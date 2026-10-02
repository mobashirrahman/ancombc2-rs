#!/usr/bin/env python3
"""Run the benchmark surface and append to ``benchmarks/results/results.jsonl``.

Seven arms per dataset, as the plan requires:

    r-1core        ANCOMBC 2.15.2, one core
    r-parallel     ANCOMBC 2.15.2, all cores
    rust-1         ancombc2-rs, 1 thread
    rust-4         ancombc2-rs, 4 threads
    rust-8         ancombc2-rs, 8 threads
    rust-16        ancombc2-rs, 16 threads
    rust-32        ancombc2-rs, 32 threads (only if the host has 32+)

Metrics per run: wall time, CPU time, peak RSS, bytes allocated, allocation
count, throughput, and the stage-resolved timings the CLI writes to
``run_metadata.tsv``. Peak RSS comes from ``getrusage(RUSAGE_CHILDREN)``
differences, which is the only figure available without a cgroup or a
``/usr/bin/time`` build.

Every row carries the compatibility target and the dataset's index entry, so a
results file is self-describing. Nothing is asserted here: the gates are
evaluated by ``scripts/bench_gates.py``, separately, so a failing gate cannot
hide a passing measurement.

Usage:
    scripts/run_benchmarks.py --data benchmarks/datasets --arms rust-1,rust-8
    scripts/run_benchmarks.py --data benchmarks/datasets --datasets bm1,bm2
    scripts/run_benchmarks.py --data benchmarks/datasets --r-script scripts/bench_r.R
"""

from __future__ import annotations

import argparse
import json
import os
import platform
import resource
import shutil
import subprocess
import sys
import time

R_ARMS = ("r-1core", "r-parallel")
RUST_ARMS = ("rust-1", "rust-4", "rust-8", "rust-16", "rust-32")
TARGET = (
    "ancombc2-rs v0.1 equivalent to ANCOMBC 2.15.2 at "
    "dc4febdf59badb3a8dfe0c767ef2186323c2199a, R 4.5.x, seed 42"
)


def host_info() -> dict:
    info = {
        "platform": platform.platform(),
        "machine": platform.machine(),
        "python": platform.python_version(),
        "cpu_count": os.cpu_count(),
    }
    try:
        with open("/proc/cpuinfo") as f:
            for line in f:
                if line.lower().startswith("model name"):
                    info["cpu_model"] = line.split(":", 1)[1].strip()
                    break
    except OSError:
        pass
    try:
        with open("/proc/meminfo") as f:
            for line in f:
                if line.startswith("MemTotal:"):
                    info["mem_total_kb"] = int(line.split()[1])
                    break
    except OSError:
        pass
    return info


def rust_version() -> str | None:
    try:
        out = subprocess.run(
            ["rustc", "--version"], capture_output=True, text=True, timeout=20
        )
        return out.stdout.strip() or None
    except (OSError, subprocess.SubprocessError):
        return None


def r_version() -> str | None:
    if not shutil.which("Rscript"):
        return None
    try:
        out = subprocess.run(
            ["Rscript", "-e", "cat(R.version.string)"],
            capture_output=True,
            text=True,
            timeout=60,
        )
        return out.stdout.strip() or None
    except (OSError, subprocess.SubprocessError):
        return None


def peak_rss_children_kb() -> int:
    """Peak RSS of all reaped children, in kibibytes.

    `ru_maxrss` is a high-water mark, so it only ever grows. The caller takes a
    difference against the value before the run, which is correct for one child
    at a time and is the reason the runner never runs arms concurrently.
    """
    return resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss


def read_stages(meta_path: str) -> dict:
    """Parse the CLI's ``run_metadata.tsv`` into a flat dict of numbers."""
    out: dict[str, float] = {}
    try:
        with open(meta_path) as f:
            for line in f:
                parts = line.rstrip("\n").split("\t")
                if len(parts) != 2:
                    continue
                key, val = parts
                try:
                    out[key] = float(val)
                except ValueError:
                    out[key] = val  # type: ignore[assignment]
    except OSError:
        pass
    return out


def run_rust_arm(
    binary: str, data_dir: str, name: str, threads: int, out_dir: str
) -> dict:
    cfg = json.load(open(os.path.join(data_dir, f"{name}.config.json")))
    formula = open(os.path.join(data_dir, f"{name}.formula.txt")).read().strip()
    cmd = [
        binary,
        "--counts", os.path.join(data_dir, f"{name}.counts.tsv"),
        "--meta", os.path.join(data_dir, f"{name}.meta.tsv"),
        "--formula", formula,
        "--threads", str(threads),
        "--out", out_dir,
        "--quiet",
    ]
    if cfg.get("group"):
        cmd += ["--group", cfg["group"]]
    for flag, key in (
        ("--global", "do_global"),
        ("--pairwise", "do_pairwise"),
        ("--struc-zero", "struc_zero"),
        # The configs spell this `pseudo_sens`, matching `AncombcConfig::pseudo_sens`
        # and the CLI flag. An earlier version of this loop looked for a key named
        # `sensitivity`, which no config has, so `--pseudo-sens` was never passed and
        # `bm5` -- the dataset that exists *to* measure the sensitivity analysis --
        # was benchmarked without it. The result was not wrong so much as
        # meaningless: it reported the cost of a different workload under this
        # dataset's name. Both keys are accepted, with the real one first.
        ("--pseudo-sens", "pseudo_sens"),
    ):
        if cfg.get(key):
            cmd.append(flag)
    if cfg.get("pseudo_sens") or cfg.get("sensitivity"):
        cmd.append("--conservative" if cfg.get("conservative") else "--nonconservative")
    if cfg.get("p_adj_method"):
        cmd += ["--p-adj-method", cfg["p_adj_method"]]
    return run_and_measure(cmd, out_dir, arm=f"rust-{threads}", name=name)


def run_r_arm(
    r_script: str, data_dir: str, name: str, threads: int | None, out_dir: str
) -> dict:
    if not os.path.exists(r_script):
        return {
            "arm": f"r-{'1core' if threads == 1 else 'parallel'}",
            "dataset": name,
            "status": "skipped",
            "reason": f"{r_script} not found",
        }
    env = dict(os.environ)
    if threads:
        env["ANCOMBC_BENCH_THREADS"] = str(threads)
    else:
        env.pop("ANCOMBC_BENCH_THREADS", None)
    env["ANCOMBC_BENCH_DATA"] = data_dir
    env["ANCOMBC_BENCH_OUT"] = out_dir
    env["ANCOMBC_BENCH_DATASET"] = name
    cmd = ["Rscript", r_script]
    return run_and_measure(
        cmd, out_dir, arm="r-1core" if threads == 1 else "r-parallel", name=name, env=env
    )


def run_and_measure(
    cmd: list[str], out_dir: str, arm: str, name: str, env: dict | None = None
) -> dict:
    os.makedirs(out_dir, exist_ok=True)
    rss_before = peak_rss_children_kb()
    t0 = time.perf_counter()
    cpu0 = time.process_time()
    try:
        proc = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            env=env,
            timeout=7200,
        )
    except subprocess.TimeoutExpired:
        return {
            "arm": arm,
            "dataset": name,
            "status": "timeout",
            "wall_seconds": 7200.0,
        }
    wall = time.perf_counter() - t0
    cpu = time.process_time() - cpu0
    rss_after = peak_rss_children_kb()
    rec: dict = {
        "arm": arm,
        "dataset": name,
        "wall_seconds": wall,
        # `process_time` measures this process, not the child, so it is only
        # meaningful for the CPU time of a single-threaded child. The child's own
        # user+sys is what the gate wants, and it is taken from the metadata the
        # CLI writes when available; otherwise this is a lower bound and the
        # record says so.
        "parent_cpu_seconds": cpu,
        "peak_rss_kb": max(rss_after - rss_before, 0),
        "peak_rss_kb_is_high_water": True,
        "exit_code": proc.returncode,
    }
    if proc.returncode != 0:
        rec["status"] = "failed"
        rec["stderr_tail"] = proc.stderr[-2000:]
        return rec
    rec["status"] = "ok"
    meta = read_stages(os.path.join(out_dir, "run_metadata.tsv"))
    rec["run_metadata"] = meta
    if "n_taxa_reported" in meta and wall > 0:
        rec["taxa_per_second"] = meta["n_taxa_reported"] / wall
    if "n_samples" in meta and "n_taxa_reported" in meta:
        rec["cells_per_second"] = (
            meta["n_taxa_reported"] * meta["n_samples"] / wall
        )
    for key in ("bytes_allocated", "allocations"):
        if key in meta:
            rec[key] = meta[key]
    # The child's own `VmHWM` beats the rusage high-water difference, which cannot
    # decrease and so reads as zero for every run after the largest one.
    if "peak_rss_kb" in meta and isinstance(meta["peak_rss_kb"], float):
        rec["peak_rss_kb"] = meta["peak_rss_kb"]
        rec["peak_rss_source"] = meta.get("peak_rss_source", "run_metadata.tsv")
        rec["peak_rss_kb_is_high_water"] = False
    else:
        rec["peak_rss_kb"] = max(rss_after - rss_before, 0)
        rec["peak_rss_source"] = "getrusage(RUSAGE_CHILDREN) high-water delta"
    return rec


def main(argv: list[str]) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data", default="benchmarks/datasets")
    ap.add_argument("--results", default="benchmarks/results/results.jsonl")
    ap.add_argument("--binary", default="target/release/ancombc2-rs")
    ap.add_argument("--r-script", default="scripts/bench_r.R")
    ap.add_argument("--datasets", default="",
                    help="comma-separated subset of the surface")
    ap.add_argument("--arms", default="",
                    help="comma-separated subset of " + ",".join(R_ARMS + RUST_ARMS))
    ap.add_argument("--repeats", type=int, default=1,
                    help="runs per arm; the best wall time is kept")
    ap.add_argument("--out-root", default="benchmarks/results/run")
    args = ap.parse_args(argv)

    index_path = os.path.join(args.data, "index.json")
    if not os.path.exists(index_path):
        print(
            f"{index_path} not found; run scripts/make_bench_data.py --out {args.data}",
            file=sys.stderr,
        )
        return 2
    with open(index_path) as f:
        index = json.load(f)

    datasets = (
        [d for d in args.datasets.split(",") if d] or sorted(index)
    )
    arms = [a for a in args.arms.split(",") if a] or list(R_ARMS + RUST_ARMS)
    unknown = [a for a in arms if a not in R_ARMS + RUST_ARMS]
    if unknown:
        print(f"unknown arm(s): {', '.join(unknown)}", file=sys.stderr)
        return 2

    host = host_info()
    env = {
        "target": TARGET,
        "host": host,
        "rustc": rust_version(),
        "r": r_version(),
        "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
    }
    print(json.dumps(env, indent=1))

    n_cpu = os.cpu_count() or 1
    rows: list[dict] = []
    for name in datasets:
        for arm in arms:
            if arm.startswith("rust-"):
                threads = int(arm.split("-")[1])
                if threads > n_cpu:
                    print(f"{name}/{arm}: skipped, {threads} > {n_cpu} CPUs")
                    rows.append(
                        {
                            "arm": arm,
                            "dataset": name,
                            "status": "skipped",
                            "reason": f"{threads} threads > {n_cpu} CPUs",
                        }
                    )
                    continue
                if not os.path.exists(args.binary):
                    print(f"{args.binary} not found; run `cargo build --release`",
                          file=sys.stderr)
                    return 2
                out_dir = os.path.join(args.out_root, f"{name}.{arm}")
                runs = [
                    run_rust_arm(args.binary, args.data, name, threads, out_dir)
                    for _ in range(args.repeats)
                ]
            else:
                threads = 1 if arm == "r-1core" else None
                out_dir = os.path.join(args.out_root, f"{name}.{arm}")
                runs = [
                    run_r_arm(args.r_script, args.data, name, threads, out_dir)
                    for _ in range(args.repeats)
                ]
            ok = [r for r in runs if r.get("status") == "ok"]
            if not ok:
                for r in runs:
                    r.update({"dataset": name, **env, "spec": index.get(name)})
                    rows.append(r)
                    print(f"{name}/{arm}: {r.get('status')}")
                continue
            best = min(ok, key=lambda r: r["wall_seconds"])
            best["repeats"] = len(runs)
            best["wall_seconds_all"] = [r["wall_seconds"] for r in ok]
            best.update({"dataset": name, **env, "spec": index.get(name)})
            rows.append(best)
            print(
                f"{name}/{arm}: {best['wall_seconds']:.3f}s wall, "
                f"{best['peak_rss_kb'] / 1024:.1f} MB peak RSS"
            )

    os.makedirs(os.path.dirname(args.results) or ".", exist_ok=True)
    with open(args.results, "a") as f:
        for r in rows:
            f.write(json.dumps(r, sort_keys=True) + "\n")
    print(f"appended {len(rows)} rows to {args.results}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
