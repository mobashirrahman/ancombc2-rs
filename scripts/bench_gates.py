#!/usr/bin/env python3
"""Evaluate the performance gates against ``benchmarks/results/results.jsonl``.

The gates are *continuation criteria*, not claims: a failure is reported as a
failure and does not stop the script, because the point of the report is to show
which gates are met, which are not, and by how much. ``--strict`` exits non-zero
on any failure for use in a gate that should block.

Gates (PLAN.md section 7):

    P1  kernel speed-up        >= 3x     rust-1 vs r-1core, best dataset
    P2  end-to-end speed-up    >= 2x     same, measured at the CLI boundary
    P3  sensitivity speed-up   >= 5x     on 8-16 threads, sensitivity datasets
    P4  resident memory        <= 70%    rust vs r-parallel
    P5  strong scaling         >= 0.70   speed-up from 1 to 16 threads / 15x threads

A gate with no data is reported as ``not measured`` and never as a pass. The
R arms are optional: without ``scripts/bench_r.R`` the R-relative gates cannot be
evaluated, and saying so is the honest outcome.

Usage:
    scripts/bench_gates.py --results benchmarks/results/results.jsonl
    scripts/bench_gates.py --results benchmarks/results/results.jsonl --strict
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from collections import defaultdict

GATES = {
    "P1": dict(desc="kernel speed-up, rust-1 vs r-1core", target=3.0, cmp=">="),
    "P2": dict(desc="end-to-end speed-up, rust-1 vs r-1core", target=2.0, cmp=">="),
    "P3": dict(desc="sensitivity speed-up, 8-16 threads", target=5.0, cmp=">="),
    "P4": dict(desc="peak RSS, rust vs r-parallel", target=0.70, cmp="<="),
    "P5": dict(desc="strong scaling efficiency, 1 -> 16 threads", target=0.70, cmp=">="),
}


def load(path: str) -> list[dict]:
    rows = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            rows.append(json.loads(line))
    return rows


def best(rows: list[dict], dataset: str, arm: str) -> dict | None:
    """The best (lowest wall time) successful run for an arm on a dataset."""
    cand = [
        r
        for r in rows
        if r.get("dataset") == dataset
        and r.get("arm") == arm
        and r.get("status") == "ok"
    ]
    if not cand:
        return None
    return min(cand, key=lambda r: r["wall_seconds"])


def latest_per_key(rows: list[dict]) -> list[dict]:
    """Keep the most recent row per (dataset, arm) pair.

    The results file is append-only, so a re-run of one arm must not be diluted
    by the earlier rows.
    """
    by_key: dict[tuple, dict] = {}
    for i, r in enumerate(rows):
        by_key[(r.get("dataset"), r.get("arm"))] = {**r, "_row": i}
    return sorted(by_key.values(), key=lambda r: r["_row"])


def ratio(num: dict | None, den: dict | None) -> float | None:
    if not num or not den:
        return None
    if num["wall_seconds"] <= 0:
        return None
    return den["wall_seconds"] / num["wall_seconds"]


def most_substantial(rows, datasets, arm="rust-1"):
    """The dataset a gate is judged on: the one Rust spends the longest on.

    Every gate uses this. Picking the *best* dataset per gate instead -- the
    smallest memory footprint for P4, the best speed-up for P3 and P5 -- makes
    each number a best case rather than a measurement, and a best case can pass a
    gate that the real workload fails. A gate is a continuation criterion, so it
    has to be evaluated where it is hard. The favourable datasets stay on the
    record as `best_case_*` fields, because they are real measurements too.
    """
    best_ds, best_wall = None, -1.0
    for d in datasets:
        rec = best(rows, d, arm)
        if not rec or not rec.get("wall_seconds"):
            continue
        if rec["wall_seconds"] > best_wall:
            best_ds, best_wall = d, rec["wall_seconds"]
    return best_ds


def evaluate(rows: list[dict]) -> list[dict]:
    rows = latest_per_key(rows)
    datasets = sorted({r.get("dataset") for r in rows if r.get("dataset")})
    out: list[dict] = []

    for name, g in GATES.items():
        rec: dict = {
            "gate": name,
            "desc": g["desc"],
            "target": f"{g['cmp']} {g['target']}",
            "measured": None,
            "verdict": "not measured",
            "detail": "",
            "evidence": [],
            # Per-dataset numbers, so one headline figure never stands in for
            # the whole surface.
            "surface": [],
        }
        if name in ("P1", "P2"):
            # The best dataset is the one where both arms ran and the measurement
            # is least dominated by fixed overhead, i.e. the longest rust run.
            best_ds, best_ratio = None, None
            for d in datasets:
                rr = ratio(best(rows, d, "rust-1"), best(rows, d, "r-1core"))
                if rr is not None:
                    best_ds, best_ratio = d, rr
            judged = most_substantial(
                rows, [d for d in datasets if ratio(best(rows, d, "rust-1"),
                                                    best(rows, d, "r-1core")) is not None]
            )
            if judged is not None:
                best_ds = judged
                best_ratio = ratio(best(rows, best_ds, "rust-1"), best(rows, best_ds, "r-1core"))
            if best_ratio is None:
                rec["detail"] = (
                    "needs both `rust-1` and `r-1core` on the same dataset; "
                    "run scripts/bench_r.R for the R arm"
                )
            else:
                rust = best(rows, best_ds, "rust-1")
                ref = best(rows, best_ds, "r-1core")
                rec["measured"] = best_ratio
                rec["verdict"] = "pass" if best_ratio >= g["target"] else "fail"
                rec["detail"] = (
                    f"most substantial dataset {best_ds} (least overhead-dominated): "
                    f"rust {rust['wall_seconds']:.3f}s vs R {ref['wall_seconds']:.3f}s"
                )
                # The whole surface, so the headline cannot hide that the small
                # datasets are an order of magnitude faster. A single number
                # taken from one dataset says nothing about the rest.
                rec["surface"] = surface_table(rows, datasets)
                rec["evidence"] = [f"{best_ds}/rust-1", f"{best_ds}/r-1core"]
        elif name == "P3":
            best_ds, best_ratio = None, None
            for d in datasets:
                spec = (best(rows, d, "rust-8") or {}).get("spec") or {}
                if not spec.get("sensitivity"):
                    continue
                rust = best(rows, d, "rust-8") or best(rows, d, "rust-16")
                arm = "rust-8" if best(rows, d, "rust-8") else "rust-16"
                ref = best(rows, d, "r-1core")
                rr = ratio(rust, ref)
                if rr is None:
                    continue
                if best_ratio is None or rr > best_ratio:
                    best_ds, best_ratio, best_arm = d, rr, arm
            # Judged on the most substantial dataset *that has a sensitivity
            # analysis*, which is what this gate is about.
            heaviest = most_substantial(rows, [d for d in datasets
                              if ((best(rows, d, "rust-8") or {}).get("spec") or {}).get("sensitivity")])
            if heaviest is not None:
                a = "rust-8" if best(rows, heaviest, "rust-8") else "rust-16"
                rr = ratio(best(rows, heaviest, a), best(rows, heaviest, "r-1core"))
                if rr is not None:
                    best_ds, best_ratio, best_arm = heaviest, rr, a
            if best_ratio is None:
                rec["detail"] = (
                    "needs a sensitivity dataset with `rust-8`/`rust-16` and `r-1core`"
                )
            else:
                rec["measured"] = best_ratio
                rec["verdict"] = "pass" if best_ratio >= g["target"] else "fail"
                rec["detail"] = (
                    f"most substantial sensitivity dataset {best_ds} on {best_arm}: "
                    f"{best_ratio:.2f}x"
                )
                rec["evidence"] = [f"{best_ds}/{best_arm}", f"{best_ds}/r-1core"]

            # A conservative sensitivity run is *three* independent refits, not
            # fifty: the grid is pseudo in {0, 0.1, 0.5, 1} and the 0 entry is the
            # main run, already computed. The outer level of the nesting order
            # parallelises those three, so 8-16 threads can at best give 3x on
            # that path no matter how much pool is available -- the reference's own
            # grid, not this implementation, is what fixes the width.
            #
            # So when the judged dataset is the conservative one, the headline ratio
            # is reported *alongside* the ceiling it is measured against, and every
            # sensitivity dataset's own stage-level scaling is listed. A reader who
            # sees "3.08x against a target of 5x" would otherwise conclude the
            # parallelisation is four-fifths short, when on this grid it is at the
            # maximum reachable value.
            ceiling = None
            spec = (best(rows, best_ds, "rust-8") or {}).get("spec") or {} if best_ds else {}
            if spec.get("sensitivity") and spec.get("conservative", True):
                ceiling = 3
                rec["measured_ceiling"] = ceiling
                rec["ceiling_note"] = (
                    f"{best_ds} is a conservative sensitivity run: the pseudo-count "
                    f"grid {{0, 0.1, 0.5, 1}} is {ceiling} refits after the main "
                    f"run, so {ceiling}x is the most the outer level can give. "
                    f"Measured {best_ratio:.2f}x is "
                    f"{100 * best_ratio / ceiling:.0f}% of that ceiling, against a "
                    f"target of {g['target']:.0f}x."
                )
            # Per-dataset stage scaling, so the conservative cap is visible rather
            # than inferred: a dataset with a 50-point grid should scale further,
            # and if it does not that is a different problem from the cap.
            stage_scaling = []
            for d in sorted(datasets):
                spec_d = (best(rows, d, "rust-8") or {}).get("spec") or {}
                if not spec_d.get("sensitivity"):
                    continue
                lo = best(rows, d, "rust-1") or {}
                hi = best(rows, d, "rust-8") or best(rows, d, "rust-16") or {}
                s_lo = (lo.get("run_metadata") or {}).get("stage_seconds.sensitivity")
                s_hi = (hi.get("run_metadata") or {}).get("stage_seconds.sensitivity")
                if not s_lo or not s_hi or s_hi <= 0:
                    continue
                stage_scaling.append({
                    "dataset": d,
                    "conservative": bool(spec_d.get("conservative", True)),
                    "sens_s_1thread": round(s_lo, 3),
                    "sens_s_8or16": round(s_hi, 3),
                    "stage_speedup": round(s_lo / s_hi, 3),
                    "refits": 3 if spec_d.get("conservative", True) else 50,
                })
            if stage_scaling:
                rec["sensitivity_stage_scaling"] = stage_scaling
        elif name == "P4":
            def rss_frac(d):
                rust = best(rows, d, "rust-1")
                ref = best(rows, d, "r-parallel") or best(rows, d, "r-1core")
                if not rust or not ref:
                    return None, None, None
                if not rust.get("peak_rss_kb") or not ref.get("peak_rss_kb"):
                    return None, None, None
                return rust["peak_rss_kb"] / ref["peak_rss_kb"], rust, ref

            judged = most_substantial(rows, [d for d in datasets if rss_frac(d)[0] is not None])
            fracs = {d: rss_frac(d)[0] for d in datasets if rss_frac(d)[0] is not None}
            if judged is None:
                rec["detail"] = "needs `rust-1` and an R arm with peak_rss_kb recorded"
            else:
                best_frac, rust, ref = rss_frac(judged)
                rec["measured"] = best_frac
                rec["verdict"] = "pass" if best_frac <= g["target"] else "fail"
                rec["detail"] = (
                    f"most substantial dataset {judged}: rust "
                    f"{rust['peak_rss_kb'] / 1024:.1f} MB vs R {ref['peak_rss_kb'] / 1024:.1f} MB"
                )
                # The most favourable dataset is recorded too. It is a real
                # measurement, and on a table small enough that process startup
                # dominates it will pass, which is exactly why it is not the gate.
                rec["best_case_dataset"] = min(fracs, key=lambda d: fracs[d])
                rec["best_case_fraction"] = min(fracs.values())
                rec["all_fractions"] = {d: round(v, 4) for d, v in sorted(fracs.items())}
                rec["surface"] = surface_table(rows, datasets)
                rec["evidence"] = [f"{judged}/rust-1", f"{judged}/r-parallel"]
        elif name == "P5":
            best_ds, best_eff = None, None
            for d in datasets:
                one = best(rows, d, "rust-1")
                many = best(rows, d, "rust-16")
                if not one or not many:
                    continue
                t = many["run_metadata"].get("threads") if isinstance(
                    many.get("run_metadata"), dict
                ) else None
                # The arm name carries the thread count that was requested; the
                # metadata carries the one the pool actually got.
                threads = float(t) if t else 16.0
                if threads <= 1.0:
                    continue
                speedup = one["wall_seconds"] / many["wall_seconds"]
                # Standard strong-scaling efficiency: the speed-up over the
                # thread count. `speedup / (threads - 1)` is also in circulation
                # but it exceeds 1 under perfect linear scaling, so a run at
                # 16/16 efficiency reports 1.067 and the gate's 0.7 target stops
                # meaning "70% of the ideal".
                eff = speedup / threads
                if best_eff is None or eff > best_eff:
                    best_ds, best_eff, best_threads, best_speed = d, eff, threads, speedup
            judged = most_substantial(
                rows,
                [d for d in datasets if best(rows, d, "rust-1") and best(rows, d, "rust-16")],
            )
            if judged is not None:
                one = best(rows, judged, "rust-1")
                many = best(rows, judged, "rust-16")
                t = (many.get("run_metadata") or {}).get("threads")
                threads = float(t) if t else 16.0
                if threads > 1.0:
                    speedup = one["wall_seconds"] / many["wall_seconds"]
                    best_ds, best_eff = judged, speedup / threads
                    best_threads, best_speed = threads, speedup
            if best_eff is None:
                rec["detail"] = "needs `rust-1` and `rust-16` on the same dataset"
            else:
                rec["measured"] = best_eff
                rec["verdict"] = "pass" if best_eff >= g["target"] else "fail"
                rec["detail"] = (
                    f"most substantial dataset {best_ds}: {best_speed:.3f}x on "
                    f"{best_threads:.0f} threads -> efficiency {best_eff:.3f}"
                )
                rec["surface"] = surface_table(rows, datasets)
                rec["evidence"] = [f"{best_ds}/rust-1", f"{best_ds}/rust-16"]
        out.append(rec)
    return out


def surface_table(rows, datasets):
    """The per-dataset numbers behind a headline figure.

    Each row is one dataset: the single-core speed-up against R, the peak-RSS
    ratio, and the strongest observed thread scaling. A reader who wants to know
    whether a gate number is representative can see every dataset here.
    """
    out = []
    for d in datasets:
        rust = best(rows, d, "rust-1")
        ref = best(rows, d, "r-1core")
        rpar = best(rows, d, "r-parallel")
        speed = ratio(rust, ref)
        rss = None
        if rust and (rpar or ref):
            base = rpar or ref
            if base.get("peak_rss_kb"):
                rss = rust["peak_rss_kb"] / base["peak_rss_kb"]
        scale = None
        one = rust
        for arm in ("rust-2", "rust-4", "rust-8", "rust-16", "rust-32"):
            m = best(rows, d, arm)
            if one and m and one.get("wall_seconds"):
                s = one["wall_seconds"] / m["wall_seconds"]
                if scale is None or s > scale:
                    scale = s
        out.append({
            "dataset": d,
            "rust_1_wall_s": round(rust["wall_seconds"], 3) if rust else None,
            "r_1core_wall_s": round(ref["wall_seconds"], 3) if ref else None,
            "speedup_vs_r": round(speed, 3) if speed else None,
            "peak_rss_ratio_vs_r": round(rss, 3) if rss else None,
            "best_thread_scaling": round(scale, 3) if scale else None,
        })
    return out


def missing_from_surface(rows, datasets):
    """Datasets the gate could not use, and why.

    A dataset whose Rust arms all failed used to vanish silently: it contributed no
    rows, the "most substantial dataset" was then decided among whatever was left,
    and the report named a smaller dataset as the headline without ever saying a
    larger one had been dropped. That is the worst way for a gate to be wrong --
    it reads as a measurement and is a selection.

    `bm5` is the live case: 5,000 taxa x 20,000 samples with the conservative
    sensitivity analysis, and the Rust arm's resident set exceeds this host's 31 GB,
    so every arm is killed. The gate must say so rather than quietly promote
    `bm6`.
    """
    out = []
    # `results.jsonl` is append-only across runs, so an arm that succeeded in an
    # earlier run and fails now has *both* rows present. Judging on the raw list
    # would keep reporting the old success and hide the current failure -- which is
    # the opposite of what this function is for. Only the latest row per
    # (dataset, arm) counts.
    latest = latest_per_key(rows)
    for d in datasets:
        rust = [r for r in latest if r.get("dataset") == d and r.get("arm", "").startswith("rust")]
        ref = [r for r in latest if r.get("dataset") == d and r.get("arm", "").startswith("r-")]
        if not rust or not ref:
            out.append({
                "dataset": d,
                "rust_rows": len(rust),
                "r_rows": len(ref),
                "statuses": sorted({str(r.get("status")) for r in rust + ref}),
            })
            continue
        ok = [r for r in rust if r.get("status") == "ok"]
        if not ok:
            out.append({
                "dataset": d,
                "rust_rows": len(rust),
                "r_rows": len(ref),
                "statuses": sorted({str(r.get("status")) for r in rust}),
                "exit_codes": sorted({str(r.get("exit_code")) for r in rust if "exit_code" in r}),
            })
    return out


def main(argv: list[str]) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--results", default="benchmarks/results/results.jsonl")
    ap.add_argument("--json", action="store_true", help="emit JSON instead of text")
    ap.add_argument("--strict", action="store_true",
                    help="exit non-zero if any gate fails")
    # Defaults to the committed report. Without this, `bench_gates.py --json`
    # printed the verdict to stdout and left `gates.json` holding whatever the
    # last `make gates` wrote -- which is how a stale "2 pass, 3 fail" survived a
    # re-run that reported "0 pass, 5 fail".
    ap.add_argument(
        "--report", default="benchmarks/results/gates.json",
        help="where to write the JSON report ('' to skip)",
    )
    args = ap.parse_args(argv)

    try:
        rows = load(args.results)
    except OSError as e:
        print(f"cannot read {args.results}: {e}", file=sys.stderr)
        return 2
    if not rows:
        print(f"{args.results} is empty; run scripts/run_benchmarks.py first",
              file=sys.stderr)
        return 2

    report = evaluate(rows)
    # Datasets with at least one *successful* Rust arm. Counting the ones that
    # merely have rows would put bm5 -- five killed arms -- in the denominator of
    # "measured", which is the number a reader would trust.
    _latest = latest_per_key(rows)
    datasets = sorted({
        r.get("dataset") for r in _latest
        if r.get("dataset") and r.get("status") == "ok"
        and str(r.get("arm", "")).startswith("rust")
    })
    # The full declared surface, so a dataset that failed to run is named rather
    # than dropped from the comparison without a word.
    index_path = os.path.join(os.path.dirname(args.results), "..", "datasets", "index.json")
    declared = []
    try:
        with open(index_path) as f:
            declared = sorted(json.load(f))
    except (OSError, ValueError):
        declared = datasets
    absent = missing_from_surface(rows, declared)
    for rec in report:
        rec["absent_datasets"] = absent
        rec["n_datasets_declared"] = len(declared)
        rec["n_datasets_measured"] = len(datasets)

    n_fail = sum(1 for r in report if r["verdict"] == "fail")
    n_pass = sum(1 for r in report if r["verdict"] == "pass")
    n_none = sum(1 for r in report if r["verdict"] == "not measured")

    # Built once, so `--json` and `--report` cannot disagree and so the file is
    # written whether the caller asked for the text report, the JSON, or both.
    payload = {
        "gates": report,
        "absent_datasets": absent,
        "pass": n_pass,
        "fail": n_fail,
        "not_measured": n_none,
        "results_file": args.results,
        "note": (
            "Gates are continuation criteria, not claims. `not measured` is "
            "reported as such and never counted as a pass."
        ),
    }

    if args.json:
        print(json.dumps(payload, indent=1))
    else:
        print("ancombc2-rs performance gates")
        print("=" * 72)
        for r in report:
            mark = {"pass": "PASS", "fail": "FAIL", "not measured": "----"}[r["verdict"]]
            m = "n/a" if r["measured"] is None else f"{r['measured']:.3f}"
            print(f"[{mark}] {r['gate']}  {r['desc']}")
            print(f"       target {r['target']:<12} measured {m}")
            if r["detail"]:
                print(f"       {r['detail']}")
            for row in r.get("surface") or []:
                if not row["rust_1_wall_s"]:
                    continue
                # A failed re-run shadows the earlier good row in
                # `latest_per_key`, so `r_1core_wall_s` can be None while
                # `rust_1_wall_s` is present -- which is exactly what a
                # container run with broken R arms produces. Print n/a rather
                # than crashing; the gate itself already reports the missing
                # reference through `speedup_vs_r is None`.
                sr = "n/a" if row["speedup_vs_r"] is None else f"{row['speedup_vs_r']:.2f}x"
                mr = "n/a" if row["peak_rss_ratio_vs_r"] is None else f"{row['peak_rss_ratio_vs_r']:.2f}x"
                sc = "n/a" if row["best_thread_scaling"] is None else f"{row['best_thread_scaling']:.2f}x"
                rr = "n/a" if row["r_1core_wall_s"] is None else f"{row['r_1core_wall_s']:>8.3f}s"
                print(
                    f"         {row['dataset']:>5}: rust-1 {row['rust_1_wall_s']:>8.3f}s"
                    f"  R {rr:>10}  speed {sr:>7}"
                    f"  rss {mr:>7}  scaling {sc:>7}"
                )
        print("=" * 72)
    if absent:
        print()
        print("datasets in the declared surface with no usable measurement:")
        for a in absent:
            print(f"  {a['dataset']}: rust rows {a['rust_rows']}, R rows {a['r_rows']}, "
                  f"statuses {a['statuses']}"
                  + (f", exit codes {a['exit_codes']}" if "exit_codes" in a else ""))
        print("  every gate above is judged on the datasets that did run, so the")
        print("  headline is NOT the most substantial dataset in the surface.")
    print()
    print(f"{n_pass} pass, {n_fail} fail, {n_none} not measured")
    if n_none:
        print(
            "A gate with no data is not a pass. `not measured` means the arm it "
            "needs was not run; see docs/compatibility.md."
        )

    if args.report:
        os.makedirs(os.path.dirname(args.report) or ".", exist_ok=True)
        with open(args.report, "w") as f:
            json.dump(payload, f, indent=1)
        print(f"wrote {args.report}")

    if args.strict and n_fail:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
