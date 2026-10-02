#!/usr/bin/env python3
"""Generate the benchmark datasets.

One dataset per point of the scaling surface in ``PLAN.md`` section 7. Each is a
counts matrix, a metadata table, a formula and a JSON config, in exactly the
layout the CLI reads and the oracle harness reads, so a benchmark measures the
same code path a user runs.

The generator is Python rather than R on purpose: the benchmark must be
reproducible from the repository without an R installation, and the *numbers*
only have to be plausible count data -- the benchmark measures throughput, not
statistical content. Parity is the golden suite's job.

Usage:
    scripts/make_bench_data.py --out benchmarks/datasets
    scripts/make_bench_data.py --out benchmarks/datasets --only bm5
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time

try:
    import numpy as np
except ImportError:  # pragma: no cover - exercised only without numpy
    np = None

# name -> (n_taxon, n_sample, n_group, n_continuous, zero_rate, sensitivity,
#          conservative, fix_eff terms)
#
# The surface is a *scaling* surface, so the shapes are deliberately far apart
# and the covariates grow with the group count, which is what a real
# microbiome study looks like. `zero_rate` is the fraction of counts forced to
# zero, which drives the missingness-pattern grouping and therefore the
# factorisation count.
SURFACE: dict[str, dict] = {
    "bm1": dict(n_taxon=50, n_sample=500, n_group=3, n_cont=1, zero_rate=0.10,
                sensitivity=False, conservative=True),
    "bm2": dict(n_taxon=100, n_sample=1000, n_group=5, n_cont=3, zero_rate=0.10,
                sensitivity=False, conservative=True),
    "bm3": dict(n_taxon=500, n_sample=5000, n_group=5, n_cont=3, zero_rate=0.50,
                sensitivity=False, conservative=True),
    "bm4": dict(n_taxon=1000, n_sample=10000, n_group=10, n_cont=8, zero_rate=0.10,
                sensitivity=False, conservative=True),
    "bm5": dict(n_taxon=5000, n_sample=20000, n_group=10, n_cont=8, zero_rate=0.10,
                sensitivity=True, conservative=True),
    "bm6": dict(n_taxon=1000, n_sample=10000, n_group=10, n_cont=8, zero_rate=0.10,
                sensitivity=True, conservative=False),
}


class Rng:
    """A splitmix64 stream, so a dataset is byte-identical on any platform.

    Python's `random` is Mersenne Twister and stable, but its `gauss` caches a
    spare value, which makes the sequence depend on call history. splitmix64 is
    a counter, so the value at step *k* does not depend on how many values were
    drawn before it -- which is what a benchmark generator needs.
    """

    def __init__(self, seed: int) -> None:
        self.state = seed & 0xFFFFFFFFFFFFFFFF

    def next_u64(self) -> int:
        self.state = (self.state + 0x9E3779B97F4A7C15) & 0xFFFFFFFFFFFFFFFF
        z = self.state
        z = ((z ^ (z >> 30)) * 0xBF58476D1CE4E5B9) & 0xFFFFFFFFFFFFFFFF
        z = ((z ^ (z >> 27)) * 0x94D049BB133111EB) & 0xFFFFFFFFFFFFFFFF
        return z ^ (z >> 31)

    def uniform(self) -> float:
        return (self.next_u64() >> 11) / float(1 << 53)

    def normal(self) -> float:
        """Box-Muller. Two uniforms per normal, no cache."""
        u1 = max(self.uniform(), 1e-300)
        u2 = self.uniform()
        return (-2.0 * __import__("math").log(u1)) ** 0.5 * __import__(
            "math"
        ).cos(2.0 * __import__("math").pi * u2)


def generate(spec: dict, seed: int, out_dir: str, name: str) -> dict:
    import math

    n_tax = spec["n_taxon"]
    n_samp = spec["n_sample"]
    n_group = spec["n_group"]
    n_cont = spec["n_cont"]
    zero_rate = spec["zero_rate"]
    rng = Rng(seed)

    os.makedirs(out_dir, exist_ok=True)
    counts_path = os.path.join(out_dir, f"{name}.counts.tsv")
    meta_path = os.path.join(out_dir, f"{name}.meta.tsv")
    config_path = os.path.join(out_dir, f"{name}.config.json")

    # A group effect on a quarter of the taxa, so the run has real signal to
    # find and the sensitivity analysis has something to be robust about.
    signal = [1.5 if (i % 4 == 0) else 0.0 for i in range(n_tax)]

    # Counts first, because the file is large and the metadata depends on nothing
    # in it -- writing counts first also means a partial run leaves the metadata
    # absent rather than stale.
    #
    # numpy's `Generator` is used when available: bm5 is 100 million cells, and
    # a scalar Python loop over them takes tens of minutes per dataset. The
    # splitmix64 stream above stays as the fallback so a machine without numpy
    # still produces a *different but valid* dataset rather than failing; the
    # dataset is byte-reproducible only with numpy present, which is recorded in
    # the config's `generator` field.
    with open(counts_path, "w", buffering=1 << 22) as f:
        f.write("\t".join(f"sample_{j + 1:06d}" for j in range(n_samp)) + "\n")
        if np is not None:
            gen = np.random.default_rng(seed)
            grp = (np.arange(n_samp) * n_group // n_samp).astype(np.float64)
            # `log(mu) = log(base) + signal * grp / (n_group - 1)`, computed
            # once per taxon, so the per-cell work is two multiplies.
            for i in range(n_tax):
                log_base = gen.normal(0.0, 0.8) + 4.0
                log_mu = log_base + signal[i] * grp / max(n_group - 1, 1)
                mu = np.exp(log_mu)
                noise = gen.normal(0.0, 1.0, n_samp)
                counts_row = (mu + noise * np.sqrt(mu)).astype(np.int64)
                np.maximum(counts_row, 0, out=counts_row)
                if zero_rate > 0:
                    drop = gen.random(n_samp) < zero_rate
                    counts_row[drop] = 0
                f.write(f"taxon_{i + 1:06d}\t")
                f.write("\t".join(map(str, counts_row.tolist())))
                f.write("\n")
        else:
            for i in range(n_tax):
                base = math.exp(rng.normal() * 0.8 + 4.0)
                row = [f"taxon_{i + 1:06d}"]
                for j in range(n_samp):
                    g = j * n_group // n_samp
                    mu = base * math.exp(signal[i] * g / max(n_group - 1, 1))
                    v = rng.normal() * math.sqrt(mu) + mu
                    c = int(v) if v > 0 else 0
                    if rng.uniform() < zero_rate:
                        c = 0
                    row.append(str(c))
                f.write("\t".join(row) + "\n")

    with open(meta_path, "w", buffering=1 << 22) as f:
        cols = ["group"] + [f"x{k + 1}" for k in range(n_cont)]
        f.write("\t".join(cols) + "\n")
        if np is not None:
            gen = np.random.default_rng(seed ^ 0x5EED)
            block = np.round(gen.normal(0.0, 1.0, (n_samp, n_cont)), 17)
            for j in range(n_samp):
                g = j * n_group // n_samp + 1
                # 17 significant digits, so the design round-trips exactly; see
                # the note in reference/R/fixtures.R.
                vals = [str(g)] + [repr(float(v)) for v in block[j]]
                f.write("\t".join([f"sample_{j + 1:06d}"] + vals) + "\n")
        else:
            for j in range(n_samp):
                g = j * n_group // n_samp + 1
                vals = [str(g)] + [
                    format(rng.normal(), ".17g") for _ in range(n_cont)
                ]
                f.write("\t".join([f"sample_{j + 1:06d}"] + vals) + "\n")

    terms = ["group"] + [f"x{k + 1}" for k in range(n_cont)]
    with open(os.path.join(out_dir, f"{name}.formula.txt"), "w") as f:
        f.write(" ".join(terms) + "\n")
    with open(config_path, "w") as f:
        json.dump(
            dict(
                name=name,
                fix_formula=" ".join(terms),
                group="group",
                p_adj_method="holm",
                pseudo=0.0,
                pseudo_sens=spec["sensitivity"],
                conservative=spec["conservative"],
                prv_cut=0.1,
                lib_cut=0.0,
                s0_perc=0.05,
                struc_zero=False,
                neg_lb=False,
                alpha=0.05,
                do_global=n_group >= 3,
                do_pairwise=n_group >= 3,
                seed=seed,
                generator="numpy" if np is not None else "python-splitmix64",
            ),
            f,
            indent=1,
        )
    return dict(
        name=name,
        n_taxon=n_tax,
        n_sample=n_samp,
        n_group=n_group,
        n_cont=n_cont,
        zero_rate=zero_rate,
        sensitivity=spec["sensitivity"],
        conservative=spec["conservative"],
        fix_eff=len(terms) + 1,  # the intercept plus one column per term
        sensitivity_runs=(4 if spec["conservative"] else 50) if spec["sensitivity"] else 0,
    )


def main(argv: list[str]) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", default="benchmarks/datasets")
    ap.add_argument("--only", action="append", default=[],
                    help="generate only these datasets (repeatable)")
    ap.add_argument("--seed", type=int, default=20240101)
    ap.add_argument("--force", action="store_true",
                    help="regenerate even if the files already exist")
    args = ap.parse_args(argv)

    names = args.only or list(SURFACE)
    unknown = [n for n in names if n not in SURFACE]
    if unknown:
        print(f"unknown dataset(s): {', '.join(unknown)}", file=sys.stderr)
        print(f"available: {', '.join(SURFACE)}", file=sys.stderr)
        return 2

    index_path = os.path.join(args.out, "index.json")
    index = {}
    if os.path.exists(index_path) and not args.force:
        with open(index_path) as f:
            index = json.load(f)

    for name in names:
        spec = SURFACE[name]
        counts_path = os.path.join(args.out, f"{name}.counts.tsv")
        if os.path.exists(counts_path) and not args.force:
            # Register the existing dataset in the index rather than skipping it
            # entirely: the index is what the benchmark runner reads, and a
            # dataset that was generated by an earlier invocation must not go
            # missing from it just because this one did not regenerate it.
            if name not in index:
                with open(os.path.join(args.out, f"{name}.config.json")) as f:
                    cfg = json.load(f)
                # The formula may be written with `+` or with whitespace, since
                # R's `terms()` accepts both.
                terms = [
                    t for t in cfg["fix_formula"].replace("~", "").replace("+", " ").split()
                    if t
                ]
                index[name] = dict(
                    name=name,
                    n_taxon=spec["n_taxon"],
                    n_sample=spec["n_sample"],
                    n_group=spec["n_group"],
                    n_cont=spec["n_cont"],
                    zero_rate=spec["zero_rate"],
                    sensitivity=spec["sensitivity"],
                    conservative=spec["conservative"],
                    fix_eff=len(terms) + 1,
                    sensitivity_runs=(
                        (4 if spec["conservative"] else 50) if spec["sensitivity"] else 0
                    ),
                    generator=cfg.get("generator", "unknown"),
                )
                print(f"{name}: exists, registered in the index")
            else:
                print(f"{name}: exists, skipping (use --force to regenerate)")
            continue
        t0 = time.perf_counter()
        entry = generate(spec, args.seed + int(name[2:]), args.out, name)
        dt = time.perf_counter() - t0
        index[name] = entry
        size_mb = os.path.getsize(counts_path) / 1e6
        print(
            f"{name}: {entry['n_taxon']} taxa x {entry['n_sample']} samples, "
            f"{entry['n_group']} groups, {entry['fix_eff']} fix_eff, "
            f"zero_rate {entry['zero_rate']}, sens {entry['sensitivity_runs']} runs, "
            f"{size_mb:.1f} MB, {dt:.1f}s"
        )

    os.makedirs(args.out, exist_ok=True)
    with open(index_path, "w") as f:
        json.dump(index, f, indent=1, sort_keys=True)
    print(f"wrote {index_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
