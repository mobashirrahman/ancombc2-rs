#!/usr/bin/env python3
"""Fail when the two simulation arms are not analysing the same thing.

Why this exists
---------------
The grid runner once passed ``struc_zero = TRUE`` and ``neg_lb = TRUE`` to the
oracle in ``scripts/sim_r.R`` and left both at ``AncombcConfig::default()``'s
``false`` in the Rust arm. The two arms were therefore running *different
analyses* on identical data. Nothing failed: every cell produced a number, and
``ancombc2-sim summarise`` compared them, because the rows it compares look the
same. The symptom only appeared when the numbers were read as a result -- at 90%
zero inflation the oracle retained 11 of 500 taxa and Rust retained all 500, and
a "FDR 0.63" cell turned out to be entirely an artefact of the mismatch.

The same class of bug appeared earlier in ``scripts/run_benchmarks.py``, which
looked for a config key named ``sensitivity`` while every config spells it
``pseudo_sens``, so the sensitivity analysis never ran in any benchmark arm.

Both were the same mistake: the analysis configuration was written down twice, in
two languages, and nothing compared the two copies. The fix is that the grid file
is now the single source and both arms read it. This script is the check that it
stays that way.

What it compares
----------------
Not everything -- the oracle's per-taxon values are not reproducible by
construction in R, and the acceptance rule for the *results* is
``ancombc2-sim summarise``'s job. This compares the quantities that must agree if
and only if both arms ran the same analysis:

* ``n_retained`` and ``n_da_retained``, which are decided entirely by the filter
  configuration (``prv_cut``, ``lib_cut``, ``struc_zero``, ``neg_lb``);
* ``n_diff_abn``, which is decided by the inference configuration on top of that.

A tolerance of zero is right here: these are counts of taxa, not estimates, and
any difference at all means the arms were configured differently.

What this does *not* check
--------------------------
It compares the arms against each other, so it detects divergence, not a
configuration that is wrong in *both* arms at once. Passing the same wrong grid
to both arms passes this check, correctly -- they really do agree. What prevents
that case is structural rather than observable: there is now one place the
configuration is written down, `grid.json`, and both arms read it, so there is no
second copy to be wrong.
"""

from __future__ import annotations

import argparse
import json
import os
import pathlib
import subprocess
import sys
import tempfile


def run_arm(cmd: list[str], label: str) -> list[dict]:
    print(f"  {label} ...", flush=True)
    env = dict(os.environ)
    # The oracle is sourced, not installed, and `jsonlite` lives in the user
    # library rather than the site library on this host. Both are the caller's to
    # provide -- this script does not guess a path.
    oracle = os.environ.get("ANCOMBC_ORACLE_DIR")
    if oracle:
        env["ANCOMBC_ORACLE_DIR"] = oracle
    proc = subprocess.run(cmd, capture_output=True, text=True, env=env)
    if proc.returncode != 0:
        sys.stderr.write(proc.stdout + proc.stderr)
        raise SystemExit(f"{label} failed with exit {proc.returncode}")
    return [json.loads(line) for line in proc.stdout.splitlines() if line.strip().startswith("{")]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--grid", default="validation/simulation/full/grid.json")
    ap.add_argument("--data", default="sim-data-full")
    ap.add_argument("--rust-bin", default="target/release/ancombc2-sim")
    ap.add_argument("--cells", default="0,60,120",
                    help="cell indices to compare; pick cells where both arms "
                         "retain taxa, so an all-zero comparison cannot pass by "
                         "both arms doing nothing")
    ap.add_argument("--reps", type=int, default=5)
    ap.add_argument("--oracle-dir", default="reference/ANCOMBC")
    args = ap.parse_args()

    repo = pathlib.Path(__file__).resolve().parent.parent
    grid = repo / args.grid
    data = repo / args.data
    if not grid.exists():
        raise SystemExit(f"no grid at {grid}")
    if not data.exists():
        raise SystemExit(
            f"no simulated data at {data}\n"
            f"  run: {args.rust_bin} generate --grid {args.grid} --out {args.data}"
        )

    tmp = pathlib.Path(tempfile.mkdtemp(prefix="sim-agree-"))
    failures: list[str] = []

    for cell in [int(c) for c in args.cells.split(",") if c.strip()]:
        rust_out = tmp / f"cell{cell}.rust.jsonl"
        r_out = tmp / f"cell{cell}.r.jsonl"
        print(f"cell {cell}:")
        run_arm([args.rust_bin, "rust", "--grid", str(grid), "--cells", str(cell),
                 "--reps", str(args.reps), "--out", str(rust_out)], "rust arm")
        run_arm(["Rscript", "--vanilla", "scripts/sim_r.R",
                 "--data", str(data), "--grid", str(grid), "--out", str(r_out),
                 "--cells", str(cell), "--reps", str(args.reps)], "R arm")

        rust = [json.loads(l) for l in rust_out.read_text().splitlines() if l.strip()]
        rr = [json.loads(l) for l in r_out.read_text().splitlines() if l.strip()]
        if not rust or not rr:
            failures.append(f"cell {cell}: an arm produced no rows")
            continue

        for key in ("n_retained", "n_da_retained", "n_diff_abn"):
            rv = sum(x[key] for x in rust if x.get(key) is not None)
            av = sum(x[key] for x in rr if x.get(key) is not None)
            if rv != av:
                failures.append(
                    f"cell {cell}: {key} rust={rv} R={av} -- the arms disagree"
                )
            else:
                print(f"  {key}: {rv} on both arms")

    print()
    if failures:
        print("ARMS DISAGREE -- they are not running the same analysis:")
        for f in failures:
            print(f"  * {f}")
        print()
        print("Every parameter both arms use must come from the grid file. A value")
        print("written down twice, in two languages, is a value that will drift.")
        return 1
    print(f"arms agree on {args.cells} (n={args.reps} reps each)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
