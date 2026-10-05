#!/usr/bin/env bash
# Run the whole benchmark surface inside the pinned image and append to the
# committed results file.
#
# The results directory is bind-mounted from the host, so the image *writes the
# repository's* `benchmarks/results/results.jsonl` rather than a copy that has to
# be merged by hand. Merging by hand is how a benchmark result goes missing, and
# `results.jsonl` is append-only on purpose: `scripts/bench_gates.py` keeps the
# last row per (dataset, arm), so a re-run supersedes rather than overwrites.
set -euo pipefail

cd /work

# The oracle is sourced from /opt/ANCOMBC. `ANCOMBC_ORACLE_DIR` is what
# `scripts/sim_r.R` and `scripts/bench_r.R` read; the vendored copy at
# `reference/ANCOMBC` is kept for the non-container path.
export ANCOMBC_ORACLE_DIR="${ANCOMBC_ORACLE_DIR:-/opt/ANCOMBC}"

# Number of R threads for the `r-parallel` arm. 0 means "use every core", which is
# what "R fully parallel" means in the plan's seven-arm list.
: "${ANCOMBC_BENCH_THREADS:=0}"
export ANCOMBC_BENCH_THREADS

# Pin BLAS threading *before* R starts, and record what it was set to. A BLAS that
# silently spawns a thread per core inside an already-parallel R arm is a
# benchmark that measures oversubscription.
: "${OMP_NUM_THREADS:=1}"
export OMP_NUM_THREADS

echo "=== ancombc2-rs benchmark surface (containerised) ==="
R --version | head -1
echo "oracle: $(git -C "$ANCOMBC_ORACLE_DIR" rev-parse HEAD 2>/dev/null || echo "vendored, not a git checkout")"
echo "rustc:  $(rustc --version)"
echo "cpus:   $(nproc)   OMP_NUM_THREADS=$OMP_NUM_THREADS"

# Both arms, then the gates. The gates are evaluated *inside* the container so the
# ratios come from rows this run wrote; evaluating them outside would compare rows
# recorded on different hosts.
# `--binary` names the image's own build explicitly. The default,
# `target/release/ancombc2-rs`, resolves under /work, which the CI job mounts
# over with the host checkout -- hiding the binaries the image built. Without
# this the harness reports "not found" and exits 2 before appending anything.
python3 scripts/run_benchmarks.py --data benchmarks/datasets \
    --binary /usr/local/bin/ancombc2-rs "$@"
python3 scripts/bench_gates.py --results benchmarks/results/results.jsonl \
    --report benchmarks/results/gates.json

echo "=== wrote benchmarks/results/results.jsonl and gates.json ==="
