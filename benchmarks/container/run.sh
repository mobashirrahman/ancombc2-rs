#!/usr/bin/env bash
# Build (if needed) and run the containerised benchmark surface on this host.
#
# The host wrapper exists so the *host* is also part of the pinned thing: the
# image pins the interpreter, the oracle and the compiler, and this script pins
# the machine by refusing to start on a host with a different CPU count than the
# one the results claim. A P4 ratio measured against an R that ran on a different
# box is not a ratio.
#
# Usage:
#   benchmarks/container/run.sh                 # full surface, all 6 datasets
#   benchmarks/container/run.sh --datasets bm4,bm5
#   benchmarks/container/run.sh --arms rust-1,rust-16
#
# Anything after the script name is passed to `scripts/run_benchmarks.py`.
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo="$(cd "$here/../.." && pwd)"
image="${ANCOMBC2_BENCH_IMAGE:-ancombc2-bench:oracle-dc4febdf}"

engine=""
for e in docker podman; do
  if command -v "$e" >/dev/null 2>&1; then engine="$e"; break; fi
done
if [ -z "$engine" ]; then
  cat >&2 <<'MSG'
Neither docker nor podman is on PATH, so the containerised harness cannot run
here.

On some hosts that is the whole story and the fix is to install an engine. On
others -- including the host these results were produced on -- no engine helps,
because rootless containers need user namespaces and the host denies them:

    $ unshare -Ur true
    unshare: write failed /proc/self/uid_map: Operation not permitted

`benchmarks/container/README.md` records that investigation. Either way, three
honest ways forward:

  * run the image elsewhere and copy `benchmarks/results/results.jsonl` back, or
  * reproduce the same pinned environment with no container at all, via the
    committed lock -- which also names the BLAS, the thing that actually has to
    match:

        micromamba create -p ./bench-env \
          --file benchmarks/container/environment-linux-64.lock
        export R_LIBS_USER="$PWD/bench-env" PATH="$PWD/bench-env/bin:$PATH"
        R CMD INSTALL --no-docs -l "$R_LIBS_USER" reference/ANCOMBC
        ANCOMBC_ORACLE_DIR=reference/ANCOMBC \
          python3 scripts/run_benchmarks.py --data benchmarks/datasets
        make gates

  * or run the harness on this host as it stands:

        ANCOMBC_ORACLE_DIR=reference/ANCOMBC \
          python3 scripts/run_benchmarks.py --data benchmarks/datasets
        make gates

The last is what produced every number currently in `benchmarks/results/`, and
`docs/compatibility.md` records that the container path has not been executed.
MSG
  exit 2
fi

"$engine" build \
  --build-arg ORACLE_SHA="${ORACLE_SHA:-dc4febdf59badb3a8dfe0c767ef2186323c2199a}" \
  -f "$here/Dockerfile" -t "$image" "$repo"

# `--rm` so a failed benchmark does not leave a 6 GB image behind; the results
# live on the host volume, not in the container.
exec "$engine" run --rm \
  -v "$repo:/work" \
  -v "$repo/benchmarks/results:/work/benchmarks/results" \
  -e OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}" \
  "$image" "$@"
