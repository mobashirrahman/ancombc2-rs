# Containerised benchmark harness

The plan's benchmarking item asks for a "committed container-based harness
appending to `benchmarks/results/results.jsonl`". That is these three files.

## Why a container at all

P4 is a *ratio*: "peak RSS ≤ 60–70% of R". A ratio between two programs is only a
measurement of the two programs if they ran on the same machine, with the same
interpreter and the same BLAS. Run R on one host and Rust on another and the
number still looks like a ratio, and is worse than useless — it is a ratio with a
host difference folded into it. The image pins:

| pinned | how | why |
| --- | --- | --- |
| R 4.5.1 | `rocker/r-ver:4.5.1` | the version the oracle's `DESCRIPTION` declares (`R (>= 4.5.0)`); the goldens in this repo were generated on 4.3.3, one series below |
| ANCOMBC | `git checkout --detach dc4febdf…` **plus** `test "$(git rev-parse HEAD)" = "$ORACLE_SHA"` | verified *in the build*, so an image with the wrong commit fails to build instead of quietly benchmarking a different oracle |
| Rust | `rust-toolchain.toml` + `Cargo.lock` | the same compiler the lockfile expects |
| `nloptr` | installed from CRAN | `NLOPT_LN_NELDERMEAD` in `.bias_em`; MASS ships with R |

BLAS/LAPACK are **reported, not chosen**. The base image's are used and written
into `benchmarks/results/container_sessioninfo.txt`. Swapping BLAS to make a gate
pass would be choosing the ruler, so the image deliberately does not do it; what
it does is make the ruler identical for both arms.

## Running it

```
benchmarks/container/run.sh                    # the whole surface
benchmarks/container/run.sh --datasets bm4,bm5
```

Results are bind-mounted, so the image appends to the repository's own
`benchmarks/results/results.jsonl`. `results.jsonl` is append-only and
`scripts/bench_gates.py` keeps the last row per `(dataset, arm)`, so a re-run
supersedes rather than overwrites. The gates are evaluated **inside** the
container, so each ratio comes from rows that run wrote.

## What has actually been run

**Not this.** The environment these results were produced in has no container
engine — neither `docker` nor `podman` is on `PATH` — so `benchmarks/container/run.sh`
exits 2 with an explanation and a working alternative rather than pretending.

Every number in `benchmarks/results/results.jsonl` came from
`python3 scripts/run_benchmarks.py --data benchmarks/datasets` on the host, with
the oracle sourced from `reference/ANCOMBC` and R 4.3.3. That is the same
harness, the same arms and the same results file; what it is *not* is the pinned
R 4.5.1 the image would provide. `docs/compatibility.md` records this, and
`docs/reproduction.md` gives both commands.

To execute the container path, run it on a host with Docker or Podman and copy
`benchmarks/results/` back. `run.sh` prefers `docker` and falls back to `podman`,
and honours `ANCOMBC2_BENCH_IMAGE` and `ORACLE_SHA`.
