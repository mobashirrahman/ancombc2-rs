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
| R 4.6 | `rocker/r-ver:4.6` | the current release, and what all five of CI's R jobs run. The oracle declares `R (>= 4.5.0)`; the R version turns out not to be what makes a number reproducible — see below |
| ANCOMBC | `git checkout --detach dc4febdf…` **plus** `test "$(git rev-parse HEAD)" = "$ORACLE_SHA"` | verified *in the build*, so an image with the wrong commit fails to build instead of quietly benchmarking a different oracle |
| Rust | `rust-toolchain.toml` + `Cargo.lock` | the same compiler the lockfile expects |
| `nloptr` | installed from CRAN | `NLOPT_LN_NELDERMEAD` in `.bias_em`; MASS ships with R |

BLAS/LAPACK are **reported, not chosen**. The base image's are used and written
into `benchmarks/results/container_sessioninfo.txt`. Swapping BLAS to make a gate
pass would be choosing the ruler, so the image deliberately does not do it; what
it does is make the ruler identical for both arms.

That "reported, not chosen" stance was written before the BLAS turned out to be
the one variable that *does* move numbers. It does: the committed goldens were
generated on Ubuntu's reference BLAS, and an OpenBLAS R shifts `beta` by up to
2.52e-06 relative — above the golden contract's `rtol 1e-8` — while every `se`
stays within 1e-11. `docs/reproduction.md` has the measurements. So the report is
load-bearing, and `environment-linux-64.lock` names the BLAS outright rather than
leaving it to be inferred.

## Running it without a container engine

`environment-linux-64.lock` is an exact `micromamba --explicit` pin for the same
environment: R 4.5.3, `libopenblas-0.3.34-pthreads`, and all 15 of the oracle's
`Imports` as binaries. It exists because the host these results were produced on
cannot run the image at all (below), and because a lock says *which* BLAS it is,
which an image digest does not.

```sh
micromamba create -p ./bench-env \
  --file benchmarks/container/environment-linux-64.lock

export R_LIBS_USER="$PWD/bench-env"
export PATH="$PWD/bench-env/bin:$PATH"
R CMD INSTALL --no-docs -l "$R_LIBS_USER" reference/ANCOMBC
ANCOMBC_ORACLE_DIR=reference/ANCOMBC \
  python3 scripts/run_benchmarks.py --data benchmarks/datasets
make gates
```

This is verified rather than aspirational: all 15 `Imports` resolve from the lock
alone, ANCOMBC installs and runs in it, and `scripts/verify_real_package.R` runs
in it.

R is pinned at 4.5.3 rather than 4.6 because that is the newest version for which
every `Imports` entry exists as a conda-forge binary — on 4.6 the solver has no
`r-rdpack` build, so `Rdpack`, `Hmisc` and `DescTools` have to come from CRAN
source and the lock stops pinning everything. Both versions are verified
bit-identical to 4.3.3 on one BLAS, 915 of 915 arrays, so the difference moves no
number.

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

**Not this, and not because a binary is missing.** The environment these results
were produced in cannot run a container at all. `benchmarks/container/run.sh`
exits 2 with an explanation and a working alternative rather than pretending.

The reason is worth stating precisely, because "docker isn't installed" is the
wrong diagnosis and would send someone off to install one:

```
$ unshare -Ur true
unshare: write failed /proc/self/uid_map: Operation not permitted
```

User namespaces are denied by the host, which is what every rootless engine
needs. A static `podman` 5.4.0 was fetched, and its bundled `conmon` and `crun`
both execute — but podman still cannot `reexec` into a user namespace. This is not
fixable by downloading anything: on this machine the container path is closed, and
`environment-linux-64.lock` above is how to get the same guarantee anyway.

Every number in `benchmarks/results/results.jsonl` came from
`python3 scripts/run_benchmarks.py --data benchmarks/datasets` on the host, with
the oracle sourced from `reference/ANCOMBC` and R 4.3.3. That is the same
harness, the same arms and the same results file; what it is *not* is the pinned R
the image would provide. `docs/compatibility.md` records this, and
`docs/reproduction.md` gives both commands.

To execute the container path, run it on a host with Docker or Podman and copy
`benchmarks/results/` back. `run.sh` prefers `docker` and falls back to `podman`,
and honours `ANCOMBC2_BENCH_IMAGE` and `ORACLE_SHA`.
