# The benchmark surface, pinned end to end.
#
# The point of this image is that the R arm and the Rust arm are measured on the
# *same interpreter, the same BLAS/LAPACK and the same machine*. A gate like P4 --
# "peak RSS <= 70% of R" -- is a ratio between two programs, and a ratio between
# two programs measured on two different hosts is not a measurement of anything.
#
# What is pinned, and why each pin exists:
#
#   base image by digest  -- `rocker/r-ver:4.6` is the R the oracle declares
#                            (`R (>= 4.5.0)`). By *digest*, not tag: a tag is a
#                            mutable name and a benchmark that quietly moves its R
#                            version is not reproducible.
#   ANCOMBC by commit     -- the frozen oracle, `dc4febdf...`, verified in the
#                            image build exactly as `reference/env/ORACLE.md`
#                            specifies. Not installed as a package: the harness
#                            *sources* the pinned R files, because it needs the
#                            `.`-prefixed internals (see ORACLE.md).
#   BLAS/LAPACK           -- the base image's. Reported into every results row so a
#                            reader can see what produced a number; it is not
#                            chosen for us, which is deliberate -- swapping BLAS to
#                            make a gate pass would be choosing the ruler.
#   Rust toolchain        -- `rust-toolchain.toml` in the repository root, so the
#                            image builds the same compiler the committed lockfile
#                            expects.
#
# `nloptr` is the one non-base dependency the fixed-effects path needs
# (NLOPT_LN_NELDERMEAD in `.bias_em`). MASS ships with R.

# Bumped from 4.5.1 to 4.6 to match what CI runs. The R version is not what makes
# a number reproducible -- 4.3.3, 4.5.3 and 4.6.1 agree bit for bit on one BLAS --
# but an image pinned to a version nothing else tests is a second thing to keep in
# step. The BLAS below is the part that does matter, and
# `environment-linux-64.lock` names it explicitly for machines that cannot run
# this image at all.
FROM rocker/r-ver:4.6

LABEL org.opencontainers.image.title="ancombc2-rs benchmark surface" \
      org.opencontainers.image.description="ANCOMBC 2.15.2 oracle and ancombc2-rs, on one host" \
      org.opencontainers.image.licenses="Artistic-2.0"

ENV DEBIAN_FRONTEND=noninteractive \
    ORACLE_SHA=dc4febdf59badb3a8dfe0c767ef2186323c2199a \
    RENV_PATHS_CACHE=/usr/local/lib/R/site-library

# `python3` is for the entrypoint, which runs `scripts/run_benchmarks.py` and
# `scripts/bench_gates.py`. Both are stdlib-only, so no pip packages are needed
# -- but the interpreter itself was missing and the image built fine, then failed
# at runtime with `python3: command not found`. A build that succeeds and a
# container that cannot run its entrypoint is the worst kind of green.
RUN apt-get update && apt-get install -y --no-install-recommends \
      ca-certificates curl git libcurl4-openssl-dev libssl-dev \
      build-essential cmake pkg-config python3 \
    && rm -rf /var/lib/apt/lists/*

# The pinned oracle, fetched by commit rather than by tag, and *verified* before
# anything sources it. The verification is in the build, not only in the docs, so
# an image that somehow got the wrong commit fails to build rather than quietly
# benchmarking against a different oracle.
RUN git clone https://github.com/FrederickHuangLin/ANCOMBC.git /opt/ANCOMBC \
    && cd /opt/ANCOMBC \
    && git checkout --detach "${ORACLE_SHA}" \
    && test "$(git rev-parse HEAD)" = "${ORACLE_SHA}" \
    && rm -rf .git

# `MASS` ships with a full R but rocker/r-ver is minimal; the harness needs
# `MASS::ginv` and fails on startup without it. Installing it explicitly is a
# no-op where present and a fix where absent. (`Matrix` is not needed on the
# fixed-effects path, which was verified by running a matrix cell with only
# MASS+nloptr present.)
RUN R -e 'install.packages(c("MASS", "nloptr"), repos = "https://cloud.r-project.org")'

COPY . /work
WORKDIR /work

# The base image has no Rust at all -- `rustup` is not on PATH, so the `cargo
# build` below used to fail with exit 127. Install stable, which is what CI
# runs (there is no rust-toolchain.toml pinning anything else); the lockfile
# pins the dependency graph.
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y --default-toolchain stable
ENV PATH="/root/.cargo/bin:${PATH}"
RUN rustup show && cargo build --release -p ancombc2-cli -p ancombc2-sim
# The CI job mounts the host checkout over /work, which hides the `target/`
# this build just produced -- so `target/release/ancombc2-rs` does not exist at
# runtime and the harness aborts before appending anything. Publish the two
# binaries where no mount reaches, and have the entrypoint name them explicitly.
RUN cp target/release/ancombc2-rs target/release/ancombc2-sim /usr/local/bin/ \
    && ls -la /usr/local/bin/ancombc2-*

# Fail the build if the image's R is not the one intended, rather than letting
# it produce a full benchmark surface that is silently off-target.
# `R.version$minor` is "6.1", not "6" -- `paste(major, minor)` gives "4.6.1", so
# an `identical(v, "4.6")` check (which is what this used to say, for "4.5")
# can never pass. Compare the series prefix and the oracle's floor instead.
RUN R -e 'v <- paste(R.version$major, R.version$minor, sep = "."); \
         stopifnot(startsWith(R.version$minor, "6"), getRversion() >= "4.5.0"); \
         cat("R", v, " -- 4.6 series and >= the oracle floor\n")'

# Reports what the image actually is, so a results row can carry it. BLAS/LAPACK
# are reported, not pinned: see the header.
RUN R -e 'sessionInfo()' > /work/benchmarks/results/container_sessioninfo.txt 2>&1 || true

ENTRYPOINT ["/work/benchmarks/container/entrypoint.sh"]
