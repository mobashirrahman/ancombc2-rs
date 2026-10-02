# ancombc2-rs
#
# Parity first, speed second. Every target here is a shortcut for something a
# developer does often; nothing here hides a failing test.
#
#   make test          unit + property + small parity (what CI runs on a push)
#   make parity        the small golden fixtures, verbosely
#   make parity-large  fx04: 9,800 x 500 with 50 sensitivity refits
#   make goldens       regenerate the goldens from the oracle (needs R)
#   make bench-data    generate the benchmark datasets
#   make bench         run the benchmark surface and evaluate the gates
#   make gates         evaluate the gates against the committed results
#   make sim           the quick simulation grid, both arms, pooled and compared
#   make sim-full      the full simulation grid (hours; the weekly CI job)
#   make realdata      Layer 4: the real datasets, both arms, compared
#   make realdata-prep re-derive validation/realdata from their sources
#   make phyloseq-shim the S4 class definitions the realdata prep needs
#   make edge-cases   the structural-zero / edge-case matrix (PLAN.md 5.5)
#   make edge-gen     regenerate those cases from the oracle
#   make fmt clippy    the CI lints
#   make clean         build artefacts only; `distclean` also removes data
#
# `make goldens` needs the pinned ANCOMBC checkout and an R with `nloptr`; see
# reference/env/ORACLE.md. Everything else needs only a Rust toolchain.

CARGO ?= cargo
PYTHON ?= python3
RSCRIPT ?= Rscript
RESULTS := benchmarks/results/results.jsonl

# The oracle needs `nloptr`, which is installed into a private library because
# the system R does not have it. Both R arms pass this through; without it the
# oracle fails to load, and the harnesses record that failure in every row rather
# than producing numbers, so a run that quietly analysed nothing still looks like
# a run. See docs/reproduction.md for the install command.
# This machine's private library is `/scratch/mdra00001/rlib`; override on the
# command line or in the environment for a different checkout.
R_LIBS_USER ?= /scratch/mdra00001/rlib

.PHONY: all build test parity parity-large properties fmt clippy lint doc \
        goldens goldens-drift smoke bench-data bench bench-container gates sim \
        sim-rust sim-r sim-agree sim-summary \
        sim-full realdata realdata-prep realdata-rust realdata-r realdata-summary \
        phyloseq-shim edge-cases edge-gen \
        miri r-install r-test distclean clean

all: lint test

build:
	$(CARGO) build --workspace --release

# --- tests -----------------------------------------------------------------

test:
	$(CARGO) test --workspace

parity:
	$(CARGO) test -p ancombc2-core --test parity -- --nocapture

# The large fixture is `#[ignore]`d so a push stays fast. It is a required gate,
# not an optional extra: the two-taxon-set behaviour only shows up on a fixture
# with structural zeros, and that is exactly what fx04 has.
parity-large:
	$(CARGO) test --release -p ancombc2-core --test parity -- --ignored --nocapture

properties:
	$(CARGO) test -p ancombc2-core --test properties -- --nocapture

# --- lints -----------------------------------------------------------------

fmt:
	$(CARGO) fmt --all -- --check

clippy:
	$(CARGO) clippy --workspace --all-targets -- -D warnings

lint: fmt clippy

# A miri smoke test catches the undefined behaviour a numerics crate is prone to:
# out-of-bounds slices in the QR and sandwich accumulation, and unaligned reads in
# the golden reader. It is slow, so only the smallest test binary runs it.
miri:
	$(CARGO) +nightly miri test -p ancombc2-stats
	$(CARGO) +nightly miri test -p ancombc2-core --lib matrix::

# --- oracle ----------------------------------------------------------------

# Regenerate the goldens. The generator verifies the oracle commit and runs a
# mirror self-check before writing, so a drifted oracle fails here rather than
# silently rewriting the contract.
goldens:
	$(RSCRIPT) scripts/smoke_oracle.R
	$(RSCRIPT) scripts/generate_goldens.R 1 2 3 4

smoke:
	$(RSCRIPT) scripts/smoke_oracle.R

# Verify the committed goldens still are what the pinned oracle produces, without
# touching the working tree.
#
# `fx01`..`fx04` are checked by **recomputing the contract from the committed
# fixtures** and comparing at a stated tolerance, not by regeneration. Two
# independent reasons, both measured:
#
#   * the fixtures predate the `set.seed(spec$seed)` in `reference/R/fixtures.R`
#     and were produced from the ambient session state, so regenerating them
#     produces different *counts* -- a different input, not a drifted golden;
#   * even taking the committed fixture files as given, `write.table` renders a
#     double at 15 significant digits, so the goldens cannot be reproduced
#     bit-for-bit. Recomputation agrees to 0 (fx01), 4.1e-12 (fx02), 1.8e-11
#     (fx03) and 0 (fx04) relative -- the text round-trip floor, three orders
#     below the Level B `rtol 1e-8` the parity suite enforces.
#
#
# The edge-case and fixture-matrix goldens *are* generated reproducibly from
# fixed seeds, so those two are compared byte-for-byte and excluding the `.rds`
# wrappers (which carry a `session` recording and a reference-compared `capture`
# environment -- see `scripts/compare_goldens.R`).
#
# An earlier version of this target regenerated in place. It destroyed the
# fixtures, reported 82 changed files that were a different input rather than
# drift, and left the repository failing parity until restored from the index.
# Each fixture is compared at its own measured floor -- 0, 4.1e-12, 1.8e-11 and
# 9.6e-04 respectively -- because `fx04` amplifies the round-off through an MLE over
# 9,800 taxa and is the only one running the 50-refit path. The oracle itself is
# deterministic: two runs over the same text agree to 0.
# `scripts/compare_goldens.R` prints each observed deviation, and passing
# `GOLDEN_TOL=--tol <x>` overrides the floors for all four at once.
GOLDEN_TOL ?=
DRIFT_TMP  ?= /tmp/ancombc2-drift
goldens-drift: oracle-sha
	@rm -rf $(DRIFT_TMP) && mkdir -p $(DRIFT_TMP)
	@cp -r validation/golden $(DRIFT_TMP)/committed-golden
	@cp -r validation/edge   $(DRIFT_TMP)/committed-edge
	@cp -r validation/matrix $(DRIFT_TMP)/committed-matrix
	R_LIBS_USER=$(R_LIBS_USER) ANCOMBC_ORACLE_DIR=reference/ANCOMBC \
	  $(RSCRIPT) --vanilla scripts/generate_goldens.R --from-committed 1 2 3 4
	@$(RSCRIPT) --vanilla scripts/compare_goldens.R \
	  $(DRIFT_TMP)/committed-golden validation/golden $(GOLDEN_TOL)
	@$(MAKE) --no-print-directory edge-gen matrix-gen
	@# `*.rds` because R's serialisation embeds the session's RNG state and the
	@# environment captured alongside the goldens, neither of which is a compared
	@# quantity. `stage_seconds.json` because it is wall-clock: it necessarily
	@# differs between two runs, which is why the four main fixtures are compared
	@# through `scripts/compare_goldens.R`, which checks the stage names and their
	@# being finite rather than the seconds themselves. Every other file in these
	@# trees -- including `convergence_trace.json` and `em_mixture.json` -- is
	@# compared byte for byte, since all of it is deterministic.
	@diff -rq --exclude='*.rds' --exclude='stage_seconds.json' \
	  $(DRIFT_TMP)/committed-edge validation/edge \
	  || { echo "ERROR: the edge-case goldens have drifted"; exit 1; }
	@diff -rq --exclude='*.rds' --exclude='stage_seconds.json' \
	  $(DRIFT_TMP)/committed-matrix validation/matrix \
	  || { echo "ERROR: the fixture-matrix goldens have drifted"; exit 1; }
	@echo "the committed goldens match the pinned oracle"

# --- benchmarks ------------------------------------------------------------

bench-data:
	$(PYTHON) scripts/make_bench_data.py --out benchmarks/datasets

# The R arms need the oracle; without it they record `skipped` and the gate
# report says `not measured` for the R-relative gates, which is the honest
# outcome rather than a pass.
bench: build bench-data
	R_LIBS_USER=$(R_LIBS_USER) ANCOMBC_ORACLE_DIR=reference/ANCOMBC \
	  $(PYTHON) scripts/run_benchmarks.py \
	    --data benchmarks/datasets --results $(RESULTS)
	$(MAKE) gates

gates:
	$(PYTHON) scripts/bench_gates.py --results $(RESULTS) \
	    --report benchmarks/results/gates.json

# The same surface, in a pinned image: one interpreter, one BLAS, one host for
# both arms. Every number currently in benchmarks/results/ was produced by `make
# bench` on the host rather than by this, because no container engine was
# available; see docs/compatibility.md, which says so rather than implying the
# container path was exercised.
bench-container:
	benchmarks/container/run.sh

# --- simulation ------------------------------------------------------------
#
# Three steps, because the two arms must analyse the same tables: the Rust
# generator writes them, the R arm reads them back, and only then are the two
# pooled. `sim` needs R and the pinned oracle for the R arm; without them the
# `sim-rust` target alone still reports the Rust arm's calibration, and
# `sim-summary` says plainly that no R arm was supplied.

SIM_GRID ?= validation/simulation/quick
SIM_DATA ?= sim-data
SIM_OUT  ?= validation/simulation/results
SIM_REPS ?=

sim-rust: build
	$(CARGO) build --release -p ancombc2-sim
	./target/release/ancombc2-sim generate --grid $(SIM_GRID) --out $(SIM_DATA)
	./target/release/ancombc2-sim rust --grid $(SIM_GRID) \
	    --out $(SIM_OUT)/quick.rust.jsonl --progress

# Two arms, one analysis. Compares the counts that the filter and inference
# configuration decide, on a few cells, and fails if they differ. Exists because
# the two arms once ran different analyses -- struc_zero/neg_lb hardcoded TRUE in
# sim_r.R and left at FALSE by the Rust config -- and nothing noticed until the
# resulting numbers were read as a result.
sim-agree:
	R_LIBS_USER=$(R_LIBS_USER) ANCOMBC_ORACLE_DIR=reference/ANCOMBC \
	  $(PYTHON) scripts/check_sim_arms_agree.py \
	    --grid $(SIM_GRID)/grid.json --data $(SIM_DATA)

sim-r:
	R_LIBS_USER=$(R_LIBS_USER) ANCOMBC_ORACLE_DIR=reference/ANCOMBC \
	  $(RSCRIPT) --vanilla scripts/sim_r.R \
	    --data $(SIM_DATA) --grid $(SIM_GRID)/grid.json \
	    --out $(SIM_OUT)/quick.r.jsonl

sim-summary:
	./target/release/ancombc2-sim summarise --grid $(SIM_GRID) \
	    --rust $(SIM_OUT)/quick.rust.jsonl --r $(SIM_OUT)/quick.r.jsonl \
	    --out $(SIM_OUT)/quick.summary.json

sim: sim-rust sim-r sim-summary

sim-full:
	$(MAKE) sim SIM_GRID=validation/simulation/full

# --- the edge-case matrix (PLAN.md section 5.5) ----------------------------
#
# Seven small cases, each a fixture plus a golden captured from the pinned
# oracle. `edge-gen` needs the oracle; `edge-cases` only reads the committed
# goldens, so it is a normal test target and runs in CI on every push.

edge-gen:
	R_LIBS_USER=$(R_LIBS_USER) ANCOMBC_ORACLE_DIR=reference/ANCOMBC \
	  $(RSCRIPT) --vanilla scripts/generate_edge_cases.R

edge-cases:
	$(CARGO) test -p ancombc2-core --test edge_cases

# --- the golden fixture matrix (PLAN.md section 5.5) ------------------------
#
# The four committed fixtures are a handful of shapes; the matrix is a sweep over
# the ten axes PLAN.md lists -- shape, sparsity, predictor type, covariate count
# and interactions, group balance, structural zeros, pseudo-count, sensitivity,
# conservative, and the adjustment method. Cells are declared in
# `reference/R/fixture_matrix.R` and generated by
# `scripts/generate_matrix_goldens.R`.
#
# `matrix-gen` needs the pinned oracle. `matrix` only reads the committed
# goldens, so it is a normal test target and runs in CI on every push -- including
# the four structural checks that assert the matrix still covers every axis level
# the plan asks for, which is the part that would otherwise rot silently.
#
# `matrix-drift` regenerates and diffs. It is the weekly CI job, and it needs the
# declared R 4.5: the committed fx01..fx04 predate a seeding fix in
# `fixtures.R` and are deliberately not regenerated, so `goldens-drift` covers
# those and `matrix-drift` covers these.

MATRIX_CELLS ?=

matrix-gen:
	R_LIBS_USER=$(R_LIBS_USER) ANCOMBC_ORACLE_DIR=reference/ANCOMBC \
	  $(RSCRIPT) --vanilla scripts/generate_matrix_goldens.R $(MATRIX_CELLS)

matrix:
	$(CARGO) test --release -p ancombc2-core --test fixture_matrix

# Runs every matrix cell through the *installed* ANCOMBC package and compares the
# result against the committed goldens, holding nothing but the interpreter
# constant. It is not a CI gate and not `matrix`: it needs an R >= 4.5.0 (the
# oracle's own `Depends`) with ANCOMBC and its Imports installed, which is a
# heavier setup than the harness in `reference/R` deliberately avoids. What it
# buys is the one thing the harness cannot check about itself -- that sourcing the
# oracle's R files reproduces what the installed package actually does.
#
#   make verify-real-package
#
# `GOLDEN_ROOT` retargets the comparison. Pointing it at a tree generated by the
# harness on the *same* interpreter is the comparison that isolates harness from
# package; the default, the committed tree, additionally measures the difference
# between this machine's BLAS and the one the goldens were generated with.
verify-real-package:
	R_LIBS_USER=$(R_LIBS_USER) GOLDEN_ROOT=$(GOLDEN_ROOT) ANCOMBC_REPO=$(CURDIR) \
	  $(RSCRIPT) --vanilla scripts/verify_real_package.R

matrix-drift: matrix-gen
	@git diff --stat validation/matrix || true
	@git diff --exit-code -- validation/matrix \
	  || (echo "the fixture matrix goldens drifted -- the oracle moved, or" \
	      "scripts/generate_matrix_goldens.R is not reproducible" && false)

# --- real-data validation (Layer 4) ---------------------------------------
#
# `realdata-prep` needs the pinned oracle for the vignette reconstruction and,
# for the two datasets that ship inside a phyloseq object, the class
# definitions described in docs/reproduction.md. It writes the datasets; the
# comparison never regenerates them, because both arms must read the same bytes.
#
# `realdata` fails the build if any dataset misses the 99.99% diff_abn
# threshold. That is the plan's Level D real-data gate, and the two arms'
# agreement on which taxa diverge is the thing being asserted.

REALDATA ?= validation/realdata
REALDATA_MB ?= microbiome/data

# The two phyloseq-sourced datasets need the S4 class definitions, which this
# machine cannot get from Bioconductor (RCurl needs libcurl headers, and there is
# no root). See docs/reproduction.md.
# The oracle pin. `reference/ANCOMBC` is committed as ordinary files -- the
# harness sources `reference/ANCOMBC/R/*.R`, so a submodule or gitlink would leave
# a clone unable to check a single golden -- which means there is no `.git` in the
# tree to ask for its own commit. This checks the recorded SHA-256 of every file
# the harness reads, and cross-checks the id against `ORACLE_SHA` in
# `crates/ancombc2-core/src/lib.rs`. Cheap, and it is the difference between a
# frozen oracle and a directory that happens to look right.
oracle-sha:
	$(RSCRIPT) --vanilla scripts/check_oracle_sha.R

phyloseq-shim:
	R_LIBS_USER=$(R_LIBS_USER) $(RSCRIPT) --vanilla scripts/install_phyloseq_shim.R \
	    --lib $(R_LIBS_USER)

realdata-prep: phyloseq-shim
	R_LIBS_USER=$(R_LIBS_USER) ANCOMBC_ORACLE_DIR=reference/ANCOMBC \
	  $(RSCRIPT) --vanilla scripts/prepare_realdata.R \
	    --out $(REALDATA) --microbiome-data $(REALDATA_MB)

realdata-rust: build
	$(CARGO) build --release -p ancombc2-sim
	./target/release/ancombc2-sim realdata --data $(REALDATA) \
	    --out $(REALDATA)/results/rust.jsonl

realdata-r:
	@mkdir -p $(REALDATA)/results
	R_LIBS_USER=$(R_LIBS_USER) ANCOMBC_ORACLE_DIR=reference/ANCOMBC \
	  $(RSCRIPT) --vanilla scripts/realdata_r.R \
	    --data $(REALDATA) --out $(REALDATA)/results/r.jsonl

realdata-summary:
	$(PYTHON) scripts/compare_realdata.py \
	    --rust $(REALDATA)/results/rust.jsonl \
	    --r $(REALDATA)/results/r.jsonl \
	    --out $(REALDATA)/results/summary.json

# The `y = x` scatters PLAN.md section 17 asks for. Committed rather than generated
# on demand, because the point of "scatter visually on y = x" as a review step is
# that a reviewer can look at the figure that was accepted.
#
# matplotlib is a build-time-only dependency and is deliberately not required for the
# numeric comparison: if it is missing, the figures are skipped with a message and
# the gate still runs, rather than the real-data layer failing to reproduce.
realdata-scatter: realdata-summary
	@$(PYTHON) -c "import matplotlib" 2>/dev/null \
	  && $(PYTHON) scripts/plot_realdata.py \
	       --data $(REALDATA)/results --out $(REALDATA)/results/scatter \
	  || echo "NOTE: matplotlib unavailable; skipping the real-data y=x scatters"

realdata: realdata-rust realdata-r realdata-summary realdata-scatter
	$(PYTHON) scripts/compare_realdata.py \
	    --rust $(REALDATA)/results/rust.jsonl \
	    --r $(REALDATA)/results/r.jsonl --fail

# --- the R wrapper ---------------------------------------------------------
#
# `r-install` builds the Rust core as a cdylib and links it, which is what
# `src/Makevars` does; `r-test` then runs the testthat suite. Both need a writable
# R library, hence ANCOMBC2RS_RLIB.

ANCOMBC2RS_RLIB ?= $(CURDIR)/.rlib

# `.Rprofile` puts the oracle's private library first on `.libPaths()`, because that
# is where nloptr lives. An `ancombc2rs` installed *there* would therefore shadow the
# one this target builds, and `r-test` would quietly exercise a stale build -- which
# it did, for a month: the 49 green tests were green against a September copy.
#
# So the oracle library is checked for a stray `ancombc2rs` and removed rather than
# left to shadow, and the test step asserts which copy it loaded.
r-install:
	mkdir -p $(ANCOMBC2RS_RLIB)
	R_LIBS_USER=$(ANCOMBC2RS_RLIB) $(RSCRIPT) -e 'dir.create("$(ANCOMBC2RS_RLIB)", showWarnings = FALSE)' \
	  -e 'install.packages(c("jsonlite", "testthat"), lib = "$(ANCOMBC2RS_RLIB)", repos = "https://cloud.r-project.org", quiet = TRUE)' \
	  -e 'if (!requireNamespace("testthat", quietly = TRUE)) stop("testthat is required for r-test")'
	R_LIBS_USER=$(ANCOMBC2RS_RLIB) $(RSCRIPT) -e \
	  'install.packages("r/ancombc2rs", repos = NULL, type = "source", lib = "$(ANCOMBC2RS_RLIB)")'
	@R_LIBS_USER=$(ANCOMBC2RS_RLIB) $(RSCRIPT) --vanilla scripts/check_r_lib.R $(ANCOMBC2RS_RLIB)
r-test: r-install
	R_LIBS_USER=$(ANCOMBC2RS_RLIB) $(RSCRIPT) -e \
	  'library(testthat); library(ancombc2rs); testthat::test_dir("r/ancombc2rs/tests/testthat", package = "ancombc2rs", stop_on_failure = TRUE)'
# `cargo doc` is a CI gate: a broken intra-doc link is a documentation failure,
# and the compatibility statement is referenced from the crate docs.
doc:
	RUSTDOCFLAGS="-D warnings" $(CARGO) doc --workspace --no-deps

# --- housekeeping ----------------------------------------------------------

clean:
	$(CARGO) clean
	rm -rf benchmarks/results/run

distclean: clean
	rm -f $(RESULTS) benchmarks/results/gates.json
	rm -rf benchmarks/datasets $(ANCOMBC2RS_RLIB) sim-data
	rm -f r/ancombc2rs/src/libancombc2_ffi.so r/ancombc2rs/src/.ancombc2_rs_built
