# Reproduction

Everything here is a command. Nothing in the repository depends on state that is
not either committed or produced by one of them.

## Prerequisites

| need | for | how to check |
| --- | --- | --- |
| Rust 1.80 or later | building and testing | `cargo --version` |
| R 4.3 or later | regenerating the goldens, the R benchmark arm, the profile | `Rscript --version` |
| `nloptr` | the oracle harness | `Rscript -e 'library(nloptr)'` |
| `MASS` | the oracle harness (`MASS::ginv`) | `Rscript -e 'library(MASS)'` |
| Python 3.9 or later | the benchmark, gate and simulation scripts | `python3 --version` |
| git | verifying the oracle checkout | `git --version` |

`nloptr` is the only non-base R dependency, and it is loaded from a private
library if the system one does not have it:

```sh
mkdir -p "$HOME/rlib"
Rscript -e 'install.packages("nloptr", lib = file.path(Sys.getenv("HOME"), "rlib"), repos = "https://cloud.r-project.org")'
export R_LIBS_USER="$HOME/rlib"
```

The oracle declares `R >= 4.5.0`. The available interpreter here is 4.3.3, and the
harness works around that by sourcing the pinned R files instead of installing the
package — see `docs/compatibility.md`.

## The oracle

`reference/ANCOMBC/` is vendored with its git history, so no network access is
needed. To check it:

```sh
git -C reference/ANCOMBC rev-parse HEAD
# dc4febdf59badb3a8dfe0c767ef2186323c2199a
```

`reference/env/ORACLE.md` records the R version, the BLAS/LAPACK
implementation, the platform, and the RNG configuration used to produce the
committed goldens. If you regenerate on a different machine, record the new
environment alongside the new goldens — a golden is only interpretable together
with the environment that produced it.

## From a clean checkout

```sh
git clone <this repository>
cd ancombc2-rust

# 1. build and run everything that does not need R
cargo test --workspace

# 2. the large golden fixture (release; it is #[ignore]d so a push stays fast)
cargo test --release -p ancombc2-core --test parity -- --ignored --nocapture

# 3. the lints CI enforces
make lint
```

## Regenerating the goldens

Only needed if the oracle changes or the fixture specifications change.

```sh
export R_LIBS_USER="$HOME/rlib"
make smoke     # mirror self-check: the instrumented harness against itself
make goldens   # fx01-fx04 -- REWRITES the fixtures; see below
```

`fx04` is 10,000 taxa x 500 samples with the non-conservative sensitivity analysis
and takes about 30 minutes.

### `make goldens` rewrites the fixtures, and for `fx01`-`fx04` that is wrong

The four committed fixtures predate the `set.seed(spec$seed)` fix in
`reference/R/fixtures.R` and cannot be reproduced. Running the generator on them
writes *different counts*, so:

* the committed goldens no longer match their inputs, and **the parity suite
  fails** — three of four fixtures, first at `y1`, Level A;
* `git diff` reports ~82 changed files, none of which is drift.

Recovery, if you have done it: `git checkout -- validation/golden
validation/fixtures`.

`docs/reference_behavior.md` section 15 has the full account.

### To confirm the goldens are still correct, without rewriting them

```sh
make goldens-drift
```

That target:

* recomputes the contract for `fx01`-`fx04` **from the committed fixtures**
  (`generate_goldens.R --from-committed`) and compares it at each fixture's
  measured text round-trip floor — 0, 4.1e-12, 1.8e-11 and 9.6e-04 relative;
* regenerates `validation/edge` and `validation/matrix` from their fixed seeds and
  compares those **byte-for-byte**, since they are reproducible.

The floors exist because `write.table` renders a double at 15 significant digits,
so the committed fixture *text* cannot represent the in-memory doubles the goldens
were computed from. Bit-exact regeneration of these four is not available. That
this is the text and not the oracle was checked rather than assumed: two runs of
`generate_goldens.R --from-committed 4` over the same text agree to **0** relative
deviation.

`fx04` needs ~30 minutes for the recomputation, so the target is slow by design;
it is a nightly/weekly check, not a per-commit one.

This is what the CI `oracle` job runs.

## Two committed files are over GitHub's 100 MB limit

```
148M  validation/golden/fx04/golden.rds
133M  validation/matrix/golden/shape-10000x500/golden.rds
```

Both are the R serialisations of the two largest golden sets, and both are
committed. **GitHub rejects a push containing a blob over 100 MB**, so adding a
remote to this repository as it stands will fail. This is recorded here rather
than left for whoever clones it to discover.

The files are kept because the golden contract is specified to be stored "in both
.rds and a canonical little-endian f64 format plus JSON metadata", and because
`make goldens-drift` regenerates and writes them. What the *tests* read is the
canonical side: `crates/ancombc2-core/tests/golden/mod.rs` loads the `.f64` and
`.json`, so dropping the `.rds` files would not weaken a single assertion --
`make goldens-drift` would report it, since it regenerates rather than skips.

Options if a remote is wanted, in the order I would pick them:

1. **Git LFS** for `*.rds`. Keeps the contract intact and the working tree
   unchanged; needs a remote with LFS, and `git lfs install` before the first
   clone.
2. **Drop the `.rds` from the tree** and have `make goldens-drift` write them to a
   gitignored directory. Cheap, keeps every assertion, at the cost of the `.rds`
   half of the stated format being generated rather than committed.
3. Ship them as release artefacts. Works, but a clone no longer has the full
   contract.

Option 1 is the one that changes the least and loses nothing.

The repository has no remote and is not published; this is a note for whoever
does that next, not a current failure.

## The golden format

Each fixture directory holds, per quantity:

* a canonical little-endian `f64` blob, so the Rust side needs no R;
* a JSON sidecar for anything that is not a number — `diff_abn`, `zero_ind`,
  `res_global`, `res_pair`;
* `manifest.json`, giving each quantity's file, kind, and shape;
* `golden.rds` and `manifest.rds`, for inspection from R.

Matrices are stored **column-major**, which is R's native order. The core is
row-major; `to_column_major` in the parity harness bridges the two and asserts the
element count.

**These files are large and that is deliberate.** `fx04` is 9,800 taxa x 500
samples x 8 coefficients across ~34 quantities, which is 404 MB of little-endian
`f64` — the contract in the format the plan specifies, uncompressed so the reader
is a plain `read_exact` and there is no decompression step that could differ
between platforms. It compresses to roughly 160 MB, and a gzip variant would
save that, but it would put a codec between the generator and the comparator,
which is the wrong place to introduce one. `benchmarks/datasets/` is a different
matter and is gitignored: 377 MB of pure function of a seed and a script.

## Benchmarks

```sh
make bench-data    # generate the datasets; needs only Python
make bench         # run the surface and evaluate the gates; needs the oracle
make gates         # re-evaluate the gates against the committed results
```

### The containerised surface

`PLAN.md` section 7 asks for a committed container-based harness. It exists —
`benchmarks/container/{Dockerfile,entrypoint.sh,run.sh}`, driven by
`make bench-container` — and pins R 4.5.1, the oracle by commit (verified *in the
build*), the Rust toolchain and `nloptr`, and reports rather than chooses BLAS. It
bind-mounts `benchmarks/results/`, so the image appends to the repository's own
`results.jsonl`, and it evaluates the gates inside the container so each ratio
comes from rows that run wrote.

**It has not been executed here.** The machine that produced the committed
results has neither `docker` nor `podman` on `PATH`, so `run.sh` exits 2 with an
explanation and the host alternative rather than proceeding. Every committed
number therefore came from the host path below, on R 4.3.3 — the same harness,
arms and results file, but not the image's pinned R 4.5.1. See
`benchmarks/container/README.md`.

To run a subset while iterating:

```sh
cargo build --release
python3 scripts/make_bench_data.py --out benchmarks/datasets --only bm1
ANCOMBC_ORACLE_DIR=reference/ANCOMBC \
  python3 scripts/run_benchmarks.py --data benchmarks/datasets \
    --results /tmp/scratch.jsonl --datasets bm1 --arms rust-1,r-1core
python3 scripts/bench_gates.py --results /tmp/scratch.jsonl
```

The results file is append-only and every row carries the host description, the
toolchain versions, the compatibility target and the dataset's index entry, so a
single line is enough to interpret it. `bench_gates.py` keeps the most recent row
per (dataset, arm), so a re-run of one arm is not diluted by earlier rows.

The full surface needs about 60 GB of disk for the datasets (`bm5` alone is
100 million cells) and takes hours with the R arm.

## Profiling

```sh
R_LIBS_USER="$HOME/rlib" Rscript scripts/profile_r.R \
    validation/fixtures/fx03/counts.tsv validation/fixtures/fx03/meta.tsv
```

The decision rule and the current outcome are in `docs/compatibility.md`. The
outcome is **A**, executed on `fx03`, and the committed profile is
`benchmarks/results/Rprof.fx03.summary.txt`. The rule requires the dominant
stage's share to fall below 35% before the optimisation is declared done; it is
still 75.8%, so that gate is red.

## The simulation grid

Results and their interpretation are in `docs/simulation_results.md`.

```sh
make sim          # the quick grid: 48 cells x 20 replicates, both arms
make sim-full     # the PLAN.md surface: 252 cells, 249,600 replicates
```

Three steps, because the two arms must analyse the same tables:

1. `ancombc2-sim generate` writes each replicate's `counts.tsv`, `meta.tsv` and
   `truth.tsv` to `sim-data/`;
2. `ancombc2-sim rust` analyses them in-process;
3. `scripts/sim_r.R` analyses the *same files* through the pinned oracle.

Then `ancombc2-sim summarise` pools the two arms and applies the acceptance rule.
Without the R rows it still reports the Rust arm's calibration and says plainly
that no R arm was supplied — a Rust-only run is a calibration report, not a
comparison, and must not be read as a pass.

To inspect a grid without running it:

```sh
./target/release/ancombc2-sim cells --grid validation/simulation/full
```

That prints the expanded cell list, the per-block rep counts and the projected
replicate total, which is how the full grid's cost can be checked before
committing a machine to it.

Per-replicate rows are committed under `validation/simulation/results/`. They are
large, because the per-taxon vectors are what the pooled metrics are computed
from.

## Real-data validation

```sh
make realdata        # both arms on every prepared dataset, then compare
make realdata-prep   # re-derive validation/realdata from their sources
```

`make realdata` needs no preparation step: `validation/realdata/` is committed, and
both arms read those files. `realdata-prep` is only for re-deriving them, and it
is the step with external requirements.

Two of the four datasets are rebuilt from the pinned oracle checkout, which
`reference/ANCOMBC/` already provides:

* `qmp-real` — the Quantitative Microbiome Project counts, `data/QMP.rda`
* `qmp-vignette` — the vignette's own construction, `sim_plnm` on the QMP
  abundance table, which needs `load_sim_helpers()` from `reference/R/oracle.R`

The other two ship inside a serialized `phyloseq` S4 object in Bioconductor's
`microbiome` package:

```sh
mkdir -p microbiome
curl -sL -o /tmp/microbiome.tar.gz \
  https://bioconductor.org/packages/3.18/bioc/src/contrib/microbiome_1.24.0.tar.gz
tar xzf /tmp/microbiome.tar.gz -C microbiome --strip-components=1
make realdata-prep REALDATA_MB=microbiome/data
```

### The `phyloseq` class-definition shim

Reading those two `.rda` files needs the S4 *class definitions*, because every
base generic on the object dispatches on its class and looks for the defining
package. `phyloseq` itself cannot be installed here: it depends on `RCurl`,
which needs libcurl development headers, and this environment has no root.

`scripts/install_phyloseq_shim.R` writes a minimal package containing only the
class definitions — no code from `phyloseq` — so that the objects can be read.
The reader in `scripts/prepare_realdata.R` then goes through `attributes()`
rather than any `phyloseq` method, so no behaviour of the real package is
relied upon: the slot payloads are ordinary matrices and data frames, and the
shape and names are already in their attributes. If the real `phyloseq` is
installed it takes precedence and the shim is never loaded.

```sh
Rscript --vanilla scripts/install_phyloseq_shim.R
```

The shim declares the union of the slot names phyloseq has used, because the
serialized objects in `microbiome` 1.24 carry `sam_data` and `phy_tree` where the
current package calls them `sample_data` and `tree`.

## Determinism

```sh
./target/release/ancombc2-rs --counts validation/fixtures/fx03/counts.tsv \
    --meta validation/fixtures/fx03/meta.tsv \
    --formula "$(cat validation/fixtures/fx03/formula.txt)" \
    --group group --global --pairwise --struc-zero --pseudo-sens --conservative \
    --out /tmp/run1
./target/release/ancombc2-rs --counts validation/fixtures/fx03/counts.tsv \
    --meta validation/fixtures/fx03/meta.tsv \
    --formula "$(cat validation/fixtures/fx03/formula.txt)" \
    --group group --global --pairwise --struc-zero --pseudo-sens --conservative \
    --threads 1 --out /tmp/run2
diff -r /tmp/run1 /tmp/run2 --exclude=run_metadata.tsv
```

The diff must be empty. `run_metadata.tsv` is excluded because it records the
thread count and the timing, which are supposed to differ.

## What is not reproducible here

* **The containerised harness has not been executed.** It is committed and
  runnable (`make bench-container`), but no container engine is available here, so
  `run.sh` declines and says why. All committed results came from the host path.
* **The `full` simulation grid has no R arm yet.** The Rust arm is complete —
  252 cells, 249,600 replicates, 0 failures, 3,998 s — and its results are written
  up in `docs/simulation_results.md`. The R arm is running and needs on the order
  of 39 hours. Until it exists there is no Rust-versus-R comparison for that grid,
  and nothing in `docs/simulation_results.md` should be read as one; the summary
  file records `"compared_against_r": false`.
* Full oracle regeneration on R >= 4.5.0 has not run here; the declared minimum
  is one series above the installed 4.3.3. `reference/env/ORACLE.md` lists the
  reasons the fixed-effects path is expected to be unaffected, and the weekly CI
  job covers it when such an interpreter exists.
* Two of the four real datasets cannot be re-derived on a machine that cannot
  install `phyloseq`; see the shim above. The datasets themselves are committed,
  so the *comparison* reproduces anywhere — only the preparation needs the shim.
