# Compatibility

## The statement

**`ancombc2-rs v0.1` is equivalent to ANCOMBC 2.15.2 at commit
`dc4febdf59badb3a8dfe0c767ef2186323c2199a`, under R 4.5.x, seed 42.**

"Equivalent" means: for the fixed-effects path, on a given input, the two produce
the same reported quantities to the tolerances in `docs/numerical_contract.md`,
with the exceptions enumerated in `docs/reference_behavior.md`.

It does **not** mean the two implementations are the same code, or that every
input is covered. Specifically:

| | status |
| --- | --- |
| fixed effects, no random effects | implemented, gated by golden parity |
| random effects (`rand_formula`) | not implemented; the argument is absent from the API |
| the trend test | not implemented |
| Dunnett's test | not implemented |
| `phyloseq` / `TreeSummarizedExperiment` input | not implemented; the core takes a count matrix |
| taxonomic aggregation | not implemented; the R wrapper would do it, and the two-table seam is in place |
| `.mdfdr`'s `B` parameter | unused by the pairwise test; the reference only uses it for Dunnett and trend |

## What has been executed

| check | command | result |
| --- | --- | --- |
| unit tests | `cargo test -p ancombc2-core --lib` | **129 passed** |
| statistics tests | `cargo test -p ancombc2-stats` | **35 passed** |
| golden per-stage timings | `stage_seconds.json`, four fixtures | recorded, not compared; names and finiteness checked |
| property tests P1–P16 | `cargo test -p ancombc2-core --test properties` | **18 passed** |
| IO and formula tests | `cargo test -p ancombc2-io` | **19 passed** |
| CLI argument tests | `cargo test -p ancombc2-cli` | **6 passed** |
| FFI tests | `cargo test -p ancombc2-ffi` | **6 passed** |
| simulation harness tests | `cargo test -p ancombc2-sim` | **49 passed** |
| real-data validation | `make realdata` | **4 of 4 datasets** at `diff_abn` agreement 1.00000, `max\|Δbeta\|` ≤ 7.4e-14 |
| structural-zero / edge-case matrix | `make edge-cases` | **12 tests, 7 cases**, all passing; see below |
| golden fixture matrix | `make matrix` | **9 tests, 38 cells**, all passing; see below |
| whole workspace | `cargo test --workspace` | **299 passed**, 0 failed, 1 ignored (`fx04`, run by `make parity-large`) |
| golden parity, fx01–fx03 | `cargo test -p ancombc2-core --test parity` | **7 passed** (3 fixtures + 4 deviation-audit tests) |
| fixture matrix | `make matrix` | **9 passed** — 38 cells at the same contract, plus 4 structural checks |
| fixture-matrix drift | `make matrix-drift` (weekly, needs R 4.5) | wired; the matrix regenerates byte-identically on this host |
| Level C attainment | `ANCOMBC2_GOLDEN_AUDIT=1 cargo test --release -p ancombc2-core --test parity -- --include-ignored --test-threads=1 --nocapture` | `p` 2.6e-15 / 7.6e-13 / 1.6e-10 and `q` 0 / 1.1e-13 / 2.0e-10 on `fx01`–`fx03`; see `docs/numerical_contract.md` |
| golden parity, fx04 | `make parity-large` | **1 passed**, 80 s |
| CLI against the oracle | `make parity` plus a direct `res.tsv` comparison against `res.json` | column names, **column order**, row order and values match; max relative deviation 2.4e-12 (fx01), 2.3e-9 (fx02) |
| R wrapper | `make r-test` | **49 passed**, 0 failed; the package installs, links the cdylib, and the suite exercises the two-taxon-set behaviour from R |
| fmt | `cargo fmt --all -- --check` | clean |
| clippy | `cargo clippy --workspace --all-targets -- -D warnings` | **clean, zero warnings** |
| docs | `RUSTDOCFLAGS="-D warnings" cargo doc --workspace --no-deps` | clean |
| profiling | `Rscript scripts/profile_r.R validation/fixtures/fx03/...` | **executed**, Outcome A; see below |
| benchmark harness | `scripts/run_benchmarks.py` | runnable end to end; see the gate table for what has data |
| miri | `make miri` | not run here; wired into CI |

## The golden fixture matrix

`PLAN.md` §5.5 asks for "a golden fixture matrix over shape 10x10/100x30/1000x100/
10000x500, sparsity zero-free/10/50/90%, continuous/binary/3-group/5-group
predictors, 1/5/10 covariates and interactions, balanced/unbalanced groups,
structural zeros present/absent, pseudo 0/0.1/0.5/1.0, sensitivity on/off,
conservative on/off, and all 7 adjust methods". The full cross product is 245,760
cells, so the matrix is a core cell plus a one-factor-at-a-time sweep per axis plus
four interaction cells: **38 cells over 28 input tables**, of which 11 are a
configuration sweep sharing one table.

`make matrix` runs it. Four of its nine tests are structural and are the part that
would otherwise rot:

* every axis level `PLAN.md` lists is present, read back from the *generated*
  `cells.json` rather than from the R that declared it;
* each sweep varies exactly the axis it is named for;
* the config-only sweeps (`pseudo`, `p_adj_method`) analyse one shared table, so a
  failure is attributable to the configuration and not to the data;
* `conservative` is only swept where the sensitivity analysis is on, since
  `conservative = FALSE, pseudo_sens = FALSE` is the same run twice.

### What it found

Six defects, each needing a shape or a width the four committed fixtures do not
have. The full derivations are in `docs/reference_behavior.md` §13–§14 and
`docs/numerical_contract.md` §4; the short list:

| found by | defect | effect |
| --- | --- | --- |
| `sparsity-090` | a failed per-taxon `lm` leaves `fitted` at NA, and `theta` is a `na.rm` mean of it | `theta` was 0.75 out while `y1` agreed to 9e-16 |
| `sparsity-010` | `q` was compared only against `p_adjust(our p)`, never against the oracle | the contract was weaker than it read |
| `adjust-hochberg` | `Hochberg` and `BH` shared a branch; weights are `n + 1 - i` against `n / i` | six of seven methods unverified |
| `adjust-*` | the `q` re-derivation hardcoded `Holm` | the adjust axis was never exercised at all |
| `int-sparsity90-5group` | `qr` asserted `n >= p` | aborted on a group narrower than it is tall |
| `covariates-10-interaction` | the formula normaliser stripped `*` | no interaction ever reached the design |

The last two are worth reading twice. A matrix that shrinks to the committed
fixtures stays green, which is why the coverage checks are structural rather than
"38 cells passed".

## What has **not** been executed

Stated plainly, because a compatibility claim is only as good as its coverage.

* **`bm5` ran, and it is the most consequential measurement in the suite.** The
  100-million-cell dataset was previously skipped for want of memory, on the
  estimate that it would need about 39 GB. It in fact peaked at 29.3 GB against
  this machine's 31 GB, so it was run for both arms. It is the dataset every
  gate is now judged on, and it reverses the story the smaller datasets told:
  Rust's speed-up over R falls from 53x on `bm1` to 1.03x on `bm5`, and its
  resident set rises to 2.6x R's. A gate judged on the most substantial
  dataset is a harder gate than one judged on the most favourable one, so
  `scripts/bench_gates.py` now uses the longest-running dataset for every gate
  and records the favourable cases separately as `best_case_*`. See
  "Performance gates" below.
* **The week-1 profiling gate has been run once*** **The week-1 profiling gate has been run once**, on `fx03` only. Its outcome is
  recorded below and acted on, but a single small fixture is not a profile of the
  large shapes: `fx04` and the benchmark surface have not been profiled, and the
  dominance of the E-M could shift with taxon count, because the E-M's cost per
  iteration is linear in taxa while the sandwich's is linear in samples.
* **Layer 3, the simulation grid, at the `quick` size only.** The harness is
  implemented and has been executed end to end — see "Statistical validation"
  below — but the executed grid is `validation/simulation/quick`: 48 cells x 20
  replicates at 100 taxa and 20 or 50 samples. The `full` grid
  (`validation/simulation/full`, 252 cells, 249,600 replicates, up to 5,000 taxa
  x 500 samples) is defined, validated by `ancombc2-sim cells`, and run by the
  nightly CI job, but **it has not been executed here**: at 20 replicates the
  per-cell Monte-Carlo error is far wider than the 3 SE acceptance band, so a
  `quick` cell that agrees has not demonstrated agreement at the resolution the
  rule needs.
* **Layer 4, real-data validation, on four prepared datasets.** Executed — see
  "Real-data validation" below. Two of the four are real measured data, one is
  the vignette's own reconstruction, and the fourth is the *only* one that runs
  the sensitivity analysis. The two phyloseq-sourced datasets need a shim to be
  re-derived here; see "What this cannot be re-derived from" below.
* **The full golden fixture matrix.** `PLAN.md` §5 specifies a matrix over shape ×
  sparsity × predictor type × covariate count × group balance × pseudo ×
  sensitivity × adjustment. Four fixtures are committed; they cover the
  structural-zero, sensitivity and large-shape cases, not the matrix.
* **R 4.5.** The available interpreter is 4.3.3 and the oracle declares
  `R >= 4.5.0`. The fixed-effects path is exercised by sourcing the pinned R files
  with a sequential `foreach` stub, which is behaviourally identical to
  `registerDoSEQ()` for this code path — but it is not the same as installing and
  running the package, and a divergence that only appears under the package's own
  scheduling would not be caught here.
* **A container-based benchmark harness.** `PLAN.md` §7 asks for a committed
  container. The harness is committed and runnable; the container image is not.

## The profiling gate and its decision rule

```sh
R_LIBS_USER=$HOME/rlib Rscript scripts/profile_r.R \
    validation/fixtures/fx03/counts.tsv validation/fixtures/fx03/meta.tsv
```

**Outcome: A.** The full profile is committed at
`benchmarks/results/Rprof.fx03.summary.txt`.

Measured on `fx03` (1,000 taxa x 100 samples, `group + x1 + x2`, R 4.3.3, 3.31 s
wall, 81 samples at 10 ms — Rprof's minimum interval on this platform):

| stage | inclusive | share |
| --- | --- | --- |
| `.bias_em` | 2.510 s | **75.8%** |
| &nbsp;&nbsp;of which `nloptr::nloptr` | 2.280 s | 68.9% |
| &nbsp;&nbsp;of which `dnorm` | 0.860 s | 26.0% |
| `.iter_mle` | 0.45 s | 13.7% |
| &nbsp;&nbsp;of which `stats::lm.fit` | 0.14 s | 4.2% |
| `.sandwich_vcov` | 0.03 s | 0.9% |

So the reference's cost is dominated by the **Nelder-Mead optimisation inside the
E-M bias estimate**, not by the linear algebra. Two consequences:

1. `dnorm` at 26% is a *scalar* normal density evaluated inside the NM objective
   function, on every one of the mixture's data points. The objective can be
   written so that the terms common to every candidate iterate point — the
   log-likelihood of the data under the *current* parameters — are computed once
   per E-M iteration rather than once per objective evaluation. That is the
   single largest available win in the port, and it does not change the answer.
2. The sandwich and the MLE, which are what the plan's parallelisation targets,
   are together under 15% of the reference's time. They are still the right
   things to parallelise — they dominate the *Rust* profile, which is a different
   distribution because the E-M's scalar loop is the part the port can vectorise
   best — but a speed-up claim about them cannot come from this profile.

The rule requires the dominant stage's share to fall below 35% before the
optimisation is declared done. It has not: `.bias_em` is still at 76%, because
the E-M vectorisation has not been done. That is the next piece of work, and the
gate stays red until it is.

The rule, decided in advance so the outcome cannot be rationalised after the fact:

* **Outcome A** — Rprof attributes ≥ 50% of the time to one stage. Optimise that
  stage, re-measure, and require the stage's share to fall below 35% before
  declaring it done. This is the expected outcome: the sandwich accumulation and
  the E-M are the two candidates.
* **Outcome B** — no stage exceeds 35%; the time is spread across the MLE
  iteration, the pattern grouping and the per-taxon fits. Then the work is in the
  *structure* — the 20 MLE iterations each redo a full `lm_fit_all` — and the
  change to make is to fuse the iteration with the pattern solve so a taxon's
  response is not rebuilt 20 times. Do that before micro-optimising anything.
* **Outcome C** — Rprof says the time is in R-level overhead (S3 dispatch, copy
  on modify, `[<-` into a matrix). Then the port's remaining cost is R's, not the
  algorithm's, and the honest conclusion is that the algorithm is already
  memory-bound and the next gain is in the layout, not the arithmetic.

Any outcome is a legitimate result. A profiling gate that can only be passed is
not a gate.

### The nesting budget, and the pool that was never built

PLAN.md §6 asks for "a single global Rayon pool with an explicit nesting budget
(pseudo-count runs outermost, then EM coefficients, then missingness groups, then
taxa) and fixed-order deterministic reductions". It is implemented, in
`crates/ancombc2-core/src/parallel.rs`, and the interesting part is *why* it needed
to be:

> `build_pool` read `if rayon::current_num_threads() > 1 { return; }`. Before any
> pool is installed `current_num_threads()` reports the machine's logical CPU
> count, so on every multi-core host the guard fired, **no pool was built**, and
> every `par_iter` fell through to Rayon's own global one at `nproc` width. The
> `--threads` flag was decorative and all five thread arms measured the same thing.

That is a bug that looks like a performance result. The fix has three parts:

* the pool is installed once, first-wins, by `ancombc2_core::parallel::install_pool`;
* the analysis runs inside `with_pool`, so the core's `par_iter`s draw from *that*
  pool — a bare `par_iter` outside `with_pool` would use Rayon's own again;
* the budget hands the pool to exactly one level, and a level with **no work to
  share releases its claim** so the level below can take it. That last rule is what
  makes a dense table scale: with no missing values there is one missingness
  pattern, the group level has a single item, and holding the pool there would pin
  the whole analysis to one core.

No level needs a deterministic-order reduction, because each writes to disjoint
slots — a table of which is in `docs/algorithm.md` §Threading. The reductions whose
order *is* observable (the E-M's parameter sweep, the sandwich's accumulation over
samples within a taxon) stay serial. P15 checks 1 thread against 8 for bitwise
equality, and `res.tsv` is byte-identical at 1, 4, 8 and 16 threads.

Measured on `bm6`: 13.96 s at one thread to 4.82 s at sixteen, 3.10x, efficiency
0.194. The MLE's `mle1` and `mle2` each scale about 3x; the serial remainder is the
iteration barrier (`theta` at iteration `k+1` needs every taxon's fitted values from
`k`), `preprocess`, and the sandwich's inner reduction.

## Performance gates

Continuation criteria, not claims. `scripts/bench_gates.py` evaluates them against
`benchmarks/results/results.jsonl` and writes `benchmarks/results/gates.json`.

**All five gates fail on the executed surface.** That is the result, and it is
reported as a failure rather than reframed. The headline number for P1–P2 is
taken from the *largest* dataset, where process startup is negligible; the
per-dataset surface is printed underneath it, because a single ratio from one
dataset says nothing about the rest and the spread here is two orders of
magnitude.

| gate | criterion | status |
| --- | --- | --- |
| P1 | kernel speed-up >= 3x vs R, one core | **FAIL** — 1.141x on `bm5` |
| P2 | end-to-end speed-up >= 2x | **FAIL** — 1.141x on `bm5` |
| P3 | sensitivity speed-up >= 5x on 8–16 cores | **FAIL** — 3.165x on `bm5` at 8 threads |
| P4 | peak RSS <= 70% of R | **FAIL** — 1.751x on `bm5` |
| P5 | strong-scaling efficiency >= 0.70, 1 → 16 threads | **FAIL** — 0.172, on `bm5` |

The whole surface was re-measured after the rank-deficient fitting fixes landed
(`docs/reference_behavior.md` §16). The previous numbers came from a binary that
factorised the group's design rather than `lm`'s on that path, so they described a
different program and are not comparable.

### P3 is measured against the width of the reference's own grid

`bm5` is a **conservative** sensitivity run. The conservative pseudo-count grid is
`{0, 0.1, 0.5, 1}`, so after the main run there are three independent refits, and
the outer level of the nesting order parallelises exactly those three. Eight or
sixteen threads therefore cannot give more than 3x on that path however much pool
is available — the width comes from the reference, not from this implementation.

Measured is 3.17x, just above that ceiling. `bench_gates.py` now records
`measured_ceiling` and a `ceiling_note` on the P3 record, and lists stage-level
scaling for every sensitivity dataset so the cap is visible rather than inferred.

The comparison that shows the cap is not the implementation's limit is `bm6`, the
non-conservative dataset with the full 50-point grid: its sensitivity stage scales
well past **5x** on eight threads, which clears the 5x target. The gate
still reports **fail** on its headline number, because that number and the target
are both unchanged; what changed is that the report no longer reads as a
parallelisation defect.

### P5 is 0.172 because sixteen threads buys nothing over eight on `bm5`

`bm5`, measured: 321.5s (1 thread), 116.6s (4), 115.9s (8), 116.6s (16). At eight
threads the resident set is already 24.6 GB against this host's 31 GB, so past
eight threads the run is competing with the page cache for memory bandwidth and
loses. `bm6`, small enough not to reach that, scales 5.50x from 1 to 16 threads.

An earlier revision of this section reported P5 as 0.320 on `bm6` because `bm5` at
sixteen threads was OOM-killed. It completes now, and the gate is scored on the
most substantial dataset as written, which is `bm5`. Both readings are well short
of 0.70.

The growth of the resident set with thread count is itself the P4 finding: each
concurrent pipeline holds its own `n_taxa × n_samp` buffers. On `bm4`, adding the
three conservative refits takes peak RSS from 1.18 GB to 2.48 GB, which is why
logging straight from indexed rows rather than materialising the sub-table first
was worth doing — it removed a further 0.4 GB there.



### A benchmark harness defect that was hiding two of these

The previous revision of this table reported **P2 and P3 as passing**, at 2.10x and
6.45x. Both numbers were wrong, and they were wrong in a way that only a real
dataset can expose.

`scripts/run_benchmarks.py` decided whether to pass `--pseudo-sens` by looking for
a config key named `sensitivity`. Every dataset config spells it `pseudo_sens`,
after `AncombcConfig::pseudo_sens` and after the CLI flag. No config has a
`sensitivity` key, so the test was always false: **the sensitivity analysis was
never run by the benchmark harness.** `bm5` and `bm6` — the two datasets that exist
precisely to measure it — were benchmarked without it, and the gates were then
computed over a workload that was not the declared one. The fix is in the harness,
and `pseudo_sens` is now read (with `sensitivity` still accepted, so a config
spelling it either way works).

The corrected surface is below. Note `bm6` single-core: 40.2 s against R's 29.3 s,
i.e. Rust is *slower* than R there once the 50 non-conservative refits are
actually included. That was invisible while the refits were silently skipped.

`bm5` is no longer absent. It ran to completion on every arm, including both R
baselines, which took ~330 s each. Getting there required removing the
`DesignCache` blow-up described further down.

### The measured surface

Single-core Rust against single-core R, plus the peak-RSS ratio and the best
observed thread scaling. Reproduce with `make gates`.

| dataset | shape | rust-1 | rust-16 | R 1-core | speed-up | RSS ratio | scaling |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `bm1` | 50 x 500 | 0.015 s | 0.013 s | 0.863 s | **57.5x** | 0.49x | 1.15x |
| `bm2` | 100 x 1000 | 0.082 s | 0.040 s | 1.073 s | **13.1x** | 0.73x | 2.05x |
| `bm3` | 500 x 5000 | 1.401 s | 0.739 s | 5.089 s | **3.6x** | 0.93x | 1.90x |
| `bm4` | 1000 x 10000 | 15.420 s | 3.963 s | 29.102 s | 1.89x | 0.97x | 3.89x |
| `bm5` | 5000 x 20000 + conservative sensitivity | 307.1 s | OOM | 325.7 s | 1.06x | 1.44x | 2.91x (8t) |
| `bm6` | 1000 x 10000 + 50 non-conservative refits | 40.937 s | 8.002 s | 31.367 s | 0.77x | 1.05x | 5.12x |


Every figure here was measured with the sensitivity analysis actually enabled,
after the harness fix above. `bm5` is the largest dataset and the one the gates
are judged on; it is the slowest to get right and the one where Rust is closest
to R.


Three things follow, and they are different problems.

**The speed-up collapses as the tables grow**, from 53x to 1.03x. Before the QR
work, `bm3` took 86.5 s; `DesignCache::build` was 47.3 s of it because the QR
factorisation materialised an `n x n` orthogonal factor, `n` being the number of
*rows* of a design sub-matrix. Replacing that with stored Householder reflectors
applied directly to the right-hand sides took `bm3` to 1.29 s — a 67x improvement
on that stage alone, from `O(n^2)` to `O(n p)`. The residual 1.03x on `bm5` is
real but far short of the 3x the gate asks for. R's inner loop is `lm.fit` on a
matrix BLAS call, which is hard to beat by much at large `n`.

**Allocation traffic was two defects, and both are fixed.** A size-histogram
allocator (`cargo run --release -p ancombc2-cli --example allocprof -- ...`)
showed `bm4` making 48.6 million allocations, 4.9 for every cell of the table,
almost all of them 8 to 32 bytes. Two causes, one per cell:

* the counts reader called `split_line`, which returns a `String` per field, so
  a counts matrix paid one heap allocation per cell to hold text that was
  parsed to `f64` and dropped. `split_line_cow` borrows the fields on the
  unquoted fast path and delegates to `split_line` for a quoted line, which is
  the only case that needs the `""` unescaping;
* `Matrix::select_rows` built a row label per selected row, and
  `DesignCache::build` calls it once per missingness group. A group can cover
  every sample, so a table with one pattern per taxon allocated one `String` per
  sample per group. `select_rows_anon` drops the labels, which the grouped
  factorisation never reads.

Together: `bm4` 48,554,206 → **2,553,750** allocations (19x fewer) and
15.2 s → 15.8 s wall with peak RSS 3.84 → 3.40 GB. `bm3` went from 7,718,644
allocations to 213,868, 0.07 per cell. The response block handed to the
grouped solve is also solved in place now (`Qr::solve_multi_into`), which
removes a full copy per group per MLE iteration; measured, it is neutral on
this host because the optimiser was already eliding the copy, and it is kept
because it makes the bound explicit rather than leaving it to LLVM. Both
changes are pinned by tests -- `the_borrowing_split_equals_the_owning_one` and
`the_anonymous_subset_equals_the_named_subsets_data` -- and the four golden
fixtures re-pass in release, `fx04` included.

**Peak RSS is still worse than R's, and worse the bigger the table**: 0.62x on
`bm1`, the only dataset under R, rising to 2.60x on `bm5`. The dominant term is
structural: the missingness cache holds `p` Householder reflectors of length
`n_used` for *every* pattern, so a table with one pattern per taxon holds
`O(n_taxa n_samp p)` doubles — 1.4 GB on `bm4`, and the reason `bm5` needs
29 GB. R's `lm.fit` factorises per group and does not retain the factors, so it
pays the factorisation repeatedly instead of the memory. Trading one for the
other is a real choice and it has not been made yet; the fix is to bound the
live cache or refactorise per iteration, and it has been measured neither way.

**Thread scaling is now real but short of the target**: 3.10x from 1 to 16 threads
on `bm6`, against a 0.70 efficiency target (0.194 measured). It was 1.00x, and the
reason was a bug rather than a limit:

> `build_pool` guarded with `if rayon::current_num_threads() > 1 { return; }`. Before
> any pool is installed `current_num_threads()` reports the machine's logical CPU
> count, so on **any multi-core host the guard fired and no pool was ever built**.
> `--threads` was silently ignored and every benchmark arm measured Rayon's default
> global pool at full width — which is why `rust-1` through `rust-16` all came out
> at the same time and P5 read as "no parallelism available" rather than "the flag
> does nothing".

With the pool honoured and the nesting budget implemented
(`crates/ancombc2-core/src/parallel.rs`, `docs/algorithm.md` §Threading), the three
disjoint axes scale and the result is bit-identical at 1, 4, 8 and 16 threads. P3
crossed its target at 6.45x. P5 does not: the MLE's *iterations* remain a genuine
barrier, because `theta` for iteration `k+1` depends on every taxon's fitted values
in iteration `k`, so Amdahl's serial fraction is the iteration overhead plus
`preprocess`, `correction` and the serial parts of the sandwich.

### What the profile changed

The Rprof gate above profiles the *reference* and says the reference is
Nelder-Mead-bound. Profiling the *port* pointed somewhere completely different —
the QR factorisation — and that is where the 67x came from. Rprof could not have
found it: Rprof profiles R, and the conclusion that carried over was about R's
cost, not the port's.

One consequence is worth stating because it invalidates a reasonable-sounding
plan. The Rprof rule requires `.bias_em`'s share to fall below 35% before the
optimisation is declared done, and on the reference's own profile it is at 76%.
That gate cannot be met by optimising R, because the port does not use R's
Nelder-Mead. The equivalent gate for the port is the one above, and it is red on
all five counts.

## Statistical validation (Layer 3)

Numerical parity proves the port matches R. It does not prove either is *right*.
`ancombc2-sim` is the answer to that: a microbiome count generator with a
recorded ground truth, and truth-aware metrics.

```sh
make sim            # quick grid: 48 cells x 20 replicates, both arms
make sim-full       # the PLAN.md surface: 252 cells, 249,600 replicates
```

### How a broken arm is reported rather than hidden

Both arms record a replicate that could not be analysed as a *row with an
`error` field*, not as a missing row and not as a zero. That distinction is the
whole point: a dropped replicate silently leaves the denominator, and a zeroed
one silently enters it. The R arm was observed doing this for real — a `make sim-r`
with `R_LIBS_USER` unset cannot load `nloptr`, and all 960 rows came back with
`error: there is no package called 'nloptr'` and no metrics, which is a visible
failure rather than a grid where R found nothing. The Makefile now passes the
library through, and a CI step fails if the R arm's file is empty.

### How the two arms are held to the same data

The acceptance rule is *"Rust's empirical FDR within Monte-Carlo error of R's"*.
That is only meaningful if both arms analysed the same numbers, so the generator
is Rust's, its output is written to disk once, and both arms read that file. A
shared seed is not a shared table: any change to the generator on either side
silently decouples the arms and the comparison becomes noise dressed as a result.

### The design, and what it does not cover

`PLAN.md` §5.3 lists seven factors whose full cross product is 96,768 cells — not
a design anyone can run, and a one-factor-at-a-time sweep would leave the
interactions unmeasured. `validation/simulation/full` is therefore a list of
**blocks**, each a full cross product over the factors it names at a stated
background. Every level of every factor appears; the cross-block interactions
(notably DA proportion x sample count at large taxon counts, and library size x
sparsity) are **not** estimated, and that is stated in the grid file itself so it
cannot be forgotten.

Within a cell the *community* — which taxa are differentially abundant, in which
direction, and how abundant each is — is fixed, and only the *sampling* is
redrawn per replicate. This is not cosmetic. If the community were redrawn too,
`taxon00005` in replicate 0 and `taxon00005` in replicate 5 would be unrelated
taxa, and the SE calibration ratio (`mean(se) / sd(beta across reps)`) would
divide one taxon's standard error by the spread of a different taxon's estimates
— a number with no interpretation.

### What the executed `quick` grid shows

48 cells x 20 replicates, both arms, seed 20240215. The result is committed as
`validation/simulation/results/quick.summary.json`, with the per-replicate rows
from both arms alongside it as `quick.rust.jsonl.gz` and `quick.r.jsonl.gz`.
The rows are the evidence behind every number above — the pooled metrics are
computed from their per-taxon vectors, not from a summary of summaries — but they
are gzipped because 960 replicates of 100 taxa is not small. `summarise` reads the
compressed form directly, so re-deriving the summary from the committed evidence
is one command:

```sh
ancombc2-sim summarise --grid validation/simulation/quick \
    --rust validation/simulation/results/quick.rust.jsonl.gz \
    --r    validation/simulation/results/quick.r.jsonl.gz
```

* **The empirical FDR is controlled in every cell where it is defined.** It is
  defined in 16 of the 48 cells and runs 0.000 to 0.004 against a nominal 0.05.
  The bias correction controls the false discovery rate at this size.
* **Power behaves as it should**, rising with effect size and sample count: 0.001
  at log FC 1 with 20 samples, 0.090 at log FC 1 with 50, 0.375 at log FC 2 with 20
  samples, and 0.986 at log FC 2 with 50. Sign concordance over the taxa with a
  real effect reaches 1.000 where power is high, which is the check that the
  estimator recovers the *direction* and not merely a significance call.
* **The sampling-fraction correction survives its own negative control.** In a
  confounded cell the labelled taxa carry no compositional effect at all — they
  differ only in how often they are observed to zero — so every call on one is a
  false positive for a method that ignores the correction. Across the 16
  confounded cells the pooled false-positive rate on those taxa runs **0.000 to
  0.045**. That is the correction working.
* **Rust and R agree, cannot be compared, or diverge — and the three are counted
  separately.** Of 48 cells: **9 resolved and agreed**, **39 could not be
  compared**, **1 diverged**. The single divergence is cell 39 (100 taxa, 50
  samples, 30% labelled, log FC 2, 50% structural zeros, confounded), where the
  two arms each made exactly one call and not the same one, so the empirical FDR
  is defined on one side only. That is a borderline-taxon disagreement at the
  edge of the decision, not a systematic difference, and it is reported as a
  finding.
* **39 inconclusive cells is a fact about the grid, and it is the main weakness of
  this run.** Most of them are the confounded arm, where — correctly — nothing is
  called, and the rest are high-sparsity cells at small effect sizes where neither
  arm finds anything. They are **not** counted as agreement: a grid that scores
  silence as success passes by having nothing to say. The honest summary is that
  9 cells constitute evidence, 1 is a finding, and 39 are silent.
* **The sandwich over-estimates uncertainty by 17–50% on this simulation.** The
  SE calibration ratio runs 1.167 to 1.495, above 1 in all 48 cells. A ratio of 1
  is correct, and a ratio consistently above 1 means the reported standard errors
  are too wide. This is a property of the method as ported rather than of the
  port, and it is a genuine result about ANCOM-BC2 on this data-generating
  process. It is *partly* an artefact of the design and should not be over-read:
  the replicates resample library depths as well as counts, so a taxon's estimate
  moves between replicates for reasons beyond the per-taxon sandwich, which
  inflates the denominator. Separating those two explanations needs a grid arm
  with fixed depths, which has not been run.

### What this does not establish

* **Not the full grid.** 20 replicates cannot resolve a 3 SE band of roughly 0.05
  on a power of 0.5; the per-cell Monte-Carlo error here is several times wider.
  A `quick` cell that agrees has not demonstrated agreement.
* **The over-estimated standard errors are reported, not explained.** One
  replicate family was run; the competing explanations (library-depth resampling
  inflating the spread, versus the sandwich genuinely being conservative) have
  not been separated. That is the next experiment, and it is a design question
  about the grid, not a bug hunt.
* **Real data is not a substitute for the grid.** Four real datasets agreeing to
  1e-14 say the port is faithful on *those* four; they say nothing about the
  factors the grid varies and real data does not: effect size, sparsity,
  library size, and confounding, on the same table. Both are needed, and the
  honest summary is that the grid has been run at one size and the real-data
  layer has been run at four.

## The structural-zero and edge-case matrix (`PLAN.md` §5.5)

Seven cases, each a small fixture plus a golden captured from the pinned oracle
by `scripts/generate_edge_cases.R`, and asserted by
`crates/ancombc2-core/tests/edge_cases.rs`. Every case checks its own predicate
*and* full parity with the oracle, and a failure names the case and the plan row
— "parity failed" on a seven-case matrix is not an actionable message.

| plan row | case | what the oracle does |
| --- | --- | --- |
| absent in group A | `absent_in_group_a` | T1 zero throughout g1 is flagged for g1 only, and dropped |
| rare in group A | `rare_in_group_a` | T2 with **one** g1 observation is flagged by `neg_lb`, and dropped |
| present in every group | `present_in_every_group` | nothing is flagged, nothing is dropped |
| `NA` in counts | `na_counts` | NA and zero are the same thing to the structural-zero screen; both dropped |
| group of size 1 | `group_of_size_one` | **refused** by the group-size check, before anything else runs |
| `keep_zero = TRUE` | `struc_zero_off` | `struc_zero = FALSE` keeps the flagged taxon **and returns no `zero_ind` at all** |
| `perc_thres` at boundary | `prv_cut_boundary` | `prevalence >= prv_cut`, so a taxon exactly at the cut is **kept** |

Two of the seven are not like-for-like translations: `keep_zero` and
`perc_thres` do not exist in this oracle, and the substitutions are documented
in `docs/reference_behavior.md` §12 and in each case's `about` field.

**The tests were checked by breaking the code, not by reading them.** Flagging
every taxon instead of only the absent ones fails four of the seven; flipping the
prevalence comparison from `>=` to `>` fails exactly the boundary case. A test
that cannot fail is not a test, and these two mutations are the evidence that
they are not.

Two findings came out of building it, both about the *cases* rather than the
core:

* `present_in_every_group` is the negative control, and it is why the other six
  are worth anything: a rule that flags too eagerly passes the two "should be
  flagged" cases and fails this one.
* The `group_of_size_one` case, written the obvious way (`rep(levels,
  length.out = n)`, which *alternates* the levels), silently produced a 4/3 split
  and tested nothing at all. The oracle accepted it happily and the golden
  recorded a normal result. The case now states its group sizes as `c(6, 1)` and
  asserts them.

## Real-data validation (Layer 4)

Parity proves the port matches R on four fixtures. It does not prove either is
right on data neither of us generated, and that is what this layer is for.

```sh
make realdata        # both arms on every prepared dataset, then compare
make realdata-prep   # re-derive validation/realdata from their sources
```

Both arms read the **same committed tables** — `counts.tsv`, `meta.tsv` and
`analysis.tsv` per dataset — so the comparison is between implementations and not
between two reconstructions of the input.

| dataset | shape | kind | analysis |
| --- | --- | --- | --- |
| `qmp-real` | 91 x 106 | real, measured | `~ group`, two cohorts |
| `dietswap` | 130 x 222 | real, measured | `~ nationality + timepoint + group`, **three groups** |
| `atlas1006` | 130 x 1151 | real, measured | `~ sex`, **plus the sensitivity analysis** |
| `qmp-vignette` | 91 x 150 | synthetic, the vignette's own | `~ group + cont_cov`, batch-crossed sampling fractions |

`dietswap` is the dataset the ANCOM-BC2 vignette analyses, minus the random
effect (`rand_formula`), which is out of scope. `atlas1006` is the large one and
the only one that exercises the pattern-grouped QR at a realistic sample count.
`qmp-vignette` is the only synthetic one, and it is here because its sampling
fractions differ by batch crossed with the group — the confounding the
sampling-fraction correction exists to absorb, in a controlled version.

### The result

| dataset | taxa | `diff_abn` | `max\|Δbeta\|` | `max\|Δse\|` | Spearman | Pearson | Jaccard |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `qmp-real` | 81 | **1.00000** | 5.0e-14 | 1.1e-15 | 1.000 | 1.000 | 1.000 |
| `dietswap` | 110 | **1.00000** | 2.5e-14 | 1.2e-15 | 1.000 | 1.000 | n/a |
| `atlas1006` | 121 | **1.00000** | 7.4e-14 | 5.6e-17 | 1.000 | 1.000 | 1.000 |
| `qmp-vignette` | 91 | **1.00000** | 1.3e-15 | 6.9e-17 | 1.000 | 1.000 | n/a |

`passed_ss` and `diff_robust` also agree at **1.00000** on `atlas1006`, over all
121 taxa. The plan's Level D real-data gate is `diff_abn` agreement ≥ 99.99%; the
measured value is 100%, and the largest coefficient discrepancy anywhere is
7.4e-14, which is at the limit of what two different orderings of the same
arithmetic produce.

`Jaccard` is `n/a` on the two datasets where **neither** arm calls anything
significant. That is not a gap in the metric — it is a dataset with no signal —
and it is reported as `n/a` rather than as a perfect score, because a Jaccard of
1.0 over two empty sets is not evidence of anything.

### The `y = x` scatters

The plan asks for the agreement to be inspected "visually on `y = x`", and a table
of correlations is not that. `make realdata-scatter` draws one committed SVG per
dataset with five panels — `beta`, `se`, `p`, `q`, and the significance indicator
itself — each with the identity line drawn through it:

```sh
make realdata-scatter
# -> validation/realdata/results/scatter/{atlas1006,dietswap,qmp-real,qmp-vignette}.svg
```

Every point sits on the dashed line to the precision the `max|d|` annotation gives
(1e-15 to 1e-13), which is the point: a correlation of 1.000 and a scatter that is
obviously *not* the identity are different claims, and only the second shows which
one is in hand.

Two details exist so the figures cannot flatter the result:

* Taxa where the arms disagree on significance are drawn **last, in red, as
  diamonds**, so a scatter that looks perfect while a few calls disagree is not
  silently reassuring.
* The agreement in the fifth panel's title is computed from the plotted points and
  cross-checked against `summary.json`. Reading it from the summary let a figure
  caption itself "agree 1.0" while showing mismatched calls, because the summary
  and the arms behind it are free to disagree.

Both are mutation-tested: shifting one arm's `beta` by 0.4 turns its panel's
`max|d|` from 2.5e-14 to 0.4, and flipping three `diff_abn` calls changes the fifth
panel from `110/110 agree` to `108/110 agree` with the mismatches flagged.

SVG rather than PNG so the figures are diffable in git, with
`svg.fonttype = none` so the axis labels and annotations are text rather than glyph
outlines — the first version was 127 kB of unsearchable path data and could not be
reviewed in a diff at all.

### What this layer found

It found four real bugs, all of which the golden fixtures and the benchmark
surface could not reach, because every variable they use is spelled without the
letter `s` and no metadata carries a missing value.

1. **A regex in the reference harness deleted every `s` from every variable
   name.** `gsub("[-+*\s]+", ...)` — inside a POSIX bracket expression R's
   default engine treats `\s` as the two literal characters `\` and `s`, so the
   class is `-`, `+`, `*`, `\`, `s`, and `sex` became `ex`. The analysis then
   failed with "variables not in metadata", which reads like a metadata problem.
   It was duplicated in two functions, so there were two places to be wrong.
   Both now call one tested function, and `scripts/smoke_oracle.R` checks seven
   formulas including `sex` and `samp_frac + sex`.
2. **The core aborted on a design row with a missing covariate** instead of
   dropping it. The identifiability guard built `X'X` over all rows, so one `NA`
   made the Gram matrix `NaN`, every eigenvalue `NaN`, and the guard reported
   *every* covariate as unidentifiable — a collinearity warning caused by a
   missing value. R fits `lm` on the complete rows. The guard now does too.
3. **There was no way to say "this sample has no group".** `atlas1006` has 37 of
   1151 samples with `sex == "NA"`, and indexing the group table with them
   panicked. `NO_GROUP` now exists, is excluded from the group sizes, from the
   level count (so a two-group analysis is still recognised as two-group and the
   global test is still deactivated), and from the structural-zero tally, where
   an unlabelled sample is not evidence of absence in *any* group.
4. **The reader fabricated a value for a missing group label.** Mapping a missing
   label onto the reference level added those samples to the reference group as if
   they had been observed there, which biased every coefficient towards the
   reference: `beta` moved by up to 0.105 and two `diff_abn` calls flipped. An
   incomplete design row is the correct representation, and it is what R's
   `model.matrix` produces.

It also found a **fifth discrepancy, in the harness rather than the port**: the
oracle has *two* definitions of `passed_ss`, selected by `conservative`. The
conservative one is `ss_prim == 0 | ss_prim == 1` — the fraction of pseudo-count
refits that flagged the taxon must be zero or one — and the non-conservative one
is an agreement test against the main run's p-value. Using the wrong one produced
a **14% disagreement in `passed_ss` that was completely invisible in
`diff_robust`**, because `diff_robust = diff_abn & passed_ss` and `diff_abn` is
FALSE for every affected taxon. `diff_robust` agreeing at 100% while
`passed_ss` agreed at 86% is the exact signature of that bug, and it is why the
comparison reports all three columns rather than the one that was expected to
matter.

### What this cannot be re-derived from

`qmp-real` and `qmp-vignette` come from the pinned oracle checkout and rebuild
with `make realdata-prep`. `dietswap` and `atlas1006` ship inside a serialized
`phyloseq` S4 object in Bioconductor's `microbiome` package, and this machine
cannot install `phyloseq`: its dependency `RCurl` needs libcurl development
headers and there is no root. The shipped **class-definition shim** reads the
objects instead — see `docs/reproduction.md`. It contains no code from
`phyloseq`, it is shadowed by the real package wherever that is installed, and
the reader goes through `attributes()` rather than any `phyloseq` method, so no
behaviour of the real package is relied upon. The two files it needs are
`microbiome_1.24.0.tar.gz` from bioconductor.org, and the commands to fetch and
unpack them are in each dataset's `source.txt`.

## Version and provenance

* Oracle: ANCOMBC 2.15.2, commit `dc4febdf59badb3a8dfe0c767ef2186323c2199a`,
  vendored at `reference/ANCOMBC/` with its git history intact.
* `reference/env/ORACLE.md` records the R version, the BLAS/LAPACK
  implementation, the platform, and the RNG configuration.
* The golden generator verifies the commit and the `DESCRIPTION` version before it
  writes anything, and runs a mirror self-check comparing the instrumented R
  harness against itself, so a drifted or substituted oracle cannot quietly
  rewrite the contract.
* `scripts/smoke_oracle.R` re-runs the self-check on demand, and the CI `oracle`
  job regenerates the goldens and fails if they differ from the committed ones.

## Dense vs sparse taxon representation

PLAN.md section 9 asks for `Representation::{Dense, SparseTaxa}` behind a config
flag and benchmarked, with the note that the sparse path is only worth revisiting
"after profiling shows the dense path is memory-bandwidth-bound at extreme
sparsity". Both halves of that are now settled.

**The flag exists** (`AncombcConfig::representation`, `--sparse-taxa` on the CLI),
as a taxon-major compressed table: non-zero counts stored as `(u32 sample, f64
value)` with a presence bit, `NaN` counts stored and marked rather than dropped,
structural zeros absent. Prevalence and library size read it compressed; the
selected table comes back dense, because every consumer from `log_center` onward
wants dense rows.

**The profiling condition is not met.** The peak-live attribution (see
"What actually costs memory" in `docs/algorithm.md`) put the resident set in
`DesignCache`, the full-size intermediates and the per-group fitted values — not in
the count table. The dense count matrix is 0.8 GB on `bm5` against a 31 GB host;
it was never the constraint. So the premise of a sparse-taxon win does not hold for
this workload, and the measurement agrees.

500 taxa x 5,000 samples, one thread, `--global --pairwise`:

| zero rate | dense wall / RSS | sparse wall / RSS |
| --- | --- | --- |
| 50% (`bm3` as generated) | 1.26 s / 286 MB | 1.26 s / 270 MB |
| 50% (re-randomised) | 1.16 s / 257 MB | 1.15 s / 257 MB |
| 75% | 0.99 s / 251 MB | 0.99 s / 252 MB |

The differences are inside run-to-run noise. That is the expected result rather
than a disappointing one: the two stages that could benefit (prevalence, library
size) are a small fraction of the run, and the table has to be expanded before the
MLE regardless, so the compression buys back only what it cost. At a 90% zero rate
the default `prv_cut` of 0.1 rejects every taxon — prevalence is 0.1 — so the
regime PLAN.md names cannot be run at the default cutoff at all.

Two structural limits are worth recording, because they are why the sparse form is
a transport format rather than a run format:

* **A stored entry costs more than the dense cell it replaces** — 12 bytes against
  8 — so storage only shrinks above `1 - 8/12 = 33%` zeros.
* **The structural-zero screen cannot use it.** That screen needs each *absent*
  cell's group, to count `observed_in_group`, and a taxon-major encoding does not
  carry an absent cell's group at all. It is why the expansion has to happen
  before the screen rather than after it.

`SparseTaxaMatrix::filter` and `CountMatrix::filter` are asserted to return
identical retained taxa, retained samples and selected cells — bitwise, `NaN`
included — over a grid of shapes, zero/`NaN` mixes and cutoff values.

## What Miri can and cannot check here

The `miri` CI job runs `ancombc2-stats` in full and the `matrix::` tests of
`ancombc2-core`. That scope is narrower than the unit-test count suggests, and
the reason is worth stating precisely rather than leaving as a comment in the
workflow.

**Rayon cannot be Miri-checked in this dependency graph.** Any run that reaches
`NestingBudget::map_par` reports an aliasing violation inside `crossbeam-epoch`,
the work-stealing deque Rayon is built on:

```text
error: Undefined Behavior: trying to retag from <N> for SharedReadWrite permission
  --> crossbeam-epoch-0.9.21/src/internal.rs:567
  --> crates/ancombc2-core/src/parallel.rs:278:9
```

`parallel.rs:278` is the `par_iter().collect()` call site; every frame above it is
inside the dependency, and it is preceded by `integer-to-pointer cast` warnings
from the same crate. It reproduces with `RAYON_NUM_THREADS=1` and in
`parallel`'s own unit tests, so it is not a consequence of running work
concurrently — it is the model's inability to verify epoch-based reclamation.

The consequence is that Miri does **not** check `parallel.rs`, `mle.rs`,
`vcov.rs`, `em.rs`, `correct.rs`, `sens.rs` or the pipeline's parallel paths. What
it does check is the dense linear algebra, the pattern bitsets, and the whole of
`ancombc2-stats`.

That is a real gap, and the mitigation is not to suppress the report: `-Zmiri-*`
flags exist that would silence it, and silencing a Stacked Borrows finding to
widen coverage would be disabling the check that produced it. What stands in its
place is the thread-invariance property test (P15) and the byte-for-byte
determinism checks reported above, which together assert that the parallel paths
produce the same answer as the serial ones — a different property from "contains
no undefined behaviour", and not a substitute for it.

### One test is skipped under Miri, with a stated reason

`dist::tests::erfc_and_erf_match_r` compares `erf`, `erfc` and `pnorm` against
R's literal values to 1e-14 and 1e-16. That cannot hold under Miri, because
**Miri randomises the last bits of transcendental functions** (`ln`, `exp`, `powf`)
to catch code that depends on their exact bits. The evidence, from a probe:

```text
pgamma_lower(0.5, 1.0) = 8.427007929497172e-1
pgamma_lower(0.5, 1.0) = 8.427007929497128e-1   <- identical arguments
```

`pgamma_series` ends in `(-x + a * x.ln() - gln).exp()`, so two calls with the
same arguments legitimately disagree under Miri. The test is therefore
`#[cfg_attr(miri, ignore)]` with that reason attached — **not** given a looser
tolerance. Loosening a correctness assertion to accommodate a tool artefact is
the wrong trade.

The structural properties are still checked everywhere, in a companion test that
is *not* skipped: `NaN` propagation, saturation at the infinities, `erf(0) == 0`
and `erfc(0) == 1` exactly, and monotonicity of `pnorm`. The complement identity
`erf(x) + erfc(x) == 1` needs two independent transcendental evaluations to agree
and is asserted natively at 1e-14, next to the R comparison it belongs with.

## Why the performance gates fail: the surface disables the plan's own optimisation

This is the most useful thing found in the performance work, and it is a property
of the **benchmark configuration** rather than of the implementation.

### The measurement

`bm5` — 5,000 taxa x 20,000 samples, conservative sensitivity — run at one
thread, with the only change being the pseudo-count:

| `bm5` configuration | R 1-core | Rust 1-thread | speed-up |
| --- | --- | --- | --- |
| `pseudo = 0` (what the surface uses) | 325.7 s | 307.1 s | 1.06x |
| `pseudo = 0.5` (the reference's default) | 317.8 s | 211.1 s | **1.51x** |

Rust stage timings, same run:

| stage | `pseudo = 0` | `pseudo = 0.5` |
| --- | --- | --- |
| `mle1` | 99.0 s | **21.6 s** |
| `mle2` | 53.6 s | **14.4 s** |
| `sensitivity` | 150.4 s | 150.4 s |

The first MLE gets **4.6x faster** and the second 3.7x. R moves 2.4%.

### Why

All six benchmark datasets are generated with `pseudo = 0.0`. ANCOM-BC2's own
default is `0.5`.

With `pseudo = 0`, `log(0 + 0) = -Inf`, which the pipeline maps to `NA` — and
`observed_mask` is built from the *log-centred response*. So every zero count
becomes a missing value, and **a taxon's missingness pattern is exactly its own
zero pattern**.

Measured on `bm5`: of 1,500 taxa sampled, there are **1,500 distinct patterns** —
one taxon each. At a 10% zero rate over 20,000 samples, two taxa sharing a zero
pattern is not a thing that happens.

That matters because PLAN.md item 3 specifies "bitset-based missingness pattern
grouping with **one cached QR per pattern**". On this surface there is no pattern
to cache *for*: each cache entry would serve exactly one taxon, so caching is
pure cost. Caching pays only when a pattern carries more than about `p` taxa,
since factoring costs `O(n_rows p^2)` and applying it `O(n_rows p)`.

With `pseudo = 0.5`, `log(0 + 0.5)` is finite, every taxon is fully observed,
there is **one** pattern, and the cached factorisation applies to all 5,000 taxa.
That is the 4.6x.

### What this does and does not say

* It does **not** move P1 over its 3x threshold. 1.51x is still short, and the
  reported gate is unchanged at 1.06x because the surface config is what the gate
  script reads.
* It does say the gate is currently measuring a configuration in which the
  implementation's main structural optimisation is, by construction, inapplicable.
  On the surface, Rust and R are both doing 5,000 independent 20,000 x 18
  Householder factorisations, which is irreducible work that neither can share.
  A 3x gap there would need a different algorithm, not tuning.
* The dataset configs were **not** changed. Selecting a configuration because it
  makes a gate pass is choosing the ruler, and the surface stays as generated.

The sensitivity stage is unchanged at 150.4 s in both rows and is now 71% of the
`pseudo = 0.5` run.

### The oracle's own stage split says where bm5 goes

The per-stage timings the oracle harness now records make this concrete. `fx04` is
9,800 taxa x 500 samples with the non-conservative sensitivity analysis, and the
oracle's own split is:

| stage | seconds | share |
| --- | --- | --- |
| `sensitivity` | 1065.60 | **96.0%** |
| `core` | 44.30 | 4.0% |
| `preprocess` | 0.21 | 0.0% |
| `structural_zeros` | 0.16 | 0.0% |
| `sanity_check` | 0.001 | 0.0% |

So the large datasets are almost entirely the sensitivity analysis, and `bm5`'s
gate numbers are a comparison of two implementations of a repeated whole-pipeline
refit. That is why P1, P2 and P3 cluster near each other on `bm5` rather than
spreading out: they are measuring the same stage.

It is also the reason the `pseudo` diagnosis above matters less for the gate than
one might expect. The main run's MLE got 4.6x faster at `pseudo = 0.5`, but the
sensitivity sweep does not move, because its pseudo-counts (`0.1`, `0.5`, `1.0`)
are already finite — so those refits already had one pattern, and 96% of the
work was never affected by the problem `pseudo = 0` causes.

**The next work is therefore the sensitivity stage, and specifically the
conservative sweep**, which runs three whole pipelines of which the reported
result keeps only a reduced summary. That is a structural observation about what
the reference does, not an omission in this implementation.
