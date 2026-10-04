# ancombc2-rs

A Rust reimplementation of ANCOM-BC2, the bias-corrected differential abundance
analysis for microbiome count data.

**Compatibility target: `ancombc2-rs v0.1` is equivalent to ANCOMBC 2.15.2 at
commit `dc4febdf59badb3a8dfe0c767ef2186323c2199a`, under R >= 4.5.0, seed 42.**

That string is a claim with a test behind it, not a slogan: the golden parity
suite compares 29 quantities against the pinned R implementation on four
fixtures, and `docs/reference_behavior.md` records every divergence it found.

## What is implemented

The **fixed-effects** path of ANCOM-BC2, in full:

* sample, taxon and prevalence filtering;
* pseudo-count addition, log transform, per-taxon centring;
* structural-zero detection;
* bitset-based missingness-pattern grouping with one cached QR per pattern;
* the iterative fixed-effect MLE that alternates `beta` and sample-specific
  `theta` (`tol = 0.01`, `max_iter = 20`);
* the blockwise HC0 sandwich covariance, with the per-sample outer products
  `x_i x_i^T` precomputed once;
* the three-component Gaussian-mixture E-M bias estimate, initialised at
  `pi = (0.75, 0.125, 0.125)`, with the variance components optimised by
  Nelder-Mead;
* `delta_wls`, the bias correction, the sampling fractions;
* SAM-style variance regularisation with `s0` at the 5th percentile and
  `W = beta / (se + s0)`;
* Wald p-values and `p.adjust` for all seven methods;
* the general contrast engine, the global quadratic Wald test, and the pairwise
  mixed-directional-FDR test implemented as published;
* the conservative and non-conservative pseudo-count sensitivity analyses.

Not implemented, and absent from the API rather than accepted and ignored:
random effects, the trend test, and Dunnett's test.

## Using it

### The command line

```sh
cargo build --release

./target/release/ancombc2-rs \
    --counts counts.tsv \
    --meta   metadata.tsv \
    --formula 'group + age + bmi' \
    --group group \
    --global --pairwise --struc-zero \
    --pseudo-sens --conservative \
    --out results
```

`counts.tsv` is a counts matrix with taxa as rows, samples as columns, and an ID
column in the corner — what R's `write.table(x, row.names = TRUE)` writes.
`metadata.tsv` is samples as rows, one column per covariate, same layout. The
delimiter is a tab unless the file ends in `.csv`.

`--help` lists every option; each one's default is the reference's default, so a
run with no algorithm flags reproduces `ancombc2(data, meta, fix_formula)` with
the reference's own defaults.

Output: `res.tsv` (the primary table, with the reference's column order and, when
the sensitivity analysis ran, `passed_ss_*` and `diff_robust_*`),
`res_global.tsv`, `res_pair.tsv`, `ss_tab.tsv`, `zero_ind.tsv`, and
`run_metadata.tsv` — the last recording the compatibility target, every option
value, the wall time, peak RSS and bytes allocated, so a result table can be
traced back to what produced it.

### As a library

```rust
use ancombc2_core::config::AncombcConfig;
use ancombc2_core::preprocess::CountMatrix;
use ancombc2_core::matrix::Matrix;

let result = ancombc2_core::ancombc2_run(&counts, &design, group_index, &cfg)?;
println!("{} taxa", result.core.taxa.len());
```

The core has no I/O and no R dependency. `ancombc2-io` reads and writes the
tables; `ancombc2-cli` is the binary above.

## Verifying the compatibility claim

```sh
make test          # unit, property, edge-case, IO, CLI, FFI and small parity
make parity-large  # fx04: 9,800 x 500, 50 sensitivity refits
make edge-cases    # the structural-zero / edge-case matrix, 7 cases
make matrix        # the golden fixture matrix, 38 cells over PLAN.md 5.5's ten axes
make sim           # the simulation grid, both arms, pooled and compared
make r-test        # the R wrapper's testthat suite, 49 tests
make lint          # fmt, clippy with -D warnings, and cargo doc
```

Everything in that table has been run on this checkout: **299 tests across the
workspace**, 0 failed — 139 core unit, 34 statistics, 18 property, 12 edge-case,
19 IO, 6 CLI, 6 FFI, 49 simulation-harness, 7 golden-parity and 9 fixture-matrix,
plus the large `fx04` fixture on request and 49 testthat tests in the R wrapper.
`fmt`, `clippy -D warnings` and `cargo doc -D warnings` are all clean.

The golden fixtures and their reference R implementation are in
`validation/` and `reference/`. To regenerate them — which requires the pinned
ANCOMBC checkout and an R with `nloptr` — run `make goldens`. The generator
verifies the oracle commit and runs a mirror self-check before writing anything,
so a drifted oracle fails there rather than quietly rewriting the contract.

`docs/reference_behavior.md` is the honest account of what does not match. In
short: the reference's quirks are reproduced and listed, and every quantity in the
contract is asserted at its own tolerance on every fixture — including the ones
whose per-taxon sub-designs are exactly rank deficient. That class used to be
carried as a twenty-entry "report rather than assert" list on the stated ground
that the reference's choice of least-squares representative was not reproducible;
it is, `lm.fit` calls `dqrls(pivot = FALSE)` so the column order is the identity
and the dropped column is the last aliased one, and §16 of that document records
the three defects that were actually behind the divergence. The rank-deficient
taxa are still counted and printed on every run — how many there are, how many
have a per-taxon `lm` that fails outright, how many are under-determined — so a
change in that population is visible, but it gates nothing.

### Threading

The plan's nesting order is implemented in `crates/ancombc2-core/src/parallel.rs`:
one global Rayon pool, claimed by exactly one level at a time — pseudo-count runs
outermost, then E-M coefficients, then missingness patterns and groups, then taxa.
A level with no work to share releases its claim, so a table with no missing values
(one missingness pattern, one item at the group level) still scales on the work
inside the group rather than pinning itself to one core.

Every parallel level writes to disjoint slots, so none of them needs a
deterministic-order reduction and none can change a result at a different thread
count. `bm6` goes from 13.96 s at one thread to 4.82 s at sixteen, and `res.tsv` is
byte-identical at 1, 4, 8 and 16.

```sh
./target/release/ancombc2-rs --counts … --meta … --formula … --threads 16 --out /tmp/r16
```

Getting `--threads` honoured at all was the substantive part: `build_pool` guarded
on `rayon::current_num_threads() > 1`, which is true before a pool exists, so on
every multi-core host **no pool was ever built** and the flag fell through to
Rayon's default global pool at full width. Every thread arm of the benchmark suite
had been measuring the same thing.

### What the fixture matrix found

The fixture matrix is the part of the testing plan that earns its keep. Six defects
came out of it, each needing a shape or a design width the four committed fixtures
do not have — among them a per-taxon `lm` fallback whose `NA`s were being written
as `0`, which moved `theta` by 0.75 on a 90%-zero table while the input table
itself agreed to 9e-16; `Hochberg` silently computing `BH`; the adjusted
p-value's re-derivation hardcoded to `Holm`, so six of the seven adjustment
methods were never exercised at all; and a formula normaliser that stripped `*`,
so no interaction ever reached the design matrix. Two of those are cases where the
contract was *weaker than it read* rather than wrong in its conclusions, which is
the worse failure, because it is read as coverage.

Four of the matrix's nine tests are structural: they assert that every axis level
`PLAN.md` lists is still covered, that each sweep varies only the axis it is named
for, and that the config-only sweeps analyse one shared table. Without them a
matrix that shrank back to the committed fixtures would still be green.

## Performance

Performance gates are **continuation criteria, not claims**. The current
outcomes are in `benchmarks/results/gates.json`; regenerate them with
`make bench`, which needs the R arm for the R-relative gates.

Every gate is judged on the **most substantial** dataset in the surface — the one
Rust spends the longest on — rather than on the most favourable one. A gate
evaluated on its best case is not a measurement, and it can pass work that the
real workload fails: P4 read 0.62x, a pass, on `bm1`, where process startup
dominates a 12 MB resident set.

| gate | criterion | measured | on |
| --- | --- | --- | --- |
| P1 | kernel speed-up ≥ 3x vs R on one core | **fail** — 1.14x | `bm5` |
| P2 | end-to-end speed-up ≥ 2x | **fail** — 1.14x | `bm5` |
| P3 | sensitivity-analysis speed-up ≥ 5x on 8–16 cores | **fail** — 3.17x | `bm5`, 8 threads |
| P4 | peak RSS ≤ 60–70% of R | **fail** — 1.75x | `bm5` |
| P5 | strong-scaling efficiency ≥ 0.70 from 1 to 16 threads | **fail** — 0.172 | `bm5` |

These are from a full re-run of all 42 arms against a single binary, re-measured
after the rank-deficient fitting fixes landed (`docs/reference_behavior.md` §16) —
the earlier surface came from a binary that factorised the wrong design on that
path, so quoting it would have been quoting a different program. `bm5` at 16
threads **completes** on this host, so P5 is scored on `bm5` rather than falling
back to `bm6`.

**P3 is capped by its own grid, not by the implementation.** `bm5` is a
*conservative* sensitivity run, and the conservative pseudo-count grid
`{0, 0.1, 0.5, 1}` is three independent refits after the main run. The outer level
of the nesting order parallelises exactly those three, so 8–16 threads can give at
most 3× however much pool is available — the reference's grid fixes the width.
Measured 3.17× is just above that ceiling. `bm6`, the non-conservative dataset with
the full 50-point grid, scales well past 5× on its sensitivity stage — clearing the
target on the path where there is enough independent work. The gate
still reports **fail**; `bench_gates.py` now records `measured_ceiling` and a
`ceiling_note` so the report cannot be misread as a parallelisation defect.

**P5 is 0.172 because 16 threads buys almost nothing over 8 on `bm5`.** At 8
threads the working set is already ~24 GB on a 31 GB host, so past 8 threads the
run competes with the page cache for memory bandwidth. `bm6`, small enough not to
hit that, scales 5.50×. The gate is scored on the most substantial dataset — the
one that hits the memory wall — and that is not counted as a pass.

The resident set grows with the thread count, because each concurrent pipeline
holds its own `n_taxa × n_samp` buffers. Measured peak RSS across the four Rust
widths:

| dataset | rust-1 | rust-4 | rust-8 | rust-16 | R parallel |
| --- | --- | --- | --- | --- | --- |
| `bm4` | 1119 MB | 1149 MB | 1192 MB | 1246 MB | 1250 MB |
| `bm6` | 1164 MB | 1467 MB | 1957 MB | 2921 MB | 1224 MB |
| `bm5` | 14928 MB | 22279 MB | 24747 MB | 24743 MB | 8526 MB |

`bm6` is the clearest: 2.5x from one thread to sixteen. P4 is scored on the most
substantial dataset, which is `bm5` — 14928 MB against R's 8526 MB, 1.75x — and
`bm5` is also the one that saturates memory, so the gate lands on the worst case
twice over. `bm4` and `bm6` are at or under R's footprint at one thread. See
`docs/performance_plan.md` for what would and would not move this.


**All five gates fail. That is the result, and it is reported as failure.** A gate
with no measurement is `not measured` and never a pass; the per-dataset surface
behind each number is in `docs/compatibility.md`, because a single ratio from one
dataset says nothing about the rest.

**`docs/performance_plan.md` is the plan for changing that**, and it starts by
diagnosing which of the five are fixable. Briefly: P1 and P2 are reachable, and
they are one bottleneck — the Householder factorisation, which is 84 % of a large
run and which this implementation and the reference execute at about the same
rate. P4 is a memory problem, not a speed one. P5 is a host limit on this machine.
P3 is capped at 3x by the reference's own three-refit grid, which makes it a
scoring question rather than an optimisation. The plan separates those rather than
promising five green gates.

An earlier revision of this table recorded **P2 and P3 as passing** (2.10x and
6.45x). Those numbers were produced by a defect in `scripts/run_benchmarks.py`:
it decided whether to pass `--pseudo-sens` by testing a config key named
`sensitivity`, while every config spells it `pseudo_sens`. The test was therefore
always false and **the sensitivity analysis never ran in any benchmark arm** — so
`bm5` and `bm6`, the two datasets that exist to measure it, were timed without it,
and the gates were computed over a workload other than the declared one. The
harness now reads `pseudo_sens`; with that fixed, single-core Rust on `bm6` is
40.2 s against R's 29.3 s, i.e. *slower* than R, which the broken harness had been
hiding. The corrected table above is the one to read.

The `quick` simulation grid has been re-run with both arms correctly configured:
**48 cells, 0 divergent** (it previously reported 1 divergent cell, which was the
configuration bug below showing up as an apparent result). 9 cells agree and 39
are inconclusive — 20 replicates cannot resolve the 3 SE band, and inconclusive is
reported as such rather than as a pass.

The **full** simulation grid has run on both arms: **252 cells x 1000 replicates,
0 divergent cells**. Every cell agrees on `lfc_bias` (252/252), and on
`empirical_fdr` and `power` wherever those are defined (161/161 and 124/124) --
the rest have no calls or no surviving DA taxa on *either* arm, so the quantity is
`0/0` and is reported as inconclusive rather than as agreement. Details in
`docs/simulation_results.md`.

The **full** simulation grid is being run and its results are written up in
`docs/simulation_results.md`, which also records a harness defect found there:
the two arms were passing different `struc_zero`/`neg_lb` settings and so were
not analysing the same thing — at 90% zero inflation the oracle retained 11 of
500 taxa where Rust retained 500. Every parameter both arms use now lives in
`grid.json`, `make sim-agree` checks the two agree before the long run, and a CI
step gates on it. The same class of bug had silently disabled the sensitivity
analysis in every benchmark arm.

**Why the gates fail is now measured, and it is partly the surface.** All six
benchmark datasets use `pseudo = 0.0`, so `log(0) = NA` and a taxon's missingness
pattern is its own zero pattern — on `bm5`, 1,500 sampled taxa have 1,500 distinct
patterns, one taxon each. The per-pattern QR cache the plan specifies therefore has
nothing to cache for, and both implementations perform 5,000 independent
20,000x18 Householder factorisations. Re-run at ANCOM-BC2's own default
`pseudo = 0.5`, where every taxon lands in one pattern, Rust's first MLE goes
99.0 s -> 21.6 s and the speed-up over R widens from 1.06x to **1.51x**. The
dataset configs were not changed; selecting one because it makes a gate pass is
choosing the ruler. See `docs/compatibility.md`.

PLAN.md section 9's dense-vs-sparse choice is now measured rather than assumed:
`Representation::{Dense, SparseTaxa}` is behind `--sparse-taxa`, and at 50% and
75% zero rates it changes neither wall time nor resident set outside run-to-run
noise. The profile says why — the count table was never where the memory went.
Details and the two structural limits are in `docs/compatibility.md`.

`bm5` (5,000 x 20,000 with the conservative sensitivity analysis) is **no longer
absent**: it completes on every arm, including both R baselines at ~330 s each.
Reaching that meant removing a `DesignCache` that held a design and its QR
reflectors for each of ~1000 missingness patterns — 5.25 GB on that dataset, and a
bad trade on time as well, since caching only pays when a pattern carries more
than about `p` taxa and real tables have roughly one pattern per taxon. See
"What actually costs memory" in `docs/algorithm.md`.


Three performance defects were found and fixed. Two came from allocation
profiling rather than inspection, each costing one heap allocation per cell of the
input table: 48.6 million allocations on `bm4`, now 2.6 million. The third was
worse than a missing optimisation — `build_pool` guarded on
`rayon::current_num_threads() > 1`, which is *always* true before a pool exists, so
**no pool was ever built on a multi-core host** and `--threads` was silently
ignored. Every benchmark arm had been measuring the same thing at full width, which
is why thread scaling read as 1.00x.

With the pool honoured and the plan's nesting budget implemented, `bm6` goes from
13.96 s to 4.82 s on 16 threads with **byte-identical output at 1, 4, 8 and 16**.
P5 still fails at 0.194 because the MLE's iterations are a genuine barrier, and P4
still fails because the per-missingness-pattern reflector cache dominates the
resident set. Both are named as open in `docs/compatibility.md`.

**Profiling outcome: A, executed on `fx03`.** 75.8% of the reference's runtime is
`.bias_em`, and 68.9% of the total is `nloptr::nloptr` inside it — the
Nelder-Mead optimisation of the mixture's variance components, not the linear
algebra. The sandwich and the MLE together are under 15%. The full profile is in
`benchmarks/results/Rprof.fx03.summary.txt` and the decision rule is in
`docs/compatibility.md`; the rule requires the dominant stage to fall below 35%
before the optimisation is declared done, and it has not, so the gate is still
red.

Profiling the *port* then found something Rprof could not: the QR
factorisation was materialising an `n x n` orthogonal factor, and replacing it
with stored Householder reflectors applied directly to the right-hand sides took
the 500x5000 dataset from 86.5 s to 1.29 s — a 67x improvement on that stage
alone, and `O(n^2)` to `O(n p)`.

## Statistical validation

Parity proves the port matches R. It does not prove either is *right*, so
`ancombc2-sim` generates microbiome counts with a recorded ground truth and
measures whether the method recovers what was planted:

```sh
make sim          # the quick grid: 48 cells x 20 replicates, both arms
make sim-full     # the PLAN.md surface: 252 cells, 249,600 replicates
make realdata     # Layer 4: four real datasets, both arms, compared
```

Both arms analyse the *same* generated tables — the generator is Rust's, its
output is written once, and R reads that file — so the comparison is between
implementations rather than between two draws from the same description.

On the executed quick grid the empirical FDR runs 0.000–0.004 against a nominal
0.05; power rises from 0.001 at log FC 1 with 20 samples to 0.986 at log FC 2
with 50; and the sampling-fraction correction keeps the false-positive rate on
the grid's own confounded negative control between 0.000 and 0.045. Of 48 cells,
9 resolved and agreed, 1 diverged (reported as a finding), and 39 could not be
compared because neither arm called anything — those are counted separately and
**not** as agreement. The sandwich also over-estimates uncertainty by 17–50% on
this data-generating process, which is reported rather than smoothed over.
`docs/compatibility.md` has the numbers, the design, and what this does not
establish.

## Edge cases

`PLAN.md` §5.5's structural-zero matrix is generated from the pinned oracle and
run as seven named cases: absent in one group, a single observation in one group,
present everywhere, `NA` against zero, a group of size 1, the
`struc_zero`-off equivalent of the plan's `keep_zero`, and a taxon sitting
exactly on the prevalence cut. Each asserts its own predicate *and* parity with
the reference, and a failure names the case and the plan row.

```sh
make edge-cases   # assert the committed goldens
make edge-gen     # re-derive them from the oracle
```

The tests are mutation-checked: flagging every taxon fails four of the seven, and
flipping the prevalence comparison from `>=` to `>` fails exactly the boundary
case.

## Real-data validation

Parity proves the port matches R on four fixtures; it does not prove either is
right on data neither of us generated. `make realdata` runs both arms over four
prepared datasets — the Quantitative Microbiome Project counts, the diet-swap
study the ANCOM-BC2 vignette analyses, a 1,151-sample atlas (the only one that
runs the sensitivity analysis), and the vignette's own synthetic
reconstruction — and compares them.

| dataset | shape | taxa | `diff_abn` | `max\|Δbeta\|` | Spearman |
| --- | --- | --- | --- | --- | --- |
| `qmp-real` | 91 x 106 | 81 | **1.00000** | 5.0e-14 | 1.000 |
| `dietswap` (3 groups) | 130 x 222 | 110 | **1.00000** | 2.5e-14 | 1.000 |
| `atlas1006` (+ sensitivity) | 130 x 1151 | 121 | **1.00000** | 7.4e-14 | 1.000 |
| `qmp-vignette` | 91 x 150 | 91 | **1.00000** | 1.3e-15 | 1.000 |

`passed_ss` and `diff_robust` also agree at 1.00000 on `atlas1006`. The plan's
Level D real-data gate is 99.99%; the measured value is 100%.

**This layer earned its keep.** It found four real bugs that four golden fixtures
and six benchmark datasets could not, because every variable they use is spelled
without the letter `s` and no fixture has a missing covariate: a regex in the
reference harness that deleted every `s` from every variable name; a core that
aborted instead of dropping a design row with a missing covariate; no way at all
to represent a sample whose group label is missing; and a reader that fabricated
a reference-level value for one. It also exposed that the oracle has *two*
definitions of `passed_ss` and that using the wrong one is invisible in
`diff_robust`. All five are written up in `docs/compatibility.md`.

## Layout

```
crates/ancombc2-core    numerics: no I/O, no R
crates/ancombc2-stats  distributions, type-7 quantile, p.adjust
crates/ancombc2-io     TSV/CSV matrices, metadata, formula parsing, result tables
crates/ancombc2-cli    the ancombc2-rs binary
crates/ancombc2-ffi    a C ABI for the R wrapper: JSON in, JSON out
crates/ancombc2-sim    the simulation grid, its metrics, and the Rust-vs-R verdict
reference/             the pinned ANCOMBC, the R harness, the environment record
validation/            fixtures, goldens, simulation and real-data results
benchmarks/            the scaling surface and its results
scripts/               fixtures, benchmarks, gates, profiling
docs/                  algorithm, numerics, reference behaviour, statistics
r/ancombc2rs/          the R wrapper package
```

## Documentation

* `docs/algorithm.md` — what the algorithm does, step by step.
* `docs/numerical_contract.md` — layouts, tolerances, and the parity levels.
* `docs/reference_behavior.md` — every divergence from the oracle, and why.
* `docs/statistical_spec.md` — the estimators, the assumptions, the tests.
* `docs/compatibility.md` — the compatibility statement, the measured benchmark
  surface, the gate outcomes, the profiling outcome, the simulation results, and
  the real-data validation.
* `docs/reproduction.md` — how to regenerate everything from scratch.

## Licence

The Rust code is licensed under the same terms as the vendored reference
(Artistic-2.0); see `LICENSE`.
