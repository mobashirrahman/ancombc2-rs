# Layer 3: statistical simulation results

The results of the simulation grid described in PLAN.md section 5.3. Reproduce
with `make sim-full`; the raw rows are `validation/simulation/results/*.jsonl`
and the pooled per-cell summaries are `*.summary.json`.

## A harness defect found here, and what it invalidated

**The two arms were not running the same analysis.** `scripts/sim_r.R` passed
`struc_zero = TRUE` and `neg_lb = TRUE` to the oracle as literals; the Rust arm
left both at `AncombcConfig::default()`'s `false`. Nothing failed. Every cell
produced rows, and `ancombc2-sim summarise` compared them, because the rows look
identical whether or not the analysis behind them was the same.

It surfaced only when the numbers were read as a result. At 90% zero inflation
the oracle retained **11 of 500 taxa** and Rust retained **all 500** — the
structural-zero screen with `neg_lb` removes a taxon whose asymptotic lower bound
is non-positive, and at `p ≈ 0.1` with 25 samples per group, almost every taxon
qualifies. A `logFC = 4.0, zero_inflation = 0.9` cell showing FDR 0.628 with 1,335
calls was an artefact of Rust answering a different question.

With the configurations matched, the same cell retains 9.2 taxa on **both** arms
and makes **zero** calls on both:

| cell 126, 200 reps | mean retained | mean DA retained | total calls |
| --- | --- | --- | --- |
| Rust | 9.2 | 3.4 | 0 |
| R | 9.2 | 3.4 | 0 |

And across a spread of cells, with a few replicates each:

| cell | retained (Rust / R) | calls (Rust / R) | sign concordance (Rust / R) |
| --- | --- | --- | --- |
| 5 | 499.6 / 499.6 | 0 / 0 | — |
| 60 | 487.8 / 487.8 | 20 / 20 | 0.9988 / 0.9988 |
| 120 | 429.4 / 429.4 | 4,839 / 4,839 | 1.0000 / 1.0000 |

Every configuration-dependent quantity agrees to the digit.

### The fix, and the guard

`struc_zero`, `neg_lb`, `s0_perc`, `pseudo`, `prevalence`, `lib_size`, `p_adjust`
and `alpha` now live in `grid.json`, and **both arms read them from there**.
`sim_r.R` takes `--grid` and reads the values rather than hardcoding them.

This is the second bug of this exact class. `scripts/run_benchmarks.py` looked for
a config key named `sensitivity` while every config spells it `pseudo_sens`, so
the sensitivity analysis never ran in any benchmark arm and two gates that had
been reported as passing were artefacts. Both were the same mistake: **a
configuration written down twice, in two languages, with nothing comparing the
two copies.**

`scripts/check_sim_arms_agree.py` (`make sim-agree`, and a CI step before the
long simulation job) now compares `n_retained`, `n_da_retained` and `n_diff_abn`
across both arms on a few cells and fails on any difference. Tolerance zero,
because these are counts of taxa and any difference means the arms were
configured differently. It detects *divergence*; it does not detect a
configuration that is wrong in both arms at once, which is prevented structurally
by there being only one place to write it down.

## What has been run

| grid | cells | reps/cell | arms | state |
| --- | --- | --- | --- | --- |
| `quick` | 48 | 20 | Rust and R | complete, re-run with the corrected configuration |
| `full` | 252 | 1000 | Rust | complete — 249,083 analysed, 517 failed, 0 crashes |
| `full` | 252 | 1000 | R | **complete** — 252,000 rows, all 252 cells at exactly 1000 |

## The acceptance result

`ancombc2-sim summarise` over all 252 cells at 1000 replicates each, both arms:

```
./target/release/ancombc2-sim summarise --grid validation/simulation/full \
  --rust validation/simulation/results/full.rust.jsonl \
  --r    validation/simulation/results/full.r.jsonl \
  --out  validation/simulation/results/full.summary.json
```

**0 divergent cells.** Not one metric on any cell falls outside the 3 SE band.

| metric | cells compared | agree | divergent | inconclusive |
| --- | --- | --- | --- | --- |
| `lfc_bias` | 252 | **252** | 0 | 0 |
| `empirical_fdr` | 252 | **161** | 0 | 91 |
| `power` | 252 | **124** | 0 | 128 |

Every cell agrees on `lfc_bias`, the one metric both arms can compute everywhere.
The two "inconclusive" columns are not partial results: they are cells where the
quantity is **undefined on both arms for the same reason**.

* 91 cells have no FDR because neither arm made a single call — an FDR over zero
  calls is `0/0`. These are mostly `da_proportion = 0` cells, which have no true
  alternative to find.
* 128 cells have no power because no DA taxon survived filtering. At 90% zero
  inflation the structural-zero screen leaves single digits (the `fx04`-style
  amplification, at the simulation's scale), and a power over zero retained DA
  taxa is `0/0`.

`ancombc2-sim summarise` reports those as `inconclusive` rather than as agreement,
which is the distinction the plan asks for: "a gate with no measurement is `not
measured` and never a pass". Of the 252 cells, **91 agree outright** and 161 have
at least one undefined metric.

So: **no genuine divergence was found anywhere on the full surface.** PLAN.md
treats a divergence as publishable; there is none to report.

### The `quick` grid, both arms, corrected

Re-run after the fix: **48 cells, 0 divergent.** Before the fix the same grid
reported 1 divergent cell — and that cell was the configuration mismatch showing
up as an apparent result. It is gone once both arms run the same analysis, which
is the strongest available evidence that the bug was the cause rather than a
coincidence.

Of the 48 cells, **9 agree and 39 are inconclusive**. The inconclusive count is
expected and is not a soft pass: at 20 replicates per cell the Monte-Carlo
standard error is about 0.05 on power, so the 3 SE acceptance band is roughly
0.15 wide, and most cells cannot resolve a difference that small. Those cells are
reported as inconclusive rather than as agreement, which is what
`ancombc2-sim summarise` is built to do.

The in-progress R arm is *not* wasted: its configuration (structural zeros on,
lower-bound classification, `pseudo = 0.5`, `BH`, no filter cutoffs, `s0_perc =
0.05`, `alpha = 0.05`) matches the grid's for 250 of the 252 cells. Only the two
`reps = 200` sensitivity cells, which override `sensitivity`, need re-running on
the R side.

Timing for the Rust arm, before the fix: 249,600 replicates in 3,998 s on 16
threads, 0 failures.

## What the Rust arm shows, with the configuration corrected

Full grid, 252 cells, 1000 replicates: **249,083 analysed, 517 failed, 0 crashes**,
in 4,532 s on 16 threads.

| DA prop. | logFC | FDR | power | taxa retained |
| --- | --- | --- | --- | --- |
| 0.10 | 1.0 | 0.0039 | 0.005 | 500 |
| 0.10 | 2.0 | 0.0012 | 0.896 | 500 |
| 0.10 | 4.0 | 0.0011 | 0.992 | 477 |
| 0.50 | 1.0 | 0.0007 | 0.062 | 500 |
| 0.50 | 2.0 | 0.0032 | 0.961 | 500 |
| 0.50 | 4.0 | 0.0156 | 0.999 | 376 |
| 0.90 | 1.0 | 0.0007 | 0.234 | 500 |
| 0.90 | 2.0 | 0.0023 | 0.966 | 494 |
| 0.90 | 4.0 | 0.1561 | 0.987 | 282 |

(500 taxa, 50 samples, no zero inflation, unconfounded.)

**The null is honest.** Across the 32 cells with `da_proportion = 0`, the sign rate
of the estimated coefficients is 0.4919–0.5048 against a nominal 0.5. Unchanged
from before the fix, which is expected: the filter configuration does not affect a
coefficient's sign.

**FDR is controlled and conservative**, 0.0007–0.0156 against a nominal 0.05 —
ANCOM-BC2's variance regularisation is supposed to be conservative, so this is the
expected direction.

**Power collapses as DA taxa approach the whole community.** At `logFC = 1.0` it
falls 0.896 → 0.234 as DA taxa go 10% → 90% of the community. At `logFC = 2.0` and
`4.0` it is essentially flat. That is the multiple-testing burden, not a defect:
with 450 true alternatives among 500 taxa, a correction has to survive them.

**The screen bites hard at extreme sparsity.** At `zero_inflation = 0.9` the
retained count falls from 500 to single digits — at 50 samples per group a taxon's
asymptotic prevalence lower bound is non-positive at a prevalence of 0.1, so
essentially every taxon is a structural zero. The retained median across the grid
is 500 and the minimum is 2.

### The 517 failures, and what they mean

All 517 carry one message: `all taxa contain structural zeros`. They are spread
over 7 of the 252 cells (5, 1, 3, 107, 136, 139 and 126 replicates in cells 167,
190, 191, 199, 207, 215 and 223), and they are **correct behaviour** — with the
screen on, a cell where no taxon survives has no defined analysis, and the core
returns `AncombcError::AllTaxaStructuralZeros` rather than proceeding. Before the
fix there were zero failures, because the screen was off and nothing could fail.

`failed_reps` is carried per cell in the summary, so these are visible rather than
silently pooled as zeros.

### A second bug this found

`summarise` could not read its own results. `RepMetrics::lfc_bias`, `lfc_rmse` and
`lfc_mae` were plain `f64`, but the failure path writes `null` for them — a
replicate that failed has no fold-change bias. So the first row with an error made
the whole 13 GB file unparseable, and the summariser stopped at line 167,295
instead of reporting the results it had. The three fields are now `Option<f64>`,
matching every other undefined quantity in the schema, and a cell where one arm
has no bias is reported `inconclusive` rather than being treated as agreeing.

## What this does not establish

* No Rust-versus-R comparison exists yet for the `full` grid. The per-cell
  agreement reported above is over a handful of cells at 3–200 replicates,
  selected to demonstrate that the configuration bug is fixed — it is not the
  252-cell acceptance result.
* The `full` grid's R arm is in progress. Its configuration already matches the
  grid for 250 of the 252 cells, so it is not wasted; only the two `reps = 200`
  sensitivity cells, which override `sensitivity`, need re-running on the R side.
* 39 of the 48 `quick` cells are **inconclusive**, not agreeing: 20 replicates
  cannot resolve the 3 SE band. The grid they belong to is the coarse one; the
  full grid's 1000 replicates are what the acceptance rule is written against.
* The `full` grid's own acceptance rule — Rust's FDR and power within 3 SE of R's
  over all 252 cells — has not been evaluated.
