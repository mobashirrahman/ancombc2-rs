# ancombc2-rs — Development Plan

A numerically compatible, high-performance Rust reimplementation of the **fixed-effects**
core of ANCOM-BC2, with a rigorous validation and benchmarking program sufficient to
support a computational-methods manuscript.

---

## 0. Scope, goals, non-goals

### 0.1 One-sentence goal

> `ancombc2-rs` reproduces ANCOM-BC2's fixed-effect inference to floating-point parity,
> at materially lower runtime and memory, for microbiome count matrices of up to
> ~20,000 taxa × ~5,000 samples, with sensitivity analysis enabled.

### 0.2 Version 1.0 feature boundary

**IN SCOPE (v1.0)**

| Area | Item |
|---|---|
| Input | dense `taxa × samples` count matrix, design matrix, group vector, taxon IDs |
| Filtering | sample filtering, prevalence filtering, library size |
| Transform | pseudo-count addition, log transform, per-taxon centering |
| Zero handling | missingness patterns, structural-zero detection (`zero_ind`) |
| Estimation | iterative MLE (`beta*`, `theta`), sandwich (HC0) variance |
| Bias | three-component Gaussian-mixture EM (`delta_em`), WLS bias (`delta_wls`) |
| Correction | bias-corrected `beta`, corrected `lfc`/`se` |
| Regularization | SAM-style `s0` (5th percentile), `W` statistic |
| Inference | Wald p-values, `p.adjust` (Holm, Hochberg, Hommel, Bonferroni, BH, BY, none) |
| Testing | global test, pairwise tests, mixed directional FDR (`mdFDR`) |
| Sensitivity | conservative and non-conservative pseudo-count sensitivity analysis |
| Integration | R package wrapper (phyloseq / TreeSummarizedExperiment / formula interface) |
| Perf | multithreading, reusable buffers, streaming/blocked memory layout |

**OUT OF SCOPE for v1.0 (deferred)**

- random intercepts / random slopes / general `lmer` compatibility
- trend (ordered-pattern) constrained testing
- Dunnett multiple comparison
- S4 / phyloseq object model in Rust
- R formula parsing, factor/contrast semantics in Rust

### 0.3 The governing principle

> **Parity first, speed second.**

ANCOM-BC2 is inference software. A 30× implementation that silently shifts
false-positive behaviour is a failed port. A 5–10× implementation with exhaustive
numerical *and* statistical equivalence is a credible, publishable result.

Every performance claim in this document is a **target with a gate**, never a promise.

---

## 1. System architecture

```
                     ┌───────────────────────────┐
                     │  R / phyloseq / TSE users │
                     │  formula, contrasts, S4   │
                     └─────────────┬─────────────┘
                                   │  model.matrix(), factor coding
                     ┌─────────────▼─────────────┐
                     │     ancombc2rs (R pkg)    │
                     └─────────────┬─────────────┘
                                   │  counts: taxa×samples (F-ordered or C)
                                   │  design: samples×p
                                   │  group: samples
                                   │  config struct
                     ┌─────────────▼─────────────┐
                     │        ancombc2-core       │
                     │                           │
                     │  preprocess               │
                     │      ↓                    │
                     │  log transform / center   │
                     │      ↓                    │
                     │  missingness pattern idx  │
                     │      ↓                    │
                     │  iterative fixed-effect   │
                     │  MLE  (QR per pattern)    │
                     │      ↓                    │
                     │  sandwich variance        │
                     │      ↓                    │
                     │  EM bias estimation       │
                     │      ↓                    │
                     │  bias correction          │
                     │      ↓                    │
                     │  sampling fractions       │
                     │      ↓                    │
                     │  variance regularization │
                     │      ↓                    │
                     │  inference (W, p, q)      │
                     │      ↓                    │
                     │  global / pairwise tests  │
                     │      ↓                    │
                     │  pseudo-count sensitivity │
                     └─────────────┬─────────────┘
                                   │  result struct (borrowed, zero-copy)
                     ┌─────────────▼─────────────┐
                     │  R tibble / TSV / CLI out │
                     └───────────────────────────┘
```

**Hard rule:** `ancombc2-core` must not depend on R, on `Rcpp`, or on any R-specific
type. It consumes plain `f64` buffers and a plain config struct. This is what makes a
Python binding, a CLI, or a WASM build possible later.

---

## 2. Repository layout

```
ancombc2-rs/
├── Cargo.toml                     # workspace
├── rust-toolchain.toml            # pinned toolchain
├── crates/
│   ├── ancombc2-core/             # pure numerics; no I/O, no R
│   │   ├── src/
│   │   │   ├── lib.rs
│   │   │   ├── config.rs          # AncombcConfig, CompatMode
│   │   │   ├── error.rs
│   │   │   ├── matrix/
│   │   │   │   ├── mod.rs
│   │   │   │   ├── dense.rs       # row-major taxa×samples
│   │   │   │   └── pattern.rs     # bitset missingness patterns
│   │   │   ├── preprocess/
│   │   │   │   ├── mod.rs
│   │   │   │   ├── filter.rs      # sample/taxon/prevalence
│   │   │   │   ├── transform.rs   # pseudo, log, center
│   │   │   │   └── zeros.rs       # structural zero_ind
│   │   │   ├── mle/
│   │   │   │   ├── mod.rs
│   │   │   │   ├── iterative.rs   # beta/theta loop
│   │   │   │   └── qr.rs          # cached per-pattern QR
│   │   │   ├── vcov/
│   │   │   │   ├── mod.rs
│   │   │   │   └── sandwich.rs    # blockwise HC0
│   │   │   ├── em/
│   │   │   │   ├── mod.rs
│   │   │   │   ├── mixture.rs     # 3-component GMM EM
│   │   │   │   ├── wls.rs         # delta_wls
│   │   │   │   └── nelder_mead.rs # variance param optimisation
│   │   │   ├── correct/
│   │   │   │   ├── mod.rs
│   │   │   │   ├── bias.rs        # beta - delta
│   │   │   │   ├── samp_frac.rs
│   │   │   │   └── regularize.rs  # s0, W
│   │   │   ├── test/
│   │   │   │   ├── mod.rs
│   │   │   │   ├── wald.rs
│   │   │   │   ├── global.rs      # quadratic Wald, contrast engine
│   │   │   │   ├── pairwise.rs
│   │   │   │   └── mdfdr.rs       # mixed directional FDR
│   │   │   ├── sens/
│   │   │   │   ├── mod.rs
│   │   │   │   ├── conservative.rs
│   │   │   │   └── nonconservative.rs
│   │   │   ├── stats/             # p.adjust, chsq, norm, quantile
│   │   │   ├── workspace.rs       # reusable buffers
│   │   │   └── rayon_plan.rs      # nesting policy
│   │   └── benches/
│   ├── ancombc2-stats/            # distributions, p-adjust, contrasts
│   ├── ancombc2-io/               # TSV/CSV/matrix IO, config parsing
│   └── ancombc2-cli/              # `ancombc2-rs` binary
├── r/ancombc2rs/                  # R wrapper package
│   ├── R/  src/  tests/  man/  DESCRIPTION  NAMESPACE
├── reference/
│   ├── R/                         # generate_r_reference.R, harness fns
│   ├── env/                       # renv.lock, sessioninfo.txt, BLAS record
│   └── README.md                  # exact oracle version + commit
├── validation/
│   ├── golden/                    # .rds + .json expected outputs
│   ├── fixtures/                  # synthetic dataset generators (R + Rust)
│   ├── simulation/                # FDR/power simulation harness
│   └── realdata/                  # dataset A/B/C adapters
├── benchmarks/
│   ├── datasets/                  # generator scripts, sizes
│   ├── run_bench.sh
│   ├── run_bench.R
│   └── results/
├── scripts/
│   ├── profile_r.R                # Rprof breakdown (Outcome A/B/C, §38)
│   ├── record_env.sh
│   └── promote_fixture.sh
├── docs/
│   ├── algorithm.md
│   ├── numerical_contract.md      # tolerances, Levels A–D
│   ├── reference_behavior.md      # quirks of the oracle
│   ├── statistical_spec.md        # what the paper says
│   ├── compatibility.md
│   └── reproduction.md
├── .github/workflows/ci.yml
├── Makefile
└── README.md
```

---

## 3. Phase 0 — Freeze the reference implementation (Week 1, days 1–2)

**Nothing is written in Rust until the oracle is frozen.**

`reference/env/` must record:

```
ANCOMBC version : 2.15.2
Git commit      : <exact 40-char SHA>
R version       : 4.5.x (full string)
Platform        : x86_64-pc-linux-gnu
BLAS/LAPACK     : e.g. OpenBLAS 0.3.27 with THREADS=32  (record sessionInfo()$BLAS)
RNG             : set.seed(42), RNGkind("Mersenne-Twister","Inversion","Rejection")
```

Compatibility target string used everywhere (CI, README, crate version):

```
ancombc2-rs v0.1 ≡ ANCOMBC 2.15.2 @ <SHA>, R >= 4.5.0, seed 42
```

> **Why this matters.** ANCOM-BC2 behaviour is *not* stable across releases. For
> example, current `conservative = TRUE` sensitivity analysis reruns the full
> algorithm for pseudo-counts {0, 0.1, 0.5, 1.0}, whereas `conservative = FALSE`
> estimates sampling fractions once and performs 50 inference refits with
> pseudo-counts 0.01…0.50. "Same as ANCOM-BC2" is therefore not a testable claim
> until it is pinned to a version + commit + RNG.

### 3.1 Reference harness

```r
# reference/R/harness.R
run_reference <- function(counts, meta, formula, cfg = list()) {
  set.seed(42)
  cfg <- utils::modifyList(list(
    p_adjust = "holm", tol = 0.01, max_iter = 20L, detail = FALSE,
    p_corr = FALSE, alpha = 0.05, s0_perc = 0.05,
    n_init = 1L, perc_thres = 0, keep_zero = FALSE,
    global = FALSE, pairwise = FALSE, pseudo_sens = FALSE,
    conservative = FALSE, protect = NULL, sample_thres = 0
  ), cfg)
  do.call(ANCOMBC::ancombc2, c(
    list(obj = phyloseq::as_phyloseq(counts, meta), formula = formula), cfg))
}
```

### 3.2 The golden contract

Validation must **never** be against final p-values alone. Capture **every
intermediate** the reference exposes, plus the ones it does not expose (obtained by
instrumenting a forked copy of the R source).

| # | Quantity | Source | Golden file |
|---|---|---|---|
| 1 | taxa retained (ids, in order) | `feature_table$taxon` | `taxa.rds` |
| 2 | samples retained (ids, in order) | `feature_table$sample` | `samples.rds` |
| 3 | structural zero flag matrix | `ss_tab$zero_ind` | `zero_ind.rds` |
| 4 | processed counts (post pseudo/log) | instrumented | `y_raw.rds` |
| 5 | per-taxon centered `Y` | instrumented | `y_centered.rds` |
| 6 | missingness pattern assignment | instrumented | `patterns.rds` |
| 7 | `beta*` (pre-correction) | instrumented | `beta_star.rds` |
| 8 | `theta` | `result$theta` | `theta.rds` |
| 9 | residuals `epsilon` | instrumented | `residuals.rds` |
| 10 | iteration count + final `epsilon` | instrumented | `conv.rds` |
| 11 | sandwich vcov block | instrumented | `vcov.rds` |
| 12 | `delta_em` per coefficient | `result$delta_em` | `delta_em.rds` |
| 13 | `delta_wls` per coefficient | `result$delta_wls` | `delta_wls.rds` |
| 14 | EM mixture params (`pi`, `mu1/2`, `l1/2`, `kappa1/2`) | instrumented | `em_params.rds` |
| 15 | bias-corrected `beta` | instrumented | `beta_corr.rds` |
| 16 | sampling fractions `samp_frac` | `result$samp_frac` | `samp_frac.rds` |
| 17 | `se_*` | `feature_table` | `se.rds` |
| 18 | `s0` | instrumented | `s0.rds` |
| 19 | `W_*` | `feature_table` | `w.rds` |
| 20 | `p_*`, `q_*` | `feature_table` | `pq.rds` |
| 21 | `diff_abn` | `feature_table` | `diff_abn.rds` |
| 22 | global test `res_global` | `result$res_global` | `global.rds` |
| 23 | pairwise `res_pair` | `result$res_pair` | `pairwise.rds` |
| 24 | `ss_tab`, `passed_ss` | `result$ss_tab` | `ss_tab.rds` |
| 25 | `diff_robust` | `feature_table` | `diff_robust.rds` |
| 26 | per-stage wall times | `Rprof` | `timing.rds` |

Each golden file stored as `.rds` (R) **and** a canonical `.npy`-equivalent
little-endian `f64` blob + JSON metadata, so the Rust test harness never needs R.

When a Rust quantity disagrees, the golden table tells you **exactly** which stage
diverged. This is the single most valuable piece of infrastructure in the project.

---

## 4. Algorithm specification (implementation checklist)

This section is the contract the Rust code must satisfy. `docs/algorithm.md` is
the long form; this is the checklist.

### 4.1 Preprocessing

1. Sample filter: drop samples with library size < `sample_thres`.
2. Taxon filter: prevalence < `perc_thres` removed (unless `keep_zero`).
3. Structural zeros: taxon with **no positive count in at least one group**
   → flagged in `ss_tab$zero_ind`, reported separately; excluded from the primary
   model path (see ANCOMBC `keep_zero = FALSE` semantics).
4. `Y = log(x + pseudo)`, `pseudo = 0` by default.
5. Center: `Y_c[i,j] = Y[i,j] - mean_j Y[i,j]` — removes taxon-specific
   sequencing efficiency.
6. `theta = 0` initially (sample-specific nuisance/bias term).

### 4.2 Missingness structure

With `pseudo = 0`, observed zeros are `log 0 = -Inf` and are **unavailable
observations** on the log scale. Each taxon therefore has its own usable sample
subset. Current ANCOMBC groups taxa by their usable-sample pattern and does one
multi-response `lm.fit` per pattern. Represent as bitsets:

```rust
pub struct MissingnessPattern { words: Vec<u64> }   // bitset over samples
// taxon A 1111011111011 ┐
// taxon B 1111011111011 ├─► P1 → {A, B, D}
// taxon D 1111011111011 ┘
// taxon C 1101110111111 ───► P2 → {C}
```

For each unique pattern: one `QR(X_pattern)`, reused for all taxa and all
iterations and all inference refits.

### 4.3 Iterative fixed-effect MLE

```
theta <- 0
repeat {
    beta  <- fit(Y_c - theta, X)          # per pattern, via cached QR
    fitted <- X %*% beta
    theta_new <- colMeans(Y_c - fitted)  # sample-specific bias estimate
    eps  <- sqrt(||beta - beta_prev||_F^2 + ||theta_new - theta||_2^2)
    if (eps <= tol) break
    theta <- theta_new
    if (iter == max_iter) break
}
```

Defaults: `tol = 0.01`, `max_iter = 20`. Residuals: `epsilon = Y_c - theta - fitted`.

**Non-negotiable:** QR decompositions are computed **once per pattern**, not per
iteration, not per taxon, not per sensitivity refit. This alone should beat the
already-optimised R path.

### 4.4 Sandwich variance (HC0)

```
V = (X'X)^{-1} ( sum_i eps_i^2 x_i x_i' ) (X'X)^{-1}
```

Precompute `XX[i] = vec(x_i x_i')` for all `i` (p(p+1)/2 numbers). Then for each
taxon accumulate `sum_i eps_ij^2 XX[i]`. Parallelise **over taxa**; accumulate
`eps^2` into a scratch `n_samples` vector in cache, then one pass over samples.

Block over taxa so the working set (`eps^2` for a block) stays in L2.

> **Compatibility hazard (must be encoded, not "fixed").** The oracle replaces
> certain `NA` entries in the sandwich accumulation with the constant `0.1`; an
> open upstream issue questions whether this is intentional. Therefore:
>
> ```rust
> pub enum CompatMode { Ancombc2_15, StrictSpec }
> ```
> `Ancombc2_15` reproduces the quirk. `StrictSpec` implements the mathematically
> intended behaviour. **Never silently diverge.** Both are tested; the default is
> the compatibility mode; the difference is documented in
> `docs/reference_behavior.md`.

### 4.5 EM bias estimation

Per fixed-effect coefficient, across taxa, ANCOM-BC2 fits a three-component
Gaussian mixture (negative-effect / null / positive-effect taxa). Initialisation
is fixed by the reference:

```
pi <- c(0.75, 0.125, 0.125)
```

EM iterates over `(pi, delta, l1, l2, kappa1, kappa2)`, with the variance
parameters optimised by **Nelder–Mead** inside the loop.

```rust
pub fn estimate_bias_em(beta: &[f64], var: &[f64], tol: f64, max_iter: usize)
    -> Result<BiasResult>;

pub struct BiasResult {
    pub delta_em: f64,
    pub delta_wls: f64,
    pub var_delta: f64,
    pub params: MixtureParams,
    pub iterations: usize,
}
```

Each coefficient is an independent task → parallelise over coefficients
(outer level of the Rayon plan, §7).

`delta_wls` is the weighted-least-squares analogue computed from the same inputs.

### 4.6 Bias correction, sampling fractions, regularization

```
beta_corr[j,k] = beta_star[j,k] - delta[k]     # except intercept / protected cols
samp_frac[i]   = mean_j(residual[i,j])         # sample-specific bias estimate
s0[k]          = quantile(se[,k], probs = s0_perc)     # default 0.05
W[j,k]         = beta_corr[j,k] / (se[j,k] + s0[k])
p[j,k]         = 2 * pnorm(-|W[j,k]|)
q[j,k]         = p.adjust(p[,k], method = p_adjust)
diff_abn[j,k]  = q[j,k] < alpha
```

**Quantile semantics are a parity trap.** R's `quantile()` defaults to
`type = 7`. Implement `type = 7` exactly and unit-test it against R for n = 1..200
and every `probs` used.

### 4.7 Global test

`H0: beta_1 = ... = beta_g = 0` as a quadratic Wald statistic with asymptotic
chi-square. Implement a **general contrast engine**, not a special case:

```rust
pub struct Contrast { pub matrix: DMatrix }

est    = A * beta
var    = A * Sigma * A'
stat   = est' * (var)^{-1} * est
p      = chisq_upper_tail(stat, df = rank(A))
```

Singular `var` must be handled by pseudo-inverse or explicit rank truncation —
match the oracle's choice, and record it in `docs/reference_behavior.md`.

### 4.8 Pairwise tests + mdFDR

Pairwise comparison of group coefficients **with direction**. ANCOM-BC2 controls
the *mixed directional* FDR across (a) multiple testing, (b) multiple group
comparisons and (c) sign/directional error. The published procedure first screens
taxa with the global test, then applies a family-wise correction to pairwise
comparisons for taxa that pass the screen.

> **Implement literally.** Do **not** "simplify" this into *run all pairs + BH*,
> because that would no longer be ANCOM-BC2.

### 4.9 Structural zeros

Computed after the core estimator (cheap, but statistically important; do not let it
distract the early numerical work). Test matrix must include: absent-in-group-A,
rare-in-group-A, present-everywhere, `NA` observations, small groups, both
`neg_lb` settings, and comparison of the **entire** `zero_ind` matrix.

### 4.10 Pseudo-count sensitivity analysis

**Conservative (`conservative = TRUE`)** — 4 fully independent runs:

```
counts ─┬─ pseudo 0.0 ─ full ANCOM-BC2 ─┐
        ├─ pseudo 0.1 ─ full ANCOM-BC2 ─┤
        ├─ pseudo 0.5 ─ full ANCOM-BC2 ─┼─► sensitivity score (ss_tab)
        └─ pseudo 1.0 ─ full ANCOM-BC2 ─┘
```

Embarrassingly parallel: `pseudo_values.par_iter().map(run_core)`.

**Non-conservative (`conservative = FALSE`)** — sampling fractions estimated once
at `pseudo = 0`, then 50 inference refits with `pseudo = 0.01 … 0.50`.
Also embarrassingly parallel; 50 jobs onto the shared Rayon pool.

---

## 5. Test strategy (rigorous, four independent layers)

Parity requires more than golden tests. Four orthogonal layers, each able to fail
independently.

### 5.1 Layer 1 — Golden parity tests (R as oracle)

Every quantity in §3.2, every fixture, tolerance per §5.6.

```
for fixture in small_ci_fixtures:
    run Rust   → actual
    run R 2.15.2 → expected (pre-generated .rds → .f64 blobs)
    compare per-quantity at its tolerance
    on failure: emit FIRST-DIVERGING-QUANTITY report with fixture id,
                coefficient index, taxon index, absolute + relative delta
```

The failure report is part of the deliverable — it is what makes the suite usable.

### 5.2 Layer 2 — Property / metamorphic tests (Rust `proptest`)

Golden tests say *"we copied R"*. Property tests say *"the implementation is
mathematically sane"*. Neither substitutes for the other.

| # | Property | Expectation |
|---|---|---|
| P1 | Permute taxa | Identical results, reordered |
| P2 | Permute samples **with** metadata | Identical statistics (up to fp reassociation) |
| P3 | Scale all counts by constant `c > 0` | Documented invariance (`pseudo=0` ⇒ exact invariance) |
| P4 | Change reference factor level | `beta`/`se`/`p` transform by the contrast map |
| P5 | Single group | Clean, typed error (not a panic) |
| P6 | Rank-deficient `X` | Deterministic, typed error naming the deficient columns |
| P7 | All-zero taxon | Filtered (or flagged as structural zero) per config |
| P8 | Constant taxon (all counts equal) | `beta == 0` for non-intercept coefficients |
| P9 | Duplicate every sample twice | `se` shrinks by `sqrt(2)`; `p` unchanged to fp |
| P10 | Add a pure-noise covariate | Existing coefficients unchanged (orthogonal design) |
| P11 | `p_adjust = "none"` | `q == p` bitwise |
| P12 | `n = 1` sample | Typed error |
| P13 | Idempotence | `run(run(x)) == run(x)` bitwise |
| P14 | `pseudo → ∞` | All `Y` equal ⇒ `beta == 0` |
| P15 | Thread-count invariance | 1 thread == 8 threads, bitwise-or-within-tol |
| P16 | `CompatMode` separation | `Ancombc2_15 ≠ StrictSpec` **only** on fixtures that exercise the quirk |

P15 is non-negotiable: any reduction-order nondeterminism in a parallel reduction
is a numerical-parity bug that will surface as nondeterministic p-values.

### 5.3 Layer 3 — Statistical validation (simulation)

Numerical parity is necessary but not sufficient: it proves the code matches R, not
that both are *correct*. Rebuild the paper's simulation design (simulated
microbiome data with real microbial mean/covariance structure) and extend it.

**Design grid** (each cell 100–1000 reps):

| Factor | Levels |
|---|---|
| DA proportion | 0, 5, 10, 30, 50, 70, 90 % |
| samples | 20, 50, 100, 500 |
| taxa | 100, 500, 1,000, 5,000 |
| zero inflation | 0, 10, 50, 90 % |
| effect size | log FC ∈ {0.5, 1, 2, 4} |
| library size | mean ∈ {1e3, 1e4, 1e5}, CV ∈ {0, 0.3, 0.6} |
| confound | sampling-fraction–DA confounded / not |

**Reported per cell:** empirical FDR (target ≤ nominal), power, sign concordance,
LFC bias vs truth, SE calibration ratio (`mean(se)/sd(beta_hat)` ≈ 1),
F1/Jaccard of DA calls.

**Acceptance:** for every cell, empirical FDR of Rust within Monte-Carlo error of
R's (3 SE of the binomial), and power within 3 SE. Report any cell where they
disagree as a **finding**, not a test failure — a genuine divergence is
publishable.

### 5.4 Layer 4 — Real-data validation

Three datasets at minimum:

- **A:** ANCOM-BC2 vignette dataset (sanity / regression-to-reference)
- **B:** QMP-like moderate dataset (performance + memory claim)
- **C:** large modern metagenomic set (tractability claim)

Metrics: Spearman + Pearson of `beta`, `max|Δbeta|`, `max|Δse|`, p-value and q-value
correlation, significant-call Jaccard, `diff_robust`/`passed_ss` agreement.

Acceptance: for a correct port the scatter is visually indistinguishable from
`y = x`, and **`diff_abn` agreement ≥ 99.99 %** (Level D, §5.6).

### 5.5 Structural-zero and edge-case matrix

| Case | Expected |
|---|---|
| Completely absent in group A | flagged, `neg_lb` semantics honoured |
| Rare in group A (< 1 obs) | flagged |
| Present in every group | not flagged |
| `NA` in counts | exact match to R's `NA` handling (incl. the sandwich quirk) |
| Group of size 1 | pairwise defined? match R (possibly dropped) |
| `keep_zero = TRUE` | included in primary fit, `zero_ind` still reported |
| `perc_thres` exactly at boundary | exact `<` vs `<=` match |

### 5.6 Numerical contract (tolerances)

`docs/numerical_contract.md`, enforced in one shared `assert_close!` used by all
suites.

**Level A — exact (bitwise or `< 1e-15`)**
taxa retained · samples retained · column names · structural-zero flags ·
`diff_abn` · `passed_ss` · `diff_robust` · pattern assignment

**Level B — floating point**
`beta` rtol `1e-8` · `theta` rtol `1e-8` · `se` rtol `1e-7` · `delta_em` rtol
`1e-7` · `delta_wls` rtol `1e-7` · `vcov` rtol `1e-7` · `s0` rtol `1e-9`

**Level C — probabilities**
`p` atol `1e-10` · `q` atol `1e-10`

**Level D — inferential**
`diff_abn` concordance > 99.99 % over 10,000-taxon runs; **100 %** on all
small deterministic fixtures.

> Rule: **never loosen a tolerance to make a test pass.** If a test fails, first
> determine whether the difference is (a) a bug, (b) a legitimate
> compiler/BLAS-order difference, or (c) an optimizer-path difference. Only
> category (b)/(c) justifies a tolerance change, and it must be recorded in the
> commit message and in `docs/reference_behavior.md`.

### 5.7 Golden fixture matrix

| Axis | Levels |
|---|---|
| shape | 10×10, 100×30, 1,000×100, 10,000×500 |
| sparsity | zero-free, 10 %, 50 %, 90 % zeros |
| predictor | continuous, binary, 3 groups, 5 groups |
| covariates | 1, 5, 10, interaction term |
| balance | balanced, unbalanced |
| structural zeros | present, absent |
| pseudo | 0, 0.1, 0.5, 1.0 |
| sensitivity | on/off × conservative on/off |
| adj. method | all 7 |
| group | NULL, 2-level, 5-level |

**CI policy:** the 10×10 / 100×30 / 1,000×100 fixtures run on **every** push
(target < 6 min). The 10,000×500 fixtures plus all simulation cells run
**nightly**. Full simulation grid runs **weekly**. A dedicated `parity` label runs
a hand-picked subset on every PR.

### 5.8 Differential testing (optional, high value)

For a few fixtures, also run the **pre-optimisation** ANCOMBC R code path (an older
tag, or the same source with pattern-grouping disabled) to prove the optimisation
itself introduced no drift. This separates "port error" from "upstream change".

---

## 6. Benchmarking design

### 6.1 Rules (to prevent the usual dishonesty)

1. **No "80× faster on toy input" claims.** Benchmark a scaling *surface*.
2. Report the R version **you** were competing against, including commit SHA.
3. Report R **with its own parallelism enabled** as a separate arm. Beating
   single-threaded R is easy and uninteresting; beating `parallel::mclapply` is not.
4. Never compare against a stale baseline. If upstream optimised, say so in the
   paper.
5. Report the *profiling outcome* (see §6.6) even if it is unflattering.

### 6.2 Scaling surface

| Run | Samples | Taxa | Covariates | Groups | Sensitivity |
|---|---|---|---|---|---|
| S1 | 50 | 500 | 3 | 2 | off |
| S2 | 100 | 1,000 | 5 | 3 | off |
| S3 | 500 | 5,000 | 5 | 3 | off |
| S4 | 1,000 | 10,000 | 10 | 5 | off |
| S5 | 5,000 | 20,000 | 10 | 5 | **on** |
| S6 | 1,000 | 10,000 | 10 | 5 | **on** (non-cons.) |

### 6.3 Arms

```
R-1        ANCOMBC 2.15.2, 1 core, no parallel sections
R-par      ANCOMBC 2.15.2, all available parallelism (mc.cores = hw)
Rs-1       ancombc2-rs, RAYON_NUM_THREADS=1
Rs-4       ancombc2-rs, 4 threads
Rs-8       ancombc2-rs, 8 threads
Rs-16      ancombc2-rs, 16 threads
Rs-32      ancombc2-rs, 32 threads
```

### 6.4 Metrics

- wall time, CPU time (user + sys)
- **peak RSS** (`/usr/bin/time -v`, `getrusage`)
- bytes allocated, allocation count (`dhat` / `jemalloc` stats)
- throughput: taxa·samples / second
- **stage-resolved** timings: preprocess, MLE, sandwich, EM, correction, tests,
  sensitivity, (de)serialisation

> Stage-resolved timings are mandatory. Without them you cannot tell a real
> algorithmic win from "R was slow at copying dataframes".

### 6.5 Benchmark harness

`benchmarks/run_bench.sh` drives one container image with pinned R + ANCOMBC +
toolchain; every run appends a row to `benchmarks/results/results.jsonl` with
`{arm, run_id, threads, hw, date, stage_times, wall, cpu, rss}`. Results are
committed so the paper's tables are reproducible from git history.

### 6.6 Pre-benchmark profiling gate (do this BEFORE writing kernels)

`scripts/profile_r.R` under `Rprof`, broken down by section:

```
iterative regression      ?
sandwich covariance       ?
EM                        ?
lmer / random effects     ?
R object handling         ?
BLAS calls                ?
```

Expected outcomes:

- **Outcome A** — computation-dominated, spread across kernels → proceed;
  expected large win.
- **Outcome B** — 75 %+ in BLAS calls → a pure Rust rewrite may yield only
  ~2–3×. Still publishable if statistical equivalence is airtight, but the paper's
  claim must change.
- **Outcome C** — >90 % in `lmer`/random effects → the *real* project is
  accelerating repeated-measure ANCOM-BC2, and the fixed-effects port becomes the
  correctness foundation.

> **This gate can save months.** Run it in week 1. Report the outcome in the
> README, whatever it is.

### 6.7 Performance gates (continuation criteria, not claims)

| Gate | Criterion | If missed |
|---|---|---|
| **P1** | fixed-effect numerical kernel ≥ 3× faster than `R-1` (or than `R-par`, if `R-par` is the honest baseline) on S3 | profile before adding features |
| **P2** | full fixed-effects end-to-end ≥ 2× vs same arm | identify the un-accelerated stage |
| **P3** | sensitivity analysis ≥ 5× on 8–16 cores (S5/S6) | revisit the parallel plan |
| **P4** | peak RSS ≤ 60–70 % of R peak on S4/S5 | revisit block/buffer strategy |
| **P5** | strong scaling efficiency ≥ 0.7 from 1→16 threads | fix reduction order / false sharing |

If P1 is missed but P2–P5 pass, **continue** — the port is still valuable. If P1
*and* P2 are missed, re-read Outcome B above and re-scope.

---

## 7. Parallelisation plan

**One** global Rayon pool. No nested pools. Explicit nesting policy:

```
Level 1  pseudo-count runs (4 or 50)          ← outermost
Level 2  EM bias estimation, per coefficient
Level 3  missingness-pattern groups
Level 4  taxa within sandwich covariance
```

Budget: with `T` hardware threads and a pseudo-count run count `R`, use
`T / R` threads per run (e.g. 32 hw threads, 4 runs → 8 each), **not** 4 × 32.
Expose this as `ancombc2_core::set_nested_threads(r)`.

**Determinism requirement.** Every parallel reduction must use a deterministic,
index-ordered accumulation (per-thread partial buffers merged in fixed index
order, never `par_iter().sum()` over floats without a fixed-order fold). Assert
P15 in CI.

---

## 8. Memory design

Target: never materialise six `taxa × samples` matrices at once.

- Reusable `Workspace { centered, fitted, residuals, beta, theta, squared_resid }`
  with in-place ops (`residuals = centered - theta - fitted`, `squared_resid[i] *= ...`).
- Sandwich covariance accumulates over a **block** of taxa sized to fit
  `n_samples × f64` in L2; the per-block `eps²` vector is reused.
- `XX[i] = vec(x_i x_i')` precomputed once: `n_samples × p(p+1)/2` doubles.
- Structural-zero `zero_ind` as a bitset, not a `bool` matrix.
- Stream outputs: write `res_global` / `res_pair` / `ss_tab` incrementally rather
  than concatenating.
- Bench RSS and allocation count, not just wall time (§6.4).

---

## 9. Dense vs sparse

Start **dense `f64` + compact missingness bitsets**. Rationale: the hot loops are
regressions on row subsets with all-numeric dense kernels, and a generic sparse
layer typically makes that slower.

Revisit only after profiling shows the dense path is memory-bandwidth-bound at
extreme sparsity (> 90 % zeros, S5/S6). Provide `Representation::{Dense, SparseTaxa}`
behind a config flag and benchmark both.

Measured 2026-10-05, on `bm3` (500 taxa x 5000 samples, 50 % zeros — the only
benchmark dataset above the ~33 % crossover where the sparse layout can in
principle win), one thread, three runs each:

| representation | wall | peak RSS | results |
| --- | --- | --- | --- |
| dense (default) | 1.20 s | ~272 MB | — |
| `--sparse-taxa` | 1.23 s | ~270 MB | byte-identical (`res`, `res_global`, `res_pair`) |

No measurable difference on either axis. The flag's own help text already warns
that only the prevalence and library-size screens benefit and only the storage;
on a 1.2 s run allocating ~2.0 GB total, the screens are too small a fraction
for the representation to show. The > 90 % condition for revisiting is met by no
dataset here, so the 42-arm surface stays dense-only by design rather than by
omission — and this paragraph is the record the code comment at
`preprocess.rs` asked for instead of an argument.

---

## 10. Release roadmap

| Release | Contents | Gate to ship |
|---|---|---|
| v0.1 | fixed-effect core, `pseudo_sens = FALSE`, `group = NULL` | P1 met; Levels A–C pass on 1,000×100 |
| v0.2 | structural zeros, prevalence filtering, primary inference | Level A exact on zero fixtures |
| v0.3 | global tests | parity on 2-, 3-, 5-group fixtures |
| v0.4 | pairwise tests + mdFDR | parity on unbalanced-group fixtures |
| v0.5 | conservative sensitivity analysis | exact `ss_tab`, `passed_ss`, `diff_robust` |
| v0.6 | non-conservative sensitivity analysis | parity on 50-refit fixtures |
| v0.7 | interaction compatibility | parity on interaction fixtures |
| v0.8 | full R wrapper (phyloseq + TSE + formula) | R-side testthat suite green |
| v0.9 | large-scale benchmarking + real-data validation | S1–S6 complete, P1–P5 reported |
| **v1.0** | **stable fixed-effects ANCOM-BC2-compatible implementation** | **Levels A–D pass; simulation FDR/power parity; all benchmarks published** |
| v1.1+ | Dunnett | |
| v1.2+ | trend (ordered-pattern) testing | |
| v2.0 | random intercept | |
| v2.1 | random slopes | |
| v2.x | general mixed models | |

Usable releases throughout, so there is always a portfolio artifact and a
citable version.

---

## 11. Compatibility vs correctness policy

Two documents, deliberately separate:

- **`docs/reference_behavior.md`** — what the oracle *does*, including quirks
  (e.g. the `NA → 0.1` substitution in the sandwich accumulation; whatever
  quantile/p-adjust edge cases the oracle exhibits). Each entry: description,
  locus in R source, upstream issue link if any, our reproduction, and status.
- **`docs/statistical_spec.md`** — what the *paper* says the method should do.

```rust
pub enum CompatMode {
    Ancombc2_15,   // default: reproduce the oracle, quirks included
    StrictSpec,    // mathematically intended behaviour
}
```

Rules:
1. Default is `Ancombc2_15`.
2. The two modes must be tested to differ **only** on fixtures that exercise a
   documented quirk.
3. Any divergence gets a `docs/reference_behavior.md` entry the same day it is
   discovered.
4. Never silently diverge. This discipline is what makes the port trustworthy.

---

## 12. Upstream engagement

Sequence: parity achieved + benchmark in hand → **then** contact the maintainer
(Huang Lin). Not before. The message carries numbers, not enthusiasm:

> We built an independent Rust implementation of ANCOM-BC2's fixed-effect core.
> Across N datasets: beta agreement = …, SE agreement = …, DA-call concordance = …
> Runtime: ANCOMBC 2.15.2 = … s (32 cores), ancombc2-rs = … s (32 cores), peak RSS …
> We would value your read on several compatibility details, in particular the
> treatment of `NA` entries in the sandwich variance accumulation (upstream issue
> #…) and the default `s0` quantile.

The package is Artistic-2.0. Regardless of licence, **implement from the paper
and verified behaviour, not by transliterating R source.** This is both the
correct engineering practice and the right posture for a methods paper.

---

## 13. Publication plan

**Title angle:** *ancombc2-rs: a scalable, numerically validated implementation of
bias-corrected differential abundance analysis for large microbiome studies.*

Not "we rewrote ANCOM-BC2 in Rust" — that is not a contribution.

**Four claims, each with a figure and a table:**

1. **Numerical equivalence** — Rust reproduces reference ANCOM-BC2 inference.
2. **Statistical equivalence** — FDR and power preserved across simulations.
3. **Computational scaling** — runtime and memory reduced; scaling laws
   characterised.
4. **Large-scale feasibility** — analyses previously impractical become tractable.

**Figures**

| # | Content |
|---|---|
| 1 | ANCOM-BC2 algorithm + `ancombc2-rs` architecture (two-panel) |
| 2 | R vs Rust parity: `beta`, `se`, `p` hexbin/scatter + parity band |
| 3 | runtime vs taxa |
| 4 | runtime vs samples |
| 5 | memory scaling |
| 6 | FDR and power simulation parity (grid heatmaps) |
| 7 | headline large real-world analysis (e.g. 20,000 × 5,000 × 8 + sensitivity) |

**The compelling demonstration:** 20,000 taxa × 5,000 samples × 8 covariates with
sensitivity analysis, where R becomes difficult and Rust remains usable — with
Figure 2 proving the answer is identical where both can run.

**Submission:** *Bioinformatics* / *Briefings in Bioinformatics* / *Genome Biology*
(Methods) / *mSystems*. Preprint on bioRxiv **after** the parity suite is public,
so reviewers can run it. Open-source the fixture generator and all validation
scripts — reproducibility is the review criterion.

---

## 14. CI matrix

```yaml
# .github/workflows/ci.yml
on: [push, pull_request]

jobs:
  fmt-clippy:      # fmt --check, clippy -D warnings, miri smoke on numerics
  unit:            # cargo test --workspace  (Layers 1 small, 2 all)
  parity-small:    # all 10x10 / 100x30 / 1,000x100 fixtures, Levels A–D
  determinism:     # threads=1 vs 8, byte-identical outputs (P15)
  nightly:
    - parity-large  # 10,000x500 fixtures
    - simulation    # full grid, 100 reps/cell
    - benchmarks    # S1–S4, results.jsonl diff
  weekly:
    - simulation-full (1000 reps/cell)
    - benchmarks-full (S1–S6)
```

`R-parity` job re-runs the R harness weekly and **fails** if the newly generated
goldens differ from the committed ones — this catches upstream drift in ANCOMBC
and forces an explicit, reviewed version bump rather than silent divergence.

---

## 15. First four weeks (concrete)

**Week 1 — reference + profiling**

- Pin ANCOMBC 2.15.2 @ SHA; record R, BLAS, RNG in `reference/env/`.
- Read and annotate: `ancombc2()`, `.iter_mle()`, `.lm_fit_all()`,
  `.sandwich_vcov()`, `.bias_em()`, `.struct_zeros()`, `.pairwise_wilcox()`/`.pw()`.
- Fork + instrument a copy of the R source to dump intermediates 4–15 of §3.2.
- Generate the full golden set for the small fixtures.
- Run `scripts/profile_r.R`; **record Outcome A/B/C in the README.**
- Write `docs/algorithm.md`, `docs/numerical_contract.md`,
  `docs/reference_behavior.md`.
- Stand up `benchmarks/` with one working S2 run per arm.

**Week 2 — preprocessing + missingness + QR**

- `CountMatrix`, filters, `log`/`center`, structural zeros.
- `MissingnessPattern` bitsets + grouping.
- Cached per-pattern `QR`.
- **Gate 1:** taxa/samples/structural zeros exact; log counts < 1e-14; centered
  < 1e-12; `beta*` rtol 1e-8. **Do not proceed until green.**

**Week 3 — iteration + sandwich**

- `theta` loop, residuals, convergence trace matching the oracle.
- Blockwise sandwich with the `NA → 0.1` compatibility quirk behind
  `CompatMode`.
- **Gate 2:** `beta*`/`theta`/residuals/`vcov`/convergence all at Level B.

**Week 4 — EM + correction + inference**

- 3-component mixture EM with Nelder-Mead variance optimisation; `delta_wls`.
- Bias correction, `samp_frac`, `s0` (R type-7 quantile), `W`, `p`, `q`, `diff_abn`.
- **Gate 3:** `delta_em`/`delta_wls`/`var_delta` rtol 1e-8 across hundreds of
  generated `beta`/variance vectors; `p`/`q` atol 1e-10.

**End of week 4 (Milestone A):**

```
counts ──► Rust ──► primary ANCOM-BC2 table
   rand_formula = NULL
   pseudo_sens  = FALSE
   group        = NULL
```

with strong numerical agreement. **Then stop adding features and benchmark**
(§6.6, gates P1–P5). That is the decision point for the whole project.

**Weeks 5–12:** global tests → pairwise + mdFDR → structural-zero edge cases →
conservative sensitivity → non-conservative sensitivity → R wrapper →
benchmarks + real-data validation + simulations → paper.

---

## 16. Risk register

| # | Risk | Likelihood | Mitigation |
|---|---|---|---|
| R1 | EM/Nelder-Mead optimiser path diverges at 1e-10 | **High** | Mirror the reference's simplex, ordering, and stopping rule exactly; if impossible, document the divergence class and set a justified tolerance from first principles |
| R2 | R `quantile` `type = 7` mismatch | Medium | Dedicated exhaustive test n = 1..200 |
| R3 | `NA → 0.1` sandwich quirk misread | Medium | Encode as `CompatMode`; two tests |
| R4 | Determinism lost under parallelism | **High** | Fixed-order reductions; P15 in CI from week 1 |
| R5 | Upstream changes the oracle mid-project | Medium | SHA-pinned goldens + weekly drift job + explicit version bump |
| R6 | Speedup is only ~2× (Outcome B) | **Medium** | Decide at week-4 gate; reposition paper as "validated, faster, memory-lean" or pivot to mixed models |
| R7 | R object overhead dominates, so Rust gains little | Medium | Same; report honestly |
| R8 | mdFDR spec under-documented | Medium | Read the paper's method section line by line; implement literally; unit-test each step; ask the maintainer |
| R9 | Simulation infrastructure becomes a second project | Medium | Reuse existing DA-benchmarking expertise; timebox to 3 weeks |
| R10 | "Just another port" framing in review | Medium | Lead with statistical validation + large-scale feasibility, not language |
| R11 | Scope creep into `lmer` | **High** | Roadmap is explicit; v2.0 only; `random effects` is a *separate* RFC |
| R12 | Finding a real upstream bug and being asked to fix R too | Medium | Report it; ship `StrictSpec`; offer a minimal upstream PR separately |

---

## 17. Definition of done (v1.0)

- [ ] Levels A–D pass on every fixture in §5.7 (small in CI, large nightly).
- [ ] All property tests P1–P16 pass, including thread-invariance.
- [ ] Simulation grid: empirical FDR/power within Monte-Carlo error of R in every
  cell; any divergence investigated and written up.
- [ ] Real datasets A/B/C: `diff_abn` concordance ≥ 99.99 %; scatter visually on
  `y = x`.
- [ ] Benchmarks S1–S6 across all 7 arms, stage-resolved, committed to git; gates
  P1–P5 reported (pass *or* explicitly not).
- [ ] `docs/` complete: algorithm, numerical contract, reference behavior,
  statistical spec, compatibility, reproduction.
- [ ] R wrapper passes its own testthat suite; `ancombc2(object, ...)` accepts
  phyloseq/TSE objects and returns the documented tibbles.
- [ ] README states the compatibility target string and the profiling outcome.
- [ ] No un-documented divergence from the oracle exists.

---

## 18. One-paragraph summary

Freeze the R oracle at an exact commit and capture **every intermediate quantity**,
not just p-values. Build the Rust port in dependency order — preprocessing,
missingness patterns, cached-QR iterative MLE, sandwich variance, EM bias
correction, regularization, inference — validating each stage against its own
golden file before moving on. Keep `core` free of R and pair it with an R wrapper
that owns `model.matrix()` and contrasts. Prove correctness with four independent
test layers: golden parity, property/metamorphic tests, FDR/power simulation, and
real-data concordance, all under an explicit tolerance contract that is never
loosened to make a test pass. Profile the R original **before** optimising, so the
performance story matches reality. Benchmark a scaling surface against both
single-threaded and parallel R, report stage-resolved times, peak RSS and
allocations, and hold to explicit gates. Only then write the paper — as numerical
equivalence, statistical equivalence, computational scaling, and large-scale
feasibility, not as a language port.
