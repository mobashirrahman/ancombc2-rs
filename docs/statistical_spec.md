# Statistical specification

What the estimators are, what they assume, and how each assumption is tested.
The reference is ANCOMBC 2.15.2 at `dc4febdf5`; where the reference's behaviour
differs from the mathematics, this document says so and points at
`docs/reference_behavior.md`.

## 1. What is being estimated

For taxon *i* and sample *j*:

```
log(count_ij + pseudo) - mean_i = x_j' beta_i + s_j + e_ij
```

* `beta_i` is the taxon's covariate effect. The group coefficients are the
  contrasts of interest.
* `s_j` is the sample's **sampling fraction**: a per-sample offset that
  confounds everything if ignored, because an unbalanced treatment changes which
  samples are sequenced deeply.
* `e_ij` is the residual. Its variance is *not* Poisson — counts are
  overdispersed — and modelling it as such is the second thing ANCOM-BC2 fixes.

The per-taxon centring (`- mean_i`) removes each taxon's overall abundance, so
`beta_i` is a contrast rather than a level and taxa of wildly different magnitude
are on a common scale.

## 2. Estimation

`beta_i` and `s` are estimated jointly by alternating least squares: fit `beta`
given `s`, then set `s` to the per-sample mean residual. The objective is not
concave jointly, so the iteration is a coordinate ascent that is not guaranteed to
find the global optimum — which is why the reference fixes `max_iter = 20` and a
loose `tol = 0.01` rather than iterating to convergence. The iteration is stopped
by an *absolute* epsilon on the sum of squares over all taxa, so for any realistic
table it always runs to the cap; a run that converged early would be a
*different* estimator, not a more accurate one.

`s` is estimated over the bias set `O1`, which includes the taxa with structural
zeros. Those taxa carry no information about their own group's effect, but they
do carry information about the sequencing depth of each sample, and excluding them
throws that away.

## 3. The variance, and why it is a sandwich

With a per-sample offset `s` estimated from the same data, the model errors within
a sample are correlated, so the OLS covariance `sigma^2 (X'X)^-` is wrong in a way
that depends on the design. The HC0 sandwich

```
V_i = (X'X)^- (sum_j eps_ij^2 x_j x_j') (X'X)^-
```

is consistent under arbitrary within-sample correlation, at the cost of being
conservative. It is the right default for count data with a nuisance offset.

**The reference's quirk.** A missing entry of the outer-product matrix is replaced
by `0.1` rather than the term being skipped, which inflates the variance slightly
for a design with missing entries. `CompatMode::Ancombc2_15` reproduces it;
`CompatMode::StrictSpec` does not. This is a genuine divergence from the intended
estimator, and the two modes are tested to differ *only* where it applies
(property test P16).

## 4. The bias correction

`beta_star` is a consistent estimator of `beta` but a biased estimate of any
*individual* taxon's contrast, because the per-sample offset is estimated from all
taxa and each taxon contributes to it. ANCOM-BC2 subtracts that bias.

The bias is modelled as a three-component Gaussian mixture over taxa, per
coefficient: a dominant component for the bulk, and two for the tails. The
posterior probability of each component gives the weight `delta_wls` and its
variance gives `var_delta`; the E-M output `delta_em` is the correction actually
subtracted. The correction's own variance is added to the coefficient's variance —
`var + var_delta + 2 sqrt(var * var_delta)`, which is the variance of a sum with a
*negative* covariance between them.

The mixture is a model, and a crude one: three normal components cannot describe a
heavy-tailed taxon distribution. Its value is that it is a smooth, symmetric
function of the taxon set, so it corrects the bulk without assuming a
distributional form for the tails. The sensitivity analysis exists because the
correction depends on which taxa are in the set.

## 5. `s0`, and why a variance is added

`se = sqrt(var_hat + s0)` with `s0` the `s0_perc` quantile (default 5th percentile,
per coefficient) of the variance distribution across taxa. This is the SAM
regulariser: without it, a taxon whose variance estimate is near zero gets an
enormous `W` and dominates the downstream mixture fit and the p-value adjustment.
The 5th percentile is a robust estimate of the variance floor, and it is taken
*per coefficient* because the intercept's variance is on a different scale from a
group contrast's.

## 6. Inference

`W = beta / se`, `p = 2 * pt(|W|, dof)`, `q = p.adjust(p)` per coefficient, and
`diff_abn = q <= alpha`.

`dof = n_used - rank` for the second MLE — the residual degrees of freedom of the
*per-taxon* regression, not of the whole table. A taxon observed on few samples
therefore has a small `dof` and a heavier-tailed p-value, which is the correct
behaviour. In the first MLE `dof` is `NULL`, and R's `pt(x, df = NULL)` evaluates
`df = Inf`, i.e. the normal limit; that path is not used for inference, only for
the variance, but it is reproduced because the reference computes it.

`p.adjust` over all seven methods is reproduced exactly, including the details that
matter and are easy to miss: `NA` values are dropped before ranking, `n` is a
double in the pairwise test, and `hommel` with `n == 2` is silently `hochberg`.
See `docs/reference_behavior.md` §6.

## 7. The multi-group tests

* **Global**: a quadratic Wald test of `H0: beta_group = 0`, using the
  off-diagonal terms of the sandwich block. Under the sandwich the statistic is
  asymptotically chi-square on the number of group coefficients, which is why it
  needs the full block and not just the diagonal.
* **Pairwise**: the published mixed-directional FDR, implemented as
  screen-then-family-wise rather than simplified to BH. The screen is the global
  test at `BH` — hard-coded in the reference, not the caller's method. The
  adjustment is *within each taxon* across its `C(n_group, 2)` contrasts, at
  `n = n_col * n_tax / R`, where `R` is the number of taxa the screen rejected.
  That `n` is deliberately inflated: the screen has already spent the
  multiplicity, and the within-taxon adjustment covers the rest.

  The division by `R` is unguarded in the reference, so `R == 0` makes every taxon
  significant. Reproduced and flagged.

## 8. Sensitivity to the pseudo-count

With `pseudo = 0`, a zero count becomes `NA` and the taxon is fitted on the samples
where it is observed. With `pseudo > 0`, no count is missing, and the fit uses
every sample with a slightly biased value. The two are different estimators, and
which one you get is a modelling choice, not a numerical detail.

`pseudo_sens` measures the disagreement: the proportion of pseudo-count runs whose
adjusted p-value exceeds `alpha`, and whether that agrees with the main run. A
taxon is `diff_robust` when it is significant under the main run and the
sensitivity runs agree.

The conservative and non-conservative variants make different trade-offs, and the
reference offers both: conservative reruns everything, which is expensive but
makes each run a complete analysis; non-conservative reuses `theta_hat` and refits
only the inference, which is 50 cheap regressions instead of 50 complete ones.
That is why non-conservative is 50x cheaper, and why it is not the same estimator.

## 9. Assumptions, and what happens when they fail

| assumption | consequence if violated | what the code does |
| --- | --- | --- |
| the sample offset is common to all taxa | the bias correction is misspecified | nothing detects this; the sensitivity analysis probes a different axis |
| residuals are conditionally heteroskedastic | HC0 is conservative, not anti-conservative | by design |
| the design is not rank deficient *per taxon* | the coefficient is not identified | a rank-deficient sub-design drops the aliased coefficient as `0`; see `docs/reference_behavior.md` §10 |
| every group has ≥ 2 samples | the variance is undefined | `GroupTooSmall`, matching the reference's error |
| ≥ 3 group levels for the multi-group tests | a global test has 1 df and is the t-test | the reference *disables* the comparisons with a warning; reproduced |
| counts are independent across samples | the sandwich is not consistent | not testable from the count table alone |

## 10. How each claim is tested

| claim | test | layer |
| --- | --- | --- |
| every quantity matches the oracle | golden parity, 4 fixtures, Levels A–D | 1 |
| taxon/sample permutation equivariance | P1 | 2 |
| sample permutation invariance | P2 | 2 |
| count scaling invariance | P3 | 2 |
| reference-level reparametrisation | P4 | 2 |
| typed errors for a single group, a rank-deficient design, one sample | P5, P6, P12 | 2 |
| an all-zero taxon is filtered, a constant taxon is not called significant | P7, P8 | 2 |
| duplicated samples halve the variance | P9 | 2 |
| an orthogonal covariate does not move the *fit* | P10 | 2 |
| `p_adjust = none` is the identity | P11 | 2 |
| idempotence, `CompatMode` separation, thread invariance | P13, P16, P15 | 2 |
| empirical FDR, power, LFC bias, SE calibration, F1 | the simulation grid | 3 |
| agreement on real data | three datasets, Spearman/Pearson of `beta`, p/q correlation, Jaccard of calls | 4 |

Layers 3 and 4 are specified in `PLAN.md` §5 and are **not yet executed**; see
`docs/compatibility.md` for the current state.
