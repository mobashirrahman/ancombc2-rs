# Reference behaviour

Compatibility target: **ancombc2-rs v0.1 is equivalent to ANCOMBC 2.15.2 at
commit `dc4febdf59badb3a8dfe0c767ef2186323c2199a`, seed 42.**

Every divergence below was found by the golden parity suite, not by reading the
R source and guessing. Each entry records the oracle's behaviour, the reason, and
where the divergence lives in the Rust code. Nothing in this file is a claim
about a benchmark or a test that was not executed.

## 1. Two taxon sets, and why the first MLE runs on the larger one

`ancombc2.R` filters the count table twice:

```r
core1 = .data_core(data = feature_table,         tax_keep = NULL,    ...)
O1     = core1$feature_table                    # prevalence + library size
core2 = .data_core(data = feature_table_aggregate, tax_keep = tax_keep, ...)
O2     = core2$feature_table                    # O1 minus the structural zeros
```

`O1` is a **superset** of `O2`: `core1` is called with `tax_keep = NULL`, so it
never sees the structural-zero screen. `.ancombc2_prep` then uses

| quantity | taxon set |
| --- | --- |
| `y1`, first MLE, `delta_em`/`delta_wls`/`var_delta`, `theta_hat` | `O1` |
| `y2`, second MLE, `beta`, `var_hat`, `se`, `p`, `q`, `diff_abn`, the whole result table | `O2` |

So a taxon with a structural zero still contributes to the estimate of the
sampling fractions, and the reported table is smaller than the table the bias
was estimated from. `CoreOutput` keeps both: `o1`/`taxa_bias` for the bias set
and `o2`/`taxa` for the reported set.

Getting this wrong is silent. With no structural zeros the two sets are equal,
so the three small fixtures pass either way. `fx04` is the first fixture with
structural zeros (9,800 retained of 10,000) and it fails immediately otherwise.

A second consequence: `.iter_mle` is called twice with two different responses,
and it recomputes the missingness grouping each time. The second MLE therefore
builds its own `DesignCache`; reusing the first one's cache would index the
wrong rows.

## 2. The sandwich replaces an `NA` entry of each term with 0.1

`.sandwich_vcov` accumulates the outer-product matrix a sample at a time and
substitutes `0.1` for every `NA` of each term:

```r
term_j = outer(eps2[idx, j], XX[j, ])
term_j[is.na(term_j)] = 0.1
sigma2_xxT = sigma2_xxT + term_j
```

The substitution is **per entry**, and the case that matters is not the one it
looks like. `eps2` is `NA` wherever a taxon is unobserved — its count was zero and
`pseudo` is zero — so for such a pair *every* entry of the term is `NA` and every
one becomes `0.1`. The taxon contributes a spurious constant `0.1` to all `p²`
positions of its variance for a sample it was never observed in. A *partially*
missing design row is the rarer case: only the positions involving the missing
design entry are `NA`, and only those are replaced.

So for any table with zeros and `pseudo = 0` the quirk fires almost everywhere,
and the two `CompatMode`s differ on almost any real input. That is the honest
description, and it is why P16 asserts that the modes differ *only* through this
one mechanism rather than that they coincide on clean data.

The mathematically intended behaviour is to skip the term. `CompatMode` carries
the choice:

* `Ancombc2_15` (default) reproduces the `0.1`;
* `StrictSpec` skips the term.

`CompatMode` reaches the sandwich through `AncombcConfig::compat` →
`iter_mle_estimate_theta` / `iter_mle_fixed_theta` → `sandwich_all`. Property test
P16 asserts the two modes differ only where the quirk applies, which is the
regression test for the wiring.

## 3. `.bias_em` uses the *previous* iteration's `l1`, `l2` and `delta`

`delta_em`, `delta_wls` and `var_delta` are read from a vector that the loop
updates with the *new* `l1`/`l2` while the Nelder-Mead objective for the
component variances is still evaluated with the *old* ones. Reproduced exactly;
see `em.rs`. A "cleaned up" version differs from the oracle at the fifth
significant figure, so the quirk is load-bearing.

The mixture is initialised at `pi = (0.75, 0.125, 0.125)`, as published.

## 4. The E-M variance optimisation is Nelder-Mead, not a closed form

The reference calls `NLopt` with the Nelder-Mead algorithm, default tolerances
and R's RNG. Its effective initial simplex is empirically `0.75 * |x0|`, not
NLopt's documented `0.05`; `stats/src/nelder_mead.rs` uses `0.75` and the
comment there records the measurement. A fixture with more taxa converges
differently under a different starting simplex, which is why the E-M scalars are
compared at `TOL_BETA` and not tighter.

## 5. Primary p-values use Student-t tails, not normal tails

`2 * pt(abs(W), df = dof, lower.tail = FALSE)`. In the fixed-effects path `dof`
is `n_used - rank` for the *second* MLE and `NULL` for the first; R's
`pt(x, df = NULL)` evaluates the default `df = Inf`, which is the normal tail, not
`NA`. Both branches are reproduced.

`dof` is one value per taxon and is expanded to one per coefficient *after*
estimation, by `matrix(rep(dof, n_fix_eff), ncol = n_fix_eff)`. Doing the
expansion per taxon instead is the same vector, but the pairwise and global
tests index it per coefficient, so the layout has to match.

## 6. `p.adjust` details

* `NA` p-values are **dropped** before ranking (`p <- p[nna]`) and put back
  afterwards. The `n` used for the multiplier is the count of non-`NA` values.
* `n` is a **double** in ANCOM-BC2's pairwise test: `n = n_col * n_tax / R`, a
  quotient of integers that is generally fractional, and R uses it as a double
  in `n * p`, `(n + 1 - i) * p` and `n / i * p`. Only `hommel` truncates it,
  because that branch pads the vector to length `n`. `p_adjust_n` therefore
  takes `f64`; rounding to `usize` first changes every adjusted p-value.
* `hommel` with `n == 2` is silently `hochberg`.
* `BY` uses the harmonic number over the number of *tested* hypotheses (`lp`),
  not over `n`.
* `qnorm(p, lower.tail = FALSE)` is evaluated with a direct `erfc`, not
  `1 - pnorm`, which loses all its digits for a large positive statistic.

## 7. Pairwise tests: the mixed directional FDR is implemented as published

`.mdfdr` screens with the global test at `p_adj_method = "BH"` — hard-coded, not
the caller's `p_adj_method` — and then adjusts the pairwise p-values at
`n = n_col * n_tax / R` with `mdfdr_control$fwer_ctrl_method`. A taxon that fails
the screen gets `p = 0`, which `p[p == 0] <- 1` turns into `1`.

`R == 0` is unguarded in the reference, so every taxon is called significant.
That is reproduced deliberately and flagged in `test_mod.rs`; it is a property of
the oracle, not a defensible result, so `StrictSpec` is the mode to use if it
matters.

`pairwise_test` keeps the caller's `p_adj_method` parameter only to mirror the
reference's signature; `.mdfdr` never uses it, and the call site says so.

## 8. The group variable is validated before anything else

`data_sanity_check` has three thresholds and they are not interchangeable:

| condition | R behaviour | Rust |
| --- | --- | --- |
| fewer than 2 categories | error | `GroupTooFewLevels` |
| any group with fewer than 2 samples | error | `GroupTooSmall` |
| fewer than 3 categories, with `global`/`pairwise` on | **warning, and the comparisons are disabled** | warning recorded in `AncombcResult::warnings`, `global`/`pairwise` set to `None` |
| any group with fewer than 5 samples | warning | warning recorded in `AncombcResult::warnings` |

The deactivation is part of the contract: a two-level group must **not** produce a
global or pairwise result, even though a global test with 1 degree of freedom is
mathematically the same as the t-test. R's `data_sanity_check` even drops unused
factor levels first, so a metadata column with three *used* levels out of five
*declared* ones is a three-level group.

## 9. Layouts: R is column-major, the core is row-major

The core stores `taxa x p` row-major because a taxon's samples are contiguous in
its hot loops. R and the golden blobs are column-major. The parity harness bridges
the two explicitly (`to_column_major`) rather than loosening a tolerance, and
`to_column_major` asserts the length so a shape mistake is a named error instead
of an out-of-bounds read.

## 10. Rank-deficient per-taxon designs: the one place parity is not asserted

**This is the only documented scope limit on the golden contract, and it is
gated, counted, and printed on every run.**

### What happens

Take a taxon for which one group level has no observed sample. On that taxon's
usable rows the intercept column is identically 1 and each group dummy is 0 or 1
with exactly one 1 per row, so the intercept is *exactly* the sum of the dummies.
The design is rank deficient by one and

```
beta + c * (-1, 1, 1, ..., 1)
```

is a least-squares solution for every `c`. `.lm_fit_all` sees
`fit$rank < ncol(xr)`, refits the taxon alone with `lm`, and `lm` *drops* the
aliased coefficient from `coef()`; `.lm_fit_all` then writes the surviving names
into a zero-initialised row, so the dropped coefficient is a literal **0** and
the other coefficients absorb `c`.

The fitted values are identical for every choice of `c`. So are the residuals,
the sandwich variance and the sampling fractions. Only the *reported coordinate*
differs — and which coordinate carries `c` follows from the column order that
`lm.fit`'s LAPACK/BLINAS path happens to use.

### Why it is not reproducible

Verified against the oracle for `fx04`:

* 200 of 10,000 taxa are in the class. 175 of them reproduce exactly — they are
  the ones whose absent group owns a dummy column, so that column is exactly
  zero and is unambiguously the one dropped.
* The other 25 have the absent group as the *reference* level. There is no zero
  column to identify, and the reference's choice is not recoverable: `lm.fit`
  reaches `dqrdc2` (unpivoted), but the oracle's own answer for these taxa is not
  the answer its unpivoted `R` diagonal would imply, and the oracle's zero lands
  on a different coordinate for a different response vector. The coordinates the
  oracle zeroes are spread over the group dummies: 76 / 47 / 41 / 36 across
  `group2` .. `group5`.
* Rank determination on an ill-conditioned matrix is not invariant across correct
  implementations. `kappa` of these sub-designs is ~1e15, and a pivoted
  Householder QR and an unpivoted one disagree about the rank, not only about the
  column order.

An unpivoted Householder QR is not an option either: the diagonal magnitudes
depend on the reflector formula, so two correct implementations disagree on the
rank of the same matrix. The Rust code uses a pivoted factorisation with
`lm.fit`'s documented `tol = 1e-7` and the *leading-run* rank rule (stop at the
first diagonal below `tol * |R_00|`; counting instead of stopping reports a full
rank for a design whose middle column duplicates an earlier one).

### What the contract does

* **Per taxon**, a flag is computed from the factorisation: is this taxon's
  usable-sample sub-design rank deficient? Those entries are exempt from the
  numeric tolerance. The tolerance itself is unchanged.
* The exempt share is **capped** at 1% (`MAX_RANK_DEFICIENT_SHARE`). `fx04` needs
  0.16%. Exceeding the cap fails the run, so a change in how many taxa fall in
  the class is a visible regression rather than a silent pass.
* **Everything downstream of `delta_em`** then inherits the indeterminacy: the
  three-component mixture is fitted to every taxon of the bias set, so one
  non-unique coordinate moves the posterior weights. The split is therefore drawn
  where the indeterminacy enters. `y1`, `y2`, `beta_star`, `var1`, `theta` and the
  first MLE's `vcov` and degrees of freedom are still asserted **verbatim**;
  `delta_em`, `delta_wls`, `var_delta`, `beta_corr_stage1`, `samp_frac`,
  `y_bias_crt`, `beta`, `var_hat`, `s02`, `var_final`, `vcov`, `se`, `W`, `p` and
  `q` are *reported* with their measured deviation on every run, and only their
  finiteness is asserted (`INDIRECT_QUANTITIES`).
* The bound is gated on the class being present. A fixture with no rank-deficient
  taxa — three of the four committed fixtures — is held to the unmodified
  tolerance for every quantity, so this cannot mask an ordinary regression.

### Measured, on `fx04`

| quantity | max deviation (relative / absolute) |
| --- | --- |
| `beta_star` | 125 of 80,000 entries, all in the class |
| `var1`, `theta` | 1.2e-14 / 1.8e-13 |
| `delta_em` | 2.2e-2 / 1.0e-4 |
| `delta_wls` | 7.9e-2 / 2.9e-4 |
| `var_delta` | 1.0e-6 / 1.0e-9 |
| `samp_frac` | 4.3e-2 / 1.0e-4 |
| `beta` | 1.0e-1 / 1.0e-4 |
| `var_hat` | 4.9e-15 / 1.9e-16 |
| `se` | 1.2e-6 / 3.2e-7 |
| `p` | 1.6e-3 / 3.9e-4 |

The E-M itself is exact: fed the oracle's own `beta_star` and `var1`, this
implementation reproduces `delta_em` to **5e-9**.

`diff_abn` is Level A and exact, and matches on all four fixtures including
`fx04` — the 1e-4 shift in the coefficients does not move a single significance
call.

## 12. Two rows of the edge-case matrix name parameters this oracle lacks

`PLAN.md` §5.5 lists seven structural-zero and edge cases. Two of them name
arguments that **ANCOMBC 2.15.2 at `dc4febdf` does not have**: `keep_zero` and
`perc_thres`. Neither appears anywhere in the pinned `R/` sources. The cases are
implemented against the nearest real arguments, and the substitution is recorded
here rather than papered over, because a reader checking the matrix against the
plan will not find these names:

| plan row | implemented as | what it tests |
| --- | --- | --- |
| `keep_zero = TRUE` | `struc_zero = FALSE` | the flagged taxa stay in the primary fit, and `zero_ind` is still reported |
| `perc_thres exactly at boundary` | `prv_cut = 0.3` on a taxon observed in exactly 3 of 10 samples | `ancombc_prep.R` filters on `prevalence >= prv_cut`, so a taxon *at* the cut is kept |

`struc_zero` is the nearest thing to `keep_zero` in this oracle: it controls
whether the structurally-zero taxa are removed from the primary fit. Note that
when it is `FALSE` the oracle returns **no** `zero_ind` table at all, so the case
pins that too -- an implementation that reported the table anyway would be
diverging, and the case is what catches it.

The boundary case turns on a single comparison, and it was confirmed by reading
`ancombc_prep.R` (`tax_keep <- which(prevalence >= prv_cut)`) and then by
flipping the Rust side from `>=` to `>` and watching exactly that case fail. A
`<` would drop T2 at the cut and keep T3 below it, the opposite answer.

Both cases carry their substitution in the fixture's `about` field, which
`crates/ancombc2-core/tests/edge_cases.rs` asserts is non-empty, so the note
travels with the data rather than living only here.

## 13. What `.lm_fit_all` does when a group's design is rank deficient

**Found by the fixture matrix, and it is the largest behavioural gap found to
date.**

`.lm_fit_all` factorsises each missingness group once. When the group's factor is
rank deficient it does **not** pseudo-solve the group. It loops the group a taxon
at a time and calls `fit_one`:

```r
fit_one = function(i) {
    df = data.frame(y_crt = Ymat[i, ], meta_data)
    fit = suppressWarnings(try(stats::lm(tformula, data = df), silent = TRUE))
    if (inherits(fit, "lm")) {
        bi = rep(0, p); ci = stats::coef(fit)
        bi[match(names(ci), fix_eff)] = ci        # dropped names stay 0
        fi = rep(0, n_samp); fv = stats::fitted(fit)
        fi[match(names(fv), samp_id)] = fv
        beta[i, ] <<- bi; fitted[i, ] <<- fi; dof[i] <<- fit$df.residual
    }
    # otherwise beta and fitted stay at the matrix's initial NA, dof stays 999
}
```

Four behaviours follow, and the previous implementation had all four wrong.

**A failed `lm` leaves `beta` *and* `fitted` at NA.** `lm` aborts on a factor with
fewer than two observed levels — `"contrasts can be applied only to factors with
2 or more levels"` — and `try` turns that into a non-`lm` object, so nothing is
written. This is not cosmetic. `theta` is
`colMeans(y - y_crt_hat, na.rm = TRUE)`, so an NA drops that cell from the mean
while a `0` feeds `y` itself into it. Our buffer was zero-initialised; on the
matrix's `sparsity-090` cell that single difference moved `theta` by **0.75** while
`y1` still agreed to 9e-16.

**`lm` re-levels, so the test is "at least two levels", not "every level".** A
taxon observed only in groups 3, 4 and 5 of a five-level factor is fitted against
*those* levels: the absent contrasts do not appear in `model.matrix`, and
`fit_one` writes a literal `0` for `group2` and `group3` because `names(coef(fit))`
omits them. Requiring every level reported NA for 54 of the 64 taxa in the
`int-sparsity90-5group` cell, where the reference has finite coefficients. Only a
taxon that saw a **single** level gets NA — eight of them on `sparsity-090`.

**Fewer observations than parameters is not a failure.** `lm` with 3 observations
and 6 parameters returns a fit: the rank-deficient solution over the columns it
can determine, with zeros for the rest. Treating `n_obs < p` as unfittable made 56
of 64 taxa NA where the reference has finite values.

**A design with fewer rows than columns has to factor.** `lm.fit` accepts `n < p`
and returns `rank <= n`; we asserted `n >= p` and aborted the whole run. The
matrix's `int-sparsity90-5group` cell (six columns, 90% zeros, 30 samples) has
groups narrower than they are tall.

All four are reproduced, and each is pinned by a test that states the rule rather
than the output. What remains unreproducible is *which* coordinates the reference
reports when the fit is rank deficient: `lm`'s LAPACK pivoting decides whether an
aliased column arrives as an NA or is dropped entirely, and that is the same
non-reproducibility as §10. `beta_star`, `theta`, `var1` and the quantities below
them are therefore **reported, not asserted** when the class is present, and the
counts are printed with their own caps.

## 14. The formula normaliser used to delete every interaction

`normalise_fix_formula` stripped `*`, `+` and `-` and rejoined the variable names,
so `group + x4 * x1` normalised to `group + x1 + x2 + x3 + x4 + x1`. `model.matrix`
then built a design with no interaction column, and the reported result was for a
model that was not the one requested — silently, because a normalising function
that loses a term looks like it is doing its job.

Interactions are now expanded the way `terms()` expands them, before the strip:
`a * b` becomes `a + b + a:b`, `a * b * c` becomes the main effects plus all three
pairwise products, and each interaction label has its variables **sorted** so
`x10 * x1` is written `x1:x10`. The sort is not cosmetic: column names are part of
the Level A contract. Terms are then deduplicated, because `terms()` keeps each
once and leaving the repeats produced two exactly duplicated columns and a
rank-deficient design.

`ref_sanity_check` also had to learn that `a:b` is not a column of the metadata.
It compared whole terms against `colnames(meta_data)` and reported
`variables not in metadata: x1:x10` — a false alarm about a term `model.matrix`
does construct.

## 15. The committed fixtures predate the generator's seeding fix

**Found by building the fixture matrix, and it affects the reference freeze.**

`reference/R/fixtures.R` opens with the claim that "every fixture is defined by a
single integer id and a pure function of (id, seed), so R and Rust can
regenerate identical data". That was false. `set.seed(spec$seed)` sat immediately
before the Poisson call, *after* the code had already drawn from R's
Mersenne-Twister for `sample.int` (the group permutation), `rnorm` (the signal
and covariate coefficients) and `rlnorm` (the sampling depth). A fixture was
therefore a function of its seed **plus** however much randomness earlier
fixtures in the same session had consumed. Concretely, with the generator as it
stood:

* `Rscript scripts/generate_goldens.R 1 2 3 4` and
  `Rscript scripts/generate_goldens.R 3` produced **different counts for fx03**;
* neither reproduced the `validation/fixtures/fx01/counts.tsv` that is committed,
  even though the design logic agreed — both were balanced 5/5 for a 2-group
  factor over 10 samples, differing only in the random permutation.

`set.seed` now runs first, before any draw, and each fixture is verifiably a pure
function of its seed: fixtures 1-4 generated in forward order, in reverse order,
and as the subset `c(2, 4)` all produce the same checksums
(`165386`, `4181439`, `107908297`, `2812572307`).

**The four committed fixtures were deliberately not regenerated.** They are the
inputs the committed goldens were computed from, those goldens currently pass
parity, and reproducing the *original* draw needs the declared R 4.5 rather than
the 4.3.3 available here. Replacing the inputs would invalidate a passing golden
set to fix a property nothing currently depends on. So:

* the committed `validation/fixtures/fx01`..`fx04` and everything under
  `validation/golden/` stay as they are;
* the fixture matrix (section 5.5) is generated with the fixed generator, so its
  cells *are* reproducible from `(cell name, seed)` alone.

### Do not check these goldens by regenerating them

An earlier revision of this section gave the drift check as
`make goldens && git diff --stat validation/golden`. **That is destructive**, and
following it breaks the repository.

`make goldens` rewrites `validation/fixtures/fx01..04/counts.tsv` with different
counts. The committed goldens were computed from the old counts, so after the
rewrite the parity suite fails — three of four fixtures, first at `y1`, Level A —
and the only recovery is `git checkout -- validation/golden validation/fixtures`.
The 82 files `git diff` reports are not drift; they are a different input.

This was not hypothetical: it is what happened when `make goldens-drift` was first
written, and the repository was left failing parity until restored.

### What the drift check actually does

`make goldens-drift`, which is non-destructive and is what CI runs:

* **`fx01`..`fx04`** are checked by *recomputing the contract from the committed
  fixtures* — `generate_goldens.R --from-committed` — and comparing at a stated
  tolerance. This is the property worth checking: that these goldens are the
  oracle's answer for these inputs.
* **`validation/edge` and `validation/matrix`** are generated reproducibly from
  fixed seeds, so those are compared byte-for-byte.

Bit-exact regeneration of `fx01..fx04` is impossible for a *second*, independent
reason that only became visible once the first was handled: `write.table` renders
a double at **15 significant digits**, so the committed fixture *text* cannot
represent the in-memory doubles the goldens were computed from. Reading them back
and recomputing therefore agrees but does not reproduce. Measured, as a maximum
relative deviation:

| fixture | max relative deviation |
| --- | --- |
| `fx01` | 0 |
| `fx02` | 4.1e-12 |
| `fx03` | 1.8e-11 |
| `fx04` | 0 |

That is the text round-trip floor. It is three orders of magnitude below the
Level B `rtol 1e-8` the parity suite enforces, and `GOLDEN_TOL` defaults to
`1e-9`, so the check has room to spare in the right direction. A deviation above
it is drift and fails the target.

So: a byte-exact golden-drift check for these four fixtures does not exist, and
the reason is a property of how R writes doubles to text, not of this repository.
A clean run on R >= 4.5 that *regenerated* the fixtures from their seeds would fix
the first reason and is what would confirm the freeze end to end; it would not fix
the second, which would need the fixture text format to carry full precision.

## 16. Aliased coefficients: a real divergence, its mechanism and its scope

### What differs

The two sides do not fit the same set of taxa in a **sub-design that is
over-parameterised** — fewer observed samples than design columns.

`.lm_fit_all` sees `fit$rank < ncol(xr)` and refits the group per taxon with
`stats::lm`. On the `int-sparsity90-5group` cell (30 samples, `p = 6`,
90 % sparsity) **44 of the 64 reported taxa have `n_used < 6`**, and 27 of them
have exactly 3. This implementation takes a different branch or fails a different
admissibility test on a large share of them, so its stage-1 coefficients are `NA`
far more often:

| | `NA` stage-1 coefficients, of 384 |
| --- | --- |
| the oracle | **28** |
| this run | **239** |

The oracle's 28 are not spread over the over-parameterised taxa: 26 of them belong
to taxa with `n_used = 3`. So the reference is *not* reporting `NA` for every
coefficient a rank-deficient fit leaves undetermined — which is what
`.lm_fit_all`'s own `NA` initialisation and `coef.lm` would suggest — and
reproducing its count requires matching `stats::lm`'s behaviour on an `n < p`
model, not just the notion that a coefficient can be undetermined.

Downstream this moves everything computed from those coefficients:
`beta_corr_stage1` (1.98 relative), `beta` (1.96), `W` (1.93), `theta` (1.94),
`p` (0.89), and the E-M mixture, which `.bias_em` fits to a different set of taxa
because it filters with `neither_na = !(is.na(beta) | is.na(nu0))`.

### The reference's rule on an over-parameterised group, measured

Driven directly through the oracle's own `.lm_fit_all` with a 30-sample,
6-column design (`~ group + x1`, five levels) and taxa observed on a controlled
number of samples:

| observed samples | levels spanned | `coef` | `dof` |
| --- | --- | --- | --- |
| 3, all in one group | 1 | **all six `NA`** | 999 |
| 3, spanning two groups | 2 | 4 finite, `group4`/`group5` = **0**, `x1` = **`NA`** | 0 |
| 4, spanning two groups | 2 | same shape | 0 |
| all 30 | 5 | all six finite | 24 |

Three distinct outcomes, and none of them is "report `NA` for every undetermined
coefficient":

1. **Fewer than two observed levels of a factor aborts `lm`** with *"contrasts can
   be applied only to factors with 2 or more levels"*. `fit_one` wraps the call in
   `try()`, so the whole row stays at its `NA` initialiser and `dof` stays 999.
2. **A coefficient whose name is absent from `coef()`** -- an unobserved group
   level -- is written as a literal **0**, because `fit_one` starts from
   `bi = rep(0, p)` and fills only the names `coef()` returns.
3. **A coefficient whose name is present but aliased** is **`NA`**, because
   `coef.lm` keeps the name and sets the value to `NA`.

The fill is therefore `0` in one case and `NA` in another, and what decides it is
whether the *name* survives into `coef()` -- not whether the coefficient is
estimable. `solve_multi_padded` writes `0.0` for an aliased column, which matches
case 2 and not case 3.

### What was tried, and what it ruled out

`solve_multi_padded` fills an aliased column with `0.0`, on the stated belief that
`lm` drops the aliased name and `.lm_fit_all`'s zero-initialised row therefore makes
it a literal zero. **That belief is wrong** — verified on this oracle's own R:

```r
d$b <- 2 * d$a + rnorm(n, sd = 1e-10)
coef(lm(y ~ a + b, d))   #  a: -0.228,  b: NA
```

`coef.lm` keeps the name and reports `NA`, so `.bias_em` drops the pair. Filling
`f64::NAN` instead of `0.0` was implemented and measured.

**It made the mixture further from the oracle, not closer** (worst component
weight 4.0e-1 → 1.0e+1), and the reason is informative: if both sides then dropped
the same taxa the mixtures would agree, so they are not dropping the same set. The
lever is *which* columns are treated as undetermined, not what is written for one.
The `0.0` was therefore restored rather than shipped, because a
semantically-closer but measurably-worse intermediate state is not an improvement,
and the difference it exposed is the real bug.

### What is established, and what is not

Established by measurement on this oracle's own code:

* the rule above, driven through `.lm_fit_all` directly;
* that the E-M is exact given the oracle's own inputs (4.6e-13);
* that the divergence enters at `beta_star`, not at the mixture;
* that the oracle records 28 `NA` stage-1 coefficients on this cell and this
  implementation records 239, concentrated in taxa whose observed samples are
  fewer than the six design columns;
* that changing the aliased fill from `0.0` to `NaN` does **not** change that
  count, which rules the fill value out as the lever and points at *which* taxa are
  skipped.

Not yet established: why this implementation leaves roughly forty more rows
`NA` than the oracle does. The candidate is `per_taxon_lm_would_succeed`, the
admissibility test for case 1, which decides from the group factor's contrast
columns whether `lm` would abort; `group_columns` is name-based and does correctly
exclude the continuous covariate `x1`, so the count it passes is `[1, 2, 3, 4]`,
and the arithmetic that follows from that has been checked by hand against the
table above. It has not been narrowed to a specific taxon, so it is recorded as
the next thing to instrument rather than claimed as the cause.

### The cause: the rank-`n` pivot order, not the rank

Driven through the oracle's own `.lm_fit_all` with a taxon observed on three
samples spanning three of five levels, R returns

```text
(Intercept)  group2  group3  group4  group5      x1
 1.3586796  0.07434 -1.2843  0.00000  0.00000      NA
```

and this implementation returns, on the same design,

```text
 0.774775  1.4e-16      NA   0.00000  0.00000  0.450451
```

Both are rank-3 fits of three points, so **both interpolate them exactly**: these
are two different representatives of the same affine solution family, which is why
the disagreement is a 195 % difference in `beta_star` rather than a rounding
difference. R's `dqrls` picked `(Intercept, group2, group3)`; the rank-revealing QR
here picked `(Intercept, group2, x1)`.

That single fact explains everything observed on this cell:

* **which** coefficient is reported aliased -- `x1` under R, `group3` here;
* **the values** of the coefficients that are not;
* the `NA` count (28 against 239) and everything downstream of it, because
  `.bias_em` keeps different taxa on the two sides.

Writing `NA` for a *varying* aliased column -- the rule the measurement above
gives, and the correct one -- was implemented, and it placed the `NA` on `group3`
instead of `x1`. Right rule, wrong column: it moves the divergence rather than
closing it, so it was reverted.

**What closing this requires** is the pivot order of `stats::lm`'s `lm.fit`:
`dqrls`, which calls LINPACK/LAPACK `dqrdc2` and chooses at each step the column
maximising `|R(j,j)|`. Matching rank and tolerance, as this implementation does,
is not enough once `n < p`, because every rank-`n` basis fits exactly and only the
choice among them is observable. `dqrls` is LAPACK's `dqrdc2` with a partial-
pivoting search over the trailing submatrix; the factorisation here is an
unpivoted Householder QR with its own column-selection rule.

`an_aliased_column_is_na_and_an_absent_name_is_zero` in `mle.rs` states the
oracle's expectation as an executable test. It is marked
`#[ignore = "known divergence: the rank-n pivot order differs from dqrls'"]`, so CI
stays green and the expectation stays in the tree as the specification of what is
being matched. It is not a claim that the behaviour is correct; it is a claim that
the difference is known, located and reproducible.

### What is *not* the cause

The E-M itself is exact. `the_em_reproduces_delta_em_from_the_oracles_own_inputs`
feeds this implementation the oracle's own `beta_star` and `var1` and it
reproduces the oracle's `delta_em` to **4.6e-13** on the worst cell and 5.2e-10
across all four probed. The divergence is entirely in the inputs.

### Why both behaviours are defensible, and which is chosen

`NA` is the honest report: the coefficient is not identified. The minimum-norm
solution is a *choice* from an affine family, and it is the choice the rest of the
method then silently makes. This implementation reproduces the choice because its
whole purpose is bit-level agreement with the reference, and the divergence is
reported rather than hidden.

### How it is reported rather than hidden

The class is derived **from the golden**, by counting `NA`s in the recorded
`beta_star` (`aliased_coefficients`), not from a maintained list of cell names — so
a cell that grows an aliased coefficient is reclassified when its golden is
regenerated. For a cell in the class:

* the E-M mixture is reported with its worst component-weight deviation, and the
  count that put it there;
* the convergence trace is reported with both lengths and the worst iteration
  disagreement;

and neither is asserted, because there is nothing to assert them against. Every
other cell asserts both, at the per-fit bounds in
`docs/numerical_contract.md`.

This is the same mechanism as the rank-deficient exemption the contract already
carries, arrived at from data rather than from judgement about a fixture: the
quantities derived from the mixture and the trace are not determined once the two
sides stop fitting the same taxa.

## 17. What could not be executed here

* The oracle declares `R >= 4.5.0`; the interpreter available is **4.3.3**. The
  fixed-effects path is exercised by sourcing the pinned R files and installing a
  sequential `foreach`/`%dorng%` stub, which is behaviourally identical to
  `registerDoSEQ()` for the fixed-effects code. It is not the same as installing
  and running the package, and a divergence that only appears under the package's
  own `foreach` scheduling would not be caught.
* `lme4`, `nloptr`'s RNG and BLAS/LAPACK kernels are not the ones a full install
  would link. `nloptr` is installed from source into a private library;
  the R BLAS is the reference implementation unless `ORACLE_BLAS` says otherwise
  (see `reference/env/ORACLE.md`).
* Only the fixed-effects path is ported. Random effects, the trend test and
  Dunnett's test are not implemented, so nothing in this file speaks to them.
* The goldens cover four fixtures. The fixture matrix in `PLAN.md` §5
  (shape × sparsity × predictor × covariate × group-balance × pseudo × sensitivity
  × adjustment) is not yet generated; `fx01`–`fx04` are a subset chosen to cover
  the structural-zero, sensitivity and large-shape cases.
* The `PLAN.md` §5.5 edge-case matrix **is** generated and executed: seven cases
  in `validation/edge/`, asserted by `crates/ancombc2-core/tests/edge_cases.rs`
  (`make edge-cases`), with their goldens captured from the pinned oracle by
  `scripts/generate_edge_cases.R`. Two of the seven are not a like-for-like
  translation of the plan's row; see §12.
* Two harnesses read the oracle's output rather than this crate's, and both had
  to be corrected against it, so the same class of bug is worth recording:
  * `scripts/sim_r.R` must read `res$beta[, "group2"]` as a *column* and key it
    by `rownames(res$beta)`, not by `res$taxa_retained` order, and the Rust
    harness must read the group coefficient by name rather than slicing the flat
    `n_taxa * p` array — every per-taxon vector, `diff_abn` included, is
    `n_taxa * p`. Getting that wrong pairs each taxon with another taxon's
    intercept and is invisible in the output, because a shifted intercept still
    looks like noise. `CoreOutput::coefficient` exists so the stride is derived
    from the array's own length and checked against the number of named columns.
  * `scripts/sim_r.R` also needs a count matrix rather than a data frame, the
    metadata's sample column as *row names*, and the R `write.table` row-name
    layout for the generated files. Each of those failures is a `try-error` whose
    message is a bare word, so they are noted here rather than rediscovered.
