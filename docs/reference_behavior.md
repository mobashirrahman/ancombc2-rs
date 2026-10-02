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

## 10. Rank-deficient per-taxon designs: parity is asserted here too

### What happens

Take a taxon for which one group level has no observed sample. On that taxon's
usable rows the intercept column is identically 1 and each group dummy is 0 or 1
with exactly one 1 per row, so the intercept is *exactly* the sum of the dummies.
The design is rank deficient by one and

```
beta + c * (-1, 1, 1, ..., 1)
```

is a least-squares solution for every `c`. `.lm_fit_all` sees
`fit$rank < ncol(xr)`, refits the taxon alone with `lm`, and `lm` re-levels the
factor to the levels that taxon actually observed, so the design it factorises is
not the design the group was fitted on.

### Why it is reproducible

An earlier revision of this section said it was not, and gave the reason that
`lm.fit`'s LAPACK/BLINAS path decides which coordinate an aliased coefficient
lands on. That reason was wrong: `lm.fit` and `lm` call `dqrls` with
`pivot = FALSE`, so the permutation is the identity — `lm.fit(X, y)$qr$pivot` is
`1 2 ... p` even when an early column is an exact combination of later ones — and
the dropped column is the last aliased one in build order. There is no pivoting to
reproduce, and the representative is fully determined.

The rule that had to be implemented instead is three-part, and all three parts are
measured against the oracle's own `.lm_fit_all`:

1. **Factorise `lm`'s design, not the group's.** `fit_one` calls `lm` on a frame
   of *all* samples with `NA` where the taxon's response is missing; `lm` drops the
   incomplete rows, then `model.matrix` applies `drop.unused.levels = TRUE`. Every
   group contrast for a level the taxon never observed is therefore absent. An
   unobserved level's contrast has no name in `coef()`, and `fit_one` writes into
   `rep(0, p)`, so it is a literal **0**; a contrast that *is* named but aliased is
   **`NA`**.
2. **The base level is the first observed one**, not the first globally. A taxon
   observed at levels {3, 4, 5} gets columns for 4 and 5 only, so `group3` is `0`
   even though `group3` is 1 on some of its rows.
3. **`stats::fitted` is the projection**, not `X %*% coef`. It stays finite where
   `coef()` is `NA`, which matters because `theta` is
   `colMeans(y - fitted, na.rm = TRUE)` and an `NA` there drops the taxon out of
   that sample's mean entirely.

Rank determination itself uses `lm.fit`'s documented `tol = 1e-7` and the
*leading-run* rule — stop at the first diagonal below `tol * |R_00|`, rather than
counting every diagonal that clears it, since counting reports a full rank for a
design whose middle column duplicates an earlier one.

### What the contract does

**Nothing special.** Every quantity in the contract is compared at its own
tolerance, and the tolerance is never widened for these taxa. The rank-deficient
class is still computed and still printed on every run — how many taxa are in it,
how many of those have a per-taxon `lm` that fails outright, how many are
under-determined — because the counts are worth seeing and a change in them should
be visible. The share caps fail the run if it moves; they gate nothing about what
is compared.

The twenty-entry report-only list that used to cover this class is gone rather than
left empty; §16 records why, and what would have to be supplied to bring a
genuine exemption back.

### Measured, on `fx04`

200 of its 10,000 taxa are in the class, and the largest deviation anywhere in the
contract is now at the tolerance the quantity itself carries. `diff_abn` is Level A
and exact, and matches on all four fixtures.

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
than the output. The last of them — "an `NA` aliased column becomes `NA` in the
reported coefficients" — was the open question for a while, on the belief that
`lm`'s pivoting made it unreproducible; it is `dqrls(pivot = FALSE)` and it is
reproducible. See §10 and §16. Every quantity is asserted, and `beta_star`,
`theta` and `var1` agree at their own tolerances on every cell including this
one.

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

## 16. Aliased coefficients: the divergence, and what actually caused it

**Status: closed.** Every quantity in the contract is asserted on every committed
fixture and every cell of the fixture matrix, including the ones whose groups are
rank deficient. Nothing in the contract is reported instead of asserted. This
section is kept because the cause was not what the previous revision of it said,
and the mistake is the kind that survives a re-read.

### The symptom

`.lm_fit_all` sees `fit$rank < ncol(xr)` for a group and refits each of its taxa
one at a time with `stats::lm`. On the `int-sparsity90-5group` cell (30 samples,
`p = 6`, 90 % sparsity, five group levels) that is 56 of 64 taxa. They came out
badly wrong: `beta_star` was off by up to **195 %** on a coefficient, and
`theta` — which is `colMeans(y - fitted, na.rm = TRUE)` over the whole set — was
off by up to 1.4 absolute.

### What the previous revision of this section claimed

That the reference's choice of least-squares representative in that path "depends
on `lm`'s LAPACK pivoting, which is the same non-reproducibility [...] here it is
not a family of equivalent least-squares solutions but an arbitrary choice of which
coordinates to report at all. [...] there is no rule to reproduce."

That claim was false, and it was falsifiable in about a minute on the oracle's own
R:

```
> lm.fit(cbind(c1, c2, c1 + c2), y)$qr$pivot
[1] 1 2 3
```

`lm.fit` and `lm` call `dqrls` with `pivot = FALSE`, so the permutation is always
the identity — there is no pivoting to reproduce at all — and the column reported
`NA` is the *last* aliased one in build order. The reference's representative is
fully determined. Everything downstream of the false claim (a twenty-entry
"report-only" list, a share cap, and the claim that the affected quantities could
not be compared at any tolerance) came from it.

### The three real causes, in the order they were found

Each was found by shrinking the disagreement, not by reasoning, and each was
verified against the oracle before the next was looked for.

**1. The fitted values were `X %*% coef`, but `stats::fitted` is the projection.**

`coef(fit)` carries `NA` for every aliased term, so multiplying it through the
design made the whole fitted row `NA`. Since `theta` is
`colMeans(y - fitted, na.rm = TRUE)`, an `NA` does not merely lose one term: it
drops the taxon out of that sample's mean entirely. Measured on
`int-sparsity90-5group`, only about a fifth of each sample's observed taxa were
contributing, three samples had no contributor at all, and `theta` came back `NaN`
— which then made the whole adjusted-response column `NaN` on the next iteration
and cascaded. Zeroing the aliased coefficients is not an approximation of the
projection, it *is* the projection: an aliased column is by definition in the span
of the columns that were kept.

**2. The rank-deficient path factorised the group's design, not `lm`'s.**

`fit_one` builds `df = data.frame(y_crt = Ymat[i, ], meta_data)` over **all**
samples and calls `lm`, which drops the incomplete rows and then has
`model.matrix` apply `drop.unused.levels = TRUE`. So the design `lm` factorises is
not the group's: every group contrast for a level the taxon never observed is
absent. That is a different matrix, and it has to be factorised as one — an
interior rejection of columns cannot be reproduced by factorising the full design
and patching the answer, because the two pick different bases.

Two details, both measured against the oracle's own `beta_star`:

* the contrast for an **unobserved** level has no name in `coef()` and is written
  as a literal `0`, from `fit_one`'s `bi = rep(0, p)` initialiser;
* the **base** level is the first *observed* one, not the first globally. A taxon
  observed at levels {3, 4, 5} gets columns for 4 and 5 only, so `group3` is `0`
  even though `group3` is 1 on some of its rows. Reading the base off the contrast
  indices instead flipped 20 of the 64 rows.

**3. The design's column names never reached the core.**

`build_design` kept them in `Design::colnames` and left `Matrix::colnames` empty,
so the pipeline's `group_columns(&x.colnames, ...)` returned an empty slice. With
no group columns, cause 2's filter is a no-op and every rank-deficient taxon is
fitted as though it had observed every level. `fix_eff` was unaffected because it
comes from `cfg.fix_eff`, which is why the defect was invisible in the output and
only showed up as a numeric divergence.

### The evidence

With all three fixed, an independent reimplementation of `.lm_fit_all` driven
directly from the golden's own `x` and `y1` reproduces the oracle's recorded
`theta` to **0.0 relative error over all 30 samples**, the first iteration's
`epsilon` to all eight published digits (`0.99801595`), and `beta_star` row by row.
The Rust side now matches that reproduction.

### What was removed with the cause

The twenty-entry `INDIRECT_QUANTITIES` list — `theta`, `beta_star`, `delta_em`,
`beta`, `se`, `p`, `q`, `vcov` and the rest of everything downstream of the first
MLE — existed to carry those quantities as *reported, not asserted*. It is now an
empty `const`, and the report-only branches in `compare_indirect` and
`compare_per_term` are gone. The rank-deficient class is still computed and still
printed, because the counts are worth seeing, but it gates nothing.

An allowance that no test can trip is worse than no allowance: it reads as a known
gap in the contract. If a future fixture reintroduces a genuine
non-reproducibility, it should be added back with the R evidence for it, not
inherited from this paragraph.

### The two branches, stated plainly

1. **Fewer than two observed levels of a factor aborts `lm`** with *"contrasts can
   be applied only to factors with 2 or more levels"*. `fit_one` wraps the call in
   `try()`, so the whole row stays at its `NA` initialiser and `dof` stays 999.
2. **A coefficient whose name is absent from `coef()`** — an unobserved group
   level — is a literal **0**.
3. **A coefficient whose name is present but aliased** is **`NA`**.

What decides 2 versus 3 is whether the *name* survives into `coef()`, not whether
the coefficient is estimable.

### A related omission, now covered

`x` itself was in the golden and in the manifest but was never compared; only the
column *names* were. Adding it at Level A immediately failed on
`shape-10x10`, which is how cause 3 surfaced — the report could not otherwise name
the design as the thing that was wrong.

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
