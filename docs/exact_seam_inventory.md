# The numerical seam: field, shape, type and layout inventory

**Task:** IMPROVED_PLAN.md S06. **Scope:** the boundary between the retained
upstream R code and anything Rust will do. This file is the inventory; nothing
here is implemented yet.

Every entry was read out of the pinned source at
`reference/ANCOMBC/R/` and, where the type is not obvious from the code, measured
against the installed original on the `linux-r453-openblas` profile. Line
references are to that pinned commit.

## 0. What is inside the seam and what is not

`ancombc2()` is one function with four stages. Only the middle one is a
numerical core; the other three are argument handling, I/O classes and table
assembly, and they stay in R.

| stage | source | disposition | why |
| --- | --- | --- | --- |
| argument and class dispatch | `ancombc2.R:394-451`, `data_sanity_check.R` | **retained R** | formula semantics, `phyloseq`/`TSE` coercion, taxonomy, factor contrasts, and eleven distinct argument guards whose exact messages are part of the contract |
| preprocessing and fitting | `.data_core`, `.get_struc_zero`, `.ancombc2_core` | **the seam** | this is what Rust will replace, one stage at a time |
| sensitivity refits | `ancombc2.R:499-780`, `.ancombc2_sens_fit`, `.ancombc2_sens_p` | **the seam** | same core, different pseudo-count; the grid is the original's |
| result assembly | `ancombc2.R:783-805`, `ancombc_prep.R:250-395` | **retained R** | `cbind` order, column prefixes, `NULL` entries, `data.frame` types |
| random effects | `.iter_remle`, `.lmer_struct_all`, `.remle_fit_all` | **retained R** | `lme4`/`lmerTest`; no exact accelerated form exists yet |
| Dunnett and trend | `.ancombc_dunn`, `.ancombc_trend`, `.constrain_est`, `.mdfdr` | **retained R** | `quadprog::solve.QP`, `multcomp`, `B`-replicate null simulation |

The seam is deliberately narrow: `.ancombc2_core` is a pure function of
`(data, aggregate_data, meta_data, fix_formula, ...)` and returns a list of
arrays and tables. Everything the original does to build its arguments, and
everything it does with the answer, stays where it is.

## 1. `.ancombc2_core` inputs

`ancombc_prep.R:92`. `O1` is `data + pseudo` and `O2` is
`aggregate_data + pseudo` (`ancombc_prep.R:111-112`), so what reaches the seam is
the *pre-filtered, pseudo-added* pair.

| name | R type | shape | layout | dimnames | notes |
| --- | --- | --- | --- | --- | --- |
| `data` | integer or double matrix | `n_tax x n_samp` | **column-major** (R native) | rows = taxa, cols = samples | `NA` counts survive; `bench_r.R`'s `counts[is.na(counts)] <- 0` is *not* applied on the exact path |
| `aggregate_data` | same as `data` | `n_tax x n_samp` | column-major | identical dimnames | equals `data` when `tax_level` is NULL |
| `meta_data` | data.frame | `n_samp x k` | column-major per column | rows = samples | `group` is a factor; covariates may be numeric; **one column is fatal** (see §5) |
| `fix_formula` | character scalar | - | - | - | a bare RHS; `model.matrix(~ <fix_formula>, meta_data)` |
| `p_adj_method` | character scalar | - | - | - | one of the seven `p.adjust` methods |
| `pseudo` | numeric scalar | - | - | - | already added to `data`; passed for the `< 50 taxa` warning only |
| `s0_perc` | numeric scalar or NULL | - | - | - | NULL means `s02 = 0` and no `quantile` call |
| `group` | character scalar or NULL | - | - | - | a *column name*, matched with `grepl(group, covariates)` |
| `alpha` | numeric scalar | - | - | - | |
| `global`/`pairwise`/`dunnet`/`trend` | logical scalar | - | - | - | `trend` forces `global <- TRUE` internally (`ancombc_prep.R:481`) |
| `iter_control` | list | `tol`, `max_iter`, `verbose` | - | - | `max_iter = 0` runs no iterations and does not error |
| `em_control` | list | `tol`, `max_iter` | - | - | `max_iter = 0` **errors**: `delta_new` is only defined inside the loop |
| `lme_control` | `lmerControl` | - | - | - | only read on the random-effects path |
| `mdfdr_control` | list | `fwer_ctrl_method`, `B` | - | - | only read by pairwise/Dunnett |
| `trend_control` | list | `contrast`, `node`, `B` | - | - | `contrast` is a list of **square matrices** |

The design matrix is built on the R main thread, once
(`ancombc_prep.R:118-120`):

```r
options(na.action = "na.pass")          # keep NAs in rows of x
x = stats::model.matrix(formula(paste0("~", fix_formula)), data = meta_data)
options(na.action = "na.omit")           # switched back
fix_eff = colnames(x)
```

`options(na.action = ...)` is **global state the original mutates**. Anything
that replaces this call has to reproduce the same `model.matrix` result with
`na.action = "na.pass"`, including the `NA` rows it leaves in `x`.

## 2. The derived arrays, in the order the original makes them

| # | line | expression | shape | type | what can be non-finite |
| --- | --- | --- | --- | --- | --- |
| 1 | `prep.R:123` | `o1 = log(O1); o1[is.infinite(o1)] = NA` | `n_tax x n_samp` | double | `NA` for a zero count and for a negative count |
| 2 | `prep.R:125` | `y1 = o1 - rowMeans(o1, na.rm = TRUE)` | `n_tax x n_samp` | double | `NA` where every sample of the taxon is missing |
| 3 | `bias:177` | `.iter_mle(x, y1, theta = NULL)` -> `beta1`, `var_hat1` | `n_tax x p`, `n_tax x p` | double | `NA` for an unfitted taxon |
| 4 | `bias:614` | `.bias_em(beta1[,i], var_hat1[,i])` per coefficient | 3 scalars per `i` | double | stops if any variance is exactly 0 |
| 5 | `prep.R:180` | `beta1 = t(t(beta1) - delta_em)` | `n_tax x p` | double | |
| 6 | `prep.R:182-185` | `theta_hat[i,] = y1[i,] - rowSums(x * rep(beta1[i,], each = n_samp), na.rm = TRUE)`, then `colMeans(..., na.rm = TRUE)` | `n_samp` | double | `NA` per sample; the warning naming the samples is part of the contract |
| 7 | `prep.R:199-201` | `o2`, `y2`, `y_bias_crt = data.frame(t(t(y2) - theta_hat))` | `n_tax x n_samp` | double | `y_bias_crt` is a **data.frame**, and it becomes the returned `bias_correct_log_table` |
| 8 | `bias:196` | `.iter_mle(x, y2, theta = theta_hat)` -> `beta_hat`, `var_hat`, `dof`, `vcov_hat` | `n_tax x p`, `n_tax x p`, `n_tax x p`, list of `n_tax` `p x p` | double | `dof` is `999` for an unfitted taxon, numeric not integer |
| 9 | `prep.R:216` | `var_hat = sweep(var_hat,2,var_delta,"+") + 2*sqrt(sweep(var_hat,2,var_delta,"*"))` | `n_tax x p` | double | `+Inf` if `var_delta` is huge |
| 10 | `prep.R:222` | `s02 = apply(var_hat, 2, function(x) quantile(x, s0_perc, na.rm = TRUE))` | `p` | double | `NA` if a whole column is `NA` |
| 11 | `prep.R:233-236` | `W = beta_hat/se_hat`; `p = 2*pt(abs(W), dof, lower.tail = FALSE)`; `p[is.na(p)] = 1`; `q = apply(p, 2, function(x) p.adjust(x, method))` | `n_tax x p` | double | `p` is forced to `1` where `NA`, `q` is not |
| 12 | `prep.R:237` | `diff_abn = q <= alpha & !is.na(q)` | `n_tax x p` | logical | |

## 3. `.lm_fit_all` -- the least-squares primitive (S09)

`ancombc_bias_correct.R:11`. This is the stage where a substitute QR is most
likely to diverge, so it is pinned exactly.

* **Grouping.** `use = is.finite(Ymat) & matrix(x_ok, byrow = TRUE)` and
  `keys = do.call(paste0, asplit(use * 1L, 2L))`; taxa sharing a usable-sample
  pattern are solved together (`bias:44-47`). When every row of `x` is complete
  and `Ymat` is all finite, there is one group of all taxa.
* **The solve is `stats::lm.fit(xr, Yr)`** (`bias:56`) -- R's own LINPACK/LAPACK
  `dqrls`-based routine. **Not** a hand-written Householder QR. The plan's
  "avoid a general QR rewrite" is not a preference here: it is what the original
  calls.
* **Rank deficiency.** `if (fit$rank < ncol(xr))` each taxon in the group falls
  back to a per-taxon `stats::lm(tformula, data = df)` (`bias:22-30`), which
  yields `NA` coefficients for the aliased ones. There is **no** rank-deficiency
  exemption in the exact contract; this fallback *is* the original's answer.
* **Outputs.** `beta` (`NA` for an unfitted taxon), `fitted` (with
  `fitted[idx, !rows] = 0` for unusable samples, `bias:65`), and
  `dof = fit$df.residual` as a **numeric** vector, `999` where unfitted.

## 4. `.sandwich_vcov` -- the covariance primitive (S10)

`ancombc_bias_correct.R:95`.

* `XTX_inv = MASS::ginv(t(x_cc) %*% x_cc)` where `x_cc` drops incomplete design
  rows. **`MASS::ginv` is the primitive**; its tolerance is `tol = 1e-15` on the
  singular values and it is applied to `t(X) %*% X`, not to `X` directly.
* Row `j` of `XX` is `as.vector(x[j,] %*% t(x[j,]))` -- the outer product in
  **column-major** order of the `p x p` block.
* `term_j = outer(eps2[idx, j], XX[j, ]); term_j[is.na(term_j)] = 0.1`. The `0.1`
  substitution for a missing entry is a documented reference behaviour, not a
  bug to fix (`docs/reference_behavior.md` section 2).
* Accumulation is `sigma2_xxT = sigma2_xxT + term_j` over `j` in order, in
  blocks of `2^17 / (p*p)` taxa. **Block size changes nothing arithmetically**
  (each block's accumulation order is the same) but a rewrite must keep the
  per-`j` order.
* `v_i = XTX_inv %*% matrix(sigma2_xxT[k,], p, p) %*% XTX_inv`, `dimnames` from
  `fix_eff`.

## 5. Things a transport must preserve

| item | why |
| --- | --- |
| `NA_real_` vs `NaN` vs `±Inf` vs `-0.0` | `log(0)` is `NA` here, not `-Inf`: `ancombc_prep.R:124` assigns `NA` to infinite logs *before* anything else. A transport that turns `NaN` into `null` loses the distinction permanently. |
| integer vs double | `data` is integer when the counts are; `O1`/`o1`/`y1` are double. A count matrix widened to double changes nothing numerically but changes the returned `feature_table`'s type, which is part of the schema. |
| column-major layout | R's native order. `as.vector(x[j,] %*% t(x[j,]))` is column-major, and a row-major bridge silently transposes a symmetric-looking `p x p` block. |
| `dof = 999` as a double | `dof` is `rep(999L, ...)` initially but `fit$df.residual` assigns a double into it, so a fitted taxon's entry is `20` and an unfitted one is `999L`. The vector's *values* differ in type. |
| dimnames at every stage | `beta`, `var_hat`, `dof`, `theta` and every `vcov_hat[[i]]` carry them, and they end up in the returned tables. |
| `options(na.action)` | the original mutates a global; a transport that computes `x` with the default `na.omit` produces a *different design matrix* whenever metadata has missing values. |
| column order of `res` | `cbind(taxon, lfc_*, se_*, W_*, p_*, q_*, diff_*)` and then `passed_ss_*`, `diff_robust_*` appended by `flag_fun` in `diff_` suffix order. |
| `NULL` entries | `zero_ind`, `ss_tab`, `res_global`, `res_pair`, `res_dunn`, `res_trend` are `NULL` depending on the options. A replacement that returns an empty data.frame where the original returns `NULL` fails the schema check. |
| one-column `meta_data` | `[.data.frame`'s `drop` collapses it, and `ancombc2()` then fails. Not a transport concern, but it is why seven committed fixtures cannot be analysed through the public API, and a transport that reshapes metadata could accidentally hide it. |

## 6. What the bridge must not do

* **No JSON.** `plan` section 2 lists it as a required closure, and
  `r/ancombc2rs/R/ancombc2.R`'s `call_core()` is exactly the defect: `NaN`
  becomes `null` and the R reconstruction changes the object's structure.
* **No R API calls from worker threads.** `model.matrix`, `p.adjust`,
  `quantile`, `pt`/`pf`/`pchisq`, `nloptr` and `MASS::ginv` are called on the R
  main thread, or replaced by an attributed exact port.
* **No reassociated reductions.** `rowMeans`, `colMeans` and `apply` have
  specific orders; `docs/performance_plan.md` sections 7.1 and 7.3 record that
  reassociating them already broke the *tolerant* contract.
* **No lookup of the oracle.** The bridge reads the caller's arguments and
  returns numbers; it never reads `validation/exact/golden/` or
  `reference/`.
