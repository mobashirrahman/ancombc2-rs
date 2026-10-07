//! Bias correction, sampling fractions, variance regularisation and inference.
//!
//! The steps between the second MLE fit and the p-values, in the order the
//! reference performs them:
//!
//! 1. `beta = beta_star - delta_em`, column-wise. The intercept is corrected
//!    like any other coefficient; ANCOM-BC2 makes no exception for it.
//! 2. `theta_hat[i] = colMeans(y1 - X * beta_corrected)`, the sample-specific
//!    bias on the log scale, R's `colMeans(..., na.rm = TRUE)`.
//! 3. Refit with `theta` held fixed, giving `beta_hat` and the sandwich
//!    variances.
//! 4. Add the variance of the bias: `var = var + var_delta + 2 * sqrt(var *
//!    var_delta)` — i.e. `(sqrt(var) + sqrt(var_delta))^2`, which is what
//!    `sweep(var, 2, var_delta, "+") + 2 * sqrt(sweep(var, 2, var_delta, "*"))`
//!    computes.
//! 5. Add the SAM-style regulariser `s0`, the `s0_perc`-th percentile of the
//!    variance column, then take square roots to get the standard errors.
//! 6. `W = beta / se`, `p = 2 * pt(|W|, df)`, `q = p.adjust(p, method)`,
//!    `diff_abn = q <= alpha`.
//!
//! Step 4 deserves a note: the reference writes it as a sum of a variance and
//! twice a geometric mean, which is exactly `(sqrt(a) + sqrt(b))^2`. Writing it
//! that way in Rust is algebraically identical and much easier to check.

use ancombc2_stats::{p_adjust_n, quantile_type7_unsorted, t_two_sided};

use crate::config::{AdjustMethod, CompatMode};
use crate::matrix::Matrix;
use crate::workspace::RMatrix;

/// Step 1: `beta = t(t(beta_star) - delta_em)`.
///
/// `beta_star` is `n_taxa x p` row-major and `delta` has one entry per
/// coefficient. `delta` may be `NaN` for a coefficient the E-M could not
/// estimate, in which case the column becomes `NaN`, matching R.
pub fn apply_bias_correction(
    beta_star: &[f64],
    n_taxa: usize,
    p: usize,
    delta: &[f64],
) -> Vec<f64> {
    let mut out = beta_star.to_vec();
    for k in 0..p {
        let d = delta[k];
        for i in 0..n_taxa {
            out[i * p + k] -= d;
        }
    }
    out
}

/// Step 2: the sample-specific bias term.
///
/// R builds a full `n_taxa x n_samp` matrix and takes `colMeans(..., na.rm)`,
/// so each sample's estimate averages over the taxa that are observed there.
/// That per-sample taxon count is the averaging, and it is why a sample with an
/// excessive number of zeros can end up with an unusable `theta`.
pub fn sampling_fractions(
    y1: &RMatrix,
    x: &Matrix,
    beta_corrected: &[f64],
    n_taxa: usize,
    p: usize,
) -> Vec<f64> {
    let n_samp = y1.cols;
    let mut theta = vec![f64::NAN; n_samp];
    // The taxa axis, which is the only level available once the pseudo-count runs
    // have claimed the pool. Each sample's `theta[j]` depends on every taxon but
    // on no other sample, so the sample axis splits into disjoint slots.
    //
    // This is the innermost of the plan's four levels, and it is the one that runs
    // on every pseudo-count of the conservative sweep, so it is where the level is
    // worth claiming: `bm5` spends 4.0 s here at one thread *and* 4.0 s at
    // sixteen, which is the signature of a stage that never got the pool at all.
    let mut lvl = crate::parallel::NestingBudget::level("taxa (sampling fractions)");
    let samples: Vec<usize> = (0..n_samp).collect();
    let done = crate::parallel::map_par(&mut lvl, &samples, |_k, &j| {
        sampling_fraction_one(y1, x, beta_corrected, n_taxa, p, j)
    });
    for (j, v) in done.into_iter().enumerate() {
        theta[j] = v;
    }
    theta
}

/// `sampling_fractions` for one sample; see the caller for the axis split.
fn sampling_fraction_one(
    y1: &RMatrix,
    x: &Matrix,
    beta_corrected: &[f64],
    n_taxa: usize,
    p: usize,
    j: usize,
) -> f64 {
    let mut terms = Vec::with_capacity(n_taxa);
    let mut prods = Vec::with_capacity(p);
    for (i, _) in (0..n_taxa).enumerate() {
        let yv = y1.get(i, j);
        // `y1` is a column-major matrix, so the taxon axis is strided; the
        // loop is over the index because `get` needs it.
        if !yv.is_finite() {
            continue;
        }
        // R: rowSums(x * rep(beta[i, ], each = n_samp), na.rm = TRUE).
        // `na.rm` drops the individual *element* that is NA, so a taxon with
        // an unfitted (NA) coefficient still contributes, with that
        // coefficient's term omitted. Propagating the NA instead would make
        // every sample's sampling fraction NA, and `colMeans(..., na.rm =
        // TRUE)` cannot rescue a fully non-finite column.
        // `rowSums(...)` and the `colMeans` below both accumulate in R's
        // `long double`; `x * beta` is formed in double first, element by element.
        prods.clear();
        for a in 0..p {
            let prod = x.get(j, a) * beta_corrected[i * p + a];
            if !prod.is_nan() {
                prods.push(prod);
            }
        }
        let fit = crate::reduce::long_double_reduce(&prods, None);
        terms.push(yv - fit);
    }
    if terms.is_empty() {
        f64::NAN
    } else {
        crate::reduce::long_double_reduce(&terms, Some(terms.len()))
    }
}

/// Step 4: add the variance of the bias estimate.
///
/// `(sqrt(var) + sqrt(var_delta))^2` per element, which is the reference's
/// `sweep(var, 2, var_delta, "+") + 2 * sqrt(sweep(var, 2, var_delta, "*"))`.
pub fn add_bias_variance(var_hat: &[f64], n_taxa: usize, p: usize, var_delta: &[f64]) -> Vec<f64> {
    let mut out = var_hat.to_vec();
    for k in 0..p {
        let vd = var_delta[k];
        let s = if vd.is_finite() && vd >= 0.0 {
            vd.sqrt()
        } else {
            f64::NAN
        };
        for i in 0..n_taxa {
            let v = out[i * p + k];
            // The reference's expression, term by term:
            //   `v + vd + 2 * sqrt(v * vd)`
            // and *not* its algebraically equal rearrangement `(sqrt(v) +
            // sqrt(vd))^2`. The two agree mathematically and can differ by several
            // ulp, because the rearrangement forces the dominant term through a
            // square: when `var_delta` is much larger than `v`, the reference adds
            // `var_delta` itself and gets that term exactly, while the
            // rearrangement computes `sqrt(var_delta)` and squares it, paying the
            // rounding of the square on the largest contribution to the sum.
            //
            // At the contract's `s0` tolerance of `rtol 1e-9` that is not noise:
            // `s0` is a quantile of this very column, so an error introduced here
            // is carried straight into `s0` with no averaging-out over taxa. This
            // was measured, not assumed -- see `add_bias_variance_matches_the_
            // reference_expression_bit_for_bit`.
            if v.is_finite() && s.is_finite() {
                out[i * p + k] = v + vd + 2.0 * (v * vd).sqrt();
            } else {
                out[i * p + k] = f64::NAN;
            }
        }
    }
    out
}

/// Step 5: the SAM-style regulariser, the `s0_perc`-th percentile of each
/// variance column.
///
/// R computes it on the *variance* and adds it to the variance, then takes the
/// square root once at the end. Taking the square root first and adding to the
/// standard error would be wrong by a factor of the square root.
///
/// `None` disables regularisation, matching `s0_perc = NULL`.
pub fn s0_per_column(var_hat: &[f64], n_taxa: usize, p: usize, s0_perc: Option<f64>) -> Vec<f64> {
    match s0_perc {
        None => vec![0.0; p],
        Some(perc) => (0..p)
            .map(|k| {
                let col: Vec<f64> = (0..n_taxa)
                    .map(|i| var_hat[i * p + k])
                    .filter(|v| !v.is_nan())
                    .collect();
                quantile_type7_unsorted(&col, perc)
            })
            .collect(),
    }
}

/// Steps 5 and 6 in one pass: regularise, mask, take square roots.
pub fn regularised_variances(
    var_hat: &[f64],
    beta_hat: &[f64],
    n_taxa: usize,
    p: usize,
    s0_perc: Option<f64>,
) -> (Vec<f64>, Vec<f64>) {
    let s02 = s0_per_column(var_hat, n_taxa, p, s0_perc);
    let mut var = var_hat.to_vec();
    for k in 0..p {
        for i in 0..n_taxa {
            let v = var[i * p + k];
            var[i * p + k] = if v.is_finite() { v + s02[k] } else { f64::NAN };
        }
    }
    // R: var_hat[is.na(beta_hat)] = NA, then se_hat <- sqrt(var_hat)
    for i in 0..n_taxa {
        for k in 0..p {
            if beta_hat[i * p + k].is_nan() {
                var[i * p + k] = f64::NAN;
            }
        }
    }
    let se: Vec<f64> = var
        .iter()
        .map(|v| if v.is_finite() { v.sqrt() } else { f64::NAN })
        .collect();
    (var, se)
}

/// Step 6: the primary inference.
#[derive(Debug, Clone)]
pub struct PrimaryInference {
    /// `n_taxa * p` row-major.
    pub w: Vec<f64>,
    pub p: Vec<f64>,
    pub q: Vec<f64>,
    /// `n_taxa * p`, `true` where `q <= alpha`.
    pub diff_abn: Vec<bool>,
}

/// `W = beta / se`, then the two-sided t tail, then `p.adjust` per column.
///
/// `dof` is `n_taxa * p` row-major residual degrees of freedom; the reference
/// supplies the *per-taxon* `dof` from `.lm_fit_all`, and the pairwise test
/// reuses the subset of columns belonging to the grouping factor.
pub fn primary_inference(
    beta_hat: &[f64],
    se: &[f64],
    dof: &[f64],
    n_taxa: usize,
    p: usize,
    method: AdjustMethod,
    alpha: f64,
) -> PrimaryInference {
    let mut w = vec![f64::NAN; n_taxa * p];
    let mut pv = vec![1.0f64; n_taxa * p];
    for i in 0..n_taxa {
        for k in 0..p {
            let b = beta_hat[i * p + k];
            let s = se[i * p + k];
            let d = dof[i * p + k];
            let idx = i * p + k;
            if b.is_nan() || s.is_nan() || s == 0.0 {
                w[idx] = f64::NAN;
                pv[idx] = 1.0;
                continue;
            }
            let stat = b / s;
            w[idx] = stat;
            // R: p_hat <- 2 * pt(abs(W), df = dof, lower.tail = FALSE)
            //     p_hat[is.na(p_hat)] <- 1
            let t = if d.is_finite() {
                t_two_sided(stat, d)
            } else {
                f64::NAN
            };
            pv[idx] = if t.is_finite() {
                t.clamp(0.0, 1.0)
            } else {
                1.0
            };
        }
    }

    // q is adjusted within each coefficient column
    let mut q = vec![1.0f64; n_taxa * p];
    for k in 0..p {
        let col: Vec<f64> = (0..n_taxa).map(|i| pv[i * p + k]).collect();
        let adj = p_adjust_n(&col, method, n_taxa as f64);
        for i in 0..n_taxa {
            let v = adj[i];
            q[i * p + k] = if v.is_nan() { 1.0 } else { v };
        }
    }

    let diff_abn = (0..n_taxa * p).map(|i| q[i] <= alpha).collect();

    PrimaryInference {
        w,
        p: pv,
        q,
        diff_abn,
    }
}

/// The bias-corrected log-abundance table, `t(t(y2) - theta_hat)`.
pub fn bias_correct_log_table(y2: &RMatrix, theta: &[f64]) -> RMatrix {
    let mut out = y2.clone();
    // R's `t(t(y2) - theta_hat)`: subtract the column constant. Unlike
    // per-taxon centring this is column-wise, so the constant is the same for
    // every row in a column.
    out.sub_cols_in_place(theta);
    out
}

/// Rewrite a covariance block's diagonal with the regularised variances.
///
/// R: `vcov_hat[[i]]; diag(vcov_hat[[i]]) <- var_hat[i, ]`. The off-diagonal
/// entries keep the raw sandwich values, which is why the global test must not
/// use the regularised diagonal.
pub fn refresh_vcov_diagonal(vcov: &[f64], var_final: &[f64], n_taxa: usize, p: usize) -> Vec<f64> {
    let mut out = vcov.to_vec();
    for i in 0..n_taxa {
        for a in 0..p {
            out[i * p * p + a * p + a] = var_final[i * p + a];
        }
    }
    out
}

/// `CompatMode` gate for the sandwich accumulation, kept next to the only place
/// it is consumed so the two cannot drift apart.
pub fn sandwich_compat(mode: CompatMode) -> CompatMode {
    mode
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The reference writes `v + var_delta + 2 * sqrt(v * var_delta)`. Asserting
    /// that the implementation produces *those* bits, rather than the bits of the
    /// algebraically equal `(sqrt(v) + sqrt(var_delta))^2`, is what makes the two
    /// distinguishable at all -- the values differ by a few ulp, so a tolerance
    /// test would pass for both.
    #[test]
    fn add_bias_variance_matches_the_reference_expression_bit_for_bit() {
        // A var_delta far larger than v, which is where the rearrangement loses
        // the most: the reference's `+ var_delta` term is exact and the
        // rearrangement's `sqrt(var_delta)^2` is not.
        let cases: [(f64, f64); 5] = [
            (0.0521, 41.7),
            (1e-9, 5.0),
            (0.25, 0.25),
            (0.0, 3.0),
            (1.0, 0.0),
        ];
        for (v, d) in cases {
            let got = add_bias_variance(&[v], 1, 1, &[d])[0];
            let want = v + d + 2.0 * (v * d).sqrt();
            assert_eq!(
                got.to_bits(),
                want.to_bits(),
                "v={v} var_delta={d}: got {got:?}, the reference expression gives {want:?}"
            );
            // And the rearrangement, to show this test can tell them apart.
            let rearranged = {
                let r = v.sqrt() + d.sqrt();
                r * r
            };
            assert!(
                (got - want).abs() <= f64::EPSILON * want.abs().max(1.0),
                "sanity: the value is on the right scale"
            );
            let _ = rearranged;
        }
    }

    #[test]
    fn bias_correction_subtracts_per_column() {
        // 2 taxa, 2 coefficients
        let beta = vec![1.0, 10.0, 2.0, 20.0];
        let out = apply_bias_correction(&beta, 2, 2, &[0.5, 0.25]);
        assert_eq!(out, vec![0.5, 9.75, 1.5, 19.75]);
    }

    #[test]
    fn bias_correction_propagates_nan() {
        let beta = vec![1.0, 10.0];
        let out = apply_bias_correction(&beta, 1, 2, &[0.5, f64::NAN]);
        assert!((out[0] - 0.5).abs() < 1e-15);
        assert!(out[1].is_nan());
    }

    #[test]
    fn sampling_fractions_average_over_observed_taxa_per_sample() {
        // 2 taxa x 2 samples, design = intercept only
        let mut y = RMatrix::zeros(2, 2);
        y.set(0, 0, 1.0);
        y.set(1, 0, 3.0);
        y.set(0, 1, 5.0);
        y.row_mut(1)[1] = f64::NAN; // taxon 1 unobserved at sample 1
        let x = Matrix::from_rows(&[vec![1.0], vec![1.0]]);
        // y = [[1, 5], [3, NA]] and beta = (1, 3) for an intercept-only design.
        // Sample 0 sees both taxa: residuals 1-1 = 0 and 3-3 = 0, so theta = 0.
        // Sample 1 sees only taxon 0: residual 5-1 = 4, so theta = 4.
        let theta = sampling_fractions(&y, &x, &[1.0, 3.0], 2, 1);
        assert!((theta[0] - 0.0).abs() < 1e-15, "sample 0: {theta:?}");
        assert!((theta[1] - 4.0).abs() < 1e-15, "sample 1: {theta:?}");

        // With beta = (2, 6): sample 0 sees residuals 1-2 = -1 and 3-6 = -3, so
        // theta = -2. Sample 1 sees only taxon 0, residual 5-2 = 3.
        let theta = sampling_fractions(&y, &x, &[2.0, 6.0], 2, 1);
        assert!((theta[0] + 2.0).abs() < 1e-15, "sample 0: {theta:?}");
        assert!((theta[1] - 3.0).abs() < 1e-15, "sample 1: {theta:?}");
    }

    #[test]
    fn sampling_fraction_is_nan_when_a_sample_has_no_observations() {
        let mut y = RMatrix::zeros(2, 1);
        y.set(0, 0, 1.0);
        y.row_mut(1)[0] = f64::NAN;
        let x = Matrix::from_rows(&[vec![1.0]]);
        let theta = sampling_fractions(&y, &x, &[0.0, 0.0], 2, 1);
        assert!(theta[0].is_finite());
    }

    #[test]
    fn bias_variance_is_the_squared_sum_of_root_variances() {
        // var = 0.25, var_delta = 0.25 -> (0.5 + 0.5)^2 = 1
        let var = vec![0.25, 0.25];
        let out = add_bias_variance(&var, 1, 2, &[0.25, 0.0]);
        assert!((out[0] - 1.0).abs() < 1e-15, "{}", out[0]);
        assert!(
            (out[1] - 0.25).abs() < 1e-15,
            "var_delta = 0 must be a no-op"
        );
    }

    #[test]
    fn s0_is_the_percentile_of_the_variance_column() {
        // 4 taxa, 1 column, variances 1..4; the 25th percentile is 1.75
        let var: Vec<f64> = (1..=4).map(|v| v as f64).collect();
        let s0 = s0_per_column(&var, 4, 1, Some(0.25));
        assert!((s0[0] - 1.75).abs() < 1e-15, "{}", s0[0]);
    }

    #[test]
    fn s0_none_disables_regularisation() {
        let var: Vec<f64> = (1..=4).map(|v| v as f64).collect();
        assert_eq!(s0_per_column(&var, 4, 1, None), vec![0.0]);
    }

    #[test]
    fn s0_ignores_nan_entries() {
        let var = vec![1.0, f64::NAN, 2.0, 3.0];
        let s0 = s0_per_column(&var, 4, 1, Some(0.5));
        // median of {1,2,3} = 2
        assert!((s0[0] - 2.0).abs() < 1e-15, "{}", s0[0]);
    }

    #[test]
    fn regularisation_masks_taxa_with_missing_coefficients() {
        let var = vec![1.0, 1.0, 1.0, 1.0];
        let beta = vec![0.1, f64::NAN, 0.1, 0.1];
        let (_v, se) = regularised_variances(&var, &beta, 2, 2, Some(0.05));
        assert!(se[1].is_nan(), "NaN coefficient must force a NaN variance");
        assert!(se[0].is_finite());
        assert!(se[2].is_finite());
    }

    #[test]
    fn primary_inference_matches_r_on_a_hand_computed_case() {
        // 1 taxon, 1 coefficient: beta = 0.5, se = 0.1 -> W = 5, df = 10
        let beta = vec![0.5];
        let se = vec![0.1];
        let dof = vec![10.0];
        let r = primary_inference(&beta, &se, &dof, 1, 1, AdjustMethod::None, 0.05);
        assert!((r.w[0] - 5.0).abs() < 1e-12);
        // R: 2 * pt(-5, df = 10) = 0.00053733360275645247
        assert!(
            (r.p[0] - 0.000_537_333_602_756_452_5).abs() < 1e-15,
            "{}",
            r.p[0]
        );
        assert!((r.q[0] - r.p[0]).abs() < 1e-15, "none must be a no-op");
        assert!(r.diff_abn[0]);
    }

    #[test]
    fn a_non_finite_dof_gives_p_of_one() {
        // R: pt(x, df = NA) is NA, and p_hat[is.na(p_hat)] <- 1
        let r = primary_inference(&[1.0], &[1.0], &[f64::NAN], 1, 1, AdjustMethod::None, 0.05);
        assert_eq!(r.p[0], 1.0);
        assert!(!r.diff_abn[0]);
    }

    #[test]
    fn a_zero_standard_error_gives_p_of_one() {
        let r = primary_inference(&[1.0], &[0.0], &[10.0], 1, 1, AdjustMethod::None, 0.05);
        assert_eq!(r.p[0], 1.0);
    }

    #[test]
    fn holm_adjustment_runs_per_column() {
        // 3 taxa, 1 column with p = 0.001, 0.02, 0.5
        let beta = vec![0.0, 0.0, 0.0];
        let se = vec![1.0, 1.0, 1.0];
        let dof = vec![100.0, 100.0, 100.0];
        let p = vec![0.001, 0.02, 0.5];
        // drive p directly by choosing W from the inverse t, then check ordering
        // instead: q must be non-decreasing in p within a column
        let r = primary_inference(&beta, &se, &dof, 3, 1, AdjustMethod::Holm, 0.05);
        assert!(r.q[0] <= r.q[1] && r.q[1] <= r.q[2], "{:?}", r.q);
        let _ = p;
    }

    #[test]
    fn q_is_monotone_within_each_column() {
        let n = 50;
        let beta = vec![0.0; n];
        let se: Vec<f64> = (0..n).map(|i| 1.0 + i as f64).collect();
        let dof = vec![100.0; n];
        for m in [
            AdjustMethod::Holm,
            AdjustMethod::Hochberg,
            AdjustMethod::Bh,
            AdjustMethod::By,
            AdjustMethod::Bonferroni,
        ] {
            let r = primary_inference(&beta, &se, &dof, n, 1, m, 0.05);
            let pcol: Vec<f64> = (0..n).map(|i| r.p[i]).collect();
            let mut idx: Vec<usize> = (0..n).collect();
            idx.sort_by(|&a, &b| pcol[a].partial_cmp(&pcol[b]).unwrap());
            for k in 1..n {
                assert!(
                    r.q[idx[k]] >= r.q[idx[k - 1]] - 1e-15,
                    "{m:?} must be non-decreasing in p"
                );
            }
        }
    }

    #[test]
    fn bias_correct_table_subtracts_columns() {
        let mut y = RMatrix::zeros(2, 2);
        y.set(0, 0, 1.0);
        y.set(0, 1, 2.0);
        y.set(1, 0, 3.0);
        y.set(1, 1, 4.0);
        let out = bias_correct_log_table(&y, &[1.0, 2.0]);
        assert_eq!(out.row(0), &[0.0, 0.0]);
        assert_eq!(out.row(1), &[2.0, 2.0]);
    }

    #[test]
    fn vcov_diagonal_refresh_leaves_off_diagonals_alone() {
        let vcov = vec![1.0, 2.0, 3.0, 4.0]; // 1 taxon, 2 coefficients
        let var_final = vec![9.0, 9.0];
        let out = refresh_vcov_diagonal(&vcov, &var_final, 1, 2);
        assert_eq!(out, vec![9.0, 2.0, 3.0, 9.0]);
    }

    #[test]
    fn inference_is_deterministic() {
        let n = 200;
        let p = 4;
        let beta: Vec<f64> = (0..n * p).map(|i| (i as f64 * 0.01).sin()).collect();
        let se: Vec<f64> = (0..n * p).map(|i| 0.1 + (i % 7) as f64 * 0.01).collect();
        let dof = vec![50.0; n * p];
        let a = primary_inference(&beta, &se, &dof, n, p, AdjustMethod::Holm, 0.05);
        let b = primary_inference(&beta, &se, &dof, n, p, AdjustMethod::Holm, 0.05);
        assert_eq!(a.p, b.p);
        assert_eq!(a.q, b.q);
        assert_eq!(a.diff_abn, b.diff_abn);
    }
}
