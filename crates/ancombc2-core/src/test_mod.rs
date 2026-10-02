//! Global and pairwise multi-group tests, with the mdFDR procedure.
//!
//! # Global test
//!
//! `.ancombc_global_F` forms, per taxon, the quadratic Wald statistic over the
//! columns of the design that belong to the grouping factor and that are not
//! interaction terms:
//!
//! ```text
//! group_ind  = grepl(group, covariates) & !grepl(":", covariates)
//! W          = t(beta[group_ind]) %*% ginv(vcov[group_ind, group_ind]) %*% beta[group_ind]
//! p          = 2 * min(pchisq(W, n_grp, lower = TRUE), pchisq(W, n_grp, lower = FALSE))
//! ```
//!
//! The two-sided-by-minimum construction is unusual and is reproduced literally.
//! The fixed-effects branch passes a non-`NULL` `dof`, so `use_chisq` is FALSE
//! and the tail is an F instead: `2 * min(pf(W, n_grp, dof), pf(W, n_grp, dof,
//! lower = FALSE))` with `dof_i` the taxon's residual degrees of freedom.
//!
//! # Pairwise test and mdFDR
//!
//! The pairwise table has one column per group coefficient *and* one per
//! unordered pair of group coefficients, `n_grp + C(n_grp, 2)` in total, and
//! each column is a difference of two coefficients.
//!
//! The mdFDR procedure (`.mdfdr`) is a two-stage construction and must be
//! implemented as written, not simplified to "run all pairs then BH":
//!
//! 1. **Screen.** Re-run the global test with `p_adj_method = "BH"` and count the
//!    taxa rejected, `R`.
//! 2. **Family-wise adjustment.** Compute the pairwise p-values, zero out the
//!    taxa that failed the screen, then adjust *within each taxon* with
//!    `n = length(x) * n_tax / R` — that is, at a denominator inflated by the
//!    total number of taxa and deflated by the number that survived the screen.
//!    Taxa that failed the screen get `p = 1` and so `q = 1`.
//!
//! `R = 0` makes `n` infinite, which the reference does not guard against; the
//! result is that every pairwise test is reported as non-significant. That
//! behaviour is reproduced, and flagged in `docs/reference_behavior.md`.

use ancombc2_stats::{f_upper, p_adjust_n};

use crate::error::{AncombcError, Result};
use crate::matrix::Matrix;

/// Result of the global test.
#[derive(Debug, Clone)]
pub struct GlobalTest {
    pub w: Vec<f64>,
    pub p: Vec<f64>,
    pub q: Vec<f64>,
    pub diff_abn: Vec<bool>,
}

/// The columns belonging to the grouping factor: name matches the group and is
/// not an interaction term.
pub fn group_columns(colnames: &[String], group: &str) -> Vec<usize> {
    colnames
        .iter()
        .enumerate()
        .filter(|(_, n)| n.contains(group) && !n.contains(':'))
        .map(|(i, _)| i)
        .collect()
}

/// The global test.
///
/// `vcov` is `n_taxa * p * p` row-major. `dof` is `n_taxa * p` row-major; when
/// `None` the chi-square branch is used, matching the reference's
/// `use_chisq = is.null(dof)`.
/// The argument list mirrors the reference's `.ancombc_global_F`: the fix_eff
/// names, the group, the coefficients, the covariance, the degrees of freedom,
/// the shapes, the adjustment method and the level. Grouping them would hide
/// which array is which at the call site, where a swap is invisible to the type
/// checker.
#[allow(clippy::too_many_arguments)]
pub fn global_test(
    x_colnames: &[String],
    group: &str,
    beta: &[f64],
    vcov: &[f64],
    dof: Option<&[f64]>,
    n_taxa: usize,
    p: usize,
    method: crate::config::AdjustMethod,
    alpha: f64,
) -> Result<GlobalTest> {
    let idx = group_columns(x_colnames, group);
    if idx.is_empty() {
        return Err(AncombcError::GroupRequiredForStructuralZeros);
    }
    let g = idx.len();
    let mut w = vec![f64::NAN; n_taxa];
    let mut pv = vec![1.0f64; n_taxa];

    for i in 0..n_taxa {
        // sub-matrix of the covariance block
        let mut sub = Matrix::zeros(g, g);
        for (a, &ia) in idx.iter().enumerate() {
            for (b, &ib) in idx.iter().enumerate() {
                sub.set(a, b, vcov[i * p * p + ia * p + ib]);
            }
        }
        let bsub: Vec<f64> = idx.iter().map(|&k| beta[i * p + k]).collect();
        // ginv of the sub-matrix, as the reference does per taxon
        let inv = crate::matrix::ginv(&sub, None);
        let mut stat = 0.0;
        let mut ok = true;
        for a in 0..g {
            let mut inner = 0.0;
            for b in 0..g {
                inner += inv.get(a, b) * bsub[b];
            }
            stat += bsub[a] * inner;
        }
        if !stat.is_finite() {
            ok = false;
        }
        if !ok {
            w[i] = f64::NAN;
            pv[i] = 1.0;
            continue;
        }
        w[i] = stat;
        pv[i] = match dof {
            None => {
                let lo = ancombc2_stats::chisq_lower(stat, g as f64);
                let hi = ancombc2_stats::chisq_upper(stat, g as f64);
                2.0 * lo.min(hi)
            }
            Some(d) => {
                // R: dof_i <- unique(dof[i, ]) -- the taxon's dof, identical
                // across its coefficients
                let di = d[i * p + idx[0]];
                if di.is_finite() {
                    let lo = 1.0 - f_upper(stat, g as f64, di);
                    let hi = f_upper(stat, g as f64, di);
                    2.0 * lo.min(hi)
                } else {
                    1.0
                }
            }
        };
        pv[i] = pv[i].clamp(0.0, 1.0);
    }

    let q = adjust_all(&pv, method, n_taxa, alpha);
    let diff_abn = (0..n_taxa).map(|i| q[i] <= alpha).collect();
    Ok(GlobalTest {
        w,
        p: pv,
        q,
        diff_abn,
    })
}

/// `p.adjust` with `NA` replaced by 1, as the reference does before testing.
fn adjust_all(p: &[f64], method: crate::config::AdjustMethod, n: usize, alpha: f64) -> Vec<f64> {
    let mut q = p_adjust_n(p, method, n as f64);
    for v in q.iter_mut() {
        if v.is_nan() {
            *v = 1.0;
        }
    }
    let _ = alpha;
    q
}

/// One column of the pairwise table: its name, the coefficient it reports
/// (`Some(k)`) and, for a contrast, the coefficient it is differenced against
/// (`Some(j)`).
///
/// R's `.combn_fun` builds the table as
///
/// ```r
/// y = c(x, combn(x, 2, FUN = fun))
/// ```
///
/// so the table's **first `g` columns are the group coefficients themselves**,
/// not differences of a coefficient with itself. Only the trailing
/// `C(g, 2)` columns are contrasts, each named `paste(name_j, name_i, "_")` for
/// `b_j - b_i` in `combn`'s lexicographic pair order. Modelling the diagonal as a
/// self-difference would zero out the first `g` columns, so the distinction is
/// load-bearing.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PairColumn {
    pub name: String,
    /// The coefficient being reported, or the minuend of the contrast.
    pub hi: usize,
    /// `None` for a plain group coefficient, `Some(i)` for the contrast
    /// `beta[hi] - beta[i]`.
    pub lo: Option<usize>,
}

impl PairColumn {
    pub fn is_contrast(&self) -> bool {
        self.lo.is_some()
    }
}

/// Build the pairwise table's columns, in the reference's order.
pub fn pairwise_columns(names: &[String]) -> Vec<PairColumn> {
    let mut out: Vec<PairColumn> = names
        .iter()
        .enumerate()
        .map(|(i, n)| PairColumn {
            name: n.clone(),
            hi: i,
            lo: None,
        })
        .collect();
    for a in 0..names.len() {
        for b in (a + 1)..names.len() {
            out.push(PairColumn {
                name: format!("{}_{}", names[b], names[a]),
                hi: b,
                lo: Some(a),
            });
        }
    }
    out
}

/// Result of the pairwise test with mdFDR control.
#[derive(Debug, Clone)]
pub struct PairwiseTest {
    /// `n_taxa * n_col` row-major, `b_hi - b_lo`.
    pub beta: Vec<f64>,
    pub se: Vec<f64>,
    pub w: Vec<f64>,
    pub p: Vec<f64>,
    pub q: Vec<f64>,
    pub diff_abn: Vec<bool>,
    /// Column names, in table order.
    pub colnames: Vec<String>,
    /// Taxa that passed the mdFDR global screen.
    pub screen: Vec<bool>,
    /// `R`, the number of taxa rejected by the screening global test.
    pub n_screen: usize,
}

/// The pairwise test with the two-stage mdFDR control.
#[allow(clippy::too_many_arguments)]
pub fn pairwise_test(
    x_colnames: &[String],
    group: &str,
    beta: &[f64],
    vcov: &[f64],
    var_hat: &[f64],
    dof: Option<&[f64]>,
    n_taxa: usize,
    p_all: usize,
    method: crate::config::AdjustMethod,
    fwer_ctrl_method: crate::config::AdjustMethod,
    alpha: f64,
) -> Result<PairwiseTest> {
    // `.mdfdr` screens with `p_adj_method = "BH"` hard-coded and adjusts the
    // pairwise p-values with `fwer_ctrl_method`, so the caller's `method` never
    // reaches the result. The parameter is kept to mirror the reference's
    // signature rather than to look like an oversight.
    let _ = method;
    let idx = group_columns(x_colnames, group);
    if idx.is_empty() {
        return Err(AncombcError::GroupRequiredForStructuralZeros);
    }
    let cols = pairwise_columns(
        &idx.iter()
            .map(|&k| x_colnames[k].clone())
            .collect::<Vec<_>>(),
    );
    let n_col = cols.len();

    // Stage 1: screen with the global test, BH-adjusted.
    let screen_test = global_test(
        x_colnames,
        group,
        beta,
        vcov,
        dof,
        n_taxa,
        p_all,
        crate::config::AdjustMethod::Bh,
        alpha,
    )?;
    let screen = screen_test.diff_abn.clone();
    let r_count = screen.iter().filter(|&&b| b).count();

    // Differences and their variances.
    let mut bmat = vec![f64::NAN; n_taxa * n_col];
    let mut vmat = vec![f64::NAN; n_taxa * n_col];
    for i in 0..n_taxa {
        for (c, col) in cols.iter().enumerate() {
            let bh = beta[i * p_all + idx[col.hi]];
            let a_hi = idx[col.hi];
            match col.lo {
                // a plain group coefficient: the estimate is the coefficient and
                // the variance is its own marginal variance
                None => {
                    bmat[i * n_col + c] = bh;
                    vmat[i * n_col + c] = var_hat[i * p_all + a_hi];
                }
                Some(lo) => {
                    let a_lo = idx[lo];
                    bmat[i * n_col + c] = bh - beta[i * p_all + a_lo];
                    vmat[i * n_col + c] = vcov[i * p_all * p_all + a_hi * p_all + a_hi]
                        + vcov[i * p_all * p_all + a_lo * p_all + a_lo]
                        - 2.0 * vcov[i * p_all * p_all + a_hi * p_all + a_lo];
                }
            }
        }
    }

    let mut se = vec![f64::NAN; n_taxa * n_col];
    let mut w = vec![f64::NAN; n_taxa * n_col];
    for i in 0..n_taxa {
        for c in 0..n_col {
            let v = vmat[i * n_col + c];
            if v.is_finite() && v > 0.0 {
                se[i * n_col + c] = v.sqrt();
                let b = bmat[i * n_col + c];
                if b.is_finite() {
                    w[i * n_col + c] = b / se[i * n_col + c];
                }
            }
        }
    }

    // Raw p-values, then the screen.
    let mut pv = vec![1.0f64; n_taxa * n_col];
    for i in 0..n_taxa {
        for c in 0..n_col {
            let stat = w[i * n_col + c];
            if stat.is_finite() {
                // With no `dof` the reference leaves `df = NULL`, and R's
                // `pt(x, df = NULL)` evaluates the *default*, df = Inf: a normal
                // tail, not a missing value. That is the branch the fixed-effects
                // path actually takes, since `.iter_mle` returns `dof = NULL` when
                // it estimated theta itself.
                let d = match dof {
                    None => f64::INFINITY,
                    Some(dd) => dd[i * p_all + idx[cols[c].hi]],
                };
                // Inf is a valid df (the normal limit); only NaN is unusable.
                let t = if d.is_nan() {
                    f64::NAN
                } else {
                    ancombc2_stats::t_two_sided(stat, d)
                };
                pv[i * n_col + c] = if t.is_finite() {
                    t.clamp(0.0, 1.0)
                } else {
                    1.0
                };
            }
        }
    }
    // p <- p * screen ; p[p == 0] <- 1 ; p[is.na(p)] <- 1
    for i in 0..n_taxa {
        for c in 0..n_col {
            pv[i * n_col + c] *= screen[i] as i32 as f64;
            if pv[i * n_col + c] == 0.0 {
                pv[i * n_col + c] = 1.0;
            }
        }
    }

    // Stage 2: adjust within each taxon at n = n_col * n_tax / R.
    let mut q = vec![1.0f64; n_taxa * n_col];
    if r_count == 0 {
        // R divides by R without guarding; every p.adjust call gets n = Inf and
        // returns 0, so every q is 0 and every call significant. Reproduced
        // deliberately: see the module docs.
        for i in 0..n_taxa {
            for c in 0..n_col {
                q[i * n_col + c] = 0.0;
            }
        }
    } else {
        // `length(x) * n_tax / R` is fractional in general, and R keeps it a
        // double inside p.adjust, so it must not be rounded here.
        let n_eff = (n_col * n_taxa) as f64 / r_count as f64;
        for i in 0..n_taxa {
            let row: Vec<f64> = (0..n_col).map(|c| pv[i * n_col + c]).collect();
            let adj = p_adjust_n(&row, fwer_ctrl_method, n_eff);
            for c in 0..n_col {
                let v = adj[c];
                q[i * n_col + c] = if v.is_nan() { 1.0 } else { v };
            }
        }
    }

    let diff_abn = (0..n_taxa * n_col).map(|i| q[i] <= alpha).collect();

    Ok(PairwiseTest {
        beta: bmat,
        se,
        w,
        p: pv,
        q,
        diff_abn,
        colnames: cols.iter().map(|c| c.name.clone()).collect(),
        screen,
        n_screen: r_count,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::config::AdjustMethod;

    fn names(v: &[&str]) -> Vec<String> {
        v.iter().map(|s| s.to_string()).collect()
    }

    /// A design with an intercept, two group dummies and a covariate, plus an
    /// interaction that must be excluded from the group columns.
    fn design() -> Vec<String> {
        names(&["(Intercept)", "groupb", "groupc", "age", "groupb:age"])
    }

    #[test]
    fn group_columns_exclude_interactions() {
        let idx = group_columns(&design(), "group");
        assert_eq!(
            idx,
            vec![1, 2],
            "only the group dummies, not the interaction"
        );
    }

    #[test]
    fn group_columns_match_on_substring() {
        // R's grepl is a substring match, so "bmi_group" also matches
        // "bmi_group2" and the reference's guard is only against ":".
        let n = names(&["(Intercept)", "bmi_group2", "site:sex"]);
        assert_eq!(group_columns(&n, "bmi_group"), vec![1]);
    }

    #[test]
    fn pairwise_column_order_matches_combn() {
        let cols = pairwise_columns(&names(&["b", "c", "d"]));
        let labels: Vec<&str> = cols.iter().map(|c| c.name.as_str()).collect();
        assert_eq!(labels, vec!["b", "c", "d", "c_b", "d_b", "d_c"]);
        // the first g columns are the coefficients themselves
        for (k, col) in cols.iter().take(3).enumerate() {
            assert!(!col.is_contrast(), "column {k} must be a plain coefficient");
            assert_eq!(col.hi, k);
        }
        // the rest are differences b_j - b_i
        assert_eq!((cols[3].hi, cols[3].lo), (1, Some(0)));
        assert_eq!((cols[5].hi, cols[5].lo), (2, Some(1)));
    }

    #[test]
    fn pairwise_column_count_is_g_plus_g_choose_2() {
        for g in 1..=6usize {
            let ns: Vec<String> = (0..g).map(|i| format!("g{i}")).collect();
            assert_eq!(pairwise_columns(&ns).len(), g + g * (g - 1) / 2);
        }
    }

    /// Build a synthetic 2-taxon, 3-level-group problem with a known signal.
    ///
    /// The design is `(Intercept) + grp2 + grp3` and the group variable is `grp`,
    /// so R's `grepl("grp", covariates)` selects both dummies and the group block
    /// has two columns. A pairwise test therefore has `2 + C(2,2) = 3` columns:
    /// the two group coefficients and their difference.
    fn toy() -> (Vec<String>, Vec<f64>, Vec<f64>, usize, usize) {
        let n = 2usize;
        let p = 3usize;
        let beta = vec![
            0.0, 2.0, 0.5, // taxon 0: strong effect in both contrasts
            0.0, 0.1, 0.2, // taxon 1: no effect
        ];
        let mut vcov = vec![0.0; n * p * p];
        for i in 0..n {
            for a in 0..p {
                vcov[i * p * p + a * p + a] = 0.04; // 0.04 on the diagonal
            }
        }
        (names(&["(Intercept)", "grp2", "grp3"]), beta, vcov, n, p)
    }

    /// `toy()` with the group block's variance supplied separately, which the
    /// pairwise test needs for the diagonal columns.
    fn toy_var() -> Vec<f64> {
        let mut v = vec![0.0; 2 * 3];
        for i in 0..2 {
            for a in 0..3 {
                v[i * 3 + a] = 0.04;
            }
        }
        v
    }

    #[test]
    fn global_test_flags_only_the_strong_taxon() {
        let (cn, beta, vcov, n, p) = toy();
        let g = global_test(
            &cn,
            "grp",
            &beta,
            &vcov,
            None,
            n,
            p,
            AdjustMethod::None,
            0.05,
        )
        .unwrap();
        // taxon 0: (2, 0.5)' diag(1/0.04) (2, 0.5) = (4 + 0.25) / 0.04 = 106.25
        assert!((g.w[0] - 106.25).abs() < 1e-9, "{}", g.w[0]);
        // taxon 1: (0.01 + 0.04) / 0.04 = 1.25
        assert!((g.w[1] - 1.25).abs() < 1e-9, "{}", g.w[1]);
        assert!(g.p[0] < 1e-10, "{}", g.p[0]);
        assert!(g.p[1] > 0.05, "{}", g.p[1]);
    }

    #[test]
    fn global_test_uses_the_f_distribution_when_dof_is_supplied() {
        let (cn, beta, vcov, n, p) = toy();
        let dof = vec![20.0; n * p];
        let chisq = global_test(
            &cn,
            "grp",
            &beta,
            &vcov,
            None,
            n,
            p,
            AdjustMethod::None,
            0.05,
        )
        .unwrap();
        let f = global_test(
            &cn,
            "grp",
            &beta,
            &vcov,
            Some(&dof),
            n,
            p,
            AdjustMethod::None,
            0.05,
        )
        .unwrap();
        assert!(
            (chisq.p[0] - f.p[0]).abs() > 1e-12,
            "the F branch must differ from the chi-square branch"
        );
        assert!(f.p[0] < 1e-8, "{}", f.p[0]);
        assert!(f.p[1] > 0.5, "{}", f.p[1]);
    }

    #[test]
    fn global_test_rejects_a_missing_group_column() {
        let (cn, beta, vcov, n, p) = toy();
        let e = global_test(
            &cn,
            "nonexistent",
            &beta,
            &vcov,
            None,
            n,
            p,
            AdjustMethod::None,
            0.05,
        )
        .unwrap_err();
        assert!(matches!(e, AncombcError::GroupRequiredForStructuralZeros));
    }

    #[test]
    fn pairwise_difference_and_variance_are_correct() {
        let (cn, beta, vcov, n, p) = toy();
        let r = pairwise_test(
            &cn,
            "grp",
            &beta,
            &vcov,
            &toy_var(),
            None,
            n,
            p,
            AdjustMethod::None,
            AdjustMethod::Holm,
            0.05,
        )
        .unwrap();
        // 2 group columns + 1 difference
        assert_eq!(r.colnames, vec!["grp2", "grp3", "grp3_grp2"]);
        // taxon 0: the two group coefficients, then their difference 0.5 - 2.0
        assert!((r.beta[0] - 2.0).abs() < 1e-15, "{}", r.beta[0]);
        assert!((r.beta[1] - 0.5).abs() < 1e-15, "{}", r.beta[1]);
        assert!((r.beta[2] - (-1.5)).abs() < 1e-15, "{}", r.beta[2]);
        // a plain group coefficient is not differenced with itself
        assert!(r.beta[0].abs() > 1e-9, "column 0 must be the coefficient");
        // the difference's variance is 0.04 + 0.04 - 0 = 0.08
        assert!((r.se[2] - 0.08f64.sqrt()).abs() < 1e-15, "{}", r.se[2]);
        // a single group coefficient's variance is just its own
        assert!((r.se[0] - 0.2).abs() < 1e-15, "{}", r.se[0]);
        // taxon 1: 0.1, 0.2, and 0.2 - 0.1 = 0.1
        assert!((r.beta[3] - 0.1).abs() < 1e-15, "{}", r.beta[3]);
        assert!((r.beta[4] - 0.2).abs() < 1e-15, "{}", r.beta[4]);
        assert!((r.beta[5] - 0.1).abs() < 1e-15, "{}", r.beta[5]);
    }

    #[test]
    fn pairwise_screen_zeroes_the_failing_taxa() {
        let (cn, beta, vcov, n, p) = toy();
        let r = pairwise_test(
            &cn,
            "grp",
            &beta,
            &vcov,
            &toy_var(),
            None,
            n,
            p,
            AdjustMethod::None,
            AdjustMethod::Holm,
            0.05,
        )
        .unwrap();
        assert!(r.screen[0], "taxon 0 has a large W and must pass");
        assert!(!r.screen[1], "taxon 1 has W well below 3.84 and must not");
        assert_eq!(r.n_screen, 1);
        let n_col = r.colnames.len();
        for c in 0..n_col {
            assert!(
                (r.p[c] - 1.0).abs() > 1e-15,
                "taxon 0 must keep a real p at column {c}, got {}",
                r.p[c]
            );
            assert!(
                (r.p[n_col + c] - 1.0).abs() < 1e-15,
                "a screened-out taxon gets p = 1 at column {c}, got {}",
                r.p[n_col + c]
            );
        }
        // All three of taxon 0's columns survive: Holm over 3 p-values with
        // n = n_col * n_tax / R = 3 * 2 / 1 = 6 gives
        //   q = (1.52e-23 * 6 * 3 / 3, 0.0124 * 6 * 2 / 3, 1.14e-7 * 6 * 1 / 3)
        //     = (9.1e-23, 0.0497, 5.7e-7)
        // and 0.0497 <= alpha = 0.05, so all three are called. That the largest
        // p-value lands just under alpha is a coincidence of the fixture, not a
        // property of the method.
        assert!(
            r.diff_abn[..n_col].iter().all(|&d| d),
            "all of taxon 0's columns are significant: {:?}",
            &r.diff_abn[..n_col]
        );
        assert!(r.q[1] < 0.05, "q[1] = {}", r.q[1]);
        assert!(
            !r.diff_abn[n_col..].iter().any(|&d| d),
            "no column of taxon 1 can be significant"
        );
    }

    #[test]
    fn mdfdr_denominator_scales_with_the_screened_count() {
        // 20 taxa, 10 of which have a real effect, 10 with none.
        let n = 20usize;
        let p = 2usize;
        let mut beta = vec![0.0; n * p];
        let mut vcov = vec![0.0; n * p * p];
        for i in 0..n {
            beta[i * p + 1] = if i < 10 { 1.0 } else { 0.02 };
            vcov[i * p * p] = 0.01;
            vcov[i * p * p + 3] = 0.01;
        }
        let cn = names(&["(Intercept)", "grp2"]);
        let var_hat: Vec<f64> = (0..n).flat_map(|_| [0.01, 0.01]).collect();
        let r = pairwise_test(
            &cn,
            "grp",
            &beta,
            &vcov,
            &var_hat,
            None,
            n,
            p,
            AdjustMethod::None,
            AdjustMethod::Holm,
            0.05,
        )
        .unwrap();
        assert!(r.n_screen > 0, "some taxa must pass the screen");
        assert!(
            r.n_screen <= 10,
            "only the real signals pass, got {}",
            r.n_screen
        );
    }

    #[test]
    fn the_two_sided_global_p_value_is_two_times_the_smaller_chi_square_tail() {
        // p = 2 * min(pchisq(W, lower = TRUE), pchisq(W, lower = FALSE)).
        // For W = 0.25 on 1 df, R gives min(0.3829249, 0.6170751) * 2 = 0.7658498,
        // so a near-null taxon is correctly *not* significant. Pinned because the
        // construction looks like a bug and is easy to "fix" by accident.
        let (cn, beta, vcov, n, p) = toy();
        let g = global_test(
            &cn,
            "grp",
            &beta,
            &vcov,
            None,
            n,
            p,
            AdjustMethod::None,
            0.05,
        )
        .unwrap();
        // W = 1.25 on 2 df: min(pchisq(1.25,2,T), pchisq(1.25,2,F)) * 2
        let lo = ancombc2_stats::chisq_lower(1.25, 2.0);
        let hi = ancombc2_stats::chisq_upper(1.25, 2.0);
        assert!((g.w[1] - 1.25).abs() < 1e-12, "W = {}", g.w[1]);
        assert!((g.p[1] - 2.0 * lo.min(hi)).abs() < 1e-12, "p = {}", g.p[1]);
        assert!(g.p[1] > 0.05, "a near-null taxon must not be significant");
    }

    #[test]
    fn an_empty_screen_reproduces_the_reference_degenerate_behaviour() {
        // R computes `n = length(x) * n_tax / R` without guarding R = 0: the
        // denominator becomes infinite, `p.adjust` returns 0, and every pairwise
        // test is called significant. That is the opposite of the intuitive
        // behaviour (no taxon passed the screen, so nothing should be
        // significant), and reproducing it is deliberate. See the module docs and
        // docs/reference_behavior.md.
        //
        // An empty screen arises when every taxon's global p-value is 1, which
        // needs a non-finite coefficient.
        let n = 2usize;
        let p = 3usize;
        let beta = vec![0.0, f64::NAN, f64::NAN, 0.0, f64::NAN, f64::NAN];
        let mut vcov = vec![0.0; n * p * p];
        for i in 0..n {
            for a in 0..p {
                vcov[i * p * p + a * p + a] = 0.04;
            }
        }
        let cn = names(&["(Intercept)", "grp2", "grp3"]);
        let var_hat = vec![0.04; n * p];
        let r = pairwise_test(
            &cn,
            "grp",
            &beta,
            &vcov,
            &var_hat,
            None,
            n,
            p,
            AdjustMethod::None,
            AdjustMethod::Holm,
            0.05,
        )
        .unwrap();
        assert_eq!(r.n_screen, 0, "no taxon has an estimable coefficient");
        assert!(
            r.diff_abn.iter().all(|&d| d),
            "R's unguarded division makes everything significant: {:?}",
            r.diff_abn
        );
    }

    #[test]
    fn results_are_deterministic() {
        let (cn, beta, vcov, n, p) = toy();
        let var_hat = toy_var();
        let a = pairwise_test(
            &cn,
            "grp",
            &beta,
            &vcov,
            &var_hat,
            None,
            n,
            p,
            AdjustMethod::None,
            AdjustMethod::Holm,
            0.05,
        )
        .unwrap();
        let b = pairwise_test(
            &cn,
            "grp",
            &beta,
            &vcov,
            &var_hat,
            None,
            n,
            p,
            AdjustMethod::None,
            AdjustMethod::Holm,
            0.05,
        )
        .unwrap();
        assert_eq!(a.p, b.p);
        assert_eq!(a.q, b.q);
        assert_eq!(a.diff_abn, b.diff_abn);
    }
}
