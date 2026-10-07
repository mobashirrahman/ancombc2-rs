//! Pseudo-count sensitivity analysis.
//!
//! ANCOM-BC2's answer to "does the result depend on the pseudo-count added to
//! the zeros?" has two modes, and they are structurally different.
//!
//! # Conservative
//!
//! The whole algorithm is re-run for `pseudo in {0.1, 0.5, 1}` on top of the
//! main run, and the agreement is measured on the **adjusted** p-values:
//!
//! ```text
//! s = (1/4) * #{ c : q_c > alpha }
//! passed_ss = (s == 0) or (s == 1)
//! ```
//!
//! The four runs share only the preprocessed design, so this is embarrassingly
//! parallel at the top level of the parallel plan.
//!
//! # Non-conservative
//!
//! The sampling fractions are estimated **once**, on the main run's data. For
//! each of the 50 pseudo-counts `0.01 … 0.50` the zeros of the bias-corrected
//! table are filled in, the table is re-logged and re-centred, the sampling
//! fraction subtracted, and the model is refitted by ordinary least squares. The
//! agreement is measured on **unadjusted** p-values:
//!
//! ```text
//! s = (1/50) * #{ c : p_c > alpha }
//! passed_ss = (s == 0 and p_0 <= alpha) or (s == 1 and p_0 > alpha)
//! ```
//!
//! The reference reports higher power *and* a higher false-positive rate for
//! this mode, which is why it is not the default.
//!
//! # Two things that are easy to get wrong
//!
//! * The refit in the non-conservative mode is a **plain `lm` per taxon**, with
//!   t-based p-values on `df.residual` — not the ANCOM-BC2 `W` statistic and not
//!   the sandwich variance. Reproducing that faithfully means reproducing `lm`.
//! * A coefficient that is not estimable gets `p = NA` from `lm`, and the
//!   reference then replaces it with 1.

use crate::config::{AdjustMethod, AncombcConfig};
use crate::error::Result;
use crate::workspace::RMatrix;

/// A sensitivity score table: one column per tested quantity.
#[derive(Debug, Clone)]
pub struct SensitivityScores {
    /// `n_taxa x n_col` row-major scores in `[0, 1]`.
    pub scores: Vec<f64>,
    /// `n_taxa x n_col`; `true` where the taxon agreed across the grid.
    pub passed: Vec<bool>,
    /// Column names, in the order of the tested quantities.
    pub colnames: Vec<String>,
    /// The pseudo-count grid that was used.
    pub pseudo: Vec<f64>,
}

impl SensitivityScores {
    pub fn n_taxa(&self) -> usize {
        if self.colnames.is_empty() {
            0
        } else {
            self.scores.len() / self.colnames.len()
        }
    }
}

/// Conservative mode: the score of one column across the runs.
///
/// `q_by_run[r]` is the adjusted p-value column of run `r`; a taxon agrees when
/// the score is exactly 0 or exactly 1, i.e. every run agreed with every other.
pub fn conservative_scores(q_by_run: &[Vec<f64>], n_taxa: usize, alpha: f64) -> Vec<f64> {
    let n_runs = q_by_run.len();
    let mut out = vec![0.0; n_taxa];
    for i in 0..n_taxa {
        let mut above = 0usize;
        for r in 0..n_runs {
            let v = q_by_run[r][i];
            if v.is_nan() || v > alpha {
                above += 1;
            }
        }
        out[i] = above as f64 / n_runs as f64;
    }
    out
}

/// `passed_ss = (s == 0) | (s == 1)`.
pub fn conservative_passed(scores: &[f64]) -> Vec<bool> {
    scores.iter().map(|&s| s == 0.0 || s == 1.0).collect()
}

/// Non-conservative mode: the score of one column across the 50 refits.
///
/// `p_by_refit[r]` is the unadjusted p-value column of refit `r`.
pub fn nonconservative_scores(p_by_refit: &[Vec<f64>], n_taxa: usize, alpha: f64) -> Vec<f64> {
    let n = p_by_refit.len();
    let mut out = vec![0.0; n_taxa];
    for i in 0..n_taxa {
        let mut above = 0usize;
        for r in 0..n {
            let v = p_by_refit[r][i];
            if v.is_nan() || v > alpha {
                above += 1;
            }
        }
        out[i] = above as f64 / n as f64;
    }
    out
}

/// `passed_ss = (s == 0 and p_0 <= alpha) or (s == 1 and p_0 > alpha)`.
///
/// The main run's p-value enters the decision, unlike the conservative rule.
pub fn nonconservative_passed(scores: &[f64], p_main: &[f64], alpha: f64) -> Vec<bool> {
    (0..scores.len())
        .map(|i| {
            let p = if p_main[i].is_nan() { 1.0 } else { p_main[i] };
            (scores[i] == 0.0 && p <= alpha) || (scores[i] == 1.0 && p > alpha)
        })
        .collect()
}

/// `.ancombc2_sens_fit`: build the response matrix for one pseudo-count.
///
/// The reference takes the bias-corrected table's *source* counts (`O2` with the
/// run's pseudo already added), replaces the zeros with `pseudo`, logs, centres
/// per taxon, and subtracts the sampling fractions. Note the sign: the response
/// is `samp_frac - y`, because the reference writes `t(t(y) - samp_frac)`, which
/// in R's column-major `t()` idiom is the *transpose* of the difference, i.e.
/// `samp_frac` minus `y`.
/// # The reduction is a parameter
///
/// The centring is `rowMeans`, and R accumulates that in `long double`. See
/// [`crate::reduce`]. So the reduction is supplied rather than defaulted: a caller
/// that reaches a result through here has to have said which accumulator produced
/// the mean, at the call site, where it is visible.
pub fn sens_response(
    counts: &crate::preprocess::CountMatrix,
    samp_frac: &[f64],
    pseudo: f64,
    red: &dyn crate::reduce::Reductions,
) -> RMatrix {
    let mut y = counts.log_center_replacing_zeros_with(pseudo, red);
    y.sub_cols_in_place(samp_frac);
    y
}

/// [`sens_response`] on a sub-table of `counts`, without materialising that
/// sub-table. See [`crate::preprocess::CountMatrix::log_center_replacing_zeros_sub`]
/// for why both axes are selected here.
/// [`sens_response`] with the reduction supplied. See that function for why it is a
/// parameter rather than a default.
pub fn sens_response_sub(
    counts: &crate::preprocess::CountMatrix,
    rows: &[usize],
    cols: &[usize],
    samp_frac: &[f64],
    pseudo: f64,
    red: &dyn crate::reduce::Reductions,
) -> RMatrix {
    let means = sens_response_means(counts, rows, cols, pseudo, red);
    sens_response_sub_with_means(counts, rows, cols, samp_frac, pseudo, &means)
}

/// The per-row centring means for one pseudo-count, and only that.
///
/// Split out for the same reason as [`crate::preprocess::CountMatrix::log_row_means`]:
/// the reduction may be R, which cannot be called from a worker thread, and the
/// non-conservative refits run in parallel. So the means are reduced once on the calling
/// thread -- `grid.len() * rows.len()` doubles -- and
/// [`sens_response_sub_with_means`] redoes only the arithmetic.
pub fn sens_response_means(
    counts: &crate::preprocess::CountMatrix,
    rows: &[usize],
    cols: &[usize],
    pseudo: f64,
    red: &dyn crate::reduce::Reductions,
) -> Vec<f64> {
    let y = counts.log_replacing_zeros_sub(rows, cols, pseudo);
    red.row_means_na_rm(&y)
}

/// [`sens_response_sub`] with the centring means supplied instead of a reduction.
///
/// Same values as [`sens_response_sub`], cell for cell; a test asserts it, so the split
/// cannot quietly become a difference.
pub fn sens_response_sub_with_means(
    counts: &crate::preprocess::CountMatrix,
    rows: &[usize],
    cols: &[usize],
    samp_frac: &[f64],
    pseudo: f64,
    means: &[f64],
) -> RMatrix {
    let mut y = counts.log_replacing_zeros_sub(rows, cols, pseudo);
    y.sub_rows_in_place(means);
    y.sub_cols_in_place(samp_frac);
    y
}

/// Run the non-conservative refits in parallel and reduce them to scores.
///
/// `fit_one` receives the pseudo-count and must return the unadjusted p-values
/// of every (taxon, fix_eff) pair as an `n_taxa * n_col` row-major block. The
/// reference's `.ancombc2_sens_fit` fits one `lm` per taxon and returns an
/// `n_tax x n_fix_eff` frame, and `ss_tab_fun` then reduces one *column* at a
/// time, so the sensitivity score is per coefficient and not per taxon. Taking
/// a single column would silently report the intercept's score for every
/// coefficient.
///
/// Splitting the refit out keeps this function testable without a design matrix,
/// and mirrors the reference's structure where `.ancombc2_sens_fit` owns the
/// transform and `.ancombc2_sens_p` owns the model.
pub fn nonconservative_run<F>(
    grid: &[f64],
    n_taxa: usize,
    n_col: usize,
    alpha: f64,
    fit_one: F,
) -> Result<SensitivityScores>
where
    F: Fn(usize, f64) -> Vec<f64> + Sync,
{
    let total = n_taxa * n_col;
    // Level 1 of the nesting order: the pseudo-count runs, outermost.
    //
    // The plan's order is pseudo-count runs outermost, taxa innermost, and the
    // budget exists so the pool is never split twice. But "outermost" is a claim
    // about *which level gets the pool*, not a statement that the outer level can
    // fill it: this grid is 50 refits on `bm6` and 3 on `bm5`, and on a 16-thread
    // host 3 items occupy three threads.
    //
    // So the level takes the pool only when it has at least as many items as the
    // pool has threads. Fifty refits do, three do not, and in the three case the
    // refits run one at a time with the whole pool inside each -- which on `bm5`
    // is the difference between three busy cores and sixteen.
    let mut lvl = crate::parallel::NestingBudget::level_for_items("pseudo-count runs", grid.len());
    let by_refit: Vec<Vec<f64>> = crate::parallel::map_par(&mut lvl, grid, |k, &pc| {
        let p = fit_one(k, pc);
        debug_assert_eq!(p.len(), total);
        p
    });
    // `ss_tab_fun`: rowMeans(p_pseudo > alpha) over the refits, column by
    // column. A fixed order over the refits keeps the reduction deterministic.
    let mut scores = vec![0.0; total];
    for c in 0..n_col {
        let col: Vec<Vec<f64>> = by_refit
            .iter()
            .map(|run| (0..n_taxa).map(|i| run[i * n_col + c]).collect())
            .collect();
        let s = nonconservative_scores(&col, n_taxa, alpha);
        for i in 0..n_taxa {
            scores[i * n_col + c] = s[i];
        }
    }
    Ok(SensitivityScores {
        scores,
        passed: Vec::new(),
        colnames: Vec::new(),
        pseudo: grid.to_vec(),
    })
}

/// Conservative mode, given the adjusted p-value columns of all four runs.
pub fn conservative_run(
    q_by_run: &[Vec<f64>],
    n_taxa: usize,
    colnames: Vec<String>,
    alpha: f64,
    grid: &[f64],
) -> SensitivityScores {
    let n_col = colnames.len();
    let mut scores = vec![0.0; n_taxa * n_col];
    let mut passed = vec![false; n_taxa * n_col];
    for c in 0..n_col {
        let col: Vec<Vec<f64>> = q_by_run
            .iter()
            .map(|run| {
                (0..n_taxa)
                    .map(|i| run[i * q_by_run[0].len() / n_taxa + c])
                    .collect()
            })
            .collect();
        let s = conservative_scores(&col, n_taxa, alpha);
        let p = conservative_passed(&s);
        for i in 0..n_taxa {
            scores[i * n_col + c] = s[i];
            passed[i * n_col + c] = p[i];
        }
    }
    SensitivityScores {
        scores,
        passed,
        colnames,
        pseudo: grid.to_vec(),
    }
}

/// The pseudo-count grids, from the configuration.
pub fn grid_for(cfg: &AncombcConfig) -> Vec<f64> {
    if cfg.conservative {
        cfg.conservative_pseudo_grid()
    } else {
        cfg.nonconservative_pseudo_grid()
    }
}

/// The adjustment method the sensitivity table is keyed on, for documentation.
pub fn default_method() -> AdjustMethod {
    AdjustMethod::None
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::preprocess::CountMatrix;
    use crate::reduce::F64Reductions;

    #[test]
    fn conservative_grid_matches_the_reference() {
        let cfg = AncombcConfig::default();
        assert_eq!(grid_for(&cfg), vec![0.0, 0.1, 0.5, 1.0]);
    }

    #[test]
    fn nonconservative_grid_matches_the_reference() {
        let cfg = AncombcConfig {
            conservative: false,
            ..Default::default()
        };
        let g = grid_for(&cfg);
        assert_eq!(g.len(), 50);
        assert!((g[0] - 0.01).abs() < 1e-15);
        assert!((g[49] - 0.50).abs() < 1e-15);
    }

    #[test]
    fn conservative_score_counts_runs_above_alpha() {
        // 2 taxa x 4 runs
        let runs = vec![
            vec![0.01, 0.9], // taxon 0 below, taxon 1 above
            vec![0.01, 0.01],
            vec![0.02, 0.01],
            vec![0.01, 0.01],
        ];
        let s = conservative_scores(&runs, 2, 0.05);
        assert!((s[0] - 0.0).abs() < 1e-15, "{}", s[0]);
        assert!((s[1] - 0.25).abs() < 1e-15, "{}", s[1]);
        let p = conservative_passed(&s);
        assert!(p[0], "s = 0 passes");
        assert!(!p[1], "s = 0.25 is neither 0 nor 1");
    }

    #[test]
    fn conservative_treats_exactly_alpha_as_agreement() {
        // The reference tests `q > alpha`, so q == alpha counts as agreement.
        let runs = vec![vec![0.05], vec![0.05]];
        let s = conservative_scores(&runs, 1, 0.05);
        assert!((s[0] - 0.0).abs() < 1e-15, "q == alpha is not 'above'");
        let runs = vec![vec![0.0500001], vec![0.05]];
        let s = conservative_scores(&runs, 1, 0.05);
        assert!((s[0] - 0.5).abs() < 1e-15);
    }

    #[test]
    fn nonconservative_passed_uses_the_main_run_p_value() {
        // s = 0 with a significant main run passes
        let p = nonconservative_passed(&[0.0], &[0.01], 0.05);
        assert!(p[0]);
        // s = 0 with a non-significant main run fails
        let p = nonconservative_passed(&[0.0], &[0.5], 0.05);
        assert!(!p[0]);
        // s = 1 with a non-significant main run passes
        let p = nonconservative_passed(&[1.0], &[0.5], 0.05);
        assert!(p[0]);
        // s = 1 with a significant main run fails
        let p = nonconservative_passed(&[1.0], &[0.01], 0.05);
        assert!(!p[0]);
        // an intermediate score never passes
        let p = nonconservative_passed(&[0.4], &[0.01], 0.05);
        assert!(!p[0]);
    }

    #[test]
    fn nonconservative_treats_na_main_p_as_one() {
        // R: p_main[is.na(p_main)] <- 1
        let p = nonconservative_passed(&[0.0], &[f64::NAN], 0.05);
        assert!(!p[0], "an NA main p behaves as 1, which is > alpha");
    }

    #[test]
    fn sens_response_subtracts_the_sampling_fractions() {
        // 1 taxon, 2 samples, counts 1 and 4
        let c = CountMatrix::new(1, 2, vec![1.0, 4.0]).unwrap();
        let out = sens_response(&c, &[0.5, 0.5], 0.01, &F64Reductions);
        // y = (log 1, log 4) centred -> (-log 2, log 2); then minus 0.5
        let l2 = 2.0f64.ln();
        assert!(
            (out.get(0, 0) - (-l2 - 0.5)).abs() < 1e-14,
            "{}",
            out.get(0, 0)
        );
        assert!(
            (out.get(0, 1) - (l2 - 0.5)).abs() < 1e-14,
            "{}",
            out.get(0, 1)
        );
    }

    #[test]
    fn sens_response_fills_zeros_only_for_positive_pseudo() {
        let c = CountMatrix::new(1, 2, vec![0.0, 4.0]).unwrap();
        // pseudo = 0 leaves the zero missing
        let out0 = sens_response(&c, &[0.0, 0.0], 0.0, &F64Reductions);
        assert!(out0.get(0, 0).is_nan());
        // a positive pseudo fills it, and the centring then uses both entries
        let out1 = sens_response(&c, &[0.0, 0.0], 0.01, &F64Reductions);
        assert!(out1.get(0, 0).is_finite());
        // and the two columns are negatives of each other after centring
        assert!((out1.get(0, 0) + out1.get(0, 1)).abs() < 1e-15);
    }

    #[test]
    fn nonconservative_run_reduces_over_the_whole_grid() {
        let n_taxa = 3usize;
        let grid: Vec<f64> = (1..=50).map(|i| i as f64 / 100.0).collect();
        // One coefficient. Taxon 0 is significant under every pseudo-count,
        // taxon 1 under none, taxon 2 under half.
        let fit = move |_k: usize, pc: f64| -> Vec<f64> {
            vec![0.001, 0.9, if pc <= 0.25 { 0.001 } else { 0.9 }]
        };
        let r = nonconservative_run(&grid, n_taxa, 1, 0.05, fit).unwrap();
        assert_eq!(r.pseudo.len(), 50);
        assert!((r.scores[0] - 0.0).abs() < 1e-15, "{}", r.scores[0]);
        assert!((r.scores[1] - 1.0).abs() < 1e-15, "{}", r.scores[1]);
        // 25 of the 50 refits are significant
        assert!((r.scores[2] - 0.5).abs() < 1e-15, "{}", r.scores[2]);
    }

    #[test]
    fn parallel_refits_agree_with_sequential_ones() {
        let n_taxa = 4usize;
        let grid: Vec<f64> = (1..=50).map(|i| i as f64 / 100.0).collect();
        let f_par = move |_k: usize, pc: f64| -> Vec<f64> {
            (0..n_taxa).map(|i| pc * (i as f64 + 1.0)).collect()
        };
        let par = nonconservative_run(&grid, n_taxa, 1, 0.05, f_par).unwrap();
        let seq: Vec<Vec<f64>> = grid
            .iter()
            .map(|&pc| (0..n_taxa).map(|i| pc * (i as f64 + 1.0)).collect())
            .collect();
        let expect = nonconservative_scores(&seq, n_taxa, 0.05);
        assert_eq!(par.scores, expect, "thread count must not matter");
    }

    #[test]
    fn conservative_run_assembles_a_matrix_of_columns() {
        // 2 taxa, 2 columns, 4 runs
        let n_taxa = 2usize;
        let n_col = 2usize;
        let runs: Vec<Vec<f64>> = (0..4)
            .map(|r| {
                (0..n_taxa * n_col)
                    .map(|k| {
                        if r % 2 == 0 {
                            0.01
                        } else {
                            if k < 2 {
                                0.9
                            } else {
                                0.01
                            }
                        }
                    })
                    .collect()
            })
            .collect();
        let s = conservative_run(
            &runs,
            n_taxa,
            vec!["q_a".into(), "q_b".into()],
            0.05,
            &[0.0, 0.1, 0.5, 1.0],
        );
        assert_eq!(s.colnames, vec!["q_a", "q_b"]);
        // column a: two of the four runs are above alpha
        assert!((s.scores[0] - 0.5).abs() < 1e-15, "{}", s.scores[0]);
        assert!((s.scores[2] - 0.0).abs() < 1e-15, "{}", s.scores[2]);
        assert!(!s.passed[0]);
        assert!(s.passed[2]);
    }

    #[test]
    fn scores_are_deterministic() {
        let runs = vec![vec![0.01, 0.9], vec![0.02, 0.8]];
        let a = conservative_scores(&runs, 2, 0.05);
        let b = conservative_scores(&runs, 2, 0.05);
        assert_eq!(a, b);
    }
}
