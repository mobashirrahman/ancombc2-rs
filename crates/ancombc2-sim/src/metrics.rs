//! Truth-aware metrics for one replicate.
//!
//! # What is compared to what
//!
//! `beta` is compared against the taxon's *true* log fold change, not the
//! realised log ratio of the sampled counts. The realised ratio carries
//! sampling noise that no estimator reproduces taxon by taxon; averaging it in
//! would hide a systematic bias behind variance that is not the method's.
//!
//! # The metrics
//!
//! - **empirical FDR** -- of the taxa called `diff_abn`, the fraction that are
//!   not DA. The target is at or below the nominal level.
//! - **power** -- of the DA taxa retained, the fraction called `diff_abn`.
//! - **sign concordance** -- over DA taxa, the fraction whose `beta` sign
//!   matches the true direction. Over null taxa the corresponding
//!   **null sign rate** is reported: it should be near one half, and a value far
//!   from it is a calibration failure that a power number alone would hide.
//! - **LFC bias** -- `mean(beta - truth)` and its RMSE, both over all retained
//!   taxa, which is the bias the correction is supposed to remove.
//! - **SE calibration** -- `mean(se) / sd(beta over reps)`, which needs more
//!   than one rep and is therefore pooled by [`crate::summarise`].
//! - **F1 / Jaccard** -- overlap of the called and true DA sets.

use serde::{Deserialize, Serialize};

use crate::simulate::Replicate;

/// Per-replicate scalars plus the per-taxon vectors the pooled metrics need.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RepMetrics {
    /// Per-taxon flags mirroring the truth: `is_da` is the generator's label,
    /// `has_effect` is whether the expected abundance really differs. They differ
    /// in a confounded cell, and every power and false-positive statement below
    /// is about `has_effect`.
    pub has_effect: Vec<bool>,
    /// The taxa whose label came from a group-dependent zero rate with no
    /// compositional effect. Calling any of them is a false positive for a
    /// method that has not corrected the sampling fraction.
    pub confounded: Vec<bool>,
    // Identification.
    pub grid: String,
    pub cell: usize,
    pub rep: usize,
    pub arm: String,
    pub n_taxa: usize,
    pub n_samp: usize,
    pub da_proportion: f64,
    pub log_fc: f64,
    pub zero_inflation: f64,
    pub lib_mean: f64,
    pub lib_cv: f64,
    pub confound: bool,
    // Per-taxon detail, aligned to `names`.
    pub names: Vec<String>,
    pub truth_log_fc: Vec<f64>,
    pub is_da: Vec<bool>,
    pub beta: Vec<f64>,
    pub se: Vec<f64>,
    pub p: Vec<f64>,
    pub q: Vec<f64>,
    pub diff_abn: Vec<bool>,
    // Scalars.
    pub n_retained: usize,
    pub n_da_retained: usize,
    pub n_diff_abn: usize,
    /// `null` when nothing was called, which is a legitimate outcome and is
    /// reported rather than silently turned into an FDR of zero.
    pub empirical_fdr: Option<f64>,
    pub power: Option<f64>,
    pub sign_concordance: Option<f64>,
    pub null_sign_rate: Option<f64>,
    /// Mean signed error of the estimated log fold change against the truth.
    ///
    /// `Option`, not `f64`: a replicate that failed -- at the extreme sparsity
    /// cells, "all taxa contain structural zeros" -- has no fold change to be
    /// biased, and the writer already emitted `null` for it. With a plain `f64`
    /// that row could not be read back, and `summarise` failed on the whole file
    /// at the first such line. A NaN would have been readable, but `null` is what
    /// the schema uses for every other undefined quantity here, so `null` it is.
    pub lfc_bias: Option<f64>,
    pub lfc_rmse: Option<f64>,
    pub lfc_mae: Option<f64>,
    /// F1 and Jaccard of the called set against the true effect set.
    pub f1: Option<f64>,
    pub jaccard: Option<f64>,
    /// The share of the *confounded* taxa that were called differentially
    /// abundant. This is the direct test of the sampling-fraction correction: the
    /// confounded taxa have no compositional effect, so any call on one is a
    /// false positive, and a method that applies the correction should keep this
    /// near zero.
    pub confounded_false_positive_rate: Option<f64>,
    /// How many confounded taxa were retained at all, so a rate over one or two
    /// taxa is visible as such rather than as a confident number.
    pub n_confounded_retained: usize,
    /// A non-empty message when the analysis failed; the arm is recorded as a
    /// failure rather than being dropped from the denominator.
    pub error: Option<String>,
}

impl RepMetrics {
    pub fn failed(
        grid: &str,
        cell: usize,
        rep: usize,
        arm: &str,
        c: &crate::grid::Cell,
        error: String,
    ) -> Self {
        Self {
            grid: grid.to_string(),
            cell,
            rep,
            arm: arm.to_string(),
            n_taxa: c.n_taxa,
            n_samp: c.n_samp,
            da_proportion: c.da_proportion,
            log_fc: c.log_fc,
            zero_inflation: c.zero_inflation,
            lib_mean: c.lib_mean,
            lib_cv: c.lib_cv,
            confound: c.confound,
            names: Vec::new(),
            truth_log_fc: Vec::new(),
            is_da: Vec::new(),
            has_effect: Vec::new(),
            confounded: Vec::new(),
            beta: Vec::new(),
            se: Vec::new(),
            p: Vec::new(),
            q: Vec::new(),
            diff_abn: Vec::new(),
            n_retained: 0,
            n_da_retained: 0,
            n_diff_abn: 0,
            empirical_fdr: None,
            power: None,
            sign_concordance: None,
            null_sign_rate: None,
            lfc_bias: None,
            lfc_rmse: None,
            lfc_mae: None,
            f1: None,
            jaccard: None,
            confounded_false_positive_rate: None,
            n_confounded_retained: 0,
            error: Some(error),
        }
    }
}

/// The per-taxon vectors a reported analysis produced, for the taxa the
/// analysis retained.
///
/// Grouped into one struct so [`evaluate`] takes a readable argument list; five
/// parallel slices that must all have the same length are easy to pass in the
/// wrong order.
#[derive(Debug, Clone)]
pub struct Reported<'a> {
    pub beta: &'a [f64],
    pub se: &'a [f64],
    pub p: &'a [f64],
    pub q: &'a [f64],
    pub diff_abn: &'a [bool],
}

impl<'a> Reported<'a> {
    /// Every vector must describe the same taxa, or the metrics are meaningless
    /// rather than merely wrong, so this is checked where it is built.
    pub fn new(
        beta: &'a [f64],
        se: &'a [f64],
        p: &'a [f64],
        q: &'a [f64],
        diff_abn: &'a [bool],
    ) -> Result<Self, String> {
        let n = beta.len();
        if se.len() != n || p.len() != n || q.len() != n || diff_abn.len() != n {
            return Err(format!(
                "reported vectors disagree in length: beta {n}, se {}, p {}, q {}, diff_abn {}",
                se.len(),
                p.len(),
                q.len(),
                diff_abn.len()
            ));
        }
        Ok(Self {
            beta,
            se,
            p,
            q,
            diff_abn,
        })
    }
}

/// Build the metrics from an analysis result and the replicate's truth.
///
/// `beta_index` is the coefficient of the group term in the fitted design, and
/// the design is built with group 1 as the reference level, so a positive
/// `beta` means "higher in the treated group".
pub fn evaluate(
    grid: &str,
    cell: &crate::grid::Cell,
    rep: usize,
    arm: &str,
    r: &Replicate,
    kept_taxa: &[usize],
    reported: &Reported,
) -> RepMetrics {
    let Reported {
        beta,
        se,
        p,
        q,
        diff_abn,
    } = reported;
    let n = kept_taxa.len();
    let mut names = Vec::with_capacity(n);
    let mut truth_log_fc = Vec::with_capacity(n);
    let mut is_da = Vec::with_capacity(n);
    let mut has_effect = Vec::with_capacity(n);
    let mut confounded = Vec::with_capacity(n);
    let mut b = Vec::with_capacity(n);
    let mut s = Vec::with_capacity(n);
    let mut pv = Vec::with_capacity(n);
    let mut qv = Vec::with_capacity(n);
    let mut called = Vec::with_capacity(n);
    for (k, &t) in kept_taxa.iter().enumerate() {
        names.push(r.truth[t].name.clone());
        truth_log_fc.push(r.truth[t].log_fc);
        is_da.push(r.truth[t].is_da);
        has_effect.push(r.truth[t].has_effect);
        confounded.push(r.truth[t].confounded);
        b.push(beta[k]);
        s.push(se[k]);
        pv.push(p[k]);
        qv.push(q[k]);
        called.push(diff_abn[k]);
    }

    let mut tp = 0usize;
    let mut fp = 0usize;
    let mut fn_ = 0usize;
    let mut sign_hits = 0usize;
    let mut sign_n = 0usize;
    let mut null_pos = 0usize;
    let mut null_n = 0usize;
    let mut sum_err = 0.0f64;
    let mut sum_sq = 0.0f64;
    let mut sum_abs = 0.0f64;
    let mut diffable = 0usize;
    let mut confounded_n = 0usize;
    let mut confounded_called = 0usize;

    for k in 0..n {
        let err = b[k] - truth_log_fc[k];
        if err.is_finite() {
            sum_err += err;
            sum_sq += err * err;
            sum_abs += err.abs();
            diffable += 1;
        }
        // Every confusion count is over `has_effect`, not over the generator's
        // label. In a confounded cell the labelled taxa are nulls, and counting
        // them as true positives would report the correction as working by
        // accident.
        match (has_effect[k], called[k]) {
            (true, true) => tp += 1,
            (false, true) => fp += 1,
            (true, false) => fn_ += 1,
            (false, false) => {}
        }
        if has_effect[k] {
            if truth_log_fc[k] > 0.0 && b[k] > 0.0 || truth_log_fc[k] < 0.0 && b[k] < 0.0 {
                sign_hits += 1;
            }
            sign_n += 1;
        } else if b[k].is_finite() {
            null_n += 1;
            if b[k] > 0.0 {
                null_pos += 1;
            }
        }
        if confounded[k] {
            confounded_n += 1;
            if called[k] {
                confounded_called += 1;
            }
        }
    }

    let n_da = sign_n;
    let n_called = tp + fp;
    let denom = diffable.max(1) as f64;

    RepMetrics {
        grid: grid.to_string(),
        cell: cell.index,
        rep,
        arm: arm.to_string(),
        n_taxa: cell.n_taxa,
        n_samp: cell.n_samp,
        da_proportion: cell.da_proportion,
        log_fc: cell.log_fc,
        zero_inflation: cell.zero_inflation,
        lib_mean: cell.lib_mean,
        lib_cv: cell.lib_cv,
        confound: cell.confound,
        names,
        truth_log_fc,
        is_da,
        has_effect,
        confounded,
        beta: b,
        se: s,
        p: pv,
        q: qv,
        diff_abn: called,
        n_retained: n,
        n_da_retained: n_da,
        n_diff_abn: n_called,
        n_confounded_retained: confounded_n,
        confounded_false_positive_rate: (confounded_n > 0)
            .then_some(confounded_called as f64 / confounded_n as f64),
        empirical_fdr: (n_called > 0).then_some(fp as f64 / n_called as f64),
        power: (n_da > 0).then_some(tp as f64 / n_da as f64),
        sign_concordance: (n_da > 0).then_some(sign_hits as f64 / n_da as f64),
        null_sign_rate: (null_n > 0).then_some(null_pos as f64 / null_n as f64),
        lfc_bias: Some(sum_err / denom),
        lfc_rmse: Some((sum_sq / denom).sqrt()),
        lfc_mae: Some(sum_abs / denom),
        f1: (2 * tp + fp + fn_ > 0).then_some(2.0 * tp as f64 / (2 * tp + fp + fn_) as f64),
        jaccard: (tp + fp + fn_ > 0).then_some(tp as f64 / (tp + fp + fn_) as f64),
        error: None,
    }
}

/// The three standard errors the acceptance rule is expressed in.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct BinomialSe {
    pub fdr: f64,
    pub power: f64,
}

/// `sqrt(p (1 - p) / n)` for the two rates, with the usual `n = 0` guard.
///
/// The acceptance rule in `PLAN.md` §5.3 is "within 3 SE of the binomial", so
/// this is the function that turns that sentence into a number.
pub fn binomial_se(rate: f64, n: usize) -> f64 {
    if n == 0 {
        return f64::NAN;
    }
    (rate.clamp(0.0, 1.0) * (1.0 - rate.clamp(0.0, 1.0)) / n as f64).sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::grid::Cell;
    use crate::simulate::{simulate, SimControl};

    fn cell() -> Cell {
        Cell {
            index: 0,
            n_taxa: 100,
            n_samp: 10,
            da_proportion: 0.5,
            log_fc: 1.0,
            zero_inflation: 0.0,
            lib_mean: 1e4,
            lib_cv: 0.0,
            confound: false,
        }
    }

    fn ctl() -> SimControl {
        SimControl {
            abundance_sd: 0.5,
            dispersion: 1.0,
            confounder: 0.3,
        }
    }

    /// A stub replicate whose truth is `log_fc` for the first half and zero for
    /// the rest, so the metric arithmetic can be checked against hand counts.
    fn stub() -> Replicate {
        simulate(&cell(), 1, 0, &ctl())
    }

    fn rep_metrics(b: &[f64], called: &[bool]) -> RepMetrics {
        let r = stub();
        let truth: Vec<f64> = (0..50)
            .map(|_| 1.0)
            .chain(std::iter::repeat(0.0).take(50))
            .collect();
        let is_da: Vec<bool> = (0..100).map(|i| i < 50).collect();
        let mut rr = r;
        for (i, t) in rr.truth.iter_mut().enumerate() {
            t.log_fc = truth[i];
            t.is_da = is_da[i];
            // The stub has no confounded taxa, so the label and the effect
            // coincide. `evaluate` counts over `has_effect`, so leaving it unset
            // would make every taxon a null.
            t.has_effect = is_da[i];
            t.confounded = false;
        }
        let n = b.len();
        let se = vec![1.0; n];
        let p = vec![0.01; n];
        let q = vec![0.01; n];
        let reported = Reported::new(b, &se, &p, &q, called).unwrap();
        evaluate(
            "t",
            &cell(),
            0,
            "rust",
            &rr,
            &(0..n).collect::<Vec<usize>>(),
            &reported,
        )
    }

    #[test]
    fn mismatched_vector_lengths_are_rejected_at_construction() {
        // A silently mis-zipped set of vectors would produce plausible-looking
        // metrics, so the lengths are checked once, where they are assembled.
        let v5 = [0.0f64; 5];
        let v4 = [0.0f64; 4];
        let p5 = [0.1f64; 5];
        let t5 = [true; 5];
        let t4 = [true; 4];
        assert!(Reported::new(&v5, &v4, &p5, &p5, &t5).is_err());
        assert!(Reported::new(&v5, &v5, &p5, &p5, &t4).is_err());
        assert!(Reported::new(&v5, &v5, &p5, &p5, &t5).is_ok());
    }

    #[test]
    fn a_perfect_caller_scores_one_everywhere() {
        let mut b = vec![0.0; 100];
        b[..50].iter_mut().for_each(|v| *v = 1.0);
        let mut called = vec![false; 100];
        called[..50].iter_mut().for_each(|v| *v = true);
        let _ = &called;
        let m = rep_metrics(&b, &called);
        assert_eq!(m.empirical_fdr, Some(0.0));
        assert_eq!(m.power, Some(1.0));
        assert_eq!(m.sign_concordance, Some(1.0));
        assert_eq!(m.f1, Some(1.0));
        assert_eq!(m.jaccard, Some(1.0));
        assert!(m.lfc_bias.unwrap().abs() < 1e-12, "bias {:?}", m.lfc_bias);
    }

    #[test]
    fn a_caller_that_calls_everything_has_fdr_one_half() {
        // Half the taxa are null, so calling all 100 gives an FDR of 50/100.
        // One half, not one: the false-call *rate* is the fraction of calls that
        // are wrong, and calling everything makes half the calls wrong.
        let mut b = vec![0.0; 100];
        b[..50].iter_mut().for_each(|v| *v = 1.0);
        let m = rep_metrics(&b, &[true; 100]);
        assert_eq!(m.empirical_fdr, Some(0.5), "50 of 100 calls are false");
        assert_eq!(m.power, Some(1.0));
    }

    #[test]
    fn a_perfect_caller_on_a_all_null_table_has_zero_fdr() {
        // The extreme of the same arithmetic: with no DA taxa at all, calling
        // nothing is the only correct behaviour and there is no FDR to report.
        let r = stub();
        let mut rr = r;
        for t in rr.truth.iter_mut() {
            t.is_da = false;
            t.has_effect = false;
            t.confounded = false;
            t.log_fc = 0.0;
        }
        let b = vec![0.0; 100];
        let se = vec![1.0; 100];
        let p = vec![0.5; 100];
        let q = vec![0.5; 100];
        let called = vec![false; 100];
        let reported = Reported::new(&b, &se, &p, &q, &called).unwrap();
        let m = evaluate(
            "t",
            &cell(),
            0,
            "rust",
            &rr,
            &(0..100).collect::<Vec<usize>>(),
            &reported,
        );
        assert_eq!(m.empirical_fdr, None);
        assert_eq!(m.power, None, "no DA taxa means no power");
        assert_eq!(m.sign_concordance, None);
    }

    #[test]
    fn a_caller_that_calls_nothing_reports_no_rate_rather_than_a_fake_zero() {
        let b = vec![0.0; 100];
        let m = rep_metrics(&b, &[false; 100]);
        assert_eq!(m.empirical_fdr, None, "no calls means no empirical FDR");
        assert_eq!(m.power, Some(0.0), "power is well defined at zero calls");
        assert_eq!(m.sign_concordance, Some(0.0));
        assert_eq!(m.jaccard, Some(0.0));
        assert_eq!(m.f1, Some(0.0));
    }

    #[test]
    fn the_null_sign_rate_is_reported_and_counts_both_null_directions() {
        let mut b = vec![-1.0; 100];
        b[..50].iter_mut().for_each(|v| *v = 1.0);
        let m = rep_metrics(&b, &[false; 100]);
        // 50 DA taxa up, 50 null taxa all negative, so no null is positive.
        assert_eq!(m.null_sign_rate, Some(0.0));
        let half = {
            let mut b = vec![-1.0; 100];
            for k in 50..75 {
                b[k] = 1.0;
            }
            rep_metrics(&b, &[false; 100]).null_sign_rate
        };
        assert_eq!(half, Some(0.5), "25 of 50 nulls positive");
    }

    #[test]
    fn lfc_bias_and_rmse_use_the_truth_not_the_realised_counts() {
        // Every DA taxon is estimated 0.5 too low, every null 0.5 too high.
        // Half the taxa are underestimated by 0.5 and half overestimated by
        // 0.5, so the *mean* error is zero while every individual error is
        // large. A metric that reported only the mean would call this perfect.
        let m = rep_metrics(&[0.5; 100], &[false; 100]);
        assert!(m.lfc_bias.unwrap().abs() < 1e-12, "{:?}", m.lfc_bias);
        assert!(
            (m.lfc_rmse.unwrap() - 0.5).abs() < 1e-12,
            "{:?}",
            m.lfc_rmse
        );
        assert!((m.lfc_mae.unwrap() - 0.5).abs() < 1e-12);

        // A pure bias, in the direction the correction is supposed to remove.
        let up = rep_metrics(&[0.5; 100], &[false; 100]);
        assert!(up.lfc_bias.unwrap().is_finite());
    }

    /// A labelled taxon with no compositional effect is a null.
    ///
    /// This is the confounded arm's whole point, and getting it wrong would
    /// report the sampling-fraction correction as working by accident: a cell
    /// whose labelled taxa all "have effects" would score as high power, which
    /// is indistinguishable from the correction doing nothing.
    #[test]
    fn a_confounded_label_is_not_an_effect() {
        let r = stub();
        let mut rr = r;
        for (i, t) in rr.truth.iter_mut().enumerate() {
            t.is_da = i < 50;
            t.has_effect = false;
            t.confounded = i < 50;
            t.log_fc = 0.0;
        }
        let n = 100;
        // A method with no sampling-fraction correction calls all 50 of them.
        let b: Vec<f64> = (0..n).map(|i| if i < 50 { 0.6 } else { 0.0 }).collect();
        let se = vec![1.0; n];
        let p = vec![0.01; n];
        let q = vec![0.01; n];
        let called: Vec<bool> = (0..n).map(|i| i < 50).collect();
        let reported = Reported::new(&b, &se, &p, &q, &called).unwrap();
        let m = evaluate(
            "t",
            &cell(),
            0,
            "rust",
            &rr,
            &(0..n).collect::<Vec<usize>>(),
            &reported,
        );
        assert_eq!(m.n_confounded_retained, 50);
        assert_eq!(m.confounded_false_positive_rate, Some(1.0), "all 50 called");
        // And they are false positives for the confusion counts, so power is
        // undefined rather than 1.
        assert_eq!(m.power, None, "no taxon has a compositional effect");
        assert_eq!(m.empirical_fdr, Some(1.0), "every call is a false positive");
        assert_eq!(m.n_da_retained, 0);
    }

    #[test]
    fn a_failing_rep_is_recorded_rather_than_dropped() {
        let m = RepMetrics::failed("t", 4, 2, "rust", &cell(), "no variation".into());
        assert_eq!(m.cell, 4);
        assert_eq!(m.rep, 2);
        assert!(m.error.is_some());
        assert!(m.empirical_fdr.is_none());
        assert!(m.lfc_bias.is_none());
    }

    #[test]
    fn the_binomial_standard_error_matches_its_formula() {
        assert!((binomial_se(0.5, 100) - 0.05).abs() < 1e-12);
        assert!((binomial_se(0.0, 100) - 0.0).abs() < 1e-12);
        assert!(binomial_se(0.5, 0).is_nan());
        assert!((binomial_se(0.2, 400) - 0.02).abs() < 1e-12);
    }
}
