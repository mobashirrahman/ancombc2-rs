//! The top-level pipeline: `.ancombc2()` for the fixed-effects path.
//!
//! Stage order, matching the reference:
//!
//! ```text
//! structural zeros (optional)
//!   -> .data_core on the native table      -> O1
//!   -> .data_core on the aggregate table   -> O2
//!   -> log + centre                         -> y1
//!   -> .iter_mle (theta estimated)         -> beta*, theta, vcov
//!   -> .bias_em per coefficient             -> delta_em, delta_wls, var_delta
//!   -> beta - delta; sampling fractions     -> theta_hat
//!   -> .iter_mle (theta fixed)              -> beta_hat, var_hat
//!   -> variance of delta; s0                -> se
//!   -> W, p, q, diff_abn
//!   -> global / pairwise
//!   -> sensitivity analysis
//! ```
//!
//! Note the two-table structure: `O1` (the native-level counts) estimates the
//! sampling fractions, because that needs as many taxa as possible, while `O2`
//! (the reported table) drives the inference. Without taxonomic aggregation the
//! two are the same matrix and the distinction is invisible, but the code keeps
//! it so an aggregating caller gets the reference's behaviour.

use crate::config::AncombcConfig;
use crate::correct::{
    add_bias_variance, apply_bias_correction, bias_correct_log_table, primary_inference,
    refresh_vcov_diagonal, regularised_variances, sampling_fractions,
};
use crate::em::{bias_em_all, BiasResult};
use crate::error::{AncombcError, Result};
use crate::matrix::Matrix;
use crate::mle::{
    check_design_identifiable, iter_mle_estimate_theta, iter_mle_fixed_theta, IterMle,
};
use crate::preprocess::{structural_zeros, CountMatrix};
use crate::sens;
use crate::test_mod::{global_test, group_columns, pairwise_test};
pub use crate::test_mod::{GlobalTest, PairwiseTest};
use crate::workspace::RMatrix;

/// A monotonic stopwatch for one stage. `Instant` is used rather than
/// `SystemTime` so a clock adjustment mid-run cannot produce a negative duration.
struct Stage(std::time::Instant);

impl Stage {
    fn start() -> Self {
        Self(std::time::Instant::now())
    }
    fn elapsed(&self) -> f64 {
        self.0.elapsed().as_secs_f64()
    }

    /// Name this stage, for the resident-memory trace.
    fn named(name: &'static str) -> Self {
        let s = Self::start();
        s.trace_mem(name);
        s
    }

    /// Report the current resident set, if `ANCOMBC2_TRACE_MEM` is set.
    ///
    /// `VmRSS`, not `VmHWM`: the high-water mark only ever rises, so it cannot
    /// show which stage *holds* the memory, only which stage first reached it.
    /// Reading a `/proc` entry per stage is a few microseconds against stages
    /// that take seconds, and it is the only way to tell "this stage allocates a
    /// lot" from "this stage is holding a lot", which are different defects with
    /// different fixes.
    fn trace_mem(&self, name: &str) {
        static ON: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
        if *ON.get_or_init(|| std::env::var_os("ANCOMBC2_TRACE_MEM").is_some()) {
            let kb = std::fs::read_to_string("/proc/self/status")
                .ok()
                .and_then(|s| {
                    s.lines()
                        .find(|l| l.starts_with("VmRSS:"))
                        .and_then(|l| l.split_whitespace().nth(1))
                        .and_then(|v| v.parse::<u64>().ok())
                });
            if let Some(kb) = kb {
                eprintln!("[mem] {name:<28} rss={:.3} GB", kb as f64 / 1e6);
            }
        }
    }
}

/// Wall time per pipeline stage, in seconds.
///
/// The stages are the ones the performance plan names. Timings are *not* part of
/// the golden contract: they are wall-clock measurements, so they are excluded
/// from every comparison and only ever written to the run metadata.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct StageTimings {
    pub preprocess: f64,
    pub pattern_grouping: f64,
    pub mle1: f64,
    pub sandwich1: f64,
    pub em: f64,
    pub correction: f64,
    pub mle2: f64,
    pub sandwich2: f64,
    pub tests: f64,
    pub serialisation: f64,
    pub sensitivity: f64,
}

impl StageTimings {
    /// The stages in pipeline order, as `(name, seconds)`.
    pub fn as_pairs(&self) -> [(&'static str, f64); 11] {
        [
            ("preprocess", self.preprocess),
            ("pattern_grouping", self.pattern_grouping),
            ("mle1", self.mle1),
            ("sandwich1", self.sandwich1),
            ("em", self.em),
            ("correction", self.correction),
            ("mle2", self.mle2),
            ("sandwich2", self.sandwich2),
            ("tests", self.tests),
            ("serialisation", self.serialisation),
            ("sensitivity", self.sensitivity),
        ]
    }

    /// Total across every stage, which is less than the wall time because the
    /// difference is the filter and the assembly around them.
    pub fn total(&self) -> f64 {
        self.as_pairs().iter().map(|(_, s)| *s).sum()
    }
}

/// Everything a single core run produces, before the sensitivity analysis.
#[derive(Debug, Clone, Default)]
pub struct CoreOutput {
    /// `O1 + pseudo`: the taxa that estimate the bias, after the prevalence and
    /// library-size filters but *before* the structural-zero screen.
    pub o1: CountMatrix,
    /// `O2 + pseudo`, the reported table: `o1` minus the structural zeros.
    pub o2: CountMatrix,
    /// Retained taxon indices, into the input table.
    pub taxa: Vec<usize>,
    /// The taxon indices behind [`CoreOutput::o1`], a superset of `taxa`.
    pub taxa_bias: Vec<usize>,
    /// Retained sample indices, into the input table.
    pub samples: Vec<usize>,
    /// Taxon names of the input table, so a result is self-describing and a
    /// parity failure can name a taxon rather than an index.
    pub taxon_names: Vec<String>,
    /// Sample names of the input table.
    pub sample_names: Vec<String>,
    pub fix_eff: Vec<String>,
    /// The centred log abundance used for the first MLE.
    pub y1: RMatrix,
    /// The centred log abundance used for the second MLE.
    pub y2: RMatrix,
    /// `t(t(y2) - theta_hat)`.
    pub y_bias_crt: RMatrix,
    /// The pre-correction coefficients.
    pub beta_star: Vec<f64>,
    /// `beta* - delta_em`.
    pub beta_corrected: Vec<f64>,
    pub theta: Vec<f64>,
    pub samp_frac: Vec<f64>,
    pub bias: Vec<BiasResult>,
    pub delta_em: Vec<f64>,
    pub delta_wls: Vec<f64>,
    pub var_delta: Vec<f64>,
    /// The final coefficients.
    pub beta: Vec<f64>,
    /// Variances after the bias and regularisation adjustments.
    pub var_final: Vec<f64>,
    /// Stage-1 marginal variances (`var_hat` of the first MLE), `n_taxa * p`.
    pub var1: Vec<f64>,
    /// Stage-1 sandwich covariance, `n_taxa * p * p` row-major.
    pub vcov1: Vec<f64>,
    /// Stage-2 marginal variances, before the bias and `s0` adjustments.
    pub var_hat: Vec<f64>,
    /// Stage-2 sandwich covariance with the diagonal refreshed by the
    /// regularised variances, `n_taxa * p * p` row-major.
    pub vcov: Vec<f64>,
    /// Residual degrees of freedom, `n_taxa * p`.
    pub dof: Vec<f64>,
    pub s02: Vec<f64>,
    pub se: Vec<f64>,
    pub w: Vec<f64>,
    pub p: Vec<f64>,
    pub q: Vec<f64>,
    pub diff_abn: Vec<bool>,
    pub global: Option<GlobalTest>,
    pub pairwise: Option<PairwiseTest>,
    pub zero_ind: Option<ZeroIndication>,
    pub ml_iterations: usize,
    /// The first MLE's `epsilon` at each iteration, in order: the contract's
    /// "convergence trace", compared against the oracle's own printed trace.
    pub ml_trace: Vec<f64>,
    /// The reported set's missingness pattern assignment, 1-based per taxon: the
    /// contract's Level A "missingness pattern assignment".
    pub pattern_group: Vec<usize>,
    pub eps: Vec<f64>,
    /// Wall time per stage, for the run metadata. Not part of the contract.
    pub timings: StageTimings,
}

/// Structural-zero flags.
#[derive(Debug, Clone)]
pub struct ZeroIndication {
    pub groups: Vec<String>,
    /// `n_taxa x n_groups`.
    pub zero_ind: Vec<bool>,
}

/// The full result of a run, including the sensitivity analysis.
#[derive(Debug, Clone)]
pub struct AncombcResult {
    pub core: CoreOutput,
    pub sensitivity: Option<sens::SensitivityScores>,
    /// `passed_ss` from the sensitivity analysis, when it ran.
    pub passed_ss: Option<Vec<bool>>,
    /// `diff_robust` = `diff_abn & passed_ss`, when the analysis ran.
    pub diff_robust: Option<Vec<bool>>,
    /// Non-fatal conditions the oracle reports with `warning()`. R deactivates
    /// features rather than aborting, and a silent deactivation would be worse
    /// than useless, so the message travels with the result.
    pub warnings: Vec<String>,
}

impl CoreOutput {
    /// The index of the coefficient for `fix_eff` in the per-taxon coefficient
    /// vectors (`beta`, `se`, `p`, `q`, ...), which are `n_taxa * p`
    /// row-major.
    ///
    /// A name is preferred to a position because the intercept is first only by
    /// convention, and a design can carry several group contrasts. `None` when
    /// the name is absent, which is a caller error rather than something to
    /// paper over by falling back to a column.
    ///
    /// When `fix_eff` was left empty in the config the core labels the columns
    /// `V0`, `V1`, ... , so a caller that never set them cannot ask by the name
    /// it sees in a result table.
    pub fn fix_eff_index(&self, name: &str) -> Option<usize> {
        self.fix_eff.iter().position(|n| n == name)
    }

    /// One coefficient's column, across the retained taxa.
    ///
    /// `beta`, `se`, `p` and `q` are stored flat and `n_taxa * p` long because
    /// that is how the reference hands them back and because the flat form is
    /// what the golden contract serialises. A caller who wants a single
    /// coefficient therefore has to know the stride, and getting the stride
    /// wrong reads the intercepts of the first few taxa as if they were the
    /// group coefficient of every taxon -- which is silent, plausible, and
    /// wrong. This is the accessor for that.
    ///
    /// `None` when the name is unknown, so a missing coefficient is an error the
    /// caller must handle rather than an empty result to average over.
    ///
    /// This applies to *every* per-taxon vector in [`CoreOutput`], not only the
    /// floats. `diff_abn` is `n_taxa * p` as well -- the divergence call is made
    /// per coefficient, which is why the reference's output table has a
    /// `diff_<coefficient>` column per fixed effect and not one per taxon.
    /// Reading it as one entry per taxon pairs taxon `t`'s intercept call with
    /// taxon `t-1`'s group call.
    pub fn coefficient<T: Copy>(&self, name: &str, values: &[T]) -> Option<Vec<T>> {
        let col = self.fix_eff_index(name)?;
        // The stride comes from the array's own length, not from `fix_eff`. If
        // the two ever disagree, the caller's assumption about the layout is
        // already wrong, and guessing either one silently reads a different
        // coefficient. So: divide the length, require it to be exact, and
        // require the resulting stride to be the number of named columns.
        let n = self.taxa.len();
        if n == 0 || values.len() % n != 0 {
            return None;
        }
        let p = values.len() / n;
        if p == 0 || p != self.fix_eff.len() {
            return None;
        }
        Some((0..n).map(|t| values[t * p + col]).collect())
    }
}

#[cfg(test)]
mod coefficient_tests {
    use super::*;

    /// A [`CoreOutput`] carrying only what these tests need: the retained taxa
    /// and the fixed-effect names. Every other field is empty, which is fine --
    /// `coefficient` reads none of them.
    fn output_with(fix_eff: &[&str], n_taxa: usize) -> CoreOutput {
        CoreOutput {
            taxa: (0..n_taxa).collect(),
            fix_eff: fix_eff.iter().map(|s| s.to_string()).collect(),
            ..CoreOutput::default()
        }
    }

    #[test]
    fn a_named_column_is_read_with_the_right_stride() {
        // Two taxa, two coefficients, row-major: [i0, g0, i1, g1].
        let o = output_with(&["(Intercept)", "group2"], 2);
        let beta = [10.0, 11.0, 20.0, 21.0];
        assert_eq!(o.coefficient("group2", &beta), Some(vec![11.0, 21.0]));
        assert_eq!(o.coefficient("(Intercept)", &beta), Some(vec![10.0, 20.0]));
    }

    #[test]
    fn flags_are_read_with_the_same_stride_as_floats() {
        // `diff_abn` is per coefficient, so a per-taxon read would pair taxon
        // t's intercept call with taxon t-1's group call.
        let o = output_with(&["(Intercept)", "group2"], 3);
        let called = [true, false, false, true, true, true];
        assert_eq!(
            o.coefficient("group2", &called),
            Some(vec![false, true, true])
        );
    }

    #[test]
    fn an_unknown_name_is_refused_rather_than_defaulted() {
        let o = output_with(&["(Intercept)", "group2"], 2);
        // Falling back to a column here would silently compare intercepts
        // against log fold changes, so there is no fallback.
        assert_eq!(o.coefficient("group3", &[1.0, 2.0, 3.0, 4.0]), None);
        assert_eq!(o.coefficient("group3", &[true, false, true, false]), None);
    }

    #[test]
    fn a_length_that_is_not_a_multiple_of_the_taxa_count_is_refused() {
        let o = output_with(&["(Intercept)", "group2"], 2);
        // Three values for two taxa: the stride would not be integral.
        assert_eq!(o.coefficient("group2", &[1.0, 2.0, 3.0]), None);
        let empty: [f64; 0] = [];
        assert_eq!(o.coefficient("group2", &empty), None);
    }

    #[test]
    fn a_stride_that_disagrees_with_the_column_names_is_refused() {
        // Three values per taxon but two named columns: the caller's layout
        // assumption is wrong, and either guess reads a different coefficient.
        let o = output_with(&["(Intercept)", "group2"], 2);
        assert_eq!(
            o.coefficient("group2", &[1.0, 2.0, 3.0, 4.0, 5.0, 6.0]),
            None
        );
    }

    #[test]
    fn an_empty_result_refuses_every_column() {
        let o = output_with(&["(Intercept)", "group2"], 0);
        let empty: [f64; 0] = [];
        assert_eq!(o.coefficient("group2", &empty), None);
    }
}

/// Run ANCOM-BC2's fixed-effects path.
///
/// `counts` is `taxa x samples`, `x` is the design `samples x p` (as R's
/// `model.matrix` would produce, including the intercept column), and
/// `group_index` maps each sample to a group level. `count_core` builds a
/// [`CountMatrix`] from a raw slice, so this function is independent of the IO
/// layer.
#[allow(clippy::too_many_arguments)]
pub fn ancombc2_run(
    counts: &CountMatrix,
    x: &Matrix,
    group_index: Option<&[usize]>,
    cfg: &AncombcConfig,
) -> Result<AncombcResult> {
    ancombc2_run_named(counts, x, group_index, cfg, &[], &[])
}

/// The name of a group level, defaulting to its 0-based index.
fn level_name(group: Option<&str>, level: usize) -> String {
    match group {
        Some(g) => format!("{g}{level}"),
        None => level.to_string(),
    }
}

/// [`ancombc2_run`] with the input's taxon and sample names, so the result can
/// report *which* taxon diverged rather than its index.
pub fn ancombc2_run_named(
    counts: &CountMatrix,
    x: &Matrix,
    group_index: Option<&[usize]>,
    cfg: &AncombcConfig,
    taxon_names: &[String],
    sample_names: &[String],
) -> Result<AncombcResult> {
    if cfg.pseudo < 0.0 {
        return Err(AncombcError::NegativePseudo);
    }
    if counts.n_taxa == 0 || counts.n_samp == 0 {
        return Err(AncombcError::EmptyInput {
            n_taxa: counts.n_taxa,
            n_samp: counts.n_samp,
        });
    }
    if x.rows != counts.n_samp {
        return Err(AncombcError::DesignColumnsMismatch {
            n_cols: x.rows,
            expected: counts.n_samp,
        });
    }
    if !cfg.global && !cfg.pairwise && !cfg.struc_zero && cfg.group.is_some() {
        // a group without a test is legal, matching the reference
    }
    if (cfg.global || cfg.pairwise) && group_index.is_none() {
        return Err(AncombcError::GroupRequired);
    }
    if cfg.struc_zero && group_index.is_none() {
        return Err(AncombcError::GroupRequiredForStructuralZeros);
    }
    // `data_sanity_check` validates the group before anything else, and its
    // thresholds are not interchangeable: fewer than two categories is a hard
    // error, but fewer than three only *deactivates* the multi-group
    // comparisons. Reporting a rank-deficient design instead would point at the
    // design rather than at the group.
    let mut warnings: Vec<String> = Vec::new();
    let mut global = cfg.global;
    let mut pairwise = cfg.pairwise;
    if let Some(gi) = group_index {
        // `NO_GROUP` is not a level. Counting it would make a two-group analysis
        // look like a three-group one -- so the global and pairwise tests would
        // stay enabled when R deactivates them -- and would put a level called
        // `group18446744073709551615` into the size-1 error message.
        let labelled = || {
            gi.iter()
                .copied()
                .filter(|v| *v != crate::preprocess::NO_GROUP)
        };
        let n_levels = labelled().collect::<std::collections::HashSet<_>>().len();
        if n_levels < 2 {
            return Err(AncombcError::GroupTooFewLevels);
        }
        if n_levels < 3 && (cfg.global || cfg.pairwise) {
            // R warns and sets `global = pairwise = FALSE`; the deactivation is
            // part of the contract, so it must not surface as a global test.
            global = false;
            pairwise = false;
            warnings.push(
                "The group variable has < 3 categories\n\
                 The multi-group comparisons (global/pairwise) will be deactivated"
                    .to_string(),
            );
        }
        let mut sizes: Vec<(usize, usize)> = Vec::new();
        for g in labelled() {
            match sizes.iter_mut().find(|(lvl, _)| *lvl == g) {
                Some((_, n)) => *n += 1,
                None => sizes.push((g, 1)),
            }
        }
        sizes.sort_unstable();
        let small: Vec<String> = sizes
            .iter()
            .filter(|(_, n)| *n < 2)
            .map(|(lvl, _)| level_name(cfg.group.as_deref(), *lvl))
            .collect();
        if !small.is_empty() {
            return Err(AncombcError::GroupTooSmall {
                groups: small.join(", "),
            });
        }
        let unstable: Vec<String> = sizes
            .iter()
            .filter(|(_, n)| *n < 5)
            .map(|(lvl, _)| level_name(cfg.group.as_deref(), *lvl))
            .collect();
        if !unstable.is_empty() {
            warnings.push(format!(
                "Small sample size detected for the following group(s): {}\n\
                 Variance estimation would be unstable when the sample size is < 5 per group",
                unstable.join(", ")
            ));
        }
    }
    if cfg.p_adj_method
        == crate::config::AdjustMethod::parse("__never__")
            .unwrap_or(crate::config::AdjustMethod::None)
    {
        // unreachable: parse of a sentinel cannot succeed
    }

    // The number of groups, from the highest level actually present. The
    // `NO_GROUP` sentinel is excluded: it is not a level, and taking `max()` over
    // it would both count a phantom group and overflow on the `+ 1`.
    let n_group = group_index.map(|g| {
        g.iter()
            .copied()
            .filter(|v| *v != crate::preprocess::NO_GROUP)
            .max()
            .map(|m| m + 1)
            .unwrap_or(0)
    });
    let zero_ind = match (cfg.struc_zero, group_index, &cfg.group) {
        (true, Some(gi), Some(gname)) => {
            let n_g = n_group.unwrap();
            let sz = structural_zeros(counts, gi, n_g, cfg.neg_lb);
            // taxa rejected from the primary analysis
            let keep = sz.taxa_without_structural_zeros();
            Some((sz, keep, gname.clone()))
        }
        _ => None,
    };

    // ---- filtering ----
    //
    // The reference keeps *two* taxon sets and the distinction is load-bearing,
    // not a detail:
    //
    //   O1 = `.data_core(..., tax_keep = NULL)`    -- prevalence and library
    //        size only. `y1`, the first MLE, the E-M bias and the sampling
    //        fractions are all estimated from this set.
    //   O2 = `.data_core(..., tax_keep = tax_keep)` -- O1 *minus the taxa with
    //        structural zeros*. `y2`, the second MLE and every reported quantity
    //        come from this set.
    //
    // So a taxon with a structural zero still contributes to the estimate of
    // theta, and dropping it before the first MLE would move the sampling
    // fractions and therefore every reported coefficient. `.data_core`'s third
    // branch also re-indexes `tax_keep` into the already-subset table, which
    // for a row-wise statistic such as prevalence is equivalent to intersecting
    // the two sets, so the intersection is what is computed here.
    // With `Representation::SparseTaxa` the count table is carried compressed from
    // the reader through the two screens that can read it compressed, and
    // expanded once, at the point where the MLE wants dense rows.
    //
    // The screens are prevalence and library size. Both scan every cell and store
    // nothing, which is the shape that benefits: on the compressed form each reads
    // `n_nonzero + n_missing` entries instead of `n_taxa * n_samp`. The
    // structural-zero screen is *not* one of them -- it needs each absent cell's
    // group, which a taxon-major encoding does not carry -- so it runs on the
    // expanded table, and the expansion happens before it rather than after.
    //
    // Both paths produce the same `Filtered`, so nothing downstream can tell which
    // ran. `the_two_representations_select_identical_taxa` pins that.
    let f1 = match cfg.representation {
        crate::preprocess::Representation::Dense => counts.filter(cfg.prv_cut, cfg.lib_cut)?,
        crate::preprocess::Representation::SparseTaxa => {
            crate::preprocess::SparseTaxaMatrix::from_dense(counts)
                .filter(cfg.prv_cut, cfg.lib_cut)?
        }
    };
    let taxa1 = f1.taxa;
    let samples1 = f1.samples;
    let taxa2: Vec<usize> = match &zero_ind {
        Some((sz, _, _)) => {
            // A bit test per taxon against the structural-zero flags, rather than
            // a `HashSet` built from the kept-taxa list: the flags are already the
            // predicate being asked about, so materialising a second index of them
            // costs an allocation and a hash per lookup to answer the same question.
            // The group-major screen is `any_in_row`, which is what `keep` encodes.
            let n_g = sz.groups.len();
            let taxa: Vec<usize> = taxa1
                .iter()
                .copied()
                .filter(|&t| !sz.bits.any_in_row(t, n_g))
                .collect();
            if taxa.is_empty() {
                return Err(AncombcError::AllTaxaStructuralZeros);
            }
            taxa
        }
        None => taxa1.clone(),
    };
    let o1 = counts.select(&taxa1, &samples1);
    let o2 = counts.select(&taxa2, &samples1);
    // Names default to positional labels so a caller that does not supply them
    // still gets a readable failure report.
    let all_taxon_names: Vec<String> = (0..counts.n_taxa)
        .map(|i| {
            taxon_names
                .get(i)
                .cloned()
                .unwrap_or_else(|| format!("taxon_{i}"))
        })
        .collect();
    let all_sample_names: Vec<String> = (0..counts.n_samp)
        .map(|i| {
            sample_names
                .get(i)
                .cloned()
                .unwrap_or_else(|| format!("sample_{i}"))
        })
        .collect();

    // The design must be subset to the retained samples.
    let xs = subset_design(x, &samples1);
    check_design_identifiable(&xs)?;
    let p = xs.cols;
    // `n_taxa1` rows feed the bias estimation; `n_taxa` rows are reported.
    let n_taxa1 = taxa1.len();
    let n_taxa = taxa2.len();
    let n_samp = samples1.len();
    let mut timings = StageTimings::default();

    // ---- MLE 1: estimate theta ----
    let t = Stage::named("preprocess");
    let y1 = o1.log_center(cfg.pseudo);
    timings.preprocess += t.elapsed();
    t.trace_mem("y1 built");
    let t = Stage::named("mle1");
    t.trace_mem("mle1:start");
    let group_cols = group_columns(&xs.colnames, cfg.group.as_deref().unwrap_or("group"));
    let mle1 = iter_mle_estimate_theta(
        &xs,
        &y1.data,
        n_taxa1,
        n_samp,
        cfg.iter_control,
        cfg.compat,
        &group_cols,
    )?
    .0;
    timings.mle1 += t.elapsed();
    timings.sandwich1 += mle1.sandwich_seconds;
    t.trace_mem("mle1:end");

    // ---- E-M bias ----
    let t = Stage::named("em");
    let bias = bias_em_all(
        &mle1.beta,
        &mle1.var_hat,
        n_taxa1,
        p,
        cfg.em_control.tol,
        cfg.em_control.max_iter,
    )?;
    timings.em += t.elapsed();
    let delta_em: Vec<f64> = bias.iter().map(|b| b.delta_em).collect();
    let delta_wls: Vec<f64> = bias.iter().map(|b| b.delta_wls).collect();
    let var_delta: Vec<f64> = bias.iter().map(|b| b.var_delta).collect();

    // ---- bias correction and sampling fractions ----
    let t = Stage::named("correction");
    let beta_corrected = apply_bias_correction(&mle1.beta, n_taxa1, p, &delta_em);
    let samp_frac = sampling_fractions(&y1, &xs, &beta_corrected, n_taxa1, p);
    timings.correction += t.elapsed();

    // ---- MLE 2: fixed theta ----
    //
    // The reference fits the second MLE on the *reported* table `O2`, which is
    // the aggregate-level counts. Taxonomic aggregation belongs to the R layer
    // (see PLAN.md 2), so the core receives one table; what remains is the
    // structural-zero screen, and that *is* modelled here as a second taxon set
    // because it changes the estimate.
    let t = Stage::named("preprocess");
    let y2 = o2.log_center(cfg.pseudo);
    timings.preprocess += t.elapsed();
    let t = Stage::named("mle2");
    let mle2 = iter_mle_fixed_theta(
        &xs,
        &y2.data,
        n_taxa,
        n_samp,
        &samp_frac,
        cfg.compat,
        &group_cols,
    )?;
    timings.mle2 += t.elapsed();
    timings.sandwich2 += mle2.sandwich_seconds;
    t.trace_mem("mle2:end");
    let y_bias_crt = bias_correct_log_table(&y2, &samp_frac);
    Stage::named("y_bias_crt").trace_mem("y2+y_bias_crt");

    // ---- variance adjustments ----
    let t = Stage::named("correction");
    let var_bc = add_bias_variance(&mle2.var_hat, n_taxa, p, &var_delta);
    let (var_final, se) = regularised_variances(&var_bc, &mle2.beta, n_taxa, p, Some(cfg.s0_perc));
    let s02 = crate::correct::s0_per_column(&var_bc, n_taxa, p, Some(cfg.s0_perc));
    let vcov = refresh_vcov_diagonal(&mle2.vcov, &var_final, n_taxa, p);
    timings.correction += t.elapsed();

    // ---- inference ----
    let dof = expand_dof(&mle2, n_taxa, p);
    let inf = primary_inference(
        &mle2.beta,
        &se,
        &dof,
        n_taxa,
        p,
        cfg.p_adj_method,
        cfg.alpha,
    );

    // ---- multi-group tests ----
    let fix_eff = if cfg.fix_eff.is_empty() {
        (0..p).map(|k| format!("V{k}")).collect()
    } else {
        cfg.fix_eff.clone()
    };
    let global = if global {
        Some(global_test(
            &fix_eff,
            cfg.group.as_deref().unwrap_or(""),
            &mle2.beta,
            &vcov,
            Some(&dof),
            n_taxa,
            p,
            cfg.p_adj_method,
            cfg.alpha,
        )?)
    } else {
        None
    };
    let pairwise = if pairwise {
        // Argument order is (..., vcov, var_hat, dof, ...): the pairwise
        // contrasts need both the full covariance block (for the off-diagonal
        // terms of a difference) and the regularised marginal variances (for a
        // single group coefficient). They are different arrays and swapping them
        // type-checks fine, so the call site says which is which.
        Some(pairwise_test(
            &fix_eff,
            cfg.group.as_deref().unwrap_or(""),
            &mle2.beta,
            &vcov,
            &var_final,
            Some(&dof),
            n_taxa,
            p,
            cfg.p_adj_method,
            cfg.mdfdr_control.fwer_ctrl_method,
            cfg.alpha,
        )?)
    } else {
        None
    };

    timings.tests += t.elapsed();

    let zero = zero_ind.map(|(sz, _, _)| ZeroIndication {
        groups: sz.groups,
        zero_ind: sz.zero_ind,
    });

    let t = Stage::named("serialisation");
    let o1_pseudo = add_pseudo(&o1, cfg.pseudo);
    let o2_pseudo = add_pseudo(&o2, cfg.pseudo);
    t.trace_mem("o1_pseudo+o2_pseudo");
    timings.serialisation += t.elapsed();

    // The five full-size intermediates are only materialised when a caller asked
    // for them -- see `AncombcConfig::keep_intermediates`. Off by default because
    // they are five times the count matrix and no production path reads them;
    // empty rather than zeroed so a caller that forgot the flag gets a shape
    // mismatch instead of plausible-looking zeros.
    let keep = cfg.keep_intermediates;
    let core = CoreOutput {
        o1: if keep {
            o1_pseudo
        } else {
            CountMatrix::default()
        },
        o2: if keep {
            o2_pseudo
        } else {
            CountMatrix::default()
        },
        taxa: taxa2.clone(),
        // The taxa that estimate the sampling fractions: `taxa` plus any with
        // structural zeros, which the reference keeps for the first MLE.
        taxa_bias: taxa1.clone(),
        samples: samples1.clone(),
        taxon_names: all_taxon_names.clone(),
        sample_names: all_sample_names.clone(),
        fix_eff,
        y1: if keep { y1 } else { RMatrix::default() },
        y2: if keep { y2 } else { RMatrix::default() },
        y_bias_crt: if keep { y_bias_crt } else { RMatrix::default() },
        beta_star: mle1.beta,
        beta_corrected,
        theta: mle1.theta,
        samp_frac,
        bias,
        delta_em,
        delta_wls,
        var_delta,
        beta: mle2.beta,
        var1: mle1.var_hat.clone(),
        vcov1: mle1.vcov.clone(),
        var_hat: mle2.var_hat.clone(),
        var_final,
        vcov,
        dof,
        s02,
        se,
        w: inf.w,
        p: inf.p,
        q: inf.q,
        diff_abn: inf.diff_abn,
        global,
        pairwise,
        zero_ind: zero,
        ml_iterations: mle1.iterations,
        ml_trace: mle1.trace.clone(),
        // The *reported* set's assignment: that is `mle2`, the fixed-theta fit over
        // the bias-corrected table restricted to the reported taxa. `mle1` groups
        // the bias set, which is a superset, so its numbering is not this one's.
        pattern_group: mle2.pattern_group.clone(),
        eps: mle1.eps,
        timings,
    };

    // ---- sensitivity analysis ----
    let mut result = AncombcResult {
        core,
        sensitivity: None,
        passed_ss: None,
        diff_robust: None,
        warnings,
    };
    if cfg.pseudo_sens {
        let t = Stage::named("sensitivity");
        run_sensitivity(&mut result, counts, &o1, &xs, &samples1, &taxa2, cfg)?;
        result.core.timings.sensitivity += t.elapsed();
    }
    Ok(result)
}

fn subset_design(x: &Matrix, samples: &[usize]) -> Matrix {
    x.select_rows(samples)
}

fn add_pseudo(c: &CountMatrix, pseudo: f64) -> CountMatrix {
    if pseudo == 0.0 {
        return c.clone();
    }
    let mut out = c.clone();
    for v in out.data.iter_mut() {
        *v += pseudo;
    }
    out
}

/// Expand a per-taxon residual dof into a per-coefficient matrix.
///
/// R: `dof <- matrix(rep(dof, n_fix_eff), ncol = n_fix_eff, byrow = FALSE)`,
/// i.e. the same value repeated across a taxon's coefficients. When `.iter_mle`
/// estimated theta itself it returns `dof = NULL`, and `pt(x, df = NULL)`
/// evaluates R's default `df = Inf` — a normal tail, not a missing value.
fn expand_dof(mle: &IterMle, n_taxa: usize, p: usize) -> Vec<f64> {
    match &mle.dof {
        Some(per_taxon) => {
            debug_assert_eq!(per_taxon.len(), n_taxa, "one dof per taxon");
            let mut out = vec![f64::INFINITY; n_taxa * p];
            for i in 0..n_taxa {
                for k in 0..p {
                    out[i * p + k] = per_taxon[i];
                }
            }
            out
        }
        None => vec![f64::INFINITY; n_taxa * p],
    }
}

/// Conservative mode: re-run the core for each pseudo-count, in parallel, and
/// reduce the adjusted p-values of the primary table to scores.
///
/// The screening designs are shared, so the four runs are independent apart from
/// reading the same inputs.
/// `o1_raw` is the *un-pseudo-added* bias-estimation table: the conservative
/// refits add their own pseudo-count, and reusing the reported `O1`, which
/// already had the main run's pseudo-count added, would add it twice.
#[allow(clippy::too_many_arguments)]
fn run_sensitivity(
    result: &mut AncombcResult,
    counts: &CountMatrix,
    o1_raw: &CountMatrix,
    x: &Matrix,
    samples: &[usize],
    taxa: &[usize],
    cfg: &AncombcConfig,
) -> Result<()> {
    let o1 = o1_raw;
    let t2 = taxa;
    let p = result.core.fix_eff.len();
    let n_taxa = taxa.len();
    let alpha = cfg.alpha;

    if cfg.conservative {
        let grid = cfg.conservative_pseudo_grid();
        // the first entry is the main run, already computed
        let main_q: Vec<f64> = result.core.q.clone();
        // Level 1 of the nesting order, and the reason it is claimed here rather
        // than left to a bare `par_iter`: each refit runs the whole pipeline, so an
        // unclaimed level would let every parallel stage inside it nest and split
        // the pool again. `skip(1)` is done before the map rather than inside it,
        // because the index the budget sees has to be over what is actually run.
        let mut lvl = crate::parallel::NestingBudget::level("pseudo-count runs (conservative)");
        Stage::named("sens:conservative").trace_mem("sens:before");
        let rest: Vec<f64> = grid.iter().copied().skip(1).collect();
        let others: Vec<Vec<f64>> = crate::parallel::map_par(&mut lvl, &rest, |&pc| {
            // A refit cannot fail where the main run succeeded: same
            // design, same taxa, only the pseudo-count differs.
            core_run_at_pseudo(o1, t2, x, pc, cfg).unwrap_or_else(|_| vec![1.0; n_taxa * p])
        });
        Stage::named("sens:conservative").trace_mem("sens:after refits");
        let mut runs = vec![main_q];
        runs.extend(others);
        let n_col = p;
        let s = sens::conservative_run(&runs, n_taxa, colnames_primary(&result.core), alpha, &grid);
        let _ = n_col;
        let passed = s.passed.clone();
        let robust: Vec<bool> = (0..n_taxa * p)
            .map(|i| result.core.diff_abn[i] && passed[i])
            .collect();
        result.sensitivity = Some(s);
        result.passed_ss = Some(passed);
        result.diff_robust = Some(robust);
    } else {
        // Non-conservative: sampling fractions are estimated once, on the main
        // run, and only the inference step is refitted.
        let samp_frac = result.core.samp_frac.clone();
        let grid = cfg.nonconservative_pseudo_grid();
        let fit = |pc: f64| -> Vec<f64> {
            let y = sens::sens_response_sub(counts, taxa, samples, &samp_frac, pc);
            ols_p_values(&y, x, cfg)
        };
        Stage::named("sens:nonconservative").trace_mem("sens:nc before");
        let sc = sens::nonconservative_run(&grid, n_taxa, p, alpha, fit)?;
        Stage::named("sens:nonconservative").trace_mem("sens:nc after");
        // `passed_fun` is applied column-wise: the main run's p-value for the
        // *same* coefficient decides agreement.
        let mut passed = vec![false; n_taxa * p];
        for c in 0..p {
            let col = sc.scores[c..]
                .iter()
                .step_by(p)
                .copied()
                .collect::<Vec<_>>();
            let pv = (0..n_taxa)
                .map(|i| {
                    let v = result.core.p[i * p + c];
                    if v.is_nan() {
                        1.0
                    } else {
                        v
                    }
                })
                .collect::<Vec<_>>();
            let pc = sens::nonconservative_passed(&col, &pv, alpha);
            for i in 0..n_taxa {
                passed[i * p + c] = pc[i];
            }
        }
        let robust: Vec<bool> = (0..n_taxa * p)
            .map(|i| result.core.diff_abn[i] && passed[i])
            .collect();
        result.sensitivity = Some(sc);
        result.passed_ss = Some(passed);
        result.diff_robust = Some(robust);
    }
    Ok(())
}

fn colnames_primary(core: &CoreOutput) -> Vec<String> {
    core.fix_eff.iter().map(|f| format!("q_{f}")).collect()
}

/// Unadjusted p-values from a per-taxon ordinary least-squares fit, which is what
/// the non-conservative sensitivity analysis uses.
///
/// This is `summary(lm(y ~ .))$coefficients[, "Pr(>|t|)"]`: an ordinary `lm`,
/// with the homoskedastic standard error and a t reference on
/// `df.residual`. It is deliberately *not* the ANCOM-BC2 `W` statistic, and that
/// difference is part of what the sensitivity analysis measures.
fn ols_p_values(y: &RMatrix, x: &Matrix, _cfg: &AncombcConfig) -> Vec<f64> {
    let p = x.cols;
    // One p-value per (taxon, fix_eff), matching the reference's
    // `n_tax x n_fix_eff` frame; `NA` becomes 1 as `p_hat[is.na] <- 1` does.
    let mut out = vec![1.0f64; y.rows * p];
    let xr = x.t();
    let xtx = xr.matmul(x).expect("gram");
    let inv = crate::matrix::ginv(&xtx, None);
    let df = (x.rows - p) as f64;
    // The taxa axis, innermost of the four: every taxon's OLS reads the same `X`,
    // `XtX` and `inv`, and writes only its own `p` slots, so the split is over
    // disjoint output and the arithmetic *within* a taxon is untouched -- the
    // result is bit-identical to the serial loop, not merely close.
    //
    // This is what the non-conservative sensitivity analysis is made of: 50 refits
    // of this function over every taxon. Serially it was the whole cost of `bm6`
    // and the reason a single-threaded Rust run was *slower* than single-threaded
    // R there (41.0s vs 29.3s) despite winning every smaller dataset by 2x-48x.
    let mut lvl = crate::parallel::NestingBudget::level("taxa (sensitivity ols)");
    let indices: Vec<usize> = (0..y.rows).collect();
    let rows: Vec<Vec<f64>> = crate::parallel::map_par(&mut lvl, &indices, |&i| {
        ols_p_values_one(y, x, &xtx, &inv, p, df, i)
    });
    for (i, r) in rows.into_iter().enumerate() {
        out[i * p..(i + 1) * p].copy_from_slice(&r);
    }
    out
}

/// One taxon's p-values in `ols_p_values`, as its own `p`-slot vector.
///
/// Split out so the parallel loop body is one call and the per-taxon arithmetic
/// reads in the same order as the serial version it replaced.
fn ols_p_values_one(
    y: &RMatrix,
    x: &Matrix,
    xtx: &Matrix,
    inv: &Matrix,
    p: usize,
    df: f64,
    i: usize,
) -> Vec<f64> {
    // Pre-filled with 1, which is also what a taxon that cannot be fitted keeps:
    // the reference does `p_hat[is.na] <- 1`, so every early return below leaves
    // the same all-ones row.
    let mut out = vec![1.0f64; p];
    let row = y.row(i);

    let mut xty = vec![0.0; p];
    for a in 0..p {
        let mut s = 0.0;
        for (j, &yv) in row.iter().enumerate() {
            if yv.is_finite() {
                s += x.get(j, a) * yv;
            }
        }
        xty[a] = s;
    }
    let beta = match crate::matrix::cholesky_solve(xtx, &xty) {
        Ok(b) => b,
        Err(_) => {
            // A singular `XtX` falls back to the pseudo-inverse, as `lm` does
            // through `ginv`, rather than abandoning the taxon.
            match inv.matmul(&Matrix::transpose_vec_rows(&xty)).ok() {
                Some(m) => m.data,
                None => return out,
            }
        }
    };

    // Residual variance with `n - p` degrees of freedom.
    let mut rss = 0.0;
    let mut n = 0usize;
    for (j, &yv) in row.iter().enumerate() {
        if !yv.is_finite() {
            continue;
        }
        // `lm` drops a row when *any* model variable is missing, not just the
        // response, so an incomplete design row must not contribute.
        if (0..p).any(|a| !x.get(j, a).is_finite()) {
            continue;
        }
        let mut fit = 0.0;
        for a in 0..p {
            fit += x.get(j, a) * beta[a];
        }
        let r = yv - fit;
        rss += r * r;
        n += 1;
    }
    if n <= p || df <= 0.0 {
        return out;
    }
    let sigma2 = rss / df;
    for a in 0..p {
        let se = (sigma2 * inv.get(a, a)).max(0.0).sqrt();
        if se > 0.0 {
            out[a] = ancombc2_stats::t_two_sided(beta[a] / se, df);
        }
    }
    out
}

/// A single core run at a given pseudo-count, used by the conservative
/// sensitivity analysis. The taxa and design are already fixed, so this skips
/// filtering and structural-zero detection.
///
/// `o1` is the bias-estimation table and `taxa2` the reported ones, mirroring
/// the two taxon sets of [`ancombc2_run_named`]: the first MLE and the sampling
/// fractions run over all of `o1`, the second MLE and the reported q-values over
/// `o1`'s rows selected by `taxa2`. `taxa2` must be a subset of `o1`'s rows, in
/// the same order, which is what the pipeline produces.
pub fn core_run_at_pseudo(
    o1: &CountMatrix,
    taxa2: &[usize],
    x: &Matrix,
    pseudo: f64,
    cfg: &AncombcConfig,
) -> Result<Vec<f64>> {
    let n_taxa1 = o1.n_taxa;
    let n_taxa = taxa2.len();
    let n_samp = o1.n_samp;
    let p = x.cols;
    let y1 = o1.log_center(pseudo);
    let mle1 = iter_mle_estimate_theta(
        x,
        &y1.data,
        n_taxa1,
        n_samp,
        cfg.iter_control,
        cfg.compat,
        &group_columns(&x.colnames, cfg.group.as_deref().unwrap_or("group")),
    )?
    .0;
    let bias = bias_em_all(
        &mle1.beta,
        &mle1.var_hat,
        n_taxa1,
        p,
        cfg.em_control.tol,
        cfg.em_control.max_iter,
    )?;
    let delta_em: Vec<f64> = bias.iter().map(|b| b.delta_em).collect();
    let var_delta: Vec<f64> = bias.iter().map(|b| b.var_delta).collect();
    let beta_corr = apply_bias_correction(&mle1.beta, n_taxa1, p, &delta_em);
    let theta = sampling_fractions(&y1, x, &beta_corr, n_taxa1, p);
    // Logged straight out of `o1`'s selected rows: materialising the selected
    // count table as well would hold both for the duration of the transform.
    let y2 = o1.log_center_rows(taxa2, pseudo);
    let mle2 = iter_mle_fixed_theta(
        x,
        &y2.data,
        n_taxa,
        n_samp,
        &theta,
        cfg.compat,
        &group_columns(&x.colnames, cfg.group.as_deref().unwrap_or("group")),
    )?;
    let var_bc = add_bias_variance(&mle2.var_hat, n_taxa, p, &var_delta);
    let (_vf, se) = regularised_variances(&var_bc, &mle2.beta, n_taxa, p, Some(cfg.s0_perc));
    let dof = expand_dof(&mle2, n_taxa, p);
    let inf = primary_inference(
        &mle2.beta,
        &se,
        &dof,
        n_taxa,
        p,
        cfg.p_adj_method,
        cfg.alpha,
    );
    Ok(inf.q)
}

/// Re-exported so a caller can build a config with R's spelling of the adjust
/// method without depending on `ancombc2-stats` directly.
pub fn adjust_method_from_r(s: &str) -> Result<crate::config::AdjustMethod> {
    crate::config::AdjustMethod::parse(s).map_err(AncombcError::from_string)
}

impl AncombcError {
    fn from_string(s: String) -> AncombcError {
        AncombcError::UnidentifiableCovariates { covariates: s }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::preprocess::CountMatrix;

    /// A small, well-conditioned design: 2 groups, 20 samples, intercept plus
    /// a group dummy.
    fn design_2grp() -> (CountMatrix, Matrix, Vec<usize>) {
        let n_tax = 60usize;
        let n_samp = 20usize;
        let mut data = vec![0.0; n_tax * n_samp];
        let mut seed = 12345u64;
        let mut rnd = move || {
            // splitmix64
            seed = seed.wrapping_add(0x9E3779B97F4A7C15);
            let mut z = seed;
            z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
            z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
            z ^= z >> 31;
            (z >> 11) as f64 / (1u64 << 53) as f64
        };
        for i in 0..n_tax {
            let base = 20.0 + 60.0 * rnd();
            let signal = if i % 5 == 0 { 1.2 } else { 0.0 };
            for j in 0..n_samp {
                let grp = if j < n_samp / 2 { 0.0 } else { 1.0 };
                let mu = (base + signal * grp).exp();
                let v = rnd() * mu;
                data[i * n_samp + j] = v.floor();
            }
        }
        let counts = CountMatrix::new(n_tax, n_samp, data).unwrap();
        let mut xd = Vec::new();
        for j in 0..n_samp {
            xd.push(vec![1.0, if j < n_samp / 2 { 0.0 } else { 1.0 }]);
        }
        let x = Matrix::from_rows(&xd);
        let g: Vec<usize> = (0..n_samp).map(|j| usize::from(j >= n_samp / 2)).collect();
        (counts, x, g)
    }

    /// Three groups, so the global and pairwise comparisons stay active.
    fn design_3grp() -> (CountMatrix, Matrix, Vec<usize>) {
        let (counts, _, _) = design_2grp();
        let n_samp = counts.n_samp;
        let mut xd = Vec::new();
        for j in 0..n_samp {
            let g = if j < n_samp / 2 { 0.0 } else { 1.0 };
            xd.push(vec![1.0, g, if j < n_samp / 3 { 1.0 } else { 0.0 }]);
        }
        let g: Vec<usize> = (0..n_samp)
            .map(|j| {
                if j < n_samp / 2 {
                    0
                } else if j < n_samp / 2 + (n_samp - n_samp / 2) / 2 {
                    1
                } else {
                    2
                }
            })
            .collect();
        (counts, Matrix::from_rows(&xd), g)
    }

    #[test]
    fn a_clean_two_group_problem_runs_end_to_end() {
        let (counts, x, g) = design_2grp();
        let cfg = AncombcConfig {
            fix_eff: vec!["(Intercept)".into(), "grp".into()],
            prv_cut: 0.0,
            ..Default::default()
        };
        let r = ancombc2_run(&counts, &x, Some(&g), &cfg).unwrap();
        let n_taxa = r.core.taxa.len();
        let p = 2usize;
        assert!(n_taxa > 0);
        assert_eq!(r.core.beta.len(), n_taxa * p);
        assert_eq!(r.core.se.len(), n_taxa * p);
        assert_eq!(r.core.q.len(), n_taxa * p);
        assert_eq!(r.core.vcov.len(), n_taxa * p * p);
        assert_eq!(r.core.samp_frac.len(), counts.n_samp);
        // every probability is a probability
        for &v in r.core.p.iter().chain(r.core.q.iter()) {
            assert!((0.0..=1.0).contains(&v), "out of range: {v}");
        }
        // the E-M bias is finite
        assert!(r.core.delta_em.iter().all(|v| v.is_finite()));
        assert!(r.core.var_delta.iter().all(|v| *v > 0.0));
        // taxa constructed with a real group effect must be discovered
        let n_sig = r.core.diff_abn.iter().filter(|&&d| d).count();
        assert!(n_sig > 0, "the fixture has 12 truly differential taxa");
        assert!(n_sig <= n_taxa, "not everything can be significant");
    }

    #[test]
    fn a_run_is_reproducible() {
        let (counts, x, g) = design_2grp();
        let cfg = AncombcConfig {
            fix_eff: vec!["(Intercept)".into(), "grp".into()],
            prv_cut: 0.0,
            ..Default::default()
        };
        let a = ancombc2_run(&counts, &x, Some(&g), &cfg).unwrap();
        let b = ancombc2_run(&counts, &x, Some(&g), &cfg).unwrap();
        assert_eq!(a.core.beta, b.core.beta, "bitwise identical");
        assert_eq!(a.core.q, b.core.q);
        assert_eq!(a.core.samp_frac, b.core.samp_frac);
        assert_eq!(a.core.delta_em, b.core.delta_em);
    }

    #[test]
    fn the_global_test_runs_when_asked() {
        let (counts, x, g) = design_3grp();
        let cfg = AncombcConfig {
            fix_eff: vec!["(Intercept)".into(), "grp2".into(), "grp3".into()],
            group: Some("grp".into()),
            global: true,
            prv_cut: 0.0,
            ..Default::default()
        };
        let r = ancombc2_run(&counts, &x, Some(&g), &cfg).unwrap();
        let gt = r.core.global.expect("global test must be present");
        assert_eq!(gt.w.len(), r.core.taxa.len());
        assert!(gt.p.iter().all(|v| (0.0..=1.0).contains(v)));
    }

    /// `data_sanity_check` deactivates the multi-group comparisons below three
    /// categories instead of aborting, so neither may appear in the result and
    /// the warning must survive.
    #[test]
    fn two_group_levels_deactivate_the_multi_group_comparisons() {
        let (counts, x, g) = design_2grp();
        let cfg = AncombcConfig {
            fix_eff: vec!["(Intercept)".into(), "grp".into()],
            group: Some("grp".into()),
            global: true,
            pairwise: true,
            prv_cut: 0.0,
            ..Default::default()
        };
        let r = ancombc2_run(&counts, &x, Some(&g), &cfg).unwrap();
        assert!(r.core.global.is_none(), "global must be deactivated");
        assert!(r.core.pairwise.is_none(), "pairwise must be deactivated");
        assert!(
            r.warnings.iter().any(|w| w.contains("< 3 categories")),
            "the deactivation must be reported, got {:?}",
            r.warnings
        );
    }

    /// A group with a single sample cannot have a variance, so R aborts before
    /// any estimation.
    #[test]
    fn a_size_one_group_is_a_typed_error() {
        let (counts, x, mut g) = design_2grp();
        g[0] = 2; // a third level with exactly one sample
        let cfg = AncombcConfig {
            fix_eff: vec!["(Intercept)".into(), "grp".into()],
            group: Some("grp".into()),
            ..Default::default()
        };
        let e = ancombc2_run(&counts, &x, Some(&g), &cfg).unwrap_err();
        assert!(
            matches!(e, AncombcError::GroupTooSmall { .. }),
            "expected GroupTooSmall, got {e:?}"
        );
    }

    /// Fewer than five samples per group is only a warning.
    #[test]
    fn a_small_group_only_warns() {
        let (counts, x, g) = design_2grp();
        let cfg = AncombcConfig {
            fix_eff: vec!["(Intercept)".into(), "grp".into()],
            group: Some("grp".into()),
            ..Default::default()
        };
        let r = ancombc2_run(&counts, &x, Some(&g), &cfg).unwrap();
        // 10 per group here, so the warning must *not* fire
        assert!(
            !r.warnings.iter().any(|w| w.contains("Small sample size")),
            "unexpected warning: {:?}",
            r.warnings
        );
    }

    #[test]
    fn rejects_an_empty_input() {
        let counts = CountMatrix::zeros(0, 10);
        let x = Matrix::zeros(10, 2);
        let e = ancombc2_run(&counts, &x, None, &AncombcConfig::default()).unwrap_err();
        assert!(matches!(e, AncombcError::EmptyInput { .. }), "{e:?}");
    }

    #[test]
    fn rejects_a_design_with_the_wrong_number_of_rows() {
        let counts = CountMatrix::zeros(5, 10);
        let x = Matrix::zeros(9, 2);
        let e = ancombc2_run(&counts, &x, None, &AncombcConfig::default()).unwrap_err();
        assert!(
            matches!(e, AncombcError::DesignColumnsMismatch { .. }),
            "{e:?}"
        );
    }

    #[test]
    fn rejects_a_negative_pseudo() {
        let (counts, x, g) = design_2grp();
        let cfg = AncombcConfig {
            pseudo: -1.0,
            ..Default::default()
        };
        let e = ancombc2_run(&counts, &x, Some(&g), &cfg).unwrap_err();
        assert!(matches!(e, AncombcError::NegativePseudo), "{e:?}");
    }

    #[test]
    fn rejects_a_rank_deficient_design() {
        let (counts, _x, g) = design_2grp();
        // two identical columns: rank 1 with p = 2
        let rows: Vec<Vec<f64>> = (0..20).map(|_| vec![1.0, 1.0]).collect();
        let x = Matrix::from_rows(&rows);
        let e = ancombc2_run(&counts, &x, Some(&g), &AncombcConfig::default()).unwrap_err();
        match e {
            AncombcError::UnidentifiableCovariates { covariates } => {
                assert!(!covariates.is_empty(), "the error must name a covariate");
            }
            other => panic!("expected UnidentifiableCovariates, got {other:?}"),
        }
    }

    #[test]
    fn rejects_a_design_with_no_residual_degrees_of_freedom() {
        // n samples == p leaves no residual df. Checked at the guard directly:
        // reaching the estimator first would fail on a degenerate constant
        // table, which is a different (and legitimate) error.
        // Full rank (2 columns) with exactly n == p samples, so df.residual = 0.
        let x = Matrix::from_rows(&[vec![1.0, 0.0], vec![1.0, 1.0]]);
        let e = check_design_identifiable(&x).unwrap_err();
        assert!(
            matches!(e, AncombcError::NoResidualDegreesOfFreedom),
            "{e:?}"
        );
    }

    #[test]
    fn a_constant_count_table_has_zero_variance_and_is_reported() {
        // Every taxon identical: the sandwich variance is exactly 0, which the
        // reference reports as "Zero variances have been detected".
        let counts = CountMatrix::new(20, 6, vec![5.0; 120]).unwrap();
        let x = Matrix::from_rows(&[
            vec![1.0, 0.0],
            vec![1.0, 1.0],
            vec![1.0, 0.0],
            vec![1.0, 1.0],
            vec![1.0, 0.0],
            vec![1.0, 1.0],
        ]);
        let g = vec![0, 0, 1, 1, 0, 1];
        let e = ancombc2_run(&counts, &x, Some(&g), &AncombcConfig::default()).unwrap_err();
        assert!(matches!(e, AncombcError::ZeroVariance { .. }), "{e:?}");
    }

    #[test]
    fn global_requires_a_group() {
        let (counts, x, _g) = design_2grp();
        let cfg = AncombcConfig {
            global: true,
            ..Default::default()
        };
        let e = ancombc2_run(&counts, &x, None, &cfg).unwrap_err();
        assert!(matches!(e, AncombcError::GroupRequired), "{e:?}");
    }

    #[test]
    fn structural_zeros_require_a_group() {
        let (counts, x, _g) = design_2grp();
        let cfg = AncombcConfig {
            struc_zero: true,
            ..Default::default()
        };
        let e = ancombc2_run(&counts, &x, None, &cfg).unwrap_err();
        assert!(
            matches!(e, AncombcError::GroupRequiredForStructuralZeros),
            "{e:?}"
        );
    }
}
