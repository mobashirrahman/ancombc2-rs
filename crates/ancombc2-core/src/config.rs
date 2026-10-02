//! Configuration mirroring the `ancombc2()` signature, restricted to the
//! fixed-effects path.
//!
//! Field names and defaults are taken from the pinned reference
//! (`ANCOMBC 2.15.2`), so a configuration written against the R API translates
//! one-for-one. Arguments the fixed-effects path ignores (`rand_formula`,
//! `dunnet`, `trend`, `lme_control`, `trend_control`) are absent rather than
//! accepted-and-ignored, so a caller cannot believe they are in effect.

pub use ancombc2_stats::AdjustMethod;

pub use crate::preprocess::Representation;

/// Which behaviour to implement where the reference is questionable.
///
/// The default is [`CompatMode::Ancombc2_15`]: reproduce the reference. A
/// divergence from the intended mathematics is a documented finding, not
/// something to paper over silently.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum CompatMode {
    /// Reproduce ANCOMBC 2.15.2 exactly, quirks included.
    #[default]
    Ancombc2_15,
    /// Implement the mathematically intended behaviour where the two differ.
    StrictSpec,
}

/// Iterative-MLE controls. Defaults `tol = 1e-2`, `max_iter = 20`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct IterControl {
    pub tol: f64,
    pub max_iter: usize,
}

impl Default for IterControl {
    fn default() -> Self {
        Self {
            tol: 1e-2,
            max_iter: 20,
        }
    }
}

/// E-M controls. Defaults `tol = 1e-5`, `max_iter = 100`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EmControl {
    pub tol: f64,
    pub max_iter: usize,
}

impl Default for EmControl {
    fn default() -> Self {
        Self {
            tol: 1e-5,
            max_iter: 100,
        }
    }
}

/// mdFDR controls. `B` is only used by the Dunnett and trend tests (out of scope
/// for v1.0); the pairwise test uses only `fwer_ctrl_method`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MdfdrControl {
    pub fwer_ctrl_method: AdjustMethod,
    pub b: usize,
}

impl Default for MdfdrControl {
    fn default() -> Self {
        Self {
            fwer_ctrl_method: AdjustMethod::Holm,
            b: 100,
        }
    }
}

/// The full configuration of one ANCOM-BC2 run.
#[derive(Debug, Clone, PartialEq)]
pub struct AncombcConfig {
    /// How the count table is stored.
    ///
    /// `Dense` (the default) is what PLAN.md section 9 concludes for the hot
    /// loops: the MLE, the sandwich and the group regressions are dense kernels
    /// over row subsets, so a sparse layer adds an indirection and removes
    /// nothing. `SparseTaxa` compresses to the stored entries and is a genuine
    /// win only for the *screening* stages -- prevalence and library size -- and
    /// only above roughly 33% zeros, below which `12` bytes per stored entry costs
    /// more than the `8` it replaces. The flag exists so that claim is measured
    /// rather than argued; see `docs/compatibility.md`.
    pub representation: Representation,
    /// Retain the full-size intermediate tables in [`crate::pipeline::CoreOutput`]: `o1`, `o2`,
    /// `y1`, `y2` and `y_bias_crt`.
    ///
    /// These five are `n_taxa x n_samp` each, so together they are five times
    /// the count matrix -- 4 GB on `bm5`. They are part of the golden contract
    /// ("processed counts, centered Y"), so the parity harness compares them, and
    /// *nothing else reads them*: not the CLI, not the FFI, and not the
    /// sensitivity analysis, which re-runs the whole pipeline and would hold five
    /// more per refit. So they are off by default and the harness turns them on.
    ///
    /// With the flag set, [`crate::pipeline::CoreOutput`]'s five matrices are as before. With it
    /// clear they are empty (`0 x 0`) placeholders, which is detectable rather
    /// than silently wrong -- a caller that needs them without asking gets a
    /// length mismatch, not zeros that look like data.
    pub keep_intermediates: bool,
    /// Fixed-effect column names, in design-matrix order. May be empty; names
    /// are carried through to the result tables for interoperability.
    pub fix_eff: Vec<String>,
    /// p-value adjustment method.
    pub p_adj_method: AdjustMethod,
    /// Added to every count before the log transform.
    pub pseudo: f64,
    /// Run the pseudo-count sensitivity analysis.
    pub pseudo_sens: bool,
    /// `true` reruns the whole algorithm per pseudo-count (0.1, 0.5, 1);
    /// `false` refits inference 50 times with pseudo-counts 0.01..0.5.
    pub conservative: bool,
    /// Prevalence cutoff; taxa below it are dropped.
    pub prv_cut: f64,
    /// Library-size cutoff; samples below it are dropped.
    pub lib_cut: f64,
    /// Quantile of the SE distribution used for the `s0` regulariser.
    pub s0_perc: f64,
    /// Name of the grouping variable, needed for structural zeros and the
    /// multi-group tests. `None` disables them.
    pub group: Option<String>,
    /// Group labels, one per sample, in sample order.
    pub group_labels: Option<Vec<String>>,
    /// Detect structural zeros.
    pub struc_zero: bool,
    /// Classify structural zeros by the asymptotic lower bound.
    pub neg_lb: bool,
    /// Significance level.
    pub alpha: f64,
    /// Run the global test.
    pub global: bool,
    /// Run the pairwise directional test.
    pub pairwise: bool,
    pub iter_control: IterControl,
    pub em_control: EmControl,
    pub mdfdr_control: MdfdrControl,
    pub compat: CompatMode,
}

impl Default for AncombcConfig {
    fn default() -> Self {
        Self {
            representation: Representation::Dense,
            keep_intermediates: false,
            fix_eff: Vec::new(),
            p_adj_method: AdjustMethod::Holm,
            pseudo: 0.0,
            pseudo_sens: false,
            conservative: true,
            prv_cut: 0.10,
            lib_cut: 0.0,
            s0_perc: 0.05,
            group: None,
            group_labels: None,
            struc_zero: false,
            neg_lb: false,
            alpha: 0.05,
            global: false,
            pairwise: false,
            iter_control: IterControl::default(),
            em_control: EmControl::default(),
            mdfdr_control: MdfdrControl::default(),
            compat: CompatMode::Ancombc2_15,
        }
    }
}

impl AncombcConfig {
    /// The pseudo-count grid used by the conservative sensitivity analysis:
    /// the main run plus 0.1, 0.5 and 1.
    pub fn conservative_pseudo_grid(&self) -> Vec<f64> {
        vec![0.0, 0.1, 0.5, 1.0]
    }

    /// The pseudo-count grid used by the non-conservative analysis: 0.01
    /// through 0.5 in steps of 0.01, as the reference builds it with
    /// `seq(0.01, 0.5, 0.01)`.
    pub fn nonconservative_pseudo_grid(&self) -> Vec<f64> {
        (1..=50).map(|i| (i as f64) / 100.0).collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn defaults_match_the_reference_signature() {
        let c = AncombcConfig::default();
        assert_eq!(c.p_adj_method, AdjustMethod::Holm);
        assert_eq!(c.pseudo, 0.0);
        assert_eq!(c.prv_cut, 0.10);
        assert_eq!(c.lib_cut, 0.0);
        assert_eq!(c.s0_perc, 0.05);
        assert_eq!(c.alpha, 0.05);
        assert!(!c.global);
        assert!(!c.pairwise);
        assert_eq!(c.iter_control.tol, 1e-2);
        assert_eq!(c.iter_control.max_iter, 20);
        assert_eq!(c.em_control.tol, 1e-5);
        assert_eq!(c.em_control.max_iter, 100);
        assert_eq!(c.mdfdr_control.b, 100);
        assert_eq!(c.compat, CompatMode::Ancombc2_15);
    }

    #[test]
    fn pseudo_grids_match_the_reference() {
        let c = AncombcConfig::default();
        assert_eq!(c.conservative_pseudo_grid(), vec![0.0, 0.1, 0.5, 1.0]);
        let g = c.nonconservative_pseudo_grid();
        assert_eq!(g.len(), 50);
        assert!((g[0] - 0.01).abs() < 1e-15);
        assert!((g[49] - 0.50).abs() < 1e-15);
    }

    #[test]
    fn strict_spec_is_not_the_default() {
        assert_ne!(CompatMode::default(), CompatMode::StrictSpec);
    }
}
