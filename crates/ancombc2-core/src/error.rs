//! Typed errors.
//!
//! Every failure mode the reference surfaces as `stop()` has a variant here, and
//! the messages name the covariate or sample involved the way the reference
//! does. A panic here would be a bug: the property suite asserts that degenerate
//! inputs produce `Err`, not a crash.

use crate::config::CompatMode;

pub type Result<T> = std::result::Result<T, AncombcError>;

#[derive(Debug, thiserror::Error)]
pub enum AncombcError {
    /// An input the core could not accept: a malformed table, an unreadable
    /// formula, an unrecognised option name. The IO layer and the CLI produce
    /// these; the numerics never do, so a `BadInput` means the problem is
    /// upstream of the algorithm rather than in it.
    #[error("{0}")]
    BadInput(String),

    #[error("shape mismatch: expected {expected} elements, got {got}")]
    Shape { expected: usize, got: usize },

    #[error("taxa x samples: {n_taxa} taxa and {n_samp} samples, at least 1 of each required")]
    EmptyInput { n_taxa: usize, n_samp: usize },

    #[error(
        "estimation failed for the following covariates: {covariates}\n\
         please ensure that these covariates do not have missing values and check \
         for multicollinearity before re-estimating the model"
    )]
    UnidentifiableCovariates { covariates: String },

    #[error(
        "no residual degrees of freedom: the model is over-parameterized; \
         please consider a more parsimonious model"
    )]
    NoResidualDegreesOfFreedom,

    #[error("zero variances have been detected for the following taxa: {taxa}\nplease remove these taxa or select a more parsimonious model")]
    ZeroVariance { taxa: String },

    #[error("no taxa remain under the current cutoff")]
    NoTaxaRemain,

    #[error("all taxa contain structural zeros")]
    AllTaxaStructuralZeros,

    #[error("no samples remain under the current cutoff")]
    NoSamplesRemain,

    #[error(
        "estimation of sampling fractions failed for the following samples: {samples}\n\
         these samples may have an excessive number of zero values"
    )]
    SamplingFractionFailed { samples: String },

    #[error("the group variable should have >= 2 categories")]
    GroupTooFewLevels,

    #[error(
        "group variable is required for the multi-group comparison: `group` is unset \
         while `global` or `pairwise` is enabled"
    )]
    GroupRequired,

    #[error(
        "group variable is required for detecting structural zeros; \
         set `struc_zero = false` to proceed"
    )]
    GroupRequiredForStructuralZeros,

    #[error("sample size per group should be >= 2; small sample size detected for: {groups}")]
    GroupTooSmall { groups: String },

    #[error("design matrix has {n_cols} columns but the expected {expected} were supplied")]
    DesignColumnsMismatch { n_cols: usize, expected: usize },

    #[error("group_labels has {got} entries but the design has {expected} samples")]
    GroupLabelCount { got: usize, expected: usize },

    #[error("pseudo must be non-negative")]
    NegativePseudo,

    #[error("compatibility mode {mode:?} does not support this operation")]
    UnsupportedInCompat { mode: CompatMode },
}

impl From<crate::matrix::MatrixError> for AncombcError {
    fn from(e: crate::matrix::MatrixError) -> Self {
        match e {
            crate::matrix::MatrixError::Shape { expected, got } => {
                AncombcError::Shape { expected, got }
            }
        }
    }
}
