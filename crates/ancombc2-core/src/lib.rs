//! `ancombc2-core`: a numerically compatible fixed-effects ANCOM-BC2 kernel.
//!
//! # Oracle
//!
//! ```text
//! ANCOMBC 2.15.2 @ dc4febdf59badb3a8dfe0c767ef2186323c2199a
//! ```
//!
//! # Scope
//!
//! Matrices, fixed effects, bias-corrected estimation, global and pairwise tests,
//! and pseudo-count sensitivity analysis. Random effects, Dunnett and trend
//! testing are out of scope for v1.0.
//!
//! # Design rules
//!
//! 1. **No R, no I/O.** This crate takes `f64` slices and a config struct and
//!    returns plain data. That is what keeps a Python binding, a CLI, or a WASM
//!    build possible.
//! 2. **Determinism.** Every parallel reduction accumulates in a fixed index
//!    order. A result must not depend on the thread count.
//! 3. **Quirks are reproduced, then labelled.** Where the reference does
//!    something mathematically questionable, [`CompatMode::Ancombc2_15`] (the
//!    default) reproduces it and [`CompatMode::StrictSpec`] implements the
//!    intended behaviour. See `docs/reference_behavior.md`.

// Several hot loops in this crate index two or three parallel buffers by index
// (`beta[t * p + a]`, `x.get(j, a)`, ...). Rewriting them as zipped iterators
// would be no faster and would hide which buffer is which, so the lint is allowed
// crate-wide rather than at each site.
#![allow(clippy::needless_range_loop)]

pub mod config;
pub mod correct;
pub mod em;
pub mod error;
pub mod matrix;
pub mod mle;
pub mod parallel;
pub mod pipeline;
pub mod preprocess;
pub mod sens;
pub mod stats;
pub mod test_mod;
pub mod vcov;
pub mod workspace;

pub use error::{AncombcError, Result};
pub use pipeline::{
    ancombc2_run, ancombc2_run_named, AncombcResult, CoreOutput, GlobalTest, PairwiseTest,
    ZeroIndication,
};
pub use preprocess::CountMatrix;
pub use sens::SensitivityScores;
pub use workspace::Workspace;

pub use config::{AdjustMethod, AncombcConfig, CompatMode, EmControl, IterControl, MdfdrControl};

/// The exact oracle this crate claims compatibility with.
///
/// Kept as a constant so the claim is greppable and testable rather than a
/// comment that rots.
pub const ORACLE_PACKAGE: &str = "ANCOMBC";
pub const ORACLE_VERSION: &str = "2.15.2";
pub const ORACLE_SHA: &str = "dc4febdf59badb3a8dfe0c767ef2186323c2199a";

/// The compatibility target string quoted in the README, the crate metadata and
/// the CI badge.
pub const COMPAT_TARGET: &str =
    "ancombc2-rs v0.1 == ANCOMBC 2.15.2 @ dc4febdf59badb3a8dfe0c767ef2186323c2199a";

/// Stage-resolved timings, for benchmarking. `Instant` is not used so the struct
/// is `PartialEq`-free and cheap to accumulate.
#[derive(Debug, Clone, Default)]
pub struct Timings {
    pub preprocess_s: f64,
    pub mle_s: f64,
    pub sandwich_s: f64,
    pub em_s: f64,
    pub correction_s: f64,
    pub inference_s: f64,
    pub global_s: f64,
    pub pairwise_s: f64,
    pub sensitivity_s: f64,
    pub sensitivity_runs: usize,
    pub ml_iterations: usize,
    pub em_iterations: Vec<usize>,
}

impl Timings {
    /// Total core time, excluding serialisation.
    pub fn total_s(&self) -> f64 {
        self.preprocess_s
            + self.mle_s
            + self.sandwich_s
            + self.em_s
            + self.correction_s
            + self.inference_s
            + self.global_s
            + self.pairwise_s
            + self.sensitivity_s
    }
}
