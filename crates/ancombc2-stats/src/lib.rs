//! R-compatible statistical primitives.
//!
//! Every routine here is a reimplementation of a specific `stats` function, and
//! each one carries a test that pins its behaviour to R. "R-compatible" is meant
//! literally: ANCOM-BC2's numerical output depends on R's exact conventions
//! (`quantile` type 7, `p.adjust` step-down ordering, `pt` upper tail), so a
//! mathematically equivalent but differently-conventioned routine is a
//! compatibility bug, not a harmless detail.

pub mod dist;
pub mod nelder_mead;
pub mod padjust;
pub mod quantile;

pub use dist::{
    chisq_lower, chisq_upper, dnorm, erf, erfc, f_upper, pgamma_lower, pgamma_upper, t_two_sided,
    t_upper,
};
pub use nelder_mead::{nelder_mead, nelder_mead_1d, nelder_mead_scalar, NmResult, NmStatus};
pub use padjust::{p_adjust, p_adjust_n, AdjustMethod};
pub use quantile::{quantile_type7, quantile_type7_unsorted, var_unbiased};
