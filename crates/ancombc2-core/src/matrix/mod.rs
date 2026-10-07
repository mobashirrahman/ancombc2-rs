//! Dense linear algebra, column-major, with R/LAPACK-compatible semantics.

mod dqrls;
mod linalg;
mod pattern;

pub use dqrls::{
    blas_kind, dqrdc2, dqrls, dqrls_multi, dqrls_multi_selected, set_blas_kind, Blas, BlasKind,
    DqrlsFit, DqrlsMulti, OpenBlasHaswell, RefBlas,
};
pub use linalg::{cholesky_solve, eigen_symmetric, ginv, qr, Matrix, MatrixError, Qr, LM_FIT_TOL};
pub use pattern::{group as group_by_observation, group_by_pattern, Bitset, PatternGroups};
