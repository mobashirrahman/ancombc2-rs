//! Dense linear algebra, column-major, with R/LAPACK-compatible semantics.

mod dqrls;
mod linalg;
mod pattern;

pub use dqrls::{dqrdc2, dqrls, Blas, DqrlsFit, OpenBlasHaswell, RefBlas};
pub use linalg::{cholesky_solve, eigen_symmetric, ginv, qr, Matrix, MatrixError, Qr};
pub use pattern::{group as group_by_observation, group_by_pattern, Bitset, PatternGroups};
