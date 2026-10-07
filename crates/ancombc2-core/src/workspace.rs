//! Reusable buffers.
//!
//! The pipeline needs several `taxa x samples` temporaries at once: the centred
//! log-abundance matrix, the theta-adjusted response, the fitted values, the
//! residuals and their squares. Allocating five dense matrices per run, and again
//! per sensitivity refit, dominates the memory profile at 10,000 x 500. A single
//! [`Workspace`] per run is allocated once and reused in place.
//!
//! `taxa x samples` matrices are stored row-major here, because every hot loop
//! walks a taxon's samples contiguously (the per-taxon sandwich accumulation, the
//! per-taxon row mean, the per-taxon residual). The `matrix` module's dense
//! `Matrix` is column-major to match R; the two coexist and conversions happen
//! once, at the boundary.

/// Row-major `rows x cols` buffer.
///
/// The default is the empty `0 x 0` matrix. A matrix is a shape plus a buffer,
/// so this is the only value `Default` can mean, and having it lets a caller
/// build the struct with `..Default::default()`.
#[derive(Debug, Clone, Default)]
pub struct RMatrix {
    pub rows: usize,
    pub cols: usize,
    pub data: Vec<f64>,
}

impl RMatrix {
    pub fn zeros(rows: usize, cols: usize) -> Self {
        Self {
            rows,
            cols,
            data: vec![0.0; rows * cols],
        }
    }

    #[inline]
    pub fn get(&self, i: usize, j: usize) -> f64 {
        self.data[i * self.cols + j]
    }

    #[inline]
    pub fn set(&mut self, i: usize, j: usize, v: f64) {
        self.data[i * self.cols + j] = v;
    }

    #[inline]
    pub fn row(&self, i: usize) -> &[f64] {
        &self.data[i * self.cols..(i + 1) * self.cols]
    }

    #[inline]
    pub fn row_mut(&mut self, i: usize) -> &mut [f64] {
        let c = self.cols;
        &mut self.data[i * c..(i + 1) * c]
    }

    /// Column `j` as a slice of row-strided values is not contiguous; this walks
    /// it. Hot loops prefer [`RMatrix::row`].
    pub fn col(&self, j: usize) -> Vec<f64> {
        (0..self.rows).map(|i| self.get(i, j)).collect()
    }

    pub fn from_row_major(rows: usize, cols: usize, data: Vec<f64>) -> Self {
        assert_eq!(data.len(), rows * cols);
        Self { rows, cols, data }
    }

    /// Mean of each row over the non-`NA` entries.
    ///
    /// **Deliberately absent.** This method used to exist and sum in `f64`, and it
    /// was wrong: R accumulates `rowMeans` in `long double` (64 bits of mantissa on
    /// x86-64), and a cancelling sum loses bits the oracle keeps. There is now no
    /// way to reach a row mean without choosing a [`Reductions`], so the
    /// approximation cannot be used by accident.
    ///
    /// See [`crate::reduce`] for the measurement.
    /// Subtract a per-row constant in place: `x[i, ] -= v[i]`.
    ///
    /// This is the per-taxon centring step, `Y[i, ] - rowMeans(Y)[i]`, which
    /// removes taxon-specific sequencing efficiency. It is *row*-wise, unlike
    /// the theta adjustment below, which is column-wise.
    ///
    /// The subtraction is [`r_sub`], not `-`, so an `NA` cell keeps R's `NA_real_`
    /// payload. `Y[i, ] - v[i]` on a matrix containing `NA_real_` produces
    /// `0x7ff8000000000002` on x86-64 where R produces `NA_real_`; the difference is
    /// invisible to every numeric check and visible to `serialize()`.
    pub fn sub_rows_in_place(&mut self, v: &[f64]) {
        debug_assert_eq!(v.len(), self.rows);
        for (i, &vv) in v.iter().enumerate() {
            let r = self.row_mut(i);
            for x in r.iter_mut() {
                *x = crate::reduce::r_sub(*x, vv);
            }
        }
    }

    /// Subtract a per-column constant in place: `x[, j] -= v[j]`.
    pub fn sub_cols_in_place(&mut self, v: &[f64]) {
        debug_assert_eq!(v.len(), self.cols);
        for i in 0..self.rows {
            let r = self.row_mut(i);
            for (j, &vv) in v.iter().enumerate() {
                r[j] = crate::reduce::r_sub(r[j], vv);
            }
        }
    }

    /// Add a per-column constant in place: `x[, j] += v[j]`.
    pub fn add_cols_in_place(&mut self, v: &[f64]) {
        debug_assert_eq!(v.len(), self.cols);
        for i in 0..self.rows {
            let r = self.row_mut(i);
            for (j, &vv) in v.iter().enumerate() {
                r[j] += vv;
            }
        }
    }

    /// `x[, j] <- v[j] - x[, j]`, which is the reference's
    /// `t(t(y) - theta)`: note the *sign* of the subtraction, which is what
    /// makes the sensitivity-analysis response `y - samp_frac`.
    pub fn sub_cols_reversed_in_place(&mut self, v: &[f64]) {
        debug_assert_eq!(v.len(), self.cols);
        for i in 0..self.rows {
            let r = self.row_mut(i);
            for (j, &vv) in v.iter().enumerate() {
                r[j] = vv - r[j];
            }
        }
    }

    /// `x^2` in place.
    pub fn square_in_place(&mut self) {
        for v in self.data.iter_mut() {
            *v *= *v;
        }
    }

    /// `x[is.infinite(x)] <- NA`, which is what the reference does after every
    /// `log`.
    ///
    /// Two things this gets right that a `!is_finite()` test does not:
    ///
    /// * the mask is `is.infinite`, so a NaN of any payload is left alone — an
    ///   `NA_real_` in the counts keeps its own bits through `log`, and a computed
    ///   `NaN` is not turned into `NA_real_`;
    /// * the replacement is `NA_real_` (`0x7ff00000000007a2`), not `f64::NAN`
    ///   (`0x7ff8000000000000`).
    pub fn replace_infinite_with_na(&mut self) {
        for v in self.data.iter_mut() {
            if v.is_infinite() {
                *v = crate::reduce::na_real();
            }
        }
    }
}

/// The set of temporaries one core run needs.
#[derive(Debug)]
pub struct Workspace {
    /// Centred log abundance, `taxa x samples`.
    pub y: RMatrix,
    /// Fitted values from the last fit, `taxa x samples`.
    pub fitted: RMatrix,
    /// Residuals `eps = (y - fitted) - theta`, `taxa x samples`.
    pub eps: RMatrix,
    /// `eps^2`, reused by the sandwich accumulation.
    pub eps2: RMatrix,
    /// Design matrix rows, `samples x p` (column-major `Matrix`).
    pub design: crate::matrix::Matrix,
    /// Per-sample outer products `vec(x_j x_j')`, `samples x p^2` row-major.
    pub xx: Vec<f64>,
    /// Cached per-pattern QR factors.
    pub qrs: Vec<crate::matrix::Qr>,
    /// Per-pattern design row indices.
    pub pattern_rows: Vec<Vec<usize>>,
    /// Per-pattern taxon indices.
    pub pattern_taxa: Vec<Vec<usize>>,
    /// Sample-specific bias term, one per sample.
    pub theta: Vec<f64>,
    /// The fixed-effect estimate, `taxa x p` row-major.
    pub beta: Vec<f64>,
}

impl Workspace {
    pub fn new(n_taxa: usize, n_samples: usize, p: usize) -> Self {
        Self {
            y: RMatrix::zeros(n_taxa, n_samples),
            fitted: RMatrix::zeros(n_taxa, n_samples),
            eps: RMatrix::zeros(n_taxa, n_samples),
            eps2: RMatrix::zeros(n_taxa, n_samples),
            design: crate::matrix::Matrix::zeros(n_samples, p),
            xx: vec![0.0; n_samples * p * p],
            qrs: Vec::new(),
            pattern_rows: Vec::new(),
            pattern_taxa: Vec::new(),
            theta: vec![0.0; n_samples],
            beta: vec![0.0; n_taxa * p],
        }
    }

    /// Bytes held by the reusable buffers, for the memory benchmark.
    pub fn bytes(&self) -> usize {
        let mats = 4 * (self.y.data.len() * 8);
        let design = self.design.data.len() * 8;
        let xx = self.xx.len() * 8;
        // `Qr` no longer holds `Q`, so the footprint is the reflectors plus `R`.
        // Reporting the *live* buffers is the point of this function: the old
        // `n^2` `Q` is exactly what used to dominate it.
        let q = self
            .qrs
            .iter()
            .map(|q| (q.reflector_bytes() + q.r.data.len()) * 8)
            .sum::<usize>();
        mats + design + xx + q + self.theta.len() * 8 + self.beta.len() * 8
    }
}

/// `x^2` in place on a plain slice, leaving `NaN` as `NaN`.
///
/// A free function rather than an `RMatrix` method because the sandwich kernel
/// consumes a flat `eps^2` buffer.
pub fn square_nan(v: &mut [f64]) {
    for x in v.iter_mut() {
        *x = *x * *x;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn row_major_round_trip() {
        let m = RMatrix::from_row_major(2, 3, vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0]);
        assert_eq!(m.row(0), &[1.0, 2.0, 3.0]);
        assert_eq!(m.row(1), &[4.0, 5.0, 6.0]);
        assert_eq!(m.get(1, 2), 6.0);
        assert_eq!(m.col(1), vec![2.0, 5.0]);
    }

    #[test]
    fn row_means_ignore_nan() {
        use crate::reduce::{F64Reductions, Reductions};
        let mut m = RMatrix::zeros(2, 3);
        m.row_mut(0).copy_from_slice(&[1.0, 2.0, f64::NAN]);
        m.row_mut(1).copy_from_slice(&[3.0, f64::NAN, 6.0]);
        let means = F64Reductions.row_means_na_rm(&m);
        assert_eq!(means, vec![1.5, 4.5]);
    }

    /// An all-`NA` row is `NA_real_`, not `NaN`, and the two have different bits.
    #[test]
    fn an_all_missing_row_is_na_real_not_nan() {
        use crate::reduce::{F64Reductions, Reductions};
        let mut m = RMatrix::zeros(1, 2);
        m.row_mut(0).copy_from_slice(&[f64::NAN, f64::NAN]);
        let got = F64Reductions.row_means_na_rm(&m)[0];
        assert_eq!(got.to_bits(), crate::reduce::NA_REAL_BITS);
        assert_ne!(got.to_bits(), f64::NAN.to_bits());
    }

    /// `x[is.infinite(x)] <- NA`: `±Inf` becomes `NA_real_`, NaN is untouched.
    #[test]
    fn only_infinity_is_replaced_by_the_na_mask() {
        let mut m = RMatrix::from_row_major(
            1,
            5,
            vec![
                f64::NEG_INFINITY,
                f64::INFINITY,
                crate::reduce::na_real(),
                f64::NAN,
                1.5,
            ],
        );
        m.replace_infinite_with_na();
        let r = m.row(0);
        assert_eq!(r[0].to_bits(), crate::reduce::NA_REAL_BITS);
        assert_eq!(r[1].to_bits(), crate::reduce::NA_REAL_BITS);
        // An `NA_real_` that was already there keeps its payload.
        assert_eq!(r[2].to_bits(), crate::reduce::NA_REAL_BITS);
        // A computed NaN keeps the plain payload and is *not* promoted to
        // `NA_real_`.
        assert_eq!(r[3].to_bits(), f64::NAN.to_bits());
        assert_eq!(r[4], 1.5);
    }

    #[test]
    fn reversed_subtraction_flips_the_sign() {
        let mut m = RMatrix::zeros(1, 2);
        m.row_mut(0).copy_from_slice(&[1.0, 2.0]);
        m.sub_cols_reversed_in_place(&[10.0, 20.0]);
        assert_eq!(m.row(0), &[9.0, 18.0]);
        m.sub_cols_in_place(&[1.0, 1.0]);
        assert_eq!(m.row(0), &[8.0, 17.0]);
        m.add_cols_in_place(&[1.0, 1.0]);
        assert_eq!(m.row(0), &[9.0, 18.0]);
    }

    #[test]
    fn workspace_allocates_once() {
        let w = Workspace::new(100, 50, 3);
        assert_eq!(w.y.data.len(), 5000);
        assert_eq!(w.eps2.data.len(), 5000);
        assert_eq!(w.design.rows, 50);
        assert_eq!(w.design.cols, 3);
        assert_eq!(w.xx.len(), 50 * 9);
        assert!(w.bytes() > 0);
    }
}
