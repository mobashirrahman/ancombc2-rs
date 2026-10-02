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

    /// Mean of each row over the *finite* entries, R's `rowMeans(x, na.rm = TRUE)`.
    pub fn row_means_na_rm(&self) -> Vec<f64> {
        (0..self.rows)
            .map(|i| {
                let mut s = 0.0;
                let mut n = 0usize;
                for &v in self.row(i) {
                    if !v.is_nan() {
                        s += v;
                        n += 1;
                    }
                }
                if n == 0 {
                    f64::NAN
                } else {
                    s / n as f64
                }
            })
            .collect()
    }

    /// Subtract a per-row constant in place: `x[i, ] -= v[i]`.
    ///
    /// This is the per-taxon centring step, `Y[i, ] - rowMeans(Y)[i]`, which
    /// removes taxon-specific sequencing efficiency. It is *row*-wise, unlike
    /// the theta adjustment below, which is column-wise.
    pub fn sub_rows_in_place(&mut self, v: &[f64]) {
        debug_assert_eq!(v.len(), self.rows);
        for (i, &vv) in v.iter().enumerate() {
            let r = self.row_mut(i);
            for x in r.iter_mut() {
                *x -= vv;
            }
        }
    }

    /// Subtract a per-column constant in place: `x[, j] -= v[j]`.
    pub fn sub_cols_in_place(&mut self, v: &[f64]) {
        debug_assert_eq!(v.len(), self.cols);
        for i in 0..self.rows {
            let r = self.row_mut(i);
            for (j, &vv) in v.iter().enumerate() {
                r[j] -= vv;
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

    /// Replace non-finite entries with `v`, R's `x[is.infinite(x)] <- NA` plus
    /// the compatibility quirk in the sandwich accumulation.
    pub fn map_nonfinite_in_place(&mut self, f: impl Fn(f64) -> f64) {
        for v in self.data.iter_mut() {
            if !v.is_finite() {
                *v = f(*v);
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
        let mut m = RMatrix::zeros(2, 3);
        m.row_mut(0).copy_from_slice(&[1.0, 2.0, f64::NAN]);
        m.row_mut(1).copy_from_slice(&[3.0, f64::NAN, 6.0]);
        let means = m.row_means_na_rm();
        assert!((means[0] - 1.5).abs() < 1e-15);
        assert!((means[1] - 4.5).abs() < 1e-15);
    }

    #[test]
    fn all_nan_row_gives_nan_mean() {
        let mut m = RMatrix::zeros(1, 2);
        m.row_mut(0).copy_from_slice(&[f64::NAN, f64::NAN]);
        assert!(m.row_means_na_rm()[0].is_nan());
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
