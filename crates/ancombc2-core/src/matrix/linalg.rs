//! Dense matrix primitives, matching R's storage order and LAPACK conventions.
//!
//! ANCOM-BC2's hot loops are small dense operations: a QR of an `n x p` design,
//! a `p x p` sandwich product per taxon, and a pseudo-inverse of `X'X`. Two
//! properties matter more than raw speed here:
//!
//!   * **Storage order.** R and LAPACK are column-major. Golden comparison and
//!     the `.sandwich_vcov` outer-product layout are both column-major, so this
//!     module is column-major to keep the correspondence auditable.
//!   * **`MASS::ginv` semantics.** The oracle computes
//!     `ginv(t(X) %*% X)`, an SVD pseudo-inverse with a relative cutoff, *not* an
//!     exact inverse. Reproducing the cutoff behaviour is required wherever a
//!     design is near-singular, which happens with a factor level present in one
//!     group only.

/// Column-major dense matrix.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Matrix {
    pub rows: usize,
    pub cols: usize,
    pub data: Vec<f64>,
    /// Column names, carried so an error message can name the offending
    /// covariate. Empty when unknown.
    pub colnames: Vec<String>,
    /// Row names, used to report retained samples and taxa in order.
    pub rownames: Vec<String>,
}

impl Matrix {
    pub fn zeros(rows: usize, cols: usize) -> Self {
        Self {
            rows,
            cols,
            data: vec![0.0; rows * cols],
            colnames: Vec::new(),
            rownames: Vec::new(),
        }
    }

    pub fn from_vec(rows: usize, cols: usize, data: Vec<f64>) -> Result<Self, MatrixError> {
        if data.len() != rows * cols {
            return Err(MatrixError::Shape {
                expected: rows * cols,
                got: data.len(),
            });
        }
        Ok(Self {
            rows,
            cols,
            data,
            colnames: Vec::new(),
            rownames: Vec::new(),
        })
    }

    /// Build from a list of **columns**, each of length `rows`.
    ///
    /// Distinct from [`Matrix::from_rows`], which takes a row-wise literal. The
    /// design matrix is naturally built one column at a time (a factor level is a
    /// whole column), so the two constructors earn their keep.
    pub fn from_cols(cols: &[Vec<f64>]) -> Self {
        let nc = cols.len();
        let nr = if nc == 0 { 0 } else { cols[0].len() };
        let mut data = vec![0.0; nr * nc];
        for (j, c) in cols.iter().enumerate() {
            assert_eq!(c.len(), nr, "all columns must have the same length");
            data[j * nr..(j + 1) * nr].copy_from_slice(c);
        }
        Self {
            rows: nr,
            cols: nc,
            data,
            colnames: Vec::new(),
            rownames: Vec::new(),
        }
    }

    /// R-style construction from a row-wise literal, for tests and fixtures.
    pub fn from_rows(rows: &[Vec<f64>]) -> Self {
        let nr = rows.len();
        let nc = if nr == 0 { 0 } else { rows[0].len() };
        let mut data = vec![0.0; nr * nc];
        for (i, r) in rows.iter().enumerate() {
            for (j, v) in r.iter().enumerate() {
                data[j * nr + i] = *v;
            }
        }
        Self {
            rows: nr,
            cols: nc,
            data,
            colnames: Vec::new(),
            rownames: Vec::new(),
        }
    }

    /// Column name `j`, or `col{j}` when the matrix carries no names.
    pub fn col_name(&self, j: usize) -> String {
        self.colnames
            .get(j)
            .cloned()
            .unwrap_or_else(|| format!("col{j}"))
    }

    pub fn with_colnames(mut self, names: Vec<String>) -> Self {
        self.colnames = names;
        self
    }

    pub fn with_rownames(mut self, names: Vec<String>) -> Self {
        self.rownames = names;
        self
    }

    #[inline]
    pub fn get(&self, i: usize, j: usize) -> f64 {
        self.data[j * self.rows + i]
    }

    #[inline]
    pub fn set(&mut self, i: usize, j: usize, v: f64) {
        let idx = j * self.rows + i;
        self.data[idx] = v;
    }

    #[inline]
    pub fn col(&self, j: usize) -> &[f64] {
        &self.data[j * self.rows..(j + 1) * self.rows]
    }

    #[inline]
    pub fn col_mut(&mut self, j: usize) -> &mut [f64] {
        let r = self.rows;
        &mut self.data[j * r..(j + 1) * r]
    }

    #[inline]
    pub fn row(&self, i: usize) -> impl Iterator<Item = f64> + '_ {
        (0..self.cols).map(move |j| self.get(i, j))
    }

    pub fn nrow(&self) -> usize {
        self.rows
    }

    pub fn ncol(&self) -> usize {
        self.cols
    }

    pub fn is_empty(&self) -> bool {
        self.rows == 0 || self.cols == 0
    }

    pub fn t(&self) -> Matrix {
        let mut out = Matrix::zeros(self.cols, self.rows);
        for j in 0..self.cols {
            for i in 0..self.rows {
                out.set(j, i, self.get(i, j));
            }
        }
        out
    }

    /// `self %*% other`
    pub fn matmul(&self, other: &Matrix) -> Result<Matrix, MatrixError> {
        if self.cols != other.rows {
            return Err(MatrixError::Shape {
                expected: self.cols,
                got: other.rows,
            });
        }
        let (n, k, m) = (self.rows, self.cols, other.cols);
        let mut out = Matrix::zeros(n, m);
        for j in 0..m {
            let oc = out.col_mut(j);
            for l in 0..k {
                let b = other.get(l, j);
                if b == 0.0 {
                    continue;
                }
                let sc = self.col(l);
                for i in 0..n {
                    oc[i] += sc[i] * b;
                }
            }
        }
        Ok(out)
    }

    pub fn add(&self, other: &Matrix) -> Result<Matrix, MatrixError> {
        if self.rows != other.rows || self.cols != other.cols {
            return Err(MatrixError::Shape {
                expected: self.rows * self.cols,
                got: other.rows * other.cols,
            });
        }
        let mut out = self.clone();
        for i in 0..out.data.len() {
            out.data[i] += other.data[i];
        }
        Ok(out)
    }

    pub fn sub(&self, other: &Matrix) -> Result<Matrix, MatrixError> {
        if self.rows != other.rows || self.cols != other.cols {
            return Err(MatrixError::Shape {
                expected: self.rows * self.cols,
                got: other.rows * other.cols,
            });
        }
        let mut out = self.clone();
        for i in 0..out.data.len() {
            out.data[i] -= other.data[i];
        }
        Ok(out)
    }

    pub fn scale(&self, s: f64) -> Matrix {
        let mut out = self.clone();
        for v in out.data.iter_mut() {
            *v *= s;
        }
        out
    }

    pub fn transpose_vec_rows(v: &[f64]) -> Matrix {
        Matrix::from_vec(v.len(), 1, v.to_vec()).expect("column vector")
    }

    /// Diagonal of the matrix, in R's order.
    pub fn diag(&self) -> Vec<f64> {
        (0..self.rows.min(self.cols))
            .map(|i| self.get(i, i))
            .collect()
    }

    pub fn set_diag(&mut self, d: &[f64]) {
        for (i, v) in d.iter().enumerate() {
            self.set(i, i, *v);
        }
    }

    /// Submatrix by column selection (R's `x[, idx]`), preserving the selected
    /// column names.
    pub fn select_cols(&self, idx: &[usize]) -> Matrix {
        let mut out = Matrix::zeros(self.rows, idx.len());
        for (k, &j) in idx.iter().enumerate() {
            let src = self.col(j);
            let dst = out.col_mut(k);
            dst.copy_from_slice(src);
        }
        out.colnames = idx.iter().map(|&j| self.col_name(j)).collect();
        out.rownames = self.rownames.clone();
        out
    }

    /// Row selection (R's `x[idx, ]`), preserving the selected row names.
    pub fn select_rows(&self, idx: &[usize]) -> Matrix {
        let mut out = Matrix::zeros(idx.len(), self.cols);
        for (k, &i) in idx.iter().enumerate() {
            for j in 0..self.cols {
                out.set(k, j, self.get(i, j));
            }
        }
        out.rownames = idx
            .iter()
            .map(|&i| {
                self.rownames
                    .get(i)
                    .cloned()
                    .unwrap_or_else(|| format!("row{i}"))
            })
            .collect();
        out.colnames = self.colnames.clone();
        out
    }

    /// [`Matrix::select_rows`] without the row names.
    ///
    /// The names cost one `String` per selected row, and the missingness-group
    /// designs are built once per group. A group can cover every sample, so on
    /// the 1000 x 10000 benchmark surface -- a distinct missingness pattern per
    /// taxon, which makes every taxon its own group -- this allocated one `String`
    /// per sample per group, 10 million in total, and that was the largest
    /// single source of allocation in the whole run.
    ///
    /// The grouped factorisation indexes `x_sub` by row number and never reads a
    /// name, so for this caller the names are pure cost. The data is identical to
    /// `select_rows`; only the labels differ.
    pub fn select_rows_anon(&self, idx: &[usize]) -> Matrix {
        let mut out = Matrix::zeros(idx.len(), self.cols);
        for (k, &i) in idx.iter().enumerate() {
            for j in 0..self.cols {
                out.set(k, j, self.get(i, j));
            }
        }
        out.colnames = self.colnames.clone();
        out
    }

    /// True when the matrix contains any NaN.
    pub fn has_nan(&self) -> bool {
        self.data.iter().any(|v| v.is_nan())
    }

    /// Symmetrise: `(A + A') / 2`
    pub fn symmetrise(&self) -> Matrix {
        let mut out = self.clone();
        for j in 0..self.cols {
            for i in 0..j {
                let v = 0.5 * (self.get(i, j) + self.get(j, i));
                out.set(i, j, v);
                out.set(j, i, v);
            }
        }
        out
    }
}

/// Operator sugar so the numeric code reads like the R it mirrors
/// (`x %*% y`, `x + y`, `x - y`). The fallible inherent methods
/// ([`Matrix::matmul`], [`Matrix::add`], [`Matrix::sub`]) remain the real API;
/// these panic on a shape mismatch, which is a programming error.
impl std::ops::Add<&Matrix> for &Matrix {
    type Output = Matrix;
    fn add(self, rhs: &Matrix) -> Matrix {
        Matrix::add(self, rhs).expect("shape mismatch in matrix addition")
    }
}

impl std::ops::Sub<&Matrix> for &Matrix {
    type Output = Matrix;
    fn sub(self, rhs: &Matrix) -> Matrix {
        Matrix::sub(self, rhs).expect("shape mismatch in matrix subtraction")
    }
}

impl std::ops::Mul<&Matrix> for &Matrix {
    type Output = Matrix;
    fn mul(self, rhs: &Matrix) -> Matrix {
        Matrix::matmul(self, rhs).expect("shape mismatch in matrix product")
    }
}

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum MatrixError {
    #[error("matrix shape mismatch: expected {expected} elements, got {got}")]
    Shape { expected: usize, got: usize },
}

// ---------------------------------------------------------------------------
// QR
// ---------------------------------------------------------------------------

/// Economy QR of an `n x p` matrix (`n >= p`), column-major Householder.
///
/// `lm.fit` uses LAPACK `dgeqp3` followed by `dormqr`; the *result* (the fitted
/// coefficients) depends only on the least-squares solution, so any accurate QR
/// agrees to rounding. This implementation is Householder with column pivoting
/// disabled, matching `lm.fit`'s no-pivoting default.
#[derive(Debug, Clone)]
pub struct Qr {
    /// `n x p` orthonormal columns of `Q`.
    ///
    /// `None`, and the reflectors applied to the right-hand sides instead, is the
    /// normal case. Materialising `Q` costs `O(n^2)` storage and time, and the
    /// only thing the pipeline ever asks of it is `Q'b` -- which is `O(n p)` per
    /// right-hand side when the reflectors are applied directly. For the design
    /// sub-matrices of a real table (`n` in the thousands, `p` under ten) that is
    /// the difference between a run that finishes and one that does not.
    pub q: Option<Matrix>,
    pub r: Matrix, // p x p, upper triangular
    pub rank: usize,
    /// Column permutation: `piv[k]` is the original column now in slot `k`.
    /// The identity, because the factorisation is unpivoted; the field exists so
    /// a future pivoted variant cannot silently change the meaning of `r`.
    pub piv: Vec<usize>,
    /// `true` for each column in `piv` order that is aliased, i.e. not among the
    /// first `rank`.
    pub aliased: Vec<bool>,
    /// The unit Householder reflectors `w_k`, so `Q = H_0 H_1 ... H_{p-1}` and
    /// `Q' b` is the reflectors applied to `b` in ascending `k`. Kept so the solve
    /// never needs `Q`.
    ///
    /// **Stored from row `k` onward**, not as a full length-`n` vector padded with
    /// zeros above row `k`. Every use reads `w[i]` only for `i` in `k..rows`, so the
    /// padding was `k` doubles per reflector that nothing could ever look at -- and
    /// this is the largest buffer in the whole analysis. The cache holds one
    /// reflector set per missingness pattern, so on `bm4` (1,000 patterns of one
    /// taxon each) the padding was about half a gigabyte, and on `bm5` (5,000
    /// patterns of ~20,000 rows) it is several gigabytes of the resident set. An
    /// all-zero reflector is still detected the same way, because the padding was
    /// zero too: `w.is_empty()` now stands in for `w.iter().all(|v| *v == 0.0)`.
    reflectors: Vec<Vec<f64>>,
}

impl Qr {
    /// Solve `min ||Ax - b||` for every column of `b` (`n x k`).
    ///
    /// Returns `None` if the design is rank deficient: the oracle detects that
    /// case (`fit$rank < ncol(xr)`) and refits taxon by taxon with `lm()`, whose
    /// behaviour is a per-taxon rank drop rather than a pseudo-solve, so the
    /// caller must be told to take a different path.
    /// Apply `f` to each response column, in parallel when the budget allows.
    ///
    /// `Matrix` is column-major, so `chunks_exact_mut(rows)` hands out one
    /// taxon's whole response column at a time and the chunks are disjoint by
    /// construction. That is what lets this be an ordinary safe `par_iter` over
    /// `&mut [f64]` rather than something needing `unsafe`: Rayon's indexed
    /// parallel iterator takes `Fn`, so a closure capturing `&mut Matrix` would not
    /// compile, and the column layout happens to make the safe form available.
    ///
    /// Determinism is unaffected. A column's dot product is a fixed ascending sum
    /// over its own rows and is never split, so no accumulation order is shared
    /// between workers.
    fn for_each_column(
        &self,
        lvl: &mut crate::parallel::Level,
        data: &mut [f64],
        rows: usize,
        f: impl Fn(&mut [f64]) + Send + Sync,
    ) {
        let n_resp = data.len() / rows.max(1);
        if lvl.may_parallelise() && n_resp > 1 {
            use rayon::prelude::*;
            data.par_chunks_mut(rows).for_each(f);
        } else {
            // A level with one column cannot split the pool, so it hands the pool
            // down to whatever is inside the column. Nothing inside a column is
            // split today; the release is what keeps a dense table from being
            // pinned to one core at the level above.
            lvl.release();
            for c in data.chunks_exact_mut(rows) {
                f(c);
            }
        }
    }

    /// [`Qr::solve_multi`] in place, so the caller does not pay for a second
    /// copy of the response block.
    ///
    /// `b` is `n_used x n_responses` and is *overwritten*: the solution is
    /// written into the first `p` rows, un-permuted, and the caller reads it from
    /// there. The rows below `p` are left holding garbage and must not be read.
    ///
    /// # Why this exists
    ///
    /// The block a caller assembles for a missingness group is
    /// `n_used x taxa_in_group`, and on the 1000 x 10000 benchmark surface that
    /// is hundreds of kilobytes *per group, per MLE iteration*. `solve_multi`
    /// used to take `&Matrix` and clone it internally, so the block was built
    /// once and copied once more -- 5.4 GB of allocation churn on that dataset,
    /// and one more live buffer at the peak. The copy exists only because the
    /// function is written against an immutable interface; taking `&mut` removes
    /// it, and the arithmetic is unchanged.
    pub fn solve_multi_into(&self, b: &mut Matrix) -> Option<()> {
        if self.rank < self.r.cols {
            return None;
        }
        let n_resp = b.cols;
        let p = self.r.rows;
        debug_assert!(b.rows >= p);
        // Apply the reflectors to `b` itself, in ascending k.
        //
        // The inner loop is over *right-hand sides*, which is one taxon's
        // response column apiece -- the innermost level of the nesting order. The
        // columns are disjoint, so splitting them cannot change a bit of the
        // result, and the reduction inside a column (the dot product) stays a
        // fixed ascending sum because it is not itself split.
        //
        // This is the level that matters on a table with no missing values: there
        // is then a single missingness pattern, the group level above has one item
        // and hands the pool down, and every expensive operation in the MLE is
        // this loop.
        let mut lvl = crate::parallel::NestingBudget::level("taxa (response columns)");
        let rows = b.rows;
        for (k, w) in self.reflectors.iter().enumerate() {
            if w.iter().all(|v| *v == 0.0) {
                continue;
            }
            self.for_each_column(&mut lvl, &mut b.data, rows, |col| {
                // `w` holds rows `k..n`, so its index is `i - k`.
                let mut dot = 0.0;
                for i in k..rows {
                    dot += w[i - k] * col[i];
                }
                let f = 2.0 * dot;
                if f == 0.0 {
                    return;
                }
                for i in k..rows {
                    col[i] -= f * w[i - k];
                }
            });
        }
        // Back-substitute in place: row `j` becomes the solved value, then the
        // rows above it are corrected. The diagonal index runs over the
        // *parameter* axis and the response axis is separate; conflating them
        // breaks for a single response.
        // The `j` loop stays sequential: row `j` of the solution depends on row
        // `j + 1`, which the iteration below has already finished. The `r` loop
        // inside each is a whole column, so it is the same split as above.
        for j in (0..self.rank).rev() {
            let rjj = self.r.get(j, j);
            self.for_each_column(&mut lvl, &mut b.data, rows, |col| {
                col[j] /= rjj;
            });
            for i in 0..j {
                let rij = self.r.get(i, j);
                if rij == 0.0 {
                    continue;
                }
                // `i < j`, so `col[j]` is read before `col[i]` is written and one
                // column slice is all the state this needs.
                self.for_each_column(&mut lvl, &mut b.data, rows, |col| {
                    col[i] -= rij * col[j];
                });
            }
        }
        // Un-permute the first `p` rows in place. The solution is currently in
        // pivoted order in rows 0..p; a cycle needs a temporary, so the rows are
        // rotated through the `x` scratch rather than permuted by an
        // allocation.
        // Un-permuting needs a temporary: the rows are in pivoted order, and
        // moving them in place would overwrite a value still to be read. The
        // scratch is one column of `p` doubles, not the whole block.
        let mut col = vec![0.0f64; p];
        for r in 0..n_resp {
            for k in 0..p {
                col[k] = b.get(k, r);
            }
            for k in 0..p {
                b.set(self.piv[k], r, col[k]);
            }
        }
        Some(())
    }

    pub fn solve_multi(&self, b: &Matrix) -> Option<Matrix> {
        let x = self.solve_pivoted(b)?;
        let p = self.r.cols;
        let mut out = Matrix::zeros(p, b.cols);
        for k in 0..p {
            for r in 0..b.cols {
                out.set(self.piv[k], r, x.get(k, r));
            }
        }
        Some(out)
    }

    /// The least-squares solution in *pivoted* order, which is the order the
    /// back-substitution naturally produces.
    ///
    /// Returns `None` for a rank-deficient design: that case needs a decision
    /// the caller makes, not a pseudo-solve. Use [`Qr::solve_multi_padded`] for
    /// the `lm`-style "drop the aliased columns" answer.
    pub fn solve_pivoted(&self, b: &Matrix) -> Option<Matrix> {
        debug_assert!(
            b.rows >= self.r.cols,
            "b has fewer rows than the design has columns"
        );
        if self.rank < self.r.cols {
            return None;
        }
        Some(self.solve_pivoted_partial(b, self.rank))
    }

    /// Back-substitute using only the first `rank` diagonals of `R`, so a
    /// rank-deficient design still yields the least-squares solution over the
    /// well-determined columns. Rows `rank..` are left untouched.
    fn solve_pivoted_partial(&self, b: &Matrix, rank: usize) -> Matrix {
        let n_resp = b.cols;
        let p = self.r.rows;
        // Q'b is p x k. Rather than forming Q, the reflectors are applied to b
        // directly: Q = H_0 H_1 ... H_{p-1} and each H_k is symmetric, so
        // Q'b = H_{p-1} ... H_0 b, which is the reflectors applied in *ascending*
        // k. The result is identical to the explicit product and costs
        // O(n p) per right-hand side instead of O(n^2) to build Q.
        let mut qtb_full = b.clone();
        for (k, w) in self.reflectors.iter().enumerate() {
            // `w` holds rows `k..n`; the all-zero test is over the whole tail.
            if w.iter().all(|v| *v == 0.0) {
                continue;
            }
            for r in 0..n_resp {
                let mut dot = 0.0;
                for i in k..b.rows {
                    dot += w[i - k] * qtb_full.get(i, r);
                }
                let f = 2.0 * dot;
                if f == 0.0 {
                    continue;
                }
                for i in k..b.rows {
                    let cur = qtb_full.get(i, r);
                    qtb_full.set(i, r, cur - f * w[i - k]);
                }
            }
        }
        let mut qtb = Matrix::zeros(p, n_resp);
        // `Q'b` has only `b.rows` rows even when `p > b.rows`: `Q` is
        // `n x min(n, p)`, so there is no row `l >= n` to read. The remaining
        // rows stay zero, which is also what `lm.fit` leaves in the trailing
        // coordinates of a rank-deficient solve.
        for l in 0..p.min(b.rows) {
            for r in 0..n_resp {
                qtb.set(l, r, qtb_full.get(l, r));
            }
        }
        // Back-substitute R x = Q'b, column by column. The solution is p x k,
        // so the diagonal index runs over the *parameter* axis and the response
        // axis is separate; conflating them breaks for a single response.
        let mut x = qtb;
        for j in (0..rank).rev() {
            let rjj = self.r.get(j, j);
            for r in 0..n_resp {
                let v = x.get(j, r) / rjj;
                x.set(j, r, v);
            }
            for i in 0..j {
                let rij = self.r.get(i, j);
                if rij == 0.0 {
                    continue;
                }
                for r in 0..n_resp {
                    let v = x.get(i, r) - rij * x.get(j, r);
                    x.set(i, r, v);
                }
            }
        }
        x
    }

    /// Bytes held by the stored reflectors, for the memory report.
    ///
    /// The tails sum to `p * n - p * (p - 1) / 2` doubles rather than `p * n`,
    /// which is the saving the truncation buys.
    pub fn reflector_bytes(&self) -> usize {
        self.reflectors.iter().map(|w| w.len() * 8).sum()
    }

    pub fn solve(&self, b: &[f64]) -> Option<Vec<f64>> {
        let bm = Matrix::from_vec(b.len(), 1, b.to_vec()).ok()?;
        self.solve_multi(&bm).map(|m| m.data)
    }

    /// The least-squares solution for every column of `b`, with an aliased column
    /// left at **zero** rather than omitted.
    ///
    /// # The `0.0` here is known to be wrong, and is kept deliberately
    ///
    /// `coef.lm` does **not** drop an aliased term: it keeps the name and reports
    /// `NA`. Verified on this oracle's own R --
    ///
    /// ```r
    /// d$b <- 2 * d$a + rnorm(n, sd = 1e-10)
    /// coef(lm(y ~ a + b, d))   #  a: -0.228,  b: NA
    /// ```
    ///
    /// so `.lm_fit_all`'s `bi[match(names(ci), fix_eff)] = ci` writes `NA` into
    /// that slot, and `.bias_em`'s `neither_na` then *drops the taxon* from the
    /// bias fit. Writing `0.0` here keeps it, with a fabricated coefficient.
    ///
    /// It was changed to `f64::NAN` and measured, and the result was **worse**:
    /// on the `int-sparsity90-5group` fixture-matrix cell the worst E-M
    /// component-weight deviation went from 4.0e-1 to 1.0e+1. That is diagnostic
    /// rather than discouraging -- if the two sides then drop the *same* taxa the
    /// mixtures should agree, so they are not dropping the same set. The remaining
    /// defect is therefore in *which* columns are declared aliased, not in what is
    /// written for one: R reaches it through `dqrls`' complete orthogonal
    /// factorisation with column pivoting, and this pivoting differs on the
    /// near-singular sub-designs that 90%-zero tables produce.
    ///
    /// So the two have to be fixed together, and until the aliased *set* matches,
    /// changing only the fill value moves the answer further from the oracle while
    /// still being semantically closer to it. A semantically-closer-but-measurably-
    /// worse intermediate state is not worth committing, so `0.0` stands, recorded
    /// here and in `docs/reference_behavior.md` section 16 as a known divergence
    /// with its mechanism, its measured scope, and what has already been ruled
    /// out.
    ///
    /// `None` is still returned for a rank-deficient design, because the grouped
    /// path must then refit per taxon -- a *separate* `lm` per taxon, each of
    /// which may have a different rank.
    pub fn solve_multi_padded(&self, b: &Matrix) -> Option<Matrix> {
        if self.rank == 0 {
            return None;
        }
        let x = self.solve_pivoted_partial(b, self.rank);
        let p = self.r.cols;
        let mut out = Matrix::zeros(p, b.cols);
        for k in 0..p {
            let orig = self.piv[k];
            for r in 0..b.cols {
                out.set(
                    orig,
                    r,
                    if self.aliased[k] {
                        f64::NAN
                    } else {
                        x.get(k, r)
                    },
                );
            }
        }
        Some(out)
    }

    /// `solve_multi_padded` for a single response, in original column order.
    pub fn solve_padded(&self, b: &[f64]) -> Option<Vec<f64>> {
        let bm = Matrix::from_vec(b.len(), 1, b.to_vec()).ok()?;
        self.solve_multi_padded(&bm).map(|m| m.data)
    }
}

/// `stats::lm.fit`'s documented default `tol`, which is the value `.lm_fit_all`
/// gets: it calls `lm.fit(xr, Yr)` and `lm(y_crt ~ .)` without a `tol` argument.
/// It is *not* LAPACK's `dlamch("epsilon")` -- that is the `tol = -1` sentinel,
/// and using it changes which columns count as aliased on an ill-conditioned
/// design.
pub const LM_FIT_TOL: f64 = 1e-7;

/// Householder QR with a rank-revealing tolerance matching `lm.fit`.
///
/// Returns the economy form: `Q` is `n x p` with orthonormal columns and `R` is
/// the `p x p` upper triangle. `lm.fit` uses LAPACK `dgeqp3` + `dormqr`; the
/// fitted coefficients depend only on the least-squares solution, so any
/// accurate QR agrees to rounding, and pivoting is off in both.
pub fn qr(a: &Matrix) -> Qr {
    let n = a.rows;
    let p = a.cols;
    // `n < p` is legal and has to be handled rather than asserted: a missingness
    // group can have fewer observed samples than the design has columns, and the
    // reference handles it -- `lm.fit` on an `n x p` matrix with `n < p` returns
    // `rank <= n < p`, and `.lm_fit_all` then refits each taxon on its own. This
    // used to `assert!(n >= p)` and abort the whole run on such a group; the
    // fixture matrix's `int-sparsity90-5group` cell has six design columns and
    // 90% zeros, so several of its groups are wider than they are tall.
    //
    // A rank-`n` factorisation still exists for such a matrix, so the reflector
    // loop simply stops at `min(n, p)`: the extra columns are unrepresented in
    // `R`, and the rank count below is bounded by `n` so they cannot be
    // reported as estimated.
    let kmax = n.min(p);
    let mut r = a.clone();
    let mut v: Vec<f64> = Vec::new();
    let mut piv: Vec<usize> = (0..p).collect();
    let mut householder: Vec<Vec<f64>> = (0..kmax).map(|k| vec![0.0f64; n - k]).collect();
    for k in 0..kmax {
        // --- no column permutation: this is `dqrls`' choice, not ours ---
        //
        // `stats::lm`'s `lm.fit` does **not** select columns. Measured on this
        // oracle's R:
        //
        // ```text
        // > lm.fit(cbind(c1, c2, c1 + c2), y)$qr$pivot
        // [1] 1 2 3
        // > coef(lm.fit(cbind(c1, c2, c1 + c2), y))
        // [1] -0.2857365  0.7701500          NA
        // ```
        //
        // The pivot is the identity even when a column is an exact combination of
        // the ones before it, and the column that comes out `NA` is the *last* one
        // in build order, not the best-conditioned one. An earlier revision of this
        // function pivoted on the largest remaining 2-norm, which is the
        // better-conditioned choice and the wrong one.
        //
        // What the rank test below then has to agree with is `dqrls`' notion of
        // rank, which is the length of the **leading run** of diagonals above
        // `tol * |R(0,0)|` -- it stops at the first diagonal that is not, rather
        // than counting every diagonal that is. That is why `aliased` below can be
        // written as `i >= rank`: the leading run is by construction a prefix, so
        // "at or after the rank" and "not in the leading run" are the same set.
        //
        // The consequence is accepted rather than worked around: an unpivoted QR
        // can report a small diagonal for a full-rank but ill-conditioned column
        // and under-state the rank. R has exactly that property -- it is the
        // behaviour being matched, not a hazard introduced here.
        let _ = &mut piv;

        // Everything below works on raw slices rather than `Matrix::get`/`set`.
        // `r` is column-major (`data[j * n + i]`), so the rows `k..n` of column
        // `j` are one contiguous run -- but `get(i, j)` hides that behind a
        // multiply and a bounds check per element, and the compiler will not
        // vectorise a loop it cannot prove is contiguous. Measured on the bm6
        // shape: 1.24 -> 2.4 GFLOP/s for the same arithmetic.
        //
        // The arithmetic is unchanged in every particular that could move a bit:
        // same operands, same order, same operations. `vnorm_sq` now sums
        // `v[k..]` rather than all of `v`, which is the same number because the
        // skipped entries are exactly `0.0` and `0.0 + x == x` for every finite
        // `x`, so the partial sums before the first real term are identical.
        let col_k = &r.data[k * n + k..(k + 1) * n];
        let norm_sq: f64 = col_k.iter().map(|v| v * v).sum();
        let norm = norm_sq.sqrt();
        if norm == 0.0 {
            continue;
        }
        // Choose the sign that avoids cancellation, as LAPACK does.
        let alpha = if r.get(k, k) >= 0.0 { -norm } else { norm };
        // One scratch buffer for the whole factorisation, not one per reflector:
        // `n` doubles were being allocated and zero-filled p times per QR.
        if v.len() < n {
            v.resize(n, 0.0);
        }
        v[..=k].fill(0.0);
        v[k..n].copy_from_slice(col_k);
        v[k] -= alpha;
        let vk = &v[k..n];
        let vnorm_sq: f64 = vk.iter().map(|x| x * x).sum();
        if vnorm_sq <= 0.0 {
            continue;
        }
        for j in k..p {
            let col = &mut r.data[j * n + k..(j + 1) * n];
            let dot: f64 = vk.iter().zip(col.iter()).map(|(a, b)| a * b).sum();
            let f = 2.0 * dot / vnorm_sq;
            for (a, b) in vk.iter().zip(col.iter_mut()) {
                *b -= f * a;
            }
        }
        // Store the unit reflector w = v / ||v|| so the Q accumulation below
        // applies H = I - 2 w w'. A division per element, not a multiply by the
        // reciprocal: `x / wnorm` and `x * (1 / wnorm)` differ in the last bit,
        // and this value is applied to every right-hand side.
        let wnorm = vnorm_sq.sqrt();
        householder[k] = vk.iter().map(|x| x / wnorm).collect();
        r.set(k, k, alpha);
    }

    // Q is *not* formed. `solve_pivoted` applies the reflectors to the right-hand
    // sides instead, which is the same arithmetic without the O(n^2) matrix.
    // `explicit_q()` builds it, for the tests that check the solve against a
    // direct Gram-Schmidt.
    // R is in *pivoted* order: R[k][j] belongs to original column piv[j].
    let mut rsmall = Matrix::zeros(p, p);
    for j in 0..p {
        // `i <= j` *and* `i < n`: with `n < p` the reflector loop only touched
        // the first `n` rows, so the lower triangle of `rsmall` below row `n` is
        // zero -- which is what makes the rank count stop at `kmax`.
        for i in 0..=j.min(n - 1) {
            rsmall.set(i, j, r.get(i, j));
        }
    }

    // `dqrls`'s rank test: a diagonal counts when it exceeds
    // `tol * |R(0,0)|`. Note there is no `max(m, n)` factor here -- that belongs
    // to `dgeqrf`, and using it would keep columns that `lm.fit` drops, turning
    // a rank-deficient group into a fit with 1e14 coefficients.
    let scale = rsmall.get(0, 0).abs().max(1e-300);
    let tol = LM_FIT_TOL * scale;
    // LAPACK counts the *leading* run of diagonals above `tol` and stops at the
    // first one that is not: a later diagonal can be large again after an exact
    // zero, and counting it would report a full rank for a design whose middle
    // column is an exact duplicate of an earlier one.
    let mut rank = 0;
    for i in 0..kmax {
        if rsmall.get(i, i).abs() > tol {
            rank = i + 1;
        } else {
            break;
        }
    }
    let aliased: Vec<bool> = (0..p).map(|i| i >= rank).collect();

    Qr {
        q: None,
        r: rsmall,
        rank,
        piv,
        aliased,
        reflectors: householder,
    }
}

impl Qr {
    /// Form `Q` explicitly, `n x p`.
    ///
    /// Only the tests need this. It is `O(n^2 p)` and allocates `n^2` doubles, so
    /// calling it in a hot path reintroduces the cost the solve avoids.
    pub fn explicit_q(&self) -> Matrix {
        // `reflectors[0]` holds rows `0..n`, so its length *is* `n`. The later
        // ones are tails, and are indexed from their own `k`.
        let n = self.reflectors.first().map_or(0, |w| w.len());
        let p = self.r.cols;
        let mut qfull = Matrix::zeros(n, n);
        for i in 0..n {
            qfull.set(i, i, 1.0);
        }
        for k in (0..p).rev() {
            let Some(w) = self.reflectors.get(k) else {
                continue;
            };
            for j in 0..n {
                let dot: f64 = w
                    .iter()
                    .enumerate()
                    .map(|(t, wi)| wi * qfull.get(k + t, j))
                    .sum();
                let f = 2.0 * dot;
                for (t, wi) in w.iter().enumerate() {
                    let i = k + t;
                    let cur = qfull.get(i, j);
                    qfull.set(i, j, cur - f * wi);
                }
            }
        }
        let mut q = Matrix::zeros(n, p);
        for j in 0..p {
            for i in 0..n {
                q.set(i, j, qfull.get(i, j));
            }
        }
        q
    }
}

// ---------------------------------------------------------------------------
// linear solve
// ---------------------------------------------------------------------------

/// Solve `A x = b` for a symmetric positive definite `A` by Cholesky.
pub fn cholesky_solve(a: &Matrix, b: &[f64]) -> Result<Vec<f64>, MatrixError> {
    let p = a.rows;
    let mut l = Matrix::zeros(p, p);
    for i in 0..p {
        for j in 0..=i {
            let mut s = a.get(i, j);
            for k in 0..j {
                s -= l.get(i, k) * l.get(j, k);
            }
            if i == j {
                if s <= 0.0 {
                    return Err(MatrixError::Shape {
                        expected: 0,
                        got: 1,
                    });
                }
                l.set(i, j, s.sqrt());
            } else {
                l.set(i, j, s / l.get(j, j));
            }
        }
    }
    // forward then back substitution
    let mut y = vec![0.0; p];
    for i in 0..p {
        let mut s = b[i];
        for k in 0..i {
            s -= l.get(i, k) * y[k];
        }
        y[i] = s / l.get(i, i);
    }
    let mut x = vec![0.0; p];
    for i in (0..p).rev() {
        let mut s = y[i];
        for k in (i + 1)..p {
            s -= l.get(k, i) * x[k];
        }
        x[i] = s / l.get(i, i);
    }
    Ok(x)
}

// ---------------------------------------------------------------------------
// eigen decomposition (symmetric, cyclic Jacobi)
// ---------------------------------------------------------------------------

/// Symmetric eigen-decomposition. Returns `(eigenvalues ascending, eigenvectors)`
/// with eigenvectors as columns.
pub fn eigen_symmetric(a: &Matrix) -> (Vec<f64>, Matrix) {
    let n = a.rows;
    let mut m = a.symmetrise();
    // V accumulates the rotations
    let mut v = Matrix::zeros(n, n);
    for i in 0..n {
        v.set(i, i, 1.0);
    }
    for _sweep in 0..100 {
        let mut off = 0.0;
        for i in 0..n {
            for j in (i + 1)..n {
                off += m.get(i, j) * m.get(i, j);
            }
        }
        if off <= 1e-30 {
            break;
        }
        for p in 0..n {
            for q in (p + 1)..n {
                let apq = m.get(p, q);
                if apq.abs() < 1e-300 {
                    continue;
                }
                let theta = (m.get(q, q) - m.get(p, p)) / (2.0 * apq);
                let t = if theta >= 0.0 {
                    1.0 / (theta + (1.0 + theta * theta).sqrt())
                } else {
                    -1.0 / (-theta + (1.0 + theta * theta).sqrt())
                };
                let c = 1.0 / (1.0 + t * t).sqrt();
                let s = t * c;
                for k in 0..n {
                    let akp = m.get(k, p);
                    let akq = m.get(k, q);
                    m.set(k, p, c * akp - s * akq);
                    m.set(k, q, s * akp + c * akq);
                }
                for k in 0..n {
                    let apk = m.get(p, k);
                    let aqk = m.get(q, k);
                    m.set(p, k, c * apk - s * aqk);
                    m.set(q, k, s * apk + c * aqk);
                }
                for k in 0..n {
                    let vkp = v.get(k, p);
                    let vkq = v.get(k, q);
                    v.set(k, p, c * vkp - s * vkq);
                    v.set(k, q, s * vkp + c * vkq);
                }
            }
        }
    }
    let mut pairs: Vec<(f64, usize)> = (0..n).map(|i| (m.get(i, i), i)).collect();
    pairs.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
    let vals: Vec<f64> = pairs.iter().map(|p| p.0).collect();
    let mut vsorted = Matrix::zeros(n, n);
    for (k, &(_, orig)) in pairs.iter().enumerate() {
        for i in 0..n {
            vsorted.set(i, k, v.get(i, orig));
        }
    }
    (vals, vsorted)
}

/// `MASS::ginv(x, tol = sqrt(.Machine$double.eps))`.
///
/// MASS computes an SVD and keeps the singular values above
/// `max(tol * d[0], 0)`. For a symmetric input the singular values are the
/// absolute eigenvalues, so an eigen-decomposition with the same relative cutoff
/// is numerically equivalent and avoids an SVD.
pub fn ginv(a: &Matrix, tol: Option<f64>) -> Matrix {
    let n = a.rows;
    if n == 0 {
        return Matrix::zeros(0, 0);
    }
    let tol = tol.unwrap_or(f64::EPSILON.sqrt());
    let s = a.symmetrise();
    let (vals, vecs) = eigen_symmetric(&s);
    let dmax = vals.iter().fold(0.0f64, |m, v| m.max(v.abs()));
    let thresh = (tol * dmax).max(0.0);
    let mut out = Matrix::zeros(n, n);
    let positive: Vec<usize> = (0..n).filter(|&k| vals[k].abs() > thresh).collect();
    if positive.is_empty() {
        return out;
    }
    for &k in &positive {
        let inv = 1.0 / vals[k];
        for i in 0..n {
            for j in 0..n {
                out.set(i, j, out.get(i, j) + vecs.get(i, k) * inv * vecs.get(j, k));
            }
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn approx(a: &Matrix, b: &Matrix, tol: f64) {
        assert_eq!(a.rows, b.rows);
        assert_eq!(a.cols, b.cols);
        for i in 0..a.data.len() {
            assert!(
                (a.data[i] - b.data[i]).abs() <= tol * (1.0 + b.data[i].abs()),
                "element {i}: {} vs {}",
                a.data[i],
                b.data[i]
            );
        }
    }

    #[test]
    fn matmul_matches_hand_computation() {
        let a = Matrix::from_rows(&[vec![1.0, 2.0], vec![3.0, 4.0]]);
        let b = Matrix::from_rows(&[vec![5.0, 6.0], vec![7.0, 8.0]]);
        let c = (&a * &b).clone();
        // R: matrix(c(1,2,3,4),2,2) %*% matrix(c(5,6,7,8),2,2) = [19,22; 43,50]
        assert_eq!(c.rows, 2);
        assert_eq!(c.cols, 2);
        assert!((c.get(0, 0) - 19.0).abs() < 1e-14);
        assert!((c.get(1, 0) - 43.0).abs() < 1e-14);
        assert!((c.get(0, 1) - 22.0).abs() < 1e-14);
        assert!((c.get(1, 1) - 50.0).abs() < 1e-14);
    }

    #[test]
    fn qr_solves_a_least_squares_problem() {
        // R: coef(lm(c(1,3,2,4,5) ~ c(0,1,2,3,4))) == c(1.2, 0.9)
        let x = Matrix::from_rows(&[
            vec![1.0, 0.0],
            vec![1.0, 1.0],
            vec![1.0, 2.0],
            vec![1.0, 3.0],
            vec![1.0, 4.0],
        ]);
        let y = vec![1.0, 3.0, 2.0, 4.0, 5.0];
        let q = qr(&x);
        assert_eq!(q.rank, 2);
        // The solve applies the reflectors to the right-hand sides instead of
        // forming Q. That is only valid if the two agree, which is the invariant
        // this test exists to protect: it is a 67x speed-up on a real table and a
        // silent wrong answer if the application order is wrong.
        let p = x.cols;
        let qm = q.explicit_q();
        let gram = qm.t().matmul(&qm).expect("Q'Q");
        for i in 0..p {
            for j in 0..p {
                let want = if i == j { 1.0 } else { 0.0 };
                assert!(
                    (gram.get(i, j) - want).abs() < 1e-12,
                    "Q must be orthonormal, but Q'Q[{i},{j}] = {}",
                    gram.get(i, j)
                );
            }
        }
        let b = q.solve(&y).unwrap();
        assert!((b[0] - 1.2).abs() < 1e-12, "intercept {}", b[0]);
        assert!((b[1] - 0.9).abs() < 1e-12, "slope {}", b[1]);
    }

    #[test]
    fn qr_handles_multiple_responses_at_once() {
        // R: lm.fit(x, cbind(y1, y2)) where y2 = 2*y1
        let x = Matrix::from_rows(&[
            vec![1.0, 0.0],
            vec![1.0, 1.0],
            vec![1.0, 2.0],
            vec![1.0, 3.0],
        ]);
        let y1 = [1.0, 3.0, 2.0, 4.0];
        let mut b = Matrix::zeros(4, 2);
        for i in 0..4 {
            b.set(i, 0, y1[i]);
            b.set(i, 1, 2.0 * y1[i]);
        }
        let q = qr(&x);
        // R: coef(lm(c(1,3,2,4) ~ c(0,1,2,3))) == c(1.3, 0.8), and y2 = 2*y1 so
        // the second response's coefficients double.
        let sol = q.solve_multi(&b).unwrap();
        assert!((sol.get(0, 0) - 1.3).abs() < 1e-12, "{}", sol.get(0, 0));
        assert!((sol.get(1, 0) - 0.8).abs() < 1e-12, "{}", sol.get(1, 0));
        assert!((sol.get(0, 1) - 2.6).abs() < 1e-12, "{}", sol.get(0, 1));
        assert!((sol.get(1, 1) - 1.6).abs() < 1e-12, "{}", sol.get(1, 1));
    }

    #[test]
    fn qr_detects_rank_deficiency() {
        // duplicated column: rank 1
        let x = Matrix::from_rows(&[vec![1.0, 2.0], vec![2.0, 4.0], vec![3.0, 6.0]]);
        let q = qr(&x);
        assert_eq!(q.rank, 1, "duplicated columns must be rank 1");
        assert!(q.solve(&[1.0, 2.0, 3.0]).is_none());
    }

    #[test]
    fn ginv_inverts_a_well_conditioned_matrix() {
        // R: ginv(diag(1, 2, 3)) == diag(1, 1/2, 1/3)
        let a = Matrix::from_rows(&[
            vec![1.0, 0.0, 0.0],
            vec![0.0, 2.0, 0.0],
            vec![0.0, 0.0, 3.0],
        ]);
        let g = ginv(&a, None);
        let expect = Matrix::from_rows(&[
            vec![1.0, 0.0, 0.0],
            vec![0.0, 0.5, 0.0],
            vec![0.0, 0.0, 1.0 / 3.0],
        ]);
        approx(&g, &expect, 1e-12);
        // and A %*% ginv(A) is the identity
        approx(
            &a.matmul(&g).unwrap(),
            &Matrix::from_rows(&[
                vec![1.0, 0.0, 0.0],
                vec![0.0, 1.0, 0.0],
                vec![0.0, 0.0, 1.0],
            ]),
            1e-12,
        );
    }

    #[test]
    fn ginv_drops_tiny_singular_values_like_mass() {
        // R:
        //   diag(1, 1e-12) ; ginv()  -> diag(1, 0)   because sqrt(eps) = 1.49e-8
        //   diag(1, 1e-4)  ; ginv()  -> diag(1, 1e4)
        let tiny = Matrix::from_rows(&[vec![1.0, 0.0], vec![0.0, 1e-12]]);
        let g = ginv(&tiny, None);
        assert!((g.get(0, 0) - 1.0).abs() < 1e-12);
        assert_eq!(g.get(1, 1), 0.0, "1e-12 < sqrt(eps), must be dropped");

        let small = Matrix::from_rows(&[vec![1.0, 0.0], vec![0.0, 1e-4]]);
        let g = ginv(&small, None);
        assert!((g.get(1, 1) - 1e4).abs() < 1e-6, "{}", g.get(1, 1));
    }

    #[test]
    fn ginv_of_a_design_gram_matrix() {
        // The oracle's exact call: ginv(t(x) %*% x) for a 1-intercept design
        let x = Matrix::from_rows(&[
            vec![1.0, 0.0],
            vec![1.0, 1.0],
            vec![1.0, 2.0],
            vec![1.0, 3.0],
            vec![1.0, 4.0],
        ]);
        let xtx = x.t().matmul(&x).unwrap();
        let g = ginv(&xtx, None);
        // R: MASS::ginv(t(x) %*% x) for x = cbind(1, 0:4) is
        //   [ 0.6 -0.2 ]
        //   [-0.2  0.1 ]
        assert!((g.get(0, 0) - 0.6).abs() < 1e-10, "{}", g.get(0, 0));
        assert!((g.get(0, 1) + 0.2).abs() < 1e-10);
        assert!((g.get(1, 1) - 0.1).abs() < 1e-10);
    }

    #[test]
    fn eigen_matches_r_symmetric() {
        // R: eigen(diag(3,1,2))$values == c(1, 2, 3)
        let a = Matrix::from_rows(&[
            vec![3.0, 0.0, 0.0],
            vec![0.0, 1.0, 0.0],
            vec![0.0, 0.0, 2.0],
        ]);
        let (vals, _) = eigen_symmetric(&a);
        assert!((vals[0] - 1.0).abs() < 1e-12);
        assert!((vals[1] - 2.0).abs() < 1e-12);
        assert!((vals[2] - 3.0).abs() < 1e-12);
    }

    #[test]
    fn eigen_of_a_general_symmetric_matrix() {
        // R: eigen(matrix(c(2,-1,0,-1,2,-1,0,-1,2), 3, 3))$values
        //   == c(2 - sqrt(2), 2, 2 + sqrt(2))
        let a = Matrix::from_rows(&[
            vec![2.0, -1.0, 0.0],
            vec![-1.0, 2.0, -1.0],
            vec![0.0, -1.0, 2.0],
        ]);
        let (vals, vecs) = eigen_symmetric(&a);
        let s2 = std::f64::consts::SQRT_2;
        assert!((vals[0] - (2.0 - s2)).abs() < 1e-12);
        assert!((vals[1] - 2.0).abs() < 1e-12);
        assert!((vals[2] - (2.0 + s2)).abs() < 1e-12);
        // reconstruct
        for i in 0..3 {
            for j in 0..3 {
                let mut s = 0.0;
                for k in 0..3 {
                    s += vecs.get(i, k) * vals[k] * vecs.get(j, k);
                }
                assert!((s - a.get(i, j)).abs() < 1e-12, "reconstruct ({i},{j})");
            }
        }
    }

    #[test]
    fn cholesky_solve_matches_r_solve() {
        // R: solve(matrix(c(4,2,2,3), 2, 2), c(1, 2)) == c(-0.125, 0.75)
        //   (matrix() is column-major, so the system is [[4,2],[2,3]] x = [1,2])
        let a = Matrix::from_rows(&[vec![4.0, 2.0], vec![2.0, 3.0]]);
        let x = cholesky_solve(&a, &[1.0, 2.0]).unwrap();
        assert!((x[0] + 0.125).abs() < 1e-14, "{}", x[0]);
        assert!((x[1] - 0.75).abs() < 1e-14, "{}", x[1]);
    }

    #[test]
    fn column_major_layout_is_r_compatible() {
        // R: as.vector(matrix(1:6, 2, 3)) is 1, 2, 3, 4, 5, 6
        let m = Matrix::from_rows(&[vec![1.0, 3.0, 5.0], vec![2.0, 4.0, 6.0]]);
        assert_eq!(m.data, vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0]);
    }

    #[test]
    fn select_helpers_match_r_subsetting() {
        let m = Matrix::from_rows(&[
            vec![1.0, 2.0, 3.0],
            vec![4.0, 5.0, 6.0],
            vec![7.0, 8.0, 9.0],
        ]);
        let c = m.select_cols(&[2, 0]);
        assert_eq!(c.data, vec![3.0, 6.0, 9.0, 1.0, 4.0, 7.0]);
        let r = m.select_rows(&[2, 0]);
        assert_eq!(r.rows, 2);
        assert!((r.get(0, 2) - 9.0).abs() < 1e-15);
    }
}
