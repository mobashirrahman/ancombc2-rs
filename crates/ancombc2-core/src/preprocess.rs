//! Preprocessing: prevalence and library-size filtering, the log transform,
//! per-taxon centring, and structural-zero detection.
//!
//! Each function corresponds to a specific step of the reference and preserves
//! its exact semantics, including the places where those semantics are
//! surprising:
//!
//! * prevalence is `rowSums(x != 0) / rowSums(!is.na(x))` — *observed* samples in
//!   the denominator, so an `NA` is neither a presence nor an absence;
//! * the cutoff is `>=`, not `>`;
//! * `log(0)` is `NA` (not `-Inf`) and the row mean is taken over the finite
//!   entries only.

use crate::error::{AncombcError, Result};
use crate::workspace::RMatrix;

/// Row-major `taxa x samples` count matrix.
#[derive(Debug, Clone, Default)]
pub struct CountMatrix {
    pub n_taxa: usize,
    pub n_samp: usize,
    /// Row-major. `NaN` marks a missing count.
    pub data: Vec<f64>,
}

/// The default is the empty table: zero taxa, zero samples, no data.
///
/// A `CountMatrix` is a shape plus a buffer, so `Default` has only one sensible
/// value, and having it lets a caller assemble the struct with `..Default::default()`
/// instead of spelling out every field.
impl CountMatrix {
    pub fn new(n_taxa: usize, n_samp: usize, data: Vec<f64>) -> Result<Self> {
        if data.len() != n_taxa * n_samp {
            return Err(AncombcError::Shape {
                expected: n_taxa * n_samp,
                got: data.len(),
            });
        }
        Ok(Self {
            n_taxa,
            n_samp,
            data,
        })
    }

    pub fn zeros(n_taxa: usize, n_samp: usize) -> Self {
        Self {
            n_taxa,
            n_samp,
            data: vec![0.0; n_taxa * n_samp],
        }
    }

    #[inline]
    pub fn get(&self, i: usize, j: usize) -> f64 {
        self.data[i * self.n_samp + j]
    }

    #[inline]
    pub fn set(&mut self, i: usize, j: usize, v: f64) {
        self.data[i * self.n_samp + j] = v;
    }

    pub fn to_rmatrix(&self) -> RMatrix {
        RMatrix::from_row_major(self.n_taxa, self.n_samp, self.data.clone())
    }

    /// `.data_core`: discard taxa below `prv_cut` and samples below `lib_cut`.
    ///
    /// The reference applies the two filters in a specific order that matters for
    /// library sizes (a taxon's contribution to a sample's library size is
    /// measured *after* the prevalence filter) and reuses the sample set computed
    /// for the unaggregated table when filtering the aggregated one.
    pub fn filter(&self, prv_cut: f64, lib_cut: f64) -> Result<Filtered> {
        // The first `.data_core` call passes `tax_keep = NULL`, so the only
        // reachable error here is "no taxa remain". The
        // AllTaxaStructuralZeros branch belongs to the aggregated pass, where
        // `tax_keep` arrives from the structural-zero screen.
        let tax_keep = self.taxa_above_prevalence(prv_cut)?;
        if tax_keep.is_empty() {
            return Err(AncombcError::NoTaxaRemain);
        }
        let samp_keep = self.samples_above_library_size(&tax_keep, lib_cut);
        if samp_keep.is_empty() {
            return Err(AncombcError::NoSamplesRemain);
        }
        Ok(Filtered {
            counts: self.select(&tax_keep, &samp_keep),
            taxa: tax_keep,
            samples: samp_keep,
        })
    }

    /// The second `.data_core` pass, as the reference calls it for the aggregated
    /// table: the structural-zero taxa are removed *first*, then the prevalence
    /// filter is recomputed on what remains, and the sample set is either reused
    /// from the first pass or recomputed.
    ///
    /// `pre_taxa` are indices into `self`; the returned `taxa` and `samples` are
    /// also indices into `self`, so the caller can align them with the unaggregated
    /// table.
    pub fn filter_aggregated(
        &self,
        pre_taxa: Option<&[usize]>,
        pre_samples: Option<&[usize]>,
        prv_cut: f64,
        lib_cut: f64,
    ) -> Result<Filtered> {
        // Step 1: drop structural-zero taxa (or take all of them).
        let after_struct: Vec<usize> = match pre_taxa {
            Some(idx) => idx.to_vec(),
            None => (0..self.n_taxa).collect(),
        };

        // Step 2: recompute prevalence on that subset and apply the cutoff.
        let prevalence_on = |taxa: &[usize]| -> Vec<f64> {
            taxa.iter()
                .map(|&t| {
                    let r = self.row(t);
                    let nonzero = r.iter().filter(|&&v| !v.is_nan() && v != 0.0).count();
                    let observed = r.iter().filter(|&&v| !v.is_nan()).count();
                    if observed == 0 {
                        f64::NAN
                    } else {
                        nonzero as f64 / observed as f64
                    }
                })
                .collect()
        };
        let prev = prevalence_on(&after_struct);
        let tax_keep: Vec<usize> = after_struct
            .iter()
            .zip(prev.iter())
            .filter(|(_, &p)| !p.is_nan() && p >= prv_cut)
            .map(|(&t, _)| t)
            .collect();

        if tax_keep.is_empty() {
            if after_struct.is_empty() {
                return Err(AncombcError::AllTaxaStructuralZeros);
            }
            return Err(AncombcError::NoTaxaRemain);
        }

        // Step 3: samples. When a sample set is supplied the reference uses it
        // verbatim and does not recompute library sizes.
        let samp_keep: Vec<usize> = match pre_samples {
            Some(idx) => idx.to_vec(),
            None => self.samples_above_library_size(&tax_keep, lib_cut),
        };
        if samp_keep.is_empty() {
            return Err(AncombcError::NoSamplesRemain);
        }

        Ok(Filtered {
            counts: self.select(&tax_keep, &samp_keep),
            taxa: tax_keep,
            samples: samp_keep,
        })
    }

    pub fn select(&self, taxa: &[usize], samples: &[usize]) -> CountMatrix {
        let mut out = CountMatrix::zeros(taxa.len(), samples.len());
        for (i, &t) in taxa.iter().enumerate() {
            for (j, &s) in samples.iter().enumerate() {
                out.set(i, j, self.get(t, s));
            }
        }
        out
    }

    /// Number of taxa having at least one non-missing count.
    pub fn taxa_with_any_observation(&self) -> usize {
        (0..self.n_taxa)
            .filter(|&i| self.row(i).iter().any(|v| !v.is_nan()))
            .count()
    }

    #[inline]
    pub fn row(&self, i: usize) -> &[f64] {
        &self.data[i * self.n_samp..(i + 1) * self.n_samp]
    }

    /// R: `rowSums(x != 0, na.rm = TRUE) / rowSums(!is.na(x))`
    pub fn prevalence(&self) -> Vec<f64> {
        (0..self.n_taxa)
            .map(|i| {
                let r = self.row(i);
                let nonzero = r.iter().filter(|&&v| !v.is_nan() && v != 0.0).count();
                let observed = r.iter().filter(|&&v| !v.is_nan()).count();
                if observed == 0 {
                    f64::NAN
                } else {
                    nonzero as f64 / observed as f64
                }
            })
            .collect()
    }

    /// `which(prevalence >= prv_cut)`, skipping NaN (R's `which` drops NA).
    fn taxa_above_prevalence(&self, prv_cut: f64) -> Result<Vec<usize>> {
        Ok(self
            .prevalence()
            .iter()
            .enumerate()
            .filter(|(_, &p)| !p.is_nan() && p >= prv_cut)
            .map(|(i, _)| i)
            .collect())
    }

    /// R: `colSums(feature_table, na.rm = TRUE)`, with `NaN` treated as missing.
    pub fn library_sizes(&self) -> Vec<f64> {
        (0..self.n_samp)
            .map(|j| {
                let mut s = 0.0;
                for i in 0..self.n_taxa {
                    let v = self.get(i, j);
                    if !v.is_nan() {
                        s += v;
                    }
                }
                s
            })
            .collect()
    }

    fn samples_above_library_size(&self, taxa: &[usize], lib_cut: f64) -> Vec<usize> {
        (0..self.n_samp)
            .filter(|&j| {
                let mut s = 0.0;
                for &i in taxa {
                    let v = self.get(i, j);
                    if !v.is_nan() {
                        s += v;
                    }
                }
                s >= lib_cut
            })
            .collect()
    }

    /// R: `O1 = data + pseudo`, then `log`, then `-Inf` becomes `NA`.
    pub fn log_center(&self, pseudo: f64) -> RMatrix {
        let mut y = RMatrix::zeros(self.n_taxa, self.n_samp);
        for i in 0..self.n_taxa {
            for j in 0..self.n_samp {
                let v = self.get(i, j) + pseudo;
                y.set(i, j, v.ln());
            }
        }
        y.map_nonfinite_in_place(|_| f64::NAN);
        let means = y.row_means_na_rm();
        y.sub_rows_in_place(&means);
        y
    }

    /// [`CountMatrix::log_center`] restricted to a subset of the rows.
    ///
    /// `log_center` then `select(taxa, all_samples)` would materialise the
    /// selected count table first and then transform it, so both are live at the
    /// same time. That is a full `n_selected x n_samp` buffer held for nothing but
    /// an index list, and on the `bm5` benchmark it is 800 MB of the run's peak
    /// resident set -- the P4 number. Reading the selected rows straight out of
    /// `self` halves it.
    ///
    /// The values and their order are identical to `log_center` on the selected
    /// table, and `.ancombc2_core` performs its own row subsetting the same way, so
    /// this is the reference's arithmetic with one fewer copy.
    pub fn log_center_rows(&self, rows: &[usize], pseudo: f64) -> RMatrix {
        let mut y = RMatrix::zeros(rows.len(), self.n_samp);
        for (i, &t) in rows.iter().enumerate() {
            let src = t * self.n_samp;
            for j in 0..self.n_samp {
                let v = self.data[src + j] + pseudo;
                y.set(i, j, v.ln());
            }
        }
        y.map_nonfinite_in_place(|_| f64::NAN);
        let means = y.row_means_na_rm();
        y.sub_rows_in_place(&means);
        y
    }

    /// The `O` table used by the non-conservative sensitivity analysis.
    ///
    /// `.ancombc2_sens_fit` operates on an `O2` that already carries the *main
    /// run's* pseudo-count, so when that pseudo is 0 the zeros are still zeros
    /// and get replaced here. It then logs and centres. Because the two
    /// pseudo-counts are different quantities, this deliberately does **not**
    /// equal [`CountMatrix::log_center`]: `log(x + c)` shifts every count, while
    /// this shifts only the zeros. They coincide exactly when the input has no
    /// zeros, which is why the equality below is asserted on such an input.
    /// [`CountMatrix::log_center_replacing_zeros`] on a sub-table, selected by row
    /// and column index, without materialising the sub-table first.
    ///
    /// Both indices are needed. The non-conservative analysis subsets the *sample*
    /// axis as well as the taxon axis, because it reads the unfiltered count table,
    /// so a rows-only variant would silently keep samples that the sample filter
    /// had already dropped -- the same taxon fitted against a different design.
    pub fn log_center_replacing_zeros_sub(
        &self,
        rows: &[usize],
        cols: &[usize],
        pseudo: f64,
    ) -> RMatrix {
        let mut y = RMatrix::zeros(rows.len(), cols.len());
        for (i, &t) in rows.iter().enumerate() {
            let src = t * self.n_samp;
            for (j, &c) in cols.iter().enumerate() {
                let v = self.data[src + c];
                let v = if v == 0.0 { pseudo } else { v };
                y.set(i, j, v.ln());
            }
        }
        y.map_nonfinite_in_place(|_| f64::NAN);
        let means = y.row_means_na_rm();
        y.sub_rows_in_place(&means);
        y
    }

    pub fn log_center_replacing_zeros(&self, pseudo: f64) -> RMatrix {
        let mut y = RMatrix::zeros(self.n_taxa, self.n_samp);
        for i in 0..self.n_taxa {
            for j in 0..self.n_samp {
                let v = self.get(i, j);
                let v = if v == 0.0 { pseudo } else { v };
                y.set(i, j, v.ln());
            }
        }
        y.map_nonfinite_in_place(|_| f64::NAN);
        let means = y.row_means_na_rm();
        y.sub_rows_in_place(&means);
        y
    }

    /// The set of usable (taxon, sample) pairs: a non-finite response and a
    /// complete design row. The design-completeness mask is folded in by the
    /// caller, which owns the design.
    pub fn observed_mask(&self, design_ok: &[bool]) -> Vec<bool> {
        let mut out = vec![false; self.n_taxa * self.n_samp];
        for i in 0..self.n_taxa {
            for j in 0..self.n_samp {
                out[i * self.n_samp + j] = self.get(i, j).is_finite() && design_ok[j];
            }
        }
        out
    }
}

/// Result of `.data_core`.
#[derive(Debug, Clone)]
pub struct Filtered {
    pub counts: CountMatrix,
    /// Taxon indices into the input table, ascending.
    pub taxa: Vec<usize>,
    /// Sample indices into the input table, ascending.
    pub samples: Vec<usize>,
}

// ---------------------------------------------------------------------------
// structural zeros
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
// sparse taxon representation
// ---------------------------------------------------------------------------

/// How the count table is stored.
///
/// `Dense` is the default and, per PLAN.md section 9, the right answer for the
/// hot loops: the MLE, the sandwich and the group regressions are dense kernels
/// over row subsets, and a sparse layer only adds an indirection to them. This
/// enum exists because that plan section makes the dense/sparse choice
/// *measurable* rather than assumed: `SparseTaxa` is here so both can be run on
/// the same input and the comparison can be recorded instead of argued.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Representation {
    #[default]
    Dense,
    /// Taxon-major compressed sparse rows: only non-zero, non-missing counts are
    /// stored, as `(sample, value)` pairs.
    SparseTaxa,
}

/// A count table stored taxon-major, with structural zeros omitted.
///
/// Each taxon holds a slice of `cols`/`vals` of equal length, ascending by sample
/// index. A cell absent from the index is a *structural zero* or a missing value;
/// the two are distinguished by the accompanying `present` bitmap, because they
/// mean different things everywhere downstream -- a missing count is dropped from
/// a taxon mean, a zero becomes `log(0) = NA` or `log(0 + pseudo)`.
///
/// Absent-and-not-missing therefore *is* the zero, and that is the entire
/// compression: at `z` the fraction of zeros, the stored size is
/// `(1 - z) * 12` bytes per cell against 8 for dense, so the representation only
/// wins below `z = 1 - 8/12 = 33%` -- and only for the *storage*. Every consumer
/// still pays to read a dense value back, which is why this is behind a flag.
#[derive(Debug, Clone, Default)]
pub struct SparseTaxaMatrix {
    pub n_taxa: usize,
    pub n_samp: usize,
    /// Row `i` spans `starts[i]..starts[i + 1]`.
    pub starts: Vec<u32>,
    /// Sample indices, ascending within a row.
    pub cols: Vec<u32>,
    /// The stored counts.
    pub vals: Vec<f64>,
    /// One bit per *stored* entry: whether the value is a real count or `NaN`.
    ///
    /// A `NaN` count is stored rather than dropped, because dropping it would
    /// make it indistinguishable from a zero and `rowSums(!is.na(x))` would come
    /// out wrong. So this is a bitmap over `vals`, not over cells.
    pub present: Vec<bool>,
}

impl SparseTaxaMatrix {
    /// Compress a dense table.
    ///
    /// Zero counts are dropped; `NaN` counts are stored and marked in `present`.
    pub fn from_dense(c: &CountMatrix) -> Self {
        let mut out = SparseTaxaMatrix {
            n_taxa: c.n_taxa,
            n_samp: c.n_samp,
            starts: Vec::with_capacity(c.n_taxa + 1),
            cols: Vec::new(),
            vals: Vec::new(),
            present: Vec::new(),
        };
        out.starts.push(0);
        for i in 0..c.n_taxa {
            let row = c.row(i);
            for (j, &v) in row.iter().enumerate() {
                if v == 0.0 {
                    continue;
                }
                out.cols.push(j as u32);
                out.vals.push(v);
                out.present.push(!v.is_nan());
            }
            out.starts.push(out.vals.len() as u32);
        }
        out
    }

    /// Bytes of storage, against `n_taxa * n_samp * 8` for the dense form.
    pub fn stored_bytes(&self) -> usize {
        self.starts.len() * 4 + self.cols.len() * 4 + self.vals.len() * 8 + self.present.len()
    }

    pub fn dense_bytes(&self) -> usize {
        self.n_taxa * self.n_samp * 8
    }

    /// The stored entries of taxon `i`, as `(sample, value)` with the sample
    /// indices ascending.
    #[inline]
    pub fn row_entries(&self, i: usize) -> (&[u32], &[f64]) {
        let (a, b) = (self.starts[i] as usize, self.starts[i + 1] as usize);
        (&self.cols[a..b], &self.vals[a..b])
    }

    /// One cell, or `0.0` when the entry is absent -- which is what a dropped
    /// structural zero is.
    #[inline]
    pub fn get(&self, i: usize, j: usize) -> f64 {
        let (a, b) = (self.starts[i] as usize, self.starts[i + 1] as usize);
        match self.cols[a..b].binary_search(&(j as u32)) {
            Ok(k) => self.vals[a + k],
            Err(_) => 0.0,
        }
    }

    /// Whether the cell is a real count, as opposed to absent (zero) or missing.
    #[inline]
    pub fn is_present(&self, i: usize, j: usize) -> bool {
        let (a, b) = (self.starts[i] as usize, self.starts[i + 1] as usize);
        match self.cols[a..b].binary_search(&(j as u32)) {
            Ok(k) => self.present[a + k],
            Err(_) => false,
        }
    }

    /// Materialise the dense table.
    ///
    /// Every consumer of a count table in the pipeline wants dense values, so this
    /// is where a sparse table has to give its memory back. It exists so the
    /// *screening* stages can read a compressed table without paying to expand it.
    pub fn to_dense(&self) -> CountMatrix {
        let mut c = CountMatrix::zeros(self.n_taxa, self.n_samp);
        for i in 0..self.n_taxa {
            let (cols, vals) = self.row_entries(i);
            for (&j, &v) in cols.iter().zip(vals.iter()) {
                c.set(i, j as usize, v);
            }
        }
        c
    }

    /// `rowSums(x != 0, na.rm = TRUE) / rowSums(!is.na(x))`, read from the sparse
    /// form without expanding it.
    ///
    /// This is one of the two stages that genuinely benefit: a prevalence screen
    /// touches every cell but stores nothing, so on the sparse form it reads
    /// `n_nonzero + n_missing` entries instead of `n_taxa * n_samp`.
    pub fn prevalence(&self) -> Vec<f64> {
        (0..self.n_taxa)
            .map(|i| {
                // Only `NaN` is dropped-as-missing; an absent entry is a *zero*,
                // and a zero is observed. So the non-NA count is
                // `n_samp - missing`, **not** the stored count -- an all-zero taxon
                // stores nothing at all, and its prevalence is 0 over `n_samp`,
                // not `NaN`. Reading `observed` off the stored entries instead
                // made every all-zero taxon `NaN`, which the prevalence cutoff
                // then silently dropped.
                let (a, b) = (self.starts[i] as usize, self.starts[i + 1] as usize);
                let stored = b - a;
                let missing = self.present[a..b].iter().filter(|&&p| !p).count();
                let nonzero = stored - missing;
                let observed = self.n_samp - missing;
                if observed == 0 {
                    f64::NAN
                } else {
                    nonzero as f64 / observed as f64
                }
            })
            .collect()
    }

    /// `colSums(feature_table, na.rm = TRUE)`.
    ///
    /// Sparse in the *other* axis: only stored entries contribute, so this is
    /// `O(nnz)` against `O(n_taxa * n_samp)` dense. It is the clearest win in the
    /// whole representation, because a library-size sum skips every zero by
    /// construction rather than by branching on it.
    pub fn library_sizes(&self) -> Vec<f64> {
        let mut out = vec![0.0f64; self.n_samp];
        for a in 0..self.n_taxa {
            let (b, e) = (self.starts[a] as usize, self.starts[a + 1] as usize);
            for k in b..e {
                if self.present[k] {
                    out[self.cols[k] as usize] += self.vals[k];
                }
            }
        }
        out
    }

    /// These are the two tallies the structural-zero screen needs from a taxon
    /// before it ever materialises a table: `observed` is `rowSums(!is.na(x))`
    /// and `non-zero` is `rowSums(x != 0, na.rm = TRUE)`.
    ///
    /// `observed` is `n_samp - missing` rather than the stored count, because an
    /// absent entry is a *zero*, and a zero is observed. A per-group version would
    /// need the group of each absent cell, which the compressed form does not
    /// carry -- so the screen reads the per-taxon totals here and still expands
    /// for the per-group split. That is the honest limit of this representation,
    /// and it is why the flag defaults to dense.
    /// `.data_core`, reading both screens from the compressed form.
    ///
    /// Deliberately the same two steps in the same order as
    /// [`CountMatrix::filter`], because the library-size sum must be taken over
    /// the *prevalence-filtered* taxa -- doing it on the compressed form changes
    /// nothing about the arithmetic, only about how many cells are touched.
    ///
    /// The selected table is returned dense, because every consumer from
    /// `log_center` onward wants dense rows. So the compressed form is a
    /// *transport and screening* format, not a run format; the memory it saves is
    /// the reader's and the two screens', not the pipeline's.
    pub fn filter(&self, prv_cut: f64, lib_cut: f64) -> Result<Filtered> {
        let taxa: Vec<usize> = self
            .prevalence()
            .iter()
            .enumerate()
            .filter(|(_, &p)| !p.is_nan() && p >= prv_cut)
            .map(|(i, _)| i)
            .collect();
        if taxa.is_empty() {
            return Err(AncombcError::NoTaxaRemain);
        }
        // Library size over the *retained* taxa only, so the column sums are
        // accumulated taxon-major from the stored entries rather than by scanning
        // `n_samp` columns of the dense table.
        let mut lib = vec![0.0f64; self.n_samp];
        for &t in &taxa {
            let (a, b) = (self.starts[t] as usize, self.starts[t + 1] as usize);
            for k in a..b {
                if self.present[k] {
                    lib[self.cols[k] as usize] += self.vals[k];
                }
            }
        }
        let samples: Vec<usize> = (0..self.n_samp).filter(|&j| lib[j] >= lib_cut).collect();
        if samples.is_empty() {
            return Err(AncombcError::NoSamplesRemain);
        }
        // Gather the selected rows/columns straight out of the compressed form,
        // so the expansion touches each retained cell once rather than the whole
        // table twice.
        let mut data = vec![0.0f64; taxa.len() * samples.len()];
        for (i, &t) in taxa.iter().enumerate() {
            let (a, b) = (self.starts[t] as usize, self.starts[t + 1] as usize);
            let mut out_j = 0usize;
            for k in a..b {
                let s = self.cols[k] as usize;
                // Advance the output cursor to this sample, leaving skipped (i.e.
                // filtered-out) samples as the zero they are.
                while out_j < samples.len() && samples[out_j] < s {
                    out_j += 1;
                }
                if out_j < samples.len() && samples[out_j] == s {
                    data[i * samples.len() + out_j] = self.vals[k];
                }
            }
        }
        Ok(Filtered {
            counts: CountMatrix {
                n_taxa: taxa.len(),
                n_samp: samples.len(),
                data,
            },
            taxa,
            samples,
        })
    }

    /// `(observed, non-zero)` counts per taxon, read without expanding.
    pub fn observed_counts(&self) -> Vec<(usize, usize)> {
        (0..self.n_taxa)
            .map(|i| {
                let (a, b) = (self.starts[i] as usize, self.starts[i + 1] as usize);
                let stored = b - a;
                let missing = self.present[a..b].iter().filter(|&&p| !p).count();
                let observed = self.n_samp - missing;
                (observed, stored - missing)
            })
            .collect()
    }
}

/// A dense bit set, one bit per `(taxon, group)` structural-zero flag.
///
/// # Why bits and not `Vec<bool>`
///
/// `Vec<bool>` is already one bit per element, so this is not about packing the
/// flags more tightly than `Vec<bool>` would. It is about the *scan*: the flags are
/// consumed by group-wise loops -- the primary screen asks whether a taxon is
/// flagged in *any* group -- and a bit set makes "is any bit set in this row" a
/// word-at-a-time
/// operation instead of a per-element branch. On `bm5` that screen is
/// `5000 x 10` flags, which is 63 words rather than 50,000 byte loads.
///
/// The layout is taxon-major (`taxon * n_groups + group`) to match the `zero_ind`
/// table the oracle writes and that the CLI serialises, so the index a caller
/// already has is the index used here.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct ZeroBitSet {
    words: Vec<u64>,
    len: usize,
}

impl ZeroBitSet {
    /// A set of `len` bits, all clear.
    pub fn new(len: usize) -> Self {
        Self {
            words: vec![0; len.div_ceil(64)],
            len,
        }
    }

    #[inline]
    pub fn len(&self) -> usize {
        self.len
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    #[inline]
    pub fn get(&self, i: usize) -> bool {
        debug_assert!(i < self.len, "bit {i} out of range of {}", self.len);
        self.words[i >> 6] >> (i & 63) & 1 == 1
    }

    #[inline]
    pub fn set(&mut self, i: usize, v: bool) {
        debug_assert!(i < self.len, "bit {i} out of range of {}", self.len);
        let w = &mut self.words[i >> 6];
        let m = 1u64 << (i & 63);
        if v {
            *w |= m;
        } else {
            *w &= !m;
        }
    }

    /// Expand to the `Vec<bool>` layout, for the serialisers.
    pub fn to_vec(&self) -> Vec<bool> {
        (0..self.len).map(|i| self.get(i)).collect()
    }

    /// True if any bit in `start..start + n_groups` is set, i.e. whether the taxon
    /// at `taxon` is structurally zero in any group.
    ///
    /// This is the `all(zero_ind[, -1] == FALSE)` screen. It short-circuits on the
    /// first flagged group, so the common case -- a taxon that is present
    /// everywhere -- costs one test.
    #[inline]
    pub fn any_in_row(&self, taxon: usize, n_groups: usize) -> bool {
        let base = taxon * n_groups;
        (0..n_groups).any(|g| self.get(base + g))
    }
}

/// Per-taxon, per-group structural-zero flags.
#[derive(Debug, Clone)]
pub struct StructuralZeros {
    /// Group labels, one per column.
    pub groups: Vec<String>,
    /// `n_taxa x n_groups` row-major, `true` where the taxon is absent from the
    /// group.
    ///
    /// A `Vec<bool>` because this is what the CLI writes and the FFI returns: both
    /// serialise it directly, and expanding the bit set is `n_taxa * n_groups`
    /// byte stores on a table that is about to be written out anyway. The
    /// group-major *screen* uses [`Self::bits`] instead, which is the hot path.
    pub zero_ind: Vec<bool>,
    /// The same flags as a bit set, for the group-wise scans.
    pub bits: ZeroBitSet,
}

impl StructuralZeros {
    #[inline]
    pub fn get(&self, taxon: usize, group: usize) -> bool {
        self.zero_ind[taxon * self.groups.len() + group]
    }

    /// The reference drops a taxon from the primary analysis when it has a
    /// structural zero in *any* group: `tax_idx = all(zero_ind[, -1] == FALSE)`.
    pub fn taxa_without_structural_zeros(&self) -> Vec<usize> {
        let n_groups = self.groups.len();
        (0..self.zero_ind.len() / n_groups)
            .filter(|&t| !self.bits.any_in_row(t, n_groups))
            .collect()
    }
}

/// `.get_struc_zero`.
///
/// Group prevalence is the fraction of a group's samples in which the taxon is
/// present. `neg_lb` additionally flags taxa whose *asymptotic lower bound*
/// `p - 1.96 * sqrt(p (1 - p) / n)` is non-positive, which is a wider net: it
/// also catches taxa that are merely very rare in a group.
///
/// # A sample with no group
///
/// `group_index` may hold [`NO_GROUP`] for a sample whose group label is
/// missing. Such a sample is excluded from the per-group tallies: it is not
/// counted in any group's size, and it is not counted as a sample of any group
/// in which the taxon is absent. That matches R, where the row is dropped from
/// the group it has no label for, and it is why the sentinel is handled here
/// rather than by the caller -- indexing a group table with it would read out of
/// bounds, and aborting the whole analysis over one missing label is not what
/// the reference does.
/// A sample's group label is missing.
///
/// A metadata table read from disk can carry a missing group for some samples --
/// `atlas1006` has `sex == "NA"` on a handful -- and R drops such a row from the
/// analysis rather than aborting. This sentinel is how a caller says so without
/// a second channel.
pub const NO_GROUP: usize = usize::MAX;

pub fn structural_zeros(
    counts: &CountMatrix,
    group_index: &[usize],
    n_groups: usize,
    neg_lb: bool,
) -> StructuralZeros {
    let n_taxa = counts.n_taxa;
    let n_samp = counts.n_samp;

    // Tallies are `group * taxon` rather than `taxon * group`, because each group
    // is scanned to completion before the next begins: that keeps a single group's
    // accumulator pair contiguous and in cache for the whole of its pass, instead
    // of touching `n_groups` separate locations per taxon.
    //
    // `u32` rather than `f64`: these are counts of samples, bounded by `n_samp`,
    // and `f64` counts exactly up to 2^53 -- so the two give identical answers and
    // the `f64` version costs 2x the bandwidth on a table this size.
    let cell = |g: usize, t: usize| g * n_taxa + t;
    let mut present = vec![0u32; n_groups * n_taxa];
    let mut observed = vec![0u32; n_groups * n_taxa];

    // One pass over the samples, accumulating every group's tallies at once.
    //
    // The reference walks `n_groups` in the outer loop and all `n_samp` samples in
    // the inner one, so it reads the count matrix `n_groups` times. Since a sample
    // belongs to exactly one group (or to none), the visits are disjoint and the
    // order of the summation does not change the counts -- so iterating samples
    // outermost and groups innermost visits each count exactly once. On `bm5` that
    // is 1e8 reads instead of 1e9.
    for j in 0..n_samp {
        let g = group_index.get(j).copied().unwrap_or(NO_GROUP);
        // A sample with no group is not evidence of absence in any group, so it
        // contributes to no tally at all.
        //
        // `NO_GROUP` is `usize::MAX`, so it is also caught by the range check
        // below; it is tested first only to make that ordering explicit.
        if g == NO_GROUP {
            continue;
        }
        // Every other index must name a real group. The reference would fail with a
        // subscript-out-of-bounds here, so a caller that gets this wrong is
        // already broken -- but the production loop reads `group_size[g]`, and a
        // silent out-of-bounds write into the tallies would be worse than a
        // panic in a debug build and corruption in a release one.
        debug_assert!(
            g < n_groups,
            "group index {g} is not one of {n_groups} groups"
        );
        if g >= n_groups {
            continue;
        }
        let base = j;
        for i in 0..n_taxa {
            let v = counts.data[i * n_samp + base];
            if !v.is_nan() {
                let c = cell(g, i);
                observed[c] += 1;
                if v != 0.0 {
                    present[c] += 1;
                }
            }
        }
    }

    let mut group_size = vec![0.0f64; n_groups];
    for &g in group_index {
        if g != NO_GROUP {
            group_size[g] += 1.0;
        }
    }

    // `p` is the prevalence *over the group*, so its denominator is the group's
    // full size and not the number of observed counts -- an NA sample is neither a
    // presence nor an absence, but it is still a sample of the group. That is the
    // distinction `neg_lb` then draws on, using `observed` for the binomial term.
    let mut bits = ZeroBitSet::new(n_taxa * n_groups);
    let mut zero_ind = vec![false; n_taxa * n_groups];
    for g in 0..n_groups {
        if group_size[g] == 0.0 {
            // An empty group is never structurally zero: there is no evidence to
            // base a flag on. The bits are already clear.
            continue;
        }
        for i in 0..n_taxa {
            let c = cell(g, i);
            let observed_in_group = observed[c];
            let p = f64::from(present[c]) / group_size[g];
            let flag = if p == 0.0 {
                true
            } else if neg_lb && observed_in_group > 0 {
                let o = f64::from(observed_in_group);
                p - 1.96 * (p * (1.0 - p) / o).sqrt() <= 0.0
            } else {
                false
            };
            if flag {
                let idx = i * n_groups + g;
                bits.set(idx, true);
                zero_ind[idx] = true;
            }
        }
    }

    StructuralZeros {
        groups: (0..n_groups).map(|g| format!("g{g}")).collect(),
        zero_ind,
        bits,
    }
}

#[cfg(test)]
mod tests {

    #[cfg(test)]
    /// Row-major construction from rows, for tests only.
    #[cfg(test)]
    fn from_test_rows(rows: &[Vec<f64>]) -> CountMatrix {
        let n_taxa = rows.len();
        let n_samp = rows[0].len();
        let mut data = Vec::with_capacity(n_taxa * n_samp);
        for r in rows {
            assert_eq!(r.len(), n_samp, "ragged test table");
            data.extend_from_slice(r);
        }
        CountMatrix::new(n_taxa, n_samp, data).expect("test table")
    }

    fn assert_same(a: &RMatrix, b: &RMatrix, what: &str) {
        assert_eq!(a.rows, b.rows, "{what}: row count");
        assert_eq!(a.cols, b.cols, "{what}: column count");
        for i in 0..a.rows {
            for j in 0..a.cols {
                assert_eq!(
                    a.get(i, j).to_bits(),
                    b.get(i, j).to_bits(),
                    "{what}: cell ({i}, {j})"
                );
            }
        }
    }

    #[test]
    fn log_center_rows_equals_select_then_log_center() {
        let c = from_test_rows(&[
            vec![1.0, 0.0, 4.0, 9.0],
            vec![0.0, 0.0, 0.0, 0.0],
            vec![7.0, 3.0, 1.0, 2.0],
            vec![5.0, 5.0, 5.0, 5.0],
        ]);
        for rows in [vec![0usize, 1, 2, 3], vec![2, 0], vec![3, 2, 1, 0, 2]] {
            for pseudo in [0.0, 0.1, 1.0] {
                let direct = c.log_center_rows(&rows, pseudo);
                let via_select = c
                    .select(&rows, &(0..c.n_samp).collect::<Vec<_>>())
                    .log_center(pseudo);
                assert_same(&direct, &via_select, "log_center_rows");
            }
        }
    }

    #[test]
    fn zero_replacing_sub_equals_select_then_transform() {
        let c = from_test_rows(&[
            vec![1.0, 0.0, 4.0, 9.0, 2.0],
            vec![0.0, 0.0, 0.0, 0.0, 0.0],
            vec![7.0, 3.0, 1.0, 2.0, 6.0],
            vec![5.0, 5.0, 5.0, 5.0, 5.0],
        ]);
        // Both axes selected, and deliberately not the identity on either: the
        // point of the test is that a rows-only variant would keep samples the
        // sample filter had dropped.
        for (rows, cols) in [
            (vec![2usize, 0], vec![1usize, 3, 4]),
            (vec![3, 2, 1, 0], vec![0, 1, 2, 3, 4]),
            (vec![1usize], vec![2usize, 0]),
        ] {
            for pseudo in [0.0, 0.01, 0.5] {
                let direct = c.log_center_replacing_zeros_sub(&rows, &cols, pseudo);
                let via_select = c.select(&rows, &cols).log_center_replacing_zeros(pseudo);
                assert_same(&direct, &via_select, "log_center_replacing_zeros_sub");
            }
        }
    }

    /// A sample with no group label must be excluded from the per-group tallies.
    ///
    /// `atlas1006` really does have samples whose `sex` is missing, and this is
    /// the code that reads it. Two things have to hold, and both are observable:
    /// an unlabelled sample must not make a taxon look structurally zero in
    /// *every* group, and the group sizes must not gain a phantom level.
    ///
    /// `group_index` covers every sample, as the caller guarantees, so the two
    /// arms of the comparison differ only in what the last two samples are
    /// labelled.
    #[test]
    fn a_sample_with_no_group_is_not_evidence_of_absence() {
        // Taxon 1 is absent from the last sample of each pair. If the two
        // unlabelled samples counted as absences, taxon 1 would be flagged zero
        // in both groups.
        let data = vec![
            10.0, 10.0, 10.0, 10.0, 10.0, 10.0, // taxon 0: present everywhere
            5.0, 5.0, 5.0, 5.0, 0.0, 0.0, // taxon 1: absent from the last of each
        ];
        let counts = CountMatrix::new(2, 6, data).unwrap();
        let labelled = [0usize, 0, 1, 1, 0, 1];
        let unlabelled = [0usize, 0, 1, 1, NO_GROUP, NO_GROUP];

        let a = structural_zeros(&counts, &labelled, 2, false);
        let b = structural_zeros(&counts, &unlabelled, 2, false);
        assert_eq!(
            a.zero_ind, b.zero_ind,
            "an unlabelled sample must not change any structural-zero flag"
        );
        // Taxon 1 is present in the first sample of each group, so it is not
        // structurally zero in either.
        // `zero_ind` is `n_taxa * n_groups` row-major, so taxon `i` in group `g`
        // is at `i * n_groups + g`.
        assert!(!a.zero_ind[2], "taxon 1 in group 0");
        assert!(!a.zero_ind[3], "taxon 1 in group 1");
    }

    /// A taxon absent from every sample *including* the unlabelled one is still
    /// flagged, so the sentinel is not simply disabling the check.
    #[test]
    fn the_sentinel_does_not_rescue_a_genuinely_absent_taxon() {
        let counts = CountMatrix::new(
            2,
            5,
            vec![10.0, 10.0, 10.0, 10.0, 10.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        )
        .unwrap();
        let sz = structural_zeros(&counts, &[0, 0, 1, 1, NO_GROUP], 2, false);
        assert!(sz.zero_ind[2], "taxon 1 in group 0");
        assert!(sz.zero_ind[3], "taxon 1 in group 1");
    }
    use super::*;

    /// Row-major count matrix from a row-wise literal.
    fn m(rows: Vec<Vec<f64>>) -> CountMatrix {
        let nr = rows.len();
        let nc = rows[0].len();
        let mut data = vec![0.0; nr * nc];
        for (i, r) in rows.iter().enumerate() {
            for (j, v) in r.iter().enumerate() {
                data[i * nc + j] = *v;
            }
        }
        CountMatrix::new(nr, nc, data).unwrap()
    }

    #[test]
    fn prevalence_uses_observed_samples_as_the_denominator() {
        // taxon 0: 1 non-zero of 3 observed -> 1/3
        // taxon 1: 3 of 3                   -> 1
        // taxon 2: 0 of 3                   -> 0
        // taxon 3: 1 non-zero of 1 observed  -> 1  (two NAs ignored)
        let c = m(vec![
            vec![1.0, 0.0, 0.0],
            vec![1.0, 2.0, 3.0],
            vec![0.0, 0.0, 0.0],
            vec![5.0, f64::NAN, f64::NAN],
        ]);
        let p = c.prevalence();
        assert!((p[0] - 1.0 / 3.0).abs() < 1e-15);
        assert!((p[1] - 1.0).abs() < 1e-15);
        assert!((p[2] - 0.0).abs() < 1e-15);
        assert!((p[3] - 1.0).abs() < 1e-15, "NA must not count as absent");
    }

    #[test]
    fn prevalence_cutoff_is_inclusive() {
        // prevalence exactly 0.5 with prv_cut = 0.5 must be kept: the reference
        // filters on `prevalence >= prv_cut`, not `>`.
        let c = m(vec![vec![1.0, 0.0, 1.0, 0.0]]);
        assert!((c.prevalence()[0] - 0.5).abs() < 1e-15);
        let f = c.filter(0.5, 0.0).unwrap();
        assert_eq!(f.taxa, vec![0], "prevalence >= prv_cut keeps the taxon");
        // one notch above and it is dropped, which surfaces as the error
        let e = c.filter(0.5000001, 0.0).unwrap_err();
        assert!(matches!(e, AncombcError::NoTaxaRemain), "got {e:?}");
    }

    #[test]
    fn library_size_filter_uses_the_taxon_filtered_table() {
        let c = m(vec![vec![10.0, 0.0, 0.0, 0.0], vec![0.0, 5.0, 5.0, 5.0]]);
        // with prv_cut = 0.25, taxon 0 has prevalence 0.25 and survives;
        // library sizes are then (10, 5, 5, 5) and lib_cut = 6 keeps sample 0
        let f = c.filter(0.25, 6.0).unwrap();
        assert_eq!(f.samples, vec![0]);
    }

    #[test]
    fn an_all_zero_table_reports_no_taxa_remain_on_the_first_pass() {
        // Every prevalence is 0, so the first `.data_core` (tax_keep = NULL)
        // fails with "No taxa remain under the current cutoff".
        let c = m(vec![vec![0.0, 0.0], vec![0.0, 0.0]]);
        let e = c.filter(0.1, 0.0).unwrap_err();
        assert!(matches!(e, AncombcError::NoTaxaRemain), "got {e:?}");
    }

    #[test]
    fn an_empty_structural_zero_set_reports_all_taxa_structural_zeros() {
        // The second `.data_core` receives `tax_keep` from the structural-zero
        // screen; when that screen rejects everything, the reference raises a
        // different error.
        let c = m(vec![vec![0.0, 0.0], vec![0.0, 0.0]]);
        let e = c.filter_aggregated(Some(&[]), None, 0.10, 0.0).unwrap_err();
        assert!(
            matches!(e, AncombcError::AllTaxaStructuralZeros),
            "got {e:?}"
        );
    }

    #[test]
    fn log_center_matches_r() {
        // R: log(c(1, e, 10)) - mean  ->  the centred log vector
        let c = m(vec![vec![1.0, std::f64::consts::E, 10.0]]);
        let y = c.log_center(0.0);
        let expect = [1.0f64.ln(), 1.0, 10.0f64.ln()];
        let mean = expect.iter().sum::<f64>() / 3.0;
        for j in 0..3 {
            assert!((y.get(0, j) - (expect[j] - mean)).abs() < 1e-15, "j={j}");
        }
    }

    #[test]
    fn log_center_makes_zeros_nan_and_uses_finite_only_means() {
        // R: o1 <- log(c(0, 1, 4)); o1[is.infinite(o1)] <- NA
        //      y1 <- o1 - rowMeans(o1, na.rm = TRUE)
        // rowMeans(na.rm = TRUE) is (log 1 + log 4) / 2 = log 2
        let c = m(vec![vec![0.0, 1.0, 4.0]]);
        let y = c.log_center(0.0);
        assert!(y.get(0, 0).is_nan(), "log 0 must become NaN");
        let mean = (1.0f64.ln() + 4.0f64.ln()) / 2.0;
        assert!((mean - 2.0f64.ln()).abs() < 1e-15);
        assert!(
            (y.get(0, 1) - (0.0 - mean)).abs() < 1e-15,
            "{}",
            y.get(0, 1)
        );
        assert!((y.get(0, 2) - (4.0f64.ln() - mean)).abs() < 1e-15);
        // the centred finite values sum to zero
        assert!((y.get(0, 1) + y.get(0, 2)).abs() < 1e-15);
    }

    #[test]
    fn log_center_with_pseudo_makes_every_entry_finite() {
        let c = m(vec![vec![0.0, 1.0, 4.0]]);
        let y = c.log_center(0.5);
        for j in 0..3 {
            assert!(y.get(0, j).is_finite(), "j={j}");
        }
    }

    #[test]
    fn replacing_zeros_and_adding_a_pseudo_are_different_transforms() {
        // The two paths are deliberately distinct. `log_center(p)` computes
        // log(x + p) for every count; `log_center_replacing_zeros(p)` computes
        // log(p) where x == 0 and log(x) elsewhere. They agree only at p = 0 on
        // a zero-free input, which is the degenerate case.
        let c = m(vec![vec![0.0, 1.0, 4.0]]);
        let zero_free = m(vec![vec![1.0, 4.0, 9.0]]);
        for pc in [0.5, 1.0, 2.0] {
            // with zeros present, both the zero cell and the non-zero cells differ
            let a = c.log_center(pc);
            let b = c.log_center_replacing_zeros(pc);
            assert!(
                (a.get(0, 1) - b.get(0, 1)).abs() > 1e-6,
                "pc={pc}: non-zero cell"
            );
            // even with no zeros, adding a pseudo-count shifts every value
            let a2 = zero_free.log_center(pc);
            let b2 = zero_free.log_center_replacing_zeros(pc);
            assert!(
                (a2.get(0, 0) - b2.get(0, 0)).abs() > 1e-6,
                "pc={pc}: zero-free"
            );
        }
        // and they agree in the degenerate case
        let a = zero_free.log_center(0.0);
        let b = zero_free.log_center_replacing_zeros(0.0);
        for j in 0..3 {
            assert!((a.get(0, j) - b.get(0, j)).abs() < 1e-15, "j={j}");
        }
    }

    #[test]
    fn replacing_zeros_with_a_pseudo_makes_zeros_finite_when_pseudo_is_positive() {
        let c = m(vec![vec![0.0, 1.0, 4.0]]);
        // pseudo = 0 keeps the zero missing, as in the main run
        let d0 = c.log_center_replacing_zeros(0.0);
        assert!(d0.get(0, 0).is_nan());
        // a positive pseudo-count brings it back
        let d1 = c.log_center_replacing_zeros(0.01);
        assert!(d1.get(0, 0).is_finite());
    }

    /// The reference definition, transcribed literally: group-major outer loop,
    /// all samples inner, `f64` tallies.
    ///
    /// The production `structural_zeros` is a single sample-major pass over `u32`
    /// tallies, which is only the same function if visiting each sample once,
    /// instead of once per group, leaves the counts unchanged. This pins that
    /// equivalence against the shape the oracle actually computes, on inputs that
    /// cover the cases where it could fail: empty groups, NA counts, NA group
    /// labels, and `neg_lb` on both sides.
    fn reference_structural_zeros(
        counts: &CountMatrix,
        group_index: &[usize],
        n_groups: usize,
        neg_lb: bool,
    ) -> Vec<bool> {
        let mut present = vec![0.0f64; counts.n_taxa * counts.n_samp];
        for i in 0..counts.n_taxa {
            for j in 0..counts.n_samp {
                let v = counts.get(i, j);
                present[i * counts.n_samp + j] = if v.is_nan() || v == 0.0 { 0.0 } else { 1.0 };
            }
        }
        let mut group_size = vec![0.0f64; n_groups];
        for &g in group_index {
            if g != NO_GROUP {
                group_size[g] += 1.0;
            }
        }
        let mut zero_ind = vec![false; counts.n_taxa * n_groups];
        for i in 0..counts.n_taxa {
            for g in 0..n_groups {
                if group_size[g] == 0.0 {
                    continue;
                }
                let mut present_in_group = 0.0;
                let mut observed_in_group = 0.0;
                for j in 0..counts.n_samp {
                    if group_index[j] != g || group_index[j] == NO_GROUP {
                        continue;
                    }
                    let v = counts.get(i, j);
                    if !v.is_nan() {
                        observed_in_group += 1.0;
                        present_in_group += present[i * counts.n_samp + j];
                    }
                }
                let p = present_in_group / group_size[g];
                if p == 0.0 {
                    zero_ind[i * n_groups + g] = true;
                } else if neg_lb && observed_in_group > 0.0 {
                    let lo = p - 1.96 * (p * (1.0 - p) / observed_in_group).sqrt();
                    if lo <= 0.0 {
                        zero_ind[i * n_groups + g] = true;
                    }
                }
            }
        }
        zero_ind
    }

    /// A deterministic spread of shapes, zeros, NAs and group labels.
    fn lcg(seed: &mut u64) -> u64 {
        *seed = seed
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        *seed >> 33
    }

    fn mixed_counts(n_taxa: usize, n_samp: usize, seed: u64) -> CountMatrix {
        let mut st = seed;
        let mut data = Vec::with_capacity(n_taxa * n_samp);
        for _ in 0..n_taxa * n_samp {
            let r = lcg(&mut st);
            data.push(match r % 7 {
                0 => 0.0,
                1 => f64::NAN,
                _ => (r % 100) as f64,
            });
        }
        CountMatrix::new(n_taxa, n_samp, data).unwrap()
    }

    #[test]
    fn the_single_pass_matches_the_reference_triple_loop() {
        let mut st = 0x5eed_1234u64;
        // Covers an empty group (a group index no sample maps to), NA group labels,
        // NA counts, and taxa that are absent from some groups and present in
        // others. Several group counts, including more groups than are populated,
        // because an unpopulated group is the branch that must stay all-false.
        for &(n_taxa, n_samp, n_groups) in &[
            (1usize, 1usize, 1usize),
            (4, 9, 2),
            (7, 12, 3),
            (5, 8, 5),
            (9, 6, 2),
            (11, 20, 4),
        ] {
            for seed in 0..6u64 {
                let counts = mixed_counts(n_taxa, n_samp, seed * 7919 + n_taxa as u64);
                let mut group_index = Vec::with_capacity(n_samp);
                for _ in 0..n_samp {
                    let r = lcg(&mut st);
                    // A third of the samples are unlabelled, and some land in groups
                    // that other samples never use.
                    // Send some samples to a group that others never use, so the
                    // "group exists but nothing maps to it" branch is covered, and
                    // keep every index in range for the current `n_groups`.
                    group_index.push(match r % 5 {
                        0 => NO_GROUP,
                        1 => (n_groups / 2 + 1) % n_groups,
                        _ => r as usize % n_groups,
                    });
                }
                for &neg_lb in &[false, true] {
                    let got = structural_zeros(&counts, &group_index, n_groups, neg_lb);
                    let want = reference_structural_zeros(&counts, &group_index, n_groups, neg_lb);
                    assert_eq!(
                        got.zero_ind, want,
                        "shape {n_taxa}x{n_samp}x{n_groups} seed {seed} neg_lb {neg_lb}"
                    );
                    // The bit set must agree with the `Vec<bool>` it was built from,
                    // since the screen uses one and the serialisers use the other.
                    assert_eq!(got.bits.to_vec(), want, "bitset disagrees at the same case");
                }
            }
        }
    }

    #[test]
    fn the_any_group_screen_agrees_with_the_reference_expression() {
        let counts = mixed_counts(6, 10, 99);
        let group_index = vec![0, 0, 1, 1, 1, NO_GROUP, 2, 2, 2, 2];
        let sz = structural_zeros(&counts, &group_index, 3, false);
        let want: Vec<usize> = (0..6).filter(|&t| (0..3).all(|g| !sz.get(t, g))).collect();
        assert_eq!(sz.taxa_without_structural_zeros(), want);
    }

    #[test]
    fn a_bitset_round_trips_through_its_vec_form() {
        let n = 200;
        let mut b = ZeroBitSet::new(n);
        for i in (0..n).step_by(7) {
            b.set(i, true);
        }
        let v = b.to_vec();
        assert_eq!(v.len(), n);
        for i in 0..n {
            assert_eq!(v[i], i % 7 == 0, "bit {i}");
            assert_eq!(b.get(i), v[i], "bit {i}");
        }
        // Clearing must clear the bit, not only fail to set it.
        b.set(0, false);
        assert!(!b.get(0));
        assert!(!b.to_vec()[0]);
        assert_eq!(b.len(), n);
        assert!(!b.is_empty());
        assert!(ZeroBitSet::new(0).is_empty());
    }

    #[test]
    fn a_row_screen_sees_flags_past_a_word_boundary() {
        // 70 groups puts the second row's flags in a second word, so a screen that
        // only inspected the first word of the row would miss them.
        let n_groups = 70;
        let mut b = ZeroBitSet::new(3 * n_groups);
        assert!(!b.any_in_row(0, n_groups));
        b.set(n_groups + 69, true);
        assert!(
            b.any_in_row(1, n_groups),
            "a flag in the last group must count"
        );
        assert!(!b.any_in_row(0, n_groups));
        assert!(!b.any_in_row(2, n_groups));
    }

    /// A table with a mix of zeros, `NaN`s and real counts, and a known answer
    /// for every screen that reads it.
    fn screen_fixture() -> CountMatrix {
        // taxon 0: all zero        -> prevalence 0,          no observation
        // taxon 1: all present     -> prevalence 1
        // taxon 2: half zero       -> prevalence 0.5, and an NA in the tail
        // taxon 3: two NAs, rest 1 -> prevalence 1 over 4 observed
        let n = 4usize;
        let mut data = vec![0.0; n * 6];
        for j in 0..6 {
            data[6 + j] = 2.0;
        }
        for j in 0..3 {
            data[12 + j] = 3.0;
        }
        data[12 + 3] = f64::NAN;
        data[12 + 4] = 0.0;
        data[12 + 5] = 3.0;
        for j in 0..6 {
            data[18 + j] = 1.0;
        }
        data[18] = f64::NAN;
        data[19] = f64::NAN;
        CountMatrix::new(n, 6, data).unwrap()
    }

    #[test]
    fn a_sparse_table_round_trips_to_the_same_dense_one() {
        let c = screen_fixture();
        let sp = SparseTaxaMatrix::from_dense(&c);
        let back = sp.to_dense();
        assert_eq!(back.n_taxa, c.n_taxa);
        assert_eq!(back.n_samp, c.n_samp);
        for i in 0..c.n_taxa {
            for j in 0..c.n_samp {
                // Bit-for-bit, `NaN` included: this is the property that lets the
                // flag be switched without changing a single reported number.
                assert!(
                    back.get(i, j).to_bits() == c.get(i, j).to_bits(),
                    "cell ({i},{j}): sparse {} vs dense {}",
                    back.get(i, j),
                    c.get(i, j)
                );
            }
        }
    }

    #[test]
    fn the_sparse_screens_agree_with_the_dense_ones() {
        let c = screen_fixture();
        let sp = SparseTaxaMatrix::from_dense(&c);
        for (i, (&d, &s)) in c
            .prevalence()
            .iter()
            .zip(sp.prevalence().iter())
            .enumerate()
        {
            assert!(
                (d.is_nan() && s.is_nan()) || d.to_bits() == s.to_bits(),
                "prevalence taxon {i}: dense {d} sparse {s}"
            );
        }
        for (j, (&d, &s)) in c
            .library_sizes()
            .iter()
            .zip(sp.library_sizes().iter())
            .enumerate()
        {
            assert_eq!(d, s, "library size sample {j}");
        }
        // Observed (non-NA) cells per taxon, which is what the structural-zero
        // screen needs from a table before it ever materialises one.
        for i in 0..c.n_taxa {
            let row = c.row(i);
            // `rowSums(!is.na(x))`: a zero counts as observed, which is exactly the
            // case an all-zero taxon exercises.
            let dense_observed = row.iter().filter(|v| !v.is_nan()).count();
            // `rowSums(x != 0, na.rm = TRUE)`.
            let dense_nonzero = row.iter().filter(|v| !v.is_nan() && **v != 0.0).count();
            let (observed, nonzero) = sp.observed_counts()[i];
            assert_eq!(
                observed, dense_observed,
                "observed taxon {i}: dense {dense_observed} sparse {observed}"
            );
            assert_eq!(
                nonzero, dense_nonzero,
                "non-zero taxon {i}: dense {dense_nonzero} sparse {nonzero}"
            );
        }
    }

    /// The break-even point is the whole reason this is behind a flag, so it is
    /// asserted rather than left as a claim in a doc comment.
    ///
    /// A stored entry costs a `u32` index plus an `f64` value plus a presence bit,
    /// i.e. more than the `8` bytes of the dense cell it replaces. At zero rate
    /// `z` the sparse form is smaller only when `(1 - z) * 12 < 8`.
    #[test]
    fn the_sparse_form_only_shrinks_above_the_break_even_zero_rate() {
        let n_taxa = 4usize;
        let n_samp = 8usize;
        for &(n_nonzero, expect_smaller) in &[(0usize, true), (8, true), (16, true), (32, false)] {
            let mut data = vec![0.0; n_taxa * n_samp];
            // `n_nonzero` of the 32 cells are non-zero, placed in the first taxa.
            for k in 0..n_nonzero {
                data[k] = 1.0;
            }
            let c = CountMatrix::new(n_taxa, n_samp, data).unwrap();
            let sp = SparseTaxaMatrix::from_dense(&c);
            let zero_rate = 1.0 - (n_nonzero as f64 / (n_taxa * n_samp) as f64);
            let break_even = 1.0 - 8.0 / 12.0;
            assert_eq!(
                sp.stored_bytes() < sp.dense_bytes(),
                expect_smaller,
                "zero rate {zero_rate:.2} (break-even {break_even:.2}): {} vs {}",
                sp.stored_bytes(),
                sp.dense_bytes()
            );
        }
    }

    #[test]
    fn an_empty_table_compresses_to_nothing_but_the_offsets() {
        let c = CountMatrix::zeros(3, 5);
        let sp = SparseTaxaMatrix::from_dense(&c);
        assert_eq!(sp.starts, vec![0, 0, 0, 0]);
        assert!(sp.cols.is_empty() && sp.vals.is_empty());
        assert_eq!(sp.stored_bytes(), (3 + 1) * 4);
        for i in 0..3 {
            for j in 0..5 {
                assert_eq!(sp.get(i, j), 0.0);
                assert!(!sp.is_present(i, j));
            }
        }
    }

    /// The flag has to be invisible from outside: both representations must
    /// return the same retained taxa, the same retained samples and the same
    /// selected table, bit for bit.
    ///
    /// Asserted on inputs that exercise the paths where they could differ -- an
    /// all-zero taxon (absent from the sparse form entirely), `NA` counts (stored,
    /// not dropped), and a library-size cutoff that has to be taken *after* the
    /// prevalence filter -- and over a grid of cutoffs, because a bug here would
    /// show up at one threshold and not another.
    #[test]
    fn both_representations_filter_identically() {
        let mut st = 0xfeed_1234u64;
        let lcg = |st: &mut u64| {
            *st = st
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            *st >> 33
        };
        for &(n_taxa, n_samp) in &[(1usize, 1usize), (3, 5), (12, 9), (40, 20)] {
            for seed in 0..5u64 {
                let mut data = Vec::with_capacity(n_taxa * n_samp);
                for _ in 0..n_taxa * n_samp {
                    let r = lcg(&mut st);
                    data.push(match (seed * 31 + r) % 11 {
                        0..=3 => 0.0,
                        4 => f64::NAN,
                        _ => (r % 50) as f64,
                    });
                }
                let c = CountMatrix::new(n_taxa, n_samp, data).unwrap();
                let sp = SparseTaxaMatrix::from_dense(&c);
                for &prv in &[0.0, 0.1, 0.5, 0.9] {
                    for &lib in &[0.0, 1.0, 100.0] {
                        let d = c.filter(prv, lib);
                        let s = sp.filter(prv, lib);
                        match (d, s) {
                            (Ok(a), Ok(b)) => {
                                assert_eq!(a.taxa, b.taxa, "taxa {n_taxa}x{n_samp} seed {seed}");
                                assert_eq!(a.samples, b.samples, "samples");
                                // Bitwise, not `==`: the tables contain `NaN`
                                // counts, and `NaN != NaN` would report a
                                // difference on every shape that has one.
                                let same = a.counts.data.len() == b.counts.data.len()
                                    && a.counts
                                        .data
                                        .iter()
                                        .zip(b.counts.data.iter())
                                        .all(|(x, y)| x.to_bits() == y.to_bits());
                                assert!(same, "selected cells differ");
                            }
                            (Err(_), Err(_)) => {}
                            (a, b) => panic!(
                                "{n_taxa}x{n_samp} seed {seed} prv {prv} lib {lib}: \
                                 one representation succeeded ({:?}) and the other did not ({:?})",
                                a.is_ok(),
                                b.is_ok()
                            ),
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn structural_zero_detects_group_absence() {
        let c = m(vec![
            vec![0.0, 0.0, 5.0, 5.0], // absent from group 0
            vec![1.0, 1.0, 5.0, 5.0], // present everywhere
        ]);
        let g = vec![0usize, 0, 1, 1];
        let sz = structural_zeros(&c, &g, 2, false);
        assert!(sz.get(0, 0), "taxon 0 absent from group 0");
        assert!(!sz.get(0, 1));
        assert!(!sz.get(1, 0) && !sz.get(1, 1));
        assert_eq!(sz.taxa_without_structural_zeros(), vec![1]);
    }

    #[test]
    fn neg_lb_catches_rare_taxa() {
        // taxon 0: 1 of 10 samples present in group 0
        //   p = 0.1, p_lo = 0.1 - 1.96*sqrt(0.1*0.9/10) = 0.1 - 0.186 < 0
        let mut rows = vec![vec![0.0; 10]];
        rows[0][0] = 1.0;
        let c = m(rows);
        let g = vec![0usize; 10];
        assert!(!structural_zeros(&c, &g, 1, false).get(0, 0));
        assert!(
            structural_zeros(&c, &g, 1, true).get(0, 0),
            "neg_lb must widen the net to rare taxa"
        );
    }

    #[test]
    fn observed_mask_folds_in_design_completeness() {
        let c = m(vec![vec![1.0, 0.0, 2.0]]);
        let mask = c.observed_mask(&[true, false, true]);
        assert_eq!(mask, vec![true, false, true]);
    }

    #[test]
    fn filter_produces_the_reference_index_sets() {
        let c = m(vec![
            vec![0.0, 1.0, 1.0, 1.0], // prevalence 0.75
            vec![0.0, 0.0, 0.0, 0.0], // prevalence 0.0
            vec![1.0, 1.0, 1.0, 1.0], // prevalence 1.0
        ]);
        let f = c.filter(0.5, 0.0).unwrap();
        assert_eq!(f.taxa, vec![0, 2]);
        assert_eq!(f.samples, vec![0, 1, 2, 3]);
        assert_eq!(f.counts.n_taxa, 2);
        assert_eq!(f.counts.get(0, 0), 0.0);
    }
}
