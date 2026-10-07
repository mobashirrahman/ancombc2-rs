//! The typed transport between R's memory and the numerical core.
//!
//! # What this module is responsible for
//!
//! Turning a set of borrowed C pointers -- which came from R -- into *owned* Rust
//! values, having first proved that every buffer's length, shape and layout agree
//! with what the caller claimed. Nothing else. No numerics happen here; a stage
//! that wanted to be a `Vec` still has to build it.
//!
//! # Why the copy is not a pessimisation to be optimised away
//!
//! Two reasons, and the second is the one that matters:
//!
//! 1. R's `REALSXP` is a `double` array with a `NA` payload that is a *specific*
//!    NaN bit pattern. Copying the 64 bits preserves `NA_real_`, `NaN` (a
//!    different NaN), `+Inf`, `-Inf` and `-0.0` exactly. Any conversion that goes
//!    through `is_nan()` and back does not.
//! 2. **No pointer that came from R may outlive this module.** Rayon jobs run on
//!    threads the R runtime knows nothing about, and R may move or collect an
//!    object at any allocation. A borrowed `&[f64]` handed to `par_iter` is a
//!    data race waiting for a GC. So the transport ends at an owned buffer and
//!    [`crate::owned::Owned`] is the only thing a parallel stage may be given.
//!
//! # Layout
//!
//! R matrices are **column-major** and so is every `Matrix` in `ancombc2-core`.
//! There is no layout conversion anywhere in this crate, and that is deliberate:
//! a conversion is a place to transpose something by accident, and the count
//! matrix and the design matrix would then disagree with no test to say so.
//! [`Layout`] exists so that a future row-major producer has to say which it is.

use std::fmt;

/// Memory layout of a 2-D buffer, in the order its elements appear.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Layout {
    /// Column-major: element `(r, c)` is at `r + c * rows`. R's own order, and
    /// the order every buffer on this boundary is in.
    ColumnMajor,
    /// Row-major: element `(r, c)` is at `c + r * cols`.
    RowMajor,
}

impl Layout {
    /// The name, for a diagnostic. Reported rather than compared: a layout
    /// mismatch is caught by the element lookup, not by the label.
    pub fn name(self) -> &'static str {
        match self {
            Layout::ColumnMajor => "column-major",
            Layout::RowMajor => "row-major",
        }
    }
}

/// Everything that can be wrong with a transport call.
///
/// Every variant carries enough to name the offending buffer and the two numbers
/// that disagreed. "invalid input" is not an error message; it is a refusal to
/// help.
#[derive(Debug, Clone, PartialEq)]
pub enum BridgeError {
    /// A pointer the caller promised was non-null was null.
    /// A pointer the caller promised was non-null was null.
    NullBuffer {
        /// Which buffer.
        field: &'static str,
    },
    /// A length disagrees with the shape it is supposed to describe.
    /// A length disagrees with the shape it is supposed to describe.
    ShapeMismatch {
        /// Which buffer.
        field: &'static str,
        /// The declared `nrow`.
        rows: usize,
        /// The declared `ncol`.
        cols: usize,
        /// The number of elements actually present.
        len: usize,
    },
    /// Two shapes that must agree do not.
    /// Two shapes that must agree do not.
    IncompatibleShapes {
        /// The first buffer.
        a: &'static str,
        /// Its `(nrow, ncol)`.
        a_shape: (usize, usize),
        /// The second buffer.
        b: &'static str,
        /// Its `(nrow, ncol)`.
        b_shape: (usize, usize),
    },
    /// A mask is not the length it claims.
    /// A buffer's length is not the one its role requires.
    LengthMismatch {
        /// Which buffer.
        field: &'static str,
        /// What the role requires.
        expected: usize,
        /// What it has.
        got: usize,
    },
    /// A mask holds a value outside its domain.
    /// A mask or a count control holds a value outside its domain.
    MaskValue {
        /// Which field.
        field: &'static str,
        /// Where.
        index: usize,
        /// What it holds.
        value: i64,
    },
    /// A group index outside `1 ..= n_levels`, other than 0 which means "absent".
    /// A group index outside `1 ..= n_levels`; `0` means absent and is legal.
    GroupIndex {
        /// Which sample.
        index: usize,
        /// What it holds.
        value: i32,
        /// How many levels there are.
        n_levels: usize,
    },
    /// The packed name buffer's offsets are not a monotone partition.
    /// The packed name buffer's offsets are not a monotone partition of its bytes.
    NameBuffer(
        /// Which buffer.
        &'static str,
    ),
    /// A scalar was not finite where a finite value is required.
    /// A scalar was not finite where a finite value is required.
    NonFiniteScalar {
        /// Which control.
        field: &'static str,
        /// What it holds.
        value: f64,
    },
    /// An argument was rejected by the numerical core.
    ///
    /// The core's errors are its own type and this crate's are this type, and the
    /// boundary is the only place the two meet. The message is carried verbatim so
    /// the diagnostic an R user sees is the one the core wrote.
    Args(String),
    /// A count is larger than any result this package can produce.
    ///
    /// Distinct from [`BridgeError::LengthMismatch`] because nothing disagrees
    /// here: the count is self-consistent, merely impossible. It is also what a
    /// caller whose struct layout disagrees with this crate's looks like, so it
    /// is reported by name rather than being left to surface as an allocation
    /// failure.
    TooLarge {
        /// Which field.
        field: &'static str,
        /// What it holds.
        got: i64,
        /// The largest value allowed.
        limit: i64,
    },
    /// The caller asked for a stage that does not exist.
    /// The caller asked for a stage or probe that does not exist.
    UnknownStage {
        /// What was asked for.
        name: &'static str,
    },
}

/// Names are cheap to compare and a wrong name is a wrong column.
impl fmt::Display for BridgeError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            BridgeError::NullBuffer { field } => {
                write!(f, "buffer `{field}` was a null pointer")
            }
            BridgeError::ShapeMismatch {
                field,
                rows,
                cols,
                len,
            } => write!(
                f,
                "buffer `{field}` is declared {rows}x{cols} = {} elements but carries {len}",
                rows * cols
            ),
            BridgeError::IncompatibleShapes {
                a,
                a_shape,
                b,
                b_shape,
            } => write!(
                f,
                "shapes disagree: `{a}` is {}x{} and `{b}` is {}x{}",
                a_shape.0, a_shape.1, b_shape.0, b_shape.1
            ),
            BridgeError::LengthMismatch {
                field,
                expected,
                got,
            } => write!(
                f,
                "buffer `{field}` has {got} element(s) where {expected} were declared"
            ),
            BridgeError::MaskValue {
                field,
                index,
                value,
            } => write!(
                f,
                "mask `{field}` holds {value} at index {index}; it is a 0/1 mask"
            ),
            BridgeError::GroupIndex {
                index,
                value,
                n_levels,
            } => write!(
                f,
                "group index {value} at position {index} is outside 1..={n_levels} (0 means absent)"
            ),
            BridgeError::NameBuffer(field) => write!(
                f,
                "packed name buffer `{field}` is not a monotone partition of its byte range"
            ),
            BridgeError::NonFiniteScalar { field, value } => {
                write!(f, "scalar `{field}` must be finite, got {value}")
            }
            BridgeError::Args(msg) => write!(f, "{msg}"),
            BridgeError::TooLarge { field, got, limit } => write!(
                f,
                "`{field}` is {got}, which is larger than any result this package can \
                 produce ({limit})"
            ),
            BridgeError::UnknownStage { name } => write!(f, "unknown stage `{name}`"),
        }
    }
}

impl std::error::Error for BridgeError {}

/// A transport result. Every failure names the buffer and the two numbers that
/// disagreed.
pub type Result<T> = std::result::Result<T, BridgeError>;

// ---------------------------------------------------------------------------
// Count matrices: integer or real, never converted between
// ---------------------------------------------------------------------------

/// A count matrix as it arrives from R, in whichever of R's two numeric types it
/// happened to be.
///
/// The two are kept apart on purpose. `data` is integer for every microbiome
/// table in this repository, and widening it to `f64` would be numerically inert
/// but would change the type of the returned `feature_table`, which is part of the
/// result schema. `NA_integer_` is `INT_MIN` and must survive as such.
#[derive(Debug, Clone, PartialEq)]
pub enum Counts {
    /// An `INTSXP` count matrix. `NA_integer_` is `INT_MIN` and is carried as
    /// such; the returned `feature_table` is this table, so its integer-ness is
    /// part of the result schema.
    Int {
        /// `nrow`.
        rows: usize,
        /// `ncol`.
        cols: usize,
        /// Element order. R's matrices are column-major.
        layout: Layout,
        /// The payload, copied.
        values: Vec<i32>,
    },
    /// A `REALSXP` count matrix, copied bit for bit so that `NA_real_`, `NaN`,
    /// `+-Inf` and `-0.0` all survive.
    Real {
        /// `nrow`.
        rows: usize,
        /// `ncol`.
        cols: usize,
        /// Element order.
        layout: Layout,
        /// The payload, copied.
        values: Vec<f64>,
    },
}

impl Counts {
    /// `nrow`.
    pub fn rows(&self) -> usize {
        match self {
            Counts::Int { rows, .. } | Counts::Real { rows, .. } => *rows,
        }
    }

    /// `ncol`.
    pub fn cols(&self) -> usize {
        match self {
            Counts::Int { cols, .. } | Counts::Real { cols, .. } => *cols,
        }
    }

    /// `(nrow, ncol)`.
    pub fn shape(&self) -> (usize, usize) {
        (self.rows(), self.cols())
    }

    /// The declared element order.
    pub fn layout(&self) -> Layout {
        match self {
            Counts::Int { layout, .. } | Counts::Real { layout, .. } => *layout,
        }
    }

    /// Whether this arrived as `INTSXP`.
    pub fn is_int(&self) -> bool {
        matches!(self, Counts::Int { .. })
    }

    /// `rows * cols`.
    pub fn len(&self) -> usize {
        match self {
            Counts::Int { rows, cols, .. } | Counts::Real { rows, cols, .. } => rows * cols,
        }
    }

    /// Whether the table has no cells at all.
    /// Whether there are no names.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Element `(r, c)`, in the declared layout.
    pub fn at(&self, r: usize, c: usize) -> f64 {
        let i = match self.layout() {
            Layout::ColumnMajor => r + c * self.rows(),
            Layout::RowMajor => c + r * self.cols(),
        };
        match self {
            Counts::Int { values, .. } => {
                let v = values[i];
                if v == i32::MIN {
                    f64::NAN
                } else {
                    v as f64
                }
            }
            Counts::Real { values, .. } => values[i],
        }
    }

    /// The raw 64 bits of element `(r, c)`.
    ///
    /// This is what a transport claim is about. `at()` maps `NA_real_` and `NaN`
    /// onto the same `f64::NAN`, so it cannot tell them apart and is not what the
    /// round-trip test compares; `bits_at()` can, and does.
    pub fn bits_at(&self, r: usize, c: usize) -> u64 {
        let i = match self.layout() {
            Layout::ColumnMajor => r + c * self.rows(),
            Layout::RowMajor => c + r * self.cols(),
        };
        match self {
            Counts::Int { values, .. } => (values[i] as i64 as u64) & 0xFFFF_FFFF,
            Counts::Real { values, .. } => values[i].to_bits(),
        }
    }

    /// Copy into owned storage, bit for bit.
    ///
    /// For the integer case this is a widening to `i64` and is exact, because
    /// `i32 -> i64` cannot lose a bit and `NA_integer_` is `INT_MIN`. It is done
    /// eagerly and once, because the alternative -- handing `&[i32]` into a
    /// closure -- borrows from R.
    pub fn to_owned_f64(&self) -> Vec<f64> {
        match self {
            Counts::Int { values, .. } => values.iter().map(|&v| v as f64).collect(),
            Counts::Real { values, .. } => values.clone(),
        }
    }
}

// ---------------------------------------------------------------------------
// Scalar controls
// ---------------------------------------------------------------------------

/// The numeric controls `.ancombc2_core` takes, as sent.
///
/// Every one of these is a *finite* double on the original's side:
/// `quantile(x, s0_perc)` with a non-finite `s0_perc` raises
/// `'probs' outside [0,1]`, and `pseudo = Inf` would make `log(Inf)` `Inf`.
/// Requiring finiteness here means a bad control is rejected at the transport with
/// a named field, rather than becoming a `NaN` three stages later.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Controls {
    /// `pseudo`, added to the counts by the stage as `ancombc_prep.R:111` does,
    /// so what crosses the boundary is the *raw* count matrix.
    pub pseudo: f64,
    /// `s0_perc`, the `quantile()` probability for `s02`. Outside `[0, 1]` the
    /// reference raises `'probs' outside [0,1]`, so a non-finite value here is
    /// refused at the transport with the field named.
    pub s0_perc: f64,
    /// `prv_cut`, the prevalence filter: `prevalence >= prv_cut` keeps a taxon.
    pub prv_cut: f64,
    /// `lib_cut`, the library-size filter.
    pub lib_cut: f64,
    /// `alpha`, the significance level.
    pub alpha: f64,
    /// `iter_control$tol`, the MLE's `epsilon` threshold. Zero makes the loop run
    /// to `max_iter` instead of converging, which is how a fixed iteration count
    /// is obtained from the reference.
    pub iter_tol: f64,
    /// `iter_control$max_iter`. Zero runs no iterations and does **not** error.
    pub iter_max: i32,
    /// `em_control$tol`, the E-M threshold.
    pub em_tol: f64,
    /// `em_control$max_iter`. Zero **does** error: `.bias_em` only defines
    /// `delta_new` inside its loop, so the worker fails with
    /// `object 'delta_new' not found`.
    pub em_max: i32,
    /// `mdfdr_control$B`, the number of null replicates in the mdFDR correction.
    pub mdfdr_b: i32,
}

impl Controls {
    /// Reject a control that would otherwise become a `NaN` three stages later.
    pub fn validate(&self) -> Result<()> {
        for (field, v) in [
            ("pseudo", self.pseudo),
            ("s0_perc", self.s0_perc),
            ("prv_cut", self.prv_cut),
            ("lib_cut", self.lib_cut),
            ("alpha", self.alpha),
            ("iter_control$tol", self.iter_tol),
            ("em_control$tol", self.em_tol),
        ] {
            if !v.is_finite() {
                return Err(BridgeError::NonFiniteScalar { field, value: v });
            }
        }
        for (field, v) in [
            ("iter_control$max_iter", self.iter_max),
            ("em_control$max_iter", self.em_max),
            ("mdfdr_control$B", self.mdfdr_b),
        ] {
            if v < 0 {
                return Err(BridgeError::MaskValue {
                    field,
                    index: 0,
                    value: v as i64,
                });
            }
        }
        Ok(())
    }
}

impl Default for Controls {
    /// The pinned original's own defaults, from `ancombc2()`'s signature and
    /// `ancombc_prep.R:92-108`. Not `AncombcConfig::default()`: this is the
    /// transport's view and the two drifting apart is exactly the bug the
    /// `defaults_match_the_original` test exists to catch.
    fn default() -> Self {
        Controls {
            pseudo: 0.0,
            s0_perc: 0.05,
            prv_cut: 0.10,
            lib_cut: 0.0,
            alpha: 0.05,
            iter_tol: 0.01,
            iter_max: 20,
            em_tol: 1e-5,
            em_max: 100,
            mdfdr_b: 100,
        }
    }
}

// ---------------------------------------------------------------------------
// Packed names
// ---------------------------------------------------------------------------

/// Names as they cross the boundary: one packed byte range plus offsets.
///
/// R's `STRSXP` is a vector of `CHARSXP`, and reaching into it from Rust means
/// dereferencing R's internal pointer table. That layout is stable in practice but
/// it is R's business, not this crate's. So the C shim -- which does speak R --
/// packs the strings into a flat buffer with `i64` offsets, and this crate only
/// ever sees bytes and integers. `unpack` below is the same partition read back,
/// and it is where a malformed buffer is caught.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PackedNames {
    /// `n + 1` offsets into `bytes`, starting at 0.
    offsets: Vec<i64>,
    /// The concatenated UTF-8.
    bytes: Vec<u8>,
}

impl PackedNames {
    /// Build from a partition, rejecting anything that is not one.
    ///
    /// Rejecting here rather than at the point of use is what stops a truncated
    /// offsets array from reading past the byte buffer.
    pub fn new(offsets: &[i64], bytes: &[u8]) -> Result<Self> {
        if offsets.is_empty() {
            return Err(BridgeError::NameBuffer("offsets"));
        }
        let n = bytes.len() as i64;
        if offsets[0] != 0 {
            return Err(BridgeError::NameBuffer("offsets"));
        }
        if *offsets.last().unwrap() != n {
            return Err(BridgeError::NameBuffer("offsets"));
        }
        if offsets.windows(2).any(|w| w[1] < w[0]) {
            return Err(BridgeError::NameBuffer("offsets"));
        }
        Ok(PackedNames {
            offsets: offsets.to_vec(),
            bytes: bytes.to_vec(),
        })
    }

    /// `rows * cols`.
    /// The number of names.
    pub fn len(&self) -> usize {
        self.offsets.len().saturating_sub(1)
    }

    /// Whether the table has no cells at all.
    /// Whether there are no names.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// String `i`, losslessly, or `None` if it is not valid UTF-8.
    ///
    /// R strings carry an encoding flag and `mkCharCE` is not this crate's
    /// business; `R CMD check` would flag a non-UTF-8 fixture long before a name
    /// mattered. An invalid byte range is a named failure rather than a lossy
    /// replacement character.
    pub fn get(&self, i: usize) -> Option<std::borrow::Cow<'_, str>> {
        let (lo, hi) = (*self.offsets.get(i)?, *self.offsets.get(i + 1)?);
        if lo < 0 || hi < lo || hi as usize > self.bytes.len() {
            return None;
        }
        std::str::from_utf8(&self.bytes[lo as usize..hi as usize])
            .ok()
            .map(std::borrow::Cow::Borrowed)
    }

    /// Pack a slice of strings. Used by the tests and by the C shim's own
    /// round-trip check.
    pub fn pack<I, S>(names: I) -> Self
    where
        I: IntoIterator<Item = S>,
        S: AsRef<str>,
    {
        let mut offsets = vec![0i64];
        let mut bytes = Vec::new();
        for s in names {
            bytes.extend_from_slice(s.as_ref().as_bytes());
            offsets.push(bytes.len() as i64);
        }
        PackedNames { offsets, bytes }
    }
}

// ---------------------------------------------------------------------------
// Validation of one transport call
// ---------------------------------------------------------------------------

/// Everything one `.ancombc2_core` call is given, already owned and already
/// checked.
#[derive(Debug, Clone, PartialEq)]
pub struct Request {
    /// `O1 + pseudo`'s base: `data`, as it arrived.
    pub data: Counts,
    /// `O2 + pseudo`'s base: `aggregate_data`.
    pub aggregate: Counts,
    /// The design, `n_samp x p`, column-major: `(nrow, ncol, layout, values)`.
    pub design: (usize, usize, Layout, Vec<f64>),
    /// `stats::complete.cases(x)` over the design: one byte per sample, 0 or 1.
    pub design_complete: Vec<u8>,
    /// Group level per sample, 1-based; `0` means the sample has no group.
    pub group_index: Vec<i32>,
    /// `levels(group)`, so a level index can be turned back into a label.
    pub group_labels: PackedNames,
    /// `colnames(x)`.
    pub fix_eff: PackedNames,
    /// `rownames(data)`.
    pub taxon_names: PackedNames,
    /// `colnames(data)`.
    pub sample_names: PackedNames,
    /// The numeric controls.
    pub controls: Controls,
    /// The logical switches: `global`, `pairwise`,
    /// `dunnet`, `trend`, `pseudo_sens`, `conservative`, `struc_zero`, `neg_lb`.
    pub flags: Flags,
}

/// The logical switches, kept out of [`Controls`] because they are not numbers
/// and a `f64` holding 0 or 1 is a worse way to say "no".
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Flags {
    /// `global`, the LRT/Wald global test.
    pub global: bool,
    /// `pairwise`, the mdFDR pairwise comparisons.
    pub pairwise: bool,
    /// `dunnet`, Dunnett's test. Retained R in the replacement; carried so the
    /// stage can refuse rather than silently ignore it.
    pub dunnet: bool,
    /// `trend`, the directional test. As above.
    pub trend: bool,
    /// `pseudo_sens`. **Defaults to TRUE**, matching `ancombc2.R:401`. The
    /// scaffold package had this wrong once (`pseudo_sens = FALSE`), and a wrong
    /// default on the bridge would hide that rather than expose it.
    pub pseudo_sens: bool,
    /// `conservative`, whether the pseudo-count grid is `{0, 0.1, 0.5, 1}` or
    /// `seq(0.01, 0.5, 0.01)`.
    pub conservative: bool,
    /// `struc_zero`, whether flagged taxa are removed from the primary fit.
    pub struc_zero: bool,
    /// `neg_lb`, which widens the structural-zero net beyond a strict zero.
    pub neg_lb: bool,
}

impl Default for Flags {
    fn default() -> Self {
        // `ancombc2()`'s own defaults. `pseudo_sens = TRUE` is the one that
        // matters here: the scaffold package got it wrong once, and a bridge that
        // defaults it to FALSE would hide that behind a default.
        Flags {
            global: false,
            pairwise: false,
            dunnet: false,
            trend: false,
            pseudo_sens: true,
            conservative: true,
            struc_zero: false,
            neg_lb: false,
        }
    }
}

impl Default for Request {
    fn default() -> Self {
        Request {
            data: Counts::Int {
                rows: 0,
                cols: 0,
                layout: Layout::ColumnMajor,
                values: Vec::new(),
            },
            aggregate: Counts::Int {
                rows: 0,
                cols: 0,
                layout: Layout::ColumnMajor,
                values: Vec::new(),
            },
            design: (0, 0, Layout::ColumnMajor, Vec::new()),
            design_complete: Vec::new(),
            group_index: Vec::new(),
            group_labels: PackedNames::pack(Vec::<String>::new()),
            fix_eff: PackedNames::pack(Vec::<String>::new()),
            taxon_names: PackedNames::pack(Vec::<String>::new()),
            sample_names: PackedNames::pack(Vec::<String>::new()),
            controls: Controls::default(),
            flags: Flags::default(),
        }
    }
}

impl Request {
    /// Every shape, length and domain check, before any value is read.
    ///
    /// Ordered so the error a caller sees is about the *first* thing that is
    /// wrong rather than about a downstream consequence of it.
    pub fn validate(&self) -> Result<()> {
        let (n_tax, n_samp) = self.data.shape();
        if n_tax == 0 || n_samp == 0 {
            return Err(BridgeError::ShapeMismatch {
                field: "data",
                rows: n_tax,
                cols: n_samp,
                len: self.data.len(),
            });
        }
        if self.aggregate.shape() != (n_tax, n_samp) {
            return Err(BridgeError::IncompatibleShapes {
                a: "data",
                a_shape: (n_tax, n_samp),
                b: "aggregate_data",
                b_shape: self.aggregate.shape(),
            });
        }
        // `.ancombc2_core` builds O1 and O2 independently and then requires
        // nothing about their agreement beyond the taxa it looks up by name, but
        // the returned `feature_table` is O2 and the retained taxa are indices
        // into O1, so a mismatch in orientation would surface as a silent
        // transposition much later. Refuse it here.
        if self.aggregate.layout() != self.data.layout() {
            return Err(BridgeError::IncompatibleShapes {
                a: "data",
                a_shape: (n_tax, n_samp),
                b: "aggregate_data",
                b_shape: self.aggregate.shape(),
            });
        }

        let (p_rows, p_cols, _, design) = &self.design;
        if *p_rows != n_samp {
            return Err(BridgeError::IncompatibleShapes {
                a: "data",
                a_shape: (n_tax, n_samp),
                b: "design",
                b_shape: (*p_rows, *p_cols),
            });
        }
        if *p_cols == 0 {
            return Err(BridgeError::ShapeMismatch {
                field: "design",
                rows: *p_rows,
                cols: *p_cols,
                len: design.len(),
            });
        }
        if design.len() != p_rows * p_cols {
            return Err(BridgeError::ShapeMismatch {
                field: "design",
                rows: *p_rows,
                cols: *p_cols,
                len: design.len(),
            });
        }

        if self.design_complete.len() != n_samp {
            return Err(BridgeError::LengthMismatch {
                field: "design_complete",
                expected: n_samp,
                got: self.design_complete.len(),
            });
        }
        for (i, &v) in self.design_complete.iter().enumerate() {
            if v > 1 {
                return Err(BridgeError::MaskValue {
                    field: "design_complete",
                    index: i,
                    value: v as i64,
                });
            }
        }

        let n_levels = self.group_labels.len();
        if self.group_index.len() != n_samp {
            return Err(BridgeError::LengthMismatch {
                field: "group_index",
                expected: n_samp,
                got: self.group_index.len(),
            });
        }
        for (i, &v) in self.group_index.iter().enumerate() {
            // 0 is the original's own "no group" marker: the candidate's R side
            // maps `NULL` group to it rather than inventing a level.
            if v < 0 || (v != 0 && v as usize > n_levels) {
                return Err(BridgeError::GroupIndex {
                    index: i,
                    value: v,
                    n_levels,
                });
            }
        }

        if self.fix_eff.len() != *p_cols {
            return Err(BridgeError::LengthMismatch {
                field: "fix_eff",
                expected: *p_cols,
                got: self.fix_eff.len(),
            });
        }
        if self.sample_names.len() != n_samp {
            return Err(BridgeError::LengthMismatch {
                field: "sample_names",
                expected: n_samp,
                got: self.sample_names.len(),
            });
        }
        // `rownames(O2)` must have one entry per taxon or the returned table has
        // a different number of rows than taxa, which no later check would catch.
        if self.taxon_names.len() != n_tax {
            return Err(BridgeError::LengthMismatch {
                field: "taxon_names",
                expected: n_tax,
                got: self.taxon_names.len(),
            });
        }

        self.controls.validate()
    }

    /// Names that came through as valid UTF-8. A name that did not is a
    /// transport failure, not something to drop: dropping it would shorten a
    /// column and shift every downstream index.
    pub fn check_names(&self) -> Result<()> {
        for (field, names, expected) in [
            ("fix_eff", &self.fix_eff, self.design.1),
            ("sample_names", &self.sample_names, self.data.cols()),
            ("taxon_names", &self.taxon_names, self.data.rows()),
            ("group_labels", &self.group_labels, self.group_labels.len()),
        ] {
            if names.len() != expected {
                return Err(BridgeError::LengthMismatch {
                    field,
                    expected,
                    got: names.len(),
                });
            }
            for i in 0..names.len() {
                if names.get(i).is_none() {
                    return Err(BridgeError::NameBuffer(field));
                }
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn real_matrix(rows: usize, cols: usize, v: Vec<f64>) -> Counts {
        Counts::Real {
            rows,
            cols,
            layout: Layout::ColumnMajor,
            values: v,
        }
    }

    fn good() -> Request {
        let n_tax = 4;
        let n_samp = 6;
        let p = 3;
        let mut r = Request::default();
        r.data = real_matrix(
            n_tax,
            n_samp,
            (0..n_tax * n_samp).map(|k| k as f64).collect(),
        );
        r.aggregate = r.data.clone();
        r.design = (
            n_samp,
            p,
            Layout::ColumnMajor,
            (0..n_samp * p).map(|k| (k % 3) as f64).collect(),
        );
        r.design_complete = vec![1; n_samp];
        r.group_index = vec![1, 1, 2, 2, 3, 3];
        r.group_labels = PackedNames::pack(["a", "b", "c"]);
        r.fix_eff = PackedNames::pack(["(Intercept)", "g2", "g3"]);
        r.taxon_names = PackedNames::pack(["T1", "T2", "T3", "T4"]);
        r.sample_names = PackedNames::pack(["S1", "S2", "S3", "S4", "S5", "S6"]);
        r
    }

    #[test]
    fn a_well_formed_request_validates() {
        let r = good();
        r.validate().unwrap();
        r.check_names().unwrap();
    }

    #[test]
    fn the_shape_error_names_the_two_numbers_that_disagree() {
        let mut r = good();
        r.aggregate = real_matrix(3, 6, vec![0.0; 18]);
        match r.validate() {
            Err(BridgeError::IncompatibleShapes {
                a,
                a_shape,
                b,
                b_shape,
            }) => {
                assert_eq!(a, "data");
                assert_eq!(a_shape, (4, 6));
                assert_eq!(b, "aggregate_data");
                assert_eq!(b_shape, (3, 6));
            }
            other => panic!("{other:?}"),
        }
    }

    #[test]
    fn a_short_buffer_is_a_named_error_not_a_read_past_the_end() {
        let mut r = good();
        r.design.3.pop();
        assert!(matches!(
            r.validate(),
            Err(BridgeError::ShapeMismatch {
                field: "design",
                ..
            })
        ));
    }

    #[test]
    fn a_non_binary_mask_is_rejected() {
        let mut r = good();
        r.design_complete[2] = 2;
        assert!(matches!(
            r.validate(),
            Err(BridgeError::MaskValue {
                field: "design_complete",
                index: 2,
                value: 2
            })
        ));
    }

    #[test]
    fn a_group_index_outside_the_levels_is_rejected_and_zero_means_absent() {
        let mut r = good();
        r.group_index[1] = 0; // absent: legal
        r.validate().unwrap();
        r.group_index[1] = 4; // one past the last level
        assert!(matches!(
            r.validate(),
            Err(BridgeError::GroupIndex {
                index: 1,
                value: 4,
                n_levels: 3
            })
        ));
        r.group_index[1] = -1;
        assert!(matches!(r.validate(), Err(BridgeError::GroupIndex { .. })));
    }

    #[test]
    fn a_name_count_that_disagrees_with_its_buffer_is_rejected() {
        let mut r = good();
        r.taxon_names = PackedNames::pack(["T1", "T2"]);
        assert!(matches!(
            r.validate(),
            Err(BridgeError::LengthMismatch {
                field: "taxon_names",
                ..
            })
        ));
    }

    #[test]
    fn a_non_finite_control_is_named_rather_than_becoming_a_nan() {
        let mut r = good();
        r.controls.s0_perc = f64::NAN;
        assert!(matches!(
            r.controls.validate(),
            Err(BridgeError::NonFiniteScalar {
                field: "s0_perc",
                ..
            })
        ));
        r.controls.s0_perc = f64::INFINITY;
        assert!(matches!(
            r.controls.validate(),
            Err(BridgeError::NonFiniteScalar {
                field: "s0_perc",
                ..
            })
        ));
    }

    #[test]
    fn a_name_buffer_must_be_a_monotone_partition() {
        assert!(PackedNames::new(&[0, 3, 3], b"abc").is_ok());
        assert!(PackedNames::new(&[0, 3], b"abcd").is_err());
        assert!(PackedNames::new(&[1, 3], b"abc").is_err());
        assert!(PackedNames::new(&[0, 3, 2], b"abc").is_err());
        assert!(PackedNames::new(&[], b"").is_err());
    }

    #[test]
    fn names_round_trip_through_the_packed_form() {
        let names = ["a", "", "with space", "unicod\u{e9}", "tab\there"];
        let p = PackedNames::pack(names);
        assert_eq!(p.len(), names.len());
        for (i, n) in names.iter().enumerate() {
            assert_eq!(p.get(i).unwrap(), *n);
        }
        assert!(p.get(names.len()).is_none());
    }

    #[test]
    fn column_major_is_the_layout_and_element_lookup_follows_it() {
        // 2x3 column-major: (0,0) (1,0) (0,1) (1,1) (0,2) (1,2)
        let m = real_matrix(2, 3, vec![10.0, 11.0, 20.0, 21.0, 30.0, 31.0]);
        assert_eq!(m.at(0, 0), 10.0);
        assert_eq!(m.at(1, 0), 11.0);
        assert_eq!(m.at(0, 1), 20.0);
        assert_eq!(m.at(1, 2), 31.0);
        // The same bytes read row-major put 11.0 at (0,1) instead of 20.0. This
        // is the assertion that matters: a transposed design does not fail
        // loudly, it produces a *different but plausible* design, so the layout
        // has to be checked rather than assumed. (The old JSON bridge's layout
        // mismatch surfaced as "estimation failed for covariate col0", which is
        // the same class of misattribution.)
        let row_major = Counts::Real {
            rows: 2,
            cols: 3,
            layout: Layout::RowMajor,
            values: m.layout_free_clone(),
        };
        assert_eq!(m.at(0, 1), 20.0);
        assert_eq!(row_major.at(0, 1), 11.0);
        assert_ne!(m.at(0, 1), row_major.at(0, 1));
    }

    impl Counts {
        fn layout_free_clone(&self) -> Vec<f64> {
            match self {
                Counts::Real { values, .. } => values.clone(),
                Counts::Int { values, .. } => values.iter().map(|&v| v as f64).collect(),
            }
        }
    }

    #[test]
    fn na_integer_survives_as_int_min_and_is_not_the_same_as_zero() {
        let m = Counts::Int {
            rows: 1,
            cols: 3,
            layout: Layout::ColumnMajor,
            values: vec![i32::MIN, 0, 1],
        };
        assert!(m.at(0, 0).is_nan());
        assert_eq!(m.at(0, 1), 0.0);
        assert_eq!(m.bits_at(0, 0), (i32::MIN as i64 as u64) & 0xFFFF_FFFF);
    }
}
