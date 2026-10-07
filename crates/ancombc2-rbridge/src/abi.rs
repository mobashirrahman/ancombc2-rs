//! The C ABI. Every `extern "C"` entry point lives here, and nowhere else.
//!
//! # Rules this layer enforces
//!
//! 1. **No unwinding across the boundary.** Rust panics cannot propagate into R;
//!    `extern "C"` functions abort the process instead of unwinding, which would
//!    take the R session with it. Every entry point is wrapped in
//!    [`catch_unwind`] and turns a caught panic into a status code plus a
//!    message, which `init.c` raises with `Rf_error`. A Rust panic inside a
//!    Rayon worker is caught by Rayon and re-raised here, so it takes the same
//!    path.
//! 2. **No allocation crosses the boundary.** Every buffer is allocated by R, in
//!    `init.c`, and R frees it. Rust writes into it and never frees it. The only
//!    thing Rust hands back that R has to release is nothing: the error message is
//!    copied into a caller-provided buffer.
//! 3. **No R pointer is stored.** Nothing here keeps a `*const` beyond the call.
//!    [`crate::owned::StageInput`] is built, used, and dropped inside one call.
//! 4. **The last error is a thread-local string**, read back through a
//!    length-then-copy pair so there is no allocation and no lifetime question.
//!    `init.c` only reads it after a non-zero status, and it is overwritten by the
//!    next call, which is documented rather than defended against.
//!
//! # Status codes
//!
//! They are deliberately distinguishable. "the arguments were wrong", "a stage
//! failed", and "the bridge panicked" are different bugs and a caller that cannot
//! tell them apart will eventually treat a panic as a validation error.

use std::cell::RefCell;
use std::panic::{catch_unwind, AssertUnwindSafe};

use crate::core_output_plan;
use crate::output::{emit_core_payloads, CorePayloads, Emitted, RBackedReductions};
use crate::owned::{Owned, StageInput, TransportValue};
use crate::transport::{BridgeError, Controls, Counts, Flags, Layout, PackedNames, Request};
use ancombc2_core::matrix::Matrix;
use ancombc2_core::mle::{lm_fit_all, DesignCache};
use ancombc2_core::preprocess::CountMatrix;

/// Everything went well.
pub const RB_OK: i32 = 0;
/// The caller's types, lengths, shapes or domains are wrong. Nothing was read.
pub const RB_ERR_ARGS: i32 = 1;
/// A numerical stage refused: rank deficiency, a zero variance, a `group` that
/// does not match. The message names the taxon or covariate where it can.
pub const RB_ERR_STAGE: i32 = 2;
/// A Rust panic was caught at the boundary. Nothing partial was written.
pub const RB_ERR_PANIC: i32 = 3;
/// An argument was well-formed but the requested stage does not exist.
pub const RB_ERR_STAGE_UNKNOWN: i32 = 4;

thread_local! {
    static LAST_ERROR: RefCell<String> = const { RefCell::new(String::new()) };
}

fn set_error(msg: impl Into<String>) {
    LAST_ERROR.with(|e| *e.borrow_mut() = msg.into());
}

/// Byte length of the last error message, **excluding** the NUL.
///
/// Read this first, allocate, then call [`ancombc2_rb_copy_last_error`]. Two calls
/// and no allocation is a little more code than returning a `malloc`'d string, and
/// it removes the only ownership question the boundary would otherwise have.
///
/// # Safety
///
/// Safe to call: it touches no caller memory. It is `unsafe` only for symmetry
/// with the rest of the surface, which is declared `unsafe` as a whole because a
/// caller cannot be trusted to have sized its buffers correctly elsewhere.
#[no_mangle]
pub unsafe extern "C" fn ancombc2_rb_last_error_len() -> i64 {
    LAST_ERROR.with(|e| e.borrow().len() as i64)
}

/// Copy the last error into `dst`, NUL-terminated.
///
/// Returns the number of bytes written **excluding** the NUL. If `cap` is too
/// small the message is truncated and the last byte written is `'\u{2026}'`, so a
/// truncated diagnostic is visibly truncated rather than quietly wrong.
///
/// # Safety
///
/// `dst` must be null, or valid for `cap` bytes. A null `dst` writes nothing and
/// returns the length.
#[no_mangle]
pub unsafe extern "C" fn ancombc2_rb_copy_last_error(dst: *mut u8, cap: i64) -> i64 {
    if dst.is_null() {
        return LAST_ERROR.with(|e| e.borrow().len() as i64);
    }
    if cap <= 0 {
        return 0;
    }
    let msg = LAST_ERROR.with(|e| e.borrow().clone());
    let bytes = msg.as_bytes();
    let cap = cap as usize;
    let n = bytes.len().min(cap.saturating_sub(1));
    std::ptr::copy_nonoverlapping(bytes.as_ptr(), dst, n);
    if n < bytes.len() && n < cap - 1 {
        // room for the ellipsis, so a truncated message says so
        let ell = "\u{2026}".as_bytes();
        let take = ell.len().min(cap - 1 - n);
        std::ptr::copy_nonoverlapping(ell.as_ptr(), dst.add(n), take);
        std::ptr::write(dst.add(n + take), 0);
        return (n + take) as i64;
    }
    std::ptr::write(dst.add(n), 0);
    n as i64
}

/// The compatibility target this bridge is written against.
///
/// Two-call protocol: call with `out = NULL` to learn the length, then call again
/// with a buffer of at least `length + 1`. Returning the length for a null
/// pointer is what makes that work without an allocation.
///
/// # Safety
///
/// `out` must be null, or valid for `cap` bytes. A null `out` writes nothing and
/// returns the length.
#[no_mangle]
pub unsafe extern "C" fn ancombc2_rb_oracle_sha(out: *mut u8, cap: i64) -> i64 {
    const SHA: &[u8] = b"dc4febdf59badb3a8dfe0c767ef2186323c2199a";
    unsafe { copy_out(SHA, out, cap) }
}

/// The bridge's own version and the pinned version it claims.
///
/// # Safety
///
/// As [`ancombc2_rb_oracle_sha`].
#[no_mangle]
pub unsafe extern "C" fn ancombc2_rb_version(out: *mut u8, cap: i64) -> i64 {
    const V: &[u8] = concat!(
        "ancombc2-rbridge 0.1.0; transport for ANCOMBC 2.15.2 @ ",
        "dc4febdf59badb3a8dfe0c767ef2186323c2199a"
    )
    .as_bytes();
    unsafe { copy_out(V, out, cap) }
}

/// Copy `src` into `out`, NUL-terminated, and return the number of bytes written
/// excluding the NUL.
///
/// A null `out` is a length query and returns `src.len()` without writing. The
/// two-call protocol is what keeps this allocation-free: the C shim asks for the
/// length with `R_alloc`, then reads.
unsafe fn copy_out(src: &[u8], out: *mut u8, cap: i64) -> i64 {
    if out.is_null() {
        return src.len() as i64;
    }
    if cap <= 0 {
        return 0;
    }
    let cap = cap as usize;
    let n = src.len().min(cap - 1);
    std::ptr::copy_nonoverlapping(src.as_ptr(), out, n);
    std::ptr::write(out.add(n), 0);
    n as i64
}

// ---------------------------------------------------------------------------
// The transport probe
// ---------------------------------------------------------------------------

/// Raw borrowed description of one call, as it arrives from `init.c`.
///
/// Everything is a pointer plus a length. No R types appear: `init.c` has already
/// checked them and extracted the dimensions, and this struct is the *validated*
/// view of the call.
#[repr(C)]
pub struct RawRequest {
    // ---- the count matrix, `data` ---------------------------------------
    /// `INTSXP` payload of `data`, or null when `data_is_int` is 0.
    pub data_int: *const i32,
    /// `REALSXP` payload of `data`, unused when `data_is_int` is 1.
    pub data_real: *const f64,
    /// `nrow(data)`.
    pub data_rows: i64,
    /// `ncol(data)`.
    pub data_cols: i64,
    /// 1 when `data` is `INTSXP`. The type is carried rather than coerced
    /// because `NA_integer_` is `INT_MIN` and the returned `feature_table` is
    /// `O2` with whatever type it arrived as.
    pub data_is_int: i32,

    // ---- the aggregate table, `aggregate_data` ---------------------------
    /// `INTSXP` payload of `aggregate_data`.
    pub aggregate_int: *const i32,
    /// `REALSXP` payload of `aggregate_data`.
    pub aggregate_real: *const f64,
    /// `nrow(aggregate_data)`.
    pub aggregate_rows: i64,
    /// `ncol(aggregate_data)`.
    pub aggregate_cols: i64,
    /// 1 when `aggregate_data` is `INTSXP`.
    pub aggregate_is_int: i32,

    // ---- the design, `x = model.matrix(~ fix_formula, meta_data)` ---------
    /// `REALSXP` payload of `x`, **column-major**, which is R's own order.
    pub design: *const f64,
    /// `nrow(x)`, the number of samples.
    pub design_rows: i64,
    /// `ncol(x)`, the number of fixed-effect columns.
    pub design_cols: i64,

    // ---- derived per-sample metadata --------------------------------------
    /// `stats::complete.cases(x)`, one byte per sample, 0 or 1. Carried as bytes
    /// rather than as `LGLSXP` because R's logical is an `int`, and a byte is
    /// what a mask is.
    pub design_complete: *const u8,
    /// Group level per sample, **1-based**; `0` means the sample has no group.
    /// 1-based because R's factor levels are 1-based and a 0-based index that
    /// meant "absent" would collide with the first level.
    pub group_index: *const i32,
    /// `length(levels(group))`; 0 when there is no group.
    pub n_group_levels: i64,

    // ---- packed names: `offsets` has `n + 1` entries ----------------------
    /// Offsets into `group_labels_bytes`, `n_group_labels + 1` of them.
    pub group_labels_offsets: *const i64,
    /// The label bytes.
    pub group_labels_bytes: *const u8,
    /// `n_group_labels`.
    pub n_group_labels: i64,

    /// Offsets into `fix_eff_bytes`, `n_fix_eff + 1` of them.
    pub fix_eff_offsets: *const i64,
    /// `colnames(x)`.
    pub fix_eff_bytes: *const u8,
    /// `ncol(x)`.
    pub n_fix_eff: i64,

    /// Offsets into `taxon_names_bytes`, `n_taxon_names + 1` of them.
    pub taxon_names_offsets: *const i64,
    /// `rownames(data)`.
    pub taxon_names_bytes: *const u8,
    /// `nrow(data)`.
    pub n_taxon_names: i64,

    /// Offsets into `sample_names_bytes`, `n_sample_names + 1` of them.
    pub sample_names_offsets: *const i64,
    /// `colnames(data)`.
    pub sample_names_bytes: *const u8,
    /// `ncol(data)`.
    pub n_sample_names: i64,

    // ---- numeric controls, from the ancombc2 formals ----------------------
    /// `pseudo`. Already added to the counts on the R side? No: the original adds
    /// it inside `.ancombc2_core` (`ancombc_prep.R:111`), and so does the stage,
    /// so what crosses here is the *raw* count matrix.
    pub pseudo: f64,
    /// `s0_perc`.
    pub s0_perc: f64,
    /// `prv_cut`.
    pub prv_cut: f64,
    /// `lib_cut`.
    pub lib_cut: f64,
    /// `alpha`.
    pub alpha: f64,
    /// `iter_control$tol`.
    pub iter_tol: f64,
    /// `iter_control$max_iter`.
    pub iter_max: i32,
    /// `em_control$tol`.
    pub em_tol: f64,
    /// `em_control$max_iter`.
    pub em_max: i32,
    /// `mdfdr_control$B`.
    pub mdfdr_b: i32,

    // ---- logical switches -------------------------------------------------
    /// `global`.
    pub global: i32,
    /// `pairwise`.
    pub pairwise: i32,
    /// `dunnet`.
    pub dunnet: i32,
    /// `trend`.
    pub trend: i32,
    /// `pseudo_sens`. The original's default is **TRUE**
    /// (`ancombc2.R:401`); the R side sends the value it actually holds rather
    /// than relying on a default here, because a default on the bridge that
    /// disagreed with the package's would be invisible.
    pub pseudo_sens: i32,
    /// `conservative`.
    pub conservative: i32,
    /// `struc_zero`.
    pub struc_zero: i32,
    /// `neg_lb`.
    pub neg_lb: i32,
}

fn as_usize(v: i64, field: &'static str) -> Result<usize, BridgeError> {
    if v < 0 {
        return Err(BridgeError::LengthMismatch {
            field,
            expected: 0,
            got: v as usize,
        });
    }
    Ok(v as usize)
}

unsafe fn slice<'a, T>(p: *const T, n: usize, field: &'static str) -> Result<&'a [T], BridgeError> {
    if p.is_null() {
        return Err(BridgeError::NullBuffer { field });
    }
    Ok(std::slice::from_raw_parts(p, n))
}

unsafe fn packed(
    offsets: *const i64,
    bytes: *const u8,
    n_names: usize,
    field: &'static str,
) -> Result<PackedNames, BridgeError> {
    let offs = slice::<i64>(offsets, n_names + 1, field)?;
    let last = *offs.last().ok_or(BridgeError::NameBuffer(field))?;
    if last < 0 {
        return Err(BridgeError::NameBuffer(field));
    }
    let buf = slice::<u8>(bytes, last as usize, field)?;
    PackedNames::new(offs, buf)
}

unsafe fn counts(
    int_p: *const i32,
    real_p: *const f64,
    rows: i64,
    cols: i64,
    is_int: i32,
    field: &'static str,
) -> Result<Counts, BridgeError> {
    let rows = as_usize(rows, field)?;
    let cols = as_usize(cols, field)?;
    let len = rows.checked_mul(cols).ok_or(BridgeError::ShapeMismatch {
        field,
        rows,
        cols,
        len: usize::MAX,
    })?;
    if is_int != 0 {
        let v = slice::<i32>(int_p, len, field)?;
        Ok(Counts::Int {
            rows,
            cols,
            layout: Layout::ColumnMajor,
            values: v.to_vec(),
        })
    } else {
        let v = slice::<f64>(real_p, len, field)?;
        Ok(Counts::Real {
            rows,
            cols,
            layout: Layout::ColumnMajor,
            values: v.to_vec(),
        })
    }
}

unsafe fn to_request(r: &RawRequest) -> Result<Request, BridgeError> {
    let n_levels = as_usize(r.n_group_levels, "group_levels")?;
    let design_rows = as_usize(r.design_rows, "design")?;
    let design_cols = as_usize(r.design_cols, "design")?;

    let req = Request {
        data: counts(
            r.data_int,
            r.data_real,
            r.data_rows,
            r.data_cols,
            r.data_is_int,
            "data",
        )?,
        aggregate: counts(
            r.aggregate_int,
            r.aggregate_real,
            r.aggregate_rows,
            r.aggregate_cols,
            r.aggregate_is_int,
            "aggregate_data",
        )?,
        design: (
            design_rows,
            design_cols,
            Layout::ColumnMajor,
            slice::<f64>(r.design, design_rows * design_cols, "design")?.to_vec(),
        ),
        design_complete: slice::<u8>(
            r.design_complete,
            as_usize(r.data_cols, "design_complete")?,
            "design_complete",
        )?
        .to_vec(),
        group_index: slice::<i32>(
            r.group_index,
            as_usize(r.data_cols, "group_index")?,
            "group_index",
        )?
        .to_vec(),
        group_labels: packed(
            r.group_labels_offsets,
            r.group_labels_bytes,
            n_levels,
            "group_labels",
        )?,
        fix_eff: packed(
            r.fix_eff_offsets,
            r.fix_eff_bytes,
            as_usize(r.n_fix_eff, "fix_eff")?,
            "fix_eff",
        )?,
        taxon_names: packed(
            r.taxon_names_offsets,
            r.taxon_names_bytes,
            as_usize(r.n_taxon_names, "taxon_names")?,
            "taxon_names",
        )?,
        sample_names: packed(
            r.sample_names_offsets,
            r.sample_names_bytes,
            as_usize(r.n_sample_names, "sample_names")?,
            "sample_names",
        )?,
        controls: Controls {
            pseudo: r.pseudo,
            s0_perc: r.s0_perc,
            prv_cut: r.prv_cut,
            lib_cut: r.lib_cut,
            alpha: r.alpha,
            iter_tol: r.iter_tol,
            iter_max: r.iter_max,
            em_tol: r.em_tol,
            em_max: r.em_max,
            mdfdr_b: r.mdfdr_b,
        },
        flags: Flags {
            global: r.global != 0,
            pairwise: r.pairwise != 0,
            dunnet: r.dunnet != 0,
            trend: r.trend != 0,
            pseudo_sens: r.pseudo_sens != 0,
            conservative: r.conservative != 0,
            struc_zero: r.struc_zero != 0,
            neg_lb: r.neg_lb != 0,
        },
    };
    req.validate()?;
    req.check_names()?;
    Ok(req)
}

/// Build the owned [`StageInput`] a parallel stage is allowed to take.
pub fn stage_input(req: &Request) -> StageInput {
    let (values_int, values_real) = match &req.data {
        Counts::Int { values, .. } => (
            Some(
                values
                    .iter()
                    .map(|&v| {
                        if v == i32::MIN {
                            TransportValue::Missing
                        } else {
                            TransportValue::Count(v as i64)
                        }
                    })
                    .collect::<Vec<_>>(),
            ),
            None,
        ),
        Counts::Real { values, .. } => (
            None,
            Some(
                values
                    .iter()
                    .map(|&v| TransportValue::Count(v.to_bits() as i64))
                    .collect::<Vec<_>>(),
            ),
        ),
    };
    // A real count matrix is carried as raw bits rather than as a double, so that
    // `NA_real_`, `NaN`, `-0.0` and the infinities all survive the trip into the
    // stage. The reinterpretation is exact in both directions.
    let data = match (values_int, values_real) {
        (Some(v), _) => Owned::from_owned(v),
        (_, Some(v)) => Owned::from_owned(v),
        _ => Owned::from_owned(Vec::new()),
    };
    let aggregate = match &req.aggregate {
        Counts::Int { values, .. } => Owned::from_owned(
            values
                .iter()
                .map(|&v| {
                    if v == i32::MIN {
                        TransportValue::Missing
                    } else {
                        TransportValue::Count(v as i64)
                    }
                })
                .collect::<Vec<_>>(),
        ),
        Counts::Real { values, .. } => Owned::from_owned(
            values
                .iter()
                .map(|&v| TransportValue::Count(v.to_bits() as i64))
                .collect::<Vec<_>>(),
        ),
    };
    StageInput {
        data,
        aggregate,
        n_tax: req.data.rows(),
        n_samp: req.data.cols(),
        design: Owned::from_owned(req.design.3.clone()),
        n_design_cols: req.design.1,
        design_complete: Owned::from_owned(req.design_complete.clone()),
        group_index: Owned::from_owned(req.group_index.clone()),
        n_group_levels: req.group_labels.len(),
        fix_eff: Owned::from_owned(
            (0..req.fix_eff.len())
                .map(|i| req.fix_eff.get(i).unwrap_or_default().into_owned())
                .collect(),
        ),
        group_labels: Owned::from_owned(
            (0..req.group_labels.len())
                .map(|i| req.group_labels.get(i).unwrap_or_default().into_owned())
                .collect(),
        ),
        taxon_names: Owned::from_owned(
            (0..req.taxon_names.len())
                .map(|i| req.taxon_names.get(i).unwrap_or_default().into_owned())
                .collect(),
        ),
        sample_names: Owned::from_owned(
            (0..req.sample_names.len())
                .map(|i| req.sample_names.get(i).unwrap_or_default().into_owned())
                .collect(),
        ),
        controls: req.controls,
        flags: req.flags,
    }
}

// ---------------------------------------------------------------------------
// The probe: a transport with no numerics behind it
// ---------------------------------------------------------------------------

/// Which buffer to echo back. One probe per buffer keeps the shim's arithmetic
/// simple and makes a mismatch name the buffer that broke.
#[repr(C)]
#[derive(Clone, Copy)]
pub enum Echo {
    /// `data` when it is `INTSXP`: `NA_integer_` must come back as `INT_MIN`.
    DataInt,
    /// `data` when it is `REALSXP`: `NA_real_`, `NaN`, `+-Inf` and `-0.0` must
    /// each come back as their own 64 bits.
    DataReal,
    /// `aggregate_data`, in whichever of R's two numeric types it arrived. The
    /// type is echoed rather than widened, because the returned `feature_table`
    /// is this table and its integer-ness is part of the schema.
    Aggregate,
    /// The design, column-major.
    Design,
    /// `stats::complete.cases(x)`.
    DesignComplete,
    /// The 1-based group index per sample, `0` meaning absent.
    GroupIndex,
}

/// What kind of R vector the caller should build for a probe's payload.
pub const KIND_NONE: i32 = 0;
/// An `INTSXP` payload.
pub const KIND_INT: i32 = 1;
/// A `REALSXP` payload.
pub const KIND_REAL: i32 = 2;
/// A `RAWSXP` payload.
pub const KIND_RAW: i32 = 3;

/// **Transport only.** Copy every buffer out of R's memory and back, bit for bit,
/// and write the derived metadata alongside.
///
/// This computes nothing. That is the point: it is the evidence that the transport
/// preserves what it is supposed to preserve, so that when a *numerical* stage
/// later disagrees with the original, the disagreement cannot be blamed on the
/// boundary.
///
/// Returns [`RB_OK`] or one of the error codes; the message is in the last-error
/// buffer.
///
/// # Safety
///
/// Every pointer in `req` must be valid for the length that accompanies it, and
/// every `out_*` pointer must be either null or valid for its stated length. The
/// caller is `init.c`, which allocated them; see `r/ANCOMBC/src/init.c`.
#[no_mangle]
pub unsafe extern "C" fn ancombc2_rb_probe(
    req: *const RawRequest,
    which: i32,
    out_int: *mut i32,
    out_real: *mut f64,
    out_u8: *mut u8,
    out_i64: *mut i64,
    out_n: i64,
    out_kind: *mut i32,
) -> i32 {
    let req = match req.as_ref() {
        Some(r) => r,
        None => {
            set_error("the request pointer was null");
            return RB_ERR_ARGS;
        }
    };
    match catch_unwind(AssertUnwindSafe(|| {
        probe_inner(
            req, which, out_int, out_real, out_u8, out_i64, out_n, out_kind,
        )
    })) {
        Ok(Ok(())) => RB_OK,
        Ok(Err((code, msg))) => {
            set_error(msg);
            code
        }
        Err(payload) => {
            set_error(format!(
                "ancombc2-rbridge panicked: {}",
                panic_message(&payload)
            ));
            RB_ERR_PANIC
        }
    }
}

fn panic_message(p: &Box<dyn std::any::Any + Send>) -> String {
    if let Some(s) = p.downcast_ref::<&str>() {
        (*s).to_string()
    } else if let Some(s) = p.downcast_ref::<String>() {
        s.clone()
    } else {
        "<non-string panic payload>".to_string()
    }
}

type Coded = Result<(), (i32, String)>;

#[allow(clippy::too_many_arguments)]
unsafe fn probe_inner(
    req: &RawRequest,
    which: i32,
    out_int: *mut i32,
    out_real: *mut f64,
    out_u8: *mut u8,
    out_i64: *mut i64,
    out_n: i64,
    out_kind: *mut i32,
) -> Coded {
    let args_err = |e: BridgeError| (RB_ERR_ARGS, e.to_string());
    let parsed = to_request(req).map_err(args_err)?;

    // The owned buffer a parallel stage would get. Built and dropped here, which
    // is the whole lifetime story: no R pointer survives this call.
    // Built and dropped inside this call. That is the entire lifetime story: no
    // R pointer survives, and nothing here is handed to Rayon.
    let staged = stage_input(&parsed);
    drop(staged);

    let n = out_n.max(0) as usize;
    // Assigned on every arm that reaches the end of the match; `KIND_NONE` is the
    // value an unknown probe leaves behind, and that arm returns an error rather
    // than writing a kind.
    #[allow(unused_assignments)]
    let mut kind = KIND_NONE;
    match which {
        x if x == Echo::DataInt as i32 => {
            let Counts::Int { values, .. } = &parsed.data else {
                return Err((
                    RB_ERR_ARGS,
                    "this case carries REALSXP counts; the integer probe does not apply".into(),
                ));
            };
            if out_int.is_null() {
                return Err((RB_ERR_ARGS, "out_int was null".into()));
            }
            let want = values.len().min(n);
            std::ptr::copy_nonoverlapping(values.as_ptr(), out_int, want);
            write_i64(out_i64, values.len() as i64);
            kind = KIND_INT;
        }
        x if x == Echo::DataReal as i32 => {
            let Counts::Real { values, .. } = &parsed.data else {
                return Err((
                    RB_ERR_ARGS,
                    "this case carries INTSXP counts; the real probe does not apply".into(),
                ));
            };
            if out_real.is_null() {
                return Err((RB_ERR_ARGS, "out_real was null".into()));
            }
            let want = values.len().min(n);
            std::ptr::copy_nonoverlapping(values.as_ptr(), out_real, want);
            write_i64(out_i64, values.len() as i64);
            kind = KIND_REAL;
        }
        x if x == Echo::Aggregate as i32 => {
            // Echoed in its own type. Widening here would be invisible in the
            // numbers and visible in the schema, which is the worse of the two.
            match &parsed.aggregate {
                Counts::Int { values, .. } => {
                    if out_int.is_null() {
                        return Err((RB_ERR_ARGS, "out_int was null".into()));
                    }
                    let want = values.len().min(n);
                    std::ptr::copy_nonoverlapping(values.as_ptr(), out_int, want);
                    write_i64(out_i64, values.len() as i64);
                    kind = KIND_INT;
                }
                Counts::Real { values, .. } => {
                    if out_real.is_null() {
                        return Err((RB_ERR_ARGS, "out_real was null".into()));
                    }
                    let want = values.len().min(n);
                    std::ptr::copy_nonoverlapping(values.as_ptr(), out_real, want);
                    write_i64(out_i64, values.len() as i64);
                    kind = KIND_REAL;
                }
            }
        }
        x if x == Echo::Design as i32 => {
            if out_real.is_null() {
                return Err((RB_ERR_ARGS, "out_real was null".into()));
            }
            let want = parsed.design.3.len().min(n);
            std::ptr::copy_nonoverlapping(parsed.design.3.as_ptr(), out_real, want);
            write_i64(out_i64, parsed.design.3.len() as i64);
            kind = KIND_REAL;
        }
        x if x == Echo::DesignComplete as i32 => {
            if out_u8.is_null() {
                return Err((RB_ERR_ARGS, "out_u8 was null".into()));
            }
            let v = &parsed.design_complete;
            let want = v.len().min(n);
            std::ptr::copy_nonoverlapping(v.as_ptr(), out_u8, want);
            write_i64(out_i64, v.len() as i64);
            kind = KIND_RAW;
        }
        x if x == Echo::GroupIndex as i32 => {
            if out_int.is_null() {
                return Err((RB_ERR_ARGS, "out_int was null".into()));
            }
            let v = &parsed.group_index;
            let want = v.len().min(n);
            std::ptr::copy_nonoverlapping(v.as_ptr(), out_int, want);
            write_i64(out_i64, v.len() as i64);
            kind = KIND_INT;
        }
        other => return Err((RB_ERR_STAGE_UNKNOWN, format!("unknown probe {other}"))),
    }
    write_i32(out_kind, kind);
    Ok(())
}

unsafe fn write_i64(out: *mut i64, v: i64) {
    if !out.is_null() {
        *out = v;
    }
}

unsafe fn write_i32(out: *mut i32, v: i32) {
    if !out.is_null() {
        *out = v;
    }
}

// ---------------------------------------------------------------------------
// The output transport
// ---------------------------------------------------------------------------

/// The core payloads as they arrive from R, for [`ancombc2_rb_emit`].
///
/// Every buffer is borrowed and validated; see [`ancombc2_rb_emit`] for the
/// contract. This struct exists so the argument list of that function is one
/// pointer rather than thirteen.
#[repr(C)]
pub struct RawPayloads {
    /// `beta_hat`, `n_tax * p` doubles.
    pub beta_hat: *const f64,
    /// `var_hat`, `n_tax * p` doubles.
    pub var_hat: *const f64,
    /// `dof`, `n_tax * p`. Doubles, because an `INTSXP` would have been widened
    /// on the R side and the *original's* type is what matters.
    pub dof: *const f64,
    /// Non-zero when the original's `dof` was integer, which is what decides
    /// whether it comes back as `INTSXP`.
    pub dof_is_int: i32,
    /// The covariances, `n_tax * p * p` doubles.
    pub vcov: *const f64,
    /// `y_bias_crt`, `n_tax * n_samp` doubles.
    pub y_bias_crt: *const f64,
    /// `theta_hat`, `n_samp` doubles.
    pub theta_hat: *const f64,
    /// `delta_em`, `delta_wls`, `var_delta`, `s02`: `p` doubles each.
    pub delta_em: *const f64,
    /// See [`Self::delta_em`].
    pub delta_wls: *const f64,
    /// See [`Self::delta_em`].
    pub var_delta: *const f64,
    /// See [`Self::delta_em`].
    pub s02: *const f64,
    /// `O2`, `n_tax * n_samp` doubles.
    pub o2: *const f64,
    /// `taxa`, `nrow(O2)` indices, 1-based as R's `which()` returns them.
    pub taxa: *const i64,
    /// `taxa_bias`, indices into the input table for the taxa in `O1`. May be
    /// longer than `taxa`: `n_tax` is `nrow(O2)`, and `struc_zero` drops rows
    /// from `O2` that `O1` keeps.
    pub taxa_bias: *const i64,
    /// Packed `rownames(O2)`.
    pub taxon_names_offsets: *const i64,
    /// See [`Self::taxon_names_offsets`].
    pub taxon_names_bytes: *const u8,
    /// `nrow(O2)`.
    pub n_taxon_names: i64,
    /// Packed `colnames(O2)`.
    pub sample_names_offsets: *const i64,
    /// See [`Self::sample_names_offsets`].
    pub sample_names_bytes: *const u8,
    /// `ncol(O2)`.
    pub n_sample_names: i64,
    /// Packed `colnames(x)`.
    pub fix_eff_offsets: *const i64,
    /// See [`Self::fix_eff_offsets`].
    pub fix_eff_bytes: *const u8,
    /// `ncol(x)`.
    pub n_fix_eff: i64,
    /// `nrow(O2)`.
    pub n_tax: i64,
    /// `ncol(O2)`.
    pub n_samp: i64,
    /// `ncol(x)`.
    pub p: i64,
    /// `length(taxa_bias)`. Its own count, because it is not `n_tax`.
    pub n_taxa_bias: i64,
}

unsafe fn raw_payloads(r: &RawPayloads) -> Result<CorePayloads, BridgeError> {
    let n_tax = as_usize(r.n_tax, "n_tax")?;
    let n_samp = as_usize(r.n_samp, "n_samp")?;
    let p = as_usize(r.p, "p")?;
    // Every product below is sized into a `Vec`, so an unchecked count turns into
    // a multi-terabyte allocation instead of an error. The real bound is `p <=
    // ncol(x)`, which R enforces; this is the belt to that braces, and it also
    // catches a caller whose struct layout disagrees with this one -- which is
    // exactly how a real bug got here, reading a pointer's low bits as a count.
    const MAX_TAX: usize = 1 << 20;
    const MAX_SAMP: usize = 1 << 24;
    const MAX_P: usize = 1 << 16;
    for (field, got, limit) in [
        ("n_tax", n_tax, MAX_TAX),
        ("n_samp", n_samp, MAX_SAMP),
        ("p", p, MAX_P),
    ] {
        if got > limit {
            return Err(BridgeError::TooLarge {
                field,
                got: got as i64,
                limit: limit as i64,
            });
        }
    }
    let nx = n_tax * p;
    let ns = n_tax * n_samp;
    let nv = n_tax * p * p;
    let nn = as_usize(r.n_taxon_names, "taxon_names")?;
    let nsn = as_usize(r.n_sample_names, "sample_names")?;
    let nf = as_usize(r.n_fix_eff, "fix_eff")?;
    Ok(CorePayloads {
        n_tax,
        n_samp,
        p,
        beta_hat: slice::<f64>(r.beta_hat, nx, "beta_hat")?.to_vec(),
        var_hat: slice::<f64>(r.var_hat, nx, "var_hat")?.to_vec(),
        dof: slice::<f64>(r.dof, nx, "dof")?.to_vec(),
        dof_is_int: r.dof_is_int != 0,
        vcov: slice::<f64>(r.vcov, nv, "vcov_hat")?.to_vec(),
        y_bias_crt: slice::<f64>(r.y_bias_crt, ns, "y_bias_crt")?.to_vec(),
        theta_hat: slice::<f64>(r.theta_hat, n_samp, "theta_hat")?.to_vec(),
        delta_em: slice::<f64>(r.delta_em, p, "delta_em")?.to_vec(),
        delta_wls: slice::<f64>(r.delta_wls, p, "delta_wls")?.to_vec(),
        var_delta: slice::<f64>(r.var_delta, p, "var_delta")?.to_vec(),
        s02: slice::<f64>(r.s02, p, "s02")?.to_vec(),
        o2: slice::<f64>(r.o2, ns, "O2")?.to_vec(),
        taxa: slice::<i64>(r.taxa, n_tax, "taxa")?
            .iter()
            .map(|&v| v as usize)
            .collect(),
        taxa_bias: slice::<i64>(
            r.taxa_bias,
            as_usize(r.n_taxa_bias, "taxa_bias")?,
            "taxa_bias",
        )?
        .iter()
        .map(|&v| v as usize)
        .collect(),
        taxon_names: packed(
            r.taxon_names_offsets,
            r.taxon_names_bytes,
            nn,
            "taxon_names",
        )?,
        sample_names: packed(
            r.sample_names_offsets,
            r.sample_names_bytes,
            nsn,
            "sample_names",
        )?,
        fix_eff: packed(r.fix_eff_offsets, r.fix_eff_bytes, nf, "fix_eff")?,
    })
}

/// How many `f64`, `i32` and `u8` the plan needs, without writing anything.
///
/// The caller allocates from these three numbers, so a plan and an allocation
/// cannot disagree about a length.
///
/// # Safety
///
/// `pay` must describe a valid payload set.
#[no_mangle]
pub unsafe extern "C" fn ancombc2_rb_emit_plan(
    pay: *const RawPayloads,
    out_f64_len: *mut i64,
    out_i32_len: *mut i64,
    out_u8_len: *mut i64,
    out_n: *mut i64,
) -> i32 {
    let r = match pay.as_ref() {
        Some(r) => r,
        None => {
            set_error("the payload pointer was null");
            return RB_ERR_ARGS;
        }
    };
    match catch_unwind(AssertUnwindSafe(|| {
        let p = to_result(raw_payloads(r))?;
        p.validate().map_err(|e| (RB_ERR_ARGS, e.to_string()))?;
        let plan = core_output_plan(p.n_tax, p.n_samp, p.p, p.dof_is_int)
            .map_err(|e| (RB_ERR_ARGS, e.to_string()))?;
        let mut e = Emitted::default();
        emit_core_payloads(&p, &plan, &mut e).map_err(|e| (RB_ERR_ARGS, e.to_string()))?;
        unsafe {
            if !out_f64_len.is_null() {
                *out_f64_len = e.f64_len as i64;
            }
            if !out_i32_len.is_null() {
                *out_i32_len = e.i32_len as i64;
            }
            if !out_u8_len.is_null() {
                *out_u8_len = e.u8_len as i64;
            }
            if !out_n.is_null() {
                *out_n = plan.len() as i64;
            }
        }
        Ok(())
    })) {
        Ok(Ok(())) => RB_OK,
        Ok(Err((code, msg))) => {
            set_error(msg);
            code
        }
        Err(payload) => {
            set_error(format!(
                "ancombc2-rbridge panicked in emit_plan: {}",
                panic_message(&payload)
            ));
            RB_ERR_PANIC
        }
    }
}

/// **Output transport.** Copy the core's payloads out of R's memory and back into
/// R's memory, bit for bit, in the types and shapes the plan declares.
///
/// Nothing is computed. This is the same discipline as
/// [`ancombc2_rb_probe`]: before any stage's output can be trusted, the path it
/// travels has to be shown to lose nothing. `NA_real_`, `NaN`, `±Inf` and `-0.0`
/// all cross as raw bits, `dof` comes back as `INTSXP` when the original's was
/// integer, and `vcov_hat` comes back as a list of `p x p` matrices with their
/// `dimnames` intact.
///
/// # Safety
///
/// Every pointer in `pay` must be valid for the length the accompanying counts
/// imply, and the `out_*` pointers must be null or valid for the lengths
/// `ancombc2_rb_emit_plan` reported. The caller is `init.c`, which allocated them.
#[no_mangle]
pub unsafe extern "C" fn ancombc2_rb_emit(
    pay: *const RawPayloads,
    out_f64: *mut f64,
    out_i32: *mut i32,
    out_offsets_f64: *mut i64,
    out_offsets_i32: *mut i64,
    out_plan_names_offsets: *mut i64,
    out_plan_names_bytes: *mut u8,
    out_plan_len: i64,
) -> i32 {
    let r = match pay.as_ref() {
        Some(r) => r,
        None => {
            set_error("the payload pointer was null");
            return RB_ERR_ARGS;
        }
    };
    match catch_unwind(AssertUnwindSafe(|| {
        let p = to_result(raw_payloads(r))?;
        p.validate().map_err(|e| (RB_ERR_ARGS, e.to_string()))?;
        let plan = core_output_plan(p.n_tax, p.n_samp, p.p, p.dof_is_int)
            .map_err(|e| (RB_ERR_ARGS, e.to_string()))?;
        let mut e = Emitted::default();
        emit_core_payloads(&p, &plan, &mut e).map_err(|e| (RB_ERR_ARGS, e.to_string()))?;
        if out_plan_len < plan.len() as i64 {
            return Err((
                RB_ERR_ARGS,
                format!(
                    "the plan has {} entries but the caller reserved {}",
                    plan.len(),
                    out_plan_len
                ),
            ));
        }

        // One owned copy, then write it out. This is where an R pointer would
        // stop being allowed to matter: `p` owns everything below this line.
        let owned = p.to_owned();

        let src_for = |name: &str| -> &[f64] {
            match name {
                "theta_hat" => &owned.theta_hat,
                "beta_hat" => &owned.beta_hat,
                "var_hat" => &owned.var_hat,
                "dof" => &owned.dof,
                "vcov_hat" => &owned.vcov,
                "y_bias_crt" => &owned.y_bias_crt,
                "O2" => &owned.o2,
                "delta_em" => &owned.delta_em,
                "delta_wls" => &owned.delta_wls,
                "var_delta" => &owned.var_delta,
                "s02" => &owned.s02,
                _ => &[],
            }
        };

        let mut w64 = 0i64;
        let mut w32 = 0i64;
        for (i, d) in plan.iter().enumerate() {
            let src = src_for(d.name);
            unsafe {
                if !out_offsets_f64.is_null() {
                    *out_offsets_f64.add(i) = e.f64_offsets[i];
                }
                if !out_offsets_i32.is_null() {
                    *out_offsets_i32.add(i) = e.i32_offsets[i];
                }
            }
            // The offset comes from the list that matches the type. Reading the
            // `f64` offset for an `INTSXP` entry gave -1, which looked like "this
            // entry writes nothing" and silently dropped `dof` -- the plan's
            // declared length then disagreed with the written one, which is how it
            // was caught.
            match d.rtype {
                crate::output::RType::Real | crate::output::RType::ListOfRealMatrix => {
                    let off = e.f64_offsets[i];
                    if off < 0 || out_f64.is_null() {
                        return Err((RB_ERR_ARGS, format!("no f64 slot for {}", d.name)));
                    }
                    let dst = unsafe { out_f64.add(off as usize) };
                    unsafe { std::ptr::copy_nonoverlapping(src.as_ptr(), dst, src.len()) };
                    w64 += d.count as i64;
                }
                crate::output::RType::Int => {
                    let off = e.i32_offsets[i];
                    if off < 0 || out_i32.is_null() {
                        return Err((RB_ERR_ARGS, format!("no i32 slot for {}", d.name)));
                    }
                    let dst32 = unsafe { out_i32.add(off as usize) };
                    for (k, &v) in src.iter().enumerate() {
                        // `dof` arrives as doubles because that is the only type
                        // the wire carries; it goes back as `int` because the
                        // original's was integer. The conversion is exact only
                        // when the values are whole, which `dof_is_int` asserts,
                        // so it is checked rather than trusted.
                        let r = v.round();
                        if (v - r).abs() > 0.0 {
                            return Err((
                                RB_ERR_ARGS,
                                format!(
                                    "dof[{}] = {v} is not whole but the caller said the \\
                                     original's dof was integer",
                                    k
                                ),
                            ));
                        }
                        unsafe { *dst32.add(k) = r as i32 };
                    }
                    w32 += d.count as i64;
                }
                _ => {
                    return Err((RB_ERR_ARGS, format!("unhandled type for {}", d.name)));
                }
            }
        }
        if w64 != e.f64_len as i64 || w32 != e.i32_len as i64 {
            return Err((
                RB_ERR_ARGS,
                "the plan's declared lengths and the written lengths disagree".into(),
            ));
        }

        // The plan's names, packed, so R can check what it got against what it
        // asked for without a second source of truth.
        unsafe {
            let mut offs: Vec<i64> = Vec::with_capacity(plan.len() + 1);
            let mut bytes: Vec<u8> = Vec::new();
            offs.push(0);
            for d in &plan {
                bytes.extend_from_slice(d.name.as_bytes());
                offs.push(bytes.len() as i64);
            }
            if !out_plan_names_offsets.is_null() {
                std::ptr::copy_nonoverlapping(offs.as_ptr(), out_plan_names_offsets, offs.len());
            }
            if !out_plan_names_bytes.is_null() {
                std::ptr::copy_nonoverlapping(bytes.as_ptr(), out_plan_names_bytes, bytes.len());
            }
        }
        Ok(())
    })) {
        Ok(Ok(())) => RB_OK,
        Ok(Err((code, msg))) => {
            set_error(msg);
            code
        }
        Err(payload) => {
            set_error(format!(
                "ancombc2-rbridge panicked in emit: {}",
                panic_message(&payload)
            ));
            RB_ERR_PANIC
        }
    }
}

#[allow(non_snake_case)]
fn to_result<T>(r: Result<T, BridgeError>) -> Result<T, (i32, String)> {
    r.map_err(|e| (RB_ERR_ARGS, e.to_string()))
}

#[cfg(test)]
mod output_layout_tests {
    //! The wire struct is shared with `r/ANCOMBC/src/init.c`, which declares a C
    //! struct of the same shape. Nothing in the type system ties the two together:
    //! a field moved in one is silently read as a different field in the other.
    //!
    //! The failure is not subtle. A pointer's low bits read as a count produced an
    //! 812 GB `Vec<i64>` and aborted the process -- loud, but only after the fact.
    //! These offsets make it fail at build time instead, and they double as the
    //! specification the C side is written against.

    use super::RawPayloads;
    use core::mem::{align_of, offset_of, size_of};

    /// The layout `r/ANCOMBC/src/init.c`'s `struct rb_payloads` must have.
    ///
    /// On 64-bit: a pointer is 8 bytes, `i32` 4, `i64` 8. `dof_is_int` is the only
    /// 4-byte field, so everything after it is 8-byte aligned and the offsets are
    /// a plain running sum -- which is what makes this list short enough to read.
    #[test]
    fn the_wire_layout_is_pinned() {
        assert_eq!(size_of::<RawPayloads>(), 216);
        assert_eq!(align_of::<RawPayloads>(), 8);
        assert_eq!(offset_of!(RawPayloads, beta_hat), 0);
        assert_eq!(offset_of!(RawPayloads, var_hat), 8);
        assert_eq!(offset_of!(RawPayloads, dof), 16);
        assert_eq!(offset_of!(RawPayloads, dof_is_int), 24);
        assert_eq!(offset_of!(RawPayloads, vcov), 32);
        assert_eq!(offset_of!(RawPayloads, y_bias_crt), 40);
        assert_eq!(offset_of!(RawPayloads, theta_hat), 48);
        assert_eq!(offset_of!(RawPayloads, delta_em), 56);
        assert_eq!(offset_of!(RawPayloads, delta_wls), 64);
        assert_eq!(offset_of!(RawPayloads, var_delta), 72);
        assert_eq!(offset_of!(RawPayloads, s02), 80);
        assert_eq!(offset_of!(RawPayloads, o2), 88);
        assert_eq!(offset_of!(RawPayloads, taxa), 96);
        assert_eq!(offset_of!(RawPayloads, taxa_bias), 104);
        assert_eq!(offset_of!(RawPayloads, taxon_names_offsets), 112);
        assert_eq!(offset_of!(RawPayloads, taxon_names_bytes), 120);
        assert_eq!(offset_of!(RawPayloads, n_taxon_names), 128);
        assert_eq!(offset_of!(RawPayloads, sample_names_offsets), 136);
        assert_eq!(offset_of!(RawPayloads, sample_names_bytes), 144);
        assert_eq!(offset_of!(RawPayloads, n_sample_names), 152);
        assert_eq!(offset_of!(RawPayloads, fix_eff_offsets), 160);
        assert_eq!(offset_of!(RawPayloads, fix_eff_bytes), 168);
        assert_eq!(offset_of!(RawPayloads, n_fix_eff), 176);
        assert_eq!(offset_of!(RawPayloads, n_tax), 184);
        assert_eq!(offset_of!(RawPayloads, n_samp), 192);
        assert_eq!(offset_of!(RawPayloads, p), 200);
        // Last, and deliberately: `n_taxa_bias` belongs to no shape above, so it is
        // the one field a reader can move by accident without noticing locally.
        assert_eq!(offset_of!(RawPayloads, n_taxa_bias), 208);
    }

    /// A count that is impossible rather than merely wrong gets a named error
    /// instead of an allocation.
    ///
    /// This is the case that produced the 812 GB abort: the count came from a
    /// misaligned field, so it was enormous, self-consistent, and impossible.
    #[test]
    fn an_impossible_count_is_refused_before_anything_is_sized_from_it() {
        let mut pay = RawPayloads {
            beta_hat: std::ptr::null(),
            var_hat: std::ptr::null(),
            dof: std::ptr::null(),
            dof_is_int: 0,
            vcov: std::ptr::null(),
            y_bias_crt: std::ptr::null(),
            theta_hat: std::ptr::null(),
            delta_em: std::ptr::null(),
            delta_wls: std::ptr::null(),
            var_delta: std::ptr::null(),
            s02: std::ptr::null(),
            o2: std::ptr::null(),
            taxa: std::ptr::null(),
            taxa_bias: std::ptr::null(),
            taxon_names_offsets: std::ptr::null(),
            taxon_names_bytes: std::ptr::null(),
            n_taxon_names: 0,
            sample_names_offsets: std::ptr::null(),
            sample_names_bytes: std::ptr::null(),
            n_sample_names: 0,
            fix_eff_offsets: std::ptr::null(),
            fix_eff_bytes: std::ptr::null(),
            n_fix_eff: 0,
            n_tax: 0,
            n_samp: 0,
            p: 0,
            n_taxa_bias: 0,
        };
        pay.n_tax = 1 << 30;
        let err = unsafe { super::raw_payloads(&pay) }.unwrap_err();
        assert!(
            matches!(err, crate::BridgeError::TooLarge { field: "n_tax", .. }),
            "got {err}"
        );
    }
}

#[cfg(test)]
mod preprocess_probe_tests {
    use crate::output::{preprocess_stages, PreprocessInputs};
    use ancombc2_core::preprocess::CountMatrix;
    use ancombc2_core::reduce::F64Reductions;

    fn counts(n_tax: usize, n_samp: usize, data: &[f64]) -> CountMatrix {
        CountMatrix::new(n_tax, n_samp, data.to_vec()).unwrap()
    }

    /// `log(0)` must become `NA_real_` and not a bare `f64::NAN`.
    ///
    /// `identical(NA_real_, NaN)` is `FALSE` in R and `serialize()` records the
    /// payload, so this is a byte difference and not a cosmetic one. It is the
    /// first stage where a lost or invented NaN payload shows up.
    #[test]
    fn the_log_zero_stage_holds_na_real_and_not_a_bare_nan() {
        let data = counts(2, 3, &[0.0, 1.0, 4.0, 9.0, 2.0, 0.0]);
        let agg = data.clone();
        let st = preprocess_stages(&PreprocessInputs {
            data: &data,
            aggregate: &agg,
            struc_zero: false,
            group_indicator: None,
            n_groups: 0,
            prv_cut: 0.0,
            lib_cut: 0.0,
            pseudo: 0.0,
            neg_lb: false,
            red: &F64Reductions,
        })
        .unwrap();
        let z = st.log1.get(0, 0);
        assert!(z.is_nan());
        assert_eq!(z.to_bits(), ancombc2_core::reduce::NA_REAL_BITS);
        assert_ne!(z.to_bits(), f64::NAN.to_bits());
    }

    /// An `NA_real_` in the counts keeps its payload through `log` and through the
    /// centring -- and is *quieted* by the arithmetic, because R quiets a signalling
    /// NaN rather than passing it through. R's `is.infinite` mask does not touch a NaN,
    /// so nothing in the stage replaces it either.
    #[test]
    fn an_na_count_survives_the_log_stage_unchanged() {
        let na = ancombc2_core::reduce::na_real();
        let data = counts(1, 3, &[na, 4.0, 9.0]);
        let agg = data.clone();
        let st = preprocess_stages(&PreprocessInputs {
            data: &data,
            aggregate: &agg,
            struc_zero: false,
            group_indicator: None,
            n_groups: 0,
            prv_cut: 0.0,
            lib_cut: 0.0,
            pseudo: 0.0,
            neg_lb: false,
            red: &F64Reductions,
        })
        .unwrap();
        // Quieted, not `NA_real_`, and the reason is the stage before: `O1 = data +
        // pseudo` is arithmetic on the `NA`, so `NA_real_ + 0` is already
        // `7ff80000000007a2` before `log` ever sees it. R agrees -- `log(NA_real_ + 0)`
        // is the quieted NA, because R's unary `log` passes a NaN through unchanged and
        // the quieting already happened one step earlier.
        let v = st.log1.get(0, 0);
        assert_eq!(v.to_bits(), ancombc2_core::reduce::QUIETED_NA_REAL_BITS);
        // And the centred row is NA there too -- with the *quieted* payload, because
        // R's arithmetic quiets a signalling NaN rather than passing it through.
        assert_eq!(
            st.y1.get(0, 0).to_bits(),
            ancombc2_core::reduce::QUIETED_NA_REAL_BITS
        );
    }

    /// `prv_cut` removes taxa, and the library size is then computed over what
    /// *remains* -- not over the whole table.
    ///
    /// The reference subsets before summing, so a taxon removed by prevalence must
    /// not contribute to `lib_size`. Summing first and filtering after would give a
    /// different `lib_size` and therefore a different `samp_keep` whenever
    /// `lib_cut` is anywhere near a column sum.
    #[test]
    fn the_library_size_is_over_the_retained_taxa_only() {
        // Taxon 0 is all zeros, so prevalence 0 and it goes at prv_cut = 0.5.
        // Taxon 1 is 1,2,3.
        let data = counts(2, 3, &[0.0, 0.0, 0.0, 1.0, 2.0, 3.0]);
        let agg = data.clone();
        let st = preprocess_stages(&PreprocessInputs {
            data: &data,
            aggregate: &agg,
            struc_zero: false,
            group_indicator: None,
            n_groups: 0,
            prv_cut: 0.5,
            lib_cut: 0.0,
            pseudo: 0.0,
            neg_lb: false,
            red: &F64Reductions,
        })
        .unwrap();
        assert_eq!(st.tax_keep1, vec![2], "only taxon 2 survives");
        assert_eq!(st.lib_size1, vec![1.0, 2.0, 3.0]);
        // And the retained table really is one taxon by three samples.
        assert_eq!(st.o1.rows, 1);
        assert_eq!(st.o1.cols, 3);
        // This input does *not* distinguish the two rules, because the removed taxon
        // contributed nothing. The next test uses one that does.
        assert_eq!(data.library_sizes(), vec![1.0, 2.0, 3.0]);
    }

    /// The same check with a removed taxon that carries a large count, so the two
    /// rules genuinely disagree and the test cannot pass by accident.
    #[test]
    fn a_removed_taxa_would_change_the_library_size_if_it_were_summed() {
        // Taxon 0 is rare -- prevalence 1/3, below prv_cut = 0.5 -- but carries the
        // 1000. Taxon 1 is present in every sample. A library size summed before the
        // taxon screen would report 1005 for the first sample; this reports 5.
        let data = counts(2, 3, &[1000.0, 0.0, 0.0, 5.0, 7.0, 8.0]);
        let agg = data.clone();
        let st = preprocess_stages(&PreprocessInputs {
            data: &data,
            aggregate: &agg,
            struc_zero: false,
            group_indicator: None,
            n_groups: 0,
            prv_cut: 0.5,
            lib_cut: 0.0,
            pseudo: 0.0,
            neg_lb: false,
            red: &F64Reductions,
        })
        .unwrap();
        assert_eq!(st.tax_keep1, vec![2]);
        assert_eq!(
            st.lib_size1,
            vec![5.0, 7.0, 8.0],
            "the 1000 must not be in the sum"
        );
        // The wrong answer, recorded so the assertion above has something to be
        // different from.
        assert_eq!(
            data.library_sizes(),
            vec![1005.0, 7.0, 8.0],
            "so the two rules genuinely differ on this input"
        );
        assert_eq!(st.o1.rows, 1);
        assert_eq!(st.o1.cols, 3);
    }

    /// `lib_cut` is a comparison against the data, not a switch, and the screen runs
    /// on the *retained* taxa.
    ///
    /// A sample whose retained library size falls below `lib_cut` is dropped. Here
    /// the first sample is below it, and the third is exactly on it -- `>=`, not `>`,
    /// so the third survives. The boundary is the whole point: `>` would drop it and
    /// change `samp_keep` by one.
    #[test]
    fn lib_cut_is_a_comparison_and_the_boundary_is_inclusive() {
        // Taxon 0: 1, 1, 10. One taxon, so it always survives the taxon screen.
        let data = counts(1, 3, &[1.0, 1.0, 10.0]);
        let agg = data.clone();
        let st = preprocess_stages(&PreprocessInputs {
            data: &data,
            aggregate: &agg,
            struc_zero: false,
            group_indicator: None,
            n_groups: 0,
            prv_cut: 0.0,
            lib_cut: 10.0,
            pseudo: 0.0,
            neg_lb: false,
            red: &F64Reductions,
        })
        .unwrap();
        assert_eq!(st.lib_size1, vec![1.0, 1.0, 10.0]);
        assert_eq!(
            st.samp_keep1,
            vec![3],
            "samples 1 and 2 are below 10; sample 3 is exactly on it and stays"
        );
        assert_eq!(st.o1.cols, 1, "one sample retained");
    }

    /// The row mean is R's `rowMeans`, so it is the one stage where the f64
    /// fallback is visibly wrong. This test pins the *fallback's* behaviour and
    /// says so; `scripts/check_preprocess_stages.R` pins R's.
    #[test]
    fn the_reduction_is_the_only_place_f64_and_r_disagree() {
        // log of a cancelling pair: 1e16 and 1 differ in log by ~36, so the row sum
        // does not cancel and both accumulators agree here. That is the point: the
        // f64 path is fine for ordinary counts and wrong only for adversarial ones.
        let data = counts(1, 2, &[1.0, 4.0]);
        let agg = data.clone();
        let st = preprocess_stages(&PreprocessInputs {
            data: &data,
            aggregate: &agg,
            struc_zero: false,
            group_indicator: None,
            n_groups: 0,
            prv_cut: 0.0,
            lib_cut: 0.0,
            pseudo: 0.0,
            neg_lb: false,
            red: &F64Reductions,
        })
        .unwrap();
        let expected = (1.0f64.ln() + 4.0f64.ln()) / 2.0;
        assert_eq!(st.means1[0], expected);
    }

    /// `pseudo` is added before the log, and the addition is part of the stage.
    #[test]
    fn pseudo_is_added_before_the_log() {
        let data = counts(1, 2, &[1.0, 3.0]);
        let agg = data.clone();
        let st = preprocess_stages(&PreprocessInputs {
            data: &data,
            aggregate: &agg,
            struc_zero: false,
            group_indicator: None,
            n_groups: 0,
            prv_cut: 0.0,
            lib_cut: 0.0,
            pseudo: 0.5,
            neg_lb: false,
            red: &F64Reductions,
        })
        .unwrap();
        assert_eq!(st.o1.row(0), &[1.5, 3.5]);
        assert_eq!(st.log1.get(0, 0), 1.5f64.ln());
        assert_eq!(st.log1.get(0, 1), 3.5f64.ln());
    }

    /// An all-missing taxon cannot reach the log stage, and that is R's rule.
    ///
    /// `which(prevalence >= prv_cut)` drops it: prevalence is `0/0`, which in R is
    /// `NaN`, and `which` drops `NaN`. So `rowMeans` never sees an all-`NA` row, the
    /// `NA_real_` mean branch in `F64Reductions` is defensive rather than
    /// reachable, and the `y1` matrix always has at least one observed entry per row.
    ///
    /// This was worth checking rather than assuming: a taxon that is missing in
    /// every sample is exactly the case a hand-written filter would plausibly let
    /// through, and it would then divide by zero here.
    #[test]
    fn an_all_missing_taxa_never_reach_the_log_stage() {
        let na = ancombc2_core::reduce::na_real();
        let data = counts(2, 2, &[na, na, 1.0, 3.0]);
        let agg = data.clone();
        let st = preprocess_stages(&PreprocessInputs {
            data: &data,
            aggregate: &agg,
            struc_zero: false,
            group_indicator: None,
            n_groups: 0,
            prv_cut: 0.0,
            lib_cut: 0.0,
            pseudo: 0.0,
            neg_lb: false,
            red: &F64Reductions,
        })
        .unwrap();
        assert_eq!(st.tax_keep1, vec![2], "the all-missing taxon is dropped");
        assert_eq!(st.o1.rows, 1);
        // And what does reach it is centred exactly: log(1), log(3) minus their mean.
        let m = (1.0f64.ln() + 3.0f64.ln()) / 2.0;
        assert_eq!(st.means1[0], m);
        assert_eq!(st.y1.row(0), &[(1.0f64.ln() - m), (3.0f64.ln() - m)]);
    }

    /// A *partly* missing taxon does reach it, and its `NA` keeps R's payload
    /// through the subtraction.
    #[test]
    fn a_partly_missing_taxa_keeps_na_real_through_centering() {
        let na = ancombc2_core::reduce::na_real();
        let data = counts(1, 3, &[na, 4.0, 9.0]);
        let agg = data.clone();
        let st = preprocess_stages(&PreprocessInputs {
            data: &data,
            aggregate: &agg,
            struc_zero: false,
            group_indicator: None,
            n_groups: 0,
            prv_cut: 0.0,
            lib_cut: 0.0,
            pseudo: 0.0,
            neg_lb: false,
            red: &F64Reductions,
        })
        .unwrap();
        assert_eq!(st.tax_keep1, vec![1]);
        let m = (4.0f64.ln() + 9.0f64.ln()) / 2.0;
        assert_eq!(st.means1[0], m, "na.rm skips the NA");
        // `NA_real_ - m` is `NA_real_` in R and a fresh NaN on the hardware.
        assert_eq!(
            st.y1.get(0, 0).to_bits(),
            ancombc2_core::reduce::QUIETED_NA_REAL_BITS
        );
        assert_ne!(st.y1.get(0, 0).to_bits(), f64::NAN.to_bits());
    }

    /// `.get_struc_zero`: an `NA` count is an *absence* for the group screen, the
    /// opposite of prevalence.
    #[test]
    fn an_na_count_is_an_absence_for_the_group_screen() {
        let na = ancombc2_core::reduce::na_real();
        // One taxon, two groups of one sample each: sample 1 present, sample 2 NA.
        // Group 1: present 1/1. Group 2: present 0/1, but its *observed* count is 0,
        // so samp_size is 0 and the lower bound is NaN.
        // Two taxa: taxon 0 is present in group 1 and NA in group 2, taxon 1 is
        // present in both. Group 1 is sample 1, group 2 is sample 2.
        let data = counts(2, 2, &[1.0, na, 7.0, 7.0]);
        let agg = data.clone();
        let g = vec![1.0, 0.0, 0.0, 1.0]; // 2 x 2, column-major
        let st = preprocess_stages(&PreprocessInputs {
            data: &data,
            aggregate: &agg,
            struc_zero: true,
            group_indicator: Some(&g),
            n_groups: 2,
            prv_cut: 0.0,
            lib_cut: 0.0,
            pseudo: 0.0,
            neg_lb: false,
            red: &F64Reductions,
        })
        .unwrap();
        assert_eq!(st.group_prevalence[0], 1.0, "group 1 has taxon 0");
        assert_eq!(st.group_prevalence[1], 0.0, "group 2 sees an NA as absent");
        assert_eq!(st.group_size[0], 1.0);
        assert_eq!(st.group_size[1], 0.0, "the NA is not an observed sample");
        assert!(
            st.group_lower[1].is_nan(),
            "dividing by a group size of 0 gives NaN"
        );
        // `zero_ind = (p_hat == 0)` flags group 2 even though neg_lb is off.
        assert_eq!(st.zero_ind[0], 0);
        assert_eq!(st.zero_ind[1], 1);
        // Taxon 1 is present in both groups and survives.
        assert_eq!(st.zero_ind[2], 0);
        assert_eq!(st.zero_ind[3], 0);
        // So taxon 0 is dropped by the screen and taxon 1 kept -- and the screen
        // runs on the *aggregate* table, so the second pass's indices come from it.
        assert_eq!(st.zero_keep, vec![2]);
        // `tax_keep2` is subset-relative, as `.data_core` reports it: taxon 2 is the
        // *first* row of the one-row structural-zero subset, so it is 1 there and 2 in
        // `aggregate`'s coordinates. Both are reported, because both are wanted.
        assert_eq!(
            st.tax_keep2,
            vec![1],
            "subset-relative, as the reference reports it"
        );
        assert_eq!(
            st.tax_keep2_absolute,
            vec![2],
            "and absolute, which a caller wants"
        );
        assert_eq!(st.o2.rows, 1);
    }

    /// Every taxon flagged is the reference's "All taxa contain structural zeros",
    /// reported by name rather than as an empty result.
    #[test]
    fn a_fully_flagged_table_stops_with_the_reference_s_message() {
        let data = counts(1, 2, &[1.0, 0.0]);
        let agg = data.clone();
        let g = vec![1.0, 0.0, 0.0, 1.0];
        let err = preprocess_stages(&PreprocessInputs {
            data: &data,
            aggregate: &agg,
            struc_zero: true,
            group_indicator: Some(&g),
            n_groups: 2,
            prv_cut: 0.0,
            lib_cut: 0.0,
            pseudo: 0.0,
            neg_lb: false,
            red: &F64Reductions,
        })
        .unwrap_err();
        assert!(err.to_string().contains("structural zeros"), "got {err}");
    }

    /// The screen reads the aggregate table, not `data`.
    ///
    /// The two have the same shape and different contents here: taxon 0 is absent in
    /// the aggregate and present in `data`. Screening `data` would keep taxon 0 and
    /// the second pass would then return a different row -- silently, because every
    /// shape still matches.
    #[test]
    fn the_structural_zero_screen_reads_the_aggregate_table() {
        let data = counts(2, 2, &[5.0, 5.0, 5.0, 5.0]); // taxon 0 present
        let agg = counts(2, 2, &[5.0, 0.0, 5.0, 5.0]); // taxon 0 absent in group 2
        let g = vec![1.0, 0.0, 0.0, 1.0];
        let st = preprocess_stages(&PreprocessInputs {
            data: &data,
            aggregate: &agg,
            struc_zero: true,
            group_indicator: Some(&g),
            n_groups: 2,
            prv_cut: 0.0,
            lib_cut: 0.0,
            pseudo: 0.0,
            neg_lb: false,
            red: &F64Reductions,
        })
        .unwrap();
        assert_eq!(st.zero_ind[0], 0);
        assert_eq!(
            st.zero_ind[1], 1,
            "absent from group 2 in the aggregate table"
        );
        assert_eq!(st.zero_keep, vec![2]);
        // The first pass keeps both taxa: it screens on `data`, where both are
        // present everywhere.
        assert_eq!(st.tax_keep1, vec![1, 2]);
        // The second keeps one, because the screen removed taxon 0.
        assert_eq!(st.tax_keep2, vec![1]);
        assert_eq!(st.tax_keep2_absolute, vec![2]);
        assert_eq!(st.o1.rows, 2);
        assert_eq!(st.o2.rows, 1);
    }

    /// `neg_lb` is a wider net: it flags a taxon that is merely rare.
    #[test]
    fn neg_lb_flags_a_taxon_that_is_merely_rare() {
        // One taxon, one group of 4 samples, present in 1. p = 0.25,
        // lower = 0.25 - 1.96 * sqrt(0.25*0.75/4) = 0.25 - 0.212 < 0.
        let data = counts(2, 4, &[5.0, 0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 4.0]);
        let agg = data.clone();
        let g = vec![1.0, 1.0, 1.0, 1.0];
        let off = preprocess_stages(&PreprocessInputs {
            data: &data,
            aggregate: &agg,
            struc_zero: true,
            group_indicator: Some(&g),
            n_groups: 1,
            prv_cut: 0.0,
            lib_cut: 0.0,
            pseudo: 0.0,
            neg_lb: false,
            red: &F64Reductions,
        })
        .unwrap();
        assert_eq!(off.zero_ind[0], 0, "p_hat != 0 so the plain test misses it");
        assert!(
            off.group_lower[0] <= 0.0,
            "the lower bound is what neg_lb reads"
        );
        // Taxon 1 is present in every sample, so its lower bound is positive.
        assert_eq!(off.zero_ind[1], 0);
        assert!(off.group_lower[1] > 0.0);

        let on = preprocess_stages(&PreprocessInputs {
            data: &data,
            aggregate: &agg,
            struc_zero: true,
            group_indicator: Some(&g),
            n_groups: 1,
            prv_cut: 0.0,
            lib_cut: 0.0,
            pseudo: 0.0,
            neg_lb: true,
            red: &F64Reductions,
        })
        .unwrap();
        assert_eq!(on.zero_ind[0], 1, "neg_lb catches it");
        assert_eq!(on.zero_keep, vec![2], "and only that taxon is dropped");
    }
}

// ---------------------------------------------------------------------------
// The preprocessing stage probe
// ---------------------------------------------------------------------------

/// Which buffer a stage's payload lives in.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum StageKind {
    /// Doubles, in the caller's `f64` buffer.
    F64,
    /// Integers, in the caller's `i32` buffer.
    I32,
}

/// The stage probe's input, in the C shim's field order.
///
/// Every array is borrowed for the call. The group indicator is `n_samp x n_groups`
/// doubles in R's column-major order; it is absent when `n_groups` is zero.
#[repr(C)]
pub struct RawPreprocess {
    /// `feature_table`, `n_tax x n_samp`.
    pub data: *const f64,
    /// `feature_table_aggregate`, `n_tax x n_samp`.
    pub aggregate: *const f64,
    /// The `n_samp x n_groups` indicator `G`, or null.
    pub group_indicator: *const f64,
    /// `nrow(feature_table)`.
    pub n_tax: i64,
    /// `ncol(feature_table)`.
    pub n_samp: i64,
    /// `nlevels(factor(meta_data[, group]))`, or 0 when there is no group.
    pub n_groups: i64,
    /// `struc_zero`.
    pub struc_zero: i32,
    /// `prv_cut`.
    pub prv_cut: f64,
    /// `lib_cut`.
    pub lib_cut: f64,
    /// `pseudo`.
    pub pseudo: f64,
    /// `neg_lb`.
    pub neg_lb: i32,
}

/// One stage array on its way out: **offsets** into the caller's buffers.
///
/// # Offsets, not pointers
///
/// An earlier version of this handed `&Vec<f64>` pointers back and let the C shim read
/// them after the Rust function returned. The Rust function's locals are dropped on
/// return, so every pointer dangled: the values came back as whatever the allocator
/// had put there next, the first two elements of every stage read as
/// plausible-looking garbage while the rest were correct, and nothing crashed, nothing
/// was `NULL`, and the shapes were all right.
///
/// So the values are copied into buffers the caller allocated and only the offsets come
/// back. Same rule as `ancombc2_rb_emit`: **nothing that crosses the boundary points
/// into Rust memory.**
#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub struct RawStage {
    /// Where the stage's `n` doubles start in the caller's `f64` buffer, or `-1`.
    pub f64_at: i64,
    /// Where the stage's `n` integers start in the caller's `i32` buffer, or `-1`.
    pub i32_at: i64,
    /// `n_rows`, or 0 for a vector.
    pub n_rows: i64,
    /// `n_cols`, or 0 for a vector.
    pub n_cols: i64,
    /// `n_rows * n_cols`, or the vector's length.
    pub n: i64,
}

/// Every preprocessing stage, in the order [`crate::output::preprocess_stages`]
/// builds it.
///
/// The field order is what `init.c` reads positionally, and
/// `preprocess_stages_have_the_documented_shapes` pins it.
#[repr(C)]
pub struct RawStages {
    /// `.data_core`'s `prevalence`, over `data`.
    pub prevalence1: RawStage,
    /// The first pass's `tax_keep`, 1-based.
    pub tax_keep1: RawStage,
    /// `colSums` of the taxon-filtered table, over every sample.
    pub lib_size1: RawStage,
    /// The first pass's `samp_keep`, 1-based.
    pub samp_keep1: RawStage,
    /// `O1 = data + pseudo`.
    pub o1: RawStage,
    /// `log(O1)` with `is.infinite` mapped to `NA_real_`.
    pub log1: RawStage,
    /// `rowMeans(log1, na.rm = TRUE)`, R's own reduction.
    pub means1: RawStage,
    /// `log1 - means1`.
    pub y1: RawStage,
    /// The second pass's `tax_keep`, as `.data_core` reports it: positions within the
    /// structural-zero-retained table.
    pub tax_keep2: RawStage,
    /// The same set in `aggregate`'s coordinates.
    pub tax_keep2_absolute: RawStage,
    /// The second pass's `samp_keep`. Reused from the first.
    pub samp_keep2: RawStage,
    /// `O2 = aggregate_data + pseudo`.
    pub o2: RawStage,
    /// `log(O2)` with `is.infinite` mapped to `NA_real_`.
    pub log2: RawStage,
    /// `rowMeans(log2, na.rm = TRUE)`.
    pub means2: RawStage,
    /// `log2 - means2`.
    pub y2: RawStage,
    /// `which(all(zero_ind[, -1] == FALSE))`, 1-based.
    pub zero_keep: RawStage,
    /// `sweep(present %*% G, 2, n_g, "/")`.
    pub group_prevalence: RawStage,
    /// `rowSums(!is.na(feature_mat)) %*% G`.
    pub group_size: RawStage,
    /// `p_hat - 1.96 * sqrt(p_hat * (1 - p_hat) / samp_size)`.
    pub group_lower: RawStage,
    /// The flags, as 0/1.
    pub zero_ind: RawStage,
}

/// **Preprocessing stage probe.** Run every preprocessing stage in the bridge and
/// hand the arrays back.
///
/// It computes the same arithmetic `ancombc2()` will, on the same inputs, with
/// nothing held back — and it hands back each stage separately, which is the point.
/// A whole-result comparison says the numbers differ; this says *which stage*.
///
/// The centring reduction is R's own (`RBackedReductions`), so `means1`, `means2`
/// and the centred matrices are what R would compute. That is not circular: R's
/// reduction is being reused, not its answer, and the arrays are compared against R
/// expressions run independently by `scripts/check_preprocess_stages.R`.
///
/// # Safety
///
/// `raw` must describe valid buffers, and the `out` pointers must be null or valid
/// for `stage_count` and `slot_count` entries respectively.
#[no_mangle]
pub unsafe extern "C" fn ancombc2_rb_preprocess_probe(
    raw: *const RawPreprocess,
    f64_out: *mut f64,
    i32_out: *mut i32,
    out: *mut RawStages,
    slot_count: i64,
    out_f64_len: *mut i64,
    out_i32_len: *mut i64,
    out_n_rows: *mut i64,
    out_n_cols: *mut i64,
    out_has_group: *mut i32,
    out_names: *mut i8,
    out_name_offsets: *mut i64,
    out_names_len: *mut i64,
) -> i32 {
    let r = match raw.as_ref() {
        Some(r) => r,
        None => {
            set_error("the preprocess input pointer was null");
            return RB_ERR_ARGS;
        }
    };
    match catch_unwind(AssertUnwindSafe(|| {
        let u = |r: Result<usize, BridgeError>| -> Result<usize, (i32, String)> { to_result(r) };
        let n_tax = u(as_usize(r.n_tax, "n_tax"))?;
        let n_samp = u(as_usize(r.n_samp, "n_samp"))?;
        let n_groups = u(as_usize(r.n_groups, "n_groups"))?;
        let n = n_tax * n_samp;
        let data = to_result(slice::<f64>(r.data, n, "feature_table"))?.to_vec();
        let agg = to_result(slice::<f64>(r.aggregate, n, "feature_table_aggregate"))?.to_vec();
        let g = if n_groups == 0 {
            Vec::new()
        } else {
            let want = n_samp * n_groups;
            let got = r.group_indicator;
            if got.is_null() {
                return Err((
                    RB_ERR_ARGS,
                    "n_groups > 0 but the group indicator is null".into(),
                ));
            }
            to_result(slice::<f64>(got, want, "group_indicator"))?.to_vec()
        };

        // R's buffers are column-major and `CountMatrix` is row-major, so the two
        // tables are transposed on the way in. Passing them straight to
        // `CountMatrix::new` builds a transposed table whose shape is still valid,
        // which is how this produced a `lib_size1` that was the column sum of the
        // wrong table and nothing else out of place.
        let data_cm = CountMatrix::from_column_major(n_tax, n_samp, &data)
            .map_err(|e| (RB_ERR_ARGS, e.to_string()))?;
        let agg_cm = CountMatrix::from_column_major(n_tax, n_samp, &agg)
            .map_err(|e| (RB_ERR_ARGS, e.to_string()))?;

        let red = RBackedReductions::linked();
        let st = to_result(crate::output::preprocess_stages(
            &crate::output::PreprocessInputs {
                data: &data_cm,
                aggregate: &agg_cm,
                struc_zero: r.struc_zero != 0,
                group_indicator: if n_groups == 0 { None } else { Some(&g) },
                n_groups,
                prv_cut: r.prv_cut,
                lib_cut: r.lib_cut,
                pseudo: r.pseudo,
                neg_lb: r.neg_lb != 0,
                red: &red,
            },
        ))?;

        // The stages, in the field order of `RawStages`. Written out rather than
        // generated so the order is readable and `init.c`'s table can be compared
        // with it by eye as well as by the length check.
        // The shape follows the payload, not the other way round. With
        // `struc_zero` off there is no group screen, so those stages are empty --
        // and an empty array that still claims `n_tax x n_groups` makes the C shim
        // set a `dim` of the wrong length on a zero-length vector, which is an R
        // error ("dims [product 6] do not match the length of object [0]") rather
        // than an empty result.
        let empty = |b: bool| usize::from(!b);
        let vec_stage = |n: usize| (n, n, 0usize);
        let int_stage = |n: usize| (n, n, 0usize);
        let mat_stage = |m: &ancombc2_core::workspace::RMatrix| {
            (
                m.data.len(),
                empty(m.data.is_empty()) * m.rows,
                empty(m.data.is_empty()) * m.cols,
            )
        };
        let group_stage = |n: usize| (n, empty(n == 0) * n_tax, empty(n == 0) * n_groups);

        // Name, shape and buffer kind for each stage, in the order the C shim's table
        // uses. The names are carried here and not in the shim's list so the two
        // cannot drift silently: the shim compares its own names against these.
        let stages: [(&str, StageKind, (usize, usize, usize)); 20] = [
            (
                "prevalence1",
                StageKind::F64,
                vec_stage(st.prevalence1.len()),
            ),
            ("tax_keep1", StageKind::I32, int_stage(st.tax_keep1.len())),
            ("lib_size1", StageKind::F64, vec_stage(st.lib_size1.len())),
            ("samp_keep1", StageKind::I32, int_stage(st.samp_keep1.len())),
            ("O1", StageKind::F64, mat_stage(&st.o1)),
            ("log1", StageKind::F64, mat_stage(&st.log1)),
            ("means1", StageKind::F64, vec_stage(st.means1.len())),
            ("y1", StageKind::F64, mat_stage(&st.y1)),
            ("tax_keep2", StageKind::I32, int_stage(st.tax_keep2.len())),
            (
                "tax_keep2_absolute",
                StageKind::I32,
                int_stage(st.tax_keep2_absolute.len()),
            ),
            ("samp_keep2", StageKind::I32, int_stage(st.samp_keep2.len())),
            ("O2", StageKind::F64, mat_stage(&st.o2)),
            ("log2", StageKind::F64, mat_stage(&st.log2)),
            ("means2", StageKind::F64, vec_stage(st.means2.len())),
            ("y2", StageKind::F64, mat_stage(&st.y2)),
            ("zero_keep", StageKind::I32, int_stage(st.zero_keep.len())),
            (
                "group_prevalence",
                StageKind::F64,
                group_stage(st.group_prevalence.len()),
            ),
            (
                "group_size",
                StageKind::F64,
                group_stage(st.group_size.len()),
            ),
            (
                "group_lower",
                StageKind::F64,
                group_stage(st.group_lower.len()),
            ),
            ("zero_ind", StageKind::I32, group_stage(st.zero_ind.len())),
        ];

        // Lay the stages out in the caller's buffers. Offsets first, so the totals are
        // known and reported before the shim reads anything.
        //
        // The two totals are reported *here*, before the copy, because the shim sizes
        // its buffers from them. Leaving this out is not a compile error and not a
        // crash: the shim allocated a one-element buffer, this wrote a hundred-odd
        // doubles past the end of it, and the symptom was a heap that R's collector
        // later tripped over with "recursive gc invocation".
        let mut descriptors = [RawStage {
            f64_at: -1,
            i32_at: -1,
            n_rows: 0,
            n_cols: 0,
            n: 0,
        }; 20];
        let mut f64_at = 0i64;
        let mut i32_at = 0i64;
        for (i, (_, kind, (n, rows, cols))) in stages.iter().enumerate() {
            descriptors[i].n = *n as i64;
            descriptors[i].n_rows = *rows as i64;
            descriptors[i].n_cols = *cols as i64;
            match kind {
                StageKind::F64 => {
                    descriptors[i].f64_at = f64_at;
                    f64_at += *n as i64;
                }
                StageKind::I32 => {
                    descriptors[i].i32_at = i32_at;
                    i32_at += *n as i64;
                }
            }
        }
        // Written on both passes: the shim reads them from the sizing pass and
        // re-checks them on the filling pass, so the two layouts cannot disagree.
        unsafe {
            if !out_f64_len.is_null() {
                *out_f64_len = f64_at;
            }
            if !out_i32_len.is_null() {
                *out_i32_len = i32_at;
            }
        }
        // The copy. Every byte that crosses the boundary lives in a buffer the caller
        // allocated, so nothing here is a pointer into Rust memory that outlives this
        // function.
        let mut w64 = 0i64;
        let mut w32 = 0i64;
        for (i, (name, kind, (n, _, _))) in stages.iter().enumerate() {
            let n = *n as i64;
            let name: &str = name;
            match *kind {
                StageKind::F64 => {
                    let src: &[f64] = match name {
                        "prevalence1" => &st.prevalence1,
                        "lib_size1" => &st.lib_size1,
                        "O1" => &st.o1.data,
                        "log1" => &st.log1.data,
                        "means1" => &st.means1,
                        "y1" => &st.y1.data,
                        "O2" => &st.o2.data,
                        "log2" => &st.log2.data,
                        "means2" => &st.means2,
                        "y2" => &st.y2.data,
                        "group_prevalence" => &st.group_prevalence,
                        "group_size" => &st.group_size,
                        "group_lower" => &st.group_lower,
                        _ => return Err((RB_ERR_ARGS, "a double stage had no source".into())),
                    };
                    if src.len() as i64 != n {
                        return Err((
                            RB_ERR_ARGS,
                            format!("double stage {name} has {} values, expected {n}", src.len()),
                        ));
                    }
                    let at = descriptors[i].f64_at;
                    let (_, rows, cols) = stages[i].2;
                    if !f64_out.is_null() {
                        unsafe {
                            if cols == 0 {
                                // A vector: element for element.
                                std::ptr::copy_nonoverlapping(
                                    src.as_ptr(),
                                    f64_out.add(at as usize),
                                    src.len(),
                                );
                            } else {
                                // A matrix. `RMatrix` and the group tables are
                                // row-major; R is column-major; so the stage crosses in
                                // R's order. Same on the way in, and the same mistake in
                                // the other direction, so both are written out.
                                if rows * cols != src.len() {
                                    return Err((
                                        RB_ERR_ARGS,
                                        format!(
                                            "stage {name} claims {rows} x {cols} but carries \
                                             {} values",
                                            src.len()
                                        ),
                                    ));
                                }
                                for r in 0..rows {
                                    for c in 0..cols {
                                        *f64_out.add(at as usize + r + c * rows) =
                                            src[r * cols + c];
                                    }
                                }
                            }
                        };
                    }
                    w64 += n;
                }
                StageKind::I32 => {
                    let src: &[i64] = match name {
                        "tax_keep1" => &st.tax_keep1,
                        "samp_keep1" => &st.samp_keep1,
                        "tax_keep2" => &st.tax_keep2,
                        "tax_keep2_absolute" => &st.tax_keep2_absolute,
                        "samp_keep2" => &st.samp_keep2,
                        "zero_keep" => &st.zero_keep,
                        "zero_ind" => &st.zero_ind,
                        _ => {
                            return Err((RB_ERR_ARGS, "an integer stage had no source".into()));
                        }
                    };
                    if src.len() as i64 != n {
                        return Err((
                            RB_ERR_ARGS,
                            format!(
                                "integer stage {name} has {} values, expected {n}",
                                src.len()
                            ),
                        ));
                    }
                    let at = descriptors[i].i32_at;
                    let (_, rows, cols) = stages[i].2;
                    if !i32_out.is_null() {
                        unsafe {
                            if cols == 0 {
                                // A vector: element for element.
                                for (k, &v) in src.iter().enumerate() {
                                    *i32_out.add(at as usize + k) = v as i32;
                                }
                            } else {
                                // A matrix, and the same transpose the double stages get.
                                // `zero_ind` is the only integer matrix today, and it was
                                // the only stage that came back transposed -- which showed
                                // up as the (1,2) and (2,1) flags being swapped and nothing
                                // else wrong.
                                if rows * cols != src.len() {
                                    return Err((
                                        RB_ERR_ARGS,
                                        format!(
                                            "stage {name} claims {rows} x {cols} but carries \
                                             {} values",
                                            src.len()
                                        ),
                                    ));
                                }
                                for r in 0..rows {
                                    for c in 0..cols {
                                        *i32_out.add(at as usize + r + c * rows) =
                                            src[r * cols + c] as i32;
                                    }
                                }
                            }
                        }
                    }
                    w32 += n;
                }
            }
        }
        if w64 != f64_at || w32 != i32_at {
            return Err((
                RB_ERR_ARGS,
                "the stage table and the bytes written disagree on a length".into(),
            ));
        }
        let want_slots = descriptors.len() as i64;
        if slot_count < want_slots {
            return Err((
                RB_ERR_ARGS,
                format!(
                    "the probe produced {want_slots} stages and the caller reserved \\
                     {slot_count}"
                ),
            ));
        }
        // The stage names, so the shim can check its own table against them.
        unsafe {
            if !out_names.is_null() {
                let mut at = 0i64;
                for (i, (name, _, _)) in stages.iter().enumerate() {
                    // The offset is recorded *before* the bytes, so it is the start
                    // of this name and the next offset is its end. Recording it after
                    // would make every offset one name too late, which the shim's
                    // monotone-partition check cannot see -- it would just report the
                    // wrong name for every stage.
                    *out_name_offsets.add(i) = at;
                    let bytes = name.as_bytes();
                    for (k, b) in bytes.iter().enumerate() {
                        *out_names.add(at as usize + k) = *b as i8;
                    }
                    at += bytes.len() as i64;
                }
                *out_name_offsets.add(stages.len()) = at;
                *out_names_len = at;
            }
        }

        unsafe {
            // `RawStages` is `#[repr(C)]` over nineteen `RawStage` fields and no
            // padding, so a `RawStage` slice is exactly the caller's buffer. The cast
            // is written out rather than hidden because that is the assumption it
            // relies on, and `the_nineteen_stages_are_in_the_documented_order` checks
            // the size.
            let slots = out as *mut RawStage;
            std::ptr::copy_nonoverlapping(descriptors.as_ptr(), slots, descriptors.len());
            if !out_n_rows.is_null() {
                *out_n_rows = n_tax as i64;
            }
            if !out_n_cols.is_null() {
                *out_n_cols = n_samp as i64;
            }
            if !out_has_group.is_null() {
                *out_has_group = i32::from(n_groups > 0);
            }
        }
        Ok(())
    })) {
        Ok(Ok(())) => RB_OK,
        Ok(Err((code, msg))) => {
            set_error(msg);
            code
        }
        Err(payload) => {
            set_error(format!(
                "ancombc2-rbridge panicked in preprocess_probe: {}",
                panic_message(&payload)
            ));
            RB_ERR_PANIC
        }
    }
}

#[cfg(test)]
mod preprocess_probe_layout_tests {
    //! `RawStages` is nineteen `RawStage` fields read positionally by
    //! `r/ANCOMBC/src/init.c`. Nothing in either type system ties the two orders
    //! together, and a swap would not crash -- it would return `y1` where `log1`
    //! belongs and every stage comparison would report a difference with no hint
    //! that the plumbing was wrong.
    //!
    //! These offsets are the specification the C side is written against.

    use super::{RawPreprocess, RawStage, RawStages};
    use core::mem::{align_of, offset_of, size_of};

    #[test]
    fn the_stage_slot_is_the_size_the_shim_assumes() {
        // Two 8-byte pointers, four i64s. No padding: every field is 8-aligned.
        assert_eq!(size_of::<RawStage>(), 40);
        assert_eq!(align_of::<RawStage>(), 8);
        assert_eq!(offset_of!(RawStage, f64_at), 0);
        assert_eq!(offset_of!(RawStage, i32_at), 8);
        assert_eq!(offset_of!(RawStage, n_rows), 16);
        assert_eq!(offset_of!(RawStage, n_cols), 24);
        assert_eq!(offset_of!(RawStage, n), 32);
    }

    #[test]
    fn the_twenty_stages_are_in_the_documented_order() {
        assert_eq!(size_of::<RawStages>(), 20 * 40);
        // A `RawStage` slice over the caller's buffer is only the same bytes if the
        // struct is a plain sequence of slots with no padding of its own.
        assert_eq!(align_of::<RawStages>(), 8);
        let names = [
            "prevalence1",
            "tax_keep1",
            "lib_size1",
            "samp_keep1",
            "o1",
            "log1",
            "means1",
            "y1",
            "tax_keep2",
            "tax_keep2_absolute",
            "samp_keep2",
            "o2",
            "log2",
            "means2",
            "y2",
            "zero_keep",
            "group_prevalence",
            "group_size",
            "group_lower",
            "zero_ind",
        ];
        let offsets = [
            offset_of!(RawStages, prevalence1),
            offset_of!(RawStages, tax_keep1),
            offset_of!(RawStages, lib_size1),
            offset_of!(RawStages, samp_keep1),
            offset_of!(RawStages, o1),
            offset_of!(RawStages, log1),
            offset_of!(RawStages, means1),
            offset_of!(RawStages, y1),
            offset_of!(RawStages, tax_keep2),
            offset_of!(RawStages, tax_keep2_absolute),
            offset_of!(RawStages, samp_keep2),
            offset_of!(RawStages, o2),
            offset_of!(RawStages, log2),
            offset_of!(RawStages, means2),
            offset_of!(RawStages, y2),
            offset_of!(RawStages, zero_keep),
            offset_of!(RawStages, group_prevalence),
            offset_of!(RawStages, group_size),
            offset_of!(RawStages, group_lower),
            offset_of!(RawStages, zero_ind),
        ];
        assert_eq!(
            offsets.len(),
            names.len(),
            "the offsets list and the names list have drifted apart"
        );
        for (i, (off, name)) in offsets.iter().zip(names.iter()).enumerate() {
            assert_eq!(*off, i * 40, "{name} is at slot {i}");
        }
    }

    #[test]
    fn the_preprocess_input_layout_is_pinned() {
        // Six 8-byte fields then six scalars.
        assert_eq!(offset_of!(RawPreprocess, data), 0);
        assert_eq!(offset_of!(RawPreprocess, aggregate), 8);
        assert_eq!(offset_of!(RawPreprocess, group_indicator), 16);
        assert_eq!(offset_of!(RawPreprocess, n_tax), 24);
        assert_eq!(offset_of!(RawPreprocess, n_samp), 32);
        assert_eq!(offset_of!(RawPreprocess, n_groups), 40);
        // `i32` then three `f64`: 4 bytes of padding before `prv_cut`.
        assert_eq!(offset_of!(RawPreprocess, struc_zero), 48);
        assert_eq!(offset_of!(RawPreprocess, prv_cut), 56);
        assert_eq!(offset_of!(RawPreprocess, lib_cut), 64);
        assert_eq!(offset_of!(RawPreprocess, pseudo), 72);
        assert_eq!(offset_of!(RawPreprocess, neg_lb), 80);
        assert_eq!(size_of::<RawPreprocess>(), 88);
    }
}

#[cfg(test)]
mod lib_size_tests {
    //! `lib_size1` came back as `13, 3, 14, 2` where the reference's `colSums` gives
    //! `6, 9, 11, 6` -- for a table whose every other stage matched. Those numbers are
    //! integers in a range no column sum of that table takes, so the suspicion is the
    //! buffer offset rather than the arithmetic, and this pins which side it is on.

    use crate::output::{preprocess_stages, PreprocessInputs};
    use ancombc2_core::preprocess::CountMatrix;
    use ancombc2_core::reduce::F64Reductions;

    fn stages_for(data: Vec<f64>, n_tax: usize, n_samp: usize) -> crate::output::PreprocessStages {
        let d = CountMatrix::new(n_tax, n_samp, data).unwrap();
        preprocess_stages(&PreprocessInputs {
            data: &d,
            aggregate: &d,
            struc_zero: false,
            group_indicator: None,
            n_groups: 0,
            prv_cut: 0.0,
            lib_cut: 0.0,
            pseudo: 0.0,
            neg_lb: false,
            red: &F64Reductions,
        })
        .unwrap()
    }

    #[test]
    fn the_library_size_is_the_column_sum_of_the_retained_taxa() {
        // [[4,0,9,1],[0,7,0,3],[2,2,2,2]] -> columns sum to 6, 9, 11, 6.
        let st = stages_for(
            vec![4.0, 0.0, 9.0, 1.0, 0.0, 7.0, 0.0, 3.0, 2.0, 2.0, 2.0, 2.0],
            3,
            4,
        );
        assert_eq!(st.tax_keep1, vec![1, 2, 3]);
        assert_eq!(st.lib_size1, vec![6.0, 9.0, 11.0, 6.0]);
        assert_eq!(st.samp_keep1, vec![1, 2, 3, 4]);
        // And the retained table, which is what those sums are over.
        assert_eq!(
            st.o1.data,
            vec![4.0, 0.0, 9.0, 1.0, 0.0, 7.0, 0.0, 3.0, 2.0, 2.0, 2.0, 2.0]
        );
    }

    /// An `NA` leaves the library size entirely: `na.rm = TRUE` skips it.
    #[test]
    fn an_na_does_not_enter_the_library_size() {
        let na = ancombc2_core::reduce::na_real();
        let st = stages_for(vec![na, 5.0, 7.0, 1.0, 0.0, 9.0], 2, 3);
        // Taxon 1 is [NA, 5, 7] and taxon 2 is [1, 0, 9]; both clear prv_cut = 0, and
        // `na.rm` skips the NA, so the first column sums to 1 and not to NA.
        assert_eq!(st.tax_keep1, vec![1, 2]);
        assert_eq!(st.lib_size1, vec![1.0, 5.0, 16.0]);
    }
}

#[cfg(test)]
mod centered_bits_tests {
    //! The stage comparison reported `y1` as `fff00000000007a2` where R gives
    //! `7ff80000000007a2`, at every position that `log(0)` had turned into `NA`.
    //! Both are "an `NA` with the sign bit set" to the eye and neither is right.

    use crate::output::{preprocess_stages, PreprocessInputs};
    use ancombc2_core::preprocess::CountMatrix;
    use ancombc2_core::reduce::{na_real, F64Reductions, NA_REAL_BITS, QUIETED_NA_REAL_BITS};

    #[test]
    fn a_centred_na_carries_the_signed_na_pattern() {
        // [[4,0,9,1],[0,7,0,3],[2,2,2,2]] with pseudo 0, so the zero at (2,1) in R's
        // 1-based numbering -- row 2, column 1 -- becomes log(0) = -Inf -> NA_real_,
        // and then the centring has to carry it.
        let d = CountMatrix::new(
            3,
            4,
            vec![4.0, 0.0, 9.0, 1.0, 0.0, 7.0, 0.0, 3.0, 2.0, 2.0, 2.0, 2.0],
        )
        .unwrap();
        let st = preprocess_stages(&PreprocessInputs {
            data: &d,
            aggregate: &d,
            struc_zero: false,
            group_indicator: None,
            n_groups: 0,
            prv_cut: 0.0,
            lib_cut: 0.0,
            pseudo: 0.0,
            neg_lb: false,
            red: &F64Reductions,
        })
        .unwrap();

        // `O + pseudo` is the zero, unchanged.
        assert_eq!(st.o1.get(1, 0), 0.0);
        // `log`, then `is.infinite` -> NA_real_ with R's exact payload.
        assert_eq!(st.log1.get(1, 0).to_bits(), NA_REAL_BITS);
        // `log1 - means1`: R quiets the NaN and keeps the payload.
        assert_eq!(
            st.y1.get(1, 0).to_bits(),
            QUIETED_NA_REAL_BITS,
            "got {:016x}",
            st.y1.get(1, 0).to_bits()
        );
        // And the *other* cells of that row are finite, so this is not a whole-array
        // mistake.
        assert!(st.y1.get(1, 1).is_finite());
        assert!(st.y1.get(1, 3).is_finite());
        assert_ne!(st.y1.get(1, 0).to_bits(), na_real().to_bits());
    }
}

// ---------------------------------------------------------------------------
// The least-squares probe
// ---------------------------------------------------------------------------

/// Inputs for the fit probe: the same `x` and `Y` `.lm_fit_all` receives.
///
/// `y` is `n_taxa x n_samp` row-major and already theta-adjusted, exactly as
/// `.iter_mle` hands it over after `sweep(y, 2, theta, "-")`.
#[repr(C)]
pub struct RawFit {
    /// `x`, `n_samp x p`, column-major.
    pub x: *const f64,
    /// `Y`, `n_taxa x n_samp`, row-major.
    pub y: *const f64,
    /// `nrow(x)`.
    pub n_samp: i64,
    /// `ncol(x)`.
    pub p: i64,
    /// `nrow(Y)`.
    pub n_taxa: i64,
    /// `1` per `(taxon, sample)` where the response is usable, `n_taxa * n_samp`.
    ///
    /// `.lm_fit_all`'s `use` is `is.finite(Y) & rep(x_ok, each row)`, and the
    /// grouping keys on that mask, so the pattern a taxon belongs to is decided by
    /// the response's finiteness as much as by the design's completeness. Passing
    /// the mask in rather than recomputing it keeps the probe from becoming a second
    /// implementation of the thing it is measuring.
    pub observed: *const u8,
}

/// **BLAS calibration.** Solve one least-squares problem the way R's `Cdqrls`
/// does, under a named BLAS, so R can compare the bits with its own `Cdqrls` and
/// learn which BLAS it is linked to.
///
/// `kind` is `0` for the netlib reference BLAS and `1` for OpenBLAS's Haswell/Zen
/// kernels. Nothing global is touched; see [`ancombc2_rb_set_blas`].
///
/// # Safety
///
/// `x` must hold `n * p` doubles, `y` and `resid_out` `n`, and `coef_out` `p`;
/// `rank_out` must be valid.
#[no_mangle]
pub unsafe extern "C" fn ancombc2_rb_dqrls(
    kind: i32,
    n: i64,
    p: i64,
    x: *const f64,
    y: *const f64,
    tol: f64,
    coef_out: *mut f64,
    resid_out: *mut f64,
    rank_out: *mut i32,
) -> i32 {
    use ancombc2_core::matrix::{dqrls_multi, OpenBlasHaswell, RefBlas};
    if x.is_null() || y.is_null() || coef_out.is_null() || resid_out.is_null() || rank_out.is_null()
    {
        set_error("a dqrls calibration pointer was null");
        return RB_ERR_ARGS;
    }
    if n <= 0 || p <= 0 {
        set_error("dqrls calibration needs n > 0 and p > 0");
        return RB_ERR_ARGS;
    }
    let (n, p) = (n as usize, p as usize);
    match catch_unwind(AssertUnwindSafe(|| {
        let xv = std::slice::from_raw_parts(x, n * p).to_vec();
        let yv = std::slice::from_raw_parts(y, n);
        let fit = match kind {
            0 => dqrls_multi::<RefBlas>(xv, n, p, &[yv], tol),
            #[cfg(target_arch = "x86_64")]
            1 => dqrls_multi::<OpenBlasHaswell<false>>(xv, n, p, &[yv], tol),
            _ => return Err(format!("unknown BLAS kind {kind} on this target")),
        };
        std::slice::from_raw_parts_mut(coef_out, p).copy_from_slice(&fit.coef[0]);
        std::slice::from_raw_parts_mut(resid_out, n).copy_from_slice(&fit.resid[0]);
        *rank_out = fit.rank as i32;
        Ok(())
    })) {
        Ok(Ok(())) => RB_OK,
        Ok(Err(m)) => {
            set_error(m);
            RB_ERR_ARGS
        }
        Err(_) => {
            set_error("panic in dqrls calibration");
            RB_ERR_PANIC
        }
    }
}

/// Select the BLAS semantics the least-squares fit reproduces: `0` reference,
/// `1` OpenBLAS Haswell/Zen. Called once by R after a calibration succeeded.
#[no_mangle]
pub extern "C" fn ancombc2_rb_set_blas(kind: i32) -> i32 {
    use ancombc2_core::matrix::{set_blas_kind, BlasKind};
    match kind {
        0 => set_blas_kind(BlasKind::Reference),
        #[cfg(target_arch = "x86_64")]
        1 => set_blas_kind(BlasKind::OpenBlasHaswell),
        _ => {
            set_error(format!("unknown BLAS kind {kind} on this target"));
            return RB_ERR_ARGS;
        }
    }
    RB_OK
}

/// **Least-squares probe.** Run the pipeline's own `lm_fit_all` and write its three
/// outputs into the caller's buffers.
///
/// It is the real function on the real inputs, so a difference is a difference in the
/// code that will run -- not in a copy made for the test.
///
/// The outputs are written into buffers the *caller* allocated, for the same reason as
/// `ancombc2_rb_preprocess_probe`: handing back pointers into Rust memory means the
/// values are gone before the caller can read them.
///
/// # Safety
///
/// `raw` must describe valid buffers for the counts it carries, and the three output
/// pointers must each have room for `n_taxa * p`, `n_taxa * n_samp` and `n_taxa`
/// doubles respectively.
#[no_mangle]
pub unsafe extern "C" fn ancombc2_rb_fit_probe(
    raw: *const RawFit,
    beta_out: *mut f64,
    fitted_out: *mut f64,
    dof_out: *mut f64,
    groups_out: *mut i32,
) -> i32 {
    let r = match raw.as_ref() {
        Some(r) => r,
        None => {
            set_error("the fit input pointer was null");
            return RB_ERR_ARGS;
        }
    };
    match catch_unwind(AssertUnwindSafe(|| {
        let u = |v: Result<usize, BridgeError>| -> Result<usize, (i32, String)> { to_result(v) };
        let n_samp = u(as_usize(r.n_samp, "n_samp"))?;
        let p = u(as_usize(r.p, "p"))?;
        let n_taxa = u(as_usize(r.n_taxa, "n_taxa"))?;
        if n_samp == 0 || p == 0 || n_taxa == 0 {
            return Err((RB_ERR_ARGS, "an empty fit has nothing to measure".into()));
        }
        for (what, ptr) in [("beta", beta_out), ("fitted", fitted_out), ("dof", dof_out)] {
            if ptr.is_null() {
                return Err((RB_ERR_ARGS, format!("the {what} output pointer was null")));
            }
        }
        let x = Matrix::from_vec(
            n_samp,
            p,
            to_result(slice::<f64>(r.x, n_samp * p, "x"))?.to_vec(),
        )
        .map_err(|e| (RB_ERR_ARGS, e.to_string()))?;
        let y = to_result(slice::<f64>(r.y, n_taxa * n_samp, "y"))?.to_vec();
        if r.observed.is_null() {
            return Err((RB_ERR_ARGS, "the observed mask was null".into()));
        }
        let raw_mask_true = (0..n_taxa * n_samp)
            .filter(|&k| unsafe { *r.observed.add(k) != 0 })
            .count();
        let mut obs: Vec<bool> = (0..n_taxa * n_samp)
            .map(|k| unsafe { *r.observed.add(k) != 0 })
            .collect();

        // `.lm_fit_all`'s `use` is `is.finite(Y) & rep(x_ok, each row)`, and it keys
        // the grouping on that mask, so the pattern a taxon belongs to is decided by
        // the response's finiteness as much as by the design's completeness. Folding
        // `x_ok` in here keeps the probe from becoming a second implementation of the
        // thing it is measuring.
        let x_ok: Vec<bool> = (0..n_samp)
            .map(|j| (0..p).all(|a| !x.get(j, a).is_nan()))
            .collect();
        for row in obs.chunks_mut(n_samp) {
            for (j, v) in row.iter_mut().enumerate() {
                *v = *v && x_ok[j];
            }
        }
        let group_cols: Vec<usize> = (1..p).collect();
        let cache = DesignCache::build(&obs, n_taxa, n_samp, &x, &group_cols);
        let out = lm_fit_all(&cache, &y, n_taxa, n_samp);

        unsafe {
            std::ptr::copy_nonoverlapping(out.beta.as_ptr(), beta_out, out.beta.len());
            std::ptr::copy_nonoverlapping(out.fitted.as_ptr(), fitted_out, out.fitted.len());
            std::ptr::copy_nonoverlapping(out.dof.as_ptr(), dof_out, out.dof.len());
        }

        // How the patterns came out, because a difference in `dof` is nearly always a
        // difference in which rows a group was given, and that is a two-line question
        // rather than a reading of the QR.
        //
        // Written in pairs: the number of rows in group `g`, then the number of taxa
        // in it. The reference's `dof[idx] <- n_used - rank` depends on the first and
        // the rank on the second.
        if !groups_out.is_null() {
            let mut at = 0usize;
            unsafe {
                // The dims first, then the groups: if a group's row count disagrees
                // with `n_samp` then the disagreement is upstream of the grouping and
                // reading the grouping would only restate it.
                *groups_out.add(at) = n_samp as i32;
                at += 1;
                *groups_out.add(at) = p as i32;
                at += 1;
                *groups_out.add(at) = n_taxa as i32;
                at += 1;
                *groups_out.add(at) = cache.rows.len() as i32;
                at += 1;
                *groups_out.add(at) = raw_mask_true as i32;
                at += 1;
                *groups_out.add(at) = x_ok.iter().filter(|&&b| b).count() as i32;
                at += 1;
                *groups_out.add(at) = x.has_nan() as i32;
                at += 1;
                for (g, taxa) in cache.taxa.iter().enumerate() {
                    *groups_out.add(at) = cache.rows[g].len() as i32;
                    at += 1;
                    *groups_out.add(at) = taxa.len() as i32;
                    at += 1;
                }
            }
        }

        Ok(())
    })) {
        Ok(Ok(())) => RB_OK,
        Ok(Err((code, msg))) => {
            set_error(msg);
            code
        }
        Err(payload) => {
            set_error(format!(
                "ancombc2-rbridge panicked in fit_probe: {}",
                panic_message(&payload)
            ));
            RB_ERR_PANIC
        }
    }
}
