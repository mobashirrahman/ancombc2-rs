//! R's `rowSums`/`colSums`/`rowMeans`/`colMeans` with `na.rm = TRUE`.
//!
//! # Why this is a trait and not a function
//!
//! R does not accumulate these in `double`. R's `do_colsum` (in `src/main/array.c`)
//! declares its running total as C `long double`:
//!
//! ```c
//! LONG_DOUBLE s = 0.;
//! for (R_xlen_t i = 0; i < nr; i++) {
//!     double x = REAL(x)[i + j * nr];
//!     if (!ISNAN(x)) { s += x; n++; }
//! }
//! if (n > 0) REAL(ans)[j] = (double)(s / n);
//! ```
//!
//! On x86-64 that is the x87 80-bit format: **64 bits of mantissa**, not 53. The
//! extra 11 bits change the answer, and not by a rounding-error-sized amount —
//! because a sum that cancels catastrophically in `double` simply retains the
//! small terms in `long double`. Measured on the primary profile, for
//! `y = c(1e16, 1, -1e16, 1e-17, 1)`:
//!
//! ```text
//! rowMeans   0x1.999999999999ap-2   (sum 2, /5)
//! f64 loop   0x1.999999999999ap-3   (sum 1, /5)
//! ```
//!
//! A factor of two, from one line of `+=`.
//!
//! # Why this is *not* a hand-rolled 80-bit accumulator
//!
//! Two reasons, and the second is the one that matters.
//!
//! 1. R's `long double` is platform-defined. On x86-64 Linux it is the x87 80-bit
//!    type with a 64-bit mantissa; on aarch64 Linux it is IEEE binary128; on
//!    Windows R has no `long double` at all and uses `double`. A hand-rolled
//!    64-bit-mantissa accumulator would be *right* on the first and *wrong* on the
//!    other two — it would encode one platform's answer as if it were the oracle's.
//! 2. A double-double (106-bit) compensated sum is not a substitute. It is more
//!    accurate, and "more accurate" is the wrong target: it rounds to 106 bits per
//!    step where R rounds to 64, so it disagrees with R on precisely the inputs
//!    where a naive `f64` sum also disagrees, just with smaller error.
//!
//! The correct target is "whatever `long double` does on the machine R is running
//! on". The only way to hit that by construction is to ask R.
//!
//! # What the no-R path gets
//!
//! [`F64Reductions`] sums in `f64`. It is exact for integer-valued inputs, which
//! covers the two reductions preprocessing performs on counts:
//!
//! * `rowSums(x != 0, na.rm = TRUE)` and `rowSums(!is.na(x))` — sums of 0/1;
//! * `colSums(feature_table, na.rm = TRUE)` — library sizes, sums of whole counts.
//!
//! Every partial sum is an exact integer and the total is below `2^53`, so no
//! rounding occurs at any step and the reduction order is irrelevant. Those three
//! are therefore *not* approximations, and `preprocess.rs` uses them as plain
//! functions for that reason.
//!
//! [`F64Reductions::row_means_na_rm`] is the exception: it reduces `log(count +
//! pseudo)`, which is not integer-valued, so it can disagree with R. It exists for
//! the CLI and the simulation harness, which have no R to ask. Anything whose
//! output must be byte-identical to the oracle uses [`Reductions`] and gets
//! [`RBackedReductions`] (in the bridge) instead.

use crate::workspace::RMatrix;

/// A sum or mean over one axis, with `na.rm = TRUE` semantics.
///
/// `na.rm` here means R's: an entry is skipped when it is a NaN of *any* payload,
/// because `is.na()` is true for `NA_real_` and for a computed `NaN` alike. A
/// `±Inf` is **not** skipped — R sums infinities, and only `NaN` is `ISNAN`.
pub trait Reductions {
    /// `rowMeans(x, na.rm = TRUE)`: one mean per row, over that row's non-NaN
    /// entries.
    ///
    /// A row with no non-NaN entry is `NA_real_`, not `NaN`: R assigns
    /// `R_NaReal` there, and the two have different bit patterns.
    fn row_means_na_rm(&self, m: &RMatrix) -> Vec<f64>;

    /// `colMeans(x, na.rm = TRUE)`.
    fn col_means_na_rm(&self, m: &RMatrix) -> Vec<f64>;

    /// `colSums(x, na.rm = TRUE)`.
    fn col_sums_na_rm(&self, m: &RMatrix) -> Vec<f64>;

    /// `rowSums(x, na.rm = TRUE)`.
    fn row_sums_na_rm(&self, m: &RMatrix) -> Vec<f64>;
}

/// R's `NA_real_`, as bits.
///
/// This is *not* `f64::NAN`. R's `NA_real_` is a NaN whose payload is
/// `0x7ff00000000007a2`, while `f64::NAN` is `0x7ff8000000000000`, and
/// `identical(NA_real_, NaN)` is `FALSE`. `serialize()` records the payload, so a
/// stage array that substitutes one for the other is not byte-identical.
///
/// Measured on the primary profile by writing the doubles out and reading them
/// back as `uint64`.
pub const NA_REAL_BITS: u64 = 0x7ff0_0000_0000_07a2;

/// `NA_real_` as an `f64`.
///
/// Not `const`: `f64::from_bits` only became usable in a `const` context in Rust
/// 1.83, and this crate's MSRV is 1.75. A `const` here would buy nothing.
#[inline]
pub fn na_real() -> f64 {
    f64::from_bits(NA_REAL_BITS)
}

/// Is this entry `NA` by R's rule?
///
/// `f64::is_nan` is exactly R's `ISNAN`, so this is `is.na(x)` — true for
/// `NA_real_` and for a computed `NaN`, false for `±Inf`. Note that
/// `is.infinite(NA_real_)` is `FALSE`, which is why the `-Inf → NA` mapping
/// below is keyed on `is_infinite` and not on "not finite".
#[inline]
pub fn is_na(v: f64) -> bool {
    v.is_nan()
}

/// Is this entry infinite? R's `is.infinite`: `±Inf` but not any NaN.
#[inline]
pub fn is_infinite(v: f64) -> bool {
    v.is_infinite()
}

/// The reductions used when there is no R to ask, and the reductions that are
/// *exactly* right for integer-valued inputs.
///
/// See the module docs for the one place this is an approximation rather than an
/// identity.
#[derive(Debug, Clone, Copy, Default)]
pub struct F64Reductions;

impl Reductions for F64Reductions {
    fn row_means_na_rm(&self, m: &RMatrix) -> Vec<f64> {
        let mut out = Vec::with_capacity(m.rows);
        for i in 0..m.rows {
            let mut s = 0.0f64;
            let mut n = 0usize;
            for &v in m.row(i) {
                if !v.is_nan() {
                    s += v;
                    n += 1;
                }
            }
            out.push(if n == 0 { na_real() } else { s / n as f64 });
        }
        out
    }

    fn col_means_na_rm(&self, m: &RMatrix) -> Vec<f64> {
        let mut sums = self.col_sums_na_rm(m);
        let counts = self.col_counts(m);
        for (s, &n) in sums.iter_mut().zip(counts.iter()) {
            *s = if n == 0 { na_real() } else { *s / n as f64 };
        }
        sums
    }

    fn col_sums_na_rm(&self, m: &RMatrix) -> Vec<f64> {
        (0..m.cols)
            .map(|j| {
                let mut s = 0.0f64;
                for i in 0..m.rows {
                    let v = m.get(i, j);
                    if !v.is_nan() {
                        s += v;
                    }
                }
                s
            })
            .collect()
    }

    fn row_sums_na_rm(&self, m: &RMatrix) -> Vec<f64> {
        (0..m.rows)
            .map(|i| {
                let mut s = 0.0f64;
                for &v in m.row(i) {
                    if !v.is_nan() {
                        s += v;
                    }
                }
                s
            })
            .collect()
    }
}

impl F64Reductions {
    /// How many non-NaN entries each column has: R's `rowSums(!is.na(x))` for the
    /// column axis, which is the denominator of a `colMeans`.
    pub fn col_counts(&self, m: &RMatrix) -> Vec<usize> {
        (0..m.cols)
            .map(|j| (0..m.rows).filter(|&i| !m.get(i, j).is_nan()).count())
            .collect()
    }

    /// How many non-NaN entries each row has.
    pub fn row_counts(&self, m: &RMatrix) -> Vec<usize> {
        (0..m.rows)
            .map(|i| m.row(i).iter().filter(|v| !v.is_nan()).count())
            .collect()
    }
}

/// Sum a slice of counts, skipping NaN: `sum(x, na.rm = TRUE)` for whole numbers.
///
/// R's `rowSums(!is.na(x))` and `rowSums(x != 0, na.rm = TRUE)` reduce 0/1 values,
/// and `colSums(feature_table, na.rm = TRUE)` reduces count integers. Every
/// partial sum is an exact integer, so the accumulator's width and the order of
/// the additions are both immaterial — which is why these three do not need a
/// `long double` and can stay plain `f64` here. The function asserts the
/// invariant it relies on rather than trusting the caller.
pub fn sum_counts(vals: impl Iterator<Item = f64>) -> f64 {
    let mut s = 0.0f64;
    for v in vals {
        if v.is_nan() {
            continue;
        }
        debug_assert!(
            v.fract() == 0.0 && v.abs() < 9.007_199_254_740_992e15,
            "a count reduction saw {v}, which is not a whole number below 2^53"
        );
        s += v;
    }
    s
}

/// `log(x)` with R's NaN-payload behaviour.
///
/// # Why this is not `x.ln()`
///
/// R's `log` is `.Internal(log(x))`, and R's own `R_log` (in `src/n/math/log.c`)
/// begins
///
/// ```c
/// if (ISNAN(x)) return x;
/// ```
///
/// so an `NA_real_` comes back as the *same* `NA_real_`, payload and all.
/// `f64::ln` is a call into the platform libm, and neither the C standard nor
/// IEEE 754 requires `log` to propagate a NaN payload. On the primary profile it
/// does not: `ln(0x7ff00000000007a2)` returns `0x7ff8000000000002` — a NaN with the
/// sign bit set and a different payload. `is.na()` is true of both, every
/// arithmetic result agrees, and `serialize()` disagrees.
///
/// Since the counts routinely contain `NA`, this is not a corner case: without the
/// NaN branch, every `NA` count acquires a different bit pattern on the way into
/// `y1`, and the whole-result byte comparison fails for a reason that has nothing to
/// do with the arithmetic.
///
/// `±0.0` and negatives are *not* special-cased here: `ln(0.0)` is `-Inf` in Rust
/// and in R alike, and the caller maps infinities to `NA_real_`.
#[inline]
pub fn r_log(x: f64) -> f64 {
    propagate_nan(x, f64::ln)
}

/// IEEE 754's *quiet* bit: bit 51 of a NaN's mantissa.
pub const QUIET_BIT: u64 = 0x0008_0000_0000_0000;

/// The IEEE sign bit, as a mask.
///
/// Not what R's arithmetic touches. It is here because the two are easy to confuse
/// when reading a hex dump, and confusing them cost an afternoon: `7ff0` becoming
/// `7ff8` looks like a sign change and is not one.
pub const SIGN_BIT: u64 = 0x8000_0000_0000_0000;

/// R's `NA_real_` as it comes back out of arithmetic: the payload with the quiet bit
/// set.
///
/// Named because it appears in every stage array and because "the quieted NA" is not a
/// thing anyone would guess. Measured: `NA_real_ - 1.5` is `7ff80000000007a2`.
pub const QUIETED_NA_REAL_BITS: u64 = NA_REAL_BITS | QUIET_BIT;

/// What R's *binary* arithmetic does to a NaN operand.
///
/// Measured on the primary profile:
///
/// ```text
/// NA_real_                  7ff00000000007a2   is.na TRUE, is.nan FALSE
/// NA_real_ - 1.5            7ff80000000007a2
/// 1.5 - NA_real_            7ff80000000007a2
/// NaN - 1.5                 7ff8000000000000
/// NA_real_ + NaN            7ff80000000007a2
/// NaN + NA_real_            7ff8000000000000
/// (NA_real_ - 1) - 2        7ff80000000007a2
/// NA_real_ * 0              7ff80000000007a2
/// NA_real_ %% 2             7ff80000000007a2
/// NA_real_ - Inf            7ff80000000007a2
/// -NA_real_                 fff00000000007a2
/// ```
///
/// Read `7ff0` -> `7ff8` carefully: bit 51 is the IEEE **quiet** bit, not the sign.
/// R's `NA_real_` is a *signalling* NaN -- quiet bit clear, payload `0x7a2` -- and any
/// arithmetic on it quiets it. The sign bit is never touched by arithmetic; the last
/// row above, R's unary `-`, is the one that sets it, and that is a sign flip rather
/// than arithmetic.
///
/// Two more things out of that table:
///
/// * the operand's **payload survives** -- `0x7a2` is still there, so the quieted NA
///   is distinguishable from a computed `NaN` by 18 bits;
/// * where **both** operands are NaN, the **left** one decides -- `NaN + NA_real_` is
///   the plain NaN and `NA_real_ + NaN` is R's.
///
/// So the rule is: take the deciding operand's bits and set the quiet bit, which is
/// `bits | QUIET_BIT`. It is idempotent, which is why `(NA_real_ - 1) - 2` is
/// unchanged by the second operation.
///
/// [`propagate_nan`] is *not* this rule, and using it for arithmetic is invisible to
/// every numeric check: `is.na()` is true either way, comparisons behave the same, and
/// arithmetic on either answer gives the same numbers. Only the bits differ, and only
/// `serialize()` looks at the bits.
#[inline]
pub fn r_nan_arith(x: f64) -> f64 {
    f64::from_bits(x.to_bits() | QUIET_BIT)
}

/// R's binary arithmetic, with the NaN rule above.
///
/// The deciding operand is the left one when both are NaN, which is what the
/// measurements show.
#[inline]
fn r_binop(x: f64, y: f64, f: impl Fn(f64, f64) -> f64) -> f64 {
    if x.is_nan() {
        r_nan_arith(x)
    } else if y.is_nan() {
        r_nan_arith(y)
    } else {
        f(x, y)
    }
}

/// Apply `f` to `x` unless `x` is a NaN, in which case return `x` itself.
///
/// This is the rule for R's *unary* math functions, and it is **not** the rule for
/// binary arithmetic. Measured, again on the primary profile:
///
/// ```text
/// log(NA_real_)   7ff00000000007a2   unchanged
/// sqrt(NA_real_)  7ff00000000007a2   unchanged
/// abs(NA_real_)   7ff00000000007a2   unchanged
/// log(NaN)        7ff8000000000000   unchanged
/// ```
///
/// R's `log` begins with `if (ISNAN(x)) return x;`, so payload and sign both survive.
/// `NA_real_ - 1.5` does not. One helper covering both would be wrong for one of them,
/// which is why there are two and the measurements sit next to each.
///
/// `abs(NA_real_)` being unchanged is consistent -- `abs` selects rather than
/// computes. R's unary `-` is arithmetic: it gives `fff00000000007a2`, a genuine sign
/// flip that neither rule here produces.
#[inline]
pub fn propagate_nan(x: f64, f: impl Fn(f64) -> f64) -> f64 {
    if x.is_nan() {
        x
    } else {
        f(x)
    }
}

/// R's `x - y`.
#[inline]
pub fn r_sub(x: f64, y: f64) -> f64 {
    r_binop(x, y, |a, b| a - b)
}

/// R's `x + y`.
#[inline]
pub fn r_add(x: f64, y: f64) -> f64 {
    r_binop(x, y, |a, b| a + b)
}

/// R's `x * y`.
#[inline]
pub fn r_mul(x: f64, y: f64) -> f64 {
    r_binop(x, y, |a, b| a * b)
}

/// R's `x / y`.
#[inline]
pub fn r_div(x: f64, y: f64) -> f64 {
    r_binop(x, y, |a, b| a / b)
}

/// Count entries that are non-zero and non-NaN: R's
/// `rowSums(x != 0, na.rm = TRUE)`.
///
/// `x != 0` is `NA` for an `NA` entry, and `na.rm = TRUE` drops it — so an `NA` is
/// neither a presence nor an absence. It is *not* counted in the numerator, and
/// it is not counted in `rowSums(!is.na(x))` either, so it leaves the ratio
/// entirely.
pub fn count_nonzero_observed(vals: impl Iterator<Item = f64>) -> usize {
    vals.filter(|&v| !v.is_nan() && v != 0.0).count()
}

/// Count non-NaN entries: R's `rowSums(!is.na(x))`.
pub fn count_observed(vals: impl Iterator<Item = f64>) -> usize {
    vals.filter(|&v| !v.is_nan()).count()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn m(rows: usize, cols: usize, data: &[f64]) -> RMatrix {
        RMatrix::from_row_major(rows, cols, data.to_vec())
    }

    /// R's `NA_real_` is `0x7ff00000000007a2`, and `f64::NAN` is
    /// `0x7ff8000000000000`. Every `is_nan()` treats them alike; `serialize()`
    /// does not.
    #[test]
    fn na_real_is_not_f64_nan() {
        assert_eq!(na_real().to_bits(), 0x7ff0_0000_0000_07a2);
        assert_ne!(na_real().to_bits(), f64::NAN.to_bits());
        assert!(is_na(na_real()));
        assert!(is_na(f64::NAN));
    }

    /// R's `is.infinite` is true for `±Inf` and false for *every* NaN, including
    /// `NA_real_`. This is the distinction `o[is.infinite(o)] <- NA` turns on:
    /// it replaces `-Inf` from `log(0)` and leaves an existing `NA_real_` alone.
    #[test]
    fn is_infinite_separates_infinity_from_every_nan() {
        assert!(is_infinite(f64::INFINITY));
        assert!(is_infinite(f64::NEG_INFINITY));
        assert!(!is_infinite(na_real()));
        assert!(!is_infinite(f64::NAN));
        assert!(!is_infinite(0.0));
        assert!(!is_infinite(-1.5));
    }

    /// `log(0)` is `-Inf`; `x[is.infinite(x)] <- NA` turns it into `NA_real_`,
    /// whose bits are `0x7ff00000000007a2`. A stage array holding a bare
    /// `f64::NAN` there is *not* byte-identical to the oracle's.
    #[test]
    fn the_log_zero_replacement_lands_on_rs_na_payload() {
        let got = if is_infinite(-1.0 / 0.0) {
            na_real()
        } else {
            f64::NAN
        };
        assert_eq!(got.to_bits(), 0x7ff0_0000_0000_07a2);
        assert_ne!(got.to_bits(), f64::NAN.to_bits());
    }

    /// `r_log` returns an `NA_real_` unchanged, payload and all, which `f64::ln`
    /// does not do.
    ///
    /// The measured libm answer on the primary profile is `0x7ff8000000000002`: a
    /// NaN, so `is.na()` is true and nothing numeric notices, but the payload
    /// differs from both the input and R's.
    #[test]
    fn r_log_preserves_the_na_payload_where_f64_ln_does_not() {
        assert_eq!(r_log(na_real()).to_bits(), NA_REAL_BITS);
        assert_eq!(r_log(f64::NAN).to_bits(), f64::NAN.to_bits());
        // And the measured divergence, so this test fails if a platform's libm
        // starts propagating and the branch looks redundant.
        assert_ne!(
            na_real().ln().to_bits(),
            NA_REAL_BITS,
            "this libm now preserves the payload, so the branch is dead code here"
        );
        // Finite values are untouched.
        assert_eq!(r_log(1.0), 0.0);
        assert_eq!(r_log(4.0), 4.0f64.ln());
        assert_eq!(r_log(0.0), f64::NEG_INFINITY);
        assert!(r_log(-1.0).is_nan());
    }

    /// R's binary arithmetic sets the sign bit on a NaN operand and keeps its payload.
    ///
    /// Every row here is a measurement from the primary profile; they are in this test
    /// rather than only in the module docs because a doc comment is not run.
    ///
    /// This is the centring step -- `y1 = o1 - rowMeans(o1, na.rm = TRUE)` -- where an
    /// `NA` count passes through, so the payload ends up in the response tensor the
    /// whole estimator consumes.
    #[test]
    fn binary_arithmetic_quiets_the_nan_and_keeps_the_payload() {
        let na = na_real();
        let nan = f64::NAN;
        let quiet_na = f64::from_bits(na.to_bits() | QUIET_BIT);
        // A computed `NaN` is already quiet, so its answer is itself.
        let quiet_nan = nan;

        // Operand order and operator do not change the answer.
        for (label, got) in [
            ("na - 1.5", r_sub(na, 1.5)),
            ("1.5 - na", r_sub(1.5, na)),
            ("na + 1.5", r_add(na, 1.5)),
            ("1.5 + na", r_add(1.5, na)),
            ("na * 2", r_mul(na, 2.0)),
            ("2 * na", r_mul(2.0, na)),
            ("na / 2", r_div(na, 2.0)),
            ("2 / na", r_div(2.0, na)),
            ("na * 0", r_mul(na, 0.0)),
            ("na / 0", r_div(na, 0.0)),
            ("0 / na", r_div(0.0, na)),
            ("na - Inf", r_sub(na, f64::INFINITY)),
            ("Inf - na", r_sub(f64::INFINITY, na)),
        ] {
            assert_eq!(got.to_bits(), quiet_na.to_bits(), "{label}: {got:e}");
        }

        // A plain NaN is already quiet, so it comes back as itself.
        assert_eq!(r_sub(nan, 1.5).to_bits(), quiet_nan.to_bits());
        assert_eq!(r_add(nan, 1.5).to_bits(), quiet_nan.to_bits());

        // Both NaN: the *left* operand decides.
        assert_eq!(r_add(nan, na).to_bits(), quiet_nan.to_bits());
        assert_eq!(r_add(na, nan).to_bits(), quiet_na.to_bits());

        // Idempotent: a second operation changes nothing.
        assert_eq!(r_sub(r_sub(na, 1.0), 2.0).to_bits(), quiet_na.to_bits());

        // The quiet bit is bit 51, not the sign. Asserting the sign separately is the
        // whole point: reading `7ff8` as a sign change is the mistake this test exists
        // to prevent.
        assert_eq!(
            quiet_na.to_bits() & SIGN_BIT,
            0,
            "arithmetic must not touch the sign"
        );
        assert_eq!(
            na.to_bits() & QUIET_BIT,
            0,
            "NA_real_ is a *signalling* NaN"
        );
        assert_ne!(
            quiet_na.to_bits() & QUIET_BIT,
            0,
            "and arithmetic quiets it"
        );

        // The hardware's own answer, and the reason this function exists: x86 returns
        // the first NaN operand with the **payload cleared**, so `NA_real_ - 1.5` is
        // `7ff8000000000000` where R is `7ff80000000007a2`. Eighteen bits of payload
        // lost, no numeric consequence, and a `serialize()` difference.
        assert_eq!((na - 1.5).to_bits(), f64::NAN.to_bits());
        assert_ne!((na - 1.5).to_bits(), quiet_na.to_bits());

        // Finite arithmetic is untouched.
        assert_eq!(r_sub(3.0, 1.5), 1.5);
        assert_eq!(r_add(3.0, 1.5), 4.5);
        assert_eq!(r_mul(3.0, 1.5), 4.5);
        assert_eq!(r_div(3.0, 1.5), 2.0);
    }

    /// A unary math function returns its NaN operand *unchanged*, which is a different
    /// rule from binary arithmetic and was the other half of the same bug.
    ///
    /// `log(NA_real_)` is `7ff00000000007a2` -- quiet bit still clear -- while
    /// `NA_real_ - 1.5` is `7ff80000000007a2`. One helper cannot do both.
    #[test]
    fn unary_math_functions_return_their_nan_operand_unchanged() {
        let na = na_real();
        assert_eq!(r_log(na).to_bits(), na.to_bits());
        assert_eq!(r_log(na).to_bits(), NA_REAL_BITS, "not the signed form");
        assert_eq!(r_log(f64::NAN).to_bits(), f64::NAN.to_bits());
        // And the two rules really do disagree.
        assert_ne!(r_log(na).to_bits(), r_sub(na, 1.5).to_bits());
    }

    /// A computed `NaN` — `Inf - Inf`, say — is not infinite, so the reference's
    /// `is.infinite` mask does not touch it and it keeps the plain payload. A
    /// transport that mapped "not finite" to `NA_real_` would turn a computed NaN
    /// into `NA_real_`, which is the S06 mistake in a different place.
    #[test]
    fn a_computed_nan_is_not_rewritten_to_na_real() {
        let computed = f64::INFINITY - f64::INFINITY;
        assert!(computed.is_nan());
        assert!(!is_infinite(computed));
        assert_eq!(computed.to_bits(), f64::NAN.to_bits());
    }

    /// A row with no observed entry is `NA_real_` in R's `rowMeans`, not `NaN`.
    #[test]
    fn an_all_missing_row_is_na_real() {
        let mm = m(2, 2, &[f64::NAN, f64::NAN, 1.0, 3.0]);
        let means = F64Reductions.row_means_na_rm(&mm);
        assert_eq!(means[0].to_bits(), 0x7ff0_0000_0000_07a2);
        assert_eq!(means[1], 2.0);
    }

    /// `na.rm = TRUE` skips NaN but sums `±Inf`, because R's mask is `ISNAN`.
    #[test]
    fn na_rm_skips_nan_and_keeps_infinity() {
        let mm = m(1, 3, &[f64::NAN, f64::INFINITY, 1.0]);
        assert_eq!(F64Reductions.row_sums_na_rm(&mm), vec![f64::INFINITY]);
        let mm2 = m(1, 2, &[f64::NAN, f64::NEG_INFINITY]);
        assert_eq!(F64Reductions.row_sums_na_rm(&mm2), vec![f64::NEG_INFINITY]);
    }

    /// The integer-valued reductions are exact, so they are not approximations.
    ///
    /// This is the property that lets `preprocess.rs` keep using plain `f64` for
    /// prevalence and library sizes while the log-scale means go through the
    /// adapter.
    #[test]
    fn count_reductions_are_exact_for_whole_numbers() {
        let vals: Vec<f64> = (0..1000)
            .map(|i| if i % 7 == 0 { 0.0 } else { 13.0 })
            .collect();
        let got = sum_counts(vals.iter().copied());
        assert_eq!(got, (1000 - (1000 + 6) / 7) as f64 * 13.0);
        assert_eq!(
            count_nonzero_observed(vals.iter().copied()),
            (1000 - (1000 + 6) / 7) as usize
        );
        // A sum of 0/1 across a thousand columns is still exactly representable.
        let flags: Vec<f64> = (0..1000).map(|i| (i % 3 == 0) as u8 as f64).collect();
        assert_eq!(sum_counts(flags.iter().copied()), ((999) / 3 + 1) as f64);
    }

    /// An `NA` is neither a presence nor an absence, so it leaves prevalence
    /// entirely rather than counting in either term.
    ///
    /// The two plausible wrong rules give 0 and 0.2 here, and the right one gives
    /// 0.25 — so the assertion distinguishes all three rather than just being true.
    #[test]
    fn an_na_leaves_prevalence_rather_than_counting_as_absent() {
        let row = [f64::NAN, 0.0, 0.0, 5.0, 0.0];
        let observed = count_observed(row.iter().copied());
        let nonzero = count_nonzero_observed(row.iter().copied());
        // 4 observed (three zeros and the 5), one of them non-zero.
        assert_eq!(observed, 4);
        assert_eq!(nonzero, 1);
        assert_eq!(nonzero as f64 / observed as f64, 0.25);
        // Counting the NA as an absence in both terms would give 1/5 = 0.2.
        assert_ne!(1.0 / 5.0, nonzero as f64 / observed as f64);
        // Dropping it from the numerator but not the denominator gives 0.
        assert_ne!(0.0, nonzero as f64 / observed as f64);
    }

    /// The case that motivates the whole module, with R's own answers.
    ///
    /// `rowMeans(c(1e16, 1, -1e16, 1e-17, 1))` is `0x1.999999999999ap-2` on the
    /// primary profile: R's `long double` kept the `+1` that an `f64` sum loses,
    /// so it sums to 2 where `f64` sums to 1. This test cannot assert R's value —
    /// it has no R — so it asserts the *direction* and magnitude of the gap, and
    /// the R-side comparator in `scripts/check_preprocess_stages.R` asserts R's
    /// exact bits.
    #[test]
    fn a_cancelling_sum_is_where_f64_and_r_disagree() {
        let y = [1e16, 1.0, -1e16, 1e-17, 1.0];
        let f64_sum: f64 = y.iter().sum();
        assert_eq!(f64_sum, 1.0, "the f64 sum loses the +1 as expected");
        // The exact sum is 2 + 1e-17, so an accumulator with more than 53 bits of
        // mantissa reports something near 0.4 rather than 0.2. `F64Reductions`
        // reports 0.2, and that is the bug this records.
        let got = F64Reductions.row_means_na_rm(&m(1, 5, &y))[0];
        assert_eq!(got, 0.2);
        assert!(
            got < 0.3,
            "this test is worthless unless f64 really loses it"
        );
    }
}

#[cfg(test)]
mod split_tests {
    //! The sensitivity refits parallelise, and the reduction may be R, which cannot be
    //! called from a worker thread. So the means are reduced on the calling thread and
    //! the workers redo the arithmetic. That is only safe if the two routes agree cell
    //! for cell, which is what these check.
    //!
    //! It is a real risk rather than a formality: the split gives the workers a `ln`
    //! per cell instead of a whole log table, and the centring then subtracts a mean it
    //! did not compute. A drift between the two would put the sensitivity scores on a
    //! different scale from the main run's q-values, which are compared directly.

    use super::*;
    use crate::preprocess::CountMatrix;
    use crate::workspace::RMatrix;

    fn table() -> CountMatrix {
        CountMatrix::new(
            4,
            5,
            vec![
                0.0, 3.0, 0.0, 7.0, 1.0, //
                4.0, 4.0, 0.0, 0.0, 0.0, //
                9.0, 1.0, 2.0, 0.0, 6.0, //
                0.0, 0.0, 0.0, 5.0, 0.0,
            ],
        )
        .unwrap()
    }

    /// A reduction that is deliberately *not* `f64`, so a mix-up cannot pass.
    ///
    /// Adding a fixed epsilon to every mean is not anything R does; it is here so the
    /// test can tell "the split carried the right numbers" from "both routes happened
    /// to agree because the difference is invisible".
    struct Shifted;

    impl Reductions for Shifted {
        fn row_means_na_rm(&self, m: &RMatrix) -> Vec<f64> {
            let mut v = F64Reductions.row_means_na_rm(m);
            for x in v.iter_mut() {
                *x += 0.125;
            }
            v
        }
        fn col_means_na_rm(&self, m: &RMatrix) -> Vec<f64> {
            F64Reductions.col_means_na_rm(m)
        }
        fn col_sums_na_rm(&self, m: &RMatrix) -> Vec<f64> {
            F64Reductions.col_sums_na_rm(m)
        }
        fn row_sums_na_rm(&self, m: &RMatrix) -> Vec<f64> {
            F64Reductions.row_sums_na_rm(m)
        }
    }

    /// Bit-for-bit, because `==` is the wrong comparison here: a `NaN` cell is never
    /// equal to itself, so `assert_eq!` on these vectors fails on the cells that matter
    /// most -- the `NA`s. And comparing bits is the point: a route that differed only
    /// in a NaN payload would pass `==` on nothing and fail here, which is the whole
    /// reason this crate has a NaN rule at all.
    fn same_bits(a: &[f64], b: &[f64]) -> bool {
        a.len() == b.len()
            && a.iter()
                .zip(b.iter())
                .all(|(x, y)| x.to_bits() == y.to_bits())
    }

    #[test]
    fn log_center_rows_with_equals_log_then_centre_by_supplied_means() {
        let c = table();
        for pseudo in [0.0, 0.5, 1.0, 7.25] {
            for rows in [vec![0usize, 1, 2, 3], vec![1usize, 3], vec![2usize]] {
                let whole = c.log_center_rows_with(&rows, pseudo, &F64Reductions);
                let means = c.log_row_means(&rows, pseudo, &F64Reductions);
                let mut split = c.log_rows(&rows, pseudo);
                split.sub_rows_in_place(&means);
                assert!(
                    same_bits(&whole.data, &split.data),
                    "rows {rows:?} at pseudo {pseudo}: the split differs"
                );
            }
        }
    }

    /// The same, with a reduction that is *not* `f64`, so the test can tell a correct
    /// split from a coincidence.
    #[test]
    fn the_split_carries_whichever_reduction_was_used() {
        let c = table();
        let rows = vec![0usize, 2, 3];
        let whole = c.log_center_rows_with(&rows, 0.5, &Shifted);
        let means = c.log_row_means(&rows, 0.5, &Shifted);
        let mut split = c.log_rows(&rows, 0.5);
        split.sub_rows_in_place(&means);
        assert!(same_bits(&whole.data, &split.data));
        // And it really is different from the `f64` answer, so the test is not vacuous.
        let plain = c.log_center_rows_with(&rows, 0.5, &F64Reductions);
        assert!(!same_bits(&whole.data, &plain.data));
    }

    /// The zero-replacing variant, which the non-conservative sensitivity path uses.
    #[test]
    fn the_zero_replacing_split_agrees_too() {
        let c = table();
        let rows = vec![0usize, 2];
        let cols = vec![1usize, 3, 4];
        for pseudo in [0.0, 0.25, 2.0] {
            let whole = c.log_center_replacing_zeros_sub_with(&rows, &cols, pseudo, &Shifted);
            let log = c.log_replacing_zeros_sub(&rows, &cols, pseudo);
            let means = Shifted.row_means_na_rm(&log);
            let mut split = log;
            split.sub_rows_in_place(&means);
            assert!(same_bits(&whole.data, &split.data), "pseudo {pseudo}");
        }
    }
}
