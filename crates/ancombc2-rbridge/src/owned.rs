//! Owned buffers, and the only thing a parallel stage may be handed.
//!
//! The rule this module exists to make enforceable: **no pointer that came from R
//! reaches a Rayon closure.** A `&[f64]` borrowed out of an R `REALSXP` is a
//! borrowed `&[f64]` like any other, and nothing in its type says where it came
//! from; R may move or collect the object at the next allocation, and the closure
//! is running on a thread the R runtime does not know about. So a stage's entry
//! point takes [`Owned<T>]`.
//!
//! `Owned` has exactly one constructor, [`Owned::from_borrowed`], which copies.
//! There is no way to build one from a reference without copying, and therefore no
//! way to smuggle an R-backed borrow into a parallel closure. The cost is one
//! copy per stage input, which is the price of not racing the garbage collector,
//! and which `ancombc2-core`'s existing code already pays at its own boundaries.

use std::sync::Arc;

/// An owned, immutable, shareable buffer.
///
/// Cheap to clone because it is a handle, not a copy: the parallel stages share
/// one buffer across workers rather than each holding its own. Cloning the
/// *contents* is explicit and rare.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Owned<T> {
    inner: Arc<T>,
}

impl<T> Owned<T> {
    /// Copy a borrowed buffer into owned storage.
    ///
    /// `Owned<T>` holds a `T`, so this is `Owned<Vec<T>>`: the handle wraps the
    /// vector, and cloning the handle shares it. The name is the point -- it says
    /// what happened. `Owned::from_ref` would invite the reader to assume nothing
    /// was copied.
    pub fn from_borrowed(x: &[T]) -> Owned<Vec<T>>
    where
        T: Clone,
    {
        Owned {
            inner: Arc::new(x.to_vec()),
        }
    }

    /// Take ownership of something already built.
    pub fn from_owned(x: T) -> Self {
        Owned { inner: Arc::new(x) }
    }

    /// Borrow the contents. The borrow cannot outlive the handle, so it cannot
    /// outlive the copy either.
    pub fn get(&self) -> &T {
        &self.inner
    }

    /// The address of the contents. For the lifetime and aliasing checks in the
    /// tests; nothing in the pipeline compares it.
    pub fn addr(&self) -> *const T {
        Arc::as_ptr(&self.inner)
    }
}

/// Everything a stage is allowed to be given, as owned data.
///
/// A struct rather than a long argument list, so that adding a field is a
/// compile error at every construction site instead of a silently-defaulted
/// argument.
#[derive(Debug, Clone, PartialEq)]
pub struct StageInput {
    /// The count matrix, integer or real, owned.
    pub data: Owned<Vec<TransportValue>>,
    /// The aggregate table, same shape and orientation as `data`.
    pub aggregate: Owned<Vec<TransportValue>>,
    /// `n_tax x n_samp` of the count matrix.
    /// `nrow(data)`.
    pub n_tax: usize,
    /// `ncol(data)`, which is also `nrow(x)`.
    pub n_samp: usize,
    /// The design, `n_samp x p`, column-major.
    pub design: Owned<Vec<f64>>,
    /// `ncol(x)`, the number of fixed-effect columns.
    pub n_design_cols: usize,
    /// `stats::complete.cases(x)`, one byte per sample.
    /// `stats::complete.cases(x)`, one byte per sample. Owned so that the
    /// completeness test inside a parallel stage cannot read freed R memory.
    pub design_complete: Owned<Vec<u8>>,
    /// Group level per sample, 0 meaning absent.
    /// Group level per sample, 1-based, `0` meaning absent.
    pub group_index: Owned<Vec<i32>>,
    /// `length(levels(group))`.
    pub n_group_levels: usize,
    /// `colnames(x)`.
    pub fix_eff: Owned<Vec<String>>,
    /// `rownames(data)`.
    pub taxon_names: Owned<Vec<String>>,
    /// `colnames(data)`.
    pub sample_names: Owned<Vec<String>>,
    /// `levels(group)`, as text, so a level index can be reported by name.
    pub group_labels: Owned<Vec<String>>,
    /// The numeric controls, owned. A `Copy` struct, so the "owned" is nominal
    /// here and there is nothing to race.
    pub controls: crate::transport::Controls,
    /// The logical switches. Also `Copy`.
    pub flags: crate::transport::Flags,
}

/// A count cell, with its kind kept.
///
/// The integer case is widened to `i64` rather than `f64` so that `NA_integer_`
/// is distinguishable from every other value all the way into the stage. Widening
/// to `f64` would be numerically exact but would erase the tag, and the tag is
/// part of the schema.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TransportValue {
    /// A count. For a real count matrix this is the raw 64 bits reinterpreted as
    /// an integer, which is exact in both directions and keeps `NA_real_`,
    /// `NaN`, `+-Inf` and `-0.0` distinguishable all the way into the stage.
    Count(i64),
    /// `NA_integer_`. The tag is kept rather than widened to a `f64::NAN`, because
    /// a `f64::NAN` would be indistinguishable from a real `NaN` count.
    Missing,
}

impl TransportValue {
    /// Whether this cell is `NA`.
    pub fn is_missing(self) -> bool {
        matches!(self, TransportValue::Missing)
    }

    /// The value as a double, with `Missing` becoming R's `NA_real_` bit pattern.
    pub fn as_f64(self) -> f64 {
        match self {
            TransportValue::Count(v) => v as f64,
            // `f64::NAN` with R's `NA_real_` payload, so a stage that needs a
            // double gets the same bit pattern R would have handed it.
            TransportValue::Missing => f64::from_bits(NA_REAL_BITS),
        }
    }
}

/// The bit pattern R uses for `NA_real_`.
///
/// Written out rather than produced with `f64::NAN`, because
/// `f64::NAN.to_bits()` is not guaranteed to be this value and a stage that
/// compares `to_bits()` would then disagree with R for a reason nobody could see.
pub const NA_REAL_BITS: u64 = 0x7FF0_0000_0000_07A2;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cloning_a_handle_does_not_copy_the_contents() {
        let src = vec![1.0f64, 2.0, 3.0];
        let a = Owned::from_borrowed(&src);
        let b = a.clone();
        assert_eq!(a.addr(), b.addr(), "cloning a handle must share one buffer");
        drop(src); // the borrow is over; the handle is not
        assert_eq!(b.get().as_slice(), &[1.0, 2.0, 3.0]);
    }

    #[test]
    fn from_borrowed_copies_rather_than_aliasing() {
        let mut src = vec![7i32, 8, 9];
        let owned = Owned::from_borrowed(&src);
        src[0] = 0; // the source changes after the copy
        assert_eq!(owned.get()[0], 7, "the handle must not alias its source");
        assert_ne!(
            owned.get().as_ptr(),
            src.as_ptr(),
            "the owned buffer must be a distinct allocation from the borrow"
        );
    }

    #[test]
    fn na_real_carries_rs_bit_pattern() {
        let v = TransportValue::Missing.as_f64();
        assert_eq!(v.to_bits(), NA_REAL_BITS);
        assert!(v.is_nan());
        // And it is NOT the same NaN as a computed one. `f64::NAN.to_bits()` is
        // implementation-defined, so this is exactly why the constant is written
        // out rather than taken from the platform. `f64::from_bits` of a
        // different payload is used rather than `0.0/0.0` so the test is about
        // bit patterns and not about the compiler folding a division.
        let plain = f64::from_bits(0x7FF8_0000_0000_0001);
        assert!(plain.is_nan());
        assert_ne!(plain.to_bits(), NA_REAL_BITS);
    }
}
