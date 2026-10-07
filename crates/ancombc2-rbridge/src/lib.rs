//! `ancombc2-rbridge`: the typed native bridge between R and the numerical core.
//!
//! # Why this crate exists
//!
//! The pre-existing `ancombc2-ffi` bridge speaks JSON. That is exactly the design
//! `IMPROVED_PLAN.md` section 2 names as a required closure: a JSON round trip
//! cannot carry `NA_real_`, `NaN`, `+Inf`, `-Inf` and `-0.0` as five distinguishable
//! things, because JSON has one "not a number" and one zero; and rebuilding the R
//! object from parsed JSON changes its structure, its attributes and its types.
//! This crate replaces it on the exact path and carries typed buffers instead.
//!
//! `ancombc2-ffi` is left in place, unmodified and unused by the replacement, so
//! the historical bridge remains buildable and the change stays reviewable.
//!
//! # Layering
//!
//! ```text
//!   r/ANCOMBC/R/*.R          public API, argument handling, result assembly
//!         |  .Call, typed
//!   r/ANCOMBC/src/init.c     R types -> borrowed pointers; allocates every buffer
//!         |  extern "C", pointers + lengths
//!   ancombc2-rbridge         validation, owned copies, the numerical seam
//!         |
//!   ancombc2-core            numerics; knows nothing about R
//! ```
//!
//! Three properties hold by construction rather than by review:
//!
//! * `ancombc2-core` has no R dependency, so `cargo build -p ancombc2-core`
//!   succeeds on a machine with no R installed. This crate is the only place the
//!   two meet.
//! * No R pointer reaches a parallel job: [`owned::StageInput`] is built from
//!   borrowed pointers, copied, and the borrow ends before anything is handed to
//!   Rayon. See [`owned`].
//! * No unwinding crosses the FFI boundary: every entry point is wrapped in
//!   `catch_unwind` and a panic becomes a status code plus a message.
//!
//! # No serialisation
//!
//! There is no `serde` and no `serde_json` dependency, and there is no string
//! format for a numeric value anywhere in the crate. The only bytes that cross the
//! boundary are packed UTF-8 *names* (see [`transport::PackedNames`]), because a
//! name is text and text is what bytes are for. Every number crosses as a number,
//! in its own type, with its bits intact.

#![deny(missing_docs)]
#![allow(clippy::needless_range_loop)]

pub mod abi;
pub mod output;
pub mod owned;
pub mod transport;

pub use abi::{
    ancombc2_rb_copy_last_error, ancombc2_rb_last_error_len, ancombc2_rb_oracle_sha,
    ancombc2_rb_probe, ancombc2_rb_version, Echo, RawRequest, KIND_INT, KIND_NONE, KIND_RAW,
    KIND_REAL, RB_ERR_ARGS, RB_ERR_PANIC, RB_ERR_STAGE, RB_ERR_STAGE_UNKNOWN, RB_OK,
};
pub use output::{
    core_output_plan, emit_core_payloads, preprocess_stages, ArrayDesc, CorePayloads, Emitted,
    PreprocessStages, RBackedReductions, RType, Reducer,
};
pub use owned::{Owned, StageInput, TransportValue, NA_REAL_BITS};
pub use transport::{BridgeError, Controls, Counts, Flags, Layout, PackedNames, Request, Result};

/// The pinned original this bridge reproduces.
pub const ORACLE_SHA: &str = "dc4febdf59badb3a8dfe0c767ef2186323c2199a";
