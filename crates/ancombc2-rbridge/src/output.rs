//! The output transport: what `.ancombc2_core` produces, and how it crosses back.
//!
//! # What crosses, and what does not
//!
//! The numerical work is not in this module. What is here is the inventory of the
//! arrays `.ancombc2_core` returns *internally*, the shapes of each, and a typed
//! copy of each into buffers R allocated. Assembly into `data.frame`s happens on
//! the R main thread, with the original's own expressions.
//!
//! # The shapes are a function of the shapes
//!
//! Every internal array's dimensions follow from `(n_tax, n_samp, p)`. Nothing
//! else determines them. So the plan can be computed on both sides from those three
//! numbers and the two sides cannot disagree about a length -- which is the whole
//! reason there is no per-array length on the wire.
//!
//! # Types, measured not assumed
//!
//! `docs/exact_output_inventory.md` records where each fact below came from. The
//! three that would have been wrong if guessed:
//!
//! * **`dof` is an integer vector when every `df.residual` is whole.** R keeps an
//!   integer vector's type when a whole-valued double is assigned into it, and
//!   `.lm_fit_all`'s `dof = rep(999L, n_tax)` is integer. Captured and confirmed:
//!   `typeof(dof) == "integer"` on the payloads. Widening it to double would be
//!   numerically inert and would change the type of a returned field.
//! * **`O2` is always double.** `.ancombc2_core` computes `O2 = aggregate_data +
//!   pseudo`, and adding a double to an integer matrix promotes it. So
//!   `feature_table` is double even when the input counts are integer.
//! * **`vcov_hat` is a *list* of `n_tax` matrices, not one `n_tax*p*p` matrix.**
//!   `.sandwich_vcov` builds it with `lapply` and gives each element `dimnames`.
//!   Flattening it into one array would lose the per-element `dimnames` that
//!   `.ancombc_global_F` and `.ancombc_trend` read via `x[group_ind, group_ind]`.
//!
//! # The bias and reported taxon sets
//!
//! `taxa_bias` (the taxa that estimated the sample-specific biases) is a superset
//! of `taxa` (the taxa in the reported table) whenever `struc_zero` flags
//! anything: `.data_core` is called twice, once with `tax_keep = NULL` for `O1`
//! and once with the structural zeros removed for `O2`. They are separate fields
//! here for exactly that reason, and merging them would silently change what
//! `theta` was estimated from.

use crate::transport::{BridgeError, Layout, PackedNames, Result};

/// The R storage type of an emitted array.
///
/// Written out rather than inferred from a Rust type, because the point of the
/// output transport is that R gets the type the original would have given it --
/// and `typeof()` is what `serialize()` records.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RType {
    /// `INTSXP`.
    Int,
    /// `REALSXP`.
    Real,
    /// `LGLSXP`. R's logical is an `int`, and the original's `diff_abn` is logical.
    Logical,
    /// `STRSXP`.
    Character,
    /// `RAWSXP`.
    Raw,
    /// A `VECSXP` of `count` `REALSXP` matrices of `rows x cols`.
    ListOfRealMatrix,
}

impl RType {
    /// The name reported to R, for a diagnostic that names a type.
    pub fn name(self) -> &'static str {
        match self {
            RType::Int => "integer",
            RType::Real => "double",
            RType::Logical => "logical",
            RType::Character => "character",
            RType::Raw => "raw",
            RType::ListOfRealMatrix => "list of matrices",
        }
    }
}

/// One emitted array.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ArrayDesc {
    /// The name the original gives it.
    pub name: &'static str,
    /// The R storage type.
    pub rtype: RType,
    /// `nrow`, or the length for a vector.
    pub rows: usize,
    /// `ncol`, or 0 for a vector.
    pub cols: usize,
    /// Element order on the wire.
    ///
    /// [`Layout::ColumnMajor`] for everything except `vcov_hat`, which is
    /// [`Layout::RowMajor`] within each taxon. It is a per-entry field rather
    /// than a module-wide fact because the covariance is the one array whose
    /// order is *not* R's, and a module-wide claim would have hidden exactly
    /// that. `vcov_hat_is_the_only_row_major_entry` pins it.
    pub layout: Layout,
    /// How many elements: `rows * cols`, or the length of the whole object for a
    /// list.
    pub count: usize,
}

impl ArrayDesc {
    fn vec3(name: &'static str, rtype: RType, rows: usize, cols: usize) -> Self {
        ArrayDesc {
            name,
            rtype,
            rows,
            cols,
            layout: Layout::ColumnMajor,
            count: rows * cols,
        }
    }

    fn vector(name: &'static str, rtype: RType, len: usize) -> Self {
        ArrayDesc {
            name,
            rtype,
            rows: len,
            cols: 0,
            layout: Layout::ColumnMajor,
            count: len,
        }
    }
}

/// Every internal array `.ancombc2_core` produces, in the order the original
/// produces them.
///
/// The order is the original's, not a convenient one: `beta_hat` before `var_hat`
/// because `var_hat` is derived from it, `vcov_hat` after `var_hat` because the
/// original overwrites its diagonal with it, `s02` last because it is added to
/// `var_hat` after the covariance is formed. Reading the plan top to bottom is
/// reading `ancombc_prep.R:216-235` top to bottom.
pub fn core_output_plan(
    n_tax: usize,
    n_samp: usize,
    p: usize,
    dof_is_int: bool,
) -> Result<Vec<ArrayDesc>> {
    if n_tax == 0 || n_samp == 0 || p == 0 {
        return Err(BridgeError::ShapeMismatch {
            field: "core outputs",
            rows: n_tax,
            cols: n_samp,
            len: n_tax * n_samp,
        });
    }
    Ok(vec![
        // prep.R:236-241 -- vcov_hat has var_hat written into its diagonal first.
        ArrayDesc::vector("theta_hat", RType::Real, n_samp),
        ArrayDesc::vec3("beta_hat", RType::Real, n_tax, p),
        ArrayDesc::vec3("var_hat", RType::Real, n_tax, p),
        ArrayDesc::vec3(
            "dof",
            if dof_is_int { RType::Int } else { RType::Real },
            n_tax,
            p,
        ),
        ArrayDesc {
            name: "vcov_hat",
            rtype: RType::ListOfRealMatrix,
            rows: n_tax,
            cols: p,
            // Row-major, unlike every other entry here. The covariance block is
            // built by `apply()`-style code upstream, and the transpose to R's
            // order belongs in exactly one place -- the list rebuild in
            // `C_ancombc2_rb_emit` -- rather than being smeared across both
            // ends of the wire.
            layout: Layout::RowMajor,
            count: n_tax * p * p,
        },
        ArrayDesc::vector("delta_em", RType::Real, p),
        ArrayDesc::vector("delta_wls", RType::Real, p),
        ArrayDesc::vector("var_delta", RType::Real, p),
        ArrayDesc::vector("s02", RType::Real, p),
        // prep.R:199-201 -- y_bias_crt is a data.frame, transported as its
        // columns plus its dimnames.
        ArrayDesc::vec3("y_bias_crt", RType::Real, n_tax, n_samp),
        // prep.R:112 -- O2 = aggregate_data + pseudo, always double.
        ArrayDesc::vec3("O2", RType::Real, n_tax, n_samp),
    ])
}

/// The payloads `.ancombc2_core` hands to its own assembly, once validated.
#[derive(Debug, Clone, PartialEq)]
pub struct CorePayloads {
    /// `nrow(O2)`, the reported taxa.
    pub n_tax: usize,
    /// `ncol(O2)`.
    pub n_samp: usize,
    /// `ncol(x)`, the fixed-effect columns.
    pub p: usize,
    /// `n_tax x p`.
    pub beta_hat: Vec<f64>,
    /// `n_tax x p`, after the `var_delta` and `s02` adjustments and after
    /// `var_hat[is.na(beta_hat)] = NA`. This is the `var_hat` the original's
    /// pairwise and Dunnett paths read, not the sandwich diagonal.
    pub var_hat: Vec<f64>,
    /// `n_tax x p`. Integer when every residual degree of freedom is whole.
    pub dof: Vec<f64>,
    /// Whether `dof` should come back as `INTSXP`.
    pub dof_is_int: bool,
    /// `n_tax` matrices of `p x p`, flattened row-major per taxon.
    ///
    /// Per element, R's own order is column-major; the flattening here is a
    /// transport detail and `emit_vcov` puts it back.
    pub vcov: Vec<f64>,
    /// `n_tax x n_samp`, the centred and theta-adjusted log abundance.
    pub y_bias_crt: Vec<f64>,
    /// `n_samp`, the estimated sampling fractions.
    pub theta_hat: Vec<f64>,
    /// `p`: the E-M estimator of the bias.
    pub delta_em: Vec<f64>,
    /// `p`: the WLS estimator of the bias.
    pub delta_wls: Vec<f64>,
    /// `p`: `1 / wls_deno`, the variance of the bias estimate.
    pub var_delta: Vec<f64>,
    /// `p`: `quantile(var_hat, s0_perc)`, added to every variance.
    pub s02: Vec<f64>,
    /// `n_tax x n_samp`, the reported count table.
    pub o2: Vec<f64>,
    /// `nrow(O2)`, the reported taxa -- *not* the bias-estimation taxa.
    pub taxon_names: PackedNames,
    /// `nrow(O2)`, the reported taxa only.
    pub taxa: Vec<usize>,
    /// The taxa that estimated the biases: a superset of [`Self::taxa`] whenever
    /// `struc_zero` flagged anything.
    pub taxa_bias: Vec<usize>,
    /// `ncol(O2)`.
    pub sample_names: PackedNames,
    /// `colnames(x)`.
    pub fix_eff: PackedNames,
}

impl CorePayloads {
    /// Every shape and length check, before a single value is read.
    ///
    /// Same rule as [`crate::transport::Request::validate`]: the first thing that
    /// is wrong is the thing reported, not a downstream consequence of it.
    pub fn validate(&self) -> Result<()> {
        let n_tax = self.n_tax;
        let n_samp = self.n_samp;
        let p = self.p;
        if n_tax == 0 || n_samp == 0 || p == 0 {
            return Err(BridgeError::ShapeMismatch {
                field: "core payloads",
                rows: n_tax,
                cols: n_samp,
                len: n_tax * n_samp,
            });
        }
        for (field, got, want) in [
            ("beta_hat", self.beta_hat.len(), n_tax * p),
            ("var_hat", self.var_hat.len(), n_tax * p),
            ("dof", self.dof.len(), n_tax * p),
            ("vcov_hat", self.vcov.len(), n_tax * p * p),
            ("y_bias_crt", self.y_bias_crt.len(), n_tax * n_samp),
            ("theta_hat", self.theta_hat.len(), n_samp),
            ("O2", self.o2.len(), n_tax * n_samp),
            ("delta_em", self.delta_em.len(), p),
            ("delta_wls", self.delta_wls.len(), p),
            ("var_delta", self.var_delta.len(), p),
            ("s02", self.s02.len(), p),
        ] {
            if got != want {
                return Err(BridgeError::ShapeMismatch {
                    field,
                    rows: n_tax,
                    cols: p,
                    len: got,
                });
            }
        }
        if self.taxa.len() != n_tax {
            return Err(BridgeError::LengthMismatch {
                field: "taxa",
                expected: n_tax,
                got: self.taxa.len(),
            });
        }
        // `taxa_bias` may be a superset -- that is the point of keeping it -- but
        // it can never be shorter, and an empty one would mean no bias was
        // estimated at all.
        if self.taxa_bias.len() < self.taxa.len() {
            return Err(BridgeError::LengthMismatch {
                field: "taxa_bias",
                expected: self.taxa.len(),
                got: self.taxa_bias.len(),
            });
        }
        if self.taxon_names.len() != n_tax {
            return Err(BridgeError::LengthMismatch {
                field: "taxon_names",
                expected: n_tax,
                got: self.taxon_names.len(),
            });
        }
        if self.sample_names.len() != n_samp {
            return Err(BridgeError::LengthMismatch {
                field: "sample_names",
                expected: n_samp,
                got: self.sample_names.len(),
            });
        }
        if self.fix_eff.len() != p {
            return Err(BridgeError::LengthMismatch {
                field: "fix_eff",
                expected: p,
                got: self.fix_eff.len(),
            });
        }
        Ok(())
    }

    /// Copy every payload into owned storage, then hand the copies out.
    ///
    /// The clones are the point, and the same reason as in
    /// [`crate::owned`]: the source buffers belong to R, and anything handed to a
    /// parallel stage has to be owned. The clones are returned rather than kept so
    /// that a caller cannot accidentally hold a second copy for longer than the
    /// call.
    pub fn to_owned(&self) -> OwnedPayloads {
        OwnedPayloads {
            n_tax: self.n_tax,
            n_samp: self.n_samp,
            p: self.p,
            beta_hat: self.beta_hat.clone(),
            var_hat: self.var_hat.clone(),
            dof: self.dof.clone(),
            dof_is_int: self.dof_is_int,
            vcov: self.vcov.clone(),
            y_bias_crt: self.y_bias_crt.clone(),
            theta_hat: self.theta_hat.clone(),
            delta_em: self.delta_em.clone(),
            delta_wls: self.delta_wls.clone(),
            var_delta: self.var_delta.clone(),
            s02: self.s02.clone(),
            o2: self.o2.clone(),
        }
    }

    /// Unflatten taxon `r`'s covariance back into R's column-major order.
    ///
    /// The wire order is row-major per taxon because it is the order that makes
    /// `vcov[r*p*p + c]` readable in a debug dump. R wants column-major, so
    /// element `(i, j)` of the matrix is at `r*p*p + i*p + j`. A wrong choice here
    /// would produce a symmetric-looking wrong answer whenever `p == 2`, which is
    /// the hardest case to notice by eye.
    pub fn vcov_element(&self, taxon: usize, i: usize, j: usize) -> f64 {
        self.vcov[taxon * self.p * self.p + i * self.p + j]
    }

    /// `true` if the covariance block is symmetric to the last bit.
    ///
    /// This is a *check*, not an assumption: the original forms
    /// `XTX_inv %*% sigma %*% XTX_inv`, which is symmetric in exact arithmetic and
    /// very nearly so in floating point. A transport bug that transposed a block
    /// would leave it symmetric and would not be caught by this test -- which is
    /// precisely why the check exists anyway, for the case where the *ordering*
    /// went wrong and the values did not.
    pub fn vcov_is_symmetric(&self) -> bool {
        for r in 0..self.n_tax {
            for i in 0..self.p {
                for j in 0..self.p {
                    if self.vcov_element(r, i, j).to_bits() != self.vcov_element(r, j, i).to_bits()
                    {
                        return false;
                    }
                }
            }
        }
        true
    }
}

/// Owned copies, with the metadata left behind because it is already owned.
///
/// Separated from [`CorePayloads`] so the numeric payload and the names cannot be
/// confused for one another at a call site.
#[derive(Debug, Clone, PartialEq)]
pub struct OwnedPayloads {
    /// `nrow(O2)`.
    pub n_tax: usize,
    /// `ncol(O2)`.
    pub n_samp: usize,
    /// `ncol(x)`.
    pub p: usize,
    /// `n_tax x p`, owned.
    pub beta_hat: Vec<f64>,
    /// `n_tax x p`, owned.
    pub var_hat: Vec<f64>,
    /// `n_tax x p`, owned.
    pub dof: Vec<f64>,
    /// Whether `dof` should be handed back as `INTSXP`.
    pub dof_is_int: bool,
    /// `n_tax * p * p`, owned.
    pub vcov: Vec<f64>,
    /// `n_tax x n_samp`, owned.
    pub y_bias_crt: Vec<f64>,
    /// `n_samp`, owned.
    pub theta_hat: Vec<f64>,
    /// `p`, owned.
    pub delta_em: Vec<f64>,
    /// `p`, owned.
    pub delta_wls: Vec<f64>,
    /// `p`, owned.
    pub var_delta: Vec<f64>,
    /// `p`, owned.
    pub s02: Vec<f64>,
    /// `n_tax x n_samp`, owned.
    pub o2: Vec<f64>,
}

/// The whole emitted payload, written into R-allocated buffers.
///
/// One flat `f64` region, one flat `i32` region and one flat `u8` region, with the
/// plan giving each array its offset. Three buffers rather than twelve because
/// R allocates in one place and the plan can be checked against the allocation
/// before anything is written.
#[derive(Debug, Default)]
pub struct Emitted {
    /// Where each plan entry starts in `f64`, or `-1` if it is not a double.
    pub f64_offsets: Vec<i64>,
    /// Where each plan entry starts in `i32`, or `-1`.
    pub i32_offsets: Vec<i64>,
    /// Where each plan entry starts in `u8`, or `-1`.
    pub u8_offsets: Vec<i64>,
    /// How many `f64` the plan needs.
    pub f64_len: usize,
    /// How many `i32` the plan needs.
    pub i32_len: usize,
    /// How many `u8` the plan needs. Zero today; the field is here so that adding
    /// a logical output does not mean changing the struct's shape and every call
    /// site with it.
    pub u8_len: usize,
}

/// Fill `out` from `payloads`, according to `plan`.
///
/// Returns the offsets. Refuses a plan and a payload set that disagree, rather
/// than writing a prefix of the data and leaving the rest of R's buffer
/// uninitialised.
pub fn emit_core_payloads(
    payloads: &CorePayloads,
    plan: &[ArrayDesc],
    out: &mut Emitted,
) -> Result<()> {
    let mut f64_at = 0i64;
    let mut i32_at = 0i64;
    // No plan entry is currently a `u8` array, so this never advances. It is kept
    // because the plan is data and the next logical output is a mask.
    let u8_at = 0i64;
    out.f64_offsets = Vec::with_capacity(plan.len());
    out.i32_offsets = Vec::with_capacity(plan.len());
    out.u8_offsets = Vec::with_capacity(plan.len());

    for d in plan {
        let src: &[f64] = match d.name {
            "theta_hat" => &payloads.theta_hat,
            "beta_hat" => &payloads.beta_hat,
            "var_hat" => &payloads.var_hat,
            "dof" => &payloads.dof,
            "vcov_hat" => &payloads.vcov,
            "y_bias_crt" => &payloads.y_bias_crt,
            "O2" => &payloads.o2,
            "delta_em" => &payloads.delta_em,
            "delta_wls" => &payloads.delta_wls,
            "var_delta" => &payloads.var_delta,
            "s02" => &payloads.s02,
            other => {
                return Err(BridgeError::UnknownStage { name: other });
            }
        };
        if src.len() != d.count {
            return Err(BridgeError::ShapeMismatch {
                field: match d.name {
                    "vcov_hat" => "vcov_hat",
                    other => other,
                },
                rows: d.rows,
                cols: d.cols,
                len: src.len(),
            });
        }
        match d.rtype {
            RType::Real | RType::ListOfRealMatrix => {
                // A list of matrices travels as its flattened doubles; R rebuilds
                // the list, and each element's `dimnames`, in `assemble.R`.
                out.f64_offsets.push(f64_at);
                f64_at += d.count as i64;
                out.i32_offsets.push(-1);
                out.u8_offsets.push(-1);
            }
            RType::Int => {
                out.f64_offsets.push(-1);
                out.i32_offsets.push(i32_at);
                i32_at += d.count as i64;
                out.u8_offsets.push(-1);
            }
            RType::Logical | RType::Raw | RType::Character => {
                return Err(BridgeError::UnknownStage { name: d.name });
            }
        }
    }
    out.f64_len = f64_at as usize;
    out.i32_len = i32_at as usize;
    out.u8_len = u8_at as usize;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn payloads(n_tax: usize, n_samp: usize, p: usize) -> CorePayloads {
        let f = |k: usize, base: f64| (0..k).map(|i| base + i as f64).collect::<Vec<f64>>();
        CorePayloads {
            n_tax,
            n_samp,
            p,
            beta_hat: f(n_tax * p, 1.0),
            var_hat: f(n_tax * p, 2.0),
            dof: f(n_tax * p, 3.0),
            dof_is_int: true,
            vcov: f(n_tax * p * p, 4.0),
            y_bias_crt: f(n_tax * n_samp, 5.0),
            theta_hat: f(n_samp, 6.0),
            delta_em: f(p, 7.0),
            delta_wls: f(p, 8.0),
            var_delta: f(p, 9.0),
            s02: f(p, 10.0),
            o2: f(n_tax * n_samp, 11.0),
            taxon_names: PackedNames::pack(["T1", "T2", "T3", "T4"]),
            taxa: vec![0, 1, 2, 3],
            taxa_bias: vec![0, 1, 2, 3],
            sample_names: PackedNames::pack(["S1", "S2", "S3", "S4", "S5", "S6"]),
            fix_eff: PackedNames::pack(["(Intercept)", "g2", "g3"]),
        }
    }

    #[test]
    fn the_plan_covers_every_payload_exactly_once() {
        let pay = payloads(4, 6, 3);
        pay.validate().unwrap();
        let plan = core_output_plan(4, 6, 3, true).unwrap();
        let names: Vec<&str> = plan.iter().map(|d| d.name).collect();
        assert_eq!(
            names,
            vec![
                "theta_hat",
                "beta_hat",
                "var_hat",
                "dof",
                "vcov_hat",
                "delta_em",
                "delta_wls",
                "var_delta",
                "s02",
                "y_bias_crt",
                "O2"
            ]
        );
        // Every plan entry's `count` must be exactly its payload's length, or
        // `emit` would refuse.
        let mut emitted = Emitted::default();
        emit_core_payloads(&pay, &plan, &mut emitted).unwrap();
        assert_eq!(emitted.i32_len, 4 * 3, "only dof is integer");
        assert_eq!(
            emitted.f64_len,
            6 + 12 + 12 + 36 + 3 + 3 + 3 + 3 + 24 + 24,
            "every other array is double"
        );
    }

    #[test]
    fn dof_is_reported_as_integer_when_it_should_be() {
        let plan = core_output_plan(4, 6, 3, true).unwrap();
        let d = plan.iter().find(|d| d.name == "dof").unwrap();
        assert_eq!(d.rtype, RType::Int);
        let plan = core_output_plan(4, 6, 3, false).unwrap();
        let d = plan.iter().find(|d| d.name == "dof").unwrap();
        assert_eq!(d.rtype, RType::Real);
    }

    /// A plan whose only integer entry is `dof`, with `dof` not integer, has *no*
    /// integer region at all.
    ///
    /// Which means `i32_len == 0`, and `i32_len == 0` used to be refused by the C
    /// shim's "the plan reserved something" check. A case with a fractional
    /// `df.residual` produces exactly this, so it is pinned here rather than
    /// discovered by a campaign.
    #[test]
    fn a_double_dof_leaves_the_integer_region_empty() {
        let mut pay = payloads(4, 6, 3);
        pay.dof_is_int = false;
        let plan = core_output_plan(pay.n_tax, pay.n_samp, pay.p, false).unwrap();
        let mut e = Emitted::default();
        emit_core_payloads(&pay, &plan, &mut e).unwrap();
        assert_eq!(e.i32_len, 0);
        assert!(e.f64_len > 0);
        assert!(plan.iter().all(|d| d.rtype != RType::Int));
        // And the double `dof` sits in the f64 region where it belongs.
        let at = plan.iter().position(|d| d.name == "dof").unwrap();
        assert!(e.f64_offsets[at] >= 0);
        assert_eq!(e.i32_offsets[at], -1);
    }

    #[test]
    fn vcov_is_the_only_row_major_entry() {
        let plan = core_output_plan(4, 6, 3, true).unwrap();
        for d in &plan {
            let want = if d.name == "vcov_hat" {
                Layout::RowMajor
            } else {
                Layout::ColumnMajor
            };
            assert_eq!(d.layout, want, "{} has the wrong layout", d.name);
        }
        // And the row-major entry is the one that is a list, so the two
        // descriptions of "unusual" cannot drift apart.
        for d in &plan {
            assert_eq!(
                d.layout == Layout::RowMajor,
                d.rtype == RType::ListOfRealMatrix
            );
        }
    }

    #[test]
    fn vcov_stays_a_list_of_matrices() {
        let plan = core_output_plan(4, 6, 3, true).unwrap();
        let v = plan.iter().find(|d| d.name == "vcov_hat").unwrap();
        assert_eq!(v.rtype, RType::ListOfRealMatrix);
        assert_eq!(v.count, 4 * 3 * 3);
        assert_eq!((v.rows, v.cols), (4, 3));
    }

    #[test]
    fn the_bias_and_reported_taxon_sets_are_separate_and_may_differ() {
        // What `struc_zero` produces: O1 keeps four taxa so theta can be estimated
        // from them, O2 reports three. `n_tax` is `nrow(O2)`, so the reported
        // payload is three rows wide while `taxa_bias` still names four input
        // rows.
        let mut pay = payloads(3, 6, 3);
        pay.taxa = vec![0, 1, 3];
        pay.taxa_bias = vec![0, 1, 2, 3];
        pay.taxon_names = PackedNames::pack(["T1", "T2", "T4"]);
        pay.validate().unwrap();
        assert!(pay.taxa_bias.len() > pay.taxa.len());

        // The reverse is impossible, and is refused with a named field.
        pay.taxa_bias = vec![2, 3];
        assert!(matches!(
            pay.validate(),
            Err(BridgeError::LengthMismatch {
                field: "taxa_bias",
                ..
            })
        ));

        // And `taxa` must have one entry per reported taxon.
        let mut pay = payloads(4, 6, 3);
        pay.taxa.pop();
        assert!(matches!(
            pay.validate(),
            Err(BridgeError::LengthMismatch { field: "taxa", .. })
        ));
    }

    #[test]
    fn a_short_payload_is_a_named_error_not_a_read_past_the_end() {
        let mut pay = payloads(4, 6, 3);
        pay.beta_hat.pop();
        assert!(matches!(
            pay.validate(),
            Err(BridgeError::ShapeMismatch {
                field: "beta_hat",
                ..
            })
        ));
    }

    #[test]
    fn a_covariance_element_reads_back_in_rs_column_major_order() {
        // 2 taxa, p = 2, filled so every element is distinguishable.
        let mut pay = payloads(2, 2, 2);
        for t in 0..2 {
            for i in 0..2 {
                for j in 0..2 {
                    pay.vcov[t * 4 + i * 2 + j] = (t * 100 + i * 10 + j) as f64;
                }
            }
        }
        assert_eq!(pay.vcov_element(0, 0, 1), 1.0);
        assert_eq!(pay.vcov_element(0, 1, 0), 10.0);
        assert_eq!(pay.vcov_element(1, 1, 1), 111.0);
    }

    #[test]
    fn a_transposed_covariance_is_not_symmetric_and_says_so() {
        let mut pay = payloads(2, 2, 3);
        for t in 0..2 {
            for i in 0..3 {
                for j in 0..3 {
                    // asymmetric on purpose
                    pay.vcov[t * 9 + i * 3 + j] =
                        (i * 3 + j) as f64 + if i == j { t as f64 } else { 0.5 };
                }
            }
        }
        assert!(!pay.vcov_is_symmetric());
        // A symmetric one passes, which is what the real payload does.
        let mut sym = payloads(2, 2, 3);
        for t in 0..2 {
            for i in 0..3 {
                for j in 0..3 {
                    sym.vcov[t * 9 + i * 3 + j] = (i + j) as f64 + t as f64;
                }
            }
        }
        assert!(sym.vcov_is_symmetric());
    }
}

// ---------------------------------------------------------------------------
// The centring reduction, asked of R rather than re-derived
// ---------------------------------------------------------------------------

use ancombc2_core::reduce::Reductions;
use ancombc2_core::workspace::RMatrix;

/// A reduction implemented by the caller.
///
/// `out` receives `rows` values for a row-wise reduction and `cols` for a
/// column-wise one. `scratch` is `rows * cols` doubles the implementation may use
/// as a flat copy; it exists so a caller that cannot build a matrix view of
/// `x` directly still has somewhere to put the bytes.
pub type RedFn =
    unsafe extern "C" fn(x: *const f64, rows: usize, cols: usize, out: *mut f64) -> i32;

/// The four reductions [`RBackedReductions`] dispatches to.
///
/// This is `r/ANCOMBC/src/init.c`'s table, and the two orders must agree — the C
/// side fills the fields positionally.
#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub struct Reducer {
    /// `rowMeans(x, na.rm = TRUE)`.
    pub row_means: RedFn,
    /// `colMeans(x, na.rm = TRUE)`.
    pub col_means: RedFn,
    /// `colSums(x, na.rm = TRUE)`.
    pub col_sums: RedFn,
    /// `rowSums(x, na.rm = TRUE)`.
    pub row_sums: RedFn,
}

/// R's own `rowSums`/`colSums`/`rowMeans`/`colMeans`, called across the boundary.
///
/// # Why the round trip
///
/// R accumulates these in C `long double`, which is 64 bits of mantissa on x86-64
/// and binary128 on aarch64 — and `double` on Windows. Re-deriving that in Rust
/// would pin one platform's answer as the oracle's, and a 106-bit compensated sum
/// would round differently at every step. Calling `base::rowMeans` is exact by
/// construction on whichever machine R runs on, and it is the same function the
/// oracle called.
///
/// The cost is a round trip per reduction. Preprocessing makes four of them, so
/// this is not a hot path; the hot loops are in S09-S11 and they reduce over
/// already-centred data.
#[derive(Debug, Clone, Copy)]
pub struct RBackedReductions {
    /// The four entry points.
    pub red: Reducer,
}

impl RBackedReductions {
    /// The reductions `init.c` installs.
    ///
    /// # Safety
    ///
    /// Must be called after `R_init_ANCOMBC` has run, i.e. after the package is
    /// loaded. Calling it earlier returns a table of stubs that fail every
    /// reduction rather than a table of null pointers, so the mistake shows up as
    /// a named error at the first reduction rather than as a jump to address zero.
    pub unsafe fn linked() -> Self {
        let mut red = Reducer {
            row_means: reducer_stub,
            col_means: reducer_stub,
            col_sums: reducer_stub,
            row_sums: reducer_stub,
        };
        if ancombc2_rb_reducer(&mut red) != crate::abi::RB_OK {
            panic!("the R reduction table could not be read; the package's init did not run");
        }
        Self { red }
    }

    fn call(&self, f: RedFn, m: &RMatrix, len: usize) -> Vec<f64> {
        let mut out = vec![0.0f64; len];
        let status = unsafe { f(m.data.as_ptr(), m.rows, m.cols, out.as_mut_ptr()) };
        if status != crate::abi::RB_OK {
            // A failed reduction must not be papered over with an `f64` fallback:
            // that would silently reintroduce the difference this type exists to
            // remove, and it would do it invisibly. The name is in the message so
            // the caller can see which primitive went wrong.
            panic!(
                "R's reduction failed with status {status}; refusing to fall back to an \\
                 f64 sum, which is what R does not do"
            );
        }
        out
    }
}

impl Reductions for RBackedReductions {
    fn row_means_na_rm(&self, m: &RMatrix) -> Vec<f64> {
        self.call(self.red.row_means, m, m.rows)
    }

    fn col_means_na_rm(&self, m: &RMatrix) -> Vec<f64> {
        self.call(self.red.col_means, m, m.cols)
    }

    fn col_sums_na_rm(&self, m: &RMatrix) -> Vec<f64> {
        self.call(self.red.col_sums, m, m.cols)
    }

    fn row_sums_na_rm(&self, m: &RMatrix) -> Vec<f64> {
        self.call(self.red.row_sums, m, m.rows)
    }
}

/// Stands in for an uninstalled reduction, so a missing table is a status code
/// and not a jump through a null pointer.
unsafe extern "C" fn reducer_stub(_x: *const f64, _r: usize, _c: usize, _o: *mut f64) -> i32 {
    crate::abi::RB_ERR_ARGS
}

extern "C" {
    /// `init.c`'s `ancombc2_rb_reducer`. Fills the table rather than returning it
    /// by value, so the two sides cannot disagree about the struct's size.
    fn ancombc2_rb_reducer(out: *mut Reducer) -> i32;
}

// ---------------------------------------------------------------------------
// The preprocessing stages, one at a time
// ---------------------------------------------------------------------------

/// Every array preprocessing produces, in the order `.ancombc2_core` produces it.
///
/// This exists because S08's acceptance is about *stage* arrays, and the only way
/// to show a stage is exact is to hand its array across the boundary and compare
/// the bytes. A whole-result comparison cannot localise a difference to a stage, and
/// "the numbers came out the same at the end" says nothing about whether the
/// centring or the covariance was at fault.
///
/// # What is here, and why
///
/// `.data_core` is `.ancombc2_core` step 0: prevalence, the taxon screen, the
/// library-size screen, and the retained indices. Then step 1's
/// `O1 = data + pseudo`, `o1 = log(O1)`, `o1[is.infinite(o1)] = NA`,
/// `y1 = o1 - rowMeans(o1, na.rm = TRUE)` — with the row means as their own stage,
/// because that reduction is the one place R's accumulator width is observable.
/// Then the same four for `O2`, and then `.get_struc_zero`'s per-group
/// prevalence and lower bound.
///
/// `theta_hat` is absent: it needs `beta1`, which is S09.
#[derive(Debug, Clone)]
pub struct PreprocessStages {
    /// `nrow(data)`.
    pub prevalence1: Vec<f64>,
    /// `which(prevalence >= prv_cut)`, 1-based.
    pub tax_keep1: Vec<i64>,
    /// `colSums(feature_table[tax_keep1, ], na.rm = TRUE)` -- over every sample of
    /// the retained taxa, because the reference computes it before it subsets the
    /// sample axis.
    pub lib_size1: Vec<f64>,
    /// `which(lib_size >= lib_cut)`, 1-based.
    pub samp_keep1: Vec<i64>,
    /// `data + pseudo` on the retained sub-table.
    pub o1: RMatrix,
    /// `log(O1)` with `is.infinite` mapped to `NA_real_`.
    pub log1: RMatrix,
    /// `rowMeans(log1, na.rm = TRUE)`.
    pub means1: Vec<f64>,
    /// `log1 - means1`.
    pub y1: RMatrix,
    /// The second `.data_core` pass's `tax_keep`, **as `.data_core` reports it**:
    /// positions within the structural-zero-retained table, not within `aggregate`.
    ///
    /// This is the reference's own surprise. When `tax_keep` is supplied,
    /// `ancombc_prep.R:19-22` subsets first and *then* calls `which()`:
    ///
    /// ```r
    /// feature_table = feature_table[tax_keep, , drop = FALSE]
    /// prevalence = prevalence_fun(feature_table)
    /// tax_keep = which(prevalence >= prv_cut)
    /// ```
    ///
    /// so `core2$tax_keep` counts from 1 within the subset. Nothing downstream reads it
    /// -- the reference only uses it to subset again -- but a stage array named after it
    /// has to hold what it holds, or the diagnostic measures something else.
    ///
    /// [`Self::tax_keep2_absolute`] is the same set in `aggregate`'s own coordinates,
    /// which is what a Rust caller wants.
    pub tax_keep2: Vec<i64>,
    /// The second pass's retained taxa in `aggregate`'s coordinates, 1-based. Equals
    /// `zero_keep[tax_keep2 - 1L]`.
    pub tax_keep2_absolute: Vec<i64>,
    /// The second pass's `samp_keep`. Reused from the first, not recomputed.
    pub samp_keep2: Vec<i64>,
    /// `aggregate_data + pseudo` on the retained sub-table.
    pub o2: RMatrix,
    /// `log(O2)` with `is.infinite` mapped to `NA_real_`.
    pub log2: RMatrix,
    /// `rowMeans(log2, na.rm = TRUE)`.
    pub means2: Vec<f64>,
    /// `log2 - means2`.
    pub y2: RMatrix,
    /// `taxa_without_structural_zeros()`: the taxa the screen kept, 1-based.
    pub zero_keep: Vec<i64>,
    /// `sweep(present %*% G, 2, n_g, "/")`, `n_tax x n_group`, row-major.
    pub group_prevalence: Vec<f64>,
    /// `rowSums(!is.na(feature_mat)) %*% G`, `n_tax x n_group`.
    pub group_size: Vec<f64>,
    /// `p_hat - 1.96 * sqrt(p_hat * (1 - p_hat) / samp_size)`.
    pub group_lower: Vec<f64>,
    /// The flags, as 0/1, `n_tax x n_group`.
    pub zero_ind: Vec<i64>,
}

/// `counts + pseudo`, the reference's `O = data + pseudo`.
///
/// A separate buffer rather than a mutation, because `counts` is the retained
/// count table and `O` is a different array; the reference has both live at once.
///
/// [`r_add`], not `+`. This is the first arithmetic the counts see, so it is where
/// the `NA_real_` payload is destroyed if plain `+` is used: `NA_real_ + 0.0` is
/// `0x7ff8000000000002` on x86-64, and every later stage then propagates *that*
/// rather than the oracle's `NA_real_`. Fixing it in `r_log` alone would have left
/// this one step earlier in the pipeline broken, which is the argument for checking
/// stage arrays rather than the result.
fn add_scalar(counts: &ancombc2_core::preprocess::CountMatrix, pseudo: f64) -> RMatrix {
    RMatrix::from_row_major(
        counts.n_taxa,
        counts.n_samp,
        counts
            .data
            .iter()
            .map(|&v| ancombc2_core::reduce::r_add(v, pseudo))
            .collect(),
    )
}

/// The inputs to [`preprocess_stages`].
///
/// A struct rather than ten positional arguments. The reference's
/// `.ancombc2_core` takes thirty formals, so a list is the only readable shape for
/// this many knobs, and a caller cannot silently swap `lib_cut` for `pseudo`.
#[derive(Clone, Copy)]
pub struct PreprocessInputs<'a> {
    /// `feature_table`: the sample-level counts.
    pub data: &'a ancombc2_core::preprocess::CountMatrix,
    /// `feature_table_aggregate`: the aggregated counts. Same shape, different
    /// contents, and the structural-zero screen reads *this* one.
    pub aggregate: &'a ancombc2_core::preprocess::CountMatrix,
    /// `struc_zero`. When `false`, every taxon is kept.
    pub struc_zero: bool,
    /// The `n_samp x n_group` 0/1 indicator `G`, column-major as R stores it.
    /// Required when [`Self::struc_zero`] is set.
    pub group_indicator: Option<&'a [f64]>,
    /// `nlevels(factor(meta_data[, group]))`.
    pub n_groups: usize,
    /// `prv_cut`, compared with `>=` against `prevalence`.
    pub prv_cut: f64,
    /// `lib_cut`, compared with `>=` against the retained-taxa column sums.
    pub lib_cut: f64,
    /// Added to every count *before* the log.
    pub pseudo: f64,
    /// `neg_lb`: also flag a taxon whose per-group lower bound is non-positive.
    pub neg_lb: bool,
    /// The centring reduction. `RBackedReductions` in the R session, so the row
    /// means are R's own `long double` accumulation.
    pub red: &'a dyn ancombc2_core::reduce::Reductions,
}

/// Compute every preprocessing stage from the same inputs `.ancombc2_core` uses.
pub fn preprocess_stages(inp: &PreprocessInputs<'_>) -> Result<PreprocessStages> {
    let PreprocessInputs {
        data,
        aggregate,
        struc_zero,
        group_indicator,
        n_groups,
        prv_cut,
        lib_cut,
        pseudo,
        neg_lb,
        red,
    } = *inp;
    // `ancombc2.R:453-462`: `.get_struc_zero` runs on the *unfiltered aggregate*
    // table against the *unfiltered* `meta_data`, before either `.data_core` pass.
    // Its `tax_keep` is then handed to the second pass. So the screen's indices are
    // into `aggregate`, not into `data` -- and the two tables have the same shape
    // but different *contents*, so using `data` here would silently screen the
    // wrong table.
    let mut zero_keep: Vec<usize> = (0..aggregate.n_taxa).collect();
    let mut group_prevalence = Vec::new();
    let mut group_size = Vec::new();
    let mut group_lower = Vec::new();
    let mut zero_ind = Vec::new();
    if struc_zero {
        let g = group_indicator.ok_or_else(|| {
            BridgeError::Args(
                "struc_zero needs a group; the reference stops with \"group must be \
                 specified\" before reaching .get_struc_zero"
                    .into(),
            )
        })?;
        let mut zs = ZeroStages::default();
        structural_zero_stages(aggregate, g, n_groups, neg_lb, &mut zs)?;
        zero_keep = zs.keep;
        group_prevalence = zs.prevalence;
        group_size = zs.size;
        group_lower = zs.lower;
        zero_ind = zs.ind;
    }

    // `ancombc2.R:466-469`: the first `.data_core` pass takes `tax_keep = NULL` and
    // `samp_keep = NULL`, so both are computed here from scratch.
    let first = data
        .filter(prv_cut, lib_cut)
        .map_err(|e| BridgeError::Args(e.to_string()))?;
    let counts1 = first.counts.clone();
    let prevalence1 = data.prevalence();
    // `lib_size` is `colSums` of the *taxon-filtered but not sample-filtered* table.
    // The reference computes it inside `if (is.null(samp_keep))` and only subsets the
    // sample axis afterwards, so it is a column sum over every sample of the retained
    // taxa. Reading it off `counts1` -- which is already sample-filtered -- gives the
    // right *decision* (the filter used the full table) and the wrong *array*, which
    // is worse: the array is what a later stage would consume.
    let counts1_all_samples = {
        let all: Vec<usize> = (0..data.n_samp).collect();
        data.select(&first.taxa, &all)
    };
    let lib_size1 = counts1_all_samples.library_sizes();

    // The retained count table, plus pseudo, logged and centred -- step 1.
    let o1 = add_scalar(&counts1, pseudo);
    let log1 = log_and_center(&o1, red);
    let means1 = red.row_means_na_rm(&log1);
    let mut y1 = log1.clone();
    y1.sub_rows_in_place(&means1);

    // `ancombc2.R:473-476`: the second pass is on the aggregate table, with
    // `tax_keep` from the structural-zero screen and `samp_keep` from the first
    // pass. The reference passes `samp_keep = colnames(O1)` -- *names*, which R
    // resolves by column name -- so this uses the first pass's indices. That is the
    // same set of columns only because `feature_table_aggregate` is derived from
    // `feature_table` and keeps its column names and order; see the note on
    // `samp_keep_is_positional_here`.
    let second = aggregate
        .filter_aggregated(
            if struc_zero { Some(&zero_keep) } else { None },
            Some(&first.samples),
            prv_cut,
            lib_cut,
        )
        .map_err(|e| BridgeError::Args(e.to_string()))?;
    let counts2 = second.counts;
    let o2 = add_scalar(&counts2, pseudo);
    let log2 = log_and_center(&o2, red);
    let means2 = red.row_means_na_rm(&log2);
    let mut y2 = log2.clone();
    y2.sub_rows_in_place(&means2);

    let st = PreprocessStages {
        prevalence1,
        tax_keep1: first.taxa.iter().map(|&t| t as i64 + 1).collect(),
        lib_size1,
        samp_keep1: first.samples.iter().map(|&s| s as i64 + 1).collect(),
        o1,
        log1,
        means1,
        y1,
        tax_keep2: {
            // `.data_core` subsets by `tax_keep` and then `which()`es, so its indices
            // count within the subset. `zero_keep` is that subset, sorted and
            // deduplicated, so an absolute index's rank in it is the reported one.
            let subset: Vec<usize> = {
                let mut v = zero_keep.clone();
                v.sort_unstable();
                v.dedup();
                v
            };
            second
                .taxa
                .iter()
                .map(|&t| {
                    subset
                        .binary_search(&t)
                        .map(|rank| rank as i64 + 1)
                        .unwrap_or_else(|_| {
                            // Cannot happen: `second.taxa` is a subset of
                            // `after_struct`, which is `zero_keep`. A miss would be a
                            // silent wrong index, so it is refused here.
                            panic!(
                                "taxon {t} was retained by the second pass but is not in \
                                 the structural-zero subset, so it has no subset-relative \
                                 index"
                            );
                        })
                })
                .collect()
        },
        tax_keep2_absolute: second.taxa.iter().map(|&t| t as i64 + 1).collect(),
        samp_keep2: second.samples.iter().map(|&s| s as i64 + 1).collect(),
        o2,
        log2,
        means2,
        y2,
        zero_keep: zero_keep.iter().map(|&t| t as i64 + 1).collect(),
        group_prevalence,
        group_size,
        group_lower,
        zero_ind,
    };

    Ok(st)
}

/// `log(O)` with `is.infinite` mapped to `NA_real_`.
///
/// [`ancombc2_core::reduce::r_log`], not `ln`: a `NaN` count keeps its payload, and
/// libm does not promise to.
fn log_and_center(o: &RMatrix, _red: &dyn ancombc2_core::reduce::Reductions) -> RMatrix {
    let mut y = o.clone();
    for v in y.data.iter_mut() {
        *v = ancombc2_core::reduce::r_log(*v);
    }
    y.replace_infinite_with_na();
    y
}

/// `.get_struc_zero` on the retained sub-table, split into its three stages.
///
/// The reference forms `p_hat = sweep(present_table %*% G, 2, n_g, "/")`. Both
/// factors are 0/1, so every product and every partial sum is an exact integer and
/// the multiplication order is immaterial — which is why a plain triple loop here
/// is bit-identical to `dgemm` and does not need to call BLAS. `sqrt` and the
/// division are then correctly-rounded IEEE operations in both.
fn structural_zero_stages(
    counts: &ancombc2_core::preprocess::CountMatrix,
    g: &[f64],
    n_groups: usize,
    neg_lb: bool,
    out: &mut ZeroStages,
) -> Result<()> {
    if n_groups == 0 || g.len() != counts.n_samp * n_groups {
        return Err(BridgeError::Args(format!(
            "the group indicator has {} values; a {} x {} matrix needs {}",
            g.len(),
            counts.n_samp,
            n_groups,
            counts.n_samp * n_groups
        )));
    }
    let n_tax = counts.n_taxa;
    let n_samp = counts.n_samp;
    // `G` is `n_samp x n_groups`, column-major as R stores it.
    let gg = |j: usize, k: usize| g[j + k * n_samp];

    let mut present = vec![0.0f64; n_tax * n_groups];
    let mut observed = vec![0.0f64; n_tax * n_groups];
    for i in 0..n_tax {
        for j in 0..n_samp {
            let v = counts.get(i, j);
            // `present_table[is.na(..)] = 0`, so an `NA` is an absence here -- the
            // opposite of prevalence, where an `NA` leaves the ratio entirely.
            let is_present = if v.is_nan() { 0.0 } else { f64::from(v != 0.0) };
            let is_observed = f64::from(!v.is_nan());
            for k in 0..n_groups {
                let m = gg(j, k);
                present[i * n_groups + k] += is_present * m;
                observed[i * n_groups + k] += is_observed * m;
            }
        }
    }

    // `n_g = colSums(G)`: a group's sample count. `G` is 0/1 so this is a count.
    let mut n_g = vec![0.0f64; n_groups];
    for k in 0..n_groups {
        n_g[k] = (0..n_samp).map(|j| gg(j, k)).sum();
    }

    let mut p_hat = vec![0.0f64; n_tax * n_groups];
    for i in 0..n_tax {
        for k in 0..n_groups {
            // `sweep(., 2, n_g, "/")`. A group of size zero divides by zero, which in
            // R is `0/0 = NaN` or `x/0 = Inf`; both are reproduced by the plain
            // division rather than special-cased.
            p_hat[i * n_groups + k] = present[i * n_groups + k] / n_g[k];
        }
    }

    let mut lower = vec![0.0f64; n_tax * n_groups];
    let mut flags = vec![0i64; n_tax * n_groups];
    for i in 0..n_tax {
        for k in 0..n_groups {
            let p = p_hat[i * n_groups + k];
            let sz = observed[i * n_groups + k];
            lower[i * n_groups + k] = p - 1.96 * (p * (1.0 - p) / sz).sqrt();
            // `zero_ind = (p_hat == 0)`, then `if (neg_lb) zero_ind[p_hat_lo <= 0]`.
            let mut z = p == 0.0;
            if neg_lb && lower[i * n_groups + k] <= 0.0 {
                z = true;
            }
            flags[i * n_groups + k] = i64::from(z);
        }
    }

    // The reference's own screen is `all(zero_ind[, -1] == FALSE)`, which is what
    // `StructuralZeros::taxa_without_structural_zeros` computes. Reusing it means
    // the taxon set this stage reports and the one the estimator will use cannot
    // come from two different definitions of "structural zero".
    let mut bits = ancombc2_core::preprocess::ZeroBitSet::new(n_tax * n_groups);
    let mut zero_ind = vec![false; n_tax * n_groups];
    for i in 0..n_tax * n_groups {
        bits.set(i, flags[i] != 0);
        zero_ind[i] = flags[i] != 0;
    }
    let zeros = ancombc2_core::preprocess::StructuralZeros {
        groups: (0..n_groups).map(|k| format!("g{k}")).collect(),
        zero_ind,
        bits,
    };

    out.keep = zeros.taxa_without_structural_zeros();
    out.prevalence = p_hat;
    out.size = observed;
    out.lower = lower;
    out.ind = flags;
    Ok(())
}

/// `.get_struc_zero`'s five arrays, as out-parameters.
///
/// A struct so the call site reads as `.get_struc_zero(&mut zero_stages)` rather than
/// as five pointers being threaded through, and so a new output is one field rather
/// than one more argument.
#[derive(Debug, Default)]
struct ZeroStages {
    /// `which(all(zero_ind[, -1] == FALSE))`.
    keep: Vec<usize>,
    /// `p_hat`, `n_tax x n_group`.
    prevalence: Vec<f64>,
    /// `samp_size`, `n_tax x n_group`.
    size: Vec<f64>,
    /// `p_hat_lo`, `n_tax x n_group`.
    lower: Vec<f64>,
    /// The flags, as 0/1, `n_tax x n_group`.
    ind: Vec<i64>,
}
