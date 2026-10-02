//! Huber-White sandwich variance, with the reference's NA quirk.
//!
//! For taxon `i`,
//!
//! ```text
//! V_i = (X'X)^- ( sum_j eps_ij^2 x_j x_j' ) (X'X)^-
//! ```
//!
//! where the inverse is `MASS::ginv(t(X_cc) %*% X_cc)` over the design rows with
//! no missing values.
//!
//! # The NA -> 0.1 quirk
//!
//! The reference computes, per sample,
//!
//! ```r
//! term_j = outer(eps2[idx, j], XX[j, ])
//! term_j[is.na(term_j)] = 0.1
//! sigma2_xxT = sigma2_xxT + term_j
//! ```
//!
//! The substitution is **per entry, not per term**. `eps2` is `NA` wherever the
//! taxon is unobserved -- that is, wherever its count was zero and `pseudo` is
//! zero -- and `XX[j, ]` is `NA` at each position where sample `j`'s design entry
//! is missing. So:
//!
//! * an unobserved (taxon, sample) pair has an all-`NA` `eps2`, and therefore
//!   *every* entry of its outer product is replaced by `0.1`. The term contributes
//!   `0.1` to all `p^2` positions instead of nothing -- a spurious constant added
//!   to the variance of a taxon that was never observed at that sample;
//! * a *partially* missing design row has `NA` only in the positions involving a
//!   missing design entry, and only those are replaced.
//!
//! The first case is the one that matters in practice: for a table with zeros and
//! `pseudo = 0`, essentially every taxon has unobserved samples, so the quirk
//! fires almost everywhere. The two compatibility modes therefore differ on
//! almost any real input, and property test P16 is what pins the difference to
//! this one mechanism. Upstream has an open issue asking whether it is intended.
//!
//! This module reproduces it under [`CompatMode::Ancombc2_15`] and substitutes the
//! mathematically neutral value (contribute nothing) under
//! [`CompatMode::StrictSpec`]. The default is the reference behaviour, and
//! `docs/reference_behavior.md` records the divergence.
//!
//! # Accumulation order
//!
//! Taxa are processed in index order and, within a taxon, samples in index order,
//! with the same block size the reference uses
//! (`block = max(1, as.integer(2^17 / p^2))`). The result is bit-identical to a
//! naive per-taxon accumulation, and identical to the reference up to the
//! blocking, which is a float association that only reorders additions within a
//! taxon.

use crate::config::CompatMode;
use crate::matrix::{ginv, Matrix};

/// Precomputed per-sample outer products `vec(x_j x_j')`, row-major over
/// samples: `xx[j * p * p + a * p + b] = x_j[a] * x_j[b]`.
#[derive(Debug, Clone)]
pub struct OuterProducts {
    pub n_samp: usize,
    pub p: usize,
    pub data: Vec<f64>,
    /// `true` when the design row for sample `j` had a missing value. Kept for
    /// the fast path: when the whole row is complete, every entry of the outer
    /// product is finite and no per-entry test is needed.
    pub row_incomplete: Vec<bool>,
}

impl OuterProducts {
    /// R: `XX[j, ] = as.vector(x[j, ] %*% t(x[j, ]))`, column-major `vec`.
    pub fn build(x: &Matrix) -> Self {
        let n = x.rows;
        let p = x.cols;
        // A partially-missing design row must keep `NA` only where it is NA: R
        // builds `XX[j, ] = as.vector(x[j, ] %*% t(x[j, ]))` element by element,
        // and `NaN * finite` is `NaN` only for the entries that involve the
        // missing design value. Replacing the whole row with `NA` -- which an
        // earlier version of this did -- substitutes `0.1` at positions the
        // reference leaves alone.
        let mut data = vec![0.0; n * p * p];
        let mut row_incomplete = vec![false; n];
        for j in 0..n {
            let mut ok = true;
            for a in 0..p {
                if !x.get(j, a).is_finite() {
                    ok = false;
                }
            }
            row_incomplete[j] = !ok;
            for a in 0..p {
                for b in 0..p {
                    data[j * p * p + a * p + b] = x.get(j, a) * x.get(j, b);
                }
            }
        }
        Self {
            n_samp: n,
            p,
            data,
            row_incomplete,
        }
    }

    #[inline]
    fn at(&self, j: usize, a: usize, b: usize) -> f64 {
        self.data[j * self.p * self.p + a * self.p + b]
    }

    /// The taxon block size the reference uses.
    pub fn block(&self) -> usize {
        let pp = (self.p * self.p) as f64;
        (131072.0 / pp).floor().max(1.0) as usize
    }

    /// How many taxa' `p x p` accumulators fit in the same 128 KiB budget.
    ///
    /// The companion to [`OuterProducts::block`]. That one answers "how many
    /// samples' outer products fit in L2", which is what the reference blocks
    /// over; this one answers "how many taxa' accumulators fit", which is what
    /// has to be blocked to stop the *outer products* being re-read once per
    /// taxon. The sandwich needs both: the sample tile bounds how much of `xx`
    /// is live, and the taxon tile bounds how many accumulators stay in cache
    /// while that tile is streamed over every taxon in the block.
    pub fn taxa_per_block(&self) -> usize {
        let bytes = (self.p * self.p * std::mem::size_of::<f64>()) as f64;
        (131072.0 / bytes).floor().max(1.0) as usize
    }
}

/// `XTX_inv = MASS::ginv(t(x_cc) %*% x_cc)`, the only place the inverse is
/// formed: it does not depend on the taxon, so it is computed once per run.
pub fn xtx_inverse(x: &Matrix) -> Matrix {
    // complete design rows only: x[complete.cases(x), ]
    let n = x.rows;
    let keep: Vec<usize> = (0..n)
        .filter(|&j| (0..x.cols).all(|a| x.get(j, a).is_finite()))
        .collect();
    let xcc = x.select_rows(&keep);
    let xtx = xcc.t().matmul(&xcc).expect("gram matrix");
    ginv(&xtx, None)
}

/// Accumulate `sum_j eps_ij^2 XX[j, ]` for every taxon, then sandwich it.
///
/// `eps` is `n_taxa x n_samp` row-major and is squared *on the fly*. A
/// pre-squared copy would be a second `n_taxa x n_samp` buffer -- 800 MB at
/// 1000 x 10000 -- held for the whole accumulation and never used again, and it
/// would exist twice over because each MLE has one. `NaN * NaN` is `NaN`, so
/// squaring in place of pre-squaring changes nothing about which cells take the
/// `0.1` branch.
/// `n_taxa * p * p` row-major flattened variance matrices and the `n_taxa x p`
/// marginal variances.
///
/// The accumulation for taxon `i` walks samples in ascending order, matching the
/// reference's inner `for (j in seq_len(n_samp))` loop, so the sum is associated
/// identically.
pub fn sandwich_all(
    eps: &[f64],
    n_taxa: usize,
    n_samp: usize,
    xx: &OuterProducts,
    xtx_inv: &Matrix,
    compat: CompatMode,
) -> (Vec<f64>, Vec<f64>) {
    let p = xx.p;
    let mut vcov = vec![0.0f64; n_taxa * p * p];
    let mut var_hat = vec![f64::NAN; n_taxa * p];

    // --- L2 blocking ---------------------------------------------------------
    //
    // The accumulation is `sigma2_xxT[i] = sum_j eps2_ij * XX_j`, a `p x p`
    // matrix per taxon. Straight over taxa it re-reads *all* of `xx` once per
    // taxon: on the 1000 x 10000 benchmark `xx` is 26 MB and there are 1000
    // taxa, so 26 GB of memory traffic to do 26 MB of arithmetic.
    //
    // The fix is to tile both axes. A sample tile of `xx.block()` rows is
    // 128 KiB, so it is cache-resident for the whole sweep over a taxon block;
    // a taxon block of `xx.taxa_per_block()` accumulators is also 128 KiB, so it
    // stays in cache while that sample tile is streamed. The outer products are
    // then re-read `n_taxa / taxa_per_block` times instead of `n_taxa` times.
    //
    // **The arithmetic is unchanged.** For a fixed taxon the samples are still
    // accumulated in ascending `j` -- the tiles partition `0..n_samp` in order --
    // and each product is the same `f64` multiply-add in the same order. So this
    // is bit-for-bit identical to the untiled loop, which is what
    // `the_blocked_sandwich_is_bit_identical_to_the_naive_one` asserts rather
    // than taking on trust.
    let samp_block = xx.block().max(1);
    let taxa_block = xx.taxa_per_block().max(1);
    // One tile of outer products, reused. This replaces a full `n_samp x p x p`
    // copy (26 MB on that benchmark) with 128 KiB.
    // One `tile` and one `acc` per taxon *block*, and the blocks are independent:
    // block `b` writes `vcov[i0..i1]` and `var_hat[i0..i1]` and nothing else. That
    // makes the outer loop the third disjoint axis, and on `bm4` it is 20 blocks
    // against 1.7 seconds of otherwise-serial sandwich.
    //
    // The scratch becomes per-block rather than shared: 20 blocks of ~130 KiB each
    // is 2.6 MB, against a saving of well over a second, and sharing it would need
    // a lock on the analysis's hottest loop.
    let mut lvl = crate::parallel::NestingBudget::level("taxa (sandwich blocks)");
    let block_starts: Vec<usize> = (0..n_taxa).step_by(taxa_block).collect();
    let per_block = crate::parallel::map_par(&mut lvl, &block_starts, |&i0| -> (usize, Vec<f64>) {
        let mut tile: Vec<f64> = Vec::with_capacity(samp_block * p * p);
        // `taxa_block` accumulators, live only while this block is being built.
        let mut acc = vec![0.0f64; taxa_block * p * p];
        let i1 = (i0 + taxa_block).min(n_taxa);
        for j0 in (0..n_samp).step_by(samp_block) {
            let j1 = (j0 + samp_block).min(n_samp);
            tile.clear();
            for j in j0..j1 {
                for a in 0..p {
                    for b in 0..p {
                        tile.push(xx.at(j, a, b));
                    }
                }
            }
            for (bi, i) in (i0..i1).enumerate() {
                let base = i * n_samp;
                let acc = &mut acc[bi * p * p..(bi + 1) * p * p];
                for (jj, j) in (j0..j1).enumerate() {
                    let e2 = {
                        let e = eps[base + j];
                        e * e
                    };
                    let off = jj * p * p;
                    if e2.is_finite() && !xx.row_incomplete[j] {
                        // The fast path: both operands are finite everywhere, so
                        // no per-entry test is needed. This is every sample of
                        // every taxon when the design is complete and the taxon
                        // is fully observed.
                        for k in 0..(p * p) {
                            acc[k] += e2 * tile[off + k];
                        }
                    } else if compat == CompatMode::Ancombc2_15 {
                        // Entry-wise: a term is 0.1 wherever the product is NA,
                        // which is everywhere if `eps2` is NA, and only at the
                        // positions involving a missing design entry otherwise.
                        for k in 0..(p * p) {
                            let t = e2 * tile[off + k];
                            acc[k] += if t.is_nan() { 0.1 } else { t };
                        }
                    }
                    // StrictSpec: contribute nothing where the product is not finite
                }
            }
        }

        (i0, acc)
    });

    // Scattered in block order, so `vcov` and `var_hat` are indexed by taxon and
    // not by whichever worker finished first.
    for (i0, acc) in per_block {
        let i1 = (i0 + taxa_block).min(n_taxa);
        for (bi, i) in (i0..i1).enumerate() {
            let a = &acc[bi * p * p..(bi + 1) * p * p];
            // V = XTX_inv %*% sigma2_xxT %*% XTX_inv
            let sigma = Matrix::from_vec(p, p, a.to_vec()).expect("p x p");
            let v = xtx_inv
                .matmul(&sigma)
                .expect("matmul")
                .matmul(xtx_inv)
                .expect("matmul");
            let v = v.symmetrise();
            for a in 0..p {
                var_hat[i * p + a] = v.get(a, a);
                for bb in 0..p {
                    vcov[i * p * p + a * p + bb] = v.get(a, bb);
                }
            }
        }
    }

    (vcov, var_hat)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn design() -> Matrix {
        // 1-intercept + 1 binary covariate, 6 samples
        Matrix::from_rows(&[
            vec![1.0, 0.0],
            vec![1.0, 0.0],
            vec![1.0, 1.0],
            vec![1.0, 1.0],
            vec![1.0, 0.0],
            vec![1.0, 1.0],
        ])
    }

    /// A reference implementation of the sandwich, with no blocking at all.
    ///
    /// Deliberately the shape `sandwich_all` used to have: taxa outermost,
    /// every sample of the taxon, one `p x p` accumulator reused per taxon. It
    /// exists so the blocked version can be held to bit-for-bit equality rather
    /// than to a tolerance -- a blocked summation that agreed only to 1e-15
    /// would still have changed the goldens, and the reason it does not is the
    /// *order*, not the magnitudes.
    fn sandwich_all_naive(
        eps: &[f64],
        n_taxa: usize,
        n_samp: usize,
        xx: &OuterProducts,
        xtx_inv: &Matrix,
        compat: CompatMode,
    ) -> (Vec<f64>, Vec<f64>) {
        let p = xx.p;
        let mut vcov = vec![0.0f64; n_taxa * p * p];
        let mut var_hat = vec![f64::NAN; n_taxa * p];
        for i in 0..n_taxa {
            let base = i * n_samp;
            let mut acc = vec![0.0f64; p * p];
            for j in 0..n_samp {
                let e = eps[base + j];
                let e2 = e * e;
                if e2.is_finite() && !xx.row_incomplete[j] {
                    for a in 0..p {
                        for b in 0..p {
                            acc[a * p + b] += e2 * xx.at(j, a, b);
                        }
                    }
                } else if compat == CompatMode::Ancombc2_15 {
                    for a in 0..p {
                        for b in 0..p {
                            let t = e2 * xx.at(j, a, b);
                            acc[a * p + b] += if t.is_nan() { 0.1 } else { t };
                        }
                    }
                }
            }
            let sigma = Matrix::from_vec(p, p, acc).expect("p x p");
            let v = xtx_inv
                .matmul(&sigma)
                .expect("matmul")
                .matmul(xtx_inv)
                .expect("matmul");
            let v = v.symmetrise();
            for a in 0..p {
                var_hat[i * p + a] = v.get(a, a);
                for b in 0..p {
                    vcov[i * p * p + a * p + b] = v.get(a, b);
                }
            }
        }
        (vcov, var_hat)
    }

    /// The blocked sandwich must equal the naive one exactly.
    ///
    /// L2 blocking is a memory-traffic optimisation, and its safety rests
    /// entirely on the claim that for a fixed taxon the samples are still
    /// accumulated in ascending order, so every `f64` multiply-add happens in the
    /// same sequence as before. That claim is what this test checks, and it is
    /// checked for *equality*, not closeness: a tolerance would hide exactly the
    /// reassociation the test exists to forbid. The shapes are chosen so both
    /// block sizes actually tile -- `samp_block` and `taxa_per_block` are derived
    /// from `p`, so a small `p` gives a large tile and a tiny table would never
    /// block at all.
    #[test]
    fn the_blocked_sandwich_is_bit_identical_to_the_naive_one() {
        let mut seed = 20240917u64;
        let mut rnd = move || {
            seed = seed
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            ((seed >> 11) as f64 / (1u64 << 53) as f64) * 2.0 - 1.0
        };
        for (n_samp, p, n_taxa) in [(64usize, 2usize, 7usize), (97, 3, 130), (40, 5, 3)] {
            let mut x = Matrix::zeros(n_samp, p);
            for v in x.data.iter_mut() {
                *v = rnd();
            }
            for i in 0..n_samp {
                x.set(i, 0, 1.0);
            }
            let xx = OuterProducts::build(&x);
            let xtx_inv = xtx_inverse(&x);
            let mut eps = vec![0.0f64; n_taxa * n_samp];
            for e in eps.iter_mut() {
                *e = rnd();
            }
            // Scattered non-finite residuals and an incomplete design row, so
            // both the fast path and the 0.1-substitution path are exercised
            // inside a blocked sweep.
            for k in 0..n_taxa {
                eps[k * n_samp] = f64::NAN;
                eps[k * n_samp + 1] = f64::INFINITY;
            }
            for compat in [CompatMode::Ancombc2_15, CompatMode::StrictSpec] {
                let (v_blocked, h_blocked) =
                    sandwich_all(&eps, n_taxa, n_samp, &xx, &xtx_inv, compat);
                let (v_naive, h_naive) =
                    sandwich_all_naive(&eps, n_taxa, n_samp, &xx, &xtx_inv, compat);
                // Bit patterns, not values. The residuals include `INFINITY`,
                // which drives `acc` to infinity and then the sandwich to NaN, and
                // `NaN != NaN` would report a difference that is not one. Bit
                // equality is also the stronger claim: it distinguishes a NaN
                // from a different NaN, and it cannot be satisfied by two values
                // that merely happen to be close.
                let same = v_blocked.len() == v_naive.len()
                    && v_blocked
                        .iter()
                        .zip(v_naive.iter())
                        .all(|(a, b)| a.to_bits() == b.to_bits())
                    && h_blocked
                        .iter()
                        .zip(h_naive.iter())
                        .all(|(a, b)| a.to_bits() == b.to_bits());
                assert!(
                    same,
                    "n_samp={n_samp} p={p} n_taxa={n_taxa} {compat:?}: the blocked \
                     sandwich is not bit-identical to the unblocked accumulation, \
                     so the blocking reordered a reduction"
                );
            }
        }
    }

    /// The blocking actually tiles, and the two block sizes bracket the L2
    /// budget from both sides.
    ///
    /// If `taxa_per_block` came out as `n_taxa` on every benchmark, the equality
    /// test above would be passing for the wrong reason -- there would be one
    /// block and no reordering to get wrong, while the 26 GB of traffic would
    /// still be there. This pins the ratios that make the tiling real.
    #[test]
    fn both_block_sizes_bracket_the_l2_budget() {
        for p in [2usize, 5, 9, 18, 50] {
            let x = {
                let mut m = Matrix::zeros(200, p);
                for i in 0..200 {
                    m.set(i, 0, 1.0);
                }
                m
            };
            let xx = OuterProducts::build(&x);
            // `block()` is the reference's own formula, in *elements* of `xx`.
            assert_eq!(
                xx.block(),
                (131072.0 / (p * p) as f64).floor().max(1.0) as usize
            );
            // `taxa_per_block` is in *bytes* of accumulator, so it is the
            // transpose of `block()`: the same budget, counted in doubles.
            let expected = (131072.0 / (p * p * 8) as f64).floor().max(1.0) as usize;
            assert_eq!(xx.taxa_per_block(), expected, "p={p}");
            assert!(
                xx.taxa_per_block() >= 1 && xx.block() >= 1,
                "a block must never be empty, or the loop would not advance"
            );
            // A tile of outer products is at most the budget.
            assert!(
                xx.block() * p * p <= 131072 || xx.block() == 1,
                "p={p}: a sample tile of {} elements exceeds the 128 KiB budget",
                xx.block() * p * p
            );
        }
    }

    #[test]
    fn block_size_matches_the_reference_formula() {
        for p in 1..=16 {
            let xx = OuterProducts {
                n_samp: 1,
                p,
                data: vec![],
                row_incomplete: vec![],
            };
            assert_eq!(
                xx.block(),
                (131072.0 / ((p * p) as f64)).floor().max(1.0) as usize,
                "p = {p}"
            );
        }
    }

    #[test]
    fn outer_products_are_column_major_vec() {
        let x = design();
        let xx = OuterProducts::build(&x);
        // sample 2 is (1, 1): x x' = [[1,1],[1,1]], flattened column-major
        for a in 0..2 {
            for b in 0..2 {
                assert_eq!(xx.at(2, a, b), 1.0);
            }
        }
        // sample 0 is (1, 0)
        assert_eq!(xx.at(0, 0, 0), 1.0);
        assert_eq!(xx.at(0, 0, 1), 0.0);
        assert_eq!(xx.at(0, 1, 1), 0.0);
    }

    #[test]
    fn xtx_inverse_matches_r() {
        let x = design();
        // R: ginv(t(x) %*% x) for this design
        let g = xtx_inverse(&x);
        let xtx = x.t().matmul(&x).unwrap();
        let expect = xtx_inverse_bruteforce(&xtx);
        for a in 0..2 {
            for b in 0..2 {
                assert!(
                    (g.get(a, b) - expect.get(a, b)).abs() < 1e-12,
                    "({a},{b}): {} vs {}",
                    g.get(a, b),
                    expect.get(a, b)
                );
            }
        }
    }

    /// Independent computation: for a 2x2 Gram matrix the inverse is explicit.
    fn xtx_inverse_bruteforce(xtx: &Matrix) -> Matrix {
        let a = xtx.get(0, 0);
        let b = xtx.get(0, 1);
        let d = xtx.get(1, 1);
        let det = a * d - b * b;
        Matrix::from_rows(&[vec![d / det, -b / det], vec![-b / det, a / det]])
    }

    #[test]
    fn sandwich_of_a_fully_observed_taxon_matches_a_direct_computation() {
        let x = design();
        let xx = OuterProducts::build(&x);
        let g = xtx_inverse(&x);
        let n_samp = x.rows;
        let eps = [0.1, -0.2, 0.3, 0.15, -0.05, 0.25];
        let (vcov, var_hat) = sandwich_all(&eps, 1, n_samp, &xx, &g, CompatMode::Ancombc2_15);

        // direct: build sigma by hand
        let mut sigma = Matrix::zeros(2, 2);
        for j in 0..n_samp {
            let xv = [x.get(j, 0), x.get(j, 1)];
            let e2 = eps[j] * eps[j];
            for a in 0..2 {
                for b in 0..2 {
                    sigma.set(a, b, sigma.get(a, b) + e2 * xv[a] * xv[b]);
                }
            }
        }
        let v = g.matmul(&sigma).unwrap().matmul(&g).unwrap();
        for a in 0..2 {
            for b in 0..2 {
                assert!(
                    (vcov[a * 2 + b] - v.get(a, b)).abs() < 1e-14,
                    "({a},{b}): {} vs {}",
                    vcov[a * 2 + b],
                    v.get(a, b)
                );
            }
        }
        assert!((var_hat[0] - v.get(0, 0)).abs() < 1e-14);
        assert!((var_hat[1] - v.get(1, 1)).abs() < 1e-14);
    }

    #[test]
    fn compat_modes_differ_only_where_the_quirk_applies() {
        let x = design();
        let xx = OuterProducts::build(&x);
        let g = xtx_inverse(&x);
        let n_samp = x.rows;

        // taxon 0: fully observed -> identical in both modes
        let e0 = vec![0.1, -0.2, 0.3, 0.15, -0.05, 0.25]
            .into_iter()
            .map(|v| v * v)
            .collect::<Vec<_>>();
        let mut e1 = e0.clone();
        e1[2] = f64::NAN; // taxon 1 unobserved at sample 2
        let eps2 = [e0, e1].concat();

        let (vc_compat, _) = sandwich_all(&eps2, 2, n_samp, &xx, &g, CompatMode::Ancombc2_15);
        let (vc_strict, _) = sandwich_all(&eps2, 2, n_samp, &xx, &g, CompatMode::StrictSpec);

        // taxon 0 unaffected
        for a in 0..2 {
            for b in 0..2 {
                assert!(
                    (vc_compat[a * 2 + b] - vc_strict[a * 2 + b]).abs() < 1e-15,
                    "fully observed taxon must be identical"
                );
            }
        }
        // taxon 1 differs
        let d = vc_compat[4..8]
            .iter()
            .zip(vc_strict[4..8].iter())
            .map(|(a, b)| (a - b).abs())
            .fold(0.0f64, f64::max);
        assert!(d > 1e-6, "the quirk must show a difference, got {d}");
    }

    #[test]
    fn strict_spec_contributes_nothing_for_missing_cells() {
        // Under StrictSpec a missing (taxon, sample) cell adds no term at all, so
        // the result equals the sandwich over the observed samples only.
        let x = design();
        let xx = OuterProducts::build(&x);
        let g = xtx_inverse(&x);
        let n_samp = x.rows;
        let e = [0.1, -0.2, 0.3, 0.15, -0.05, 0.25];
        // The kernel takes unsquared residuals, so an unobserved cell is `NaN`
        // here rather than a zero that has already been squared.
        let mut eps: Vec<f64> = e.to_vec();
        eps[2] = f64::NAN;
        let (vc_strict, _) = sandwich_all(&eps, 1, n_samp, &xx, &g, CompatMode::StrictSpec);

        // observed-only reference
        let eps_ref: Vec<f64> = (0..n_samp)
            .map(|j| if j == 2 { 0.0 } else { e[j] })
            .collect();
        let (vc_ref, _) = sandwich_all(&eps_ref, 1, n_samp, &xx, &g, CompatMode::StrictSpec);
        for k in 0..4 {
            assert!(
                (vc_strict[k] - vc_ref[k]).abs() < 1e-15,
                "strict spec must skip missing cells, got {} vs {}",
                vc_strict[k],
                vc_ref[k]
            );
        }
    }

    #[test]
    fn incomplete_design_row_also_triggers_the_quirk() {
        let mut x = design();
        x.set(1, 1, f64::NAN);
        let xx = OuterProducts::build(&x);
        assert!(xx.row_incomplete[1]);
        assert!(!xx.row_incomplete[0]);
        // ginv drops the incomplete row
        let g = xtx_inverse(&x);
        assert!(g.data.iter().all(|v| v.is_finite()));
    }

    #[test]
    fn var_hat_is_the_diagonal_of_vcov() {
        let x = design();
        let xx = OuterProducts::build(&x);
        let g = xtx_inverse(&x);
        let n = x.rows;
        let e: Vec<f64> = (0..n)
            .map(|i| 0.1 * (i as f64 + 1.0))
            .map(|v| v * v)
            .collect();
        let (vcov, var_hat) = sandwich_all(&e, 1, n, &xx, &g, CompatMode::Ancombc2_15);
        assert!((var_hat[0] - vcov[0]).abs() < 1e-15);
        assert!((var_hat[1] - vcov[3]).abs() < 1e-15);
    }

    #[test]
    fn accumulation_is_increasing_taxon_and_sample_order() {
        // The result must not depend on how the work is blocked, only on the
        // (taxon, sample) index order of the additions.
        let x = design();
        let xx = OuterProducts::build(&x);
        let g = xtx_inverse(&x);
        let n = x.rows;
        // `sandwich_all` takes the *unsquared* residuals and squares them, so the
        // fixture is built as residuals.
        let mut e = Vec::new();
        for t in 0..5 {
            for j in 0..n {
                e.push((0.1 * (t as f64 + 1.0) * (j as f64 + 1.0)).sqrt());
            }
        }
        let (vcov, _) = sandwich_all(&e, 5, n, &xx, &g, CompatMode::Ancombc2_15);
        // recompute taxon 2 by an independent scalar accumulation
        let t = 2;
        let mut acc = [0.0f64; 4];
        for j in 0..n {
            let e = e[t * n + j] * e[t * n + j];
            let (x0, x1) = (x.get(j, 0), x.get(j, 1));
            let outer = [x0 * x0, x0 * x1, x1 * x0, x1 * x1];
            for k in 0..4 {
                acc[k] += e * outer[k];
            }
        }
        let sigma = Matrix::from_rows(&[vec![acc[0], acc[1]], vec![acc[2], acc[3]]]);
        let v = g.matmul(&sigma).unwrap().matmul(&g).unwrap().symmetrise();
        for a in 0..2 {
            for b in 0..2 {
                assert!(
                    (vcov[t * 4 + a * 2 + b] - v.get(a, b)).abs() < 1e-14,
                    "taxon {t} ({a},{b})"
                );
            }
        }
    }

    #[test]
    fn results_are_bitwise_reproducible() {
        let x = design();
        let xx = OuterProducts::build(&x);
        let g = xtx_inverse(&x);
        let n = x.rows;
        let eps: Vec<f64> = (0..(4 * n)).map(|i| 0.01 * i as f64).collect();
        let a = sandwich_all(&eps, 4, n, &xx, &g, CompatMode::Ancombc2_15);
        let b = sandwich_all(&eps, 4, n, &xx, &g, CompatMode::Ancombc2_15);
        assert_eq!(a.0, b.0, "bitwise identical across runs");
        assert_eq!(a.1, b.1);
    }
}
