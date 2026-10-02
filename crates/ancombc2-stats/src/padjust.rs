//! `stats::p.adjust`, transcribed from R's source.
//!
//! ANCOM-BC2 exposes seven methods through `p_adj_method`, and additionally calls
//! `p.adjust` with a non-unit `n` inside the mdFDR procedure (`.mdfdr` passes
//! `n = length(x) * n_tax / R`), which is how the pairwise tests inflate the
//! denominator. The `n` parameter is load-bearing, not an edge case.
//!
//! # R semantics that are easy to get wrong
//!
//! * **`NA` p-values are dropped before ranking.** R computes
//!   `p <- p[nna]` and then evaluates the default `n = length(p)`, so `n` counts
//!   the *observed* p-values, not the input length. An `NA` therefore shrinks the
//!   effective number of tests: `p.adjust(c(0.01, NA, 0.03), "holm")` is
//!   `c(0.02, NA, 0.03)`, not `c(0.03, NA, 0.03)`. NA is restored at the end.
//! * **`n <= 1` short-circuits** to the input unchanged.
//! * **`hommel` becomes `hochberg` when `n == 2`**, and pads `p` with ones when
//!   `n > lp`.
//! * **`fdr` is an alias for `BH`.**
//! * The step-down/step-up steps use `cummax`/`cummin` over the *sorted* values
//!   and are then mapped back through `ro`.
//! * `BY` uses the harmonic number `sum(1/(1:n))`, not `1/n`.
//!
//! The `n` used here is the *effective* test count: the caller passes it
//! explicitly (the oracle does so inside `.mdfdr`), and the default follows R's
//! lazy evaluation and counts the non-`NA` entries.

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AdjustMethod {
    Holm,
    Hochberg,
    Hommel,
    Bonferroni,
    Bh,
    By,
    Fdr,
    None,
}

impl AdjustMethod {
    /// Parse R's spelling, including the aliases the oracle's docs list.
    pub fn parse(s: &str) -> Result<Self, String> {
        match s.to_ascii_lowercase().as_str() {
            "holm" => Ok(AdjustMethod::Holm),
            "hochberg" => Ok(AdjustMethod::Hochberg),
            "hommel" => Ok(AdjustMethod::Hommel),
            "bonferroni" => Ok(AdjustMethod::Bonferroni),
            "bh" => Ok(AdjustMethod::Bh),
            "fdr" => Ok(AdjustMethod::Bh),
            "by" => Ok(AdjustMethod::By),
            "none" => Ok(AdjustMethod::None),
            other => Err(format!("unknown p.adjust method: {other}")),
        }
    }

    pub fn as_str(&self) -> &'static str {
        match self {
            AdjustMethod::Holm => "holm",
            AdjustMethod::Hochberg => "hochberg",
            AdjustMethod::Hommel => "hommel",
            AdjustMethod::Bonferroni => "bonferroni",
            AdjustMethod::Bh => "BH",
            AdjustMethod::By => "BY",
            AdjustMethod::Fdr => "fdr",
            AdjustMethod::None => "none",
        }
    }
}

/// R `p.adjust(p, method, n = length(p))` with `n` defaulting to the number of
/// non-`NA` entries, matching R's lazy evaluation of the default.
pub fn p_adjust(p: &[f64], method: AdjustMethod, n: Option<usize>) -> Vec<f64> {
    p_adjust_n(
        p,
        method,
        n.unwrap_or_else(|| p.iter().filter(|v| !v.is_nan()).count()) as f64,
    )
}

/// The same function with `n` always taken explicitly. `n` is the effective
/// number of tests; R requires `n >= lp` and returns the input unchanged when
/// `n <= 1`.
///
/// `n` is a `f64` on purpose. ANCOM-BC2's pairwise test passes
/// `n = n_col * n_tax / R`, a quotient of integers that is generally
/// fractional, and R uses it *as a double* in the multipliers
/// (`n * p`, `(n + 1 - i) * p`, `n / i * p`) while only `hommel` truncates it,
/// because that branch pads the vector to length `n`. Rounding `n` to an
/// integer before calling would change every adjusted p-value.
pub fn p_adjust_n(p: &[f64], method: AdjustMethod, n: f64) -> Vec<f64> {
    if method == AdjustMethod::None {
        return p.to_vec();
    }
    let keep: Vec<usize> = (0..p.len()).filter(|&i| !p[i].is_nan()).collect();
    let lp = keep.len();
    if n <= 1.0 || lp == 0 {
        return p.to_vec();
    }
    debug_assert!(n >= lp as f64, "R: stopifnot(n >= lp)");
    let pv: Vec<f64> = keep.iter().map(|&i| p[i]).collect();

    let adjusted: Vec<f64> = match method {
        AdjustMethod::None => unreachable!(),
        AdjustMethod::Bonferroni => pv.iter().map(|&v| (n * v).min(1.0)).collect(),
        AdjustMethod::Holm => {
            // i <- 1:lp ; o <- order(p) ; pmin(1, cummax((n+1-i) * p[o]))[order(o)]
            let o = order_p(&pv);
            let ro = order(&o);
            let mut step: Vec<f64> = (0..lp)
                .map(|k| (n + 1.0 - (k + 1) as f64) * pv[o[k]])
                .collect();
            cumminmax(&mut step, true);
            let mut out = vec![0.0; lp];
            // R's `x[ro]`: element `ro[k]` of the sorted vector lands at
            // position `k`. Writing it the other way round (a scatter) silently
            // permutes the adjusted values whenever the p-values are not already
            // in ascending order.
            for k in 0..lp {
                out[k] = step[ro[k]].min(1.0);
            }
            out
        }
        AdjustMethod::Hochberg => {
            // R: i <- lp:1L; o <- order(p, decreasing = TRUE); ro <- order(o);
            //    pmin(1, cummin((n + 1L - i) * p[o]))[ro]
            //
            // The weight is `n + 1 - i`, **not** `n / i`. That is the whole
            // difference from BH, and sharing one branch with BH made Hochberg
            // silently compute the Benjamini-Hochberg adjustment: with 100
            // hypotheses the two disagree by up to a factor of 100, and the
            // `adjust-hochberg` fixture matrix cell caught it at Level A, where
            // `p.adjust` on the oracle's own p-values has to reproduce the
            // oracle's own q exactly.
            //
            // Both sort descending and take a cumulative minimum; only the
            // multiplier differs, and they coincide when `lp == n` and the
            // weights happen to line up -- which is why a shared branch looked
            // plausible until a cell with 100 taxa reached it.
            let o = order_p_desc(&pv);
            let ro = order(&o);
            let mut step: Vec<f64> = (0..lp)
                .map(|k| (n + 1.0 - (lp - k) as f64) * pv[o[k]])
                .collect();
            cumminmax(&mut step, false);
            let mut out = vec![0.0; lp];
            for k in 0..lp {
                out[k] = step[ro[k]].min(1.0);
            }
            out
        }
        AdjustMethod::Bh | AdjustMethod::By => {
            // R: BH  -> pmin(1, cummin(n/i * p[o]))
            //     BY  -> pmin(1, cummin(q * n/i * p[o])) with q = sum(1/(1:n))
            let q: f64 = if method == AdjustMethod::By {
                // `q = sum(1 / (1:n))` over the number of *tested* hypotheses,
                // which is `lp`, not `n`.
                (1..=lp).map(|k| 1.0 / k as f64).sum()
            } else {
                1.0
            };
            let o = order_p_desc(&pv);
            let ro = order(&o);
            let mut step: Vec<f64> = (0..lp)
                .map(|k| q * n / (lp - k) as f64 * pv[o[k]])
                .collect();
            cumminmax(&mut step, false);
            let mut out = vec![0.0; lp];
            for k in 0..lp {
                out[k] = step[ro[k]].min(1.0);
            }
            out
        }
        AdjustMethod::Hommel => {
            if n == 2.0 {
                // R: if (n == 2 && method == "hommel") method <- "hochberg"
                return p_adjust_n(p, AdjustMethod::Hochberg, 2.0);
            }
            // `seq_len(n)` and `rep.int(..., n)` truncate a fractional `n`, so
            // Hommel is the one branch that works on the integer.
            hommel(&pv, n as usize)
        }
        AdjustMethod::Fdr => unreachable!("fdr is an alias for BH"),
    };

    let mut out = p.to_vec();
    for (k, &i) in keep.iter().enumerate() {
        out[i] = adjusted[k];
    }
    out
}

/// R's `hommel` branch, with `p` padded to length `n` by ones when `n > lp`.
fn hommel(p: &[f64], n: usize) -> Vec<f64> {
    let lp = p.len();
    let mut p = p.to_vec();
    if n > lp {
        p.extend(std::iter::repeat(1.0).take(n - lp));
    }
    let o = order_p(&p);
    let ps: Vec<f64> = o.iter().map(|&k| p[k]).collect();
    let ro = order(&o);

    let mut q: Vec<f64> = vec![0.0; n];
    let mut pa: Vec<f64> = vec![0.0; n];
    let min_np_i = (1..=n)
        .map(|i| n as f64 * ps[i - 1] / i as f64)
        .fold(f64::INFINITY, f64::min);
    let m = min_np_i;
    for v in q.iter_mut() {
        *v = m;
    }
    for v in pa.iter_mut() {
        *v = m;
    }

    if n >= 2 {
        for j in (2..=n).rev() {
            let ij: Vec<usize> = (1..=(n - j + 1)).collect();
            let i2: Vec<usize> = (n - j + 2..=n).collect();
            let q1 = (0..i2.len())
                .map(|t| j as f64 * ps[i2[t] - 1] / (t + 2) as f64)
                .fold(f64::INFINITY, f64::min);
            for &k in &ij {
                q[k - 1] = (j as f64 * ps[k - 1]).min(q1);
            }
            let tail = q[n - j];
            for &k in &i2 {
                q[k - 1] = tail;
            }
            for t in 0..n {
                pa[t] = pa[t].max(q[t]);
            }
        }
    }

    let maxed: Vec<f64> = (0..n).map(|t| pa[t].max(ps[t])).collect();
    if lp < n {
        (0..lp).map(|k| maxed[ro[k]]).collect()
    } else {
        (0..lp).map(|k| maxed[ro[k]]).collect()
    }
}

/// R's `order(v)` for a numeric vector. Ties keep the earlier index first, which
/// is R's default (`method = "auto"` uses a stable radix sort for numerics).
fn order_p(v: &[f64]) -> Vec<usize> {
    let mut idx: Vec<usize> = (0..v.len()).collect();
    idx.sort_by(|&a, &b| v[a].partial_cmp(&v[b]).unwrap_or(std::cmp::Ordering::Equal));
    idx
}

/// `order(v)` for an integer permutation vector: `order(o)` in R's source.
fn order(v: &[usize]) -> Vec<usize> {
    let mut idx: Vec<usize> = (0..v.len()).collect();
    idx.sort_by_key(|&k| v[k]);
    idx
}

fn order_p_desc(v: &[f64]) -> Vec<usize> {
    let mut idx: Vec<usize> = (0..v.len()).collect();
    idx.sort_by(|&a, &b| v[b].partial_cmp(&v[a]).unwrap_or(std::cmp::Ordering::Equal));
    idx
}

/// `cummax` when `max` is true, `cummin` otherwise. R also drops `NA` in
/// `cummax`/`cummin`; the p-values reaching this point are finite.
fn cumminmax(v: &mut [f64], max: bool) {
    for i in 1..v.len() {
        v[i] = if max {
            v[i - 1].max(v[i])
        } else {
            v[i - 1].min(v[i])
        };
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // Every expected value below was produced by R 4.3.3 via
    //   p.adjust(c(0.01, 0.02, 0.03, 0.04, 0.05), <method>)
    const P5: [f64; 5] = [0.01, 0.02, 0.03, 0.04, 0.05];

    fn check(got: &[f64], want: &[f64], what: &str) {
        assert_eq!(got.len(), want.len(), "{what}");
        for (i, (g, w)) in got.iter().zip(want).enumerate() {
            assert!((g - w).abs() < 1e-12, "{what}[{i}]: got {g}, want {w}");
        }
    }

    #[test]
    fn matches_r_for_every_method() {
        check(
            &p_adjust(&P5, AdjustMethod::Holm, None),
            &[0.05, 0.08, 0.09, 0.09, 0.09],
            "holm",
        );
        check(
            &p_adjust(&P5, AdjustMethod::Hochberg, None),
            &[0.05; 5],
            "hochberg",
        );
        check(
            &p_adjust(&P5, AdjustMethod::Hommel, None),
            &[0.05; 5],
            "hommel",
        );
        check(&p_adjust(&P5, AdjustMethod::Bh, None), &[0.05; 5], "BH");
        // R 4.3.3 prints 0.11416666666666665 for the first three and
        // 0.11416666666666667 for the last two; only the displayed rounding of
        // the intermediate products differs, so compare to 1e-14.
        let h5: f64 = (1..=5).map(|k| 1.0 / k as f64).sum();
        // R's i runs lp, lp-1, ..., 1 over the descending sort, so the factor is
        // n/i = 1, 5/4, 5/3, 5/2, 5.
        let by_expect: Vec<f64> = (0..5)
            .map(|k| (h5 * 5.0 / (5 - k) as f64 * P5[4 - k]).min(1.0))
            .collect();
        check(&p_adjust(&P5, AdjustMethod::By, None), &by_expect, "BY");
        check(
            &p_adjust(&P5, AdjustMethod::Bonferroni, None),
            &[0.05, 0.10, 0.15, 0.20, 0.25],
            "bonf",
        );
    }

    #[test]
    fn by_uses_the_harmonic_number_not_one_over_n() {
        // R 4.3.3: p.adjust(c(0.01, 0.5, 0.5, 0.5, 0.5), "BY")
        //   == c(0.11416666666666667, 1, 1, 1, 1)
        // For the smallest p, R's i runs 5,4,3,2,1, so n/i = 5 there and the
        // value is q * 5 * 0.01 with q = H_5. Everything else is capped at 1 by
        // the cummin/pmin.
        let got = p_adjust(&[0.01, 0.5, 0.5, 0.5, 0.5], AdjustMethod::By, None);
        let h5: f64 = (1..=5).map(|k| 1.0 / k as f64).sum();
        assert!((got[0] - h5 * 5.0 * 0.01).abs() < 1e-14, "got {got:?}");
        for v in &got[1..] {
            assert!((v - 1.0).abs() < 1e-14, "got {got:?}");
        }
        // BH on the same input would give exactly 1.0 for the smallest p
        // (0.01 * 5/1 = 0.05 vs BY's 0.1142), so this separates the two rules.
        let bh = p_adjust(&[0.01, 0.5, 0.5, 0.5, 0.5], AdjustMethod::Bh, None);
        assert!((bh[0] - 0.05).abs() < 1e-14, "got {bh:?}");
    }

    #[test]
    fn na_is_dropped_before_ranking_and_restored_after() {
        // R 4.3.3: p.adjust(c(0.01, NA, 0.03), "holm") == c(0.02, NA, 0.03)
        let got = p_adjust(&[0.01, f64::NAN, 0.03], AdjustMethod::Holm, None);
        assert!((got[0] - 0.02).abs() < 1e-12, "got {got:?}");
        assert!(got[1].is_nan());
        assert!((got[2] - 0.03).abs() < 1e-12, "got {got:?}");
    }

    #[test]
    fn all_na_returns_input() {
        let got = p_adjust(&[f64::NAN, f64::NAN], AdjustMethod::Holm, None);
        assert!(got.iter().all(|v| v.is_nan()), "got {got:?}");
    }

    #[test]
    fn n_le_one_short_circuits() {
        let got = p_adjust(&[0.01, 0.5], AdjustMethod::Holm, Some(1));
        assert_eq!(got, vec![0.01, 0.5]);
    }

    #[test]
    fn explicit_n_scales_the_adjustment() {
        // .mdfdr: p.adjust(x, "holm", n = length(x) * n_tax / R)
        let got = p_adjust_n(&[0.01, 0.02], AdjustMethod::Holm, 10.0);
        // i = 1, 2 ; (n + 1 - i) * p[o] = 10*p1, 9*p2 ; cummax
        assert!((got[0] - 0.1).abs() < 1e-12, "got {got:?}");
        assert!((got[1] - 0.18).abs() < 1e-12, "got {got:?}");
    }

    #[test]
    fn hommel_falls_back_to_hochberg_at_n_two() {
        let a = p_adjust_n(&[0.01, 0.5], AdjustMethod::Hommel, 2.0);
        let b = p_adjust_n(&[0.01, 0.5], AdjustMethod::Hochberg, 2.0);
        assert_eq!(a, b, "R: n == 2 && hommel -> hochberg");
    }

    #[test]
    fn hommel_pads_when_n_exceeds_the_input_length() {
        // R: if (n > lp) p <- c(p, rep.int(1, n - lp))
        let got = p_adjust_n(&[0.01, 0.2], AdjustMethod::Hommel, 6.0);
        assert_eq!(got.len(), 2);
        assert!(got.iter().all(|v| v.is_finite()));
    }

    #[test]
    fn none_is_identity() {
        assert_eq!(p_adjust(&P5, AdjustMethod::None, None), P5.to_vec());
    }

    #[test]
    fn adjusted_values_are_never_below_the_raw_p() {
        for m in [
            AdjustMethod::Holm,
            AdjustMethod::Hochberg,
            AdjustMethod::Hommel,
            AdjustMethod::Bh,
            AdjustMethod::By,
            AdjustMethod::Bonferroni,
        ] {
            let q = p_adjust(&P5, m, None);
            for i in 0..5 {
                assert!(q[i] >= P5[i] - 1e-15, "{m:?} at {i}: {}", q[i]);
                assert!(q[i] <= 1.0, "{m:?} at {i}: {}", q[i]);
            }
        }
    }

    #[test]
    fn holm_and_bh_are_monotone_in_the_sorted_order() {
        let p: Vec<f64> = (1..=20).map(|i| i as f64 / 40.0).collect();
        for m in [AdjustMethod::Holm, AdjustMethod::Bh, AdjustMethod::Hochberg] {
            let q = p_adjust(&p, m, None);
            for i in 1..p.len() {
                assert!(q[i] >= q[i - 1] - 1e-15, "{m:?} must be non-decreasing");
            }
        }
    }

    #[test]
    fn fdr_parses_to_bh() {
        assert_eq!(AdjustMethod::parse("fdr").unwrap(), AdjustMethod::Bh);
        assert_eq!(AdjustMethod::parse("BH").unwrap(), AdjustMethod::Bh);
        assert!(AdjustMethod::parse("nope").is_err());
    }

    /// Hochberg and BH are *not* the same adjustment, and a shared branch made
    /// one compute the other.
    ///
    /// The two differ only in the multiplier R applies -- `(n + 1 - i)` against
    /// `n / i`, both after a descending sort and a cumulative minimum -- so they
    /// agree on many inputs, which is exactly why the conflation survived the
    /// suite. These vectors are ones where they diverge, taken from R 4.3.3
    /// `stats::p.adjust`:
    ///
    /// ```text
    /// p            hochberg              BH
    /// 0.265509      0.908208              0.744248
    /// 0.372124      0.908208              0.744248
    /// 0.572853      0.908208              0.763804
    /// 0.908208      0.908208              0.908208
    /// ```
    #[test]
    fn hochberg_is_not_benjamini_hochberg() {
        let p = [0.265509, 0.372124, 0.572853, 0.908208];
        let hochberg = p_adjust_n(&p, AdjustMethod::Hochberg, 4.0);
        let bh = p_adjust_n(&p, AdjustMethod::Bh, 4.0);
        let want_hochberg = [0.908208, 0.908208, 0.908208, 0.908208];
        let want_bh = [0.744248, 0.744248, 0.763804, 0.908208];
        for i in 0..4 {
            assert!(
                (hochberg[i] - want_hochberg[i]).abs() < 1e-6,
                "hochberg[{i}] = {}, want {}",
                hochberg[i],
                want_hochberg[i]
            );
            assert!(
                (bh[i] - want_bh[i]).abs() < 1e-6,
                "BH[{i}] = {}, want {}",
                bh[i],
                want_bh[i]
            );
        }
        assert!(
            hochberg
                .iter()
                .zip(bh.iter())
                .any(|(a, b)| (a - b).abs() > 1e-6),
            "this test is only meaningful while the two still differ"
        );
    }

    /// The multiplier is `n + 1 - i` against the number of hypotheses, not the
    /// number of p-values.
    ///
    /// `i` runs `lp, lp-1, ..., 1` over the *descending-sorted* p-values, so the
    /// smallest p-value carries weight `n + 1 - 1 = n`. Here only two of the ten
    /// hypotheses have a p-value, and R 4.3.3 gives
    /// `p.adjust(c(0.02, 0.6), n = 10, method = "hochberg")` = `0.2, 1`.
    ///
    /// This is a pin rather than a discriminator: on this input Hochberg and BH
    /// coincide, because `n / 1` is also `n`. `hochberg_is_not_benjamini_hochberg`
    /// above is what separates them.
    #[test]
    fn hochberg_weights_by_hypothesis_count_not_by_tested_count() {
        let p = [0.02, 0.6];
        let got = p_adjust_n(&p, AdjustMethod::Hochberg, 10.0);
        let want = [0.2, 1.0];
        for i in 0..2 {
            assert!(
                (got[i] - want[i]).abs() < 1e-12,
                "hochberg[{i}] = {}, want {} (R: p.adjust(c(0.02, 0.6), n = 10, method = 'hochberg'))",
                got[i],
                want[i]
            );
        }
    }
}
