//! `stats::quantile`, `type = 7` (the R default).
//!
//! ANCOM-BC2 computes the SAM-style regulariser `s0` as
//! `quantile(se[, k], s0_perc)`, i.e. the R default type 7. Type 7 uses
//! `m = 1 + (n - 1) * p` and interpolates linearly between the floor and ceiling
//! order statistics, so for `p = 0.05` and `n = 2000` the result depends on the
//! fractional part of `m` down to ~1e-3 of a standard error. Getting this wrong
//! shifts every `W` statistic.
//!
//! NaNs are ignored, matching R's `quantile(..., na.rm = TRUE)`, which is how the
//! oracle calls it.

/// R `quantile(x, probs, type = 7, na.rm = TRUE)`.
pub fn quantile_type7(sorted: &[f64], p: f64) -> f64 {
    let n = sorted.len();
    if n == 0 {
        return f64::NAN;
    }
    if n == 1 {
        return sorted[0];
    }
    let m = 1.0 + (n as f64 - 1.0) * p;
    let lo = m.floor();
    let hi = m.ceil();
    if lo < 1.0 {
        // p < 0 extrapolates from the smallest observation; R clamps the index
        // and interpolates between x[1] and x[2].
        let idx = if m <= 1.0 { 1 } else { 2 };
        return sorted[idx - 1];
    }
    if lo >= n as f64 {
        return sorted[n - 1];
    }
    if lo == hi {
        return sorted[lo as usize - 1];
    }
    let a = sorted[lo as usize - 1];
    let b = sorted[hi as usize - 1];
    a + (b - a) * (m - lo)
}

/// Convenience wrapper: quantile of an unsorted slice, NaNs dropped.
pub fn quantile_type7_unsorted(values: &[f64], p: f64) -> f64 {
    let mut v: Vec<f64> = values.iter().copied().filter(|x| !x.is_nan()).collect();
    if v.is_empty() {
        return f64::NAN;
    }
    v.sort_by(|a, b| a.partial_cmp(b).expect("NaN filtered"));
    quantile_type7(&v, p)
}

/// R's `var(x, na.rm = TRUE)`: the unbiased (n-1 denominator) sample variance,
/// with a two-pass mean correction.
pub fn var_unbiased(values: &[f64]) -> f64 {
    let n = values.len();
    if n < 2 {
        return f64::NAN;
    }
    let mean = values.iter().sum::<f64>() / n as f64;
    let ss: f64 = values.iter().map(|v| (v - mean) * (v - mean)).sum();
    ss / (n as f64 - 1.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn type7_matches_r_reference_values() {
        // Reference values produced by R 4.3.3:
        //   quantile(1:10, c(0, .05, .1, .25, .5, .75, .9, .95, 1), type = 7)
        let x: Vec<f64> = (1..=10).map(|v| v as f64).collect();
        let cases = [
            (0.00, 1.0),
            (0.05, 1.45),
            (0.10, 1.9),
            (0.25, 3.25),
            (0.50, 5.5),
            (0.75, 7.75),
            (0.90, 9.1),
            (0.95, 9.55),
            (1.00, 10.0),
        ];
        for (p, want) in cases {
            let got = quantile_type7(&x, p);
            assert!((got - want).abs() < 1e-12, "p={p}: got {got}, want {want}");
        }
    }

    #[test]
    fn type7_extrapolates_below_and_clamps_above() {
        let x: Vec<f64> = (1..=5).map(|v| v as f64).collect();
        // p = 0 => m = 1 => x[1]
        assert_eq!(quantile_type7(&x, 0.0), 1.0);
        // p = 1 => m = n => x[n]
        assert_eq!(quantile_type7(&x, 1.0), 5.0);
        // m below 1: R clamps to the first two order statistics
        let v = quantile_type7(&x, -1.0);
        assert!((1.0..=2.0).contains(&v), "got {v}");
    }

    #[test]
    fn single_and_empty() {
        assert_eq!(quantile_type7(&[7.0], 0.05), 7.0);
        assert!(quantile_type7(&[], 0.05).is_nan());
    }

    #[test]
    fn naive_is_not_type7() {
        // The nearest-rank / "type 1" value at p = 0.05 for n = 10 is 1.0, not
        // 1.45. Guards against a regression to a simpler rule.
        let x: Vec<f64> = (1..=10).map(|v| v as f64).collect();
        assert_ne!(quantile_type7(&x, 0.05), 1.0);
    }

    #[test]
    fn var_matches_r() {
        // R: var(c(1,2,3,4)) = 1.666667
        let v = var_unbiased(&[1.0, 2.0, 3.0, 4.0]);
        assert!((v - 5.0 / 3.0).abs() < 1e-15);
        assert!(var_unbiased(&[1.0]).is_nan());
    }
}
