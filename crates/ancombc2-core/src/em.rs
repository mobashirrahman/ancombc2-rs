//! The E-M bias estimator: `.bias_em`.
//!
//! ANCOM-BC2 models a coefficient across all taxa as a three-component Gaussian
//! mixture — negative-effect taxa, null taxa, positive-effect taxa — and
//! estimates the sampling-fraction bias `delta` as the mixture's central
//! location. Per component `k` with prior `pi_k`,
//!
//! ```text
//! pdf0 = dnorm(beta,           delta,     sqrt(nu0))
//! pdf1 = dnorm(beta, delta + l1, sqrt(nu0 + kappa1))
//! pdf2 = dnorm(beta, delta + l2, sqrt(nu0 + kappa2))
//! r_ki = pi_k pdf_k / sum_j pi_j pdf_j
//! ```
//!
//! The M-step is closed-form for the weights, the locations and the variance
//! offsets, **except** for `kappa1` and `kappa2`: the reference optimises those
//! two with Nelder-Mead at every E-M iteration, via
//! `nloptr::nloptr(algorithm = "NLOPT_LN_NELDERMEAD", lb = 0)`. That is the
//! single largest source of divergence risk in the whole port, so it is
//! reproduced rather than replaced (see `ancombc2_stats::nelder_mead`).
//!
//! Two details are easy to get wrong and both are reproduced here:
//!
//! * `l1` is clamped with `min(., 0)` and `l2` with `max(., 0)`, so component 1
//!   can only sit at or below the centre and component 2 at or above it.
//! * the variance used per component is `nu0 + kappa_k`, with `nu0` the taxon's
//!   sandwich variance; the WLS estimator adds `var_delta = 1 / sum(1/nu)` to
//!   that later.
//!
//! The initial values are quantiles of `beta`: the centre from the interquartile
//! range, the component locations from the 12.5th and 87.5th percentiles, and
//! the variance offsets from the corresponding sample variances, each with a
//! documented fallback.

use ancombc2_stats::{dnorm, nelder_mead_scalar, quantile_type7_unsorted, var_unbiased};

use crate::error::{AncombcError, Result};

/// Priors fixed by the reference.
const PI0_0: f64 = 0.75;
const PI1_0: f64 = 0.125;
const PI2_0: f64 = 0.125;

/// The estimated bias and the variance of that estimate.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BiasResult {
    /// E-M estimate of the sample-specific bias.
    pub delta_em: f64,
    /// Weighted-least-squares estimate of the same quantity.
    pub delta_wls: f64,
    /// `1 / sum(1/nu)`, the variance of the bias estimate.
    pub var_delta: f64,
    /// Final mixture parameters, for the golden contract.
    pub params: MixtureParams,
    pub iterations: usize,
}

/// The fitted three-component mixture.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct MixtureParams {
    pub pi0: f64,
    pub pi1: f64,
    pub pi2: f64,
    pub delta: f64,
    pub l1: f64,
    pub l2: f64,
    pub kappa1: f64,
    pub kappa2: f64,
}

impl MixtureParams {
    pub fn initial() -> Self {
        Self {
            pi0: PI0_0,
            pi1: PI1_0,
            pi2: PI2_0,
            delta: f64::NAN,
            l1: f64::NAN,
            l2: f64::NAN,
            kappa1: f64::NAN,
            kappa2: f64::NAN,
        }
    }
}

/// `.bias_em` for one fixed effect.
pub fn bias_em(beta: &[f64], var_hat: &[f64], tol: f64, max_iter: usize) -> Result<BiasResult> {
    // Drop taxa with a missing coefficient or variance; R filters them with
    // `neither_na`, keeping the *paired* entries.
    let mut b: Vec<f64> = Vec::with_capacity(beta.len());
    let mut nu0: Vec<f64> = Vec::with_capacity(beta.len());
    let mut zero_var_taxa: Vec<usize> = Vec::new();
    for i in 0..beta.len() {
        if beta[i].is_nan() || var_hat[i].is_nan() {
            continue;
        }
        if var_hat[i] == 0.0 {
            zero_var_taxa.push(i);
        }
        b.push(beta[i]);
        nu0.push(var_hat[i]);
    }
    if !zero_var_taxa.is_empty() {
        return Err(AncombcError::ZeroVariance {
            taxa: zero_var_taxa
                .iter()
                .map(|i| format!("taxon_{i}"))
                .collect::<Vec<_>>()
                .join(", "),
        });
    }
    if b.is_empty() {
        return Ok(BiasResult {
            delta_em: f64::NAN,
            delta_wls: f64::NAN,
            var_delta: 0.0,
            params: MixtureParams::initial(),
            iterations: 0,
        });
    }

    let init = initial_params(&b);
    let mut cur = init;
    let mut iterations = 0usize;

    while iterations < max_iter {
        // ---- E-step ----
        let (r0, r1, r2) = e_step(&b, &nu0, &cur);
        // ---- M-step ----
        let next = m_step(&b, &nu0, &r0, &r1, &r2, &cur);

        let epsilon = (next.pi0 - cur.pi0).powi(2)
            + (next.pi1 - cur.pi1).powi(2)
            + (next.pi2 - cur.pi2).powi(2)
            + (next.delta - cur.delta).powi(2)
            + (next.l1 - cur.l1).powi(2)
            + (next.l2 - cur.l2).powi(2)
            + (next.kappa1 - cur.kappa1).powi(2)
            + (next.kappa2 - cur.kappa2).powi(2);
        iterations += 1;
        cur = next;
        if epsilon.sqrt() <= tol {
            break;
        }
    }

    let (delta_em, delta_wls, var_delta) = wls_estimators(&b, &nu0, &cur);
    Ok(BiasResult {
        delta_em,
        delta_wls,
        var_delta,
        params: cur,
        iterations,
    })
}

/// Run `bias_em` for every fixed effect, in parallel.
///
/// `beta` and `var_hat` are `n_taxa x n_eff` **row-major**, so a coefficient's
/// values are strided rather than contiguous and must be gathered before the
/// E-M can sweep over them. Gathering costs `n_eff` small buffers, which is
/// negligible next to the O(n_eff * max_iter * n_taxa) work the fit itself does.
///
/// Each coefficient is an independent problem, so this is the second level of
/// the parallel plan (inside a pseudo-count run, above the missingness groups).
/// Results are collected by index, so the thread count cannot affect them.
pub fn bias_em_all(
    beta: &[f64],
    var_hat: &[f64],
    n_taxa: usize,
    n_eff: usize,
    tol: f64,
    max_iter: usize,
) -> Result<Vec<BiasResult>> {
    debug_assert_eq!(beta.len(), n_taxa * n_eff);
    debug_assert_eq!(var_hat.len(), n_taxa * n_eff);
    let mut err: Option<AncombcError> = None;
    // Level 2 of the nesting order: the E-M coefficients.
    //
    // Claiming the pool here is what stops the E-M from nesting inside a
    // pseudo-count sweep. The `p` coefficients are independent of each other --
    // each is a separate three-component mixture fitted over all the taxa -- so
    // splitting them is deterministic with no reduction to order.
    let mut lvl = crate::parallel::NestingBudget::level("E-M coefficients");
    let coefs: Vec<usize> = (0..n_eff).collect();
    let results: Vec<Result<BiasResult>> = crate::parallel::map_par(&mut lvl, &coefs, |&k| {
        let b: Vec<f64> = (0..n_taxa).map(|i| beta[i * n_eff + k]).collect();
        let v: Vec<f64> = (0..n_taxa).map(|i| var_hat[i * n_eff + k]).collect();
        bias_em(&b, &v, tol, max_iter)
    });
    let mut out = Vec::with_capacity(n_eff);
    for r in results {
        match r {
            Ok(x) => out.push(x),
            Err(e) => {
                if err.is_none() {
                    err = Some(e);
                }
            }
        }
    }
    match err {
        Some(e) => Err(e),
        None => Ok(out),
    }
}

/// Exposed for the EM trajectory test, which pins each M-step against the
/// reference's printed trace. Not part of the public API surface.
#[doc(hidden)]
pub fn m_step_public(
    beta: &[f64],
    nu0: &[f64],
    r0: &[f64],
    r1: &[f64],
    r2: &[f64],
    prev: &MixtureParams,
) -> MixtureParams {
    m_step(beta, nu0, r0, r1, r2, prev)
}

/// The quantile-based starting values, with the reference's fallbacks.
fn initial_params(beta: &[f64]) -> MixtureParams {
    let q = |p: f64| quantile_type7_unsorted(beta, p);
    let q25 = q(0.25);
    let q75 = q(0.75);
    let q125 = q(0.125);
    let q875 = q(0.875);

    // delta_0: mean of beta inside [q25, q75]
    let mid: Vec<f64> = beta
        .iter()
        .copied()
        .filter(|v| *v >= q25 && *v <= q75)
        .collect();
    let mut delta = mean(&mid);
    if delta.is_nan() {
        delta = mean(beta);
    }

    // l1_0: mean of beta below q125, falling back to min(beta)
    let lo: Vec<f64> = beta.iter().copied().filter(|v| *v < q125).collect();
    let mut l1 = mean(&lo);
    if l1.is_nan() {
        l1 = beta.iter().copied().fold(f64::INFINITY, f64::min);
    }

    // l2_0: mean of beta above q875, falling back to max(beta)
    let hi: Vec<f64> = beta.iter().copied().filter(|v| *v > q875).collect();
    let mut l2 = mean(&hi);
    if l2.is_nan() {
        l2 = beta.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    }

    // kappa: the sample variances of those same tails, with 1 as the fallback
    let mut kappa1 = var_unbiased(&lo);
    if kappa1.is_nan() || kappa1 == 0.0 {
        kappa1 = 1.0;
    }
    let mut kappa2 = var_unbiased(&hi);
    if kappa2.is_nan() || kappa2 == 0.0 {
        kappa2 = 1.0;
    }

    MixtureParams {
        pi0: PI0_0,
        pi1: PI1_0,
        pi2: PI2_0,
        delta,
        l1,
        l2,
        kappa1,
        kappa2,
    }
}

fn mean(v: &[f64]) -> f64 {
    if v.is_empty() {
        return f64::NAN;
    }
    let s: f64 = v.iter().sum();
    s / v.len() as f64
}

/// Responsibilities `r_ki`. `NaN` densities become 0, as `r[is.na(r)] = 0`.
fn e_step(beta: &[f64], nu0: &[f64], p: &MixtureParams) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    let n = beta.len();
    let mut r0 = vec![0.0; n];
    let mut r1 = vec![0.0; n];
    let mut r2 = vec![0.0; n];
    for i in 0..n {
        let pdf0 = dnorm(beta[i], p.delta, nu0[i].sqrt());
        let pdf1 = dnorm(beta[i], p.delta + p.l1, (nu0[i] + p.kappa1).sqrt());
        let pdf2 = dnorm(beta[i], p.delta + p.l2, (nu0[i] + p.kappa2).sqrt());
        let denom = p.pi0 * pdf0 + p.pi1 * pdf1 + p.pi2 * pdf2;
        if denom.is_nan() || denom == 0.0 {
            continue;
        }
        r0[i] = p.pi0 * pdf0 / denom;
        r1[i] = p.pi1 * pdf1 / denom;
        r2[i] = p.pi2 * pdf2 / denom;
    }
    (r0, r1, r2)
}

fn m_step(
    beta: &[f64],
    nu0: &[f64],
    r0: &[f64],
    r1: &[f64],
    r2: &[f64],
    prev: &MixtureParams,
) -> MixtureParams {
    let n = beta.len() as f64;
    let pi0 = mean(r0);
    let pi1 = mean(r1);
    let pi2 = mean(r2);

    let mut num = 0.0;
    let mut den = 0.0;
    for i in 0..beta.len() {
        num += r0[i] * beta[i] / nu0[i]
            + r1[i] * (beta[i] - prev.l1) / (nu0[i] + prev.kappa1)
            + r2[i] * (beta[i] - prev.l2) / (nu0[i] + prev.kappa2);
        den += r0[i] / nu0[i] + r1[i] / (nu0[i] + prev.kappa1) + r2[i] / (nu0[i] + prev.kappa2);
    }
    let mut delta = num / den;
    if delta.is_nan() {
        delta = prev.delta;
    }

    // l1 is constrained to be <= 0, l2 to be >= 0: component 1 sits at or below
    // the centre, component 2 at or above.
    //
    // As with the variance objectives below, the centring term is `beta - delta`
    // with the **previous** iteration's `delta`, not the `delta_new` computed a
    // few lines above. The reference names the new value `delta_new` and never
    // reaches for it here, so the location update and the component-location
    // update are based on different centres. See docs/reference_behavior.md.
    let mut n1 = 0.0;
    let mut d1 = 0.0;
    let mut n2 = 0.0;
    let mut d2 = 0.0;
    for i in 0..beta.len() {
        n1 += r1[i] * (beta[i] - prev.delta) / (nu0[i] + prev.kappa1);
        d1 += r1[i] / (nu0[i] + prev.kappa1);
        n2 += r2[i] * (beta[i] - prev.delta) / (nu0[i] + prev.kappa2);
        d2 += r2[i] / (nu0[i] + prev.kappa2);
    }
    let mut l1 = if d1 == 0.0 { 0.0 } else { (n1 / d1).min(0.0) };
    if l1.is_nan() {
        l1 = 0.0;
    }
    let mut l2 = if d2 == 0.0 { 0.0 } else { (n2 / d2).max(0.0) };
    if l2.is_nan() {
        l2 = 0.0;
    }

    // The two variance offsets are Nelder-Mead problems. R's objective negates
    // the weighted log-likelihood and maps an infinite log-pdf to 0:
    //   -sum(r_k * log_pdf),  log_pdf[is.infinite(log_pdf)] = 0
    // Reproduced literally, because the `is.infinite` substitution changes the
    // objective on the boundary kappa = 0.
    //
    // **The location parameters in that objective are the previous iteration's,
    // not the ones just computed.** The reference closes over `delta` and `l1` /
    // `l2` — the values read at the top of the loop — while the freshly updated
    // ones are named `delta_new`, `l1_new`, `l2_new` and go unused here. Using
    // the new values gives a *different objective* and moves `kappa1` off the
    // boundary, which visibly changes `delta_em` (on the fx01 fixture, 9.5e-3
    // versus 0). It looks like a bug in the reference; it is reproduced
    // regardless, and recorded in docs/reference_behavior.md.
    let kappa1 = nelder_mead_scalar(
        |k| {
            let mut s = 0.0;
            for i in 0..beta.len() {
                let lp = dnorm(beta[i], prev.delta + prev.l1, (nu0[i] + k).sqrt()).ln();
                s += r1[i] * if lp.is_infinite() { 0.0 } else { lp };
            }
            -s
        },
        prev.kappa1,
        Some(0.0),
    );
    let kappa2 = nelder_mead_scalar(
        |k| {
            let mut s = 0.0;
            for i in 0..beta.len() {
                let lp = dnorm(beta[i], prev.delta + prev.l2, (nu0[i] + k).sqrt()).ln();
                s += r2[i] * if lp.is_infinite() { 0.0 } else { lp };
            }
            -s
        },
        prev.kappa2,
        Some(0.0),
    );

    let _ = n;
    MixtureParams {
        pi0,
        pi1,
        pi2,
        delta,
        l1,
        l2,
        kappa1,
        kappa2,
    }
}

/// The E-M and WLS bias estimates, plus the variance of the bias.
///
/// R assigns each taxon to a cluster using the *final* mixture weights as
/// quantiles, then forms the inverse-variance weighted mean. The cluster
/// assignment is what makes this a "mixed" estimator rather than a plain
/// weighted mean.
fn wls_estimators(beta: &[f64], nu0: &[f64], p: &MixtureParams) -> (f64, f64, f64) {
    let delta_em = p.delta;

    let q_lo = quantile_type7_unsorted(beta, p.pi1);
    let q_hi = quantile_type7_unsorted(beta, 1.0 - p.pi2);

    let mut nu = nu0.to_vec();
    for i in 0..beta.len() {
        if beta[i] < q_lo {
            nu[i] += p.kappa1;
        } else if beta[i] >= q_hi {
            nu[i] += p.kappa2;
        }
    }

    let mut wls_deno = 0.0;
    let mut wls_nume = 0.0;
    for i in 0..beta.len() {
        let w = 1.0 / nu[i];
        wls_deno += w;
        // R computes `wls_nume <- 1/nu` and then *overwrites* the cluster
        // members with `(wls_nume * (beta - l_k))`, i.e. the weight times the
        // component-centred value.
        let c = if beta[i] < q_lo {
            beta[i] - p.l1
        } else if beta[i] >= q_hi {
            beta[i] - p.l2
        } else {
            beta[i]
        };
        wls_nume += w * c;
    }
    let delta_wls = wls_nume / wls_deno;
    let var_delta = 1.0 / wls_deno;
    let var_delta = if var_delta.is_nan() { 0.0 } else { var_delta };
    (delta_em, delta_wls, var_delta)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn counts(n: usize) -> (Vec<f64>, Vec<f64>) {
        let beta: Vec<f64> = (0..n)
            .map(|i| {
                let t = i as f64 / n as f64;
                -1.5 + 3.0 * t + 0.05 * ((i * 7) as f64).sin()
            })
            .collect();
        let var: Vec<f64> = (0..n).map(|i| 0.2 + 0.3 * ((i % 5) as f64) / 5.0).collect();
        (beta, var)
    }

    #[test]
    fn estimates_are_finite_and_in_range() {
        let (b, v) = counts(300);
        let r = bias_em(&b, &v, 1e-5, 100).unwrap();
        assert!(r.delta_em.is_finite(), "delta_em {}", r.delta_em);
        assert!(r.delta_wls.is_finite(), "delta_wls {}", r.delta_wls);
        assert!(r.var_delta > 0.0, "var_delta {}", r.var_delta);
        // the bias is small relative to the spread of beta
        assert!(r.delta_em.abs() < 1.0, "delta_em {}", r.delta_em);
    }

    #[test]
    fn initial_priors_match_the_reference() {
        // After one iteration the weights are the mean responsibilities, so they
        // are close to but no longer exactly the starting priors. What the
        // iteration must do is move them *from* (0.75, 0.125, 0.125) towards the
        // data, and stay a probability vector.
        let (b, v) = counts(200);
        let r = bias_em(&b, &v, 1e-5, 1).unwrap();
        let s = r.params.pi0 + r.params.pi1 + r.params.pi2;
        assert!((s - 1.0).abs() < 1e-12, "weights must sum to 1, got {s}");
        assert!((r.params.pi0 - 0.75).abs() > 1e-6, "pi0 must have moved");
        assert!(r.params.pi0 > 0.5, "pi0 {}", r.params.pi0);
        assert!(r.params.pi1 > 0.0 && r.params.pi1 < 0.5);
        assert!(r.params.pi2 > 0.0 && r.params.pi2 < 0.5);
    }

    #[test]
    fn initial_params_are_the_documented_quantile_starting_values() {
        let beta = vec![-2.0, -1.0, -0.5, 0.0, 0.5, 1.0, 2.0];
        let p = initial_params(&beta);
        assert!((p.pi0 - 0.75).abs() < 1e-15);
        assert!((p.pi1 - 0.125).abs() < 1e-15);
        assert!((p.pi2 - 0.125).abs() < 1e-15);
        // the tails are small here, so the mean falls back to min / max
        assert_eq!(p.l1, -2.0, "l1 falls back to min(beta)");
        assert_eq!(p.l2, 2.0, "l2 falls back to max(beta)");
        // and the variance fallbacks are both 1
        assert_eq!(p.kappa1, 1.0);
        assert_eq!(p.kappa2, 1.0);
        // delta is the mean of the interquartile range
        let mid: Vec<f64> = beta
            .iter()
            .copied()
            .filter(|v| *v >= -0.5 && *v <= 0.5)
            .collect();
        assert!((p.delta - mean(&mid)).abs() < 1e-15);
    }

    #[test]
    fn component_locations_stay_on_their_side_of_the_centre() {
        let (b, v) = counts(400);
        let r = bias_em(&b, &v, 1e-5, 100).unwrap();
        assert!(r.params.l1 <= 1e-12, "l1 must be <= 0, got {}", r.params.l1);
        assert!(
            r.params.l2 >= -1e-12,
            "l2 must be >= 0, got {}",
            r.params.l2
        );
    }

    #[test]
    fn a_pure_null_mixture_recovers_a_near_zero_bias() {
        // beta centred on 0: the estimated bias should be close to 0.
        let b: Vec<f64> = (0..500)
            .map(|i| 0.05 * (((i * 37) % 101) as f64 - 50.0) / 50.0)
            .collect();
        let v = vec![0.25; 500];
        let r = bias_em(&b, &v, 1e-5, 100).unwrap();
        assert!(r.delta_em.abs() < 0.05, "delta_em {}", r.delta_em);
    }

    #[test]
    fn a_shifted_mixture_recovers_the_shift() {
        // beta centred on 0.7: the bias estimate should land near 0.7.
        let b: Vec<f64> = (0..500)
            .map(|i| 0.7 + 0.05 * (((i * 37) % 101) as f64 - 50.0) / 50.0)
            .collect();
        let v = vec![0.25; 500];
        let r = bias_em(&b, &v, 1e-5, 100).unwrap();
        assert!((r.delta_em - 0.7).abs() < 0.1, "delta_em {}", r.delta_em);
    }

    #[test]
    fn missing_entries_are_dropped_as_a_pair() {
        let (mut b, mut v) = counts(100);
        b[3] = f64::NAN;
        v[3] = f64::NAN;
        b[7] = f64::NAN;
        // v[7] is finite: the pair is still dropped
        let r = bias_em(&b, &v, 1e-5, 100).unwrap();
        assert!(r.delta_em.is_finite());
    }

    #[test]
    fn zero_variance_is_an_error_naming_the_taxon() {
        let (b, mut v) = counts(50);
        v[10] = 0.0;
        let e = bias_em(&b, &v, 1e-5, 100).unwrap_err();
        match e {
            AncombcError::ZeroVariance { taxa } => {
                assert!(taxa.contains("taxon_10"), "{taxa}")
            }
            other => panic!("expected ZeroVariance, got {other:?}"),
        }
    }

    #[test]
    fn all_missing_input_returns_nan_rather_than_panicking() {
        let b = vec![f64::NAN; 10];
        let v = vec![f64::NAN; 10];
        let r = bias_em(&b, &v, 1e-5, 100).unwrap();
        assert!(r.delta_em.is_nan());
        assert_eq!(r.var_delta, 0.0);
    }

    #[test]
    fn a_constant_beta_yields_a_finite_estimate() {
        // Every taxon identical: the quantiles collapse and the fallbacks fire.
        let b = vec![0.3; 200];
        let v = vec![0.1; 200];
        let r = bias_em(&b, &v, 1e-5, 100).unwrap();
        assert!(r.delta_em.is_finite(), "delta_em {}", r.delta_em);
    }

    #[test]
    fn convergence_tolerance_is_respected() {
        let (b, v) = counts(300);
        let loose = bias_em(&b, &v, 1e-2, 100).unwrap();
        let tight = bias_em(&b, &v, 1e-8, 200).unwrap();
        assert!(
            tight.iterations >= loose.iterations,
            "tighter tol must iterate at least as long"
        );
    }

    #[test]
    fn max_iter_is_respected() {
        let (b, v) = counts(300);
        let r = bias_em(&b, &v, 1e-12, 3).unwrap();
        assert!(r.iterations <= 3, "iterations {}", r.iterations);
    }

    #[test]
    fn all_results_are_deterministic() {
        let (b, v) = counts(257);
        let a = bias_em(&b, &v, 1e-5, 100).unwrap();
        let c = bias_em(&b, &v, 1e-5, 100).unwrap();
        assert_eq!(a.delta_em, c.delta_em, "bitwise identical");
        assert_eq!(a.delta_wls, c.delta_wls);
        assert_eq!(a.var_delta, c.var_delta);
        assert_eq!(a.params, c.params);
    }

    #[test]
    fn parallel_and_sequential_agree_bitwise() {
        // Two distinct coefficients over 400 taxa, laid out row-major: entry
        // (i, k) is at i * n_eff + k, so each column is strided.
        let n = 400usize;
        let n_eff = 2usize;
        let mut beta = vec![0.0; n * n_eff];
        let mut var = vec![0.0; n * n_eff];
        for i in 0..n {
            let t = i as f64 / n as f64;
            beta[i * n_eff] = -1.5 + 3.0 * t + 0.05 * ((i * 7) as f64).sin();
            beta[i * n_eff + 1] = 0.3 - 1.0 * t;
            var[i * n_eff] = 0.2 + 0.3 * ((i % 5) as f64) / 5.0;
            var[i * n_eff + 1] = 0.5;
        }
        let par = bias_em_all(&beta, &var, n, n_eff, 1e-5, 100).unwrap();
        let seq: Vec<_> = (0..n_eff)
            .map(|k| {
                let b: Vec<f64> = (0..n).map(|i| beta[i * n_eff + k]).collect();
                let v: Vec<f64> = (0..n).map(|i| var[i * n_eff + k]).collect();
                bias_em(&b, &v, 1e-5, 100).unwrap()
            })
            .collect();
        for k in 0..n_eff {
            assert_eq!(par[k].delta_em, seq[k].delta_em, "coefficient {k}");
        }
    }

    #[test]
    fn each_coefficient_gets_its_own_column() {
        // The regression guard for the row-major layout: entry (i, k) lives at
        // `i * n_eff + k`, so slicing contiguously would feed coefficient 1 the
        // wrong values. Each parallel result must equal a direct call on the
        // strided column.
        let n = 300usize;
        let n_eff = 3usize;
        let mut beta = vec![0.0; n * n_eff];
        let mut var = vec![0.0; n * n_eff];
        for i in 0..n {
            let t = i as f64 / n as f64;
            for k in 0..n_eff {
                beta[i * n_eff + k] = (k as f64 + 1.0) * t - (k as f64) * 0.3;
                var[i * n_eff + k] = 0.2 + 0.1 * k as f64;
            }
        }
        let par = bias_em_all(&beta, &var, n, n_eff, 1e-5, 100).unwrap();
        for k in 0..n_eff {
            let b: Vec<f64> = (0..n).map(|i| beta[i * n_eff + k]).collect();
            let v: Vec<f64> = (0..n).map(|i| var[i * n_eff + k]).collect();
            let direct = bias_em(&b, &v, 1e-5, 100).unwrap();
            assert_eq!(
                par[k].delta_em, direct.delta_em,
                "coefficient {k} did not get its own column"
            );
        }
    }

    #[test]
    fn wls_var_delta_is_exactly_one_over_the_weight_sum() {
        let (b, v) = counts(120);
        let r = bias_em(&b, &v, 1e-5, 100).unwrap();
        // Recompute the weighted sum independently from the final parameters.
        let p = r.params;
        let q_lo = quantile_type7_unsorted(&b, p.pi1);
        let q_hi = quantile_type7_unsorted(&b, 1.0 - p.pi2);
        let mut denom = 0.0;
        for i in 0..b.len() {
            let mut nu = v[i];
            if b[i] < q_lo {
                nu += p.kappa1;
            } else if b[i] >= q_hi {
                nu += p.kappa2;
            }
            denom += 1.0 / nu;
        }
        assert!(
            (r.var_delta - 1.0 / denom).abs() < 1e-15,
            "var_delta {} vs {}",
            r.var_delta,
            1.0 / denom
        );
        // and the magnitude is sensible: 1 / sum(1/nu) for 120 taxa of var ~0.4
        assert!(r.var_delta > 0.0 && r.var_delta < 0.05, "{}", r.var_delta);
    }
}
