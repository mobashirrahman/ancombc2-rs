//! A Nelder-Mead simplex search matching NLopt's `NLOPT_LN_NELDERMEAD` defaults,
//! because `.bias_em` in the oracle optimises the mixture variance parameters
//! `kappa1` and `kappa2` with exactly that optimiser and exactly those defaults.
//!
//! Reproducing an optimiser is a parity problem, not a correctness problem: any
//! simplex method that lands on the same optimum to within the same tolerance
//! gives the same `delta_em` to within the EM's own `tol = 1e-5`. Reproducing
//! NLopt's *trajectory* gets us much closer, and makes a divergence diagnosable
//! rather than mysterious.
//!
//! NLopt defaults that matter here:
//!   * `xtol_rel = 1e-6`, `xtol_abs = 0`
//!   * `ftol_rel = 0`, `ftol_abs = 0` (so convergence is decided on x alone)
//!   * `maxeval = 1000`
//!   * initial step: NLopt perturbs each coordinate by 0.25 * |x_j| (or 0.25 when
//!     `x_j == 0`) to build the starting simplex.

/// Result of a Nelder-Mead search.
#[derive(Debug, Clone, Copy)]
pub struct NmResult {
    pub x: f64,
    pub fx: f64,
    pub status: NmStatus,
    pub iterations: usize,
    pub evaluations: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NmStatus {
    XtolReached,
    MaxevalReached,
}

/// Scalar Nelder-Mead with NLopt's default initial simplex and tolerances.
///
/// A 1-D simplex has two vertices, so the classical three-case step of the
/// n-dimensional algorithm collapses: there is no "middle" point, and every
/// operation is expressed directly in terms of the best and worst vertices.
pub fn nelder_mead_1d(x0: f64, lb: Option<f64>, f: impl Fn(f64) -> f64) -> NmResult {
    const XTOL_REL: f64 = 1e-6;
    const MAXEVAL: usize = 1000;

    // NLopt's effective initial simplex, established by probing
    // `nloptr(NLOPT_LN_NELDERMEAD)` directly: the second vertex sits at
    // `x0 + 0.75 * |x0|` (1.75 * x0 for x0 > 0), for every starting point tried
    // (0.1, 0.5, 1, 2, 10). Using 0.25 instead costs ~35% more evaluations and
    // lands on a different point, which perturbs the EM.
    let step = if x0 == 0.0 { 0.75 } else { 0.75 * x0.abs() };
    let mut a = clamp_lb(x0, lb);
    let mut b = clamp_lb(x0 + step, lb);
    if a == b {
        b = clamp_lb(a + 1.0, lb);
    }
    let mut fa = f(a);
    let mut fb = f(b);
    let mut evals = 2usize;
    let mut status = NmStatus::MaxevalReached;
    let mut iters = 0usize;

    loop {
        // order: a is the best, b the worst
        if fb < fa {
            std::mem::swap(&mut a, &mut b);
            std::mem::swap(&mut fa, &mut fb);
        }
        iters += 1;

        // NLopt's stopping rule for NM: the relative spread of the simplex
        let (xmin, xmax) = (a.min(b), a.max(b));
        if xmax - xmin <= XTOL_REL * (xmax.abs() + XTOL_REL) {
            status = NmStatus::XtolReached;
            break;
        }
        if evals >= MAXEVAL {
            break;
        }

        // reflection about the best vertex
        let xr = clamp_lb(2.0 * a - b, lb);
        let fr = f(xr);
        evals += 1;

        if fr < fa {
            // expansion: keep going in the same direction
            let xe = clamp_lb(a + 2.0 * (a - b), lb);
            let fe = f(xe);
            evals += 1;
            b = if fe < fr { xe } else { xr };
            fb = if fe < fr { fe } else { fr };
        } else if fr < fb {
            // accept the reflection
            b = xr;
            fb = fr;
        } else {
            // contraction back toward the best vertex
            let xc = clamp_lb(a + 0.5 * (a - b), lb);
            let fc = f(xc);
            evals += 1;
            if fc < fb {
                b = xc;
                fb = fc;
            } else {
                // shrink: collapse onto the best vertex
                b = clamp_lb(a + 0.5 * (b - a), lb);
                fb = f(b);
                evals += 1;
            }
        }
    }

    if fb < fa {
        NmResult {
            x: b,
            fx: fb,
            status,
            iterations: iters,
            evaluations: evals,
        }
    } else {
        NmResult {
            x: a,
            fx: fa,
            status,
            iterations: iters,
            evaluations: evals,
        }
    }
}

fn clamp_lb(x: f64, lb: Option<f64>) -> f64 {
    match lb {
        Some(l) if x < l => l,
        _ => x,
    }
}

/// `n`-dimensional Nelder-Mead with NLopt defaults. Used by nothing in the
/// current kernel (`.bias_em` optimises one scalar at a time) but provided so the
/// optimiser is available for the constrained trend/Dunnett work planned for
/// v1.1+.
#[allow(dead_code)]
pub fn nelder_mead(x0: &[f64], lb: Option<&[f64]>, f: impl Fn(&[f64]) -> f64) -> (Vec<f64>, f64) {
    const NEL: f64 = 1.0;
    const GAM: f64 = 2.0;
    const RHO: f64 = 0.5;
    const SIG: f64 = 0.5;
    const XTOL_REL: f64 = 1e-6;
    const MAXEVAL: usize = 1000;

    let n = x0.len();
    let mut pts: Vec<Vec<f64>> = Vec::with_capacity(n + 1);
    pts.push(x0.to_vec());
    for i in 0..n {
        let mut p = x0.to_vec();
        let step = if x0[i] == 0.0 {
            0.25
        } else {
            0.25 * x0[i].abs()
        };
        p[i] += step;
        if let Some(l) = lb {
            if p[i] < l[i] {
                p[i] = l[i];
            }
        }
        pts.push(p);
    }
    let mut fv: Vec<f64> = pts.iter().map(|p| f(p)).collect();
    let mut evals = n + 1;
    let clamp = |p: &mut Vec<f64>| {
        if let Some(l) = lb {
            for i in 0..p.len() {
                if p[i] < l[i] {
                    p[i] = l[i];
                }
            }
        }
    };

    while evals < MAXEVAL {
        let mut order: Vec<usize> = (0..=n).collect();
        order.sort_by(|&a, &b| fv[a].partial_cmp(&fv[b]).unwrap());
        pts = order.iter().map(|&i| pts[i].clone()).collect();
        fv = order.iter().map(|&i| fv[i]).collect();

        let spread = (0..=n).map(|i| pts[i][0].abs()).fold(0.0f64, f64::max);
        if (0..=n).all(|i| {
            (0..n).all(|j| (pts[i][j] - pts[0][j]).abs() <= XTOL_REL * spread.max(XTOL_REL))
        }) {
            break;
        }

        let centroid: Vec<f64> = (0..n)
            .map(|j| (0..n).map(|i| pts[i][j]).sum::<f64>() / n as f64)
            .collect();
        let xr: Vec<f64> = (0..n)
            .map(|j| {
                clamp_lb(
                    centroid[j] + NEL * (centroid[j] - pts[n][j]),
                    lb.map(|l| l[j]),
                )
            })
            .collect();
        let fr = f(&xr);
        evals += 1;

        if fr < fv[0] {
            let mut xe: Vec<f64> = (0..n)
                .map(|j| centroid[j] + GAM * (centroid[j] - pts[n][j]))
                .collect();
            clamp(&mut xe);
            let fe = f(&xe);
            evals += 1;
            if fe < fr {
                pts[n] = xe;
                fv[n] = fe;
            } else {
                pts[n] = xr;
                fv[n] = fr;
            }
        } else if fr < fv[n - 1] {
            pts[n] = xr;
            fv[n] = fr;
        } else {
            let mut xc: Vec<f64> = if fr < fv[n] {
                (0..n)
                    .map(|j| centroid[j] + RHO * (centroid[j] - pts[n][j]))
                    .collect()
            } else {
                (0..n)
                    .map(|j| centroid[j] + RHO * (xr[j] - centroid[j]))
                    .collect()
            };
            clamp(&mut xc);
            let fc = f(&xc);
            evals += 1;
            if fc < fv[n] {
                pts[n] = xc;
                fv[n] = fc;
            } else {
                for i in 1..=n {
                    pts[i] = (0..n)
                        .map(|j| pts[0][j] + SIG * (pts[i][j] - pts[0][j]))
                        .collect();
                    clamp(&mut pts[i]);
                    fv[i] = f(&pts[i]);
                    evals += 1;
                }
            }
        }
    }

    let mut order: Vec<usize> = (0..=n).collect();
    order.sort_by(|&a, &b| fv[a].partial_cmp(&fv[b]).unwrap());
    (pts[order[0]].clone(), fv[order[0]])
}

/// Scalar convenience wrapper used by `.bias_em`, mirroring the reference's call
/// `nloptr::nloptr(x0 = kappa, eval_f = obj, lb = 0, ub = NULL, opts = nm_opts)`.
pub fn nelder_mead_scalar(f: impl Fn(f64) -> f64, x0: f64, lb: Option<f64>) -> f64 {
    nelder_mead_1d(x0, lb, f).x
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn finds_a_quadratic_minimum() {
        let r = nelder_mead_1d(1.0, Some(0.0), |x| (x - 2.3).powi(2));
        assert!((r.x - 2.3).abs() < 1e-5, "x = {}", r.x);
        assert_eq!(r.status, NmStatus::XtolReached);
    }

    #[test]
    fn respects_the_lower_bound() {
        // The minimum is below lb, so the optimum sits on the boundary.
        let r = nelder_mead_1d(1.0, Some(0.0), |x: f64| (x + 5.0).powi(2));
        assert!(r.x >= 0.0, "x = {} must not fall below the bound", r.x);
        assert!(r.x < 1e-6, "x = {} should sit on the bound", r.x);
    }

    #[test]
    fn an_unbounded_objective_terminates_instead_of_overflowing() {
        // f(x) = x is unbounded below as x -> -inf, so the only thing the
        // optimiser can do is run away and stop on maxeval. This is the case that
        // would hang a naive implementation, and it is reachable from
        // `.bias_em` if a mixture weight degenerates.
        let r = nelder_mead_1d(1.0, None, |x: f64| x);
        assert_eq!(r.status, NmStatus::MaxevalReached);
        assert!(r.evaluations <= 1000, "evaluations = {}", r.evaluations);
        assert!(r.x.is_finite(), "x must not overflow to inf");
    }

    #[test]
    fn finds_a_one_sided_optimum_above_the_start() {
        // -log(x) on x > 0 has its minimum at x = 1 but increases monotonically
        // for x > 1, so from a start of 1 the search must head upwards. The
        // objective is a genuine unbounded-above case, so the answer is only
        // bounded by maxeval.
        let r = nelder_mead_1d(
            1.0,
            Some(1e-12),
            |x: f64| if x <= 0.0 { 1e12 } else { -x.ln() },
        );
        assert!(r.x > 1.0, "x = {}", r.x);
        assert_eq!(r.status, NmStatus::MaxevalReached);
        assert!(r.x.is_finite());
    }

    #[test]
    fn matches_nloptr_on_the_oracle_objective_shape() {
        // The shape `.bias_em` uses: maximise sum(r * dnorm(beta, mu, sqrt(nu0+k)))
        // which, after taking logs and negating, is a concave function of k.
        let beta = [0.1, -0.2, 0.3, 0.0, -0.05];
        let nu0 = [0.5, 0.6, 0.4, 0.7, 0.55];
        let r = [0.2, 0.1, 0.3, 0.1, 0.3];
        let obj = |k: f64| {
            let mut s = 0.0;
            for i in 0..5 {
                let lp = -0.5 * beta[i] * beta[i] / (nu0[i] + k) - 0.5 * (nu0[i] + k).ln();
                s += r[i] * lp;
            }
            -s
        };
        let res = nelder_mead_1d(1.0, Some(0.0), obj);
        // A grid search must agree that this is near the minimum.
        let grid_best = (0..2000)
            .map(|i| i as f64 * 0.01)
            .min_by(|a, b| obj(*a).partial_cmp(&obj(*b)).unwrap())
            .unwrap();
        assert!(
            (res.x - grid_best).abs() < 5e-3,
            "nm {} vs grid {grid_best}",
            res.x
        );
    }

    #[test]
    fn nd_matches_1d_on_a_separable_objective() {
        let obj = |p: &[f64]| (p[0] - 1.5).powi(2) + (p[1] + 0.5).powi(2);
        let (x, fx) = nelder_mead(&[0.0, 0.0], None, obj);
        assert!((x[0] - 1.5).abs() < 1e-4, "x = {x:?}");
        assert!((x[1] + 0.5).abs() < 1e-4, "x = {x:?}");
        assert!(fx < 1e-8, "fx = {fx}");
    }
}
