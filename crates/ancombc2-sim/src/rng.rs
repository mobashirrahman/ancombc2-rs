//! A self-contained, reproducible generator for the simulation harness.
//!
//! # Why not reuse the oracle's RNG
//!
//! The simulation grid does not need R's `Mersenne-Twister` and the bit stream
//! to agree: the Rust and R arms of a cell are handed *the same* generated
//! count table, so they see identical data by construction. What the grid does
//! need is that a given `(grid, cell, rep)` regenerates bit-for-bit on any
//! machine, months later, so that a reported FDR can be re-derived. That is what
//! [`Rng`] provides: a fixed algorithm, seeded from the run seed and the cell
//! and rep indices, with no dependence on thread count, platform, or the
//! standard library's internals.
//!
//! Every distribution is sampled by a named algorithm so the number is
//! reproducible from this file alone.

/// SplitMix64. Small, fast, and with a well-defined bit stream.
#[derive(Clone)]
pub struct Rng {
    state: u64,
}

impl Rng {
    /// Seed from the run seed, the cell index, and the rep index.
    ///
    /// Mixing the three means two cells cannot produce correlated streams even
    /// with adjacent indices, and that a rep is reproducible without replaying
    /// the reps before it.
    pub fn for_rep(run_seed: u64, cell: usize, rep: usize) -> Self {
        let mut s = run_seed ^ 0x9e37_79b9_7f4a_7c15;
        s = mix(s ^ (cell as u64).wrapping_mul(0xff51_afd7_ed55_8ccd));
        s = mix(s ^ (rep as u64).wrapping_mul(0xc4ce_b9fe_1a85_ec53));
        Self { state: mix(s) }
    }

    /// The stream for a *cell*, with no rep component.
    ///
    /// The community -- which taxa exist, which are differentially abundant, in
    /// which direction, and how abundant each is -- is a property of the cell,
    /// not of the replicate. A replicate resamples the *sampling*: library
    /// depths, which cells are structurally zero, and the counts.
    ///
    /// Keeping the two apart is what makes the per-taxon metrics mean anything.
    /// If the community were redrawn per replicate, then `taxon00005` in rep 0
    /// and `taxon00005` in rep 1 would be unrelated taxa, and the SE
    /// calibration ratio -- `mean(se) / sd(beta across reps)` -- would be
    /// comparing a taxon's standard error to the spread of a *different*
    /// taxon's estimates, which is a number with no interpretation. It would
    /// also inflate the between-replicate variance of every rate, so the Monte
    /// Carlo error would be far wider than the design needs.
    pub fn for_cell(run_seed: u64, cell: usize) -> Self {
        let mut s = run_seed ^ 0x9e37_79b9_7f4a_7c15;
        s = mix(s ^ (cell as u64).wrapping_mul(0xff51_afd7_ed55_8ccd));
        // A distinct constant from the rep mixer, so a cell stream cannot
        // coincide with a rep stream.
        s = mix(s ^ 0x5851_f42d_4c95_7f2d);
        Self { state: mix(s) }
    }

    pub fn from_seed(seed: u64) -> Self {
        Self { state: mix(seed) }
    }

    /// Next `u64`.
    pub fn next_u64(&mut self) -> u64 {
        self.state = self.state.wrapping_add(0x9e37_79b9_7f4a_7c15);
        let mut z = self.state;
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
        z ^ (z >> 31)
    }

    /// Uniform on `[0, 1)`, using 53 significant bits.
    pub fn uniform(&mut self) -> f64 {
        ((self.next_u64() >> 11) as f64) * (1.0 / ((1u64 << 53) as f64))
    }

    /// Uniform on `[lo, hi)`.
    pub fn uniform_range(&mut self, lo: f64, hi: f64) -> f64 {
        lo + (hi - lo) * self.uniform()
    }

    /// Standard normal, Box-Muller.
    ///
    /// The polar form is deliberately avoided: it discards a variate per draw,
    /// which makes the stream depend on a rejection loop rather than only on
    /// the number of draws.
    pub fn normal(&mut self) -> f64 {
        // 1 - uniform() keeps the argument of `ln` strictly positive even when
        // `uniform()` returns exactly 0.
        let u1 = 1.0 - self.uniform();
        let u2 = self.uniform();
        (-2.0 * u1.ln()).sqrt() * (std::f64::consts::TAU * u2).cos()
    }

    /// `Normal(mean, sd)`.
    pub fn normal_sd(&mut self, mean: f64, sd: f64) -> f64 {
        mean + sd * self.normal()
    }

    /// `Poisson(lambda)` for `lambda < 30`, by Knuth's product method.
    fn poisson_knuth(&mut self, lambda: f64) -> u64 {
        let l = (-lambda).exp();
        let mut k = 0u64;
        let mut p = 1.0f64;
        loop {
            p *= self.uniform();
            if p <= l {
                return k;
            }
            k += 1;
        }
    }

    /// `Poisson(lambda)` for `lambda >= 30`, by the normal approximation with
    /// a continuity correction and rejection back into the lower tail.
    fn poisson_normal(&mut self, lambda: f64) -> u64 {
        let sd = lambda.sqrt();
        for _ in 0..64 {
            let v = self.normal_sd(lambda, sd);
            if v < 0.0 {
                continue;
            }
            let v = v.floor();
            // Attenuate back towards the mean where the normal approximation
            // is poor, which is why the guard below never rejects for large
            // lambda.
            if v > lambda + 6.0 * sd {
                continue;
            }
            return v as u64;
        }
        lambda.round() as u64
    }

    /// `Poisson(lambda)`.
    pub fn poisson(&mut self, lambda: f64) -> u64 {
        if !lambda.is_finite() || lambda <= 0.0 {
            return 0;
        }
        if lambda < 30.0 {
            self.poisson_knuth(lambda)
        } else {
            self.poisson_normal(lambda)
        }
    }

    /// `Gamma(shape, scale)` for `shape >= 1`, by Marsaglia-Tsang squeeze.
    pub fn gamma(&mut self, shape: f64, scale: f64) -> f64 {
        debug_assert!(shape >= 1.0, "gamma expects shape >= 1; use gamma_lt_one");
        if shape == 1.0 {
            // Exponential.
            return -scale * (1.0 - self.uniform()).ln();
        }
        // Marsaglia-Tsang: with `d = shape - 1/3` the returned variate is
        // `d (1 + c X)^3` for `X` standard normal and `c = 1 / sqrt(9 d)`. The
        // moment check is `E[(1 + cX)^3] = 1 + 3c^2 = 1 + 1/(3d)`, so
        // `E[d(1 + cX)^3] = d + 1/3 = shape`. Note the reciprocal: with
        // `c = sqrt(9 d)` the mean is wrong by a factor of `9d`, which is
        // invisible for `shape = 1` because that case returns early.
        let d = shape - 1.0 / 3.0;
        let c = 1.0 / (9.0 * d).sqrt();
        loop {
            let x = self.normal();
            let v = 1.0 + c * x;
            if v <= 0.0 {
                continue;
            }
            let v = v * v * v;
            let u = self.uniform();
            if u < 1.0 - 0.0331 * x * x * x * x {
                return d * v * scale;
            }
            if u.ln() < 0.5 * x * x + d * (1.0 - v + v.ln()) {
                return d * v * scale;
            }
        }
    }

    /// `Gamma(shape, scale)` for `shape < 1`, by the boost trick: sample
    /// `Gamma(shape + 1)` and rescale, which is the standard construction
    /// because the Marsaglia-Tsang method above needs `shape >= 1`.
    pub fn gamma_lt_one(&mut self, shape: f64, scale: f64) -> f64 {
        let g = self.gamma(shape + 1.0, 1.0);
        let u = self.uniform();
        g * (u.max(f64::MIN_POSITIVE)).powf(1.0 / shape) * scale
    }

    /// `Gamma(shape, scale)` for any positive `shape`.
    pub fn gamma_any(&mut self, shape: f64, scale: f64) -> f64 {
        if shape >= 1.0 {
            self.gamma(shape, scale)
        } else {
            self.gamma_lt_one(shape, scale)
        }
    }

    /// A Bernoulli draw.
    pub fn bernoulli(&mut self, p: f64) -> bool {
        self.uniform() < p
    }

    /// `LogNormal(log(mean), cv)`, parameterised by the *mean* and the
    /// coefficient of variation rather than by a log-scale `sd`, because that
    /// is how the grid specifies library sizes.
    ///
    /// For `L ~ LogNormal(mu, s^2)` the mean is `exp(mu + s^2/2)` and the CV is
    /// `sqrt(exp(s^2) - 1)`, so `s^2 = log(1 + cv^2)` and
    /// `mu = log(mean) - s^2/2`.
    pub fn lognormal_mean_cv(&mut self, mean: f64, cv: f64) -> f64 {
        if mean <= 0.0 {
            return 0.0;
        }
        if cv <= 0.0 {
            return mean;
        }
        let s2 = (1.0 + cv * cv).ln();
        let mu = mean.ln() - 0.5 * s2;
        // Rounding to an integer is the library's *depth*; the draw is a depth
        // count, not a rate. The exponentiation is the whole point of the
        // log-normal: `mu` is a location on the log scale, not the value.
        self.normal_sd(mu, s2.sqrt()).exp().max(1.0).round()
    }
}

fn mix(mut z: u64) -> u64 {
    z = (z ^ (z >> 33)).wrapping_mul(0xff51_afd7_ed55_8ccd);
    z = (z ^ (z >> 33)).wrapping_mul(0xc4ce_b9fe_1a85_ec53);
    z ^ (z >> 33)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_seed_reproduces_the_stream() {
        let mut a = Rng::for_rep(42, 3, 7);
        let mut b = Rng::for_rep(42, 3, 7);
        for _ in 0..100 {
            assert_eq!(a.next_u64(), b.next_u64());
        }
    }

    #[test]
    fn a_cell_stream_is_the_same_for_every_rep_of_that_cell() {
        let mut a = Rng::for_cell(42, 3);
        let mut b = Rng::for_cell(42, 3);
        for _ in 0..100 {
            assert_eq!(a.next_u64(), b.next_u64());
        }
        let d: Vec<u64> = (0..4).map(|_| Rng::for_cell(42, 4).next_u64()).collect();
        let e: Vec<u64> = (0..4).map(|_| Rng::for_cell(42, 3).next_u64()).collect();
        assert_ne!(d, e, "different cells must not share a community");
    }

    #[test]
    fn a_cell_stream_cannot_collide_with_a_rep_stream() {
        let cell: Vec<u64> = (0..8).map(|_| Rng::for_cell(7, 5).next_u64()).collect();
        for rep in 0..4 {
            let r: Vec<u64> = (0..8).map(|_| Rng::for_rep(7, 5, rep).next_u64()).collect();
            assert_ne!(cell, r, "rep {rep} collides with the cell stream");
        }
    }

    #[test]
    fn cells_and_reps_do_not_share_a_stream() {
        let a: Vec<u64> = (0..8).map(|_| Rng::for_rep(42, 1, 1).next_u64()).collect();
        let b: Vec<u64> = (0..8).map(|_| Rng::for_rep(42, 2, 1).next_u64()).collect();
        assert_ne!(a, b, "adjacent cells must not share a stream");
    }

    #[test]
    fn uniform_stays_in_range_and_centres() {
        let mut r = Rng::from_seed(1);
        let mut sum = 0.0;
        let n = 200_000;
        for _ in 0..n {
            let u = r.uniform();
            assert!((0.0..1.0).contains(&u));
            sum += u;
        }
        assert!(
            (sum / n as f64 - 0.5).abs() < 0.01,
            "mean {:.4}",
            sum / n as f64
        );
    }

    #[test]
    fn normal_has_the_requested_moments() {
        let mut r = Rng::from_seed(2);
        let n = 400_000;
        let (mut s, mut s2) = (0.0, 0.0);
        for _ in 0..n {
            let v = r.normal();
            s += v;
            s2 += v * v;
        }
        let m = s / n as f64;
        let v = s2 / n as f64 - m * m;
        assert!(m.abs() < 0.02, "mean {m}");
        assert!((v.sqrt() - 1.0).abs() < 0.02, "sd {}", v.sqrt());
    }

    #[test]
    fn poisson_matches_its_mean_and_variance() {
        let mut r = Rng::from_seed(3);
        // For a Poisson, mean == variance == lambda. Both are estimated from
        // the same draws, which makes this a joint check on both.
        for lambda in [0.5f64, 5.0, 50.0, 5000.0] {
            let n = 200_000;
            let (mut s, mut s2) = (0.0, 0.0);
            for _ in 0..n {
                let v = r.poisson(lambda) as f64;
                s += v;
                s2 += v * v;
            }
            let m = s / n as f64;
            let v = s2 / n as f64 - m * m;
            assert!(
                (m - lambda).abs() < 0.03 * lambda.max(1.0) + 0.05,
                "lambda {lambda}: mean {m}"
            );
            assert!(
                (v - lambda).abs() < 0.05 * lambda.max(1.0) + 0.1,
                "lambda {lambda}: var {v}"
            );
        }
    }

    #[test]
    fn gamma_matches_its_shape_and_scale() {
        let mut r = Rng::from_seed(4);
        for (shape, scale) in [(0.5f64, 1.0), (1.0, 2.0), (2.5, 3.0), (10.0, 0.5)] {
            let n = 200_000;
            let (mut s, mut s2) = (0.0, 0.0);
            for _ in 0..n {
                let v = r.gamma_any(shape, scale);
                s += v;
                s2 += v * v;
            }
            let m = s / n as f64;
            let v = s2 / n as f64 - m * m;
            assert!(
                (m - shape * scale).abs() < 0.02 * shape * scale + 0.01,
                "shape {shape} scale {scale}: mean {m}"
            );
            assert!(
                (v.sqrt() - (shape * scale * scale).sqrt()).abs()
                    < 0.03 * (shape * scale * scale).sqrt() + 0.01,
                "shape {shape} scale {scale}: sd {}",
                v.sqrt()
            );
        }
    }

    #[test]
    fn lognormal_mean_cv_hits_its_specification() {
        let mut r = Rng::from_seed(5);
        for (mean, cv) in [(1e3f64, 0.0f64), (1e4, 0.3), (1e5, 0.6)] {
            let n = 100_000;
            let xs: Vec<f64> = (0..n).map(|_| r.lognormal_mean_cv(mean, cv)).collect();
            let m = xs.iter().sum::<f64>() / n as f64;
            let v = xs.iter().map(|x| (x - m) * (x - m)).sum::<f64>() / n as f64;
            let want_cv = if cv == 0.0 { 0.0 } else { v.sqrt() / m };
            assert!(
                (m - mean).abs() < 0.02 * mean,
                "mean {mean} cv {cv}: got {m}"
            );
            assert!(
                (want_cv - cv).abs() < 0.03 * cv + 0.005,
                "cv {cv}: got {want_cv}"
            );
        }
    }
}
