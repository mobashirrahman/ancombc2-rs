//! Shared builders for the property suite: deterministic synthetic count tables,
//! designs built the way `model.matrix` would, and thin wrappers over the
//! pipeline.
//!
//! Everything here is a pure function of `(n, m, seed)`, so a failure is
//! reproducible by construction.

// This module is compiled once per test binary, and each binary uses a
// different part of it, so a helper that one of them never calls is dead code
// *in that binary* and nothing else. The alternative -- a shared crate -- would
// mean publishing these builders as an API, which is a much larger commitment
// than the test suite warrants.
#![allow(dead_code)]

use ancombc2_core::config::AncombcConfig;
use ancombc2_core::matrix::Matrix;
use ancombc2_core::preprocess::CountMatrix;
use ancombc2_core::{ancombc2_run_named, AncombcResult};

/// splitmix64, so a "seed" gives the same stream in any language and any run.
fn splitmix64(state: &mut u64) -> u64 {
    *state = state.wrapping_add(0x9E3779B97F4A7C15);
    let mut z = *state;
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
    z ^ (z >> 31)
}

struct Rng(u64);

impl Rng {
    fn new(seed: u64) -> Self {
        Rng(seed.wrapping_mul(0x2545F4914F6CDD1D).wrapping_add(1))
    }

    /// Uniform on [0, 1).
    fn f64(&mut self) -> f64 {
        (splitmix64(&mut self.0) >> 11) as f64 / (1u64 << 53) as f64
    }

    /// Standard normal, via Box-Muller.
    fn normal(&mut self) -> f64 {
        let u1 = self.f64().max(1e-300);
        let u2 = self.f64();
        (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos()
    }

    fn below(&mut self, p: f64) -> bool {
        self.f64() < p
    }
}

/// A `taxa x samples` count table with a known differential-abundance structure.
///
/// A deterministic fraction of taxa carry a real group effect of known
/// magnitude, the rest are null, and every sample carries a lognormal depth
/// offset so the sampling-fraction correction has something to remove.
pub fn counts(n_tax: usize, n_samp: usize, seed: u64) -> CountMatrix {
    let mut rng = Rng::new(seed);
    let mut data = vec![0.0; n_tax * n_samp];
    for i in 0..n_tax {
        let base = 20.0 + 80.0 * rng.f64();
        let is_signal = i % 5 == 0;
        let effect = if is_signal {
            1.0 + 0.5 * rng.f64()
        } else {
            0.0
        };
        let sign = if rng.below(0.5) { -1.0 } else { 1.0 };
        let spread = 0.3 + 0.4 * rng.f64();
        for j in 0..n_samp {
            let grp = (j * 2) as f64 / n_samp as f64;
            let mu = (base + sign * effect * grp + spread * rng.normal()).exp();
            data[i * n_samp + j] = (rng.f64() * mu).floor();
        }
    }
    CountMatrix::new(n_tax, n_samp, data).expect("counts")
}

/// A small table with no differential abundance, for error paths.
pub fn tiny_counts(n_tax: usize, n_samp: usize, seed: u64) -> CountMatrix {
    let mut rng = Rng::new(seed);
    let mut data = vec![0.0; n_tax * n_samp];
    for v in data.iter_mut() {
        *v = (rng.f64() * 50.0).floor() + 1.0;
    }
    CountMatrix::new(n_tax, n_samp, data).expect("counts")
}

/// `k` group labels, balanced as far as possible, in a shuffled order.
pub fn labels(n_samp: usize, k: usize) -> Vec<usize> {
    let mut v: Vec<usize> = (0..n_samp).map(|i| i % k).collect();
    // deterministic Fisher-Yates
    let mut rng = Rng::new(n_samp as u64 * 7919 + k as u64);
    for i in (1..v.len()).rev() {
        let j = (rng.f64() * (i + 1) as f64) as usize;
        v.swap(i, j.min(i));
    }
    v
}

/// A two-group problem: counts, design, labels.
pub fn two_group(n_tax: usize, n_samp: usize, seed: u64) -> (CountMatrix, Matrix, Vec<usize>) {
    let c = counts(n_tax, n_samp, seed);
    let g = labels(n_samp, 2);
    let (x, _) = design_k_groups(&g, 2, 0.0);
    (c, x, g)
}

/// A `k`-group problem.
pub fn k_group(
    n_tax: usize,
    n_samp: usize,
    k: usize,
    seed: u64,
) -> (CountMatrix, Matrix, Vec<usize>) {
    let c = counts(n_tax, n_samp, seed);
    let g = labels(n_samp, k);
    let (x, _) = design_k_groups(&g, k, 0.0);
    (c, x, g)
}

/// `(Intercept) + grp2 + ... + grpK`, with the group levels in ascending order.
///
/// Names match R's `model.matrix` output for a factor under the default
/// `contr.treatment`, because the multi-group tests locate the group columns by
/// substring and the result tables carry these names.
pub fn design_k_groups(g: &[usize], k: usize, extra: f64) -> (Matrix, Vec<String>) {
    let n = g.len();
    let mut cols: Vec<Vec<f64>> = vec![vec![1.0; n]];
    let mut names: Vec<String> = vec!["(Intercept)".into()];
    for lev in 2..=k {
        cols.push(
            (0..n)
                .map(|j| if g[j] + 1 == lev { 1.0 } else { 0.0 })
                .collect(),
        );
        names.push(format!("grp{lev}"));
    }
    if extra > 0.0 {
        // a constant column makes the design rank deficient, which several
        // properties need
        cols.push(vec![extra; n]);
        names.push("extra".into());
    }
    let m = Matrix::from_cols(&cols);
    (m.with_colnames(names.clone()), names)
}

/// `(Intercept) + grp2`, the common two-group design.
pub fn design_binary(g: &[usize], extra: f64) -> Matrix {
    design_k_groups(g, 2, extra).0
}

/// Repeat every sample twice, keeping the design and labels aligned.
pub fn duplicate_samples(
    c: &CountMatrix,
    x: &Matrix,
    g: &[usize],
) -> (CountMatrix, Matrix, Vec<usize>) {
    let m = c.n_samp;
    let mut samples: Vec<usize> = Vec::with_capacity(2 * m);
    for j in 0..m {
        samples.push(j);
        samples.push(j);
    }
    let c2 = c.select(&(0..c.n_taxa).collect::<Vec<_>>(), &samples);
    let x2 = x.select_rows(&samples);
    let g2: Vec<usize> = samples.iter().map(|&j| g[j]).collect();
    (c2, x2, g2)
}

/// The default configuration: no filtering that would drop rows, no tests that
/// need a group, and the reference's tolerances.
pub fn default_config(p: usize, group: &str) -> AncombcConfig {
    AncombcConfig {
        fix_eff: (0..p)
            .map(|k| {
                if k == 0 {
                    "(Intercept)".to_string()
                } else {
                    format!("{group}{k}")
                }
            })
            .collect(),
        prv_cut: 0.0,
        lib_cut: 0.0,
        pseudo_sens: false,
        ..Default::default()
    }
}

pub fn run(c: &CountMatrix, x: &Matrix, g: &[usize], cfg: &AncombcConfig) -> AncombcResult {
    run_with(c, x, g, cfg, cfg.pseudo)
}

pub fn run_with(
    c: &CountMatrix,
    x: &Matrix,
    g: &[usize],
    cfg: &AncombcConfig,
    pseudo: f64,
) -> AncombcResult {
    let mut cfg = cfg.clone();
    cfg.pseudo = pseudo;
    ancombc2_run_named(c, x, Some(g), &cfg, &[], &[], &ancombc2_core::F64Reductions)
        .expect("the run must succeed")
}
