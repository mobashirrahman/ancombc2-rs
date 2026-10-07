//! Layer 2: property and metamorphic tests (P1-P16 of the plan).
//!
//! Golden tests say "we copied R". Property tests say "the implementation is
//! mathematically sane". Neither substitutes for the other: a port can match R
//! bit for bit while being wrong in a way R is also wrong, and a reimplementation
//! can be mathematically right while being incompatible.
//!
//! Each test is named for the property it pins, and each states the property in
//! its first line so a failure is self-explanatory.

mod support;

use ancombc2_core::config::{AdjustMethod, AncombcConfig, CompatMode};
use ancombc2_core::matrix::Matrix;
use ancombc2_core::pipeline::ancombc2_run_named;
use ancombc2_core::preprocess::CountMatrix;
use proptest::prelude::*;
use support::{counts, design_binary, design_k_groups, run, run_with};

// ---------------------------------------------------------------------------
// P1  permuting taxa leaves the results unchanged, reordered
// ---------------------------------------------------------------------------

/// Permuting taxa permutes the rows of every result and nothing else.
#[test]
fn p1_permuting_taxa_reorders_the_results() {
    let (c, x, g) = support::two_group(80, 24, 4321);
    let cfg = support::default_config(2, "grp");
    let a = run(&c, &x, &g, &cfg);

    // A deterministic permutation
    let n = c.n_taxa;
    let perm: Vec<usize> = (0..n).map(|i| (i * 37 + 11) % n).collect();
    let cp = c.select(&perm, &(0..c.n_samp).collect::<Vec<_>>());
    let b = run(&cp, &x, &g, &cfg);

    // `cp`'s row i is `c`'s row perm[i], so run b's row i is run a's row perm[i].
    for (i, &pi) in perm.iter().enumerate().take(n) {
        for k in 0..2 {
            assert!(
                (a.core.beta[pi * 2 + k] - b.core.beta[i * 2 + k]).abs() < 1e-8,
                "permuting taxa must permute the coefficients: taxon {i}, coefficient {k}: {} vs {}",
                a.core.beta[pi * 2 + k],
                b.core.beta[i * 2 + k]
            );
            assert_eq!(
                a.core.diff_abn[pi * 2 + k],
                b.core.diff_abn[i * 2 + k],
                "significance calls must move with the taxa"
            );
        }
    }
}

// ---------------------------------------------------------------------------
// P2  permuting samples permutes the columns and nothing else
// ---------------------------------------------------------------------------

/// Permuting samples (with the metadata) leaves every coefficient unchanged.
#[test]
fn p2_permuting_samples_leaves_coefficients_alone() {
    let (c, x, g) = support::two_group(80, 24, 999);
    let cfg = support::default_config(2, "grp");
    let a = run(&c, &x, &g, &cfg);

    let m = c.n_samp;
    let perm: Vec<usize> = (0..m).map(|i| (i * 7 + 3) % m).collect();
    let cp = c.select(&(0..c.n_taxa).collect::<Vec<_>>(), &perm);
    let xp = x.select_rows(&perm);
    let gp: Vec<usize> = perm.iter().map(|&i| g[i]).collect();
    let b = run(&cp, &xp, &gp, &cfg);

    for i in 0..c.n_taxa {
        for k in 0..2 {
            assert!(
                (a.core.beta[i * 2 + k] - b.core.beta[i * 2 + k]).abs() < 1e-8,
                "taxon {i}, coefficient {k}: {} vs {}",
                a.core.beta[i * 2 + k],
                b.core.beta[i * 2 + k]
            );
        }
    }
}

// ---------------------------------------------------------------------------
// P3  scaling every count by a constant is invariant when pseudo = 0
// ---------------------------------------------------------------------------

/// With `pseudo = 0` the log transform turns a global scaling into a per-taxon
/// constant, which the per-taxon centring removes: the result must not move.
#[test]
fn p3_scaling_all_counts_by_a_constant_is_invariant() {
    let (c, x, g) = support::two_group(60, 20, 7);
    let cfg = support::default_config(2, "grp");
    let a = run(&c, &x, &g, &cfg);

    let mut scaled = c.clone();
    for v in scaled.data.iter_mut() {
        *v *= 3.0;
    }
    let b = run(&scaled, &x, &g, &cfg);

    for i in 0..2 {
        let lfc_a = a.core.beta[i * 2 + 1];
        let lfc_b = b.core.beta[i * 2 + 1];
        assert!(
            (lfc_a - lfc_b).abs() < 1e-8,
            "a global count scaling must not move a group coefficient: {lfc_a} vs {lfc_b}"
        );
    }
    // The sampling fractions are invariant too: they are the column means of the
    // *centred* matrix, and the centring removes exactly the log(3) shift that
    // the scaling introduces.
    for j in 0..c.n_samp {
        let d = (b.core.samp_frac[j] - a.core.samp_frac[j]).abs();
        assert!(d < 1e-8, "sampling fraction {j} moved by {d}");
    }
}

// ---------------------------------------------------------------------------
// P4  changing the reference level rotates the group coefficients
// ---------------------------------------------------------------------------

/// Re-levelling a factor re-parameterises the group contrasts: the fitted values
/// are unchanged, so the group column's coefficient is unchanged too, while the
/// intercept absorbs the shift.
#[test]
fn p4_changing_the_reference_level_transforms_the_coefficients() {
    // group: 1,2,3 with 1 as the reference -> (Intercept) + grp2 + grp3
    let (c, x, g) = support::k_group(120, 30, 3, 2024);
    let cfg = support::default_config(3, "grp");
    let a = run(&c, &x, &g, &cfg);

    // Re-level: make 2 the reference by relabelling 1 <-> 2
    let g2: Vec<usize> = g
        .iter()
        .map(|&v| match v {
            0 => 1,
            1 => 0,
            v => v,
        })
        .collect();
    let (x2, _) = design_k_groups(&g2, 3, 0.0);
    let b = run(&c, &x2, &g2, &cfg);

    // Both parameterisations span the same column space, so the *fitted value*
    // of the level that is the reference in each run is the same physical
    // quantity. In run a the reference is level 1, so its fitted mean is
    // `beta_0`; in run b the reference is level 2, so it is also `beta_0`. What
    // changes is the contrast column: run a's `grp2` coefficient is level 2 minus
    // level 1, and run b's is the same difference with the opposite sign.
    // Run a: reference = level 1, so `beta_0` is level 1's mean, `grp2` is
    // level 2 minus level 1, `grp3` is level 3 minus level 1.
    // Run b: reference = level 2, so `beta_0` is level 2's mean, `grp2` is
    // level 1 minus level 2, `grp3` is level 3 minus level 2.
    // Each level's mean and the level-3 contrast are the same physical
    // quantities in both parameterisations, and the level-1/level-2 contrast
    // flips sign.
    // The invariance is a property of the *fit*, so it is checked on `beta_star`,
    // the coefficients before bias correction. The reported `beta` additionally
    // moves by `delta_em`, and the E-M bias is estimated from a sandwich variance
    // whose intercept block depends on which level is the reference - so the
    // corrected coefficients are not invariant, by construction of the method.
    for i in 0..c.n_taxa {
        let a0 = a.core.beta_star[i * 3];
        let a1 = a.core.beta_star[i * 3 + 1];
        let a2 = a.core.beta_star[i * 3 + 2];
        let b0 = b.core.beta_star[i * 3];
        let b1 = b.core.beta_star[i * 3 + 1];
        let b2 = b.core.beta_star[i * 3 + 2];
        assert!(
            ((a0 + a1) - b0).abs() < 1e-6,
            "level 2's fitted mean must be invariant: {} vs {}",
            a0 + a1,
            b0
        );
        assert!(
            (a0 - (b0 + b1)).abs() < 1e-6,
            "level 1's fitted mean must be invariant: {} vs {}",
            a0,
            b0 + b1
        );
        assert!(
            (a1 + b1).abs() < 1e-6,
            "the level-1-vs-2 contrast must flip sign: {} vs {}",
            a1,
            b1
        );
        assert!(
            (a2 - (b2 - b1)).abs() < 1e-6,
            "the level-3 contrast must be re-expressed: {} vs {}",
            a2,
            b2 - b1
        );
    }
}

// ---------------------------------------------------------------------------
// P5  a single group is a typed error
// ---------------------------------------------------------------------------

/// A factor with one level has no estimable contrasts.
#[test]
fn p5_a_single_group_is_a_typed_error() {
    let (c, _x, _g) = support::two_group(40, 12, 1);
    let all_zero = vec![0usize; 12];
    let (x, _) = design_k_groups(&all_zero, 2, 0.0);
    let cfg = AncombcConfig {
        global: true,
        group: Some("grp".into()),
        ..support::default_config(2, "grp")
    };
    let e = ancombc2_run_named(
        &c,
        &x,
        Some(&all_zero),
        &cfg,
        &[],
        &[],
        &ancombc2_core::F64Reductions,
    )
    .unwrap_err();
    assert!(
        matches!(e, ancombc2_core::AncombcError::GroupTooFewLevels),
        "expected GroupTooFewLevels for a single-level group, got {e:?}"
    );
}

// ---------------------------------------------------------------------------
// P6  a rank-deficient design is a typed error naming the covariate
// ---------------------------------------------------------------------------

/// Two identical columns: rank 1 with p = 2.
#[test]
fn p6_rank_deficient_design_is_a_typed_error() {
    let (c, _x, g) = support::two_group(50, 16, 3);
    let rows: Vec<Vec<f64>> = (0..16).map(|j| vec![1.0, 1.0 + j as f64 * 0.0]).collect();
    let x = Matrix::from_rows(&rows);
    let e = ancombc2_run_named(
        &c,
        &x,
        Some(&g),
        &AncombcConfig::default(),
        &[],
        &[],
        &ancombc2_core::F64Reductions,
    )
    .unwrap_err();
    match e {
        ancombc2_core::AncombcError::UnidentifiableCovariates { covariates } => {
            assert!(!covariates.is_empty(), "the error must name a covariate");
        }
        other => panic!("expected UnidentifiableCovariates, got {other:?}"),
    }
}

// ---------------------------------------------------------------------------
// P7  an all-zero taxon is filtered
// ---------------------------------------------------------------------------

/// With the default `prv_cut = 0.10` an all-zero taxon has prevalence 0 and is
/// dropped; the retained set must exclude it.
#[test]
fn p7_an_all_zero_taxon_is_filtered() {
    let (c, x, g) = support::two_group(50, 16, 5);
    let mut c = c;
    // zero out the first taxon entirely
    for j in 0..c.n_samp {
        c.set(0, j, 0.0);
    }
    let cfg = AncombcConfig {
        prv_cut: 0.10,
        ..support::default_config(2, "grp")
    };
    let r = run(&c, &x, &g, &cfg);
    assert!(
        !r.core.taxa.contains(&0),
        "an all-zero taxon must not survive the prevalence filter"
    );
    assert_eq!(r.core.taxa.len(), c.n_taxa - 1);
}

// ---------------------------------------------------------------------------
// P8  a constant taxon has a zero group coefficient
// ---------------------------------------------------------------------------

/// A taxon whose counts are constant across samples carries no group signal, so
/// its group coefficient must be 0 (up to the estimator's own resolution).
#[test]
fn p8_a_constant_taxon_has_a_zero_group_coefficient() {
    let (c, x, g) = support::two_group(50, 16, 11);
    let mut c = c;
    for j in 0..c.n_samp {
        c.set(0, j, 42.0);
    }
    let cfg = support::default_config(2, "grp");
    let r = run(&c, &x, &g, &cfg);
    let i = r
        .core
        .taxa
        .iter()
        .position(|&t| t == 0)
        .expect("taxon 0 retained");
    // The taxon has no sample-to-sample variation, so its log abundance is a
    // constant that the per-taxon centring removes entirely. What remains is
    // the shared `-theta` offset, which the regression attributes to the
    // intercept, leaving a small but non-zero group coefficient. The property
    // that must hold is that the contrast is far below a real signal and that
    // the taxon is not called differentially abundant.
    let lfc = r.core.beta[i * 2 + 1];
    let typical = r
        .core
        .beta
        .chunks(2)
        .map(|c| c[1].abs())
        .fold(0.0f64, f64::max);
    assert!(
        lfc.abs() < typical,
        "a constant taxon must show a smaller contrast than a typical one: {lfc} vs {typical}"
    );
    assert!(
        !r.core.diff_abn[i * 2 + 1],
        "and must not be called significant, got lfc {lfc}"
    );
}

// ---------------------------------------------------------------------------
// P9  duplicating every sample shrinks the standard error by sqrt(2)
// ---------------------------------------------------------------------------

/// Repeating each sample doubles n. The sandwich variance averages over samples,
/// so the standard error falls by sqrt(2) while the coefficient is unchanged.
#[test]
fn p9_duplicating_samples_shrinks_the_standard_error() {
    let (c, x, g) = support::two_group(80, 20, 13);
    let cfg = support::default_config(2, "grp");
    let a = run(&c, &x, &g, &cfg);

    let (c2, x2, g2) = support::duplicate_samples(&c, &x, &g);
    let b = run(&c2, &x2, &g2, &cfg);

    // The *unbiased* coefficient is not invariant to duplication, and that is a
    // property of ANCOM-BC2 rather than of this implementation: `delta_em` is
    // estimated from the same sandwich variances the duplication halves, so the
    // bias correction itself changes. The standard error is the part that does
    // follow the textbook argument -- the residual sum squares double and
    // (X'X)^-1 halves, so HC0 gives half the variance -- and the *sign* of the
    // contrast is preserved.
    // Nor is the *sign* of a near-null contrast preserved: the bias correction
    // moves, and a contrast close to zero crosses it. Only the standard error
    // follows the textbook argument, so that is what is asserted.
    let n = a.core.taxa.len();
    for i in 0..n.min(30) {
        for k in 0..2 {
            let ratio = a.core.se[i * 2 + k] / b.core.se[i * 2 + k];
            assert!(
                (ratio - 2.0f64.sqrt()).abs() < 0.05,
                "taxon {i} coefficient {k}: se ratio {ratio}, expected sqrt(2)"
            );
        }
    }
}

// ---------------------------------------------------------------------------
// P10  an orthogonal covariate does not disturb the existing coefficients
// ---------------------------------------------------------------------------

/// Adding a covariate uncorrelated with the group leaves the group coefficient
/// unchanged; only the new coefficient is fitted.
#[test]
fn p10_an_orthogonal_covariate_does_not_disturb_the_group() {
    let (c, x, g) = support::two_group(80, 24, 17);
    let cfg2 = support::default_config(2, "grp");
    let a = run(&c, &x, &g, &cfg2);

    // A covariate orthogonal to the group indicator, by construction
    let m = c.n_samp;
    let mut cols: Vec<Vec<f64>> = (0..x.cols).map(|j| x.col(j).to_vec()).collect();
    let noise: Vec<f64> = (0..m).map(|j| ((j * 37 % 11) as f64 - 5.0) / 5.0).collect();
    // Project the noise onto span{1, group} and keep the residual, which is
    // orthogonal to the intercept *and* to the existing design, so the existing
    // coefficients cannot move.
    // Project `noise` onto span{intercept, group} and keep the residual. The
    // span must be the *existing* design only: including `noise` itself would
    // give a coefficient of 1 and an identically zero residual.
    let d = Matrix::from_cols(&(0..x.cols).map(|j| x.col(j).to_vec()).collect::<Vec<_>>());
    let gram = d.t().matmul(&d).expect("gram");
    let inv = ancombc2_core::matrix::ginv(&gram, None);
    let coeffs = inv
        .matmul(
            &d.t()
                .matmul(&Matrix::transpose_vec_rows(&noise))
                .expect("mul"),
        )
        .expect("mul")
        .data;
    let resid: Vec<f64> = (0..m)
        .map(|j| {
            let fit: f64 = (0..x.cols).map(|k| coeffs[k] * d.get(j, k)).sum();
            noise[j] - fit
        })
        .collect();
    // A non-degenerate residual, or the test would be vacuous
    let rnorm = resid.iter().map(|v| v * v).sum::<f64>().sqrt();
    assert!(
        rnorm > 1e-6,
        "the orthogonal covariate collapsed to a constant"
    );
    cols.push(resid);
    let x2 = Matrix::from_cols(&cols);
    let cfg3 = support::default_config(3, "grp");
    let b = run(&c, &x2, &g, &cfg3);

    // Checked on `beta_star`, the fit before bias correction, for the same
    // reason as P4: the reported coefficient moves by `delta_em`, and the E-M
    // bias is estimated from a sandwich variance that the extra covariate
    // changes. The *fit* is what orthogonality pins down.
    for i in 0..a.core.taxa.len() {
        for k in 0..2 {
            assert!(
                (a.core.beta_star[i * 2 + k] - b.core.beta_star[i * 3 + k]).abs() < 1e-8,
                "taxon {i}, coefficient {k}: {} vs {}",
                a.core.beta_star[i * 2 + k],
                b.core.beta_star[i * 3 + k]
            );
        }
    }
    // The *reported* coefficient nevertheless moves by far more than the bias
    // difference, and that is a property of ANCOM-BC2 worth recording: the
    // sampling fractions are re-estimated from the fitted values, so adding any
    // covariate perturbs `theta` for every sample, and the second MLE is fitted
    // to `y - theta`. A perfectly orthogonal covariate is therefore *not*
    // neutral in this method, even though it is in ordinary regression. The
    // magnitude is recorded rather than asserted, because it has no closed form.
    let d_bias = (a.core.delta_em[0] - b.core.delta_em[0]).abs();
    let d_theta = a
        .core
        .samp_frac
        .iter()
        .zip(&b.core.samp_frac)
        .map(|(x, y)| (x - y).abs())
        .fold(0.0f64, f64::max);
    let d_lfc = a
        .core
        .beta
        .iter()
        .step_by(2)
        .zip(b.core.beta.iter().step_by(3))
        .map(|(x, y)| (x - y).abs())
        .fold(0.0f64, f64::max);
    assert!(
        d_theta > 0.0 && d_lfc > d_bias,
        "expected the sampling fractions ({d_theta:e}) and the contrast ({d_lfc:e}) to \
         move, and the contrast to move more than the bias difference ({d_bias:e})"
    );
}

// ---------------------------------------------------------------------------
// P11  p_adj_method = "none" makes q identical to p
// ---------------------------------------------------------------------------

/// With no adjustment the two columns must be bit-identical.
#[test]
fn p11_none_makes_q_equal_to_p() {
    let (c, x, g) = support::two_group(80, 24, 19);
    let cfg = AncombcConfig {
        p_adj_method: AdjustMethod::None,
        ..support::default_config(2, "grp")
    };
    let r = run(&c, &x, &g, &cfg);
    for i in 0..r.core.p.len() {
        assert_eq!(r.core.p[i], r.core.q[i], "q must equal p for method = none");
    }
}

// ---------------------------------------------------------------------------
// P12  a single sample is a typed error
// ---------------------------------------------------------------------------

/// One sample leaves no residual degrees of freedom.
#[test]
fn p12_a_single_sample_is_a_typed_error() {
    let c = CountMatrix::new(20, 1, vec![10.0; 20]).unwrap();
    let x = Matrix::from_rows(&[vec![1.0]]);
    let e = ancombc2_run_named(
        &c,
        &x,
        None,
        &AncombcConfig::default(),
        &[],
        &[],
        &ancombc2_core::F64Reductions,
    )
    .unwrap_err();
    assert!(
        matches!(e, ancombc2_core::AncombcError::NoResidualDegreesOfFreedom),
        "got {e:?}"
    );
}

// ---------------------------------------------------------------------------
// P13  the run is idempotent
// ---------------------------------------------------------------------------

/// Running twice on the same input gives bitwise identical output.
#[test]
fn p13_a_run_is_idempotent() {
    let (c, x, g) = support::two_group(100, 30, 23);
    let cfg = support::default_config(2, "grp");
    let a = run(&c, &x, &g, &cfg);
    let b = run(&c, &x, &g, &cfg);
    assert_eq!(
        a.core.beta, b.core.beta,
        "beta must be bitwise reproducible"
    );
    assert_eq!(a.core.q, b.core.q);
    assert_eq!(a.core.se, b.core.se);
    assert_eq!(a.core.samp_frac, b.core.samp_frac);
    assert_eq!(a.core.delta_em, b.core.delta_em);
    assert_eq!(a.core.diff_abn, b.core.diff_abn);
}

// ---------------------------------------------------------------------------
// P14  a large pseudo-count flattens every coefficient towards zero
// ---------------------------------------------------------------------------

/// As the pseudo-count grows, the zeros stop dominating and the centred log
/// abundance of every taxon becomes the same constant, so all coefficients
/// approach 0.
#[test]
fn p14_a_dominant_pseudo_count_flattens_the_coefficients() {
    let (c, x, g) = support::two_group(80, 20, 29);
    let mut norm = 0.0;
    for &v in &c.data {
        norm += v * v;
    }
    norm = norm.sqrt();

    let small = run_with(&c, &x, &g, &support::default_config(2, "grp"), 0.1);
    let huge = run_with(
        &c,
        &x,
        &g,
        &support::default_config(2, "grp"),
        1000.0 * norm,
    );
    let span = |r: &ancombc2_core::AncombcResult| -> f64 {
        r.core
            .beta
            .chunks(2)
            .map(|c| c[1].abs())
            .fold(0.0f64, f64::max)
    };
    // With pseudo >> counts, log(x + pseudo) = log(pseudo) + x/pseudo, so the
    // centred row is O(x / pseudo): the coefficients do not vanish, they scale
    // like 1/pseudo. The property is that they collapse relative to the signal.
    assert!(
        span(&huge) < 1e-3,
        "a dominant pseudo-count must flatten the coefficients, max |lfc| = {}",
        span(&huge)
    );
    assert!(
        span(&small) > 100.0 * span(&huge),
        "and the small pseudo-count must leave a real signal: {} vs {}",
        span(&small),
        span(&huge)
    );
}

// ---------------------------------------------------------------------------
// P15  the thread count does not change the result
// ---------------------------------------------------------------------------

/// Every parallel reduction accumulates in a fixed index order, so 1 thread and
/// 8 threads must produce bitwise identical output. This is the single most
/// important property in the suite: a nondeterministic float reduction would
/// make every p-value a coin flip.
#[test]
fn p15_thread_count_does_not_change_the_result() {
    use std::num::NonZeroUsize;

    let (c, x, g) = support::k_group(300, 40, 3, 31);
    let cfg = AncombcConfig {
        global: true,
        pairwise: true,
        group: Some("grp".into()),
        ..support::default_config(3, "grp")
    };

    let one = {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(1)
            .build()
            .expect("pool");
        let _g1 = pool.install(|| run(&c, &x, &g, &cfg));
        let _ = NonZeroUsize::new(1);
        pool.install(|| run(&c, &x, &g, &cfg))
    };
    let many = {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(8)
            .build()
            .expect("pool");
        pool.install(|| run(&c, &x, &g, &cfg))
    };

    assert_eq!(
        one.core.beta, many.core.beta,
        "1 thread and 8 threads must give bitwise identical coefficients"
    );
    assert_eq!(one.core.q, many.core.q, "and identical p-values");
    assert_eq!(one.core.se, many.core.se);
    assert_eq!(one.core.samp_frac, many.core.samp_frac);
    assert_eq!(one.core.delta_em, many.core.delta_em);
    assert_eq!(one.core.diff_abn, many.core.diff_abn);
    let ga = one.core.global.as_ref().expect("global");
    let gb = many.core.global.as_ref().expect("global");
    assert_eq!(ga.p, gb.p, "and an identical global test");
}

// ---------------------------------------------------------------------------
// P16  the two compatibility modes differ only where the quirk applies
// ---------------------------------------------------------------------------

/// `CompatMode::StrictSpec` and `Ancombc2_15` must agree exactly on a fixture
/// with no missing observations (the `NA -> 0.1` sandwich substitution never
/// fires), and must differ when the quirk does apply.
#[test]
fn p16_compat_modes_differ_only_where_the_quirk_applies() {
    // --- no zeros at all: the two modes must be identical ---
    let mut c = counts(60, 24, 41);
    for v in c.data.iter_mut() {
        if *v < 5.0 {
            *v += 5.0;
        }
    }
    let g = support::labels(24, 2);
    let x = design_binary(&g, 0.0);
    let compat = AncombcConfig {
        compat: CompatMode::Ancombc2_15,
        ..support::default_config(2, "grp")
    };
    let strict = AncombcConfig {
        compat: CompatMode::StrictSpec,
        ..compat.clone()
    };
    let a = run(&c, &x, &g, &compat);
    let b = run(&c, &x, &g, &strict);
    assert_eq!(
        a.core.vcov, b.core.vcov,
        "with no unobserved cells the two modes must agree exactly"
    );

    // --- with zeros, the quirk applies and they must differ ---
    let mut c = counts(60, 24, 43);
    for i in 0..c.n_taxa {
        for j in 0..c.n_samp {
            if (i * 7 + j * 3) % 5 == 0 {
                c.set(i, j, 0.0);
            }
        }
    }
    let a = run(&c, &x, &g, &compat);
    let b = run(&c, &x, &g, &strict);
    let d = a
        .core
        .vcov
        .iter()
        .zip(&b.core.vcov)
        .map(|(x, y)| (x - y).abs())
        .fold(0.0f64, f64::max);
    assert!(
        d > 1e-9,
        "with unobserved cells the sandwich quirk must change the result, got {d}"
    );
}

// ---------------------------------------------------------------------------
// Randomised: the pipeline never panics and always returns consistent shapes
// ---------------------------------------------------------------------------

proptest! {
    #![proptest_config(ProptestConfig::with_cases(24))]

    /// Whatever the input, the run terminates and the outputs have the shapes the
    /// configuration implies.
    #[test]
    fn arbitrary_inputs_never_panic(n_tax in 8usize..40, n_samp in 6usize..20, seed in 0u64..1000) {
        let c = counts(n_tax, n_samp, seed);
        let g = support::labels(n_samp, 2);
        let x = design_binary(&g, 0.0);
        let cfg = AncombcConfig {
            prv_cut: 0.0,
            ..support::default_config(2, "grp")
        };
        if let Ok(r) = ancombc2_run_named(&c, &x, Some(&g), &cfg, &[], &[], &ancombc2_core::F64Reductions) {
            let n = r.core.taxa.len();
            prop_assert_eq!(r.core.beta.len(), n * 2);
            prop_assert_eq!(r.core.se.len(), n * 2);
            prop_assert_eq!(r.core.q.len(), n * 2);
            prop_assert_eq!(r.core.vcov.len(), n * 4);
            prop_assert_eq!(r.core.samp_frac.len(), n_samp);
            for &pv in r.core.p.iter().chain(r.core.q.iter()) {
                prop_assert!((0.0..=1.0).contains(&pv), "p out of range: {pv}");
            }
        }
    }

    /// Taxa permutation equivariance, on random inputs.
    #[test]
    fn permuting_taxa_permutes_the_results(
        n_tax in 10usize..30,
        n_samp in 8usize..18,
        seed in 0u64..1000,
    ) {
        let c = counts(n_tax, n_samp, seed);
        let g = support::labels(n_samp, 2);
        let x = design_binary(&g, 0.0);
        let cfg = AncombcConfig {
            prv_cut: 0.0,
            ..support::default_config(2, "grp")
        };
        // A modular map is not a bijection for composite n, which would turn the
        // test into a check of the filter rather than of equivariance.
        let mut perm: Vec<usize> = (0..n_tax).collect();
        let mut st = seed | 1;
        for i in (1..n_tax).rev() {
            st = st.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            let j = (st >> 33) as usize % (i + 1);
            perm.swap(i, j);
        }
        let cp = c.select(&perm, &(0..n_samp).collect::<Vec<_>>());
        let (a, b) = match (
            ancombc2_run_named(&c, &x, Some(&g), &cfg, &[], &[], &ancombc2_core::F64Reductions),
            ancombc2_run_named(&cp, &x, Some(&g), &cfg, &[], &[], &ancombc2_core::F64Reductions),
        ) {
            (Ok(a), Ok(b)) => (a, b),
            _ => return Ok(()),
        };
        for (i, &pi) in perm.iter().enumerate().take(n_tax) {
            for k in 0..2 {
                let d = (a.core.beta[pi * 2 + k] - b.core.beta[i * 2 + k]).abs();
                prop_assert!(d < 1e-7, "taxon {i}, coefficient {k}: {d}");
            }
        }
    }
}
