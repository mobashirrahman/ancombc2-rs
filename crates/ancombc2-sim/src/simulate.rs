//! The microbiome count generator.
//!
//! # Model
//!
//! For taxon `i` and sample `j`:
//!
//! 1. A baseline log abundance `mu_i ~ Normal(0, abundance_sd)`, giving a
//!    log-normal spread of taxa across the community -- the "some taxa are
//!    common, some are rare" structure the abundance filter acts on.
//! 2. A group effect: `0` for a non-DA taxon, and `s_i * log_fc` for a DA
//!    taxon, where `s_i` is `+1` or `-1` with equal probability. Both
//!    directions are present so a method that only detects upregulation is
//!    caught.
//! 3. A library depth `L_j ~ LogNormal(mean = lib_mean, cv = lib_cv)`.
//! 4. Composition normalised over taxa:
//!    `p_ij = exp(mu_i + s_i log_fc 1{DA} 1{j in group 2}) / sum_k exp(...)`.
//! 5. A structural zero: with probability `zero_inflation` the cell is zero
//!    regardless of `p_ij`. These are the cells that create *missingness
//!    patterns*, so the grid exercises the pattern-grouped QR path rather than
//!    only the all-observed one. With `confound`, DA taxa are additionally
//!    zeroed at `zero_inflation + confounder` in group 1 only, so the
//!    per-sample sampling fraction correlates with group and with DA status.
//!    That is the confounding the sampling-fraction correction removes, and
//!    without it the correction is untested.
//! 6. `x_ij ~ NegBin(mean = L_j p_ij, size = dispersion)`, gamma-Poisson.
//!
//! Every step is a named function on a seeded [`Rng`], so a rep is a pure
//! function of `(grid seed, cell index, rep index)`.

use ancombc2_core::CountMatrix;
use serde::{Deserialize, Serialize};

use crate::grid::Cell;
use crate::rng::Rng;

/// One generated replicate: the data, the design, and the truth.
#[derive(Debug, Clone)]
pub struct Replicate {
    pub counts: CountMatrix,
    /// Group index per sample, 0 or 1. Group 1 is the one with the effect.
    pub group: Vec<usize>,
    pub taxon_names: Vec<String>,
    pub sample_names: Vec<String>,
    /// Per-taxon ground truth.
    pub truth: Vec<TaxonTruth>,
}

/// Ground truth for one taxon.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TaxonTruth {
    pub name: String,
    /// The taxon was *labelled* differentially abundant by the generator.
    ///
    /// This is not the same as having a compositional effect. In a confounded
    /// cell the labelled taxa have no effect on their expected abundance at all
    /// -- they differ only in how often they are observed to zero -- so `is_da`
    /// marks the label and [`TaxonTruth::has_effect`] marks the truth. The two
    /// are equal in an unconfounded cell.
    pub is_da: bool,
    /// The taxon's expected abundance actually differs between the groups.
    ///
    /// False for every labelled taxon in a confounded cell, which is what makes
    /// that arm a negative control: a method that ignores the sampling fraction
    /// will call them differentially abundant, and the correction should not.
    pub has_effect: bool,
    /// The label came from a group-dependent zero rate rather than a
    /// compositional effect, so these taxa are false positives for a method that
    /// does not correct the sampling fraction.
    pub confounded: bool,
    /// Sign of the true effect, `+1` or `-1`, or `0` for a null taxon.
    pub direction: i32,
    /// The true log fold change of the *expected count ratio* between the two
    /// groups. This is the estimand ANCOM-BC2 targets, and it is what `beta` is
    /// compared against -- not the realised log ratio of the sampled counts,
    /// which carries sampling noise that the estimator is not expected to
    /// reproduce taxon by taxon.
    pub log_fc: f64,
    /// The realised fraction of zeros over all samples: the per-sample
    /// analogue is what ANCOM-BC2 reports as the sampling fraction.
    pub zero_fraction: f64,
    /// The realised fraction of zeros in the control group.
    pub zero_fraction_control: f64,
    /// The realised fraction of zeros in the treated group.
    pub zero_fraction_treated: f64,
}

/// Extra factors that are not grid axes.
#[derive(Debug, Clone, Copy)]
pub struct SimControl {
    pub abundance_sd: f64,
    pub dispersion: f64,
    /// Extra zero probability applied to DA taxa in group 1 when
    /// `confound` is set.
    pub confounder: f64,
}

pub fn simulate(cell: &Cell, seed: u64, rep: usize, ctl: &SimControl) -> Replicate {
    // Two streams. The community is fixed for the whole cell; the sampling is
    // redrawn per replicate. See [`Rng::for_cell`] for why.
    let mut comm = Rng::for_cell(seed, cell.index);
    let mut rng = Rng::for_rep(seed, cell.index, rep);
    let n_taxa = cell.n_taxa;
    let n_samp = cell.n_samp;

    // Group 1 is the treated group. The split is as balanced as possible with
    // the odd sample going to the control, so the two groups differ in size by
    // at most one and the truth is never trivially a group-size artefact.
    let n_g1 = n_samp / 2;
    let group: Vec<usize> = (0..n_samp).map(|j| if j < n_g1 { 0 } else { 1 }).collect();

    // Which taxa are DA, and in which direction.
    let n_da = ((cell.da_proportion * n_taxa as f64).round() as usize).min(n_taxa);
    let mut is_da = vec![false; n_taxa];
    for k in 0..n_da {
        is_da[k] = true;
    }
    // Shuffle so DA taxa are not the first `n_da` rows, which would correlate
    // the truth with the taxon's position in the file and with the bitset
    // pattern numbering. From the *cell* stream, so every replicate of this
    // cell has the same taxa in the same roles.
    for k in (1..n_taxa).rev() {
        let j = (comm.next_u64() % (k as u64 + 1)) as usize;
        is_da.swap(k, j);
    }
    let direction: Vec<i32> = is_da
        .iter()
        .map(|d| {
            if !*d {
                0
            } else if comm.bernoulli(0.5) {
                1
            } else {
                -1
            }
        })
        .collect();

    // Baseline log abundances, also fixed for the cell: the same taxon has the
    // same expected relative abundance in every replicate.
    let mu: Vec<f64> = (0..n_taxa)
        .map(|_| comm.normal_sd(0.0, ctl.abundance_sd))
        .collect();

    // Library depths.
    let lib: Vec<f64> = (0..n_samp)
        .map(|_| rng.lognormal_mean_cv(cell.lib_mean, cell.lib_cv))
        .collect();

    // Structural-zero mask, drawn before the counts so the zero fraction is
    // known exactly rather than inferred from a negative binomial's own zeros.
    let mut zero_mask = vec![false; n_taxa * n_samp];
    for i in 0..n_taxa {
        for j in 0..n_samp {
            let mut p = cell.zero_inflation;
            if cell.confound && is_da[i] && group[j] == 0 {
                p += ctl.confounder;
            }
            zero_mask[i * n_samp + j] = p > 0.0 && rng.bernoulli(p.min(1.0));
        }
    }

    // Expected counts, normalised across taxa within a sample.
    let mut data = vec![0.0f64; n_taxa * n_samp];
    for j in 0..n_samp {
        let mut shifted: Vec<f64> = Vec::with_capacity(n_taxa);
        let mut total = 0.0f64;
        for i in 0..n_taxa {
            // A confounded cell withholds the compositional effect entirely;
            // see the module docs. The confounded taxa are therefore nulls that
            // differ only in detection rate.
            let effect = if is_da[i] && !cell.confound && group[j] == 1 {
                cell.log_fc * direction[i] as f64
            } else {
                0.0
            };
            let e = (mu[i] + effect).exp();
            shifted.push(e);
            total += e;
        }
        for i in 0..n_taxa {
            let p_ij = shifted[i] / total;
            let lambda = lib[j] * p_ij;
            if zero_mask[i * n_samp + j] {
                data[i * n_samp + j] = 0.0;
                continue;
            }
            // Negative binomial by gamma-Poisson. `size` is the NB's `r`, so
            // the gamma shape is `r` and its scale is `mu / r`.
            let g = rng.gamma_any(ctl.dispersion, lambda / ctl.dispersion);
            data[i * n_samp + j] = rng.poisson(g) as f64;
        }
    }

    let taxon_names: Vec<String> = (0..n_taxa).map(|i| format!("taxon{i:05}")).collect();
    let sample_names: Vec<String> = (0..n_samp).map(|j| format!("sample{j:04}")).collect();

    let n_g1 = n_g1.max(1);
    let n_g2 = (n_samp - n_g1).max(1);
    let truth: Vec<TaxonTruth> = (0..n_taxa)
        .map(|i| {
            let frac = |pred: &dyn Fn(usize) -> bool| {
                (0..n_samp)
                    .filter(|j| pred(*j) && data[i * n_samp + j] == 0.0)
                    .count() as f64
            };
            let zeros = (0..n_samp).filter(|j| data[i * n_samp + j] == 0.0).count();
            TaxonTruth {
                name: taxon_names[i].clone(),
                is_da: is_da[i],
                direction: direction[i],
                log_fc: if is_da[i] && !cell.confound {
                    cell.log_fc * direction[i] as f64
                } else {
                    0.0
                },
                has_effect: is_da[i] && !cell.confound,
                confounded: is_da[i] && cell.confound,
                zero_fraction: zeros as f64 / n_samp as f64,
                zero_fraction_control: frac(&|j| group[j] == 0) / n_g1 as f64,
                zero_fraction_treated: frac(&|j| group[j] == 1) / n_g2 as f64,
            }
        })
        .collect();

    Replicate {
        // A count of zero is a legitimate observed zero, so the table is
        // complete here; NA handling is a separate concern covered by the
        // golden and edge-case fixtures.
        counts: CountMatrix::new(n_taxa, n_samp, data).expect("shape matches"),
        group,
        taxon_names,
        sample_names,
        truth,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::grid::Cell;

    fn cell() -> Cell {
        Cell {
            index: 0,
            n_taxa: 200,
            n_samp: 20,
            da_proportion: 0.3,
            log_fc: 1.0,
            zero_inflation: 0.2,
            lib_mean: 1e4,
            lib_cv: 0.3,
            confound: false,
        }
    }

    fn ctl() -> SimControl {
        SimControl {
            abundance_sd: 0.5,
            dispersion: 1.0,
            confounder: 0.3,
        }
    }

    /// The community is a property of the cell, not the replicate: the same
    /// taxa, the same DA status, the same directions, the same baseline
    /// abundances, in every replicate.
    ///
    /// This is what makes the per-taxon metrics valid. If the community were
    /// redrawn per replicate, `taxon00005` in rep 0 and `taxon00005` in rep 1
    /// would be unrelated taxa, and the SE calibration ratio would divide one
    /// taxon's standard error by the spread of a different taxon's estimates.
    #[test]
    fn the_community_is_fixed_across_a_cell_and_the_sampling_is_not() {
        let a = simulate(&cell(), 42, 0, &ctl());
        let b = simulate(&cell(), 42, 1, &ctl());
        let c = simulate(&cell(), 42, 2, &ctl());
        for r in [&a, &b, &c] {
            assert_eq!(
                r.truth.iter().map(|t| t.is_da).collect::<Vec<_>>(),
                a.truth.iter().map(|t| t.is_da).collect::<Vec<_>>()
            );
            assert_eq!(
                r.truth.iter().map(|t| t.log_fc).collect::<Vec<_>>(),
                a.truth.iter().map(|t| t.log_fc).collect::<Vec<_>>()
            );
        }
        // And the counts do differ, or the replicates would be pointless.
        assert_ne!(a.counts.data, b.counts.data);
        assert_ne!(b.counts.data, c.counts.data);
    }

    /// Two different cells must not share a community.
    #[test]
    fn a_different_cell_gets_a_different_community() {
        let mut other = cell();
        other.index = 1;
        let a = simulate(&cell(), 42, 0, &ctl());
        let b = simulate(&other, 42, 0, &ctl());
        assert_ne!(
            a.truth.iter().map(|t| t.log_fc).collect::<Vec<_>>(),
            b.truth.iter().map(|t| t.log_fc).collect::<Vec<_>>()
        );
    }

    #[test]
    fn a_rep_is_a_pure_function_of_its_seed_cell_and_rep() {
        let a = simulate(&cell(), 42, 3, &ctl());
        let b = simulate(&cell(), 42, 3, &ctl());
        assert_eq!(a.counts.data, b.counts.data);
        let c = simulate(&cell(), 42, 4, &ctl());
        assert_ne!(a.counts.data, c.counts.data, "reps must differ");
        let d = simulate(&cell(), 43, 3, &ctl());
        assert_ne!(a.counts.data, d.counts.data, "the run seed must matter");
    }

    #[test]
    fn the_da_taxa_are_not_the_first_rows() {
        let r = simulate(&cell(), 7, 0, &ctl());
        let first_half_da = r.truth.iter().take(100).filter(|t| t.is_da).count();
        assert!(
            first_half_da < 60,
            "DA taxa should be spread through the table, got {first_half_da} in the first half"
        );
    }

    #[test]
    fn both_directions_occur_and_match_the_truth() {
        let r = simulate(&cell(), 11, 0, &ctl());
        let ups = r.truth.iter().filter(|t| t.log_fc > 0.0).count();
        let downs = r.truth.iter().filter(|t| t.log_fc < 0.0).count();
        assert!(ups > 0 && downs > 0, "both directions must be present");
        for t in &r.truth {
            assert_eq!(
                t.is_da,
                t.log_fc != 0.0,
                "is_da and log_fc must agree for {}",
                t.name
            );
        }
    }

    #[test]
    fn the_da_proportion_is_respected() {
        let mut c = cell();
        c.da_proportion = 0.3;
        c.n_taxa = 1000;
        let r = simulate(&c, 5, 0, &ctl());
        let n_da = r.truth.iter().filter(|t| t.is_da).count();
        assert!(n_da.abs_diff(300) <= 1, "got {n_da} DA taxa");
    }

    #[test]
    fn zero_inflation_raises_the_observed_zero_fraction() {
        let mut c = cell();
        c.zero_inflation = 0.0;
        let low = simulate(&c, 9, 0, &ctl());
        c.zero_inflation = 0.6;
        let high = simulate(&c, 9, 0, &ctl());
        let mean = |r: &Replicate| {
            r.truth.iter().map(|t| t.zero_fraction).sum::<f64>() / r.truth.len() as f64
        };
        assert!(
            mean(&high) > mean(&low) + 0.4,
            "zero inflation must dominate: {} vs {}",
            mean(&high),
            mean(&low)
        );
    }

    #[test]
    fn the_count_of_ones_factor_is_visibly_overdispersed() {
        // A negative binomial with size 1 has far more near-zero cells than a
        // Poisson of the same mean. This is what makes the sampler a real
        // microbiome simulator rather than a Poisson with extra steps.
        let mut c = cell();
        c.zero_inflation = 0.0;
        c.lib_cv = 0.0;
        let r = simulate(&c, 13, 0, &ctl());
        let zeros = r.counts.data.iter().filter(|v| **v == 0.0).count();
        let frac = zeros as f64 / r.counts.data.len() as f64;
        assert!(frac > 0.02, "expected visible zeros, got {frac}");
    }

    #[test]
    fn confounding_ties_the_sampling_fraction_to_the_group() {
        // The premise of the confounded arm, stated as a property of the
        // generated data: DA taxa must be markedly sparser in the control
        // group than in the treated group, and much more so than in the
        // unconfounded arm. If this does not hold, the confounded cells are not
        // actually confounded and the arm measures nothing.
        let control_deficit = |confound: bool| {
            let mut c = cell();
            c.confound = confound;
            c.zero_inflation = 0.1;
            let r = simulate(&c, 17, 0, &ctl());
            let da: Vec<&TaxonTruth> = r.truth.iter().filter(|t| t.is_da).collect();
            // Control minus treated: positive means the control is sparser,
            // which is the direction the confounder acts in. It is not as large
            // as the confounder itself, because the treated group is also the
            // one whose DA taxa carry the effect and so have higher expected
            // counts and fewer incidental zeros.
            da.iter()
                .map(|t| t.zero_fraction_control - t.zero_fraction_treated)
                .sum::<f64>()
                / da.len() as f64
        };
        let plain = control_deficit(false);
        let confounded = control_deficit(true);
        assert!(
            confounded > plain + 0.1,
            "the confounded control-group deficit {confounded} must exceed the plain {plain}"
        );
    }

    #[test]
    fn library_sizes_follow_the_grid() {
        let mut c = cell();
        c.lib_mean = 1e4;
        c.lib_cv = 0.6;
        c.zero_inflation = 0.5;
        let r = simulate(&c, 19, 0, &ctl());
        // With 50% of cells zeroed the depth is only recoverable up to that
        // factor, so the check is a loose bracket.
        let depth: f64 = r.counts.data.iter().sum();
        let per_sample = depth / r.counts.n_samp as f64;
        assert!(
            per_sample > 0.25 * c.lib_mean && per_sample < 0.85 * c.lib_mean,
            "depth per sample {per_sample} is not near lib_mean {}",
            c.lib_mean
        );
    }
}
