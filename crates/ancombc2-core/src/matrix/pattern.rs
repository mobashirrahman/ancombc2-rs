//! Missing-value patterns: the grouping that makes the fixed-effects fit cheap.
//!
//! With `pseudo = 0` a zero count becomes `log 0 = -Inf`, which the oracle
//! replaces by `NA`. Each taxon is therefore observed on its own subset of
//! samples, and two taxa can only share a least-squares solve if their observed
//! samples coincide.
//!
//! `.lm_fit_all` in the reference groups taxa by exactly this pattern and
//! solves each group with one multi-response `lm.fit` — a `QR` per group, shared
//! by every taxon in the group and by every iteration of the bias loop. Grouping
//! with bitsets keeps the grouping step itself off the critical path: it is
//! `n_taxa * n_samples / 64` word operations plus a hash of the bit patterns,
//! rather than the oracle's `do.call(paste0, asplit(use * 1L, 2L))`, which
//! materialises a string per taxon.

use std::collections::HashMap;

/// A bitset over samples.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Bitset {
    words: Vec<u64>,
    n_bits: usize,
}

impl Bitset {
    pub fn zeros(n_bits: usize) -> Self {
        Self {
            words: vec![0u64; n_bits.div_ceil(64)],
            n_bits,
        }
    }

    #[inline]
    pub fn set(&mut self, i: usize) {
        self.words[i / 64] |= 1u64 << (i % 64);
    }

    #[inline]
    pub fn get(&self, i: usize) -> bool {
        self.words[i / 64] & (1u64 << (i % 64)) != 0
    }

    pub fn len(&self) -> usize {
        self.n_bits
    }

    pub fn is_empty(&self) -> bool {
        self.n_bits == 0
    }

    pub fn count_ones(&self) -> usize {
        self.words.iter().map(|w| w.count_ones() as usize).sum()
    }

    pub fn as_slice(&self) -> &[u64] {
        &self.words
    }
}

/// Taxa grouped by their pattern of observed samples.
///
/// The group order is the order of first appearance, matching
/// `split(seq_len(n_tax), factor(keys, levels = unique(keys)))` in the oracle,
/// which makes the accumulation order (and therefore the floating-point result)
/// identical.
#[derive(Debug, Clone)]
pub struct PatternGroups {
    /// Taxon indices per group, in ascending order.
    pub groups: Vec<Vec<usize>>,
    /// Observed-sample indices per group, in ascending order.
    pub rows: Vec<Vec<usize>>,
}

impl PatternGroups {
    pub fn n_groups(&self) -> usize {
        self.groups.len()
    }

    /// Total observed (taxon, sample) pairs; a proxy for the work the fit will do.
    pub fn nnz(&self) -> usize {
        self.groups
            .iter()
            .zip(self.rows.iter())
            .map(|(g, r)| g.len() * r.len())
            .sum()
    }
}

/// Group taxa by their observed-sample pattern.
///
/// `observed` is `n_taxa * n_samples`, row-major, where a non-finite entry marks
/// an unusable observation (R's `is.finite(Ymat) & x_ok`).
pub fn group_by_pattern(observed: &[bool], n_taxa: usize, n_samp: usize) -> PatternGroups {
    debug_assert_eq!(observed.len(), n_taxa * n_samp);
    let mut map: HashMap<Bitset, usize> = HashMap::new();
    let mut groups: Vec<Vec<usize>> = Vec::new();
    let mut rows: Vec<Vec<usize>> = Vec::new();

    for t in 0..n_taxa {
        let base = t * n_samp;
        let mut bits = Bitset::zeros(n_samp);
        for s in 0..n_samp {
            if observed[base + s] {
                bits.set(s);
            }
        }
        match map.get(&bits) {
            Some(&g) => groups[g].push(t),
            None => {
                let idx: Vec<usize> = (0..n_samp).filter(|&s| bits.get(s)).collect();
                map.insert(bits, groups.len());
                groups.push(vec![t]);
                rows.push(idx);
            }
        }
    }
    PatternGroups { groups, rows }
}

/// Fast path: when every taxon is observed on every sample there is exactly one
/// group and no bitset work is needed. The oracle takes the same branch
/// (`if (all(x_ok) && all(is.finite(Ymat)))`).
pub fn single_group(n_taxa: usize, n_samp: usize) -> PatternGroups {
    PatternGroups {
        groups: vec![(0..n_taxa).collect()],
        rows: vec![(0..n_samp).collect()],
    }
}

/// Decide whether the complete-observation fast path applies, and group.
pub fn group(observed: &[bool], n_taxa: usize, n_samp: usize) -> PatternGroups {
    if observed.iter().all(|&b| b) {
        return single_group(n_taxa, n_samp);
    }
    group_by_pattern(observed, n_taxa, n_samp)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn complete_data_is_a_single_group() {
        let obs = vec![true; 12];
        let g = group(&obs, 4, 3);
        assert_eq!(g.n_groups(), 1);
        assert_eq!(g.groups[0], vec![0, 1, 2, 3]);
        assert_eq!(g.rows[0], vec![0, 1, 2]);
    }

    #[test]
    fn groups_taxa_sharing_a_pattern() {
        // The example from the design notes:
        //   A 1111011111011
        //   B 1111011111011
        //   C 1101110111111
        //   D 1111011111011
        let pats = [
            vec![1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1],
            vec![1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1],
            vec![1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1],
            vec![1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1],
        ];
        let n_s = 13;
        let mut obs = vec![false; 4 * n_s];
        for (t, p) in pats.iter().enumerate() {
            for (s, &v) in p.iter().enumerate() {
                obs[t * n_s + s] = v != 0;
            }
        }
        let g = group(&obs, 4, n_s);
        assert_eq!(g.n_groups(), 2);
        assert_eq!(g.groups[0], vec![0, 1, 3], "first-appearance order");
        assert_eq!(g.groups[1], vec![2]);
        assert_eq!(g.rows[0].len(), 11);
        assert_eq!(g.rows[1].len(), 11);
        assert!(!g.rows[0].contains(&4));
        assert!(!g.rows[1].contains(&2));
    }

    #[test]
    fn nnz_counts_the_fit_work() {
        // 2 taxa x 3 samples: taxon 0 observed on {0,1}, taxon 1 on all three.
        let obs = vec![true, true, false, true, true, true];
        let g = group(&obs, 2, 3);
        assert_eq!(g.groups[0], vec![0]);
        assert_eq!(g.rows[0], vec![0, 1]);
        assert_eq!(g.groups[1], vec![1]);
        assert_eq!(g.rows[1], vec![0, 1, 2]);
        assert_eq!(g.nnz(), 2 + 3);
    }

    #[test]
    fn all_zero_rows_still_form_a_group() {
        // The oracle refits taxa with no usable sample one at a time; the group
        // must exist (with an empty row set) so the caller takes that path.
        let obs = vec![false, false, true, true];
        let g = group(&obs, 2, 2);
        assert_eq!(g.n_groups(), 2);
        assert!(g.rows[0].is_empty());
        assert_eq!(g.rows[1], vec![0, 1]);
    }

    #[test]
    fn bitset_basics() {
        let mut b = Bitset::zeros(130);
        assert_eq!(b.count_ones(), 0);
        b.set(0);
        b.set(64);
        b.set(129);
        assert_eq!(b.count_ones(), 3);
        assert!(b.get(0) && b.get(64) && b.get(129));
        assert!(!b.get(1));
        assert_eq!(b.len(), 130);
        assert_eq!(b.as_slice().len(), 3);
    }

    #[test]
    fn grouping_is_order_stable() {
        // 3 taxa x 2 samples: taxon 0 seen in both, taxon 1 only in sample 0,
        // taxon 2 only in sample 1. Group ids follow first appearance, so the
        // accumulation order of the sandwich accumulation matches the oracle's.
        let obs = vec![true, true, true, false, false, true];
        let g = group(&obs, 3, 2);
        assert_eq!(g.groups[0], vec![0]);
        assert_eq!(g.rows[0], vec![0, 1]);
        assert_eq!(g.groups[1], vec![1]);
        assert_eq!(g.rows[1], vec![0]);
        assert_eq!(g.groups[2], vec![2]);
        assert_eq!(g.rows[2], vec![1]);
    }

    #[test]
    fn large_input_is_consistent() {
        let n_tax = 500;
        let n_samp = 40;
        let mut obs = vec![true; n_tax * n_samp];
        for t in 0..n_tax {
            if t % 7 == 0 {
                obs[t * n_samp + (t % n_samp)] = false;
            }
        }
        let g = group(&obs, n_tax, n_samp);
        let total: usize = g.groups.iter().map(|v| v.len()).sum();
        assert_eq!(total, n_tax, "every taxon must appear in exactly one group");
        // reconstructing the patterns from the groups must reproduce the input
        let mut recon = vec![false; n_tax * n_samp];
        for (gi, taxa) in g.groups.iter().enumerate() {
            for &t in taxa {
                for &s in &g.rows[gi] {
                    recon[t * n_samp + s] = true;
                }
            }
        }
        assert_eq!(recon, obs);
    }
}
