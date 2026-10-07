//! The iterative MLE: `.iter_mle` and its two supporting routines.
//!
//! The estimator alternates between a per-taxon least-squares fit of
//! `Y - theta` on the design `X`, and an update of the sample-specific bias term
//! `theta` given the fitted values:
//!
//! ```text
//! theta <- 0
//! repeat {
//!     beta  <- fit(Y - theta, X)
//!     theta <- colMeans(Y - fitted, na.rm = TRUE)
//! } until  sqrt(||dbeta||^2 + ||dtheta||^2) <= tol   or  iter == max_iter
//! ```
//!
//! Two structural properties carry the performance, and both are already present
//! in the reference (which groups by missingness pattern and reuses one QR per
//! group):
//!
//! 1. `X` does not change across iterations, so each pattern's QR is factored
//!    **once** and reused for every iteration, every taxon in the group, and
//!    every sensitivity refit.
//! 2. Only the right-hand side changes between iterations, so each iteration is
//!    a triangular solve against a cached factorisation, not a refit.
//!
//! # Residuals and the sandwich input
//!
//! `eps = (Y - fitted) - theta`, i.e. the reference's `t(t(y - y_crt_hat) - theta)`.
//! The two subtractions are in that order and the sign matters for the
//! compatibility mode, because the `0.1` substitution in `.sandwich_vcov` keys
//! off `is.finite(eps)`.
//!
//! # Rank deficiency
//!
//! `.lm_fit_all` detects a rank-deficient sub-design (`fit$rank < ncol(xr)`) and
//! falls back to fitting each taxon individually with `lm()`, which drops an
//! absent factor level instead of failing. The unidentifiable column is reported
//! as `NA` in `beta` and the taxon keeps `dof = 999`. A taxon with no usable
//! sample at all is also reported unfitted. Both are represented here as
//! `beta = NaN`, matching the reference's `NA`.

use crate::config::IterControl;
use crate::error::{AncombcError, Result};
use crate::matrix::{dqrls_multi_selected, Matrix, LM_FIT_TOL};
use crate::vcov::{sandwich_all, OuterProducts};

/// Cached per-pattern factorisations, built once and reused across iterations.
///
/// The design sub-matrices are retained alongside the factors: the fitted values
/// are `X_sub * beta`, and recomputing `X_sub * beta` costs the same as one
/// `Q` pass, so there is no reason to keep the factors and discard the input.
/// Retaining it also makes the group membership auditable, which matters when a
/// parity failure has to be localised to a pattern.
#[derive(Debug, Clone)]
pub struct DesignCache {
    /// The full design, retained so each pattern's sub-design and factorisation
    /// can be built on demand.
    ///
    /// Holding the pattern's own `designs` and `qrs` arrays instead was
    /// a `O(n_patterns * n_rows * p)` structure: on `bm5` that is ~1000 patterns
    /// at ~18000 rows each, and the two arrays came to **5.25 GB** -- measured,
    /// with the design and its reflectors alone. That is more than the entire
    /// count matrix times thirty, on a run that is supposed to be the memory
    /// reference.
    ///
    /// It was also a poor trade on time. Building a pattern's factorisation costs
    /// `O(n_rows p^2)` while applying it to one taxon costs `O(n_rows p)`, so
    /// caching only pays when a pattern carries more than about `p` taxa. The
    /// real tables do not: `bm4` and `bm5` each have roughly one pattern per
    /// taxon, because a taxon's missing counts give it its own pattern. The cache
    /// was paying `p` times the work to reuse nothing.
    ///
    /// So the per-pattern factors are built where they are used -- inside the
    /// group loop, which is already parallel over patterns -- and freed when that
    /// group finishes. Peak is then `n_threads` patterns' worth rather than
    /// `n_patterns`, and the factorisation is built once per `lm_fit_all` pass
    /// instead of once per MLE run. [`Self::build_pattern`] is the single place
    /// that knows how.
    pub x: Matrix,
    pub rows: Vec<Vec<usize>>,
    pub taxa: Vec<Vec<usize>>,
    pub p: usize,
    pub n_samp: usize,
    /// Indices of the treatment-contrast columns for the group factor.
    ///
    /// Needed only by the rank-deficient path, and only because of what the
    /// reference does there: it refits each taxon of a rank-deficient group with
    /// `stats::lm`, and `lm` *fails* when a factor has fewer than two levels
    /// among the rows it actually uses. Knowing which columns are group
    /// contrasts is what lets that failure be reproduced without a formula
    /// parser. Empty when the caller does not know, in which case the per-taxon
    /// level test is skipped and the older group-rank rule applies.
    pub group_cols: Vec<usize>,
}

impl DesignCache {
    /// Group taxa by observed-sample pattern.
    ///
    /// Only the grouping is done here; the per-pattern designs and
    /// factorisations are built by [`Self::build_pattern`] where they are used.
    pub fn build(
        observed: &[bool],
        n_taxa: usize,
        n_samp: usize,
        x: &Matrix,
        group_cols: &[usize],
    ) -> Self {
        let groups = crate::matrix::group_by_observation(observed, n_taxa, n_samp);
        DesignCache {
            x: x.clone(),
            rows: groups.rows,
            taxa: groups.groups,
            p: x.cols,
            n_samp,
            group_cols: group_cols.to_vec(),
        }
    }

    /// The sub-design for pattern `g`.
    ///
    /// Name-free: the fit indexes the sub-design by row number, and one `String`
    /// per row per pattern was the run's largest allocation source long before the
    /// retention was questioned. See [`Matrix::select_rows_anon`]. The
    /// factorisation is R's `dqrdc2`, run by [`dqrls_multi_selected`] where the
    /// responses are known.
    pub fn build_pattern(&self, g: usize) -> Matrix {
        self.x.select_rows_anon(&self.rows[g])
    }
}

/// Result of one `.lm_fit_all` pass: coefficients, fitted values, residual
/// degrees of freedom.
#[derive(Debug, Clone)]
pub struct FitAll {
    /// `n_taxa * p` row-major; `NaN` for an unfitted taxon.
    pub beta: Vec<f64>,
    /// `n_taxa * n_samp` row-major; `0` for a sample excluded from the fit.
    pub fitted: Vec<f64>,
    /// Residual degrees of freedom per taxon; `999` when not estimable.
    pub dof: Vec<f64>,
}

/// One `.lm_fit_all` call: fit `X` to every taxon's theta-adjusted response.
///
/// `y` is `n_taxa x n_samp` row-major, already theta-adjusted.
/// One missingness group's contribution to a `lm_fit_all` pass.
///
/// Returned by value rather than written through `&mut`, for two reasons. It is
/// what lets the group loop be a `par_iter`: the writes become disjoint by
/// construction rather than by an argument about which taxa each group owns. And
/// it gives each worker its own scratch instead of sharing one `used` buffer,
/// which would otherwise be the only thing in the loop two threads could race on.
#[derive(Debug, Clone, Default)]
struct GroupWrites {
    /// `(taxon, coefficients)`, row-major within a taxon.
    beta: Vec<(usize, Vec<f64>)>,
    /// `(taxon, fitted)` where the fitted values are **parallel to the group's
    /// `rows`**, so a taxon costs one `f64` per used sample and nothing else.
    ///
    /// This was `Vec<(usize, f64)>` -- a `(sample, value)` pair per entry -- and
    /// that `usize` is pure redundancy: every taxon of a group is fitted on the
    /// same row set, so the sample index is already known from `cache.rows[g]`.
    /// Storing it per value doubled the largest structure the group loop retains
    /// until its `collect` returns. Measured on `bm5` that array was **1.4 GB**,
    /// and the conservative sensitivity analysis runs three of these pipelines at
    /// once, so it was paid three times over.
    fitted: Vec<(usize, Vec<f64>)>,
    /// `(taxon, residual degrees of freedom)`.
    dof: Vec<(usize, f64)>,
}

/// Fit every taxon of one missingness group.
///
/// `g` indexes `cache`'s per-pattern arrays. An empty `rows` yields nothing, which
/// leaves the caller-initialised `NA`/`999` in place -- `.lm_fit_all`'s
/// "no usable sample" branch, where `fit_one` fails and the taxon stays at NA.
fn fit_one_group(cache: &DesignCache, y: &[f64], n_samp: usize, g: usize) -> GroupWrites {
    let p = cache.p;
    let taxa = &cache.taxa[g];
    let rows = &cache.rows[g];
    let mut out = GroupWrites::default();
    if rows.is_empty() {
        // No usable sample: `.lm_fit_all` falls back to `fit_one`, whose `lm`
        // fails for a single-level factor and leaves the taxon at NA with dof 999.
        return out;
    }
    // Built here, not cached: see `DesignCache::x`. The factorisation is
    // `O(n_rows p^2)` against `O(n_rows p)` to apply it, so for the one-taxon
    // patterns that dominate real tables this is a win on time as well as memory.
    let xsub = cache.build_pattern(g);
    let xsub = &xsub;
    // `stats::lm.fit(xr, Yr)`: R's `dqrls`, one factorisation of the design and one
    // solve per response column. Run before the rank test because the reference's
    // rank test *is* this fit's rank.
    let yr: Vec<Vec<f64>> = taxa
        .iter()
        .map(|&t| rows.iter().map(|&s| y[t * n_samp + s]).collect())
        .collect();
    let ycols: Vec<&[f64]> = yr.iter().map(Vec::as_slice).collect();
    let fit = dqrls_multi_selected(xsub.data.clone(), rows.len(), p, &ycols, LM_FIT_TOL);
    if fit.rank < p {
        // Rank-deficient sub-design. The reference does not pseudo-solve; it
        // refits every taxon of the group on its own with `lm`, and `lm`
        // *drops* an aliased column from `coef()`. `.lm_fit_all` then writes
        // the surviving names into a zero-initialised row, so an aliased
        // coefficient is a literal 0 -- not NA, and not the 1e14 value a
        // pseudo-inverse would produce. Reproduced here taxon by taxon,
        // because a per-taxon `lm` can have a different rank from the group.
        for &t in taxa {
            // `fit_one`'s `lm` can still fail -- see
            // `per_taxon_lm_would_succeed` -- and a failed taxon keeps NA for both
            // its coefficients and its fitted values.
            if !per_taxon_lm_would_succeed(cache, xsub, rows) {
                continue;
            }
            let resp: Vec<f64> = rows.iter().map(|&s| y[t * n_samp + s]).collect();
            // `fit_one` does **not** factorise the group's design. It builds
            // `df = data.frame(y_crt = Ymat[i, ], meta_data)` -- all `n_samp`
            // rows, with `NA` where this taxon's response is missing -- and calls
            // `lm(tformula, data = df)`. `lm` drops the incomplete rows first,
            // then `model.matrix` applies `drop.unused.levels = TRUE` to what is
            // left, so **the design it factorises is not `xsub`**: every group
            // contrast for a level this taxon never observed is gone.
            //
            // Measured on the oracle, for the first retained taxon of the
            // `int-sparsity90-5group` cell (3 usable samples, 6 design columns,
            // levels {1, 2, 4} of 5):
            //
            // ```text
            // > colnames(model.matrix(lm(y_crt ~ group + x1, df)))
            // [1] "(Intercept)" "groupg2"  "groupg4"  "x1"
            // > fit_one writes
            // (Intercept)=1.40771  group2=-1.30597  group3=0  group4=-2.86441
            // group5=0  x1=NA
            // ```
            //
            // `group3` and `group5` are literal `0` because they have no name in
            // `coef()` and `fit_one` writes into `rep(0, p)`; `x1` is `NA`
            // because it *is* named and `coef.lm` marks it aliased.
            //
            // This has to be a real factorisation of the reduced design rather
            // than the full one with its unobserved columns patched to zero. The
            // two differ in which columns form the basis: the full design's
            // columns are `[1,1,1]`, `[0,0,1]`, `[0,0,0]`, `[0,1,0]`, `[0,0,0]`,
            // `[·]`, so an unpivoted greedy QR admits 1, 2 and 4 and rejects 3,
            // 5 and 6 -- an *interior* rejection, which the `aliased = i >= rank`
            // rule in [`qr`] cannot express. That rule is correct for
            // `lm.fit`'s contiguous leading run and wrong here, so patching after
            // the fact gave `beta_star` 166% wrong on this cell.
            //
            // `xsub` is the group's sub-design: its rows are `rows` in order, so
            // it is indexed by position within `rows`, never by the global sample
            // number -- which reads past the end for any group that is not the
            // first block of samples.
            //
            // A treatment contrast is 1 exactly on the rows of its own level, so
            // "the taxon observed that level" is "the column is ever 1". Dropping
            // the columns that are not is then just filtering, and the surviving
            // contrasts keep their original names and values: R re-bases against
            // the first *observed* level, which is the same 0/1 coding shifted,
            // not a new column.
            //
            // The base level is the first **observed** one, not the first globally.
            // `model.matrix` re-levels the factor to what the na.omit-reduced frame
            // contains and then applies `contr.treatment`, which drops the first
            // level it sees. So a taxon observed at levels {3, 4, 5} gets columns
            // for 4 and 5 only: `group3` is the base and gets no column at all,
            // even though `group3` is 1 on some of its rows and is therefore
            // "observed" in every other sense. Measured against the oracle's own
            // `beta_star` on `int-sparsity90-5group`, keeping the smallest
            // observed level's column instead flipped 20 of the 64 rows -- `F`
            // where the oracle has `0`, and `0` where it has `F`.
            let g0 = cache.group_cols.first().copied();
            // A treatment contrast is 1 exactly on its own level's rows, so
            // "observed" is "ever 1". `group_cols` is ascending in column order
            // and `model.matrix` emits the contrasts in level order, so the
            // contrast for level L sits at `g0 + (L - 2)` -- which means the level
            // has to be recovered from the column index, not compared with it:
            // the *first* level has no contrast at all, so the smallest contrast
            // index is the second level, not the base.
            let level_of = |j: usize| g0.map(|g| 2 + j - g);
            //
            // The level has to come from the *rows*, because the base level is the
            // one with no contrast of its own: a row with every group contrast 0
            // is level 1, and level 1 being the base is exactly the case that
            // makes `group2` the first column rather than the dropped one. Reading
            // the base off the contrast indices instead -- "the lowest contrast
            // index observed is the base" -- drops `group2` for a taxon observed
            // at levels {1, 2, 4}, which is 20 of the 64 rows on
            // `int-sparsity90-5group`.
            let row_level = |li: usize| {
                cache
                    .group_cols
                    .iter()
                    .copied()
                    .find(|&j| xsub.get(li, j) != 0.0)
                    .and_then(level_of)
                    .unwrap_or(1)
            };
            let base = (0..rows.len()).map(&row_level).min();
            let cols: Vec<usize> = (0..p)
                .filter(|&j| match level_of(j) {
                    Some(l) if cache.group_cols.contains(&j) => {
                        l != base.unwrap_or(1) && (0..rows.len()).any(|li| row_level(li) == l)
                    }
                    _ => !cache.group_cols.contains(&j),
                })
                .collect();
            let red = xsub.select_cols(&cols);
            // `stats::lm` -> `lm.fit`: the same `dqrls`, on the reduced design.
            let red_fit = dqrls_multi_selected(
                red.data.clone(),
                rows.len(),
                cols.len(),
                &[resp.as_slice()],
                LM_FIT_TOL,
            );
            if red_fit.rank == 0 {
                continue;
            }
            // `fit_one`'s `bi = rep(0, p)`: a column with no name in `coef()` keeps
            // the initialiser, so it is a literal 0 rather than `NA`. A named but
            // aliased column is `NA` (`lm.fit` fills those slots after unpivoting).
            let mut coef = vec![0.0f64; p];
            for (pos, &pc) in red_fit.pivot.iter().enumerate() {
                coef[cols[pc]] = if pos < red_fit.rank {
                    red_fit.coef[0][pos]
                } else {
                    f64::NAN
                };
            }
            // `stats::fitted` is `lm.fit`'s *projection*, not `X %*% coef(fit)`.
            // That distinction is invisible for a full-rank fit and decisive for a
            // rank-deficient one: `coef` carries `NA` for every aliased term, so
            // the product would be `NA` along the whole row, and `theta` is
            // `colMeans(y - fitted, na.rm = TRUE)` -- an `NA` there does not merely
            // lose one term, it drops the taxon out of that sample's mean entirely.
            //
            // Measured on `int-sparsity90-5group`: propagating the `NA` left only
            // about a fifth of each sample's observed taxa contributing to `theta`,
            // three samples with *no* contributor at all (so `theta` came back
            // `NaN`, which then made the whole adjusted-response column `NaN` on
            // the next iteration and cascaded), and the first iteration's epsilon
            // came out 2.23 where the oracle's is 0.998.
            //
            // Zeroing the aliased coefficients is not an approximation of the
            // projection, it *is* the projection: an aliased column is by
            // definition in the span of the ones that were kept, so dropping it
            // leaves the fitted space unchanged and `X %*% coef_kept` is the unique
            // least-squares fit in that space. What `coef()` reports as `NA` is
            // only the fact that the coordinate is not identified; the fitted
            // values are.
            // And it is computed the way `lm.fit` computes it: `y - residuals`,
            // not `X %*% coef`, which differ in the low bits.
            let fitted_row: Vec<f64> = resp
                .iter()
                .zip(&red_fit.resid[0])
                .map(|(yv, r)| yv - r)
                .collect();
            out.fitted.push((t, fitted_row));
            out.beta.push((t, coef));
            // `fit$df.residual` is measured on the *reduced* model, so the rank
            // that goes with it is the reduced rank.
            out.dof.push((t, (rows.len() - red_fit.rank) as f64));
        }
        return out;
    }
    // Full rank, so `dqrdc2` moved nothing and `pivot` is the identity: the
    // coefficients are already in design-column order.
    debug_assert!(fit.pivot.iter().enumerate().all(|(i, &j)| i == j));
    for (k, &t) in taxa.iter().enumerate() {
        // `fitted.values <- y - z$residuals`
        let fitted_row: Vec<f64> = yr[k]
            .iter()
            .zip(&fit.resid[k])
            .map(|(yv, r)| yv - r)
            .collect();
        out.fitted.push((t, fitted_row));
        out.beta.push((t, fit.coef[k].clone()));
    }
    // `df.residual = n - rank`
    let df = (rows.len() - fit.rank) as f64;
    out.dof.extend(taxa.iter().map(|&t| (t, df)));
    out
}

/// One `.lm_fit_all` call: fit `X` to every taxon's theta-adjusted response.
///
/// `y` is `n_taxa x n_samp` row-major, already theta-adjusted.
///
/// # Level 3 of the nesting order
///
/// The missingness groups are fitted in parallel. Each group writes only its own
/// taxa's slots of `beta`, `fitted` and `dof`, and a taxon belongs to exactly one
/// group, so the writes are disjoint by construction: **there is no reduction
/// here whose order is observable**, and splitting the pool cannot change a bit of
/// the result. That is the whole reason this level can be parallel while the
/// MLE's *iterations* stay sequential -- `theta` depends on every taxon's fitted
/// values, so iteration `k+1` is genuinely a barrier, and only the work within one
/// iteration is shared out.
///
/// The levels above it (pseudo-count runs, E-M coefficients) take the pool instead
/// when they are present, because [`crate::parallel::NestingBudget`] hands it to exactly one level.
/// See `parallel.rs`.
pub fn lm_fit_all(cache: &DesignCache, y: &[f64], n_taxa: usize, n_samp: usize) -> FitAll {
    let p = cache.p;
    let mut beta = vec![f64::NAN; n_taxa * p];
    // NaN-initialised, matching the reference: `.lm_fit_all` builds
    // `fitted = matrix(NA_real_, n_tax, n_samp)` and only writes what a fit
    // produced. A taxon whose per-taxon `lm` *fails* therefore keeps NA across
    // its whole row, and that is not cosmetic -- `theta` is
    // `colMeans(y - fitted, na.rm = TRUE)`, so NA drops the cell from the mean
    // while a 0 would feed `y` itself into it. On the 90%-zero fixture-matrix
    // cell that single difference moved `theta` by 0.75 while `y1` still agreed
    // to 9e-16.
    //
    // Both successful paths in `fit_one_group` write the zeros explicitly,
    // because the reference writes them explicitly too: `fitted[idx, !rows] = 0`
    // in the grouped branch, and `fit_one`'s `fi = rep(0, n_samp)` over the whole
    // row.
    let mut fitted = vec![f64::NAN; n_taxa * n_samp];
    let mut dof = vec![999.0; n_taxa];

    let mut lvl = crate::parallel::NestingBudget::level("missingness groups");
    let indices: Vec<usize> = (0..cache.taxa.len()).collect();
    let results = crate::parallel::map_par(&mut lvl, &indices, |_k, &g| {
        fit_one_group(cache, y, n_samp, g)
    });
    // `map_par` preserves index order, so zipping against `cache.rows` recovers
    // which pattern each result belongs to.
    for (rows, g) in cache.rows.iter().zip(results) {
        for (t, coef) in g.beta {
            beta[t * p..t * p + p].copy_from_slice(&coef);
        }
        // `fitted[idx, !rows] = 0`: a successful fit writes a literal 0 at every
        // sample the group's fit did not use, and the buffer is NaN-initialised,
        // so that zero is an explicit write and not something the initialiser
        // could be left to do. It is done here, once per group, rather than
        // inside the parallel worker: the mask depends only on the pattern, so
        // computing it per taxon was `n_taxa` redundant `n_samp`-long allocations
        // and the zero positions were then carried back as retained pairs.
        //
        // This is not cosmetic. `theta` is `colMeans(y - fitted, na.rm = TRUE)`,
        // so an NA there drops the cell from the mean while a 0 feeds `y` into
        // it -- the same distinction the NaN initialisation above exists for.
        let mut used = vec![false; n_samp];
        for &s in rows {
            used[s] = true;
        }
        // The fitted values come back parallel to the group's `rows`, so the
        // sample index is recovered from the cache rather than carried per value.
        for (t, row) in g.fitted {
            let base = t * n_samp;
            for (ri, &s) in rows.iter().enumerate() {
                fitted[base + s] = row[ri];
            }
            for (s, &u) in used.iter().enumerate() {
                if !u {
                    fitted[base + s] = 0.0;
                }
            }
        }
        for (t, d) in g.dof {
            dof[t] = d;
        }
    }

    FitAll { beta, fitted, dof }
}

/// Whether `.lm_fit_all`'s per-taxon `fit_one` would produce a fit for this
/// taxon, or abort and leave it at NA.
///
/// `fit_one` builds `data.frame(y_crt = Ymat[i, ], meta_data)` and calls
/// `lm(tformula, ...)` over *all* samples; `lm`'s default `na.action` drops the
/// rows where the response is NA, so the fit uses exactly the taxon's usable
/// samples -- the same rows as `xsub`. Two things make it fail there:
///
/// * **a factor with fewer than two observed levels.** A taxon seen in a handful
///   of samples can have all of them inside one group, and `lm` then aborts with
///   "contrasts can be applied only to factors with 2 or more levels". This is
///   the case that occurs on a heavily zero-inflated table, and it is why the
///   reference reports eight all-NA taxa on a 90%-zero cell.
/// * **no usable sample at all.** Already handled by the caller's empty-rows
///   branch, but stated here so the function is complete on its own.
///
/// A rank-deficient *numeric* column is deliberately **not** a failure here: `lm`
/// returns NA for an aliased coefficient and `fit_one` then writes zeros for the
/// dropped names, which is a finite result. Treating that as a failure would turn
/// every aliased taxon into an NA one, which is a different behaviour from the
/// reference's.
///
/// `cache.group_cols` empty means the caller does not know which columns are
/// group contrasts; the level test is then skipped rather than guessed, so the
/// pre-existing behaviour applies.
fn per_taxon_lm_would_succeed(cache: &DesignCache, xsub: &Matrix, rows: &[usize]) -> bool {
    if rows.is_empty() {
        // `lm` on a frame with no complete cases aborts, and `fit_one` leaves the
        // taxon at NA.
        return false;
    }
    if cache.group_cols.is_empty() {
        // The caller does not know which columns are group contrasts. Skipping the
        // test is the conservative choice for parity: it keeps the pre-existing
        // behaviour rather than guessing a level structure.
        return true;
    }
    // How many levels of the group factor this taxon actually saw.
    //
    // Level 1 is present when some usable sample lies outside every contrast
    // column -- under treatment contrasts the level-1 column *is* the intercept
    // column -- and each contrast column that is non-zero somewhere contributes
    // one more level.
    let has_level_1 = rows
        .iter()
        .enumerate()
        .any(|(ri, _)| cache.group_cols.iter().all(|&c| xsub.get(ri, c) == 0.0));
    let contrasted = cache
        .group_cols
        .iter()
        .filter(|&&c| {
            rows.iter()
                .enumerate()
                .any(|(ri, _)| xsub.get(ri, c) != 0.0)
        })
        .count();
    // The test is "at least two levels", not "every level present", and that
    // distinction is the whole rule. `lm` **re-levels**: a taxon observed only in
    // groups 3, 4 and 5 of a five-level factor is fitted against *those* levels,
    // the absent contrasts simply do not appear, and `fit_one` writes a literal 0
    // for `group2` and `group3` because `names(coef(fit))` omits them. Requiring
    // every level reported NA for 56 of 64 taxa in the `int-sparsity90-5group`
    // matrix cell, where the reference has finite coefficients.
    //
    // Only a taxon that saw a *single* level makes `lm` abort, with "contrasts can
    // be applied only to factors with 2 or more levels". That is the case on the
    // `sparsity-090` cell, where the reference reports eight all-NA taxa.
    usize::from(has_level_1) + contrasted >= 2
}

// `rows.len() < p` is deliberately *not* a failure condition. `lm` with fewer
// observations than parameters still returns a fit: it reports the rank-deficient
// solution over the columns it can determine, and `fit_one` writes those names
// and a literal 0 for the aliased ones.

/// Output of `.iter_mle`.
#[derive(Debug, Clone)]
pub struct IterMle {
    /// Wall time of the sandwich accumulation, reported separately because it is
    /// the one stage whose cost is linear in samples rather than in cells.
    pub sandwich_seconds: f64,
    pub beta: Vec<f64>,
    pub theta: Vec<f64>,
    /// `n_taxa * n_samp` row-major residuals `eps`.
    pub eps: Vec<f64>,
    /// `n_taxa * p * p` row-major sandwich covariance blocks.
    pub vcov: Vec<f64>,
    /// `n_taxa * p` row-major marginal variances.
    pub var_hat: Vec<f64>,
    /// Residual degrees of freedom, **one per taxon** (not per coefficient);
    /// the reference expands it with `matrix(rep(dof, n_fix_eff), ncol = ...)`.
    /// `None` when `theta` was estimated, matching the reference's `dof = NULL`.
    pub dof: Option<Vec<f64>>,
    pub iterations: usize,
    pub epsilon: f64,
    /// The missingness pattern each taxon was assigned to, 1-based.
    ///
    /// The contract's "missingness pattern assignment", and a Level A quantity:
    /// the assignment decides which taxa share one QR, so a different assignment
    /// is a different factorisation and a different `beta` for every taxon in the
    /// affected patterns. Numbering is first-appearance order, matching the
    /// reference's `factor(keys, levels = unique(keys))`.
    pub pattern_group: Vec<usize>,
    /// `epsilon` at the end of each iteration, in order: the contract's
    /// "convergence trace".
    ///
    /// The reference prints one line per iteration when `verbose = TRUE` and
    /// returns only the final `epsilon`, so this is the Rust half of a quantity
    /// that has to be captured from both sides -- the oracle from its printed
    /// trace, this from the loop.
    ///
    /// `iterations` and `epsilon` are kept as well because they are the two fields
    /// the reference does return, and a caller that only wants the summary should
    /// not have to know about the vector.
    pub trace: Vec<f64>,
}

/// `.iter_mle` with `theta = NULL`: estimate the sampling fractions iteratively.
pub fn iter_mle_estimate_theta(
    x: &Matrix,
    y: &[f64],
    n_taxa: usize,
    n_samp: usize,
    control: IterControl,
    compat: crate::config::CompatMode,
    group_cols: &[usize],
) -> Result<(IterMle, DesignCache)> {
    let observed = observed_mask(y, n_taxa, n_samp, x);
    trace_mem("mle:observed_mask");
    let cache = DesignCache::build(&observed, n_taxa, n_samp, x, group_cols);
    trace_mem("mle:design_cache");
    let p = x.cols;

    let mut theta = vec![0.0; n_samp];
    let mut beta = vec![f64::NAN; n_taxa * p];
    let mut fitted = vec![0.0; n_taxa * n_samp];
    let mut iterations = 0usize;
    let mut epsilon = 100.0f64;
    // The convergence trace, one entry per completed iteration.
    let mut trace: Vec<f64> = Vec::new();
    let mut adjusted = vec![0.0; n_taxa * n_samp];
    adjusted_response_into(y, n_taxa, n_samp, &theta, &mut adjusted);
    trace_mem("mle:adjusted");

    while epsilon > control.tol && iterations < control.max_iter {
        let fit = lm_fit_all(&cache, &adjusted, n_taxa, n_samp);
        trace_mem("mle:lm_fit_all");
        let new_theta = theta_new(y, &fit.fitted, n_taxa, n_samp);
        epsilon = delta_norm(&fit.beta, &beta, &new_theta, &theta);
        beta = fit.beta;
        fitted = fit.fitted;
        theta = new_theta;
        iterations += 1;
        trace.push(epsilon);
        adjusted_response_into(y, n_taxa, n_samp, &theta, &mut adjusted);
        trace_mem("mle:iter end");
    }

    let eps = residuals(y, &fitted, n_taxa, n_samp, &theta);
    trace_mem("mle:residuals");
    let t_sandwich = std::time::Instant::now();
    let xx = OuterProducts::build(x);
    let xtx_inv = crate::vcov::xtx_inverse(x);
    // The sandwich squares the residuals on the fly, so no second
    // `n_taxa x n_samp` buffer is allocated. See `sandwich_all`.
    trace_mem("mle:outer_products");
    let (vcov, var_hat) = sandwich_all(&eps, n_taxa, n_samp, &xx, &xtx_inv, compat);
    trace_mem("mle:sandwich_all");

    Ok((
        IterMle {
            beta,
            theta,
            eps,
            vcov,
            var_hat,
            dof: None,
            iterations,
            epsilon,
            trace,
            pattern_group: group_of_taxon(&cache),
            sandwich_seconds: t_sandwich.elapsed().as_secs_f64(),
        },
        cache,
    ))
}

/// Report the resident set at a point inside the MLE, if tracing is enabled.
///
/// Duplicated from `pipeline::Stage::trace_mem` rather than shared, because this
/// is the module where the peak actually is, and the trace is a debugging aid
/// that should not force a dependency from the numerics onto the pipeline.
fn trace_mem(name: &str) {
    static ON: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    if !*ON.get_or_init(|| std::env::var_os("ANCOMBC2_TRACE_MEM").is_some()) {
        return;
    }
    let kb = std::fs::read_to_string("/proc/self/status")
        .ok()
        .and_then(|s| {
            s.lines()
                .find(|l| l.starts_with("VmRSS:"))
                .and_then(|l| l.split_whitespace().nth(1))
                .and_then(|v| v.parse::<u64>().ok())
        });
    if let Some(kb) = kb {
        eprintln!("[mem]   {name:<24} rss={:.3} GB", kb as f64 / 1e6);
    }
}

/// `.iter_mle` with `theta` supplied: a single fit, no iteration.
/// The second MLE, with `theta` held fixed.
///
/// The missingness grouping is re-derived from *this* response rather than
/// reused from the first MLE. The reference calls `.iter_mle` twice, once per
/// response, and `.iter_mle` recomputes the groups each time; the two responses
/// can have different taxon sets (the second is `O2`, the aggregate-level table
/// with the structural zeros removed) and even the same set can have a different
/// missingness pattern, so a cached grouping would silently index the wrong
/// rows. Keeping the parameter would invite exactly that mistake.
pub fn iter_mle_fixed_theta(
    x: &Matrix,
    y: &[f64],
    n_taxa: usize,
    n_samp: usize,
    theta: &[f64],
    compat: crate::config::CompatMode,
    group_cols: &[usize],
) -> Result<IterMle> {
    let adjusted = adjusted_response(y, n_taxa, n_samp, theta);
    let observed = observed_mask(y, n_taxa, n_samp, x);
    trace_mem("mle:observed_mask");
    let cache = DesignCache::build(&observed, n_taxa, n_samp, x, group_cols);
    trace_mem("mle:design_cache");
    let fit = lm_fit_all(&cache, &adjusted, n_taxa, n_samp);
    let eps = residuals(y, &fit.fitted, n_taxa, n_samp, theta);
    let t_sandwich = std::time::Instant::now();
    let xx = OuterProducts::build(x);
    let xtx_inv = crate::vcov::xtx_inverse(x);
    // The sandwich squares the residuals on the fly, so no second
    // `n_taxa x n_samp` buffer is allocated. See `sandwich_all`.
    trace_mem("mle:outer_products");
    let (vcov, var_hat) = sandwich_all(&eps, n_taxa, n_samp, &xx, &xtx_inv, compat);
    trace_mem("mle:sandwich_all");
    Ok(IterMle {
        // `.iter_mle` with `theta` supplied does not iterate, so the reference
        // prints no convergence trace for this call and neither does this. One
        // entry would imply an iteration that did not happen.
        trace: Vec::new(),
        pattern_group: group_of_taxon(&cache),
        sandwich_seconds: t_sandwich.elapsed().as_secs_f64(),
        beta: fit.beta,
        theta: theta.to_vec(),
        eps,
        vcov,
        var_hat,
        dof: Some(fit.dof),
        iterations: 0,
        epsilon: 0.0,
    })
}

/// The 1-based pattern id of every taxon, in taxon order.
///
/// `cache.taxa[g]` holds the taxa of pattern `g` in ascending order, so this is
/// an inverse permutation of the grouping rather than a second grouping -- there
/// is no second place for the assignment to be computed, and therefore no second
/// place for it to be wrong.
fn group_of_taxon(cache: &DesignCache) -> Vec<usize> {
    let n_taxa: usize = cache.taxa.iter().map(|g| g.len()).sum();
    let mut out = vec![0usize; n_taxa];
    for (g, taxa) in cache.taxa.iter().enumerate() {
        for &t in taxa {
            out[t] = g + 1;
        }
    }
    out
}

fn observed_mask(y: &[f64], n_taxa: usize, n_samp: usize, x: &Matrix) -> Vec<bool> {
    let x_ok: Vec<bool> = (0..x.rows)
        .map(|j| (0..x.cols).all(|a| x.get(j, a).is_finite()))
        .collect();
    let mut out = vec![false; n_taxa * n_samp];
    for i in 0..n_taxa {
        for j in 0..n_samp {
            out[i * n_samp + j] = y[i * n_samp + j].is_finite() && x_ok[j];
        }
    }
    out
}

fn adjusted_response(y: &[f64], n_taxa: usize, n_samp: usize, theta: &[f64]) -> Vec<f64> {
    let mut out = vec![0.0; n_taxa * n_samp];
    adjusted_response_into(y, n_taxa, n_samp, theta, &mut out);
    out
}

/// `y - theta` per sample, written into `out`.
///
/// The iteration reuses one buffer rather than allocating a fresh
/// `n_taxa x n_samp` response on every one of its (up to twenty) passes. A
/// `Vec` of that size is a large allocation and a large memset, repeated; a
/// borrowed buffer is neither. `NaN` entries stay `NaN`, as `sweep` leaves them.
fn adjusted_response_into(y: &[f64], n_taxa: usize, n_samp: usize, theta: &[f64], out: &mut [f64]) {
    debug_assert_eq!(out.len(), n_taxa * n_samp);
    for i in 0..n_taxa {
        let row = &y[i * n_samp..(i + 1) * n_samp];
        let dst = &mut out[i * n_samp..(i + 1) * n_samp];
        for j in 0..n_samp {
            dst[j] = if row[j].is_finite() {
                row[j] - theta[j]
            } else {
                f64::NAN
            };
        }
    }
}

/// R: `theta_new = colMeans(y - y_crt_hat, na.rm = TRUE)`.
fn theta_new(y: &[f64], fitted: &[f64], n_taxa: usize, n_samp: usize) -> Vec<f64> {
    // Level 4 of the nesting order: the taxa axis, reached as independent
    // per-sample means.
    //
    // This is a `colMeans` over taxa, so the *outer* loop is over samples and each
    // sample's mean is over a fixed ascending run of taxa. Splitting the samples
    // therefore gives each worker a whole, self-contained reduction whose summation
    // order is the one the serial loop used -- the result is bit-identical, and no
    // deterministic-order accumulator is needed because there is no cross-sample
    // reduction to reorder.
    //
    // Note what is *not* parallelised here: the inner taxon loop. Splitting it
    // would need a tree or pairwise sum to keep the order, which changes the last
    // bit for no gain, and it would nest inside the sample split besides.
    let mut lvl = crate::parallel::NestingBudget::level("taxa (per-sample means)");
    let samples: Vec<usize> = (0..n_samp).collect();
    crate::parallel::map_par(&mut lvl, &samples, |_k, &j| {
        // `colMeans(y - fitted, na.rm = TRUE)`: the difference is formed in double,
        // then summed and divided in R's `long double`. `na.rm` skips NaN only.
        let d = (0..n_taxa).map(|i| y[i * n_samp + j] - fitted[i * n_samp + j]);
        let mean = crate::reduce::mean_na_rm(d);
        // An all-NaN column is `NaN` here, as before, not `NA_real_`: the iteration
        // tests it with `is_nan`, and the payload is not part of any stage output.
        if mean.is_nan() {
            f64::NAN
        } else {
            mean
        }
    })
}

/// R: `sqrt(sum((beta_new - beta)^2, na.rm) + sum((theta_new - theta)^2, na.rm))`.
fn delta_norm(beta_new: &[f64], beta: &[f64], theta_new: &[f64], theta: &[f64]) -> f64 {
    let mut s = 0.0;
    for (a, b) in beta_new.iter().zip(beta) {
        if a.is_finite() && b.is_finite() {
            let d = a - b;
            s += d * d;
        }
    }
    for (a, b) in theta_new.iter().zip(theta) {
        if a.is_finite() && b.is_finite() {
            let d = a - b;
            s += d * d;
        }
    }
    s.sqrt()
}

/// R: `eps = t(t(y - y_crt_hat) - theta)`.
fn residuals(y: &[f64], fitted: &[f64], n_taxa: usize, n_samp: usize, theta: &[f64]) -> Vec<f64> {
    let mut out = vec![0.0; n_taxa * n_samp];
    for i in 0..n_taxa {
        for j in 0..n_samp {
            let a = y[i * n_samp + j];
            let b = fitted[i * n_samp + j];
            let t = theta[j];
            out[i * n_samp + j] = if a.is_finite() && b.is_finite() && t.is_finite() {
                (a - b) - t
            } else {
                f64::NAN
            };
        }
    }
    out
}

/// The reference's over-parameterisation guard.
///
/// `.iter_mle` fits `lm(rnorm(n) ~ fix_formula)` as a smoke test and stops if any
/// coefficient is `NA` (an unidentifiable covariate) or the residual degrees of
/// freedom are zero. Reproduced here so the error surfaces at the same point,
/// with the same message.
///
/// # Incomplete design rows
///
/// The smoke test is an `lm` fit, and `lm` drops the rows whose design entry is
/// `NA` -- a factor level that does not exist, a covariate with a missing value.
/// The Gram matrix is therefore built from the *complete* rows only. Including
/// the incomplete ones makes `X'X` NaN, every eigenvalue NaN, and the guard then
/// reports every covariate as unidentifiable: a warning that looks like a
/// collinearity problem and is really a missing value. `atlas1006` has 37 of
/// 1151 samples with an unknown `sex`, and that is how they surfaced.
pub fn check_design_identifiable(x: &Matrix) -> Result<()> {
    let n = x.rows;
    let p = x.cols;
    if p == 0 {
        return Err(AncombcError::UnidentifiableCovariates {
            covariates: "(none)".into(),
        });
    }
    // `lm`'s complete-case rows, and the sample count they leave behind, since
    // that is also the `n` in the residual-degrees-of-freedom test below.
    let complete: Vec<usize> = (0..n)
        .filter(|j| (0..p).all(|a| x.get(*j, a).is_finite()))
        .collect();
    let n_used = complete.len();
    // The rank test is over the complete rows. `xtx` is accumulated by hand
    // rather than by `x.t().matmul(x)`, because a matrix product cannot skip
    // rows and the filtered sub-matrix would have to be materialised first.
    let mut xtx = vec![0.0f64; p * p];
    for &j in &complete {
        for a in 0..p {
            let xa = x.get(j, a);
            for b in 0..p {
                xtx[a * p + b] += xa * x.get(j, b);
            }
        }
    }
    let (vals, _) = crate::matrix::eigen_symmetric(&Matrix::from_vec(p, p, xtx)?);
    let dmax = vals.iter().fold(0.0f64, |m, v| m.max(v.abs()));
    let tol = (n_used.max(p) as f64) * f64::EPSILON * dmax.max(f64::MIN_POSITIVE);
    // `!(a > b)` rather than `a <= b`, stated explicitly because it is
    // deliberate: an eigenvalue may be NaN, and every comparison involving NaN
    // is false, so `!(|v| > tol)` reports the column while `|v| <= tol` would
    // silently pass it.

    #[allow(
        clippy::neg_cmp_op_on_partial_ord,
        reason = "a NaN eigenvalue must be reported, not pass the rank test"
    )]
    let bad: Vec<String> = vals
        .iter()
        .take(p)
        .enumerate()
        .filter(|(_, v)| !(v.abs() > tol))
        .map(|(k, _)| x.col_name(k))
        .collect();
    if !bad.is_empty() {
        return Err(AncombcError::UnidentifiableCovariates {
            covariates: bad.join(", "),
        });
    }
    if n_used == p {
        // R's guard is `lm_smoke$df.residual == 0`, i.e. exactly zero rather than
        // non-positive. With n == p that is the case; with n < p the rank check
        // above has already fired, because n < p makes X'X singular. The count
        // is of *complete* rows, because that is the `n` `lm` fits on.
        return Err(AncombcError::NoResidualDegreesOfFreedom);
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::{lm_fit_all, DesignCache};

    /// `~ group + x1` on 30 samples with five levels: an intercept, four
    /// treatment contrasts, and one continuous column. Column order matches
    /// `model.matrix`.
    fn design() -> crate::matrix::Matrix {
        let n_samp = 30usize;
        let p = 6usize;
        let mut x = crate::matrix::Matrix::zeros(n_samp, p);
        let names = ["(Intercept)", "group2", "group3", "group4", "group5", "x1"]
            .iter()
            .map(|s| s.to_string())
            .collect::<Vec<_>>();
        for j in 0..n_samp {
            x.set(j, 0, 1.0);
            // five levels of six samples each; the reference level's contrast is
            // absent, absorbed into the intercept.
            let level = j / 6;
            for a in 1..=4usize {
                x.set(j, a, if level == a { 1.0 } else { 0.0 });
            }
            // Non-zero and distinct per sample, so it is never confusable with a
            // level indicator.
            x.set(j, 5, 0.5 + j as f64 * 0.37);
        }
        x.colnames = names;
        x
    }

    /// One taxon observed on exactly `keep`, everything else missing.
    fn observed_on(keep: &[usize], n_samp: usize) -> Vec<bool> {
        let mut v = vec![false; n_samp];
        for &j in keep {
            v[j] = true;
        }
        v
    }

    /// The base level is the first **observed** one, not the first globally, so a
    /// taxon that never saw level 1 still reports a coefficient for it.
    ///
    /// `model.matrix` re-levels the factor to what the na.omit-reduced frame
    /// contains and then applies `contr.treatment`, which drops the first level
    /// *it sees*. The test case is samples 12, 18 and 24 -- levels 3, 4 and 5 of
    /// the five, with levels 1 and 2 unobserved -- so the base is level 3, `group3`
    /// is the one whose name never appears in `coef()`, and `group4`/`group5` both
    /// do.
    ///
    /// This is the case that distinguishes the two candidate rules. Keeping the
    /// smallest *contrast index* instead of the smallest *level* drops `group3`
    /// here, which is 20 of the 64 rows on `int-sparsity90-5group`.
    #[test]
    fn the_base_level_is_the_first_observed_one() {
        let x = design();
        let n_samp = x.rows;
        let p = x.cols;
        let obs = observed_on(&[12, 18, 24], n_samp);
        let cache = DesignCache::build(&obs, 1, n_samp, &x, &[1usize, 2, 3, 4]);

        let mut y = vec![f64::NAN; n_samp];
        for (k, &j) in [12usize, 18, 24].iter().enumerate() {
            y[j] = 1.0 + k as f64;
        }
        let fit = lm_fit_all(&cache, &y, 1, n_samp);
        let got = &fit.beta[0..p];

        let zeros: Vec<usize> = (0..p).filter(|&a| got[a] == 0.0).collect();
        let nan: Vec<usize> = (0..p).filter(|&a| got[a].is_nan()).collect();
        let solved: Vec<usize> = (0..p)
            .filter(|&a| got[a].is_finite() && got[a] != 0.0)
            .collect();

        // Two zeros, for two different reasons, and keeping them apart is the
        // point of the test. `group2` (column 1) is a level this taxon never
        // observed, so it has no column in `lm`'s model matrix at all.
        // `group3` (column 2) *is* observed -- it is 1 on sample 12 -- and is
        // still zeroed, because it is the first level the reduced factor has and
        // `contr.treatment` drops that one. Under the rule "the lowest contrast
        // index observed is the base" it would be column 1 that is dropped and
        // column 2 that is fitted, which is the bug this pins.
        assert_eq!(
            zeros,
            vec![1, 2],
            "group2 is unobserved and group3 is the re-based base; both are 0, for \
             different reasons; got {got:?}"
        );
        assert_eq!(
            solved,
            vec![0, 3, 4],
            "levels 4 and 5 are observed and are not the base, so both are fitted; \
             got {got:?}"
        );
        assert_eq!(
            nan,
            vec![5],
            "three samples and four columns, so `x1` is the aliased one; got {got:?}"
        );
    }

    /// `stats::fitted` is `lm.fit`'s projection, so it stays finite when `coef()`
    /// is `NA`, and it interpolates when `n == rank`.
    ///
    /// The distinction is invisible for a full-rank fit and decisive for this one:
    /// `theta` is `colMeans(y - fitted, na.rm = TRUE)`, so an `NA` in `fitted` does
    /// not cost one term, it drops the taxon out of that sample's mean. On
    /// `int-sparsity90-5group` that left about a fifth of each sample's observed
    /// taxa contributing and three samples with no contributor at all, so `theta`
    /// came back `NaN` and cascaded on the next iteration.
    #[test]
    fn fitted_is_the_projection_and_interpolates_when_n_equals_rank() {
        let x = design();
        let n_samp = x.rows;
        let obs = observed_on(&[12, 18, 24], n_samp);
        let cache = DesignCache::build(&obs, 1, n_samp, &x, &[1usize, 2, 3, 4]);

        let mut y = vec![f64::NAN; n_samp];
        for (k, &j) in [12usize, 18, 24].iter().enumerate() {
            y[j] = 1.0 + k as f64;
        }
        let fit = lm_fit_all(&cache, &y, 1, n_samp);

        assert!(
            fit.beta[5].is_nan(),
            "the premise: `x1` is aliased, so `coef` carries NA"
        );
        for &j in [12usize, 18, 24].iter() {
            let f = fit.fitted[j];
            assert!(
                f.is_finite(),
                "sample {j}: `fitted` must be the projection, not `X %*% coef`, so \
                 an aliased coefficient cannot make it NA; got {f}"
            );
            assert!(
                (f - y[j]).abs() < 1e-9,
                "sample {j}: three observations against a rank-3 model interpolates \
                 exactly, so the fitted value is the response; got {f} against {}",
                y[j]
            );
        }
    }

    /// A taxon observed on three samples of a five-level factor, spanning three
    /// levels, must come back as **three finite coefficients, two literal zeros
    /// and one `NA`** -- not three zeros.
    ///
    /// Measured on the oracle's own `.lm_fit_all` with this exact design (30
    /// samples, samples 1/7/13 observed, groups 1-3 spanned):
    ///
    /// ```text
    /// (Intercept)  group2  group3  group4  group5      x1
    ///  1.3586796  0.07434 -1.2843  0.00000  0.00000      NA
    /// ```
    ///
    /// The three outcomes are distinct and each is reproduced by a different rule
    /// in the reference:
    ///
    /// * the three spanned levels get coefficients, `lm` having re-levelled to
    ///   them;
    /// * `group4` and `group5` are **absent from `coef()`** because their levels
    ///   were not observed, and `fit_one` starts from `bi = rep(0, p)`, so they are
    ///   literal zeros;
    /// * `x1` is **present but aliased**, and `coef.lm` keeps an aliased name with
    ///   the value `NA`.
    ///
    /// What distinguishes the last two is whether the *name survives into
    /// `coef()`*, not whether the coefficient is estimable. Both are aliased here;
    /// only one is reported as `NA`.
    /// **Ignored, and that is the finding.** The expectation below is the
    /// oracle's, measured through its own `.lm_fit_all`. This implementation does
    /// not meet it, and the reason is narrower and more specific than "aliased
    /// coefficients are handled differently": on a group with fewer observed
    /// samples than design columns, both sides pick a rank-`n` basis for the fit,
    /// and **the bases differ**. R's `dqrls` chose
    /// `(Intercept, group2, group3)`; this implementation's rank-revealing QR
    /// chose `(Intercept, group2, x1)`. Both interpolate the same three points
    /// exactly, so the coefficients are different representatives of the same
    /// affine family -- which is why the divergence is a 195% difference in
    /// `beta_star` rather than a rounding difference.
    ///
    /// That also decides *which* column is reported aliased: `x1` under R, `group3`
    /// here. Writing `NA` for a varying aliased column -- the right rule, taken
    /// from the measurement -- was implemented and produced `NA` on `group3` instead,
    /// so it moved the divergence rather than closing it.
    ///
    /// Closing this needs the pivot order of `stats::lm`'s `lm.fit` (`dqrls` /
    /// `dqrdc2`, which maximises `|R(j,j)|` at each step) reproduced, not just its
    /// rank and tolerance. Until then the divergence is reported and scoped in
    /// `docs/reference_behavior.md` section 16 rather than papered over, and this
    /// test stands as the executable statement of what is being matched.
    #[test]
    fn an_aliased_column_is_na_and_an_absent_name_is_zero() {
        let x = design();
        let n_samp = x.rows;
        let p = x.cols;
        let group_cols = vec![1usize, 2, 3, 4];
        // samples 0, 6, 12 -> levels 0, 1, 2 of the five
        let obs = observed_on(&[0, 6, 12], n_samp);
        let cache = DesignCache::build(&obs, 1, n_samp, &x, &group_cols);

        // The response is `n_taxa x n_samp`; the unobserved samples are `NaN`,
        // which is how the reference's own table reaches `.lm_fit_all`.
        let mut y = vec![f64::NAN; n_samp];
        for (k, &j) in [0usize, 6, 12].iter().enumerate() {
            y[j] = 1.0 + k as f64;
        }
        let fit = lm_fit_all(&cache, &y, 1, n_samp);
        let got = &fit.beta[0..p];

        // Three outcomes, so three disjoint classes -- and `0.0` is one of them,
        // which is why this cannot be a count of "finite" values: the two
        // re-levelled zeros are finite. `|coefficient|` separates them, and the
        // spanned levels' coefficients are not zero on this response.
        let solved: Vec<usize> = (0..p)
            .filter(|&a| got[a].is_finite() && got[a] != 0.0)
            .collect();
        let zeros: Vec<usize> = (0..p).filter(|&a| got[a] == 0.0).collect();
        let nan: Vec<usize> = (0..p).filter(|&a| got[a].is_nan()).collect();

        assert_eq!(
            nan,
            vec![5],
            "x1 varies and sits beyond the rank, so `coef.lm` names it and \
             reports NA; got {got:?}"
        );
        assert_eq!(
            zeros,
            vec![3, 4],
            "group4 and group5 are constant over the fitted samples, so `lm` \
             re-levels them away, their names never appear, and `fit_one`'s \
             `rep(0, p)` leaves them at zero; got {got:?}"
        );
        assert_eq!(
            solved,
            vec![0, 1, 2],
            "the first `rank` columns are the ones kept, in the order they were \
             built -- that is `dqrls`' rule, with no column selection; got {got:?}"
        );
        assert_eq!(
            solved.len() + nan.len() + zeros.len(),
            p,
            "every column must be in exactly one class; got {got:?}"
        );
    }
}
