# The algorithm

What ANCOM-BC2 computes, in the order the code computes it. The reference is
ANCOMBC 2.15.2 at `dc4febdf5`; every "the reference does X" below was checked
against that commit, and `docs/reference_behavior.md` records where the check
found something other than what the source appears to say.

## The problem

Count data from a microbiome experiment is compositional and sparse. Two things
corrupt a naive differential-abundance test:

* **sampling fractions.** A taxon that happens to be sequenced shallowly in some
  samples has its log abundance shifted, and that shift looks exactly like a
  treatment effect if the treatment is unbalanced.
* **sampling variation.** Even at a fixed fraction, the counts are overdispersed
  relative to the Poisson variance a regression assumes.

ANCOM-BC2 addresses the first with an explicit bias term and the second with a
sandwich estimator that does not assume the variance is Poisson.

## Step 0: two tables

The reference filters twice and keeps both results, and this is not an
implementation detail — see `docs/reference_behavior.md` §1.

* `O1` = prevalence and library-size filtered. Used for `y1`, the first MLE, the
  E-M bias, and the sampling fractions.
* `O2` = `O1` minus the taxa with structural zeros. Used for `y2`, the second MLE,
  and the entire reported table.

A taxon with a structural zero therefore still contributes to the estimate of
`theta`. Dropping it earlier changes `theta`, and `theta` changes every reported
coefficient.

## Step 1: filtering and structural zeros

* **Prevalence.** A taxon is dropped if its prevalence — non-zero samples over
  observed samples — is below `prv_cut` (default 0.1).
* **Library size.** A sample is dropped if its total count is below `lib_cut`
  (default 0). The cutoff is computed on the taxon-filtered table.
* **Structural zeros.** `.get_struc_zero` marks a taxon as a structural zero in a
  group when it is (nearly) absent there by construction, not by sampling. A taxon
  with a structural zero anywhere is removed from the primary analysis. With
  `neg_lb = TRUE` the classification uses the asymptotic lower bound instead of a
  fixed threshold.

## Step 2: the first MLE, and `theta`

`y1 = log(O1 + pseudo)`, with `log(0)` mapped to `NA`, then centred per taxon:
`y1[i, ] - mean(y1[i, ], na.rm = TRUE)`. The centring is what makes the
coefficients contrasts rather than abundances.

The MLE alternates:

```
theta = 0
repeat up to max_iter times, while epsilon > tol:
    fit X beta ~ y1 - theta          (one lm per missingness pattern, see below)
    theta_new = colMeans(y1 - fitted, na.rm = TRUE)
    epsilon = sqrt(sum((beta_new - beta)^2) + sum((theta_new - theta)^2))
    beta, theta = beta_new, theta_new
```

`tol = 0.01`, `max_iter = 20`. `epsilon` is a sum of squares over every taxon and
every sample, so for a realistic table it never falls below `tol` and the loop
runs to `max_iter`. That is the reference's behaviour and it is load-bearing: the
reported `beta` is the twentieth fit, and the reported `theta` is the twentieth
update, which was computed from the *nineteenth* fit's fitted values. The two are
deliberately not consistent with each other, because the reference does not make
them consistent.

The residual is `eps = t(t(y1 - fitted) - theta)`, the same mismatched pair.

## Step 3: the missingness-pattern grouping

`.lm_fit_all` solves every taxon at once, but only the *factorisation* can be
shared: two taxa can be solved together only if they are observed on the same
samples, because only then is the sub-design the same.

* A taxon's pattern is the set of samples where its response is finite and the
  design row is complete.
* Taxa are grouped by pattern; the patterns are keyed in the order they first
  appear, matching `split(seq_len(n_tax), factor(keys, levels = unique(keys)))`.
* One QR per pattern, cached and reused across all 20 iterations. The
  factorisation depends only on the design, not on the response, so this is the
  single largest win in the algorithm.
* The response matrix for a group is `n_used x n_group`, and `lm.fit` solves it
  with one call, so the cost is one triangular solve per coefficient, not per
  taxon.
* **`Q` is never formed.** `Qr` stores the unit Householder reflectors and
  applies them to the right-hand sides instead: with `Q = H_0 H_1 ... H_{k-1}`
  and each `H` symmetric, `Q' b` is the reflectors applied to `b` in *ascending*
  order, at `O(n p)` per response rather than `O(n^2)` to build `Q`. This matters
  because `n` here is the number of *rows* of a pattern's sub-design — the number
  of samples in which a group's taxa are all observed — which reaches into the
  thousands on a real table while `p` stays under ten. An earlier version
  materialised the `n x n` factor and the 500x5000 benchmark spent 47.3 s of its
  47.6 s first-MLE budget in the cache build alone; applying the reflectors
  directly took the same dataset from 86.5 s to 1.29 s end to end.
* The group response block is solved **in place** (`Qr::solve_multi_into`). It is
  built once, per group, per MLE iteration, and the out-of-place solve copied it
  in full to apply the reflectors; solving in place removes that copy and one
  more live buffer at the peak. On this host the compiler was already eliding the
  copy, so the measured effect is nil — it is kept because the bound is then
  explicit rather than left to an optimisation pass.
* Neither the design sub-matrices nor the count reader allocates per cell. A
  counts matrix is parsed through `split_line_cow`, which borrows each field on
  the unquoted path (one `String` per field is one heap allocation per cell to
  hold text that is immediately parsed to `f64`), and the pattern designs are
  built by `Matrix::select_rows_anon`, which omits row labels the grouped
  factorisation never reads (one `String` per row, per group). Together these
  were 4.9 allocations per input cell on the 1000x10000 benchmark, 48.6 million
  in total, and are now 2.6 million. See `docs/compatibility.md` for the
  measurement.
  `Qr::explicit_q` still builds `Q`, and exists only so the test that checks the
  reflector application against the direct product can run.

Two edge cases, both reproduced:

* A group with no usable sample refits each taxon alone with `lm`, which fails for
  a single-level factor and leaves the taxon `NA` with `dof = 999`.
* A rank-deficient sub-design also refits per taxon, and `lm` *drops* the aliased
  coefficient, which `.lm_fit_all` then writes as a literal `0`. See
  `docs/reference_behavior.md` §10.

## Step 4: the sandwich covariance

`V_i = (X'X)^- (sum_j eps_ij^2 x_j x_j') (X'X)^-`, the HC0 estimator, per taxon
and per coefficient block.

The per-sample outer products `x_j x_j'` do not depend on the taxon, so they are
computed once into an `n_samp x p^2` matrix and the accumulation runs over all
taxa together — a single pass that is cache-friendly in the taxon axis.

The reference fills a missing entry of the outer-product matrix with `0.1` rather
than skipping the term. `CompatMode` chooses.

`X'X` is computed once for the whole run, not per taxon, because the design is the
same for every taxon.

## Step 5: the E-M bias estimate

For each coefficient, the reference fits a three-component Gaussian mixture to
`beta` across taxa, with component variances `nu1, nu2, nu3` and weights
`pi = (0.75, 0.125, 0.125)`, initialised at the observed variance. The posterior
probabilities give `delta_wls` and `var_delta`; `delta_em` is the E-step output.

The two things that are easy to get wrong, and are reproduced because they are
load-bearing:

* The loop reads `l1`, `l2` and `delta` from the *previous* iteration when writing
  `l1_new`, `l2_new` and the Nelder-Mead objective. A "clean" version differs from
  the reference at the fifth significant figure.
* The variance components are optimised with NLopt's Nelder-Mead, whose effective
  initial simplex is empirically `0.75 * |x0|`, not the documented `0.05`.

## Step 6: the bias correction and the sampling fractions

`beta1 = t(t(beta1) - delta_em)`, and then

```
theta_hat[i, ] = y1[i, ] - rowSums(x * rep(beta1[i, ], each = n_samp), na.rm = TRUE)
theta_hat       = colMeans(theta_hat, na.rm = TRUE)
```

so `theta_hat` is the per-sample bias, estimated over *every* taxon of the bias
set. A sample whose estimate is `NA` — an excessive number of zero values — makes
the reference warn, and the warning is reproduced.

The `na.rm` in the inner `rowSums` is why an unfitted taxon's coefficient is
dropped rather than propagated; see `docs/numerical_contract.md` §2.

## Step 7: the second MLE

`y2 = log(O2 + pseudo)`, centred, minus `theta_hat`, and `lm` fitted with
`theta` held fixed. The missingness grouping is recomputed, because the response
is a different matrix with a different taxon set.

`dof` is `n_used - rank` for this stage, one value per taxon, expanded to one per
coefficient afterwards.

## Step 8: variance, `s0`, and inference

```
var_hat = var_hat + var_delta + 2 * sqrt(var_hat * var_delta)   # the delta variance
s02     = quantile(var_hat, s0_perc, na.rm = TRUE)              # per coefficient
var_hat = var_hat + s02
se      = sqrt(var_hat)
W       = beta / se
p       = 2 * pt(abs(W), df = dof, lower.tail = FALSE)
q       = p.adjust(p, method = p_adj_method)                   # per coefficient
diff_abn = q <= alpha
```

`s0` is the SAM regulariser: a small positive constant added to every variance so
that a coefficient estimated from a handful of observations cannot dominate the
mixture fit downstream. It is the 5th percentile of the variance distribution, per
coefficient, not a global constant.

The primary p-value uses **Student-t** tails with the stage-2 `dof`, not normal
tails. In the first stage `dof` is `NULL` and R's `pt(x, df = NULL)` evaluates the
default `df = Inf` — the normal limit. Both branches are reproduced.

## Step 9: the multi-group tests

* **Global.** A quadratic Wald test of `H0: beta_group = 0` over the group
  coefficients, with the off-diagonal terms of the sandwich block. The null
  distribution is the chi-square on the number of group coefficients.
* **Pairwise.** The published mixed-directional FDR: screen with the global test
  at `BH`, then adjust the pairwise p-values *within each taxon* at
  `n = n_col * n_tax / R` with `mdfdr_control$fwer_ctrl_method`, where `R` is the
  number of taxa the screen rejected. A taxon that fails the screen gets `p = 1`.

  The reference divides by `R` without guarding against `R == 0`, so every taxon is
  called significant when the screen rejects nothing. That is reproduced
  deliberately and flagged in the source, because it is a property of the oracle
  and not a defensible result.

## Step 10: the sensitivity analysis

`pseudo_sens` asks whether a taxon's significance call survives a different
pseudo-count.

* **Conservative** (`conservative = TRUE`, the default): the whole algorithm is
  rerun at pseudo-counts 0.1, 0.5 and 1 in addition to the main run's, and
  `ss_tab` is the proportion of runs whose adjusted p-value exceeds `alpha`.
* **Non-conservative**: `theta_hat` is estimated **once**, on the main run, and
  only the inference is refitted — 50 times, at pseudo-counts 0.01 to 0.50 in
  steps of 0.01.

Either way `passed_ss` is `ss_tab == 0 & p_main <= alpha` or
`ss_tab == 1 & p_main > alpha`, and `diff_robust = diff_abn & passed_ss`.

The non-conservative refit is a per-taxon `lm` on the theta-adjusted data, so it
is embarrassingly parallel and dominates the runtime of a large analysis. That is
where the thread count matters, and it is the outer loop of the nesting budget.

## The rank-deficient path

One missingness group per distinct observation pattern, one cached factorisation
each. When a group's factor is rank deficient the group is **not** pseudo-solved:
every taxon of it is refit on its own, which is what `.lm_fit_all`'s `fit_one`
does. Three rules there are easy to get wrong, and all three are reproduced:

* **The per-taxon fit is over the taxon's own usable samples**, not the group's.
* **It can abort.** `lm` fails on a factor with fewer than two observed levels, and
  the taxon is then left at `NA` for both its coefficients and its fitted values.
  `theta` is a `na.rm` mean of `y - fitted`, so an `NA` and a `0` are different
  answers, not the same answer written differently.
* **It re-levels.** A taxon that saw only groups 3, 4 and 5 of a five-level factor
  is fitted against those, and the two absent contrasts are reported as literal
  zeros rather than as `NA`.

A group can also be narrower than the design is tall, which is legal and has to
factor: `rank <= min(n, p)`, and the caller refits per taxon. See
`docs/reference_behavior.md` §13 for the reference's source and the measurements.

## Threading

**Implemented.** One global Rayon pool and an explicit nesting budget; `parallel.rs`
is the whole mechanism.

### The budget

The pool is claimed by exactly one level at a time, and the nesting order is the
plan's:

```text
pseudo-count runs        sens.rs, pipeline.rs     outermost
  E-M coefficients       em.rs
    missingness patterns  DesignCache::build      <- the axis that matters
    missingness groups    lm_fit_all
      taxa                 solve_multi_into, theta_new, the sandwich's blocks
```

A level parallelises only if nothing above it has claimed the pool, and a level
with no work to share **releases its claim** so the level below can take it. That
release is not a nicety: a table with no missing values has exactly one missingness
pattern, so the group level has a single item, and if it held the pool against the
work inside the group the whole analysis would run on one core with `--threads 16`.

`map_par` takes the *same* closure on both paths. A hand-written serial fallback
would be a second implementation that could drift from the parallel one, which is a
bug waiting for the inputs that take one branch and not the other.

### Why no axis needs a deterministic-order reduction

Each of the parallel levels writes to **disjoint** slots:

| level | writes | why it is disjoint |
| --- | --- | --- |
| missingness patterns | `qrs[g]`, `designs[g]` | one entry per pattern |
| missingness groups | `beta[t*p..]`, `fitted[t*n_samp..]`, `dof[t]` | a taxon is in exactly one group |
| taxa, in the solve | one column of the response block | a column is one taxon |
| taxa, in `theta_new` | `theta[j]` | each sample's mean is self-contained |
| taxa, in the sandwich | `vcov[i*p*p..]`, `var_hat[i*p..]` | one taxon block each |

So none of them has an order that is observable in the result, and splitting them
cannot change a bit. The reductions whose order *is* observable -- the E-M's
parameter sweep, and the sandwich's accumulation over samples within a taxon --
stay serial. The P15 property test checks 1 thread against 8 for bitwise equality
of the coefficients, q, se, samp_frac, delta_em, diff_abn and both multi-group
tests, and the CLI's `res.tsv` is byte-identical across 1, 4, 8 and 16 threads.

### A level that cannot fill the pool hands it down

"At most one level consumes the pool" is right for avoiding oversubscription and
wrong for *utilisation* when the outermost level has too few items to occupy the
pool. The conservative sensitivity analysis is the case in the executed
benchmarks: it runs the pseudo-count grid `0.1, 0.5, 1.0` alongside the main run,
so the outermost level has **three** items. On a 16-core host that is three busy
threads and thirteen idle for the whole stage — and because the claim is *held*,
every level inside each refit is also told it may not parallelise, even though each
of those levels has thousands of items.

`NestingBudget::level_for_items(name, n)` therefore takes the pool only when
`n >= rayon::current_num_threads()`; otherwise it declines, `map_par` runs its items
one at a time, and the level inside each item claims the pool and uses all of it.
On `bm5` that is the difference between three occupied cores and sixteen.

This is not nesting — the budget's rule is unchanged, the pool is still split
once. It is only a decision about *which* level gets it, and it only changes
behaviour when the outer level provably cannot fill the pool.

The non-conservative grid (50 refits) has more items than threads, so it keeps the
claim and runs as before; the split is per-grid, not hard-coded.

### `sampling_fractions` is the taxa axis, and it was serial

The plan's innermost level is `taxa`, and `sampling_fractions` is the one stage
that had no `map_par` at all: a serial loop over samples, each summing over every
taxon, with the taxon axis strided because `y1` is column-major. On `bm5` that is
10^8 strided reads, and the stage took **4.0 s at one thread and 4.0 s at
sixteen** — the signature of a stage that never received the pool.

Each sample's `theta[j]` depends on every taxon and on no other sample, so the
sample axis splits into disjoint slots. Parallelised, the same stage is **0.47 s**:
8.6x. It is a small share of the total (4%), but it was a stage that could not
scale at all.

### The group's `used` mask is per group

`fit_one_group` allocates its own "which samples did this fit use" mask rather than
sharing one. It is the only mutable state the group loop would otherwise share, and
it is the thing two workers would race on.

### `--threads` has to be honoured

The analysis runs inside `with_pool`, so the core's `par_iter`s draw from the pool
`install_pool` built rather than from Rayon's own global one. A bare `par_iter`
outside `with_pool` would use Rayon's default at `nproc` width, and `--threads`
would be decorative.

## Group results are carried in tiles

`lm_fit_all`'s groups are handled `GROUP_TILE` at a time. Each group's fitted
values have to travel from its worker to the scatter loop, so holding every group's
at once holds a second copy of the whole `n_taxa x n_samp` table — 100 million
entries on `bm5`, which is the difference between finishing at 29 GB and being
killed.


## `pseudo = 0` turns every zero into a missing value

This is not a performance footnote; it is the single largest structural fact
about the benchmark surface, and it inverts what the missingness-pattern machinery
is for.

`log_center` maps `-Inf` to `NA`, and `observed_mask` is built from the
log-centred response, so with `pseudo = 0` a zero count is *missing*, not zero. A
taxon's missingness pattern is then precisely the set of samples where its count is
zero.

On `bm5` — a 10% zero rate over 20,000 samples — 1,500 sampled taxa have **1,500
distinct patterns**. One taxon each. So:

* the per-pattern QR cache, which PLAN.md item 3 specifies, has nothing to cache
  *for*; each entry would serve a single taxon, and factoring costs
  `O(n_rows p^2)` against `O(n_rows p)` to apply;
* the `structural_zeros` screen still runs, but on a table where the reference
  treats zeros as absences anyway, so it agrees with the prevalence filter.

At ANCOM-BC2's own default `pseudo = 0.5`, `log(0 + 0.5)` is finite, every taxon
is fully observed, there is **one** pattern, and one factorisation serves all
5,000 taxa. Measured on `bm5` at one thread, that is `mle1` 99.0 s -> 21.6 s.

All six benchmark datasets are generated with `pseudo = 0.0`, so the surface runs
in the regime where that optimisation cannot apply. `docs/compatibility.md` has
the full comparison, including the R arm measured the same way.

## What actually costs memory

Every number below is measured on this host (AMD Ryzen 7 3700X, 16 CPUs, 31 GB)
by `crates/ancombc2-cli/examples/allocprof`, which installs a counting *and*
resident-set global allocator. Total allocated bytes say "this run asks for 35 GB",
which does not say which structure is responsible; the peak *live* set by size
class says "2.6 GB of 256 KB allocations are live at once", which localises it
immediately. The gap between the two is also worth knowing: on `bm5` peak live
tracks resident RSS to within 1%, so there is no allocator fragmentation or
arena-retention story to tell — the memory really is in use.

Four structures dominated, in the order they were found:

| structure | was | now | why |
| --- | --- | --- | --- |
| `DesignCache` per-pattern designs + QR reflectors | **5.25 GB** | 0.14 GB | cached per pattern, retained for the whole run |
| `o1`, `o2`, `y1`, `y2`, `y_bias_crt` in `CoreOutput` | 5 x count matrix | opt-in | part of the golden contract, read by nothing else |
| `GroupWrites::fitted` as `(sample, value)` pairs | 1.4 GB | 0.7 GB | a per-value `usize` that `cache.rows[g]` already knows |
| `structural_zeros`' `present` matrix + triple loop | 800 MB + `O(n_tax n_grp n_smp)` | bitset + one pass | re-scanned every sample once per group |

### The `DesignCache` was the wrong shape

`DesignCache::build` factorised every missingness pattern up front and held the
result for the whole run. On `bm5` that is ~1000 patterns at ~18,000 rows each,
so the designs and their reflectors came to 5.25 GB — over thirty times the count
matrix, on the dataset meant to be the memory reference.

It was also a bad trade on time, which is the part worth generalising. Factoring a
pattern costs `O(n_rows p^2)`; *applying* the factorisation to one taxon costs
`O(n_rows p)`. Caching therefore only pays when a pattern carries more than about
`p` taxa. Real tables do not have that: `bm4` and `bm5` each have roughly one
pattern per taxon, because a taxon's missing counts give it a private pattern. The
cache was paying `p` times the work — 18 here — to reuse nothing.

The per-pattern design and factorisation are now built inside the group loop,
which is already parallel over patterns, and freed when the group finishes. Peak
becomes `n_threads` patterns' worth rather than `n_patterns`. The cost is one
extra factorisation per additional MLE iteration; measured, that is +16% on the
`bm5` MLE, against a 37x reduction in that stage's resident set.

### The five intermediates are a contract, not a use

`o1`, `o2`, `y1`, `y2` and `y_bias_crt` are all `n_taxa x n_samp` — five times the
count matrix between them, 4 GB on `bm5`. They are part of the 26-quantity golden
contract ("processed counts, centered Y"), so the parity harness compares them.
Nothing else does: not the CLI, not the FFI, and not the sensitivity analysis,
which re-runs the whole pipeline and would otherwise hold five more *per refit*.
They are behind `AncombcConfig::keep_intermediates`, default off, and empty rather
than zeroed when off so a caller that forgets the flag gets a shape mismatch
instead of plausible-looking zeros.

### Structural zeros: one pass, and a bit set

`structural_zeros` used to build an `n_taxa x n_samp` `f64` presence matrix and
then walk all samples once *per group*, re-reading the count matrix `n_groups`
times. Since a sample belongs to exactly one group the visits are disjoint and the
summation order does not matter, so a single sample-major pass computes the same
tallies: `O(n_taxa * n_samp)` instead of `O(n_taxa * n_groups * n_samp)`, and the
800 MB presence matrix is gone. `crates/ancombc2-core/src/preprocess.rs` keeps the
literal reference triple loop in its test module as an oracle for the rewrite, and
the two are compared across shapes that include empty groups, `NA` counts, `NA`
group labels and both settings of `neg_lb`.

The flags themselves moved from `Vec<bool>` to a `ZeroBitSet`. The screen that
consumes them — "does this taxon have a structural zero in *any* group" — is a
per-row scan, which is what a bit set is for. The `Vec<bool>` stays on
`StructuralZeros` because that is what the CLI writes and the FFI returns; both
serialise it directly, and expanding the bits is `n_taxa * n_groups` byte stores
on a table that is about to be written out anyway.
