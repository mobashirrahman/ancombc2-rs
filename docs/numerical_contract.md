# Numerical contract

What the core guarantees, and how each guarantee is checked.

## 1. Layouts

| array | shape | layout | why |
| --- | --- | --- | --- |
| `CountMatrix::data` | `taxa x samples` | row-major | a taxon's samples are contiguous, which is what the hot loops walk |
| `RMatrix` (`y1`, `y2`, `y_bias_crt`) | `taxa x samples` | row-major | as above |
| `beta_star`, `beta`, `se`, `w`, `p`, `q`, `var1`, `var_hat`, `var_final`, `dof` | `taxa x p` | row-major | as above |
| `vcov`, `vcov1` | `taxa x p x p` | row-major, `p x p` blocks | the block is one taxon's covariance, and its rows are the coefficient axis |
| golden `.f64` blobs | any | **column-major** | R's native order |
| `Matrix` (`x`, the design) | `samples x p` | column-major | matches `model.matrix`, and `matmul` is written for it |

R and the core disagree about layout, so the parity harness bridges it explicitly
with `to_column_major` rather than loosening a tolerance. That function asserts
the element count, so a shape mistake is a named error rather than an
out-of-bounds read — the failure mode that surfaced when the two taxon sets were
first conflated.

Every array crossing the `ancombc2-io` boundary is converted explicitly. The
`CoreOutput` fields are documented one by one with their shape and layout.

**The FFI boundary carries both layouts at once**, and that is a trap worth naming:
`ancombc2_ffi`'s JSON request takes `counts` row-major (`CountMatrix`) and `design`
column-major (`Matrix`), because the two types differ upstream. Sending a
row-major design produces two *identical* columns rather than an error, and the
resulting message — "estimation failed for the following covariates: col0" — points
at the intercept when the real fault is the layout. The R wrapper therefore builds
the design column-major explicitly and `the_ffi_design_layout_is_column_major`
asserts it.

## 2. Missing values

`NaN` is a count's or a log's "not observed", and it is *not* interchangeable with
a number anywhere:

* `observed_mask` treats a non-finite response cell as unobserved, and a taxon's
  missingness pattern is the set of its observed samples. A taxon with all cells
  unobserved keeps `NaN` coefficients and `dof = 999`.
* The sampling-fraction sum uses `na.rm` semantics, which drop the individual
  non-finite *term* rather than the whole row. A taxon with an unfitted
  coefficient therefore still contributes, with that term omitted. Propagating the
  `NaN` instead makes every sample's fraction `NaN`, and `colMeans(..., na.rm =
  TRUE)` cannot rescue a fully non-finite column.
* `Inf` in a result is always a bug and is rejected by the parity harness; `NaN`
  is legitimate (`y_bias_crt` marks a zero count with it, exactly as R's `log(0)`
  does).

## 3. Reduction order and determinism

The result must not depend on the thread count. Two mechanisms:

* **Fixed-order reductions.** Every `par_iter` in the core collects into a
  pre-sized buffer indexed by the item, never by completion order. There is no
  `fold` whose associativity could change, and no atomic accumulation.
* **One global pool**, sized once by `--threads`. Nested parallelism is avoided
  by the nesting budget: pseudo-count runs are the outer `par_iter`, the E-M
  coefficients the next, and the missingness groups the innermost, so an inner
  loop never has more items than there are workers.

Property test P15 asserts that 1 and 8 threads agree, and CI runs the CLI twice
with different `--threads` and diffs the output byte for byte.

Floating-point addition is not associative, so a *parallel* sum can differ from a
*sequential* one in the last bits even with a fixed order. The design keeps every
reduction in taxon order regardless of thread count, which makes the result
bit-identical, and the byte-for-byte CI diff is what proves it.

## 4. Tolerances and the parity levels

The levels are in `PLAN.md` §5. In `crates/ancombc2-core/tests/golden/mod.rs`:

| level | quantities | tolerance |
| --- | --- | --- |
| A | retained taxa and samples, column names, structural-zero flags, `diff_abn`, `passed_ss`, `diff_robust`, pattern assignment | exact |
| B | `beta*`, `theta` | `rtol 1e-8` |
| B | `se`, `delta_em`, `delta_wls`, `vcov` | `rtol 1e-7` |
| B | `s0` | `rtol 1e-9` |
| C | `p`, `q` | `atol 1e-10` + error propagation |
| B | `delta_em`, `delta_wls`, `var_delta` | `rtol 1e-7` |
| D | `diff_abn` concordance | 100% on the small fixtures |

### The fixture matrix, and what it changed

The four committed fixtures are a handful of shapes. PLAN.md §5.5 asks for a
*matrix* over ten axes, and building it
(`reference/R/fixture_matrix.R`, `make matrix`) found five defects that no single
fixture could reach, because each needs a shape, a sparsity, or a design width the
others do not have:

| found by | defect |
| --- | --- |
| `sparsity-090` | the per-taxon `lm` fallback's NA is written into `fitted`, and `theta` is a `na.rm` mean of it — `theta` was 0.75 out |
| `sparsity-090` | `q` had never been compared against the oracle, only against `p_adjust(our p)` |
| `adjust-hochberg` | `Hochberg` and `BH` shared a branch; the weights are `n + 1 - i` against `n / i` |
| `adjust-*` | the `q` re-derivation hardcoded `Holm`, so six of the seven methods were never exercised |
| `int-sparsity90-5group` | `qr` asserted `n >= p` and aborted on a group narrower than it is tall |
| `covariates-10-interaction` | the formula normaliser stripped `*`, so no interaction ever reached the design |

Two of those — the `q` comparison and the hardcoded method — mean the contract was
*weaker than it read* rather than wrong in its conclusions. A tolerance or a check
that silently covers less than it claims is worse than one that is absent, because
it is read as coverage.

Three rules govern how a tolerance may be used:

1. **Never loosen a tolerance to make a test pass.** A divergence is fixed in the
   algorithm or documented as a scope limit. Every tolerance change in the history
   of this file came with a change to the implementation.
2. **A derived quantity's tolerance is set by its operands, not by itself.**
   `samp_frac` is a difference of quantities known to `1e-8`, so it is compared at
   `TOL_BETA.with_operand_scale(y1)` rather than at a tolerance keyed to its own
   small magnitude. `W = beta / se` likewise. Comparing a ratio at a tolerance
   keyed to the ratio would pass a computation whose inputs were each wrong.
3. **An exemption is not a relaxed tolerance.** The one place the contract is
   scoped rather than held — the rank-deficient class — is enumerated, counted,
   capped and printed. See `docs/reference_behavior.md` §10.

### What Level C actually achieves

The `1e-10` in the table above is a *floor*, and the check that runs is the floor
plus the **observed** discrepancy in `W`, because a p-value is
`2 * (1 - F_t(|W|))` and `|dp/dW| = 2 f_t(|W|) <= 0.8` — an error in `W`
propagates into `p` at no more than unit gain. That is error propagation, not a
loosened bound, and it is what lets the floor be the plan's `1e-10` while the
largest `p` disagreement seen is larger than `1e-10`.

`ANCOMBC2_GOLDEN_AUDIT=1` records the largest `|rust - oracle|` for every
compared quantity, so the tolerance is a measurement rather than an assumption.
Over the three fixtures where parity is asserted — no rank-deficient taxon:

| fixture | shape | `p` | `q` |
| --- | --- | --- | --- |
| `fx01` | 10 x 10 | 2.6e-15 | 0 (exact) |
| `fx02` | 100 x 30 | 7.6e-13 | 1.1e-13 |
| `fx03` | 1000 x 100 | 1.6e-10 | 2.0e-10 |

`fx04` reaches 2.0e-4 for `p` and 7.0e-4 for `q`, with 200 of its 10,000 taxa on
an exactly singular sub-design. Those quantities were once reported rather than
asserted for that reason; they are asserted now, and the remaining deviation is
inherited from the E-M rather than new — see `docs/reference_behavior.md` §16 for
why the premise behind the exemption was false. `fx03`'s 1.6e-10 traces to
`delta_em` at 1.7e-11, which is
inside its own Level B contract — the residual is summation order inside the E-M
sweep.

This constant was `atol = 1e-8` until the audit was added, justified by the E-M
iteration tolerance of 1e-5. That reasoning was wrong and is worth recording: the
E-M tolerance says how accurately the *algorithm* converges, not how closely two
implementations of the same algorithm agree when both run it to the same stopping
rule. The measurement is two orders of magnitude tighter than the value it
replaces, and the audit that produced it is itself covered by tests, so the
numbers in this section cannot silently rot.

`beta_star`, `theta` and `var1` are reported rather than asserted whenever the
rank-deficient class is present, for the reason in
`docs/reference_behavior.md` §13: which coordinates the reference's per-taxon `lm`
reports is decided by its LAPACK pivoting, and `theta` is a mean over those
coordinates. The class is counted in three separate buckets — singular design,
`lm` aborted, and under-determined (`n_obs < p`) — each with its own cap, because
they have different causes and different correct rates.

### Two checks on `q`

`q` is checked twice, because either check alone is insufficient:

1. **Re-derivation.** `p_adjust` is applied to the *golden* p-values and must
   equal the *golden* `q` exactly (Level A). This tests the adjustment in
   isolation, and it is what catches a wrong adjustment method.
2. **The pipeline.** The Rust `q` is compared against the *golden* `q`, and
   separately against `p_adjust(our p)`.

   `q` needs a **larger** error-propagation term than `p`, and that is not a
   loosened bound. `q = p_adjust(p)` is not 1-Lipschitz: Bonferroni multiplies by
   `n`, and BH/BY divide by a rank while taking a running minimum, so a
   perturbation of `dp` moves `q` by up to `n * dp`. The check uses the worst case
   any of the seven methods can produce, `atol + n * observed_dp`. Using `p`'s
   tolerance for `q` is not a tighter check, it is a wrong one — on the
   `sparsity-010` matrix cell `p` agrees to 3.5e-11 and `q` to 3.5e-9, exactly a
   factor of the 100 taxa.

Check 2 against the golden `q` is the one that was missing. Comparing our `q` with
`p_adjust(our p)` alone proves the pipeline is internally coherent, but a `q` that
is a deterministic function of a wrong `p` passes it. It is listed in
`INDIRECT_QUANTITIES`, so under a rank-deficient fixture it routes through
`compare_indirect` and is reported rather than asserted, exactly as `p` is.

## 5. Failures report the first diverging quantity

`compare_core` walks the quantities in pipeline order and returns the first
divergence, with the level, the flat index, the taxon name, the coefficient name,
both values, and the absolute and relative difference. The order is the pipeline
order, so a sandwich divergence is reported as a sandwich divergence rather than as
a downstream p-value mismatch.

## 6. The linear algebra

There is no BLAS dependency. `Matrix::matmul`, `qr`, `cholesky_solve`, `ginv` and
`eigen_symmetric` are in `crates/ancombc2-core/src/matrix/linalg.rs`, and they are
covered by the Miri job, because out-of-bounds indexing in a numerics crate is
silent in release and catastrophic in the worst case.

`qr` is a Householder factorisation with **column pivoting on the remaining
2-norm**, and the rank test is `lm.fit`'s: a diagonal counts when it exceeds
`1e-7 * |R_00|`, counting the *leading run* and stopping at the first failure.

Two of those three choices are load-bearing and were found by the parity suite:

* **`tol = 1e-7`**, `lm.fit`'s documented default, not LAPACK's
  `dlamch("epsilon")`. They differ by nine orders of magnitude and change which
  columns count as aliased on an ill-conditioned design.
* **Stop, do not skip.** A later diagonal can be large again after an exact zero;
  counting it reports a full rank for a design whose middle column duplicates an
  earlier one, and the fit then produces coefficients of order `1e14`.

Pivoting is on because the alternative is not an option: the diagonal magnitudes
of an unpivoted Householder QR depend on the reflector formula, so two correct
implementations disagree on the rank of the same ill-conditioned matrix. The
cost is that *which* column is dropped on a rank-deficient design is not
reproducible against the reference; that is the documented scope limit.

`ginv` is used only where the reference uses `MASS::ginv`. The sandwich's
`(X'X)^-1` is computed once per run, not per taxon, because it does not depend on
the taxon.

## 7. The sandwich quirk

`.sandwich_vcov` writes `XX[rowSums(XX, na.rm = TRUE) == 0] = 0.1`, so a missing
entry of a term contributes `0.1` to the sum rather than being dropped.
`CompatMode::Ancombc2_15` reproduces it and `CompatMode::StrictSpec` skips the
term. `CompatMode` reaches the kernel through
`AncombcConfig::compat → iter_mle_* → sandwich_all`; property test P16 is the
regression test for that wiring, and it was written *after* the wiring was found
to be broken — the configuration field existed and was silently ignored.

## EM mixture parameters

The three-component Gaussian mixture that `.bias_em` fits, captured per fixed
effect. Like the convergence trace, it is computed by the reference but not
returned: `.bias_em` returns only `c(delta_em, delta_wls, var_delta)`.

It is captured by tracing the **exit** of the oracle's own `.bias_em` and reading
the loop's final `*_new` values — the parameters it actually used. Re-implementing
the routine in the harness instead would be a ~150-line transcription of a
Nelder-Mead EM fit into the one number the whole method rests on, where a
transcription slip would be a silent parity risk that no assertion could see.
Reading them off the function that produced them cannot drift from the oracle.

Per term, the mixture carries `pi` (three component weights), `l` (`l1`, `l2` — the
component offsets, constrained `l1 <= 0 <= l2` by `min`/`max`), `kappa` (`kappa1`,
`kappa2` — the component variance offsets) and `delta`.

| fixture | terms | iterations |
| --- | --- | --- |
| `fx01` | 2 | 14, 100 |
| `fx02` | 3 | 27, 100, 39 |
| `fx03` | 5 | 100 ×5 |
| `fx04` | 8 | 24, 100, 100, 100, 100, 74, 78, 75 |

Two properties worth recording:

* Many fits stop at `max_iter = 100` rather than at `tol = 1e-5`, including every
  term of `fx03`. The bias estimate is the reference's own output in that case, so
  parity requires reproducing the same cap, not converging harder.
* `kappa` is frequently exactly `0`. `nloptr` optimises it with `lb = 0` and the
  likelihood is monotone in the variance offset at the boundary for a small taxon
  set, so the optimum sits on the bound. It is a real fitted value, not a
  placeholder.

Compared at Level B `rtol 1e-7`, with the iteration count exact as for the
convergence trace. `delta` is held to `rtol 1e-7` and not tighter because it is the
same number the `.f64` contract carries as `delta_em`, and it cannot be held
exactly: R sums the E-step accumulator in a single `sum()` over a vector this
implementation accumulates in blocks, so the two agree to about four units in the
last place (`-0.11310136764609605` against `-0.1131013676460965` on `fx01`).
Bit equality here would demand more than the existing `delta_em` parity gets. The
*identity* of the fit is pinned on the oracle side instead, where
`reference/R/harness.R` asserts the traced `delta_new` equals the `delta_em` the
run returned — an exact equality between two values from the same function, which
is what rules out a trace having captured the wrong call.

Both failure modes are exercised: an altered `pi` and an altered iteration count
each fail `fx01_tiny_two_group` with the specific message.

## Convergence trace

The contract's "convergence trace" is the first MLE's `epsilon` at each
iteration. The reference **prints** it — `ML iteration = k, epsilon = e`, via
`message()` when `verbose = TRUE` — and returns only the final epsilon, so the
oracle's half is captured from that printed trace and the Rust half from the loop
(`IterMle::trace`, surfaced as `CoreOutput::ml_trace`).

Observed traces:

| fixture | epsilons | iterations |
| --- | --- | --- |
| `fx01` | 2.4, 8.8e-16 | 2 |
| `fx02` | 2.3, 0.043, 5.9e-05 | 3 |
| `fx03` | 4.7, 2, 0.012, 9e-05 | 4 |
| `fx04` | 11, 17, 0.28, 0.0052 | 4 |

Two levels of check, and the second is not the usual one:

* **Level A, exact** on the iteration count. Two runs taking a different number of
  steps converged differently, and no tolerance on the values would make that
  comparable.
* On the values, **at the precision the oracle publishes**. The reference prints
  `signif(epsilon, 2)` — two significant figures — so the check is that this run's
  epsilon rounds to the same two figures. That is an exact comparison at the
  available precision, not a hand-picked tolerance.

Holding the trace to Level B's `rtol 1e-8`, which the rest of the numerics is held
to, would look stricter and be wrong: the oracle's own printed `2.4` is
`2.3946148989178186` as computed, so an `rtol 1e-8` assertion would fail against an
exact reimplementation. The reference does not carry more precision here.

Below `tol * 1e-12` an epsilon is accumulated double rounding rather than a
property of anything — the reference's final iteration reports `8.8e-16` where an
exact reimplementation computes `6.6e-16`, from the same arithmetic. Those
iterations are checked only for still being under the tolerance, which is the
property that holds. The floor is derived from the tolerance rather than chosen.

Both failure modes are exercised: an altered epsilon and an altered iteration
count each fail `fx01_tiny_two_group` with the specific message.

## `s0` is bounded by the E-M, not by itself

PLAN.md §5 specifies `rtol 1e-9` for `s0` and `rtol 1e-7` for
`se`/`delta_em`/`delta_wls`/`vcov`. `s0` was being compared at `1e-7` — two
orders looser than the contract, and passing. Tightening it to `1e-9` exposed a
real conflict between the two numbers.

`s0` is `quantile(var_hat[, k], s0_perc)` over R type 7, and the column it
quantiles is

```text
var_hat = var + var_delta + 2*sqrt(var*var_delta) + s0
```

so `s0` inherits `var_delta` and the E-M's `pi`/`kappa` estimates through two
terms. `var_delta = 1/sum(1/nu)` is a **single scalar per coefficient** — there is
no averaging-out over taxa, so nothing dilutes an error in it. The bound `s0` can
actually meet is therefore the E-M's bound.

That only bites when the E-M does not converge. Measured on the
`predictor-5group` matrix cell (6 terms, `max_iter = 100`, `tol = 1e-5`):

| term | E-M iterations | `delta_em` relative deviation | `s0` relative deviation |
| --- | --- | --- | --- |
| `(Intercept)` | 44 | 2.7e-12 | 9e-16 |
| `group2` | 57 | 3.7e-16 | 2e-16 |
| `group3` | 53 | 4.7e-12 | 2e-16 |
| `group4` | 100 (capped) | 4.8e-14 | 4e-16 |
| **`group5`** | **100 (capped)** | **2.2e-9** | **1.2e-9** |
| `x1` | 92 | 1.2e-15 | 2e-16 |

One term out of six. `delta_em` for `group5` is 2.2e-9 from the oracle — inside
the E-M's own `rtol 1e-7` by four orders of magnitude — and that term is the one
bracketing the 5% quantile, so `s0` inherits it. `group4` also hit the cap and
agrees to 5e-14; hitting the cap makes divergence *possible*, not certain.

When an iteration stops at `max_iter` instead of at `tol`, the returned value is
wherever the iteration happened to be — a continuous function of the trajectory.
Two implementations following the same stopping rule land 1e-9-ish apart there.
This is the same phenomenon already documented for `delta_wls` on the
`covariates-10` cell, arriving through `s0` instead.

### How the contract handles it

`s0` uses the same rule as the mixture: `em_tolerance_for`, which widens a term's
bound to the E-M's **recorded final epsilon** — the size of its last parameter step
— when that is coarser than the contract, and otherwise leaves it alone. Every
comparison prints the bound it applied, so the slack is visible on each run rather
than being an absence of evidence.

Two earlier rules were tried here and both were wrong:

* **the iteration count.** A fit that stopped at 61 iterations with a final epsilon
  of 9.1e-6 has "converged" by the recorded test, and still leaves `s0` at
  1.07e-9 on the `covariates-10-interaction` cell's third coefficient — outside
  `rtol 1e-9`. Converged is not the same as accurate.
* **the achieved `delta_em`.** Measured, which made it the obvious candidate, but it
  did not predict that case: `delta_em` for the term in question agrees to better
  than 1e-9 while `s0` does not, so the disagreement is downstream of the bias
  estimate and a bound keyed on it never engaged at all.

The final epsilon over-states the error by roughly three orders on that cell
(9.1e-6 recorded, 1.07e-9 observed), so the bound it produces is loose. That is
the trade being made deliberately: a loose bound that is derived from the recorded
run beats a tight one that is a guess about which fixtures are hard, and a bound
that tightens automatically as the E-M converges is still doing real work on the
terms where it matters.

`em_convergence_is_read_from_the_recorded_mixture` pins the parsing, and
`the_s0_bound_follows_the_em_convergence` pins that the bound moves with the
recorded state rather than with a maintained list of cell names.

## Per-stage timings

The contract's last quantity is the oracle's per-stage wall time, recorded by
`reference/R/harness.R` into `stage_seconds.json` beside the canonical files for
each fixture.

**They are recorded, not compared.** Wall-clock is not reproducible: asserting
that the oracle spent 44.3 s in `core` would assert that it ran at the same speed
on a different machine, which is not a property of this contract. The parity suite
checks the things that can actually break, and only those:

* every oracle stage name is recognised, so a stage that silently stopped being
  timed is caught;
* every recorded value, on both sides, is finite and non-negative;
* the Rust side produced timings at all, so the check cannot pass vacuously.

Both failure modes are exercised rather than assumed: injecting a negative duration
and injecting an unrecognised stage name each fail `fx01_tiny_two_group` with the
specific message.

The stage granularity differs by design. The oracle harness times five regions
(`sanity_check`, `structural_zeros`, `preprocess`, `core`, `sensitivity`); the Rust
`StageTimings` splits the core into `preprocess`, `pattern_grouping`, `mle1`,
`sandwich1`, `em`, `correction`, `mle2`, `sandwich2`, `tests` and
`serialisation`, plus `sensitivity`. The check accounts for each oracle region
against the Rust set rather than requiring the two lists to be equal, because the
difference is instrumentation, not algorithm.

The measurement is worth having for its own sake: on `fx04` it shows the oracle
spending **96% of its time in `sensitivity`** (1065.6 s of 1110.3 s) and 4% in
`core`. The large datasets are a comparison of two implementations of a repeated
whole-pipeline refit, which is not obvious from a wall-clock number alone.

