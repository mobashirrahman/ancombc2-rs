# Preprocessing and reduction, stage by stage (IMPROVED_PLAN.md S08)

What preprocessing computes, what R actually computes at each step, and the nine
places the Rust and R versions disagreed before this task. Every claim below is a
measurement from the primary profile, and the measurements are reproducible with
`make preprocess-stages`.

The claim under test is narrow: **for the same input, every preprocessing stage the
bridge computes is bit-identical to the stage R's own expressions produce.** That is
not Rust numerical parity — nothing here computes a fit — but it is the thing that
has to be true before any of it can be.

## 1. Why stages, and not results

`.ancombc2_core()` has one return value. Comparing it against the reference's tells
you *that* something differs; it cannot tell you whether the filtering, the log, the
centring or the covariance was at fault, and by the time the difference reaches the
result it has been through a least-squares fit, an EM iteration and a sandwich
estimator.

So the stage probe runs preprocessing in the bridge and hands **each stage's array**
back separately: `prevalence1`, `tax_keep1`, `lib_size1`, `samp_keep1`, `O1`, `log1`,
`means1`, `y1`, and the same eight for the aggregate pass, plus `zero_keep`,
`group_prevalence`, `group_size`, `group_lower` and `zero_ind`. Twenty in all.

`scripts/check_preprocess_stages.R` then runs the *reference's own expressions* —
transcribed from `ancombc_prep.R:2-84` and `:110-130`, not imported from the
installed original, so a divergence cannot hide inside the package under test — and
compares with `identical(serialize(x, ...), serialize(y, ...))`.

**274/274 stage comparisons pass** across thirteen fixtures: tiny, sparse,
missing-value (double and integer), boundary-prevalence, boundary-library,
boundary-zero-pseudo, one-sample, one-group, no-group, aggregate-differs, neg-lb and
neg-lb-off.

## 2. The nine differences

Each of these was found by the comparison, not by reading the code. None of them
produces a wrong *number*; every one of them is invisible to `all.equal()` and to
every numeric check, and visible only to `serialize()`.

### 2.1 `f64::ln` does not preserve a NaN payload

`log(NA_real_)` is `NA_real_` in R, because R's `R_log` begins `if (ISNAN(x)) return
x;`. `f64::ln` calls libm, and neither C nor IEEE 754 requires a NaN payload to
survive `log`. On this platform:

```text
f64::ln(0x7ff00000000007a2)  ->  0x7ff8000000000002
```

A different NaN. `is.na()` is true of both, every arithmetic result agrees, and
`serialize()` disagrees — on every `NA` count in the table. [`r_log`][reduce] is now a
one-line `if x.is_nan() { x } else { x.ln() }`.

### 2.2 R's `is.infinite` mask is not "not finite"

`o[is.infinite(o)] <- NA` replaces `±Inf` and **not** a NaN of any payload. The first
version mapped every non-finite value to `f64::NAN`, which both destroyed an existing
`NA_real_`'s payload and rewrote a computed `NaN` into `NA_real_`. Now
[`replace_infinite_with_na`][workspace] keys on `is_infinite` and writes
[`NA_REAL_BITS`][reduce].

### 2.3 `NA_real_` is a *signalling* NaN, and arithmetic quiets it

```text
NA_real_                7ff00000000007a2   quiet bit clear
NA_real_ - 1.5          7ff80000000007a2   quiet bit set, payload kept
1.5 - NA_real_          7ff80000000007a2
NaN - 1.5               7ff8000000000000   already quiet, so unchanged
NaN + NA_real_          7ff8000000000000   the *left* operand decides
-NA_real_               fff00000000007a2   a genuine sign flip, not arithmetic
```

Read `7ff0` → `7ff8` carefully: that is IEEE's **quiet** bit, bit 51, not the sign.
Arithmetic never touches the sign. The rule is `bits | QUIET_BIT`, it keeps the
operand's payload, it is idempotent, and where both operands are NaN the left one
decides. This is [`r_nan_arith`][reduce], and it applies to `O = data + pseudo` as
much as to the centring — so an `NA` count is already quieted by the time `log` sees
it.

### 2.4 …but the *hardware* loses the payload entirely

```text
na_real() - 1.5   ->  0x7ff8000000000000     (payload cleared)
R: NA_real_ - 1.5 ->  0x7ff80000000007a2     (payload kept)
```

Eighteen bits, no numeric consequence, and a `serialize()` difference. This is the
whole reason for #2.3.

### 2.5 A unary math function and a binary operator have *different* NaN rules

`log(NA_real_)`, `sqrt(NA_real_)` and `abs(NA_real_)` all return the operand
**unchanged**; `NA_real_ - 1.5` does not. One "propagate the NaN" helper covering both
would be wrong for one of them, so there are two — [`propagate_nan`][reduce] and
[`r_binop`][reduce] — with the measurements written next to each.

### 2.6 `rowSums`/`rowMeans` accumulate in `long double`, not `double`

R's `do_colsum` declares its running total as C `long double`. On x86-64 that is the
x87 80-bit format: **64 bits of mantissa**, not 53.

```text
y = c(1e16, 1, -1e16, 1e-17, 1)
rowMeans   0x1.999999999999ap-2   (sum 2, /5)
f64 loop   0x1.999999999999ap-3   (sum 1, /5)
```

A factor of two, from one line of `+=`. A hand-rolled 80-bit accumulator would be
right on x86-64 and wrong on aarch64 (binary128) and on Windows (`double`), so it
would encode one platform's answer as the oracle's. A 106-bit compensated sum is not a
substitute: it rounds differently at every step.

So the reduction is asked of R: [`RBackedReductions`][bridge] calls
`base::rowMeans`/`colMeans`/`colSums`/`rowSums` across the boundary. Four calls per
run, none of them in a hot loop.

The three reductions preprocessing does over *counts* — `rowSums(x != 0)`,
`rowSums(!is.na(x))`, `colSums(feature_table)` — are exact in `f64`, because every
partial sum is a whole number below `2^53`. Those stay plain `f64`, and
[`sum_counts`][reduce] asserts the invariant rather than trusting the caller.

### 2.7 The library size is over the *retained* taxa, and over *all* samples

`.data_core` subsets the taxon axis, then computes `colSums`, then subsets the sample
axis. So `lib_size` is the column sum of the taxon-filtered table with every sample
still present. The first version read it off the sample-filtered table: the right
*decision* (the filter used the full table) and the wrong *array*, which is worse,
because the array is what a later stage would consume.

### 2.8 The structural-zero screen reads the **aggregate** table

`ancombc2.R:453-462` runs `.get_struc_zero` on `feature_table_aggregate` against the
*unfiltered* `meta_data`, before either `.data_core` pass. The two tables have the
same shape and different contents, so screening `feature_table` instead would produce
a different `tax_keep` — silently, because every shape still matches. The
`aggregate-differs` fixture exists only to catch this.

### 2.9 `.data_core`'s second pass reports subset-relative indices

```r
feature_table = feature_table[tax_keep, , drop = FALSE]
prevalence = prevalence_fun(feature_table)
tax_keep = which(prevalence >= prv_cut)     # counts within the subset
```

So `core2$tax_keep` is not in `aggregate`'s coordinates. Nothing downstream reads it —
the reference only uses it to subset again — but a stage array named after it has to
hold what it holds. Both are reported: `tax_keep2` as the reference reports it, and
`tax_keep2_absolute` for a Rust caller.

### 2.10 …and three more that are not arithmetic at all

* `colSums` returns a **named** vector, so `lib_size1` carries the sample names;
  `rowMeans` does not, so `means1` carries none. Getting either backwards is invisible
  to a numeric comparison and visible to `serialize()`.
* `Rf_findNamespace` is not in `Rinternals.h` — only the internal `R_FindNamespace` is
  — so the reduction wrappers are passed *from* R as arguments rather than looked up
  by name from C. The first attempt resolved the name and segfaulted on the first call.
* A `n x 1` stage is still a matrix. The reference's
  `feature_table[tax_keep, , drop = FALSE]` keeps its `dim`, so a shape is emitted for
  every declared shape, not only for those with more than one column.

## 3. Two ways to cross the boundary wrongly, both found here

Both are the same mistake as the one S06's acceptance test exists to prevent, in a new
place.

**Pointers into Rust memory.** The first version of the probe filled a
`[RawStage; 20]` with `&Vec<f64>` pointers and let the C shim read them after
`ancombc2_rb_preprocess_probe` returned. The Rust function's locals are dropped on
return, so every pointer dangled: nothing crashed, nothing was `NULL`, every shape was
right, and the first two elements of each stage came back as whatever the allocator
had left there. Now the values are copied into buffers the caller allocated and only
offsets cross.

**Sizing the caller's buffer from the wrong side.** `out_f64_len` was never written, so
C allocated a one-element `REALSXP` and Rust wrote a hundred-odd doubles past the end
of it. The symptom was not a crash but `"recursive gc invocation"` some allocations
later. Both lengths are now reported before the copy, on both passes, and the C shim
checks them against what it wrote.

A third, smaller one: the payload buffers are real R vectors rather than `R_alloc`,
because the bridge calls back into R for the row means and holding an `R_alloc` buffer
across an R allocation is not safe.

## 4. Closing the gap the reduction opened

Threading `Reductions` through the pipeline turned up a constraint worth writing down,
because it was not obvious until the compiler said so: **`dyn Reductions` is not
`Sync`, and R cannot be called from a thread the main thread did not start.** The
sensitivity analysis refits once per pseudo-count *in parallel*, so a reduction inside
that loop is either unsound or serial.

It is avoidable, and the fix is structural rather than a comment. The centring means
depend only on the taxon subset and the pseudo-count, both fixed before the loop, so:

```text
calling thread   for each pseudo-count:  means = log_row_means(...)     <- may call R
                                          cost: grid.len() * taxa.len() doubles
worker threads   for each pseudo-count:  y = log_rows(...)               <- a ln per cell
                                          y -= means[k]
```

Three things follow, and each is pinned by a test in `reduce::split_tests`:

* the two routes are **bit-identical**, compared with `to_bits()` rather than `==` --
  `==` fails on every `NA` cell, which are the cells that matter;
* the test uses a reduction that is deliberately *not* `f64` (`Shifted`, which adds
  0.125 to every mean) so it distinguishes "carried the right numbers" from "agreed by
  coincidence";
* `map_par` now passes the **index** as well as the item, because the pseudo-count grids
  repeat entries and a mean looked up by value would be ambiguous.

The entry points say which accumulator they use, out loud:

```rust
ancombc2_run(...)            // f64 wrapper, for the CLI and the simulation harness
ancombc2_run_named(..., red) // the reduction is a parameter, never a default
```

A default argument would have been the smaller change and would have hidden exactly the
thing that has to stay visible. `ancombc2-io`, `ancombc2-ffi` and `ancombc2-sim` each
pass `&F64Reductions` with a comment saying why.

## 5. What this does not establish

* **No Rust numerical parity.** Nothing in preprocessing is a fit. The reductions are
  R's; the arithmetic is checked stage by stage. C3 is still open.
* **The pipeline still uses `f64`.** The threading above makes it *possible* for the R
  session to pass an R-backed reduction, and `ancombc2_run_named` takes one -- but
  nothing in R calls the pipeline yet, so `ancombc2-io` and the CLI still pass
  `F64Reductions`. S09 makes that call.
* **`theta_hat` is absent.** It needs `beta1`, which is S09.
* **The design matrix is an input, not a stage.** `x = model.matrix(...)` arrives over
  the S06 transport and is not recomputed here.
* **Thirteen fixtures, not the 114-case campaign.** These are hand-built to hit the
  boundaries S08 names — tiny, sparse, missing-value, boundary. They are deliberately
  small enough that a difference is readable in a hex dump; the 114-case campaign is
  what checks that nothing else moved.

## 6. Reproducing

```sh
make replacement-build
make preprocess-stages       # 274/274 stage comparisons
make r-bridge-selftest       # S06 regression, 32 checks
make r-output-selftest       # S07 regression, 43 checks
```

Evidence: `validation/exact/evidence/s08_preprocess_stages.txt`.

[reduce]: ../crates/ancombc2-core/src/reduce.rs
[workspace]: ../crates/ancombc2-core/src/workspace.rs
[bridge]: ../crates/ancombc2-rbridge/src/output.rs
