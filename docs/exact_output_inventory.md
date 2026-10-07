# The output inventory, and the output transport (IMPROVED_PLAN.md S07)

What `.ancombc2_core()` produces, what shapes and types it has, and how each of
those travels across the native bridge and back. This is the mirror of
[`exact_seam_inventory.md`](exact_seam_inventory.md), which does the same for the
inputs.

The measured numbers come from `validation/exact/payloads/tiny-defaults.rds`,
captured by `scripts/capture_core_payloads.R` from the pinned original's
`.ancombc2_core()` with `n_tax = 5`, `n_samp = 8`, `p = 3`. **Those payloads are
test fixtures. A candidate run may not read them, and using them is not evidence
of numerical parity** — they contain the original's numbers, so agreeing with them
proves only that the transport loses nothing.

## 1. The inventory

Ten of the twelve internals are reachable from `.ancombc2_core()`'s return value;
`O1` and `x` are inputs to the result rather than outputs, and `O1` in particular
is only needed by the bias-estimation side.

| object | R type | shape | in the returned result? | names carried |
| --- | --- | --- | --- | --- |
| `beta_hat` | `REALSXP` | `n_tax x p` | via `res$lfc_*` | `dimnames`: taxa, fixed effects |
| `var_hat` | `REALSXP` | `n_tax x p` | via `res$se_*` | `dimnames`: taxa, fixed effects |
| `dof` | **`INTSXP`** | `n_tax x p` | no | `dimnames`: taxa, fixed effects |
| `vcov_hat` | `VECSXP` of `REALSXP` | list of `n_tax` `p x p` matrices | no | each element's `dimnames`: fixed effects, fixed effects |
| `delta_em` | `REALSXP` | length `p` | yes, verbatim | **none** |
| `delta_wls` | `REALSXP` | length `p` | yes, verbatim | **none** |
| `var_delta` | `REALSXP` | length `p` | no | **none** |
| `s02` | `REALSXP` | length `p` | no | fixed effects (`apply()` named it) |
| `y_bias_crt` | **`VECSXP`, classed `data.frame`** | `n_tax x n_samp` | yes, as `bias_correct_log_table` | row names: taxa; column names: samples |
| `theta_hat` | `REALSXP` | length `n_samp` | yes, as `samp_frac` | samples |
| `O2` | `REALSXP` | `n_tax x n_samp` | yes, as `feature_table` | `dimnames`: taxa, samples |

Five details in that table are the reason this file exists. Each was found by
capturing the original rather than by reading it, and each would produce a result
that is *nearly* right.

### `dof` is integer, and that has to survive

`dof` is the only integer among the numbers. The original computes it from
`sapply` over degrees of freedom, and R's `sapply` simplifies a length-`p` integer
result to `INTSXP`. `serialize()` records that type, so a `dof` that came back as
`REALSXP` fails the byte comparison even with every value identical.

The wire carries one numeric representation — `f64` — so the type is recorded
separately as `dof_is_int` and the plan declares `dof` as `INTSXP` only when the
original's was. The conversion is checked rather than trusted: `ancombc2_rb_emit`
refuses a non-integral value in an entry the caller declared integer, rather than
truncating it.

### `var_hat` is transported *after* the original finishes with it

`ancombc_prep.R:232-236` runs

```r
s02 = apply(var_hat, 2, function(x) stats::quantile(x, s0_perc, na.rm = TRUE))
var_hat = t(t(var_hat) + s02)
var_hat[is.na(beta_hat)] = NA
```

and writes `var_hat[i, ]` into the diagonal of `vcov_hat[[i]]`. The captured
`var_hat` and `vcov_hat` are therefore the **post**-transformation values, and
`sqrt(var_hat)` is already `se_hat`. The first version of `R/assemble.R` re-applied
the `s02` addition, which inflated every standard error and shifted every p-value
— and would still have looked right on any case with `s0_perc = 0`. The payloads
are a capture of the core's outputs, not of its inputs.

### `vcov_hat` is a list, not a matrix

`.sandwich_vcov` builds it with `lapply` and gives each element `dimnames`. It is
transported as one flat `n_tax * p * p` block and rebuilt as a list. Collapsing it
to a single `n_tax*p` x `p` matrix would change `typeof()` and every attribute, and
the byte comparison would fail on the class alone.

The wire is **row-major within each taxon** while R is column-major, so the
transpose happens in exactly two visible places: `rb_flatten_vcov()` in
`R/bridge.R` and the list rebuild in `C_ancombc2_rb_emit`. A covariance is
symmetric, so a transposed covariance is still a covariance and passes every
symmetry check; the acceptance test therefore compares a deliberately
non-symmetric `p = 3` block element by element rather than relying on shape.

### `y_bias_crt` is a data.frame on the way out and a matrix on the wire

The original returns it as `bias_correct_log_table`, a `data.frame` whose row names
are the taxa. It travels as a matrix, so the round trip is `as.matrix()` on the way
in and `data.frame()` on the way out, with both margins taken from the restored
`dimnames`. Its rownames are **not** reset: `res` is the table the original sets to
`NULL` at line 266, and `bias_correct_log_table` is not.

### The names are not uniform

`theta_hat` is named by the samples and `s02` by the fixed effects, because
`apply()` gives it `var_hat`'s column names — while `delta_em`, `delta_wls` and
`var_delta` carry no names at all, and both of the first two are returned verbatim
in the result. A rule of the form "vectors get the column names" would add names to
three objects that must not have them.

`rb_restore_names()` therefore writes out the rule per entry rather than inferring
it from a shape.

## 2. `taxa` and `taxa_bias` are two different sets

`struc_zero = TRUE` makes the original flag taxa whose bias-corrected value is zero
or negative. Those taxa are dropped from the *reported* table (`O2`) but retained in
`O1`, because the sampling fractions `theta_hat` are estimated from `O1`.

So:

- `n_tax` is `nrow(O2)`, the **reported** taxa.
- `taxa` has `n_tax` entries: indices into the *input* table.
- `taxa_bias` has `nrow(O1)` entries and can be **longer** than `n_tax`.

Neither is an index into `O2`. Collapsing them into one vector would report the
wrong number of taxa on every case with a structural zero, and would be invisible on
every case without one.

## 3. The transport

`ancombc2_rb_emit` is `rb_probe`'s mirror. It computes nothing: it copies the core's
payloads out of R's memory, into the bridge, and back into R's memory, and the
acceptance test is that the bytes either side are equal.

```
R payload  --(C_ancombc2_rb_emit_plan)-->  lengths
            --(C_ancombc2_rb_emit)------->  one flat f64 region, one flat i32 region
            --(rb_restore_names)-------->  dimnames, names
            --(ancombc2_assemble_core)->>  the original's result
```

The plan is asked for before anything is allocated, so the allocation and the plan
cannot disagree about a length. `r/ANCOMBC/src/init.c` keeps its own literal table
of the eleven shapes and compares it against the bridge's plan entry by entry: two
sources of truth is the point, and a disagreement is a loud error rather than a
quiet reshape.

`RawPayloads` in `crates/ancombc2-rbridge/src/abi.rs` and `struct rb_payloads` in
`r/ANCOMBC/src/init.c` are the same 216-byte layout in two languages, and nothing in
either type system ties them together. Every field offset is pinned by
`abi::output_layout_tests::the_wire_layout_is_pinned`, so moving a field is a
build failure instead of a misread. This is not hypothetical: an earlier version had
`n_taxa_bias` in a different position in each, and the low bits of a pointer were
read as a count, producing an 812 GB allocation and an abort. The counts are now
also bounded, so an impossible count is a named error rather than an allocation.

## 4. What is verified, and what is not

`make r-output-selftest` — 43 checks, all passing. Among them:

- each of the eleven payloads is `serialize()`-identical to the original's, `dof`
  included;
- the covariance list survives as a list, with each element's `dim` and `dimnames`,
  and a non-symmetric `p = 3` block comes back with `(1,2)`, `(2,1)` and `(1,3)` in
  their original places;
- `NA_real_`, `NaN`, `Inf`, `-Inf` and `-0.0` all survive, and `-0.0` keeps its sign
  bit — checked with `1/x`, because `-0.0 == 0.0` and equality cannot see it;
- `a reported set shorter than the bias set` round-trips, and a *shorter* bias set
  is refused;
- the assembled `.ancombc2_core()` result is byte-identical to the original's, all
  ten names and four `NULL`s included.

**Not** verified, and not claimed:

- The bridge computes nothing. The numbers in the fixture came from the pinned
  original. This is transport parity, not Rust numerical parity; C3 is still open.
- `global`, `pairwise`, `dunnet` and `trend` are not transported. Asking for one is
  a named error, not a `NULL`. Their assemblies need `group`, `para2$fits` and the
  LME4 fits, which S08/S09 supply.
- The single fixture is `tiny-defaults`, where all four test switches are `FALSE`.
  A case with `struc_zero = TRUE` would exercise section 2 end to end; the
  `taxa`/`taxa_bias` handling is currently covered by a synthetic case.

## 5. Reproducing

```sh
make replacement-build
make r-output-selftest          # 43 checks
make r-bridge-selftest          # S06 regression, 32 checks
make upstream-manifest-verify   # 53 files, 1 locally modified (NAMESPACE)
```

Evidence: `validation/exact/evidence/s07_output_transport.txt`.
