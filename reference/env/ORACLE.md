# Oracle freeze for ancombc2-rs

## Compatibility target

```
ancombc2-rs v0.1  ==  ANCOMBC 2.15.2 @ dc4febdf59badb3a8dfe0c767ef2186323c2199a
```

This string is asserted in `crates/ancombc2-core/src/compat.rs` (`ORACLE_SHA`) and in
`docs/numerical_contract.md`. Any change to it is a breaking change to the
compatibility claim and requires an explicit version bump in `CHANGELOG.md`.

## Oracle

| Item | Value |
|---|---|
| Package | ANCOMBC |
| Version | 2.15.2 |
| Repository | https://github.com/FrederickHuangLin/ANCOMBC |
| Commit | `dc4febdf59badb3a8dfe0c767ef2186323c2199a` |
| Commit date | 2026-09-18 12:18:29 -0400 |
| Commit subject | `Add the ANCOMBC-concepts help topic and bump version to 2.15.2` |
| License | Artistic-2.0 |
| Declared R requirement | `R (>= 4.5.0)` |

## Execution environment actually used to generate goldens

| Item | Value |
|---|---|
| R version | 4.3.3 (2024-02-29) "Angel Food Cake" |
| Platform | x86_64-pc-linux-gnu (64-bit) |
| Extra library path | `/scratch/mdra00001/rlib` |
| RNG | `set.seed(42)`, default RNGkind |
| Rust | 1.98.1 (797e8a9bc 2026-08-05) |

The oracle's DESCRIPTION declares `R (>= 4.5.0)`. The golden fixtures in this
repository were generated on **R 4.3.3**, one minor series below the declared
minimum. The fixed-effects code path exercised here (`R/ancombc_prep.R`,
`R/ancombc_bias_correct.R`, `R/ancombc_mult.R`, `R/utils.R`) uses no R-4.4+-only
feature, and every quantity is a `stats`/base-R computation whose semantics
(`p.adjust`, `quantile` type 7, `pt`, `pf`, `pchisq`, `lm.fit`, `MASS::ginv`,
`nloptr` NLOPT_LN_NELDERMEAD) are stable across 4.3/4.5. `sessioninfo.txt` in
this directory records the full session.

**Re-run the golden generator on R >= 4.5.0 before publication and diff the
result.** The `nightly` CI job does this when a >= 4.5 interpreter is available.

## Packages the fixed-effects path needs

Only these are loaded. Every other ANCOMBC import is reached exclusively from
`rand_formula`, `dunnet`, or `trend` code, all of which are out of scope for v1.0.

| Package | Version | Used for |
|---|---|---|
| stats | base | `p.adjust`, `quantile`, `pt`, `pf`, `pchisq`, `dnorm`, `lm.fit`, `model.matrix`, `lm`, `anova`, `var`, `quantile` |
| MASS | 7.3-x (shipped) | `ginv` in `.sandwich_vcov` and `.ancombc_global_F` |
| nloptr | 2.x | `nloptr` / `NLOPT_LN_NELDERMEAD` in `.bias_em` |

Deliberately **not** loaded: `lme4`, `lmerTest`, `multcomp`, `quadprog`,
`doParallel`, `doRNG`, `foreach`, `DescTools`, `Hmisc`, `Rdpack`, `energy`,
`gtools`, `Matrix`. `reference/R/oracle.R` stubs `foreach`/`%dorng%`/
`registerDoSEQ` sequentially, which is exactly the `n_cl = 1` behaviour of the
real package.

## How the oracle is loaded

`reference/R/oracle.R` **sources** the pinned R files rather than installing the
package. Reasons:

1. `R (>= 4.5.0)` is declared but only 4.3.3 is available.
2. Sourcing gives the harness access to `.iter_mle`, `.lm_fit_all`,
   `.sandwich_vcov` and `.bias_em` — all `.`-prefixed internal functions that are
   not exported — which is required to capture golden quantities 4–15 of the
   contract in `docs/numerical_contract.md`.
3. Instrumenting a copy is the only way to obtain `beta*`, `theta` traces,
   residuals, and the sandwich block, none of which the package returns.

The source is verified against the pinned commit before loading:

```r
stopifnot(identical(
  system("git -C reference/ANCOMBC rev-parse HEAD", intern = TRUE),
  "dc4febdf59badb3a8dfe0c767ef2186323c2199a"))
```

The vendored copy lives at `reference/ANCOMBC/` (see `reference/README.md`).

**Where the oracle is fetched from.** The canonical repository is
`https://github.com/FrederickHuangLin/ANCOMBC`. Both CI jobs that need the oracle
(`oracle` and `benchmarks`) clone that URL and then run the `test` above, so a
fetch that lands anywhere else fails the job instead of silently benchmarking a
different package. The `bench-container` job applies the same check inside the
image build, so an image with the wrong commit cannot be built.

This was not always so: both jobs pointed at `https://github.com/zdk123/ANCOMBC`,
a different repository. The vendored copy, and this document, are from
`FrederickHuangLin/ANCOMBC` at `dc4febdf...`.
