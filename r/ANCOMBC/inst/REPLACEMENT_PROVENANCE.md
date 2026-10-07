# Replacement build provenance

This directory is a **replacement build** of ANCOMBC. It is not the upstream
tarball and does not claim to be. This file says exactly what was carried over,
what was changed, and what the change is for, so that a reader never has to
guess which is which.

## Upstream

| item | value |
| --- | --- |
| package | ANCOMBC |
| version | 2.15.2 |
| commit | `dc4febdf59badb3a8dfe0c767ef2186323c2199a` |
| commit date | 2026-09-18 12:18:29 -0400 |
| commit subject | `Add the ANCOMBC-concepts help topic and bump version to 2.15.2` |
| repository | <https://github.com/FrederickHuangLin/ANCOMBC> |
| licence | Artistic-2.0 |
| declared R requirement | `R (>= 4.5.0)` |
| retained upstream source | `reference/ANCOMBC/`, an immutable oracle checkout |

## What is carried over, byte for byte

`R/`, `NAMESPACE`, `NEWS`, `README.md`, `data/`, `man/`, `inst/CITATION`,
`inst/REFERENCES.bib`, `tests/` and `vignettes/`, copied from the pinned commit
with no edits. `inst/UPSTREAM_MANIFEST.sha256` lists every file in this package
with the SHA-256 of the copy here and of the pinned checkout, and marks which
files upstream's own digest list covers.

Attribution is unchanged: `Authors@R`, `License` and `inst/CITATION` are
upstream's. This build adds no licence terms and removes none.

## What this build adds

| file | why |
| --- | --- |
| `DESCRIPTION` | the upstream fields verbatim, plus `X-Replacement-Of`, `X-Build-Provenance` and `X-Upstream-License`, so the installed package does not misrepresent its build provenance. |
| `inst/REPLACEMENT_PROVENANCE.md` | this file. Its presence is what lets `scripts/exact_runner.R` refuse to run the original under the candidate's name. |
| `inst/UPSTREAM_MANIFEST.sha256` | the per-file provenance record, with a verification script. |
| `NAMESPACE` | upstream's, plus one `useDynLib(ANCOMBC, .registration = TRUE)` line. The only change to a retained upstream file, and the reason is in the next section. |
| `R/bridge.R` | the R side of the typed native bridge: the input transport (S06) and the output transport (S07). Reached from nowhere in the retained upstream code. |
| `R/assemble.R` | rebuilds the original's `.ancombc2_core()` result from transported payloads. Reached from nowhere in the retained upstream code. |
| `src/init.c`, `src/Makevars` | the native bridge: R type checks, buffer ownership, `.Call` registration, and the output emit shim. |
| `tests/bridge_selftest.R` | the input transport's acceptance checks (IMPROVED_PLAN.md S06). |
| `tests/output_selftest.R` | the output transport's acceptance checks (IMPROVED_PLAN.md S07). |
| `tests/bridge-check.R` | runs both selftests under `R CMD check`. |

## Local modifications to retained upstream files

`scripts/upstream_manifest.py --verify` lists every file in this package whose
digest differs from `reference/ANCOMBC/`. There is exactly one:

| file | change | why |
| --- | --- | --- |
| `NAMESPACE` | one `useDynLib(ANCOMBC, .registration = TRUE)` directive added, above `import(stats, except = filter)` | `.registration = TRUE` is what makes `R_init_ANCOMBC` in `src/init.c` run, and it is what makes a renamed Rust entry point a **load-time** error instead of a run-time one. It adds no import and no export, so the package's exported names, their signatures and their defaults are untouched. Nothing else in the retained R code changed: `R/ancombc2.R` is still byte-identical to the pinned commit, which is the point -- the bridge is reached from `R/bridge.R`, not by editing `ancombc2()`. |

Verified by `python3 scripts/upstream_manifest.py --verify`: `45 files, 0 locally
modified` before the bridge existed, and after S06 and S07, still exactly the one
`NAMESPACE` row above, which the same script refuses to accept unless this file
names it. The count grew to 56 files as `R/assemble.R` and
`tests/output_selftest.R` were added; both are new files with no upstream
counterpart, so neither appears as a modification.

## What is deliberately not here yet

This is the **compatibility scaffold**, not the accelerated implementation.
As of this commit:

* the entire `ancombc2()` interface runs from retained upstream code, including
  random effects (`rand_formula`), the global test, pairwise comparisons,
  Dunnett's test, the trend test, both sensitivity modes and taxonomic
  aggregation;
* no Rust code is called from `ancombc2()`;
* `r/ancombc2rs/` -- the pre-existing JSON-bridge prototype -- is **not** part of
  this package and is not on its search path.

Compatibility here means the public API, the exported names, the returned
schema and the results are the original's. It does not mean this build is
upstream's build, and it does not by itself satisfy any performance gate.

## Why the scaffold is the right first step

The numerical stages that will move to Rust are guarded by `lme4`, `MASS::ginv`,
`nloptr`, `multcomp`, `quadprog`, R's formula semantics and R's own QR and
solve. Re-deriving those from a function name is how a "compatible"
reimplementation quietly diverges. Retaining the original code, and proving
byte-identical results against the installed original first, means every later
Rust change is measured against a baseline that is known to be exact rather than
against a baseline that was already approximate.

## How to verify this file

```sh
python3 scripts/upstream_manifest.py --verify
make exact-env-check PROFILE=linux-r453-openblas
```

The first proves every retained file is still byte-identical to the pinned
checkout and that no local modification is undocumented. The second proves the
runtime, the pinned digests and the installed original.
