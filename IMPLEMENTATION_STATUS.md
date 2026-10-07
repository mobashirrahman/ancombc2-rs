# Implementation status

Authoritative execution record for `IMPROVED_PLAN.md` (tasks S00-S21, gates
C1-C4/B1/P1-P5). `PLAN.md`, `docs/performance_plan.md` and every file listed in
the audit manifest are **history**: they are preserved, not corrected in place,
and their numbers are not the release baseline.

## Objective and baseline

- Plan: `IMPROVED_PLAN.md` (dated 2026-10-05, status "implementation has not started")
- Starting commit: `e8ec02e16e7b7cf0cb0d62ef3efbf93bc135c59b`
- Pre-existing dirty files (preserved, never reset):
  - `M PLAN.md`
  - `M scripts/bench_gates.py` (user change: `n/a`-safe printing of
    `r_1core_wall_s` when a re-run shadows an earlier good row)
  - `M scripts/__pycache__/bench_gates.cpython-312.pyc` (tracked cache of the above)
  - `?? IMPLEMENTATION_PROMPT.md`, `?? IMPROVED_PLAN.md` (untracked at session start)
  - `git diff --binary | sha256sum` = `fc84c15d71c2597c8961fd59b549ab5f0ae6e88eb48e32de567202590801343a`
- Certified runtime profiles: **none certified yet** (S01 in progress). Two
  candidate profiles exist on this host; see the table below.
- Original library / package hash: the pinned source is vendored at
  `reference/ANCOMBC/` and its per-file SHA-256 list is
  `reference/env/oracle.sha256` (22 files). No *installed* original exists yet.
- Candidate library / artifact hash: **none**. `r/ancombc2rs` is the pre-existing
  scaffold (package name `ancombc2rs`, JSON FFI); `r/ANCOMBC/` does not exist.
- Active task and exact next substep: **S04** -> build the rest of the exact
  fixtures under `validation/exact/` (fx04 plus one input per case the small set
  will grow into), each from its committed file read once through
  `scripts/make_exact_input.R`, and add a regeneration check that reproduces the
  inputs and both arms' outputs byte for byte in a temporary directory. Then
  wire the observed-case list into `run_exact.py --json-out` so
  `make exact-cases` consumes real observations rather than a hand-written
  file. S00-S03 are `passed`.

### Host and runtime inventory (measured 2026-10-05, session 1)

| item | value | how measured |
| --- | --- | --- |
| kernel / arch | Linux 5.x, `x86_64` | `uname -a` |
| CPU | AMD Ryzen 7 3700X, **8 physical cores**, 16 logical (2 threads/core), 1 socket | `lscpu` |
| RAM | 32,805,716 kB total | `/proc/meminfo` |
| container engine | **absent** (`docker`, `podman` not on PATH) | `which` |
| cgroup v2 | mounted, `memory` controller delegated at `/sys/fs/cgroup`, **no write permission** (`mkdir` denied at root, `user.slice`, and `user-7059086.slice`) | `mkdir -p` probes |
| sudo | password required, not available | `sudo -n true` |
| network | HTTPS to `cloud.r-project.org`, `micro.mamba.pm` reachable | `curl -sSI` |
| system R | **4.3.3** "Angel Food Cake", reference BLAS `/usr/lib/x86_64-linux-gnu/blas/libblas.so.3.12.0` | `Rscript -e 'extSoftVersion()'` |
| user R library | `$ANCOMBC_RLIB` (phyloseq, nloptr, foreach, ...; **no** lme4/lmerTest/multcomp/quadprog/doParallel/doRNG/energy/gtools/DescTools/Hmisc/Rdpack) | `installed.packages()` |
| cargo / rustc | present on PATH (`/home/mdra00001/.cargo/bin`) | `which` |

**Provisioned during session 1 (S01 groundwork):** `micromamba 2.9.0` at
`$MM/bin/micromamba`, and the committed explicit
lock `benchmarks/container/environment-linux-64.lock` realised at
`$MM/benchenv`:

| item | value |
| --- | --- |
| R | **4.5.3** (2026-03-11) "Reassessed Reassurer", `x86_64-conda-linux-gnu` |
| BLAS | `$MM/benchenv/lib/libopenblasp-r0.3.34.so` (OpenBLAS 0.3.34, pthreads) |
| LAPACK | `liblapack.so.3 -> libopenblasp-r0.3.34.so` (same OpenBLAS) |
| all 18 `Imports` of the pinned DESCRIPTION | present, versions recorded in S01 evidence |
| invoked as | `micromamba run -p $MM/benchenv Rscript ...` |

Candidate profiles, both on this host, neither certified until S01's verifier
runs:

| profile id | R | BLAS | role |
| --- | --- | --- | --- |
| `linux-r453-openblas` | 4.5.3 | OpenBLAS 0.3.34 | **primary**: matches the committed lock and the pinned `Depends: R (>= 4.5.0)` |
| `linux-r433-refblas` | 4.3.3 | reference netlib libblas 3.12.0 | secondary: the runtime the legacy goldens were drawn on; below the pinned floor, so it cannot certify the *package*'s `Depends` |

## Tasks

| ID | State | Changed files | Acceptance evidence | Remaining work |
| --- | --- | --- | --- | --- |
| S00 | passed | `IMPLEMENTATION_STATUS.md`, `validation/exact/audit/s00_claims_audit.json`, `README.md`, `docs/compatibility.md`, `docs/release_status.md` | `git status --porcelain` still shows only the three pre-existing modifications plus this session's additions; the three pre-existing files were not touched (`git diff --stat` unchanged for them). `python3 -c json.load(validation/exact/audit/s00_claims_audit.json)` -> 12 claims. Banners added to `README.md`, `docs/compatibility.md`, `docs/release_status.md` naming `IMPROVED_PLAN.md`/`IMPLEMENTATION_STATUS.md` and stating exact replacement + corrected performance are unfinished | none |
| S01 | passed | `scripts/probe_runtime.R`, `scripts/verify_profile.py`, `Makefile` (`exact-env-check`, `exact-env-selftest`), `validation/exact/profiles/linux-r453-openblas.json`, `validation/exact/profiles/linux-r433-refblas.json`, `.rlib/original/` (installed original) | `make exact-env-check PROFILE=linux-r453-openblas` -> **22/22 pass, exit 0**. `make exact-env-selftest` -> **11/11 mutations non-zero, each tripping its designed check; unmutated profile 22/22, exit 0**. `make exact-env-check PROFILE=linux-r433-refblas` -> **exit 1** (`runtime.satisfies_pinned_depends`, `original.installed`), correctly identifying it as uncertifiable before analysis. Log: `validation/exact/evidence/s01_profile_preflight.txt`. Original public API dumped to `validation/exact/evidence/s01_original_public_api.txt` | none for the primary profile. `certified` stays `false` in the JSON until S14 has run every required case against this installed original |
| S04b | passed | `scripts/make_case_manifest.py`, `validation/exact/cases.json`, `Makefile` (`exact-case-manifest`, `exact-case-manifest-check`) | `validation/exact/cases.json` is **generated** from `validation/exact/fixtures.json`; `make exact-case-manifest-check` -> "up to date (25 small, 114 required)", exit 0. One source of truth, so the manifest and the suite cannot describe different case sets. `run_exact.py --observed-out` writes the observed list in the shape `check_exact.R --cases` reads, so `make exact-cases-required` consumes real observations | the 114-case `required` set is the whole fixture manifest; S14 widens it with the matrix cells, the real datasets, 100+ generated seeds and a small valid/invalid case per public parameter |
| S02 | passed | `scripts/exact_runner.R`, `scripts/run_exact.py`, `scripts/with_profile_r.py`, `scripts/make_exact_input.R`, `Makefile` (`exact-input`, `exact-run`, `exact-runner-selftest`), `validation/exact/inputs/{fx01,fx02,fx03}-basic.rds` | Log: `validation/exact/evidence/s02_exact_runner.txt`. Two runs of the **installed original** per fixture, `pkg_identity=fd42b62f1464f428` matching the S01 pin every time: fx02 58245 B `79c7da78139554f2` twice; fx03 1848544 B `0175f5d238f4bae1` twice; fx01 error payload 744 B `25b8aae4edb88c04` twice. `make exact-runner-selftest` -> **8/8 guard checks non-zero**, including `original_refused_as_candidate` and `replacement_refused_as_original` | candidate arm: none to compare until S05; `run_exact.py` reports `candidate: not-installed-yet` rather than a pass |
| S03 | passed | `scripts/check_exact.R`, `validation/exact/cases.json`, `Makefile` (`exact-comparator-selftest`, `exact-compare`, `exact-cases`) | Log: `validation/exact/evidence/s03_comparator.txt`. `make exact-comparator-selftest` -> **30/30 pass, exit 0**: identical objects accepted; one flipped mantissa bit, `-0`<->`+0` both ways, `NA_real_`<->`NaN` both ways, `Inf`->finite, `42L`->`42`, a changed name, data.frame and matrix row order, `NULL`->`list()`, a dropped NULL field, swapped list elements, a changed and an added attribute, a factor value swap, an extra factor level, a changed warning, a changed `.Random.seed`, a reordering of conditions, and a 1e-12 relative difference all **rejected**. `make exact-compare` on two original runs -> `VERDICT: identical`, exit 0; on a capture with byte 5000 of `result.bytes` flipped -> exit nonzero, naming `$$feature_table` index 1241 (1 of 3000). Manifest: complete set -> exit 0; missing / duplicated / absent-observed / unknown set -> nonzero in all four | the observed-case file is hand-built here; S04 wires `run_exact.py` to emit it |
| S04 | passed | `validation/exact/fixtures.json` (95 cases), `scripts/exact_input_lib.R`, `scripts/make_exact_inputs.R`, `scripts/make_exact_source_fixtures.R`, `scripts/check_exact_inputs.py`, `scripts/make_exact_input.R`, `validation/exact/inputs/` (95 RDS + `built.json` + `manifest.json`), `validation/exact/fixtures-src/` (4 new fixtures), `Makefile` (`exact-inputs`, `exact-inputs-verify`, `exact-source-fixtures`) | `python3 scripts/check_exact_inputs.py --verify` -> **"REGENERATION: byte-identical for every input and every source file", exit 0** over 95 cases: it rebuilds every RDS in a temporary directory and compares SHA-256 of the input *and* of each committed `counts.tsv`/`meta.tsv`/config it was read from, and asserts the manifest's `formal_order` against the installed `formals(ANCOMBC::ancombc2)`. `scripts/make_exact_source_fixtures.R` writes each new fixture, reads it back through the same reader the exact inputs use, and compares with `tolerance 0` before writing the manifest | the expectation table was **wrong on 19 of 72** first-draft cases; corrected against measured reference behaviour and re-measured. See the S04 findings below |
| S05 | passed | `r/ANCOMBC/` (DESCRIPTION, NAMESPACE, R/, data/, man/, inst/, tests/, vignettes/, NEWS, README.md), `r/ANCOMBC/inst/REPLACEMENT_PROVENANCE.md`, `r/ANCOMBC/inst/UPSTREAM_MANIFEST.sha256`, `scripts/upstream_manifest.py`, `scripts/check_replacement.R`, `Makefile` (`replacement-build`, `replacement-check`, `upstream-manifest*`), profile `replacement` block | `python3 scripts/upstream_manifest.py --verify` -> **45 files, 0 locally modified, "UPSTREAM MANIFEST: consistent"**; it also re-checks that every `reference/env/oracle.sha256` digest still holds. `R CMD INSTALL --library=.rlib/replacement r/ANCOMBC` -> `* DONE (ANCOMBC)`. `make replacement-check` -> **24/24 pass, exit 0** (log below): identical 32-formal names and order, byte-identical `args()` text, per-formal default expressions identical, `pseudo_sens` TRUE, no `...`, the original's 7 exports present and none added, every export's signature identical, attribution and licence retained, `X-Replacement-Of` present, named and positional calls succeed, and **no namespace loaded from the original's library and the original's library absent from the candidate's search path** | the whole `ancombc2()` interface currently runs from retained upstream R; **no Rust is called**. This is the compatibility scaffold, exactly as `inst/REPLACEMENT_PROVENANCE.md` states. Acceleration is S08-S13 (S06-S07 moved the data, not the arithmetic) |
| S04a | passed | `validation/exact/fixtures.json` corrections after measurement, `scripts/exact_input_lib.R` (`parse_override` `__OMITTED__`/`__NULL__`, `as_r_literal` tagged object, `parse_count_cell`, `args_from` modes), `scripts/make_case_manifest.py` | **final campaign: 114/114 byte-identical, 0 outcome mismatches, 0 reference nondeterminism, 0 infrastructure failures, 0 declared-vs-measured disagreements.** 82 success / 32 error (reference errors are comparable outcomes, each with its own serialized condition compared byte for byte). `make exact-cases-required` -> "MANIFEST: complete", exit 0. The 25-case `small` set runs and checks separately -> "MANIFEST: complete", exit 0. Evidence: `validation/exact/evidence/s04_two_arm_summary.txt`, `s04_two_arm_114.txt`, `s04_two_arm_small_25.txt`, `s04_exact_inputs.txt`; machine-readable `validation/exact/two_arm.json`, `validation/exact/observed.json`, `validation/exact/inputs/manifest.json` | none for this pass. S14 still has to add the matrix cells, the four real datasets, 100+ generated seeds, a small valid/invalid case per public parameter, and thread-comparison cases |
| S06 | passed | `docs/exact_seam_inventory.md`, `crates/ancombc2-rbridge/` (Cargo.toml, `src/{lib,transport,owned,abi}.rs`), `r/ANCOMBC/src/{init.c,Makevars}`, `r/ANCOMBC/R/bridge.R`, `r/ANCOMBC/tests/{bridge_selftest.R,bridge-check.R}`, `r/ANCOMBC/NAMESPACE` (+1 line), `Makefile` (`r-bridge-selftest`, `r-bridge-rust-tests`, `r-core-without-r`), profile `replacement.native` | Log: `validation/exact/evidence/s06_native_transport.txt`. **`make r-bridge-selftest` -> 32/32 pass, exit 0** against the *installed* replacement: every one of the 64 bits of every cell survives; `NA_integer_` returns as `INT_MIN`; `NA_real_` returns with R's own NaN payload `0x7ff00000000007a2` and a computed `NaN` does **not** become it; `-0.0` and `+0.0` stay distinguishable where `identical()` and `sprintf()` cannot tell them; `±Inf` survive; dims, four name vectors, the `complete.cases` mask and the 1-based group index (with `0` = absent) all survive; the returned buffer does not alias the input; 11 malformed inputs each produce an R error naming the argument. `make r-bridge-rust-tests` -> 14/14 + clippy `-D warnings` clean. `make r-core-without-r` -> **154 core tests pass with no R present**. `python3 scripts/upstream_manifest.py --verify` -> 54 files, **1 locally modified** (`NAMESPACE`, one `useDynLib` line, documented in `inst/REPLACEMENT_PROVENANCE.md`). The 25-case small set still **25/25 byte-identical** with the bridge linked | `ancombc2()` does not call the bridge yet: a loaded-but-uncalled bridge is not Rust execution and is not reported as such. The numerical seam is S08-S13 |
| S07 | passed | `docs/exact_output_inventory.md`, `scripts/capture_core_payloads.R`, `validation/exact/payloads/tiny-defaults.rds`, `crates/ancombc2-rbridge/src/output.rs`, `crates/ancombc2-rbridge/src/{lib,abi,transport}.rs`, `r/ANCOMBC/src/init.c` (`wire_doubles`, `C_ancombc2_rb_emit`), `r/ANCOMBC/R/bridge.R` (`rb_emit_payload`, `rb_restore_names`, `rb_flatten_vcov`), `r/ANCOMBC/R/assemble.R`, `r/ANCOMBC/tests/output_selftest.R`, `Makefile` (`r-output-selftest`), `r/ANCOMBC/inst/REPLACEMENT_PROVENANCE.md` | Log: `validation/exact/evidence/s07_output_transport.txt`. **`make r-output-selftest` -> 43/43 pass, exit 0** against the *installed* replacement: the eleven payloads come back `serialize()`-identical to the pinned original's, `dof` included as `INTSXP`; `vcov_hat` stays a list of `p x p` matrices with `dimnames`, and a deliberately non-symmetric `p = 3` block returns `(1,2)`/`(2,1)`/`(1,3)` to their original places (symmetry cannot catch a transpose, so it is not used as the check); `NA_real_`, `NaN`, `Inf`, `-Inf` and `-0.0` all survive, the sign of zero checked with `1/x` because `-0.0 == 0.0`; a reported taxon set shorter than the bias set round-trips and a shorter bias set is refused; the **assembled `.ancombc2_core()` result is byte-identical to the original's**, all ten names and four `NULL`s. `make r-bridge-rust-tests` -> **25/25** + clippy `-D warnings` clean and `cargo fmt --check` clean, including `the_wire_layout_is_pinned`, which pins all 28 offsets of the 216-byte `RawPayloads` that `src/init.c` mirrors, and `vcov_is_the_only_row_major_entry`. `make r-core-without-r` -> 154 core tests still pass. `make r-bridge-selftest` -> 32/32 (S06 regression). `python3 scripts/upstream_manifest.py --verify` -> **53 files, 1 locally modified** (`NAMESPACE` only), stable across rebuilds because compiler output is excluded. `make exact` re-run on the final build: **114/114** byte-identical (82 results + 32 errors, 0 problems), `make exact-cases-required` -> MANIFEST complete; `make exact-small` **25/25** | **This is transport parity, not Rust numerical parity.** The bridge computes nothing: the numbers came from the pinned original, so agreeing with them proves the path loses nothing and nothing more. C3 stays open and `ancombc2()` still does not call the bridge. `global`/`pairwise`/`dunnet`/`trend` are not transported and asking for one is a named error, not a `NULL`; their assemblies need `group` and `para2$fits`, which S08/S09 supply. The single fixture has all four switches `FALSE`, so `struc_zero`'s effect on `taxa`/`taxa_bias` is covered by a synthetic case rather than by a captured one. |
| S08 | passed | `docs/exact_preprocess_stages.md`, `crates/ancombc2-core/src/reduce.rs` (new), `crates/ancombc2-core/src/preprocess.rs`, `crates/ancombc2-core/src/workspace.rs`, `crates/ancombc2-core/src/lib.rs`, `crates/ancombc2-rbridge/src/output.rs` (`preprocess_stages`, `PreprocessStages`, `PreprocessInputs`, `RBackedReductions`, `Reducer`), `crates/ancombc2-rbridge/src/abi.rs` (`ancombc2_rb_preprocess_probe`, `RawStage`, `RawStages`, `RawPreprocess`), `r/ANCOMBC/src/init.c` (`eval_na_rm_reduction`, `reduce_rows`, `C_ancombc2_rb_preprocess_probe`), `r/ANCOMBC/R/bridge.R` (`rb_preprocess_stages`, the four reduction wrappers), `scripts/check_preprocess_stages.R` (new), `Makefile` (`preprocess-stages`), `r/ANCOMBC/tests/bridge-check.R` | Log: `validation/exact/evidence/s08_preprocess_stages.txt`. **`make preprocess-stages` -> 274/274 stage comparisons pass** across 13 fixtures (tiny, sparse, missing-value double and integer, boundary-prevalence, boundary-library, boundary-zero-pseudo, one-sample, one-group, no-group, aggregate-differs, neg-lb, neg-lb-off), each compared with `identical(serialize(x, NULL, ascii = FALSE, xdr = TRUE, version = 3), ...)` against the reference's own expressions transcribed into the script rather than imported from the installed original. All twenty stages exact: `prevalence1`, `tax_keep1`, `lib_size1`, `samp_keep1`, `O1`/`log1`/`means1`/`y1`, the aggregate eight, `zero_keep`, `group_prevalence`, `group_size`, `group_lower`, `zero_ind`. `make exact` on the final build -> **114/114** byte-identical (82 results + 32 errors, 0 problems) and `make exact-cases-required` -> MANIFEST complete; `make exact-small` -> **25/25** -- both with the four library digests identical before and after, which is the check that no install happened mid-campaign. `ancombc2_run_named` and `core_run_at_pseudo` now take the reduction as a **parameter, never a default**, so every call site says which accumulator it used; `parallel::map_par` passes the item **index** as well as the item, because the pseudo-count grids repeat entries and a mean looked up by value would be ambiguous. The reduction is also unreachable from the parallel sensitivity refits by construction -- `dyn Reductions` is not `Sync` and R cannot be called off the main thread -- so the centring means are reduced on the calling thread (`grid.len() * taxa.len()` doubles) and the workers redo only the `ln`; `reduce::split_tests` asserts that split is **bit-identical** to the unsplit route, using a reduction that deliberately is not `f64` so it cannot pass by coincidence. `cargo clippy --release --workspace --all-targets -- -D warnings` and `cargo fmt --check` clean; whole workspace green (core 170, rbridge 44, integration 88); `make r-bridge-selftest` -> 32/32, `make r-output-selftest` -> 43/43, `make replacement-check` -> 24/24; `upstream-manifest-verify` -> 53 files, **1 locally modified** (`NAMESPACE` only) | **Ten measured divergences, none of which produced a wrong number** -- every one is invisible to `all.equal()` and visible only to `serialize()`. `f64::ln` does not preserve a NaN payload (`0x7ff00000000007a2` -> `0x7ff8000000000002`); R's `is.infinite` mask is not "not finite"; `NA_real_` is a *signalling* NaN that arithmetic **quiets** (bit 51, not the sign) while the hardware clears the payload entirely; unary math functions propagate a NaN unchanged and binary operators quiet it, so they need different rules; `rowSums`/`rowMeans` accumulate in `long double` (64-bit mantissa -- `rowMeans(c(1e16,1,-1e16,1e-17,1))` is `0x1.999999999999ap-2` against an `f64` loop's `0x1.999999999999ap-3`, a factor of two), so the reduction is asked of R rather than re-derived, which is right on every platform by construction; `lib_size` is over the retained taxa and *all* samples; the structural-zero screen reads the **aggregate** table; `.data_core`'s second pass reports subset-relative indices (both forms now reported); `colSums` is named and `rowMeans` is not. Two boundary errors also: the probe first returned **pointers into Rust memory** that dangled on return -- no crash, no NULL, right shapes, first two elements of each stage garbage -- and then never wrote `out_f64_len`, so C allocated one element and Rust wrote a hundred past it, surfacing later as `"recursive gc invocation"`. Buffers are now caller-owned with offsets, and both lengths are reported before the copy | **Stage exactness, not Rust numerical parity.** No fit is computed here: the reductions are R's and the arithmetic is checked step by step. `theta_hat` is absent because it needs `beta1` (S09); the design matrix is an S06 input, not a recomputed stage. Thirteen hand-built boundary fixtures, deliberately small enough that a difference is readable in a hex dump; the 114-case campaign is re-run on the final build and confirms nothing else moved. **The pipeline still uses `F64Reductions`** -- `ancombc2()` does not call it yet -- so an `f64` row mean can still reach a *CLI* result; S09 makes the R-backed reduction the one the R session passes. C3 stays open |
| S09 | **partial** (least squares done; iterative MLE not yet measured) | `crates/ancombc2-core/src/matrix/dqrls.rs` (new), `crates/ancombc2-core/src/mle.rs`, `crates/ancombc2-core/src/reduce.rs` (`long_double_reduce`, `mean_na_rm`), `crates/ancombc2-core/src/correct.rs`, `crates/ancombc2-rbridge/src/abi.rs` (`ancombc2_rb_dqrls`, `ancombc2_rb_set_blas`), `r/ANCOMBC/src/init.c`, `r/ANCOMBC/R/bridge.R` (`rb_blas_kind`), `crates/ancombc2-core/tests/dqrls_vs_r.rs`, `scripts/make_dqrls_cases.R`, `scripts/check_fit_stages.R`, `Makefile` (`dqrls-vs-r`, `fit-stages`) | Logs: `validation/exact/evidence/s09_fit_stages.txt`, `s09_dqrls_vs_r_openblas.txt`. **`make fit-stages` -> 191/191 pass** (was 37 pass / 53 fail): every fixture's `beta`, `fitted` and `dof` are `serialize()`-identical to `.lm_fit_all` transcribed from the pinned source, including `NA_real_` payloads, an incomplete design row, a taxon with no usable sample, and a factor level absent for a taxon (the per-taxon `lm()` fallback, with `drop.unused.levels`). **`make dqrls-vs-r` -> 440/440 problems bit-identical** to R's own `Cdqrls` under OpenBLAS 0.3.34 (R 4.5.3), and 440/440 under the netlib reference BLAS (R 4.3.3), sizes 3..20000, aliased/near-collinear/badly scaled/zero columns included. Whole workspace green, clippy `-D warnings` and `fmt` clean; `make r-bridge-selftest` 32/32, `r-output-selftest` 43/43, `replacement-check` 24/24, `preprocess-stages` pass, `r-bridge-rust-tests` pass | `.iter_mle`'s loop (theta update, convergence, the rank smoke test) is not yet compared stage by stage with R's; the golden `fx04` deviation (4e-9, rank-deficient taxa) is unchanged and unexplained by the fit. `ancombc2()` still does not call Rust. **Limits that are properties of R, not of this port:** (1) the BLAS rounding is part of the result, so byte identity is defined *per BLAS family* -- netlib and OpenBLAS Haswell/Zen are implemented and calibrated at run time against R's own `Cdqrls` (`rb_blas_kind()` is `NA` when neither matches, and then nothing may claim byte identity); other CPUs' OpenBLAS kernels (SkylakeX/AVX-512, etc.) are not implemented; (2) OpenBLAS threads `ddot`/`daxpy` above n = 10000, so R's own answer there depends on its thread count -- identical only with `OPENBLAS_NUM_THREADS=1`, and the profile pins that |
| S10 | pending | | | |
| S11 | pending | | | |
| S12 | pending | | | |
| S13 | pending | | | |
| S14 | pending | | | |
| S15 | pending | | | |
| S16 | pending | | | |
| S17 | pending | | | |
| S18 | pending | | | |
| S19 | pending | | | |
| S20 | pending | | | |
| S21 | pending | | | |

## Last checkpoint

- What is implemented: **S00-S08 are all `passed`.** The headline measurement is
  still the scaffold one: **114 of 114 exact cases return byte-identical complete
  results from the installed original and from the installed replacement**, with
  zero outcome mismatches, zero reference nondeterminism and zero declared-vs-
  measured disagreements. Both arms still run the retained upstream R, so that
  proves the harness, the inputs, the comparator and the schema. It is not Rust
  parity and it does not pass C3.
- What S06, S07 and S08 add on top of that is the seam itself, both directions and
  then the arithmetic in between. S06 (`docs/exact_seam_inventory.md`) put the core's
  inputs on a typed native wire: 32/32 checks, every one of the 64 bits of every cell
  surviving, `NA_integer_` as `INT_MIN`, R's own NaN payload kept distinct from a
  computed `NaN`, `-0.0` still distinguishable where `identical()` cannot tell them.
  S07 (`docs/exact_output_inventory.md`) did the same for the outputs: 43/43 checks,
  all eleven payloads `serialize()`-identical, `dof` back as `INTSXP`, `vcov_hat` back
  as a *list*, and the assembled `.ancombc2_core()` result byte-identical.
- S08 (`docs/exact_preprocess_stages.md`) is where the arithmetic starts, and it found
  **ten measured divergences, none of which produced a wrong number**: every one is
  invisible to `all.equal()` and visible only to `serialize()`. The largest is that
  R accumulates `rowSums`/`rowMeans` in `long double`, so a cancelling sum keeps bits
  an `f64` loop drops -- a factor of two on `rowMeans(c(1e16, 1, -1e16, 1e-17, 1))`.
  The reduction is now asked of R rather than re-derived, which is correct on every
  platform by construction. The rest are `ln` not preserving NaN payloads, R quieting a
  `NA` where the hardware clears it, `is.infinite` not meaning "not finite", and four
  reference semantics read off the code rather than assumed. **274/274 stage
  comparisons now pass.**
- **None of this is Rust numerical parity.** S06 and S07 move data; S08 checks
  preprocessing, with R supplying the reductions. The first stage that *computes* is
  S09. C3 stays open.
- One consequence worth flagging for S09: **the reduction cannot be reached from a
  worker thread.** The sensitivity analysis refits in parallel, so its centring means
  are now reduced on the calling thread and the workers redo only the arithmetic --
  a split that `reduce::split_tests` shows is bit-identical. Any future primitive that
  needs R has to be hoisted the same way, and cannot be hidden behind a default.
- Earlier (superseded): **S00-S03 are done.** Two execution profiles are
  pinned and machine-verified; the unmodified pinned original is installed in an
  isolated library under the primary profile. A real two-arm exact runner exists
  and has established that the installed original is byte-repeatable on three
  fixtures. Every claim this repository makes about byte identity or
  R-relative performance is now marked unvalidated, with the reason, in
  `validation/exact/audit/s00_claims_audit.json`.
- Exact commands run and outcomes:
  - `git rev-parse HEAD` -> `e8ec02e16e7b7cf0cb0d62ef3efbf93bc135c59b`
  - `git diff --binary | sha256sum` -> `fc84c15d...`
  - `micromamba create -y -p .../benchenv --file benchmarks/container/environment-linux-64.lock` -> success (779 log lines, `Transaction finished`)
  - `micromamba run -p .../benchenv R --version` -> `R version 4.5.3 (2026-03-11) -- "Reassured Reassurer"`
  - `micromamba run -p .../benchenv Rscript -e 'extSoftVersion()...'` -> BLAS `libopenblasp-r0.3.34.so`, `MISSING: <none>` for all 18 Imports
  - `mkdir -p /sys/fs/cgroup/{ancombc_probe,user.slice/probe_test,user.slice/user-7059086.slice/probe_test}` -> all `Permission denied`
  - `R CMD INSTALL --library=.rlib/original reference/ANCOMBC` (R 4.5.3 env) -> `* DONE (ANCOMBC)`; `packageVersion` -> `2.15.2`; all 18 Imports resolve at the recorded versions; 7 exports
  - `R CMD INSTALL --library=<tmp> reference/ANCOMBC` (system R 4.3.3) -> `ERROR: this R is version 4.3.3, package 'ANCOMBC' requires R >= 4.5.0`, exit 1
  - `make exact-env-check PROFILE=linux-r453-openblas` -> 22/22, exit 0
  - `make exact-env-selftest PROFILE=linux-r453-openblas` -> 11/11 negative mutations non-zero, exit 0
  - `make exact-env-check PROFILE=linux-r433-refblas` -> exit 1 (unsupported, by design)
  - `make exact-input PROFILE=linux-r453-openblas FIXTURE=validation/fixtures/fx0{1,2,3} ...` -> three input RDS written
  - `make exact-run PROFILE=linux-r453-openblas CASE=... REPEATS=2` for fx01/fx02/fx03 -> original self-deterministic on all three (2 identical payload hashes each)
  - `make exact-runner-selftest PROFILE=linux-r453-openblas CASE=...` -> 8/8 non-zero
  - `make exact-comparator-selftest PROFILE=linux-r453-openblas` -> 30/30, exit 0
  - `make exact-compare LEFT=<r1> RIGHT=<r2>` -> identical (exit 0); with one flipped byte -> nonzero, first difference located
  - `make exact-cases CASESET=small OBSERVED=...` -> complete exit 0; partial, duplicated, absent-observed and unknown-set all nonzero
  - `make exact-inputs` / `exact-inputs-verify` -> 114 inputs built, then "REGENERATION: byte-identical for every input and every source file", exit 0
  - `make exact-case-manifest-check` -> "up to date (25 small, 114 required)", exit 0
  - `python3 scripts/upstream_manifest.py --verify` -> "45 files, 0 locally modified, UPSTREAM MANIFEST: consistent", exit 0
  - `make replacement-check` -> 24/24, exit 0
  - full two-arm campaign (114 cases, both installed packages, `--repeats 1`): **114/114 byte-identical**, exit 0; `make exact-cases-required` -> "MANIFEST: complete", exit 0
  - small set (25 cases) run and checked separately -> 25/25 byte-identical, "MANIFEST: complete", exit 0
- S01 findings that change later work:
  - The original's `ancombc2` has **32 formals**, in this order: `data`,
    `taxa_are_rows`, `assay.type`, `assay_name`, `rank`, `tax_level`,
    `aggregate_data`, `meta_data`, `fix_formula`, `rand_formula`,
    `p_adj_method`, `pseudo`, `pseudo_sens`, `conservative`, `prv_cut`,
    `lib_cut`, `s0_perc`, `group`, `struc_zero`, `neg_lb`, `alpha`, `n_cl`,
    `verbose`, `global`, `pairwise`, `dunnet`, `trend`, `iter_control`,
    `em_control`, `lme_control`, `mdfdr_control`, `trend_control`. There is **no
    `...`**. Defaults that are lazy: `assay.type = assay_name`,
    `rank = tax_level`, `lme_control = lme4::lmerControl()`. Control-list
    defaults are `iter_control = list(tol = 0.01, max_iter = 20, verbose = FALSE)`,
    `em_control = list(tol = 1e-05, max_iter = 100)`,
    `mdfdr_control = list(fwer_ctrl_method = "holm", B = 100)`,
    `trend_control = list(contrast = NULL, node = NULL, B = 100)`.
    **`pseudo_sens = TRUE`**, confirming the defect recorded against
    `r/ancombc2rs/R/ancombc2.R` (`pseudo_sens = FALSE`).
    Full dump: `validation/exact/evidence/s01_original_public_api.txt`.
  - The package exports exactly 7 names: `ancom`, `ancombc`, `ancombc2`,
    `data_sanity_check`, `secom_dist`, `secom_linear`, `sim_plnm`. C2 requires
    all 7 plus the 32-formal `ancombc2` signature.
  - The pinned original has **no `src/` directory**: it is a pure-R package, so
    `R CMD INSTALL` needs no toolchain and the replacement has no compiled
    original code to preserve.
  - `reference/R/oracle.R`'s `verify_oracle()` guards its commit check with
    `dir.exists(file.path(dir, ".git"))`, and `reference/ANCOMBC/.git` **does not
    exist**, so on this checkout the commit check is silently skipped. The
    per-file `reference/env/oracle.sha256` digests are the real pin; S01's
    verifier enforces those instead. Recorded as `CLAIM-ORACLE-GIT-REV`.
- S02 findings that change later work -- all three are **pinned-original
  behaviour that must be preserved**, not candidate defects:
  - **`fix_formula` is a character string, not a formula object.**
    `data_sanity_check.R:253-254` does `gsub("\\*", "+", fix_formula)` and
    `strsplit(fix_formula, "\\s*\\+\\s*")` on it, and `.ancombc2_sens_fit`
    builds `stats::formula(paste0("y ~ ", fix_formula))`. Passing a formula
    object takes a different code path in the original, so `make_exact_input.R`
    accepts only a string.
  - **A one-column `meta_data` data.frame makes `ancombc2` fail**, and this is
    real: `data_sanity_check.R` ends its matrix branch with
    `meta_data = meta_data[colnames(feature_table), ]`, and `[.data.frame`'s
    `drop` defaults to TRUE when there is exactly one column, so a single-column
    metadata frame collapses to a vector. The next check then reports
    `The following variables specified are not in the meta data: group`.
    Measured on `fx01-basic`: `simpleError`, 744-byte serialized condition,
    identical across two runs. The package's own vignette always passes a
    multi-column `smd`. `fx01-basic` is therefore kept as an **error case**, not
    deleted: reference errors are comparable outcomes (S14).
  - **`verbose` is not forwarded to `data_sanity_check`** (`ancombc2.R:425-441`
    omits it), so that function always uses its own default `verbose = TRUE` and
    its progress messages appear whatever the caller asked for. 19-20
    conditions are recorded per successful run and are part of the comparison.
- First remaining mismatch or blocker: **no mismatch remains in the 114 cases that
  exist.** The blocker is now scope, not correctness: the candidate runs the
  retained upstream R, so nothing has been accelerated and no Rust stage has
  produced a byte. Three environment limits are confirmed and recorded rather
  than worked around:
  1. **8 physical cores.** P5 is defined "on the declared scaling host with at
     least 16 available physical cores". This host has 8. P5 cannot be certified
     here; it needs a runner with >= 16 physical cores. Do not reinterpret the
     16-*logical* figure as satisfying it.
  2. **No cgroup v2 write delegation.** P4 requires
     `peak_cgroup_memory_bytes` from a dedicated cgroup per arm including worker
     processes. `mkdir` is denied at every visible level, and there is no sudo.
     Per S16.4 the fallback is a labelled process-tree sampling *diagnostic* that
     explicitly cannot pass the gate. P4 needs cgroup delegation or a configured
     runner.
- Long-running processes and log paths: none still running. Provisioning log:
  `$MM/logs/mm-create.log`. S01 acceptance log:
  `validation/exact/evidence/s01_profile_preflight.txt`.
- S03 findings worth carrying forward:
  - **`identical(-0, 0)` is TRUE in R.** Any acceptance rule built on `identical`
    or `==` on the object cannot see signed zero; only the serialized bytes
    (`0x8000000000000000` vs `0x0000000000000000`) do. The comparator therefore
    never calls `identical` on a number.
  - **R's parser folds a negated numeric literal inconsistently.** Measured here:
    the literal `-0` inside `check_exact.R`'s `selftest()` folded to `+0`, while
    the same literal in a two-line script did not. A self-test written with the
    literal would have been comparing `+0` to `+0` and passing for the wrong
    reason. Negative zero is now built as `0 * -1`, and
    `signed_zero_fixture_is_sane` asserts `1/base$nz == -Inf` before the two
    signed-zero checks are allowed to mean anything.
  - **`isTRUE(all.equal(1, 1 + 1e-12))` is TRUE.** Recorded as an explicit
    passing check named `tiny_relative_difference_is_invisible_to_all_equal`, so
    the document states what the old comparator would have let through rather
    than merely asserting that the new one is stricter.
- Next bounded action: **S04** -- extend `validation/exact/inputs/` with fx04
  and the parameter/feature cases, each built by `scripts/make_exact_input.R`
  from its committed file read once (never regenerated in place), and add a
  regeneration check that rebuilds inputs and both arms' outputs byte for byte
  in a temporary directory. Also make `run_exact.py --json-out` emit the
  observed-case list in the shape `scripts/check_exact.R --cases` reads, so the
  manifest check consumes real observations.

## S09: least squares and the iterative MLE

**Update (end of session 2):** the least-squares half is done -- see the S09 row above
(191/191 stage checks, 440/440 `Cdqrls` problems on two BLAS families). The text below
is the investigation that led there and is kept as the record of *why* the transport
and the QR had to change; its "Where S09 stands" figures (37 pass, 53 fail) are
superseded. The iterative-MLE half remains open.

### What the reference actually calls

`stats::lm.fit` is not a LAPACK call. It is

```r
z <- .Call(C_Cdqrls, x, y, tol, FALSE)   # tol = 1e-7, LINPACK dqrls / dqrlss
```

and the wrapper then does three things that a least-squares port does not do by
default: it unpermutes the coefficients against the pivots and fills every aliased
slot with `NA_real_`, it sets `df.residual = n - rank`, and it computes
`fitted.values <- y - z$residuals` rather than `x %*% coef`. That last one is not
cosmetic: on a hard case the two disagree bitwise, and on a rank-deficient one they
differ structurally, because `coef` carries `NA`.

Two measurements that shape the port:

- **Rank does not depend on `y`.** Verified by fitting the same `x` against several
  different `y`. That matters because `.iter_mle`'s rank test is a `rnorm()`-seeded
  smoke test on `x %*% beta`; if rank could move with `y`, the iteration's branch
  would be nondeterministic and no care elsewhere would make the trace reproducible.
- **Pivoting is not cosmetic.** On a design with a column scaled by `1e12`, `dqrls`
  returns different coefficients in the low bits than an unpivoted solve, and it
  reports full rank where the current Rust QR reports rank deficiency.

### The three known divergences in the current Rust fit

1. `fit_one_group` computes `dof` as `rows.len() - p`; the reference uses
   `n - rank`. On a full-rank design these agree, which is why `dof` already passes
   on the fixtures that reach this path.
2. It computes fitted values as `x %*% coef`, not as the solver's residuals.
3. Its QR is unpivoted Householder. `crates/ancombc2-core/src/matrix/linalg.rs`
   carries a comment claiming LAPACK `dgeqp3`; that comment is wrong and is corrected
   as part of the port.

### The instrument

`scripts/check_fit_stages.R` compares the pipeline's own `lm_fit_all` against
`.lm_fit_all` transcribed from `reference/ANCOMBC/R/ancombc_bias_correct.R:11-84`,
plus `stats::lm.fit` directly for a single group. It reaches Rust through
`ANCOMBC:::rb_fit_probe()`, which runs the real `lm_fit_all` on the real inputs, so a
difference is a difference in the code that will run. Eight fixtures: full rank, two
factors, a scaled column, a missing count, a taxon with all counts missing, a huge
response, an incomplete design row, and both together.

`rb_fit_probe` also returns the pattern grouping it computed -- the group count and
each group's row and taxon counts -- because a `dof` difference is nearly always
either a row-count difference or a rank difference, and that is a two-line question
rather than a reading of the QR.

### Four bugs in the instrument, and why none of them showed up as an error

Every one of these produced plausible numbers rather than a crash or a NaN. Each was
found by narrowing a *wrong answer* down to a single layer, not by reading code and
spotting it.

**R's `LGLSXP` is four bytes per cell, not one.** The `observed` mask was handed to
the core by casting `LOGICAL()` to `uint8_t *`, on the strength of a comment
asserting that "`uint8_t` and R's logical are both one byte with 0/1". They are not:
`LGLSXP` is an `int32_t` array. The core therefore read the first *quarter* of the
buffer, in a stride of one -- and because `TRUE` is `1` stored as `01 00 00 00`, each
cell read as `TRUE, FALSE, FALSE, FALSE`. An all-`TRUE` 72-cell mask arrived with 18
cells set, which is exactly the number observed. The mask is now narrowed into a real
byte buffer, with `NA` normalised rather than cast (`INT_MIN` would set every bit).

**R's column-major index of `(t, j)` is `t + j * n_taxa`, and both flattenings had it
backwards.** The shim built the row-major response buffer as `ymat[j + t * n_samp]`.
That is not a permutation that looks wrong -- it leaves the buffer in R's own
column-major order, so taxon 1's response arrives interleaved with taxa 2, 3 and so
on, and *every taxon gets a plausible-looking response while every fit is quietly the
wrong one*. `dof` was unaffected, the fitted values were constant within each design
half as a correct fit requires, and nothing errored. The sentinel `y[t, j] = 100t + j`
is what made it visible: the core reported a first-half mean of 351, which is
`mean(101, 201, 301, 401, 501, 601)`. The same mistake was in the `observed`
flattening, where it is invisible on an all-`TRUE` mask and only appears once a cell
is missing.

**The outputs were never transposed.** `beta` and `fitted` come back from the core
**row-major**, because that is the core's own layout and `lm_fit_all` writes a taxon's
row at a time. The shim gave them a `dim` without transposing, so R read down the
taxon axis instead of along the sample axis: `fitted[1, ]` was the interleaved
responses of taxa 1, 3, 5, .... The fix transposes in the shim rather than making the
core lay its output out for R, which would put a presentation concern in the numerics.

**`matrix_arg` hands back `NULL` for a `REALSXP` out-parameter when the matrix is
`INTSXP`.** The design pointer was then NULL and every design value read as whatever
was at address zero. `rb_fit_probe` coerces in R, so in practice `x` is already
`REALSXP`; the shim now widens it explicitly and refuses a null design pointer,
because a NULL design is a silent wrong answer rather than a crash.

The grouping itself was never at fault: `group_by_observation` returns 12 rows for one
group on these inputs, confirmed in a unit test before the transport bugs were found.

### Where S09 stands

With the transport correct, `beta[1]` is now **bit-identical** to `dqrls` on the
full-rank fixture and `beta[2]` differs by 4.2e-16 relative -- one or two ulp. The
fitted values differ by 1.5e-15, and the sum of squared residuals is identical to
R's to all 15 printed digits, with the core's answer the marginally better of the two.
So the remaining `beta`/`fitted` differences on the full-rank path are **rounding, not
structure**: a different pivoting and update order in the QR, which is divergence 3
above and exactly what the LINPACK `dqrls` port is for. `dof` agrees on every fixture
that reaches the grouped path.

Current instrument reading: **37 checks pass, 53 fail**, and every remaining failure is
in `beta` or `fitted` at the ulp level. Nothing structural is left in the transport.

Two things are still open, and they are the real work of S09:

1. **Bit-exactness on the full-rank path.** The ulp differences need the actual
   `dqrls` algorithm -- LINPACK's `dqrlss`, with its own pivoting, its own
   `d1mach`-style tolerance and its own accumulation order -- not a more accurate QR.
   Any correct QR agrees with any other to rounding; byte-identical output requires the
   *same* rounding. This is the port, and it also has to settle the corrected `dgeqp3`
   comment and the rank-deficient cases where `dqrls` reports full rank and the current
   QR does not.
2. **The rank-deficient fallback.** The Rust side has a per-taxon reduced-design
   fallback in `fit_one_group`; the reference falls back to
   `stats::lm(tformula, data = meta_data)`, which brings `model.matrix`,
   `contr.treatment` and dropped factor levels with it. The transcription refuses
   fixtures that reach that path without a formula rather than silently agreeing with
   them, so `missing-all-in-a-taxon` and `incomplete-row-and-missing` are not yet
   comparable at all.

## Release gates

| Gate | State | Evidence | Remaining work |
| --- | --- | --- | --- |
| C1 original fidelity | **in_progress** | installed original at `.rlib/original`, identity digest `fd42b62f1464f428...` matching the S01 pin on every one of 114 runs; `R CMD INSTALL` succeeded; all 18 Imports at the recorded versions; 22/22 preflight checks; per-case repeatability established on fx01/fx02/fx03 (2 runs each, identical payload digests) | S01/S02 done. Outstanding for a full C1: a repeatability sweep over all 114 cases (the campaign ran one repetition each; `--repeats 2` exists) and the dependency manifest as a committed artefact rather than a profile field |
| C2 package replacement | **in_progress** | `make replacement-check` -> **24/24, exit 0**: identical 32-formal names/order, byte-identical `args()` text and per-formal default expressions, `pseudo_sens = TRUE`, no `...`, the original's 7 exports and nothing more, every export's signature identical, attribution/licence retained, `X-Replacement-Of` present; **no namespace loaded from the original's library** and the original's library absent from the candidate's search path | standalone `R CMD build` + `R CMD INSTALL` from the archive without the checkout is still S20. The `required` case set is not yet the one C2 names (no real datasets, no phyloseq/TSE input) |
| C3 exact results | **in_progress** | 114/114 whole-result **byte-identical** against the installed original on the primary profile, via `serialize(x, NULL, ascii = FALSE, xdr = TRUE, version = 3)` with no normalisation; comparator's 30 negative checks all reject; missing/duplicated cases rejected | This is the *scaffold* agreement -- both arms currently run retained upstream R, so it proves the harness, the inputs and the schema, not Rust parity. C3 is not passed until the same 114 cases (plus the S14 additions) pass with the accelerated stages executing. S08-S13 remain |
| C4 execution integrity | **in_progress** | runner refuses a library that does not exist, one with no ANCOMBC, the original run as the candidate and a replacement run as the original (8/8 guard checks); the candidate arm resolves only its own package; each arm gets `R_LIBS = <its own library>:<site>`; mid-campaign package-identity changes are detected | no Rust backend evidence yet, so "backend evidence proves the intended Rust path" is not satisfiable. S08 checks preprocessing stages and still computes no fit; S09-S13, then S20 |
| B1 benchmark validity | pending | none; `scripts/bench_r.R` omits `cfg`, `reference/R/harness.R` stubs `%dorng%` sequentially | S15, S16 |
| P1 kernel performance | pending | historical 1.141x (invalid for comparison: sequential R stub, kernel/total boundary conflated) | S17 onward |
| P2 complete-call performance | pending | historical 1.141x / 2.00x (same invalidity) | S17 onward |
| P3 sensitivity performance | pending | historical 3.165x on bm5, plus a 5.50x bm6 scaling figure (invalid) | S17 onward |
| P4 memory | pending | historical 1.751x peak-RSS ratio; **metric is incompatible** (`gc()` vs `VmHWM`) and the required cgroup metric is unmeasurable here | S16.4 infrastructure + S18; needs cgroup delegation |
| P5 scaling | pending | historical 0.172 on bm5 / 0.172-5.50x range (invalid) | needs >= 16 physical cores; host has 8 |

## Environment blockers (external changes required)

| id | blocker | affects | precise external change needed | workaround permitted as evidence? |
| --- | --- | --- | --- | --- |
| ENV-1 | 8 physical cores, 16 logical | P5 | provision a runner with >= 16 physical cores and record `lscpu` in the campaign manifest | no; a 16-logical-core number is not P5 |
| ENV-2 | no cgroup v2 write delegation, no sudo | P4 | delegate `memory` on a cgroup subtree to this account, or run P4 on a configured runner | yes, as a labelled process-tree sampling *diagnostic* only (S16.4) |
| ENV-3 | no container engine (`docker`/`podman`) | the `benchmarks/container/Dockerfile` path | install an engine, or use the realised lock environment directly | yes; the lock environment is what the Dockerfile realises |
