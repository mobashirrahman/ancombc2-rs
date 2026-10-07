# Plan: exact ANCOMBC replacement with measured Rust acceleration

Date: 2026-10-05. Status: implementation has not started.

This plan supersedes the scope and release gates in `PLAN.md` and the optimization
recommendations in `docs/performance_plan.md`. Preserve those documents and their
results as history. Historical measurements are not the new release baseline.

## 1. Required outcome

Install a replacement package and run existing code such as
`ANCOMBC::ancombc2(...)` without changing the call or its downstream consumers.
The replacement must return exactly the same result as the pinned original on
the same input and runtime. Acceleration must be measured on that installed
replacement, including its R adapter and native bridge.

The original is ANCOMBC 2.15.2, commit
`dc4febdf59badb3a8dfe0c767ef2186323c2199a`.

The full `ancombc2()` interface is in scope: fixed and random effects, global,
pairwise, Dunnett and trend tests, both sensitivity modes, taxonomic aggregation,
matrix/data.frame/phyloseq/TSE input, aliases, controls and default behavior.
Replacing a package must also preserve its other exports. Keep their licensed
original implementations; acceleration work focuses on ANCOM-BC2.

### Architecture choice

Build an **exact hybrid implementation**: retain the original R-facing code and
dependency routines where needed, and replace validated numerical stages with
Rust. This avoids requiring a small model to independently recreate `lmer`, R's
formula semantics, numerical optimizers and the entire package namespace.

The shipped replacement must be self-contained apart from its declared ordinary
R dependencies. It must not load the original ANCOMBC installation, source files
from `reference/`, fetch upstream code at runtime, or run an oracle to correct its
answers. Bundled upstream code must retain its Artistic-2.0 license, attribution,
provenance and checksums. Retain each imported native routine's own license and
provenance as well; prefer calling an installed dependency when appropriate.
Describe retained R stages accurately.

Keeping original code is a functional milestone. It is not completion of the
acceleration objective. Required benchmark cases must demonstrably execute Rust
numerical work, and performance gates must pass for the shipped execution path.

### Meaning of exact

For each certified runtime profile and identical initial state:

1. Match public formal arguments, their order, aliases and default expressions.
   Include positional calls and omitted arguments. Do not add an accepting `...`
   to the replacement entrypoint if the original does not have one.
2. Match every returned field, including `NULL` fields, list order, classes,
   attributes, row/column names, shapes, element order and storage types.
3. Match the raw representation of numbers, including signed zero and the
   distinction between `NA_real_`, `NaN`, infinities and ordinary doubles.
4. Match warnings/errors/messages and their order for the exercised behavior
   cases. Check reference-observable RNG state and options after the call, input
   mutation and worker cleanup. Runtime durations are inherently variable;
   compare progress-message templates and numerical content, never use a broad
   message filter to hide differences.
5. For successful calls, require identical complete-result serialization:

   ```r
   result_bytes <- function(x) serialize(x, NULL, ascii = FALSE,
                                        xdr = TRUE, version = 3)
   stopifnot(identical(result_bytes(reference_result),
                       result_bytes(candidate_result)))
   ```

Do not strip fields, reorder rows, round values, coerce types, remove attributes,
or replace missing values before this check. A diagnostic normalized comparison
may help find a mismatch; it cannot determine acceptance. `all.equal()` alone is
insufficient. A matching hash is convenient evidence; also compare raw bytes.

The byte identity contract concerns returned analysis objects and explicitly
defined analysis exports. Package binaries, installation metadata, benchmark
telemetry and elapsed-time text are not analysis results. Preserve all original
result fields. Keep new backend/timing telemetry in separate files.

Exactness is comparison against the original **in the same runtime**, not against
one universal answer across different BLAS/R/OS versions. Certify explicit
profiles. Do not infer support for every `R >= 4.5.0` from two tested versions.
Each additional supported profile gets its own exact checks against its original.

### Rules that protect the requirement

- Never accept an epsilon, ULP allowance, rounding rule or rank-deficiency
  exemption in an exact release gate.
- Preserve the pinned original even where its behavior is mathematically unusual.
  Statistical improvements belong in a separately named experimental API.
- Fix the first differing numerical stage. Do not compensate in a later table.
- Do not populate candidate results from golden files or consult the reference
  during candidate execution.
- Keep `ancombc2-core` usable without R. Put R-specific exact primitives and
  adapters behind a separate native integration boundary.
- Never call the R API from Rayon workers. Use R on its main thread, an audited
  thread-safe native primitive, or the original process backend as appropriate.
- Parallelize independent work and preserve each result's operation order. A
  reduction rewrite or FMA can change bits even when the algebra is equivalent.
- A missing prerequisite, skipped required case, timeout, mismatched result or
  invalid measurement prevents the corresponding gate from passing.

## 2. Findings to resolve

These findings were checked on this checkout. Reproduce them before changing
their relevant files; several older documents contradict the current code.

| Finding | Current evidence | Required closure |
|---|---|---|
| Numerical equality is approximate | `crates/ancombc2-core/tests/golden/mod.rs` uses tolerances; fresh fx03 differences: beta `2.84e-11`, p `1.58e-10`, q `3.89e-10` | Whole returned objects pass strict bytes against the installed original |
| Wrong sensitivity default | `r/ancombc2rs/R/ancombc2.R`: `pseudo_sens = FALSE`; pinned original: `TRUE` | Exact public signature and omitted-argument tests |
| Missing features | Random effects, Dunnett/trend and aggregation lack implementation in the Rust wrapper | All original options execute through the replacement and pass exact checks |
| Different namespace and return schema | Package is `ancombc2rs`; FFI exposes `global`, `pairwise`, `sensitivity` | `ANCOMBC::ancombc2` works with the original complete schema |
| JSON bridge is unsuitable for a raw-bit contract | `call_core()` uses JSON; nonfinite values become null and R reconstruction changes structure | A typed native bridge preserves values and result attributes |
| Old fixtures cannot establish exact input equality | `Makefile` documents goldens computed before lossy TSV round trips | New goldens originate from the exact input read by both arms |
| R benchmark ignores configuration | `scripts/bench_r.R` calls `H$ref_run(counts, meta, fix_formula)` without `cfg` | Same argument object is consumed and recorded by both installed packages |
| R parallel benchmark is sequential | `reference/R/harness.R` installs `foreach_stub` and a sequential `%dorng%` | Benchmark the real installed original with verified workers |
| Memory ratios use incompatible metrics | R records current `gc()` usage; Rust records process `VmHWM` | Measure both with the same external metric including workers |
| Gate records can be mixed and unrepresentative | One repetition; latest rows keyed only by dataset/arm; kernel and total gates share wall time | Immutable campaigns, repeated paired measurements and distinct boundaries |
| Source installation relies on the workspace | `r/ancombc2rs/src/Makevars` uses `../../..` | Standalone source archive installs with the repository unavailable |
| Current claims overstate evidence | Documents report full completion and bit identity alongside tolerant or limited checks | All public claims are generated from the new evidence manifest |

## 3. Release gates

Freeze these definitions before the corrected baseline. Changing a definition
requires a recorded rationale and a new campaign; never reinterpret old rows.

| Gate | Required result |
|---|---|
| C1: original fidelity | Original is an installed, checksum-verified package; full dependency/runtime manifest; original repeatability is checked per case |
| C2: package replacement | Same namespace, exports, ANCOM-BC2 formals/defaults, supported features and complete returned-object schema; standalone installation succeeds |
| C3: exact results | Zero differing result bytes across required fixtures, behavior cases, real datasets and seeded generated cases; no exemptions or skipped required cases |
| C4: execution integrity | Candidate is independent of the original installation; backend evidence proves the intended Rust path; requested thread counts are honored and matching reference configurations agree |
| B1: benchmark validity | Same inputs/configuration/runtime, verified workers, matching output bytes, complete repetitions, comparable timing/memory boundaries |
| P1: kernel performance | At least 3x one-worker speedup on bm5's matched fixed-effects numerical kernel, with equal measured boundaries |
| P2: complete-call performance | At least 2x one-worker speedup on each of bm4, bm5 and bm6, including R adaptation and returned-object construction |
| P3: sensitivity performance | At least 5x sensitivity-stage speedup on bm6 at 16 workers compared with the original at 16 workers; report conservative bm5 separately |
| P4: memory | At most 70% of the original's peak externally measured process-group memory on bm5, at both one and 16 workers |
| P5: scaling | `T(1)/(16*T(16)) >= 0.70` on both bm5 and bm6 on the declared scaling host with at least 16 available physical cores |

P1-P5 are engineering targets, not predicted outcomes. Performance is measured
after exactness. The conservative grid's three additional refits limit that
particular outer loop; they do not prove a universal 3x speedup ceiling against R
or forbid parallelism inside a refit. Do not use that claim to dismiss a deficit.

Report cold process launch, warm complete calls, kernel time, sensitivity time,
one-worker versus one-worker, and matched-worker comparisons separately. P1 is
kernel time; P2 is the median warm complete-call time. Also report regressions on
small inputs and on retained R feature paths. An overall performance completion
claim additionally requires no reproducible complete-call regression over 5% on
those declared workloads, assessed with paired measurements and uncertainty.

Use at least seven measured paired repetitions after two warmups. Keep every
sample, failures included. Use medians and interquartile ranges; report paired
bootstrap confidence intervals using a fixed reporting seed. A release speed
gate uses the lower 95% confidence bound; a memory upper-bound gate uses the
upper 95% bound. A point estimate that crosses a target without that evidence is
inconclusive. More repetitions may resolve uncertainty, but never keep only the
best runs.

A functional release and an accelerated release are separate named milestones.
The overall task is complete only when C1-C4, B1 and P1-P5 pass. If a performance
target proves unattainable with exact execution, report the evidence and the
remaining work. Do not declare completion or silently lower the target.

## 4. Execution protocol for a small model

Implement tasks in the order below. Work on one task or one numbered substep at
a time. Read the relevant original function before editing its counterpart.
Prefer retaining a proven routine over inventing equivalent numerics.

Create `IMPLEMENTATION_STATUS.md` in S00. Give each task a state (`pending`,
`in_progress`, `passed`, `blocked`), changed files, exact commands, evidence paths
and the next substep. `passed` requires its acceptance checks. A checkpoint is
not completion. Do not change this plan into a report that asserts unrun checks.

For each task:

1. Inspect its named files and current implementation status.
2. Write one narrowly scoped failing check for the actual defect or contract.
3. Make the smallest implementation change that resolves it.
4. Run that check and the existing checks for affected code.
5. Record results, including failures, and continue to the next unpassed task.

Run a broad suite when a shared numerical stage changes or at a phase boundary.
Do not repeatedly rerun hours of simulations after unrelated documentation edits.
When computation must run for hours, keep the process/logs and checkpoint its
state; a timeout is not permission to accept fewer required cases.

### S00 — Preserve the baseline and expose unproven claims

**Read:** `git status`, `README.md`, `docs/release_status.md`,
`docs/compatibility.md`, `docs/numerical_contract.md`.

**Do:** record the starting commit and dirty files. Pre-existing modifications
include `scripts/bench_gates.py` and its tracked Python cache; preserve them.
Create the status file and an audit manifest. Mark current R-relative benchmark
and byte-identity claims as unvalidated for this contract. Keep historical raw
results. Add a clear distinction between Rust self-determinism and equality to
the original. Do not repair every old paragraph in this task.

**Accept:** no pre-existing changes are lost; public status links to this plan
and accurately states that exact replacement and corrected performance are
unfinished. All task states begin pending except tasks actually executed.

### S01 — Pin and verify execution profiles

**Read:** `reference/env/ORACLE.md`, `reference/env/oracle.sha256`,
`benchmarks/container/environment-linux-64.lock`, container scripts and CI.

**Do:** add `validation/exact/profiles/` and a machine-readable profile verifier.
Start with the existing Linux R 4.5.3/OpenBLAS lock, after checking its complete
dependency set. Record exact R, BLAS/LAPACK, compiler, Rust flags, architecture,
locale, thread environment, RNGkind, dependency versions and package/source
hashes. Install the unmodified pinned original in an isolated library. Add a
reference-BLAS profile and any additional supported OS/R profiles individually;
never reuse goldens across profiles. Validate required dependencies before runs.

**Accept:** wrong oracle SHA, R version, BLAS, package hash or missing required
dependency causes a nonzero preflight result. The installed original reports
the pinned version and resolves every declared dependency. Unsupported profiles
are identified before analysis rather than treated as exact-certified.

### S02 — Create a real original/candidate runner

**Read:** `scripts/verify_real_package.R`,
`reference/ANCOMBC/R/ancombc2.R`, `reference/ANCOMBC/NAMESPACE`.

**Do:** create `scripts/exact_runner.R` and `scripts/run_exact.py`. Launch the
original and candidate in separate fresh R processes with isolated library paths;
they will share the package name ANCOMBC. Load one identical arguments/state RDS
file in both. Call the installed exported function. Capture the full result,
conditions, initial/final RNG state and selected observable state. Record package
paths/hashes. Do not source `reference/R/harness.R` or install loop stubs in this
runner. Create only the runner now; candidate installation follows in S05.

**Accept:** two runs of the installed original over a deterministic fixture
produce equal full-result bytes. The runner rejects a missing/wrong library and
cannot accidentally resolve the original as the candidate. Unexpected reference
nondeterminism becomes a diagnosed failure rather than an ignored case.

### S03 — Implement the strict comparator and its negative checks

**Do:** add `scripts/check_exact.R`. Use complete serialization equality for
acceptance, plus a recursive diagnostic naming the first field, attribute and
element that differs. For numeric mismatches report raw hex bytes and values;
do not use the diagnostic comparator to allow discrepancies. Capture required
case counts in `validation/exact/cases.json` and reject absent/duplicate cases.

**Accept:** deliberately modifying one mantissa bit, signed zero, NA/NaN kind,
integer/double type, a name, row order, attribute, NULL field or list order fails
the check. Altering a warning or RNG state fails its behavior check. Identical
objects pass. A skipped or missing required case produces nonzero exit status.

### S04 — Establish exact inputs and independent goldens

**Read:** `scripts/generate_goldens.R`, `reference/R/fixtures.R`,
`reference/R/serialize.R`, the drift comments in `Makefile`.

**Do:** create exact fixtures under `validation/exact/`, leaving legacy goldens
untouched. For existing fixtures, read their committed input files once through
the declared R path, then save the argument object/state as RDS. The installed
original must compute its exact golden from that same saved input. For new
fixtures use explicit seeds and save the input before computing either result.
For CLI cases, compare parser output bytes before comparing analyses; use a
round-trip-safe input format or writer. Store input/configuration/profile hashes.

**Accept:** regeneration in a temporary directory reproduces exact fixture
inputs and outputs. Both arms consume the same bytes. No test compares a Rust
run over rounded TSV input to an original run over earlier unrounded values.
Candidate execution has no access to the oracle output directory.

### S05 — Build a self-contained replacement namespace

**Read:** pinned `DESCRIPTION`, `NAMESPACE`, public R code, package data and
license; `r/ancombc2rs/src/Makevars`.

**Do:** create a replacement package under `r/ANCOMBC/` from the licensed pinned
source. Preserve original exports, dependencies, public preprocessing, argument
handling, result construction and unaccelerated methods. Add provenance notices
and a tracked upstream file/hash manifest. Keep `reference/ANCOMBC` immutable.
Use separate installation libraries for original and replacement. Keep candidate
distribution metadata truthful; API compatibility does not mean pretending its
build provenance is the original's.

**Accept:** existing named and positional `ANCOMBC::ancombc2()` calls work after
switching only the library installation. Its `formals()` match the original,
including `pseudo_sens = TRUE`, lazy aliases and control-list defaults. Every
original export is present. Original versus candidate calls using retained code
pass strict bytes. This is explicitly the compatibility scaffold, not Rust
acceleration. No dependency on an installed original ANCOMBC exists.

### S06 — Inventory the numerical seam and build native input transport

**Read:** `.ancombc2_core` in the pinned R source,
`crates/ancombc2-ffi/src/lib.rs`, `r/ancombc2rs/src/init.c`,
`crates/ancombc2-core/src/pipeline.rs` and `workspace.rs`.

**Do:** write a field/shape/type/layout inventory for the actual core seam.
Create a separate R-native adapter crate/C shim as needed. Keep the public
entrypoint and its R semantics intact; replace a numerical internal call rather
than manually reconstructing the entire public workflow. Pass typed numeric
buffers and explicit metadata through `.Call`, with documented ownership,
protection and layout conversions. Preserve R NA/NaN tags. Keep R pointers out
of Rayon jobs and convert failures into normal R conditions safely.

**Accept:** transport-only round trips preserve every input value bit, dimensions,
names, masks and integer/group metadata. Include NA, NaN, infinities and signed
zero. Layout and lifetime checks pass. No JSON/null conversion is used on the
exact package path. The pure Rust library still builds without R.

### S07 — Build native output transport and exact result assembly

**Read:** pinned core and public return lists; current FFI output and
`r/ancombc2rs/R/tables.R`.

**Do:** return precisely the internal objects the retained public R code expects.
Construct typed matrices/data.frames on the R thread with the original storage
types and attributes. Preserve bias-estimation and reported taxon sets separately.
Use original assembly code for `feature_table`, `bias_correct_log_table`,
`ss_tab`, `zero_ind`, `samp_frac`, `delta_em`, `delta_wls`, `res`, `res_global`,
`res_pair`, `res_dunn` and `res_trend`, including absent-test NULL entries. Keep
extra diagnostic fields outside the public result.

**Accept:** a transport/assembly test with original numerical payloads passes
whole-result bytes. A consumer reading every original field needs no changes.
This test proves transport/schema only; it must not be reported as Rust
numerical parity. Do not use original payloads in production candidate runs.

### S08 — Fix preprocessing and reduction differences first

**Read:** `crates/ancombc2-core/src/preprocess.rs`, pinned `.data_core` and R's
actual log/centering/reduction call sites.

**Do:** compare identical stage inputs before comparing outputs. Add raw-bit
stage diagnostics for filtering, transformed counts, centering, design and
structural zeros. Match scalar math, missing-value rules, reduction order and
accumulator precision. Reuse R's operations on the main thread if the current
Rust calculation cannot match; put any shared primitive behind the adapter.
Do not assume Rust `sum()` matches R `rowMeans` or BLAS.

**Accept:** stage arrays and retained identifiers are exact on tiny, sparse,
missing-value and boundary fixtures. A failure is reduced to an explicit small
input and fixed here before continuing downstream.

### S09 — Match least squares and iterative MLE

**Read:** `mle.rs`, `matrix/linalg.rs`, pinned `.lm_fit_all`/`.iter_mle`, and
`docs/reference_behavior.md` section 16.

**Do:** compare one fit, one iteration, then the complete iteration trace.
Match the original QR/solve primitive and operation order, rank thresholds,
dropped factor levels, coefficient names, fitted projections, aliased NA versus
absent-column zero, degrees of freedom and stopping behavior. The current custom
Householder QR is not presumed bit-equivalent to the original primitive. Use
the original runtime primitive or an attributed exact port of its actual routine.
Avoid a general QR rewrite. Keep per-pattern cache keys faithful to the complete
design and solver controls.

**Accept:** full-rank, singular, underdetermined, one-level and interaction cases
agree at each iteration and in their final arrays. There are zero rank-deficient
exemptions. Unknown behavior is resolved by a minimized call to the installed
original, not by a tolerance or a mathematically preferred answer.

### S10 — Match covariance and generalized inverses

**Read:** `vcov.rs`, pinned `.sandwich_vcov`, its `MASS::ginv` call sites,
and actual block/reduction ordering.

**Do:** compare residuals, each taxon's covariance block, inverse decisions and
variance outputs. Keep BLAS/LAPACK choice, singular-value cutoff, matrix layout
and accumulation semantics equal to the original. A direct reuse of the original
primitive is acceptable. Preserve exactness before adding cached outer products
or blocked scratch buffers.

**Accept:** ordinary and singular covariance arrays are raw-bit identical on the
required stage cases. The final object retains the reference's missing values
and types. Serial and blocked implementations agree exactly when both exist.

### S11 — Match EM, optimizer behavior and correction

**Read:** `em.rs`, `correct.rs`, `ancombc2-stats/src/nelder_mead.rs`, pinned
`.bias_em` and the WLS/correction call sites.

**Do:** compare initialization, one E step, one M step, optimizer objective/trace,
stopping decisions, component parameters, variance, delta values and correction.
Use the same installed `nloptr` routine and options where necessary; a Rust
Nelder-Mead implementation is not automatically equivalent. Match density
functions and reductions before changing the optimizer. Preserve the original
behavior on nonconvergence and unusual mixtures.

**Accept:** exact stage arrays and convergence decisions on both sensitivity
modes, sparse cases and explicit optimizer-boundary cases. No observed-epsilon
or nonconvergence allowance is permitted in the exact gate.

### S12 — Match regularization, distributions and multiple testing

**Read:** `stats.rs`, `test_mod.rs`, `ancombc2-stats/src/dist.rs`,
`padjust.rs`, `quantile.rs` and the corresponding original call sites.

**Do:** compare `s0`, W statistics, distribution calls and all seven adjustments.
Use the original Rmath/statistics routines or an exact attributed port; do not
substitute algebraically different survival-function formulas. Match global and
pairwise contrast order, covariance calculations and mdFDR control defaults.
Include values next to significance thresholds and extreme tails.

**Accept:** exact raw p/q/W/SE/regularization values and every discrete flag.
Passing significance decisions while p/q bits differ fails. Test all adjustment
methods and both sides of each comparison boundary.

### S13 — Complete full feature coverage and dispatch integration

**Read:** retained public R code, `sens.rs`, original random-effects,
Dunnett/trend, aggregation and input dispatch paths.

**Do:** integrate exact stages into the replacement's actual workflow. Keep
original dependent routines for random models and constrained tests until an
exact accelerated implementation exists. Preserve input coercion, taxonomy,
factor contrasts, aliases and all controls. Produce separate backend evidence
showing stages executed in Rust versus retained R. Document dispatch conditions
as deterministic feature checks; never run both implementations and choose the
answer that matches. The default package cannot accept an option and ignore it.

**Accept:** strict tests cover random intercepts/slopes, Dunnett, trend with each
supported control combination, aggregation, all input classes, aliases,
positional/default calls and combinations with sensitivity. Fixed-effects
benchmark cases execute real Rust numerical stages. Retained R features are
supported and labeled; they are not mislabeled as pure Rust ports.

### S14 — Expand exact coverage and check thread behavior

**Do:** extend the case manifest and runner to include all 38 existing matrix
cells, seven edge-case inputs, fx01-fx04 and four real datasets. Add deterministic
generated cases (at least 100 seeds), missing metadata, custom contrasts,
rank/conditioning boundaries, library/prevalence cutoffs, alpha boundaries and
nonconvergence. Add small valid and invalid cases for every public parameter.
Store each failure as a small permanent regression input. The manifest must name
every case and expected outcome; reference errors are comparable outcomes,
not skipped cases. Add other-export smoke cases to protect the package facade.

Compare candidate versus original separately at `n_cl = 1, 4, 8, 16` on the
declared host. Check candidate thread invariance only where the original is
invariant. Respect original thread-dependent RNG behavior when present. Verify
the worker count through observed workers/execution, not a requested-count field.
Exercise different thread counts in consecutive calls within one R session so
a once-created global Rayon pool cannot silently ignore later requests.

**Accept:** all required differential and worker checks pass for every required
profile/case. Zero required omissions, zero result-byte mismatches, zero
unhandled feature options. Packaging and final runtime-isolation portions of
C2/C4 remain pending until S20. The existing tolerance suite remains a diagnostic
suite. Statistical FDR/power agreement is additional evidence and cannot
substitute for exact per-call results.

### S15 — Repair benchmark configuration and execution paths

**Read:** `scripts/bench_r.R`, `scripts/run_benchmarks.py`,
`scripts/make_bench_data.py`, `scripts/bench_gates.py` and container entrypoints.

**Do:** use installed original and candidate packages from S02 for the primary
comparison. Normalize existing dataset config keys once into a complete R
argument object, including `group`, `global`, `pairwise`, `pseudo_sens`,
`conservative`, `struc_zero` and all numerical controls. Decode configuration
with a JSON parser, preserving strings; remove the regex scalar parser. Send
the same argument object to both packages. Distinguish omitted defaults from
explicit values in default-behavior cases. Record input/config hashes and the
effective requested/observed workers. Use the real original backend. Explicitly
pin nested BLAS/OpenMP threads to the declared budget.

**Accept:** changing sensitivity, conservative mode or requested tests changes
the effective configuration of both arms. Stage-presence/refit-count checks prove
the work occurred. The pre-benchmark full-result equality check passes. R
parallel runs prove multiple workers. The installed candidate proves Rust work.
Any mismatch prevents a ratio from being scored. CLI-only timings are a separately
labeled surface and cannot stand in for package performance.

### S16 — Repair timing, memory and campaign bookkeeping

**Do in separate substeps:**

1. Add a versioned benchmark schema and immutable campaign ID. Record git SHA,
   dirty diff hash, binary/package hashes, dataset/input/config hashes, runtime
   profile, CPU affinity, physical/logical cores, memory limit and worker budget.
   All paired arms must match these fields except the intended implementation.
2. Store each repetition individually. Alternate or deterministically randomize
   paired arm order. Record warmups separately. Reset identical initial RNG/state
   before each call. Separate cold startup from warm complete-call timing.
3. Add equal kernel/sensitivity stage timing boundaries. Prefer external or
   minimal stage instrumentation whose results have been exact-checked against
   the uninstrumented installed package. Aggregate thread/process work by elapsed
   stage interval, not a sum of worker wall times.
4. Measure memory externally for both arms. For single-process RSS, record each
   fresh process's true high-water mark. For a parallel package use a dedicated
   cgroup v2 per arm and its reset `memory.peak`, including worker processes, as
   the P4 process-group metric. Record cache/accounting policy and cgroup baseline.
   Name it `peak_cgroup_memory_bytes`, not RSS. Run both arms under identical
   policy. If cgroup delegation is unavailable, retain a labeled process-tree
   sampling diagnostic and execute the required P4 campaign on a configured
   runner; sampled peaks or parent-only peaks cannot pass that gate.
5. Replace latest-dataset/arm selection with explicitly selecting one complete
   campaign. Reject incompatible schema/environment/build/config records. A
   failed repetition remains failed and cannot be replaced by an older success.

**Accept:** synthetic allocating/sleeping parent-and-worker programs demonstrate
the timer boundaries and process-group accounting. Repeating a smaller workload
after a larger one does not produce zero memory. Mixed campaigns, missing
samples, unverified workers and output mismatches are rejected. Kernel and total
gates read different measured fields. No `gc()` or cumulative-rusage subtraction
is accepted as comparable peak process-group memory.

### S17 — Run the corrected baseline and choose work from profiles

**Do:** run a small paired smoke campaign, then the complete bm1-bm6 surface using
the installed packages, declared one/matched-worker budgets and repetition
protocol. Include declared retained R feature workloads. Save raw measurements,
exactness evidence and environment manifests under a new campaign directory.
Regenerate P1-P5 from that campaign only. Keep old rows as invalid-for-comparison
historical evidence. Profile the dominant candidate stages with reproducible
instrumentation and matched boundaries; commit the instrumentation or script.

**Accept:** all benchmark-validity checks pass; the report states each gate's
status and uncertainty. A reproducible profile identifies the largest elapsed
time and live allocation contributors. Do not carry forward the old unsupported
claim that 84% of a matched run is QR without reproducing it here.

### S18 — Reduce live memory without changing the calculation

**Read:** `pipeline.rs`, `workspace.rs`, `sens.rs`, `mle.rs`, `parallel.rs` and
the corrected allocation profile.

**Do one change per patch:** share immutable counts/design; remove redundant
transport copies; shorten buffer lifetimes; reuse scratch after its last reader;
avoid retaining complete refit pipelines when the sensitivity result needs only
their q/p decisions; retain only required cached design/QR state. Bound concurrent
sensitivity pipelines by a declared memory budget and exploit independent work
inside refits when appropriate. Apply results in original order. Do not change
grids, omit tests, reduce iterations, approximate covariance or drop returned data.

**Accept per change:** affected exact cases, large/sparse exact cases and the
targeted external memory measurement pass. Record peak memory and complete-call
time before/after. Re-run the complete campaign when the optimized implementation
is ready. P4 must pass on the required workload/budgets; memory reduction cannot
be claimed using a changed metric or less requested work.

### S19 — Accelerate the measured dominant work while preserving bits

**Do one experiment per patch:** remove R/JSON marshalling overhead already
identified; cache only identical repeated factorizations with complete keys;
batch independent responses against the same original factorization; move
per-taxon orchestration into Rust; parallelize independent taxa/patterns/refits
with per-item scratch. Preserve each scalar sequence and original solver choices.
Vectorize across independent outputs when that preserves each reduction.

Do not start with four-accumulator sums, fast-math, FMA changes, approximate
distributions, early termination or a generic optimized QR backend. Earlier
blocked-QR work measured no improvement and reassociated reductions broke even
the tolerant contract (`docs/performance_plan.md`, sections 7.1 and 7.3).
Revisit an experiment only with a new profile or a changed workload that justifies
it. Mathematical equivalence is not evidence of byte equality.

**Accept per experiment:** exact stage checks and complete-result checks pass;
paired measurements show an improvement on its target workload; relevant retained
R paths do not regress. Reject an experiment that changes bits. Continue profiling
and implementation until the declared performance gates pass or there is a
concrete, reproducible engineering blocker; never invent a speed projection.

### S20 — Verify the distributable package and runtime isolation

**Do:** make source archives include the required Rust sources/lock/configuration,
license and native bridge. Eliminate `../../..` workspace assumptions. Build and
install the actual archive in a temporary directory without the checkout or
original library on its search path. Verify native library lookup and declared
dependencies for each certified platform. Do not silently skip native compilation.
Run existing user scripts against the installed artifact. Repeat exact and
benchmark smoke checks from that artifact and record its checksum.

**Accept:** `R CMD build`, standalone `R CMD INSTALL` and applicable
`R CMD check` checks succeed. The built replacement resolves only its own
ANCOMBC namespace and declared dependencies. Existing ANCOM-BC2 consumer scripts
run with only the library installation changed. Exact results and backend
evidence match the tested release artifact.

### S21 — Wire enforceable CI and publish one consistent status

**Do:** add explicit Makefile targets for profile preflight, exact comparison,
replacement installation, benchmark validity and release checking. Wire strict
small cases into normal CI; full fixtures, feature/profile matrix and packaging
into required release CI; long simulations into scheduled/release evidence where
appropriate. Run actual benchmark infrastructure on the declared runner, including
cgroup memory and the required physical-core budget. Missing tools make a job
unmet/failed, not green. Upload failure logs and partial raw records.

Update `README.md`, `CHANGELOG.md`, `docs/compatibility.md`,
`docs/numerical_contract.md`, `docs/reference_behavior.md`,
`docs/release_status.md`, `docs/reproduction.md`, `docs/performance_plan.md`,
`reference/env/ORACLE.md` and container documentation from the new manifests.
Replace stale claims about test counts, simulations, R versions, incomplete
features and performance ceilings. Cite exact commands/artifact hashes instead
of saying every functional criterion is met from a tolerant test run.

**Accept:** one release manifest contains C1-C4/B1/P1-P5, all required case
counts, runtime profiles, artifact hashes and campaign ID. CI fails if any
required gate is absent, invalid, inconclusive or failed. The shipped artifact,
documentation and measured implementation agree. All gates must pass before
describing the objective as completed.

## 5. Proposed command interface

These targets do not exist yet. Implement them in the tasks above; do not run
them and treat a missing command as an implementation result.

```sh
make exact-env-check PROFILE=linux-r45-openblas
make replacement-build PROFILE=linux-r45-openblas
make exact PROFILE=linux-r45-openblas CASESET=small
make exact PROFILE=linux-r45-openblas CASESET=required
make bench-validate PROFILE=linux-r45-openblas
make bench-exact PROFILE=linux-r45-openblas CAMPAIGN=baseline-v2
make release-check PROFILE=linux-r45-openblas CAMPAIGN=baseline-v2
```

The implementation must record the actual interpreter, library paths, commands
and output locations in `IMPLEMENTATION_STATUS.md`. All generated oracle outputs
go to a temporary directory or the new exact tree; legacy fixtures remain intact.

## 6. Completion checklist

- [ ] Every required original feature and namespace/export compatibility works.
- [ ] The replacement runs independently from an original ANCOMBC installation.
- [ ] Actual candidate numerical execution is traced without changing results.
- [ ] Complete result bytes match on every required case/profile.
- [ ] Defaults, positional calls, attributes, errors and state are checked.
- [ ] No tolerance/exemption/rounding/report-only acceptance path remains.
- [ ] Both benchmark arms analyze the same inputs and options with real workers.
- [ ] Equal timing boundaries and external process-group memory are measured.
- [ ] Repeated paired campaigns pass validation and all performance gates.
- [ ] Standalone package artifacts pass exact and installation checks.
- [ ] Documentation cites the same tested artifacts and campaign.

Finite testing cannot prove every possible input. The engineering basis for a
strong replacement claim is retaining the original semantics/primitives,
restricting accelerated changes to operations with an exactness argument, and
enforcing the strict differential suite. No model should promise universal
cross-platform byte identity or arbitrary speedups without that evidence.
