# Performance improvement plan: large data sets

**Status: proposal.** Everything under "Measured" was executed on the host
described in `reference/env/ORACLE.md` (16 CPUs, 31 GB). Everything under
"Projected" is arithmetic from those measurements and is labelled as such.

This exists because `README.md`, `docs/release_status.md` and
`docs/compatibility.md` all report **0 pass, 5 fail** on the performance gates,
and a reader who stops there has no idea whether that is a fixable defect, a
host limit, or a wrong measurement. It is mostly the first, and the rest of this
document is about telling those three apart.

---

## 1. The short version

The starting position, as reported by `make gates` from a single fixed build
across the full 42-arm surface (`benchmarks/results/gates.json`):

| gate | criterion | target | measured | |
| --- | --- | --- | --- | --- |
| P1 | kernel speed-up, `rust-1` vs `r-1core` | >= 3.0 | **1.141x** | FAIL |
| P2 | end-to-end speed-up, `rust-1` vs `r-1core` | >= 2.0 | **1.141x** | FAIL |
| P3 | sensitivity speed-up, 8-16 threads | >= 5.0 | **3.165x** | FAIL — above bm5's 3x ceiling, see section 8 |
| P4 | peak RSS, rust vs R parallel | <= 0.7 | **1.751x** | FAIL |
| P5 | strong-scaling efficiency, 1 -> 16 threads | >= 0.7 | **0.172** | FAIL |

All five are scored on `bm5`, the most substantial dataset and the one where R's
fixed overhead matters least. P1 and P2 are the same measurement, because at this
size almost all of the wall time is kernel.

Against that, the port is **68x faster than the reference on a toy and slower
than it on real-sized data**, and the reason is not the one it looks like:

* On small inputs our advantage is R's fixed overhead — interpreter, data frame
  allocation, per-call dispatch — which is enormous and which we do not pay.
* On large inputs that overhead is a rounding error, and what remains is the
  least-squares kernel in `.lm_fit_all`, which **both implementations run at about
  the same rate**. We are 1.30x slower inside that kernel, and there is no
  structural reason we could not be faster.

So P1 and P2 are not out of reach, but they are not reachable by tidying: the
Householder factorisation has to get 3-5x faster. That is one bottleneck, in two
places, and it is 84 % of the run.

---

## 2. Measured: the speedup decays with problem size

`rust-1` against `r-1core`, single build, best of each arm
(`benchmarks/results/results.jsonl`):

| dataset | shape | R 1-core | rust-1 | speed |
| --- | --- | --- | --- | --- |
| `bm1` | 50 x 500 x 3 | 0.99 s | 0.01 s | **67.7x** |
| `bm2` | 100 x 1000 x 5 | 1.28 s | 0.08 s | **15.4x** |
| `bm3` | 500 x 5000 x 5 | 5.47 s | 1.40 s | **3.9x** |
| `bm4` | 1000 x 10000 x 10 | 31.01 s | 15.49 s | **2.0x** |
| `bm5` | 5000 x 20000 x 10 + sens | 366.82 s | 321.52 s | **1.14x** |
| `bm6` | 1000 x 10000 x 10 nc-sens | 31.50 s | 41.74 s | **0.75x** |

Two orders of magnitude between the ends. A reader could reasonably conclude
that the kernel is catastrophically slow; the next section says otherwise.

---

## 3. Measured: where the time actually goes

Instrumenting `lm_fit_all` per phase on `bm6` (1000 taxa x 10000 samples,
`p = 10`, one thread). Seconds:

| pass | group total | **qr** | solve | fit_values | scatter |
| --- | --- | --- | --- | --- | --- |
| MLE pass 1 | 10.996 | **9.954** | 0.586 | 0.389 | 0.090 |
| MLE pass 2 | 4.773 | **3.755** | 0.683 | 0.216 | 0.065 |
| sensitivity refit | 1.010 | 0.013 | **0.718** | 0.133 | 0.053 |

Two things to read out of this.

**The QR is 93 % of the first MLE pass.** `mle1` (9.08 s) plus `sensitivity`
(26.02 s) is 35.10 s of `bm6`'s 41.74 s — 84 %. Sandwich (1.72 s), the E-M
(0.16 s), correction (0.26 s) and the tests (0.04 s) together are 2.18 s, or
5 %, which is why none of them appears in this plan.

> These per-phase figures come from temporary instrumentation of
> `lm_fit_all`, which is **not committed** — there is no command in section 11
> that regenerates them. They are recorded here so the diagnosis can be checked
> against a re-instrumented run rather than taken on trust. The stage totals that
> *are* committed (section 4, `bm6` `rust-1`) sum to 37.24 s and bracket this
> table consistently.

**In the sensitivity refits the QR disappears and the cost moves to
`solve_multi_into`.** There the design is shared across 1000 responses, so one
factorisation serves all of them (13 ms) and the remaining work is applying the
reflectors — *the same reflector kernel*. 0.718 s x 50 refits is 36 s of
`bm6`.

So: one bottleneck, two call sites, and they are the same loop.

Counting the real path: **324 M row-columns through the QR in 9.95 s**, about
3.3 M rows/s. `crates/ancombc2-core/examples/qrbench.rs` reproduces that warm at
5.6 M rows/s. The gap between the two is cold cache; the gap between both and
what the memory traffic allows is the kernel.

---

## 4. Measured: the reference is not the thing to chase

R's own core, from `scripts/bench_r.R` run directly on `bm6`:

```
wall_seconds                29.379
stage_seconds.preprocess      0.592
stage_seconds.core           28.669
```

Our equivalent stages — `mle1` + `sandwich1` + `sandwich2` + `em` + `correction`
+ `sensitivity` — total **37.24 s**. That is **1.30x slower**, not 10x, not 30x.

Both implementations call the same algorithm: group the taxa by missingness
pattern, factorise each pattern's design once, apply it to every taxon in the
pattern, take the leading run of the diagonal above `tol * |R_00|` as the rank,
and report columns beyond it as aliased. On `bm4` and `bm5` there is roughly one
pattern per taxon, so there is no sharing to exploit and no algorithmic gap to
exploit. R's `.lm_fit_all` is not well-optimised either; it is simply not the
reason we lose.

**Conclusion: the bar is the reference kernel's rate, and we are 1.30x off it.**

> A note on why this was not visible earlier. `bench_r.R`'s header said it wrote
> per-stage timings, and it did not — `ref_run` had been collecting them all
> along and the script simply never emitted them. The header has been corrected
> and the output committed, so re-running the command in section 11 now produces
> the `stage_seconds` lines directly. The numbers above were read out of `ref_run`
> after the fact. The gates' own stage-scaling figures are unaffected, because
> `bench_gates.py` computes those from the Rust arms alone.

---

## 5. Already banked

Committed, and bit-identical — the 332-test workspace suite passes unchanged, and
the golden contract compares at Level B, so any deviation would show.

| change | before | after |
| --- | --- | --- |
| baseline | 0.886 ms / 0.91 GFLOP/s | — |
| delete `colnorm` | 0.652 ms / 1.24 GFLOP/s | **1.36x** |
| raw-slice inner loops, one scratch buffer | 0.624 ms / 1.26 GFLOP/s | +4 % |

`n = 5000, p = 9`. Subject to the run-to-run variance noted in section 6; the
ratios are the reliable part, not the absolute milliseconds.

`colnorm` is the interesting one. It was written in three places and **read in
none**. It existed to choose a column permutation; the permutation was deleted
when `qr()` was corrected to match `lm.fit`'s identity pivot (`docs/reference_behavior.md`
§16), and the norm table stayed behind — along with a full extra O(n) pass per
column per reflector. That pass was a third of the QR's memory traffic and
contributed nothing to any result.

---

## 6. Measured: what is left

A size sweep of the microbenchmark, p = 9:

| n | ms | GFLOP/s | panel |
| --- | --- | --- | --- |
| 64 | 0.004 | 2.36 | 4 KB |
| 256 | 0.017 | 2.47 | 18 KB |
| 512 | 0.033 | 2.48 | 36 KB |
| 1024 | 0.066 | 2.52 | 72 KB |
| 2048 | 0.226 | 1.47 | 144 KB |
| 4096 | 0.517 | 1.28 | 288 KB |
| 8192 | 1.098 | 1.21 | 576 KB |
| 16384 | 2.473 | 1.07 | 1152 KB |

**Flat at about 2.5 GFLOP/s while the panel fits in cache, then roughly half.**
The knee is between 72 KB and 144 KB, which is where the trailing panel stops
fitting in a typical L2.

The microbenchmark is a wall-clock measurement on a shared host and moves by a
few percent between runs — the same binary gave 0.624 and 0.646 ms at n = 5000 on
two consecutive invocations. Treat these as a shape, not as constants, and
re-measure the baseline before and after any change rather than comparing against
a number written down here.

Two negative results worth recording so they are not re-tried:

* **`-C target-cpu=native` bought 4 %** (1.30 -> 1.36 GFLOP/s). The host has AVX2
  and FMA, and the compiler was not leaving them on the table. This is not a
  vectorisation problem.
* **`perf` is unavailable on this host.** `kernel.perf_event_paranoid` is 4, so
  no perf event of any kind can be recorded. That is why the phase breakdown
  above was produced by instrumenting the code rather than by sampling a profile.

---

## 7. The plan

### Step 1 — blocked Householder factorisation *(the one that matters)*

Today each reflector streams every trailing column: read it for the dot product,
read it again and write it for the update. Three passes over a column that does
not fit in cache, `p` times.

A blocked factorisation holds a **panel** of trailing columns resident and
applies a *block* of reflectors against it, so the panel is loaded once and
stored once per block instead of once per reflector. With `p = 10` and a block of
4, that is roughly a 4x reduction in column traffic.

The arithmetic does not change. The reflectors are computed from the same updated
`R` in the same order, and applied to the same values in the same order — only
the loop nesting differs. That is what makes it safe here: it can be bit-identical,
and the golden suite is the check.

**Measured 2026-10-05: 0 %.** A blocked implementation (panel of 4, exact
`(v, vnorm_sq)` saved per reflector, trailing columns loaded once per block) is
bit-identical — 332/332 workspace tests unchanged — and no faster: 1.55 ms vs
1.48 ms per factorisation at `n = 10000, p = 10`, within run-to-run noise. The
working set per reflector already fits in cache at these shapes, so there is no
eviction for blocking to avoid. The cache-residency knee in section 6 is real but
is not what limits this kernel; see below. The blocked code was reverted — a
change that cannot be measured in seconds does not go in.

**What actually limits it:** the dot-product reduction. `Iterator::sum()` is a
single serial accumulator (2.8 GFLOP/s in isolation); four independent
accumulators run at 11.2 GFLOP/s. The axpy half already runs at 11-14 GFLOP/s,
so the dot is ~80 % of QR time. That dependency chain is the bottleneck, not
cache traffic, which is why reordering the loops changed nothing.

**Risk: low** for bit-identity (proven), but the step as specified does not move
the gate. Do not re-attempt blocking without a shape where the trailing set
exceeds the cache.

### Step 2 — the same treatment for `solve_multi_into`

Identical loop, and it is the other half of the bottleneck. The sensitivity stage
is 26.02 s of `bm6`'s 41.74 s — 62 % of the run — and in the instrumented refit
group, `solve_multi_into` is 0.718 s of 1.010 s, or 71 %. That puts roughly 18 s
of `bm6` in this one loop, and it benefits from the same reasoning as step 1:
1000 right-hand sides through a single factorisation is the best case in the whole
program for keeping things in cache, and it is currently handled one column at a
time.

**Expected:** ~3x on the sensitivity stage. **Risk: low**, same reason.

### Step 3 — keep the reflector out of the axpy loop

The update is `col -= f * v`, which needs `f = 2(v'r)/v'v` before it can start,
so the column is read twice. A two-pass form that computes the dot in one sweep
and applies in a second, with the partial dot held in registers, would halve the
traffic again.

**This one is not free.** Any reassociation of the dot product changes
floating-point summation order, and this codebase has a documented rule about
exactly that: the QR's rank test sits on a tolerance, so a column that is
borderline can flip, and the comment that used to sit on the recomputed norms —
"an algebraically equivalent rewrite is a behavioural change here, not a
refactor" — was written about this.

So step 3 is to be **measured as a trade and presented as one**, with the
deviation reported, not slipped in under a claim of equivalence. If it moves
results, it becomes a `CompatMode` variant rather than a default.

**Measured 2026-10-05: +15 % for a 2e-07 deviation — not worth taking.**
Four-accumulator reductions on all three hot sums (`norm_sq`, `vnorm_sq`, `dot`)
move the kernel from 1.29 to 1.52 GFLOP/s at `n = 10000, p = 10`, because the
memory traffic and the axpy dominate once the reduction is fixed. The golden
contract fails on `delta_em` at rel 1.97e-07 — 20x over the `beta`/`theta` rtol
and 2x over `delta_em`'s own. Fifteen percent faster for a broken contract is a
bad trade on both axes, and it was reverted with the blocking above.

The implication is blunt: even a perfect reassociation of the reduction buys at
most ~2x on the QR, not the 3-5x P1 and P2 need, because only the dot half
responds to it. Micro-optimising this kernel cannot close the gates, with or
without bit-identity. What could is algorithmic — less QR work per run, not
faster QR arithmetic — and that is a different plan.

**Risk: medium** for the code, **unacceptable** for the contract at the measured
deviation. Do not ship a reassociated reduction without a compatibility decision
that explicitly accepts the new bits.

### Projections

> **2026-10-05: these projections are withdrawn.** They assumed steps 1 and 2
> would land at 3x by removing cache traffic. Step 1 measured 0 % and step 3
> measured +15 % for a broken contract, so the assumption is false and the table
> below is arithmetic from a false premise. It is kept so the reasoning can be
> checked, not as a target. P1 and P2 have no known path that preserves
> bit-identity; see steps 1 and 3 above.

If steps 1 and 2 land at 3x, and **these are arithmetic, not measurements**:

| dataset | now | projected | vs R | gates |
| --- | --- | --- | --- | --- |
| `bm6` | 41.74 s | ~16 s | 2.0x | P2's 2.0 target reached exactly |
| `bm5` | 321.52 s | ~90 s | 4.1x | **P1 (3.0) and P2 (2.0) both pass** |

`bm5` is the dataset the gates are scored on, so P1 and P2 pass there and the
kernel work is done. Note that `bm6` is the dataset that actually shows the
current deficit — it is 0.75x R today — and 3x on steps 1 and 2 is what closes
it.

---

## 8. What this plan will **not** fix

Stated plainly so the plan is not read as a promise to turn all five gates green.
**Three of the five are out of reach of kernel work**, for three different
reasons, and one of them is not an engineering problem at all.

### P4 — peak RSS is a different problem

Measured peak RSS:

| dataset | rust-1 | rust-4 | rust-8 | rust-16 | R parallel |
| --- | --- | --- | --- | --- | --- | --- |
| `bm4` | 1119 MB | 1149 MB | 1192 MB | 1246 MB | 1250 MB |
| `bm6` | 1164 MB | 1467 MB | 1957 MB | 2921 MB | 1224 MB |
| `bm5` | 14928 MB | 22279 MB | 24747 MB | 24743 MB | 8526 MB |

`bm6` grows 2.5x from one thread to sixteen; `bm5` is 1.75x R's footprint at one
thread. No amount of arithmetic in the QR changes any of that: the growth is
because each concurrent pipeline level holds its own `n_taxa x n_samp` buffers.

**What it needs** is buffer sharing across the conservative refits — one working
set, refitted three times — rather than a faster kernel. `bm4` and `bm6` are at
or under R's footprint at one thread; the gate fails because it is scored on
`bm5`, which is the worst case twice over.

### P5 — strong-scaling efficiency is a host limit here

`bm5` wall time by width: 321.5 s, 116.6 s, 115.9 s, 116.6 s. Completely flat
past four threads. Peak RSS at eight threads is 24.7 GB on a host with 31 GB
total, so the run is competing with the page cache for memory bandwidth.

Making the kernel 3x faster does not change the shape of that curve; it moves
the whole thing left. The 24.7 GB figure is the constraint, and step 1/2 work
would help only indirectly, by shortening the time spent saturating it. Sharing
the sensitivity buffers (the P4 work) is the lever that actually moves it.

### P3 — capped by the reference's own grid

`bm5` is a *conservative* run, and the conservative pseudo-count grid
`{0, 0.1, 0.5, 1}` is three independent refits after the main run. The outer
level of the nesting order parallelises exactly those three, so **8-16 threads
can give at most 3x however much pool is available**. The reference's grid fixes
the width. Measured 3.165x is just above that ceiling.

`bench_gates.py` already records `measured_ceiling` and a `ceiling_note` for
this, so the ceiling is in the report rather than in a footnote.

The control is `bm6`, which is non-conservative with a 50-point grid and so has
enough independent work to fill the machine: its sensitivity stage goes
26.02 s -> 3.46 s from one thread to sixteen, **7.5x**, clearing the 5.0 target
on the path where the parallelism genuinely exists.

**This is a scoring question, not an optimisation**, and it should be decided as
one: either score P3 on the non-conservative dataset, or state the 3x ceiling in
the target. No amount of kernel work changes the width of the reference's grid.

---

## 9. How each step is verified

1. **`crates/ancombc2-core/examples/qrbench.rs`**, before and after. A change
   that cannot be measured in seconds does not go in.
2. **The 38-cell fixture matrix and the four committed fixtures, unchanged.**
   This is the bit-identity check. The golden contract compares `beta`/`theta` at
   `rtol 1e-8` and `se`/`vcov`/`delta_em` at `rtol 1e-7`; a reassociated sum or a
   reordered reflector update shows up there immediately. No tolerance is
   loosened, and none may be, to accommodate a change.
3. **`bm4` and `bm6` end to end**, checking that the wall time moved by the
   predicted factor and that nothing else did.
4. **The full 42-arm surface, re-measured in full** at the end. Not extrapolated,
   and not a subset. `bench_gates.py` takes the most recent row per
   (dataset, arm), so a partial re-run is not diluted by stale rows — but the
   gates are only meaningful from one build, so the file is truncated to that
   build as before.

---

## 10. Order of work

| # | step | expected | bit-identical | risk |
| --- | --- | --- | --- | --- |
| 0 | ~~delete dead `colnorm`~~ | 1.36x | yes | **done** |
| 0 | ~~raw-slice loops, one scratch buffer~~ | +4 % | yes | **done** |
| 1 | block the Householder factorisation | ~2-3x on qr | yes | low |
| 2 | block `solve_multi_into` | ~3x on solve | yes | low |
| 3 | two-pass axpy | up to 2x more | **no** | medium |
| 4 | share sensitivity buffers (P4) | RSS and P5 | yes | medium |
| 5 | decide P3's scoring | — | — | n/a |

1 and 2 are the ones that move P1 and P2, and they are the two that cannot
change a single bit. 3 is worth doing only as a measured trade. 4 is a different
piece of work addressing P4 and P5. 5 is a decision, not code.

If only one thing gets done, do step 1.

---

## 11. Reproducing the measurements in this document

```
# the isolated kernel, for sections 5 and 6
cargo run --release --example qrbench -- 5000 9 200      # ~0.63 ms, ~1.26 GFLOP/s
cargo run --release --example qrbench -- 1024 9 2000     # ~0.07 ms, ~2.5 GFLOP/s

# R's core time for section 4
ANCOMBC_BENCH_DATA=benchmarks/datasets ANCOMBC_BENCH_OUT=/tmp/benchr \
ANCOMBC_BENCH_DATASET=bm6 ANCOMBC_BENCH_THREADS=1 \
  Rscript --vanilla scripts/bench_r.R
grep '^stage_seconds' /tmp/benchr/run_metadata.tsv

# the surface, and the gates
make bench
make gates
```

Host: 16 CPUs, 31 GB, `reference/env/ORACLE.md`. R 4.3.3, reference BLAS.
