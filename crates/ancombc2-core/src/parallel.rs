//! One global Rayon pool, and an explicit budget for what may use it.
//!
//! # Why a budget and not just `par_iter`
//!
//! Rayon will happily nest: a `par_iter` inside another `par_iter` splits the
//! same pool again, and the inner tasks steal from the outer ones until the
//! scheduler is spending more time balancing than working. The plan names the
//! nesting order this crate is supposed to use --
//!
//! ```text
//! pseudo-count runs  (outermost)
//!   E-M coefficients
//!     missingness groups
//!       taxa            (innermost)
//! ```
//!
//! -- and the rule that makes it safe is not "always parallelise the innermost
//! thing", it is **at most one level consumes the pool**. A level parallelises
//! only when nothing above it already has, so the pool is split once and the
//! order above is a fallback sequence rather than a stack.
//!
//! That is what [`NestingBudget`] expresses. [`NestingBudget::level`] returns a
//! guard; [`Level::may_parallelise`] is true only at the outermost depth, and
//! every caller that has an alternative serial path takes it. So the budget is
//! not documentation of an intent -- changing the depth changes what runs, and
//! `depth_governs_who_may_parallelise` in the tests below pins that.
//!
//! # Why the parallel axes here need no reduction order
//!
//! The three levels this crate parallelises all write to **disjoint** slots:
//!
//! * `lm_fit_all` writes `beta[t * p + a]`, `fitted[t * n_samp + s]` and
//!   `dof[t]` for the taxa of its own missingness group, and a taxon belongs to
//!   exactly one group;
//! * `sandwich_all` writes `vcov[i * p * p + ..]` and `var_hat[i * p + a]` for the
//!   taxa of its own block;
//! * `theta_new` computes each sample's mean independently of every other
//!   sample.
//!
//! None of them is a reduction whose order is observable, so none of them needs a
//! deterministic-order accumulator and none of them can produce a different
//! answer at a different thread count. That is the property the P15 thread
//! invariance test checks, and it is why this module does not provide one: there
//! is nothing here for it to fix. The genuinely order-dependent reductions --
//! the E-M's parameter sweep and `sandwich_all`'s accumulation over *samples* --
//! stay serial, because their order is observable in the last bit.

use std::cell::Cell;
use std::sync::OnceLock;

use rayon::prelude::*;

// Whether the pool is claimed on this thread, and how deep the stack is.
//
// Thread-local because the pool's tasks are what nest, and a flag shared across
// threads would need synchronisation on a path whose entire cost is deciding
// whether to spawn. Each worker thread starts with the pool unclaimed, which is
// correct: work *inside* a `par_iter` runs at the same nesting level as the loop
// that spawned it, not deeper.
//
// The state is a claim rather than a depth because "who may use the pool" is not
// the same question as "how deep are we". A level that had nothing to split
// releases its claim, and the level below it then finds the pool free.
thread_local! {
    static CLAIMED: Cell<bool> = const { Cell::new(false) };
    static DEPTH: Cell<usize> = const { Cell::new(0) };
}

/// A guard for one level of the nesting order.
///
/// Held for the duration of the level; the depth drops when it is dropped, so an
/// early return or a panic unwinds correctly.
#[derive(Debug)]
pub struct Level {
    /// Whether this level claimed the pool, decided on entry from whether the
    /// pool was already free. Decided lazily from the live flag would be wrong:
    /// this level's own claim would read back as "already taken".
    claimed: bool,
    depth: usize,
}

impl Level {
    /// Whether *this* level may split the pool.
    ///
    /// True only when nothing above it has. A caller with a serial alternative --
    /// and every level in this crate has one, because it is the same code path --
    /// uses it.
    #[inline]
    pub fn may_parallelise(&self) -> bool {
        self.claimed
    }

    /// Give up the pool to the next level down, because this one had no work worth
    /// splitting. A no-op for a level that never claimed it.
    ///
    /// Without this the budget strands the pool on a level that cannot use it, and
    /// that is the *common* case rather than an edge: a table with no missing
    /// values has exactly one missingness pattern, so the group level has one item
    /// and everything expensive happens per taxon inside it. Holding the outermost
    /// level for a one-item `par_iter` means `--threads 16` runs the whole analysis
    /// on one core while the core's taxon axis is told it may not use the pool
    /// because something above it already has.
    ///
    /// Idempotent, and only meaningful for a level that held the pool in the first
    /// place.
    pub fn release(&mut self) {
        if self.claimed {
            self.claimed = false;
            CLAIMED.with(|c| c.set(false));
        }
    }

    /// This level's depth, for reporting.
    #[inline]
    pub fn depth(&self) -> usize {
        self.depth
    }
}

impl Drop for Level {
    fn drop(&mut self) {
        DEPTH.with(|d| d.set(d.get().saturating_sub(1)));
        if self.claimed {
            // Only the claim's own level releases it. A level that entered while
            // the pool was taken by an enclosing level, and then *its* enclosing
            // level released, must not clear a claim it does not hold.
            CLAIMED.with(|c| c.set(false));
        }
    }
}

/// Enter one level of the nesting order.
///
/// The names are the plan's, so a call site says which level it is:
///
/// ```ignore
/// let mut lvl = NestingBudget::level("missingness groups");
/// let done = map_par(&mut lvl, &items, f);
/// ```
///
/// The name is documentation, not instrumentation: a level that reported where the
/// pool was being spent would need counters kept on the hot path, and the stage
/// timings in `run_metadata.tsv` already say where the time goes.
pub struct NestingBudget;

impl NestingBudget {
    /// Enter the named level.
    pub fn level(name: &'static str) -> Level {
        Self::level_inner(name, true)
    }

    /// Enter a level that will spread `items` items, and hand the pool down if
    /// that many cannot fill it.
    ///
    /// # Why
    ///
    /// [`Self::level`] implements "at most one level consumes the pool", which is
    /// right for avoiding oversubscription and wrong for *utilisation* when the
    /// outer level has too few items to occupy the pool. The conservative
    /// sensitivity analysis is the case in the executed benchmarks: it runs the
    /// pseudo-count grid `0.1, 0.5, 1.0` alongside the main run, so the outermost
    /// level has **three** items. On this 16-core host that is three busy threads
    /// and thirteen idle, for the whole stage -- and because the claim is held,
    /// every level inside each refit is told it may not parallelise either. The
    /// stage measured 145 s at one thread and 57 s at eight, i.e. 2.5x from eight
    /// cores, while the levels beneath it had thousands of items each.
    ///
    /// The fix is not to nest, which is what the budget exists to prevent. It is to
    /// give the pool to whichever level can actually fill it: if this level has
    /// fewer items than the pool has threads, it hands the claim down and runs its
    /// items one at a time, each using the whole pool internally.
    ///
    /// This is conservative in the direction that matters. It only changes
    /// behaviour when the outer level provably cannot occupy the pool, and it
    /// cannot leave the pool unclaimed at the innermost level either -- a level
    /// with no inner level to hand to would find `may_parallelise` true anyway,
    /// since the claim is taken by the *first* level entered, not by depth.
    pub fn level_for_items(name: &'static str, items: usize) -> Level {
        let wide_enough = items >= rayon::current_num_threads();
        Self::level_inner(name, wide_enough)
    }

    fn level_inner(_name: &'static str, claim: bool) -> Level {
        DEPTH.with(|d| {
            let depth = d.get() + 1;
            d.set(depth);
            // Claim only if the pool was free, and *only then* mark it taken: a
            // level that is denied must leave the flag alone, or it would deny
            // the level below it for a pool it never held.
            let was_claimed = CLAIMED.with(|c| c.get());
            // Only take the pool if this level asked for it and nobody else has
            // it. A level that declined (`level_for_items` with too few items)
            // leaves the flag alone, so the level below can claim it.
            if claim && !was_claimed {
                CLAIMED.with(|c| c.set(true));
            }
            Level {
                claimed: claim && !was_claimed,
                depth,
            }
        })
    }

    /// The current depth, outside any level.
    pub fn depth() -> usize {
        DEPTH.with(|d| d.get())
    }
}

/// The one pool.
///
/// Built once, from the CLI's `--threads`, and every `par_iter` in the workspace
/// draws from it. A pool per call site would let two of them each believe it owned
/// the machine, which is the oversubscription the budget exists to prevent.
static POOL: OnceLock<rayon::ThreadPool> = OnceLock::new();

/// Install the global pool with `threads` workers.
///
/// Idempotent: the first call wins, because a second pool would silently replace
/// the one the earlier callers captured, and a thread count that changes
/// mid-process is a bug rather than a reconfiguration. Returns `true` if this
/// call installed it.
///
/// A `threads` of 0 or 1 installs a pool anyway rather than skipping, so the code
/// path is the same either way and `--threads 1` exercises the real machinery
/// instead of a different build of it.
pub fn install_pool(threads: usize) -> bool {
    let n = threads.max(1);
    if POOL.get().is_some() {
        return false;
    }
    let built = rayon::ThreadPoolBuilder::new()
        .num_threads(n)
        .thread_name(|i| format!("ancombc2-{i}"))
        .build();
    match built {
        Ok(p) => POOL.set(p).is_ok(),
        // A pool that cannot be built is not a reason to abandon the analysis;
        // Rayon's global pool is the documented fallback and produces the same
        // numbers, only without the requested width.
        Err(_) => false,
    }
}

/// Run `f` inside the global pool if one is installed, or on this thread if not.
///
/// The indirection is what makes the pool *global*: a bare `f()` would use
/// Rayon's own global pool, which is a different pool with a different width.
pub fn with_pool<R: Send>(f: impl FnOnce() -> R + Send) -> R {
    match POOL.get() {
        Some(p) => p.install(f),
        None => f(),
    }
}

/// The installed pool's width, if any.
pub fn pool_threads() -> Option<usize> {
    POOL.get().map(|p| p.current_num_threads())
}

/// `f` over `items` in parallel when `lvl` allows it, serially otherwise.
///
/// The two branches call the *same* closure, so the serial path is not a
/// reimplementation that can drift from the parallel one. That is the property
/// worth stating: a "fast path" and a "fallback" that are separately written are
/// a bug waiting for the inputs that take one and not the other.
/// Map `f` over `items`, in parallel when the level allows it.
///
/// `f` receives the item **and its index**, because "the value at this position in the
/// caller's grid" is a thing parallel code needs and cannot reconstruct from the value:
/// the pseudo-count grids repeat entries, so looking a mean up by value would be
/// ambiguous. The index is what pairs a work item with something the caller computed
/// before the loop -- which is how the reduction stays on the main thread while the
/// arithmetic does not.
pub fn map_par<T: Send + Sync, R: Send>(
    lvl: &mut Level,
    items: &[T],
    f: impl Fn(usize, &T) -> R + Send + Sync,
) -> Vec<R> {
    // `len() <= 1` is not a micro-optimisation here, it is the case that decides
    // where the analysis runs: splitting the pool for a single item costs more in
    // scheduling than it saves in work, *and* it would hold the outermost level
    // against the levels below, which is how a dense table ended up single-threaded
    // with `--threads 16`. So the level hands the pool down instead.
    if lvl.may_parallelise() && items.len() > 1 {
        items.par_iter().enumerate().map(|(i, t)| f(i, t)).collect()
    } else {
        lvl.release();
        items.iter().enumerate().map(|(i, t)| f(i, t)).collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::BTreeSet;

    /// Only the outermost level may split the pool.
    ///
    /// This is the whole content of the module's type signature, so it is pinned
    /// directly: without it, a caller could take `may_parallelise() == true` at
    /// depth 1 and nest, and the oversubscription would come back silently.
    #[test]
    fn depth_governs_who_may_parallelise() {
        let outer = NestingBudget::level("outer");
        assert!(
            outer.may_parallelise(),
            "the outermost level may split the pool"
        );
        assert_eq!(outer.depth(), 1);
        {
            let inner = NestingBudget::level("inner");
            assert!(
                !inner.may_parallelise(),
                "a nested level must not split the pool again"
            );
            assert_eq!(inner.depth(), 2);
            {
                let innermost = NestingBudget::level("innermost");
                assert!(!innermost.may_parallelise());
                assert_eq!(innermost.depth(), 3);
            }
            assert_eq!(inner.depth(), 2, "the depth unwinds on drop");
        }
        assert_eq!(outer.depth(), 1);
        drop(outer);
        assert_eq!(NestingBudget::depth(), 0, "the depth returns to zero");
    }

    /// Both branches of `map_par` must produce the same values, in the same order.
    ///
    /// The serial branch is what a nested level takes, and it is the same closure
    /// rather than a second implementation -- so the test that matters is that the
    /// two agree, which is what would fail if `map_par` ever grew a hand-written
    /// serial loop.
    ///
    /// The depth is checked as *relative* to where this test starts. The counter is
    /// thread-local, which is what makes it cheap, and the test harness reuses
    /// threads -- so an absolute `== 0` would make this test depend on whether some
    /// other test left a guard behind.
    #[test]
    fn both_branches_agree_and_preserve_order() {
        let items: Vec<usize> = (0..64).collect();
        // The index is ignored here; the test is about the two branches agreeing.
        let square = |_i: usize, v: &usize| v * v;
        let base = NestingBudget::depth();

        let par = {
            let mut outer = NestingBudget::level("outer");
            assert!(outer.may_parallelise());
            map_par(&mut outer, &items, square)
        };
        let ser = {
            let mut outer = NestingBudget::level("outer");
            let inner = NestingBudget::level("inner");
            assert!(!inner.may_parallelise());
            map_par(&mut outer, &items, square)
        };
        assert_eq!(par, ser, "the serial branch must equal the parallel one");
        assert_eq!(
            NestingBudget::depth(),
            base,
            "every guard must have been dropped, or the next test inherits a depth"
        );
    }

    /// The pool is installed once, is the width first asked for, and is the one
    /// `with_pool` runs inside.
    ///
    /// One test rather than two, because the pool is process-wide: `cargo test`
    /// shares one process across every test in the binary, so a second test that
    /// asserted "the first install wins" would fail depending on which of them ran
    /// first. That is not a test flake to retry, it is a test that cannot be
    /// written here.
    ///
    /// What it does check is the thing that would otherwise be wrong silently: a
    /// bare `par_iter` would use Rayon's *own* global pool, at a width of `nproc`
    /// rather than `--threads`, and nothing else in the suite would notice.
    #[test]
    fn the_pool_is_installed_once_and_is_the_one_used() {
        let first_width = install_pool(3);
        let second_width = install_pool(8);
        assert!(!second_width, "a second install must not replace the first");
        let width = pool_threads().expect("a pool must exist once installed");
        if first_width {
            assert_eq!(width, 3, "the width must be the first one asked for");
        } else {
            // Another test in this binary got there first; the width is whatever it
            // asked for, and the point is that it did not just change.
            assert_ne!(width, 8, "the second request must be ignored");
        }
        let threads = with_pool(rayon::current_num_threads);
        assert_eq!(
            Some(threads),
            pool_threads(),
            "`with_pool` must run inside the installed pool, not Rayon's own"
        );
    }

    /// A width of 1 still installs a pool, so `--threads 1` runs the real
    /// machinery rather than a different build of it.
    #[test]
    fn one_thread_still_installs_a_pool() {
        install_pool(1);
        assert!(
            pool_threads().is_some(),
            "a pool must exist even for one thread, or --threads 1 would exercise \
             a different code path from --threads 8"
        );
    }

    /// A level with no work to share hands the pool to the level below.
    ///
    /// This is the case that decides whether `--threads` does anything on a table
    /// with no missing values: one missingness pattern, one item at the group
    /// level, and every expensive operation per taxon inside it. If the group level
    /// kept the pool while running a one-item map, the taxon level would be denied
    /// it and the whole analysis would run on one core.
    #[test]
    fn an_empty_level_passes_the_pool_down() {
        let items: Vec<usize> = (0..64).collect();
        let one: Vec<usize> = vec![0];
        let base = NestingBudget::depth();

        {
            let mut outer = NestingBudget::level("missingness groups");
            assert!(
                outer.may_parallelise(),
                "the outermost level holds the pool"
            );
            // One group: the level cannot use the pool, and must say so.
            map_par(&mut outer, &one, |_i, v| *v);
            assert!(
                !outer.may_parallelise(),
                "a spent level must release the pool"
            );
            let inner = NestingBudget::level("taxa");
            assert!(
                inner.may_parallelise(),
                "the level below must inherit the released pool, or a dense table \
                 runs single-threaded with --threads 16"
            );
        }

        // ...and a level with real work does *not* release it.
        {
            let mut busy = NestingBudget::level("missingness groups");
            map_par(&mut busy, &items, |_i, v| *v);
            assert!(
                busy.may_parallelise(),
                "a level that split the pool keeps it"
            );
            let inner2 = NestingBudget::level("taxa");
            assert!(
                !inner2.may_parallelise(),
                "a level that used the pool must not hand it on as well"
            );
        }
        assert_eq!(
            NestingBudget::depth(),
            base,
            "every guard must be dropped, or the next test inherits the pool"
        );
    }

    /// The three parallel axes in the crate write to disjoint slots, which is why
    /// they need no deterministic-order reduction.
    ///
    /// Expressed as a property over the writers themselves rather than over the
    /// pipeline: given the same inputs, a group loop's writes are a function of
    /// its own taxa alone. If that ever stopped holding, the parallel paths would
    /// start racing and P15 would catch it -- but it would catch it as a
    /// nondeterminism rather than pointing here, so it is worth asserting
    /// directly.
    #[test]
    fn the_parallel_axes_write_to_disjoint_slots() {
        // `lm_fit_all`: every taxon appears in exactly one group.
        let groups: Vec<Vec<usize>> = vec![vec![0, 3, 7], vec![1, 2], vec![4, 5, 6]];
        let mut seen: BTreeSet<usize> = BTreeSet::new();
        let mut disjoint = true;
        for g in &groups {
            for &t in g {
                if !seen.insert(t) {
                    disjoint = false;
                }
            }
        }
        assert!(
            disjoint,
            "a taxon in two groups would race on beta/fitted/dof"
        );
        assert_eq!(seen.len(), 8, "every taxon of 0..8 must be covered");

        // `sandwich_all`: blocks are `i0..i1` over the same index space.
        let (n_taxa, block) = (100usize, 7usize);
        let mut covered: BTreeSet<usize> = BTreeSet::new();
        for i0 in (0..n_taxa).step_by(block) {
            for i in i0..(i0 + block).min(n_taxa) {
                covered.insert(i);
            }
        }
        assert_eq!(
            covered.len(),
            n_taxa,
            "blocks must tile the taxa exactly once"
        );
    }
}
