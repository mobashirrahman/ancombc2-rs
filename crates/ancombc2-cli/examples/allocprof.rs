//! Attribute a run's allocation traffic by size.
//!
//! The benchmark contract reports `bytes_allocated` and `allocation count` as
//! totals, which is enough to gate on but not enough to act on: a total says
//! *that* a run allocates a lot, not *where*. This example installs a counting
//! global allocator and reports a log2(size) histogram, which is what turned
//! "bm4 allocates 48 million times" into two specific defects -- the counts
//! reader building a `String` per cell, and the missingness-group designs
//! building a row label per row per group -- each of which was one allocation
//! per cell of the table.
//!
//! ```text
//! cargo run --release -p ancombc2-cli --example allocprof -- benchmarks/datasets/bm4
//! ```
//!
//! Two caveats, so the numbers are not over-read:
//!
//! * The per-bucket `MB` column is `count * 2^(bucket-1)`, the *lower* bound of
//!   the bucket, so it systematically under-counts. The `bytes` total is the
//!   authoritative figure; the histogram is for locating, not for accounting.
//! * The counts are this process's, including reading the input and writing
//!   nothing. Reading is a real part of a benchmark arm, so it is left in.

use ancombc2_core::AncombcConfig;
use std::alloc::{GlobalAlloc, Layout, System};
use std::sync::atomic::{AtomicI64, AtomicU64, Ordering::Relaxed};

/// Note a live allocation of `size`, and update the peak live total.
#[inline]
fn charge_live(size: usize) {
    let b = (usize::BITS - size.leading_zeros()) as usize;
    let b = b.min(31);
    if b as u64 == BT_BUCKET.load(Ordering::Relaxed) && BT_LEFT.load(Ordering::Relaxed) > 0 {
        let ok = IN_HOOK.with(|f| {
            if f.get() {
                false
            } else {
                f.set(true);
                true
            }
        });
        if ok {
            BT_LEFT.fetch_sub(1, Ordering::Relaxed);
            if let Ok(mut v) = BT_DONE.lock() {
                v.push(format!(
                    "--- size {size} ---\n{}",
                    std::backtrace::Backtrace::force_capture()
                ));
            }
            IN_HOOK.with(|f| f.set(false));
        }
    }
    let n = size as i64;
    LIVE[b].fetch_add(n, Ordering::Relaxed);
    // The peak has to be the *sum* over buckets. Tracking the max of any single
    // bucket reports one structure's size, not the process's live set, which
    // understates it by however many large structures are live at once.
    let now = LIVE_TOTAL.fetch_add(n, Ordering::Relaxed) + n;
    // Claim the peak first, then snapshot: only the thread that actually set the
    // new maximum records, and the snapshot is consistent because the loser does
    // not touch it.
    if PEAK.load(Ordering::Relaxed) < now {
        PEAK.store(now, Ordering::Relaxed);
        for (slot, live) in PEAK_SNAP.iter().zip(LIVE.iter()) {
            slot.store(live.load(Ordering::Relaxed), Ordering::Relaxed);
        }
        // Capture on *new peak*, not on the first allocation of this size: the
        // point is to attribute the allocation that made the peak, which is by
        // definition not one of the first three.
        if b as u64 == BT_BUCKET.load(Ordering::Relaxed) && now >= BT_MIN.load(Ordering::Relaxed) {
            let ok = IN_HOOK.with(|f| {
                if f.get() {
                    false
                } else {
                    f.set(true);
                    true
                }
            });
            if ok {
                BT_LEFT.fetch_sub(1, Ordering::Relaxed);
                if let Ok(mut v) = BT_DONE.lock() {
                    v.push(format!(
                        "--- new peak {:.3} GB by a {size}-byte allocation ---\n{}",
                        now as f64 / 1e9,
                        std::backtrace::Backtrace::force_capture()
                    ));
                }
                IN_HOOK.with(|f| f.set(false));
            }
        }
    }
}

#[inline]
fn release_live(size: usize) {
    let b = ((usize::BITS - size.leading_zeros()) as usize).min(31);
    LIVE[b].fetch_sub(size as i64, Ordering::Relaxed);
    LIVE_TOTAL.fetch_sub(size as i64, Ordering::Relaxed);
}

/// Live bytes per log2 bucket, for a snapshot.
fn live_histogram() -> Vec<(usize, i64)> {
    (0..32)
        .map(|b| (b, LIVE[b].load(Ordering::Relaxed)))
        .filter(|(_, v)| *v != 0)
        .collect()
}

use std::sync::atomic::Ordering;

static BYTES: AtomicU64 = AtomicU64::new(0);
static COUNT: AtomicU64 = AtomicU64::new(0);
/// Allocation count by `floor(log2(size))`, so two neighbouring buckets never
/// differ by more than a factor of two.
static HIST: [AtomicU64; 32] = [const { AtomicU64::new(0) }; 32];
/// Live bytes, by log2 size bucket, so *peak* memory can be attributed to a size
/// class rather than only total traffic.
///
/// Traffic says "bm5 asks for 40 GB", which does not say which of the ten-odd
/// multi-hundred-megabyte structures is responsible. Peak says "the live set has
/// 12 GiB in the 2^30 bucket", which narrows it immediately -- and a large live
/// footprint that total traffic never shows is exactly the P4 failure mode, since
/// P4 is about resident set size.
static LIVE: [AtomicI64; 32] = [const { AtomicI64::new(0) }; 32];
/// Live bytes across *all* buckets, so the peak is a real total.
static LIVE_TOTAL: AtomicI64 = AtomicI64::new(0);
/// The largest live total seen at any single allocation boundary.
static PEAK: AtomicI64 = AtomicI64::new(0);
/// The per-bucket composition at the moment the peak was reached, so the peak
/// can be attributed to sizes rather than only totalled. Recorded only when the
/// total is a new maximum, which is rare enough that the copy is free.
static PEAK_SNAP: [AtomicI64; 32] = [const { AtomicI64::new(0) }; 32];

/// Backtraces captured for allocations in one log2 bucket, for attributing a
/// size class to a call site.
///
/// Opt-in via `ALLOCPROF_BT=<bucket>` and capped, because capturing inside the
/// allocator will itself allocate: the guard below makes the hook reentrant-safe
/// by declining to record while a capture is in progress, and the cap bounds the
/// work. An earlier version of this tool recursed instead and overflowed the
/// stack, which is why this is a flag and not a default.
static BT_BUCKET: AtomicU64 = AtomicU64::new(u64::MAX);
static BT_LEFT: AtomicU64 = AtomicU64::new(0);
/// Only capture once the live total passes this, so the records are the
/// allocations that build the *final* peak rather than the first ones of a size.
static BT_MIN: AtomicI64 = AtomicI64::new(0);
static BT_DONE: std::sync::Mutex<Vec<String>> = std::sync::Mutex::new(Vec::new());

thread_local! {
    /// True while a capture is in flight on this thread, so the allocations that
    /// the capture itself makes are not recorded.
    static IN_HOOK: std::cell::Cell<bool> = const { std::cell::Cell::new(false) };
}

struct Prof;

unsafe impl GlobalAlloc for Prof {
    unsafe fn alloc(&self, l: Layout) -> *mut u8 {
        BYTES.fetch_add(l.size() as u64, Relaxed);
        COUNT.fetch_add(1, Relaxed);
        let b = (usize::BITS - l.size().leading_zeros()) as usize;
        HIST[b.min(31)].fetch_add(1, Relaxed);
        charge_live(l.size());
        System.alloc(l)
    }

    unsafe fn alloc_zeroed(&self, l: Layout) -> *mut u8 {
        // `vec![0.0; n]` goes through here, not `alloc`, and in this crate that is
        // the dominant way a large buffer is created -- the `n_taxa x n_samp`
        // response, the observed mask, the per-group solved block. Hooking only
        // `alloc` misses every one of them, which is how an earlier version of
        // this profiler reported a peak live set well below the process RSS.
        BYTES.fetch_add(l.size() as u64, Relaxed);
        COUNT.fetch_add(1, Relaxed);
        let b = (usize::BITS - l.size().leading_zeros()) as usize;
        HIST[b.min(31)].fetch_add(1, Relaxed);
        charge_live(l.size());
        System.alloc_zeroed(l)
    }

    unsafe fn dealloc(&self, p: *mut u8, l: Layout) {
        release_live(l.size());
        System.dealloc(p, l)
    }

    unsafe fn realloc(&self, p: *mut u8, l: Layout, n: usize) -> *mut u8 {
        // Charged at the new size: a `realloc` is the allocation growing, and the
        // `bytes` total is meant to be "how much memory the run asked for".
        BYTES.fetch_add(n as u64, Relaxed);
        COUNT.fetch_add(1, Relaxed);
        let b = (usize::BITS - n.leading_zeros()) as usize;
        HIST[b.min(31)].fetch_add(1, Relaxed);
        release_live(l.size());
        charge_live(n);
        System.realloc(p, l, n)
    }
}

#[global_allocator]
static ALLOC: Prof = Prof;

fn main() {
    let d = std::env::args()
        .nth(1)
        .expect("usage: allocprof DATASET_PREFIX  (e.g. benchmarks/datasets/bm4)");
    if let Ok(v) = std::env::var("ALLOCPROF_BT") {
        if let Ok(b) = v.parse::<u64>() {
            BT_BUCKET.store(b, Ordering::Relaxed);
            BT_LEFT.store(3, Ordering::Relaxed);
        }
        if let Ok(v) = std::env::var("ALLOCPROF_BT_MIN_GB") {
            if let Ok(g) = v.parse::<f64>() {
                BT_MIN.store((g * 1e9) as i64, Ordering::Relaxed);
            }
        }
    }
    use std::path::PathBuf;
    let counts =
        ancombc2_io::read_counts(&PathBuf::from(format!("{d}.counts.tsv"))).expect("counts");
    let meta =
        ancombc2_io::read_metadata(&PathBuf::from(format!("{d}.meta.tsv"))).expect("metadata");
    let formula = std::fs::read_to_string(format!("{d}.formula.txt")).expect("formula");
    let f = ancombc2_io::formula::parse(formula.trim()).expect("formula");
    let design = ancombc2_io::build_design(&meta, &f, Some("group")).expect("design");

    let cfg = AncombcConfig {
        fix_eff: design.colnames.clone(),
        group: Some("group".to_string()),
        group_labels: design.group.clone(),
        ..Default::default()
    };
    let run = ancombc2_io::run(&counts, &design, &cfg).expect("run");

    let cells = counts.n_taxa() as f64 * counts.n_samp() as f64;
    println!("taxa retained: {}", run.core.taxa.len());
    println!("bytes         {:.3} GB", BYTES.load(Relaxed) as f64 / 1e9);
    println!("allocations   {}", COUNT.load(Relaxed));
    println!(
        "per input cell {:.2}",
        COUNT.load(Relaxed) as f64 / cells.max(1.0)
    );
    println!("peak live total {:.3} GB", PEAK.load(Relaxed) as f64 / 1e9);
    println!("log2(size)     count         MB");
    for (b, c) in HIST.iter().enumerate() {
        let n = c.load(Relaxed);
        if n > 0 {
            println!(
                "  {b:>2}  {n:>12}  {:>8.1}",
                (n as f64 * 2f64.powi(b as i32 - 1)) / 1e6
            );
        }
    }
    if let Ok(v) = BT_DONE.lock() {
        if !v.is_empty() {
            println!(
                "\n=== backtraces for log2 bucket {} ===",
                BT_BUCKET.load(Ordering::Relaxed)
            );
            for t in v.iter() {
                println!("{t}");
            }
        }
    }
    println!();
    println!("live bytes by size bucket AT PEAK (the P4 number, by size class):");
    for (b, slot) in PEAK_SNAP.iter().enumerate() {
        let v = slot.load(Ordering::Relaxed);
        if v != 0 {
            println!(
                "  {b:>2}  {:>8.3} GB live   (~{:.0} allocs of 2^{})",
                v as f64 / 1e9,
                v as f64 / 2f64.powi(b as i32),
                b
            );
        }
    }
    println!();
    println!("live bytes by size bucket at exit (a large *live* figure here is");
    println!("what P4 measures; a large count above with a small figure is churn):");
    for (b, v) in live_histogram() {
        println!("  {b:>2}  {:>8.3} GB live", v as f64 / 1e9);
    }
}
