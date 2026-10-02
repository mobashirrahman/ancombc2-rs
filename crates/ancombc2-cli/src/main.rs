//! `ancombc2-rs`: the command-line front end.
//!
//! One run, one output directory, no hidden state. Every option that changes a
//! number has a flag, and every flag's default is the reference's default, so
//! `ancombc2-rs` with no options beyond the inputs reproduces
//! `ancombc2(..., fix_formula = <formula>)` with the reference's own defaults.
//!
//! Threading: the algorithm parallelises over the pseudo-count runs, the E-M
//! coefficients and the missingness groups, in that nesting order. `--threads`
//! sizes the single global pool; the reduction order inside each parallel loop is
//! fixed, so the result does not depend on the thread count. Property test P15
//! asserts that.

use std::alloc::{GlobalAlloc, Layout, System};
use std::path::{Path, PathBuf};
use std::process::ExitCode;
use std::sync::atomic::{AtomicU64, Ordering};
use std::time::Instant;

use ancombc2_core::config::{AncombcConfig, CompatMode};
use ancombc2_io as io;
use ancombc2_io::formula;

/// A `GlobalAlloc` that counts bytes and allocation calls.
///
/// The plan asks for "bytes allocated, allocation count" as benchmark metrics, and
/// they are the only figures that say whether a speed-up came from doing less
/// work or from doing the same work more cleverly. Two relaxed atomics on the
/// allocation path are the price, and it is paid only in the binary -- the
/// library crates used by the parity tests do not carry it.
struct Counting;

static BYTES: AtomicU64 = AtomicU64::new(0);
static ALLOCS: AtomicU64 = AtomicU64::new(0);

unsafe impl GlobalAlloc for Counting {
    unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
        BYTES.fetch_add(layout.size() as u64, Ordering::Relaxed);
        ALLOCS.fetch_add(1, Ordering::Relaxed);
        // SAFETY: forwarding to the system allocator with the caller's layout.
        unsafe { System.alloc(layout) }
    }

    unsafe fn dealloc(&self, ptr: *mut u8, layout: Layout) {
        // SAFETY: forwarding to the system allocator with the caller's layout.
        unsafe { System.dealloc(ptr, layout) }
    }

    unsafe fn realloc(&self, ptr: *mut u8, layout: Layout, new_size: usize) -> *mut u8 {
        BYTES.fetch_add(
            new_size.saturating_sub(layout.size()) as u64,
            Ordering::Relaxed,
        );
        ALLOCS.fetch_add(1, Ordering::Relaxed);
        // SAFETY: forwarding to the system allocator with both layouts.
        unsafe { System.realloc(ptr, layout, new_size) }
    }

    unsafe fn alloc_zeroed(&self, layout: Layout) -> *mut u8 {
        BYTES.fetch_add(layout.size() as u64, Ordering::Relaxed);
        ALLOCS.fetch_add(1, Ordering::Relaxed);
        // SAFETY: forwarding to the system allocator with the caller's layout.
        unsafe { System.alloc_zeroed(layout) }
    }
}

#[global_allocator]
static ALLOC: Counting = Counting;

const USAGE: &str = "\
ancombc2-rs -- ANCOM-BC2, a Rust reimplementation of ANCOMBC 2.15.2

USAGE:
    ancombc2-rs --counts <file> --meta <file> --formula <expr> [OPTIONS]

INPUTS
    --counts <file>        counts matrix, TSV or CSV, taxa as rows and samples
                           as columns, with an ID column in the corner
    --meta <file>          sample metadata, TSV or CSV, samples as rows
    --formula <expr>       fixed-effect formula, e.g. '~ group + x1 + x2'
    --group <name>         metadata column to use as the group; required for
                           --global, --pairwise and structural zeros

ALGORITHM
    --p-adj-method <m>     holm | hochberg | hommel | bonferroni | BH | BY | none
                           [default: holm]
    --pseudo <f>           pseudo-count added before the log [default: 0]
    --prv-cut <f>          prevalence cutoff [default: 0.1]
    --lib-cut <f>          library-size cutoff [default: 0]
    --s0-perc <f>          quantile of the SE distribution used for s0
                           [default: 0.05]
    --alpha <f>            significance level [default: 0.05]
    --global               run the global quadratic Wald test
    --pairwise             run the pairwise mixed-directional-FDR test
    --mdfdr-fwer <m>       family-wise method for the pairwise test [default: holm]
    --struc-zero           detect structural zeros
    --sparse-taxa          carry the count table as a compressed sparse-taxon
                           table through the prevalence and library-size screens.
                           Off by default: it only shrinks storage above ~33%
                           zeros, and only the screens benefit. See
                           docs/compatibility.md
    --neg-lb               classify structural zeros by the lower bound
    --iter-tol <f>         MLE convergence tolerance [default: 0.01]
    --iter-max <n>         MLE iteration cap [default: 20]
    --em-tol <f>           E-M convergence tolerance [default: 1e-5]
    --em-max <n>           E-M iteration cap [default: 100]
    --pseudo-sens          run the pseudo-count sensitivity analysis
    --conservative         rerun everything per pseudo-count; without it,
                           inference is refitted 50 times with pseudo 0.01..0.5
    --compat <mode>        ancombc2-2.15 (default) or strict

OUTPUT
    --out <dir>            output directory [default: .]
    --quiet                do not print the progress summary
    --threads <n>          size of the global worker pool [default: all cores]
    --version              print the compatibility target and exit
    --help                 print this message
";

fn main() -> ExitCode {
    match real_main() {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("ancombc2-rs: {e}");
            ExitCode::FAILURE
        }
    }
}

/// The compatibility target, printed by `--version` and recorded in the run
/// metadata so a result table can be traced back to what produced it.
const TARGET: &str = concat!(
    "ancombc2-rs v0.1 equivalent to ANCOMBC 2.15.2 at ",
    "dc4febdf59badb3a8dfe0c767ef2186323c2199a, R 4.5.x, seed 42"
);

fn real_main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let opts = match Options::parse(&args)? {
        Parsed::Help => {
            print!("{USAGE}");
            return Ok(());
        }
        Parsed::Version => {
            println!("{TARGET}");
            return Ok(());
        }
        Parsed::Run(o) => *o,
    };

    // One global pool, sized once. Every parallel loop in the core draws from
    // it, and the nesting budget (pseudo-count outermost, then the E-M
    // coefficients, then the missingness groups, then taxa) keeps the inner
    // loops from oversubscribing.
    let threads = opts.threads.unwrap_or_else(|| {
        std::thread::available_parallelism()
            .map(|n| n.get())
            .unwrap_or(1)
    });
    build_pool(threads);
    let t_read = Instant::now();
    let counts = io::read_counts(&opts.counts)?;
    let meta = io::read_metadata(&opts.meta)?;
    let f = formula::parse(&opts.formula)?;
    let design = io::build_design(&meta, &f, opts.group.as_deref())?;
    if !opts.quiet {
        eprintln!(
            "ancombc2-rs: read {} taxa x {} samples, design {} x {} ({:.2}s)",
            counts.n_taxa(),
            counts.n_samp(),
            design.matrix.rows,
            design.matrix.cols,
            t_read.elapsed().as_secs_f64()
        );
    }

    let mut cfg = AncombcConfig {
        fix_eff: design.colnames.clone(),
        p_adj_method: io::parse_adjust(&opts.p_adj_method)?,
        pseudo: opts.pseudo,
        pseudo_sens: opts.pseudo_sens,
        conservative: opts.conservative,
        prv_cut: opts.prv_cut,
        lib_cut: opts.lib_cut,
        s0_perc: opts.s0_perc,
        group: opts.group.clone(),
        group_labels: design.group.clone(),
        struc_zero: opts.struc_zero,
        representation: opts.representation,
        neg_lb: opts.neg_lb,
        alpha: opts.alpha,
        global: opts.global,
        pairwise: opts.pairwise,
        compat: opts.compat,
        ..Default::default()
    };
    cfg.iter_control.tol = opts.iter_tol;
    cfg.iter_control.max_iter = opts.iter_max;
    cfg.em_control.tol = opts.em_tol;
    cfg.em_control.max_iter = opts.em_max;
    cfg.mdfdr_control.fwer_ctrl_method = io::parse_adjust(&opts.mdfdr_fwer)?;

    // The sensitivity analysis writes a second set of tables, and both sets get
    // the `passed_ss`/`diff_robust` columns, so the flag is required to run it.
    if cfg.pseudo_sens && !cfg.conservative {
        // 50 refits is the reference's grid; nothing to configure.
    }

    // Inside the pool, so the `par_iter`s in the core draw from *our* pool rather
    // than from Rayon's own global one. `with_pool` is what makes the width
    // `--threads` asked for, rather than the machine's CPU count.
    let t_run = Instant::now();
    let result = ancombc2_core::parallel::with_pool(|| io::run(&counts, &design, &cfg))?;
    let elapsed = t_run.elapsed().as_secs_f64();

    let out = &opts.out;
    let primary = io::primary_table(&result);
    primary.write_to(&out.join("res.tsv"))?;
    if let Some(t) = io::global_table(&result.core) {
        t.write_to(&out.join("res_global.tsv"))?;
    }
    if let Some(t) = io::pairwise_table(&result.core) {
        t.write_to(&out.join("res_pair.tsv"))?;
    }
    if let (Some(sc), Some(passed)) = (&result.sensitivity, &result.passed_ss) {
        write_sensitivity(&out.join("ss_tab.tsv"), &result, sc, passed)?;
    }
    if let Some(zi) = &result.core.zero_ind {
        write_zero_ind(&out.join("zero_ind.tsv"), zi, &result.core)?;
    }
    write_metadata(
        &out.join("run_metadata.tsv"),
        &cfg,
        &result,
        elapsed,
        threads,
    )?;

    if !opts.quiet {
        let n_taxa = result.core.taxa.len();
        let n_sig = result.core.diff_abn.iter().filter(|b| **b).count();
        eprintln!(
            "ancombc2-rs: {n_taxa} taxa retained, {n_sig} of {} coefficient-wise calls, \
             {elapsed:.2}s on {threads} threads",
            result.core.diff_abn.len()
        );
        for w in &result.warnings {
            eprintln!("ancombc2-rs: warning: {w}");
        }
        if result.sensitivity.is_some() {
            let n_robust = result
                .diff_robust
                .as_ref()
                .map(|r| r.iter().filter(|b| **b).count())
                .unwrap_or(0);
            eprintln!("ancombc2-rs: {n_robust} robust calls after the sensitivity analysis");
        }
        eprintln!("ancombc2-rs: wrote {}", out.display());
    }
    Ok(())
}

fn write_sensitivity(
    path: &Path,
    result: &ancombc2_core::AncombcResult,
    sc: &ancombc2_core::SensitivityScores,
    passed: &[bool],
) -> std::io::Result<()> {
    use std::fmt::Write as _;
    let core = &result.core;
    let p = core.fix_eff.len();
    let n_taxa = core.taxa.len();
    let mut s = String::new();
    let _ = writeln!(s, "taxon");
    for name in &core.fix_eff {
        let _ = writeln!(s, "ss_prim_{name}");
    }
    for name in sc.colnames.iter().skip(p) {
        let _ = writeln!(s, "ss_pair_{name}");
    }
    for i in 0..n_taxa {
        let name = core
            .taxa
            .get(i)
            .and_then(|&t| core.taxon_names.get(t))
            .cloned()
            .unwrap_or_else(|| format!("taxon_{i}"));
        let _ = writeln!(s, "{name}");
        for c in 0..sc.colnames.len() {
            let v = sc
                .scores
                .get(i * sc.colnames.len() + c)
                .copied()
                .unwrap_or(f64::NAN);
            let _ = writeln!(
                s,
                "\t{}\t{}",
                ancombc2_io::table::format_double(v),
                passed[i * p + c.min(p - 1)]
            );
        }
    }
    std::fs::write(path, s)
}

fn write_zero_ind(
    path: &Path,
    zi: &ancombc2_core::ZeroIndication,
    core: &ancombc2_core::CoreOutput,
) -> std::io::Result<()> {
    use std::fmt::Write as _;
    let n_g = zi.groups.len();
    let mut s = String::new();
    let _ = writeln!(s, "taxon");
    for g in &zi.groups {
        let _ = writeln!(s, "zero_{g}");
    }
    for i in 0..zi.zero_ind.len().min(core.taxa_bias.len()) {
        let name = core
            .taxa_bias
            .get(i)
            .and_then(|&t| core.taxon_names.get(t))
            .cloned()
            .unwrap_or_else(|| format!("taxon_{i}"));
        let _ = writeln!(s, "{name}");
        for g in 0..n_g {
            let _ = write!(
                s,
                "\t{}",
                if zi.zero_ind[i * n_g + g] {
                    "TRUE"
                } else {
                    "FALSE"
                }
            );
        }
        let _ = writeln!(s);
    }
    std::fs::write(path, s)
}

fn write_metadata(
    path: &Path,
    cfg: &AncombcConfig,
    result: &ancombc2_core::AncombcResult,
    elapsed: f64,
    threads: usize,
) -> std::io::Result<()> {
    use std::fmt::Write as _;
    let core = &result.core;
    let mut s = String::new();
    let _ = writeln!(s, "key\tvalue");
    let _ = writeln!(s, "target\t{TARGET}");
    let _ = writeln!(s, "threads\t{threads}");
    let _ = writeln!(s, "wall_seconds\t{elapsed:.6}");
    let _ = writeln!(s, "peak_rss_kb\t{}", peak_rss_kb());
    let _ = writeln!(s, "bytes_allocated\t{}", BYTES.load(Ordering::Relaxed));
    let _ = writeln!(s, "allocations\t{}", ALLOCS.load(Ordering::Relaxed));
    let _ = writeln!(s, "peak_rss_source\t{}", rss_source());
    let _ = writeln!(s, "n_taxa_reported\t{}", core.taxa.len());
    let _ = writeln!(s, "n_taxa_bias\t{}", core.taxa_bias.len());
    let _ = writeln!(s, "n_samples\t{}", core.samples.len());
    let _ = writeln!(s, "n_fix_eff\t{}", core.fix_eff.len());
    let _ = writeln!(s, "fix_eff\t{}", core.fix_eff.join(";"));
    let _ = writeln!(s, "p_adj_method\t{:?}", cfg.p_adj_method);
    let _ = writeln!(s, "pseudo\t{}", cfg.pseudo);
    let _ = writeln!(s, "prv_cut\t{}", cfg.prv_cut);
    let _ = writeln!(s, "lib_cut\t{}", cfg.lib_cut);
    let _ = writeln!(s, "s0_perc\t{}", cfg.s0_perc);
    let _ = writeln!(s, "alpha\t{}", cfg.alpha);
    let _ = writeln!(s, "global\t{}", cfg.global);
    let _ = writeln!(s, "pairwise\t{}", cfg.pairwise);
    let _ = writeln!(s, "struc_zero\t{}", cfg.struc_zero);
    let _ = writeln!(s, "representation\t{:?}", cfg.representation);
    let _ = writeln!(s, "neg_lb\t{}", cfg.neg_lb);
    let _ = writeln!(s, "pseudo_sens\t{}", cfg.pseudo_sens);
    let _ = writeln!(s, "conservative\t{}", cfg.conservative);
    let _ = writeln!(s, "compat\t{:?}", cfg.compat);
    let _ = writeln!(s, "iter_tol\t{}", cfg.iter_control.tol);
    let _ = writeln!(s, "iter_max\t{}", cfg.iter_control.max_iter);
    let _ = writeln!(s, "em_tol\t{}", cfg.em_control.tol);
    let _ = writeln!(s, "em_max\t{}", cfg.em_control.max_iter);
    let _ = writeln!(s, "ml_iterations\t{}", core.ml_iterations);
    // The variance chain, one value per fixed effect: the E-M's bias and its
    // variance, the regulariser, and the post-regularisation variance. These are
    // length-`p`, not length-`n_taxa * p`, so they cost nothing to record -- and
    // they are the only way to attribute an `se` or `s0` disagreement to a link
    // rather than to "the numerics". `s0` in particular is a quantile of a column
    // built from these, so an error here reaches `s0` with no averaging-out over
    // taxa.
    let mut per_term = |name: &str, v: &[f64]| {
        let joined: Vec<String> = v.iter().map(|x| format!("{x:?}")).collect();
        let _ = writeln!(s, "{name}\t{}", joined.join(","));
    };
    per_term("delta_em", &core.delta_em);
    per_term("delta_wls", &core.delta_wls);
    per_term("var_delta", &core.var_delta);
    per_term("s02", &core.s02);
    per_term("convergence_trace", &core.ml_trace);
    for (name, secs) in core.timings.as_pairs() {
        if secs > 0.0 {
            let _ = writeln!(s, "stage_seconds.{name}\t{secs:.6}");
        }
    }
    let _ = writeln!(s, "stage_seconds_total\t{:.6}", core.timings.total());
    let _ = writeln!(
        s,
        "n_diff_abn\t{}",
        core.diff_abn.iter().filter(|b| **b).count()
    );
    if let Some(sc) = &result.sensitivity {
        let _ = writeln!(s, "sensitivity_runs\t{}", sc.pseudo.len());
        let _ = writeln!(
            s,
            "sensitivity_pseudo_min\t{}",
            sc.pseudo.first().copied().unwrap_or(f64::NAN)
        );
        let _ = writeln!(
            s,
            "sensitivity_pseudo_max\t{}",
            sc.pseudo.last().copied().unwrap_or(f64::NAN)
        );
    }
    let _ = writeln!(s, "n_warnings\t{}", result.warnings.len());
    for (i, w) in result.warnings.iter().enumerate() {
        let _ = writeln!(s, "warning_{i}\t{}", w.replace('\n', " "));
    }
    std::fs::write(path, s)
}

/// This process's peak resident set size, in kibibytes.
///
/// Read from `/proc/self/status` `VmHWM`, which is the kernel's own high-water
/// mark and therefore exact. The alternative -- a `getrusage(RUSAGE_CHILDREN)`
/// difference around the child -- does not work: `ru_maxrss` is a high-water mark
/// over *all* reaped children, so it never decreases and a second, smaller run
/// measures as zero.
fn peak_rss_kb() -> u64 {
    match std::fs::read_to_string("/proc/self/status") {
        Ok(text) => {
            for line in text.lines() {
                if let Some(rest) = line.strip_prefix("VmHWM:") {
                    let kb: String = rest
                        .trim()
                        .chars()
                        .take_while(|c| c.is_ascii_digit())
                        .collect();
                    return kb.parse().unwrap_or(0);
                }
            }
            0
        }
        Err(_) => 0,
    }
}

/// Where [`peak_rss_kb`] got its number, recorded so a reader knows whether it is
/// a measurement or an absence of one.
fn rss_source() -> &'static str {
    if Path::new("/proc/self/status").exists() {
        "/proc/self/status:VmHWM"
    } else {
        "unavailable"
    }
}

fn build_pool(threads: usize) {
    // Delegates to the core's single global pool, which is first-wins so a
    // library user who installed their own is not overridden. The analysis itself
    // runs inside it -- see the `with_pool` call below -- so this is the pool the
    // core's `par_iter`s actually draw from.
    //
    // This used to guard with `if rayon::current_num_threads() > 1 { return; }`,
    // which is backwards: before any pool is installed `current_num_threads`
    // reports the machine's logical CPU count, so on any multi-core host the guard
    // fired and **no pool was ever built**. `--threads` was silently ignored and
    // every benchmark arm measured Rayon's default global pool at full width --
    // which is why rust-1 through rust-16 all came out at the same time and P5
    // read as 1.00x scaling rather than as "the flag does nothing".
    ancombc2_core::parallel::install_pool(threads);
}

#[derive(Debug)]
struct Options {
    counts: PathBuf,
    meta: PathBuf,
    formula: String,
    group: Option<String>,
    p_adj_method: String,
    pseudo: f64,
    prv_cut: f64,
    lib_cut: f64,
    s0_perc: f64,
    alpha: f64,
    global: bool,
    pairwise: bool,
    mdfdr_fwer: String,
    struc_zero: bool,
    representation: ancombc2_core::preprocess::Representation,
    neg_lb: bool,
    iter_tol: f64,
    iter_max: usize,
    em_tol: f64,
    em_max: usize,
    pseudo_sens: bool,
    conservative: bool,
    compat: CompatMode,
    out: PathBuf,
    quiet: bool,
    threads: Option<usize>,
}

#[derive(Debug)]
enum Parsed {
    Help,
    Version,
    /// Boxed because `Options` is much larger than the other two variants and
    /// the enum is returned by value from a function that runs once.
    Run(Box<Options>),
}

impl Options {
    fn parse(args: &[String]) -> Result<Parsed, Box<dyn std::error::Error>> {
        let mut o = Options {
            counts: PathBuf::new(),
            meta: PathBuf::new(),
            formula: String::new(),
            group: None,
            p_adj_method: "holm".into(),
            pseudo: 0.0,
            prv_cut: 0.1,
            lib_cut: 0.0,
            s0_perc: 0.05,
            alpha: 0.05,
            global: false,
            pairwise: false,
            mdfdr_fwer: "holm".into(),
            struc_zero: false,
            representation: ancombc2_core::preprocess::Representation::Dense,
            neg_lb: false,
            iter_tol: 0.01,
            iter_max: 20,
            em_tol: 1e-5,
            em_max: 100,
            pseudo_sens: false,
            conservative: true,
            compat: CompatMode::Ancombc2_15,
            out: PathBuf::from("."),
            quiet: false,
            threads: None,
        };
        let mut have_formula = false;
        let mut i = 0;
        while i < args.len() {
            let a = args[i].as_str();
            let mut next = |what: &str| -> Result<String, Box<dyn std::error::Error>> {
                i += 1;
                args.get(i)
                    .cloned()
                    .ok_or_else(|| format!("{what} needs a value").into())
            };
            match a {
                "--help" | "-h" => return Ok(Parsed::Help),
                "--version" | "-V" => return Ok(Parsed::Version),
                "--counts" => o.counts = PathBuf::from(next("--counts")?),
                "--meta" => o.meta = PathBuf::from(next("--meta")?),
                "--formula" => {
                    o.formula = next("--formula")?;
                    have_formula = true;
                }
                "--group" => o.group = Some(next("--group")?),
                "--p-adj-method" => o.p_adj_method = next("--p-adj-method")?,
                "--pseudo" => o.pseudo = num(&next("--pseudo")?, "--pseudo")?,
                "--prv-cut" => o.prv_cut = num(&next("--prv-cut")?, "--prv-cut")?,
                "--lib-cut" => o.lib_cut = num(&next("--lib-cut")?, "--lib-cut")?,
                "--s0-perc" => o.s0_perc = num(&next("--s0-perc")?, "--s0-perc")?,
                "--alpha" => o.alpha = num(&next("--alpha")?, "--alpha")?,
                "--global" => o.global = true,
                "--pairwise" => o.pairwise = true,
                "--mdfdr-fwer" => o.mdfdr_fwer = next("--mdfdr-fwer")?,
                "--struc-zero" => o.struc_zero = true,
                "--sparse-taxa" => {
                    o.representation = ancombc2_core::preprocess::Representation::SparseTaxa
                }
                "--neg-lb" => o.neg_lb = true,
                "--iter-tol" => o.iter_tol = num(&next("--iter-tol")?, "--iter-tol")?,
                "--iter-max" => o.iter_max = next("--iter-max")?.parse()?,
                "--em-tol" => o.em_tol = num(&next("--em-tol")?, "--em-tol")?,
                "--em-max" => o.em_max = next("--em-max")?.parse()?,
                "--pseudo-sens" => o.pseudo_sens = true,
                "--conservative" => o.conservative = true,
                "--nonconservative" => o.conservative = false,
                "--compat" => o.compat = io::parse_compat(&next("--compat")?)?,
                "--out" => o.out = PathBuf::from(next("--out")?),
                "--quiet" => o.quiet = true,
                "--threads" => o.threads = Some(next("--threads")?.parse()?),
                other => return Err(format!("unknown option {other:?}; try --help").into()),
            }
            i += 1;
        }
        if o.counts.as_os_str().is_empty() {
            return Err("--counts is required".into());
        }
        if o.meta.as_os_str().is_empty() {
            return Err("--meta is required".into());
        }
        if !have_formula || o.formula.trim().is_empty() {
            return Err("--formula is required".into());
        }
        if (o.global || o.pairwise || o.struc_zero) && o.group.is_none() {
            return Err(
                "--group is required for --global, --pairwise and --struc-zero, because \
                 the reference requires it for those features"
                    .into(),
            );
        }
        if o.prv_cut < 0.0 || o.prv_cut > 1.0 {
            return Err("--prv-cut must be in [0, 1]".into());
        }
        if o.alpha <= 0.0 || o.alpha >= 1.0 {
            return Err("--alpha must be in (0, 1)".into());
        }
        if o.pseudo < 0.0 {
            return Err("--pseudo must be non-negative".into());
        }
        if !(0.0..=1.0).contains(&o.s0_perc) {
            return Err("--s0-perc must be in [0, 1]".into());
        }
        Ok(Parsed::Run(Box::new(o)))
    }
}

fn num(s: &str, what: &str) -> Result<f64, Box<dyn std::error::Error>> {
    s.parse::<f64>()
        .map_err(|_| format!("{what}: {s:?} is not a number").into())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn args(v: &[&str]) -> Vec<String> {
        v.iter().map(|s| s.to_string()).collect()
    }

    #[test]
    fn help_and_version_short_circuit() {
        assert!(matches!(
            Options::parse(&args(&["--help"])).unwrap(),
            Parsed::Help
        ));
        assert!(matches!(
            Options::parse(&args(&["--version"])).unwrap(),
            Parsed::Version
        ));
    }

    #[test]
    fn requires_the_three_inputs() {
        let e = Options::parse(&args(&["--meta", "m.tsv", "--formula", "~ a"])).unwrap_err();
        assert!(e.to_string().contains("--counts"));
        let e = Options::parse(&args(&["--counts", "c.tsv", "--formula", "~ a"])).unwrap_err();
        assert!(e.to_string().contains("--meta"));
        let e = Options::parse(&args(&["--counts", "c.tsv", "--meta", "m.tsv"])).unwrap_err();
        assert!(e.to_string().contains("--formula"));
    }

    #[test]
    fn defaults_match_the_reference() {
        let o = Options::parse(&args(&[
            "--counts",
            "c.tsv",
            "--meta",
            "m.tsv",
            "--formula",
            "~ a",
        ]))
        .map(|p| match p {
            Parsed::Run(o) => *o,
            _ => panic!("expected a run"),
        })
        .unwrap();
        assert_eq!(o.p_adj_method, "holm");
        assert_eq!(o.pseudo, 0.0);
        assert_eq!(o.prv_cut, 0.1);
        assert_eq!(o.s0_perc, 0.05);
        assert_eq!(o.alpha, 0.05);
        assert_eq!(o.iter_max, 20);
        assert_eq!(o.em_max, 100);
        assert!(o.conservative, "the reference's default is conservative");
        assert!(!o.pseudo_sens);
        assert!(!o.global && !o.pairwise && !o.struc_zero);
    }

    #[test]
    fn multi_group_features_require_a_group() {
        for flag in ["--global", "--pairwise", "--struc-zero"] {
            let e = Options::parse(&args(&[
                "--counts",
                "c.tsv",
                "--meta",
                "m.tsv",
                "--formula",
                "~ a",
                flag,
            ]))
            .unwrap_err();
            assert!(e.to_string().contains("--group is required"), "{flag}: {e}");
        }
    }

    #[test]
    fn rejects_out_of_range_numbers() {
        let base = ["--counts", "c.tsv", "--meta", "m.tsv", "--formula", "~ a"];
        for (flag, v) in [
            ("--prv-cut", "2"),
            ("--alpha", "0"),
            ("--s0-perc", "1.5"),
            ("--pseudo", "-1"),
        ] {
            let mut a = base.to_vec();
            a.push(flag);
            a.push(v);
            assert!(
                Options::parse(&args(&a)).is_err(),
                "{flag} {v} should be rejected"
            );
        }
    }

    #[test]
    fn rejects_unknown_options() {
        let e = Options::parse(&args(&["--nope"])).unwrap_err();
        assert!(e.to_string().contains("unknown option"));
    }
}
