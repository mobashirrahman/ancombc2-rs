//! `ancombc2-sim`: Layer 3 statistical validation.
//!
//! Numerical parity with the oracle proves the port matches R. It does not
//! prove either is *right*. This binary rebuilds the simulation design from
//! `PLAN.md` §5.3, runs the Rust arm in-process, and summarises the cells
//! against the R arm that `scripts/sim_r.R` produces over the same generated
//! tables.
//!
//! ```text
//! ancombc2-sim generate --grid validation/simulation/full --out data
//! ancombc2-sim rust     --grid validation/simulation/full --out rust.jsonl
//! ancombc2-sim r        --grid validation/simulation/full --out r.jsonl
//! ancombc2-sim summarise --grid validation/simulation/full \
//!     --rust rust.jsonl --r r.jsonl --out summary.json
//! ```
//!
//! Every subcommand prints what it did, and the summary states plainly whether
//! the run had an R arm to compare against, so a Rust-only run can never be
//! read as a passing comparison.

// Several loops here walk an index and two or three parallel buffers by that
// index (`is_da[i]`, `group[j]`, `data[i * n_samp + j]`). Zipping them would
// hide which buffer is which, and the generator is clearer when the row and
// column indices are explicit.
#![allow(clippy::needless_range_loop)]

pub mod grid;
pub mod metrics;
pub mod rng;
pub mod run;
pub mod simulate;
pub mod summarise;

use std::collections::HashMap;
use std::io::Write;
use std::path::{Path, PathBuf};

use clap::{Parser, Subcommand};
use grid::Grid;
use metrics::RepMetrics;
use summarise::GridSummary;

#[derive(Parser, Debug)]
#[command(
    name = "ancombc2-sim",
    about = "Layer 3 statistical validation for ancombc2-rs",
    long_about = "Simulation grid, truth-aware metrics, and the Rust-vs-R comparison. \
                  Parity with the oracle proves the port matches R; this proves whether \
                  the method, as ported, recovers what was planted."
)]
struct Cli {
    #[command(subcommand)]
    cmd: Cmd,
}

#[derive(Subcommand, Debug)]
enum Cmd {
    /// Expand the grid and write every replicate's counts, metadata, and truth.
    Generate {
        #[arg(long)]
        grid: std::path::PathBuf,
        /// Root for the per-replicate directories, one per cell and rep.
        #[arg(long)]
        out: std::path::PathBuf,
        /// Write only this cell, for running a single cell through R.
        #[arg(long)]
        cell: Option<usize>,
        /// Write only this rep, for debugging a single replicate.
        #[arg(long)]
        rep: Option<usize>,
    },
    /// Run the Rust arm over the grid and write one JSONL row per replicate.
    Rust {
        #[arg(long)]
        grid: std::path::PathBuf,
        #[arg(long)]
        out: std::path::PathBuf,
        /// Restrict to these cell indices.
        #[arg(long, value_delimiter = ',')]
        cells: Option<Vec<usize>>,
        /// Restrict to this many reps per cell, overriding the grid.
        #[arg(long)]
        reps: Option<usize>,
        /// Report per-cell progress on stderr.
        #[arg(long)]
        progress: bool,
    },
    /// Run the R arm over data already written by `generate`.
    ///
    /// The R script is `scripts/sim_r.R`; this subcommand is a thin wrapper that
    /// locates Rscript and passes the arguments through, so a CI job and a local
    /// run take the same path.
    R {
        #[arg(long)]
        data: std::path::PathBuf,
        #[arg(long)]
        out: std::path::PathBuf,
        /// Override the R script path.
        #[arg(long, default_value = "scripts/sim_r.R")]
        script: std::path::PathBuf,
        /// An explicit Rscript, for an R outside `PATH`.
        #[arg(long, default_value = "Rscript")]
        rscript: String,
        /// Library search paths, `R_LIBS` style, for an out-of-tree install.
        #[arg(long, value_delimiter = ':')]
        r_libs: Option<Vec<String>>,
    },
    /// Pool the arms into per-cell summaries and apply the acceptance rule.
    Summarise {
        #[arg(long)]
        grid: std::path::PathBuf,
        #[arg(long)]
        rust: std::path::PathBuf,
        /// The R arm's rows. Without it the summary is Rust-only and says so.
        #[arg(long)]
        r: Option<std::path::PathBuf>,
        #[arg(long)]
        out: Option<std::path::PathBuf>,
        /// Exit non-zero when any cell is divergent. Off by default, because
        /// `PLAN.md` treats a divergence as a finding rather than a failure.
        #[arg(long)]
        fail_on_divergence: bool,
    },
    /// Print a grid's cells without running anything, to check a spec.
    Cells {
        #[arg(long)]
        grid: std::path::PathBuf,
    },
    /// Layer 4: run the Rust arm over the committed real datasets.
    ///
    /// Reads each `validation/realdata/<name>/` directory -- the same
    /// `counts.tsv`, `meta.tsv` and `analysis.json` the R arm reads -- so the
    /// two arms see identical input and the comparison is between
    /// implementations rather than between two reconstructions of the dataset.
    Realdata {
        /// The dataset root, holding one directory per dataset.
        #[arg(long, default_value = "validation/realdata")]
        data: std::path::PathBuf,
        #[arg(long)]
        out: std::path::PathBuf,
        /// Only these datasets; all of them by default.
        #[arg(long, value_delimiter = ',')]
        datasets: Option<Vec<String>>,
    },
}

fn main() -> std::process::ExitCode {
    let cli = Cli::parse();
    let r = match cli.cmd {
        Cmd::Generate {
            grid,
            out,
            cell,
            rep,
        } => cmd_generate(&grid, &out, cell, rep),
        Cmd::Rust {
            grid,
            out,
            cells,
            reps,
            progress,
        } => cmd_rust(&grid, &out, cells, reps, progress),
        Cmd::R {
            data,
            out,
            script,
            rscript,
            r_libs,
        } => cmd_r(&data, &out, &script, &rscript, r_libs),
        Cmd::Summarise {
            grid,
            rust,
            r,
            out,
            fail_on_divergence,
        } => cmd_summarise(
            &grid,
            &rust,
            r.as_deref(),
            out.as_deref(),
            fail_on_divergence,
        ),
        Cmd::Cells { grid } => cmd_cells(&grid),
        Cmd::Realdata {
            data,
            out,
            datasets,
        } => cmd_realdata(&data, &out, datasets),
    };
    match r {
        Ok(()) => std::process::ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("ancombc2-sim: {e}");
            std::process::ExitCode::FAILURE
        }
    }
}

fn load(path: &Path) -> Result<Grid, String> {
    Grid::load(path)
}

fn cmd_generate(
    grid_path: &Path,
    out: &Path,
    cell: Option<usize>,
    rep: Option<usize>,
) -> Result<(), String> {
    let g = load(grid_path)?;
    let cells = g.cells();
    let want: Vec<grid::Cell> = match cell {
        Some(i) => cells.into_iter().filter(|c| c.index == i).collect(),
        None => cells,
    };
    if want.is_empty() {
        return Err(format!("no cell with index {cell:?} in this grid"));
    }
    let reps: Vec<usize> = match rep {
        Some(r) => vec![r],
        None => (0..g.reps).collect(),
    };
    let ctl = run::sim_control(&g);
    let mut n = 0usize;
    let mut bytes = 0u64;
    for c in &want {
        for r_idx in &reps {
            let rep_data = simulate::simulate(c, g.seed, *r_idx, &ctl);
            let dir = out
                .join(format!("cell{:04}", c.index))
                .join(format!("rep{r_idx:04}"));
            run::write_replicate(&dir, &rep_data, c, *r_idx)
                .map_err(|e| format!("writing {}: {e}", dir.display()))?;
            if let Ok(md) = std::fs::metadata(dir.join("counts.tsv")) {
                bytes += md.len();
            }
            n += 1;
        }
    }
    println!(
        "wrote {n} replicate(s) for {} cell(s) to {} ({:.1} MiB of counts)",
        want.len(),
        out.display(),
        bytes as f64 / (1024.0 * 1024.0)
    );
    Ok(())
}

fn write_jsonl(path: &Path, rows: &[RepMetrics]) -> Result<(), String> {
    if let Some(p) = path.parent() {
        if !p.as_os_str().is_empty() {
            std::fs::create_dir_all(p).map_err(|e| format!("creating {}: {e}", p.display()))?;
        }
    }
    let f = std::fs::File::create(path)
        .map_err(|e| format!("cannot create {}: {e}", path.display()))?;
    let mut w = std::io::BufWriter::new(f);
    for r in rows {
        let line = serde_json::to_string(r).map_err(|e| format!("serialising a replicate: {e}"))?;
        writeln!(w, "{line}").map_err(|e| format!("writing {}: {e}", path.display()))?;
    }
    w.flush()
        .map_err(|e| format!("flushing {}: {e}", path.display()))?;
    Ok(())
}

/// Read a JSONL file, transparently decompressing a `.gz`.
///
/// The committed per-replicate rows are gzipped — 960 replicates of 100 taxa is
/// 23 MB raw — and a reviewer re-running `summarise` against them is exactly the
/// person who should not have to gunzip by hand first. The decompression shells
/// out to `gzip` rather than pulling in a compression dependency, which is the
/// right trade for a test harness: one external binary, already a prerequisite of
/// the benchmark scripts, versus a crate in the dependency graph of the
/// validator.
fn read_jsonl(path: &Path) -> Result<Vec<RepMetrics>, String> {
    let rows: Vec<RepMetrics> = read_jsonl_streaming(path)?.collect::<Result<Vec<_>, _>>()?;
    Ok(rows)
}

/// Stream a JSONL file a row at a time.
///
/// This was `fs::read_to_string` followed by a split into lines, which holds the
/// whole file as one `String` *and* then all the parsed rows at once. The `full`
/// grid's Rust arm is 249,600 rows and **17 GB**, so that needed more memory than
/// the machine has: `summarise` was killed by the OOM reaper rather than
/// reporting anything.
///
/// A line-oriented reader costs nothing extra and makes the summariser's memory
/// a function of the largest row rather than of the file, which is the only
/// property that lets a 17 GB artefact be summarised on a 31 GB host.
///
/// A `.gz` path still shells out to `gzip -dc` (see below for why), but pipes it
/// rather than buffering the whole decompression.
fn read_jsonl_streaming(
    path: &Path,
) -> Result<Box<dyn Iterator<Item = Result<RepMetrics, String>> + Send>, String> {
    let gz = path
        .extension()
        .map(|e| e.eq_ignore_ascii_case("gz"))
        .unwrap_or(false);
    use std::io::BufRead as _;
    if gz {
        let child = std::process::Command::new("gzip")
            .arg("-dc")
            .arg(path)
            .stdout(std::process::Stdio::piped())
            .spawn()
            .map_err(|e| {
                format!(
                    "cannot decompress {}: {e}; install gzip, or pass the \
                     uncompressed file",
                    path.display()
                )
            })?;
        let reader = std::io::BufReader::with_capacity(
            1 << 20,
            child
                .stdout
                .ok_or_else(|| "gzip produced no stdout".to_string())?,
        );
        let name = path.display().to_string();
        let n = name.clone();
        return Ok(Box::new(
            reader
                .lines()
                .enumerate()
                .filter(|(_, l)| !matches!(l, Ok(l) if l.trim().is_empty()))
                .map(move |(i, line)| {
                    line.map_err(|e| format!("{n}:{}: {e}", i + 1))
                        .and_then(|l| {
                            serde_json::from_str(&l).map_err(|e| format!("{n}:{}: {e}", i + 1))
                        })
                }),
        ));
    }
    let file =
        std::fs::File::open(path).map_err(|e| format!("cannot read {}: {e}", path.display()))?;
    let reader = std::io::BufReader::with_capacity(1 << 20, file);
    let name = path.display().to_string();
    Ok(Box::new(
        reader
            .lines()
            .enumerate()
            .filter(|(_, l)| !matches!(l, Ok(l) if l.trim().is_empty()))
            .map(move |(i, line)| {
                line.map_err(|e| format!("{name}:{}: {e}", i + 1))
                    .and_then(|l| {
                        serde_json::from_str(&l).map_err(|e| format!("{name}:{}: {e}", i + 1))
                    })
            }),
    ))
}

fn cmd_rust(
    grid_path: &Path,
    out: &Path,
    cells: Option<Vec<usize>>,
    reps: Option<usize>,
    progress: bool,
) -> Result<(), String> {
    let mut g = load(grid_path)?;
    if let Some(r) = reps {
        // A rep override is a *ceiling*, not a replacement: a block that asked
        // for fewer reps keeps its own count, so lowering the override cannot
        // silently lengthen the sensitivity blocks.
        g.reps = r;
        for b in &mut g.blocks {
            if let Some(br) = b.reps {
                b.reps = Some(br.min(r));
            }
        }
    }
    let all = g.cells();
    let selected: Vec<grid::Cell> = match cells {
        Some(want) => {
            for i in &want {
                if !all.iter().any(|c| c.index == *i) {
                    return Err(format!("no cell with index {i} in this grid"));
                }
            }
            all.into_iter()
                .filter(|c| want.contains(&c.index))
                .collect()
        }
        None => all,
    };
    let total: usize = selected.iter().map(|c| g.reps_for(c)).sum();
    println!(
        "rust arm: {} cell(s), {total} replicate(s), seed {}",
        selected.len(),
        g.seed
    );
    let t0 = std::time::Instant::now();
    let mut rows = Vec::with_capacity(total);
    for c in &selected {
        let got: Vec<RepMetrics> = run::run_cell(&g, c, "rust");
        let failed = got.iter().filter(|r| r.error.is_some()).count();
        if progress {
            eprintln!(
                "cell {:>4} n_taxa={} n_samp={} da={} zinf={} conf={} reps={} failed={}",
                c.index,
                c.n_taxa,
                c.n_samp,
                c.da_proportion,
                c.zero_inflation,
                c.confound,
                got.len(),
                failed
            );
        }
        rows.extend(got);
    }
    write_jsonl(out, &rows)?;
    let (ok, failed) = summarise::count_rows(&rows, "rust");
    println!(
        "wrote {} row(s) to {}: {ok} analysed, {failed} failed, {:.1}s",
        rows.len(),
        out.display(),
        t0.elapsed().as_secs_f64()
    );
    Ok(())
}

fn cmd_r(
    data: &Path,
    out: &Path,
    script: &Path,
    rscript: &str,
    r_libs: Option<Vec<String>>,
) -> Result<(), String> {
    if !script.exists() {
        return Err(format!("R script {} not found", script.display()));
    }
    let mut cmd = std::process::Command::new(rscript);
    cmd.arg("--vanilla")
        .arg(script)
        .arg("--data")
        .arg(data)
        .arg("--out")
        .arg(out);
    if let Some(libs) = &r_libs {
        cmd.env("R_LIBS", libs.join(":"));
    }
    let status = cmd
        .status()
        .map_err(|e| format!("cannot run {rscript}: {e}; is R installed and on PATH?"))?;
    if !status.success() {
        return Err(format!("{} exited with {status}", script.display()));
    }
    println!("r arm: wrote {}", out.display());
    Ok(())
}

fn cmd_summarise(
    grid_path: &Path,
    rust_path: &Path,
    r_path: Option<&Path>,
    out: Option<&Path>,
    fail_on_divergence: bool,
) -> Result<(), String> {
    let g = load(grid_path)?;
    let mut rows = read_jsonl(rust_path)?;
    let have_r = r_path.is_some();
    if let Some(p) = r_path {
        rows.extend(read_jsonl(p)?);
    }
    let summary: GridSummary = summarise::summarise(&g, &rows);
    print!("{}", summarise::table(&summary));
    if let Some(p) = out {
        if let Some(dir) = p.parent() {
            if !dir.as_os_str().is_empty() {
                std::fs::create_dir_all(dir)
                    .map_err(|e| format!("creating {}: {e}", dir.display()))?;
            }
        }
        let text = serde_json::to_string_pretty(&summary)
            .map_err(|e| format!("serialising the summary: {e}"))?;
        std::fs::write(p, format!("{text}\n"))
            .map_err(|e| format!("writing {}: {e}", p.display()))?;
        println!("wrote summary to {}", p.display());
    }
    if !have_r {
        println!(
            "\nno R arm was supplied, so no cell was compared. This is a Rust-only run:\n\
             it reports the Rust arm's calibration and power, and says nothing about whether\n\
             the two implementations agree."
        );
    }
    if fail_on_divergence && summary.n_divergent > 0 {
        return Err(format!(
            "{} divergent cell(s); see findings[] in the summary",
            summary.n_divergent
        ));
    }
    Ok(())
}

/// One prepared real dataset's analysis configuration.
///
/// Read from a two-column key/value table rather than JSON, because the file is
/// written by `scripts/prepare_realdata.R` without a JSON dependency on either
/// side and a formula containing `~` or `+` is then just a line. An unknown key
/// is an error rather than a silent default: a typo in a dataset's
/// configuration that quietly reverts to the default would produce a comparison
/// of something other than what the dataset says.
fn read_analysis(path: &Path) -> Result<std::collections::HashMap<String, String>, String> {
    let text = std::fs::read_to_string(path)
        .map_err(|e| format!("cannot read {}: {e}", path.display()))?;
    let mut map = std::collections::HashMap::new();
    for (i, line) in text.lines().enumerate() {
        if line.trim().is_empty() || line.starts_with('#') {
            continue;
        }
        if i == 0 && line.starts_with("key\t") {
            continue;
        }
        let mut f = line.splitn(2, '\t');
        let k = f.next().unwrap_or("").trim().to_string();
        let v = f.next().unwrap_or("").trim().to_string();
        if k.is_empty() {
            continue;
        }
        if map.insert(k.clone(), v).is_some() {
            return Err(format!("{}: duplicate key `{k}`", path.display()));
        }
    }
    Ok(map)
}

fn as_bool(m: &std::collections::HashMap<String, String>, k: &str) -> Option<bool> {
    m.get(k).map(|v| v == "true")
}

fn as_f64(
    m: &std::collections::HashMap<String, String>,
    k: &str,
    path: &Path,
) -> Result<Option<f64>, String> {
    match m.get(k) {
        None => Ok(None),
        Some(v) => v
            .parse::<f64>()
            .map(Some)
            .map_err(|e| format!("{}: `{k}` is {v:?}, not a number: {e}", path.display())),
    }
}

/// Minimal reader for a real dataset's `counts.tsv` and `meta.tsv`.
///
/// Deliberately not `ancombc2_io`: the IO crate's reader owns the filtering and
/// formula handling that a real run wants, and this needs only the raw table
/// plus a design matrix built from the same formula string. Duplicating it would
/// mean two readers whose disagreement shows up as a data difference.
/// One prepared dataset: its counts, the metadata as name-keyed rows, and the
/// metadata's variable names.
struct Dataset {
    counts: ancombc2_core::preprocess::CountMatrix,
    taxa: Vec<String>,
    meta: Vec<HashMap<String, String>>,
    samples: Vec<String>,
    vars: Vec<String>,
}

fn read_dataset(dir: &Path) -> Result<Dataset, String> {
    let text = std::fs::read_to_string(dir.join("counts.tsv"))
        .map_err(|e| format!("cannot read {}/counts.tsv: {e}", dir.display()))?;
    let mut lines = text.lines().filter(|l| !l.trim().is_empty());
    let header: Vec<&str> = lines
        .next()
        .ok_or_else(|| format!("{}/counts.tsv is empty", dir.display()))?
        .split('\t')
        .collect();
    let samples: Vec<String> = header.iter().map(|s| s.to_string()).collect();
    let mut taxa = Vec::new();
    let mut data = Vec::new();
    for line in lines {
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() != samples.len() + 1 {
            return Err(format!(
                "{}/counts.tsv: row has {} fields, expected {}",
                dir.display(),
                f.len(),
                samples.len() + 1
            ));
        }
        taxa.push(f[0].to_string());
        for v in &f[1..] {
            data.push(v.parse::<f64>().map_err(|e| {
                format!("{}/counts.tsv: {:?} is not a count: {e}", dir.display(), v)
            })?);
        }
    }
    let n_taxa = taxa.len();

    let mtext = std::fs::read_to_string(dir.join("meta.tsv"))
        .map_err(|e| format!("cannot read {}/meta.tsv: {e}", dir.display()))?;
    let mut mlines = mtext.lines().filter(|l| !l.trim().is_empty());
    let mheader: Vec<&str> = mlines
        .next()
        .ok_or_else(|| format!("{}/meta.tsv is empty", dir.display()))?
        .split('\t')
        .collect();
    let mut meta_names = Vec::new();
    let mut meta_rows: Vec<HashMap<String, String>> = Vec::new();
    for line in mlines {
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() != mheader.len() + 1 {
            return Err(format!(
                "{}/meta.tsv: row has {} fields, expected {}",
                dir.display(),
                f.len(),
                mheader.len() + 1
            ));
        }
        meta_names.push(f[0].to_string());
        let mut row = HashMap::new();
        for (j, k) in mheader.iter().enumerate() {
            row.insert((*k).to_string(), f[j + 1].to_string());
        }
        meta_rows.push(row);
    }
    if meta_names != samples {
        return Err(format!(
            "{}/meta.tsv samples do not match counts.tsv: {} vs {}",
            dir.display(),
            meta_names.len(),
            samples.len()
        ));
    }

    let counts = ancombc2_core::preprocess::CountMatrix::new(n_taxa, samples.len(), data)
        .map_err(|e| format!("{}: {e}", dir.display()))?;
    // The design is built by the caller from the formula.
    Ok(Dataset {
        counts,
        taxa,
        meta: meta_rows,
        samples,
        vars: mheader.iter().map(|s| s.to_string()).collect(),
    })
}

fn cmd_realdata(root: &Path, out: &Path, only: Option<Vec<String>>) -> Result<(), String> {
    use ancombc2_core::config::{AdjustMethod, AncombcConfig};
    use ancombc2_core::matrix::Matrix;
    use ancombc2_core::pipeline::ancombc2_run_named;

    let mut dirs: Vec<PathBuf> = std::fs::read_dir(root)
        .map_err(|e| format!("cannot read {}: {e}", root.display()))?
        .filter_map(|e| e.ok())
        .map(|e| e.path())
        .filter(|p| p.is_dir() && p.join("counts.tsv").exists())
        .collect();
    dirs.sort();
    if let Some(want) = &only {
        dirs.retain(|d| {
            want.iter().any(|w| {
                d.file_name()
                    .map(|n| n.to_string_lossy() == w.as_str())
                    .unwrap_or(false)
            })
        });
        if dirs.is_empty() {
            return Err(format!(
                "no dataset matched {want:?} under {}",
                root.display()
            ));
        }
    }
    if dirs.is_empty() {
        return Err(format!("no dataset directories under {}", root.display()));
    }

    if let Some(p) = out.parent() {
        if !p.as_os_str().is_empty() {
            std::fs::create_dir_all(p).map_err(|e| format!("creating {}: {e}", p.display()))?;
        }
    }
    let mut w = std::io::BufWriter::new(
        std::fs::File::create(out).map_err(|e| format!("cannot create {}: {e}", out.display()))?,
    );
    let mut wrote = 0usize;

    for dir in &dirs {
        let name = dir
            .file_name()
            .map(|n| n.to_string_lossy().into_owned())
            .unwrap_or_else(|| dir.display().to_string());
        let t0 = std::time::Instant::now();
        let row = (|| -> Result<String, String> {
            let cfg_path = dir.join("analysis.tsv");
            let aj = read_analysis(&cfg_path)?;
            let aj_formula = aj
                .get("formula")
                .cloned()
                .ok_or_else(|| format!("{}: no `formula` key", cfg_path.display()))?;
            let aj_group = aj.get("group").cloned();
            // The factor's level order, reference first. This is stated in the
            // dataset rather than inferred from the order the labels appear in,
            // because R's `factor()` sorts alphabetically and would otherwise
            // build the *opposite* contrast: the two arms would report
            // coefficients with the same meaning under different names, and
            // every sign would be inverted.
            let aj_group_levels: Option<Vec<String>> = aj.get("group_levels").map(|v| {
                v.split(',')
                    .map(|t| t.trim().to_string())
                    .filter(|t| !t.is_empty())
                    .collect()
            });
            let ds = read_dataset(dir)?;
            let Dataset {
                counts,
                taxa,
                meta,
                samples,
                vars,
            } = ds;

            // The design: an intercept, one column per formula term. A term that
            // names the group variable becomes a factor contrast against its
            // first level, exactly as R's `model.matrix` would, and any other
            // term is read as a numeric column.
            let terms: Vec<String> = aj_formula
                .split('+')
                .map(|t| t.trim().trim_start_matches('~').trim().to_string())
                .filter(|t| !t.is_empty())
                .collect();
            if terms.is_empty() {
                return Err(format!(
                    "{}/analysis.json has an empty formula",
                    dir.display()
                ));
            }
            let mut colnames = vec!["(Intercept)".to_string()];
            let mut columns: Vec<Vec<f64>> = vec![vec![1.0; samples.len()]];
            for t in &terms {
                if !vars.contains(t) {
                    return Err(format!(
                        "{}: the formula names `{t}`, which is not a \
                         metadata column (have: {})",
                        dir.display(),
                        vars.join(", ")
                    ));
                }
                let raw: Vec<String> = meta
                    .iter()
                    .map(|r| r.get(t).cloned().unwrap_or_default())
                    .collect();
                if aj_group.as_deref() == Some(t.as_str()) {
                    // A factor: one column per level past the first. The level
                    // order is the dataset's `group_levels` when it gives one,
                    // and otherwise the order the labels first appear in --
                    // which is checked against the data rather than trusted.
                    let mut levels: Vec<String> = Vec::new();
                    for v in &raw {
                        if v.is_empty() || v == "NA" {
                            continue;
                        }
                        if !levels.contains(v) {
                            levels.push(v.clone());
                        }
                    }
                    // A sample whose group label is missing gets an *incomplete*
                    // design row, not a fabricated reference-level value.
                    //
                    // R's `model.matrix` gives NA for a factor with no level and
                    // `lm.fit` then drops that row, so those samples are excluded
                    // from the fit. Mapping a missing label onto the reference
                    // level instead adds the sample to the reference group as if
                    // it had been observed there -- which biases every
                    // coefficient toward the reference. On `atlas1006`, 37 of
                    // 1151 samples have `sex == "NA"`, and doing it the other way
                    // moved `beta` by up to 0.105 and flipped two `diff_abn`
                    // calls. The core already skips incomplete rows in the
                    // sandwich; it needs the NaN to know.
                    let missing: Vec<bool> =
                        raw.iter().map(|v| v.is_empty() || v == "NA").collect();
                    if let Some(declared) = &aj_group_levels {
                        // `NA` and the empty string are missing values here too,
                        // for the same reason as above: they are not levels.
                        let missing: Vec<&String> =
                            levels.iter().filter(|l| !declared.contains(l)).collect();
                        let extra: Vec<&String> =
                            declared.iter().filter(|l| !levels.contains(l)).collect();
                        if !missing.is_empty() || !extra.is_empty() {
                            return Err(format!(
                                "{}: `group_levels` is {declared:?} but the data has \
                                 {levels:?} (unlisted: {missing:?}, absent: {extra:?})",
                                cfg_path.display()
                            ));
                        }
                        levels = declared.clone();
                    }
                    if levels.len() < 2 {
                        return Err(format!(
                            "{}/{t} has {} level(s); a group needs at least 2",
                            dir.display(),
                            levels.len()
                        ));
                    }
                    for lv in levels.iter().skip(1) {
                        columns.push(
                            raw.iter()
                                .zip(&missing)
                                .map(|(v, m)| {
                                    if *m {
                                        f64::NAN
                                    } else if v == lv {
                                        1.0
                                    } else {
                                        0.0
                                    }
                                })
                                .collect(),
                        );
                        colnames.push(format!("{t}{lv}"));
                    }
                } else {
                    // Numeric if every present value parses, a factor otherwise.
                    // This is what R's `model.matrix` does with a data frame
                    // column: a character column becomes a factor, and the
                    // contrast columns are then named after its levels. Reading
                    // every non-group term as numeric would fail on any
                    // categorical covariate, and *silently* treating a
                    // categorical term as numeric is not an option at all.
                    //
                    // A bare `NA` is a missing value, because that is
                    // `read.delim`'s default `na.strings`. Treating it as a
                    // level would add a spurious factor level -- which is a real
                    // hazard: `atlas1006` has `sex == "NA"` on some samples.
                    let present: Vec<&String> =
                        raw.iter().filter(|v| !v.is_empty() && *v != "NA").collect();
                    let all_numeric =
                        !present.is_empty() && present.iter().all(|v| v.parse::<f64>().is_ok());
                    if all_numeric {
                        // A missing value becomes NaN, and the core drops
                        // incomplete design rows the way `lm.fit` does -- which
                        // is what R's `model.matrix` produces for a factor level
                        // that does not exist. Fabricating a value here instead
                        // (mapping a missing group onto the reference level)
                        // biases every coefficient towards the reference.
                        let col: Vec<f64> = raw
                            .iter()
                            .map(|v| v.parse::<f64>().unwrap_or(f64::NAN))
                            .collect();
                        columns.push(col);
                        colnames.push(t.clone());
                    } else {
                        let mut levels: Vec<String> = Vec::new();
                        for v in &present {
                            if !levels.contains(v) {
                                levels.push((*v).clone());
                            }
                        }
                        if levels.len() < 2 {
                            return Err(format!(
                                "{}: column `{t}` is neither numeric nor a factor with \
                                 at least two levels (present values: {:?})",
                                dir.display(),
                                levels
                            ));
                        }
                        for lv in levels.iter().skip(1) {
                            columns.push(
                                raw.iter()
                                    .map(|v| if v == lv { 1.0 } else { 0.0 })
                                    .collect(),
                            );
                            colnames.push(format!("{t}{lv}"));
                        }
                    }
                }
            }
            let n_samp = samples.len();
            let mut x = Matrix::zeros(n_samp, columns.len());
            for (j, col) in columns.iter().enumerate() {
                for i in 0..n_samp {
                    x.set(i, j, col[i]);
                }
            }

            // The group index, for the test machinery and the structural zeros.
            let group_index: Option<Vec<usize>> = aj_group.as_ref().map(|g| {
                // The same level order the design used, so the test machinery's
                // group numbering matches the coefficient columns.
                let mut levels: Vec<String> = aj_group_levels.clone().unwrap_or_default();
                if levels.is_empty() {
                    for v in meta.iter().map(|r| r.get(g).cloned().unwrap_or_default()) {
                        if v.is_empty() || v == "NA" {
                            continue;
                        }
                        if !levels.contains(&v) {
                            levels.push(v);
                        }
                    }
                }
                // A missing group label has no level. It becomes `usize::MAX`,
                // which the core treats as "not in any group": the sample is
                // still a row of the design, and the per-taxon fits drop it the
                // way `lm.fit` drops an incomplete row. Panicking here instead
                // would make one missing label abort the whole dataset, which is
                // not what R does.
                meta.iter()
                    .map(|r| {
                        let v = r.get(g).cloned().unwrap_or_default();
                        levels.iter().position(|l| l == &v).unwrap_or(usize::MAX)
                    })
                    .collect()
            });

            let unknown: Vec<&String> = aj
                .keys()
                .filter(|k| {
                    !matches!(
                        k.as_str(),
                        "formula"
                            | "group"
                            | "group_levels"
                            | "p_adj_method"
                            | "pseudo"
                            | "pseudo_sens"
                            | "conservative"
                            | "prv_cut"
                            | "lib_cut"
                            | "s0_perc"
                            | "struc_zero"
                            | "neg_lb"
                            | "alpha"
                            | "global"
                            | "pairwise"
                            | "coefficient"
                    )
                })
                .collect();
            if !unknown.is_empty() {
                return Err(format!(
                    "{}: unknown key(s) {:?}; a typo would silently fall back to a \
                     default and compare something other than what the dataset says",
                    cfg_path.display(),
                    unknown
                ));
            }
            let p_adj = aj
                .get("p_adj_method")
                .cloned()
                .unwrap_or_else(|| "BH".to_string());
            let cfg = AncombcConfig {
                fix_eff: colnames.clone(),
                p_adj_method: AdjustMethod::parse(&p_adj)
                    .map_err(|e| format!("{}: {e}", cfg_path.display()))?,
                pseudo: as_f64(&aj, "pseudo", &cfg_path)?.unwrap_or(0.5),
                pseudo_sens: as_bool(&aj, "pseudo_sens").unwrap_or(false),
                conservative: as_bool(&aj, "conservative").unwrap_or(true),
                prv_cut: as_f64(&aj, "prv_cut", &cfg_path)?.unwrap_or(0.0),
                lib_cut: as_f64(&aj, "lib_cut", &cfg_path)?.unwrap_or(0.0),
                s0_perc: as_f64(&aj, "s0_perc", &cfg_path)?.unwrap_or(0.05),
                group: aj_group.clone(),
                group_labels: None,
                struc_zero: as_bool(&aj, "struc_zero").unwrap_or(false),
                neg_lb: as_bool(&aj, "neg_lb").unwrap_or(false),
                alpha: as_f64(&aj, "alpha", &cfg_path)?.unwrap_or(0.05),
                global: as_bool(&aj, "global").unwrap_or(false),
                pairwise: as_bool(&aj, "pairwise").unwrap_or(false),
                ..AncombcConfig::default()
            };

            let res =
                ancombc2_run_named(&counts, &x, group_index.as_deref(), &cfg, &taxa, &samples)
                    .map_err(|e| e.to_string())?;
            let c = &res.core;

            // The coefficient under comparison, by name. The group contrast is
            // preferred -- `<variable><level>` for the level that is not the
            // reference -- and only if the design has no group does this fall
            // back to the first non-intercept column. A *positional* fallback
            // would be wrong on a design with a covariate: one arm would take
            // the group and the other the covariate, and every sign would
            // disagree by construction while both arms reported success.
            let want = match aj.get("coefficient").cloned() {
                Some(c) => c,
                None => {
                    let group_contrast = aj_group
                        .as_ref()
                        .zip(aj_group_levels.as_ref())
                        .and_then(|(g, lv)| (lv.len() >= 2).then(|| format!("{g}{}", lv[1])))
                        .filter(|c| colnames.contains(c));
                    match group_contrast {
                        Some(c) => c,
                        None => colnames
                            .iter()
                            .find(|n| n.as_str() != "(Intercept)")
                            .cloned()
                            .unwrap_or_else(|| "(Intercept)".to_string()),
                    }
                }
            };
            let take = |v: &Vec<f64>| -> Result<Vec<f64>, String> {
                c.coefficient(&want, v).ok_or_else(|| {
                    format!(
                        "cannot extract `{want}`: fix_eff = {:?}, {} taxa, {} values",
                        c.fix_eff,
                        c.taxa.len(),
                        v.len()
                    )
                })
            };
            let j = |v: &Vec<f64>| -> String {
                if v.is_empty() {
                    "[]".to_string()
                } else {
                    let parts: Vec<String> = v
                        .iter()
                        .map(|x| {
                            if x.is_finite() {
                                format!("{x:.17e}")
                            } else {
                                "null".to_string()
                            }
                        })
                        .collect();
                    format!("[{}]", parts.join(","))
                }
            };
            let jb = |v: &[bool]| -> String {
                if v.is_empty() {
                    "[]".to_string()
                } else {
                    let parts: Vec<&str> = v
                        .iter()
                        .map(|b| if *b { "true" } else { "false" })
                        .collect();
                    format!("[{}]", parts.join(","))
                }
            };
            let jn = |v: Option<f64>| match v {
                Some(x) if x.is_finite() => format!("{x:.17e}"),
                _ => "null".to_string(),
            };
            let names: Vec<String> = c.taxa.iter().map(|i| taxa[*i].clone()).collect();
            let jnames = |v: &[String]| -> String {
                let parts: Vec<String> = v
                    .iter()
                    .map(|s| format!("\"{}\"", s.replace('\\', "\\\\").replace('"', "\\\"")))
                    .collect();
                format!("[{}]", parts.join(","))
            };

            // `passed_ss` and `diff_robust` are absent when the sensitivity
            // analysis did not run, which is `null` rather than `false` -- and
            // when it did, they are `n_taxa * p` like every other per-taxon
            // vector, so they need the same column extraction as `diff_abn`.
            let flag_col = |v: &Option<Vec<bool>>| -> String {
                match v {
                    None => "null".to_string(),
                    Some(flags) => match c.coefficient(&want, flags) {
                        Some(col) => jb(&col),
                        None => "null".to_string(),
                    },
                }
            };
            let jss = flag_col(&res.passed_ss);
            let jdr = flag_col(&res.diff_robust);

            let elapsed = t0.elapsed().as_secs_f64();
            Ok(format!(
                "{{\"arm\":\"rust\",\"dataset\":\"{name}\",\"n_taxa_in\":{},\"n_samp_in\":{n_samp},\"n_taxa_retained\":{},\"coefficient\":\"{want}\",\"n_coefficients\":{},\"elapsed_seconds\":{},\"taxa\":{},\"beta\":{},\"se\":{},\"p\":{},\"q\":{},\"diff_abn\":{},\"passed_ss\":{jss},\"diff_robust\":{jdr},\"error\":null}}",
                counts.n_taxa,
                names.len(),
                columns.len(),
                jn(Some(elapsed)),
                jnames(&names),
                j(&take(&c.beta)?),
                j(&take(&c.se)?),
                j(&take(&c.p)?),
                j(&take(&c.q)?),
                jb(&take_flags(&c.diff_abn, &want, c)?),
            ))
        })();
        match row {
            Ok(text) => {
                writeln!(w, "{text}").map_err(|e| format!("writing {}: {e}", out.display()))?;
                eprintln!(
                    "rust arm on {name}: ok ({:.1}s)",
                    t0.elapsed().as_secs_f64()
                );
            }
            Err(e) => {
                let msg = e.replace('\n', " ").replace('"', "'");
                writeln!(
                    w,
                    "{{\"arm\":\"rust\",\"dataset\":\"{name}\",\"error\":\"{msg}\"}}"
                )
                .map_err(|e| format!("writing {}: {e}", out.display()))?;
                eprintln!("rust arm on {name}: FAILED -- {e}");
            }
        }
        wrote += 1;
    }
    writeln!(w).ok();
    w.flush()
        .map_err(|e| format!("flushing {}: {e}", out.display()))?;
    println!("wrote {wrote} row(s) to {}", out.display());
    Ok(())
}

/// Extract one coefficient's flags. `diff_abn` is `n_taxa * p`, so it needs the
/// same treatment as the numeric per-taxon vectors.
fn take_flags(
    values: &[bool],
    name: &str,
    c: &ancombc2_core::pipeline::CoreOutput,
) -> Result<Vec<bool>, String> {
    c.coefficient(name, values).ok_or_else(|| {
        format!(
            "cannot extract `{name}`: fix_eff = {:?}, {} taxa, {} flags",
            c.fix_eff,
            c.taxa.len(),
            values.len()
        )
    })
}

fn cmd_cells(grid_path: &Path) -> Result<(), String> {
    let g = load(grid_path)?;
    let cells = g.cells();
    println!(
        "grid {}: {} cells, {} replicates in total, seed {}",
        g.name,
        cells.len(),
        g.total_reps(),
        g.seed
    );
    if !g.note.is_empty() {
        println!("note: {}", g.note);
    }
    if !g.blocks.is_empty() {
        println!("blocks:");
        // Per-block spans are derived from each block's own expanded levels, so
        // the printed span cannot drift from the cell indices the rest of the
        // harness uses.
        let mut start = 0usize;
        for (i, b) in g.blocks.iter().enumerate() {
            let levels = b.levels(&g);
            let end = start + levels.len();
            println!(
                "  [{i}] {:<28} cells {}..{}  ({} cells x {} reps)",
                b.name,
                start,
                end - 1,
                levels.len(),
                b.reps.unwrap_or(g.reps)
            );
            start = end;
        }
    }
    println!("  cell  taxa  samp   DA%  LFC  zinf%   lib_mean  lib_cv  confound   reps  sens");
    for c in &cells {
        let a = g.analysis_for(c);
        println!(
            "{:>5}  {:>4}  {:>4}  {:>4}  {:>3}  {:>4}  {:>9}  {:>6}  {:>8}  {:>5}  {}",
            c.index,
            c.n_taxa,
            c.n_samp,
            (c.da_proportion * 100.0).round() as i64,
            c.log_fc,
            (c.zero_inflation * 100.0).round() as i64,
            c.lib_mean,
            c.lib_cv,
            if c.confound { "yes" } else { "no" },
            a.reps,
            a.sensitivity,
        );
    }
    Ok(())
}
