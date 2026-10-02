//! Running one arm of the grid.
//!
//! The Rust arm runs in-process: generate, analyse with
//! [`ancombc2_core::ancombc2_run`], evaluate. The R arm is driven by
//! `scripts/sim_r.R` over the same generated tables, so both arms see identical
//! data. Neither arm's result is trusted over the other's -- the comparison in
//! [`crate::summarise`] is what decides whether they agree.

use ancombc2_core::config::{AdjustMethod, AncombcConfig};
use ancombc2_core::matrix::Matrix;
use ancombc2_core::pipeline::ancombc2_run_named;
use rayon::prelude::*;

use crate::grid::{Cell, Grid};
use crate::metrics::{evaluate, RepMetrics, Reported};
use crate::simulate::{simulate, Replicate, SimControl};

pub fn sim_control(grid: &Grid) -> SimControl {
    SimControl {
        abundance_sd: grid.abundance_sd,
        dispersion: grid.dispersion,
        // The confounder is sized so a DA taxon is roughly 30 points sparser in
        // the control group on top of the cell's own zero inflation. That is
        // enough to decorrelate the sampling fraction from the group without
        // emptying the control group entirely.
        confounder: 0.3,
    }
}

/// The analysis config for a cell, from the grid's and the cell block's
/// analysis-level factors.
pub fn analysis_config(grid: &Grid, cell: &Cell) -> AncombcConfig {
    let analysis = grid.analysis_for(cell);
    let group2_name = group_name(cell.n_taxa, cell.n_samp);
    let mut cfg = AncombcConfig {
        // The design's column names. Without them the core falls back to `V0`,
        // `V1`, ... and a caller that wants "the group coefficient" cannot ask
        // for it by name -- which is how a flat `n_taxa * p` array gets indexed
        // as one value per taxon and every taxon ends up wearing another
        // taxon's intercept. The CLI sets these from the parsed formula, so the
        // sim does the same from its own design.
        fix_eff: vec![INTERCEPT_COEF.to_string(), group2_name.to_string()],
        group: Some(GROUP_VARIABLE.to_string()),
        group_labels: Some(vec![
            format!("{GROUP_VARIABLE}1"),
            format!("{GROUP_VARIABLE}2"),
        ]),
        pseudo: grid.pseudo,
        prv_cut: analysis.prevalence,
        lib_cut: analysis.lib_size,
        p_adj_method: AdjustMethod::parse(&analysis.p_adjust).unwrap_or(AdjustMethod::Bh),
        s0_perc: analysis.s0_perc,
        // Read from the grid rather than left at the `AncombcConfig` default.
        // They were previously defaults-on-this-side / literals-in-sim_r.R, which
        // meant the two arms ran different analyses; see `grid.rs`.
        struc_zero: analysis.struc_zero,
        neg_lb: analysis.neg_lb,
        ..AncombcConfig::default()
    };
    match analysis.sensitivity.as_str() {
        // The conservative analysis sweeps the pseudo-count and reports
        // `passed_ss`; the non-conservative one re-estimates the sampling
        // fractions once and refits the inference 50 times, which the grid
        // enables with `pseudo_sens` plus `conservative = false`.
        "conservative" => {
            cfg.pseudo_sens = true;
            cfg.conservative = true;
        }
        "nonconservative" => {
            cfg.pseudo_sens = true;
            cfg.conservative = false;
        }
        _ => cfg.pseudo_sens = false,
    }
    cfg.alpha = grid.alpha;
    cfg
}

/// The design: an intercept and a group term, group 1 as the reference.
///
/// The group term is the estimand, so its coefficient is the reported `beta`.
fn design(group: &[usize]) -> Matrix {
    let n = group.len();
    let mut m = Matrix::zeros(n, 2);
    for j in 0..n {
        m.set(j, 0, 1.0);
        m.set(j, 1, if group[j] == 1 { 1.0 } else { 0.0 });
    }
    m
}

/// Analyse one replicate. Errors are returned, not swallowed: an arm that
/// cannot analyse a rep records the failure with the rep's identity so the
/// denominator of every pooled metric stays honest.
/// The group variable, as the design and the metadata name it.
pub const GROUP_VARIABLE: &str = "group";

/// The intercept column's name, matching R's `model.matrix`.
pub const INTERCEPT_COEF: &str = "(Intercept)";

/// The name of the group coefficient in the design [`design`] builds.
///
/// The design is `~ group` with the control as the reference level, so the
/// treatment contrast is the second column, and R names it after the *treated*
/// level -- `group2`, not `group`.
fn group_name(_n_taxa: usize, _n_samp: usize) -> String {
    format!("{GROUP_VARIABLE}2")
}

pub fn run_replicate(grid: &Grid, cell: &Cell, rep: usize, arm: &str) -> RepMetrics {
    let ctl = sim_control(grid);
    let r = simulate(cell, grid.seed, rep, &ctl);
    let cfg = analysis_config(grid, cell);
    let x = design(&r.group);
    // The *named* entry point, not `ancombc2_run`. `ancombc2_run` passes empty
    // taxon and sample name vectors, and the analysis uses the sample names to
    // align the design's rows with the count table's columns; with no names the
    // group term is estimated against the wrong rows and every coefficient
    // collapses towards zero. The names are also what let a result name a taxon
    // rather than an index.
    match ancombc2_run_named(
        &r.counts,
        &x,
        Some(&r.group),
        &cfg,
        &r.taxon_names,
        &r.sample_names,
    ) {
        Ok(res) => {
            let out = &res.core;
            // The group coefficient, by *name*. `beta`, `se`, `p` and `q` are
            // `n_taxa * p` row-major, so passing them straight through as one
            // value per taxon reads the intercept of taxon `t` and the group
            // coefficient of taxon `t-1` as if they were every taxon's effect.
            // The result is plausible, sign-bearing, correlated with the truth,
            // and wrong -- and it is exactly what a per-taxon metric cannot
            // detect on its own, because a shifted intercept looks like noise.
            let coef = group_name(cell.n_taxa, cell.n_samp);
            let columns = [
                ("beta", out.beta.as_slice()),
                ("se", out.se.as_slice()),
                ("p", out.p.as_slice()),
                ("q", out.q.as_slice()),
            ];
            let mut extracted: Vec<Vec<f64>> = Vec::with_capacity(columns.len());
            for (label, values) in &columns {
                match out.coefficient(&coef, values) {
                    Some(v) => extracted.push(v),
                    None => {
                        return RepMetrics::failed(
                            &grid.name,
                            cell.index,
                            rep,
                            arm,
                            cell,
                            format!(
                                "cannot extract the `{coef}` column of {label}: fix_eff = {:?}, \
                                 {} taxa, {} values",
                                out.fix_eff,
                                out.taxa.len(),
                                values.len()
                            ),
                        )
                    }
                }
            }
            let (beta, se, p, q) = (
                extracted[0].clone(),
                extracted[1].clone(),
                extracted[2].clone(),
                extracted[3].clone(),
            );
            // `diff_abn` is per coefficient too, so it needs the same treatment.
            let called = match out.coefficient(&coef, &out.diff_abn) {
                Some(v) => v,
                None => {
                    return RepMetrics::failed(
                        &grid.name,
                        cell.index,
                        rep,
                        arm,
                        cell,
                        format!(
                            "cannot extract the `{coef}` column of diff_abn: fix_eff = {:?}, \
                             {} taxa, {} flags",
                            out.fix_eff,
                            out.taxa.len(),
                            out.diff_abn.len()
                        ),
                    )
                }
            };
            let reported = match Reported::new(&beta, &se, &p, &q, &called) {
                Ok(v) => v,
                Err(e) => {
                    return RepMetrics::failed(
                        &grid.name,
                        cell.index,
                        rep,
                        arm,
                        cell,
                        format!("inconsistent result vectors: {e}"),
                    )
                }
            };
            evaluate(&grid.name, cell, rep, arm, &r, &out.taxa, &reported)
        }
        Err(e) => RepMetrics::failed(&grid.name, cell.index, rep, arm, cell, e.to_string()),
    }
}

/// Run every rep of one cell, at the cell block's rep count.
pub fn run_cell(grid: &Grid, cell: &Cell, arm: &str) -> Vec<RepMetrics> {
    (0..grid.reps_for(cell))
        .into_par_iter()
        .map(|rep| run_replicate(grid, cell, rep, arm))
        .collect()
}

/// Run every cell of the grid.
pub fn run_grid(grid: &Grid, arm: &str, progress: bool) -> Vec<RepMetrics> {
    let cells = grid.cells();
    let mut out = Vec::with_capacity(grid.total_reps());
    for cell in &cells {
        let rows = run_cell(grid, cell, arm);
        if progress {
            let failed = rows.iter().filter(|r| r.error.is_some()).count();
            eprintln!(
                "cell {:>4}/{} n_taxa={} n_samp={} da={} zinf={} conf={} reps={} failed={}",
                cell.index + 1,
                cells.len(),
                cell.n_taxa,
                cell.n_samp,
                cell.da_proportion,
                cell.zero_inflation,
                cell.confound,
                rows.len(),
                failed
            );
        }
        out.extend(rows);
    }
    out
}

// ---------------------------------------------------------------------------
// Writing a rep out, for the R arm.
// ---------------------------------------------------------------------------

/// Write a replicate's counts, metadata, and truth for the R arm to read.
///
/// The R arm must analyse *exactly* the data the Rust arm analysed, so the
/// generator's output is written verbatim rather than regenerated: a shared
/// seed is not a shared table if the generator is ever changed.
pub fn write_replicate(
    dir: &std::path::Path,
    r: &Replicate,
    cell: &Cell,
    rep: usize,
) -> std::io::Result<()> {
    std::fs::create_dir_all(dir)?;
    let n_taxa = r.counts.n_taxa;
    let n_samp = r.counts.n_samp;

    // R's `write.table(x, row.names = TRUE)` layout: an *unnamed* first header
    // field, then one field per sample, and one field per row for the taxon name
    // plus the counts. The project's own reader requires exactly this shape, and
    // a named first column is rejected as ragged -- so writing a named one would
    // produce a file the Rust CLI cannot read, while `read.delim(row.names = 1)`
    // accepts both. The unnamed form is what the oracle and the benchmark
    // harness both read.
    let mut counts = String::with_capacity(n_taxa * n_samp * 4);
    for (j, s) in r.sample_names.iter().enumerate() {
        if j > 0 {
            counts.push('\t');
        }
        counts.push_str(s);
    }
    counts.push('\n');
    for i in 0..n_taxa {
        counts.push_str(&r.taxon_names[i]);
        for j in 0..n_samp {
            counts.push('\t');
            // Counts are integers; writing them as such is what R's `read.table`
            // needs to keep them out of a numeric-precision discussion.
            counts.push_str(&format!("{}", r.counts.data[i * n_samp + j] as i64));
        }
        counts.push('\n');
    }
    std::fs::write(dir.join("counts.tsv"), counts)?;

    // The same unnamed-first-field layout as the counts: the header names only
    // the variables, and the row labels are the sample names. This is what both
    // the project's metadata reader and the benchmark datasets use, and a named
    // first column is read as a ragged row.
    let mut meta = String::from("group\n");
    for j in 0..n_samp {
        meta.push_str(&format!("{}\t{}\n", r.sample_names[j], r.group[j] + 1));
    }
    std::fs::write(dir.join("meta.tsv"), meta)?;

    // `is_da` is the generator's label and `has_effect` is whether the expected
    // abundance really differs. The two differ in a confounded cell, and the
    // metrics count over `has_effect`, so both columns are written: a truth file
    // with only the label would make the R arm score the confounded taxa as true
    // positives and disagree with the Rust arm for a reason that is in the file
    // rather than in the code.
    let mut truth = String::from("taxon\tis_da\thas_effect\tconfounded\tdirection\tlog_fc\tzero_fraction\tzero_fraction_control\tzero_fraction_treated\n");
    for t in &r.truth {
        truth.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            t.name,
            if t.is_da { 1 } else { 0 },
            if t.has_effect { 1 } else { 0 },
            if t.confounded { 1 } else { 0 },
            t.direction,
            t.log_fc,
            t.zero_fraction,
            t.zero_fraction_control,
            t.zero_fraction_treated
        ));
    }
    std::fs::write(dir.join("truth.tsv"), truth)?;

    let cfg = serde_json::json!({
        "cell": cell.index,
        "rep": rep,
        "n_taxa": n_taxa,
        "n_samp": n_samp,
        "da_proportion": cell.da_proportion,
        "log_fc": cell.log_fc,
        "zero_inflation": cell.zero_inflation,
        "lib_mean": cell.lib_mean,
        "lib_cv": cell.lib_cv,
        "confound": cell.confound,
    });
    std::fs::write(dir.join("replicate.json"), format!("{cfg}\n"))?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::grid::Grid;

    fn grid() -> Grid {
        Grid {
            name: "t".into(),
            seed: 4242,
            n_taxa: vec![60],
            n_samp: vec![10],
            da_proportion: vec![0.5],
            log_fc: vec![1.0],
            zero_inflation: vec![0.0],
            lib_mean: vec![1e4],
            lib_cv: vec![0.0],
            confound: vec![false],
            dispersion: 1.0,
            abundance_sd: 0.5,
            struc_zero: true,
            neg_lb: true,
            s0_perc: 0.05,
            reps: 2,
            alpha: 0.05,
            pseudo: 0.5,
            sensitivity: "none".into(),
            prevalence: 0.0,
            lib_size: 0.0,
            p_adjust: "BH".into(),
            note: String::new(),
            blocks: Vec::new(),
        }
    }

    #[test]
    fn a_replicate_analyses_and_reports_its_truth() {
        let g = grid();
        let cells = g.cells();
        let m = run_replicate(&g, &cells[0], 0, "rust");
        assert!(m.error.is_none(), "{:?}", m.error);
        assert_eq!(m.n_retained, m.beta.len());
        assert_eq!(m.beta.len(), m.se.len());
        assert_eq!(m.diff_abn.len(), m.is_da.len());
        assert!(m.n_retained > 0, "some taxa must survive filtering");
        assert!(m.lfc_rmse.map(f64::is_finite).unwrap_or(false));
    }

    #[test]
    fn the_analysis_is_reproducible() {
        let g = grid();
        let cells = g.cells();
        let a = run_replicate(&g, &cells[0], 1, "rust");
        let b = run_replicate(&g, &cells[0], 1, "rust");
        assert_eq!(a.beta, b.beta);
        assert_eq!(a.p, b.p);
        assert_eq!(a.diff_abn, b.diff_abn);
    }

    #[test]
    fn the_design_puts_the_effect_on_the_second_column() {
        let x = design(&[0, 0, 1, 1]);
        assert_eq!(x.rows, 4);
        assert_eq!(x.cols, 2);
        assert_eq!(x.get(0, 1), 0.0);
        assert_eq!(x.get(2, 1), 1.0);
        assert_eq!(x.get(3, 0), 1.0);
    }

    #[test]
    fn the_group_term_recovers_the_planted_effect_direction() {
        // A large, clean effect on a well-sampled cell must be recovered with
        // the right sign. This is the analysis-side end-to-end sanity check:
        // it fails if the design, the group coding, or the reported coefficient
        // is wrong, none of which the parity fixtures would catch on simulated
        // data.
        let mut g = grid();
        g.n_taxa = vec![80];
        g.da_proportion = vec![0.5];
        g.log_fc = vec![2.0];
        let cells = g.cells();
        let m = run_replicate(&g, &cells[0], 0, "rust");
        assert!(m.error.is_none());
        let ups: Vec<f64> = m
            .beta
            .iter()
            .zip(m.truth_log_fc.iter())
            .filter(|(_, t)| **t > 0.0)
            .map(|(b, _)| *b)
            .collect();
        let downs: Vec<f64> = m
            .beta
            .iter()
            .zip(m.truth_log_fc.iter())
            .filter(|(_, t)| **t < 0.0)
            .map(|(b, _)| *b)
            .collect();
        assert!(
            !ups.is_empty() && !downs.is_empty(),
            "both directions planted"
        );
        let mean_up = ups.iter().sum::<f64>() / ups.len() as f64;
        let mean_down = downs.iter().sum::<f64>() / downs.len() as f64;
        assert!(
            mean_up > 0.0,
            "up-regulated taxa must be positive: {mean_up}"
        );
        assert!(
            mean_down < 0.0,
            "down-regulated taxa must be negative: {mean_down}"
        );
    }

    #[test]
    fn a_written_replicate_round_trips_through_disk() {
        let g = grid();
        let cells = g.cells();
        let r = simulate(&cells[0], g.seed, 0, &sim_control(&g));
        let dir = std::env::temp_dir().join(format!("simtest-{}", std::process::id()));
        write_replicate(&dir, &r, &cells[0], 0).unwrap();
        let counts = std::fs::read_to_string(dir.join("counts.tsv")).unwrap();
        let lines: Vec<&str> = counts.trim().lines().collect();
        assert_eq!(lines.len(), r.counts.n_taxa + 1);
        // One field per sample in the header, one more per row for the taxon.
        assert_eq!(lines[0].split('\t').count(), r.counts.n_samp);
        assert_eq!(lines[1].split('\t').count(), r.counts.n_samp + 1);
        assert_eq!(lines[1].split('\t').next(), Some(r.taxon_names[0].as_str()));
        let meta = std::fs::read_to_string(dir.join("meta.tsv")).unwrap();
        let meta_lines: Vec<&str> = meta.trim().lines().collect();
        assert_eq!(
            meta_lines.len(),
            r.counts.n_samp + 1,
            "a header plus one row per sample"
        );
        assert_eq!(
            meta_lines[0], "group",
            "only the group variable, no name column"
        );
        assert_eq!(
            meta_lines[1].split('\t').count(),
            2,
            "a label and one value"
        );
        let truth = std::fs::read_to_string(dir.join("truth.tsv")).unwrap();
        let truth_lines: Vec<&str> = truth.trim().lines().collect();
        assert_eq!(truth_lines.len(), r.counts.n_taxa + 1);
        // The label, the effect, and the confounded flag are all separate
        // columns, and the R arm needs all three.
        let header: Vec<&str> = truth_lines[0].split('\t').collect();
        for want in ["is_da", "has_effect", "confounded", "log_fc"] {
            assert!(
                header.contains(&want),
                "truth.tsv is missing the `{want}` column"
            );
        }
        assert!(dir.join("replicate.json").exists());
        std::fs::remove_dir_all(&dir).ok();
    }
}
