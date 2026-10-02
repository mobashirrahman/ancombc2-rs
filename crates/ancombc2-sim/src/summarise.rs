//! Pooling reps into per-cell summaries, and the Rust-vs-R acceptance rule.
//!
//! # The acceptance rule
//!
//! `PLAN.md` §5.3: *"for every cell, empirical FDR of Rust within Monte-Carlo
//! error of R's (3 SE of the binomial), and power within 3 SE."*
//!
//! Three things are done here that the sentence leaves open, each recorded in
//! the output so a reader can disagree with it:
//!
//! 1. **Pooling.** A cell's FDR and power are the ratios of *pooled* counts
//!    over its reps, not the mean of the per-rep ratios. The per-rep ratio is
//!    undefined whenever a rep calls nothing, and averaging ratios also
//!    over-weights a rep that made one call. The pooled ratio is the
//!    quantity whose binomial standard error the rule is written against.
//! 2. **The standard error.** `sqrt(p (1 - p) / n)` on the *pooled* counts.
//!    For FDR, `n` is the total number of calls and `p` the pooled false-call
//!    fraction; for power, `n` is the total number of DA taxa and `p` the
//!    pooled hit rate. Reps within a cell are independent, so pooling is what
//!    makes the binomial formula applicable at all.
//! 3. **The verdict.** `agree` when the gap is within `3 * sqrt(se_r^2 +
//!    se_o^2)`, the standard error of a difference of two independent
//!    proportions. `rust_only` / `r_only` when one arm is outside its own
//!    3 SE band but the two are not outside each other -- a one-armed
//!    excursion. `divergent` otherwise.
//!
//! A `divergent` cell is a *finding*, not a test failure: `PLAN.md` says a
//! genuine divergence is publishable.

use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

use crate::grid::Grid;
use crate::metrics::{binomial_se, RepMetrics};

/// A ratio pooled over a cell's reps, with the SE the acceptance rule uses.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Pooled {
    pub numerator: f64,
    pub denominator: f64,
    pub estimate: f64,
    pub se: f64,
    /// 3 SE, the width of the acceptance band.
    pub band: f64,
    pub low: f64,
    pub high: f64,
    /// Replicates that contributed.
    pub reps: usize,
    /// Reps whose analysis failed, and which therefore did *not* contribute.
    pub failed_reps: usize,
}

impl Pooled {
    fn ratio(num: f64, den: f64, reps: usize, failed: usize) -> Self {
        if den <= 0.0 {
            return Self {
                numerator: num,
                denominator: den,
                estimate: f64::NAN,
                se: f64::NAN,
                band: f64::NAN,
                low: f64::NAN,
                high: f64::NAN,
                reps,
                failed_reps: failed,
            };
        }
        let estimate = num / den;
        let se = binomial_se(estimate, den.round() as usize);
        let band = 3.0 * se;
        Self {
            numerator: num,
            denominator: den,
            estimate,
            se,
            band,
            // A binomial rate is bounded; the band is clipped to [0, 1] so a
            // small-n cell cannot have an acceptance band that leaves the
            // parameter space.
            low: (estimate - band).max(0.0),
            high: (estimate + band).min(1.0),
            reps,
            failed_reps: failed,
        }
    }

    /// The mean of the finite values, or `NaN` when there are none.
    fn mean(values: &[f64]) -> f64 {
        let finite: Vec<f64> = values.iter().copied().filter(|v| v.is_finite()).collect();
        if finite.is_empty() {
            f64::NAN
        } else {
            finite.iter().sum::<f64>() / finite.len() as f64
        }
    }
}

/// A mean over reps where an undefined per-rep value poisons the mean.
///
/// A rate like sign concordance is undefined for a rep that retained no DA
/// taxa, and averaging over the reps that do have a value would quietly
/// overstate the cell's support. The accumulator below sums the per-rep values,
/// so a single `NaN` propagates and the cell reports `NaN` -- the honest
/// answer for a cell with nothing to say.
#[derive(Default)]
struct MaybeMean {
    sum: f64,
    n: usize,
}

impl MaybeMean {
    fn new() -> Self {
        Self { sum: 0.0, n: 0 }
    }
    fn add(&mut self, v: f64) {
        self.sum += v;
        self.n += 1;
    }
    fn get(&self) -> f64 {
        if self.n == 0 {
            f64::NAN
        } else {
            self.sum / self.n as f64
        }
    }
}

/// How the two arms compare on one metric.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Verdict {
    /// Within three SE of the difference.
    Agree,
    /// One arm outside its own 3 SE band, the two consistent with each other.
    RustOnly,
    ROnly,
    /// The arms disagree beyond the combined 3 SE.
    Divergent,
    /// Neither arm has a value for this metric -- a cell where nothing was
    /// called has no empirical FDR on either side, and two absences are not a
    /// disagreement.
    ///
    /// Reported separately because "we cannot tell" and "they disagree" are
    /// different findings, and a grid that scored them the same would either
    /// bury the real divergences or cry wolf on every sparse cell.
    Inconclusive,
}

/// One metric's comparison.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Comparison {
    pub metric: String,
    pub rust: f64,
    pub r: f64,
    pub rust_se: f64,
    pub r_se: f64,
    pub difference: f64,
    /// `3 * sqrt(se_rust^2 + se_r^2)`.
    pub tolerance: f64,
    pub verdict: Verdict,
}

/// A cell's full summary, both arms.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CellSummary {
    pub cell: usize,
    pub n_taxa: usize,
    pub n_samp: usize,
    pub da_proportion: f64,
    pub log_fc: f64,
    pub zero_inflation: f64,
    pub lib_mean: f64,
    pub lib_cv: f64,
    pub confound: bool,
    pub rust: ArmSummary,
    pub r: Option<ArmSummary>,
    pub comparisons: Vec<Comparison>,
    /// True when every comparison reached a verdict and none was `divergent`.
    pub agrees: bool,
    /// True when at least one comparison had no value on one or both sides.
    pub inconclusive: bool,
}

/// One arm's per-cell summary.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ArmSummary {
    pub arm: String,
    pub reps: usize,
    pub failed_reps: usize,
    pub n_retained: f64,
    pub fdr: Pooled,
    pub power: Pooled,
    pub sign_concordance: f64,
    pub null_sign_rate: f64,
    /// Pooled over the cell's successful replicates. `None` only if every
    /// replicate of the cell failed, which is the case for the extreme-sparsity
    /// cells where no taxon survives the structural-zero screen at all.
    pub lfc_bias: Option<f64>,
    pub lfc_bias_se: f64,
    pub lfc_rmse: Option<f64>,
    pub f1: f64,
    pub jaccard: f64,
    /// `mean(se) / sd(beta over reps)`, per taxon then averaged. `NaN` when
    /// there are too few reps for a per-taxon standard deviation, which is the
    /// honest answer for a one-rep cell.
    pub se_calibration: f64,
    /// The same statistic from the R arm, for reference.
    pub mean_se: f64,
    /// Calls among the confounded taxa, pooled over the cell's replicates.
    ///
    /// These taxa carry no compositional effect -- they differ only in detection
    /// rate -- so this is the sampling-fraction correction's false-positive rate
    /// on the grid's own negative control. A correction that works keeps it near
    /// zero; one that does not work produces the same number as an unadjusted
    /// method's.
    pub confounded_false_positive: Pooled,
    /// How many confounded taxa the cell retained in total, so a rate over one
    /// or two taxa is visible as such.
    pub n_confounded_retained: f64,
}

/// A grid's full summary.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GridSummary {
    pub grid: String,
    pub note: String,
    pub seed: u64,
    pub alpha: f64,
    pub cells: Vec<CellSummary>,
    pub n_cells: usize,
    pub n_divergent: usize,
    /// Cells where at least one metric had no value on one or both sides.
    pub n_inconclusive: usize,
    /// Every divergent cell's description, so a finding is quotable without
    /// re-deriving the cell from its index.
    pub findings: Vec<String>,
    /// A run with no R rows is not a comparison, and the summary says so rather
    /// than reporting every cell as agreeing.
    pub compared_against_r: bool,
}

#[derive(Default)]
struct Acc {
    reps: usize,
    failed: usize,
    retained: f64,
    calls: f64,
    false_calls: f64,
    da: f64,
    hits: f64,
    sign_conc: MaybeMean,
    null_sign: MaybeMean,
    confounded_taxa: f64,
    confounded_called: f64,
    /// Per-taxon `beta` across reps, for the SE calibration ratio.
    beta: BTreeMap<String, Vec<f64>>,
    se: BTreeMap<String, Vec<f64>>,
}

/// Pool one arm's replicates for one cell.
///
/// `cell` is not optional. Pooling every row an arm produced into every cell
/// produces a table of identical numbers -- each cell reporting the whole run's
/// aggregate -- which reads as a suspiciously clean result rather than as the
/// bug it is.
fn summarise_arm(rows: &[RepMetrics], cell: usize, arm: &str) -> ArmSummary {
    let mut a = Acc {
        sign_conc: MaybeMean::new(),
        null_sign: MaybeMean::new(),
        ..Acc::default()
    };
    let mut sum_bias = 0.0;
    let mut sum_rmse = 0.0;
    let mut sum_f1 = 0.0;
    let mut sum_jac = 0.0;
    let mine: Vec<&RepMetrics> = rows.iter().filter(|r| r.cell == cell).collect();
    let keep: Vec<&RepMetrics> = mine
        .iter()
        .copied()
        .filter(|r| r.arm == arm)
        .filter(|r| r.error.is_none())
        .collect();
    a.reps = keep.len();
    a.failed = mine
        .iter()
        .copied()
        .filter(|r| r.arm == arm && r.error.is_some())
        .count();
    for r in &keep {
        a.retained += r.n_retained as f64;
        a.calls += r.n_diff_abn as f64;
        a.da += r.n_da_retained as f64;
        a.hits += (0..r.n_retained)
            .filter(|k| r.is_da[*k] && r.diff_abn[*k])
            .count() as f64;
        a.false_calls += (0..r.n_retained)
            .filter(|k| !r.is_da[*k] && r.diff_abn[*k])
            .count() as f64;
        a.sign_conc.add(r.sign_concordance.unwrap_or(f64::NAN));
        a.confounded_taxa += r.n_confounded_retained as f64;
        a.confounded_called += r
            .confounded_false_positive_rate
            .map(|v| v * r.n_confounded_retained as f64)
            .unwrap_or(0.0);
        a.null_sign.add(r.null_sign_rate.unwrap_or(f64::NAN));
        // A failed replicate contributes nothing: it has no fold change. It was
        // previously `NaN`, which is also excluded by the `is_finite` filter on
        // the pooled mean, so the pooled value is unchanged -- but `NaN` did not
        // survive the round trip, and `null` does.
        if let Some(v) = r.lfc_bias {
            sum_bias += v;
        }
        if let Some(v) = r.lfc_rmse {
            sum_rmse += v;
        }
        sum_f1 += r.f1.unwrap_or(f64::NAN);
        sum_jac += r.jaccard.unwrap_or(f64::NAN);
        for k in 0..r.n_retained {
            a.beta
                .entry(r.names[k].clone())
                .or_default()
                .push(r.beta[k]);
            a.se.entry(r.names[k].clone()).or_default().push(r.se[k]);
        }
    }
    let n = a.reps.max(1) as f64;

    // SE calibration: for each taxon, `mean(se) / sd(beta across reps)`. A taxon
    // retained in only one rep has no `sd`, and contributes nothing rather than
    // an infinite ratio.
    let mut ratios = Vec::new();
    for (name, betas) in &a.beta {
        let ses = match a.se.get(name) {
            Some(s) if s.len() == betas.len() && !s.is_empty() => s,
            _ => continue,
        };
        if betas.len() < 2 {
            continue;
        }
        let mean_se = ses.iter().sum::<f64>() / ses.len() as f64;
        let mb = betas.iter().sum::<f64>() / betas.len() as f64;
        let var = betas.iter().map(|b| (b - mb) * (b - mb)).sum::<f64>() / (betas.len() - 1) as f64;
        let sd = var.sqrt();
        if sd > 0.0 && mean_se.is_finite() {
            ratios.push(mean_se / sd);
        }
    }

    ArmSummary {
        arm: arm.to_string(),
        reps: a.reps,
        failed_reps: a.failed,
        n_retained: a.retained / n,
        fdr: Pooled::ratio(a.false_calls, a.calls, a.reps, a.failed),
        power: Pooled::ratio(a.hits, a.da, a.reps, a.failed),
        sign_concordance: a.sign_conc.get(),
        null_sign_rate: a.null_sign.get(),
        lfc_bias: if n > 0.0 { Some(sum_bias / n) } else { None },
        lfc_bias_se: lfc_bias_se(rows, cell, arm),
        lfc_rmse: if n > 0.0 { Some(sum_rmse / n) } else { None },
        f1: sum_f1 / n,
        jaccard: sum_jac / n,
        se_calibration: Pooled::mean(&ratios),
        mean_se: f64::NAN,
        confounded_false_positive: Pooled::ratio(
            a.confounded_called,
            a.confounded_taxa,
            a.reps,
            a.failed,
        ),
        n_confounded_retained: a.confounded_taxa,
    }
}

/// The standard error of the per-rep LFC bias, across reps.
fn lfc_bias_se(rows: &[RepMetrics], cell: usize, arm: &str) -> f64 {
    let v: Vec<f64> = rows
        .iter()
        .filter(|r| r.cell == cell && r.arm == arm && r.error.is_none())
        .filter_map(|r| r.lfc_bias)
        .filter(|v| v.is_finite())
        .collect();
    if v.len() < 2 {
        return f64::NAN;
    }
    let m = v.iter().sum::<f64>() / v.len() as f64;
    (v.iter().map(|x| (x - m) * (x - m)).sum::<f64>() / (v.len() - 1) as f64).sqrt()
}

fn compare(metric: &str, rust: &Pooled, r: &Pooled) -> Comparison {
    let difference = rust.estimate - r.estimate;
    let tolerance = 3.0 * (rust.se * rust.se + r.se * r.se).sqrt();
    let rust_has = rust.estimate.is_finite();
    let r_has = r.estimate.is_finite();
    let verdict = match (rust_has, r_has) {
        // Neither arm had any calls, so neither has a rate. Nothing to compare.
        (false, false) => Verdict::Inconclusive,
        // One arm produced calls and the other produced none. That is not an
        // inability to compare -- it is a difference in what the two
        // implementations *did*, and it is exactly the asymmetry that a
        // "cannot tell" verdict would hide.
        (true, false) | (false, true) => Verdict::Divergent,
        (true, true) if !tolerance.is_finite() => Verdict::Inconclusive,
        (true, true) if difference.abs() <= tolerance => Verdict::Agree,
        (true, true) if rust.estimate < rust.low || rust.estimate > rust.high => Verdict::RustOnly,
        (true, true) if r.estimate < r.low || r.estimate > r.high => Verdict::ROnly,
        (true, true) => Verdict::Divergent,
    };
    Comparison {
        metric: metric.to_string(),
        rust: rust.estimate,
        r: r.estimate,
        rust_se: rust.se,
        r_se: r.se,
        difference,
        tolerance,
        verdict,
    }
}

/// `NaN`-safe equality of the 3 SE rule.
fn compare_scalar(metric: &str, rust: f64, r: f64, n: f64) -> Comparison {
    let se = n.max(1.0).sqrt();
    let difference = rust - r;
    let tolerance = 3.0 * se;
    let verdict = if !difference.is_finite() {
        Verdict::Inconclusive
    } else if difference.abs() <= tolerance {
        Verdict::Agree
    } else {
        Verdict::Divergent
    };
    Comparison {
        metric: metric.to_string(),
        rust,
        r,
        rust_se: se,
        r_se: se,
        difference,
        tolerance,
        verdict,
    }
}

/// Summarise the grid from both arms' per-rep rows.
pub fn summarise(grid: &Grid, rows: &[RepMetrics]) -> GridSummary {
    let cells = grid.cells();
    let has_r = rows.iter().any(|r| r.arm == "r");
    let mut out = Vec::with_capacity(cells.len());
    let mut findings = Vec::new();
    let mut n_divergent = 0usize;
    let mut n_inconclusive = 0usize;

    for cell in &cells {
        let rust = summarise_arm(rows, cell.index, "rust");
        let r = if has_r {
            Some(summarise_arm(rows, cell.index, "r"))
        } else {
            None
        };
        let mut comparisons = Vec::new();
        if let Some(ra) = &r {
            comparisons.push(compare("empirical_fdr", &rust.fdr, &ra.fdr));
            comparisons.push(compare("power", &rust.power, &ra.power));
            // Compared only when the cell actually has a negative control, so a
            // grid without confounded cells reports `not measured` for it rather
            // than a vacuous agreement on two absences.
            if rust.n_confounded_retained > 0.0 || ra.n_confounded_retained > 0.0 {
                comparisons.push(compare(
                    "confounded_false_positive",
                    &rust.confounded_false_positive,
                    &ra.confounded_false_positive,
                ));
            }
            // A bias comparison is a difference of means over reps, so its
            // tolerance is the SE of that difference rather than a binomial.
            let se = (rust.lfc_bias_se * rust.lfc_bias_se + ra.lfc_bias_se * ra.lfc_bias_se).sqrt();
            // A cell where every replicate failed on one arm has no bias to
            // compare. That is `inconclusive`, not agreement: the difference does
            // not exist, and treating a missing value as a zero difference is how
            // a missing value becomes a pass.
            let d = match (rust.lfc_bias, ra.lfc_bias) {
                (Some(x), Some(y)) => x - y,
                _ => f64::NAN,
            };
            comparisons.push(Comparison {
                metric: "lfc_bias".into(),
                rust: rust.lfc_bias.unwrap_or(f64::NAN),
                r: ra.lfc_bias.unwrap_or(f64::NAN),
                rust_se: rust.lfc_bias_se,
                r_se: ra.lfc_bias_se,
                difference: d,
                tolerance: 3.0 * se,
                verdict: if !d.is_finite() || !se.is_finite() {
                    Verdict::Inconclusive
                } else if d.abs() <= 3.0 * se {
                    Verdict::Agree
                } else {
                    Verdict::Divergent
                },
            });
        }
        let cell_divergences: Vec<&Comparison> = comparisons
            .iter()
            .filter(|c| c.verdict == Verdict::Divergent)
            .collect();
        let this_cell_inconclusive = comparisons
            .iter()
            .filter(|c| c.verdict == Verdict::Inconclusive)
            .count();
        // A cell is only *resolved* when every comparison reached a verdict. A
        // cell with nothing called has no FDR to agree on, and calling that
        // "agreement" would let a sparse grid pass by having nothing to say.
        let agrees = cell_divergences.is_empty() && this_cell_inconclusive == 0;
        let inconclusive = this_cell_inconclusive > 0;
        // Three distinct outcomes, and conflating them is how a sparse grid
        // ends up reporting dozens of findings: a cell can have a real
        // divergence, be unresolved, or neither. `agrees` requires full
        // resolution, so it is the strictest of the three.
        if !cell_divergences.is_empty() {
            n_divergent += 1;
            for c in &cell_divergences {
                let detail = if !c.rust.is_finite() && !c.r.is_finite() {
                    format!(
                        "neither arm produced a value (rust {} calls, r {} calls)",
                        c.rust_se, c.r_se
                    )
                } else if !c.rust.is_finite() || !c.r.is_finite() {
                    format!(
                        "one arm produced no value: rust {:.4} vs r {:.4}; the arms made \
                         different calls, so the rate is defined on one side only",
                        c.rust, c.r
                    )
                } else {
                    format!(
                        "gap {:.4} exceeds 3 SE {:.4}",
                        c.difference.abs(),
                        c.tolerance
                    )
                };
                findings.push(format!(
                    "cell {} (taxa {}, samples {}, DA {}%, LFC {}, zero-infl {}%, lib {} CV {}, confound {}): {} rust {:.4} vs r {:.4}, {}",
                    cell.index,
                    cell.n_taxa,
                    cell.n_samp,
                    (cell.da_proportion * 100.0).round(),
                    cell.log_fc,
                    (cell.zero_inflation * 100.0).round(),
                    cell.lib_mean,
                    cell.lib_cv,
                    cell.confound,
                    c.metric,
                    c.rust,
                    c.r,
                    detail
                ));
            }
        }
        if inconclusive {
            n_inconclusive += 1;
        }
        out.push(CellSummary {
            cell: cell.index,
            n_taxa: cell.n_taxa,
            n_samp: cell.n_samp,
            da_proportion: cell.da_proportion,
            log_fc: cell.log_fc,
            zero_inflation: cell.zero_inflation,
            lib_mean: cell.lib_mean,
            lib_cv: cell.lib_cv,
            confound: cell.confound,
            rust,
            r,
            comparisons,
            agrees,
            inconclusive,
        });
    }

    GridSummary {
        grid: grid.name.clone(),
        note: grid.note.clone(),
        seed: grid.seed,
        alpha: grid.alpha,
        n_cells: out.len(),
        n_divergent,
        n_inconclusive,
        findings,
        cells: out,
        compared_against_r: has_r,
    }
}

/// A short human-readable table of the Rust arm, for a terminal or a log.
pub fn table(summary: &GridSummary) -> String {
    let mut s = String::new();
    s.push_str(&format!(
        "grid {} ({} cells, seed {}, alpha {}, compared against R: {})\n",
        summary.grid,
        summary.n_cells,
        summary.seed,
        summary.alpha,
        if summary.compared_against_r {
            "yes"
        } else {
            "no"
        }
    ));
    s.push_str("  cell  taxa  samp   DA%  zinf%  conf |    FDR +-3se     power +-3se   sign  null  bias   se_cal  confFP  calls\n");
    for c in &summary.cells {
        s.push_str(&format!(
            "{:>5}  {:>4}  {:>4}  {:>4}  {:>4}  {:>4} | {}  {}  {}  {}  {}  {}  {}  {}\n",
            c.cell,
            c.n_taxa,
            c.n_samp,
            (c.da_proportion * 100.0).round() as i64,
            (c.zero_inflation * 100.0).round() as i64,
            if c.confound { "y" } else { "n" },
            fmt(&c.rust.fdr),
            fmt(&c.rust.power),
            f3(c.rust.sign_concordance),
            f3(c.rust.null_sign_rate),
            c.rust.lfc_bias.map(f4).unwrap_or_else(|| "-".into()),
            f3(c.rust.se_calibration),
            fmt(&c.rust.confounded_false_positive),
            c.rust.fdr.denominator.round() as i64,
        ));
    }
    if summary.n_inconclusive > 0 {
        s.push_str(&format!(
            "{} cell(s) had at least one metric with no value on one or both arms; \
             those are not counted as agreement\n",
            summary.n_inconclusive
        ));
    }
    if summary.n_divergent > 0 {
        s.push_str(&format!(
            "\n{} divergent cell(s) -- these are findings, see findings[]:\n",
            summary.n_divergent
        ));
        for f in &summary.findings {
            s.push_str(&format!("  {f}\n"));
        }
    }
    s
}

fn f3(v: f64) -> String {
    if v.is_finite() {
        format!("{v:.3}")
    } else {
        "-".into()
    }
}
fn f4(v: f64) -> String {
    if v.is_finite() {
        format!("{v:+.4}")
    } else {
        "-".into()
    }
}
fn fmt(p: &Pooled) -> String {
    if p.estimate.is_finite() {
        format!("{:.3}+-{:.3}", p.estimate, p.band)
    } else {
        "-".into()
    }
}

/// Count how many reps each arm contributed, for a run manifest.
pub fn count_rows(rows: &[RepMetrics], arm: &str) -> (usize, usize) {
    let r: Vec<&RepMetrics> = rows.iter().filter(|x| x.arm == arm).collect();
    (
        r.iter().filter(|x| x.error.is_none()).count(),
        r.iter().filter(|x| x.error.is_some()).count(),
    )
}

/// Convenience for the R-arm script check: does every cell have both arms?
pub fn coverage(grid: &Grid, rows: &[RepMetrics], arm: &str) -> Vec<usize> {
    let n_cells = grid.cells().len();
    (0..n_cells)
        .filter(|c| rows.iter().any(|r| r.cell == *c && r.arm == arm))
        .collect()
}

/// The `mean(se) / sd(beta)` calibration ratio for one arm over the whole grid,
/// as a single number. This is the statistic that would be near 1 for a
/// correctly calibrated sandwich.
///
/// Pooled over every cell, so it is a grid-wide average rather than any cell's
/// calibration; the per-cell value is the one the summary table reports.
pub fn overall_se_calibration(rows: &[RepMetrics], arm: &str) -> f64 {
    let ratios: Vec<f64> = rows
        .iter()
        .filter(|r| r.arm == arm && r.error.is_none())
        .map(|r| r.se.iter().sum::<f64>() / r.se.len().max(1) as f64)
        .collect();
    Pooled::mean(&ratios)
}

/// Not used by the harness directly, but keeps the scalar comparison helper
/// exercised: the same 3 SE rule applied to a non-binomial metric.
pub fn scalar_verdict(rust: f64, r: f64, se: f64) -> Verdict {
    compare_scalar("x", rust, r, se * se).verdict
}

/// A tiny helper the CLI uses when writing a summary for a cell with only one
/// arm: the "agreement" is then vacuous, and the summary must say so.
pub fn single_arm_agrees(comparisons: &[Comparison]) -> bool {
    comparisons.is_empty() || comparisons.iter().all(|c| c.verdict != Verdict::Divergent)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::grid::Cell;
    use crate::metrics::RepMetrics;

    fn cell() -> Cell {
        Cell {
            index: 0,
            n_taxa: 100,
            n_samp: 10,
            da_proportion: 0.5,
            log_fc: 1.0,
            zero_inflation: 0.0,
            lib_mean: 1e4,
            lib_cv: 0.0,
            confound: false,
        }
    }

    fn grid() -> Grid {
        Grid {
            name: "t".into(),
            seed: 1,
            n_taxa: vec![100],
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

    /// A rep whose contribution is controlled by `hits` (DA taxa called) and
    /// `false_calls` (null taxa called), built consistently with the vectors:
    /// the DA hits are the first taxa, the false calls the last, so
    /// `summarise` can re-derive both from the vectors alone. Reps differ by a
    /// deterministic per-rep jitter on `beta`, which is what gives a per-taxon
    /// spread for the SE calibration ratio.
    fn rep(arm: &str, rep: usize, hits: usize, false_calls: usize, da: usize) -> RepMetrics {
        let n = 100;
        let mut m = RepMetrics::failed("t", 0, rep, arm, &cell(), String::new());
        m.error = None;
        m.n_retained = n;
        m.n_da_retained = da;
        m.n_diff_abn = hits + false_calls;
        m.names = (0..n).map(|i| format!("t{i:03}")).collect();
        m.truth_log_fc = (0..n).map(|i| if i < da { 1.0 } else { 0.0 }).collect();
        m.is_da = (0..n).map(|i| i < da).collect();
        m.has_effect = m.is_da.clone();
        m.confounded = vec![false; n];
        let jitter = rep as f64 * 0.1;
        m.beta = (0..n)
            .map(|i| {
                let truth = if i < da { 1.0 } else { 0.0 };
                truth + 0.1 * jitter * (1.0 + i as f64 / n as f64)
            })
            .collect();
        m.se = (0..n).map(|i| 0.5 + 0.01 * i as f64).collect();
        m.p = vec![0.01; n];
        m.q = vec![0.01; n];
        // The first `hits` DA taxa and the last `false_calls` null taxa.
        m.diff_abn = (0..n)
            .map(|i| i < hits || (i >= n - false_calls && i >= da))
            .collect();
        // The scalar fields are what the harness would have written; they are
        // recomputed from the vectors by `evaluate`, and `summarise` uses only
        // the vectors, so these are here to keep the row shape complete.
        m.empirical_fdr = (m.n_diff_abn > 0).then_some(false_calls as f64 / m.n_diff_abn as f64);
        m.power = (da > 0).then_some(hits as f64 / da as f64);
        m.sign_concordance = Some(0.9);
        m.null_sign_rate = Some(0.5);
        m.lfc_bias = Some(0.05 * jitter);
        m.lfc_rmse = Some(0.2);
        m.f1 = Some(0.7);
        m.jaccard = Some(0.6);
        m
    }

    #[test]
    fn pooling_is_on_the_counts_not_the_mean_of_ratios() {
        // One rep makes 1 call of 0 false, the next makes 99 calls of 49
        // false. The pooled FDR is 49/100; the mean of ratios is (0 + 0.495)/2,
        // which would be a different and wrong summary.
        // Rep 0 makes one correct call; rep 1 makes all 50 DA taxa and 49 nulls
        // as well. Pooled: 100 calls of which 49 are false, and 51 hits.
        let rows = vec![rep("rust", 0, 1, 0, 50), rep("rust", 1, 50, 49, 50)];
        let s = summarise_arm(&rows, 0, "rust");
        assert!((s.fdr.estimate - 0.49).abs() < 1e-12, "{}", s.fdr.estimate);
        assert_eq!(s.fdr.numerator, 49.0);
        assert_eq!(s.fdr.denominator, 100.0);
        assert!(
            (s.power.estimate - 0.51).abs() < 1e-12,
            "power is pooled hits/da = 51/100, got {}",
            s.power.estimate
        );
    }

    #[test]
    fn a_cell_with_no_calls_has_no_fdr_rather_than_zero() {
        let rows = vec![rep("rust", 0, 0, 0, 50)];
        let s = summarise_arm(&rows, 0, "rust");
        assert!(
            s.fdr.estimate.is_nan(),
            "an undefined FDR must not read as 0"
        );
        assert!(s.fdr.denominator == 0.0);
        assert!(
            (s.power.estimate - 0.0).abs() < 1e-12,
            "power is 0, not NaN"
        );
    }

    #[test]
    fn failed_reps_are_counted_and_excluded_from_the_pool() {
        let mut bad = rep("rust", 1, 5, 5, 50);
        bad.error = Some("boom".into());
        let rows = vec![rep("rust", 0, 5, 5, 50), bad];
        let s = summarise_arm(&rows, 0, "rust");
        assert_eq!(s.reps, 1);
        assert_eq!(s.failed_reps, 1);
        assert!((s.fdr.estimate - 0.5).abs() < 1e-12);
    }

    #[test]
    fn se_calibration_needs_two_reps_and_lands_near_one_for_a_correct_sandwich() {
        let one = vec![rep("rust", 0, 5, 5, 50)];
        assert!(
            summarise_arm(&one, 0, "rust").se_calibration.is_nan(),
            "one rep has no per-taxon sd"
        );
        // Two reps with a known beta spread and a known se: the ratio must come
        // out finite and positive, which is the property the metric needs to be
        // reportable. A calibrated sandwich lands near 1 on real data, and that
        // is what the executed grid reports.
        let two = vec![rep("rust", 0, 5, 5, 50), rep("rust", 1, 5, 5, 50)];
        let cal = summarise_arm(&two, 0, "rust").se_calibration;
        assert!(cal.is_finite() && cal > 0.0, "got {cal}");
    }

    /// A two-cell grid, for the tests that need more than one.
    fn two_cell_grid() -> Grid {
        Grid {
            n_taxa: vec![100, 100],
            log_fc: vec![1.0, 2.0],
            ..grid()
        }
    }

    /// A regression test for a bug that produced a table of identical numbers:
    /// every cell's summary was pooled from *all* of the arm's rows.
    #[test]
    fn a_cell_summarises_only_its_own_replicates() {
        let mut rows = vec![rep("rust", 0, 50, 0, 50), rep("rust", 1, 50, 0, 50)];
        // A second cell with the opposite behaviour.
        for mut r in [rep("rust", 0, 1, 0, 50), rep("rust", 1, 1, 0, 50)] {
            r.cell = 1;
            rows.push(r);
        }
        let g = summarise(&two_cell_grid(), &rows);
        assert!((g.cells[0].rust.fdr.estimate - 0.0).abs() < 1e-12);
        assert!(
            (g.cells[0].rust.fdr.numerator - 0.0).abs() < 1e-12,
            "cell 0 made no false calls"
        );
        assert!(
            (g.cells[1].rust.power.estimate - 0.02).abs() < 1e-12,
            "cell 1 called 1 of 50 in each of two reps, so pooled power is 0.02, got {}",
            g.cells[1].rust.power.estimate
        );
        assert_eq!(g.cells[0].rust.reps, 2);
        assert_eq!(g.cells[1].rust.reps, 2);
    }

    #[test]
    fn identical_arms_agree_and_wildly_different_arms_diverge() {
        let rows = vec![
            rep("rust", 0, 5, 1, 50),
            rep("rust", 1, 5, 1, 50),
            rep("r", 0, 5, 1, 50),
            rep("r", 1, 5, 1, 50),
        ];
        let g = summarise(&grid(), &rows);
        assert!(g.compared_against_r);
        assert!(g.cells[0].agrees, "{:?}", g.cells[0].comparisons);

        // Rust calls 5 of 50 DA taxa with a 50% false rate; R calls all 50 with
        // none false. Both the FDR and the power gaps are outside the combined
        // 3 SE, so both are reported as findings.
        let rows = vec![
            rep("rust", 0, 5, 5, 50),
            rep("rust", 1, 5, 5, 50),
            rep("r", 0, 50, 0, 50),
            rep("r", 1, 50, 0, 50),
        ];
        let g = summarise(&grid(), &rows);
        assert!(!g.cells[0].agrees);
        assert_eq!(g.n_divergent, 1);
        assert!(!g.findings.is_empty());
        let joined = g.findings.join("\n");
        assert!(joined.contains("empirical_fdr"), "{joined}");
        assert!(joined.contains("power"), "{joined}");
    }

    /// Two absent rates are not a disagreement.
    ///
    /// A cell where neither arm called anything has no empirical FDR on either
    /// side. Scoring that as `divergent` -- which is what comparing `NaN` to
    /// `NaN` naturally does -- turns a sparse grid into a page of findings and
    /// buries the real ones.
    #[test]
    fn a_metric_with_no_value_on_either_side_is_inconclusive_not_divergent() {
        let rows = vec![rep("rust", 0, 0, 0, 50), rep("r", 0, 0, 0, 50)];
        let g = summarise(&grid(), &rows);
        let fdr = g.cells[0]
            .comparisons
            .iter()
            .find(|c| c.metric == "empirical_fdr")
            .expect("fdr is compared");
        assert_eq!(fdr.verdict, Verdict::Inconclusive);
        assert!(!fdr.rust.is_finite() && !fdr.r.is_finite());
        assert_eq!(g.n_divergent, 0, "an undefined rate is not a divergence");
        assert_eq!(g.n_inconclusive, 1);
        assert!(!g.cells[0].agrees, "an unresolved cell has not agreed");
        assert!(g.findings.is_empty());
    }

    /// The converse: one arm silent and the other not is still a divergence,
    /// because that asymmetry is exactly the bug a sparse cell would hide.
    #[test]
    fn a_rate_missing_on_only_one_side_is_still_a_divergence() {
        let rows = vec![rep("rust", 0, 0, 0, 50), rep("r", 0, 10, 1, 50)];
        let g = summarise(&grid(), &rows);
        let fdr = g.cells[0]
            .comparisons
            .iter()
            .find(|c| c.metric == "empirical_fdr")
            .expect("fdr is compared");
        assert_eq!(fdr.verdict, Verdict::Divergent);
        assert_eq!(g.n_divergent, 1);
        assert!(!g.findings.is_empty());
    }

    #[test]
    fn a_rust_only_run_never_claims_agreement_with_r() {
        let rows = vec![rep("rust", 0, 5, 1, 50)];
        let g = summarise(&grid(), &rows);
        assert!(!g.compared_against_r);
        assert!(g.cells[0].comparisons.is_empty());
        assert!(single_arm_agrees(&g.cells[0].comparisons));
    }

    #[test]
    fn the_verdict_follows_the_three_se_rule() {
        let a = Pooled::ratio(10.0, 100.0, 1, 0); // 0.1, se ~ 0.03
        let b = Pooled::ratio(11.0, 100.0, 1, 0); // 0.11
        assert_eq!(compare("fdr", &a, &b).verdict, Verdict::Agree);
        let c = Pooled::ratio(50.0, 100.0, 1, 0); // 0.5
        let d = compare("fdr", &a, &c);
        assert_eq!(d.verdict, Verdict::Divergent);
        assert!(d.tolerance > 0.0);
    }

    #[test]
    fn coverage_reports_which_cells_an_arm_filled() {
        let rows = vec![rep("rust", 0, 1, 0, 1)];
        assert_eq!(coverage(&grid(), &rows, "rust"), vec![0]);
        assert!(coverage(&grid(), &rows, "r").is_empty());
        assert_eq!(count_rows(&rows, "rust"), (1, 0));
    }
}
