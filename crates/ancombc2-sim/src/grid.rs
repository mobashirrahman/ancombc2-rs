//! The simulation grid: a declarative spec, expanded into cells.
//!
//! The grid file is JSON so the same spec drives the Rust arm, the R arm, and
//! the nightly CI job. `validation/simulation/quick` is what CI runs on a
//! pull request and `validation/simulation/full` is the weekly surface from
//! `PLAN.md` §5.3. Nothing about a cell is decided in code: a factor added to the
//! spec appears in both arms automatically.

use serde::{Deserialize, Serialize};
use std::path::Path;

/// One factor combination from the grid.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct Cell {
    /// Index of this cell in the expanded grid, stable for a given spec.
    pub index: usize,
    pub n_taxa: usize,
    pub n_samp: usize,
    /// Fraction of taxa that are differentially abundant, `(0, 1]`.
    pub da_proportion: f64,
    /// Magnitude of the log fold change for a DA taxon. The sign is random per
    /// taxon, so both directions are represented.
    pub log_fc: f64,
    /// Extra probability of a structural zero, on top of the zeros the
    /// negative binomial produces from low abundance.
    pub zero_inflation: f64,
    /// Mean sequencing depth.
    pub lib_mean: f64,
    /// Coefficient of variation of the depth.
    pub lib_cv: f64,
    /// Correlate the sampling fraction with the group, which is the
    /// confounding ANCOM-BC2's sampling-fraction correction exists to absorb.
    pub confound: bool,
}

/// A grid spec, as read from `grid.json`.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Grid {
    pub name: String,
    /// Master seed. Every cell and rep derives its own stream from this.
    pub seed: u64,
    pub n_taxa: Vec<usize>,
    pub n_samp: Vec<usize>,
    pub da_proportion: Vec<f64>,
    pub log_fc: Vec<f64>,
    pub zero_inflation: Vec<f64>,
    pub lib_mean: Vec<f64>,
    pub lib_cv: Vec<f64>,
    pub confound: Vec<bool>,
    /// Negative binomial dispersion, the `size` parameter.
    #[serde(default = "default_dispersion")]
    pub dispersion: f64,
    /// Overdispersion of the taxa's baseline abundances on the log scale.
    #[serde(default = "default_abundance_sd")]
    pub abundance_sd: f64,
    /// Replications per cell.
    pub reps: usize,
    /// Significance level for the reported FDR and power.
    #[serde(default = "default_alpha")]
    pub alpha: f64,
    /// Pseudo-count for the analysis.
    #[serde(default = "default_pseudo")]
    pub pseudo: f64,
    /// Which sensitivity analysis to run: `none`, `conservative` (0/0.1/0.5/1.0)
    /// or `nonconservative` (50 refits over 0.01..0.50).
    #[serde(default = "default_sens")]
    pub sensitivity: String,
    /// Prevalence filter.
    #[serde(default = "default_prevalence")]
    pub prevalence: f64,
    /// Library-size filter, as a fraction of the minimum.
    #[serde(default = "default_lib_size")]
    pub lib_size: f64,
    /// Multiplicity correction for the reported `q`.
    #[serde(default = "default_adjust")]
    pub p_adjust: String,
    /// Detect structural zeros.
    ///
    /// **Both arms read this from here.** It was previously hardcoded `TRUE` in
    /// `scripts/sim_r.R` and left at its `AncombcConfig::default()` of `false` in
    /// the Rust arm, so the two arms were not analysing the same thing: at 90%
    /// zero inflation the oracle retained 11 of 500 taxa and Rust retained all
    /// 500. Every Rust-versus-R cell comparison produced up to that point was
    /// therefore between two different analyses. The value is here, and
    /// `sim_r.R` reads it from the grid rather than from a literal, so the two
    /// arms cannot drift apart again.
    #[serde(default = "default_true")]
    pub struc_zero: bool,
    /// Classify structural zeros by the asymptotic lower bound rather than by
    /// outright absence. Read from the grid by both arms, as `struc_zero` is.
    #[serde(default = "default_true")]
    pub neg_lb: bool,
    /// Quantile of the SE distribution used for `s0`.
    #[serde(default = "default_s0_perc")]
    pub s0_perc: f64,
    /// A human note carried into every result row, e.g. what the grid omits.
    #[serde(default)]
    pub note: String,
    /// Named blocks of the grid.
    ///
    /// `PLAN.md` §5.3 lists the *levels* of seven factors. Their full cross
    /// product is 96,768 cells, which is not a design anyone can run, and a
    /// one-factor-at-a-time sweep would leave the interactions that matter
    /// unmeasured. A block list is the middle path: each block is a full cross
    /// product over a few factors at a fixed background, so every interaction
    /// the block covers is measured and every level of every factor appears
    /// somewhere. When `blocks` is absent the grid is the single top-level
    /// cross product, which is what a small smoke grid wants.
    #[serde(default)]
    pub blocks: Vec<Block>,
}

/// One block: a named, fully crossed set of factors at a fixed background.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Block {
    pub name: String,
    /// Overrides for the block's background. A field that is `None` inherits
    /// the grid's value, so a block only states what it changes.
    #[serde(default)]
    pub n_taxa: Option<Vec<usize>>,
    #[serde(default)]
    pub n_samp: Option<Vec<usize>>,
    #[serde(default)]
    pub da_proportion: Option<Vec<f64>>,
    #[serde(default)]
    pub log_fc: Option<Vec<f64>>,
    #[serde(default)]
    pub zero_inflation: Option<Vec<f64>>,
    #[serde(default)]
    pub lib_mean: Option<Vec<f64>>,
    #[serde(default)]
    pub lib_cv: Option<Vec<f64>>,
    #[serde(default)]
    pub confound: Option<Vec<bool>>,
    /// A block may override the analysis-level factors, which is how the
    /// sensitivity and filter blocks are expressed.
    #[serde(default)]
    pub reps: Option<usize>,
    #[serde(default)]
    pub sensitivity: Option<String>,
    #[serde(default)]
    pub prevalence: Option<f64>,
    #[serde(default)]
    pub lib_size: Option<f64>,
    #[serde(default)]
    pub p_adjust: Option<String>,
    /// Per-block overrides of the structural-zero screen. See the top-level
    /// `struc_zero` for why both arms read it from the grid.
    pub struc_zero: Option<bool>,
    pub neg_lb: Option<bool>,
    pub s0_perc: Option<f64>,
}

fn default_dispersion() -> f64 {
    1.0
}
fn default_abundance_sd() -> f64 {
    0.5
}
fn default_alpha() -> f64 {
    0.05
}
fn default_true() -> bool {
    true
}

fn default_s0_perc() -> f64 {
    0.05
}

fn default_pseudo() -> f64 {
    0.5
}
fn default_sens() -> String {
    "none".to_string()
}
fn default_prevalence() -> f64 {
    0.0
}
fn default_lib_size() -> f64 {
    0.0
}
fn default_adjust() -> String {
    "BH".to_string()
}

impl Grid {
    /// Load a grid from a spec file, or from a directory containing
    /// `grid.json`.
    ///
    /// Both forms are accepted because a grid is naturally addressed by its
    /// directory -- the Makefile, the CI job and `make sim-full` all name
    /// `validation/simulation/quick` rather than the file inside it -- and
    /// requiring the caller to know the file name is a chance for the same grid
    /// to be loaded from two paths.
    pub fn load(path: &Path) -> Result<Self, String> {
        let file = if path.is_dir() {
            path.join("grid.json")
        } else {
            path.to_path_buf()
        };
        let text = std::fs::read_to_string(&file)
            .map_err(|e| format!("cannot read grid {}: {e}", file.display()))?;
        let grid: Grid = serde_json::from_str(&text)
            .map_err(|e| format!("cannot parse grid {}: {e}", file.display()))?;
        grid.validate()?;
        Ok(grid)
    }

    /// Fail loudly on a spec that cannot produce a meaningful result, rather
    /// than silently running a degenerate grid.
    pub fn validate(&self) -> Result<(), String> {
        if self.n_taxa.is_empty() || self.n_samp.is_empty() {
            return Err("n_taxa and n_samp must be non-empty".into());
        }
        if self.reps == 0 {
            return Err("reps must be at least 1".into());
        }
        if !(0.0..=1.0).contains(&self.alpha) || self.alpha == 0.0 {
            return Err(format!("alpha must be in (0, 1], got {}", self.alpha));
        }
        if self.dispersion <= 0.0 {
            return Err(format!(
                "dispersion must be positive, got {}",
                self.dispersion
            ));
        }
        if self.pseudo < 0.0 {
            return Err("pseudo must be non-negative".into());
        }
        for p in &self.da_proportion {
            if !(0.0..=1.0).contains(p) {
                return Err(format!("da_proportion out of range: {p}"));
            }
        }
        for z in &self.zero_inflation {
            if !(0.0..1.0).contains(z) {
                return Err(format!("zero_inflation out of range: {z}"));
            }
        }
        for m in &self.lib_mean {
            if *m <= 0.0 {
                return Err(format!("lib_mean must be positive, got {m}"));
            }
        }
        for c in &self.lib_cv {
            if *c < 0.0 {
                return Err(format!("lib_cv must be non-negative, got {c}"));
            }
        }
        for n in &self.n_taxa {
            if *n < 2 {
                return Err(format!("n_taxa must be at least 2, got {n}"));
            }
        }
        for n in &self.n_samp {
            // ANCOM-BC2 needs at least two samples per group for the sample
            // intercepts and the within-group contrasts to be identified.
            if *n < 4 {
                return Err(format!("n_samp must be at least 4, got {n}"));
            }
        }
        for (i, b) in self.blocks.iter().enumerate() {
            if let Some(s) = &b.sensitivity {
                if !matches!(s.as_str(), "none" | "conservative" | "nonconservative") {
                    return Err(format!("block {i} ({}) unknown sensitivity: {s}", b.name));
                }
            }
            for lv in [
                b.n_taxa.clone().unwrap_or_default(),
                b.n_samp.clone().unwrap_or_default(),
            ] {
                if lv.iter().any(|v| *v < 2) {
                    return Err(format!("block {i} ({}) has a size below 2", b.name));
                }
            }
        }
        self.check_sensitivities()?;
        Ok(())
    }

    fn check_sensitivities(&self) -> Result<(), String> {
        match self.sensitivity.as_str() {
            "none" | "conservative" | "nonconservative" => {}
            other => return Err(format!("unknown sensitivity mode: {other}")),
        }
        Ok(())
    }

    /// Every factor combination, in a fixed order so `Cell::index` is stable.
    ///
    /// The order is the block order, and within a block the nested loop order of
    /// [`Block::levels`]. `Cell::index` therefore identifies a cell in a given
    /// spec file and is stable across runs of that file, which is what the
    /// per-cell reporting relies on.
    pub fn cells(&self) -> Vec<Cell> {
        let mut out = Vec::new();
        if self.blocks.is_empty() {
            for (i, lv) in self.top_levels().into_iter().enumerate() {
                out.push(Cell { index: i, ..lv });
            }
            return out;
        }
        for b in &self.blocks {
            for lv in b.levels(self) {
                out.push(Cell {
                    index: out.len(),
                    ..lv
                });
            }
        }
        out
    }

    /// Replications for a given cell: the block's override, else the grid's.
    pub fn reps_for(&self, cell: &Cell) -> usize {
        if self.blocks.is_empty() {
            return self.reps;
        }
        for b in &self.blocks {
            let start = self.block_start(b);
            if cell.index >= start && cell.index < start + block_cells(b, self).len() {
                return b.reps.unwrap_or(self.reps);
            }
        }
        self.reps
    }

    /// The grid's own factor levels, as a single cross product.
    fn top_levels(&self) -> Vec<Cell> {
        let lv = Levels {
            n_taxa: self.n_taxa.clone(),
            n_samp: self.n_samp.clone(),
            da_proportion: self.da_proportion.clone(),
            log_fc: self.log_fc.clone(),
            zero_inflation: self.zero_inflation.clone(),
            lib_mean: self.lib_mean.clone(),
            lib_cv: self.lib_cv.clone(),
            confound: self.confound.clone(),
        };
        cross(&lv)
    }

    fn block_start(&self, target: &Block) -> usize {
        let mut at = 0usize;
        for b in &self.blocks {
            if std::ptr::eq(b, target) {
                return at;
            }
            at += block_cells(b, self).len();
        }
        at
    }

    pub fn total_reps(&self) -> usize {
        if self.blocks.is_empty() {
            return self.cells().len() * self.reps;
        }
        self.blocks
            .iter()
            .map(|b| b.reps.unwrap_or(self.reps) * block_cells(b, self).len())
            .sum()
    }

    /// The analysis factors a cell runs with, including any block override.
    ///
    /// The Rust arm and the R arm both read these, so a block that changes the
    /// sensitivity mode cannot have one arm silently left on the default.
    pub fn analysis_for(&self, cell: &Cell) -> Analysis {
        let mut a = Analysis {
            reps: self.reps,
            sensitivity: self.sensitivity.clone(),
            prevalence: self.prevalence,
            lib_size: self.lib_size,
            p_adjust: self.p_adjust.clone(),
            struc_zero: self.struc_zero,
            neg_lb: self.neg_lb,
            s0_perc: self.s0_perc,
        };
        if self.blocks.is_empty() {
            return a;
        }
        for b in &self.blocks {
            let start = self.block_start(b);
            if cell.index >= start && cell.index < start + block_cells(b, self).len() {
                if let Some(r) = b.reps {
                    a.reps = r;
                }
                if let Some(s) = &b.sensitivity {
                    a.sensitivity = s.clone();
                }
                if let Some(p) = b.prevalence {
                    a.prevalence = p;
                }
                if let Some(l) = b.lib_size {
                    a.lib_size = l;
                }
                if let Some(p) = &b.p_adjust {
                    a.p_adjust = p.clone();
                }
                if let Some(v) = b.struc_zero {
                    a.struc_zero = v;
                }
                if let Some(v) = b.neg_lb {
                    a.neg_lb = v;
                }
                if let Some(v) = b.s0_perc {
                    a.s0_perc = v;
                }
            }
        }
        a
    }
}

/// The analysis-level factors for one cell.
#[derive(Debug, Clone, PartialEq)]
pub struct Analysis {
    pub reps: usize,
    pub sensitivity: String,
    pub prevalence: f64,
    pub lib_size: f64,
    pub p_adjust: String,
    pub struc_zero: bool,
    pub neg_lb: bool,
    pub s0_perc: f64,
}

struct Levels {
    n_taxa: Vec<usize>,
    n_samp: Vec<usize>,
    da_proportion: Vec<f64>,
    log_fc: Vec<f64>,
    zero_inflation: Vec<f64>,
    lib_mean: Vec<f64>,
    lib_cv: Vec<f64>,
    confound: Vec<bool>,
}

impl Block {
    fn resolved(&self, g: &Grid) -> Levels {
        Levels {
            n_taxa: self.n_taxa.clone().unwrap_or_else(|| g.n_taxa.clone()),
            n_samp: self.n_samp.clone().unwrap_or_else(|| g.n_samp.clone()),
            da_proportion: self
                .da_proportion
                .clone()
                .unwrap_or_else(|| g.da_proportion.clone()),
            log_fc: self.log_fc.clone().unwrap_or_else(|| g.log_fc.clone()),
            zero_inflation: self
                .zero_inflation
                .clone()
                .unwrap_or_else(|| g.zero_inflation.clone()),
            lib_mean: self.lib_mean.clone().unwrap_or_else(|| g.lib_mean.clone()),
            lib_cv: self.lib_cv.clone().unwrap_or_else(|| g.lib_cv.clone()),
            confound: self.confound.clone().unwrap_or_else(|| g.confound.clone()),
        }
    }

    pub fn levels(&self, g: &Grid) -> Vec<Cell> {
        cross(&self.resolved(g))
    }
}

fn block_cells(b: &Block, g: &Grid) -> Vec<Cell> {
    b.levels(g)
}

fn cross(lv: &Levels) -> Vec<Cell> {
    let mut out = Vec::new();
    for &n_taxa in &lv.n_taxa {
        for &n_samp in &lv.n_samp {
            for &da_proportion in &lv.da_proportion {
                for &log_fc in &lv.log_fc {
                    for &zero_inflation in &lv.zero_inflation {
                        for &lib_mean in &lv.lib_mean {
                            for &lib_cv in &lv.lib_cv {
                                for &confound in &lv.confound {
                                    out.push(Cell {
                                        index: out.len(),
                                        n_taxa,
                                        n_samp,
                                        da_proportion,
                                        log_fc,
                                        zero_inflation,
                                        lib_mean,
                                        lib_cv,
                                        confound,
                                    });
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn spec() -> Grid {
        Grid {
            name: "t".into(),
            seed: 1,
            n_taxa: vec![10, 20],
            n_samp: vec![8],
            da_proportion: vec![0.0, 0.5],
            log_fc: vec![1.0],
            zero_inflation: vec![0.0],
            lib_mean: vec![1e4],
            lib_cv: vec![0.3],
            confound: vec![false, true],
            dispersion: 1.0,
            abundance_sd: 0.5,
            struc_zero: true,
            neg_lb: true,
            s0_perc: 0.05,
            reps: 3,
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
    fn the_cartesian_product_is_the_full_cross_product() {
        let g = spec();
        let cells = g.cells();
        // 2 taxa x 1 sample count x 2 DA proportions x 2 confound levels; the
        // single-element factors are omitted.
        let expected = 2 * 2 * 2;
        assert_eq!(cells.len(), expected);
        assert_eq!(g.total_reps(), cells.len() * 3);
        assert!(cells.windows(2).all(|w| w[0].index + 1 == w[1].index));
    }

    #[test]
    fn a_degenerate_spec_is_rejected() {
        let mut g = spec();
        g.alpha = 0.0;
        assert!(g.validate().is_err());
        g = spec();
        g.dispersion = -1.0;
        assert!(g.validate().is_err());
        g = spec();
        g.zero_inflation = vec![1.5];
        assert!(g.validate().is_err());
        g = spec();
        g.sensitivity = "wild".into();
        assert!(g.validate().is_err());
        g = spec();
        g.n_samp = vec![2];
        assert!(
            g.validate().is_err(),
            "two samples cannot be split into two groups"
        );
    }

    #[test]
    fn a_directory_resolves_to_the_grid_json_inside_it() {
        let dir = std::env::temp_dir().join(format!("gridload-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let text = serde_json::to_string(&spec()).unwrap();
        std::fs::write(dir.join("grid.json"), &text).unwrap();
        let from_dir = Grid::load(&dir).unwrap();
        let from_file = Grid::load(&dir.join("grid.json")).unwrap();
        assert_eq!(from_dir.cells(), from_file.cells());
        assert_eq!(from_dir.name, from_file.name);
        // A directory with no grid.json says so rather than reading nothing.
        let empty = std::env::temp_dir().join(format!("gridempty-{}", std::process::id()));
        std::fs::create_dir_all(&empty).unwrap();
        assert!(Grid::load(&empty).is_err());
        std::fs::remove_dir_all(&dir).ok();
        std::fs::remove_dir_all(&empty).ok();
    }

    #[test]
    fn cell_expansion_is_deterministic() {
        assert_eq!(spec().cells(), spec().cells());
    }
}
