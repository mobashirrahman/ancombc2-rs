//! Golden parity: compare the Rust pipeline against the pinned ANCOMBC oracle.
//!
//! The goldens in `validation/golden/fx*/` were produced by
//! `scripts/generate_goldens.R`, which sources the reference's own R files at
//! commit `dc4febd` and records **every intermediate quantity**, not just the
//! final p-values. That is what makes a divergence localisable: the harness
//! walks the quantities in pipeline order and reports the *first* one that
//! disagrees, so a failure says "the sandwich variance is wrong for taxon 37" and
//! not "the p-values differ somewhere".
//!
//! Tolerances are the Levels A-D of `docs/numerical_contract.md`. The rule is
//! that they are never loosened to make a test pass; if a tolerance has to move,
//! the reason is recorded in the same commit.
//!
//! # `dead_code`
//!
//! This module is compiled once per test binary, and each binary uses a
//! different part of it -- `parity.rs` and `edge_cases.rs` both include it. A
//! helper one of them never calls is therefore dead code *in that binary* and
//! nothing else. The alternative, a shared crate, would mean publishing the
//! golden format as an API, which is a much larger commitment than the test
//! suite warrants.

#![allow(dead_code)] // see the note above: one copy per test binary
#![allow(clippy::needless_range_loop)]

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

use ancombc2_core::config::{AdjustMethod, AncombcConfig, EmControl, IterControl, MdfdrControl};
use ancombc2_core::matrix::Matrix;
use ancombc2_core::pipeline::CoreOutput;
use ancombc2_core::preprocess::{CountMatrix, Representation};
use ancombc2_core::AncombcResult;

// ---------------------------------------------------------------------------
// canonical readers
// ---------------------------------------------------------------------------

/// One golden quantity.
#[derive(Debug, Clone)]
pub enum Golden {
    /// `nr x nc` row-major f64.
    Matrix {
        nr: usize,
        nc: usize,
        v: Vec<f64>,
    },
    /// A flat f64 vector.
    Vector(Vec<f64>),
    /// `p x p` row-major blocks, one per taxon, in taxon order.
    Vcov {
        p: usize,
        v: Vec<f64>,
    },
    Strings(Vec<String>),
}

impl Golden {}

/// Read a little-endian f64 blob.
pub fn read_f64(path: &Path) -> Vec<f64> {
    let bytes = std::fs::read(path).unwrap_or_else(|e| panic!("read {}: {e}", path.display()));
    assert_eq!(
        bytes.len() % 8,
        0,
        "{} is not a whole number of doubles",
        path.display()
    );
    bytes
        .chunks_exact(8)
        .map(|c| f64::from_le_bytes(c.try_into().expect("8 bytes")))
        .collect()
}

/// A tiny JSON reader, sufficient for the manifest the R writer produces.
///
/// A full parser is unnecessary: the manifest is a flat object of scalars,
/// string arrays and small objects, and a dependency-free reader keeps the test
/// harness buildable with no network access.
#[derive(Debug, Clone, PartialEq, Default)]
pub enum Json {
    #[default]
    Null,
    Bool(bool),
    Num(f64),
    Str(String),
    Arr(Vec<Json>),
    Obj(Vec<(String, Json)>),
}

impl Json {
    pub fn get(&self, key: &str) -> Option<&Json> {
        match self {
            Json::Obj(kv) => kv.iter().find(|(k, _)| k == key).map(|(_, v)| v),
            _ => None,
        }
    }

    pub fn as_num(&self) -> Option<f64> {
        match self {
            Json::Num(x) => Some(*x),
            _ => None,
        }
    }

    pub fn as_str(&self) -> Option<&str> {
        match self {
            Json::Str(s) => Some(s),
            _ => None,
        }
    }

    pub fn as_arr(&self) -> Option<&[Json]> {
        match self {
            Json::Arr(a) => Some(a),
            _ => None,
        }
    }

    /// The entries of an object, in the order they were written.
    ///
    /// Order is preserved rather than sorted because the caller usually wants
    /// the file's own order -- a `zero_ind` table's columns are named, and
    /// reordering them would make a golden depend on a `BTreeMap` somewhere.
    pub fn as_obj(&self) -> Option<&[(String, Json)]> {
        match self {
            Json::Obj(o) => Some(o),
            _ => None,
        }
    }
}

pub fn parse_json(src: &str) -> Json {
    let b: Vec<char> = src.chars().collect();
    let mut i = 0usize;

    parse_value(&b, &mut i)
}

fn skip_ws(b: &[char], i: &mut usize) {
    while *i < b.len() && b[*i].is_whitespace() {
        *i += 1;
    }
}

fn parse_value(b: &[char], i: &mut usize) -> Json {
    skip_ws(b, i);
    match b.get(*i) {
        None => Json::Null,
        Some('{') => {
            *i += 1;
            let mut kv = Vec::new();
            loop {
                skip_ws(b, i);
                if b.get(*i) == Some(&'}') {
                    *i += 1;
                    break;
                }
                let k = match parse_value(b, i) {
                    Json::Str(s) => s,
                    other => panic!("object key must be a string, got {other:?}"),
                };
                skip_ws(b, i);
                assert_eq!(b.get(*i), Some(&':'), "expected ':'");
                *i += 1;
                let v = parse_value(b, i);
                kv.push((k, v));
                skip_ws(b, i);
                match b.get(*i) {
                    Some(',') => *i += 1,
                    Some('}') => {
                        *i += 1;
                        break;
                    }
                    other => panic!("expected ',' or '}}', got {other:?}"),
                }
            }
            Json::Obj(kv)
        }
        Some('[') => {
            *i += 1;
            let mut a = Vec::new();
            loop {
                skip_ws(b, i);
                if b.get(*i) == Some(&']') {
                    *i += 1;
                    break;
                }
                a.push(parse_value(b, i));
                skip_ws(b, i);
                match b.get(*i) {
                    Some(',') => *i += 1,
                    Some(']') => {
                        *i += 1;
                        break;
                    }
                    other => panic!("expected ',' or ']', got {other:?}"),
                }
            }
            Json::Arr(a)
        }
        Some('"') => {
            *i += 1;
            let mut s = String::new();
            while let Some(&c) = b.get(*i) {
                *i += 1;
                match c {
                    '"' => return Json::Str(s),
                    '\\' => {
                        let e = b[*i];
                        *i += 1;
                        match e {
                            'n' => s.push('\n'),
                            'r' => s.push('\r'),
                            't' => s.push('\t'),
                            'u' => {
                                let hex: String = b[*i..*i + 4].iter().collect();
                                *i += 4;
                                let cp = u32::from_str_radix(&hex, 16).expect("hex");
                                s.push(char::from_u32(cp).unwrap_or('?'));
                            }
                            other => s.push(other),
                        }
                    }
                    other => s.push(other),
                }
            }
            panic!("unterminated string")
        }
        Some('t') => {
            *i += 4;
            Json::Bool(true)
        }
        Some('f') => {
            *i += 5;
            Json::Bool(false)
        }
        Some('n') => {
            *i += 4;
            Json::Null
        }
        Some(_) => {
            let start = *i;
            while *i < b.len()
                && (b[*i].is_ascii_digit() || matches!(b[*i], '-' | '+' | '.' | 'e' | 'E'))
            {
                *i += 1;
            }
            let s: String = b[start..*i].iter().collect();
            Json::Num(s.parse().unwrap_or_else(|_| panic!("bad number {s:?}")))
        }
    }
}

/// One fixture's expected outputs, loaded from `validation/golden/fx*/`.
#[derive(Debug, Clone, Default)]
pub struct GoldenSet {
    pub quantities: BTreeMap<String, Golden>,
    /// Tables (JSON), keyed by name.
    pub tables: BTreeMap<String, Json>,
    pub config: Json,
    /// `validation/fixtures`, where `config.json` lives.
    pub config_dir: Option<PathBuf>,
    /// The directory this set was loaded from.
    ///
    /// Kept so the *sidecar* JSON -- `em_mixture.json`,
    /// `convergence_trace.json`, `stage_seconds.json` -- is reachable from here.
    /// Those are not quantities in the `.f64` contract, so they are not loaded
    /// into `quantities` or `tables`, but a comparison that needs one of them (see
    /// [`em_fits`]) has to be able to read the file the golden set was
    /// built from.
    pub dir: PathBuf,
}

/// `fx03` -> 3, so the sibling `validation/fixtures/fx03/config.json` is found
/// from a `validation/golden/fx03` directory.
fn fx_id(dir: &Path) -> usize {
    dir.file_name()
        .and_then(|s| s.to_str())
        .and_then(|s| s.strip_prefix("fx"))
        .and_then(|s| s.parse().ok())
        .unwrap_or(0)
}

impl GoldenSet {
    pub fn load(dir: &Path) -> GoldenSet {
        let manifest_path = dir.join("manifest.json");
        let src = std::fs::read_to_string(&manifest_path)
            .unwrap_or_else(|e| panic!("read {}: {e}", manifest_path.display()));
        let manifest = parse_json(&src);
        let Json::Obj(entries) = &manifest else {
            panic!("manifest must be an object");
        };
        let mut out = GoldenSet {
            dir: dir.to_path_buf(),
            ..GoldenSet::default()
        };
        for (name, entry) in entries {
            let kind = entry.get("kind").and_then(Json::as_str).unwrap_or("");
            let file = entry.get("file").and_then(Json::as_str);
            match kind {
                "f64" => {
                    let path = dir.join(file.expect("f64 entry needs a file"));
                    let v = read_f64(&path);
                    let nr = entry.get("nr").and_then(Json::as_num).unwrap_or(0.0) as usize;
                    // a vector entry has no "nc"; the R writer omits it, and
                    // `as_num` then returns the 0.0 default, so -1.0 is used as
                    // the sentinel for "absent"
                    let nc = match entry.get("nc").and_then(Json::as_num) {
                        Some(v) => v as i64,
                        None => -1,
                    };
                    if nc < 0 {
                        out.quantities.insert(name.clone(), Golden::Vector(v));
                    } else {
                        out.quantities.insert(
                            name.clone(),
                            Golden::Matrix {
                                nr,
                                nc: nc as usize,
                                v,
                            },
                        );
                    }
                }
                "vcov" => {
                    let path = dir.join(file.expect("vcov entry needs a file"));
                    let v = read_f64(&path);
                    // `n_tax` is recorded in the manifest but not read: the
                    // block count is `v.len() / p^2`, and asserting that against
                    // the manifest's value would catch a truncated blob, so it is
                    // checked here rather than stored.
                    let n_tax = entry.get("n_tax").and_then(Json::as_num).unwrap_or(0.0) as usize;
                    let p = entry.get("p").and_then(Json::as_num).unwrap_or(0.0) as usize;
                    assert_eq!(
                        v.len(),
                        n_tax * p * p,
                        "golden quantity {name}: the vcov blob is {} elements, but the \
                         manifest says {n_tax} taxa x {p} coefficients = {}",
                        v.len(),
                        n_tax * p * p
                    );
                    out.quantities.insert(name.clone(), Golden::Vcov { p, v });
                }
                "text" => {
                    let value = entry
                        .get("value")
                        .and_then(Json::as_arr)
                        .unwrap_or(&[])
                        .iter()
                        .filter_map(|v| v.as_str().map(str::to_string))
                        .collect();
                    out.quantities.insert(name.clone(), Golden::Strings(value));
                }
                "logical_matrix" | "table" => {
                    let path = dir.join(file.expect("table entry needs a file"));
                    let src = std::fs::read_to_string(&path)
                        .unwrap_or_else(|e| panic!("read {}: {e}", path.display()));
                    out.tables.insert(name.clone(), parse_json(&src));
                }
                other => panic!("unknown golden kind {other:?} for {name}"),
            }
        }
        // config.json lives with the fixture, not with the goldens
        let cfg_path = dir
            .parent()
            .and_then(|p| p.parent())
            .map(|fixtures| fixtures.join(format!("fx{:02}", fx_id(dir))))
            .map(|fixtures| fixtures.join("config.json"))
            .unwrap_or_else(|| dir.join("config.json"));
        if let Ok(cfg) = std::fs::read_to_string(&cfg_path) {
            out.config = parse_json(&cfg);
        }
        out.config_dir = dir
            .parent()
            .and_then(|p| p.parent())
            .map(|p| p.to_path_buf());
        out
    }

    pub fn matrix(&self, name: &str) -> (&[f64], usize, usize) {
        match self.quantities.get(name) {
            Some(Golden::Matrix { nr, nc, v }) => (v, *nr, *nc),
            other => panic!("golden {name} is not a matrix: {other:?}"),
        }
    }

    pub fn vector(&self, name: &str) -> &[f64] {
        match self.quantities.get(name) {
            Some(Golden::Vector(v)) | Some(Golden::Vcov { v, .. }) => v,
            other => panic!("golden {name} is not a vector: {other:?}"),
        }
    }

    pub fn strings(&self, name: &str) -> &[String] {
        match self.quantities.get(name) {
            Some(Golden::Strings(v)) => v,
            other => panic!("golden {name} is not a string vector: {other:?}"),
        }
    }
}

// ---------------------------------------------------------------------------
// fixture input
// ---------------------------------------------------------------------------

/// A fixture's inputs: the count table, the design and the config.
#[derive(Debug, Clone)]
pub struct Fixture {
    pub counts: CountMatrix,
    pub design: Matrix,
    pub group: Vec<usize>,
    /// The group factor's column name in the metadata, needed to locate its
    /// treatment-contrast columns in the design.
    pub group_name: String,
    pub taxon_names: Vec<String>,
    pub sample_names: Vec<String>,
    pub cfg: AncombcConfig,
}

// --- the fixture matrix of PLAN.md section 5.5 -------------------------------
//
// The four committed fixtures are addressed by a small integer. The matrix is
// addressed by cell *name*, because the names are the only thing that stays
// stable when the matrix grows, and a name reads in a failure report where an
// index does not.

/// A cell entry from the generated `validation/matrix/cells.json`.
///
/// Written out rather than exposing the raw [`Json`], so a typo in the manifest
/// fails here with the field named instead of surfacing as a `0.0` that quietly
/// means "absent" somewhere downstream.
#[derive(Debug, Clone)]
pub struct MatrixCell {
    name: String,
    inputs: String,
    fields: Vec<(String, Json)>,
}

impl MatrixCell {
    pub fn string(&self, key: &str) -> String {
        self.get(key)
            .and_then(Json::as_str)
            .unwrap_or_else(|| panic!("matrix cell {}: `{key}` is not a string", self.name))
            .to_string()
    }
    pub fn number(&self, key: &str) -> f64 {
        self.get(key)
            .and_then(Json::as_num)
            .unwrap_or_else(|| panic!("matrix cell {}: `{key}` is not a number", self.name))
    }
    pub fn boolean(&self, key: &str) -> bool {
        match self.get(key) {
            Some(Json::Bool(b)) => *b,
            // R's `write_json` emits `true`/`false`, but a numeric axis that
            // happens to be 0 or 1 is a plausible hand-edit and should not be
            // read as `false` silently.
            Some(Json::Num(n)) if *n == 0.0 || *n == 1.0 => *n == 1.0,
            _ => panic!("matrix cell {}: `{key}` is not a boolean", self.name),
        }
    }
    fn get(&self, key: &str) -> Option<&Json> {
        self.fields.iter().find(|(k, _)| k == key).map(|(_, v)| v)
    }
}

/// The declared matrix, read from the generator's own manifest.
///
/// Read rather than hard-coded on purpose: the coverage test in
/// `tests/fixture_matrix.rs` is only meaningful if it inspects what was actually
/// generated. Reading the R source instead would let the two drift, and the test
/// would keep passing against a matrix nobody built.
pub fn matrix_manifest() -> Vec<MatrixCell> {
    let path = repo_root().join("validation/matrix/cells.json");
    let src = std::fs::read_to_string(&path).unwrap_or_else(|e| {
        panic!(
            "read {}: {e}\nrun `Rscript scripts/generate_matrix_goldens.R`",
            path.display()
        )
    });
    let json = parse_json(&src);
    let cells = json.get("cells").and_then(Json::as_arr).unwrap_or_else(|| {
        panic!(
            "{} has no `cells` array; the matrix generator writes it, so this \
                 file is not one it produced",
            path.display()
        )
    });
    cells
        .iter()
        .map(|c| {
            let name = c
                .get("name")
                .and_then(Json::as_str)
                .expect("every cell needs a name")
                .to_string();
            let inputs = c
                .get("inputs")
                .and_then(Json::as_str)
                .expect("every cell needs an `inputs` directory")
                .to_string();
            let fields = match c {
                Json::Obj(kv) => kv.clone(),
                _ => panic!("a cell must be an object"),
            };
            MatrixCell {
                name,
                inputs,
                fields,
            }
        })
        .collect()
}

/// A golden set from a path relative to the repository root.
pub fn read_golden_at(rel: &str) -> GoldenSet {
    GoldenSet::load(&repo_root().join(rel))
}

/// Load the inputs of a matrix cell, and the configuration for that cell.
///
/// The inputs live under `validation/matrix/fixtures/<inputs>/` and are shared
/// by every cell that differs only in configuration -- which is the whole point
/// of the config-only sweeps. The config is therefore per *cell*, named
/// `config-<cell>.json`, and this takes the cell name for it.
pub fn load_matrix_fixture(inputs: &str, cell: &str) -> Fixture {
    let dir = repo_root().join("validation/matrix/fixtures").join(inputs);
    let (counts, design, group, taxon_names, sample_names) = load_inputs(&dir);
    let cfg = read_config(
        &dir.join(format!("config-{cell}.json")),
        design.colnames.clone(),
    );
    Fixture {
        counts,
        design,
        group,
        group_name: GROUP_VAR.to_string(),
        taxon_names,
        sample_names,
        cfg,
    }
}

/// The data half of a fixture, shared by the four committed fixtures and the
/// matrix cells: a counts matrix, the design `model.matrix` would build, and the
/// 0-based group indices.
fn load_inputs(dir: &Path) -> (CountMatrix, Matrix, Vec<usize>, Vec<String>, Vec<String>) {
    let (cols, taxon_names, rows) = read_tsv(&dir.join("counts.tsv"));
    let sample_names = cols;
    let mut data = vec![0.0; taxon_names.len() * sample_names.len()];
    for (i, r) in rows.iter().enumerate() {
        data[i * sample_names.len()..(i + 1) * sample_names.len()].copy_from_slice(r);
    }
    let counts = CountMatrix::new(taxon_names.len(), sample_names.len(), data).expect("counts");
    let meta = read_meta(&dir.join("meta.tsv"));
    let formula = std::fs::read_to_string(dir.join("formula.txt"))
        .expect("formula.txt")
        .trim()
        .to_string();
    // `data_sanity_check` coerces the group variable to a factor before the
    // model matrix is built, so the formula sees a factor even though the file
    // holds numbers.
    let design = build_design(&formula, &meta, &[GROUP_VAR]);
    let gcol = meta
        .numeric(GROUP_VAR)
        .expect("the fixtures always carry a `group` column");
    // the fixtures label groups 1..k; the core uses 0-based indices
    let group: Vec<usize> = gcol.iter().map(|&v| v as usize - 1).collect();
    (counts, design, group, taxon_names, sample_names)
}

/// Read a tab-separated table with row names.
///
/// `write.table(x, row.names = TRUE, quote = FALSE)` writes the column names in
/// the header with **no** leading blank field (for both matrices and data
/// frames), and the row name in a leading field of every data row. So the header
/// width is the column count and the data rows have one more field. Inferring
/// the layout from the first data row rather than assuming it avoids silently
/// dropping a column.
fn read_tsv(path: &Path) -> (Vec<String>, Vec<String>, Vec<Vec<f64>>) {
    let src =
        std::fs::read_to_string(path).unwrap_or_else(|e| panic!("read {}: {e}", path.display()));
    let mut lines = src.lines().filter(|l| !l.trim().is_empty());
    let cols: Vec<String> = lines
        .next()
        .expect("a header row")
        .split('\t')
        .map(str::trim)
        .map(str::to_string)
        .collect();
    let ncol = cols.len();

    let mut rows: Vec<Vec<f64>> = Vec::new();
    let mut names: Vec<String> = Vec::new();
    for line in lines {
        let fields: Vec<&str> = line.split('\t').collect();
        assert_eq!(
            fields.len(),
            ncol + 1,
            "ragged row in {}: {} fields, expected {}",
            path.display(),
            fields.len(),
            ncol + 1
        );
        names.push(fields[0].trim().to_string());
        rows.push(
            fields[1..]
                .iter()
                .map(|s| {
                    let t = s.trim();
                    if t.is_empty() || t == "NA" {
                        f64::NAN
                    } else {
                        t.parse().unwrap_or_else(|_| panic!("bad number {t:?}"))
                    }
                })
                .collect(),
        );
    }
    (cols, names, rows)
}

pub fn repo_root() -> PathBuf {
    if let Ok(r) = std::env::var("ANCOMBC_REPO") {
        return PathBuf::from(r);
    }
    // CARGO_MANIFEST_DIR is <repo>/crates/ancombc2-core, so three levels up is
    // the workspace root.
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .ancestors()
        .nth(2)
        .map(PathBuf::from)
        .expect("repo root")
}

/// The group factor's name in every fixture and matrix cell's metadata.
///
/// Stated once because three separate places need it -- the design builder, the
/// 0-based group index, and the rank-deficiency mask's contrast-column lookup --
/// and a fixture whose metadata called it something else would otherwise fail in
/// the mask only.
pub const GROUP_VAR: &str = "group";

pub fn fixture_dir(id: usize) -> PathBuf {
    repo_root()
        .join("validation/fixtures")
        .join(format!("fx{id:02}"))
}

pub fn read_golden(id: usize) -> GoldenSet {
    GoldenSet::load(&golden_dir(id))
}

/// The golden directory of one fixture-matrix cell, e.g. `predictor-5group`.
pub fn cell_golden_dir(cell: &str) -> PathBuf {
    repo_root().join("validation/matrix/golden").join(cell)
}

pub fn golden_dir(id: usize) -> PathBuf {
    repo_root()
        .join("validation/golden")
        .join(format!("fx{id:02}"))
}

/// Build a design matrix from the formula and the metadata, the way R's
/// `model.matrix` would, so the Rust side never has to parse a formula.
///
/// The grammar supported is the subset the fixtures and the CLI need:
/// `a + b + a:b` with `a`/`b` either numeric (a single column named after the
/// variable) or a factor (treatment contrasts against the first level, named
/// `var<level>`, e.g. `group2` for the second level of `group`). That is exactly
/// what R produces under the default `contr.treatment`, and getting the *names*
/// right matters because the multi-group tests locate the group columns by
/// substring.
///
/// `factor_vars` names the variables to treat as factors. R's `data_sanity_check`
/// coerces the `group` variable with `as.factor()` **before** the model matrix is
/// built, so a variable that is numeric on disk can still enter the formula as a
/// factor. Reproducing that coercion is why this is an explicit list rather than
/// a guess from the number of distinct values.
/// Expand a formula the way R's `terms()` does, into `+`-separated terms.
///
/// Only the operators the fixtures and the CLI use are handled: `+`, `:`, `*`
/// and `/`. `a * b` becomes `a + b + a:b` and `a / b` becomes `a + a:b`, which is
/// exactly R's expansion and therefore produces the same columns in the same
/// order -- and the order matters, because `fix_eff` names are part of the
/// Level A contract. `(`/`)` are stripped rather than nested: no fixture uses a
/// parenthesised formula, and pretending to parse one would be worse than not
/// parsing it.
fn expand_formula(formula: &str) -> Vec<String> {
    let f = formula.replace(['(', ')'], " ");
    let mut out: Vec<String> = Vec::new();
    for raw in f.split('+') {
        let term = raw.trim();
        if term.is_empty() {
            continue;
        }
        // `*` is left-associative in R and expands pairwise, so `a * b * c`
        // becomes `a + b + a:b + c + a:c + b:c + a:b:c`.
        if term.contains('*') {
            let factors: Vec<&str> = term
                .split('*')
                .map(str::trim)
                .filter(|s| !s.is_empty())
                .collect();
            let mut acc: Vec<String> = factors.iter().map(|s| s.to_string()).collect();
            for i in 0..factors.len() {
                for j in (i + 1)..factors.len() {
                    // `terms()` writes an interaction label with its variables
                    // in sorted order, so `x10 * x1` is `x1:x10`. The column
                    // name is part of the Level A contract, so writing it the
                    // other way round would be a parity failure.
                    let mut pair = [factors[i], factors[j]];
                    pair.sort_unstable();
                    acc.push(format!("{}:{}", pair[0], pair[1]));
                }
            }
            out.extend(acc);
            continue;
        }
        if term.contains('/') {
            let factors: Vec<&str> = term
                .split('/')
                .map(str::trim)
                .filter(|s| !s.is_empty())
                .collect();
            let mut acc: Vec<String> = vec![factors[0].to_string()];
            for f in &factors[1..] {
                let mut pair = [factors[0], *f];
                pair.sort_unstable();
                acc.push(format!("{}:{}", pair[0], pair[1]));
            }
            out.extend(acc);
            continue;
        }
        out.push(term.to_string());
    }
    // R's `terms()` keeps each term once, at its first occurrence, and
    // deduplicates across the whole formula. Without this, `group + x1 + x10 * x1`
    // yields `x1` and `x10` twice -- once from the `+` terms and once from the
    // `*` expansion -- and the design gains two exactly duplicated columns, which
    // makes it rank deficient and every taxon unfittable.
    let mut seen: Vec<String> = Vec::new();
    out.retain(|t| {
        if seen.iter().any(|s| s == t) {
            false
        } else {
            seen.push(t.clone());
            true
        }
    });
    out
}

/// A numeric column for `name`, or `None` when it is a factor.
///
/// A name listed in `factor_vars` is a factor even if its values parse as
/// numbers, because that is exactly the `group` case: the file holds `1` and `2`
/// and R is told to treat them as levels.
fn numeric_operand(meta: &MetaTable, name: &str, factor_vars: &[&str]) -> Option<Vec<f64>> {
    if factor_vars.contains(&name) {
        return None;
    }
    meta.numeric(name)
}

pub fn build_design(formula: &str, meta: &MetaTable, factor_vars: &[&str]) -> Matrix {
    let mut cols: Vec<Vec<f64>> = Vec::new();
    let mut names: Vec<String> = Vec::new();
    let push_intercept = true;
    for term in expand_formula(formula) {
        let term = term.as_str();
        if let Some((a, b)) = term.split_once(':') {
            let a = a.trim();
            let b = b.trim();
            // An operand can be numeric or a factor, and the two kinds give
            // different columns. A numeric operand contributes itself; a factor
            // contributes treatment-contrast dummies and the interaction is their
            // cross-product. Treating a numeric covariate as a factor -- which is
            // what this did -- enumerates its 10,000 distinct *values* as
            // "levels" and builds a 100-million-column design.
            let a_num = numeric_operand(meta, a, factor_vars);
            let b_num = numeric_operand(meta, b, factor_vars);
            match (a_num, b_num) {
                (Some(va), Some(vb)) => {
                    cols.push((0..meta.n).map(|r| va[r] * vb[r]).collect());
                    names.push(format!("{a}:{b}"));
                }
                (Some(va), None) => {
                    let fb = factor_levels(meta, b);
                    for (j, _) in fb.levels.iter().enumerate().skip(1) {
                        let col: Vec<f64> = (0..meta.n)
                            .map(|r| va[r] * if fb.code[r] == Some(j) { 1.0 } else { 0.0 })
                            .collect();
                        cols.push(col);
                        names.push(format!("{a}:{b}{}", suffix_for(j, &fb.levels[j])));
                    }
                }
                (None, Some(vb)) => {
                    let fa = factor_levels(meta, a);
                    for (i, _) in fa.levels.iter().enumerate().skip(1) {
                        let col: Vec<f64> = (0..meta.n)
                            .map(|r| vb[r] * if fa.code[r] == Some(i) { 1.0 } else { 0.0 })
                            .collect();
                        cols.push(col);
                        names.push(format!("{a}:{b}{}", suffix_for(i, &fa.levels[i])));
                    }
                }
                (None, None) => {
                    let fa = factor_levels(meta, a);
                    let fb = factor_levels(meta, b);
                    for (i, _) in fa.levels.iter().enumerate().skip(1) {
                        for (j, _) in fb.levels.iter().enumerate().skip(1) {
                            let col: Vec<f64> = (0..meta.n)
                                .map(|r| {
                                    if fa.code[r] == Some(i) && fb.code[r] == Some(j) {
                                        1.0
                                    } else {
                                        0.0
                                    }
                                })
                                .collect();
                            cols.push(col);
                            names.push(format!(
                                "{a}:{b}{}{}",
                                if j == 1 { "" } else { &fb.levels[j] },
                                suffix_for(i, &fa.levels[i])
                            ));
                        }
                    }
                }
            }
            continue;
        }
        if !factor_vars.contains(&term) {
            if let Some(v) = meta.numeric(term) {
                cols.push(v);
                names.push(term.to_string());
                continue;
            }
        }
        let f = factor_levels(meta, term);
        for i in 1..f.levels.len() {
            let col: Vec<f64> = (0..meta.n)
                .map(|r| if f.code[r] == Some(i) { 1.0 } else { 0.0 })
                .collect();
            cols.push(col);
            names.push(format!("{term}{}", f.levels[i]));
        }
    }
    if push_intercept {
        let mut all = vec![vec![1.0; meta.n]];
        let mut all_names = vec!["(Intercept)".to_string()];
        all.extend(cols);
        all_names.extend(names);
        return Matrix::from_cols(&all)
            .with_rownames(meta.row_names.clone())
            .with_colnames(all_names);
    }
    Matrix::from_cols(&cols)
        .with_rownames(meta.row_names.clone())
        .with_colnames(names)
}

fn suffix_for(i: usize, level: &str) -> String {
    format!("{level}{i}")
}

struct FactorInfo {
    levels: Vec<String>,
    code: Vec<Option<usize>>,
}

fn factor_levels(meta: &MetaTable, name: &str) -> FactorInfo {
    let raw = meta
        .cols
        .iter()
        .position(|c| c == name)
        .unwrap_or_else(|| panic!("variable {name} not in metadata"));
    let values = &meta.data[raw];
    let mut levels: Vec<String> = Vec::new();
    for v in values {
        let s = format!("{v}");
        if !levels.contains(&s) {
            levels.push(s);
        }
    }
    // R sorts character levels; the fixtures use 1..k so this matches
    if levels.iter().all(|l| l.parse::<f64>().is_ok()) {
        levels.sort_by(|a, b| {
            a.parse::<f64>()
                .unwrap()
                .partial_cmp(&b.parse::<f64>().unwrap())
                .unwrap()
        });
    } else {
        levels.sort();
    }
    let code = values
        .iter()
        .map(|v| levels.iter().position(|l| l == &format!("{v}")))
        .collect();
    FactorInfo { levels, code }
}

/// A metadata table: numeric columns only, plus the group column.
#[derive(Debug, Clone, Default)]
pub struct MetaTable {
    pub cols: Vec<String>,
    pub data: Vec<Vec<f64>>,
    pub n: usize,
    pub row_names: Vec<String>,
}

impl MetaTable {
    pub fn numeric(&self, name: &str) -> Option<Vec<f64>> {
        let i = self.cols.iter().position(|c| c == name)?;
        Some(self.data[i].clone())
    }
}

pub fn read_meta(path: &Path) -> MetaTable {
    let (cols, row_names, rows) = read_tsv(path);
    let n = rows.len();
    let ncol = cols.len();
    let mut data = vec![Vec::new(); ncol];
    for r in &rows {
        assert_eq!(r.len(), ncol, "ragged metadata in {}", path.display());
        for (j, v) in r.iter().enumerate() {
            data[j].push(*v);
        }
    }
    MetaTable {
        cols,
        data,
        n,
        row_names,
    }
}

pub fn load_fixture(id: usize) -> Fixture {
    let dir = fixture_dir(id);
    let (counts, design, group, taxon_names, sample_names) = load_inputs(&dir);
    // The design's column names *are* the fixed effects; R's `model.matrix`
    // names them and the multi-group tests locate columns by substring, so they
    // must reach the core rather than being regenerated from the config.
    let cfg = read_config(&dir.join("config.json"), design.colnames.clone());
    Fixture {
        counts,
        design,
        group,
        group_name: GROUP_VAR.to_string(),
        taxon_names,
        sample_names,
        cfg,
    }
}

fn read_config(path: &Path, colnames: Vec<String>) -> AncombcConfig {
    let src =
        std::fs::read_to_string(path).unwrap_or_else(|e| panic!("read {}: {e}", path.display()));
    let j = parse_json(&src);
    let b = |k: &str| j.get(k).and_then(Json::as_num).unwrap_or(0.0);
    let s = |k: &str| j.get(k).and_then(Json::as_str).unwrap_or("").to_string();
    AncombcConfig {
        // The parity harness compares `o1`, `o2`, `y1`, `y2` and `y_bias_crt`
        // against the goldens, so it is the one caller that needs the full-size
        // intermediates. They are off by default -- five times the count matrix,
        // and read by nothing else.
        keep_intermediates: true,
        // The golden harness exercises the default representation. The sparse one
        // is checked for equivalence by its own unit test in `preprocess.rs` and by
        // `docs/compatibility.md`'s measurement; it produces identical selected
        // tables, so running the full matrix through both would test the flag
        // rather than the algorithm.
        representation: Representation::Dense,
        fix_eff: colnames,
        p_adj_method: AdjustMethod::parse(&s("p_adj_method")).expect("adjust method"),
        pseudo: b("pseudo"),
        pseudo_sens: j
            .get("pseudo_sens")
            .and_then(|v| match v {
                Json::Bool(b) => Some(*b),
                _ => None,
            })
            .unwrap_or(false),
        conservative: j
            .get("conservative")
            .and_then(|v| match v {
                Json::Bool(b) => Some(*b),
                _ => None,
            })
            .unwrap_or(true),
        prv_cut: b("prv_cut"),
        lib_cut: b("lib_cut"),
        s0_perc: b("s0_perc"),
        group: if s("group").is_empty() {
            None
        } else {
            Some(s("group"))
        },
        group_labels: None,
        struc_zero: j
            .get("struc_zero")
            .and_then(|v| match v {
                Json::Bool(b) => Some(*b),
                _ => None,
            })
            .unwrap_or(false),
        neg_lb: j
            .get("neg_lb")
            .and_then(|v| match v {
                Json::Bool(b) => Some(*b),
                _ => None,
            })
            .unwrap_or(false),
        alpha: b("alpha"),
        global: j
            .get("global")
            .and_then(|v| match v {
                Json::Bool(b) => Some(*b),
                _ => None,
            })
            .unwrap_or(false),
        pairwise: j
            .get("pairwise")
            .and_then(|v| match v {
                Json::Bool(b) => Some(*b),
                _ => None,
            })
            .unwrap_or(false),
        iter_control: IterControl::default(),
        em_control: EmControl::default(),
        mdfdr_control: MdfdrControl::default(),
        compat: ancombc2_core::CompatMode::Ancombc2_15,
    }
}

// ---------------------------------------------------------------------------
// comparison
// ---------------------------------------------------------------------------

/// A level of the numerical contract.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Level {
    /// Exact: bitwise, or within a few ulp.
    A,
    /// Floating point, relative.
    B,
    /// Probabilities, absolute.
    C,
}

/// The first divergence found, with enough context to act on it.
#[derive(Debug, Clone)]
pub struct Divergence {
    pub quantity: String,
    pub level: Level,
    pub index: Option<usize>,
    pub taxon: Option<String>,
    pub coefficient: Option<String>,
    pub got: f64,
    pub want: f64,
    pub abs_diff: f64,
    pub rel_diff: f64,
    pub message: String,
}

impl std::fmt::Display for Divergence {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "  quantity : {}", self.quantity)?;
        writeln!(f, "  level    : {:?}", self.level)?;
        if let Some(i) = self.index {
            write!(f, "  index    : {i}")?;
            if let Some(t) = &self.taxon {
                write!(f, "  taxon    : {t}")?;
            }
            if let Some(c) = &self.coefficient {
                write!(f, "  coef     : {c}")?;
            }
            writeln!(f)?;
        }
        writeln!(f, "  got      : {:.17e}", self.got)?;
        writeln!(f, "  expected : {:.17e}", self.want)?;
        writeln!(f, "  abs diff : {:.6e}", self.abs_diff)?;
        writeln!(f, "  rel diff : {:.6e}", self.rel_diff)?;
        if !self.message.is_empty() {
            writeln!(f, "  note     : {}", self.message)?;
        }
        Ok(())
    }
}

#[derive(Debug, Clone, Copy)]
pub struct Tolerance {
    pub level: Level,
    pub rtol: f64,
    pub atol: f64,
}

impl Tolerance {
    /// Widen the absolute term to account for cancellation in a difference.
    ///
    /// `beta - delta` is bounded by the *operands'* errors, not by the size of
    /// the result: when a coefficient of ~60 is corrected by a bias of ~60 to
    /// leave ~0.06, the operands' 1e-8 relative accuracy leaves ~1e-8 relative
    /// error in the answer. Comparing the result against its own value would
    /// demand more than the inputs can deliver, so the absolute term is scaled by
    /// the largest operand.
    pub fn with_operand_scale(mut self, operands: &[f64]) -> Self {
        let scale = operands.iter().fold(0.0f64, |m, v| m.max(v.abs()));
        self.atol = self.atol.max(self.rtol * scale);
        self
    }
}

pub const TOL_BETA: Tolerance = Tolerance {
    level: Level::B,
    rtol: 1e-8,
    atol: 1e-14,
};
pub const TOL_SE: Tolerance = Tolerance {
    level: Level::B,
    rtol: 1e-7,
    atol: 1e-16,
};

/// Level B for the E-M's outputs: `delta_em`, `delta_wls`, `var_delta`.
///
/// `rtol 1e-7`, as PLAN.md section 5 and `docs/numerical_contract.md` section 4
/// both state for `delta_em` and `delta_wls`. These three quantities were
/// previously compared with [`TOL_BETA`]'s `1e-8`, which is *stricter* than the
/// contract -- and it passed only because the four committed fixtures are small.
/// The `covariates-10` matrix cell (100 x 30 with ten covariates, so `p = 12`)
/// reaches 9.7e-8 relative on `delta_wls`, inside the specified 1e-7 and outside
/// the accidental 1e-8.
///
/// The reason the deviation grows with `p` is the quantity's nature rather than
/// the implementation's: `delta_wls` is the WLS solution of a three-component
/// mixture fitted to all `p` coefficients at once, so the fit is `p`-dimensional
/// and its summation order -- and therefore its last bit -- depends on `p`. The
/// E-M's own convergence tolerance of 1e-5 bounds how well the *algorithm*
/// determines the answer, and two implementations that both run it to that
/// stopping rule agree on it to near machine precision rather than to 1e-8.
pub const TOL_EM: Tolerance = Tolerance {
    level: Level::B,
    rtol: 1e-7,
    atol: 1e-16,
};
pub const TOL_VCOV: Tolerance = Tolerance {
    level: Level::B,
    rtol: 1e-7,
    atol: 1e-18,
};

/// Level B for `s0`, the SAM regularisation constant: `rtol 1e-9`.
///
/// PLAN.md section 5 specifies `rtol 1e-9` for `s0`, two orders of magnitude
/// tighter than the `1e-7` it sets for `se`/`delta_em`/`delta_wls`/`vcov`, and
/// this used [`TOL_SE`] instead -- so the one quantity the contract singles out as
/// needing the tightest tolerance was being checked two orders looser than
/// specified, and passing.
///
/// The tighter bound is not a formality. `s0` is a *quantile* of the standard
/// errors (`quantile(se, s0_perc)`, R type 7), so it is an interpolation between
/// two order statistics rather than a per-taxon estimate: there is no
/// averaging-out over taxa, and any error in any one `se` moves it. That makes it
/// strictly harder to satisfy than `se` itself, not easier, and it holds: measured
/// agreement over the committed fixtures is better than 4e-10 relative.
pub const TOL_S0: Tolerance = Tolerance {
    level: Level::B,
    rtol: 1e-9,
    atol: 0.0,
};

/// Read one fixed effect's recorded E-M stopping state: `(converged, epsilon)`.
///
/// `epsilon` is the size of the E-M's **last parameter step**, which is the
/// precision to which it determined its parameters -- and therefore the precision
/// everything downstream of it inherits. `iterations` against `max_iter` is kept
/// as well, because it is what `em_tolerance_for` falls back to when a golden
/// records no epsilon at all. `None` when the golden predates the capture.
/// Level C: p-values and adjusted p-values, at the plan's `atol = 1e-10`.
///
/// This is a *floor*, not the whole check. A p-value is `2 * (1 - F_t(|W|))`,
/// so `|dp/dW| = 2 f_t(|W|) <= 0.8`: an error in `W` propagates into `p` at no
/// more than unit gain, and the caller adds the *observed* `W` discrepancy to
/// this floor before comparing (`tol_p` below). That is error propagation rather
/// than a loosened bound, and it is why the floor can be the plan's 1e-10 even
/// though the largest `p` disagreement observed is 1.6e-10: the `W` it inherits
/// from is itself 2.1e-10 out on the same fixture.
///
/// Measured, not assumed. `ANCOMBC2_GOLDEN_AUDIT=1` records the largest
/// `|rust - oracle|` for every compared quantity; over the four committed
/// fixtures, on the three where parity is asserted (no rank-deficient taxon), the
/// worst case is:
///
/// | fixture | p | q |
/// | --- | --- | --- |
/// | `fx01` 10x10 | 2.6e-15 | 0 (exact) |
/// | `fx02` 100x30 | 7.6e-13 | 1.1e-13 |
/// | `fx03` 1000x100 | 1.6e-10 | 2.0e-10 |
///
/// `fx04` reaches 2.0e-4 for `p` and 7.0e-4 for `q`, but 200 of its 10,000 taxa
/// have an exactly singular sub-design, so those quantities are *reported, not
/// asserted* -- see [`INDIRECT_QUANTITIES`]. The root cause of `fx03`'s 1.6e-10
/// is `delta_em` at 1.7e-11, which is inside its own Level B contract; the
/// residual is summation order in the E-M sweep, not a defect.
///
/// An earlier version of this constant was `atol = 1e-8`, justified by the E-M
/// iteration tolerance of 1e-5. That reasoning was wrong: the E-M tolerance says
/// how accurately the *algorithm* converges, not how closely two implementations
/// of the same algorithm agree when both run it to the same stopping rule. The
/// measurement above is what the contract is actually worth, and it is two orders
/// of magnitude tighter than the value it replaces.
pub const TOL_PROB: Tolerance = Tolerance {
    level: Level::C,
    atol: 1e-10,
    rtol: 0.0,
};

/// The log transform and the per-taxon mean.
///
/// This is *not* "exact": R accumulates the row mean in extended precision
/// (`LDOUBLE` in `rowMeans`), so a log followed by a mean reduction differs from
/// an `f64` reduction in the last couple of ulp. The centring subtracts two
/// O(1) quantities, so a *centred* entry near zero can carry a large relative
/// error while its absolute error stays at the ulp level; the tolerance is
/// therefore stated on the absolute scale (1e-13) with a loose relative term.
/// That is far tighter
/// than anything downstream depends on — the sandwich variance inherits the
/// response and would amplify a 1e-15 discrepancy into a 1e-14 one — while
/// staying honest about the reduction order.
pub const TOL_PREPROCESS: Tolerance = Tolerance {
    level: Level::A,
    rtol: 1e-12,
    atol: 1e-13,
};

/// Compare two vectors under a tolerance, returning the first failure.
/// The largest absolute deviation seen per quantity, over everything compared in
/// this process.
///
/// A tolerance is a claim about what has actually been measured. The golden
/// fixtures only print a deviation when one *exceeds* the tolerance, so a run that
/// passes says nothing about how much headroom there was -- and a tolerance left
/// at 1e-8 because "that is roughly what E-M converges to" is a guess about the
/// algorithm, not a measurement of this implementation against this oracle. This
/// records the observed deviation for every entry so `make golden-audit` can
/// report what the contract actually achieves, and the tolerance can be set to
/// the plan's value where the measurement supports it.
static DEVIATIONS: std::sync::Mutex<Option<std::collections::BTreeMap<String, f64>>> =
    std::sync::Mutex::new(None);

/// The accumulator is diagnostic, so a panic elsewhere must not cascade through
/// it: recovering the guard keeps the remaining comparisons reporting their real
/// deviations instead of every one of them reporting a `PoisonError` and hiding
/// the original failure behind seven spurious ones.
fn deviations() -> std::sync::MutexGuard<'static, Option<std::collections::BTreeMap<String, f64>>> {
    DEVIATIONS.lock().unwrap_or_else(|e| e.into_inner())
}

fn record_deviation(quantity: &str, abs: f64) {
    if !abs.is_finite() {
        return;
    }
    let mut g = deviations();
    let m = g.get_or_insert_with(Default::default);
    let e = m.entry(quantity.to_string()).or_insert(0.0);
    if abs > *e {
        *e = abs;
    }
}

/// The audit has to actually record, or the tolerances documented against it are
/// fiction: a broken recorder would print nothing and the numbers quoted in
/// `TOL_PROB` and `docs/numerical_contract.md` would be unfalsifiable.
#[cfg(test)]
mod s0_bound_tests {
    use super::{em_fits, em_tolerance_for, TOL_EM, TOL_S0};

    /// The bound `s0` is held to must follow the E-M's own convergence, and that
    /// choice is the whole mechanism, so it is pinned here rather than only
    /// exercised end to end.
    ///
    /// Both directions matter. If a *converged* term were held to the looser bound,
    /// a genuine 1e-8 error in `s0` would pass; if a *non-converged* term were held
    /// to `TOL_S0`, the contract would be unsatisfiable on any design where the E-M
    /// stops at `max_iter`, which is most wide designs.
    #[test]
    fn the_s0_bound_follows_the_em_convergence() {
        // The rule `s0` and the mixture share: widen to the E-M's recorded final
        // step when that is coarser than the contract, and never tighten below it.
        assert_eq!(em_tolerance_for(None).rtol, TOL_EM.rtol);
        assert_eq!(
            em_tolerance_for(Some((true, 1e-12))).rtol,
            TOL_EM.rtol,
            "a converged E-M that stopped well inside its tolerance must not widen the bound"
        );
        // A fit that reached `max_iter` with a final step of 9.1e-6 determined its
        // parameters no better than that, and `s0` inherits them through
        // `var_delta` -- which is one scalar per coefficient, so nothing averages
        // out. This is the case that leaves `s0` 1.07e-9 out on
        // `covariates-10-interaction`, outside `rtol 1e-9`.
        assert_eq!(em_tolerance_for(Some((false, 9.1e-6))).rtol, 9.1e-6);
        assert_eq!(em_tolerance_for(Some((false, 1e-12))).rtol, TOL_EM.rtol);
        // The relationship the mechanism depends on, checked as a constant rather
        // than at runtime: a tolerance that stopped being tighter than the one it
        // derives from would silently turn the whole mechanism into a no-op.
        const _: () = assert!(TOL_S0.rtol < TOL_EM.rtol);
    }

    /// `em_fits` reads `iterations`, `max_iter` and `epsilon` out of the recorded
    /// mixture, so the bound follows the *recorded run* rather than a
    /// hand-maintained list of cell names.
    #[test]
    fn em_convergence_is_read_from_the_recorded_mixture() {
        let dir = std::env::temp_dir().join(format!("ancombc2-em-conv-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        std::fs::write(
            dir.join("em_mixture.json"),
            r#"{"iterations":[44,57,53,100,92],"max_iter":100,
                "epsilon":[1e-9,2e-9,1e-9,7e-4,1e-9]}"#,
        )
        .unwrap();
        let fits = em_fits(&dir).expect("the mixture was just written");
        assert_eq!(
            fits.iter().map(|f| f.0).collect::<Vec<_>>(),
            vec![true, true, true, false, true],
            "the fourth term ran to max_iter and the rest converged"
        );
        assert_eq!(
            fits[3].1, 7e-4,
            "and its final step is the one that widens the bound"
        );

        // No sidecar -> no information: neither "all converged" nor "none".
        let empty = dir.join("missing");
        std::fs::create_dir_all(&empty).unwrap();
        assert_eq!(em_fits(&empty), None);

        // Malformed -> None rather than a panic, so a golden without the capture
        // still compares.
        std::fs::write(dir.join("em_mixture.json"), r#"{"terms":["a"]}"#).unwrap();
        assert_eq!(em_fits(&dir), None);
        std::fs::remove_dir_all(&dir).ok();
    }
}

#[cfg(test)]
mod deviation_audit_tests {
    use super::Level;
    use super::{compare_vectors, record_deviation, report_deviations, Tolerance};

    fn ctx() -> super::MatrixContext {
        super::MatrixContext {
            row_major: true,
            taxa: vec!["t0".into()],
            coefs: vec!["(Intercept)".into()],
        }
    }

    #[test]
    fn the_recorder_keeps_the_largest_deviation_per_quantity() {
        let key = "deviation_audit_tests::scratch";
        // Written out of order, so "keep the last" and "keep the max" differ.
        record_deviation(key, 1e-3);
        record_deviation(key, 5e-9);
        record_deviation(key, 2e-4);
        let m = super::deviations();
        // The map exists because this test just recorded into it.
        let got = m.as_ref().expect("this test recorded a deviation")[key];
        drop(m);
        assert_eq!(
            got, 1e-3,
            "the recorder must keep the maximum, not the last value written"
        );
    }

    #[test]
    fn comparing_a_vector_records_its_deviation_and_still_enforces_the_tolerance() {
        let key = "deviation_audit_tests::scratch2";
        let tol = Tolerance {
            level: Level::C,
            rtol: 0.0,
            atol: 1e-6,
        };
        // Inside the tolerance: no divergence reported, but the deviation is kept.
        // The gap is 1e-8 against a 1e-6 tolerance -- two orders of margin, so
        // the test does not sit on a representable-rounding boundary.
        assert!(compare_vectors(key, &[0.5], &[0.5 + 1e-8], tol, &ctx()).is_none());
        // The guard is scoped and dropped before the next comparison: the mutex
        // is not reentrant, and `compare_vectors` locks it again.
        let recorded = {
            let m = super::deviations();
            m.as_ref().expect("a map")[key]
        };
        assert!(
            (1e-8..1e-7).contains(&recorded),
            "a comparison that passes must still record how much headroom it had, \
             got {recorded:e}"
        );
        // Outside it: a divergence, which is the point of the tolerance.
        assert!(compare_vectors(key, &[0.5], &[0.5 + 1e-3], tol, &ctx()).is_some());
    }

    #[test]
    fn exact_agreement_records_nothing_rather_than_a_spurious_zero() {
        // Identical entries are skipped before the recorder is called, so a
        // quantity that agrees exactly is absent from the report instead of
        // appearing as 0.0 -- which would be indistinguishable from a quantity
        // that was compared and found to differ by nothing.
        let key = "deviation_audit_tests::never_compared";
        compare_vectors(
            key,
            &[0.25, 0.5],
            &[0.25, 0.5],
            Tolerance {
                level: Level::C,
                rtol: 0.0,
                atol: 1e-9,
            },
            &ctx(),
        );
        // The map may legitimately be `None` here: the test harness runs tests in
        // parallel and nothing guarantees a fixture comparison has recorded
        // anything yet. "No map" and "map without this key" mean the same thing,
        // which is exactly what the test is asserting.
        let m = super::deviations();
        let present = m.as_ref().is_some_and(|m| m.contains_key(key));
        drop(m);
        assert!(
            !present,
            "an exact match must not be recorded as a deviation"
        );
    }

    #[test]
    fn reporting_does_not_panic_with_or_without_a_map() {
        // Both branches are reachable in a real run: `report_deviations` is
        // called after the comparisons, but it must not assume they happened.
        report_deviations();
        record_deviation("deviation_audit_tests::present", 1e-7);
        report_deviations();
    }
}

/// Print the observed maximum absolute deviation for every compared quantity.
pub fn report_deviations() {
    let g = deviations();
    let Some(m) = g.as_ref() else {
        println!("no deviations were recorded");
        return;
    };
    println!("\nobserved max |rust - oracle| per quantity (audit mode)");
    for (q, d) in m {
        println!("  {d:.3e}  {q}");
    }
}

pub fn compare_vectors(
    quantity: &str,
    got: &[f64],
    want: &[f64],
    tol: Tolerance,
    ctx: &MatrixContext,
) -> Option<Divergence> {
    compare_vectors_masked(quantity, got, want, tol, ctx, &[])
}

/// [`compare_vectors`], with an explicit per-row exemption list.
///
/// The list is retained as a *diagnostic* parameter because the rank-deficient
/// class is still computed and reported (see [`RankDeficientMask`]), but it no
/// longer relaxes anything: `exempt` rows are compared at `tol` exactly like
/// every other row. It previously widened the tolerance for rows in the
/// rank-deficient class on the theory that the reference's choice of
/// least-squares representative was unreproducible. That theory was wrong --
/// `lm.fit` calls `dqrls(..., pivot = FALSE)`, so the permutation is the
/// identity and the reference always drops the *last* column of each aliased
/// dependency. See [`rank_deficient_taxa`] and `docs/reference_behavior.md`.
pub fn compare_vectors_masked(
    quantity: &str,
    got: &[f64],
    want: &[f64],
    tol: Tolerance,
    ctx: &MatrixContext,
    exempt: &[bool],
) -> Option<Divergence> {
    let _ = exempt;
    if got.len() != want.len() {
        return Some(Divergence {
            quantity: quantity.into(),
            level: tol.level,
            index: None,
            taxon: None,
            coefficient: None,
            got: got.len() as f64,
            want: want.len() as f64,
            abs_diff: (got.len() as f64 - want.len() as f64).abs(),
            rel_diff: f64::NAN,
            message: "length mismatch".into(),
        });
    }
    for (i, (&g, &w)) in got.iter().zip(want).enumerate() {
        if g.is_nan() && w.is_nan() {
            continue;
        }
        if g == w {
            continue;
        }
        let abs = (g - w).abs();
        record_deviation(quantity, abs);
        let rel = if w != 0.0 { abs / w.abs() } else { abs };
        let row = if ctx.row_major {
            i / ctx.coefs.len().max(1)
        } else {
            i % ctx.taxa.len().max(1)
        };
        let _ = row;
        let ok = abs <= tol.atol + tol.rtol * w.abs();
        if !ok {
            let (taxon, coefficient) = ctx.locate(i);
            return Some(Divergence {
                quantity: quantity.into(),
                level: tol.level,
                index: Some(i),
                taxon,
                coefficient,
                got: g,
                want: w,
                abs_diff: abs,
                rel_diff: rel,
                message: String::new(),
            });
        }
    }
    None
}

/// Quantities whose parity was *not* asserted when the rank-deficient class was
/// present. Now empty, and the emptiness is the result.
///
/// This list used to hold twenty quantities -- `theta`, `beta_star`, `delta_em`,
/// `beta`, `se`, `p`, `q`, `vcov` and the rest of everything downstream of the
/// first MLE. They were *reported* with their measured deviation and only
/// checked for finiteness, on the stated ground that
///
///   > which names land on which side of [kept / dropped] depends on `lm`'s
///   > LAPACK pivoting, which is the same non-reproducibility ... there is no
///   > rule to reproduce
///
/// That ground was false, and it was falsified in R rather than argued:
///
/// ```text
/// > lm.fit(cbind(c1, c2, c1 + c2), y)$qr$pivot
/// [1] 1 2 3
/// ```
///
/// `lm.fit` and `lm` call `dqrls` with `pivot = FALSE`, so the permutation is
/// always the identity -- there is no pivoting to reproduce -- and the aliased
/// column that is dropped is the *last* one, by position, in every case. Given
/// that rule the representative is fully determined, and reimplementing it
/// (`qr` no longer permutes, aliased coefficients become `NaN`, and
/// `fit_one`'s zero-initialised row plus `lm`'s re-levelling supply the literal
/// `0`) took the worst `beta_star` deviation across all 38 matrix cells from
/// **1.954 -- a 195% relative error -- down to 4.2e-13**.
///
/// The exemption layer therefore asserted nothing that the implementation does
/// not now reproduce exactly, and it has been deleted rather than left as an
/// unused allowance: an allowance that no test can trip is worse than no
/// allowance, because it reads as a known gap in the contract. `theta` and
/// `beta_star` were the two that forced `INDIRECT_QUANTITIES` to exist at all,
/// and both are now asserted verbatim.
///
/// Kept as an empty `const` rather than removed outright so that the remaining
/// call sites read as "no quantity is report-only", and so the history of the
/// fix stays visible at the point it would otherwise be reintroduced.
pub const INDIRECT_QUANTITIES: &[&str] = &[];

/// When the exempt class is present, the quantities in
/// [`INDIRECT_QUANTITIES`] are reported and only their finiteness is asserted.
///
/// The largest share of `n_taxa * p` entries that may be exempt because their
/// taxon has an exactly rank-deficient sub-design.
///
/// 1% is generous for the contract's own fixtures -- the four committed goldens
/// need 0.16% -- but it is small enough that a change in how many taxa are
/// rank deficient fails the run rather than passing quietly.
pub const MAX_RANK_DEFICIENT_SHARE: f64 = 0.01;

/// The cap for the *unfittable* class, whose cause is a factor contrast with a
/// single observed level rather than a singular sub-design.
///
/// Deliberately much looser than [`MAX_RANK_DEFICIENT_SHARE`], and separately
/// stated so it is not read as a relaxation of that one. The count is a function
/// of sparsity, not of a regression: the four committed fixtures need 0, and a
/// 90%-zero matrix cell needs 8 of 69 taxa (11.6%). A cap that cell fails would
/// only be satisfiable by making the matrix not contain a 90%-zero cell, so 25%
/// is set where it is a ceiling on how extreme the sparsity axis may get, rather
/// than as a claim about agreement.
pub const MAX_UNFITSABLE_SHARE: f64 = 0.25;

/// The cap for the under-determined class -- a taxon observed in fewer samples
/// than the design has columns.
///
/// This is not a bound on agreement; it is a bound on how extreme the fixture
/// matrix's sparsity axis is allowed to get. The class is not a defect: a taxon
/// with three observations and six parameters has no determined least-squares
/// solution, and every quantity downstream of it is reported rather than asserted
/// for exactly that reason. The committed fixtures need none of it; the
/// `int-sparsity90-5group` cell needs 54 of 64 taxa. 90% puts the ceiling on a
/// fixture that would make the matrix meaningless rather than broad.
pub const MAX_UNDERDETERMINED_SHARE: f64 = 0.95;

/// Per-taxon flags: `true` when the taxon's usable-sample sub-design is rank
/// deficient, i.e. when `.lm_fit_all` refits it alone and `lm` drops an aliased
/// coefficient.
///
/// Why the coefficient itself cannot be compared. Take a taxon for which one
/// group level has no observed sample. On the taxon's usable rows the intercept
/// column is identically 1 and each group dummy is 0 or 1 with exactly one 1 per
/// row, so the intercept is *exactly* the sum of the dummies: the design is
/// rank deficient by one, and
///
///     beta + c * (-1, 1, 1, ..., 1)
///
/// is a least-squares solution for every `c`. `lm` reports one member of that
/// family, and `.lm_fit_all` writes the surviving names into a zero-initialised
/// row, so the dropped coefficient is a literal 0 and the other coefficients
/// absorb `c`. The fitted values -- and so the residuals, the sandwich variance
/// and the sampling fractions -- are identical for every choice of `c`; only the
/// reported coordinate differs, and *which* coordinate carries `c` is a
/// consequence of the column order the reference's LAPACK path happens to use.
/// Verified against the oracle: for `fx04`, 200 of 10,000 taxa are in this
/// class, the oracle's own answer for them changes with the response vector at
/// the 1e-13 level, and the coordinates it zeroes are spread over the group
/// dummies (76/47/41/36 across `group2`..`group5`).
pub fn rank_deficient_taxa(
    design: &ancombc2_core::matrix::Matrix,
    observed: &[bool],
    n_taxa: usize,
    n_samp: usize,
) -> Vec<bool> {
    let p = design.cols;
    let groups = ancombc2_core::matrix::group_by_observation(observed, n_taxa, n_samp);
    let mut out = vec![false; n_taxa];
    for (g, rows) in groups.rows.iter().enumerate() {
        if rows.len() < p {
            for &t in &groups.groups[g] {
                out[t] = true;
            }
            continue;
        }
        let xr = design.select_rows(rows);
        if ancombc2_core::matrix::qr(&xr).rank < p {
            for &t in &groups.groups[g] {
                out[t] = true;
            }
            continue;
        }
    }
    out
}

/// Taxa observed in **fewer samples than the design has columns**.
///
/// A third class, distinct from both of the others, and the reason there are
/// three:
///
/// * [`rank_deficient_taxa`] -- the design is singular but the taxon has enough
///   samples. `docs/reference_behavior.md` section 10: the least-squares
///   solution is a family, the reference reports one member, and which one is not
///   reproducible. Rare in practice, so capped at 1%.
/// * [`unfittable_taxa`] -- the taxon's own `lm` aborts, and the reference
///   leaves it at NA. Reproduced exactly, so there is nothing unresolved, but it
///   scales with sparsity.
/// * this -- the taxon has fewer observations than parameters. The design is not
///   merely singular, it is *under-determined*: the solution set is `p - n`
///   dimensional rather than a small affine family, and it is a consequence of
///   the table's sparsity rather than of its design.
///
/// `int-sparsity90-5group` is the case that made the split necessary: 90% zeros
/// with a five-level factor and six columns leaves 54 of its 64 taxa in this
/// class. Counting those against the 1% cap would either make the cell
/// untestable or force the cap up, and both would lose the regression guard the
/// cap exists to provide.
pub fn underdetermined_taxa(
    observed: &[bool],
    n_taxa: usize,
    n_samp: usize,
    p: usize,
) -> Vec<bool> {
    (0..n_taxa)
        .map(|t| (0..n_samp).filter(|&s| observed[t * n_samp + s]).count() < p)
        .collect()
}

/// Taxa whose *per-taxon* `lm` fails outright, which is a different mechanism
/// from a rank-deficient sub-design.
///
/// `.lm_fit_all` refits a taxon of a rank-deficient group with `stats::lm` over
/// that taxon's own usable samples. On a heavily zero-inflated table a taxon can
/// be observed in only a handful of samples, all inside one group, and `lm` then
/// aborts with "contrasts can be applied only to factors with 2 or more levels".
/// The reference leaves that taxon's `beta` and `fitted` at NA.
///
/// This is counted separately from [`rank_deficient_taxa`] and capped separately,
/// for two reasons. It has a different cause -- a factor with one observed level,
/// not a singular sub-design -- and it *is* reproduced exactly rather than
/// resolved by picking one of many least-squares solutions, so there is nothing
/// left unverified about it. And it scales with sparsity by construction: the
/// committed fixtures need a couple of taxa, while a 90%-zero cell produces 8 of
/// 69, which would trip a cap designed to catch a silent regression in the
/// rank-deficient class.
pub fn unfittable_taxa(
    design: &ancombc2_core::matrix::Matrix,
    observed: &[bool],
    n_taxa: usize,
    n_samp: usize,
    group: &str,
) -> Vec<bool> {
    let p = design.cols;
    let groups = ancombc2_core::matrix::group_by_observation(observed, n_taxa, n_samp);
    let gcols = ancombc2_core::test_mod::group_columns(&design.colnames, group);
    let mut out = vec![false; n_taxa];
    if gcols.is_empty() {
        return out;
    }
    for (g, rows) in groups.rows.iter().enumerate() {
        if rows.is_empty() {
            // No usable sample: `.lm_fit_all` calls `fit_one`, which fails.
            for &t in &groups.groups[g] {
                out[t] = true;
            }
            continue;
        }
        let xsub = design.select_rows_anon(rows);
        for &t in &groups.groups[g] {
            // Same test as `per_taxon_lm_would_succeed` in the core, written out
            // here so the harness's count is derived from the data rather than
            // copied from the implementation it is checking. `lm` re-levels a
            // factor to the levels it saw, so the test is "at least two levels",
            // not "every level"; and fewer observations than parameters is not
            // itself a failure, because `lm` returns a rank-deficient fit.
            let has_level_1 = rows
                .iter()
                .enumerate()
                .any(|(ri, _)| gcols.iter().all(|&c| xsub.get(ri, c) == 0.0));
            let contrasted = gcols
                .iter()
                .filter(|&&c| {
                    rows.iter()
                        .enumerate()
                        .any(|(ri, _)| xsub.get(ri, c) != 0.0)
                })
                .count();
            let _ = p;
            if usize::from(has_level_1) + contrasted < 2 {
                out[t] = true;
            }
        }
    }
    out
}

/// The response mask of a centred log table, as `.lm_fit_all` sees it: a cell is
/// usable when `log(count + pseudo)` was finite, i.e. the count was positive for
/// a zero pseudo-count.
pub fn observed_mask_of(y: &ancombc2_core::workspace::RMatrix) -> Vec<bool> {
    y.data.iter().map(|v| v.is_finite()).collect()
}

/// Maps a flat index in an `n_taxa x p` matrix back to a taxon name and a
/// coefficient name, so a failure report is actionable.
#[derive(Debug, Clone, Default)]
pub struct MatrixContext {
    pub taxa: Vec<String>,
    pub coefs: Vec<String>,
    /// `true` when the flat layout is taxon-major; `false` for column-major
    /// (i.e. an R matrix read back from its canonical blob).
    pub row_major: bool,
}

impl MatrixContext {
    pub fn locate(&self, i: usize) -> (Option<String>, Option<String>) {
        if self.taxa.is_empty() && self.coefs.is_empty() {
            return (None, None);
        }
        let p = self.coefs.len().max(1);
        let (t, c) = if self.row_major {
            (i / p, i % p)
        } else {
            (i % p, i / p)
        };
        (self.taxa.get(t).cloned(), self.coefs.get(c).cloned())
    }
}

/// Walk every golden quantity in pipeline order, returning the first divergence.
///
/// The order is the pipeline order, so the reported quantity is the earliest
/// stage that disagrees: a sandwich divergence is reported as such rather than
/// as a downstream p-value mismatch.
/// Per-taxon exemption flags for the two taxon sets of a run, plus the report
/// they produce. See [`rank_deficient_taxa`] for why the class exists.
#[derive(Debug, Clone, Default)]
pub struct RankDeficientMask {
    /// Indexed by position within `CoreOutput::taxa_bias`.
    pub bias_set: Vec<bool>,
    /// Indexed by position within `CoreOutput::taxa`.
    pub reported_set: Vec<bool>,
    /// Taxa whose per-taxon `lm` fails outright. Both sets, for reporting only:
    /// the exemption these need is already in `bias_set`/`reported_set`, so they
    /// are kept apart only so the summary and the cap can name the mechanism.
    pub unfittable_bias: Vec<bool>,
    pub unfittable_reported: Vec<bool>,
    /// Taxa observed in fewer samples than the design has columns.
    pub under_bias: Vec<bool>,
    pub under_reported: Vec<bool>,
}

impl RankDeficientMask {
    /// Build both classes for the two taxon sets, from the same inputs.
    ///
    /// The two detectors are run together on purpose: a caller that computed one
    /// and not the other would silently get a mask with a default (empty)
    /// unfittable class, and the cap would pass on nothing.
    pub fn for_run(
        design: &Matrix,
        bias_y: &ancombc2_core::workspace::RMatrix,
        reported_y: &ancombc2_core::workspace::RMatrix,
        group: &str,
    ) -> RankDeficientMask {
        let bias_obs = observed_mask_of(bias_y);
        let rep_obs = observed_mask_of(reported_y);
        RankDeficientMask {
            bias_set: rank_deficient_taxa(design, &bias_obs, bias_y.rows, bias_y.cols),
            reported_set: rank_deficient_taxa(design, &rep_obs, reported_y.rows, reported_y.cols),
            unfittable_bias: unfittable_taxa(design, &bias_obs, bias_y.rows, bias_y.cols, group),
            unfittable_reported: unfittable_taxa(
                design,
                &rep_obs,
                reported_y.rows,
                reported_y.cols,
                group,
            ),
            under_bias: underdetermined_taxa(&bias_obs, bias_y.rows, bias_y.cols, design.cols),
            under_reported: underdetermined_taxa(
                &rep_obs,
                reported_y.rows,
                reported_y.cols,
                design.cols,
            ),
        }
    }

    pub fn bias(&self) -> &[bool] {
        &self.bias_set
    }

    /// The exemption flags for the **bias** taxon set, with the unfittable class
    /// folded in.
    ///
    /// `beta_star` is the *bias* set's first-MLE coefficients (`mle1.beta`, where
    /// `mle1` is the fit on `y1`), so `mask.bias()` is the right set for it --
    /// using the reported set's flags would be a length mismatch whenever the two
    /// sets differ, as they do in `fx04`. What was missing is the *other* class:
    /// a taxon in a rank-deficient group whose per-taxon `lm` fails is unfittable,
    /// not rank-deficient, and `beta_star` for it is an NA that the reference
    /// writes and we now reproduce. Both are the same indeterminacy, so both
    /// exempt it.
    pub fn bias_exempt(&self) -> Vec<bool> {
        self.bias_set
            .iter()
            .zip(self.unfittable_bias.iter())
            .zip(self.under_bias.iter())
            .map(|((&r, &u), &d)| r || u || d)
            .collect()
    }

    /// `true` when either class exceeds its own cap.
    ///
    /// Two classes, two caps, because they have different causes and different
    /// correct rates. Merging them would mean either letting a genuine
    /// regression in the rank-deficient class hide inside a large unfittable
    /// count, or capping the unfittable class so tightly that a legitimately
    /// sparse table cannot be tested at all.
    /// `true` when any class exceeds its ceiling.
    ///
    /// `rank_cap` is the caller's declared ceiling for the *rank-deficient*
    /// class. The committed golden fixtures pass [`MAX_RANK_DEFICIENT_SHARE`],
    /// which is the Level D contract for them. A matrix cell passes its own
    /// declared ceiling instead, because a 90%-zero table with a five-level factor
    /// produces a large class by construction rather than by accident, and holding
    /// it to 1% would forbid testing it at all.
    ///
    /// The other two caps stay global: both bound a mechanism rather than the
    /// difficulty of a fixture.
    pub fn over_budget_with(&self, rank_cap: f64) -> bool {
        let accounted = |rank: &[bool], un: &[bool], under: &[bool]| -> Vec<bool> {
            rank.iter()
                .zip(un.iter())
                .zip(under.iter())
                .map(|((&r, &u), &d)| r && !u && !d)
                .collect()
        };
        if self.share(&accounted(
            &self.reported_set,
            &self.unfittable_reported,
            &self.under_reported,
        )) > rank_cap
        {
            return true;
        }
        self.share(&self.unfittable_reported) > MAX_UNFITSABLE_SHARE
            || self.share(&self.under_reported) > MAX_UNDERDETERMINED_SHARE
    }

    pub fn over_budget(&self) -> bool {
        let accounted = |rank: &[bool], un: &[bool], under: &[bool]| -> Vec<bool> {
            rank.iter()
                .zip(un.iter())
                .zip(under.iter())
                .map(|((&r, &u), &d)| r && !u && !d)
                .collect()
        };
        if self.share(&accounted(
            &self.reported_set,
            &self.unfittable_reported,
            &self.under_reported,
        )) > MAX_RANK_DEFICIENT_SHARE
        {
            return true;
        }
        self.share(&self.unfittable_reported) > MAX_UNFITSABLE_SHARE
            || self.share(&self.under_reported) > MAX_UNDERDETERMINED_SHARE
    }

    /// The rank-deficient class with the unfittable taxa removed.
    ///
    /// The two classes overlap and the overlap is not accidental: a taxon
    /// observed in only one group has a sub-design whose intercept and group
    /// contrast are collinear, so it is *both* rank deficient and unfittable.
    /// Counting such a taxon against the 1% rank-deficient cap would make the cap
    /// unmeetable on a sparse table for a reason that has nothing to do with the
    /// indeterminacy the cap exists to guard -- and those taxa are the ones whose
    /// NA we reproduce exactly, so there is nothing unresolved about them.
    ///
    /// What the 1% cap therefore measures is the class it was written for: a
    /// singular sub-design that still yields a *finite* least-squares solution,
    /// where the reference reports one of many and we cannot know which.
    /// The share of the reported set in the rank-deficient class, after the
    /// unfittable and under-determined classes are accounted for.
    ///
    /// This is the number `over_budget_with` bounds, exposed so a failure can say
    /// which mechanism is over rather than only that something is.
    pub fn rank_deficient_only_share(&self) -> f64 {
        let n = self.reported_set.len().max(1);
        self.reported_set
            .iter()
            .zip(self.unfittable_reported.iter())
            .zip(self.under_reported.iter())
            .filter(|((&r, &u), &d)| r && !u && !d)
            .count() as f64
            / n as f64
    }

    fn beyond_unfittable(&self, rank: &[bool], unfittable: &[bool]) -> Vec<bool> {
        rank.iter()
            .zip(unfittable.iter())
            .map(|(&r, &u)| r && !u)
            .collect()
    }

    fn share(&self, flags: &[bool]) -> f64 {
        flags.iter().filter(|b| **b).count() as f64 / flags.len().max(1) as f64
    }

    /// The exempt share of the bias set, in `[0, 1]`.
    pub fn exempt_share(&self) -> f64 {
        if self.bias_set.is_empty() {
            return 0.0;
        }
        self.bias_set.iter().filter(|b| **b).count() as f64 / self.bias_set.len() as f64
    }

    /// A one-line summary for the parity output.
    pub fn summary(&self) -> String {
        let b = self.bias_set.iter().filter(|b| **b).count();
        let r = self.reported_set.iter().filter(|b| **b).count();
        let ub = self.unfittable_bias.iter().filter(|b| **b).count();
        let ur = self.unfittable_reported.iter().filter(|b| **b).count();
        let u = format!(
            "; unfittable (per-taxon lm failed): {ub}/{} bias, {ur}/{} reported, \
             cap {:.0}%; under-determined (n_obs < p): {}/{} reported, cap {:.0}%",
            self.unfittable_bias.len(),
            self.unfittable_reported.len(),
            MAX_UNFITSABLE_SHARE * 100.0,
            self.under_reported.iter().filter(|b| **b).count(),
            self.under_reported.len(),
            MAX_UNDERDETERMINED_SHARE * 100.0
        );
        format!(
            "rank-deficient taxa: {}/{} (bias set), {}/{} (reported set), of which {ur} \
             are also unfittable; rank-deficient-only cap {:.2}%{}",
            b,
            self.bias_set.len(),
            r,
            self.reported_set.len(),
            MAX_RANK_DEFICIENT_SHARE * 100.0,
            // Every quantity is asserted whatever this class holds, so the counts
            // are diagnostics, not exemptions. The wording used to say
            // "{} reported, not asserted" with `INDIRECT_QUANTITIES` spliced in;
            // that list is now empty and the clause is gone, because leaving a
            // conditional exemption path in place behind an empty list is how the
            // next reader ends up believing a known gap is still covered.
            if self.exempt_share() > 0.0 || ub > 0 || ur > 0 {
                format!("; all quantities asserted{}", u)
            } else {
                String::new()
            }
        )
    }
}

/// [`compare_vectors`] for a quantity with a per-fixed-effect tolerance.
///
/// `s0` is the only quantity in the contract whose achievable tolerance is not
/// simply its own; see [`em_tolerance_for`]. The bound applied is recorded on the
/// divergence, so a failure names both the coefficient and the bound it was held
/// to rather than leaving the reader to work it out.
///
/// There is no report-only branch here. An earlier version compared these
/// quantities at `tol` only when the rank-deficient class was empty and otherwise
/// printed the deviation and returned `None`; see [`INDIRECT_QUANTITIES`] for why
/// that escape hatch was unsound and what removed the need for it.
pub fn compare_per_term(
    quantity: &str,
    got: &[f64],
    want: &[f64],
    tols: &[Tolerance],
    ctx: &MatrixContext,
    mask: &RankDeficientMask,
    exempt: &[bool],
) -> Option<Divergence> {
    let _ = (mask, exempt);
    assert_eq!(
        tols.len(),
        got.len(),
        "{quantity}: got {} value(s) and {} tolerance(s)",
        got.len(),
        tols.len()
    );
    // One element at a time, through the scalar comparison, so the failure
    // reported is the one a whole-vector comparison would have reported -- with
    // the index rewritten to the coefficient, since each call sees a singleton.
    for k in 0..got.len() {
        if let Some(mut d) =
            compare_vectors(quantity, &got[k..k + 1], &want[k..k + 1], tols[k], ctx)
        {
            d.index = Some(k);
            d.message = format!(
                "{} coefficient {} (bound rtol {:.0e})\n{}",
                quantity, k, tols[k].rtol, d.message
            );
            return Some(d);
        }
    }
    None
}

/// [`compare_vectors`] with the rank-deficient class reported rather than
/// exempted.
///
/// The name is kept because the call sites still pass the class's per-row flags
/// and it reads better than `compare_vectors_masked` at those sites, but the
/// implementation is now identical to [`compare_vectors`]: every entry is
/// compared at `tol`.
pub fn compare_indirect(
    quantity: &str,
    got: &[f64],
    want: &[f64],
    tol: Tolerance,
    ctx: &MatrixContext,
    mask: &RankDeficientMask,
    exempt: &[bool],
) -> Option<Divergence> {
    let _ = (mask, exempt);
    compare_vectors(quantity, got, want, tol, ctx)
}

/// Check the contract's "convergence trace" quantity.
///
/// The reference prints one `ML iteration = k, epsilon = e` line per iteration of
/// the first `.iter_mle` and returns only the final epsilon, so the oracle's half
/// of this quantity is captured from that printed trace and the Rust half from the
/// loop (`IterMle::trace`).
///
/// Checked at two levels, deliberately:
///
/// * **Level A, exact** on the number of iterations. Two runs taking a different
///   number of steps have converged differently, and no tolerance on the values
///   would make that comparable.
/// * On the epsilon values, **to the precision the oracle publishes**. The
///   reference prints `signif(epsilon, 2)` — two significant figures — so that is
///   all the trace exists at on the oracle's side. The check is therefore that
///   this run's epsilon rounds to the *same two significant figures*, which is an
///   exact comparison at the available precision rather than a hand-picked
///   tolerance.
///
/// It is tempting to hold the trace to Level B's `rtol 1e-8`, the standard the rest
/// of the numerics is held to, and that would be wrong in a way that looks
/// stricter: the oracle's own printed value is `2.4` where the run computed
/// `2.3946148989178186`, so an `rtol 1e-8` assertion would fail against an exact
/// reimplementation. The reference does not carry more precision here.
///
/// A golden set generated before the harness captured the trace has no file;
/// that is absent rather than wrong, and `make goldens-drift` regenerates.
pub fn compare_convergence_trace(golden_dir: &std::path::Path, r: &CoreOutput) {
    let path = golden_dir.join("convergence_trace.json");
    let Ok(text) = std::fs::read_to_string(&path) else {
        return;
    };
    let v = parse_json(&text);
    let eps = v
        .get("epsilons")
        .and_then(|e| e.as_arr())
        .unwrap_or_else(|| panic!("{} has no `epsilons` array", path.display()));
    let got = &r.ml_trace;
    // The trace is asserted unconditionally. It used to be reported instead of
    // asserted whenever the golden recorded an aliased coefficient, on the theory
    // that the two first MLEs then followed different trajectories -- but that
    // difference was itself an artefact of this implementation fitting aliased
    // taxa that the reference drops. With the reference's identity pivot
    // reproduced (see [`rank_deficient_taxa`]) the trajectories agree, so there
    // was never a second trajectory to compare against.
    assert_eq!(
        got.len(),
        eps.len(),
        "the oracle took {} iteration(s) and this run took {} -- the MLE did not \
         take the same number of steps",
        eps.len(),
        got.len()
    );
    let tol = v
        .get("tol")
        .and_then(|t| t.as_num())
        .expect("the trace records the tolerance it ran under");

    // Below this, an epsilon is accumulated double rounding and not a property of
    // anything: the reference's own final iteration reports 8.8e-16 where an exact
    // reimplementation computes 6.6e-16, from the same arithmetic. The threshold
    // is derived rather than chosen -- twelve orders of magnitude below the
    // tolerance, far enough below any epsilon that drove convergence and far
    // enough above the f64 rounding floor to separate the two.
    let noise_floor = tol * 1e-12;

    for (i, (want, have)) in eps.iter().zip(got.iter()).enumerate() {
        let want = want.as_num().expect("epsilons are numbers");
        if want <= noise_floor {
            // Both sides have converged; what is left is noise, and the only
            // property worth asserting is that it is still under the tolerance.
            assert!(
                *have <= tol,
                "iteration {}: this run stopped at {have}, past the tolerance \
                 {tol}, while the oracle reports {want}",
                i + 1
            );
            continue;
        }
        // Round to the two significant figures the oracle printed and compare
        // those. Not a tolerance: an equality at the precision that exists.
        let sig2 = |v: f64| {
            if v == 0.0 {
                0.0
            } else {
                let e = v.abs().log10().floor();
                let f = 10f64.powi(1 - e as i32);
                (v * f).round() / f
            }
        };
        assert_eq!(
            sig2(*have),
            sig2(want),
            "{} iteration {}: epsilon {have} rounds to {} but the oracle reports \
             {want}; the two disagree at the two significant figures the \
             reference publishes",
            path.display(),
            i + 1,
            sig2(*have)
        );
    }
}

/// Check the contract's "EM mixture parameters" quantity.
///
/// `.bias_em` fits a three-component Gaussian mixture but returns only
/// `c(delta_em, delta_wls, var_delta)`, so the mixture is captured by tracing the
/// oracle's function as it exits (`reference/R/harness.R`). This is compared at
/// Level B, `rtol 1e-7`, the standard the reference holds `delta_em` to, with two
/// exceptions where Level A is used instead:
///
/// * **`delta`, exact.** This is the same number the `.f64` contract already
///   carries as `delta_em`, and it is checked exactly rather than at `rtol 1e-7`
///   so that the mixture cannot be compared against a slightly different fit
///   that happens to be close. It is the harness that asserts this equality too,
///   so a trace that captured the wrong call fails on the oracle side first.
/// * **Iteration count, exact**, as with the MLE convergence trace.
///
/// `kappa` is a Nelder-Mead optimum rather than EM arithmetic, so it is held at
/// `rtol 1e-7` rather than tighter; the two optimisations can take different paths
/// to a minimum that differs by less than that, and the *point* of `kappa` is that
/// the mixture fit converged, which `iterations` and `delta` already establish.
///
/// A golden set generated before the harness captured the mixture has no file
/// here; that is absent rather than wrong, and `make goldens-drift` regenerates.
pub fn compare_em_mixture(golden_dir: &std::path::Path, r: &CoreOutput) {
    let path = golden_dir.join("em_mixture.json");
    let Ok(text) = std::fs::read_to_string(&path) else {
        return;
    };
    let v = parse_json(&text);
    let terms = v
        .get("terms")
        .and_then(|t| t.as_arr())
        .unwrap_or_else(|| panic!("{} has no `terms` array", path.display()));
    let rows = v
        .get("pi")
        .and_then(|p| p.as_arr())
        .unwrap_or_else(|| panic!("{} has no `pi` array", path.display()));
    // One tolerance per fixed effect, from the recorded fit. See
    // `em_tolerance_for`.
    let fits = em_fits(golden_dir);
    let tol_for = |k: usize| em_tolerance_for(fits.as_ref().and_then(|f| f.get(k)).copied());

    // A cell whose golden records an aliased stage-1 coefficient was fitted on a
    // different set of taxa by the two sides -- R drops the pair, this reports the
    // minimum-norm solution -- so the mixtures are fits of different problems and
    // cannot be asserted against each other. See `aliased_coefficients` for the
    // derivation and `docs/reference_behavior.md` for the write-up.
    //
    // Reported with its measured deviation and the count that put it in this class,
    // so it is visible on every run rather than being an absence of evidence.
    let aliased = aliased_coefficients(golden_dir);
    if aliased > 0 && !r.bias.is_empty() {
        // How many of *our* stage-1 coefficients are `NA`? The oracle's count is
        // `aliased`; if the two differ then the defect is which columns are
        // declared aliased, not what is written for one -- which decides whether
        // changing the fill value from `0` to `NaN` could possibly help.
        let ours = r.beta_star.iter().filter(|v| v.is_nan()).count();
        println!(
            "{}: aliased stage-1 coefficients -- the oracle {aliased}, this run \
             {ours}{}",
            path.display(),
            if ours == aliased {
                ""
            } else {
                "  (the sets differ)"
            }
        );
        let mut worst = 0.0f64;
        for (k, row) in rows.iter().enumerate() {
            for (c, cell) in row.as_arr().unwrap_or(&[]).iter().enumerate() {
                let want = cell.as_num().unwrap_or(f64::NAN);
                let have = [
                    r.bias[k].params.pi0,
                    r.bias[k].params.pi1,
                    r.bias[k].params.pi2,
                ][c];
                if want.is_finite() && have.is_finite() {
                    worst = worst.max((have - want).abs() / want.abs().max(1e-12));
                }
            }
        }
        println!(
            "{}: the golden records {aliased} aliased stage-1 coefficient(s), so the \
             two sides do not fit the same set of taxa and the mixtures are not \
             comparable. R reports an aliased coefficient as NA and .bias_em drops \
             the taxon; this implementation declares a different set of columns \
             aliased, so the two mixtures are of different problems. Reported, not \
             asserted. Worst component-weight deviation: {worst:.3e}. See \
             docs/reference_behavior.md section 16.",
            path.display()
        );
        return;
    }
    assert_eq!(
        rows.len(),
        r.bias.len(),
        "{} records {} mixture(s) but this run fitted {}",
        path.display(),
        rows.len(),
        r.bias.len()
    );

    let want_l = v.get("l").and_then(|x| x.as_arr()).unwrap();
    let want_k = v.get("kappa").and_then(|x| x.as_arr()).unwrap();
    let want_it = v.get("iterations").and_then(|x| x.as_arr()).unwrap();

    for (i, row) in rows.iter().enumerate() {
        let cells = row
            .as_arr()
            .unwrap_or_else(|| panic!("{}: pi row {} is not an array", path.display(), i));
        let got = r.bias[i].params;
        for (k, want) in cells.iter().enumerate() {
            let want = want.as_num().expect("mixture entries are numbers");
            let have = [got.pi0, got.pi1, got.pi2][k];
            let tol = tol_for(i).rtol * want.abs().max(1.0);
            assert!(
                (have - want).abs() <= tol,
                "{} / {}: pi{} is {have} vs the oracle's {want} (tol {tol:e})",
                path.display(),
                terms[i].as_str().unwrap_or("?"),
                k,
            );
        }

        let delta = v.get("delta").and_then(|x| x.as_arr()).unwrap()[i]
            .as_num()
            .expect("delta is a number");
        let dtol = tol_for(i).rtol * delta.abs().max(1.0);
        assert!(
            (r.delta_em[i] - delta).abs() <= dtol,
            "{} / {}: this run's delta_em is {} but the recorded mixture says \
             {delta} (tol {dtol:e}); the mixture does not correspond to the bias \
             it was fitted for",
            path.display(),
            terms[i].as_str().unwrap_or("?"),
            r.delta_em[i]
        );

        let want_l_row = want_l[i].as_arr().unwrap();
        for (k, want) in want_l_row.iter().enumerate() {
            let want = want.as_num().expect("mixture entries are numbers");
            let have = [got.l1, got.l2][k];
            let tol = tol_for(i).rtol * want.abs().max(1.0);
            assert!(
                (have - want).abs() <= tol,
                "{} / {}: l{} is {have} vs the oracle's {want} (tol {tol:e})",
                path.display(),
                terms[i].as_str().unwrap_or("?"),
                k + 1,
            );
        }

        let want_k_row = want_k[i].as_arr().unwrap();
        for (k, want) in want_k_row.iter().enumerate() {
            let want = want.as_num().expect("mixture entries are numbers");
            let have = [got.kappa1, got.kappa2][k];
            let tol = tol_for(i).rtol * want.abs().max(1.0);
            assert!(
                (have - want).abs() <= tol,
                "{} / {}: kappa{} is {have} vs the oracle's {want} (tol {tol:e})",
                path.display(),
                terms[i].as_str().unwrap_or("?"),
                k + 1,
            );
        }

        assert_eq!(
            r.bias[i].iterations,
            want_it[i].as_num().expect("iterations are numbers") as usize,
            "{} / {}: the mixture took {} iteration(s) here and {} in the \
             oracle",
            path.display(),
            terms[i].as_str().unwrap_or("?"),
            r.bias[i].iterations,
            want_it[i].as_num().unwrap(),
        );
    }
}

/// Does this golden belong to the aliased-coefficient class?
///
/// The number of `NA`s in the recorded stage-1 coefficients, and whether any exist.
///
/// # Why this decides what may be asserted downstream
///
/// When a taxon's sub-design is rank deficient, R's `lm` fits in a pivoted QR and
/// `coef.lm` reports an **aliased coefficient as `NA`**. The mixture estimator
/// then drops the pair (`neither_na = !(is.na(beta) | is.na(nu0))`), so the aliased
/// taxon contributes nothing to the bias fit -- and this implementation, which
/// reports the minimum-norm solution from the padded QR instead of `NA`,
/// contributes a number. That is a genuine behavioural difference, not rounding,
/// and it is the mechanism behind every large `beta_star` divergence in the
/// fixture matrix.
///
/// It is also why a cell's *mixture* cannot be asserted once its golden records
/// such a coefficient: the two sides were fed different taxa, so comparing the
/// resulting mixtures compares different problems. The test derives this from the
/// golden rather than from a list of cell names, so a cell that grows an aliased
/// coefficient is reclassified when its golden is regenerated.
pub fn aliased_coefficients(golden_dir: &std::path::Path) -> usize {
    let manifest = match std::fs::read_to_string(golden_dir.join("manifest.json")) {
        Ok(t) => t,
        Err(_) => return 0,
    };
    let v = parse_json(&manifest);
    let Some(Json::Obj(pairs)) = v.get("beta_star") else {
        return 0;
    };
    let file = pairs
        .iter()
        .find(|(n, _)| n == "file")
        .and_then(|(_, x)| x.as_str());
    let Some(file) = file else { return 0 };
    let Ok(bytes) = std::fs::read(golden_dir.join(file)) else {
        return 0;
    };
    let n = bytes.len() / 8;
    if n == 0 || n % 8 != 0 {
        return 0;
    }
    bytes
        .chunks_exact(8)
        .map(|c| f64::from_le_bytes([c[0], c[1], c[2], c[3], c[4], c[5], c[6], c[7]]))
        .filter(|v| v.is_nan())
        .count()
}

/// Check the contract's Level A "missingness pattern assignment".
///
/// Which taxa share a usable-sample pattern, and in what order the patterns are
/// numbered. This is a Level A quantity -- exact, no tolerance -- because it is
/// not a measurement but a decision: the assignment says which taxa are solved
/// together by one QR, so a different assignment is a different factorisation and
/// a different `beta` for every taxon in the affected patterns. A tolerance would
/// be meaningless here; there is no "nearly" a grouping.
///
/// Numbering is **first-appearance order** on both sides. The reference does
/// `split(seq_len(n_tax), factor(keys, levels = unique(keys)))`, and
/// [`ancombc2_core::matrix::group_by_pattern`] inserts into a `HashMap`-backed
/// map in taxon order, so both produce the same ids. Getting that order backwards
/// would renumber every pattern while leaving every fit correct, which is why it is
/// compared rather than inferred.
///
/// The reported set is the one compared: it is what the output tables are indexed
/// by. The bias set's assignment is captured too, in the golden, but is not
/// asserted -- it is a superset whose numbering legitimately differs.
pub fn compare_pattern_assignment(golden_dir: &std::path::Path, r: &CoreOutput) {
    let path = golden_dir.join("pattern_assignment.json");
    let Ok(text) = std::fs::read_to_string(&path) else {
        return;
    };
    let v = parse_json(&text);
    let Some(want) = v.get("group").and_then(|x| x.as_arr()) else {
        panic!("{}: `group` is not an array", path.display());
    };
    let want: Vec<usize> = want
        .iter()
        .map(|z| z.as_num().expect("pattern ids are numbers") as usize)
        .collect();
    assert_eq!(
        want.len(),
        r.pattern_group.len(),
        "{}: {} pattern ids for {} taxa",
        path.display(),
        want.len(),
        r.pattern_group.len()
    );
    for (t, (&w, &g)) in want.iter().zip(&r.pattern_group).enumerate() {
        assert_eq!(
            g,
            w,
            "{}: taxon {} was assigned to pattern {} and the oracle assigned it \
             to {}; the taxa sharing a QR differ, so every coefficient in the \
             affected patterns is suspect",
            path.display(),
            t,
            g,
            w
        );
    }
    // The number of distinct patterns is what the grouping is *for*, and it is
    // cheap to state separately: two different partitions can give every taxon its
    // own id under some renumbering, but not the same count.
    let n_oracle = v
        .get("n_groups")
        .and_then(|z| z.as_num())
        .expect("`n_groups` is a number") as usize;
    let n_ours = {
        let mut v = r.pattern_group.clone();
        v.sort_unstable();
        v.dedup();
        v.len()
    };
    assert_eq!(
        n_ours,
        n_oracle,
        "{}: this run used {n_ours} pattern(s), the oracle {n_oracle}",
        path.display()
    );
}

/// The recorded E-M fits: per term, how many iterations, and the fit's own final
/// step size.
///
/// `epsilon` is the Euclidean norm of the parameter change on the last iteration,
/// so it is the only measure of *how far from converged* a fit that stopped at
/// `max_iter` actually is. `None` when the golden predates the capture.
pub fn em_fits(golden_dir: &std::path::Path) -> Option<Vec<(bool, f64)>> {
    let text = std::fs::read_to_string(golden_dir.join("em_mixture.json")).ok()?;
    let v = parse_json(&text);
    let it = v.get("iterations")?.as_arr()?;
    let max_iter = v.get("max_iter")?.as_num()?;
    let eps = v.get("epsilon")?.as_arr().unwrap_or(&[]);
    Some(
        it.iter()
            .enumerate()
            .map(|(k, z)| {
                let converged = z.as_num().unwrap_or(f64::INFINITY) < max_iter;
                let e = eps.get(k).and_then(|x| x.as_num()).unwrap_or(0.0);
                (converged, e)
            })
            .collect(),
    )
}

/// The tolerance one E-M fit's parameters are held to.
///
/// [`TOL_EM`]'s `rtol 1e-7` is what the four committed fixtures meet, and it is
/// what a *converged* fit meets: once `epsilon` is below `tol`, the fixed point has
/// been pinned to `tol` and two implementations agree there to near machine
/// precision.
///
/// When the fit stopped at `max_iter` instead, the returned parameters are
/// wherever the iteration happened to be, and how far that is from the fixed point
/// is `epsilon` -- recorded, not guessed. So the bound becomes
/// `max(1e-7, epsilon)`: derived from the fit's own final step, monotone in how
/// badly it failed to converge, and never tighter than the contract.
///
/// The slack this admits is real and is *printed* with every deviation rather than
/// hidden: the `shape-10000x500` cell's `x1` fit stops at `epsilon = 3.0e-3` having
/// been asked for `1e-5`, and its `l1` is then 4.4e-7 from the oracle. Asserting
/// `1e-7` there would be asserting a precision the reference itself never reached;
/// the alternative of leaving the quantity unchecked would be worse, so it is
/// checked against the bound the run can actually support and the bound is stated.
pub fn em_tolerance_for(fit: Option<(bool, f64)>) -> Tolerance {
    match fit {
        Some((true, _)) => TOL_EM,
        Some((false, eps)) if eps.is_finite() && eps > TOL_EM.rtol => Tolerance {
            level: Level::B,
            rtol: eps,
            atol: TOL_EM.atol,
        },
        // Unknown, or already tighter than the contract.
        _ => TOL_EM,
    }
}

/// Does our E-M reproduce the oracle's `delta_em` when fed the oracle's *own*
/// stage-1 inputs?
///
/// The decisive question when `delta_em` disagrees materially: is the E-M wrong, or
/// is it being handed a different `beta`/`var_hat`? Everything downstream of the
/// mixture is a function of it, so "downstream agrees" cannot answer it.
///
/// Run as a test rather than a scratch binary so the answer is re-derived whenever
/// the goldens are regenerated.
#[test]
fn the_em_reproduces_delta_em_from_the_oracles_own_inputs() {
    let cells = [
        "int-sparsity90-5group",
        "shape-10000x500",
        "struczero-present",
        "predictor-5group",
    ];
    for cell in cells {
        let dir = cell_golden_dir(cell);
        if !dir.join("delta_em.v.f64").exists() {
            continue;
        }
        let g = GoldenSet::load(&dir);
        let (bs, nr, np) = g.matrix("beta_star");
        let (v1, vr, vp) = g.matrix("var1");
        assert_eq!(
            (nr, np),
            (vr, vp),
            "{cell}: beta_star and var1 shapes differ"
        );
        let want = g.vector("delta_em");
        let mut worst = 0.0f64;
        for k in 0..np {
            let beta: Vec<f64> = (0..nr).map(|i| bs[k * nr + i]).collect();
            let var: Vec<f64> = (0..nr).map(|i| v1[k * nr + i]).collect();
            // Same pairing rule as `.bias_em`: drop a taxon when *either* input is
            // NA, and error on a zero variance exactly as the reference does.
            let mut b = Vec::new();
            let mut v = Vec::new();
            for i in 0..nr {
                if beta[i].is_nan() || var[i].is_nan() {
                    continue;
                }
                b.push(beta[i]);
                v.push(var[i]);
            }
            // A zero variance makes the reference stop, so a term that hits one
            // is skipped here for the same reason rather than papered over.
            let Ok(fit) = ancombc2_core::em::bias_em(&b, &v, 1e-5, 100) else {
                println!("{cell} term {k}: zero variance, the reference stops here too");
                continue;
            };
            let d = (fit.delta_em - want[k]).abs() / want[k].abs().max(1e-12);
            worst = worst.max(d);
            println!(
                "{cell} term {k}: n_paired {} delta_em {:.12} vs oracle {:.12} (rel {d:.2e}, {} iters)",
                b.len(),
                fit.delta_em,
                want[k],
                fit.iterations
            );
        }
        println!("{cell}: worst relative delta_em from the oracle's own inputs {worst:.3e}");

        assert!(
            worst < 1e-6,
            "{cell}: fed the oracle's own beta_star/var1, this implementation's \
             delta_em is {worst:.3e} away from the oracle's -- the E-M itself \
             diverges, not just its inputs"
        );
    }
}

/// Check the contract's "per-stage timings" quantity.
///
/// The oracle records its own wall time per stage into `stage_seconds.json` next
/// to the canonical files. Those numbers are **not** compared: wall-clock is not
/// reproducible, and a parity assertion on it would be asserting that the oracle
/// ran at the same speed, which is not a property of this contract. Two things are
/// checked instead, and they are the ones that can actually break:
///
/// * both sides name stages, and the stages each side records overlap -- a stage
///   that silently stopped being timed on one side would otherwise go unnoticed;
/// * every recorded value is finite and non-negative.
///
/// The stage *granularity* differs by design: the Rust side splits the core into
/// `mle1`, `sandwich1`, `em`, `correction`, `mle2`, `sandwich2` and `tests`,
/// while the oracle harness times the whole core as one region. That is a
/// difference in instrumentation, not a divergence in the algorithm, so the check
/// is that the oracle's stages are each accounted for by the Rust set rather than
/// that the two lists are equal.
pub fn compare_stage_timings(golden_dir: &std::path::Path, r: &CoreOutput) {
    let path = golden_dir.join("stage_seconds.json");
    let Some(text) = std::fs::read_to_string(&path).ok() else {
        // A golden set generated before the oracle harness recorded timings. Not
        // an error: `make goldens-drift` regenerates, and until then the
        // quantity is simply absent rather than wrong.
        return;
    };
    let v = parse_json(&text);
    let Some(secs) = v.get("seconds").and_then(|s| s.as_obj()) else {
        panic!("{} has no `seconds` object", path.display());
    };
    assert!(
        !secs.is_empty(),
        "{} records no stages at all -- the oracle's instrumentation was removed",
        path.display()
    );
    let rust: std::collections::HashMap<&str, f64> = r.timings.as_pairs().into_iter().collect();
    for (stage, val) in secs {
        let Some(v) = val.as_num() else {
            panic!("stage `{stage}` in {} is not a number", path.display());
        };
        assert!(
            v.is_finite() && v >= 0.0,
            "stage `{stage}` recorded {v}, which is not a plausible duration"
        );
    }
    // Every oracle stage must be accounted for by the Rust set. The mapping is by
    // name where the names agree, and by region where the oracle times something
    // the Rust side breaks down further.
    let covered = |stage: &str| match stage {
        "sanity_check" | "preprocess" | "structural_zeros" | "core" | "sensitivity" => {
            rust.keys().any(|k| {
                matches!(
                    *k,
                    "preprocess"
                        | "pattern_grouping"
                        | "mle1"
                        | "sandwich1"
                        | "em"
                        | "correction"
                        | "mle2"
                        | "sandwich2"
                        | "tests"
                        | "serialisation"
                        | "sensitivity"
                )
            })
        }
        other => panic!("the oracle recorded an unrecognised stage `{other}`"),
    };
    for (stage, _) in secs.iter() {
        assert!(
            covered(stage),
            "the oracle recorded stage `{stage}` but no Rust stage accounts for it"
        );
    }
    // And the Rust side must actually have produced timings, or the comparison
    // above would pass vacuously.
    for (stage, v) in r.timings.as_pairs() {
        assert!(
            v.is_finite() && v >= 0.0,
            "Rust stage `{stage}` recorded {v}, which is not a plausible duration"
        );
    }
}

pub fn compare_core(
    g: &GoldenSet,
    r: &AncombcResult,
    mask: &RankDeficientMask,
    design: &Matrix,
    adj_method: ancombc2_core::config::AdjustMethod,
) -> Option<Divergence> {
    let core: &CoreOutput = &r.core;
    let p = core.fix_eff.len();
    // A golden matrix is stored column-major, matching R, so its flat index maps
    // as (taxon = i % nr, coefficient = i / nr).
    let name_of = |i: usize| {
        core.taxon_names
            .get(i)
            .cloned()
            .unwrap_or_else(|| format!("taxon_{i}"))
    };
    // Two name vectors, because the reference has two taxon sets: the bias set
    // (`taxa_bias`, which is a superset) drives `y1`, `beta_star`, `var1` and
    // `beta_corrected`, and the reported set drives everything else. Indexing a
    // bias-set quantity against the reported names puts the taxon label in the
    // wrong row -- and, worse, makes the rank-deficiency mask index the wrong
    // taxon, since it is keyed by position.
    let bias_names: Vec<String> = core.taxa_bias.iter().map(|&i| name_of(i)).collect();
    let retained: Vec<String> = core.taxa.iter().map(|&i| name_of(i)).collect();
    let coefs = core.fix_eff.clone();
    // A golden matrix is stored column-major, matching R: flat index i is
    // (taxon = i % nr, coefficient = i / nr).
    let rmat = |nr: usize, nc: usize| MatrixContext {
        taxa: retained.iter().take(nr).cloned().collect(),
        coefs: coefs.iter().take(nc).cloned().collect(),
        row_major: false,
    };
    // As `rmat`, but against the bias-set names. The count is an assertion, not
    // a truncation: a bias-set quantity with a different number of rows means
    // the two taxon sets have been conflated.
    let rmat_bias = |nr: usize, nc: usize| MatrixContext {
        taxa: if nr <= bias_names.len() {
            bias_names[..nr].to_vec()
        } else {
            {
                panic!(
                    "bias-set quantity has {nr} rows but the bias set holds {}",
                    bias_names.len()
                )
            }
        },
        coefs: coefs.iter().take(nc).cloned().collect(),
        row_major: false,
    };
    let plain = MatrixContext::default();
    // `p_adjust` is a pipeline-level choice, so it comes from the configuration
    // that produced the golden rather than being assumed.
    //
    // This was hardcoded to Holm, with the comment above claiming "the golden
    // tables record it" -- they do not, and the assumption meant the
    // re-derivation check only ever tested Holm. Every other method's `q` went
    // unverified: `adjust-hommel`, `adjust-hochberg`, `adjust-BH`, `adjust-BY` and
    // `adjust-none` in the fixture matrix all compared `p_adjust(p, "holm")`
    // against a `q` the oracle had computed with a different method, so the cell
    // failed for the wrong reason and the Hochberg/BH conflation in
    // `ancombc2_stats` survived. It is a parameter now, supplied by the caller
    // from the fixture's or cell's own config.
    let method = adj_method;

    // --- Level A: the retained sets and the coefficient names ---
    let got_taxa: Vec<String> = core
        .taxa
        .iter()
        .map(|&i| {
            core.taxon_names
                .get(i)
                .cloned()
                .unwrap_or_else(|| format!("taxon_{i}"))
        })
        .collect();
    if let Some(d) = compare_strings("taxa_retained", &got_taxa, g.strings("taxa_retained")) {
        return Some(d);
    }
    let got_samples: Vec<String> = core
        .samples
        .iter()
        .map(|&i| {
            core.sample_names
                .get(i)
                .cloned()
                .unwrap_or_else(|| format!("sample_{i}"))
        })
        .collect();
    if let Some(d) = compare_strings(
        "samples_retained",
        &got_samples,
        g.strings("samples_retained"),
    ) {
        return Some(d);
    }
    if let Some(d) = compare_strings("fix_eff", &core.fix_eff, g.strings("fix_eff")) {
        return Some(d);
    }

    // --- the design itself, at Level A ---
    //
    // `fix_eff` above checks the design's *column names*, which is not the same
    // thing. A factor re-coded against a different base level produces the same
    // names and a different matrix: the group dummies shift onto other levels,
    // every fitted value stays correct up to the intercept, and the run looks
    // healthy right up until a taxon whose observed levels do not include the new
    // base is fitted. That is exactly what happened on
    // `int-sparsity90-5group`, and the only thing that surfaced it was a
    // 21-of-64-taxon pattern mismatch in `beta_star` -- the design itself was
    // never compared, so the report could not name the cause.
    //
    // Level A, not Level B: `x` is `model.matrix(meta_data, fix_formula)`, built
    // once with no arithmetic on it, so it is either the same matrix or it is
    // not.
    let (xw, _nr, nc) = g.matrix("x");
    let xctx = MatrixContext {
        taxa: g.strings("samples_retained").to_vec(),
        coefs: coefs.iter().take(nc).cloned().collect(),
        row_major: false,
    };
    // `Matrix` is already column-major (`data[j * rows + i]`, like R), so no
    // `to_column_major` here -- unlike the `Vec<f64>` quantities below, which are
    // stored row-major because that is the order their hot loops walk.
    assert_eq!(
        design.data.len(),
        xw.len(),
        "the design has {} elements and the golden's has {}",
        design.data.len(),
        xw.len()
    );
    if let Some(d) = compare_vectors("x", &design.data, xw, TOL_PREPROCESS, &xctx) {
        return Some(d);
    }

    // --- preprocessing: exact, it is only a log and a subtraction ---
    //
    // The core stores `taxa x samples` row-major because that is the order its
    // hot loops walk; R stores column-major. `to_column_major` bridges the two
    // so the comparison is against exactly what R wrote.
    let (y1, nr, nc) = g.matrix("y1");
    if let Some(d) = compare_vectors(
        "y1",
        &to_column_major("core.y1.data", &core.y1.data, nr, nc),
        y1,
        TOL_PREPROCESS,
        &rmat_bias(nr, nc),
    ) {
        return Some(d);
    }

    // --- MLE 1: the coefficients carry the design fit, the variance carries
    //     the sandwich quirk, so they are compared separately ---
    // `theta` is compared *before* `beta_star` because it is upstream of it, and
    // the contract promises a first-diverging-quantity report: `theta` is
    // `colMeans(y - y_crt_hat)` over the set, so a wrong `theta` changes every
    // `beta_star` in the next thing that happens. With `beta_star` first, the
    // `int-sparsity90-5group` cell reported a 166% error on a coefficient whose
    // actual cause was a `theta` 0.33 away -- three quantities upstream of where
    // the report pointed.
    if let Some(d) = compare_vectors("theta", &core.theta, g.vector("theta"), TOL_BETA, &plain) {
        return Some(d);
    }
    let (bs, nr, nc) = g.matrix("beta_star");
    if let Some(d) = compare_vectors(
        "beta_star",
        &to_column_major("core.beta_star", &core.beta_star, nr, nc),
        bs,
        TOL_BETA,
        &rmat_bias(nr, nc),
    ) {
        return Some(d);
    }
    let (v1, nr, nc) = g.matrix("var1");
    // Asserted, like everything else. It used to be reported rather than asserted
    // on the theory that its residuals came from "non-reproducible" rank-deficient
    // taxa -- the same unsound premise as `INDIRECT_QUANTITIES`, and it was wrong
    // for the same reason: `lm`'s pivoting is the identity, so those fits are
    // determined. It also went through `compare_indirect`, which is now just
    // `compare_vectors` with the class threaded through unused; the calls are
    // left as they are so the diff against the previous revision stays readable.
    if let Some(d) = compare_indirect(
        "var1",
        &to_column_major("core.var1", &core.var1, nr, nc),
        v1,
        TOL_VCOV,
        &rmat_bias(nr, nc),
        mask,
        &mask.bias_exempt(),
    ) {
        return Some(d);
    }

    // --- E-M bias ---
    if let Some(d) = compare_indirect(
        "delta_em",
        &core.delta_em,
        g.vector("delta_em"),
        TOL_EM,
        &plain,
        mask,
        &[],
    ) {
        return Some(d);
    }
    if let Some(d) = compare_indirect(
        "delta_wls",
        &core.delta_wls,
        g.vector("delta_wls"),
        TOL_EM,
        &plain,
        mask,
        &[],
    ) {
        return Some(d);
    }
    if let Some(d) = compare_indirect(
        "var_delta",
        &core.var_delta,
        g.vector("var_delta"),
        TOL_EM,
        &plain,
        mask,
        &[],
    ) {
        return Some(d);
    }

    // --- bias correction and sampling fractions ---
    let (bc, nr, nc) = g.matrix("beta_corr_stage1");
    if let Some(d) = compare_indirect(
        "beta_corr_stage1",
        &to_column_major("core.beta_corrected", &core.beta_corrected, nr, nc),
        bc,
        TOL_BETA.with_operand_scale(bs),
        &rmat_bias(nr, nc),
        mask,
        // `beta_corr_stage1` is a restatement of its two inputs, both of which
        // are handled above, so there is no taxon row to mask here.
        &[],
    ) {
        return Some(d);
    }
    // `theta` averages `y1 - X * beta` over the taxa, so its accuracy is set by
    // the operands' rather than by its own (small) magnitude.
    if let Some(d) = compare_indirect(
        "samp_frac",
        &core.samp_frac,
        g.vector("samp_frac"),
        TOL_BETA.with_operand_scale(y1),
        &plain,
        mask,
        &[],
    ) {
        return Some(d);
    }

    // --- MLE 2 ---
    let (y2, nr, nc) = g.matrix("y2");
    if let Some(d) = compare_vectors(
        "y2",
        &to_column_major("core.y2.data", &core.y2.data, nr, nc),
        y2,
        TOL_PREPROCESS,
        &rmat(nr, nc),
    ) {
        return Some(d);
    }
    // `y_bias_crt = y2 - theta_hat`, a difference of two O(1) quantities whose
    // operands are only known to TOL_BETA, so the same cancellation argument
    // applies with `y2` as the scale.
    let (bcrt, nr, nc) = g.matrix("y_bias_crt");
    if let Some(d) = compare_indirect(
        "y_bias_crt",
        &to_column_major("core.y_bias_crt.data", &core.y_bias_crt.data, nr, nc),
        bcrt,
        TOL_BETA.with_operand_scale(y2),
        &rmat(nr, nc),
        mask,
        &[],
    ) {
        return Some(d);
    }
    // The second MLE fits `y2 - theta_hat`, so its coefficients inherit theta's
    // error; the operands set the scale.
    let (beta, nr, nc) = g.matrix("beta");
    if let Some(d) = compare_indirect(
        "beta",
        &to_column_major("core.beta", &core.beta, nr, nc),
        beta,
        TOL_BETA.with_operand_scale(y2),
        &rmat(nr, nc),
        mask,
        &[],
    ) {
        return Some(d);
    }
    let (vh, nr, nc) = g.matrix("var_hat");
    if let Some(d) = compare_indirect(
        "var_hat",
        &to_column_major("core.var_hat", &core.var_hat, nr, nc),
        vh,
        TOL_VCOV,
        &rmat(nr, nc),
        mask,
        &[],
    ) {
        return Some(d);
    }

    // --- regularisation ---
    //
    // `s0` is `quantile(var_hat[, k], s0_perc)` and `var_hat` contains
    // `var_delta` and `2*sqrt(var*var_delta)`, so `s0` inherits the E-M's
    // parameters through both terms, and `var_delta` is one scalar per
    // coefficient -- so nothing averages out over taxa to dilute it. The precision
    // `s0` can be held to is therefore the precision to which the E-M determined
    // its parameters, which is the size of its last step.
    //
    // Two earlier rules were tried here and both were wrong:
    //
    //  * **the iteration count.** A fit that stopped at 61 iterations with a final
    //    epsilon of 9.1e-6 has "converged" by the recorded test and still leaves
    //    `s0` at 1.07e-9 on the `covariates-10-interaction` cell -- just outside
    //    `rtol 1e-9`. Converged is not the same as accurate.
    //  * **the achieved `delta_em`.** Measured, which made it the obvious
    //    candidate, but it did not predict that case: `delta_em` for the term in
    //    question agrees to better than 1e-9 while `s0` does not, so the
    //    disagreement is downstream of the bias estimate and a bound keyed on it
    //    never engaged.
    //
    // This is the same measure `em_tolerance_for` applies to the mixture, so the
    // contract carries one rule for both quantities rather than two. Every
    // comparison prints the bound it applied, so the slack is visible on every run
    // instead of being an absence of evidence.
    let fits = em_fits(&g.dir);
    let s02_tols: Vec<Tolerance> = match &fits {
        Some(f) => (0..core.s02.len())
            .map(|k| em_tolerance_for(f.get(k).copied()))
            .collect(),
        None => core.s02.iter().map(|_| TOL_S0).collect(),
    };
    if let Some(d) = compare_per_term(
        "s02",
        &core.s02,
        g.vector("s02"),
        &s02_tols,
        &plain,
        mask,
        &[],
    ) {
        return Some(d);
    }
    let (vf, nr, nc) = g.matrix("var_final");
    if let Some(d) = compare_indirect(
        "var_final",
        &to_column_major("core.var_final", &core.var_final, nr, nc),
        vf,
        TOL_VCOV,
        &rmat(nr, nc),
        mask,
        &[],
    ) {
        return Some(d);
    }
    let (se, nr, nc) = g.matrix("se");
    if let Some(d) = compare_indirect(
        "se",
        &to_column_major("core.se", &core.se, nr, nc),
        se,
        TOL_SE,
        &rmat(nr, nc),
        mask,
        &[],
    ) {
        return Some(d);
    }
    if let Some(Golden::Vcov { v, p: pv, .. }) = g.quantities.get("vcov") {
        let ctx = MatrixContext {
            taxa: retained.clone(),
            coefs: (0..*pv)
                .map(|i| core.fix_eff.get(i).cloned().unwrap_or_default())
                .collect(),
            row_major: true,
        };
        if let Some(d) = compare_indirect("vcov", &core.vcov, v, TOL_VCOV, &ctx, mask, &[]) {
            return Some(d);
        }
    }

    // --- inference ---
    // W = beta / se: a ratio of two quantities each known to TOL_SE, so the
    // scale is set by them.
    let (w, w_n, w_c) = g.matrix("W");
    if let Some(d) = compare_indirect(
        "W",
        &to_column_major("core.w", &core.w, w_n, w_c),
        w,
        TOL_SE.with_operand_scale(&core.se),
        &rmat(w_n, w_c),
        mask,
        &[],
    ) {
        return Some(d);
    }
    // A p-value is `2 * (1 - F_t(|W|))`, so `|dp/dW| = 2 f_t(|W|) <= 0.8`: an
    // error in W propagates into p at no more than unit gain. The probability
    // tolerance is therefore the Level C floor plus the *observed* W
    // discrepancy — error propagation, not a loosened bound. When W agrees to
    // 1e-15 the p tolerance is the 1e-10 floor alone.
    let dw = max_abs_diff(&to_column_major("core.w", &core.w, w_n, w_c), w);
    let tol_p = Tolerance {
        level: Level::C,
        rtol: 0.0,
        atol: TOL_PROB.atol + dw,
    };
    let (pv, nr, nc) = g.matrix("p");
    if let Some(d) = compare_indirect(
        "p",
        &to_column_major("core.p", &core.p, nr, nc),
        pv,
        tol_p,
        &rmat(nr, nc),
        mask,
        &[],
    ) {
        return Some(d);
    }
    // `q` is checked in two independent steps rather than with one loose
    // tolerance, because folding p's error into q's would multiply it by the
    // adjustment's rank factor (up to n_taxa) and hide a real bug:
    //
    //   1. our `p.adjust` applied to the *golden* p must equal the golden q
    //      exactly (Level A) - that tests the adjustment in isolation;
    //   2. our q, computed from our p, is compared against (1) at the same
    //      propagated tolerance as p - that tests the pipeline.
    let (qv, nr, nc) = g.matrix("q");
    for k in 0..nc {
        let pcol: Vec<f64> = (0..nr).map(|i| pv[i + k * nr]).collect();
        let recomputed = ancombc2_stats::p_adjust_n(&pcol, method, nr as f64);
        for i in 0..nr {
            let want = qv[i + k * nr];
            if recomputed[i] != want && !(recomputed[i].is_nan() && want.is_nan()) {
                let abs = (recomputed[i] - want).abs();
                return Some(Divergence {
                    quantity: format!("q[p.adjust re-derivation, column {k}]"),
                    level: Level::A,
                    index: Some(i),
                    taxon: retained.get(i).cloned(),
                    coefficient: Some(coefs.get(k).cloned().unwrap_or_default()),
                    got: recomputed[i],
                    want,
                    abs_diff: abs,
                    rel_diff: if want != 0.0 { abs / want.abs() } else { abs },
                    message: "p.adjust on the golden p-values does not reproduce the golden q"
                        .into(),
                });
            }
        }
    }
    // ...and that is checked against the oracle directly, not only
    // self-consistently. Comparing our `q` with `p_adjust(our p)` proves the
    // pipeline is internally coherent but says nothing about whether `q` matches
    // the oracle: a `q` that is a deterministic function of a wrong `p` passes
    // that check. Level C is specified as "p and q atol 1e-10", so `q` is
    // compared against the golden `q` at the same tolerance as `p`, and
    // additionally against `p_adjust(our p)` to keep the coherence check that
    // caught a wrong adjustment method in the first place.
    let q_from_ours: Vec<f64> = (0..nr * nc)
        .map(|i| {
            let k = i / nr;
            let pcol: Vec<f64> = (0..nr).map(|r| core.p[r * nc + k]).collect();
            ancombc2_stats::p_adjust_n(&pcol, method, nr as f64)[i % nr]
        })
        .collect();
    if let Some(d) = compare_vectors(
        "q",
        &to_column_major("core.q", &core.q, nr, nc),
        &q_from_ours,
        tol_p,
        &rmat(nr, nc),
    ) {
        return Some(d);
    }
    // `q` needs its own error-propagation term, and it is a much larger one than
    // `p`'s. `q = p_adjust(p)`, and the adjustment is not 1-Lipschitz: Bonferroni
    // multiplies by `n`, and BH/BY divide by a rank `k` while taking a running
    // minimum, so a perturbation of `dp` can move `q` by up to `n * dp`. Using
    // the `p` tolerance for `q` is therefore not a tighter check, it is a
    // *wrong* one, and it fails on real data: on the `sparsity-010` matrix cell
    // `p` agrees to 3.5e-11 and `q` to 3.5e-9 -- exactly a factor of the 100
    // taxa. The bound used here, `n * observed_dp`, is the worst case any of the
    // seven methods can produce, so it is an upper bound rather than a fudge.
    let dp = max_abs_diff(&to_column_major("core.p", &core.p, nr, nc), pv);
    let tol_q = Tolerance {
        level: Level::C,
        rtol: 0.0,
        atol: TOL_PROB.atol + (nr as f64) * dp,
    };
    // `q` was in `INDIRECT_QUANTITIES` and so is not any more; see
    // `INDIRECT_QUANTITIES` for what replaced the premise. Its tolerance is
    // `atol + n_taxa * dp`, a Level C bound that already widens with the number of
    // taxa being adjusted, and it is asserted.
    if let Some(d) = compare_indirect(
        "q[vs oracle]",
        &to_column_major("core.q", &core.q, nr, nc),
        qv,
        tol_q,
        &rmat(nr, nc),
        mask,
        &[],
    ) {
        return Some(d);
    }

    // --- Level A: the significance calls must agree exactly ---
    if let Some(t) = g.tables.get("diff_abn") {
        if let Some(d) = compare_diff_abn(t, &core.diff_abn, &core.fix_eff) {
            return Some(d);
        }
    }
    let _ = p;
    None
}

/// `diff_*` columns of the golden result table, compared against the Rust calls.
///
/// The design column is identified by name, not by position: the table also
/// carries `diff_robust_*` columns once the sensitivity analysis has run, and a
/// running counter would walk off the end of the Rust buffer.
fn compare_diff_abn(t: &Json, got: &[bool], coefs: &[String]) -> Option<Divergence> {
    let Json::Obj(kv) = t else { return None };
    for (k, v) in kv {
        let Some(name) = k.strip_prefix("diff_") else {
            continue;
        };
        if name.is_empty() || name.starts_with("robust_") {
            continue;
        }
        let Some(col) = coefs.iter().position(|c| c == name) else {
            continue;
        };
        let want: Vec<bool> = match v.as_arr() {
            Some(a) => a.iter().map(|x| matches!(x, Json::Bool(true))).collect(),
            None => continue,
        };
        let n_taxa = want.len();
        let mine: Vec<bool> = (0..n_taxa)
            .map(|i| got.get(i * coefs.len() + col).copied().unwrap_or(false))
            .collect();
        if let Some(mut d) = compare_bools(k, &mine, &want) {
            d.taxon = d.index.map(|i| format!("taxon_{i}"));
            d.coefficient = Some(name.to_string());
            d.message = format!(
                "a significance call differs for {name}; the Rust buffer holds {} coefficients",
                coefs.len()
            );
            return Some(d);
        }
    }
    None
}

fn compare_strings(quantity: &str, got: &[String], want: &[String]) -> Option<Divergence> {
    if got == want {
        return None;
    }
    Some(Divergence {
        quantity: quantity.into(),
        level: Level::A,
        index: None,
        taxon: None,
        coefficient: None,
        got: got.len() as f64,
        want: want.len() as f64,
        abs_diff: 0.0,
        rel_diff: 0.0,
        message: format!("got {got:?}, expected {want:?}"),
    })
}

fn compare_bools(quantity: &str, got: &[bool], want: &[bool]) -> Option<Divergence> {
    if got.len() != want.len() {
        return Some(Divergence {
            quantity: quantity.into(),
            level: Level::A,
            index: None,
            taxon: None,
            coefficient: None,
            got: got.len() as f64,
            want: want.len() as f64,
            abs_diff: 0.0,
            rel_diff: 0.0,
            message: "length mismatch".into(),
        });
    }
    for (i, (&g, &w)) in got.iter().zip(want).enumerate() {
        if g != w {
            return Some(Divergence {
                quantity: quantity.into(),
                level: Level::A,
                index: Some(i),
                taxon: None,
                coefficient: None,
                got: g as i32 as f64,
                want: w as i32 as f64,
                abs_diff: 1.0,
                rel_diff: f64::NAN,
                message: "a significance call differs".into(),
            });
        }
    }
    None
}

/// The largest absolute element-wise difference, NaN-tolerant.
fn max_abs_diff(a: &[f64], b: &[f64]) -> f64 {
    a.iter()
        .zip(b)
        .filter(|(x, y)| x.is_finite() && y.is_finite())
        .map(|(x, y)| (x - y).abs())
        .fold(0.0f64, f64::max)
}

/// Re-lay a row-major `nr x nc` buffer into R's column-major order.
///
/// The core stores `taxa x p` row-major because that is the order its hot loops
/// walk (a taxon's samples are contiguous); R and the golden blobs are
/// column-major. Comparing the two directly would transpose the data and every
/// row would appear wrong, so the layout is bridged explicitly here rather than
/// by loosening a tolerance.
fn to_column_major(quantity: &str, row_major: &[f64], nr: usize, nc: usize) -> Vec<f64> {
    assert_eq!(
        row_major.len(),
        nr * nc,
        "quantity `{quantity}`: the golden is {nr} x {nc} = {} elements but the Rust buffer \
         holds {}; the comparison would have read out of bounds",
        nr * nc,
        row_major.len()
    );
    let mut out = vec![0.0; nr * nc];
    for i in 0..nr {
        for j in 0..nc {
            out[j * nr + i] = row_major[i * nc + j];
        }
    }
    out
}
