//! Reading ANCOM-BC2 inputs and writing its result tables.
//!
//! The formats are the ones the R package and its users actually exchange:
//!
//! * a **counts matrix** in TSV or CSV, taxa as rows and samples as columns,
//!   with an ID column in the corner -- what `read.delim` reads into a data frame
//!   and what `ANCOMBC2` receives as `data`;
//! * a **metadata table** in TSV or CSV, samples as rows, with a sample-name
//!   column;
//! * a **model formula** in R's formula syntax, which is *parsed*, not
//!   interpreted as a string. `~ group + x1 + x2:x3` builds the same columns in
//!   the same order as `stats::model.matrix`, because the design matrix has to
//!   match the oracle's exactly for the parity contract to mean anything.
//!
//! Every rejection carries the line and column, because a counts file with a
//! stray comma in an ID is the single most common way this goes wrong.

use std::collections::BTreeSet;
use std::fmt;
use std::path::Path;

use ancombc2_core::config::{AdjustMethod, AncombcConfig, CompatMode};
use ancombc2_core::error::{AncombcError, Result};
use ancombc2_core::matrix::Matrix;
use ancombc2_core::pipeline::{ancombc2_run_named, AncombcResult, CoreOutput, PairwiseTest};
use ancombc2_core::preprocess::CountMatrix;

pub mod formula;
pub mod table;

pub use formula::{Formula, Term};
pub use table::ResultTable;

/// The delimiter a file uses, decided by its extension.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Delimiter {
    Tab,
    Comma,
}

impl Delimiter {
    /// TSV unless the path ends in `.csv` or `.csv.gz`.
    pub fn of_path(path: &Path) -> Self {
        let name = path
            .file_name()
            .and_then(|s| s.to_str())
            .unwrap_or_default()
            .to_ascii_lowercase();
        if name.ends_with(".csv") || name.ends_with(".csv.gz") {
            Delimiter::Comma
        } else {
            Delimiter::Tab
        }
    }

    fn byte(self) -> u8 {
        match self {
            Delimiter::Tab => b'\t',
            Delimiter::Comma => b',',
        }
    }
}

/// An I/O failure, with enough context to fix the input.
#[derive(Debug)]
pub enum IoError {
    Read {
        path: String,
        source: std::io::Error,
    },
    /// A row whose field count does not match the header.
    Ragged {
        path: String,
        line: usize,
        expected: usize,
        got: usize,
    },
    /// A value that is not a number, or is `NA`/`NaN` where a count is required.
    BadValue {
        path: String,
        line: usize,
        column: String,
        value: String,
        what: &'static str,
    },
    /// An empty file, or a header with no fields.
    Empty { path: String },
    /// A duplicate identifier.
    Duplicate {
        path: String,
        kind: &'static str,
        name: String,
    },
    /// A name referenced by the formula that the metadata does not define.
    UnknownVariable {
        name: String,
        available: Vec<String>,
    },
    /// A structural problem in the formula.
    Formula(String),
    /// Sample names in the counts header and the metadata rows do not match.
    SampleMismatch {
        only_in_counts: Vec<String>,
        only_in_metadata: Vec<String>,
    },
    Write {
        path: String,
        source: std::io::Error,
    },
}

impl fmt::Display for IoError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            IoError::Read { path, source } => write!(f, "cannot read {path}: {source}"),
            IoError::Ragged {
                path,
                line,
                expected,
                got,
            } => write!(
                f,
                "{path} line {line}: {got} fields, but the header has {expected}"
            ),
            IoError::BadValue {
                path,
                line,
                column,
                value,
                what,
            } => write!(
                f,
                "{path} line {line}, column `{column}`: {value:?} is not {what}"
            ),
            IoError::Empty { path } => write!(f, "{path} is empty"),
            IoError::Duplicate { path, kind, name } => {
                write!(f, "{path}: duplicate {kind} name {name:?}")
            }
            IoError::UnknownVariable { name, available } => write!(
                f,
                "the formula uses {name:?}, which is not a metadata column; available: {}",
                available.join(", ")
            ),
            IoError::Formula(m) => write!(f, "cannot parse the formula: {m}"),
            IoError::SampleMismatch {
                only_in_counts,
                only_in_metadata,
            } => write!(
                f,
                "sample names differ: only in the counts header [{}], only in the metadata [{}]",
                only_in_counts.join(", "),
                only_in_metadata.join(", ")
            ),
            IoError::Write { path, source } => write!(f, "cannot write {path}: {source}"),
        }
    }
}

impl std::error::Error for IoError {}

impl From<IoError> for AncombcError {
    fn from(e: IoError) -> Self {
        AncombcError::BadInput(e.to_string())
    }
}

/// Split one delimited line, honouring `"` quoting and doubled `""` escapes.
///
/// `read.table`'s default `quote = "\""`, so a quoted field may contain the
/// delimiter and a doubled quote is a literal one. Both matter: a taxon ID with a
/// comma in it is why this exists.
/// [`split_line`] without allocating a `String` per field.
///
/// `split_line` returns owned text, which is what a metadata column wants -- it
/// has to keep the text. A counts matrix does not: every field is parsed to an
/// `f64` and the text is dropped. Routing the counts reader through `split_line`
/// cost one heap allocation *per cell*, which on the 1000 x 10000 benchmark
/// surface was 10 million allocations, most of them 8 to 32 bytes, and dominated
/// the allocation count for the whole run.
///
/// The fast path is the one a counts matrix actually takes: an unquoted line,
/// where each field is a borrow of the line and the split allocates one `Vec`
/// instead of one `String` per field. Quoting needs the `""` unescaping
/// `split_line` does and cannot borrow, so a quoted line is delegated to it
/// unchanged -- quoted counts are rare, and correctness there is worth more than
/// the allocation.
pub fn split_line_cow<'a>(line: &'a str, delim: u8) -> Vec<std::borrow::Cow<'a, str>> {
    if line.as_bytes().contains(&b'"') {
        return split_line(line, delim)
            .into_iter()
            .map(std::borrow::Cow::Owned)
            .collect();
    }
    line.split(delim as char)
        .map(std::borrow::Cow::Borrowed)
        .collect()
}

pub fn split_line(line: &str, delim: u8) -> Vec<String> {
    let bytes = line.as_bytes();
    let mut out = Vec::new();
    let mut cur = String::new();
    let mut in_quotes = false;
    let mut i = 0;
    while i < bytes.len() {
        let b = bytes[i];
        if in_quotes {
            if b == b'"' {
                if i + 1 < bytes.len() && bytes[i + 1] == b'"' {
                    cur.push('"');
                    i += 2;
                    continue;
                }
                in_quotes = false;
            } else {
                cur.push(b as char);
            }
        } else if b == b'"' && cur.is_empty() {
            in_quotes = true;
        } else if b == delim {
            out.push(std::mem::take(&mut cur));
        } else {
            // Non-ASCII identifiers are common in microbiome tables; pushing
            // byte-by-byte through `char` would mangle them, so the byte runs are
            // copied as UTF-8 where possible.
            if b < 0x80 || cur.is_empty() {
                cur.push(b as char);
            } else {
                // fall through to the char branch below for multi-byte sequences
                let ch = line[i..].chars().next().unwrap_or('\u{fffd}');
                cur.push(ch);
                i += ch.len_utf8();
                continue;
            }
        }
        i += 1;
    }
    out.push(cur);
    out
}

/// A table with an ID column and named numeric columns.
#[derive(Debug, Clone)]
struct RawTable {
    ids: Vec<String>,
    header: Vec<String>,
    /// One row per ID, one value per header column, `NaN` for a missing value.
    rows: Vec<Vec<f64>>,
    path: String,
}

fn read_raw(path: &Path) -> std::result::Result<RawTable, IoError> {
    let text = std::fs::read_to_string(path).map_err(|source| IoError::Read {
        path: path.display().to_string(),
        source,
    })?;
    let delim = Delimiter::of_path(path).byte();
    let mut lines = text.lines().filter(|l| !l.trim().is_empty()).peekable();
    let header_line = lines.next().ok_or_else(|| IoError::Empty {
        path: path.display().to_string(),
    })?;
    let header = split_line(header_line, delim);
    if header.len() < 2 {
        return Err(IoError::Empty {
            path: path.display().to_string(),
        });
    }
    // R's `write.table(x, row.names = TRUE)` writes the *column* names as the
    // header and puts the row names in an unnamed first column, so the header has
    // exactly `ncol` fields and every data row has `ncol + 1`. Assuming a header
    // for the ID column would drop the first sample's name and shift every
    // column after it.
    let ncol = header.len();
    let pstr = path.display().to_string();
    let mut ids = Vec::new();
    let mut rows = Vec::new();
    for (n, line) in lines.enumerate() {
        let lineno = n + 2;
        let fields = split_line_cow(line, delim);
        if fields.len() != ncol + 1 {
            return Err(IoError::Ragged {
                path: pstr.clone(),
                line: lineno,
                expected: ncol + 1,
                got: fields.len(),
            });
        }
        // The id is kept, so it has to be owned. That is one allocation per
        // *taxon*, not per cell, and the row name is the API's output.
        ids.push(fields[0].to_string());
        let mut row = Vec::with_capacity(ncol);
        for (j, raw) in fields[1..].iter().enumerate() {
            row.push(parse_count(raw).map_err(|value| IoError::BadValue {
                path: pstr.clone(),
                line: lineno,
                column: header[j].clone(),
                value: value.to_string(),
                what: "a count",
            })?);
        }
        rows.push(row);
    }
    Ok(RawTable {
        ids,
        header,
        rows,
        path: pstr,
    })
}

/// A count, or `NaN` for the missing-value spellings R accepts.
fn parse_count(s: &str) -> std::result::Result<f64, String> {
    let t = s.trim();
    if t.is_empty() || t == "NA" || t == "NaN" {
        return Ok(f64::NAN);
    }
    t.parse::<f64>().map_err(|_| s.to_string())
}

/// A metadata table: samples as rows, one column per covariate.
#[derive(Debug, Clone)]
pub struct Metadata {
    pub sample_names: Vec<String>,
    /// Column names in file order, excluding the sample-name column.
    pub columns: Vec<String>,
    /// `values[column][row]`, as text. Kept as text because a covariate can be a
    /// factor label, and the design matrix needs the level order the *file*
    /// implies, not a numeric reading.
    pub values: Vec<Vec<String>>,
}

impl Metadata {
    pub fn column(&self, name: &str) -> Option<&[String]> {
        self.columns
            .iter()
            .position(|c| c == name)
            .and_then(|i| self.values.get(i))
            .map(|v| v.as_slice())
    }
}

/// Read a metadata table. Values stay textual; see [`Metadata`].
pub fn read_metadata(path: &Path) -> std::result::Result<Metadata, IoError> {
    let raw = read_raw_text(path)?;
    if raw.header.is_empty() {
        return Err(IoError::Empty {
            path: path.display().to_string(),
        });
    }
    Ok(Metadata {
        sample_names: raw.ids,
        columns: raw.header,
        values: raw.cells,
    })
}

/// Like [`read_raw`] but keeping the cells as text, for metadata.
fn read_raw_text(path: &Path) -> std::result::Result<RawText, IoError> {
    let text = std::fs::read_to_string(path).map_err(|source| IoError::Read {
        path: path.display().to_string(),
        source,
    })?;
    let delim = Delimiter::of_path(path).byte();
    let mut lines = text.lines().filter(|l| !l.trim().is_empty()).peekable();
    let header_line = lines.next().ok_or_else(|| IoError::Empty {
        path: path.display().to_string(),
    })?;
    let header = split_line(header_line, delim);
    if header.is_empty() {
        return Err(IoError::Empty {
            path: path.display().to_string(),
        });
    }
    // A one-covariate metadata table is legal, so only an empty header is an error.
    let ncol = header.len();
    let pstr = path.display().to_string();
    let mut ids = Vec::new();
    let mut cells: Vec<Vec<String>> = Vec::with_capacity(ncol);
    for _ in 0..ncol {
        cells.push(Vec::new());
    }
    for (n, line) in lines.enumerate() {
        let lineno = n + 2;
        // Owned text, unlike the counts reader: a metadata column *is* the text,
        // because a covariate can be a factor label and the design matrix needs
        // the level order the file implies rather than a numeric reading.
        let fields = split_line(line, delim);
        if fields.len() != ncol + 1 {
            return Err(IoError::Ragged {
                path: pstr.clone(),
                line: lineno,
                expected: ncol + 1,
                got: fields.len(),
            });
        }
        ids.push(fields[0].clone());
        for j in 0..ncol {
            cells[j].push(fields[j + 1].clone());
        }
    }
    Ok(RawText { ids, header, cells })
}

struct RawText {
    ids: Vec<String>,
    header: Vec<String>,
    cells: Vec<Vec<String>>,
}

/// A counts matrix, taxa as rows, with the sample names in the header.
#[derive(Debug, Clone)]
pub struct Counts {
    pub taxon_names: Vec<String>,
    pub sample_names: Vec<String>,
    pub counts: CountMatrix,
}

impl Counts {
    pub fn n_taxa(&self) -> usize {
        self.taxon_names.len()
    }
    pub fn n_samp(&self) -> usize {
        self.sample_names.len()
    }
}

/// Read a counts matrix, taxa as rows.
pub fn read_counts(path: &Path) -> std::result::Result<Counts, IoError> {
    let raw = read_raw(path)?;
    if raw.ids.is_empty() {
        return Err(IoError::Empty {
            path: path.display().to_string(),
        });
    }
    let mut seen = BTreeSet::new();
    for id in &raw.ids {
        if !seen.insert(id.clone()) {
            return Err(IoError::Duplicate {
                path: raw.path.clone(),
                kind: "taxon",
                name: id.clone(),
            });
        }
    }
    let mut data = Vec::with_capacity(raw.rows.len() * raw.header.len());
    for row in &raw.rows {
        data.extend_from_slice(row);
    }
    let counts =
        CountMatrix::new(raw.ids.len(), raw.header.len(), data).map_err(|e| IoError::BadValue {
            path: raw.path.clone(),
            line: 0,
            column: "<matrix>".into(),
            value: e.to_string(),
            what: "a rectangular count matrix",
        })?;
    Ok(Counts {
        taxon_names: raw.ids,
        sample_names: raw.header,
        counts,
    })
}

/// A parsed design: the model matrix R's `model.matrix` would produce, plus the
/// group levels.
#[derive(Debug, Clone)]
pub struct Design {
    pub matrix: Matrix,
    /// Column names, in matrix order, including `(Intercept)` first.
    pub colnames: Vec<String>,
    /// Group labels per sample, in metadata row order. `None` when no group was
    /// requested.
    pub group: Option<Vec<String>>,
}

/// Build the design matrix for `formula` from `meta`.
///
/// The column order and naming follow `stats::model.matrix`: the intercept, then
/// the terms of the formula left to right; a character or factor column becomes
/// `name` plus one column per level but the first, named `name<level>`, in the
/// order the levels first appear -- which is the order `factor()` gives, i.e.
/// sorted for text and numeric for numbers. Interactions expand to
/// `name1:name2` products.
pub fn build_design(
    meta: &Metadata,
    formula: &Formula,
    group: Option<&str>,
) -> std::result::Result<Design, IoError> {
    for name in formula.variables() {
        let _ = &name;
        if meta.column(&name).is_none() {
            return Err(IoError::UnknownVariable {
                name,
                available: meta.columns.clone(),
            });
        }
    }
    if let Some(g) = group {
        if meta.column(g).is_none() {
            return Err(IoError::UnknownVariable {
                name: g.to_string(),
                available: meta.columns.clone(),
            });
        }
    }
    let n = meta.sample_names.len();
    // The intercept column is all ones. `Matrix::zeros` would leave it all
    // zeros, which is collinear with nothing and yet makes the design rank
    // deficient -- and the failure surfaces as "covariate col0 is not
    // identifiable", which points at the wrong thing entirely.
    let mut matrix = Matrix::zeros(n, 1);
    for i in 0..n {
        matrix.set(i, 0, 1.0);
    }
    let mut colnames = vec!["(Intercept)".to_string()];
    for term in &formula.terms {
        match term {
            Term::Intercept => {}
            Term::Variable(v) => {
                let col = meta.column(v).expect("checked above");
                // The group column is a factor by construction. ANCOM-BC2's `group`
                // argument *is* a grouping variable, and the reference harness
                // coerces it before `model.matrix`, so a metadata column whose
                // labels happen to be numbers still becomes `group2`, `group3`, ...
                // rather than one numeric column.
                let as_group = Some(v.as_str()) == group;
                if as_group {
                    let levels = as_factor_levels(col).unwrap_or_else(|| numeric_levels(col));
                    if levels.len() >= 2 {
                        for lv in &levels[1..] {
                            let name = format!("{v}{lv}");
                            let mut m = Matrix::zeros(n, 1);
                            for (i, ci) in col.iter().enumerate().take(n) {
                                m.set(i, 0, if ci == lv { 1.0 } else { 0.0 });
                            }
                            matrix = append_column(matrix, m);
                            colnames.push(name);
                        }
                        continue;
                    }
                }
                // R's rule, from `model.matrix` and `data_sanity_check`: a column
                // that is *numeric* in the data frame is used numerically, and only
                // a non-numeric one is coerced with `as.factor`. A continuous
                // covariate read as text is therefore numeric, and treating its
                // 100 distinct values as 100 levels would build a 200-column
                // design and lose the intercept to collinearity.
                match as_factor_levels(col) {
                    None => add_numeric_column(&mut matrix, &mut colnames, v.to_string(), col),
                    Some(levels) if levels.len() < 2 => {
                        // A one-level factor: R errors on it, and a constant numeric
                        // column is collinear with the intercept. Let the core
                        // report the unidentifiable covariate.
                        add_numeric_column(&mut matrix, &mut colnames, v.to_string(), col);
                    }
                    Some(levels) => {
                        // `model.matrix` drops the first level of a factor.
                        for lv in &levels[1..] {
                            let name = format!("{v}{lv}");
                            let mut m = Matrix::zeros(n, 1);
                            for (i, ci) in col.iter().enumerate().take(n) {
                                m.set(i, 0, if ci == lv { 1.0 } else { 0.0 });
                            }
                            matrix = append_column(matrix, m);
                            colnames.push(name);
                        }
                    }
                }
            }
            Term::Interaction(a, b) => {
                let (ca, cb) = (
                    meta.column(a).expect("checked above"),
                    meta.column(b).expect("checked above"),
                );
                let name = format!("{a}:{b}");
                let mut m = Matrix::zeros(n, 1);
                for i in 0..n {
                    m.set(
                        i,
                        0,
                        parse_count(&ca[i]).unwrap_or(f64::NAN)
                            * parse_count(&cb[i]).unwrap_or(f64::NAN),
                    );
                }
                matrix = append_column(matrix, m);
                colnames.push(name);
            }
        }
    }
    let group_labels = group.map(|g| meta.column(g).expect("checked above").to_vec());
    // The names go on the matrix as well as in `Design::colnames`, and they have
    // to be in both places rather than one.
    //
    // The pipeline reads `x.colnames` to work out *which columns are the group
    // contrasts* (`group_columns`), and that is what tells the rank-deficient path
    // how `lm` re-levels the factor per taxon. With the names living only on
    // `Design`, the matrix handed to the core had none, `group_columns` returned
    // empty, and every taxon in a rank-deficient group was fitted as if it had
    // observed every level: on `int-sparsity90-5group` the first iteration's
    // epsilon came out 1.88 against the oracle's 0.998, and three samples had no
    // contributor to `theta` at all.
    //
    // `Design::colnames` stays as a field because it is the public accessor, but
    // it is now a copy of the matrix's own rather than the only record.
    Ok(Design {
        matrix: matrix.with_colnames(colnames.clone()),
        colnames,
        group: group_labels,
    })
}

/// The levels of a *non-numeric* column, in the order `factor()` would order
/// them, or `None` when every value parses as a number -- because then R keeps the
/// column numeric and `model.matrix` gives it one column, not one per value.
pub fn as_factor_levels(col: &[String]) -> Option<Vec<String>> {
    // A missing value is a legal level, and R spells it `NA`.
    let numeric = col
        .iter()
        .all(|v| v.trim().is_empty() || v == "NA" || v.parse::<f64>().is_ok());
    if numeric {
        return None;
    }
    let mut levels: Vec<String> = Vec::new();
    let mut seen = BTreeSet::new();
    for v in col {
        if seen.insert(v.clone()) {
            levels.push(v.clone());
        }
    }
    // `factor()` on text sorts with the collation locale; byte order is the
    // closest available and agrees for the ASCII labels a group column has.
    levels.sort();
    Some(levels)
}

fn add_numeric_column(
    matrix: &mut Matrix,
    colnames: &mut Vec<String>,
    name: String,
    col: &[String],
) {
    let n = col.len();
    let mut m = Matrix::zeros(n, 1);
    for (i, ci) in col.iter().enumerate() {
        m.set(i, 0, parse_count(ci).unwrap_or(f64::NAN));
    }
    // `append_column` copies, so the existing intercept and earlier columns
    // survive; replacing the matrix instead would drop them.
    let left = std::mem::replace(matrix, Matrix::zeros(n, 1));
    *matrix = append_column(left, m);
    colnames.push(name);
}

fn append_column(left: Matrix, right: Matrix) -> Matrix {
    let rows = left.rows;
    let p = left.cols + right.cols;
    let mut out = Matrix::zeros(rows, p);
    for i in 0..rows {
        for j in 0..left.cols {
            out.set(i, j, left.get(i, j));
        }
        for j in 0..right.cols {
            out.set(i, left.cols + j, right.get(i, j));
        }
    }
    out
}

/// The levels of a numeric column, ordered as numbers. Used for the group
/// column, which is a factor even when its labels are numeric.
fn numeric_levels(col: &[String]) -> Vec<String> {
    let mut levels: Vec<String> = Vec::new();
    let mut seen = BTreeSet::new();
    for v in col {
        if seen.insert(v.clone()) {
            levels.push(v.clone());
        }
    }
    levels.sort_by(|a, b| {
        a.parse::<f64>()
            .unwrap_or(f64::NAN)
            .partial_cmp(&b.parse::<f64>().unwrap_or(f64::NAN))
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    levels
}

/// Map group labels to level indices, in the order `factor()` would order them.
pub fn group_indices(labels: &[String]) -> Vec<usize> {
    let levels = as_factor_levels(labels).unwrap_or_else(|| numeric_levels(labels));
    labels
        .iter()
        .map(|l| levels.iter().position(|x| x == l).unwrap_or(0))
        .collect()
}

/// Run ANCOM-BC2 on a counts matrix, a design and a configuration.
pub fn run(counts: &Counts, design: &Design, cfg: &AncombcConfig) -> Result<AncombcResult> {
    run_with_names(
        counts,
        design,
        cfg,
        counts.taxon_names.clone(),
        counts.sample_names.clone(),
    )
}

/// [`run`], with the names to report in the result. When `design.group` holds
/// labels for the *subset* of samples, `sample_names` must match it.
pub fn run_with_names(
    counts: &Counts,
    design: &Design,
    cfg: &AncombcConfig,
    taxon_names: Vec<String>,
    sample_names: Vec<String>,
) -> Result<AncombcResult> {
    let group_index = design.group.as_ref().map(|g| group_indices(g));
    ancombc2_run_named(
        &counts.counts,
        &design.matrix,
        group_index.as_deref(),
        cfg,
        &taxon_names,
        &sample_names,
    )
}

/// The primary result table: one row per retained taxon.
///
/// The column *order* is the reference's, because it is part of the contract:
/// `.ancombc2_prep` builds the table with
/// `cbind(taxon, beta_prim, se_prim, W_prim, p_prim, q_prim, diff_prim)`, so all
/// the `lfc_` columns come first, then all the `se_`, and so on -- not the
/// per-coefficient grouping `lfc_, se_, W_, p_, q_, diff_` repeated. When the
/// sensitivity analysis ran, `flag_fun` appends `passed_ss_*` and then
/// `diff_robust_*`.
pub fn primary_table(result: &AncombcResult) -> ResultTable {
    let core = &result.core;
    let n_taxa = core.taxa.len();
    let p = core.fix_eff.len();
    let has_ss = result.passed_ss.is_some();
    let mut t = ResultTable::new("res", "taxon");
    // Header
    for name in &core.fix_eff {
        t.push_num(&format!("lfc_{name}"));
    }
    for name in &core.fix_eff {
        t.push_num(&format!("se_{name}"));
    }
    for name in &core.fix_eff {
        t.push_num(&format!("W_{name}"));
    }
    for name in &core.fix_eff {
        t.push_num(&format!("p_{name}"));
    }
    for name in &core.fix_eff {
        t.push_num(&format!("q_{name}"));
    }
    for name in &core.fix_eff {
        t.push_bool(&format!("diff_{name}"));
    }
    if has_ss {
        for name in &core.fix_eff {
            t.push_bool(&format!("passed_ss_{name}"));
        }
        for name in &core.fix_eff {
            t.push_bool(&format!("diff_robust_{name}"));
        }
    }
    for i in 0..n_taxa {
        t.write_str(&taxon_name(core, i));
        for a in 0..p {
            t.write_num(core.beta[i * p + a]);
        }
        for a in 0..p {
            t.write_num(core.se[i * p + a]);
        }
        for a in 0..p {
            t.write_num(core.w[i * p + a]);
        }
        for a in 0..p {
            t.write_num(core.p[i * p + a]);
        }
        for a in 0..p {
            t.write_num(core.q[i * p + a]);
        }
        for a in 0..p {
            t.write_bool(core.diff_abn[i * p + a]);
        }
        if let Some(passed) = &result.passed_ss {
            for a in 0..p {
                t.write_bool(passed[i * p + a]);
            }
        }
        if let Some(robust) = &result.diff_robust {
            for a in 0..p {
                t.write_bool(robust[i * p + a]);
            }
        }
        t.end_row();
    }
    t
}

/// The global test table, if the global test ran.
pub fn global_table(core: &CoreOutput) -> Option<ResultTable> {
    let g = core.global.as_ref()?;
    let mut t = ResultTable::new("res_global", "taxon");
    for c in ["W", "p_val", "q_val", "diff_abn"] {
        t.push_str(c);
    }
    for i in 0..g.w.len() {
        t.write_str(&taxon_name(core, i));
        t.write_num(g.w[i]);
        t.write_num(g.p[i]);
        t.write_num(g.q[i]);
        t.write_bool(g.diff_abn[i]);
        t.end_row();
    }
    Some(t)
}

/// The pairwise table, if the pairwise test ran.
pub fn pairwise_table(core: &CoreOutput) -> Option<ResultTable> {
    let pt: &PairwiseTest = core.pairwise.as_ref()?;
    let names = &pt.colnames;
    let n_taxa = core.taxa.len();
    let n_col = names.len();
    let mut t = ResultTable::new("res_pair", "taxon");
    for name in names {
        t.push_num(&format!("lfc_{name}"));
        t.push_num(&format!("se_{name}"));
        t.push_num(&format!("W_{name}"));
        t.push_num(&format!("p_{name}"));
        t.push_num(&format!("q_{name}"));
        t.push_bool(&format!("diff_{name}"));
    }
    for i in 0..n_taxa {
        t.write_str(&taxon_name(core, i));
        for c in 0..n_col {
            let k = i * n_col + c;
            t.write_num(pt.beta[k]);
            t.write_num(pt.se[k]);
            t.write_num(pt.w[k]);
            t.write_num(pt.p[k]);
            t.write_num(pt.q[k]);
            t.write_bool(pt.diff_abn[k]);
        }
        t.end_row();
    }
    Some(t)
}

fn taxon_name(core: &CoreOutput, i: usize) -> String {
    core.taxa
        .get(i)
        .and_then(|&t| core.taxon_names.get(t))
        .cloned()
        .unwrap_or_else(|| format!("taxon_{i}"))
}

/// Parse an adjustment-method name the way `p.adjust` spells it.
pub fn parse_adjust(name: &str) -> Result<AdjustMethod> {
    AdjustMethod::parse(name).map_err(|_| {
        AncombcError::BadInput(format!(
            "unknown p-value adjustment method {name:?}; expected one of holm, \
             hochberg, hommel, bonferroni, BH, BY, none"
        ))
    })
}

/// Parse a compatibility mode.
pub fn parse_compat(name: &str) -> Result<CompatMode> {
    match name {
        "ancombc2-2.15" | "ancombc2_15" | "compat" => Ok(CompatMode::Ancombc2_15),
        "strict" | "strict-spec" | "strict_spec" => Ok(CompatMode::StrictSpec),
        other => Err(AncombcError::BadInput(format!(
            "unknown compatibility mode {other:?}; expected ancombc2-2.15 or strict"
        ))),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The design matrix must carry its own column names, not leave them only on
    /// `Design::colnames`.
    ///
    /// The core reads `x.colnames` to work out which columns are the group
    /// contrasts, and that is what tells the rank-deficient fitting path how `lm`
    /// re-levels the factor for each taxon. With the names living only on the
    /// `Design`, `group_columns` came back empty, every rank-deficient taxon was
    /// fitted as if it had observed every level, and `theta` came out with three
    /// samples having no contributor at all -- while `fix_eff`, which is taken
    /// from `cfg.fix_eff`, stayed perfectly correct, so nothing in the output
    /// hinted at it.
    #[test]
    fn the_design_matrix_carries_its_own_colnames() {
        let meta = Metadata {
            sample_names: (0..6).map(|i| format!("s{i}")).collect(),
            columns: vec!["group".into(), "x1".into()],
            values: vec![
                vec![
                    "1".into(),
                    "1".into(),
                    "2".into(),
                    "2".into(),
                    "3".into(),
                    "3".into(),
                ],
                vec![
                    "0.5".into(),
                    "-1.5".into(),
                    "2.0".into(),
                    "0.25".into(),
                    "-0.75".into(),
                    "1.125".into(),
                ],
            ],
        };
        let f = crate::formula::parse("group + x1").expect("formula parses");
        let design = build_design(&meta, &f, Some("group")).expect("design builds");
        assert_eq!(
            design.matrix.colnames, design.colnames,
            "the matrix the core receives has to carry the same names as the \
             `Design` that produced it"
        );
        assert_eq!(
            design.matrix.colnames,
            vec!["(Intercept)", "group2", "group3", "x1"],
            "an integer-valued grouping column is still a factor: it becomes one \
             column per level but the first, not one numeric column"
        );
        assert_eq!(
            ancombc2_core::test_mod::group_columns(&design.matrix.colnames, "group"),
            vec![1usize, 2],
            "and the group contrasts are recoverable from the matrix alone, which \
             is the only place the core can look"
        );
    }

    #[test]
    fn splits_tabs_commas_and_quotes() {
        assert_eq!(split_line("a\tb\tc", b'\t'), vec!["a", "b", "c"]);
        assert_eq!(split_line("a,b,c", b','), vec!["a", "b", "c"]);
        assert_eq!(split_line("a,\"b,c\",d", b','), vec!["a", "b,c", "d"]);
        assert_eq!(split_line("a,\"b\"\"q\",d", b','), vec!["a", "b\"q", "d"]);
        assert_eq!(split_line("a,,c", b','), vec!["a", "", "c"]);
    }

    #[test]
    fn a_numeric_column_is_not_a_factor() {
        // R keeps a numeric column numeric, so this must not become 100 levels.
        let col: Vec<String> = (0..100).map(|i| format!("{}.{}", i, i * 7 % 13)).collect();
        assert!(as_factor_levels(&col).is_none());
    }

    #[test]
    fn factor_levels_are_sorted_and_deduplicated() {
        let col = vec!["b".to_string(), "a".into(), "b".into()];
        assert_eq!(as_factor_levels(&col).unwrap(), vec!["a", "b"]);
    }

    /// The borrowing split must be indistinguishable from the owning one.
    ///
    /// `split_line_cow` is a performance change to the counts reader, and it is
    /// only safe if it returns exactly what `split_line` returns -- including on
    /// a quoted line, where the `""` unescaping has to happen. A test that only
    /// covered the unquoted fast path would let a divergence through on quoted
    /// input, which is precisely the path that falls back.
    #[test]
    fn the_borrowing_split_equals_the_owning_one() {
        for line in [
            "a\tb\tc",
            "",
            "\t",
            "a\t\tb",
            "one",
            "1\t2.5\tNA\tNaN\t-3e2",
            "a\t\"q\"\tb",
            "\"a,b\"\t\"c\"\"d\"\te",
            "\"\"\"\tx",
        ] {
            for delim in *b"\t," {
                let want = split_line(line, delim);
                let got: Vec<String> = split_line_cow(line, delim)
                    .into_iter()
                    .map(|c| c.into_owned())
                    .collect();
                assert_eq!(want, got, "line {line:?} delim {delim}");
            }
        }
    }

    /// A counts file with quoted fields still reads, which is the branch that
    /// makes the fast path safe to add.
    #[test]
    fn a_quoted_counts_field_is_read_through_the_fallback() {
        let dir = std::env::temp_dir().join("ancombc2_io_quoted_counts");
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("q.counts.tsv");
        // The header carries the sample names only; the row names sit in an
        // unnamed first column, as `write.table(x, row.names = TRUE)` writes it.
        std::fs::write(&path, "s1\ts2\nt1\t1\t\"2\"\nt2\t3\t4\n").unwrap();
        let t = read_counts(&path).expect("a quoted count must parse");
        assert_eq!(t.n_taxa(), 2);
        assert_eq!(t.n_samp(), 2);
        assert_eq!(t.counts.get(0, 0), 1.0);
        assert_eq!(
            t.counts.get(0, 1),
            2.0,
            "the quoted field must parse as a count"
        );
        assert_eq!(t.counts.get(1, 1), 4.0);
        std::fs::remove_file(&path).ok();
    }

    #[test]
    fn a_factor_with_numbers_is_still_a_factor_when_one_value_is_not() {
        // `group` is a numeric-looking column of labels; if any value does not
        // parse, R would read the whole column as character.
        let col = vec!["1".to_string(), "2".into(), "control".into()];
        assert_eq!(as_factor_levels(&col).unwrap(), vec!["1", "2", "control"]);
    }

    #[test]
    fn parse_accepts_rs_missing_spellings() {
        assert!(parse_count("NA").unwrap().is_nan());
        assert!(parse_count("NaN").unwrap().is_nan());
        assert!(parse_count("").unwrap().is_nan());
        assert_eq!(parse_count(" 12 ").unwrap(), 12.0);
        assert!(parse_count("x").is_err());
    }
}
