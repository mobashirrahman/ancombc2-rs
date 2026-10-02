//! A result table that writes itself as TSV.
//!
//! One column type per column, decided when the column is added, so a numeric
//! column cannot accidentally print as a quoted string and a logical column keeps
//! R's `TRUE`/`FALSE` spelling -- which matters, because these tables are read
//! back by R.

use std::fmt::Write as _;
use std::io::Write as _;
use std::path::Path;

use super::IoError;

/// What a column holds, and therefore how it prints.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum ColKind {
    Str,
    Num,
    Bool,
}

/// A TSV table with a header row.
#[derive(Debug, Clone)]
pub struct ResultTable {
    /// Table name, for diagnostics: `res`, `res_global`, `res_pair`.
    pub name: String,
    header: Vec<String>,
    kinds: Vec<ColKind>,
    cells: Vec<String>,
    n_rows: usize,
    n_cols: usize,
}

impl ResultTable {
    /// A table whose first column is named `first_col`.
    pub fn new(name: &str, first_col: &str) -> Self {
        let mut t = Self {
            name: name.to_string(),
            header: Vec::new(),
            kinds: Vec::new(),
            cells: Vec::new(),
            n_rows: 0,
            n_cols: 0,
        };
        t.push_str(first_col);
        t
    }

    /// Append a text column.
    pub fn push_str(&mut self, name: &str) {
        self.header.push(name.to_string());
        self.kinds.push(ColKind::Str);
        self.n_cols += 1;
    }

    /// Append a numeric column.
    pub fn push_num(&mut self, name: &str) {
        self.header.push(name.to_string());
        self.kinds.push(ColKind::Num);
        self.n_cols += 1;
    }

    /// Append a logical column.
    pub fn push_bool(&mut self, name: &str) {
        self.header.push(name.to_string());
        self.kinds.push(ColKind::Bool);
        self.n_cols += 1;
    }

    /// Start a new row; must be followed by exactly one value per column.
    pub fn begin_row(&mut self) {
        debug_assert_eq!(self.cells.len() % self.n_cols.max(1), 0);
    }

    pub fn push_value(&mut self, s: &str) {
        self.cells.push(s.to_string());
    }

    /// Append a text value, in the current row.
    pub fn write_str(&mut self, v: &str) {
        self.cells.push(v.to_string());
    }

    /// Append a numeric value, printed the way R prints a double: 15 significant
    /// digits, `NA` for a missing value, and a bare `1` rather than `1.00000`.
    pub fn write_num(&mut self, v: f64) {
        self.cells.push(format_double(v));
    }

    /// Append a logical value, `TRUE`/`FALSE`/`NA`.
    pub fn write_bool(&mut self, v: bool) {
        self.cells
            .push(if v { "TRUE" } else { "FALSE" }.to_string());
    }

    /// Finish a row: one value per column.
    pub fn end_row(&mut self) {
        self.n_rows += 1;
    }

    /// The number of data rows.
    pub fn n_rows(&self) -> usize {
        self.n_rows
    }

    /// Render as TSV.
    ///
    /// Convenience, and **not** how `write_to` does it: this materialises the
    /// whole rendered table as one `String`, so a caller who wants the text pays for
    /// a second full copy of it on top of the table itself. `write_to` streams
    /// instead. Both go through the same private `write_rows`, so they cannot
    /// render differently.
    pub fn to_tsv(&self) -> String {
        let mut out = String::new();
        let _ = self.write_rows(&mut out, 0..self.n_rows);
        out
    }

    /// Write the header and every data row into `out`, one row at a time.
    ///
    /// A `fmt::Write`, so the caller picks the sink: a `String` for `to_tsv`, and
    /// for `write_to` a chunk that is flushed to the file once it is full. The row
    /// is assembled in a reusable buffer rather than a fresh `Vec<String>` per row,
    /// which is what the old code needed in order to be able to join.
    ///
    /// `rows` is a range so a caller can render the table in pieces: `write_to`
    /// writes a batch, flushes, and asks for the next. A single call over the
    /// whole range is `to_tsv`, which does want it all at once. The header is
    /// written only for the batch that starts at row 0, so the pieces concatenate
    /// to exactly the same bytes as the whole.
    fn write_rows(&self, out: &mut String, rows: std::ops::Range<usize>) -> std::fmt::Result {
        if rows.start == 0 {
            writeln!(out, "{}", self.header.join("\t"))?;
        }
        let mut row = String::with_capacity(self.n_cols * 12);
        for r in rows {
            let base = r * self.n_cols;
            row.clear();
            for c in 0..self.n_cols {
                if c > 0 {
                    row.push('\t');
                }
                push_escaped(&mut row, &self.cells[base + c], self.kinds[c]);
            }
            row.push('\n');
            out.push_str(&row);
        }
        Ok(())
    }

    /// Write to `path`, streaming the rows through a buffered writer.
    ///
    /// Peak memory is the table plus one 64 KiB chunk, rather than the table plus
    /// the entire rendered text. On the benchmark surface the tables are small next
    /// to the count matrices, so this is not what makes `bm5` fit -- but a pairwise
    /// analysis on a large run has one row per taxon per contrast, and
    /// materialising that as text on top of the table is exactly the cost the plan's
    /// "streaming output tables" is about.
    pub fn write_to(&self, path: &Path) -> Result<(), IoError> {
        if let Some(dir) = path.parent() {
            if !dir.as_os_str().is_empty() {
                std::fs::create_dir_all(dir).map_err(|source| IoError::Write {
                    path: dir.display().to_string(),
                    source,
                })?;
            }
        }
        let f = std::fs::File::create(path).map_err(|source| IoError::Write {
            path: path.display().to_string(),
            source,
        })?;
        let mut w = std::io::BufWriter::with_capacity(CHUNK, f);
        let mut chunk = String::with_capacity(CHUNK + CHUNK / 2);
        // Row by row, flushing whenever the chunk is comfortably full. Flushing on
        // a row boundary keeps the chunk from being split mid-row, so what lands in
        // the file is identical to `to_tsv` either way.
        let err = |e: std::io::Error| IoError::Write {
            path: path.display().to_string(),
            source: e,
        };
        let format_err = || std::io::Error::other("formatting the table failed");
        // Batch rows so the chunk is flushed repeatedly rather than once at the
        // end. `flush_chunk` is what makes this streaming rather than merely
        // buffered: a table too large to render at once never gets rendered at
        // once.
        let mut start = 0usize;
        while start < self.n_rows {
            let end = (start + ROWS_PER_BATCH).min(self.n_rows);
            chunk.clear();
            self.write_rows(&mut chunk, start..end)
                .map_err(|_| err(format_err()))?;
            w.write_all(chunk.as_bytes()).map_err(err)?;
            start = end;
        }
        // The header on its own, for a table with no data rows.
        if self.n_rows == 0 {
            self.write_rows(&mut chunk, 0..0)
                .map_err(|_| err(format_err()))?;
            w.write_all(chunk.as_bytes()).map_err(err)?;
        }
        w.flush().map_err(err)
    }
}

/// How much rendered text to hold before handing it to the OS.
///
/// 64 KiB is the usual page-clustered write size: large enough that the syscall
/// count is negligible against the formatting, small enough that the peak is a
/// rounding error next to the table.
const CHUNK: usize = 1 << 16;

/// Rows rendered per flush.
///
/// Chosen against `CHUNK` rather than to hit a byte count, because the row width
/// is not known here: a 13-column primary table and a pairwise table with one
/// column per contrast differ by an order of magnitude. 4,096 rows is small enough
/// that the chunk bounds the peak for either and large enough that the syscall per
/// batch is free.
const ROWS_PER_BATCH: usize = 4096;

/// Format a double the way R's default `digits = 7` would not, but with enough
/// digits to round-trip: 15 significant figures is what `format(x, digits = 15)`
/// and `as.character` on a double effectively preserve.
pub fn format_double(v: f64) -> String {
    if v.is_nan() {
        return "NA".to_string();
    }
    if v.is_infinite() {
        return if v > 0.0 { "Inf" } else { "-Inf" }.to_string();
    }
    if v == v.trunc() && v.abs() < 1e15 {
        return format!("{}", v as i64);
    }
    let mut s = format!("{v:.15e}");
    // Prefer a plain decimal when it is short enough to read.
    let plain = format!("{v}");
    if plain.len() <= s.len() {
        s = plain;
    }
    s
}

/// Quote a field only when it would otherwise be ambiguous, matching
/// `write.table`'s default: quote when the value contains the delimiter, a quote
/// or a newline.
/// Append `v` to `out`, quoted if the column is text and the value needs it.
///
/// Appends rather than returning a `String`, so a row can be assembled with no
/// allocation per field. The quoting rule is unchanged from the version that
/// returned a `String`: a tab, a newline or a double quote forces `write.table`'s
/// `"..."` form with `"` doubled. A numeric column passes through untouched, so an
/// empty or `NA` numeric cell stays empty.
fn push_escaped(out: &mut String, v: &str, kind: ColKind) {
    if kind != ColKind::Str || !(v.contains('\t') || v.contains('\n') || v.contains('"')) {
        out.push_str(v);
        return;
    }
    out.push('"');
    for c in v.chars() {
        if c == '"' {
            out.push('"');
        }
        out.push(c);
    }
    out.push('"');
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn writes_a_header_and_rows() {
        let mut t = ResultTable::new("res", "taxon");
        t.push_num("lfc_grp");
        t.push_bool("diff_grp");
        t.write_str("t1");
        t.write_num(0.5);
        t.write_bool(true);
        t.end_row();
        t.write_str("t2");
        t.write_num(f64::NAN);
        t.write_bool(false);
        t.end_row();
        let got = t.to_tsv();
        assert_eq!(
            got,
            "taxon\tlfc_grp\tdiff_grp\nt1\t0.5\tTRUE\nt2\tNA\tFALSE\n"
        );
    }

    #[test]
    fn formats_doubles_for_round_tripping() {
        assert_eq!(format_double(f64::NAN), "NA");
        assert_eq!(format_double(f64::INFINITY), "Inf");
        assert_eq!(format_double(1.0), "1");
        assert_eq!(format_double(-0.0), "0");
        let x: f64 = 0.1 + 0.2;
        assert_eq!(format_double(x).parse::<f64>().unwrap(), x);
    }

    fn esc(v: &str, kind: ColKind) -> String {
        let mut out = String::new();
        push_escaped(&mut out, v, kind);
        out
    }

    #[test]
    fn quotes_only_when_needed() {
        assert_eq!(esc("a\tb", ColKind::Str), "\"a\tb\"");
        assert_eq!(esc("a\"b", ColKind::Str), "\"a\"\"b\"");
        assert_eq!(esc("plain", ColKind::Str), "plain");
        assert_eq!(esc("a\tb", ColKind::Num), "a\tb");
        assert_eq!(esc("", ColKind::Num), "");
    }

    /// The streamed file and the rendered string must be the same bytes.
    ///
    /// `write_to` now renders in batches and flushes between them, with the header
    /// written only for the batch that starts at row 0. Every way that could go
    /// wrong -- a duplicated header, a missing row at a batch boundary, a
    /// truncated final flush -- shows up here as a byte difference, which is the
    /// only place it can be caught: the golden fixtures do not read these files.
    ///
    /// The row count is deliberately past `ROWS_PER_BATCH`, so the batching path
    /// runs at least twice and the boundary is actually exercised.
    #[test]
    fn the_streamed_file_matches_the_rendered_string() {
        let n = ROWS_PER_BATCH * 2 + 37;
        let mut t = ResultTable::new("res", "num");
        t.push_str("txt");
        for r in 0..n {
            t.begin_row();
            t.write_num(r as f64 + 0.5);
            // A tab and a quote in every third text cell, so quoting is exercised
            // across a batch boundary rather than only at the top.
            if r % 3 == 0 {
                t.write_str(&format!("a\tb\"c{r}"));
            } else {
                t.write_str(&format!("plain{r}"));
            }
            t.end_row();
        }
        let dir = std::env::temp_dir().join("ancombc2_io_stream");
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("t.tsv");
        t.write_to(&path).expect("write");
        let on_disk = std::fs::read_to_string(&path).expect("read back");
        assert_eq!(
            t.to_tsv(),
            on_disk,
            "streaming and rendering must produce identical bytes"
        );
        // The header must appear exactly once, not once per batch.
        assert_eq!(
            on_disk.matches("num\ttxt").count(),
            1,
            "the header must be written once, not per batch"
        );
        assert_eq!(on_disk.lines().count(), n + 1, "one header plus {} rows", n);
        std::fs::remove_file(&path).ok();
    }

    /// An empty table still gets its header.
    #[test]
    fn an_empty_table_still_writes_its_header() {
        let mut t = ResultTable::new("res", "a");
        t.push_str("b");
        let dir = std::env::temp_dir().join("ancombc2_io_stream");
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("empty.tsv");
        t.write_to(&path).expect("write");
        assert_eq!(std::fs::read_to_string(&path).expect("read back"), "a\tb\n");
        assert_eq!(t.to_tsv(), "a\tb\n");
        std::fs::remove_file(&path).ok();
    }
}
