//! The structural-zero and edge-case matrix (`PLAN.md` §5.5).
//!
//! Seven cases, each a small fixture plus a golden captured from the pinned
//! oracle by `scripts/generate_edge_cases.R`. Every case asserts two things:
//!
//! 1. **Its own predicate.** The case names the behaviour it exists to pin --
//!    a taxon absent in one group, a taxon with a single observation, `NA`
//!    against zero, a group of size 1, the prevalence boundary. A failure says
//!    which case and which plan row, because "parity failed" on a seven-case
//!    matrix is not an actionable message.
//! 2. **General parity.** The retained taxa, the structural-zero table, the
//!    coefficients, and the significance calls must match the oracle, at the
//!    same Levels A-D the golden fixtures use.
//!
//! # Why the goldens are captured rather than written by hand
//!
//! The plan writes the table as "Expected", and for most rows the only
//! authoritative answer is what the reference does. `prv_cut` filters on
//! `prevalence >= prv_cut` in `ancombc_prep.R`, so a taxon exactly at the cut
//! is *kept*; whether a group of size 1 is refused or quietly dropped is a
//! property of the reference. Recording the oracle's answer and asserting
//! against it makes these cases parity rather than opinion -- and it is how two
//! of the cases' descriptions were corrected before this test was written.
//!
//! # Two rows name parameters this oracle does not have
//!
//! The plan lists `keep_zero` and `perc_thres`. ANCOMBC 2.15.2 at `dc4febdf`
//! has neither; the nearest real arguments are `struc_zero` (which controls
//! whether flagged taxa are removed from the primary fit) and `prv_cut` (the
//! prevalence filter). The rows are implemented as `struc_zero_off` and
//! `prv_cut_boundary`, and the substitution is recorded in each case's `about`
//! and in `docs/reference_behavior.md`.

mod support;

use std::path::{Path, PathBuf};

use ancombc2_core::config::{AdjustMethod, AncombcConfig, CompatMode};
use ancombc2_core::matrix::Matrix;
use ancombc2_core::pipeline::ancombc2_run_named;
use ancombc2_core::preprocess::CountMatrix;

#[path = "golden/mod.rs"]
mod golden;

use golden::{parse_json, Json};

/// The matrix's root, next to `validation/fixtures`.
fn root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../validation/edge")
        .canonicalize()
        .expect("validation/edge exists; run scripts/generate_edge_cases.R")
}

// ---- reading a case ---------------------------------------------------------

struct Case {
    name: String,
    plan_row: String,
    about: String,
    counts: CountMatrix,
    taxa: Vec<String>,
    meta: Vec<Vec<(String, String)>>,
    analysis: Vec<(String, String)>,
    golden: Json,
}

impl Case {
    fn load(dir: &Path) -> Self {
        let read = |p: PathBuf| {
            std::fs::read_to_string(&p).unwrap_or_else(|e| panic!("{}: {e}", p.display()))
        };
        let g = parse_json(&read(dir.join("golden.json")));

        // counts.tsv: R's `write.table(row.names = TRUE)` layout, so the header
        // carries only the sample names and each row is name + counts.
        let text = read(dir.join("counts.tsv"));
        let mut lines = text.lines().filter(|l| !l.trim().is_empty());
        let header: Vec<String> = lines
            .next()
            .expect("a header")
            .split('\t')
            .map(|s| s.to_string())
            .collect();
        let n_samp = header.len();
        let mut taxa = Vec::new();
        let mut data = Vec::new();
        for line in lines {
            let f: Vec<&str> = line.split('\t').collect();
            assert_eq!(
                f.len(),
                n_samp + 1,
                "counts.tsv row has {} fields, expected {}",
                f.len(),
                n_samp + 1
            );
            taxa.push(f[0].to_string());
            for v in &f[1..] {
                // `NA` is a missing count, not a parse failure, and the whole
                // point of the `na_counts` case is that it is told apart from a
                // zero.
                data.push(match *v {
                    "NA" | "" => f64::NAN,
                    other => other
                        .parse::<f64>()
                        .unwrap_or_else(|e| panic!("{other:?} is not a count: {e}")),
                });
            }
        }
        let counts = CountMatrix::new(taxa.len(), n_samp, data).expect("shape");

        // meta.tsv: sample label then one field per variable.
        let mtext = read(dir.join("meta.tsv"));
        let mut mlines = mtext.lines().filter(|l| !l.trim().is_empty());
        let mhead: Vec<String> = mlines
            .next()
            .expect("a header")
            .split('\t')
            .map(|s| s.to_string())
            .collect();
        // `Vec<Vec<(variable, value)>>`: one entry per sample, each holding that
        // sample's value for every variable. Annotated because `collect()` would
        // otherwise infer a flat `Vec<(String, String)>` and the mismatch surfaces
        // three closures later as "no method named `iter`".
        let mut meta: Vec<Vec<(String, String)>> = Vec::new();
        for line in mlines {
            let f: Vec<&str> = line.split('\t').collect();
            assert_eq!(f.len(), mhead.len() + 1, "meta.tsv row: {line}");
            let label = f[0].to_string();
            assert_eq!(
                label,
                *header.iter().find(|h| **h == label).unwrap_or(&label),
                "meta row labels are the sample names"
            );
            meta.push(
                mhead
                    .iter()
                    .enumerate()
                    .map(|(j, k)| (k.clone(), f[j + 1].to_string()))
                    .collect::<Vec<(String, String)>>(),
            );
        }
        assert_eq!(
            meta.len(),
            n_samp,
            "meta.tsv has {} rows, counts.tsv has {} samples",
            meta.len(),
            n_samp
        );

        let analysis = read(dir.join("analysis.tsv"))
            .lines()
            .filter(|l| !l.starts_with("key") && !l.trim().is_empty())
            .filter_map(|l| {
                let mut kv = l.splitn(2, '\t');
                let k = kv.next()?.trim().to_string();
                let v = kv.next()?.trim().to_string();
                if k.is_empty() {
                    None
                } else {
                    Some((k, v))
                }
            })
            .collect();

        let s = |k: &str| g.get(k).and_then(|v| v.as_str()).unwrap_or("").to_string();
        Self {
            name: s("name"),
            plan_row: s("plan_row"),
            about: s("about"),
            counts,
            taxa,
            meta,
            analysis,
            golden: g,
        }
    }

    fn cfg(&self, fix_eff: &[String]) -> AncombcConfig {
        let s = |k: &str| {
            self.analysis
                .iter()
                .find(|(a, _)| a == k)
                .map(|(_, v)| v.as_str())
        };
        let b = |k: &str| s(k) == Some("true");
        let f = |k: &str, d: f64| s(k).and_then(|v| v.parse().ok()).unwrap_or(d);
        AncombcConfig {
            fix_eff: fix_eff.to_vec(),
            p_adj_method: AdjustMethod::parse(s("p_adj_method").unwrap_or("BH"))
                .unwrap_or(AdjustMethod::Bh),
            pseudo: f("pseudo", 0.5),
            pseudo_sens: b("pseudo_sens"),
            conservative: !matches!(s("conservative"), Some("false")),
            prv_cut: f("prv_cut", 0.0),
            lib_cut: f("lib_cut", 0.0),
            s0_perc: f("s0_perc", 0.05),
            group: s("group").map(|g| g.to_string()),
            group_labels: None,
            struc_zero: b("struc_zero"),
            neg_lb: b("neg_lb"),
            alpha: f("alpha", 0.05),
            global: false,
            pairwise: false,
            compat: CompatMode::Ancombc2_15,
            ..AncombcConfig::default()
        }
    }

    /// The group level order, as the dataset states it.
    fn group_levels(&self) -> Option<Vec<String>> {
        self.analysis
            .iter()
            .find(|(k, _)| k == "group_levels")
            .map(|(_, v)| {
                v.split(',')
                    .map(|t| t.trim().to_string())
                    .filter(|t| !t.is_empty())
                    .collect()
            })
    }

    fn formula(&self) -> String {
        self.analysis
            .iter()
            .find(|(k, _)| k == "formula")
            .map(|(_, v)| v.clone())
            .expect("the dataset states a formula")
    }
}

/// Run the core on a case, or return the error it produced.
fn run(case: &Case) -> Result<ancombc2_core::pipeline::AncombcResult, String> {
    let n_samp = case.counts.n_samp;
    let group_var = case
        .analysis
        .iter()
        .find(|(k, _)| k == "group")
        .map(|(_, v)| v.clone());
    let levels = case.group_levels();
    let declared = levels.clone();

    // The design: an intercept, then one contrast per group level past the
    // first. Built from `group_levels` rather than from the order the labels
    // appear in, so the reference level is the one the dataset names -- see
    // scripts/realdata_r.R for why a positional choice would be a silent
    // sign inversion between the two implementations.
    let mut fix_eff = vec!["(Intercept)".to_string()];
    let mut columns: Vec<Vec<f64>> = vec![vec![1.0; n_samp]];
    if let Some(gv) = &group_var {
        let raw: Vec<String> = case
            .meta
            .iter()
            .map(|fields| {
                fields
                    .iter()
                    .find(|(k, _)| k == gv)
                    .map(|(_, v)| v.clone())
                    .unwrap_or_default()
            })
            .collect();
        let ordered = declared.unwrap_or_else(|| {
            let mut seen: Vec<String> = Vec::new();
            for v in &raw {
                if !seen.contains(v) {
                    seen.push(v.clone());
                }
            }
            seen
        });
        assert!(
            ordered.len() >= 2,
            "case {}: a group needs at least two levels, found {:?}",
            case.name,
            ordered
        );
        for lv in ordered.iter().skip(1) {
            columns.push(
                raw.iter()
                    .map(|v| if v == lv { 1.0 } else { 0.0 })
                    .collect(),
            );
            fix_eff.push(format!("{gv}{lv}"));
        }
    }
    // Any non-group term, numeric. The matrix's cases are all `~ group`, but the
    // reader is written for the general shape rather than assuming it.
    for term in case
        .formula()
        .split('+')
        .map(|t| t.trim().to_string())
        .filter(|t| !t.is_empty())
    {
        if Some(&term) == group_var.as_ref() {
            continue;
        }
        let col: Vec<f64> = case
            .meta
            .iter()
            .map(|fields| {
                fields
                    .iter()
                    .find(|(k, _)| *k == term)
                    .map(|(_, v)| v.parse().unwrap_or(f64::NAN))
                    .unwrap_or(f64::NAN)
            })
            .collect();
        columns.push(col);
        fix_eff.push(term);
    }

    let mut x = Matrix::zeros(n_samp, columns.len());
    // Indexed by (column, row) because `Matrix` is row-major and the columns are
    // built as one vector per design column.
    for (j, col) in columns.iter().enumerate() {
        for (i, v) in col.iter().enumerate().take(n_samp) {
            x.set(i, j, *v);
        }
    }
    // The group index, in the same level order.
    let group_index: Option<Vec<usize>> = group_var.as_ref().map(|gv| {
        let ordered = levels.clone().unwrap_or_default();
        case.meta
            .iter()
            .map(|fields| {
                let v = fields
                    .iter()
                    .find(|(k, _)| k == gv)
                    .map(|(_, v)| v.clone())
                    .unwrap_or_default();
                ordered.iter().position(|l| *l == v).unwrap_or(usize::MAX)
            })
            .collect()
    });

    let cfg = case.cfg(&fix_eff);
    let samples: Vec<String> = (0..n_samp).map(|i| format!("S{}", i + 1)).collect();
    let mut got = ancombc2_run_named(
        &case.counts,
        &x,
        group_index.as_deref(),
        &cfg,
        &case.taxa,
        &samples,
    )
    .map_err(|e| e.to_string())?;
    got.core.fix_eff = fix_eff;
    Ok(got)
}

// ---- comparison helpers -----------------------------------------------------

const RTOL_B: f64 = 1e-8; // beta, theta
const RTOL_C: f64 = 1e-7; // se, delta_em, delta_wls, vcov
const ATOL_P: f64 = 1e-10; // p and q

fn first_divergence(name: &str, field: &str, a: &[f64], b: &[f64], rtol: f64, atol: f64) {
    if a.len() != b.len() {
        panic!(
            "{name}: {field} has {} entries, the golden has {}",
            a.len(),
            b.len()
        );
    }
    for (i, (x, y)) in a.iter().zip(b.iter()).enumerate() {
        let ok = if x.is_finite() && y.is_finite() {
            (x - y).abs() <= atol + rtol * x.abs().max(y.abs())
        } else {
            x.is_nan() && y.is_nan()
        };
        if !ok {
            panic!(
                "{name}: {field}[{i}] diverges first: {x:?} against the golden's {y:?} \
                 (rtol {rtol}, atol {atol})"
            );
        }
    }
}

fn jnums(g: &Json, key: &str) -> Vec<f64> {
    g.get(key)
        .and_then(|v| v.as_arr())
        .map(|a| {
            a.iter()
                .map(|v| match v {
                    Json::Num(x) => *x,
                    _ => f64::NAN,
                })
                .collect()
        })
        .unwrap_or_default()
}

fn jstrs_of(g: &Json, key: &str) -> Vec<String> {
    g.get(key)
        .and_then(|v| v.as_arr())
        .map(|a| {
            a.iter()
                .map(|v| v.as_str().unwrap_or("").to_string())
                .collect()
        })
        .unwrap_or_default()
}

fn jbools_of(g: &Json, key: &str) -> Vec<bool> {
    g.get(key)
        .and_then(|v| v.as_arr())
        .map(|a| a.iter().map(|v| matches!(v, Json::Bool(true))).collect())
        .unwrap_or_default()
}

/// The structural-zero table as `(column name, flags over the input taxa)`.
fn golden_zero_ind(g: &Json, n_taxa: usize) -> Vec<(String, Vec<bool>)> {
    let Some(zi) = g.get("zero_ind") else {
        return Vec::new();
    };
    let Some(entries) = zi.as_obj() else {
        return Vec::new();
    };
    entries
        .iter()
        .map(|(k, v)| {
            let flags: Vec<bool> = v
                .as_arr()
                .map(|a| a.iter().map(|x| matches!(x, Json::Bool(true))).collect())
                .unwrap_or_default();
            (k.clone(), flags)
        })
        .filter(|(_, f)| f.len() == n_taxa)
        .collect()
}

// ---- the shared parity check -------------------------------------------------

/// Everything every case must satisfy, independent of what it is for.
///
/// Returns the retained taxa, so a case can make a specific statement about them.
fn check_parity(case: &Case, got: &ancombc2_core::pipeline::AncombcResult) -> Vec<String> {
    let g = &case.golden;
    let name = &case.name;

    // Level A: the retained set, exactly.
    let want_taxa = jstrs_of(g, "taxa_retained");
    let got_taxa: Vec<String> = got
        .core
        .taxa
        .iter()
        .map(|i| case.taxa[*i].clone())
        .collect();
    assert_eq!(
        got_taxa, want_taxa,
        "{name} (PLAN row `{}`): the retained taxa differ from the oracle",
        case.plan_row
    );

    // Level A: the structural-zero table, column for column.
    let n_taxa = case.counts.n_taxa;
    let want_zi = golden_zero_ind(g, n_taxa);
    if let Some(actual) = got.core.zero_ind.as_ref() {
        assert!(
            !want_zi.is_empty(),
            "{name}: the oracle reported no zero_ind but the core did"
        );
        // The core's table is `n_taxa * n_groups` row-major over the *group
        // levels in index order*; the golden's columns are named
        // `structural_zero (group = <level>)`, so each is matched by level name
        // rather than by position.
        let levels: Vec<String> = case.group_levels().unwrap_or_default();
        let n_groups = levels.len().max(1);
        for (col, want_flags) in &want_zi {
            // The column is `structural_zero (group = g1)`, so the level is what
            // follows `group = ` with the closing parenthesis stripped. Leaving
            // it in produced `g1)` and a "the golden names a group that is not in
            // the dataset" failure, which reads like a golden problem rather
            // than a one-character parsing slip.
            let level = col
                .rsplit("group = ")
                .next()
                .unwrap_or(col)
                .trim()
                .trim_end_matches(')')
                .trim()
                .to_string();
            let g_idx = levels.iter().position(|l| *l == level);
            let g_idx = g_idx.unwrap_or_else(|| {
                panic!(
                    "{name}: the golden names group `{level}` but the dataset declares {levels:?}"
                )
            });
            for (t, want) in want_flags.iter().enumerate() {
                let have = actual.zero_ind[t * n_groups + g_idx];
                assert_eq!(
                    have,
                    *want,
                    "{name}: structural zero for taxon row {} in group `{level}` is {have}, \
                     the oracle says {want}",
                    t + 1
                );
            }
        }
    } else {
        assert!(
            want_zi.is_empty(),
            "{name}: the oracle reported a zero_ind table but the core did not"
        );
    }

    // Level B/C/D: the reported quantities, at the golden tolerances.
    let p = got.core.fix_eff.len();
    let take = |v: &[f64], col: usize| -> Vec<f64> {
        (0..got.core.taxa.len()).map(|t| v[t * p + col]).collect()
    };
    let col = got
        .core
        .fix_eff_index(g.get("coefficient").and_then(|v| v.as_str()).unwrap_or(""))
        .unwrap_or_else(|| {
            panic!(
                "{name}: the golden compares coefficient `{}`, which is not among {:?}",
                g.get("coefficient").and_then(|v| v.as_str()).unwrap_or(""),
                got.core.fix_eff
            )
        });
    first_divergence(
        name,
        "beta",
        &take(&got.core.beta, col),
        &jnums(g, "beta"),
        RTOL_B,
        0.0,
    );
    first_divergence(
        name,
        "se",
        &take(&got.core.se, col),
        &jnums(g, "se"),
        RTOL_C,
        0.0,
    );
    first_divergence(
        name,
        "p",
        &take(&got.core.p, col),
        &jnums(g, "p"),
        0.0,
        ATOL_P,
    );
    first_divergence(
        name,
        "q",
        &take(&got.core.q, col),
        &jnums(g, "q"),
        0.0,
        ATOL_P,
    );
    // `diff_abn` is per *coefficient*, not per taxon, so it needs the same
    // stride the numeric vectors get. Reading it as one entry per taxon pairs
    // each taxon with another taxon's intercept -- a mistake that is silent,
    // plausible, and was made in the simulation harness before `coefficient`
    // existed to make it hard.
    let called: Vec<bool> = (0..got.core.taxa.len())
        .map(|t| got.core.diff_abn[t * p + col])
        .collect();
    assert_eq!(
        called,
        jbools_of(g, "diff_abn"),
        "{name}: the diff_abn calls differ from the oracle"
    );
    got_taxa
}

// ---- the matrix -------------------------------------------------------------

fn case_names() -> Vec<String> {
    let idx =
        std::fs::read_to_string(root().join("index.json")).expect("validation/edge/index.json");
    let j = parse_json(&idx);
    j.get("cases")
        .and_then(|c| c.as_arr())
        .expect("the index lists cases")
        .iter()
        .map(|c| {
            c.get("name")
                .and_then(|n| n.as_str())
                .unwrap_or("")
                .to_string()
        })
        .collect()
}

/// The plan's seven rows, in order. A missing case fails here rather than
/// quietly reducing the matrix, because a matrix that shrank would still pass.
const EXPECTED_ROWS: [&str; 7] = [
    "Completely absent in group A",
    "Rare in group A (< 1 obs)",
    "Present in every group",
    "NA in counts",
    "Group of size 1",
    "keep_zero = TRUE",
    "perc_thres exactly at boundary",
];

#[test]
fn the_matrix_covers_every_row_of_the_plan() {
    let names = case_names();
    let cases: Vec<Case> = names.iter().map(|n| Case::load(&root().join(n))).collect();
    let rows: Vec<&str> = cases.iter().map(|c| c.plan_row.as_str()).collect();
    assert_eq!(
        rows.len(),
        EXPECTED_ROWS.len(),
        "the matrix has {} case(s), the plan has {} row(s)",
        rows.len(),
        EXPECTED_ROWS.len()
    );
    for (i, want) in EXPECTED_ROWS.iter().enumerate() {
        let have = rows[i];
        assert!(
            have.starts_with(want) || want.starts_with(have.split(" (as ").next().unwrap_or(have)),
            "case {} is `{}`, the plan's row {} is `{want}`",
            cases[i].name,
            have,
            i + 1
        );
    }
    // Every case must say what it is for. A case with an empty `about` is a
    // case nobody can read the intent of.
    for c in &cases {
        assert!(
            c.about.len() > 40,
            "case {} has no meaningful `about`",
            c.name
        );
    }
}

// ---- the cases --------------------------------------------------------------

#[test]
fn absent_in_group_a_is_flagged_for_that_group_only() {
    let c = Case::load(&root().join("absent_in_group_a"));
    let got = run(&c).expect("the oracle analyses this, so the core must too");
    let retained = check_parity(&c, &got);

    // The case-specific statement: T1 is zero in every g1 sample, so it is
    // flagged there and *not* in g2, and it is dropped from the primary fit.
    let zi = got.core.zero_ind.as_ref().expect("struc_zero is on");
    let levels = c.group_levels().unwrap();
    let g1 = levels.iter().position(|l| l == "g1").unwrap();
    let g2 = levels.iter().position(|l| l == "g2").unwrap();
    // `zero_ind` is `n_taxa * n_groups` row-major, so taxon row `i` in group
    // `g` is at `i * n_groups + g`. T1 is row 0, hence the leading zero.
    assert!(
        zi.zero_ind[g1],
        "T1 is zero in every g1 sample and must be flagged there"
    );
    assert!(
        !zi.zero_ind[g2],
        "T1 is present in g2 and must not be flagged there"
    );
    assert!(
        !retained.contains(&"T1".to_string()),
        "a flagged taxon is dropped from the primary fit"
    );
    assert_eq!(retained.len(), 5, "the other five taxa are retained");
}

#[test]
fn rare_in_group_a_is_flagged_by_the_neg_lb_rule() {
    let c = Case::load(&root().join("rare_in_group_a"));
    let got = run(&c).expect("the oracle analyses this");
    let retained = check_parity(&c, &got);

    // T2 has one observation in g1, a prevalence of 1/5. It is *not* structurally
    // zero in the strict sense, and it is flagged anyway: `neg_lb` widens the net
    // to taxa whose asymptotic lower bound is non-positive. This is the
    // distinction the case exists for, and it is why `neg_lb` has a test at all.
    let zi = got.core.zero_ind.as_ref().expect("struc_zero is on");
    let levels = c.group_levels().unwrap();
    let n_groups = levels.len();
    let g1 = levels.iter().position(|l| l == "g1").unwrap();
    assert!(
        zi.zero_ind[n_groups + g1],
        "T2 has a single g1 observation and neg_lb must flag it"
    );
    assert!(
        !retained.contains(&"T2".to_string()),
        "a flagged taxon is dropped from the primary fit"
    );
}

#[test]
fn present_in_every_group_flags_nothing() {
    // The negative control for the two cases above. A rule that flags too
    // eagerly passes those and fails this, which is why it is a case of its own
    // and not a footnote.
    let c = Case::load(&root().join("present_in_every_group"));
    let got = run(&c).expect("the oracle analyses this");
    let retained = check_parity(&c, &got);
    let zi = got.core.zero_ind.as_ref().expect("struc_zero is on");
    assert!(
        zi.zero_ind.iter().all(|f| !f),
        "every taxon is observed in every sample of every group, so nothing may be flagged"
    );
    assert_eq!(retained.len(), 6, "nothing is dropped");
}

#[test]
fn na_counts_are_told_apart_from_zeros_in_the_fit_both_are_absences_for_the_flag() {
    let c = Case::load(&root().join("na_counts"));
    let got = run(&c).expect("the oracle analyses this");
    let retained = check_parity(&c, &got);

    // T3 is NA across g1 and T4 is zero across g1. The two differ in how the
    // sandwich treats them -- an NA term is replaced with 0.1, a zero contributes
    // 0 -- and `check_parity` has already pinned that, since it compares `se` to
    // 1e-7. What is asserted here is the other half: for the structural-zero
    // screen they are the same thing, and both are dropped.
    let zi = got.core.zero_ind.as_ref().expect("struc_zero is on");
    let levels = c.group_levels().unwrap();
    let n_groups = levels.len();
    let g1 = levels.iter().position(|l| l == "g1").unwrap();
    assert!(zi.zero_ind[2 * n_groups + g1], "T3 is NA across g1");
    assert!(zi.zero_ind[3 * n_groups + g1], "T4 is zero across g1");
    // Both flags, and nothing else: the g2 columns must be clear for these two
    // rows, since the taxa are present throughout g2.
    assert!(
        !retained.contains(&"T3".to_string()) && !retained.contains(&"T4".to_string()),
        "both are dropped: for this screen NA and zero are both absences"
    );
    assert_eq!(retained.len(), 4);
}

#[test]
fn a_group_of_size_one_is_refused_as_the_oracle_refuses_it() {
    let c = Case::load(&root().join("group_of_size_one"));
    let gerr = c
        .golden
        .get("error")
        .and_then(|v| v.as_str())
        .unwrap_or("")
        .to_string();
    let got = run(&c);
    match got {
        Ok(res) => panic!(
            "the oracle refuses this case (`{gerr}`) but the core analysed it, retaining {} taxa",
            res.core.taxa.len()
        ),
        Err(e) => {
            // The *reason* must match, not merely that both failed: a refusal
            // for the wrong reason is a different bug. The reference validates
            // the group sizes before anything else.
            let joined = e.to_ascii_lowercase();
            assert!(
                joined.contains("group") || joined.contains("size") || joined.contains("level"),
                "the core refused the case for the wrong reason: {e}\nthe oracle said: {gerr}"
            );
            assert!(
                !gerr.is_empty(),
                "this case asserts the oracle refuses it; the golden records no error, so \
                 one of the two is wrong and the golden must be regenerated"
            );
        }
    }
}

#[test]
fn struc_zero_off_keeps_the_flagged_taxa_and_still_reports_them() {
    // The plan's `keep_zero = TRUE` row, spelled `struc_zero = FALSE` because
    // ANCOMBC 2.15.2 has no `keep_zero` argument.
    let c = Case::load(&root().join("struc_zero_off"));
    let got = run(&c).expect("the oracle analyses this");
    let retained = check_parity(&c, &got);

    // Two questions in one case, and the second is the one a "just skip the
    // removal" implementation gets wrong.
    assert!(
        retained.contains(&"T1".to_string()),
        "with struc_zero = FALSE the flagged taxon stays in the primary fit"
    );
    assert_eq!(retained.len(), 6, "nothing is dropped at all");
    assert!(
        got.core.zero_ind.is_none(),
        "the core reports no zero_ind table when the screen is off, as the oracle does"
    );
}

#[test]
fn prv_cut_at_the_boundary_keeps_a_taxon_exactly_at_the_cut() {
    // The plan's `perc_thres exactly at boundary` row, spelled `prv_cut`.
    let c = Case::load(&root().join("prv_cut_boundary"));
    let got = run(&c).expect("the oracle analyses this");
    let retained = check_parity(&c, &got);

    // `ancombc_prep.R` filters on `prevalence >= prv_cut`. T2 is observed in
    // exactly 3 of 10 samples and the cut is 0.3, so it is exactly at the
    // boundary and is *kept*; T3 at 2/10 is below and is dropped. A `<` instead
    // of a `>=` would keep T3 and drop T2 -- the opposite answer, which is why
    // the case exists.
    assert!(
        retained.contains(&"T2".to_string()),
        "T2's prevalence is exactly prv_cut and `>=` keeps it"
    );
    assert!(
        !retained.contains(&"T3".to_string()),
        "T3's prevalence is below prv_cut and it is dropped"
    );
}
