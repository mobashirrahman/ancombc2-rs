//! The `dqrls` transcription against R's own `Cdqrls`, bit for bit.
//!
//! Needs a case file from `scripts/make_dqrls_cases.R`, produced by the *same* R
//! (and BLAS) the comparison is meant to hold against. Skipped when
//! `DQRLS_CASES` is unset.
//!
//!     DQRLS_CASES=cases.bin cargo test --release -p ancombc2-core --test dqrls_vs_r

use ancombc2_core::matrix::{dqrls, Blas, OpenBlasHaswell, RefBlas};

fn f64s(path: &str) -> Vec<f64> {
    std::fs::read(path)
        .expect("read DQRLS_CASES")
        .chunks_exact(8)
        .map(|c| f64::from_le_bytes(c.try_into().unwrap()))
        .collect()
}

fn same(a: &[f64], b: &[f64]) -> bool {
    a.len() == b.len() && a.iter().zip(b).all(|(x, y)| x.to_bits() == y.to_bits())
}

/// Mismatching cases for BLAS `B` against the case file.
fn run<B: Blas>(d: &[f64]) -> (usize, Vec<String>) {
    let (mut at, mut cases, mut bad) = (0usize, 0usize, Vec::<String>::new());
    while at < d.len() {
        let (n, p) = (d[at] as usize, d[at + 1] as usize);
        at += 2;
        let mut take = |len: usize| {
            let s = &d[at..at + len];
            at += len;
            s.to_vec()
        };
        let (x, y) = (take(n * p), take(n));
        let (ec, er, eq, eqr, epv, erk) = (take(p), take(n), take(n), take(p), take(p), take(1));
        let f = dqrls::<B>(x, n, p, &y, 1e-7);
        let piv: Vec<f64> = f.pivot.iter().map(|&j| (j + 1) as f64).collect();
        let mut what = vec![];
        if f.rank as f64 != erk[0] {
            what.push("rank");
        }
        if !same(&piv, &epv) {
            what.push("pivot");
        }
        if !same(&f.qraux, &eqr) {
            what.push("qraux");
        }
        if !same(&f.qty, &eq) {
            what.push("effects");
        }
        if !same(&f.coef, &ec) {
            what.push("coef");
        }
        if !same(&f.resid, &er) {
            what.push("resid");
        }
        if !what.is_empty() {
            bad.push(format!("case {cases} (n={n}, p={p}): {what:?}"));
        }
        cases += 1;
    }
    (cases, bad)
}

#[test]
fn dqrls_matches_r_bit_for_bit() {
    let Ok(path) = std::env::var("DQRLS_CASES") else {
        eprintln!("DQRLS_CASES unset; skipping");
        return;
    };
    let d = f64s(&path);
    // DQRLS_BLAS selects the BLAS the case file was drawn under.
    let which = std::env::var("DQRLS_BLAS").unwrap_or_else(|_| "ref".into());
    let (cases, bad) = match which.as_str() {
        "ref" => run::<RefBlas>(&d),
        "haswell" => run::<OpenBlasHaswell<false>>(&d),
        "haswell-fma" => run::<OpenBlasHaswell<true>>(&d),
        other => panic!("unknown DQRLS_BLAS {other}"),
    };
    assert!(
        bad.is_empty(),
        "{} of {cases} cases differ ({which}):\n{}",
        bad.len(),
        bad[..bad.len().min(12)].join("\n")
    );
    eprintln!("{cases} cases bit-identical ({which})");
}
