//! Layer 1: golden parity against the pinned ANCOMBC oracle.
//!
//! Each fixture is run through both the R oracle (offline, into
//! `validation/golden/`) and this pipeline, and every quantity in the contract of
//! `docs/numerical_contract.md` is compared. A failure prints the *first*
//! divergent quantity, in pipeline order, with the taxon, the coefficient, and
//! the three numbers.
//!
//! The small fixtures (1, 2, 3) run on every build. Fixture 4 (9,800 x 500) is
//! marked `ignore` and runs nightly; the rationale is in the CI matrix.

mod golden;

use ancombc2_core::ancombc2_run_named;
use golden::{compare_core, fixture_dir, read_golden};

fn run_fixture(id: usize) {
    let dir = fixture_dir(id);
    if !dir.join("counts.tsv").exists() {
        panic!("fixture {id} is missing: {}", dir.display());
    }
    let f = golden::load_fixture(id);
    let g = read_golden(id);
    let r = ancombc2_run_named(
        &f.counts,
        &f.design,
        Some(&f.group),
        &f.cfg,
        &f.taxon_names,
        &f.sample_names,
        &ancombc2_core::F64Reductions,
    )
    .unwrap_or_else(|e| {
        panic!("fixture {id} failed to run: {e}");
    });
    // The rank-deficient class: a taxon whose usable-sample sub-design is
    // exactly singular has a non-unique least-squares solution, so the
    // coordinate the reference reports is not reproducible. See
    // `golden::rank_deficient_taxa` and docs/reference_behavior.md.
    // The two taxon sets are factorised separately, exactly as `.iter_mle` is
    // called twice with two different responses, so each gets its own flags.
    let design = f.design.select_rows(&r.core.samples);
    let mask = golden::RankDeficientMask::for_run(&design, &r.core.y1, &r.core.y2, &f.group_name);
    if let Some(d) = compare_core(&g, &r, &mask, &design, f.cfg.p_adj_method) {
        panic!(
            "fixture {id} diverges from the ANCOMBC 2.15.2 oracle:\n{d}\n({})",
            mask.summary()
        );
    }
    assert!(
        !mask.over_budget(),
        "fixture {id}: {} -- more taxa became rank deficient than the documented \
         class allows; investigate rather than raise the cap",
        mask.summary()
    );
    println!("fx{id:02}: {}", mask.summary());
    println!(
        "fx{id:02}: {} taxa x {} samples, {} coefficients, {} sensitivity runs -- parity OK",
        r.core.taxa.len(),
        r.core.samples.len(),
        r.core.fix_eff.len(),
        r.sensitivity.as_ref().map(|s| s.pseudo.len()).unwrap_or(0),
    );
    // The contract's per-stage timings. The oracle's numbers are recorded but not
    // compared -- wall-clock is not reproducible -- so this checks that both sides
    // name stages, that every oracle stage is accounted for, and that both sets of
    // values are finite and non-negative.
    golden::compare_pattern_assignment(&golden::golden_dir(id), &r.core);
    golden::compare_em_mixture(&golden::golden_dir(id), &r.core);
    golden::compare_convergence_trace(&golden::golden_dir(id), &r.core);
    golden::compare_stage_timings(&golden::golden_dir(id), &r.core);

    // Audit mode prints what the contract actually achieves, so a tolerance is
    // never a guess about the algorithm instead of a measurement of this
    // implementation against this oracle. Off by default: it is diagnostic, and
    // on a four-test binary with a global accumulator the output would interleave.
    if std::env::var_os("ANCOMBC2_GOLDEN_AUDIT").is_some() {
        golden::report_deviations();
    }
}

#[test]
fn fx01_tiny_two_group() {
    run_fixture(1);
}

#[test]
fn fx02_small_two_group() {
    run_fixture(2);
}

#[test]
fn fx03_three_group_with_sensitivity() {
    run_fixture(3);
}

#[test]
#[ignore = "nightly: 9,800 x 500 with the 50-refit sensitivity analysis"]
fn fx04_large_five_group() {
    run_fixture(4);
}
