//! A C ABI for the R wrapper.
//!
//! The boundary is deliberately narrow: one entry point that takes a JSON
//! request and returns a JSON response, both as NUL-terminated UTF-8. JSON is
//! used because R has a first-class JSON reader in every modern installation and
//! the payloads are small — a count matrix goes in as an array of row *indices*
//! into a caller-supplied dense buffer, not as text.
//!
//! Why a C ABI at all, rather than the wrapper shelling out to the CLI: a shell-out
//! pays process start-up and serialises the whole count matrix through a file, and
//! it makes the wrapper untestable without building the binary first. Here the R
//! package links one shared object and the numbers go straight across.
//!
//! # Ownership
//!
//! Every string returned by this module is allocated with libc's `malloc` and must
//! be released with [`ancombc2_rs_free_string`]. Rust's allocator is not R's, so
//! `R_free` would be wrong.
//!
//! # Errors
//!
//! A failure is returned as a JSON object with an `error` key, not as a null
//! pointer and not by aborting. The wrapper turns that into an R condition.

use std::ffi::{c_char, CStr, CString};
use std::panic::{catch_unwind, AssertUnwindSafe};

use ancombc2_core::config::{AdjustMethod, AncombcConfig, CompatMode};
use ancombc2_core::matrix::Matrix;
use ancombc2_core::pipeline::ancombc2_run_named;
use ancombc2_core::preprocess::CountMatrix;

/// The request the wrapper sends.
///
/// ```json
/// {
///   "counts": [12, 0, 3, ...],        // row-major, taxa x samples
///   "n_taxa": 1000, "n_samp": 100,
///   "taxon_names": ["t1", ...],       // optional; positional labels otherwise
///   "sample_names": ["s1", ...],      // optional
///   "design": [1, 0, 1, 0, ...],      // COLUMN-major, samples x p
///
/// The design is **column-major**, because `ancombc2_core::matrix::Matrix` is --
/// it is R's layout, and the design matrix comes from R's `model.matrix`. The
/// count matrix is row-major, because `CountMatrix` is. One FFI boundary with two
/// layouts is a trap, so the mismatch is called out in
/// `docs/numerical_contract.md` and asserted by the parity test
/// `the_ffi_design_layout_is_column_major`.
///   "fix_eff": ["(Intercept)", "grp2", ...],
///   "group": [0, 0, 1, ...] | null,   // level index per sample
///   "p_adj_method": "holm",
///   "pseudo": 0.0, "prv_cut": 0.1, "lib_cut": 0.0, "s0_perc": 0.05,
///   "alpha": 0.05, "global": false, "pairwise": false,
///   "struc_zero": false, "neg_lb": false,
///   "pseudo_sens": false, "conservative": true,
///   "group_name": "grp",
///   "compat": "ancombc2-2.15"
/// }
/// ```
///
/// Every field except `counts`, `n_taxa`, `n_samp` and `design` is optional and
/// falls back to the reference's default, so a minimal request is
/// `{"counts": [...], "n_taxa": n, "n_samp": m, "design": [...]}`.
#[derive(serde::Deserialize)]
#[serde(default)]
struct Request {
    counts: Vec<f64>,
    n_taxa: usize,
    n_samp: usize,
    design: Vec<f64>,
    fix_eff: Vec<String>,
    taxon_names: Vec<String>,
    sample_names: Vec<String>,
    group: Option<Vec<usize>>,
    group_name: Option<String>,
    p_adj_method: String,
    pseudo: f64,
    pseudo_sens: bool,
    conservative: bool,
    prv_cut: f64,
    lib_cut: f64,
    s0_perc: f64,
    struc_zero: bool,
    neg_lb: bool,
    alpha: f64,
    global: bool,
    pairwise: bool,
    mdfdr_fwer: String,
    iter_tol: f64,
    iter_max: usize,
    em_tol: f64,
    em_max: usize,
    compat: String,
}

impl Default for Request {
    fn default() -> Self {
        // The reference's own defaults, from `ancombc2()`'s signature. Duplicated
        // here rather than derived, because the JSON layer cannot import a Rust
        // `Default` for a struct it is deserialising into field by field; the
        // parity test `defaults_match_the_json_defaults` is what keeps the two in
        // step.
        Self {
            counts: Vec::new(),
            n_taxa: 0,
            n_samp: 0,
            design: Vec::new(),
            fix_eff: Vec::new(),
            taxon_names: Vec::new(),
            sample_names: Vec::new(),
            group: None,
            group_name: None,
            p_adj_method: "holm".into(),
            pseudo: 0.0,
            pseudo_sens: false,
            conservative: true,
            prv_cut: 0.10,
            lib_cut: 0.0,
            s0_perc: 0.05,
            struc_zero: false,
            neg_lb: false,
            alpha: 0.05,
            global: false,
            pairwise: false,
            mdfdr_fwer: "holm".into(),
            iter_tol: 0.01,
            iter_max: 20,
            em_tol: 1e-5,
            em_max: 100,
            compat: "ancombc2-2.15".into(),
        }
    }
}

fn error_json(msg: &str) -> String {
    serde_json::json!({ "error": msg }).to_string()
}

fn run(req: Request) -> Result<String, String> {
    if req.n_taxa == 0 || req.n_samp == 0 {
        return Err("n_taxa and n_samp must both be positive".into());
    }
    let want = req
        .n_taxa
        .checked_mul(req.n_samp)
        .ok_or("n_taxa * n_samp overflows")?;
    if req.counts.len() != want {
        return Err(format!(
            "counts has {} elements but n_taxa * n_samp = {want}",
            req.counts.len()
        ));
    }
    let n_samp = req.n_samp;
    if req.design.len() % n_samp != 0 || req.design.is_empty() {
        return Err(format!(
            "design has {} elements, which is not a positive multiple of n_samp = {n_samp}",
            req.design.len()
        ));
    }
    let p = req.design.len() / n_samp;
    let design = Matrix::from_vec(n_samp, p, req.design).map_err(|e| e.to_string())?;

    let fix_eff = if req.fix_eff.is_empty() {
        (0..p).map(|k| format!("V{k}")).collect()
    } else {
        if req.fix_eff.len() != p {
            return Err(format!(
                "fix_eff has {} names but the design has {p} columns",
                req.fix_eff.len()
            ));
        }
        req.fix_eff
    };

    if let Some(g) = &req.group {
        if g.len() != n_samp {
            return Err(format!(
                "group has {} entries but there are {n_samp} samples",
                g.len()
            ));
        }
    }

    let adjust = AdjustMethod::parse(&req.p_adj_method)
        .map_err(|_| format!("unknown p_adj_method {:?}", req.p_adj_method))?;
    let fwer = AdjustMethod::parse(&req.mdfdr_fwer)
        .map_err(|_| format!("unknown mdfdr fwer method {:?}", req.mdfdr_fwer))?;
    let compat = match req.compat.as_str() {
        "ancombc2-2.15" | "ancombc2_15" => CompatMode::Ancombc2_15,
        "strict" | "strict-spec" => CompatMode::StrictSpec,
        other => return Err(format!("unknown compat mode {other:?}")),
    };

    let mut cfg = AncombcConfig {
        fix_eff: fix_eff.clone(),
        p_adj_method: adjust,
        pseudo: req.pseudo,
        pseudo_sens: req.pseudo_sens,
        conservative: req.conservative,
        prv_cut: req.prv_cut,
        lib_cut: req.lib_cut,
        s0_perc: req.s0_perc,
        group: req.group_name.clone(),
        group_labels: None,
        struc_zero: req.struc_zero,
        neg_lb: req.neg_lb,
        alpha: req.alpha,
        global: req.global,
        pairwise: req.pairwise,
        compat,
        ..Default::default()
    };
    cfg.iter_control.tol = req.iter_tol;
    cfg.iter_control.max_iter = req.iter_max;
    cfg.em_control.tol = req.em_tol;
    cfg.em_control.max_iter = req.em_max;
    cfg.mdfdr_control.fwer_ctrl_method = fwer;

    let counts = CountMatrix::new(req.n_taxa, req.n_samp, req.counts).map_err(|e| e.to_string())?;
    let result = ancombc2_run_named(
        &counts,
        &design,
        req.group.as_deref(),
        &cfg,
        &req.taxon_names,
        &req.sample_names,
    )
    .map_err(|e| e.to_string())?;

    let core = &result.core;
    let taxa: Vec<String> = core
        .taxa
        .iter()
        .map(|&i| {
            core.taxon_names
                .get(i)
                .cloned()
                .unwrap_or_else(|| format!("taxon_{i}"))
        })
        .collect();
    let taxa_bias: Vec<String> = core
        .taxa_bias
        .iter()
        .map(|&i| {
            core.taxon_names
                .get(i)
                .cloned()
                .unwrap_or_else(|| format!("taxon_{i}"))
        })
        .collect();

    let mut out = serde_json::json!({
        "fix_eff": fix_eff,
        "taxa": taxa,
        "taxa_bias": taxa_bias,
        "samples": core.samples,
        "beta": core.beta,
        "se": core.se,
        "w": core.w,
        "p": core.p,
        "q": core.q,
        "diff_abn": core.diff_abn,
        "var_hat": core.var_hat,
        "var_final": core.var_final,
        "s02": core.s02,
        "theta": core.theta,
        "samp_frac": core.samp_frac,
        "delta_em": core.delta_em,
        "delta_wls": core.delta_wls,
        "var_delta": core.var_delta,
        "ml_iterations": core.ml_iterations,
        "warnings": result.warnings,
    });
    if let Some(g) = &core.global {
        out["global"] = serde_json::json!({
            "w": g.w, "p": g.p, "q": g.q, "diff_abn": g.diff_abn,
        });
    }
    if let Some(pt) = &core.pairwise {
        out["pairwise"] = serde_json::json!({
            "colnames": pt.colnames,
            "beta": pt.beta, "se": pt.se, "w": pt.w,
            "p": pt.p, "q": pt.q, "diff_abn": pt.diff_abn,
        });
    }
    if let Some(zi) = &core.zero_ind {
        out["zero_ind"] = serde_json::json!({
            "groups": zi.groups, "taxa": taxa_bias, "zero_ind": zi.zero_ind,
        });
    }
    if let (Some(sc), Some(passed)) = (&result.sensitivity, &result.passed_ss) {
        out["sensitivity"] = serde_json::json!({
            "pseudo": sc.pseudo,
            "colnames": sc.colnames,
            "scores": sc.scores,
            "passed_ss": passed,
            "diff_robust": result.diff_robust,
        });
    }
    Ok(out.to_string())
}

/// # Safety
///
/// `request` must be a valid NUL-terminated UTF-8 string. The returned pointer
/// is owned by the caller and must be released with [`ancombc2_rs_free_string`].
#[no_mangle]
pub unsafe extern "C" fn ancombc2_rs_run(request: *const c_char) -> *mut c_char {
    if request.is_null() {
        return to_c_string(&error_json("the request pointer is null"));
    }
    // SAFETY: the caller promises a valid NUL-terminated string; `to_str` is
    // checked, and the panic guard below catches anything else.
    let text = match unsafe { CStr::from_ptr(request) }.to_str() {
        Ok(t) => t,
        Err(e) => return to_c_string(&error_json(&format!("the request is not UTF-8: {e}"))),
    };
    let req: Request = match serde_json::from_str(text) {
        Ok(r) => r,
        Err(e) => return to_c_string(&error_json(&format!("cannot parse the request: {e}"))),
    };
    // A panic across the FFI boundary is undefined behaviour, and an
    // index-out-of-bounds in the numerics would abort the R session rather than
    // raise an R error. Catching it here turns it into a normal condition.
    let out = match catch_unwind(AssertUnwindSafe(|| run(req))) {
        Ok(Ok(s)) => s,
        Ok(Err(e)) => error_json(&e),
        Err(_) => error_json("the Rust core panicked; this is a bug in ancombc2-rs"),
    };
    to_c_string(&out)
}

/// The version of the compatibility target this library implements.
///
/// # Safety
///
/// The returned pointer is owned by the caller.
#[no_mangle]
pub unsafe extern "C" fn ancombc2_rs_version() -> *mut c_char {
    to_c_string(concat!(
        "ancombc2-rs v0.1 equivalent to ANCOMBC 2.15.2 at ",
        "dc4febdf59badb3a8dfe0c767ef2186323c2199a, R 4.5.x, seed 42"
    ))
}

/// Release a string returned by this library.
///
/// # Safety
///
/// `s` must have come from this library and must not have been freed already.
#[no_mangle]
pub unsafe extern "C" fn ancombc2_rs_free_string(s: *mut c_char) {
    if s.is_null() {
        return;
    }
    // SAFETY: `CString::from_raw` takes back the allocation this module made.
    drop(unsafe { CString::from_raw(s) });
}

fn to_c_string(s: &str) -> *mut c_char {
    // A JSON response never contains an interior NUL, so the lossy conversion
    // cannot actually lose anything; the fallback keeps the signature infallible.
    CString::new(s)
        .unwrap_or_else(|_| CString::new(error_json("the response contained a NUL byte")).unwrap())
        .into_raw()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn call(json: &str) -> serde_json::Value {
        let c = CString::new(json).unwrap();
        let ptr = unsafe { ancombc2_rs_run(c.as_ptr()) };
        assert!(!ptr.is_null());
        let text = unsafe { CStr::from_ptr(ptr) }.to_str().unwrap().to_string();
        unsafe { ancombc2_rs_free_string(ptr) };
        serde_json::from_str(&text).expect("the response must be JSON")
    }

    #[test]
    fn reports_an_error_instead_of_aborting() {
        let v = call("{}");
        assert!(v["error"].is_string(), "{v}");
        let v = call("not json");
        assert!(v["error"].is_string(), "{v}");
    }

    #[test]
    fn checks_the_shapes() {
        let v = call(r#"{"counts":[1,2,3],"n_taxa":2,"n_samp":2,"design":[1,1,1,1]}"#);
        assert!(
            v["error"].as_str().unwrap().contains("n_taxa * n_samp"),
            "{v}"
        );
        let v = call(r#"{"counts":[1,2,3,4],"n_taxa":2,"n_samp":2,"design":[1,1,1]}"#);
        assert!(v["error"].as_str().unwrap().contains("design has"), "{v}");
        let v = call(
            r#"{"counts":[1,2,3,4],"n_taxa":2,"n_samp":2,"design":[1,1,1,1],"fix_eff":["a"]}"#,
        );
        assert!(v["error"].as_str().unwrap().contains("fix_eff has"), "{v}");
        let v =
            call(r#"{"counts":[1,2,3,4],"n_taxa":2,"n_samp":2,"design":[1,1,1,1],"group":[0]}"#);
        assert!(v["error"].as_str().unwrap().contains("group has"), "{v}");
    }

    #[test]
    fn rejects_an_unknown_adjust_method() {
        let v = call(
            r#"{"counts":[1,2,3,4],"n_taxa":2,"n_samp":2,"design":[1,1,1,1],
                "p_adj_method":"nonesuch"}"#,
        );
        assert!(v["error"].as_str().unwrap().contains("p_adj_method"), "{v}");
    }

    /// The JSON defaults and the Rust `AncombcConfig` defaults must agree, or a
    /// minimal request from R would run with different settings than the CLI.
    #[test]
    fn defaults_match_the_json_defaults() {
        let d = Request::default();
        let c = AncombcConfig::default();
        assert_eq!(d.p_adj_method, "holm");
        assert_eq!(d.pseudo, c.pseudo);
        assert_eq!(d.prv_cut, c.prv_cut);
        assert_eq!(d.lib_cut, c.lib_cut);
        assert_eq!(d.s0_perc, c.s0_perc);
        assert_eq!(d.alpha, c.alpha);
        assert_eq!(d.conservative, c.conservative);
        assert!(!d.pseudo_sens);
        assert!(!d.global && !d.pairwise && !d.struc_zero && !d.neg_lb);
        assert_eq!(d.iter_tol, c.iter_control.tol);
        assert_eq!(d.iter_max, c.iter_control.max_iter);
        assert_eq!(d.em_tol, c.em_control.tol);
        assert_eq!(d.em_max, c.em_control.max_iter);
        assert_eq!(d.mdfdr_fwer, "holm");
        assert_eq!(d.compat, "ancombc2-2.15");
    }

    /// The design arrives column-major. A two-column design that is *row*-major
    /// transposes into two identical columns for a 2x2-shaped input, which the
    /// core reports as an unidentifiable intercept; this test is what stops that
    /// being mistaken for a modelling failure.
    #[test]
    fn the_ffi_design_layout_is_column_major() {
        // 4 samples x 2 columns, column-major: col0 all ones, col1 = 0,0,1,1.
        let design = vec![1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0];
        let m = Matrix::from_vec(4, 2, design).unwrap();
        assert!(m.get(0, 0) == 1.0 && m.get(0, 1) == 0.0);
        assert!(m.get(2, 0) == 1.0 && m.get(2, 1) == 1.0);
        // The same bytes read as row-major would give two identical columns.
        let wrong = Matrix::from_vec(4, 2, vec![1.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0]).unwrap();
        assert!(wrong.get(0, 0) == wrong.get(0, 1));
    }

    #[test]
    fn version_reports_the_target() {
        let ptr = unsafe { ancombc2_rs_version() };
        let text = unsafe { CStr::from_ptr(ptr) }.to_str().unwrap().to_string();
        unsafe { ancombc2_rs_free_string(ptr) };
        assert!(
            text.contains("dc4febdf59badb3a8dfe0c767ef2186323c2199a"),
            "{text}"
        );
    }
}
