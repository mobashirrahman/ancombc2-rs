//! Time `qr()` in isolation, so the MLE's dominant cost can be measured and then
//! optimised against a number rather than against a 42-minute benchmark.
//!
//!     cargo run --release --example qrbench -- 5000 9 200
//!
//! reports nanoseconds per factorisation and the effective rate in
//! `2 * n * p^2` flops per second -- the Householder count for one QR.

use ancombc2_core::matrix::{qr, Matrix};
use std::time::Instant;

fn main() {
    let mut a = std::env::args().skip(1);
    let n: usize = a.next().map_or(5000, |s| s.parse().unwrap());
    let p: usize = a.next().map_or(9, |s| s.parse().unwrap());
    let reps: usize = a.next().map_or(200, |s| s.parse().unwrap());

    // A fixed, reproducible matrix: deterministic values with a spread of
    // magnitudes so the diagonal is not degenerate.
    let mut seed = 0x243f_6a88_85a3_08d3u64;
    let mut x = Matrix::zeros(n, p);
    for j in 0..p {
        for i in 0..n {
            seed = seed
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            let u = ((seed >> 11) as f64) / ((1u64 << 53) as f64);
            x.set(i, j, u - 0.5);
        }
    }
    x.set(0, 0, 1.0); // the intercept, as the real design has

    // warm up
    let q = qr(&x);
    assert!(q.rank >= 1);

    let t = Instant::now();
    let mut rank_sum = 0usize;
    for _ in 0..reps {
        let q = qr(&x);
        rank_sum += q.rank;
    }
    let per = t.elapsed().as_secs_f64() / reps as f64;
    let flops = 2.0 * n as f64 * p as f64 * p as f64;
    let per_ms = per * 1e3;
    let gflops = flops / per / 1e9;
    let mean_rank = rank_sum / reps;
    println!(
        "qr  n={n:6} p={p:3}  {per_ms:9.3} ms/factorisation   {gflops:6.2} GFLOP/s   (mean rank {mean_rank})"
    );
}
