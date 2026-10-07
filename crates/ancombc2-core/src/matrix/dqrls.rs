//! A transcription of R's `Cdqrls` path: `dqrdc2` + `dqrsl` + `dqrls`.
//!
//! `stats::lm.fit` is not a LAPACK call. It is
//! `.Call(C_Cdqrls, x, y, tol, FALSE)`, which runs R's modified LINPACK
//! `dqrdc2` (Householder with *limited* column pivoting: a column is moved to the
//! end only when its norm has fallen below `tol` times its original norm) and
//! then `dqrsl` to produce `qty`, the coefficients and the residuals. Any other
//! correct least-squares routine agrees with it to rounding; byte-identical output
//! needs the *same* rounding, so the algorithm is transcribed line for line from
//! `R-4.5.3/src/appl/{dqrdc2,dqrsl,dqrls}.f`, including its order of operations.
//!
//! # The BLAS is part of the algorithm
//!
//! LINPACK calls `dnrm2`, `ddot`, `daxpy` and `dscal` in whatever BLAS R is linked
//! against, and those routines' rounding is part of the result. The [`Blas`] trait
//! names the three that round (`dscal`, `dswap` and `dcopy` are exact); [`RefBlas`]
//! is the netlib reference BLAS bundled with R (LAPACK 3.12's `dnrm2`, which is
//! Blue's algorithm, and sequential `ddot`/`daxpy`). Other BLAS implementations
//! are other implementations of this trait.
//!
//! All matrices here are **column-major** with leading dimension `n`, exactly as
//! the Fortran has them.

/// The BLAS level-1 routines whose rounding reaches the result.
pub trait Blas {
    /// `dnrm2(n, x, 1)`.
    fn nrm2(x: &[f64]) -> f64;
    /// `ddot(n, x, 1, y, 1)`.
    fn dot(x: &[f64], y: &[f64]) -> f64;
    /// `daxpy(n, a, x, 1, y, 1)`: `y += a * x`.
    fn axpy(a: f64, x: &[f64], y: &mut [f64]);
}

/// Netlib reference BLAS as bundled with R 4.5 (`src/extra/blas`).
pub struct RefBlas;

// Constants of the reference `dnrm2` (Blue's algorithm, LAPACK >= 3.10), for IEEE
// binary64: `radix = 2`, `minexponent = -1021`, `maxexponent = 1024`, `digits = 53`.

/// `2^e` for a normal exponent, exactly.
fn pow2(e: i64) -> f64 {
    f64::from_bits(((1023 + e) as u64) << 52)
}

impl Blas for RefBlas {
    fn nrm2(x: &[f64]) -> f64 {
        if x.is_empty() {
            return 0.0;
        }
        let (tsml, tbig, ssml, sbig) = (pow2(-511), pow2(486), pow2(537), pow2(-538));
        let (mut notbig, mut asml, mut amed, mut abig) = (true, 0.0f64, 0.0f64, 0.0f64);
        for &v in x {
            let ax = v.abs();
            if ax > tbig {
                let s = ax * sbig;
                abig += s * s;
                notbig = false;
            } else if ax < tsml {
                if notbig {
                    let s = ax * ssml;
                    asml += s * s;
                }
            } else {
                amed += ax * ax;
            }
        }
        let (scl, sumsq);
        if abig > 0.0 {
            if amed > 0.0 || amed > f64::MAX || amed.is_nan() {
                abig += (amed * sbig) * sbig;
            }
            scl = 1.0 / sbig;
            sumsq = abig;
        } else if asml > 0.0 {
            if amed > 0.0 || amed > f64::MAX || amed.is_nan() {
                let amed = amed.sqrt();
                let asml = asml.sqrt() / ssml;
                let (ymin, ymax) = if asml > amed {
                    (amed, asml)
                } else {
                    (asml, amed)
                };
                let r = ymin / ymax;
                scl = 1.0;
                sumsq = ymax * ymax * (1.0 + r * r);
            } else {
                scl = 1.0 / ssml;
                sumsq = asml;
            }
        } else {
            scl = 1.0;
            sumsq = amed;
        }
        scl * sumsq.sqrt()
    }

    fn dot(x: &[f64], y: &[f64]) -> f64 {
        let mut s = 0.0;
        for (a, b) in x.iter().zip(y) {
            s += a * b;
        }
        s
    }

    fn axpy(a: f64, x: &[f64], y: &mut [f64]) {
        if a == 0.0 {
            return;
        }
        for (yi, xi) in y.iter_mut().zip(x) {
            *yi += a * xi;
        }
    }
}

/// OpenBLAS's Haswell/Zen kernels, as used by `libopenblas` on this class of CPU
/// (`kernel/x86_64/{nrm2.S,ddot.c,daxpy.c}` plus the `*_microk_haswell-2.c` files).
///
/// * `dnrm2` is the x87 kernel: squares and sums are formed in 80-bit extended
///   precision (64-bit mantissa) in **four interleaved accumulators** (element `i`
///   of each block of eight goes to accumulator `i % 4`, the tail to the first), the
///   accumulators are combined, the square root is taken *in extended precision*
///   and only then rounded to `double`. That double rounding is the point.
/// * `ddot` over 16 or more elements runs sixteen fused-multiply-add lanes
///   (four `ymm` registers of four) and a fixed reduction tree; the `n % 16` tail is
///   a scalar loop.
/// * `daxpy` is fused (`vfmadd231pd`) in blocks of 16, with a scalar tail.
///
/// Whether the scalar tails are themselves contracted to FMA depends on how the
/// kernel was compiled, so it is a const parameter that the tests pin down against
/// R rather than assume.
pub struct OpenBlasHaswell<const TAIL_FMA: bool>;

#[cfg(target_arch = "x86_64")]
fn nrm2_x87(x: &[f64]) -> f64 {
    use core::arch::asm;
    // Four 80-bit accumulators held in memory (10 bytes used of 16).
    let mut acc = [[0u8; 16]; 4];
    let blocks = x.len() / 8 * 8;
    let add = |acc: &mut [u8; 16], v: f64| unsafe {
        asm!(
            "fld qword ptr [{v}]",
            "fmul st(0), st(0)",
            "fld tbyte ptr [{a}]",
            "faddp st(1), st(0)",
            "fstp tbyte ptr [{a}]",
            v = in(reg) &v,
            a = in(reg) acc.as_mut_ptr(),
            out("st(0)") _, out("st(1)") _,
            options(nostack)
        );
    };
    for (slot, &v) in x[..blocks].iter().enumerate() {
        add(&mut acc[slot % 4], v);
    }
    for &v in &x[blocks..] {
        add(&mut acc[0], v);
    }
    let mut out = 0.0f64;
    unsafe {
        // [A,B,C,D] -> C += A; C += B; D += C; sqrt in extended; round to double.
        asm!(
            "fld tbyte ptr [{d}]",
            "fld tbyte ptr [{c}]",
            "fld tbyte ptr [{b}]",
            "fld tbyte ptr [{a}]",
            "faddp st(2), st(0)",
            "faddp st(1), st(0)",
            "faddp st(1), st(0)",
            "fsqrt",
            "fstp qword ptr [{o}]",
            a = in(reg) acc[0].as_ptr(),
            b = in(reg) acc[1].as_ptr(),
            c = in(reg) acc[2].as_ptr(),
            d = in(reg) acc[3].as_ptr(),
            o = in(reg) &mut out,
            out("st(0)") _, out("st(1)") _, out("st(2)") _, out("st(3)") _,
            options(nostack)
        );
    }
    out
}

#[cfg(target_arch = "x86_64")]
impl<const TAIL_FMA: bool> Blas for OpenBlasHaswell<TAIL_FMA> {
    fn nrm2(x: &[f64]) -> f64 {
        nrm2_x87(x)
    }

    fn dot(x: &[f64], y: &[f64]) -> f64 {
        let n = x.len();
        let n1 = n & !15;
        let mut dot = 0.0f64;
        if n1 > 0 {
            // acc[k][l]: register k (0..4), lane l (0..4).
            let mut acc = [[0.0f64; 4]; 4];
            let mut i = 0;
            while i < n1 {
                for (k, a) in acc.iter_mut().enumerate() {
                    for (l, al) in a.iter_mut().enumerate() {
                        let e = i + 4 * k + l;
                        *al = x[e].mul_add(y[e], *al);
                    }
                }
                i += 16;
            }
            // vextractf128 + vaddpd: lane pair (0,2) and (1,3).
            let mut r = [[0.0f64; 2]; 4];
            for k in 0..4 {
                r[k] = [acc[k][0] + acc[k][2], acc[k][1] + acc[k][3]];
            }
            let s45 = [r[0][0] + r[1][0], r[0][1] + r[1][1]];
            let s67 = [r[2][0] + r[3][0], r[2][1] + r[3][1]];
            let s = [s45[0] + s67[0], s45[1] + s67[1]];
            dot = s[0] + s[1]; // vhaddpd
        }
        for i in n1..n {
            dot = if TAIL_FMA {
                y[i].mul_add(x[i], dot)
            } else {
                dot + y[i] * x[i]
            };
        }
        dot
    }

    fn axpy(a: f64, x: &[f64], y: &mut [f64]) {
        if a == 0.0 {
            return; // interface/axpy.c
        }
        let n = x.len();
        let n1 = n & !15;
        for i in 0..n1 {
            y[i] = a.mul_add(x[i], y[i]);
        }
        for i in n1..n {
            y[i] = if TAIL_FMA {
                a.mul_add(x[i], y[i])
            } else {
                y[i] + a * x[i]
            };
        }
    }
}

/// The outputs of `dqrls` for one response column.
#[derive(Debug, Clone)]
pub struct DqrlsFit {
    /// Coefficients in *pivoted* order, length `p`; slots `rank..p` are `0`
    /// (`dqrls` zeroes them; `lm.fit` later turns them into `NA`).
    pub coef: Vec<f64>,
    pub resid: Vec<f64>,
    pub qty: Vec<f64>,
    pub qraux: Vec<f64>,
    /// 0-based: `pivot[j]` is the original column now in position `j`.
    pub pivot: Vec<usize>,
    pub rank: usize,
}

/// R's `dqrdc2`. `x` is `n x p` column-major and is overwritten by the packed QR.
/// Returns the rank `k`.
pub fn dqrdc2<B: Blas>(
    x: &mut [f64],
    n: usize,
    p: usize,
    tol: f64,
    qraux: &mut [f64],
    jpvt: &mut [usize],
    work: &mut [f64],
) -> usize {
    // work is p x 2, column-major: work1 = work[..p], work2 = work[p..].
    let (w1, w2) = work.split_at_mut(p);
    if n > 0 {
        for j in 0..p {
            qraux[j] = B::nrm2(&x[j * n..j * n + n]);
            w1[j] = qraux[j];
            w2[j] = qraux[j];
            if w2[j] == 0.0 {
                w2[j] = 1.0;
            }
        }
    }
    let lup = n.min(p);
    // Fortran `k` is 1-based "one past the last live column".
    let mut k = p + 1;
    for l in 1..=lup {
        let lz = l - 1; // 0-based column/row of the pivot
        loop {
            if l >= k || qraux[lz] >= w2[lz] * tol {
                break;
            }
            for i in 0..n {
                let t = x[lz * n + i];
                for j in l..p {
                    x[(j - 1) * n + i] = x[j * n + i];
                }
                x[(p - 1) * n + i] = t;
            }
            let i = jpvt[lz];
            let t = qraux[lz];
            let tt = w1[lz];
            let ttt = w2[lz];
            for j in l..p {
                jpvt[j - 1] = jpvt[j];
                qraux[j - 1] = qraux[j];
                w1[j - 1] = w1[j];
                w2[j - 1] = w2[j];
            }
            jpvt[p - 1] = i;
            qraux[p - 1] = t;
            w1[p - 1] = tt;
            w2[p - 1] = ttt;
            k -= 1;
        }
        if l != n {
            let len = n - lz;
            let mut nrmxl = B::nrm2(&x[lz * n + lz..lz * n + n]);
            if nrmxl != 0.0 {
                if x[lz * n + lz] != 0.0 {
                    nrmxl = nrmxl.copysign(x[lz * n + lz]);
                }
                let inv = 1.0 / nrmxl;
                for v in &mut x[lz * n + lz..lz * n + n] {
                    *v *= inv; // dscal
                }
                x[lz * n + lz] += 1.0; // 1 + x == x + 1 exactly
                for j in l..p {
                    let (head, tail) = x.split_at_mut(j * n);
                    let xl = &head[lz * n + lz..lz * n + n];
                    let xj = &mut tail[lz..n];
                    let t = -B::dot(xl, xj) / xl[0];
                    B::axpy(t, xl, xj);
                    if qraux[j] != 0.0 {
                        let r = xj[0].abs() / qraux[j];
                        let mut tt = 1.0 - r * r;
                        tt = tt.max(0.0);
                        let t = tt;
                        if t.abs() >= 1e-6 {
                            qraux[j] *= t.sqrt();
                        } else {
                            qraux[j] = B::nrm2(&xj[1..len]);
                            w1[j] = qraux[j];
                        }
                    }
                }
                qraux[lz] = x[lz * n + lz];
                x[lz * n + lz] = -nrmxl;
            }
        }
    }
    (k - 1).min(n)
}

/// `dqrsl` with `job = 1110` (`qty`, `b` and `rsd`), as `dqrls` calls it.
/// `x` is the packed QR from [`dqrdc2`], restored on return.
#[allow(clippy::too_many_arguments)] // mirrors the Fortran signature
fn dqrsl_1110<B: Blas>(
    x: &mut [f64],
    n: usize,
    k: usize,
    qraux: &[f64],
    y: &[f64],
    qty: &mut [f64],
    b: &mut [f64],
    rsd: &mut [f64],
) {
    let ju = k.min(n - 1);
    if ju == 0 {
        qty[0] = y[0];
        if x[0] != 0.0 {
            b[0] = y[0] / x[0];
        }
        rsd[0] = 0.0;
        return;
    }
    qty[..n].copy_from_slice(&y[..n]);
    for j in 0..ju {
        if qraux[j] != 0.0 {
            let temp = x[j * n + j];
            x[j * n + j] = qraux[j];
            let t = -B::dot(&x[j * n + j..j * n + n], &qty[j..n]) / x[j * n + j];
            B::axpy(t, &x[j * n + j..j * n + n], &mut qty[j..n]);
            x[j * n + j] = temp;
        }
    }
    b[..k].copy_from_slice(&qty[..k]);
    if k < n {
        rsd[k..n].copy_from_slice(&qty[k..n]);
    }
    for r in rsd.iter_mut().take(k) {
        *r = 0.0;
    }
    for jj in 0..k {
        let j = k - 1 - jj;
        if x[j * n + j] == 0.0 {
            break; // info = j; dqrls ignores it
        }
        b[j] /= x[j * n + j];
        if j > 0 {
            let t = -b[j];
            let (bh, _) = b.split_at_mut(j);
            B::axpy(t, &x[j * n..j * n + j], bh);
        }
    }
    for jj in 0..ju {
        let j = ju - 1 - jj;
        if qraux[j] != 0.0 {
            let temp = x[j * n + j];
            x[j * n + j] = qraux[j];
            let t = -B::dot(&x[j * n + j..j * n + n], &rsd[j..n]) / x[j * n + j];
            B::axpy(t, &x[j * n + j..j * n + n], &mut rsd[j..n]);
            x[j * n + j] = temp;
        }
    }
}

/// R's `dqrls` for a single response. `x` (`n x p`, column-major) is consumed as
/// the QR workspace, as R's own copy is.
pub fn dqrls<B: Blas>(mut x: Vec<f64>, n: usize, p: usize, y: &[f64], tol: f64) -> DqrlsFit {
    let mut qraux = vec![0.0; p];
    let mut jpvt: Vec<usize> = (0..p).collect();
    let mut work = vec![0.0; 2 * p];
    let k = dqrdc2::<B>(&mut x, n, p, tol, &mut qraux, &mut jpvt, &mut work);
    let mut coef = vec![0.0; p];
    let mut qty = vec![0.0; n];
    let mut resid = vec![0.0; n];
    if k > 0 {
        dqrsl_1110::<B>(&mut x, n, k, &qraux, y, &mut qty, &mut coef, &mut resid);
    } else {
        resid.copy_from_slice(&y[..n]);
    }
    for c in coef.iter_mut().skip(k) {
        *c = 0.0;
    }
    DqrlsFit {
        coef,
        resid,
        qty,
        qraux,
        pivot: jpvt,
        rank: k,
    }
}
