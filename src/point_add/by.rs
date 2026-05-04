//! Bernstein–Yang divsteps: classical reference harness and moonshot data.
//!
//! References:
//! - D. J. Bernstein, B.-Y. Yang, "Fast constant-time gcd computation and
//!   modular inversion", IACR ePrint 2019/266, TCHES 2019(3).
//!   https://eprint.iacr.org/2019/266
//!
//! This module is analysis-only. It does not change the quantum circuit.
//! It is here so future sessions can keep the moonshot work self-contained
//! inside `src/point_add/`.
//!
//! ## Scope of the classical work here
//! 1. `divstep2` reference for secp256k1.
//! 2. Empirical survey of actual iteration counts on random secp256k1 inputs.
//! 3. Empirical survey of `jumpdivsteps2` matrix-entry magnitudes, to tighten
//!    the reversible cost model for jumped B-Y.
//!
//! ## Key takeaway so far
//! Plain B-Y (`w = 1`) is still worse than Kaliski on raw iteration count.
//! I initially believed jumped B-Y might be re-opened if the empirical
//! transition-matrix entries were much smaller than the worst-case `2^w`
//! bound. After correcting a bug in the matrix-survey code, the updated
//! survey shows the opposite: the low-word jump matrices frequently hit the
//! full `2^w` growth. So the original pessimistic reversible cost model was
//! basically right.

use std::time::Instant;

use alloy_primitives::{U256, U512};
use sha3::digest::{ExtendableOutput, Update, XofReader};

use super::test_timeout::{check_deadline, two_min_deadline};

/// secp256k1 prime: p = 2^256 − 2^32 − 977.
pub const SECP256K1_P: U256 = U256::from_limbs([
    0xFFFFFFFEFFFFFC2F,
    0xFFFFFFFFFFFFFFFF,
    0xFFFFFFFFFFFFFFFF,
    0xFFFFFFFFFFFFFFFF,
]);

/// Theoretical safegcd iteration bound (Bernstein–Yang 2019/266,
/// Theorem 11.2 linearized bound used in the paper's constant-time recip2):
///
///     N_bound(n) = ceil((49 n + 57) / 17)
///
/// For n = 256, this is 742.
pub fn safegcd_iters(n_bits: usize) -> usize {
    (49 * n_bits + 57 + 16) / 17
}

// ─────────────────────────────────────────────────────────────────────────
// Signed integer helper (257-bit via sign + U256 magnitude)
// ─────────────────────────────────────────────────────────────────────────

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct SInt {
    pub neg: bool,
    pub mag: U256,
}

impl SInt {
    pub fn zero() -> Self {
        Self {
            neg: false,
            mag: U256::ZERO,
        }
    }

    pub fn from_u(x: U256) -> Self {
        Self { neg: false, mag: x }
    }

    pub fn negate(self) -> Self {
        if self.mag.is_zero() {
            self
        } else {
            Self {
                neg: !self.neg,
                mag: self.mag,
            }
        }
    }

    pub fn bit0(&self) -> bool {
        // Parity is the same for ±x.
        self.mag.bit(0)
    }

    pub fn is_zero(&self) -> bool {
        self.mag.is_zero()
    }

    pub fn is_one_pos(&self) -> bool {
        !self.neg && self.mag == U256::from(1)
    }

    pub fn is_one_neg(&self) -> bool {
        self.neg && self.mag == U256::from(1)
    }

    pub fn add(a: Self, b: Self) -> Self {
        match (a.neg, b.neg) {
            (false, false) => Self {
                neg: false,
                mag: a.mag.wrapping_add(b.mag),
            },
            (true, true) => Self {
                neg: true,
                mag: a.mag.wrapping_add(b.mag),
            },
            (false, true) => sub_mag(a.mag, b.mag),
            (true, false) => sub_mag(b.mag, a.mag),
        }
    }

    pub fn sub(a: Self, b: Self) -> Self {
        Self::add(a, b.negate())
    }

    pub fn shr1_even(self) -> Self {
        debug_assert!(!self.bit0(), "shr1_even on odd integer");
        Self {
            neg: self.neg,
            mag: self.mag >> 1,
        }
    }
}

fn sub_mag(a: U256, b: U256) -> SInt {
    if a >= b {
        SInt {
            neg: false,
            mag: a.wrapping_sub(b),
        }
    } else {
        SInt {
            neg: true,
            mag: b.wrapping_sub(a),
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────
// Classical modular arithmetic for coefficient tracking
// ─────────────────────────────────────────────────────────────────────────

fn addm(a: U256, b: U256, p: U256) -> U256 {
    a.add_mod(b, p)
}

fn subm(a: U256, b: U256, p: U256) -> U256 {
    let (r, borrow) = a.overflowing_sub(b);
    if borrow {
        r.wrapping_add(p)
    } else {
        r
    }
}

fn negm(a: U256, p: U256) -> U256 {
    if a.is_zero() {
        a
    } else {
        p.wrapping_sub(a)
    }
}

fn mulm(a: U256, b: U256, p: U256) -> U256 {
    a.mul_mod(b, p)
}

// ─────────────────────────────────────────────────────────────────────────
// divstep2 classical reference
// ─────────────────────────────────────────────────────────────────────────

#[derive(Clone, Copy, Debug)]
pub struct Coeffs {
    pub uu: U256,
    pub vv: U256,
    pub qq: U256,
    pub rr: U256,
}

impl Coeffs {
    pub fn initial() -> Self {
        Self {
            uu: U256::from(1),
            vv: U256::ZERO,
            qq: U256::ZERO,
            rr: U256::from(1),
        }
    }
}

#[derive(Debug, Clone)]
pub struct DivstepsRun {
    pub converged: bool,
    pub iters_done: usize,
    pub max_abs_delta: i64,
    pub final_f: SInt,
    pub final_g: SInt,
    pub final_coeffs: Coeffs,
}

/// Run one-step-at-a-time `divstep2` until convergence or until max_iters.
///
/// This follows the integer `divsteps2` of BY 2019/266 Figure 10.1,
/// specialized to modular-inverse tracking over an odd prime modulus p.
pub fn run_divsteps(g0: U256, p: U256, max_iters: usize) -> DivstepsRun {
    assert!(p.bit(0), "p must be odd");
    assert!(g0 < p && !g0.is_zero(), "g0 must lie in [1, p)");

    let mut delta: i64 = 1;
    let mut f = SInt::from_u(p);
    let mut g = SInt::from_u(g0);
    let mut coeffs = Coeffs::initial();
    let mut max_abs_delta = 1i64;
    let mut converged_iter = None;

    for i in 0..max_iters {
        if g.is_zero() {
            converged_iter = Some(i);
            break;
        }

        let g_odd = g.bit0();
        if delta > 0 && g_odd {
            // Case A:
            //   (δ, f, g) ← (1 − δ, g, (g − f) / 2)
            //   (U,V,Q,R) ← (2Q, 2R, Q−U, R−V)
            let nf = g;
            let ng = SInt::sub(g, f).shr1_even();
            let nu = addm(coeffs.qq, coeffs.qq, p);
            let nv = addm(coeffs.rr, coeffs.rr, p);
            let nq = subm(coeffs.qq, coeffs.uu, p);
            let nr = subm(coeffs.rr, coeffs.vv, p);
            delta = 1 - delta;
            f = nf;
            g = ng;
            coeffs = Coeffs {
                uu: nu,
                vv: nv,
                qq: nq,
                rr: nr,
            };
        } else if g_odd {
            // Case B:
            //   (δ, f, g) ← (1 + δ, f, (g + f) / 2)
            //   (U,V,Q,R) ← (2U, 2V, Q+U, R+V)
            let ng = SInt::add(g, f).shr1_even();
            let nu = addm(coeffs.uu, coeffs.uu, p);
            let nv = addm(coeffs.vv, coeffs.vv, p);
            let nq = addm(coeffs.qq, coeffs.uu, p);
            let nr = addm(coeffs.rr, coeffs.vv, p);
            delta = 1 + delta;
            g = ng;
            coeffs = Coeffs {
                uu: nu,
                vv: nv,
                qq: nq,
                rr: nr,
            };
        } else {
            // Case C:
            //   (δ, f, g) ← (1 + δ, f, g / 2)
            //   (U,V,Q,R) ← (2U, 2V, Q, R)
            let ng = g.shr1_even();
            let nu = addm(coeffs.uu, coeffs.uu, p);
            let nv = addm(coeffs.vv, coeffs.vv, p);
            delta = 1 + delta;
            g = ng;
            coeffs = Coeffs {
                uu: nu,
                vv: nv,
                qq: coeffs.qq,
                rr: coeffs.rr,
            };
        }

        let abs_delta = delta.unsigned_abs() as i64;
        if abs_delta > max_abs_delta {
            max_abs_delta = abs_delta;
        }
    }

    let iters_done = converged_iter.unwrap_or(max_iters);
    DivstepsRun {
        converged: converged_iter.is_some(),
        iters_done,
        max_abs_delta,
        final_f: f,
        final_g: g,
        final_coeffs: coeffs,
    }
}

/// Run exactly `iters` divsteps, continuing after convergence with the
/// `g = 0` even branch. Constant-time BY recip does this: once `g` is zero,
/// later steps only double the top coefficient row, preserving the fixed
/// invariant `2^iters f = U p + V g0`.
///
/// This is the right model for an approximate fixed-cap circuit: convergence
/// before the cap yields a valid inverse scaled by the public `2^-iters`; lack
/// of convergence is the permitted failure event.
pub fn run_divsteps_fixed(g0: U256, p: U256, iters: usize) -> DivstepsRun {
    assert!(p.bit(0), "p must be odd");
    assert!(g0 < p && !g0.is_zero(), "g0 must lie in [1, p)");

    let mut delta: i64 = 1;
    let mut f = SInt::from_u(p);
    let mut g = SInt::from_u(g0);
    let mut coeffs = Coeffs::initial();
    let mut max_abs_delta = 1i64;

    for _ in 0..iters {
        let g_odd = g.bit0();
        if delta > 0 && g_odd {
            let nf = g;
            let ng = SInt::sub(g, f).shr1_even();
            let nu = addm(coeffs.qq, coeffs.qq, p);
            let nv = addm(coeffs.rr, coeffs.rr, p);
            let nq = subm(coeffs.qq, coeffs.uu, p);
            let nr = subm(coeffs.rr, coeffs.vv, p);
            delta = 1 - delta;
            f = nf;
            g = ng;
            coeffs = Coeffs {
                uu: nu,
                vv: nv,
                qq: nq,
                rr: nr,
            };
        } else if g_odd {
            let ng = SInt::add(g, f).shr1_even();
            let nu = addm(coeffs.uu, coeffs.uu, p);
            let nv = addm(coeffs.vv, coeffs.vv, p);
            let nq = addm(coeffs.qq, coeffs.uu, p);
            let nr = addm(coeffs.rr, coeffs.vv, p);
            delta = 1 + delta;
            g = ng;
            coeffs = Coeffs {
                uu: nu,
                vv: nv,
                qq: nq,
                rr: nr,
            };
        } else {
            let ng = g.shr1_even();
            let nu = addm(coeffs.uu, coeffs.uu, p);
            let nv = addm(coeffs.vv, coeffs.vv, p);
            delta = 1 + delta;
            g = ng;
            coeffs = Coeffs {
                uu: nu,
                vv: nv,
                qq: coeffs.qq,
                rr: coeffs.rr,
            };
        }

        let abs_delta = delta.unsigned_abs() as i64;
        if abs_delta > max_abs_delta {
            max_abs_delta = abs_delta;
        }
    }

    DivstepsRun {
        converged: g.is_zero(),
        iters_done: iters,
        max_abs_delta,
        final_f: f,
        final_g: g,
        final_coeffs: coeffs,
    }
}

/// Recover `g0^{-1} mod p` from a converged divsteps run.
///
/// From the invariant `2^k f_k = U p + V g0`, with final `f_k = ±1`:
///
///     g0^{-1} ≡ sign(f_k) · V · 2^{-k}  (mod p)
pub fn recover_modinv(run: &DivstepsRun, p: U256) -> Option<U256> {
    if !run.converged {
        return None;
    }
    if !(run.final_f.is_one_pos() || run.final_f.is_one_neg()) {
        return None;
    }

    // 2^{-1} mod p = (p+1)/2 for odd p.
    let two_inv = (p.wrapping_add(U256::from(1))) >> 1;
    let mut two_inv_k = U256::from(1);
    let mut base = two_inv;
    let mut e = run.iters_done as u64;
    while e > 0 {
        if e & 1 == 1 {
            two_inv_k = mulm(two_inv_k, base, p);
        }
        e >>= 1;
        if e > 0 {
            base = mulm(base, base, p);
        }
    }
    let v_scaled = mulm(run.final_coeffs.vv, two_inv_k, p);
    if run.final_f.is_one_pos() {
        Some(v_scaled)
    } else {
        Some(negm(v_scaled, p))
    }
}

/// Fermat-little-theorem inverse for cross-checking.
pub fn fermat_modinv(a: U256, p: U256) -> U256 {
    assert!(!a.is_zero());
    let exp = p.wrapping_sub(U256::from(2));
    let mut result = U256::from(1);
    let mut base = a % p;
    for i in 0..256 {
        if exp.bit(i) {
            result = mulm(result, base, p);
        }
        base = mulm(base, base, p);
    }
    result
}

// ─────────────────────────────────────────────────────────────────────────
// Deterministic sampler for surveys
// ─────────────────────────────────────────────────────────────────────────

pub struct Sampler {
    reader: Box<dyn XofReader>,
    p: U256,
}

impl Sampler {
    pub fn new(seed: &[u8], p: U256) -> Self {
        let mut hasher = sha3::Shake128::default();
        hasher.update(seed);
        Self {
            reader: Box::new(hasher.finalize_xof()),
            p,
        }
    }

    pub fn next(&mut self) -> U256 {
        loop {
            let mut buf = [0u8; 32];
            self.reader.read(&mut buf);
            let x = U256::from_le_slice(&buf);
            if x < self.p && !x.is_zero() {
                return x;
            }
        }
    }
}

#[derive(Debug, Default, Clone)]
pub struct SurveyStats {
    pub samples: usize,
    pub all_converged: bool,
    pub min_iters: usize,
    pub max_iters: usize,
    pub sum_iters: u128,
    pub max_abs_delta: i64,
    pub modinv_matches: usize,
    pub modinv_mismatches: usize,
}

impl SurveyStats {
    pub fn mean_iters(&self) -> f64 {
        if self.samples == 0 {
            0.0
        } else {
            self.sum_iters as f64 / self.samples as f64
        }
    }
}

pub fn survey(sampler: &mut Sampler, n_samples: usize, p: U256, max_iters: usize) -> SurveyStats {
    let mut stats = SurveyStats {
        samples: 0,
        all_converged: true,
        min_iters: usize::MAX,
        max_iters: 0,
        sum_iters: 0,
        max_abs_delta: 0,
        modinv_matches: 0,
        modinv_mismatches: 0,
    };

    let deadline = two_min_deadline();
    for i in 0..n_samples {
        if (i & 127) == 0 {
            check_deadline(deadline, "by::survey");
        }
        let x = sampler.next();
        let run = run_divsteps(x, p, max_iters);
        if !run.converged {
            stats.all_converged = false;
        }
        let k = run.iters_done;
        stats.samples += 1;
        if k < stats.min_iters {
            stats.min_iters = k;
        }
        if k > stats.max_iters {
            stats.max_iters = k;
        }
        stats.sum_iters += k as u128;
        if run.max_abs_delta > stats.max_abs_delta {
            stats.max_abs_delta = run.max_abs_delta;
        }

        let expected = fermat_modinv(x, p);
        match recover_modinv(&run, p) {
            Some(v) if v == expected => stats.modinv_matches += 1,
            _ => stats.modinv_mismatches += 1,
        }
    }
    stats
}

// ─────────────────────────────────────────────────────────────────────────
// jumpdivsteps2 matrix survey
// ─────────────────────────────────────────────────────────────────────────
//
// BY 2019/266 Fig. 10.2 defines jumpdivsteps2 recursively. The returned
// matrix P satisfies
//
//     (f_n, g_n)^T = (1 / 2^n) · P · (f, g)^T
//
// and entries of P are bounded by 2^n in the worst case.
//
// For reversible quantum cost, what matters is the ACTUAL entry bit-width,
// because applying `a·f + b·g` costs roughly `(bitlen(a)+bitlen(b)) · n` in
// conditional-add/sub operations. So we measure the empirical distribution of
// entry sizes on random low-word inputs.

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct TransitionMatrix {
    pub m00: i128,
    pub m01: i128,
    pub m10: i128,
    pub m11: i128,
    pub delta_final: i64,
}

/// Truncate a signed integer to `t` bits as in BY Fig. 10.1:
///
///     truncate(f, t) = ((f + 2^{t-1}) mod 2^t) - 2^{t-1}
///
/// Here we operate on ordinary signed i128 for the low-word survey only.
pub fn truncate_i128(f: i128, t: usize) -> i128 {
    if t == 0 {
        return 0;
    }
    let two_t_minus_1: i128 = 1i128 << (t - 1);
    ((f + two_t_minus_1) & ((two_t_minus_1 << 1) - 1)) - two_t_minus_1
}

/// Classical Fig. 10.1 `divsteps2(n, t, delta, f, g)` on low-word signed ints.
/// Returns `(delta_n, f_n, g_n, matrix)`.
pub fn divsteps2_lowword(
    mut n: usize,
    mut t: usize,
    mut delta: i64,
    mut f: i128,
    mut g: i128,
) -> (i64, i128, i128, TransitionMatrix) {
    assert!(t >= n && n >= 1);
    f = truncate_i128(f, t);
    g = truncate_i128(g, t);
    let (mut u, mut v, mut q, mut r) = (1i128, 0i128, 0i128, 1i128);
    while n > 0 {
        f = truncate_i128(f, t);
        if delta > 0 && (g & 1) != 0 {
            let (ndelta, nf, ng, nu, nv, nq, nr) = (-delta, g, -f, q, r, -u, -v);
            delta = ndelta;
            f = nf;
            g = ng;
            u = nu;
            v = nv;
            q = nq;
            r = nr;
        }
        let g0 = g & 1;
        delta = 1 + delta;
        g = (g + g0 * f) / 2;
        q = (q + g0 * u) / 2;
        r = (r + g0 * v) / 2;
        n -= 1;
        t -= 1;
        g = truncate_i128(g, t);
    }
    (
        delta,
        f,
        g,
        TransitionMatrix {
            m00: u,
            m01: v,
            m10: q,
            m11: r,
            delta_final: delta,
        },
    )
}

/// Directly accumulate the integer 2×2 transition matrix over `w` divsteps.
///
/// If `P_w` is the returned matrix, then
///
///     (f_w, g_w)^T = (1 / 2^w) · P_w · (f_0, g_0)^T
///
/// where `(f_i, g_i)` are the states produced by BY `divstep` on the low-word
/// approximation. This is the quantity relevant to reversible cost: applying
/// `P_w` to the full-width quantum registers costs proportional to the bit-width
/// of the entries of `P_w`.
///
/// The low-word state evolution follows Fig. 10.1's `divsteps2`: after each
/// step, `t` shrinks by 1 and `g` is truncated to the new `t` bits; `f` is
/// truncated at the start of the next step. We mirror that behavior.
pub fn jump_matrix_direct_lowword(
    w: usize,
    mut t: usize,
    mut delta: i64,
    mut f: i128,
    mut g: i128,
) -> (i64, i128, i128, TransitionMatrix) {
    assert!(t >= w && w >= 1);
    // Integer matrices corresponding to the three branch cases, with the
    // common 1/2 factor pulled out:
    //  A: (f', g') = (g, (g-f)/2)     = (1/2) * [[0,2],[-1,1]] [f,g]
    //  B: (f', g') = (f, (g+f)/2)     = (1/2) * [[2,0],[ 1,1]] [f,g]
    //  C: (f', g') = (f, g/2)         = (1/2) * [[2,0],[ 0,1]] [f,g]
    let (mut p00, mut p01, mut p10, mut p11) = (1i128, 0i128, 0i128, 1i128);
    let mut n = w;
    f = truncate_i128(f, t);
    g = truncate_i128(g, t);
    while n > 0 {
        f = truncate_i128(f, t);
        if delta > 0 && (g & 1) != 0 {
            // Case A
            let (np00, np01, np10, np11) = (
                0 * p00 + 2 * p10,
                0 * p01 + 2 * p11,
                -1 * p00 + 1 * p10,
                -1 * p01 + 1 * p11,
            );
            let new_f = g;
            let new_g = (g - f) / 2;
            delta = 1 - delta;
            f = new_f;
            g = new_g;
            p00 = np00;
            p01 = np01;
            p10 = np10;
            p11 = np11;
        } else if (g & 1) != 0 {
            // Case B
            let (np00, np01, np10, np11) = (
                2 * p00 + 0 * p10,
                2 * p01 + 0 * p11,
                1 * p00 + 1 * p10,
                1 * p01 + 1 * p11,
            );
            let new_g = (g + f) / 2;
            delta = 1 + delta;
            g = new_g;
            p00 = np00;
            p01 = np01;
            p10 = np10;
            p11 = np11;
        } else {
            // Case C
            let (np00, np01, np10, np11) = (2 * p00, 2 * p01, p10, p11);
            let new_g = g / 2;
            delta = 1 + delta;
            g = new_g;
            p00 = np00;
            p01 = np01;
            p10 = np10;
            p11 = np11;
        }
        n -= 1;
        t -= 1;
        g = truncate_i128(g, t);
    }
    let f_out = truncate_i128(f, t + 1); // after n=w steps, f known to t-w+1 bits
    let g_out = truncate_i128(g, t); // and g to t-w bits. Here `t` already decremented.
    (
        delta,
        f_out,
        g_out,
        TransitionMatrix {
            m00: p00,
            m01: p01,
            m10: p10,
            m11: p11,
            delta_final: delta,
        },
    )
}

#[derive(Clone, Debug, Default)]
pub struct JumpStats {
    pub samples: usize,
    pub w: usize,
    pub max_entry_abs: i128,
    pub sum_log2_entry_abs: f64,
    pub nonzero_entries: usize,
}

pub fn jump_matrix_entry_survey(seed: &[u8], n_samples: usize, w: usize) -> JumpStats {
    let mut hasher = sha3::Shake128::default();
    hasher.update(seed);
    let mut reader = hasher.finalize_xof();
    let mut stats = JumpStats {
        samples: 0,
        w,
        max_entry_abs: 0,
        sum_log2_entry_abs: 0.0,
        nonzero_entries: 0,
    };
    let deadline = two_min_deadline();
    let mut buf = [0u8; 24];
    for i in 0..n_samples {
        if (i & 1023) == 0 {
            check_deadline(deadline, "by::jump_matrix_entry_survey");
        }
        reader.read(&mut buf);
        let mut f_low = u64::from_le_bytes(buf[0..8].try_into().unwrap()) as i128;
        let g_low = u64::from_le_bytes(buf[8..16].try_into().unwrap()) as i128;
        let delta = (u64::from_le_bytes(buf[16..24].try_into().unwrap()) % 41) as i64 - 20;
        f_low |= 1; // ensure odd
        let (_, _, _, m) = jump_matrix_direct_lowword(w, w, delta, f_low, g_low);
        for &e in &[m.m00, m.m01, m.m10, m.m11] {
            let abs = e.wrapping_abs();
            if abs > stats.max_entry_abs {
                stats.max_entry_abs = abs;
            }
            if abs > 0 {
                stats.sum_log2_entry_abs += (abs as f64).log2();
                stats.nonzero_entries += 1;
            }
        }
        stats.samples += 1;
    }
    stats
}

#[derive(Clone, Debug, Default)]
pub struct JumpHistogram {
    pub samples: usize,
    pub distinct_matrices: usize,
    pub most_common_count: usize,
    pub most_common_matrix: Option<TransitionMatrix>,
    pub total_unique_rows: usize,
}

/// Enumerate all possible low-word states for a given w and record how many
/// distinct transition matrices actually occur.
///
/// State space:
///   - delta in [-20, 20] (empirical |delta| cap from the 10k secp256k1 survey)
///   - f_low odd w-bit value
///   - g_low arbitrary w-bit value
///
/// This is the exact state space a fixed-width jumped-BY step would need to
/// handle if we bound delta to the observed range.
pub fn jump_matrix_histogram_all_states(w: usize) -> JumpHistogram {
    use std::collections::BTreeMap;
    let mut counts: BTreeMap<TransitionMatrix, usize> = BTreeMap::new();
    let f_states: usize = 1usize << (w - 1); // odd w-bit values
    let g_states: usize = 1usize << w;
    let mut samples = 0usize;
    for delta in -20i64..=20i64 {
        for f_odd in 0..f_states {
            let f_low: i128 = ((f_odd << 1) | 1) as i128;
            for g_raw in 0..g_states {
                let g_low: i128 = g_raw as i128;
                let (_, _, _, m) = jump_matrix_direct_lowword(w, w, delta, f_low, g_low);
                *counts.entry(m).or_insert(0) += 1;
                samples += 1;
            }
        }
    }
    let distinct_matrices = counts.len();
    let mut most_common_count = 0usize;
    let mut most_common_matrix = None;
    for (m, c) in &counts {
        if *c > most_common_count {
            most_common_count = *c;
            most_common_matrix = Some(*m);
        }
    }
    JumpHistogram {
        samples,
        distinct_matrices,
        most_common_count,
        most_common_matrix,
        total_unique_rows: counts.values().filter(|&&c| c == 1).count(),
    }
}

/// Count how many distinct low-w states can reach the *same* jump matrix.
///
/// If the number of distinct matrices is dramatically smaller than the state
/// space, a reversible implementation can use a QROM indexed by a compressed
/// class rather than by all (delta, f_low, g_low) tuples.

/// Env-gated smoke output used by `src/point_add/mod.rs` when BY_TEST=1.
pub fn run_classical_test() {
    let p = SECP256K1_P;
    let theoretical_bound = safegcd_iters(256);
    let max_iters = theoretical_bound + 100;
    let mut sampler = Sampler::new(b"divstep2-survey-seed-v1", p);
    let stats = survey(&mut sampler, 10_000, p, max_iters);

    eprintln!("=== B-Y divstep2 empirical survey on secp256k1 ===");
    eprintln!("samples            : {}", stats.samples);
    eprintln!("all_converged      : {}", stats.all_converged);
    eprintln!("theoretical bound  : {}", theoretical_bound);
    eprintln!("min iters observed : {}", stats.min_iters);
    eprintln!("max iters observed : {}", stats.max_iters);
    eprintln!("mean iters         : {:.2}", stats.mean_iters());
    eprintln!("max |δ| observed   : {}", stats.max_abs_delta);
    eprintln!("modinv matches     : {}", stats.modinv_matches);
    eprintln!("modinv mismatches  : {}", stats.modinv_mismatches);
    eprintln!("=================================================");

    for &w in &[4usize, 8, 12, 16] {
        let js = jump_matrix_entry_survey(b"jumpdivstep-matrix-seed-v1", 100_000, w);
        let mean_log2 = if js.nonzero_entries == 0 {
            0.0
        } else {
            js.sum_log2_entry_abs / (js.nonzero_entries as f64)
        };
        eprintln!("=== jumpdivstep matrix-entry survey (w={}) ===", w);
        eprintln!("samples                 : {}", js.samples);
        eprintln!("max |entry| observed    : {}", js.max_entry_abs);
        eprintln!(
            "max log2 |entry|        : {:.3}",
            (js.max_entry_abs as f64).log2()
        );
        eprintln!("mean log2 |entry|       : {:.3}", mean_log2);
        eprintln!("theoretical max log2    : {}", w);
        eprintln!("===========================================");
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn divstep_smoke() {
        let p = SECP256K1_P;
        let inputs: &[U256] = &[
            U256::from(1),
            U256::from(2),
            U256::from(3),
            U256::from(0xDEADBEEFu64),
            U256::from_limbs([
                0x0123456789ABCDEF,
                0xFEDCBA9876543210,
                0x0F0F0F0F0F0F0F0F,
                0x1234567890ABCDEF,
            ]),
            p.wrapping_sub(U256::from(1)),
        ];
        let max_iters = safegcd_iters(256);
        for &x in inputs {
            let run = run_divsteps(x, p, max_iters);
            assert!(run.converged, "did not converge for x={}", x);
            let got = recover_modinv(&run, p).expect("recovery");
            let expected = fermat_modinv(x, p);
            assert_eq!(got, expected, "modinv mismatch x={}", x);
        }
    }

    #[test]
    fn survey_10k() {
        let p = SECP256K1_P;
        let theoretical_bound = safegcd_iters(256);
        let max_iters = theoretical_bound + 100;
        let mut sampler = Sampler::new(b"divstep2-survey-seed-v1", p);
        let stats = survey(&mut sampler, 10_000, p, max_iters);

        eprintln!("=== B-Y divstep2 empirical survey on secp256k1 ===");
        eprintln!("samples            : {}", stats.samples);
        eprintln!("all_converged      : {}", stats.all_converged);
        eprintln!("theoretical bound  : {}", theoretical_bound);
        eprintln!("min iters observed : {}", stats.min_iters);
        eprintln!("max iters observed : {}", stats.max_iters);
        eprintln!("mean iters         : {:.2}", stats.mean_iters());
        eprintln!("max |δ| observed   : {}", stats.max_abs_delta);
        eprintln!("modinv matches     : {}", stats.modinv_matches);
        eprintln!("modinv mismatches  : {}", stats.modinv_mismatches);
        eprintln!("=================================================");

        assert!(stats.all_converged);
        assert_eq!(stats.modinv_mismatches, 0);
        assert!(
            stats.max_iters <= theoretical_bound,
            "observed max iters {} exceeds theoretical bound {}",
            stats.max_iters,
            theoretical_bound
        );
    }

    fn row_popcount_adds_i128(row: (i128, i128)) -> usize {
        let terms = row.0.unsigned_abs().count_ones() as usize
            + row.1.unsigned_abs().count_ones() as usize;
        terms.saturating_sub(1)
    }

    fn matrix_popcount_adds_i128(m: TransitionMatrix) -> usize {
        row_popcount_adds_i128((m.m00, m.m01)) + row_popcount_adds_i128((m.m10, m.m11))
    }

    #[test]
    fn approximate_divstep_cutoff_survey() {
        // With approximate failure tolerance, BY's empirical convergence tail
        // is much shorter than the 742-step proof bound. This matters because
        // jump windows scale directly with the iteration cap. Keep this as a
        // distributional fact, not as an exact-circuit claim.
        let p = SECP256K1_P;
        let samples = 20_000usize;
        let mut sampler = Sampler::new(b"by-approx-cutoff-v1", p);
        let mut iters = Vec::with_capacity(samples);
        for _ in 0..samples {
            let x = sampler.next();
            let run = run_divsteps(x, p, safegcd_iters(256));
            assert!(run.converged);
            iters.push(run.iters_done);
        }
        iters.sort_unstable();
        let q99 = iters[(samples * 99) / 100];
        let q999 = iters[(samples * 999) / 1000];
        let fail_550 = iters.iter().filter(|&&k| k > 550).count();
        let fail_560 = iters.iter().filter(|&&k| k > 560).count();
        eprintln!(
            "BY divstep cutoff: q99={q99}, q999={q999}, fail>550={:.4}, fail>560={:.4}, max={}",
            fail_550 as f64 / samples as f64,
            fail_560 as f64 / samples as f64,
            iters[samples - 1]
        );
        assert!(fail_550 as f64 / samples as f64 <= 0.01, "550-step approximate cutoff exceeded 1% on sample");
    }

    #[test]
    fn approximate_one_percent_cutoff_does_not_fund_lowword_by_near_miss() {
        // User-approved approximation budget: up to ~1% classical wrong
        // outputs may be tolerable, but phase and ancilla cleanup must remain
        // exact.  A fixed shorter BY cap is compatible with that: the circuit
        // can still be a clean reversible fixed-length circuit, and the only
        // failures are non-converged denominators.  Check whether this tolerance
        // by itself closes the old fully charged scratch600 BY near-miss.
        let p = SECP256K1_P;
        let samples = 50_000usize;
        let mut sampler = Sampler::new(b"by-approx-1pct-budget-v1", p);
        let mut iters = Vec::with_capacity(samples);
        for _ in 0..samples {
            let x = sampler.next();
            let run = run_divsteps(x, p, safegcd_iters(256));
            assert!(run.converged);
            iters.push(run.iters_done);
        }
        iters.sort_unstable();
        let fail_count = |cutoff: usize| -> usize { iters.iter().filter(|&&k| k > cutoff).count() };
        let fail_ppm = |cutoff: usize| -> usize { fail_count(cutoff) * 1_000_000 / samples };
        let mut best_cutoff = 560usize;
        for cutoff in 520usize..=560 {
            if fail_count(cutoff) * 100 <= samples {
                best_cutoff = cutoff;
                break;
            }
        }
        let mut aligned_cutoff = 560usize;
        for cutoff in (16usize..=560).step_by(16) {
            if fail_count(cutoff) * 100 <= samples {
                aligned_cutoff = cutoff;
                break;
            }
        }

        // Fully charged old low-scratch BY accounting decomposes exactly as:
        // scaffold_after_div + streamed replay + pattern-delta decoder +
        // lowword pattern oracle = 642716 + 3308*K + 111*K + 372*K.
        // K=560 reproduces the recorded 2,765,676 near-miss.  The any-cutoff
        // number grants a partial final window; the aligned number requires
        // whole 16-step windows.
        let scaffold_after_div = 642_716i64;
        let replay_per_step = 3_308i64;
        let decoder_per_step = 62_160i64 / 560;
        let lowword_selector_per_step = 208_320i64 / 560;
        let projected = |k: usize| -> i64 {
            scaffold_after_div + (replay_per_step + decoder_per_step + lowword_selector_per_step) * k as i64
        };
        let best_projected = projected(best_cutoff);
        let aligned_projected = projected(aligned_cutoff);
        let best_gap = best_projected - 2_700_000;
        let aligned_gap = aligned_projected - 2_700_000;
        println!("METRIC by_approx_1pct_best_cutoff_steps={best_cutoff}");
        println!("METRIC by_approx_1pct_best_fail_ppm={}", fail_ppm(best_cutoff));
        println!("METRIC by_approx_1pct_best_projected_gap_ccx={best_gap}");
        println!("METRIC by_approx_1pct_aligned_cutoff_steps={aligned_cutoff}");
        println!("METRIC by_approx_1pct_aligned_projected_gap_ccx={aligned_gap}");
        println!("METRIC by_approx_fail544_ppm={}", fail_ppm(544));
        println!("METRIC by_approx_fail550_ppm={}", fail_ppm(550));
        println!("METRIC by_approx_fail560_ppm={}", fail_ppm(560));
        eprintln!(
            "BY 1% approximate cutoff budget: best_cutoff={best_cutoff}, fail_ppm={}, projected={best_projected}, gap={best_gap}; aligned_cutoff={aligned_cutoff}, aligned_gap={aligned_gap}; fail544={}ppm fail550={}ppm fail560={}ppm",
            fail_ppm(best_cutoff),
            fail_ppm(544),
            fail_ppm(550),
            fail_ppm(560)
        );
        assert!(fail_ppm(best_cutoff) <= 10_000, "selected cutoff violates 1% classical-mismatch budget");
    }

    #[test]
    fn approximate_one_percent_pointadd_budget_needs_two_denominators() {
        // The benchmark-level tolerance is on point-add outputs, not on a
        // single denominator.  A two-DIV affine add therefore needs each fixed
        // cutoff failure rate to be roughly <=0.5% (union-bound accounting),
        // unless a later route proves strong correlation or an error detector.
        let p = SECP256K1_P;
        let samples = 50_000usize;
        let mut sampler = Sampler::new(b"by-approx-1pct-pointadd-budget-v1", p);
        let mut iters = Vec::with_capacity(samples);
        for _ in 0..samples {
            let x = sampler.next();
            let run = run_divsteps(x, p, safegcd_iters(256));
            assert!(run.converged);
            iters.push(run.iters_done);
        }
        iters.sort_unstable();
        let fail_count = |cutoff: usize| -> usize { iters.iter().filter(|&&k| k > cutoff).count() };
        let fail_ppm = |cutoff: usize| -> usize { fail_count(cutoff) * 1_000_000 / samples };
        let mut cutoff_half_percent = 560usize;
        for cutoff in 520usize..=560 {
            if fail_count(cutoff) * 200 <= samples {
                cutoff_half_percent = cutoff;
                break;
            }
        }
        let mut cutoff_one_percent_two_den_union = 560usize;
        for cutoff in 520usize..=560 {
            // Exact independent estimate is 1-(1-p)^2; use integer ppm for a
            // diagnostic while selecting by the safer union-bound above.
            let ppm = fail_ppm(cutoff) as u128;
            let two_ppm = 2 * ppm - (ppm * ppm) / 1_000_000;
            if two_ppm <= 10_000 {
                cutoff_one_percent_two_den_union = cutoff;
                break;
            }
        }
        let scaffold_after_div = 642_716i64;
        let replay_per_step = 3_308i64;
        let decoder_per_step = 62_160i64 / 560;
        let lowword_selector_per_step = 208_320i64 / 560;
        let projected = |k: usize| -> i64 {
            scaffold_after_div + (replay_per_step + decoder_per_step + lowword_selector_per_step) * k as i64
        };
        let pointadd_gap = projected(cutoff_half_percent) - 2_700_000;
        let independent_cutoff_gap = projected(cutoff_one_percent_two_den_union) - 2_700_000;
        let per_den_ppm = fail_ppm(cutoff_half_percent);
        let union_ppm = 2 * per_den_ppm;
        let indep_ppm = {
            let ppm = per_den_ppm as u128;
            (2 * ppm - (ppm * ppm) / 1_000_000) as usize
        };
        println!("METRIC by_approx_pointadd_cutoff_steps={cutoff_half_percent}");
        println!("METRIC by_approx_pointadd_per_den_fail_ppm={per_den_ppm}");
        println!("METRIC by_approx_pointadd_union_fail_ppm={union_ppm}");
        println!("METRIC by_approx_pointadd_independent_fail_ppm={indep_ppm}");
        println!("METRIC by_approx_pointadd_projected_gap_ccx={pointadd_gap}");
        println!("METRIC by_approx_pointadd_independent_cutoff_steps={cutoff_one_percent_two_den_union}");
        println!("METRIC by_approx_pointadd_independent_gap_ccx={independent_cutoff_gap}");
        println!("METRIC by_approx_pointadd_fail550_ppm={}", fail_ppm(550));
        println!("METRIC by_approx_pointadd_fail552_ppm={}", fail_ppm(552));
        eprintln!(
            "BY 1% point-add cutoff budget: union cutoff={cutoff_half_percent}, per_den_fail={per_den_ppm}ppm, union_fail={union_ppm}ppm, gap={pointadd_gap}; independent cutoff={cutoff_one_percent_two_den_union}, independent_gap={independent_cutoff_gap}"
        );
        assert!(union_ppm <= 10_000, "two-denominator union-bound failure exceeds 1%");
    }

    #[test]
    fn harness_scale_approx_cutoff_leaves_by_lowword_over_budget() {
        // The real benchmark executes 9024 random point-add cases and expects
        // zero wrong classical outputs.  Percent-level approximate failures are
        // therefore unusable.  Use a harness-scale tail target and verify that
        // the old lowword-selector BY near-miss no longer gets meaningful step
        // savings once this stricter requirement is respected.
        let p = SECP256K1_P;
        let samples = 100_000usize;
        let mut sampler = Sampler::new(b"by-harness-scale-cutoff-v1", p);
        let mut iters = Vec::with_capacity(samples);
        for _ in 0..samples {
            let x = sampler.next();
            let run = run_divsteps(x, p, safegcd_iters(256));
            assert!(run.converged);
            iters.push(run.iters_done);
        }
        iters.sort_unstable();
        let fail_count = |cutoff: usize| -> usize { iters.iter().filter(|&&k| k > cutoff).count() };
        let fail_ppm = |cutoff: usize| -> usize { fail_count(cutoff) * 1_000_000 / samples };
        let pointadd_union_ppm = |cutoff: usize| -> usize { 2 * fail_ppm(cutoff) };
        let mut cutoff_100ppm_pointadd = 576usize;
        let mut cutoff_20ppm_pointadd = 576usize;
        for cutoff in 540usize..=576 {
            let union_ppm = pointadd_union_ppm(cutoff);
            if union_ppm <= 100 && cutoff_100ppm_pointadd == 576 {
                cutoff_100ppm_pointadd = cutoff;
            }
            if union_ppm <= 20 {
                cutoff_20ppm_pointadd = cutoff;
                break;
            }
        }
        let scaffold_after_div = 642_716i64;
        let replay_per_step = 3_308i64;
        let decoder_per_step = 62_160i64 / 560;
        let lowword_selector_per_step = 208_320i64 / 560;
        let projected = |k: usize| -> i64 {
            scaffold_after_div + (replay_per_step + decoder_per_step + lowword_selector_per_step) * k as i64
        };
        let gap_100ppm = projected(cutoff_100ppm_pointadd) - 2_700_000;
        let gap_20ppm = projected(cutoff_20ppm_pointadd) - 2_700_000;
        println!("METRIC by_approx_harness_100ppm_cutoff_steps={cutoff_100ppm_pointadd}");
        println!("METRIC by_approx_harness_100ppm_union_fail_ppm={}", pointadd_union_ppm(cutoff_100ppm_pointadd));
        println!("METRIC by_approx_harness_100ppm_projected_gap_ccx={gap_100ppm}");
        println!("METRIC by_approx_harness_20ppm_cutoff_steps={cutoff_20ppm_pointadd}");
        println!("METRIC by_approx_harness_20ppm_union_fail_ppm={}", pointadd_union_ppm(cutoff_20ppm_pointadd));
        println!("METRIC by_approx_harness_projected_gap_ccx={gap_20ppm}");
        println!("METRIC by_approx_harness_fail556_ppm={}", fail_ppm(556));
        println!("METRIC by_approx_harness_fail560_ppm={}", fail_ppm(560));
        println!("METRIC by_approx_harness_fail564_ppm={}", fail_ppm(564));
        println!("METRIC by_approx_harness_fail568_ppm={}", fail_ppm(568));
        eprintln!(
            "BY harness-scale cutoff: cutoff100ppm={cutoff_100ppm_pointadd}, gap100ppm={gap_100ppm}, cutoff20ppm={cutoff_20ppm_pointadd}, gap20ppm={gap_20ppm}, fail560={}ppm, fail564={}ppm",
            fail_ppm(560),
            fail_ppm(564)
        );
        assert!(gap_20ppm > 0, "harness-scale approximate cutoff unexpectedly funds old lowword BY route");
    }

    fn approx_odd_pattern_from_truncated_lowbits_for_test(
        w: usize,
        mut t: usize,
        mut delta: i64,
        f: SInt,
        g: SInt,
    ) -> u16 {
        let mut ff = sint_low_i128(f, t);
        let mut gg = sint_low_i128(g, t);
        let mut pattern = 0u16;
        for i in 0..w {
            ff = truncate_i128(ff, t);
            gg = truncate_i128(gg, t);
            let odd = (gg & 1) != 0;
            if odd {
                pattern |= 1u16 << i;
            }
            if delta > 0 && odd {
                let nf = gg;
                let ng = (gg - ff) / 2;
                delta = 1 - delta;
                ff = nf;
                gg = ng;
            } else if odd {
                let ng = (gg + ff) / 2;
                delta = 1 + delta;
                gg = ng;
            } else {
                let ng = gg / 2;
                delta = 1 + delta;
                gg = ng;
            }
            t = t.saturating_sub(1);
        }
        pattern
    }

    #[test]
    fn truncated_lowword_patterns_are_not_harness_clean() {
        // Gidney-style truncation works for high-bit residue accumulation
        // because the algorithm masks approximation error.  BY branch patterns
        // are the opposite: they are 2-adic low-bit controls.  If we copy fewer
        // than 16 low bits and guess the rest of a 16-step window, the selected
        // microprogram changes too often to survive a 9024-case harness.
        const W: usize = 16;
        const WINDOWS: usize = 35;
        let samples = 20_000usize;
        let mut sampler = Sampler::new(b"by-truncated-lowword-pattern-v1", SECP256K1_P);
        let ts = [8usize, 10, 12, 14, 15];
        let mut window_mismatches = vec![0usize; ts.len()];
        let mut trajectory_mismatches = vec![0usize; ts.len()];
        for _ in 0..samples {
            let x = sampler.next();
            let mut delta = 1i64;
            let mut f = SInt::from_u(SECP256K1_P);
            let mut g = SInt::from_u(x);
            let mut any = vec![false; ts.len()];
            for _ in 0..WINDOWS {
                let start_delta = delta;
                let start_f = f;
                let start_g = g;
                let mut true_pattern = 0u16;
                for i in 0..W {
                    if g.bit0() {
                        true_pattern |= 1u16 << i;
                    }
                    divstep_sint_state(&mut delta, &mut f, &mut g);
                }
                for (j, &t) in ts.iter().enumerate() {
                    let approx = approx_odd_pattern_from_truncated_lowbits_for_test(W, t, start_delta, start_f, start_g);
                    if approx != true_pattern {
                        window_mismatches[j] += 1;
                        any[j] = true;
                    }
                }
            }
            for (j, hit) in any.into_iter().enumerate() {
                if hit {
                    trajectory_mismatches[j] += 1;
                }
            }
        }
        let ppm = |count: usize, total: usize| -> usize { count * 1_000_000 / total };
        let total_windows = samples * WINDOWS;
        let window_t8 = ppm(window_mismatches[0], total_windows);
        let window_t12 = ppm(window_mismatches[2], total_windows);
        let window_t14 = ppm(window_mismatches[3], total_windows);
        let window_t15 = ppm(window_mismatches[4], total_windows);
        let traj_t14 = ppm(trajectory_mismatches[3], samples);
        let traj_t15 = ppm(trajectory_mismatches[4], samples);
        let selector_t14 = (5_952usize * 14).div_ceil(16) * WINDOWS;
        let selector_t15 = (5_952usize * 15).div_ceil(16) * WINDOWS;
        let selector_saving_t14 = 5_952usize * WINDOWS - selector_t14;
        let selector_saving_t15 = 5_952usize * WINDOWS - selector_t15;
        println!("METRIC by_trunc_lowword_window_mismatch_t8_ppm={window_t8}");
        println!("METRIC by_trunc_lowword_window_mismatch_t12_ppm={window_t12}");
        println!("METRIC by_trunc_lowword_window_mismatch_t14_ppm={window_t14}");
        println!("METRIC by_trunc_lowword_window_mismatch_t15_ppm={window_t15}");
        println!("METRIC by_trunc_lowword_trajectory_mismatch_t14_ppm={traj_t14}");
        println!("METRIC by_trunc_lowword_trajectory_mismatch_t15_ppm={traj_t15}");
        println!("METRIC by_trunc_lowword_selector_saving_t14_ccx={selector_saving_t14}");
        println!("METRIC by_trunc_lowword_selector_saving_t15_ccx={selector_saving_t15}");
        eprintln!(
            "BY truncated lowword patterns: window mismatch ppm t8={window_t8}, t12={window_t12}, t14={window_t14}, t15={window_t15}; trajectory ppm t14={traj_t14}, t15={traj_t15}; selector savings t14={selector_saving_t14}, t15={selector_saving_t15}"
        );
        assert!(traj_t15 > 10_000, "15-bit truncated lowword pattern would be harness-clean enough to revisit");
        assert!(selector_saving_t15 < 20_000, "15-bit truncation unexpectedly saves enough selector cost");
    }

    fn two_inv_pow(p: U256, iters: usize) -> U256 {
        let two_inv = (p.wrapping_add(U256::from(1))) >> 1;
        let mut acc = U256::from(1);
        let mut base = two_inv;
        let mut e = iters as u64;
        while e > 0 {
            if (e & 1) != 0 {
                acc = mulm(acc, base, p);
            }
            e >>= 1;
            if e != 0 {
                base = mulm(base, base, p);
            }
        }
        acc
    }

    #[test]
    fn fixed_by_coeff_channel_is_tagged_div_when_converged() {
        // Structural algebra for replacing Kaliski tagged-DIV with BY:
        // after fixed K divsteps, if f=±1 and g=0, the top coefficient V obeys
        //     V*x = sign(f)*2^K  (mod p),
        // and the bottom coefficient R obeys
        //     R*x = 0            (mod p)  -> R=0 for nonzero x.
        // Therefore carrying a tagged numerator y+x through the same
        // coefficient channel gives V*(y+x); multiplying by sign(f)*2^-K and
        // subtracting 1 recovers y/x, while the bottom channel is zero. This is
        // the BY analogue of the Kaliski y+x tagged DIV transform.
        let p = SECP256K1_P;
        let k = 550usize;
        let two_inv_k = two_inv_pow(p, k);
        let samples = 5_000usize;
        let mut sx = Sampler::new(b"by-fixed-tagged-div-x-v1", p);
        let mut sy = Sampler::new(b"by-fixed-tagged-div-y-v1", p);
        let mut failures = 0usize;
        for _ in 0..samples {
            let x = sx.next();
            let y = sy.next();
            let run = run_divsteps_fixed(x, p, k);
            if !run.converged || !(run.final_f.is_one_pos() || run.final_f.is_one_neg()) {
                failures += 1;
                continue;
            }
            let tag = addm(y, x, p);
            assert_eq!(mulm(run.final_coeffs.rr, tag, p), U256::ZERO, "bottom BY tagged channel did not self-zero");
            let raw = mulm(run.final_coeffs.vv, tag, p);
            let scaled = mulm(raw, two_inv_k, p);
            let plus_one = if run.final_f.is_one_pos() { scaled } else { negm(scaled, p) };
            let quotient = subm(plus_one, U256::from(1), p);
            let expected = mulm(y, fermat_modinv(x, p), p);
            assert_eq!(quotient, expected, "BY tagged quotient mismatch");
        }
        let fail_rate = failures as f64 / samples as f64;
        eprintln!(
            "fixed BY tagged-DIV algebra at K={k}: failures={failures}/{samples} ({fail_rate:.4})"
        );
        assert!(fail_rate <= 0.01, "550-step fixed BY tagged DIV exceeded 1% failure tolerance");
    }

    fn sint_low_i128(x: SInt, w: usize) -> i128 {
        let mask = if w == 64 { u64::MAX } else { (1u64 << w) - 1 };
        let low = (x.mag.as_limbs()[0] & mask) as i128;
        let signed = if x.neg { -low } else { low };
        truncate_i128(signed, w)
    }

    fn divstep_sint_state(delta: &mut i64, f: &mut SInt, g: &mut SInt) {
        let g_odd = g.bit0();
        if *delta > 0 && g_odd {
            let nf = *g;
            let ng = SInt::sub(*g, *f).shr1_even();
            *delta = 1 - *delta;
            *f = nf;
            *g = ng;
        } else if g_odd {
            let ng = SInt::add(*g, *f).shr1_even();
            *delta = 1 + *delta;
            *g = ng;
        } else {
            let ng = g.shr1_even();
            *delta = 1 + *delta;
            *g = ng;
        }
    }

    fn divstep_i128_once_for_ambiguity(delta: i64, f: i128, g: i128) -> (i64, i128, i128, u8) {
        let odd = (g & 1) != 0;
        if delta > 0 && odd {
            (1 - delta, g, (g - f) / 2, b'A')
        } else if odd {
            (1 + delta, f, (g + f) / 2, b'B')
        } else {
            (1 + delta, f, g / 2, b'C')
        }
    }

    fn inverse_divstep_candidates_for_ambiguity(delta: i64, f: i128, g: i128) -> Vec<(i64, i128, i128, u8)> {
        // Invert one BY denominator divstep from a poststate.  The point of the
        // last-shot consumed-denominator idea was to avoid storing branch bits
        // and later recover them from the consumed denominator state.  These are
        // the three mathematically valid predecessor formulas.
        let mut out = Vec::with_capacity(3);
        let candidates = [
            (delta - 1, f, 2 * g, b'C'),
            (delta - 1, f, 2 * g - f, b'B'),
            (1 - delta, f - 2 * g, f, b'A'),
        ];
        for &(d0, f0, g0, case) in &candidates {
            if (f0 & 1) == 0 {
                continue;
            }
            let valid = match case {
                b'C' => (g0 & 1) == 0,
                b'B' => (g0 & 1) != 0 && d0 <= 0,
                b'A' => (g0 & 1) != 0 && d0 > 0,
                _ => false,
            };
            if valid {
                let (d1, f1, g1, got_case) = divstep_i128_once_for_ambiguity(d0, f0, g0);
                assert_eq!((d1, f1, g1, got_case), (delta, f, g, case));
                out.push((d0, f0, g0, case));
            }
        }
        out
    }

    #[test]
    fn consumed_denominator_window_branchless_recovery_is_exponentially_ambiguous() {
        // Ground-up last-shot BY invalidation: if we try to consume the
        // denominator state in-place and clear the window branch controls from
        // the post-window denominator alone, the inverse relation is massively
        // many-to-one.  Even the tiny poststate (delta,f,g)=(0,1,0) has more
        // than half a million valid 16-step predecessor branch patterns.  Thus
        // a SOTA BY denominator cannot be branchless-poststate cleanup; it must
        // carry compressed history or consume q/pattern inside a reversible
        // fixed-matrix update with an explicit local inverse.
        let mut states: std::collections::BTreeMap<(i64, i128, i128), usize> = std::collections::BTreeMap::new();
        states.insert((0, 1, 0), 1);
        let mut totals = Vec::new();
        for _depth in 1..=16 {
            let mut next: std::collections::BTreeMap<(i64, i128, i128), usize> = std::collections::BTreeMap::new();
            for (&(d, f, g), &count) in states.iter() {
                for (d0, f0, g0, _case) in inverse_divstep_candidates_for_ambiguity(d, f, g) {
                    assert!(d0.abs() < 64);
                    assert!(f0.abs() < (1i128 << 32));
                    assert!(g0.abs() < (1i128 << 32));
                    *next.entry((d0, f0, g0)).or_insert(0) += count;
                }
            }
            let total: usize = next.values().copied().sum();
            totals.push(total);
            states = next;
        }
        let total4 = totals[3];
        let total8 = totals[7];
        let total16 = *totals.last().unwrap();
        let bits4 = usize::BITS as usize - (total4 - 1).leading_zeros() as usize;
        let bits8 = usize::BITS as usize - (total8 - 1).leading_zeros() as usize;
        let bits16 = usize::BITS as usize - (total16 - 1).leading_zeros() as usize;
        eprintln!(
            "BY consumed-denominator branchless poststate ambiguity: totals_by_depth={totals:?}, depth4_patterns={total4}, depth8_patterns={total8}, depth16_states={}, depth16_patterns={total16}",
            states.len()
        );
        println!("METRIC scratch600_consumed_w4_poststate_patterns={total4}");
        println!("METRIC scratch600_consumed_w4_poststate_bits={bits4}");
        println!("METRIC scratch600_consumed_w8_poststate_patterns={total8}");
        println!("METRIC scratch600_consumed_w8_poststate_bits={bits8}");
        println!("METRIC scratch600_consumed_w16_poststate_states={}", states.len());
        println!("METRIC scratch600_consumed_w16_poststate_patterns={total16}");
        println!("METRIC scratch600_consumed_w16_poststate_bits={bits16}");
        assert!(states.len() > 500_000, "poststate unexpectedly identifies most predecessor states");
        assert!(total16 > 500_000, "poststate unexpectedly identifies most branch patterns");
    }

    #[test]
    fn windowed_scaled_by_tagged_division_matches_microstep_algebra() {
        // Full classical model of the intended w=16 BY tagged-DIV route:
        // denominator evolves by exact 16 divsteps/window, while the tagged
        // modular channel applies 2^-16 P each window. After 35 windows (560
        // steps), convergence failures are far below 1%, and output recovery is
        // simply sign(f)*r - 1 because the 2^-K scaling has been paid per window.
        let p = SECP256K1_P;
        let w = 16usize;
        let windows = 35usize;
        let inv_scale = two_inv_pow(p, w);
        let samples = 3_000usize;
        let mut sx = Sampler::new(b"by-windowed-scaled-div-x-v1", p);
        let mut sy = Sampler::new(b"by-windowed-scaled-div-y-v1", p);
        let mut failures = 0usize;
        for _ in 0..samples {
            let x = sx.next();
            let y = sy.next();
            let mut delta = 1i64;
            let mut f = SInt::from_u(p);
            let mut g = SInt::from_u(x);
            let mut r = U256::ZERO;
            let mut s = addm(y, x, p);
            for _ in 0..windows {
                let f_low = sint_low_i128(f, w);
                let g_low = sint_low_i128(g, w);
                let (_, _, _, m) = jump_matrix_direct_lowword(w, w, delta, f_low, g_low);
                let nr = mulm(
                    addm(
                        mulm(signed_i128_mod_p(m.m00, p), r, p),
                        mulm(signed_i128_mod_p(m.m01, p), s, p),
                        p,
                    ),
                    inv_scale,
                    p,
                );
                let ns = mulm(
                    addm(
                        mulm(signed_i128_mod_p(m.m10, p), r, p),
                        mulm(signed_i128_mod_p(m.m11, p), s, p),
                        p,
                    ),
                    inv_scale,
                    p,
                );
                r = nr;
                s = ns;
                for _ in 0..w {
                    divstep_sint_state(&mut delta, &mut f, &mut g);
                }
            }
            if !g.is_zero() || !(f.is_one_pos() || f.is_one_neg()) {
                failures += 1;
                continue;
            }
            assert_eq!(s, U256::ZERO, "scaled BY bottom tagged channel did not zero");
            let plus_one = if f.is_one_pos() { r } else { negm(r, p) };
            let quotient = subm(plus_one, U256::from(1), p);
            let expected = mulm(y, fermat_modinv(x, p), p);
            assert_eq!(quotient, expected, "windowed scaled BY quotient mismatch");
        }
        let fail_rate = failures as f64 / samples as f64;
        eprintln!(
            "windowed scaled BY tagged DIV: windows={windows}, steps={}, failures={failures}/{samples} ({fail_rate:.4})",
            windows * w
        );
        assert!(fail_rate <= 0.01);
    }

    #[test]
    fn jumpdivstep_matrix_arithmetic_intensity_model() {
        // BY/jumpdivsteps is attractive because branch selection is local to
        // low words + delta, not a full-width u>v comparator. The price is a
        // selected signed 2x2 matrix. This row-popcount model estimates the
        // shifted add/sub terms needed to apply that matrix to one full-width
        // pair. It is not a complete circuit cost, but it is the right first
        // lower-bound for deciding if BY deserves a live prototype.
        let samples = 50_000usize;
        let mut hasher = sha3::Shake128::default();
        hasher.update(b"by-jump-matrix-popcount-v1");
        let mut reader = hasher.finalize_xof();
        let mut buf = [0u8; 24];
        for &w in &[4usize, 8, 12, 16] {
            let mut total = 0usize;
            let mut max_cost = 0usize;
            let mut costs = Vec::with_capacity(samples);
            for _ in 0..samples {
                reader.read(&mut buf);
                let f_low = (u64::from_le_bytes(buf[0..8].try_into().unwrap()) as i128) | 1;
                let g_low = u64::from_le_bytes(buf[8..16].try_into().unwrap()) as i128;
                let delta = (u64::from_le_bytes(buf[16..24].try_into().unwrap()) % 41) as i64 - 20;
                let (_, _, _, m) = jump_matrix_direct_lowword(w, w, delta, f_low, g_low);
                let c = matrix_popcount_adds_i128(m);
                total += c;
                max_cost = max_cost.max(c);
                costs.push(c);
            }
            costs.sort_unstable();
            let mean = total as f64 / samples as f64;
            let p90 = costs[(samples * 90) / 100];
            let exact_windows = safegcd_iters(256).div_ceil(w);
            let mean_terms_per_pair = mean * exact_windows as f64;
            eprintln!(
                "BY jump w={w}: mean row-add terms/window={mean:.2}, p90={p90}, max={max_cost}, exact_windows={}, mean_terms_per_pair={mean_terms_per_pair:.1}",
                exact_windows
            );
            assert!(mean_terms_per_pair < 600.0, "BY row-add intensity unexpectedly high");
        }
    }

    #[test]
    fn jumpdivstep_budget_model_suggests_live_prototype() {
        // Very optimistic but actionable budget model for a BY jump inversion:
        // apply the selected 2x2 matrix to three full-width pairs:
        //   (f,g) plus the two coefficient columns. Each row-popcount term is
        // charged as one n-bit add/sub. This ignores reversible matrix synthesis,
        // sign handling, reductions, and cleanup, so it is a lower bound; still,
        // if this were already > Kaliski there would be no reason to prototype.
        const N: usize = 256;
        const PAIRS_PER_WINDOW: usize = 3;
        let samples = 50_000usize;
        let mut hasher = sha3::Shake128::default();
        hasher.update(b"by-jump-budget-v1");
        let mut reader = hasher.finalize_xof();
        let mut buf = [0u8; 24];
        let w = 16usize;
        let mut total_terms = 0usize;
        for _ in 0..samples {
            reader.read(&mut buf);
            let f_low = (u64::from_le_bytes(buf[0..8].try_into().unwrap()) as i128) | 1;
            let g_low = u64::from_le_bytes(buf[8..16].try_into().unwrap()) as i128;
            let delta = (u64::from_le_bytes(buf[16..24].try_into().unwrap()) % 41) as i64 - 20;
            let (_, _, _, m) = jump_matrix_direct_lowword(w, w, delta, f_low, g_low);
            total_terms += matrix_popcount_adds_i128(m);
        }
        let mean_terms_per_window = total_terms as f64 / samples as f64;
        let exact_windows = safegcd_iters(256).div_ceil(w);
        let approx_windows_1pct = 550usize.div_ceil(w);
        let exact_toffoli_lb = mean_terms_per_window * exact_windows as f64 * PAIRS_PER_WINDOW as f64 * N as f64;
        let approx_toffoli_lb = mean_terms_per_window * approx_windows_1pct as f64 * PAIRS_PER_WINDOW as f64 * N as f64;
        eprintln!(
            "BY w=16 budget lower-bound: mean_terms/window={mean_terms_per_window:.2}, exact_windows={exact_windows}, exact≈{exact_toffoli_lb:.0} Toffoli, approx_windows={approx_windows_1pct}, approx≈{approx_toffoli_lb:.0} Toffoli"
        );
        assert!(exact_toffoli_lb < 600_000.0, "BY lower bound no longer beats Kaliski enough to prototype");
        assert!(approx_toffoli_lb < 500_000.0, "Approx BY lower bound too high");
    }

    fn count_ccx(ops: &[crate::circuit::Op]) -> usize {
        ops.iter()
            .filter(|o| matches!(o.kind, crate::circuit::OperationType::CCX | crate::circuit::OperationType::CCZ))
            .count()
    }

    fn add_shifted_term_for_cost(
        b: &mut super::super::B,
        src: &[super::super::QubitId],
        dst: &[super::super::QubitId],
        shift: usize,
        subtract: bool,
    ) {
        if shift >= dst.len() {
            return;
        }
        let len = dst.len() - shift;
        let addend = b.alloc_qubits(len);
        let copy_len = src.len().min(len);
        for i in 0..copy_len {
            b.cx(src[i], addend[i]);
        }
        let dst_slice: Vec<_> = dst[shift..shift + len].to_vec();
        if subtract {
            super::super::sub_nbit_qq_fast(b, &addend, &dst_slice);
        } else {
            super::super::add_nbit_qq_fast(b, &addend, &dst_slice);
        }
        for i in 0..copy_len {
            b.cx(src[i], addend[i]);
        }
        b.free_vec(&addend);
    }

    fn add_coeff_times_for_cost(
        b: &mut super::super::B,
        coeff: i128,
        src: &[super::super::QubitId],
        dst: &[super::super::QubitId],
    ) {
        let subtract = coeff < 0;
        let mut mag = coeff.unsigned_abs();
        let mut shift = 0usize;
        while mag != 0 {
            if (mag & 1) != 0 {
                add_shifted_term_for_cost(b, src, dst, shift, subtract);
            }
            mag >>= 1;
            shift += 1;
        }
    }

    fn emit_constant_matrix_apply_for_cost(b: &mut super::super::B, m: TransitionMatrix, width: usize) {
        let f = b.alloc_qubits(width);
        let g = b.alloc_qubits(width);
        let out0 = b.alloc_qubits(width);
        let out1 = b.alloc_qubits(width);
        add_coeff_times_for_cost(b, m.m00, &f, &out0);
        add_coeff_times_for_cost(b, m.m01, &g, &out0);
        add_coeff_times_for_cost(b, m.m10, &f, &out1);
        add_coeff_times_for_cost(b, m.m11, &g, &out1);
        // This is only a forward cost/peak probe for row formation; outputs are
        // not freed because the full BY state update would swap/use them.
        let _ = (f, g, out0, out1);
    }

    fn det_sign_pow2(m: TransitionMatrix, w: usize) -> i128 {
        let det = m.m00 * m.m11 - m.m01 * m.m10;
        let scale = 1i128 << w;
        assert!(det == scale || det == -scale, "unexpected jump determinant {det}, expected ±{scale}");
        det / scale
    }

    fn scaled_inverse_matrix(m: TransitionMatrix, w: usize) -> TransitionMatrix {
        // For new = P old / 2^w and det(P)=s·2^w, old = s·adj(P) new.
        let s = det_sign_pow2(m, w);
        TransitionMatrix {
            m00: s * m.m11,
            m01: -s * m.m01,
            m10: -s * m.m10,
            m11: s * m.m00,
            delta_final: m.delta_final,
        }
    }

    fn emit_scaled_pair_update_with_cleanup_for_cost(
        b: &mut super::super::B,
        m: TransitionMatrix,
        width: usize,
        w: usize,
    ) {
        // More faithful BY jump pair update cost:
        //   temp = P·old is accumulated at width+w bits;
        //   temp low w bits are mathematically zero;
        //   new is the high `width` bits, i.e. P·old / 2^w;
        //   old is cleaned using old = (2^w/det(P)) adj(P) new.
        let f = b.alloc_qubits(width);
        let g = b.alloc_qubits(width);
        let tmp0 = b.alloc_qubits(width + w);
        let tmp1 = b.alloc_qubits(width + w);

        add_coeff_times_for_cost(b, m.m00, &f, &tmp0);
        add_coeff_times_for_cost(b, m.m01, &g, &tmp0);
        add_coeff_times_for_cost(b, m.m10, &f, &tmp1);
        add_coeff_times_for_cost(b, m.m11, &g, &tmp1);

        let new0 = tmp0[w..w + width].to_vec();
        let new1 = tmp1[w..w + width].to_vec();
        let inv = scaled_inverse_matrix(m, w);
        add_coeff_times_for_cost(b, -inv.m00, &new0, &f);
        add_coeff_times_for_cost(b, -inv.m01, &new1, &f);
        add_coeff_times_for_cost(b, -inv.m10, &new0, &g);
        add_coeff_times_for_cost(b, -inv.m11, &new1, &g);

        let _ = (f, g, tmp0, tmp1);
    }

    #[test]
    fn constant_jump_matrix_apply_cost_probe() {
        // Build actual circuits for constant selected BY matrices to calibrate
        // the row-popcount lower bound. This is still not a full reversible BY
        // update, but it includes the real n-bit add/sub primitive cost and
        // scratch peak for forming the two output rows.
        const WIDTH: usize = 256 + 16 + 2;
        let mut hasher = sha3::Shake128::default();
        hasher.update(b"by-constant-matrix-apply-cost-v1");
        let mut reader = hasher.finalize_xof();
        let mut buf = [0u8; 24];
        let mut total_ccx = 0usize;
        let mut total_terms = 0usize;
        let mut max_peak = 0u32;
        let samples = 24usize;
        for _ in 0..samples {
            reader.read(&mut buf);
            let f_low = (u64::from_le_bytes(buf[0..8].try_into().unwrap()) as i128) | 1;
            let g_low = u64::from_le_bytes(buf[8..16].try_into().unwrap()) as i128;
            let delta = (u64::from_le_bytes(buf[16..24].try_into().unwrap()) % 41) as i64 - 20;
            let (_, _, _, m) = jump_matrix_direct_lowword(16, 16, delta, f_low, g_low);
            let mut b = super::super::B::new();
            let start = b.ops.len();
            emit_constant_matrix_apply_for_cost(&mut b, m, WIDTH);
            let ccx = count_ccx(&b.ops[start..]);
            total_ccx += ccx;
            total_terms += matrix_popcount_adds_i128(m);
            max_peak = max_peak.max(b.peak_qubits);
        }
        let mean_ccx = total_ccx as f64 / samples as f64;
        let mean_terms = total_terms as f64 / samples as f64;
        eprintln!(
            "constant BY w=16 matrix apply cost probe: mean_ccx={mean_ccx:.1}, mean_terms={mean_terms:.2}, ccx_per_term={:.1}, max_peak={max_peak}",
            mean_ccx / mean_terms
        );
        assert!(mean_ccx < 10_000.0, "constant matrix row formation too costly to prototype");
    }

    #[test]
    fn scaled_pair_update_cleanup_cost_probe() {
        // Circuit-level calibration for the reversible replacement step, not
        // just row formation. It forms P·old in width+w bits, interprets the
        // high bits as (P·old)/2^w, then cleans old with the scaled adjugate.
        // This is the core operation a jumped-BY inversion would repeat.
        const WIDTH: usize = 256 + 16 + 2;
        const W: usize = 16;
        let mut hasher = sha3::Shake128::default();
        hasher.update(b"by-scaled-pair-update-cleanup-v1");
        let mut reader = hasher.finalize_xof();
        let mut buf = [0u8; 24];
        let samples = 24usize;
        let mut total_ccx = 0usize;
        let mut max_peak = 0u32;
        for _ in 0..samples {
            reader.read(&mut buf);
            let f_low = (u64::from_le_bytes(buf[0..8].try_into().unwrap()) as i128) | 1;
            let g_low = u64::from_le_bytes(buf[8..16].try_into().unwrap()) as i128;
            let delta = (u64::from_le_bytes(buf[16..24].try_into().unwrap()) % 41) as i64 - 20;
            let (_, _, _, m) = jump_matrix_direct_lowword(W, W, delta, f_low, g_low);
            let mut b = super::super::B::new();
            emit_scaled_pair_update_with_cleanup_for_cost(&mut b, m, WIDTH, W);
            total_ccx += count_ccx(&b.ops);
            max_peak = max_peak.max(b.peak_qubits);
        }
        let mean_ccx = total_ccx as f64 / samples as f64;
        eprintln!(
            "scaled BY w=16 pair update+cleanup probe: mean_ccx={mean_ccx:.1}, max_peak={max_peak}"
        );
        assert!(mean_ccx < 9_000.0, "scaled pair replacement too expensive");
        assert!(max_peak < 1_800, "single-pair replacement peak unexpectedly high");
    }

    fn cadd_qq_fast_for_cost(
        b: &mut super::super::B,
        acc: &[super::super::QubitId],
        a: &[super::super::QubitId],
        ctrl: super::super::QubitId,
    ) {
        let n = acc.len();
        let f = b.alloc_qubits(n);
        for i in 0..n {
            b.ccx(ctrl, a[i], f[i]);
        }
        super::super::add_nbit_qq_fast(b, &f, acc);
        for i in 0..n {
            let m = b.alloc_bit();
            b.hmr(f[i], m);
            b.cz_if(ctrl, a[i], m);
        }
        b.free_vec(&f);
    }

    fn csub_qq_fast_for_cost(
        b: &mut super::super::B,
        acc: &[super::super::QubitId],
        a: &[super::super::QubitId],
        ctrl: super::super::QubitId,
    ) {
        let n = acc.len();
        let f = b.alloc_qubits(n);
        for i in 0..n {
            b.ccx(ctrl, a[i], f[i]);
        }
        super::super::sub_nbit_qq_fast(b, &f, acc);
        for i in 0..n {
            let m = b.alloc_bit();
            b.hmr(f[i], m);
            b.cz_if(ctrl, a[i], m);
        }
        b.free_vec(&f);
    }

    fn inv_odd_mod_pow2_u64(a: u64, w: usize) -> u64 {
        assert!(w > 0 && w <= 63 && (a & 1) == 1);
        let mask = (1u64 << w) - 1;
        let mut x = 1u64;
        // Hensel/Newton doubling; enough rounds for w<=63.
        for _ in 0..6 {
            x = x.wrapping_mul(2u64.wrapping_sub(a.wrapping_mul(x))) & mask;
        }
        x & mask
    }

    #[test]
    fn jump_matrix_depends_on_delta_and_g_over_f_ratio() {
        // BY low-word jumps do not really depend on both low f and low g.
        // Since f is always odd, normalizing by f shows the transition matrix
        // is a function of (delta, h=g/f mod 2^w). Exact enumeration for
        // w<=8 matches the earlier histogram law: distinct matrices = 41*2^w.
        use std::collections::BTreeMap;
        for &w in &[4usize, 6, 8] {
            let mask = (1u64 << w) - 1;
            let mut by_key: BTreeMap<(i64, u64), TransitionMatrix> = BTreeMap::new();
            for delta in -20i64..=20i64 {
                for f_odd in 0..(1usize << (w - 1)) {
                    let f_low = ((f_odd << 1) | 1) as u64;
                    let inv_f = inv_odd_mod_pow2_u64(f_low, w);
                    for g_raw in 0..(1usize << w) {
                        let h = (g_raw as u64).wrapping_mul(inv_f) & mask;
                        let (_, _, _, m) = jump_matrix_direct_lowword(
                            w,
                            w,
                            delta,
                            f_low as i128,
                            g_raw as i128,
                        );
                        match by_key.insert((delta, h), m) {
                            Some(prev) => assert_eq!(prev, m, "matrix not determined by delta,h for w={w}"),
                            None => {}
                        }
                    }
                }
            }
            eprintln!(
                "BY normalized jump keys w={w}: keys={}, expected={}",
                by_key.len(),
                41usize * (1usize << w)
            );
            assert_eq!(by_key.len(), 41usize * (1usize << w));
        }
    }

    fn ratio_window_next_with_pattern(w: usize, delta: i64, h: u64) -> (i64, u64, u64) {
        let signed_h = if (h & (1u64 << (w - 1))) != 0 {
            (h as i128) - (1i128 << w)
        } else {
            h as i128
        };
        let mut d = delta;
        let mut f = SInt::from_u(U256::from(1));
        let mag = U256::from(signed_h.unsigned_abs());
        let mut g = SInt { neg: signed_h < 0, mag };
        let mut pattern = 0u64;
        for i in 0..w {
            if g.bit0() {
                pattern |= 1u64 << i;
            }
            divstep_sint_state(&mut d, &mut f, &mut g);
        }
        let mask = (1u64 << w) - 1;
        let f_low = if f.neg {
            ((!f.mag.as_limbs()[0]).wrapping_add(1)) & mask
        } else {
            f.mag.as_limbs()[0] & mask
        };
        let g_low = if g.neg {
            ((!g.mag.as_limbs()[0]).wrapping_add(1)) & mask
        } else {
            g.mag.as_limbs()[0] & mask
        };
        let inv_f = inv_odd_mod_pow2_u64(f_low, w);
        (d, g_low.wrapping_mul(inv_f) & mask, pattern)
    }

    fn ratio_window_next(w: usize, delta: i64, h: u64) -> (i64, u64) {
        let (d, h_next, _) = ratio_window_next_with_pattern(w, delta, h);
        (d, h_next)
    }

    fn low_ratio_of_sints(f: SInt, g: SInt, w: usize) -> u64 {
        let mask = (1u64 << w) - 1;
        let f_low = if f.neg {
            ((!f.mag.as_limbs()[0]).wrapping_add(1)) & mask
        } else {
            f.mag.as_limbs()[0] & mask
        };
        let g_low = if g.neg {
            ((!g.mag.as_limbs()[0]).wrapping_add(1)) & mask
        } else {
            g.mag.as_limbs()[0] & mask
        };
        g_low.wrapping_mul(inv_odd_mod_pow2_u64(f_low, w)) & mask
    }

    #[test]
    fn pattern_augmented_low_ratio_state_still_not_forward_complete() {
        // Sharper invalidation of the tempting h-only generator.  The current
        // 16 branch bits are determined by (delta, h=g/f mod 2^16), but the
        // *next* window's h is not.  The missing information is high 2-adic
        // data of the denominator pair, not merely the branch pattern.  So a
        // savings-capable generator cannot keep only a 16-bit h register plus
        // the pattern history; it needs a sliding higher-precision state, a
        // rank payload, or a consumed-denominator schedule.
        const W: usize = 16;
        let samples = 2_000usize;
        let mut sampler = Sampler::new(b"by-pattern-augmented-forward-dead-v1", SECP256K1_P);
        let mut windows = 0usize;
        let mut h_next_mismatches = 0usize;
        let mut first: Option<(i64, u64, u64, u64)> = None;
        for _ in 0..samples {
            let x = sampler.next();
            let mut delta = 1i64;
            let mut f = SInt::from_u(SECP256K1_P);
            let mut g = SInt::from_u(x);
            for _ in 0..35 {
                let h = low_ratio_of_sints(f, g, W);
                let (_d_model, h_model, pat) = ratio_window_next_with_pattern(W, delta, h);
                let mut delta_actual = delta;
                let mut f_actual = f;
                let mut g_actual = g;
                for _ in 0..W {
                    divstep_sint_state(&mut delta_actual, &mut f_actual, &mut g_actual);
                }
                let h_actual = low_ratio_of_sints(f_actual, g_actual, W);
                if h_actual != h_model {
                    h_next_mismatches += 1;
                    first.get_or_insert((delta, h, pat, h_actual ^ h_model));
                }
                windows += 1;
                delta = delta_actual;
                f = f_actual;
                g = g_actual;
            }
        }
        let mismatch_rate = h_next_mismatches as f64 / windows as f64;
        eprintln!(
            "BY h16+pattern forward incompleteness: mismatches={h_next_mismatches}/{windows} ({mismatch_rate:.4}), first={first:?}"
        );
        assert!(mismatch_rate > 0.50, "h16+pattern unexpectedly predicts next-window h often enough");
    }

    #[test]
    fn low_ratio_window_state_needs_large_rank_history() {
        // Tempting idea: keep only h=g/f mod 2^w and delta to select BY jump
        // matrices, instead of a full denominator pair. But the h-update is
        // many-to-one. On actual 35-window secp256k1 trajectories, recovering
        // the previous h from (delta',h') needs up to 17 bits of rank per
        // window in this sample. That is about 595 history bits before any
        // arithmetic scratch, so h-only state is not the 600-scratch escape.
        use std::collections::HashMap;
        const W: usize = 16;
        let mut counts: HashMap<(i64, u64), u32> = HashMap::new();
        for delta in -24i64..=64i64 {
            for h in 0u64..(1u64 << W) {
                let out = ratio_window_next(W, delta, h);
                *counts.entry(out).or_insert(0) += 1;
            }
        }

        let samples = 2_000usize;
        let mut sampler = Sampler::new(b"by-low-ratio-rank-history-v1", SECP256K1_P);
        let mut max_rank = 0u32;
        let mut over_16 = 0usize;
        for _ in 0..samples {
            let x = sampler.next();
            let mut delta = 1i64;
            let mut f = SInt::from_u(SECP256K1_P);
            let mut g = SInt::from_u(x);
            let mut sample_max = 0u32;
            for _ in 0..35 {
                assert!((-24..=64).contains(&delta), "delta out of modeled range: {delta}");
                let h = low_ratio_of_sints(f, g, W);
                let out = ratio_window_next(W, delta, h);
                let rank_space = *counts.get(&out).expect("counted output");
                sample_max = sample_max.max(rank_space);
                for _ in 0..W {
                    divstep_sint_state(&mut delta, &mut f, &mut g);
                }
            }
            max_rank = max_rank.max(sample_max);
            if sample_max > (1u32 << 16) {
                over_16 += 1;
            }
        }
        let fail_16 = over_16 as f64 / samples as f64;
        eprintln!(
            "BY low-ratio reversible-state rank: max_rank={max_rank}, bits={}, fail_if_16bit_rank={fail_16:.4}",
            32 - (max_rank - 1).leading_zeros()
        );
        assert!(max_rank > (1u32 << 16), "16-bit/window rank unexpectedly sufficient");
        assert!(fail_16 > 0.01, "16-bit/window rank would meet 1% tolerance; revisit h-only path");
    }

    fn low_signed_sint_for_streaming_test(x: SInt, bits: usize) -> i128 {
        assert!((1..=63).contains(&bits));
        let mask = (1u64 << bits) - 1;
        let low = x.mag.as_limbs()[0] & mask;
        let residue = if x.neg { ((!low).wrapping_add(1)) & mask } else { low };
        if (residue & (1u64 << (bits - 1))) != 0 {
            residue as i128 - (1i128 << bits)
        } else {
            residue as i128
        }
    }

    fn low_signed_sint16_for_streaming_test(x: SInt) -> i128 {
        low_signed_sint_for_streaming_test(x, 16)
    }

    fn u256_limb16_for_streaming_test(x: U256, win: usize) -> u64 {
        let shift = 16 * win;
        if shift >= 256 {
            0
        } else {
            ((x >> shift) & U256::from(0xffffu64)).to::<u64>()
        }
    }

    fn mask_u512_bits_for_streaming_test(bits: usize) -> U512 {
        if bits == 512 {
            !U512::ZERO
        } else {
            (U512::from(1u64) << bits) - U512::from(1u64)
        }
    }

    fn inv_odd_u512_pow2_for_streaming_test(a: U512, bits: usize) -> U512 {
        let modulus = U512::from(1u64) << bits;
        let mask = modulus - U512::from(1u64);
        let mut x = U512::from(1u64);
        // Newton iteration doubles the number of correct low bits each round.
        for _ in 0..10 {
            let ax = (a * x) & mask;
            let two_minus_ax = (U512::from(2u64) + modulus - ax) & mask;
            x = (x * two_minus_ax) & mask;
        }
        x
    }

    fn sign_extend_u512_for_streaming_test(x: U512, from_bits: usize, to_bits: usize) -> U512 {
        debug_assert!(from_bits <= to_bits);
        let mut y = x & mask_u512_bits_for_streaming_test(from_bits);
        if from_bits > 0 && y.bit(from_bits - 1) {
            for i in from_bits..to_bits {
                y.set_bit(i, true);
            }
        }
        y & mask_u512_bits_for_streaming_test(to_bits)
    }

    fn mul_i128_mod_width_for_streaming_test(x: U512, coeff: i128, bits: usize) -> U512 {
        let c = signed_coeff_mod_width_for_test(coeff, bits);
        (x * c) & mask_u512_bits_for_streaming_test(bits)
    }

    fn add_i128_term_mod_width_for_streaming_test(acc: U512, x: U512, coeff: i128, bits: usize) -> U512 {
        (acc + mul_i128_mod_width_for_streaming_test(x, coeff, bits)) & mask_u512_bits_for_streaming_test(bits)
    }

    type Big9 = [u64; 9];

    fn big9_mask(bits: usize) -> Big9 {
        let mut out = [0u64; 9];
        let full = bits / 64;
        let rem = bits % 64;
        for limb in out.iter_mut().take(full.min(9)) {
            *limb = u64::MAX;
        }
        if full < 9 && rem != 0 {
            out[full] = (1u64 << rem) - 1;
        }
        out
    }

    fn big9_apply_mask(x: &mut Big9, bits: usize) {
        let mask = big9_mask(bits);
        for i in 0..9 {
            x[i] &= mask[i];
        }
    }

    fn big9_from_u256(x: U256) -> Big9 {
        let mut out = [0u64; 9];
        let limbs = x.as_limbs();
        out[..4].copy_from_slice(limbs);
        out
    }

    fn big9_bit(x: &Big9, bit: usize) -> bool {
        if bit >= 576 {
            false
        } else {
            ((x[bit / 64] >> (bit % 64)) & 1) != 0
        }
    }

    fn big9_set_bit(x: &mut Big9, bit: usize) {
        x[bit / 64] |= 1u64 << (bit % 64);
    }

    fn big9_add_word(tmp: &mut [u64; 19], mut idx: usize, mut val: u64) {
        while val != 0 {
            let (sum, carry) = tmp[idx].overflowing_add(val);
            tmp[idx] = sum;
            val = carry as u64;
            idx += 1;
        }
    }

    fn big9_mul_mod(a: &Big9, b: &Big9, bits: usize) -> Big9 {
        let mut tmp = [0u64; 19];
        let need_limbs = ((bits + 63) / 64).min(9);
        for i in 0..need_limbs {
            for j in 0..need_limbs {
                if i + j >= need_limbs {
                    continue;
                }
                let prod = (a[i] as u128) * (b[j] as u128);
                big9_add_word(&mut tmp, i + j, prod as u64);
                big9_add_word(&mut tmp, i + j + 1, (prod >> 64) as u64);
            }
        }
        let mut out = [0u64; 9];
        out.copy_from_slice(&tmp[..9]);
        big9_apply_mask(&mut out, bits);
        out
    }

    fn big9_sub1(mut x: Big9, bits: usize) -> Big9 {
        let mut i = 0usize;
        loop {
            let (v, borrow) = x[i].overflowing_sub(1);
            x[i] = v;
            if !borrow {
                break;
            }
            i += 1;
        }
        big9_apply_mask(&mut x, bits);
        x
    }

    fn big9_add1(mut x: Big9, bits: usize) -> Big9 {
        let mut i = 0usize;
        loop {
            let (v, carry) = x[i].overflowing_add(1);
            x[i] = v;
            if !carry {
                break;
            }
            i += 1;
        }
        big9_apply_mask(&mut x, bits);
        x
    }

    fn big9_shr1(mut x: Big9, bits: usize) -> Big9 {
        let mut carry = 0u64;
        for limb in x.iter_mut().rev() {
            let next = *limb << 63;
            *limb = (*limb >> 1) | carry;
            carry = next;
        }
        big9_apply_mask(&mut x, bits.saturating_sub(1));
        x
    }

    fn big9_inv_odd_mod_pow2(a: &Big9, bits: usize) -> Big9 {
        if bits == 0 {
            return [0u64; 9];
        }
        assert!(big9_bit(a, 0));
        let mut inv = [0u64; 9];
        inv[0] = 1;
        for i in 1..bits {
            let prod = big9_mul_mod(a, &inv, i + 1);
            if big9_bit(&prod, i) {
                big9_set_bit(&mut inv, i);
            }
        }
        big9_apply_mask(&mut inv, bits);
        inv
    }

    fn h_ratio_step_big9_for_streaming_test(delta: i64, h: Big9, t: usize) -> (i64, Big9, bool) {
        let odd = big9_bit(&h, 0);
        if t == 1 {
            let next_delta = if delta > 0 && odd { 1 - delta } else { 1 + delta };
            return (next_delta, [0u64; 9], odd);
        }
        let next_bits = t - 1;
        if delta > 0 && odd {
            let half_num = big9_shr1(big9_sub1(h, t), t);
            let inv_h = big9_inv_odd_mod_pow2(&h, next_bits);
            (1 - delta, big9_mul_mod(&half_num, &inv_h, next_bits), odd)
        } else if odd {
            (1 + delta, big9_shr1(big9_add1(h, t), t), odd)
        } else {
            (1 + delta, big9_shr1(h, t), odd)
        }
    }

    fn sint_residue_u512_for_streaming_test(x: SInt, bits: usize) -> U512 {
        let mask = mask_u512_bits_for_streaming_test(bits);
        let mag = U512::from(x.mag) & mask;
        if x.neg {
            ((!mag) + U512::from(1u64)) & mask
        } else {
            mag
        }
    }

    fn h_ratio_step_u512_for_streaming_test(delta: i64, h: U512, t: usize) -> (i64, U512, bool) {
        let odd = h.bit(0);
        if t == 1 {
            let next_delta = if delta > 0 && odd { 1 - delta } else { 1 + delta };
            return (next_delta, U512::ZERO, odd);
        }
        let next_bits = t - 1;
        let mask = mask_u512_bits_for_streaming_test(next_bits);
        if delta > 0 && odd {
            let half_num = (h - U512::from(1u64)) >> 1usize;
            let inv_h = inv_odd_u512_pow2_for_streaming_test(h, next_bits);
            (1 - delta, (half_num * inv_h) & mask, odd)
        } else if odd {
            (1 + delta, ((h + U512::from(1u64)) >> 1usize) & mask, odd)
        } else {
            (1 + delta, (h >> 1usize) & mask, odd)
        }
    }

    fn low_i16_from_residue_for_streaming_test(x: U512) -> i128 {
        let low = x.as_limbs()[0] & 0xffff;
        if (low & 0x8000) != 0 {
            low as i128 - (1i128 << 16)
        } else {
            low as i128
        }
    }

    fn streaming_selector_model_matches_for_test(x: U256, state_limbs: usize) -> bool {
        const W: usize = 16;
        const WINDOWS: usize = 35;
        let bits = state_limbs * W;
        let prod_bits = bits + W;
        let mut a = [[U512::ZERO; 2]; 2];
        a[0][0] = U512::from(1u64);
        a[1][1] = U512::from(1u64);
        let mut c = [U512::ZERO; 2];
        let p = SECP256K1_P;
        let mut delta = 1i64;
        let mut actual_delta = 1i64;
        let mut f = SInt::from_u(p);
        let mut g = SInt::from_u(x);
        for win in 0..WINDOWS {
            let limb_f = U512::from(u256_limb16_for_streaming_test(p, win));
            let limb_g = U512::from(u256_limb16_for_streaming_test(x, win));
            let v0_bits = (a[0][0] * limb_f + a[0][1] * limb_g + c[0])
                & mask_u512_bits_for_streaming_test(bits);
            let v1_bits = (a[1][0] * limb_f + a[1][1] * limb_g + c[1])
                & mask_u512_bits_for_streaming_test(bits);
            let low0 = low_i16_from_residue_for_streaming_test(v0_bits);
            let low1 = low_i16_from_residue_for_streaming_test(v1_bits);
            let exp0 = low_signed_sint16_for_streaming_test(f);
            let exp1 = low_signed_sint16_for_streaming_test(g);
            if low0 != exp0 || low1 != exp1 {
                return false;
            }
            let bits_vec = branch_bits_for_lowword_window(W, delta, low0, low1);
            let m = matrix_from_branch_bits(delta, &bits_vec);

            let v0_ext = sign_extend_u512_for_streaming_test(v0_bits, bits, prod_bits);
            let v1_ext = sign_extend_u512_for_streaming_test(v1_bits, bits, prod_bits);
            let mut n0 = U512::ZERO;
            n0 = add_i128_term_mod_width_for_streaming_test(n0, v0_ext, m.m00, prod_bits);
            n0 = add_i128_term_mod_width_for_streaming_test(n0, v1_ext, m.m01, prod_bits);
            let mut n1 = U512::ZERO;
            n1 = add_i128_term_mod_width_for_streaming_test(n1, v0_ext, m.m10, prod_bits);
            n1 = add_i128_term_mod_width_for_streaming_test(n1, v1_ext, m.m11, prod_bits);
            // The low-word window matrix guarantees divisibility by 2^W.
            if (n0.as_limbs()[0] & 0xffff) != 0 || (n1.as_limbs()[0] & 0xffff) != 0 {
                return false;
            }
            c[0] = sign_extend_u512_for_streaming_test(n0 >> W, bits, bits);
            c[1] = sign_extend_u512_for_streaming_test(n1 >> W, bits, bits);

            let old_a = a;
            a[0][0] = U512::ZERO;
            a[0][0] = add_i128_term_mod_width_for_streaming_test(a[0][0], old_a[0][0], m.m00, bits);
            a[0][0] = add_i128_term_mod_width_for_streaming_test(a[0][0], old_a[1][0], m.m01, bits);
            a[0][1] = U512::ZERO;
            a[0][1] = add_i128_term_mod_width_for_streaming_test(a[0][1], old_a[0][1], m.m00, bits);
            a[0][1] = add_i128_term_mod_width_for_streaming_test(a[0][1], old_a[1][1], m.m01, bits);
            a[1][0] = U512::ZERO;
            a[1][0] = add_i128_term_mod_width_for_streaming_test(a[1][0], old_a[0][0], m.m10, bits);
            a[1][0] = add_i128_term_mod_width_for_streaming_test(a[1][0], old_a[1][0], m.m11, bits);
            a[1][1] = U512::ZERO;
            a[1][1] = add_i128_term_mod_width_for_streaming_test(a[1][1], old_a[0][1], m.m10, bits);
            a[1][1] = add_i128_term_mod_width_for_streaming_test(a[1][1], old_a[1][1], m.m11, bits);
            delta = m.delta_final;

            for _ in 0..W {
                divstep_sint_state(&mut actual_delta, &mut f, &mut g);
            }
            if actual_delta != delta {
                return false;
            }
        }
        true
    }

    fn streaming_selector_constant_folded_separate_matches_for_test(
        x: U256,
        b_limbs: usize,
        c0_limbs: usize,
        c1_limbs: usize,
    ) -> bool {
        const W: usize = 16;
        const WINDOWS: usize = 35;
        let b_bits = b_limbs * W;
        let c0_bits = c0_limbs * W;
        let c1_bits = c1_limbs * W;
        let prod_bits = b_bits.max(c0_bits).max(c1_bits) + W;
        let p = SECP256K1_P;
        let mut bx = [U512::ZERO, U512::from(1u64)];
        let mut c = [U512::from(p) & mask_u512_bits_for_streaming_test(c0_bits), U512::ZERO];
        let mut delta = 1i64;
        let mut actual_delta = 1i64;
        let mut f = SInt::from_u(p);
        let mut g = SInt::from_u(x);
        for win in 0..WINDOWS {
            let limb_x = U512::from(u256_limb16_for_streaming_test(x, win));
            let b0_ext = sign_extend_u512_for_streaming_test(bx[0], b_bits, prod_bits);
            let b1_ext = sign_extend_u512_for_streaming_test(bx[1], b_bits, prod_bits);
            let c0_ext = sign_extend_u512_for_streaming_test(c[0], c0_bits, prod_bits);
            let c1_ext = sign_extend_u512_for_streaming_test(c[1], c1_bits, prod_bits);
            let v0_ext = (b0_ext * limb_x + c0_ext) & mask_u512_bits_for_streaming_test(prod_bits);
            let v1_ext = (b1_ext * limb_x + c1_ext) & mask_u512_bits_for_streaming_test(prod_bits);
            let low0 = low_i16_from_residue_for_streaming_test(v0_ext);
            let low1 = low_i16_from_residue_for_streaming_test(v1_ext);
            if low0 != low_signed_sint16_for_streaming_test(f)
                || low1 != low_signed_sint16_for_streaming_test(g)
            {
                return false;
            }
            let bits_vec = branch_bits_for_lowword_window(W, delta, low0, low1);
            let m = matrix_from_branch_bits(delta, &bits_vec);
            let mut n0 = U512::ZERO;
            n0 = add_i128_term_mod_width_for_streaming_test(n0, v0_ext, m.m00, prod_bits);
            n0 = add_i128_term_mod_width_for_streaming_test(n0, v1_ext, m.m01, prod_bits);
            let mut n1 = U512::ZERO;
            n1 = add_i128_term_mod_width_for_streaming_test(n1, v0_ext, m.m10, prod_bits);
            n1 = add_i128_term_mod_width_for_streaming_test(n1, v1_ext, m.m11, prod_bits);
            if (n0.as_limbs()[0] & 0xffff) != 0 || (n1.as_limbs()[0] & 0xffff) != 0 {
                return false;
            }
            c[0] = sign_extend_u512_for_streaming_test(n0 >> W, prod_bits - W, c0_bits);
            c[1] = sign_extend_u512_for_streaming_test(n1 >> W, prod_bits - W, c1_bits);

            let old_bx = bx;
            bx[0] = U512::ZERO;
            bx[0] = add_i128_term_mod_width_for_streaming_test(bx[0], old_bx[0], m.m00, b_bits);
            bx[0] = add_i128_term_mod_width_for_streaming_test(bx[0], old_bx[1], m.m01, b_bits);
            bx[1] = U512::ZERO;
            bx[1] = add_i128_term_mod_width_for_streaming_test(bx[1], old_bx[0], m.m10, b_bits);
            bx[1] = add_i128_term_mod_width_for_streaming_test(bx[1], old_bx[1], m.m11, b_bits);
            delta = m.delta_final;

            for _ in 0..W {
                divstep_sint_state(&mut actual_delta, &mut f, &mut g);
            }
            if actual_delta != delta {
                return false;
            }
        }
        true
    }

    fn streaming_selector_tail_c_only_matches_for_test(x: U256, b_limbs: usize, c0_limbs: usize, c1_limbs: usize) -> bool {
        const W: usize = 16;
        let b_bits = b_limbs * W;
        let c0_bits = c0_limbs * W;
        let c1_bits = c1_limbs * W;
        let prod_bits = b_bits.max(c0_bits).max(c1_bits) + W;
        let p = SECP256K1_P;
        let mut bx = [U512::ZERO, U512::from(1u64)];
        let mut c = [U512::from(p) & mask_u512_bits_for_streaming_test(c0_bits), U512::ZERO];
        let mut delta = 1i64;
        let mut actual_delta = 1i64;
        let mut f = SInt::from_u(p);
        let mut g = SInt::from_u(x);
        for win in 0..35 {
            let limb_x = U512::from(u256_limb16_for_streaming_test(x, win));
            let use_b = win < 16;
            let b0_ext = if use_b { sign_extend_u512_for_streaming_test(bx[0], b_bits, prod_bits) } else { U512::ZERO };
            let b1_ext = if use_b { sign_extend_u512_for_streaming_test(bx[1], b_bits, prod_bits) } else { U512::ZERO };
            let c0_ext = sign_extend_u512_for_streaming_test(c[0], c0_bits, prod_bits);
            let c1_ext = sign_extend_u512_for_streaming_test(c[1], c1_bits, prod_bits);
            let v0_ext = (b0_ext * limb_x + c0_ext) & mask_u512_bits_for_streaming_test(prod_bits);
            let v1_ext = (b1_ext * limb_x + c1_ext) & mask_u512_bits_for_streaming_test(prod_bits);
            let low0 = low_i16_from_residue_for_streaming_test(v0_ext);
            let low1 = low_i16_from_residue_for_streaming_test(v1_ext);
            if low0 != low_signed_sint16_for_streaming_test(f)
                || low1 != low_signed_sint16_for_streaming_test(g)
            {
                return false;
            }
            let bits_vec = branch_bits_for_lowword_window(W, delta, low0, low1);
            let m = matrix_from_branch_bits(delta, &bits_vec);
            let mut n0 = U512::ZERO;
            n0 = add_i128_term_mod_width_for_streaming_test(n0, v0_ext, m.m00, prod_bits);
            n0 = add_i128_term_mod_width_for_streaming_test(n0, v1_ext, m.m01, prod_bits);
            let mut n1 = U512::ZERO;
            n1 = add_i128_term_mod_width_for_streaming_test(n1, v0_ext, m.m10, prod_bits);
            n1 = add_i128_term_mod_width_for_streaming_test(n1, v1_ext, m.m11, prod_bits);
            if (n0.as_limbs()[0] & 0xffff) != 0 || (n1.as_limbs()[0] & 0xffff) != 0 {
                return false;
            }
            c[0] = sign_extend_u512_for_streaming_test(n0 >> W, prod_bits - W, c0_bits);
            c[1] = sign_extend_u512_for_streaming_test(n1 >> W, prod_bits - W, c1_bits);
            if use_b {
                let old_bx = bx;
                bx[0] = U512::ZERO;
                bx[0] = add_i128_term_mod_width_for_streaming_test(bx[0], old_bx[0], m.m00, b_bits);
                bx[0] = add_i128_term_mod_width_for_streaming_test(bx[0], old_bx[1], m.m01, b_bits);
                bx[1] = U512::ZERO;
                bx[1] = add_i128_term_mod_width_for_streaming_test(bx[1], old_bx[0], m.m10, b_bits);
                bx[1] = add_i128_term_mod_width_for_streaming_test(bx[1], old_bx[1], m.m11, b_bits);
            }
            delta = m.delta_final;
            for _ in 0..W {
                divstep_sint_state(&mut actual_delta, &mut f, &mut g);
            }
            if actual_delta != delta {
                return false;
            }
        }
        true
    }

    fn streaming_selector_constant_folded_matches_for_test(x: U256, state_limbs: usize) -> bool {
        const W: usize = 16;
        const WINDOWS: usize = 35;
        let bits = state_limbs * W;
        let prod_bits = bits + W;
        let mut bx = [U512::ZERO; 2];
        bx[1] = U512::from(1u64);
        let p = SECP256K1_P;
        let mut c = [U512::from(p) & mask_u512_bits_for_streaming_test(bits), U512::ZERO];
        let mut delta = 1i64;
        let mut actual_delta = 1i64;
        let mut f = SInt::from_u(p);
        let mut g = SInt::from_u(x);
        for win in 0..WINDOWS {
            let limb_x = U512::from(u256_limb16_for_streaming_test(x, win));
            let v0_bits = (bx[0] * limb_x + c[0]) & mask_u512_bits_for_streaming_test(bits);
            let v1_bits = (bx[1] * limb_x + c[1]) & mask_u512_bits_for_streaming_test(bits);
            let low0 = low_i16_from_residue_for_streaming_test(v0_bits);
            let low1 = low_i16_from_residue_for_streaming_test(v1_bits);
            if low0 != low_signed_sint16_for_streaming_test(f)
                || low1 != low_signed_sint16_for_streaming_test(g)
            {
                return false;
            }
            let bits_vec = branch_bits_for_lowword_window(W, delta, low0, low1);
            let m = matrix_from_branch_bits(delta, &bits_vec);
            let v0_ext = sign_extend_u512_for_streaming_test(v0_bits, bits, prod_bits);
            let v1_ext = sign_extend_u512_for_streaming_test(v1_bits, bits, prod_bits);
            let mut n0 = U512::ZERO;
            n0 = add_i128_term_mod_width_for_streaming_test(n0, v0_ext, m.m00, prod_bits);
            n0 = add_i128_term_mod_width_for_streaming_test(n0, v1_ext, m.m01, prod_bits);
            let mut n1 = U512::ZERO;
            n1 = add_i128_term_mod_width_for_streaming_test(n1, v0_ext, m.m10, prod_bits);
            n1 = add_i128_term_mod_width_for_streaming_test(n1, v1_ext, m.m11, prod_bits);
            if (n0.as_limbs()[0] & 0xffff) != 0 || (n1.as_limbs()[0] & 0xffff) != 0 {
                return false;
            }
            c[0] = sign_extend_u512_for_streaming_test(n0 >> W, bits, bits);
            c[1] = sign_extend_u512_for_streaming_test(n1 >> W, bits, bits);

            let old_bx = bx;
            bx[0] = U512::ZERO;
            bx[0] = add_i128_term_mod_width_for_streaming_test(bx[0], old_bx[0], m.m00, bits);
            bx[0] = add_i128_term_mod_width_for_streaming_test(bx[0], old_bx[1], m.m01, bits);
            bx[1] = U512::ZERO;
            bx[1] = add_i128_term_mod_width_for_streaming_test(bx[1], old_bx[0], m.m10, bits);
            bx[1] = add_i128_term_mod_width_for_streaming_test(bx[1], old_bx[1], m.m11, bits);
            delta = m.delta_final;

            for _ in 0..W {
                divstep_sint_state(&mut actual_delta, &mut f, &mut g);
            }
            if actual_delta != delta {
                return false;
            }
        }
        true
    }

    #[test]
    fn streaming_limb_selector_is_exact_but_state_heavy() {
        // Architectural pivot test for a real BY-DIV selector generator.  The
        // h-only state is not forward-complete, but a streaming base-2^16 limb
        // recurrence is: keep the current affine carry
        //
        //     (f_j,g_j) = A_j · (p>>16j, x>>16j) + c_j
        //
        // modulo a fixed number of 16-bit limbs.  This generates all 35 branch
        // windows without a full 560-bit denominator pair.  The catch is state:
        // the first robust sampled setting needs roughly A(4 entries)+c(2
        // entries) at 12 limbs each = 1152 logical bits before compression.  This
        // is a plausible Toffoli architecture, but not yet a 1175q Google-shape
        // circuit unless A is recomputed/compressed from pattern history.
        let samples = 64usize;
        let mut sampler = Sampler::new(b"by-streaming-limb-selector-v1", SECP256K1_P);
        let mut k8_failures = 0usize;
        for _ in 0..samples {
            let x = sampler.next();
            if !streaming_selector_model_matches_for_test(x, 8) {
                k8_failures += 1;
            }
            assert!(
                streaming_selector_model_matches_for_test(x, 12),
                "12-limb streaming selector lost exactness"
            );
        }
        let state_bits = 6 * 12 * 16;
        eprintln!(
            "BY streaming limb selector: samples={samples}, k8_failures={k8_failures}, k12_state_bits={state_bits}"
        );
        assert!(k8_failures > 0, "8-limb streaming state unexpectedly exact on all samples");
        assert!(state_bits > 600, "streaming selector would already fit the low-scratch target");
    }

    fn streaming_selector_projective_normalized_matches_for_test(x: U256, state_limbs: usize) -> bool {
        const W: usize = 16;
        const WINDOWS: usize = 35;
        let bits = state_limbs * W;
        let prod_bits = bits + W;
        let p = SECP256K1_P;
        let inv_p = inv_odd_u512_pow2_for_streaming_test(U512::from(p), bits);
        // Scaled state: f = 1 + b0*x_tail, g = c1 + b1*x_tail.
        // The missing common odd scale is intentionally not stored.
        let mut b0 = U512::ZERO;
        let mut b1 = inv_p;
        let mut c1 = U512::ZERO;
        let mut delta = 1i64;
        let mut actual_delta = 1i64;
        let mut f = SInt::from_u(p);
        let mut g = SInt::from_u(x);
        for win in 0..WINDOWS {
            let limb_x = U512::from(u256_limb16_for_streaming_test(x, win));
            let v0_bits = (U512::from(1u64) + b0 * limb_x) & mask_u512_bits_for_streaming_test(bits);
            let v1_bits = (c1 + b1 * limb_x) & mask_u512_bits_for_streaming_test(bits);
            let low0 = low_i16_from_residue_for_streaming_test(v0_bits);
            let low1 = low_i16_from_residue_for_streaming_test(v1_bits);
            let predicted_bits = branch_bits_for_lowword_window(W, delta, low0, low1);

            let mut d_check = actual_delta;
            let mut f_check = f;
            let mut g_check = g;
            for &pred_odd in &predicted_bits {
                if pred_odd != g_check.bit0() {
                    return false;
                }
                divstep_sint_state(&mut d_check, &mut f_check, &mut g_check);
            }
            let m = matrix_from_branch_bits(delta, &predicted_bits);

            let v0_ext = sign_extend_u512_for_streaming_test(v0_bits, bits, prod_bits);
            let v1_ext = sign_extend_u512_for_streaming_test(v1_bits, bits, prod_bits);
            let mut n0 = U512::ZERO;
            n0 = add_i128_term_mod_width_for_streaming_test(n0, v0_ext, m.m00, prod_bits);
            n0 = add_i128_term_mod_width_for_streaming_test(n0, v1_ext, m.m01, prod_bits);
            let mut n1 = U512::ZERO;
            n1 = add_i128_term_mod_width_for_streaming_test(n1, v0_ext, m.m10, prod_bits);
            n1 = add_i128_term_mod_width_for_streaming_test(n1, v1_ext, m.m11, prod_bits);
            if (n0.as_limbs()[0] & 0xffff) != 0 || (n1.as_limbs()[0] & 0xffff) != 0 {
                return false;
            }
            let cc0 = sign_extend_u512_for_streaming_test(n0 >> W, bits, bits);
            let cc1 = sign_extend_u512_for_streaming_test(n1 >> W, bits, bits);
            if !cc0.bit(0) {
                return false;
            }
            let mut nb0 = U512::ZERO;
            nb0 = add_i128_term_mod_width_for_streaming_test(nb0, b0, m.m00, bits);
            nb0 = add_i128_term_mod_width_for_streaming_test(nb0, b1, m.m01, bits);
            let mut nb1 = U512::ZERO;
            nb1 = add_i128_term_mod_width_for_streaming_test(nb1, b0, m.m10, bits);
            nb1 = add_i128_term_mod_width_for_streaming_test(nb1, b1, m.m11, bits);
            let inv_cc0 = inv_odd_u512_pow2_for_streaming_test(cc0, bits);
            b0 = (nb0 * inv_cc0) & mask_u512_bits_for_streaming_test(bits);
            b1 = (nb1 * inv_cc0) & mask_u512_bits_for_streaming_test(bits);
            c1 = (cc1 * inv_cc0) & mask_u512_bits_for_streaming_test(bits);
            delta = m.delta_final;

            for _ in 0..W {
                divstep_sint_state(&mut actual_delta, &mut f, &mut g);
            }
            if actual_delta != delta {
                return false;
            }
        }
        true
    }

    #[test]
    fn streaming_limb_selector_folds_constant_p_column_but_still_too_large() {
        // First compression of the streaming selector: the first denominator
        // column multiplies the constant p, so it can be folded into the carry
        // vector.  The state drops from six entries (full A plus c) to four
        // entries (two x-column coefficients plus two carries).  Exactness now
        // needs 17 16-bit limbs because p itself is a 256-bit positive constant;
        // 16 limbs aliases it as a negative two's-complement value.
        let samples = 64usize;
        let mut sampler = Sampler::new(b"by-streaming-limb-selector-fold-p-v1", SECP256K1_P);
        let mut k16_failures = 0usize;
        for _ in 0..samples {
            let x = sampler.next();
            if !streaming_selector_constant_folded_matches_for_test(x, 16) {
                k16_failures += 1;
            }
            assert!(
                streaming_selector_constant_folded_matches_for_test(x, 17),
                "17-limb constant-folded streaming selector lost exactness"
            );
        }
        let state_bits = 4 * 17 * 16;
        eprintln!(
            "BY streaming selector with constant-p fold: samples={samples}, k16_failures={k16_failures}, k17_state_bits={state_bits}"
        );
        assert!(k16_failures > 0, "16-limb folded selector unexpectedly exact on all samples");
        assert!(state_bits < 1152, "constant-p fold did not improve the naive A/c state");
        assert!(state_bits > 600, "folded selector would already fit the low-scratch target");
    }

    #[test]
    fn separate_width_streaming_selector_reaches_816_bits() {
        // The previous folded selector used one width for all four entries.
        // In fact the x-column coefficients only need 9 limbs, while the carry
        // rows need 17 and 16 limbs respectively.  This is the first exact
        // selector state below 1k bits.
        let samples = 64usize;
        let mut sampler = Sampler::new(b"by-streaming-selector-separate-width-v1", SECP256K1_P);
        let mut b8_failures = 0usize;
        let mut c1_15_failures = 0usize;
        for _ in 0..samples {
            let x = sampler.next();
            if !streaming_selector_constant_folded_separate_matches_for_test(x, 8, 17, 16) {
                b8_failures += 1;
            }
            if !streaming_selector_constant_folded_separate_matches_for_test(x, 9, 17, 15) {
                c1_15_failures += 1;
            }
            assert!(
                streaming_selector_constant_folded_separate_matches_for_test(x, 9, 17, 16),
                "separate-width folded streaming selector lost exactness"
            );
        }
        let state_bits = (2 * 9 + 17 + 16) * 16;
        eprintln!(
            "BY separate-width streaming selector: samples={samples}, b8_failures={b8_failures}, c1_15_failures={c1_15_failures}, state_bits={state_bits}"
        );
        assert!(b8_failures > 0, "8-limb x-column unexpectedly exact");
        assert!(c1_15_failures > 0, "15-limb c1 unexpectedly exact");
        assert!(state_bits == 816, "unexpected separate-width state size");
    }

    #[test]
    fn tail_exhausted_streaming_selector_has_528_bit_carry_core() {
        // After 16 windows all bits of the 256-bit quantum denominator x have
        // been consumed.  From that point on the x-column coefficients no
        // longer contribute to branch selection; only the two carry rows remain.
        // This identifies a 528-bit rolling selector core for windows 16..34.
        let samples = 64usize;
        let mut sampler = Sampler::new(b"by-streaming-selector-tail-c-only-v1", SECP256K1_P);
        for _ in 0..samples {
            let x = sampler.next();
            assert!(
                streaming_selector_tail_c_only_matches_for_test(x, 9, 17, 16),
                "tail c-only selector did not match full divsteps"
            );
        }
        let c_core_bits = (17 + 16) * 16;
        eprintln!(
            "BY tail-exhausted streaming selector: samples={samples}, c_core_bits={c_core_bits}, b_workspace_bits={}",
            2 * 9 * 16
        );
        assert_eq!(c_core_bits, 528);
        assert!(c_core_bits < 600, "post-tail selector core misses low-scratch target");
    }

    fn inv_odd_u32_pow2_for_ratio_test(a: u32, bits: usize) -> u32 {
        if bits == 0 {
            return 0;
        }
        assert_eq!(a & 1, 1);
        let mut inv = 1u32;
        for i in 1..bits {
            if (((a.wrapping_mul(inv)) >> i) & 1) != 0 {
                inv |= 1u32 << i;
            }
        }
        inv & ((1u32 << bits) - 1)
    }

    fn ratio_a_step_anf_stats_for_test(n: usize, out_bit: usize) -> (usize, usize) {
        let vars = n - 1; // h is odd: h = 1 + 2*z, z has n-1 bits.
        let size = 1usize << vars;
        let mask = (1u32 << (n - 1)) - 1;
        let mut anf = vec![0u8; size];
        for (z, slot) in anf.iter_mut().enumerate() {
            let h = 1u32 + ((z as u32) << 1);
            let half_num = (h - 1) >> 1;
            let inv_h = inv_odd_u32_pow2_for_ratio_test(h, n - 1);
            let h_prime = half_num.wrapping_mul(inv_h) & mask;
            *slot = ((h_prime >> out_bit) & 1) as u8;
        }
        for bit in 0..vars {
            for idx in 0..size {
                if (idx & (1usize << bit)) != 0 {
                    anf[idx] ^= anf[idx ^ (1usize << bit)];
                }
            }
        }
        let density = anf.iter().filter(|&&c| c != 0).count();
        let degree = anf
            .iter()
            .enumerate()
            .filter_map(|(i, &c)| if c != 0 { Some(i.count_ones() as usize) } else { None })
            .max()
            .unwrap_or(0);
        (degree, density)
    }

    fn v2_i128_for_ratio_test(x: i128) -> usize {
        if x == 0 {
            128
        } else {
            (x.unsigned_abs()).trailing_zeros() as usize
        }
    }

    fn ratio_window_v2_hist_for_test(w: usize, samples: usize) -> Vec<usize> {
        let mut sampler = Sampler::new(b"by-ratio-window-size-v2-v1", SECP256K1_P);
        let mut hist = vec![0usize; w + 1];
        for _ in 0..samples {
            let x = sampler.next();
            let mut delta = 1i64;
            let mut f = SInt::from_u(SECP256K1_P);
            let mut g = SInt::from_u(x);
            let mut done = 0usize;
            while done + w <= 35 * 16 {
                let f_low = low_signed_sint_for_streaming_test(f, w);
                let g_low = low_signed_sint_for_streaming_test(g, w);
                let bits = branch_bits_for_lowword_window(w, delta, f_low, g_low);
                let m = matrix_from_branch_bits(delta, &bits);
                let v = v2_i128_for_ratio_test(m.m01).min(w);
                hist[v] += 1;
                for _ in 0..w {
                    divstep_sint_state(&mut delta, &mut f, &mut g);
                }
                done += w;
            }
        }
        hist
    }

    #[test]
    fn wider_ratio_windows_do_not_remove_mobius_inverse_problem() {
        // Larger windows reduce the number of Möbius updates, but make the
        // decoder/table larger.  They would be attractive if m01 gained enough
        // powers of two to make D=q0+m01*H almost constant.  Sampling W=32 still
        // leaves many low-valuation denominators.
        let hist32 = ratio_window_v2_hist_for_test(32, 64);
        let windows32: usize = hist32.iter().sum();
        let weak32: usize = hist32.iter().take(5).sum();
        eprintln!(
            "BY ratio W=32 denominator v2(m01): windows={windows32}, hist={hist32:?}, v<=4={weak32}"
        );
        assert_eq!(windows32, 64 * 17);
        assert!(weak32 * 3 > windows32, "W=32 unexpectedly makes denominators almost constant");
    }

    #[test]
    fn full_ratio_initial_constant_multiply_is_not_the_main_blocker() {
        // h0 = x * p^-1 mod 2^560 is a rectangular constant multiply from a
        // 256-bit input.  A naive controlled-add implementation needs one
        // shifted 560-bit-ish add per x bit; this is large but still in the
        // few-10^5 Toffoli range for compute+uncompute, unlike per-A variable
        // inversions.
        const TOTAL_BITS: usize = 35 * 16;
        const X_BITS: usize = 256;
        let p_big = big9_from_u256(SECP256K1_P);
        let p_inv = big9_inv_odd_mod_pow2(&p_big, TOTAL_BITS);
        let inv_pop = (0..TOTAL_BITS).filter(|&i| big9_bit(&p_inv, i)).count();
        let add_widths: usize = (0..X_BITS).map(|i| TOTAL_BITS - i).sum();
        let one_toffoli_per_bit_roundtrip = 2 * add_widths;
        let two_toffoli_per_bit_roundtrip = 4 * add_widths;
        eprintln!(
            "BY full-ratio h0 constmul budget: p_inv_pop={inv_pop}, add_widths={add_widths}, roundtrip_1t={one_toffoli_per_bit_roundtrip}, roundtrip_2t={two_toffoli_per_bit_roundtrip}"
        );
        assert_eq!(add_widths, 110_720);
        assert!(one_toffoli_per_bit_roundtrip < 250_000);
        assert!(two_toffoli_per_bit_roundtrip < 500_000);
    }

    #[test]
    fn ratio_window_mobius_denominators_are_not_near_constant() {
        // Windowing the ratio update gives
        //   h' = (m10 + m11*h) / (m00 + m01*h)
        // and, after consuming the low 16 bits, an odd denominator
        //   D = q0 + m01*H.
        // If m01 had high 2-adic valuation in most windows, D^{-1} could be a
        // short series.  In real traces m01 is often odd or only weakly even,
        // so windowing does not by itself remove the variable-inverse problem.
        const W: usize = 16;
        let samples = 64usize;
        let mut sampler = Sampler::new(b"by-ratio-window-den-v2-v1", SECP256K1_P);
        let mut hist = [0usize; 17];
        let mut windows = 0usize;
        for _ in 0..samples {
            let x = sampler.next();
            let mut delta = 1i64;
            let mut f = SInt::from_u(SECP256K1_P);
            let mut g = SInt::from_u(x);
            for _ in 0..35 {
                let f_low = low_signed_sint16_for_streaming_test(f);
                let g_low = low_signed_sint16_for_streaming_test(g);
                let bits = branch_bits_for_lowword_window(W, delta, f_low, g_low);
                let m = matrix_from_branch_bits(delta, &bits);
                let v = v2_i128_for_ratio_test(m.m01).min(16);
                hist[v] += 1;
                windows += 1;
                for _ in 0..W {
                    divstep_sint_state(&mut delta, &mut f, &mut g);
                }
            }
        }
        let weak = hist[0] + hist[1] + hist[2];
        eprintln!(
            "BY ratio window denominator v2(m01): windows={windows}, hist={hist:?}, weak_v1_v2={weak}"
        );
        assert!(weak > windows / 3, "m01 high-v2 enough to merit a separate follow-up");
        assert_eq!(hist[0], 0, "m01 should be even after scaled 16-step windows");
    }

    fn bitlen_u256_for_compact_pair_test(x: U256) -> usize {
        let limbs = x.as_limbs();
        for i in (0..4).rev() {
            if limbs[i] != 0 {
                return i * 64 + (64 - limbs[i].leading_zeros() as usize);
            }
        }
        0
    }

    fn bitlen_sint_for_compact_pair_test(x: SInt) -> usize {
        bitlen_u256_for_compact_pair_test(x.mag)
    }

    fn emit_variable_coeff_times_acc_lower_bound_for_test(
        b: &mut super::super::B,
        src: &[super::super::QubitId],
        coeff: &[super::super::QubitId],
        acc: &[super::super::QubitId],
    ) {
        // Lower-bound-ish implementation for selected matrix rows: for each
        // possible coefficient bit, controlled-add the shifted source into the
        // accumulator.  This omits coefficient lookup, signs, overflow cleanup,
        // modular corrections, and old-register cleanup, so if this is already
        // too large the selected fixed-matrix BY route is dead.
        for (shift, &ctrl) in coeff.iter().enumerate() {
            super::super::cucc_add_ctrl(b, src, &acc[shift..shift + src.len()], ctrl);
        }
    }

    #[test]
    fn selected_fixed_matrix_window_variable_coeff_lower_bound_kills_by() {
        // Hard gate from the user's instruction: do not build BY unless the
        // selected-window primitive is already SOTA-shaped.  The fixed-matrix
        // budget used classical coefficients.  In a real circuit those 16-step
        // matrix coefficients are pattern-selected quantum data.  Even the
        // optimistic variable-coefficient multiply-accumulate lower bound for a
        // single window is far above the ~10k/window target before lookup,
        // signs, q corrections, cleanup, and history allocation.
        const WIDTH: usize = 274;
        const CBITS: usize = 17;
        let mut b = super::super::B::new();
        let src = b.alloc_qubits(WIDTH);
        let coeff = b.alloc_qubits(CBITS);
        let acc = b.alloc_qubits(WIDTH + CBITS);
        let start_one = b.ops.len();
        emit_variable_coeff_times_acc_lower_bound_for_test(&mut b, &src, &coeff, &acc);
        let one_coeff_ccx = count_ccx(&b.ops[start_one..]);

        let mut b_win = super::super::B::new();
        let f = b_win.alloc_qubits(WIDTH);
        let g = b_win.alloc_qubits(WIDTH);
        let coeffs = (0..8).map(|_| b_win.alloc_qubits(CBITS)).collect::<Vec<_>>();
        let accs = (0..4).map(|_| b_win.alloc_qubits(WIDTH + CBITS)).collect::<Vec<_>>();
        let start = b_win.ops.len();
        // Four coeff*source terms form two new rows; four more are the minimal
        // adjugate old-clean terms.  This is still missing table lookup and q.
        emit_variable_coeff_times_acc_lower_bound_for_test(&mut b_win, &f, &coeffs[0], &accs[0]);
        emit_variable_coeff_times_acc_lower_bound_for_test(&mut b_win, &g, &coeffs[1], &accs[0]);
        emit_variable_coeff_times_acc_lower_bound_for_test(&mut b_win, &f, &coeffs[2], &accs[1]);
        emit_variable_coeff_times_acc_lower_bound_for_test(&mut b_win, &g, &coeffs[3], &accs[1]);
        emit_variable_coeff_times_acc_lower_bound_for_test(&mut b_win, &accs[0][..WIDTH], &coeffs[4], &accs[2]);
        emit_variable_coeff_times_acc_lower_bound_for_test(&mut b_win, &accs[1][..WIDTH], &coeffs[5], &accs[2]);
        emit_variable_coeff_times_acc_lower_bound_for_test(&mut b_win, &accs[0][..WIDTH], &coeffs[6], &accs[3]);
        emit_variable_coeff_times_acc_lower_bound_for_test(&mut b_win, &accs[1][..WIDTH], &coeffs[7], &accs[3]);
        let window_lower_ccx = count_ccx(&b_win.ops[start..]);
        eprintln!(
            "BY selected fixed-matrix variable-coeff lower bound: one_coeff_ccx={one_coeff_ccx}, window_lower_ccx={window_lower_ccx}, peak={}q",
            b_win.peak_qubits
        );
        assert!(one_coeff_ccx > 5_000, "controlled variable coeff multiply unexpectedly cheap");
        assert!(window_lower_ccx > 40_000, "selected matrix window may still be SOTA-shaped; revisit BY");
    }

    #[test]
    fn denominator_pair_fixed_slack_schedule_50_sidecar_on_samples() {
        // A fixed allocator is much easier than a data-dependent compactor: at
        // step k reserve enough low bits for the worst sampled |f|,|g| sizes and
        // use only the predeclared high-zone as history storage.  This test
        // checks that the same 49-bit sidecar works with per-step worst-case
        // sampled bitlengths, not just per-trace adaptive slack.  The fixed
        // schedule needs 50 bits on an independent seed because some traces
        // reach the full 560-step endpoint.
        let samples = 8192usize;
        let mut max_used_by_step = vec![0usize; 35 * 16 + 1];
        let mut sampler = Sampler::new(b"by-pair-fixed-slack-v1", SECP256K1_P);
        for _ in 0..samples {
            let x = sampler.next();
            let mut delta = 1i64;
            let mut f = SInt::from_u(SECP256K1_P);
            let mut g = SInt::from_u(x);
            for step in 0..=35 * 16 {
                let used = bitlen_sint_for_compact_pair_test(f) + bitlen_sint_for_compact_pair_test(g);
                max_used_by_step[step] = max_used_by_step[step].max(used);
                if step == 35 * 16 || g.is_zero() {
                    break;
                }
                divstep_sint_state(&mut delta, &mut f, &mut g);
            }
        }
        let mut worst_deficit = 0isize;
        let mut worst_step = 0usize;
        for (step, &used) in max_used_by_step.iter().enumerate() {
            if used == 0 && step != 0 {
                continue;
            }
            let slack = 512isize - used as isize;
            let deficit = step as isize - slack;
            if deficit > worst_deficit {
                worst_deficit = deficit;
                worst_step = step;
            }
        }
        eprintln!(
            "BY fixed slack schedule: samples={samples}, worst_step={worst_step}, worst_sidecar={worst_deficit} bits"
        );
        assert!(worst_deficit <= 50, "fixed slack schedule exceeded 50-bit sidecar");
    }

    #[test]
    fn denominator_pair_plus_50_sidecar_can_hold_raw_history_on_samples() {
        // Better low-qubit idea than the scalar ratio: keep the 256-bit f/g
        // denominator pair and stash consumed branch bits into high zero slack
        // as the pair shrinks, with a tiny sidecar for the initial lag.  Real
        // traces need at most 50 extra raw-history bits over the pair slack in
        // sampled checks, matching the 560-step vs 510-slack fixed-schedule
        // endpoint intuition.
        let samples = 8192usize;
        let mut sampler = Sampler::new(b"by-pair-slack-history-v1", SECP256K1_P);
        let mut worst_deficit = 0isize;
        let mut max_converge = 0usize;
        for _ in 0..samples {
            let x = sampler.next();
            let mut delta = 1i64;
            let mut f = SInt::from_u(SECP256K1_P);
            let mut g = SInt::from_u(x);
            for step in 0..(35 * 16) {
                if g.is_zero() {
                    max_converge = max_converge.max(step);
                    break;
                }
                divstep_sint_state(&mut delta, &mut f, &mut g);
                let used = bitlen_sint_for_compact_pair_test(f) + bitlen_sint_for_compact_pair_test(g);
                let slack = 512isize - used as isize;
                let deficit = (step + 1) as isize - slack;
                worst_deficit = worst_deficit.max(deficit);
                if step == 35 * 16 - 1 {
                    max_converge = max_converge.max(35 * 16);
                }
            }
        }
        eprintln!(
            "BY denominator-pair slack+sidecar history: samples={samples}, max_converge={max_converge}, worst_raw_history_deficit={worst_deficit} bits"
        );
        assert!(worst_deficit > 0, "pair slack unexpectedly stores all raw branch bits");
        assert!(worst_deficit <= 50, "sidecar exceeded the 50-bit target");
    }

    #[test]
    fn denominator_pair_slack_history_toy_exhaustive_warns_exact_tail() {
        // Exhaustive toy counterpart to the sampled secp slack-history check.
        // For a low-qubit BY route, the raw branch stream can be stored in the
        // shrinking f/g denominator pair only if the initial sidecar deficit
        // stays a small fraction of n, not another field-sized history bank.
        // Use a BY-like 560/256 step ratio on toy pseudo-Mersenne primes.
        let cases = [
            (8usize, 251u16),
            (10usize, 1021u16),
            (12usize, 4093u16),
            (14usize, 16381u16),
        ];
        let mut max_scaled = 0.0f64;
        for &(n, p16) in &cases {
            let steps = (560usize * n).div_ceil(256);
            let mut worst_deficit = 0isize;
            let mut worst_x = 0u16;
            let mut max_converge = 0usize;
            for x in 1..p16 {
                let mut delta = 1i64;
                let mut f = SInt::from_u(U256::from(p16 as u64));
                let mut g = SInt::from_u(U256::from(x as u64));
                for step in 0..steps {
                    if g.is_zero() {
                        max_converge = max_converge.max(step);
                        break;
                    }
                    divstep_sint_state(&mut delta, &mut f, &mut g);
                    let used = bitlen_sint_for_compact_pair_test(f) + bitlen_sint_for_compact_pair_test(g);
                    let slack = (2 * n) as isize - used as isize;
                    let deficit = (step + 1) as isize - slack;
                    if deficit > worst_deficit {
                        worst_deficit = deficit;
                        worst_x = x;
                    }
                    if step + 1 == steps {
                        max_converge = max_converge.max(steps);
                    }
                }
            }
            let scaled_to_256 = (worst_deficit.max(0) as f64) * 256.0 / n as f64;
            max_scaled = max_scaled.max(scaled_to_256);
            eprintln!(
                "BY toy slack-history: n={n}, steps={steps}, worst_deficit={worst_deficit}, scaled_to_256={scaled_to_256:.1}, worst_x={worst_x}, max_converge={max_converge}"
            );
            if n == 14 {
                println!("METRIC by_slack_history_toy_n14_steps={steps}");
                println!("METRIC by_slack_history_toy_n14_worst_deficit={worst_deficit}");
                println!("METRIC by_slack_history_toy_n14_scaled_to_256={scaled_to_256:.3}");
            }
            assert!(worst_deficit >= 0, "deficit accounting should not stay negative");
        }
        assert!(
            max_scaled > 120.0,
            "toy exhaustive sidecar stayed small; sampled secp slack route is stronger than expected"
        );
    }

    #[test]
    fn ratio_a_step_serial_inverse_budget_is_too_large() {
        // A bit-serial triangular inverse avoids large scratch, but each A step
        // still needs convolution work across the remaining active width.  Use
        // sum(t^2/2) over real A-step positions as a crude lower-order gate
        // proxy.  It is already multi-million before reversible cleanup and
        // before the modular replay itself.
        let samples = 64usize;
        let mut sampler = Sampler::new(b"by-ratio-a-serial-budget-v1", SECP256K1_P);
        let mut total_proxy = 0f64;
        let mut max_proxy = 0f64;
        for _ in 0..samples {
            let x = sampler.next();
            let mut delta = 1i64;
            let mut f = SInt::from_u(SECP256K1_P);
            let mut g = SInt::from_u(x);
            let mut proxy = 0f64;
            for step in 0..(35 * 16) {
                let t = (35 * 16 - step) as f64;
                if delta > 0 && g.bit0() {
                    proxy += 0.5 * t * t;
                }
                divstep_sint_state(&mut delta, &mut f, &mut g);
            }
            total_proxy += proxy;
            max_proxy = max_proxy.max(proxy);
        }
        let mean_proxy = total_proxy / samples as f64;
        let mean_proxy_rounded = mean_proxy.round() as usize;
        let max_proxy_rounded = max_proxy.round() as usize;
        let replay_body_projection = 2_645_196usize;
        let projected_toffoli = replay_body_projection + mean_proxy_rounded;
        let gap_to_2700k = projected_toffoli as isize - 2_700_000isize;
        eprintln!(
            "BY ratio serial A-inverse proxy: mean_t2_over2={mean_proxy:.0}, max_t2_over2={max_proxy:.0}, projected_toffoli={projected_toffoli}, gap={gap_to_2700k}"
        );
        println!("METRIC scratch600_full_ratio_a_serial_proxy_mean_ccx={mean_proxy_rounded}");
        println!("METRIC scratch600_full_ratio_a_serial_proxy_max_ccx={max_proxy_rounded}");
        println!("METRIC scratch600_full_ratio_projected_toffoli={projected_toffoli}");
        println!("METRIC scratch600_full_ratio_gap_to_2700k={gap_to_2700k}");
        assert!(mean_proxy > 2_000_000.0);
    }

    #[test]
    fn ratio_a_step_is_inverse_dense_and_common() {
        // The 560-bit full-ratio selector solves state size, but not
        // automatically gate cost.  Its A update contains a modular inverse of
        // the current odd h.  On a toy 16-bit odd input, even one output bit has
        // almost maximal ANF degree and a dense monomial set.  A steps are also
        // common in real 560-step secp256k1 traces, so a naive per-A inverse is
        // not a SOTA route.
        let (degree, density) = ratio_a_step_anf_stats_for_test(16, 14);
        let samples = 64usize;
        let mut sampler = Sampler::new(b"by-ratio-a-count-v1", SECP256K1_P);
        let mut total_a = 0usize;
        let mut max_a = 0usize;
        for _ in 0..samples {
            let x = sampler.next();
            let mut delta = 1i64;
            let mut f = SInt::from_u(SECP256K1_P);
            let mut g = SInt::from_u(x);
            let mut a_count = 0usize;
            for _ in 0..(35 * 16) {
                if delta > 0 && g.bit0() {
                    a_count += 1;
                }
                divstep_sint_state(&mut delta, &mut f, &mut g);
            }
            total_a += a_count;
            max_a = max_a.max(a_count);
        }
        let mean_a = total_a as f64 / samples as f64;
        eprintln!(
            "BY ratio A-step obstruction: anf_degree={degree}, density={density}/32768, mean_a_steps={mean_a:.1}, max_a_steps={max_a}"
        );
        assert!(degree >= 14);
        assert!(density > 7_000);
        assert!(mean_a > 100.0);
    }

    #[test]
    fn full_ratio_state_streams_all_branches_in_560_bits() {
        // Stronger selector compression: since BY branch decisions depend only
        // on delta and the 2-adic ratio h=g/f, start with
        // h0 = x / p mod 2^560 and never materialize the denominator pair or
        // the affine x-column/carry decomposition.  The active ratio width
        // shrinks by one bit per divstep; the vacated bits can be the branch
        // history in a reversible implementation.  This makes the whole
        // denominator selector/history a 560-bit object, inside the low-qubit
        // target's ~600-bit allowance.
        const TOTAL_BITS: usize = 35 * 16;
        let samples = 64usize;
        let p_big = big9_from_u256(SECP256K1_P);
        let p_inv = big9_inv_odd_mod_pow2(&p_big, TOTAL_BITS);
        let mut sampler = Sampler::new(b"by-full-ratio-state-v1", SECP256K1_P);
        for _ in 0..samples {
            let x = sampler.next();
            let x_big = big9_from_u256(x);
            let mut h = big9_mul_mod(&x_big, &p_inv, TOTAL_BITS);
            let mut delta = 1i64;
            let mut f = SInt::from_u(SECP256K1_P);
            let mut g = SInt::from_u(x);
            for t in (1..=TOTAL_BITS).rev() {
                let (next_d, next_h, odd_h) = h_ratio_step_big9_for_streaming_test(delta, h, t);
                assert_eq!(odd_h, g.bit0(), "full-ratio odd mismatch at t={t}");
                divstep_sint_state(&mut delta, &mut f, &mut g);
                assert_eq!(next_d, delta, "full-ratio delta mismatch at t={t}");
                delta = next_d;
                h = next_h;
            }
            assert!(g.is_zero(), "560 BY ratio steps did not converge g");
            assert_eq!(h, [0u64; 9], "full ratio did not taper to zero");
        }
        eprintln!("BY full ratio selector: samples={samples}, total_bits={TOTAL_BITS}");
        assert_eq!(TOTAL_BITS, 560);
    }

    #[test]
    fn tail_ratio_state_streams_remaining_branches_in_304_bits() {
        // After the first 16 windows, the 256-bit x tail is exhausted.  Instead
        // of keeping both carry rows (528 bits), keep the 2-adic ratio
        // h=g/f mod 2^304.  The closed BY ratio update streams the remaining
        // 304 branch bits exactly while the active h width tapers by one bit per
        // divstep.  The same 304-qubit register can also hold the consumed tail
        // branch history in its vacated high bits in a reversible circuit.
        const PREFIX_WINDOWS: usize = 16;
        const REM_BITS: usize = (35 - PREFIX_WINDOWS) * 16;
        let samples = 64usize;
        let mut sampler = Sampler::new(b"by-tail-ratio-state-v1", SECP256K1_P);
        for _ in 0..samples {
            let x = sampler.next();
            let mut delta = 1i64;
            let mut f = SInt::from_u(SECP256K1_P);
            let mut g = SInt::from_u(x);
            for _ in 0..(PREFIX_WINDOWS * 16) {
                divstep_sint_state(&mut delta, &mut f, &mut g);
            }
            let f_res = sint_residue_u512_for_streaming_test(f, REM_BITS);
            assert!(f_res.bit(0), "f must remain odd for BY ratio state");
            let g_res = sint_residue_u512_for_streaming_test(g, REM_BITS);
            let mut h = (g_res * inv_odd_u512_pow2_for_streaming_test(f_res, REM_BITS))
                & mask_u512_bits_for_streaming_test(REM_BITS);
            let mut d_h = delta;
            for t in (1..=REM_BITS).rev() {
                let (next_d, next_h, odd_h) = h_ratio_step_u512_for_streaming_test(d_h, h, t);
                assert_eq!(odd_h, g.bit0(), "ratio odd mismatch at t={t}");
                divstep_sint_state(&mut delta, &mut f, &mut g);
                assert_eq!(next_d, delta, "ratio delta mismatch at t={t}");
                d_h = next_d;
                h = next_h;
            }
            assert!(g.is_zero(), "remaining BY steps did not converge g");
            assert_eq!(h, U512::ZERO, "ratio state did not taper to zero");
        }
        eprintln!(
            "BY tail ratio selector: samples={samples}, prefix_windows={PREFIX_WINDOWS}, tail_bits={REM_BITS}"
        );
        assert_eq!(REM_BITS, 304);
    }

    #[test]
    fn truncated_x_column_selector_state_is_not_locally_reversible() {
        // The 9-limb x-column in the separate-width selector is a low-residue
        // scratch model, not a standalone reversible state.  Even the C branch
        // update contains b0 <- 2*b0 (mod 2^k), which is two-to-one.  Therefore
        // a real circuit must recompute this x-column from retained pattern
        // history, keep wider exact coefficients, or use a nontrivial MBUC
        // cleanup.  Do not wire the 816-bit model as an in-place rolling
        // register without solving this.
        let k = 9 * 16;
        let modulus = U512::from(1u64) << k;
        let mask = modulus - U512::from(1u64);
        let x = U512::from(0x1234u64);
        let y0: U512 = (x << 1usize) & mask;
        let y1: U512 = ((x + (U512::from(1u64) << (k - 1))) << 1usize) & mask;
        eprintln!(
            "BY truncated x-column non-injectivity: width_bits={k}, example_output_low16={}",
            y0.as_limbs()[0] & 0xffff
        );
        assert_eq!(y0, y1, "doubling modulo 2^k unexpectedly injective");
    }

    #[test]
    fn first16_pattern_history_entropy_is_low_gate_but_not_low_qubit() {
        // If the 288-bit x-column workspace is recomputed from early pattern
        // history instead of carried, only the first 16 window patterns are
        // needed after the x tail is exhausted.  They compress well enough for
        // a 1425q-style low-gate point, but not for the 1175q low-qubit target
        // when added to the 528-bit carry core.
        use std::collections::HashMap;
        const W: usize = 16;
        const WINDOWS: usize = 16;
        let samples = 5_000usize;
        let mut counts: Vec<HashMap<u16, usize>> = (0..WINDOWS).map(|_| HashMap::new()).collect();
        let mut sampler = Sampler::new(b"by-first16-pattern-entropy-v1", SECP256K1_P);
        for _ in 0..samples {
            let x = sampler.next();
            let mut delta = 1i64;
            let mut f = SInt::from_u(SECP256K1_P);
            let mut g = SInt::from_u(x);
            for win in 0..WINDOWS {
                let mut pat = 0u16;
                for i in 0..W {
                    if g.bit0() {
                        pat |= 1u16 << i;
                    }
                    divstep_sint_state(&mut delta, &mut f, &mut g);
                }
                *counts[win].entry(pat).or_insert(0) += 1;
            }
        }
        let mut entropy = 0.0f64;
        let mut fixed_bits = 0usize;
        for c in &counts {
            fixed_bits += usize::BITS as usize - (c.len() - 1).leading_zeros() as usize;
            for &n in c.values() {
                let p = n as f64 / samples as f64;
                entropy -= p * p.log2();
            }
        }
        let carry_core = 528usize;
        let low_gate_persistent = carry_core + fixed_bits;
        let low_gate_scratch = 1425usize - 512usize;
        let low_qubit_scratch = 1175usize - 512usize;
        let gap_to_low_gate = low_gate_persistent as isize - low_gate_scratch as isize;
        let gap_to_low_qubit = low_gate_persistent as isize - low_qubit_scratch as isize;
        eprintln!(
            "BY first16 pattern history entropy: H≈{entropy:.1}, fixed_bits={fixed_bits}, carry_plus_fixed={low_gate_persistent}, gap_low_gate={gap_to_low_gate}, gap_low_qubit={gap_to_low_qubit}"
        );
        println!("METRIC by_first16_pattern_entropy_bits={entropy:.3}");
        println!("METRIC by_first16_pattern_fixed_bits={fixed_bits}");
        println!("METRIC by_tail_carry_core_bits={carry_core}");
        println!("METRIC by_first16_carry_plus_fixed_bits={low_gate_persistent}");
        println!("METRIC by_first16_gap_to_lowgate_scratch_bits={gap_to_low_gate}");
        println!("METRIC by_first16_gap_to_lowqubit_scratch_bits={gap_to_low_qubit}");
        assert!(entropy < 220.0, "first16 pattern history too large for low-gate selector plan");
        assert!(low_gate_persistent < low_gate_scratch, "carry+first16 history misses 1425q low-gate budget");
        assert!(low_gate_persistent > low_qubit_scratch, "carry+first16 history would already hit 1175q low-qubit budget");
    }

    #[test]
    fn first16_tail_selector_lowgate_toffoli_margin_is_tight() {
        // Go/no-go budget for the only BY streaming-selector variant that
        // currently fits a Google memory regime: first-16 pattern history plus
        // the 528-bit tail carry core fits the 1425q low-gate scratch allowance
        // but misses the 1175q low-qubit allowance.  The low-gate Toffoli
        // target is stricter, so the selector generator must be almost free.
        // Use measured local oracles from prior tests:
        //   lowword pattern oracle:       5,952 CCX / 16-step window
        //   lowword pattern+q oracle:     9,408 CCX / 16-step window
        // and the fast scaled-BY whole-point budget from the scratch model:
        //   scaffold_after_div + fast BY body ~= 1,938,476 CCX
        let google_low_gate = 2_100_000isize;
        let fast_scaled_by_projected = 1_938_476isize;
        let margin = google_low_gate - fast_scaled_by_projected;
        let first16_windows = 16isize;
        let pattern_oracle_per_window = 5_952isize;
        let pattern_q_oracle_per_window = 9_408isize;
        let one_pass_pattern = first16_windows * pattern_oracle_per_window;
        let compute_uncompute_pattern = 2 * one_pass_pattern;
        let one_pass_pattern_q = first16_windows * pattern_q_oracle_per_window;
        let rem_after_pattern = margin - one_pass_pattern;
        let rem_after_pattern_q = margin - one_pass_pattern_q;
        let gap_after_pattern_compute_uncompute = compute_uncompute_pattern - margin;
        eprintln!("BY first16/tail low-gate budget: margin={margin}, one_pass_pattern={one_pass_pattern}, rem_after_pattern={rem_after_pattern}, one_pass_pattern_q={one_pass_pattern_q}, rem_after_pattern_q={rem_after_pattern_q}, cu_pattern_gap={gap_after_pattern_compute_uncompute}");
        println!("METRIC by_first16_lowgate_margin_ccx={margin}");
        println!("METRIC by_first16_pattern_onepass_ccx={one_pass_pattern}");
        println!("METRIC by_first16_pattern_remaining_ccx={rem_after_pattern}");
        println!("METRIC by_first16_pattern_q_onepass_ccx={one_pass_pattern_q}");
        println!("METRIC by_first16_pattern_q_remaining_ccx={rem_after_pattern_q}");
        println!("METRIC by_first16_pattern_compute_uncompute_gap_ccx={gap_after_pattern_compute_uncompute}");
        assert!(rem_after_pattern > 0, "even one-pass first16 pattern generation misses low-gate budget");
        assert!(rem_after_pattern_q > 0, "pattern+q first16 oracle would already miss low-gate budget");
        assert!(gap_after_pattern_compute_uncompute > 0, "compute+uncompute pattern unexpectedly fits low-gate margin");
    }

    fn naf_weight_i128_for_by_budget(x: i128) -> usize {
        let mut n = x.unsigned_abs();
        let mut w = 0usize;
        while n != 0 {
            if (n & 1) == 0 {
                n >>= 1;
            } else {
                w += 1;
                if (n & 3) == 1 || n == 1 {
                    n = (n - 1) >> 1;
                } else {
                    n = (n + 1) >> 1;
                }
            }
        }
        w
    }

    fn row_csd_update_lb_for_by_budget(a: i128, b: i128, width_a: usize, width_b: usize) -> usize {
        let wa = naf_weight_i128_for_by_budget(a);
        let wb = naf_weight_i128_for_by_budget(b);
        let raw = wa * width_a + wb * width_b;
        let free_first = [if wa > 0 { width_a } else { 0 }, if wb > 0 { width_b } else { 0 }]
            .into_iter()
            .max()
            .unwrap_or(0);
        raw.saturating_sub(free_first)
    }

    #[test]
    fn first16_tail_carry_update_lower_bound_is_tight_but_not_dead() {
        // Hard-piece-first check for the first16/tail BY low-gate subpath. The
        // previous test left only 66,292 CCX after one-pass first16 pattern
        // generation. This gives the tail carry update an unrealistically
        // generous lower bound: fixed-pattern signed-CSD constant row updates,
        // one source term per row copied for free, shifts by 2^16 free, and no
        // charge for controls, cleanup, low-word extraction, or reversibility.
        // The bound fits, but barely; any real circuit overhead can kill it.
        const W: usize = 16;
        const PREFIX_WINDOWS: usize = 16;
        const WINDOWS: usize = 35;
        const C0_WIDTH: usize = 17 * 16;
        const C1_WIDTH: usize = 16 * 16;
        let samples = 512usize;
        let mut sampler = Sampler::new(b"by-tail-carry-update-lb-v1", SECP256K1_P);
        let mut costs = Vec::with_capacity(samples);
        for _ in 0..samples {
            let x = sampler.next();
            let mut delta = 1i64;
            let mut f = SInt::from_u(SECP256K1_P);
            let mut g = SInt::from_u(x);
            let mut cost = 0usize;
            for win in 0..WINDOWS {
                let f_low = low_signed_sint16_for_streaming_test(f);
                let g_low = low_signed_sint16_for_streaming_test(g);
                let bits = branch_bits_for_lowword_window(W, delta, f_low, g_low);
                let m = matrix_from_branch_bits(delta, &bits);
                if win >= PREFIX_WINDOWS {
                    cost += row_csd_update_lb_for_by_budget(m.m00, m.m01, C0_WIDTH, C1_WIDTH);
                    cost += row_csd_update_lb_for_by_budget(m.m10, m.m11, C0_WIDTH, C1_WIDTH);
                }
                for _ in 0..W {
                    divstep_sint_state(&mut delta, &mut f, &mut g);
                }
            }
            costs.push(cost);
        }
        costs.sort_unstable();
        let mean = costs.iter().sum::<usize>() as f64 / samples as f64;
        let p90 = costs[(samples * 90) / 100];
        let p99 = costs[(samples * 99) / 100];
        let max = *costs.last().unwrap();
        let remaining_after_first16_pattern = 66_292isize;
        let p99_gap = p99 as isize - remaining_after_first16_pattern;
        let max_gap = max as isize - remaining_after_first16_pattern;
        eprintln!("BY tail carry update optimistic CSD lower bound: mean={mean:.1}, p90={p90}, p99={p99}, max={max}, p99_gap={p99_gap}, max_gap={max_gap}");
        println!("METRIC by_tail_carry_update_lb_mean_ccx={mean:.3}");
        println!("METRIC by_tail_carry_update_lb_p90_ccx={p90}");
        println!("METRIC by_tail_carry_update_lb_p99_ccx={p99}");
        println!("METRIC by_tail_carry_update_lb_max_ccx={max}");
        println!("METRIC by_tail_carry_update_lb_gap_to_remaining_ccx={p99_gap}");
        println!("METRIC by_tail_carry_update_lb_max_gap_to_remaining_ccx={max_gap}");
        assert!(p99_gap < 0, "even optimistic p99 tail carry lower bound misses low-gate margin");
        assert!(max_gap < 0, "sampled max tail carry lower bound misses low-gate margin");
    }

    fn tail_carry_fresh_update_cost_for_by_budget(mtx: TransitionMatrix) -> usize {
        const C0_WIDTH: usize = 17 * 16;
        const C1_WIDTH: usize = 16 * 16;
        const TMP_WIDTH: usize = C0_WIDTH + 16;
        let mut b = super::super::B::new();
        let c0 = b.alloc_qubits(C0_WIDTH);
        let c1 = b.alloc_qubits(C1_WIDTH);
        let y0 = b.alloc_qubits(TMP_WIDTH);
        let y1 = b.alloc_qubits(TMP_WIDTH);
        add_signed_coeff_times_for_cost(&mut b, mtx.m00, &c0, &y0);
        add_signed_coeff_times_for_cost(&mut b, mtx.m01, &c1, &y0);
        add_signed_coeff_times_for_cost(&mut b, mtx.m10, &c0, &y1);
        add_signed_coeff_times_for_cost(&mut b, mtx.m11, &c1, &y1);
        arith_shift_right_inplace_for_cost(&mut b, &y0, 16);
        arith_shift_right_inplace_for_cost(&mut b, &y1, 16);
        count_ccx(&b.ops)
    }

    #[test]
    fn first16_tail_carry_fresh_update_circuit_kills_lowgate_margin() {
        // The previous CSD arithmetic lower bound fit with only ~10k CCX to
        // spare. This is the promised actual-circuit checkpoint. It is still
        // generous: it computes each tail carry update into fresh rows, does not
        // clean the old carry rows, does not uncompute low-word/pattern logic,
        // and uses a fixed known matrix rather than quantum-selected controls.
        // If even this forward-only fresh update exceeds the 66,292 CCX margin,
        // the first16/tail low-gate subpath is not worth deeper engineering.
        use std::collections::HashMap;
        const W: usize = 16;
        const PREFIX_WINDOWS: usize = 16;
        const WINDOWS: usize = 35;
        const SAMPLES: usize = 64;
        let mut sampler = Sampler::new(b"by-tail-carry-fresh-circuit-v1", SECP256K1_P);
        let mut cache: HashMap<(i128, i128, i128, i128), usize> = HashMap::new();
        let mut costs = Vec::with_capacity(SAMPLES);
        for _ in 0..SAMPLES {
            let x = sampler.next();
            let mut delta = 1i64;
            let mut f = SInt::from_u(SECP256K1_P);
            let mut g = SInt::from_u(x);
            let mut cost = 0usize;
            for win in 0..WINDOWS {
                let f_low = low_signed_sint16_for_streaming_test(f);
                let g_low = low_signed_sint16_for_streaming_test(g);
                let bits = branch_bits_for_lowword_window(W, delta, f_low, g_low);
                let m = matrix_from_branch_bits(delta, &bits);
                if win >= PREFIX_WINDOWS {
                    let key = (m.m00, m.m01, m.m10, m.m11);
                    let c = *cache.entry(key).or_insert_with(|| tail_carry_fresh_update_cost_for_by_budget(m));
                    cost += c;
                }
                for _ in 0..W {
                    divstep_sint_state(&mut delta, &mut f, &mut g);
                }
            }
            costs.push(cost);
        }
        costs.sort_unstable();
        let mean = costs.iter().sum::<usize>() as f64 / SAMPLES as f64;
        let p90 = costs[(SAMPLES * 90) / 100];
        let p99 = costs[(SAMPLES * 99) / 100];
        let max = costs[SAMPLES - 1];
        let remaining_after_first16_pattern = 66_292isize;
        let p90_gap = p90 as isize - remaining_after_first16_pattern;
        let min_cached_window = cache.values().copied().min().unwrap_or(0);
        let max_cached_window = cache.values().copied().max().unwrap_or(0);
        eprintln!("BY tail carry fresh-update circuit: unique_matrices={}, min_window={min_cached_window}, max_window={max_cached_window}, mean={mean:.1}, p90={p90}, p99={p99}, max={max}, p90_gap={p90_gap}", cache.len());
        println!("METRIC by_tail_carry_fresh_unique_matrices={}", cache.len());
        println!("METRIC by_tail_carry_fresh_window_min_ccx={min_cached_window}");
        println!("METRIC by_tail_carry_fresh_window_max_ccx={max_cached_window}");
        println!("METRIC by_tail_carry_fresh_mean_ccx={mean:.3}");
        println!("METRIC by_tail_carry_fresh_p90_ccx={p90}");
        println!("METRIC by_tail_carry_fresh_p99_ccx={p99}");
        println!("METRIC by_tail_carry_fresh_max_ccx={max}");
        println!("METRIC by_tail_carry_fresh_gap_to_remaining_ccx={p90_gap}");
        assert!(p90_gap > 0, "fresh forward-only tail carry update fits low-gate margin; build real reversible version");
    }

    #[test]
    fn projective_normalized_streaming_selector_loses_high_bits() {
        // Tempting compression: since BY branch choices are invariant under a
        // common odd scale, normalize the folded selector so c0=1 and keep only
        // three entries `(b0,b1,c1)`.  This would be 3*17*16=816 bits, much
        // closer to the target.  It fails because repeated normalization throws
        // away high 2-adic information needed by later windows; this is the
        // same obstruction as the earlier h-only state, now in affine-streaming
        // coordinates.
        let samples = 64usize;
        let mut sampler = Sampler::new(b"by-projective-streaming-selector-v1", SECP256K1_P);
        let mut failures = 0usize;
        for _ in 0..samples {
            let x = sampler.next();
            if !streaming_selector_projective_normalized_matches_for_test(x, 17) {
                failures += 1;
            }
        }
        let projected_state_bits = 3 * 17 * 16;
        eprintln!(
            "BY projective-normalized streaming selector: samples={samples}, failures={failures}, projected_state_bits={projected_state_bits}"
        );
        assert!(failures > 0, "projective normalized selector unexpectedly exact on all samples");
        assert!(projected_state_bits < 1088, "projective normalization would not reduce state");
    }

    #[test]
    fn naive_variable_coefficient_jump_apply_is_too_expensive() {
        // If we synthesize the w-bit matrix entries into quantum coefficient
        // registers and then multiply each full-width row by every possible
        // coefficient bit, cost scales with bit-width rather than popcount.
        // This quantifies that dead end: selected matrices must be applied via
        // a better decomposition/control scheme than generic variable small ×
        // wide multiplication.
        const WIDTH: usize = 274;
        const W: usize = 16;
        let mut b = super::super::B::new();
        let src = b.alloc_qubits(WIDTH);
        let dst = b.alloc_qubits(WIDTH + W);
        let coeff_bits = b.alloc_qubits(W + 1);
        let start = b.ops.len();
        for shift in 0..=W {
            let len = src.len().min(dst.len() - shift);
            let src_slice = src[..len].to_vec();
            let dst_slice = dst[shift..shift + len].to_vec();
            cadd_qq_fast_for_cost(&mut b, &dst_slice, &src_slice, coeff_bits[shift]);
        }
        let one_coeff_ccx = count_ccx(&b.ops[start..]);
        let pair_update_cleanup_ccx = one_coeff_ccx * 8; // 4 P entries + 4 scaled-adjugate entries.
        let approx_two_pair_35 = pair_update_cleanup_ccx as f64 * 2.0 * 35.0;
        eprintln!(
            "naive variable BY coefficient apply: one_coeff_ccx={one_coeff_ccx}, pair_update_cleanup_ccx≈{pair_update_cleanup_ccx}, two_pair_35_windows≈{approx_two_pair_35:.0}"
        );
        assert!(approx_two_pair_35 > 3_000_000.0, "naive variable coefficient apply unexpectedly viable");
    }

    #[test]
    fn by_microstep_inplace_cost_model_is_not_the_jump_win() {
        // Low-scratch in-place BY microsteps are algebraically clean but they
        // pay controlled full-width additions every bit. This test keeps us
        // honest: the SOTA-shaped path needs jumped/selected matrices, not 550
        // raw coherent microsteps, unless the controlled-add implementation is
        // radically improved.
        const N: usize = 256;
        const WIDTH: usize = 274;
        let p = SECP256K1_P;
        let mut b = super::super::B::new();
        let a_ctrl = b.alloc_qubit(); // A branch: delta>0 && odd
        let b_ctrl = b.alloc_qubit(); // B branch: odd && !A
        let f = b.alloc_qubits(WIDTH);
        let g = b.alloc_qubits(WIDTH);
        let r = b.alloc_qubits(N);
        let s = b.alloc_qubits(N);

        let start = b.ops.len();
        // Denominator pair: g +=/-= f on odd, then f += g on A.
        cadd_qq_fast_for_cost(&mut b, &g, &f, b_ctrl);
        csub_qq_fast_for_cost(&mut b, &g, &f, a_ctrl);
        cadd_qq_fast_for_cost(&mut b, &f, &g, a_ctrl);
        // Tagged modular channel mirrors the same shears, then doubles top.
        super::super::cmod_add_qq(&mut b, &s, &r, b_ctrl, p);
        super::super::cmod_sub_qq(&mut b, &s, &r, a_ctrl, p);
        super::super::cmod_add_qq(&mut b, &r, &s, a_ctrl, p);
        super::super::mod_double_inplace_fast(&mut b, &r, p);
        let ccx = count_ccx(&b.ops[start..]);
        let approx_total = ccx as f64 * 550.0;
        eprintln!(
            "BY raw microstep in-place cost model: ccx_per_step={ccx}, approx_550≈{approx_total:.0}, peak={}q",
            b.peak_qubits
        );
        assert!(approx_total > 1_500_000.0, "raw microsteps unexpectedly competitive; revisit jump need");
    }

    fn signed_i128_mod_p(x: i128, p: U256) -> U256 {
        if x >= 0 {
            U256::from(x as u128) % p
        } else {
            let mag = U256::from(x.unsigned_abs());
            if mag.is_zero() { U256::ZERO } else { p.wrapping_sub(mag % p) }
        }
    }

    fn popcount_u256(x: U256) -> usize {
        (0..256).filter(|&i| x.bit(i)).count()
    }

    fn u256_to_u512_for_by_tests(x: U256) -> U512 {
        U512::from_limbs([
            x.as_limbs()[0],
            x.as_limbs()[1],
            x.as_limbs()[2],
            x.as_limbs()[3],
            0,
            0,
            0,
            0,
        ])
    }

    fn mod_mul_two_small_coeffs_acc_for_cost(
        b: &mut super::super::B,
        src: &[super::super::QubitId],
        c0: i128,
        acc0: &[super::super::QubitId],
        c1: i128,
        acc1: &[super::super::QubitId],
        p: U256,
    ) {
        if c0 == 0 && c1 == 0 {
            return;
        }
        let n = src.len();
        let tmp = b.alloc_qubits(n);
        for i in 0..n {
            b.cx(src[i], tmp[i]);
        }
        let mag0 = c0.unsigned_abs();
        let mag1 = c1.unsigned_abs();
        let top0 = if mag0 == 0 { 0 } else { 127 - mag0.leading_zeros() as usize };
        let top1 = if mag1 == 0 { 0 } else { 127 - mag1.leading_zeros() as usize };
        let top = top0.max(top1);
        for i in 0..=top {
            if ((mag0 >> i) & 1) != 0 {
                if c0 < 0 {
                    super::super::mod_sub_qq_fast(b, acc0, &tmp, p);
                } else {
                    super::super::mod_add_qq_fast(b, acc0, &tmp, p);
                }
            }
            if ((mag1 >> i) & 1) != 0 {
                if c1 < 0 {
                    super::super::mod_sub_qq_fast(b, acc1, &tmp, p);
                } else {
                    super::super::mod_add_qq_fast(b, acc1, &tmp, p);
                }
            }
            if i < top {
                super::super::mod_double_inplace_fast(b, &tmp, p);
            }
        }
        for _ in 0..top {
            super::super::mod_halve_inplace_fast(b, &tmp, p);
        }
        for i in 0..n {
            b.cx(src[i], tmp[i]);
        }
        b.free_vec(&tmp);
    }

    fn emit_scaled_modular_pair_update_with_sparse_cleanup_for_cost(
        b: &mut super::super::B,
        m: TransitionMatrix,
        w: usize,
        p: U256,
    ) {
        // Coefficient convention: C' = 2^-w · P · C (mod p). Forward rows use
        // sparse P followed by w modular halvings; cleanup uses sparse adj(P),
        // avoiding the dense 2^-w inverse constants. The row former shares one
        // doubling walk of each source across both destination rows.
        let x0 = b.alloc_qubits(256);
        let x1 = b.alloc_qubits(256);
        let y0 = b.alloc_qubits(256);
        let y1 = b.alloc_qubits(256);

        mod_mul_two_small_coeffs_acc_for_cost(b, &x0, m.m00, &y0, m.m10, &y1, p);
        mod_mul_two_small_coeffs_acc_for_cost(b, &x1, m.m01, &y0, m.m11, &y1, p);
        for _ in 0..w {
            super::super::mod_halve_inplace_fast(b, &y0, p);
            super::super::mod_halve_inplace_fast(b, &y1, p);
        }

        let inv = scaled_inverse_matrix(m, w); // sparse adjugate with det sign.
        mod_mul_two_small_coeffs_acc_for_cost(b, &y0, -inv.m00, &x0, -inv.m10, &x1, p);
        mod_mul_two_small_coeffs_acc_for_cost(b, &y1, -inv.m01, &x0, -inv.m11, &x1, p);
        let _ = (x0, x1, y0, y1);
    }

    #[test]
    fn modular_primitive_cost_breakdown_for_by_rows() {
        let p = SECP256K1_P;
        let mut b = super::super::B::new();
        let a = b.alloc_qubits(256);
        let acc = b.alloc_qubits(256);
        let start_add = b.ops.len();
        super::super::mod_add_qq_fast(&mut b, &acc, &a, p);
        let add_ccx = count_ccx(&b.ops[start_add..]);
        let start_sub = b.ops.len();
        super::super::mod_sub_qq_fast(&mut b, &acc, &a, p);
        let sub_ccx = count_ccx(&b.ops[start_sub..]);
        let start_double = b.ops.len();
        super::super::mod_double_inplace_fast(&mut b, &acc, p);
        let double_ccx = count_ccx(&b.ops[start_double..]);
        let start_halve = b.ops.len();
        super::super::mod_halve_inplace_fast(&mut b, &acc, p);
        let halve_ccx = count_ccx(&b.ops[start_halve..]);
        eprintln!(
            "mod primitive costs for BY rows: add={add_ccx}, sub={sub_ccx}, double={double_ccx}, halve={halve_ccx}, peak={}q",
            b.peak_qubits
        );
        assert!(add_ccx > 0 && halve_ccx > 0);
    }

    fn add_shifted_small_reg_for_cost(
        b: &mut super::super::B,
        small: &[super::super::QubitId],
        acc: &[super::super::QubitId],
        shift: usize,
        subtract: bool,
    ) {
        if shift >= acc.len() {
            return;
        }
        let len = acc.len() - shift;
        let tmp = b.alloc_qubits(len);
        let copy_len = small.len().min(len);
        for i in 0..copy_len {
            b.cx(small[i], tmp[i]);
        }
        let acc_slice = acc[shift..].to_vec();
        if subtract {
            super::super::sub_nbit_qq_fast(b, &tmp, &acc_slice);
        } else {
            super::super::add_nbit_qq_fast(b, &tmp, &acc_slice);
        }
        for i in 0..copy_len {
            b.cx(small[i], tmp[i]);
        }
        b.free_vec(&tmp);
    }

    fn emit_approx_batched_halve16_canonical(b: &mut super::super::B, v: &[super::super::QubitId]) {
        assert!(v.len() >= 274);
        const W: usize = 16;
        let m = b.alloc_qubits(W);
        let pinv = 51_919u64;
        let neg_pinv = ((!pinv).wrapping_add(1)) & ((1u64 << W) - 1);
        for bit_i in 0..W {
            if ((neg_pinv >> bit_i) & 1) != 0 {
                let len = W - bit_i;
                let src = v[..len].to_vec();
                let dst = m[bit_i..W].to_vec();
                super::super::add_nbit_qq_fast(b, &src, &dst);
            }
        }
        for &sh in &[0usize, 4, 6, 7, 8, 9, 32] {
            add_shifted_small_reg_for_cost(b, &m, v, sh, true);
        }
        add_shifted_small_reg_for_cost(b, &m, v, 256, false);
        for i in 0..(v.len() - W) {
            b.swap(v[i], v[i + W]);
        }
        for i in 0..W {
            b.cx(v[240 + i], m[i]);
        }
        b.free_vec(&m);
    }

    fn emit_approx_batched_halve16_for_cost(b: &mut super::super::B, v: &[super::super::QubitId]) {
        // Approximate canonical modular division by 2^16 for secp256k1:
        //   m = -v_low * p^{-1} mod 2^16;
        //   v <- (v + m*p) >> 16.
        // Since p=2^256-c, adding m*p is adding m at bit 256 and subtracting
        // m*c with c=2^32+977 (bits 0,4,6,7,8,9,32). For almost all inputs,
        // m is recovered from the top 16 output bits; rare small-input borrow
        // cases are a negligible approximate-DIV exception.
        assert!(v.len() >= 274);
        const W: usize = 16;
        let m = b.alloc_qubits(W);
        let pinv = 51_919u64; // p^{-1} mod 2^16 for secp256k1.
        let neg_pinv = ((!pinv).wrapping_add(1)) & ((1u64 << W) - 1);
        for bit_i in 0..W {
            if ((neg_pinv >> bit_i) & 1) != 0 {
                let len = W - bit_i;
                let src = v[..len].to_vec();
                let dst = m[bit_i..W].to_vec();
                super::super::add_nbit_qq_fast(b, &src, &dst);
            }
        }
        for &sh in &[0usize, 4, 6, 7, 8, 9, 32] {
            add_shifted_small_reg_for_cost(b, &m, v, sh, true);
        }
        add_shifted_small_reg_for_cost(b, &m, v, 256, false);
        // Right shift by 16 is a wire/swap layer. For this cost probe we only
        // model Toffoli, so no gates are needed. Approx-uncompute m from the
        // top output bits (v[256..272] before the conceptual reindexing).
        for i in 0..W {
            b.cx(v[256 + i], m[i]);
        }
        b.free_vec(&m);
    }

    fn set_slice_u512_by<R: sha3::digest::XofReader>(sim: &mut crate::sim::Simulator<R>, qs: &[super::super::QubitId], val: U512) {
        for (i, &q) in qs.iter().enumerate() {
            if val.bit(i) {
                *sim.qubit_mut(q) |= 1;
            } else {
                *sim.qubit_mut(q) &= !1;
            }
        }
    }

    fn get_slice_u512_by<R: sha3::digest::XofReader>(sim: &crate::sim::Simulator<R>, qs: &[super::super::QubitId]) -> U512 {
        let mut bytes = [0u8; 64];
        for (i, &q) in qs.iter().enumerate() {
            if (sim.qubit(q) & 1) != 0 {
                bytes[i / 8] |= 1u8 << (i % 8);
            }
        }
        U512::from_le_slice(&bytes)
    }

    #[test]
    fn approximate_batched_halve16_canonical_circuit_matches_classical() {
        let mut b = super::super::B::new();
        let v = b.alloc_qubits(274);
        emit_approx_batched_halve16_canonical(&mut b, &v);
        let num_qubits = b.next_qubit as usize;
        let num_bits = b.next_bit as usize;
        let ops = b.ops;
        let p = u256_to_u512_for_by_tests(SECP256K1_P);
        let pinv = 51_919u64;
        let mask = (1u64 << 16) - 1;
        let mut sampler = Sampler::new(b"by-batched-halve16-circuit-v1", SECP256K1_P);
        for _ in 0..64 {
            let t = sampler.next();
            let low = t.as_limbs()[0] & mask;
            let m = low.wrapping_mul((!pinv).wrapping_add(1)) & mask;
            let expected: U512 = (u256_to_u512_for_by_tests(t) + U512::from(m) * p) >> 16usize;
            let mut hasher = sha3::Shake128::default();
            hasher.update(b"by-batched-halve16-sim-xof-v1");
            let mut xof = hasher.finalize_xof();
            let mut sim = crate::sim::Simulator::new(num_qubits, num_bits, &mut xof);
            set_slice_u512_by(&mut sim, &v, u256_to_u512_for_by_tests(t));
            sim.apply(&ops);
            let got = get_slice_u512_by(&sim, &v);
            assert_eq!(got, expected, "batched halve16 circuit mismatch for T={t}");
        }
    }

    #[test]
    fn batched_halve16_top_bits_recover_correction_with_negligible_exception() {
        // Classical validation of the approximate uncompute used by the cost
        // model above. For canonical T, m = -T*p^{-1} mod 2^16. After
        // q=(T+m*p)/2^16, the top 16 bits of q equal m except when T < m*c,
        // a tiny O(2^48/p) set. That is far below the user's 1% allowance.
        let p_u = u256_to_u512_for_by_tests(SECP256K1_P);
        let modulus = 1u64 << 16;
        let pinv = 51_919u64;
        let mut failures = 0usize;
        let samples = 20_000usize;
        let mut sampler = Sampler::new(b"by-batched-halve16-topbits-v1", SECP256K1_P);
        for _ in 0..samples {
            let t = sampler.next();
            let low = t.as_limbs()[0] & (modulus - 1);
            let m = low.wrapping_mul((!pinv).wrapping_add(1)) & (modulus - 1);
            let t_u = u256_to_u512_for_by_tests(t);
            let q: U512 = (t_u + U512::from(m) * p_u) >> 16usize;
            let q_top: U512 = q >> 240usize;
            let top = q_top.to::<u64>() & (modulus - 1);
            if top != m {
                failures += 1;
            }
        }
        // Exhibit the known rare exception shape.
        let t_one = U512::from(1u64);
        let m_one = (1u64.wrapping_mul((!pinv).wrapping_add(1))) & (modulus - 1);
        let q_one: U512 = (t_one + U512::from(m_one) * p_u) >> 16usize;
        let q_one_top: U512 = q_one >> 240usize;
        let top_one = q_one_top.to::<u64>() & (modulus - 1);
        eprintln!(
            "batched halve16 top-bit correction: sample_failures={failures}/{samples}, T=1 has m={m_one}, top={top_one}"
        );
        assert_eq!(failures, 0);
        assert_ne!(top_one, m_one, "expected rare small-T exception disappeared; revisit proof");
    }

    fn emit_approx_highfold_p_for_cost(b: &mut super::super::B, v: &[super::super::QubitId]) {
        // Approximate T <- T - k*p with k = signed high bits T>>256.
        // Cost model treats k as an 18-bit magnitude/control slice; sign handling
        // would add a small constant amount and does not change the conclusion.
        assert!(v.len() >= 274);
        let k = v[256..274].to_vec();
        for &sh in &[0usize, 4, 6, 7, 8, 9, 32] {
            add_shifted_small_reg_for_cost(b, &k, v, sh, false);
        }
        add_shifted_small_reg_for_cost(b, &k, v, 256, true);
    }

    fn add_low_coeff_mod16_for_cost(
        b: &mut super::super::B,
        coeff_mod: u64,
        src: &[super::super::QubitId],
        dst: &[super::super::QubitId],
        inverse: bool,
    ) {
        assert_eq!(dst.len(), 16);
        if inverse {
            for sh in (0..16usize).rev() {
                if ((coeff_mod >> sh) & 1) != 0 {
                    add_shifted_term_for_cost(b, src, dst, sh, true);
                }
            }
        } else {
            for sh in 0..16usize {
                if ((coeff_mod >> sh) & 1) != 0 {
                    add_shifted_term_for_cost(b, src, dst, sh, false);
                }
            }
        }
    }

    fn compute_row_correction_m_from_sources(
        b: &mut super::super::B,
        coeff0: i128,
        src0: &[super::super::QubitId],
        coeff1: i128,
        src1: &[super::super::QubitId],
        m: &[super::super::QubitId],
        inverse: bool,
    ) {
        const W: u64 = 1u64 << 16;
        let neg_pinv = ((!51_919u64).wrapping_add(1)) & (W - 1);
        let c0 = ((coeff0.rem_euclid(W as i128) as u64).wrapping_mul(neg_pinv)) & (W - 1);
        let c1 = ((coeff1.rem_euclid(W as i128) as u64).wrapping_mul(neg_pinv)) & (W - 1);
        if inverse {
            add_low_coeff_mod16_for_cost(b, c1, src1, m, true);
            add_low_coeff_mod16_for_cost(b, c0, src0, m, true);
        } else {
            add_low_coeff_mod16_for_cost(b, c0, src0, m, false);
            add_low_coeff_mod16_for_cost(b, c1, src1, m, false);
        }
    }

    fn arith_shift_right_mod_width_for_test(v: U512, width: usize, shift: usize) -> U512 {
        let mut q = v >> shift;
        if v.bit(width - 1) {
            for i in (width - shift)..width {
                q.set_bit(i, true);
            }
        }
        q
    }

    fn signed_coeff_mod_width_for_test(c: i128, width: usize) -> U512 {
        if c >= 0 {
            U512::from(c as u128)
        } else {
            (U512::from(1u64) << width) - U512::from(c.unsigned_abs())
        }
    }

    fn add_signed_shifted_term_for_cost(
        b: &mut super::super::B,
        src: &[super::super::QubitId],
        dst: &[super::super::QubitId],
        shift: usize,
        subtract: bool,
    ) {
        if shift >= dst.len() {
            return;
        }
        let len = dst.len() - shift;
        let addend = b.alloc_qubits(len);
        let copy_len = src.len().min(len);
        for i in 0..copy_len {
            b.cx(src[i], addend[i]);
        }
        if !src.is_empty() {
            let sign = src[src.len() - 1];
            for i in copy_len..len {
                b.cx(sign, addend[i]);
            }
        }
        let dst_slice = dst[shift..shift + len].to_vec();
        if subtract {
            super::super::sub_nbit_qq_fast(b, &addend, &dst_slice);
        } else {
            super::super::add_nbit_qq_fast(b, &addend, &dst_slice);
        }
        if !src.is_empty() {
            let sign = src[src.len() - 1];
            for i in copy_len..len {
                b.cx(sign, addend[i]);
            }
        }
        for i in 0..copy_len {
            b.cx(src[i], addend[i]);
        }
        b.free_vec(&addend);
    }

    fn add_signed_coeff_times_for_cost(
        b: &mut super::super::B,
        coeff: i128,
        src: &[super::super::QubitId],
        dst: &[super::super::QubitId],
    ) {
        let subtract = coeff < 0;
        let mut mag = coeff.unsigned_abs();
        let mut shift = 0usize;
        while mag != 0 {
            if (mag & 1) != 0 {
                add_signed_shifted_term_for_cost(b, src, dst, shift, subtract);
            }
            mag >>= 1;
            shift += 1;
        }
    }

    fn arith_shift_right_inplace_for_cost(b: &mut super::super::B, v: &[super::super::QubitId], shift: usize) {
        let n = v.len();
        let sign = b.alloc_qubit();
        b.cx(v[n - 1], sign);
        for i in 0..(n - shift) {
            b.swap(v[i], v[i + shift]);
        }
        for i in (n - shift)..n {
            b.cx(sign, v[i]);
        }
        b.cx(v[n - 1], sign);
        b.free(sign);
    }

    fn compute_signed_q_from_m_for_matrix(
        b: &mut super::super::B,
        mtx: TransitionMatrix,
        m0: &[super::super::QubitId],
        m1: &[super::super::QubitId],
    ) -> (Vec<super::super::QubitId>, Vec<super::super::QubitId>) {
        let sgn = det_sign_pow2(mtx, 16);
        let q0 = b.alloc_qubits(34);
        let q1 = b.alloc_qubits(34);
        add_coeff_times_for_cost(b, sgn * mtx.m11, m0, &q0);
        add_coeff_times_for_cost(b, -sgn * mtx.m01, m1, &q0);
        add_coeff_times_for_cost(b, -sgn * mtx.m10, m0, &q1);
        add_coeff_times_for_cost(b, sgn * mtx.m00, m1, &q1);
        arith_shift_right_inplace_for_cost(b, &q0, 16);
        arith_shift_right_inplace_for_cost(b, &q1, 16);
        for i in 18..34 {
            b.cx(q0[17], q0[i]);
            b.cx(q1[17], q1[i]);
        }
        let q0_live = q0[..18].to_vec();
        let q1_live = q1[..18].to_vec();
        b.free_vec(&q0[18..]);
        b.free_vec(&q1[18..]);
        (q0_live, q1_live)
    }

    fn subtract_signed_q_times_solinas_c_for_cost(
        b: &mut super::super::B,
        q: &[super::super::QubitId],
        x: &[super::super::QubitId],
    ) {
        for &sh in &[0usize, 4, 6, 7, 8, 9, 32] {
            add_signed_shifted_term_for_cost(b, q, x, sh, true);
        }
    }

    fn clear_signed_q_from_z_high_for_cost(
        b: &mut super::super::B,
        q: &[super::super::QubitId],
        z: &[super::super::QubitId],
    ) {
        for i in 18..q.len() {
            b.cx(q[17], q[i]);
        }
        for i in 0..18 {
            b.cx(z[256 + i], q[i]);
        }
    }

    fn emit_signed_row_scaled_from_sources_with_m_for_test(
        b: &mut super::super::B,
        coeff0: i128,
        src0: &[super::super::QubitId],
        coeff1: i128,
        src1: &[super::super::QubitId],
        out: &[super::super::QubitId],
        m: &[super::super::QubitId],
    ) {
        add_coeff_times_for_cost(b, coeff0, src0, out);
        add_coeff_times_for_cost(b, coeff1, src1, out);
        for &sh in &[0usize, 4, 6, 7, 8, 9, 32] {
            add_shifted_small_reg_for_cost(b, m, out, sh, true);
        }
        add_shifted_small_reg_for_cost(b, m, out, 256, false);
        let sign = b.alloc_qubit();
        b.cx(out[out.len() - 1], sign);
        for i in 0..(out.len() - 16) {
            b.swap(out[i], out[i + 16]);
        }
        for i in (out.len() - 16)..out.len() {
            b.cx(sign, out[i]);
        }
        b.cx(out[out.len() - 1], sign);
        b.free(sign);
    }

    fn emit_signed_row_scaled_from_sources_for_test(
        b: &mut super::super::B,
        coeff0: i128,
        src0: &[super::super::QubitId],
        coeff1: i128,
        src1: &[super::super::QubitId],
        out: &[super::super::QubitId],
    ) {
        let m = b.alloc_qubits(16);
        compute_row_correction_m_from_sources(b, coeff0, src0, coeff1, src1, &m, false);
        emit_signed_row_scaled_from_sources_with_m_for_test(b, coeff0, src0, coeff1, src1, out, &m);
        compute_row_correction_m_from_sources(b, coeff0, src0, coeff1, src1, &m, true);
        b.free_vec(&m);
    }

    fn emit_positive_row_scaled_from_sources_for_test(
        b: &mut super::super::B,
        coeff0: i128,
        src0: &[super::super::QubitId],
        coeff1: i128,
        src1: &[super::super::QubitId],
        out: &[super::super::QubitId],
    ) {
        add_coeff_times_for_cost(b, coeff0, src0, out);
        add_coeff_times_for_cost(b, coeff1, src1, out);
        let m = b.alloc_qubits(16);
        compute_row_correction_m_from_sources(b, coeff0, src0, coeff1, src1, &m, false);
        for &sh in &[0usize, 4, 6, 7, 8, 9, 32] {
            add_shifted_small_reg_for_cost(b, &m, out, sh, true);
        }
        add_shifted_small_reg_for_cost(b, &m, out, 256, false);
        for i in 0..(out.len() - 16) {
            b.swap(out[i], out[i + 16]);
        }
        compute_row_correction_m_from_sources(b, coeff0, src0, coeff1, src1, &m, true);
        b.free_vec(&m);
    }

    #[test]
    fn signed_matrix_forward_rows_clean_m_and_match_twos_complement() {
        const WIDTH: usize = 274;
        let mtx = jump_matrix_direct_lowword(16, 16, 1, 1, 3).3;
        assert!(mtx.m00 < 0 || mtx.m01 < 0 || mtx.m10 < 0 || mtx.m11 < 0, "sample matrix should exercise signs: {:?}", mtx);
        let mut b = super::super::B::new();
        let x0 = b.alloc_qubits(256);
        let x1 = b.alloc_qubits(256);
        let y0 = b.alloc_qubits(WIDTH);
        let y1 = b.alloc_qubits(WIDTH);
        emit_signed_row_scaled_from_sources_for_test(&mut b, mtx.m00, &x0, mtx.m01, &x1, &y0);
        emit_signed_row_scaled_from_sources_for_test(&mut b, mtx.m10, &x0, mtx.m11, &x1, &y1);
        let ccx = count_ccx(&b.ops);
        let peak = b.peak_qubits;
        let num_qubits = b.next_qubit as usize;
        let num_bits = b.next_bit as usize;
        let ops = b.ops;
        let width_mod = U512::from(1u64) << WIDTH;
        let width_mask = width_mod - U512::from(1u64);
        let p512 = u256_to_u512_for_by_tests(SECP256K1_P);
        let pinv = 51_919u64;
        let low_mask = (1u64 << 16) - 1;
        let mut sx = Sampler::new(b"by-signed-forward-row-x0-v1", SECP256K1_P);
        let mut sy = Sampler::new(b"by-signed-forward-row-x1-v1", SECP256K1_P);
        for _ in 0..32 {
            let a = sx.next();
            let c = sy.next();
            let x0w = u256_to_u512_for_by_tests(a);
            let x1w = u256_to_u512_for_by_tests(c);
            let expected_row = |c0: i128, c1: i128| -> U512 {
                let t = (x0w * signed_coeff_mod_width_for_test(c0, WIDTH)
                    + x1w * signed_coeff_mod_width_for_test(c1, WIDTH)) & width_mask;
                let corr = (t.as_limbs()[0] & low_mask).wrapping_mul((!pinv).wrapping_add(1)) & low_mask;
                let v = (t + U512::from(corr) * p512) & width_mask;
                arith_shift_right_mod_width_for_test(v, WIDTH, 16)
            };
            let exp0 = expected_row(mtx.m00, mtx.m01);
            let exp1 = expected_row(mtx.m10, mtx.m11);
            let mut hasher = sha3::Shake128::default();
            hasher.update(b"by-signed-forward-row-sim-v1");
            let mut xof = hasher.finalize_xof();
            let mut sim = crate::sim::Simulator::new(num_qubits, num_bits, &mut xof);
            set_slice_u512_by(&mut sim, &x0, x0w);
            set_slice_u512_by(&mut sim, &x1, x1w);
            sim.apply(&ops);
            assert_eq!(get_slice_u512_by(&sim, &y0), exp0, "signed row0 mismatch for {:?}", mtx);
            assert_eq!(get_slice_u512_by(&sim, &y1), exp1, "signed row1 mismatch for {:?}", mtx);
        }
        eprintln!("signed BY matrix forward rows: ccx={ccx}, peak={peak}q, matrix={:?}", mtx);
        assert!(ccx < 35_000, "signed forward rows too costly");
        assert!(peak < 2_200, "signed forward row peak too high");
    }

    #[test]
    fn row_correction_m_from_sources_circuit_matches_classical() {
        let mut b = super::super::B::new();
        let x0 = b.alloc_qubits(256);
        let x1 = b.alloc_qubits(256);
        let m = b.alloc_qubits(16);
        compute_row_correction_m_from_sources(&mut b, 65535, &x0, 1, &x1, &m, false);
        let ops = b.ops;
        let num_qubits = b.next_qubit as usize;
        let num_bits = b.next_bit as usize;
        let a = U256::from_limbs([
            0x7f7df51fc0ad69fa,
            0x79422d087c39ea56,
            0x00a59c1897e6d50a,
            0xfc2ad18cfe76cc7f,
        ]) % SECP256K1_P;
        let c = U256::from_limbs([
            0x96e72f29e7c30894,
            0x4ae30ac8953f8e71,
            0xc9ab887a528b640a,
            0x9d92bbd5d05a25ba,
        ]) % SECP256K1_P;
        let pinv = 51_919u64;
        let low = ((a.as_limbs()[0].wrapping_neg()).wrapping_add(c.as_limbs()[0])) & 0xffff;
        let expected = low.wrapping_mul((!pinv).wrapping_add(1)) & 0xffff;
        let mut hasher = sha3::Shake128::default();
        hasher.update(b"by-row-correction-m-sim-v1");
        let mut xof = hasher.finalize_xof();
        let mut sim = crate::sim::Simulator::new(num_qubits, num_bits, &mut xof);
        set_slice_u512_by(&mut sim, &x0, u256_to_u512_for_by_tests(a));
        set_slice_u512_by(&mut sim, &x1, u256_to_u512_for_by_tests(c));
        sim.apply(&ops);
        let got = get_slice_u512_by(&sim, &m).to::<u64>();
        assert_eq!(got, expected, "m mismatch got={got} exp={expected}");
    }

    #[test]
    fn fixed_positive_matrix_forward_rows_clean_m_and_match_classical() {
        // First actual noncanonical row circuit: m is computed from the original
        // row sources and uncomputed from those same sources after the shift,
        // so no top-bit recovery or quotient history is needed for the forward
        // rows. This is only the forward half for a positive sampled matrix;
        // old-row cleanup is still the open problem.
        const WIDTH: usize = 274;
        let mtx = jump_matrix_direct_lowword(16, 16, -20, 1, 1).3;
        assert_eq!((mtx.m00, mtx.m01, mtx.m10, mtx.m11), (65536, 0, 65535, 1));
        let mut b = super::super::B::new();
        let x0 = b.alloc_qubits(256);
        let x1 = b.alloc_qubits(256);
        let y0 = b.alloc_qubits(WIDTH);
        let y1 = b.alloc_qubits(WIDTH);
        emit_positive_row_scaled_from_sources_for_test(&mut b, mtx.m00, &x0, mtx.m01, &x1, &y0);
        emit_positive_row_scaled_from_sources_for_test(&mut b, mtx.m10, &x0, mtx.m11, &x1, &y1);
        let ccx = count_ccx(&b.ops);
        let num_qubits = b.next_qubit as usize;
        let num_bits = b.next_bit as usize;
        let ops = b.ops;
        let p512 = u256_to_u512_for_by_tests(SECP256K1_P);
        let pinv = 51_919u64;
        let mask = (1u64 << 16) - 1;
        let mut sx = Sampler::new(b"by-fixed-positive-row-x0-v1", SECP256K1_P);
        let mut sy = Sampler::new(b"by-fixed-positive-row-x1-v1", SECP256K1_P);
        for _ in 0..32 {
            let a = sx.next();
            let c = sy.next();
            let t0 = u256_to_u512_for_by_tests(a) * U512::from(mtx.m00 as u128);
            let low0 = t0.as_limbs()[0] & mask;
            let corr0 = low0.wrapping_mul((!pinv).wrapping_add(1)) & mask;
            let exp0: U512 = (t0 + U512::from(corr0) * p512) >> 16usize;
            let t1 = u256_to_u512_for_by_tests(a) * U512::from(mtx.m10 as u128)
                + u256_to_u512_for_by_tests(c) * U512::from(mtx.m11 as u128);
            let low1 = t1.as_limbs()[0] & mask;
            let corr1 = low1.wrapping_mul((!pinv).wrapping_add(1)) & mask;
            let exp1: U512 = (t1 + U512::from(corr1) * p512) >> 16usize;
            let mut hasher = sha3::Shake128::default();
            hasher.update(b"by-fixed-positive-row-sim-v1");
            let mut xof = hasher.finalize_xof();
            let mut sim = crate::sim::Simulator::new(num_qubits, num_bits, &mut xof);
            set_slice_u512_by(&mut sim, &x0, u256_to_u512_for_by_tests(a));
            set_slice_u512_by(&mut sim, &x1, u256_to_u512_for_by_tests(c));
            sim.apply(&ops);
            let got0 = get_slice_u512_by(&sim, &y0);
            let got1 = get_slice_u512_by(&sim, &y1);
            assert_eq!(got0, exp0, "row0 mismatch a={a} got={got0} exp={exp0}");
            assert_eq!(got1, exp1, "row1 mismatch a={a} c={c} got={got1} exp={exp1}");
        }
        eprintln!(
            "fixed positive BY matrix forward rows: ccx={ccx}, peak={}q, matrix={:?}",
            num_qubits, mtx
        );
        assert!(ccx < 20_000, "forward rows too expensive for BY window budget");
    }

    fn emit_positive_triangular_old_cleanup_for_test(
        b: &mut super::super::B,
        x0: &[super::super::QubitId],
        x1: &[super::super::QubitId],
        y0: &[super::super::QubitId],
        y1: &[super::super::QubitId],
    ) -> (Vec<super::super::QubitId>, Vec<super::super::QubitId>) {
        // Matrix [[2^16,0],[2^16-1,1]]. Outputs satisfy:
        //   y0 = x0
        //   2^16*y1 = (2^16-1)x0 + x1 + m*p
        // Therefore z = 2^16*y1 - (2^16-1)y0 = x1 + m*p.
        // To zero x1, subtract z low bits (leaving m*c) and then subtract m*c.
        let m = b.alloc_qubits(16);
        compute_row_correction_m_from_sources(b, 65535, x0, 1, x1, &m, false);
        let z = b.alloc_qubits(274);
        add_coeff_times_for_cost(b, 65536, y1, &z);
        add_coeff_times_for_cost(b, -65535, y0, &z);
        let z_low = z[..256].to_vec();
        super::super::sub_nbit_qq_fast(b, &z_low, x1);
        for &sh in &[0usize, 4, 6, 7, 8, 9, 32] {
            add_shifted_small_reg_for_cost(b, &m, x1, sh, true);
        }
        // Approximate m cleanup from z's high bits. For z=x1+m*p with x1<p,
        // top bits equal m except the same tiny x1<m*c exception as before.
        for i in 0..16 {
            b.cx(z[256 + i], m[i]);
        }
        // Uncompute z from y.
        add_coeff_times_for_cost(b, 65535, y0, &z);
        add_coeff_times_for_cost(b, -65536, y1, &z);
        // x0 cleanup is exact: y0=x0 for this triangular matrix.
        let y0_low = y0[..256].to_vec();
        super::super::sub_nbit_qq_fast(b, &y0_low, x0);
        (m, z)
    }

    fn emit_fixed_matrix_old_cleanup_with_mq_for_test(
        b: &mut super::super::B,
        mtx: TransitionMatrix,
        x0: &[super::super::QubitId],
        x1: &[super::super::QubitId],
        y0: &[super::super::QubitId],
        y1: &[super::super::QubitId],
        m0: &[super::super::QubitId],
        m1: &[super::super::QubitId],
        q0: &[super::super::QubitId],
        q1: &[super::super::QubitId],
    ) -> (Vec<super::super::QubitId>, Vec<super::super::QubitId>) {
        let sgn = det_sign_pow2(mtx, 16);
        let z0 = b.alloc_qubits(274);
        let z1 = b.alloc_qubits(274);
        add_signed_coeff_times_for_cost(b, sgn * mtx.m11, y0, &z0);
        add_signed_coeff_times_for_cost(b, -sgn * mtx.m01, y1, &z0);
        add_signed_coeff_times_for_cost(b, -sgn * mtx.m10, y0, &z1);
        add_signed_coeff_times_for_cost(b, sgn * mtx.m00, y1, &z1);

        let z0_low = z0[..256].to_vec();
        let z1_low = z1[..256].to_vec();
        super::super::sub_nbit_qq_fast(b, &z0_low, x0);
        super::super::sub_nbit_qq_fast(b, &z1_low, x1);
        subtract_signed_q_times_solinas_c_for_cost(b, q0, x0);
        subtract_signed_q_times_solinas_c_for_cost(b, q1, x1);

        // Clear m using P*q = m (mod 2^16).
        add_low_coeff_mod16_for_cost(b, mtx.m00.rem_euclid(1 << 16) as u64, q0, m0, true);
        add_low_coeff_mod16_for_cost(b, mtx.m01.rem_euclid(1 << 16) as u64, q1, m0, true);
        add_low_coeff_mod16_for_cost(b, mtx.m10.rem_euclid(1 << 16) as u64, q0, m1, true);
        add_low_coeff_mod16_for_cost(b, mtx.m11.rem_euclid(1 << 16) as u64, q1, m1, true);

        clear_signed_q_from_z_high_for_cost(b, q0, &z0);
        clear_signed_q_from_z_high_for_cost(b, q1, &z1);

        add_signed_coeff_times_for_cost(b, -sgn * mtx.m11, y0, &z0);
        add_signed_coeff_times_for_cost(b, sgn * mtx.m01, y1, &z0);
        add_signed_coeff_times_for_cost(b, sgn * mtx.m10, y0, &z1);
        add_signed_coeff_times_for_cost(b, -sgn * mtx.m00, y1, &z1);
        (z0, z1)
    }

    fn emit_fixed_matrix_old_cleanup_for_test(
        b: &mut super::super::B,
        mtx: TransitionMatrix,
        x0: &[super::super::QubitId],
        x1: &[super::super::QubitId],
        y0: &[super::super::QubitId],
        y1: &[super::super::QubitId],
    ) -> (
        Vec<super::super::QubitId>,
        Vec<super::super::QubitId>,
        Vec<super::super::QubitId>,
        Vec<super::super::QubitId>,
        Vec<super::super::QubitId>,
        Vec<super::super::QubitId>,
    ) {
        let m0 = b.alloc_qubits(16);
        let m1 = b.alloc_qubits(16);
        compute_row_correction_m_from_sources(b, mtx.m00, x0, mtx.m01, x1, &m0, false);
        compute_row_correction_m_from_sources(b, mtx.m10, x0, mtx.m11, x1, &m1, false);
        let (q0, q1) = compute_signed_q_from_m_for_matrix(b, mtx, &m0, &m1);
        let (z0, z1) = emit_fixed_matrix_old_cleanup_with_mq_for_test(b, mtx, x0, x1, y0, y1, &m0, &m1, &q0, &q1);
        (m0, m1, q0, q1, z0, z1)
    }

    #[test]
    fn signed_sample_fixed_matrix_replacement_cleans_old_rows() {
        const WIDTH: usize = 274;
        let mtx = jump_matrix_direct_lowword(16, 16, 1, 1, 3).3;
        assert_eq!((mtx.m00, mtx.m01, mtx.m10, mtx.m11), (-8192, 24576, -3, 1));
        let mut b = super::super::B::new();
        let x0 = b.alloc_qubits(256);
        let x1 = b.alloc_qubits(256);
        let y0 = b.alloc_qubits(WIDTH);
        let y1 = b.alloc_qubits(WIDTH);
        emit_signed_row_scaled_from_sources_for_test(&mut b, mtx.m00, &x0, mtx.m01, &x1, &y0);
        emit_signed_row_scaled_from_sources_for_test(&mut b, mtx.m10, &x0, mtx.m11, &x1, &y1);
        let (m0, m1, q0, q1, z0, z1) = emit_fixed_matrix_old_cleanup_for_test(&mut b, mtx, &x0, &x1, &y0, &y1);
        let ccx = count_ccx(&b.ops);
        let peak = b.peak_qubits;
        let num_qubits = b.next_qubit as usize;
        let num_bits = b.next_bit as usize;
        let ops = b.ops;
        let width_mod = U512::from(1u64) << WIDTH;
        let width_mask = width_mod - U512::from(1u64);
        let p512 = u256_to_u512_for_by_tests(SECP256K1_P);
        let pinv = 51_919u64;
        let low_mask = (1u64 << 16) - 1;
        let mut sx = Sampler::new(b"by-signed-repl-x0-v1", SECP256K1_P);
        let mut sy = Sampler::new(b"by-signed-repl-x1-v1", SECP256K1_P);
        for _ in 0..32 {
            let a = sx.next();
            let c = sy.next();
            let x0w = u256_to_u512_for_by_tests(a);
            let x1w = u256_to_u512_for_by_tests(c);
            let expected_row = |c0: i128, c1: i128| -> U512 {
                let t = (x0w * signed_coeff_mod_width_for_test(c0, WIDTH)
                    + x1w * signed_coeff_mod_width_for_test(c1, WIDTH)) & width_mask;
                let corr = (t.as_limbs()[0] & low_mask).wrapping_mul((!pinv).wrapping_add(1)) & low_mask;
                let v = (t + U512::from(corr) * p512) & width_mask;
                arith_shift_right_mod_width_for_test(v, WIDTH, 16)
            };
            let exp0 = expected_row(mtx.m00, mtx.m01);
            let exp1 = expected_row(mtx.m10, mtx.m11);
            let mut hasher = sha3::Shake128::default();
            hasher.update(b"by-signed-repl-sim-v1");
            let mut xof = hasher.finalize_xof();
            let mut sim = crate::sim::Simulator::new(num_qubits, num_bits, &mut xof);
            set_slice_u512_by(&mut sim, &x0, x0w);
            set_slice_u512_by(&mut sim, &x1, x1w);
            sim.apply(&ops);
            assert_eq!(get_slice_u512_by(&sim, &x0), U512::ZERO, "x0 not zero");
            assert_eq!(get_slice_u512_by(&sim, &x1), U512::ZERO, "x1 not zero");
            assert_eq!(get_slice_u512_by(&sim, &m0), U512::ZERO, "m0 not zero");
            assert_eq!(get_slice_u512_by(&sim, &m1), U512::ZERO, "m1 not zero");
            assert_eq!(get_slice_u512_by(&sim, &q0), U512::ZERO, "q0 not zero");
            assert_eq!(get_slice_u512_by(&sim, &q1), U512::ZERO, "q1 not zero");
            assert_eq!(get_slice_u512_by(&sim, &z0), U512::ZERO, "z0 not zero");
            assert_eq!(get_slice_u512_by(&sim, &z1), U512::ZERO, "z1 not zero");
            assert_eq!(get_slice_u512_by(&sim, &y0), exp0, "y0 mismatch");
            assert_eq!(get_slice_u512_by(&sim, &y1), exp1, "y1 mismatch");
        }
        eprintln!("signed sample BY fixed-matrix replacement: ccx={ccx}, peak={peak}q");
        assert!(ccx < 45_000, "signed fixed-matrix replacement too costly");
        assert!(peak < 2_700, "signed fixed-matrix replacement peak too high");
    }

    #[test]
    fn fixed_matrix_replacement_sample_cost_distribution() {
        const WIDTH: usize = 274;
        let mut hasher = sha3::Shake128::default();
        hasher.update(b"by-fixed-matrix-replacement-cost-dist-v1");
        let mut reader = hasher.finalize_xof();
        let mut buf = [0u8; 24];
        let samples = 32usize;
        let mut costs = Vec::with_capacity(samples);
        let mut peaks = Vec::with_capacity(samples);
        for _ in 0..samples {
            reader.read(&mut buf);
            let f_low = (u64::from_le_bytes(buf[0..8].try_into().unwrap()) as i128) | 1;
            let g_low = u64::from_le_bytes(buf[8..16].try_into().unwrap()) as i128;
            let delta = (u64::from_le_bytes(buf[16..24].try_into().unwrap()) % 41) as i64 - 20;
            let (_, _, _, mtx) = jump_matrix_direct_lowword(16, 16, delta, f_low, g_low);
            let mut b = super::super::B::new();
            let x0 = b.alloc_qubits(256);
            let x1 = b.alloc_qubits(256);
            let y0 = b.alloc_qubits(WIDTH);
            let y1 = b.alloc_qubits(WIDTH);
            emit_signed_row_scaled_from_sources_for_test(&mut b, mtx.m00, &x0, mtx.m01, &x1, &y0);
            emit_signed_row_scaled_from_sources_for_test(&mut b, mtx.m10, &x0, mtx.m11, &x1, &y1);
            let _regs = emit_fixed_matrix_old_cleanup_for_test(&mut b, mtx, &x0, &x1, &y0, &y1);
            costs.push(count_ccx(&b.ops));
            peaks.push(b.peak_qubits);
        }
        costs.sort_unstable();
        peaks.sort_unstable();
        let mean_cost = costs.iter().sum::<usize>() as f64 / samples as f64;
        let p90_cost = costs[(samples * 90) / 100];
        let max_cost = costs[samples - 1];
        let p90_peak = peaks[(samples * 90) / 100];
        let max_peak = peaks[samples - 1];
        eprintln!(
            "BY fixed-matrix replacement cost distribution: mean_ccx={mean_cost:.1}, p90_ccx={p90_cost}, max_ccx={max_cost}, p90_peak={p90_peak}q, max_peak={max_peak}q"
        );
        assert!(p90_cost < 45_000, "fixed-matrix replacement p90 too costly");
        assert!(max_peak < 2_800, "fixed-matrix replacement sample exceeds cap");
    }

    #[test]
    fn fixed_matrix_window_with_free_mq_history_still_misses_target() {
        // Grant the selected-window denominator route a very generous oracle:
        // the 16-bit row corrections m0/m1 and 18-bit adjugate quotients q0/q1
        // are already present as clean history bits.  The window only has to
        // form the two scaled rows, consume the old rows, and clear those m/q
        // bits from the output-side residuals.  If even this misses ~10k/window,
        // then the current fixed-matrix arithmetic cannot be rescued merely by
        // making the selector oracle smarter.
        const WIDTH: usize = 274;
        const SAMPLES: usize = 24;
        const W: i128 = 1 << 16;
        let p_mod = SECP256K1_P;
        let p512 = u256_to_u512_for_by_tests(p_mod);
        let pinv = 51_919u64;
        let neg_pinv = ((!pinv).wrapping_add(1)) & ((1u64 << 16) - 1);
        let low_mask = (1u64 << 16) - 1;
        let width_mod = U512::from(1u64) << WIDTH;
        let width_mask = width_mod - U512::from(1u64);
        let row_m = |c0: i128, x0: U256, c1: i128, x1: U256| -> u64 {
            let a = (c0.rem_euclid(W) as u64).wrapping_mul(x0.as_limbs()[0] & low_mask);
            let b = (c1.rem_euclid(W) as u64).wrapping_mul(x1.as_limbs()[0] & low_mask);
            a.wrapping_add(b).wrapping_mul(neg_pinv) & low_mask
        };
        let row_expected = |c0: i128, x0w: U512, c1: i128, x1w: U512| -> U512 {
            let t = (x0w * signed_coeff_mod_width_for_test(c0, WIDTH)
                + x1w * signed_coeff_mod_width_for_test(c1, WIDTH)) & width_mask;
            let corr = (t.as_limbs()[0] & low_mask).wrapping_mul(neg_pinv) & low_mask;
            let v = (t + U512::from(corr) * p512) & width_mask;
            arith_shift_right_mod_width_for_test(v, WIDTH, 16)
        };
        let q_pair = |mtx: TransitionMatrix, m0: u64, m1: u64| -> (i128, i128) {
            let sgn = det_sign_pow2(mtx, 16);
            let q0_num = sgn * mtx.m11 * m0 as i128 - sgn * mtx.m01 * m1 as i128;
            let q1_num = -sgn * mtx.m10 * m0 as i128 + sgn * mtx.m00 * m1 as i128;
            assert_eq!(q0_num % W, 0, "q0 integrality failed for {mtx:?}");
            assert_eq!(q1_num % W, 0, "q1 integrality failed for {mtx:?}");
            let q0 = q0_num / W;
            let q1 = q1_num / W;
            assert!((-131_072..131_072).contains(&q0), "q0 needs more than 18 bits: {q0}");
            assert!((-131_072..131_072).contains(&q1), "q1 needs more than 18 bits: {q1}");
            (q0, q1)
        };

        let mut hasher = sha3::Shake128::default();
        hasher.update(b"by-free-mq-window-budget-v1");
        let mut reader = hasher.finalize_xof();
        let mut buf = [0u8; 24];
        let mut costs = Vec::with_capacity(SAMPLES);
        let mut peaks = Vec::with_capacity(SAMPLES);
        for sample_idx in 0..SAMPLES {
            reader.read(&mut buf);
            let f_low = (u64::from_le_bytes(buf[0..8].try_into().unwrap()) as i128) | 1;
            let g_low = u64::from_le_bytes(buf[8..16].try_into().unwrap()) as i128;
            let delta = (u64::from_le_bytes(buf[16..24].try_into().unwrap()) % 41) as i64 - 20;
            let (_, _, _, mtx) = jump_matrix_direct_lowword(16, 16, delta, f_low, g_low);
            let mut b = super::super::B::new();
            let x0 = b.alloc_qubits(256);
            let x1 = b.alloc_qubits(256);
            let y0 = b.alloc_qubits(WIDTH);
            let y1 = b.alloc_qubits(WIDTH);
            let m0 = b.alloc_qubits(16);
            let m1 = b.alloc_qubits(16);
            let q0 = b.alloc_qubits(18);
            let q1 = b.alloc_qubits(18);
            emit_signed_row_scaled_from_sources_with_m_for_test(&mut b, mtx.m00, &x0, mtx.m01, &x1, &y0, &m0);
            emit_signed_row_scaled_from_sources_with_m_for_test(&mut b, mtx.m10, &x0, mtx.m11, &x1, &y1, &m1);
            let (z0, z1) = emit_fixed_matrix_old_cleanup_with_mq_for_test(
                &mut b, mtx, &x0, &x1, &y0, &y1, &m0, &m1, &q0, &q1,
            );
            let ccx = count_ccx(&b.ops);
            costs.push(ccx);
            peaks.push(b.peak_qubits);

            if sample_idx < 3 {
                let num_qubits = b.next_qubit as usize;
                let num_bits = b.next_bit as usize;
                let ops = b.ops.clone();
                let mut sx = Sampler::new(b"by-free-mq-x0-v1", p_mod);
                let mut sy = Sampler::new(b"by-free-mq-x1-v1", p_mod);
                for _ in 0..8 {
                    let a = sx.next();
                    let c = sy.next();
                    let m0v = row_m(mtx.m00, a, mtx.m01, c);
                    let m1v = row_m(mtx.m10, a, mtx.m11, c);
                    let (q0v, q1v) = q_pair(mtx, m0v, m1v);
                    let x0w = u256_to_u512_for_by_tests(a);
                    let x1w = u256_to_u512_for_by_tests(c);
                    let exp0 = row_expected(mtx.m00, x0w, mtx.m01, x1w);
                    let exp1 = row_expected(mtx.m10, x0w, mtx.m11, x1w);
                    let mut h = sha3::Shake128::default();
                    h.update(b"by-free-mq-window-sim-v1");
                    let mut xof = h.finalize_xof();
                    let mut sim = crate::sim::Simulator::new(num_qubits, num_bits, &mut xof);
                    set_slice_u512_by(&mut sim, &x0, x0w);
                    set_slice_u512_by(&mut sim, &x1, x1w);
                    set_slice_u512_by(&mut sim, &m0, U512::from(m0v));
                    set_slice_u512_by(&mut sim, &m1, U512::from(m1v));
                    set_slice_u512_by(&mut sim, &q0, twos_u512_for_delta(q0v as i64, 18));
                    set_slice_u512_by(&mut sim, &q1, twos_u512_for_delta(q1v as i64, 18));
                    sim.apply(&ops);
                    assert_eq!(get_slice_u512_by(&sim, &x0), U512::ZERO, "x0 not zero");
                    assert_eq!(get_slice_u512_by(&sim, &x1), U512::ZERO, "x1 not zero");
                    assert_eq!(get_slice_u512_by(&sim, &m0), U512::ZERO, "m0 not zero");
                    assert_eq!(get_slice_u512_by(&sim, &m1), U512::ZERO, "m1 not zero");
                    assert_eq!(get_slice_u512_by(&sim, &q0), U512::ZERO, "q0 not zero");
                    assert_eq!(get_slice_u512_by(&sim, &q1), U512::ZERO, "q1 not zero");
                    assert_eq!(get_slice_u512_by(&sim, &z0), U512::ZERO, "z0 not zero");
                    assert_eq!(get_slice_u512_by(&sim, &z1), U512::ZERO, "z1 not zero");
                    assert_eq!(get_slice_u512_by(&sim, &y0), exp0, "y0 mismatch");
                    assert_eq!(get_slice_u512_by(&sim, &y1), exp1, "y1 mismatch");
                }
            }
        }
        costs.sort_unstable();
        peaks.sort_unstable();
        let mean = costs.iter().sum::<usize>() as f64 / SAMPLES as f64;
        let p90 = costs[(SAMPLES * 90) / 100];
        let max = costs[SAMPLES - 1];
        let max_peak = peaks[SAMPLES - 1];
        let target = 10_000.0;
        let gap = mean - target;
        println!("METRIC by_fixed_window_free_mq_mean_ccx={mean:.3}");
        println!("METRIC by_fixed_window_free_mq_p90_ccx={p90}");
        println!("METRIC by_fixed_window_free_mq_max_ccx={max}");
        println!("METRIC by_fixed_window_free_mq_gap_to_target_ccx={gap:.3}");
        println!("METRIC by_fixed_window_free_mq_max_peak={max_peak}");
        eprintln!(
            "BY fixed-matrix window with free m/q history: mean_ccx={mean:.1}, p90_ccx={p90}, max_ccx={max}, gap_to_10k={gap:.1}, max_peak={max_peak}q"
        );
        assert!(mean > target * 1.15, "free m/q window unexpectedly meets the selected-window target");
    }

    #[test]
    fn last_shot_fixed_matrix_window_consumption_misses_sota_budget() {
        // Final BY SOTA gate: after product-clean replay, the only credible
        // remaining denominator path is a 16-step selected fixed-matrix/q
        // consumption update.  Measure the actual reversible one-window object
        // we currently know how to synthesize (form scaled rows, clean old rows,
        // clear m/q/z).  A SOTA-shaped two-denominator plan needs roughly
        // <10k CCX/window; this existing arithmetic is about 2× too expensive,
        // so wiring it into point-add would be predictably late/dead.
        const WIDTH: usize = 274;
        const SAMPLES: usize = 24;
        const WINDOWS: usize = 36; // 576-step exact setting used by the harness
        let mut hasher = sha3::Shake128::default();
        hasher.update(b"by-last-shot-fixed-window-budget-v1");
        let mut reader = hasher.finalize_xof();
        let mut buf = [0u8; 24];
        let mut costs = Vec::with_capacity(SAMPLES);
        let mut max_peak = 0u32;
        for _ in 0..SAMPLES {
            reader.read(&mut buf);
            let f_low = (u64::from_le_bytes(buf[0..8].try_into().unwrap()) as i128) | 1;
            let g_low = u64::from_le_bytes(buf[8..16].try_into().unwrap()) as i128;
            let delta = (u64::from_le_bytes(buf[16..24].try_into().unwrap()) % 41) as i64 - 20;
            let (_, _, _, mtx) = jump_matrix_direct_lowword(16, 16, delta, f_low, g_low);
            let mut b = super::super::B::new();
            let x0 = b.alloc_qubits(256);
            let x1 = b.alloc_qubits(256);
            let y0 = b.alloc_qubits(WIDTH);
            let y1 = b.alloc_qubits(WIDTH);
            emit_signed_row_scaled_from_sources_for_test(&mut b, mtx.m00, &x0, mtx.m01, &x1, &y0);
            emit_signed_row_scaled_from_sources_for_test(&mut b, mtx.m10, &x0, mtx.m11, &x1, &y1);
            let _regs = emit_fixed_matrix_old_cleanup_for_test(&mut b, mtx, &x0, &x1, &y0, &y1);
            costs.push(count_ccx(&b.ops));
            max_peak = max_peak.max(b.peak_qubits);
        }
        costs.sort_unstable();
        let mean = costs.iter().sum::<usize>() as f64 / SAMPLES as f64;
        let p90 = costs[(SAMPLES * 90) / 100];
        let max = costs[SAMPLES - 1];
        let two_denominators = 2.0 * WINDOWS as f64 * mean;
        let target_per_window = 10_000.0;
        eprintln!(
            "BY last-shot fixed-matrix/q window budget: mean_ccx={mean:.1}, p90_ccx={p90}, max_ccx={max}, max_peak={max_peak}q, two_denominators≈{two_denominators:.0}, target_per_window≈{target_per_window:.0}"
        );
        assert!(mean > target_per_window * 1.8, "fixed-matrix window unexpectedly reached SOTA target; wire it immediately");
        assert!(two_denominators > 1_300_000.0, "two-denominator fixed-window budget unexpectedly fits SOTA margin");
    }

    #[test]
    fn by_fixed_window_three_million_budget_still_misses_with_current_arithmetic() {
        // Relaxing the gate target from 2.7M to 3.0M does not make the known
        // BY fixed-window denominator object safe.  With two centered 873.6k
        // replays and the deleted-Kaliski scaffold, the two denominator windows
        // can spend about 13.4k CCX/window.  The current fixed-window arithmetic
        // misses even when m/q history is granted for free.
        let target_3m = 3_000_000.0;
        let two_replay_scaffold = 285_766.0;
        let centered_replay = 873_600.0;
        let windows = 36.0;
        let budget_per_window = (target_3m - two_replay_scaffold - 2.0 * centered_replay) / (2.0 * windows);
        let free_mq_mean = 18_138.583;
        let last_shot_mean = 20_323.5;
        let free_mq_projected = two_replay_scaffold + 2.0 * centered_replay + 2.0 * windows * free_mq_mean;
        let last_shot_projected = two_replay_scaffold + 2.0 * centered_replay + 2.0 * windows * last_shot_mean;
        let free_mq_gap = free_mq_projected - target_3m;
        let last_shot_gap = last_shot_projected - target_3m;
        let needed_free_mq_reduction = free_mq_mean - budget_per_window;
        println!("METRIC by_3m_fixed_window_budget_per_window_ccx={budget_per_window:.3}");
        println!("METRIC by_3m_free_mq_projected_gap_ccx={free_mq_gap:.3}");
        println!("METRIC by_3m_last_shot_projected_gap_ccx={last_shot_gap:.3}");
        println!("METRIC by_3m_free_mq_needed_reduction_per_window_ccx={needed_free_mq_reduction:.3}");
        eprintln!(
            "BY 3M fixed-window budget: budget_per_window={budget_per_window:.1}, free_mq_projected={free_mq_projected:.0} gap={free_mq_gap:.0}, last_shot_projected={last_shot_projected:.0} gap={last_shot_gap:.0}"
        );
        assert!(free_mq_gap > 300_000.0, "free m/q fixed-window arithmetic is close enough to revisit for 3M");
        assert!(needed_free_mq_reduction > 4_000.0, "only a small per-window reduction is needed; try optimizing fixed-window arithmetic");
    }

    #[test]
    fn controlled_dirty_qoffset_adder_small_basis_check() {
        let n = 8usize;
        let mask = (1u64 << n) - 1;
        let mut b = super::super::B::new();
        let ctrl = b.alloc_qubit();
        let target = b.alloc_qubits(n);
        let offset = b.alloc_qubits(n);
        let dirty = b.alloc_qubits(n - 2);
        let clean2 = [b.alloc_qubit(), b.alloc_qubit()];
        let clean3 = b.alloc_qubit();
        super::super::venting::ciadd_dirty_3clean_qoffset(&mut b, &target, &dirty, &clean2, clean3, &offset, ctrl);
        let ccx = count_ccx(&b.ops);
        let peak = b.peak_qubits;
        let num_qubits = b.next_qubit as usize;
        let num_bits = b.next_bit as usize;
        let ops = b.ops;
        for ctrl_v in [false, true] {
            for target_v in [0x00u64, 0x35, 0xf1] {
                for offset_v in [0x00u64, 0x17, 0x80] {
                    let dirty_v = 0x2du64 & ((1u64 << (n - 2)) - 1);
                    let mut hasher = sha3::Shake128::default();
                    hasher.update(b"by-ciadd-qoffset-small-v1");
                    let mut xof = hasher.finalize_xof();
                    let mut sim = crate::sim::Simulator::new(num_qubits, num_bits, &mut xof);
                    if ctrl_v {
                        *sim.qubit_mut(ctrl) |= 1;
                    }
                    set_slice_u512_by(&mut sim, &target, U512::from(target_v));
                    set_slice_u512_by(&mut sim, &offset, U512::from(offset_v));
                    set_slice_u512_by(&mut sim, &dirty, U512::from(dirty_v));
                    sim.apply(&ops);
                    let expected = if ctrl_v { target_v.wrapping_add(offset_v) & mask } else { target_v & mask };
                    assert_eq!(get_slice_u512_by(&sim, &target).to::<u64>() & mask, expected, "target mismatch");
                    assert_eq!(get_slice_u512_by(&sim, &offset).to::<u64>() & mask, offset_v & mask, "offset changed");
                    assert_eq!(get_slice_u512_by(&sim, &dirty).to::<u64>() & ((1u64 << (n - 2)) - 1), dirty_v, "dirty changed");
                    assert_eq!(sim.global_phase() & 1, 0, "phase changed");
                }
            }
        }
        eprintln!("controlled dirty qoffset small check: n={n}, ccx={ccx}, peak={peak}q");
    }

    #[test]
    fn masked_controlled_qoffset_borrows_offset_as_dirty_gate_good_scratch_short() {
        // New control trick for the dirty qoffset adder: compute a clean masked
        // offset m_i = ctrl & offset_i, then run the *uncontrolled* dirty
        // qoffset adder using the original offset register as dirty workspace.
        // This restores offset and mask, and is far cheaper than controlling
        // every vented-adder primitive.  The catch is the clean n-bit mask: gate
        // count is SOTA-shaped, but compressed-history + mask still misses the
        // user's ~600 scratch cap unless the mask can be overlapped/streamed.
        let n = 8usize;
        let maskv = (1u64 << n) - 1;
        let mut b = super::super::B::new();
        let ctrl = b.alloc_qubit();
        let target = b.alloc_qubits(n);
        let offset = b.alloc_qubits(n);
        let mask = b.alloc_qubits(n);
        let clean2 = [b.alloc_qubit(), b.alloc_qubit()];
        for k in 0..n {
            b.ccx(ctrl, offset[k], mask[k]);
        }
        super::super::venting::iadd_dirty_2clean_qoffset(&mut b, &target, &offset[..n - 2], &clean2, &mask, false);
        for k in 0..n {
            b.ccx(ctrl, offset[k], mask[k]);
        }
        let ccx8 = count_ccx(&b.ops);
        let num_qubits = b.next_qubit as usize;
        let num_bits = b.next_bit as usize;
        let ops = b.ops;
        for ctrl_v in [false, true] {
            for target_v in [0x00u64, 0x35, 0xf1] {
                for offset_v in [0x00u64, 0x17, 0x80] {
                    let mut hasher = sha3::Shake128::default();
                    hasher.update(b"by-masked-borrow-qoffset-small-v1");
                    let mut xof = hasher.finalize_xof();
                    let mut sim = crate::sim::Simulator::new(num_qubits, num_bits, &mut xof);
                    if ctrl_v { *sim.qubit_mut(ctrl) |= 1; }
                    set_slice_u512_by(&mut sim, &target, U512::from(target_v));
                    set_slice_u512_by(&mut sim, &offset, U512::from(offset_v));
                    sim.apply(&ops);
                    let expected = if ctrl_v { target_v.wrapping_add(offset_v) & maskv } else { target_v & maskv };
                    assert_eq!(get_slice_u512_by(&sim, &target).to::<u64>() & maskv, expected, "target mismatch");
                    assert_eq!(get_slice_u512_by(&sim, &offset).to::<u64>() & maskv, offset_v & maskv, "offset changed");
                    assert_eq!(get_slice_u512_by(&sim, &mask).to::<u64>() & maskv, 0, "mask dirty");
                    assert_eq!(sim.global_phase() & 1, 0, "phase changed");
                }
            }
        }

        let mut b = super::super::B::new();
        let ctrl = b.alloc_qubit();
        let target = b.alloc_qubits(256);
        let offset = b.alloc_qubits(256);
        let mask = b.alloc_qubits(256);
        let clean2 = [b.alloc_qubit(), b.alloc_qubit()];
        for k in 0..256 {
            b.ccx(ctrl, offset[k], mask[k]);
        }
        super::super::venting::iadd_dirty_2clean_qoffset(&mut b, &target, &offset[..254], &clean2, &mask, false);
        for k in 0..256 {
            b.ccx(ctrl, offset[k], mask[k]);
        }
        let ccx = count_ccx(&b.ops);
        let peak = b.peak_qubits as usize;
        let scaled_microstep_with_this_add = ccx + 256 + 255 + 255;
        let div560 = scaled_microstep_with_this_add as f64 * 560.0;
        let scratch_with_compressed_history = 481usize + 26 + 256 + 3;
        eprintln!(
            "masked-borrow controlled qoffset: ccx8={ccx8}, ccx256={ccx}, peak={peak}q, div560≈{div560:.0}, scratch_with_history≈{scratch_with_compressed_history}q"
        );
        println!("METRIC masked_borrow_qoffset_ccx={ccx}");
        println!("METRIC masked_borrow_qoffset_peak={peak}");
        println!("METRIC masked_borrow_qoffset_div560={div560:.0}");
        println!("METRIC masked_borrow_qoffset_scratch_with_history={scratch_with_compressed_history}");
        assert!(ccx < 1_400, "masked-borrow control failed to hit gate target");
        assert!(div560 < 1_200_000.0, "masked-borrow controlled qoffset not replay-shaped");
        assert!(scratch_with_compressed_history > 600, "mask unexpectedly fits the 600 scratch target; revisit BY integration");
    }

    #[test]
    fn streamed_mask_controlled_qoffset_fits_scratch_and_hits_lowqubit_target() {
        // Real version of the mask-streaming idea: keep only a one-qubit
        // ctrl&offset[k] mask at a time, while using an independent dirty bank.
        // Simple q_offset->dst broadcasts are emitted as direct controlled
        // toggles instead of materializing the mask bit.  That small structural
        // change is enough to beat the linear partial-mask tradeoff while
        // keeping the 600-scratch model intact.
        let n = 8usize;
        let maskv = (1u64 << n) - 1;
        let mut b = super::super::B::new();
        let ctrl = b.alloc_qubit();
        let target = b.alloc_qubits(n);
        let offset = b.alloc_qubits(n);
        let dirty = b.alloc_qubits(n - 2);
        let clean2 = [b.alloc_qubit(), b.alloc_qubit()];
        let mask = b.alloc_qubit();
        super::super::venting::ciadd_dirty_3clean_qoffset_stream_mask(&mut b, &target, &dirty, &clean2, mask, &offset, ctrl);
        let num_qubits = b.next_qubit as usize;
        let num_bits = b.next_bit as usize;
        let ops = b.ops;
        for ctrl_v in [false, true] {
            for target_v in [0x00u64, 0x35, 0xf1] {
                for offset_v in [0x00u64, 0x17, 0x80] {
                    let dirty_v = 0x2du64 & ((1u64 << (n - 2)) - 1);
                    let mut hasher = sha3::Shake128::default();
                    hasher.update(b"by-stream-mask-qoffset-small-v1");
                    let mut xof = hasher.finalize_xof();
                    let mut sim = crate::sim::Simulator::new(num_qubits, num_bits, &mut xof);
                    if ctrl_v { *sim.qubit_mut(ctrl) |= 1; }
                    set_slice_u512_by(&mut sim, &target, U512::from(target_v));
                    set_slice_u512_by(&mut sim, &offset, U512::from(offset_v));
                    set_slice_u512_by(&mut sim, &dirty, U512::from(dirty_v));
                    sim.apply(&ops);
                    let expected = if ctrl_v { target_v.wrapping_add(offset_v) & maskv } else { target_v & maskv };
                    assert_eq!(get_slice_u512_by(&sim, &target).to::<u64>() & maskv, expected, "target mismatch ctrl={ctrl_v} target={target_v:x} offset={offset_v:x}");
                    assert_eq!(get_slice_u512_by(&sim, &offset).to::<u64>() & maskv, offset_v & maskv, "offset changed");
                    assert_eq!(get_slice_u512_by(&sim, &dirty).to::<u64>() & ((1u64 << (n - 2)) - 1), dirty_v, "dirty changed");
                    assert_eq!(sim.qubit(mask) & 1, 0, "mask dirty");
                    assert_eq!(sim.global_phase() & 1, 0, "phase changed");
                }
            }
        }
        let mut b = super::super::B::new();
        let ctrl = b.alloc_qubit();
        let target = b.alloc_qubits(256);
        let offset = b.alloc_qubits(256);
        let dirty = b.alloc_qubits(254);
        let clean2 = [b.alloc_qubit(), b.alloc_qubit()];
        let mask = b.alloc_qubit();
        super::super::venting::ciadd_dirty_3clean_qoffset_stream_mask(&mut b, &target, &dirty, &clean2, mask, &offset, ctrl);
        let ccx = count_ccx(&b.ops);
        let peak = b.peak_qubits as usize;
        let scaled_microstep_with_this_add = ccx + 256 + 255 + 255;
        let div560 = scaled_microstep_with_this_add as f64 * 560.0;
        let scaffold_after_div = 642_716usize;
        let branch_decode_margin = 150_000usize;
        let projected_point_add = div560 as usize + scaffold_after_div + branch_decode_margin;
        let scratch_with_history = 481usize + 26 + 3;
        eprintln!(
            "streamed-mask controlled qoffset: ccx={ccx}, peak={peak}q, div560≈{div560:.0}, projected≈{projected_point_add}, scratch_with_history≈{scratch_with_history}q"
        );
        println!("METRIC streamed_mask_qoffset_ccx={ccx}");
        println!("METRIC streamed_mask_qoffset_peak={peak}");
        println!("METRIC streamed_mask_qoffset_div560={div560:.0}");
        println!("METRIC streamed_mask_qoffset_projected_point_add={projected_point_add}");
        println!("METRIC streamed_mask_qoffset_scratch_with_history={scratch_with_history}");
        assert!(ccx < 2_650, "streamed-mask qoffset no longer beats the 600q low-qubit gate target");
        assert!(projected_point_add < 2_700_000, "streamed-mask qoffset misses the Google low-qubit point-add target");
        assert!(scratch_with_history < 600, "streamed-mask qoffset does not fit 600q scratch model");
    }

    #[test]
    fn streamed_mask_qoffset_still_has_no_selector_margin_for_integration() {
        // Early-invalidation guardrail after the scratch-good streamed-mask add:
        // do not hook the BY replay into point-add unless the denominator
        // branch-pattern source is also budgeted.  The 2.645M low-qubit model
        // already spends a 150k selector/decode allowance.  The measured
        // reversible pattern+delta decoder is ~62k, leaving <90k for producing
        // the branch patterns from the quantum denominator.  Even the cheap
        // 16-step lowword pattern oracle costs 5952/window ≈208k for one
        // 560-step denominator before exact full-state plumbing; the known
        // tapered exact generator is far worse.  Therefore the qoffset replay
        // primitive is not, by itself, an integration plan.
        let streamed_projected_with_allowance = 2_645_196usize;
        let allowance = 150_000usize;
        let decoder = 62_160usize;
        let remaining_selector_margin = allowance - decoder;
        let lowword_pattern_oracle_per_window = 5_952usize;
        let windows = 35usize;
        let lowword_one_denominator = lowword_pattern_oracle_per_window * windows;
        let projected_with_lowword_one_denominator =
            streamed_projected_with_allowance - allowance + decoder + lowword_one_denominator;
        let tapered_exact_one_denominator = 2_008_160usize;
        let projected_with_tapered_exact =
            streamed_projected_with_allowance - allowance + decoder + tapered_exact_one_denominator;
        eprintln!(
            "streamed qoffset selector guardrail: remaining_margin={remaining_selector_margin}, lowword_one_den={lowword_one_denominator}, projected_lowword≈{projected_with_lowword_one_denominator}, projected_tapered≈{projected_with_tapered_exact}"
        );
        println!("METRIC streamed_qoffset_selector_margin={remaining_selector_margin}");
        println!("METRIC streamed_qoffset_lowword_selector_ccx={lowword_one_denominator}");
        println!("METRIC streamed_qoffset_projected_with_lowword_selector={projected_with_lowword_one_denominator}");
        println!("METRIC streamed_qoffset_projected_with_tapered_selector={projected_with_tapered_exact}");
        assert!(remaining_selector_margin < lowword_one_denominator, "cheap lowword selector now fits; revisit BY integration");
        assert!(projected_with_lowword_one_denominator > 2_700_000, "lowword selector no longer invalidates integration");
        assert!(projected_with_tapered_exact > 4_000_000, "tapered exact selector unexpectedly SOTA-shaped");
    }

    #[test]
    fn partial_prefix_mask_qoffset_closes_lowqubit_by_near_miss() {
        // A hybrid between the scratch-good streamed-mask qoffset adder and the
        // gate-good full masked-offset adder: keep only a small prefix of
        // `ctrl & offset[k]` masks live, and stream the rest.  Each stored mask
        // bit is reused several times inside the vented carry-xor cleanup, so a
        // few dozen masks can save enough Toffolis while still staying under
        // the 600-scratch BY history budget.
        let n = 8usize;
        let maskv = (1u64 << n) - 1;
        let mut b8 = super::super::B::new();
        let ctrl8 = b8.alloc_qubit();
        let target8 = b8.alloc_qubits(n);
        let offset8 = b8.alloc_qubits(n);
        let dirty8 = b8.alloc_qubits(n - 2);
        let clean8 = [b8.alloc_qubit(), b8.alloc_qubit()];
        let stream8 = b8.alloc_qubit();
        let prefix8 = b8.alloc_qubits(3);
        super::super::venting::ciadd_dirty_3clean_qoffset_partial_mask(
            &mut b8, &target8, &dirty8, &clean8, stream8, &prefix8, &offset8, ctrl8,
        );
        let num_qubits = b8.next_qubit as usize;
        let num_bits = b8.next_bit as usize;
        let ops = b8.ops;
        for ctrl_v in [false, true] {
            for target_v in [0x00u64, 0x35, 0xf1] {
                for offset_v in [0x00u64, 0x17, 0x80] {
                    let dirty_v = 0x2du64 & ((1u64 << (n - 2)) - 1);
                    let mut hasher = sha3::Shake128::default();
                    hasher.update(b"by-partial-prefix-qoffset-small-v1");
                    let mut xof = hasher.finalize_xof();
                    let mut sim = crate::sim::Simulator::new(num_qubits, num_bits, &mut xof);
                    if ctrl_v { *sim.qubit_mut(ctrl8) |= 1; }
                    set_slice_u512_by(&mut sim, &target8, U512::from(target_v));
                    set_slice_u512_by(&mut sim, &offset8, U512::from(offset_v));
                    set_slice_u512_by(&mut sim, &dirty8, U512::from(dirty_v));
                    sim.apply(&ops);
                    let expected = if ctrl_v { target_v.wrapping_add(offset_v) & maskv } else { target_v & maskv };
                    assert_eq!(get_slice_u512_by(&sim, &target8).to::<u64>() & maskv, expected, "target mismatch");
                    assert_eq!(get_slice_u512_by(&sim, &offset8).to::<u64>() & maskv, offset_v & maskv, "offset changed");
                    assert_eq!(get_slice_u512_by(&sim, &dirty8).to::<u64>() & ((1u64 << (n - 2)) - 1), dirty_v, "dirty changed");
                    assert_eq!(get_slice_u512_by(&sim, &prefix8), U512::ZERO, "prefix masks not cleared");
                    assert_eq!(sim.qubit(stream8) & 1, 0, "stream mask not cleared");
                    assert_eq!(sim.global_phase() & 1, 0, "phase changed");
                }
            }
        }

        let mut costs = Vec::new();
        for &prefix_len in &[0usize, 16, 32, 48, 64, 80] {
            let mut b = super::super::B::new();
            let ctrl = b.alloc_qubit();
            let target = b.alloc_qubits(256);
            let offset = b.alloc_qubits(256);
            let dirty = b.alloc_qubits(254);
            let clean2 = [b.alloc_qubit(), b.alloc_qubit()];
            let stream = b.alloc_qubit();
            let prefix = b.alloc_qubits(prefix_len);
            super::super::venting::ciadd_dirty_3clean_qoffset_partial_mask(
                &mut b, &target, &dirty, &clean2, stream, &prefix, &offset, ctrl,
            );
            costs.push((prefix_len, count_ccx(&b.ops), b.peak_qubits as usize));
        }
        let cost0 = costs.iter().find(|&&(m, _, _)| m == 0).unwrap().1;
        let cost32 = costs.iter().find(|&&(m, _, _)| m == 32).unwrap().1;
        let cost48 = costs.iter().find(|&&(m, _, _)| m == 48).unwrap().1;
        let cost64 = costs.iter().find(|&&(m, _, _)| m == 64).unwrap().1;
        let cost80 = costs.iter().find(|&&(m, _, _)| m == 80).unwrap().1;
        let peak32 = costs.iter().find(|&&(m, _, _)| m == 32).unwrap().2;
        let scratch_base = 510usize;
        let scratch32 = scratch_base + 32;
        let scratch48 = scratch_base + 48;
        let scratch64 = scratch_base + 64;
        let scratch80 = scratch_base + 80;
        let harness_cutoff_steps = 564usize;
        let harness_windows = harness_cutoff_steps.div_ceil(16);
        let scaffold_after_div = 642_716usize;
        let lowword_selector = 5_952usize * harness_windows;
        let decoder = 1_776usize * harness_windows;
        let projected_gap = |cost: usize| -> isize {
            let scaled_step = cost + 256 + 255 + 255;
            (scaffold_after_div + lowword_selector + decoder + scaled_step * harness_cutoff_steps) as isize - 2_700_000
        };
        let gap32 = projected_gap(cost32);
        let gap48 = projected_gap(cost48);
        let gap64 = projected_gap(cost64);
        let gap80 = projected_gap(cost80);
        println!("METRIC by_partial_prefix_qoffset_cost0_ccx={cost0}");
        println!("METRIC by_partial_prefix_qoffset_cost32_ccx={cost32}");
        println!("METRIC by_partial_prefix_qoffset_cost48_ccx={cost48}");
        println!("METRIC by_partial_prefix_qoffset_cost64_ccx={cost64}");
        println!("METRIC by_partial_prefix_qoffset_cost80_ccx={cost80}");
        println!("METRIC by_partial_prefix_qoffset_peak32={peak32}");
        println!("METRIC by_partial_prefix_qoffset_windows={harness_windows}");
        println!("METRIC by_partial_prefix_qoffset_selector_ccx={lowword_selector}");
        println!("METRIC by_partial_prefix_qoffset_decoder_ccx={decoder}");
        println!("METRIC by_partial_prefix_qoffset_scratch32={scratch32}");
        println!("METRIC by_partial_prefix_qoffset_scratch48={scratch48}");
        println!("METRIC by_partial_prefix_qoffset_scratch64={scratch64}");
        println!("METRIC by_partial_prefix_qoffset_scratch80={scratch80}");
        println!("METRIC by_partial_prefix_qoffset_projected32_gap_ccx={gap32}");
        println!("METRIC by_partial_prefix_qoffset_projected48_gap_ccx={gap48}");
        println!("METRIC by_partial_prefix_qoffset_projected64_gap_ccx={gap64}");
        println!("METRIC by_partial_prefix_qoffset_projected80_gap_ccx={gap80}");
        eprintln!(
            "BY partial-prefix qoffset costs: {costs:?}, windows={harness_windows}, selector={lowword_selector}, decoder={decoder}, scratch32={scratch32}, gap32={gap32}, scratch48={scratch48}, gap48={gap48}, scratch64={scratch64}, gap64={gap64}, scratch80={scratch80}, gap80={gap80}"
        );
        assert!(scratch32 < 600, "32 prefix masks no longer fit scratch budget");
        assert!(gap32 < 0, "32 prefix masks do not close the harness-scale BY lowword near-miss");
    }

    #[test]
    fn partial_prefix_qoffset_validates_across_widths_and_finds_scratch_limit() {
        // Broaden the previous single n=8/prefix=3 smoke test before treating
        // the partial-prefix qoffset adder as a real low-scratch BY substrate.
        // This is still local primitive validation, not a point-add hook-up.
        for &n in &[8usize, 10, 12, 16] {
            let prefix_cases = [0usize, 1, n / 4, n / 2, n];
            let maskv = if n == 64 { u64::MAX } else { (1u64 << n) - 1 };
            for &prefix_len in &prefix_cases {
                let mut b = super::super::B::new();
                let ctrl = b.alloc_qubit();
                let target = b.alloc_qubits(n);
                let offset = b.alloc_qubits(n);
                let dirty = b.alloc_qubits(n - 2);
                let clean2 = [b.alloc_qubit(), b.alloc_qubit()];
                let stream = b.alloc_qubit();
                let prefix = b.alloc_qubits(prefix_len);
                super::super::venting::ciadd_dirty_3clean_qoffset_partial_mask(
                    &mut b, &target, &dirty, &clean2, stream, &prefix, &offset, ctrl,
                );
                let num_qubits = b.next_qubit as usize;
                let num_bits = b.next_bit as usize;
                let ops = b.ops;
                for ctrl_v in [false, true] {
                    for target_v in [0u64, 0x35, maskv.wrapping_sub(3) & maskv] {
                        for offset_v in [0u64, 0x17, (1u64 << (n - 1)) & maskv] {
                            let dirty_v = 0x5a5au64 & ((1u64 << (n - 2)) - 1);
                            let mut hasher = sha3::Shake128::default();
                            hasher.update(b"by-partial-prefix-qoffset-wide-v1");
                            hasher.update(&(n as u64).to_le_bytes());
                            hasher.update(&(prefix_len as u64).to_le_bytes());
                            let mut xof = hasher.finalize_xof();
                            let mut sim = crate::sim::Simulator::new(num_qubits, num_bits, &mut xof);
                            if ctrl_v { *sim.qubit_mut(ctrl) |= 1; }
                            set_slice_u512_by(&mut sim, &target, U512::from(target_v));
                            set_slice_u512_by(&mut sim, &offset, U512::from(offset_v));
                            set_slice_u512_by(&mut sim, &dirty, U512::from(dirty_v));
                            sim.apply(&ops);
                            let expected = if ctrl_v { target_v.wrapping_add(offset_v) & maskv } else { target_v & maskv };
                            assert_eq!(get_slice_u512_by(&sim, &target).to::<u64>() & maskv, expected, "target mismatch n={n} prefix={prefix_len}");
                            assert_eq!(get_slice_u512_by(&sim, &offset).to::<u64>() & maskv, offset_v & maskv, "offset changed n={n} prefix={prefix_len}");
                            assert_eq!(get_slice_u512_by(&sim, &dirty).to::<u64>() & ((1u64 << (n - 2)) - 1), dirty_v, "dirty changed n={n} prefix={prefix_len}");
                            assert_eq!(get_slice_u512_by(&sim, &prefix), U512::ZERO, "prefix masks dirty n={n} prefix={prefix_len}");
                            assert_eq!(sim.qubit(stream) & 1, 0, "stream mask dirty n={n} prefix={prefix_len}");
                            assert_eq!(sim.global_phase() & 1, 0, "phase changed n={n} prefix={prefix_len}");
                        }
                    }
                }
            }
        }

        let harness_cutoff_steps = 564usize;
        let harness_windows = harness_cutoff_steps.div_ceil(16);
        let scaffold_after_div = 642_716usize;
        let lowword_selector = 5_952usize * harness_windows;
        let decoder = 1_776usize * harness_windows;
        let scratch_base = 510usize;
        let mut best_prefix = 0usize;
        let mut best_cost = usize::MAX;
        let mut best_gap = isize::MAX;
        for prefix_len in (0usize..=90).step_by(10) {
            let mut b = super::super::B::new();
            let ctrl = b.alloc_qubit();
            let target = b.alloc_qubits(256);
            let offset = b.alloc_qubits(256);
            let dirty = b.alloc_qubits(254);
            let clean2 = [b.alloc_qubit(), b.alloc_qubit()];
            let stream = b.alloc_qubit();
            let prefix = b.alloc_qubits(prefix_len);
            super::super::venting::ciadd_dirty_3clean_qoffset_partial_mask(
                &mut b, &target, &dirty, &clean2, stream, &prefix, &offset, ctrl,
            );
            let cost = count_ccx(&b.ops);
            let step = cost + 256 + 255 + 255;
            let gap = (scaffold_after_div + lowword_selector + decoder + step * harness_cutoff_steps) as isize - 2_700_000;
            if scratch_base + prefix_len <= 600 && gap < best_gap {
                best_prefix = prefix_len;
                best_cost = cost;
                best_gap = gap;
            }
        }
        let scratch_best = scratch_base + best_prefix;
        println!("METRIC by_partial_prefix_qoffset_best_prefix_under600={best_prefix}");
        println!("METRIC by_partial_prefix_qoffset_best_cost_under600_ccx={best_cost}");
        println!("METRIC by_partial_prefix_qoffset_best_scratch_under600={scratch_best}");
        println!("METRIC by_partial_prefix_qoffset_best_gap_under600_ccx={best_gap}");
        eprintln!(
            "BY partial-prefix qoffset wide validation passed; best_under600 prefix={best_prefix}, cost={best_cost}, scratch={scratch_best}, gap={best_gap}"
        );
        assert!(best_gap < -100_000, "best under-600 prefix no longer has robust margin");
    }

    #[test]
    fn partial_prefix_qoffset_scratch_schedule_requires_temp_reuse_but_fits() {
        // The prefix90 budget only fits strict scratch if the clean prefix-mask
        // row is not allocated simultaneously with the lowword selector's local
        // simulator. This executable ledger makes that reuse assumption
        // explicit and keeps the route from silently double-counting scratch.
        const W: usize = 16;
        const DBITS: usize = 10;
        let mut b = super::super::B::new();
        let f = b.alloc_qubits(W);
        let g = b.alloc_qubits(W);
        let delta = b.alloc_qubits(DBITS);
        let pattern_tmp = b.alloc_qubits(W);
        let a_tmp = b.alloc_qubits(W);
        let pattern_hist = b.alloc_qubits(W);
        for i in 0..W {
            emit_2adic_by_branch_step_for_test(&mut b, &f, &g, &delta, pattern_tmp[i], a_tmp[i]);
        }
        for i in 0..W {
            b.cx(pattern_tmp[i], pattern_hist[i]);
        }
        for i in (0..W).rev() {
            emit_2adic_by_branch_step_reverse_for_test(&mut b, &f, &g, &delta, pattern_tmp[i], a_tmp[i]);
        }
        let lowword_oracle_peak = b.peak_qubits as usize;
        let lowword_oracle_ccx = count_ccx(&b.ops);

        let compressed_history = 481usize;
        let decoder_control = 26usize;
        let small_clean = 3usize;
        let prefix = 90usize;
        let pattern_bits_in_oracle_peak = W;
        let selector_extra = lowword_oracle_peak - pattern_bits_in_oracle_peak;
        let selector_peak_with_reuse = compressed_history + selector_extra;
        let selector_peak_without_reuse = selector_peak_with_reuse + prefix;
        let replay_peak = compressed_history + decoder_control + small_clean + prefix;
        let scheduled_peak = selector_peak_with_reuse.max(replay_peak);
        println!("METRIC by_partial_prefix_schedule_lowword_oracle_peak={lowword_oracle_peak}");
        println!("METRIC by_partial_prefix_schedule_lowword_oracle_ccx={lowword_oracle_ccx}");
        println!("METRIC by_partial_prefix_schedule_selector_peak_reuse={selector_peak_with_reuse}");
        println!("METRIC by_partial_prefix_schedule_selector_peak_no_reuse={selector_peak_without_reuse}");
        println!("METRIC by_partial_prefix_schedule_replay_peak={replay_peak}");
        println!("METRIC by_partial_prefix_schedule_peak={scheduled_peak}");
        eprintln!(
            "BY partial-prefix scratch schedule: lowword_peak={lowword_oracle_peak}, selector_reuse={selector_peak_with_reuse}, selector_no_reuse={selector_peak_without_reuse}, replay={replay_peak}, scheduled={scheduled_peak}"
        );
        assert!(selector_peak_without_reuse > 600, "temp reuse is no longer required; simplify the schedule");
        assert!(scheduled_peak <= 600, "partial-prefix qoffset scratch schedule no longer fits strict 600");
    }

    #[test]
    fn partial_prefix_qoffset_two_denominator_ledger_blocks_naive_promotion() {
        // Adversarial-accountant response to the prior BY blow-up: the one-DIV
        // budget is not enough evidence.  If the architecture needs separate
        // pair1 tagged-DIV and pair2 product-clean replay, the partial-prefix
        // qoffset win is overwhelmed.  Keep this guardrail before any hook-up.
        let prefix90_qoffset = 2_094usize;
        let step = prefix90_qoffset + 256 + 255 + 255;
        let steps = 564usize;
        let windows = steps.div_ceil(16);
        let selector_per_den = 5_952usize * windows;
        let decoder_per_den = 1_776usize * windows;

        let one_div_scaffold = 642_716usize;
        let one_div_total = one_div_scaffold + step * steps + selector_per_den + decoder_per_den;
        let one_div_gap = one_div_total as isize - 2_700_000;

        let current_total = 4_132_750isize;
        let current_two_kaliski = 3_190_000isize;
        let deleted_pair1_muls = 149_889isize + 150_145isize;
        let pair1_scale_loop = 407isize * 255;
        let pair2_scale_loop = 404isize * 255;
        let pair2_product_mul = 150_145isize;
        let two_replay_scaffold = current_total
            - current_two_kaliski
            - deleted_pair1_muls
            - pair1_scale_loop
            - pair2_scale_loop
            - pair2_product_mul;
        let two_den_total = two_replay_scaffold
            + 2 * (step * steps + selector_per_den + decoder_per_den) as isize;
        let two_den_gap = two_den_total - 2_700_000;
        let missing_assumption_swing = two_den_total - one_div_total as isize;
        println!("METRIC by_partial_prefix_one_div_total={one_div_total}");
        println!("METRIC by_partial_prefix_one_div_gap_ccx={one_div_gap}");
        println!("METRIC by_partial_prefix_two_den_scaffold={two_replay_scaffold}");
        println!("METRIC by_partial_prefix_two_den_total={two_den_total}");
        println!("METRIC by_partial_prefix_two_den_gap_ccx={two_den_gap}");
        println!("METRIC by_partial_prefix_missing_second_den_swing_ccx={missing_assumption_swing}");
        eprintln!(
            "BY partial-prefix adversarial two-den ledger: one_div_total={one_div_total} gap={one_div_gap}, two_den_total={two_den_total} gap={two_den_gap}, swing={missing_assumption_swing}"
        );
        assert!(one_div_gap < 0, "one-DIV optimistic model no longer fits; update prior tests");
        assert!(two_den_gap > 500_000, "two-denominator promotion might fit; revisit before demoting");
    }

    #[test]
    fn partial_prefix_one_div_is_too_expensive_for_strategy_e_escape() {
        // After the two-denominator ledger, the obvious escape is to use the
        // one-DIV partial-prefix primitive in a slope-coordinate/Strategy-E
        // point-add map.  That route still needs an in-place variable multiply
        // to convert the carried slope to affine y.  Even granting a future
        // schoolbook-like product-clean multiply, the partial-prefix DIV body is
        // too expensive for the Strategy-E budget.
        let partial_prefix_one_div_total = 2_533_964usize;
        let one_div_scaffold_from_by_model = 642_716usize;
        let partial_prefix_div_body = partial_prefix_one_div_total - one_div_scaffold_from_by_model;
        let strategy_e_non_div_scaffold = 942_750usize;
        let known_product_clean = 1_145_760usize;
        let hypothetical_schoolbook_product_clean = 180_000usize;
        let current_product_total = strategy_e_non_div_scaffold + partial_prefix_div_body + known_product_clean;
        let hypothetical_product_total = strategy_e_non_div_scaffold + partial_prefix_div_body + hypothetical_schoolbook_product_clean;
        let current_gap = current_product_total as isize - 2_700_000;
        let hypothetical_gap = hypothetical_product_total as isize - 2_700_000;
        println!("METRIC by_partial_prefix_strategy_e_div_body_ccx={partial_prefix_div_body}");
        println!("METRIC by_partial_prefix_strategy_e_current_product_gap_ccx={current_gap}");
        println!("METRIC by_partial_prefix_strategy_e_hyp_product_gap_ccx={hypothetical_gap}");
        eprintln!(
            "partial-prefix one-DIV in Strategy E: div_body={partial_prefix_div_body}, current_total={current_product_total} gap={current_gap}, hypothetical_total={hypothetical_product_total} gap={hypothetical_gap}"
        );
        assert!(hypothetical_gap > 0, "partial-prefix DIV plus schoolbook product would make Strategy E viable; revisit");
    }

    #[test]
    fn strategy_e_three_million_still_needs_unavailable_div_selector_budget() {
        // The user's relaxed question is whether anything can reach ~3M while
        // staying in the low-qubit regime.  Strategy E is the cleanest algebraic
        // way to delete BY's second denominator, but with the current known
        // product-clean multiply it leaves only a tiny budget for the single
        // DIV's denominator controls.  This guardrail prevents counting the
        // fixed-control replay body as a solved low-scratch DIV.
        let target_3m = 3_000_000usize;
        let strategy_e_non_div_scaffold = 942_750usize;
        let known_product_clean = 1_145_760usize;
        let best_fixed_control_replay = 873_600usize;
        let measured_decoder = 1_776usize * 36usize;
        let measured_lowword_selector = 5_952usize * 36usize;
        let max_div_body_for_3m = target_3m - strategy_e_non_div_scaffold - known_product_clean;
        let selector_budget_after_replay = max_div_body_for_3m as isize - best_fixed_control_replay as isize;
        let decoder_only_gap = (strategy_e_non_div_scaffold
            + known_product_clean
            + best_fixed_control_replay
            + measured_decoder) as isize
            - target_3m as isize;
        let lowword_gap = (strategy_e_non_div_scaffold
            + known_product_clean
            + best_fixed_control_replay
            + measured_decoder
            + measured_lowword_selector) as isize
            - target_3m as isize;
        let partial_prefix_div_body = 1_891_248usize;
        let partial_prefix_gap = (strategy_e_non_div_scaffold
            + known_product_clean
            + partial_prefix_div_body) as isize
            - target_3m as isize;
        println!("METRIC strategy_e_3m_max_div_body_ccx={max_div_body_for_3m}");
        println!("METRIC strategy_e_3m_selector_budget_after_fixed_replay_ccx={selector_budget_after_replay}");
        println!("METRIC strategy_e_3m_decoder_only_gap_ccx={decoder_only_gap}");
        println!("METRIC strategy_e_3m_lowword_selector_gap_ccx={lowword_gap}");
        println!("METRIC strategy_e_3m_partial_prefix_gap_ccx={partial_prefix_gap}");
        eprintln!(
            "Strategy E 3M risk ledger: max_div_body={max_div_body_for_3m}, selector_budget_after_replay={selector_budget_after_replay}, decoder_only_gap={decoder_only_gap}, lowword_gap={lowword_gap}, partial_prefix_gap={partial_prefix_gap}"
        );
        assert!(selector_budget_after_replay < measured_decoder as isize, "fixed replay plus measured decoder now fits under 3M; revisit Strategy E");
        assert!(lowword_gap > 0, "measured lowword selector would make Strategy E <3M; wire a toy schedule next");
        assert!(partial_prefix_gap > 0, "partial-prefix DIV is enough for Strategy E under 3M; revisit");
    }

    #[test]
    fn strategy_e_current_product_clean_is_still_a_second_denominator() {
        // Stronger risk check: the "known product-clean multiply" used in old
        // Strategy-E budgets is not a generic near-schoolbook in-place multiply;
        // it is itself a denominator-controlled BY/product-clean replay on
        // c=Rx-Qx.  Therefore Strategy E only deletes the second denominator if
        // a new non-DIV product-clean multiply exists.  With current primitives
        // the second selector/parser must be charged too.
        let target_3m = 3_000_000isize;
        let strategy_e_non_div_scaffold = 942_750isize;
        let one_div_fixed_replay = 873_600isize;
        let current_product_clean_div_replay = 1_145_760isize;
        let optimistic_centered_product_replay = 873_600isize;
        let selector_decoder_per_den = (5_952isize + 1_776isize) * 36isize;
        let forgotten_second_selector_total = strategy_e_non_div_scaffold
            + one_div_fixed_replay
            + optimistic_centered_product_replay
            + selector_decoder_per_den;
        let optimistic_two_den_total = forgotten_second_selector_total + selector_decoder_per_den;
        let current_two_den_total = strategy_e_non_div_scaffold
            + one_div_fixed_replay
            + current_product_clean_div_replay
            + 2 * selector_decoder_per_den;
        let forgotten_second_selector_gap = forgotten_second_selector_total - target_3m;
        let optimistic_two_den_gap = optimistic_two_den_total - target_3m;
        let current_two_den_gap = current_two_den_total - target_3m;
        let missing_second_selector_swing = selector_decoder_per_den;
        let partial_prefix_div_body = 1_891_248isize;
        let product_clean_budget_after_partial_prefix = target_3m - strategy_e_non_div_scaffold - partial_prefix_div_body;
        println!("METRIC strategy_e_second_den_selector_decoder_ccx={selector_decoder_per_den}");
        println!("METRIC strategy_e_forgotten_second_selector_gap_ccx={forgotten_second_selector_gap}");
        println!("METRIC strategy_e_optimistic_two_den_gap_to_3m_ccx={optimistic_two_den_gap}");
        println!("METRIC strategy_e_current_two_den_gap_to_3m_ccx={current_two_den_gap}");
        println!("METRIC strategy_e_missing_second_selector_swing_ccx={missing_second_selector_swing}");
        println!("METRIC strategy_e_product_clean_budget_after_partial_prefix_ccx={product_clean_budget_after_partial_prefix}");
        eprintln!(
            "Strategy E current product-clean risk: forgotten_second_selector_total={forgotten_second_selector_total} gap={forgotten_second_selector_gap}, optimistic_two_den_total={optimistic_two_den_total} gap={optimistic_two_den_gap}, current_two_den_total={current_two_den_total} gap={current_two_den_gap}"
        );
        assert!(forgotten_second_selector_gap < 0, "missing-second-selector trap no longer looks tempting; update note");
        assert!(optimistic_two_den_gap > 0, "even charging two selectors, optimistic Strategy E fits 3M; revisit");
        assert!(product_clean_budget_after_partial_prefix < 180_000, "partial-prefix Strategy E can afford a schoolbook-like product-clean multiply under 3M");
    }

    fn secp_curve_for_strategy_e_share_test() -> crate::weierstrass_elliptic_curve::WeierstrassEllipticCurve {
        crate::weierstrass_elliptic_curve::WeierstrassEllipticCurve {
            modulus: SECP256K1_P,
            a: U256::from(0),
            b: U256::from(7),
            gx: U256::from_str_radix(
                "79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798",
                16,
            ).unwrap(),
            gy: U256::from_str_radix(
                "483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8",
                16,
            ).unwrap(),
            order: U256::from_str_radix(
                "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141",
                16,
            ).unwrap(),
        }
    }

    fn rand_scalar_for_strategy_e_share_test(rng: &mut u64, order: U256) -> U256 {
        let mut limbs = [0u64; 4];
        for limb in &mut limbs {
            *rng = rng.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            *limb = *rng;
        }
        let x = U256::from_limbs(limbs) % order;
        if x.is_zero() { U256::from(1u64) } else { x }
    }

    fn by_branch_cases_for_strategy_e_share_test(g0: U256, p: U256, iters: usize) -> Vec<u8> {
        let mut delta: i64 = 1;
        let mut f = SInt::from_u(p);
        let mut g = SInt::from_u(g0);
        let mut out = Vec::with_capacity(iters);
        for _ in 0..iters {
            let g_odd = g.bit0();
            if delta > 0 && g_odd {
                out.push(2); // A
                let nf = g;
                let ng = SInt::sub(g, f).shr1_even();
                delta = 1 - delta;
                f = nf;
                g = ng;
            } else if g_odd {
                out.push(1); // B
                let ng = SInt::add(g, f).shr1_even();
                delta = 1 + delta;
                g = ng;
            } else {
                out.push(0); // C
                let ng = g.shr1_even();
                delta = 1 + delta;
                g = ng;
            }
        }
        out
    }

    fn mutual_information_millibits_for_strategy_e_share_test<const N: usize>(counts: [[usize; N]; N]) -> f64 {
        let total: usize = counts.iter().flatten().sum();
        if total == 0 { return 0.0; }
        let mut rows = [0usize; N];
        let mut cols = [0usize; N];
        for i in 0..N {
            for j in 0..N {
                rows[i] += counts[i][j];
                cols[j] += counts[i][j];
            }
        }
        let total_f = total as f64;
        let mut mi = 0.0f64;
        for i in 0..N {
            for j in 0..N {
                let c = counts[i][j];
                if c == 0 { continue; }
                let pxy = c as f64 / total_f;
                let px = rows[i] as f64 / total_f;
                let py = cols[j] as f64 / total_f;
                mi += pxy * (pxy / (px * py)).log2();
            }
        }
        1000.0 * mi
    }

    #[test]
    fn strategy_e_product_denominator_controls_do_not_share_pair1_branch_stream() {
        // If Strategy E could derive the product-clean denominator controls
        // from the slope-DIV branch stream, the "second denominator" risk above
        // might be avoidable.  A direct secp sample shows no such simple sharing:
        // BY branch cases for dx and for c=Rx-Qx are essentially independent.
        // This does not prove every algebra impossible, but it blocks counting
        // branch-control reuse without a new concrete invariant.
        let curve = secp_curve_for_strategy_e_share_test();
        let p = SECP256K1_P;
        let iters = 576usize;
        let target_samples = 256usize;
        let mut rng = 0x57aa_7e61_5bad_e001u64;
        let mut samples = 0usize;
        let mut odd_counts = [[0usize; 2]; 2];
        let mut case_counts = [[0usize; 3]; 3];
        let mut odd_match = 0usize;
        let mut case_match = 0usize;
        while samples < target_samples {
            let k1 = rand_scalar_for_strategy_e_share_test(&mut rng, curve.order);
            let k2 = rand_scalar_for_strategy_e_share_test(&mut rng, curve.order);
            let (px, py) = curve.mul(curve.gx, curve.gy, k1);
            let (qx, qy) = curve.mul(curve.gx, curve.gy, k2);
            if (px.is_zero() && py.is_zero()) || (qx.is_zero() && qy.is_zero()) || px == qx {
                continue;
            }
            let (rx, ry) = curve.add(px, py, qx, qy);
            if rx.is_zero() && ry.is_zero() { continue; }
            let dx = subm(px, qx, p);
            let prod_den = subm(rx, qx, p);
            if dx.is_zero() || prod_den.is_zero() { continue; }
            let a = by_branch_cases_for_strategy_e_share_test(dx, p, iters);
            let b = by_branch_cases_for_strategy_e_share_test(prod_den, p, iters);
            for i in 0..iters {
                let oa = (a[i] != 0) as usize;
                let ob = (b[i] != 0) as usize;
                odd_counts[oa][ob] += 1;
                case_counts[a[i] as usize][b[i] as usize] += 1;
                if oa == ob { odd_match += 1; }
                if a[i] == b[i] { case_match += 1; }
            }
            samples += 1;
        }
        let total = samples * iters;
        let odd_match_ppm = odd_match * 1_000_000usize / total;
        let case_match_ppm = case_match * 1_000_000usize / total;
        let odd_mi_millibits = mutual_information_millibits_for_strategy_e_share_test::<2>(odd_counts);
        let case_mi_millibits = mutual_information_millibits_for_strategy_e_share_test::<3>(case_counts);
        println!("METRIC strategy_e_branch_share_samples={samples}");
        println!("METRIC strategy_e_branch_odd_match_ppm={odd_match_ppm}");
        println!("METRIC strategy_e_branch_case_match_ppm={case_match_ppm}");
        println!("METRIC strategy_e_branch_odd_mi_millibits={odd_mi_millibits:.3}");
        println!("METRIC strategy_e_branch_case_mi_millibits={case_mi_millibits:.3}");
        eprintln!(
            "Strategy E branch sharing probe: samples={samples}, odd_match_ppm={odd_match_ppm}, case_match_ppm={case_match_ppm}, odd_mi_millibits={odd_mi_millibits:.3}, case_mi_millibits={case_mi_millibits:.3}, odd_counts={odd_counts:?}, case_counts={case_counts:?}"
        );
        assert!(odd_mi_millibits < 10.0, "pair1/product odd controls show exploitable correlation; investigate sharing");
        assert!(case_mi_millibits < 10.0, "pair1/product branch cases show exploitable correlation; investigate sharing");
    }

    #[test]
    fn endomorphism_y_denominators_do_not_share_branch_stream() {
        // The j=0 endomorphism identity rewrites the input slope as
        //   lambda = ((beta*Px-Qx)(beta^2*Px-Qx))/(Py+Qy).
        // Applying the same identity to the output subtraction uses denominator
        // Ry-Qy.  If these y-denominator branch streams were correlated, this
        // could be a real control-sharing escape from the two-denominator BY
        // obstruction.  Sample the exact branch cases directly before granting
        // any selector reuse budget.
        let curve = secp_curve_for_strategy_e_share_test();
        let p = SECP256K1_P;
        let iters = 576usize;
        let target_samples = 256usize;
        let mut rng = 0x0e5d_0d15_c0fe_0001u64;
        let mut samples = 0usize;
        let mut odd_counts = [[0usize; 2]; 2];
        let mut case_counts = [[0usize; 3]; 3];
        let mut odd_match = 0usize;
        let mut case_match = 0usize;
        while samples < target_samples {
            let k1 = rand_scalar_for_strategy_e_share_test(&mut rng, curve.order);
            let k2 = rand_scalar_for_strategy_e_share_test(&mut rng, curve.order);
            let (px, py) = curve.mul(curve.gx, curve.gy, k1);
            let (qx, qy) = curve.mul(curve.gx, curve.gy, k2);
            if (px.is_zero() && py.is_zero()) || (qx.is_zero() && qy.is_zero()) || px == qx {
                continue;
            }
            let (rx, ry) = curve.add(px, py, qx, qy);
            if rx.is_zero() && ry.is_zero() { continue; }
            let in_y_den = addm(py, qy, p);
            let out_y_den = subm(ry, qy, p);
            if in_y_den.is_zero() || out_y_den.is_zero() { continue; }
            let a = by_branch_cases_for_strategy_e_share_test(in_y_den, p, iters);
            let b = by_branch_cases_for_strategy_e_share_test(out_y_den, p, iters);
            for i in 0..iters {
                let oa = (a[i] != 0) as usize;
                let ob = (b[i] != 0) as usize;
                odd_counts[oa][ob] += 1;
                case_counts[a[i] as usize][b[i] as usize] += 1;
                if oa == ob { odd_match += 1; }
                if a[i] == b[i] { case_match += 1; }
            }
            samples += 1;
        }
        let total = samples * iters;
        let odd_match_ppm = odd_match * 1_000_000usize / total;
        let case_match_ppm = case_match * 1_000_000usize / total;
        let odd_mi_millibits = mutual_information_millibits_for_strategy_e_share_test::<2>(odd_counts);
        let case_mi_millibits = mutual_information_millibits_for_strategy_e_share_test::<3>(case_counts);
        println!("METRIC endomorphism_y_branch_share_samples={samples}");
        println!("METRIC endomorphism_y_branch_odd_match_ppm={odd_match_ppm}");
        println!("METRIC endomorphism_y_branch_case_match_ppm={case_match_ppm}");
        println!("METRIC endomorphism_y_branch_odd_mi_millibits={odd_mi_millibits:.3}");
        println!("METRIC endomorphism_y_branch_case_mi_millibits={case_mi_millibits:.3}");
        eprintln!(
            "Endomorphism y-denominator branch sharing probe: samples={samples}, odd_match_ppm={odd_match_ppm}, case_match_ppm={case_match_ppm}, odd_mi_millibits={odd_mi_millibits:.3}, case_mi_millibits={case_mi_millibits:.3}, odd_counts={odd_counts:?}, case_counts={case_counts:?}"
        );
        assert!(odd_mi_millibits < 10.0, "endomorphism y-denominator odd controls show exploitable correlation; investigate sharing");
        assert!(case_mi_millibits < 10.0, "endomorphism y-denominator branch cases show exploitable correlation; investigate sharing");
    }

    #[test]
    fn partial_mask_controlled_qoffset_linear_tradeoff_just_misses_600q_target() {
        // First-order model after the masked-borrow primitive: full mask gives
        // good gates but 766q scratch with compressed history; no mask gives
        // good scratch but 3557-CCX add.  The 600q cap leaves only ~90 clean
        // mask bits after 481 pattern-history bits and 26 decoder bits.  Even a
        // generous linear interpolation between the measured endpoints lands
        // just above the 2.7M low-qubit point-add target once branch/decode
        // margin is included, so mask streaming must beat linear interpolation.
        let no_mask_ccx = 3557usize;
        let full_mask_ccx = 1274usize;
        let full_mask_bits = 256usize;
        let compressed_history = 481usize;
        let decoder = 26usize;
        let small_clean = 3usize;
        let scratch_cap = 600usize;
        let mask_budget = scratch_cap - compressed_history - decoder - small_clean;
        let saved = (no_mask_ccx - full_mask_ccx) * mask_budget / full_mask_bits;
        let interpolated_add_ccx = no_mask_ccx - saved;
        let non_add_step_overhead = 256usize + 255 + 255;
        let replay560 = (interpolated_add_ccx + non_add_step_overhead) * 560;
        let scaffold_after_div = 642_716usize;
        let branch_decode_margin = 150_000usize;
        let projected = replay560 + scaffold_after_div + branch_decode_margin;
        eprintln!(
            "partial mask qoffset model: mask_budget={mask_budget}, add_ccx≈{interpolated_add_ccx}, replay560≈{replay560}, projected≈{projected}"
        );
        println!("METRIC partial_mask_qoffset_mask_budget={mask_budget}");
        println!("METRIC partial_mask_qoffset_interpolated_add_ccx={interpolated_add_ccx}");
        println!("METRIC partial_mask_qoffset_projected_point_add={projected}");
        assert!(mask_budget < 100);
        assert!(projected > 2_700_000, "linear partial-mask model would already hit low-qubit SOTA");
    }

    #[test]
    fn naive_controlled_dirty_qoffset_is_scratch_right_but_toffoli_heavy() {
        let mut b = super::super::B::new();
        let ctrl = b.alloc_qubit();
        let target = b.alloc_qubits(256);
        let offset = b.alloc_qubits(256);
        let dirty = b.alloc_qubits(254);
        let clean2 = [b.alloc_qubit(), b.alloc_qubit()];
        let clean3 = b.alloc_qubit();
        super::super::venting::ciadd_dirty_3clean_qoffset(&mut b, &target, &dirty, &clean2, clean3, &offset, ctrl);
        let ccx = count_ccx(&b.ops);
        let peak = b.peak_qubits as usize;
        let scaled_microstep_with_this_add = ccx + 256 + 255 + 255;
        let div560 = scaled_microstep_with_this_add as f64 * 560.0;
        eprintln!(
            "naive controlled dirty qoffset add: ccx={ccx}, peak={peak}q, dirty={}q, clean=3q, scaled_step≈{scaled_microstep_with_this_add}, div560≈{div560:.0}",
            dirty.len()
        );
        assert_eq!(peak, 1 + target.len() + offset.len() + dirty.len() + 3, "unexpected hidden clean workspace");
        assert!(ccx > 3_000, "naive controlled dirty qoffset unexpectedly cheap; revisit BY cmod_add rewrite");
        assert!(div560 > 2_000_000.0, "naive controlled dirty qoffset would be SOTA-shaped after all");
    }

    #[test]
    fn dirty_quantum_offset_adder_is_plausible_cmod_add_substrate() {
        // Existing venting code already has the right primitive shape for the
        // missing no-clean-temp add: add a quantum offset using only two clean
        // qubits plus n-2 dirty qubits. It is not controlled and not modular by
        // p, so it is not a drop-in replacement for cmod_add_qq yet, but it
        // shows the target scratch/cost scale is realistic rather than magical.
        let mut b = super::super::B::new();
        let target = b.alloc_qubits(256);
        let offset = b.alloc_qubits(256);
        let dirty = b.alloc_qubits(254);
        let clean2 = [b.alloc_qubit(), b.alloc_qubit()];
        super::super::venting::iadd_dirty_2clean_qoffset(&mut b, &target, &dirty, &clean2, &offset, false);
        let ccx = count_ccx(&b.ops);
        let peak = b.peak_qubits as usize;
        eprintln!(
            "dirty quantum-offset add substrate: ccx={ccx}, peak={peak}q, dirty={}q, clean=2q",
            dirty.len()
        );
        assert!(ccx < 1_000, "dirty qoffset add too expensive for cmod_add rewrite substrate");
        assert_eq!(peak, target.len() + offset.len() + dirty.len() + 2, "unexpected hidden clean workspace");
    }

    #[test]
    fn compressed_pattern_history_scratch_model_is_600q_if_add_workspace_is_removed() {
        // Peak accounting for the scaled BY DIV core. The current cmod_add
        // implementation uses clean 256-bit addend/carry workspaces, so history
        // and arithmetic temp add. If controlled modular add is implemented with
        // no clean n-bit temp (or uses the history bank as dirty workspace), the
        // compressed branch-pattern history becomes the dominant scratch item.
        let pair_regs = 512usize;
        let current_microstep_peak = 1287usize;
        let current_local_workspace = current_microstep_peak - pair_regs;
        let compressed_pattern_bits = 481usize; // fixed per-window distinct-pattern IDs from entropy test.
        let decoder_bits = 16usize + 10usize; // A scratch for one window + signed delta.
        let current_scratch = compressed_pattern_bits + decoder_bits + current_local_workspace;
        let no_clean_temp_workspace = 90usize; // target: controls/ext/carry only, no 256-bit f register.
        let target_scratch = compressed_pattern_bits + decoder_bits + no_clean_temp_workspace;
        eprintln!(
            "BY compressed-pattern scratch model: current_scratch≈{current_scratch}, target_no_clean_temp≈{target_scratch}, local_workspace_now={current_local_workspace}"
        );
        assert!(current_scratch > 1_000, "current clean-temp implementation already fits 600 scratch; update model");
        assert!(target_scratch <= 620, "no-clean-temp target no longer near 600 scratch");
    }

    #[test]
    fn clean_two_replay_by_budget_requires_replay_or_phase_breakthrough() {
        // Updated after the clean pattern-decoder scaffold: if both BY replays
        // use the current all-A-history clean schedule (decode, replay, reverse
        // decoder), the point-add budget is slightly above 2.7M before any
        // denominator branch generation. This keeps the moonshot target honest:
        // we need a cheaper replay (fused average / fixed-control window) or a
        // phase-safe way to clear A locally, not just denominator plumbing.
        let current_total = 4_111_918.0;
        let non_inv_scaffold = 942_750.0;
        let deleted_pair1_muls = 149_889.0 + 150_145.0;
        let two_replay_scaffold = non_inv_scaffold
            - deleted_pair1_muls
            - 407.0 * 255.0
            - 404.0 * 255.0
            - 150_145.0;
        let fast_raw_replay = 1_145_760.0;
        let decoded_forward_replay = 1_207_920.0;
        let clean_decoded_replay = 1_270_080.0;
        let clean_two_replay_total = two_replay_scaffold + 2.0 * clean_decoded_replay;
        let forward_only_two_replay_total = two_replay_scaffold + 2.0 * decoded_forward_replay;
        let fixed_control_replay = 800_900.0; // measured fixed branch-numerator lower target.
        let fixed_control_two_replay_with_300k_branch = two_replay_scaffold + 2.0 * fixed_control_replay + 300_000.0;
        eprintln!(
            "clean BY two-replay budget: scaffold≈{two_replay_scaffold:.0}, forward_only≈{forward_only_two_replay_total:.0}, clean≈{clean_two_replay_total:.0}, fixed_control_plus300k_branch≈{fixed_control_two_replay_with_300k_branch:.0}"
        );
        assert!(forward_only_two_replay_total > 2_690_000.0, "current forward-only decoded replay has large hidden margin; update model");
        assert!(clean_two_replay_total > 2_700_000.0, "current clean decoded replay already beats 2.7M; denominator generation is the only blocker");
        assert!(fixed_control_two_replay_with_300k_branch < 2_200_000.0, "fixed-control replay target no longer has low-gate margin");
    }

    #[test]
    fn centered_signed_replay_budget_hits_google_if_parity_and_denominator_are_folded() {
        let non_inv_scaffold = 942_750.0;
        let deleted_pair1_muls = 149_889.0 + 150_145.0;
        let two_replay_scaffold = non_inv_scaffold
            - deleted_pair1_muls
            - 407.0 * 255.0
            - 404.0 * 255.0
            - 150_145.0;
        let centered_replay = 873_600.0;
        let tapered_den_compute_per_denominator = 303_828.0;
        let naive_parity_cleanup_per_replay = 725_760.0;
        let folded_parity_selected_denominator = two_replay_scaffold
            + 2.0 * centered_replay
            + 2.0 * tapered_den_compute_per_denominator;
        let naive_parity_total = folded_parity_selected_denominator + 2.0 * naive_parity_cleanup_per_replay;
        eprintln!(
            "centered signed BY budget: scaffold≈{two_replay_scaffold:.0}, two_replay≈{:.0}, +selected_den≈{folded_parity_selected_denominator:.0}, naive_parity≈{naive_parity_total:.0}",
            2.0 * centered_replay
        );
        assert!(folded_parity_selected_denominator < 2_700_000.0, "centered replay no longer has low-qubit SOTA margin with selected denominator");
        assert!(naive_parity_total > 3_500_000.0, "naive parity cleanup unexpectedly acceptable");
    }

    #[test]
    fn low_scratch_scaled_by_budget_still_beats_27m_after_pair1_mul_deletion() {
        // Important refinement: a true tagged DIV does not merely replace the
        // two Kaliski bodies. It also deletes pair1's two schoolbook
        // multiplications because the DIV itself maps (0,dy+dx)->(lambda+1,0).
        // That ~300k saving gives enough margin for the low-scratch vented
        // modular-add variant, whose measured step cost is higher but peak is
        // much closer to the 600-scratch target.
        let current_total = 4_132_750.0;
        let current_two_kaliski = 3_190_000.0;
        let non_inv_scaffold = current_total - current_two_kaliski;
        let deleted_pair1_muls = 149_889.0 + 150_145.0;
        let scaffold_after_div = non_inv_scaffold - deleted_pair1_muls;
        let fast_by_div = 2_046.0 * 560.0;
        let low_scratch_vented_by_div = 3_318.0 * 560.0; // measured with KAL_VENT_MODADD=1.
        let branch_decode_margin = 150_000.0;
        let fast_projected = scaffold_after_div + fast_by_div + branch_decode_margin;
        let low_scratch_projected = scaffold_after_div + low_scratch_vented_by_div + branch_decode_margin;
        let two_replay_scaffold = current_total
            - current_two_kaliski
            - deleted_pair1_muls
            - 407.0 * 255.0 // pair1 scale loop
            - 404.0 * 255.0 // pair2 double/scale loop
            - 150_145.0; // pair2 schoolbook product add
        let two_replay_decode_margin = 60_000.0; // two pattern decoders are ~46k by pattern_decoder_budget.
        let two_fast_replays = two_replay_scaffold + 2.0 * fast_by_div + two_replay_decode_margin;
        eprintln!(
            "BY DIV with pair1 mul deletion: scaffold_after_div≈{scaffold_after_div:.0}, fast_projected≈{fast_projected:.0}, low_scratch_vented_projected≈{low_scratch_projected:.0}, two_fast_replays≈{two_fast_replays:.0}"
        );
        assert!(fast_projected < 2_100_000.0, "fast BY DIV no longer reaches low-gate target band");
        assert!(low_scratch_projected < 2_700_000.0, "low-scratch vented BY DIV no longer beats 2.7M");
        assert!(two_fast_replays < 2_700_000.0, "two fast BY replays no longer beat 2.7M");
    }

    #[test]
    fn naive_quantum_branch_generator_would_erase_scaled_by_savings() {
        // The benchmark-path blocker is not the modular replay any more; it is
        // generating the 560 data-dependent branch bits reversibly from the
        // quantum denominator.  A direct 2-adic f/g generator keeps enough
        // precision to emit all branch bits, but per step it still performs a
        // full-width controlled swap, controlled negation, and controlled add.
        // This executable budget is the guardrail: wiring this naive generator
        // would be a real circuit, but not a SOTA circuit.
        let steps = 560.0;
        let width = 560.0; // 2-adic precision needed to stream 560 branch bits.
        let cswap_fg = width; // one Fredkin per bit.
        let cneg_g = width; // controlled bit flips plus one controlled increment (Toffoli-scale).
        let cadd_fg = 3.0 * width; // load ctrl&f, ripple add, unload (optimistic fast-adder count).
        let delta_logic = 80.0; // signed-delta positivity/update, deliberately small.
        let branch_gen_step = cswap_fg + cneg_g + cadd_fg + delta_logic;
        let one_generator_compute_uncompute = 2.0 * steps * branch_gen_step;
        let two_generators = 2.0 * one_generator_compute_uncompute;
        let two_replays = 2.0 * 2_046.0 * steps;
        let current_total = 4_132_750.0;
        let current_two_kaliski = 3_190_000.0;
        let deleted_pair1_muls = 149_889.0 + 150_145.0;
        let two_replay_scaffold = current_total
            - current_two_kaliski
            - deleted_pair1_muls
            - 407.0 * 255.0
            - 404.0 * 255.0
            - 150_145.0;
        let projected_with_naive_generators = two_replay_scaffold + two_replays + two_generators;
        eprintln!(
            "naive BY 2-adic branch generator: step≈{branch_gen_step:.0}, one_compute_uncompute≈{one_generator_compute_uncompute:.0}, two_generators≈{two_generators:.0}, projected≈{projected_with_naive_generators:.0}"
        );
        assert!(one_generator_compute_uncompute > 2_000_000.0, "naive generator unexpectedly cheap; revisit integration");
        assert!(projected_with_naive_generators > current_total, "naive branch generation would already be a saving; model is stale");
    }

    #[test]
    fn scaled_by_div_point_add_budget_has_sota_margin_if_history_workspace_solved() {
        // The structural point of the scaled controlled microstep is that it
        // replaces both Kaliski invocations by one in-place tagged DIV. This is
        // a budget model, not an implementation: it assumes the remaining
        // branch-history/workspace problem is solved without changing the
        // measured 2046-CCX microstep arithmetic.
        let current_total = 4_132_750.0;
        let current_two_kaliski = 3_190_000.0; // measured ~1.60M + ~1.59M from point-add decomposition.
        let non_inv_scaffold = current_total - current_two_kaliski;
        let scaled_by_div = 2_046.0 * 560.0;
        let branch_decode_margin = 150_000.0;
        let projected = non_inv_scaffold + scaled_by_div + branch_decode_margin;
        eprintln!(
            "scaled-BY DIV point-add budget: scaffold≈{non_inv_scaffold:.0}, div≈{scaled_by_div:.0}, margin≈{branch_decode_margin:.0}, projected≈{projected:.0}"
        );
        assert!(projected < 2_700_000.0, "scaled BY DIV would not beat Google low-qubit Toffoli target");
        assert!(projected < current_total - 1_500_000.0, "scaled BY DIV lacks inversion-sized saving");
    }

    #[test]
    fn branch_pattern_entropy_supports_compressed_history_target() {
        // Raw branch history is 560 bits. Encoding each 16-step window as a
        // branch pattern gives a concrete compression target that is closer to
        // a reversible implementation than matrix IDs: the pattern itself is
        // the control microprogram for scaled replay.
        use std::collections::HashMap;
        const W: usize = 16;
        const WINDOWS: usize = 35;
        let samples = 10_000usize;
        let mut sampler = Sampler::new(b"by-branch-pattern-entropy-v1", SECP256K1_P);
        let mut seqs: Vec<Vec<u16>> = Vec::with_capacity(samples);
        let mut counts: Vec<HashMap<u16, usize>> = (0..WINDOWS).map(|_| HashMap::new()).collect();
        for _ in 0..samples {
            let x = sampler.next();
            let mut delta = 1i64;
            let mut f = SInt::from_u(SECP256K1_P);
            let mut g = SInt::from_u(x);
            let mut seq = Vec::with_capacity(WINDOWS);
            for j in 0..WINDOWS {
                let mut pat = 0u16;
                for i in 0..W {
                    if g.bit0() {
                        pat |= 1u16 << i;
                    }
                    divstep_sint_state(&mut delta, &mut f, &mut g);
                }
                *counts[j].entry(pat).or_insert(0) += 1;
                seq.push(pat);
            }
            seqs.push(seq);
        }
        let mut entropy_sum = 0.0f64;
        let mut fixed_bits = 0usize;
        for c in &counts {
            fixed_bits += ((c.len() + 1) as f64).log2().ceil() as usize;
            for &n in c.values() {
                let p = n as f64 / samples as f64;
                entropy_sum -= p * p.log2();
            }
        }
        let mut code_lengths = Vec::with_capacity(samples);
        for seq in &seqs {
            let mut len = 0.0f64;
            for (j, &pat) in seq.iter().enumerate() {
                let n = *counts[j].get(&pat).unwrap();
                let p = n as f64 / samples as f64;
                len -= p.log2();
            }
            code_lengths.push(len);
        }
        code_lengths.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let p99 = code_lengths[samples * 99 / 100];
        let p999 = code_lengths[samples * 999 / 1000];
        let fail_520 = code_lengths.iter().filter(|&&l| l > 520.0).count() as f64 / samples as f64;
        eprintln!(
            "BY branch-pattern entropy: H≈{entropy_sum:.1} bits, p99≈{p99:.1}, p999≈{p999:.1}, fail>520≈{fail_520:.4}, fixed_distinct_bits={fixed_bits}"
        );
        assert!(entropy_sum < 500.0, "branch-pattern entropy too high for compressed history");
        assert!(p99 < 510.0, "branch-pattern p99 too high for 600-scratch target");
        assert!(fixed_bits < 560, "fixed per-window pattern IDs do not compress raw history");
    }

    #[test]
    fn branch_pattern_fixed_id_decoder_floor_blocks_generic_compression() {
        // The entropy result above is only a storage target.  If each window is
        // stored as a fixed dictionary ID, replay still has to recover the 16
        // odd-pattern controls.  A generic row decoder has to at least touch
        // every distinct pattern row in the codebook, before equality tests,
        // output-bit toggles, cleanup, or a second replay are charged.  That
        // already exceeds the remaining margin of the two-fast-replay budget,
        // so the compressed BY route needs a structured parser/consumer, not a
        // plain fixed-ID QROM decode.
        use std::collections::HashMap;
        const W: usize = 16;
        const WINDOWS: usize = 35;
        let samples = 10_000usize;
        let mut sampler =
            Sampler::new(b"by-branch-pattern-fixed-id-decoder-v1", SECP256K1_P);
        let mut counts: Vec<HashMap<u16, usize>> =
            (0..WINDOWS).map(|_| HashMap::new()).collect();
        for _ in 0..samples {
            let x = sampler.next();
            let mut delta = 1i64;
            let mut f = SInt::from_u(SECP256K1_P);
            let mut g = SInt::from_u(x);
            for j in 0..WINDOWS {
                let mut pat = 0u16;
                for i in 0..W {
                    if g.bit0() {
                        pat |= 1u16 << i;
                    }
                    divstep_sint_state(&mut delta, &mut f, &mut g);
                }
                *counts[j].entry(pat).or_insert(0) += 1;
            }
        }

        let fixed_id_bits = counts
            .iter()
            .map(|c| ((c.len() + 1) as f64).log2().ceil() as usize)
            .sum::<usize>();
        let distinct_rows = counts.iter().map(|c| c.len()).sum::<usize>();
        let nonzero_pattern_bits = counts
            .iter()
            .flat_map(|c| c.keys())
            .map(|pat| pat.count_ones() as usize)
            .sum::<usize>();
        let max_window_rows = counts.iter().map(|c| c.len()).max().unwrap_or(0);
        let two_fast_replay_with_pattern_decode = 2_637_286usize;
        let remaining_to_2700k = 2_700_000usize - two_fast_replay_with_pattern_decode;
        let row_floor_gap =
            two_fast_replay_with_pattern_decode as isize + distinct_rows as isize
                - 2_700_000isize;
        let bit_floor_gap =
            two_fast_replay_with_pattern_decode as isize + nonzero_pattern_bits as isize
                - 2_700_000isize;
        println!("METRIC by_branch_pattern_fixed_id_samples={samples}");
        println!("METRIC by_branch_pattern_fixed_id_bits={fixed_id_bits}");
        println!("METRIC by_branch_pattern_fixed_id_distinct_rows={distinct_rows}");
        println!("METRIC by_branch_pattern_fixed_id_max_window_rows={max_window_rows}");
        println!("METRIC by_branch_pattern_fixed_id_nonzero_table_bits={nonzero_pattern_bits}");
        println!("METRIC by_branch_pattern_fixed_id_remaining_to_2700k={remaining_to_2700k}");
        println!("METRIC by_branch_pattern_fixed_id_row_floor_gap_ccx={row_floor_gap}");
        println!("METRIC by_branch_pattern_fixed_id_bit_floor_gap_ccx={bit_floor_gap}");
        eprintln!(
            "BY fixed-ID pattern decoder floor: fixed_bits={fixed_id_bits}, rows={distinct_rows}, max_window_rows={max_window_rows}, nonzero_bits={nonzero_pattern_bits}, remaining={remaining_to_2700k}, row_gap={row_floor_gap}, bit_gap={bit_floor_gap}"
        );
        assert!(
            fixed_id_bits < 560,
            "fixed pattern IDs no longer compress raw history; update the entropy target"
        );
        assert!(
            distinct_rows > remaining_to_2700k,
            "generic fixed-ID row decoder now fits the two-replay margin; build it"
        );
        assert!(
            bit_floor_gap > row_floor_gap,
            "pattern-bit table floor should be stricter than row-touch floor"
        );
    }

    #[test]
    fn raw_pattern_streaming_parser_still_needs_reversible_delta_history() {
        // Escape hatch checked after the fixed-ID decoder miss: keep the raw
        // 560 odd bits and stream only one A control at a time.  The scratch
        // arithmetic is tempting, but the delta transition for an odd step is
        // not injective without an A/history bit: delta=1 and delta=-1 both map
        // to delta'=0.  A separate checkpoint-width probe finds that even the
        // sampled 5-bit delta checkpoint breaks the 663-bit Google scratch
        // model, and the current exact clean pattern decoder already consumes
        // the whole two-fast-replay margin before any compressed-ID expansion.
        use std::collections::BTreeMap;

        let mut images: BTreeMap<i64, Vec<i64>> = BTreeMap::new();
        for delta in -20i64..=20 {
            let next = if delta > 0 { 1 - delta } else { 1 + delta };
            images.entry(next).or_default().push(delta);
        }
        let collision_images = images.values().filter(|srcs| srcs.len() > 1).count();
        let worst_collision = images.values().map(Vec::len).max().unwrap_or(0);
        let zero_sources = images.get(&0).cloned().unwrap_or_default();

        const GOOGLE_LOW_QUBIT_SCRATCH: usize = 663;
        const STRICT_SCRATCH: usize = 600;
        const RAW_PATTERN_BITS: usize = 560;
        const DELTA_BITS: usize = 10;
        const SINGLE_A_BIT: usize = 1;
        const WINDOW_A_BITS: usize = 16;
        const TARGET_LOCAL_WORKSPACE: usize = 90;
        const EXACT_PATTERN_DECODER_PER_REPLAY: usize = 62_160;
        const TWO_FAST_REPLAY_BEFORE_BRANCH_DECODE: usize = 2_577_286;

        let raw_single_a_scratch =
            RAW_PATTERN_BITS + DELTA_BITS + SINGLE_A_BIT + TARGET_LOCAL_WORKSPACE;
        let raw_with_one_delta_checkpoint = raw_single_a_scratch + DELTA_BITS;
        let raw_window_a_scratch =
            RAW_PATTERN_BITS + DELTA_BITS + WINDOW_A_BITS + TARGET_LOCAL_WORKSPACE;
        let raw_full_a_history_scratch =
            RAW_PATTERN_BITS + DELTA_BITS + RAW_PATTERN_BITS + TARGET_LOCAL_WORKSPACE;
        let exact_two_decoder_projected =
            TWO_FAST_REPLAY_BEFORE_BRANCH_DECODE + 2 * EXACT_PATTERN_DECODER_PER_REPLAY;
        let exact_two_decoder_gap = exact_two_decoder_projected as isize - 2_700_000isize;

        println!("METRIC by_raw_pattern_single_a_scratch={raw_single_a_scratch}");
        println!("METRIC by_raw_pattern_one_checkpoint_scratch={raw_with_one_delta_checkpoint}");
        println!("METRIC by_raw_pattern_window_a_scratch={raw_window_a_scratch}");
        println!("METRIC by_raw_pattern_full_a_history_scratch={raw_full_a_history_scratch}");
        println!("METRIC by_raw_pattern_delta_collision_images={collision_images}");
        println!("METRIC by_raw_pattern_delta_worst_collision={worst_collision}");
        println!("METRIC by_raw_pattern_exact_two_decoder_projected={exact_two_decoder_projected}");
        println!("METRIC by_raw_pattern_exact_two_decoder_gap_to_2700k={exact_two_decoder_gap}");
        eprintln!(
            "BY raw-pattern streaming parser: singleA_scratch={raw_single_a_scratch}, one_checkpoint={raw_with_one_delta_checkpoint}, windowA={raw_window_a_scratch}, fullA={raw_full_a_history_scratch}, zero_sources={zero_sources:?}, exact_two_decoder_projected={exact_two_decoder_projected}, gap={exact_two_decoder_gap}"
        );

        assert!(
            raw_single_a_scratch <= GOOGLE_LOW_QUBIT_SCRATCH,
            "single-A raw-pattern schedule is no longer the relevant tempting case"
        );
        assert!(
            raw_single_a_scratch > STRICT_SCRATCH,
            "raw pattern bits unexpectedly fit strict 600 scratch without compression"
        );
        assert!(
            zero_sources.contains(&-1) && zero_sources.contains(&1),
            "odd-step delta update no longer exhibits the expected local collision"
        );
        assert!(
            collision_images > 0 && worst_collision > 1,
            "delta update became locally injective; revisit streaming parser"
        );
        assert!(
            raw_with_one_delta_checkpoint > GOOGLE_LOW_QUBIT_SCRATCH,
            "one delta checkpoint now fits raw-pattern Google scratch; revisit raw streaming"
        );
        assert!(
            raw_window_a_scratch > GOOGLE_LOW_QUBIT_SCRATCH,
            "window-local A scratch now fits raw-pattern Google scratch"
        );
        assert!(
            exact_two_decoder_projected > 2_700_000,
            "current exact pattern decoder now leaves compressed-ID margin"
        );
    }

    #[test]
    fn raw_pattern_delta_checkpoint_exceeds_google_scratch_slack() {
        // The single-A raw-pattern schedule has only two scratch bits left
        // under the 663-bit Google model.  A reversible streaming parser needs
        // either enough old-delta checkpoint state to reverse a window, or it
        // must retain the ambiguous A decisions.  Sample the actual 35-window
        // BY trajectories: the window-start delta alphabet is already wider
        // than two bits, and the per-step A ambiguity is far larger still.
        use std::collections::HashMap;

        const W: usize = 16;
        const WINDOWS: usize = 35;
        const GOOGLE_LOW_QUBIT_SCRATCH: usize = 663;
        const RAW_SINGLE_A_SCRATCH: usize = 661;
        let samples = 10_000usize;
        let mut sampler = Sampler::new(b"by-raw-pattern-delta-checkpoint-v1", SECP256K1_P);
        let mut window_delta_counts: Vec<HashMap<i64, usize>> =
            (0..WINDOWS).map(|_| HashMap::new()).collect();
        let mut ambiguous_a_counts = Vec::with_capacity(samples);
        for _ in 0..samples {
            let x = sampler.next();
            let mut delta = 1i64;
            let mut f = SInt::from_u(SECP256K1_P);
            let mut g = SInt::from_u(x);
            let mut ambiguous = 0usize;
            for counts in &mut window_delta_counts {
                *counts.entry(delta).or_insert(0) += 1;
                for _ in 0..W {
                    let old_delta = delta;
                    let odd = g.bit0();
                    divstep_sint_state(&mut delta, &mut f, &mut g);
                    if odd && delta <= 0 && old_delta != 0 {
                        ambiguous += 1;
                    }
                }
            }
            ambiguous_a_counts.push(ambiguous);
        }
        ambiguous_a_counts.sort_unstable();
        let max_window_delta_rows = window_delta_counts
            .iter()
            .map(HashMap::len)
            .max()
            .unwrap_or(0);
        let checkpoint_bits =
            usize::BITS as usize - max_window_delta_rows.saturating_sub(1).leading_zeros() as usize;
        let checkpoint_scratch = RAW_SINGLE_A_SCRATCH + checkpoint_bits;
        let scratch_slack = GOOGLE_LOW_QUBIT_SCRATCH - RAW_SINGLE_A_SCRATCH;
        let ambiguous_mean =
            ambiguous_a_counts.iter().sum::<usize>() as f64 / ambiguous_a_counts.len() as f64;
        let ambiguous_p99 = ambiguous_a_counts[ambiguous_a_counts.len() * 99 / 100];
        let ambiguous_max = *ambiguous_a_counts.last().unwrap();
        println!("METRIC by_raw_pattern_delta_checkpoint_samples={samples}");
        println!("METRIC by_raw_pattern_delta_checkpoint_max_rows={max_window_delta_rows}");
        println!("METRIC by_raw_pattern_delta_checkpoint_bits={checkpoint_bits}");
        println!("METRIC by_raw_pattern_delta_checkpoint_scratch={checkpoint_scratch}");
        println!("METRIC by_raw_pattern_delta_checkpoint_scratch_slack={scratch_slack}");
        println!("METRIC by_raw_pattern_ambiguous_a_bits_mean={ambiguous_mean:.3}");
        println!("METRIC by_raw_pattern_ambiguous_a_bits_p99={ambiguous_p99}");
        println!("METRIC by_raw_pattern_ambiguous_a_bits_max={ambiguous_max}");
        eprintln!(
            "BY raw-pattern delta checkpoint: max_rows={max_window_delta_rows}, bits={checkpoint_bits}, scratch={checkpoint_scratch}, slack={scratch_slack}, ambiguous_a_mean={ambiguous_mean:.1}, p99={ambiguous_p99}, max={ambiguous_max}"
        );
        assert!(
            checkpoint_bits > scratch_slack,
            "delta checkpoint now fits the raw-pattern scratch slack; revisit streaming parser"
        );
        assert!(
            checkpoint_scratch > GOOGLE_LOW_QUBIT_SCRATCH,
            "raw-pattern delta checkpoint now fits Google scratch"
        );
        assert!(
            ambiguous_p99 > 32,
            "ambiguous A history became sparse enough to revisit retained-A cleanup"
        );
    }

    #[test]
    fn actual_matrix_sequence_entropy_supports_sub600_history_target() {
        // Storing raw 22-bit (delta,h) keys costs 770 bits for 35 windows, but
        // actual secp256k1 trajectories are highly non-uniform, especially near
        // convergence. An entropy-coded matrix history is not a circuit yet, but
        // this shows the information-theoretic target is below the user's
        // ~600-scratch budget.
        use std::collections::HashMap;
        const W: usize = 16;
        const WINDOWS: usize = 35;
        let samples = 10_000usize;
        let mut sampler = Sampler::new(b"by-matrix-sequence-entropy-v1", SECP256K1_P);
        let mut seqs: Vec<Vec<TransitionMatrix>> = Vec::with_capacity(samples);
        let mut counts: Vec<HashMap<TransitionMatrix, usize>> = (0..WINDOWS).map(|_| HashMap::new()).collect();
        for _ in 0..samples {
            let x = sampler.next();
            let mut delta = 1i64;
            let mut f = SInt::from_u(SECP256K1_P);
            let mut g = SInt::from_u(x);
            let mut seq = Vec::with_capacity(WINDOWS);
            for j in 0..WINDOWS {
                let f_low = sint_low_i128(f, W);
                let g_low = sint_low_i128(g, W);
                let (_, _, _, mtx) = jump_matrix_direct_lowword(W, W, delta, f_low, g_low);
                *counts[j].entry(mtx).or_insert(0) += 1;
                seq.push(mtx);
                for _ in 0..W {
                    divstep_sint_state(&mut delta, &mut f, &mut g);
                }
            }
            seqs.push(seq);
        }
        let mut entropy_sum = 0.0f64;
        for c in &counts {
            for &n in c.values() {
                let p = n as f64 / samples as f64;
                entropy_sum -= p * p.log2();
            }
        }
        let mut code_lengths = Vec::with_capacity(samples);
        for seq in &seqs {
            let mut len = 0.0f64;
            for (j, mtx) in seq.iter().enumerate() {
                let n = *counts[j].get(mtx).unwrap();
                let p = n as f64 / samples as f64;
                len -= p.log2();
            }
            code_lengths.push(len);
        }
        code_lengths.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let p99 = code_lengths[(samples * 99) / 100];
        let p999 = code_lengths[(samples * 999) / 1000];
        let fail_550 = code_lengths.iter().filter(|&&l| l > 550.0).count() as f64 / samples as f64;
        eprintln!(
            "BY matrix-sequence entropy: H≈{entropy_sum:.1} bits, p99_len≈{p99:.1}, p999_len≈{p999:.1}, fail>550≈{fail_550:.4}"
        );
        assert!(entropy_sum < 520.0, "matrix history entropy too high for sub600 target");
        assert!(p99 < 540.0, "p99 matrix history code length too high");
        assert!(fail_550 < 0.01, "550-bit matrix history would exceed 1% failure tolerance");
    }

    fn toy_by_matrix_sequence(n: usize, p: u16, x: u16, w: usize) -> Vec<TransitionMatrix> {
        let windows = (2 * n).div_ceil(w);
        let mut delta = 1i64;
        let mut f = SInt::from_u(U256::from(p as u64));
        let mut g = SInt::from_u(U256::from(x as u64));
        let mut seq = Vec::with_capacity(windows);
        for _ in 0..windows {
            let f_low = sint_low_i128(f, w);
            let g_low = sint_low_i128(g, w);
            let (_, _, _, mtx) = jump_matrix_direct_lowword(w, w, delta, f_low, g_low);
            seq.push(mtx);
            for _ in 0..w {
                divstep_sint_state(&mut delta, &mut f, &mut g);
            }
        }
        seq
    }

    fn by_matrix_sequence_rank_parity_anf_stats_for_test(n: usize, p: u16, w: usize) -> (usize, usize, usize, usize) {
        use std::collections::{BTreeMap, BTreeSet};
        let size = 1usize << n;
        let mut seqs = Vec::with_capacity(p as usize - 1);
        let mut distinct = BTreeSet::<Vec<TransitionMatrix>>::new();
        for x in 1..p {
            let seq = toy_by_matrix_sequence(n, p, x, w);
            distinct.insert(seq.clone());
            seqs.push((x as usize, seq));
        }
        let distinct_count = distinct.len();
        let rank_bits = if distinct_count <= 1 {
            0
        } else {
            usize::BITS as usize - (distinct_count - 1).leading_zeros() as usize
        };
        let ranks: BTreeMap<Vec<TransitionMatrix>, usize> = distinct
            .into_iter()
            .enumerate()
            .map(|(rank, seq)| (seq, rank))
            .collect();
        let mut anf = vec![0u8; size];
        for (x, seq) in seqs {
            anf[x] = (ranks.get(&seq).expect("ranked matrix sequence") & 1) as u8;
        }
        for bit in 0..n {
            for idx in 0..size {
                if (idx & (1usize << bit)) != 0 {
                    anf[idx] ^= anf[idx ^ (1usize << bit)];
                }
            }
        }
        let density = anf.iter().filter(|&&c| c != 0).count();
        let degree = anf
            .iter()
            .enumerate()
            .filter_map(|(i, &c)| if c != 0 { Some(i.count_ones() as usize) } else { None })
            .max()
            .unwrap_or(0);
        (degree, density, distinct_count, rank_bits)
    }

    #[test]
    fn by_matrix_sequence_rank_decoder_is_dense_on_toys() {
        // The matrix-history entropy result is only information-theoretic.
        // If the compressed payload is represented by a canonical rank of the
        // exact matrix sequence, even one rank bit is already a high-degree,
        // dense global function of the initial denominator on exhaustive toy
        // fields. That keeps the compressed BY route in the "needs a real
        // reversible parser/consumer" bucket rather than a cheap MBUC cleanup.
        let cases = [(8usize, 251u16), (10, 1021), (12, 4093), (14, 16381)];
        for &(n, p) in &cases {
            let (degree, density, distinct, rank_bits) =
                by_matrix_sequence_rank_parity_anf_stats_for_test(n, p, 4);
            let table = 1usize << n;
            eprintln!(
                "BY matrix-sequence rank ANF: n={n}, degree={degree}, density={density}/{table}, distinct={distinct}, rank_bits={rank_bits}"
            );
            if n == 14 {
                println!("METRIC by_matrix_sequence_rank_degree_n14={degree}");
                println!("METRIC by_matrix_sequence_rank_density_n14={density}");
                println!("METRIC by_matrix_sequence_rank_distinct_n14={distinct}");
                println!("METRIC by_matrix_sequence_rank_bits_n14={rank_bits}");
            }
            assert!(degree + 1 >= n);
            assert!(density > table / 4);
            assert!(rank_bits + 1 >= n);
        }
    }

    fn mat2_mul_i128(a: [[i128; 2]; 2], b: [[i128; 2]; 2]) -> [[i128; 2]; 2] {
        [
            [a[0][0] * b[0][0] + a[0][1] * b[1][0], a[0][0] * b[0][1] + a[0][1] * b[1][1]],
            [a[1][0] * b[0][0] + a[1][1] * b[1][0], a[1][0] * b[0][1] + a[1][1] * b[1][1]],
        ]
    }

    fn mat2_det_i128(a: [[i128; 2]; 2]) -> i128 {
        a[0][0] * a[1][1] - a[0][1] * a[1][0]
    }

    fn mat2_max_abs_i128(a: [[i128; 2]; 2]) -> i128 {
        a.iter().flatten().map(|x| x.abs()).max().unwrap_or(0)
    }

    fn mat2_inv_unimodular_for_test(a: [[i128; 2]; 2]) -> [[i128; 2]; 2] {
        let det = mat2_det_i128(a);
        assert!(det == 1 || det == -1, "not unimodular: det={det}");
        let s = det;
        [[s * a[1][1], -s * a[0][1]], [-s * a[1][0], s * a[0][0]]]
    }

    fn snf2_for_test(input: [[i128; 2]; 2]) -> ([[i128; 2]; 2], [[i128; 2]; 2], [[i128; 2]; 2]) {
        // Returns (U,D,V) with U*input*V = D. This tiny 2x2 Smith-normal-form
        // helper is test-only; it lets us reason about in-place scaled BY
        // windows without introducing a dependency on a CAS.
        let mut m = input;
        let mut u = [[1i128, 0], [0, 1]];
        let mut v = [[1i128, 0], [0, 1]];
        for _ in 0..10_000 {
            if m[0][0] == 0 {
                let mut pos = None;
                for i in 0..2 {
                    for j in 0..2 {
                        if m[i][j] != 0 {
                            pos = Some((i, j));
                            break;
                        }
                    }
                    if pos.is_some() {
                        break;
                    }
                }
                let Some((i, j)) = pos else { return (u, m, v); };
                if i != 0 {
                    m.swap(0, i);
                    u.swap(0, i);
                }
                if j != 0 {
                    for r in 0..2 {
                        m[r].swap(0, j);
                        v[r].swap(0, j);
                    }
                }
            }

            let mut changed = false;
            if m[1][0] != 0 {
                let q = m[1][0] / m[0][0];
                for c in 0..2 {
                    m[1][c] -= q * m[0][c];
                    u[1][c] -= q * u[0][c];
                }
                changed = true;
                if m[1][0] != 0 && m[1][0].abs() < m[0][0].abs() {
                    m.swap(0, 1);
                    u.swap(0, 1);
                }
            }
            if m[0][1] != 0 {
                let q = m[0][1] / m[0][0];
                for r in 0..2 {
                    m[r][1] -= q * m[r][0];
                    v[r][1] -= q * v[r][0];
                }
                changed = true;
                if m[0][1] != 0 && m[0][1].abs() < m[0][0].abs() {
                    for r in 0..2 {
                        m[r].swap(0, 1);
                        v[r].swap(0, 1);
                    }
                }
            }
            if changed {
                continue;
            }
            if m[1][0] != 0 {
                assert_eq!(m[1][0] % m[0][0], 0);
                let q = m[1][0] / m[0][0];
                for c in 0..2 {
                    m[1][c] -= q * m[0][c];
                    u[1][c] -= q * u[0][c];
                }
                continue;
            }
            if m[0][1] != 0 {
                assert_eq!(m[0][1] % m[0][0], 0);
                let q = m[0][1] / m[0][0];
                for r in 0..2 {
                    m[r][1] -= q * m[r][0];
                    v[r][1] -= q * v[r][0];
                }
                continue;
            }
            if m[1][1] == 0 || m[1][1] % m[0][0] == 0 {
                for i in 0..2 {
                    if m[i][i] < 0 {
                        for c in 0..2 {
                            m[i][c] = -m[i][c];
                            u[i][c] = -u[i][c];
                        }
                    }
                }
                return (u, m, v);
            }
            // Mix the lower-right entry back into the pivot block and keep
            // reducing until the diagonal divisibility condition holds.
            for r in 0..2 {
                m[r][0] += m[r][1];
                v[r][0] += v[r][1];
            }
        }
        panic!("2x2 SNF did not converge for {input:?}");
    }

    #[test]
    fn smith_factorization_reduces_by_window_to_inplace_shifts_and_unimodular_maps() {
        const W: usize = 16;
        let mut hasher = sha3::Shake128::default();
        hasher.update(b"by-smith-factorization-v1");
        let mut reader = hasher.finalize_xof();
        let mut buf = [0u8; 24];
        let mut max_uv = 0i128;
        let mut max_uvi = 0i128;
        let mut diag_hist = std::collections::BTreeMap::<(i128, i128), usize>::new();
        for _ in 0..4096 {
            reader.read(&mut buf);
            let f_low = (u64::from_le_bytes(buf[0..8].try_into().unwrap()) as i128) | 1;
            let g_low = u64::from_le_bytes(buf[8..16].try_into().unwrap()) as i128;
            let delta = (u64::from_le_bytes(buf[16..24].try_into().unwrap()) % 41) as i64 - 20;
            let (_, _, _, pmat) = jump_matrix_direct_lowword(W, W, delta, f_low, g_low);
            let p = [[pmat.m00, pmat.m01], [pmat.m10, pmat.m11]];
            let (u, d, v) = snf2_for_test(p);
            assert_eq!(mat2_mul_i128(mat2_mul_i128(u, p), v), d);
            assert_eq!(d[0][1], 0);
            assert_eq!(d[1][0], 0);
            assert_eq!(d[0][0] * d[1][1], 1i128 << W);
            assert!((d[0][0] as u128).is_power_of_two());
            assert!((d[1][1] as u128).is_power_of_two());
            *diag_hist.entry((d[0][0], d[1][1])).or_default() += 1;
            let ui = mat2_inv_unimodular_for_test(u);
            let vi = mat2_inv_unimodular_for_test(v);
            max_uv = max_uv.max(mat2_max_abs_i128(u)).max(mat2_max_abs_i128(v));
            max_uvi = max_uvi.max(mat2_max_abs_i128(ui)).max(mat2_max_abs_i128(vi));
        }
        eprintln!(
            "BY Smith factorization: diag_hist={diag_hist:?}, max_UV_entry={max_uv}, max_inverse_entry={max_uvi}"
        );
        assert!(max_uvi > (1i128 << W), "naive SNF unexpectedly gave uniformly small factors");
    }

    fn gcd_i128_for_test(mut a: i128, mut b: i128) -> i128 {
        a = a.abs();
        b = b.abs();
        while b != 0 {
            let r = a % b;
            a = b;
            b = r;
        }
        a
    }

    fn egcd_i128_for_test(a: i128, b: i128) -> (i128, i128, i128) {
        if b == 0 {
            let g = a.abs();
            let x = if a < 0 { -1 } else { 1 };
            return (x, 0, g);
        }
        let (x1, y1, g) = egcd_i128_for_test(b, a % b);
        (y1, x1 - (a / b) * y1, g)
    }

    fn centered_mod_i128(x: i128, modulus: i128) -> i128 {
        let m = modulus.abs();
        ((x + m / 2).rem_euclid(m)) - m / 2
    }

    fn hermite_scaled_window_factor_for_test(
        p: [[i128; 2]; 2],
    ) -> Option<([[i128; 2]; 2], [[i128; 2]; 2], [[i128; 2]; 2])> {
        // Find small U,V,H with U*P*V = H = [[1,e],[0,65536]], |e|<=32768.
        // Then P/65536 = U^-1 * [[2^-16, e*2^-16],[0,1]] * V^-1,
        // i.e. an in-place scaled window can be implemented as two unimodular
        // maps plus a single batched divide-by-2^16 row.
        let mut best = None;
        for radius in [2i128, 8] {
            for r in -radius..=radius {
                for s in -radius..=radius {
                    if r == 0 && s == 0 {
                        continue;
                    }
                    if gcd_i128_for_test(r, s) != 1 {
                        continue;
                    }
                    let y0 = p[0][0] * r + p[0][1] * s;
                    let y1 = p[1][0] * r + p[1][1] * s;
                    if gcd_i128_for_test(y0, y1) != 1 {
                        continue;
                    }
                    let (alpha, beta, gy) = egcd_i128_for_test(y0, y1);
                    if gy != 1 {
                        continue;
                    }
                    let mut u = [[alpha, beta], [-y1, y0]];
                    let (aa, bb, grs) = egcd_i128_for_test(r, s);
                    if grs != 1 {
                        continue;
                    }
                    let mut v = [[r, -bb], [s, aa]];
                    let mut h = mat2_mul_i128(mat2_mul_i128(u, p), v);
                    if h[0][0] == -1 && h[1][0] == 0 {
                        for c in 0..2 {
                            u[0][c] = -u[0][c];
                        }
                        h = mat2_mul_i128(mat2_mul_i128(u, p), v);
                    }
                    if h[0][0] != 1 || h[1][0] != 0 || h[0][1].abs() > (1i128 << 80) {
                        continue;
                    }
                    if h[1][1] < 0 {
                        for c in 0..2 {
                            u[1][c] = -u[1][c];
                        }
                        h = mat2_mul_i128(mat2_mul_i128(u, p), v);
                    }
                    if h[1][1] != (1i128 << 16) {
                        continue;
                    }
                    let e_reduced = centered_mod_i128(h[0][1], h[1][1]);
                    let k = e_reduced - h[0][1];
                    for row in 0..2 {
                        v[row][1] += k * v[row][0];
                    }
                    h = mat2_mul_i128(mat2_mul_i128(u, p), v);
                    if h[0][0] != 1 || h[1][0] != 0 || h[1][1] != (1i128 << 16) {
                        continue;
                    }
                    if h[0][1].abs() > (1i128 << 15) {
                        continue;
                    }
                    let score = mat2_max_abs_i128(u).max(mat2_max_abs_i128(v)).max(h[0][1].abs());
                    if best.as_ref().map_or(true, |(best_score, _, _, _): &(i128, _, _, _)| score < *best_score) {
                        best = Some((score, u, h, v));
                    }
                }
            }
            if best.is_some() {
                break;
            }
        }
        best.map(|(_, u, h, v)| (u, h, v))
    }

    #[derive(Clone, Copy, Debug)]
    enum RowOp2ForTest {
        Add { dst: usize, src: usize, k: i128 },
        Swap,
        Neg { row: usize },
    }

    fn apply_row_op_to_mat_for_test(m: &mut [[i128; 2]; 2], op: RowOp2ForTest) {
        match op {
            RowOp2ForTest::Add { dst, src, k } => {
                for c in 0..2 {
                    m[dst][c] += k * m[src][c];
                }
            }
            RowOp2ForTest::Swap => m.swap(0, 1),
            RowOp2ForTest::Neg { row } => {
                for c in 0..2 {
                    m[row][c] = -m[row][c];
                }
            }
        }
    }

    fn inverse_row_op_for_test(op: RowOp2ForTest) -> RowOp2ForTest {
        match op {
            RowOp2ForTest::Add { dst, src, k } => RowOp2ForTest::Add { dst, src, k: -k },
            RowOp2ForTest::Swap => RowOp2ForTest::Swap,
            RowOp2ForTest::Neg { row } => RowOp2ForTest::Neg { row },
        }
    }

    fn reduce_unimodular_to_identity_ops_for_test(mut m: [[i128; 2]; 2]) -> Vec<RowOp2ForTest> {
        let det = mat2_det_i128(m);
        assert!(det == 1 || det == -1, "not unimodular: {m:?}, det={det}");
        let mut ops = Vec::new();
        for _ in 0..256 {
            let a = m[0][0];
            let c = m[1][0];
            if c == 0 {
                break;
            }
            if a.abs() >= c.abs() {
                let q = a / c;
                assert_ne!(q, 0, "Euclid quotient unexpectedly zero for {m:?}");
                let op = RowOp2ForTest::Add { dst: 0, src: 1, k: -q };
                apply_row_op_to_mat_for_test(&mut m, op);
                ops.push(op);
            } else {
                let op = RowOp2ForTest::Swap;
                apply_row_op_to_mat_for_test(&mut m, op);
                ops.push(op);
            }
        }
        assert_eq!(m[1][0], 0, "Euclid reduction failed: {m:?}");
        assert!(m[0][0] == 1 || m[0][0] == -1, "bad pivot after reduction: {m:?}");
        if m[0][0] == -1 {
            let op = RowOp2ForTest::Neg { row: 0 };
            apply_row_op_to_mat_for_test(&mut m, op);
            ops.push(op);
        }
        assert_eq!(m[0][0], 1);
        let d = m[1][1];
        assert!(d == 1 || d == -1, "bad lower diagonal after reduction: {m:?}");
        if m[0][1] != 0 {
            let k = -m[0][1] / d;
            let op = RowOp2ForTest::Add { dst: 0, src: 1, k };
            apply_row_op_to_mat_for_test(&mut m, op);
            ops.push(op);
        }
        if m[1][1] == -1 {
            let op = RowOp2ForTest::Neg { row: 1 };
            apply_row_op_to_mat_for_test(&mut m, op);
            ops.push(op);
        }
        assert_eq!(m, [[1, 0], [0, 1]], "did not reduce to identity");
        ops
    }

    fn emit_mod_shear_small_coeff_for_test(
        b: &mut super::super::B,
        dst: &[super::super::QubitId],
        src: &[super::super::QubitId],
        k: i128,
        p: U256,
    ) {
        if k == 0 {
            return;
        }
        mod_mul_two_small_coeffs_acc_for_cost(b, src, k, dst, 0, dst, p);
    }

    fn emit_row_op_mod_for_test(
        b: &mut super::super::B,
        x0: &[super::super::QubitId],
        x1: &[super::super::QubitId],
        op: RowOp2ForTest,
        p: U256,
    ) {
        match op {
            RowOp2ForTest::Add { dst: 0, src: 1, k } => emit_mod_shear_small_coeff_for_test(b, x0, x1, k, p),
            RowOp2ForTest::Add { dst: 1, src: 0, k } => emit_mod_shear_small_coeff_for_test(b, x1, x0, k, p),
            RowOp2ForTest::Add { .. } => unreachable!("invalid 2-row op"),
            RowOp2ForTest::Swap => {
                for i in 0..x0.len() {
                    b.swap(x0[i], x1[i]);
                }
            }
            RowOp2ForTest::Neg { row: 0 } => super::super::mod_neg_inplace_fast(b, x0, p),
            RowOp2ForTest::Neg { row: 1 } => super::super::mod_neg_inplace_fast(b, x1, p),
            RowOp2ForTest::Neg { .. } => unreachable!("invalid 2-row op"),
        }
    }

    fn emit_unimodular_matrix_mod_inplace_for_test(
        b: &mut super::super::B,
        m: [[i128; 2]; 2],
        x0: &[super::super::QubitId],
        x1: &[super::super::QubitId],
        p: U256,
    ) -> (usize, i128) {
        let ops = reduce_unimodular_to_identity_ops_for_test(m);
        let max_k = ops
            .iter()
            .filter_map(|op| match op {
                RowOp2ForTest::Add { k, .. } => Some(k.abs()),
                _ => None,
            })
            .max()
            .unwrap_or(0);
        let count = ops.len();
        for op in ops.into_iter().rev().map(inverse_row_op_for_test) {
            emit_row_op_mod_for_test(b, x0, x1, op, p);
        }
        (count, max_k)
    }

    #[test]
    fn hermite_factorization_keeps_scaled_by_window_in_place_with_small_coefficients() {
        const W: usize = 16;
        let mut hasher = sha3::Shake128::default();
        hasher.update(b"by-hermite-factorization-v1");
        let mut reader = hasher.finalize_xof();
        let mut buf = [0u8; 24];
        let mut max_coeff = 0i128;
        let mut p99_scores = Vec::new();
        for _ in 0..4096 {
            reader.read(&mut buf);
            let f_low = (u64::from_le_bytes(buf[0..8].try_into().unwrap()) as i128) | 1;
            let g_low = u64::from_le_bytes(buf[8..16].try_into().unwrap()) as i128;
            let delta = (u64::from_le_bytes(buf[16..24].try_into().unwrap()) % 41) as i64 - 20;
            let (_, _, _, pmat) = jump_matrix_direct_lowword(W, W, delta, f_low, g_low);
            let p = [[pmat.m00, pmat.m01], [pmat.m10, pmat.m11]];
            let (u, h, v) = hermite_scaled_window_factor_for_test(p).expect("small Hermite factor");
            assert_eq!(mat2_mul_i128(mat2_mul_i128(u, p), v), h);
            assert_eq!(h[0][0], 1);
            assert_eq!(h[1][0], 0);
            assert_eq!(h[1][1], 1i128 << W);
            assert!(h[0][1].abs() <= (1i128 << 15));
            let ui = mat2_inv_unimodular_for_test(u);
            let vi = mat2_inv_unimodular_for_test(v);
            let score = mat2_max_abs_i128(u)
                .max(mat2_max_abs_i128(v))
                .max(mat2_max_abs_i128(ui))
                .max(mat2_max_abs_i128(vi))
                .max(h[0][1].abs());
            max_coeff = max_coeff.max(score);
            p99_scores.push(score);
        }
        p99_scores.sort_unstable();
        let p99 = p99_scores[p99_scores.len() * 99 / 100];
        eprintln!(
            "BY Hermite in-place factorization: p99_coeff={p99}, max_coeff={max_coeff}"
        );
        assert!(max_coeff <= (1i128 << W), "Hermite factors exceeded w-bit coefficient scale");
    }

    fn inv_pow2_mod_p_for_test(w: usize, p: U256) -> U256 {
        let inv2 = (p.wrapping_add(U256::from(1u64))) >> 1usize;
        let mut acc = U256::from(1u64);
        for _ in 0..w {
            acc = mulm(acc, inv2, p);
        }
        acc
    }

    fn emit_fixed_hermite_inplace_window_for_test(
        b: &mut super::super::B,
        pmat: TransitionMatrix,
        x0: &[super::super::QubitId],
        x1: &[super::super::QubitId],
        p_mod: U256,
    ) -> (usize, i128, i128) {
        const W: usize = 16;
        let p = [[pmat.m00, pmat.m01], [pmat.m10, pmat.m11]];
        let (u, h, v) = hermite_scaled_window_factor_for_test(p).expect("Hermite factor");
        assert_eq!(mat2_mul_i128(mat2_mul_i128(u, p), v), h);
        let ui = mat2_inv_unimodular_for_test(u);
        let vi = mat2_inv_unimodular_for_test(v);
        let e = h[0][1];
        let (v_ops, v_max_k) = emit_unimodular_matrix_mod_inplace_for_test(b, vi, x0, x1, p_mod);
        emit_mod_shear_small_coeff_for_test(b, x0, x1, e, p_mod);
        for _ in 0..W {
            super::super::mod_halve_inplace_fast(b, x0, p_mod);
        }
        let (u_ops, u_max_k) = emit_unimodular_matrix_mod_inplace_for_test(b, ui, x0, x1, p_mod);
        (v_ops + u_ops + 1, v_max_k.max(u_max_k).max(e.abs()), e)
    }

    #[test]
    fn fixed_hermite_inplace_modular_window_matches_scaled_by_matrix() {
        const W: usize = 16;
        let p_mod = SECP256K1_P;
        let pmat = jump_matrix_direct_lowword(W, W, 1, 1, 3).3;
        assert_eq!((pmat.m00, pmat.m01, pmat.m10, pmat.m11), (-8192, 24576, -3, 1));
        let p = [[pmat.m00, pmat.m01], [pmat.m10, pmat.m11]];
        let (u, h, v) = hermite_scaled_window_factor_for_test(p).expect("Hermite factor");
        assert_eq!(mat2_mul_i128(mat2_mul_i128(u, p), v), h);
        let e = h[0][1];
        assert_eq!(h[0][0], 1);
        assert_eq!(h[1][0], 0);
        assert_eq!(h[1][1], 1i128 << W);

        let mut b = super::super::B::new();
        let x0 = b.alloc_qubits(256);
        let x1 = b.alloc_qubits(256);
        let (factor_ops, max_shear, _) = emit_fixed_hermite_inplace_window_for_test(&mut b, pmat, &x0, &x1, p_mod);
        let ccx = count_ccx(&b.ops);
        let peak = b.peak_qubits;
        let num_qubits = b.next_qubit as usize;
        let num_bits = b.next_bit as usize;
        let ops = b.ops;

        let inv2w = inv_pow2_mod_p_for_test(W, p_mod);
        let row_expected = |a: U256, c: U256, c0: i128, c1: i128| -> U256 {
            let t0 = mulm(signed_i128_mod_p(c0, p_mod), a, p_mod);
            let t1 = mulm(signed_i128_mod_p(c1, p_mod), c, p_mod);
            mulm(addm(t0, t1, p_mod), inv2w, p_mod)
        };
        let mut sx = Sampler::new(b"by-hermite-inplace-x0-v1", p_mod);
        let mut sy = Sampler::new(b"by-hermite-inplace-x1-v1", p_mod);
        for _ in 0..32 {
            let a = sx.next();
            let c = sy.next();
            let exp0 = row_expected(a, c, pmat.m00, pmat.m01);
            let exp1 = row_expected(a, c, pmat.m10, pmat.m11);
            let mut hasher = sha3::Shake128::default();
            hasher.update(b"by-hermite-inplace-sim-v1");
            let mut xof = hasher.finalize_xof();
            let mut sim = crate::sim::Simulator::new(num_qubits, num_bits, &mut xof);
            set_slice_u512_by(&mut sim, &x0, u256_to_u512_for_by_tests(a));
            set_slice_u512_by(&mut sim, &x1, u256_to_u512_for_by_tests(c));
            sim.apply(&ops);
            assert_eq!(get_slice_u512_by(&sim, &x0), u256_to_u512_for_by_tests(exp0), "row0 mismatch");
            assert_eq!(get_slice_u512_by(&sim, &x1), u256_to_u512_for_by_tests(exp1), "row1 mismatch");
        }
        eprintln!(
            "BY fixed Hermite in-place modular window: ccx={ccx}, peak={peak}q, e={e}, factor_ops={factor_ops}, max_shear={max_shear}"
        );
        assert!(peak < 1_600, "in-place Hermite window lost scratch advantage");
        assert!(ccx < 80_000, "fixed Hermite in-place sample window too costly");
    }

    fn emit_fixed_branch_numerator_scaled_window_for_test(
        b: &mut super::super::B,
        mut delta: i64,
        bits: &[bool],
        x0: &[super::super::QubitId],
        x1: &[super::super::QubitId],
        p: U256,
    ) -> (usize, usize, usize) {
        // Apply the 16 numerator microstep matrices directly, then apply the
        // common 2^-16 scaling to both rows. This uses the branch bits as the
        // circuit description and avoids Hermite-factor synthesis entirely in
        // the fixed-control model.
        let mut a_cases = 0usize;
        let mut b_cases = 0usize;
        let mut c_cases = 0usize;
        for &odd in bits {
            if delta > 0 && odd {
                // A: (x0,x1) -> (2*x1, x1-x0)
                super::super::mod_sub_qq_fast(b, x1, x0, p);
                super::super::mod_add_qq_fast(b, x0, x1, p);
                super::super::mod_double_inplace_fast(b, x0, p);
                delta = 1 - delta;
                a_cases += 1;
            } else if odd {
                // B: (x0,x1) -> (2*x0, x0+x1)
                super::super::mod_add_qq_fast(b, x1, x0, p);
                super::super::mod_double_inplace_fast(b, x0, p);
                delta = 1 + delta;
                b_cases += 1;
            } else {
                // C: (x0,x1) -> (2*x0, x1)
                super::super::mod_double_inplace_fast(b, x0, p);
                delta = 1 + delta;
                c_cases += 1;
            }
        }
        for _ in 0..bits.len() {
            super::super::mod_halve_inplace_fast(b, x0, p);
            super::super::mod_halve_inplace_fast(b, x1, p);
        }
        (a_cases, b_cases, c_cases)
    }

    #[test]
    fn fixed_branch_numerator_window_matches_scaled_by_matrix() {
        const W: usize = 16;
        let p_mod = SECP256K1_P;
        let delta = 1i64;
        let f_low = 1i128;
        let g_low = 3i128;
        let bits = branch_bits_for_lowword_window(W, delta, f_low, g_low);
        let (_, _, _, pmat) = jump_matrix_direct_lowword(W, W, delta, f_low, g_low);
        let mut b = super::super::B::new();
        let x0 = b.alloc_qubits(256);
        let x1 = b.alloc_qubits(256);
        let cases = emit_fixed_branch_numerator_scaled_window_for_test(&mut b, delta, &bits, &x0, &x1, p_mod);
        let ccx = count_ccx(&b.ops);
        let peak = b.peak_qubits;
        let num_qubits = b.next_qubit as usize;
        let num_bits = b.next_bit as usize;
        let ops = b.ops;
        let inv2w = inv_pow2_mod_p_for_test(W, p_mod);
        let row_expected = |a: U256, c: U256, c0: i128, c1: i128| -> U256 {
            let t0 = mulm(signed_i128_mod_p(c0, p_mod), a, p_mod);
            let t1 = mulm(signed_i128_mod_p(c1, p_mod), c, p_mod);
            mulm(addm(t0, t1, p_mod), inv2w, p_mod)
        };
        let mut sx = Sampler::new(b"by-branch-num-x0-v1", p_mod);
        let mut sy = Sampler::new(b"by-branch-num-x1-v1", p_mod);
        for _ in 0..32 {
            let a = sx.next();
            let c = sy.next();
            let exp0 = row_expected(a, c, pmat.m00, pmat.m01);
            let exp1 = row_expected(a, c, pmat.m10, pmat.m11);
            let mut hasher = sha3::Shake128::default();
            hasher.update(b"by-branch-num-sim-v1");
            let mut xof = hasher.finalize_xof();
            let mut sim = crate::sim::Simulator::new(num_qubits, num_bits, &mut xof);
            set_slice_u512_by(&mut sim, &x0, u256_to_u512_for_by_tests(a));
            set_slice_u512_by(&mut sim, &x1, u256_to_u512_for_by_tests(c));
            sim.apply(&ops);
            assert_eq!(get_slice_u512_by(&sim, &x0), u256_to_u512_for_by_tests(exp0), "row0 mismatch");
            assert_eq!(get_slice_u512_by(&sim, &x1), u256_to_u512_for_by_tests(exp1), "row1 mismatch");
        }
        eprintln!(
            "BY fixed branch-numerator scaled window: ccx={ccx}, peak={peak}q, cases={cases:?}, matrix={pmat:?}"
        );
        assert!(peak < 1_200, "branch-numerator fixed window lost scratch advantage");
        assert!(ccx < 35_000, "branch-numerator fixed window not cheaper than Hermite sample");
    }

    #[test]
    fn quantum_controlled_branch_numerator_replay_is_too_expensive_naively() {
        // The fixed branch-numerator window is the right arithmetic lower
        // bound, but if every BY case is selected by live quantum controls with
        // today's cmod_add/cmod_sub primitives, the control tax reverts to the
        // old microstep bottleneck. Keep this as a guardrail: the selected
        // implementation needs a structural trick, not generic controlled
        // modular additions for all possible cases.
        const W: usize = 16;
        let p = SECP256K1_P;
        let mut b = super::super::B::new();
        let a_ctrl = b.alloc_qubit();
        let b_ctrl = b.alloc_qubit();
        let r = b.alloc_qubits(256);
        let s = b.alloc_qubits(256);
        let start = b.ops.len();
        for _ in 0..W {
            emit_tagged_modular_microstep_for_cost(&mut b, &r, &s, a_ctrl, b_ctrl, p);
        }
        for _ in 0..W {
            super::super::mod_halve_inplace_fast(&mut b, &r, p);
            super::super::mod_halve_inplace_fast(&mut b, &s, p);
        }
        let ccx = count_ccx(&b.ops[start..]);
        let approx35 = ccx as f64 * 35.0;
        eprintln!(
            "BY naive quantum-controlled branch-numerator replay: window_ccx={ccx}, approx35≈{approx35:.0}, peak={}q",
            b.peak_qubits
        );
        assert!(approx35 > 2_500_000.0, "naive controlled branch replay unexpectedly SOTA-shaped");
    }

    fn emit_cmod_neg_for_test(
        b: &mut super::super::B,
        v: &[super::super::QubitId],
        ctrl: super::super::QubitId,
        p: U256,
    ) {
        // ctrl ? (p-v) : v. Like mod_neg_inplace_fast, this maps v=0 to the
        // noncanonical representative p; good enough for this structural probe
        // and the same exceptional shape as existing fast negation.
        for &q in v {
            b.cx(ctrl, q);
        }
        super::super::cadd_nbit_const_fast(b, v, p.wrapping_add(U256::from(1u64)), ctrl);
    }

    fn emit_scaled_by_controlled_microstep_inverse_negr_for_test(
        b: &mut super::super::B,
        u_neg_r: &[super::super::QubitId],
        s: &[super::super::QubitId],
        odd_ctrl: super::super::QubitId,
        a_ctrl: super::super::QubitId,
        p: U256,
    ) {
        // Inverse scaled step in the sign-flipped representation u=-r:
        //   C: (u,s) -> (u, 2s)
        //   B: (u,s) -> (u, 2s+u)
        //   A: (u,s) -> (u+2s, -u)
        // This replaces cmod_sub by cmod_add and has the same cost shape as
        // the forward scaled microstep.
        super::super::mod_double_inplace_fast(b, s, p);
        super::super::cmod_add_qq(b, s, u_neg_r, odd_ctrl, p);
        for i in 0..u_neg_r.len() {
            super::super::cswap(b, a_ctrl, u_neg_r[i], s[i]);
        }
        emit_cmod_neg_for_test(b, s, a_ctrl, p);
    }

    fn emit_scaled_by_controlled_microstep_inverse_for_test(
        b: &mut super::super::B,
        r: &[super::super::QubitId],
        s: &[super::super::QubitId],
        odd_ctrl: super::super::QubitId,
        a_ctrl: super::super::QubitId,
        p: U256,
    ) {
        // Inverse of the scaled step below:
        //   undo halve(s), undo s += r under odd, undo controlled neg(s), undo A swap.
        super::super::mod_double_inplace_fast(b, s, p);
        super::super::cmod_sub_qq(b, s, r, odd_ctrl, p);
        emit_cmod_neg_for_test(b, s, a_ctrl, p);
        for i in 0..r.len() {
            super::super::cswap(b, a_ctrl, r[i], s[i]);
        }
    }

    fn emit_twos_complement_cneg_for_test(
        b: &mut super::super::B,
        v: &[super::super::QubitId],
        ctrl: super::super::QubitId,
    ) {
        for &q in v {
            b.cx(ctrl, q);
        }
        super::super::cadd_nbit_const_fast(b, v, U256::from(1u64), ctrl);
    }

    fn emit_twos_complement_cneg_exact_for_test(
        b: &mut super::super::B,
        v: &[super::super::QubitId],
        ctrl: super::super::QubitId,
    ) {
        for &q in v {
            b.cx(ctrl, q);
        }
        super::super::cadd_nbit_const(b, v, U256::from(1u64), ctrl);
    }

    fn emit_logical_shift_right_even_for_test(b: &mut super::super::B, v: &[super::super::QubitId]) {
        // Reversible rotation that equals logical `/2` on the promised even
        // subspace because the old low bit is zero and rotates into the top.
        for i in 0..v.len() - 1 {
            b.swap(v[i], v[i + 1]);
        }
    }

    fn emit_arithmetic_shift_right_even_for_test(b: &mut super::super::B, v: &[super::super::QubitId]) {
        // Reversible arithmetic `/2` on sign-extended even signed values. The
        // logical rotate puts the promised zero low bit in the top position;
        // duplicate the old sign bit (now at n-2) into that known-zero slot.
        emit_logical_shift_right_even_for_test(b, v);
        if v.len() >= 2 {
            b.cx(v[v.len() - 2], v[v.len() - 1]);
        }
    }

    fn emit_arithmetic_shift_left_even_inverse_for_test(b: &mut super::super::B, v: &[super::super::QubitId]) {
        if v.len() >= 2 {
            b.cx(v[v.len() - 2], v[v.len() - 1]);
        }
        emit_logical_shift_left_even_inverse_for_test(b, v);
    }

    fn emit_delta_positive_into_for_test(
        b: &mut super::super::B,
        delta: &[super::super::QubitId],
        flag: super::super::QubitId,
    ) {
        let nz = b.alloc_qubit();
        super::super::cmp_neq_zero_into(b, delta, nz);
        let sign = delta[delta.len() - 1];
        b.x(sign);
        b.ccx(nz, sign, flag);
        b.x(sign);
        super::super::cmp_neq_zero_into(b, delta, nz);
        b.free(nz);
    }

    fn emit_logical_shift_left_even_inverse_for_test(b: &mut super::super::B, v: &[super::super::QubitId]) {
        for i in (0..v.len() - 1).rev() {
            b.swap(v[i], v[i + 1]);
        }
    }

    fn emit_2adic_by_branch_step_reverse_for_test(
        b: &mut super::super::B,
        f: &[super::super::QubitId],
        g: &[super::super::QubitId],
        delta: &[super::super::QubitId],
        odd_hist: super::super::QubitId,
        a_hist: super::super::QubitId,
    ) {
        // Reverse of emit_2adic_by_branch_step_for_test, then clear the two
        // branch-history bits from the restored pre-step denominator state.
        super::super::sub_nbit_const_fast(b, delta, U256::from(1u64));
        emit_twos_complement_cneg_for_test(b, delta, a_hist);
        emit_logical_shift_left_even_inverse_for_test(b, g);
        super::super::cucc_sub_ctrl(b, f, g, odd_hist);
        emit_twos_complement_cneg_for_test(b, g, a_hist);
        for i in 0..f.len() {
            super::super::cswap(b, a_hist, f[i], g[i]);
        }

        let positive = b.alloc_qubit();
        emit_delta_positive_into_for_test(b, delta, positive);
        b.ccx(odd_hist, positive, a_hist);
        emit_delta_positive_into_for_test(b, delta, positive);
        b.free(positive);
        b.cx(g[0], odd_hist);
    }

    fn emit_2adic_denominator_step_with_controls_for_test(
        b: &mut super::super::B,
        f: &[super::super::QubitId],
        g: &[super::super::QubitId],
        odd: super::super::QubitId,
        a: super::super::QubitId,
    ) {
        for i in 0..f.len() {
            super::super::cswap(b, a, f[i], g[i]);
        }
        emit_twos_complement_cneg_for_test(b, g, a);
        super::super::cucc_add_ctrl(b, f, g, odd);
        emit_logical_shift_right_even_for_test(b, g);
    }

    fn emit_2adic_by_branch_step_for_test(
        b: &mut super::super::B,
        f: &[super::super::QubitId],
        g: &[super::super::QubitId],
        delta: &[super::super::QubitId],
        odd_out: super::super::QubitId,
        a_out: super::super::QubitId,
    ) {
        b.cx(g[0], odd_out);
        let positive = b.alloc_qubit();
        emit_delta_positive_into_for_test(b, delta, positive);
        b.ccx(odd_out, positive, a_out);
        emit_delta_positive_into_for_test(b, delta, positive);
        b.free(positive);

        for i in 0..f.len() {
            super::super::cswap(b, a_out, f[i], g[i]);
        }
        emit_twos_complement_cneg_for_test(b, g, a_out);
        super::super::cucc_add_ctrl(b, f, g, odd_out);
        emit_logical_shift_right_even_for_test(b, g);

        emit_twos_complement_cneg_for_test(b, delta, a_out);
        super::super::add_nbit_const_fast(b, delta, U256::from(1u64));
    }

    fn emit_signed_by_branch_step_for_test(
        b: &mut super::super::B,
        f: &[super::super::QubitId],
        g: &[super::super::QubitId],
        delta: &[super::super::QubitId],
        odd_out: super::super::QubitId,
        a_out: super::super::QubitId,
    ) {
        b.cx(g[0], odd_out);
        let positive = b.alloc_qubit();
        emit_delta_positive_into_for_test(b, delta, positive);
        b.ccx(odd_out, positive, a_out);
        emit_delta_positive_into_for_test(b, delta, positive);
        b.free(positive);

        for i in 0..f.len() {
            super::super::cswap(b, a_out, f[i], g[i]);
        }
        emit_twos_complement_cneg_for_test(b, g, a_out);
        super::super::cucc_add_ctrl(b, f, g, odd_out);
        emit_arithmetic_shift_right_even_for_test(b, g);

        emit_twos_complement_cneg_for_test(b, delta, a_out);
        super::super::add_nbit_const_fast(b, delta, U256::from(1u64));
    }

    fn emit_signed_by_branch_step_reverse_for_test(
        b: &mut super::super::B,
        f: &[super::super::QubitId],
        g: &[super::super::QubitId],
        delta: &[super::super::QubitId],
        odd_hist: super::super::QubitId,
        a_hist: super::super::QubitId,
    ) {
        super::super::sub_nbit_const_fast(b, delta, U256::from(1u64));
        emit_twos_complement_cneg_for_test(b, delta, a_hist);
        emit_arithmetic_shift_left_even_inverse_for_test(b, g);
        super::super::cucc_sub_ctrl(b, f, g, odd_hist);
        emit_twos_complement_cneg_for_test(b, g, a_hist);
        for i in 0..f.len() {
            super::super::cswap(b, a_hist, f[i], g[i]);
        }

        let positive = b.alloc_qubit();
        emit_delta_positive_into_for_test(b, delta, positive);
        b.ccx(odd_hist, positive, a_hist);
        emit_delta_positive_into_for_test(b, delta, positive);
        b.free(positive);
        b.cx(g[0], odd_hist);
    }

    fn emit_cmod_neg_exact_for_test(
        b: &mut super::super::B,
        v: &[super::super::QubitId],
        ctrl: super::super::QubitId,
        p: U256,
    ) {
        for &q in v {
            b.cx(ctrl, q);
        }
        super::super::cadd_nbit_const(b, v, p.wrapping_add(U256::from(1u64)), ctrl);
    }

    fn cmod_add_qq_exact_for_test(
        b: &mut super::super::B,
        acc: &[super::super::QubitId],
        a: &[super::super::QubitId],
        ctrl: super::super::QubitId,
        p: U256,
    ) {
        let f = b.alloc_qubits(acc.len());
        for i in 0..acc.len() {
            b.ccx(ctrl, a[i], f[i]);
        }
        super::super::mod_add_qq(b, acc, &f, p);
        for i in 0..acc.len() {
            b.ccx(ctrl, a[i], f[i]);
        }
        b.free_vec(&f);
    }

    fn emit_scaled_by_controlled_microstep_exact_for_test(
        b: &mut super::super::B,
        r: &[super::super::QubitId],
        s: &[super::super::QubitId],
        odd_ctrl: super::super::QubitId,
        a_ctrl: super::super::QubitId,
        p: U256,
    ) {
        for i in 0..r.len() {
            super::super::cswap(b, a_ctrl, r[i], s[i]);
        }
        emit_cmod_neg_exact_for_test(b, s, a_ctrl, p);
        cmod_add_qq_exact_for_test(b, s, r, odd_ctrl, p);
        super::super::mod_halve_inplace(b, s, p);
    }

    fn emit_scaled_by_controlled_microstep_exact_cneg_for_test(
        b: &mut super::super::B,
        r: &[super::super::QubitId],
        s: &[super::super::QubitId],
        odd_ctrl: super::super::QubitId,
        a_ctrl: super::super::QubitId,
        p: U256,
    ) {
        for i in 0..r.len() {
            super::super::cswap(b, a_ctrl, r[i], s[i]);
        }
        // The A-controlled negation is the only MBU operation directly
        // controlled by A in the fast microstep. Making just this exact may be
        // enough to allow window-local A clearing while preserving most of the
        // fast cmod_add+halve savings.
        emit_cmod_neg_exact_for_test(b, s, a_ctrl, p);
        super::super::cmod_add_qq(b, s, r, odd_ctrl, p);
        super::super::mod_halve_inplace_fast(b, s, p);
    }

    fn mod_add_qq_fast_keep_reduction_flag_for_test(
        b: &mut super::super::B,
        acc: &[super::super::QubitId],
        a: &[super::super::QubitId],
        p: U256,
        reduction_flag: super::super::QubitId,
    ) {
        // Same forward modular addition as mod_add_qq_fast, but keep the
        // reduction flag live instead of paying cmp_lt_into_fast to uncompute
        // it immediately. This is not a complete primitive; it quantifies the
        // exact flag-cleaning obstacle for a fused modular average.
        let n = acc.len();
        assert_eq!(n, a.len());
        let (acc_ext, acc_ovf) = super::super::ext_reg(b, acc);
        let (a_ext, a_ovf) = super::super::ext_reg(b, a);
        super::super::add_nbit_qq_fast(b, &a_ext, &acc_ext);
        let c = U256::MAX.wrapping_sub(p).wrapping_add(U256::from(1u64));
        super::super::add_nbit_const_fast(b, &acc_ext, c);
        b.cx(acc_ovf, reduction_flag);
        b.x(reduction_flag);
        super::super::csub_nbit_const_fast(b, &acc_ext, c, reduction_flag);
        b.x(reduction_flag);
        b.cx(reduction_flag, acc_ovf);
        super::super::unext_reg(b, a_ovf);
        super::super::unext_reg(b, acc_ovf);
        let _ = (acc_ext, a_ext);
    }

    fn cmod_add_qq_keep_reduction_flag_for_test(
        b: &mut super::super::B,
        acc: &[super::super::QubitId],
        a: &[super::super::QubitId],
        ctrl: super::super::QubitId,
        p: U256,
        reduction_flag: super::super::QubitId,
    ) {
        let f = b.alloc_qubits(acc.len());
        for i in 0..acc.len() {
            b.ccx(ctrl, a[i], f[i]);
        }
        mod_add_qq_fast_keep_reduction_flag_for_test(b, acc, &f, p, reduction_flag);
        for i in 0..acc.len() {
            let m = b.alloc_bit();
            b.hmr(f[i], m);
            b.cz_if(ctrl, a[i], m);
        }
        b.free_vec(&f);
    }

    fn emit_scaled_by_controlled_microstep_live_addflag_for_test(
        b: &mut super::super::B,
        r: &[super::super::QubitId],
        s: &[super::super::QubitId],
        odd_ctrl: super::super::QubitId,
        a_ctrl: super::super::QubitId,
        reduction_flag: super::super::QubitId,
        p: U256,
    ) {
        for i in 0..r.len() {
            super::super::cswap(b, a_ctrl, r[i], s[i]);
        }
        emit_cmod_neg_for_test(b, s, a_ctrl, p);
        cmod_add_qq_keep_reduction_flag_for_test(b, s, r, odd_ctrl, p, reduction_flag);
        super::super::mod_halve_inplace_fast(b, s, p);
    }

    fn emit_signed_controlled_add_for_test(
        b: &mut super::super::B,
        acc: &[super::super::QubitId],
        a: &[super::super::QubitId],
        ctrl: super::super::QubitId,
    ) {
        let f = b.alloc_qubits(acc.len());
        for i in 0..acc.len() {
            b.ccx(ctrl, a[i], f[i]);
        }
        super::super::add_nbit_qq_fast(b, &f, acc);
        for i in 0..acc.len() {
            let m = b.alloc_bit();
            b.hmr(f[i], m);
            b.cz_if(ctrl, a[i], m);
        }
        b.free_vec(&f);
    }

    fn emit_signed_controlled_sub_for_test(
        b: &mut super::super::B,
        acc: &[super::super::QubitId],
        a: &[super::super::QubitId],
        ctrl: super::super::QubitId,
    ) {
        let f = b.alloc_qubits(acc.len());
        for i in 0..acc.len() {
            b.ccx(ctrl, a[i], f[i]);
        }
        super::super::sub_nbit_qq_fast(b, &f, acc);
        for i in 0..acc.len() {
            let m = b.alloc_bit();
            b.hmr(f[i], m);
            b.cz_if(ctrl, a[i], m);
        }
        b.free_vec(&f);
    }

    fn emit_signed_controlled_add_negcopy_for_test(
        b: &mut super::super::B,
        acc: &[super::super::QubitId],
        a: &[super::super::QubitId],
        ctrl: super::super::QubitId,
    ) {
        let f = b.alloc_qubits(acc.len());
        for i in 0..acc.len() {
            b.ccx(ctrl, a[i], f[i]);
        }
        super::super::add_nbit_qq_fast(b, &f, acc);
        for i in 0..acc.len() {
            let m = b.alloc_bit();
            b.hmr(f[i], m);
            b.cz_if(ctrl, a[i], m);
            b.neg_if(m);
        }
        b.free_vec(&f);
    }

    fn emit_signed_controlled_sub_negcopy_for_test(
        b: &mut super::super::B,
        acc: &[super::super::QubitId],
        a: &[super::super::QubitId],
        ctrl: super::super::QubitId,
    ) {
        let f = b.alloc_qubits(acc.len());
        for i in 0..acc.len() {
            b.ccx(ctrl, a[i], f[i]);
        }
        super::super::sub_nbit_qq_fast(b, &f, acc);
        for i in 0..acc.len() {
            let m = b.alloc_bit();
            b.hmr(f[i], m);
            b.cz_if(ctrl, a[i], m);
            b.neg_if(m);
        }
        b.free_vec(&f);
    }

    fn cuccaro_add_fast_negphase_for_test(
        b: &mut super::super::B,
        a: &[super::super::QubitId],
        acc: &[super::super::QubitId],
        c_in: super::super::QubitId,
    ) {
        let n = a.len();
        assert_eq!(n, acc.len());
        if n == 0 { return; }
        if n == 1 {
            b.cx(c_in, acc[0]);
            b.cx(a[0], acc[0]);
            return;
        }
        let carries = b.alloc_qubits(n - 1);
        b.cx(a[0], acc[0]);
        b.cx(a[0], c_in);
        b.ccx(c_in, acc[0], carries[0]);
        b.cx(carries[0], a[0]);
        for i in 1..n - 1 {
            b.cx(a[i], acc[i]);
            b.cx(a[i], a[i - 1]);
            b.ccx(a[i - 1], acc[i], carries[i]);
            b.cx(carries[i], a[i]);
        }
        b.cx(a[n - 2], acc[n - 1]);
        b.cx(a[n - 1], acc[n - 1]);
        for i in (1..n - 1).rev() {
            b.cx(carries[i], a[i]);
            let m = b.alloc_bit();
            b.hmr(carries[i], m);
            b.cz_if(a[i - 1], acc[i], m);
            b.neg_if(m);
            b.cx(a[i], a[i - 1]);
            b.cx(a[i - 1], acc[i]);
        }
        b.cx(carries[0], a[0]);
        let m0 = b.alloc_bit();
        b.hmr(carries[0], m0);
        b.cz_if(c_in, acc[0], m0);
        b.neg_if(m0);
        b.cx(a[0], c_in);
        b.cx(c_in, acc[0]);
        b.free_vec(&carries);
    }

    fn cuccaro_sub_fast_negphase_for_test(
        b: &mut super::super::B,
        a: &[super::super::QubitId],
        acc: &[super::super::QubitId],
        c_in: super::super::QubitId,
    ) {
        let n = a.len();
        assert_eq!(n, acc.len());
        if n == 0 { return; }
        if n == 1 {
            b.cx(a[0], acc[0]);
            b.cx(c_in, acc[0]);
            return;
        }
        let carries = b.alloc_qubits(n - 1);
        b.cx(c_in, acc[0]);
        b.cx(a[0], c_in);
        b.ccx(c_in, acc[0], carries[0]);
        b.cx(carries[0], a[0]);
        for i in 1..n - 1 {
            b.cx(a[i - 1], acc[i]);
            b.cx(a[i], a[i - 1]);
            b.ccx(a[i - 1], acc[i], carries[i]);
            b.cx(carries[i], a[i]);
        }
        b.cx(a[n - 1], acc[n - 1]);
        b.cx(a[n - 2], acc[n - 1]);
        for i in (1..n - 1).rev() {
            b.cx(carries[i], a[i]);
            let m = b.alloc_bit();
            b.hmr(carries[i], m);
            b.cz_if(a[i - 1], acc[i], m);
            b.neg_if(m);
            b.cx(a[i], a[i - 1]);
            b.cx(a[i], acc[i]);
        }
        b.cx(carries[0], a[0]);
        let m0 = b.alloc_bit();
        b.hmr(carries[0], m0);
        b.cz_if(c_in, acc[0], m0);
        b.neg_if(m0);
        b.cx(a[0], c_in);
        b.cx(a[0], acc[0]);
        b.free_vec(&carries);
    }

    fn add_nbit_qq_fast_negphase_for_test(
        b: &mut super::super::B,
        a: &[super::super::QubitId],
        acc: &[super::super::QubitId],
    ) {
        let c_in = b.alloc_qubit();
        cuccaro_add_fast_negphase_for_test(b, a, acc, c_in);
        b.free(c_in);
    }

    fn sub_nbit_qq_fast_negphase_for_test(
        b: &mut super::super::B,
        a: &[super::super::QubitId],
        acc: &[super::super::QubitId],
    ) {
        let c_in = b.alloc_qubit();
        cuccaro_sub_fast_negphase_for_test(b, a, acc, c_in);
        b.free(c_in);
    }

    fn emit_signed_controlled_add_negcarry_for_test(
        b: &mut super::super::B,
        acc: &[super::super::QubitId],
        a: &[super::super::QubitId],
        ctrl: super::super::QubitId,
    ) {
        let f = b.alloc_qubits(acc.len());
        for i in 0..acc.len() { b.ccx(ctrl, a[i], f[i]); }
        add_nbit_qq_fast_negphase_for_test(b, &f, acc);
        for i in 0..acc.len() {
            let m = b.alloc_bit();
            b.hmr(f[i], m);
            b.cz_if(ctrl, a[i], m);
        }
        b.free_vec(&f);
    }

    fn emit_signed_controlled_sub_negcarry_for_test(
        b: &mut super::super::B,
        acc: &[super::super::QubitId],
        a: &[super::super::QubitId],
        ctrl: super::super::QubitId,
    ) {
        let f = b.alloc_qubits(acc.len());
        for i in 0..acc.len() { b.ccx(ctrl, a[i], f[i]); }
        sub_nbit_qq_fast_negphase_for_test(b, &f, acc);
        for i in 0..acc.len() {
            let m = b.alloc_bit();
            b.hmr(f[i], m);
            b.cz_if(ctrl, a[i], m);
        }
        b.free_vec(&f);
    }

    fn emit_signed_controlled_add_negboth_for_test(
        b: &mut super::super::B,
        acc: &[super::super::QubitId],
        a: &[super::super::QubitId],
        ctrl: super::super::QubitId,
    ) {
        let f = b.alloc_qubits(acc.len());
        for i in 0..acc.len() { b.ccx(ctrl, a[i], f[i]); }
        add_nbit_qq_fast_negphase_for_test(b, &f, acc);
        for i in 0..acc.len() {
            let m = b.alloc_bit();
            b.hmr(f[i], m);
            b.cz_if(ctrl, a[i], m);
            b.neg_if(m);
        }
        b.free_vec(&f);
    }

    fn emit_signed_controlled_sub_negboth_for_test(
        b: &mut super::super::B,
        acc: &[super::super::QubitId],
        a: &[super::super::QubitId],
        ctrl: super::super::QubitId,
    ) {
        let f = b.alloc_qubits(acc.len());
        for i in 0..acc.len() { b.ccx(ctrl, a[i], f[i]); }
        sub_nbit_qq_fast_negphase_for_test(b, &f, acc);
        for i in 0..acc.len() {
            let m = b.alloc_bit();
            b.hmr(f[i], m);
            b.cz_if(ctrl, a[i], m);
            b.neg_if(m);
        }
        b.free_vec(&f);
    }

    fn emit_signed_controlled_add_exact_for_test(
        b: &mut super::super::B,
        acc: &[super::super::QubitId],
        a: &[super::super::QubitId],
        ctrl: super::super::QubitId,
    ) {
        let f = b.alloc_qubits(acc.len());
        for i in 0..acc.len() {
            b.ccx(ctrl, a[i], f[i]);
        }
        super::super::add_nbit_qq(b, &f, acc);
        for i in 0..acc.len() {
            b.ccx(ctrl, a[i], f[i]);
        }
        b.free_vec(&f);
    }

    fn emit_signed_controlled_sub_exact_for_test(
        b: &mut super::super::B,
        acc: &[super::super::QubitId],
        a: &[super::super::QubitId],
        ctrl: super::super::QubitId,
    ) {
        let f = b.alloc_qubits(acc.len());
        for i in 0..acc.len() {
            b.ccx(ctrl, a[i], f[i]);
        }
        super::super::sub_nbit_qq(b, &f, acc);
        for i in 0..acc.len() {
            b.ccx(ctrl, a[i], f[i]);
        }
        b.free_vec(&f);
    }

    fn emit_signed_redundant_halve_live_parity_for_test(
        b: &mut super::super::B,
        v: &[super::super::QubitId],
        parity_hist: super::super::QubitId,
        p: U256,
    ) {
        b.cx(v[0], parity_hist);
        super::super::cadd_nbit_const_fast(b, v, p, parity_hist);
        emit_arithmetic_shift_right_even_for_test(b, v);
    }

    fn emit_signed_redundant_halve_centered_live_parity_for_test(
        b: &mut super::super::B,
        v: &[super::super::QubitId],
        parity_hist: super::super::QubitId,
        p: U256,
    ) {
        // Centered no-reduction halve: if signed pre-halve T is odd, add p for
        // T<0 and subtract p for T>=0, then arithmetic shift. The copied old
        // sign is cleared after the shift using old_sign = new_sign XOR parity
        // on the promised centered range.
        let sign_hist = b.alloc_qubit();
        let add_ctrl = b.alloc_qubit();
        let sub_ctrl = b.alloc_qubit();
        b.cx(v[0], parity_hist);
        b.cx(v[v.len() - 1], sign_hist);
        b.ccx(parity_hist, sign_hist, add_ctrl);
        b.x(sign_hist);
        b.ccx(parity_hist, sign_hist, sub_ctrl);
        b.x(sign_hist);
        super::super::cadd_nbit_const_fast(b, v, p, add_ctrl);
        super::super::csub_nbit_const_fast(b, v, p, sub_ctrl);
        b.x(sign_hist);
        b.ccx(parity_hist, sign_hist, sub_ctrl);
        b.x(sign_hist);
        b.ccx(parity_hist, sign_hist, add_ctrl);
        b.free(sub_ctrl);
        b.free(add_ctrl);
        emit_arithmetic_shift_right_even_for_test(b, v);
        b.cx(v[v.len() - 1], sign_hist);
        b.cx(parity_hist, sign_hist);
        b.free(sign_hist);
    }

    fn emit_signed_redundant_halve_centered_live_parity_exact_const_for_test(
        b: &mut super::super::B,
        v: &[super::super::QubitId],
        parity_hist: super::super::QubitId,
        p: U256,
    ) {
        let sign_hist = b.alloc_qubit();
        let add_ctrl = b.alloc_qubit();
        let sub_ctrl = b.alloc_qubit();
        b.cx(v[0], parity_hist);
        b.cx(v[v.len() - 1], sign_hist);
        b.ccx(parity_hist, sign_hist, add_ctrl);
        b.x(sign_hist);
        b.ccx(parity_hist, sign_hist, sub_ctrl);
        b.x(sign_hist);
        super::super::cadd_nbit_const(b, v, p, add_ctrl);
        super::super::csub_nbit_const(b, v, p, sub_ctrl);
        b.x(sign_hist);
        b.ccx(parity_hist, sign_hist, sub_ctrl);
        b.x(sign_hist);
        b.ccx(parity_hist, sign_hist, add_ctrl);
        b.free(sub_ctrl);
        b.free(add_ctrl);
        emit_arithmetic_shift_right_even_for_test(b, v);
        b.cx(v[v.len() - 1], sign_hist);
        b.cx(parity_hist, sign_hist);
        b.free(sign_hist);
    }

    fn emit_scaled_by_redundant_signed_microstep_live_parity_for_test(
        b: &mut super::super::B,
        r: &[super::super::QubitId],
        s: &[super::super::QubitId],
        odd_ctrl: super::super::QubitId,
        a_ctrl: super::super::QubitId,
        parity_hist: super::super::QubitId,
        p: U256,
    ) {
        for i in 0..r.len() {
            super::super::cswap(b, a_ctrl, r[i], s[i]);
        }
        emit_twos_complement_cneg_for_test(b, s, a_ctrl);
        emit_signed_controlled_add_for_test(b, s, r, odd_ctrl);
        emit_signed_redundant_halve_live_parity_for_test(b, s, parity_hist, p);
    }

    fn emit_scaled_by_centered_signed_microstep_live_parity_for_test(
        b: &mut super::super::B,
        r: &[super::super::QubitId],
        s: &[super::super::QubitId],
        odd_ctrl: super::super::QubitId,
        a_ctrl: super::super::QubitId,
        parity_hist: super::super::QubitId,
        p: U256,
    ) {
        for i in 0..r.len() {
            super::super::cswap(b, a_ctrl, r[i], s[i]);
        }
        emit_twos_complement_cneg_for_test(b, s, a_ctrl);
        emit_signed_controlled_add_for_test(b, s, r, odd_ctrl);
        emit_signed_redundant_halve_centered_live_parity_for_test(b, s, parity_hist, p);
    }

    fn emit_signed_redundant_unhalve_centered_with_parity_for_test(
        b: &mut super::super::B,
        v: &[super::super::QubitId],
        parity_hist: super::super::QubitId,
        p: U256,
    ) {
        // Inverse of the centered halve with a stored parity bit: arithmetic
        // shift left, then if parity was set undo the sign-conditioned ±p.
        // The correction controls must be uncomputed from the pre-correction
        // sign.  When parity=1 the correction flips the sign, so using the
        // post-correction sign leaves dirty controls and R-phase garbage.
        emit_arithmetic_shift_left_even_inverse_for_test(b, v);
        let sign_hist = b.alloc_qubit();
        let add_ctrl = b.alloc_qubit();
        let sub_ctrl = b.alloc_qubit();
        let sign = v[v.len() - 1];
        b.cx(sign, sign_hist);
        b.ccx(parity_hist, sign_hist, add_ctrl);
        b.x(sign_hist);
        b.ccx(parity_hist, sign_hist, sub_ctrl);
        b.x(sign_hist);
        super::super::cadd_nbit_const_fast(b, v, p, add_ctrl);
        super::super::csub_nbit_const_fast(b, v, p, sub_ctrl);
        b.x(sign_hist);
        b.ccx(parity_hist, sign_hist, sub_ctrl);
        b.x(sign_hist);
        b.ccx(parity_hist, sign_hist, add_ctrl);
        b.free(sub_ctrl);
        b.free(add_ctrl);
        b.cx(sign, sign_hist);
        b.cx(parity_hist, sign_hist);
        b.free(sign_hist);
    }

    fn emit_signed_redundant_unhalve_centered_with_parity_exact_const_for_test(
        b: &mut super::super::B,
        v: &[super::super::QubitId],
        parity_hist: super::super::QubitId,
        p: U256,
    ) {
        emit_arithmetic_shift_left_even_inverse_for_test(b, v);
        let sign_hist = b.alloc_qubit();
        let add_ctrl = b.alloc_qubit();
        let sub_ctrl = b.alloc_qubit();
        let sign = v[v.len() - 1];
        b.cx(sign, sign_hist);
        b.ccx(parity_hist, sign_hist, add_ctrl);
        b.x(sign_hist);
        b.ccx(parity_hist, sign_hist, sub_ctrl);
        b.x(sign_hist);
        super::super::cadd_nbit_const(b, v, p, add_ctrl);
        super::super::csub_nbit_const(b, v, p, sub_ctrl);
        b.x(sign_hist);
        b.ccx(parity_hist, sign_hist, sub_ctrl);
        b.x(sign_hist);
        b.ccx(parity_hist, sign_hist, add_ctrl);
        b.free(sub_ctrl);
        b.free(add_ctrl);
        b.cx(sign, sign_hist);
        b.cx(parity_hist, sign_hist);
        b.free(sign_hist);
    }

    fn emit_scaled_by_centered_signed_microstep_inverse_live_parity_for_test(
        b: &mut super::super::B,
        r: &[super::super::QubitId],
        s: &[super::super::QubitId],
        odd_ctrl: super::super::QubitId,
        a_ctrl: super::super::QubitId,
        parity_hist: super::super::QubitId,
        p: U256,
    ) {
        emit_signed_redundant_unhalve_centered_with_parity_for_test(b, s, parity_hist, p);
        emit_signed_controlled_sub_for_test(b, s, r, odd_ctrl);
        emit_twos_complement_cneg_for_test(b, s, a_ctrl);
        for i in 0..r.len() {
            super::super::cswap(b, a_ctrl, r[i], s[i]);
        }
    }

    fn emit_scaled_by_centered_signed_microstep_live_parity_variant_for_test(
        b: &mut super::super::B,
        r: &[super::super::QubitId],
        s: &[super::super::QubitId],
        odd_ctrl: super::super::QubitId,
        a_ctrl: super::super::QubitId,
        parity_hist: super::super::QubitId,
        p: U256,
        exact_signed_add: bool,
        exact_parity_const: bool,
    ) {
        for i in 0..r.len() {
            super::super::cswap(b, a_ctrl, r[i], s[i]);
        }
        if exact_parity_const {
            emit_twos_complement_cneg_exact_for_test(b, s, a_ctrl);
        } else {
            emit_twos_complement_cneg_for_test(b, s, a_ctrl);
        }
        if exact_signed_add {
            emit_signed_controlled_add_exact_for_test(b, s, r, odd_ctrl);
        } else {
            emit_signed_controlled_add_for_test(b, s, r, odd_ctrl);
        }
        if exact_parity_const {
            emit_signed_redundant_halve_centered_live_parity_exact_const_for_test(b, s, parity_hist, p);
        } else {
            emit_signed_redundant_halve_centered_live_parity_for_test(b, s, parity_hist, p);
        }
    }

    fn emit_scaled_by_centered_signed_microstep_inverse_live_parity_variant_for_test(
        b: &mut super::super::B,
        r: &[super::super::QubitId],
        s: &[super::super::QubitId],
        odd_ctrl: super::super::QubitId,
        a_ctrl: super::super::QubitId,
        parity_hist: super::super::QubitId,
        p: U256,
        exact_signed_add: bool,
        exact_parity_const: bool,
    ) {
        if exact_parity_const {
            emit_signed_redundant_unhalve_centered_with_parity_exact_const_for_test(b, s, parity_hist, p);
        } else {
            emit_signed_redundant_unhalve_centered_with_parity_for_test(b, s, parity_hist, p);
        }
        if exact_signed_add {
            emit_signed_controlled_sub_exact_for_test(b, s, r, odd_ctrl);
        } else {
            emit_signed_controlled_sub_for_test(b, s, r, odd_ctrl);
        }
        if exact_parity_const {
            emit_twos_complement_cneg_exact_for_test(b, s, a_ctrl);
        } else {
            emit_twos_complement_cneg_for_test(b, s, a_ctrl);
        }
        for i in 0..r.len() {
            super::super::cswap(b, a_ctrl, r[i], s[i]);
        }
    }

    fn emit_scaled_by_centered_signed_microstep_live_parity_negboth_exact_parity_for_test(
        b: &mut super::super::B,
        r: &[super::super::QubitId],
        s: &[super::super::QubitId],
        odd_ctrl: super::super::QubitId,
        a_ctrl: super::super::QubitId,
        parity_hist: super::super::QubitId,
        p: U256,
    ) {
        for i in 0..r.len() {
            super::super::cswap(b, a_ctrl, r[i], s[i]);
        }
        emit_twos_complement_cneg_exact_for_test(b, s, a_ctrl);
        emit_signed_controlled_add_negboth_for_test(b, s, r, odd_ctrl);
        emit_signed_redundant_halve_centered_live_parity_exact_const_for_test(b, s, parity_hist, p);
    }

    fn emit_scaled_by_centered_signed_microstep_inverse_live_parity_negboth_exact_parity_for_test(
        b: &mut super::super::B,
        r: &[super::super::QubitId],
        s: &[super::super::QubitId],
        odd_ctrl: super::super::QubitId,
        a_ctrl: super::super::QubitId,
        parity_hist: super::super::QubitId,
        p: U256,
    ) {
        emit_signed_redundant_unhalve_centered_with_parity_exact_const_for_test(b, s, parity_hist, p);
        emit_signed_controlled_sub_negboth_for_test(b, s, r, odd_ctrl);
        emit_twos_complement_cneg_exact_for_test(b, s, a_ctrl);
        for i in 0..r.len() {
            super::super::cswap(b, a_ctrl, r[i], s[i]);
        }
    }

    fn emit_scaled_by_centered_signed_microstep_live_parity_negcarry_exact_parity_for_test(
        b: &mut super::super::B,
        r: &[super::super::QubitId],
        s: &[super::super::QubitId],
        odd_ctrl: super::super::QubitId,
        a_ctrl: super::super::QubitId,
        parity_hist: super::super::QubitId,
        p: U256,
    ) {
        for i in 0..r.len() {
            super::super::cswap(b, a_ctrl, r[i], s[i]);
        }
        emit_twos_complement_cneg_exact_for_test(b, s, a_ctrl);
        emit_signed_controlled_add_negcarry_for_test(b, s, r, odd_ctrl);
        emit_signed_redundant_halve_centered_live_parity_exact_const_for_test(b, s, parity_hist, p);
    }

    fn emit_scaled_by_centered_signed_microstep_inverse_live_parity_negcarry_exact_parity_for_test(
        b: &mut super::super::B,
        r: &[super::super::QubitId],
        s: &[super::super::QubitId],
        odd_ctrl: super::super::QubitId,
        a_ctrl: super::super::QubitId,
        parity_hist: super::super::QubitId,
        p: U256,
    ) {
        emit_signed_redundant_unhalve_centered_with_parity_exact_const_for_test(b, s, parity_hist, p);
        emit_signed_controlled_sub_negcarry_for_test(b, s, r, odd_ctrl);
        emit_twos_complement_cneg_exact_for_test(b, s, a_ctrl);
        for i in 0..r.len() {
            super::super::cswap(b, a_ctrl, r[i], s[i]);
        }
    }

    fn emit_scaled_by_centered_signed_microstep_live_parity_negcopy_exact_parity_for_test(
        b: &mut super::super::B,
        r: &[super::super::QubitId],
        s: &[super::super::QubitId],
        odd_ctrl: super::super::QubitId,
        a_ctrl: super::super::QubitId,
        parity_hist: super::super::QubitId,
        p: U256,
    ) {
        for i in 0..r.len() {
            super::super::cswap(b, a_ctrl, r[i], s[i]);
        }
        emit_twos_complement_cneg_exact_for_test(b, s, a_ctrl);
        emit_signed_controlled_add_negcopy_for_test(b, s, r, odd_ctrl);
        emit_signed_redundant_halve_centered_live_parity_exact_const_for_test(b, s, parity_hist, p);
    }

    fn emit_scaled_by_centered_signed_microstep_inverse_live_parity_negcopy_exact_parity_for_test(
        b: &mut super::super::B,
        r: &[super::super::QubitId],
        s: &[super::super::QubitId],
        odd_ctrl: super::super::QubitId,
        a_ctrl: super::super::QubitId,
        parity_hist: super::super::QubitId,
        p: U256,
    ) {
        emit_signed_redundant_unhalve_centered_with_parity_exact_const_for_test(b, s, parity_hist, p);
        emit_signed_controlled_sub_negcopy_for_test(b, s, r, odd_ctrl);
        emit_twos_complement_cneg_exact_for_test(b, s, a_ctrl);
        for i in 0..r.len() {
            super::super::cswap(b, a_ctrl, r[i], s[i]);
        }
    }

    fn emit_centered_signed_clear_parity_after_inverse_for_test(
        b: &mut super::super::B,
        r: &[super::super::QubitId],
        s: &[super::super::QubitId],
        odd_ctrl: super::super::QubitId,
        parity_hist: super::super::QubitId,
    ) {
        b.cx(s[0], parity_hist);
        b.ccx(odd_ctrl, r[0], parity_hist);
    }

    fn emit_scaled_by_controlled_microstep_for_test(
        b: &mut super::super::B,
        r: &[super::super::QubitId],
        s: &[super::super::QubitId],
        odd_ctrl: super::super::QubitId,
        a_ctrl: super::super::QubitId,
        p: U256,
    ) {
        // Direct scaled BY microstep for the modular tagged pair:
        //   C: (r,s) -> (r, s/2)
        //   B: (r,s) -> (r, (s+r)/2)
        //   A: (r,s) -> (s, (s-r)/2)
        // Implement A by a controlled physical swap, then row1 <- -row1 + row0.
        // This removes the branch-numerator A-only r+=s correction and replaces
        // the per-step double+final-scale convention by one immediate halve.
        for i in 0..r.len() {
            super::super::cswap(b, a_ctrl, r[i], s[i]);
        }
        emit_cmod_neg_for_test(b, s, a_ctrl, p);
        super::super::cmod_add_qq(b, s, r, odd_ctrl, p);
        super::super::mod_halve_inplace_fast(b, s, p);
    }

    fn a_mask_from_pattern_and_delta_for_test(pattern: u16, w: usize, mut delta: i64) -> (u16, i64) {
        let mut a_mask = 0u16;
        for i in 0..w {
            let odd = ((pattern >> i) & 1) != 0;
            let a = delta > 0 && odd;
            if a {
                a_mask |= 1u16 << i;
                delta = 1 - delta;
            } else {
                delta = 1 + delta;
            }
        }
        (a_mask, delta)
    }

    fn emit_pattern_delta_decode_window_for_test(
        b: &mut super::super::B,
        pattern: &[super::super::QubitId],
        delta: &[super::super::QubitId],
        a_bits: &[super::super::QubitId],
    ) {
        assert_eq!(pattern.len(), a_bits.len());
        for i in 0..pattern.len() {
            let positive = b.alloc_qubit();
            emit_delta_positive_into_for_test(b, delta, positive);
            b.ccx(pattern[i], positive, a_bits[i]);
            emit_delta_positive_into_for_test(b, delta, positive);
            b.free(positive);
            emit_twos_complement_cneg_exact_for_test(b, delta, a_bits[i]);
            // Use exact Cuccaro add here. The decoder is run forward, then its
            // inverse after the modular replay has used A bits; measurement-
            // based fast add/sub can leave phase tied to those intervening
            // controls even when the basis data is restored.
            super::super::add_nbit_const(b, delta, U256::from(1u64));
        }
    }

    fn emit_pattern_delta_decode_window_reverse_for_test(
        b: &mut super::super::B,
        pattern: &[super::super::QubitId],
        delta: &[super::super::QubitId],
        a_bits: &[super::super::QubitId],
    ) {
        assert_eq!(pattern.len(), a_bits.len());
        for i in (0..pattern.len()).rev() {
            super::super::sub_nbit_const(b, delta, U256::from(1u64));
            emit_twos_complement_cneg_exact_for_test(b, delta, a_bits[i]);
            let positive = b.alloc_qubit();
            emit_delta_positive_into_for_test(b, delta, positive);
            b.ccx(pattern[i], positive, a_bits[i]);
            emit_delta_positive_into_for_test(b, delta, positive);
            b.free(positive);
        }
    }

    fn twos_u512_for_delta(delta: i64, width: usize) -> U512 {
        let modulus = 1i128 << width;
        let v = if delta < 0 { modulus + delta as i128 } else { delta as i128 };
        U512::from(v as u128)
    }

    #[test]
    fn reversible_pattern_delta_decoder_matches_and_cleans() {
        const W: usize = 16;
        const DBITS: usize = 10;
        let mut b = super::super::B::new();
        let pattern = b.alloc_qubits(W);
        let delta = b.alloc_qubits(DBITS);
        let a_bits = b.alloc_qubits(W);
        emit_pattern_delta_decode_window_for_test(&mut b, &pattern, &delta, &a_bits);
        let forward_ops_len = b.ops.len();
        emit_pattern_delta_decode_window_reverse_for_test(&mut b, &pattern, &delta, &a_bits);
        let total_ccx = count_ccx(&b.ops);
        let forward_ccx = count_ccx(&b.ops[..forward_ops_len]);
        let peak = b.peak_qubits;
        let num_qubits = b.next_qubit as usize;
        let num_bits = b.next_bit as usize;
        let ops = b.ops;
        let cases = [
            (0x0000u16, 1i64),
            (0xffffu16, 1i64),
            (0xa55au16, 7i64),
            (0x5015u16, -3i64),
            (0x8c31u16, 19i64),
        ];
        for &(pat, d0) in &cases {
            let (exp_a, exp_delta) = a_mask_from_pattern_and_delta_for_test(pat, W, d0);
            let mut hasher = sha3::Shake128::default();
            hasher.update(b"by-pattern-delta-decoder-v1");
            let mut xof = hasher.finalize_xof();
            let mut sim = crate::sim::Simulator::new(num_qubits, num_bits, &mut xof);
            set_slice_u512_by(&mut sim, &pattern, U512::from(pat));
            set_slice_u512_by(&mut sim, &delta, twos_u512_for_delta(d0, DBITS));
            sim.apply(&ops[..forward_ops_len]);
            assert_eq!(get_slice_u512_by(&sim, &a_bits), U512::from(exp_a), "A-mask mismatch for {pat:#x}/{d0}");
            assert_eq!(get_slice_u512_by(&sim, &delta), twos_u512_for_delta(exp_delta, DBITS), "delta mismatch for {pat:#x}/{d0}");
            sim.apply(&ops[forward_ops_len..]);
            assert_eq!(get_slice_u512_by(&sim, &pattern), U512::from(pat), "pattern changed");
            assert_eq!(get_slice_u512_by(&sim, &delta), twos_u512_for_delta(d0, DBITS), "delta did not restore");
            assert_eq!(get_slice_u512_by(&sim, &a_bits), U512::ZERO, "A scratch not cleaned");
            assert_eq!(sim.global_phase() & 1, 0, "phase garbage");
        }
        eprintln!(
            "BY reversible pattern+delta decoder: forward_ccx={forward_ccx}, roundtrip_ccx={total_ccx}, peak={peak}q"
        );
        assert!(forward_ccx < 2_000, "pattern decoder synthesis exceeds budget");
    }

    #[test]
    fn pattern_decoder_budget_fits_branch_decode_margin() {
        // Rough reversible budget for pattern+delta -> A-mask+next-delta. Delta
        // stays tiny (empirically |delta|<~32), so a 10-bit signed register is
        // ample. Each of 560 microsteps needs: sign/nonzero predicate, one AND
        // with odd to write A, and a controlled add/negate of a 10-bit delta.
        // This is deliberately pessimistic and still small versus modular
        // replay.
        let delta_bits = 10usize;
        let steps = 560usize;
        let gt_zero_ccx = delta_bits; // sign + nonzero tree, rounded up.
        let write_a_ccx = 1usize;
        let update_delta_ccx = 3 * delta_bits; // controlled +/- and optional negate.
        let per_step = gt_zero_ccx + write_a_ccx + update_delta_ccx;
        let total = per_step * steps;
        eprintln!(
            "BY pattern decoder budget: delta_bits={delta_bits}, per_step≈{per_step}, total≈{total} CCX"
        );
        assert!(total < 30_000, "pattern decoder exceeds reserved branch/decode margin");
    }

    #[test]
    fn window_pattern_and_delta_reconstruct_a_controls() {
        const W: usize = 16;
        let mut hasher = sha3::Shake128::default();
        hasher.update(b"by-pattern-a-mask-v1");
        let mut reader = hasher.finalize_xof();
        let mut buf = [0u8; 24];
        for _ in 0..20_000 {
            reader.read(&mut buf);
            let f_low = (u64::from_le_bytes(buf[0..8].try_into().unwrap()) as i128) | 1;
            let g_low = u64::from_le_bytes(buf[8..16].try_into().unwrap()) as i128;
            let delta = (u64::from_le_bytes(buf[16..24].try_into().unwrap()) % 41) as i64 - 20;
            let controls = branch_controls_for_lowword_window_for_test(W, delta, f_low, g_low);
            let mut pattern = 0u16;
            let mut expected_a = 0u16;
            for (i, &(odd, a)) in controls.iter().enumerate() {
                if odd {
                    pattern |= 1u16 << i;
                }
                if a {
                    expected_a |= 1u16 << i;
                }
            }
            let (got_a, got_delta) = a_mask_from_pattern_and_delta_for_test(pattern, W, delta);
            let direct_delta = jump_matrix_direct_lowword(W, W, delta, f_low, g_low).3.delta_final;
            assert_eq!(got_a, expected_a, "A-mask reconstruction failed");
            assert_eq!(got_delta, direct_delta, "final delta reconstruction failed");
        }
    }

    fn branch_controls_for_lowword_window_for_test(
        w: usize,
        mut delta: i64,
        mut f: i128,
        mut g: i128,
    ) -> Vec<(bool, bool)> {
        let mut out = Vec::with_capacity(w);
        f = truncate_i128(f, w);
        g = truncate_i128(g, w);
        for t in (1..=w).rev() {
            f = truncate_i128(f, t);
            let odd = (g & 1) != 0;
            let a_case = delta > 0 && odd;
            out.push((odd, a_case));
            if a_case {
                let nf = g;
                let ng = (g - f) / 2;
                delta = 1 - delta;
                f = nf;
                g = ng;
            } else if odd {
                g = (g + f) / 2;
                delta = 1 + delta;
            } else {
                g /= 2;
                delta = 1 + delta;
            }
            g = truncate_i128(g, t - 1);
        }
        out
    }

    #[test]
    fn two_adic_branch_generator_matches_classical_prefix_on_small_width() {
        // A real denominator-control generator can be built 2-adically: keep
        // f,g modulo 2^W, use the same swap/neg/add/halve skeleton as the
        // scaled numerator replay, and update the small delta register. This
        // test proves the branch bits are the BY branch bits on a small exact
        // instance. The following budget test explains why this direct version
        // is not the final SOTA implementation.
        const W: usize = 96;
        const STEPS: usize = 64;
        const DBITS: usize = 12;
        let p = U256::from((1u128 << 61) - 1); // odd, tiny relative to 2^W.
        let x = U256::from(0x1234_5678_9abc_defu64) % p;
        let mut delta_c = 1i64;
        let mut f_c = SInt::from_u(p);
        let mut g_c = SInt::from_u(x);
        let mut expected = Vec::with_capacity(STEPS);
        for _ in 0..STEPS {
            let odd = g_c.bit0();
            let a = delta_c > 0 && odd;
            expected.push((odd, a));
            divstep_sint_state(&mut delta_c, &mut f_c, &mut g_c);
        }

        let mut b = super::super::B::new();
        let f = b.alloc_qubits(W);
        let g = b.alloc_qubits(W);
        let delta = b.alloc_qubits(DBITS);
        let odd = b.alloc_qubits(STEPS);
        let a = b.alloc_qubits(STEPS);
        for i in 0..STEPS {
            emit_2adic_by_branch_step_for_test(&mut b, &f, &g, &delta, odd[i], a[i]);
        }
        let ccx = count_ccx(&b.ops);
        let num_qubits = b.next_qubit as usize;
        let num_bits = b.next_bit as usize;
        let ops = b.ops;
        let mut hasher = sha3::Shake128::default();
        hasher.update(b"by-2adic-branch-generator-small-v1");
        let mut xof = hasher.finalize_xof();
        let mut sim = crate::sim::Simulator::new(num_qubits, num_bits, &mut xof);
        set_slice_u512_by(&mut sim, &f, U512::from(p));
        set_slice_u512_by(&mut sim, &g, U512::from(x));
        set_slice_u512_by(&mut sim, &delta, U512::from(1u64));
        sim.apply(&ops);
        for (i, &(odd_e, a_e)) in expected.iter().enumerate() {
            assert_eq!((sim.qubit(odd[i]) & 1) != 0, odd_e, "odd mismatch at step {i}");
            assert_eq!((sim.qubit(a[i]) & 1) != 0, a_e, "A mismatch at step {i}");
        }
        eprintln!(
            "2-adic BY branch generator small prefix: steps={STEPS}, width={W}, ccx={ccx}, peak={}q",
            b.peak_qubits
        );
    }

    fn low_mask_u256_for_test(bits: usize) -> U256 {
        if bits >= 256 {
            U256::MAX
        } else if bits == 0 {
            U256::ZERO
        } else {
            (U256::from(1u64) << bits).wrapping_sub(U256::from(1u64))
        }
    }

    fn truncate_sint_bits_for_test(x: SInt, bits: usize) -> SInt {
        if bits == 0 || x.mag.is_zero() {
            return SInt::zero();
        }
        let mask = low_mask_u256_for_test(bits);
        let mag_low = x.mag & mask;
        let residue = if x.neg {
            if mag_low.is_zero() {
                U256::ZERO
            } else if bits >= 256 {
                U256::ZERO.wrapping_sub(mag_low)
            } else {
                ((U256::from(1u64) << bits).wrapping_sub(mag_low)) & mask
            }
        } else {
            mag_low
        };
        if residue.is_zero() {
            return SInt::zero();
        }
        let sign_bit = residue.bit(bits - 1);
        if sign_bit {
            let mag = if bits >= 256 {
                U256::ZERO.wrapping_sub(residue)
            } else {
                (U256::from(1u64) << bits).wrapping_sub(residue)
            };
            SInt { neg: true, mag }
        } else {
            SInt { neg: false, mag: residue }
        }
    }

    fn divstep_sint_state_truncated_for_test(delta: &mut i64, f: &mut SInt, g: &mut SInt, bits: usize) {
        divstep_sint_state(delta, f, g);
        *f = truncate_sint_bits_for_test(*f, bits);
        *g = truncate_sint_bits_for_test(*g, bits);
    }

    #[test]
    fn fixed_precision_2adic_denominator_branch_curve() {
        // A possible escape from the full 560-bit branch generator is to keep a
        // fixed truncated 2-adic denominator state and accept a small branch
        // mismatch rate. This probe measures the precision curve directly on
        // secp256k1 samples. Result: fixed precision predicts roughly that many
        // initial steps and then loses essentially every 560-step trajectory;
        // field-width denominator state is not an approximate shortcut.
        const STEPS: usize = 560;
        const SAMPLES: usize = 800;
        let precisions = [64usize, 96, 128, 160, 192, 224, 256];
        let mut sampler = Sampler::new(b"by-fixed-precision-branch-curve-v1", SECP256K1_P);
        let mut failure_rates = Vec::with_capacity(precisions.len());
        for &bits in &precisions {
            let mut failures = 0usize;
            let mut first_mismatch_sum = 0usize;
            for _ in 0..SAMPLES {
                let x = sampler.next();
                let mut d_full = 1i64;
                let mut f_full = SInt::from_u(SECP256K1_P);
                let mut g_full = SInt::from_u(x);
                let mut d_local = 1i64;
                let mut f_local = truncate_sint_bits_for_test(f_full, bits);
                let mut g_local = truncate_sint_bits_for_test(g_full, bits);
                let mut first_bad = None;
                for step in 0..STEPS {
                    let odd_full = g_full.bit0();
                    let a_full = d_full > 0 && odd_full;
                    let odd_local = g_local.bit0();
                    let a_local = d_local > 0 && odd_local;
                    if first_bad.is_none() && (odd_full != odd_local || a_full != a_local) {
                        first_bad = Some(step);
                    }
                    divstep_sint_state(&mut d_full, &mut f_full, &mut g_full);
                    divstep_sint_state_truncated_for_test(&mut d_local, &mut f_local, &mut g_local, bits);
                }
                if let Some(step) = first_bad {
                    failures += 1;
                    first_mismatch_sum += step;
                }
            }
            let rate = failures as f64 / SAMPLES as f64;
            let mean_first = if failures == 0 { STEPS as f64 } else { first_mismatch_sum as f64 / failures as f64 };
            failure_rates.push(rate);
            eprintln!(
                "BY fixed-precision branch curve: bits={bits}, mismatch_rate={rate:.4}, mean_first_mismatch={mean_first:.1}"
            );
        }
        assert!(failure_rates[0] > 0.99, "64-bit branch state unexpectedly accurate");
        assert!(failure_rates[6] > 0.99, "256-bit branch state would be an approximate shortcut; revisit");
    }

    fn divstep_i128_exact_for_test(delta: &mut i64, f: &mut i128, g: &mut i128) {
        let odd = (*g & 1) != 0;
        if *delta > 0 && odd {
            let nf = *g;
            let ng = (*g - *f) / 2;
            *delta = 1 - *delta;
            *f = nf;
            *g = ng;
        } else if odd {
            *g = (*g + *f) / 2;
            *delta = 1 + *delta;
        } else {
            *g /= 2;
            *delta = 1 + *delta;
        }
    }

    #[test]
    fn consumed_lowword_window_has_exact_quotient_update_and_pattern_inverse() {
        // This is the algebraic shape of a self-cleaning denominator window.
        // The low W bits choose the branch pattern / matrix. After W divsteps
        // the full signed denominator has been divided by 2^W, so the active
        // high quotient state is
        //
        //     high' = P·high + (P·low)/2^W.
        //
        // The W-bit pattern plus the quotient state is reversible because
        // old = sign(det(P))·adj(P)·new. Thus a production denominator window
        // does not need to preserve the consumed low words separately; the
        // branch pattern is exactly the determinant-history payload.
        const W: usize = 16;
        const SAMPLES: usize = 5_000;
        let scale = 1i128 << W;
        let mut hasher = sha3::Shake128::default();
        hasher.update(b"by-consumed-lowword-window-v1");
        let mut reader = hasher.finalize_xof();
        let mut buf = [0u8; 32];
        let mut max_abs_q = 0i128;
        let mut max_abs_new = 0i128;
        for _ in 0..SAMPLES {
            reader.read(&mut buf);
            let mut f0 = (i128::from_le_bytes(buf[0..16].try_into().unwrap()) >> 40) | 1;
            let mut g0 = i128::from_le_bytes(buf[16..32].try_into().unwrap()) >> 40;
            // Keep values well inside i128 so products by 16-bit coefficients
            // cannot overflow; the identity itself is width-independent.
            f0 %= 1i128 << 86;
            g0 %= 1i128 << 86;
            f0 |= 1;
            let d0 = ((buf[0] as i64) % 41) - 20;

            let low_f = truncate_i128(f0, W);
            let low_g = truncate_i128(g0, W);
            let high_f = (f0 - low_f) / scale;
            let high_g = (g0 - low_g) / scale;
            let bits = branch_bits_for_lowword_window(W, d0, low_f, low_g);
            let m = matrix_from_branch_bits(d0, &bits);

            let q0_num = m.m00 * low_f + m.m01 * low_g;
            let q1_num = m.m10 * low_f + m.m11 * low_g;
            assert_eq!(q0_num % scale, 0, "row0 low correction not integral");
            assert_eq!(q1_num % scale, 0, "row1 low correction not integral");
            let q0 = q0_num / scale;
            let q1 = q1_num / scale;
            max_abs_q = max_abs_q.max(q0.abs()).max(q1.abs());

            let split0 = m.m00 * high_f + m.m01 * high_g + q0;
            let split1 = m.m10 * high_f + m.m11 * high_g + q1;
            let direct0_num = m.m00 * f0 + m.m01 * g0;
            let direct1_num = m.m10 * f0 + m.m11 * g0;
            assert_eq!(direct0_num % scale, 0, "row0 direct update not divisible");
            assert_eq!(direct1_num % scale, 0, "row1 direct update not divisible");
            assert_eq!(split0, direct0_num / scale, "split row0 mismatch");
            assert_eq!(split1, direct1_num / scale, "split row1 mismatch");

            let mut d_run = d0;
            let mut f_run = f0;
            let mut g_run = g0;
            for _ in 0..W {
                divstep_i128_exact_for_test(&mut d_run, &mut f_run, &mut g_run);
            }
            assert_eq!(d_run, m.delta_final, "delta mismatch");
            assert_eq!(f_run, split0, "full divstep f mismatch");
            assert_eq!(g_run, split1, "full divstep g mismatch");
            max_abs_new = max_abs_new.max(f_run.abs()).max(g_run.abs());

            let inv = scaled_inverse_matrix(m, W);
            let rec_f = inv.m00 * f_run + inv.m01 * g_run;
            let rec_g = inv.m10 * f_run + inv.m11 * g_run;
            assert_eq!(rec_f, f0, "pattern+quotient did not reconstruct f");
            assert_eq!(rec_g, g0, "pattern+quotient did not reconstruct g");
        }
        eprintln!(
            "BY consumed lowword window algebra: samples={SAMPLES}, max_abs_q={max_abs_q}, max_abs_new={max_abs_new}"
        );
        assert!(max_abs_q < (1i128 << 17), "lowword quotient correction larger than expected");
    }

    #[test]
    fn tapered_fixed_matrix_denominator_budget_is_sota_shaped_if_selection_solved() {
        // Positive moonshot target: if each 16-step denominator window is
        // applied as the already-known fixed scaled matrix replacement, and
        // the active 2-adic width drops by 16 bits per window, the denominator
        // side is far cheaper than per-bit replay. This intentionally assumes
        // the matrix/pattern selector is already available; it quantifies the
        // arithmetic target for a selected-window implementation.
        const W: usize = 16;
        const START_WIDTH: usize = 560;
        let mut sampler = Sampler::new(b"by-tapered-fixed-den-budget-v1", SECP256K1_P);
        let samples = 8usize;
        let mut total_compute = 0usize;
        let mut max_compute = 0usize;
        let mut max_peak = 0u32;
        for _ in 0..samples {
            let x = sampler.next();
            let mut delta = 1i64;
            let mut f = SInt::from_u(SECP256K1_P);
            let mut g = SInt::from_u(x);
            let mut sample_ccx = 0usize;
            for win in 0..35 {
                let mut bits = Vec::with_capacity(W);
                let d0 = delta;
                for _ in 0..W {
                    bits.push(g.bit0());
                    divstep_sint_state(&mut delta, &mut f, &mut g);
                }
                let m = matrix_from_branch_bits(d0, &bits);
                let width = START_WIDTH - win * W;
                let mut b = super::super::B::new();
                emit_scaled_pair_update_with_cleanup_for_cost(&mut b, m, width, W);
                sample_ccx += count_ccx(&b.ops);
                max_peak = max_peak.max(b.peak_qubits);
            }
            total_compute += sample_ccx;
            max_compute = max_compute.max(sample_ccx);
        }
        let mean_compute = total_compute as f64 / samples as f64;
        let mean_compute_uncompute = 2.0 * mean_compute;
        eprintln!(
            "BY tapered fixed-matrix denominator budget: mean_compute≈{mean_compute:.0}, max_compute={max_compute}, mean_compute_uncompute≈{mean_compute_uncompute:.0}, max_peak={max_peak}q"
        );
        assert!(mean_compute < 450_000.0, "fixed-matrix tapered denominator no longer beats per-bit replay");
        assert!(mean_compute_uncompute < 900_000.0, "denominator compute+uncompute too large for BY SOTA budget");
    }

    #[test]
    fn full_width_denominator_microstep_window_replay_is_not_enough() {
        // Given branch controls, the full denominator can be updated by the
        // same swap/neg/add/halve skeleton as the lowword generator. Measure a
        // 16-step window at the real signed width. This is the straightforward
        // self-cleaning denominator body; it is useful as a target, but too
        // expensive compared with the fixed-matrix/window replacement lower
        // bounds (~8k/window).
        const W: usize = 16;
        const WIDTH: usize = 274;
        let mut b = super::super::B::new();
        let f = b.alloc_qubits(WIDTH);
        let g = b.alloc_qubits(WIDTH);
        let odd = b.alloc_qubits(W);
        let a = b.alloc_qubits(W);
        let start = b.ops.len();
        for i in 0..W {
            emit_2adic_denominator_step_with_controls_for_test(&mut b, &f, &g, odd[i], a[i]);
        }
        let window_ccx = count_ccx(&b.ops[start..]);
        let compute_35 = window_ccx as f64 * 35.0;
        let compute_uncompute = compute_35 * 2.0;
        eprintln!(
            "BY full-width denominator controlled microstep window: window_ccx={window_ccx}, compute35≈{compute_35:.0}, compute_uncompute≈{compute_uncompute:.0}, peak={}q",
            b.peak_qubits
        );
        assert!(window_ccx > 20_000, "full-width denominator replay unexpectedly beats fixed-window target");
        assert!(compute_uncompute > 1_500_000.0, "direct denominator replay might be SOTA-shaped; revisit");
    }

    #[test]
    fn lowword_pattern_oracle_is_cheap_and_clean() {
        // Window-level branch generation component: copy only the low 16 bits
        // of f,g into a scratch 2-adic simulator, run 16 BY steps to produce
        // the odd-pattern, CNOT the pattern to persistent history, then reverse
        // the simulator and clear its local A/pattern scratch. This is the
        // right oracle shape for a windowed DIV. It does not update the full
        // denominator; that selected/window update remains the hard part.
        const W: usize = 16;
        const DBITS: usize = 10;
        let mut b = super::super::B::new();
        let f = b.alloc_qubits(W);
        let g = b.alloc_qubits(W);
        let delta = b.alloc_qubits(DBITS);
        let pattern_tmp = b.alloc_qubits(W);
        let a_tmp = b.alloc_qubits(W);
        let pattern_hist = b.alloc_qubits(W);
        for i in 0..W {
            emit_2adic_by_branch_step_for_test(&mut b, &f, &g, &delta, pattern_tmp[i], a_tmp[i]);
        }
        for i in 0..W {
            b.cx(pattern_tmp[i], pattern_hist[i]);
        }
        for i in (0..W).rev() {
            emit_2adic_by_branch_step_reverse_for_test(&mut b, &f, &g, &delta, pattern_tmp[i], a_tmp[i]);
        }
        let ccx = count_ccx(&b.ops);
        let peak = b.peak_qubits;
        let num_qubits = b.next_qubit as usize;
        let num_bits = b.next_bit as usize;
        let ops = b.ops;
        let cases = [
            (1u64, 3u64, 1i64),
            (0xffffu64, 0x1234u64, -5i64),
            (0x9d31u64 | 1, 0xbeefu64, 17i64),
            (0x8001u64, 0x7fffu64, 0i64),
        ];
        for &(f0, g0, d0) in &cases {
            let bits = branch_bits_for_lowword_window(W, d0, f0 as i128, g0 as i128);
            let mut exp_pat = 0u16;
            for (i, bit) in bits.iter().enumerate() {
                if *bit {
                    exp_pat |= 1u16 << i;
                }
            }
            let mut hasher = sha3::Shake128::default();
            hasher.update(b"by-lowword-pattern-oracle-v1");
            let mut xof = hasher.finalize_xof();
            let mut sim = crate::sim::Simulator::new(num_qubits, num_bits, &mut xof);
            set_slice_u512_by(&mut sim, &f, U512::from(f0));
            set_slice_u512_by(&mut sim, &g, U512::from(g0));
            set_slice_u512_by(&mut sim, &delta, twos_u512_for_delta(d0, DBITS));
            sim.apply(&ops);
            assert_eq!(get_slice_u512_by(&sim, &f), U512::from(f0), "f changed");
            assert_eq!(get_slice_u512_by(&sim, &g), U512::from(g0), "g changed");
            assert_eq!(get_slice_u512_by(&sim, &delta), twos_u512_for_delta(d0, DBITS), "delta changed");
            assert_eq!(get_slice_u512_by(&sim, &pattern_tmp), U512::ZERO, "pattern tmp dirty");
            assert_eq!(get_slice_u512_by(&sim, &a_tmp), U512::ZERO, "A tmp dirty");
            assert_eq!(get_slice_u512_by(&sim, &pattern_hist), U512::from(exp_pat), "pattern mismatch");
            assert_eq!(sim.global_phase() & 1, 0, "phase garbage");
        }
        eprintln!(
            "BY lowword 16-step pattern oracle: ccx={ccx}, peak={peak}q"
        );
        assert!(ccx < 25_000, "lowword pattern oracle too expensive for windowed DIV");
        assert!(peak < 150, "lowword pattern oracle unexpectedly wide");
    }

    #[test]
    fn lowword_pattern_selector_subwindow_budget_probe() {
        // The scratch-600 frontier says the streamed-mask BY replay is only
        // ~65k CCX over the 2.7M target after charging the old W=16 lowword
        // selector.  Before looking for exotic selectors, check whether simply
        // using narrower lowword pattern windows lowers the branch-source cost.
        // This charges both the reversible pattern oracle and the pattern->A
        // delta decoder; the modular replay body is still the measured
        // streamed-mask one, so this is an honest selector-only probe.
        const DBITS: usize = 10;
        let streamed_projected_with_allowance = 2_645_196usize;
        let selector_allowance = 150_000usize;
        let mut totals = Vec::new();
        for &w in &[4usize, 8, 10, 14, 16] {
            assert_eq!(560 % w, 0);
            let mut b = super::super::B::new();
            let f = b.alloc_qubits(w);
            let g = b.alloc_qubits(w);
            let delta = b.alloc_qubits(DBITS);
            let pattern_tmp = b.alloc_qubits(w);
            let a_tmp = b.alloc_qubits(w);
            let pattern_hist = b.alloc_qubits(w);
            for i in 0..w {
                emit_2adic_by_branch_step_for_test(&mut b, &f, &g, &delta, pattern_tmp[i], a_tmp[i]);
            }
            for i in 0..w {
                b.cx(pattern_tmp[i], pattern_hist[i]);
            }
            for i in (0..w).rev() {
                emit_2adic_by_branch_step_reverse_for_test(&mut b, &f, &g, &delta, pattern_tmp[i], a_tmp[i]);
            }
            let oracle_ccx = count_ccx(&b.ops);

            let mut b = super::super::B::new();
            let pattern = b.alloc_qubits(w);
            let delta = b.alloc_qubits(DBITS);
            let a_bits = b.alloc_qubits(w);
            emit_pattern_delta_decode_window_for_test(&mut b, &pattern, &delta, &a_bits);
            let decode_ccx = count_ccx(&b.ops);

            let windows = 560 / w;
            let selector_total = (oracle_ccx + decode_ccx) * windows;
            let projected = streamed_projected_with_allowance - selector_allowance + selector_total;
            let gap = projected as isize - 2_700_000isize;
            eprintln!(
                "BY lowword selector subwindow w={w}: oracle_ccx={oracle_ccx}, decode_ccx={decode_ccx}, windows={windows}, selector_total={selector_total}, projected={projected}, gap={gap}"
            );
            totals.push((w, oracle_ccx, decode_ccx, selector_total, projected, gap));
        }
        let best = *totals.iter().min_by_key(|entry| entry.4).unwrap();
        println!("METRIC scratch600_lowword_selector_best_w={}", best.0);
        println!("METRIC scratch600_lowword_selector_best_oracle_ccx={}", best.1);
        println!("METRIC scratch600_lowword_selector_best_decode_ccx={}", best.2);
        println!("METRIC scratch600_lowword_selector_best_total_ccx={}", best.3);
        println!("METRIC scratch600_lowword_selector_best_projected_toffoli={}", best.4);
        println!("METRIC scratch600_lowword_selector_best_gap_to_2700k={}", best.5);
        for (w, oracle, decode, total, projected, gap) in totals {
            if w == 8 {
                println!("METRIC scratch600_lowword_selector_w8_oracle_ccx={oracle}");
                println!("METRIC scratch600_lowword_selector_w8_decode_ccx={decode}");
                println!("METRIC scratch600_lowword_selector_w8_total_ccx={total}");
                println!("METRIC scratch600_lowword_selector_w8_projected_toffoli={projected}");
                println!("METRIC scratch600_lowword_selector_w8_gap_to_2700k={gap}");
            }
            if w == 16 {
                println!("METRIC scratch600_lowword_selector_w16_oracle_ccx={oracle}");
                println!("METRIC scratch600_lowword_selector_w16_decode_ccx={decode}");
                println!("METRIC scratch600_lowword_selector_w16_total_ccx={total}");
                println!("METRIC scratch600_lowword_selector_w16_projected_toffoli={projected}");
                println!("METRIC scratch600_lowword_selector_w16_gap_to_2700k={gap}");
            }
        }
        assert!(best.3 > 0, "selector accounting should be nonzero");
    }

    #[test]
    fn w4_lowword_selector_naive_full_pair_plumbing_is_too_expensive() {
        // The w=4 lowword oracle+decoder closes the old selector-margin gap
        // only if the missing full-denominator/selector state update costs
        // almost nothing.  Charge the naive exact plumbing: update the signed
        // full-width denominator pair for the same four divsteps under already
        // generated controls.  This is the obvious reversible state source for
        // subsequent lowword windows, and it immediately spends far more than
        // the ~15k slack created by the w=4 selector oracle.  Therefore the
        // w=4 result is not an integration plan by itself; it needs a compact
        // ratio/carry or fixed-matrix state update.
        const W: usize = 4;
        const WIDTH: usize = 274;
        let mut b = super::super::B::new();
        let f = b.alloc_qubits(WIDTH);
        let g = b.alloc_qubits(WIDTH);
        let odd = b.alloc_qubits(W);
        let a = b.alloc_qubits(W);
        let start = b.ops.len();
        for i in 0..W {
            emit_2adic_denominator_step_with_controls_for_test(&mut b, &f, &g, odd[i], a[i]);
        }
        let window_ccx = count_ccx(&b.ops[start..]);
        let windows = 560 / W;
        let compute_ccx = window_ccx * windows;
        let compute_uncompute_ccx = 2 * compute_ccx;
        let w4_selector_gap_surplus = 14_964usize; // projected 2.685036M is this much under 2.7M.
        let compute_excess = compute_ccx - w4_selector_gap_surplus;
        let roundtrip_excess = compute_uncompute_ccx - w4_selector_gap_surplus;
        eprintln!(
            "BY w4 naive full-pair plumbing: window_ccx={window_ccx}, windows={windows}, compute={compute_ccx}, compute_uncompute={compute_uncompute_ccx}, compute_excess={compute_excess}, peak={}q",
            b.peak_qubits
        );
        println!("METRIC scratch600_w4_full_pair_window_ccx={window_ccx}");
        println!("METRIC scratch600_w4_full_pair_compute_ccx={compute_ccx}");
        println!("METRIC scratch600_w4_full_pair_compute_uncompute_ccx={compute_uncompute_ccx}");
        println!("METRIC scratch600_w4_full_pair_compute_excess_ccx={compute_excess}");
        println!("METRIC scratch600_w4_full_pair_roundtrip_excess_ccx={roundtrip_excess}");
        println!("METRIC scratch600_w4_full_pair_peak_q={}", b.peak_qubits);
        assert!(compute_excess > 800_000, "naive full-pair plumbing might fit the w4 slack; revisit BY integration");
    }

    #[test]
    fn w4_fixed_matrix_denominator_update_still_spends_selector_slack() {
        // Next compact-plumbing candidate after the per-bit full-pair replay:
        // assume the 4-step matrix is already selected for free and apply it as
        // one scaled fixed-matrix replacement.  This is an optimistic lower
        // bound for fixed-window denominator plumbing because it ignores QROM /
        // coefficient selection and still uses old+new buffers.  Even this
        // spends far more than the 14,964-CCX slack opened by the w=4 lowword
        // selector oracle, so a viable selector needs a genuinely consumed or
        // ratio-specific update rather than fixed-matrix pair replacement.
        const W: usize = 4;
        const WIDTH: usize = 274;
        const WINDOWS: usize = 560 / W;
        let mut costs = Vec::new();
        let mut max_peak = 0usize;
        let mut sampler = Sampler::new(b"by-w4-fixed-matrix-den-update-v1", SECP256K1_P);
        for _ in 0..8 {
            let x = sampler.next();
            let mut delta = 1i64;
            let mut f = SInt::from_u(SECP256K1_P);
            let mut g = SInt::from_u(x);
            for _ in 0..WINDOWS {
                let f_low = low_signed_sint_for_streaming_test(f, W);
                let g_low = low_signed_sint_for_streaming_test(g, W);
                let (_, _, _, mtx) = jump_matrix_direct_lowword(W, W, delta, f_low, g_low);
                let mut b = super::super::B::new();
                emit_scaled_pair_update_with_cleanup_for_cost(&mut b, mtx, WIDTH, W);
                costs.push(count_ccx(&b.ops));
                max_peak = max_peak.max(b.peak_qubits as usize);
                for _ in 0..W {
                    divstep_sint_state(&mut delta, &mut f, &mut g);
                }
            }
        }
        costs.sort_unstable();
        let mean_window = costs.iter().sum::<usize>() as f64 / costs.len() as f64;
        let p90_window = costs[(costs.len() * 90) / 100];
        let max_window = *costs.last().unwrap();
        let compute_ccx = (mean_window * WINDOWS as f64).round() as usize;
        let compute_uncompute_ccx = 2 * compute_ccx;
        let w4_selector_gap_surplus = 14_964usize;
        let compute_excess = compute_ccx - w4_selector_gap_surplus;
        let roundtrip_excess = compute_uncompute_ccx - w4_selector_gap_surplus;
        eprintln!(
            "BY w4 fixed-matrix denominator update: mean_window={mean_window:.1}, p90={p90_window}, max={max_window}, compute={compute_ccx}, compute_excess={compute_excess}, peak={max_peak}q"
        );
        println!("METRIC scratch600_w4_fixed_matrix_mean_window_ccx={mean_window:.3}");
        println!("METRIC scratch600_w4_fixed_matrix_p90_window_ccx={p90_window}");
        println!("METRIC scratch600_w4_fixed_matrix_max_window_ccx={max_window}");
        println!("METRIC scratch600_w4_fixed_matrix_compute_ccx={compute_ccx}");
        println!("METRIC scratch600_w4_fixed_matrix_compute_uncompute_ccx={compute_uncompute_ccx}");
        println!("METRIC scratch600_w4_fixed_matrix_compute_excess_ccx={compute_excess}");
        println!("METRIC scratch600_w4_fixed_matrix_roundtrip_excess_ccx={roundtrip_excess}");
        println!("METRIC scratch600_w4_fixed_matrix_peak_q={max_peak}");
        assert!(compute_excess > 250_000, "fixed-matrix w4 denominator update might fit; revisit selector integration");
    }

    #[test]
    fn tiny_lowword_selector_windows_still_cannot_fund_fixed_matrix_update() {
        // The w=4 selector oracle alone gets under 2.7M, but only by 14,964
        // CCX and without updating the full denominator state.  Check the even
        // smaller w=1,2 subwindows: if their selector slack still cannot pay
        // for the optimistic fixed-matrix denominator replacement, the simple
        // "smaller window" selector escape hatch is closed.
        const DBITS: usize = 10;
        const WIDTH: usize = 274;
        const GOOGLE_LOW_QUBIT_TOFFOLI: usize = 2_700_000;
        let streamed_projected_with_allowance = 2_645_196usize;
        let selector_allowance = 150_000usize;
        let mut rows = Vec::new();
        for &w in &[1usize, 2, 4] {
            assert_eq!(560 % w, 0);
            let mut b = super::super::B::new();
            let f = b.alloc_qubits(w);
            let g = b.alloc_qubits(w);
            let delta = b.alloc_qubits(DBITS);
            let pattern_tmp = b.alloc_qubits(w);
            let a_tmp = b.alloc_qubits(w);
            let pattern_hist = b.alloc_qubits(w);
            for i in 0..w {
                emit_2adic_by_branch_step_for_test(&mut b, &f, &g, &delta, pattern_tmp[i], a_tmp[i]);
            }
            for i in 0..w {
                b.cx(pattern_tmp[i], pattern_hist[i]);
            }
            for i in (0..w).rev() {
                emit_2adic_by_branch_step_reverse_for_test(&mut b, &f, &g, &delta, pattern_tmp[i], a_tmp[i]);
            }
            let oracle_ccx = count_ccx(&b.ops);

            let mut b = super::super::B::new();
            let pattern = b.alloc_qubits(w);
            let delta = b.alloc_qubits(DBITS);
            let a_bits = b.alloc_qubits(w);
            emit_pattern_delta_decode_window_for_test(&mut b, &pattern, &delta, &a_bits);
            let decode_ccx = count_ccx(&b.ops);

            let windows = 560 / w;
            let selector_total = (oracle_ccx + decode_ccx) * windows;
            let projected = streamed_projected_with_allowance - selector_allowance + selector_total;
            let selector_slack = GOOGLE_LOW_QUBIT_TOFFOLI.saturating_sub(projected);

            let mut sampler = Sampler::new(b"by-tiny-fixed-matrix-update-v1", SECP256K1_P);
            let mut costs = Vec::new();
            for _ in 0..2 {
                let x = sampler.next();
                let mut delta = 1i64;
                let mut f = SInt::from_u(SECP256K1_P);
                let mut g = SInt::from_u(x);
                for _ in 0..windows {
                    let f_low = low_signed_sint_for_streaming_test(f, w);
                    let g_low = low_signed_sint_for_streaming_test(g, w);
                    let (_, _, _, mtx) = jump_matrix_direct_lowword(w, w, delta, f_low, g_low);
                    let mut b = super::super::B::new();
                    emit_scaled_pair_update_with_cleanup_for_cost(&mut b, mtx, WIDTH, w);
                    costs.push(count_ccx(&b.ops));
                    for _ in 0..w {
                        divstep_sint_state(&mut delta, &mut f, &mut g);
                    }
                }
            }
            let mean_window = costs.iter().sum::<usize>() as f64 / costs.len() as f64;
            let compute_ccx = (mean_window * windows as f64).round() as usize;
            let update_excess = compute_ccx as isize - selector_slack as isize;
            eprintln!(
                "BY tiny lowword selector w={w}: selector={selector_total}, projected={projected}, selector_slack={selector_slack}, fixed_matrix_compute={compute_ccx}, update_excess={update_excess}"
            );
            rows.push((w, selector_total, projected, selector_slack, compute_ccx, update_excess));
        }

        let best = *rows.iter().min_by_key(|row| row.5).unwrap();
        println!("METRIC scratch600_tiny_lowword_best_w={}", best.0);
        println!("METRIC scratch600_tiny_lowword_best_selector_ccx={}", best.1);
        println!("METRIC scratch600_tiny_lowword_best_projected={}", best.2);
        println!("METRIC scratch600_tiny_lowword_best_selector_slack={}", best.3);
        println!("METRIC scratch600_tiny_lowword_best_fixed_matrix_compute={}", best.4);
        println!("METRIC scratch600_tiny_lowword_best_update_excess={}", best.5);
        for (w, selector, projected, slack, compute, excess) in rows {
            if w == 1 {
                println!("METRIC scratch600_tiny_lowword_w1_selector_ccx={selector}");
                println!("METRIC scratch600_tiny_lowword_w1_projected={projected}");
                println!("METRIC scratch600_tiny_lowword_w1_selector_slack={slack}");
                println!("METRIC scratch600_tiny_lowword_w1_fixed_matrix_compute={compute}");
                println!("METRIC scratch600_tiny_lowword_w1_update_excess={excess}");
            }
            if w == 2 {
                println!("METRIC scratch600_tiny_lowword_w2_selector_ccx={selector}");
                println!("METRIC scratch600_tiny_lowword_w2_projected={projected}");
                println!("METRIC scratch600_tiny_lowword_w2_selector_slack={slack}");
                println!("METRIC scratch600_tiny_lowword_w2_fixed_matrix_compute={compute}");
                println!("METRIC scratch600_tiny_lowword_w2_update_excess={excess}");
            }
        }
        assert!(best.5 > 0, "tiny lowword windows create enough slack for fixed-matrix update; revisit BY integration");
    }

    #[test]
    fn lowword_pattern_and_q_oracle_is_still_cheap_and_clean() {
        // Strengthen the lowword oracle into the consumed-window primitive:
        // sign-extend the W-bit low words into a slightly wider local simulator,
        // run W divsteps, copy out both the branch pattern and the small
        // quotient correction q=(P·low)/2^W, then reverse the simulator.  This
        // is the local side information needed by the selected quotient update.
        const W: usize = 16;
        const QBITS: usize = 34;
        const DBITS: usize = 10;
        let mut b = super::super::B::new();
        let f = b.alloc_qubits(QBITS);
        let g = b.alloc_qubits(QBITS);
        let delta = b.alloc_qubits(DBITS);
        let pattern_tmp = b.alloc_qubits(W);
        let a_tmp = b.alloc_qubits(W);
        let pattern_hist = b.alloc_qubits(W);
        let q0_hist = b.alloc_qubits(QBITS);
        let q1_hist = b.alloc_qubits(QBITS);
        for i in 0..W {
            emit_signed_by_branch_step_for_test(&mut b, &f, &g, &delta, pattern_tmp[i], a_tmp[i]);
        }
        for i in 0..W {
            b.cx(pattern_tmp[i], pattern_hist[i]);
        }
        for i in 0..QBITS {
            b.cx(f[i], q0_hist[i]);
            b.cx(g[i], q1_hist[i]);
        }
        for i in (0..W).rev() {
            emit_signed_by_branch_step_reverse_for_test(&mut b, &f, &g, &delta, pattern_tmp[i], a_tmp[i]);
        }
        let ccx = count_ccx(&b.ops);
        let peak = b.peak_qubits;
        let num_qubits = b.next_qubit as usize;
        let num_bits = b.next_bit as usize;
        let ops = b.ops;
        let cases = [
            (1i128, 3i128, 1i64),
            (-1i128, 0x1234i128, -5i64),
            (-25323i128, -0x4111i128, 17i64),
            (0x7fffi128, -0x1234i128, 0i64),
        ];
        for &(f0, g0, d0) in &cases {
            let low_f = truncate_i128(f0, W);
            let low_g = truncate_i128(g0, W);
            let bits = branch_bits_for_lowword_window(W, d0, low_f, low_g);
            let m = matrix_from_branch_bits(d0, &bits);
            let scale = 1i128 << W;
            let q0 = (m.m00 * low_f + m.m01 * low_g) / scale;
            let q1 = (m.m10 * low_f + m.m11 * low_g) / scale;
            let mut exp_pat = 0u16;
            for (i, bit) in bits.iter().enumerate() {
                if *bit { exp_pat |= 1u16 << i; }
            }
            let mut hasher = sha3::Shake128::default();
            hasher.update(b"by-lowword-pattern-q-oracle-v1");
            let mut xof = hasher.finalize_xof();
            let mut sim = crate::sim::Simulator::new(num_qubits, num_bits, &mut xof);
            set_slice_u512_by(&mut sim, &f, twos_u512_for_delta(low_f as i64, QBITS));
            set_slice_u512_by(&mut sim, &g, twos_u512_for_delta(low_g as i64, QBITS));
            set_slice_u512_by(&mut sim, &delta, twos_u512_for_delta(d0, DBITS));
            sim.apply(&ops);
            assert_eq!(get_slice_u512_by(&sim, &f), twos_u512_for_delta(low_f as i64, QBITS), "f not restored");
            assert_eq!(get_slice_u512_by(&sim, &g), twos_u512_for_delta(low_g as i64, QBITS), "g not restored");
            assert_eq!(get_slice_u512_by(&sim, &delta), twos_u512_for_delta(d0, DBITS), "delta not restored");
            assert_eq!(get_slice_u512_by(&sim, &pattern_tmp), U512::ZERO, "pattern tmp dirty");
            assert_eq!(get_slice_u512_by(&sim, &a_tmp), U512::ZERO, "A tmp dirty");
            assert_eq!(get_slice_u512_by(&sim, &pattern_hist), U512::from(exp_pat), "pattern mismatch");
            assert_eq!(get_slice_u512_by(&sim, &q0_hist), twos_u512_for_delta(q0 as i64, QBITS), "q0 mismatch");
            assert_eq!(get_slice_u512_by(&sim, &q1_hist), twos_u512_for_delta(q1 as i64, QBITS), "q1 mismatch");
            assert_eq!(sim.global_phase() & 1, 0, "phase garbage");
        }
        eprintln!(
            "BY lowword pattern+q oracle: ccx={ccx}, peak={peak}q, qbits={QBITS}"
        );
        assert!(ccx < 18_000, "pattern+q oracle too expensive for denominator window");
        assert!(peak < 300, "pattern+q oracle unexpectedly wide");
    }

    fn emit_consumed_high_cleanup_row_for_cost(
        b: &mut super::super::B,
        coeff0: i128,
        src0: &[super::super::QubitId],
        coeff1: i128,
        src1: &[super::super::QubitId],
        old_high: &[super::super::QubitId],
        w: usize,
    ) {
        let z_width = old_high.len() + w;
        let z = b.alloc_qubits(z_width);
        add_signed_coeff_times_for_cost(b, coeff0, src0, &z);
        add_signed_coeff_times_for_cost(b, coeff1, src1, &z);
        super::super::add_nbit_const_fast(b, &z, U256::from(1u64 << (w - 1)));
        let z_high = z[w..].to_vec();
        super::super::sub_nbit_qq_fast(b, &z_high, old_high);
        super::super::sub_nbit_const_fast(b, &z, U256::from(1u64 << (w - 1)));
        add_signed_coeff_times_for_cost(b, -coeff1, src1, &z);
        add_signed_coeff_times_for_cost(b, -coeff0, src0, &z);
        b.free_vec(&z);
    }

    fn emit_consumed_high_pair_update_with_q_for_cost(
        b: &mut super::super::B,
        mtx: TransitionMatrix,
        high_width: usize,
        w: usize,
        q_bits: usize,
    ) {
        assert!(high_width > 0);
        let h0 = b.alloc_qubits(high_width);
        let h1 = b.alloc_qubits(high_width);
        let q0 = b.alloc_qubits(q_bits);
        let q1 = b.alloc_qubits(q_bits);
        let y0 = b.alloc_qubits(high_width);
        let y1 = b.alloc_qubits(high_width);

        add_signed_coeff_times_for_cost(b, mtx.m00, &h0, &y0);
        add_signed_coeff_times_for_cost(b, mtx.m01, &h1, &y0);
        add_signed_shifted_term_for_cost(b, &q0, &y0, 0, false);
        add_signed_coeff_times_for_cost(b, mtx.m10, &h0, &y1);
        add_signed_coeff_times_for_cost(b, mtx.m11, &h1, &y1);
        add_signed_shifted_term_for_cost(b, &q1, &y1, 0, false);

        let sgn = det_sign_pow2(mtx, w);
        emit_consumed_high_cleanup_row_for_cost(
            b,
            sgn * mtx.m11,
            &y0,
            -sgn * mtx.m01,
            &y1,
            &h0,
            w,
        );
        emit_consumed_high_cleanup_row_for_cost(
            b,
            -sgn * mtx.m10,
            &y0,
            sgn * mtx.m00,
            &y1,
            &h1,
            w,
        );

        let _ = (h0, h1, q0, q1, y0, y1);
    }

    #[test]
    fn consumed_high_state_window_update_spends_the_selector_margin() {
        // The consumed-lowword identity avoids preserving the low W bits, but
        // the rolling high state still needs a reversible update:
        //
        //   high' = P·high + q,  q = (P·low)/2^W.
        //
        // To clean old high, reconstruct old = s·adj(P)·high' and drop the
        // balanced low W bits.  This charges that high-state body separately
        // from the modular qoffset replay.  It is intentionally optimistic:
        // matrix selection and final q-history cleanup are still uncharged.
        const W: usize = 16;
        const TOTAL_BITS: usize = 35 * W;
        const QBITS: usize = 18;
        const SAMPLES: usize = 4;
        let mut sampler = Sampler::new(b"by-consumed-high-state-update-v1", SECP256K1_P);
        let mut sample_costs = Vec::with_capacity(SAMPLES);
        let mut sample_peaks = Vec::with_capacity(SAMPLES);
        for _ in 0..SAMPLES {
            let x = sampler.next();
            let mut delta = 1i64;
            let mut f = SInt::from_u(SECP256K1_P);
            let mut g = SInt::from_u(x);
            let mut sample_ccx = 0usize;
            let mut sample_peak = 0usize;
            for win in 0..35 {
                let high_width = TOTAL_BITS - (win + 1) * W;
                let f_low = low_signed_sint_for_streaming_test(f, W);
                let g_low = low_signed_sint_for_streaming_test(g, W);
                let (_, _, _, mtx) = jump_matrix_direct_lowword(W, W, delta, f_low, g_low);
                if high_width > 0 {
                    let mut b = super::super::B::new();
                    emit_consumed_high_pair_update_with_q_for_cost(&mut b, mtx, high_width, W, QBITS);
                    sample_ccx += count_ccx(&b.ops);
                    sample_peak = sample_peak.max(b.peak_qubits as usize);
                }
                for _ in 0..W {
                    divstep_sint_state(&mut delta, &mut f, &mut g);
                }
            }
            sample_costs.push(sample_ccx);
            sample_peaks.push(sample_peak);
        }
        sample_costs.sort_unstable();
        sample_peaks.sort_unstable();
        let mean_compute = sample_costs.iter().sum::<usize>() as f64 / SAMPLES as f64;
        let max_compute = sample_costs[SAMPLES - 1];
        let max_peak = sample_peaks[SAMPLES - 1];

        let mut bq = super::super::B::new();
        let f = bq.alloc_qubits(34);
        let g = bq.alloc_qubits(34);
        let delta = bq.alloc_qubits(10);
        let pattern_tmp = bq.alloc_qubits(W);
        let a_tmp = bq.alloc_qubits(W);
        let pattern_hist = bq.alloc_qubits(W);
        let q0_hist = bq.alloc_qubits(34);
        let q1_hist = bq.alloc_qubits(34);
        for i in 0..W {
            emit_signed_by_branch_step_for_test(&mut bq, &f, &g, &delta, pattern_tmp[i], a_tmp[i]);
        }
        for i in 0..W {
            bq.cx(pattern_tmp[i], pattern_hist[i]);
        }
        for i in 0..34 {
            bq.cx(f[i], q0_hist[i]);
            bq.cx(g[i], q1_hist[i]);
        }
        for i in (0..W).rev() {
            emit_signed_by_branch_step_reverse_for_test(&mut bq, &f, &g, &delta, pattern_tmp[i], a_tmp[i]);
        }
        let q_oracle_window = count_ccx(&bq.ops);
        let q_oracle_total = q_oracle_window * 35;
        let optimistic_compute_uncompute = (2.0 * mean_compute).round() as usize;
        let streamed_replay_body_with_allowance = 2_645_196usize;
        let selector_allowance = 150_000usize;
        let decoder_total = 62_160usize;
        let optimistic_pointadd = streamed_replay_body_with_allowance - selector_allowance
            + decoder_total
            + q_oracle_total
            + optimistic_compute_uncompute;
        let optimistic_gap = optimistic_pointadd as isize - 2_700_000isize;
        println!("METRIC by_consumed_high_update_mean_compute_ccx={mean_compute:.3}");
        println!("METRIC by_consumed_high_update_max_compute_ccx={max_compute}");
        println!("METRIC by_consumed_high_update_compute_uncompute_ccx={optimistic_compute_uncompute}");
        println!("METRIC by_consumed_high_update_max_peak_q={max_peak}");
        println!("METRIC by_consumed_high_q_oracle_window_ccx={q_oracle_window}");
        println!("METRIC by_consumed_high_q_oracle_total_ccx={q_oracle_total}");
        println!("METRIC by_consumed_high_optimistic_pointadd_toffoli={optimistic_pointadd}");
        println!("METRIC by_consumed_high_optimistic_gap_to_2700k={optimistic_gap}");
        eprintln!(
            "BY consumed high-state update: mean_compute={mean_compute:.0}, max_compute={max_compute}, compute_uncompute={optimistic_compute_uncompute}, q_oracle_total={q_oracle_total}, optimistic_pointadd={optimistic_pointadd}, gap={optimistic_gap}, max_peak={max_peak}q"
        );
        assert!(
            optimistic_gap > 300_000,
            "consumed high-state update may fit the BY selector margin; revisit integration"
        );
    }

    #[test]
    fn tiny_lowword_consumed_high_update_still_spends_selector_slack() {
        // The w=1/2 selector oracles leave more target slack than w=4, but the
        // consumed-lowword escape still has to roll the full high denominator
        // state.  Charge an optimistic consumed-high update with sampled
        // matrices, current lowword pattern+q oracle, and the same streamed
        // replay body.  This is deliberately generous: matrix selection and
        // q-history cleanup are still not charged.  If it misses here, smaller
        // lowword windows do not rescue the streamed qoffset row.
        const TOTAL_BITS: usize = 560;
        const SAMPLES: usize = 1;
        let streamed_replay_body_with_allowance = 2_645_196usize;
        let selector_allowance = 150_000usize;
        let decoder_total = 62_160usize;
        let mut rows = Vec::new();
        for &w in &[1usize, 2, 4] {
            let windows = TOTAL_BITS / w;
            let q_update_bits = w + 2;
            let q_oracle_bits = 2 * w + 2;
            let mut bq = super::super::B::new();
            let f = bq.alloc_qubits(q_oracle_bits);
            let g = bq.alloc_qubits(q_oracle_bits);
            let delta = bq.alloc_qubits(10);
            let pattern_tmp = bq.alloc_qubits(w);
            let a_tmp = bq.alloc_qubits(w);
            let pattern_hist = bq.alloc_qubits(w);
            let q0_hist = bq.alloc_qubits(q_oracle_bits);
            let q1_hist = bq.alloc_qubits(q_oracle_bits);
            for i in 0..w {
                emit_signed_by_branch_step_for_test(
                    &mut bq,
                    &f,
                    &g,
                    &delta,
                    pattern_tmp[i],
                    a_tmp[i],
                );
            }
            for i in 0..w {
                bq.cx(pattern_tmp[i], pattern_hist[i]);
            }
            for i in 0..q_oracle_bits {
                bq.cx(f[i], q0_hist[i]);
                bq.cx(g[i], q1_hist[i]);
            }
            for i in (0..w).rev() {
                emit_signed_by_branch_step_reverse_for_test(
                    &mut bq,
                    &f,
                    &g,
                    &delta,
                    pattern_tmp[i],
                    a_tmp[i],
                );
            }
            let q_oracle_window = count_ccx(&bq.ops);
            let q_oracle_total = q_oracle_window * windows;

            let mut sampler = Sampler::new(b"by-tiny-consumed-high-update-v1", SECP256K1_P);
            let mut sample_costs = Vec::with_capacity(SAMPLES);
            let mut sample_peaks = Vec::with_capacity(SAMPLES);
            for _ in 0..SAMPLES {
                let x = sampler.next();
                let mut delta = 1i64;
                let mut f = SInt::from_u(SECP256K1_P);
                let mut g = SInt::from_u(x);
                let mut sample_ccx = 0usize;
                let mut sample_peak = 0usize;
                for win in 0..windows {
                    let high_width = TOTAL_BITS - (win + 1) * w;
                    let f_low = low_signed_sint_for_streaming_test(f, w);
                    let g_low = low_signed_sint_for_streaming_test(g, w);
                    let (_, _, _, mtx) = jump_matrix_direct_lowword(w, w, delta, f_low, g_low);
                    if high_width > 0 {
                        let mut b = super::super::B::new();
                        emit_consumed_high_pair_update_with_q_for_cost(
                            &mut b,
                            mtx,
                            high_width,
                            w,
                            q_update_bits,
                        );
                        sample_ccx += count_ccx(&b.ops);
                        sample_peak = sample_peak.max(b.peak_qubits as usize);
                    }
                    for _ in 0..w {
                        divstep_sint_state(&mut delta, &mut f, &mut g);
                    }
                }
                sample_costs.push(sample_ccx);
                sample_peaks.push(sample_peak);
            }
            let mean_compute = sample_costs.iter().sum::<usize>() as f64 / SAMPLES as f64;
            let max_peak = *sample_peaks.iter().max().unwrap();
            let optimistic_compute_uncompute = (2.0 * mean_compute).round() as usize;
            let optimistic_pointadd = streamed_replay_body_with_allowance - selector_allowance
                + decoder_total
                + q_oracle_total
                + optimistic_compute_uncompute;
            let optimistic_gap = optimistic_pointadd as isize - 2_700_000isize;
            eprintln!(
                "BY tiny consumed-high w={w}: q_oracle_window={q_oracle_window}, q_oracle_total={q_oracle_total}, compute={mean_compute:.0}, compute_uncompute={optimistic_compute_uncompute}, pointadd={optimistic_pointadd}, gap={optimistic_gap}, peak={max_peak}"
            );
            rows.push((
                w,
                q_oracle_window,
                q_oracle_total,
                mean_compute.round() as usize,
                optimistic_compute_uncompute,
                optimistic_pointadd,
                optimistic_gap,
                max_peak,
            ));
        }
        let best = *rows.iter().min_by_key(|row| row.6).unwrap();
        println!("METRIC by_tiny_consumed_high_best_w={}", best.0);
        println!("METRIC by_tiny_consumed_high_best_q_oracle_window={}", best.1);
        println!("METRIC by_tiny_consumed_high_best_q_oracle_total={}", best.2);
        println!("METRIC by_tiny_consumed_high_best_compute_ccx={}", best.3);
        println!("METRIC by_tiny_consumed_high_best_compute_uncompute_ccx={}", best.4);
        println!("METRIC by_tiny_consumed_high_best_pointadd={}", best.5);
        println!("METRIC by_tiny_consumed_high_best_gap_to_2700k={}", best.6);
        println!("METRIC by_tiny_consumed_high_best_peak_q={}", best.7);
        for row in rows {
            if row.0 == 1 {
                println!("METRIC by_tiny_consumed_high_w1_pointadd={}", row.5);
                println!("METRIC by_tiny_consumed_high_w1_gap_to_2700k={}", row.6);
                println!("METRIC by_tiny_consumed_high_w1_peak_q={}", row.7);
            }
            if row.0 == 2 {
                println!("METRIC by_tiny_consumed_high_w2_pointadd={}", row.5);
                println!("METRIC by_tiny_consumed_high_w2_gap_to_2700k={}", row.6);
                println!("METRIC by_tiny_consumed_high_w2_peak_q={}", row.7);
            }
        }
        assert!(
            best.6 > 250_000,
            "tiny lowword consumed-high update may fit the streamed qoffset selector budget"
        );
        assert!(
            best.7 > 600,
            "tiny consumed-high update unexpectedly fits the scratch-600 state gate"
        );
    }

    #[test]
    fn by_denominator_branch_history_self_cleans_on_reverse() {
        // This is the constructive counterpart to the dead compute/copy/uncompute
        // branch generators: if the denominator state itself is part of the DIV
        // primitive, its branch history can be cleared while running the
        // denominator backward. That is the self-cleaning shape needed for a
        // real integrated DIV. This small 2-adic circuit proves the mechanism.
        const W: usize = 96;
        const STEPS: usize = 64;
        const DBITS: usize = 12;
        let p = U256::from((1u128 << 61) - 1);
        let x = U256::from(0x0fed_cba9_8765_4321u64) % p;
        let mut b = super::super::B::new();
        let f = b.alloc_qubits(W);
        let g = b.alloc_qubits(W);
        let delta = b.alloc_qubits(DBITS);
        let odd = b.alloc_qubits(STEPS);
        let a = b.alloc_qubits(STEPS);
        for i in 0..STEPS {
            emit_2adic_by_branch_step_for_test(&mut b, &f, &g, &delta, odd[i], a[i]);
        }
        for i in (0..STEPS).rev() {
            emit_2adic_by_branch_step_reverse_for_test(&mut b, &f, &g, &delta, odd[i], a[i]);
        }
        let ccx = count_ccx(&b.ops);
        let peak = b.peak_qubits;
        let num_qubits = b.next_qubit as usize;
        let num_bits = b.next_bit as usize;
        let ops = b.ops;
        let mut hasher = sha3::Shake128::default();
        hasher.update(b"by-den-branch-self-clean-v1");
        let mut xof = hasher.finalize_xof();
        let mut sim = crate::sim::Simulator::new(num_qubits, num_bits, &mut xof);
        set_slice_u512_by(&mut sim, &f, U512::from(p));
        set_slice_u512_by(&mut sim, &g, U512::from(x));
        set_slice_u512_by(&mut sim, &delta, U512::from(1u64));
        sim.apply(&ops);
        assert_eq!(get_slice_u512_by(&sim, &f), U512::from(p));
        assert_eq!(get_slice_u512_by(&sim, &g), U512::from(x));
        assert_eq!(get_slice_u512_by(&sim, &delta), U512::from(1u64));
        for i in 0..STEPS {
            assert_eq!(sim.qubit(odd[i]) & 1, 0, "odd history not clean at {i}");
            assert_eq!(sim.qubit(a[i]) & 1, 0, "A history not clean at {i}");
        }
        assert_eq!(sim.global_phase() & 1, 0, "phase garbage in self-cleaning branch history");
        eprintln!(
            "BY denominator branch history self-clean: steps={STEPS}, width={W}, ccx={ccx}, peak={peak}q"
        );
    }

    #[test]
    fn tapered_2adic_branch_generator_cost_is_still_too_high() {
        // The correct high-precision ratio/denominator state can be tapered:
        // after each branch bit one 2-adic bit is consumed, so the active
        // f/g width drops from 560 to 1. This is the principled version of the
        // direct branch generator and roughly halves the uniform-W cost, but it
        // is still too expensive to sit next to two 1.145M modular replays.
        const STEPS: usize = 560;
        const DBITS: usize = 12;
        let mut b = super::super::B::new();
        let f = b.alloc_qubits(STEPS);
        let g = b.alloc_qubits(STEPS);
        let delta = b.alloc_qubits(DBITS);
        let odd = b.alloc_qubits(STEPS);
        let a = b.alloc_qubits(STEPS);
        let start = b.ops.len();
        for i in 0..STEPS {
            let rem = STEPS - i;
            emit_2adic_by_branch_step_for_test(&mut b, &f[..rem], &g[..rem], &delta, odd[i], a[i]);
        }
        let ccx = count_ccx(&b.ops[start..]);
        let compute_uncompute = 2 * ccx;
        let two_denominators = 2 * compute_uncompute;
        eprintln!(
            "tapered 2-adic BY branch generator: compute_ccx={ccx}, compute_uncompute={compute_uncompute}, two_denominators={two_denominators}, peak={}q",
            b.peak_qubits
        );
        assert!(ccx < 2_000_000, "tapering failed to reduce the direct generator");
        assert!(two_denominators > 2_000_000, "tapered direct generator unexpectedly fits SOTA margin");
    }

    #[test]
    fn inverse_scaled_by_560_negr_frame_recovers_fast_cost() {
        let p = SECP256K1_P;
        let mut sx = Sampler::new(b"by-inverse-negr-560-x-v1", p);
        let mut sq = Sampler::new(b"by-inverse-negr-560-q-v1", p);
        let (x, q, controls, f_final) = loop {
            let x = sx.next();
            let q = sq.next();
            let mut delta = 1i64;
            let mut f = SInt::from_u(p);
            let mut g = SInt::from_u(x);
            let mut controls = Vec::with_capacity(560);
            for _ in 0..560 {
                let odd = g.bit0();
                let a = delta > 0 && odd;
                controls.push((odd, a));
                divstep_sint_state(&mut delta, &mut f, &mut g);
            }
            if g.is_zero() && (f.is_one_pos() || f.is_one_neg()) {
                break (x, q, controls, f);
            }
        };
        let mut b = super::super::B::new();
        let odd = b.alloc_qubits(560);
        let a_ctrl = b.alloc_qubits(560);
        let u = b.alloc_qubits(256);
        let s = b.alloc_qubits(256);
        for i in (0..560).rev() {
            emit_scaled_by_controlled_microstep_inverse_negr_for_test(&mut b, &u, &s, odd[i], a_ctrl[i], p);
        }
        let ccx = count_ccx(&b.ops);
        let peak = b.peak_qubits;
        let num_qubits = b.next_qubit as usize;
        let num_bits = b.next_bit as usize;
        let ops = b.ops;
        let u0 = if f_final.is_one_pos() { negm(q, p) } else { q };
        let expected_s = mulm(q, x, p);
        let mut hasher = sha3::Shake128::default();
        hasher.update(b"by-inverse-negr-560-sim-v1");
        let mut xof = hasher.finalize_xof();
        let mut sim = crate::sim::Simulator::new(num_qubits, num_bits, &mut xof);
        for (i, &(odd_v, a_v)) in controls.iter().enumerate() {
            if odd_v {
                *sim.qubit_mut(odd[i]) |= 1;
            }
            if a_v {
                *sim.qubit_mut(a_ctrl[i]) |= 1;
            }
        }
        set_slice_u512_by(&mut sim, &u, u256_to_u512_for_by_tests(u0));
        set_slice_u512_by(&mut sim, &s, U512::ZERO);
        sim.apply(&ops);
        let got_u = get_slice_u512_by(&sim, &u);
        let p512 = u256_to_u512_for_by_tests(p);
        assert!(got_u == U512::ZERO || got_u == p512, "u=-r was not logical zero: {got_u}");
        assert_eq!(get_slice_u512_by(&sim, &s), u256_to_u512_for_by_tests(expected_s), "product output mismatch");
        eprintln!(
            "BY inverse scaled 560-step neg-r product-clean scaffold: ccx={ccx}, peak={peak}q"
        );
        assert!(ccx <= 1_145_760, "neg-r inverse should match forward cost");
    }

    #[test]
    fn inverse_scaled_by_560_cleans_lam_and_writes_product() {
        // This is the missing pair2 cleanup schedule. If forward scaled BY maps
        // (0, q*x) -> (sign*q, 0), then the inverse map sends (sign*q, 0) ->
        // (0, q*x). Thus pair2 can clean lam and write Ry+Qy without a
        // Kaliski-style inversion or a separate q*x multiplication.
        let p = SECP256K1_P;
        let mut sx = Sampler::new(b"by-inverse-560-x-v1", p);
        let mut sq = Sampler::new(b"by-inverse-560-q-v1", p);
        let (x, q, controls, f_final) = loop {
            let x = sx.next();
            let q = sq.next();
            let mut delta = 1i64;
            let mut f = SInt::from_u(p);
            let mut g = SInt::from_u(x);
            let mut controls = Vec::with_capacity(560);
            for _ in 0..560 {
                let odd = g.bit0();
                let a = delta > 0 && odd;
                controls.push((odd, a));
                divstep_sint_state(&mut delta, &mut f, &mut g);
            }
            if g.is_zero() && (f.is_one_pos() || f.is_one_neg()) {
                break (x, q, controls, f);
            }
        };
        let mut b = super::super::B::new();
        let odd = b.alloc_qubits(560);
        let a_ctrl = b.alloc_qubits(560);
        let r = b.alloc_qubits(256);
        let s = b.alloc_qubits(256);
        for i in (0..560).rev() {
            emit_scaled_by_controlled_microstep_inverse_for_test(&mut b, &r, &s, odd[i], a_ctrl[i], p);
        }
        let ccx = count_ccx(&b.ops);
        let peak = b.peak_qubits;
        let num_qubits = b.next_qubit as usize;
        let num_bits = b.next_bit as usize;
        let ops = b.ops;
        let r0 = if f_final.is_one_pos() { q } else { negm(q, p) };
        let expected_s = mulm(q, x, p);
        let mut hasher = sha3::Shake128::default();
        hasher.update(b"by-inverse-560-sim-v1");
        let mut xof = hasher.finalize_xof();
        let mut sim = crate::sim::Simulator::new(num_qubits, num_bits, &mut xof);
        for (i, &(odd_v, a_v)) in controls.iter().enumerate() {
            if odd_v {
                *sim.qubit_mut(odd[i]) |= 1;
            }
            if a_v {
                *sim.qubit_mut(a_ctrl[i]) |= 1;
            }
        }
        set_slice_u512_by(&mut sim, &r, u256_to_u512_for_by_tests(r0));
        set_slice_u512_by(&mut sim, &s, U512::ZERO);
        sim.apply(&ops);
        let got_r = get_slice_u512_by(&sim, &r);
        let p512 = u256_to_u512_for_by_tests(p);
        assert!(got_r == U512::ZERO || got_r == p512, "lam/r was not logical zero: {got_r}");
        assert_eq!(get_slice_u512_by(&sim, &s), u256_to_u512_for_by_tests(expected_s), "product output mismatch");
        eprintln!(
            "BY inverse scaled 560-step product-clean scaffold: ccx={ccx}, peak={peak}q"
        );
        assert!(ccx < 1_400_000, "inverse scaled BY cleanup too costly");
    }

    fn controls_for_560_sample_for_test(x: U256, p: U256) -> Option<(Vec<(bool, bool)>, SInt)> {
        let mut delta = 1i64;
        let mut f = SInt::from_u(p);
        let mut g = SInt::from_u(x);
        let mut controls = Vec::with_capacity(560);
        for _ in 0..560 {
            let odd = g.bit0();
            let a = delta > 0 && odd;
            controls.push((odd, a));
            divstep_sint_state(&mut delta, &mut f, &mut g);
        }
        if g.is_zero() && (f.is_one_pos() || f.is_one_neg()) {
            Some((controls, f))
        } else {
            None
        }
    }

    #[test]
    fn quantum_branch_values_do_not_reduce_replay_toffoli_accounting() {
        // Important accounting check: the benchmark averages over classical
        // shot masks, not over quantum control values. A CCX controlled by a
        // branch qubit is still a Toffoli gate issued to the quantum computer
        // on every live shot. Therefore branch sparsity does not lower replay
        // Toffoli. This kills a tempting "executed controls are sparse" margin.
        let p = SECP256K1_P;
        let mut b = super::super::B::new();
        let odd = b.alloc_qubits(560);
        let a = b.alloc_qubits(560);
        let r = b.alloc_qubits(256);
        let s = b.alloc_qubits(256);
        for i in 0..560 {
            emit_scaled_by_controlled_microstep_for_test(&mut b, &r, &s, odd[i], a[i], p);
        }
        let static_ccx = count_ccx(&b.ops);
        let num_qubits = b.next_qubit as usize;
        let num_bits = b.next_bit as usize;
        let ops = b.ops;
        let mut sampler = Sampler::new(b"by-executed-fast-replay-v1", p);
        let mut samples = 0usize;
        let mut sum_exec = 0u64;
        let mut max_exec = 0u64;
        while samples < 32 {
            let x = sampler.next();
            let Some((controls, _)) = controls_for_560_sample_for_test(x, p) else { continue; };
            let mut hasher = sha3::Shake128::default();
            hasher.update(b"by-executed-fast-replay-sim-v1");
            let mut xof = hasher.finalize_xof();
            let mut sim = crate::sim::Simulator::new(num_qubits, num_bits, &mut xof);
            for (i, &(odd_v, a_v)) in controls.iter().enumerate() {
                if odd_v { *sim.qubit_mut(odd[i]) |= 1; }
                if a_v { *sim.qubit_mut(a[i]) |= 1; }
            }
            set_slice_u512_by(&mut sim, &s, u256_to_u512_for_by_tests(x));
            sim.apply(&ops);
            sum_exec += sim.stats.toffoli_gates;
            max_exec = max_exec.max(sim.stats.toffoli_gates);
            samples += 1;
        }
        let mean_exec = sum_exec as f64 / samples as f64;
        let mean_per_shot = mean_exec / 64.0;
        let max_per_shot = (max_exec as f64) / 64.0;
        eprintln!(
            "BY fast replay Toffoli accounting: static={static_ccx}, mean_per_shot={mean_per_shot:.1}, max_per_shot={max_per_shot:.1}, samples={samples}"
        );
        assert!((mean_per_shot - static_ccx as f64).abs() < 1.0, "quantum branch controls unexpectedly reduce Toffoli accounting");
    }

    #[test]
    fn exact_cneg_scaled_microstep_may_enable_window_local_a_clearing() {
        let p = SECP256K1_P;
        let inv2 = (p.wrapping_add(U256::from(1u64))) >> 1usize;
        let mut sx = Sampler::new(b"by-window-local-exact-cneg-x-v1", p);
        let mut sy = Sampler::new(b"by-window-local-exact-cneg-y-v1", p);
        let (x, y, controls, boundary_delta, exp_r, exp_s, f_final) = loop {
            let x = sx.next();
            let y = sy.next();
            let mut delta = 1i64;
            let mut f = SInt::from_u(p);
            let mut g = SInt::from_u(x);
            let mut r_exp = U256::ZERO;
            let mut s_exp = addm(y, x, p);
            let mut controls = Vec::with_capacity(560);
            let mut boundary_delta = Vec::with_capacity(35);
            for step in 0..560 {
                if step % 16 == 0 {
                    boundary_delta.push(delta);
                }
                let odd = g.bit0();
                let a = delta > 0 && odd;
                controls.push((odd, a));
                if a {
                    let nr = s_exp;
                    let ns = mulm(subm(s_exp, r_exp, p), inv2, p);
                    r_exp = nr;
                    s_exp = ns;
                } else if odd {
                    s_exp = mulm(addm(s_exp, r_exp, p), inv2, p);
                } else {
                    s_exp = mulm(s_exp, inv2, p);
                }
                divstep_sint_state(&mut delta, &mut f, &mut g);
            }
            if g.is_zero() && (f.is_one_pos() || f.is_one_neg()) {
                break (x, y, controls, boundary_delta, r_exp, s_exp, f);
            }
        };
        let mut b = super::super::B::new();
        let pattern = b.alloc_qubits(560);
        let delta_starts: Vec<Vec<super::super::QubitId>> = (0..35).map(|_| b.alloc_qubits(10)).collect();
        let delta_work = b.alloc_qubits(10);
        let a_window = b.alloc_qubits(16);
        let r = b.alloc_qubits(256);
        let s = b.alloc_qubits(256);
        for win in 0..35 {
            for i in 0..10 {
                b.cx(delta_starts[win][i], delta_work[i]);
            }
            emit_pattern_delta_decode_window_for_test(
                &mut b,
                &pattern[win * 16..win * 16 + 16],
                &delta_work,
                &a_window,
            );
            for i in 0..16 {
                emit_scaled_by_controlled_microstep_exact_cneg_for_test(
                    &mut b,
                    &r,
                    &s,
                    pattern[win * 16 + i],
                    a_window[i],
                    p,
                );
            }
            emit_pattern_delta_decode_window_reverse_for_test(
                &mut b,
                &pattern[win * 16..win * 16 + 16],
                &delta_work,
                &a_window,
            );
            for i in 0..10 {
                b.cx(delta_starts[win][i], delta_work[i]);
            }
        }
        let ccx = count_ccx(&b.ops);
        let peak = b.peak_qubits;
        let num_qubits = b.next_qubit as usize;
        let num_bits = b.next_bit as usize;
        let ops = b.ops;
        let mut hasher = sha3::Shake128::default();
        hasher.update(b"by-window-local-exact-cneg-sim-v1");
        let mut xof = hasher.finalize_xof();
        let mut sim = crate::sim::Simulator::new(num_qubits, num_bits, &mut xof);
        for (i, &(odd_v, _)) in controls.iter().enumerate() {
            if odd_v {
                *sim.qubit_mut(pattern[i]) |= 1;
            }
        }
        for (win, &d) in boundary_delta.iter().enumerate() {
            set_slice_u512_by(&mut sim, &delta_starts[win], twos_u512_for_delta(d, 10));
        }
        set_slice_u512_by(&mut sim, &r, U512::ZERO);
        set_slice_u512_by(&mut sim, &s, u256_to_u512_for_by_tests(addm(y, x, p)));
        sim.apply(&ops);
        assert_eq!(get_slice_u512_by(&sim, &r), u256_to_u512_for_by_tests(exp_r), "r mismatch");
        assert_eq!(get_slice_u512_by(&sim, &s), u256_to_u512_for_by_tests(exp_s), "s mismatch");
        let plus_one = if f_final.is_one_pos() { exp_r } else { negm(exp_r, p) };
        let quotient = subm(plus_one, U256::from(1u64), p);
        assert_eq!(quotient, mulm(y, fermat_modinv(x, p), p), "tagged quotient mismatch");
        assert_eq!(get_slice_u512_by(&sim, &a_window), U512::ZERO, "A window scratch not clean");
        assert_eq!(get_slice_u512_by(&sim, &delta_work), U512::ZERO, "delta work not clean");
        let phase = sim.global_phase() & 1;
        eprintln!(
            "BY window-local exact-cneg replay: ccx={ccx}, peak={peak}q, phase={phase}"
        );
        assert_eq!(phase, 0, "exact A-controlled negation did not fix early A-clear phase");
        assert!(ccx < 1_550_000, "exact-cneg window-local replay too expensive to keep alive");
    }

    #[test]
    fn exact_scaled_microstep_is_phase_safe_but_too_expensive_for_window_local_clear() {
        // Exact reversible modular arithmetic should allow controls to be
        // cleared locally, unlike the MBU fast path. Quantify the tax before
        // pursuing a full exact replay.
        let p = SECP256K1_P;
        let mut b = super::super::B::new();
        let odd = b.alloc_qubit();
        let a = b.alloc_qubit();
        let r = b.alloc_qubits(256);
        let s = b.alloc_qubits(256);
        emit_scaled_by_controlled_microstep_exact_for_test(&mut b, &r, &s, odd, a, p);
        let exact_ccx = count_ccx(&b.ops);
        let exact_peak = b.peak_qubits;
        let mut bf = super::super::B::new();
        let oddf = bf.alloc_qubit();
        let af = bf.alloc_qubit();
        let rf = bf.alloc_qubits(256);
        let sf = bf.alloc_qubits(256);
        emit_scaled_by_controlled_microstep_for_test(&mut bf, &rf, &sf, oddf, af, p);
        let fast_ccx = count_ccx(&bf.ops);
        let approx560 = exact_ccx as f64 * 560.0;
        eprintln!(
            "BY exact scaled microstep: exact_ccx={exact_ccx}, fast_ccx={fast_ccx}, approx560≈{approx560:.0}, peak={exact_peak}q"
        );
        assert!(exact_ccx > fast_ccx + 1_000, "exact microstep tax unexpectedly small");
        assert!(approx560 > 2_000_000.0, "exact replay might be SOTA-shaped; revisit window-local clear");
    }

    #[test]
    fn window_local_a_clear_fails_phase_with_mbu_microsteps() {
        // Tempting low-scratch schedule: decode one 16-step A window, use it
        // immediately, then reverse the decoder and clear A before later
        // windows. Classical data comes out right, but current scaled-BY
        // microsteps use measurement-based modular arithmetic whose phase
        // corrections still depend on those controls. Clearing A early leaves
        // phase garbage. Keep this as an executable invalidation: with the
        // current MBU microsteps, A controls must remain live until the replay
        // is complete, or the modular primitives must be made exact/phase-safe
        // under early control clearing.
        let p = SECP256K1_P;
        let inv2 = (p.wrapping_add(U256::from(1u64))) >> 1usize;
        let mut sx = Sampler::new(b"by-window-local-decoder-560-x-v1", p);
        let mut sy = Sampler::new(b"by-window-local-decoder-560-y-v1", p);
        let (x, y, controls, boundary_delta, exp_r, exp_s, f_final) = loop {
            let x = sx.next();
            let y = sy.next();
            let mut delta = 1i64;
            let mut f = SInt::from_u(p);
            let mut g = SInt::from_u(x);
            let mut r_exp = U256::ZERO;
            let mut s_exp = addm(y, x, p);
            let mut controls = Vec::with_capacity(560);
            let mut boundary_delta = Vec::with_capacity(35);
            for step in 0..560 {
                if step % 16 == 0 {
                    boundary_delta.push(delta);
                }
                let odd = g.bit0();
                let a = delta > 0 && odd;
                controls.push((odd, a));
                if a {
                    let nr = s_exp;
                    let ns = mulm(subm(s_exp, r_exp, p), inv2, p);
                    r_exp = nr;
                    s_exp = ns;
                } else if odd {
                    s_exp = mulm(addm(s_exp, r_exp, p), inv2, p);
                } else {
                    s_exp = mulm(s_exp, inv2, p);
                }
                divstep_sint_state(&mut delta, &mut f, &mut g);
            }
            if g.is_zero() && (f.is_one_pos() || f.is_one_neg()) {
                break (x, y, controls, boundary_delta, r_exp, s_exp, f);
            }
        };

        let mut b = super::super::B::new();
        let pattern = b.alloc_qubits(560);
        let delta_starts: Vec<Vec<super::super::QubitId>> = (0..35).map(|_| b.alloc_qubits(10)).collect();
        let delta_work = b.alloc_qubits(10);
        let a_window = b.alloc_qubits(16);
        let r = b.alloc_qubits(256);
        let s = b.alloc_qubits(256);
        for win in 0..35 {
            for i in 0..10 {
                b.cx(delta_starts[win][i], delta_work[i]);
            }
            emit_pattern_delta_decode_window_for_test(
                &mut b,
                &pattern[win * 16..win * 16 + 16],
                &delta_work,
                &a_window,
            );
            for i in 0..16 {
                emit_scaled_by_controlled_microstep_for_test(
                    &mut b,
                    &r,
                    &s,
                    pattern[win * 16 + i],
                    a_window[i],
                    p,
                );
            }
            emit_pattern_delta_decode_window_reverse_for_test(
                &mut b,
                &pattern[win * 16..win * 16 + 16],
                &delta_work,
                &a_window,
            );
            for i in 0..10 {
                b.cx(delta_starts[win][i], delta_work[i]);
            }
        }
        let ccx = count_ccx(&b.ops);
        let peak = b.peak_qubits;
        let num_qubits = b.next_qubit as usize;
        let num_bits = b.next_bit as usize;
        let ops = b.ops;
        let mut hasher = sha3::Shake128::default();
        hasher.update(b"by-window-local-decoder-560-sim-v1");
        let mut xof = hasher.finalize_xof();
        let mut sim = crate::sim::Simulator::new(num_qubits, num_bits, &mut xof);
        for (i, &(odd_v, _)) in controls.iter().enumerate() {
            if odd_v {
                *sim.qubit_mut(pattern[i]) |= 1;
            }
        }
        for (win, &d) in boundary_delta.iter().enumerate() {
            set_slice_u512_by(&mut sim, &delta_starts[win], twos_u512_for_delta(d, 10));
        }
        set_slice_u512_by(&mut sim, &r, U512::ZERO);
        set_slice_u512_by(&mut sim, &s, u256_to_u512_for_by_tests(addm(y, x, p)));
        sim.apply(&ops);
        assert_eq!(get_slice_u512_by(&sim, &r), u256_to_u512_for_by_tests(exp_r), "r mismatch");
        assert_eq!(get_slice_u512_by(&sim, &s), u256_to_u512_for_by_tests(exp_s), "s mismatch");
        let plus_one = if f_final.is_one_pos() { exp_r } else { negm(exp_r, p) };
        let quotient = subm(plus_one, U256::from(1u64), p);
        assert_eq!(quotient, mulm(y, fermat_modinv(x, p), p), "tagged quotient mismatch");
        assert_eq!(get_slice_u512_by(&sim, &a_window), U512::ZERO, "A window scratch not clean");
        assert_eq!(get_slice_u512_by(&sim, &delta_work), U512::ZERO, "delta work not clean");
        for (win, &d) in boundary_delta.iter().enumerate() {
            assert_eq!(get_slice_u512_by(&sim, &delta_starts[win]), twos_u512_for_delta(d, 10), "boundary delta changed");
        }
        let phase = sim.global_phase() & 1;
        eprintln!(
            "BY window-local A-clear phase failure: ccx={ccx}, peak={peak}q, boundary_delta_bits=350, phase={phase}"
        );
        assert_ne!(phase, 0, "window-local A clearing unexpectedly phase-clean; revisit low-scratch schedule");
        assert!(peak < 2_300, "window-local schedule did not reduce peak enough to be worth revisiting");
    }

    #[test]
    fn live_reduction_flags_make_window_local_a_clear_phase_safe_candidate() {
        // Keep the modular-add reduction flags live instead of immediately
        // uncomputing them with measurement-based cmp_lt. This tests whether
        // those flags are the phase dependency that made early A-clearing fail.
        // The flags are deliberate garbage here; the next problem is a cheap
        // way to clean or absorb them.
        let p = SECP256K1_P;
        let inv2 = (p.wrapping_add(U256::from(1u64))) >> 1usize;
        let mut sx = Sampler::new(b"by-live-flags-window-local-x-v1", p);
        let mut sy = Sampler::new(b"by-live-flags-window-local-y-v1", p);
        let (x, y, controls, boundary_delta, exp_r, exp_s, f_final) = loop {
            let x = sx.next();
            let y = sy.next();
            let mut delta = 1i64;
            let mut f = SInt::from_u(p);
            let mut g = SInt::from_u(x);
            let mut r_exp = U256::ZERO;
            let mut s_exp = addm(y, x, p);
            let mut controls = Vec::with_capacity(560);
            let mut boundary_delta = Vec::with_capacity(35);
            for step in 0..560 {
                if step % 16 == 0 { boundary_delta.push(delta); }
                let odd = g.bit0();
                let a = delta > 0 && odd;
                controls.push((odd, a));
                if a {
                    let nr = s_exp;
                    let ns = mulm(subm(s_exp, r_exp, p), inv2, p);
                    r_exp = nr;
                    s_exp = ns;
                } else if odd {
                    s_exp = mulm(addm(s_exp, r_exp, p), inv2, p);
                } else {
                    s_exp = mulm(s_exp, inv2, p);
                }
                divstep_sint_state(&mut delta, &mut f, &mut g);
            }
            if g.is_zero() && (f.is_one_pos() || f.is_one_neg()) {
                break (x, y, controls, boundary_delta, r_exp, s_exp, f);
            }
        };

        let mut b = super::super::B::new();
        let pattern = b.alloc_qubits(560);
        let red_flags = b.alloc_qubits(560);
        let delta_starts: Vec<Vec<super::super::QubitId>> = (0..35).map(|_| b.alloc_qubits(10)).collect();
        let delta_work = b.alloc_qubits(10);
        let a_window = b.alloc_qubits(16);
        let r = b.alloc_qubits(256);
        let s = b.alloc_qubits(256);
        for win in 0..35 {
            for i in 0..10 { b.cx(delta_starts[win][i], delta_work[i]); }
            emit_pattern_delta_decode_window_for_test(
                &mut b,
                &pattern[win * 16..win * 16 + 16],
                &delta_work,
                &a_window,
            );
            for i in 0..16 {
                let step = win * 16 + i;
                emit_scaled_by_controlled_microstep_live_addflag_for_test(
                    &mut b,
                    &r,
                    &s,
                    pattern[step],
                    a_window[i],
                    red_flags[step],
                    p,
                );
            }
            emit_pattern_delta_decode_window_reverse_for_test(
                &mut b,
                &pattern[win * 16..win * 16 + 16],
                &delta_work,
                &a_window,
            );
            for i in 0..10 { b.cx(delta_starts[win][i], delta_work[i]); }
        }
        let ccx = count_ccx(&b.ops);
        let peak = b.peak_qubits;
        let num_qubits = b.next_qubit as usize;
        let num_bits = b.next_bit as usize;
        let ops = b.ops;
        let mut hasher = sha3::Shake128::default();
        hasher.update(b"by-live-flags-window-local-sim-v1");
        let mut xof = hasher.finalize_xof();
        let mut sim = crate::sim::Simulator::new(num_qubits, num_bits, &mut xof);
        for (i, &(odd_v, _)) in controls.iter().enumerate() {
            if odd_v { *sim.qubit_mut(pattern[i]) |= 1; }
        }
        for (win, &d) in boundary_delta.iter().enumerate() {
            set_slice_u512_by(&mut sim, &delta_starts[win], twos_u512_for_delta(d, 10));
        }
        set_slice_u512_by(&mut sim, &r, U512::ZERO);
        set_slice_u512_by(&mut sim, &s, u256_to_u512_for_by_tests(addm(y, x, p)));
        sim.apply(&ops);
        assert_eq!(get_slice_u512_by(&sim, &r), u256_to_u512_for_by_tests(exp_r), "r mismatch");
        assert_eq!(get_slice_u512_by(&sim, &s), u256_to_u512_for_by_tests(exp_s), "s mismatch");
        let plus_one = if f_final.is_one_pos() { exp_r } else { negm(exp_r, p) };
        let quotient = subm(plus_one, U256::from(1u64), p);
        assert_eq!(quotient, mulm(y, fermat_modinv(x, p), p), "tagged quotient mismatch");
        assert_eq!(get_slice_u512_by(&sim, &a_window), U512::ZERO, "A window scratch not clean");
        assert_eq!(get_slice_u512_by(&sim, &delta_work), U512::ZERO, "delta work not clean");
        let red_nonzero = red_flags.iter().any(|&q| (sim.qubit(q) & 1) != 0);
        let phase = sim.global_phase() & 1;
        eprintln!(
            "BY live reduction flags + window-local A clear: ccx={ccx}, peak={peak}q, phase={phase}, flags_nonzero={red_nonzero}"
        );
        assert!(red_nonzero, "live reduction flags unexpectedly clean already");
        assert_eq!(phase, 0, "live reduction flags did not fix early A-clear phase");
        assert!(ccx < 1_140_000, "live-flag window-local replay lost target Toffoli band");
        assert!(peak < 2_850, "live-flag window-local replay exceeds current cap");
    }

    #[test]
    fn scaled_by_pattern_decoder_560_tagged_div_scaffold_is_clean() {
        // Clean version of the raw-pattern scaffold: expand 560 odd-pattern
        // bits into A controls using the reversible pattern+delta decoder,
        // run the scaled-BY replay, then reverse the decoders to clean A and
        // restore delta. This is deliberately not the final low-scratch
        // schedule (it keeps all A controls during replay), but it proves the
        // decoder integrates with the 560-step arithmetic and quantifies the
        // exact overhead.
        let p = SECP256K1_P;
        let inv2 = (p.wrapping_add(U256::from(1u64))) >> 1usize;
        let mut sx = Sampler::new(b"by-pattern-decoder-560-x-v1", p);
        let mut sy = Sampler::new(b"by-pattern-decoder-560-y-v1", p);
        let (x, y, controls, exp_r, exp_s, f_final) = loop {
            let x = sx.next();
            let y = sy.next();
            let mut delta = 1i64;
            let mut f = SInt::from_u(p);
            let mut g = SInt::from_u(x);
            let mut r_exp = U256::ZERO;
            let mut s_exp = addm(y, x, p);
            let mut controls = Vec::with_capacity(560);
            for _ in 0..560 {
                let odd = g.bit0();
                let a = delta > 0 && odd;
                controls.push((odd, a));
                if a {
                    let nr = s_exp;
                    let ns = mulm(subm(s_exp, r_exp, p), inv2, p);
                    r_exp = nr;
                    s_exp = ns;
                } else if odd {
                    s_exp = mulm(addm(s_exp, r_exp, p), inv2, p);
                } else {
                    s_exp = mulm(s_exp, inv2, p);
                }
                divstep_sint_state(&mut delta, &mut f, &mut g);
            }
            if g.is_zero() && (f.is_one_pos() || f.is_one_neg()) {
                break (x, y, controls, r_exp, s_exp, f);
            }
        };

        let mut b = super::super::B::new();
        let pattern = b.alloc_qubits(560);
        let a_hist = b.alloc_qubits(560);
        let delta = b.alloc_qubits(10);
        let r = b.alloc_qubits(256);
        let s = b.alloc_qubits(256);
        for win in 0..35 {
            emit_pattern_delta_decode_window_for_test(
                &mut b,
                &pattern[win * 16..win * 16 + 16],
                &delta,
                &a_hist[win * 16..win * 16 + 16],
            );
        }
        let decode_forward_ccx = count_ccx(&b.ops);
        for i in 0..560 {
            emit_scaled_by_controlled_microstep_for_test(&mut b, &r, &s, pattern[i], a_hist[i], p);
        }
        let replay_plus_decode_ccx = count_ccx(&b.ops);
        for win in (0..35).rev() {
            emit_pattern_delta_decode_window_reverse_for_test(
                &mut b,
                &pattern[win * 16..win * 16 + 16],
                &delta,
                &a_hist[win * 16..win * 16 + 16],
            );
        }
        let ccx = count_ccx(&b.ops);
        let peak = b.peak_qubits;
        let num_qubits = b.next_qubit as usize;
        let num_bits = b.next_bit as usize;
        let ops = b.ops;
        let mut hasher = sha3::Shake128::default();
        hasher.update(b"by-pattern-decoder-560-sim-v1");
        let mut xof = hasher.finalize_xof();
        let mut sim = crate::sim::Simulator::new(num_qubits, num_bits, &mut xof);
        for (i, &(odd_v, _)) in controls.iter().enumerate() {
            if odd_v {
                *sim.qubit_mut(pattern[i]) |= 1;
            }
        }
        set_slice_u512_by(&mut sim, &delta, U512::from(1u64));
        set_slice_u512_by(&mut sim, &r, U512::ZERO);
        set_slice_u512_by(&mut sim, &s, u256_to_u512_for_by_tests(addm(y, x, p)));
        sim.apply(&ops);
        assert_eq!(get_slice_u512_by(&sim, &r), u256_to_u512_for_by_tests(exp_r), "r mismatch");
        assert_eq!(get_slice_u512_by(&sim, &s), u256_to_u512_for_by_tests(exp_s), "s mismatch");
        assert_eq!(exp_s, U256::ZERO, "bottom tagged channel did not zero");
        let plus_one = if f_final.is_one_pos() { exp_r } else { negm(exp_r, p) };
        let quotient = subm(plus_one, U256::from(1u64), p);
        assert_eq!(quotient, mulm(y, fermat_modinv(x, p), p), "tagged quotient mismatch");
        assert_eq!(get_slice_u512_by(&sim, &a_hist), U512::ZERO, "A history not cleaned");
        assert_eq!(get_slice_u512_by(&sim, &delta), U512::from(1u64), "delta not restored");
        assert_eq!(sim.global_phase() & 1, 0, "phase garbage");
        eprintln!(
            "BY pattern-decoder 560-step tagged-DIV scaffold: decode_ccx={decode_forward_ccx}, replay_plus_decode_ccx={replay_plus_decode_ccx}, roundtrip_ccx={ccx}, peak={peak}q"
        );
        assert!(ccx < 1_300_000, "decoded pattern scaffold exceeded integration margin");
        assert!(peak < 2_700, "decoded pattern scaffold too wide for current cap");
    }

    #[test]
    fn scaled_by_pattern_history_560_tagged_div_scaffold_reduces_peak() {
        let p = SECP256K1_P;
        let inv2 = (p.wrapping_add(U256::from(1u64))) >> 1usize;
        let mut sx = Sampler::new(b"by-pattern-560-sim-x-v1", p);
        let mut sy = Sampler::new(b"by-pattern-560-sim-y-v1", p);
        let (x, y, controls, exp_r, exp_s, f_final) = loop {
            let x = sx.next();
            let y = sy.next();
            let mut delta = 1i64;
            let mut f = SInt::from_u(p);
            let mut g = SInt::from_u(x);
            let mut r_exp = U256::ZERO;
            let mut s_exp = addm(y, x, p);
            let mut controls = Vec::with_capacity(560);
            for _ in 0..560 {
                let odd = g.bit0();
                let a = delta > 0 && odd;
                controls.push((odd, a));
                if a {
                    let nr = s_exp;
                    let ns = mulm(subm(s_exp, r_exp, p), inv2, p);
                    r_exp = nr;
                    s_exp = ns;
                } else if odd {
                    s_exp = mulm(addm(s_exp, r_exp, p), inv2, p);
                } else {
                    s_exp = mulm(s_exp, inv2, p);
                }
                divstep_sint_state(&mut delta, &mut f, &mut g);
            }
            if g.is_zero() && (f.is_one_pos() || f.is_one_neg()) {
                break (x, y, controls, r_exp, s_exp, f);
            }
        };

        let mut b = super::super::B::new();
        let pattern = b.alloc_qubits(560); // raw pattern history; compressed IDs are a later decoder layer.
        let a_scratch = b.alloc_qubits(16);
        let r = b.alloc_qubits(256);
        let s = b.alloc_qubits(256);
        for win in 0..35 {
            for i in 0..16 {
                if controls[win * 16 + i].1 {
                    b.x(a_scratch[i]);
                }
            }
            for i in 0..16 {
                emit_scaled_by_controlled_microstep_for_test(
                    &mut b,
                    &r,
                    &s,
                    pattern[win * 16 + i],
                    a_scratch[i],
                    p,
                );
            }
            for i in 0..16 {
                if controls[win * 16 + i].1 {
                    b.x(a_scratch[i]);
                }
            }
        }
        let ccx = count_ccx(&b.ops);
        let peak = b.peak_qubits;
        let num_qubits = b.next_qubit as usize;
        let num_bits = b.next_bit as usize;
        let ops = b.ops;
        let mut hasher = sha3::Shake128::default();
        hasher.update(b"by-pattern-560-sim-xof-v1");
        let mut xof = hasher.finalize_xof();
        let mut sim = crate::sim::Simulator::new(num_qubits, num_bits, &mut xof);
        for (i, &(odd_v, _)) in controls.iter().enumerate() {
            if odd_v {
                *sim.qubit_mut(pattern[i]) |= 1;
            }
        }
        set_slice_u512_by(&mut sim, &r, U512::ZERO);
        set_slice_u512_by(&mut sim, &s, u256_to_u512_for_by_tests(addm(y, x, p)));
        sim.apply(&ops);
        assert_eq!(get_slice_u512_by(&sim, &r), u256_to_u512_for_by_tests(exp_r), "r mismatch");
        assert_eq!(get_slice_u512_by(&sim, &s), u256_to_u512_for_by_tests(exp_s), "s mismatch");
        assert_eq!(exp_s, U256::ZERO, "bottom tagged channel did not zero");
        let plus_one = if f_final.is_one_pos() { exp_r } else { negm(exp_r, p) };
        let quotient = subm(plus_one, U256::from(1u64), p);
        assert_eq!(quotient, mulm(y, fermat_modinv(x, p), p), "tagged quotient mismatch");
        eprintln!(
            "BY pattern-history 560-step tagged-DIV scaffold: ccx={ccx}, peak={peak}q, raw_pattern_bits=560"
        );
        assert!(peak < 1_900, "raw pattern-history scaffold peak too high");
    }

    #[test]
    fn scaled_by_controlled_560_tagged_div_basis_simulation() {
        let p = SECP256K1_P;
        let inv2 = (p.wrapping_add(U256::from(1u64))) >> 1usize;
        let mut sx = Sampler::new(b"by-scaled-560-sim-x-v1", p);
        let mut sy = Sampler::new(b"by-scaled-560-sim-y-v1", p);
        let (x, y, controls, exp_r, exp_s, f_final) = loop {
            let x = sx.next();
            let y = sy.next();
            let mut delta = 1i64;
            let mut f = SInt::from_u(p);
            let mut g = SInt::from_u(x);
            let mut r_exp = U256::ZERO;
            let mut s_exp = addm(y, x, p);
            let mut controls = Vec::with_capacity(560);
            for _ in 0..560 {
                let odd = g.bit0();
                let a = delta > 0 && odd;
                controls.push((odd, a));
                if a {
                    let nr = s_exp;
                    let ns = mulm(subm(s_exp, r_exp, p), inv2, p);
                    r_exp = nr;
                    s_exp = ns;
                } else if odd {
                    s_exp = mulm(addm(s_exp, r_exp, p), inv2, p);
                } else {
                    s_exp = mulm(s_exp, inv2, p);
                }
                divstep_sint_state(&mut delta, &mut f, &mut g);
            }
            if g.is_zero() && (f.is_one_pos() || f.is_one_neg()) {
                break (x, y, controls, r_exp, s_exp, f);
            }
        };

        let mut b = super::super::B::new();
        let odd = b.alloc_qubits(560);
        let a_ctrl = b.alloc_qubits(560);
        let r = b.alloc_qubits(256);
        let s = b.alloc_qubits(256);
        for i in 0..560 {
            emit_scaled_by_controlled_microstep_for_test(&mut b, &r, &s, odd[i], a_ctrl[i], p);
        }
        let ccx = count_ccx(&b.ops);
        let peak = b.peak_qubits;
        let num_qubits = b.next_qubit as usize;
        let num_bits = b.next_bit as usize;
        let ops = b.ops;
        let mut hasher = sha3::Shake128::default();
        hasher.update(b"by-scaled-560-sim-xof-v1");
        let mut xof = hasher.finalize_xof();
        let mut sim = crate::sim::Simulator::new(num_qubits, num_bits, &mut xof);
        for (i, &(odd_v, a_v)) in controls.iter().enumerate() {
            if odd_v {
                *sim.qubit_mut(odd[i]) |= 1;
            }
            if a_v {
                *sim.qubit_mut(a_ctrl[i]) |= 1;
            }
        }
        set_slice_u512_by(&mut sim, &r, U512::ZERO);
        set_slice_u512_by(&mut sim, &s, u256_to_u512_for_by_tests(addm(y, x, p)));
        sim.apply(&ops);
        assert_eq!(get_slice_u512_by(&sim, &r), u256_to_u512_for_by_tests(exp_r), "r mismatch");
        assert_eq!(get_slice_u512_by(&sim, &s), u256_to_u512_for_by_tests(exp_s), "s mismatch");
        assert_eq!(exp_s, U256::ZERO, "bottom tagged channel did not zero");
        let plus_one = if f_final.is_one_pos() { exp_r } else { negm(exp_r, p) };
        let quotient = subm(plus_one, U256::from(1u64), p);
        assert_eq!(quotient, mulm(y, fermat_modinv(x, p), p), "tagged quotient mismatch");
        eprintln!(
            "BY scaled controlled 560-step tagged-DIV basis sim: ccx={ccx}, peak={peak}q"
        );
    }

    #[test]
    fn scaled_by_controlled_560_scaffold_cost_model_fits_current_cap() {
        let p = SECP256K1_P;
        let mut b = super::super::B::new();
        let odd = b.alloc_qubits(560);
        let a_ctrl = b.alloc_qubits(560);
        let r = b.alloc_qubits(256);
        let s = b.alloc_qubits(256);
        for i in 0..560 {
            emit_scaled_by_controlled_microstep_for_test(&mut b, &r, &s, odd[i], a_ctrl[i], p);
        }
        let ccx = count_ccx(&b.ops);
        let peak = b.peak_qubits;
        eprintln!(
            "BY scaled controlled 560-step scaffold: ccx={ccx}, peak={peak}q, raw_control_bits=1120"
        );
        assert!(ccx < 1_160_000, "scaled 560-step scaffold cost drifted");
        assert!(peak < 2_500, "raw-control scaled scaffold exceeds current cap too much");
    }

    #[test]
    fn scaled_by_controlled_window_matches_jump_matrix() {
        const W: usize = 16;
        let p = SECP256K1_P;
        let delta = 1i64;
        let f_low = 1i128;
        let g_low = 3i128;
        let controls = branch_controls_for_lowword_window_for_test(W, delta, f_low, g_low);
        let (_, _, _, pmat) = jump_matrix_direct_lowword(W, W, delta, f_low, g_low);
        let mut b = super::super::B::new();
        let odd = b.alloc_qubits(W);
        let a_ctrl = b.alloc_qubits(W);
        let r = b.alloc_qubits(256);
        let s = b.alloc_qubits(256);
        for i in 0..W {
            emit_scaled_by_controlled_microstep_for_test(&mut b, &r, &s, odd[i], a_ctrl[i], p);
        }
        let ccx = count_ccx(&b.ops);
        let peak = b.peak_qubits;
        let num_qubits = b.next_qubit as usize;
        let num_bits = b.next_bit as usize;
        let ops = b.ops;
        let inv2w = inv_pow2_mod_p_for_test(W, p);
        let row_expected = |x0: U256, x1: U256, c0: i128, c1: i128| -> U256 {
            let t0 = mulm(signed_i128_mod_p(c0, p), x0, p);
            let t1 = mulm(signed_i128_mod_p(c1, p), x1, p);
            mulm(addm(t0, t1, p), inv2w, p)
        };
        let mut sx = Sampler::new(b"by-scaled-window-r-v1", p);
        let mut sy = Sampler::new(b"by-scaled-window-s-v1", p);
        for _ in 0..16 {
            let rv = sx.next();
            let sv = sy.next();
            let exp_r = row_expected(rv, sv, pmat.m00, pmat.m01);
            let exp_s = row_expected(rv, sv, pmat.m10, pmat.m11);
            let mut hasher = sha3::Shake128::default();
            hasher.update(b"by-scaled-window-sim-v1");
            let mut xof = hasher.finalize_xof();
            let mut sim = crate::sim::Simulator::new(num_qubits, num_bits, &mut xof);
            for (i, &(odd_v, a_v)) in controls.iter().enumerate() {
                if odd_v {
                    *sim.qubit_mut(odd[i]) |= 1;
                }
                if a_v {
                    *sim.qubit_mut(a_ctrl[i]) |= 1;
                }
            }
            set_slice_u512_by(&mut sim, &r, u256_to_u512_for_by_tests(rv));
            set_slice_u512_by(&mut sim, &s, u256_to_u512_for_by_tests(sv));
            sim.apply(&ops);
            assert_eq!(get_slice_u512_by(&sim, &r), u256_to_u512_for_by_tests(exp_r), "r mismatch");
            assert_eq!(get_slice_u512_by(&sim, &s), u256_to_u512_for_by_tests(exp_s), "s mismatch");
        }
        eprintln!(
            "BY scaled controlled 16-step window: ccx={ccx}, peak={peak}q, matrix={pmat:?}"
        );
        assert!(ccx < 35_000, "scaled controlled window too costly");
        assert!(peak < 1_350, "scaled controlled window peak too high");
    }

    #[derive(Clone, Copy, Debug)]
    struct SignedWideForTest {
        neg: bool,
        mag: U512,
    }

    fn sw_zero_for_test() -> SignedWideForTest {
        SignedWideForTest { neg: false, mag: U512::ZERO }
    }

    fn sw_from_u512_for_test(x: U512) -> SignedWideForTest {
        SignedWideForTest { neg: false, mag: x }
    }

    fn sw_add_for_test(a: SignedWideForTest, b: SignedWideForTest) -> SignedWideForTest {
        if a.mag.is_zero() { return b; }
        if b.mag.is_zero() { return a; }
        if a.neg == b.neg {
            SignedWideForTest { neg: a.neg, mag: a.mag.wrapping_add(b.mag) }
        } else if a.mag >= b.mag {
            SignedWideForTest { neg: a.neg, mag: a.mag.wrapping_sub(b.mag) }
        } else {
            SignedWideForTest { neg: b.neg, mag: b.mag.wrapping_sub(a.mag) }
        }
    }

    fn sw_sub_for_test(a: SignedWideForTest, b: SignedWideForTest) -> SignedWideForTest {
        sw_add_for_test(a, SignedWideForTest { neg: !b.neg && !b.mag.is_zero(), mag: b.mag })
    }

    fn sw_half_modp_no_reduce_for_test(mut x: SignedWideForTest, p512: U512) -> SignedWideForTest {
        if x.mag.bit(0) {
            x = sw_add_for_test(x, sw_from_u512_for_test(p512));
        }
        assert!(!x.mag.bit(0), "wide representative not even before halve");
        SignedWideForTest { neg: x.neg, mag: x.mag >> 1usize }
    }

    fn low_u256_from_u512_for_test(x: U512) -> U256 {
        let limbs = x.as_limbs();
        U256::from_limbs([limbs[0], limbs[1], limbs[2], limbs[3]])
    }

    fn sw_mod_p_for_test(x: SignedWideForTest, p: U256) -> U256 {
        let p512 = u256_to_u512_for_by_tests(p);
        let mut r = x.mag;
        while r >= p512 {
            r = r.wrapping_sub(p512);
        }
        if x.neg && !r.is_zero() {
            low_u256_from_u512_for_test(p512.wrapping_sub(r))
        } else {
            low_u256_from_u512_for_test(r)
        }
    }

    fn sw_centered_from_u256_for_test(x: U256, p: U256) -> SignedWideForTest {
        if x > (p >> 1usize) {
            SignedWideForTest { neg: true, mag: u256_to_u512_for_by_tests(p.wrapping_sub(x)) }
        } else {
            sw_from_u512_for_test(u256_to_u512_for_by_tests(x))
        }
    }

    fn sw_half_modp_centered_for_test(mut x: SignedWideForTest, p512: U512) -> SignedWideForTest {
        if x.mag.bit(0) {
            if x.neg {
                x = sw_add_for_test(x, sw_from_u512_for_test(p512));
            } else {
                x = sw_sub_for_test(x, sw_from_u512_for_test(p512));
            }
        }
        assert!(!x.mag.bit(0), "centered representative not even before halve");
        SignedWideForTest { neg: x.neg, mag: x.mag >> 1usize }
    }

    #[test]
    fn centered_signed_redundant_replay_stays_within_half_modulus() {
        // Better representative discipline for the no-red-flag replay: keep
        // signed values centered and, on an odd pre-halve representative, add
        // or subtract p according to the sign.  This keeps the whole tagged
        // channel inside about ±p/2 in all sampled trajectories.  It does not
        // by itself solve reversibility (the parity branch is still needed),
        // but it makes a narrow signed circuit and range-based parity cleanup
        // plausible enough to pursue.
        let p = SECP256K1_P;
        let p512 = u256_to_u512_for_by_tests(p);
        let samples = 2_000usize;
        let mut sx = Sampler::new(b"by-centered-signed-x-v1", p);
        let mut sy = Sampler::new(b"by-centered-signed-y-v1", p);
        let mut failures = 0usize;
        let mut max_mag = U512::ZERO;
        let mut parity_true_total = 0usize;
        for _ in 0..samples {
            let x = sx.next();
            let y = sy.next();
            let mut delta = 1i64;
            let mut f = SInt::from_u(p);
            let mut g = SInt::from_u(x);
            let mut r = sw_zero_for_test();
            let mut s = sw_centered_from_u256_for_test(addm(y, x, p), p);
            for _ in 0..560 {
                let odd = g.bit0();
                let a = delta > 0 && odd;
                if a {
                    let nr = s;
                    let t = sw_sub_for_test(s, r);
                    if t.mag.bit(0) { parity_true_total += 1; }
                    s = sw_half_modp_centered_for_test(t, p512);
                    r = nr;
                } else if odd {
                    let t = sw_add_for_test(s, r);
                    if t.mag.bit(0) { parity_true_total += 1; }
                    s = sw_half_modp_centered_for_test(t, p512);
                } else {
                    if s.mag.bit(0) { parity_true_total += 1; }
                    s = sw_half_modp_centered_for_test(s, p512);
                }
                max_mag = max_mag.max(r.mag).max(s.mag);
                divstep_sint_state(&mut delta, &mut f, &mut g);
            }
            if !g.is_zero() || !(f.is_one_pos() || f.is_one_neg()) {
                failures += 1;
                continue;
            }
            assert_eq!(sw_mod_p_for_test(s, p), U256::ZERO, "centered bottom channel not zero mod p");
            let r_mod = sw_mod_p_for_test(r, p);
            let plus_one = if f.is_one_pos() { r_mod } else { negm(r_mod, p) };
            let quotient = subm(plus_one, U256::from(1u64), p);
            assert_eq!(quotient, mulm(y, fermat_modinv(x, p), p), "centered signed quotient mismatch");
        }
        let fail_rate = failures as f64 / samples as f64;
        let parity_mean = parity_true_total as f64 / samples as f64;
        eprintln!(
            "BY centered signed replay: samples={samples}, failures={failures} ({fail_rate:.4}), max_mag_bits={}, parity_mean={parity_mean:.1}"
            , 512 - max_mag.leading_zeros()
        );
        assert!(fail_rate <= 0.01, "centered signed replay convergence tail too high");
        assert!(max_mag < p512, "centered representatives escaped one modulus");
    }

    #[test]
    fn redundant_signed_scaled_by_replay_avoids_reduction_flags_algebraically() {
        // Moonshot representation rewrite: skip modular reduction before the
        // per-step halve.  For any integer representative T of the intended
        // field value, (T + (T&1)*p)/2 is a valid representative of T/2 mod p.
        // This deletes the modular-add reduction flag algebraically.  The price
        // is noncanonical signed representatives and a parity-cleaning problem
        // for a circuit implementation.
        let p = SECP256K1_P;
        let p512 = u256_to_u512_for_by_tests(p);
        let samples = 2_000usize;
        let mut sx = Sampler::new(b"by-redundant-signed-x-v1", p);
        let mut sy = Sampler::new(b"by-redundant-signed-y-v1", p);
        let mut max_multiple = 0usize;
        let mut failures = 0usize;
        let mut parity_true_total = 0usize;
        let mut parity_step_counts = vec![0usize; 560];
        for _ in 0..samples {
            let x = sx.next();
            let y = sy.next();
            let mut delta = 1i64;
            let mut f = SInt::from_u(p);
            let mut g = SInt::from_u(x);
            let mut r = sw_zero_for_test();
            let mut s = sw_from_u512_for_test(u256_to_u512_for_by_tests(addm(y, x, p)));
            for step in 0..560 {
                let odd = g.bit0();
                let a = delta > 0 && odd;
                if a {
                    let nr = s;
                    let t = sw_sub_for_test(s, r);
                    if t.mag.bit(0) { parity_true_total += 1; parity_step_counts[step] += 1; }
                    s = sw_half_modp_no_reduce_for_test(t, p512);
                    r = nr;
                } else if odd {
                    let t = sw_add_for_test(s, r);
                    if t.mag.bit(0) { parity_true_total += 1; parity_step_counts[step] += 1; }
                    s = sw_half_modp_no_reduce_for_test(t, p512);
                } else {
                    if s.mag.bit(0) { parity_true_total += 1; parity_step_counts[step] += 1; }
                    s = sw_half_modp_no_reduce_for_test(s, p512);
                }
                for k in 0..8 {
                    if r.mag < p512 * U512::from((k + 1) as u64) && s.mag < p512 * U512::from((k + 1) as u64) {
                        max_multiple = max_multiple.max(k + 1);
                        break;
                    }
                }
                divstep_sint_state(&mut delta, &mut f, &mut g);
            }
            if !g.is_zero() || !(f.is_one_pos() || f.is_one_neg()) {
                failures += 1;
                continue;
            }
            assert_eq!(sw_mod_p_for_test(s, p), U256::ZERO, "bottom tagged channel not zero mod p");
            let r_mod = sw_mod_p_for_test(r, p);
            let plus_one = if f.is_one_pos() { r_mod } else { negm(r_mod, p) };
            let quotient = subm(plus_one, U256::from(1u64), p);
            assert_eq!(quotient, mulm(y, fermat_modinv(x, p), p), "redundant signed quotient mismatch");
        }
        let fail_rate = failures as f64 / samples as f64;
        let parity_mean = parity_true_total as f64 / samples as f64;
        let parity_entropy: f64 = parity_step_counts.iter().map(|&c| {
            let q = c as f64 / samples as f64;
            if q <= 0.0 || q >= 1.0 { 0.0 } else { -q * q.log2() - (1.0 - q) * (1.0 - q).log2() }
        }).sum();
        eprintln!(
            "BY redundant signed replay: samples={samples}, failures={failures} ({fail_rate:.4}), max_seen_multiple<={max_multiple}p, parity_mean={parity_mean:.1}, parity_entropy≈{parity_entropy:.1} bits"
        );
        assert!(fail_rate <= 0.01, "redundant signed replay convergence tail too high");
        assert!(max_multiple <= 4, "redundant representatives grow too large for a narrow circuit rewrite");
    }

    fn sw_twos_for_width_for_test(x: SignedWideForTest, width: usize) -> U512 {
        if x.mag.is_zero() {
            U512::ZERO
        } else if x.neg {
            (U512::from(1u64) << width).wrapping_sub(x.mag)
        } else {
            x.mag
        }
    }

    #[test]
    fn redundant_signed_microstep_is_cheap_if_parity_history_can_be_cleaned() {
        // Circuit-level probe for the redundant signed replay representation.
        // It has no modular-reduction comparator and therefore no red flag; it
        // leaves only the pre-halve parity bit. This is not a clean primitive
        // yet, but it quantifies the payoff if parity can be cleaned by a range
        // discipline or fused window schedule.
        const WIDE: usize = 260;
        let p = SECP256K1_P;
        let p512 = u256_to_u512_for_by_tests(p);
        let mut b = super::super::B::new();
        let odd = b.alloc_qubit();
        let a_ctrl = b.alloc_qubit();
        let parity = b.alloc_qubit();
        let r = b.alloc_qubits(WIDE);
        let s = b.alloc_qubits(WIDE);
        emit_scaled_by_redundant_signed_microstep_live_parity_for_test(&mut b, &r, &s, odd, a_ctrl, parity, p);
        let ccx = count_ccx(&b.ops);
        let peak = b.peak_qubits;
        let num_qubits = b.next_qubit as usize;
        let num_bits = b.next_bit as usize;
        let ops = b.ops;
        let cases = [(false, false, "C"), (true, false, "B"), (true, true, "A")];
        let mut sx = Sampler::new(b"by-redundant-step-r-v1", p);
        let mut sy = Sampler::new(b"by-redundant-step-s-v1", p);
        for &(odd_v, a_v, name) in &cases {
            for _ in 0..12 {
                let rv = sw_from_u512_for_test(u256_to_u512_for_by_tests(sx.next()));
                let sv = sw_from_u512_for_test(u256_to_u512_for_by_tests(sy.next()));
                let (exp_r, exp_s) = match name {
                    "A" => (sv, sw_half_modp_no_reduce_for_test(sw_sub_for_test(sv, rv), p512)),
                    "B" => (rv, sw_half_modp_no_reduce_for_test(sw_add_for_test(sv, rv), p512)),
                    "C" => (rv, sw_half_modp_no_reduce_for_test(sv, p512)),
                    _ => unreachable!(),
                };
                let mut hasher = sha3::Shake128::default();
                hasher.update(b"by-redundant-step-sim-v1");
                let mut xof = hasher.finalize_xof();
                let mut sim = crate::sim::Simulator::new(num_qubits, num_bits, &mut xof);
                if odd_v { *sim.qubit_mut(odd) |= 1; }
                if a_v { *sim.qubit_mut(a_ctrl) |= 1; }
                set_slice_u512_by(&mut sim, &r, sw_twos_for_width_for_test(rv, WIDE));
                set_slice_u512_by(&mut sim, &s, sw_twos_for_width_for_test(sv, WIDE));
                sim.apply(&ops);
                assert_eq!(get_slice_u512_by(&sim, &r), sw_twos_for_width_for_test(exp_r, WIDE), "r mismatch {name}");
                assert_eq!(get_slice_u512_by(&sim, &s), sw_twos_for_width_for_test(exp_s, WIDE), "s mismatch {name}");
            }
        }
        let replay560 = ccx * 560;
        eprintln!(
            "BY redundant signed live-parity microstep: ccx={ccx}, replay560≈{replay560}, peak={peak}q, width={WIDE}"
        );
        assert!(ccx < 1_400, "redundant signed microstep lost the no-reduction payoff");
        assert!(replay560 < 800_000, "redundant signed replay would not beat fixed-control target even if cleaned");
    }

    #[test]
    fn centered_parity_is_recoverable_from_poststate_range_for_add_cases() {
        // Small exact model of the centered redundant step.  Surprise/good
        // news: with both input registers promised centered, the pre-halve
        // parity is recoverable from poststate+case by testing whether the
        // even preimage is centered.  For B, parity=0 iff (2*s_out - r_out) is
        // centered; for A, parity=0 iff (r_out - 2*s_out) is centered.  C has
        // the analogous test on 2*s_out.  This turns parity cleanup from an
        // information-theoretic blocker into a range-test synthesis problem.
        use std::collections::BTreeMap;
        let p = 31i64;
        let centered = |x: i64| (-p / 2..=p / 2).contains(&x);
        let vals: Vec<i64> = (-(p / 2)..=(p / 2)).collect();
        let centered_halve = |mut t: i64| -> (i64, bool) {
            let parity = (t & 1) != 0;
            if parity {
                if t < 0 { t += p; } else { t -= p; }
            }
            assert_eq!(t & 1, 0);
            (t / 2, parity)
        };
        let mut b_map: BTreeMap<(i64, i64), bool> = BTreeMap::new();
        let mut a_map: BTreeMap<(i64, i64), bool> = BTreeMap::new();
        let mut c_map: BTreeMap<i64, bool> = BTreeMap::new();
        let mut mismatches = 0usize;
        for &r_old in &vals {
            for &s_old in &vals {
                let (c_s, c_par) = centered_halve(s_old);
                let c_rec = !centered(2 * c_s);
                mismatches += (c_rec != c_par) as usize;
                c_map.insert(c_s, c_par);

                let (b_s, b_par) = centered_halve(s_old + r_old);
                let b_rec = !centered(2 * b_s - r_old);
                mismatches += (b_rec != b_par) as usize;
                if let Some(prev) = b_map.insert((r_old, b_s), b_par) {
                    assert_eq!(prev, b_par, "B parity collision for r={r_old}, s_out={b_s}");
                }

                let (a_s, a_par) = centered_halve(s_old - r_old);
                let a_rec = !centered(s_old - 2 * a_s);
                mismatches += (a_rec != a_par) as usize;
                if let Some(prev) = a_map.insert((s_old, a_s), a_par) {
                    assert_eq!(prev, a_par, "A parity collision for r_out={s_old}, s_out={a_s}");
                }
            }
        }
        eprintln!(
            "BY centered parity recovery by range: states={}, B_keys={}, A_keys={}, C_keys={}, mismatches={mismatches}",
            vals.len() * vals.len(), b_map.len(), a_map.len(), c_map.len()
        );
        assert_eq!(mismatches, 0, "centered range parity recovery formula failed");
    }

    fn sw_highbits_outside_center_for_test(x: SignedWideForTest) -> bool {
        x.mag.bit(255)
    }

    #[test]
    fn centered_parity_highbits_recovery_is_too_approximate_without_boundary_fix() {
        // Tempting shortcut: exact centered recovery tests |E|>p/2, and p/2 is
        // close to 2^255, so maybe testing only signed-255 overflow is enough.
        // It is not: actual centered BY states hit the boundary band often
        // enough to exceed the user's 1% approximate-failure tolerance. Any
        // cheap range recovery needs an exact/special boundary correction, not
        // just high bits.
        let p = SECP256K1_P;
        let p512 = u256_to_u512_for_by_tests(p);
        let samples = 2_000usize;
        let mut sx = Sampler::new(b"by-centered-highbits-parity-x-v1", p);
        let mut sy = Sampler::new(b"by-centered-highbits-parity-y-v1", p);
        let mut mismatches = 0usize;
        let mut checks = 0usize;
        for _ in 0..samples {
            let x = sx.next();
            let y = sy.next();
            let mut delta = 1i64;
            let mut f = SInt::from_u(p);
            let mut g = SInt::from_u(x);
            let mut r = sw_zero_for_test();
            let mut s = sw_centered_from_u256_for_test(addm(y, x, p), p);
            for _ in 0..560 {
                let odd = g.bit0();
                let a = delta > 0 && odd;
                if a {
                    let old_r = r;
                    let old_s = s;
                    let t = sw_sub_for_test(old_s, old_r);
                    let par = t.mag.bit(0);
                    let ns = sw_half_modp_centered_for_test(t, p512);
                    let even_preimage = sw_sub_for_test(old_s, SignedWideForTest { neg: ns.neg, mag: ns.mag << 1usize });
                    let rec = sw_highbits_outside_center_for_test(even_preimage);
                    if rec != par { mismatches += 1; }
                    r = old_s;
                    s = ns;
                } else if odd {
                    let t = sw_add_for_test(s, r);
                    let par = t.mag.bit(0);
                    let ns = sw_half_modp_centered_for_test(t, p512);
                    let even_preimage = sw_sub_for_test(SignedWideForTest { neg: ns.neg, mag: ns.mag << 1usize }, r);
                    let rec = sw_highbits_outside_center_for_test(even_preimage);
                    if rec != par { mismatches += 1; }
                    s = ns;
                } else {
                    let par = s.mag.bit(0);
                    let ns = sw_half_modp_centered_for_test(s, p512);
                    let even_preimage = SignedWideForTest { neg: ns.neg, mag: ns.mag << 1usize };
                    let rec = sw_highbits_outside_center_for_test(even_preimage);
                    if rec != par { mismatches += 1; }
                    s = ns;
                }
                checks += 1;
                divstep_sint_state(&mut delta, &mut f, &mut g);
            }
        }
        let fail_rate = mismatches as f64 / checks as f64;
        eprintln!(
            "BY centered parity highbits recovery dead end: mismatches={mismatches}/{checks} ({fail_rate:.6})"
        );
        assert!(fail_rate > 0.01, "high-bit shortcut might satisfy approximate tolerance; revisit");
    }

    fn emit_highbits_outside_center_into_for_test(
        b: &mut super::super::B,
        v: &[super::super::QubitId],
        flag: super::super::QubitId,
    ) {
        let sign = v[v.len() - 1];
        let scratch = b.alloc_qubits(v.len() - 1 - 255);
        for (j, i) in (255..v.len() - 1).enumerate() {
            b.cx(v[i], scratch[j]);
            b.cx(sign, scratch[j]);
        }
        super::super::cmp_neq_zero_into(b, &scratch, flag);
        for (j, i) in (255..v.len() - 1).enumerate().rev() {
            b.cx(sign, scratch[j]);
            b.cx(v[i], scratch[j]);
        }
        b.free_vec(&scratch);
    }

    fn emit_centered_b_parity_highbits_recovery_for_cost(
        b: &mut super::super::B,
        r_out: &[super::super::QubitId],
        s_out: &[super::super::QubitId],
        flag: super::super::QubitId,
    ) {
        let n = s_out.len();
        let tmp = b.alloc_qubits(n);
        for i in 0..n - 1 { b.cx(s_out[i], tmp[i + 1]); }
        super::super::sub_nbit_qq_fast(b, r_out, &tmp);
        emit_highbits_outside_center_into_for_test(b, &tmp, flag);
        super::super::add_nbit_qq_fast(b, r_out, &tmp);
        for i in (0..n - 1).rev() { b.cx(s_out[i], tmp[i + 1]); }
        b.free_vec(&tmp);
    }

    #[test]
    fn highbits_centered_parity_recovery_cost_is_plausible_if_folded() {
        const WIDE: usize = 260;
        let mut b = super::super::B::new();
        let r = b.alloc_qubits(WIDE);
        let s = b.alloc_qubits(WIDE);
        let flag = b.alloc_qubit();
        let start = b.ops.len();
        emit_centered_b_parity_highbits_recovery_for_cost(&mut b, &r, &s, flag);
        let ccx = count_ccx(&b.ops[start..]);
        let cleanup560 = ccx as f64 * 560.0;
        eprintln!(
            "BY highbits centered parity recovery: ccx_per_Bflag={ccx}, cleanup560_worst≈{cleanup560:.0}, peak={}q",
            b.peak_qubits
        );
        assert!(ccx < 600, "high-bit parity recovery no better than full range comparator");
    }

    fn emit_centered_b_parity_recovery_for_cost(
        b: &mut super::super::B,
        r_out: &[super::super::QubitId],
        s_out: &[super::super::QubitId],
        flag: super::super::QubitId,
        p: U256,
    ) {
        // Recover B-case parity by testing whether the even preimage
        // s_old = 2*s_out - r_out is centered.  This intentionally ignores
        // shift gates in the cost (they are Clifford/swaps); it includes the
        // full-width add/sub and comparator stack a naive cleanup would pay.
        let n = s_out.len();
        let tmp = b.alloc_qubits(n);
        for i in 0..n { b.cx(s_out[i], tmp[i]); }
        // Conceptual tmp <- 2*tmp by wire permutation (0 Toffoli).
        super::super::sub_nbit_qq_fast(b, r_out, &tmp);
        let bias = p >> 1usize;
        super::super::add_nbit_const_fast(b, &tmp, bias);
        let p_reg = super::super::load_const(b, n, p);
        super::super::cmp_lt_into_fast(b, &tmp, &p_reg, flag); // flag ^= centered
        super::super::unload_const(b, &p_reg, p);
        super::super::sub_nbit_const_fast(b, &tmp, bias);
        super::super::add_nbit_qq_fast(b, r_out, &tmp);
        for i in 0..n { b.cx(s_out[i], tmp[i]); }
        b.free_vec(&tmp);
        b.x(flag); // parity = !centered(even preimage)
    }

    #[test]
    fn naive_centered_parity_recovery_cost_would_erase_redundant_replay_win() {
        const WIDE: usize = 260;
        let p = SECP256K1_P;
        let mut b = super::super::B::new();
        let r = b.alloc_qubits(WIDE);
        let s = b.alloc_qubits(WIDE);
        let flag = b.alloc_qubit();
        let start = b.ops.len();
        emit_centered_b_parity_recovery_for_cost(&mut b, &r, &s, flag, p);
        let ccx = count_ccx(&b.ops[start..]);
        let replay_cleanup = ccx as f64 * 560.0;
        eprintln!(
            "BY naive centered parity recovery: ccx_per_flag={ccx}, cleanup560≈{replay_cleanup:.0}, peak={}q",
            b.peak_qubits
        );
        assert!(ccx > 900, "naive range recovery unexpectedly cheap; synthesize it for real");
        assert!(replay_cleanup > 500_000.0, "naive parity cleanup would preserve redundant replay win");
    }

    #[test]
    fn centered_signed_microstep_keeps_narrow_reps_at_submillion_cost() {
        const WIDE: usize = 260;
        let p = SECP256K1_P;
        let p512 = u256_to_u512_for_by_tests(p);
        let mut b = super::super::B::new();
        let odd = b.alloc_qubit();
        let a_ctrl = b.alloc_qubit();
        let parity = b.alloc_qubit();
        let r = b.alloc_qubits(WIDE);
        let s = b.alloc_qubits(WIDE);
        emit_scaled_by_centered_signed_microstep_live_parity_for_test(&mut b, &r, &s, odd, a_ctrl, parity, p);
        let ccx = count_ccx(&b.ops);
        let peak = b.peak_qubits;
        let num_qubits = b.next_qubit as usize;
        let num_bits = b.next_bit as usize;
        let ops = b.ops;
        let cases = [(false, false, "C"), (true, false, "B"), (true, true, "A")];
        let mut sx = Sampler::new(b"by-centered-step-r-v1", p);
        let mut sy = Sampler::new(b"by-centered-step-s-v1", p);
        for &(odd_v, a_v, name) in &cases {
            for _ in 0..12 {
                let rv = sw_centered_from_u256_for_test(sx.next(), p);
                let sv = sw_centered_from_u256_for_test(sy.next(), p);
                let (exp_r, exp_s) = match name {
                    "A" => (sv, sw_half_modp_centered_for_test(sw_sub_for_test(sv, rv), p512)),
                    "B" => (rv, sw_half_modp_centered_for_test(sw_add_for_test(sv, rv), p512)),
                    "C" => (rv, sw_half_modp_centered_for_test(sv, p512)),
                    _ => unreachable!(),
                };
                let mut hasher = sha3::Shake128::default();
                hasher.update(b"by-centered-step-sim-v1");
                let mut xof = hasher.finalize_xof();
                let mut sim = crate::sim::Simulator::new(num_qubits, num_bits, &mut xof);
                if odd_v { *sim.qubit_mut(odd) |= 1; }
                if a_v { *sim.qubit_mut(a_ctrl) |= 1; }
                set_slice_u512_by(&mut sim, &r, sw_twos_for_width_for_test(rv, WIDE));
                set_slice_u512_by(&mut sim, &s, sw_twos_for_width_for_test(sv, WIDE));
                sim.apply(&ops);
                assert_eq!(get_slice_u512_by(&sim, &r), sw_twos_for_width_for_test(exp_r, WIDE), "r mismatch {name}");
                assert_eq!(get_slice_u512_by(&sim, &s), sw_twos_for_width_for_test(exp_s, WIDE), "s mismatch {name}");
            }
        }
        let replay560 = ccx * 560;
        eprintln!(
            "BY centered signed live-parity microstep: ccx={ccx}, replay560≈{replay560}, peak={peak}q, width={WIDE}"
        );
        assert!(ccx < 1_700, "centered signed microstep too costly to stay SOTA-shaped");
        assert!(replay560 < 1_000_000, "centered signed replay loses the sub-1M target");
    }

    #[test]
    fn centered_signed_560_scaffold_hits_submillion_replay_with_live_parity() {
        const WIDE: usize = 260;
        let p = SECP256K1_P;
        let p512 = u256_to_u512_for_by_tests(p);
        let mut sx = Sampler::new(b"by-centered-560-x-v1", p);
        let mut sy = Sampler::new(b"by-centered-560-y-v1", p);
        let (x, y, controls, f_final, exp_r, exp_s) = loop {
            let x = sx.next();
            let y = sy.next();
            let mut delta = 1i64;
            let mut f = SInt::from_u(p);
            let mut g = SInt::from_u(x);
            let mut controls = Vec::with_capacity(560);
            let mut r = sw_zero_for_test();
            let mut s = sw_centered_from_u256_for_test(addm(y, x, p), p);
            for _ in 0..560 {
                let odd = g.bit0();
                let a = delta > 0 && odd;
                controls.push((odd, a));
                if a {
                    let nr = s;
                    s = sw_half_modp_centered_for_test(sw_sub_for_test(s, r), p512);
                    r = nr;
                } else if odd {
                    s = sw_half_modp_centered_for_test(sw_add_for_test(s, r), p512);
                } else {
                    s = sw_half_modp_centered_for_test(s, p512);
                }
                divstep_sint_state(&mut delta, &mut f, &mut g);
            }
            if g.is_zero() && (f.is_one_pos() || f.is_one_neg()) {
                break (x, y, controls, f, r, s);
            }
        };
        let mut b = super::super::B::new();
        let odd = b.alloc_qubits(560);
        let a_ctrl = b.alloc_qubits(560);
        let parity = b.alloc_qubits(560);
        let r = b.alloc_qubits(WIDE);
        let s = b.alloc_qubits(WIDE);
        for i in 0..560 {
            emit_scaled_by_centered_signed_microstep_live_parity_for_test(&mut b, &r, &s, odd[i], a_ctrl[i], parity[i], p);
        }
        let ccx = count_ccx(&b.ops);
        let peak = b.peak_qubits;
        let num_qubits = b.next_qubit as usize;
        let num_bits = b.next_bit as usize;
        let ops = b.ops;
        let mut hasher = sha3::Shake128::default();
        hasher.update(b"by-centered-560-sim-v1");
        let mut xof = hasher.finalize_xof();
        let mut sim = crate::sim::Simulator::new(num_qubits, num_bits, &mut xof);
        for (i, &(odd_v, a_v)) in controls.iter().enumerate() {
            if odd_v { *sim.qubit_mut(odd[i]) |= 1; }
            if a_v { *sim.qubit_mut(a_ctrl[i]) |= 1; }
        }
        set_slice_u512_by(&mut sim, &r, U512::ZERO);
        set_slice_u512_by(&mut sim, &s, sw_twos_for_width_for_test(sw_centered_from_u256_for_test(addm(y, x, p), p), WIDE));
        sim.apply(&ops);
        assert_eq!(get_slice_u512_by(&sim, &r), sw_twos_for_width_for_test(exp_r, WIDE), "centered r mismatch");
        assert_eq!(get_slice_u512_by(&sim, &s), sw_twos_for_width_for_test(exp_s, WIDE), "centered s mismatch");
        assert_eq!(sw_mod_p_for_test(exp_s, p), U256::ZERO, "bottom channel not zero mod p");
        let r_mod = sw_mod_p_for_test(exp_r, p);
        let plus_one = if f_final.is_one_pos() { r_mod } else { negm(r_mod, p) };
        let quotient = subm(plus_one, U256::from(1u64), p);
        assert_eq!(quotient, mulm(y, fermat_modinv(x, p), p), "centered scaffold quotient mismatch");
        let parity_nonzero = parity.iter().any(|&q| (sim.qubit(q) & 1) != 0);
        eprintln!(
            "BY centered signed 560 scaffold: ccx={ccx}, peak={peak}q, parity_nonzero={parity_nonzero}"
        );
        assert!(parity_nonzero, "centered parity history unexpectedly clean");
        assert!(ccx < 900_000, "centered signed scaffold lost sub-million replay target");
        assert!(peak < 2_800, "centered signed scaffold exceeds current cap");
        assert_eq!(sim.global_phase() & 1, 0, "centered signed scaffold phase garbage");
    }

    #[test]
    fn centered_signed_inverse_560_product_clean_scaffold_matches_forward_cost() {
        const WIDE: usize = 260;
        let p = SECP256K1_P;
        let p512 = u256_to_u512_for_by_tests(p);
        let mut sx = Sampler::new(b"by-centered-inv-560-x-v1", p);
        let mut sq = Sampler::new(b"by-centered-inv-560-q-v1", p);
        let (x, q, controls, parity_hist, f_final, final_r, final_s, start_s) = loop {
            let x = sx.next();
            let q = sq.next();
            let start_s = sw_centered_from_u256_for_test(mulm(q, x, p), p);
            let mut delta = 1i64;
            let mut f = SInt::from_u(p);
            let mut g = SInt::from_u(x);
            let mut controls = Vec::with_capacity(560);
            let mut parity_hist = Vec::with_capacity(560);
            let mut r = sw_zero_for_test();
            let mut s = start_s;
            for _ in 0..560 {
                let odd = g.bit0();
                let a = delta > 0 && odd;
                controls.push((odd, a));
                if a {
                    let nr = s;
                    let t = sw_sub_for_test(s, r);
                    parity_hist.push(t.mag.bit(0));
                    s = sw_half_modp_centered_for_test(t, p512);
                    r = nr;
                } else if odd {
                    let t = sw_add_for_test(s, r);
                    parity_hist.push(t.mag.bit(0));
                    s = sw_half_modp_centered_for_test(t, p512);
                } else {
                    parity_hist.push(s.mag.bit(0));
                    s = sw_half_modp_centered_for_test(s, p512);
                }
                divstep_sint_state(&mut delta, &mut f, &mut g);
            }
            if g.is_zero() && (f.is_one_pos() || f.is_one_neg()) {
                break (x, q, controls, parity_hist, f, r, s, start_s);
            }
        };
        assert_eq!(sw_mod_p_for_test(final_s, p), U256::ZERO, "forward centered final s not zero");
        let r_mod = sw_mod_p_for_test(final_r, p);
        let recovered_q = if f_final.is_one_pos() { r_mod } else { negm(r_mod, p) };
        assert_eq!(recovered_q, q, "forward centered q frame mismatch before inverse test");

        let mut b = super::super::B::new();
        let odd = b.alloc_qubits(560);
        let a_ctrl = b.alloc_qubits(560);
        let parity = b.alloc_qubits(560);
        let r = b.alloc_qubits(WIDE);
        let s = b.alloc_qubits(WIDE);
        for i in (0..560).rev() {
            emit_scaled_by_centered_signed_microstep_inverse_live_parity_for_test(&mut b, &r, &s, odd[i], a_ctrl[i], parity[i], p);
        }
        let ccx = count_ccx(&b.ops);
        let peak = b.peak_qubits;
        let num_qubits = b.next_qubit as usize;
        let num_bits = b.next_bit as usize;
        let ops = b.ops;
        let mut hasher = sha3::Shake128::default();
        hasher.update(b"by-centered-inv-560-sim-v1");
        let mut xof = hasher.finalize_xof();
        let mut sim = crate::sim::Simulator::new(num_qubits, num_bits, &mut xof);
        for (i, &(odd_v, a_v)) in controls.iter().enumerate() {
            if odd_v { *sim.qubit_mut(odd[i]) |= 1; }
            if a_v { *sim.qubit_mut(a_ctrl[i]) |= 1; }
            if parity_hist[i] { *sim.qubit_mut(parity[i]) |= 1; }
        }
        set_slice_u512_by(&mut sim, &r, sw_twos_for_width_for_test(final_r, WIDE));
        set_slice_u512_by(&mut sim, &s, sw_twos_for_width_for_test(final_s, WIDE));
        sim.apply(&ops);
        assert_eq!(get_slice_u512_by(&sim, &r), U512::ZERO, "inverse centered r not zero");
        assert_eq!(get_slice_u512_by(&sim, &s), sw_twos_for_width_for_test(start_s, WIDE), "inverse centered product mismatch");
        assert_eq!(sw_mod_p_for_test(start_s, p), mulm(q, x, p), "start product frame mismatch");
        eprintln!(
            "BY centered signed inverse 560 scaffold: ccx={ccx}, peak={peak}q"
        );
        assert!(ccx < 900_000, "centered inverse scaffold lost sub-million target");
        assert!(peak < 2_800, "centered inverse scaffold exceeds current cap");
        assert_eq!(sim.global_phase() & 1, 0, "centered inverse scaffold phase garbage");
    }

    #[test]
    fn centered_signed_roundtrip_parity_clear_is_phase_clean_with_exact_controls() {
        const WIDE: usize = 260;
        const STEPS: usize = 96;
        let p = SECP256K1_P;
        let p512 = u256_to_u512_for_by_tests(p);
        let mut sx = Sampler::new(b"by-centered-roundtrip-exact-x-v1", p);
        let mut sy = Sampler::new(b"by-centered-roundtrip-exact-y-v1", p);
        let x = sx.next();
        let y = sy.next();
        let start_s = sw_centered_from_u256_for_test(addm(y, x, p), p);
        let mut delta = 1i64;
        let mut f = SInt::from_u(p);
        let mut g = SInt::from_u(x);
        let mut controls = Vec::with_capacity(STEPS);
        let mut class_r = sw_zero_for_test();
        let mut class_s = start_s;
        for _ in 0..STEPS {
            let odd = g.bit0();
            let a = delta > 0 && odd;
            controls.push((odd, a));
            if a {
                let nr = class_s;
                class_s = sw_half_modp_centered_for_test(sw_sub_for_test(class_s, class_r), p512);
                class_r = nr;
            } else if odd {
                class_s = sw_half_modp_centered_for_test(sw_add_for_test(class_s, class_r), p512);
            } else {
                class_s = sw_half_modp_centered_for_test(class_s, p512);
            }
            divstep_sint_state(&mut delta, &mut f, &mut g);
        }

        let variants = [
            ("fast_mbu", false, false),
            ("exact_signed_only", true, false),
            ("exact_parity_controls", false, true),
            ("all_exact", true, true),
        ];
        let mut phases = Vec::new();
        for &(name, exact_signed_add, exact_parity_const) in &variants {
            let mut b = super::super::B::new();
            let odd = b.alloc_qubits(STEPS);
            let a_ctrl = b.alloc_qubits(STEPS);
            let parity = b.alloc_qubits(STEPS);
            let r = b.alloc_qubits(WIDE);
            let s = b.alloc_qubits(WIDE);
            for i in 0..STEPS {
                emit_scaled_by_centered_signed_microstep_live_parity_variant_for_test(
                    &mut b,
                    &r,
                    &s,
                    odd[i],
                    a_ctrl[i],
                    parity[i],
                    p,
                    exact_signed_add,
                    exact_parity_const,
                );
            }
            for i in (0..STEPS).rev() {
                emit_scaled_by_centered_signed_microstep_inverse_live_parity_variant_for_test(
                    &mut b,
                    &r,
                    &s,
                    odd[i],
                    a_ctrl[i],
                    parity[i],
                    p,
                    exact_signed_add,
                    exact_parity_const,
                );
                emit_centered_signed_clear_parity_after_inverse_for_test(&mut b, &r, &s, odd[i], parity[i]);
            }
            let ccx = count_ccx(&b.ops);
            let peak = b.peak_qubits;
            let num_qubits = b.next_qubit as usize;
            let num_bits = b.next_bit as usize;
            let ops = b.ops;
            let mut hasher = sha3::Shake128::default();
            hasher.update(b"by-centered-roundtrip-exact-sim-v1");
            hasher.update(name.as_bytes());
            let mut xof = hasher.finalize_xof();
            let mut sim = crate::sim::Simulator::new(num_qubits, num_bits, &mut xof);
            for (i, &(odd_v, a_v)) in controls.iter().enumerate() {
                if odd_v { *sim.qubit_mut(odd[i]) |= 1; }
                if a_v { *sim.qubit_mut(a_ctrl[i]) |= 1; }
            }
            set_slice_u512_by(&mut sim, &r, U512::ZERO);
            set_slice_u512_by(&mut sim, &s, sw_twos_for_width_for_test(start_s, WIDE));
            sim.apply(&ops);
            assert_eq!(get_slice_u512_by(&sim, &r), U512::ZERO, "{name}: r not restored");
            assert_eq!(get_slice_u512_by(&sim, &s), sw_twos_for_width_for_test(start_s, WIDE), "{name}: s not restored");
            for (j, &q) in parity.iter().enumerate() {
                assert_eq!(sim.qubit(q) & 1, 0, "{name}: parity[{j}] not cleared");
            }
            let phase = sim.global_phase() & 1;
            phases.push((name, phase, ccx, peak));
        }
        for (name, phase, ccx, peak) in &phases {
            eprintln!("BY centered parity-clear roundtrip variant {name}: phase={phase}, ccx={ccx}, peak={peak}q");
        }
        let fast_phase = phases.iter().find(|(name, _, _, _)| *name == "fast_mbu").unwrap().1;
        let exact_signed_phase = phases.iter().find(|(name, _, _, _)| *name == "exact_signed_only").unwrap().1;
        let exact_parity_phase = phases.iter().find(|(name, _, _, _)| *name == "exact_parity_controls").unwrap().1;
        let all_exact_phase = phases.iter().find(|(name, _, _, _)| *name == "all_exact").unwrap().1;
        assert_eq!(fast_phase, 0, "fixed unhalve sign history should make fast centered parity clear phase clean");
        assert_eq!(exact_signed_phase, 0, "fixed unhalve sign history should also clean exact-signed variant");
        assert_eq!(exact_parity_phase, 0, "exact parity-controlled ±p corrections should remain phase clean");
        assert_eq!(all_exact_phase, 0, "all-exact centered parity-clear roundtrip should be phase clean");
    }

    #[test]
    fn centered_signed_560_parity_can_be_cleaned_phase_safely_only_with_all_exact_controls() {
        const WIDE: usize = 260;
        const STEPS: usize = 560;
        let p = SECP256K1_P;
        let mut sx = Sampler::new(b"by-centered-clean560-x-v1", p);
        let mut sy = Sampler::new(b"by-centered-clean560-y-v1", p);
        let x = sx.next();
        let y = sy.next();
        let start_s = sw_centered_from_u256_for_test(addm(y, x, p), p);
        let mut delta = 1i64;
        let mut f = SInt::from_u(p);
        let mut g = SInt::from_u(x);
        let mut controls = Vec::with_capacity(STEPS);
        for _ in 0..STEPS {
            let odd = g.bit0();
            let a = delta > 0 && odd;
            controls.push((odd, a));
            divstep_sint_state(&mut delta, &mut f, &mut g);
        }

        let mut b = super::super::B::new();
        let odd = b.alloc_qubits(STEPS);
        let a_ctrl = b.alloc_qubits(STEPS);
        let parity = b.alloc_qubits(STEPS);
        let r = b.alloc_qubits(WIDE);
        let s = b.alloc_qubits(WIDE);
        for i in 0..STEPS {
            emit_scaled_by_centered_signed_microstep_live_parity_variant_for_test(
                &mut b,
                &r,
                &s,
                odd[i],
                a_ctrl[i],
                parity[i],
                p,
                true,
                true,
            );
        }
        for i in (0..STEPS).rev() {
            emit_scaled_by_centered_signed_microstep_inverse_live_parity_variant_for_test(
                &mut b,
                &r,
                &s,
                odd[i],
                a_ctrl[i],
                parity[i],
                p,
                true,
                true,
            );
            emit_centered_signed_clear_parity_after_inverse_for_test(&mut b, &r, &s, odd[i], parity[i]);
        }
        let ccx = count_ccx(&b.ops);
        let peak = b.peak_qubits;
        let num_qubits = b.next_qubit as usize;
        let num_bits = b.next_bit as usize;
        let ops = b.ops;
        let mut hasher = sha3::Shake128::default();
        hasher.update(b"by-centered-clean560-sim-v1");
        let mut xof = hasher.finalize_xof();
        let mut sim = crate::sim::Simulator::new(num_qubits, num_bits, &mut xof);
        for (i, &(odd_v, a_v)) in controls.iter().enumerate() {
            if odd_v { *sim.qubit_mut(odd[i]) |= 1; }
            if a_v { *sim.qubit_mut(a_ctrl[i]) |= 1; }
        }
        set_slice_u512_by(&mut sim, &r, U512::ZERO);
        set_slice_u512_by(&mut sim, &s, sw_twos_for_width_for_test(start_s, WIDE));
        sim.apply(&ops);
        assert_eq!(get_slice_u512_by(&sim, &r), U512::ZERO, "clean560 r not restored");
        assert_eq!(get_slice_u512_by(&sim, &s), sw_twos_for_width_for_test(start_s, WIDE), "clean560 s not restored");
        for (j, &q) in parity.iter().enumerate() {
            assert_eq!(sim.qubit(q) & 1, 0, "clean560 parity[{j}] not cleared");
        }
        eprintln!(
            "BY centered signed clean 560 roundtrip with all exact controls: ccx={ccx}, peak={peak}q, phase={}",
            sim.global_phase() & 1
        );
        assert_eq!(sim.global_phase() & 1, 0, "clean560 all-exact-control roundtrip phase garbage");
        assert!(ccx > 3_000_000, "all-exact-control fallback unexpectedly SOTA-shaped");
        assert!(peak < 2_800, "all-exact-control clean roundtrip exceeds current cap");
    }

    #[test]
    fn centered_signed_fast_signed_phase_after_exact_parity_controls_is_clean_after_unhalve_fix() {
        const WIDE: usize = 260;
        const STEPS: usize = 560;
        let p = SECP256K1_P;
        let mut b = super::super::B::new();
        let odd = b.alloc_qubits(STEPS);
        let a_ctrl = b.alloc_qubits(STEPS);
        let parity = b.alloc_qubits(STEPS);
        let r = b.alloc_qubits(WIDE);
        let s = b.alloc_qubits(WIDE);
        for i in 0..STEPS {
            emit_scaled_by_centered_signed_microstep_live_parity_variant_for_test(
                &mut b, &r, &s, odd[i], a_ctrl[i], parity[i], p, false, true,
            );
        }
        for i in (0..STEPS).rev() {
            emit_scaled_by_centered_signed_microstep_inverse_live_parity_variant_for_test(
                &mut b, &r, &s, odd[i], a_ctrl[i], parity[i], p, false, true,
            );
            emit_centered_signed_clear_parity_after_inverse_for_test(&mut b, &r, &s, odd[i], parity[i]);
        }
        let ccx = count_ccx(&b.ops);
        let peak = b.peak_qubits;
        const NON_DIV_SCAFFOLD: usize = 642_716;
        const GOOGLE_LOW_QUBIT_TOFFOLI: usize = 2_700_000;
        let per_div_budget = (GOOGLE_LOW_QUBIT_TOFFOLI - NON_DIV_SCAFFOLD) / 2;
        let two_clean_div_projection = NON_DIV_SCAFFOLD + 2 * ccx;
        eprintln!(
            "BY centered exact-parity/fast-signed full clean budget: ccx={ccx}, peak={peak}q, per_div_budget={per_div_budget}, two_clean_div_projection={two_clean_div_projection}"
        );
        assert!(peak < 2_800, "centered exact-parity clean replay no longer fits the current qubit cap");
        assert!(
            ccx > per_div_budget,
            "exact-parity clean replay now fits a one-DIV low-qubit budget; revisit BY promotion"
        );
        assert!(
            two_clean_div_projection > GOOGLE_LOW_QUBIT_TOFFOLI,
            "two clean centered BY DIVs now fit the low-qubit target; revisit BY promotion"
        );
        let num_qubits = b.next_qubit as usize;
        let num_bits = b.next_bit as usize;
        let ops = b.ops;
        let mut sx = Sampler::new(b"by-centered-exactparity-phase-x-v1", p);
        let mut sy = Sampler::new(b"by-centered-exactparity-phase-y-v1", p);
        let mut saw = [false; 2];
        for sample in 0..12 {
            let x = sx.next();
            let y = sy.next();
            let start_s = sw_centered_from_u256_for_test(addm(y, x, p), p);
            let mut delta = 1i64;
            let mut f = SInt::from_u(p);
            let mut g = SInt::from_u(x);
            let mut controls = Vec::with_capacity(STEPS);
            for _ in 0..STEPS {
                let odd_v = g.bit0();
                let a_v = delta > 0 && odd_v;
                controls.push((odd_v, a_v));
                divstep_sint_state(&mut delta, &mut f, &mut g);
            }
            let mut hasher = sha3::Shake128::default();
            hasher.update(b"by-centered-exactparity-phase-sim-v1");
            hasher.update(&(sample as u64).to_le_bytes());
            let mut xof = hasher.finalize_xof();
            let mut sim = crate::sim::Simulator::new(num_qubits, num_bits, &mut xof);
            for (i, &(odd_v, a_v)) in controls.iter().enumerate() {
                if odd_v { *sim.qubit_mut(odd[i]) |= 1; }
                if a_v { *sim.qubit_mut(a_ctrl[i]) |= 1; }
            }
            set_slice_u512_by(&mut sim, &r, U512::ZERO);
            set_slice_u512_by(&mut sim, &s, sw_twos_for_width_for_test(start_s, WIDE));
            sim.apply(&ops);
            assert_eq!(get_slice_u512_by(&sim, &r), U512::ZERO, "sample {sample}: r not restored");
            assert_eq!(get_slice_u512_by(&sim, &s), sw_twos_for_width_for_test(start_s, WIDE), "sample {sample}: s not restored");
            for &q in &parity { assert_eq!(sim.qubit(q) & 1, 0, "sample {sample}: parity not clean"); }
            saw[(sim.global_phase() & 1) as usize] = true;
        }
        eprintln!("BY centered exact-parity/fast-signed full clean phases after unhalve fix: saw0={}, saw1={}", saw[0], saw[1]);
        assert!(saw[0] && !saw[1], "fixed unhalve sign history should remove data-dependent fast signed-control phase");
    }

    #[test]
    fn negating_signed_copy_measurements_does_not_fix_centered_control_phase() {
        const WIDE: usize = 260;
        const STEPS: usize = 560;
        let p = SECP256K1_P;
        let mut b = super::super::B::new();
        let odd = b.alloc_qubits(STEPS);
        let a_ctrl = b.alloc_qubits(STEPS);
        let parity = b.alloc_qubits(STEPS);
        let r = b.alloc_qubits(WIDE);
        let s = b.alloc_qubits(WIDE);
        for i in 0..STEPS {
            emit_scaled_by_centered_signed_microstep_live_parity_negcopy_exact_parity_for_test(
                &mut b, &r, &s, odd[i], a_ctrl[i], parity[i], p,
            );
        }
        for i in (0..STEPS).rev() {
            emit_scaled_by_centered_signed_microstep_inverse_live_parity_negcopy_exact_parity_for_test(
                &mut b, &r, &s, odd[i], a_ctrl[i], parity[i], p,
            );
            emit_centered_signed_clear_parity_after_inverse_for_test(&mut b, &r, &s, odd[i], parity[i]);
        }
        let ccx = count_ccx(&b.ops);
        let num_qubits = b.next_qubit as usize;
        let num_bits = b.next_bit as usize;
        let ops = b.ops;
        let mut sx = Sampler::new(b"by-centered-negcopy-phase-x-v1", p);
        let mut sy = Sampler::new(b"by-centered-negcopy-phase-y-v1", p);
        let mut saw = [false; 2];
        for sample in 0..12 {
            let x = sx.next();
            let y = sy.next();
            let start_s = sw_centered_from_u256_for_test(addm(y, x, p), p);
            let mut delta = 1i64;
            let mut f = SInt::from_u(p);
            let mut g = SInt::from_u(x);
            let mut controls = Vec::with_capacity(STEPS);
            for _ in 0..STEPS {
                let odd_v = g.bit0();
                let a_v = delta > 0 && odd_v;
                controls.push((odd_v, a_v));
                divstep_sint_state(&mut delta, &mut f, &mut g);
            }
            let mut hasher = sha3::Shake128::default();
            hasher.update(b"by-centered-negcopy-phase-sim-v1");
            hasher.update(&(sample as u64).to_le_bytes());
            let mut xof = hasher.finalize_xof();
            let mut sim = crate::sim::Simulator::new(num_qubits, num_bits, &mut xof);
            for (i, &(odd_v, a_v)) in controls.iter().enumerate() {
                if odd_v { *sim.qubit_mut(odd[i]) |= 1; }
                if a_v { *sim.qubit_mut(a_ctrl[i]) |= 1; }
            }
            set_slice_u512_by(&mut sim, &r, U512::ZERO);
            set_slice_u512_by(&mut sim, &s, sw_twos_for_width_for_test(start_s, WIDE));
            sim.apply(&ops);
            assert_eq!(get_slice_u512_by(&sim, &r), U512::ZERO, "sample {sample}: r not restored");
            assert_eq!(get_slice_u512_by(&sim, &s), sw_twos_for_width_for_test(start_s, WIDE), "sample {sample}: s not restored");
            saw[(sim.global_phase() & 1) as usize] = true;
        }
        eprintln!("BY centered negcopy signed-control phase test: ccx={ccx}, saw0={}, saw1={}", saw[0], saw[1]);
        assert!(saw[0] && saw[1], "neg_if on copied signed-add controls unexpectedly made phase a global constant");
    }

    #[test]
    fn negating_cuccaro_carry_measurements_does_not_fix_centered_control_phase() {
        const WIDE: usize = 260;
        const STEPS: usize = 560;
        let p = SECP256K1_P;
        let mut b = super::super::B::new();
        let odd = b.alloc_qubits(STEPS);
        let a_ctrl = b.alloc_qubits(STEPS);
        let parity = b.alloc_qubits(STEPS);
        let r = b.alloc_qubits(WIDE);
        let s = b.alloc_qubits(WIDE);
        for i in 0..STEPS {
            emit_scaled_by_centered_signed_microstep_live_parity_negcarry_exact_parity_for_test(
                &mut b, &r, &s, odd[i], a_ctrl[i], parity[i], p,
            );
        }
        for i in (0..STEPS).rev() {
            emit_scaled_by_centered_signed_microstep_inverse_live_parity_negcarry_exact_parity_for_test(
                &mut b, &r, &s, odd[i], a_ctrl[i], parity[i], p,
            );
            emit_centered_signed_clear_parity_after_inverse_for_test(&mut b, &r, &s, odd[i], parity[i]);
        }
        let ccx = count_ccx(&b.ops);
        let num_qubits = b.next_qubit as usize;
        let num_bits = b.next_bit as usize;
        let ops = b.ops;
        let mut sx = Sampler::new(b"by-centered-negcarry-phase-x-v1", p);
        let mut sy = Sampler::new(b"by-centered-negcarry-phase-y-v1", p);
        let mut saw = [false; 2];
        for sample in 0..12 {
            let x = sx.next();
            let y = sy.next();
            let start_s = sw_centered_from_u256_for_test(addm(y, x, p), p);
            let mut delta = 1i64;
            let mut f = SInt::from_u(p);
            let mut g = SInt::from_u(x);
            let mut controls = Vec::with_capacity(STEPS);
            for _ in 0..STEPS {
                let odd_v = g.bit0();
                let a_v = delta > 0 && odd_v;
                controls.push((odd_v, a_v));
                divstep_sint_state(&mut delta, &mut f, &mut g);
            }
            let mut hasher = sha3::Shake128::default();
            hasher.update(b"by-centered-negcarry-phase-sim-v1");
            hasher.update(&(sample as u64).to_le_bytes());
            let mut xof = hasher.finalize_xof();
            let mut sim = crate::sim::Simulator::new(num_qubits, num_bits, &mut xof);
            for (i, &(odd_v, a_v)) in controls.iter().enumerate() {
                if odd_v { *sim.qubit_mut(odd[i]) |= 1; }
                if a_v { *sim.qubit_mut(a_ctrl[i]) |= 1; }
            }
            set_slice_u512_by(&mut sim, &r, U512::ZERO);
            set_slice_u512_by(&mut sim, &s, sw_twos_for_width_for_test(start_s, WIDE));
            sim.apply(&ops);
            assert_eq!(get_slice_u512_by(&sim, &r), U512::ZERO, "sample {sample}: r not restored");
            assert_eq!(get_slice_u512_by(&sim, &s), sw_twos_for_width_for_test(start_s, WIDE), "sample {sample}: s not restored");
            saw[(sim.global_phase() & 1) as usize] = true;
        }
        eprintln!("BY centered negcarry signed-control phase test: ccx={ccx}, saw0={}, saw1={}", saw[0], saw[1]);
        assert!(saw[0] && saw[1], "neg_if on Cuccaro carry measurements unexpectedly fixed the cleaned-control phase");
    }

    #[test]
    fn negating_all_signed_mbu_measurements_does_not_fix_centered_control_phase() {
        const WIDE: usize = 260;
        const STEPS: usize = 560;
        let p = SECP256K1_P;
        let mut b = super::super::B::new();
        let odd = b.alloc_qubits(STEPS);
        let a_ctrl = b.alloc_qubits(STEPS);
        let parity = b.alloc_qubits(STEPS);
        let r = b.alloc_qubits(WIDE);
        let s = b.alloc_qubits(WIDE);
        for i in 0..STEPS {
            emit_scaled_by_centered_signed_microstep_live_parity_negboth_exact_parity_for_test(
                &mut b, &r, &s, odd[i], a_ctrl[i], parity[i], p,
            );
        }
        for i in (0..STEPS).rev() {
            emit_scaled_by_centered_signed_microstep_inverse_live_parity_negboth_exact_parity_for_test(
                &mut b, &r, &s, odd[i], a_ctrl[i], parity[i], p,
            );
            emit_centered_signed_clear_parity_after_inverse_for_test(&mut b, &r, &s, odd[i], parity[i]);
        }
        let ccx = count_ccx(&b.ops);
        let num_qubits = b.next_qubit as usize;
        let num_bits = b.next_bit as usize;
        let ops = b.ops;
        let mut sx = Sampler::new(b"by-centered-negboth-phase-x-v1", p);
        let mut sy = Sampler::new(b"by-centered-negboth-phase-y-v1", p);
        let mut saw = [false; 2];
        for sample in 0..12 {
            let x = sx.next();
            let y = sy.next();
            let start_s = sw_centered_from_u256_for_test(addm(y, x, p), p);
            let mut delta = 1i64;
            let mut f = SInt::from_u(p);
            let mut g = SInt::from_u(x);
            let mut controls = Vec::with_capacity(STEPS);
            for _ in 0..STEPS {
                let odd_v = g.bit0();
                let a_v = delta > 0 && odd_v;
                controls.push((odd_v, a_v));
                divstep_sint_state(&mut delta, &mut f, &mut g);
            }
            let mut hasher = sha3::Shake128::default();
            hasher.update(b"by-centered-negboth-phase-sim-v1");
            hasher.update(&(sample as u64).to_le_bytes());
            let mut xof = hasher.finalize_xof();
            let mut sim = crate::sim::Simulator::new(num_qubits, num_bits, &mut xof);
            for (i, &(odd_v, a_v)) in controls.iter().enumerate() {
                if odd_v { *sim.qubit_mut(odd[i]) |= 1; }
                if a_v { *sim.qubit_mut(a_ctrl[i]) |= 1; }
            }
            set_slice_u512_by(&mut sim, &r, U512::ZERO);
            set_slice_u512_by(&mut sim, &s, sw_twos_for_width_for_test(start_s, WIDE));
            sim.apply(&ops);
            assert_eq!(get_slice_u512_by(&sim, &r), U512::ZERO, "sample {sample}: r not restored");
            assert_eq!(get_slice_u512_by(&sim, &s), sw_twos_for_width_for_test(start_s, WIDE), "sample {sample}: s not restored");
            saw[(sim.global_phase() & 1) as usize] = true;
        }
        eprintln!("BY centered negboth signed-control phase test: ccx={ccx}, saw0={}, saw1={}", saw[0], saw[1]);
        assert!(saw[0] && saw[1], "neg_if on all signed MBU measurements unexpectedly fixed the cleaned-control phase");
    }

    #[test]
    fn centered_clean_roundtrip_fixed_trace_for_benchmark_hook_is_phase_clean() {
        const WIDE: usize = 260;
        const STEPS: usize = 560;
        const ODD_WORDS: [u64; 9] = [
            0x9f0102a4a879b9a7,
            0x39950f607ecb1db3,
            0xefaf7e99e64fb43a,
            0x6f3857abf7ed1f44,
            0x5b90e29f6d3d3b0c,
            0xb9f3f86e0ff7143e,
            0xb54e3a746addb473,
            0xd88e00e18c323864,
            0x00000000066e560a,
        ];
        const A_WORDS: [u64; 9] = [
            0x9501008408488925,
            0x0881002054411510,
            0x2525548924450402,
            0x2508548955211544,
            0x4910209521111104,
            0x8911080205550412,
            0x9542124422548410,
            0x4802002104120824,
            0x0000000002220202,
        ];
        let p = SECP256K1_P;
        let mut sx = Sampler::new(b"by-centered-clean560-x-v1", p);
        let mut sy = Sampler::new(b"by-centered-clean560-y-v1", p);
        let x = sx.next();
        let y = sy.next();
        let start_s = sw_centered_from_u256_for_test(addm(y, x, p), p);
        let mut delta = 1i64;
        let mut f = SInt::from_u(p);
        let mut g = SInt::from_u(x);
        let mut controls = Vec::with_capacity(STEPS);
        for _ in 0..STEPS {
            let odd_v = g.bit0();
            let a_v = delta > 0 && odd_v;
            controls.push((odd_v, a_v));
            divstep_sint_state(&mut delta, &mut f, &mut g);
        }
        for i in 0..STEPS {
            assert_eq!(controls[i].0, ((ODD_WORDS[i / 64] >> (i % 64)) & 1) != 0, "odd bit mismatch {i}");
            assert_eq!(controls[i].1, ((A_WORDS[i / 64] >> (i % 64)) & 1) != 0, "A bit mismatch {i}");
        }
        let mut b = super::super::B::new();
        let odd = b.alloc_qubits(STEPS);
        let a_ctrl = b.alloc_qubits(STEPS);
        let parity = b.alloc_qubits(STEPS);
        let r = b.alloc_qubits(WIDE);
        let s = b.alloc_qubits(WIDE);
        for i in 0..STEPS {
            if controls[i].0 { b.x(odd[i]); }
            if controls[i].1 { b.x(a_ctrl[i]); }
        }
        for i in 0..STEPS {
            emit_scaled_by_centered_signed_microstep_live_parity_variant_for_test(
                &mut b, &r, &s, odd[i], a_ctrl[i], parity[i], p, true, true,
            );
        }
        for i in (0..STEPS).rev() {
            emit_scaled_by_centered_signed_microstep_inverse_live_parity_variant_for_test(
                &mut b, &r, &s, odd[i], a_ctrl[i], parity[i], p, true, true,
            );
            emit_centered_signed_clear_parity_after_inverse_for_test(&mut b, &r, &s, odd[i], parity[i]);
        }
        for i in (0..STEPS).rev() {
            if controls[i].1 { b.x(a_ctrl[i]); }
            if controls[i].0 { b.x(odd[i]); }
        }
        let num_qubits = b.next_qubit as usize;
        let num_bits = b.next_bit as usize;
        let ops = b.ops;
        let mut hasher = sha3::Shake128::default();
        hasher.update(b"by-centered-clean560-fixed-hook-sim-v1");
        let mut xof = hasher.finalize_xof();
        let mut sim = crate::sim::Simulator::new(num_qubits, num_bits, &mut xof);
        set_slice_u512_by(&mut sim, &r, U512::ZERO);
        set_slice_u512_by(&mut sim, &s, sw_twos_for_width_for_test(start_s, WIDE));
        sim.apply(&ops);
        assert_eq!(get_slice_u512_by(&sim, &r), U512::ZERO, "fixed hook r not restored");
        assert_eq!(get_slice_u512_by(&sim, &s), sw_twos_for_width_for_test(start_s, WIDE), "fixed hook s not restored");
        for &q in odd.iter().chain(a_ctrl.iter()).chain(parity.iter()) {
            assert_eq!(sim.qubit(q) & 1, 0, "fixed hook control/parity not clean");
        }
        eprintln!("BY centered clean fixed trace hook test: phase={}, qubits={num_qubits}", sim.global_phase() & 1);
        assert_eq!(sim.global_phase() & 1, 0, "fixed hook all-exact roundtrip phase garbage");
    }

    #[test]
    fn live_reduction_flag_history_is_dense_and_high_entropy() {
        // If we cannot clean live reduction flags, can we at least treat them
        // as a small sparse history? No: they are arithmetic carry/borrow bits
        // of the tagged numerator channel and are common enough to resist a
        // simple position list. This does not kill compression, but it kills a
        // sparse-red-flag escape.
        let p = SECP256K1_P;
        let inv2 = (p.wrapping_add(U256::from(1u64))) >> 1usize;
        let samples = 3_000usize;
        let mut sx = Sampler::new(b"by-live-flag-density-x-v1", p);
        let mut sy = Sampler::new(b"by-live-flag-density-y-v1", p);
        let mut per_step_true = vec![0usize; 560];
        let mut counts = Vec::with_capacity(samples);
        let mut accepted = 0usize;
        while accepted < samples {
            let x = sx.next();
            let y = sy.next();
            let mut delta = 1i64;
            let mut f = SInt::from_u(p);
            let mut g = SInt::from_u(x);
            let mut controls = Vec::with_capacity(560);
            for _ in 0..560 {
                let odd = g.bit0();
                let a = delta > 0 && odd;
                controls.push((odd, a));
                divstep_sint_state(&mut delta, &mut f, &mut g);
            }
            if !g.is_zero() || !(f.is_one_pos() || f.is_one_neg()) { continue; }
            let mut r = U256::ZERO;
            let mut s = addm(y, x, p);
            let mut sample_true = 0usize;
            for (i, &(odd, a)) in controls.iter().enumerate() {
                let flag = if !odd {
                    false
                } else if a {
                    // After A swap+cneg, the add is (p-r_old)+s_old.  Fast
                    // cneg maps r=0 to p, which also sets the live flag.
                    r.is_zero() || s >= r
                } else {
                    s >= p.wrapping_sub(r)
                };
                if flag {
                    per_step_true[i] += 1;
                    sample_true += 1;
                }
                if a {
                    let nr = s;
                    let ns = mulm(subm(s, r, p), inv2, p);
                    r = nr;
                    s = ns;
                } else if odd {
                    s = mulm(addm(s, r, p), inv2, p);
                } else {
                    s = mulm(s, inv2, p);
                }
            }
            counts.push(sample_true);
            accepted += 1;
        }
        counts.sort_unstable();
        let mean_true = counts.iter().sum::<usize>() as f64 / samples as f64;
        let p90 = counts[(samples * 90) / 100];
        let p99 = counts[(samples * 99) / 100];
        let entropy: f64 = per_step_true
            .iter()
            .map(|&c| {
                let q = c as f64 / samples as f64;
                if q <= 0.0 || q >= 1.0 { 0.0 } else { -q * q.log2() - (1.0 - q) * (1.0 - q).log2() }
            })
            .sum();
        eprintln!(
            "BY live reduction flag history: mean_true={mean_true:.1}, p90={p90}, p99={p99}, independent_entropy≈{entropy:.1} bits"
        );
        assert!(mean_true > 80.0, "live flags unexpectedly sparse enough for a position-list escape");
        assert!(entropy > 250.0, "live flag history unexpectedly low entropy");
    }

    #[test]
    fn live_reduction_flag_is_recoverable_from_doubled_output_but_cleanup_is_costly() {
        // Algebra for cleaning a live modular-add reduction flag after the
        // following halve: for canonical inputs, z = 2*out_s mod p is the
        // modular-add result, so the reduction flag is (odd && z < addend).
        // Recovering it directly would require a modular double/copy plus a
        // full comparator, which costs more than the flag-uncompute we skipped.
        let p = SECP256K1_P;
        let mut b = super::super::B::new();
        let odd = b.alloc_qubit();
        let a_ctrl = b.alloc_qubit();
        let red_flag = b.alloc_qubit();
        let r = b.alloc_qubits(256);
        let s = b.alloc_qubits(256);
        emit_scaled_by_controlled_microstep_live_addflag_for_test(&mut b, &r, &s, odd, a_ctrl, red_flag, p);
        let num_qubits = b.next_qubit as usize;
        let num_bits = b.next_bit as usize;
        let ops = b.ops;
        let inv2 = (p.wrapping_add(U256::from(1u64))) >> 1usize;
        let cases = [(false, false, "C"), (true, false, "B"), (true, true, "A")];
        let mut sx = Sampler::new(b"by-live-flag-recover-r-v1", p);
        let mut sy = Sampler::new(b"by-live-flag-recover-s-v1", p);
        let mut mismatches = 0usize;
        let mut nonzero_mismatches = 0usize;
        let mut checked = 0usize;
        for &(odd_v, a_v, name) in &cases {
            let mut samples = vec![(U256::ZERO, U256::ZERO), (U256::ZERO, sy.next()), (sx.next(), U256::ZERO)];
            for _ in 0..16 { samples.push((sx.next(), sy.next())); }
            for (rv, sv) in samples {
                let (exp_r, exp_s) = match name {
                    "A" => (sv, mulm(subm(sv, rv, p), inv2, p)),
                    "B" => (rv, mulm(addm(sv, rv, p), inv2, p)),
                    "C" => (rv, mulm(sv, inv2, p)),
                    _ => unreachable!(),
                };
                let mut hasher = sha3::Shake128::default();
                hasher.update(b"by-live-flag-recover-sim-v1");
                let mut xof = hasher.finalize_xof();
                let mut sim = crate::sim::Simulator::new(num_qubits, num_bits, &mut xof);
                if odd_v { *sim.qubit_mut(odd) |= 1; }
                if a_v { *sim.qubit_mut(a_ctrl) |= 1; }
                set_slice_u512_by(&mut sim, &r, u256_to_u512_for_by_tests(rv));
                set_slice_u512_by(&mut sim, &s, u256_to_u512_for_by_tests(sv));
                sim.apply(&ops);
                let got_flag = (sim.qubit(red_flag) & 1) != 0;
                let z = addm(exp_s, exp_s, p);
                let recovered = odd_v && z < exp_r;
                if got_flag != recovered {
                    mismatches += 1;
                    if !(a_v && rv.is_zero()) {
                        nonzero_mismatches += 1;
                    }
                }
                checked += 1;
            }
        }
        let direct_cleanup_lb = 255 + 256 + 255; // copy/double output, compare, uncompute doubled copy.
        eprintln!(
            "BY live reduction flag recovery: checked={checked}, mismatches={mismatches}, nonzero_mismatches={nonzero_mismatches}, direct_cleanup_lb≈{direct_cleanup_lb} CCX/flag"
        );
        assert_eq!(nonzero_mismatches, 0, "canonical nonzero flag recovery relation failed");
        assert!(direct_cleanup_lb > 2 * 256, "direct cleanup lower bound no longer explains the blocker");
    }

    #[test]
    fn live_reduction_flag_microstep_hits_replay_target_but_needs_cleanup() {
        let p = SECP256K1_P;
        let mut b = super::super::B::new();
        let odd = b.alloc_qubit();
        let a_ctrl = b.alloc_qubit();
        let red_flag = b.alloc_qubit();
        let r = b.alloc_qubits(256);
        let s = b.alloc_qubits(256);
        emit_scaled_by_controlled_microstep_live_addflag_for_test(&mut b, &r, &s, odd, a_ctrl, red_flag, p);
        let ccx = count_ccx(&b.ops);
        let peak = b.peak_qubits;
        let num_qubits = b.next_qubit as usize;
        let num_bits = b.next_bit as usize;
        let ops = b.ops;
        let inv2 = (p.wrapping_add(U256::from(1u64))) >> 1usize;
        let cases = [(false, false, "C"), (true, false, "B"), (true, true, "A")];
        let mut sx = Sampler::new(b"by-live-flag-step-r-v1", p);
        let mut sy = Sampler::new(b"by-live-flag-step-s-v1", p);
        let mut saw_flag = false;
        for &(odd_v, a_v, name) in &cases {
            for _ in 0..12 {
                let rv = sx.next();
                let sv = sy.next();
                let (exp_r, exp_s) = match name {
                    "A" => (sv, mulm(subm(sv, rv, p), inv2, p)),
                    "B" => (rv, mulm(addm(sv, rv, p), inv2, p)),
                    "C" => (rv, mulm(sv, inv2, p)),
                    _ => unreachable!(),
                };
                let mut hasher = sha3::Shake128::default();
                hasher.update(b"by-live-flag-step-sim-v1");
                let mut xof = hasher.finalize_xof();
                let mut sim = crate::sim::Simulator::new(num_qubits, num_bits, &mut xof);
                if odd_v { *sim.qubit_mut(odd) |= 1; }
                if a_v { *sim.qubit_mut(a_ctrl) |= 1; }
                set_slice_u512_by(&mut sim, &r, u256_to_u512_for_by_tests(rv));
                set_slice_u512_by(&mut sim, &s, u256_to_u512_for_by_tests(sv));
                sim.apply(&ops);
                assert_eq!(get_slice_u512_by(&sim, &r), u256_to_u512_for_by_tests(exp_r), "r mismatch {name}");
                assert_eq!(get_slice_u512_by(&sim, &s), u256_to_u512_for_by_tests(exp_s), "s mismatch {name}");
                saw_flag |= (sim.qubit(red_flag) & 1) != 0;
            }
        }
        let replay560 = ccx * 560;
        eprintln!(
            "BY live-reduction-flag microstep: ccx={ccx}, replay560≈{replay560}, peak={peak}q"
        );
        assert!(saw_flag, "reduction flag never set in samples; test is not exercising live garbage");
        assert!(ccx < 1_850, "live-flag microstep does not recover replay margin");
        assert!(replay560 < 1_040_000, "live-flag replay not near target band");
    }

    #[test]
    fn scaled_by_controlled_microstep_matches_all_cases_and_hits_target_cost() {
        let p = SECP256K1_P;
        let mut b = super::super::B::new();
        let odd = b.alloc_qubit();
        let a_ctrl = b.alloc_qubit();
        let r = b.alloc_qubits(256);
        let s = b.alloc_qubits(256);
        emit_scaled_by_controlled_microstep_for_test(&mut b, &r, &s, odd, a_ctrl, p);
        let ccx = count_ccx(&b.ops);
        let peak = b.peak_qubits;
        let num_qubits = b.next_qubit as usize;
        let num_bits = b.next_bit as usize;
        let ops = b.ops;
        let inv2 = (p.wrapping_add(U256::from(1u64))) >> 1usize;
        let cases = [
            (false, false, "C"),
            (true, false, "B"),
            (true, true, "A"),
        ];
        let mut sx = Sampler::new(b"by-scaled-step-r-v1", p);
        let mut sy = Sampler::new(b"by-scaled-step-s-v1", p);
        for &(odd_v, a_v, name) in &cases {
            let mut samples = vec![(U256::ZERO, U256::ZERO), (U256::ZERO, sy.next()), (sx.next(), U256::ZERO)];
            for _ in 0..16 {
                samples.push((sx.next(), sy.next()));
            }
            for (rv, sv) in samples {
                let (exp_r, exp_s) = match name {
                    "A" => (sv, mulm(subm(sv, rv, p), inv2, p)),
                    "B" => (rv, mulm(addm(sv, rv, p), inv2, p)),
                    "C" => (rv, mulm(sv, inv2, p)),
                    _ => unreachable!(),
                };
                let mut hasher = sha3::Shake128::default();
                hasher.update(b"by-scaled-step-sim-v1");
                let mut xof = hasher.finalize_xof();
                let mut sim = crate::sim::Simulator::new(num_qubits, num_bits, &mut xof);
                if odd_v {
                    *sim.qubit_mut(odd) |= 1;
                }
                if a_v {
                    *sim.qubit_mut(a_ctrl) |= 1;
                }
                set_slice_u512_by(&mut sim, &r, u256_to_u512_for_by_tests(rv));
                set_slice_u512_by(&mut sim, &s, u256_to_u512_for_by_tests(sv));
                sim.apply(&ops);
                assert_eq!(get_slice_u512_by(&sim, &r), u256_to_u512_for_by_tests(exp_r), "r mismatch case {name}");
                assert_eq!(get_slice_u512_by(&sim, &s), u256_to_u512_for_by_tests(exp_s), "s mismatch case {name}");
            }
        }
        let total_560 = ccx as f64 * 560.0;
        eprintln!(
            "BY scaled controlled microstep: ccx={ccx}, total560≈{total_560:.0}, peak={peak}q"
        );
        assert!(total_560 < 1_250_000.0, "scaled controlled microsteps no longer SOTA-shaped");
        assert!(peak < 1_350, "scaled controlled microstep peak drifted too high");
    }

    fn emit_cmod_signed_mux_add_for_test(
        b: &mut super::super::B,
        acc: &[super::super::QubitId],
        a: &[super::super::QubitId],
        odd_ctrl: super::super::QubitId,
        neg_ctrl: super::super::QubitId,
        p: U256,
    ) {
        // Valid when neg_ctrl => odd_ctrl. Computes
        //   acc += odd_ctrl ? (neg_ctrl ? -a : a) : 0  (mod p).
        // It shares the ctrl&a addend for the add/sub cases, instead of paying
        // separate cmod_add and cmod_sub bodies.
        let n = acc.len();
        let f = b.alloc_qubits(n);
        for i in 0..n {
            b.ccx(odd_ctrl, a[i], f[i]);
        }
        for &q in &f {
            b.cx(neg_ctrl, q);
        }
        super::super::cadd_nbit_const_fast(b, &f, p.wrapping_add(U256::from(1u64)), neg_ctrl);
        super::super::mod_add_qq_fast(b, acc, &f, p);
        super::super::csub_nbit_const_fast(b, &f, p.wrapping_add(U256::from(1u64)), neg_ctrl);
        for &q in &f {
            b.cx(neg_ctrl, q);
        }
        for i in 0..n {
            let m = b.alloc_bit();
            b.hmr(f[i], m);
            b.cz_if(odd_ctrl, a[i], m);
        }
        b.free_vec(&f);
    }

    #[test]
    fn signed_mux_controlled_modular_add_works_but_not_enough() {
        let p = SECP256K1_P;
        let mut b = super::super::B::new();
        let odd = b.alloc_qubit();
        let neg = b.alloc_qubit();
        let acc = b.alloc_qubits(256);
        let a = b.alloc_qubits(256);
        emit_cmod_signed_mux_add_for_test(&mut b, &acc, &a, odd, neg, p);
        let ccx = count_ccx(&b.ops);
        let peak = b.peak_qubits;
        let num_qubits = b.next_qubit as usize;
        let num_bits = b.next_bit as usize;
        let ops = b.ops;
        let cases = [(false, false), (true, false), (true, true)];
        let mut sx = Sampler::new(b"by-signed-mux-acc-v1", p);
        let mut sy = Sampler::new(b"by-signed-mux-a-v1", p);
        for &(odd_v, neg_v) in &cases {
            for _ in 0..16 {
                let x = sx.next();
                let y = sy.next();
                let expected = if !odd_v {
                    x
                } else if neg_v {
                    subm(x, y, p)
                } else {
                    addm(x, y, p)
                };
                let mut hasher = sha3::Shake128::default();
                hasher.update(b"by-signed-mux-sim-v1");
                let mut xof = hasher.finalize_xof();
                let mut sim = crate::sim::Simulator::new(num_qubits, num_bits, &mut xof);
                if odd_v {
                    *sim.qubit_mut(odd) |= 1;
                }
                if neg_v {
                    *sim.qubit_mut(neg) |= 1;
                }
                set_slice_u512_by(&mut sim, &acc, u256_to_u512_for_by_tests(x));
                set_slice_u512_by(&mut sim, &a, u256_to_u512_for_by_tests(y));
                sim.apply(&ops);
                assert_eq!(get_slice_u512_by(&sim, &acc), u256_to_u512_for_by_tests(expected));
                assert_eq!(get_slice_u512_by(&sim, &a), u256_to_u512_for_by_tests(y));
            }
        }
        let static_a_total = 560.0 * (ccx as f64 + 1280.0 + 255.0) + 2.0 * 560.0 * 255.0;
        eprintln!(
            "BY signed mux controlled add/sub: ccx={ccx}, peak={peak}q, static_A_total≈{static_a_total:.0}"
        );
        assert!(ccx < 2_000, "signed mux failed to beat separate cmod add+sub");
        assert!(static_a_total > 2_000_000.0, "signed mux alone unexpectedly solves selected replay");
    }

    #[test]
    fn enumerated_branch_block_select_explodes_beyond_single_step() {
        // Another tempting idea is to group b divsteps and SELECT one fixed
        // branch-numerator block for each possible case sequence. Even ignoring
        // equality-control and QROM overhead, the sum of all fixed block bodies
        // grows too quickly once b>1.
        let mod_add = 1024usize;
        let mod_sub = 1277usize;
        let dbl = 255usize;
        let halve = 255usize;
        let case_cost = |c: char| match c {
            'A' => mod_sub + mod_add + dbl,
            'B' => mod_add + dbl,
            'C' => dbl,
            _ => unreachable!(),
        };
        let mut summaries = Vec::new();
        for block in 1usize..=4 {
            let mut seqs = std::collections::BTreeSet::<Vec<char>>::new();
            for delta0 in -80i64..=80 {
                for pat in 0usize..(1usize << block) {
                    let mut delta = delta0;
                    let mut seq = Vec::with_capacity(block);
                    for i in 0..block {
                        let odd = ((pat >> i) & 1) != 0;
                        if delta > 0 && odd {
                            seq.push('A');
                            delta = 1 - delta;
                        } else if odd {
                            seq.push('B');
                            delta = 1 + delta;
                        } else {
                            seq.push('C');
                            delta = 1 + delta;
                        }
                    }
                    seqs.insert(seq);
                }
            }
            let body_sum: usize = seqs.iter().map(|s| s.iter().map(|&c| case_cost(c)).sum::<usize>()).sum();
            let total = 560usize.div_ceil(block) * body_sum + 2 * 560 * halve;
            summaries.push((block, seqs.len(), body_sum, total));
        }
        eprintln!("BY enumerated branch-block SELECT lower bounds: {summaries:?}");
        assert!(summaries[1].3 > 5_000_000, "2-step enumerated SELECT unexpectedly viable");
        assert!(summaries[2].3 > 10_000_000, "3-step enumerated SELECT unexpectedly viable");
    }

    #[test]
    fn selected_replay_budget_requires_more_than_a_signed_mux() {
        // Quantify the remaining gap after the obvious primitive improvement.
        // A signed add/sub mux can combine the A-first (s-=r) and B-first
        // (s+=r) updates into one controlled modular operation per divstep.
        // But the extra A-only r+=s update is too common to pay at every step,
        // and not sparse enough for a naive position list. This budget tells us
        // what a real block-specialization scheme must beat.
        let p = SECP256K1_P;
        let mut b = super::super::B::new();
        let ctrl = b.alloc_qubit();
        let r = b.alloc_qubits(256);
        let s = b.alloc_qubits(256);
        let start = b.ops.len();
        super::super::cmod_add_qq(&mut b, &s, &r, ctrl, p);
        let cmod_add = count_ccx(&b.ops[start..]);
        let start = b.ops.len();
        super::super::mod_add_qq_fast(&mut b, &s, &r, p);
        let mod_add = count_ccx(&b.ops[start..]);
        let start = b.ops.len();
        super::super::mod_double_inplace_fast(&mut b, &r, p);
        let dbl = count_ccx(&b.ops[start..]);
        let start = b.ops.len();
        super::super::mod_halve_inplace_fast(&mut b, &r, p);
        let halve = count_ccx(&b.ops[start..]);

        let steps = 560.0;
        let scale_halves = 2.0 * steps * halve as f64;
        let ideal_signed_mux_static_a = steps * (2.0 * cmod_add as f64 + dbl as f64) + scale_halves;
        let mean_a = 133.5; // measured by actual_branch_cases_are_not_sparse_enough_for_a_correction_list.
        let signed_mux_with_value_proportional_a =
            steps * (cmod_add as f64 + dbl as f64) + mean_a * mod_add as f64 + scale_halves;
        eprintln!(
            "BY selected replay budget targets: cmod_add={cmod_add}, mod_add={mod_add}, dbl={dbl}, halve={halve}, signed_mux_static_A≈{ideal_signed_mux_static_a:.0}, signed_mux_value_A_lb≈{signed_mux_with_value_proportional_a:.0}"
        );
        assert!(ideal_signed_mux_static_a > 1_700_000.0, "static A mux would already be enough; revisit selected replay");
        assert!(signed_mux_with_value_proportional_a < 1_500_000.0, "even value-proportional A corrections are too costly");
    }

    #[test]
    fn actual_branch_cases_are_not_sparse_enough_for_a_correction_list() {
        // Check a tempting escape hatch: handle the odd add/sub stream with a
        // single signed mux per divstep, then encode the extra A-only r+=s
        // updates as a sparse correction list. Actual secp256k1 trajectories
        // kill this: A-cases are not rare after all, so a simple A-position
        // payload would be larger than raw branch history.
        const W: usize = 16;
        const WINDOWS: usize = 35;
        let samples = 10_000usize;
        let mut sampler = Sampler::new(b"by-actual-branch-case-dist-v1", SECP256K1_P);
        let mut a_counts = Vec::with_capacity(samples);
        let mut b_counts = Vec::with_capacity(samples);
        let mut c_counts = Vec::with_capacity(samples);
        for _ in 0..samples {
            let x = sampler.next();
            let mut delta = 1i64;
            let mut f = SInt::from_u(SECP256K1_P);
            let mut g = SInt::from_u(x);
            let mut ac = 0usize;
            let mut bc = 0usize;
            let mut cc = 0usize;
            for _ in 0..WINDOWS {
                for _ in 0..W {
                    let odd = g.bit0();
                    if delta > 0 && odd {
                        ac += 1;
                    } else if odd {
                        bc += 1;
                    } else {
                        cc += 1;
                    }
                    divstep_sint_state(&mut delta, &mut f, &mut g);
                }
            }
            a_counts.push(ac);
            b_counts.push(bc);
            c_counts.push(cc);
        }
        a_counts.sort_unstable();
        b_counts.sort_unstable();
        c_counts.sort_unstable();
        let mean_a = a_counts.iter().sum::<usize>() as f64 / samples as f64;
        let mean_b = b_counts.iter().sum::<usize>() as f64 / samples as f64;
        let mean_c = c_counts.iter().sum::<usize>() as f64 / samples as f64;
        let p99_a = a_counts[samples * 99 / 100];
        let p999_a = a_counts[samples * 999 / 1000];
        let sparse_a_bits_p99 = p99_a * 10; // 10 bits address one of 560 steps, loose fixed-list encoding.
        eprintln!(
            "BY actual branch cases over 560 steps: mean(A,B,C)=({mean_a:.1},{mean_b:.1},{mean_c:.1}), p99_A={p99_a}, p999_A={p999_a}, p99_A_position_bits≈{sparse_a_bits_p99}"
        );
        assert!(mean_a > 100.0, "A cases unexpectedly sparse; revisit correction-list idea");
        assert!(p99_a * 10 > 1_000, "A-position list unexpectedly compact");
    }

    #[test]
    fn fixed_branch_numerator_window_cost_distribution() {
        const W: usize = 16;
        let p_mod = SECP256K1_P;
        let mut hasher = sha3::Shake128::default();
        hasher.update(b"by-branch-numerator-cost-dist-v1");
        let mut reader = hasher.finalize_xof();
        let mut buf = [0u8; 24];
        let samples = 64usize;
        let mut costs = Vec::with_capacity(samples);
        let mut peaks = Vec::with_capacity(samples);
        let mut a_total = 0usize;
        let mut b_total = 0usize;
        let mut c_total = 0usize;
        for _ in 0..samples {
            reader.read(&mut buf);
            let f_low = (u64::from_le_bytes(buf[0..8].try_into().unwrap()) as i128) | 1;
            let g_low = u64::from_le_bytes(buf[8..16].try_into().unwrap()) as i128;
            let delta = (u64::from_le_bytes(buf[16..24].try_into().unwrap()) % 41) as i64 - 20;
            let bits = branch_bits_for_lowword_window(W, delta, f_low, g_low);
            let mut b = super::super::B::new();
            let x0 = b.alloc_qubits(256);
            let x1 = b.alloc_qubits(256);
            let (ac, bc, cc) = emit_fixed_branch_numerator_scaled_window_for_test(&mut b, delta, &bits, &x0, &x1, p_mod);
            a_total += ac;
            b_total += bc;
            c_total += cc;
            costs.push(count_ccx(&b.ops));
            peaks.push(b.peak_qubits as usize);
        }
        costs.sort_unstable();
        peaks.sort_unstable();
        let mean = costs.iter().sum::<usize>() as f64 / samples as f64;
        let p90 = costs[samples * 90 / 100];
        let max = costs[samples - 1];
        let approx35 = mean * 35.0;
        eprintln!(
            "BY fixed branch-numerator cost distribution: mean_ccx={mean:.1}, p90={p90}, max={max}, approx35≈{approx35:.0}, max_peak={}q, avg_cases=({:.2},{:.2},{:.2})",
            peaks[samples - 1],
            a_total as f64 / samples as f64,
            b_total as f64 / samples as f64,
            c_total as f64 / samples as f64
        );
        assert!(peaks[samples - 1] < 1_200, "branch-numerator distribution lost scratch advantage");
        assert!(approx35 < 900_000.0, "fixed branch-numerator arithmetic not SOTA-shaped");
    }

    #[test]
    fn fixed_hermite_inplace_window_cost_distribution() {
        const W: usize = 16;
        let p_mod = SECP256K1_P;
        let mut hasher = sha3::Shake128::default();
        hasher.update(b"by-hermite-inplace-cost-dist-v1");
        let mut reader = hasher.finalize_xof();
        let mut buf = [0u8; 24];
        let samples = 24usize;
        let mut costs = Vec::with_capacity(samples);
        let mut peaks = Vec::with_capacity(samples);
        let mut shears = Vec::with_capacity(samples);
        let mut ops_counts = Vec::with_capacity(samples);
        for _ in 0..samples {
            reader.read(&mut buf);
            let f_low = (u64::from_le_bytes(buf[0..8].try_into().unwrap()) as i128) | 1;
            let g_low = u64::from_le_bytes(buf[8..16].try_into().unwrap()) as i128;
            let delta = (u64::from_le_bytes(buf[16..24].try_into().unwrap()) % 41) as i64 - 20;
            let (_, _, _, pmat) = jump_matrix_direct_lowword(W, W, delta, f_low, g_low);
            let mut b = super::super::B::new();
            let x0 = b.alloc_qubits(256);
            let x1 = b.alloc_qubits(256);
            let (ops_count, max_shear, _) = emit_fixed_hermite_inplace_window_for_test(&mut b, pmat, &x0, &x1, p_mod);
            costs.push(count_ccx(&b.ops));
            peaks.push(b.peak_qubits as usize);
            shears.push(max_shear);
            ops_counts.push(ops_count);
        }
        costs.sort_unstable();
        peaks.sort_unstable();
        shears.sort_unstable();
        ops_counts.sort_unstable();
        let mean = costs.iter().sum::<usize>() as f64 / samples as f64;
        let p90 = costs[samples * 90 / 100];
        let max = costs[samples - 1];
        let peak_max = peaks[samples - 1];
        let approx35 = mean * 35.0;
        eprintln!(
            "BY Hermite in-place window cost distribution: mean_ccx={mean:.1}, p90={p90}, max={max}, approx35≈{approx35:.0}, max_peak={peak_max}q, max_shear={}, max_factor_ops={}",
            shears[samples - 1], ops_counts[samples - 1]
        );
        assert!(peak_max < 1_600, "Hermite in-place distribution lost scratch advantage");
        assert!(approx35 < 2_000_000.0, "naive fixed Hermite arithmetic is too costly to be SOTA-shaped");
    }

    fn inv_odd_mod_pow2_for_test(a: i128, bits: usize) -> i128 {
        if bits == 0 {
            return 0;
        }
        let modulus = 1i128 << bits;
        let (x, _, g) = egcd_i128_for_test(a.rem_euclid(modulus), modulus);
        assert_eq!(g, 1);
        x.rem_euclid(modulus)
    }

    fn h_ratio_step_for_test(delta: i64, h: i128, t: usize) -> (i64, i128, bool) {
        assert!(t >= 1);
        let odd = (h & 1) != 0;
        let next_bits = t - 1;
        if next_bits == 0 {
            let next_delta = if delta > 0 && odd { 1 - delta } else { 1 + delta };
            return (next_delta, 0, odd);
        }
        let next_mod = 1i128 << next_bits;
        if delta > 0 && odd {
            // h' = ((g-f)/2)/g = (h-1)/(2h) mod 2^(t-1).
            let inv_h = inv_odd_mod_pow2_for_test(h, next_bits);
            let next_h = (((h - 1) / 2) * inv_h).rem_euclid(next_mod);
            (1 - delta, next_h, odd)
        } else if odd {
            // h' = (g+f)/(2f) = (h+1)/2 mod 2^(t-1).
            (1 + delta, ((h + 1) / 2).rem_euclid(next_mod), odd)
        } else {
            // h' = g/(2f) = h/2 mod 2^(t-1).
            (1 + delta, (h / 2).rem_euclid(next_mod), odd)
        }
    }

    #[test]
    fn low_ratio_microstep_update_generates_branch_bits_without_full_denominator() {
        // BY branch generation does not need a full 256-bit denominator pair.
        // With h=g/f mod 2^t and odd f, the next branch bit is h&1 and h has
        // a closed 2-adic update. This keeps the selector generator at
        // O(w)-bit state; the hard part is selected modular replay, not finding
        // the branch bits.
        const W: usize = 16;
        let mut hasher = sha3::Shake128::default();
        hasher.update(b"by-low-ratio-step-v1");
        let mut reader = hasher.finalize_xof();
        let mut buf = [0u8; 24];
        for _ in 0..20_000 {
            reader.read(&mut buf);
            let mut f = truncate_i128((u64::from_le_bytes(buf[0..8].try_into().unwrap()) as i128) | 1, W);
            let mut g = truncate_i128(u64::from_le_bytes(buf[8..16].try_into().unwrap()) as i128, W);
            let mut delta = (u64::from_le_bytes(buf[16..24].try_into().unwrap()) % 41) as i64 - 20;
            for t in (1..=W).rev() {
                f = truncate_i128(f, t);
                g = truncate_i128(g, t);
                let modulus = 1i128 << t;
                let h = (g.rem_euclid(modulus) * inv_odd_mod_pow2_for_test(f, t)).rem_euclid(modulus);
                let (next_delta_h, next_h, odd_h) = h_ratio_step_for_test(delta, h, t);
                let odd_g = (g & 1) != 0;
                assert_eq!(odd_h, odd_g, "h parity did not match g parity");
                if delta > 0 && odd_g {
                    let nf = g;
                    let ng = (g - f) / 2;
                    delta = 1 - delta;
                    f = nf;
                    g = ng;
                } else if odd_g {
                    g = (g + f) / 2;
                    delta = 1 + delta;
                } else {
                    g /= 2;
                    delta = 1 + delta;
                }
                assert_eq!(delta, next_delta_h, "delta update mismatch");
                if t > 1 {
                    g = truncate_i128(g, t - 1);
                    f = truncate_i128(f, t); // next loop truncates f to t-1 at top.
                    let next_mod = 1i128 << (t - 1);
                    let f_next = truncate_i128(f, t - 1);
                    let g_next = truncate_i128(g, t - 1);
                    let h_from_fg = (g_next.rem_euclid(next_mod) * inv_odd_mod_pow2_for_test(f_next, t - 1)).rem_euclid(next_mod);
                    assert_eq!(next_h, h_from_fg, "h update mismatch at t={t}");
                }
            }
        }
        eprintln!("BY low-ratio branch generator: 16-bit h+delta state suffices per window; branch history is the reversibility payload");
    }

    fn branch_bits_for_lowword_window(w: usize, mut delta: i64, mut f: i128, mut g: i128) -> Vec<bool> {
        let mut bits = Vec::with_capacity(w);
        f = truncate_i128(f, w);
        g = truncate_i128(g, w);
        for t in (1..=w).rev() {
            f = truncate_i128(f, t);
            let odd = (g & 1) != 0;
            bits.push(odd);
            if delta > 0 && odd {
                let nf = g;
                let ng = (g - f) / 2;
                delta = 1 - delta;
                f = nf;
                g = ng;
            } else if odd {
                g = (g + f) / 2;
                delta = 1 + delta;
            } else {
                g /= 2;
                delta = 1 + delta;
            }
            g = truncate_i128(g, t - 1);
        }
        bits
    }

    fn matrix_from_branch_bits(mut delta: i64, bits: &[bool]) -> TransitionMatrix {
        let (mut p00, mut p01, mut p10, mut p11) = (1i128, 0i128, 0i128, 1i128);
        for &odd in bits {
            if delta > 0 && odd {
                let (np00, np01, np10, np11) = (2 * p10, 2 * p11, -p00 + p10, -p01 + p11);
                delta = 1 - delta;
                p00 = np00;
                p01 = np01;
                p10 = np10;
                p11 = np11;
            } else if odd {
                let (np00, np01, np10, np11) = (2 * p00, 2 * p01, p00 + p10, p01 + p11);
                delta = 1 + delta;
                p00 = np00;
                p01 = np01;
                p10 = np10;
                p11 = np11;
            } else {
                p00 *= 2;
                p01 *= 2;
                delta = 1 + delta;
            }
        }
        TransitionMatrix { m00: p00, m01: p01, m10: p10, m11: p11, delta_final: delta }
    }

    #[test]
    fn branch_bits_reconstruct_by_jump_matrix() {
        const W: usize = 16;
        let mut hasher = sha3::Shake128::default();
        hasher.update(b"by-branch-bits-reconstruct-v1");
        let mut reader = hasher.finalize_xof();
        let mut buf = [0u8; 24];
        for _ in 0..10_000 {
            reader.read(&mut buf);
            let f_low = (u64::from_le_bytes(buf[0..8].try_into().unwrap()) as i128) | 1;
            let g_low = u64::from_le_bytes(buf[8..16].try_into().unwrap()) as i128;
            let delta = (u64::from_le_bytes(buf[16..24].try_into().unwrap()) % 41) as i64 - 20;
            let (_, _, _, direct) = jump_matrix_direct_lowword(W, W, delta, f_low, g_low);
            let bits = branch_bits_for_lowword_window(W, delta, f_low, g_low);
            let rebuilt = matrix_from_branch_bits(delta, &bits);
            assert_eq!(rebuilt, direct);
        }
    }

    #[test]
    fn branch_bit_history_by_tagged_div_budget_model() {
        // Each w=16 BY matrix is exactly determined by the 16 odd/even branch
        // bits plus the starting delta. This gives a concrete 560-bit selector
        // history for 35 windows, unlike empirical entropy codes. It does not
        // solve generation of those bits from x; it only makes the matrix
        // selector representation compatible with the 2800q cap.
        const WINDOWS: usize = 35;
        const BRANCH_HISTORY_BITS: usize = WINDOWS * 16;
        const DELTA_AND_CONTROL: usize = 16;
        const MOD_PEAK: usize = 2224;
        let peak_model = MOD_PEAK + BRANCH_HISTORY_BITS + DELTA_AND_CONTROL;
        eprintln!(
            "BY branch-bit history budget: branch_bits={BRANCH_HISTORY_BITS}, peak_model≈{peak_model}q"
        );
        assert!(peak_model <= 2_800, "branch-bit selector history does not fit cap");
    }

    #[test]
    fn h_only_compressed_history_by_tagged_div_budget_model() {
        // Structural target model: delete the full integer denominator pair and
        // keep only the 16-bit low ratio h plus delta. Matrix/history is stored
        // in a compressed code (empirical p99 below), and arithmetic is just the
        // modular fixed-matrix replacement per window. This is not a circuit,
        // but it is the first BY model that is simultaneously sub-MToffoli and
        // under the 2800q cap.
        const WINDOWS: usize = 35;
        const WIDTH: usize = 274;
        let mut hasher = sha3::Shake128::default();
        hasher.update(b"by-h-only-compressed-budget-v1");
        let mut reader = hasher.finalize_xof();
        let mut buf = [0u8; 24];
        let samples = 24usize;
        let mut mod_costs = Vec::with_capacity(samples);
        let mut mod_peaks = Vec::with_capacity(samples);
        for _ in 0..samples {
            reader.read(&mut buf);
            let f_low = (u64::from_le_bytes(buf[0..8].try_into().unwrap()) as i128) | 1;
            let g_low = u64::from_le_bytes(buf[8..16].try_into().unwrap()) as i128;
            let delta = (u64::from_le_bytes(buf[16..24].try_into().unwrap()) % 41) as i64 - 20;
            let (_, _, _, mtx) = jump_matrix_direct_lowword(16, 16, delta, f_low, g_low);
            let mut b_mod = super::super::B::new();
            let x0 = b_mod.alloc_qubits(256);
            let x1 = b_mod.alloc_qubits(256);
            let y0 = b_mod.alloc_qubits(WIDTH);
            let y1 = b_mod.alloc_qubits(WIDTH);
            emit_signed_row_scaled_from_sources_for_test(&mut b_mod, mtx.m00, &x0, mtx.m01, &x1, &y0);
            emit_signed_row_scaled_from_sources_for_test(&mut b_mod, mtx.m10, &x0, mtx.m11, &x1, &y1);
            let _regs = emit_fixed_matrix_old_cleanup_for_test(&mut b_mod, mtx, &x0, &x1, &y0, &y1);
            mod_costs.push(count_ccx(&b_mod.ops));
            mod_peaks.push(b_mod.peak_qubits as usize);
        }
        mod_costs.sort_unstable();
        mod_peaks.sort_unstable();
        let mean_mod_window = mod_costs.iter().sum::<usize>() as f64 / samples as f64;
        let approx_arith = mean_mod_window * WINDOWS as f64;

        // Conservative p99 code length from the 10k entropy experiment, rounded
        // up with margin; h/delta/control allowance covers live low-ratio state.
        let history_bits = 480usize;
        let h_delta_control = 32usize;
        let peak_model = mod_peaks[samples - 1] + history_bits + h_delta_control;
        eprintln!(
            "BY h-only compressed-history budget: mean_mod_window≈{mean_mod_window:.1}, approx35≈{approx_arith:.0}, mod_peak={}q, history_bits={history_bits}, peak_model≈{peak_model}q",
            mod_peaks[samples - 1]
        );
        assert!(approx_arith < 900_000.0, "h-only BY arithmetic no longer sub-MToffoli");
        assert!(peak_model < 2_800, "h-only compressed BY model exceeds cap");
    }

    #[test]
    fn by_tagged_div_stored_matrix_upper_bound_model() {
        // Upper-bound architecture with per-window matrix history already known:
        // update the integer denominator pair with sparse scaled rows, and the
        // modular tagged pair with the fixed-matrix replacement developed above.
        // This separates arithmetic viability from the remaining matrix-selection
        // / history-compression problem.
        const WIDTH: usize = 274;
        let mut hasher = sha3::Shake128::default();
        hasher.update(b"by-tagged-div-stored-matrix-upper-bound-v1");
        let mut reader = hasher.finalize_xof();
        let mut buf = [0u8; 24];
        let samples = 32usize;
        let mut window_costs = Vec::with_capacity(samples);
        let mut mod_peaks = Vec::with_capacity(samples);
        let mut den_peaks = Vec::with_capacity(samples);
        for _ in 0..samples {
            reader.read(&mut buf);
            let f_low = (u64::from_le_bytes(buf[0..8].try_into().unwrap()) as i128) | 1;
            let g_low = u64::from_le_bytes(buf[8..16].try_into().unwrap()) as i128;
            let delta = (u64::from_le_bytes(buf[16..24].try_into().unwrap()) % 41) as i64 - 20;
            let (_, _, _, mtx) = jump_matrix_direct_lowword(16, 16, delta, f_low, g_low);

            let mut b_den = super::super::B::new();
            emit_scaled_pair_update_with_cleanup_for_cost(&mut b_den, mtx, WIDTH, 16);
            let den_ccx = count_ccx(&b_den.ops);
            den_peaks.push(b_den.peak_qubits as usize);

            let mut b_mod = super::super::B::new();
            let x0 = b_mod.alloc_qubits(256);
            let x1 = b_mod.alloc_qubits(256);
            let y0 = b_mod.alloc_qubits(WIDTH);
            let y1 = b_mod.alloc_qubits(WIDTH);
            emit_signed_row_scaled_from_sources_for_test(&mut b_mod, mtx.m00, &x0, mtx.m01, &x1, &y0);
            emit_signed_row_scaled_from_sources_for_test(&mut b_mod, mtx.m10, &x0, mtx.m11, &x1, &y1);
            let _regs = emit_fixed_matrix_old_cleanup_for_test(&mut b_mod, mtx, &x0, &x1, &y0, &y1);
            let mod_ccx = count_ccx(&b_mod.ops);
            mod_peaks.push(b_mod.peak_qubits as usize);
            window_costs.push(den_ccx + mod_ccx);
        }
        window_costs.sort_unstable();
        mod_peaks.sort_unstable();
        den_peaks.sort_unstable();
        let mean_window = window_costs.iter().sum::<usize>() as f64 / samples as f64;
        let p90_window = window_costs[(samples * 90) / 100];
        let max_window = window_costs[samples - 1];
        let approx_total = mean_window * 35.0;
        let stored_key_bits = 35 * 22; // delta plus h=g/f mod 2^16 upper-bound selector.
        let scheduled_peak_model = mod_peaks[samples - 1] + 2 * WIDTH;
        eprintln!(
            "BY tagged-DIV stored-matrix upper bound: mean_window_ccx={mean_window:.1}, p90={p90_window}, max={max_window}, approx35≈{approx_total:.0}, den_peak={}q, mod_peak={}q, scheduled_peak≈{scheduled_peak_model}q, selector_bits={stored_key_bits}",
            den_peaks[samples - 1],
            mod_peaks[samples - 1]
        );
        assert!(approx_total < 1_200_000.0, "stored-matrix BY arithmetic no longer cheaper than Kaliski");
        assert!(scheduled_peak_model < 2_900, "stored-matrix BY upper-bound peak drifted too high");
    }

    #[test]
    fn qcorr_roundtrip_recovers_m_for_sampled_by_matrices() {
        // If q = s*adj(P)*m / 2^w, then P*q = m. This is the missing
        // reversibility hook for general old-row cleanup: after q has been used
        // to remove q*c from the old rows, m can be uncomputed from q even
        // though the old sources have been zeroed.
        const W: usize = 16;
        let pinv = 51_919i128;
        let mask = (1i128 << W) - 1;
        let mut hasher = sha3::Shake128::default();
        hasher.update(b"by-qcorr-roundtrip-v1");
        let mut reader = hasher.finalize_xof();
        let mut buf = [0u8; 24];
        for _ in 0..5_000 {
            reader.read(&mut buf);
            let f_low = (u64::from_le_bytes(buf[0..8].try_into().unwrap()) as i128) | 1;
            let g_low = u64::from_le_bytes(buf[8..16].try_into().unwrap()) as i128;
            let delta = (u64::from_le_bytes(buf[16..24].try_into().unwrap()) % 41) as i64 - 20;
            let (_, _, _, mtx) = jump_matrix_direct_lowword(W, W, delta, f_low, g_low);
            // Use deterministic low row values; only low words matter for m.
            let x0_low = (f_low * 17 + 3) & mask;
            let x1_low = (g_low * 19 - 5) & mask;
            let t0_low = (mtx.m00 * x0_low + mtx.m01 * x1_low) & mask;
            let t1_low = (mtx.m10 * x0_low + mtx.m11 * x1_low) & mask;
            let m0 = (-t0_low * pinv) & mask;
            let m1 = (-t1_low * pinv) & mask;
            let sgn = det_sign_pow2(mtx, W);
            let q0_num = sgn * (mtx.m11 * m0 - mtx.m01 * m1);
            let q1_num = sgn * (-mtx.m10 * m0 + mtx.m00 * m1);
            assert_eq!(q0_num & mask, 0);
            assert_eq!(q1_num & mask, 0);
            let q0 = q0_num >> W;
            let q1 = q1_num >> W;
            assert_eq!(mtx.m00 * q0 + mtx.m01 * q1, m0, "P*q did not recover m0");
            assert_eq!(mtx.m10 * q0 + mtx.m11 * q1, m1, "P*q did not recover m1");
        }
    }

    #[test]
    fn adjugate_m_correction_is_integral_for_sampled_by_matrices() {
        // General cleanup formula behind the triangular prototype. If
        // 2^w*y = P*x + p*m and det(P)=s*2^w, then
        //   s*adj(P)*y = x + p * (s*adj(P)*m / 2^w).
        // The correction vector is integral for valid BY rows, so old-row
        // cleanup can in principle use the same low-word m values computed
        // from the original sources, not a dense modular inverse matrix.
        const W: usize = 16;
        let p = SECP256K1_P;
        let pinv = 51_919i128;
        let mask = (1i128 << W) - 1;
        let mut hasher = sha3::Shake128::default();
        hasher.update(b"by-adjugate-m-integral-v1");
        let mut reader = hasher.finalize_xof();
        let mut buf = [0u8; 88];
        for _ in 0..2_000 {
            reader.read(&mut buf);
            let f_low = (u64::from_le_bytes(buf[0..8].try_into().unwrap()) as i128) | 1;
            let g_low = u64::from_le_bytes(buf[8..16].try_into().unwrap()) as i128;
            let delta = (u64::from_le_bytes(buf[16..24].try_into().unwrap()) % 41) as i64 - 20;
            let x0 = U256::from_le_slice(&buf[24..56]) % p;
            let x1 = U256::from_le_slice(&buf[56..88]) % p;
            let (_, _, _, mtx) = jump_matrix_direct_lowword(W, W, delta, f_low, g_low);
            let t0_low = (mtx.m00 * (x0.as_limbs()[0] as i128) + mtx.m01 * (x1.as_limbs()[0] as i128)) & mask;
            let t1_low = (mtx.m10 * (x0.as_limbs()[0] as i128) + mtx.m11 * (x1.as_limbs()[0] as i128)) & mask;
            let m0 = (-t0_low * pinv) & mask;
            let m1 = (-t1_low * pinv) & mask;
            let sgn = det_sign_pow2(mtx, W);
            let c0_num = sgn * (mtx.m11 * m0 - mtx.m01 * m1);
            let c1_num = sgn * (-mtx.m10 * m0 + mtx.m00 * m1);
            assert_eq!(c0_num & mask, 0, "adjugate m correction 0 not divisible by 2^w");
            assert_eq!(c1_num & mask, 0, "adjugate m correction 1 not divisible by 2^w");
        }
    }

    #[test]
    fn positive_triangular_fixed_matrix_replacement_cleans_old_rows() {
        const WIDTH: usize = 274;
        let mtx = jump_matrix_direct_lowword(16, 16, -20, 1, 1).3;
        assert_eq!((mtx.m00, mtx.m01, mtx.m10, mtx.m11), (65536, 0, 65535, 1));
        let mut b = super::super::B::new();
        let x0 = b.alloc_qubits(256);
        let x1 = b.alloc_qubits(256);
        let y0 = b.alloc_qubits(WIDTH);
        let y1 = b.alloc_qubits(WIDTH);
        emit_positive_row_scaled_from_sources_for_test(&mut b, mtx.m00, &x0, mtx.m01, &x1, &y0);
        emit_positive_row_scaled_from_sources_for_test(&mut b, mtx.m10, &x0, mtx.m11, &x1, &y1);
        let (m_reg, z_reg) = emit_positive_triangular_old_cleanup_for_test(&mut b, &x0, &x1, &y0, &y1);
        let ccx = count_ccx(&b.ops);
        let peak = b.peak_qubits;
        let num_qubits = b.next_qubit as usize;
        let num_bits = b.next_bit as usize;
        let ops = b.ops;
        let p512 = u256_to_u512_for_by_tests(SECP256K1_P);
        let pinv = 51_919u64;
        let mask = (1u64 << 16) - 1;
        let mut sx = Sampler::new(b"by-positive-tri-repl-x0-v1", SECP256K1_P);
        let mut sy = Sampler::new(b"by-positive-tri-repl-x1-v1", SECP256K1_P);
        for _ in 0..32 {
            let a = sx.next();
            let c = sy.next();
            let exp0 = u256_to_u512_for_by_tests(a);
            let t1 = u256_to_u512_for_by_tests(a) * U512::from(65535u64)
                + u256_to_u512_for_by_tests(c);
            let corr1 = (t1.as_limbs()[0] & mask).wrapping_mul((!pinv).wrapping_add(1)) & mask;
            let exp1: U512 = (t1 + U512::from(corr1) * p512) >> 16usize;
            let mut hasher = sha3::Shake128::default();
            hasher.update(b"by-positive-tri-repl-sim-v1");
            let mut xof = hasher.finalize_xof();
            let mut sim = crate::sim::Simulator::new(num_qubits, num_bits, &mut xof);
            set_slice_u512_by(&mut sim, &x0, u256_to_u512_for_by_tests(a));
            set_slice_u512_by(&mut sim, &x1, u256_to_u512_for_by_tests(c));
            sim.apply(&ops);
            assert_eq!(get_slice_u512_by(&sim, &x0), U512::ZERO, "x0 not zero");
            assert_eq!(get_slice_u512_by(&sim, &x1), U512::ZERO, "x1 not zero");
            assert_eq!(get_slice_u512_by(&sim, &m_reg), U512::ZERO, "m not zero");
            assert_eq!(get_slice_u512_by(&sim, &z_reg), U512::ZERO, "z not zero");
            assert_eq!(get_slice_u512_by(&sim, &y0), exp0, "y0 changed");
            assert_eq!(get_slice_u512_by(&sim, &y1), exp1, "y1 mismatch");
        }
        eprintln!(
            "positive triangular BY fixed-matrix replacement: ccx={ccx}, peak={peak}q"
        );
        assert!(ccx < 35_000, "fixed positive replacement too costly");
        assert!(peak < 2_500, "fixed positive replacement peak too high");
    }

    #[test]
    fn noncanonical_scaled_pair_map_is_injective_on_canonical_domain() {
        // Row scaling alone loses representative quotient (T and T+p collide),
        // but the TWO-row matrix map can still be injective on canonical input
        // pairs because det(P)=±2^w and p is odd. This is the algebraic reason
        // a fixed-matrix pair replacement might clean quotient bits using both
        // rows/sources instead of storing m histories.
        use std::collections::HashSet;
        let p_small: i128 = 251;
        let w = 4usize;
        let mask = (1i128 << w) - 1;
        let pinv = 3i128; // 251^{-1} mod 16.
        let matrices = [
            jump_matrix_direct_lowword(w, w, 1, 1, 3).3,
            jump_matrix_direct_lowword(w, w, -3, 1, 5).3,
            jump_matrix_direct_lowword(w, w, 7, 1, -2).3,
            jump_matrix_direct_lowword(w, w, 0, 1, 6).3,
        ];
        for mtx in matrices {
            det_sign_pow2(mtx, w);
            let mut seen = HashSet::new();
            for x0 in 0..p_small {
                for x1 in 0..p_small {
                    let t0 = mtx.m00 * x0 + mtx.m01 * x1;
                    let t1 = mtx.m10 * x0 + mtx.m11 * x1;
                    let c0 = (-(t0 & mask) * pinv) & mask;
                    let c1 = (-(t1 & mask) * pinv) & mask;
                    let q0 = (t0 + c0 * p_small) >> w;
                    let q1 = (t1 + c1 * p_small) >> w;
                    assert!(seen.insert((q0, q1)), "collision for matrix {:?} at ({x0},{x1})", mtx);
                }
            }
        }
    }

    #[test]
    fn noncanonical_batched_shift_needs_quotient_uncompute() {
        // Important caveat for the highfold idea: for noncanonical T, the final
        // scaled residue does not uniquely encode the quotient k such that
        // T=k*p+R. T and T+p represent the same residue and produce the same
        // scaled output, but their low-word correction m differs by one. A
        // reversible circuit must therefore either keep k, recover it from the
        // row sources, or fuse reduction with cleanup; it cannot just erase k
        // from the output row alone.
        let p = SECP256K1_P;
        let p512 = u256_to_u512_for_by_tests(p);
        let pinv = 51_919u64;
        let mask = (1u64 << 16) - 1;
        let t = U256::from(123456789u64);
        let low0 = t.as_limbs()[0] & mask;
        let m0 = low0.wrapping_mul((!pinv).wrapping_add(1)) & mask;
        let q0: U512 = (u256_to_u512_for_by_tests(t) + U512::from(m0) * p512) >> 16usize;
        let t1 = u256_to_u512_for_by_tests(t) + p512;
        let low1 = t1.as_limbs()[0] & mask;
        let m1 = low1.wrapping_mul((!pinv).wrapping_add(1)) & mask;
        let q1: U512 = (t1 + U512::from(m1) * p512) >> 16usize;
        assert_eq!(q0, q1, "scaled residue should ignore representative quotient");
        assert_ne!(m0, m1, "correction m should change with representative quotient");
    }

    #[test]
    fn highfold_then_batched_halve16_matches_row_distribution() {
        // For actual BY row values T=a*x+b*y with signed w=16 matrix entries,
        // first folding k=T>>256 copies of p brings T into canonical range, and
        // then the batched halve's top-bit m recovery succeeds on samples.
        let p_u = u256_to_u512_for_by_tests(SECP256K1_P);
        let pinv = 51_919u64;
        let mask = (1u64 << 16) - 1;
        let mut hasher = sha3::Shake128::default();
        hasher.update(b"by-row-highfold-batched-halve-v1");
        let mut reader = hasher.finalize_xof();
        let mut buf = [0u8; 88];
        let samples = 20_000usize;
        let mut failures = 0usize;
        for _ in 0..samples {
            reader.read(&mut buf);
            let f_low = (u64::from_le_bytes(buf[0..8].try_into().unwrap()) as i128) | 1;
            let g_low = u64::from_le_bytes(buf[8..16].try_into().unwrap()) as i128;
            let delta = (u64::from_le_bytes(buf[16..24].try_into().unwrap()) % 41) as i64 - 20;
            let x = U256::from_le_slice(&buf[24..56]) % SECP256K1_P;
            let y = U256::from_le_slice(&buf[56..88]) % SECP256K1_P;
            let (_, _, _, mtx) = jump_matrix_direct_lowword(16, 16, delta, f_low, g_low);
            for &(a, bb) in &[(mtx.m00, mtx.m01), (mtx.m10, mtx.m11)] {
                // Use i128 for the small high quotient and U512 for positive
                // magnitude arithmetic; sampled signs are handled by checking
                // both row signs through signed_i128_mod_p equivalence.
                let ax = if a >= 0 { u256_to_u512_for_by_tests(x) * U512::from(a as u128) } else { U512::ZERO };
                let by = if bb >= 0 { u256_to_u512_for_by_tests(y) * U512::from(bb as u128) } else { U512::ZERO };
                if a < 0 || bb < 0 {
                    // Fall back to modular representative for signed rows in
                    // this distribution test; the circuit cost model below is
                    // sign-symmetric.
                    let row_mod = addm(mulm(signed_i128_mod_p(a, SECP256K1_P), x, SECP256K1_P), mulm(signed_i128_mod_p(bb, SECP256K1_P), y, SECP256K1_P), SECP256K1_P);
                    let low = row_mod.as_limbs()[0] & mask;
                    let corr = low.wrapping_mul((!pinv).wrapping_add(1)) & mask;
                    let q: U512 = (u256_to_u512_for_by_tests(row_mod) + U512::from(corr) * p_u) >> 16usize;
                    let q_top: U512 = q >> 240usize;
                    let top = q_top.to::<u64>() & mask;
                    if top != corr { failures += 1; }
                } else {
                    let t = ax + by;
                    let k: U512 = t >> 256usize;
                    let folded = t - k * p_u;
                    let low = folded.as_limbs()[0] & mask;
                    let corr = low.wrapping_mul((!pinv).wrapping_add(1)) & mask;
                    let q: U512 = (folded + U512::from(corr) * p_u) >> 16usize;
                    let q_top: U512 = q >> 240usize;
                    let top = q_top.to::<u64>() & mask;
                    if top != corr { failures += 1; }
                }
            }
        }
        eprintln!("BY row highfold+halve16 sampled failures={failures}/{}", samples * 2);
        assert_eq!(failures, 0);
    }

    #[test]
    fn approximate_batched_shift_reopens_scaled_by_jump_budget() {
        const WIDTH: usize = 274;
        const W: usize = 16;
        let mut b = super::super::B::new();
        let v = b.alloc_qubits(WIDTH);
        let start = b.ops.len();
        emit_approx_highfold_p_for_cost(&mut b, &v);
        let highfold_ccx = count_ccx(&b.ops[start..]);
        let start_shift = b.ops.len();
        emit_approx_batched_halve16_for_cost(&mut b, &v);
        let shift_ccx = count_ccx(&b.ops[start_shift..]);

        let mut hasher = sha3::Shake128::default();
        hasher.update(b"by-approx-batched-shift-budget-v1");
        let mut reader = hasher.finalize_xof();
        let mut buf = [0u8; 24];
        let samples = 24usize;
        let mut total_pair_ccx = 0usize;
        for _ in 0..samples {
            reader.read(&mut buf);
            let f_low = (u64::from_le_bytes(buf[0..8].try_into().unwrap()) as i128) | 1;
            let g_low = u64::from_le_bytes(buf[8..16].try_into().unwrap()) as i128;
            let delta = (u64::from_le_bytes(buf[16..24].try_into().unwrap()) % 41) as i64 - 20;
            let (_, _, _, m) = jump_matrix_direct_lowword(W, W, delta, f_low, g_low);
            let mut b2 = super::super::B::new();
            emit_scaled_pair_update_with_cleanup_for_cost(&mut b2, m, WIDTH, W);
            total_pair_ccx += count_ccx(&b2.ops);
        }
        let mean_integer_pair = total_pair_ccx as f64 / samples as f64;
        let row_scale_ccx = highfold_ccx + shift_ccx;
        // Two forward rows need highfold+shift. Two old rows cleaned by the
        // sparse adjugate need a highfold to turn the residual small multiple
        // of p into zero. The base integer_pair already includes the sparse
        // row additions/subtractions themselves.
        let modular_pair_window = mean_integer_pair + 2.0 * row_scale_ccx as f64 + 2.0 * highfold_ccx as f64;
        let approx35 = modular_pair_window * 35.0;
        eprintln!(
            "approx batched-shift BY scaled modular budget: highfold_ccx={highfold_ccx}, shift16_ccx={shift_ccx}, integer_pair≈{mean_integer_pair:.1}, modular_pair/window≈{modular_pair_window:.1}, approx35≈{approx35:.0}, shift_peak={}q",
            b.peak_qubits
        );
        assert!(approx35 < 800_000.0, "batched shift no longer gives a SOTA-shaped BY modular pair");
    }

    #[test]
    fn scaled_modular_jump_sparse_cleanup_is_too_expensive_with_current_primitives() {
        // Tried repair after discovering dense unscaled inverses: keep the
        // coefficient/tagged channel in the scaled BY convention. A window then
        // costs sparse forward P rows, public halvings by w, and sparse
        // adjugate cleanup. With the current constant-multiply/halve primitives
        // this is still too expensive; keep the result as an invalidation and
        // as a target for a better small-constant modular row former.
        const W: usize = 16;
        let p = SECP256K1_P;
        let mut hasher = sha3::Shake128::default();
        hasher.update(b"by-scaled-modular-sparse-cleanup-v1");
        let mut reader = hasher.finalize_xof();
        let mut buf = [0u8; 24];
        let samples = 12usize;
        let mut total_ccx = 0usize;
        let mut max_peak = 0u32;
        for _ in 0..samples {
            reader.read(&mut buf);
            let f_low = (u64::from_le_bytes(buf[0..8].try_into().unwrap()) as i128) | 1;
            let g_low = u64::from_le_bytes(buf[8..16].try_into().unwrap()) as i128;
            let delta = (u64::from_le_bytes(buf[16..24].try_into().unwrap()) % 41) as i64 - 20;
            let (_, _, _, m) = jump_matrix_direct_lowword(W, W, delta, f_low, g_low);
            let mut b = super::super::B::new();
            emit_scaled_modular_pair_update_with_sparse_cleanup_for_cost(&mut b, m, W, p);
            total_ccx += count_ccx(&b.ops);
            max_peak = max_peak.max(b.peak_qubits);
        }
        let mean_ccx = total_ccx as f64 / samples as f64;
        let approx_35 = mean_ccx * 35.0;
        eprintln!(
            "scaled modular BY pair update sparse-cleanup: mean_ccx/window={mean_ccx:.1}, approx_35≈{approx_35:.0}, max_peak={max_peak}q"
        );
        assert!(approx_35 > 2_000_000.0, "scaled modular sparse cleanup unexpectedly competitive; revisit BY path");
    }

    fn emit_tagged_modular_microstep_for_cost(
        b: &mut super::super::B,
        r: &[super::super::QubitId],
        s: &[super::super::QubitId],
        a_ctrl: super::super::QubitId,
        b_ctrl: super::super::QubitId,
        p: U256,
    ) {
        // A: s -= r; r += s; r *= 2.  B: s += r; r *= 2.  C: r *= 2.
        super::super::cmod_add_qq(b, s, r, b_ctrl, p);
        super::super::cmod_sub_qq(b, s, r, a_ctrl, p);
        super::super::cmod_add_qq(b, r, s, a_ctrl, p);
        super::super::mod_double_inplace_fast(b, r, p);
    }

    #[test]
    fn hybrid_jump_denominator_with_microstep_tag_channel_still_too_costly() {
        // Valid hybrid after the dense-inverse correction: use jumped sparse
        // scaled updates only for the integer denominator pair, but update the
        // modular tagged channel by raw in-place BY microsteps to avoid dense
        // 2^-w inverse matrices. This is coherent and low-scratch, but the
        // modular microsteps dominate.
        const N: usize = 256;
        const W: usize = 16;
        const WIDTH: usize = N + W + 2;
        let p = SECP256K1_P;
        let approx_windows = 550usize.div_ceil(W);

        let mut hasher = sha3::Shake128::default();
        hasher.update(b"by-hybrid-den-jump-mod-micro-v1");
        let mut reader = hasher.finalize_xof();
        let mut buf = [0u8; 24];
        let samples = 24usize;
        let mut total_den_pair_ccx = 0usize;
        for _ in 0..samples {
            reader.read(&mut buf);
            let f_low = (u64::from_le_bytes(buf[0..8].try_into().unwrap()) as i128) | 1;
            let g_low = u64::from_le_bytes(buf[8..16].try_into().unwrap()) as i128;
            let delta = (u64::from_le_bytes(buf[16..24].try_into().unwrap()) % 41) as i64 - 20;
            let (_, _, _, m) = jump_matrix_direct_lowword(W, W, delta, f_low, g_low);
            let mut b = super::super::B::new();
            emit_scaled_pair_update_with_cleanup_for_cost(&mut b, m, WIDTH, W);
            total_den_pair_ccx += count_ccx(&b.ops);
        }
        let mean_den_pair_ccx = total_den_pair_ccx as f64 / samples as f64;

        let mut b = super::super::B::new();
        let a_ctrl = b.alloc_qubit();
        let b_ctrl = b.alloc_qubit();
        let r = b.alloc_qubits(N);
        let s = b.alloc_qubits(N);
        let start = b.ops.len();
        emit_tagged_modular_microstep_for_cost(&mut b, &r, &s, a_ctrl, b_ctrl, p);
        let mod_micro_ccx = count_ccx(&b.ops[start..]);

        let approx_total = mean_den_pair_ccx * approx_windows as f64 + mod_micro_ccx as f64 * 550.0;
        eprintln!(
            "BY hybrid denom-jump + tagged-micro budget: den_pair/window≈{mean_den_pair_ccx:.1}, mod_micro/step={mod_micro_ccx}, approx_total≈{approx_total:.0}"
        );
        assert!(approx_total > 1_800_000.0, "hybrid unexpectedly beats Kaliski; revisit implementation path");
    }

    #[test]
    fn modular_jump_inverse_cleanup_is_dense_dead_end() {
        // Correct an important over-optimism: scaled adjugate cleanup is sparse
        // for the INTEGER denominator pair because the update is P/2^w. The
        // modular coefficient/tagged channel is updated by P, whose inverse is
        // 2^-w * adj(P) mod p. The 2^-w factor makes the constants dense.
        // Therefore per-window modular row replacement cannot use sparse
        // adjugate cleanup; it needs either raw microsteps or a new structural
        // trick.
        const W: usize = 16;
        let p = SECP256K1_P;
        let inv_scale = two_inv_pow(p, W);
        let mut hasher = sha3::Shake128::default();
        hasher.update(b"by-modular-inverse-density-v1");
        let mut reader = hasher.finalize_xof();
        let mut buf = [0u8; 24];
        let samples = 2_000usize;
        let mut total_pop = 0usize;
        let mut min_pop = usize::MAX;
        for _ in 0..samples {
            reader.read(&mut buf);
            let f_low = (u64::from_le_bytes(buf[0..8].try_into().unwrap()) as i128) | 1;
            let g_low = u64::from_le_bytes(buf[8..16].try_into().unwrap()) as i128;
            let delta = (u64::from_le_bytes(buf[16..24].try_into().unwrap()) % 41) as i64 - 20;
            let (_, _, _, m) = jump_matrix_direct_lowword(W, W, delta, f_low, g_low);
            let s = if det_sign_pow2(m, W) >= 0 { 1i128 } else { -1i128 };
            let inv_entries = [s * m.m11, -s * m.m01, -s * m.m10, s * m.m00];
            let pop: usize = inv_entries
                .iter()
                .map(|&e| popcount_u256(mulm(signed_i128_mod_p(e, p), inv_scale, p)))
                .sum();
            total_pop += pop;
            min_pop = min_pop.min(pop);
        }
        let mean_pop = total_pop as f64 / samples as f64;
        eprintln!(
            "BY modular inverse cleanup density: mean_popcount_4entries={mean_pop:.1}, min_popcount_4entries={min_pop}"
        );
        assert!(mean_pop > 450.0, "modular inverse cleanup unexpectedly sparse");
    }

    #[test]
    fn optimistic_two_pair_integer_cleanup_lower_bound() {
        // Optimistic lower bound for the tagged-DIV shape if BOTH pairs could
        // use the sparse integer scaled-adjugate cleanup. Later tests show the
        // modular coefficient/tag pair cannot use this directly (unscaled
        // inverse is dense; scaled modular row formation is currently costly),
        // so this is a floor, not an implementation forecast.
        const N: usize = 256;
        const W: usize = 16;
        const WIDTH: usize = N + W + 2;
        const PAIRS: usize = 2;
        let exact_windows = safegcd_iters(N).div_ceil(W);
        let approx_windows = 550usize.div_ceil(W);

        let mut hasher = sha3::Shake128::default();
        hasher.update(b"by-tagged-div-two-pair-budget-v1");
        let mut reader = hasher.finalize_xof();
        let mut buf = [0u8; 24];
        let samples = 24usize;
        let mut total_pair_ccx = 0usize;
        let mut single_pair_peak = 0u32;
        for _ in 0..samples {
            reader.read(&mut buf);
            let f_low = (u64::from_le_bytes(buf[0..8].try_into().unwrap()) as i128) | 1;
            let g_low = u64::from_le_bytes(buf[8..16].try_into().unwrap()) as i128;
            let delta = (u64::from_le_bytes(buf[16..24].try_into().unwrap()) % 41) as i64 - 20;
            let (_, _, _, m) = jump_matrix_direct_lowword(W, W, delta, f_low, g_low);
            let mut b = super::super::B::new();
            emit_scaled_pair_update_with_cleanup_for_cost(&mut b, m, WIDTH, W);
            total_pair_ccx += count_ccx(&b.ops);
            single_pair_peak = single_pair_peak.max(b.peak_qubits);
        }
        let mean_pair_ccx = total_pair_ccx as f64 / samples as f64;
        let exact_ccx = mean_pair_ccx * PAIRS as f64 * exact_windows as f64;
        let approx_ccx = mean_pair_ccx * PAIRS as f64 * approx_windows as f64;
        let other_persistent_pair = 2 * WIDTH;
        let lowword_control = 2 * W + 16;
        let scheduled_peak = single_pair_peak as usize + other_persistent_pair + lowword_control;
        let scratch_beyond_two_field_regs = scheduled_peak.saturating_sub(2 * N);
        eprintln!(
            "BY optimistic 2-pair integer-cleanup lower bound: width={WIDTH}, mean_pair_ccx={mean_pair_ccx:.1}, exact≈{exact_ccx:.0}, approx≈{approx_ccx:.0}, scheduled_peak≈{scheduled_peak}q, scratch_beyond_2n≈{scratch_beyond_two_field_regs}q"
        );
        assert!(approx_ccx < 600_000.0, "approx tagged-DIV BY budget not SOTA-shaped");
        assert!(scheduled_peak < 2_400, "two-pair BY tagged-DIV lower-bound peak too high");
    }

    #[test]
    fn jumpdivstep_full_state_cleanup_budget_exceeds_2800_after_carry_fix() {
        // Stronger model than row-only: use the measured replacement+cleanup
        // pair cost and schedule the three BY pairs sequentially. This is the
        // best current proxy for a real jumped-BY inversion before low-word
        // matrix synthesis is included. After fixing shifted-row carry slack,
        // this full 3-pair state no longer fits the current 2800q cap.
        const N: usize = 256;
        const W: usize = 16;
        const WIDTH: usize = N + W + 2;
        const PAIRS: usize = 3;
        let exact_windows = safegcd_iters(N).div_ceil(W);
        let approx_windows = 550usize.div_ceil(W);

        let mut hasher = sha3::Shake128::default();
        hasher.update(b"by-full-state-cleanup-budget-v1");
        let mut reader = hasher.finalize_xof();
        let mut buf = [0u8; 24];
        let samples = 24usize;
        let mut total_pair_ccx = 0usize;
        let mut single_pair_peak = 0u32;
        for _ in 0..samples {
            reader.read(&mut buf);
            let f_low = (u64::from_le_bytes(buf[0..8].try_into().unwrap()) as i128) | 1;
            let g_low = u64::from_le_bytes(buf[8..16].try_into().unwrap()) as i128;
            let delta = (u64::from_le_bytes(buf[16..24].try_into().unwrap()) % 41) as i64 - 20;
            let (_, _, _, m) = jump_matrix_direct_lowword(W, W, delta, f_low, g_low);
            let mut b = super::super::B::new();
            emit_scaled_pair_update_with_cleanup_for_cost(&mut b, m, WIDTH, W);
            total_pair_ccx += count_ccx(&b.ops);
            single_pair_peak = single_pair_peak.max(b.peak_qubits);
        }
        let mean_pair_ccx = total_pair_ccx as f64 / samples as f64;
        let exact_ccx = mean_pair_ccx * PAIRS as f64 * exact_windows as f64;
        let approx_ccx = mean_pair_ccx * PAIRS as f64 * approx_windows as f64;
        let other_persistent_pairs = (PAIRS - 1) * 2 * WIDTH;
        let lowword_control = 2 * W + 16;
        let scheduled_peak = single_pair_peak as usize + other_persistent_pairs + lowword_control;
        eprintln!(
            "BY full-state cleanup budget: width={WIDTH}, mean_pair_ccx={mean_pair_ccx:.1}, exact≈{exact_ccx:.0}, approx≈{approx_ccx:.0}, scheduled_peak≈{scheduled_peak}q"
        );
        assert!(exact_ccx < 1_250_000.0, "exact BY cleanup budget no longer competitive");
        assert!(scheduled_peak > 2_800, "3-pair BY cleanup unexpectedly fits again; revisit full inverse path");
    }

    #[test]
    fn jumpdivstep_full_state_budget_model() {
        // Ground-up BY jump inversion budget from the calibrated row-former.
        // State model for one inversion:
        //   (f,g) signed pair + two coefficient columns = 6 wide registers.
        // Row application is sequential with two shared output rows and one
        // Cuccaro carry strip. This is the first budget that includes both
        // Toffoli and qubits in the same model.
        const N: usize = 256;
        const W: usize = 16;
        const WIDTH: usize = N + W + 2;
        const PAIRS: usize = 3;
        let exact_windows = safegcd_iters(N).div_ceil(W);
        let approx_windows = 550usize.div_ceil(W);

        let mut hasher = sha3::Shake128::default();
        hasher.update(b"by-full-state-budget-v1");
        let mut reader = hasher.finalize_xof();
        let mut buf = [0u8; 24];
        let samples = 24usize;
        let mut total_pair_ccx = 0usize;
        for _ in 0..samples {
            reader.read(&mut buf);
            let f_low = (u64::from_le_bytes(buf[0..8].try_into().unwrap()) as i128) | 1;
            let g_low = u64::from_le_bytes(buf[8..16].try_into().unwrap()) as i128;
            let delta = (u64::from_le_bytes(buf[16..24].try_into().unwrap()) % 41) as i64 - 20;
            let (_, _, _, m) = jump_matrix_direct_lowword(W, W, delta, f_low, g_low);
            let mut b = super::super::B::new();
            emit_constant_matrix_apply_for_cost(&mut b, m, WIDTH);
            total_pair_ccx += count_ccx(&b.ops);
        }
        let mean_pair_ccx = total_pair_ccx as f64 / samples as f64;
        let exact_row_ccx = mean_pair_ccx * PAIRS as f64 * exact_windows as f64;
        let approx_row_ccx = mean_pair_ccx * PAIRS as f64 * approx_windows as f64;

        let persistent_state = PAIRS * 2 * WIDTH; // six wide registers.
        let shared_outputs = 2 * WIDTH;
        let carry_strip = WIDTH;
        let lowword_control = 2 * W + 16; // f_low,g_low,delta/misc rough allowance.
        let peak_model = persistent_state + shared_outputs + carry_strip + lowword_control;
        eprintln!(
            "BY full-state budget model: width={WIDTH}, mean_pair_ccx={mean_pair_ccx:.1}, exact_row≈{exact_row_ccx:.0}, approx_row≈{approx_row_ccx:.0}, peak_model≈{peak_model}q"
        );
        assert!(exact_row_ccx < 700_000.0, "exact BY row budget too high");
        assert!(peak_model < 2_800, "BY modeled peak exceeds current cap");
    }

    #[test]
    fn jumpdivstep_matrix_entry_survey_test() {
        let samples = 100_000;
        for &w in &[4usize, 8, 12, 16] {
            let stats = jump_matrix_entry_survey(b"jumpdivstep-matrix-seed-v1", samples, w);
            let mean_log2 = if stats.nonzero_entries == 0 {
                0.0
            } else {
                stats.sum_log2_entry_abs / (stats.nonzero_entries as f64)
            };
            eprintln!("=== jumpdivstep matrix-entry survey (w={}) ===", w);
            eprintln!("samples                 : {}", stats.samples);
            eprintln!("max |entry| observed    : {}", stats.max_entry_abs);
            eprintln!(
                "max log2 |entry|        : {:.3}",
                (stats.max_entry_abs as f64).log2()
            );
            eprintln!("mean log2 |entry|       : {:.3}", mean_log2);
            eprintln!("theoretical max log2    : {}", w);
            eprintln!("===========================================");
            assert!(
                stats.max_entry_abs <= (1i128 << w),
                "w={} entry {} exceeded 2^w",
                w,
                stats.max_entry_abs
            );
        }
    }

    #[test]
    fn jumpdivstep_matrix_histogram() {
        // New moonshot stress-test: even if entries hit 2^w, maybe the NUMBER
        // of distinct matrices is tiny, allowing a heavily-compressed QROM.
        // This keeps the moonshot alive only if strong collapse occurs.
        for &w in &[4usize, 6, 8] {
            let hist = jump_matrix_histogram_all_states(w);
            eprintln!("=== jumpdivstep matrix histogram (w={}) ===", w);
            eprintln!("samples              : {}", hist.samples);
            eprintln!("distinct matrices    : {}", hist.distinct_matrices);
            eprintln!("most common count    : {}", hist.most_common_count);
            eprintln!("unique singleton mats: {}", hist.total_unique_rows);
            if let Some(m) = hist.most_common_matrix {
                eprintln!(
                    "most common matrix   : [[{}, {}], [{}, {}]]",
                    m.m00, m.m01, m.m10, m.m11
                );
            }
            eprintln!(
                "compression factor   : {:.2}",
                hist.samples as f64 / hist.distinct_matrices as f64
            );
            eprintln!("============================================");
            assert!(hist.distinct_matrices >= 1);
        }
    }
}
