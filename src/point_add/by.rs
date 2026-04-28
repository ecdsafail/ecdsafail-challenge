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

    fn ratio_window_next(w: usize, delta: i64, h: u64) -> (i64, u64) {
        let signed_h = if (h & (1u64 << (w - 1))) != 0 {
            (h as i128) - (1i128 << w)
        } else {
            h as i128
        };
        let mut d = delta;
        let mut f = SInt::from_u(U256::from(1));
        let mag = U256::from(signed_h.unsigned_abs());
        let mut g = SInt { neg: signed_h < 0, mag };
        for _ in 0..w {
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
        (d, g_low.wrapping_mul(inv_f) & mask)
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

    fn emit_signed_row_scaled_from_sources_for_test(
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
        let sgn = det_sign_pow2(mtx, 16);
        let m0 = b.alloc_qubits(16);
        let m1 = b.alloc_qubits(16);
        compute_row_correction_m_from_sources(b, mtx.m00, x0, mtx.m01, x1, &m0, false);
        compute_row_correction_m_from_sources(b, mtx.m10, x0, mtx.m11, x1, &m1, false);
        let (q0, q1) = compute_signed_q_from_m_for_matrix(b, mtx, &m0, &m1);
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
        subtract_signed_q_times_solinas_c_for_cost(b, &q0, x0);
        subtract_signed_q_times_solinas_c_for_cost(b, &q1, x1);

        // Clear m using P*q = m (mod 2^16).
        add_low_coeff_mod16_for_cost(b, mtx.m00.rem_euclid(1 << 16) as u64, &q0, &m0, true);
        add_low_coeff_mod16_for_cost(b, mtx.m01.rem_euclid(1 << 16) as u64, &q1, &m0, true);
        add_low_coeff_mod16_for_cost(b, mtx.m10.rem_euclid(1 << 16) as u64, &q0, &m1, true);
        add_low_coeff_mod16_for_cost(b, mtx.m11.rem_euclid(1 << 16) as u64, &q1, &m1, true);

        clear_signed_q_from_z_high_for_cost(b, &q0, &z0);
        clear_signed_q_from_z_high_for_cost(b, &q1, &z1);

        add_signed_coeff_times_for_cost(b, -sgn * mtx.m11, y0, &z0);
        add_signed_coeff_times_for_cost(b, sgn * mtx.m01, y1, &z0);
        add_signed_coeff_times_for_cost(b, sgn * mtx.m10, y0, &z1);
        add_signed_coeff_times_for_cost(b, -sgn * mtx.m00, y1, &z1);
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
