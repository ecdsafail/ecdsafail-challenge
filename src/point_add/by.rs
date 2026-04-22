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

use alloy_primitives::U256;
use sha3::digest::{ExtendableOutput, Update, XofReader};

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
        Self { neg: false, mag: U256::ZERO }
    }

    pub fn from_u(x: U256) -> Self {
        Self { neg: false, mag: x }
    }

    pub fn negate(self) -> Self {
        if self.mag.is_zero() {
            self
        } else {
            Self { neg: !self.neg, mag: self.mag }
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
            (false, false) => Self { neg: false, mag: a.mag.wrapping_add(b.mag) },
            (true, true) => Self { neg: true, mag: a.mag.wrapping_add(b.mag) },
            (false, true) => sub_mag(a.mag, b.mag),
            (true, false) => sub_mag(b.mag, a.mag),
        }
    }

    pub fn sub(a: Self, b: Self) -> Self {
        Self::add(a, b.negate())
    }

    pub fn shr1_even(self) -> Self {
        debug_assert!(!self.bit0(), "shr1_even on odd integer");
        Self { neg: self.neg, mag: self.mag >> 1 }
    }
}

fn sub_mag(a: U256, b: U256) -> SInt {
    if a >= b {
        SInt { neg: false, mag: a.wrapping_sub(b) }
    } else {
        SInt { neg: true, mag: b.wrapping_sub(a) }
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
    if borrow { r.wrapping_add(p) } else { r }
}

fn negm(a: U256, p: U256) -> U256 {
    if a.is_zero() { a } else { p.wrapping_sub(a) }
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
            coeffs = Coeffs { uu: nu, vv: nv, qq: nq, rr: nr };
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
            coeffs = Coeffs { uu: nu, vv: nv, qq: nq, rr: nr };
        } else {
            // Case C:
            //   (δ, f, g) ← (1 + δ, f, g / 2)
            //   (U,V,Q,R) ← (2U, 2V, Q, R)
            let ng = g.shr1_even();
            let nu = addm(coeffs.uu, coeffs.uu, p);
            let nv = addm(coeffs.vv, coeffs.vv, p);
            delta = 1 + delta;
            g = ng;
            coeffs = Coeffs { uu: nu, vv: nv, qq: coeffs.qq, rr: coeffs.rr };
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

/// Recover `g0^{-1} mod p` from a converged divsteps run.
///
/// From the invariant `2^k f_k = U p + V g0`, with final `f_k = ±1`:
///
///     g0^{-1} ≡ sign(f_k) · V · 2^{-k}  (mod p)
pub fn recover_modinv(run: &DivstepsRun, p: U256) -> Option<U256> {
    if !run.converged { return None; }
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
        Self { reader: Box::new(hasher.finalize_xof()), p }
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
        if self.samples == 0 { 0.0 } else { self.sum_iters as f64 / self.samples as f64 }
    }
}

pub fn survey(
    sampler: &mut Sampler,
    n_samples: usize,
    p: U256,
    max_iters: usize,
) -> SurveyStats {
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

    for _ in 0..n_samples {
        let x = sampler.next();
        let run = run_divsteps(x, p, max_iters);
        if !run.converged {
            stats.all_converged = false;
        }
        let k = run.iters_done;
        stats.samples += 1;
        if k < stats.min_iters { stats.min_iters = k; }
        if k > stats.max_iters { stats.max_iters = k; }
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
    if t == 0 { return 0; }
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
        TransitionMatrix { m00: u, m01: v, m10: q, m11: r, delta_final: delta },
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
            p00 = np00; p01 = np01; p10 = np10; p11 = np11;
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
            p00 = np00; p01 = np01; p10 = np10; p11 = np11;
        } else {
            // Case C
            let (np00, np01, np10, np11) = (
                2 * p00,
                2 * p01,
                p10,
                p11,
            );
            let new_g = g / 2;
            delta = 1 + delta;
            g = new_g;
            p00 = np00; p01 = np01; p10 = np10; p11 = np11;
        }
        n -= 1;
        t -= 1;
        g = truncate_i128(g, t);
    }
    let f_out = truncate_i128(f, t + 1); // after n=w steps, f known to t-w+1 bits
    let g_out = truncate_i128(g, t);     // and g to t-w bits. Here `t` already decremented.
    (delta, f_out, g_out, TransitionMatrix {
        m00: p00,
        m01: p01,
        m10: p10,
        m11: p11,
        delta_final: delta,
    })
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
    let mut buf = [0u8; 24];
    for _ in 0..n_samples {
        reader.read(&mut buf);
        let mut f_low = u64::from_le_bytes(buf[0..8].try_into().unwrap()) as i128;
        let g_low = u64::from_le_bytes(buf[8..16].try_into().unwrap()) as i128;
        let delta = (u64::from_le_bytes(buf[16..24].try_into().unwrap()) % 41) as i64 - 20;
        f_low |= 1; // ensure odd
        let (_, _, _, m) = jump_matrix_direct_lowword(w, w, delta, f_low, g_low);
        for &e in &[m.m00, m.m01, m.m10, m.m11] {
            let abs = e.wrapping_abs();
            if abs > stats.max_entry_abs { stats.max_entry_abs = abs; }
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
        eprintln!("max log2 |entry|        : {:.3}", (js.max_entry_abs as f64).log2());
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
        assert!(stats.max_iters <= theoretical_bound,
            "observed max iters {} exceeds theoretical bound {}",
            stats.max_iters, theoretical_bound);
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
            eprintln!("max log2 |entry|        : {:.3}", (stats.max_entry_abs as f64).log2());
            eprintln!("mean log2 |entry|       : {:.3}", mean_log2);
            eprintln!("theoretical max log2    : {}", w);
            eprintln!("===========================================");
            assert!(stats.max_entry_abs <= (1i128 << w),
                "w={} entry {} exceeded 2^w", w, stats.max_entry_abs);
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
                eprintln!("most common matrix   : [[{}, {}], [{}, {}]]", m.m00, m.m01, m.m10, m.m11);
            }
            eprintln!("compression factor   : {:.2}", hist.samples as f64 / hist.distinct_matrices as f64);
            eprintln!("============================================");
            assert!(hist.distinct_matrices >= 1);
        }
    }
}
