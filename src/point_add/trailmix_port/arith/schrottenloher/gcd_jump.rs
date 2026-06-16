//! "Jump" GCD dialog: Stein-style binary GCD that removes up to `jump`
//! trailing zeros of `v` per step (vs the divstep convention of exactly 1).
//! Fewer steps => the per-step adders (csub/cswap in the GCD, the
//! controlled mod-add in `apply_bv`) fire fewer times => big Toffoli win, at
//! ~unchanged dialog length (the recorded shift count costs ~1 bit/step,
//! which the step reduction offsets). See `scripts/schr_aggressive.py`.
//!
//! This module currently holds the CLASSICAL reference + roundtrip tests
//! that validate the algebra (forward GCD -> per-step dialog -> apply
//! reconstructs z*x mod q). The quantum forward-GCD / `apply_bv` / reverse
//! build on top once the algebra is confirmed.
//!
//! Per-step dialog symbol = (b0, `b0_and_b1`, j):
//!   b0        = v was odd (a subtract happened)
//!   `b0_and_b1` = a swap happened (only when b0)
//!   j         = number of trailing zeros removed this step, in 1..=jump
//!
//! Reconstruction mirror (apply, reverse iter order): v = 2^j * v mod q;
//! if b0: v += u mod q; if `b0_and_b1`: swap(u,v).

use crate::point_add::trailmix_port::compat::num_bigint::BigUint;
use crate::point_add::trailmix_port::compat::num_traits::{One, Zero};

/// One recorded dialog step.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct JumpSym {
    pub b0: u8,
    pub b0_and_b1: u8,
    pub j: u8,
}

/// Forward jump-GCD on `(u, v)` (u must be odd, e.g. the modulus q). Runs
/// until `(u, v) == (1, 0)` or `cap` steps. Returns the per-step dialog and
/// whether it converged. `jump >= 1`; jump=1 reproduces the divstep dialog.
#[must_use]
pub fn to_bitvector_classical_jump(
    mut u: BigUint,
    mut v: BigUint,
    jump: usize,
    cap: usize,
) -> (Vec<JumpSym>, bool) {
    assert!(u.bit(0), "u must be odd");
    assert!(jump >= 1);
    let mut dialog = Vec::with_capacity(cap);
    let one = BigUint::one();
    for _ in 0..cap {
        let b0 = u8::from(v.bit(0));
        let mut b0_and_b1 = 0u8;
        if b0 == 1 {
            if u > v {
                std::mem::swap(&mut u, &mut v);
                b0_and_b1 = 1;
            }
            // v >= u now; v - u is even.
            v -= &u;
        }
        // Remove up to `jump` trailing zeros of v (v is even here unless 0).
        let mut j = 0u8;
        while (j as usize) < jump && !v.is_zero() && !v.bit(0) {
            v >>= 1u32;
            j += 1;
        }
        dialog.push(JumpSym { b0, b0_and_b1, j });
        if u == one && v.is_zero() {
            return (dialog, true);
        }
    }
    (dialog, u == one && v.is_zero())
}

/// Replay the jump dialog (reverse iter order): the inverse of the forward
/// reduction, reconstructing the modular product. On input `(u, v) = (z, 0)`
/// and a dialog from `to_bitvector_classical_jump(q, x, ...)`, returns
/// `(0, z * x mod q)` — exactly the `IPModMul` contract, generalized to jumps.
#[must_use]
pub fn apply_bitvector_classical_jump(
    mut u: BigUint,
    mut v: BigUint,
    dialog: &[JumpSym],
    q: &BigUint,
) -> (BigUint, BigUint) {
    for s in dialog.iter().rev() {
        for _ in 0..s.j {
            v = (&v * 2u32) % q; // v *= 2^j mod q
        }
        if s.b0 == 1 {
            v = (&v + &u) % q;
        }
        if s.b0_and_b1 == 1 {
            std::mem::swap(&mut u, &mut v);
        }
    }
    (u, v)
}

// ---------------------------------------------------------------------------
// Jump-BEFORE-swap variant (cleaner dialog alphabet).
//
// Reorder the step so the trailing-zero removal happens FIRST: shift v to odd
// (up to `jump`), THEN conditionally swap + subtract. Because v is even at the
// start of every step (post-subtract), there is no independent parity flag --
// the symbol decouples into INDEPENDENT (j = zeros removed, swap). For jump=2
// the per-step alphabet is exactly 5: (j=1,swap), (j=2,swap), and the overflow
// (ctz>jump -> v still even -> no subtract, so swap is absent). 5^3=125 < 2^7,
// so 3 symbols pack into 7 bits (~2.33 b/step) vs base-6's 2.585. This avoids
// the cross-step b0<->s_2 correlation of the swap-before-jump ordering.
// ---------------------------------------------------------------------------

/// One jump-before-swap dialog step. `subtracted=0` is the overflow symbol
/// (ctz>jump: v stayed even, no swap, no subtract).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct JumpBsSym {
    pub j: u8,          // trailing zeros removed this step (0..=jump)
    pub swap: u8,       // swap bit (only meaningful when subtracted==1)
    pub subtracted: u8, // 1 = v became odd after the shift -> subtract happened
}

/// Forward jump-before-swap GCD on `(u, v)` (u odd). Per step: remove up to
/// `jump` trailing zeros of v; if v is now odd, (swap if u>v) and v -= u; else
/// (ctz>jump) record the overflow symbol. Converges to (1, 0).
#[must_use]
pub fn to_bitvector_classical_jump_bs(
    mut u: BigUint,
    mut v: BigUint,
    jump: usize,
    cap: usize,
) -> (Vec<JumpBsSym>, bool) {
    assert!(u.bit(0), "u must be odd");
    assert!(jump >= 1);
    let one = BigUint::one();
    let mut dialog = Vec::with_capacity(cap);
    for _ in 0..cap {
        // Shift up to `jump` trailing zeros (v never reaches 0 by shifting: a
        // power of two shifts down to the odd 1).
        let mut j = 0u8;
        while (j as usize) < jump && !v.bit(0) {
            v >>= 1u32;
            j += 1;
        }
        if v.bit(0) {
            // v odd -> swap if u>v, subtract.
            let mut swap = 0u8;
            if u > v {
                std::mem::swap(&mut u, &mut v);
                swap = 1;
            }
            v -= &u;
            dialog.push(JumpBsSym {
                j,
                swap,
                subtracted: 1,
            });
        } else {
            // ctz > jump: v still even, no subtract this step (overflow).
            dialog.push(JumpBsSym {
                j,
                swap: 0,
                subtracted: 0,
            });
        }
        if u == one && v.is_zero() {
            return (dialog, true);
        }
    }
    (dialog, u == one && v.is_zero())
}

/// Replay a jump-before-swap dialog (reverse iter order). Inverse of the
/// forward step (shift; swap; sub) is (v += u; swap; v *= 2^j). On input
/// (u, v) = (z, 0) with a dialog from `to_bitvector_classical_jump_bs(q, x)`,
/// reconstructs the modular product in v.
#[must_use]
pub fn apply_bitvector_classical_jump_bs(
    mut u: BigUint,
    mut v: BigUint,
    dialog: &[JumpBsSym],
    q: &BigUint,
) -> (BigUint, BigUint) {
    for s in dialog.iter().rev() {
        if s.subtracted == 1 {
            v = (&v + &u) % q;
            if s.swap == 1 {
                std::mem::swap(&mut u, &mut v);
            }
        }
        for _ in 0..s.j {
            v = (&v * 2u32) % q;
        }
    }
    (u, v)
}

/// Monte-Carlo register-shrink schedule + iteration budget for the jump-GCD.
/// Instead of the closed-form `ceil(N - i*slope + pad)`, simulate `samples`
/// random inversions and, at each iteration, measure the bit-width actually
/// needed (`max(bitlen(u), bitlen(v))`). Set `current_n[i] = mean_i +
/// margin_sigma * std_i` (per-iter margin, tight where the variance is small).
/// The iteration budget is the max convergence step over the samples plus a
/// small slack. Returns (schedule, `iters_budget`, `sum_of_current_n`).
///
/// `sum_of_current_n` is the proxy for the GCD csub/cswap Toffoli (those
/// scale with the per-iter working width), so a tighter MC schedule directly
/// lowers Toffoli vs the closed form.
#[must_use]
pub fn mc_schedule(
    q: &BigUint,
    jump: usize,
    samples: usize,
    margin_sigma: f64,
    cap: usize,
) -> (Vec<usize>, usize, usize) {
    use rand::{thread_rng, Rng};
    let one = BigUint::one();
    let mut sum = vec![0f64; cap];
    let mut sumsq = vec![0f64; cap];
    let mut conv_max = 0usize;
    let mut rng = thread_rng();
    for _ in 0..samples {
        let x = (BigUint::from_bytes_le(&rng.gen::<[u8; 32]>()) % (q - 1u32)) + 1u32;
        let (mut u, mut v) = (q.clone(), x);
        let mut conv = cap;
        for i in 0..cap {
            // width needed at the START of iter i (before truncation).
            let w = u.bits().max(v.bits()) as f64;
            sum[i] += w;
            sumsq[i] += w * w;
            // one jump step
            if v.bit(0) {
                if u > v {
                    std::mem::swap(&mut u, &mut v);
                }
                v -= &u;
            }
            let mut j = 0;
            while j < jump && !v.is_zero() && !v.bit(0) {
                v >>= 1u32;
                j += 1;
            }
            if u == one && v.is_zero() {
                conv = i + 1;
                break;
            }
        }
        conv_max = conv_max.max(conv);
    }
    let iters_budget = conv_max + 5; // small slack above the observed max
    let mut schedule = Vec::with_capacity(iters_budget);
    let mut total = 0usize;
    for i in 0..iters_budget {
        let mean = sum[i] / samples as f64;
        let var = (sumsq[i] / samples as f64 - mean * mean).max(0.0);
        let cn = (mean + margin_sigma * var.sqrt()).ceil() as usize;
        let cn = cn.clamp(1, 256);
        schedule.push(cn);
        total += cn;
    }
    (schedule, iters_budget, total)
}

// ---------------------------------------------------------------------------
// Base-5 dialog packing (jump=2): the 5 steady jump-before-swap symbols pack
// 3-per-7-bits (5^3 = 125 < 2^7 = 128, 97.7% efficient), ~2.33 bits/step vs
// the 6-symbol base-6's 2.585. Classical reference for the quantum packer.
// ---------------------------------------------------------------------------

/// Map a jump-before-swap symbol `(subtracted, swap, s_2)` to a base-5 digit
/// in 0..=4. The 5 steady (jump=2) symbols: (1,0,0)->0, (1,1,0)->1, (1,0,1)->2,
/// (1,1,1)->3, (0,0,1)->4 (overflow). The step-0 boundary symbol (0,1,1) and
/// the invalid combos are not in the steady alphabet (handled by the t1 prefix).
#[must_use]
pub fn jump_bs_sym_to_digit(subtracted: u8, swap: u8, s2: u8) -> u8 {
    if subtracted == 1 {
        swap + 2 * s2
    } else {
        4
    }
}

/// Inverse of [`jump_bs_sym_to_digit`].
#[must_use]
pub fn jump_bs_digit_to_sym(d: u8) -> (u8, u8, u8) {
    match d {
        0 => (1, 0, 0),
        1 => (1, 1, 0),
        2 => (1, 0, 1),
        3 => (1, 1, 1),
        4 => (0, 0, 1),
        _ => panic!("jump_bs digit {d} out of range 0..=4"),
    }
}

/// Pack 3 base-5 digits into a 7-bit value `d0 + 5*d1 + 25*d2` in 0..125.
#[must_use]
pub fn jump_bs_pack3(d: [u8; 3]) -> u8 {
    d[0] + 5 * d[1] + 25 * d[2]
}

/// Inverse of [`jump_bs_pack3`].
#[must_use]
pub fn jump_bs_unpack3(v: u8) -> [u8; 3] {
    [v % 5, (v / 5) % 5, v / 25]
}

// ---------------------------------------------------------------------------
// Quantum jump-GCD (RAW dialog; base-5 packing is the qubit-reduction pass).
// ---------------------------------------------------------------------------

use crate::point_add::trailmix_port::circuit::{Circuit, QReg};

/// Controlled right-shift-by-1 (LSB-first logical shift): when `ctrl` is set,
/// moves v[i+1] -> v[i] (v[0], known |0> in the GCD, wraps to the top). Gated
/// analogue of `right_shift`, each swap a Fredkin. Kept LOCAL to this module
/// to stay decoupled from the in-flight arith/shift.rs reorg.
fn controlled_right_shift(circ: &mut Circuit, ctrl: &QReg, v: &[QReg]) {
    for i in 0..v.len().saturating_sub(1) {
        circ.cswap(ctrl, &v[i], &v[i + 1]);
    }
}

/// Exact gate-inverse of `controlled_right_shift` (cswaps in reverse order).
fn controlled_left_shift(circ: &mut Circuit, ctrl: &QReg, v: &[QReg]) {
    for i in (1..v.len()).rev() {
        circ.cswap(ctrl, &v[i], &v[i - 1]);
    }
}

/// Controlled pseudo-Mersenne modular doubling (secp256k1): if `ctrl`,
/// a := 2*a mod q; else a unchanged. `a.len() == n+1`, `a[n] = |0>`.
///
/// Same shape as `pm_prims::mod_double_pm` but gating ONLY the shift: after a
/// controlled left-shift, `a[n] = (ctrl AND overflow)`, so the borrowed-dirty
/// +f window stays unconditional on `a[n]` (adds nothing when ctrl=0). Only
/// the ancilla cleanup needs the control (`a[0] == a[n]` post-add holds iff
/// ctrl shifted), so the final clean is a `ccx`. Inlined here (public
/// `controlled_add_const_gidney`) to avoid editing the in-flight `pm_prims.rs`.
fn controlled_mod_double_secp256k1(circ: &mut Circuit, ctrl: &QReg, a: &[QReg]) {
    const LSBS: usize = 30 + 33; // padding + f_bitlen = 63
    let dh = LSBS + (LSBS - 1); // dirty-borrow window: a[LSBS..125]
    let top = a.len() - 1;
    let f_bytes = super::F_SECP256K1.to_le_bytes();
    controlled_left_shift(circ, ctrl, a);
    crate::point_add::trailmix_port::arith::gidney_const_adder::controlled_add_const_gidney(
        circ,
        &a[top],
        &a[..LSBS],
        &f_bytes,
        &a[LSBS..dh],
    );
    circ.ccx(ctrl, &a[0], &a[top]);
}

/// Exact gate-inverse of `controlled_mod_double_secp256k1`: if `ctrl`,
/// a := a/2 mod q (doubling-representation halve); else unchanged.
fn controlled_mod_double_secp256k1_reverse(circ: &mut Circuit, ctrl: &QReg, a: &[QReg]) {
    const LSBS: usize = 30 + 33;
    let dh = LSBS + (LSBS - 1);
    let top = a.len() - 1;
    let f_bytes = super::F_SECP256K1.to_le_bytes();
    circ.ccx(ctrl, &a[0], &a[top]);
    // -f = X-sandwich of the controlled +f.
    for q in &a[..LSBS] {
        circ.x(q);
    }
    crate::point_add::trailmix_port::arith::gidney_const_adder::controlled_add_const_gidney(
        circ,
        &a[top],
        &a[..LSBS],
        &f_bytes,
        &a[LSBS..dh],
    );
    for q in &a[..LSBS] {
        circ.x(q);
    }
    controlled_right_shift(circ, ctrl, a);
}

/// Venting budget for the GCD csub/cadd. In the FULL circuit the GCD runs
/// below the `apply_bv` qubit peak, so full venting (`usize::MAX` -> Gidney 2n)
/// is free; the standalone roundtrip test has no `apply_bv` to dominate, so a
/// vent-free Cuccaro (0) keeps that test's peak honest. Forward and reverse
/// must agree.
const JUMP_GCD_VENTS: usize = usize::MAX;

/// Top-bits window for the swap-decision comparator. The schedule is the EXTREME
/// width, so a fast shot's value sits below `current_n`; the swap is decided by
/// the highest bit where u,v differ, which is at most `gap` bits below
/// `current_n`. Measured (`gen_jump_schedule` `[cmp]` sweep, 40M held-out shots):
/// max gap = 58, P(gap>=64)=0/40M, P(gap>=56)=1.25e-7. trunc=72 gives a 14-bit
/// margin over the observed max => per-pass mis-swap ~1e-12, whole-circuit
/// reliability ~1.0. The precision contract `jump_gcd_subtract_precond` guards
/// every shot against a tie mis-swap. Forward and reverse must agree.
const JUMP_CMP_TRUNC: usize = 72;

/// Conservative iteration budget for the quantum jump-GCD: classical max
/// over 2000 random inversions + slack. Forward and reverse MUST agree, so
/// this is a pure function of `jump`. Refined later via the MC schedule.
#[must_use]
pub fn jump_iters_budget(jump: usize) -> usize {
    super::jump_schedule::jump_schedule(jump).1
}

/// q = 2^256 - `F_SECP256K1` as a `BigUint` (the GCD's `u_init`).
fn q_secp256k1() -> BigUint {
    (BigUint::one() << 256u32) - BigUint::from(super::F_SECP256K1)
}

/// Quantum forward jump-GCD on (u=q, v=`v_full`). Phase-A validation build:
/// full-width 256 recurrence (no register shrink). Per step the variable
/// shift removes up to `jump` trailing zeros of v via nested controlled
/// rotations, recording shift flags `s_2..s_jump`. The per-step dialog symbol
/// is stored RAW as (jump+1) bits: (b0, `b0_and_b1`, `s_2`, ..., `s_jump`). Returns
/// the raw garbage tape (caller frees it via the reverse pass).
pub fn forward_gcd_jump_quantum_secp256k1(
    circ: &mut Circuit,
    v_full: &mut Vec<QReg>,
    jump: usize,
) -> Vec<QReg> {
    use crate::point_add::trailmix_port::arith::poc_arith;
    use crate::point_add::trailmix_port::arith::schrottenloher::msb_compare;

    let n = 256usize;
    assert!(jump >= 1);
    assert!(v_full.len() >= n, "v_full must be at least n=256 bits");
    let (sched, iters) = super::jump_schedule::jump_schedule(jump);
    let sym_bits = jump + 1; // subtracted, swap, then s_2..s_jump
    let garbage_len = iters * sym_bits + 1; // +1: garbage[0] = step-0 shift1 flag (t1)

    let prev = circ.push_section("fwd_gcd_jump");

    // u_full = q (n bits); shrinks per the schedule (tail bits return to the
    // allocator pool for dialog reuse once the value no longer needs them).
    let mut u_full = circ.alloc_qreg_bits("uj_full", n);
    let q_bytes = q_secp256k1().to_bytes_le();
    for i in 0..n {
        if (q_bytes.get(i / 8).copied().unwrap_or(0) >> (i % 8)) & 1 == 1 {
            circ.x(&u_full[i]);
        }
    }

    let subtracted = circ.alloc_qreg("gcdj.subtracted"); // post-shift parity (1 => subtract)
    let swap_flag = circ.alloc_qreg("gcdj.swap");
    // shift flags s_2..s_jump (jump-1 of them), reused across steps.
    let s_flags: Vec<QReg> = (2..=jump).map(|_| circ.alloc_qreg("gcdj.s")).collect();
    // Step-0 shift1-fired flag (= x even); used only at i==0, inserted at
    // garbage[0] after the loop so apply/reverse find it at the front.
    let t1 = circ.alloc_qreg("gcdj.t1");

    let mut garbage_vec: Vec<QReg> = Vec::with_capacity(garbage_len);

    for i in 0..iters {
        let iter_prev = circ.push_section(&format!("fwdj_iter_{i:04}"));

        // Shrink u_full/v_full to the scheduled width: free tail bits that are
        // |0> for schedule-fitting inputs (a non-fitting tail input would have
        // a nonzero high bit here -> zero_and_free fails; that is the schedule
        // miss rate the generator bounds).
        let current_n = (sched[i] as usize).max(1);
        while u_full.len() > current_n {
            let q = u_full.pop().expect("u_full nonempty");
            circ.zero_and_free(q);
        }
        while v_full.len() > current_n {
            let q = v_full.pop().expect("v_full nonempty");
            circ.zero_and_free(q);
        }
        // Narrow swap-decision comparator: scan only the top JUMP_CMP_TRUNC bits
        // of the active width. A naive top-k window can sit ABOVE a fast shot's
        // differing bit (the EXTREME schedule makes current_n >> the shot's true
        // width) -> tie -> mis-swap; JUMP_CMP_TRUNC is sized 14 bits above the
        // measured max gap (58), so the window always covers the differing bit
        // (mis-swap ~1e-12). The forward pass's precision contract guards every
        // shot; forward and reverse use the identical window.
        let cmp_lo = current_n.saturating_sub(JUMP_CMP_TRUNC);
        let cmp_eff = current_n - cmp_lo;

        // 1) SHIFT-FIRST: remove up to `jump` trailing zeros of v. Steps>=1
        //    start with v even (post-subtract) so shift1 is unconditional (free
        //    swaps); step 0 gates shift1 on (v even) and records it in `t1`.
        let s_sh = circ.push_section("jshift");
        if i == 0 {
            circ.cx(&v_full[0], &t1); // t1 = v[0]
            circ.x(&t1); // t1 = NOT(v[0]) = (x even)
            controlled_right_shift(circ, &t1, &v_full[..current_n]);
        } else {
            poc_arith::right_shift(circ, &v_full[..current_n]);
        }
        // s_2..s_jump: shift again while still even (nesting is automatic --
        // once an odd bit reaches LSB the shift stops, so s_k goes 0 after).
        for s_k in &s_flags {
            circ.cx(&v_full[0], s_k);
            circ.x(s_k);
            controlled_right_shift(circ, s_k, &v_full[..current_n]);
        }
        circ.pop_section(&s_sh);
        // 2) subtracted = v[0] (post-shift parity): 1 => v odd => swap+subtract;
        //    0 => ctz>jump overflow (v still even, no subtract this step).
        circ.cx(&v_full[0], &subtracted);
        // 3) swap_flag ^= subtracted AND (v < u) on the active width.
        msb_compare::controlled_lt_msbs_gidney(
            circ,
            &subtracted,
            &v_full[cmp_lo..current_n],
            &u_full[cmp_lo..current_n],
            cmp_eff,
            &swap_flag,
        );
        // 4) cswap(swap_flag, u, v)
        let s_cs = circ.push_section("cswap");
        for j in 0..current_n {
            circ.cswap(&swap_flag, &u_full[j], &v_full[j]);
        }
        circ.pop_section(&s_cs);
        // Precision contract: post-swap, subtracted=1 must leave v >= u (the
        // v-=u precondition). Exact comparator => always holds; this guards the
        // JC narrow-window variant against precision-tie mis-swaps.
        circ.contract_check("jump_gcd_subtract_precond", |view, shot| {
            if !view.read_bit_shot(&subtracted, shot) {
                return Ok(());
            }
            let uu = view.read_u256_shot(&u_full[..current_n], shot);
            let vv = view.read_u256_shot(&v_full[..current_n], shot);
            if vv < uu {
                return Err(format!("v < u post-swap (approx mis-swap): v={vv} u={uu}"));
            }
            Ok(())
        });
        // 5) v -= subtracted * u  (X-sandwich two's-complement subtract)
        let s_sub = circ.push_section("csub");
        for q in &v_full[..current_n] {
            circ.x(q);
        }
        crate::point_add::trailmix_port::arith::gidney_const_adder::controlled_hybrid_add(
            circ,
            &subtracted,
            &v_full[..current_n],
            &u_full[..current_n],
            JUMP_GCD_VENTS,
        );
        for q in &v_full[..current_n] {
            circ.x(q);
        }
        circ.pop_section(&s_sub);
        // 6) record the symbol (subtracted, swap, s_2..s_jump) into fresh |0>
        //    slots (returning the ancilla to |0> for the next step).
        let mut slots: Vec<QReg> = (0..sym_bits).map(|_| circ.alloc_qreg("jdlg")).collect();
        circ.swap(&subtracted, &slots[0]);
        circ.swap(&swap_flag, &slots[1]);
        for (idx, s_k) in s_flags.iter().enumerate() {
            circ.swap(s_k, &slots[2 + idx]);
        }
        garbage_vec.append(&mut slots);
        // 7) strict-dealloc touch: subtracted/swap/s_k are |0> now (swapped
        //    out), so these CX gates are control-0 no-ops that advance
        //    u_full's touch edge past the dialog/comparator allocs.
        for k in 0..current_n {
            circ.cx(&subtracted, &u_full[k]);
            circ.cx(&swap_flag, &u_full[k]);
        }
        for s_k in &s_flags {
            circ.cx(s_k, &u_full[0]);
        }
        circ.pop_section(&iter_prev);
    }

    // Convergence contract: the binary GCD must reach (u, v) == (1, 0) on every
    // shot. Catch divergence HERE at the GCD boundary -- not at the late free of
    // u_full. u_full/v_full are now shrunk to the last scheduled width (~1 bit).
    let uw = u_full.len();
    let vw = v_full.len();
    circ.contract_check("jump_gcd_converged", |view, shot| {
        let u = view.read_u256_shot(&u_full[..uw], shot);
        let v = view.read_u256_shot(&v_full[..vw], shot);
        if u != BigUint::one() {
            return Err(format!("u={u} != 1 (GCD did not converge)"));
        }
        if !v.is_zero() {
            return Err(format!("v={v} != 0 (GCD did not converge)"));
        }
        Ok(())
    });

    // u converged to 1; X to 0 then free. Shrink v_full to one |0> bit so the
    // caller's register is not dead weight during the apply_bv peak (the
    // reverse regrows it symmetrically).
    circ.x(&u_full[0]);
    while v_full.len() > 1 {
        let q = v_full.pop().expect("v_full nonempty");
        circ.zero_and_free(q);
    }
    drop(subtracted);
    drop(swap_flag);
    drop(s_flags);
    drop(u_full);
    // Move t1 to the front: garbage = [t1, sym_0, ..., sym_{iters-1}].
    garbage_vec.insert(0, t1);
    assert_eq!(garbage_vec.len(), garbage_len);
    circ.pop_section(&prev);
    garbage_vec
}

/// Reverse of `forward_gcd_jump_quantum_secp256k1`: restores `v_full[..256]`
/// to the original `x` and drains `garbage` to all-|0>. Mirrors the forward
/// ops in reverse order with each step inverted.
pub fn forward_gcd_jump_quantum_secp256k1_reverse(
    circ: &mut Circuit,
    v_full: &mut Vec<QReg>,
    garbage: &mut Vec<QReg>,
    jump: usize,
) {
    use crate::point_add::trailmix_port::arith::poc_arith;
    use crate::point_add::trailmix_port::arith::schrottenloher::msb_compare;

    let n = 256usize;
    assert!(jump >= 1);
    let (sched, iters) = super::jump_schedule::jump_schedule(jump);
    let sym_bits = jump + 1;
    assert!(garbage.len() > iters * sym_bits); // +1 = t1 prefix at garbage[0]

    let prev = circ.push_section("fwd_gcd_jump_rev");

    // u_full regrows from 1 bit, symmetric to the forward's shrink to ~1.
    let mut u_full = circ.alloc_qreg_bits("uj_full_re", 1);
    // forward ended u_full=0 (post-X); reverse re-inits u_final=1 via X.
    circ.x(&u_full[0]);

    let subtracted = circ.alloc_qreg("gcdj_re.subtracted");
    let swap_flag = circ.alloc_qreg("gcdj_re.swap");
    let s_flags: Vec<QReg> = (2..=jump).map(|_| circ.alloc_qreg("gcdj_re.s")).collect();

    for i in (0..iters).rev() {
        let iter_prev = circ.push_section(&format!("revj_iter_{i:04}"));

        // Grow u_full/v_full back to step i's scheduled width (fresh |0> bits,
        // symmetric to the forward's shrink).
        let current_n = (sched[i] as usize).max(1);
        while u_full.len() < current_n {
            u_full.push(circ.alloc_qreg("uj_full_re"));
        }
        while v_full.len() < current_n {
            v_full.push(circ.alloc_qreg("vj_full_re"));
        }
        // Narrow swap-decision comparator: scan only the top JUMP_CMP_TRUNC bits
        // of the active width. A naive top-k window can sit ABOVE a fast shot's
        // differing bit (the EXTREME schedule makes current_n >> the shot's true
        // width) -> tie -> mis-swap; JUMP_CMP_TRUNC is sized 14 bits above the
        // measured max gap (58), so the window always covers the differing bit
        // (mis-swap ~1e-12). The forward pass's precision contract guards every
        // shot; forward and reverse use the identical window.
        let cmp_lo = current_n.saturating_sub(JUMP_CMP_TRUNC);
        let cmp_eff = current_n - cmp_lo;

        // Pull sym_i off the tape end into the ancilla (subtracted, swap, s_k).
        let glen = garbage.len();
        let cur: Vec<QReg> = garbage.split_off(glen - sym_bits);
        circ.swap(&subtracted, &cur[0]);
        circ.swap(&swap_flag, &cur[1]);
        for (idx, s_k) in s_flags.iter().enumerate() {
            circ.swap(s_k, &cur[2 + idx]);
        }

        // Inverse of forward (shift1, s_2+, subtracted, cmp, cswap, sub) in
        // reverse order: sub^-1, cswap^-1, cmp^-1, subtracted^-1, s_2+^-1, shift1^-1.
        // a) sub^-1: v += subtracted*u (X-sandwich cancels).
        let s_add = circ.push_section("rev_cadd");
        crate::point_add::trailmix_port::arith::gidney_const_adder::controlled_hybrid_add(
            circ,
            &subtracted,
            &v_full[..current_n],
            &u_full[..current_n],
            JUMP_GCD_VENTS,
        );
        circ.pop_section(&s_add);
        // b) cswap^-1 (involutory).
        let s_cs = circ.push_section("rev_cswap");
        for j in 0..current_n {
            circ.cswap(&swap_flag, &u_full[j], &v_full[j]);
        }
        circ.pop_section(&s_cs);
        // c) comparator^-1 (involutory) -> uncomputes swap_flag (v,u are in
        //    the post-shift pre-cswap state here, matching the forward cmp).
        msb_compare::controlled_lt_msbs_gidney(
            circ,
            &subtracted,
            &v_full[cmp_lo..current_n],
            &u_full[cmp_lo..current_n],
            cmp_eff,
            &swap_flag,
        );
        // d) subtracted^-1: cx(v[0], subtracted) -- v[0] is still the post-shift
        //    parity here (shifts not yet undone), so this clears subtracted.
        circ.cx(&v_full[0], &subtracted);
        // e) s_2+ inverse (reverse order): controlled left-shift, uncompute s_k.
        let s_sh = circ.push_section("rev_jshift");
        for s_k in s_flags.iter().rev() {
            controlled_left_shift(circ, s_k, &v_full[..current_n]);
            circ.x(s_k);
            circ.cx(&v_full[0], s_k);
        }
        // f) shift1 inverse: i>=1 unconditional left-shift; i==0 gated on t1
        //    (garbage[0]), then uncompute t1.
        if i == 0 {
            controlled_left_shift(circ, &garbage[0], &v_full[..current_n]);
            circ.x(&garbage[0]);
            circ.cx(&v_full[0], &garbage[0]);
        } else {
            poc_arith::left_shift(circ, &v_full[..current_n]);
        }
        circ.pop_section(&s_sh);
        // strict-dealloc touch (subtracted/swap/s_k are |0> here).
        for k in 0..current_n {
            circ.cx(&subtracted, &u_full[k]);
            circ.cx(&swap_flag, &u_full[k]);
        }
        for s_k in &s_flags {
            circ.cx(s_k, &u_full[0]);
        }
        // Drain the (now |0>) symbol slots; at i==0 also drain t1 (garbage[0]).
        for q in cur {
            circ.zero_and_free(q);
        }
        if i == 0 {
            let t1 = garbage.pop().expect("t1 prefix at garbage[0]");
            circ.zero_and_free(t1);
        }
        circ.pop_section(&iter_prev);
    }
    debug_assert!(garbage.is_empty(), "garbage tape not fully drained");

    // u_full grew back to sched[0]=256 bits holding q; deinit and free.
    let q_bytes = q_secp256k1().to_bytes_le();
    for i in 0..n.min(u_full.len()) {
        if (q_bytes.get(i / 8).copied().unwrap_or(0) >> (i % 8)) & 1 == 1 {
            circ.x(&u_full[i]);
        }
    }
    drop(subtracted);
    drop(swap_flag);
    drop(s_flags);
    drop(u_full);
    circ.pop_section(&prev);
}

/// Replay a RAW jump-before-swap dialog (reverse iter order) to reconstruct
/// the modular product. Mirrors `apply_bitvector_classical_jump_bs`: per step
/// `if subtracted: y += x; if swap: swap(x,y); then y := 2^j * y mod q`. The
/// x2^j is shift1 (i>=1 unconditional `mod_double`; i==0 gated on the t1 prefix
/// = garbage[0]) plus (jump-1) controlled doublings on `s_2..s_jump`.
///
/// On input (`x_reg=z`, `y_reg=0`) with a dialog from `forward_gcd_jump(q, x)`,
/// `y_reg` -> z*x mod q. `x_reg/y_reg` are n+1 = 257 bits (extra MSB = `mod_double`
/// overflow slot). `garbage` = [t1, `sym_0`, ..., sym_{iters-1}] with each sym =
/// (subtracted, swap, `s_2`, ..., `s_jump`).
pub fn apply_bitvector_jump_quantum_secp256k1(
    circ: &mut Circuit,
    garbage: &[QReg],
    x_reg: &[QReg],
    y_reg: &[QReg],
    jump: usize,
) {
    use crate::point_add::trailmix_port::arith::schrottenloher::pm_prims::{
        controlled_mod_add_pm_secp256k1_vents, mod_double_pm_secp256k1_vents,
    };
    let n = 256usize;
    assert!(jump >= 1);
    assert_eq!(x_reg.len(), n + 1, "x_reg must be n+1 = 257 bits");
    assert_eq!(y_reg.len(), n + 1, "y_reg must be n+1 = 257 bits");
    let iters = jump_iters_budget(jump);
    let sym_bits = jump + 1;
    assert!(garbage.len() > iters * sym_bits); // +1 = t1 prefix at garbage[0]
    let vents = 0usize; // apply-phase venting knob (0 for now; raised in JC)

    let prev = circ.push_section("apply_bv_jump");
    for i in (0..iters).rev() {
        let off = 1 + sym_bits * i; // garbage[0] = t1 prefix
                                    // 1) if subtracted: y += x mod q.
        controlled_mod_add_pm_secp256k1_vents(circ, &garbage[off], x_reg, y_reg, vents);
        // 2) if swap: swap(x, y).
        for j in 0..=n {
            circ.cswap(&garbage[off + 1], &x_reg[j], &y_reg[j]);
        }
        // 3) y := 2^j * y mod q: shift1 (i>=1 uncond; i==0 gated on t1) + s_2+.
        if i == 0 {
            controlled_mod_double_secp256k1(circ, &garbage[0], y_reg);
        } else {
            mod_double_pm_secp256k1_vents(circ, y_reg, vents);
        }
        for idx in 0..(jump - 1) {
            controlled_mod_double_secp256k1(circ, &garbage[off + 2 + idx], y_reg);
        }
    }
    circ.pop_section(&prev);
}

/// Compress a RAW jump=2 dialog `[t1, sym_0, ..., sym_{iters-1}]` (3 bits/sym)
/// in place into `[t1, code_0, ..., code_{nwin-1}]` (7 bits/window of 3 syms).
/// iters must be a multiple of 3 (the schedule rounds it). Frees 2 qubits per
/// window. Off-peak (called after `forward_gcd`, before the apply peak).
pub fn compress_dialog_jump(circ: &mut Circuit, garbage: &mut Vec<QReg>) {
    use crate::point_add::trailmix_port::arith::schrottenloher::gcd_compress_jump::compress_3sym_qrom_refs;
    let nsym = garbage.len() - 1;
    assert_eq!(
        nsym % 9,
        0,
        "iters must be a multiple of 3 for base-5 packing"
    );
    let nwin = nsym / 9;
    let prev = circ.push_section("b5pack_dialog");
    let dirty: Vec<QReg> = (0..13).map(|_| circ.alloc_qreg("b5d")).collect();
    let dref: Vec<&QReg> = dirty.iter().collect();
    for w in 0..nwin {
        let win: [&QReg; 9] = std::array::from_fn(|k| &garbage[1 + 9 * w + k]);
        compress_3sym_qrom_refs(circ, &win, &dref);
    }
    // compact: keep [t1] + 7/window, free the 2 (now |0>) high bits per window.
    let old = std::mem::take(garbage);
    let mut it = old.into_iter();
    garbage.push(it.next().expect("t1"));
    for _ in 0..nwin {
        for _ in 0..7 {
            garbage.push(it.next().expect("code bit"));
        }
        for _ in 0..2 {
            circ.zero_and_free(it.next().expect("freed bit"));
        }
    }
    debug_assert!(it.next().is_none());
    for q in dirty {
        circ.zero_and_free(q);
    }
    circ.pop_section(&prev);
}

/// Inverse of `compress_dialog_jump`: `[t1, code_*]` (7/window) -> raw
/// `[t1, sym_*]` (3/sym). Re-allocates 2 |0> bits per window, then decompresses.
pub fn decompress_dialog_jump(circ: &mut Circuit, garbage: &mut Vec<QReg>) {
    use crate::point_add::trailmix_port::arith::schrottenloher::gcd_compress_jump::compress_3sym_qrom_reverse_refs;
    let ncode = garbage.len() - 1;
    assert_eq!(ncode % 7, 0, "packed dialog must be 7 bits/window");
    let nwin = ncode / 7;
    let prev = circ.push_section("b5pack_dialog");
    // grow: [t1] + (7 code + 2 fresh |0>)/window.
    let old = std::mem::take(garbage);
    let mut it = old.into_iter();
    garbage.push(it.next().expect("t1"));
    for _ in 0..nwin {
        for _ in 0..7 {
            garbage.push(it.next().expect("code bit"));
        }
        garbage.push(circ.alloc_qreg("b5d_re"));
        garbage.push(circ.alloc_qreg("b5d_re"));
    }
    let dirty: Vec<QReg> = (0..13).map(|_| circ.alloc_qreg("b5d")).collect();
    let dref: Vec<&QReg> = dirty.iter().collect();
    for w in 0..nwin {
        let win: [&QReg; 9] = std::array::from_fn(|k| &garbage[1 + 9 * w + k]);
        compress_3sym_qrom_reverse_refs(circ, &win, &dref);
    }
    for q in dirty {
        circ.zero_and_free(q);
    }
    circ.pop_section(&prev);
}

/// Inverse apply (`y -> y * x^-1 mod q`): replays the dialog in FORWARD iter
/// order with each apply step inverted -- per step `y := 2^-j*y; if swap:
/// swap(x,y); if subtracted: y -= x`. Used by the in-place `mod_mul` reverse
/// (division). Same garbage layout [t1, `sym_0`, ...].
pub fn apply_bitvector_jump_quantum_secp256k1_inv(
    circ: &mut Circuit,
    garbage: &[QReg],
    x_reg: &[QReg],
    y_reg: &[QReg],
    jump: usize,
) {
    use crate::point_add::trailmix_port::arith::schrottenloher::pm_prims::{
        controlled_mod_sub_pm_secp256k1_vents, mod_halve_pm_secp256k1,
    };
    let n = 256usize;
    assert!(jump >= 1);
    assert_eq!(x_reg.len(), n + 1, "x_reg must be n+1 = 257 bits");
    assert_eq!(y_reg.len(), n + 1, "y_reg must be n+1 = 257 bits");
    let iters = jump_iters_budget(jump);
    let sym_bits = jump + 1;
    assert!(garbage.len() > iters * sym_bits);
    let vents = 0usize;

    let prev = circ.push_section("apply_bv_jump_inv");
    for i in 0..iters {
        let off = 1 + sym_bits * i;
        // 1) y := 2^-j * y mod q: shift1 halve (i>=1 uncond; i==0 gated on t1)
        //    + s_2+ controlled halves (inverse of the forward apply's x2^j).
        if i == 0 {
            controlled_mod_double_secp256k1_reverse(circ, &garbage[0], y_reg);
        } else {
            mod_halve_pm_secp256k1(circ, y_reg);
        }
        for idx in 0..(jump - 1) {
            controlled_mod_double_secp256k1_reverse(circ, &garbage[off + 2 + idx], y_reg);
        }
        // 2) if swap: swap(x, y).
        for j in 0..=n {
            circ.cswap(&garbage[off + 1], &x_reg[j], &y_reg[j]);
        }
        // 3) if subtracted: y -= x mod q.
        controlled_mod_sub_pm_secp256k1_vents(circ, &garbage[off], x_reg, y_reg, vents);
    }
    circ.pop_section(&prev);
}

/// Packed jump=2 apply (`y -> y*x mod q`): the dialog is base-5 PACKED
/// (`[t1, code_0(7), ...]`). Per 3-symbol window, decompress (re-alloc 2
/// ancilla, borrow `x_reg`[0..13] dirty) -> read the 3 symbols across the 3
/// reverse-order steps -> recompress -> free the ancilla. Same step as
/// `apply_bitvector_jump_quantum_secp256k1`.
/// `vents`/`coupled` tune the apply's q-q adds: `coupled=false` (low-qubit)
/// vents ONLY the register add (`f_vents=0`, +f stays borrowed-dirty, peak-safe);
/// `coupled=true` (low-tof) materializes + vents the +f window too, spending the
/// qubit headroom for a lower Toffoli count.
pub fn apply_bitvector_jump_packed_secp256k1(
    circ: &mut Circuit,
    garbage: &[QReg],
    x_reg: &[QReg],
    y_reg: &[QReg],
    vents: usize,
    coupled: bool,
) {
    use crate::point_add::trailmix_port::arith::schrottenloher::gcd_compress_jump::{
        compress_3sym_qrom_refs, compress_3sym_qrom_reverse_refs,
    };
    use crate::point_add::trailmix_port::arith::schrottenloher::pm_prims::{
        controlled_mod_add_pm_secp256k1_regvents, controlled_mod_add_pm_secp256k1_vents,
        mod_double_pm_secp256k1_vents,
    };
    let n = 256usize;
    let iters = jump_iters_budget(2);
    let nwin = iters / 3;
    assert_eq!(garbage.len(), 1 + 7 * nwin, "packed dialog length");
    let prev = circ.push_section("apply_bv_jump_pk");
    let mut cur_anc: Vec<QReg> = Vec::new();
    for i in (0..iters).rev() {
        let w = i / 3;
        let pos = i % 3;
        let cbase = 1 + 7 * w;
        if pos == 2 {
            cur_anc = vec![circ.alloc_qreg("jwin_a"), circ.alloc_qreg("jwin_a")];
            let win9: [&QReg; 9] = std::array::from_fn(|k| {
                if k < 7 {
                    &garbage[cbase + k]
                } else {
                    &cur_anc[k - 7]
                }
            });
            let dirty: Vec<&QReg> = (0..13).map(|k| &x_reg[k]).collect();
            compress_3sym_qrom_reverse_refs(circ, &win9, &dirty);
        }
        let resolve = |k: usize| -> &QReg {
            if k < 7 {
                &garbage[cbase + k]
            } else {
                &cur_anc[k - 7]
            }
        };
        let sub = resolve(3 * pos);
        let swp = resolve(3 * pos + 1);
        let s2 = resolve(3 * pos + 2);
        if coupled {
            controlled_mod_add_pm_secp256k1_vents(circ, sub, x_reg, y_reg, vents);
        } else {
            controlled_mod_add_pm_secp256k1_regvents(circ, sub, x_reg, y_reg, vents);
        }
        // bit n (overflow slot) is |0> on both x and y here (mod-add/sub clears
        // y[n]; mod_double's precondition keeps x[n]=0), so the n-th cswap is
        // 0<->0 -- skip it.
        for j in 0..n {
            circ.cswap(swp, &x_reg[j], &y_reg[j]);
        }
        if i == 0 {
            controlled_mod_double_secp256k1(circ, &garbage[0], y_reg);
        } else {
            mod_double_pm_secp256k1_vents(circ, y_reg, if coupled { vents } else { 0 });
        }
        controlled_mod_double_secp256k1(circ, s2, y_reg);
        if pos == 0 {
            let win9: [&QReg; 9] = std::array::from_fn(|k| {
                if k < 7 {
                    &garbage[cbase + k]
                } else {
                    &cur_anc[k - 7]
                }
            });
            let dirty: Vec<&QReg> = (0..13).map(|k| &x_reg[k]).collect();
            compress_3sym_qrom_refs(circ, &win9, &dirty);
            for q in std::mem::take(&mut cur_anc) {
                circ.zero_and_free(q);
            }
        }
    }
    circ.pop_section(&prev);
}

/// Packed jump=2 inverse apply (`y -> y*x^-1 mod q`): per-window decompress as
/// above, but forward iter order with the inverse step (2^-j -> swap -> -u).
pub fn apply_bitvector_jump_packed_inv_secp256k1(
    circ: &mut Circuit,
    garbage: &[QReg],
    x_reg: &[QReg],
    y_reg: &[QReg],
    vents: usize,
    coupled: bool,
) {
    use crate::point_add::trailmix_port::arith::schrottenloher::gcd_compress_jump::{
        compress_3sym_qrom_refs, compress_3sym_qrom_reverse_refs,
    };
    use crate::point_add::trailmix_port::arith::schrottenloher::pm_prims::{
        controlled_mod_sub_pm_secp256k1_regvents, controlled_mod_sub_pm_secp256k1_vents,
        mod_halve_pm_secp256k1,
    };
    let n = 256usize;
    let iters = jump_iters_budget(2);
    let nwin = iters / 3;
    assert_eq!(garbage.len(), 1 + 7 * nwin, "packed dialog length");
    let prev = circ.push_section("apply_bv_jump_pk_inv");
    let mut cur_anc: Vec<QReg> = Vec::new();
    for i in 0..iters {
        let w = i / 3;
        let pos = i % 3;
        let cbase = 1 + 7 * w;
        if pos == 0 {
            cur_anc = vec![circ.alloc_qreg("jwin_a"), circ.alloc_qreg("jwin_a")];
            let win9: [&QReg; 9] = std::array::from_fn(|k| {
                if k < 7 {
                    &garbage[cbase + k]
                } else {
                    &cur_anc[k - 7]
                }
            });
            let dirty: Vec<&QReg> = (0..13).map(|k| &x_reg[k]).collect();
            compress_3sym_qrom_reverse_refs(circ, &win9, &dirty);
        }
        let resolve = |k: usize| -> &QReg {
            if k < 7 {
                &garbage[cbase + k]
            } else {
                &cur_anc[k - 7]
            }
        };
        let sub = resolve(3 * pos);
        let swp = resolve(3 * pos + 1);
        let s2 = resolve(3 * pos + 2);
        if i == 0 {
            controlled_mod_double_secp256k1_reverse(circ, &garbage[0], y_reg);
        } else {
            mod_halve_pm_secp256k1(circ, y_reg);
        }
        controlled_mod_double_secp256k1_reverse(circ, s2, y_reg);
        // bit n (overflow slot) is |0> on both x and y here (mod-add/sub clears
        // y[n]; mod_double's precondition keeps x[n]=0), so the n-th cswap is
        // 0<->0 -- skip it.
        for j in 0..n {
            circ.cswap(swp, &x_reg[j], &y_reg[j]);
        }
        if coupled {
            controlled_mod_sub_pm_secp256k1_vents(circ, sub, x_reg, y_reg, vents);
        } else {
            controlled_mod_sub_pm_secp256k1_regvents(circ, sub, x_reg, y_reg, vents);
        }
        if pos == 2 {
            let win9: [&QReg; 9] = std::array::from_fn(|k| {
                if k < 7 {
                    &garbage[cbase + k]
                } else {
                    &cur_anc[k - 7]
                }
            });
            let dirty: Vec<&QReg> = (0..13).map(|k| &x_reg[k]).collect();
            compress_3sym_qrom_refs(circ, &win9, &dirty);
            for q in std::mem::take(&mut cur_anc) {
                circ.zero_and_free(q);
            }
        }
    }
    circ.pop_section(&prev);
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::{thread_rng, Rng};

    fn q_secp() -> BigUint {
        (BigUint::one() << 256u32) - BigUint::from(super::super::F_SECP256K1)
    }

    /// Does `x` fit the register-shrink schedule? (exact bitlen never exceeds
    /// current_n[i] at any step). Construction-time zero_and_free asserts the
    /// freed tail bits are |0>, so the test shots must all fit (the rare
    /// non-fitting tail is the generator's bounded miss rate). Mirrors the
    /// generator's `fits`.
    fn fits_schedule(x: &BigUint, jump: usize) -> bool {
        let (sched, iters) = super::super::jump_schedule::jump_schedule(jump);
        let q = q_secp();
        let one = BigUint::one();
        let (mut u, mut v) = (q, x.clone());
        for i in 0..iters {
            if u.bits().max(v.bits()) as usize > sched[i] as usize {
                return false;
            }
            // jump-before-swap step (matches the schedule generator).
            let mut j = 0usize;
            while j < jump && !v.is_zero() && !v.bit(0) {
                v >>= 1u32;
                j += 1;
            }
            if v.bit(0) {
                if u > v {
                    std::mem::swap(&mut u, &mut v);
                }
                v -= &u;
            }
            if u == one && v.is_zero() {
                return true;
            }
        }
        false
    }

    /// For random x in [1,q) and z in [0,q): forward jump-GCD on (q, x) must
    /// converge, and replaying its dialog on (z, 0) must give z*x mod q.
    /// Swept over jump in {1,2,3,4} (1 = the divstep baseline).
    #[test]
    fn jump_dialog_roundtrip_recovers_product() {
        let q = q_secp();
        let mut rng = thread_rng();
        for &jump in &[1usize, 2, 3, 4] {
            let mut converged = 0usize;
            let trials = 2000;
            let mut max_iters = 0usize;
            for _ in 0..trials {
                let x = (BigUint::from_bytes_le(&rng.gen::<[u8; 32]>()) % (&q - 1u32)) + 1u32;
                let z = BigUint::from_bytes_le(&rng.gen::<[u8; 32]>()) % &q;
                let (dialog, ok) = to_bitvector_classical_jump(q.clone(), x.clone(), jump, 1000);
                if !ok {
                    continue;
                }
                converged += 1;
                max_iters = max_iters.max(dialog.len());
                let (uo, vo) =
                    apply_bitvector_classical_jump(z.clone(), BigUint::zero(), &dialog, &q);
                let want = (&z * &x) % &q;
                assert!(
                    uo.is_zero(),
                    "jump={jump}: u should be 0 after apply, got {uo}"
                );
                assert_eq!(vo, want, "jump={jump}: apply must give z*x mod q");
            }
            eprintln!("jump={jump}: {converged}/{trials} converged, max_iters={max_iters}");
            assert_eq!(
                converged, trials,
                "jump={jump}: all should converge within cap"
            );
        }
    }

    /// Jump-BEFORE-swap roundtrip: forward GCD(q,x) -> dialog, apply on (z,0)
    /// recovers z*x mod q. Also checks the 5-symbol invariant (overflow symbol
    /// has swap=0; j in 1..=jump for non-step-0 symbols).
    #[test]
    fn jump_bs_dialog_roundtrip_recovers_product() {
        let q = q_secp();
        let mut rng = thread_rng();
        for &jump in &[1usize, 2, 3, 4] {
            let trials = 2000;
            let mut max_iters = 0usize;
            let mut converged = 0usize;
            let mut distinct: std::collections::HashSet<(u8, u8, u8)> = Default::default();
            for _ in 0..trials {
                let x = (BigUint::from_bytes_le(&rng.gen::<[u8; 32]>()) % (&q - 1u32)) + 1u32;
                let z = BigUint::from_bytes_le(&rng.gen::<[u8; 32]>()) % &q;
                let (dialog, ok) = to_bitvector_classical_jump_bs(q.clone(), x.clone(), jump, 1000);
                if !ok {
                    continue;
                }
                converged += 1;
                max_iters = max_iters.max(dialog.len());
                for s in &dialog {
                    distinct.insert((s.j, s.swap, s.subtracted));
                    if s.subtracted == 0 {
                        assert_eq!(s.swap, 0, "jump={jump}: overflow symbol must have swap=0");
                    }
                    assert!(s.j as usize <= jump, "jump={jump}: j out of range");
                }
                let (_uo, vo) =
                    apply_bitvector_classical_jump_bs(z.clone(), BigUint::zero(), &dialog, &q);
                assert_eq!(vo, (&z * &x) % &q, "jump={jump}: apply must give z*x mod q");
            }
            let mut syms: Vec<_> = distinct.iter().copied().collect();
            syms.sort();
            eprintln!(
                "jump-bs jump={jump}: {converged}/{trials} conv, max_iters={max_iters}, \
                 {} distinct symbols: {syms:?}",
                syms.len()
            );
            assert_eq!(
                converged, trials,
                "jump-bs jump={jump}: all should converge"
            );
        }
    }

    /// Base-5 packing arithmetic: sym<->digit roundtrip, pack3/unpack3 bijective
    /// into 7 bits, and every symbol the classical bs GCD (jump=2) emits maps to
    /// a valid 0..4 digit (the t1 prefix resolves the step-0 j=0 vs j=1 digit-1
    /// ambiguity at the apply).
    #[test]
    fn jump_bs_base5_packing_bijective() {
        for d in 0..5u8 {
            let (sub, sw, s2) = jump_bs_digit_to_sym(d);
            assert_eq!(jump_bs_sym_to_digit(sub, sw, s2), d);
        }
        let mut seen = std::collections::HashSet::new();
        for d0 in 0..5u8 {
            for d1 in 0..5u8 {
                for d2 in 0..5u8 {
                    let v = jump_bs_pack3([d0, d1, d2]);
                    assert!(v < 128, "packed {v} >= 128");
                    assert!(seen.insert(v), "pack collision at {v}");
                    assert_eq!(jump_bs_unpack3(v), [d0, d1, d2]);
                }
            }
        }
        assert_eq!(seen.len(), 125, "expect 125 distinct packed values");

        // Every classical bs symbol (jump=2) -> a valid 0..4 digit.
        let q = q_secp();
        let mut rng = thread_rng();
        for _ in 0..500 {
            let x = (BigUint::from_bytes_le(&rng.gen::<[u8; 32]>()) % (&q - 1u32)) + 1u32;
            let (dialog, ok) = to_bitvector_classical_jump_bs(q.clone(), x, 2, 1000);
            if !ok {
                continue;
            }
            for s in &dialog {
                // quantum view: subtracted=s.subtracted, swap=s.swap, s_2=(j==2).
                let s2 = u8::from(s.j == 2);
                let d = jump_bs_sym_to_digit(s.subtracted, s.swap, s2);
                assert!(
                    d < 5,
                    "sym (sub={},sw={},s2={s2}) -> {d}",
                    s.subtracted,
                    s.swap
                );
                assert_eq!(jump_bs_digit_to_sym(d), (s.subtracted, s.swap, s2));
            }
        }
    }

    /// MC-derived register-shrink schedule + iter budget per jump. Prints the
    /// iteration budget and `sum(current_n)` (the GCD csub/cswap Toffoli proxy)
    /// so we can see the jump win on BOTH steps and per-iter width. jump=1 is
    /// cross-checked against the closed-form sum (~59,974, gcd_pack.rs:257).
    #[test]
    fn mc_schedule_per_jump() {
        let q = q_secp();
        let samples = 40000;
        for &jump in &[1usize, 2, 3, 4] {
            let (sched, iters, total) = mc_schedule(&q, jump, samples, 6.0, 460);
            let head: Vec<usize> = sched.iter().take(4).copied().collect();
            let tail: Vec<usize> = sched.iter().rev().take(4).rev().copied().collect();
            eprintln!(
                "jump={jump}: iters_budget={iters}  sum(current_n)={total}  \
                 schedule head={head:?} ... tail={tail:?}"
            );
        }
        // closed-form jump-1 sum for reference (pad=37, 405 iters): ~59974.
    }

    /// Controlled mod-double: ctrl=1 shots give 2x mod q, ctrl=0 leave x;
    /// a[256] and ctrl clean afterward.
    #[test]
    fn controlled_mod_double_secp256k1_gated_value() {
        use crate::point_add::trailmix_port::compat::num_traits::Zero;
        let q = q_secp();
        let mut circ = Circuit::new();
        let a = circ.alloc_qreg_bits("a", 257);
        let ctrl = circ.alloc_qreg("ctrl");
        let mut rng = thread_rng();
        let mut expected: Vec<BigUint> = Vec::with_capacity(64);
        let mut ctrls: Vec<bool> = Vec::with_capacity(64);
        for shot in 0..64 {
            let x = BigUint::from_bytes_le(&rng.gen::<[u8; 32]>()) % &q;
            let mut buf = [0u8; 32];
            for (i, b) in x.to_bytes_le().iter().take(32).enumerate() {
                buf[i] = *b;
            }
            circ.sim_load_reg_bytes_shot(&a[..256], &buf, shot);
            let c = shot % 2 == 0;
            circ.sim_load_reg_bytes_shot(std::slice::from_ref(&ctrl), &[c as u8], shot);
            expected.push(if c { (&x * 2u32) % &q } else { x });
            ctrls.push(c);
        }

        super::controlled_mod_double_secp256k1(&mut circ, &ctrl, &a);

        let mut outs = a;
        outs.push(ctrl);
        let (sim, detached) = circ.destroy_sim(outs);
        for (shot, want) in expected.iter().enumerate() {
            let mut got = BigUint::zero();
            for i in 0..256 {
                if sim.read_bit_shot(&detached[i], shot) == 1 {
                    got |= BigUint::one() << i;
                }
            }
            assert_eq!(&got, want, "shot {shot} (ctrl={}): mismatch", ctrls[shot]);
            assert_eq!(
                sim.read_bit_shot(&detached[256], shot),
                0,
                "shot {shot}: a[256] dirty"
            );
            assert_eq!(
                sim.read_bit_shot(&detached[257], shot),
                ctrls[shot] as u8,
                "shot {shot}: ctrl perturbed"
            );
        }
    }

    /// Quantum forward jump-GCD then reverse must restore v_full=x and drain
    /// the garbage tape, phase-clean. Validates the variable-shift mechanism
    /// and its exact inverse in-circuit (jump=2, raw dialog, full width).
    #[test]
    fn jump_quantum_gcd_roundtrip_secp256k1() {
        use crate::point_add::trailmix_port::compat::num_traits::Zero;
        let jump = 2usize;
        let n = 256usize;
        let q = q_secp();
        let mut rng = thread_rng();

        // 64 random x in [1, q) that converge under the budget.
        let mut shot_data: Vec<BigUint> = Vec::new();
        while shot_data.len() < 64 {
            let mut x = BigUint::from_bytes_le(&rng.gen::<[u8; 32]>()) % &q;
            if x.is_zero() {
                x = BigUint::one();
            }
            let (_d, ok) =
                to_bitvector_classical_jump_bs(q.clone(), x.clone(), jump, jump_iters_budget(jump));
            if !ok || !fits_schedule(&x, jump) {
                continue;
            }
            shot_data.push(x);
        }

        let mut circ = Circuit::new();
        // Shrinking GCD + raw (jump+1)-bit dialog (~813q for jump=2). Peak is
        // dialog-tape-driven (~820q) since u+v shrink as the tape grows.
        circ.set_max_qubit_peak(950);
        let mut v_full = circ.alloc_qreg_bits("v_full", n);
        for (shot, x) in shot_data.iter().enumerate() {
            let mut buf = [0u8; 32];
            for (i, b) in x.to_bytes_le().iter().take(32).enumerate() {
                buf[i] = *b;
            }
            circ.sim_load_reg_bytes_shot(&v_full[..n], &buf, shot);
        }

        let ccx0 = circ.ccx_emitted;
        let ccz0 = circ.ccz_emitted;
        let mut garbage = forward_gcd_jump_quantum_secp256k1(&mut circ, &mut v_full, jump);
        forward_gcd_jump_quantum_secp256k1_reverse(&mut circ, &mut v_full, &mut garbage, jump);
        assert!(garbage.is_empty(), "roundtrip should drain garbage tape");
        let rt_tof = (circ.ccx_emitted - ccx0) + (circ.ccz_emitted - ccz0);
        eprintln!(
            "  jump={jump} fwd+rev GCD: tof={rt_tof} peak={} iters={}",
            circ.peak_qubits,
            jump_iters_budget(jump)
        );

        circ.assert_phase_clean();
        let (sim, detached) = circ.destroy_sim(v_full);
        for (shot, x) in shot_data.iter().enumerate() {
            let mut got = BigUint::zero();
            for i in 0..n {
                if sim.read_bit_shot(&detached[i], shot) == 1 {
                    got |= BigUint::one() << i;
                }
            }
            assert_eq!(got, *x, "shot {shot}: v not restored");
        }
    }

    /// Compress/decompress dialog passes are inverse: forward -> compress (tape
    /// shrinks) -> decompress (tape grows back) -> reverse restores v.
    #[test]
    fn jump_packed_dialog_roundtrip_secp256k1() {
        let jump = 2usize;
        let n = 256usize;
        let q = q_secp();
        let mut rng = thread_rng();
        let mut shot_data: Vec<BigUint> = Vec::new();
        while shot_data.len() < 64 {
            let mut x = BigUint::from_bytes_le(&rng.gen::<[u8; 32]>()) % &q;
            if x.is_zero() {
                x = BigUint::one();
            }
            let (_d, ok) =
                to_bitvector_classical_jump_bs(q.clone(), x.clone(), jump, jump_iters_budget(jump));
            if !ok || !fits_schedule(&x, jump) {
                continue;
            }
            shot_data.push(x);
        }
        let mut circ = Circuit::new();
        circ.set_max_qubit_peak(950);
        let mut v_full = circ.alloc_qreg_bits("v_full", n);
        for (shot, x) in shot_data.iter().enumerate() {
            let mut buf = [0u8; 32];
            for (i, b) in x.to_bytes_le().iter().take(32).enumerate() {
                buf[i] = *b;
            }
            circ.sim_load_reg_bytes_shot(&v_full[..n], &buf, shot);
        }
        let mut garbage = forward_gcd_jump_quantum_secp256k1(&mut circ, &mut v_full, jump);
        let raw_len = garbage.len();
        compress_dialog_jump(&mut circ, &mut garbage);
        let packed_len = garbage.len();
        assert!(
            packed_len < raw_len,
            "compress should shrink ({raw_len}->{packed_len})"
        );
        decompress_dialog_jump(&mut circ, &mut garbage);
        assert_eq!(garbage.len(), raw_len, "decompress should restore length");
        eprintln!(
            "  jump=2 dialog: raw={raw_len} packed={packed_len} (save {})",
            raw_len - packed_len
        );
        forward_gcd_jump_quantum_secp256k1_reverse(&mut circ, &mut v_full, &mut garbage, jump);
        assert!(garbage.is_empty(), "tape drained");
        circ.assert_phase_clean();
        let (sim, detached) = circ.destroy_sim(v_full);
        for (shot, x) in shot_data.iter().enumerate() {
            let mut got = BigUint::zero();
            for i in 0..n {
                if sim.read_bit_shot(&detached[i], shot) == 1 {
                    got |= BigUint::one() << i;
                }
            }
            assert_eq!(got, *x, "shot {shot}: v not restored");
        }
    }

    /// In-circuit multiply via the jump dialog: forward_gcd_jump(x) -> dialog,
    /// apply on (x_reg=z, y_reg=0) -> y_reg = z*x mod q (and x_reg = 0), then
    /// reverse_gcd_jump restores v=x and drains the tape, phase-clean.
    #[test]
    fn jump_apply_multiplies_secp256k1() {
        use crate::point_add::trailmix_port::compat::num_traits::Zero;
        let jump = 2usize;
        let n = 256usize;
        let q = q_secp();
        let mut rng = thread_rng();

        let mut xs: Vec<BigUint> = Vec::new();
        let mut zs: Vec<BigUint> = Vec::new();
        while xs.len() < 64 {
            let mut x = BigUint::from_bytes_le(&rng.gen::<[u8; 32]>()) % &q;
            if x.is_zero() {
                x = BigUint::one();
            }
            let (_d, ok) =
                to_bitvector_classical_jump_bs(q.clone(), x.clone(), jump, jump_iters_budget(jump));
            if !ok || !fits_schedule(&x, jump) {
                continue;
            }
            let z = BigUint::from_bytes_le(&rng.gen::<[u8; 32]>()) % &q;
            xs.push(x);
            zs.push(z);
        }

        let mut circ = Circuit::new();
        // Shrinking GCD; peak is the apply_bv phase: raw (jump+1)-bit tape
        // (~813q) + 2 mul regs (514q) + working. Packing (JC) brings the tape
        // (and thus the peak) down toward ~1175.
        circ.set_max_qubit_peak(1450);
        let mut v_reg = circ.alloc_qreg_bits("v_reg", n); // GCD input (x)
        let x_reg = circ.alloc_qreg_bits("x_reg", n + 1); // multiplicand z
        let y_reg = circ.alloc_qreg_bits("y_reg", n + 1); // accumulator -> z*x

        let load = |circ: &mut Circuit, reg: &[QReg], val: &BigUint, shot: usize| {
            let mut buf = [0u8; 32];
            for (i, b) in val.to_bytes_le().iter().take(32).enumerate() {
                buf[i] = *b;
            }
            circ.sim_load_reg_bytes_shot(&reg[..n], &buf, shot);
        };
        for (shot, (x, z)) in xs.iter().zip(zs.iter()).enumerate() {
            load(&mut circ, &v_reg, x, shot);
            load(&mut circ, &x_reg, z, shot);
        }

        let mut garbage = forward_gcd_jump_quantum_secp256k1(&mut circ, &mut v_reg, jump);
        apply_bitvector_jump_quantum_secp256k1(&mut circ, &garbage, &x_reg, &y_reg, jump);
        forward_gcd_jump_quantum_secp256k1_reverse(&mut circ, &mut v_reg, &mut garbage, jump);
        assert!(garbage.is_empty(), "tape drained");
        circ.assert_phase_clean();

        let mut outs: Vec<QReg> = Vec::new();
        outs.extend(v_reg);
        outs.extend(x_reg);
        outs.extend(y_reg);
        let (sim, detached) = circ.destroy_sim(outs);
        let read = |off: usize, bits: usize, shot: usize| -> BigUint {
            let mut v = BigUint::zero();
            for i in 0..bits {
                if sim.read_bit_shot(&detached[off + i], shot) == 1 {
                    v |= BigUint::one() << i;
                }
            }
            v
        };
        for (shot, (x, z)) in xs.iter().zip(zs.iter()).enumerate() {
            let v_got = read(0, n, shot);
            let x_got = read(n, n + 1, shot);
            let y_got = read(2 * n + 1, n + 1, shot);
            let want = (z * x) % &q;
            // Alg-7 reductions keep values in [0, 2^256) (residue-class, not
            // strictly [0,q)), so check mod q. x_reg is the cofactor (= a
            // multiple of q, ≡0 mod q), cleaned by the reverse division in the
            // full in-place mod_mul.
            assert_eq!(&y_got % &q, want, "shot {shot}: y != z*x mod q (y={y_got})");
            assert!(
                (&x_got % &q).is_zero(),
                "shot {shot}: x_reg not ≡0 mod q (x={x_got})"
            );
            assert_eq!(v_got, *x, "shot {shot}: v not restored");
        }
    }
}
