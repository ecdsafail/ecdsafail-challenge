//! The product-min jump-GCD modular inversion (`y -> y * x^-1 mod q`
//! and the forward `y -> y * x` direction), built on this crate's `B`
//! builder and the `schedule` / `arith` / `comparator` / `codec` modules.
//!
//! ## Fused apply with incremental codec
//! Schedule `schedule::SCHED_J2` / `GAP_J2`, jump=2, `ITERS = 258`. The dialog
//! apply is fused into the forward/reverse GCD: each divstep symbol is applied to
//! the coordinate pair the instant it is computed, so the apply adders run in the
//! GCD passes' headroom and the full raw tape is never materialized.
//!
//!   1. `forward_gcd_jump` runs the jump=2 divstep loop, computing the per-step
//!      dialog `(subtracted, swap, s_2)`. When `apply_inv` is set it applies the
//!      inverse step to `[y, tmp]` (`Direction::Inverse`, producing `y*x^-1`)
//!      before the symbol is swapped to the tape; otherwise it only records.
//!   2. each codec window is compressed inline the instant its symbols are
//!      recorded, so the resident tape stays compressed.
//!   3. `reverse_gcd_jump` is the exact gate-inverse of the divstep loop: it
//!      restores `x` and drains the tape to empty, decompressing one window at a
//!      time. When `apply_fwd` is set it applies the forward step to `[tmp, y]`
//!      (`Direction::Forward`, producing `y*x`) in reverse iter order.
//!
//! The compressed tape held live across a pass is `dialog_tape_qubits` qubits
//! (603 for the len-258 all-triple tiling -- see `codec::dialog_tape_qubits`).
//!
//! ## Algorithm structure
//! - The classical jump-before-swap divstep is the algebra `forward_gcd_jump`
//!   realizes in the circuit.
//! - The forward / reverse GCD: the per-step
//!   shift / swap-decision / cswap / subtract structure and the schedule-driven
//!   register shrink/regrow, with the comparator/subtract routed through this
//!   crate's `comparator` / `arith`.
//! - The reverse-iter apply (`if subtracted: y+=x; if swap: swap(x,y);
//!   y *= 2^j mod q`): the multiply is the apply, the inversion is the
//!   reverse-of-multiply composition.
//! - The per-step variable shift (`controlled_right_shift` /
//!   `controlled_left_shift`) and the apply's controlled doubling.

use super::arith::{self, F_SECP256K1};
use super::schedule::{GAP_J2, ITERS, JUMP, SCHED_J2};
use super::{B, BExt};
use crate::circuit::{QubitId};

/// Which product the apply pass reconstructs into `y`.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Direction {
    /// `y -> y * x^-1 mod q` (the modular inversion: apply the dialog reversed
    /// vs the multiply, i.e. forward iter order with the inverse step).
    Inverse,
    /// `y -> y * x   mod q` (the multiply: reverse iter order, forward step).
    Forward,
}

/// secp256k1 prime q = 2^256 - f = 0xFFFF...FFFFFC2F, little-endian bytes.
/// (Build-time constant; avoids a bignum dependency in the submission.)
#[must_use]
pub fn q_secp256k1_le() -> [u8; 32] {
    let mut b = [0xFFu8; 32];
    b[0] = 0x2F;
    b[1] = 0xFC;
    b[4] = 0xFE;
    b
}

/// All-triple group count the product-min selector picks for this
/// schedule: `n3 = iters/3` (clamped inside `codec::jump_dialog_regions`).
#[must_use]
pub fn n3_for_iters(iters: usize) -> usize {
    iters / 3
}

// =====================================================================
// MBU AND-uncompute (HMR + conditional-Z), measurement-vented. `t` holds
// `a AND b` and is returned to |0>.
// =====================================================================

fn clear_and(circ: &mut B, t: &QubitId, a: &QubitId, b: &QubitId) {
    let bit = circ.alloc_bit();
    circ.hmr(*t, bit);
    circ.cz_if_bit(*a, *b, bit);
}

// =====================================================================
// Per-step variable shift (controlled right/left shift by 1).
// =====================================================================

/// Controlled logical right-shift-by-1 (LSB-first): when `ctrl`, v[i+1] -> v[i].
fn controlled_right_shift(circ: &mut B, ctrl: &QubitId, v: &[QubitId]) {
    for i in 0..v.len().saturating_sub(1) {
        circ.cswap(*ctrl, v[i], v[i + 1]);
    }
}

/// Exact gate-inverse of [`controlled_right_shift`] (cswaps reversed).
fn controlled_left_shift(circ: &mut B, ctrl: &QubitId, v: &[QubitId]) {
    for i in (1..v.len()).rev() {
        circ.cswap(*ctrl, v[i], v[i - 1]);
    }
}

/// Unconditional logical right-shift-by-1: v[i+1] -> v[i], top bit -> |0>.
/// A swap chain from the bottom up; v[0] is |0> in the GCD post-subtract so it
/// wraps harmlessly to the top.
fn right_shift(circ: &mut B, v: &[QubitId]) {
    for i in 0..v.len().saturating_sub(1) {
        circ.swap(v[i], v[i + 1]);
    }
}

/// Exact gate-inverse of [`right_shift`].
fn left_shift(circ: &mut B, v: &[QubitId]) {
    for i in (1..v.len()).rev() {
        circ.swap(v[i], v[i - 1]);
    }
}

// =====================================================================
// Controlled mod_double for the apply (gated `y := 2*y mod q`), built on
// this crate's arith::add_f_window (re-using add_f via a ctrl).
// =====================================================================

/// If `ctrl`: `a := 2*a mod q`; else `a` unchanged. `a.len() == n = 256`; the n+1
/// shift view is `a ++ ovf` with a transient `ovf` (not a slot in `a`), restored to
/// |0>. Normal controlled shift + the (MBU) `+f` fold gated on the shifted-out
/// overflow + an MBU ancilla clean.
fn controlled_mod_double(circ: &mut B, ctrl: &QubitId, a: &[QubitId]) {
    let n = a.len();
    assert_eq!(n, 256, "controlled_mod_double expects 256-bit a");
    let f_bytes = F_SECP256K1.to_le_bytes();
    let ovf = circ.alloc_qubit();
    // 1) controlled left-shift by 1 over the n+1 view a ++ ovf: MSB -> ovf.
    let w: Vec<&QubitId> = a.iter().chain(std::iter::once(&ovf)).collect();
    for i in (0..n).rev() {
        circ.cswap(*ctrl, *w[i], *w[i + 1]);
    }
    // 2) if ovf (= ctrl AND old-MSB), add f to the low LSBS bits.
    arith::add_f_window_pub(circ, &ovf, a, arith::LSBS, &f_bytes, None);
    // 3) clear ovf (= ctrl AND a[0], post-fold) by MBU.
    clear_and(circ, &ovf, ctrl, &a[0]);
    circ.zero_and_free(ovf);
}

/// Exact gate-inverse of [`controlled_mod_double`]: if `ctrl`, `a := a/2 mod q`.
fn controlled_mod_double_reverse(circ: &mut B, ctrl: &QubitId, a: &[QubitId]) {
    let n = a.len();
    assert_eq!(n, 256, "controlled_mod_double_reverse expects 256-bit a");
    let f_bytes = F_SECP256K1.to_le_bytes();
    let ovf = circ.alloc_qubit();
    // inverse of 3): rebuild ovf = ctrl AND a[0] with a CCX.
    circ.ccx(*ctrl, a[0], ovf);
    // inverse of 2): subtract f (X-sandwich of add_f_window) gated on ovf.
    for q in &a[..arith::LSBS] {
        circ.x(*q);
    }
    arith::add_f_window_pub(circ, &ovf, a, arith::LSBS, &f_bytes, None);
    for q in &a[..arith::LSBS] {
        circ.x(*q);
    }
    // inverse of 1): controlled right-shift over a ++ ovf; ovf shifted back to |0>.
    let w: Vec<&QubitId> = a.iter().chain(std::iter::once(&ovf)).collect();
    for i in 0..n {
        circ.cswap(*ctrl, *w[i], *w[i + 1]);
    }
    circ.zero_and_free(ovf);
}

// =====================================================================
// Forward jump-GCD: record the dialog onto the compressed tape.
// =====================================================================

/// Forward jump=2 GCD on `(u = q, v = x)` over the baked schedule. `v` is the
/// input register (>= 256 bits; bits beyond the active width must be |0> for
/// schedule-fitting inputs). Computes the per-step dialog `(subtracted, swap,
/// s_2)`, compresses each codec window inline, and returns the compressed tape
/// (`[t1, win_0 code bits, ..]`, `codec::dialog_tape_qubits` long). On a fitting
/// input, ends with `u = 0` (X'd from 1), `v = 0`, both shrunk away; the tape
/// holds the whole dialog.
///
/// When `apply_inv` is `Some((x_reg, y_reg))`, the inverse apply step is fused in
/// (applied to the coordinate pair each iter before the symbol is taped);
/// otherwise the coordinate registers are untouched.
#[must_use]
pub fn forward_gcd_jump(circ: &mut B, v: &mut Vec<QubitId>, apply_inv: Option<(&[QubitId], &[QubitId])>) -> Vec<QubitId> {
    let n = 256usize;
    assert_eq!(JUMP, 2, "ludicrous apply/codec are jump=2 specific");
    assert!(v.len() >= n, "v must be at least n=256 bits");
    let iters = ITERS;
    let sym_bits = 3; // (subtracted, swap, s_2)

    // u_full = q (n bits); shrinks per the schedule.
    let mut u: Vec<QubitId> = (0..n).map(|_| circ.alloc_qubit()).collect();
    let q_bytes = q_secp256k1_le();
    for (i, qb) in u.iter().enumerate() {
        if (q_bytes.get(i / 8).copied().unwrap_or(0) >> (i % 8)) & 1 == 1 {
            circ.x(*qb);
        }
    }

    let subtracted = circ.alloc_qubit();
    let swap_flag = circ.alloc_qubit();
    let s2 = circ.alloc_qubit(); // the single jump=2 extra-shift flag
    let t1 = circ.alloc_qubit(); // step-0 shift1-fired flag (= x even)

    // Incremental all-triple codec: compress each codec window inline the instant
    // its symbols are recorded, so the resident tape grows compressed (~603) not
    // raw (775), bounding the peak. `tape` accumulates the t1 prefix
    // followed by each window's compressed code bits, in iter order.
    let n3 = n3_for_iters(iters);
    let mut window_plan: Vec<super::codec::DialogCodec> = Vec::new();
    for (codec, count) in super::codec::jump_dialog_regions(n3, iters) {
        for _ in 0..count {
            window_plan.push(codec);
        }
    }
    let mut tape: Vec<QubitId> = Vec::with_capacity(super::codec::dialog_tape_qubits(n3, iters));
    let mut win_idx = 0usize; // current window
    let mut pending: Vec<QubitId> = Vec::new(); // raw symbol slots of the current window
    for i in 0..iters {
        let current_n = (SCHED_J2[i] as usize).max(1);
        while u.len() > current_n {
            let q = u.pop().expect("u nonempty");
            circ.zero_and_free(q);
        }
        while v.len() > current_n {
            let q = v.pop().expect("v nonempty");
            circ.zero_and_free(q);
        }
        // swap-decision window: top GAP_J2[i] bits of the active width.
        let cmp_eff = (GAP_J2[i] as usize).min(current_n).max(1);

        // 1) Shift-first: remove up to jump=2 trailing zeros of v.
        //    i==0 gates shift1 on (v even) and records it in t1; i>=1 is
        //    unconditional (v even post-subtract).
        if i == 0 {
            circ.cx(v[0], t1); // t1 = v[0]
            circ.x(t1); // t1 = NOT(v[0]) = (x even)
            controlled_right_shift(circ, &t1, &v[..current_n]);
        } else {
            right_shift(circ, &v[..current_n]);
        }
        // s_2: shift again while still even. Nesting is automatic (once an odd
        // bit reaches LSB the shift stops).
        circ.cx(v[0], s2);
        circ.x(s2);
        controlled_right_shift(circ, &s2, &v[..current_n]);

        // 2) subtracted = v[0] (post-shift parity): 1 => odd => swap+subtract.
        circ.cx(v[0], subtracted);

        // 3) swap decision. Step 0: u=q, v=x<q always, so (v<u)=1 deterministically,
        //    so swap_flag = subtracted (one CX). Else the narrow top-k comparator
        //    decides swap_flag ^= subtracted AND (v < u).
        if i == 0 {
            circ.cx(subtracted, swap_flag);
        } else {
            controlled_swap_decision_v_lt_u(
                circ,
                &subtracted,
                &v[..current_n],
                &u[..current_n],
                cmp_eff,
                &swap_flag,
            );
        }
        // 4) cswap(swap_flag, u, v).
        for j in 0..current_n {
            circ.cswap(swap_flag, u[j], v[j]);
        }
        // 5) v -= subtracted * u (controlled mod-free subtract on the active width;
        //    post-swap v >= u so no borrow on a fitting input). X-sandwich add.
        for q in &v[..current_n] {
            circ.x(*q);
        }
        controlled_add_active(circ, &subtracted, &u[..current_n], &v[..current_n]);
        for q in &v[..current_n] {
            circ.x(*q);
        }

        // Fused inverse apply (Direction::Inverse): apply this divstep's symbol to the
        // coordinate pair using the live symbol bits, before they are swapped to the
        // tape. Forward iter order matches apply_step_reverse's order, so the tape is
        // never materialized for the apply -> the apply adders run in the GCD's headroom.
        if let Some((xr, yr)) = apply_inv {
            apply_step_reverse(circ, i, &subtracted, &swap_flag, &s2, &t1, xr, yr);
        }

        // 6) record the symbol into fresh |0> slots (returning the ancilla to |0>).
        let slots: Vec<QubitId> = (0..sym_bits).map(|_| circ.alloc_qubit()).collect();
        circ.swap(subtracted, slots[0]);
        circ.swap(swap_flag, slots[1]);
        circ.swap(s2, slots[2]);
        pending.extend(slots);

        // When the current window's symbols are complete, compress it inline.
        let codec = window_plan[win_idx];
        if pending.len() == codec.syms() * sym_bits {
            let data = codec.compress_window(circ, &pending);
            tape.extend(data);
            pending.clear();
            win_idx += 1;
        }
    }
    assert_eq!(win_idx, window_plan.len(), "all windows compressed");
    assert!(pending.is_empty(), "no leftover symbols");

    // u converged to 1; X to 0 then free. v=0; shrink away.
    circ.x(u[0]);
    while let Some(q) = v.pop() {
        circ.zero_and_free(q);
    }
    for q in u {
        circ.zero_and_free(q);
    }
    circ.zero_and_free(subtracted);
    circ.zero_and_free(swap_flag);
    circ.zero_and_free(s2);
    // Prepend t1: tape = [t1, win_0 code bits, ..]  (compressed layout).
    tape.insert(0, t1);
    assert_eq!(tape.len(), super::codec::dialog_tape_qubits(n3, iters));
    tape
}

/// Reverse of [`forward_gcd_jump`]: restores `v[..256]` to the original `x` and
/// drains the compressed `tape` to empty (all |0>, freed). Exact gate-inverse,
/// step by step. Decompresses one codec window at a time (in reverse iter order),
/// consumes its symbols, and frees the raw slots -- so the resident tape stays
/// compressed.
pub fn reverse_gcd_jump(circ: &mut B, v: &mut Vec<QubitId>, tape: &mut Vec<QubitId>, apply_fwd: Option<(&[QubitId], &[QubitId])>) {
    let n = 256usize;
    let iters = ITERS;
    let n3 = n3_for_iters(iters);
    assert_eq!(
        tape.len(),
        super::codec::dialog_tape_qubits(n3, iters),
        "tape must be the compressed dialog"
    );

    // Window plan (iter order); we consume windows from the last one back.
    let mut window_plan: Vec<super::codec::DialogCodec> = Vec::new();
    for (codec, count) in super::codec::jump_dialog_regions(n3, iters) {
        for _ in 0..count {
            window_plan.push(codec);
        }
    }
    let mut win_idx = window_plan.len(); // next window to decompress (from the end)
    // `pending` holds the raw symbol slots of the currently-decompressed window,
    // consumed symbol-by-symbol from the end (reverse symbol order).
    let mut pending: Vec<QubitId> = Vec::new();

    // u regrows from 1 bit (forward ended u=0 post-X; re-init u_final=1 via X).
    let mut u: Vec<QubitId> = vec![circ.alloc_qubit()];
    circ.x(u[0]);

    let subtracted = circ.alloc_qubit();
    let swap_flag = circ.alloc_qubit();
    let s2 = circ.alloc_qubit();

    for i in (0..iters).rev() {
        let current_n = (SCHED_J2[i] as usize).max(1);
        while u.len() < current_n {
            u.push(circ.alloc_qubit());
        }
        while v.len() < current_n {
            v.push(circ.alloc_qubit());
        }
        let cmp_eff = (GAP_J2[i] as usize).min(current_n).max(1);

        // If the current window is exhausted, decompress the next one (from the
        // tape end) into raw symbol slots.
        if pending.is_empty() {
            win_idx -= 1;
            let codec = window_plan[win_idx];
            let cb = codec.code_bits();
            let tlen = tape.len();
            let data: Vec<QubitId> = tape.split_off(tlen - cb);
            pending = codec.decompress_window(circ, &data);
        }
        // Pull the last symbol (3 bits) off `pending` into the ancilla.
        let plen = pending.len();
        let cur: Vec<QubitId> = pending.split_off(plen - 3);
        circ.swap(subtracted, cur[0]);
        circ.swap(swap_flag, cur[1]);
        circ.swap(s2, cur[2]);
        // Free the 3 now-|0> symbol slots here -- before the apply -- so they are not
        // carried as dead live qubits through the apply + GCD-undo, where they would
        // inflate the peak by sym_bits=3. The next window's decompress re-allocs
        // from these freed slots.
        for q in cur {
            circ.zero_and_free(q);
        }

        // Fused forward apply (Direction::Forward multiply): apply this divstep's
        // symbol to the coordinate pair using the live symbol bits, in reverse iter
        // order matching apply_step_forward's order, so the tape is never materialized
        // for the apply -> the apply adders run in the (reverse) GCD's headroom.
        if let Some((xr, yr)) = apply_fwd {
            let t1 = tape[0];
            apply_step_forward(circ, i, &subtracted, &swap_flag, &s2, &t1, xr, yr);
        }

        // Inverse step (reverse op order): sub^-1, cswap^-1, cmp^-1, subtracted^-1,
        // s_2^-1, shift1^-1.
        // a) sub^-1: v += subtracted*u (X-sandwich cancels).
        controlled_add_active(circ, &subtracted, &u[..current_n], &v[..current_n]);
        // b) cswap^-1 (involutory).
        for j in 0..current_n {
            circ.cswap(swap_flag, u[j], v[j]);
        }
        // c) uncompute swap_flag. Step 0 is a CNOT (swap_flag == subtracted). For
        //    i>=1 the flag holds `subtracted AND (v<u)`; clear it by measurement-vent
        //    (hmr + capped Z-recompute under push_condition) -- ~half the normal
        //    comparator's Toffoli, scorer-discounted ~0.5 more.
        if i == 0 {
            circ.cx(subtracted, swap_flag);
        } else {
            super::comparator::swap_decision_uncompute_vented(
                circ,
                &subtracted,
                &v[..current_n],
                &u[..current_n],
                cmp_eff,
                &swap_flag,
            );
        }
        // d) subtracted^-1: cx(v[0], subtracted) (v[0] still post-shift parity).
        circ.cx(v[0], subtracted);
        // e) s_2 inverse: controlled left-shift, uncompute s_2.
        controlled_left_shift(circ, &s2, &v[..current_n]);
        circ.x(s2);
        circ.cx(v[0], s2);
        // f) shift1 inverse: i>=1 unconditional left-shift; i==0 gated on t1.
        if i == 0 {
            controlled_left_shift(circ, &tape[0], &v[..current_n]);
            circ.x(tape[0]);
            circ.cx(v[0], tape[0]);
        } else {
            left_shift(circ, &v[..current_n]);
        }

        // (symbol slots already freed before the apply, above). At i==0 drain t1.
        if i == 0 {
            let t1 = tape.pop().expect("t1 prefix");
            circ.zero_and_free(t1);
        }
    }
    assert!(tape.is_empty(), "tape not fully drained");

    // u grew back to 256 bits holding q; deinit and free.
    let q_bytes = q_secp256k1_le();
    for (i, qb) in u.iter().enumerate().take(n) {
        if (q_bytes.get(i / 8).copied().unwrap_or(0) >> (i % 8)) & 1 == 1 {
            circ.x(*qb);
        }
    }
    for q in u {
        circ.zero_and_free(q);
    }
    circ.zero_and_free(subtracted);
    circ.zero_and_free(swap_flag);
    circ.zero_and_free(s2);
}

// =====================================================================
// Active-width helpers (the GCD body runs on `current_n` bits, NOT n=256).
// =====================================================================

/// `target ^= ctrl AND (v_top < u_top)` on the top `k` MSBs of the active-width
/// slices.
fn controlled_swap_decision_v_lt_u(
    circ: &mut B,
    ctrl: &QubitId,
    v: &[QubitId],
    u: &[QubitId],
    k: usize,
    target: &QubitId,
) {
    super::comparator::controlled_swap_decision_lt_truncated(circ, ctrl, v, u, k, target);
}

/// `y += ctrl * x (mod 2^width)` over the active width (no carry-out captured).
/// The GCD subtract uses this inside an X-sandwich; post-swap `v >= u` so the
/// two's-complement subtract never borrows out on a fitting input.
///
/// Measurement-vented threaded add: ~2 Toffoli/bit (one forward carry +
/// one gated sum), with the carry-uncompute vented by `hmr` (0 Toffoli) and no
/// gated-addend materialization -- vs the single-ancilla Cuccaro's ~3 Toffoli/bit.
/// The per-bit carry qubits (allocated + freed inside) fit in the GCD passes'
/// headroom below the global (apply) peak. mod-2^m: the carry-out is dropped, as
/// the GCD subtract never carries out on a fitting input.
fn controlled_add_active(circ: &mut B, ctrl: &QubitId, x: &[QubitId], y: &[QubitId]) {
    // The GCD subtract `v -= u` (here `y += x` inside the X-sandwich):
    // a = target = y (= v), b = addend = x (= u). cap = PAD.
    // Schedule-driven GCD subtract: pull the baked carry-cap `k` for this call
    // (next_gcd_k) and emit the capped varchunk/adaptive adder. The baked `k`
    // fixes the chunk decomposition deterministically.
    // Kaliski bit-0 optimization: under the odd-u invariant the subtrahend bit0
    // `x[0] == 1`, and the accumulator bit0 is the control (forward, inside the
    // X-sandwich `~v[0] == NOT ctrl`) or 0 (reverse-add). Either way the bit-0
    // carry-out `acc[0] AND ctrl AND x[0]` is provably 0, so the bit-0 result is the
    // known `cx(ctrl, y[0])` with no carry into bit 1 -- emit it directly and run the
    // capped adder on bits 1.. with carry-in 0. Saves the bit-0 carry CCX (~1-2 tof)
    // per GCD conditional sub/add * 2 (fwd+rev) * ITERS ~= 1000+ tof. Not the apply.
    let k = super::next_gcd_k();
    let branch = super::next_gcd_branch();
    circ.cx(*ctrl, y[0]); // bit-0 sum (x[0] == 1 under the odd-u invariant)
    if x.len() > 1 {
        let yr: Vec<&QubitId> = y[1..].iter().collect();
        let xr: Vec<&QubitId> = x[1..].iter().collect();
        super::gidney::controlled_hybrid_add_capped_branch(circ, ctrl, &yr, &xr, k, super::PAD, branch);
    }
}

// =====================================================================
// Apply: read the (decompressed) tape and reconstruct the modular product.
// =====================================================================

/// One forward apply step (the multiply body) at iter `i` with the symbol bits
/// `(sub, swp, s2)` and the t1 prefix `t1`:
///   if sub: y += x mod q; if swp: swap(x,y); y *= 2 (shift1) ; if s2: y *= 2.
fn apply_step_forward(
    circ: &mut B,
    i: usize,
    sub: &QubitId,
    swp: &QubitId,
    s2: &QubitId,
    t1: &QubitId,
    x_reg: &[QubitId],
    y_reg: &[QubitId],
) {
    let n = 256usize;

    // 1) if subtracted: y += x mod q. Apply cofactor add -> cout adder at the
    // schedule's k.
    let k = super::next_cout_k();
    let ffg = super::next_ffg();
    arith::controlled_mod_add_k(circ, sub, &x_reg[..n], &y_reg[..n], Some(k), Some(ffg));
    // 2) if swap: swap(x, y).
    for j in 0..n {
        circ.cswap(*swp, x_reg[j], y_reg[j]);
    }
    // 3) y := 2*(1+s2)*y mod q. i==0: two separate controlled doublings (shift1 is
    //    t1-gated). i>0: the fused double+cdouble -- one combined (e+2d)*f fold
    //    (the unfused form costs extra inversion Toffoli).
    if i == 0 {
        controlled_mod_double(circ, t1, y_reg);
        controlled_mod_double(circ, s2, y_reg);
    } else {
        super::fused::fused_double_cdouble(circ, s2, y_reg);
    }
}

/// One inverse apply step (the divide body, gate-inverse of [`apply_step_forward`])
/// at iter `i`:  s_2 halve; shift1 halve; if swp: swap(x,y); if sub: y -= x mod q.
fn apply_step_reverse(
    circ: &mut B,
    i: usize,
    sub: &QubitId,
    swp: &QubitId,
    s2: &QubitId,
    t1: &QubitId,
    x_reg: &[QubitId],
    y_reg: &[QubitId],
) {
    let n = 256usize;

    // inverse of 3): i==0 two separate reverse-doublings (s2 halve then t1 halve);
    // i>0 the fused inverse double+cdouble.
    if i == 0 {
        controlled_mod_double_reverse(circ, s2, y_reg);
        controlled_mod_double_reverse(circ, t1, y_reg);
    } else {
        super::fused::fused_double_cdouble_reverse(circ, s2, y_reg);
    }
    // inverse of 2): swap (involutory).
    for j in 0..n {
        circ.cswap(*swp, x_reg[j], y_reg[j]);
    }
    // inverse of 1): y -= x mod q. The apply-path operands carry pseudo-Mersenne
    // representation drift, so the borrow clean uses the MBU form.
    // Apply cofactor sub -> schedule cout k.
    let k = super::next_cout_k();
    controlled_mod_sub_vented(circ, sub, &x_reg[..n], &y_reg[..n], Some(k));
}

/// Apply-path `y := y - ctrl * x (mod q)` whose borrow ancilla is cleaned by
/// MBU (hmr + cz_if_bit) rather than the normal comparator
/// `zero_and_free`. Tolerates pseudo-Mersenne representation drift (operands in
/// `[0, 2^256)` not strictly `< q`), which the apply pipeline produces. A gated
/// add over the full width with an MSBs phase-correction borrow clean.
fn controlled_mod_sub_vented(circ: &mut B, ctrl: &QubitId, x: &[QubitId], y: &[QubitId], sched_k: Option<usize>) {
    let n = x.len();
    assert_eq!(y.len(), n, "x,y equal width");
    let f_bytes = F_SECP256K1.to_le_bytes();
    let anc = circ.alloc_qubit();
    // X-sandwich: ~y += ctrl*x => y -= ctrl*x; cout = ctrl AND borrow.
    for q in y {
        circ.x(*q);
    }
    controlled_add_active_cout(circ, ctrl, x, y, &anc, sched_k);
    for q in y {
        circ.x(*q);
    }
    // gated -f fold on the borrow.
    for q in &y[..arith::LSBS] {
        circ.x(*q);
    }
    let ffg = super::next_ffg();
    arith::add_f_window_pub(circ, &anc, y, arith::LSBS, &f_bytes, Some(ffg));
    for q in &y[..arith::LSBS] {
        circ.x(*q);
    }
    // Measurement-vented borrow clean by MBU over the top `msbs` bits (a vented
    // adaptive msbs comparator, not a normal 2k MAJ/un-MAJ). The borrow == ctrl AND
    // carryout(y_final + x + 1); the -f fold only touched the low LSBS bits, so the
    // predicate is recomputed on the top msbs, which the fold left intact (msbs =
    // PAD). HMR the borrow, then on the fired shots recompute the carry as a deferred
    // Z via the Gidney `a + ~b + ~cin` carry chain (~1 Toffoli/bit); cz_if_bit(ctrl,
    // carry) cancels the phase. This never asserts |0>, so operand drift is tolerated.
    // Flip x_top so the comparator's internal `~b` yields `+x_top` (carryout(y+x+1)).
    let k = arith::MSBS.min(n);
    let lo = n - k;
    let ctrl = *ctrl;
    let bit = circ.alloc_bit();
    circ.hmr(anc, bit);
    circ.zero_and_free(anc);
    circ.push_condition(bit);
    let yt: Vec<QubitId> = y[lo..n].to_vec();
    let xt: Vec<QubitId> = x[lo..n].to_vec();
    for q in &xt {
        circ.x(*q); // ~b = ~(~x_top) = x_top inside the comparator -> carryout(y+x+1)
    }
    // Recompute the borrow through the headroom-adaptive (chunked) comparator
    // backend (`compare_geq_chunked_middle`); held-carry count = the full window
    // `k`. `flag = [yt >= xt] = carryout(yt + x + 1)` (xt was complemented to ~x
    // above) = the borrow predicate; deposit Z^(ctrl AND flag). The chunked
    // comparator builds the [>=] carry with an implicit +1 carry-in, so no separate
    // `zcin` is needed (this also drops a qubit vs the explicit carry-in form).
    let flag = circ.alloc_qubit();
    super::comparator::compare_geq_chunked_middle(circ, &yt, &xt, &flag, |c, fl| {
        c.cz(ctrl, *fl);
    }, k);
    circ.zero_and_free(flag);
    for q in &xt {
        circ.x(*q);
    }
    circ.pop_condition();
}

/// `y += ctrl * x` over `n` bits depositing carry-out into `cout` (= ctrl AND
/// overflow). Measurement-vented chunked-gated adder (~2.5 Toffoli/bit, bounded
/// peak) -- the inverse apply's mod-sub register add, on the peak-bound hot path.
/// `cout` is |0> on entry; `x` restored.
fn controlled_add_active_cout(circ: &mut B, ctrl: &QubitId, x: &[QubitId], y: &[QubitId], cout: &QubitId, sched_k: Option<usize>) {
    match sched_k {
        Some(k) => {
            // cout adder at the schedule's k (a = target = y, b = addend = x).
            let yr: Vec<&QubitId> = y.iter().collect();
            let xr: Vec<&QubitId> = x.iter().collect();
            super::gidney::controlled_hybrid_add_cout_refs(circ, ctrl, &yr, &xr, cout, k);
        }
        None => arith::controlled_add_vented_chunked_cout(circ, ctrl, x, y, arith::APPLY_CHUNK, Some(cout)),
    }
}


/// Jump-GCD in-place modular inversion / multiply: `(xv, y) -> (xv, y*x^{±1})`.
/// `Direction::Inverse` => `y := y * xv^-1 mod q`; `Direction::Forward` => `y := y * xv mod q`.
/// xv restored. Fused both directions (apply folded into the GCD passes).
#[must_use]
pub fn mod_mul_inverse_in_place(
    circ: &mut B,
    mut xv: Vec<QubitId>,
    y: &[QubitId],
    dir: Direction,
) -> Vec<QubitId> {
    let n = 256usize;
    assert_eq!(xv.len(), n, "xv must be 256 bits");
    assert_eq!(y.len(), n, "y must be 256 bits");

    match dir {
        Direction::Inverse => {
            // Fused inverse apply folds into the forward GCD. apply pair (x_reg=y,
            // y_reg=tmp); result y = z*x^-1, tmp = 0.
            let tmp: Vec<QubitId> = (0..n).map(|_| circ.alloc_qubit()).collect();
            for j in 0..n {
                circ.swap(y[j], tmp[j]); // tmp = z, y = 0
            }
            let mut tape = forward_gcd_jump(circ, &mut xv, Some((y, &tmp)));
            for q in tmp {
                circ.zero_and_free(q);
            }
            reverse_gcd_jump(circ, &mut xv, &mut tape, None);
            xv
        }
        Direction::Forward => {
            // Fused forward apply folds into the reverse GCD. apply pair (x_reg=tmp=z,
            // y_reg=y=0); result y = z*x, tmp = canonical 0 (drift cleared).
            let mut tape = forward_gcd_jump(circ, &mut xv, None);
            let tmp: Vec<QubitId> = (0..n).map(|_| circ.alloc_qubit()).collect();
            for j in 0..n {
                circ.swap(y[j], tmp[j]); // tmp = z, y = 0
            }
            reverse_gcd_jump(circ, &mut xv, &mut tape, Some((&tmp, y)));
            clear_zeroed_drift(circ, &tmp[..n]);
            for q in tmp {
                circ.zero_and_free(q);
            }
            xv
        }
    }
}

/// Clear a register holding a value `== 0 (mod q)` whose representative is `0` or
/// `q`: XOR q maps `q -> 0` (caught by the caller's zero_and_free if not).
fn clear_zeroed_drift(circ: &mut B, reg: &[QubitId]) {
    let q_bytes = q_secp256k1_le();
    for (i, qb) in reg.iter().enumerate() {
        if (q_bytes.get(i / 8).copied().unwrap_or(0) >> (i % 8)) & 1 == 1 {
            circ.x(*qb);
        }
    }
}

