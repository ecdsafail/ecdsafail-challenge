//! Safegcd v3f: MSB-wrap detection + wrap_log tracking for d/e.
//!
//! Bug in v3e/v2-s5: NOT(e)+d+1 = d-e (mod 2^NX). When d < e, wraps to
//! 2^NX+d-e, and arithmetic shift + div2_mod p-add leaves a 2^(NX-1) offset.
//!
//! Fix: after NOT(e)+d+1, check MSB=e[NX-1]. If 1 (wrap), add p:
//!   (2^NX+d-e+p) mod 2^NX = d-e+p ∈ (0,p).  After fix, both cases have e∈[0,p).
//!
//! Wrap flag stored per-step in wrap_log.  In reverse: undo wrap fix,
//! recompute e[NX-1] at same phase, XOR to clean.  Carry ancillae from
//! the wrap add/sub are cleaned via CX(wrap, carry) since carry==wrap
//! deterministically (the add wraps, the sub borrows).
use super::*;
use crate::point_add::N;
use alloy_primitives::U256;

const DEFAULT_SAFEGCD_ITERS: usize = 741;
const DB: usize = 12;
const NX: usize = N + 14;           // 270 bits
const LOG_BITS: usize = 3;          // b0, b1, e_odd per step

// ── public config accessors ──────────────────────────────────────────

pub(crate) fn safegcd_iters() -> usize {
    std::env::var("SAFEGCD_ITERS").ok().and_then(|s| s.parse().ok()).unwrap_or(DEFAULT_SAFEGCD_ITERS)
}
pub(crate) fn safegcd_enabled() -> bool { std::env::var("SAFEGCD").ok().as_deref() == Some("1") }

// ── low-level helpers ────────────────────────────────────────────────

fn cswap(b: &mut B, ctrl: QubitId, a: QubitId, bb: QubitId) {
    if a == bb { return; }
    b.cx(bb, a);
    b.ccx(ctrl, a, bb);
    b.cx(bb, a);
}
fn load_u256_ext(b: &mut B, reg: &[QubitId], val: U256) {
    for i in 0..N { if val.bit(i) { b.x(reg[i]); } }
}
fn cnot_reg(b: &mut B, ctrl: QubitId, reg: &[QubitId]) {
    for &q in reg { b.cx(ctrl, q); }
}

/// Arithmetic shift right (preserves MSB). sign_anc: scratch qubit, |0⟩ initially.
fn arith_shift_right_w(b: &mut B, reg: &[QubitId], sign_anc: QubitId, w: usize) {
    b.cx(reg[w - 1], sign_anc);
    for i in 0..w - 1 { b.swap(reg[i], reg[i + 1]); }
    b.swap(sign_anc, reg[w - 1]);
}
fn arith_unshift_right_w(b: &mut B, reg: &[QubitId], sign_anc: QubitId, w: usize) {
    b.swap(sign_anc, reg[w - 1]);
    for i in (0..w - 1).rev() { b.swap(reg[i], reg[i + 1]); }
    b.cx(reg[w - 1], sign_anc);
}

// ── delta control ────────────────────────────────────────────────────

fn delta_inc(b: &mut B, d: &[QubitId]) {
    let c = b.alloc_qubits(DB - 1);
    b.x(d[0]); b.cx(d[0], c[0]);
    for i in 1..DB - 1 { b.ccx(d[i], c[i - 1], c[i]); }
    for i in 1..DB { b.cx(c[i - 1], d[i]); }
    for i in (1..DB - 1).rev() { b.ccx(d[i], c[i - 1], c[i]); }
    b.cx(d[0], c[0]); b.free_vec(&c);
}
fn delta_dec(b: &mut B, d: &[QubitId]) {
    for i in 0..DB { b.x(d[i]); } delta_inc(b, d); for i in 0..DB { b.x(d[i]); }
}
fn cdelta_inc(b: &mut B, d: &[QubitId], ctrl: QubitId) {
    let anc = b.alloc_qubit();
    let c = b.alloc_qubits(DB - 1);
    b.cx(ctrl, d[0]);
    b.ccx(ctrl, d[0], c[0]);
    for i in 1..DB - 1 {
        b.ccx(d[i], c[i - 1], anc);
        b.ccx(ctrl, anc, c[i]);
        b.ccx(d[i], c[i - 1], anc);
    }
    for i in 1..DB { b.ccx(ctrl, c[i - 1], d[i]); }
    for i in (1..DB - 1).rev() {
        b.ccx(d[i], c[i - 1], anc);
        b.ccx(ctrl, anc, c[i]);
        b.ccx(d[i], c[i - 1], anc);
    }
    b.ccx(ctrl, d[0], c[0]);
    b.free(anc); b.free_vec(&c);
}
fn cdelta_dec(b: &mut B, d: &[QubitId], ctrl: QubitId) {
    for i in 0..DB { b.cx(ctrl, d[i]); }
    cdelta_inc(b, d, ctrl);
    for i in 0..DB { b.cx(ctrl, d[i]); }
}
fn delta_nz_compute(b: &mut B, delta: &[QubitId]) -> (QubitId, Vec<QubitId>) {
    let s = b.alloc_qubits(12); let nz = s[11];
    for p in 0..6 {
        let a = delta[2 * p]; let bb = delta[2 * p + 1]; let c = s[p];
        b.cx(a, c); b.cx(bb, c); b.ccx(a, bb, c);
    }
    for p in 0..3 {
        let a = s[2 * p]; let bb = s[2 * p + 1]; let c = s[6 + p];
        b.cx(a, c); b.cx(bb, c); b.ccx(a, bb, c);
    }
    b.cx(s[6], s[9]); b.cx(s[7], s[9]); b.ccx(s[6], s[7], s[9]);
    b.cx(s[8], s[10]);
    b.cx(s[9], nz); b.cx(s[10], nz); b.ccx(s[9], s[10], nz);
    (nz, s)
}
fn delta_nz_uncompute(b: &mut B, delta: &[QubitId], s: &[QubitId]) {
    let nz = s[11];
    b.ccx(s[9], s[10], nz); b.cx(s[10], nz); b.cx(s[9], nz);
    b.cx(s[8], s[10]);
    b.ccx(s[6], s[7], s[9]); b.cx(s[7], s[9]); b.cx(s[6], s[9]);
    for p in (0..3).rev() {
        let a = s[2 * p]; let bb = s[2 * p + 1]; let c = s[6 + p];
        b.ccx(a, bb, c); b.cx(bb, c); b.cx(a, c);
    }
    for p in (0..6).rev() {
        let a = delta[2 * p]; let bb = delta[2 * p + 1]; let c = s[p];
        b.ccx(a, bb, c); b.cx(bb, c); b.cx(a, c);
    }
}

// ═══════════════════════════════════════════════════════════════════════
//  Safegcd divsteps: forward pass
// ═══════════════════════════════════════════════════════════════════════

pub(crate) fn emit_safegcd_tobitvector_steps(
    b: &mut B,
    u: &[QubitId], v: &[QubitId],
    d: &[QubitId], e: &[QubitId],
    dialog_log: &[QubitId],
    wrap_log: &[QubitId],
    have_de: bool,
    p: U256,
    sign_log: &[QubitId],
) -> Vec<QubitId> {
    let iters = safegcd_iters();
    assert!(dialog_log.len() >= LOG_BITS * iters);
    assert_eq!(wrap_log.len(), iters);
    if have_de { assert_eq!(d.len(), NX); assert_eq!(e.len(), NX); }
    assert_eq!(u.len(), NX); assert_eq!(v.len(), NX);

    let one = b.alloc_qubits(NX); b.x(one[0]);
    let delta = b.alloc_qubits(DB); b.x(delta[0]);

    for step in 0..iters {
        let b0 = dialog_log[LOG_BITS * step];
        let b1 = dialog_log[LOG_BITS * step + 1];
        let e_odd = if have_de { dialog_log[LOG_BITS * step + 2] } else { QubitId(0) };

        // ── b0 = v[0] ──
        b.cx(v[0], b0);

        // ── b1 = (delta > 0) ∧ b0 ──
        let (dnz, dscr) = delta_nz_compute(b, &delta);
        let ds = b.alloc_qubit(); b.cx(delta[DB - 1], ds); b.x(ds);
        let t1 = b.alloc_qubit();
        b.ccx(ds, dnz, t1); b.ccx(t1, b0, b1); b.ccx(ds, dnz, t1);
        b.free(t1);
        b.x(ds); b.cx(delta[DB - 1], ds); b.free(ds);
        delta_nz_uncompute(b, &delta, &dscr);

        // ── CSWAP u↔v, d↔e when b1 ──
        for i in 0..NX { cswap(b, b1, u[i], v[i]); }
        if have_de { for i in 0..NX { cswap(b, b1, d[i], e[i]); } }

        // ── CNOT when b1 ──
        cnot_reg(b, b1, &v[..NX]);
        if have_de { cnot_reg(b, b1, &e[..NX]); }

        // ── b1b0 = b1 ∧ b0 (Branch 1 control) ──
        let b1b0 = b.alloc_qubit(); b.ccx(b1, b0, b1b0);

        // ── Add 1 when b1∧b0 (part of NOT+add+1) ──
        let ci1 = b.alloc_qubit(); let s1 = b.alloc_qubit();
        cuccaro_add_ctrl_lowq(b, &one, &v[..NX], b1b0, ci1, s1);
        b.free(s1); b.free(ci1);
        if have_de {
            let cd = b.alloc_qubit(); let sd = b.alloc_qubit();
            cuccaro_add_ctrl_lowq(b, &one, &e[..NX], b1b0, cd, sd);
            b.free(sd); b.free(cd);
        }

        // ── Add u to v, d to e when b0 ──
        let ci = b.alloc_qubit(); let sc = b.alloc_qubit();
        cuccaro_add_ctrl_lowq(b, &u[..NX], &v[..NX], b0, ci, sc);
        b.free(sc); b.free(ci);
        if have_de {
            let cd = b.alloc_qubit(); let sd = b.alloc_qubit();
            cuccaro_add_ctrl_lowq(b, &d[..NX], &e[..NX], b0, cd, sd);
            b.free(sd); b.free(cd);
        }

        // ── v3f: Wrap detection and fix for e ──
        // After NOT(e)+1+d: e = d_old - e_old (mod 2^NX).
        // e[NX-1] = 1 iff wrap (d_old < e_old in the sub).  Only meaningful for Branch 1.
        if have_de {
            let wrap_bit = wrap_log[step];
            b.ccx(e[NX - 1], b1b0, wrap_bit);

            // If wrap: e += p.  The add wraps mod 2^NX, carry=1 deterministically
            // because e ≈ 2^NX and p > 0.  We clean the carry with CX(wrap, carry).
            let p_reg = b.alloc_qubits(NX); load_u256_ext(b, &p_reg, p);
            let cw = b.alloc_qubit(); let sw = b.alloc_qubit();
            cuccaro_add_ctrl_lowq(b, &p_reg, &e[..NX], wrap_bit, cw, sw);
            b.cx(wrap_bit, cw);  // cw ← cw⊕wrap = carry⊕wrap = wrap⊕wrap = 0
            b.free(sw); b.free(cw);
            load_u256_ext(b, &p_reg, p); b.free_vec(&p_reg);
        }

        b.free(b1b0);

        // ── e: div2_mod ──
        if have_de {
            b.cx(e[0], e_odd);

            let p_reg = b.alloc_qubits(NX); load_u256_ext(b, &p_reg, p);
            let ce = b.alloc_qubit(); let se = b.alloc_qubit();
            cuccaro_add_ctrl_lowq(b, &p_reg, &e[..NX], e_odd, ce, se);
            b.free(se); b.free(ce);
            load_u256_ext(b, &p_reg, p); b.free_vec(&p_reg);

            let e_sign = b.alloc_qubit();
            arith_shift_right_w(b, e, e_sign, NX);
            b.free(e_sign);
        }

        // ── v: arithmetic shift right ──
        arith_shift_right_w(b, v, sign_log[step], NX);

        // ── Delta update ──
        let b1s = b.alloc_qubit(); b.cx(b1, b1s);
        for i in 0..DB { b.cx(b1s, delta[i]); }
        delta_inc(b, &delta);
        cdelta_inc(b, &delta, b1s);
        b.cx(b1, b1s); b.free(b1s);
    }

    b.x(one[0]); b.free_vec(&one);
    delta
}

// ═══════════════════════════════════════════════════════════════════════
//  Safegcd divsteps: reverse (uncompute) pass
// ═══════════════════════════════════════════════════════════════════════

pub(crate) fn emit_safegcd_tobitvector_steps_reverse(
    b: &mut B,
    u: &[QubitId], v: &[QubitId],
    d: &[QubitId], e: &[QubitId],
    dialog_log: &[QubitId],
    wrap_log: &[QubitId],
    have_de: bool,
    p: U256,
    sign_log: &[QubitId],
    delta: &[QubitId],
) {
    let iters = safegcd_iters();
    assert!(dialog_log.len() >= LOG_BITS * iters);
    assert_eq!(wrap_log.len(), iters);
    if have_de { assert_eq!(d.len(), NX); assert_eq!(e.len(), NX); }
    assert_eq!(u.len(), NX); assert_eq!(v.len(), NX);
    assert_eq!(delta.len(), DB);

    let one = b.alloc_qubits(NX); b.x(one[0]);

    for step in (0..iters).rev() {
        let b0 = dialog_log[LOG_BITS * step];
        let b1 = dialog_log[LOG_BITS * step + 1];
        let e_odd = if have_de { dialog_log[LOG_BITS * step + 2] } else { QubitId(0) };

        // ── Compute b1b0 early (needed for wrap cleanup AND 1-sub undo) ──
        let b1b0 = b.alloc_qubit(); b.ccx(b1, b0, b1b0);

        // ── Reverse: v unshift ──
        arith_unshift_right_w(b, v, sign_log[step], NX);

        // ── Reverse: delta ──
        let b1s = b.alloc_qubit(); b.cx(b1, b1s);
        cdelta_dec(b, &delta, b1s);
        delta_dec(b, &delta);
        for i in 0..DB { b.cx(b1s, delta[i]); }
        b.cx(b1, b1s); b.free(b1s);

        // ── Reverse: e ──
        if have_de {
            // Unshift
            let e_sign = b.alloc_qubit();
            b.cx(e[NX - 1], e_sign);
            arith_unshift_right_w(b, e, e_sign, NX);
            b.free(e_sign);

            // Undo e_odd p-add
            let p_reg = b.alloc_qubits(NX); load_u256_ext(b, &p_reg, p);
            let ce = b.alloc_qubit(); let se = b.alloc_qubit();
            cuccaro_sub_ctrl_lowq(b, &p_reg, &e[..NX], e_odd, ce, se);
            b.free(se); b.free(ce);
            load_u256_ext(b, &p_reg, p); b.free_vec(&p_reg);

            // Uncompute e_odd BEFORE wrap undo (e[0] must match forward's
            // e[0] which was sampled AFTER the wrap fix, because p is odd).
            b.cx(e[0], e_odd);

            // v3f: Undo wrap fix (subtract p if wrap_log[step])
            let p_reg2 = b.alloc_qubits(NX); load_u256_ext(b, &p_reg2, p);
            let cw = b.alloc_qubit(); let sw = b.alloc_qubit();
            cuccaro_sub_ctrl_lowq(b, &p_reg2, &e[..NX], wrap_log[step], cw, sw);
            b.cx(wrap_log[step], cw);  // carry←0 (borrow from wrap-fix undo)
            b.free(sw); b.free(cw);
            load_u256_ext(b, &p_reg2, p); b.free_vec(&p_reg2);

            // v3f: Clean wrap_log[step] by recomputing e[NX-1] at same phase.
            // After undoing wrap fix: e = d_old - e_old (mod 2^NX).
            // e[NX-1] = (e_after_swap < d_after_swap) = forward's wrap condition.
            // CCX(e[NX-1], b1b0, wrap_log) XORs the same AND back, clearing it.
            b.ccx(e[NX - 1], b1b0, wrap_log[step]);
        }

        // ── Reverse: sub u from v (undo add with b0) ──
        let ci = b.alloc_qubit(); let sc = b.alloc_qubit();
        cuccaro_sub_ctrl_lowq(b, &u[..NX], &v[..NX], b0, ci, sc);
        b.free(sc); b.free(ci);
        if have_de {
            let cd = b.alloc_qubit(); let sd = b.alloc_qubit();
            cuccaro_sub_ctrl_lowq(b, &d[..NX], &e[..NX], b0, cd, sd);
            b.free(sd); b.free(cd);
        }

        // ── Reverse: sub 1 (undo add with b1b0) ──
        let ci1 = b.alloc_qubit(); let s1 = b.alloc_qubit();
        cuccaro_sub_ctrl_lowq(b, &one, &v[..NX], b1b0, ci1, s1);
        b.free(s1); b.free(ci1);
        if have_de {
            let cd = b.alloc_qubit(); let sd = b.alloc_qubit();
            cuccaro_sub_ctrl_lowq(b, &one, &e[..NX], b1b0, cd, sd);
            b.free(sd); b.free(cd);
        }

        // Clean and free b1b0
        b.ccx(b1, b0, b1b0); b.free(b1b0);

        // ── Reverse: CNOT (undo cnot_reg when b1) ──
        if have_de { cnot_reg(b, b1, &e[..NX]); }
        cnot_reg(b, b1, &v[..NX]);

        // ── Reverse: CSWAP (undo cswap when b1) ──
        for i in 0..NX { cswap(b, b1, u[i], v[i]); }
        if have_de { for i in 0..NX { cswap(b, b1, d[i], e[i]); } }

        // ── Reverse: b1 uncompute ──
        let (dnz, dscr) = delta_nz_compute(b, &delta);
        let ds = b.alloc_qubit(); b.cx(delta[DB - 1], ds); b.x(ds);
        let t1 = b.alloc_qubit();
        b.ccx(ds, dnz, t1); b.ccx(t1, b0, b1); b.ccx(ds, dnz, t1);
        b.free(t1);
        b.x(ds); b.cx(delta[DB - 1], ds); b.free(ds);
        delta_nz_uncompute(b, &delta, &dscr);

        // ── Reverse: b0 uncompute ──
        b.cx(v[0], b0);
    }

    b.x(delta[0]);
    b.x(one[0]); b.free_vec(&one);
}

// ═══════════════════════════════════════════════════════════════════════
//  Apply inverse: extract d[0..255] → target (with u_neg correction)
// ═══════════════════════════════════════════════════════════════════════

fn safegcd_apply_inverse(
    b: &mut B, target: &[QubitId], d: &[QubitId], u: &[QubitId], p: U256,
) {
    let result = b.alloc_qubits(N);
    for i in 0..N { b.cx(d[i], result[i]); }

    let u_neg = b.alloc_qubit(); b.cx(u[NX - 1], u_neg);
    cnot_reg(b, u_neg, &result[..N]);
    let p1c = b.alloc_qubits(N); load_u256_ext(b, &p1c, p + U256::from(1u64));
    let ci = b.alloc_qubit(); let sc = b.alloc_qubit();
    cuccaro_add_ctrl_lowq(b, &p1c, &result[..N], u_neg, ci, sc);
    b.free(sc); b.free(ci);
    load_u256_ext(b, &p1c, p + U256::from(1u64)); b.free_vec(&p1c);
    b.cx(u[NX - 1], u_neg); b.free(u_neg);

    for i in 0..N { b.swap(result[i], target[i]); }
    b.free_vec(&result);
}

// ═══════════════════════════════════════════════════════════════════════
//  Main quotient: safegcd inverse of factor, placed into target
// ═══════════════════════════════════════════════════════════════════════

pub(crate) fn emit_safegcd_raw_quotient(
    b: &mut B, factor: &[QubitId], target: &[QubitId], p: U256,
) {
    let u = b.alloc_qubits(NX); load_u256_ext(b, &u, p);
    let v = b.alloc_qubits(NX);
    for i in 0..N { b.swap(v[i], factor[i]); }
    let have_de = true;
    let d = b.alloc_qubits(NX);
    let e = b.alloc_qubits(NX); b.x(e[0]);
    let dialog_log = b.alloc_qubits(LOG_BITS * safegcd_iters());
    let wrap_log = b.alloc_qubits(safegcd_iters());
    let sign_log = b.alloc_qubits(safegcd_iters());

    let delta = emit_safegcd_tobitvector_steps(
        b, &u, &v, &d, &e, &dialog_log, &wrap_log, have_de, p, &sign_log,
    );

    safegcd_apply_inverse(b, target, &d, &u, p);

    emit_safegcd_tobitvector_steps_reverse(
        b, &u, &v, &d, &e, &dialog_log, &wrap_log, have_de, p, &sign_log, &delta,
    );

    b.free_vec(&delta);
    b.free_vec(&wrap_log);
    b.free_vec(&sign_log);
    b.x(e[0]); b.free_vec(&e); b.free_vec(&d);
    for i in 0..N { b.swap(v[i], factor[i]); } b.free_vec(&v);
    load_u256_ext(b, &u, p); b.free_vec(&u);
    b.free_vec(&dialog_log);
}

// ═══════════════════════════════════════════════════════════════════════
//  Point-addition: safegcd quotient + dialog IPMUL + fused-round84
// ═══════════════════════════════════════════════════════════════════════

pub(crate) fn emit_safegcd_raw_pa(
    b: &mut B, tx: &[QubitId], ty: &[QubitId], ox: &[BitId], oy: &[BitId], p: U256,
) {
    mod_sub_qb(b, tx, ox, p);
    mod_sub_qb(b, ty, oy, p);
    emit_safegcd_raw_quotient(b, tx, ty, p);

    round84_emit_fused_square_xtail(b, tx, ty, ox, p);

    if dialog_fuse_c_form_enabled() {
        mod_add_triple_qb(b, tx, ox, p);
    } else {
        mod_sub_qb(b, tx, ox, p);
        mod_neg_inplace_fast(b, tx, p);
    }

    emit_dialog_gcd_raw_ipmul(b, tx, ty, p);

    mod_sub_qb(b, ty, oy, p);

    if dialog_fuse_x_restore_enabled() {
        mod_const_minus_reg_qb(b, tx, ox, p);
    } else {
        mod_neg_inplace_fast(b, tx, p);
        mod_add_qb(b, tx, ox, p);
    }
}
