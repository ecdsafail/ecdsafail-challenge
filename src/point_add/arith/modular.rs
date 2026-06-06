//! Modular arithmetic: add/sub/neg (qq and qb), doubling/halving, shifts, and
//! the controlled (cmod_*) variants. All operate over secp256k1's prime via
//! Solinas reduction on "extended" (n+1)-wide registers; bit n is a transient
//! overflow/sign ancilla allocated for the duration of a mod-op.
use super::*;

/// `acc := (acc + a) mod p`. Both `acc` and `a` are n-bit quantum registers
/// with value in [0, p). Solinas reduction using c = 2^n - p: sum ∈ [0, 2p),
/// then add c, branch on top bit to either clear it (reduction) or undo
/// the add (no reduction). Saves one full (n+1)-wide Cuccaro compared to
/// the sub-p/add-p/csub-p pattern.
pub(crate) fn mod_add_qq(b: &mut B, acc: &[QubitId], a: &[QubitId], p: U256) {
    let n = acc.len();
    assert_eq!(n, a.len());
    debug_assert_eq!(n, 256);

    let (acc_ext, acc_ovf) = ext_reg(b, acc);
    let (a_ext, a_ovf) = ext_reg(b, a);

    // Step 1: (n+1)-bit add. acc_ext ∈ [0, 2p).
    add_nbit_qq(b, &a_ext, &acc_ext);

    // Step 2: add c. If sum was >= p, the top bit of (sum + c) becomes 1.
    let c = U256::MAX.wrapping_sub(p).wrapping_add(U256::from(1));
    add_nbit_const(b, &acc_ext, c);

    // Step 3: flag := acc_ovf (= top bit of sum + c).
    let flag = b.alloc_qubit();
    b.cx(acc_ovf, flag);

    // Step 4: if flag=0 (no reduction needed), undo the add of c.
    b.x(flag);
    csub_nbit_const(b, &acc_ext, c, flag);
    b.x(flag);

    // Step 5: if flag=1, clear the top bit (drops 2^n → yields sum - p).
    b.cx(flag, acc_ovf);

    // Step 6: uncompute flag. Same identity as the old version:
    //   flag == (acc_final < a_orig)
    // because in the flag=1 case acc_final = acc_orig + a - p < a (since acc_orig < p),
    // and in the flag=0 case acc_final = acc_orig + a ≥ a.
    cmp_lt_into(b, &acc_ext[..n], &a_ext[..n], flag);
    b.free(flag);

    unext_reg(b, a_ovf);
    unext_reg(b, acc_ovf);
    let _ = (acc_ext, a_ext);
}

pub(crate) fn mod_sub_qq(b: &mut B, acc: &[QubitId], a: &[QubitId], p: U256) {
    // mod_add_qq is a bijection on (acc, a): (acc, a) ↦ (acc + a mod p, a).
    // Its gate-level inverse therefore acts as (acc, a) ↦ (acc - a mod p, a),
    // which is exactly what we want. emit_inverse replays the forward's gates
    // reversed, skipping R markers — valid because mod_add_qq is clean
    // (every ancilla is driven to |0⟩ before its R).
    let a_copy: Vec<QubitId> = a.to_vec();
    emit_inverse(b, move |b| mod_add_qq(b, acc, &a_copy, p));
}

/// Ancilla-light copy of [`mod_add_qq`]. The only difference: the two Solinas
/// constant corrections (`+c` in step 2 and the conditional `-c` in step 4) use
/// the extended-carry clean adders, which load `c = 2^256 - p` into a 256-qubit
/// register (= `acc_ext.len() - 1`) and fold the overflow into `acc_ext[n]` via a
/// measurement-free Cuccaro. The stock `mod_add_qq` instead calls
/// `add_nbit_const`/`csub_nbit_const`, which materialize a full 257-qubit loaded
/// constant — the sole +1 transient that pins the round84 mid-sub peak at 1309.
/// Replacing it with the 256-wide load drops that peak to 1308. Value- and
/// phase-identical to `mod_add_qq`; clean (no `b.hmr`), so `emit_inverse`-safe.
pub(crate) fn mod_add_qq_lowq(b: &mut B, acc: &[QubitId], a: &[QubitId], p: U256) {
    let n = acc.len();
    assert_eq!(n, a.len());
    debug_assert_eq!(n, 256);

    let (acc_ext, acc_ovf) = ext_reg(b, acc);
    let (a_ext, a_ovf) = ext_reg(b, a);

    // Step 1: (n+1)-bit add. acc_ext ∈ [0, 2p).
    add_nbit_qq(b, &a_ext, &acc_ext);

    // The addend `a_ext` is preserved by every step below (the only later read of
    // it, the step-6 compare, touches a_ext[..n] only), so `a_ovf = a_ext[n]` is
    // provably |0> and idle during the two const corrections (steps 2 & 4). Under
    // the borrow flag we lend it to the Cuccaro as the carry-in slot, removing the
    // sole fresh +1 transient that pins the round84-lowq mid-sub peak at 1308.
    let borrow = if r84_lowq_cin_borrow_enabled() {
        Some(a_ovf)
    } else {
        None
    };

    // Step 2: add c (ancilla-light: 256-wide const load + clean carry capture).
    let c = U256::MAX.wrapping_sub(p).wrapping_add(U256::from(1));
    add_nbit_const_extcarry_clean_with_cin(b, &acc_ext, c, borrow);

    // Step 3: flag := acc_ovf (= top bit of sum + c).
    let flag = b.alloc_qubit();
    b.cx(acc_ovf, flag);

    // Step 4: if flag=0 (no reduction needed), undo the add of c.
    b.x(flag);
    csub_nbit_const_extcarry_clean_with_cin(b, &acc_ext, c, flag, borrow);
    b.x(flag);

    // Step 5: if flag=1, clear the top bit (drops 2^n → yields sum - p).
    b.cx(flag, acc_ovf);

    // Step 6: uncompute flag (same identity as mod_add_qq).
    cmp_lt_into(b, &acc_ext[..n], &a_ext[..n], flag);
    b.free(flag);

    unext_reg(b, a_ovf);
    unext_reg(b, acc_ovf);
    let _ = (acc_ext, a_ext);
}

/// Ancilla-light `acc := (acc - a) mod p`. Exact gate-level inverse of
/// [`mod_add_qq_lowq`] (which is clean), so `emit_inverse` replays it as
/// `(acc, a) ↦ (acc - a mod p, a)` with operand preserved and zero phase.
pub(crate) fn mod_sub_qq_lowq(b: &mut B, acc: &[QubitId], a: &[QubitId], p: U256) {
    let a_copy: Vec<QubitId> = a.to_vec();
    emit_inverse(b, move |b| mod_add_qq_lowq(b, acc, &a_copy, p));
}

/// Fast `acc := (acc - a) mod p`. Direct sub + conditional add-p + flag
/// uncompute via neg+cmp_lt+neg. All ops use measurement-based Cuccaro.
pub(crate) fn mod_sub_qq_fast(b: &mut B, acc: &[QubitId], a: &[QubitId], p: U256) {
    let n = acc.len();
    assert_eq!(n, a.len());
    debug_assert_eq!(n, 256);

    let (acc_ext, acc_ovf) = ext_reg(b, acc);
    let (a_ext, a_ovf) = ext_reg(b, a);

    // Step 1: (n+1)-bit sub.
    sub_nbit_qq_fast(b, &a_ext, &acc_ext);

    // Step 2: flag = acc_ovf (=1 iff underflow, i.e. acc < a).
    let flag = b.alloc_qubit();
    b.cx(acc_ovf, flag);
    // We only need the borrow as a separate flag; the low register is
    // corrected modulo 2^n, so clear the extension bit immediately.
    b.cx(flag, acc_ovf);

    // Step 3: underflow correction. With p = 2^n - c, the wrapped 256-bit
    // subtraction needs only a conditional subtract of c on the low register.
    let c = U256::MAX.wrapping_sub(p).wrapping_add(U256::from(1));
    if kal_vent_modadd_enabled() {
        // Use venting cisub with a_ext as dirty qubits.
        let c_low = c.as_limbs()[0];
        let q_clean2: [QubitId; 2] = [b.alloc_qubit(), b.alloc_qubit()];
        venting::cisub_dirty_2clean_classical(
            b,
            &acc_ext[..n],
            &a_ext[..n - 2],
            &q_clean2,
            c_low,
            flag,
        );
        b.free(q_clean2[0]);
        b.free(q_clean2[1]);
    } else if secp_direct_const_arith_enabled() {
        csub_nbit_const_direct_fast(b, &acc_ext[..n], c, flag);
    } else {
        csub_nbit_const_fast(b, &acc_ext[..n], c, flag);
    }

    // Step 4: uncompute flag. Identity: flag = NOT(acc_final < (p - a)).
    // Negate a in place, compare, un-negate.
    b.x(flag);
    mod_neg_inplace_fast(b, &a_ext[..n], p);
    cmp_lt_into_fast(b, &acc_ext[..n], &a_ext[..n], flag);
    mod_neg_inplace_fast(b, &a_ext[..n], p);
    b.free(flag);

    unext_reg(b, a_ovf);
    unext_reg(b, acc_ovf);
    let _ = (acc_ext, a_ext);
}

/// Low-peak `acc := (acc + a) mod p`. Identical structure to `mod_add_qq` but
/// the two Solinas-constant corrections (`+c`, conditional `-c`) are vented onto
/// the operand `a_ext` as dirty scratch (2 clean qubits) instead of a fresh
/// n-qubit loaded-constant register. The main add and the flag-uncompute compare
/// stay ancilla-free (Cuccaro / cmp_lt_into), so the only transient is +2 clean.
/// Used inside the round84 Solinas reduction where the materialized `load_const`
/// coexisting with tmp_ext + z1_reg was the peak binder. `c = 2^256 - p` fits in
/// 64 bits, so `c_low` carries the whole constant.
pub(crate) fn mod_add_qq_vent(b: &mut B, acc: &[QubitId], a: &[QubitId], p: U256) {
    let n = acc.len();
    assert_eq!(n, a.len());
    debug_assert_eq!(n, 256);

    let (acc_ext, acc_ovf) = ext_reg(b, acc);
    let (a_ext, a_ovf) = ext_reg(b, a);

    add_nbit_qq(b, &a_ext, &acc_ext);

    let c = U256::MAX.wrapping_sub(p).wrapping_add(U256::from(1));
    let c_low = c.as_limbs()[0];
    let n1 = acc_ext.len();
    {
        let q_clean2: [QubitId; 2] = [b.alloc_qubit(), b.alloc_qubit()];
        venting::iadd_dirty_2clean_classical(
            b,
            &acc_ext,
            &a_ext[..n1 - 2],
            &q_clean2,
            c_low,
            false,
        );
        b.free(q_clean2[0]);
        b.free(q_clean2[1]);
    }

    let flag = b.alloc_qubit();
    b.cx(acc_ovf, flag);

    b.x(flag);
    {
        let q_clean2: [QubitId; 2] = [b.alloc_qubit(), b.alloc_qubit()];
        venting::cisub_dirty_2clean_classical(
            b,
            &acc_ext,
            &a_ext[..n1 - 2],
            &q_clean2,
            c_low,
            flag,
        );
        b.free(q_clean2[0]);
        b.free(q_clean2[1]);
    }
    b.x(flag);

    b.cx(flag, acc_ovf);

    cmp_lt_into(b, &acc_ext[..n], &a_ext[..n], flag);
    b.free(flag);

    unext_reg(b, a_ovf);
    unext_reg(b, acc_ovf);
    let _ = (acc_ext, a_ext);
}

/// `acc := (acc - a) mod p`, low-peak. Explicit gate-reverse of
/// `mod_add_qq_vent` (the venting protocols use measurement, so `emit_inverse`
/// cannot reverse them; each venting step is undone by its matched dual:
/// iadd↔isub, cisub↔ciadd). The flag-uncompute is `cmp_lt_into` (self-inverse,
/// no materialized neg), so no n-wide const register is ever live.
pub(crate) fn mod_sub_qq_vent(b: &mut B, acc: &[QubitId], a: &[QubitId], p: U256) {
    let n = acc.len();
    assert_eq!(n, a.len());
    debug_assert_eq!(n, 256);

    let (acc_ext, acc_ovf) = ext_reg(b, acc);
    let (a_ext, a_ovf) = ext_reg(b, a);

    let c = U256::MAX.wrapping_sub(p).wrapping_add(U256::from(1));
    let c_low = c.as_limbs()[0];
    let n1 = acc_ext.len();

    // Reverse of forward step 6: cmp_lt_into is its own inverse (XOR into flag).
    let flag = b.alloc_qubit();
    cmp_lt_into(b, &acc_ext[..n], &a_ext[..n], flag);

    // Reverse of step 5.
    b.cx(flag, acc_ovf);

    // Reverse of step 4: forward applied (cisub c) under !flag; undo with ciadd.
    b.x(flag);
    {
        let q_clean2: [QubitId; 2] = [b.alloc_qubit(), b.alloc_qubit()];
        venting::ciadd_dirty_2clean_classical(
            b,
            &acc_ext,
            &a_ext[..n1 - 2],
            &q_clean2,
            c_low,
            flag,
            false,
        );
        b.free(q_clean2[0]);
        b.free(q_clean2[1]);
    }
    b.x(flag);

    // Reverse of step 3.
    b.cx(acc_ovf, flag);
    b.free(flag);

    // Reverse of step 2: undo the unconditional (iadd c) with a cisub under an
    // always-on control.
    {
        let one = b.alloc_qubit();
        b.x(one);
        let q_clean2: [QubitId; 2] = [b.alloc_qubit(), b.alloc_qubit()];
        venting::cisub_dirty_2clean_classical(b, &acc_ext, &a_ext[..n1 - 2], &q_clean2, c_low, one);
        b.free(q_clean2[0]);
        b.free(q_clean2[1]);
        b.x(one);
        b.free(one);
    }

    // Reverse of step 1.
    sub_nbit_qq(b, &a_ext, &acc_ext);

    unext_reg(b, a_ovf);
    unext_reg(b, acc_ovf);
    let _ = (acc_ext, a_ext);
}

/// Fast mod_neg using measurement-based Cuccaro for the addition.
pub(crate) fn mod_neg_inplace_fast(b: &mut B, v: &[QubitId], p: U256) {
    for &q in v {
        b.x(q);
    }
    let n = v.len();
    let ca = load_const(b, n, p.wrapping_add(U256::from(1)));
    add_nbit_qq_fast(b, &ca, v);
    unload_const(b, &ca, p.wrapping_add(U256::from(1)));
}


pub(crate) fn mod_add_qb(b: &mut B, acc: &[QubitId], bits: &[BitId], p: U256) {
    // acc := (acc + bits) mod p. `bits` is a classical bit register.
    let a = load_bits(b, bits);
    mod_add_qq_fast(b, acc, &a, p);
    unload_bits(b, &a, bits);
}

pub(crate) fn mod_add_double_qb(b: &mut B, acc: &[QubitId], bits: &[BitId], p: U256) {
    // acc := acc + 2*bits mod p. Reuse a single loaded copy of the classical
    // point and walk it through the cheap secp256k1 double/halve pair.
    let a = load_bits(b, bits);
    mod_double_inplace_fast(b, &a, p);
    mod_add_qq_fast(b, acc, &a, p);
    mod_halve_inplace_fast(b, &a, p);
    unload_bits(b, &a, bits);
}

pub(crate) fn mod_sub_qb(b: &mut B, acc: &[QubitId], bits: &[BitId], p: U256) {
    // acc -= bits mod p. Uses fast mod_sub_qq via neg+add+neg.
    let a = load_bits(b, bits);
    mod_sub_qq_fast(b, acc, &a, p);
    unload_bits(b, &a, bits);
}

// ═══════════════════════════════════════════════════════════════════════════
//  Non-modular n-bit primitives
// ═══════════════════════════════════════════════════════════════════════════

/// Fast Cuccaro sub: `acc -= a mod 2^n` with measurement UMA (0 Toffoli
/// for UMA sweep). Exact gate-level inverse of `cuccaro_add_fast`.
/// Fast `acc += a mod 2^n` using measurement-based Cuccaro.

pub(crate) fn mod_double_inplace_fast(b: &mut B, v: &[QubitId], p: U256) {
    mod_double_inplace_fast_with_dirty(b, v, p, None)
}

pub(crate) fn mod_double_inplace_fast_with_dirty(
    b: &mut B,
    v: &[QubitId],
    p: U256,
    dirty_src: Option<&[QubitId]>,
) {
    let n = v.len();
    let ovf = b.alloc_qubit();
    b.swap(v[n - 1], ovf);
    for i in (0..n - 1).rev() {
        b.swap(v[i], v[i + 1]);
    }
    debug_assert_eq!(n, 256);
    // For secp256k1, p = 2^n - c. After the shift, the old top bit is in
    // `ovf` and the low register holds T mod 2^n for T = 2*v. If ovf=1 then
    // T = 2^n + low and T mod p = low + c; otherwise T mod p = low.
    let c = U256::MAX.wrapping_sub(p).wrapping_add(U256::from(1));
    let use_venting = std::env::var("KAL_VENT_DOUBLE").ok().as_deref() == Some("1")
        && dirty_src.map_or(false, |d| d.len() >= n - 2);
    if let Some(w) = double_carry_trunc_window() {
        // Carry-tail-truncated sparse-constant add (default OFF).
        cadd_nbit_const_direct_trunc_fast(b, v, c, ovf, w);
    } else if use_venting {
        let dirty = dirty_src.unwrap();
        let q_clean2: [QubitId; 2] = [b.alloc_qubit(), b.alloc_qubit()];
        venting::ciadd_dirty_2clean_classical(
            b,
            v,
            &dirty[..n - 2],
            &q_clean2,
            c.as_limbs()[0],
            ovf,
            false,
        );
        b.free(q_clean2[0]);
        b.free(q_clean2[1]);
    } else if direct_const_walks_enabled()
        || std::env::var("KAL_DIRECT_CONST_DOUBLE").ok().as_deref() == Some("1")
    {
        cadd_nbit_const_direct_fast(b, v, c, ovf);
    } else {
        cadd_nbit_const_fast(b, v, c, ovf);
    }
    // Result parity equals the old top bit: even if ovf=0, odd if ovf=1.
    b.cx(v[0], ovf);
    b.free(ovf);
}

pub(crate) fn mod_double_inplace_direct_const_fast(b: &mut B, v: &[QubitId], p: U256) {
    let n = v.len();
    let ovf = b.alloc_qubit();
    b.swap(v[n - 1], ovf);
    for i in (0..n - 1).rev() {
        b.swap(v[i], v[i + 1]);
    }
    debug_assert_eq!(n, 256);
    let c = U256::MAX.wrapping_sub(p).wrapping_add(U256::from(1));
    cadd_nbit_const_direct_fast(b, v, c, ovf);
    b.cx(v[0], ovf);
    b.free(ovf);
}

/// Shift v left by k bits mod p. Returns (spill, flag_inv, ovf) which MUST
/// be passed to mod_shift_right_by_k for cleanup. Bennett-pattern: flags
/// stay alive across the body so the inverse can cleanly cancel them.
///
/// k must be small enough that spill·c < p. For k≤22 with secp256k1 this holds.
pub(crate) fn lowq_shift22() -> bool {
    if d1_phase_corrected_product_core_active() {
        return true;
    }
    // Default OFF: on the current scaffold it no longer reduces every global
    // peak, but it is the measured phase-corrected low-Q shift core for D1.
    // Keep the historical standalone knob for qubit-first experiments.
    match std::env::var("LOWQ_SHIFT22") {
        Ok(v) => v != "0",
        Err(_) => false,
    }
}

pub(crate) fn mod_shift_left_by_k(
    b: &mut B,
    v: &[QubitId],
    p: U256,
    k: usize,
) -> (Vec<QubitId>, QubitId, QubitId) {
    let n = v.len();
    debug_assert_eq!(n, 256);
    let c = U256::MAX.wrapping_sub(p).wrapping_add(U256::from(1));

    let spill = b.alloc_qubits(k);
    let ovf = b.alloc_qubit();
    let flag_inv = b.alloc_qubit();

    // Step 1: k rounds of shift-by-1, capturing top bits into spill.
    for shift_i in 0..k {
        b.swap(v[n - 1], spill[k - 1 - shift_i]);
        for i in (0..n - 1).rev() {
            b.swap(v[i], v[i + 1]);
        }
    }

    // Step 2: add spill · c to v_ext (using ovf as bit n).
    // c = 2^32 + 977 = 2^32 + 2^10 - 2^6 + 2^4 + 2^0.
    // Consolidate 4 bits (6,7,8,9) of 977 into 2^10 - 2^6: saves 2 Cuccaros per shift.
    // Op list: ADD at 0, 4, 10, 32; SUB at 6. Total 5 ops instead of 7.
    let mut v_ext = v.to_vec();
    v_ext.push(ovf);
    let cuccaro_op = |b: &mut B, pos: usize, is_sub: bool| {
        let pad_width = n + 1 - pos;
        let padded = b.alloc_qubits(pad_width);
        for i in 0..k.min(pad_width) {
            b.cx(spill[i], padded[i]);
        }
        let v_slice: Vec<QubitId> = v_ext[pos..n + 1].to_vec();
        let c_in = b.alloc_qubit();
        if lowq_shift22() {
            if is_sub {
                cuccaro_sub(b, &padded, &v_slice, c_in);
            } else {
                cuccaro_add(b, &padded, &v_slice, c_in);
            }
        } else if is_sub {
            // Fast cuccaro: saves ~n CCX per op. Peak during this op (~514
            // transient) is still below the mod_add_qq_fast peak (517) inside
            // the enclosing Solinas, so no global peak increase.
            cuccaro_sub_fast(b, &padded, &v_slice, c_in);
        } else {
            cuccaro_add_fast(b, &padded, &v_slice, c_in);
        }
        b.free(c_in);
        for i in 0..k.min(pad_width) {
            b.cx(spill[i], padded[i]);
        }
        b.free_vec(&padded);
    };
    b.set_phase("shift22_cuccaro_op_0");
    cuccaro_op(b, 0, false);
    b.set_phase("shift22_cuccaro_op_4");
    cuccaro_op(b, 4, false);
    b.set_phase("shift22_cuccaro_op_6");
    cuccaro_op(b, 6, true);
    b.set_phase("shift22_cuccaro_op_10");
    cuccaro_op(b, 10, false);
    b.set_phase("shift22_cuccaro_op_32");
    cuccaro_op(b, 32, false);

    // Step 3: const add.
    b.set_phase("shift22_step3");
    if lowq_shift22() {
        add_nbit_const(b, &v_ext, c);
    } else {
        add_nbit_const_fast(b, &v_ext, c);
    }
    b.x(ovf);
    b.cx(ovf, flag_inv); // flag_inv = NOT(top_bit_after_add) = (value < p)
    b.x(ovf);

    // Step 4: conditional const sub.
    b.set_phase("shift22_step4");
    if lowq_shift22() {
        csub_nbit_const(b, &v_ext, c, flag_inv);
    } else {
        csub_nbit_const_fast(b, &v_ext, c, flag_inv);
    }
    b.x(flag_inv);
    b.cx(flag_inv, ovf);
    b.x(flag_inv);

    (spill, flag_inv, ovf)
}

/// Gate-level inverse of mod_shift_left_by_k.
pub(crate) fn mod_shift_right_by_k(
    b: &mut B,
    v: &[QubitId],
    p: U256,
    k: usize,
    spill: Vec<QubitId>,
    flag_inv: QubitId,
    ovf: QubitId,
) {
    let n = v.len();
    debug_assert_eq!(n, 256);
    let c = U256::MAX.wrapping_sub(p).wrapping_add(U256::from(1));

    let mut v_ext = v.to_vec();
    v_ext.push(ovf);

    // Reverse step 4.
    b.x(flag_inv);
    b.cx(flag_inv, ovf);
    b.x(flag_inv);
    b.set_phase("rshift22_rev_step4");
    if lowq_shift22() {
        cadd_nbit_const(b, &v_ext, c, flag_inv);
    } else {
        cadd_nbit_const_fast(b, &v_ext, c, flag_inv);
    }

    // Reverse step 3.
    b.x(ovf);
    b.cx(ovf, flag_inv);
    b.x(ovf);
    b.set_phase("rshift22_rev_step3");
    if lowq_shift22() {
        sub_nbit_const(b, &v_ext, c);
    } else {
        sub_nbit_const_fast(b, &v_ext, c);
    }
    b.free(flag_inv);
    b.set_phase("rshift22_rev_step2");

    // Reverse step 2: inverse of the consolidated op list (5 ops, in reverse order, flipped signs).
    let cuccaro_op = |b: &mut B, pos: usize, is_sub: bool| {
        let pad_width = n + 1 - pos;
        let padded = b.alloc_qubits(pad_width);
        for i in 0..k.min(pad_width) {
            b.cx(spill[i], padded[i]);
        }
        let v_slice: Vec<QubitId> = v_ext[pos..n + 1].to_vec();
        let c_in = b.alloc_qubit();
        if lowq_shift22() {
            if is_sub {
                cuccaro_sub(b, &padded, &v_slice, c_in);
            } else {
                cuccaro_add(b, &padded, &v_slice, c_in);
            }
        } else if is_sub {
            cuccaro_sub_fast(b, &padded, &v_slice, c_in);
        } else {
            cuccaro_add_fast(b, &padded, &v_slice, c_in);
        }
        b.free(c_in);
        for i in 0..k.min(pad_width) {
            b.cx(spill[i], padded[i]);
        }
        b.free_vec(&padded);
    };
    // Reverse: undo ADD at 32, 10; undo SUB at 6; undo ADD at 4, 0.
    cuccaro_op(b, 32, true); // undo +spill·2^32
    cuccaro_op(b, 10, true); // undo +spill·2^10
    cuccaro_op(b, 6, false); // undo -spill·2^6
    cuccaro_op(b, 4, true); // undo +spill·2^4
    cuccaro_op(b, 0, true); // undo +spill·2^0

    // Reverse step 1: reverse swap cascades.
    for shift_i in (0..k).rev() {
        for i in 0..n - 1 {
            b.swap(v[i], v[i + 1]);
        }
        b.swap(v[n - 1], spill[k - 1 - shift_i]);
    }

    b.free(ovf);
    b.free_vec(&spill);
}

pub(crate) fn mod_shift_left_by_k_lowq(
    b: &mut B,
    v: &[QubitId],
    p: U256,
    k: usize,
) -> (Vec<QubitId>, QubitId, QubitId) {
    let n = v.len();
    debug_assert_eq!(n, 256);
    let c = U256::MAX.wrapping_sub(p).wrapping_add(U256::from(1));

    let spill = b.alloc_qubits(k);
    let ovf = b.alloc_qubit();
    let flag_inv = b.alloc_qubit();

    for shift_i in 0..k {
        b.swap(v[n - 1], spill[k - 1 - shift_i]);
        for i in (0..n - 1).rev() {
            b.swap(v[i], v[i + 1]);
        }
    }

    let mut v_ext = v.to_vec();
    v_ext.push(ovf);
    let cuccaro_op = |b: &mut B, pos: usize, is_sub: bool| {
        let pad_width = n + 1 - pos;
        let padded = b.alloc_qubits(pad_width);
        for i in 0..k.min(pad_width) {
            b.cx(spill[i], padded[i]);
        }
        let v_slice: Vec<QubitId> = v_ext[pos..n + 1].to_vec();
        let c_in = b.alloc_qubit();
        if is_sub {
            cuccaro_sub(b, &padded, &v_slice, c_in);
        } else {
            cuccaro_add(b, &padded, &v_slice, c_in);
        }
        b.free(c_in);
        for i in 0..k.min(pad_width) {
            b.cx(spill[i], padded[i]);
        }
        b.free_vec(&padded);
    };
    cuccaro_op(b, 0, false);
    cuccaro_op(b, 4, false);
    cuccaro_op(b, 6, true);
    cuccaro_op(b, 10, false);
    cuccaro_op(b, 32, false);

    add_nbit_const(b, &v_ext, c);
    b.x(ovf);
    b.cx(ovf, flag_inv);
    b.x(ovf);
    csub_nbit_const(b, &v_ext, c, flag_inv);
    b.x(flag_inv);
    b.cx(flag_inv, ovf);
    b.x(flag_inv);

    (spill, flag_inv, ovf)
}

pub(crate) fn mod_shift_right_by_k_lowq(
    b: &mut B,
    v: &[QubitId],
    p: U256,
    k: usize,
    spill: Vec<QubitId>,
    flag_inv: QubitId,
    ovf: QubitId,
) {
    let n = v.len();
    debug_assert_eq!(n, 256);
    let c = U256::MAX.wrapping_sub(p).wrapping_add(U256::from(1));

    let mut v_ext = v.to_vec();
    v_ext.push(ovf);

    b.x(flag_inv);
    b.cx(flag_inv, ovf);
    b.x(flag_inv);
    cadd_nbit_const(b, &v_ext, c, flag_inv);

    b.x(ovf);
    b.cx(ovf, flag_inv);
    b.x(ovf);
    sub_nbit_const(b, &v_ext, c);
    b.free(flag_inv);

    let cuccaro_op = |b: &mut B, pos: usize, is_sub: bool| {
        let pad_width = n + 1 - pos;
        let padded = b.alloc_qubits(pad_width);
        for i in 0..k.min(pad_width) {
            b.cx(spill[i], padded[i]);
        }
        let v_slice: Vec<QubitId> = v_ext[pos..n + 1].to_vec();
        let c_in = b.alloc_qubit();
        if is_sub {
            cuccaro_sub(b, &padded, &v_slice, c_in);
        } else {
            cuccaro_add(b, &padded, &v_slice, c_in);
        }
        b.free(c_in);
        for i in 0..k.min(pad_width) {
            b.cx(spill[i], padded[i]);
        }
        b.free_vec(&padded);
    };
    cuccaro_op(b, 32, true);
    cuccaro_op(b, 10, true);
    cuccaro_op(b, 6, false);
    cuccaro_op(b, 4, true);
    cuccaro_op(b, 0, true);

    for shift_i in (0..k).rev() {
        for i in 0..n - 1 {
            b.swap(v[i], v[i + 1]);
        }
        b.swap(v[n - 1], spill[k - 1 - shift_i]);
    }

    b.free(ovf);
    b.free_vec(&spill);
}

/// Fast `v := v/2 mod p`. Explicit reverse of `mod_double_inplace` with
/// measurement-based Cuccaro (not emit_inverse).
pub(crate) fn mod_halve_inplace_fast(b: &mut B, v: &[QubitId], p: U256) {
    mod_halve_inplace_fast_with_dirty(b, v, p, None)
}

pub(crate) fn mod_halve_inplace_direct_const_fast(b: &mut B, v: &[QubitId], p: U256) {
    let n = v.len();
    let ovf = b.alloc_qubit();
    debug_assert_eq!(n, 256);
    let c = U256::MAX.wrapping_sub(p).wrapping_add(U256::from(1));
    b.cx(v[0], ovf);
    csub_nbit_const_direct_fast(b, v, c, ovf);
    for i in 0..n - 1 {
        b.swap(v[i], v[i + 1]);
    }
    b.swap(v[n - 1], ovf);
    b.free(ovf);
}

/// Variant of `mod_halve_inplace_fast` that optionally borrows `dirty_src`
/// qubits for the controlled-sub step, using Gidney's venting
/// `cisub_dirty_2clean_classical`. Saves n transient qubits at the peak
/// when dirty qubits are available from the caller.
pub(crate) fn mod_halve_inplace_fast_with_dirty(
    b: &mut B,
    v: &[QubitId],
    p: U256,
    dirty_src: Option<&[QubitId]>,
) {
    let n = v.len();
    let ovf = b.alloc_qubit();
    debug_assert_eq!(n, 256);
    let c = U256::MAX.wrapping_sub(p).wrapping_add(U256::from(1));
    b.cx(v[0], ovf);
    // If caller provided enough dirty qubits AND c fits in u64 (it does
    // for secp256k1: c = 2^32 + 977), use the venting variant.
    let use_venting = kal_vent_halve_enabled() && dirty_src.map_or(false, |d| d.len() >= n - 2);
    if let Some(w) = double_carry_trunc_window() {
        // Carry-tail-truncated sparse-constant sub (inverse of the truncated
        // double; default OFF; same window so double/halve stay exact inverses).
        csub_nbit_const_direct_trunc_fast(b, v, c, ovf, w);
    } else if use_venting {
        // c as u64 (it fits: c = 0x1000003D1).
        // For n=256, we still need to pass the full 256-bit constant via u64.
        // Since c only has 33 bits, u64 is fine.
        let c_u64: u64 = c.as_limbs()[0] | (c.as_limbs()[1] << 32); // hack for U256
                                                                    // Actually, U256 limbs are u64[4]. Bit 32 of U256 is limbs[0] bit 32.
                                                                    // limbs[0] holds bits 0..64. So just take limbs[0] for bits < 64.
        let c_low = c.as_limbs()[0];
        let dirty = dirty_src.unwrap();
        let dirty_slice = &dirty[..n - 2];
        // We need 2 clean ancilla. Alloc them fresh.
        let q_clean2: [QubitId; 2] = [b.alloc_qubit(), b.alloc_qubit()];
        venting::cisub_dirty_2clean_classical(b, v, dirty_slice, &q_clean2, c_low, ovf);
        b.free(q_clean2[0]);
        b.free(q_clean2[1]);
        let _ = c_u64; // unused, c_low is the right value
    } else if direct_const_walks_enabled()
        || std::env::var("KAL_DIRECT_CONST_HALVE").ok().as_deref() == Some("1")
    {
        csub_nbit_const_direct_fast(b, v, c, ovf);
    } else {
        csub_nbit_const_fast(b, v, c, ovf);
    }
    for i in 0..n - 1 {
        b.swap(v[i], v[i + 1]);
    }
    b.swap(v[n - 1], ovf);
    b.free(ovf);
}

/// Controlled lazy mod-double: EXACT controlled form of `mod_double_inplace_fast`
/// (Solinas reduction, lazy [0,2^n) coset rep, same carry-trunc window). Identity
/// when ctrl=0. Used by the K=2 prototype's conditional 2nd double so it composes
/// correctly with the uncontrolled `mod_double_inplace_fast` in the apply.
pub(crate) fn cmod_double_inplace_lazy(b: &mut B, v: &[QubitId], p: U256, ctrl: QubitId) {
    let n = v.len();
    let ovf = b.alloc_qubit();
    cswap(b, ctrl, v[n - 1], ovf);
    for i in (0..n - 1).rev() {
        cswap(b, ctrl, v[i], v[i + 1]);
    }
    let c = U256::MAX.wrapping_sub(p).wrapping_add(U256::from(1));
    if let Some(w) = double_carry_trunc_window() {
        cadd_nbit_const_direct_trunc_fast(b, v, c, ovf, w);
    } else if direct_const_walks_enabled()
        || std::env::var("KAL_DIRECT_CONST_DOUBLE").ok().as_deref() == Some("1")
    {
        cadd_nbit_const_direct_fast(b, v, c, ovf);
    } else {
        cadd_nbit_const_fast(b, v, c, ovf);
    }
    // Clear ovf: result parity == old top bit == ovf (gated by ctrl).
    b.ccx(ctrl, v[0], ovf);
    b.free(ovf);
}

/// Controlled lazy mod-halve: EXACT controlled form of `mod_halve_inplace_fast`
/// (inverse of `cmod_double_inplace_lazy`, same window). Identity when ctrl=0.
pub(crate) fn cmod_halve_inplace_lazy(b: &mut B, v: &[QubitId], p: U256, ctrl: QubitId) {
    let n = v.len();
    let ovf = b.alloc_qubit();
    let c = U256::MAX.wrapping_sub(p).wrapping_add(U256::from(1));
    b.ccx(ctrl, v[0], ovf);
    if let Some(w) = double_carry_trunc_window() {
        csub_nbit_const_direct_trunc_fast(b, v, c, ovf, w);
    } else if direct_const_walks_enabled()
        || std::env::var("KAL_DIRECT_CONST_HALVE").ok().as_deref() == Some("1")
    {
        csub_nbit_const_direct_fast(b, v, c, ovf);
    } else {
        csub_nbit_const_fast(b, v, c, ovf);
    }
    for i in 0..n - 1 {
        cswap(b, ctrl, v[i], v[i + 1]);
    }
    cswap(b, ctrl, v[n - 1], ovf);
    b.free(ovf);
}

// ═══════════════════════════════════════════════════════════════════════════
//  Conditional modular add/sub helpers
// ═══════════════════════════════════════════════════════════════════════════
//
// Used by the multipliers. Each variant loads `(ctrl ? a : 0)` into a
// fresh temporary via CCX or CX_if, runs the unconditional mod_add_qq /
// mod_sub_qq, then unloads.

/// Like `cmp_lt_into` but uses carry-ancilla + measurement-based uncompute
/// for the inv_MAJ sweep. Saves n CCX. NOT emit_inverse-safe.

/// Like `mod_add_qq` but uses `cmp_lt_into_fast` for the flag uncompute.
/// NOT safe inside emit_inverse blocks.
pub(crate) fn mod_add_qq_fast(b: &mut B, acc: &[QubitId], a: &[QubitId], p: U256) {
    let n = acc.len();
    assert_eq!(n, a.len());
    debug_assert_eq!(n, 256);

    let (acc_ext, acc_ovf) = ext_reg(b, acc);
    let (a_ext, a_ovf) = ext_reg(b, a);

    // Use fast (measurement-based) Cuccaro everywhere.
    add_nbit_qq_fast(b, &a_ext, &acc_ext);
    let c = U256::MAX.wrapping_sub(p).wrapping_add(U256::from(1));
    // add_nbit_const with fast Cuccaro OR venting (using `a` as dirty).
    let use_vent = kal_vent_modadd_enabled();
    if use_vent {
        let n1 = acc_ext.len();
        // Use `a_ext` as dirty qubits (it was just used as add operand,
        // its value is preserved through the venting sub-protocol).
        let c_low = c.as_limbs()[0];
        let q_clean2: [QubitId; 2] = [b.alloc_qubit(), b.alloc_qubit()];
        venting::iadd_dirty_2clean_classical(
            b,
            &acc_ext,
            &a_ext[..n1 - 2],
            &q_clean2,
            c_low,
            false,
        );
        b.free(q_clean2[0]);
        b.free(q_clean2[1]);
    } else if secp_direct_const_arith_enabled() {
        add_nbit_const_direct_uncontrolled_fast(b, &acc_ext, c);
    } else {
        let n1 = acc_ext.len();
        let ca = load_const(b, n1, c);
        add_nbit_qq_fast(b, &ca, &acc_ext);
        unload_const(b, &ca, c);
    }
    let flag = b.alloc_qubit();
    b.cx(acc_ovf, flag);
    b.x(flag);
    // csub_nbit_const with fast Cuccaro OR venting.
    if use_vent {
        let c_low = c.as_limbs()[0];
        let n1 = acc_ext.len();
        let q_clean2: [QubitId; 2] = [b.alloc_qubit(), b.alloc_qubit()];
        venting::cisub_dirty_2clean_classical(
            b,
            &acc_ext,
            &a_ext[..n1 - 2],
            &q_clean2,
            c_low,
            flag,
        );
        b.free(q_clean2[0]);
        b.free(q_clean2[1]);
    } else if secp_direct_const_arith_enabled() {
        csub_nbit_const_direct_fast(b, &acc_ext, c, flag);
    } else {
        let n1 = acc_ext.len();
        let ca = b.alloc_qubits(n1);
        for i in 0..n1 {
            if bit(c, i) {
                b.cx(flag, ca[i]);
            }
        }
        sub_nbit_qq_fast(b, &ca, &acc_ext);
        for i in 0..n1 {
            if bit(c, i) {
                b.cx(flag, ca[i]);
            }
        }
        b.free_vec(&ca);
    }
    b.x(flag);
    b.cx(flag, acc_ovf);
    cmp_lt_into_fast(b, &acc_ext[..n], &a_ext[..n], flag);
    b.free(flag);

    unext_reg(b, a_ovf);
    unext_reg(b, acc_ovf);
    let _ = (acc_ext, a_ext);
}

/// Specialization of mod_add_qq_fast when acc = 0 on entry. Replaces the
/// initial Cuccaro add with CX-copy (0 CCX instead of n-1 CCX).
/// Saves 255 CCX per call.
pub(crate) fn mod_add_qq_fast_from_zero(b: &mut B, acc: &[QubitId], a: &[QubitId], p: U256) {
    let n = acc.len();
    assert_eq!(n, a.len());
    debug_assert_eq!(n, 256);

    let (acc_ext, acc_ovf) = ext_reg(b, acc);
    let (a_ext, a_ovf) = ext_reg(b, a);

    // acc is 0 on entry. CX-copy a into acc (0 CCX). Top bits both 0.
    for i in 0..n {
        b.cx(a[i], acc[i]);
    }
    // acc_ovf and a_ovf are both 0 (both freshly allocated as 0 by ext_reg).

    let c = U256::MAX.wrapping_sub(p).wrapping_add(U256::from(1));
    let use_vent = kal_vent_modadd_enabled();
    if use_vent {
        let n1 = acc_ext.len();
        let c_low = c.as_limbs()[0];
        let q_clean2: [QubitId; 2] = [b.alloc_qubit(), b.alloc_qubit()];
        venting::iadd_dirty_2clean_classical(
            b,
            &acc_ext,
            &a_ext[..n1 - 2],
            &q_clean2,
            c_low,
            false,
        );
        b.free(q_clean2[0]);
        b.free(q_clean2[1]);
    } else {
        let n1 = acc_ext.len();
        let ca = load_const(b, n1, c);
        add_nbit_qq_fast(b, &ca, &acc_ext);
        unload_const(b, &ca, c);
    }
    let flag = b.alloc_qubit();
    b.cx(acc_ovf, flag);
    b.x(flag);
    if use_vent {
        let c_low = c.as_limbs()[0];
        let n1 = acc_ext.len();
        let q_clean2: [QubitId; 2] = [b.alloc_qubit(), b.alloc_qubit()];
        venting::cisub_dirty_2clean_classical(
            b,
            &acc_ext,
            &a_ext[..n1 - 2],
            &q_clean2,
            c_low,
            flag,
        );
        b.free(q_clean2[0]);
        b.free(q_clean2[1]);
    } else {
        let n1 = acc_ext.len();
        let ca = b.alloc_qubits(n1);
        for i in 0..n1 {
            if bit(c, i) {
                b.cx(flag, ca[i]);
            }
        }
        sub_nbit_qq_fast(b, &ca, &acc_ext);
        for i in 0..n1 {
            if bit(c, i) {
                b.cx(flag, ca[i]);
            }
        }
        b.free_vec(&ca);
    }
    b.x(flag);
    b.cx(flag, acc_ovf);
    cmp_lt_into_fast(b, &acc_ext[..n], &a_ext[..n], flag);
    b.free(flag);

    unext_reg(b, a_ovf);
    unext_reg(b, acc_ovf);
    let _ = (acc_ext, a_ext);
}

pub(crate) fn cmod_add_qq(b: &mut B, acc: &[QubitId], a: &[QubitId], ctrl: QubitId, p: U256) {
    let n = acc.len();
    let f = b.alloc_qubits(n);
    for i in 0..n {
        b.ccx(ctrl, a[i], f[i]);
    }
    mod_add_qq_fast(b, acc, &f, p);
    // Gidney measurement-based AND uncomputation: f[i] = ctrl AND a[i],
    // which is unchanged by mod_add_qq (Cuccaro restores the addend).
    // HMR + classically-conditioned CZ costs 0 Toffoli vs 256 CCX.
    for i in 0..n {
        let m = b.alloc_bit();
        b.hmr(f[i], m);
        b.cz_if(ctrl, a[i], m);
    }
    b.free_vec(&f);
}

pub(crate) fn cmod_sub_qq(b: &mut B, acc: &[QubitId], a: &[QubitId], ctrl: QubitId, p: U256) {
    let n = acc.len();
    let f = b.alloc_qubits(n);
    for i in 0..n {
        b.ccx(ctrl, a[i], f[i]);
    }
    mod_sub_qq_fast(b, acc, &f, p);
    for i in 0..n {
        let m = b.alloc_bit();
        b.hmr(f[i], m);
        b.cz_if(ctrl, a[i], m);
    }
    b.free_vec(&f);
}

pub(crate) fn cmod_add_qq_lowq(b: &mut B, acc: &[QubitId], a: &[QubitId], ctrl: QubitId, p: U256) {
    let n = acc.len();
    let f = b.alloc_qubits(n);
    for i in 0..n {
        b.ccx(ctrl, a[i], f[i]);
    }
    mod_add_qq(b, acc, &f, p);
    for i in 0..n {
        let m = b.alloc_bit();
        b.hmr(f[i], m);
        b.cz_if(ctrl, a[i], m);
    }
    b.free_vec(&f);
}

pub(crate) fn cmod_sub_qq_lowq(b: &mut B, acc: &[QubitId], a: &[QubitId], ctrl: QubitId, p: U256) {
    let n = acc.len();
    let f = b.alloc_qubits(n);
    for i in 0..n {
        b.ccx(ctrl, a[i], f[i]);
    }
    mod_sub_qq(b, acc, &f, p);
    for i in 0..n {
        let m = b.alloc_bit();
        b.hmr(f[i], m);
        b.cz_if(ctrl, a[i], m);
    }
    b.free_vec(&f);
}

