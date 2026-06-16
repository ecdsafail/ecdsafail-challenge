use alloy_primitives::U256;

use crate::circuit::{BitId, QubitId};

use super::super::{
    bit, cadd_nbit_const_direct_fast, cmp_lt_into, cswap, csub_nbit_const_direct_fast, B,
    SECP256K1_P,
};
use super::{primitives, rfold, schedule};

fn p_plus_1() -> U256 {
    SECP256K1_P.wrapping_add(U256::from(1u64))
}

fn half_p_ceil() -> U256 {
    (SECP256K1_P >> 1) + U256::from(1u64)
}

fn xor_const(b: &mut B, reg: &[QubitId], value: usize) {
    for (i, &q) in reg.iter().enumerate() {
        if (value >> i) & 1 == 1 {
            b.x(q);
        }
    }
}

fn and_into_clean(b: &mut B, ctrls: &[QubitId], target: QubitId) {
    match ctrls {
        [] => b.x(target),
        [a] => b.cx(*a, target),
        [a, c] => b.ccx(*a, *c, target),
        _ => {
            let mid = ctrls.len() / 2;
            let left = b.alloc_qubit();
            let right = b.alloc_qubit();
            and_into_clean(b, &ctrls[..mid], left);
            and_into_clean(b, &ctrls[mid..], right);
            b.ccx(left, right, target);
            and_into_clean(b, &ctrls[mid..], right);
            b.release_zeroed(right);
            and_into_clean(b, &ctrls[..mid], left);
            b.release_zeroed(left);
        }
    }
}

fn mcx_linear_into(b: &mut B, ctrls: &[QubitId], target: QubitId) {
    match ctrls {
        [] => b.x(target),
        [a] => b.cx(*a, target),
        [a, c] => b.ccx(*a, *c, target),
        _ => {
            let work = b.alloc_qubits(ctrls.len() - 2);
            b.ccx(ctrls[0], ctrls[1], work[0]);
            for i in 2..ctrls.len() - 1 {
                b.ccx(work[i - 2], ctrls[i], work[i - 1]);
            }
            b.ccx(work[work.len() - 1], ctrls[ctrls.len() - 1], target);
            for i in (2..ctrls.len() - 1).rev() {
                b.ccx(work[i - 2], ctrls[i], work[i - 1]);
            }
            b.ccx(ctrls[0], ctrls[1], work[0]);
            b.free_vec(&work);
        }
    }
}

fn mcx_into(b: &mut B, ctrls: &[QubitId], target: QubitId) {
    mcx_linear_into(b, ctrls, target);
}

fn eq_const_into(b: &mut B, reg: &[QubitId], value: usize, target: QubitId) {
    for (i, &q) in reg.iter().enumerate() {
        if (value >> i) & 1 == 0 {
            b.x(q);
        }
    }
    mcx_into(b, reg, target);
    for (i, &q) in reg.iter().enumerate() {
        if (value >> i) & 1 == 0 {
            b.x(q);
        }
    }
}

fn or_is_zero(b: &mut B, reg: &[QubitId], out: QubitId) {
    for &q in reg {
        b.x(q);
    }
    and_into_clean(b, reg, out);
    for &q in reg {
        b.x(q);
    }
}

fn or_nonzero(b: &mut B, reg: &[QubitId], out: QubitId) {
    or_is_zero(b, reg, out);
    b.x(out);
}

fn controlled_field_neg(b: &mut B, ctrl: QubitId, reg: &[QubitId]) {
    for &q in reg {
        b.cx(ctrl, q);
    }
    cadd_nbit_const_direct_fast(b, reg, p_plus_1(), ctrl);
}

fn field_neg(b: &mut B, reg: &[QubitId]) {
    let one = b.alloc_qubit();
    b.x(one);
    controlled_field_neg(b, one, reg);
    b.x(one);
    b.release_zeroed(one);
}

fn load_creg_temp_257(b: &mut B, bits: &[BitId]) -> Vec<QubitId> {
    assert_eq!(bits.len(), 256);
    let temp = b.alloc_qubits(257);
    for (i, &bit) in bits.iter().enumerate() {
        b.x_if(temp[i], bit);
    }
    temp
}

fn unload_creg_temp_257(b: &mut B, temp: &[QubitId], bits: &[BitId]) {
    assert_eq!(bits.len(), 256);
    assert_eq!(temp.len(), 257);
    for (i, &bit) in bits.iter().enumerate() {
        b.x_if(temp[i], bit);
    }
}

fn mod_add_creg_qload(b: &mut B, acc: &[QubitId], bits: &[BitId]) {
    b.set_phase("shrunken_pz_creg_add_qload");
    assert_eq!(acc.len(), 257);
    let temp = load_creg_temp_257(b, bits);
    let one = b.alloc_qubit();
    b.x(one);
    rfold::controlled_mod_add_rfold_mbu(b, one, acc, &temp);
    b.x(one);
    b.release_zeroed(one);
    unload_creg_temp_257(b, &temp, bits);
    b.free_vec(&temp);
}

fn mod_sub_creg_qload(b: &mut B, acc: &[QubitId], bits: &[BitId]) {
    b.set_phase("shrunken_pz_creg_sub_qload");
    assert_eq!(acc.len(), 257);
    let temp = load_creg_temp_257(b, bits);
    let one = b.alloc_qubit();
    b.x(one);
    rfold::controlled_mod_sub_rfold_mbu(b, one, acc, &temp);
    b.x(one);
    b.release_zeroed(one);
    unload_creg_temp_257(b, &temp, bits);
    b.free_vec(&temp);
}

fn mod_sub_from_creg_qload(b: &mut B, acc: &[QubitId], bits: &[BitId]) {
    field_neg(b, acc);
    mod_add_creg_qload(b, acc, bits);
}

fn r_const() -> U256 {
    U256::from(977u64) + (U256::from(1u64) << 32)
}

fn xor_or_into(b: &mut B, p: QubitId, q: QubitId, target: QubitId) {
    b.cx(p, target);
    b.cx(q, target);
    b.ccx(p, q, target);
}

fn compare_geq_p_secp256k1_into(b: &mut B, a: &[QubitId], flag: QubitId) {
    assert_eq!(a.len(), 257);
    b.set_phase("shrunken_pz_cmp_p");

    let low4_all_ones = b.alloc_qubit();
    and_into_clean(b, &a[..4], low4_all_ones);

    let low_tail_or = b.alloc_qubit();
    xor_or_into(b, a[4], low4_all_ones, low_tail_or);

    let low6_ge = b.alloc_qubit();
    b.ccx(a[5], low_tail_or, low6_ge);

    let hi_or_67 = b.alloc_qubit();
    xor_or_into(b, a[6], a[7], hi_or_67);

    let hi_or_89 = b.alloc_qubit();
    xor_or_into(b, a[8], a[9], hi_or_89);

    let hi4_nonzero = b.alloc_qubit();
    xor_or_into(b, hi_or_67, hi_or_89, hi4_nonzero);

    let low10_ge = b.alloc_qubit();
    xor_or_into(b, low6_ge, hi4_nonzero, low10_ge);

    let mid_all_ones = b.alloc_qubit();
    and_into_clean(b, &a[10..32], mid_all_ones);

    let mid_and_low = b.alloc_qubit();
    b.ccx(mid_all_ones, low10_ge, mid_and_low);

    let tail_or = b.alloc_qubit();
    xor_or_into(b, a[32], mid_and_low, tail_or);

    let high_all_ones = b.alloc_qubit();
    and_into_clean(b, &a[33..256], high_all_ones);

    let high_and_tail = b.alloc_qubit();
    b.ccx(high_all_ones, tail_or, high_and_tail);

    xor_or_into(b, a[256], high_and_tail, flag);

    b.ccx(high_all_ones, tail_or, high_and_tail);
    b.release_zeroed(high_and_tail);

    and_into_clean(b, &a[33..256], high_all_ones);
    b.release_zeroed(high_all_ones);

    xor_or_into(b, a[32], mid_and_low, tail_or);
    b.release_zeroed(tail_or);

    b.ccx(mid_all_ones, low10_ge, mid_and_low);
    b.release_zeroed(mid_and_low);

    and_into_clean(b, &a[10..32], mid_all_ones);
    b.release_zeroed(mid_all_ones);

    xor_or_into(b, low6_ge, hi4_nonzero, low10_ge);
    b.release_zeroed(low10_ge);

    xor_or_into(b, hi_or_67, hi_or_89, hi4_nonzero);
    b.release_zeroed(hi4_nonzero);

    xor_or_into(b, a[8], a[9], hi_or_89);
    b.release_zeroed(hi_or_89);

    xor_or_into(b, a[6], a[7], hi_or_67);
    b.release_zeroed(hi_or_67);

    b.ccx(a[5], low_tail_or, low6_ge);
    b.release_zeroed(low6_ge);

    xor_or_into(b, a[4], low4_all_ones, low_tail_or);
    b.release_zeroed(low_tail_or);

    and_into_clean(b, &a[..4], low4_all_ones);
    b.release_zeroed(low4_all_ones);
}

fn reduce_once_secp256k1_from_rfold_lowq(b: &mut B, acc: &[QubitId], flag: QubitId) {
    compare_geq_p_secp256k1_into(b, acc, flag);
    cadd_nbit_const_direct_fast(b, &acc[..256], r_const(), flag);
}

fn reduce_once_secp256k1_from_rfold_reverse_lowq(b: &mut B, acc: &[QubitId], flag: QubitId) {
    csub_nbit_const_direct_fast(b, &acc[..256], r_const(), flag);
    compare_geq_p_secp256k1_into(b, acc, flag);
}

fn ctrl_aliases(ctrl: QubitId, left: &[QubitId], right: &[QubitId]) -> bool {
    left.contains(&ctrl) || right.contains(&ctrl)
}

fn controlled_mod_add_rfold_safe(
    b: &mut B,
    ctrl: QubitId,
    acc: &[QubitId],
    addend: &[QubitId],
) {
    if ctrl_aliases(ctrl, acc, addend) {
        let copied = b.alloc_qubit();
        b.cx(ctrl, copied);
        rfold::controlled_mod_add_rfold_mbu(b, copied, acc, addend);
        b.cx(ctrl, copied);
        b.release_zeroed(copied);
    } else {
        rfold::controlled_mod_add_rfold_mbu(b, ctrl, acc, addend);
    }
}

fn controlled_mod_sub_rfold_safe(
    b: &mut B,
    ctrl: QubitId,
    acc: &[QubitId],
    subtrahend: &[QubitId],
) {
    if ctrl_aliases(ctrl, acc, subtrahend) {
        let copied = b.alloc_qubit();
        b.cx(ctrl, copied);
        rfold::controlled_mod_sub_rfold_mbu(b, copied, acc, subtrahend);
        b.cx(ctrl, copied);
        b.release_zeroed(copied);
    } else {
        rfold::controlled_mod_sub_rfold_mbu(b, ctrl, acc, subtrahend);
    }
}

fn horner_forward(b: &mut B, result: &[QubitId], a: &[QubitId], x: &[QubitId]) -> QubitId {
    let n = result.len();
    assert_eq!(n, 257);
    assert_eq!(a.len(), n);
    assert_eq!(x.len(), n);
    b.set_phase("shrunken_pz_horner_add");
    controlled_mod_add_rfold_safe(b, x[n - 1], result, a);
    for i in (0..n - 1).rev() {
        b.set_phase("shrunken_pz_horner_double");
        rfold::mod_double_rfold_mbu(b, result);
        b.set_phase("shrunken_pz_horner_add");
        controlled_mod_add_rfold_safe(b, x[i], result, a);
    }
    b.set_phase("shrunken_pz_horner_reduce");
    let flag = b.alloc_qubit();
    reduce_once_secp256k1_from_rfold_lowq(b, result, flag);
    flag
}

fn horner_reverse(b: &mut B, result: &[QubitId], a: &[QubitId], x: &[QubitId], flag: QubitId) {
    let n = result.len();
    assert_eq!(n, 257);
    assert_eq!(a.len(), n);
    assert_eq!(x.len(), n);
    b.set_phase("shrunken_pz_horner_unreduce");
    reduce_once_secp256k1_from_rfold_reverse_lowq(b, result, flag);
    b.free(flag);
    for i in 0..n - 1 {
        b.set_phase("shrunken_pz_horner_sub");
        controlled_mod_sub_rfold_safe(b, x[i], result, a);
        b.set_phase("shrunken_pz_horner_halve");
        rfold::mod_halve_rfold_mbu(b, result);
    }
    b.set_phase("shrunken_pz_horner_sub");
    controlled_mod_sub_rfold_safe(b, x[n - 1], result, a);
}

fn horner_canonical_flag_consume(b: &mut B, result: &[QubitId], flag: QubitId) {
    b.set_phase("shrunken_pz_horner_flag_consume");
    compare_geq_p_secp256k1_into(b, result, flag);
    b.free(flag);
}

fn mod_mac_inplace(b: &mut B, acc: &[QubitId], a: &[QubitId], x: &[QubitId]) {
    b.set_phase("shrunken_pz_mod_mac_preshift");
    assert_eq!(acc.len(), 257);
    assert_eq!(a.len(), 257);
    assert_eq!(x.len(), 257);
    for _ in 0..256 {
        rfold::mod_halve_rfold_mbu(b, acc);
    }
    let flag = horner_forward(b, acc, a, x);
    horner_canonical_flag_consume(b, acc, flag);
}

fn mod_msc_inplace(b: &mut B, acc: &[QubitId], a: &[QubitId], x: &[QubitId]) {
    b.set_phase("shrunken_pz_mod_msc");
    assert_eq!(acc.len(), 257);
    assert_eq!(a.len(), 257);
    assert_eq!(x.len(), 257);
    let pre_flag = b.alloc_qubit();
    horner_reverse(b, acc, a, x, pre_flag);
    b.set_phase("shrunken_pz_mod_msc_postshift");
    for _ in 0..256 {
        rfold::mod_double_rfold_mbu(b, acc);
    }
    let post_flag = b.alloc_qubit();
    reduce_once_secp256k1_from_rfold_lowq(b, acc, post_flag);
    horner_canonical_flag_consume(b, acc, post_flag);
}

fn compare_geq_const(b: &mut B, reg: &[QubitId], c: U256, out: QubitId) {
    let creg = b.alloc_qubits(reg.len());
    for i in 0..reg.len().min(256) {
        if bit(c, i) {
            b.x(creg[i]);
        }
    }
    cmp_lt_into(b, reg, &creg, out);
    b.x(out);
    for i in 0..reg.len().min(256) {
        if bit(c, i) {
            b.x(creg[i]);
        }
    }
    b.free_vec(&creg);
}

fn ctrl_inc(b: &mut B, ctrl: QubitId, reg: &[QubitId]) {
    let one = b.alloc_qubits(reg.len());
    b.x(one[0]);
    primitives::ctrl_add(b, ctrl, reg, &one);
    b.x(one[0]);
    b.free_vec(&one);
}

fn ctrl_dec(b: &mut B, ctrl: QubitId, reg: &[QubitId]) {
    let one = b.alloc_qubits(reg.len());
    b.x(one[0]);
    primitives::ctrl_sub(b, ctrl, reg, &one);
    b.x(one[0]);
    b.free_vec(&one);
}

fn add_refs(b: &mut B, acc: &[QubitId], addend: &[QubitId]) {
    let one = b.alloc_qubit();
    b.x(one);
    primitives::ctrl_add(b, one, acc, addend);
    b.x(one);
    b.release_zeroed(one);
}

fn sub_refs(b: &mut B, acc: &[QubitId], subtrahend: &[QubitId]) {
    let one = b.alloc_qubit();
    b.x(one);
    primitives::ctrl_sub(b, one, acc, subtrahend);
    b.x(one);
    b.release_zeroed(one);
}

fn deposit_bitlen_position(b: &mut B, src: &[QubitId], pos: &[QubitId]) {
    let n = src.len();
    for k in (0..n).rev() {
        let flag = b.alloc_qubit();
        let mut ctrls = Vec::with_capacity(n - k);
        ctrls.push(src[k]);
        for &q in &src[k + 1..] {
            b.x(q);
            ctrls.push(q);
        }
        mcx_into(b, &ctrls, flag);
        let gd = k ^ (k + 1);
        for (j, &q) in pos.iter().enumerate() {
            if (gd >> j) & 1 == 1 {
                b.cx(flag, q);
            }
        }
        mcx_into(b, &ctrls, flag);
        for &q in &src[k + 1..] {
            b.x(q);
        }
        b.release_zeroed(flag);
    }
}

fn bit_length_lean_middle(
    b: &mut B,
    src: &[QubitId],
    pos: &[QubitId],
    body: impl FnOnce(&mut B) -> bool,
) {
    if src.is_empty() {
        body(b);
        return;
    }
    b.set_phase("shrunken_pz_bitlen");
    deposit_bitlen_position(b, src, pos);
    let clean = body(b);
    if clean {
        deposit_bitlen_position(b, src, pos);
    }
}

fn bit_length_lean(b: &mut B, src: &[QubitId], s: &[QubitId], dec: bool) {
    if src.is_empty() {
        return;
    }
    let n = src.len();
    let pos = b.alloc_qubits(s.len());
    xor_const(b, &pos, n);
    bit_length_lean_middle(b, src, &pos, |b| {
        if dec {
            for &q in s {
                b.x(q);
            }
        }
        add_refs(b, s, &pos);
        let one = b.alloc_qubit();
        b.x(one);
        ctrl_inc(b, one, s);
        b.x(one);
        b.release_zeroed(one);
        if dec {
            for &q in s {
                b.x(q);
            }
        }
        true
    });
    xor_const(b, &pos, n);
    b.free_vec(&pos);
}

fn add_const_uncontrolled(b: &mut B, reg: &[QubitId], value: i64) {
    let one = b.alloc_qubit();
    b.x(one);
    let w = reg.len();
    let val = i128::from(value).rem_euclid(1i128 << w) as u128;
    cadd_nbit_const_direct_fast(b, reg, U256::from(val), one);
    b.x(one);
    b.release_zeroed(one);
}

fn clz_diff_body_middle(
    b: &mut B,
    a: &[QubitId],
    c: &[QubitId],
    w: usize,
    lo_a: usize,
    lo_b: usize,
    body: impl FnOnce(&mut B, &[QubitId]),
) {
    let aw = &a[lo_a.min(a.len())..];
    let cw = &c[lo_b.min(c.len())..];
    let pa = b.alloc_qubits(w);
    let pb = b.alloc_qubits(w);
    xor_const(b, &pa, aw.len());
    bit_length_lean_middle(b, aw, &pa, |_| false);
    xor_const(b, &pb, cw.len());
    bit_length_lean_middle(b, cw, &pb, |_| false);

    add_const_uncontrolled(b, &pa, 1 + lo_a as i64);
    sub_refs(b, &pa, &pb);
    add_const_uncontrolled(b, &pa, -(1 + lo_b as i64));

    body(b, &pa);

    add_const_uncontrolled(b, &pa, 1 + lo_b as i64);
    add_refs(b, &pa, &pb);
    add_const_uncontrolled(b, &pa, -(1 + lo_a as i64));

    bit_length_lean_middle(b, cw, &pb, |_| false);
    xor_const(b, &pb, cw.len());
    bit_length_lean_middle(b, aw, &pa, |_| false);
    xor_const(b, &pa, aw.len());
    b.free_vec(&pb);
    b.free_vec(&pa);
}

fn rotate_left_pow2(b: &mut B, reg: &[QubitId], ctrl: QubitId, shift: usize) {
    if shift == 0 || shift >= reg.len() {
        return;
    }
    for i in (0..reg.len() - shift).rev() {
        cswap(b, ctrl, reg[i], reg[i + shift]);
    }
}

fn rotate_right_pow2(b: &mut B, reg: &[QubitId], ctrl: QubitId, shift: usize) {
    if shift == 0 || shift >= reg.len() {
        return;
    }
    for i in 0..reg.len() - shift {
        cswap(b, ctrl, reg[i], reg[i + shift]);
    }
}

fn rotate_left(b: &mut B, reg: &[QubitId], s: &[QubitId]) {
    b.set_phase("shrunken_pz_shift_left");
    for (j, &ctrl) in s.iter().enumerate() {
        rotate_left_pow2(b, reg, ctrl, 1usize << j);
    }
}

fn rotate_right(b: &mut B, reg: &[QubitId], s: &[QubitId]) {
    b.set_phase("shrunken_pz_shift_right");
    for (j, &ctrl) in s.iter().enumerate().rev() {
        rotate_right_pow2(b, reg, ctrl, 1usize << j);
    }
}

fn narrow_lt(b: &mut B, a: &[QubitId], c: &[QubitId], out: QubitId, lo: usize) {
    let hi = a.len().min(c.len());
    let lo = lo.min(hi.saturating_sub(1));
    primitives::borrow_compare_refs(b, &a[lo..hi], &c[lo..hi], out);
}

fn set_bit_at_s_gated(b: &mut B, q: &[QubitId], s: &[QubitId], active: QubitId) {
    for (i, &target) in q.iter().enumerate() {
        let eq = b.alloc_qubit();
        eq_const_into(b, s, i, eq);
        b.ccx(active, eq, target);
        eq_const_into(b, s, i, eq);
        b.release_zeroed(eq);
    }
}

fn division_substep_windowed(
    b: &mut B,
    a: &[QubitId],
    c: &[QubitId],
    q_div: &[QubitId],
    s_rot: &[QubitId],
    offset: QubitId,
    active: QubitId,
    lo_a: usize,
    lo_b: usize,
    rot_bits: usize,
) {
    let rb = rot_bits.min(s_rot.len());
    let w = s_rot.len();
    clz_diff_body_middle(b, a, c, w, lo_a, lo_b, |b, diff| {
        for j in 0..w {
            b.ccx(active, diff[j], s_rot[j]);
        }
    });
    rotate_left(b, c, &s_rot[..rb]);
    let or = b.alloc_qubit();
    narrow_lt(b, a, c, or, lo_a);
    b.ccx(active, or, offset);
    narrow_lt(b, a, c, or, lo_a);
    b.release_zeroed(or);
    rotate_right(b, c, std::slice::from_ref(&offset));
    ctrl_dec(b, offset, s_rot);
    clz_diff_body_middle(b, a, c, w, lo_a, lo_a, |b, diff| {
        b.ccx(active, diff[0], offset);
    });
    primitives::ctrl_sub(b, active, a, c);
    set_bit_at_s_gated(b, q_div, s_rot, active);
    rotate_right(b, c, &s_rot[..rb]);
    let t = b.alloc_qubits(w);
    xor_const(b, &t, q_div.len());
    let rev: Vec<_> = q_div.iter().rev().copied().collect();
    bit_length_lean(b, &rev, &t, true);
    primitives::ctrl_sub(b, active, s_rot, &t);
    bit_length_lean(b, &rev, &t, false);
    xor_const(b, &t, q_div.len());
    b.free_vec(&t);
}

fn division_substep_windowed_inv(
    b: &mut B,
    a: &[QubitId],
    c: &[QubitId],
    q_div: &[QubitId],
    s_rot: &[QubitId],
    offset: QubitId,
    active: QubitId,
    lo_a: usize,
    lo_b: usize,
    rot_bits: usize,
) {
    let rb = rot_bits.min(s_rot.len());
    let w = s_rot.len();
    let t = b.alloc_qubits(w);
    xor_const(b, &t, q_div.len());
    let rev: Vec<_> = q_div.iter().rev().copied().collect();
    bit_length_lean(b, &rev, &t, true);
    primitives::ctrl_add(b, active, s_rot, &t);
    bit_length_lean(b, &rev, &t, false);
    xor_const(b, &t, q_div.len());
    b.free_vec(&t);
    rotate_left(b, c, &s_rot[..rb]);
    set_bit_at_s_gated(b, q_div, s_rot, active);
    primitives::ctrl_add(b, active, a, c);
    clz_diff_body_middle(b, a, c, w, lo_a, lo_a, |b, diff| {
        b.ccx(active, diff[0], offset);
    });
    ctrl_inc(b, offset, s_rot);
    rotate_left(b, c, std::slice::from_ref(&offset));
    let or = b.alloc_qubit();
    narrow_lt(b, a, c, or, lo_a);
    b.ccx(active, or, offset);
    narrow_lt(b, a, c, or, lo_a);
    b.release_zeroed(or);
    rotate_right(b, c, &s_rot[..rb]);
    clz_diff_body_middle(b, a, c, w, lo_a, lo_b, |b, diff| {
        for j in 0..w {
            b.ccx(active, diff[j], s_rot[j]);
        }
    });
}

fn multiply_substep_windowed(
    b: &mut B,
    a: &[QubitId],
    c: &[QubitId],
    q_mul: &[QubitId],
    s_rot: &[QubitId],
    off: QubitId,
    active: QubitId,
    ca_window: usize,
    cb_window: usize,
    rot_bits: usize,
) {
    let rb = rot_bits.min(s_rot.len());
    let w = s_rot.len();
    let t = b.alloc_qubits(w);
    let rev: Vec<_> = q_mul.iter().rev().copied().collect();
    xor_const(b, &t, q_mul.len());
    bit_length_lean(b, &rev, &t, true);
    for j in 0..w {
        b.ccx(active, t[j], s_rot[j]);
    }
    bit_length_lean(b, &rev, &t, false);
    xor_const(b, &t, q_mul.len());
    b.free_vec(&t);
    set_bit_at_s_gated(b, q_mul, s_rot, active);
    rotate_left(b, c, &s_rot[..rb]);
    primitives::ctrl_add(b, active, a, c);
    clz_diff_body_middle(b, a, c, w, ca_window, ca_window, |b, diff| {
        b.ccx(active, diff[0], off);
    });
    rotate_left(b, c, std::slice::from_ref(&off));
    ctrl_inc(b, off, s_rot);
    let lt = b.alloc_qubit();
    narrow_lt(b, a, c, lt, ca_window);
    b.ccx(active, lt, off);
    narrow_lt(b, a, c, lt, ca_window);
    b.release_zeroed(lt);
    rotate_right(b, c, &s_rot[..rb]);
    clz_diff_body_middle(b, c, a, w, cb_window, ca_window, |b, diff| {
        primitives::ctrl_add(b, active, s_rot, diff);
    });
}

fn multiply_substep_windowed_inv(
    b: &mut B,
    a: &[QubitId],
    c: &[QubitId],
    q_mul: &[QubitId],
    s_rot: &[QubitId],
    off: QubitId,
    active: QubitId,
    ca_window: usize,
    cb_window: usize,
    rot_bits: usize,
) {
    let rb = rot_bits.min(s_rot.len());
    let w = s_rot.len();
    clz_diff_body_middle(b, c, a, w, cb_window, ca_window, |b, diff| {
        primitives::ctrl_sub(b, active, s_rot, diff);
    });
    rotate_left(b, c, &s_rot[..rb]);
    let lt = b.alloc_qubit();
    narrow_lt(b, a, c, lt, ca_window);
    b.ccx(active, lt, off);
    narrow_lt(b, a, c, lt, ca_window);
    b.release_zeroed(lt);
    ctrl_dec(b, off, s_rot);
    rotate_right(b, c, std::slice::from_ref(&off));
    clz_diff_body_middle(b, a, c, w, ca_window, ca_window, |b, diff| {
        b.ccx(active, diff[0], off);
    });
    primitives::ctrl_sub(b, active, a, c);
    rotate_right(b, c, &s_rot[..rb]);
    set_bit_at_s_gated(b, q_mul, s_rot, active);
    let t = b.alloc_qubits(w);
    let rev: Vec<_> = q_mul.iter().rev().copied().collect();
    xor_const(b, &t, q_mul.len());
    bit_length_lean(b, &rev, &t, true);
    for j in 0..w {
        b.ccx(active, t[j], s_rot[j]);
    }
    bit_length_lean(b, &rev, &t, false);
    xor_const(b, &t, q_mul.len());
    b.free_vec(&t);
}

fn gate_hold(
    b: &mut B,
    x: &[QubitId],
    y: &[QubitId],
    active: QubitId,
    g: QubitId,
    body: impl FnOnce(&mut B, QubitId),
) {
    let lt = b.alloc_qubit();
    primitives::borrow_compare_refs(b, x, y, lt);
    b.ccx(lt, active, g);
    body(b, g);
    b.ccx(lt, active, g);
    primitives::borrow_compare_refs(b, x, y, lt);
    b.release_zeroed(lt);
}

fn done_counter_fn(b: &mut B, aa: &[QubitId], qq: &[QubitId], counter: &[QubitId], inverse: bool) {
    let done = b.alloc_qubit();
    let conv = |b: &mut B, done: QubitId| {
        let az = b.alloc_qubit();
        let qz = b.alloc_qubit();
        or_is_zero(b, aa, az);
        or_is_zero(b, qq, qz);
        b.ccx(az, qz, done);
        or_is_zero(b, qq, qz);
        or_is_zero(b, aa, az);
        b.release_zeroed(qz);
        b.release_zeroed(az);
    };
    let cnz = |b: &mut B, done: QubitId| {
        let z = b.alloc_qubit();
        or_nonzero(b, counter, z);
        b.cx(z, done);
        or_nonzero(b, counter, z);
        b.release_zeroed(z);
    };
    if inverse {
        cnz(b, done);
        ctrl_dec(b, done, counter);
        conv(b, done);
    } else {
        conv(b, done);
        ctrl_inc(b, done, counter);
        cnz(b, done);
    }
    b.release_zeroed(done);
}

fn rb(bits: usize) -> usize {
    if bits == 0 {
        1
    } else {
        usize::BITS as usize - bits.leading_zeros() as usize
    }
}

fn pass_step(
    b: &mut B,
    aa: &[QubitId],
    bb: &[QubitId],
    ca: &[QubitId],
    cb: &[QubitId],
    q: &[QubitId],
    counter: &[QubitId],
    parity: QubitId,
    s_rot: &[QubitId],
    off: QubitId,
    i: usize,
    inverse: bool,
) {
    let (lo_a, lo_b, ca_window, cb_window, _) = schedule::reg_los(i);
    let (sdb, s2b) = schedule::shift_bounds(i);
    let swap = |b: &mut B, active: QubitId| {
        let qz = b.alloc_qubit();
        let anz = b.alloc_qubit();
        let t = b.alloc_qubit();
        let g = b.alloc_qubit();
        or_is_zero(b, q, qz);
        or_nonzero(b, aa, anz);
        b.ccx(qz, anz, t);
        b.ccx(t, active, g);
        for j in 0..aa.len() {
            cswap(b, g, aa[j], bb[j]);
        }
        for j in 0..ca.len() {
            cswap(b, g, ca[j], cb[j]);
        }
        b.cx(g, parity);
        b.ccx(t, active, g);
        b.ccx(qz, anz, t);
        or_nonzero(b, aa, anz);
        or_is_zero(b, q, qz);
        b.release_zeroed(g);
        b.release_zeroed(t);
        b.release_zeroed(anz);
        b.release_zeroed(qz);
    };
    if inverse {
        done_counter_fn(b, aa, q, counter, true);
        let active = b.alloc_qubit();
        or_is_zero(b, counter, active);
        swap(b, active);
        let g_div = b.alloc_qubit();
        gate_hold(b, ca, cb, active, g_div, |b, g| {
            division_substep_windowed_inv(b, aa, bb, q, s_rot, off, g, lo_a, lo_b, rb(sdb));
        });
        b.release_zeroed(g_div);
        let g_mul = b.alloc_qubit();
        gate_hold(b, aa, bb, active, g_mul, |b, g| {
            multiply_substep_windowed_inv(
                b,
                ca,
                cb,
                q,
                s_rot,
                off,
                g,
                ca_window,
                cb_window,
                rb(s2b),
            );
        });
        b.release_zeroed(g_mul);
        or_is_zero(b, counter, active);
        b.release_zeroed(active);
    } else {
        let active = b.alloc_qubit();
        or_is_zero(b, counter, active);
        let g_mul = b.alloc_qubit();
        gate_hold(b, aa, bb, active, g_mul, |b, g| {
            multiply_substep_windowed(b, ca, cb, q, s_rot, off, g, ca_window, cb_window, rb(s2b));
        });
        b.release_zeroed(g_mul);
        let g_div = b.alloc_qubit();
        gate_hold(b, ca, cb, active, g_div, |b, g| {
            division_substep_windowed(b, aa, bb, q, s_rot, off, g, lo_a, lo_b, rb(sdb));
        });
        b.release_zeroed(g_div);
        swap(b, active);
        or_is_zero(b, counter, active);
        b.release_zeroed(active);
        done_counter_fn(b, aa, q, counter, false);
    }
}

fn resize(b: &mut B, reg: &mut Vec<QubitId>, target: usize) {
    while reg.len() > target {
        let q = reg.pop().expect("nonempty dynamic register");
        b.free(q);
    }
    while reg.len() < target {
        reg.push(b.alloc_qubit());
    }
}

fn invert_forward(
    b: &mut B,
    aa: &mut Vec<QubitId>,
    bb: &mut Vec<QubitId>,
    ca: &mut Vec<QubitId>,
    cb: &mut Vec<QubitId>,
    q: &mut Vec<QubitId>,
    counter: &[QubitId],
    parity: QubitId,
    s_rot: &[QubitId],
    off: QubitId,
) {
    for i in 0..schedule::SHRUNKEN_PZ_NSTEPS {
        let (wa, wb, wca, wcb, wq) = schedule::reg_widths(i);
        resize(b, aa, wa.max(wb));
        resize(b, bb, wa.max(wb));
        resize(b, ca, wca.max(wcb));
        resize(b, cb, wca.max(wcb));
        resize(b, q, wq.max(1));
        pass_step(b, aa, bb, ca, cb, q, counter, parity, s_rot, off, i, false);
    }
}

fn invert_backward(
    b: &mut B,
    aa: &mut Vec<QubitId>,
    bb: &mut Vec<QubitId>,
    ca: &mut Vec<QubitId>,
    cb: &mut Vec<QubitId>,
    q: &mut Vec<QubitId>,
    counter: &[QubitId],
    parity: QubitId,
    s_rot: &[QubitId],
    off: QubitId,
) {
    for i in (0..schedule::SHRUNKEN_PZ_NSTEPS).rev() {
        pass_step(b, aa, bb, ca, cb, q, counter, parity, s_rot, off, i, true);
        if i > 0 {
            let (wa, wb, wca, wcb, wq) = schedule::reg_widths(i - 1);
            resize(b, aa, wa.max(wb));
            resize(b, bb, wa.max(wb));
            resize(b, ca, wca.max(wcb));
            resize(b, cb, wca.max(wcb));
            resize(b, q, wq.max(1));
        }
    }
}

fn load_p(b: &mut B, reg: &[QubitId]) {
    for (i, &q) in reg.iter().enumerate().take(256) {
        if bit(SECP256K1_P, i) {
            b.x(q);
        }
    }
}

pub(crate) fn divide_forward(
    b: &mut B,
    mut dx: Vec<QubitId>,
    dy: Vec<QubitId>,
) -> (Vec<QubitId>, Vec<QubitId>, Vec<QubitId>) {
    assert_eq!(dx.len(), 257);
    assert_eq!(dy.len(), 257);
    let sgn = b.alloc_qubit();
    compare_geq_const(b, &dx, half_p_ceil(), sgn);
    controlled_field_neg(b, sgn, &dx);
    let (a0, b0, ca0, cb0, q0) = schedule::reg_widths(0);
    let (wg0, wc0) = (a0.max(b0), ca0.max(cb0));
    resize(b, &mut dx, wg0);
    let mut aa = b.alloc_qubits(wg0);
    let mut ca = b.alloc_qubits(wc0);
    let mut cb = b.alloc_qubits(wc0);
    let mut q = b.alloc_qubits(q0.max(1));
    let s_rot = b.alloc_qubits(9);
    let off = b.alloc_qubit();
    let parity = b.alloc_qubit();
    let counter = b.alloc_qubits(10);
    load_p(b, &aa);
    b.x(cb[0]);
    b.x(parity);
    invert_forward(
        b, &mut aa, &mut dx, &mut ca, &mut cb, &mut q, &counter, parity, &s_rot, off,
    );
    let (ta, tb, tca, tq) = (aa.len(), dx.len(), ca.len(), q.len());
    load_p(b, &ca);
    b.x(dx[0]);
    b.free_vec(&aa);
    b.free_vec(&dx);
    b.free_vec(&ca);
    b.free_vec(&q);
    aa.clear();
    dx.clear();
    ca.clear();
    q.clear();
    let cb_w = cb.len();
    resize(b, &mut cb, 257);
    let lambda = b.alloc_qubits(257);
    rfold::mod_mul_rfold_mbu(b, &lambda, &cb[..257], &dy);
    resize(b, &mut cb, cb_w);
    let f = b.alloc_qubit();
    b.cx(sgn, f);
    b.cx(parity, f);
    b.x(f);
    controlled_field_neg(b, f, &lambda);
    b.x(f);
    b.cx(parity, f);
    b.cx(sgn, f);
    b.release_zeroed(f);
    let mut ghosts = Vec::with_capacity(dy.len());
    for &qbit in &dy {
        ghosts.push(b.hmr_ghost(qbit));
    }
    for qbit in dy {
        b.release_zeroed(qbit);
    }
    aa = b.alloc_qubits(ta);
    dx = b.alloc_qubits(tb);
    b.x(dx[0]);
    ca = b.alloc_qubits(tca);
    load_p(b, &ca);
    q = b.alloc_qubits(tq);
    invert_backward(
        b, &mut aa, &mut dx, &mut ca, &mut cb, &mut q, &counter, parity, &s_rot, off,
    );
    b.x(parity);
    b.release_zeroed(parity);
    b.x(cb[0]);
    load_p(b, &aa);
    b.free_vec(&aa);
    b.free_vec(&ca);
    b.free_vec(&cb);
    b.free_vec(&q);
    b.free_vec(&s_rot);
    b.free_vec(&counter);
    b.release_zeroed(off);
    resize(b, &mut dx, 257);
    controlled_field_neg(b, sgn, &dx);
    compare_geq_const(b, &dx, half_p_ceil(), sgn);
    b.release_zeroed(sgn);
    let dy_new = b.alloc_qubits(257);
    rfold::mod_mul_rfold_mbu(b, &dy_new, &lambda, &dx);
    for (ghost, &qbit) in ghosts.into_iter().zip(dy_new.iter()) {
        b.resolve_ghost(ghost, qbit);
    }
    (dx, dy_new, lambda)
}

pub(crate) fn divide_cancel(
    b: &mut B,
    mut dx: Vec<QubitId>,
    dy: Vec<QubitId>,
    lambda: Vec<QubitId>,
) -> (Vec<QubitId>, Vec<QubitId>) {
    assert_eq!(dx.len(), 257);
    assert_eq!(dy.len(), 257);
    assert_eq!(lambda.len(), 257);
    let mut lam_ghosts = Vec::with_capacity(lambda.len());
    for &qbit in &lambda {
        lam_ghosts.push(b.hmr_ghost(qbit));
    }
    for qbit in lambda {
        b.release_zeroed(qbit);
    }
    let sgn = b.alloc_qubit();
    compare_geq_const(b, &dx, half_p_ceil(), sgn);
    controlled_field_neg(b, sgn, &dx);
    let (a0, b0, ca0, cb0, q0) = schedule::reg_widths(0);
    let (wg0, wc0) = (a0.max(b0), ca0.max(cb0));
    resize(b, &mut dx, wg0);
    let mut aa = b.alloc_qubits(wg0);
    let mut ca = b.alloc_qubits(wc0);
    let mut cb = b.alloc_qubits(wc0);
    let mut q = b.alloc_qubits(q0.max(1));
    let s_rot = b.alloc_qubits(9);
    let off = b.alloc_qubit();
    let parity = b.alloc_qubit();
    let counter = b.alloc_qubits(10);
    load_p(b, &aa);
    b.x(cb[0]);
    b.x(parity);
    invert_forward(
        b, &mut aa, &mut dx, &mut ca, &mut cb, &mut q, &counter, parity, &s_rot, off,
    );
    let (ta, tb, tca, tq) = (aa.len(), dx.len(), ca.len(), q.len());
    load_p(b, &ca);
    b.x(dx[0]);
    b.free_vec(&aa);
    b.free_vec(&dx);
    b.free_vec(&ca);
    b.free_vec(&q);
    aa.clear();
    dx.clear();
    ca.clear();
    q.clear();
    let cb_w = cb.len();
    resize(b, &mut cb, 257);
    let temp = b.alloc_qubits(257);
    rfold::mod_mul_rfold_mbu(b, &temp, &cb[..257], &dy);
    let f = b.alloc_qubit();
    b.cx(sgn, f);
    b.cx(parity, f);
    b.x(f);
    controlled_field_neg(b, f, &temp);
    for (ghost, &qbit) in lam_ghosts.into_iter().zip(temp.iter()) {
        b.resolve_ghost(ghost, qbit);
    }
    controlled_field_neg(b, f, &temp);
    b.x(f);
    b.cx(parity, f);
    b.cx(sgn, f);
    b.release_zeroed(f);
    rfold::mod_mul_rfold_mbu_undo(b, &temp, &cb[..257], &dy);
    b.free_vec(&temp);
    resize(b, &mut cb, cb_w);
    aa = b.alloc_qubits(ta);
    dx = b.alloc_qubits(tb);
    b.x(dx[0]);
    ca = b.alloc_qubits(tca);
    load_p(b, &ca);
    q = b.alloc_qubits(tq);
    invert_backward(
        b, &mut aa, &mut dx, &mut ca, &mut cb, &mut q, &counter, parity, &s_rot, off,
    );
    b.x(parity);
    b.release_zeroed(parity);
    b.x(cb[0]);
    load_p(b, &aa);
    b.free_vec(&aa);
    b.free_vec(&ca);
    b.free_vec(&cb);
    b.free_vec(&q);
    b.free_vec(&s_rot);
    b.free_vec(&counter);
    b.release_zeroed(off);
    resize(b, &mut dx, 257);
    controlled_field_neg(b, sgn, &dx);
    compare_geq_const(b, &dx, half_p_ceil(), sgn);
    b.release_zeroed(sgn);
    (dx, dy)
}

pub(crate) fn component_peak_probe() -> (u32, &'static str) {
    let mut b = B::new_count_only();
    let dx = b.alloc_qubits(257);
    let dy = b.alloc_qubits(257);
    let (dx, dy, lambda) = divide_forward(&mut b, dx, dy);
    let (_dx, _dy) = divide_cancel(&mut b, dx, dy, lambda);
    (b.peak_qubits, b.peak_phase)
}

pub(crate) fn compare_p_selftest() -> Result<(), String> {
    use sha3::{
        digest::{ExtendableOutput, Update},
        Shake256,
    };

    let mut b = B::new();
    let a = b.alloc_qubits(257);
    let flag = b.alloc_qubit();
    compare_geq_p_secp256k1_into(&mut b, &a, flag);

    let mut seed = Shake256::default();
    seed.update(b"shrunken-pz-compare-p");
    let mut xof = seed.finalize_xof();
    let mut sim = crate::sim::Simulator::new(b.next_qubit as usize, b.next_bit as usize, &mut xof);

    let p = SECP256K1_P;
    let cases = [
        (U256::ZERO, false),
        (U256::from(1u64), false),
        (p.wrapping_sub(U256::from(2u64)), false),
        (p.wrapping_sub(U256::from(1u64)), false),
        (p, false),
        (p.wrapping_add(U256::from(1u64)), false),
        (U256::MAX, false),
        (U256::ZERO, true),
    ];
    for shot in 0..64usize {
        let (low, high) = cases[shot % cases.len()];
        for i in 0..256 {
            if bit(low, i) {
                *sim.qubit_mut(a[i]) |= 1u64 << shot;
            }
        }
        if high {
            *sim.qubit_mut(a[256]) |= 1u64 << shot;
        }
    }

    sim.apply_iter(b.ops.iter());
    if sim.phase != 0 {
        return Err(format!("compare_p_selftest phase garbage: 0x{:016x}", sim.phase));
    }
    for shot in 0..64usize {
        let (low, high) = cases[shot % cases.len()];
        let expect = high || low >= p;
        let got = ((sim.qubit(flag) >> shot) & 1) != 0;
        if got != expect {
            return Err(format!(
                "shot {shot}: compare_p got {}, expected {} for high={} low=0x{low:x}",
                u8::from(got),
                u8::from(expect),
                u8::from(high),
            ));
        }
    }
    for q in 258..b.next_qubit as usize {
        if sim.qubit(QubitId(q as u64)) != 0 {
            return Err(format!("compare_p dirty ancilla q{q}"));
        }
    }
    Ok(())
}

pub(crate) fn ec_add_inplace_shrunken_pz(
    b: &mut B,
    tx: &mut Vec<QubitId>,
    ty: &mut Vec<QubitId>,
    ox: &[BitId],
    oy: &[BitId],
) {
    assert_eq!(tx.len(), 256);
    assert_eq!(ty.len(), 256);
    assert_eq!(ox.len(), 256);
    assert_eq!(oy.len(), 256);

    b.set_phase("shrunken_pz_ec_pad");
    tx.push(b.alloc_qubit());
    ty.push(b.alloc_qubit());

    b.set_phase("shrunken_pz_ec_dy_build");
    mod_sub_from_creg_qload(b, &ty, oy);
    b.set_phase("shrunken_pz_ec_dx_build");
    mod_sub_from_creg_qload(b, &tx, ox);

    b.set_phase("shrunken_pz_ec_inv_forward");
    let tx_inner = std::mem::take(tx);
    let ty_inner = std::mem::take(ty);
    let (tx_inner, ty_inner, lambda) = divide_forward(b, tx_inner, ty_inner);
    *tx = tx_inner;
    *ty = ty_inner;

    b.set_phase("shrunken_pz_ec_dx_clean");
    mod_sub_from_creg_qload(b, &tx, ox);

    b.set_phase("shrunken_pz_ec_new_x");
    field_neg(b, &tx);
    mod_mac_inplace(b, &tx, &lambda, &lambda);
    mod_sub_creg_qload(b, &tx, ox);

    b.set_phase("shrunken_pz_ec_dx_diff");
    field_neg(b, &tx);
    rfold::mod_double_rfold_mbu(b, &tx);
    mod_mac_inplace(b, &tx, &lambda, &lambda);
    mod_sub_creg_qload(b, &tx, ox);

    b.set_phase("shrunken_pz_ec_new_y");
    mod_mac_inplace(b, &ty, &lambda, &tx);
    mod_sub_creg_qload(b, &ty, oy);

    b.set_phase("shrunken_pz_ec_dx_diff_clean");
    mod_add_creg_qload(b, &tx, ox);
    mod_msc_inplace(b, &tx, &lambda, &lambda);
    rfold::mod_halve_rfold_mbu(b, &tx);
    field_neg(b, &tx);

    b.set_phase("shrunken_pz_ec_alt_new_dy");
    mod_add_creg_qload(b, &ty, oy);
    b.set_phase("shrunken_pz_ec_alt_new_dx");
    mod_sub_from_creg_qload(b, &tx, ox);
    b.set_phase("shrunken_pz_ec_inv_cancel");
    let tx_inner = std::mem::take(tx);
    let ty_inner = std::mem::take(ty);
    let (tx_inner, ty_inner) = divide_cancel(b, tx_inner, ty_inner, lambda);
    *tx = tx_inner;
    *ty = ty_inner;

    b.set_phase("shrunken_pz_ec_restore_x");
    mod_sub_from_creg_qload(b, &tx, ox);
    b.set_phase("shrunken_pz_ec_restore_y");
    mod_sub_creg_qload(b, &ty, oy);

    b.set_phase("shrunken_pz_ec_unpad");
    b.release_zeroed(ty.pop().expect("ty padded"));
    b.release_zeroed(tx.pop().expect("tx padded"));
    b.set_phase("shrunken_pz_ec_done");
}
