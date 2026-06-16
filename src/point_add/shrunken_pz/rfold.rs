use alloy_primitives::U256;

use crate::circuit::{BitId, QubitId};

use super::super::{
    add_nbit_qq, cadd_nbit_const_direct_fast, cmp_lt_into, csub_nbit_const_direct_fast, B,
};

const RFOLD_WINDOW: usize = 73;
const COMPARE_TOPK: usize = 64;

fn r_const() -> U256 {
    U256::from(977u64) + (U256::from(1u64) << 32)
}

fn top_window<'a>(xs: &'a [QubitId], width: usize) -> &'a [QubitId] {
    let n = xs.len().min(256);
    let w = width.min(n);
    &xs[n - w..n]
}

fn phase_correct_lt(b: &mut B, u: &[QubitId], v: &[QubitId], ctrl: Option<QubitId>, phase: BitId) {
    let flag = b.alloc_qubit();
    b.push_condition(phase);
    cmp_lt_into(b, u, v, flag);
    match ctrl {
        Some(ctrl) => b.cz(ctrl, flag),
        None => b.cz(flag, flag),
    }
    cmp_lt_into(b, u, v, flag);
    b.pop_condition();
    b.release_zeroed(flag);
}

fn left_shift_one(b: &mut B, a: &[QubitId]) {
    for i in (1..a.len()).rev() {
        b.swap(a[i], a[i - 1]);
    }
}

fn right_shift_one(b: &mut B, a: &[QubitId]) {
    for i in 0..a.len() - 1 {
        b.swap(a[i], a[i + 1]);
    }
}

fn add_with_carry_to_high(b: &mut B, a: &[QubitId], addend: &[QubitId]) {
    assert_eq!(a.len(), addend.len() + 1);
    let pad = b.alloc_qubit();
    let mut add_ext = addend.to_vec();
    add_ext.push(pad);
    add_nbit_qq(b, &add_ext, a);
    b.release_zeroed(pad);
}

pub(crate) fn controlled_mod_add_rfold_mbu(
    b: &mut B,
    ctrl: QubitId,
    acc: &[QubitId],
    addend: &[QubitId],
) {
    assert_eq!(acc.len(), 257);
    assert!(addend.len() == 256 || addend.len() == 257);
    b.set_phase("shrunken_pz_rfold_cadd_int");
    if addend.len() == 257 {
        super::primitives::ctrl_add(b, ctrl, acc, addend);
    } else {
        let pad = b.alloc_qubit();
        let mut add_ext = addend.to_vec();
        add_ext.push(pad);
        super::primitives::ctrl_add(b, ctrl, acc, &add_ext);
        b.release_zeroed(pad);
    }

    b.set_phase("shrunken_pz_rfold_cadd_rfold");
    cadd_nbit_const_direct_fast(b, &acc[..RFOLD_WINDOW], r_const(), acc[256]);

    b.set_phase("shrunken_pz_rfold_cadd_phase");
    let phase = b.alloc_bit();
    b.hmr(acc[256], phase);
    phase_correct_lt(
        b,
        top_window(acc, COMPARE_TOPK),
        top_window(addend, COMPARE_TOPK),
        Some(ctrl),
        phase,
    );
}

pub(crate) fn controlled_mod_sub_rfold_mbu(
    b: &mut B,
    ctrl: QubitId,
    acc: &[QubitId],
    subtrahend: &[QubitId],
) {
    assert_eq!(acc.len(), 257);
    b.set_phase("shrunken_pz_rfold_csub");
    for &q in &acc[..256] {
        b.x(q);
    }
    controlled_mod_add_rfold_mbu(b, ctrl, acc, subtrahend);
    for &q in &acc[..256] {
        b.x(q);
    }
}

pub(crate) fn mod_double_rfold_mbu(b: &mut B, acc: &[QubitId]) {
    assert_eq!(acc.len(), 257);
    b.set_phase("shrunken_pz_rfold_double_shift");
    left_shift_one(b, acc);
    b.set_phase("shrunken_pz_rfold_double_rfold");
    cadd_nbit_const_direct_fast(b, &acc[..RFOLD_WINDOW], r_const(), acc[256]);
    let phase = b.alloc_bit();
    b.hmr(acc[256], phase);
    b.z_if(acc[0], phase);
}

pub(crate) fn mod_halve_rfold_mbu(b: &mut B, acc: &[QubitId]) {
    assert_eq!(acc.len(), 257);
    b.set_phase("shrunken_pz_rfold_halve");
    b.cx(acc[0], acc[256]);
    for &q in &acc[..RFOLD_WINDOW] {
        b.x(q);
    }
    cadd_nbit_const_direct_fast(b, &acc[..RFOLD_WINDOW], r_const(), acc[256]);
    for &q in &acc[..RFOLD_WINDOW] {
        b.x(q);
    }
    right_shift_one(b, acc);
}

pub(crate) fn reduce_once_secp256k1_from_rfold(b: &mut B, acc: &[QubitId], flag: QubitId) {
    assert_eq!(acc.len(), 257);
    b.set_phase("shrunken_pz_rfold_reduce");
    let p = super::super::SECP256K1_P;
    let p_reg = b.alloc_qubits(257);
    for i in 0..256 {
        if super::super::bit(p, i) {
            b.x(p_reg[i]);
        }
    }
    cmp_lt_into(b, acc, &p_reg, flag);
    b.x(flag);
    cadd_nbit_const_direct_fast(b, &acc[..256], r_const(), flag);
    b.x(flag);
    for i in 0..256 {
        if super::super::bit(p, i) {
            b.x(p_reg[i]);
        }
    }
    b.free_vec(&p_reg);
}

pub(crate) fn reduce_once_secp256k1_from_rfold_reverse(b: &mut B, acc: &[QubitId], flag: QubitId) {
    assert_eq!(acc.len(), 257);
    b.set_phase("shrunken_pz_rfold_reduce_rev");
    csub_nbit_const_direct_fast(b, &acc[..256], r_const(), flag);
    let p = super::super::SECP256K1_P;
    let p_reg = b.alloc_qubits(257);
    for i in 0..256 {
        if super::super::bit(p, i) {
            b.x(p_reg[i]);
        }
    }
    cmp_lt_into(b, acc, &p_reg, flag);
    b.x(flag);
    for i in 0..256 {
        if super::super::bit(p, i) {
            b.x(p_reg[i]);
        }
    }
    b.free_vec(&p_reg);
}

pub(crate) fn mod_mul_rfold_mbu(b: &mut B, result: &[QubitId], a: &[QubitId], x: &[QubitId]) {
    assert_eq!(result.len(), 257);
    assert_eq!(a.len(), 257);
    assert_eq!(x.len(), 257);
    controlled_mod_add_rfold_mbu(b, x[255], result, a);
    for i in (0..255).rev() {
        mod_double_rfold_mbu(b, result);
        controlled_mod_add_rfold_mbu(b, x[i], result, a);
    }
}

pub(crate) fn mod_mul_rfold_mbu_undo(b: &mut B, result: &[QubitId], a: &[QubitId], x: &[QubitId]) {
    assert_eq!(result.len(), 257);
    assert_eq!(a.len(), 257);
    assert_eq!(x.len(), 257);
    for i in 0..255 {
        controlled_mod_sub_rfold_mbu(b, x[i], result, a);
        mod_halve_rfold_mbu(b, result);
    }
    controlled_mod_sub_rfold_mbu(b, x[255], result, a);
}
