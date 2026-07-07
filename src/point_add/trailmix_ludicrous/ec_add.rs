
use super::arith::{
    mod_add, mod_add_exact, mod_neg, mod_rsub_vented_loaded, mod_sub_classical_low3,
    mod_sub_shifted_low, mod_sub_vented,
};
use super::gcd::{mod_mul_inverse_in_place, Direction};
use super::square::mod_square_sub_pm_secp256k1_symmetric;
use super::{B, BExt};
use crate::point_add::{arith::mod_const_minus_reg_qb, SECP256K1_P};
use crate::circuit::{BitId, QubitId};

const N: usize = 256;

fn coord_addsub(circ: &mut B, dst: &[QubitId], coord: &[BitId], subtract: bool) {
    debug_assert_eq!(dst.len(), N);
    debug_assert_eq!(coord.len(), N);
    let split_low3 = subtract
        && std::env::var("TLM_COORD_SPLIT_LOW3")
            .ok()
            .as_deref()
            .unwrap_or("0")
            != "0";
    if split_low3 {
        let temp = circ.alloc_qubits(N - 3);
        for i in 3..N {
            circ.x_if_bit(temp[i - 3], coord[i]);
        }
        mod_sub_shifted_low(circ, &temp, dst, 3);
        for i in 3..N {
            circ.x_if_bit(temp[i - 3], coord[i]);
        }
        for q in temp {
            circ.zero_and_free(q);
        }
        mod_sub_classical_low3(circ, dst, &coord[..3]);
        return;
    }
    let temp = circ.alloc_qubits(N);
    for i in 0..N {
        circ.x_if_bit(temp[i], coord[i]);
    }

    if subtract {
        mod_sub_vented(circ, &temp, dst);
    } else {
        mod_add(circ, &temp, dst);
    }
    for i in 0..N {
        circ.x_if_bit(temp[i], coord[i]);
    }
    for q in temp {
        circ.zero_and_free(q);
    }
}

fn coord_add3x(circ: &mut B, dst: &[QubitId], coord: &[BitId]) {
    debug_assert_eq!(dst.len(), N);
    debug_assert_eq!(coord.len(), N);

    let three_coord = classical_times3_mod_q(circ, coord);

    let temp = circ.alloc_qubits(N);
    for i in 0..N {
        circ.x_if_bit(temp[i], three_coord[i]);
    }

    if std::env::var("TLM_COORD_ADD3X_TRUNC").ok().as_deref() == Some("1") {
        mod_add(circ, &temp, dst);
    } else {
        mod_add_exact(circ, &temp, dst);
    }
    for i in 0..N {
        circ.x_if_bit(temp[i], three_coord[i]);
    }
    for q in temp {
        circ.zero_and_free(q);
    }

    for &b in &three_coord {
        circ.bit_store0(b);
    }
}

fn classical_times3_mod_q(circ: &mut B, coord: &[BitId]) -> Vec<BitId> {
    debug_assert_eq!(coord.len(), N);
    const C: u128 = (1u128 << 32) + 977;

    let s: Vec<BitId> = circ.alloc_bits(N + 2);
    for &b in &s {
        circ.bit_store0(b);
    }
    classical_add_into(circ, &s, coord);
    classical_add_into(circ, &s, coord);
    classical_add_into(circ, &s, coord);

    let r: Vec<BitId> = circ.alloc_bits(N + 1);
    for i in 0..N {
        circ.bit_copy(r[i], s[i]);
    }
    circ.bit_store0(r[N]);

    let av_bits = 35usize;
    let av: Vec<BitId> = circ.alloc_bits(av_bits);
    classical_set_const_times_bit(circ, &av, C, s[N], false);
    classical_add_const_times_bit(circ, &av, 2 * C, s[N + 1]);
    classical_add_into(circ, &r, &av);

    let tmp: Vec<BitId> = circ.alloc_bits(N + 2);
    for i in 0..(N + 1) {
        circ.bit_copy(tmp[i], r[i]);
    }
    circ.bit_store0(tmp[N + 1]);
    {

        let cbits: Vec<BitId> = circ.alloc_bits(av_bits);
        classical_set_const(circ, &cbits, C);
        classical_add_into(circ, &tmp, &cbits);
        for &b in &cbits {
            circ.bit_store0(b);
        }
    }

    let geflag = circ.alloc_bit();
    circ.bit_store0(geflag);
    circ.push_condition(tmp[N]);
    circ.bit_store1(geflag);
    circ.pop_condition();
    circ.push_condition(tmp[N + 1]);
    circ.bit_store1(geflag);
    circ.pop_condition();

    let result: Vec<BitId> = circ.alloc_bits(N);
    for i in 0..N {
        circ.bit_store0(result[i]);
        circ.push_condition(geflag);
        circ.push_condition(tmp[i]);
        circ.bit_store1(result[i]);
        circ.pop_condition();
        circ.pop_condition();
        circ.bit_invert(geflag);
        circ.push_condition(geflag);
        circ.push_condition(r[i]);
        circ.bit_store1(result[i]);
        circ.pop_condition();
        circ.pop_condition();
        circ.bit_invert(geflag);
    }

    circ.bit_store0(geflag);
    for &b in tmp.iter().chain(av.iter()).chain(r.iter()).chain(s.iter()) {
        circ.bit_store0(b);
    }
    result
}

fn classical_set_const(circ: &mut B, dst: &[BitId], k: u128) {
    for (i, &b) in dst.iter().enumerate() {
        let bit = i < 128 && ((k >> i) & 1) == 1;
        if bit {
            circ.bit_store0(b);
            circ.bit_invert(b);
        } else {
            circ.bit_store0(b);
        }
    }
}

fn classical_set_const_times_bit(circ: &mut B, dst: &[BitId], k: u128, gate: BitId, _accumulate: bool) {
    for (i, &b) in dst.iter().enumerate() {
        circ.bit_store0(b);
        let bit = i < 128 && ((k >> i) & 1) == 1;
        if bit {
            circ.push_condition(gate);
            circ.bit_store1(b);
            circ.pop_condition();
        }
    }
}

fn classical_add_const_times_bit(circ: &mut B, dst: &[BitId], k: u128, gate: BitId) {
    let w = dst.len();
    let addend: Vec<BitId> = circ.alloc_bits(w);
    classical_set_const_times_bit(circ, &addend, k, gate, false);
    classical_add_into(circ, dst, &addend);
    for &b in &addend {
        circ.bit_store0(b);
    }
}

fn classical_add_into(circ: &mut B, acc: &[BitId], addend: &[BitId]) {
    let carry = circ.alloc_bit();
    circ.bit_store0(carry);
    let newcarry = circ.alloc_bit();
    for i in 0..acc.len() {
        let a_i = addend.get(i).copied();

        circ.bit_store0(newcarry);
        if let Some(a) = a_i {
            circ.bit_and_xor_into(newcarry, acc[i], a);
            circ.bit_and_xor_into(newcarry, acc[i], carry);
            circ.bit_and_xor_into(newcarry, a, carry);
        } else {
            circ.bit_and_xor_into(newcarry, acc[i], carry);
        }

        if let Some(a) = a_i {
            circ.bit_xor_into(acc[i], a);
        }
        circ.bit_xor_into(acc[i], carry);

        circ.bit_copy(carry, newcarry);
    }
    circ.bit_store0(newcarry);
    circ.bit_store0(carry);
}

fn classical_plus1_mod_2n(circ: &mut B, coord: &[BitId]) -> Vec<BitId> {
    debug_assert_eq!(coord.len(), N);
    let s: Vec<BitId> = circ.alloc_bits(N);
    for i in 0..N {
        circ.bit_copy(s[i], coord[i]);
    }
    let one: Vec<BitId> = circ.alloc_bits(1);
    circ.bit_store0(one[0]);
    circ.bit_invert(one[0]);
    classical_add_into(circ, &s, &one);
    circ.bit_store0(one[0]);
    s
}

fn coord_rsub(circ: &mut B, x: &[QubitId], coord: &[BitId]) {
    debug_assert_eq!(x.len(), N);
    debug_assert_eq!(coord.len(), N);

    if std::env::var("TLM_COORD_RSUB_FUSED").ok().as_deref() == Some("1") {
        let coord_p1 = classical_plus1_mod_2n(circ, coord);
        let t: Vec<QubitId> = (0..N).map(|_| circ.alloc_qubit()).collect();
        for i in 0..N {
            circ.x_if_bit(t[i], coord_p1[i]);
        }
        mod_rsub_vented_loaded(circ, &t, x);
        for i in 0..N {
            circ.x_if_bit(t[i], coord_p1[i]);
        }
        for q in t {
            circ.zero_and_free(q);
        }
        for &b in &coord_p1 {
            circ.bit_store0(b);
        }
        return;
    }
    if std::env::var("TLM_FUSE_X_RESTORE")
        .ok()
        .as_deref()
        == Some("1")
    {
        mod_const_minus_reg_qb(circ, x, coord, SECP256K1_P);
        return;
    }
    let t: Vec<QubitId> = (0..N).map(|_| circ.alloc_qubit()).collect();
    for i in 0..N {
        circ.x_if_bit(t[i], coord[i]);
    }
    mod_sub_vented(circ, &t, x);
    for i in 0..N {
        circ.x_if_bit(t[i], coord[i]);
    }
    for q in t {
        circ.zero_and_free(q);
    }
    mod_neg(circ, x);
}

pub fn ec_add(
    circ: &mut B,
    x2: &mut Vec<QubitId>,
    y2: &[QubitId],
    ox: &[BitId],
    oy: &[BitId],
) {
    assert_eq!(x2.len(), N, "x2 is 256 bits");
    assert_eq!(y2.len(), N, "y2 is 256 bits");
    assert_eq!(ox.len(), N, "ox is 256 classical bits");
    assert_eq!(oy.len(), N, "oy is 256 classical bits");

    circ.set_phase("tlm_coord_x_sub");
    coord_addsub(circ, x2, ox, true);
    circ.set_phase("tlm_coord_y_sub");
    coord_addsub(circ, &y2[..N], oy, true);

    circ.set_phase("tlm_inverse");
    let xv = std::mem::take(x2);
    *x2 = mod_mul_inverse_in_place(circ, xv, y2, Direction::Inverse);

    circ.set_phase("tlm_coord_add3x");
    coord_add3x(circ, x2, ox);

    circ.set_phase("tlm_square");
    mod_square_sub_pm_secp256k1_symmetric(circ, &y2[..N], x2);

    circ.set_phase("tlm_forward_multiply");
    let xv = std::mem::take(x2);
    *x2 = mod_mul_inverse_in_place(circ, xv, y2, Direction::Forward);

    circ.set_phase("tlm_coord_y_sub_final");
    coord_addsub(circ, &y2[..N], oy, true);
    circ.set_phase("tlm_coord_rsub_final");
    coord_rsub(circ, x2, ox);
}

pub fn build_times3_test() -> (Vec<crate::circuit::Op>, Vec<BitId>, Vec<BitId>) {
    let mut circ = B::new_for_test();
    let ox = circ.alloc_bits(N);
    let t = classical_times3_mod_q(&mut circ, &ox);

    let tout = circ.alloc_bits(N);
    for i in 0..N {
        circ.bit_copy(tout[i], t[i]);
    }
    circ.declare_bit_register(&ox);
    circ.declare_bit_register(&tout);
    (circ.take_ops(), ox, tout)
}

pub fn build_add3x_test() -> (Vec<crate::circuit::Op>, Vec<BitId>, Vec<QubitId>) {
    let mut circ = B::new_for_test();
    let dst: Vec<QubitId> = circ.alloc_qubits(N);
    let ox = circ.alloc_bits(N);
    coord_add3x(&mut circ, &dst, &ox);
    circ.declare_qubit_register(&dst);
    circ.declare_bit_register(&ox);
    (circ.take_ops(), ox, dst)
}

fn coord_add3x_orig(circ: &mut B, dst: &[QubitId], coord: &[BitId]) {
    let temp: Vec<QubitId> = (0..=N).map(|_| circ.alloc_qubit()).collect();
    for i in 0..N {
        circ.x_if_bit(temp[i], coord[i]);
    }
    mod_add(circ, &temp[..N], dst);
    super::arith::mod_double(circ, &temp);
    mod_add(circ, &temp[..N], dst);
    super::arith::mod_double_reverse(circ, &temp);
    for i in 0..N {
        circ.x_if_bit(temp[i], coord[i]);
    }
    for q in temp {
        circ.zero_and_free(q);
    }
}

pub fn build_add3x_test_orig() -> (Vec<crate::circuit::Op>, Vec<BitId>, Vec<QubitId>) {
    let mut circ = B::new_for_test();
    let dst: Vec<QubitId> = circ.alloc_qubits(N);
    let ox = circ.alloc_bits(N);
    coord_add3x_orig(&mut circ, &dst, &ox);
    circ.declare_qubit_register(&dst);
    circ.declare_bit_register(&ox);
    (circ.take_ops(), ox, dst)
}
