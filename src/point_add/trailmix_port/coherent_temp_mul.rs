//! Gate-only canonical multiplication for short-lived products.
//!
//! Unlike the canonical MBU multiplier, this route coherently clears every
//! reduction predicate. Its forward block therefore contains only
//! self-inverse quantum gates. The matching inverse replays the exact emitted
//! `Op` records in reverse order; it does not rebuild an algebraic undo.

use std::collections::BTreeSet;

use crate::circuit::{Op, OperationType, QubitId, NO_BIT};
use crate::point_add::trailmix_port::circuit::{Circuit, QReg};

const OP_KIND_COUNT: usize = 18;

#[derive(Debug)]
pub struct CanonicalCoherentMulInverse {
    count_only: bool,
    op_start: usize,
    gate_count: usize,
    kind_counts: [usize; OP_KIND_COUNT],
    workspace_ids: Vec<u32>,
    result_ids: Vec<u32>,
    a_ids: Vec<u32>,
    b_ids: Vec<u32>,
}

impl CanonicalCoherentMulInverse {
    #[must_use]
    pub fn gate_count(&self) -> usize {
        self.gate_count
    }

    #[must_use]
    pub fn toffoli_count(&self) -> usize {
        self.kind_counts[OperationType::CCX as usize]
            + self.kind_counts[OperationType::CCZ as usize]
    }

    #[must_use]
    pub fn workspace_qubits(&self) -> usize {
        self.workspace_ids.len()
    }
}

fn ids(reg: &[QReg]) -> Vec<u32> {
    reg.iter().map(QReg::id).collect()
}

fn assert_register_alias_rules(result: &[QReg], a: &[QReg], b: &[QReg]) {
    let result_ids: BTreeSet<_> = ids(result).into_iter().collect();
    let a_ids: BTreeSet<_> = ids(a).into_iter().collect();
    let b_ids: BTreeSet<_> = ids(b).into_iter().collect();
    assert!(
        result_ids.is_disjoint(&a_ids) && result_ids.is_disjoint(&b_ids),
        "coherent multiplier result must not alias an operand"
    );
    assert!(
        a_ids.is_disjoint(&b_ids) || a_ids == b_ids,
        "coherent multiplier operands may be disjoint or exact aliases only"
    );
}

fn assert_modulus(width: usize, modulus_le: &[u8]) {
    assert!(width > 0, "coherent multiplier width must be positive");
    assert_eq!(
        modulus_le.first().copied().unwrap_or(0) & 1,
        1,
        "coherent multiplier requires an odd modulus"
    );
    assert!(
        (width..modulus_le.len() * 8)
            .all(|i| ((modulus_le[i / 8] >> (i % 8)) & 1) == 0),
        "coherent multiplier modulus does not fit its value width"
    );
}

fn is_literal_gate(op: &Op) -> bool {
    matches!(
        op.kind,
        OperationType::X
            | OperationType::Z
            | OperationType::CX
            | OperationType::CZ
            | OperationType::Swap
            | OperationType::CCX
            | OperationType::CCZ
    ) && op.c_target == NO_BIT
        && op.c_condition == NO_BIT
}

fn assert_gate_alphabet(kind_counts: &[usize; OP_KIND_COUNT]) {
    for kind in 0..OP_KIND_COUNT {
        let allowed = matches!(
            kind,
            x if x == OperationType::X as usize
                || x == OperationType::Z as usize
                || x == OperationType::CX as usize
                || x == OperationType::CZ as usize
                || x == OperationType::Swap as usize
                || x == OperationType::CCX as usize
                || x == OperationType::CCZ as usize
        );
        assert!(
            allowed || kind_counts[kind] == 0,
            "coherent multiplier emitted forbidden operation kind {kind}"
        );
    }
}

fn coherent_compare_geq_p_secp256k1(circ: &mut Circuit, value: &[QReg], flag: &QReg) {
    use crate::point_add::trailmix_port::arith::khattar_gidney::xor_and_of_khattar_gidney;

    assert_eq!(value.len(), 257);
    let low4_all_ones = circ.alloc_qreg("coherent-cmp-p.low4-all-ones");
    xor_and_of_khattar_gidney(circ, &value[..4], &low4_all_ones);

    let low_tail_or = circ.alloc_qreg("coherent-cmp-p.low-tail-or");
    circ.cx(&value[4], &low_tail_or);
    circ.cx(&low4_all_ones, &low_tail_or);
    circ.ccx(&value[4], &low4_all_ones, &low_tail_or);

    let low6_ge = circ.alloc_qreg("coherent-cmp-p.low6-ge");
    circ.ccx(&value[5], &low_tail_or, &low6_ge);

    let hi_or_67 = circ.alloc_qreg("coherent-cmp-p.hi-or-67");
    circ.cx(&value[6], &hi_or_67);
    circ.cx(&value[7], &hi_or_67);
    circ.ccx(&value[6], &value[7], &hi_or_67);

    let hi_or_89 = circ.alloc_qreg("coherent-cmp-p.hi-or-89");
    circ.cx(&value[8], &hi_or_89);
    circ.cx(&value[9], &hi_or_89);
    circ.ccx(&value[8], &value[9], &hi_or_89);

    let hi4_nonzero = circ.alloc_qreg("coherent-cmp-p.hi4-nonzero");
    circ.cx(&hi_or_67, &hi4_nonzero);
    circ.cx(&hi_or_89, &hi4_nonzero);
    circ.ccx(&hi_or_67, &hi_or_89, &hi4_nonzero);

    let low10_ge = circ.alloc_qreg("coherent-cmp-p.low10-ge");
    circ.cx(&low6_ge, &low10_ge);
    circ.cx(&hi4_nonzero, &low10_ge);
    circ.ccx(&low6_ge, &hi4_nonzero, &low10_ge);

    let mid_all_ones = circ.alloc_qreg("coherent-cmp-p.mid-all-ones");
    xor_and_of_khattar_gidney(circ, &value[10..32], &mid_all_ones);

    let mid_and_low = circ.alloc_qreg("coherent-cmp-p.mid-and-low");
    circ.ccx(&mid_all_ones, &low10_ge, &mid_and_low);

    let tail_or = circ.alloc_qreg("coherent-cmp-p.tail-or");
    circ.cx(&value[32], &tail_or);
    circ.cx(&mid_and_low, &tail_or);
    circ.ccx(&value[32], &mid_and_low, &tail_or);

    let high_all_ones = circ.alloc_qreg("coherent-cmp-p.high-all-ones");
    xor_and_of_khattar_gidney(circ, &value[33..256], &high_all_ones);

    let high_and_tail = circ.alloc_qreg("coherent-cmp-p.high-and-tail");
    circ.ccx(&high_all_ones, &tail_or, &high_and_tail);

    circ.cx(&value[256], flag);
    circ.cx(&high_and_tail, flag);
    circ.ccx(&value[256], &high_and_tail, flag);

    circ.ccx(&high_all_ones, &tail_or, &high_and_tail);
    circ.zero_and_free(high_and_tail);
    xor_and_of_khattar_gidney(circ, &value[33..256], &high_all_ones);
    circ.zero_and_free(high_all_ones);

    circ.ccx(&value[32], &mid_and_low, &tail_or);
    circ.cx(&mid_and_low, &tail_or);
    circ.cx(&value[32], &tail_or);
    circ.zero_and_free(tail_or);

    circ.ccx(&mid_all_ones, &low10_ge, &mid_and_low);
    circ.zero_and_free(mid_and_low);
    xor_and_of_khattar_gidney(circ, &value[10..32], &mid_all_ones);
    circ.zero_and_free(mid_all_ones);

    circ.ccx(&low6_ge, &hi4_nonzero, &low10_ge);
    circ.cx(&hi4_nonzero, &low10_ge);
    circ.cx(&low6_ge, &low10_ge);
    circ.zero_and_free(low10_ge);

    circ.ccx(&hi_or_67, &hi_or_89, &hi4_nonzero);
    circ.cx(&hi_or_89, &hi4_nonzero);
    circ.cx(&hi_or_67, &hi4_nonzero);
    circ.zero_and_free(hi4_nonzero);

    circ.ccx(&value[8], &value[9], &hi_or_89);
    circ.cx(&value[9], &hi_or_89);
    circ.cx(&value[8], &hi_or_89);
    circ.zero_and_free(hi_or_89);

    circ.ccx(&value[6], &value[7], &hi_or_67);
    circ.cx(&value[7], &hi_or_67);
    circ.cx(&value[6], &hi_or_67);
    circ.zero_and_free(hi_or_67);

    circ.ccx(&value[5], &low_tail_or, &low6_ge);
    circ.zero_and_free(low6_ge);

    circ.ccx(&value[4], &low4_all_ones, &low_tail_or);
    circ.cx(&low4_all_ones, &low_tail_or);
    circ.cx(&value[4], &low_tail_or);
    circ.zero_and_free(low_tail_or);
    xor_and_of_khattar_gidney(circ, &value[..4], &low4_all_ones);
    circ.zero_and_free(low4_all_ones);
}

fn coherent_compare_geq_modulus(
    circ: &mut Circuit,
    value: &[QReg],
    width: usize,
    modulus_le: &[u8],
    flag: &QReg,
) {
    if width == 256 && modulus_le == crate::point_add::trailmix_port::mod_arith::SECP256K1_P_LE {
        coherent_compare_geq_p_secp256k1(circ, value, flag);
    } else {
        crate::point_add::trailmix_port::arith::compare::compare_geq_const(
            circ, value, modulus_le, flag,
        );
    }
}

fn coherent_controlled_sub_modulus(
    circ: &mut Circuit,
    ctrl: &QReg,
    value: &[QReg],
    width: usize,
    modulus_le: &[u8],
) {
    if width == 256 && modulus_le == crate::point_add::trailmix_port::mod_arith::SECP256K1_P_LE {
        let mut r = [0u8; 32];
        r[0] = 0xd1;
        r[1] = 0x03;
        r[4] = 0x01;
        crate::point_add::trailmix_port::arith::const_add::controlled_add_const(
            circ, ctrl, value, &r,
        );
        circ.cx(ctrl, &value[256]);
    } else {
        crate::point_add::trailmix_port::arith::const_add::controlled_sub_const(
            circ,
            ctrl,
            value,
            modulus_le,
        );
    }
}

fn coherent_controlled_mod_add(
    circ: &mut Circuit,
    ctrl: &QReg,
    acc: &[QReg],
    addend: &[QReg],
    width: usize,
    modulus_le: &[u8],
) {
    crate::point_add::trailmix_port::arith::ripple_add::controlled_add(
        circ, ctrl, acc, addend,
    );

    let reduced = circ.alloc_qreg("coherent-cma.reduced");
    coherent_compare_geq_modulus(circ, acc, width, modulus_le, &reduced);
    coherent_controlled_sub_modulus(circ, &reduced, acc, width, modulus_le);

    // reduced = ctrl AND (acc_post < addend). Copy ctrl because the physical
    // comparator temporarily scrambles both compared registers; this is also
    // what makes exact squaring (ctrl aliases addend) legal.
    let ctrl_copy = circ.alloc_qreg("coherent-cma.ctrl-copy");
    circ.cx(ctrl, &ctrl_copy);
    let acc_ge_addend = circ.alloc_qreg("coherent-cma.acc-ge-addend");
    crate::point_add::trailmix_port::arith::compare::compare_geq_physical_middle(
        circ,
        &acc[..width],
        &addend[..width],
        &acc_ge_addend,
        |circ, geq| {
            circ.cx(&ctrl_copy, &reduced);
            circ.ccx(&ctrl_copy, geq, &reduced);
        },
    );
    circ.cx(ctrl, &ctrl_copy);

    circ.zero_and_free(acc_ge_addend);
    circ.zero_and_free(ctrl_copy);
    circ.zero_and_free(reduced);
}

fn coherent_mod_double(
    circ: &mut Circuit,
    acc: &[QReg],
    width: usize,
    modulus_le: &[u8],
) {
    crate::point_add::trailmix_port::arith::shift::left_shift(circ, acc);
    let reduced = circ.alloc_qreg("coherent-double.reduced");
    coherent_compare_geq_modulus(circ, acc, width, modulus_le, &reduced);
    coherent_controlled_sub_modulus(circ, &reduced, acc, width, modulus_le);

    // The unreduced doubled value is even and the modulus is odd, so the
    // canonical output parity is exactly the reduction predicate.
    circ.cx(&acc[0], &reduced);
    circ.zero_and_free(reduced);
}

fn emit_coherent_horner(
    circ: &mut Circuit,
    result: &[QReg],
    a: &[QReg],
    b: &[QReg],
    width: usize,
    modulus_le: &[u8],
) {
    coherent_controlled_mod_add(circ, &b[width - 1], result, a, width, modulus_le);
    for i in (0..width - 1).rev() {
        coherent_mod_double(circ, result, width, modulus_le);
        coherent_controlled_mod_add(circ, &b[i], result, a, width, modulus_le);
    }
}

/// Emit a canonical gate-only product and return its exact literal inverse.
///
/// `result`, `a`, and `b` have `width + 1` lanes, with the final lane clean.
/// `a` and `b` may be the same register (squaring), but partial aliasing and
/// result/operand aliasing are rejected.
#[doc(hidden)]
pub fn mod_mul_canonical_coherent_temp_with_modulus(
    circ: &mut Circuit,
    result: &[QReg],
    a: &[QReg],
    b: &[QReg],
    width: usize,
    modulus_le: &[u8],
) -> CanonicalCoherentMulInverse {
    assert_eq!(result.len(), width + 1);
    assert_eq!(a.len(), width + 1);
    assert_eq!(b.len(), width + 1);
    assert_modulus(width, modulus_le);
    assert_register_alias_rules(result, a, b);

    let prev = circ.push_section("mul_canonical_coherent_temp");
    let op_start = circ.b.ops.len();
    let counted_start = circ.b.counted_ops;
    let kind_start = circ.b.counted_kind_ops;
    let active_start = circ.b.active_qubits;
    let bits_start = circ.b.next_bit;
    let registers_start = circ.b.next_register;

    circ.begin_gate_only_block();
    emit_coherent_horner(circ, result, a, b, width, modulus_le);
    let workspace_ids = circ.finish_gate_only_block();

    assert_eq!(
        circ.b.active_qubits, active_start,
        "coherent multiplier retained workspace"
    );
    assert_eq!(circ.b.next_bit, bits_start, "coherent multiplier allocated a bit");
    assert_eq!(
        circ.b.next_register, registers_start,
        "coherent multiplier emitted register metadata"
    );
    assert!(
        workspace_ids
            .iter()
            .all(|id| circ.b.free_qubits.contains(id)),
        "coherent multiplier workspace is not clean-released"
    );

    let gate_count = circ.b.counted_ops - counted_start;
    let mut kind_counts = [0usize; OP_KIND_COUNT];
    for (kind, count) in kind_counts.iter_mut().enumerate() {
        *count = circ.b.counted_kind_ops[kind] - kind_start[kind];
    }
    assert_gate_alphabet(&kind_counts);
    if !circ.b.count_only {
        assert_eq!(circ.b.ops.len() - op_start, gate_count);
        assert!(
            circ.b.ops[op_start..].iter().all(is_literal_gate),
            "coherent multiplier emitted a conditioned or non-gate operation"
        );
    }
    circ.pop_section(&prev);

    CanonicalCoherentMulInverse {
        count_only: circ.b.count_only,
        op_start,
        gate_count,
        kind_counts,
        workspace_ids,
        result_ids: ids(result),
        a_ids: ids(a),
        b_ids: ids(b),
    }
}

/// secp256k1 specialization of
/// [`mod_mul_canonical_coherent_temp_with_modulus`].
pub fn mod_mul_canonical_coherent_temp(
    circ: &mut Circuit,
    result: &[QReg],
    a: &[QReg],
    b: &[QReg],
) -> CanonicalCoherentMulInverse {
    mod_mul_canonical_coherent_temp_with_modulus(
        circ,
        result,
        a,
        b,
        256,
        &crate::point_add::trailmix_port::mod_arith::SECP256K1_P_LE,
    )
}

/// Consume a matching token and append the exact forward `Op` records in
/// reverse order. Every allowed gate is self-inverse, so no semantic inverse
/// arithmetic, reset, measurement, or classical-bit operation is involved.
pub fn mod_mul_canonical_coherent_temp_reverse(
    circ: &mut Circuit,
    inverse: CanonicalCoherentMulInverse,
    result: &[QReg],
    a: &[QReg],
    b: &[QReg],
) {
    assert_eq!(ids(result), inverse.result_ids, "coherent inverse result mismatch");
    assert_eq!(ids(a), inverse.a_ids, "coherent inverse operand a mismatch");
    assert_eq!(ids(b), inverse.b_ids, "coherent inverse operand b mismatch");
    assert_eq!(
        circ.b.count_only, inverse.count_only,
        "coherent inverse builder mode mismatch"
    );

    let prev = circ.push_section("mul_canonical_coherent_temp_reverse");
    circ.flush_pending_frees();
    let active_start = circ.b.active_qubits;
    let free_pool = circ.b.free_qubits.clone();
    for &id in &inverse.workspace_ids {
        circ.b.reacquire(QubitId(u64::from(id)));
    }

    if circ.b.count_only {
        assert!(
            circ.b.clone_fiat_hash().is_none(),
            "count-only literal reversal cannot synthesize a Fiat-Shamir byte stream"
        );
        for kind in [
            OperationType::X,
            OperationType::Z,
            OperationType::CX,
            OperationType::CZ,
            OperationType::Swap,
            OperationType::CCX,
            OperationType::CCZ,
        ] {
            circ.b
                .add_counted_kind(kind, inverse.kind_counts[kind as usize]);
        }
    } else {
        let end = inverse.op_start + inverse.gate_count;
        assert!(end <= circ.b.ops.len(), "coherent forward gate slice is missing");
        for offset in 0..inverse.gate_count {
            let op = circ.b.ops[end - 1 - offset];
            assert!(is_literal_gate(&op));
            circ.b.push_op(op);
        }
    }

    for &id in &inverse.workspace_ids {
        circ.b.free_clean(QubitId(u64::from(id)));
    }
    circ.b.free_qubits = free_pool;
    assert_eq!(
        circ.b.active_qubits, active_start,
        "coherent inverse retained workspace"
    );
    circ.pop_section(&prev);
}
