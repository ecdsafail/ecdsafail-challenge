
use super::{BitId, QubitId, B};
use crate::circuit::{Op, OperationType};

#[allow(dead_code)]
pub(crate) fn xor_right_shifted_carries_into_classical(
    b: &mut B,
    q_src: &[QubitId],
    offset_bits: u64,
    q_dst: &[QubitId],
    carry_in: bool,
) {
    let n = q_dst.len();
    assert!(n <= q_src.len() && q_src.len() <= n + 1, "len mismatch");
    if n == 0 {
        return;
    }

    let bit = |k: usize| -> bool {
        if k >= 64 {
            false
        } else {
            (offset_bits >> k) & 1 != 0
        }
    };

    let ccx_inv =
        |b: &mut B, ctrl_a: QubitId, inv_a: bool, ctrl_b: QubitId, inv_b: bool, target: QubitId| {
            if inv_a {
                b.x(ctrl_a);
            }
            if inv_b {
                b.x(ctrl_b);
            }
            b.ccx(ctrl_a, ctrl_b, target);
            if inv_b {
                b.x(ctrl_b);
            }
            if inv_a {
                b.x(ctrl_a);
            }
        };

    for k in (1..n).rev() {
        ccx_inv(b, q_src[k], bit(k), q_dst[k - 1], false, q_dst[k]);
    }

    for k in 0..n {
        if bit(k) {
            b.x(q_dst[k]);
        }
    }

    let carry_in_xor_offset0 = carry_in ^ bit(0);
    if carry_in_xor_offset0 {

        if bit(0) {
            b.x(q_src[0]);
        }
        b.cx(q_src[0], q_dst[0]);
        if bit(0) {
            b.x(q_src[0]);
        }
    }

    for k in 1..n {
        ccx_inv(b, q_src[k], bit(k), q_dst[k - 1], bit(k), q_dst[k]);
    }
}

pub(crate) fn add_vented_2clean_classical(
    b: &mut B,
    q_target: &[QubitId],
    q_clean2: &[QubitId; 2],
    offset_bits: u64,
    carry_in: bool,
    vent_keys: &[BitId],
) {
    add_vented_2clean_classical_cxt(
        b,
        q_target,
        q_clean2,
        offset_bits,
        carry_in,
        vent_keys,
        None,
    );
}

pub(crate) fn add_vented_2clean_classical_cxt(
    b: &mut B,
    q_target: &[QubitId],
    q_clean2: &[QubitId; 2],
    offset_bits: u64,
    carry_in: bool,
    vent_keys: &[BitId],
    carry_xor_target: Option<&[Option<QubitId>]>,
) {
    let n = q_target.len();
    if n == 0 {
        return;
    }
    let bit = |k: usize| -> bool {
        if k >= 64 {
            false
        } else {
            (offset_bits >> k) & 1 != 0
        }
    };

    if n == 1 {
        if carry_in {
            b.x(q_target[0]);
        }
        if bit(0) {
            b.x(q_target[0]);
        }
        return;
    }

    for k in 0..n {
        if bit(k) {
            b.x(q_target[k]);
        }
    }

    let get_carry_qubit = |k: usize| -> Option<QubitId> {
        if k == 0 {
            None
        } else if k == n - 1 {
            Some(q_target[n - 1])
        } else {
            Some(q_clean2[k % 2])
        }
    };

    for k in 0..n - 1 {

        if k < n - 2 {
            if let Some(q) = get_carry_qubit(k + 1) {

                let mut op = Op::empty();
                op.kind = OperationType::R;
                op.q_target = q;
                b.ops.push(op);
            }
        }

        if k == 0 {
            let eff_carry = carry_in ^ bit(0);
            if eff_carry {

                if let Some(q) = get_carry_qubit(1) {
                    b.cx(q_target[0], q);
                }
            }
        } else {
            let carry_q = get_carry_qubit(k).expect("non-boundary carry");
            let carry_next = get_carry_qubit(k + 1).expect("non-boundary next carry");
            if bit(k) {
                b.x(carry_q);
                b.ccx(q_target[k], carry_q, carry_next);
                b.x(carry_q);
            } else {
                b.ccx(q_target[k], carry_q, carry_next);
            }
        }

        if k == 0 {
            if carry_in {
                b.x(q_target[0]);
            }
        } else {
            let carry_q = get_carry_qubit(k).expect("non-boundary carry");
            b.cx(carry_q, q_target[k]);
        }

        if let Some(cxt) = carry_xor_target {
            if k < cxt.len() {
                if let Some(dst) = cxt[k] {
                    if k == 0 {
                        if carry_in {
                            b.x(dst);
                        }
                    } else {
                        let carry_q = get_carry_qubit(k).expect("non-boundary carry");
                        b.cx(carry_q, dst);
                    }
                }
            }
        }

        if k > 0 {
            let carry_q = get_carry_qubit(k).expect("non-boundary carry");
            b.hmr(carry_q, vent_keys[k]);
        }

        if bit(k) {
            if let Some(q) = get_carry_qubit(k + 1) {
                b.x(q);
            }
        }
    }
}

pub(crate) fn iadd_linear_clean_classical(
    b: &mut B,
    q_target: &[QubitId],
    q_clean: &[QubitId],
    offset_bits: u64,
    carry_in: bool,
) {
    let n = q_target.len();
    if n == 0 {
        return;
    }
    assert!(q_clean.len() >= n.saturating_sub(2), "need n-2 clean");
    let q_clean = &q_clean[..n.saturating_sub(2)];

    let bit = |k: usize| -> bool {
        if k >= 64 {
            false
        } else {
            (offset_bits >> k) & 1 != 0
        }
    };

    if n == 1 {
        if bit(0) {
            b.x(q_target[0]);
        }
        if carry_in {
            b.x(q_target[0]);
        }
        return;
    }

    if n == 2 {

        if bit(0) {
            b.x(q_target[1]);
        }

        for k in 0..2 {
            if bit(k) {
                b.x(q_target[k]);
            }
        }

        let eff0 = carry_in ^ bit(0);
        if eff0 {
            b.cx(q_target[0], q_target[1]);
        }

        if carry_in {
            b.x(q_target[0]);
        }
        return;
    }

    for &q in q_clean.iter() {
        let mut op = Op::empty();
        op.kind = OperationType::R;
        op.q_target = q;
        b.ops.push(op);
    }

    let get_carry = |k: usize| -> Option<QubitId> {
        if k == 0 {
            None
        } else if k == n - 1 {
            Some(q_target[n - 1])
        } else {
            Some(q_clean[k - 1])
        }
    };

    for k in 0..n - 1 {
        if bit(k) {
            if let Some(q) = get_carry(k + 1) {
                b.x(q);
            }
        }
    }

    for k in 0..n {
        if bit(k) {
            b.x(q_target[k]);
        }
    }

    for k in 0..n - 1 {

        let next = get_carry(k + 1).expect("k+1 in bounds");
        if k == 0 {

            let eff = carry_in ^ bit(0);
            if eff {
                b.cx(q_target[0], next);
            }
        } else {
            let cur = get_carry(k).expect("k in bounds");
            if bit(k) {
                b.x(cur);
                b.ccx(q_target[k], cur, next);
                b.x(cur);
            } else {
                b.ccx(q_target[k], cur, next);
            }
        }
    }

    for k in (0..n - 2).rev() {

        let next = get_carry(k + 1).expect("k+1 in bounds");
        b.cx(next, q_target[k + 1]);

        let m = b.alloc_bit();
        b.hmr(next, m);

        if bit(k) {
            let mut op = Op::empty();
            op.kind = OperationType::Neg;
            op.c_condition = m;
            b.ops.push(op);
        }

        if k == 0 {

            let eff = carry_in ^ bit(0);
            if eff {

                let mut op = Op::empty();
                op.kind = OperationType::Z;
                op.q_target = q_target[k];
                op.c_condition = m;
                b.ops.push(op);
            }
        } else {
            let cur = get_carry(k).expect("k in bounds");

            if bit(k) {
                b.x(cur);
                b.cz_if(q_target[k], cur, m);
                b.x(cur);
            } else {
                b.cz_if(q_target[k], cur, m);
            }
        }
    }

    if carry_in {
        b.x(q_target[0]);
    }
}

#[allow(dead_code)]
pub(crate) fn iadd_dirty_2clean_classical(
    b: &mut B,
    q_target: &[QubitId],
    q_dirty: &[QubitId],
    q_clean2: &[QubitId; 2],
    offset_bits: u64,
    carry_in: bool,
) {
    let n = q_target.len();
    if n == 0 {
        return;
    }

    if n <= 4 {
        iadd_linear_clean_classical(b, q_target, q_clean2, offset_bits, carry_in);
        return;
    }
    assert!(q_dirty.len() >= n - 2, "need n-2 dirty qubits");
    let q_dirty = &q_dirty[..n - 2];

    let vent_keys: Vec<BitId> = (0..n).map(|_| b.alloc_bit()).collect();

    let cxt: Vec<Option<QubitId>> = (0..n)
        .map(|k| {
            if k == 0 {
                None
            } else {
                q_dirty.get(k - 1).copied()
            }
        })
        .collect();

    add_vented_2clean_classical_cxt(
        b,
        q_target,
        q_clean2,
        offset_bits,
        carry_in,
        &vent_keys,
        Some(&cxt),
    );

    for k in 0..n {
        b.x(q_target[k]);
    }

    for k in 0..n - 2 {
        let mut op = Op::empty();
        op.kind = OperationType::Z;
        op.q_target = q_dirty[k];
        op.c_condition = vent_keys[k + 1];
        b.ops.push(op);
    }

    xor_right_shifted_carries_into_classical(b, &q_target[..n - 1], offset_bits, q_dirty, carry_in);
    for k in 0..n - 2 {
        let mut op = Op::empty();
        op.kind = OperationType::Z;
        op.q_target = q_dirty[k];
        op.c_condition = vent_keys[k + 1];
        b.ops.push(op);
    }
    for k in 0..n {
        b.x(q_target[k]);
    }
}

pub(crate) fn ciadd_dirty_2clean_classical(
    b: &mut B,
    q_target: &[QubitId],
    q_dirty: &[QubitId],
    q_clean2: &[QubitId; 2],
    offset_bits: u64,
    ctrl: QubitId,
    carry_in: bool,
) {

    assert!(
        !carry_in,
        "ciadd_dirty_2clean_classical requires carry_in=false; pre-process if needed"
    );
    let n = q_target.len();
    if n == 0 {
        return;
    }
    if n <= 4 {

        let a: Vec<QubitId> = (0..n).map(|_| b.alloc_qubit()).collect();
        for i in 0..n {
            if (offset_bits >> i) & 1 != 0 {
                b.cx(ctrl, a[i]);
            }
        }

        for i in 0..n {
            if (offset_bits >> i) & 1 != 0 {
                b.cx(ctrl, a[i]);
            }
        }
        for q in a {
            b.free(q);
        }
        panic!("ciadd_dirty_2clean: n<=4 fallback not implemented; use uncontrolled path");
    }
    assert!(q_dirty.len() >= n - 2, "need n-2 dirty qubits");
    let q_dirty = &q_dirty[..n - 2];

    let vent_keys: Vec<BitId> = (0..n).map(|_| b.alloc_bit()).collect();

    let cxt: Vec<Option<QubitId>> = (0..n)
        .map(|k| {
            if k == 0 {
                None
            } else {
                q_dirty.get(k - 1).copied()
            }
        })
        .collect();

    c_add_vented_2clean_inline(
        b,
        q_target,
        q_clean2,
        offset_bits,
        ctrl,
        carry_in,
        &vent_keys,
        &cxt,
    );

    for k in 0..n {

        b.cx(ctrl, q_target[k]);
    }
    for k in 0..n - 2 {

        let mut op = Op::empty();
        op.kind = OperationType::Z;
        op.q_target = q_dirty[k];
        op.c_condition = vent_keys[k + 1];
        b.ops.push(op);
    }

    c_xor_right_shifted_carries_into_classical(
        b,
        &q_target[..n - 1],
        offset_bits,
        ctrl,
        q_dirty,
        carry_in,
    );
    for k in 0..n - 2 {
        let mut op = Op::empty();
        op.kind = OperationType::Z;
        op.q_target = q_dirty[k];
        op.c_condition = vent_keys[k + 1];
        b.ops.push(op);
    }
    for k in 0..n {
        b.cx(ctrl, q_target[k]);
    }
}

fn c_add_vented_2clean_inline(
    b: &mut B,
    q_target: &[QubitId],
    q_clean2: &[QubitId; 2],
    offset_bits: u64,
    ctrl: QubitId,
    carry_in: bool,
    vent_keys: &[BitId],
    carry_xor_target: &[Option<QubitId>],
) {
    let n = q_target.len();
    if n < 2 {

        if n == 1 {
            if carry_in {
                b.cx(ctrl, q_target[0]);
            }
            if (offset_bits & 1) != 0 {
                b.cx(ctrl, q_target[0]);
            }
        }
        return;
    }

    let bit = |k: usize| -> bool {
        if k >= 64 {
            false
        } else {
            (offset_bits >> k) & 1 != 0
        }
    };

    for k in 0..n {
        if bit(k) {
            b.cx(ctrl, q_target[k]);
        }
    }

    let get_carry_qubit = |k: usize| -> Option<QubitId> {
        if k == 0 {
            None
        } else if k == n - 1 {
            Some(q_target[n - 1])
        } else {
            Some(q_clean2[k % 2])
        }
    };

    for k in 0..n - 1 {

        if k < n - 2 {
            if let Some(q) = get_carry_qubit(k + 1) {
                let mut op = Op::empty();
                op.kind = OperationType::R;
                op.q_target = q;
                b.ops.push(op);
            }
        }

        if k == 0 {
            let next = get_carry_qubit(1);
            if let Some(next_q) = next {
                if bit(0) {

                    if carry_in {

                        b.x(ctrl);
                        b.ccx(q_target[0], ctrl, next_q);
                        b.x(ctrl);
                    } else {
                        b.ccx(q_target[0], ctrl, next_q);
                    }
                } else if carry_in {

                    b.cx(q_target[0], next_q);
                }

            }
        } else {
            let cur = get_carry_qubit(k).expect("non-boundary carry");
            let next = get_carry_qubit(k + 1).expect("non-boundary next carry");
            if bit(k) {

                b.cx(ctrl, cur);
                b.ccx(q_target[k], cur, next);
                b.cx(ctrl, cur);
            } else {
                b.ccx(q_target[k], cur, next);
            }
        }

        if k == 0 {
            if carry_in {
                b.x(q_target[0]);
            }
        } else {
            let cur = get_carry_qubit(k).expect("non-boundary carry");
            b.cx(cur, q_target[k]);
        }

        if k < carry_xor_target.len() {
            if let Some(dst) = carry_xor_target[k] {
                if k == 0 {
                    if carry_in {
                        b.x(dst);
                    }
                } else {
                    let cur = get_carry_qubit(k).expect("non-boundary carry");
                    b.cx(cur, dst);
                }
            }
        }

        if k > 0 {
            let cur = get_carry_qubit(k).expect("non-boundary carry");
            b.hmr(cur, vent_keys[k]);
        }

        if bit(k) {
            if let Some(q) = get_carry_qubit(k + 1) {
                b.cx(ctrl, q);
            }
        }
    }
}

pub(crate) fn add_vented_2clean_qoffset(
    b: &mut B,
    q_target: &[QubitId],
    q_clean2: &[QubitId; 2],
    q_offset: &[QubitId],
    carry_in: bool,
    vent_keys: &[BitId],
    carry_xor_target: Option<&[Option<QubitId>]>,
) {
    let n = q_target.len();
    assert_eq!(q_offset.len(), n, "q_offset length must match q_target");
    if n == 0 {
        return;
    }
    if n == 1 {
        if carry_in {
            b.x(q_target[0]);
        }
        b.cx(q_offset[0], q_target[0]);
        return;
    }

    for k in 0..n {
        b.cx(q_offset[k], q_target[k]);
    }

    let get_carry_qubit = |k: usize| -> Option<QubitId> {
        if k == 0 {
            None
        } else if k == n - 1 {
            Some(q_target[n - 1])
        } else {
            Some(q_clean2[k % 2])
        }
    };

    for k in 0..n - 1 {
        if k < n - 2 {
            if let Some(q) = get_carry_qubit(k + 1) {
                let mut op = Op::empty();
                op.kind = OperationType::R;
                op.q_target = q;
                b.ops.push(op);
            }
        }

        if k == 0 {
            let next = get_carry_qubit(1);
            if let Some(next_q) = next {
                if carry_in {
                    b.x(q_offset[0]);
                    b.ccx(q_target[0], q_offset[0], next_q);
                    b.x(q_offset[0]);
                } else {
                    b.ccx(q_target[0], q_offset[0], next_q);
                }
            }
        } else {
            let cur = get_carry_qubit(k).expect("non-boundary carry");
            let next = get_carry_qubit(k + 1).expect("non-boundary next carry");

            b.cx(q_offset[k], cur);
            b.ccx(q_target[k], cur, next);
            b.cx(q_offset[k], cur);
        }

        if k == 0 {
            if carry_in {
                b.x(q_target[0]);
            }
        } else {
            let cur = get_carry_qubit(k).expect("non-boundary carry");
            b.cx(cur, q_target[k]);
        }

        if let Some(cxt) = carry_xor_target {
            if k < cxt.len() {
                if let Some(dst) = cxt[k] {
                    if k == 0 {
                        if carry_in {
                            b.x(dst);
                        }
                    } else {
                        let cur = get_carry_qubit(k).expect("non-boundary carry");
                        b.cx(cur, dst);
                    }
                }
            }
        }

        if k > 0 {
            let cur = get_carry_qubit(k).expect("non-boundary carry");
            b.hmr(cur, vent_keys[k]);
        }

        if let Some(q) = get_carry_qubit(k + 1) {
            b.cx(q_offset[k], q);
        }
    }
}

pub(crate) fn xor_right_shifted_carries_into_qoffset(
    b: &mut B,
    q_src: &[QubitId],
    q_offset: &[QubitId],
    q_dst: &[QubitId],
    carry_in: bool,
) {
    let n = q_dst.len();
    assert!(n <= q_src.len() && q_src.len() <= n + 1, "len mismatch");
    if n == 0 {
        return;
    }

    let ccx_with_qxor = |b: &mut B,
                         ctrl_a: QubitId,
                         xor_a: Option<QubitId>,
                         ctrl_b: QubitId,
                         xor_b: Option<QubitId>,
                         target: QubitId| {
        if let Some(x) = xor_a {
            b.cx(x, ctrl_a);
        }
        if let Some(x) = xor_b {
            b.cx(x, ctrl_b);
        }
        b.ccx(ctrl_a, ctrl_b, target);
        if let Some(x) = xor_b {
            b.cx(x, ctrl_b);
        }
        if let Some(x) = xor_a {
            b.cx(x, ctrl_a);
        }
    };

    for k in (1..n).rev() {
        ccx_with_qxor(b, q_src[k], Some(q_offset[k]), q_dst[k - 1], None, q_dst[k]);
    }

    for k in 0..n {
        b.cx(q_offset[k], q_dst[k]);
    }

    b.cx(q_offset[0], q_src[0]);
    if carry_in {
        b.x(q_offset[0]);
    }
    b.ccx(q_src[0], q_offset[0], q_dst[0]);
    if carry_in {
        b.x(q_offset[0]);
    }
    b.cx(q_offset[0], q_src[0]);

    for k in 1..n {
        ccx_with_qxor(
            b,
            q_src[k],
            Some(q_offset[k]),
            q_dst[k - 1],
            Some(q_offset[k]),
            q_dst[k],
        );
    }
}

pub(crate) fn iadd_dirty_2clean_qoffset(
    b: &mut B,
    q_target: &[QubitId],
    q_dirty: &[QubitId],
    q_clean2: &[QubitId; 2],
    q_offset: &[QubitId],
    carry_in: bool,
) {
    let n = q_target.len();
    assert_eq!(q_offset.len(), n);
    if n == 0 {
        return;
    }
    if n <= 4 {
        panic!("iadd_dirty_2clean_qoffset: n<=4 not supported yet, use cuccaro_add");
    }
    assert!(q_dirty.len() >= n - 2, "need n-2 dirty qubits");
    let q_dirty = &q_dirty[..n - 2];

    let vent_keys: Vec<BitId> = (0..n).map(|_| b.alloc_bit()).collect();
    let cxt: Vec<Option<QubitId>> = (0..n)
        .map(|k| {
            if k == 0 {
                None
            } else {
                q_dirty.get(k - 1).copied()
            }
        })
        .collect();

    add_vented_2clean_qoffset(
        b,
        q_target,
        q_clean2,
        q_offset,
        carry_in,
        &vent_keys,
        Some(&cxt),
    );

    for k in 0..n {
        b.x(q_target[k]);
    }
    for k in 0..n - 2 {
        let mut op = Op::empty();
        op.kind = OperationType::Z;
        op.q_target = q_dirty[k];
        op.c_condition = vent_keys[k + 1];
        b.ops.push(op);
    }
    xor_right_shifted_carries_into_qoffset(b, &q_target[..n - 1], q_offset, q_dirty, carry_in);
    for k in 0..n - 2 {
        let mut op = Op::empty();
        op.kind = OperationType::Z;
        op.q_target = q_dirty[k];
        op.c_condition = vent_keys[k + 1];
        b.ops.push(op);
    }
    for k in 0..n {
        b.x(q_target[k]);
    }
}

pub(crate) fn isub_dirty_2clean_qoffset(
    b: &mut B,
    q_target: &[QubitId],
    q_dirty: &[QubitId],
    q_clean2: &[QubitId; 2],
    q_offset: &[QubitId],
) {
    let n = q_target.len();
    for k in 0..n {
        b.x(q_target[k]);
    }
    iadd_dirty_2clean_qoffset(b, q_target, q_dirty, q_clean2, q_offset, false);
    for k in 0..n {
        b.x(q_target[k]);
    }
}

fn c_xor_right_shifted_carries_into_classical(
    b: &mut B,
    q_src: &[QubitId],
    offset_bits: u64,
    ctrl: QubitId,
    q_dst: &[QubitId],
    carry_in: bool,
) {
    let n = q_dst.len();
    assert!(n <= q_src.len() && q_src.len() <= n + 1);
    if n == 0 {
        return;
    }
    let bit = |k: usize| -> bool {
        if k >= 64 {
            false
        } else {
            (offset_bits >> k) & 1 != 0
        }
    };

    let ccx_ctrl_mix = |b: &mut B,
                        ctrl_a: QubitId,
                        a_xor_ctrl: bool,
                        ctrl_b: QubitId,
                        b_xor_ctrl: bool,
                        target: QubitId| {
        if a_xor_ctrl {
            b.cx(ctrl, ctrl_a);
        }
        if b_xor_ctrl {
            b.cx(ctrl, ctrl_b);
        }
        b.ccx(ctrl_a, ctrl_b, target);
        if b_xor_ctrl {
            b.cx(ctrl, ctrl_b);
        }
        if a_xor_ctrl {
            b.cx(ctrl, ctrl_a);
        }
    };

    for k in (1..n).rev() {
        ccx_ctrl_mix(b, q_src[k], bit(k), q_dst[k - 1], false, q_dst[k]);
    }

    for k in 0..n {
        if bit(k) {
            b.cx(ctrl, q_dst[k]);
        }
    }

    let cin_eff_uses_ctrl = bit(0);
    let cin_classical_part = carry_in ^ false;
    if cin_eff_uses_ctrl {

        if carry_in {

            b.cx(ctrl, q_src[0]);
            b.x(ctrl);
            b.ccx(q_src[0], ctrl, q_dst[0]);
            b.x(ctrl);
            b.cx(ctrl, q_src[0]);
        } else {

            b.cx(ctrl, q_src[0]);
            b.ccx(q_src[0], ctrl, q_dst[0]);
            b.cx(ctrl, q_src[0]);
        }
    } else {

        if cin_classical_part {

            b.cx(q_src[0], q_dst[0]);
        }

    }
    for k in 1..n {
        ccx_ctrl_mix(b, q_src[k], bit(k), q_dst[k - 1], bit(k), q_dst[k]);
    }
}

pub(crate) fn cisub_dirty_2clean_classical(
    b: &mut B,
    q_target: &[QubitId],
    q_dirty: &[QubitId],
    q_clean2: &[QubitId; 2],
    c_bits: u64,
    ctrl: QubitId,
) {
    let n = q_target.len();

    for k in 0..n {
        b.cx(ctrl, q_target[k]);
    }
    ciadd_dirty_2clean_classical(
        b, q_target, q_dirty, q_clean2, c_bits, ctrl,
        false,
    );
    for k in 0..n {
        b.cx(ctrl, q_target[k]);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sim::Simulator;
    use sha3::{
        digest::{ExtendableOutput, Update},
        Shake256,
    };

    fn anf_degree_density_from_truth_table(mut table: Vec<u8>, vars: usize) -> (usize, usize) {
        let states = 1usize << vars;

        for bit in 0..vars {
            let step = 1usize << bit;
            for mask in 0..states {
                if (mask & step) != 0 {
                    table[mask] ^= table[mask ^ step];
                }
            }
        }

        let mut degree = 0usize;
        let mut density = 0usize;
        for (mask, &coeff) in table.iter().enumerate() {
            if coeff != 0 {
                density += 1;
                degree = degree.max(mask.count_ones() as usize);
            }
        }
        (degree, density)
    }

    fn product_phase_anf_degree_density(n: usize, phase_mask: u64) -> (usize, usize) {
        assert!(n > 0 && n <= 10, "test keeps exhaustive table small");
        let vars = 2 * n;
        let states = 1usize << vars;
        let x_mask = (1u64 << n) - 1;
        let mut table = vec![0u8; states];
        for state in 0..states {
            let x = (state as u64) & x_mask;
            let y = ((state as u64) >> n) & x_mask;
            let prod = x * y;
            table[state] = ((prod & phase_mask).count_ones() & 1) as u8;
        }
        anf_degree_density_from_truth_table(table, vars)
    }

    fn carry_save_product_bits_for_phase_test(n: usize, x: u64, y: u64) -> Vec<u8> {

        let mut cols = vec![Vec::<u8>::new(); 2 * n + 8];
        for i in 0..n {
            for j in 0..n {
                let bit = (((x >> i) & 1) & ((y >> j) & 1)) as u8;
                cols[i + j].push(bit);
            }
        }
        for k in 0..cols.len() - 1 {
            while cols[k].len() > 2 {
                let a = cols[k].pop().unwrap();
                let b = cols[k].pop().unwrap();
                let c = cols[k].pop().unwrap();
                let sum = a ^ b ^ c;
                let carry = (a & b) ^ (a & c) ^ (b & c);
                cols[k].push(sum);
                cols[k + 1].push(carry);
            }
        }
        let mut out = Vec::with_capacity(4 * n + 4);
        for col in cols.iter().take(2 * n + 2) {
            out.push(*col.get(0).unwrap_or(&0));
            out.push(*col.get(1).unwrap_or(&0));
        }
        out
    }

    fn carry_save_product_phase_anf_degree_density(
        n: usize,
        top_column_only: bool,
    ) -> (usize, usize) {
        assert!(
            n > 0 && n <= 8,
            "test keeps exhaustive carry-save table small"
        );
        let vars = 2 * n;
        let states = 1usize << vars;
        let x_mask = (1u64 << n) - 1;
        let mut table = vec![0u8; states];
        for state in 0..states {
            let x = (state as u64) & x_mask;
            let y = ((state as u64) >> n) & x_mask;
            let bits = carry_save_product_bits_for_phase_test(n, x, y);
            table[state] = if top_column_only {
                let k = 2 * (2 * n - 2);
                bits[k] ^ bits[k + 1]
            } else {
                bits.iter().fold(0u8, |acc, &b| acc ^ b)
            };
        }
        anf_degree_density_from_truth_table(table, vars)
    }

    #[test]
    fn raw_product_measurement_phase_is_dense_not_free_kickmix() {

        for &n in &[4usize, 6, 8, 10] {
            let full_mask = if 2 * n == 64 {
                u64::MAX
            } else {
                (1u64 << (2 * n)) - 1
            };
            let high_mask = 1u64 << (2 * n - 2);
            let (deg_full, dens_full) = product_phase_anf_degree_density(n, full_mask);
            let (deg_high, dens_high) = product_phase_anf_degree_density(n, high_mask);
            eprintln!(
                "raw_product_phase n={n} full_deg={deg_full} full_density={dens_full} high_deg={deg_high} high_density={dens_high}"
            );
            if n == 10 {
                println!("METRIC raw_product_mbu_fullmask_degree_n10={deg_full}");
                println!("METRIC raw_product_mbu_fullmask_density_n10={dens_full}");
                println!("METRIC raw_product_mbu_highbit_degree_n10={deg_high}");
                println!("METRIC raw_product_mbu_highbit_density_n10={dens_high}");
            }
        }

        let (deg_full, dens_full) = product_phase_anf_degree_density(10, (1u64 << 20) - 1);
        let (deg_high, dens_high) = product_phase_anf_degree_density(10, 1u64 << 18);
        assert_eq!(deg_full, 19);
        assert_eq!(dens_full, 427_812);
        assert_eq!(deg_high, 19);
        assert_eq!(dens_high, 120_581);
    }

    #[test]
    fn carry_save_product_scratch_mbu_still_has_dense_phases() {

        for &n in &[4usize, 6, 8] {
            let (deg_all, dens_all) = carry_save_product_phase_anf_degree_density(n, false);
            let (deg_top, dens_top) = carry_save_product_phase_anf_degree_density(n, true);
            eprintln!(
                "carry_save_product_phase n={n} all_deg={deg_all} all_density={dens_all} top_deg={deg_top} top_density={dens_top}"
            );
            if n == 8 {
                println!("METRIC carry_save_product_mbu_all_degree_n8={deg_all}");
                println!("METRIC carry_save_product_mbu_all_density_n8={dens_all}");
                println!("METRIC carry_save_product_mbu_top_degree_n8={deg_top}");
                println!("METRIC carry_save_product_mbu_top_density_n8={dens_top}");
            }
        }
        let (deg_all, dens_all) = carry_save_product_phase_anf_degree_density(8, false);
        let (deg_top, dens_top) = carry_save_product_phase_anf_degree_density(8, true);
        assert_eq!(deg_all, 16);
        assert_eq!(dens_all, 20_440);
        assert_eq!(deg_top, 15);
        assert_eq!(dens_top, 3_602);
    }

    fn classical_carry(x: u64, d: u64, cin: bool, n: usize) -> u64 {

        let mut c: u64 = 0;
        let mut prev = cin;
        for k in 0..n {
            let xk = (x >> k) & 1 != 0;
            let dk = (d >> k) & 1 != 0;

            let new_carry = (prev && xk) || (prev && dk) || (xk && dk);
            if new_carry {
                c |= 1 << (k + 1);
            }
            prev = new_carry;
        }

        if cin {
            c |= 1;
        }
        c
    }

    fn run_xor_rsh_carries(n: usize, trials: usize) -> bool {
        let mut hasher = Shake256::default();
        hasher.update(&[n as u8, trials as u8, 42]);
        use sha3::digest::XofReader;
        let mut xof = <sha3::Shake256 as sha3::digest::ExtendableOutput>::finalize_xof(hasher);
        for _trial in 0..trials {
            let mut buf = [0u8; 32];
            xof.read(&mut buf);
            let src_raw = u64::from_le_bytes(buf[0..8].try_into().unwrap());
            let dst_raw = u64::from_le_bytes(buf[8..16].try_into().unwrap());
            let offset_raw = u64::from_le_bytes(buf[16..24].try_into().unwrap());
            let cin_raw = buf[24];
            let src = if n < 64 {
                src_raw & ((1u64 << n) - 1)
            } else {
                src_raw
            };
            let dst = if n < 64 {
                dst_raw & ((1u64 << n) - 1)
            } else {
                dst_raw
            };
            let offset = if n < 64 {
                offset_raw & ((1u64 << n) - 1)
            } else {
                offset_raw
            };
            let cin = (cin_raw & 1) != 0;

            let mut bb = B::new();
            let q_src: Vec<QubitId> = bb.alloc_qubits(n);
            let q_dst: Vec<QubitId> = bb.alloc_qubits(n);

            xor_right_shifted_carries_into_classical(&mut bb, &q_src, offset, &q_dst, cin);

            let ops = bb.ops.clone();
            let num_qubits = bb.next_qubit as usize;
            let num_bits = 0usize;
            let mut inner_hasher = Shake256::default();
            inner_hasher.update(&[77u8]);
            let mut inner_xof =
                <sha3::Shake256 as sha3::digest::ExtendableOutput>::finalize_xof(inner_hasher);
            let mut sim = Simulator::new(num_qubits, num_bits, &mut inner_xof);
            sim.clear_for_shot();

            for k in 0..n {
                if (src >> k) & 1 != 0 {
                    *sim.qubit_mut(q_src[k]) = 1;
                }
                if (dst >> k) & 1 != 0 {
                    *sim.qubit_mut(q_dst[k]) = 1;
                }
            }
            sim.apply(&ops);

            let expected_carries = classical_carry(src, offset, cin, n + 1);
            let expected_rsh = expected_carries >> 1;
            let expected_dst = (dst ^ expected_rsh) & ((1u64 << n) - 1);

            let mut got_dst: u64 = 0;
            for k in 0..n {
                if sim.qubit(q_dst[k]) & 1 != 0 {
                    got_dst |= 1 << k;
                }
            }
            if got_dst != expected_dst {
                eprintln!(
                    "n={} src={:#x} dst={:#x} offset={:#x} cin={} got={:#x} exp={:#x}",
                    n, src, dst, offset, cin, got_dst, expected_dst
                );
                return false;
            }
        }
        true
    }

    #[test]
    fn test_xor_rsh_carries_small() {
        for n in 1..=8 {
            assert!(run_xor_rsh_carries(n, 20), "failed at n={n}");
        }
    }

    fn run_vented_add_2clean(n: usize, trials: usize) -> (usize, usize) {
        let mut hasher = Shake256::default();
        hasher.update(&[n as u8, trials as u8, 51]);
        use sha3::digest::XofReader;
        let mut xof = <sha3::Shake256 as sha3::digest::ExtendableOutput>::finalize_xof(hasher);
        let mut ok = 0;
        let mut bad = 0;
        for _trial in 0..trials {
            let mut buf = [0u8; 24];
            xof.read(&mut buf);
            let target_raw = u64::from_le_bytes(buf[0..8].try_into().unwrap());
            let offset_raw = u64::from_le_bytes(buf[8..16].try_into().unwrap());
            let cin_raw = buf[16];
            let mask = if n < 64 { (1u64 << n) - 1 } else { u64::MAX };
            let target = target_raw & mask;
            let offset = offset_raw & mask;
            let cin = (cin_raw & 1) != 0;

            let mut bb = B::new();
            let q_target: Vec<QubitId> = bb.alloc_qubits(n);
            let q_clean2: [QubitId; 2] = [bb.alloc_qubit(), bb.alloc_qubit()];
            let vent_keys: Vec<BitId> = (0..n).map(|_| bb.alloc_bit()).collect();

            add_vented_2clean_classical(&mut bb, &q_target, &q_clean2, offset, cin, &vent_keys);

            let ops = bb.ops.clone();
            let num_qubits = bb.next_qubit as usize;
            let num_bits = bb.next_bit as usize;
            let mut inner_hasher = Shake256::default();
            inner_hasher.update(&[101u8]);
            let mut inner_xof =
                <sha3::Shake256 as sha3::digest::ExtendableOutput>::finalize_xof(inner_hasher);
            let mut sim = Simulator::new(num_qubits, num_bits, &mut inner_xof);
            sim.clear_for_shot();
            for k in 0..n {
                if (target >> k) & 1 != 0 {
                    *sim.qubit_mut(q_target[k]) = 1;
                }
            }
            sim.apply(&ops);

            let expected_sum = (target.wrapping_add(offset).wrapping_add(cin as u64)) & mask;
            let mut got: u64 = 0;
            for k in 0..n {
                if sim.qubit(q_target[k]) & 1 != 0 {
                    got |= 1 << k;
                }
            }
            if got == expected_sum {
                ok += 1;
            } else {
                bad += 1;
                if bad < 3 {
                    eprintln!(
                        "vented add FAIL n={} t={:#x} o={:#x} cin={} got={:#x} exp={:#x}",
                        n, target, offset, cin, got, expected_sum
                    );
                }
            }
        }
        (ok, bad)
    }

    #[test]
    fn test_vented_add_2clean_small() {
        for n in 2..=8 {
            let (ok, bad) = run_vented_add_2clean(n, 20);
            assert_eq!(bad, 0, "n={n}: {ok}/{} passed", ok + bad);
        }
    }

    fn run_linear_clean_add(n: usize, trials: usize) -> (usize, usize) {
        let mut hasher = Shake256::default();
        hasher.update(&[n as u8, trials as u8, 73]);
        use sha3::digest::XofReader;
        let mut xof = <sha3::Shake256 as sha3::digest::ExtendableOutput>::finalize_xof(hasher);
        let mut ok = 0;
        let mut bad = 0;
        for _trial in 0..trials {
            let mut buf = [0u8; 24];
            xof.read(&mut buf);
            let target_raw = u64::from_le_bytes(buf[0..8].try_into().unwrap());
            let offset_raw = u64::from_le_bytes(buf[8..16].try_into().unwrap());
            let cin_raw = buf[16];
            let mask = if n < 64 { (1u64 << n) - 1 } else { u64::MAX };
            let target = target_raw & mask;
            let offset = offset_raw & mask;
            let cin = (cin_raw & 1) != 0;

            let mut bb = B::new();
            let q_target: Vec<QubitId> = bb.alloc_qubits(n);
            let n_clean = n.saturating_sub(2).max(2);
            let q_clean: Vec<QubitId> = bb.alloc_qubits(n_clean);

            iadd_linear_clean_classical(&mut bb, &q_target, &q_clean, offset, cin);

            let ops = bb.ops.clone();
            let num_qubits = bb.next_qubit as usize;
            let num_bits = bb.next_bit as usize;
            let mut inner_hasher = Shake256::default();
            inner_hasher.update(&[151u8]);
            let mut inner_xof =
                <sha3::Shake256 as sha3::digest::ExtendableOutput>::finalize_xof(inner_hasher);
            let mut sim = Simulator::new(num_qubits, num_bits, &mut inner_xof);
            sim.clear_for_shot();
            for k in 0..n {
                if (target >> k) & 1 != 0 {
                    *sim.qubit_mut(q_target[k]) = 1;
                }
            }
            sim.apply(&ops);

            let expected_sum = (target.wrapping_add(offset).wrapping_add(cin as u64)) & mask;
            let mut got: u64 = 0;
            for k in 0..n {
                if sim.qubit(q_target[k]) & 1 != 0 {
                    got |= 1 << k;
                }
            }
            if got == expected_sum {
                ok += 1;
            } else {
                bad += 1;
                if bad < 3 {
                    eprintln!(
                        "HRS FAIL n={} t={:#x} o={:#x} cin={} got={:#x} exp={:#x}",
                        n, target, offset, cin, got, expected_sum
                    );
                }
            }
        }
        (ok, bad)
    }

    #[test]
    fn test_iadd_linear_clean_small() {
        for n in 1..=8 {
            let (ok, bad) = run_linear_clean_add(n, 20);
            assert_eq!(bad, 0, "n={n}: {ok}/{} passed", ok + bad);
        }
    }

    fn run_iadd_dirty_2clean(n: usize, trials: usize) -> (usize, usize) {
        let mut hasher = Shake256::default();
        hasher.update(&[n as u8, trials as u8, 97]);
        use sha3::digest::XofReader;
        let mut xof = <sha3::Shake256 as sha3::digest::ExtendableOutput>::finalize_xof(hasher);
        let mut ok = 0;
        let mut bad = 0;
        for _trial in 0..trials {
            let mut buf = [0u8; 32];
            xof.read(&mut buf);
            let target_raw = u64::from_le_bytes(buf[0..8].try_into().unwrap());
            let offset_raw = u64::from_le_bytes(buf[8..16].try_into().unwrap());
            let dirty_raw = u64::from_le_bytes(buf[16..24].try_into().unwrap());
            let cin_raw = buf[24];
            let mask = if n < 64 { (1u64 << n) - 1 } else { u64::MAX };
            let target = target_raw & mask;
            let offset = offset_raw & mask;
            let dirty_init = dirty_raw & mask;
            let cin = (cin_raw & 1) != 0;

            let mut bb = B::new();
            let q_target: Vec<QubitId> = bb.alloc_qubits(n);
            let q_dirty: Vec<QubitId> = bb.alloc_qubits(n.saturating_sub(2).max(1));
            let q_clean2: [QubitId; 2] = [bb.alloc_qubit(), bb.alloc_qubit()];

            iadd_dirty_2clean_classical(&mut bb, &q_target, &q_dirty, &q_clean2, offset, cin);

            let ops = bb.ops.clone();
            let num_qubits = bb.next_qubit as usize;
            let num_bits = bb.next_bit as usize;
            let mut inner_hasher = Shake256::default();
            inner_hasher.update(&[201u8]);
            let mut inner_xof =
                <sha3::Shake256 as sha3::digest::ExtendableOutput>::finalize_xof(inner_hasher);
            let mut sim = Simulator::new(num_qubits, num_bits, &mut inner_xof);
            sim.clear_for_shot();
            for k in 0..n {
                if (target >> k) & 1 != 0 {
                    *sim.qubit_mut(q_target[k]) = 1;
                }
            }

            for (k, &q) in q_dirty.iter().enumerate() {
                if (dirty_init >> k) & 1 != 0 {
                    *sim.qubit_mut(q) = 1;
                }
            }
            sim.apply(&ops);

            let expected_sum = (target.wrapping_add(offset).wrapping_add(cin as u64)) & mask;
            let mut got: u64 = 0;
            for k in 0..n {
                if sim.qubit(q_target[k]) & 1 != 0 {
                    got |= 1 << k;
                }
            }

            let mut got_dirty: u64 = 0;
            for (k, &q) in q_dirty.iter().enumerate() {
                if sim.qubit(q) & 1 != 0 {
                    got_dirty |= 1 << k;
                }
            }
            let dirty_ok = if n > 4 {
                got_dirty == (dirty_init & ((1u64 << q_dirty.len()) - 1).min(mask))
            } else {
                true
            };

            let phase = sim.global_phase() & 1;

            if got == expected_sum && dirty_ok && phase == 0 {
                ok += 1;
            } else {
                bad += 1;
                if bad < 3 {
                    eprintln!(
                        "iadd_dirty_2clean FAIL n={} t={:#x} o={:#x} d={:#x} cin={} got={:#x} exp={:#x} dirty_ok={} phase={}",
                        n, target, offset, dirty_init, cin, got, expected_sum, dirty_ok, phase
                    );
                }
            }
        }
        (ok, bad)
    }

    #[test]
    fn test_iadd_dirty_2clean_small() {
        for n in 2..=8 {
            let (ok, bad) = run_iadd_dirty_2clean(n, 10);
            assert_eq!(bad, 0, "n={n}: {ok}/{} passed", ok + bad);
        }
    }

    fn run_ciadd_dirty_2clean(n: usize, trials: usize) -> (usize, usize) {
        let mut hasher = Shake256::default();
        hasher.update(&[n as u8, trials as u8, 113]);
        use sha3::digest::XofReader;
        let mut xof = <sha3::Shake256 as sha3::digest::ExtendableOutput>::finalize_xof(hasher);
        let mut ok = 0;
        let mut bad = 0;
        for _trial in 0..trials {
            let mut buf = [0u8; 40];
            xof.read(&mut buf);
            let target_raw = u64::from_le_bytes(buf[0..8].try_into().unwrap());
            let offset_raw = u64::from_le_bytes(buf[8..16].try_into().unwrap());
            let dirty_raw = u64::from_le_bytes(buf[16..24].try_into().unwrap());
            let cin_raw = buf[24];
            let ctrl_raw = buf[25];
            let mask = if n < 64 { (1u64 << n) - 1 } else { u64::MAX };
            let target = target_raw & mask;
            let offset = offset_raw & mask;
            let dirty_init = dirty_raw & mask;
            let cin = false;
            let _ = cin_raw;
            let ctrl_val = (ctrl_raw & 1) != 0;

            let mut bb = B::new();
            let q_target: Vec<QubitId> = bb.alloc_qubits(n);
            let q_dirty: Vec<QubitId> = bb.alloc_qubits(n.saturating_sub(2).max(1));
            let q_clean2: [QubitId; 2] = [bb.alloc_qubit(), bb.alloc_qubit()];
            let q_ctrl = bb.alloc_qubit();

            ciadd_dirty_2clean_classical(
                &mut bb, &q_target, &q_dirty, &q_clean2, offset, q_ctrl, cin,
            );

            let ops = bb.ops.clone();
            let num_qubits = bb.next_qubit as usize;
            let num_bits = bb.next_bit as usize;
            let mut inner_hasher = Shake256::default();
            inner_hasher.update(&[211u8]);
            let mut inner_xof =
                <sha3::Shake256 as sha3::digest::ExtendableOutput>::finalize_xof(inner_hasher);
            let mut sim = Simulator::new(num_qubits, num_bits, &mut inner_xof);
            sim.clear_for_shot();
            for k in 0..n {
                if (target >> k) & 1 != 0 {
                    *sim.qubit_mut(q_target[k]) = 1;
                }
            }
            for (k, &q) in q_dirty.iter().enumerate() {
                if (dirty_init >> k) & 1 != 0 {
                    *sim.qubit_mut(q) = 1;
                }
            }
            if ctrl_val {
                *sim.qubit_mut(q_ctrl) = 1;
            }
            sim.apply(&ops);

            let expected_sum = if ctrl_val {
                (target.wrapping_add(offset).wrapping_add(cin as u64)) & mask
            } else {
                target
            };
            let mut got: u64 = 0;
            for k in 0..n {
                if sim.qubit(q_target[k]) & 1 != 0 {
                    got |= 1 << k;
                }
            }
            let mut got_dirty: u64 = 0;
            for (k, &q) in q_dirty.iter().enumerate() {
                if sim.qubit(q) & 1 != 0 {
                    got_dirty |= 1 << k;
                }
            }
            let dirty_ok = got_dirty == (dirty_init & ((1u64 << q_dirty.len()) - 1).min(mask));
            let phase = sim.global_phase() & 1;
            let ctrl_preserved = sim.qubit(q_ctrl) & 1 == (ctrl_val as u64);

            if got == expected_sum && dirty_ok && phase == 0 && ctrl_preserved {
                ok += 1;
            } else {
                bad += 1;
                if bad < 3 {
                    eprintln!(
                        "ciadd_dirty FAIL n={} t={:#x} o={:#x} d={:#x} cin={} ctrl={} got={:#x} exp={:#x} d_ok={} phase={} ctrl_preserved={}",
                        n, target, offset, dirty_init, cin, ctrl_val, got, expected_sum, dirty_ok, phase, ctrl_preserved
                    );
                }
            }
        }
        (ok, bad)
    }

    #[test]
    fn test_ciadd_dirty_2clean_small() {
        for n in 5..=10 {
            let (ok, bad) = run_ciadd_dirty_2clean(n, 10);
            assert_eq!(bad, 0, "n={n}: {ok}/{} passed", ok + bad);
        }
    }

    fn run_cisub_dirty(n: usize, trials: usize) -> (usize, usize) {
        let mut hasher = Shake256::default();
        hasher.update(&[n as u8, trials as u8, 179]);
        use sha3::digest::XofReader;
        let mut xof = <sha3::Shake256 as sha3::digest::ExtendableOutput>::finalize_xof(hasher);
        let mut ok = 0;
        let mut bad = 0;
        for _trial in 0..trials {
            let mut buf = [0u8; 40];
            xof.read(&mut buf);
            let target_raw = u64::from_le_bytes(buf[0..8].try_into().unwrap());
            let c_raw = u64::from_le_bytes(buf[8..16].try_into().unwrap());
            let dirty_raw = u64::from_le_bytes(buf[16..24].try_into().unwrap());
            let ctrl_raw = buf[25];
            let mask = if n < 64 { (1u64 << n) - 1 } else { u64::MAX };
            let target = target_raw & mask;
            let c = c_raw & mask;
            let dirty_init = dirty_raw & mask;
            let ctrl_val = (ctrl_raw & 1) != 0;

            let mut bb = B::new();
            let q_target: Vec<QubitId> = bb.alloc_qubits(n);
            let q_dirty: Vec<QubitId> = bb.alloc_qubits(n.saturating_sub(2).max(1));
            let q_clean2: [QubitId; 2] = [bb.alloc_qubit(), bb.alloc_qubit()];
            let q_ctrl = bb.alloc_qubit();

            cisub_dirty_2clean_classical(&mut bb, &q_target, &q_dirty, &q_clean2, c, q_ctrl);

            let ops = bb.ops.clone();
            let num_qubits = bb.next_qubit as usize;
            let num_bits = bb.next_bit as usize;
            let mut inner_hasher = Shake256::default();
            inner_hasher.update(&[221u8]);
            let mut inner_xof =
                <sha3::Shake256 as sha3::digest::ExtendableOutput>::finalize_xof(inner_hasher);
            let mut sim = Simulator::new(num_qubits, num_bits, &mut inner_xof);
            sim.clear_for_shot();
            for k in 0..n {
                if (target >> k) & 1 != 0 {
                    *sim.qubit_mut(q_target[k]) = 1;
                }
            }
            for (k, &q) in q_dirty.iter().enumerate() {
                if (dirty_init >> k) & 1 != 0 {
                    *sim.qubit_mut(q) = 1;
                }
            }
            if ctrl_val {
                *sim.qubit_mut(q_ctrl) = 1;
            }
            sim.apply(&ops);

            let expected = if ctrl_val {
                target.wrapping_sub(c) & mask
            } else {
                target
            };
            let mut got: u64 = 0;
            for k in 0..n {
                if sim.qubit(q_target[k]) & 1 != 0 {
                    got |= 1 << k;
                }
            }
            let mut got_dirty: u64 = 0;
            for (k, &q) in q_dirty.iter().enumerate() {
                if sim.qubit(q) & 1 != 0 {
                    got_dirty |= 1 << k;
                }
            }
            let dirty_ok = got_dirty == (dirty_init & ((1u64 << q_dirty.len()) - 1).min(mask));
            let phase = sim.global_phase() & 1;

            if got == expected && dirty_ok && phase == 0 {
                ok += 1;
            } else {
                bad += 1;
                if bad < 3 {
                    eprintln!(
                        "cisub FAIL n={} t={:#x} c={:#x} ctrl={} got={:#x} exp={:#x} d_ok={} phase={}",
                        n, target, c, ctrl_val, got, expected, dirty_ok, phase
                    );
                }
            }
        }
        (ok, bad)
    }

    #[test]
    fn test_cisub_dirty_small() {
        for n in 5..=10 {
            let (ok, bad) = run_cisub_dirty(n, 10);
            assert_eq!(bad, 0, "n={n}: {ok}/{} passed", ok + bad);
        }
    }

    #[test]
    fn test_cisub_dirty_large() {
        let n = 256;
        let mut hasher = Shake256::default();
        hasher.update(&[n as u8, 50u8, 17]);
        use sha3::digest::XofReader;
        let mut xof = <sha3::Shake256 as sha3::digest::ExtendableOutput>::finalize_xof(hasher);
        let trials = 50;
        let c_low = 0x1_0000_03D1u64;
        let mut ok = 0;
        let mut bad = 0;
        for _trial in 0..trials {
            let mut buf = [0u8; 40];
            xof.read(&mut buf);
            let target = u64::from_le_bytes(buf[0..8].try_into().unwrap());
            let dirty_init = u64::from_le_bytes(buf[8..16].try_into().unwrap());
            let ctrl_val = (buf[16] & 1) != 0;

            let mut bb = B::new();
            let q_target: Vec<QubitId> = bb.alloc_qubits(n);
            let q_dirty: Vec<QubitId> = bb.alloc_qubits(n - 2);
            let q_clean2: [QubitId; 2] = [bb.alloc_qubit(), bb.alloc_qubit()];
            let q_ctrl = bb.alloc_qubit();

            cisub_dirty_2clean_classical(&mut bb, &q_target, &q_dirty, &q_clean2, c_low, q_ctrl);

            let ops = bb.ops.clone();
            let num_qubits = bb.next_qubit as usize;
            let num_bits = bb.next_bit as usize;
            let mut inner_hasher = Shake256::default();
            inner_hasher.update(&[19u8]);
            let mut inner_xof =
                <sha3::Shake256 as sha3::digest::ExtendableOutput>::finalize_xof(inner_hasher);
            let mut sim = Simulator::new(num_qubits, num_bits, &mut inner_xof);
            sim.clear_for_shot();
            for k in 0..64 {
                if (target >> k) & 1 != 0 {
                    *sim.qubit_mut(q_target[k]) = 1;
                }
            }
            for (k, &q) in q_dirty.iter().enumerate().take(64) {
                if (dirty_init >> k) & 1 != 0 {
                    *sim.qubit_mut(q) = 1;
                }
            }
            if ctrl_val {
                *sim.qubit_mut(q_ctrl) = 1;
            }
            sim.apply(&ops);

            let expected = if ctrl_val {
                target.wrapping_sub(c_low)
            } else {
                target
            };
            let mut got: u64 = 0;
            for k in 0..64 {
                if sim.qubit(q_target[k]) & 1 != 0 {
                    got |= 1 << k;
                }
            }
            let mut got_dirty: u64 = 0;
            for (k, &q) in q_dirty.iter().enumerate().take(64) {
                if sim.qubit(q) & 1 != 0 {
                    got_dirty |= 1 << k;
                }
            }
            let dirty_ok = got_dirty == dirty_init;
            let phase = sim.global_phase() & 1;
            let ctrl_preserved = sim.qubit(q_ctrl) & 1 == (ctrl_val as u64);
            if got == expected && dirty_ok && phase == 0 && ctrl_preserved {
                ok += 1;
            } else {
                bad += 1;
                if bad < 3 {
                    eprintln!(
                        "cisub n=256 FAIL t={:#x} d={:#x} ctrl={} got={:#x} exp={:#x} d_ok={} phase={} ctrl_ok={}",
                        target, dirty_init, ctrl_val, got, expected, dirty_ok, phase, ctrl_preserved
                    );
                }
            }
        }
        assert_eq!(bad, 0, "n=256 cisub: {ok}/{trials} passed");
    }

    #[test]
    fn test_cisub_dirty_kaliski_pattern() {

        let n = 256;
        let c_low = 0x1_0000_03D1u64;
        let trials = 50;
        let mut hasher = Shake256::default();
        hasher.update(&[99u8]);
        use sha3::digest::XofReader;
        let mut xof = <sha3::Shake256 as sha3::digest::ExtendableOutput>::finalize_xof(hasher);
        let mut ok = 0;
        let mut bad = 0;
        for _trial in 0..trials {
            let mut buf = [0u8; 16];
            xof.read(&mut buf);
            let target = u64::from_le_bytes(buf[0..8].try_into().unwrap());
            let dirty_u_lsb = (buf[8] & 1) != 0;
            let ctrl_val = (buf[9] & 1) != 0;

            let mut bb = B::new();
            let q_target: Vec<QubitId> = bb.alloc_qubits(n);
            let q_dirty: Vec<QubitId> = bb.alloc_qubits(n - 2);
            let q_clean2: [QubitId; 2] = [bb.alloc_qubit(), bb.alloc_qubit()];
            let q_ctrl = bb.alloc_qubit();

            cisub_dirty_2clean_classical(&mut bb, &q_target, &q_dirty, &q_clean2, c_low, q_ctrl);

            let ops = bb.ops.clone();
            let num_qubits = bb.next_qubit as usize;
            let num_bits = bb.next_bit as usize;
            let mut inner_hasher = Shake256::default();
            inner_hasher.update(&[77u8]);
            let mut inner_xof =
                <sha3::Shake256 as sha3::digest::ExtendableOutput>::finalize_xof(inner_hasher);
            let mut sim = Simulator::new(num_qubits, num_bits, &mut inner_xof);
            sim.clear_for_shot();
            for k in 0..64 {
                if (target >> k) & 1 != 0 {
                    *sim.qubit_mut(q_target[k]) = 1;
                }
            }

            if dirty_u_lsb {
                *sim.qubit_mut(q_dirty[0]) = 1;
            }
            if ctrl_val {
                *sim.qubit_mut(q_ctrl) = 1;
            }
            sim.apply(&ops);

            let expected = if ctrl_val {
                target.wrapping_sub(c_low)
            } else {
                target
            };
            let mut got: u64 = 0;
            for k in 0..64 {
                if sim.qubit(q_target[k]) & 1 != 0 {
                    got |= 1 << k;
                }
            }
            let got_dirty0 = sim.qubit(q_dirty[0]) & 1 != 0;
            let dirty_ok = got_dirty0 == dirty_u_lsb;
            let phase = sim.global_phase() & 1;
            let ctrl_preserved = sim.qubit(q_ctrl) & 1 == (ctrl_val as u64);
            if got == expected && dirty_ok && phase == 0 && ctrl_preserved {
                ok += 1;
            } else {
                bad += 1;
                if bad < 3 {
                    eprintln!(
                        "cisub kaliski FAIL t={:#x} d0={} ctrl={} got={:#x} exp={:#x} d_ok={} phase={} ctrl_ok={}",
                        target, dirty_u_lsb, ctrl_val, got, expected, dirty_ok, phase, ctrl_preserved
                    );
                }
            }
        }
        assert_eq!(bad, 0, "kaliski pattern cisub: {ok}/{trials} passed");
    }

    fn run_iadd_qoffset_dirty(n: usize, trials: usize) -> (usize, usize) {
        let mut hasher = Shake256::default();
        hasher.update(&[n as u8, trials as u8, 199]);
        use sha3::digest::XofReader;
        let mut xof = <sha3::Shake256 as sha3::digest::ExtendableOutput>::finalize_xof(hasher);
        let mut ok = 0;
        let mut bad = 0;
        for _trial in 0..trials {
            let mut buf = [0u8; 40];
            xof.read(&mut buf);
            let target_raw = u64::from_le_bytes(buf[0..8].try_into().unwrap());
            let offset_raw = u64::from_le_bytes(buf[8..16].try_into().unwrap());
            let dirty_raw = u64::from_le_bytes(buf[16..24].try_into().unwrap());
            let cin_raw = buf[24];
            let mask = if n < 64 { (1u64 << n) - 1 } else { u64::MAX };
            let target = target_raw & mask;
            let offset = offset_raw & mask;
            let dirty_init = dirty_raw & mask;
            let cin = (cin_raw & 1) != 0;

            let mut bb = B::new();
            let q_target: Vec<QubitId> = bb.alloc_qubits(n);
            let q_offset: Vec<QubitId> = bb.alloc_qubits(n);
            let q_dirty: Vec<QubitId> = bb.alloc_qubits(n.saturating_sub(2).max(1));
            let q_clean2: [QubitId; 2] = [bb.alloc_qubit(), bb.alloc_qubit()];

            iadd_dirty_2clean_qoffset(&mut bb, &q_target, &q_dirty, &q_clean2, &q_offset, cin);

            let ops = bb.ops.clone();
            let num_qubits = bb.next_qubit as usize;
            let num_bits = bb.next_bit as usize;
            let mut inner_hasher = Shake256::default();
            inner_hasher.update(&[231u8]);
            let mut inner_xof =
                <sha3::Shake256 as sha3::digest::ExtendableOutput>::finalize_xof(inner_hasher);
            let mut sim = Simulator::new(num_qubits, num_bits, &mut inner_xof);
            sim.clear_for_shot();
            for k in 0..n {
                if (target >> k) & 1 != 0 {
                    *sim.qubit_mut(q_target[k]) = 1;
                }
                if (offset >> k) & 1 != 0 {
                    *sim.qubit_mut(q_offset[k]) = 1;
                }
            }
            for (k, &q) in q_dirty.iter().enumerate() {
                if (dirty_init >> k) & 1 != 0 {
                    *sim.qubit_mut(q) = 1;
                }
            }
            sim.apply(&ops);

            let expected = (target.wrapping_add(offset).wrapping_add(cin as u64)) & mask;
            let mut got: u64 = 0;
            for k in 0..n {
                if sim.qubit(q_target[k]) & 1 != 0 {
                    got |= 1 << k;
                }
            }

            let mut got_offset: u64 = 0;
            for k in 0..n {
                if sim.qubit(q_offset[k]) & 1 != 0 {
                    got_offset |= 1 << k;
                }
            }
            let mut got_dirty: u64 = 0;
            for (k, &q) in q_dirty.iter().enumerate() {
                if sim.qubit(q) & 1 != 0 {
                    got_dirty |= 1 << k;
                }
            }
            let dirty_ok = got_dirty == (dirty_init & ((1u64 << q_dirty.len()) - 1).min(mask));
            let offset_ok = got_offset == offset;
            let phase = sim.global_phase() & 1;

            if got == expected && dirty_ok && offset_ok && phase == 0 {
                ok += 1;
            } else {
                bad += 1;
                if bad < 3 {
                    eprintln!(
                        "iadd_qoffset FAIL n={} t={:#x} o={:#x} d={:#x} cin={} got={:#x} exp={:#x} d_ok={} o_ok={} phase={}",
                        n, target, offset, dirty_init, cin, got, expected, dirty_ok, offset_ok, phase
                    );
                }
            }
        }
        (ok, bad)
    }

    #[test]
    fn test_iadd_qoffset_dirty_small() {
        for n in 5..=10 {
            let (ok, bad) = run_iadd_qoffset_dirty(n, 10);
            assert_eq!(bad, 0, "n={n}: {ok}/{} passed", ok + bad);
        }
    }

    fn run_iadd_qoffset_narrow(n: usize, m: usize, trials: usize) -> (usize, usize) {
        let mut hasher = Shake256::default();
        hasher.update(&[n as u8, m as u8, trials as u8, 211]);
        use sha3::digest::XofReader;
        let mut xof = <sha3::Shake256 as sha3::digest::ExtendableOutput>::finalize_xof(hasher);
        let mut ok = 0;
        let mut bad = 0;
        let nmask = if n < 64 { (1u64 << n) - 1 } else { u64::MAX };
        let mmask = if m < 64 { (1u64 << m) - 1 } else { u64::MAX };
        for _trial in 0..trials {
            let mut buf = [0u8; 40];
            xof.read(&mut buf);
            let target = u64::from_le_bytes(buf[0..8].try_into().unwrap()) & nmask;
            let offset = u64::from_le_bytes(buf[8..16].try_into().unwrap()) & mmask;
            let dirty_init = u64::from_le_bytes(buf[16..24].try_into().unwrap()) & nmask;
            let cin = (buf[24] & 1) != 0;

            let mut bb = B::new();
            let q_target: Vec<QubitId> = bb.alloc_qubits(n);
            let q_offset: Vec<QubitId> = bb.alloc_qubits(m);
            let q_dirty: Vec<QubitId> = bb.alloc_qubits(n.saturating_sub(2).max(1));
            let q_clean2: [QubitId; 2] = [bb.alloc_qubit(), bb.alloc_qubit()];

            iadd_dirty_2clean_qoffset_narrow(
                &mut bb, &q_target, &q_dirty, &q_clean2, &q_offset, cin,
            );

            let ops = bb.ops.clone();
            let num_qubits = bb.next_qubit as usize;
            let num_bits = bb.next_bit as usize;
            let mut inner_hasher = Shake256::default();
            inner_hasher.update(&[233u8]);
            let mut inner_xof =
                <sha3::Shake256 as sha3::digest::ExtendableOutput>::finalize_xof(inner_hasher);
            let mut sim = Simulator::new(num_qubits, num_bits, &mut inner_xof);
            sim.clear_for_shot();
            for k in 0..n {
                if k < 64 && (target >> k) & 1 != 0 {
                    *sim.qubit_mut(q_target[k]) = 1;
                }
            }
            for k in 0..m {
                if k < 64 && (offset >> k) & 1 != 0 {
                    *sim.qubit_mut(q_offset[k]) = 1;
                }
            }
            for (k, &q) in q_dirty.iter().enumerate() {
                if k < 64 && (dirty_init >> k) & 1 != 0 {
                    *sim.qubit_mut(q) = 1;
                }
            }
            sim.apply(&ops);

            let expected: u128 = (target as u128) + (offset as u128) + (cin as u128);
            let readbits = n.min(127);
            let expected = expected & ((1u128 << readbits) - 1);
            let mut got: u128 = 0;
            for k in 0..readbits {
                if sim.qubit(q_target[k]) & 1 != 0 {
                    got |= 1u128 << k;
                }
            }
            let mut got_offset: u64 = 0;
            for k in 0..m {
                if sim.qubit(q_offset[k]) & 1 != 0 {
                    got_offset |= 1 << k;
                }
            }
            let mut got_dirty: u64 = 0;
            for (k, &q) in q_dirty.iter().enumerate() {
                if sim.qubit(q) & 1 != 0 {
                    got_dirty |= 1 << k;
                }
            }
            let dmask = if q_dirty.len() < 64 {
                (1u64 << q_dirty.len()) - 1
            } else {
                u64::MAX
            };
            let dirty_ok = got_dirty == (dirty_init & dmask & nmask);
            let offset_ok = got_offset == offset;
            let phase = sim.global_phase() & 1;

            if got == expected && dirty_ok && offset_ok && phase == 0 {
                ok += 1;
            } else {
                bad += 1;
                if bad < 3 {
                    eprintln!(
                        "narrow FAIL n={} m={} t={:#x} o={:#x} d={:#x} cin={} got={:#x} exp={:#x} d_ok={} o_ok={} phase={}",
                        n, m, target, offset, dirty_init, cin, got, expected, dirty_ok, offset_ok, phase
                    );
                }
            }
        }
        (ok, bad)
    }

    #[test]
    fn test_iadd_qoffset_narrow_small() {

        for n in 5..=12 {
            for m in 1..n {
                let (ok, bad) = run_iadd_qoffset_narrow(n, m, 12);
                assert_eq!(bad, 0, "n={n} m={m}: {ok}/{} passed", ok + bad);
            }
        }
    }

    #[test]
    fn test_iadd_qoffset_narrow_wide() {

        let (ok, bad) = run_iadd_qoffset_narrow(256, 22, 40);
        assert_eq!(bad, 0, "n=256 m=22: {ok}/{} passed", ok + bad);
    }
}
