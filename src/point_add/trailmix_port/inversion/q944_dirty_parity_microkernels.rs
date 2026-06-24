//! Standalone allocation-free arithmetic using one arbitrary dirty carry lane.
//!
//! This module is deliberately not connected to the PZ route. It isolates the
//! gate identities needed to decide whether a dirty parity lane could replace
//! the clean Q945 arithmetic carry. Every lender is restored exactly.

use crate::circuit::{Op, OperationType};
use crate::point_add::trailmix_port::arith::mcx::mcx_dirty_ladder;
use crate::point_add::trailmix_port::circuit::{Circuit, QReg};
use crate::point_add::trailmix_port::inversion::q945_local_hosts::Q945_NON_HCLZ_ROWS;
use crate::point_add::trailmix_port::inversion::q949_robust_envelope::{
    q949_robust_clz_lows, q949_robust_pair_symmetric_widths,
};
use crate::point_add::B;

fn assert_distinct(roles: &[(&str, &QReg)]) {
    for (i, (left_name, left)) in roles.iter().enumerate() {
        for (right_name, right) in &roles[..i] {
            assert!(
                !std::ptr::eq(*left, *right),
                "Q944 dirty-parity alias: {left_name} aliases {right_name}"
            );
        }
    }
}

fn assert_arithmetic_layout(gate: &QReg, carry: &QReg, a: &[&QReg], b: &[&QReg]) {
    assert_eq!(a.len(), b.len(), "Q944 dirty-parity width mismatch");
    let mut roles = Vec::with_capacity(2 + a.len() + b.len());
    roles.push(("gate", gate));
    roles.push(("carry", carry));
    roles.extend(a.iter().map(|&q| ("a", q)));
    roles.extend(b.iter().map(|&q| ("b", q)));
    assert_distinct(&roles);
}

/// Literal controlled Cuccaro ripple with the caller's arbitrary carry `d`.
///
/// Gate-by-gate, one bit maps incoming carry `k` to
/// `MAJ(a_i,b_i,k)`. The reverse pass restores `d`, `b`, and `gate`, while
/// producing
///
/// `a -> a + gate * (b + d) (mod 2^n)`.
///
/// This raw identity is exposed only so the independent proof binary can test
/// the research claim directly. Full-route code must use the corrected forms.
#[doc(hidden)]
pub fn controlled_add_dirty_carry_raw_refs(
    circ: &mut Circuit,
    gate: &QReg,
    carry: &QReg,
    a: &[&QReg],
    b: &[&QReg],
) {
    assert_arithmetic_layout(gate, carry, a, b);
    for i in 0..a.len() {
        circ.cx(carry, b[i]);
        circ.cx(carry, a[i]);
        circ.ccx(a[i], b[i], carry);
    }
    for i in (0..a.len()).rev() {
        circ.ccx(a[i], b[i], carry);
        circ.cx(carry, a[i]);
        circ.ccx(gate, b[i], a[i]);
        circ.cx(carry, b[i]);
    }
}

/// Toggle `a -= gate*carry (mod 2^n)` with `b` as restored dirty lenders.
fn controlled_decrement_dirty_lenders(
    circ: &mut Circuit,
    gate: &QReg,
    carry: &QReg,
    a: &[&QReg],
    b: &[&QReg],
) {
    // High-to-low order leaves every lower bit at its pre-decrement value.
    for j in (0..a.len()).rev() {
        for &q in &a[..j] {
            circ.x(q);
        }
        let mut controls = Vec::with_capacity(j + 2);
        controls.extend([gate, carry]);
        controls.extend_from_slice(&a[..j]);
        mcx_dirty_ladder(circ, &controls, a[j], &b[..j]);
        for &q in a[..j].iter().rev() {
            circ.x(q);
        }
    }
}

/// Toggle `a += gate*carry (mod 2^n)` with `b` as restored dirty lenders.
fn controlled_increment_dirty_lenders(
    circ: &mut Circuit,
    gate: &QReg,
    carry: &QReg,
    a: &[&QReg],
    b: &[&QReg],
) {
    // High-to-low order leaves every lower bit at its pre-increment value.
    for j in (0..a.len()).rev() {
        let mut controls = Vec::with_capacity(j + 2);
        controls.extend([gate, carry]);
        controls.extend_from_slice(&a[..j]);
        mcx_dirty_ladder(circ, &controls, a[j], &b[..j]);
    }
}

/// Allocation-free corrected controlled addition with arbitrary dirty carry.
///
/// Semantics: `a -> a + gate*b (mod 2^n)`. `gate`, `carry`, and every `b`
/// lender are restored. The raw ripple contributes `gate*carry`, which the
/// final controlled decrement removes.
pub fn controlled_add_dirty_carry_refs(
    circ: &mut Circuit,
    gate: &QReg,
    carry: &QReg,
    a: &[&QReg],
    b: &[&QReg],
) {
    let section = circ.push_section("q944.dirty-carry-add");
    controlled_add_dirty_carry_raw_refs(circ, gate, carry, a, b);
    controlled_decrement_dirty_lenders(circ, gate, carry, a, b);
    circ.pop_section(&section);
}

/// Allocation-free corrected controlled subtraction with arbitrary dirty carry.
///
/// The X-bracketed raw ripple gives `a - gate*(b+carry)`. The final controlled
/// increment removes the `-gate*carry` term. All non-target lanes are restored.
pub fn controlled_sub_dirty_carry_refs(
    circ: &mut Circuit,
    gate: &QReg,
    carry: &QReg,
    a: &[&QReg],
    b: &[&QReg],
) {
    let section = circ.push_section("q944.dirty-carry-sub");
    assert_arithmetic_layout(gate, carry, a, b);
    for &q in a {
        circ.x(q);
    }
    controlled_add_dirty_carry_raw_refs(circ, gate, carry, a, b);
    for &q in a {
        circ.x(q);
    }
    controlled_increment_dirty_lenders(circ, gate, carry, a, b);
    circ.pop_section(&section);
}

fn assert_comparator_layout(
    gate: &QReg,
    carry: &QReg,
    v: &[&QReg],
    u: &[&QReg],
    out: &QReg,
) {
    assert!(!v.is_empty(), "Q944 dirty-parity comparator requires n >= 1");
    assert_eq!(v.len(), u.len(), "Q944 dirty-parity comparator width mismatch");
    let mut roles = Vec::with_capacity(3 + v.len() + u.len());
    roles.extend([("gate", gate), ("carry", carry), ("out", out)]);
    roles.extend(v.iter().map(|&q| ("v", q)));
    roles.extend(u.iter().map(|&q| ("u", q)));
    assert_distinct(&roles);
}

/// Raw gated strict comparator with arbitrary dirty carry.
///
/// The final carry of `u + !v + d` is
/// `[v<u] XOR d*[v=u]`. The operation therefore toggles
/// `out` by `gate * ([v<u] XOR d*[v=u])`, then restores all other lanes.
#[doc(hidden)]
pub fn strict_compare_gated_dirty_carry_raw_refs(
    circ: &mut Circuit,
    v: &[&QReg],
    u: &[&QReg],
    gate: &QReg,
    out: &QReg,
    carry: &QReg,
) {
    assert_comparator_layout(gate, carry, v, u, out);
    for &q in v {
        circ.x(q);
    }
    circ.cx(v[0], u[0]);
    circ.cx(v[0], carry);
    circ.ccx(carry, u[0], v[0]);
    for i in 1..v.len() {
        circ.cx(v[i], u[i]);
        circ.cx(v[i], v[i - 1]);
        circ.ccx(v[i - 1], u[i], v[i]);
    }
    circ.ccx(gate, v[v.len() - 1], out);
    for i in (1..v.len()).rev() {
        circ.ccx(v[i - 1], u[i], v[i]);
        circ.cx(v[i], v[i - 1]);
        circ.cx(v[i], u[i]);
    }
    circ.ccx(carry, u[0], v[0]);
    circ.cx(v[0], carry);
    circ.cx(v[0], u[0]);
    for &q in v {
        circ.x(q);
    }
}

/// Allocation-free corrected gated strict comparator.
///
/// Semantics: `out ^= gate * [v<u]`. `v`, `u`, `gate`, and the arbitrary
/// dirty `carry` are restored. Equality is encoded temporarily into `v`; `u`
/// supplies exactly `n` restored dirty lenders for the `(n+2)`-controlled
/// correction `gate*carry*[v=u]`.
pub fn strict_compare_gated_dirty_carry_refs(
    circ: &mut Circuit,
    v: &[&QReg],
    u: &[&QReg],
    gate: &QReg,
    out: &QReg,
    carry: &QReg,
) {
    let section = circ.push_section("q944.dirty-carry-strict-compare");
    strict_compare_gated_dirty_carry_raw_refs(circ, v, u, gate, out, carry);

    for i in 0..v.len() {
        circ.cx(u[i], v[i]);
        circ.x(v[i]);
    }
    let mut controls = Vec::with_capacity(v.len() + 2);
    controls.extend([gate, carry]);
    controls.extend_from_slice(v);
    mcx_dirty_ladder(circ, &controls, out, u);
    for i in (0..v.len()).rev() {
        circ.x(v[i]);
        circ.cx(u[i], v[i]);
    }
    circ.pop_section(&section);
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct Q944GateCounts {
    pub x: usize,
    pub cx: usize,
    pub ccx: usize,
    pub total: usize,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct Q944WidthCounts {
    pub width: usize,
    pub add: Q944GateCounts,
    pub sub: Q944GateCounts,
    pub strict_compare: Q944GateCounts,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct Q944PeakRowCounts {
    pub row: usize,
    pub division_arithmetic_width: usize,
    pub multiply_arithmetic_width: usize,
    pub division_compare_width: usize,
    pub multiply_compare_width: usize,
    pub division_add: Q944GateCounts,
    pub division_sub: Q944GateCounts,
    pub multiply_add: Q944GateCounts,
    pub multiply_sub: Q944GateCounts,
    pub division_compare: Q944GateCounts,
    pub multiply_compare: Q944GateCounts,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Q944ExhaustiveReport {
    pub widths_checked: usize,
    pub raw_compare_gate_shape_checks: usize,
    pub raw_add_states_checked: usize,
    pub corrected_add_states_checked: usize,
    pub corrected_sub_states_checked: usize,
    pub raw_compare_states_checked: usize,
    pub corrected_compare_states_checked: usize,
    pub phase_clean_gate_streams_checked: usize,
    pub allocation_free_gate_streams_checked: usize,
    pub width_counts: Vec<Q944WidthCounts>,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Q944StructuralReport {
    pub rows_checked: usize,
    pub emitted_streams_checked: usize,
    pub phase_clean_gate_streams_checked: usize,
    pub allocation_free_gate_streams_checked: usize,
    pub rows: Vec<Q944PeakRowCounts>,
}

#[must_use]
pub const fn q944_correction_toffoli(width: usize) -> usize {
    2 * width * width - 2 * width + 1
}

#[must_use]
pub const fn q944_add_counts(width: usize) -> Q944GateCounts {
    let ccx = 3 * width + q944_correction_toffoli(width);
    let x = width * (width - 1);
    let cx = 4 * width;
    Q944GateCounts {
        x,
        cx,
        ccx,
        total: x + cx + ccx,
    }
}

#[must_use]
pub const fn q944_sub_counts(width: usize) -> Q944GateCounts {
    let ccx = 3 * width + q944_correction_toffoli(width);
    let x = 2 * width;
    let cx = 4 * width;
    Q944GateCounts {
        x,
        cx,
        ccx,
        total: x + cx + ccx,
    }
}

#[must_use]
pub const fn q944_strict_compare_counts(width: usize) -> Q944GateCounts {
    let x = 4 * width;
    let cx = 6 * width;
    let ccx = 6 * width + 1;
    Q944GateCounts {
        x,
        cx,
        ccx,
        total: x + cx + ccx,
    }
}

fn gate_counts(ops: &[Op]) -> Q944GateCounts {
    let mut counts = Q944GateCounts::default();
    for op in ops {
        match op.kind {
            OperationType::X => counts.x += 1,
            OperationType::CX => counts.cx += 1,
            OperationType::CCX => counts.ccx += 1,
            other => panic!("Q944 microkernel emitted phase/non-classical gate {other:?}"),
        }
    }
    counts.total = ops.len();
    assert_eq!(counts.total, counts.x + counts.cx + counts.ccx);
    counts
}

#[derive(Clone, Copy)]
enum ArithmeticKind {
    RawAdd,
    Add,
    Sub,
}

fn build_arithmetic(width: usize, kind: ArithmeticKind) -> B {
    let mut circ = Circuit::new();
    let gate = circ.alloc_qreg("q944.gate");
    let carry = circ.alloc_qreg("q944.dirty-parity");
    let a = circ.alloc_qreg_bits("q944.a", width);
    let b = circ.alloc_qreg_bits("q944.b", width);
    let ar: Vec<&QReg> = a.iter().collect();
    let br: Vec<&QReg> = b.iter().collect();
    match kind {
        ArithmeticKind::RawAdd => {
            controlled_add_dirty_carry_raw_refs(&mut circ, &gate, &carry, &ar, &br)
        }
        ArithmeticKind::Add => {
            controlled_add_dirty_carry_refs(&mut circ, &gate, &carry, &ar, &br)
        }
        ArithmeticKind::Sub => {
            controlled_sub_dirty_carry_refs(&mut circ, &gate, &carry, &ar, &br)
        }
    }
    let builder = circ.into_builder();
    let input_qubits = 2 * width + 2;
    assert_eq!(builder.next_qubit as usize, input_qubits);
    assert_eq!(builder.active_qubits as usize, input_qubits);
    assert_eq!(builder.peak_qubits as usize, input_qubits);
    drop((ar, br));
    drop((gate, carry, a, b));
    builder
}

fn build_comparator(width: usize, raw: bool) -> B {
    let mut circ = Circuit::new();
    let gate = circ.alloc_qreg("q944.gate");
    let carry = circ.alloc_qreg("q944.dirty-parity");
    let v = circ.alloc_qreg_bits("q944.v", width);
    let u = circ.alloc_qreg_bits("q944.u", width);
    let out = circ.alloc_qreg("q944.out");
    let vr: Vec<&QReg> = v.iter().collect();
    let ur: Vec<&QReg> = u.iter().collect();
    if raw {
        strict_compare_gated_dirty_carry_raw_refs(
            &mut circ, &vr, &ur, &gate, &out, &carry,
        );
    } else {
        strict_compare_gated_dirty_carry_refs(&mut circ, &vr, &ur, &gate, &out, &carry);
    }
    let builder = circ.into_builder();
    let input_qubits = 2 * width + 3;
    assert_eq!(builder.next_qubit as usize, input_qubits);
    assert_eq!(builder.active_qubits as usize, input_qubits);
    assert_eq!(builder.peak_qubits as usize, input_qubits);
    drop((vr, ur));
    drop((gate, carry, v, u, out));
    builder
}

fn apply_scalar(ops: &[Op], mut state: u64) -> u64 {
    let bit = |word: u64, id: u64| ((word >> id) & 1) != 0;
    for op in ops {
        match op.kind {
            OperationType::X => state ^= 1u64 << op.q_target.0,
            OperationType::CX => {
                if bit(state, op.q_control1.0) {
                    state ^= 1u64 << op.q_target.0;
                }
            }
            OperationType::CCX => {
                if bit(state, op.q_control1.0) && bit(state, op.q_control2.0) {
                    state ^= 1u64 << op.q_target.0;
                }
            }
            other => panic!("Q944 scalar proof saw unexpected gate {other:?}"),
        }
    }
    state
}

/// Bind the raw comparator to the exact existing MAJ/un-MAJ cascade. This
/// shape check runs before any semantic interpretation and rejects accidental
/// reuse of the physical carry beyond bit zero.
fn assert_raw_comparator_gate_shape(ops: &[Op], width: usize) {
    let gate = 0u64;
    let carry = 1u64;
    let v = |i: usize| 2 + i as u64;
    let u = |i: usize| 2 + width as u64 + i as u64;
    let out = 2 + 2 * width as u64;
    let mut cursor = 0usize;

    let next_x = |cursor: &mut usize, target: u64| {
        let op = &ops[*cursor];
        assert_eq!((op.kind, op.q_target.0), (OperationType::X, target));
        *cursor += 1;
    };
    let next_cx = |cursor: &mut usize, control: u64, target: u64| {
        let op = &ops[*cursor];
        assert_eq!(
            (op.kind, op.q_control1.0, op.q_target.0),
            (OperationType::CX, control, target)
        );
        *cursor += 1;
    };
    let next_ccx = |cursor: &mut usize, control1: u64, control2: u64, target: u64| {
        let op = &ops[*cursor];
        assert_eq!(
            (
                op.kind,
                op.q_control2.0,
                op.q_control1.0,
                op.q_target.0,
            ),
            (OperationType::CCX, control1, control2, target)
        );
        *cursor += 1;
    };

    for i in 0..width {
        next_x(&mut cursor, v(i));
    }
    next_cx(&mut cursor, v(0), u(0));
    next_cx(&mut cursor, v(0), carry);
    next_ccx(&mut cursor, carry, u(0), v(0));
    for i in 1..width {
        next_cx(&mut cursor, v(i), u(i));
        next_cx(&mut cursor, v(i), v(i - 1));
        next_ccx(&mut cursor, v(i - 1), u(i), v(i));
    }
    next_ccx(&mut cursor, gate, v(width - 1), out);
    for i in (1..width).rev() {
        next_ccx(&mut cursor, v(i - 1), u(i), v(i));
        next_cx(&mut cursor, v(i), v(i - 1));
        next_cx(&mut cursor, v(i), u(i));
    }
    next_ccx(&mut cursor, carry, u(0), v(0));
    next_cx(&mut cursor, v(0), carry);
    next_cx(&mut cursor, v(0), u(0));
    for i in 0..width {
        next_x(&mut cursor, v(i));
    }
    assert_eq!(cursor, ops.len());
}

fn check_arithmetic_basis(width: usize, builder: &B, kind: ArithmeticKind) -> usize {
    let states = 1usize << (2 * width + 2);
    let mask = (1u64 << width) - 1;
    for input in 0..states as u64 {
        let gate = input & 1;
        let carry = (input >> 1) & 1;
        let a = (input >> 2) & mask;
        let b = (input >> (width + 2)) & mask;
        let output = apply_scalar(&builder.ops, input);
        let got_gate = output & 1;
        let got_carry = (output >> 1) & 1;
        let got_a = (output >> 2) & mask;
        let got_b = (output >> (width + 2)) & mask;
        let want_a = match kind {
            ArithmeticKind::RawAdd if gate != 0 => a.wrapping_add(b + carry) & mask,
            ArithmeticKind::Add if gate != 0 => a.wrapping_add(b) & mask,
            ArithmeticKind::Sub if gate != 0 => a.wrapping_sub(b) & mask,
            _ => a,
        };
        assert_eq!((got_gate, got_carry), (gate, carry));
        assert_eq!(got_b, b, "width={width}: dirty lender changed");
        assert_eq!(got_a, want_a, "width={width} input={input:#x}");
    }
    states
}

fn check_comparator_basis(width: usize, builder: &B, raw: bool) -> usize {
    let states = 1usize << (2 * width + 3);
    let mask = (1u64 << width) - 1;
    for input in 0..states as u64 {
        let gate = input & 1;
        let carry = (input >> 1) & 1;
        let v = (input >> 2) & mask;
        let u = (input >> (width + 2)) & mask;
        let out = (input >> (2 * width + 2)) & 1;
        let output = apply_scalar(&builder.ops, input);
        let got_gate = output & 1;
        let got_carry = (output >> 1) & 1;
        let got_v = (output >> 2) & mask;
        let got_u = (output >> (width + 2)) & mask;
        let got_out = (output >> (2 * width + 2)) & 1;
        let strict = u64::from(v < u);
        let equal_error = carry & u64::from(v == u);
        let predicate = if raw { strict ^ equal_error } else { strict };
        let want_out = out ^ (gate & predicate);
        assert_eq!((got_gate, got_carry), (gate, carry));
        assert_eq!((got_v, got_u), (v, u), "width={width}: operand changed");
        assert_eq!(got_out, want_out, "width={width} input={input:#x}");
    }
    states
}

/// Exhaustive scalar interpretation over all basis states for widths 1..=5.
#[must_use]
pub fn exhaustive_q944_dirty_parity_microkernels_check() -> Q944ExhaustiveReport {
    let mut raw_add_states_checked = 0;
    let mut corrected_add_states_checked = 0;
    let mut corrected_sub_states_checked = 0;
    let mut raw_compare_states_checked = 0;
    let mut corrected_compare_states_checked = 0;
    let mut width_counts = Vec::new();

    for width in 1..=5 {
        let raw_add = build_arithmetic(width, ArithmeticKind::RawAdd);
        let add = build_arithmetic(width, ArithmeticKind::Add);
        let sub = build_arithmetic(width, ArithmeticKind::Sub);
        let raw_compare = build_comparator(width, true);
        let strict_compare = build_comparator(width, false);

        assert_raw_comparator_gate_shape(&raw_compare.ops, width);
        raw_add_states_checked += check_arithmetic_basis(width, &raw_add, ArithmeticKind::RawAdd);
        corrected_add_states_checked += check_arithmetic_basis(width, &add, ArithmeticKind::Add);
        corrected_sub_states_checked += check_arithmetic_basis(width, &sub, ArithmeticKind::Sub);
        raw_compare_states_checked += check_comparator_basis(width, &raw_compare, true);
        corrected_compare_states_checked += check_comparator_basis(width, &strict_compare, false);

        let add_counts = gate_counts(&add.ops);
        let sub_counts = gate_counts(&sub.ops);
        let compare_counts = gate_counts(&strict_compare.ops);
        assert_eq!(gate_counts(&raw_add.ops), Q944GateCounts {
            x: 0,
            cx: 4 * width,
            ccx: 3 * width,
            total: 7 * width,
        });
        assert_eq!(gate_counts(&raw_compare.ops), Q944GateCounts {
            x: 2 * width,
            cx: 4 * width,
            ccx: 2 * width + 1,
            total: 8 * width + 1,
        });
        assert_eq!(add_counts, q944_add_counts(width));
        assert_eq!(sub_counts, q944_sub_counts(width));
        assert_eq!(compare_counts, q944_strict_compare_counts(width));
        width_counts.push(Q944WidthCounts {
            width,
            add: add_counts,
            sub: sub_counts,
            strict_compare: compare_counts,
        });
    }

    Q944ExhaustiveReport {
        widths_checked: width_counts.len(),
        raw_compare_gate_shape_checks: width_counts.len(),
        raw_add_states_checked,
        corrected_add_states_checked,
        corrected_sub_states_checked,
        raw_compare_states_checked,
        corrected_compare_states_checked,
        phase_clean_gate_streams_checked: 5 * 5,
        allocation_free_gate_streams_checked: 5 * 5,
        width_counts,
    }
}

fn checked_arithmetic_counts(width: usize, kind: ArithmeticKind) -> Q944GateCounts {
    let builder = build_arithmetic(width, kind);
    let counts = gate_counts(&builder.ops);
    match kind {
        ArithmeticKind::RawAdd => unreachable!("structural check does not use raw add"),
        ArithmeticKind::Add => assert_eq!(counts, q944_add_counts(width)),
        ArithmeticKind::Sub => assert_eq!(counts, q944_sub_counts(width)),
    }
    counts
}

fn checked_compare_counts(width: usize) -> Q944GateCounts {
    let builder = build_comparator(width, false);
    let counts = gate_counts(&builder.ops);
    assert_eq!(counts, q944_strict_compare_counts(width));
    counts
}

/// Emit and count every arithmetic width occurring at the seven Q945 peak rows.
/// This is a structural microkernel check only; it does not integrate a route.
#[must_use]
pub fn q944_q945_peak_width_emission_check() -> Q944StructuralReport {
    const EXPECTED: [(usize, usize, usize, usize, usize); 7] = [
        (363, 81, 248, 80, 77),
        (364, 81, 248, 80, 74),
        (374, 74, 254, 74, 76),
        (375, 73, 255, 73, 76),
        (376, 73, 255, 73, 76),
        (379, 72, 256, 72, 77),
        (380, 72, 256, 72, 77),
    ];
    assert_eq!(Q945_NON_HCLZ_ROWS, EXPECTED.map(|entry| entry.0));

    let mut rows = Vec::new();
    for &(row, div_width, mul_width, div_cmp_width, mul_cmp_width) in &EXPECTED {
        let widths = q949_robust_pair_symmetric_widths(row);
        let lows = q949_robust_clz_lows(row);
        assert_eq!((widths[0], widths[2]), (div_width, mul_width));
        assert_eq!(widths[0] - lows[0], div_cmp_width);
        assert_eq!(widths[2] - lows[2], mul_cmp_width);

        rows.push(Q944PeakRowCounts {
            row,
            division_arithmetic_width: div_width,
            multiply_arithmetic_width: mul_width,
            division_compare_width: div_cmp_width,
            multiply_compare_width: mul_cmp_width,
            division_add: checked_arithmetic_counts(div_width, ArithmeticKind::Add),
            division_sub: checked_arithmetic_counts(div_width, ArithmeticKind::Sub),
            multiply_add: checked_arithmetic_counts(mul_width, ArithmeticKind::Add),
            multiply_sub: checked_arithmetic_counts(mul_width, ArithmeticKind::Sub),
            division_compare: checked_compare_counts(div_cmp_width),
            multiply_compare: checked_compare_counts(mul_cmp_width),
        });
    }

    let emitted_streams_checked = rows.len() * 6;
    Q944StructuralReport {
        rows_checked: rows.len(),
        emitted_streams_checked,
        phase_clean_gate_streams_checked: emitted_streams_checked,
        allocation_free_gate_streams_checked: emitted_streams_checked,
        rows,
    }
}
