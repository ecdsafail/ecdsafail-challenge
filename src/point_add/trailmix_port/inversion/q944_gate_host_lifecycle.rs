//! Standalone proof for a zero lane hosting the complete active-and-less gate.
//!
//! The independently verified dirty-parity comparator remains unchanged. This
//! module composes it as compute/body/uncompute and checks the composition over
//! every basis state at small widths before any production integration.

use crate::circuit::{Op, OperationType};
use crate::point_add::trailmix_port::circuit::{Circuit, QReg};
use crate::point_add::trailmix_port::inversion::q944_dirty_parity_microkernels::
    strict_compare_gated_dirty_carry_refs;
use crate::point_add::B;

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct Q944GateHostLifecycleCounts {
    pub x: usize,
    pub cx: usize,
    pub ccx: usize,
    pub total: usize,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Q944GateHostLifecycleReport {
    pub widths_checked: usize,
    pub basis_states_checked: usize,
    pub zero_host_entry_checks: usize,
    pub zero_host_exit_checks: usize,
    pub counter_restoration_checks: usize,
    pub parity_restoration_checks: usize,
    pub operand_restoration_checks: usize,
    pub body_semantic_checks: usize,
    pub allocation_free_streams_checked: usize,
    pub phase_clean_streams_checked: usize,
    pub width_counts: Vec<(usize, Q944GateHostLifecycleCounts)>,
}

fn build_hosted_lifecycle(width: usize) -> B {
    let mut circ = Circuit::new();
    let counter = circ.alloc_qreg("q944.hosted.counter");
    let parity = circ.alloc_qreg("q944.hosted.dirty-parity");
    let host = circ.alloc_qreg("q944.hosted.gate");
    let body_control = circ.alloc_qreg("q944.hosted.body-control");
    let body_target = circ.alloc_qreg("q944.hosted.body-target");
    let v = circ.alloc_qreg_bits("q944.hosted.v", width);
    let u = circ.alloc_qreg_bits("q944.hosted.u", width);
    let vr: Vec<&QReg> = v.iter().collect();
    let ur: Vec<&QReg> = u.iter().collect();

    let toggle_gate = |circ: &mut Circuit| {
        circ.x(&counter);
        strict_compare_gated_dirty_carry_refs(
            circ,
            &vr,
            &ur,
            &counter,
            &host,
            &parity,
        );
        circ.x(&counter);
    };
    toggle_gate(&mut circ);
    circ.ccx(&host, &body_control, &body_target);
    toggle_gate(&mut circ);
    let builder = circ.into_builder();
    let inputs = 2 * width + 5;
    assert_eq!(builder.next_qubit as usize, inputs);
    assert_eq!(builder.active_qubits as usize, inputs);
    assert_eq!(builder.peak_qubits as usize, inputs);
    drop((vr, ur));
    drop((counter, parity, host, body_control, body_target, v, u));
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
            other => panic!("Q944 hosted lifecycle emitted phase gate {other:?}"),
        }
    }
    state
}

fn gate_counts(ops: &[Op]) -> Q944GateHostLifecycleCounts {
    let mut counts = Q944GateHostLifecycleCounts::default();
    for op in ops {
        match op.kind {
            OperationType::X => counts.x += 1,
            OperationType::CX => counts.cx += 1,
            OperationType::CCX => counts.ccx += 1,
            other => panic!("Q944 hosted lifecycle emitted phase gate {other:?}"),
        }
    }
    counts.total = ops.len();
    assert_eq!(counts.total, counts.x + counts.cx + counts.ccx);
    counts
}

/// Exhaustively prove the complete hosted-gate lifecycle for widths 1 through
/// 5. The host is constrained to zero on entry; every other input, including
/// the comparator carry, is arbitrary.
#[must_use]
pub fn exhaustive_q944_gate_host_lifecycle_check() -> Q944GateHostLifecycleReport {
    let mut basis_states_checked = 0usize;
    let mut zero_host_entry_checks = 0usize;
    let mut zero_host_exit_checks = 0usize;
    let mut counter_restoration_checks = 0usize;
    let mut parity_restoration_checks = 0usize;
    let mut operand_restoration_checks = 0usize;
    let mut body_semantic_checks = 0usize;
    let mut width_counts = Vec::new();

    for width in 1..=5 {
        let builder = build_hosted_lifecycle(width);
        let counts = gate_counts(&builder.ops);
        assert_eq!(counts.x, 8 * width + 4);
        assert_eq!(counts.cx, 12 * width);
        assert_eq!(counts.ccx, 12 * width + 3);
        assert_eq!(counts.total, 32 * width + 7);

        let qubits = 2 * width + 5;
        let mask = (1u64 << width) - 1;
        for input in 0..(1u64 << qubits) {
            let host = (input >> 2) & 1;
            if host != 0 {
                continue;
            }
            basis_states_checked += 1;
            zero_host_entry_checks += 1;

            let counter = input & 1;
            let parity = (input >> 1) & 1;
            let body_control = (input >> 3) & 1;
            let body_target = (input >> 4) & 1;
            let v = (input >> 5) & mask;
            let u = (input >> (5 + width)) & mask;
            let expected_toggle = (counter ^ 1) & u64::from(v < u) & body_control;
            let expected = input ^ (expected_toggle << 4);
            let output = apply_scalar(&builder.ops, input);

            assert_eq!(output, expected, "width={width} input={input:#x}");
            zero_host_exit_checks += usize::from(((output >> 2) & 1) == 0);
            counter_restoration_checks += usize::from((output & 1) == counter);
            parity_restoration_checks += usize::from(((output >> 1) & 1) == parity);
            operand_restoration_checks += usize::from(
                ((output >> 5) & mask) == v
                    && ((output >> (5 + width)) & mask) == u,
            );
            body_semantic_checks += usize::from(((output >> 4) & 1) == (body_target ^ expected_toggle));
        }
        width_counts.push((width, counts));
    }

    assert_eq!(zero_host_entry_checks, basis_states_checked);
    assert_eq!(zero_host_exit_checks, basis_states_checked);
    assert_eq!(counter_restoration_checks, basis_states_checked);
    assert_eq!(parity_restoration_checks, basis_states_checked);
    assert_eq!(operand_restoration_checks, basis_states_checked);
    assert_eq!(body_semantic_checks, basis_states_checked);
    Q944GateHostLifecycleReport {
        widths_checked: width_counts.len(),
        basis_states_checked,
        zero_host_entry_checks,
        zero_host_exit_checks,
        counter_restoration_checks,
        parity_restoration_checks,
        operand_restoration_checks,
        body_semantic_checks,
        allocation_free_streams_checked: width_counts.len(),
        phase_clean_streams_checked: width_counts.len(),
        width_counts,
    }
}
