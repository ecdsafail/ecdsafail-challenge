//! Allocation-free quotient witness handoff for the five blocked Q944 rows.
//!
//! The 25-bit quotient register is zero before division. Its top lane `q[24]`
//! can therefore host the outer division predicate while the shift is built.
//! Every quotient bit except the sentinel is deposited before subtraction.
//! The deposited one-hot bit clears the binary shift inside the hosted body.
//! Only `s=24` remains parked in high shift lanes while the outer lifecycle
//! clears `q[24]`, after which that sentinel moves back into the quotient.
//! Reverse execution is the exact gate inverse, ordered so `s[0]` and `s[1]`
//! remain clean for the outer comparator.

use crate::circuit::{Op, OperationType};
use crate::point_add::trailmix_port::arith::mcx::mcx_dirty_ladder;
use crate::point_add::trailmix_port::circuit::{Circuit, QReg};
use crate::point_add::B;

pub const Q944_QUOTIENT_WIDTH: usize = 25;
pub const Q944_QUOTIENT_SENTINEL: usize = 24;
pub const Q944_SHIFT_WIDTH: usize = 5;

fn assert_unique(roles: &[(&str, &QReg)]) {
    for (index, (name, lane)) in roles.iter().enumerate() {
        for (other_name, other) in &roles[..index] {
            assert!(
                lane.id() != other.id(),
                "Q944 quotient witness alias: {name} aliases {other_name}"
            );
        }
    }
}

fn assert_layout(q: &[QReg], s: &[&QReg], active: &QReg, lenders: &[&QReg]) {
    assert_eq!(q.len(), Q944_QUOTIENT_WIDTH);
    assert_eq!(s.len(), Q944_SHIFT_WIDTH);
    assert_eq!(active.id(), q[Q944_QUOTIENT_SENTINEL].id());
    assert!(lenders.len() >= Q944_SHIFT_WIDTH - 1);
    let mut roles = Vec::with_capacity(q.len() + s.len() + lenders.len());
    roles.extend(q.iter().map(|lane| ("q", lane)));
    roles.extend(s.iter().map(|lane| ("s", *lane)));
    roles.extend(lenders.iter().map(|lane| ("lender", *lane)));
    assert_unique(&roles);
}

/// Deposit `q[s] ^= active`, except when `s=24`, where the active host already
/// occupies the correct quotient lane. Every dirty lender is restored.
pub fn q944_partial_demux_excluding_sentinel(
    circ: &mut Circuit,
    q: &[QReg],
    s: &[&QReg],
    active: &QReg,
    lenders: &[&QReg],
) {
    assert_layout(q, s, active, lenders);
    let allocation_serial = circ.b.allocation_serial;
    let next_qubit = circ.b.next_qubit;
    let active_qubits = circ.b.active_qubits;
    let free_qubits = circ.b.free_qubits.clone();
    let previous = circ.push_section("q944.qw.partial-demux");
    for (index, target) in q.iter().enumerate() {
        if index == Q944_QUOTIENT_SENTINEL {
            continue;
        }
        for (bit, lane) in s.iter().enumerate() {
            if (index >> bit) & 1 == 0 {
                circ.x(lane);
            }
        }
        let mut controls = Vec::with_capacity(1 + s.len());
        controls.push(active);
        controls.extend(s.iter().copied());
        mcx_dirty_ladder(circ, &controls, target, &lenders[..controls.len() - 2]);
        for (bit, lane) in s.iter().enumerate().rev() {
            if (index >> bit) & 1 == 0 {
                circ.x(lane);
            }
        }
    }
    circ.pop_section(&previous);
    assert_eq!(circ.b.allocation_serial, allocation_serial);
    assert_eq!(circ.b.next_qubit, next_qubit);
    assert_eq!(circ.b.active_qubits, active_qubits);
    assert_eq!(circ.b.free_qubits, free_qubits);
}

/// Clear the binary shift from the deposited non-sentinel one-hot quotient.
/// For `s=24`, this deliberately leaves `s[4:3]=11` and `s[2:0]=000`: the
/// outer predicate lifecycle may therefore reuse `s[0]` and `s[1]` while the
/// sentinel remains parked in disjoint high lanes.
pub fn q944_clear_non_sentinel_shift(circ: &mut Circuit, q: &[QReg], s: &[&QReg]) {
    assert_eq!(q.len(), Q944_QUOTIENT_WIDTH);
    assert_eq!(s.len(), Q944_SHIFT_WIDTH);
    let previous = circ.push_section("q944.qw.clear-non-sentinel-shift");
    for (index, source) in q.iter().enumerate() {
        if index == Q944_QUOTIENT_SENTINEL {
            continue;
        }
        for (bit, target) in s.iter().enumerate() {
            if (index >> bit) & 1 == 1 {
                circ.cx(source, target);
            }
        }
    }
    circ.pop_section(&previous);
}

/// After the outer predicate has cleared `q[24]`, move the parked `s=24`
/// sentinel from `s[4:3]=11` into `q[24]` and clear both high shift lanes.
pub fn q944_commit_parked_sentinel(circ: &mut Circuit, q: &[QReg], s: &[&QReg]) {
    assert_eq!(q.len(), Q944_QUOTIENT_WIDTH);
    assert_eq!(s.len(), Q944_SHIFT_WIDTH);
    let previous = circ.push_section("q944.qw.commit-parked-sentinel");
    circ.cx(s[4], s[3]);
    circ.cx(s[4], &q[Q944_QUOTIENT_SENTINEL]);
    circ.cx(&q[Q944_QUOTIENT_SENTINEL], s[4]);
    circ.pop_section(&previous);
}

/// Complete the forward one-hot handoff when no outer lifecycle needs the low
/// shift lanes between the two stages.
pub fn q944_forward_finalize_quotient(circ: &mut Circuit, q: &[QReg], s: &[&QReg]) {
    q944_clear_non_sentinel_shift(circ, q, s);
    q944_commit_parked_sentinel(circ, q, s);
}

/// First part of the exact inverse of `q944_forward_finalize_quotient`.
/// It consumes a possible `q[24]` witness into `s[4:3]` while deliberately
/// leaving `s[0:2]=0` for `gate_hold_counter_zero` scratch.
pub fn q944_reverse_park_sentinel(circ: &mut Circuit, q: &[QReg], s: &[&QReg]) {
    assert_eq!(q.len(), Q944_QUOTIENT_WIDTH);
    assert_eq!(s.len(), Q944_SHIFT_WIDTH);
    let previous = circ.push_section("q944.qw.reverse-park-sentinel");
    circ.cx(&q[Q944_QUOTIENT_SENTINEL], s[4]);
    circ.cx(s[4], &q[Q944_QUOTIENT_SENTINEL]);
    circ.cx(s[4], s[3]);
    circ.pop_section(&previous);
}

/// Complete the reverse quotient-to-shift decode after the outer predicate has
/// been computed into the now-clean `q[24]` host.
pub fn q944_reverse_materialize_non_sentinel_index(
    circ: &mut Circuit,
    q: &[QReg],
    s: &[&QReg],
) {
    assert_eq!(q.len(), Q944_QUOTIENT_WIDTH);
    assert_eq!(s.len(), Q944_SHIFT_WIDTH);
    let previous = circ.push_section("q944.qw.reverse-materialize-index");
    for (index, source) in q.iter().enumerate() {
        if index == Q944_QUOTIENT_SENTINEL {
            continue;
        }
        for (bit, target) in s.iter().enumerate() {
            if (index >> bit) & 1 == 1 {
                circ.cx(source, target);
            }
        }
    }
    circ.pop_section(&previous);
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct Q944QuotientWitnessGateCounts {
    pub x: usize,
    pub cx: usize,
    pub ccx: usize,
    pub total: usize,
    pub toffoli_class: usize,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Q944QuotientWitnessProofReport {
    pub quotient_width: usize,
    pub shift_width: usize,
    pub sentinel: usize,
    pub classical_equivalence_cases: usize,
    pub circuit_forward_cases: usize,
    pub circuit_roundtrip_cases: usize,
    pub gate_scratch_boundary_checks: usize,
    pub phase_checks: usize,
    pub predicate_restoration_checks: usize,
    pub lender_restoration_checks: usize,
    pub allocation_free_streams_checked: usize,
    pub forward_counts: Q944QuotientWitnessGateCounts,
    pub roundtrip_counts: Q944QuotientWitnessGateCounts,
}

fn counts(builder: &B, end: usize) -> Q944QuotientWitnessGateCounts {
    let mut x = 0;
    let mut cx = 0;
    let mut ccx = 0;
    for op in &builder.ops[..end] {
        match op.kind {
            OperationType::X => x += 1,
            OperationType::CX => cx += 1,
            OperationType::CCX => ccx += 1,
            other => panic!("Q944 quotient witness emitted unsupported gate {other:?}"),
        }
    }
    Q944QuotientWitnessGateCounts {
        x,
        cx,
        ccx,
        total: x + cx + ccx,
        toffoli_class: ccx,
    }
}

fn apply_scalar(ops: &[Op], mut state: u64) -> (u64, bool) {
    let bit = |word: u64, id: u64| ((word >> id) & 1) != 0;
    let mut phase = false;
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
            OperationType::CCZ => {
                if bit(state, op.q_control1.0)
                    && bit(state, op.q_control2.0)
                    && bit(state, op.q_target.0)
                {
                    phase = !phase;
                }
            }
            other => panic!("Q944 quotient witness scalar saw unsupported gate {other:?}"),
        }
    }
    (state, phase)
}

struct Harness {
    builder: B,
    predicate: usize,
    q: Vec<usize>,
    s: Vec<usize>,
    lenders: Vec<usize>,
    target: usize,
    hosted_body_end: usize,
    forward_end: usize,
}

fn build_harness() -> Harness {
    let mut circ = Circuit::new();
    let predicate = circ.alloc_qreg("q944.qw.predicate");
    let q = circ.alloc_qreg_bits("q944.qw.q", Q944_QUOTIENT_WIDTH);
    let s = circ.alloc_qreg_bits("q944.qw.s", Q944_SHIFT_WIDTH);
    let lenders = circ.alloc_qreg_bits("q944.qw.lender", Q944_SHIFT_WIDTH - 1);
    let target = circ.alloc_qreg("q944.qw.target");
    let s_refs: Vec<&QReg> = s.iter().collect();
    let lender_refs: Vec<&QReg> = lenders.iter().collect();
    let active = &q[Q944_QUOTIENT_SENTINEL];

    circ.cx(&predicate, active);
    q944_partial_demux_excluding_sentinel(&mut circ, &q, &s_refs, active, &lender_refs);
    circ.cx(active, &target);
    q944_clear_non_sentinel_shift(&mut circ, &q, &s_refs);
    let hosted_body_end = circ.b.ops.len();
    circ.cx(&predicate, active);
    q944_commit_parked_sentinel(&mut circ, &q, &s_refs);
    let forward_end = circ.b.ops.len();

    q944_reverse_park_sentinel(&mut circ, &q, &s_refs);
    circ.cx(&predicate, active);
    q944_reverse_materialize_non_sentinel_index(&mut circ, &q, &s_refs);
    circ.cx(active, &target);
    q944_partial_demux_excluding_sentinel(&mut circ, &q, &s_refs, active, &lender_refs);
    circ.cx(&predicate, active);

    let expected_qubits = 1 + Q944_QUOTIENT_WIDTH + Q944_SHIFT_WIDTH
        + (Q944_SHIFT_WIDTH - 1)
        + 1;
    let builder = circ.into_builder();
    assert_eq!(builder.next_qubit as usize, expected_qubits);
    assert_eq!(builder.active_qubits as usize, expected_qubits);
    assert_eq!(builder.peak_qubits as usize, expected_qubits);
    let report = Harness {
        builder,
        predicate: predicate.id() as usize,
        q: q.iter().map(|lane| lane.id() as usize).collect(),
        s: s.iter().map(|lane| lane.id() as usize).collect(),
        lenders: lenders.iter().map(|lane| lane.id() as usize).collect(),
        target: target.id() as usize,
        hosted_body_end,
        forward_end,
    };
    drop((s_refs, lender_refs));
    drop((predicate, q, s, lenders, target));
    report
}

fn bit(state: u64, index: usize) -> bool {
    ((state >> index) & 1) != 0
}

/// Exhaustive abstract arithmetic equivalence and exact 25-bit circuit proof.
#[must_use]
pub fn exhaustive_q944_quotient_witness_check() -> Q944QuotientWitnessProofReport {
    let mut classical_equivalence_cases = 0usize;
    for width in 1..=5 {
        let modulus = 1usize << width;
        for a in 0..modulus {
            for b in 0..modulus {
                for predicate in 0..=1usize {
                    for shift in 0..Q944_QUOTIENT_WIDTH {
                        let baseline_a = (a + modulus - predicate * b) % modulus;
                        let baseline_q = predicate << shift;

                        // Candidate order: retain f in q[24], deposit every
                        // non-sentinel q bit, subtract once, clear the binary
                        // shift from the one-hot witness, clear f, then commit
                        // the parked sentinel. Production guarantees s=0 when
                        // f=0 because the bit-length difference is f-masked.
                        let mut candidate_q = 0usize;
                        let mut candidate_shift = predicate * shift;
                        let mut candidate_host = predicate;
                        if candidate_host != 0 && shift != Q944_QUOTIENT_SENTINEL {
                            candidate_q ^= 1usize << shift;
                        }
                        let candidate_a = (a + modulus - predicate * b) % modulus;
                        if candidate_host != 0 && shift != Q944_QUOTIENT_SENTINEL {
                            candidate_shift ^= shift;
                        }
                        assert_eq!(candidate_shift & 0b11, 0);
                        candidate_host ^= predicate;
                        if candidate_shift == Q944_QUOTIENT_SENTINEL {
                            assert_eq!(candidate_host, 0);
                            candidate_q ^= 1usize << Q944_QUOTIENT_SENTINEL;
                            candidate_shift = 0;
                        }
                        assert_eq!((candidate_a, candidate_q), (baseline_a, baseline_q));
                        assert_eq!(candidate_shift, 0);
                        assert_eq!(candidate_host, 0);
                        let restored = (candidate_a + predicate * b) % modulus;
                        assert_eq!(restored, a);
                        classical_equivalence_cases += 1;
                    }
                }
            }
        }
    }

    let harness = build_harness();
    let mut circuit_forward_cases = 0usize;
    let mut circuit_roundtrip_cases = 0usize;
    let mut gate_scratch_boundary_checks = 0usize;
    let mut phase_checks = 0usize;
    let mut predicate_restoration_checks = 0usize;
    let mut lender_restoration_checks = 0usize;
    for predicate in 0..=1u64 {
        for shift in 0..Q944_QUOTIENT_WIDTH {
            for target in 0..=1u64 {
                for lender_word in 0..(1u64 << harness.lenders.len()) {
                    let mut initial = predicate << harness.predicate;
                    initial |= target << harness.target;
                    let effective_shift = if predicate == 0 { 0 } else { shift };
                    for (bit_index, lane) in harness.s.iter().enumerate() {
                        initial |= (((effective_shift >> bit_index) & 1) as u64) << lane;
                    }
                    for (bit_index, lane) in harness.lenders.iter().enumerate() {
                        initial |= ((lender_word >> bit_index) & 1) << lane;
                    }

                    let (hosted_boundary, hosted_boundary_phase) = apply_scalar(
                        &harness.builder.ops[..harness.hosted_body_end],
                        initial,
                    );
                    assert!(!hosted_boundary_phase);
                    for (shift_bit, &lane) in harness.s.iter().enumerate() {
                        let expected = predicate != 0
                            && shift == Q944_QUOTIENT_SENTINEL
                            && (shift_bit == 3 || shift_bit == 4);
                        assert_eq!(bit(hosted_boundary, lane), expected);
                    }
                    assert_eq!(
                        bit(hosted_boundary, harness.q[Q944_QUOTIENT_SENTINEL]),
                        predicate != 0,
                    );
                    gate_scratch_boundary_checks += 1;

                    let (forward, forward_phase) =
                        apply_scalar(&harness.builder.ops[..harness.forward_end], initial);
                    let mut expected_forward = initial;
                    for lane in &harness.s {
                        expected_forward &= !(1u64 << lane);
                    }
                    expected_forward ^= predicate << harness.target;
                    if predicate != 0 {
                        expected_forward |= 1u64 << harness.q[shift];
                    }
                    assert_eq!(forward, expected_forward);
                    assert!(!forward_phase);
                    circuit_forward_cases += 1;
                    phase_checks += 1;
                    predicate_restoration_checks +=
                        usize::from(bit(forward, harness.predicate) == (predicate != 0));
                    lender_restoration_checks += usize::from(
                        harness
                            .lenders
                            .iter()
                            .all(|&lane| bit(forward, lane) == bit(initial, lane)),
                    );

                    let (roundtrip, roundtrip_phase) =
                        apply_scalar(&harness.builder.ops[harness.forward_end..], forward);
                    assert_eq!(roundtrip, initial);
                    assert!(!roundtrip_phase);
                    circuit_roundtrip_cases += 1;
                    phase_checks += 1;
                    predicate_restoration_checks +=
                        usize::from(bit(roundtrip, harness.predicate) == (predicate != 0));
                    lender_restoration_checks += usize::from(
                        harness
                            .lenders
                            .iter()
                            .all(|&lane| bit(roundtrip, lane) == bit(initial, lane)),
                    );
                }
            }
        }
    }
    assert_eq!(predicate_restoration_checks, 2 * circuit_forward_cases);
    assert_eq!(lender_restoration_checks, 2 * circuit_forward_cases);
    Q944QuotientWitnessProofReport {
        quotient_width: Q944_QUOTIENT_WIDTH,
        shift_width: Q944_SHIFT_WIDTH,
        sentinel: Q944_QUOTIENT_SENTINEL,
        classical_equivalence_cases,
        circuit_forward_cases,
        circuit_roundtrip_cases,
        gate_scratch_boundary_checks,
        phase_checks,
        predicate_restoration_checks,
        lender_restoration_checks,
        allocation_free_streams_checked: 2,
        forward_counts: counts(&harness.builder, harness.forward_end),
        roundtrip_counts: counts(&harness.builder, harness.builder.ops.len()),
    }
}

#[must_use]
pub const fn q944_dirty_arithmetic_toffoli(width: usize) -> usize {
    2 * width * width + width + 1
}

#[must_use]
pub const fn q944_distributed_arithmetic_toffoli(width: usize) -> usize {
    Q944_QUOTIENT_WIDTH * q944_dirty_arithmetic_toffoli(width)
}
