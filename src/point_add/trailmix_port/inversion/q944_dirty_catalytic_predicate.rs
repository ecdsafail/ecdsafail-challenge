//! Arbitrary-dirty catalytic control for elementary involutions.
//!
//! For an involution `G`, desired predicate `f`, and arbitrary dirty bit `d`,
//!
//! `G^d; d ^= f; G^d; d ^= f = G^f`.
//!
//! The circuit proxy below applies this identity gate by gate. Its predicate
//! toggle may be a direct CNOT in the small algebra proof or the independently
//! proved arbitrary-dirty-carry strict comparator in the production-shaped
//! proof. No production PZ route is enabled by this module.

use crate::circuit::{Op, OperationType, QubitId, NO_BIT, NO_QUBIT};
use crate::point_add::trailmix_port::circuit::{Circuit, QReg};
use crate::point_add::trailmix_port::inversion::q944_dirty_parity_microkernels::strict_compare_gated_dirty_carry_refs;
use crate::point_add::B;

#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd)]
pub enum Q944CatalyticKind {
    X,
    CX,
    CCX,
    CCZ,
}

impl Q944CatalyticKind {
    pub const ALL: [Self; 4] = [Self::X, Self::CX, Self::CCX, Self::CCZ];

    pub const fn label(self) -> &'static str {
        match self {
            Self::X => "X",
            Self::CX => "CX",
            Self::CCX => "CCX",
            Self::CCZ => "CCZ",
        }
    }

    pub const fn operand_count(self) -> usize {
        match self {
            Self::X => 1,
            Self::CX => 2,
            Self::CCX | Self::CCZ => 3,
        }
    }
}

#[derive(Clone, Copy)]
pub enum Q944CatalyticPrimitive<'a> {
    X(&'a QReg),
    CX(&'a QReg, &'a QReg),
    CCX(&'a QReg, &'a QReg, &'a QReg),
    CCZ(&'a QReg, &'a QReg, &'a QReg),
}

impl<'a> Q944CatalyticPrimitive<'a> {
    pub const fn kind(self) -> Q944CatalyticKind {
        match self {
            Self::X(_) => Q944CatalyticKind::X,
            Self::CX(_, _) => Q944CatalyticKind::CX,
            Self::CCX(_, _, _) => Q944CatalyticKind::CCX,
            Self::CCZ(_, _, _) => Q944CatalyticKind::CCZ,
        }
    }

    fn operands(self) -> Vec<&'a QReg> {
        match self {
            Self::X(target) => vec![target],
            Self::CX(control, target) => vec![control, target],
            Self::CCX(a, b, target) | Self::CCZ(a, b, target) => vec![a, b, target],
        }
    }
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct Q944CatalyticGateCounts {
    pub x: usize,
    pub cx: usize,
    pub ccx: usize,
    pub ccz: usize,
    pub total: usize,
    pub toffoli_class: usize,
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct Q944CatalyticProxyReport {
    pub primitives_rewritten: usize,
    pub predicate_toggles: usize,
    pub allocation_free_checks: usize,
    pub dirty_restoration_obligations: usize,
    pub lender_restoration_obligations: usize,
    pub emitted: Q944CatalyticGateCounts,
}

pub struct Q944CatalyticProxy<'a, 'q, F>
where
    F: FnMut(&mut Circuit, &QReg),
{
    circ: &'a mut Circuit,
    dirty: &'q QReg,
    lender: &'q QReg,
    toggle_predicate: F,
    report: Q944CatalyticProxyReport,
}

fn assert_distinct(roles: &[(&str, &QReg)]) {
    for (index, (name, lane)) in roles.iter().enumerate() {
        for (other_name, other) in &roles[..index] {
            assert!(
                lane.id() != other.id(),
                "Q944 catalytic alias: {name} aliases {other_name}"
            );
        }
    }
}

fn emit_controlled_by_dirty(
    circ: &mut Circuit,
    primitive: Q944CatalyticPrimitive<'_>,
    dirty: &QReg,
    lender: &QReg,
) {
    match primitive {
        Q944CatalyticPrimitive::X(target) => circ.cx(dirty, target),
        Q944CatalyticPrimitive::CX(control, target) => circ.ccx(dirty, control, target),
        Q944CatalyticPrimitive::CCX(a, b, target) => {
            // Dirty-lender surrounded C^3X. The lender-dependent term occurs
            // twice and cancels; d*a*b reaches the target exactly once.
            circ.ccx(dirty, a, lender);
            circ.ccx(lender, b, target);
            circ.ccx(dirty, a, lender);
            circ.ccx(lender, b, target);
        }
        Q944CatalyticPrimitive::CCZ(a, b, c) => {
            // Phase analogue of the surrounded C^3X construction. The two
            // lender*b*c phase terms cancel, leaving d*a*b*c.
            circ.ccx(dirty, a, lender);
            circ.ccz(lender, b, c);
            circ.ccx(dirty, a, lender);
            circ.ccz(lender, b, c);
        }
    }
}

fn gate_counts_delta(before: [usize; 18], after: [usize; 18]) -> Q944CatalyticGateCounts {
    let x = after[OperationType::X as usize] - before[OperationType::X as usize];
    let cx = after[OperationType::CX as usize] - before[OperationType::CX as usize];
    let ccx = after[OperationType::CCX as usize] - before[OperationType::CCX as usize];
    let ccz = after[OperationType::CCZ as usize] - before[OperationType::CCZ as usize];
    Q944CatalyticGateCounts {
        x,
        cx,
        ccx,
        ccz,
        total: x + cx + ccx + ccz,
        toffoli_class: ccx + ccz,
    }
}

impl<'a, 'q, F> Q944CatalyticProxy<'a, 'q, F>
where
    F: FnMut(&mut Circuit, &QReg),
{
    pub fn new(
        circ: &'a mut Circuit,
        dirty: &'q QReg,
        lender: &'q QReg,
        toggle_predicate: F,
    ) -> Self {
        assert_distinct(&[("dirty", dirty), ("lender", lender)]);
        Self {
            circ,
            dirty,
            lender,
            toggle_predicate,
            report: Q944CatalyticProxyReport::default(),
        }
    }

    pub fn rewrite(&mut self, primitive: Q944CatalyticPrimitive<'_>) {
        let operands = primitive.operands();
        let mut roles = vec![("dirty", self.dirty), ("lender", self.lender)];
        roles.extend(operands.iter().map(|lane| ("operand", *lane)));
        assert_distinct(&roles);

        let allocation_serial = self.circ.b.allocation_serial;
        let next_qubit = self.circ.b.next_qubit;
        let active_qubits = self.circ.b.active_qubits;
        let free_qubits = self.circ.b.free_qubits.clone();
        let before = self.circ.b.counted_kind_ops;

        emit_controlled_by_dirty(self.circ, primitive, self.dirty, self.lender);
        (self.toggle_predicate)(self.circ, self.dirty);
        emit_controlled_by_dirty(self.circ, primitive, self.dirty, self.lender);
        (self.toggle_predicate)(self.circ, self.dirty);

        assert_eq!(self.circ.b.allocation_serial, allocation_serial);
        assert_eq!(self.circ.b.next_qubit, next_qubit);
        assert_eq!(self.circ.b.active_qubits, active_qubits);
        assert_eq!(self.circ.b.free_qubits, free_qubits);
        let delta = gate_counts_delta(before, self.circ.b.counted_kind_ops);
        self.report.primitives_rewritten += 1;
        self.report.predicate_toggles += 2;
        self.report.allocation_free_checks += 1;
        self.report.dirty_restoration_obligations += 1;
        self.report.lender_restoration_obligations += usize::from(matches!(
            primitive,
            Q944CatalyticPrimitive::CCX(_, _, _) | Q944CatalyticPrimitive::CCZ(_, _, _)
        ));
        self.report.emitted.x += delta.x;
        self.report.emitted.cx += delta.cx;
        self.report.emitted.ccx += delta.ccx;
        self.report.emitted.ccz += delta.ccz;
        self.report.emitted.total += delta.total;
        self.report.emitted.toffoli_class += delta.toffoli_class;
    }

    pub const fn report(&self) -> Q944CatalyticProxyReport {
        self.report
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct Q944ClassifiedPrimitive {
    pub kind: Q944CatalyticKind,
    pub operands: [QubitId; 3],
    pub operand_count: usize,
    /// Classical target changed by this involution. `CCZ` has no mutable
    /// target and therefore records `NO_QUBIT`.
    pub mutable_target: QubitId,
}

impl Q944ClassifiedPrimitive {
    pub fn operand_slice(&self) -> &[QubitId] {
        &self.operands[..self.operand_count]
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Q944CatalyticReject {
    ClassicalCondition,
    PredicateIsTarget,
    PredicateControlledCczNeedsCzBase,
    UnsupportedOperation(OperationType),
}

impl Q944CatalyticReject {
    pub const fn label(self) -> &'static str {
        match self {
            Self::ClassicalCondition => "classical-condition",
            Self::PredicateIsTarget => "predicate-is-target",
            Self::PredicateControlledCczNeedsCzBase => {
                "predicate-controlled-ccz-needs-unsupported-cz-base"
            }
            Self::UnsupportedOperation(_) => "unsupported-operation",
        }
    }
}

fn classified(
    kind: Q944CatalyticKind,
    operands: &[QubitId],
    mutable_target: QubitId,
) -> Q944ClassifiedPrimitive {
    let mut out = [NO_QUBIT; 3];
    out[..operands.len()].copy_from_slice(operands);
    Q944ClassifiedPrimitive {
        kind,
        operands: out,
        operand_count: operands.len(),
        mutable_target,
    }
}

/// Convert one operation from the `active=1` body template into a supported
/// elementary base involution. The formal predicate may occur only as a
/// control; operations not mentioning it are still controlled catalytically.
pub fn q944_classify_template_op(
    op: &Op,
    formal_predicate: QubitId,
) -> Result<Q944ClassifiedPrimitive, Q944CatalyticReject> {
    if op.c_condition != NO_BIT {
        return Err(Q944CatalyticReject::ClassicalCondition);
    }
    match op.kind {
        OperationType::X => {
            if op.q_target == formal_predicate {
                Err(Q944CatalyticReject::PredicateIsTarget)
            } else {
                Ok(classified(
                    Q944CatalyticKind::X,
                    &[op.q_target],
                    op.q_target,
                ))
            }
        }
        OperationType::CX => {
            if op.q_target == formal_predicate {
                Err(Q944CatalyticReject::PredicateIsTarget)
            } else if op.q_control1 == formal_predicate {
                Ok(classified(
                    Q944CatalyticKind::X,
                    &[op.q_target],
                    op.q_target,
                ))
            } else {
                Ok(classified(
                    Q944CatalyticKind::CX,
                    &[op.q_control1, op.q_target],
                    op.q_target,
                ))
            }
        }
        OperationType::CCX => {
            if op.q_target == formal_predicate {
                return Err(Q944CatalyticReject::PredicateIsTarget);
            }
            if op.q_control1 == formal_predicate {
                Ok(classified(
                    Q944CatalyticKind::CX,
                    &[op.q_control2, op.q_target],
                    op.q_target,
                ))
            } else if op.q_control2 == formal_predicate {
                Ok(classified(
                    Q944CatalyticKind::CX,
                    &[op.q_control1, op.q_target],
                    op.q_target,
                ))
            } else {
                Ok(classified(
                    Q944CatalyticKind::CCX,
                    &[op.q_control2, op.q_control1, op.q_target],
                    op.q_target,
                ))
            }
        }
        OperationType::CCZ => {
            if [op.q_control2, op.q_control1, op.q_target].contains(&formal_predicate) {
                Err(Q944CatalyticReject::PredicateControlledCczNeedsCzBase)
            } else {
                Ok(classified(
                    Q944CatalyticKind::CCZ,
                    &[op.q_control2, op.q_control1, op.q_target],
                    NO_QUBIT,
                ))
            }
        }
        other => Err(Q944CatalyticReject::UnsupportedOperation(other)),
    }
}

/// Choose a catalytic dirty bit and a second restored dirty lender. Both are
/// long-lived candidates and must be disjoint from the primitive and every
/// predicate-comparator lane supplied in `forbidden`.
pub fn q944_select_catalytic_dirty_pair(
    primitive: &Q944ClassifiedPrimitive,
    candidates: &[QubitId],
    forbidden: &[QubitId],
) -> Option<(QubitId, QubitId)> {
    let usable = |candidate: QubitId| {
        candidate != NO_QUBIT
            && !primitive.operand_slice().contains(&candidate)
            && !forbidden.contains(&candidate)
    };
    for &dirty in candidates {
        if !usable(dirty) {
            continue;
        }
        for &lender in candidates {
            if lender != dirty && usable(lender) {
                return Some((dirty, lender));
            }
        }
    }
    None
}

/// Exact per-primitive cost when `f = !done AND (v<u)` is toggled by the
/// corrected dirty-carry comparator at full width `n`.
#[must_use]
pub const fn q944_catalytic_cost(width: usize, kind: Q944CatalyticKind) -> Q944CatalyticGateCounts {
    let x = 8 * width + 4;
    let mut cx = 12 * width;
    let mut ccx = 12 * width + 2;
    let mut ccz = 0;
    match kind {
        Q944CatalyticKind::X => cx += 2,
        Q944CatalyticKind::CX => ccx += 2,
        Q944CatalyticKind::CCX => ccx += 8,
        Q944CatalyticKind::CCZ => {
            ccx += 4;
            ccz += 4;
        }
    }
    Q944CatalyticGateCounts {
        x,
        cx,
        ccx,
        ccz,
        total: x + cx + ccx + ccz,
        toffoli_class: ccx + ccz,
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Q944CatalyticProofReport {
    pub algebra_cases_checked: usize,
    pub direct_kinds_checked: usize,
    pub direct_basis_states_checked: usize,
    pub comparator_widths_checked: usize,
    pub comparator_kind_width_pairs_checked: usize,
    pub comparator_basis_states_checked: usize,
    pub classical_output_checks: usize,
    pub phase_checks: usize,
    pub dirty_restoration_checks: usize,
    pub lender_restoration_checks: usize,
    pub operand_restoration_checks: usize,
    pub allocation_free_streams_checked: usize,
    pub phase_clean_classical_streams_checked: usize,
    pub phase_sensitive_streams_checked: usize,
    pub direct_counts: Vec<(Q944CatalyticKind, Q944CatalyticGateCounts)>,
    pub comparator_counts: Vec<(usize, Q944CatalyticKind, Q944CatalyticGateCounts)>,
}

fn apply_scalar_phase(ops: &[Op], mut state: u64) -> (u64, bool) {
    let bit = |word: u64, id: QubitId| ((word >> id.0) & 1) != 0;
    let mut phase = false;
    for op in ops {
        match op.kind {
            OperationType::X => state ^= 1u64 << op.q_target.0,
            OperationType::CX => {
                if bit(state, op.q_control1) {
                    state ^= 1u64 << op.q_target.0;
                }
            }
            OperationType::CCX => {
                if bit(state, op.q_control1) && bit(state, op.q_control2) {
                    state ^= 1u64 << op.q_target.0;
                }
            }
            OperationType::CCZ => {
                if bit(state, op.q_control1) && bit(state, op.q_control2) && bit(state, op.q_target)
                {
                    phase = !phase;
                }
            }
            other => panic!("Q944 catalytic proof saw unsupported gate {other:?}"),
        }
    }
    (state, phase)
}

fn total_gate_counts(builder: &B) -> Q944CatalyticGateCounts {
    gate_counts_delta([0; 18], builder.counted_kind_ops)
}

struct DirectHarness {
    builder: B,
    proxy: Q944CatalyticProxyReport,
}

fn build_direct_harness(kind: Q944CatalyticKind) -> DirectHarness {
    let mut circ = Circuit::new();
    let predicate = circ.alloc_qreg("q944.catalytic.f");
    let dirty = circ.alloc_qreg("q944.catalytic.d");
    let lender = circ.alloc_qreg("q944.catalytic.psi");
    let lanes = circ.alloc_qreg_bits("q944.catalytic.g", 3);
    let primitive = match kind {
        Q944CatalyticKind::X => Q944CatalyticPrimitive::X(&lanes[0]),
        Q944CatalyticKind::CX => Q944CatalyticPrimitive::CX(&lanes[0], &lanes[1]),
        Q944CatalyticKind::CCX => Q944CatalyticPrimitive::CCX(&lanes[0], &lanes[1], &lanes[2]),
        Q944CatalyticKind::CCZ => Q944CatalyticPrimitive::CCZ(&lanes[0], &lanes[1], &lanes[2]),
    };
    let mut proxy = Q944CatalyticProxy::new(&mut circ, &dirty, &lender, |circ, dirty| {
        circ.cx(&predicate, dirty);
    });
    proxy.rewrite(primitive);
    let proxy_report = proxy.report();
    drop(proxy);
    let builder = circ.into_builder();
    assert_eq!(builder.next_qubit, 6);
    assert_eq!(builder.active_qubits, 6);
    assert_eq!(builder.peak_qubits, 6);
    drop((predicate, dirty, lender, lanes));
    DirectHarness {
        builder,
        proxy: proxy_report,
    }
}

struct ComparatorHarness {
    builder: B,
    proxy: Q944CatalyticProxyReport,
}

fn build_comparator_harness(width: usize, kind: Q944CatalyticKind) -> ComparatorHarness {
    let mut circ = Circuit::new();
    let done = circ.alloc_qreg("q944.catalytic.done");
    let dirty = circ.alloc_qreg("q944.catalytic.d");
    let parity = circ.alloc_qreg("q944.catalytic.parity");
    let lender = circ.alloc_qreg("q944.catalytic.psi");
    let v = circ.alloc_qreg_bits("q944.catalytic.v", width);
    let u = circ.alloc_qreg_bits("q944.catalytic.u", width);
    let lanes = circ.alloc_qreg_bits("q944.catalytic.g", 3);
    let vr: Vec<&QReg> = v.iter().collect();
    let ur: Vec<&QReg> = u.iter().collect();
    let primitive = match kind {
        Q944CatalyticKind::X => Q944CatalyticPrimitive::X(&lanes[0]),
        Q944CatalyticKind::CX => Q944CatalyticPrimitive::CX(&lanes[0], &lanes[1]),
        Q944CatalyticKind::CCX => Q944CatalyticPrimitive::CCX(&lanes[0], &lanes[1], &lanes[2]),
        Q944CatalyticKind::CCZ => Q944CatalyticPrimitive::CCZ(&lanes[0], &lanes[1], &lanes[2]),
    };
    let mut proxy = Q944CatalyticProxy::new(&mut circ, &dirty, &lender, |circ, dirty| {
        circ.x(&done);
        strict_compare_gated_dirty_carry_refs(circ, &vr, &ur, &done, dirty, &parity);
        circ.x(&done);
    });
    proxy.rewrite(primitive);
    let proxy_report = proxy.report();
    drop(proxy);
    let inputs = 2 * width + 7;
    let builder = circ.into_builder();
    assert_eq!(builder.next_qubit as usize, inputs);
    assert_eq!(builder.active_qubits as usize, inputs);
    assert_eq!(builder.peak_qubits as usize, inputs);
    drop((vr, ur));
    drop((done, dirty, parity, lender, v, u, lanes));
    ComparatorHarness {
        builder,
        proxy: proxy_report,
    }
}

fn expected_effect(
    kind: Q944CatalyticKind,
    state: u64,
    predicate: bool,
    lane_base: usize,
) -> (u64, bool) {
    let bit = |index: usize| ((state >> index) & 1) != 0;
    let a = bit(lane_base);
    let b = bit(lane_base + 1);
    let c = bit(lane_base + 2);
    let mut output = state;
    let mut phase = false;
    if predicate {
        match kind {
            Q944CatalyticKind::X => output ^= 1u64 << lane_base,
            Q944CatalyticKind::CX if a => output ^= 1u64 << (lane_base + 1),
            Q944CatalyticKind::CCX if a && b => output ^= 1u64 << (lane_base + 2),
            Q944CatalyticKind::CCZ if a && b && c => phase = true,
            _ => {}
        }
    }
    (output, phase)
}

/// Exhaustive algebra, classical-output, phase, dirty-lane, lender, operand,
/// and allocation checks for direct and comparator-derived predicates.
#[must_use]
pub fn exhaustive_q944_dirty_catalytic_predicate_check() -> Q944CatalyticProofReport {
    let mut algebra_cases_checked = 0;
    for f in [false, true] {
        for d in [false, true] {
            assert_eq!(d ^ (d ^ f), f);
            assert_eq!(d ^ f ^ f, d);
            algebra_cases_checked += 2;
        }
    }

    let mut direct_basis_states_checked = 0usize;
    let mut comparator_basis_states_checked = 0usize;
    let mut classical_output_checks = 0usize;
    let mut phase_checks = 0usize;
    let mut dirty_restoration_checks = 0usize;
    let mut lender_restoration_checks = 0usize;
    let mut operand_restoration_checks = 0usize;
    let mut direct_counts = Vec::new();
    let mut comparator_counts = Vec::new();

    for kind in Q944CatalyticKind::ALL {
        let harness = build_direct_harness(kind);
        assert_eq!(harness.proxy.primitives_rewritten, 1);
        assert_eq!(harness.proxy.predicate_toggles, 2);
        for state in 0..(1u64 << 6) {
            let predicate = state & 1 != 0;
            let expected = expected_effect(kind, state, predicate, 3);
            let actual = apply_scalar_phase(&harness.builder.ops, state);
            assert_eq!(actual, expected, "direct kind={kind:?} state={state:#x}");
            direct_basis_states_checked += 1;
            classical_output_checks += 1;
            phase_checks += 1;
            dirty_restoration_checks += usize::from(((actual.0 >> 1) & 1) == ((state >> 1) & 1));
            lender_restoration_checks += usize::from(((actual.0 >> 2) & 1) == ((state >> 2) & 1));
            // The full `(state, phase)` equality above checks controls,
            // targets, unused lanes, and both dirty lanes simultaneously.
            operand_restoration_checks += 1;
        }
        direct_counts.push((kind, total_gate_counts(&harness.builder)));
    }

    for width in 1..=3 {
        for kind in Q944CatalyticKind::ALL {
            let harness = build_comparator_harness(width, kind);
            assert_eq!(harness.proxy.primitives_rewritten, 1);
            assert_eq!(harness.proxy.predicate_toggles, 2);
            let inputs = 2 * width + 7;
            let mask = (1u64 << width) - 1;
            let lane_base = 4 + 2 * width;
            for state in 0..(1u64 << inputs) {
                let done = state & 1 != 0;
                let v = (state >> 4) & mask;
                let u = (state >> (4 + width)) & mask;
                let predicate = !done && v < u;
                let expected = expected_effect(kind, state, predicate, lane_base);
                let actual = apply_scalar_phase(&harness.builder.ops, state);
                assert_eq!(
                    actual, expected,
                    "comparator width={width} kind={kind:?} state={state:#x}"
                );
                comparator_basis_states_checked += 1;
                classical_output_checks += 1;
                phase_checks += 1;
                dirty_restoration_checks +=
                    usize::from(((actual.0 >> 1) & 1) == ((state >> 1) & 1));
                lender_restoration_checks +=
                    usize::from(((actual.0 >> 3) & 1) == ((state >> 3) & 1));
                operand_restoration_checks += 1;
            }
            let counts = total_gate_counts(&harness.builder);
            assert_eq!(counts, q944_catalytic_cost(width, kind));
            comparator_counts.push((width, kind, counts));
        }
    }

    let all_states = direct_basis_states_checked + comparator_basis_states_checked;
    assert_eq!(classical_output_checks, all_states);
    assert_eq!(phase_checks, all_states);
    assert_eq!(dirty_restoration_checks, all_states);
    assert_eq!(lender_restoration_checks, all_states);
    assert_eq!(operand_restoration_checks, all_states);
    Q944CatalyticProofReport {
        algebra_cases_checked,
        direct_kinds_checked: direct_counts.len(),
        direct_basis_states_checked,
        comparator_widths_checked: 3,
        comparator_kind_width_pairs_checked: comparator_counts.len(),
        comparator_basis_states_checked,
        classical_output_checks,
        phase_checks,
        dirty_restoration_checks,
        lender_restoration_checks,
        operand_restoration_checks,
        allocation_free_streams_checked: direct_counts.len() + comparator_counts.len(),
        phase_clean_classical_streams_checked: 3 * 4,
        phase_sensitive_streams_checked: 4,
        direct_counts,
        comparator_counts,
    }
}
