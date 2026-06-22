//! Reversible gate port of the public register-sharing EEA reference.
//!
//! This module is intentionally staged. The arithmetic primitives and their
//! reduced-width exhaustive proofs land before the complete 1,479-step circuit.

use crate::point_add::trailmix_port::circuit::{Circuit, QReg};

use super::register_shared_eea_microkernels::{
    apply_scalar, controlled_add_one, controlled_decrement_mod_2n, controlled_increment_mod_2n,
    controlled_sub_one, controlled_swap_registers, decrement_mod_2n, gate_counts, increment_mod_2n,
    multi_controlled_x_vchain, variable_rotate_high, variable_rotate_high_refs,
    variable_rotate_low, variable_rotate_low_refs, RegisterSharedGateCounts,
};
use crate::point_add::B;

pub const REFERENCE_LENGTH_WIDTH: usize = 9;
pub const REFERENCE_STEPS: usize = 1_479;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct ReferenceActiveWindows {
    pub r_add_sub: (usize, usize),
    pub quotient_swap: (usize, usize),
    pub t_add_sub: (usize, usize),
    pub length_update_t: (usize, usize),
    pub length_update_r: (usize, usize),
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ReferenceCuccaroProofReport {
    pub widths_checked: usize,
    pub basis_states_checked: usize,
    pub carry_clean_checks: usize,
    pub inverse_pair_checks: usize,
    pub add9: RegisterSharedGateCounts,
    pub sub9: RegisterSharedGateCounts,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ReferenceLocationSwapProofReport {
    pub work_widths_checked: usize,
    pub basis_states_checked: usize,
    pub scratch_clean_checks: usize,
    pub inverse_pair_checks: usize,
    pub full259_length9: RegisterSharedGateCounts,
    pub full259_length9_inverse: RegisterSharedGateCounts,
    pub reference_steps: usize,
    pub full_window_toffoli_upper_bound: usize,
    pub scheduled_window_sum: usize,
    pub scheduled_toffoli: usize,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ReferenceScheduleProofReport {
    pub steps_checked: usize,
    pub r_window_sum: usize,
    pub quotient_swap_window_sum: usize,
    pub t_window_sum: usize,
    pub length_t_window_sum: usize,
    pub length_r_window_sum: usize,
    pub maximum_r_window: usize,
    pub maximum_quotient_swap_window: usize,
    pub maximum_t_window: usize,
    pub maximum_length_t_window: usize,
    pub maximum_length_r_window: usize,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct ReferenceInitializerProofReport {
    pub cases_checked: usize,
    pub reflected_cases_checked: usize,
    pub non_reflected_cases_checked: usize,
    pub input_qubits: usize,
    pub initialization_transient_peak_qubits: usize,
    pub packed_peak_qubits: usize,
    pub packed_active_qubits: usize,
    pub final_active_qubits: usize,
    pub emitted_ops: usize,
    pub emitted_toffoli: usize,
    pub emitted_hmr: usize,
    pub emitted_resets: usize,
    pub classical_roundtrip_checks: usize,
    pub phase_cleanup_checks: usize,
    pub ancilla_cleanup_checks: usize,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ReferenceCoefficientArithmeticProofReport {
    pub work_widths_checked: usize,
    pub basis_states_checked: usize,
    pub inverse_pair_checks: usize,
    pub scratch_clean_checks: usize,
    pub control_off_identity_checks: usize,
    pub length_restore_checks: usize,
    pub t_add257: RegisterSharedGateCounts,
    pub t_sub257: RegisterSharedGateCounts,
    pub reference_steps: usize,
    pub scheduled_window_sum: usize,
    pub scheduled_add_toffoli: usize,
    pub scheduled_sub_toffoli: usize,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ReferenceRemainderArithmeticProofReport {
    pub work_widths_checked: usize,
    pub basis_states_checked: usize,
    pub inverse_pair_checks: usize,
    pub scratch_clean_checks: usize,
    pub control_off_identity_checks: usize,
    pub zero_remainder_identity_checks: usize,
    pub length_restore_checks: usize,
    pub r_add257: RegisterSharedGateCounts,
    pub r_sub257: RegisterSharedGateCounts,
    pub reference_steps: usize,
    pub scheduled_window_sum: usize,
    pub scheduled_add_toffoli: usize,
    pub scheduled_sub_toffoli: usize,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ReferenceNormalizedControlProofReport {
    pub length_widths_checked: usize,
    pub basis_states_checked: usize,
    pub inverse_pair_checks: usize,
    pub scratch_clean_checks: usize,
    pub oracle_transition_checks: usize,
    pub phase_update9: RegisterSharedGateCounts,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ReferenceLengthSwapProofReport {
    pub work_widths_checked: usize,
    pub basis_states_checked: usize,
    pub quadratic_oracle_equivalence_checks: usize,
    pub inverse_pair_checks: usize,
    pub scratch_clean_checks: usize,
    pub control_off_identity_checks: usize,
    pub conditional_work_length_swap259: RegisterSharedGateCounts,
    pub conditional_steps: usize,
    pub standalone_active_qubits: usize,
    pub standalone_peak_qubits: usize,
    pub temporary_peak_qubits: usize,
    pub projected_reference_peak_qubits: usize,
    pub scheduled_toffoli: usize,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ReferenceWholeStepProofReport {
    pub modulus: u64,
    pub nonzero_inputs_checked: usize,
    pub steps_per_input: usize,
    pub boundary_transitions_checked: usize,
    pub inverse_transition_checks: usize,
    pub scratch_clean_checks: usize,
    pub data_qubits: usize,
    pub step_active_qubits: usize,
    pub step_peak_qubits: usize,
    pub temporary_peak_qubits: usize,
    pub emitted_ops: usize,
    pub emitted_toffoli: usize,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct ReferenceScheduledInversionProfile {
    pub steps: usize,
    pub inversion_state_qubits: usize,
    pub passenger_qubits: usize,
    pub point_add_state_qubits: usize,
    pub inversion_peak_qubits: usize,
    pub projected_point_add_peak_qubits: usize,
    pub emitted_ops: usize,
    pub emitted_toffoli: usize,
    pub emitted_hmr: usize,
    pub emitted_resets: usize,
}

fn measurement_classical_gate_counts(ops: &[crate::circuit::Op]) -> RegisterSharedGateCounts {
    use crate::circuit::OperationType;

    let mut counts = RegisterSharedGateCounts::default();
    for operation in ops {
        match operation.kind {
            OperationType::X => counts.x += 1,
            OperationType::CX => counts.cx += 1,
            OperationType::CCX | OperationType::CCZ => counts.ccx += 1,
            OperationType::Neg
            | OperationType::Z
            | OperationType::CZ
            | OperationType::R
            | OperationType::Hmr
            | OperationType::PushCondition
            | OperationType::PopCondition => {}
            other => panic!("length-swap proof emitted unsupported operation {other:?}"),
        }
    }
    counts.total = counts.x + counts.cx + counts.ccx;
    counts
}

fn ceil_safe(value: f64) -> isize {
    (value - 1e-12).ceil() as isize
}

fn floor_safe(value: f64) -> isize {
    (value + 1e-12).floor() as isize
}

#[must_use]
pub fn reference_active_windows(n: usize, step: usize) -> ReferenceActiveWindows {
    assert!(step > 0);
    let n = n as isize;
    let step = step as isize;
    let phi = (5.0f64.sqrt() + 1.0) / 2.0;
    let c = 1.0 / phi.log2();
    let k1 = ceil_safe((step as f64 - (n + 2) as f64) / (4.0 * c - 1.0)).max(1) + 2;
    let upper1 = n + 3;
    let k2 = ceil_safe((step as f64 - 3.0 * (n + 2) as f64) / (4.0 * c - 3.0)).max(1) + 1;
    let upper2 = floor_safe(step as f64 / 2.0).min(n) + 2;
    let upper3 = ceil_safe(step as f64 / 4.0).min(n) + 1;
    let k4 = ceil_safe((step as f64 - 4.0 * (n + 2) as f64) / (4.0 * c - 4.0)).max(1);
    let upper4 = floor_safe(step as f64 / 4.0 + 3.0).min(n + 3);
    let k5 = ceil_safe(step as f64 / (4.0 * c));
    let upper5 = floor_safe(step as f64 / 4.0 + 4.0).min(n + 3);
    assert!(k1 <= upper1 && k2 <= upper2 && 1 <= upper3 && k4 <= upper4 && k5 <= upper5);
    ReferenceActiveWindows {
        r_add_sub: (k1 as usize, upper1 as usize),
        quotient_swap: (k2 as usize, upper2 as usize),
        t_add_sub: (1, upper3 as usize),
        length_update_t: (k4 as usize, upper4 as usize),
        length_update_r: (k5 as usize, upper5 as usize),
    }
}

fn inclusive_width(window: (usize, usize)) -> usize {
    window.1 - window.0 + 1
}

#[must_use]
pub fn exhaustive_reference_schedule_check() -> ReferenceScheduleProofReport {
    let mut report = ReferenceScheduleProofReport {
        steps_checked: 0,
        r_window_sum: 0,
        quotient_swap_window_sum: 0,
        t_window_sum: 0,
        length_t_window_sum: 0,
        length_r_window_sum: 0,
        maximum_r_window: 0,
        maximum_quotient_swap_window: 0,
        maximum_t_window: 0,
        maximum_length_t_window: 0,
        maximum_length_r_window: 0,
    };
    for step in 1..=REFERENCE_STEPS {
        let windows = reference_active_windows(256, step);
        let r = inclusive_width(windows.r_add_sub);
        let swap = inclusive_width(windows.quotient_swap);
        let t = inclusive_width(windows.t_add_sub);
        let length_t = inclusive_width(windows.length_update_t);
        let length_r = inclusive_width(windows.length_update_r);
        report.steps_checked += 1;
        report.r_window_sum += r;
        report.quotient_swap_window_sum += swap;
        report.t_window_sum += t;
        report.length_t_window_sum += length_t;
        report.length_r_window_sum += length_r;
        report.maximum_r_window = report.maximum_r_window.max(r);
        report.maximum_quotient_swap_window = report.maximum_quotient_swap_window.max(swap);
        report.maximum_t_window = report.maximum_t_window.max(t);
        report.maximum_length_t_window = report.maximum_length_t_window.max(length_t);
        report.maximum_length_r_window = report.maximum_length_r_window.max(length_r);
    }
    report
}

fn majority(circ: &mut Circuit, a: &QReg, b: &QReg, carry: &QReg) {
    circ.cx(a, b);
    circ.cx(a, carry);
    circ.ccx(carry, b, a);
}

fn unmajority_add(circ: &mut Circuit, a: &QReg, b: &QReg, carry: &QReg) {
    circ.ccx(carry, b, a);
    circ.cx(a, carry);
    circ.cx(carry, b);
}

fn majority_inverse(circ: &mut Circuit, a: &QReg, b: &QReg, carry: &QReg) {
    circ.ccx(carry, b, a);
    circ.cx(a, carry);
    circ.cx(a, b);
}

fn unmajority_add_inverse(circ: &mut Circuit, a: &QReg, b: &QReg, carry: &QReg) {
    circ.cx(carry, b);
    circ.cx(a, carry);
    circ.ccx(carry, b, a);
}

/// `b += a mod 2^n`; `carry` is restored and `overflow` receives carry-out.
pub fn cuccaro_add_mod_2n(
    circ: &mut Circuit,
    a: &[QReg],
    b: &[QReg],
    carry: &QReg,
    overflow: &QReg,
) {
    assert_eq!(a.len(), b.len());
    assert!(!a.is_empty());
    majority(circ, &a[0], &b[0], carry);
    for index in 1..a.len() {
        majority(circ, &a[index], &b[index], &a[index - 1]);
    }
    circ.cx(&a[a.len() - 1], overflow);
    for index in (1..a.len()).rev() {
        unmajority_add(circ, &a[index], &b[index], &a[index - 1]);
    }
    unmajority_add(circ, &a[0], &b[0], carry);
}

/// Exact inverse of [`cuccaro_add_mod_2n`].
pub fn cuccaro_sub_mod_2n(
    circ: &mut Circuit,
    a: &[QReg],
    b: &[QReg],
    carry: &QReg,
    overflow: &QReg,
) {
    assert_eq!(a.len(), b.len());
    assert!(!a.is_empty());
    unmajority_add_inverse(circ, &a[0], &b[0], carry);
    for index in 1..a.len() {
        unmajority_add_inverse(circ, &a[index], &b[index], &a[index - 1]);
    }
    circ.cx(&a[a.len() - 1], overflow);
    for index in (1..a.len()).rev() {
        majority_inverse(circ, &a[index], &b[index], &a[index - 1]);
    }
    majority_inverse(circ, &a[0], &b[0], carry);
}

pub fn toggle_constant(circ: &mut Circuit, register: &[QReg], value: usize) {
    let reduced = value % (1usize << register.len());
    for (index, bit) in register.iter().enumerate() {
        if ((reduced >> index) & 1) != 0 {
            circ.x(bit);
        }
    }
}

pub fn add_const_mod_2n(circ: &mut Circuit, register: &[QReg], value: usize, scratch: &[QReg]) {
    assert!(scratch.len() >= register.len() + 2);
    let (constant, tail) = scratch.split_at(register.len());
    toggle_constant(circ, constant, value);
    cuccaro_add_mod_2n(circ, constant, register, &tail[0], &tail[1]);
    toggle_constant(circ, constant, value);
}

pub fn sub_const_mod_2n(circ: &mut Circuit, register: &[QReg], value: usize, scratch: &[QReg]) {
    assert!(scratch.len() >= register.len() + 2);
    let (constant, tail) = scratch.split_at(register.len());
    toggle_constant(circ, constant, value);
    cuccaro_sub_mod_2n(circ, constant, register, &tail[0], &tail[1]);
    toggle_constant(circ, constant, value);
}

fn equality_flag(
    circ: &mut Circuit,
    control: &QReg,
    register: &[QReg],
    value: usize,
    flag: &QReg,
    chain: &[QReg],
) {
    let reduced = value % (1usize << register.len());
    for (index, bit) in register.iter().enumerate() {
        if ((reduced >> index) & 1) == 0 {
            circ.x(bit);
        }
    }
    let mut controls = Vec::with_capacity(register.len() + 1);
    controls.push(control);
    controls.extend(register.iter());
    multi_controlled_x_vchain(circ, &controls, flag, chain);
    for (index, bit) in register.iter().enumerate() {
        if ((reduced >> index) & 1) == 0 {
            circ.x(bit);
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn controlled_dynamic_swaps(
    circ: &mut Circuit,
    active: &QReg,
    sign: &QReg,
    work1_window: &[QReg],
    sum: &[QReg],
    first_global_index: usize,
    flag: &QReg,
    chain: &[QReg],
    reverse: bool,
) {
    if reverse {
        for (local_index, work_bit) in work1_window.iter().enumerate().rev() {
            let global_index = first_global_index + local_index;
            equality_flag(circ, active, sum, global_index, flag, chain);
            circ.cswap(flag, work_bit, sign);
            equality_flag(circ, active, sum, global_index, flag, chain);
        }
    } else {
        for (local_index, work_bit) in work1_window.iter().enumerate() {
            let global_index = first_global_index + local_index;
            equality_flag(circ, active, sum, global_index, flag, chain);
            circ.cswap(flag, work_bit, sign);
            equality_flag(circ, active, sum, global_index, flag, chain);
        }
    }
}

/// Exact dynamic swap used by the quotient phase.
///
/// The active predicate is `phase1 XOR phase2`. The selected packed Work1 bit
/// is indexed by `l_t + l_q + 1` while the quotient grows and by `l_t + l_q`
/// while it shrinks. This accounts for the zero delimiter between little-endian
/// `t` and big-endian `q`. The construction is linear in the active window and
/// uses exactly `length_width + 2` clean scratch lanes.
#[allow(clippy::too_many_arguments)]
pub fn location_controlled_swap_one_hot(
    circ: &mut Circuit,
    phase1: &QReg,
    phase2: &QReg,
    sign: &QReg,
    work1_window: &[QReg],
    first_global_index: usize,
    l_t: &[QReg],
    l_q: &[QReg],
    scratch: &[QReg],
) {
    assert_eq!(l_t.len(), l_q.len());
    assert!(!l_t.is_empty());
    assert!(scratch.len() >= l_t.len() + 2);
    let carry = &scratch[0];
    let overflow = &scratch[1];
    let flag = &scratch[2];
    let chain = &scratch[3..];
    assert!(chain.len() >= l_t.len().saturating_sub(1));

    circ.cx(phase1, phase2);
    let controls = [phase2, phase1];
    circ.x(phase1);
    controlled_add_one(circ, &controls, l_q, chain);
    circ.x(phase1);
    cuccaro_add_mod_2n(circ, l_t, l_q, carry, overflow);
    controlled_dynamic_swaps(
        circ,
        phase2,
        sign,
        work1_window,
        l_q,
        first_global_index,
        flag,
        chain,
        false,
    );
    cuccaro_sub_mod_2n(circ, l_t, l_q, carry, overflow);
    controlled_sub_one(circ, &controls, l_q, chain);
    circ.cx(phase1, phase2);
}

/// Exact inverse of [`location_controlled_swap_one_hot`].
#[allow(clippy::too_many_arguments)]
pub fn location_controlled_swap_one_hot_inverse(
    circ: &mut Circuit,
    phase1: &QReg,
    phase2: &QReg,
    sign: &QReg,
    work1_window: &[QReg],
    first_global_index: usize,
    l_t: &[QReg],
    l_q: &[QReg],
    scratch: &[QReg],
) {
    assert_eq!(l_t.len(), l_q.len());
    assert!(!l_t.is_empty());
    assert!(scratch.len() >= l_t.len() + 2);
    let carry = &scratch[0];
    let overflow = &scratch[1];
    let flag = &scratch[2];
    let chain = &scratch[3..];
    assert!(chain.len() >= l_t.len().saturating_sub(1));

    circ.cx(phase1, phase2);
    let controls = [phase2, phase1];
    controlled_add_one(circ, &controls, l_q, chain);
    cuccaro_add_mod_2n(circ, l_t, l_q, carry, overflow);
    controlled_dynamic_swaps(
        circ,
        phase2,
        sign,
        work1_window,
        l_q,
        first_global_index,
        flag,
        chain,
        true,
    );
    cuccaro_sub_mod_2n(circ, l_t, l_q, carry, overflow);
    circ.x(phase1);
    controlled_sub_one(circ, &controls, l_q, chain);
    circ.x(phase1);
    circ.cx(phase1, phase2);
}

fn toggle_control_and_nonnegative(
    circ: &mut Circuit,
    control: &QReg,
    signed_length: &[QReg],
    active: &QReg,
) {
    let sign = signed_length.last().expect("nonempty signed length");
    circ.x(sign);
    circ.ccx(control, sign, active);
    circ.x(sign);
}

fn compute_t_sub_enable(
    circ: &mut Circuit,
    phase1: &QReg,
    phase2: &QReg,
    sign: &QReg,
    condition: &QReg,
    enable: &QReg,
) {
    // condition = phase2 OR !sign; enable = phase1 AND condition.
    circ.cx(phase2, condition);
    circ.x(sign);
    circ.cx(sign, condition);
    circ.ccx(phase2, sign, condition);
    circ.x(sign);
    circ.ccx(phase1, condition, enable);
}

fn uncompute_t_sub_enable(
    circ: &mut Circuit,
    phase1: &QReg,
    phase2: &QReg,
    sign: &QReg,
    condition: &QReg,
    enable: &QReg,
) {
    circ.ccx(phase1, condition, enable);
    circ.x(sign);
    circ.ccx(phase2, sign, condition);
    circ.cx(sign, condition);
    circ.x(sign);
    circ.cx(phase2, condition);
}

/// Corrected little-endian coefficient-add variant derived from
/// `location_controlled_add_gate_single` at pinned commit b836f59.
///
/// The public reference iterates the packed coefficient lanes in the opposite
/// direction and omits the standalone sign update used by this route.
#[allow(clippy::too_many_arguments)]
pub fn coefficient_add_single(
    circ: &mut Circuit,
    phase1: &QReg,
    sign: &QReg,
    work1: &[QReg],
    work2: &[QReg],
    l_t: &[QReg],
    scratch: &[QReg],
) {
    assert_eq!(work1.len(), work2.len());
    assert!(!work1.is_empty() && !l_t.is_empty());
    assert!(scratch.len() >= l_t.len() + 2);
    let carry = &scratch[0];
    let active = &scratch[1];
    let tmp = &scratch[2];
    let length_chain = &scratch[3..];

    for index in 0..work1.len() {
        toggle_control_and_nonnegative(circ, phase1, l_t, active);
        circ.ccx(active, carry, &work2[index]);
        circ.ccx(active, carry, &work1[index]);
        multi_controlled_x_vchain(
            circ,
            &[active, &work1[index], &work2[index]],
            carry,
            std::slice::from_ref(tmp),
        );
        toggle_control_and_nonnegative(circ, phase1, l_t, active);
        if index + 1 != work1.len() {
            decrement_mod_2n(circ, l_t, length_chain);
        }
    }

    circ.ccx(phase1, carry, sign);
    circ.cx(phase1, sign);

    for index in (0..work1.len()).rev() {
        toggle_control_and_nonnegative(circ, phase1, l_t, active);
        multi_controlled_x_vchain(
            circ,
            &[active, &work1[index], &work2[index]],
            carry,
            std::slice::from_ref(tmp),
        );
        circ.ccx(active, carry, &work1[index]);
        circ.ccx(active, &work1[index], &work2[index]);
        toggle_control_and_nonnegative(circ, phase1, l_t, active);
        if index != 0 {
            increment_mod_2n(circ, l_t, length_chain);
        }
    }
}

/// Exact inverse of [`coefficient_add_single`].
#[allow(clippy::too_many_arguments)]
pub fn coefficient_add_single_inverse(
    circ: &mut Circuit,
    phase1: &QReg,
    sign: &QReg,
    work1: &[QReg],
    work2: &[QReg],
    l_t: &[QReg],
    scratch: &[QReg],
) {
    assert_eq!(work1.len(), work2.len());
    assert!(!work1.is_empty() && !l_t.is_empty());
    assert!(scratch.len() >= l_t.len() + 2);
    let carry = &scratch[0];
    let active = &scratch[1];
    let tmp = &scratch[2];
    let length_chain = &scratch[3..];

    for index in 0..work1.len() {
        if index != 0 {
            decrement_mod_2n(circ, l_t, length_chain);
        }
        toggle_control_and_nonnegative(circ, phase1, l_t, active);
        circ.ccx(active, &work1[index], &work2[index]);
        circ.ccx(active, carry, &work1[index]);
        multi_controlled_x_vchain(
            circ,
            &[active, &work1[index], &work2[index]],
            carry,
            std::slice::from_ref(tmp),
        );
        toggle_control_and_nonnegative(circ, phase1, l_t, active);
    }

    circ.cx(phase1, sign);
    circ.ccx(phase1, carry, sign);

    for index in (0..work1.len()).rev() {
        if index + 1 != work1.len() {
            increment_mod_2n(circ, l_t, length_chain);
        }
        toggle_control_and_nonnegative(circ, phase1, l_t, active);
        multi_controlled_x_vchain(
            circ,
            &[active, &work1[index], &work2[index]],
            carry,
            std::slice::from_ref(tmp),
        );
        circ.ccx(active, carry, &work1[index]);
        circ.ccx(active, carry, &work2[index]);
        toggle_control_and_nonnegative(circ, phase1, l_t, active);
    }
}

/// Corrected little-endian coefficient-subtract variant derived from
/// `location_controlled_sub_gate_single` at pinned commit b836f59.
#[allow(clippy::too_many_arguments)]
pub fn coefficient_sub_single(
    circ: &mut Circuit,
    phase1: &QReg,
    phase2: &QReg,
    sign: &QReg,
    work1: &[QReg],
    work2: &[QReg],
    l_t: &[QReg],
    scratch: &[QReg],
) {
    assert_eq!(work1.len(), work2.len());
    assert!(!work1.is_empty() && !l_t.is_empty());
    assert!(scratch.len() >= l_t.len() + 4);
    let carry = &scratch[0];
    let active = &scratch[1];
    let tmp = &scratch[2];
    let condition = &scratch[3];
    let enable = &scratch[4];
    let length_chain = &scratch[5..];

    compute_t_sub_enable(circ, phase1, phase2, sign, condition, enable);
    for index in 0..work1.len() {
        toggle_control_and_nonnegative(circ, enable, l_t, active);
        circ.ccx(active, &work1[index], &work2[index]);
        circ.ccx(active, carry, &work1[index]);
        multi_controlled_x_vchain(
            circ,
            &[active, &work1[index], &work2[index]],
            carry,
            std::slice::from_ref(tmp),
        );
        toggle_control_and_nonnegative(circ, enable, l_t, active);
        if index + 1 != work1.len() {
            decrement_mod_2n(circ, l_t, length_chain);
        }
    }
    for index in (0..work1.len()).rev() {
        toggle_control_and_nonnegative(circ, enable, l_t, active);
        multi_controlled_x_vchain(
            circ,
            &[active, &work1[index], &work2[index]],
            carry,
            std::slice::from_ref(tmp),
        );
        circ.ccx(active, carry, &work1[index]);
        circ.ccx(active, carry, &work2[index]);
        toggle_control_and_nonnegative(circ, enable, l_t, active);
        if index != 0 {
            increment_mod_2n(circ, l_t, length_chain);
        }
    }
    uncompute_t_sub_enable(circ, phase1, phase2, sign, condition, enable);
}

/// Exact inverse of [`coefficient_sub_single`].
#[allow(clippy::too_many_arguments)]
pub fn coefficient_sub_single_inverse(
    circ: &mut Circuit,
    phase1: &QReg,
    phase2: &QReg,
    sign: &QReg,
    work1: &[QReg],
    work2: &[QReg],
    l_t: &[QReg],
    scratch: &[QReg],
) {
    assert_eq!(work1.len(), work2.len());
    assert!(!work1.is_empty() && !l_t.is_empty());
    assert!(scratch.len() >= l_t.len() + 4);
    let carry = &scratch[0];
    let active = &scratch[1];
    let tmp = &scratch[2];
    let condition = &scratch[3];
    let enable = &scratch[4];
    let length_chain = &scratch[5..];

    compute_t_sub_enable(circ, phase1, phase2, sign, condition, enable);
    for index in 0..work1.len() {
        if index != 0 {
            decrement_mod_2n(circ, l_t, length_chain);
        }
        toggle_control_and_nonnegative(circ, enable, l_t, active);
        circ.ccx(active, carry, &work2[index]);
        circ.ccx(active, carry, &work1[index]);
        multi_controlled_x_vchain(
            circ,
            &[active, &work1[index], &work2[index]],
            carry,
            std::slice::from_ref(tmp),
        );
        toggle_control_and_nonnegative(circ, enable, l_t, active);
    }
    for index in (0..work1.len()).rev() {
        if index + 1 != work1.len() {
            increment_mod_2n(circ, l_t, length_chain);
        }
        toggle_control_and_nonnegative(circ, enable, l_t, active);
        multi_controlled_x_vchain(
            circ,
            &[active, &work1[index], &work2[index]],
            carry,
            std::slice::from_ref(tmp),
        );
        circ.ccx(active, carry, &work1[index]);
        circ.ccx(active, &work1[index], &work2[index]);
        toggle_control_and_nonnegative(circ, enable, l_t, active);
    }
    uncompute_t_sub_enable(circ, phase1, phase2, sign, condition, enable);
}

struct RemainderScratch<'a> {
    carry: &'a QReg,
    active: &'a QReg,
    tmp: &'a QReg,
    length_carry: &'a QReg,
    length_overflow: &'a QReg,
    constant: &'a [QReg],
    walk: &'a [QReg],
    nonzero: &'a QReg,
    operation: &'a QReg,
    phase_sign: &'a QReg,
    enable: &'a QReg,
    nonzero_chain: &'a [QReg],
}

fn remainder_scratch_width(length_width: usize, remainder_length_width: usize) -> usize {
    5 + (length_width + 2)
        + length_width.saturating_sub(1)
        + 4
        + remainder_length_width.saturating_sub(2)
}

fn split_remainder_scratch<'a>(
    scratch: &'a [QReg],
    length_width: usize,
    remainder_length_width: usize,
) -> RemainderScratch<'a> {
    assert!(scratch.len() >= remainder_scratch_width(length_width, remainder_length_width));
    let mut offset = 0usize;
    let carry = &scratch[offset];
    offset += 1;
    let active = &scratch[offset];
    offset += 1;
    let tmp = &scratch[offset];
    offset += 1;
    let length_carry = &scratch[offset];
    offset += 1;
    let length_overflow = &scratch[offset];
    offset += 1;
    let constant = &scratch[offset..offset + length_width + 2];
    offset += length_width + 2;
    let walk = &scratch[offset..offset + length_width.saturating_sub(1)];
    offset += length_width.saturating_sub(1);
    let nonzero = &scratch[offset];
    offset += 1;
    let operation = &scratch[offset];
    offset += 1;
    let phase_sign = &scratch[offset];
    offset += 1;
    let enable = &scratch[offset];
    offset += 1;
    let nonzero_chain = &scratch[offset..offset + remainder_length_width.saturating_sub(2)];
    RemainderScratch {
        carry,
        active,
        tmp,
        length_carry,
        length_overflow,
        constant,
        walk,
        nonzero,
        operation,
        phase_sign,
        enable,
        nonzero_chain,
    }
}

fn compute_nonzero(circ: &mut Circuit, register: &[QReg], nonzero: &QReg, chain: &[QReg]) {
    assert!(!register.is_empty());
    for bit in register {
        circ.x(bit);
    }
    let controls: Vec<&QReg> = register.iter().collect();
    multi_controlled_x_vchain(circ, &controls, nonzero, chain);
    for bit in register {
        circ.x(bit);
    }
    circ.x(nonzero);
}

fn uncompute_nonzero(circ: &mut Circuit, register: &[QReg], nonzero: &QReg, chain: &[QReg]) {
    circ.x(nonzero);
    for bit in register {
        circ.x(bit);
    }
    let controls: Vec<&QReg> = register.iter().collect();
    multi_controlled_x_vchain(circ, &controls, nonzero, chain);
    for bit in register {
        circ.x(bit);
    }
}

fn compute_zero(circ: &mut Circuit, register: &[QReg], zero: &QReg, chain: &[QReg]) {
    assert!(!register.is_empty());
    for bit in register {
        circ.x(bit);
    }
    let controls: Vec<&QReg> = register.iter().collect();
    multi_controlled_x_vchain(circ, &controls, zero, chain);
    for bit in register {
        circ.x(bit);
    }
}

fn uncompute_zero(circ: &mut Circuit, register: &[QReg], zero: &QReg, chain: &[QReg]) {
    compute_zero(circ, register, zero, chain);
}

fn normalized_phase_scratch_width(length_width: usize) -> usize {
    5 + length_width.saturating_sub(2)
}

#[allow(clippy::too_many_arguments)]
pub fn normalized_phase_update(
    circ: &mut Circuit,
    phase1: &QReg,
    phase2: &QReg,
    sign: &QReg,
    l_q: &[QReg],
    l_r_prime: &[QReg],
    l_s: &[QReg],
    scratch: &[QReg],
) {
    assert_eq!(l_q.len(), l_r_prime.len());
    assert_eq!(l_q.len(), l_s.len());
    assert!(scratch.len() >= normalized_phase_scratch_width(l_q.len()));
    let zero_q = &scratch[0];
    let nonzero_r = &scratch[1];
    let zero_s = &scratch[2];
    let condition = &scratch[3];
    let temporary = &scratch[4];
    let chain = &scratch[5..];

    compute_zero(circ, l_q, zero_q, chain);
    compute_nonzero(circ, l_r_prime, nonzero_r, chain);
    compute_zero(circ, l_s, zero_s, chain);
    circ.ccx(zero_q, nonzero_r, condition);

    circ.cx(sign, temporary);
    circ.cx(phase1, temporary);
    circ.ccx(condition, temporary, phase2);
    circ.cx(phase1, temporary);
    circ.cx(sign, temporary);
    circ.ccx(condition, phase2, sign);
    circ.cx(zero_s, phase1);
    circ.cx(zero_s, phase2);

    circ.ccx(zero_q, nonzero_r, condition);
    uncompute_zero(circ, l_s, zero_s, chain);
    uncompute_nonzero(circ, l_r_prime, nonzero_r, chain);
    uncompute_zero(circ, l_q, zero_q, chain);
}

#[allow(clippy::too_many_arguments)]
pub fn normalized_phase_update_inverse(
    circ: &mut Circuit,
    phase1: &QReg,
    phase2: &QReg,
    sign: &QReg,
    l_q: &[QReg],
    l_r_prime: &[QReg],
    l_s: &[QReg],
    scratch: &[QReg],
) {
    assert_eq!(l_q.len(), l_r_prime.len());
    assert_eq!(l_q.len(), l_s.len());
    assert!(scratch.len() >= normalized_phase_scratch_width(l_q.len()));
    let zero_q = &scratch[0];
    let nonzero_r = &scratch[1];
    let zero_s = &scratch[2];
    let condition = &scratch[3];
    let temporary = &scratch[4];
    let chain = &scratch[5..];

    compute_zero(circ, l_q, zero_q, chain);
    compute_nonzero(circ, l_r_prime, nonzero_r, chain);
    compute_zero(circ, l_s, zero_s, chain);
    circ.ccx(zero_q, nonzero_r, condition);

    circ.cx(zero_s, phase2);
    circ.cx(zero_s, phase1);
    circ.ccx(condition, phase2, sign);
    circ.cx(sign, temporary);
    circ.cx(phase1, temporary);
    circ.ccx(condition, temporary, phase2);
    circ.cx(phase1, temporary);
    circ.cx(sign, temporary);

    circ.ccx(zero_q, nonzero_r, condition);
    uncompute_zero(circ, l_s, zero_s, chain);
    uncompute_nonzero(circ, l_r_prime, nonzero_r, chain);
    uncompute_zero(circ, l_q, zero_q, chain);
}

fn toggle_or_latch(circ: &mut Circuit, latch: &QReg, condition: &QReg, temporary: &QReg) {
    circ.x(latch);
    circ.ccx(condition, latch, temporary);
    circ.x(latch);
    circ.cx(temporary, latch);
    circ.x(latch);
    circ.ccx(condition, latch, temporary);
    circ.x(latch);
}

fn controlled_add_no_overflow(
    circ: &mut Circuit,
    control: &QReg,
    source: &[QReg],
    target: &[QReg],
    carry: &QReg,
    temporary: &QReg,
) {
    assert_eq!(source.len(), target.len());
    assert!(!source.is_empty());
    let controlled_cx = |circ: &mut Circuit, source: &QReg, target: &QReg| {
        circ.ccx(control, source, target);
    };
    let controlled_ccx = |circ: &mut Circuit, left: &QReg, right: &QReg, target: &QReg| {
        multi_controlled_x_vchain(
            circ,
            &[control, left, right],
            target,
            std::slice::from_ref(temporary),
        );
    };

    controlled_cx(circ, carry, &target[0]);
    controlled_cx(circ, carry, &source[0]);
    controlled_ccx(circ, &source[0], &target[0], carry);
    for index in 1..source.len() {
        controlled_cx(circ, &source[index - 1], &target[index]);
        controlled_cx(circ, &source[index - 1], &source[index]);
        controlled_ccx(circ, &source[index], &target[index], &source[index - 1]);
    }
    for index in (1..source.len()).rev() {
        controlled_ccx(circ, &source[index], &target[index], &source[index - 1]);
        controlled_cx(circ, &source[index - 1], &source[index]);
        controlled_cx(circ, &source[index], &target[index]);
    }
    controlled_ccx(circ, &source[0], &target[0], carry);
    controlled_cx(circ, carry, &source[0]);
    controlled_cx(circ, &source[0], &target[0]);
}

fn controlled_add_no_overflow_inverse(
    circ: &mut Circuit,
    control: &QReg,
    source: &[QReg],
    target: &[QReg],
    carry: &QReg,
    temporary: &QReg,
) {
    assert_eq!(source.len(), target.len());
    assert!(!source.is_empty());
    let controlled_cx = |circ: &mut Circuit, source: &QReg, target: &QReg| {
        circ.ccx(control, source, target);
    };
    let controlled_ccx = |circ: &mut Circuit, left: &QReg, right: &QReg, target: &QReg| {
        multi_controlled_x_vchain(
            circ,
            &[control, left, right],
            target,
            std::slice::from_ref(temporary),
        );
    };

    controlled_cx(circ, &source[0], &target[0]);
    controlled_cx(circ, carry, &source[0]);
    controlled_ccx(circ, &source[0], &target[0], carry);
    for index in 1..source.len() {
        controlled_cx(circ, &source[index], &target[index]);
        controlled_cx(circ, &source[index - 1], &source[index]);
        controlled_ccx(circ, &source[index], &target[index], &source[index - 1]);
    }
    for index in (1..source.len()).rev() {
        controlled_ccx(circ, &source[index], &target[index], &source[index - 1]);
        controlled_cx(circ, &source[index - 1], &source[index]);
        controlled_cx(circ, &source[index - 1], &target[index]);
    }
    controlled_ccx(circ, &source[0], &target[0], carry);
    controlled_cx(circ, carry, &source[0]);
    controlled_cx(circ, carry, &target[0]);
}

fn length_scan_scratch_width(length_width: usize) -> usize {
    4 + (length_width + 2) + length_width.saturating_sub(1)
}

#[allow(clippy::too_many_arguments)]
fn length_scan_compute(
    circ: &mut Circuit,
    total_work_width: usize,
    window_lower: usize,
    window_upper: usize,
    work1: &[QReg],
    work2: &[QReg],
    l_s: &[QReg],
    l_q: &[QReg],
    l_r_prime: &[QReg],
    scratch: &[QReg],
) {
    assert_eq!(work1.len(), total_work_width);
    assert_eq!(work2.len(), total_work_width);
    assert!(1 <= window_lower && window_lower <= window_upper);
    assert!(window_upper <= total_work_width);
    assert_eq!(l_s.len(), l_q.len());
    assert_eq!(l_s.len(), l_r_prime.len());
    assert!(scratch.len() >= length_scan_scratch_width(l_s.len()));
    let latch_u = &scratch[0];
    let latch_v = &scratch[1];
    let condition = &scratch[2];
    let temporary = &scratch[3];
    let constant = &scratch[4..4 + l_s.len() + 2];
    let pool = &scratch[4 + l_s.len() + 2..];
    let sign_s = l_s.last().expect("nonempty l_s");
    let sign_q = l_q.last().expect("nonempty l_q");
    let sign_r = l_r_prime.last().expect("nonempty l_r_prime");

    sub_const_mod_2n(
        circ,
        l_r_prime,
        total_work_width + 1 - window_upper,
        constant,
    );
    for position in (window_lower..=window_upper).rev() {
        let index = position - 1;
        multi_controlled_x_vchain(circ, &[&work1[index], sign_s, sign_r], condition, pool);
        toggle_or_latch(circ, latch_u, condition, temporary);
        multi_controlled_x_vchain(circ, &[&work1[index], sign_s, sign_r], condition, pool);
        controlled_increment_mod_2n(circ, latch_u, l_s, pool);

        multi_controlled_x_vchain(circ, &[&work2[index], sign_q, sign_r], condition, pool);
        toggle_or_latch(circ, latch_v, condition, temporary);
        multi_controlled_x_vchain(circ, &[&work2[index], sign_q, sign_r], condition, pool);
        controlled_increment_mod_2n(circ, latch_v, l_q, pool);
        controlled_decrement_mod_2n(circ, sign_r, l_r_prime, pool);
    }
    let carry = &constant[l_s.len()];
    let overflow = &constant[l_s.len() + 1];
    cuccaro_sub_mod_2n(circ, l_q, l_s, carry, overflow);
}

#[allow(clippy::too_many_arguments)]
fn length_scan_uncompute(
    circ: &mut Circuit,
    total_work_width: usize,
    window_lower: usize,
    window_upper: usize,
    work1: &[QReg],
    work2: &[QReg],
    l_s: &[QReg],
    l_q: &[QReg],
    l_r_prime: &[QReg],
    scratch: &[QReg],
) {
    assert!(scratch.len() >= length_scan_scratch_width(l_s.len()));
    let latch_u = &scratch[0];
    let latch_v = &scratch[1];
    let condition = &scratch[2];
    let temporary = &scratch[3];
    let constant = &scratch[4..4 + l_s.len() + 2];
    let pool = &scratch[4 + l_s.len() + 2..];
    let sign_s = l_s.last().expect("nonempty l_s");
    let sign_q = l_q.last().expect("nonempty l_q");
    let sign_r = l_r_prime.last().expect("nonempty l_r_prime");
    let carry = &constant[l_s.len()];
    let overflow = &constant[l_s.len() + 1];

    cuccaro_add_mod_2n(circ, l_q, l_s, carry, overflow);
    for position in window_lower..=window_upper {
        let index = position - 1;
        controlled_increment_mod_2n(circ, sign_r, l_r_prime, pool);

        controlled_decrement_mod_2n(circ, latch_v, l_q, pool);
        multi_controlled_x_vchain(circ, &[&work2[index], sign_q, sign_r], condition, pool);
        toggle_or_latch(circ, latch_v, condition, temporary);
        multi_controlled_x_vchain(circ, &[&work2[index], sign_q, sign_r], condition, pool);

        controlled_decrement_mod_2n(circ, latch_u, l_s, pool);
        multi_controlled_x_vchain(circ, &[&work1[index], sign_s, sign_r], condition, pool);
        toggle_or_latch(circ, latch_u, condition, temporary);
        multi_controlled_x_vchain(circ, &[&work1[index], sign_s, sign_r], condition, pool);
    }
    add_const_mod_2n(
        circ,
        l_r_prime,
        total_work_width + 1 - window_upper,
        constant,
    );
}

fn length_update_scratch_width(length_width: usize) -> usize {
    length_scan_scratch_width(length_width) + 2
}

#[allow(clippy::too_many_arguments)]
pub fn conditional_length_update(
    circ: &mut Circuit,
    total_work_width: usize,
    window_lower: usize,
    window_upper: usize,
    control: &QReg,
    work1: &[QReg],
    work2: &[QReg],
    l_s: &[QReg],
    l_q: &[QReg],
    l_r_prime: &[QReg],
    target: &[QReg],
    scratch: &[QReg],
) {
    let scan_width = length_scan_scratch_width(l_s.len());
    assert!(scratch.len() >= scan_width + 2);
    let scan = &scratch[..scan_width];
    let carry = &scratch[scan_width];
    let temporary = &scratch[scan_width + 1];
    length_scan_compute(
        circ,
        total_work_width,
        window_lower,
        window_upper,
        work1,
        work2,
        l_s,
        l_q,
        l_r_prime,
        scan,
    );
    controlled_add_no_overflow(circ, control, l_s, target, carry, temporary);
    length_scan_uncompute(
        circ,
        total_work_width,
        window_lower,
        window_upper,
        work1,
        work2,
        l_s,
        l_q,
        l_r_prime,
        scan,
    );
}

#[allow(clippy::too_many_arguments)]
pub fn conditional_length_update_inverse(
    circ: &mut Circuit,
    total_work_width: usize,
    window_lower: usize,
    window_upper: usize,
    control: &QReg,
    work1: &[QReg],
    work2: &[QReg],
    l_s: &[QReg],
    l_q: &[QReg],
    l_r_prime: &[QReg],
    target: &[QReg],
    scratch: &[QReg],
) {
    let scan_width = length_scan_scratch_width(l_s.len());
    assert!(scratch.len() >= scan_width + 2);
    let scan = &scratch[..scan_width];
    let carry = &scratch[scan_width];
    let temporary = &scratch[scan_width + 1];
    length_scan_compute(
        circ,
        total_work_width,
        window_lower,
        window_upper,
        work1,
        work2,
        l_s,
        l_q,
        l_r_prime,
        scan,
    );
    controlled_add_no_overflow_inverse(circ, control, l_s, target, carry, temporary);
    length_scan_uncompute(
        circ,
        total_work_width,
        window_lower,
        window_upper,
        work1,
        work2,
        l_s,
        l_q,
        l_r_prime,
        scan,
    );
}

#[allow(clippy::too_many_arguments)]
fn toggle_static_zero_flag(circ: &mut Circuit, source: &[&QReg], flag: &QReg) {
    use crate::point_add::trailmix_port::arith::khattar_gidney::{
        kg_prefix_ancilla_count, KgPrefixAnd,
    };

    if source.is_empty() {
        circ.x(flag);
        return;
    }
    for bit in source {
        circ.x(bit);
    }
    let ancillae = circ.alloc_qreg_bits(
        "rs.dynamic-bitlen.zero-prefix",
        kg_prefix_ancilla_count(source.len()),
    );
    let ancilla_refs: Vec<&QReg> = ancillae.iter().collect();
    let width = source.len();
    let done = KgPrefixAnd::new(source, &ancilla_refs).forward(circ, |circ, index, controls| {
        if index != width {
            return;
        }
        match controls {
            [] => circ.x(flag),
            [control] => circ.cx(control, flag),
            [left, right] => circ.ccx(left, right, flag),
            _ => unreachable!("KG prefix controls have width at most two"),
        }
    });
    done.reverse(circ, |_, _, _| {});
    for lane in ancillae {
        circ.zero_and_free(lane);
    }
    for bit in source {
        circ.x(bit);
    }
}

fn bit_length_lean_allow_zero_legacy(
    circ: &mut Circuit,
    source: &[&QReg],
    output: &[QReg],
    decrement: bool,
) {
    use super::shrunken_pz_state_machine::bit_length_lean;

    if source.is_empty() {
        return;
    }
    let zero = circ.alloc_qreg("rs.dynamic-bitlen.zero");
    let carries = circ.alloc_qreg_bits(
        "rs.dynamic-bitlen.zero-carries",
        output.len().saturating_sub(1),
    );
    toggle_static_zero_flag(circ, source, &zero);
    if decrement {
        controlled_increment_mod_2n(circ, &zero, output, &carries);
        bit_length_lean(circ, source, output, true);
    } else {
        bit_length_lean(circ, source, output, false);
        controlled_decrement_mod_2n(circ, &zero, output, &carries);
    }
    toggle_static_zero_flag(circ, source, &zero);
    for lane in carries {
        circ.zero_and_free(lane);
    }
    circ.zero_and_free(zero);
}

fn bit_length_lean_allow_zero(
    circ: &mut Circuit,
    source: &[&QReg],
    output: &[QReg],
    decrement: bool,
) {
    use super::shrunken_pz_state_machine::{
        bit_length_lean, lowq_fused_zero_prefix_bitlen_requested,
    };

    if !lowq_fused_zero_prefix_bitlen_requested() {
        bit_length_lean_allow_zero_legacy(circ, source, output, decrement);
        return;
    }
    if source.is_empty() {
        return;
    }
    bit_length_lean(circ, source, output, decrement);
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct FusedZeroPrefixLocalResources {
    pub extra_qubits: usize,
    pub emitted_ops: usize,
    pub emitted_toffoli: usize,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct FusedZeroPrefixBitLengthProofReport {
    pub widths_checked: usize,
    pub accumulator_width: usize,
    pub update_cases_checked: usize,
    pub controlled_cases_checked: usize,
    pub default_stream_equivalence_checks: usize,
    pub baseline_equivalence_checks: usize,
    pub phase_clean_checks: usize,
    pub ancilla_clean_checks: usize,
    pub maximum_baseline: FusedZeroPrefixLocalResources,
    pub maximum_fused: FusedZeroPrefixLocalResources,
    pub local_width: usize,
    pub local_accumulator_width: usize,
    pub local_add_baseline: FusedZeroPrefixLocalResources,
    pub local_add_fused: FusedZeroPrefixLocalResources,
    pub local_sub_baseline: FusedZeroPrefixLocalResources,
    pub local_sub_fused: FusedZeroPrefixLocalResources,
    pub scheduled_steps: usize,
    pub scheduled_baseline_inversion_peak_qubits: usize,
    pub scheduled_fused_inversion_peak_qubits: usize,
    pub scheduled_baseline_emitted_ops: usize,
    pub scheduled_fused_emitted_ops: usize,
    pub scheduled_baseline_emitted_toffoli: usize,
    pub scheduled_fused_emitted_toffoli: usize,
    pub scheduled_baseline_emitted_hmr: usize,
    pub scheduled_fused_emitted_hmr: usize,
    pub scheduled_baseline_emitted_resets: usize,
    pub scheduled_fused_emitted_resets: usize,
}

struct FusedZeroPrefixProofCircuit {
    builder: B,
    source_ids: Vec<u32>,
    accumulator_ids: Vec<u32>,
    control_id: Option<u32>,
}

fn configure_fused_zero_prefix_proof(fused: bool) {
    if fused {
        std::env::set_var("LOWQ_FUSED_ZERO_PREFIX_BITLEN", "1");
    } else {
        std::env::remove_var("LOWQ_FUSED_ZERO_PREFIX_BITLEN");
    }
}

fn build_fused_zero_prefix_update(
    source_width: usize,
    accumulator_width: usize,
    decrement: bool,
    legacy_entry: bool,
) -> FusedZeroPrefixProofCircuit {
    let mut circuit = Circuit::new();
    let source = circuit.alloc_qreg_bits("fused-zero-proof.source", source_width);
    let accumulator =
        circuit.alloc_qreg_bits("fused-zero-proof.accumulator", accumulator_width);
    let source_refs: Vec<&QReg> = source.iter().collect();
    if legacy_entry {
        bit_length_lean_allow_zero_legacy(
            &mut circuit,
            &source_refs,
            &accumulator,
            decrement,
        );
    } else {
        bit_length_lean_allow_zero(&mut circuit, &source_refs, &accumulator, decrement);
    }
    FusedZeroPrefixProofCircuit {
        source_ids: source.iter().map(QReg::id).collect(),
        accumulator_ids: accumulator.iter().map(QReg::id).collect(),
        control_id: None,
        builder: circuit.into_builder(),
    }
}

fn build_fused_zero_prefix_controlled_xor(
    source_width: usize,
    accumulator_width: usize,
) -> FusedZeroPrefixProofCircuit {
    let mut circuit = Circuit::new();
    let source = circuit.alloc_qreg_bits("fused-zero-control-proof.source", source_width);
    let control = circuit.alloc_qreg("fused-zero-control-proof.control");
    let accumulator = circuit.alloc_qreg_bits(
        "fused-zero-control-proof.accumulator",
        accumulator_width,
    );
    let temporary =
        circuit.alloc_qreg_bits("fused-zero-control-proof.temporary", accumulator_width);
    let source_refs: Vec<&QReg> = source.iter().collect();
    bit_length_lean_allow_zero(&mut circuit, &source_refs, &temporary, false);
    for (length_bit, accumulator_bit) in temporary.iter().zip(&accumulator) {
        circuit.ccx(&control, length_bit, accumulator_bit);
    }
    bit_length_lean_allow_zero(&mut circuit, &source_refs, &temporary, true);
    for lane in temporary {
        circuit.zero_and_free(lane);
    }
    FusedZeroPrefixProofCircuit {
        source_ids: source.iter().map(QReg::id).collect(),
        accumulator_ids: accumulator.iter().map(QReg::id).collect(),
        control_id: Some(control.id()),
        builder: circuit.into_builder(),
    }
}

fn fused_zero_prefix_resources(
    circuit: &FusedZeroPrefixProofCircuit,
) -> FusedZeroPrefixLocalResources {
    use crate::circuit::OperationType;

    let external_qubits = circuit.source_ids.len()
        + circuit.accumulator_ids.len()
        + usize::from(circuit.control_id.is_some());
    FusedZeroPrefixLocalResources {
        extra_qubits: circuit.builder.peak_qubits as usize - external_qubits,
        emitted_ops: circuit.builder.ops.len(),
        emitted_toffoli: circuit
            .builder
            .ops
            .iter()
            .filter(|operation| {
                matches!(operation.kind, OperationType::CCX | OperationType::CCZ)
            })
            .count(),
    }
}

fn assert_fused_zero_prefix_default_stream(
    default: &FusedZeroPrefixProofCircuit,
    legacy: &FusedZeroPrefixProofCircuit,
) {
    assert_eq!(default.source_ids, legacy.source_ids);
    assert_eq!(default.accumulator_ids, legacy.accumulator_ids);
    assert_eq!(default.builder.ops, legacy.builder.ops);
    assert_eq!(default.builder.next_qubit, legacy.builder.next_qubit);
    assert_eq!(default.builder.next_bit, legacy.builder.next_bit);
    assert_eq!(default.builder.active_qubits, legacy.builder.active_qubits);
    assert_eq!(default.builder.peak_qubits, legacy.builder.peak_qubits);
    assert_eq!(default.builder.free_qubits, legacy.builder.free_qubits);
    assert_eq!(
        default.builder.allocation_serial,
        legacy.builder.allocation_serial
    );
}

fn read_fused_zero_prefix_register<R: sha3::digest::XofReader>(
    simulator: &crate::sim::Simulator<'_, R>,
    ids: &[u32],
    shot: usize,
) -> u64 {
    use crate::circuit::QubitId;

    ids.iter().enumerate().fold(0u64, |value, (bit, &id)| {
        value | (((simulator.qubit(QubitId(u64::from(id))) >> shot) & 1) << bit)
    })
}

fn fused_zero_prefix_external_ids(circuit: &FusedZeroPrefixProofCircuit) -> Vec<u32> {
    circuit
        .source_ids
        .iter()
        .chain(circuit.control_id.iter())
        .chain(circuit.accumulator_ids.iter())
        .copied()
        .collect()
}

fn assert_fused_zero_prefix_internal_clean<R: sha3::digest::XofReader>(
    simulator: &crate::sim::Simulator<'_, R>,
    circuit: &FusedZeroPrefixProofCircuit,
    live: u64,
    shots: usize,
    context: &str,
) -> usize {
    use crate::circuit::QubitId;

    let external = fused_zero_prefix_external_ids(circuit);
    let mut internal_lanes = 0usize;
    for id in 0..circuit.builder.next_qubit {
        if !external.contains(&id) {
            assert_eq!(
                simulator.qubit(QubitId(u64::from(id))) & live,
                0,
                "{context} left internal q{id} dirty"
            );
            internal_lanes += 1;
        }
    }
    internal_lanes * shots
}

fn verify_fused_zero_prefix_update_equivalence(
    baseline: &FusedZeroPrefixProofCircuit,
    fused: &FusedZeroPrefixProofCircuit,
    source_width: usize,
    decrement: bool,
) -> (usize, usize, usize) {
    use crate::circuit::QubitId;
    use crate::sim::Simulator;
    use sha3::{
        digest::{ExtendableOutput, Update},
        Shake128,
    };

    assert_eq!(baseline.source_ids, fused.source_ids);
    assert_eq!(baseline.accumulator_ids, fused.accumulator_ids);
    let accumulator_width = fused.accumulator_ids.len();
    let accumulator_modulus = 1u64 << accumulator_width;
    let cases: Vec<(u64, u64)> = (0..(1u64 << source_width))
        .flat_map(|source| (0..accumulator_modulus).map(move |accumulator| (source, accumulator)))
        .collect();
    let mut phase_clean_checks = 0usize;
    let mut ancilla_clean_checks = 0usize;

    for (batch, chunk) in cases.chunks(64).enumerate() {
        let mut baseline_seed = Shake128::default();
        baseline_seed.update(b"q883-fused-zero-update-baseline");
        baseline_seed.update(&(source_width as u64).to_le_bytes());
        baseline_seed.update(&(decrement as u64).to_le_bytes());
        baseline_seed.update(&(batch as u64).to_le_bytes());
        let mut baseline_xof = baseline_seed.finalize_xof();
        let mut baseline_simulator = Simulator::new(
            baseline.builder.next_qubit as usize,
            baseline.builder.next_bit as usize,
            &mut baseline_xof,
        );

        let mut fused_seed = Shake128::default();
        fused_seed.update(b"q883-fused-zero-update-fused");
        fused_seed.update(&(source_width as u64).to_le_bytes());
        fused_seed.update(&(decrement as u64).to_le_bytes());
        fused_seed.update(&(batch as u64).to_le_bytes());
        let mut fused_xof = fused_seed.finalize_xof();
        let mut fused_simulator = Simulator::new(
            fused.builder.next_qubit as usize,
            fused.builder.next_bit as usize,
            &mut fused_xof,
        );

        for (shot, &(source_value, accumulator_value)) in chunk.iter().enumerate() {
            for (bit, (&baseline_id, &fused_id)) in baseline
                .source_ids
                .iter()
                .zip(&fused.source_ids)
                .enumerate()
            {
                if (source_value >> bit) & 1 == 1 {
                    *baseline_simulator.qubit_mut(QubitId(u64::from(baseline_id))) |=
                        1u64 << shot;
                    *fused_simulator.qubit_mut(QubitId(u64::from(fused_id))) |= 1u64 << shot;
                }
            }
            for (bit, (&baseline_id, &fused_id)) in baseline
                .accumulator_ids
                .iter()
                .zip(&fused.accumulator_ids)
                .enumerate()
            {
                if (accumulator_value >> bit) & 1 == 1 {
                    *baseline_simulator.qubit_mut(QubitId(u64::from(baseline_id))) |=
                        1u64 << shot;
                    *fused_simulator.qubit_mut(QubitId(u64::from(fused_id))) |= 1u64 << shot;
                }
            }
        }
        baseline_simulator.apply_iter(baseline.builder.ops.iter());
        fused_simulator.apply_iter(fused.builder.ops.iter());
        let live = if chunk.len() == 64 {
            u64::MAX
        } else {
            (1u64 << chunk.len()) - 1
        };
        assert_eq!(baseline_simulator.phase & live, 0);
        assert_eq!(fused_simulator.phase & live, 0);
        phase_clean_checks += 2 * chunk.len();

        for (shot, &(source_value, accumulator_value)) in chunk.iter().enumerate() {
            let bit_length = if source_value == 0 {
                0
            } else {
                u64::from(64 - source_value.leading_zeros())
            };
            let expected = if decrement {
                accumulator_value.wrapping_sub(bit_length) % accumulator_modulus
            } else {
                accumulator_value.wrapping_add(bit_length) % accumulator_modulus
            };
            let baseline_source = read_fused_zero_prefix_register(
                &baseline_simulator,
                &baseline.source_ids,
                shot,
            );
            let fused_source =
                read_fused_zero_prefix_register(&fused_simulator, &fused.source_ids, shot);
            let baseline_accumulator = read_fused_zero_prefix_register(
                &baseline_simulator,
                &baseline.accumulator_ids,
                shot,
            );
            let fused_accumulator = read_fused_zero_prefix_register(
                &fused_simulator,
                &fused.accumulator_ids,
                shot,
            );
            assert_eq!(baseline_source, source_value);
            assert_eq!(fused_source, source_value);
            assert_eq!(baseline_accumulator, expected);
            assert_eq!(fused_accumulator, expected);
        }
        ancilla_clean_checks += assert_fused_zero_prefix_internal_clean(
            &baseline_simulator,
            baseline,
            live,
            chunk.len(),
            "legacy update",
        );
        ancilla_clean_checks += assert_fused_zero_prefix_internal_clean(
            &fused_simulator,
            fused,
            live,
            chunk.len(),
            "fused update",
        );
    }

    (cases.len(), phase_clean_checks, ancilla_clean_checks)
}

fn verify_fused_zero_prefix_controlled_equivalence(
    baseline: &FusedZeroPrefixProofCircuit,
    fused: &FusedZeroPrefixProofCircuit,
    source_width: usize,
) -> (usize, usize, usize) {
    use crate::circuit::QubitId;
    use crate::sim::Simulator;
    use sha3::{
        digest::{ExtendableOutput, Update},
        Shake128,
    };

    assert_eq!(baseline.source_ids, fused.source_ids);
    assert_eq!(baseline.accumulator_ids, fused.accumulator_ids);
    let baseline_control = baseline.control_id.expect("baseline control");
    let fused_control = fused.control_id.expect("fused control");
    let accumulator_width = fused.accumulator_ids.len();
    let accumulator_modulus = 1u64 << accumulator_width;
    let cases: Vec<(u64, u64, u64)> = (0..(1u64 << source_width))
        .flat_map(|source| {
            (0..=1).flat_map(move |control| {
                (0..accumulator_modulus)
                    .map(move |accumulator| (source, control, accumulator))
            })
        })
        .collect();
    let mut phase_clean_checks = 0usize;
    let mut ancilla_clean_checks = 0usize;

    for (batch, chunk) in cases.chunks(64).enumerate() {
        let mut baseline_seed = Shake128::default();
        baseline_seed.update(b"q883-fused-zero-control-baseline");
        baseline_seed.update(&(source_width as u64).to_le_bytes());
        baseline_seed.update(&(batch as u64).to_le_bytes());
        let mut baseline_xof = baseline_seed.finalize_xof();
        let mut baseline_simulator = Simulator::new(
            baseline.builder.next_qubit as usize,
            baseline.builder.next_bit as usize,
            &mut baseline_xof,
        );

        let mut fused_seed = Shake128::default();
        fused_seed.update(b"q883-fused-zero-control-fused");
        fused_seed.update(&(source_width as u64).to_le_bytes());
        fused_seed.update(&(batch as u64).to_le_bytes());
        let mut fused_xof = fused_seed.finalize_xof();
        let mut fused_simulator = Simulator::new(
            fused.builder.next_qubit as usize,
            fused.builder.next_bit as usize,
            &mut fused_xof,
        );

        for (shot, &(source_value, control_value, accumulator_value)) in
            chunk.iter().enumerate()
        {
            for (bit, (&baseline_id, &fused_id)) in baseline
                .source_ids
                .iter()
                .zip(&fused.source_ids)
                .enumerate()
            {
                if (source_value >> bit) & 1 == 1 {
                    *baseline_simulator.qubit_mut(QubitId(u64::from(baseline_id))) |=
                        1u64 << shot;
                    *fused_simulator.qubit_mut(QubitId(u64::from(fused_id))) |= 1u64 << shot;
                }
            }
            if control_value == 1 {
                *baseline_simulator.qubit_mut(QubitId(u64::from(baseline_control))) |=
                    1u64 << shot;
                *fused_simulator.qubit_mut(QubitId(u64::from(fused_control))) |= 1u64 << shot;
            }
            for (bit, (&baseline_id, &fused_id)) in baseline
                .accumulator_ids
                .iter()
                .zip(&fused.accumulator_ids)
                .enumerate()
            {
                if (accumulator_value >> bit) & 1 == 1 {
                    *baseline_simulator.qubit_mut(QubitId(u64::from(baseline_id))) |=
                        1u64 << shot;
                    *fused_simulator.qubit_mut(QubitId(u64::from(fused_id))) |= 1u64 << shot;
                }
            }
        }
        baseline_simulator.apply_iter(baseline.builder.ops.iter());
        fused_simulator.apply_iter(fused.builder.ops.iter());
        let live = if chunk.len() == 64 {
            u64::MAX
        } else {
            (1u64 << chunk.len()) - 1
        };
        assert_eq!(baseline_simulator.phase & live, 0);
        assert_eq!(fused_simulator.phase & live, 0);
        phase_clean_checks += 2 * chunk.len();

        for (shot, &(source_value, control_value, accumulator_value)) in
            chunk.iter().enumerate()
        {
            let bit_length = if source_value == 0 {
                0
            } else {
                u64::from(64 - source_value.leading_zeros())
            };
            let expected = accumulator_value
                ^ if control_value == 1 { bit_length } else { 0 };
            assert_eq!(
                read_fused_zero_prefix_register(
                    &baseline_simulator,
                    &baseline.source_ids,
                    shot,
                ),
                source_value
            );
            assert_eq!(
                read_fused_zero_prefix_register(&fused_simulator, &fused.source_ids, shot),
                source_value
            );
            assert_eq!(
                (baseline_simulator.qubit(QubitId(u64::from(baseline_control))) >> shot) & 1,
                control_value
            );
            assert_eq!(
                (fused_simulator.qubit(QubitId(u64::from(fused_control))) >> shot) & 1,
                control_value
            );
            assert_eq!(
                read_fused_zero_prefix_register(
                    &baseline_simulator,
                    &baseline.accumulator_ids,
                    shot,
                ),
                expected
            );
            assert_eq!(
                read_fused_zero_prefix_register(
                    &fused_simulator,
                    &fused.accumulator_ids,
                    shot,
                ),
                expected
            );
        }
        ancilla_clean_checks += assert_fused_zero_prefix_internal_clean(
            &baseline_simulator,
            baseline,
            live,
            chunk.len(),
            "legacy controlled",
        );
        ancilla_clean_checks += assert_fused_zero_prefix_internal_clean(
            &fused_simulator,
            fused,
            live,
            chunk.len(),
            "fused controlled",
        );
    }

    (cases.len(), phase_clean_checks, ancilla_clean_checks)
}

fn maximize_fused_zero_prefix_resources(
    maximum: &mut FusedZeroPrefixLocalResources,
    candidate: FusedZeroPrefixLocalResources,
) {
    maximum.extra_qubits = maximum.extra_qubits.max(candidate.extra_qubits);
    maximum.emitted_ops = maximum.emitted_ops.max(candidate.emitted_ops);
    maximum.emitted_toffoli = maximum.emitted_toffoli.max(candidate.emitted_toffoli);
}

/// Exhaustively prove that the feature-gated final complemented-prefix term
/// replaces the historical zero-flag correction without changing semantics.
/// Widths zero through eight cover every source, every five-bit accumulator,
/// both arithmetic directions, and both values of an external control in the
/// production compute/XOR/uncompute usage.
#[doc(hidden)]
pub fn fused_zero_prefix_bit_length_roundtrip_check(
) -> FusedZeroPrefixBitLengthProofReport {
    const MAX_SOURCE_WIDTH: usize = 8;
    const ACCUMULATOR_WIDTH: usize = 5;
    const LOCAL_WIDTH: usize = 259;
    const LOCAL_ACCUMULATOR_WIDTH: usize = 10;

    let mut update_cases_checked = 0usize;
    let mut controlled_cases_checked = 0usize;
    let mut default_stream_equivalence_checks = 0usize;
    let mut baseline_equivalence_checks = 0usize;
    let mut phase_clean_checks = 0usize;
    let mut ancilla_clean_checks = 0usize;
    let mut maximum_baseline = FusedZeroPrefixLocalResources::default();
    let mut maximum_fused = FusedZeroPrefixLocalResources::default();

    for source_width in 0..=MAX_SOURCE_WIDTH {
        for decrement in [false, true] {
            configure_fused_zero_prefix_proof(false);
            let baseline = build_fused_zero_prefix_update(
                source_width,
                ACCUMULATOR_WIDTH,
                decrement,
                false,
            );
            let legacy = build_fused_zero_prefix_update(
                source_width,
                ACCUMULATOR_WIDTH,
                decrement,
                true,
            );
            assert_fused_zero_prefix_default_stream(&baseline, &legacy);
            default_stream_equivalence_checks += 1;

            configure_fused_zero_prefix_proof(true);
            let fused = build_fused_zero_prefix_update(
                source_width,
                ACCUMULATOR_WIDTH,
                decrement,
                false,
            );
            let baseline_resources = fused_zero_prefix_resources(&baseline);
            let fused_resources = fused_zero_prefix_resources(&fused);
            assert!(fused_resources.extra_qubits <= baseline_resources.extra_qubits);
            assert!(fused_resources.emitted_ops <= baseline_resources.emitted_ops);
            assert!(fused_resources.emitted_toffoli <= baseline_resources.emitted_toffoli);
            maximize_fused_zero_prefix_resources(&mut maximum_baseline, baseline_resources);
            maximize_fused_zero_prefix_resources(&mut maximum_fused, fused_resources);

            let (cases, phase_checks, ancilla_checks) =
                verify_fused_zero_prefix_update_equivalence(
                    &baseline,
                    &fused,
                    source_width,
                    decrement,
                );
            update_cases_checked += cases;
            baseline_equivalence_checks += cases;
            phase_clean_checks += phase_checks;
            ancilla_clean_checks += ancilla_checks;
        }

        configure_fused_zero_prefix_proof(false);
        let baseline_controlled =
            build_fused_zero_prefix_controlled_xor(source_width, ACCUMULATOR_WIDTH);
        configure_fused_zero_prefix_proof(true);
        let fused_controlled =
            build_fused_zero_prefix_controlled_xor(source_width, ACCUMULATOR_WIDTH);
        let (cases, phase_checks, ancilla_checks) =
            verify_fused_zero_prefix_controlled_equivalence(
                &baseline_controlled,
                &fused_controlled,
                source_width,
            );
        controlled_cases_checked += cases;
        baseline_equivalence_checks += cases;
        phase_clean_checks += phase_checks;
        ancilla_clean_checks += ancilla_checks;
    }

    configure_fused_zero_prefix_proof(false);
    let local_add_baseline = fused_zero_prefix_resources(&build_fused_zero_prefix_update(
        LOCAL_WIDTH,
        LOCAL_ACCUMULATOR_WIDTH,
        false,
        false,
    ));
    let local_sub_baseline = fused_zero_prefix_resources(&build_fused_zero_prefix_update(
        LOCAL_WIDTH,
        LOCAL_ACCUMULATOR_WIDTH,
        true,
        false,
    ));
    configure_fused_zero_prefix_proof(true);
    let local_add_fused = fused_zero_prefix_resources(&build_fused_zero_prefix_update(
        LOCAL_WIDTH,
        LOCAL_ACCUMULATOR_WIDTH,
        false,
        false,
    ));
    let local_sub_fused = fused_zero_prefix_resources(&build_fused_zero_prefix_update(
        LOCAL_WIDTH,
        LOCAL_ACCUMULATOR_WIDTH,
        true,
        false,
    ));
    assert!(local_add_fused.extra_qubits < local_add_baseline.extra_qubits);
    assert!(local_add_fused.emitted_ops < local_add_baseline.emitted_ops);
    assert!(local_add_fused.emitted_toffoli < local_add_baseline.emitted_toffoli);
    assert!(local_sub_fused.extra_qubits < local_sub_baseline.extra_qubits);
    assert!(local_sub_fused.emitted_ops < local_sub_baseline.emitted_ops);
    assert!(local_sub_fused.emitted_toffoli < local_sub_baseline.emitted_toffoli);

    configure_fused_zero_prefix_proof(false);
    let scheduled_baseline = profile_reference_scheduled_inversion();
    configure_fused_zero_prefix_proof(true);
    let scheduled_fused = profile_reference_scheduled_inversion();
    assert_eq!(scheduled_fused.steps, scheduled_baseline.steps);
    assert!(scheduled_fused.inversion_peak_qubits <= scheduled_baseline.inversion_peak_qubits);
    assert!(scheduled_fused.emitted_ops < scheduled_baseline.emitted_ops);
    assert!(scheduled_fused.emitted_toffoli < scheduled_baseline.emitted_toffoli);

    configure_fused_zero_prefix_proof(true);
    FusedZeroPrefixBitLengthProofReport {
        widths_checked: MAX_SOURCE_WIDTH + 1,
        accumulator_width: ACCUMULATOR_WIDTH,
        update_cases_checked,
        controlled_cases_checked,
        default_stream_equivalence_checks,
        baseline_equivalence_checks,
        phase_clean_checks,
        ancilla_clean_checks,
        maximum_baseline,
        maximum_fused,
        local_width: LOCAL_WIDTH,
        local_accumulator_width: LOCAL_ACCUMULATOR_WIDTH,
        local_add_baseline,
        local_add_fused,
        local_sub_baseline,
        local_sub_fused,
        scheduled_steps: scheduled_baseline.steps,
        scheduled_baseline_inversion_peak_qubits: scheduled_baseline.inversion_peak_qubits,
        scheduled_fused_inversion_peak_qubits: scheduled_fused.inversion_peak_qubits,
        scheduled_baseline_emitted_ops: scheduled_baseline.emitted_ops,
        scheduled_fused_emitted_ops: scheduled_fused.emitted_ops,
        scheduled_baseline_emitted_toffoli: scheduled_baseline.emitted_toffoli,
        scheduled_fused_emitted_toffoli: scheduled_fused.emitted_toffoli,
        scheduled_baseline_emitted_hmr: scheduled_baseline.emitted_hmr,
        scheduled_fused_emitted_hmr: scheduled_fused.emitted_hmr,
        scheduled_baseline_emitted_resets: scheduled_baseline.emitted_resets,
        scheduled_fused_emitted_resets: scheduled_fused.emitted_resets,
    }
}

fn controlled_xor_dynamic_prefix_bit_length_quadratic(
    circ: &mut Circuit,
    control: &QReg,
    right_length: &[QReg],
    work: &[QReg],
    output: &[QReg],
) {
    let flag = circ.alloc_qreg("rs.dynamic-prefix.flag");
    let chain = circ.alloc_qreg_bits(
        "rs.dynamic-prefix.chain",
        right_length.len().saturating_sub(1),
    );
    let temporary = circ.alloc_qreg_bits("rs.dynamic-prefix.length", output.len());
    for boundary in 0..=work.len() {
        equality_flag(circ, control, right_length, boundary, &flag, &chain);
        let source: Vec<&QReg> = work[..work.len() - boundary].iter().collect();
        bit_length_lean_allow_zero(circ, &source, &temporary, false);
        for (source_bit, output_bit) in temporary.iter().zip(output) {
            circ.ccx(&flag, source_bit, output_bit);
        }
        bit_length_lean_allow_zero(circ, &source, &temporary, true);
        equality_flag(circ, control, right_length, boundary, &flag, &chain);
    }
    for lane in temporary {
        circ.zero_and_free(lane);
    }
    for lane in chain {
        circ.zero_and_free(lane);
    }
    circ.zero_and_free(flag);
}

fn controlled_xor_dynamic_suffix_bit_length_quadratic(
    circ: &mut Circuit,
    control: &QReg,
    left_length: &[QReg],
    work: &[QReg],
    output: &[QReg],
) {
    let flag = circ.alloc_qreg("rs.dynamic-suffix.flag");
    let chain = circ.alloc_qreg_bits(
        "rs.dynamic-suffix.chain",
        left_length.len().saturating_sub(1),
    );
    let temporary = circ.alloc_qreg_bits("rs.dynamic-suffix.length", output.len());
    for boundary in 0..=work.len() {
        equality_flag(circ, control, left_length, boundary, &flag, &chain);
        let start = (boundary + 1).min(work.len());
        let source: Vec<&QReg> = work[start..].iter().rev().collect();
        bit_length_lean_allow_zero(circ, &source, &temporary, false);
        for (source_bit, output_bit) in temporary.iter().zip(output) {
            circ.ccx(&flag, source_bit, output_bit);
        }
        bit_length_lean_allow_zero(circ, &source, &temporary, true);
        equality_flag(circ, control, left_length, boundary, &flag, &chain);
    }
    for lane in temporary {
        circ.zero_and_free(lane);
    }
    for lane in chain {
        circ.zero_and_free(lane);
    }
    circ.zero_and_free(flag);
}

/// XOR `control AND max(bitlen(source) - boundary, 0)` into `output`.
///
/// The extra high lane makes the modular subtraction a signed subtraction. A
/// boundary occupies `output.len()` bits and `source.len()` is constrained to
/// the positive signed range, so the top lane is exactly the negative
/// predicate. The following add is the exact inverse of the subtract and
/// restores all arithmetic scratch, including the overflow lane.
fn controlled_xor_saturating_bit_length_difference(
    circ: &mut Circuit,
    control: &QReg,
    boundary: &[QReg],
    source: &[&QReg],
    output: &[QReg],
) {
    assert_eq!(boundary.len(), output.len());
    assert!(!output.is_empty());
    let signed_width = output.len() + 1;
    assert!(
        source.len() <= (1usize << (signed_width - 1)) - 1,
        "signed bit-length workspace is too narrow"
    );

    let length = circ.alloc_qreg_bits("rs.rotated-bitlen.length", signed_width);
    bit_length_lean_allow_zero(circ, source, &length, false);

    let boundary_copy = circ.alloc_qreg_bits("rs.rotated-bitlen.boundary", signed_width);
    for (source_bit, target_bit) in boundary.iter().zip(&boundary_copy) {
        circ.cx(source_bit, target_bit);
    }
    let carry = circ.alloc_qreg("rs.rotated-bitlen.carry");
    let overflow = circ.alloc_qreg("rs.rotated-bitlen.overflow");
    cuccaro_sub_mod_2n(circ, &boundary_copy, &length, &carry, &overflow);

    let sign = &length[signed_width - 1];
    let enabled = circ.alloc_qreg("rs.rotated-bitlen.enabled");
    circ.x(sign);
    circ.ccx(control, sign, &enabled);
    for (difference_bit, output_bit) in length.iter().zip(output) {
        circ.ccx(&enabled, difference_bit, output_bit);
    }
    circ.ccx(control, sign, &enabled);
    circ.x(sign);
    circ.zero_and_free(enabled);

    cuccaro_add_mod_2n(circ, &boundary_copy, &length, &carry, &overflow);
    circ.zero_and_free(overflow);
    circ.zero_and_free(carry);
    for (source_bit, target_bit) in boundary.iter().zip(&boundary_copy) {
        circ.cx(source_bit, target_bit);
    }
    for lane in boundary_copy {
        circ.zero_and_free(lane);
    }

    bit_length_lean_allow_zero(circ, source, &length, true);
    for lane in length {
        circ.zero_and_free(lane);
    }
}

/// The rightmost `right_length` lanes contain the packed `r'` component.
/// Rotating them to the low end places `t'` immediately above that boundary,
/// so a single full-register bit length reveals `right_length + bitlen(t')`
/// whenever `t'` is nonzero. Saturation maps the zero case to zero exactly.
fn controlled_xor_rotated_prefix_bit_length(
    circ: &mut Circuit,
    control: &QReg,
    right_length: &[QReg],
    work: &[QReg],
    output: &[QReg],
) {
    let view: Vec<&QReg> = work.iter().collect();
    variable_rotate_high_refs(circ, right_length, &view);
    controlled_xor_saturating_bit_length_difference(circ, control, right_length, &view, output);
    variable_rotate_low_refs(circ, right_length, &view);
}

/// XOR the raw `t'` bit length into `output` while Work2 is in its shifted
/// coefficient-update layout. The physical rotation includes `l_s`, but the
/// packed remainder boundary does not, so those two lengths are intentionally
/// distinct.
fn controlled_xor_raw_t_prime_bit_length(
    circ: &mut Circuit,
    control: &QReg,
    l_s: &[QReg],
    l_r_prime: &[QReg],
    work2: &[QReg],
    output: &[QReg],
) {
    assert_eq!(l_s.len(), l_r_prime.len());
    assert_eq!(l_s.len(), output.len());
    let rotation = circ.alloc_qreg_bits("rs.raw-t-prime.rotation", output.len());
    for (source, destination) in l_r_prime.iter().zip(&rotation) {
        circ.cx(source, destination);
    }
    let carry = circ.alloc_qreg("rs.raw-t-prime.rotation-carry");
    let overflow = circ.alloc_qreg("rs.raw-t-prime.rotation-overflow");
    cuccaro_add_mod_2n(circ, l_s, &rotation, &carry, &overflow);

    let view: Vec<&QReg> = work2.iter().collect();
    variable_rotate_high_refs(circ, &rotation, &view);
    controlled_xor_saturating_bit_length_difference(circ, control, l_r_prime, &view, output);
    variable_rotate_low_refs(circ, &rotation, &view);

    cuccaro_sub_mod_2n(circ, l_s, &rotation, &carry, &overflow);
    circ.zero_and_free(overflow);
    circ.zero_and_free(carry);
    for (source, destination) in l_r_prime.iter().zip(&rotation) {
        circ.cx(source, destination);
    }
    free_clean(circ, rotation);
}

/// In the reversed Work1 view, the packed `t` component occupies the
/// rightmost `left_length` lanes. The separator is zero and remains above the
/// rotated remainder, hence `bitlen(rotated) - left_length = bitlen(r)`.
fn controlled_xor_rotated_suffix_bit_length(
    circ: &mut Circuit,
    control: &QReg,
    left_length: &[QReg],
    work: &[QReg],
    output: &[QReg],
) {
    let reversed: Vec<&QReg> = work.iter().rev().collect();
    variable_rotate_high_refs(circ, left_length, &reversed);
    controlled_xor_saturating_bit_length_difference(circ, control, left_length, &reversed, output);
    variable_rotate_low_refs(circ, left_length, &reversed);
}

#[allow(clippy::too_many_arguments)]
fn conditional_work_and_length_swap_quadratic_oracle(
    circ: &mut Circuit,
    control: &QReg,
    iteration_parity: &QReg,
    work1: &[QReg],
    work2: &[QReg],
    l_t: &[QReg],
    l_t_prime: &[QReg],
    l_q: &[QReg],
    l_s: &[QReg],
    l_r_prime: &[QReg],
) {
    assert_eq!(work1.len(), work2.len());
    assert_eq!(l_t.len(), l_t_prime.len());
    assert_eq!(l_t.len(), l_r_prime.len());
    let old_r_length = circ.alloc_qreg_bits("rs.swap-length.old-lrp", l_r_prime.len());

    controlled_xor_dynamic_suffix_bit_length_quadratic(circ, control, l_t, work1, &old_r_length);
    for (current, next) in l_t.iter().zip(l_t_prime) {
        circ.cswap(control, current, next);
    }
    for (current, next) in l_r_prime.iter().zip(&old_r_length) {
        circ.cswap(control, current, next);
    }
    controlled_swap_registers(circ, control, work1, work2);
    controlled_xor_dynamic_suffix_bit_length_quadratic(circ, control, l_t, work1, &old_r_length);
    for lane in old_r_length {
        circ.zero_and_free(lane);
    }
    circ.cx(control, iteration_parity);
    let _ = (l_q, l_s);
}

pub fn conditional_work_and_length_swap(
    circ: &mut Circuit,
    control: &QReg,
    iteration_parity: &QReg,
    work1: &[QReg],
    work2: &[QReg],
    l_t: &[QReg],
    l_t_prime: &[QReg],
    l_q: &[QReg],
    l_s: &[QReg],
    l_r_prime: &[QReg],
    _t_window: (usize, usize),
    _r_window: (usize, usize),
    _scratch: &[QReg],
) {
    assert_eq!(work1.len(), work2.len());
    assert_eq!(l_t.len(), l_r_prime.len());
    assert_eq!(l_t.len(), l_t_prime.len());
    let old_r_length = circ.alloc_qreg_bits("rs.swap-length.old-lrp", l_r_prime.len());

    controlled_xor_rotated_suffix_bit_length(circ, control, l_t, work1, &old_r_length);

    for (current, next) in l_t.iter().zip(l_t_prime) {
        circ.cswap(control, current, next);
    }
    for (current, next) in l_r_prime.iter().zip(&old_r_length) {
        circ.cswap(control, current, next);
    }
    controlled_swap_registers(circ, control, work1, work2);

    controlled_xor_rotated_suffix_bit_length(circ, control, l_t, work1, &old_r_length);
    for lane in old_r_length {
        circ.zero_and_free(lane);
    }
    circ.cx(control, iteration_parity);

    let _ = (l_q, l_s);
}

#[allow(clippy::too_many_arguments)]
pub fn conditional_work_and_length_swap_inverse(
    circ: &mut Circuit,
    control: &QReg,
    iteration_parity: &QReg,
    work1: &[QReg],
    work2: &[QReg],
    l_t: &[QReg],
    l_t_prime: &[QReg],
    l_q: &[QReg],
    l_s: &[QReg],
    l_r_prime: &[QReg],
    t_window: (usize, usize),
    r_window: (usize, usize),
    scratch: &[QReg],
) {
    conditional_work_and_length_swap(
        circ,
        control,
        iteration_parity,
        work1,
        work2,
        l_t,
        l_t_prime,
        l_q,
        l_s,
        l_r_prime,
        t_window,
        r_window,
        scratch,
    );
}

fn compute_remainder_operation(
    circ: &mut Circuit,
    phase1: &QReg,
    l_r_prime: &[QReg],
    scratch: &RemainderScratch<'_>,
) {
    compute_nonzero(circ, l_r_prime, scratch.nonzero, scratch.nonzero_chain);
    circ.x(phase1);
    circ.ccx(phase1, scratch.nonzero, scratch.operation);
    circ.x(phase1);
}

fn uncompute_remainder_operation(
    circ: &mut Circuit,
    phase1: &QReg,
    l_r_prime: &[QReg],
    scratch: &RemainderScratch<'_>,
) {
    circ.x(phase1);
    circ.ccx(phase1, scratch.nonzero, scratch.operation);
    circ.x(phase1);
    uncompute_nonzero(circ, l_r_prime, scratch.nonzero, scratch.nonzero_chain);
}

fn toggle_remainder_add_enable(
    circ: &mut Circuit,
    operation: &QReg,
    phase2: &QReg,
    sign: &QReg,
    phase_sign: &QReg,
    enable: &QReg,
) {
    // enable ^= operation AND !(phase2 AND sign), restoring phase_sign.
    circ.ccx(phase2, sign, phase_sign);
    circ.x(phase_sign);
    circ.ccx(operation, phase_sign, enable);
    circ.x(phase_sign);
    circ.ccx(phase2, sign, phase_sign);
}

fn prepare_remainder_range(
    circ: &mut Circuit,
    total_work_width: usize,
    window_upper: usize,
    l_t: &[QReg],
    l_q: &[QReg],
    l_s: &[QReg],
    scratch: &RemainderScratch<'_>,
) {
    assert_eq!(l_t.len(), l_q.len());
    assert_eq!(l_t.len(), l_s.len());
    assert!(window_upper <= total_work_width);
    cuccaro_add_mod_2n(
        circ,
        l_t,
        l_q,
        scratch.length_carry,
        scratch.length_overflow,
    );
    // At global position j=K, sign(l_q) is one iff l_t+l_q <= K-2.
    sub_const_mod_2n(circ, l_q, window_upper - 1, scratch.constant);
    // sign(l_s) is one iff l_s <= N-K.
    sub_const_mod_2n(
        circ,
        l_s,
        total_work_width - window_upper + 1,
        scratch.constant,
    );
}

fn unprepare_remainder_range(
    circ: &mut Circuit,
    total_work_width: usize,
    window_upper: usize,
    l_t: &[QReg],
    l_q: &[QReg],
    l_s: &[QReg],
    scratch: &RemainderScratch<'_>,
) {
    add_const_mod_2n(
        circ,
        l_s,
        total_work_width - window_upper + 1,
        scratch.constant,
    );
    add_const_mod_2n(circ, l_q, window_upper - 1, scratch.constant);
    cuccaro_sub_mod_2n(
        circ,
        l_t,
        l_q,
        scratch.length_carry,
        scratch.length_overflow,
    );
}

fn toggle_remainder_range_active(
    circ: &mut Circuit,
    control: &QReg,
    l_q: &[QReg],
    l_s: &[QReg],
    scratch: &RemainderScratch<'_>,
) {
    let controls = [
        control,
        l_q.last().expect("nonempty l_q"),
        l_s.last().expect("nonempty l_s"),
    ];
    multi_controlled_x_vchain(
        circ,
        &controls,
        scratch.active,
        std::slice::from_ref(scratch.tmp),
    );
}

fn walk_remainder_range_down(
    circ: &mut Circuit,
    l_q: &[QReg],
    l_s: &[QReg],
    scratch: &RemainderScratch<'_>,
) {
    decrement_mod_2n(circ, l_s, scratch.walk);
    increment_mod_2n(circ, l_q, scratch.walk);
}

fn walk_remainder_range_up(
    circ: &mut Circuit,
    l_q: &[QReg],
    l_s: &[QReg],
    scratch: &RemainderScratch<'_>,
) {
    increment_mod_2n(circ, l_s, scratch.walk);
    decrement_mod_2n(circ, l_q, scratch.walk);
}

#[allow(clippy::too_many_arguments)]
pub fn remainder_sub_window(
    circ: &mut Circuit,
    total_work_width: usize,
    window_upper: usize,
    phase1: &QReg,
    sign: &QReg,
    work1: &[QReg],
    work2: &[QReg],
    l_t: &[QReg],
    l_q: &[QReg],
    l_s: &[QReg],
    l_r_prime: &[QReg],
    scratch: &[QReg],
) {
    assert_eq!(work1.len(), work2.len());
    let scratch = split_remainder_scratch(scratch, l_t.len(), l_r_prime.len());
    compute_remainder_operation(circ, phase1, l_r_prime, &scratch);
    prepare_remainder_range(
        circ,
        total_work_width,
        window_upper,
        l_t,
        l_q,
        l_s,
        &scratch,
    );
    for index in (0..work1.len()).rev() {
        toggle_remainder_range_active(circ, scratch.operation, l_q, l_s, &scratch);
        circ.ccx(scratch.active, &work2[index], &work1[index]);
        circ.ccx(scratch.active, scratch.carry, &work2[index]);
        multi_controlled_x_vchain(
            circ,
            &[scratch.active, &work2[index], &work1[index]],
            scratch.carry,
            std::slice::from_ref(scratch.tmp),
        );
        toggle_remainder_range_active(circ, scratch.operation, l_q, l_s, &scratch);
        if index != 0 {
            walk_remainder_range_down(circ, l_q, l_s, &scratch);
        }
    }
    circ.ccx(scratch.operation, scratch.carry, sign);
    for index in 0..work1.len() {
        toggle_remainder_range_active(circ, scratch.operation, l_q, l_s, &scratch);
        multi_controlled_x_vchain(
            circ,
            &[scratch.active, &work2[index], &work1[index]],
            scratch.carry,
            std::slice::from_ref(scratch.tmp),
        );
        circ.ccx(scratch.active, scratch.carry, &work2[index]);
        circ.ccx(scratch.active, scratch.carry, &work1[index]);
        toggle_remainder_range_active(circ, scratch.operation, l_q, l_s, &scratch);
        if index + 1 != work1.len() {
            walk_remainder_range_up(circ, l_q, l_s, &scratch);
        }
    }
    unprepare_remainder_range(
        circ,
        total_work_width,
        window_upper,
        l_t,
        l_q,
        l_s,
        &scratch,
    );
    uncompute_remainder_operation(circ, phase1, l_r_prime, &scratch);
}

#[allow(clippy::too_many_arguments)]
pub fn remainder_sub_window_inverse(
    circ: &mut Circuit,
    total_work_width: usize,
    window_upper: usize,
    phase1: &QReg,
    sign: &QReg,
    work1: &[QReg],
    work2: &[QReg],
    l_t: &[QReg],
    l_q: &[QReg],
    l_s: &[QReg],
    l_r_prime: &[QReg],
    scratch: &[QReg],
) {
    assert_eq!(work1.len(), work2.len());
    let scratch = split_remainder_scratch(scratch, l_t.len(), l_r_prime.len());
    compute_remainder_operation(circ, phase1, l_r_prime, &scratch);
    prepare_remainder_range(
        circ,
        total_work_width,
        window_upper,
        l_t,
        l_q,
        l_s,
        &scratch,
    );
    for index in (0..work1.len()).rev() {
        if index + 1 != work1.len() {
            walk_remainder_range_down(circ, l_q, l_s, &scratch);
        }
        toggle_remainder_range_active(circ, scratch.operation, l_q, l_s, &scratch);
        circ.ccx(scratch.active, scratch.carry, &work1[index]);
        circ.ccx(scratch.active, scratch.carry, &work2[index]);
        multi_controlled_x_vchain(
            circ,
            &[scratch.active, &work2[index], &work1[index]],
            scratch.carry,
            std::slice::from_ref(scratch.tmp),
        );
        toggle_remainder_range_active(circ, scratch.operation, l_q, l_s, &scratch);
    }
    circ.ccx(scratch.operation, scratch.carry, sign);
    for index in 0..work1.len() {
        if index != 0 {
            walk_remainder_range_up(circ, l_q, l_s, &scratch);
        }
        toggle_remainder_range_active(circ, scratch.operation, l_q, l_s, &scratch);
        multi_controlled_x_vchain(
            circ,
            &[scratch.active, &work2[index], &work1[index]],
            scratch.carry,
            std::slice::from_ref(scratch.tmp),
        );
        circ.ccx(scratch.active, scratch.carry, &work2[index]);
        circ.ccx(scratch.active, &work2[index], &work1[index]);
        toggle_remainder_range_active(circ, scratch.operation, l_q, l_s, &scratch);
    }
    unprepare_remainder_range(
        circ,
        total_work_width,
        window_upper,
        l_t,
        l_q,
        l_s,
        &scratch,
    );
    uncompute_remainder_operation(circ, phase1, l_r_prime, &scratch);
}

#[allow(clippy::too_many_arguments)]
pub fn remainder_add_window(
    circ: &mut Circuit,
    total_work_width: usize,
    window_upper: usize,
    phase1: &QReg,
    phase2: &QReg,
    sign: &QReg,
    work1: &[QReg],
    work2: &[QReg],
    l_t: &[QReg],
    l_q: &[QReg],
    l_s: &[QReg],
    l_r_prime: &[QReg],
    scratch: &[QReg],
) {
    assert_eq!(work1.len(), work2.len());
    let scratch = split_remainder_scratch(scratch, l_t.len(), l_r_prime.len());
    compute_remainder_operation(circ, phase1, l_r_prime, &scratch);
    toggle_remainder_add_enable(
        circ,
        scratch.operation,
        phase2,
        sign,
        scratch.phase_sign,
        scratch.enable,
    );
    prepare_remainder_range(
        circ,
        total_work_width,
        window_upper,
        l_t,
        l_q,
        l_s,
        &scratch,
    );
    for index in (0..work1.len()).rev() {
        toggle_remainder_range_active(circ, scratch.enable, l_q, l_s, &scratch);
        circ.ccx(scratch.active, scratch.carry, &work1[index]);
        circ.ccx(scratch.active, scratch.carry, &work2[index]);
        multi_controlled_x_vchain(
            circ,
            &[scratch.active, &work2[index], &work1[index]],
            scratch.carry,
            std::slice::from_ref(scratch.tmp),
        );
        toggle_remainder_range_active(circ, scratch.enable, l_q, l_s, &scratch);
        if index != 0 {
            walk_remainder_range_down(circ, l_q, l_s, &scratch);
        }
    }
    for index in 0..work1.len() {
        toggle_remainder_range_active(circ, scratch.enable, l_q, l_s, &scratch);
        multi_controlled_x_vchain(
            circ,
            &[scratch.active, &work2[index], &work1[index]],
            scratch.carry,
            std::slice::from_ref(scratch.tmp),
        );
        circ.ccx(scratch.active, scratch.carry, &work2[index]);
        circ.ccx(scratch.active, &work2[index], &work1[index]);
        toggle_remainder_range_active(circ, scratch.enable, l_q, l_s, &scratch);
        if index + 1 != work1.len() {
            walk_remainder_range_up(circ, l_q, l_s, &scratch);
        }
    }
    unprepare_remainder_range(
        circ,
        total_work_width,
        window_upper,
        l_t,
        l_q,
        l_s,
        &scratch,
    );
    toggle_remainder_add_enable(
        circ,
        scratch.operation,
        phase2,
        sign,
        scratch.phase_sign,
        scratch.enable,
    );
    uncompute_remainder_operation(circ, phase1, l_r_prime, &scratch);
}

#[allow(clippy::too_many_arguments)]
pub fn remainder_add_window_inverse(
    circ: &mut Circuit,
    total_work_width: usize,
    window_upper: usize,
    phase1: &QReg,
    phase2: &QReg,
    sign: &QReg,
    work1: &[QReg],
    work2: &[QReg],
    l_t: &[QReg],
    l_q: &[QReg],
    l_s: &[QReg],
    l_r_prime: &[QReg],
    scratch: &[QReg],
) {
    assert_eq!(work1.len(), work2.len());
    let scratch = split_remainder_scratch(scratch, l_t.len(), l_r_prime.len());
    compute_remainder_operation(circ, phase1, l_r_prime, &scratch);
    toggle_remainder_add_enable(
        circ,
        scratch.operation,
        phase2,
        sign,
        scratch.phase_sign,
        scratch.enable,
    );
    prepare_remainder_range(
        circ,
        total_work_width,
        window_upper,
        l_t,
        l_q,
        l_s,
        &scratch,
    );
    for index in (0..work1.len()).rev() {
        if index + 1 != work1.len() {
            walk_remainder_range_down(circ, l_q, l_s, &scratch);
        }
        toggle_remainder_range_active(circ, scratch.enable, l_q, l_s, &scratch);
        circ.ccx(scratch.active, &work2[index], &work1[index]);
        circ.ccx(scratch.active, scratch.carry, &work2[index]);
        multi_controlled_x_vchain(
            circ,
            &[scratch.active, &work2[index], &work1[index]],
            scratch.carry,
            std::slice::from_ref(scratch.tmp),
        );
        toggle_remainder_range_active(circ, scratch.enable, l_q, l_s, &scratch);
    }
    for index in 0..work1.len() {
        if index != 0 {
            walk_remainder_range_up(circ, l_q, l_s, &scratch);
        }
        toggle_remainder_range_active(circ, scratch.enable, l_q, l_s, &scratch);
        multi_controlled_x_vchain(
            circ,
            &[scratch.active, &work2[index], &work1[index]],
            scratch.carry,
            std::slice::from_ref(scratch.tmp),
        );
        circ.ccx(scratch.active, scratch.carry, &work2[index]);
        circ.ccx(scratch.active, scratch.carry, &work1[index]);
        toggle_remainder_range_active(circ, scratch.enable, l_q, l_s, &scratch);
    }
    unprepare_remainder_range(
        circ,
        total_work_width,
        window_upper,
        l_t,
        l_q,
        l_s,
        &scratch,
    );
    toggle_remainder_add_enable(
        circ,
        scratch.operation,
        phase2,
        sign,
        scratch.phase_sign,
        scratch.enable,
    );
    uncompute_remainder_operation(circ, phase1, l_r_prime, &scratch);
}

fn remainder_phase_sign_flip(
    circ: &mut Circuit,
    phase1: &QReg,
    phase2: &QReg,
    sign: &QReg,
    l_r_prime: &[QReg],
    scratch: &[QReg],
) {
    assert!(scratch.len() >= l_r_prime.len().max(2));
    let nonzero = &scratch[0];
    let operation = &scratch[1];
    let chain = &scratch[2..];
    compute_nonzero(circ, l_r_prime, nonzero, chain);
    circ.x(phase1);
    circ.ccx(phase1, nonzero, operation);
    circ.x(phase1);
    circ.ccx(operation, phase2, sign);
    circ.x(phase1);
    circ.ccx(phase1, nonzero, operation);
    circ.x(phase1);
    uncompute_nonzero(circ, l_r_prime, nonzero, chain);
}

fn free_clean(circ: &mut Circuit, register: Vec<QReg>) {
    for lane in register {
        circ.zero_and_free(lane);
    }
}

fn toggle_initial_coefficient_enable(
    circ: &mut Circuit,
    phase1: &QReg,
    phase2: &QReg,
    sign: &QReg,
    enable: &QReg,
    chain: &[QReg],
) {
    // phase1 AND (phase2 OR !sign), split into disjoint terms.
    circ.ccx(phase1, phase2, enable);
    circ.x(phase2);
    circ.x(sign);
    multi_controlled_x_vchain(circ, &[phase1, phase2, sign], enable, chain);
    circ.x(sign);
    circ.x(phase2);
}

fn toggle_output_coefficient_enable(
    circ: &mut Circuit,
    phase1: &QReg,
    phase2: &QReg,
    sign: &QReg,
    output_less_than: &QReg,
    enable: &QReg,
    chain: &[QReg],
) {
    // On the promised coefficient transition, the input branch bit is
    // phase1 AND (phase2 OR sign_out OR (output < t)).
    circ.ccx(phase1, phase2, enable);
    circ.x(phase2);
    multi_controlled_x_vchain(circ, &[phase1, phase2, sign], enable, chain);
    circ.x(sign);
    multi_controlled_x_vchain(
        circ,
        &[phase1, phase2, sign, output_less_than],
        enable,
        chain,
    );
    circ.x(sign);
    circ.x(phase2);
}

fn toggle_coefficient_target_length(
    circ: &mut Circuit,
    control: &QReg,
    work2: &[QReg],
    l_s: &[QReg],
    l_r_prime: &[QReg],
    target_length: &[QReg],
) {
    assert_eq!(target_length.len(), l_s.len());
    assert_eq!(target_length.len(), l_r_prime.len());
    let length_width = target_length.len();
    let boundary = circ.alloc_qreg_bits("rs.coeff-compare.boundary", length_width);
    for (source, destination) in l_r_prime.iter().zip(&boundary) {
        circ.cx(source, destination);
    }
    let boundary_carry = circ.alloc_qreg("rs.coeff-compare.boundary-carry");
    let boundary_overflow = circ.alloc_qreg("rs.coeff-compare.boundary-overflow");
    cuccaro_add_mod_2n(circ, l_s, &boundary, &boundary_carry, &boundary_overflow);

    controlled_xor_rotated_prefix_bit_length(circ, control, &boundary, work2, target_length);

    cuccaro_sub_mod_2n(circ, l_s, &boundary, &boundary_carry, &boundary_overflow);
    circ.zero_and_free(boundary_overflow);
    circ.zero_and_free(boundary_carry);
    for (source, destination) in l_r_prime.iter().zip(&boundary) {
        circ.cx(source, destination);
    }
    free_clean(circ, boundary);
}

fn toggle_coefficient_length_above_boundary(
    circ: &mut Circuit,
    target_length: &[QReg],
    l_t: &[QReg],
    l_s: &[QReg],
    target: &QReg,
) {
    assert_eq!(target_length.len(), l_t.len());
    assert_eq!(target_length.len(), l_s.len());
    let length_width = l_t.len();
    let wide = length_width + 1;
    let difference = circ.alloc_qreg_bits("rs.coeff-compare.difference", wide);
    let t_copy = circ.alloc_qreg_bits("rs.coeff-compare.t-copy", wide);
    for (source, destination) in target_length.iter().zip(&difference) {
        circ.cx(source, destination);
    }
    for (source, destination) in l_t.iter().zip(&t_copy) {
        circ.cx(source, destination);
    }
    let boundary_carry = circ.alloc_qreg("rs.coeff-compare.boundary-carry");
    let boundary_overflow = circ.alloc_qreg("rs.coeff-compare.boundary-overflow");
    cuccaro_add_mod_2n(
        circ,
        l_s,
        &t_copy[..length_width],
        &boundary_carry,
        &boundary_overflow,
    );
    let compare_carry = circ.alloc_qreg("rs.coeff-compare.carry");
    let compare_overflow = circ.alloc_qreg("rs.coeff-compare.overflow");
    cuccaro_sub_mod_2n(
        circ,
        &t_copy,
        &difference,
        &compare_carry,
        &compare_overflow,
    );
    let nonzero = circ.alloc_qreg("rs.coeff-compare.nonzero");
    let nonzero_chain =
        circ.alloc_qreg_bits("rs.coeff-compare.nonzero-chain", wide.saturating_sub(2));
    compute_nonzero(circ, &difference, &nonzero, &nonzero_chain);
    let sign = &difference[wide - 1];
    circ.x(sign);
    circ.ccx(&nonzero, sign, target);
    circ.x(sign);
    uncompute_nonzero(circ, &difference, &nonzero, &nonzero_chain);
    free_clean(circ, nonzero_chain);
    circ.zero_and_free(nonzero);
    cuccaro_add_mod_2n(
        circ,
        &t_copy,
        &difference,
        &compare_carry,
        &compare_overflow,
    );
    circ.zero_and_free(compare_overflow);
    circ.zero_and_free(compare_carry);
    cuccaro_sub_mod_2n(
        circ,
        l_s,
        &t_copy[..length_width],
        &boundary_carry,
        &boundary_overflow,
    );
    circ.zero_and_free(boundary_overflow);
    circ.zero_and_free(boundary_carry);
    for (source, destination) in l_t.iter().zip(&t_copy) {
        circ.cx(source, destination);
    }
    for lane in t_copy {
        circ.zero_and_free(lane);
    }
    for (source, destination) in target_length.iter().zip(&difference) {
        circ.cx(source, destination);
    }
    for lane in difference {
        circ.zero_and_free(lane);
    }
}

#[allow(clippy::too_many_arguments)]
fn toggle_coefficient_less_than(
    circ: &mut Circuit,
    control: &QReg,
    work1: &[QReg],
    work2: &[QReg],
    l_t: &[QReg],
    l_s: &[QReg],
    target_length: &[QReg],
    target: &QReg,
) {
    let above_t = circ.alloc_qreg("rs.coeff-compare.above-t");
    toggle_coefficient_length_above_boundary(circ, target_length, l_t, l_s, &above_t);

    let cursor = circ.alloc_qreg_bits("rs.coeff-compare.cursor", l_t.len());
    for (source, destination) in l_t.iter().zip(&cursor) {
        circ.cx(source, destination);
    }
    let cursor_scratch = circ.alloc_qreg_bits(
        "rs.coeff-compare.cursor-scratch",
        l_t.len().saturating_sub(1),
    );
    decrement_mod_2n(circ, &cursor, &cursor_scratch);
    let carry = circ.alloc_qreg("rs.coeff-compare.local-carry");
    let active = circ.alloc_qreg("rs.coeff-compare.active");
    let tmp = circ.alloc_qreg("rs.coeff-compare.tmp");

    for index in 0..work1.len() {
        toggle_control_and_nonnegative(circ, control, &cursor, &active);
        circ.ccx(&active, &work1[index], &work2[index]);
        circ.ccx(&active, &carry, &work1[index]);
        multi_controlled_x_vchain(
            circ,
            &[&active, &work1[index], &work2[index]],
            &carry,
            std::slice::from_ref(&tmp),
        );
        toggle_control_and_nonnegative(circ, control, &cursor, &active);
        if index + 1 != work1.len() {
            decrement_mod_2n(circ, &cursor, &cursor_scratch);
        }
    }

    circ.x(&above_t);
    circ.ccx(&carry, &above_t, target);
    circ.x(&above_t);

    for index in (0..work1.len()).rev() {
        if index + 1 != work1.len() {
            increment_mod_2n(circ, &cursor, &cursor_scratch);
        }
        toggle_control_and_nonnegative(circ, control, &cursor, &active);
        multi_controlled_x_vchain(
            circ,
            &[&active, &work1[index], &work2[index]],
            &carry,
            std::slice::from_ref(&tmp),
        );
        circ.ccx(&active, &carry, &work1[index]);
        circ.ccx(&active, &work1[index], &work2[index]);
        toggle_control_and_nonnegative(circ, control, &cursor, &active);
    }
    increment_mod_2n(circ, &cursor, &cursor_scratch);
    circ.zero_and_free(tmp);
    circ.zero_and_free(active);
    circ.zero_and_free(carry);
    free_clean(circ, cursor_scratch);
    for (source, destination) in l_t.iter().zip(&cursor) {
        circ.cx(source, destination);
    }
    free_clean(circ, cursor);

    toggle_coefficient_length_above_boundary(circ, target_length, l_t, l_s, &above_t);
    circ.zero_and_free(above_t);
}

fn coefficient_add_data_only(
    circ: &mut Circuit,
    control: &QReg,
    work1: &[QReg],
    work2: &[QReg],
    l_t: &[QReg],
    inverse: bool,
) {
    let cursor = circ.alloc_qreg_bits("rs.coeff-add.cursor", l_t.len());
    for (source, destination) in l_t.iter().zip(&cursor) {
        circ.cx(source, destination);
    }
    let cursor_scratch =
        circ.alloc_qreg_bits("rs.coeff-add.cursor-scratch", l_t.len().saturating_sub(1));
    let carry = circ.alloc_qreg("rs.coeff-add.carry");
    let active = circ.alloc_qreg("rs.coeff-add.active");
    let tmp = circ.alloc_qreg("rs.coeff-add.tmp");

    if inverse {
        for index in 0..work1.len() {
            if index != 0 {
                decrement_mod_2n(circ, &cursor, &cursor_scratch);
            }
            toggle_control_and_nonnegative(circ, control, &cursor, &active);
            circ.ccx(&active, &work1[index], &work2[index]);
            circ.ccx(&active, &carry, &work1[index]);
            multi_controlled_x_vchain(
                circ,
                &[&active, &work1[index], &work2[index]],
                &carry,
                std::slice::from_ref(&tmp),
            );
            toggle_control_and_nonnegative(circ, control, &cursor, &active);
        }
        for index in (0..work1.len()).rev() {
            if index + 1 != work1.len() {
                increment_mod_2n(circ, &cursor, &cursor_scratch);
            }
            toggle_control_and_nonnegative(circ, control, &cursor, &active);
            multi_controlled_x_vchain(
                circ,
                &[&active, &work1[index], &work2[index]],
                &carry,
                std::slice::from_ref(&tmp),
            );
            circ.ccx(&active, &carry, &work1[index]);
            circ.ccx(&active, &carry, &work2[index]);
            toggle_control_and_nonnegative(circ, control, &cursor, &active);
        }
    } else {
        for index in 0..work1.len() {
            toggle_control_and_nonnegative(circ, control, &cursor, &active);
            circ.ccx(&active, &carry, &work2[index]);
            circ.ccx(&active, &carry, &work1[index]);
            multi_controlled_x_vchain(
                circ,
                &[&active, &work1[index], &work2[index]],
                &carry,
                std::slice::from_ref(&tmp),
            );
            toggle_control_and_nonnegative(circ, control, &cursor, &active);
            if index + 1 != work1.len() {
                decrement_mod_2n(circ, &cursor, &cursor_scratch);
            }
        }
        for index in (0..work1.len()).rev() {
            toggle_control_and_nonnegative(circ, control, &cursor, &active);
            multi_controlled_x_vchain(
                circ,
                &[&active, &work1[index], &work2[index]],
                &carry,
                std::slice::from_ref(&tmp),
            );
            circ.ccx(&active, &carry, &work1[index]);
            circ.ccx(&active, &work1[index], &work2[index]);
            toggle_control_and_nonnegative(circ, control, &cursor, &active);
            if index != 0 {
                increment_mod_2n(circ, &cursor, &cursor_scratch);
            }
        }
    }

    circ.zero_and_free(tmp);
    circ.zero_and_free(active);
    circ.zero_and_free(carry);
    free_clean(circ, cursor_scratch);
    for (source, destination) in l_t.iter().zip(&cursor) {
        circ.cx(source, destination);
    }
    free_clean(circ, cursor);
}

#[allow(clippy::too_many_arguments)]
fn coefficient_phase_block(
    circ: &mut Circuit,
    phase1: &QReg,
    phase2: &QReg,
    sign: &QReg,
    work1: &[QReg],
    work2: &[QReg],
    full_work2: &[QReg],
    l_t: &[QReg],
    l_t_prime: &[QReg],
    l_s: &[QReg],
    l_r_prime: &[QReg],
    inverse: bool,
) {
    let enable = circ.alloc_qreg("rs.coeff-block.enable");
    let less_than = circ.alloc_qreg("rs.coeff-block.less-than");
    let add_only = circ.alloc_qreg("rs.coeff-block.add-only");
    let chain = circ.alloc_qreg_bits("rs.coeff-block.chain", 2);
    assert_eq!(l_t.len(), l_t_prime.len());

    if inverse {
        toggle_coefficient_less_than(circ, phase1, work1, work2, l_t, l_s, l_t_prime, &less_than);
        toggle_output_coefficient_enable(circ, phase1, phase2, sign, &less_than, &enable, &chain);

        circ.x(&enable);
        circ.ccx(phase1, &enable, &add_only);
        circ.x(&enable);
        controlled_xor_raw_t_prime_bit_length(
            circ, &add_only, l_s, l_r_prime, full_work2, l_t_prime,
        );
        coefficient_add_data_only(circ, &add_only, work1, work2, l_t, true);
        controlled_xor_raw_t_prime_bit_length(
            circ, &add_only, l_s, l_r_prime, full_work2, l_t_prime,
        );
        circ.x(&enable);
        circ.ccx(phase1, &enable, &add_only);
        circ.x(&enable);

        circ.cx(&less_than, sign);
        circ.cx(phase1, sign);
        toggle_coefficient_less_than(circ, &enable, work1, work2, l_t, l_s, l_t_prime, &less_than);
        toggle_initial_coefficient_enable(circ, phase1, phase2, sign, &enable, &chain);
    } else {
        toggle_initial_coefficient_enable(circ, phase1, phase2, sign, &enable, &chain);
        toggle_coefficient_less_than(circ, &enable, work1, work2, l_t, l_s, l_t_prime, &less_than);
        circ.cx(phase1, sign);
        circ.cx(&less_than, sign);

        circ.x(&enable);
        circ.ccx(phase1, &enable, &add_only);
        circ.x(&enable);
        controlled_xor_raw_t_prime_bit_length(
            circ, &add_only, l_s, l_r_prime, full_work2, l_t_prime,
        );
        coefficient_add_data_only(circ, &add_only, work1, work2, l_t, false);
        controlled_xor_raw_t_prime_bit_length(
            circ, &add_only, l_s, l_r_prime, full_work2, l_t_prime,
        );
        circ.x(&enable);
        circ.ccx(phase1, &enable, &add_only);
        circ.x(&enable);

        toggle_output_coefficient_enable(circ, phase1, phase2, sign, &less_than, &enable, &chain);
        toggle_coefficient_less_than(circ, phase1, work1, work2, l_t, l_s, l_t_prime, &less_than);
    }

    free_clean(circ, chain);
    circ.zero_and_free(add_only);
    circ.zero_and_free(less_than);
    circ.zero_and_free(enable);
}

#[allow(clippy::too_many_arguments)]
pub fn register_shared_full_window_step(
    circ: &mut Circuit,
    phase1: &QReg,
    phase2: &QReg,
    iteration_parity: &QReg,
    sign: &QReg,
    work1: &[QReg],
    work2: &[QReg],
    l_t: &[QReg],
    l_t_prime: &[QReg],
    l_q: &[QReg],
    l_s: &[QReg],
    l_r_prime: &[QReg],
    emit_length_swap: bool,
) {
    assert_eq!(work1.len(), work2.len());
    assert_eq!(l_t.len(), l_t_prime.len());
    assert_eq!(l_t.len(), l_q.len());
    assert_eq!(l_t.len(), l_s.len());
    assert_eq!(l_t.len(), l_r_prime.len());
    let work_width = work1.len();
    let length_width = l_t.len();

    let pre_scratch = circ.alloc_qreg_bits("rs.step.pre-scratch", length_width + 4);
    super::register_shared_eea_microkernels::pre_shift(
        circ,
        phase1,
        phase2,
        work2,
        l_s,
        &pre_scratch,
    );
    free_clean(circ, pre_scratch);

    let remainder_scratch = circ.alloc_qreg_bits(
        "rs.step.remainder-scratch",
        remainder_scratch_width(length_width, length_width),
    );
    remainder_sub_window(
        circ,
        work_width,
        work_width,
        phase1,
        sign,
        work1,
        work2,
        l_t,
        l_q,
        l_s,
        l_r_prime,
        &remainder_scratch,
    );
    remainder_phase_sign_flip(circ, phase1, phase2, sign, l_r_prime, &remainder_scratch);
    remainder_add_window(
        circ,
        work_width,
        work_width,
        phase1,
        phase2,
        sign,
        work1,
        work2,
        l_t,
        l_q,
        l_s,
        l_r_prime,
        &remainder_scratch,
    );
    free_clean(circ, remainder_scratch);

    let location_scratch = circ.alloc_qreg_bits("rs.step.location-scratch", length_width + 2);
    location_controlled_swap_one_hot(
        circ,
        phase1,
        phase2,
        sign,
        work1,
        0,
        l_t,
        l_q,
        &location_scratch,
    );
    free_clean(circ, location_scratch);

    coefficient_phase_block(
        circ, phase1, phase2, sign, work1, work2, work2, l_t, l_t_prime, l_s, l_r_prime, false,
    );

    let post_scratch = circ.alloc_qreg_bits("rs.step.post-scratch", length_width + 4);
    super::register_shared_eea_microkernels::post_shift(
        circ,
        phase1,
        phase2,
        work2,
        l_s,
        &post_scratch,
    );
    free_clean(circ, post_scratch);

    let phase_scratch = circ.alloc_qreg_bits(
        "rs.step.phase-scratch",
        normalized_phase_scratch_width(length_width),
    );
    normalized_phase_update(
        circ,
        phase1,
        phase2,
        sign,
        l_q,
        l_r_prime,
        l_s,
        &phase_scratch,
    );
    free_clean(circ, phase_scratch);

    if emit_length_swap {
        let condition_scratch = circ.alloc_qreg_bits("rs.step.swap-condition", length_width + 1);
        let zero_q = &condition_scratch[0];
        let zero_s = &condition_scratch[1];
        let control = &condition_scratch[2];
        let chain = &condition_scratch[3..];
        compute_zero(circ, l_q, zero_q, chain);
        compute_zero(circ, l_s, zero_s, chain);
        circ.ccx(zero_q, zero_s, control);
        conditional_work_and_length_swap(
            circ,
            control,
            iteration_parity,
            work1,
            work2,
            l_t,
            l_t_prime,
            l_q,
            l_s,
            l_r_prime,
            (1, work_width),
            (1, work_width),
            &[],
        );
        circ.ccx(zero_q, zero_s, control);
        uncompute_zero(circ, l_s, zero_s, chain);
        uncompute_zero(circ, l_q, zero_q, chain);
        free_clean(circ, condition_scratch);
    }
}

#[allow(clippy::too_many_arguments)]
pub fn register_shared_scheduled_step(
    circ: &mut Circuit,
    step: usize,
    phase1: &QReg,
    phase2: &QReg,
    iteration_parity: &QReg,
    sign: &QReg,
    work1: &[QReg],
    work2: &[QReg],
    l_t: &[QReg],
    l_t_prime: &[QReg],
    l_q: &[QReg],
    l_s: &[QReg],
    l_r_prime: &[QReg],
) {
    assert_eq!(work1.len(), 259);
    assert_eq!(work2.len(), 259);
    assert_eq!(l_t.len(), REFERENCE_LENGTH_WIDTH);
    assert_eq!(l_t_prime.len(), REFERENCE_LENGTH_WIDTH);
    assert_eq!(l_q.len(), REFERENCE_LENGTH_WIDTH);
    assert_eq!(l_s.len(), REFERENCE_LENGTH_WIDTH);
    assert_eq!(l_r_prime.len(), REFERENCE_LENGTH_WIDTH);
    let windows = reference_active_windows(256, step);

    let pre_scratch = circ.alloc_qreg_bits("rs.scheduled.pre-scratch", REFERENCE_LENGTH_WIDTH + 4);
    super::register_shared_eea_microkernels::pre_shift(
        circ,
        phase1,
        phase2,
        work2,
        l_s,
        &pre_scratch,
    );
    free_clean(circ, pre_scratch);

    let r_start = windows.r_add_sub.0 - 1;
    let r_end = windows.r_add_sub.1;
    let remainder_scratch = circ.alloc_qreg_bits(
        "rs.scheduled.remainder-scratch",
        remainder_scratch_width(REFERENCE_LENGTH_WIDTH, REFERENCE_LENGTH_WIDTH),
    );
    remainder_sub_window(
        circ,
        259,
        windows.r_add_sub.1,
        phase1,
        sign,
        &work1[r_start..r_end],
        &work2[r_start..r_end],
        l_t,
        l_q,
        l_s,
        l_r_prime,
        &remainder_scratch,
    );
    remainder_phase_sign_flip(circ, phase1, phase2, sign, l_r_prime, &remainder_scratch);
    remainder_add_window(
        circ,
        259,
        windows.r_add_sub.1,
        phase1,
        phase2,
        sign,
        &work1[r_start..r_end],
        &work2[r_start..r_end],
        l_t,
        l_q,
        l_s,
        l_r_prime,
        &remainder_scratch,
    );
    free_clean(circ, remainder_scratch);

    let swap_start = windows.quotient_swap.0 - 1;
    let swap_end = windows.quotient_swap.1;
    let location_scratch =
        circ.alloc_qreg_bits("rs.scheduled.location-scratch", REFERENCE_LENGTH_WIDTH + 2);
    location_controlled_swap_one_hot(
        circ,
        phase1,
        phase2,
        sign,
        &work1[swap_start..swap_end],
        swap_start,
        l_t,
        l_q,
        &location_scratch,
    );
    free_clean(circ, location_scratch);

    let t_end = windows.t_add_sub.1;
    coefficient_phase_block(
        circ,
        phase1,
        phase2,
        sign,
        &work1[..t_end],
        &work2[..t_end],
        work2,
        l_t,
        l_t_prime,
        l_s,
        l_r_prime,
        false,
    );

    let post_scratch =
        circ.alloc_qreg_bits("rs.scheduled.post-scratch", REFERENCE_LENGTH_WIDTH + 4);
    super::register_shared_eea_microkernels::post_shift(
        circ,
        phase1,
        phase2,
        work2,
        l_s,
        &post_scratch,
    );
    free_clean(circ, post_scratch);

    let phase_scratch = circ.alloc_qreg_bits(
        "rs.scheduled.phase-scratch",
        normalized_phase_scratch_width(REFERENCE_LENGTH_WIDTH),
    );
    normalized_phase_update(
        circ,
        phase1,
        phase2,
        sign,
        l_q,
        l_r_prime,
        l_s,
        &phase_scratch,
    );
    free_clean(circ, phase_scratch);

    if step % 4 == 0 {
        let condition_scratch =
            circ.alloc_qreg_bits("rs.scheduled.swap-condition", REFERENCE_LENGTH_WIDTH + 1);
        let zero_q = &condition_scratch[0];
        let zero_s = &condition_scratch[1];
        let control = &condition_scratch[2];
        let chain = &condition_scratch[3..];
        compute_zero(circ, l_q, zero_q, chain);
        compute_zero(circ, l_s, zero_s, chain);
        circ.ccx(zero_q, zero_s, control);
        conditional_work_and_length_swap(
            circ,
            control,
            iteration_parity,
            work1,
            work2,
            l_t,
            l_t_prime,
            l_q,
            l_s,
            l_r_prime,
            windows.length_update_t,
            windows.length_update_r,
            &[],
        );
        circ.ccx(zero_q, zero_s, control);
        uncompute_zero(circ, l_s, zero_s, chain);
        uncompute_zero(circ, l_q, zero_q, chain);
        free_clean(circ, condition_scratch);
    }
}

#[allow(clippy::too_many_arguments)]
pub fn register_shared_scheduled_step_inverse(
    circ: &mut Circuit,
    step: usize,
    phase1: &QReg,
    phase2: &QReg,
    iteration_parity: &QReg,
    sign: &QReg,
    work1: &[QReg],
    work2: &[QReg],
    l_t: &[QReg],
    l_t_prime: &[QReg],
    l_q: &[QReg],
    l_s: &[QReg],
    l_r_prime: &[QReg],
) {
    assert_eq!(work1.len(), 259);
    assert_eq!(work2.len(), 259);
    assert_eq!(l_t.len(), REFERENCE_LENGTH_WIDTH);
    assert_eq!(l_t_prime.len(), REFERENCE_LENGTH_WIDTH);
    assert_eq!(l_q.len(), REFERENCE_LENGTH_WIDTH);
    assert_eq!(l_s.len(), REFERENCE_LENGTH_WIDTH);
    assert_eq!(l_r_prime.len(), REFERENCE_LENGTH_WIDTH);
    let windows = reference_active_windows(256, step);

    if step % 4 == 0 {
        let condition_scratch = circ.alloc_qreg_bits(
            "rs.scheduled-inverse.swap-condition",
            REFERENCE_LENGTH_WIDTH + 1,
        );
        let zero_q = &condition_scratch[0];
        let zero_s = &condition_scratch[1];
        let control = &condition_scratch[2];
        let chain = &condition_scratch[3..];
        compute_zero(circ, l_q, zero_q, chain);
        compute_zero(circ, l_s, zero_s, chain);
        circ.ccx(zero_q, zero_s, control);
        conditional_work_and_length_swap_inverse(
            circ,
            control,
            iteration_parity,
            work1,
            work2,
            l_t,
            l_t_prime,
            l_q,
            l_s,
            l_r_prime,
            windows.length_update_t,
            windows.length_update_r,
            &[],
        );
        circ.ccx(zero_q, zero_s, control);
        uncompute_zero(circ, l_s, zero_s, chain);
        uncompute_zero(circ, l_q, zero_q, chain);
        free_clean(circ, condition_scratch);
    }

    let phase_scratch = circ.alloc_qreg_bits(
        "rs.scheduled-inverse.phase-scratch",
        normalized_phase_scratch_width(REFERENCE_LENGTH_WIDTH),
    );
    normalized_phase_update_inverse(
        circ,
        phase1,
        phase2,
        sign,
        l_q,
        l_r_prime,
        l_s,
        &phase_scratch,
    );
    free_clean(circ, phase_scratch);

    let post_scratch = circ.alloc_qreg_bits(
        "rs.scheduled-inverse.post-scratch",
        REFERENCE_LENGTH_WIDTH + 4,
    );
    super::register_shared_eea_microkernels::post_shift_inverse(
        circ,
        phase1,
        phase2,
        work2,
        l_s,
        &post_scratch,
    );
    free_clean(circ, post_scratch);

    let t_end = windows.t_add_sub.1;
    coefficient_phase_block(
        circ,
        phase1,
        phase2,
        sign,
        &work1[..t_end],
        &work2[..t_end],
        work2,
        l_t,
        l_t_prime,
        l_s,
        l_r_prime,
        true,
    );

    let swap_start = windows.quotient_swap.0 - 1;
    let swap_end = windows.quotient_swap.1;
    let location_scratch = circ.alloc_qreg_bits(
        "rs.scheduled-inverse.location-scratch",
        REFERENCE_LENGTH_WIDTH + 2,
    );
    location_controlled_swap_one_hot_inverse(
        circ,
        phase1,
        phase2,
        sign,
        &work1[swap_start..swap_end],
        swap_start,
        l_t,
        l_q,
        &location_scratch,
    );
    free_clean(circ, location_scratch);

    let r_start = windows.r_add_sub.0 - 1;
    let r_end = windows.r_add_sub.1;
    let remainder_scratch = circ.alloc_qreg_bits(
        "rs.scheduled-inverse.remainder-scratch",
        remainder_scratch_width(REFERENCE_LENGTH_WIDTH, REFERENCE_LENGTH_WIDTH),
    );
    remainder_add_window_inverse(
        circ,
        259,
        windows.r_add_sub.1,
        phase1,
        phase2,
        sign,
        &work1[r_start..r_end],
        &work2[r_start..r_end],
        l_t,
        l_q,
        l_s,
        l_r_prime,
        &remainder_scratch,
    );
    remainder_phase_sign_flip(circ, phase1, phase2, sign, l_r_prime, &remainder_scratch);
    remainder_sub_window_inverse(
        circ,
        259,
        windows.r_add_sub.1,
        phase1,
        sign,
        &work1[r_start..r_end],
        &work2[r_start..r_end],
        l_t,
        l_q,
        l_s,
        l_r_prime,
        &remainder_scratch,
    );
    free_clean(circ, remainder_scratch);

    let pre_scratch = circ.alloc_qreg_bits(
        "rs.scheduled-inverse.pre-scratch",
        REFERENCE_LENGTH_WIDTH + 4,
    );
    super::register_shared_eea_microkernels::pre_shift_inverse(
        circ,
        phase1,
        phase2,
        work2,
        l_s,
        &pre_scratch,
    );
    free_clean(circ, pre_scratch);
}

struct RegisterSharedCore {
    phase1: QReg,
    phase2: QReg,
    iteration_parity: QReg,
    sign: QReg,
    work1: Vec<QReg>,
    work2: Vec<QReg>,
    l_t: Vec<QReg>,
    l_t_prime: Vec<QReg>,
    l_q: Vec<QReg>,
    l_s: Vec<QReg>,
    l_r_prime: Vec<QReg>,
}

struct RegisterSharedTerminal {
    iteration_parity: QReg,
    work2: Vec<QReg>,
    l_t_prime: Vec<QReg>,
    l_s: Vec<QReg>,
}

const REGISTER_SHARED_WORK_WIDTH: usize = 259;
const REGISTER_SHARED_FIELD_WIDTH: usize = 257;
const REGISTER_SHARED_HALF_BYTES: [u8; 33] = [
    0x18, 0xfe, 0xff, 0x7f, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0x7f,
    0x00,
];

fn toggle_initial_work1(circ: &mut Circuit, work1: &[QReg]) {
    use crate::point_add::trailmix_port::mod_arith::SECP256K1_P_LE;

    assert_eq!(work1.len(), REGISTER_SHARED_WORK_WIDTH);
    circ.x(&work1[0]);
    for bit in 0..256 {
        if ((SECP256K1_P_LE[bit / 8] >> (bit % 8)) & 1) != 0 {
            circ.x(&work1[REGISTER_SHARED_WORK_WIDTH - 1 - bit]);
        }
    }
}

fn toggle_terminal_work1(circ: &mut Circuit, work1: &[QReg]) {
    use crate::point_add::trailmix_port::mod_arith::SECP256K1_P_LE;

    assert_eq!(work1.len(), REGISTER_SHARED_WORK_WIDTH);
    for bit in 0..256 {
        if ((SECP256K1_P_LE[bit / 8] >> (bit % 8)) & 1) != 0 {
            circ.x(&work1[bit]);
        }
    }
    circ.x(&work1[REGISTER_SHARED_WORK_WIDTH - 1]);
}

fn register_shared_initialize(circ: &mut Circuit, mut dx: Vec<QReg>) -> RegisterSharedCore {
    use super::shrunken_pz_state_machine::{bit_length_lean, controlled_field_neg};
    use crate::point_add::trailmix_port::arith::compare::compare_geq_const;

    assert_eq!(dx.len(), REGISTER_SHARED_FIELD_WIDTH);
    let iteration_parity = circ.alloc_qreg("rs.divider.iteration-parity");
    compare_geq_const(circ, &dx, &REGISTER_SHARED_HALF_BYTES, &iteration_parity);
    controlled_field_neg(circ, &iteration_parity, &dx);

    let l_r_prime = circ.alloc_qreg_bits("rs.divider.l-r-prime", REFERENCE_LENGTH_WIDTH);
    let source: Vec<&QReg> = dx.iter().take(256).collect();
    bit_length_lean(circ, &source, &l_r_prime, false);

    dx.push(circ.alloc_qreg("rs.divider.work2-pad0"));
    dx.push(circ.alloc_qreg("rs.divider.work2-pad1"));
    dx.reverse();
    let work2 = dx;

    let work1 = circ.alloc_qreg_bits("rs.divider.work1", REGISTER_SHARED_WORK_WIDTH);
    toggle_initial_work1(circ, &work1);
    let l_t = circ.alloc_qreg_bits("rs.divider.l-t", REFERENCE_LENGTH_WIDTH);
    let l_t_prime = circ.alloc_qreg_bits("rs.divider.l-t-prime", REFERENCE_LENGTH_WIDTH);
    let l_q = circ.alloc_qreg_bits("rs.divider.l-q", REFERENCE_LENGTH_WIDTH);
    let l_s = circ.alloc_qreg_bits("rs.divider.l-s", REFERENCE_LENGTH_WIDTH);
    circ.x(&l_t[0]);
    let phase1 = circ.alloc_qreg("rs.divider.phase1");
    let phase2 = circ.alloc_qreg("rs.divider.phase2");
    let sign = circ.alloc_qreg("rs.divider.sign");

    RegisterSharedCore {
        phase1,
        phase2,
        iteration_parity,
        sign,
        work1,
        work2,
        l_t,
        l_t_prime,
        l_q,
        l_s,
        l_r_prime,
    }
}

fn register_shared_forward(circ: &mut Circuit, core: &RegisterSharedCore) {
    for step in 1..=REFERENCE_STEPS {
        register_shared_scheduled_step(
            circ,
            step,
            &core.phase1,
            &core.phase2,
            &core.iteration_parity,
            &core.sign,
            &core.work1,
            &core.work2,
            &core.l_t,
            &core.l_t_prime,
            &core.l_q,
            &core.l_s,
            &core.l_r_prime,
        );
    }
}

fn register_shared_release_terminal(
    circ: &mut Circuit,
    core: RegisterSharedCore,
) -> RegisterSharedTerminal {
    toggle_terminal_work1(circ, &core.work1);
    free_clean(circ, core.work1);
    toggle_constant(circ, &core.l_t, 256);
    free_clean(circ, core.l_t);
    free_clean(circ, core.l_q);
    free_clean(circ, core.l_r_prime);
    circ.zero_and_free(core.phase1);
    circ.zero_and_free(core.phase2);
    circ.zero_and_free(core.sign);
    RegisterSharedTerminal {
        iteration_parity: core.iteration_parity,
        work2: core.work2,
        l_t_prime: core.l_t_prime,
        l_s: core.l_s,
    }
}

fn register_shared_rebuild_terminal(
    circ: &mut Circuit,
    terminal: RegisterSharedTerminal,
) -> RegisterSharedCore {
    let work1 = circ.alloc_qreg_bits("rs.divider.work1.rebuilt", REGISTER_SHARED_WORK_WIDTH);
    toggle_terminal_work1(circ, &work1);
    let l_t = circ.alloc_qreg_bits("rs.divider.l-t.rebuilt", REFERENCE_LENGTH_WIDTH);
    let l_q = circ.alloc_qreg_bits("rs.divider.l-q.rebuilt", REFERENCE_LENGTH_WIDTH);
    let l_r_prime = circ.alloc_qreg_bits("rs.divider.l-r-prime.rebuilt", REFERENCE_LENGTH_WIDTH);
    toggle_constant(circ, &l_t, 256);
    RegisterSharedCore {
        phase1: circ.alloc_qreg("rs.divider.phase1.rebuilt"),
        phase2: circ.alloc_qreg("rs.divider.phase2.rebuilt"),
        iteration_parity: terminal.iteration_parity,
        sign: circ.alloc_qreg("rs.divider.sign.rebuilt"),
        work1,
        work2: terminal.work2,
        l_t,
        l_t_prime: terminal.l_t_prime,
        l_q,
        l_s: terminal.l_s,
        l_r_prime,
    }
}

fn register_shared_reverse(circ: &mut Circuit, core: &RegisterSharedCore) {
    for step in (1..=REFERENCE_STEPS).rev() {
        register_shared_scheduled_step_inverse(
            circ,
            step,
            &core.phase1,
            &core.phase2,
            &core.iteration_parity,
            &core.sign,
            &core.work1,
            &core.work2,
            &core.l_t,
            &core.l_t_prime,
            &core.l_q,
            &core.l_s,
            &core.l_r_prime,
        );
    }
}

fn register_shared_finish(circ: &mut Circuit, mut core: RegisterSharedCore) -> Vec<QReg> {
    use super::shrunken_pz_state_machine::{bit_length_lean, controlled_field_neg};
    use crate::point_add::trailmix_port::arith::compare::compare_geq_const;

    circ.zero_and_free(core.phase1);
    circ.zero_and_free(core.phase2);
    circ.zero_and_free(core.sign);
    circ.x(&core.l_t[0]);
    free_clean(circ, core.l_t);
    free_clean(circ, core.l_t_prime);
    free_clean(circ, core.l_q);
    free_clean(circ, core.l_s);
    toggle_initial_work1(circ, &core.work1);
    free_clean(circ, core.work1);

    core.work2.reverse();
    let pad1 = core.work2.pop().expect("register-shared Work2 pad1");
    let pad0 = core.work2.pop().expect("register-shared Work2 pad0");
    circ.zero_and_free(pad1);
    circ.zero_and_free(pad0);
    assert_eq!(core.work2.len(), REGISTER_SHARED_FIELD_WIDTH);

    let source: Vec<&QReg> = core.work2.iter().take(256).collect();
    bit_length_lean(circ, &source, &core.l_r_prime, true);
    free_clean(circ, core.l_r_prime);
    controlled_field_neg(circ, &core.iteration_parity, &core.work2);
    compare_geq_const(
        circ,
        &core.work2,
        &REGISTER_SHARED_HALF_BYTES,
        &core.iteration_parity,
    );
    circ.zero_and_free(core.iteration_parity);
    core.work2
}

fn toggle_terminal_inverse_sign(circ: &mut Circuit, terminal: &RegisterSharedTerminal) {
    use super::shrunken_pz_state_machine::controlled_field_neg;

    circ.x(&terminal.iteration_parity);
    controlled_field_neg(
        circ,
        &terminal.iteration_parity,
        &terminal.work2[..REGISTER_SHARED_FIELD_WIDTH],
    );
    circ.x(&terminal.iteration_parity);
}

/// Experimental complete register-shared divider lifecycle.
///
/// The Q883 candidate selects this API in source. It remains unsubmitted until
/// source-bound support and trusted evaluator gates are complete.
pub fn register_shared_divide_forward(
    circ: &mut Circuit,
    dx: Vec<QReg>,
    dy: Vec<QReg>,
) -> (Vec<QReg>, Vec<QReg>, Vec<QReg>) {
    use crate::point_add::trailmix_port::arith::rfold_mbu::mod_mul_canonical_mbu;

    assert_eq!(dx.len(), REGISTER_SHARED_FIELD_WIDTH);
    assert_eq!(dy.len(), REGISTER_SHARED_FIELD_WIDTH);
    let core = register_shared_initialize(circ, dx);
    register_shared_forward(circ, &core);
    let terminal = register_shared_release_terminal(circ, core);

    variable_rotate_high(circ, &terminal.l_s, &terminal.work2);
    toggle_terminal_inverse_sign(circ, &terminal);
    let mut lambda = circ.alloc_qreg_bits("rs.divider.lambda", REGISTER_SHARED_FIELD_WIDTH);
    mod_mul_canonical_mbu(
        circ,
        &lambda,
        &terminal.work2[..REGISTER_SHARED_FIELD_WIDTH],
        &dy,
    );
    toggle_terminal_inverse_sign(circ, &terminal);
    variable_rotate_low(circ, &terminal.l_s, &terminal.work2);

    let lambda_top = lambda.pop().expect("canonical lambda top lane");
    circ.zero_and_free(lambda_top);
    let dy_ghosts: Vec<_> = dy.iter().map(|lane| circ.hmr_ghost(lane)).collect();
    free_clean(circ, dy);

    let core = register_shared_rebuild_terminal(circ, terminal);
    register_shared_reverse(circ, &core);
    let dx = register_shared_finish(circ, core);

    lambda.push(circ.alloc_qreg("rs.divider.lambda-top-restored"));
    let dy = circ.alloc_qreg_bits("rs.divider.dy-restored", REGISTER_SHARED_FIELD_WIDTH);
    mod_mul_canonical_mbu(circ, &dy, &lambda, &dx);
    for (ghost, lane) in dy_ghosts.into_iter().zip(&dy) {
        circ.resolve_ghost(ghost, lane);
    }
    (dx, dy, lambda)
}

/// Experimental inverse-witness cleanup for [`register_shared_divide_forward`].
pub fn register_shared_divide_cancel(
    circ: &mut Circuit,
    dx: Vec<QReg>,
    dy: Vec<QReg>,
    lambda: Vec<QReg>,
) -> (Vec<QReg>, Vec<QReg>) {
    use crate::point_add::trailmix_port::arith::rfold_mbu::{
        mod_mul_canonical_mbu, mod_mul_canonical_mbu_undo,
    };

    assert_eq!(dx.len(), REGISTER_SHARED_FIELD_WIDTH);
    assert_eq!(dy.len(), REGISTER_SHARED_FIELD_WIDTH);
    assert_eq!(lambda.len(), REGISTER_SHARED_FIELD_WIDTH);
    let lambda_ghosts: Vec<_> = lambda.iter().map(|lane| circ.hmr_ghost(lane)).collect();
    free_clean(circ, lambda);

    let core = register_shared_initialize(circ, dx);
    register_shared_forward(circ, &core);
    let terminal = register_shared_release_terminal(circ, core);
    variable_rotate_high(circ, &terminal.l_s, &terminal.work2);
    toggle_terminal_inverse_sign(circ, &terminal);

    let quotient = circ.alloc_qreg_bits("rs.divider.quotient-check", REGISTER_SHARED_FIELD_WIDTH);
    mod_mul_canonical_mbu(
        circ,
        &quotient,
        &terminal.work2[..REGISTER_SHARED_FIELD_WIDTH],
        &dy,
    );
    for (ghost, lane) in lambda_ghosts.into_iter().zip(&quotient) {
        circ.resolve_ghost(ghost, lane);
    }
    mod_mul_canonical_mbu_undo(
        circ,
        &quotient,
        &terminal.work2[..REGISTER_SHARED_FIELD_WIDTH],
        &dy,
    );
    free_clean(circ, quotient);

    toggle_terminal_inverse_sign(circ, &terminal);
    variable_rotate_low(circ, &terminal.l_s, &terminal.work2);
    let core = register_shared_rebuild_terminal(circ, terminal);
    register_shared_reverse(circ, &core);
    let dx = register_shared_finish(circ, core);
    (dx, dy)
}

#[allow(clippy::too_many_arguments)]
pub fn register_shared_full_window_step_inverse(
    circ: &mut Circuit,
    phase1: &QReg,
    phase2: &QReg,
    iteration_parity: &QReg,
    sign: &QReg,
    work1: &[QReg],
    work2: &[QReg],
    l_t: &[QReg],
    l_t_prime: &[QReg],
    l_q: &[QReg],
    l_s: &[QReg],
    l_r_prime: &[QReg],
    emit_length_swap: bool,
) {
    assert_eq!(work1.len(), work2.len());
    assert_eq!(l_t.len(), l_t_prime.len());
    let work_width = work1.len();
    let length_width = l_t.len();

    if emit_length_swap {
        let condition_scratch = circ.alloc_qreg_bits("rs.step.swap-condition", length_width + 1);
        let zero_q = &condition_scratch[0];
        let zero_s = &condition_scratch[1];
        let control = &condition_scratch[2];
        let chain = &condition_scratch[3..];
        compute_zero(circ, l_q, zero_q, chain);
        compute_zero(circ, l_s, zero_s, chain);
        circ.ccx(zero_q, zero_s, control);
        conditional_work_and_length_swap_inverse(
            circ,
            control,
            iteration_parity,
            work1,
            work2,
            l_t,
            l_t_prime,
            l_q,
            l_s,
            l_r_prime,
            (1, work_width),
            (1, work_width),
            &[],
        );
        circ.ccx(zero_q, zero_s, control);
        uncompute_zero(circ, l_s, zero_s, chain);
        uncompute_zero(circ, l_q, zero_q, chain);
        free_clean(circ, condition_scratch);
    }

    let phase_scratch = circ.alloc_qreg_bits(
        "rs.step.phase-scratch",
        normalized_phase_scratch_width(length_width),
    );
    normalized_phase_update_inverse(
        circ,
        phase1,
        phase2,
        sign,
        l_q,
        l_r_prime,
        l_s,
        &phase_scratch,
    );
    free_clean(circ, phase_scratch);

    let post_scratch = circ.alloc_qreg_bits("rs.step.post-scratch", length_width + 4);
    super::register_shared_eea_microkernels::post_shift_inverse(
        circ,
        phase1,
        phase2,
        work2,
        l_s,
        &post_scratch,
    );
    free_clean(circ, post_scratch);

    coefficient_phase_block(
        circ, phase1, phase2, sign, work1, work2, work2, l_t, l_t_prime, l_s, l_r_prime, true,
    );

    let location_scratch = circ.alloc_qreg_bits("rs.step.location-scratch", length_width + 2);
    location_controlled_swap_one_hot_inverse(
        circ,
        phase1,
        phase2,
        sign,
        work1,
        0,
        l_t,
        l_q,
        &location_scratch,
    );
    free_clean(circ, location_scratch);

    let remainder_scratch = circ.alloc_qreg_bits(
        "rs.step.remainder-scratch",
        remainder_scratch_width(length_width, length_width),
    );
    remainder_add_window_inverse(
        circ,
        work_width,
        work_width,
        phase1,
        phase2,
        sign,
        work1,
        work2,
        l_t,
        l_q,
        l_s,
        l_r_prime,
        &remainder_scratch,
    );
    remainder_phase_sign_flip(circ, phase1, phase2, sign, l_r_prime, &remainder_scratch);
    remainder_sub_window_inverse(
        circ,
        work_width,
        work_width,
        phase1,
        sign,
        work1,
        work2,
        l_t,
        l_q,
        l_s,
        l_r_prime,
        &remainder_scratch,
    );
    free_clean(circ, remainder_scratch);

    let pre_scratch = circ.alloc_qreg_bits("rs.step.pre-scratch", length_width + 4);
    super::register_shared_eea_microkernels::pre_shift_inverse(
        circ,
        phase1,
        phase2,
        work2,
        l_s,
        &pre_scratch,
    );
    free_clean(circ, pre_scratch);
}

#[derive(Clone, Copy)]
enum CoefficientKernel {
    Add,
    AddInverse,
    Sub,
    SubInverse,
}

fn build_coefficient_kernel(
    work_width: usize,
    length_width: usize,
    kernel: CoefficientKernel,
) -> B {
    let mut circ = Circuit::new();
    let phase1 = circ.alloc_qreg("rs.coeff.phase1");
    let phase2 = circ.alloc_qreg("rs.coeff.phase2");
    let sign = circ.alloc_qreg("rs.coeff.sign");
    let work1 = circ.alloc_qreg_bits("rs.coeff.work1", work_width);
    let work2 = circ.alloc_qreg_bits("rs.coeff.work2", work_width);
    let l_t = circ.alloc_qreg_bits("rs.coeff.l-t", length_width);
    let scratch_width = match kernel {
        CoefficientKernel::Add | CoefficientKernel::AddInverse => length_width + 2,
        CoefficientKernel::Sub | CoefficientKernel::SubInverse => length_width + 4,
    };
    let scratch = circ.alloc_qreg_bits("rs.coeff.scratch", scratch_width);
    match kernel {
        CoefficientKernel::Add => {
            coefficient_add_single(&mut circ, &phase1, &sign, &work1, &work2, &l_t, &scratch)
        }
        CoefficientKernel::AddInverse => coefficient_add_single_inverse(
            &mut circ, &phase1, &sign, &work1, &work2, &l_t, &scratch,
        ),
        CoefficientKernel::Sub => coefficient_sub_single(
            &mut circ, &phase1, &phase2, &sign, &work1, &work2, &l_t, &scratch,
        ),
        CoefficientKernel::SubInverse => coefficient_sub_single_inverse(
            &mut circ, &phase1, &phase2, &sign, &work1, &work2, &l_t, &scratch,
        ),
    }
    circ.into_builder()
}

#[must_use]
pub fn exhaustive_reference_coefficient_arithmetic_check(
) -> ReferenceCoefficientArithmeticProofReport {
    const TEST_LENGTH_WIDTH: usize = 3;
    let mut basis_states_checked = 0usize;
    let mut inverse_pair_checks = 0usize;
    let mut scratch_clean_checks = 0usize;
    let mut control_off_identity_checks = 0usize;
    let mut length_restore_checks = 0usize;

    for work_width in 1..=4 {
        let add = build_coefficient_kernel(work_width, TEST_LENGTH_WIDTH, CoefficientKernel::Add);
        let add_inverse =
            build_coefficient_kernel(work_width, TEST_LENGTH_WIDTH, CoefficientKernel::AddInverse);
        let sub = build_coefficient_kernel(work_width, TEST_LENGTH_WIDTH, CoefficientKernel::Sub);
        let sub_inverse =
            build_coefficient_kernel(work_width, TEST_LENGTH_WIDTH, CoefficientKernel::SubInverse);
        let data_width = 3 + 2 * work_width + TEST_LENGTH_WIDTH;
        let length_offset = 3 + 2 * work_width;
        let length_mask = (1u64 << TEST_LENGTH_WIDTH) - 1;
        for input in 0..(1u64 << data_width) {
            for (forward, inverse) in [(&add, &add_inverse), (&sub, &sub_inverse)] {
                let output = apply_scalar(&forward.ops, input);
                assert_eq!(
                    output >> data_width,
                    0,
                    "coefficient kernel left scratch dirty"
                );
                assert_eq!(
                    (output >> length_offset) & length_mask,
                    (input >> length_offset) & length_mask,
                    "coefficient kernel changed l_t"
                );
                if input & 1 == 0 {
                    assert_eq!(output, input, "phase1=0 must disable coefficient kernel");
                    control_off_identity_checks += 1;
                }
                assert_eq!(apply_scalar(&inverse.ops, output), input);
                basis_states_checked += 1;
                inverse_pair_checks += 1;
                scratch_clean_checks += 1;
                length_restore_checks += 1;
            }
        }
    }

    let t_add257 = gate_counts(
        &build_coefficient_kernel(257, REFERENCE_LENGTH_WIDTH, CoefficientKernel::Add).ops,
    );
    let t_sub257 = gate_counts(
        &build_coefficient_kernel(257, REFERENCE_LENGTH_WIDTH, CoefficientKernel::Sub).ops,
    );
    let add1 = gate_counts(
        &build_coefficient_kernel(1, REFERENCE_LENGTH_WIDTH, CoefficientKernel::Add).ops,
    )
    .ccx;
    let add2 = gate_counts(
        &build_coefficient_kernel(2, REFERENCE_LENGTH_WIDTH, CoefficientKernel::Add).ops,
    )
    .ccx;
    let sub1 = gate_counts(
        &build_coefficient_kernel(1, REFERENCE_LENGTH_WIDTH, CoefficientKernel::Sub).ops,
    )
    .ccx;
    let sub2 = gate_counts(
        &build_coefficient_kernel(2, REFERENCE_LENGTH_WIDTH, CoefficientKernel::Sub).ops,
    )
    .ccx;
    let add_slope = add2 - add1;
    let sub_slope = sub2 - sub1;
    for width in 1..=257 {
        let observed_add = gate_counts(
            &build_coefficient_kernel(width, REFERENCE_LENGTH_WIDTH, CoefficientKernel::Add).ops,
        )
        .ccx;
        let observed_sub = gate_counts(
            &build_coefficient_kernel(width, REFERENCE_LENGTH_WIDTH, CoefficientKernel::Sub).ops,
        )
        .ccx;
        assert_eq!(observed_add, add1 + (width - 1) * add_slope);
        assert_eq!(observed_sub, sub1 + (width - 1) * sub_slope);
    }
    let schedule = exhaustive_reference_schedule_check();
    let scheduled_add_toffoli =
        add_slope * schedule.t_window_sum + (add1 - add_slope) * REFERENCE_STEPS;
    let scheduled_sub_toffoli =
        sub_slope * schedule.t_window_sum + (sub1 - sub_slope) * REFERENCE_STEPS;

    ReferenceCoefficientArithmeticProofReport {
        work_widths_checked: 4,
        basis_states_checked,
        inverse_pair_checks,
        scratch_clean_checks,
        control_off_identity_checks,
        length_restore_checks,
        t_add257,
        t_sub257,
        reference_steps: REFERENCE_STEPS,
        scheduled_window_sum: schedule.t_window_sum,
        scheduled_add_toffoli,
        scheduled_sub_toffoli,
    }
}

#[derive(Clone, Copy)]
enum RemainderKernel {
    Add,
    AddInverse,
    Sub,
    SubInverse,
}

fn build_remainder_kernel(
    work_width: usize,
    total_work_width: usize,
    window_upper: usize,
    length_width: usize,
    kernel: RemainderKernel,
) -> B {
    assert!(work_width > 0);
    assert!(work_width <= window_upper && window_upper <= total_work_width);
    let mut circ = Circuit::new();
    let phase1 = circ.alloc_qreg("rs.remainder.phase1");
    let phase2 = circ.alloc_qreg("rs.remainder.phase2");
    let sign = circ.alloc_qreg("rs.remainder.sign");
    let work1 = circ.alloc_qreg_bits("rs.remainder.work1", work_width);
    let work2 = circ.alloc_qreg_bits("rs.remainder.work2", work_width);
    let l_t = circ.alloc_qreg_bits("rs.remainder.l-t", length_width);
    let l_q = circ.alloc_qreg_bits("rs.remainder.l-q", length_width);
    let l_s = circ.alloc_qreg_bits("rs.remainder.l-s", length_width);
    let l_r_prime = circ.alloc_qreg_bits("rs.remainder.l-r-prime", length_width);
    let scratch = circ.alloc_qreg_bits(
        "rs.remainder.scratch",
        remainder_scratch_width(length_width, length_width),
    );
    match kernel {
        RemainderKernel::Add => remainder_add_window(
            &mut circ,
            total_work_width,
            window_upper,
            &phase1,
            &phase2,
            &sign,
            &work1,
            &work2,
            &l_t,
            &l_q,
            &l_s,
            &l_r_prime,
            &scratch,
        ),
        RemainderKernel::AddInverse => remainder_add_window_inverse(
            &mut circ,
            total_work_width,
            window_upper,
            &phase1,
            &phase2,
            &sign,
            &work1,
            &work2,
            &l_t,
            &l_q,
            &l_s,
            &l_r_prime,
            &scratch,
        ),
        RemainderKernel::Sub => remainder_sub_window(
            &mut circ,
            total_work_width,
            window_upper,
            &phase1,
            &sign,
            &work1,
            &work2,
            &l_t,
            &l_q,
            &l_s,
            &l_r_prime,
            &scratch,
        ),
        RemainderKernel::SubInverse => remainder_sub_window_inverse(
            &mut circ,
            total_work_width,
            window_upper,
            &phase1,
            &sign,
            &work1,
            &work2,
            &l_t,
            &l_q,
            &l_s,
            &l_r_prime,
            &scratch,
        ),
    }
    circ.into_builder()
}

#[must_use]
pub fn exhaustive_reference_remainder_arithmetic_check() -> ReferenceRemainderArithmeticProofReport
{
    const TEST_LENGTH_WIDTH: usize = 2;
    let mut basis_states_checked = 0usize;
    let mut inverse_pair_checks = 0usize;
    let mut scratch_clean_checks = 0usize;
    let mut control_off_identity_checks = 0usize;
    let mut zero_remainder_identity_checks = 0usize;
    let mut length_restore_checks = 0usize;

    for work_width in 1..=3 {
        let add = build_remainder_kernel(
            work_width,
            work_width,
            work_width,
            TEST_LENGTH_WIDTH,
            RemainderKernel::Add,
        );
        let add_inverse = build_remainder_kernel(
            work_width,
            work_width,
            work_width,
            TEST_LENGTH_WIDTH,
            RemainderKernel::AddInverse,
        );
        let sub = build_remainder_kernel(
            work_width,
            work_width,
            work_width,
            TEST_LENGTH_WIDTH,
            RemainderKernel::Sub,
        );
        let sub_inverse = build_remainder_kernel(
            work_width,
            work_width,
            work_width,
            TEST_LENGTH_WIDTH,
            RemainderKernel::SubInverse,
        );
        let data_width = 3 + 2 * work_width + 4 * TEST_LENGTH_WIDTH;
        let lengths_offset = 3 + 2 * work_width;
        let l_r_prime_offset = lengths_offset + 3 * TEST_LENGTH_WIDTH;
        let lengths_mask = (1u64 << (4 * TEST_LENGTH_WIDTH)) - 1;
        let one_length_mask = (1u64 << TEST_LENGTH_WIDTH) - 1;
        for input in 0..(1u64 << data_width) {
            for (forward, inverse) in [(&add, &add_inverse), (&sub, &sub_inverse)] {
                let output = apply_scalar(&forward.ops, input);
                assert_eq!(
                    output >> data_width,
                    0,
                    "remainder kernel left scratch dirty"
                );
                assert_eq!(
                    (output >> lengths_offset) & lengths_mask,
                    (input >> lengths_offset) & lengths_mask,
                    "remainder kernel changed a length register"
                );
                if input & 1 != 0 {
                    assert_eq!(output, input, "phase1=1 must disable remainder kernel");
                    control_off_identity_checks += 1;
                }
                if ((input >> l_r_prime_offset) & one_length_mask) == 0 {
                    assert_eq!(output, input, "l_r_prime=0 must disable remainder kernel");
                    zero_remainder_identity_checks += 1;
                }
                assert_eq!(apply_scalar(&inverse.ops, output), input);
                basis_states_checked += 1;
                inverse_pair_checks += 1;
                scratch_clean_checks += 1;
                length_restore_checks += 1;
            }
        }
    }

    let r_add257 = gate_counts(
        &build_remainder_kernel(257, 259, 259, REFERENCE_LENGTH_WIDTH, RemainderKernel::Add).ops,
    );
    let r_sub257 = gate_counts(
        &build_remainder_kernel(257, 259, 259, REFERENCE_LENGTH_WIDTH, RemainderKernel::Sub).ops,
    );
    let add1 = gate_counts(
        &build_remainder_kernel(1, 259, 259, REFERENCE_LENGTH_WIDTH, RemainderKernel::Add).ops,
    )
    .ccx;
    let add2 = gate_counts(
        &build_remainder_kernel(2, 259, 259, REFERENCE_LENGTH_WIDTH, RemainderKernel::Add).ops,
    )
    .ccx;
    let sub1 = gate_counts(
        &build_remainder_kernel(1, 259, 259, REFERENCE_LENGTH_WIDTH, RemainderKernel::Sub).ops,
    )
    .ccx;
    let sub2 = gate_counts(
        &build_remainder_kernel(2, 259, 259, REFERENCE_LENGTH_WIDTH, RemainderKernel::Sub).ops,
    )
    .ccx;
    let add_slope = add2 - add1;
    let sub_slope = sub2 - sub1;
    for width in 1..=257 {
        let observed_add = gate_counts(
            &build_remainder_kernel(
                width,
                259,
                259,
                REFERENCE_LENGTH_WIDTH,
                RemainderKernel::Add,
            )
            .ops,
        )
        .ccx;
        let observed_sub = gate_counts(
            &build_remainder_kernel(
                width,
                259,
                259,
                REFERENCE_LENGTH_WIDTH,
                RemainderKernel::Sub,
            )
            .ops,
        )
        .ccx;
        assert_eq!(observed_add, add1 + (width - 1) * add_slope);
        assert_eq!(observed_sub, sub1 + (width - 1) * sub_slope);
    }
    let schedule = exhaustive_reference_schedule_check();
    let scheduled_add_toffoli =
        add_slope * schedule.r_window_sum + (add1 - add_slope) * REFERENCE_STEPS;
    let scheduled_sub_toffoli =
        sub_slope * schedule.r_window_sum + (sub1 - sub_slope) * REFERENCE_STEPS;

    ReferenceRemainderArithmeticProofReport {
        work_widths_checked: 3,
        basis_states_checked,
        inverse_pair_checks,
        scratch_clean_checks,
        control_off_identity_checks,
        zero_remainder_identity_checks,
        length_restore_checks,
        r_add257,
        r_sub257,
        reference_steps: REFERENCE_STEPS,
        scheduled_window_sum: schedule.r_window_sum,
        scheduled_add_toffoli,
        scheduled_sub_toffoli,
    }
}

fn build_normalized_phase_update(length_width: usize, inverse: bool) -> B {
    let mut circ = Circuit::new();
    let phase1 = circ.alloc_qreg("rs.normalized-phase.phase1");
    let phase2 = circ.alloc_qreg("rs.normalized-phase.phase2");
    let sign = circ.alloc_qreg("rs.normalized-phase.sign");
    let l_q = circ.alloc_qreg_bits("rs.normalized-phase.l-q", length_width);
    let l_r_prime = circ.alloc_qreg_bits("rs.normalized-phase.l-r-prime", length_width);
    let l_s = circ.alloc_qreg_bits("rs.normalized-phase.l-s", length_width);
    let scratch = circ.alloc_qreg_bits(
        "rs.normalized-phase.scratch",
        normalized_phase_scratch_width(length_width),
    );
    if inverse {
        normalized_phase_update_inverse(
            &mut circ, &phase1, &phase2, &sign, &l_q, &l_r_prime, &l_s, &scratch,
        );
    } else {
        normalized_phase_update(
            &mut circ, &phase1, &phase2, &sign, &l_q, &l_r_prime, &l_s, &scratch,
        );
    }
    circ.into_builder()
}

#[must_use]
pub fn exhaustive_normalized_phase_update_check() -> ReferenceNormalizedControlProofReport {
    let mut basis_states_checked = 0usize;
    let mut inverse_pair_checks = 0usize;
    let mut scratch_clean_checks = 0usize;
    let mut oracle_transition_checks = 0usize;
    for length_width in 1..=3 {
        let forward = build_normalized_phase_update(length_width, false);
        let inverse = build_normalized_phase_update(length_width, true);
        let mask = (1u64 << length_width) - 1;
        let data_width = 3 + 3 * length_width;
        for input in 0..(1u64 << data_width) {
            let mut phase1 = input & 1 != 0;
            let mut phase2 = (input >> 1) & 1 != 0;
            let mut sign = (input >> 2) & 1 != 0;
            let l_q = (input >> 3) & mask;
            let l_r_prime = (input >> (3 + length_width)) & mask;
            let l_s = (input >> (3 + 2 * length_width)) & mask;
            if l_q == 0 && l_r_prime != 0 {
                phase2 ^= sign ^ phase1;
                sign ^= phase2;
            }
            if l_s == 0 {
                phase1 ^= true;
                phase2 ^= true;
            }
            let expected = u64::from(phase1)
                | (u64::from(phase2) << 1)
                | (u64::from(sign) << 2)
                | (l_q << 3)
                | (l_r_prime << (3 + length_width))
                | (l_s << (3 + 2 * length_width));
            let output = apply_scalar(&forward.ops, input);
            assert_eq!(output, expected);
            assert_eq!(output >> data_width, 0);
            assert_eq!(apply_scalar(&inverse.ops, output), input);
            basis_states_checked += 1;
            inverse_pair_checks += 1;
            scratch_clean_checks += 1;
            oracle_transition_checks += 1;
        }
    }
    ReferenceNormalizedControlProofReport {
        length_widths_checked: 3,
        basis_states_checked,
        inverse_pair_checks,
        scratch_clean_checks,
        oracle_transition_checks,
        phase_update9: gate_counts(&build_normalized_phase_update(9, false).ops),
    }
}

fn build_conditional_length_update(
    work_width: usize,
    length_width: usize,
    window: (usize, usize),
    target_r_prime: bool,
    inverse: bool,
) -> B {
    let mut circ = Circuit::new();
    let control = circ.alloc_qreg("rs.length-update.control");
    let work1 = circ.alloc_qreg_bits("rs.length-update.work1", work_width);
    let work2 = circ.alloc_qreg_bits("rs.length-update.work2", work_width);
    let l_t = circ.alloc_qreg_bits("rs.length-update.l-t", length_width);
    let l_q = circ.alloc_qreg_bits("rs.length-update.l-q", length_width);
    let l_s = circ.alloc_qreg_bits("rs.length-update.l-s", length_width);
    let l_r_prime = circ.alloc_qreg_bits("rs.length-update.l-r-prime", length_width);
    let scratch = circ.alloc_qreg_bits(
        "rs.length-update.scratch",
        length_update_scratch_width(length_width),
    );
    let target = if target_r_prime { &l_r_prime } else { &l_t };
    if inverse {
        conditional_length_update_inverse(
            &mut circ, work_width, window.0, window.1, &control, &work1, &work2, &l_s, &l_q,
            &l_r_prime, target, &scratch,
        );
    } else {
        conditional_length_update(
            &mut circ, work_width, window.0, window.1, &control, &work1, &work2, &l_s, &l_q,
            &l_r_prime, target, &scratch,
        );
    }
    circ.into_builder()
}

fn build_work_and_length_swap(work_width: usize, length_width: usize, inverse: bool) -> B {
    let mut circ = Circuit::new();
    let control = circ.alloc_qreg("rs.swap-length.control");
    let iteration = circ.alloc_qreg("rs.swap-length.iteration");
    let work1 = circ.alloc_qreg_bits("rs.swap-length.work1", work_width);
    let work2 = circ.alloc_qreg_bits("rs.swap-length.work2", work_width);
    let l_t = circ.alloc_qreg_bits("rs.swap-length.l-t", length_width);
    let l_t_prime = circ.alloc_qreg_bits("rs.swap-length.l-t-prime", length_width);
    let l_q = circ.alloc_qreg_bits("rs.swap-length.l-q", length_width);
    let l_s = circ.alloc_qreg_bits("rs.swap-length.l-s", length_width);
    let l_r_prime = circ.alloc_qreg_bits("rs.swap-length.l-r-prime", length_width);
    let scratch = Vec::new();
    let window = (1, work_width);
    if inverse {
        conditional_work_and_length_swap_inverse(
            &mut circ, &control, &iteration, &work1, &work2, &l_t, &l_t_prime, &l_q, &l_s,
            &l_r_prime, window, window, &scratch,
        );
    } else {
        conditional_work_and_length_swap(
            &mut circ, &control, &iteration, &work1, &work2, &l_t, &l_t_prime, &l_q, &l_s,
            &l_r_prime, window, window, &scratch,
        );
    }
    circ.into_builder()
}

fn build_work_and_length_swap_quadratic_oracle(work_width: usize, length_width: usize) -> B {
    let mut circ = Circuit::new();
    let control = circ.alloc_qreg("rs.swap-length.control");
    let iteration = circ.alloc_qreg("rs.swap-length.iteration");
    let work1 = circ.alloc_qreg_bits("rs.swap-length.work1", work_width);
    let work2 = circ.alloc_qreg_bits("rs.swap-length.work2", work_width);
    let l_t = circ.alloc_qreg_bits("rs.swap-length.l-t", length_width);
    let l_t_prime = circ.alloc_qreg_bits("rs.swap-length.l-t-prime", length_width);
    let l_q = circ.alloc_qreg_bits("rs.swap-length.l-q", length_width);
    let l_s = circ.alloc_qreg_bits("rs.swap-length.l-s", length_width);
    let l_r_prime = circ.alloc_qreg_bits("rs.swap-length.l-r-prime", length_width);
    conditional_work_and_length_swap_quadratic_oracle(
        &mut circ, &control, &iteration, &work1, &work2, &l_t, &l_t_prime, &l_q, &l_s, &l_r_prime,
    );
    circ.into_builder()
}

fn build_conditional_work_swap(work_width: usize) -> B {
    let mut circ = Circuit::new();
    let control = circ.alloc_qreg("rs.work-swap.control");
    let work1 = circ.alloc_qreg_bits("rs.work-swap.work1", work_width);
    let work2 = circ.alloc_qreg_bits("rs.work-swap.work2", work_width);
    controlled_swap_registers(&mut circ, &control, &work1, &work2);
    circ.into_builder()
}

#[must_use]
pub fn exhaustive_reference_length_swap_check() -> ReferenceLengthSwapProofReport {
    const TEST_LENGTH_WIDTH: usize = 2;

    fn bit_length(value: u64) -> usize {
        if value == 0 {
            0
        } else {
            64 - value.leading_zeros() as usize
        }
    }

    fn pack_work1(width: usize, t: u64, r: u64) -> u64 {
        let l_t = bit_length(t);
        assert!(l_t + 1 + bit_length(r) <= width);
        let mut packed = t;
        for bit in 0..bit_length(r) {
            if (r >> bit) & 1 != 0 {
                packed |= 1u64 << (width - 1 - bit);
            }
        }
        packed
    }

    fn pack_work2(width: usize, t_prime: u64, r_prime: u64) -> u64 {
        let l_r_prime = bit_length(r_prime);
        assert!(bit_length(t_prime) + l_r_prime <= width);
        let mut packed = t_prime;
        for bit in 0..l_r_prime {
            if (r_prime >> bit) & 1 != 0 {
                packed |= 1u64 << (width - 1 - bit);
            }
        }
        packed
    }

    let mut basis_states_checked = 0usize;
    let mut quadratic_oracle_equivalence_checks = 0usize;
    let mut inverse_pair_checks = 0usize;
    let mut scratch_clean_checks = 0usize;
    let mut control_off_identity_checks = 0usize;
    for work_width in 1..=3 {
        let forward = build_work_and_length_swap(work_width, TEST_LENGTH_WIDTH, false);
        let inverse = build_work_and_length_swap(work_width, TEST_LENGTH_WIDTH, true);
        let quadratic_oracle =
            build_work_and_length_swap_quadratic_oracle(work_width, TEST_LENGTH_WIDTH);
        let data_width = 2 + 2 * work_width + 5 * TEST_LENGTH_WIDTH;
        let lengths_offset = 2 + 2 * work_width;
        for input in 0..(1u64 << data_width) {
            if input & 1 != 0 {
                continue;
            }
            let output = apply_scalar(&forward.ops, input);
            assert_eq!(
                output,
                apply_scalar(&quadratic_oracle.ops, input),
                "rotated length swap disagrees with the quadratic oracle"
            );
            quadratic_oracle_equivalence_checks += 1;
            assert_eq!(
                output >> data_width,
                0,
                "length swap left scratch dirty: work_width={work_width} input={input:#x} output={output:#x}"
            );
            assert_eq!(output, input, "control=0 must disable work/length swap");
            control_off_identity_checks += 1;
            assert_eq!(apply_scalar(&inverse.ops, output), input);
            basis_states_checked += 1;
            inverse_pair_checks += 1;
            scratch_clean_checks += 1;
        }

        let limit = 1u64 << work_width;
        for t in 0..limit {
            for t_prime in 0..limit {
                for r in 1..limit {
                    for r_prime in 1..limit {
                        let l_t = bit_length(t);
                        let l_r_prime = bit_length(r_prime);
                        if l_t + 1 + bit_length(r) > work_width
                            || bit_length(t_prime) + 1 + l_r_prime > work_width
                        {
                            continue;
                        }
                        let work1 = pack_work1(work_width, t, r);
                        let work2 = pack_work2(work_width, t_prime, r_prime);
                        for iteration in 0..=1u64 {
                            let input = 1
                                | (iteration << 1)
                                | (work1 << 2)
                                | (work2 << (2 + work_width))
                                | ((l_t as u64) << lengths_offset)
                                | ((bit_length(t_prime) as u64)
                                    << (lengths_offset + TEST_LENGTH_WIDTH))
                                | ((l_r_prime as u64) << (lengths_offset + 4 * TEST_LENGTH_WIDTH));
                            let expected = 1
                                | ((iteration ^ 1) << 1)
                                | (work2 << 2)
                                | (work1 << (2 + work_width))
                                | ((bit_length(t_prime) as u64) << lengths_offset)
                                | ((l_t as u64) << (lengths_offset + TEST_LENGTH_WIDTH))
                                | ((bit_length(r) as u64)
                                    << (lengths_offset + 4 * TEST_LENGTH_WIDTH));
                            let output = apply_scalar(&forward.ops, input);
                            assert_eq!(
                                output,
                                apply_scalar(&quadratic_oracle.ops, input),
                                "rotated length swap disagrees with the quadratic oracle on packed support"
                            );
                            quadratic_oracle_equivalence_checks += 1;
                            assert_eq!(
                                output,
                                expected,
                                "valid packed length swap mismatch: work_width={work_width} t={t} t_prime={t_prime} r={r} r_prime={r_prime}"
                            );
                            assert_eq!(
                                apply_scalar(&inverse.ops, output),
                                input,
                                "valid packed inverse mismatch: work_width={work_width} t={t} t_prime={t_prime} r={r} r_prime={r_prime} iteration={iteration}"
                            );
                            basis_states_checked += 1;
                            inverse_pair_checks += 1;
                            scratch_clean_checks += 1;
                        }
                    }
                }
            }
        }
    }

    let full = build_work_and_length_swap(259, 9, false);
    let standalone_active_qubits = full.active_qubits as usize;
    let standalone_peak_qubits = full.peak_qubits as usize;
    let temporary_peak_qubits = standalone_peak_qubits - standalone_active_qubits;
    let projected_reference_peak_qubits = 815 + temporary_peak_qubits;
    let conditional_work_length_swap259 = measurement_classical_gate_counts(&full.ops);
    let conditional_steps = REFERENCE_STEPS / 4;
    ReferenceLengthSwapProofReport {
        work_widths_checked: 3,
        basis_states_checked,
        quadratic_oracle_equivalence_checks,
        inverse_pair_checks,
        scratch_clean_checks,
        control_off_identity_checks,
        conditional_work_length_swap259,
        conditional_steps,
        standalone_active_qubits,
        standalone_peak_qubits,
        temporary_peak_qubits,
        projected_reference_peak_qubits,
        scheduled_toffoli: conditional_steps * conditional_work_length_swap259.ccx,
    }
}

fn build_full_window_step(n: usize, length_width: usize, inverse: bool) -> B {
    let work_width = n + 3;
    let mut circ = Circuit::new();
    let phase1 = circ.alloc_qreg("rs.whole-step.phase1");
    let phase2 = circ.alloc_qreg("rs.whole-step.phase2");
    let iteration_parity = circ.alloc_qreg("rs.whole-step.iteration");
    let sign = circ.alloc_qreg("rs.whole-step.sign");
    let work1 = circ.alloc_qreg_bits("rs.whole-step.work1", work_width);
    let work2 = circ.alloc_qreg_bits("rs.whole-step.work2", work_width);
    let l_t = circ.alloc_qreg_bits("rs.whole-step.l-t", length_width);
    let l_t_prime = circ.alloc_qreg_bits("rs.whole-step.l-t-prime", length_width);
    let l_q = circ.alloc_qreg_bits("rs.whole-step.l-q", length_width);
    let l_s = circ.alloc_qreg_bits("rs.whole-step.l-s", length_width);
    let l_r_prime = circ.alloc_qreg_bits("rs.whole-step.l-r-prime", length_width);
    if inverse {
        register_shared_full_window_step_inverse(
            &mut circ,
            &phase1,
            &phase2,
            &iteration_parity,
            &sign,
            &work1,
            &work2,
            &l_t,
            &l_t_prime,
            &l_q,
            &l_s,
            &l_r_prime,
            true,
        );
    } else {
        register_shared_full_window_step(
            &mut circ,
            &phase1,
            &phase2,
            &iteration_parity,
            &sign,
            &work1,
            &work2,
            &l_t,
            &l_t_prime,
            &l_q,
            &l_s,
            &l_r_prime,
            true,
        );
    }
    circ.into_builder()
}

fn apply_basis_vector(ops: &[crate::circuit::Op], mut state: Vec<bool>) -> Vec<bool> {
    use crate::circuit::OperationType;

    let index = |id: crate::circuit::QubitId| id.0 as usize;
    for operation in ops {
        match operation.kind {
            OperationType::X => state[index(operation.q_target)] ^= true,
            OperationType::CX => {
                if state[index(operation.q_control1)] {
                    state[index(operation.q_target)] ^= true;
                }
            }
            OperationType::CCX => {
                if state[index(operation.q_control1)] && state[index(operation.q_control2)] {
                    state[index(operation.q_target)] ^= true;
                }
            }
            OperationType::Swap => {
                state.swap(index(operation.q_control1), index(operation.q_target));
            }
            OperationType::R | OperationType::Hmr => {
                state[index(operation.q_target)] = false;
            }
            OperationType::Neg
            | OperationType::Z
            | OperationType::CZ
            | OperationType::CCZ
            | OperationType::PushCondition
            | OperationType::PopCondition
            | OperationType::DebugPrint => {}
            other => panic!("whole-step basis replay saw unsupported operation {other:?}"),
        }
    }
    state
}

fn packed_snapshot_bits(
    snapshot: super::register_shared_eea::RegisterSharedPackedSnapshot,
    length_width: usize,
    total_qubits: usize,
) -> Vec<bool> {
    let work_width = snapshot.n + 3;
    let mut state = vec![false; total_qubits];
    state[0] = snapshot.phase1;
    state[1] = snapshot.phase2;
    state[2] = snapshot.iteration_parity;
    state[3] = snapshot.sign;
    let mut offset = 4usize;
    for bit in 0..work_width {
        state[offset + bit] = ((snapshot.work1 >> bit) & 1) != 0;
    }
    offset += work_width;
    for bit in 0..work_width {
        state[offset + bit] = ((snapshot.work2 >> bit) & 1) != 0;
    }
    offset += work_width;
    let t_prime_length = if snapshot.t_prime == 0 {
        0
    } else {
        64 - snapshot.t_prime.leading_zeros() as usize
    };
    for value in [
        snapshot.l_t,
        t_prime_length,
        snapshot.l_q,
        snapshot.l_s,
        snapshot.l_r_prime,
    ] {
        for bit in 0..length_width {
            state[offset + bit] = ((value >> bit) & 1) != 0;
        }
        offset += length_width;
    }
    state
}

#[must_use]
pub fn exhaustive_reference_whole_step_check() -> ReferenceWholeStepProofReport {
    use super::register_shared_eea::register_shared_small_packed_trace;
    use crate::circuit::OperationType;

    const N: usize = 6;
    const MODULUS: u64 = 37;
    const LENGTH_WIDTH: usize = 5;
    const STEPS: usize = 36;

    let forward = build_full_window_step(N, LENGTH_WIDTH, false);
    let inverse = build_full_window_step(N, LENGTH_WIDTH, true);
    let data_qubits = 4 + 2 * (N + 3) + 5 * LENGTH_WIDTH;
    assert_eq!(forward.active_qubits as usize, data_qubits);
    assert_eq!(inverse.active_qubits as usize, data_qubits);
    let total_qubits = (forward.next_qubit as usize).max(inverse.next_qubit as usize);
    let mut boundary_transitions_checked = 0usize;
    let mut inverse_transition_checks = 0usize;
    let mut scratch_clean_checks = 0usize;

    for input in 1..MODULUS {
        let trace = register_shared_small_packed_trace(input, MODULUS, N, STEPS);
        assert_eq!(trace.len(), STEPS + 1);
        let mut state = packed_snapshot_bits(trace[0], LENGTH_WIDTH, total_qubits);
        for step in 1..=STEPS {
            let before = state.clone();
            let output = apply_basis_vector(&forward.ops, state);
            let expected = packed_snapshot_bits(trace[step], LENGTH_WIDTH, total_qubits);
            assert_eq!(
                &output[..data_qubits],
                &expected[..data_qubits],
                "whole-step mismatch: input={input} step={step} before={:?} expected={:?}",
                trace[step - 1],
                trace[step]
            );
            assert!(
                output[data_qubits..].iter().all(|&bit| !bit),
                "whole-step scratch dirty: input={input} step={step}"
            );
            let restored = apply_basis_vector(&inverse.ops, output.clone());
            assert_eq!(
                restored, before,
                "whole-step inverse mismatch: input={input} step={step}"
            );
            boundary_transitions_checked += 1;
            inverse_transition_checks += 1;
            scratch_clean_checks += 1;
            state = output;
        }
    }

    let emitted_toffoli = forward
        .ops
        .iter()
        .filter(|operation| matches!(operation.kind, OperationType::CCX | OperationType::CCZ))
        .count();
    ReferenceWholeStepProofReport {
        modulus: MODULUS,
        nonzero_inputs_checked: (MODULUS - 1) as usize,
        steps_per_input: STEPS,
        boundary_transitions_checked,
        inverse_transition_checks,
        scratch_clean_checks,
        data_qubits,
        step_active_qubits: forward.active_qubits as usize,
        step_peak_qubits: forward.peak_qubits as usize,
        temporary_peak_qubits: forward.peak_qubits as usize - data_qubits,
        emitted_ops: forward.ops.len(),
        emitted_toffoli,
    }
}

#[must_use]
pub fn profile_reference_scheduled_inversion() -> ReferenceScheduledInversionProfile {
    use crate::circuit::OperationType;

    std::env::set_var("POINT_ADD_COUNT_ONLY", "1");
    let mut circ = Circuit::new();
    let _passenger = circ.alloc_input_qreg_bits("rs.profile.passenger", 257);
    let phase1 = circ.alloc_qreg("rs.profile.phase1");
    let phase2 = circ.alloc_qreg("rs.profile.phase2");
    let iteration_parity = circ.alloc_qreg("rs.profile.iteration");
    let sign = circ.alloc_qreg("rs.profile.sign");
    let work1 = circ.alloc_qreg_bits("rs.profile.work1", 259);
    let work2 = circ.alloc_qreg_bits("rs.profile.work2", 259);
    let l_t = circ.alloc_qreg_bits("rs.profile.l-t", REFERENCE_LENGTH_WIDTH);
    let l_t_prime = circ.alloc_qreg_bits("rs.profile.l-t-prime", REFERENCE_LENGTH_WIDTH);
    let l_q = circ.alloc_qreg_bits("rs.profile.l-q", REFERENCE_LENGTH_WIDTH);
    let l_s = circ.alloc_qreg_bits("rs.profile.l-s", REFERENCE_LENGTH_WIDTH);
    let l_r_prime = circ.alloc_qreg_bits("rs.profile.l-r-prime", REFERENCE_LENGTH_WIDTH);
    let inversion_state_qubits = 2 * 259 + 5 * REFERENCE_LENGTH_WIDTH + 4;
    let passenger_qubits = 257;
    let point_add_state_qubits = inversion_state_qubits + passenger_qubits;
    assert_eq!(circ.b.active_qubits as usize, point_add_state_qubits);

    for step in 1..=REFERENCE_STEPS {
        register_shared_scheduled_step(
            &mut circ,
            step,
            &phase1,
            &phase2,
            &iteration_parity,
            &sign,
            &work1,
            &work2,
            &l_t,
            &l_t_prime,
            &l_q,
            &l_s,
            &l_r_prime,
        );
    }
    let builder = circ.into_builder();
    std::env::remove_var("POINT_ADD_COUNT_ONLY");
    let emitted_toffoli = builder.counted_kind_ops[OperationType::CCX as usize]
        + builder.counted_kind_ops[OperationType::CCZ as usize];
    ReferenceScheduledInversionProfile {
        steps: REFERENCE_STEPS,
        inversion_state_qubits,
        passenger_qubits,
        point_add_state_qubits,
        inversion_peak_qubits: builder.peak_qubits as usize - passenger_qubits,
        projected_point_add_peak_qubits: builder.peak_qubits as usize,
        emitted_ops: builder.counted_kind_ops.iter().sum(),
        emitted_toffoli,
        emitted_hmr: builder.counted_kind_ops[OperationType::Hmr as usize],
        emitted_resets: builder.counted_kind_ops[OperationType::R as usize],
    }
}

fn build_cuccaro(width: usize, inverse: bool) -> B {
    let mut circ = Circuit::new();
    let a = circ.alloc_qreg_bits("rs.cuccaro.a", width);
    let b = circ.alloc_qreg_bits("rs.cuccaro.b", width);
    let carry = circ.alloc_qreg("rs.cuccaro.carry");
    let overflow = circ.alloc_qreg("rs.cuccaro.overflow");
    if inverse {
        cuccaro_sub_mod_2n(&mut circ, &a, &b, &carry, &overflow);
    } else {
        cuccaro_add_mod_2n(&mut circ, &a, &b, &carry, &overflow);
    }
    circ.into_builder()
}

#[must_use]
pub fn exhaustive_reference_cuccaro_check() -> ReferenceCuccaroProofReport {
    let mut basis_states_checked = 0usize;
    let mut carry_clean_checks = 0usize;
    let mut inverse_pair_checks = 0usize;
    for width in 1..=7 {
        let add = build_cuccaro(width, false);
        let sub = build_cuccaro(width, true);
        let mask = (1u64 << width) - 1;
        for a in 0..=mask {
            for b in 0..=mask {
                let input = a | (b << width);
                let sum = a + b;
                let expected = a | ((sum & mask) << width) | ((sum >> width) << (2 * width + 1));
                let output = apply_scalar(&add.ops, input);
                assert_eq!(output, expected);
                assert_eq!((output >> (2 * width)) & 1, 0);
                assert_eq!(apply_scalar(&sub.ops, output), input);
                basis_states_checked += 1;
                carry_clean_checks += 1;
                inverse_pair_checks += 1;
            }
        }
    }
    ReferenceCuccaroProofReport {
        widths_checked: 7,
        basis_states_checked,
        carry_clean_checks,
        inverse_pair_checks,
        add9: gate_counts(&build_cuccaro(REFERENCE_LENGTH_WIDTH, false).ops),
        sub9: gate_counts(&build_cuccaro(REFERENCE_LENGTH_WIDTH, true).ops),
    }
}

fn build_location_swap(work_width: usize, length_width: usize, inverse: bool) -> B {
    let mut circ = Circuit::new();
    let phase1 = circ.alloc_qreg("rs.location.phase1");
    let phase2 = circ.alloc_qreg("rs.location.phase2");
    let sign = circ.alloc_qreg("rs.location.sign");
    let work1 = circ.alloc_qreg_bits("rs.location.work1", work_width);
    let l_t = circ.alloc_qreg_bits("rs.location.lt", length_width);
    let l_q = circ.alloc_qreg_bits("rs.location.lq", length_width);
    let scratch = circ.alloc_qreg_bits("rs.location.scratch", length_width + 2);
    if inverse {
        location_controlled_swap_one_hot_inverse(
            &mut circ, &phase1, &phase2, &sign, &work1, 0, &l_t, &l_q, &scratch,
        );
    } else {
        location_controlled_swap_one_hot(
            &mut circ, &phase1, &phase2, &sign, &work1, 0, &l_t, &l_q, &scratch,
        );
    }
    circ.into_builder()
}

#[must_use]
pub fn exhaustive_reference_location_swap_check() -> ReferenceLocationSwapProofReport {
    const TEST_LENGTH_WIDTH: usize = 3;
    let mut basis_states_checked = 0usize;
    let mut scratch_clean_checks = 0usize;
    let mut inverse_pair_checks = 0usize;
    for work_width in 1..=6 {
        let forward = build_location_swap(work_width, TEST_LENGTH_WIDTH, false);
        let inverse = build_location_swap(work_width, TEST_LENGTH_WIDTH, true);
        let work_mask = (1u64 << work_width) - 1;
        let length_mask = (1u64 << TEST_LENGTH_WIDTH) - 1;
        let data_width = 3 + work_width + 2 * TEST_LENGTH_WIDTH;
        for input in 0..(1u64 << data_width) {
            let phase1 = (input & 1) != 0;
            let phase2 = ((input >> 1) & 1) != 0;
            let mut sign = ((input >> 2) & 1) != 0;
            let mut work = (input >> 3) & work_mask;
            let l_t = (input >> (3 + work_width)) & length_mask;
            let mut l_q = (input >> (3 + work_width + TEST_LENGTH_WIDTH)) & length_mask;
            let active = phase1 ^ phase2;
            if active && !phase1 {
                l_q = l_q.wrapping_add(1) & length_mask;
            }
            let location = (l_t + l_q) & length_mask;
            if active && location < work_width as u64 {
                let work_bit = ((work >> location) & 1) != 0;
                if work_bit != sign {
                    work ^= 1u64 << location;
                    sign = work_bit;
                }
            }
            if active && phase1 {
                l_q = l_q.wrapping_sub(1) & length_mask;
            }
            let expected = u64::from(phase1)
                | (u64::from(phase2) << 1)
                | (u64::from(sign) << 2)
                | (work << 3)
                | (l_t << (3 + work_width))
                | (l_q << (3 + work_width + TEST_LENGTH_WIDTH));
            let output = apply_scalar(&forward.ops, input);
            assert_eq!(output, expected);
            assert_eq!(output >> data_width, 0);
            assert_eq!(apply_scalar(&inverse.ops, output), input);
            basis_states_checked += 1;
            scratch_clean_checks += 1;
            inverse_pair_checks += 1;
        }
    }

    let full259_length9 = gate_counts(&build_location_swap(259, REFERENCE_LENGTH_WIDTH, false).ops);
    let full259_length9_inverse =
        gate_counts(&build_location_swap(259, REFERENCE_LENGTH_WIDTH, true).ops);
    let schedule = exhaustive_reference_schedule_check();
    let fixed_toffoli = full259_length9.ccx - 35 * 259;
    ReferenceLocationSwapProofReport {
        work_widths_checked: 6,
        basis_states_checked,
        scratch_clean_checks,
        inverse_pair_checks,
        full259_length9,
        full259_length9_inverse,
        reference_steps: REFERENCE_STEPS,
        full_window_toffoli_upper_bound: REFERENCE_STEPS * full259_length9.ccx,
        scheduled_window_sum: schedule.quotient_swap_window_sum,
        scheduled_toffoli: 35 * schedule.quotient_swap_window_sum + fixed_toffoli * REFERENCE_STEPS,
    }
}

/// Emit and simulate the exact reversible boundary map between canonical
/// challenge inputs and the public register-shared EEA initial layout.
///
/// This deliberately excludes the 1,479 EEA steps. It proves that reflection,
/// bit-length extraction, register packing, and their inverse fit beneath the
/// 912-qubit reference allocation and clean every non-input lane.
#[doc(hidden)]
#[must_use]
pub fn reference_initializer_roundtrip_check() -> ReferenceInitializerProofReport {
    use super::shrunken_pz_state_machine::{bit_length_lean, controlled_field_neg};
    use crate::circuit::{OperationType, QubitId};
    use crate::point_add::trailmix_port::arith::compare::compare_geq_const;
    use crate::point_add::trailmix_port::mod_arith::SECP256K1_P_LE;
    use crate::point_add::SECP256K1_P;
    use crate::sim::Simulator;
    use ruint::aliases::U256;
    use sha3::{
        digest::{ExtendableOutput, Update, XofReader},
        Shake128,
    };

    const WORK_BITS: usize = 259;
    const FIELD_BITS: usize = 257;
    const PASSENGER_BITS: usize = 257;
    const CONTROL_BITS: usize = 4;
    const SCRATCH_BITS: usize = 97;
    const PACKED_QUBITS: usize = 912;
    const HALF_BYTES: [u8; 33] = [
        0x18, 0xfe, 0xff, 0x7f, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
        0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
        0xff, 0x7f, 0x00,
    ];

    fn ids(register: &[QReg]) -> Vec<u32> {
        register.iter().map(QReg::id).collect()
    }

    fn load<R: XofReader>(
        simulator: &mut Simulator<'_, R>,
        register: &[u32],
        value: U256,
        shot: usize,
    ) {
        for (bit, &id) in register.iter().take(256).enumerate() {
            if value.bit(bit) {
                *simulator.qubit_mut(QubitId(u64::from(id))) |= 1u64 << shot;
            }
        }
    }

    fn read<R: XofReader>(simulator: &Simulator<'_, R>, register: &[u32], shot: usize) -> U256 {
        let mut value = U256::ZERO;
        for (bit, &id) in register.iter().take(256).enumerate() {
            if ((simulator.qubit(QubitId(u64::from(id))) >> shot) & 1) != 0 {
                value.set_bit(bit, true);
            }
        }
        value
    }

    fn free_register(circ: &mut Circuit, register: Vec<QReg>) {
        for lane in register {
            circ.zero_and_free(lane);
        }
    }

    std::env::remove_var("POINT_ADD_COUNT_ONLY");
    let mut circ = Circuit::new();
    assert!(!circ.b.count_only);
    let mut dx = circ.alloc_input_qreg_bits("q925-init.dx", FIELD_BITS);
    let dy = circ.alloc_input_qreg_bits("q925-init.dy", PASSENGER_BITS);
    let dx_ids = ids(&dx);
    let dy_ids = ids(&dy);
    assert_eq!(circ.b.active_qubits as usize, 2 * FIELD_BITS);

    // The paper initializes the iteration parity with the reflection bit and
    // runs EEA on |dx| <= p/2.
    let iteration_parity = circ.alloc_qreg("q925-init.iteration-parity");
    compare_geq_const(&mut circ, &dx, &HALF_BYTES, &iteration_parity);
    controlled_field_neg(&mut circ, &iteration_parity, &dx);

    let l_r_prime = circ.alloc_qreg_bits("q925-init.l-r-prime", REFERENCE_LENGTH_WIDTH);
    let source: Vec<&QReg> = dx.iter().take(256).collect();
    bit_length_lean(&mut circ, &source, &l_r_prime, false);
    let initialization_transient_peak_qubits = circ.b.peak_qubits as usize;

    // dx is little-endian. Two clean high pads followed by handle reversal
    // produce the public Work2 layout: zero-padded t' followed by big-endian r'.
    dx.push(circ.alloc_qreg("q925-init.work2-pad0"));
    dx.push(circ.alloc_qreg("q925-init.work2-pad1"));
    dx.reverse();
    let mut work2 = dx;
    assert_eq!(work2.len(), WORK_BITS);

    // Work1 = [t_le | q_be | r_be] for (t,q,r)=(1,0,p). Position 1 is the
    // delimiter; position 2 is the leading zero of the 257-bit p field.
    let work1 = circ.alloc_qreg_bits("q925-init.work1", WORK_BITS);
    circ.x(&work1[0]);
    for bit in 0..256 {
        if ((SECP256K1_P_LE[bit / 8] >> (bit % 8)) & 1) != 0 {
            circ.x(&work1[WORK_BITS - 1 - bit]);
        }
    }

    let l_t = circ.alloc_qreg_bits("q925-init.l-t", REFERENCE_LENGTH_WIDTH);
    let l_q = circ.alloc_qreg_bits("q925-init.l-q", REFERENCE_LENGTH_WIDTH);
    let l_s = circ.alloc_qreg_bits("q925-init.l-s", REFERENCE_LENGTH_WIDTH);
    circ.x(&l_t[0]);
    let phase1 = circ.alloc_qreg("q925-init.phase1");
    let phase2 = circ.alloc_qreg("q925-init.phase2");
    let sign = circ.alloc_qreg("q925-init.sign");
    let scratch = circ.alloc_qreg_bits("q925-init.reference-scratch", SCRATCH_BITS);
    let packed_active_qubits = circ.b.active_qubits as usize;
    let packed_peak_qubits = circ.b.peak_qubits as usize;
    assert_eq!(CONTROL_BITS, 1 + 3);
    assert_eq!(packed_active_qubits, PACKED_QUBITS);
    assert_eq!(packed_peak_qubits, PACKED_QUBITS);

    // Inverse boundary map. The untouched all-zero registers are released,
    // constants are toggled away, and the reflected input is restored.
    free_register(&mut circ, scratch);
    circ.zero_and_free(phase1);
    circ.zero_and_free(phase2);
    circ.zero_and_free(sign);
    circ.x(&l_t[0]);
    free_register(&mut circ, l_t);
    free_register(&mut circ, l_q);
    free_register(&mut circ, l_s);

    circ.x(&work1[0]);
    for bit in 0..256 {
        if ((SECP256K1_P_LE[bit / 8] >> (bit % 8)) & 1) != 0 {
            circ.x(&work1[WORK_BITS - 1 - bit]);
        }
    }
    free_register(&mut circ, work1);

    work2.reverse();
    let pad1 = work2.pop().expect("Work2 pad1");
    let pad0 = work2.pop().expect("Work2 pad0");
    circ.zero_and_free(pad1);
    circ.zero_and_free(pad0);
    assert_eq!(work2.len(), FIELD_BITS);
    dx = work2;

    let source: Vec<&QReg> = dx.iter().take(256).collect();
    bit_length_lean(&mut circ, &source, &l_r_prime, true);
    free_register(&mut circ, l_r_prime);
    controlled_field_neg(&mut circ, &iteration_parity, &dx);
    compare_geq_const(&mut circ, &dx, &HALF_BYTES, &iteration_parity);
    circ.zero_and_free(iteration_parity);
    circ.flush_pending_frees();
    let final_active_qubits = circ.b.active_qubits as usize;
    assert_eq!(final_active_qubits, 2 * FIELD_BITS);

    let builder = circ.into_builder();
    let emitted_toffoli = builder.counted_kind_ops[OperationType::CCX as usize]
        + builder.counted_kind_ops[OperationType::CCZ as usize];
    let emitted_hmr = builder.counted_kind_ops[OperationType::Hmr as usize];
    let emitted_resets = builder.counted_kind_ops[OperationType::R as usize];

    let half = SECP256K1_P >> 1;
    let mut cases = Vec::with_capacity(64);
    for shot in 0..64usize {
        let dx_value = match shot {
            0 => U256::from(1u64),
            1 => half,
            2 => half + U256::from(1u64),
            3 => SECP256K1_P - U256::from(1u64),
            _ => {
                let high_bit = 1 + ((37 * shot) % 255);
                let low = (U256::from(1u64) << high_bit) | U256::from((2 * shot + 1) as u64);
                if shot & 1 == 0 {
                    low
                } else {
                    SECP256K1_P - low
                }
            }
        };
        let dy_small = U256::from((5 * shot + 3) as u64);
        let dy_value = if shot % 3 == 0 {
            SECP256K1_P - dy_small
        } else {
            dy_small
        };
        assert!(dx_value != U256::ZERO && dx_value < SECP256K1_P);
        assert!(dy_value < SECP256K1_P);
        cases.push((dx_value, dy_value));
    }

    let mut seed = Shake128::default();
    seed.update(b"q925-register-shared-initializer-roundtrip");
    let mut xof = seed.finalize_xof();
    let mut simulator = Simulator::new(
        builder.next_qubit as usize,
        builder.next_bit as usize,
        &mut xof,
    );
    simulator.clear_for_shot();
    for (shot, &(dx_value, dy_value)) in cases.iter().enumerate() {
        load(&mut simulator, &dx_ids, dx_value, shot);
        load(&mut simulator, &dy_ids, dy_value, shot);
    }
    simulator.apply_iter(builder.ops.iter());

    let mut reflected_cases_checked = 0usize;
    for (shot, &(dx_value, dy_value)) in cases.iter().enumerate() {
        assert_eq!(read(&simulator, &dx_ids, shot), dx_value);
        assert_eq!(read(&simulator, &dy_ids, shot), dy_value);
        assert_eq!(
            (simulator.qubit(QubitId(u64::from(dx_ids[256]))) >> shot) & 1,
            0
        );
        assert_eq!(
            (simulator.qubit(QubitId(u64::from(dy_ids[256]))) >> shot) & 1,
            0
        );
        reflected_cases_checked += usize::from(dx_value > half);
    }
    assert_eq!(
        simulator.phase, 0,
        "initializer roundtrip left phase garbage"
    );

    for id in dx_ids.iter().chain(dy_ids.iter()) {
        *simulator.qubit_mut(QubitId(u64::from(*id))) = 0;
    }
    for id in 0..builder.next_qubit {
        assert_eq!(
            simulator.qubit(QubitId(u64::from(id))),
            0,
            "initializer roundtrip left q{id} dirty"
        );
    }

    ReferenceInitializerProofReport {
        cases_checked: cases.len(),
        reflected_cases_checked,
        non_reflected_cases_checked: cases.len() - reflected_cases_checked,
        input_qubits: 2 * FIELD_BITS,
        initialization_transient_peak_qubits,
        packed_peak_qubits,
        packed_active_qubits,
        final_active_qubits,
        emitted_ops: builder.ops.len(),
        emitted_toffoli,
        emitted_hmr,
        emitted_resets,
        classical_roundtrip_checks: cases.len(),
        phase_cleanup_checks: cases.len(),
        ancilla_cleanup_checks: cases.len(),
    }
}
