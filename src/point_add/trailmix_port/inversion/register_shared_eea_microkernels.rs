//! Gate-level primitives and liveness checks for the register-sharing EEA.
//!
//! The module contains the allocation-free rotations, length arithmetic, and
//! control primitives used by the staged port. No complete inversion or
//! challenge-width claim follows from these component proofs.

use crate::circuit::{Op, OperationType};
use crate::point_add::trailmix_port::circuit::{Circuit, QReg};
use crate::point_add::B;

pub const WORK_BITS: usize = 259;
pub const FIELD_PASSENGER_BITS: usize = 257;
pub const LENGTH_BITS: [usize; 4] = [10, 10, 10, 11];
pub const REFERENCE_LENGTH_BITS: [usize; 4] = [9, 9, 9, 9];
pub const REFERENCE_SCRATCH_POOLS: [usize; 11] = [11, 13, 13, 2, 22, 1, 22, 1, 3, 5, 4];
pub const CONTROL_BITS: usize = 4;
pub const PROJECTED_INVERSION_PEAK: usize =
    2 * WORK_BITS + FIELD_PASSENGER_BITS + 41 + CONTROL_BITS;
pub const REFERENCE_PORT_PEAK: usize =
    2 * WORK_BITS + FIELD_PASSENGER_BITS + 36 + CONTROL_BITS + 97;

/// Explicit proof/profile override for the register-shared decrement stream.
pub const REGISTER_SHARED_REVERSE_DECREMENT_STREAM_ENV: &str =
    "LOWQ_REGISTER_SHARED_REVERSE_DECREMENT_STREAM";

/// Source-bake point for the exact reverse-decrement stream.
///
/// Keep this off until the focused proof and the count-only production profile
/// have both passed. An explicit `0` or `1` environment value overrides this
/// fallback for isolated proof/profile processes.
pub const REGISTER_SHARED_REVERSE_DECREMENT_STREAM_SOURCE_BAKE: bool = false;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ReversibleDecrementStream {
    LegacyComplementedBorrow,
    ExactReverseIncrement,
}

#[must_use]
pub fn production_decrement_stream() -> ReversibleDecrementStream {
    static STREAM: std::sync::OnceLock<ReversibleDecrementStream> = std::sync::OnceLock::new();

    *STREAM.get_or_init(|| {
        let enabled = match std::env::var(REGISTER_SHARED_REVERSE_DECREMENT_STREAM_ENV) {
            Ok(value) => match value.as_str() {
                "0" => false,
                "1" => true,
                _ => panic!(
                    "{REGISTER_SHARED_REVERSE_DECREMENT_STREAM_ENV} must be exactly 0 or 1, got {value:?}"
                ),
            },
            Err(std::env::VarError::NotPresent) => {
                REGISTER_SHARED_REVERSE_DECREMENT_STREAM_SOURCE_BAKE
            }
            Err(error) => {
                panic!("invalid {REGISTER_SHARED_REVERSE_DECREMENT_STREAM_ENV}: {error}")
            }
        };
        if enabled {
            ReversibleDecrementStream::ExactReverseIncrement
        } else {
            ReversibleDecrementStream::LegacyComplementedBorrow
        }
    })
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct RegisterSharedGateCounts {
    pub x: usize,
    pub cx: usize,
    pub ccx: usize,
    pub total: usize,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct RegisterSharedShiftWidthReport {
    pub width: usize,
    pub basis_states_checked: usize,
    pub forward: RegisterSharedGateCounts,
    pub reverse: RegisterSharedGateCounts,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct RegisterSharedShiftProofReport {
    pub widths_checked: usize,
    pub basis_states_checked: usize,
    pub allocation_free_streams_checked: usize,
    pub phase_clean_streams_checked: usize,
    pub widths: Vec<RegisterSharedShiftWidthReport>,
    pub work259_shift_toffoli: usize,
    pub legacy_schematic_four_rotations_toffoli: usize,
    pub reference_steps: usize,
    pub legacy_schematic_four_rotations_total_toffoli: usize,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct RegisterSharedLengthProofReport {
    pub widths_checked: usize,
    pub basis_states_checked: usize,
    pub scratch_clean_checks: usize,
    pub inverse_pair_checks: usize,
    pub increment9: RegisterSharedGateCounts,
    pub decrement9: RegisterSharedGateCounts,
    pub controlled_increment9: RegisterSharedGateCounts,
    pub controlled_decrement9: RegisterSharedGateCounts,
    pub two_control_add_one9: RegisterSharedGateCounts,
    pub two_control_sub_one9: RegisterSharedGateCounts,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct RegisterSharedReverseDecrementWidthReport {
    pub width: usize,
    pub direct_basis_states: usize,
    pub controlled_basis_states: usize,
    pub increment: RegisterSharedGateCounts,
    pub legacy_decrement: RegisterSharedGateCounts,
    pub exact_reverse_decrement: RegisterSharedGateCounts,
    pub controlled_increment: RegisterSharedGateCounts,
    pub legacy_controlled_decrement: RegisterSharedGateCounts,
    pub exact_reverse_controlled_decrement: RegisterSharedGateCounts,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct RegisterSharedReverseDecrementProofReport {
    pub widths_checked: usize,
    pub direct_basis_states_checked: usize,
    pub controlled_basis_states_checked: usize,
    pub exact_reverse_stream_checks: usize,
    pub scalar_forward_checks: usize,
    pub inverse_pair_checks: usize,
    pub phase_clean_streams_checked: usize,
    pub ancilla_clean_checks: usize,
    pub allocation_profile_checks: usize,
    pub toffoli_preservation_checks: usize,
    pub local_legacy_ops: usize,
    pub local_exact_reverse_ops: usize,
    pub local_ops_removed: usize,
    pub widths: Vec<RegisterSharedReverseDecrementWidthReport>,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct RegisterSharedControlProofReport {
    pub work_swap_widths_checked: usize,
    pub work_swap_basis_states_checked: usize,
    pub work_swap_inverse_checks: usize,
    pub work259_swap: RegisterSharedGateCounts,
    pub phase_update_basis_states_checked: usize,
    pub phase_update_scratch_clean_checks: usize,
    pub phase_update: RegisterSharedGateCounts,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct RegisterSharedPrePostProofReport {
    pub work_widths_checked: usize,
    pub basis_states_checked: usize,
    pub scratch_clean_checks: usize,
    pub inverse_pair_checks: usize,
    pub pre_shift259_length9: RegisterSharedGateCounts,
    pub pre_shift259_length9_inverse: RegisterSharedGateCounts,
    pub post_shift259_length9: RegisterSharedGateCounts,
    pub post_shift259_length9_inverse: RegisterSharedGateCounts,
    pub emitted_rotations_per_step: usize,
    pub reference_steps: usize,
    pub total_pre_post_toffoli: usize,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct RegisterSharedBarrelProofReport {
    pub widths_checked: usize,
    pub basis_states_checked: usize,
    pub inverse_pair_checks: usize,
    pub work259_amount9_high: RegisterSharedGateCounts,
    pub work259_amount9_low: RegisterSharedGateCounts,
    pub toffoli_per_direction: usize,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct RegisterSharedAllocationReport {
    pub paper_core_peak_qubits: usize,
    pub paper_core_final_active_qubits: usize,
    pub paper_core_reset_operations: usize,
    pub reference_port_peak_qubits: usize,
    pub reference_port_final_active_qubits: usize,
    pub reference_port_reset_operations: usize,
    pub reference_port_scratch_bits: usize,
    pub input_dx_bits: usize,
    pub passenger_bits: usize,
    pub work1_bits: usize,
    pub work2_bits: usize,
    pub paper_length_bits: usize,
    pub reference_length_bits: usize,
    pub control_bits: usize,
}

/// Controlled cyclic rotation toward lower wire indices.
pub fn controlled_rotate_low(circ: &mut Circuit, control: &QReg, register: &[QReg]) {
    for index in 0..register.len().saturating_sub(1) {
        circ.cswap(control, &register[index], &register[index + 1]);
    }
}

/// Exact inverse of [`controlled_rotate_low`].
pub fn controlled_rotate_high(circ: &mut Circuit, control: &QReg, register: &[QReg]) {
    for index in (0..register.len().saturating_sub(1)).rev() {
        circ.cswap(control, &register[index], &register[index + 1]);
    }
}

fn controlled_rotate_high_two(circ: &mut Circuit, control: &QReg, register: &[QReg]) {
    for index in (0..register.len().saturating_sub(2)).rev() {
        circ.cswap(control, &register[index], &register[index + 1]);
        circ.cswap(control, &register[index + 1], &register[index + 2]);
    }
}

fn controlled_rotate_high_two_inverse(circ: &mut Circuit, control: &QReg, register: &[QReg]) {
    for index in 0..register.len().saturating_sub(2) {
        circ.cswap(control, &register[index + 1], &register[index + 2]);
        circ.cswap(control, &register[index], &register[index + 1]);
    }
}

fn gcd(mut left: usize, mut right: usize) -> usize {
    while right != 0 {
        let remainder = left % right;
        left = right;
        right = remainder;
    }
    left
}

/// Controlled cyclic rotation toward higher wire indices by an arbitrary offset.
/// Each permutation cycle is emitted with a pivot and `cycle_len - 1` Fredkins.
pub fn controlled_rotate_high_by(
    circ: &mut Circuit,
    control: &QReg,
    register: &[QReg],
    offset: usize,
) {
    let width = register.len();
    if width < 2 {
        return;
    }
    let offset = offset % width;
    if offset == 0 {
        return;
    }
    let cycles = gcd(width, offset);
    for start in 0..cycles {
        let mut current = (start + offset) % width;
        while current != start {
            circ.cswap(control, &register[start], &register[current]);
            current = (current + offset) % width;
        }
    }
}

/// Exact inverse of [`controlled_rotate_high_by`].
pub fn controlled_rotate_low_by(
    circ: &mut Circuit,
    control: &QReg,
    register: &[QReg],
    offset: usize,
) {
    let width = register.len();
    if width < 2 {
        return;
    }
    controlled_rotate_high_by(circ, control, register, width - (offset % width));
}

/// Reference-view variant of [`controlled_rotate_high_by`].
///
/// This permits the same allocation-free permutation to act on a reversed or
/// otherwise borrowed view of a packed register.
pub fn controlled_rotate_high_by_refs(
    circ: &mut Circuit,
    control: &QReg,
    register: &[&QReg],
    offset: usize,
) {
    let width = register.len();
    if width < 2 {
        return;
    }
    let offset = offset % width;
    if offset == 0 {
        return;
    }
    let cycles = gcd(width, offset);
    for start in 0..cycles {
        let mut current = (start + offset) % width;
        while current != start {
            circ.cswap(control, register[start], register[current]);
            current = (current + offset) % width;
        }
    }
}

/// Exact inverse of [`controlled_rotate_high_by_refs`].
pub fn controlled_rotate_low_by_refs(
    circ: &mut Circuit,
    control: &QReg,
    register: &[&QReg],
    offset: usize,
) {
    let width = register.len();
    if width < 2 {
        return;
    }
    controlled_rotate_high_by_refs(circ, control, register, width - (offset % width));
}

/// Rotate by the little-endian quantum amount without allocating scratch.
pub fn variable_rotate_high(circ: &mut Circuit, amount: &[QReg], register: &[QReg]) {
    if register.len() < 2 {
        return;
    }
    let mut offset = 1 % register.len();
    for control in amount {
        controlled_rotate_high_by(circ, control, register, offset);
        offset = (2 * offset) % register.len();
    }
}

/// Exact inverse of [`variable_rotate_high`].
pub fn variable_rotate_low(circ: &mut Circuit, amount: &[QReg], register: &[QReg]) {
    if register.len() < 2 {
        return;
    }
    let mut offsets = Vec::with_capacity(amount.len());
    let mut offset = 1 % register.len();
    for _ in amount {
        offsets.push(offset);
        offset = (2 * offset) % register.len();
    }
    for (control, offset) in amount.iter().zip(offsets).rev() {
        controlled_rotate_low_by(circ, control, register, offset);
    }
}

/// Reference-view variant of [`variable_rotate_high`].
pub fn variable_rotate_high_refs(circ: &mut Circuit, amount: &[QReg], register: &[&QReg]) {
    if register.len() < 2 {
        return;
    }
    let mut offset = 1 % register.len();
    for control in amount {
        controlled_rotate_high_by_refs(circ, control, register, offset);
        offset = (2 * offset) % register.len();
    }
}

/// Exact inverse of [`variable_rotate_high_refs`].
pub fn variable_rotate_low_refs(circ: &mut Circuit, amount: &[QReg], register: &[&QReg]) {
    if register.len() < 2 {
        return;
    }
    let mut offsets = Vec::with_capacity(amount.len());
    let mut offset = 1 % register.len();
    for _ in amount {
        offsets.push(offset);
        offset = (2 * offset) % register.len();
    }
    for (control, offset) in amount.iter().zip(offsets).rev() {
        controlled_rotate_low_by_refs(circ, control, register, offset);
    }
}

/// Multi-controlled X using a clean v-chain. Every ancilla is restored to zero.
pub fn multi_controlled_x_vchain(
    circ: &mut Circuit,
    controls: &[&QReg],
    target: &QReg,
    ancillas: &[QReg],
) {
    match controls.len() {
        0 => circ.x(target),
        1 => circ.cx(controls[0], target),
        2 => circ.ccx(controls[0], controls[1], target),
        count => {
            assert!(ancillas.len() >= count - 2);
            circ.ccx(controls[0], controls[1], &ancillas[0]);
            for index in 2..count - 1 {
                circ.ccx(controls[index], &ancillas[index - 2], &ancillas[index - 1]);
            }
            circ.ccx(controls[count - 1], &ancillas[count - 3], target);
            for index in (2..count - 1).rev() {
                circ.ccx(controls[index], &ancillas[index - 2], &ancillas[index - 1]);
            }
            circ.ccx(controls[0], controls[1], &ancillas[0]);
        }
    }
}

/// Add one modulo `2^n`, restoring the `n-1` clean carry lanes.
pub fn increment_mod_2n(circ: &mut Circuit, register: &[QReg], carries: &[QReg]) {
    let width = register.len();
    if width == 0 {
        return;
    }
    if width == 1 {
        circ.x(&register[0]);
        return;
    }
    assert!(carries.len() >= width - 1);
    circ.cx(&register[0], &carries[0]);
    for index in 1..width - 1 {
        circ.ccx(&register[index], &carries[index - 1], &carries[index]);
    }
    circ.cx(&carries[width - 2], &register[width - 1]);
    for index in (1..width - 1).rev() {
        circ.ccx(&register[index], &carries[index - 1], &carries[index]);
        circ.cx(&carries[index - 1], &register[index]);
    }
    circ.cx(&register[0], &carries[0]);
    circ.x(&register[0]);
}

/// Subtract one modulo `2^n`, restoring the `n-1` clean borrow lanes.
pub fn decrement_mod_2n(circ: &mut Circuit, register: &[QReg], borrows: &[QReg]) {
    let width = register.len();
    if width == 0 {
        return;
    }
    if width == 1 {
        circ.x(&register[0]);
        return;
    }
    assert!(borrows.len() >= width - 1);
    circ.x(&register[0]);
    circ.cx(&register[0], &borrows[0]);
    circ.x(&register[0]);
    for index in 1..width - 1 {
        circ.x(&register[index]);
        circ.ccx(&register[index], &borrows[index - 1], &borrows[index]);
        circ.x(&register[index]);
    }
    circ.cx(&borrows[width - 2], &register[width - 1]);
    for index in (1..width - 1).rev() {
        circ.x(&register[index]);
        circ.ccx(&register[index], &borrows[index - 1], &borrows[index]);
        circ.x(&register[index]);
        circ.cx(&borrows[index - 1], &register[index]);
    }
    circ.x(&register[0]);
    circ.cx(&register[0], &borrows[0]);
    circ.x(&register[0]);
    circ.x(&register[0]);
}

/// Subtract one by emitting the exact reversed operation stream of
/// [`increment_mod_2n`].
pub fn decrement_mod_2n_exact_reverse_increment(
    circ: &mut Circuit,
    register: &[QReg],
    carries: &[QReg],
) {
    let width = register.len();
    if width == 0 {
        return;
    }
    if width == 1 {
        circ.x(&register[0]);
        return;
    }
    assert!(carries.len() >= width - 1);
    circ.x(&register[0]);
    circ.cx(&register[0], &carries[0]);
    for index in 1..width - 1 {
        circ.cx(&carries[index - 1], &register[index]);
        circ.ccx(&register[index], &carries[index - 1], &carries[index]);
    }
    circ.cx(&carries[width - 2], &register[width - 1]);
    for index in (1..width - 1).rev() {
        circ.ccx(&register[index], &carries[index - 1], &carries[index]);
    }
    circ.cx(&register[0], &carries[0]);
}

pub fn decrement_mod_2n_with_stream(
    circ: &mut Circuit,
    register: &[QReg],
    scratch: &[QReg],
    stream: ReversibleDecrementStream,
) {
    match stream {
        ReversibleDecrementStream::LegacyComplementedBorrow => {
            decrement_mod_2n(circ, register, scratch);
        }
        ReversibleDecrementStream::ExactReverseIncrement => {
            decrement_mod_2n_exact_reverse_increment(circ, register, scratch);
        }
    }
}

/// Production dispatcher used by the complete register-shared EEA path.
pub fn production_decrement_mod_2n(circ: &mut Circuit, register: &[QReg], scratch: &[QReg]) {
    decrement_mod_2n_with_stream(circ, register, scratch, production_decrement_stream());
}

/// Controlled add one modulo `2^n`, restoring the clean carry lanes.
pub fn controlled_increment_mod_2n(
    circ: &mut Circuit,
    control: &QReg,
    register: &[QReg],
    carries: &[QReg],
) {
    let width = register.len();
    if width == 0 {
        return;
    }
    if width == 1 {
        circ.cx(control, &register[0]);
        return;
    }
    assert!(carries.len() >= width - 1);
    circ.ccx(control, &register[0], &carries[0]);
    circ.cx(control, &register[0]);
    for index in 1..width - 1 {
        circ.ccx(&register[index], &carries[index - 1], &carries[index]);
    }
    circ.cx(&carries[width - 2], &register[width - 1]);
    for index in (1..width - 1).rev() {
        circ.ccx(&register[index], &carries[index - 1], &carries[index]);
        circ.cx(&carries[index - 1], &register[index]);
    }
    circ.cx(control, &register[0]);
    circ.ccx(control, &register[0], &carries[0]);
    circ.cx(control, &register[0]);
}

/// Controlled subtract one modulo `2^n`, restoring the clean borrow lanes.
pub fn controlled_decrement_mod_2n(
    circ: &mut Circuit,
    control: &QReg,
    register: &[QReg],
    borrows: &[QReg],
) {
    let width = register.len();
    if width == 0 {
        return;
    }
    if width == 1 {
        circ.cx(control, &register[0]);
        return;
    }
    assert!(borrows.len() >= width - 1);
    circ.x(&register[0]);
    circ.ccx(control, &register[0], &borrows[0]);
    circ.x(&register[0]);
    circ.cx(control, &register[0]);
    for index in 1..width - 1 {
        circ.x(&register[index]);
        circ.ccx(&register[index], &borrows[index - 1], &borrows[index]);
        circ.x(&register[index]);
    }
    circ.cx(&borrows[width - 2], &register[width - 1]);
    for index in (1..width - 1).rev() {
        circ.x(&register[index]);
        circ.ccx(&register[index], &borrows[index - 1], &borrows[index]);
        circ.x(&register[index]);
        circ.cx(&borrows[index - 1], &register[index]);
    }
    circ.cx(control, &register[0]);
    circ.x(&register[0]);
    circ.ccx(control, &register[0], &borrows[0]);
    circ.x(&register[0]);
    circ.cx(control, &register[0]);
}

/// Controlled subtract one emitted as the exact reversed operation stream of
/// [`controlled_increment_mod_2n`].
pub fn controlled_decrement_mod_2n_exact_reverse_increment(
    circ: &mut Circuit,
    control: &QReg,
    register: &[QReg],
    carries: &[QReg],
) {
    let width = register.len();
    if width == 0 {
        return;
    }
    if width == 1 {
        circ.cx(control, &register[0]);
        return;
    }
    assert!(carries.len() >= width - 1);
    circ.cx(control, &register[0]);
    circ.ccx(control, &register[0], &carries[0]);
    circ.cx(control, &register[0]);
    for index in 1..width - 1 {
        circ.cx(&carries[index - 1], &register[index]);
        circ.ccx(&register[index], &carries[index - 1], &carries[index]);
    }
    circ.cx(&carries[width - 2], &register[width - 1]);
    for index in (1..width - 1).rev() {
        circ.ccx(&register[index], &carries[index - 1], &carries[index]);
    }
    circ.cx(control, &register[0]);
    circ.ccx(control, &register[0], &carries[0]);
}

pub fn controlled_decrement_mod_2n_with_stream(
    circ: &mut Circuit,
    control: &QReg,
    register: &[QReg],
    scratch: &[QReg],
    stream: ReversibleDecrementStream,
) {
    match stream {
        ReversibleDecrementStream::LegacyComplementedBorrow => {
            controlled_decrement_mod_2n(circ, control, register, scratch);
        }
        ReversibleDecrementStream::ExactReverseIncrement => {
            controlled_decrement_mod_2n_exact_reverse_increment(circ, control, register, scratch);
        }
    }
}

/// Production dispatcher used by the complete register-shared EEA path.
pub fn production_controlled_decrement_mod_2n(
    circ: &mut Circuit,
    control: &QReg,
    register: &[QReg],
    scratch: &[QReg],
) {
    controlled_decrement_mod_2n_with_stream(
        circ,
        control,
        register,
        scratch,
        production_decrement_stream(),
    );
}

/// Add one when all controls are set. Higher bits are updated before lower bits.
pub fn controlled_add_one(
    circ: &mut Circuit,
    controls: &[&QReg],
    register: &[QReg],
    scratch: &[QReg],
) {
    for index in (1..register.len()).rev() {
        let mut bit_controls = Vec::with_capacity(controls.len() + index);
        bit_controls.extend_from_slice(controls);
        bit_controls.extend(register[..index].iter());
        multi_controlled_x_vchain(circ, &bit_controls, &register[index], scratch);
    }
    if let Some(low) = register.first() {
        multi_controlled_x_vchain(circ, controls, low, scratch);
    }
}

/// Subtract one when all controls are set.
///
/// Each v-chain is self-reversing, so the low-to-high traversal is the exact
/// operation-stream reverse of [`controlled_add_one`].
pub fn controlled_sub_one(
    circ: &mut Circuit,
    controls: &[&QReg],
    register: &[QReg],
    scratch: &[QReg],
) {
    if let Some(low) = register.first() {
        multi_controlled_x_vchain(circ, controls, low, scratch);
    }
    for index in 1..register.len() {
        let mut bit_controls = Vec::with_capacity(controls.len() + index);
        bit_controls.extend_from_slice(controls);
        bit_controls.extend(register[..index].iter());
        multi_controlled_x_vchain(circ, &bit_controls, &register[index], scratch);
    }
}

pub fn controlled_swap_registers(
    circ: &mut Circuit,
    control: &QReg,
    left: &[QReg],
    right: &[QReg],
) {
    assert_eq!(left.len(), right.len());
    for (left_bit, right_bit) in left.iter().zip(right) {
        circ.cswap(control, left_bit, right_bit);
    }
}

fn toggle_phase1_is_zero(circ: &mut Circuit, phase1: &QReg, target: &QReg) {
    circ.x(phase1);
    circ.cx(phase1, target);
    circ.x(phase1);
}

/// Published pre-shift block, including its length update and clean controls.
pub fn pre_shift(
    circ: &mut Circuit,
    phase1: &QReg,
    phase2: &QReg,
    work2: &[QReg],
    shift_length: &[QReg],
    scratch: &[QReg],
) {
    assert!(scratch.len() >= shift_length.len() + 4);
    let phase1_is_zero = &scratch[0];
    let both = &scratch[1];
    let carries = &scratch[2..2 + shift_length.len().saturating_sub(1)];
    toggle_phase1_is_zero(circ, phase1, phase1_is_zero);
    controlled_rotate_low(circ, phase1_is_zero, work2);
    controlled_increment_mod_2n(circ, phase1_is_zero, shift_length, carries);
    circ.ccx(phase1_is_zero, phase2, both);
    controlled_rotate_high_two(circ, both, work2);
    production_controlled_decrement_mod_2n(circ, both, shift_length, carries);
    production_controlled_decrement_mod_2n(circ, both, shift_length, carries);
    circ.ccx(phase1_is_zero, phase2, both);
    toggle_phase1_is_zero(circ, phase1, phase1_is_zero);
}

/// Exact inverse of [`pre_shift`].
pub fn pre_shift_inverse(
    circ: &mut Circuit,
    phase1: &QReg,
    phase2: &QReg,
    work2: &[QReg],
    shift_length: &[QReg],
    scratch: &[QReg],
) {
    assert!(scratch.len() >= shift_length.len() + 4);
    let phase1_is_zero = &scratch[0];
    let both = &scratch[1];
    let carries = &scratch[2..2 + shift_length.len().saturating_sub(1)];
    toggle_phase1_is_zero(circ, phase1, phase1_is_zero);
    circ.ccx(phase1_is_zero, phase2, both);
    controlled_increment_mod_2n(circ, both, shift_length, carries);
    controlled_increment_mod_2n(circ, both, shift_length, carries);
    controlled_rotate_high_two_inverse(circ, both, work2);
    circ.ccx(phase1_is_zero, phase2, both);
    production_controlled_decrement_mod_2n(circ, phase1_is_zero, shift_length, carries);
    controlled_rotate_high(circ, phase1_is_zero, work2);
    toggle_phase1_is_zero(circ, phase1, phase1_is_zero);
}

/// Published post-shift block, including its length update and clean control.
pub fn post_shift(
    circ: &mut Circuit,
    phase1: &QReg,
    phase2: &QReg,
    work2: &[QReg],
    shift_length: &[QReg],
    scratch: &[QReg],
) {
    assert!(scratch.len() >= shift_length.len() + 4);
    let both = &scratch[0];
    let carries = &scratch[1..1 + shift_length.len().saturating_sub(1)];
    controlled_rotate_low(circ, phase1, work2);
    controlled_increment_mod_2n(circ, phase1, shift_length, carries);
    circ.ccx(phase1, phase2, both);
    controlled_rotate_high_two(circ, both, work2);
    production_controlled_decrement_mod_2n(circ, both, shift_length, carries);
    production_controlled_decrement_mod_2n(circ, both, shift_length, carries);
    circ.ccx(phase1, phase2, both);
}

/// Exact inverse of [`post_shift`].
pub fn post_shift_inverse(
    circ: &mut Circuit,
    phase1: &QReg,
    phase2: &QReg,
    work2: &[QReg],
    shift_length: &[QReg],
    scratch: &[QReg],
) {
    assert!(scratch.len() >= shift_length.len() + 4);
    let both = &scratch[0];
    let carries = &scratch[1..1 + shift_length.len().saturating_sub(1)];
    circ.ccx(phase1, phase2, both);
    controlled_increment_mod_2n(circ, both, shift_length, carries);
    controlled_increment_mod_2n(circ, both, shift_length, carries);
    controlled_rotate_high_two_inverse(circ, both, work2);
    circ.ccx(phase1, phase2, both);
    production_controlled_decrement_mod_2n(circ, phase1, shift_length, carries);
    controlled_rotate_high(circ, phase1, work2);
}

/// Published phase-state update. Only the three length-register sign bits enter.
pub fn phase_update(
    circ: &mut Circuit,
    phase1: &QReg,
    phase2: &QReg,
    sign: &QReg,
    lq_sign: &QReg,
    lrp_sign: &QReg,
    ls_sign: &QReg,
    condition: &QReg,
    temporary: &QReg,
) {
    circ.x(lrp_sign);
    circ.ccx(lq_sign, lrp_sign, condition);
    circ.x(lrp_sign);
    circ.cx(sign, temporary);
    circ.cx(phase1, temporary);
    circ.ccx(condition, temporary, phase2);
    circ.cx(phase1, temporary);
    circ.cx(sign, temporary);
    circ.ccx(condition, phase2, sign);
    circ.x(lrp_sign);
    circ.ccx(lq_sign, lrp_sign, condition);
    circ.x(lrp_sign);
    circ.cx(ls_sign, phase1);
    circ.cx(ls_sign, phase2);
}

pub(crate) fn gate_counts(ops: &[Op]) -> RegisterSharedGateCounts {
    let mut counts = RegisterSharedGateCounts::default();
    for op in ops {
        match op.kind {
            OperationType::X => counts.x += 1,
            OperationType::CX => counts.cx += 1,
            OperationType::CCX => counts.ccx += 1,
            other => panic!("register-sharing shift emitted unexpected gate {other:?}"),
        }
    }
    counts.total = ops.len();
    assert_eq!(counts.total, counts.x + counts.cx + counts.ccx);
    counts
}

fn build_shift(width: usize, low: bool) -> B {
    assert!(width >= 2);
    let mut circ = Circuit::new();
    let control = circ.alloc_qreg("rs.shift.control");
    let register = circ.alloc_qreg_bits("rs.shift.work", width);
    if low {
        controlled_rotate_low(&mut circ, &control, &register);
    } else {
        controlled_rotate_high(&mut circ, &control, &register);
    }
    let builder = circ.into_builder();
    assert_eq!(builder.next_qubit as usize, width + 1);
    assert_eq!(builder.active_qubits as usize, width + 1);
    assert_eq!(builder.peak_qubits as usize, width + 1);
    builder
}

pub(crate) fn apply_scalar(ops: &[Op], mut state: u64) -> u64 {
    let bit = |word: u64, id: u64| ((word >> id) & 1) != 0;
    for op in ops {
        match op.kind {
            OperationType::X => {
                state ^= 1u64 << op.q_target.0;
            }
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
            OperationType::R | OperationType::Hmr => {
                state &= !(1u64 << op.q_target.0);
            }
            OperationType::Neg
            | OperationType::Z
            | OperationType::CZ
            | OperationType::CCZ
            | OperationType::PushCondition
            | OperationType::PopCondition => {}
            other => panic!("register-sharing scalar shift saw {other:?}"),
        }
    }
    state
}

fn rotate_low(value: u64, width: usize) -> u64 {
    let mask = (1u64 << width) - 1;
    ((value >> 1) | ((value & 1) << (width - 1))) & mask
}

fn rotate_high(value: u64, width: usize) -> u64 {
    let mask = (1u64 << width) - 1;
    ((value << 1) | (value >> (width - 1))) & mask
}

#[must_use]
pub fn exhaustive_register_shared_shift_check() -> RegisterSharedShiftProofReport {
    const REFERENCE_STEPS: usize = 1_479;
    let mut reports = Vec::new();
    let mut total_states = 0usize;
    for width in 2usize..=8 {
        let forward = build_shift(width, true);
        let reverse = build_shift(width, false);
        let states = 1usize << (width + 1);
        let mask = (1u64 << width) - 1;
        for input in 0..states as u64 {
            let control = input & 1;
            let value = (input >> 1) & mask;
            let forward_output = apply_scalar(&forward.ops, input);
            let reverse_output = apply_scalar(&reverse.ops, input);
            let expected_forward = if control == 0 {
                value
            } else {
                rotate_low(value, width)
            };
            let expected_reverse = if control == 0 {
                value
            } else {
                rotate_high(value, width)
            };
            assert_eq!(forward_output & 1, control);
            assert_eq!(reverse_output & 1, control);
            assert_eq!((forward_output >> 1) & mask, expected_forward);
            assert_eq!((reverse_output >> 1) & mask, expected_reverse);
            assert_eq!(apply_scalar(&reverse.ops, forward_output), input);
            assert_eq!(apply_scalar(&forward.ops, reverse_output), input);
        }
        let forward_counts = gate_counts(&forward.ops);
        let reverse_counts = gate_counts(&reverse.ops);
        let expected = RegisterSharedGateCounts {
            x: 0,
            cx: 2 * (width - 1),
            ccx: width - 1,
            total: 3 * (width - 1),
        };
        assert_eq!(forward_counts, expected);
        assert_eq!(reverse_counts, expected);
        total_states += states;
        reports.push(RegisterSharedShiftWidthReport {
            width,
            basis_states_checked: states,
            forward: forward_counts,
            reverse: reverse_counts,
        });
    }
    let work259_shift_toffoli = WORK_BITS - 1;
    let legacy_schematic_four_rotations_toffoli = 4 * work259_shift_toffoli;
    RegisterSharedShiftProofReport {
        widths_checked: reports.len(),
        basis_states_checked: total_states,
        allocation_free_streams_checked: 2 * reports.len(),
        phase_clean_streams_checked: 2 * reports.len(),
        widths: reports,
        work259_shift_toffoli,
        legacy_schematic_four_rotations_toffoli,
        reference_steps: REFERENCE_STEPS,
        legacy_schematic_four_rotations_total_toffoli: REFERENCE_STEPS
            * legacy_schematic_four_rotations_toffoli,
    }
}

#[derive(Clone, Copy)]
enum LengthPrimitive {
    Increment,
    Decrement,
    ExactReverseDecrement,
    ControlledIncrement,
    ControlledDecrement,
    ExactReverseControlledDecrement,
    TwoControlAddOne,
    TwoControlSubOne,
}

fn build_length_primitive(width: usize, primitive: LengthPrimitive) -> B {
    assert!(width > 0);
    let mut circ = Circuit::new();
    match primitive {
        LengthPrimitive::Increment
        | LengthPrimitive::Decrement
        | LengthPrimitive::ExactReverseDecrement => {
            let register = circ.alloc_qreg_bits("rs.length", width);
            let scratch = circ.alloc_qreg_bits("rs.length.scratch", width.saturating_sub(1));
            match primitive {
                LengthPrimitive::Increment => increment_mod_2n(&mut circ, &register, &scratch),
                LengthPrimitive::Decrement => decrement_mod_2n(&mut circ, &register, &scratch),
                LengthPrimitive::ExactReverseDecrement => {
                    decrement_mod_2n_exact_reverse_increment(&mut circ, &register, &scratch);
                }
                _ => unreachable!(),
            }
            circ.into_builder()
        }
        LengthPrimitive::ControlledIncrement
        | LengthPrimitive::ControlledDecrement
        | LengthPrimitive::ExactReverseControlledDecrement => {
            let control = circ.alloc_qreg("rs.length.control");
            let register = circ.alloc_qreg_bits("rs.length", width);
            let scratch = circ.alloc_qreg_bits("rs.length.scratch", width.saturating_sub(1));
            match primitive {
                LengthPrimitive::ControlledIncrement => {
                    controlled_increment_mod_2n(&mut circ, &control, &register, &scratch);
                }
                LengthPrimitive::ControlledDecrement => {
                    controlled_decrement_mod_2n(&mut circ, &control, &register, &scratch);
                }
                LengthPrimitive::ExactReverseControlledDecrement => {
                    controlled_decrement_mod_2n_exact_reverse_increment(
                        &mut circ, &control, &register, &scratch,
                    );
                }
                _ => unreachable!(),
            }
            circ.into_builder()
        }
        LengthPrimitive::TwoControlAddOne | LengthPrimitive::TwoControlSubOne => {
            let controls = circ.alloc_qreg_bits("rs.length.controls", 2);
            let register = circ.alloc_qreg_bits("rs.length", width);
            let scratch = circ.alloc_qreg_bits("rs.length.scratch", width.saturating_sub(1));
            let control_refs = [&controls[0], &controls[1]];
            match primitive {
                LengthPrimitive::TwoControlAddOne => {
                    controlled_add_one(&mut circ, &control_refs, &register, &scratch);
                }
                LengthPrimitive::TwoControlSubOne => {
                    controlled_sub_one(&mut circ, &control_refs, &register, &scratch);
                }
                _ => unreachable!(),
            }
            circ.into_builder()
        }
    }
}

#[must_use]
pub fn exhaustive_register_shared_length_check() -> RegisterSharedLengthProofReport {
    let mut basis_states_checked = 0usize;
    let mut scratch_clean_checks = 0usize;
    let mut inverse_pair_checks = 0usize;
    for width in 1..=8 {
        let mask = (1u64 << width) - 1;
        let increment = build_length_primitive(width, LengthPrimitive::Increment);
        let decrement = build_length_primitive(width, LengthPrimitive::Decrement);
        for value in 0..=mask {
            let incremented = apply_scalar(&increment.ops, value);
            let decremented = apply_scalar(&decrement.ops, value);
            assert_eq!(incremented & mask, value.wrapping_add(1) & mask);
            assert_eq!(decremented & mask, value.wrapping_sub(1) & mask);
            assert_eq!(incremented >> width, 0);
            assert_eq!(decremented >> width, 0);
            assert_eq!(apply_scalar(&decrement.ops, incremented), value);
            assert_eq!(apply_scalar(&increment.ops, decremented), value);
            basis_states_checked += 2;
            scratch_clean_checks += 2;
            inverse_pair_checks += 2;
        }

        let controlled_increment =
            build_length_primitive(width, LengthPrimitive::ControlledIncrement);
        let controlled_decrement =
            build_length_primitive(width, LengthPrimitive::ControlledDecrement);
        for input in 0..(1u64 << (width + 1)) {
            let control = input & 1;
            let value = (input >> 1) & mask;
            let expected_increment = if control == 0 {
                value
            } else {
                value.wrapping_add(1) & mask
            };
            let expected_decrement = if control == 0 {
                value
            } else {
                value.wrapping_sub(1) & mask
            };
            let incremented = apply_scalar(&controlled_increment.ops, input);
            let decremented = apply_scalar(&controlled_decrement.ops, input);
            assert_eq!(incremented & 1, control);
            assert_eq!(decremented & 1, control);
            assert_eq!((incremented >> 1) & mask, expected_increment);
            assert_eq!((decremented >> 1) & mask, expected_decrement);
            assert_eq!(incremented >> (width + 1), 0);
            assert_eq!(decremented >> (width + 1), 0);
            assert_eq!(apply_scalar(&controlled_decrement.ops, incremented), input);
            assert_eq!(apply_scalar(&controlled_increment.ops, decremented), input);
            basis_states_checked += 2;
            scratch_clean_checks += 2;
            inverse_pair_checks += 2;
        }

        let add_one = build_length_primitive(width, LengthPrimitive::TwoControlAddOne);
        let sub_one = build_length_primitive(width, LengthPrimitive::TwoControlSubOne);
        for input in 0..(1u64 << (width + 2)) {
            let controls = input & 3;
            let value = (input >> 2) & mask;
            let active = controls == 3;
            let expected_add = if active {
                value.wrapping_add(1) & mask
            } else {
                value
            };
            let expected_sub = if active {
                value.wrapping_sub(1) & mask
            } else {
                value
            };
            let added = apply_scalar(&add_one.ops, input);
            let subtracted = apply_scalar(&sub_one.ops, input);
            assert_eq!(added & 3, controls);
            assert_eq!(subtracted & 3, controls);
            assert_eq!((added >> 2) & mask, expected_add);
            assert_eq!((subtracted >> 2) & mask, expected_sub);
            assert_eq!(added >> (width + 2), 0);
            assert_eq!(subtracted >> (width + 2), 0);
            assert_eq!(apply_scalar(&sub_one.ops, added), input);
            assert_eq!(apply_scalar(&add_one.ops, subtracted), input);
            basis_states_checked += 2;
            scratch_clean_checks += 2;
            inverse_pair_checks += 2;
        }
    }

    RegisterSharedLengthProofReport {
        widths_checked: 8,
        basis_states_checked,
        scratch_clean_checks,
        inverse_pair_checks,
        increment9: gate_counts(&build_length_primitive(9, LengthPrimitive::Increment).ops),
        decrement9: gate_counts(&build_length_primitive(9, LengthPrimitive::Decrement).ops),
        controlled_increment9: gate_counts(
            &build_length_primitive(9, LengthPrimitive::ControlledIncrement).ops,
        ),
        controlled_decrement9: gate_counts(
            &build_length_primitive(9, LengthPrimitive::ControlledDecrement).ops,
        ),
        two_control_add_one9: gate_counts(
            &build_length_primitive(9, LengthPrimitive::TwoControlAddOne).ops,
        ),
        two_control_sub_one9: gate_counts(
            &build_length_primitive(9, LengthPrimitive::TwoControlSubOne).ops,
        ),
    }
}

fn assert_same_length_allocation_profile(left: &B, right: &B) {
    assert_eq!(left.next_qubit, right.next_qubit);
    assert_eq!(left.next_bit, right.next_bit);
    assert_eq!(left.active_qubits, right.active_qubits);
    assert_eq!(left.peak_qubits, right.peak_qubits);
    assert_eq!(left.free_qubits, right.free_qubits);
    assert_eq!(left.allocation_serial, right.allocation_serial);
}

fn assert_phase_neutral_classical_stream(ops: &[Op]) {
    use crate::circuit::NO_BIT;

    for operation in ops {
        assert!(
            matches!(
                operation.kind,
                OperationType::X | OperationType::CX | OperationType::CCX
            ),
            "reverse-decrement proof saw phase-capable operation {:?}",
            operation.kind
        );
        assert_eq!(operation.c_condition, NO_BIT);
    }
}

/// Exhaustively prove the exact reverse-decrement streams through width 16.
///
/// The operation-stream equalities are stronger than scalar equivalence: every
/// exact decrement operation must equal the corresponding increment operation
/// visited in reverse order. The scalar sweep separately checks arithmetic,
/// inverse behavior, controls, and clean scratch lanes for every basis input.
#[must_use]
pub fn exhaustive_exact_reverse_decrement_check() -> RegisterSharedReverseDecrementProofReport {
    const MAX_WIDTH: usize = 16;

    assert!(!REGISTER_SHARED_REVERSE_DECREMENT_STREAM_SOURCE_BAKE);
    let mut widths = Vec::with_capacity(MAX_WIDTH);
    let mut direct_basis_states_checked = 0usize;
    let mut controlled_basis_states_checked = 0usize;
    let mut exact_reverse_stream_checks = 0usize;
    let mut scalar_forward_checks = 0usize;
    let mut inverse_pair_checks = 0usize;
    let mut phase_clean_streams_checked = 0usize;
    let mut ancilla_clean_checks = 0usize;
    let mut allocation_profile_checks = 0usize;
    let mut toffoli_preservation_checks = 0usize;
    let mut local_legacy_ops = 0usize;
    let mut local_exact_reverse_ops = 0usize;

    for width in 1..=MAX_WIDTH {
        let increment = build_length_primitive(width, LengthPrimitive::Increment);
        let legacy_decrement = build_length_primitive(width, LengthPrimitive::Decrement);
        let exact_reverse_decrement =
            build_length_primitive(width, LengthPrimitive::ExactReverseDecrement);
        assert!(exact_reverse_decrement
            .ops
            .iter()
            .eq(increment.ops.iter().rev()));
        exact_reverse_stream_checks += 1;

        let controlled_increment =
            build_length_primitive(width, LengthPrimitive::ControlledIncrement);
        let legacy_controlled_decrement =
            build_length_primitive(width, LengthPrimitive::ControlledDecrement);
        let exact_reverse_controlled_decrement =
            build_length_primitive(width, LengthPrimitive::ExactReverseControlledDecrement);
        assert!(exact_reverse_controlled_decrement
            .ops
            .iter()
            .eq(controlled_increment.ops.iter().rev()));
        exact_reverse_stream_checks += 1;

        let two_control_add = build_length_primitive(width, LengthPrimitive::TwoControlAddOne);
        let two_control_sub = build_length_primitive(width, LengthPrimitive::TwoControlSubOne);
        assert!(two_control_sub
            .ops
            .iter()
            .eq(two_control_add.ops.iter().rev()));
        exact_reverse_stream_checks += 1;

        for builder in [
            &increment,
            &legacy_decrement,
            &exact_reverse_decrement,
            &controlled_increment,
            &legacy_controlled_decrement,
            &exact_reverse_controlled_decrement,
            &two_control_add,
            &two_control_sub,
        ] {
            assert_phase_neutral_classical_stream(&builder.ops);
            phase_clean_streams_checked += 1;
        }

        assert_same_length_allocation_profile(&legacy_decrement, &exact_reverse_decrement);
        assert_same_length_allocation_profile(
            &legacy_controlled_decrement,
            &exact_reverse_controlled_decrement,
        );
        allocation_profile_checks += 2;

        let increment_counts = gate_counts(&increment.ops);
        let legacy_decrement_counts = gate_counts(&legacy_decrement.ops);
        let exact_reverse_decrement_counts = gate_counts(&exact_reverse_decrement.ops);
        let controlled_increment_counts = gate_counts(&controlled_increment.ops);
        let legacy_controlled_decrement_counts = gate_counts(&legacy_controlled_decrement.ops);
        let exact_reverse_controlled_decrement_counts =
            gate_counts(&exact_reverse_controlled_decrement.ops);

        assert_eq!(exact_reverse_decrement_counts, increment_counts);
        assert_eq!(
            legacy_decrement_counts.cx,
            exact_reverse_decrement_counts.cx
        );
        assert_eq!(
            legacy_decrement_counts.ccx,
            exact_reverse_decrement_counts.ccx
        );
        assert_eq!(
            legacy_decrement_counts
                .total
                .checked_sub(exact_reverse_decrement_counts.total)
                .expect("exact reverse decrement grew the direct stream"),
            4 * width.saturating_sub(1)
        );
        toffoli_preservation_checks += 1;

        assert_eq!(
            exact_reverse_controlled_decrement_counts,
            controlled_increment_counts
        );
        assert_eq!(
            legacy_controlled_decrement_counts.cx,
            exact_reverse_controlled_decrement_counts.cx
        );
        assert_eq!(
            legacy_controlled_decrement_counts.ccx,
            exact_reverse_controlled_decrement_counts.ccx
        );
        assert_eq!(
            legacy_controlled_decrement_counts
                .total
                .checked_sub(exact_reverse_controlled_decrement_counts.total)
                .expect("exact reverse decrement grew the controlled stream"),
            4 * width.saturating_sub(1)
        );
        toffoli_preservation_checks += 1;

        local_legacy_ops +=
            legacy_decrement_counts.total + legacy_controlled_decrement_counts.total;
        local_exact_reverse_ops +=
            exact_reverse_decrement_counts.total + exact_reverse_controlled_decrement_counts.total;

        let mask = (1u64 << width) - 1;
        let direct_basis_states = 1usize << width;
        for value in 0..direct_basis_states as u64 {
            let incremented = apply_scalar(&increment.ops, value);
            let decremented = apply_scalar(&exact_reverse_decrement.ops, value);
            assert_eq!(incremented & mask, value.wrapping_add(1) & mask);
            assert_eq!(decremented & mask, value.wrapping_sub(1) & mask);
            assert_eq!(incremented >> width, 0);
            assert_eq!(decremented >> width, 0);
            assert_eq!(
                apply_scalar(&exact_reverse_decrement.ops, incremented),
                value
            );
            assert_eq!(apply_scalar(&increment.ops, decremented), value);
            scalar_forward_checks += 2;
            inverse_pair_checks += 2;
            ancilla_clean_checks += 2;
        }
        direct_basis_states_checked += direct_basis_states;

        let controlled_basis_states = 1usize << (width + 1);
        for input in 0..controlled_basis_states as u64 {
            let control = input & 1;
            let value = (input >> 1) & mask;
            let expected_increment = if control == 0 {
                value
            } else {
                value.wrapping_add(1) & mask
            };
            let expected_decrement = if control == 0 {
                value
            } else {
                value.wrapping_sub(1) & mask
            };
            let incremented = apply_scalar(&controlled_increment.ops, input);
            let decremented = apply_scalar(&exact_reverse_controlled_decrement.ops, input);
            assert_eq!(incremented & 1, control);
            assert_eq!(decremented & 1, control);
            assert_eq!((incremented >> 1) & mask, expected_increment);
            assert_eq!((decremented >> 1) & mask, expected_decrement);
            assert_eq!(incremented >> (width + 1), 0);
            assert_eq!(decremented >> (width + 1), 0);
            assert_eq!(
                apply_scalar(&exact_reverse_controlled_decrement.ops, incremented),
                input
            );
            assert_eq!(apply_scalar(&controlled_increment.ops, decremented), input);
            scalar_forward_checks += 2;
            inverse_pair_checks += 2;
            ancilla_clean_checks += 2;
        }
        controlled_basis_states_checked += controlled_basis_states;

        widths.push(RegisterSharedReverseDecrementWidthReport {
            width,
            direct_basis_states,
            controlled_basis_states,
            increment: increment_counts,
            legacy_decrement: legacy_decrement_counts,
            exact_reverse_decrement: exact_reverse_decrement_counts,
            controlled_increment: controlled_increment_counts,
            legacy_controlled_decrement: legacy_controlled_decrement_counts,
            exact_reverse_controlled_decrement: exact_reverse_controlled_decrement_counts,
        });
    }

    RegisterSharedReverseDecrementProofReport {
        widths_checked: widths.len(),
        direct_basis_states_checked,
        controlled_basis_states_checked,
        exact_reverse_stream_checks,
        scalar_forward_checks,
        inverse_pair_checks,
        phase_clean_streams_checked,
        ancilla_clean_checks,
        allocation_profile_checks,
        toffoli_preservation_checks,
        local_legacy_ops,
        local_exact_reverse_ops,
        local_ops_removed: local_legacy_ops
            .checked_sub(local_exact_reverse_ops)
            .expect("exact reverse decrement grew the local proof streams"),
        widths,
    }
}

fn build_work_swap(width: usize) -> B {
    let mut circ = Circuit::new();
    let control = circ.alloc_qreg("rs.swap.control");
    let left = circ.alloc_qreg_bits("rs.swap.left", width);
    let right = circ.alloc_qreg_bits("rs.swap.right", width);
    controlled_swap_registers(&mut circ, &control, &left, &right);
    circ.into_builder()
}

fn build_phase_update() -> B {
    let mut circ = Circuit::new();
    let phase1 = circ.alloc_qreg("rs.phase1");
    let phase2 = circ.alloc_qreg("rs.phase2");
    let sign = circ.alloc_qreg("rs.sign");
    let lq_sign = circ.alloc_qreg("rs.lq.sign");
    let lrp_sign = circ.alloc_qreg("rs.lrp.sign");
    let ls_sign = circ.alloc_qreg("rs.ls.sign");
    let condition = circ.alloc_qreg("rs.phase.condition");
    let temporary = circ.alloc_qreg("rs.phase.temporary");
    phase_update(
        &mut circ, &phase1, &phase2, &sign, &lq_sign, &lrp_sign, &ls_sign, &condition, &temporary,
    );
    circ.into_builder()
}

#[must_use]
pub fn exhaustive_register_shared_control_check() -> RegisterSharedControlProofReport {
    let mut work_swap_basis_states_checked = 0usize;
    let mut work_swap_inverse_checks = 0usize;
    for width in 1..=8 {
        let swap = build_work_swap(width);
        let mask = (1u64 << width) - 1;
        for input in 0..(1u64 << (2 * width + 1)) {
            let control = input & 1;
            let left = (input >> 1) & mask;
            let right = (input >> (width + 1)) & mask;
            let output = apply_scalar(&swap.ops, input);
            let (expected_left, expected_right) = if control == 0 {
                (left, right)
            } else {
                (right, left)
            };
            assert_eq!(output & 1, control);
            assert_eq!((output >> 1) & mask, expected_left);
            assert_eq!((output >> (width + 1)) & mask, expected_right);
            assert_eq!(apply_scalar(&swap.ops, output), input);
            work_swap_basis_states_checked += 1;
            work_swap_inverse_checks += 1;
        }
    }

    let phase = build_phase_update();
    let phase_update_basis_states_checked = 1usize << 6;
    for input in 0..phase_update_basis_states_checked as u64 {
        let mut phase1 = (input & 1) != 0;
        let mut phase2 = ((input >> 1) & 1) != 0;
        let mut sign = ((input >> 2) & 1) != 0;
        let lq_sign = ((input >> 3) & 1) != 0;
        let lrp_sign = ((input >> 4) & 1) != 0;
        let ls_sign = ((input >> 5) & 1) != 0;
        let condition = lq_sign && !lrp_sign;
        let temporary = sign ^ phase1;
        phase2 ^= condition && temporary;
        sign ^= condition && phase2;
        phase1 ^= ls_sign;
        phase2 ^= ls_sign;
        let expected = u64::from(phase1)
            | (u64::from(phase2) << 1)
            | (u64::from(sign) << 2)
            | (u64::from(lq_sign) << 3)
            | (u64::from(lrp_sign) << 4)
            | (u64::from(ls_sign) << 5);
        let output = apply_scalar(&phase.ops, input);
        assert_eq!(output, expected);
        assert_eq!(output >> 6, 0);
    }

    RegisterSharedControlProofReport {
        work_swap_widths_checked: 8,
        work_swap_basis_states_checked,
        work_swap_inverse_checks,
        work259_swap: gate_counts(&build_work_swap(WORK_BITS).ops),
        phase_update_basis_states_checked,
        phase_update_scratch_clean_checks: phase_update_basis_states_checked,
        phase_update: gate_counts(&phase.ops),
    }
}

#[derive(Clone, Copy)]
enum PrePostPrimitive {
    Pre,
    PreInverse,
    Post,
    PostInverse,
}

fn build_pre_post_shift(work_width: usize, length_width: usize, primitive: PrePostPrimitive) -> B {
    assert!(work_width >= 3);
    assert!(length_width > 0);
    let mut circ = Circuit::new();
    let phase1 = circ.alloc_qreg("rs.prepost.phase1");
    let phase2 = circ.alloc_qreg("rs.prepost.phase2");
    let work2 = circ.alloc_qreg_bits("rs.prepost.work2", work_width);
    let shift_length = circ.alloc_qreg_bits("rs.prepost.length", length_width);
    let scratch = circ.alloc_qreg_bits("rs.prepost.scratch", length_width + 4);
    match primitive {
        PrePostPrimitive::Pre => {
            pre_shift(&mut circ, &phase1, &phase2, &work2, &shift_length, &scratch)
        }
        PrePostPrimitive::PreInverse => {
            pre_shift_inverse(&mut circ, &phase1, &phase2, &work2, &shift_length, &scratch)
        }
        PrePostPrimitive::Post => {
            post_shift(&mut circ, &phase1, &phase2, &work2, &shift_length, &scratch)
        }
        PrePostPrimitive::PostInverse => {
            post_shift_inverse(&mut circ, &phase1, &phase2, &work2, &shift_length, &scratch)
        }
    }
    circ.into_builder()
}

fn expected_pre_shift(
    phase1: bool,
    phase2: bool,
    mut work: u64,
    work_width: usize,
    mut shift_length: u64,
    length_mask: u64,
) -> (u64, u64) {
    if !phase1 {
        work = rotate_low(work, work_width);
        shift_length = shift_length.wrapping_add(1) & length_mask;
        if phase2 {
            work = rotate_high(rotate_high(work, work_width), work_width);
            shift_length = shift_length.wrapping_sub(2) & length_mask;
        }
    }
    (work, shift_length)
}

fn expected_post_shift(
    phase1: bool,
    phase2: bool,
    mut work: u64,
    work_width: usize,
    mut shift_length: u64,
    length_mask: u64,
) -> (u64, u64) {
    if phase1 {
        work = rotate_low(work, work_width);
        shift_length = shift_length.wrapping_add(1) & length_mask;
        if phase2 {
            work = rotate_high(rotate_high(work, work_width), work_width);
            shift_length = shift_length.wrapping_sub(2) & length_mask;
        }
    }
    (work, shift_length)
}

#[must_use]
pub fn exhaustive_register_shared_pre_post_check() -> RegisterSharedPrePostProofReport {
    const TEST_LENGTH_WIDTH: usize = 4;
    const REFERENCE_LENGTH_WIDTH: usize = 9;
    const REFERENCE_STEPS: usize = 1_479;
    let mut basis_states_checked = 0usize;
    let mut scratch_clean_checks = 0usize;
    let mut inverse_pair_checks = 0usize;
    for work_width in 3..=8 {
        let pre = build_pre_post_shift(work_width, TEST_LENGTH_WIDTH, PrePostPrimitive::Pre);
        let pre_inverse =
            build_pre_post_shift(work_width, TEST_LENGTH_WIDTH, PrePostPrimitive::PreInverse);
        let post = build_pre_post_shift(work_width, TEST_LENGTH_WIDTH, PrePostPrimitive::Post);
        let post_inverse =
            build_pre_post_shift(work_width, TEST_LENGTH_WIDTH, PrePostPrimitive::PostInverse);
        let work_mask = (1u64 << work_width) - 1;
        let length_mask = (1u64 << TEST_LENGTH_WIDTH) - 1;
        let data_width = 2 + work_width + TEST_LENGTH_WIDTH;
        for input in 0..(1u64 << data_width) {
            let phase1 = (input & 1) != 0;
            let phase2 = ((input >> 1) & 1) != 0;
            let work = (input >> 2) & work_mask;
            let shift_length = (input >> (2 + work_width)) & length_mask;
            let (pre_work, pre_length) =
                expected_pre_shift(phase1, phase2, work, work_width, shift_length, length_mask);
            let (post_work, post_length) =
                expected_post_shift(phase1, phase2, work, work_width, shift_length, length_mask);
            let expected_pre = (input & 3) | (pre_work << 2) | (pre_length << (2 + work_width));
            let expected_post = (input & 3) | (post_work << 2) | (post_length << (2 + work_width));
            let pre_output = apply_scalar(&pre.ops, input);
            let post_output = apply_scalar(&post.ops, input);
            assert_eq!(pre_output, expected_pre);
            assert_eq!(post_output, expected_post);
            assert_eq!(pre_output >> data_width, 0);
            assert_eq!(post_output >> data_width, 0);
            assert_eq!(apply_scalar(&pre_inverse.ops, pre_output), input);
            assert_eq!(apply_scalar(&post_inverse.ops, post_output), input);
            basis_states_checked += 2;
            scratch_clean_checks += 2;
            inverse_pair_checks += 2;
        }
    }

    let pre_shift259_length9 = gate_counts(
        &build_pre_post_shift(WORK_BITS, REFERENCE_LENGTH_WIDTH, PrePostPrimitive::Pre).ops,
    );
    let pre_shift259_length9_inverse = gate_counts(
        &build_pre_post_shift(
            WORK_BITS,
            REFERENCE_LENGTH_WIDTH,
            PrePostPrimitive::PreInverse,
        )
        .ops,
    );
    let post_shift259_length9 = gate_counts(
        &build_pre_post_shift(WORK_BITS, REFERENCE_LENGTH_WIDTH, PrePostPrimitive::Post).ops,
    );
    let post_shift259_length9_inverse = gate_counts(
        &build_pre_post_shift(
            WORK_BITS,
            REFERENCE_LENGTH_WIDTH,
            PrePostPrimitive::PostInverse,
        )
        .ops,
    );
    RegisterSharedPrePostProofReport {
        work_widths_checked: 6,
        basis_states_checked,
        scratch_clean_checks,
        inverse_pair_checks,
        pre_shift259_length9,
        pre_shift259_length9_inverse,
        post_shift259_length9,
        post_shift259_length9_inverse,
        emitted_rotations_per_step: 6,
        reference_steps: REFERENCE_STEPS,
        total_pre_post_toffoli: REFERENCE_STEPS
            * (pre_shift259_length9.ccx + post_shift259_length9.ccx),
    }
}

fn build_variable_rotation(width: usize, amount_width: usize, high: bool) -> B {
    let mut circ = Circuit::new();
    let amount = circ.alloc_qreg_bits("rs.barrel.amount", amount_width);
    let register = circ.alloc_qreg_bits("rs.barrel.work", width);
    if high {
        variable_rotate_high(&mut circ, &amount, &register);
    } else {
        variable_rotate_low(&mut circ, &amount, &register);
    }
    circ.into_builder()
}

fn rotate_high_by_value(mut value: u64, width: usize, offset: usize) -> u64 {
    for _ in 0..offset % width {
        value = rotate_high(value, width);
    }
    value
}

#[must_use]
pub fn exhaustive_register_shared_barrel_check() -> RegisterSharedBarrelProofReport {
    let mut basis_states_checked = 0usize;
    let mut inverse_pair_checks = 0usize;
    for width in 2usize..=8 {
        let amount_width = (usize::BITS - (width - 1).leading_zeros()) as usize;
        let high = build_variable_rotation(width, amount_width, true);
        let low = build_variable_rotation(width, amount_width, false);
        let value_mask = (1u64 << width) - 1;
        let amount_mask = (1u64 << amount_width) - 1;
        for input in 0..(1u64 << (width + amount_width)) {
            let amount = (input & amount_mask) as usize;
            let value = (input >> amount_width) & value_mask;
            let expected = (input & amount_mask)
                | (rotate_high_by_value(value, width, amount) << amount_width);
            let output = apply_scalar(&high.ops, input);
            assert_eq!(output, expected);
            assert_eq!(apply_scalar(&low.ops, output), input);
            basis_states_checked += 1;
            inverse_pair_checks += 1;
        }
    }
    let work259_amount9_high = gate_counts(&build_variable_rotation(WORK_BITS, 9, true).ops);
    let work259_amount9_low = gate_counts(&build_variable_rotation(WORK_BITS, 9, false).ops);
    assert_eq!(work259_amount9_high, work259_amount9_low);
    RegisterSharedBarrelProofReport {
        widths_checked: 7,
        basis_states_checked,
        inverse_pair_checks,
        work259_amount9_high,
        work259_amount9_low,
        toffoli_per_direction: work259_amount9_high.ccx,
    }
}

fn allocation_skeleton(length_widths: &[usize], scratch_widths: &[usize]) -> B {
    let mut circ = Circuit::new();
    let mut dx = circ.alloc_qreg_bits("rs.dx", FIELD_PASSENGER_BITS);
    let dy = circ.alloc_qreg_bits("rs.passenger", FIELD_PASSENGER_BITS);
    let work2_padding = WORK_BITS - dx.len();
    dx.extend(circ.alloc_qreg_bits("rs.work2.pad", work2_padding));
    let work1 = circ.alloc_qreg_bits("rs.work1", WORK_BITS);
    let lengths: Vec<Vec<QReg>> = length_widths
        .iter()
        .enumerate()
        .map(|(index, &width)| circ.alloc_qreg_bits(&format!("rs.length.{index}"), width))
        .collect();
    let controls = circ.alloc_qreg_bits("rs.controls", CONTROL_BITS);
    let scratch: Vec<Vec<QReg>> = scratch_widths
        .iter()
        .enumerate()
        .map(|(index, &width)| circ.alloc_qreg_bits(&format!("rs.scratch.{index}"), width))
        .collect();

    for lane in dx
        .into_iter()
        .chain(dy)
        .chain(work1)
        .chain(lengths.into_iter().flatten())
        .chain(controls)
        .chain(scratch.into_iter().flatten())
    {
        circ.zero_and_free(lane);
    }
    circ.flush_pending_frees();
    circ.into_builder()
}

fn assert_reset_only_skeleton(builder: &B, expected_peak: usize) {
    assert_eq!(builder.peak_qubits as usize, expected_peak);
    assert_eq!(builder.active_qubits, 0);
    assert_eq!(builder.ops.len(), expected_peak);
    assert!(builder
        .ops
        .iter()
        .all(|operation| operation.kind == OperationType::R));
}

#[must_use]
pub fn register_shared_allocation_skeleton_check() -> RegisterSharedAllocationReport {
    let paper_core = allocation_skeleton(&LENGTH_BITS, &[]);
    assert_reset_only_skeleton(&paper_core, PROJECTED_INVERSION_PEAK);
    let reference_port = allocation_skeleton(&REFERENCE_LENGTH_BITS, &REFERENCE_SCRATCH_POOLS);
    assert_reset_only_skeleton(&reference_port, REFERENCE_PORT_PEAK);
    RegisterSharedAllocationReport {
        paper_core_peak_qubits: paper_core.peak_qubits as usize,
        paper_core_final_active_qubits: paper_core.active_qubits as usize,
        paper_core_reset_operations: paper_core.ops.len(),
        reference_port_peak_qubits: reference_port.peak_qubits as usize,
        reference_port_final_active_qubits: reference_port.active_qubits as usize,
        reference_port_reset_operations: reference_port.ops.len(),
        reference_port_scratch_bits: REFERENCE_SCRATCH_POOLS.iter().sum(),
        input_dx_bits: FIELD_PASSENGER_BITS,
        passenger_bits: FIELD_PASSENGER_BITS,
        work1_bits: WORK_BITS,
        work2_bits: WORK_BITS,
        paper_length_bits: LENGTH_BITS.iter().sum(),
        reference_length_bits: REFERENCE_LENGTH_BITS.iter().sum(),
        control_bits: CONTROL_BITS,
    }
}
