//! Reversible gate port of the public register-sharing EEA reference.
//!
//! This module is intentionally staged. The arithmetic primitives and their
//! reduced-width exhaustive proofs land before the complete 1,479-step circuit.

use crate::point_add::trailmix_port::circuit::{Circuit, QReg};
use std::sync::atomic::{AtomicUsize, Ordering};

use super::register_shared_eea_microkernels::{
    apply_scalar, controlled_add_one, controlled_increment_mod_2n, controlled_sub_one,
    controlled_swap_registers, gate_counts, increment_mod_2n, multi_controlled_x_vchain,
    production_controlled_decrement_mod_2n as controlled_decrement_mod_2n,
    production_decrement_mod_2n as decrement_mod_2n, variable_rotate_high,
    variable_rotate_high_refs, variable_rotate_low, variable_rotate_low_refs,
    RegisterSharedGateCounts,
};
use super::shrunken_pz_state_machine::with_bit_length_callsite;
use crate::point_add::trailmix_port::arith::mcx::mcx_dirty_ladder;
use crate::point_add::B;

pub const REFERENCE_LENGTH_WIDTH: usize = 9;
pub const REFERENCE_R_LENGTH_WIDTH: usize = 8;
pub const SUB820_L_Q_WIDTH: usize = 8;
pub const REFERENCE_STEPS: usize = 1_479;
pub const COEFFICIENT_RAW_BITLEN_LOAN_FLAG: &str = "LOWQ_REUSE_COEFFICIENT_RAW_BITLEN_LOAN";
pub const INPLACE_ROTATED_BITLEN_BOUNDARY_FLAG: &str = "LOWQ_INPLACE_ROTATED_BITLEN_BOUNDARY";
pub const FUSED_PREFIX_SCRATCH_LOAN_FLAG: &str = "LOWQ_LOAN_FUSED_PREFIX_SCRATCH";
pub const PROMISED_LQ_SWAP_BORROW_FLAG: &str = "LOWQ_REUSE_LQ_AS_SWAP_OLD_R_LENGTH";
pub const SPLIT_COEFFICIENT_ROTATION_LIFETIME_FLAG: &str =
    "LOWQ_SPLIT_COEFFICIENT_ROTATION_LIFETIME";
pub const COEFFICIENT_LESS_THAN_LANE_REUSE_FLAG: &str = "LOWQ_REUSE_COEFFICIENT_LESS_THAN_LANES";
pub const CLEAN_CHAIN_COEFFICIENT_ADD_LENDER_FLAG: &str =
    "LOWQ_REUSE_CLEAN_CHAIN_FOR_COEFFICIENT_ADD";
pub const PRESERVED_DY_TOP_PREFIX_LOAN_FLAG: &str = "LOWQ_REUSE_PRESERVED_DY_TOP_FOR_PREFIX";
pub const MIXED_WIDTH_L_R_PRIME_FLAG: &str = "LOWQ_MIXED_WIDTH_L_R_PRIME";
pub const PAIRED_BITLEN_SOURCE_COMPLEMENT_FLAG: &str =
    "LOWQ_PAIRED_BITLEN_SOURCE_COMPLEMENT";
pub const COEFFICIENT_NONNEGATIVE_X_CANCEL_FLAG: &str =
    "LOWQ_COEFFICIENT_NONNEGATIVE_X_CANCEL";
pub const Q845_LIFETIME_COEFFICIENT_FUSION_FLAG: &str =
    "LOWQ_Q845_LIFETIME_COEFFICIENT_FUSION";
pub const PROMISED_SWAP_SUPPORT_LIFETIME_FUSION_FLAG: &str =
    "LOWQ_FUSE_PROMISED_SWAP_SUPPORT_LIFETIME";
pub const Q845_SWAP_ONLY_T_PRIME_LENGTH_FLAG: &str =
    "LOWQ_Q845_SWAP_ONLY_T_PRIME_LENGTH";
pub const Q851_TRUNCATED_SWAP_ONLY_GUARD_FLAG: &str =
    "LOWQ_Q851_TRUNCATED_SWAP_ONLY_GUARD";
pub const Q851_FIXED_SIGN_EVENT_FLAG: &str = "LOWQ_Q851_FIXED_SIGN_EVENT";
pub const Q830_DIRTY_FIXED_SIGN_EVENT_FLAG: &str =
    "LOWQ_Q830_DIRTY_FIXED_SIGN_EVENT";
pub const Q830_DIRECT_SWAP_METADATA_FLAG: &str =
    "LOWQ_Q830_DIRECT_SWAP_METADATA";
pub const Q830_COEFFICIENT_COUNTER_RELOCATION_FLAG: &str =
    "LOWQ_Q830_COEFFICIENT_COUNTER_RELOCATION";
pub const SUB800_INPLACE_GUARD_ADDRESS_FLAG: &str =
    "LOWQ_SUB800_INPLACE_GUARD_ADDRESS";
pub const SUB800_RAW_PREFIX_PRESERVED_LENDER_FLAG: &str =
    "LOWQ_SUB800_RAW_PREFIX_PRESERVED_LENDER";
pub const SUB800_RAW_PREFIX_PREDICATE_LENDER_FLAG: &str =
    "LOWQ_SUB800_RAW_PREFIX_PREDICATE_LENDER";
pub const SUB800_MIXED_BOUNDARY_SCRATCH_EXTENSION_FLAG: &str =
    "LOWQ_SUB800_MIXED_BOUNDARY_SCRATCH_EXTENSION";
pub const SUB800_BORROWED_ROTATED_UNDERFLOW_FLAG: &str =
    "LOWQ_SUB800_BORROWED_ROTATED_UNDERFLOW";
pub const SUB800_SPLIT_MIXED_ROTATED_LENGTH_FLAG: &str =
    "LOWQ_SUB800_SPLIT_MIXED_ROTATED_LENGTH";
pub const SUB800_SPLIT_SAME_ROTATED_LENGTH_FLAG: &str =
    "LOWQ_SUB800_SPLIT_SAME_ROTATED_LENGTH";
pub const SUB800_SPLIT_TWO_HIGH_ROTATED_LENGTH_FLAG: &str =
    "LOWQ_SUB800_SPLIT_TWO_HIGH_ROTATED_LENGTH";
pub const SUB800_SPLIT_THREE_HIGH_ROTATED_LENGTH_FLAG: &str =
    "LOWQ_SUB800_SPLIT_THREE_HIGH_ROTATED_LENGTH";
pub const SUB800_SPLIT_FOUR_HIGH_ROTATED_LENGTH_FLAG: &str =
    "LOWQ_SUB800_SPLIT_FOUR_HIGH_ROTATED_LENGTH";
pub const SUB800_ULS_CLEAN_LENDER_FLAG: &str = "LOWQ_SUB800_ULS_CLEAN_LENDER";
pub const SUB800_ULS_FUSED_TARGET_FLAG: &str = "LOWQ_SUB800_ULS_FUSED_TARGET";
pub const SUB800_ULS_DIRECT_SELECTOR_FLAG: &str = "LOWQ_SUB800_ULS_DIRECT_SELECTOR";
pub const Q839_SEVEN_PLATEAU_LENDERS_FLAG: &str =
    "LOWQ_Q839_SEVEN_PLATEAU_LENDERS";

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct Sub800Q839RouteCoverage {
    pub uls_fused_forward_calls: usize,
    pub uls_fused_reverse_calls: usize,
    pub split_three_same_calls: usize,
    pub split_three_mixed_calls: usize,
    pub uls_direct_forward_calls: usize,
    pub uls_direct_reverse_calls: usize,
    pub split_four_same_calls: usize,
    pub split_four_mixed_calls: usize,
    pub seven_plateau_uls_forward_loans: usize,
    pub seven_plateau_uls_reverse_loans: usize,
    pub seven_plateau_short_increments: usize,
    pub seven_plateau_short_decrements: usize,
    pub seven_plateau_support_loans: usize,
    pub seven_plateau_support_fallbacks: usize,
}

static SUB800_Q839_ULS_FUSED_FORWARD_CALLS: AtomicUsize = AtomicUsize::new(0);
static SUB800_Q839_ULS_FUSED_REVERSE_CALLS: AtomicUsize = AtomicUsize::new(0);
static SUB800_Q839_SPLIT_THREE_SAME_CALLS: AtomicUsize = AtomicUsize::new(0);
static SUB800_Q839_SPLIT_THREE_MIXED_CALLS: AtomicUsize = AtomicUsize::new(0);
static SUB800_Q838_ULS_DIRECT_FORWARD_CALLS: AtomicUsize = AtomicUsize::new(0);
static SUB800_Q838_ULS_DIRECT_REVERSE_CALLS: AtomicUsize = AtomicUsize::new(0);
static SUB800_Q838_SPLIT_FOUR_SAME_CALLS: AtomicUsize = AtomicUsize::new(0);
static SUB800_Q838_SPLIT_FOUR_MIXED_CALLS: AtomicUsize = AtomicUsize::new(0);
static Q839_SEVEN_PLATEAU_ULS_FORWARD_LOANS: AtomicUsize = AtomicUsize::new(0);
static Q839_SEVEN_PLATEAU_ULS_REVERSE_LOANS: AtomicUsize = AtomicUsize::new(0);
static Q839_SEVEN_PLATEAU_SHORT_INCREMENTS: AtomicUsize = AtomicUsize::new(0);
static Q839_SEVEN_PLATEAU_SHORT_DECREMENTS: AtomicUsize = AtomicUsize::new(0);
static Q839_SEVEN_PLATEAU_SUPPORT_LOANS: AtomicUsize = AtomicUsize::new(0);
static Q839_SEVEN_PLATEAU_SUPPORT_FALLBACKS: AtomicUsize = AtomicUsize::new(0);

pub fn reset_sub800_q839_route_coverage() {
    SUB800_Q839_ULS_FUSED_FORWARD_CALLS.store(0, Ordering::Relaxed);
    SUB800_Q839_ULS_FUSED_REVERSE_CALLS.store(0, Ordering::Relaxed);
    SUB800_Q839_SPLIT_THREE_SAME_CALLS.store(0, Ordering::Relaxed);
    SUB800_Q839_SPLIT_THREE_MIXED_CALLS.store(0, Ordering::Relaxed);
    SUB800_Q838_ULS_DIRECT_FORWARD_CALLS.store(0, Ordering::Relaxed);
    SUB800_Q838_ULS_DIRECT_REVERSE_CALLS.store(0, Ordering::Relaxed);
    SUB800_Q838_SPLIT_FOUR_SAME_CALLS.store(0, Ordering::Relaxed);
    SUB800_Q838_SPLIT_FOUR_MIXED_CALLS.store(0, Ordering::Relaxed);
    Q839_SEVEN_PLATEAU_ULS_FORWARD_LOANS.store(0, Ordering::Relaxed);
    Q839_SEVEN_PLATEAU_ULS_REVERSE_LOANS.store(0, Ordering::Relaxed);
    Q839_SEVEN_PLATEAU_SHORT_INCREMENTS.store(0, Ordering::Relaxed);
    Q839_SEVEN_PLATEAU_SHORT_DECREMENTS.store(0, Ordering::Relaxed);
    Q839_SEVEN_PLATEAU_SUPPORT_LOANS.store(0, Ordering::Relaxed);
    Q839_SEVEN_PLATEAU_SUPPORT_FALLBACKS.store(0, Ordering::Relaxed);
}

#[must_use]
pub fn sub800_q839_route_coverage() -> Sub800Q839RouteCoverage {
    Sub800Q839RouteCoverage {
        uls_fused_forward_calls: SUB800_Q839_ULS_FUSED_FORWARD_CALLS.load(Ordering::Relaxed),
        uls_fused_reverse_calls: SUB800_Q839_ULS_FUSED_REVERSE_CALLS.load(Ordering::Relaxed),
        split_three_same_calls: SUB800_Q839_SPLIT_THREE_SAME_CALLS.load(Ordering::Relaxed),
        split_three_mixed_calls: SUB800_Q839_SPLIT_THREE_MIXED_CALLS.load(Ordering::Relaxed),
        uls_direct_forward_calls: SUB800_Q838_ULS_DIRECT_FORWARD_CALLS.load(Ordering::Relaxed),
        uls_direct_reverse_calls: SUB800_Q838_ULS_DIRECT_REVERSE_CALLS.load(Ordering::Relaxed),
        split_four_same_calls: SUB800_Q838_SPLIT_FOUR_SAME_CALLS.load(Ordering::Relaxed),
        split_four_mixed_calls: SUB800_Q838_SPLIT_FOUR_MIXED_CALLS.load(Ordering::Relaxed),
        seven_plateau_uls_forward_loans: Q839_SEVEN_PLATEAU_ULS_FORWARD_LOANS
            .load(Ordering::Relaxed),
        seven_plateau_uls_reverse_loans: Q839_SEVEN_PLATEAU_ULS_REVERSE_LOANS
            .load(Ordering::Relaxed),
        seven_plateau_short_increments: Q839_SEVEN_PLATEAU_SHORT_INCREMENTS
            .load(Ordering::Relaxed),
        seven_plateau_short_decrements: Q839_SEVEN_PLATEAU_SHORT_DECREMENTS
            .load(Ordering::Relaxed),
        seven_plateau_support_loans: Q839_SEVEN_PLATEAU_SUPPORT_LOANS
            .load(Ordering::Relaxed),
        seven_plateau_support_fallbacks: Q839_SEVEN_PLATEAU_SUPPORT_FALLBACKS
            .load(Ordering::Relaxed),
    }
}

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
pub struct Sub800InplaceGuardAddressProofReport {
    pub stages_checked: usize,
    pub basis_states_checked: usize,
    pub coordinate_checks: usize,
    pub inverse_pair_checks: usize,
    pub scratch_clean_checks: usize,
    pub phase_clean_checks: usize,
    pub ancilla_clean_checks: usize,
    pub allocated_address_lanes: usize,
    pub prepare: RegisterSharedGateCounts,
    pub reverse_boundary: RegisterSharedGateCounts,
    pub full_roundtrip: RegisterSharedGateCounts,
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
pub struct ReferencePhaseOverlaidRemainderScratchProofReport {
    pub length_widths_checked: usize,
    pub window_configurations_checked: usize,
    pub basis_states_checked: usize,
    pub baseline_equivalence_checks: usize,
    pub inverse_pair_checks: usize,
    pub phase_boundary_clean_checks: usize,
    pub ancilla_clean_checks: usize,
    pub control_combinations_checked: usize,
    pub zero_remainder_checks: usize,
    pub range_equality_checks: usize,
    pub range_boundary_checks: usize,
    pub baseline_reference9_scratch_lanes: usize,
    pub overlaid_reference9_scratch_lanes: usize,
    pub scratch_lanes_saved: usize,
    pub baseline_reference9_peak_qubits: usize,
    pub overlaid_reference9_peak_qubits: usize,
    pub baseline_add257: RegisterSharedGateCounts,
    pub overlaid_add257: RegisterSharedGateCounts,
    pub baseline_sub257: RegisterSharedGateCounts,
    pub overlaid_sub257: RegisterSharedGateCounts,
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
pub struct ReferenceBorrowedCoefficientComparatorProofReport {
    pub widths_checked: usize,
    pub basis_states_checked: usize,
    pub control_off_identity_checks: usize,
    pub equality_boundary_checks: usize,
    pub subtraction_underflow_checks: usize,
    pub addition_overflow_checks: usize,
    pub oracle_equivalence_checks: usize,
    pub inverse_pair_checks: usize,
    pub operand_restore_checks: usize,
    pub ancilla_clean_checks: usize,
    pub baseline_reference9: RegisterSharedGateCounts,
    pub borrowed_reference9: RegisterSharedGateCounts,
    pub baseline_reference9_active_qubits: usize,
    pub baseline_reference9_peak_qubits: usize,
    pub baseline_reference9_temporary_qubits: usize,
    pub borrowed_reference9_active_qubits: usize,
    pub borrowed_reference9_peak_qubits: usize,
    pub borrowed_reference9_temporary_qubits: usize,
    pub borrowed_caller_lanes: usize,
    pub borrowed_reference9_fresh_qubits: usize,
    pub reference9_toffoli_reduction: usize,
    pub reference9_standalone_peak_reduction: usize,
    pub production_incremental_caller_qubits: usize,
    pub reference9_production_local_peak_reduction: usize,
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct Q847LifetimeLocalResources {
    pub active_qubits: usize,
    pub peak_qubits: usize,
    pub temporary_qubits: usize,
    pub emitted_ops: usize,
    pub emitted_toffoli: usize,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct PromisedLqSwapBorrowProofReport {
    pub configurations_checked: usize,
    pub basis_states_checked: usize,
    pub scalar_equivalence_checks: usize,
    pub simulator_equivalence_checks: usize,
    pub control_off_checks: usize,
    pub control_off_nonzero_l_q_checks: usize,
    pub control_off_invalid_support_checks: usize,
    pub control_on_promised_support_checks: usize,
    pub lender_restore_checks: usize,
    pub inverse_pair_checks: usize,
    pub phase_clean_checks: usize,
    pub ancilla_clean_checks: usize,
    pub default_stream_identity_checks: usize,
    pub reference_lender_lanes: usize,
    pub reset_ops_removed_per_invocation: usize,
    pub whole_point_add_invocations: usize,
    pub whole_point_add_ops_delta: i64,
    pub whole_point_add_toffoli_delta: i64,
    pub baseline_local: Q847LifetimeLocalResources,
    pub candidate_local: Q847LifetimeLocalResources,
    pub local_qubit_delta: i64,
    pub local_ops_delta: i64,
    pub local_toffoli_delta: i64,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct PromisedSwapSupportLifetimeFusionProofReport {
    pub configurations_checked: usize,
    pub basis_states_checked: usize,
    pub scalar_equivalence_checks: usize,
    pub simulator_equivalence_checks: usize,
    pub control_off_checks: usize,
    pub control_off_nonzero_l_q_checks: usize,
    pub control_on_promised_support_checks: usize,
    pub excluded_mixed_width_overflow_states: usize,
    pub lender_restore_checks: usize,
    pub inverse_pair_checks: usize,
    pub phase_clean_checks: usize,
    pub ancilla_clean_checks: usize,
    pub default_stream_identity_checks: usize,
    pub baseline_local: Q847LifetimeLocalResources,
    pub candidate_local: Q847LifetimeLocalResources,
    pub local_qubit_delta: i64,
    pub local_ops_delta: i64,
    pub local_toffoli_delta: i64,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct Q845SwapOnlyCoefficientProofReport {
    pub dependency_checks: usize,
    pub full_feature_environment_checks: usize,
    pub truncated_guard_default_off_stream_identity_checks: usize,
    pub truncated_range_comparator_cases_checked: usize,
    pub truncated_guard_layout_cases_checked: usize,
    pub truncated_guard_trace_checks: usize,
    pub production_truncated_address_cases_checked: usize,
    pub layout_cases_checked: usize,
    pub internal_boundary_layout_cases_checked: usize,
    pub promised_basis_states_checked: usize,
    pub oracle_transition_checks: usize,
    pub inverse_pair_checks: usize,
    pub scratch_clean_checks: usize,
    pub cursor_restore_checks: usize,
    pub count_restore_checks: usize,
    pub residue_preservation_checks: usize,
    pub excluded_suffix_preservation_checks: usize,
    pub default_off_stream_identity_checks: usize,
    pub ephemeral_swap_cases_checked: usize,
    pub ephemeral_control_on_checks: usize,
    pub ephemeral_control_off_checks: usize,
    pub persistent_lifecycle_equivalence_checks: usize,
    pub ephemeral_inverse_pair_checks: usize,
    pub ephemeral_l_t_prime_zero_checks: usize,
    pub ephemeral_phase_clean_checks: usize,
    pub ephemeral_ancilla_clean_checks: usize,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct Q851FixedSignEventProofReport {
    pub identity_pairs_checked: usize,
    pub body_256_checks: usize,
    pub production_schedule_widths_checked: usize,
    pub domain_fallback_stream_identity_checks: usize,
    pub route_branch_cases_checked: usize,
    pub route_transition_index_checks: usize,
    pub route_body_sign_observation_checks: usize,
    pub route_cursor_restore_checks: usize,
    pub transition_events_checked: usize,
    pub transition_basis_states_checked: usize,
    pub direction_stream_identity_checks: usize,
    pub sequence_widths_checked: usize,
    pub forward_sequence_cases_checked: usize,
    pub reverse_sequence_cases_checked: usize,
    pub exact_reverse_stream_checks: usize,
    pub cursor_restore_checks: usize,
    pub scratch_clean_checks: usize,
    pub phase_clean_checks: usize,
    pub ancilla_clean_checks: usize,
    pub default_off_stream_identity_checks: usize,
    pub allocation_free_microkernels_checked: usize,
    pub allocation_free_sequence_transitions_checked: usize,
    pub transition_toffoli: usize,
    pub transition_x_min: usize,
    pub transition_x_max: usize,
    pub transition_ops_min: usize,
    pub transition_ops_max: usize,
    pub baseline_transition9: RegisterSharedGateCounts,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct SplitCoefficientRotationLifetimeProofReport {
    pub configurations_checked: usize,
    pub lender_modes_checked: usize,
    pub basis_states_checked: usize,
    pub scalar_equivalence_checks: usize,
    pub simulator_equivalence_checks: usize,
    pub control_off_checks: usize,
    pub lender_restore_checks: usize,
    pub inverse_pair_checks: usize,
    pub phase_clean_checks: usize,
    pub ancilla_clean_checks: usize,
    pub default_stream_identity_checks: usize,
    pub split_release_checks: usize,
    pub split_recompute_checks: usize,
    pub reference_rotation_lanes_released: usize,
    pub whole_point_add_invocations: usize,
    pub whole_point_add_ops_delta: i64,
    pub whole_point_add_toffoli_delta: i64,
    pub baseline_local: Q847LifetimeLocalResources,
    pub candidate_local: Q847LifetimeLocalResources,
    pub local_qubit_delta: i64,
    pub local_ops_delta: i64,
    pub local_toffoli_delta: i64,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct CoefficientLessThanLaneReuseProofReport {
    pub configurations_checked: usize,
    pub basis_states_checked: usize,
    pub scalar_equivalence_checks: usize,
    pub simulator_equivalence_checks: usize,
    pub control_off_checks: usize,
    pub boundary_restore_checks: usize,
    pub lender_restore_checks: usize,
    pub inverse_pair_checks: usize,
    pub phase_clean_checks: usize,
    pub ancilla_clean_checks: usize,
    pub default_stream_identity_checks: usize,
    pub alias_rejections: usize,
    pub caller_lanes_reused: usize,
    pub reset_ops_removed_per_invocation: usize,
    pub whole_point_add_invocations: usize,
    pub whole_point_add_ops_delta: i64,
    pub baseline_local: Q847LifetimeLocalResources,
    pub candidate_local: Q847LifetimeLocalResources,
    pub local_qubit_delta: i64,
    pub local_ops_delta: i64,
    pub local_toffoli_delta: i64,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct CleanChainCoefficientAddLenderProofReport {
    pub configurations_checked: usize,
    pub directions_checked: usize,
    pub basis_states_checked: usize,
    pub scalar_equivalence_checks: usize,
    pub simulator_equivalence_checks: usize,
    pub control_off_identity_checks: usize,
    pub length_preservation_checks: usize,
    pub lender_clean_entry_checks: usize,
    pub lender_restore_checks: usize,
    pub roundtrip_checks: usize,
    pub lender_window_phase_clean_checks: usize,
    pub ancilla_clean_checks: usize,
    pub default_stream_identity_checks: usize,
    pub distinctness_rejections: usize,
    pub composition_basis_states_checked: usize,
    pub composition_equivalence_checks: usize,
    pub composition_lender_clean_entry_checks: usize,
    pub composition_lender_restore_checks: usize,
    pub composition_lender_window_phase_clean_checks: usize,
    pub composition_ancilla_clean_checks: usize,
    pub caller_lanes_reused: usize,
    pub reset_ops_removed_per_invocation: usize,
    pub whole_point_add_invocations: usize,
    pub whole_point_add_ops_delta: i64,
    pub baseline_local: Q847LifetimeLocalResources,
    pub candidate_local: Q847LifetimeLocalResources,
    pub local_qubit_delta: i64,
    pub local_ops_delta: i64,
    pub local_toffoli_delta: i64,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct InPlaceLtCursorLegacyDifferentialProofReport {
    pub less_than_configurations_checked: usize,
    pub coefficient_add_configurations_checked: usize,
    pub coefficient_add_directions_checked: usize,
    pub basis_states_checked: usize,
    pub scalar_equivalence_checks: usize,
    pub simulator_equivalence_checks: usize,
    pub control_off_checks: usize,
    pub length_restore_checks: usize,
    pub inverse_pair_checks: usize,
    pub phase_clean_checks: usize,
    pub ancilla_clean_checks: usize,
    pub legacy_streams_checked: usize,
    pub copied_cursor_lanes_removed: usize,
    pub local_ops_removed: usize,
    pub local_cx_removed: usize,
    pub local_resets_removed: usize,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct PreservedDyTopPrefixLoanProofReport {
    pub configurations_checked: usize,
    pub directions_checked: usize,
    pub basis_states_checked: usize,
    pub scalar_equivalence_checks: usize,
    pub simulator_equivalence_checks: usize,
    pub control_off_checks: usize,
    pub lender_clean_entry_checks: usize,
    pub lender_restore_checks: usize,
    pub inverse_pair_checks: usize,
    pub phase_clean_checks: usize,
    pub ancilla_clean_checks: usize,
    pub borrow_windows_checked: usize,
    pub alias_rejections: usize,
    pub baseline_owned_prefix_lanes: usize,
    pub candidate_owned_prefix_lanes: usize,
    pub local_qubit_delta: i64,
    pub local_ops_delta: i64,
    pub local_toffoli_delta: i64,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct MixedLrPrimeProofReport {
    pub boundary_configurations_checked: usize,
    pub swap_configurations_checked: usize,
    pub directions_checked: usize,
    pub boundary_basis_states_checked: usize,
    pub swap_basis_states_checked: usize,
    pub scalar_equivalence_checks: usize,
    pub simulator_equivalence_checks: usize,
    pub control_off_checks: usize,
    pub unsupported_control_on_states: usize,
    pub omitted_high_lane_clean_checks: usize,
    pub inverse_pair_checks: usize,
    pub phase_clean_checks: usize,
    pub ancilla_clean_checks: usize,
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

/// `b += a mod 2^n` without materializing the discarded carry-out.
fn cuccaro_add_mod_2n_no_overflow(circ: &mut Circuit, a: &[QReg], b: &[QReg], carry: &QReg) {
    assert_eq!(a.len(), b.len());
    assert!(!a.is_empty());
    majority(circ, &a[0], &b[0], carry);
    for index in 1..a.len() {
        majority(circ, &a[index], &b[index], &a[index - 1]);
    }
    for index in (1..a.len()).rev() {
        unmajority_add(circ, &a[index], &b[index], &a[index - 1]);
    }
    unmajority_add(circ, &a[0], &b[0], carry);
}

/// Exact inverse of [`cuccaro_add_mod_2n_no_overflow`].
fn cuccaro_sub_mod_2n_no_overflow(circ: &mut Circuit, a: &[QReg], b: &[QReg], carry: &QReg) {
    assert_eq!(a.len(), b.len());
    assert!(!a.is_empty());
    unmajority_add_inverse(circ, &a[0], &b[0], carry);
    for index in 1..a.len() {
        unmajority_add_inverse(circ, &a[index], &b[index], &a[index - 1]);
    }
    for index in (1..a.len()).rev() {
        majority_inverse(circ, &a[index], &b[index], &a[index - 1]);
    }
    majority_inverse(circ, &a[0], &b[0], carry);
}

fn cuccaro_add_mod_2n_no_overflow_refs(
    circ: &mut Circuit,
    a: &[&QReg],
    b: &[QReg],
    carry: &QReg,
) {
    assert_eq!(a.len(), b.len());
    assert!(!a.is_empty());
    majority(circ, a[0], &b[0], carry);
    for index in 1..a.len() {
        majority(circ, a[index], &b[index], a[index - 1]);
    }
    for index in (1..a.len()).rev() {
        unmajority_add(circ, a[index], &b[index], a[index - 1]);
    }
    unmajority_add(circ, a[0], &b[0], carry);
}

fn cuccaro_sub_mod_2n_no_overflow_refs(
    circ: &mut Circuit,
    a: &[&QReg],
    b: &[QReg],
    carry: &QReg,
) {
    assert_eq!(a.len(), b.len());
    assert!(!a.is_empty());
    unmajority_add_inverse(circ, a[0], &b[0], carry);
    for index in 1..a.len() {
        unmajority_add_inverse(circ, a[index], &b[index], a[index - 1]);
    }
    for index in (1..a.len()).rev() {
        majority_inverse(circ, a[index], &b[index], a[index - 1]);
    }
    majority_inverse(circ, a[0], &b[0], carry);
}

/// Add an `(n-1)`-bit value into an `n`-bit target using a caller-owned clean
/// zero extension. The extension and carry are restored.
fn cuccaro_add_zero_extended_no_overflow(
    circ: &mut Circuit,
    short: &[QReg],
    target: &[QReg],
    carry: &QReg,
    zero_extension: &QReg,
) {
    assert_eq!(short.len() + 1, target.len());
    assert!(short.iter().all(|lane| lane.id() != zero_extension.id()));
    assert!(target.iter().all(|lane| lane.id() != zero_extension.id()));
    let source: Vec<&QReg> = short.iter().chain(std::iter::once(zero_extension)).collect();
    cuccaro_add_mod_2n_no_overflow_refs(circ, &source, target, carry);
}

fn cuccaro_sub_zero_extended_no_overflow(
    circ: &mut Circuit,
    short: &[QReg],
    target: &[QReg],
    carry: &QReg,
    zero_extension: &QReg,
) {
    assert_eq!(short.len() + 1, target.len());
    assert!(short.iter().all(|lane| lane.id() != zero_extension.id()));
    assert!(target.iter().all(|lane| lane.id() != zero_extension.id()));
    let source: Vec<&QReg> = short.iter().chain(std::iter::once(zero_extension)).collect();
    cuccaro_sub_mod_2n_no_overflow_refs(circ, &source, target, carry);
}

/// Controlled increment with a caller-composed reference slice of clean carries.
fn controlled_increment_mod_2n_carry_refs(
    circ: &mut Circuit,
    control: &QReg,
    register: &[QReg],
    carries: &[&QReg],
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
    circ.ccx(control, &register[0], carries[0]);
    circ.cx(control, &register[0]);
    for index in 1..width - 1 {
        circ.ccx(&register[index], carries[index - 1], carries[index]);
    }
    circ.cx(carries[width - 2], &register[width - 1]);
    for index in (1..width - 1).rev() {
        circ.ccx(&register[index], carries[index - 1], carries[index]);
        circ.cx(carries[index - 1], &register[index]);
    }
    circ.cx(control, &register[0]);
    circ.ccx(control, &register[0], carries[0]);
    circ.cx(control, &register[0]);
}

/// Controlled decrement with a caller-composed reference slice of clean borrows.
fn controlled_decrement_mod_2n_carry_refs(
    circ: &mut Circuit,
    control: &QReg,
    register: &[QReg],
    borrows: &[&QReg],
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
    circ.ccx(control, &register[0], borrows[0]);
    circ.x(&register[0]);
    circ.cx(control, &register[0]);
    for index in 1..width - 1 {
        circ.x(&register[index]);
        circ.ccx(&register[index], borrows[index - 1], borrows[index]);
        circ.x(&register[index]);
    }
    circ.cx(borrows[width - 2], &register[width - 1]);
    for index in (1..width - 1).rev() {
        circ.x(&register[index]);
        circ.ccx(&register[index], borrows[index - 1], borrows[index]);
        circ.x(&register[index]);
        circ.cx(borrows[index - 1], &register[index]);
    }
    circ.cx(control, &register[0]);
    circ.x(&register[0]);
    circ.ccx(control, &register[0], borrows[0]);
    circ.x(&register[0]);
    circ.cx(control, &register[0]);
}

/// Borrowed-view Cuccaro add used to append a clean zero-extension lane.
fn cuccaro_add_mod_2n_refs(
    circ: &mut Circuit,
    a: &[&QReg],
    b: &[QReg],
    carry: &QReg,
    overflow: &QReg,
) {
    assert_eq!(a.len(), b.len());
    assert!(!a.is_empty());
    majority(circ, a[0], &b[0], carry);
    for index in 1..a.len() {
        majority(circ, a[index], &b[index], a[index - 1]);
    }
    circ.cx(a[a.len() - 1], overflow);
    for index in (1..a.len()).rev() {
        unmajority_add(circ, a[index], &b[index], a[index - 1]);
    }
    unmajority_add(circ, a[0], &b[0], carry);
}

/// Exact inverse of [`cuccaro_add_mod_2n_refs`].
fn cuccaro_sub_mod_2n_refs(
    circ: &mut Circuit,
    a: &[&QReg],
    b: &[QReg],
    carry: &QReg,
    overflow: &QReg,
) {
    assert_eq!(a.len(), b.len());
    assert!(!a.is_empty());
    unmajority_add_inverse(circ, a[0], &b[0], carry);
    for index in 1..a.len() {
        unmajority_add_inverse(circ, a[index], &b[index], a[index - 1]);
    }
    circ.cx(a[a.len() - 1], overflow);
    for index in (1..a.len()).rev() {
        majority_inverse(circ, a[index], &b[index], a[index - 1]);
    }
    majority_inverse(circ, a[0], &b[0], carry);
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
    assert_eq!(l_t.len(), l_q.len() + 1);
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
    cuccaro_add_zero_extended_no_overflow(circ, l_q, l_t, carry, overflow);
    controlled_dynamic_swaps(
        circ,
        phase2,
        sign,
        work1_window,
        l_t,
        first_global_index,
        flag,
        chain,
        false,
    );
    cuccaro_sub_zero_extended_no_overflow(circ, l_q, l_t, carry, overflow);
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
    assert_eq!(l_t.len(), l_q.len() + 1);
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
    cuccaro_add_zero_extended_no_overflow(circ, l_q, l_t, carry, overflow);
    controlled_dynamic_swaps(
        circ,
        phase2,
        sign,
        work1_window,
        l_t,
        first_global_index,
        flag,
        chain,
        true,
    );
    cuccaro_sub_zero_extended_no_overflow(circ, l_q, l_t, carry, overflow);
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

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum CoefficientNonnegativeBracket {
    Legacy,
    CancelMiddleX,
}

fn production_coefficient_nonnegative_bracket() -> CoefficientNonnegativeBracket {
    if coefficient_nonnegative_x_cancel_requested() {
        CoefficientNonnegativeBracket::CancelMiddleX
    } else {
        CoefficientNonnegativeBracket::Legacy
    }
}

fn begin_coefficient_nonnegative_bracket(
    circ: &mut Circuit,
    control: &QReg,
    signed_length: &[QReg],
    active: &QReg,
    bracket: CoefficientNonnegativeBracket,
) {
    let sign = signed_length.last().expect("nonempty signed length");
    begin_coefficient_nonnegative_sign_bracket(circ, control, sign, active, bracket);
}

fn begin_coefficient_nonnegative_sign_bracket(
    circ: &mut Circuit,
    control: &QReg,
    sign: &QReg,
    active: &QReg,
    bracket: CoefficientNonnegativeBracket,
) {
    if bracket == CoefficientNonnegativeBracket::Legacy {
        circ.x(sign);
        circ.ccx(control, sign, active);
        circ.x(sign);
        return;
    }
    circ.x(sign);
    circ.ccx(control, sign, active);
}

fn end_coefficient_nonnegative_bracket(
    circ: &mut Circuit,
    control: &QReg,
    signed_length: &[QReg],
    active: &QReg,
    bracket: CoefficientNonnegativeBracket,
) {
    let sign = signed_length.last().expect("nonempty signed length");
    end_coefficient_nonnegative_sign_bracket(circ, control, sign, active, bracket);
}

fn end_coefficient_nonnegative_sign_bracket(
    circ: &mut Circuit,
    control: &QReg,
    sign: &QReg,
    active: &QReg,
    bracket: CoefficientNonnegativeBracket,
) {
    if bracket == CoefficientNonnegativeBracket::Legacy {
        circ.x(sign);
        circ.ccx(control, sign, active);
        circ.x(sign);
        return;
    }
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
    constant_overflow: &'a QReg,
    walk: &'a [QReg],
    nonzero: &'a QReg,
    operation: &'a QReg,
    phase_sign: &'a QReg,
    enable: &'a QReg,
    nonzero_chain: &'a [QReg],
}

pub const PHASE_OVERLAID_REMAINDER_SCRATCH_FLAG: &str = "LOWQ_Q839_PHASE_REMAINDER_SCRATCH";

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum RemainderScratchLayout {
    Baseline,
    PhaseOverlaid,
}

fn phase_overlaid_remainder_scratch_requested() -> bool {
    std::env::var(PHASE_OVERLAID_REMAINDER_SCRATCH_FLAG)
        .ok()
        .as_deref()
        == Some("1")
}

fn selected_remainder_scratch_layout() -> RemainderScratchLayout {
    if phase_overlaid_remainder_scratch_requested() {
        RemainderScratchLayout::PhaseOverlaid
    } else {
        RemainderScratchLayout::Baseline
    }
}

fn baseline_remainder_scratch_width(length_width: usize, remainder_length_width: usize) -> usize {
    5 + (length_width + 2)
        + length_width.saturating_sub(1)
        + 4
        + remainder_length_width.saturating_sub(2)
}

fn phase_overlaid_remainder_scratch_width(
    length_width: usize,
    remainder_length_width: usize,
) -> usize {
    (length_width + 6).max(remainder_length_width).max(8)
}

fn remainder_scratch_width_for_layout(
    length_width: usize,
    remainder_length_width: usize,
    layout: RemainderScratchLayout,
) -> usize {
    match layout {
        RemainderScratchLayout::Baseline => {
            baseline_remainder_scratch_width(length_width, remainder_length_width)
        }
        RemainderScratchLayout::PhaseOverlaid => {
            phase_overlaid_remainder_scratch_width(length_width, remainder_length_width)
        }
    }
}

fn remainder_scratch_width(length_width: usize, remainder_length_width: usize) -> usize {
    remainder_scratch_width_for_layout(
        length_width,
        remainder_length_width,
        selected_remainder_scratch_layout(),
    )
}

fn assert_unique_qreg_ids(label: &str, lanes: &[&QReg]) {
    let mut ids = std::collections::HashSet::with_capacity(lanes.len());
    for (index, lane) in lanes.iter().enumerate() {
        assert!(
            ids.insert(lane.id()),
            "{label}: lane {index} aliases another simultaneously live lane"
        );
    }
}

fn assert_remainder_scratch_disjoint_from_data(
    label: &str,
    controls: &[&QReg],
    registers: &[&[QReg]],
    scratch: &[QReg],
) {
    let mut lanes = Vec::with_capacity(
        controls.len()
            + registers
                .iter()
                .map(|register| register.len())
                .sum::<usize>()
            + scratch.len(),
    );
    lanes.extend_from_slice(controls);
    for register in registers {
        lanes.extend(register.iter());
    }
    lanes.extend(scratch.iter());
    assert_unique_qreg_ids(label, &lanes);
}

fn split_baseline_remainder_scratch<'a>(
    scratch: &'a [QReg],
    length_width: usize,
    remainder_length_width: usize,
) -> RemainderScratch<'a> {
    assert!(
        scratch.len() >= baseline_remainder_scratch_width(length_width, remainder_length_width)
    );
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
    let constant_overflow = constant.last().expect("constant scratch is nonempty");
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
        constant_overflow,
        walk,
        nonzero,
        operation,
        phase_sign,
        enable,
        nonzero_chain,
    }
}

fn assert_phase_overlaid_remainder_scratch_layout(
    scratch: &[QReg],
    view: &RemainderScratch<'_>,
    length_width: usize,
    remainder_length_width: usize,
) {
    assert_eq!(
        scratch.len(),
        phase_overlaid_remainder_scratch_width(length_width, remainder_length_width),
        "phase-overlaid remainder scratch requires its exact lane count"
    );

    assert_eq!(view.phase_sign.id(), view.length_overflow.id());
    assert_eq!(view.carry.id(), view.length_carry.id());
    assert_eq!(view.carry.id(), view.constant[0].id());
    assert_eq!(
        view.constant_overflow.id(),
        view.constant[length_width + 1].id()
    );
    assert_eq!(
        view.walk.iter().map(QReg::id).collect::<Vec<_>>(),
        view.constant[1..length_width]
            .iter()
            .map(QReg::id)
            .collect::<Vec<_>>()
    );
    match length_width {
        1 => {
            assert_eq!(view.active.id(), view.constant[1].id());
            assert_eq!(view.tmp.id(), scratch[7].id());
        }
        2 => {
            assert_eq!(view.active.id(), view.walk[0].id());
            assert_eq!(view.tmp.id(), view.constant[2].id());
        }
        _ => {
            assert_eq!(view.active.id(), view.walk[0].id());
            assert_eq!(view.tmp.id(), view.walk[1].id());
        }
    }

    let mut nonzero_phase = vec![view.nonzero, view.operation];
    nonzero_phase.extend(view.nonzero_chain.iter());
    assert_unique_qreg_ids("remainder nonzero phase", &nonzero_phase);

    assert_unique_qreg_ids(
        "remainder add-enable phase",
        &[view.nonzero, view.operation, view.enable, view.phase_sign],
    );

    assert_unique_qreg_ids(
        "remainder prepare-add phase",
        &[
            view.nonzero,
            view.operation,
            view.enable,
            view.length_carry,
            view.length_overflow,
        ],
    );

    let mut constant_phase = vec![
        view.nonzero,
        view.operation,
        view.enable,
        view.length_overflow,
    ];
    constant_phase.extend(view.constant.iter());
    assert_unique_qreg_ids("remainder constant phase", &constant_phase);

    assert_unique_qreg_ids(
        "remainder arithmetic phase",
        &[
            view.nonzero,
            view.operation,
            view.enable,
            view.length_overflow,
            view.constant_overflow,
            view.carry,
            view.active,
            view.tmp,
        ],
    );

    let mut walk_phase = vec![
        view.nonzero,
        view.operation,
        view.enable,
        view.length_overflow,
        view.constant_overflow,
        view.carry,
    ];
    walk_phase.extend(view.walk.iter());
    assert_unique_qreg_ids("remainder walk phase", &walk_phase);
}

fn split_phase_overlaid_remainder_scratch<'a>(
    scratch: &'a [QReg],
    length_width: usize,
    remainder_length_width: usize,
) -> RemainderScratch<'a> {
    assert!(!scratch.is_empty());
    assert!(length_width > 0 && remainder_length_width > 0);
    assert_eq!(
        scratch.len(),
        phase_overlaid_remainder_scratch_width(length_width, remainder_length_width)
    );

    // Lanes 0..=3 hold predicates that cross phases. The constant workspace is
    // reinterpreted only after each clean handoff; width one needs one extra
    // traversal lane because it has no walk-chain lanes to host active/tmp.
    let constant = &scratch[4..4 + length_width + 2];
    let walk = &constant[1..length_width];
    let (active, tmp) = match length_width {
        1 => (&constant[1], &scratch[7]),
        2 => (&walk[0], &constant[2]),
        _ => (&walk[0], &walk[1]),
    };
    let view = RemainderScratch {
        nonzero: &scratch[0],
        operation: &scratch[1],
        enable: &scratch[2],
        phase_sign: &scratch[3],
        length_overflow: &scratch[3],
        carry: &constant[0],
        length_carry: &constant[0],
        walk,
        active,
        tmp,
        constant,
        constant_overflow: &constant[length_width + 1],
        nonzero_chain: &scratch[2..2 + remainder_length_width.saturating_sub(2)],
    };
    assert_phase_overlaid_remainder_scratch_layout(
        scratch,
        &view,
        length_width,
        remainder_length_width,
    );
    view
}

fn split_remainder_scratch_for_layout<'a>(
    scratch: &'a [QReg],
    length_width: usize,
    remainder_length_width: usize,
    layout: RemainderScratchLayout,
) -> RemainderScratch<'a> {
    match layout {
        RemainderScratchLayout::Baseline => {
            split_baseline_remainder_scratch(scratch, length_width, remainder_length_width)
        }
        RemainderScratchLayout::PhaseOverlaid => {
            split_phase_overlaid_remainder_scratch(scratch, length_width, remainder_length_width)
        }
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
    assert_eq!(l_q.len() + 1, l_s.len());
    assert_l_r_prime_metadata_width(l_s.len(), l_r_prime.len());
    assert!(scratch.len() >= normalized_phase_scratch_width(l_s.len()));
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
    assert_eq!(l_q.len() + 1, l_s.len());
    assert_l_r_prime_metadata_width(l_s.len(), l_r_prime.len());
    assert!(scratch.len() >= normalized_phase_scratch_width(l_s.len()));
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

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
struct DynamicBitLengthZeroAllocationTrace {
    flag_allocations: usize,
    carry_allocations: usize,
    prefix_allocations: usize,
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct FusedPrefixScratchLoanAllocationTrace {
    pub calls: usize,
    pub owned_lanes: usize,
    pub borrowed_lanes: usize,
    pub maximum_owned_lanes: usize,
    pub maximum_borrowed_lanes: usize,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
struct PreservedDyTopBorrowWindow {
    entry_ops_idx: usize,
    restore_ops_idx: usize,
    lender_id: u32,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
struct Q839SupportLenderBorrowWindow {
    entry_ops_idx: usize,
    restore_ops_idx: usize,
    lender_id: u32,
}

thread_local! {
    static DYNAMIC_BIT_LENGTH_ZERO_ALLOCATION_TRACE: std::cell::Cell<(
        bool,
        DynamicBitLengthZeroAllocationTrace,
    )> = std::cell::Cell::new((false, DynamicBitLengthZeroAllocationTrace {
        flag_allocations: 0,
        carry_allocations: 0,
        prefix_allocations: 0,
    }));
    static FUSED_PREFIX_SCRATCH_LOAN_ALLOCATION_TRACE: std::cell::Cell<(
        bool,
        FusedPrefixScratchLoanAllocationTrace,
    )> = std::cell::Cell::new((false, FusedPrefixScratchLoanAllocationTrace {
        calls: 0,
        owned_lanes: 0,
        borrowed_lanes: 0,
        maximum_owned_lanes: 0,
        maximum_borrowed_lanes: 0,
    }));
    static PRESERVED_DY_TOP_BORROW_WINDOWS: std::cell::RefCell<(
        bool,
        Vec<PreservedDyTopBorrowWindow>,
    )> = std::cell::RefCell::new((false, Vec::new()));
    static Q839_SUPPORT_LENDER_BORROW_WINDOWS: std::cell::RefCell<(
        bool,
        Vec<Q839SupportLenderBorrowWindow>,
    )> = std::cell::RefCell::new((false, Vec::new()));
}

fn begin_dynamic_bit_length_zero_allocation_trace() {
    DYNAMIC_BIT_LENGTH_ZERO_ALLOCATION_TRACE.with(|trace| {
        trace.set((true, DynamicBitLengthZeroAllocationTrace::default()));
    });
}

fn finish_dynamic_bit_length_zero_allocation_trace() -> DynamicBitLengthZeroAllocationTrace {
    DYNAMIC_BIT_LENGTH_ZERO_ALLOCATION_TRACE.with(|trace| {
        let (_, snapshot) = trace.get();
        trace.set((false, snapshot));
        snapshot
    })
}

fn trace_dynamic_bit_length_zero_allocations(
    flag_allocations: usize,
    carry_allocations: usize,
    prefix_allocations: usize,
) {
    DYNAMIC_BIT_LENGTH_ZERO_ALLOCATION_TRACE.with(|trace| {
        let (enabled, mut snapshot) = trace.get();
        if enabled {
            snapshot.flag_allocations += flag_allocations;
            snapshot.carry_allocations += carry_allocations;
            snapshot.prefix_allocations += prefix_allocations;
            trace.set((enabled, snapshot));
        }
    });
}

fn begin_fused_prefix_scratch_loan_allocation_trace() {
    FUSED_PREFIX_SCRATCH_LOAN_ALLOCATION_TRACE.with(|trace| {
        trace.set((true, FusedPrefixScratchLoanAllocationTrace::default()));
    });
}

fn finish_fused_prefix_scratch_loan_allocation_trace() -> FusedPrefixScratchLoanAllocationTrace {
    FUSED_PREFIX_SCRATCH_LOAN_ALLOCATION_TRACE.with(|trace| {
        let (_, snapshot) = trace.get();
        trace.set((false, snapshot));
        snapshot
    })
}

fn trace_fused_prefix_scratch_loan_allocation(owned_lanes: usize, borrowed_lanes: usize) {
    FUSED_PREFIX_SCRATCH_LOAN_ALLOCATION_TRACE.with(|trace| {
        let (enabled, mut snapshot) = trace.get();
        if enabled {
            snapshot.calls += 1;
            snapshot.owned_lanes += owned_lanes;
            snapshot.borrowed_lanes += borrowed_lanes;
            snapshot.maximum_owned_lanes = snapshot.maximum_owned_lanes.max(owned_lanes);
            snapshot.maximum_borrowed_lanes = snapshot.maximum_borrowed_lanes.max(borrowed_lanes);
            trace.set((enabled, snapshot));
        }
    });
}

fn begin_preserved_dy_top_borrow_trace() {
    PRESERVED_DY_TOP_BORROW_WINDOWS.with(|trace| {
        *trace.borrow_mut() = (true, Vec::new());
    });
}

fn finish_preserved_dy_top_borrow_trace() -> Vec<PreservedDyTopBorrowWindow> {
    PRESERVED_DY_TOP_BORROW_WINDOWS.with(|trace| {
        let mut trace = trace.borrow_mut();
        trace.0 = false;
        std::mem::take(&mut trace.1)
    })
}

fn record_preserved_dy_top_borrow_window(
    lender: &QReg,
    entry_ops_idx: usize,
    restore_ops_idx: usize,
) {
    PRESERVED_DY_TOP_BORROW_WINDOWS.with(|trace| {
        let mut trace = trace.borrow_mut();
        if trace.0 {
            assert!(restore_ops_idx > entry_ops_idx);
            trace.1.push(PreservedDyTopBorrowWindow {
                entry_ops_idx,
                restore_ops_idx,
                lender_id: lender.id(),
            });
        }
    });
}

fn begin_q839_support_lender_borrow_trace() {
    Q839_SUPPORT_LENDER_BORROW_WINDOWS.with(|trace| {
        *trace.borrow_mut() = (true, Vec::new());
    });
}

fn finish_q839_support_lender_borrow_trace() -> Vec<Q839SupportLenderBorrowWindow> {
    Q839_SUPPORT_LENDER_BORROW_WINDOWS.with(|trace| {
        let mut trace = trace.borrow_mut();
        trace.0 = false;
        std::mem::take(&mut trace.1)
    })
}

fn record_q839_support_lender_borrow_window(
    lender: &QReg,
    entry_ops_idx: usize,
    restore_ops_idx: usize,
) {
    Q839_SUPPORT_LENDER_BORROW_WINDOWS.with(|trace| {
        let mut trace = trace.borrow_mut();
        if trace.0 {
            assert!(restore_ops_idx > entry_ops_idx);
            trace.1.push(Q839SupportLenderBorrowWindow {
                entry_ops_idx,
                restore_ops_idx,
                lender_id: lender.id(),
            });
        }
    });
}

fn toggle_static_zero_flag(
    circ: &mut Circuit,
    source: &[&QReg],
    flag: &QReg,
    prefix_scratch: Option<&[&QReg]>,
) {
    use super::shrunken_pz_state_machine::DIRECT_PREFIX_KG_SCRATCH_LEN;
    use crate::point_add::trailmix_port::arith::khattar_gidney::{
        kg_prefix_ancilla_count, KgPrefixAnd,
    };

    if let Some(scratch) = prefix_scratch {
        assert_eq!(
            scratch.len(),
            DIRECT_PREFIX_KG_SCRATCH_LEN,
            "zero-flag prefix scratch requires exactly {DIRECT_PREFIX_KG_SCRATCH_LEN} lanes"
        );
        assert!(
            kg_prefix_ancilla_count(source.len()) <= scratch.len(),
            "zero-flag source width {} exceeds the five-lane KG scratch budget",
            source.len()
        );
        for (index, lane) in scratch.iter().enumerate() {
            assert!(
                scratch[..index].iter().all(|other| other.id() != lane.id()),
                "zero-flag prefix scratch lane {index} aliases an earlier scratch lane"
            );
            assert!(
                source.iter().all(|source| source.id() != lane.id()),
                "zero-flag prefix scratch lane {index} aliases the source"
            );
            assert_ne!(
                lane.id(),
                flag.id(),
                "zero-flag prefix scratch lane {index} aliases the output flag"
            );
        }
    }

    if source.is_empty() {
        circ.x(flag);
        return;
    }
    let prefix_allocation_serial = prefix_scratch.map(|_| circ.b.allocation_serial);
    for bit in source {
        circ.x(bit);
    }
    let ancillae = if prefix_scratch.is_some() {
        Vec::new()
    } else {
        trace_dynamic_bit_length_zero_allocations(0, 0, kg_prefix_ancilla_count(source.len()));
        circ.alloc_qreg_bits(
            "rs.dynamic-bitlen.zero-prefix",
            kg_prefix_ancilla_count(source.len()),
        )
    };
    let ancilla_refs: Vec<&QReg> =
        prefix_scratch.map_or_else(|| ancillae.iter().collect(), |scratch| scratch.to_vec());
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
    if let Some(allocation_serial) = prefix_allocation_serial {
        assert_eq!(
            circ.b.allocation_serial, allocation_serial,
            "caller-supplied zero-flag prefix traversal allocated an internal qubit"
        );
    }
}

fn full_zero_carry_prefix_scratch_requested() -> bool {
    std::env::var("LOWQ_REUSE_ZERO_CARRIES_FOR_FULL_PREFIX_SCRATCH")
        .ok()
        .as_deref()
        == Some("1")
}

fn assert_full_zero_carry_scratch(
    source: &[&QReg],
    output: &[QReg],
    zero: &QReg,
    carries: &[&QReg],
) {
    use super::shrunken_pz_state_machine::DIRECT_PREFIX_FULL_SCRATCH_LEN;

    assert_eq!(
        carries.len(),
        DIRECT_PREFIX_FULL_SCRATCH_LEN,
        "full zero-carry prefix reuse requires exactly {DIRECT_PREFIX_FULL_SCRATCH_LEN} lanes"
    );
    for (index, lane) in carries.iter().enumerate() {
        assert!(
            carries[..index].iter().all(|other| other.id() != lane.id()),
            "zero-carry scratch lane {index} aliases an earlier scratch lane"
        );
        assert!(
            source.iter().all(|source| source.id() != lane.id()),
            "zero-carry scratch lane {index} aliases the source"
        );
        assert!(
            output.iter().all(|output| output.id() != lane.id()),
            "zero-carry scratch lane {index} aliases the output"
        );
        assert_ne!(
            lane.id(),
            zero.id(),
            "zero-carry scratch lane {index} aliases the zero flag"
        );
    }
}

fn bit_length_lean_allow_zero_legacy_with_borrowed_carry(
    circ: &mut Circuit,
    source: &[&QReg],
    output: &[QReg],
    decrement: bool,
    borrowed_carry: Option<&QReg>,
) {
    use super::shrunken_pz_state_machine::{
        bit_length_lean, bit_length_lean_with_full_prefix_scratch,
        bit_length_lean_with_increment_scratch, dirty_controlled_inc_suffix,
        DIRECT_PREFIX_FULL_SCRATCH_LEN, DIRECT_PREFIX_KG_SCRATCH_LEN,
    };

    if source.is_empty() {
        return;
    }
    if let Some(borrowed_carry) = borrowed_carry {
        assert!(
            source
                .iter()
                .all(|source_bit| source_bit.id() != borrowed_carry.id()),
            "borrowed zero-correction carry aliases the source"
        );
        assert!(
            output
                .iter()
                .all(|output_bit| output_bit.id() != borrowed_carry.id()),
            "borrowed zero-correction carry aliases the output"
        );
    }
    trace_dynamic_bit_length_zero_allocations(1, 0, 0);
    let zero = circ.alloc_qreg("rs.dynamic-bitlen.zero");
    let dirty_correction = std::env::var("LOWQ_RS_DIRTY_ZERO_CORRECTION")
        .ok()
        .as_deref()
        == Some("1")
        && source.len() >= output.len().saturating_sub(2);
    let reuse_carries = !dirty_correction
        && std::env::var("LOWQ_REUSE_ZERO_CARRIES_FOR_PREFIX")
            .ok()
            .as_deref()
            == Some("1");
    let full_prefix_scratch = full_zero_carry_prefix_scratch_requested();
    if full_prefix_scratch {
        assert!(
            reuse_carries,
            "full zero-carry prefix scratch requires LOWQ_REUSE_ZERO_CARRIES_FOR_PREFIX=1"
        );
        assert_eq!(
            std::env::var("LOWQ_DIRECT_PREFIX_BITLEN").ok().as_deref(),
            Some("1"),
            "full zero-carry prefix scratch requires the direct-prefix route"
        );
        assert_ne!(
            std::env::var("LOWQ_DIRECT_PREFIX_DIRTY_UPDATE")
                .ok()
                .as_deref(),
            Some("1"),
            "full zero-carry prefix scratch forbids dirty prefix updates"
        );
        assert_ne!(
            std::env::var("LOWQ_DIRECT_PREFIX_NO_FLAG").ok().as_deref(),
            Some("1"),
            "full zero-carry prefix scratch requires a materialized flag"
        );
        assert_ne!(
            std::env::var("LOWQ_RS_DIRTY_ZERO_CORRECTION")
                .ok()
                .as_deref(),
            Some("1"),
            "full zero-carry prefix scratch forbids dirty zero correction"
        );
        assert!(
            output.len().saturating_sub(1) <= DIRECT_PREFIX_FULL_SCRATCH_LEN,
            "full zero-carry prefix scratch supports outputs of at most ten bits"
        );
    }
    let carries = if dirty_correction {
        Vec::new()
    } else {
        let carry_count = if full_prefix_scratch {
            DIRECT_PREFIX_FULL_SCRATCH_LEN
        } else {
            output.len().saturating_sub(1)
        };
        let borrowed_count = usize::from(borrowed_carry.is_some() && carry_count != 0);
        let owned_count = carry_count - borrowed_count;
        trace_dynamic_bit_length_zero_allocations(0, owned_count, 0);
        circ.alloc_qreg_bits("rs.dynamic-bitlen.zero-carries", owned_count)
    };
    let output_refs: Vec<&QReg> = output.iter().collect();
    let mut carry_refs: Vec<&QReg> = carries.iter().collect();
    // Keep caller-owned storage at the tail. The zero correction consumes it
    // before the bit-length update, which may then reuse the same clean lane.
    if !dirty_correction && (full_prefix_scratch || output.len() > 1) {
        if let Some(borrowed_carry) = borrowed_carry {
            carry_refs.push(borrowed_carry);
        }
    }
    if full_prefix_scratch {
        assert_full_zero_carry_scratch(source, output, &zero, &carry_refs);
    }
    let zero_prefix_scratch =
        full_prefix_scratch.then(|| &carry_refs[..DIRECT_PREFIX_KG_SCRATCH_LEN]);
    toggle_static_zero_flag(circ, source, &zero, zero_prefix_scratch);
    if decrement {
        if dirty_correction {
            dirty_controlled_inc_suffix(circ, &[&zero], &output_refs, 0, false, source);
        } else {
            controlled_increment_mod_2n_carry_refs(circ, &zero, output, &carry_refs);
        }
        if full_prefix_scratch {
            bit_length_lean_with_full_prefix_scratch(circ, source, output, true, &carry_refs);
        } else if reuse_carries {
            bit_length_lean_with_increment_scratch(circ, source, output, true, &carry_refs);
        } else {
            bit_length_lean(circ, source, output, true);
        }
    } else {
        if full_prefix_scratch {
            bit_length_lean_with_full_prefix_scratch(circ, source, output, false, &carry_refs);
        } else if reuse_carries {
            bit_length_lean_with_increment_scratch(circ, source, output, false, &carry_refs);
        } else {
            bit_length_lean(circ, source, output, false);
        }
        if dirty_correction {
            dirty_controlled_inc_suffix(circ, &[&zero], &output_refs, 0, true, source);
        } else {
            controlled_decrement_mod_2n_carry_refs(circ, &zero, output, &carry_refs);
        }
    }
    toggle_static_zero_flag(circ, source, &zero, zero_prefix_scratch);
    for lane in carries {
        circ.zero_and_free(lane);
    }
    circ.zero_and_free(zero);
}

fn bit_length_lean_allow_zero_legacy(
    circ: &mut Circuit,
    source: &[&QReg],
    output: &[QReg],
    decrement: bool,
) {
    bit_length_lean_allow_zero_legacy_with_borrowed_carry(circ, source, output, decrement, None);
}

fn fused_prefix_scratch_loan_requested() -> bool {
    std::env::var(FUSED_PREFIX_SCRATCH_LOAN_FLAG)
        .ok()
        .as_deref()
        == Some("1")
}

fn preserved_dy_top_prefix_loan_requested() -> bool {
    std::env::var(PRESERVED_DY_TOP_PREFIX_LOAN_FLAG)
        .ok()
        .as_deref()
        == Some("1")
}

fn mixed_width_l_r_prime_requested() -> bool {
    std::env::var(MIXED_WIDTH_L_R_PRIME_FLAG).ok().as_deref() == Some("1")
}

fn production_l_r_prime_width() -> usize {
    if mixed_width_l_r_prime_requested() {
        REFERENCE_R_LENGTH_WIDTH
    } else {
        REFERENCE_LENGTH_WIDTH
    }
}

fn assert_l_r_prime_metadata_width(full_width: usize, r_width: usize) {
    assert!(r_width == full_width || r_width + 1 == full_width);
}

fn assert_fused_prefix_scratch_lenders(source: &[&QReg], output: &[QReg], lenders: &[&QReg]) {
    use super::shrunken_pz_state_machine::DIRECT_PREFIX_FULL_SCRATCH_LEN;

    assert!(
        lenders.len() <= DIRECT_PREFIX_FULL_SCRATCH_LEN,
        "fused-prefix scratch loan provides {} lanes but the layout has only {DIRECT_PREFIX_FULL_SCRATCH_LEN}",
        lenders.len()
    );
    for (index, lane) in lenders.iter().enumerate() {
        assert!(
            lenders[..index].iter().all(|other| other.id() != lane.id()),
            "fused-prefix scratch lender {index} aliases an earlier lender"
        );
        assert!(
            source.iter().all(|other| other.id() != lane.id()),
            "fused-prefix scratch lender {index} aliases the source"
        );
        assert!(
            output.iter().all(|other| other.id() != lane.id()),
            "fused-prefix scratch lender {index} aliases the output"
        );
    }
}

fn bit_length_lean_allow_zero_with_borrowed_scratch_impl(
    circ: &mut Circuit,
    source: &[&QReg],
    output: &[QReg],
    decrement: bool,
    borrowed_carry: Option<&QReg>,
    borrowed_prefix_scratch: &[&QReg],
    source_is_complemented: bool,
    omitted_high_bits: usize,
) {
    use super::shrunken_pz_state_machine::{
        bit_length_lean, bit_length_lean_complemented_source,
        bit_length_lean_with_full_prefix_scratch,
        bit_length_lean_with_full_prefix_scratch_complemented_source,
        bit_length_lean_with_full_prefix_scratch_split_four_high,
        bit_length_lean_with_full_prefix_scratch_split_high,
        bit_length_lean_with_full_prefix_scratch_split_three_high,
        bit_length_lean_with_full_prefix_scratch_split_two_high,
        lowq_fused_zero_prefix_bitlen_requested, DIRECT_PREFIX_FULL_SCRATCH_LEN,
    };

    let prefix_loan = fused_prefix_scratch_loan_requested();
    if prefix_loan {
        assert!(
            lowq_fused_zero_prefix_bitlen_requested(),
            "fused-prefix scratch loan requires LOWQ_FUSED_ZERO_PREFIX_BITLEN=1"
        );
        assert!(
            full_zero_carry_prefix_scratch_requested(),
            "fused-prefix scratch loan requires LOWQ_REUSE_ZERO_CARRIES_FOR_FULL_PREFIX_SCRATCH=1"
        );
        assert_fused_prefix_scratch_lenders(source, output, borrowed_prefix_scratch);
    }
    if !lowq_fused_zero_prefix_bitlen_requested() {
        assert!(
            !source_is_complemented,
            "pre-complemented source requires fused zero-prefix bit length"
        );
        bit_length_lean_allow_zero_legacy_with_borrowed_carry(
            circ,
            source,
            output,
            decrement,
            borrowed_carry,
        );
        return;
    }
    if let Some(borrowed_carry) = borrowed_carry {
        assert!(
            source
                .iter()
                .all(|source_bit| source_bit.id() != borrowed_carry.id()),
            "borrowed zero-correction carry aliases the source"
        );
        assert!(
            output
                .iter()
                .all(|output_bit| output_bit.id() != borrowed_carry.id()),
            "borrowed zero-correction carry aliases the output"
        );
    }
    assert_eq!(
        std::env::var("LOWQ_DIRECT_PREFIX_BITLEN").ok().as_deref(),
        Some("1"),
        "LOWQ_FUSED_ZERO_PREFIX_BITLEN requires LOWQ_DIRECT_PREFIX_BITLEN=1"
    );
    if source.is_empty() {
        return;
    }
    assert!(
        omitted_high_bits == 0 || full_zero_carry_prefix_scratch_requested(),
        "split-high bit length requires full prefix scratch"
    );
    if full_zero_carry_prefix_scratch_requested() {
        assert_eq!(
            std::env::var("LOWQ_REUSE_ZERO_CARRIES_FOR_PREFIX")
                .ok()
                .as_deref(),
            Some("1"),
            "fused full prefix scratch requires LOWQ_REUSE_ZERO_CARRIES_FOR_PREFIX=1"
        );
        assert!(
            output.len().saturating_sub(1) <= DIRECT_PREFIX_FULL_SCRATCH_LEN,
            "fused full prefix scratch supports outputs of at most ten bits"
        );
        let borrowed_lanes = if prefix_loan {
            borrowed_prefix_scratch.len()
        } else {
            0
        };
        let owned_lanes = DIRECT_PREFIX_FULL_SCRATCH_LEN - borrowed_lanes;
        trace_fused_prefix_scratch_loan_allocation(owned_lanes, borrowed_lanes);
        let scratch = circ.alloc_qreg_bits("rs.dynamic-bitlen.fused-prefix-scratch", owned_lanes);
        let scratch_refs: Vec<&QReg> = scratch
            .iter()
            .chain(borrowed_prefix_scratch.iter().copied().take(borrowed_lanes))
            .collect();
        assert_eq!(scratch_refs.len(), DIRECT_PREFIX_FULL_SCRATCH_LEN);
        if omitted_high_bits == 1 {
            bit_length_lean_with_full_prefix_scratch_split_high(
                circ,
                source,
                output,
                decrement,
                &scratch_refs,
                source_is_complemented,
            );
        } else if omitted_high_bits == 2 {
            bit_length_lean_with_full_prefix_scratch_split_two_high(
                circ,
                source,
                output,
                decrement,
                &scratch_refs,
                source_is_complemented,
            );
        } else if omitted_high_bits == 3 {
            bit_length_lean_with_full_prefix_scratch_split_three_high(
                circ,
                source,
                output,
                decrement,
                &scratch_refs,
                source_is_complemented,
            );
        } else if omitted_high_bits == 4 {
            bit_length_lean_with_full_prefix_scratch_split_four_high(
                circ,
                source,
                output,
                decrement,
                &scratch_refs,
                source_is_complemented,
            );
        } else if source_is_complemented {
            bit_length_lean_with_full_prefix_scratch_complemented_source(
                circ,
                source,
                output,
                decrement,
                &scratch_refs,
            );
        } else {
            bit_length_lean_with_full_prefix_scratch(
                circ,
                source,
                output,
                decrement,
                &scratch_refs,
            );
        }
        for lane in scratch {
            circ.zero_and_free(lane);
        }
    } else if source_is_complemented {
        bit_length_lean_complemented_source(circ, source, output, decrement);
    } else {
        bit_length_lean(circ, source, output, decrement);
    }
}

fn bit_length_lean_allow_zero_with_borrowed_scratch(
    circ: &mut Circuit,
    source: &[&QReg],
    output: &[QReg],
    decrement: bool,
    borrowed_carry: Option<&QReg>,
    borrowed_prefix_scratch: &[&QReg],
) {
    bit_length_lean_allow_zero_with_borrowed_scratch_impl(
        circ,
        source,
        output,
        decrement,
        borrowed_carry,
        borrowed_prefix_scratch,
        false,
        0,
    );
}

fn bit_length_lean_allow_zero_with_borrowed_scratch_complemented_source(
    circ: &mut Circuit,
    source: &[&QReg],
    output: &[QReg],
    decrement: bool,
    borrowed_carry: Option<&QReg>,
    borrowed_prefix_scratch: &[&QReg],
) {
    bit_length_lean_allow_zero_with_borrowed_scratch_impl(
        circ,
        source,
        output,
        decrement,
        borrowed_carry,
        borrowed_prefix_scratch,
        true,
        0,
    );
}

fn bit_length_lean_allow_zero_with_borrowed_scratch_split_high(
    circ: &mut Circuit,
    source: &[&QReg],
    output: &[QReg],
    decrement: bool,
    borrowed_carry: Option<&QReg>,
    borrowed_prefix_scratch: &[&QReg],
    source_is_complemented: bool,
) {
    bit_length_lean_allow_zero_with_borrowed_scratch_impl(
        circ,
        source,
        output,
        decrement,
        borrowed_carry,
        borrowed_prefix_scratch,
        source_is_complemented,
        1,
    );
}

fn bit_length_lean_allow_zero_with_borrowed_scratch_split_two_high(
    circ: &mut Circuit,
    source: &[&QReg],
    output: &[QReg],
    decrement: bool,
    borrowed_carry: Option<&QReg>,
    borrowed_prefix_scratch: &[&QReg],
    source_is_complemented: bool,
) {
    bit_length_lean_allow_zero_with_borrowed_scratch_impl(
        circ,
        source,
        output,
        decrement,
        borrowed_carry,
        borrowed_prefix_scratch,
        source_is_complemented,
        2,
    );
}

fn bit_length_lean_allow_zero_with_borrowed_scratch_split_three_high(
    circ: &mut Circuit,
    source: &[&QReg],
    output: &[QReg],
    decrement: bool,
    borrowed_carry: Option<&QReg>,
    borrowed_prefix_scratch: &[&QReg],
    source_is_complemented: bool,
) {
    bit_length_lean_allow_zero_with_borrowed_scratch_impl(
        circ,
        source,
        output,
        decrement,
        borrowed_carry,
        borrowed_prefix_scratch,
        source_is_complemented,
        3,
    );
}

fn bit_length_lean_allow_zero_with_borrowed_scratch_split_four_high(
    circ: &mut Circuit,
    source: &[&QReg],
    output: &[QReg],
    decrement: bool,
    borrowed_carry: Option<&QReg>,
    borrowed_prefix_scratch: &[&QReg],
    source_is_complemented: bool,
) {
    bit_length_lean_allow_zero_with_borrowed_scratch_impl(
        circ,
        source,
        output,
        decrement,
        borrowed_carry,
        borrowed_prefix_scratch,
        source_is_complemented,
        4,
    );
}

fn bit_length_lean_allow_zero_with_borrowed_carry(
    circ: &mut Circuit,
    source: &[&QReg],
    output: &[QReg],
    decrement: bool,
    borrowed_carry: Option<&QReg>,
) {
    bit_length_lean_allow_zero_with_borrowed_scratch(
        circ,
        source,
        output,
        decrement,
        borrowed_carry,
        &[],
    );
}

fn bit_length_lean_allow_zero(
    circ: &mut Circuit,
    source: &[&QReg],
    output: &[QReg],
    decrement: bool,
) {
    bit_length_lean_allow_zero_with_borrowed_carry(circ, source, output, decrement, None);
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct DirtyZeroBitLengthProofReport {
    pub cases_checked: usize,
    pub directions_checked: usize,
    pub maximum_extra_qubits: usize,
    pub maximum_emitted_ops: usize,
    pub maximum_emitted_toffoli: usize,
}

/// Exhaustively verify the dirty zero-correction composition for every
/// eight-bit source (including zero), every five-bit accumulator, and both
/// update directions. The caller enables the direct-prefix route before entry.
fn zero_bit_length_roundtrip_check_mode(
    dirty_correction: bool,
    reuse_carries: bool,
    full_prefix_scratch: bool,
) -> DirtyZeroBitLengthProofReport {
    use crate::circuit::{OperationType, QubitId};
    use crate::sim::Simulator;
    use sha3::{
        digest::{ExtendableOutput, Update},
        Shake128,
    };

    assert_eq!(
        std::env::var("LOWQ_DIRECT_PREFIX_BITLEN").ok().as_deref(),
        Some("1")
    );
    assert_eq!(
        std::env::var("LOWQ_DIRECT_PREFIX_DIRTY_UPDATE")
            .ok()
            .as_deref(),
        dirty_correction.then_some("1")
    );
    assert_eq!(
        std::env::var("LOWQ_DIRECT_PREFIX_NO_FLAG").ok().as_deref(),
        dirty_correction.then_some("1")
    );
    assert_eq!(
        std::env::var("LOWQ_RS_DIRTY_ZERO_CORRECTION")
            .ok()
            .as_deref(),
        dirty_correction.then_some("1")
    );
    assert_eq!(
        std::env::var("LOWQ_REUSE_ZERO_CARRIES_FOR_PREFIX")
            .ok()
            .as_deref(),
        reuse_carries.then_some("1")
    );
    assert_eq!(
        std::env::var("LOWQ_REUSE_ZERO_CARRIES_FOR_FULL_PREFIX_SCRATCH")
            .ok()
            .as_deref(),
        full_prefix_scratch.then_some("1")
    );

    let mut cases_checked = 0usize;
    let mut maximum_extra_qubits = 0usize;
    let mut maximum_emitted_ops = 0usize;
    let mut maximum_emitted_toffoli = 0usize;
    for decrement in [false, true] {
        let mut circuit = Circuit::new();
        let source = circuit.alloc_qreg_bits("dirty-zero-proof.source", 8);
        let output = circuit.alloc_qreg_bits("dirty-zero-proof.output", 5);
        let source_refs: Vec<&QReg> = source.iter().collect();
        bit_length_lean_allow_zero(&mut circuit, &source_refs, &output, decrement);

        let source_ids: Vec<u32> = source.iter().map(QReg::id).collect();
        let output_ids: Vec<u32> = output.iter().map(QReg::id).collect();
        let external: Vec<u32> = source_ids
            .iter()
            .chain(output_ids.iter())
            .copied()
            .collect();
        let builder = circuit.into_builder();
        maximum_extra_qubits =
            maximum_extra_qubits.max(builder.peak_qubits as usize - external.len());
        maximum_emitted_ops = maximum_emitted_ops.max(builder.ops.len());
        maximum_emitted_toffoli = maximum_emitted_toffoli.max(
            builder
                .ops
                .iter()
                .filter(|op| matches!(op.kind, OperationType::CCX | OperationType::CCZ))
                .count(),
        );

        let cases: Vec<(u64, u64)> = (0u64..=255)
            .flat_map(|source_value| {
                (0u64..32).map(move |output_value| (source_value, output_value))
            })
            .collect();
        for (batch, chunk) in cases.chunks(64).enumerate() {
            let mut seed = Shake128::default();
            seed.update(if decrement {
                b"dirty-zero-bitlen-sub"
            } else {
                b"dirty-zero-bitlen-add"
            });
            seed.update(&(batch as u64).to_le_bytes());
            let mut xof = seed.finalize_xof();
            let mut simulator = Simulator::new(
                builder.next_qubit as usize,
                builder.next_bit as usize,
                &mut xof,
            );
            for (shot, &(source_value, output_value)) in chunk.iter().enumerate() {
                for (bit, &id) in source_ids.iter().enumerate() {
                    if (source_value >> bit) & 1 == 1 {
                        *simulator.qubit_mut(QubitId(u64::from(id))) |= 1u64 << shot;
                    }
                }
                for (bit, &id) in output_ids.iter().enumerate() {
                    if (output_value >> bit) & 1 == 1 {
                        *simulator.qubit_mut(QubitId(u64::from(id))) |= 1u64 << shot;
                    }
                }
            }
            simulator.apply_iter(builder.ops.iter());
            let live = if chunk.len() == 64 {
                u64::MAX
            } else {
                (1u64 << chunk.len()) - 1
            };
            assert_eq!(simulator.phase & live, 0, "phase failure in batch {batch}");

            for (shot, &(source_value, output_value)) in chunk.iter().enumerate() {
                let bit_length = if source_value == 0 {
                    0
                } else {
                    64 - source_value.leading_zeros() as u64
                };
                let expected = if decrement {
                    output_value.wrapping_sub(bit_length) & 31
                } else {
                    output_value.wrapping_add(bit_length) & 31
                };
                let read = |ids: &[u32]| {
                    ids.iter().enumerate().fold(0u64, |value, (bit, &id)| {
                        value | (((simulator.qubit(QubitId(u64::from(id))) >> shot) & 1) << bit)
                    })
                };
                assert_eq!(read(&source_ids), source_value);
                assert_eq!(read(&output_ids), expected);
            }
            for id in 0..builder.next_qubit {
                if !external.contains(&id) {
                    assert_eq!(
                        simulator.qubit(QubitId(u64::from(id))) & live,
                        0,
                        "ancilla q{id} dirty in batch {batch}"
                    );
                }
            }
            cases_checked += chunk.len();
        }
    }

    DirtyZeroBitLengthProofReport {
        cases_checked,
        directions_checked: 2,
        maximum_extra_qubits,
        maximum_emitted_ops,
        maximum_emitted_toffoli,
    }
}

#[doc(hidden)]
pub fn dirty_zero_bit_length_roundtrip_check() -> DirtyZeroBitLengthProofReport {
    zero_bit_length_roundtrip_check_mode(true, false, false)
}

#[doc(hidden)]
pub fn reused_zero_carry_bit_length_roundtrip_check() -> DirtyZeroBitLengthProofReport {
    zero_bit_length_roundtrip_check_mode(false, true, false)
}

#[doc(hidden)]
pub fn fully_reused_zero_carry_bit_length_roundtrip_check() -> DirtyZeroBitLengthProofReport {
    zero_bit_length_roundtrip_check_mode(false, true, true)
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
    pub accumulator_values_per_source: usize,
    pub update_cases_checked: usize,
    pub controlled_cases_checked: usize,
    pub directions_checked: usize,
    pub control_values_checked: usize,
    pub baseline_equivalence_checks: usize,
    pub default_stream_equivalence_checks: usize,
    pub phase_clean_checks: usize,
    pub ancilla_clean_checks: usize,
    pub precondition_rejections: usize,
    pub trace_sensitivity_flag_allocations: usize,
    pub trace_sensitivity_carry_allocations: usize,
    pub trace_sensitivity_prefix_allocations: usize,
    pub fused_zero_flag_allocations: usize,
    pub fused_zero_carry_allocations: usize,
    pub fused_zero_prefix_allocations: usize,
    pub maximum_baseline_extra_qubits: usize,
    pub maximum_fused_extra_qubits: usize,
    pub maximum_baseline_emitted_ops: usize,
    pub maximum_fused_emitted_ops: usize,
    pub maximum_baseline_emitted_toffoli: usize,
    pub maximum_fused_emitted_toffoli: usize,
    pub width_direction_toffoli_increase_cases: usize,
    pub maximum_width_direction_toffoli_increase: usize,
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

fn configure_fused_zero_prefix_proof(full_scratch_baseline: bool, fused: bool) {
    std::env::set_var("LOWQ_DIRECT_PREFIX_BITLEN", "1");
    for name in [
        "LOWQ_DIRECT_PREFIX_DIRTY_UPDATE",
        "LOWQ_DIRECT_PREFIX_NO_FLAG",
        "LOWQ_RS_DIRTY_ZERO_CORRECTION",
    ] {
        std::env::remove_var(name);
    }
    if full_scratch_baseline {
        std::env::set_var("LOWQ_REUSE_ZERO_CARRIES_FOR_PREFIX", "1");
        std::env::set_var("LOWQ_REUSE_ZERO_CARRIES_FOR_FULL_PREFIX_SCRATCH", "1");
    } else {
        std::env::remove_var("LOWQ_REUSE_ZERO_CARRIES_FOR_PREFIX");
        std::env::remove_var("LOWQ_REUSE_ZERO_CARRIES_FOR_FULL_PREFIX_SCRATCH");
    }
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
    let accumulator = circuit.alloc_qreg_bits("fused-zero-proof.accumulator", accumulator_width);
    let source_refs: Vec<&QReg> = source.iter().collect();
    if legacy_entry {
        bit_length_lean_allow_zero_legacy(&mut circuit, &source_refs, &accumulator, decrement);
    } else {
        bit_length_lean_allow_zero(&mut circuit, &source_refs, &accumulator, decrement);
    }
    let source_ids = source.iter().map(QReg::id).collect();
    let accumulator_ids = accumulator.iter().map(QReg::id).collect();
    FusedZeroPrefixProofCircuit {
        builder: circuit.into_builder(),
        source_ids,
        accumulator_ids,
        control_id: None,
    }
}

fn build_fused_zero_prefix_controlled_xor(source_width: usize) -> FusedZeroPrefixProofCircuit {
    const ACCUMULATOR_WIDTH: usize = 5;

    let mut circuit = Circuit::new();
    let source = circuit.alloc_qreg_bits("fused-zero-control-proof.source", source_width);
    let control = circuit.alloc_qreg("fused-zero-control-proof.control");
    let accumulator =
        circuit.alloc_qreg_bits("fused-zero-control-proof.accumulator", ACCUMULATOR_WIDTH);
    let temporary =
        circuit.alloc_qreg_bits("fused-zero-control-proof.temporary", ACCUMULATOR_WIDTH);
    let source_refs: Vec<&QReg> = source.iter().collect();
    bit_length_lean_allow_zero(&mut circuit, &source_refs, &temporary, false);
    for (length_bit, accumulator_bit) in temporary.iter().zip(&accumulator) {
        circuit.ccx(&control, length_bit, accumulator_bit);
    }
    bit_length_lean_allow_zero(&mut circuit, &source_refs, &temporary, true);
    for lane in temporary {
        circuit.zero_and_free(lane);
    }
    let source_ids = source.iter().map(QReg::id).collect();
    let accumulator_ids = accumulator.iter().map(QReg::id).collect();
    FusedZeroPrefixProofCircuit {
        builder: circuit.into_builder(),
        source_ids,
        accumulator_ids,
        control_id: Some(control.id()),
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
            .filter(|operation| matches!(operation.kind, OperationType::CCX | OperationType::CCZ))
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

fn assert_zero_dynamic_bit_length_trace(trace: DynamicBitLengthZeroAllocationTrace) {
    assert_eq!(
        trace,
        DynamicBitLengthZeroAllocationTrace::default(),
        "fused traversal allocated an rs.dynamic-bitlen.zero* lane"
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
        baseline_seed.update(b"fused-zero-prefix-update-baseline");
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
        fused_seed.update(b"fused-zero-prefix-update-fused");
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
                    *baseline_simulator.qubit_mut(QubitId(u64::from(baseline_id))) |= 1u64 << shot;
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
                    *baseline_simulator.qubit_mut(QubitId(u64::from(baseline_id))) |= 1u64 << shot;
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
        assert_eq!(
            baseline_simulator.phase & live,
            0,
            "baseline phase failure at width {source_width}, batch {batch}"
        );
        assert_eq!(
            fused_simulator.phase & live,
            0,
            "fused phase failure at width {source_width}, batch {batch}"
        );
        phase_clean_checks += 2 * chunk.len();

        for (shot, &(source_value, accumulator_value)) in chunk.iter().enumerate() {
            let bit_length = if source_value == 0 {
                0
            } else {
                64 - u64::from(source_value.leading_zeros())
            };
            let expected = if decrement {
                accumulator_value.wrapping_sub(bit_length) % accumulator_modulus
            } else {
                accumulator_value.wrapping_add(bit_length) % accumulator_modulus
            };
            let baseline_source =
                read_fused_zero_prefix_register(&baseline_simulator, &baseline.source_ids, shot);
            let fused_source =
                read_fused_zero_prefix_register(&fused_simulator, &fused.source_ids, shot);
            let baseline_accumulator = read_fused_zero_prefix_register(
                &baseline_simulator,
                &baseline.accumulator_ids,
                shot,
            );
            let fused_accumulator =
                read_fused_zero_prefix_register(&fused_simulator, &fused.accumulator_ids, shot);
            assert_eq!(baseline_source, source_value);
            assert_eq!(fused_source, source_value);
            assert_eq!(baseline_accumulator, expected);
            assert_eq!(fused_accumulator, expected);
            assert_eq!(fused_source, baseline_source);
            assert_eq!(fused_accumulator, baseline_accumulator);
        }
        ancilla_clean_checks += assert_fused_zero_prefix_internal_clean(
            &baseline_simulator,
            baseline,
            live,
            chunk.len(),
            "baseline update",
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
    let accumulator_modulus = 1u64 << fused.accumulator_ids.len();
    let cases: Vec<(u64, u64, u64)> = (0..(1u64 << source_width))
        .flat_map(|source| {
            (0..=1u64).flat_map(move |control| {
                (0..accumulator_modulus).map(move |accumulator| (source, control, accumulator))
            })
        })
        .collect();
    let mut phase_clean_checks = 0usize;
    let mut ancilla_clean_checks = 0usize;

    for (batch, chunk) in cases.chunks(64).enumerate() {
        let mut baseline_seed = Shake128::default();
        baseline_seed.update(b"fused-zero-prefix-control-baseline");
        baseline_seed.update(&(source_width as u64).to_le_bytes());
        baseline_seed.update(&(batch as u64).to_le_bytes());
        let mut baseline_xof = baseline_seed.finalize_xof();
        let mut baseline_simulator = Simulator::new(
            baseline.builder.next_qubit as usize,
            baseline.builder.next_bit as usize,
            &mut baseline_xof,
        );

        let mut fused_seed = Shake128::default();
        fused_seed.update(b"fused-zero-prefix-control-fused");
        fused_seed.update(&(source_width as u64).to_le_bytes());
        fused_seed.update(&(batch as u64).to_le_bytes());
        let mut fused_xof = fused_seed.finalize_xof();
        let mut fused_simulator = Simulator::new(
            fused.builder.next_qubit as usize,
            fused.builder.next_bit as usize,
            &mut fused_xof,
        );

        for (shot, &(source_value, control_value, accumulator_value)) in chunk.iter().enumerate() {
            for (bit, (&baseline_id, &fused_id)) in baseline
                .source_ids
                .iter()
                .zip(&fused.source_ids)
                .enumerate()
            {
                if (source_value >> bit) & 1 == 1 {
                    *baseline_simulator.qubit_mut(QubitId(u64::from(baseline_id))) |= 1u64 << shot;
                    *fused_simulator.qubit_mut(QubitId(u64::from(fused_id))) |= 1u64 << shot;
                }
            }
            if control_value == 1 {
                *baseline_simulator.qubit_mut(QubitId(u64::from(baseline_control))) |= 1u64 << shot;
                *fused_simulator.qubit_mut(QubitId(u64::from(fused_control))) |= 1u64 << shot;
            }
            for (bit, (&baseline_id, &fused_id)) in baseline
                .accumulator_ids
                .iter()
                .zip(&fused.accumulator_ids)
                .enumerate()
            {
                if (accumulator_value >> bit) & 1 == 1 {
                    *baseline_simulator.qubit_mut(QubitId(u64::from(baseline_id))) |= 1u64 << shot;
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
        assert_eq!(
            baseline_simulator.phase & live,
            0,
            "baseline controlled phase failure at width {source_width}, batch {batch}"
        );
        assert_eq!(
            fused_simulator.phase & live,
            0,
            "fused controlled phase failure at width {source_width}, batch {batch}"
        );
        phase_clean_checks += 2 * chunk.len();

        for (shot, &(source_value, control_value, accumulator_value)) in chunk.iter().enumerate() {
            let bit_length = if source_value == 0 {
                0
            } else {
                64 - u64::from(source_value.leading_zeros())
            };
            let expected = accumulator_value ^ if control_value == 1 { bit_length } else { 0 };
            let baseline_source =
                read_fused_zero_prefix_register(&baseline_simulator, &baseline.source_ids, shot);
            let fused_source =
                read_fused_zero_prefix_register(&fused_simulator, &fused.source_ids, shot);
            let baseline_accumulator = read_fused_zero_prefix_register(
                &baseline_simulator,
                &baseline.accumulator_ids,
                shot,
            );
            let fused_accumulator =
                read_fused_zero_prefix_register(&fused_simulator, &fused.accumulator_ids, shot);
            assert_eq!(baseline_source, source_value);
            assert_eq!(fused_source, source_value);
            assert_eq!(
                (baseline_simulator.qubit(QubitId(u64::from(baseline_control))) >> shot) & 1,
                control_value
            );
            assert_eq!(
                (fused_simulator.qubit(QubitId(u64::from(fused_control))) >> shot) & 1,
                control_value
            );
            assert_eq!(baseline_accumulator, expected);
            assert_eq!(fused_accumulator, expected);
            assert_eq!(fused_source, baseline_source);
            assert_eq!(fused_accumulator, baseline_accumulator);
        }
        ancilla_clean_checks += assert_fused_zero_prefix_internal_clean(
            &baseline_simulator,
            baseline,
            live,
            chunk.len(),
            "baseline controlled context",
        );
        ancilla_clean_checks += assert_fused_zero_prefix_internal_clean(
            &fused_simulator,
            fused,
            live,
            chunk.len(),
            "fused controlled context",
        );
    }

    (cases.len(), phase_clean_checks, ancilla_clean_checks)
}

fn assert_fused_zero_prefix_rejection(
    result: Result<(), Box<dyn std::any::Any + Send>>,
    expected: &str,
) {
    let payload = result.expect_err("invalid fused zero-prefix call was accepted");
    let message = payload
        .downcast_ref::<String>()
        .map(String::as_str)
        .or_else(|| payload.downcast_ref::<&str>().copied())
        .unwrap_or("non-string panic payload");
    assert!(
        message.contains(expected),
        "unexpected fused zero-prefix rejection: {message}"
    );
}

fn fused_zero_prefix_precondition_rejections() -> usize {
    use super::shrunken_pz_state_machine::{
        bit_length_lean, bit_length_lean_with_increment_scratch,
    };
    use std::panic::{catch_unwind, AssertUnwindSafe};

    configure_fused_zero_prefix_proof(true, true);
    let previous_hook = std::panic::take_hook();
    std::panic::set_hook(Box::new(|_| {}));

    let source_target_alias = catch_unwind(AssertUnwindSafe(|| {
        let mut circuit = Circuit::new();
        let source = circuit.alloc_qreg_bits("fused-zero-reject.alias", 8);
        let source_refs: Vec<&QReg> = source.iter().collect();
        bit_length_lean(&mut circuit, &source_refs, &source, false);
    }));
    let duplicate_source = catch_unwind(AssertUnwindSafe(|| {
        let mut circuit = Circuit::new();
        let source = circuit.alloc_qreg_bits("fused-zero-reject.duplicate-source", 2);
        let output = circuit.alloc_qreg_bits("fused-zero-reject.output", 5);
        let source_refs = vec![&source[0], &source[0]];
        bit_length_lean(&mut circuit, &source_refs, &output, false);
    }));
    let narrow_target = catch_unwind(AssertUnwindSafe(|| {
        let mut circuit = Circuit::new();
        let source = circuit.alloc_qreg_bits("fused-zero-reject.wide-source", 8);
        let output = circuit.alloc_qreg_bits("fused-zero-reject.narrow-target", 3);
        let source_refs: Vec<&QReg> = source.iter().collect();
        bit_length_lean(&mut circuit, &source_refs, &output, false);
    }));
    let scratch_alias = catch_unwind(AssertUnwindSafe(|| {
        let mut circuit = Circuit::new();
        let source = circuit.alloc_qreg_bits("fused-zero-reject.scratch-source", 8);
        let output = circuit.alloc_qreg_bits("fused-zero-reject.scratch-output", 5);
        let source_refs: Vec<&QReg> = source.iter().collect();
        let scratch = vec![&source[0]];
        bit_length_lean_with_increment_scratch(
            &mut circuit,
            &source_refs,
            &output,
            false,
            &scratch,
        );
    }));
    std::env::remove_var("LOWQ_DIRECT_PREFIX_BITLEN");
    let missing_direct_route = catch_unwind(AssertUnwindSafe(|| {
        let mut circuit = Circuit::new();
        let source = circuit.alloc_qreg_bits("fused-zero-reject.no-direct-source", 1);
        let output = circuit.alloc_qreg_bits("fused-zero-reject.no-direct-output", 1);
        let source_refs: Vec<&QReg> = source.iter().collect();
        bit_length_lean_allow_zero(&mut circuit, &source_refs, &output, false);
    }));

    std::panic::set_hook(previous_hook);
    assert_fused_zero_prefix_rejection(source_target_alias, "aliases the target");
    assert_fused_zero_prefix_rejection(duplicate_source, "aliases an earlier source lane");
    assert_fused_zero_prefix_rejection(narrow_target, "cannot represent source width");
    assert_fused_zero_prefix_rejection(scratch_alias, "scratch lane 0 aliases the source");
    assert_fused_zero_prefix_rejection(
        missing_direct_route,
        "requires LOWQ_DIRECT_PREFIX_BITLEN=1",
    );
    configure_fused_zero_prefix_proof(true, true);
    5
}

/// Exhaustively prove the feature-gated `i=n` prefix fusion against the
/// historical zero-flag composition. Widths zero through eight cover every
/// source, every five-bit accumulator, both update directions, and both values
/// of an external control in the compute/XOR/uncompute usage.
#[doc(hidden)]
pub fn fused_zero_prefix_bit_length_roundtrip_check() -> FusedZeroPrefixBitLengthProofReport {
    const MAX_SOURCE_WIDTH: usize = 8;
    const ACCUMULATOR_WIDTH: usize = 5;
    const LOCAL_WIDTH: usize = 259;
    const LOCAL_ACCUMULATOR_WIDTH: usize = 10;

    configure_fused_zero_prefix_proof(false, false);
    begin_dynamic_bit_length_zero_allocation_trace();
    let _trace_sensitivity =
        build_fused_zero_prefix_update(MAX_SOURCE_WIDTH, ACCUMULATOR_WIDTH, false, false);
    let trace_sensitivity = finish_dynamic_bit_length_zero_allocation_trace();
    assert!(trace_sensitivity.flag_allocations > 0);
    assert!(trace_sensitivity.carry_allocations > 0);
    assert!(trace_sensitivity.prefix_allocations > 0);

    let mut update_cases_checked = 0usize;
    let mut controlled_cases_checked = 0usize;
    let mut baseline_equivalence_checks = 0usize;
    let mut default_stream_equivalence_checks = 0usize;
    let mut phase_clean_checks = 0usize;
    let mut ancilla_clean_checks = 0usize;
    let mut fused_trace = DynamicBitLengthZeroAllocationTrace::default();
    let mut maximum_baseline = FusedZeroPrefixLocalResources::default();
    let mut maximum_fused = FusedZeroPrefixLocalResources::default();
    let mut width_direction_toffoli_increase_cases = 0usize;
    let mut maximum_width_direction_toffoli_increase = 0usize;

    for source_width in 0..=MAX_SOURCE_WIDTH {
        for decrement in [false, true] {
            configure_fused_zero_prefix_proof(true, false);
            let baseline =
                build_fused_zero_prefix_update(source_width, ACCUMULATOR_WIDTH, decrement, false);
            let legacy =
                build_fused_zero_prefix_update(source_width, ACCUMULATOR_WIDTH, decrement, true);
            assert_fused_zero_prefix_default_stream(&baseline, &legacy);
            default_stream_equivalence_checks += 1;

            configure_fused_zero_prefix_proof(true, true);
            begin_dynamic_bit_length_zero_allocation_trace();
            let fused =
                build_fused_zero_prefix_update(source_width, ACCUMULATOR_WIDTH, decrement, false);
            let trace = finish_dynamic_bit_length_zero_allocation_trace();
            assert_zero_dynamic_bit_length_trace(trace);
            fused_trace.flag_allocations += trace.flag_allocations;
            fused_trace.carry_allocations += trace.carry_allocations;
            fused_trace.prefix_allocations += trace.prefix_allocations;

            let baseline_resources = fused_zero_prefix_resources(&baseline);
            let fused_resources = fused_zero_prefix_resources(&fused);
            assert!(
                fused_resources.extra_qubits <= baseline_resources.extra_qubits,
                "fused width {source_width} decrement={decrement} increased the qubit peak"
            );
            if fused_resources.emitted_toffoli > baseline_resources.emitted_toffoli {
                width_direction_toffoli_increase_cases += 1;
                maximum_width_direction_toffoli_increase = maximum_width_direction_toffoli_increase
                    .max(fused_resources.emitted_toffoli - baseline_resources.emitted_toffoli);
            }
            maximum_baseline.extra_qubits = maximum_baseline
                .extra_qubits
                .max(baseline_resources.extra_qubits);
            maximum_baseline.emitted_ops = maximum_baseline
                .emitted_ops
                .max(baseline_resources.emitted_ops);
            maximum_baseline.emitted_toffoli = maximum_baseline
                .emitted_toffoli
                .max(baseline_resources.emitted_toffoli);
            maximum_fused.extra_qubits =
                maximum_fused.extra_qubits.max(fused_resources.extra_qubits);
            maximum_fused.emitted_ops = maximum_fused.emitted_ops.max(fused_resources.emitted_ops);
            maximum_fused.emitted_toffoli = maximum_fused
                .emitted_toffoli
                .max(fused_resources.emitted_toffoli);

            let (cases, phase_checks, ancilla_checks) = verify_fused_zero_prefix_update_equivalence(
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

        configure_fused_zero_prefix_proof(true, false);
        let baseline_controlled = build_fused_zero_prefix_controlled_xor(source_width);
        configure_fused_zero_prefix_proof(true, true);
        begin_dynamic_bit_length_zero_allocation_trace();
        let fused_controlled = build_fused_zero_prefix_controlled_xor(source_width);
        let trace = finish_dynamic_bit_length_zero_allocation_trace();
        assert_zero_dynamic_bit_length_trace(trace);
        fused_trace.flag_allocations += trace.flag_allocations;
        fused_trace.carry_allocations += trace.carry_allocations;
        fused_trace.prefix_allocations += trace.prefix_allocations;
        let (cases, phase_checks, ancilla_checks) = verify_fused_zero_prefix_controlled_equivalence(
            &baseline_controlled,
            &fused_controlled,
            source_width,
        );
        controlled_cases_checked += cases;
        baseline_equivalence_checks += cases;
        phase_clean_checks += phase_checks;
        ancilla_clean_checks += ancilla_checks;
    }

    configure_fused_zero_prefix_proof(true, false);
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
    configure_fused_zero_prefix_proof(true, true);
    begin_dynamic_bit_length_zero_allocation_trace();
    let local_add_fused = fused_zero_prefix_resources(&build_fused_zero_prefix_update(
        LOCAL_WIDTH,
        LOCAL_ACCUMULATOR_WIDTH,
        false,
        false,
    ));
    let local_add_trace = finish_dynamic_bit_length_zero_allocation_trace();
    assert_zero_dynamic_bit_length_trace(local_add_trace);
    begin_dynamic_bit_length_zero_allocation_trace();
    let local_sub_fused = fused_zero_prefix_resources(&build_fused_zero_prefix_update(
        LOCAL_WIDTH,
        LOCAL_ACCUMULATOR_WIDTH,
        true,
        false,
    ));
    let local_sub_trace = finish_dynamic_bit_length_zero_allocation_trace();
    assert_zero_dynamic_bit_length_trace(local_sub_trace);
    assert!(local_add_fused.extra_qubits <= local_add_baseline.extra_qubits);
    assert!(local_add_fused.emitted_ops <= local_add_baseline.emitted_ops);
    assert!(local_add_fused.emitted_toffoli <= local_add_baseline.emitted_toffoli);
    assert!(local_sub_fused.extra_qubits <= local_sub_baseline.extra_qubits);
    assert!(local_sub_fused.emitted_ops <= local_sub_baseline.emitted_ops);
    assert!(local_sub_fused.emitted_toffoli <= local_sub_baseline.emitted_toffoli);

    configure_fused_zero_prefix_proof(true, false);
    let scheduled_baseline = profile_reference_scheduled_inversion();
    configure_fused_zero_prefix_proof(true, true);
    let scheduled_fused = profile_reference_scheduled_inversion();
    assert_eq!(scheduled_fused.steps, scheduled_baseline.steps);
    assert!(scheduled_fused.inversion_peak_qubits <= scheduled_baseline.inversion_peak_qubits);
    assert!(scheduled_fused.emitted_ops <= scheduled_baseline.emitted_ops);
    assert!(scheduled_fused.emitted_toffoli <= scheduled_baseline.emitted_toffoli);

    let precondition_rejections = fused_zero_prefix_precondition_rejections();
    configure_fused_zero_prefix_proof(true, true);

    FusedZeroPrefixBitLengthProofReport {
        widths_checked: MAX_SOURCE_WIDTH + 1,
        accumulator_width: ACCUMULATOR_WIDTH,
        accumulator_values_per_source: 1usize << ACCUMULATOR_WIDTH,
        update_cases_checked,
        controlled_cases_checked,
        directions_checked: 2,
        control_values_checked: 2,
        baseline_equivalence_checks,
        default_stream_equivalence_checks,
        phase_clean_checks,
        ancilla_clean_checks,
        precondition_rejections,
        trace_sensitivity_flag_allocations: trace_sensitivity.flag_allocations,
        trace_sensitivity_carry_allocations: trace_sensitivity.carry_allocations,
        trace_sensitivity_prefix_allocations: trace_sensitivity.prefix_allocations,
        fused_zero_flag_allocations: fused_trace.flag_allocations,
        fused_zero_carry_allocations: fused_trace.carry_allocations,
        fused_zero_prefix_allocations: fused_trace.prefix_allocations,
        maximum_baseline_extra_qubits: maximum_baseline.extra_qubits,
        maximum_fused_extra_qubits: maximum_fused.extra_qubits,
        maximum_baseline_emitted_ops: maximum_baseline.emitted_ops,
        maximum_fused_emitted_ops: maximum_fused.emitted_ops,
        maximum_baseline_emitted_toffoli: maximum_baseline.emitted_toffoli,
        maximum_fused_emitted_toffoli: maximum_fused.emitted_toffoli,
        width_direction_toffoli_increase_cases,
        maximum_width_direction_toffoli_increase,
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
fn rotated_bitlen_scratch_reuse_requested() -> bool {
    std::env::var("LOWQ_REUSE_ROTATED_BITLEN_SCRATCH")
        .ok()
        .as_deref()
        == Some("1")
}

fn coefficient_raw_bitlen_loan_requested() -> bool {
    std::env::var(COEFFICIENT_RAW_BITLEN_LOAN_FLAG)
        .ok()
        .as_deref()
        == Some("1")
}

fn promised_l_q_swap_borrow_requested() -> bool {
    std::env::var(PROMISED_LQ_SWAP_BORROW_FLAG).ok().as_deref() == Some("1")
}

fn promised_swap_support_lifetime_fusion_requested() -> bool {
    std::env::var(PROMISED_SWAP_SUPPORT_LIFETIME_FUSION_FLAG)
        .ok()
        .as_deref()
        == Some("1")
}

fn split_coefficient_rotation_lifetime_requested() -> bool {
    std::env::var(SPLIT_COEFFICIENT_ROTATION_LIFETIME_FLAG)
        .ok()
        .as_deref()
        == Some("1")
}

fn coefficient_less_than_lane_reuse_requested() -> bool {
    std::env::var(COEFFICIENT_LESS_THAN_LANE_REUSE_FLAG)
        .ok()
        .as_deref()
        == Some("1")
}

fn clean_chain_coefficient_add_lender_requested() -> bool {
    std::env::var(CLEAN_CHAIN_COEFFICIENT_ADD_LENDER_FLAG)
        .ok()
        .as_deref()
        == Some("1")
}

fn paired_bitlen_source_complement_requested() -> bool {
    std::env::var(PAIRED_BITLEN_SOURCE_COMPLEMENT_FLAG)
        .ok()
        .as_deref()
        == Some("1")
}

fn coefficient_nonnegative_x_cancel_requested() -> bool {
    std::env::var(COEFFICIENT_NONNEGATIVE_X_CANCEL_FLAG)
        .ok()
        .as_deref()
        == Some("1")
}

fn q845_lifetime_coefficient_fusion_requested() -> bool {
    std::env::var(Q845_LIFETIME_COEFFICIENT_FUSION_FLAG)
        .ok()
        .as_deref()
        == Some("1")
}

fn q845_swap_only_t_prime_length_requested() -> bool {
    std::env::var(Q845_SWAP_ONLY_T_PRIME_LENGTH_FLAG)
        .ok()
        .as_deref()
        == Some("1")
}

fn q851_truncated_swap_only_guard_requested() -> bool {
    std::env::var(Q851_TRUNCATED_SWAP_ONLY_GUARD_FLAG)
        .ok()
        .as_deref()
        == Some("1")
}

fn q851_fixed_sign_event_requested() -> bool {
    std::env::var(Q851_FIXED_SIGN_EVENT_FLAG).ok().as_deref() == Some("1")
}

fn q830_dirty_fixed_sign_event_requested() -> bool {
    std::env::var(Q830_DIRTY_FIXED_SIGN_EVENT_FLAG)
        .ok()
        .as_deref()
        == Some("1")
}

fn q830_direct_swap_metadata_requested() -> bool {
    std::env::var(Q830_DIRECT_SWAP_METADATA_FLAG)
        .ok()
        .as_deref()
        == Some("1")
}

fn q830_coefficient_counter_relocation_requested() -> bool {
    std::env::var(Q830_COEFFICIENT_COUNTER_RELOCATION_FLAG)
        .ok()
        .as_deref()
        == Some("1")
}

fn sub800_inplace_guard_address_requested() -> bool {
    std::env::var(SUB800_INPLACE_GUARD_ADDRESS_FLAG)
        .ok()
        .as_deref()
        == Some("1")
}

fn sub800_raw_prefix_preserved_lender_requested() -> bool {
    std::env::var(SUB800_RAW_PREFIX_PRESERVED_LENDER_FLAG)
        .ok()
        .as_deref()
        == Some("1")
}

fn sub800_raw_prefix_predicate_lender_requested() -> bool {
    std::env::var(SUB800_RAW_PREFIX_PREDICATE_LENDER_FLAG)
        .ok()
        .as_deref()
        == Some("1")
}

fn sub800_mixed_boundary_scratch_extension_requested() -> bool {
    std::env::var(SUB800_MIXED_BOUNDARY_SCRATCH_EXTENSION_FLAG)
        .ok()
        .as_deref()
        == Some("1")
}

fn sub800_borrowed_rotated_underflow_requested() -> bool {
    std::env::var(SUB800_BORROWED_ROTATED_UNDERFLOW_FLAG)
        .ok()
        .as_deref()
        == Some("1")
}
fn sub800_split_mixed_rotated_length_requested() -> bool {
    std::env::var(SUB800_SPLIT_MIXED_ROTATED_LENGTH_FLAG)
        .ok()
        .as_deref()
        == Some("1")
}
fn sub800_split_same_rotated_length_requested() -> bool {
    std::env::var(SUB800_SPLIT_SAME_ROTATED_LENGTH_FLAG)
        .ok()
        .as_deref()
        == Some("1")
}
fn sub800_split_two_high_rotated_length_requested() -> bool {
    std::env::var(SUB800_SPLIT_TWO_HIGH_ROTATED_LENGTH_FLAG)
        .ok()
        .as_deref()
        == Some("1")
}
fn sub800_split_three_high_rotated_length_requested() -> bool {
    std::env::var(SUB800_SPLIT_THREE_HIGH_ROTATED_LENGTH_FLAG)
        .ok()
        .as_deref()
        == Some("1")
}
fn sub800_split_four_high_rotated_length_requested() -> bool {
    std::env::var(SUB800_SPLIT_FOUR_HIGH_ROTATED_LENGTH_FLAG)
        .ok()
        .as_deref()
        == Some("1")
}
fn sub800_uls_clean_lender_requested() -> bool {
    q845_swap_only_t_prime_length_requested()
        && std::env::var(SUB800_ULS_CLEAN_LENDER_FLAG)
        .ok()
        .as_deref()
        == Some("1")
}
fn sub800_uls_fused_target_requested() -> bool {
    q845_swap_only_t_prime_length_requested()
        && std::env::var(SUB800_ULS_FUSED_TARGET_FLAG)
        .ok()
        .as_deref()
        == Some("1")
}
fn sub800_uls_direct_selector_requested() -> bool {
    q845_swap_only_t_prime_length_requested()
        && sub800_uls_fused_target_requested()
        && std::env::var(SUB800_ULS_DIRECT_SELECTOR_FLAG)
            .ok()
        .as_deref()
        == Some("1")
}
fn q839_seven_plateau_lenders_requested() -> bool {
    if std::env::var(Q839_SEVEN_PLATEAU_LENDERS_FLAG)
        .ok()
        .as_deref()
        != Some("1")
    {
        return false;
    }
    assert!(
        q845_swap_only_t_prime_length_requested(),
        "seven-plateau lending requires the Q845 swap-only t-prime route"
    );
    assert!(
        q845_lifetime_coefficient_fusion_requested() && sub800_uls_fused_target_requested(),
        "seven-plateau ULS lending requires the fused Q845 coefficient target"
    );
    assert!(
        q851_fixed_sign_event_requested(),
        "seven-plateau ULS lending requires the fixed-sign cursor event"
    );
    assert!(
        !sub800_uls_direct_selector_requested(),
        "seven-plateau ULS lending and the direct selector are mutually exclusive"
    );
    true
}
fn q845_swap_only_coefficient_dependencies_satisfied() -> bool {
    (!q845_swap_only_t_prime_length_requested() || q845_lifetime_coefficient_fusion_requested())
        && (!q851_truncated_swap_only_guard_requested()
            || q845_swap_only_t_prime_length_requested())
        && (!q851_fixed_sign_event_requested() || q845_swap_only_t_prime_length_requested())
        && (!sub800_inplace_guard_address_requested()
            || q845_swap_only_t_prime_length_requested())
}

fn q845_swap_only_swap_dependencies_satisfied() -> bool {
    !q845_swap_only_t_prime_length_requested() || promised_l_q_swap_borrow_requested()
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
struct RawBitLengthAllocationTrace {
    raw_rotation_carry_allocations: usize,
    raw_rotation_overflow_allocations: usize,
    raw_rotation_split_releases: usize,
    raw_rotation_split_recomputes: usize,
    raw_rotation_lanes_released: usize,
    rotated_boundary_qubits_allocated: usize,
    rotated_inplace_boundary_uses: usize,
    rotated_carry_allocations: usize,
    rotated_overflow_allocations: usize,
    rotated_enabled_allocations: usize,
}

thread_local! {
    static RAW_BIT_LENGTH_ALLOCATION_TRACE: std::cell::Cell<(
        bool,
        RawBitLengthAllocationTrace,
    )> = std::cell::Cell::new((false, RawBitLengthAllocationTrace {
        raw_rotation_carry_allocations: 0,
        raw_rotation_overflow_allocations: 0,
        raw_rotation_split_releases: 0,
        raw_rotation_split_recomputes: 0,
        raw_rotation_lanes_released: 0,
        rotated_boundary_qubits_allocated: 0,
        rotated_inplace_boundary_uses: 0,
        rotated_carry_allocations: 0,
        rotated_overflow_allocations: 0,
        rotated_enabled_allocations: 0,
    }));
}

fn begin_raw_bit_length_allocation_trace() {
    RAW_BIT_LENGTH_ALLOCATION_TRACE.with(|trace| {
        trace.set((true, RawBitLengthAllocationTrace::default()));
    });
}

fn finish_raw_bit_length_allocation_trace() -> RawBitLengthAllocationTrace {
    RAW_BIT_LENGTH_ALLOCATION_TRACE.with(|trace| {
        let (_, snapshot) = trace.get();
        trace.set((false, snapshot));
        snapshot
    })
}

fn trace_raw_bit_length_allocations(delta: RawBitLengthAllocationTrace) {
    RAW_BIT_LENGTH_ALLOCATION_TRACE.with(|trace| {
        let (enabled, mut snapshot) = trace.get();
        if enabled {
            snapshot.raw_rotation_carry_allocations += delta.raw_rotation_carry_allocations;
            snapshot.raw_rotation_overflow_allocations += delta.raw_rotation_overflow_allocations;
            snapshot.raw_rotation_split_releases += delta.raw_rotation_split_releases;
            snapshot.raw_rotation_split_recomputes += delta.raw_rotation_split_recomputes;
            snapshot.raw_rotation_lanes_released += delta.raw_rotation_lanes_released;
            snapshot.rotated_boundary_qubits_allocated += delta.rotated_boundary_qubits_allocated;
            snapshot.rotated_inplace_boundary_uses += delta.rotated_inplace_boundary_uses;
            snapshot.rotated_carry_allocations += delta.rotated_carry_allocations;
            snapshot.rotated_overflow_allocations += delta.rotated_overflow_allocations;
            snapshot.rotated_enabled_allocations += delta.rotated_enabled_allocations;
            trace.set((enabled, snapshot));
        }
    });
}

fn inplace_rotated_bitlen_boundary_requested() -> bool {
    std::env::var(INPLACE_ROTATED_BITLEN_BOUNDARY_FLAG)
        .ok()
        .as_deref()
        == Some("1")
}

fn assert_rotated_bitlen_scratch_disjoint(
    control: &QReg,
    boundary: &[QReg],
    source: &[&QReg],
    output: &[QReg],
    scratch: &[&QReg],
) -> bool {
    let borrowed = scratch.len() >= 3;
    if borrowed {
        for (index, lane) in scratch[..3].iter().enumerate() {
            assert!(
                scratch[..index].iter().all(|other| other.id() != lane.id()),
                "rotated bit-length scratch lane {index} aliases an earlier lane"
            );
            assert_ne!(lane.id(), control.id());
            assert!(boundary.iter().all(|other| other.id() != lane.id()));
            assert!(source.iter().all(|other| other.id() != lane.id()));
            assert!(output.iter().all(|other| other.id() != lane.id()));
        }
    }
    borrowed
}

fn toggle_subtraction_borrow_anf(
    circ: &mut Circuit,
    x: &QReg,
    y: &QReg,
    borrow: &QReg,
    target: &QReg,
) {
    // br(x,y,b) = y XOR b XOR yb XOR xy XOR xb.
    circ.cx(y, target);
    circ.cx(borrow, target);
    circ.ccx(y, borrow, target);
    circ.ccx(x, y, target);
    circ.ccx(x, borrow, target);
}

fn toggle_controlled_subtraction_borrow_anf(
    circ: &mut Circuit,
    control: &QReg,
    x: &QReg,
    y: Option<&QReg>,
    borrow: &QReg,
    target: &QReg,
    dirty: &QReg,
) {
    use crate::point_add::trailmix_port::arith::mcx::mcx_dirty_any_k;

    if let Some(y) = y {
        circ.ccx(control, y, target);
        circ.ccx(control, borrow, target);
        mcx_dirty_any_k(circ, &[control, y, borrow], target, dirty);
        mcx_dirty_any_k(circ, &[control, x, y], target, dirty);
        mcx_dirty_any_k(circ, &[control, x, borrow], target, dirty);
    } else {
        // br(x,0,b) = b XOR xb.
        circ.ccx(control, borrow, target);
        mcx_dirty_any_k(circ, &[control, x, borrow], target, dirty);
    }
}

fn controlled_xor_saturating_difference_materialized_boundary(
    circ: &mut Circuit,
    control: &QReg,
    boundary: &[QReg],
    source: &[&QReg],
    length: &[QReg],
    output: &[QReg],
    scratch: &[&QReg],
) {
    trace_raw_bit_length_allocations(RawBitLengthAllocationTrace {
        rotated_boundary_qubits_allocated: length.len(),
        ..RawBitLengthAllocationTrace::default()
    });
    let boundary_copy = circ.alloc_qreg_bits("rs.rotated-bitlen.boundary", length.len());
    for (source_bit, target_bit) in boundary.iter().zip(&boundary_copy) {
        circ.cx(source_bit, target_bit);
    }
    let borrowed =
        assert_rotated_bitlen_scratch_disjoint(control, boundary, source, output, scratch);
    if !borrowed {
        trace_raw_bit_length_allocations(RawBitLengthAllocationTrace {
            rotated_carry_allocations: 1,
            rotated_overflow_allocations: 1,
            ..RawBitLengthAllocationTrace::default()
        });
    }
    let carry_owned = (!borrowed).then(|| circ.alloc_qreg("rs.rotated-bitlen.carry"));
    let overflow_owned = (!borrowed).then(|| circ.alloc_qreg("rs.rotated-bitlen.overflow"));
    let carry = if borrowed {
        scratch[0]
    } else {
        carry_owned.as_ref().expect("owned rotated carry")
    };
    let overflow = if borrowed {
        scratch[1]
    } else {
        overflow_owned.as_ref().expect("owned rotated overflow")
    };
    cuccaro_sub_mod_2n(circ, &boundary_copy, length, carry, overflow);

    let sign = &length[length.len() - 1];
    if !borrowed {
        trace_raw_bit_length_allocations(RawBitLengthAllocationTrace {
            rotated_enabled_allocations: 1,
            ..RawBitLengthAllocationTrace::default()
        });
    }
    let enabled_owned = (!borrowed).then(|| circ.alloc_qreg("rs.rotated-bitlen.enabled"));
    let enabled = if borrowed {
        scratch[2]
    } else {
        enabled_owned.as_ref().expect("owned rotated enable")
    };
    circ.x(sign);
    circ.ccx(control, sign, enabled);
    for (difference_bit, output_bit) in length.iter().zip(output) {
        circ.ccx(enabled, difference_bit, output_bit);
    }
    circ.ccx(control, sign, enabled);
    circ.x(sign);
    if let Some(enabled) = enabled_owned {
        circ.zero_and_free(enabled);
    }

    cuccaro_add_mod_2n(circ, &boundary_copy, length, carry, overflow);
    if let Some(overflow) = overflow_owned {
        circ.zero_and_free(overflow);
    }
    if let Some(carry) = carry_owned {
        circ.zero_and_free(carry);
    }
    for (source_bit, target_bit) in boundary.iter().zip(&boundary_copy) {
        circ.cx(source_bit, target_bit);
    }
    for lane in boundary_copy {
        circ.zero_and_free(lane);
    }
}

fn controlled_xor_saturating_difference_inplace_boundary(
    circ: &mut Circuit,
    control: &QReg,
    boundary: &[QReg],
    source: &[&QReg],
    length: &[QReg],
    output: &[QReg],
    scratch: &[&QReg],
) {
    let split_four_high_length = length.len() + 4 == boundary.len();
    let split_three_high_length = length.len() + 3 == boundary.len();
    let split_two_high_length = length.len() + 2 == boundary.len();
    let split_same_length = length.len() + 1 == boundary.len();
    let borrowed_underflow = split_four_high_length
        || split_three_high_length
        || split_two_high_length
        || split_same_length
        || length.len() == boundary.len();
    assert!(
        split_four_high_length
            || split_three_high_length
            || split_two_high_length
            || split_same_length
            || length.len() == boundary.len()
            || length.len() == boundary.len() + 1
    );
    trace_raw_bit_length_allocations(RawBitLengthAllocationTrace {
        rotated_inplace_boundary_uses: 1,
        ..RawBitLengthAllocationTrace::default()
    });
    let borrowed =
        assert_rotated_bitlen_scratch_disjoint(control, boundary, source, output, scratch);
    assert!(
        !borrowed_underflow || borrowed,
        "borrowed rotated underflow requires three disjoint clean lenders"
    );
    if !borrowed {
        trace_raw_bit_length_allocations(RawBitLengthAllocationTrace {
            rotated_carry_allocations: 1,
            rotated_enabled_allocations: 1,
            ..RawBitLengthAllocationTrace::default()
        });
    }
    let carry_owned = (!borrowed).then(|| circ.alloc_qreg("rs.rotated-bitlen.carry"));
    let enabled_owned = (!borrowed).then(|| circ.alloc_qreg("rs.rotated-bitlen.enabled"));
    let carry = if borrowed {
        scratch[0]
    } else {
        carry_owned.as_ref().expect("owned rotated carry")
    };
    let enabled = if borrowed_underflow {
        carry
    } else if borrowed {
        scratch[2]
    } else {
        enabled_owned.as_ref().expect("owned rotated enable")
    };

    if split_four_high_length {
        assert!(borrowed);
        assert!(scratch.len() >= 9);
        assert_eq!(length.len() + 4, boundary.len());
        assert_eq!(output.len(), boundary.len());
        assert!(!length.is_empty());
        assert!(!source.is_empty());
        let boundary_low = &boundary[..length.len()];
        let boundary_high5 = &boundary[length.len()];
        let boundary_high6 = &boundary[length.len() + 1];
        let boundary_high7 = &boundary[length.len() + 2];
        let boundary_high8 = &boundary[length.len() + 3];
        let borrow = scratch[1];
        let high5_borrow = scratch[2];
        let high6_borrow = scratch[3];
        let high7_borrow = scratch[4];
        let length_high5 = scratch[5];
        let length_high6 = scratch[6];
        let length_high7 = scratch[7];
        let length_high8 = scratch[8];
        let dirty = source[0];
        for (index, lane) in scratch[..9].iter().enumerate() {
            assert!(scratch[..index].iter().all(|other| other.id() != lane.id()));
            assert_ne!(lane.id(), control.id());
            assert!(boundary.iter().all(|other| other.id() != lane.id()));
            assert!(source.iter().all(|other| other.id() != lane.id()));
            assert!(output.iter().all(|other| other.id() != lane.id()));
            assert!(length.iter().all(|other| other.id() != lane.id()));
        }

        cuccaro_sub_mod_2n(circ, boundary_low, length, carry, borrow);
        toggle_subtraction_borrow_anf(circ, length_high5, boundary_high5, borrow, high5_borrow);
        toggle_subtraction_borrow_anf(
            circ,
            length_high6,
            boundary_high6,
            high5_borrow,
            high6_borrow,
        );
        toggle_subtraction_borrow_anf(
            circ,
            length_high7,
            boundary_high7,
            high6_borrow,
            high7_borrow,
        );

        circ.cx(control, enabled);
        toggle_controlled_subtraction_borrow_anf(
            circ,
            control,
            length_high8,
            Some(boundary_high8),
            high7_borrow,
            enabled,
            dirty,
        );

        for (high, boundary_high, incoming) in [
            (length_high5, boundary_high5, borrow),
            (length_high6, boundary_high6, high5_borrow),
            (length_high7, boundary_high7, high6_borrow),
            (length_high8, boundary_high8, high7_borrow),
        ] {
            circ.cx(boundary_high, high);
            circ.cx(incoming, high);
        }
        for (difference_bit, output_bit) in length.iter().zip(&output[..length.len()]) {
            circ.ccx(enabled, difference_bit, output_bit);
        }
        for (offset, high) in [length_high5, length_high6, length_high7, length_high8]
            .into_iter()
            .enumerate()
        {
            circ.ccx(enabled, high, &output[length.len() + offset]);
        }
        for (high, boundary_high, incoming) in [
            (length_high8, boundary_high8, high7_borrow),
            (length_high7, boundary_high7, high6_borrow),
            (length_high6, boundary_high6, high5_borrow),
            (length_high5, boundary_high5, borrow),
        ] {
            circ.cx(incoming, high);
            circ.cx(boundary_high, high);
        }

        toggle_controlled_subtraction_borrow_anf(
            circ,
            control,
            length_high8,
            Some(boundary_high8),
            high7_borrow,
            enabled,
            dirty,
        );
        circ.cx(control, enabled);
        toggle_subtraction_borrow_anf(
            circ,
            length_high7,
            boundary_high7,
            high6_borrow,
            high7_borrow,
        );
        toggle_subtraction_borrow_anf(
            circ,
            length_high6,
            boundary_high6,
            high5_borrow,
            high6_borrow,
        );
        toggle_subtraction_borrow_anf(circ, length_high5, boundary_high5, borrow, high5_borrow);
        cuccaro_add_mod_2n(circ, boundary_low, length, carry, borrow);
        return;
    }

    if split_three_high_length {
        assert!(borrowed);
        assert!(scratch.len() >= 8);
        assert_eq!(length.len() + 3, boundary.len());
        assert_eq!(output.len(), boundary.len());
        assert!(!length.is_empty());
        assert!(!source.is_empty());
        let boundary_low = &boundary[..length.len()];
        let boundary_high6 = &boundary[length.len()];
        let boundary_high7 = &boundary[length.len() + 1];
        let boundary_high8 = &boundary[length.len() + 2];
        let borrow = scratch[1];
        let high6_borrow = scratch[2];
        let high7_borrow = scratch[3];
        let length_high6 = scratch[5];
        let length_high7 = scratch[6];
        let length_high8 = scratch[7];
        let dirty = source[0];
        for (index, lane) in scratch[..8].iter().enumerate() {
            assert!(scratch[..index].iter().all(|other| other.id() != lane.id()));
            assert_ne!(lane.id(), control.id());
            assert!(boundary.iter().all(|other| other.id() != lane.id()));
            assert!(source.iter().all(|other| other.id() != lane.id()));
            assert!(output.iter().all(|other| other.id() != lane.id()));
            assert!(length.iter().all(|other| other.id() != lane.id()));
        }
        assert_ne!(dirty.id(), control.id());
        assert!(boundary.iter().all(|lane| lane.id() != dirty.id()));
        assert!(output.iter().all(|lane| lane.id() != dirty.id()));
        assert!(length.iter().all(|lane| lane.id() != dirty.id()));

        cuccaro_sub_mod_2n(circ, boundary_low, length, carry, borrow);
        toggle_subtraction_borrow_anf(
            circ,
            length_high6,
            boundary_high6,
            borrow,
            high6_borrow,
        );
        toggle_subtraction_borrow_anf(
            circ,
            length_high7,
            boundary_high7,
            high6_borrow,
            high7_borrow,
        );

        // Start from c, then toggle c*br(h8,k8,v) to obtain c AND !borrow.
        circ.cx(control, enabled);
        toggle_controlled_subtraction_borrow_anf(
            circ,
            control,
            length_high8,
            Some(boundary_high8),
            high7_borrow,
            enabled,
            dirty,
        );

        circ.cx(boundary_high6, length_high6);
        circ.cx(borrow, length_high6);
        circ.cx(boundary_high7, length_high7);
        circ.cx(high6_borrow, length_high7);
        circ.cx(boundary_high8, length_high8);
        circ.cx(high7_borrow, length_high8);
        for (difference_bit, output_bit) in length.iter().zip(&output[..length.len()]) {
            circ.ccx(enabled, difference_bit, output_bit);
        }
        circ.ccx(enabled, length_high6, &output[length.len()]);
        circ.ccx(enabled, length_high7, &output[length.len() + 1]);
        circ.ccx(enabled, length_high8, &output[length.len() + 2]);
        circ.cx(high7_borrow, length_high8);
        circ.cx(boundary_high8, length_high8);
        circ.cx(high6_borrow, length_high7);
        circ.cx(boundary_high7, length_high7);
        circ.cx(borrow, length_high6);
        circ.cx(boundary_high6, length_high6);

        toggle_controlled_subtraction_borrow_anf(
            circ,
            control,
            length_high8,
            Some(boundary_high8),
            high7_borrow,
            enabled,
            dirty,
        );
        circ.cx(control, enabled);
        toggle_subtraction_borrow_anf(
            circ,
            length_high7,
            boundary_high7,
            high6_borrow,
            high7_borrow,
        );
        toggle_subtraction_borrow_anf(
            circ,
            length_high6,
            boundary_high6,
            borrow,
            high6_borrow,
        );
        cuccaro_add_mod_2n(circ, boundary_low, length, carry, borrow);
        return;
    }

    if split_two_high_length {
        use crate::point_add::trailmix_port::arith::mcx::mcx_dirty_any_k;

        assert!(borrowed);
        assert!(scratch.len() >= 7);
        assert_eq!(length.len() + 2, boundary.len());
        assert_eq!(output.len(), boundary.len());
        assert!(!length.is_empty());
        let boundary_low = &boundary[..length.len()];
        let boundary_high7 = &boundary[length.len()];
        let boundary_high8 = &boundary[length.len() + 1];
        let borrow = scratch[1];
        let high_borrow = scratch[2];
        let length_high7 = scratch[5];
        let length_high8 = scratch[6];
        let dirty = &length[0];
        for (index, lane) in scratch[..7].iter().enumerate() {
            assert!(scratch[..index].iter().all(|other| other.id() != lane.id()));
            assert_ne!(lane.id(), control.id());
            assert!(boundary.iter().all(|other| other.id() != lane.id()));
            assert!(source.iter().all(|other| other.id() != lane.id()));
            assert!(output.iter().all(|other| other.id() != lane.id()));
            assert!(length.iter().all(|other| other.id() != lane.id()));
        }

        cuccaro_sub_mod_2n(circ, boundary_low, length, carry, borrow);

        // v = br(a,g,u) = g XOR u XOR gu XOR ag XOR au.
        circ.cx(boundary_high7, high_borrow);
        circ.cx(borrow, high_borrow);
        circ.ccx(boundary_high7, borrow, high_borrow);
        circ.ccx(length_high7, boundary_high7, high_borrow);
        circ.ccx(length_high7, borrow, high_borrow);

        // U = br(h,k,v). Starting from c, toggle c*U to obtain c AND !U.
        circ.cx(control, enabled);
        circ.ccx(control, boundary_high8, enabled);
        circ.ccx(control, high_borrow, enabled);
        mcx_dirty_any_k(
            circ,
            &[control, boundary_high8, high_borrow],
            enabled,
            dirty,
        );
        mcx_dirty_any_k(
            circ,
            &[control, length_high8, boundary_high8],
            enabled,
            dirty,
        );
        mcx_dirty_any_k(
            circ,
            &[control, length_high8, high_borrow],
            enabled,
            dirty,
        );

        // d7 = a XOR g XOR u and d8 = h XOR k XOR v.
        circ.cx(boundary_high7, length_high7);
        circ.cx(borrow, length_high7);
        circ.cx(boundary_high8, length_high8);
        circ.cx(high_borrow, length_high8);
        for (difference_bit, output_bit) in length.iter().zip(&output[..length.len()]) {
            circ.ccx(enabled, difference_bit, output_bit);
        }
        circ.ccx(enabled, length_high7, &output[length.len()]);
        circ.ccx(enabled, length_high8, &output[length.len() + 1]);
        circ.cx(high_borrow, length_high8);
        circ.cx(boundary_high8, length_high8);
        circ.cx(borrow, length_high7);
        circ.cx(boundary_high7, length_high7);

        mcx_dirty_any_k(
            circ,
            &[control, length_high8, high_borrow],
            enabled,
            dirty,
        );
        mcx_dirty_any_k(
            circ,
            &[control, length_high8, boundary_high8],
            enabled,
            dirty,
        );
        mcx_dirty_any_k(
            circ,
            &[control, boundary_high8, high_borrow],
            enabled,
            dirty,
        );
        circ.ccx(control, high_borrow, enabled);
        circ.ccx(control, boundary_high8, enabled);
        circ.cx(control, enabled);

        circ.ccx(length_high7, borrow, high_borrow);
        circ.ccx(length_high7, boundary_high7, high_borrow);
        circ.ccx(boundary_high7, borrow, high_borrow);
        circ.cx(borrow, high_borrow);
        circ.cx(boundary_high7, high_borrow);
        cuccaro_add_mod_2n(circ, boundary_low, length, carry, borrow);
        return;
    }

    if split_same_length {
        use crate::point_add::trailmix_port::arith::mcx::mcx_dirty_any_k;

        assert!(borrowed);
        assert_eq!(length.len() + 1, boundary.len());
        assert_eq!(output.len(), boundary.len());
        assert!(!length.is_empty());
        let boundary_high = &boundary[boundary.len() - 1];
        let boundary_low = &boundary[..length.len()];
        let borrow = scratch[1];
        let length_high = scratch[2];
        let dirty = &length[0];

        cuccaro_sub_mod_2n(circ, boundary_low, length, carry, borrow);

        // For L=(h,l), B=(g,b), and u=[l<b], the final borrow is
        // U = g XOR u XOR gu XOR hg XOR hu over F_2. Starting from the
        // control and toggling control*U materializes control AND NOT U.
        circ.cx(control, enabled);
        circ.ccx(control, boundary_high, enabled);
        circ.ccx(control, borrow, enabled);
        mcx_dirty_any_k(circ, &[control, boundary_high, borrow], enabled, dirty);
        mcx_dirty_any_k(
            circ,
            &[control, length_high, boundary_high],
            enabled,
            dirty,
        );
        mcx_dirty_any_k(circ, &[control, length_high, borrow], enabled, dirty);

        circ.cx(boundary_high, length_high);
        circ.cx(borrow, length_high);
        for (difference_bit, output_bit) in
            length.iter().zip(&output[..length.len()])
        {
            circ.ccx(enabled, difference_bit, output_bit);
        }
        circ.ccx(enabled, length_high, &output[output.len() - 1]);
        circ.cx(borrow, length_high);
        circ.cx(boundary_high, length_high);

        mcx_dirty_any_k(circ, &[control, length_high, borrow], enabled, dirty);
        mcx_dirty_any_k(
            circ,
            &[control, length_high, boundary_high],
            enabled,
            dirty,
        );
        mcx_dirty_any_k(circ, &[control, boundary_high, borrow], enabled, dirty);
        circ.ccx(control, borrow, enabled);
        circ.ccx(control, boundary_high, enabled);
        circ.cx(control, enabled);

        cuccaro_add_mod_2n(circ, boundary_low, length, carry, borrow);
        return;
    }

    let (low_length, sign) = if borrowed_underflow {
        (length, scratch[1])
    } else {
        let (low_length, sign_lane) = length.split_at(boundary.len());
        (low_length, &sign_lane[0])
    };

    // The inverse-Cuccaro carry is restored after subtraction. In the compact
    // route, keep the underflow in a borrowed lane and reuse the restored carry
    // as the enable predicate before running the exact inverse addition.
    cuccaro_sub_mod_2n(circ, boundary, low_length, carry, sign);
    circ.x(sign);
    circ.ccx(control, sign, enabled);
    for (difference_bit, output_bit) in low_length.iter().zip(output) {
        circ.ccx(enabled, difference_bit, output_bit);
    }
    circ.ccx(control, sign, enabled);
    circ.x(sign);
    cuccaro_add_mod_2n(circ, boundary, low_length, carry, sign);

    if let Some(enabled) = enabled_owned {
        circ.zero_and_free(enabled);
    }
    if let Some(carry) = carry_owned {
        circ.zero_and_free(carry);
    }
}

fn controlled_xor_saturating_difference_inplace_mixed_boundary(
    circ: &mut Circuit,
    control: &QReg,
    boundary: &[QReg],
    source: &[&QReg],
    length: &[QReg],
    output: &[QReg],
    scratch: &[&QReg],
) {
    assert_eq!(output.len(), boundary.len() + 1);
    let split_four_high_length = length.len() + 4 == output.len();
    let split_three_high_length = length.len() + 3 == output.len();
    let split_two_high_length = length.len() + 2 == output.len();
    let split_mixed_length = length.len() + 1 == output.len();
    let borrowed_underflow = split_four_high_length
        || split_three_high_length
        || split_two_high_length
        || split_mixed_length
        || length.len() == output.len();
    assert!(
        split_four_high_length
            || split_three_high_length
            || split_two_high_length
            || split_mixed_length
            || length.len() == output.len()
            || length.len() == output.len() + 1
    );
    trace_raw_bit_length_allocations(RawBitLengthAllocationTrace {
        rotated_inplace_boundary_uses: 1,
        ..RawBitLengthAllocationTrace::default()
    });
    let borrowed =
        assert_rotated_bitlen_scratch_disjoint(control, boundary, source, output, scratch);
    assert!(
        !borrowed_underflow || borrowed,
        "borrowed mixed-width underflow requires three disjoint clean lenders"
    );
    if !borrowed {
        trace_raw_bit_length_allocations(RawBitLengthAllocationTrace {
            rotated_carry_allocations: 1,
            rotated_enabled_allocations: 1,
            ..RawBitLengthAllocationTrace::default()
        });
    }
    let carry_owned = (!borrowed).then(|| circ.alloc_qreg("rs.rotated-bitlen.carry"));
    let enabled_owned = (!borrowed).then(|| circ.alloc_qreg("rs.rotated-bitlen.enabled"));
    let carry = if borrowed {
        scratch[0]
    } else {
        carry_owned.as_ref().expect("owned rotated carry")
    };
    let enabled = if borrowed_underflow {
        carry
    } else if borrowed {
        scratch[2]
    } else {
        enabled_owned.as_ref().expect("owned rotated enable")
    };

    if split_four_high_length {
        assert!(borrowed);
        assert!(scratch.len() >= 9);
        assert_eq!(length.len() + 3, boundary.len());
        assert!(!length.is_empty());
        assert!(!source.is_empty());
        let boundary_low = &boundary[..length.len()];
        let boundary_high5 = &boundary[length.len()];
        let boundary_high6 = &boundary[length.len() + 1];
        let boundary_high7 = &boundary[length.len() + 2];
        let borrow = scratch[1];
        let high5_borrow = scratch[2];
        let high6_borrow = scratch[3];
        let high7_borrow = scratch[4];
        let length_high5 = scratch[5];
        let length_high6 = scratch[6];
        let length_high7 = scratch[7];
        let length_high8 = scratch[8];
        let dirty = source[0];
        for (index, lane) in scratch[..9].iter().enumerate() {
            assert!(scratch[..index].iter().all(|other| other.id() != lane.id()));
            assert_ne!(lane.id(), control.id());
            assert!(boundary.iter().all(|other| other.id() != lane.id()));
            assert!(source.iter().all(|other| other.id() != lane.id()));
            assert!(output.iter().all(|other| other.id() != lane.id()));
            assert!(length.iter().all(|other| other.id() != lane.id()));
        }

        cuccaro_sub_mod_2n(circ, boundary_low, length, carry, borrow);
        toggle_subtraction_borrow_anf(circ, length_high5, boundary_high5, borrow, high5_borrow);
        toggle_subtraction_borrow_anf(
            circ,
            length_high6,
            boundary_high6,
            high5_borrow,
            high6_borrow,
        );
        toggle_subtraction_borrow_anf(
            circ,
            length_high7,
            boundary_high7,
            high6_borrow,
            high7_borrow,
        );

        circ.cx(control, enabled);
        toggle_controlled_subtraction_borrow_anf(
            circ,
            control,
            length_high8,
            None,
            high7_borrow,
            enabled,
            dirty,
        );

        for (high, boundary_high, incoming) in [
            (length_high5, Some(boundary_high5), borrow),
            (length_high6, Some(boundary_high6), high5_borrow),
            (length_high7, Some(boundary_high7), high6_borrow),
            (length_high8, None, high7_borrow),
        ] {
            if let Some(boundary_high) = boundary_high {
                circ.cx(boundary_high, high);
            }
            circ.cx(incoming, high);
        }
        for (difference_bit, output_bit) in length.iter().zip(&output[..length.len()]) {
            circ.ccx(enabled, difference_bit, output_bit);
        }
        for (offset, high) in [length_high5, length_high6, length_high7, length_high8]
            .into_iter()
            .enumerate()
        {
            circ.ccx(enabled, high, &output[length.len() + offset]);
        }
        for (high, boundary_high, incoming) in [
            (length_high8, None, high7_borrow),
            (length_high7, Some(boundary_high7), high6_borrow),
            (length_high6, Some(boundary_high6), high5_borrow),
            (length_high5, Some(boundary_high5), borrow),
        ] {
            circ.cx(incoming, high);
            if let Some(boundary_high) = boundary_high {
                circ.cx(boundary_high, high);
            }
        }

        toggle_controlled_subtraction_borrow_anf(
            circ,
            control,
            length_high8,
            None,
            high7_borrow,
            enabled,
            dirty,
        );
        circ.cx(control, enabled);
        toggle_subtraction_borrow_anf(
            circ,
            length_high7,
            boundary_high7,
            high6_borrow,
            high7_borrow,
        );
        toggle_subtraction_borrow_anf(
            circ,
            length_high6,
            boundary_high6,
            high5_borrow,
            high6_borrow,
        );
        toggle_subtraction_borrow_anf(circ, length_high5, boundary_high5, borrow, high5_borrow);
        cuccaro_add_mod_2n(circ, boundary_low, length, carry, borrow);
        return;
    }

    if split_three_high_length {
        assert!(borrowed);
        assert!(scratch.len() >= 8);
        assert_eq!(length.len() + 2, boundary.len());
        assert!(!length.is_empty());
        assert!(!source.is_empty());
        let boundary_low = &boundary[..length.len()];
        let boundary_high6 = &boundary[length.len()];
        let boundary_high7 = &boundary[length.len() + 1];
        let borrow = scratch[1];
        let high6_borrow = scratch[2];
        let high7_borrow = scratch[3];
        let length_high6 = scratch[5];
        let length_high7 = scratch[6];
        let length_high8 = scratch[7];
        let dirty = source[0];
        for (index, lane) in scratch[..8].iter().enumerate() {
            assert!(scratch[..index].iter().all(|other| other.id() != lane.id()));
            assert_ne!(lane.id(), control.id());
            assert!(boundary.iter().all(|other| other.id() != lane.id()));
            assert!(source.iter().all(|other| other.id() != lane.id()));
            assert!(output.iter().all(|other| other.id() != lane.id()));
            assert!(length.iter().all(|other| other.id() != lane.id()));
        }
        assert_ne!(dirty.id(), control.id());
        assert!(boundary.iter().all(|lane| lane.id() != dirty.id()));
        assert!(output.iter().all(|lane| lane.id() != dirty.id()));
        assert!(length.iter().all(|lane| lane.id() != dirty.id()));

        cuccaro_sub_mod_2n(circ, boundary_low, length, carry, borrow);
        toggle_subtraction_borrow_anf(
            circ,
            length_high6,
            boundary_high6,
            borrow,
            high6_borrow,
        );
        toggle_subtraction_borrow_anf(
            circ,
            length_high7,
            boundary_high7,
            high6_borrow,
            high7_borrow,
        );

        // The missing ninth boundary bit is zero, so br(h8,0,v)=v XOR h8*v.
        circ.cx(control, enabled);
        toggle_controlled_subtraction_borrow_anf(
            circ,
            control,
            length_high8,
            None,
            high7_borrow,
            enabled,
            dirty,
        );

        circ.cx(boundary_high6, length_high6);
        circ.cx(borrow, length_high6);
        circ.cx(boundary_high7, length_high7);
        circ.cx(high6_borrow, length_high7);
        circ.cx(high7_borrow, length_high8);
        for (difference_bit, output_bit) in length.iter().zip(&output[..length.len()]) {
            circ.ccx(enabled, difference_bit, output_bit);
        }
        circ.ccx(enabled, length_high6, &output[length.len()]);
        circ.ccx(enabled, length_high7, &output[length.len() + 1]);
        circ.ccx(enabled, length_high8, &output[length.len() + 2]);
        circ.cx(high7_borrow, length_high8);
        circ.cx(high6_borrow, length_high7);
        circ.cx(boundary_high7, length_high7);
        circ.cx(borrow, length_high6);
        circ.cx(boundary_high6, length_high6);

        toggle_controlled_subtraction_borrow_anf(
            circ,
            control,
            length_high8,
            None,
            high7_borrow,
            enabled,
            dirty,
        );
        circ.cx(control, enabled);
        toggle_subtraction_borrow_anf(
            circ,
            length_high7,
            boundary_high7,
            high6_borrow,
            high7_borrow,
        );
        toggle_subtraction_borrow_anf(
            circ,
            length_high6,
            boundary_high6,
            borrow,
            high6_borrow,
        );
        cuccaro_add_mod_2n(circ, boundary_low, length, carry, borrow);
        return;
    }

    if split_two_high_length {
        use crate::point_add::trailmix_port::arith::mcx::mcx_dirty_any_k;

        assert!(borrowed);
        assert!(scratch.len() >= 7);
        assert_eq!(length.len() + 1, boundary.len());
        assert!(!length.is_empty());
        let boundary_low = &boundary[..length.len()];
        let boundary_high7 = &boundary[length.len()];
        let borrow = scratch[1];
        let high_borrow = scratch[2];
        let length_high7 = scratch[5];
        let length_high8 = scratch[6];
        let dirty = &length[0];
        for (index, lane) in scratch[..7].iter().enumerate() {
            assert!(scratch[..index].iter().all(|other| other.id() != lane.id()));
            assert_ne!(lane.id(), control.id());
            assert!(boundary.iter().all(|other| other.id() != lane.id()));
            assert!(source.iter().all(|other| other.id() != lane.id()));
            assert!(output.iter().all(|other| other.id() != lane.id()));
            assert!(length.iter().all(|other| other.id() != lane.id()));
        }

        cuccaro_sub_mod_2n(circ, boundary_low, length, carry, borrow);

        // v = br(a,g,u). The missing ninth boundary bit is the constant zero.
        circ.cx(boundary_high7, high_borrow);
        circ.cx(borrow, high_borrow);
        circ.ccx(boundary_high7, borrow, high_borrow);
        circ.ccx(length_high7, boundary_high7, high_borrow);
        circ.ccx(length_high7, borrow, high_borrow);

        // U = br(h,0,v) = v XOR hv.
        circ.cx(control, enabled);
        circ.ccx(control, high_borrow, enabled);
        mcx_dirty_any_k(
            circ,
            &[control, length_high8, high_borrow],
            enabled,
            dirty,
        );

        circ.cx(boundary_high7, length_high7);
        circ.cx(borrow, length_high7);
        circ.cx(high_borrow, length_high8);
        for (difference_bit, output_bit) in length.iter().zip(&output[..length.len()]) {
            circ.ccx(enabled, difference_bit, output_bit);
        }
        circ.ccx(enabled, length_high7, &output[length.len()]);
        circ.ccx(enabled, length_high8, &output[length.len() + 1]);
        circ.cx(high_borrow, length_high8);
        circ.cx(borrow, length_high7);
        circ.cx(boundary_high7, length_high7);

        mcx_dirty_any_k(
            circ,
            &[control, length_high8, high_borrow],
            enabled,
            dirty,
        );
        circ.ccx(control, high_borrow, enabled);
        circ.cx(control, enabled);

        circ.ccx(length_high7, borrow, high_borrow);
        circ.ccx(length_high7, boundary_high7, high_borrow);
        circ.ccx(boundary_high7, borrow, high_borrow);
        circ.cx(borrow, high_borrow);
        circ.cx(boundary_high7, high_borrow);
        cuccaro_add_mod_2n(circ, boundary_low, length, carry, borrow);
        return;
    }

    if split_mixed_length {
        use crate::point_add::trailmix_port::arith::mcx::mcx_dirty_any_k;

        assert!(borrowed);
        assert_eq!(length.len(), boundary.len());
        assert!(!length.is_empty());
        let borrow = scratch[1];
        let high = scratch[2];
        let dirty = &length[0];

        // Write the low subtraction borrow into `borrow`. The full nine-bit
        // difference has high bit `high XOR borrow`, while underflow is
        // `borrow AND NOT high`. Materialize the enabled predicate in the
        // restored carry with one restored dirty lender.
        cuccaro_sub_mod_2n(circ, boundary, length, carry, borrow);
        circ.cx(control, enabled);
        circ.x(high);
        mcx_dirty_any_k(circ, &[control, borrow, high], enabled, dirty);
        circ.x(high);

        circ.cx(borrow, high);
        for (difference_bit, output_bit) in
            length.iter().zip(&output[..length.len()])
        {
            circ.ccx(enabled, difference_bit, output_bit);
        }
        circ.ccx(enabled, high, &output[output.len() - 1]);
        circ.cx(borrow, high);

        circ.x(high);
        mcx_dirty_any_k(circ, &[control, borrow, high], enabled, dirty);
        circ.x(high);
        circ.cx(control, enabled);
        cuccaro_add_mod_2n(circ, boundary, length, carry, borrow);
        return;
    }

    // The support theorem excludes 256, so the persistent unsigned boundary
    // needs only eight lanes. The prefix scratch is clean again at this point;
    // the opt-in route reuses one of those lanes as the ninth zero.
    let scratch_extension = sub800_mixed_boundary_scratch_extension_requested();
    assert!(
        !scratch_extension || borrowed,
        "mixed-boundary scratch extension requires three disjoint clean lenders"
    );
    let zero_extension_owned =
        (!scratch_extension).then(|| circ.alloc_qreg("rs.l-r-prime.zero-extension"));
    let zero_extension = if scratch_extension {
        scratch[1]
    } else {
        zero_extension_owned
            .as_ref()
            .expect("owned mixed-boundary zero extension")
    };
    let mut extended_boundary: Vec<&QReg> = boundary.iter().collect();
    extended_boundary.push(zero_extension);
    for lane in scratch.iter().take(3) {
        if !scratch_extension {
            assert_ne!(lane.id(), zero_extension.id());
        }
    }
    let (low_length, sign) = if borrowed_underflow {
        let sign = if scratch_extension {
            scratch[2]
        } else {
            scratch[1]
        };
        (length, sign)
    } else {
        let (low_length, sign_lane) = length.split_at(output.len());
        (low_length, &sign_lane[0])
    };
    cuccaro_sub_mod_2n_refs(circ, &extended_boundary, low_length, carry, sign);
    circ.x(sign);
    circ.ccx(control, sign, enabled);
    for (difference_bit, output_bit) in low_length.iter().zip(output) {
        circ.ccx(enabled, difference_bit, output_bit);
    }
    circ.ccx(control, sign, enabled);
    circ.x(sign);
    cuccaro_add_mod_2n_refs(circ, &extended_boundary, low_length, carry, sign);
    drop(extended_boundary);
    if let Some(zero_extension) = zero_extension_owned {
        circ.zero_and_free(zero_extension);
    }

    if let Some(enabled) = enabled_owned {
        circ.zero_and_free(enabled);
    }
    if let Some(carry) = carry_owned {
        circ.zero_and_free(carry);
    }
}

#[derive(Clone, Copy)]
enum SaturatingDifferenceBoundaryRoute {
    Configured,
    Materialized,
    Inplace,
}

fn assert_paired_bitlen_source_complement_preconditions(
    control: &QReg,
    boundary: &[QReg],
    source: &[&QReg],
    output: &[QReg],
    scratch: &[&QReg],
    borrowed_zero_correction_carry: Option<&QReg>,
) {
    use super::shrunken_pz_state_machine::lowq_fused_zero_prefix_bitlen_requested;

    assert!(
        lowq_fused_zero_prefix_bitlen_requested(),
        "paired source complement requires fused zero-prefix bit length"
    );
    assert!(
        source.len() > 1,
        "paired source complement requires a source wider than one bit"
    );
    for (index, source_lane) in source.iter().enumerate() {
        assert!(
            source[..index]
                .iter()
                .all(|other| other.id() != source_lane.id()),
            "paired source-complement lane {index} aliases an earlier source lane"
        );
        assert_ne!(
            source_lane.id(),
            control.id(),
            "paired source-complement lane {index} aliases the control"
        );
        assert!(
            boundary
                .iter()
                .all(|other| other.id() != source_lane.id()),
            "paired source-complement lane {index} aliases the boundary"
        );
        assert!(
            output
                .iter()
                .all(|other| other.id() != source_lane.id()),
            "paired source-complement lane {index} aliases the output"
        );
        assert!(
            scratch
                .iter()
                .all(|other| other.id() != source_lane.id()),
            "paired source-complement lane {index} aliases borrowed scratch"
        );
        if let Some(carry) = borrowed_zero_correction_carry {
            assert_ne!(
                source_lane.id(),
                carry.id(),
                "paired source-complement lane {index} aliases the borrowed carry"
            );
        }
    }
}

fn toggle_split_bit_length_high(
    circ: &mut Circuit,
    source: &[&QReg],
    output_width: usize,
    high: &QReg,
    dirty: &QReg,
    source_is_complemented: bool,
) {
    use crate::point_add::trailmix_port::arith::mcx::mcx_dirty_any_k;

    assert!(output_width >= 2);
    let threshold = 1usize << (output_width - 1);
    if source.len() < threshold {
        return;
    }
    let high_source = &source[threshold - 1..];
    assert!(!high_source.is_empty());
    assert!(high_source.iter().all(|lane| lane.id() != high.id()));
    assert!(high_source.iter().all(|lane| lane.id() != dirty.id()));
    assert_ne!(high.id(), dirty.id());

    if !source_is_complemented {
        for lane in high_source {
            circ.x(lane);
        }
    }
    circ.x(high);
    mcx_dirty_any_k(circ, high_source, high, dirty);
    if !source_is_complemented {
        for lane in high_source {
            circ.x(lane);
        }
    }
}

fn toggle_split_bit_length_two_high(
    circ: &mut Circuit,
    source: &[&QReg],
    high7: &QReg,
    high8: &QReg,
    scratch: &[&QReg],
    source_is_complemented: bool,
) {
    use super::shrunken_pz_state_machine::DIRECT_PREFIX_KG_SCRATCH_LEN;
    use crate::point_add::trailmix_port::arith::khattar_gidney::{
        kg_prefix_ancilla_count, KgPrefixAnd,
    };

    assert_eq!(source.len(), 259);
    assert!(scratch.len() >= DIRECT_PREFIX_KG_SCRATCH_LEN + 2);
    assert_ne!(high7.id(), high8.id());
    assert!(source.iter().all(|lane| lane.id() != high7.id()));
    assert!(source.iter().all(|lane| lane.id() != high8.id()));
    let anc = &scratch[..DIRECT_PREFIX_KG_SCRATCH_LEN];
    assert!(anc.iter().all(|lane| lane.id() != high7.id()));
    assert!(anc.iter().all(|lane| lane.id() != high8.id()));
    assert!(kg_prefix_ancilla_count(source.len()) <= anc.len());

    if !source_is_complemented {
        for lane in source {
            circ.x(lane);
        }
    }
    let qbits: Vec<&QReg> = source.iter().rev().copied().collect();
    let z128_prefix = source.len() - 127;
    let z256_prefix = source.len() - 255;

    // B[7] = [B<128] XOR [B<256], while B[8] = 1 XOR [B<256].
    circ.x(high8);
    let done = KgPrefixAnd::new(&qbits, anc).forward(circ, |_, _, _| {});
    done.reverse(circ, |circ, prefix_len, controls| {
        let toggle = |circ: &mut Circuit, target: &QReg| match controls {
            [control] => circ.cx(control, target),
            [left, right] => circ.ccx(left, right, target),
            _ => unreachable!("KG prefix controls must contain one or two qubits"),
        };
        if prefix_len == z128_prefix {
            toggle(circ, high7);
        }
        if prefix_len == z256_prefix {
            toggle(circ, high7);
            toggle(circ, high8);
        }
    });

    if !source_is_complemented {
        for lane in source {
            circ.x(lane);
        }
    }
}

fn toggle_split_bit_length_three_high(
    circ: &mut Circuit,
    source: &[&QReg],
    high6: &QReg,
    high7: &QReg,
    high8: &QReg,
    scratch: &[&QReg],
    source_is_complemented: bool,
) {
    use super::shrunken_pz_state_machine::DIRECT_PREFIX_KG_SCRATCH_LEN;
    use crate::point_add::trailmix_port::arith::khattar_gidney::{
        kg_prefix_ancilla_count, KgPrefixAnd,
    };

    const SOURCE_WIDTH: usize = 259;
    const Z64_PREFIX_LEN: usize = 196;
    const Z128_PREFIX_LEN: usize = 132;
    const Z192_PREFIX_LEN: usize = 68;
    const Z256_PREFIX_LEN: usize = 4;

    assert_eq!(source.len(), SOURCE_WIDTH);
    assert!(scratch.len() >= DIRECT_PREFIX_KG_SCRATCH_LEN + 3);
    assert_ne!(high6.id(), high7.id());
    assert_ne!(high6.id(), high8.id());
    assert_ne!(high7.id(), high8.id());
    for high in [high6, high7, high8] {
        assert!(source.iter().all(|lane| lane.id() != high.id()));
    }
    let anc = &scratch[..DIRECT_PREFIX_KG_SCRATCH_LEN];
    for high in [high6, high7, high8] {
        assert!(anc.iter().all(|lane| lane.id() != high.id()));
    }
    assert!(kg_prefix_ancilla_count(source.len()) <= anc.len());

    if !source_is_complemented {
        for lane in source {
            circ.x(lane);
        }
    }
    let qbits: Vec<&QReg> = source.iter().rev().copied().collect();

    // b6=Z64^Z128^Z192^Z256, b7=Z128^Z256, b8=1^Z256.
    circ.x(high8);
    let done = KgPrefixAnd::new(&qbits, anc).forward(circ, |_, _, _| {});
    done.reverse(circ, |circ, prefix_len, controls| {
        let toggle = |circ: &mut Circuit, target: &QReg| match controls {
            [control] => circ.cx(control, target),
            [left, right] => circ.ccx(left, right, target),
            _ => unreachable!("KG prefix controls must contain one or two qubits"),
        };
        match prefix_len {
            Z64_PREFIX_LEN => toggle(circ, high6),
            Z128_PREFIX_LEN => {
                toggle(circ, high6);
                toggle(circ, high7);
            }
            Z192_PREFIX_LEN => toggle(circ, high6),
            Z256_PREFIX_LEN => {
                toggle(circ, high6);
                toggle(circ, high7);
                toggle(circ, high8);
            }
            _ => {}
        }
    });

    if !source_is_complemented {
        for lane in source {
            circ.x(lane);
        }
    }
}

fn toggle_split_bit_length_four_high(
    circ: &mut Circuit,
    source: &[&QReg],
    high5: &QReg,
    high6: &QReg,
    high7: &QReg,
    high8: &QReg,
    scratch: &[&QReg],
    source_is_complemented: bool,
) {
    use super::shrunken_pz_state_machine::DIRECT_PREFIX_KG_SCRATCH_LEN;
    use crate::point_add::trailmix_port::arith::khattar_gidney::{
        kg_prefix_ancilla_count, KgPrefixAnd,
    };

    const SOURCE_WIDTH: usize = 259;
    const THRESHOLDS: [(usize, usize); 8] = [
        (32, 228),
        (64, 196),
        (96, 164),
        (128, 132),
        (160, 100),
        (192, 68),
        (224, 36),
        (256, 4),
    ];

    assert_eq!(source.len(), SOURCE_WIDTH);
    assert!(scratch.len() >= DIRECT_PREFIX_KG_SCRATCH_LEN + 4);
    let highs = [high5, high6, high7, high8];
    for (index, high) in highs.iter().enumerate() {
        assert!(highs[..index].iter().all(|other| other.id() != high.id()));
        assert!(source.iter().all(|lane| lane.id() != high.id()));
    }
    let anc = &scratch[..DIRECT_PREFIX_KG_SCRATCH_LEN];
    for high in highs {
        assert!(anc.iter().all(|lane| lane.id() != high.id()));
    }
    assert!(kg_prefix_ancilla_count(source.len()) <= anc.len());

    if !source_is_complemented {
        for lane in source {
            circ.x(lane);
        }
    }
    let qbits: Vec<&QReg> = source.iter().rev().copied().collect();
    circ.x(high8);
    let done = KgPrefixAnd::new(&qbits, anc).forward(circ, |_, _, _| {});
    done.reverse(circ, |circ, prefix_len, controls| {
        let toggle = |circ: &mut Circuit, target: &QReg| match controls {
            [control] => circ.cx(control, target),
            [left, right] => circ.ccx(left, right, target),
            _ => unreachable!("KG prefix controls must contain one or two qubits"),
        };
        let threshold = THRESHOLDS
            .iter()
            .find_map(|(threshold, length)| (*length == prefix_len).then_some(*threshold));
        if let Some(threshold) = threshold {
            toggle(circ, high5);
            if threshold % 64 == 0 {
                toggle(circ, high6);
            }
            if threshold % 128 == 0 {
                toggle(circ, high7);
            }
            if threshold == 256 {
                toggle(circ, high8);
            }
        }
    });

    if !source_is_complemented {
        for lane in source {
            circ.x(lane);
        }
    }
}

fn controlled_xor_saturating_bit_length_difference_with_route(
    circ: &mut Circuit,
    control: &QReg,
    boundary: &[QReg],
    source: &[&QReg],
    output: &[QReg],
    scratch: &[&QReg],
    borrowed_zero_correction_carry: Option<&QReg>,
    route: SaturatingDifferenceBoundaryRoute,
) {
    assert!(boundary.len() == output.len() || boundary.len() + 1 == output.len());
    assert!(!output.is_empty());
    let signed_width = output.len() + 1;
    assert!(
        source.len() <= (1usize << (signed_width - 1)) - 1,
        "signed bit-length workspace is too narrow"
    );
    let inplace = match route {
        SaturatingDifferenceBoundaryRoute::Configured => {
            inplace_rotated_bitlen_boundary_requested()
        }
        SaturatingDifferenceBoundaryRoute::Materialized => false,
        SaturatingDifferenceBoundaryRoute::Inplace => true,
    };
    let borrowed_underflow = inplace && sub800_borrowed_rotated_underflow_requested();
    let production_split_shape = source.len() == 259 && output.len() == 9;
    let split_four_high_length = borrowed_underflow
        && sub800_split_four_high_rotated_length_requested()
        && production_split_shape;
    let split_three_high_length = !split_four_high_length
        && borrowed_underflow
        && sub800_split_three_high_rotated_length_requested()
        && production_split_shape;
    let split_two_high_length = !split_four_high_length
        && !split_three_high_length
        && borrowed_underflow
        && sub800_split_two_high_rotated_length_requested()
        && production_split_shape;
    let split_mixed_length = !split_four_high_length
        && !split_three_high_length
        && !split_two_high_length
        && borrowed_underflow
        && boundary.len() + 1 == output.len()
        && sub800_split_mixed_rotated_length_requested();
    let split_same_length = !split_four_high_length
        && !split_three_high_length
        && !split_two_high_length
        && borrowed_underflow
        && boundary.len() == output.len()
        && sub800_split_same_rotated_length_requested();
    let split_length = split_four_high_length
        || split_three_high_length
        || split_two_high_length
        || split_mixed_length
        || split_same_length;
    if split_four_high_length {
        if boundary.len() == output.len() {
            SUB800_Q838_SPLIT_FOUR_SAME_CALLS.fetch_add(1, Ordering::Relaxed);
        } else {
            assert_eq!(boundary.len() + 1, output.len());
            SUB800_Q838_SPLIT_FOUR_MIXED_CALLS.fetch_add(1, Ordering::Relaxed);
        }
    } else if split_three_high_length {
        if boundary.len() == output.len() {
            SUB800_Q839_SPLIT_THREE_SAME_CALLS.fetch_add(1, Ordering::Relaxed);
        } else {
            assert_eq!(boundary.len() + 1, output.len());
            SUB800_Q839_SPLIT_THREE_MIXED_CALLS.fetch_add(1, Ordering::Relaxed);
        }
    }
    assert!(
        !borrowed_underflow || scratch.len() >= 3,
        "borrowed rotated underflow requires three clean scratch lanes"
    );
    assert!(
        !split_length || output.len() >= 2,
        "split rotated length requires at least two output lanes"
    );
    assert!(
        !split_four_high_length
            || (signed_width == 10 && scratch.len() >= 9 && production_split_shape),
        "split-four-high rotated length requires signed width ten, nine scratch lanes, and a 259-bit source"
    );
    assert!(
        !split_three_high_length
            || (signed_width == 10 && scratch.len() >= 8 && production_split_shape),
        "split-three-high rotated length requires signed width ten, eight scratch lanes, and a 259-bit source"
    );
    assert!(
        !split_two_high_length
            || (scratch.len() >= 7 && production_split_shape),
        "split-two-high rotated length requires nine output lanes, seven scratch lanes, and a 259-bit source"
    );
    let omitted_split_bits = if split_four_high_length {
        4
    } else if split_three_high_length {
        3
    } else if split_two_high_length {
        2
    } else {
        usize::from(split_length)
    };
    let length_width = signed_width - usize::from(borrowed_underflow) - omitted_split_bits;

    // Baseline is X_S K_add X_S V X_S K_sub X_S. The intervening USE block
    // V is source-disjoint, so the middle X_S pair commutes through V and
    // cancels. Keep S complemented across the pair and emit only the outer
    // brackets: X_S K_add V K_sub X_S.
    let paired_source_complement = paired_bitlen_source_complement_requested();
    if paired_source_complement {
        assert_paired_bitlen_source_complement_preconditions(
            control,
            boundary,
            source,
            output,
            scratch,
            borrowed_zero_correction_carry,
        );
        for lane in source {
            circ.x(lane);
        }
    }

    let length = circ.alloc_qreg_bits("rs.rotated-bitlen.length", length_width);
    if split_four_high_length {
        bit_length_lean_allow_zero_with_borrowed_scratch_split_four_high(
            circ,
            source,
            &length,
            false,
            borrowed_zero_correction_carry,
            scratch,
            paired_source_complement,
        );
    } else if split_three_high_length {
        bit_length_lean_allow_zero_with_borrowed_scratch_split_three_high(
            circ,
            source,
            &length,
            false,
            borrowed_zero_correction_carry,
            scratch,
            paired_source_complement,
        );
    } else if split_two_high_length {
        bit_length_lean_allow_zero_with_borrowed_scratch_split_two_high(
            circ,
            source,
            &length,
            false,
            borrowed_zero_correction_carry,
            scratch,
            paired_source_complement,
        );
    } else if split_length {
        bit_length_lean_allow_zero_with_borrowed_scratch_split_high(
            circ,
            source,
            &length,
            false,
            borrowed_zero_correction_carry,
            scratch,
            paired_source_complement,
        );
    } else if paired_source_complement {
        bit_length_lean_allow_zero_with_borrowed_scratch_complemented_source(
            circ,
            source,
            &length,
            false,
            borrowed_zero_correction_carry,
            scratch,
        );
    } else {
        bit_length_lean_allow_zero_with_borrowed_scratch(
            circ,
            source,
            &length,
            false,
            borrowed_zero_correction_carry,
            scratch,
        );
    }
    if split_four_high_length {
        toggle_split_bit_length_four_high(
            circ,
            source,
            scratch[5],
            scratch[6],
            scratch[7],
            scratch[8],
            scratch,
            paired_source_complement,
        );
    } else if split_three_high_length {
        toggle_split_bit_length_three_high(
            circ,
            source,
            scratch[5],
            scratch[6],
            scratch[7],
            scratch,
            paired_source_complement,
        );
    } else if split_two_high_length {
        toggle_split_bit_length_two_high(
            circ,
            source,
            scratch[5],
            scratch[6],
            scratch,
            paired_source_complement,
        );
    } else if split_length {
        toggle_split_bit_length_high(
            circ,
            source,
            output.len(),
            scratch[2],
            source[0],
            paired_source_complement,
        );
    }

    if inplace {
        if boundary.len() == output.len() {
            controlled_xor_saturating_difference_inplace_boundary(
                circ, control, boundary, source, &length, output, scratch,
            );
        } else {
            controlled_xor_saturating_difference_inplace_mixed_boundary(
                circ, control, boundary, source, &length, output, scratch,
            );
        }
    } else {
        controlled_xor_saturating_difference_materialized_boundary(
            circ, control, boundary, source, &length, output, scratch,
        );
    }

    if split_four_high_length {
        toggle_split_bit_length_four_high(
            circ,
            source,
            scratch[5],
            scratch[6],
            scratch[7],
            scratch[8],
            scratch,
            paired_source_complement,
        );
    } else if split_three_high_length {
        toggle_split_bit_length_three_high(
            circ,
            source,
            scratch[5],
            scratch[6],
            scratch[7],
            scratch,
            paired_source_complement,
        );
    } else if split_two_high_length {
        toggle_split_bit_length_two_high(
            circ,
            source,
            scratch[5],
            scratch[6],
            scratch,
            paired_source_complement,
        );
    } else if split_length {
        toggle_split_bit_length_high(
            circ,
            source,
            output.len(),
            scratch[2],
            source[0],
            paired_source_complement,
        );
    }
    if split_four_high_length {
        bit_length_lean_allow_zero_with_borrowed_scratch_split_four_high(
            circ,
            source,
            &length,
            true,
            borrowed_zero_correction_carry,
            scratch,
            paired_source_complement,
        );
        if paired_source_complement {
            for lane in source {
                circ.x(lane);
            }
        }
    } else if split_three_high_length {
        bit_length_lean_allow_zero_with_borrowed_scratch_split_three_high(
            circ,
            source,
            &length,
            true,
            borrowed_zero_correction_carry,
            scratch,
            paired_source_complement,
        );
        if paired_source_complement {
            for lane in source {
                circ.x(lane);
            }
        }
    } else if split_two_high_length {
        bit_length_lean_allow_zero_with_borrowed_scratch_split_two_high(
            circ,
            source,
            &length,
            true,
            borrowed_zero_correction_carry,
            scratch,
            paired_source_complement,
        );
        if paired_source_complement {
            for lane in source {
                circ.x(lane);
            }
        }
    } else if split_length {
        bit_length_lean_allow_zero_with_borrowed_scratch_split_high(
            circ,
            source,
            &length,
            true,
            borrowed_zero_correction_carry,
            scratch,
            paired_source_complement,
        );
        if paired_source_complement {
            for lane in source {
                circ.x(lane);
            }
        }
    } else if paired_source_complement {
        bit_length_lean_allow_zero_with_borrowed_scratch_complemented_source(
            circ,
            source,
            &length,
            true,
            borrowed_zero_correction_carry,
            scratch,
        );
        for lane in source {
            circ.x(lane);
        }
    } else {
        bit_length_lean_allow_zero_with_borrowed_scratch(
            circ,
            source,
            &length,
            true,
            borrowed_zero_correction_carry,
            scratch,
        );
    }
    for lane in length {
        circ.zero_and_free(lane);
    }
}

fn controlled_xor_saturating_bit_length_difference(
    circ: &mut Circuit,
    control: &QReg,
    boundary: &[QReg],
    source: &[&QReg],
    output: &[QReg],
    scratch: &[&QReg],
    borrowed_zero_correction_carry: Option<&QReg>,
) {
    controlled_xor_saturating_bit_length_difference_with_route(
        circ,
        control,
        boundary,
        source,
        output,
        scratch,
        borrowed_zero_correction_carry,
        SaturatingDifferenceBoundaryRoute::Configured,
    );
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
    scratch: &[QReg],
    borrowed_prefix_scratch: &[&QReg],
) {
    assert!(borrowed_prefix_scratch.len() <= 2);
    for lender in borrowed_prefix_scratch {
        assert_ne!(lender.id(), control.id());
        assert!(right_length.iter().all(|lane| lane.id() != lender.id()));
        assert!(work.iter().all(|lane| lane.id() != lender.id()));
        assert!(output.iter().all(|lane| lane.id() != lender.id()));
        assert!(scratch.iter().all(|lane| lane.id() != lender.id()));
    }
    let view: Vec<&QReg> = work.iter().collect();
    let scratch_refs: Vec<&QReg> = scratch
        .iter()
        .chain(borrowed_prefix_scratch.iter().copied())
        .collect();
    variable_rotate_high_refs(circ, right_length, &view);
    controlled_xor_saturating_bit_length_difference(
        circ,
        control,
        right_length,
        &view,
        output,
        &scratch_refs,
        None,
    );
    variable_rotate_low_refs(circ, right_length, &view);
}

#[allow(clippy::too_many_arguments)]
fn controlled_xor_rotated_prefix_with_predicate_lenders(
    circ: &mut Circuit,
    zero_q: &QReg,
    zero_s: &QReg,
    control: &QReg,
    l_q: &[QReg],
    l_s: &[QReg],
    right_length: &[QReg],
    work: &[QReg],
    output: &[QReg],
    scratch: &[QReg],
) {
    assert!(scratch.len() >= l_q.len().saturating_sub(2));
    circ.ccx(zero_q, zero_s, control);
    uncompute_zero(circ, l_s, zero_s, scratch);
    uncompute_zero(circ, l_q, zero_q, scratch);
    controlled_xor_rotated_prefix_bit_length(
        circ,
        control,
        right_length,
        work,
        output,
        scratch,
        &[zero_q, zero_s],
    );
    compute_zero(circ, l_q, zero_q, scratch);
    compute_zero(circ, l_s, zero_s, scratch);
    circ.ccx(zero_q, zero_s, control);
}

/// XOR the raw `t'` bit length into `output` while Work2 is in its shifted
/// coefficient-update layout. The physical rotation includes `l_s`, but the
/// packed remainder boundary does not, so those two lengths are intentionally
/// distinct.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum RawTPrimeRotationLifetime {
    Continuous,
    Split,
}

fn controlled_xor_raw_t_prime_bit_length_allocated_with_lifetime(
    circ: &mut Circuit,
    control: &QReg,
    l_s: &[QReg],
    l_r_prime: &[QReg],
    work2: &[QReg],
    output: &[QReg],
    lifetime: RawTPrimeRotationLifetime,
) {
    assert_eq!(l_s.len(), output.len());
    assert_l_r_prime_metadata_width(l_s.len(), l_r_prime.len());
    let mut rotation = Some(circ.alloc_qreg_bits("rs.raw-t-prime.rotation", output.len()));
    for (source, destination) in l_r_prime
        .iter()
        .zip(rotation.as_ref().expect("live raw rotation"))
    {
        circ.cx(source, destination);
    }
    trace_raw_bit_length_allocations(RawBitLengthAllocationTrace {
        raw_rotation_carry_allocations: 1,
        raw_rotation_overflow_allocations: 1,
        ..RawBitLengthAllocationTrace::default()
    });
    let carry = circ.alloc_qreg("rs.raw-t-prime.rotation-carry");
    let overflow = circ.alloc_qreg("rs.raw-t-prime.rotation-overflow");
    cuccaro_add_mod_2n(
        circ,
        l_s,
        rotation.as_ref().expect("live raw rotation"),
        &carry,
        &overflow,
    );

    let view: Vec<&QReg> = work2.iter().collect();
    variable_rotate_high_refs(circ, rotation.as_ref().expect("live raw rotation"), &view);
    if lifetime == RawTPrimeRotationLifetime::Split {
        let released = rotation.take().expect("split raw rotation release");
        cuccaro_sub_mod_2n(circ, l_s, &released, &carry, &overflow);
        for (source, destination) in l_r_prime.iter().zip(&released) {
            circ.cx(source, destination);
        }
        trace_raw_bit_length_allocations(RawBitLengthAllocationTrace {
            raw_rotation_split_releases: 1,
            raw_rotation_lanes_released: released.len(),
            ..RawBitLengthAllocationTrace::default()
        });
        free_clean(circ, released);
    }
    controlled_xor_saturating_bit_length_difference(
        circ,
        control,
        l_r_prime,
        &view,
        output,
        &[],
        None,
    );
    if lifetime == RawTPrimeRotationLifetime::Split {
        let recomputed = circ.alloc_qreg_bits("rs.raw-t-prime.rotation", output.len());
        for (source, destination) in l_r_prime.iter().zip(&recomputed) {
            circ.cx(source, destination);
        }
        cuccaro_add_mod_2n(circ, l_s, &recomputed, &carry, &overflow);
        trace_raw_bit_length_allocations(RawBitLengthAllocationTrace {
            raw_rotation_split_recomputes: 1,
            ..RawBitLengthAllocationTrace::default()
        });
        rotation = Some(recomputed);
    }
    variable_rotate_low_refs(circ, rotation.as_ref().expect("live raw rotation"), &view);

    let rotation = rotation.take().expect("final raw rotation release");
    cuccaro_sub_mod_2n(circ, l_s, &rotation, &carry, &overflow);
    circ.zero_and_free(overflow);
    circ.zero_and_free(carry);
    for (source, destination) in l_r_prime.iter().zip(&rotation) {
        circ.cx(source, destination);
    }
    free_clean(circ, rotation);
}

fn controlled_xor_raw_t_prime_bit_length_allocated(
    circ: &mut Circuit,
    control: &QReg,
    l_s: &[QReg],
    l_r_prime: &[QReg],
    work2: &[QReg],
    output: &[QReg],
) {
    controlled_xor_raw_t_prime_bit_length_allocated_with_lifetime(
        circ,
        control,
        l_s,
        l_r_prime,
        work2,
        output,
        RawTPrimeRotationLifetime::Continuous,
    );
}

fn assert_coefficient_raw_bitlen_loan_preconditions(
    control: &QReg,
    l_s: &[QReg],
    l_r_prime: &[QReg],
    work2: &[QReg],
    output: &[QReg],
    coefficient_chain: &[QReg],
) {
    assert_eq!(l_s.len(), output.len());
    assert_l_r_prime_metadata_width(l_s.len(), l_r_prime.len());
    assert!(!output.is_empty());
    assert_eq!(
        coefficient_chain.len(),
        2,
        "coefficient raw bit-length loan requires exactly two chain lanes"
    );

    let mut ids = Vec::with_capacity(
        1 + l_s.len() + l_r_prime.len() + work2.len() + output.len() + coefficient_chain.len(),
    );
    for (index, lane) in std::iter::once(control)
        .chain(l_s)
        .chain(l_r_prime)
        .chain(work2)
        .chain(output)
        .chain(coefficient_chain)
        .enumerate()
    {
        assert!(
            !ids.contains(&lane.id()),
            "coefficient raw bit-length loan lane {index} aliases an operand or lender"
        );
        ids.push(lane.id());
    }
}

fn controlled_xor_raw_t_prime_bit_length_loaned_with_lifetime(
    circ: &mut Circuit,
    control: &QReg,
    l_s: &[QReg],
    l_r_prime: &[QReg],
    work2: &[QReg],
    output: &[QReg],
    coefficient_chain: &[QReg],
    lifetime: RawTPrimeRotationLifetime,
) {
    assert_coefficient_raw_bitlen_loan_preconditions(
        control,
        l_s,
        l_r_prime,
        work2,
        output,
        coefficient_chain,
    );
    let mut rotation = Some(circ.alloc_qreg_bits("rs.raw-t-prime.rotation", output.len()));
    for (source, destination) in l_r_prime
        .iter()
        .zip(rotation.as_ref().expect("live raw rotation"))
    {
        circ.cx(source, destination);
    }
    trace_raw_bit_length_allocations(RawBitLengthAllocationTrace {
        raw_rotation_carry_allocations: 1,
        ..RawBitLengthAllocationTrace::default()
    });
    let carry = circ.alloc_qreg("rs.raw-t-prime.rotation-carry");
    cuccaro_add_mod_2n_no_overflow(
        circ,
        l_s,
        rotation.as_ref().expect("live raw rotation"),
        &carry,
    );

    let view: Vec<&QReg> = work2.iter().collect();
    let inner_scratch = vec![&coefficient_chain[0], &coefficient_chain[1], &carry];
    variable_rotate_high_refs(circ, rotation.as_ref().expect("live raw rotation"), &view);
    if lifetime == RawTPrimeRotationLifetime::Split {
        let released = rotation.take().expect("split raw rotation release");
        cuccaro_sub_mod_2n_no_overflow(circ, l_s, &released, &carry);
        for (source, destination) in l_r_prime.iter().zip(&released) {
            circ.cx(source, destination);
        }
        trace_raw_bit_length_allocations(RawBitLengthAllocationTrace {
            raw_rotation_split_releases: 1,
            raw_rotation_lanes_released: released.len(),
            ..RawBitLengthAllocationTrace::default()
        });
        free_clean(circ, released);
    }
    controlled_xor_saturating_bit_length_difference(
        circ,
        control,
        l_r_prime,
        &view,
        output,
        &inner_scratch,
        Some(&carry),
    );
    if lifetime == RawTPrimeRotationLifetime::Split {
        let recomputed = circ.alloc_qreg_bits("rs.raw-t-prime.rotation", output.len());
        for (source, destination) in l_r_prime.iter().zip(&recomputed) {
            circ.cx(source, destination);
        }
        cuccaro_add_mod_2n_no_overflow(circ, l_s, &recomputed, &carry);
        trace_raw_bit_length_allocations(RawBitLengthAllocationTrace {
            raw_rotation_split_recomputes: 1,
            ..RawBitLengthAllocationTrace::default()
        });
        rotation = Some(recomputed);
    }
    variable_rotate_low_refs(circ, rotation.as_ref().expect("live raw rotation"), &view);

    let rotation = rotation.take().expect("final raw rotation release");
    cuccaro_sub_mod_2n_no_overflow(circ, l_s, &rotation, &carry);
    drop(inner_scratch);
    circ.zero_and_free(carry);
    for (source, destination) in l_r_prime.iter().zip(&rotation) {
        circ.cx(source, destination);
    }
    free_clean(circ, rotation);
}

fn controlled_xor_raw_t_prime_bit_length_loaned(
    circ: &mut Circuit,
    control: &QReg,
    l_s: &[QReg],
    l_r_prime: &[QReg],
    work2: &[QReg],
    output: &[QReg],
    coefficient_chain: &[QReg],
) {
    controlled_xor_raw_t_prime_bit_length_loaned_with_lifetime(
        circ,
        control,
        l_s,
        l_r_prime,
        work2,
        output,
        coefficient_chain,
        RawTPrimeRotationLifetime::Continuous,
    );
}

fn controlled_xor_raw_t_prime_bit_length(
    circ: &mut Circuit,
    control: &QReg,
    l_s: &[QReg],
    l_r_prime: &[QReg],
    work2: &[QReg],
    output: &[QReg],
    coefficient_chain: &[QReg],
) {
    let lifetime = if split_coefficient_rotation_lifetime_requested() {
        RawTPrimeRotationLifetime::Split
    } else {
        RawTPrimeRotationLifetime::Continuous
    };
    if coefficient_raw_bitlen_loan_requested() {
        controlled_xor_raw_t_prime_bit_length_loaned_with_lifetime(
            circ,
            control,
            l_s,
            l_r_prime,
            work2,
            output,
            coefficient_chain,
            lifetime,
        );
    } else {
        controlled_xor_raw_t_prime_bit_length_allocated_with_lifetime(
            circ, control, l_s, l_r_prime, work2, output, lifetime,
        );
    }
}

/// In the reversed Work1 view, the packed `t` component occupies the
/// rightmost `left_length` lanes. The separator is zero and remains above the
/// rotated remainder, hence `bitlen(rotated) - left_length = bitlen(r)`.
fn controlled_xor_rotated_suffix_bit_length_with_prefix_scratch(
    circ: &mut Circuit,
    control: &QReg,
    left_length: &[QReg],
    work: &[QReg],
    output: &[QReg],
    scratch: &[QReg],
    borrowed_prefix_scratch: &[&QReg],
) {
    assert!(borrowed_prefix_scratch.len() <= 3);
    for (index, lender) in borrowed_prefix_scratch.iter().enumerate() {
        assert!(
            borrowed_prefix_scratch[..index]
                .iter()
                .all(|other| other.id() != lender.id())
        );
        assert_ne!(lender.id(), control.id());
        assert!(left_length.iter().all(|lane| lane.id() != lender.id()));
        assert!(work.iter().all(|lane| lane.id() != lender.id()));
        assert!(output.iter().all(|lane| lane.id() != lender.id()));
        assert!(scratch.iter().all(|lane| lane.id() != lender.id()));
    }
    let reversed: Vec<&QReg> = work.iter().rev().collect();
    let scratch_refs: Vec<&QReg> = scratch
        .iter()
        .chain(borrowed_prefix_scratch.iter().copied())
        .collect();
    variable_rotate_high_refs(circ, left_length, &reversed);
    controlled_xor_saturating_bit_length_difference(
        circ,
        control,
        left_length,
        &reversed,
        output,
        &scratch_refs,
        None,
    );
    variable_rotate_low_refs(circ, left_length, &reversed);
}

fn controlled_xor_rotated_suffix_bit_length(
    circ: &mut Circuit,
    control: &QReg,
    left_length: &[QReg],
    work: &[QReg],
    output: &[QReg],
    scratch: &[QReg],
    preserved_dy_top_scratch: &[&QReg],
) {
    assert!(preserved_dy_top_scratch.len() <= 1);
    let entry_ops_idx = circ.total_ops() as usize;
    controlled_xor_rotated_suffix_bit_length_with_prefix_scratch(
        circ,
        control,
        left_length,
        work,
        output,
        scratch,
        preserved_dy_top_scratch,
    );
    if let Some(lender) = preserved_dy_top_scratch.first() {
        record_preserved_dy_top_borrow_window(lender, entry_ops_idx, circ.total_ops() as usize);
    }
}

#[allow(clippy::too_many_arguments)]
fn controlled_xor_rotated_suffix_with_zero_s_lender(
    circ: &mut Circuit,
    zero_s: &QReg,
    clean_lender: &QReg,
    l_s: &[QReg],
    control: &QReg,
    left_length: &[QReg],
    work: &[QReg],
    output: &[QReg],
    scratch: &[QReg],
    additional_prefix_scratch: &[&QReg],
) {
    assert!(
        scratch.len() + additional_prefix_scratch.len() >= l_s.len().saturating_sub(2)
    );
    let mut zero_predicate_scratch: Vec<QReg> =
        scratch.iter().map(QReg::borrowed_alias).collect();
    if zero_predicate_scratch.len() < l_s.len().saturating_sub(2) {
        zero_predicate_scratch.push(
            output
                .last()
                .expect("hosted high output lane")
                .borrowed_alias(),
        );
    }
    uncompute_zero(circ, l_s, zero_s, &zero_predicate_scratch);
    let borrowed_prefix_scratch: Vec<&QReg> = [zero_s, clean_lender]
        .into_iter()
        .chain(additional_prefix_scratch.iter().copied())
        .collect();
    controlled_xor_rotated_suffix_bit_length_with_prefix_scratch(
        circ,
        control,
        left_length,
        work,
        output,
        scratch,
        &borrowed_prefix_scratch,
    );
    compute_zero(circ, l_s, zero_s, &zero_predicate_scratch);
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
    assert_l_r_prime_metadata_width(l_t.len(), l_r_prime.len());
    let old_r_length = circ.alloc_qreg_bits("rs.swap-length.old-lrp", l_t.len());

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
    scratch: &[QReg],
) {
    assert_eq!(work1.len(), work2.len());
    assert_l_r_prime_metadata_width(l_t.len(), l_r_prime.len());
    assert_eq!(l_t.len(), l_t_prime.len());
    let old_r_length = circ.alloc_qreg_bits("rs.swap-length.old-lrp", l_t.len());

    with_bit_length_callsite(
        "p.bitlen.rs.swap-old-r.pre.deposit",
        "p.bitlen.rs.swap-old-r.pre.erase",
        || {
            controlled_xor_rotated_suffix_bit_length(
                circ,
                control,
                l_t,
                work1,
                &old_r_length,
                scratch,
                &[],
            );
        },
    );

    for (current, next) in l_t.iter().zip(l_t_prime) {
        circ.cswap(control, current, next);
    }
    for (current, next) in l_r_prime.iter().zip(&old_r_length) {
        circ.cswap(control, current, next);
    }
    controlled_swap_registers(circ, control, work1, work2);

    with_bit_length_callsite(
        "p.bitlen.rs.swap-old-r.post.deposit",
        "p.bitlen.rs.swap-old-r.post.erase",
        || {
            controlled_xor_rotated_suffix_bit_length(
                circ,
                control,
                l_t,
                work1,
                &old_r_length,
                scratch,
                &[],
            );
        },
    );
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

fn multi_controlled_x_vchain_borrowed(
    circ: &mut Circuit,
    controls: &[&QReg],
    target: &QReg,
    ancillas: &[&QReg],
) {
    match controls.len() {
        0 => circ.x(target),
        1 => circ.cx(controls[0], target),
        2 => circ.ccx(controls[0], controls[1], target),
        count => {
            assert!(ancillas.len() >= count - 2);
            circ.ccx(controls[0], controls[1], ancillas[0]);
            for index in 2..count - 1 {
                circ.ccx(controls[index], ancillas[index - 2], ancillas[index - 1]);
            }
            circ.ccx(controls[count - 1], ancillas[count - 3], target);
            for index in (2..count - 1).rev() {
                circ.ccx(controls[index], ancillas[index - 2], ancillas[index - 1]);
            }
            circ.ccx(controls[0], controls[1], ancillas[0]);
        }
    }
}

fn support_control_uses_preserved_dy_top(
    support_control: &QReg,
    preserved_dy_top_scratch: &[&QReg],
) -> bool {
    assert!(preserved_dy_top_scratch.len() <= 1);
    preserved_dy_top_scratch
        .first()
        .map(|lender| lender.id() == support_control.id())
        .unwrap_or(false)
}

/// Toggle the support-qualified swap predicate into `control` while restoring
/// `l_q` and every workspace lane. The zero predicate alone does not imply the
/// packed-length invariant needed to clear an old-r lender. Materializing the
/// discrepancy in the promised-zero `l_q` makes that invariant explicit:
///
/// `l_q = l_r_prime XOR suffix_bitlen(work2, l_t_prime)`.
///
/// `zero_s AND zero(l_q)` is then exactly the reversible support predicate.
#[allow(clippy::too_many_arguments)]
fn materialize_promised_l_q_swap_discrepancy(
    circ: &mut Circuit,
    zero_q: &QReg,
    zero_s: &QReg,
    control: &QReg,
    l_s: &[QReg],
    work2: &[QReg],
    l_t_prime: &[QReg],
    l_q: &[QReg],
    l_r_prime: &[QReg],
    scratch: &[QReg],
    preserved_dy_top_scratch: &[&QReg],
    predicate_lending: bool,
    support_workspace: &[&QReg],
) {
    assert!(!l_q.is_empty());
    assert_eq!(l_q.len() + 1, l_t_prime.len());
    assert_eq!(l_q.len(), l_r_prime.len());
    assert!(!support_workspace.is_empty());
    let base_control = support_workspace[0];
    assert_ne!(control.id(), base_control.id());

    circ.ccx(zero_q, zero_s, base_control);
    let conditional_zero_q_lender =
        support_control_uses_preserved_dy_top(base_control, preserved_dy_top_scratch);
    if conditional_zero_q_lender {
        // On the active branch base_control implies zero_q=1. Make zero_q a
        // clean hosted high output there. The controlled route preserves its
        // arbitrary inactive-branch value without measuring or resetting it.
        circ.cx(base_control, zero_q);
    }
    with_bit_length_callsite(
        "p.bitlen.rs.swap-support.compute.deposit",
        "p.bitlen.rs.swap-support.compute.erase",
        || {
            if predicate_lending {
                assert_eq!(preserved_dy_top_scratch.len(), 1);
                let mut extended_l_q: Vec<QReg> =
                    l_q.iter().map(QReg::borrowed_alias).collect();
                let (bit_length_scratch, additional_prefix_scratch) =
                    if conditional_zero_q_lender {
                        extended_l_q.push(zero_q.borrowed_alias());
                        (scratch, &[][..])
                    } else {
                        let (bit_length_scratch, high_lane) =
                            scratch.split_at(scratch.len() - 1);
                        extended_l_q.push(high_lane[0].borrowed_alias());
                        (bit_length_scratch, preserved_dy_top_scratch)
                    };
                controlled_xor_rotated_suffix_with_zero_s_lender(
                    circ,
                    zero_s,
                    control,
                    l_s,
                    base_control,
                    l_t_prime,
                    work2,
                    &extended_l_q,
                    bit_length_scratch,
                    additional_prefix_scratch,
                );
            } else {
                panic!("eight-bit l_q requires predicate and preserved-prefix lenders");
            }
        },
    );
    if conditional_zero_q_lender {
        circ.cx(base_control, zero_q);
    }
    for (source, destination) in l_r_prime.iter().zip(l_q) {
        circ.ccx(base_control, source, destination);
    }
    circ.ccx(zero_q, zero_s, base_control);
}

fn toggle_materialized_promised_l_q_swap_control(
    circ: &mut Circuit,
    zero_s: &QReg,
    control: &QReg,
    l_q: &[QReg],
    scratch: &[QReg],
    support_workspace: &[&QReg],
) {
    assert!(!l_q.is_empty());
    assert!(!support_workspace.is_empty());

    for lane in l_q {
        circ.x(lane);
    }
    let controls: Vec<&QReg> = std::iter::once(zero_s).chain(l_q).collect();
    let required_ancillas = controls.len().saturating_sub(2);
    let support_ancillas: Vec<&QReg> = scratch
        .iter()
        .chain(support_workspace.iter().copied())
        .take(required_ancillas)
        .collect();
    assert_eq!(support_ancillas.len(), required_ancillas);
    multi_controlled_x_vchain_borrowed(circ, &controls, control, &support_ancillas);
    for lane in l_q {
        circ.x(lane);
    }
}

#[allow(clippy::too_many_arguments)]
fn unmaterialize_promised_l_q_swap_discrepancy(
    circ: &mut Circuit,
    zero_q: &QReg,
    zero_s: &QReg,
    control: &QReg,
    l_s: &[QReg],
    work2: &[QReg],
    l_t_prime: &[QReg],
    l_q: &[QReg],
    l_r_prime: &[QReg],
    scratch: &[QReg],
    preserved_dy_top_scratch: &[&QReg],
    predicate_lending: bool,
    support_workspace: &[&QReg],
) {
    assert!(!l_q.is_empty());
    assert_eq!(l_q.len() + 1, l_t_prime.len());
    assert_eq!(l_q.len(), l_r_prime.len());
    assert!(!support_workspace.is_empty());
    let base_control = support_workspace[0];
    assert_ne!(control.id(), base_control.id());

    circ.ccx(zero_q, zero_s, base_control);
    for (source, destination) in l_r_prime.iter().zip(l_q) {
        circ.ccx(base_control, source, destination);
    }
    let conditional_zero_q_lender =
        support_control_uses_preserved_dy_top(base_control, preserved_dy_top_scratch);
    if conditional_zero_q_lender {
        circ.cx(base_control, zero_q);
    }
    with_bit_length_callsite(
        "p.bitlen.rs.swap-support.uncompute.deposit",
        "p.bitlen.rs.swap-support.uncompute.erase",
        || {
            if predicate_lending {
                assert_eq!(preserved_dy_top_scratch.len(), 1);
                let mut extended_l_q: Vec<QReg> =
                    l_q.iter().map(QReg::borrowed_alias).collect();
                let (bit_length_scratch, additional_prefix_scratch) =
                    if conditional_zero_q_lender {
                        extended_l_q.push(zero_q.borrowed_alias());
                        (scratch, &[][..])
                    } else {
                        let (bit_length_scratch, high_lane) =
                            scratch.split_at(scratch.len() - 1);
                        extended_l_q.push(high_lane[0].borrowed_alias());
                        (bit_length_scratch, preserved_dy_top_scratch)
                    };
                controlled_xor_rotated_suffix_with_zero_s_lender(
                    circ,
                    zero_s,
                    control,
                    l_s,
                    base_control,
                    l_t_prime,
                    work2,
                    &extended_l_q,
                    bit_length_scratch,
                    additional_prefix_scratch,
                );
            } else {
                panic!("eight-bit l_q requires predicate and preserved-prefix lenders");
            }
        },
    );
    if conditional_zero_q_lender {
        circ.cx(base_control, zero_q);
    }
    circ.ccx(zero_q, zero_s, base_control);
}

#[allow(clippy::too_many_arguments)]
fn toggle_promised_l_q_swap_control(
    circ: &mut Circuit,
    zero_q: &QReg,
    zero_s: &QReg,
    control: &QReg,
    l_s: &[QReg],
    work2: &[QReg],
    l_t_prime: &[QReg],
    l_q: &[QReg],
    l_r_prime: &[QReg],
    scratch: &[QReg],
    preserved_dy_top_scratch: &[&QReg],
    support_workspace: &[&QReg],
) {
    materialize_promised_l_q_swap_discrepancy(
        circ,
        zero_q,
        zero_s,
        control,
        l_s,
        work2,
        l_t_prime,
        l_q,
        l_r_prime,
        scratch,
        preserved_dy_top_scratch,
        false,
        support_workspace,
    );
    toggle_materialized_promised_l_q_swap_control(
        circ,
        zero_s,
        control,
        l_q,
        scratch,
        support_workspace,
    );
    unmaterialize_promised_l_q_swap_discrepancy(
        circ,
        zero_q,
        zero_s,
        control,
        l_s,
        work2,
        l_t_prime,
        l_q,
        l_r_prime,
        scratch,
        preserved_dy_top_scratch,
        false,
        support_workspace,
    );
}

/// Specialized swap used only while the support-qualified zero control is
/// materialized. When `control` is one, `l_q` is a clean lender and the old
/// `l_r_prime` agrees with the swapped-in packed suffix; when it is zero, every
/// use of the lender is controlled and arbitrary `l_q` data is preserved.
fn controlled_xor_rotated_suffix_lq8_hosted(
    circ: &mut Circuit,
    control: &QReg,
    left_length: &[QReg],
    work: &[QReg],
    l_q: &[QReg],
    scratch: &[QReg],
    hosted_high_output: Option<&QReg>,
    preserved_dy_top_scratch: &[&QReg],
    predicate_prefix_scratch: &[&QReg],
) {
    assert_eq!(left_length.len(), l_q.len() + 1);
    assert_eq!(preserved_dy_top_scratch.len(), 1);
    assert_eq!(
        predicate_prefix_scratch.len(),
        2 - usize::from(hosted_high_output.is_some())
    );
    let mut extended_l_q: Vec<QReg> = l_q.iter().map(QReg::borrowed_alias).collect();
    let bit_length_scratch = if let Some(hosted_high_output) = hosted_high_output {
        extended_l_q.push(hosted_high_output.borrowed_alias());
        scratch
    } else {
        let (bit_length_scratch, high_lane) = scratch.split_at(scratch.len() - 1);
        extended_l_q.push(high_lane[0].borrowed_alias());
        bit_length_scratch
    };
    let borrowed_prefix_scratch: Vec<&QReg> = predicate_prefix_scratch
        .iter()
        .copied()
        .chain(preserved_dy_top_scratch.iter().copied())
        .collect();
    controlled_xor_rotated_suffix_bit_length_with_prefix_scratch(
        circ,
        control,
        left_length,
        work,
        &extended_l_q,
        bit_length_scratch,
        &borrowed_prefix_scratch,
    );
}

#[allow(clippy::too_many_arguments)]
fn conditional_work_and_length_swap_promised_l_q_zero(
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
    scratch: &[QReg],
    hosted_high_output: Option<&QReg>,
    preserved_dy_top_scratch: &[&QReg],
    predicate_prefix_scratch: &[&QReg],
) {
    assert_eq!(work1.len(), work2.len());
    assert_eq!(l_t.len(), l_t_prime.len());
    assert_eq!(l_t.len(), l_q.len() + 1);
    assert_eq!(l_q.len(), l_r_prime.len());
    assert_l_r_prime_metadata_width(l_t.len(), l_r_prime.len());

    with_bit_length_callsite(
        "p.bitlen.rs.promised-swap.pre.deposit",
        "p.bitlen.rs.promised-swap.pre.erase",
        || {
            controlled_xor_rotated_suffix_lq8_hosted(
                circ,
                control,
                l_t,
                work1,
                l_q,
                scratch,
                hosted_high_output,
                preserved_dy_top_scratch,
                predicate_prefix_scratch,
            );
        },
    );

    for (current, next) in l_t.iter().zip(l_t_prime) {
        circ.cswap(control, current, next);
    }
    for (current, next) in l_r_prime.iter().zip(l_q) {
        circ.cswap(control, current, next);
    }
    controlled_swap_registers(circ, control, work1, work2);

    with_bit_length_callsite(
        "p.bitlen.rs.promised-swap.post.deposit",
        "p.bitlen.rs.promised-swap.post.erase",
        || {
            controlled_xor_rotated_suffix_lq8_hosted(
                circ,
                control,
                l_t,
                work1,
                l_q,
                scratch,
                hosted_high_output,
                preserved_dy_top_scratch,
                predicate_prefix_scratch,
            );
        },
    );
    circ.cx(control, iteration_parity);

    let _ = l_s;
}

#[allow(clippy::too_many_arguments)]
fn conditional_work_and_length_swap_promised_l_q_zero_inverse(
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
    hosted_high_output: Option<&QReg>,
    preserved_dy_top_scratch: &[&QReg],
    predicate_prefix_scratch: &[&QReg],
) {
    conditional_work_and_length_swap_promised_l_q_zero(
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
        hosted_high_output,
        preserved_dy_top_scratch,
        predicate_prefix_scratch,
    );
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum PromisedLqSwapRoute {
    Configured,
    Allocated,
    PromisedAllocated,
    BorrowLq,
}

/// Exchange the reachable packed Work1/Work2 metadata without materializing
/// `l_t_prime`. On the active branch, `l_t` follows
/// `A -> A xor B -> B`; `l_q` similarly carries the old Work1 suffix length
/// across the remainder-length swap and is erased from the post-swap word.
#[allow(clippy::too_many_arguments)]
fn conditional_work_and_length_swap_direct_metadata(
    circ: &mut Circuit,
    zero_q: &QReg,
    zero_s: &QReg,
    control: &QReg,
    iteration_parity: &QReg,
    work1: &[QReg],
    work2: &[QReg],
    l_t: &[QReg],
    l_q: &[QReg],
    l_s: &[QReg],
    l_r_prime: &[QReg],
    scratch: &[QReg],
    preserved_dy_top_scratch: &[&QReg],
) {
    assert_eq!(work1.len(), work2.len());
    assert_eq!(l_t.len(), l_q.len() + 1);
    assert_eq!(l_q.len(), l_r_prime.len());
    assert_l_r_prime_metadata_width(l_t.len(), l_r_prime.len());
    assert_eq!(preserved_dy_top_scratch.len(), 1);
    assert!(scratch.len() >= l_t.len().saturating_sub(2));

    // Reachable scheduled states satisfy the packed-length invariant whenever
    // q=s=0. Materialize that base control, then release both zero-predicate
    // targets as clean decoder lenders while retaining the control.
    circ.ccx(zero_q, zero_s, control);
    uncompute_zero(circ, l_s, zero_s, scratch);
    uncompute_zero(circ, l_q, zero_q, scratch);

    let predicate_prefix_scratch = [zero_s];
    with_bit_length_callsite(
        "p.bitlen.rs.direct-swap-old-r.pre.deposit",
        "p.bitlen.rs.direct-swap-old-r.pre.erase",
        || {
            controlled_xor_rotated_suffix_lq8_hosted(
                circ,
                control,
                l_t,
                work1,
                l_q,
                scratch,
                Some(zero_q),
                preserved_dy_top_scratch,
                &predicate_prefix_scratch,
            );
        },
    );
    with_bit_length_callsite(
        "p.bitlen.rs.direct-swap-new-t.pre.deposit",
        "p.bitlen.rs.direct-swap-new-t.pre.erase",
        || {
            controlled_xor_rotated_prefix_bit_length(
                circ,
                control,
                l_r_prime,
                work2,
                l_t,
                scratch,
                &[zero_q, zero_s],
            );
        },
    );

    for (current, next) in l_r_prime.iter().zip(l_q) {
        circ.cswap(control, current, next);
    }
    controlled_swap_registers(circ, control, work1, work2);

    with_bit_length_callsite(
        "p.bitlen.rs.direct-swap-old-t.post.deposit",
        "p.bitlen.rs.direct-swap-old-t.post.erase",
        || {
            controlled_xor_rotated_prefix_bit_length(
                circ,
                control,
                l_r_prime,
                work2,
                l_t,
                scratch,
                &[zero_q, zero_s],
            );
        },
    );
    with_bit_length_callsite(
        "p.bitlen.rs.direct-swap-old-r-prime.post.deposit",
        "p.bitlen.rs.direct-swap-old-r-prime.post.erase",
        || {
            controlled_xor_rotated_suffix_lq8_hosted(
                circ,
                control,
                l_t,
                work1,
                l_q,
                scratch,
                Some(zero_q),
                preserved_dy_top_scratch,
                &predicate_prefix_scratch,
            );
        },
    );
    circ.cx(control, iteration_parity);

    compute_zero(circ, l_q, zero_q, scratch);
    compute_zero(circ, l_s, zero_s, scratch);
    circ.ccx(zero_q, zero_s, control);
}

/// Keep the zero predicates live across the swap. This is the only production
/// entry point that can select the promised `l_q` lender route.
#[allow(clippy::too_many_arguments)]
fn conditional_work_and_length_swap_under_zero_predicate(
    circ: &mut Circuit,
    zero_q: &QReg,
    zero_s: &QReg,
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
    preserved_dy_top_scratch: &[&QReg],
    inverse: bool,
    route: PromisedLqSwapRoute,
) {
    let route = match route {
        PromisedLqSwapRoute::Configured if promised_l_q_swap_borrow_requested() => {
            PromisedLqSwapRoute::BorrowLq
        }
        PromisedLqSwapRoute::Configured => PromisedLqSwapRoute::Allocated,
        route => route,
    };
    assert!(
        q845_swap_only_swap_dependencies_satisfied(),
        "Q845 swap-only t-prime lifecycle requires the promised l_q swap route"
    );
    if q845_swap_only_t_prime_length_requested() {
        assert_eq!(
            route,
            PromisedLqSwapRoute::BorrowLq,
            "Q845 swap-only t-prime lifecycle cannot use an allocated swap route"
        );
    }
    if q830_direct_swap_metadata_requested() {
        assert!(
            q845_swap_only_t_prime_length_requested(),
            "direct swap metadata requires the swap-only t-prime lifecycle"
        );
        assert_eq!(
            route,
            PromisedLqSwapRoute::BorrowLq,
            "direct swap metadata requires promised l_q lending"
        );
        conditional_work_and_length_swap_direct_metadata(
            circ,
            zero_q,
            zero_s,
            control,
            iteration_parity,
            work1,
            work2,
            l_t,
            l_q,
            l_s,
            l_r_prime,
            scratch,
            preserved_dy_top_scratch,
        );
        return;
    }
    if route == PromisedLqSwapRoute::Allocated {
        circ.ccx(zero_q, zero_s, control);
        if inverse {
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
                t_window,
                r_window,
                scratch,
            );
        } else {
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
        circ.ccx(zero_q, zero_s, control);
        return;
    }

    let required_support_ancillas = l_q.len().saturating_sub(1);
    let support_workspace_width =
        1usize.max(required_support_ancillas.saturating_sub(scratch.len()));
    let borrow_support_workspace = q839_seven_plateau_lenders_requested()
        && !preserved_dy_top_scratch.is_empty();
    if borrow_support_workspace {
        assert_eq!(support_workspace_width, 1);
        assert_eq!(preserved_dy_top_scratch.len(), 1);
        assert!(
            promised_swap_support_lifetime_fusion_requested(),
            "q837 support lending requires fused promised-swap support lifetime"
        );
        assert!(
            sub800_raw_prefix_predicate_lender_requested(),
            "q837 support lending requires predicate-prefix lending"
        );
        let lender = preserved_dy_top_scratch[0];
        assert_ne!(lender.id(), zero_q.id());
        assert_ne!(lender.id(), zero_s.id());
        assert_ne!(lender.id(), control.id());
        assert!(scratch.iter().all(|lane| lane.id() != lender.id()));
        Q839_SEVEN_PLATEAU_SUPPORT_LOANS.fetch_add(1, Ordering::Relaxed);
    } else if q839_seven_plateau_lenders_requested() {
        Q839_SEVEN_PLATEAU_SUPPORT_FALLBACKS.fetch_add(1, Ordering::Relaxed);
    }
    let support_lender_entry_ops_idx =
        borrow_support_workspace.then(|| circ.total_ops() as usize);
    let support_workspace_owned = (!borrow_support_workspace).then(|| {
        circ.alloc_qreg_bits("rs.swap-length.promised-support", support_workspace_width)
    });
    let support_workspace: Vec<&QReg> = if borrow_support_workspace {
        vec![preserved_dy_top_scratch[0]]
    } else {
        support_workspace_owned
            .as_ref()
            .expect("owned support workspace")
            .iter()
            .collect()
    };
    let raw_prefix_scratch = if sub800_raw_prefix_preserved_lender_requested() {
        preserved_dy_top_scratch
    } else {
        &[]
    };
    if q845_swap_only_t_prime_length_requested() {
        // `l_t_prime` is zero between swap sites. The base zero predicate is
        // sufficient to materialize the packed Work2 coefficient length used
        // by the support test; the same operation erases the swapped-in length
        // after the controlled exchange.
        if sub800_raw_prefix_predicate_lender_requested() {
            controlled_xor_rotated_prefix_with_predicate_lenders(
                circ,
                zero_q,
                zero_s,
                control,
                l_q,
                l_s,
                l_r_prime,
                work2,
                l_t_prime,
                scratch,
            );
        } else {
            circ.ccx(zero_q, zero_s, control);
            controlled_xor_rotated_prefix_bit_length(
                circ,
                control,
                l_r_prime,
                work2,
                l_t_prime,
                scratch,
                raw_prefix_scratch,
            );
            circ.ccx(zero_q, zero_s, control);
        }
    }
    let fuse_support_lifetime = promised_swap_support_lifetime_fusion_requested();
    if fuse_support_lifetime {
        assert_eq!(
            route,
            PromisedLqSwapRoute::BorrowLq,
            "support-lifetime fusion requires the promised l_q lender route"
        );
        // Keep D = suffix_bitlen(work2, l_t_prime) XOR l_r_prime in l_q
        // across the controlled swap. D is zero exactly on the control-on
        // branch, while every swap operation preserves arbitrary D when the
        // control is off. This removes one complete support scan pair.
        materialize_promised_l_q_swap_discrepancy(
            circ,
            zero_q,
            zero_s,
            control,
            l_s,
            work2,
            l_t_prime,
            l_q,
            l_r_prime,
            scratch,
            preserved_dy_top_scratch,
            sub800_raw_prefix_predicate_lender_requested(),
            &support_workspace,
        );
        toggle_materialized_promised_l_q_swap_control(
            circ,
            zero_s,
            control,
            l_q,
            scratch,
            &support_workspace,
        );
    } else {
        toggle_promised_l_q_swap_control(
            circ,
            zero_q,
            zero_s,
            control,
            l_s,
            work2,
            l_t_prime,
            l_q,
            l_r_prime,
            scratch,
            preserved_dy_top_scratch,
            &support_workspace,
        );
    }
    let conditional_zero_q_lender =
        borrow_support_workspace && sub800_raw_prefix_predicate_lender_requested();
    let hosted_swap_high_output = conditional_zero_q_lender.then_some(zero_q);
    let predicate_swap_lenders = if sub800_raw_prefix_predicate_lender_requested() {
        assert_eq!(
            route,
            PromisedLqSwapRoute::BorrowLq,
            "predicate lending requires the promised l_q swap route"
        );
        uncompute_zero(circ, l_s, zero_s, scratch);
        if conditional_zero_q_lender {
            // The support predicate implies zero_q=1. Its controlled toggle
            // makes it a clean hosted high output on the active branch. The
            // inactive branch preserves its arbitrary value without lending it
            // to measurement-based prefix scratch.
            circ.cx(control, zero_q);
            vec![zero_s]
        } else {
            vec![zero_s, support_workspace[0]]
        }
    } else {
        Vec::new()
    };
    match (inverse, route) {
        (false, PromisedLqSwapRoute::PromisedAllocated) => conditional_work_and_length_swap(
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
        ),
        (true, PromisedLqSwapRoute::PromisedAllocated) => conditional_work_and_length_swap_inverse(
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
        ),
        (false, PromisedLqSwapRoute::BorrowLq) => {
            conditional_work_and_length_swap_promised_l_q_zero(
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
                hosted_swap_high_output,
                preserved_dy_top_scratch,
                &predicate_swap_lenders,
            )
        }
        (true, PromisedLqSwapRoute::BorrowLq) => {
            conditional_work_and_length_swap_promised_l_q_zero_inverse(
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
                hosted_swap_high_output,
                preserved_dy_top_scratch,
                &predicate_swap_lenders,
            )
        }
        (_, PromisedLqSwapRoute::Allocated | PromisedLqSwapRoute::Configured) => {
            unreachable!("promised swap route resolved")
        }
    }
    if conditional_zero_q_lender {
        circ.cx(control, zero_q);
    }
    if !predicate_swap_lenders.is_empty() {
        compute_zero(circ, l_s, zero_s, scratch);
    }
    if fuse_support_lifetime {
        // The promised swap restores the retained discrepancy in l_q: zero on
        // the swap branch and the untouched arbitrary value off branch.
        toggle_materialized_promised_l_q_swap_control(
            circ,
            zero_s,
            control,
            l_q,
            scratch,
            &support_workspace,
        );
        unmaterialize_promised_l_q_swap_discrepancy(
            circ,
            zero_q,
            zero_s,
            control,
            l_s,
            work2,
            l_t_prime,
            l_q,
            l_r_prime,
            scratch,
            preserved_dy_top_scratch,
            sub800_raw_prefix_predicate_lender_requested(),
            &support_workspace,
        );
    } else {
        toggle_promised_l_q_swap_control(
            circ,
            zero_q,
            zero_s,
            control,
            l_s,
            work2,
            l_t_prime,
            l_q,
            l_r_prime,
            scratch,
            preserved_dy_top_scratch,
            &support_workspace,
        );
    }
    if q845_swap_only_t_prime_length_requested() {
        if sub800_raw_prefix_predicate_lender_requested() {
            controlled_xor_rotated_prefix_with_predicate_lenders(
                circ,
                zero_q,
                zero_s,
                control,
                l_q,
                l_s,
                l_r_prime,
                work2,
                l_t_prime,
                scratch,
            );
        } else {
            circ.ccx(zero_q, zero_s, control);
            controlled_xor_rotated_prefix_bit_length(
                circ,
                control,
                l_r_prime,
                work2,
                l_t_prime,
                scratch,
                raw_prefix_scratch,
            );
            circ.ccx(zero_q, zero_s, control);
        }
    }
    if let Some(entry_ops_idx) = support_lender_entry_ops_idx {
        record_q839_support_lender_borrow_window(
            preserved_dy_top_scratch[0],
            entry_ops_idx,
            circ.total_ops() as usize,
        );
    }
    if let Some(support_workspace) = support_workspace_owned {
        free_clean(circ, support_workspace);
    }
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
    assert_eq!(l_t.len(), l_q.len() + 1);
    assert_eq!(l_t.len(), l_s.len());
    assert!(window_upper <= total_work_width);
    cuccaro_add_zero_extended_no_overflow(
        circ,
        l_q,
        l_t,
        scratch.length_carry,
        scratch.length_overflow,
    );
    // At global position j=K, sign(l_t) is one iff l_t+l_q <= K-2.
    sub_const_mod_2n(circ, l_t, window_upper - 1, scratch.constant);
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
    add_const_mod_2n(circ, l_t, window_upper - 1, scratch.constant);
    cuccaro_sub_zero_extended_no_overflow(
        circ,
        l_q,
        l_t,
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
    let layout = selected_remainder_scratch_layout();
    if layout == RemainderScratchLayout::PhaseOverlaid {
        assert_remainder_scratch_disjoint_from_data(
            "remainder sub",
            &[phase1, sign],
            &[work1, work2, l_t, l_q, l_s, l_r_prime],
            scratch,
        );
    }
    let scratch = split_remainder_scratch_for_layout(scratch, l_t.len(), l_r_prime.len(), layout);
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
        toggle_remainder_range_active(circ, scratch.operation, l_t, l_s, &scratch);
        circ.ccx(scratch.active, &work2[index], &work1[index]);
        circ.ccx(scratch.active, scratch.carry, &work2[index]);
        multi_controlled_x_vchain(
            circ,
            &[scratch.active, &work2[index], &work1[index]],
            scratch.carry,
            std::slice::from_ref(scratch.tmp),
        );
        toggle_remainder_range_active(circ, scratch.operation, l_t, l_s, &scratch);
        if index != 0 {
            walk_remainder_range_down(circ, l_t, l_s, &scratch);
        }
    }
    circ.ccx(scratch.operation, scratch.carry, sign);
    for index in 0..work1.len() {
        toggle_remainder_range_active(circ, scratch.operation, l_t, l_s, &scratch);
        multi_controlled_x_vchain(
            circ,
            &[scratch.active, &work2[index], &work1[index]],
            scratch.carry,
            std::slice::from_ref(scratch.tmp),
        );
        circ.ccx(scratch.active, scratch.carry, &work2[index]);
        circ.ccx(scratch.active, scratch.carry, &work1[index]);
        toggle_remainder_range_active(circ, scratch.operation, l_t, l_s, &scratch);
        if index + 1 != work1.len() {
            walk_remainder_range_up(circ, l_t, l_s, &scratch);
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
    let layout = selected_remainder_scratch_layout();
    if layout == RemainderScratchLayout::PhaseOverlaid {
        assert_remainder_scratch_disjoint_from_data(
            "remainder sub inverse",
            &[phase1, sign],
            &[work1, work2, l_t, l_q, l_s, l_r_prime],
            scratch,
        );
    }
    let scratch = split_remainder_scratch_for_layout(scratch, l_t.len(), l_r_prime.len(), layout);
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
            walk_remainder_range_down(circ, l_t, l_s, &scratch);
        }
        toggle_remainder_range_active(circ, scratch.operation, l_t, l_s, &scratch);
        circ.ccx(scratch.active, scratch.carry, &work1[index]);
        circ.ccx(scratch.active, scratch.carry, &work2[index]);
        multi_controlled_x_vchain(
            circ,
            &[scratch.active, &work2[index], &work1[index]],
            scratch.carry,
            std::slice::from_ref(scratch.tmp),
        );
        toggle_remainder_range_active(circ, scratch.operation, l_t, l_s, &scratch);
    }
    circ.ccx(scratch.operation, scratch.carry, sign);
    for index in 0..work1.len() {
        if index != 0 {
            walk_remainder_range_up(circ, l_t, l_s, &scratch);
        }
        toggle_remainder_range_active(circ, scratch.operation, l_t, l_s, &scratch);
        multi_controlled_x_vchain(
            circ,
            &[scratch.active, &work2[index], &work1[index]],
            scratch.carry,
            std::slice::from_ref(scratch.tmp),
        );
        circ.ccx(scratch.active, scratch.carry, &work2[index]);
        circ.ccx(scratch.active, &work2[index], &work1[index]);
        toggle_remainder_range_active(circ, scratch.operation, l_t, l_s, &scratch);
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
    let layout = selected_remainder_scratch_layout();
    if layout == RemainderScratchLayout::PhaseOverlaid {
        assert_remainder_scratch_disjoint_from_data(
            "remainder add",
            &[phase1, phase2, sign],
            &[work1, work2, l_t, l_q, l_s, l_r_prime],
            scratch,
        );
    }
    let scratch = split_remainder_scratch_for_layout(scratch, l_t.len(), l_r_prime.len(), layout);
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
        toggle_remainder_range_active(circ, scratch.enable, l_t, l_s, &scratch);
        circ.ccx(scratch.active, scratch.carry, &work1[index]);
        circ.ccx(scratch.active, scratch.carry, &work2[index]);
        multi_controlled_x_vchain(
            circ,
            &[scratch.active, &work2[index], &work1[index]],
            scratch.carry,
            std::slice::from_ref(scratch.tmp),
        );
        toggle_remainder_range_active(circ, scratch.enable, l_t, l_s, &scratch);
        if index != 0 {
            walk_remainder_range_down(circ, l_t, l_s, &scratch);
        }
    }
    for index in 0..work1.len() {
        toggle_remainder_range_active(circ, scratch.enable, l_t, l_s, &scratch);
        multi_controlled_x_vchain(
            circ,
            &[scratch.active, &work2[index], &work1[index]],
            scratch.carry,
            std::slice::from_ref(scratch.tmp),
        );
        circ.ccx(scratch.active, scratch.carry, &work2[index]);
        circ.ccx(scratch.active, &work2[index], &work1[index]);
        toggle_remainder_range_active(circ, scratch.enable, l_t, l_s, &scratch);
        if index + 1 != work1.len() {
            walk_remainder_range_up(circ, l_t, l_s, &scratch);
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
    let layout = selected_remainder_scratch_layout();
    if layout == RemainderScratchLayout::PhaseOverlaid {
        assert_remainder_scratch_disjoint_from_data(
            "remainder add inverse",
            &[phase1, phase2, sign],
            &[work1, work2, l_t, l_q, l_s, l_r_prime],
            scratch,
        );
    }
    let scratch = split_remainder_scratch_for_layout(scratch, l_t.len(), l_r_prime.len(), layout);
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
            walk_remainder_range_down(circ, l_t, l_s, &scratch);
        }
        toggle_remainder_range_active(circ, scratch.enable, l_t, l_s, &scratch);
        circ.ccx(scratch.active, &work2[index], &work1[index]);
        circ.ccx(scratch.active, scratch.carry, &work2[index]);
        multi_controlled_x_vchain(
            circ,
            &[scratch.active, &work2[index], &work1[index]],
            scratch.carry,
            std::slice::from_ref(scratch.tmp),
        );
        toggle_remainder_range_active(circ, scratch.enable, l_t, l_s, &scratch);
    }
    for index in 0..work1.len() {
        if index != 0 {
            walk_remainder_range_up(circ, l_t, l_s, &scratch);
        }
        toggle_remainder_range_active(circ, scratch.enable, l_t, l_s, &scratch);
        multi_controlled_x_vchain(
            circ,
            &[scratch.active, &work2[index], &work1[index]],
            scratch.carry,
            std::slice::from_ref(scratch.tmp),
        );
        circ.ccx(scratch.active, scratch.carry, &work2[index]);
        circ.ccx(scratch.active, scratch.carry, &work1[index]);
        toggle_remainder_range_active(circ, scratch.enable, l_t, l_s, &scratch);
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

    with_bit_length_callsite(
        "p.bitlen.rs.coeff-target.deposit",
        "p.bitlen.rs.coeff-target.erase",
        || {
            controlled_xor_rotated_prefix_bit_length(
                circ,
                control,
                &boundary,
                work2,
                target_length,
                &[],
                &[],
            );
        },
    );

    cuccaro_sub_mod_2n(circ, l_s, &boundary, &boundary_carry, &boundary_overflow);
    circ.zero_and_free(boundary_overflow);
    circ.zero_and_free(boundary_carry);
    for (source, destination) in l_r_prime.iter().zip(&boundary) {
        circ.cx(source, destination);
    }
    free_clean(circ, boundary);
}

fn borrowed_coefficient_comparator_requested() -> bool {
    std::env::var("LOWQ_REUSE_COEFFICIENT_COMPARATOR_SCRATCH")
        .ok()
        .as_deref()
        == Some("1")
}

fn toggle_coefficient_length_above_boundary_allocated(
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

fn assert_borrowed_coefficient_comparator_preconditions(
    target_length: &[QReg],
    l_t: &[QReg],
    l_s: &[QReg],
    target: &QReg,
    scratch: &[&QReg],
) {
    assert!(!target_length.is_empty());
    assert_eq!(target_length.len(), l_t.len());
    assert_eq!(target_length.len(), l_s.len());
    assert_eq!(
        scratch.len(),
        3,
        "borrow-only coefficient comparator requires exactly three clean caller lanes"
    );

    let mut ids = Vec::with_capacity(3 * target_length.len() + scratch.len() + 1);
    for (index, lane) in target_length
        .iter()
        .chain(l_t)
        .chain(l_s)
        .chain(std::iter::once(target))
        .chain(scratch.iter().copied())
        .enumerate()
    {
        assert!(
            !ids.contains(&lane.id()),
            "borrow-only coefficient comparator lane {index} aliases an operand or scratch lane"
        );
        ids.push(lane.id());
    }
}

/// Toggle `target_length > (l_t + l_s mod 2^n)` into `target`.
///
/// `scratch` is caller-owned clean storage in the order `(zero pad, carry,
/// borrow)`. The add/subtract pairs restore all three lanes exactly. Cleanliness
/// is a caller precondition; the exhaustive proof below checks it on every
/// reduced-width basis state.
fn toggle_coefficient_length_above_boundary_borrowed(
    circ: &mut Circuit,
    target_length: &[QReg],
    l_t: &[QReg],
    l_s: &[QReg],
    target: &QReg,
    scratch: &[&QReg],
) {
    assert_borrowed_coefficient_comparator_preconditions(target_length, l_t, l_s, target, scratch);
    let length_width = l_t.len();
    let wide = length_width + 1;
    let zero_pad = scratch[0];
    let carry = scratch[1];
    let borrow = scratch[2];

    let boundary = circ.alloc_qreg_bits("rs.coeff-compare.borrow-boundary", wide);
    for (source, destination) in l_t.iter().zip(&boundary) {
        circ.cx(source, destination);
    }
    cuccaro_add_mod_2n_no_overflow(circ, l_s, &boundary[..length_width], carry);

    let mut zero_extended_target: Vec<&QReg> = target_length.iter().collect();
    zero_extended_target.push(zero_pad);
    cuccaro_sub_mod_2n_refs(circ, &zero_extended_target, &boundary, carry, borrow);
    circ.cx(borrow, target);
    cuccaro_add_mod_2n_refs(circ, &zero_extended_target, &boundary, carry, borrow);

    cuccaro_sub_mod_2n_no_overflow(circ, l_s, &boundary[..length_width], carry);
    for (source, destination) in l_t.iter().zip(&boundary) {
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
    borrowed_scratch: &[&QReg],
) {
    if borrowed_coefficient_comparator_requested() {
        toggle_coefficient_length_above_boundary_borrowed(
            circ,
            target_length,
            l_t,
            l_s,
            target,
            borrowed_scratch,
        );
    } else {
        toggle_coefficient_length_above_boundary_allocated(circ, target_length, l_t, l_s, target);
    }
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
struct CoefficientLessThanLifetimeTrace {
    after_first_boundary: usize,
    before_second_boundary: usize,
}

thread_local! {
    static COEFFICIENT_LESS_THAN_LIFETIME_TRACE: std::cell::Cell<(
        bool,
        CoefficientLessThanLifetimeTrace,
    )> = std::cell::Cell::new((false, CoefficientLessThanLifetimeTrace {
        after_first_boundary: 0,
        before_second_boundary: 0,
    }));
}

fn begin_coefficient_less_than_lifetime_trace() {
    COEFFICIENT_LESS_THAN_LIFETIME_TRACE.with(|trace| {
        trace.set((true, CoefficientLessThanLifetimeTrace::default()));
    });
}

fn finish_coefficient_less_than_lifetime_trace() -> CoefficientLessThanLifetimeTrace {
    COEFFICIENT_LESS_THAN_LIFETIME_TRACE.with(|trace| {
        let (_, snapshot) = trace.get();
        trace.set((false, snapshot));
        snapshot
    })
}

fn record_coefficient_less_than_lifetime_boundary(circ: &Circuit, first: bool) {
    COEFFICIENT_LESS_THAN_LIFETIME_TRACE.with(|trace| {
        let (enabled, mut snapshot) = trace.get();
        if enabled {
            if first {
                snapshot.after_first_boundary = circ.total_ops() as usize;
            } else {
                snapshot.before_second_boundary = circ.total_ops() as usize;
            }
            trace.set((enabled, snapshot));
        }
    });
}

#[allow(clippy::too_many_arguments)]
fn assert_coefficient_less_than_lane_reuse_preconditions(
    control: &QReg,
    work1: &[QReg],
    work2: &[QReg],
    l_t: &[QReg],
    l_s: &[QReg],
    target_length: &[QReg],
    target: &QReg,
    scratch: &[&QReg],
) {
    assert_eq!(work1.len(), work2.len());
    assert_eq!(l_t.len(), l_s.len());
    assert_eq!(l_t.len(), target_length.len());
    assert_eq!(
        scratch.len(),
        3,
        "coefficient less-than lane reuse requires carry, active, and tmp lenders"
    );
    for (index, lane) in scratch.iter().enumerate() {
        assert!(
            scratch[..index].iter().all(|other| other.id() != lane.id()),
            "coefficient less-than lender {index} aliases an earlier lender"
        );
        assert_ne!(lane.id(), control.id());
        assert_ne!(lane.id(), target.id());
        assert!(work1.iter().all(|other| other.id() != lane.id()));
        assert!(work2.iter().all(|other| other.id() != lane.id()));
        assert!(l_t.iter().all(|other| other.id() != lane.id()));
        assert!(l_s.iter().all(|other| other.id() != lane.id()));
        assert!(target_length.iter().all(|other| other.id() != lane.id()));
    }
}

#[allow(clippy::too_many_arguments)]
fn toggle_coefficient_less_than_with_lane_route(
    circ: &mut Circuit,
    control: &QReg,
    work1: &[QReg],
    work2: &[QReg],
    l_t: &[QReg],
    l_s: &[QReg],
    target_length: &[QReg],
    target: &QReg,
    compare_scratch: &[&QReg],
    reuse_caller_lanes: bool,
) {
    if reuse_caller_lanes {
        assert_coefficient_less_than_lane_reuse_preconditions(
            control,
            work1,
            work2,
            l_t,
            l_s,
            target_length,
            target,
            compare_scratch,
        );
    }
    let above_t = circ.alloc_qreg("rs.coeff-compare.above-t");
    toggle_coefficient_length_above_boundary(
        circ,
        target_length,
        l_t,
        l_s,
        &above_t,
        compare_scratch,
    );
    record_coefficient_less_than_lifetime_boundary(circ, true);

    let cursor_scratch = circ.alloc_qreg_bits(
        "rs.coeff-compare.cursor-scratch",
        l_t.len().saturating_sub(1),
    );
    // Borrow l_t as the cursor; the reverse scan and final increment restore it.
    decrement_mod_2n(circ, l_t, &cursor_scratch);
    let carry_owned =
        (!reuse_caller_lanes).then(|| circ.alloc_qreg("rs.coeff-compare.local-carry"));
    let active_owned = (!reuse_caller_lanes).then(|| circ.alloc_qreg("rs.coeff-compare.active"));
    let tmp_owned = (!reuse_caller_lanes).then(|| circ.alloc_qreg("rs.coeff-compare.tmp"));
    let carry = if reuse_caller_lanes {
        compare_scratch[0]
    } else {
        carry_owned.as_ref().expect("owned coefficient carry")
    };
    let active = if reuse_caller_lanes {
        compare_scratch[1]
    } else {
        active_owned.as_ref().expect("owned coefficient active")
    };
    let tmp = if reuse_caller_lanes {
        compare_scratch[2]
    } else {
        tmp_owned.as_ref().expect("owned coefficient tmp")
    };
    let bracket = production_coefficient_nonnegative_bracket();

    for index in 0..work1.len() {
        begin_coefficient_nonnegative_bracket(circ, control, l_t, active, bracket);
        circ.ccx(&active, &work1[index], &work2[index]);
        circ.ccx(&active, &carry, &work1[index]);
        multi_controlled_x_vchain(
            circ,
            &[&active, &work1[index], &work2[index]],
            &carry,
            std::slice::from_ref(&tmp),
        );
        end_coefficient_nonnegative_bracket(circ, control, l_t, active, bracket);
        if index + 1 != work1.len() {
            decrement_mod_2n(circ, l_t, &cursor_scratch);
        }
    }

    circ.x(&above_t);
    circ.ccx(&carry, &above_t, target);
    circ.x(&above_t);

    for index in (0..work1.len()).rev() {
        if index + 1 != work1.len() {
            increment_mod_2n(circ, l_t, &cursor_scratch);
        }
        begin_coefficient_nonnegative_bracket(circ, control, l_t, active, bracket);
        multi_controlled_x_vchain(
            circ,
            &[&active, &work1[index], &work2[index]],
            &carry,
            std::slice::from_ref(&tmp),
        );
        circ.ccx(&active, &carry, &work1[index]);
        circ.ccx(&active, &work1[index], &work2[index]);
        end_coefficient_nonnegative_bracket(circ, control, l_t, active, bracket);
    }
    increment_mod_2n(circ, l_t, &cursor_scratch);
    if let Some(tmp) = tmp_owned {
        circ.zero_and_free(tmp);
    }
    if let Some(active) = active_owned {
        circ.zero_and_free(active);
    }
    if let Some(carry) = carry_owned {
        circ.zero_and_free(carry);
    }
    free_clean(circ, cursor_scratch);

    record_coefficient_less_than_lifetime_boundary(circ, false);
    toggle_coefficient_length_above_boundary(
        circ,
        target_length,
        l_t,
        l_s,
        &above_t,
        compare_scratch,
    );
    circ.zero_and_free(above_t);
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
    compare_scratch: &[&QReg],
) {
    toggle_coefficient_less_than_with_lane_route(
        circ,
        control,
        work1,
        work2,
        l_t,
        l_s,
        target_length,
        target,
        compare_scratch,
        coefficient_less_than_lane_reuse_requested(),
    );
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
struct CoefficientAddLenderTrace {
    calls: usize,
    entry_ops_idx: usize,
    restore_ops_idx: usize,
    lender_mask: u64,
}

thread_local! {
    static COEFFICIENT_ADD_LENDER_TRACE: std::cell::Cell<(
        bool,
        CoefficientAddLenderTrace,
    )> = std::cell::Cell::new((false, CoefficientAddLenderTrace {
        calls: 0,
        entry_ops_idx: 0,
        restore_ops_idx: 0,
        lender_mask: 0,
    }));
}

fn begin_coefficient_add_lender_trace() {
    COEFFICIENT_ADD_LENDER_TRACE.with(|trace| {
        trace.set((true, CoefficientAddLenderTrace::default()));
    });
}

fn finish_coefficient_add_lender_trace() -> CoefficientAddLenderTrace {
    COEFFICIENT_ADD_LENDER_TRACE.with(|trace| {
        let (_, snapshot) = trace.get();
        trace.set((false, snapshot));
        snapshot
    })
}

fn record_coefficient_add_lender_entry(circ: &Circuit, lenders: &[&QReg]) {
    COEFFICIENT_ADD_LENDER_TRACE.with(|trace| {
        let (enabled, mut snapshot) = trace.get();
        if enabled {
            assert_eq!(snapshot.calls, 0, "coefficient-add trace supports one call");
            snapshot.calls = 1;
            snapshot.entry_ops_idx = circ.total_ops() as usize;
            snapshot.lender_mask = lenders.iter().fold(0u64, |mask, lane| {
                mask | 1u64.checked_shl(lane.id()).unwrap_or(0)
            });
            trace.set((enabled, snapshot));
        }
    });
}

fn record_coefficient_add_lender_restore(circ: &Circuit) {
    COEFFICIENT_ADD_LENDER_TRACE.with(|trace| {
        let (enabled, mut snapshot) = trace.get();
        if enabled {
            assert_eq!(snapshot.calls, 1, "coefficient-add restore without entry");
            snapshot.restore_ops_idx = circ.total_ops() as usize;
            trace.set((enabled, snapshot));
        }
    });
}

fn assert_clean_chain_coefficient_add_lender_preconditions(
    control: &QReg,
    work1: &[QReg],
    work2: &[QReg],
    l_t: &[QReg],
    lenders: &[&QReg],
) {
    assert_eq!(work1.len(), work2.len());
    assert!(!work1.is_empty());
    assert!(!l_t.is_empty());
    assert_eq!(
        lenders.len(),
        2,
        "clean-chain coefficient add requires exactly active and v-chain tmp lenders"
    );
    assert_ne!(
        lenders[0].id(),
        lenders[1].id(),
        "clean-chain coefficient-add lenders must be distinct"
    );
    for (index, lane) in lenders.iter().enumerate() {
        assert_ne!(lane.id(), control.id(), "lender {index} aliases control");
        assert!(
            work1.iter().all(|other| other.id() != lane.id()),
            "lender {index} aliases work1"
        );
        assert!(
            work2.iter().all(|other| other.id() != lane.id()),
            "lender {index} aliases work2"
        );
        assert!(
            l_t.iter().all(|other| other.id() != lane.id()),
            "lender {index} aliases l_t"
        );
    }
}

/// Apply only the coefficient data update.
///
/// When `reuse_clean_chain` is true, `lenders` must be two distinct clean
/// caller lanes ordered as `(active, v-chain tmp)`. Both lenders and the
/// in-place `l_t` cursor are restored exactly; the focused exhaustive proof
/// checks the entry and restore cut points.
fn coefficient_add_data_only_with_lane_route(
    circ: &mut Circuit,
    control: &QReg,
    work1: &[QReg],
    work2: &[QReg],
    l_t: &[QReg],
    inverse: bool,
    lenders: &[&QReg],
    reuse_clean_chain: bool,
) {
    if reuse_clean_chain {
        assert_clean_chain_coefficient_add_lender_preconditions(
            control, work1, work2, l_t, lenders,
        );
        record_coefficient_add_lender_entry(circ, lenders);
    }
    let cursor_scratch =
        circ.alloc_qreg_bits("rs.coeff-add.cursor-scratch", l_t.len().saturating_sub(1));
    let carry = circ.alloc_qreg("rs.coeff-add.carry");
    let active_owned = (!reuse_clean_chain).then(|| circ.alloc_qreg("rs.coeff-add.active"));
    let tmp_owned = (!reuse_clean_chain).then(|| circ.alloc_qreg("rs.coeff-add.tmp"));
    let active = if reuse_clean_chain {
        lenders[0]
    } else {
        active_owned.as_ref().expect("owned coefficient-add active")
    };
    let tmp = if reuse_clean_chain {
        lenders[1]
    } else {
        tmp_owned
            .as_ref()
            .expect("owned coefficient-add v-chain tmp")
    };
    let bracket = production_coefficient_nonnegative_bracket();

    if inverse {
        for index in 0..work1.len() {
            if index != 0 {
                decrement_mod_2n(circ, l_t, &cursor_scratch);
            }
            begin_coefficient_nonnegative_bracket(circ, control, l_t, active, bracket);
            circ.ccx(active, &work1[index], &work2[index]);
            circ.ccx(active, &carry, &work1[index]);
            multi_controlled_x_vchain(
                circ,
                &[active, &work1[index], &work2[index]],
                &carry,
                std::slice::from_ref(tmp),
            );
            end_coefficient_nonnegative_bracket(circ, control, l_t, active, bracket);
        }
        for index in (0..work1.len()).rev() {
            if index + 1 != work1.len() {
                increment_mod_2n(circ, l_t, &cursor_scratch);
            }
            begin_coefficient_nonnegative_bracket(circ, control, l_t, active, bracket);
            multi_controlled_x_vchain(
                circ,
                &[active, &work1[index], &work2[index]],
                &carry,
                std::slice::from_ref(tmp),
            );
            circ.ccx(active, &carry, &work1[index]);
            circ.ccx(active, &carry, &work2[index]);
            end_coefficient_nonnegative_bracket(circ, control, l_t, active, bracket);
        }
    } else {
        for index in 0..work1.len() {
            begin_coefficient_nonnegative_bracket(circ, control, l_t, active, bracket);
            circ.ccx(active, &carry, &work2[index]);
            circ.ccx(active, &carry, &work1[index]);
            multi_controlled_x_vchain(
                circ,
                &[active, &work1[index], &work2[index]],
                &carry,
                std::slice::from_ref(tmp),
            );
            end_coefficient_nonnegative_bracket(circ, control, l_t, active, bracket);
            if index + 1 != work1.len() {
                decrement_mod_2n(circ, l_t, &cursor_scratch);
            }
        }
        for index in (0..work1.len()).rev() {
            begin_coefficient_nonnegative_bracket(circ, control, l_t, active, bracket);
            multi_controlled_x_vchain(
                circ,
                &[active, &work1[index], &work2[index]],
                &carry,
                std::slice::from_ref(tmp),
            );
            circ.ccx(active, &carry, &work1[index]);
            circ.ccx(active, &work1[index], &work2[index]);
            end_coefficient_nonnegative_bracket(circ, control, l_t, active, bracket);
            if index != 0 {
                increment_mod_2n(circ, l_t, &cursor_scratch);
            }
        }
    }

    if reuse_clean_chain {
        record_coefficient_add_lender_restore(circ);
    }
    if let Some(tmp) = tmp_owned {
        circ.zero_and_free(tmp);
    }
    if let Some(active) = active_owned {
        circ.zero_and_free(active);
    }
    circ.zero_and_free(carry);
    free_clean(circ, cursor_scratch);
}

fn coefficient_add_data_only(
    circ: &mut Circuit,
    control: &QReg,
    work1: &[QReg],
    work2: &[QReg],
    l_t: &[QReg],
    inverse: bool,
    lenders: &[&QReg],
) {
    coefficient_add_data_only_with_lane_route(
        circ,
        control,
        work1,
        work2,
        l_t,
        inverse,
        lenders,
        clean_chain_coefficient_add_lender_requested(),
    );
}

#[allow(clippy::too_many_arguments)]
fn assert_q845_fused_guard_preconditions(
    control: &QReg,
    target_length: &[QReg],
    l_t: &[QReg],
    l_s: &[QReg],
    target: &QReg,
    scratch: &[&QReg],
) {
    assert!(!target_length.is_empty());
    assert_eq!(target_length.len(), l_t.len());
    assert_eq!(target_length.len(), l_s.len());
    assert_eq!(
        scratch.len(),
        3,
        "Q845 fused guard requires zero-pad, carry, and borrow lenders"
    );
    let mut ids = Vec::with_capacity(3 * target_length.len() + scratch.len() + 2);
    for (index, lane) in std::iter::once(control)
        .chain(std::iter::once(target))
        .chain(target_length)
        .chain(l_t)
        .chain(l_s)
        .chain(scratch.iter().copied())
        .enumerate()
    {
        assert!(
            !ids.contains(&lane.id()),
            "Q845 fused guard lane {index} aliases another operand or lender"
        );
        ids.push(lane.id());
    }
}

/// Toggle `target` when `control` is set and
/// `target_length > l_t + l_s + 1`.
///
/// The extra one is injected as the Cuccaro carry-in, while the zero-extension
/// lane receives the true carry-out. This avoids the `n - 1` clean lanes used
/// by a standalone wide increment. All three caller lenders and the allocated
/// boundary are restored exactly.
#[allow(clippy::too_many_arguments)]
fn toggle_q845_coefficient_length_above_guarded_boundary(
    circ: &mut Circuit,
    control: &QReg,
    target_length: &[QReg],
    l_t: &[QReg],
    l_s: &[QReg],
    target: &QReg,
    scratch: &[&QReg],
) {
    assert_q845_fused_guard_preconditions(control, target_length, l_t, l_s, target, scratch);
    let length_width = l_t.len();
    let zero_pad = scratch[0];
    let carry = scratch[1];
    let compare_borrow = scratch[2];
    let boundary = circ.alloc_qreg_bits("rs.q845-coeff-fused.boundary", length_width + 1);
    for (source, destination) in l_t.iter().zip(&boundary) {
        circ.cx(source, destination);
    }

    // boundary = zero_extend(l_t + l_s + 1). The Cuccaro carry lane is
    // restored to one by the adder and then returned to zero explicitly.
    circ.x(carry);
    cuccaro_add_mod_2n(
        circ,
        l_s,
        &boundary[..length_width],
        carry,
        &boundary[length_width],
    );
    circ.x(carry);

    let mut zero_extended_target: Vec<&QReg> = target_length.iter().collect();
    zero_extended_target.push(zero_pad);
    cuccaro_sub_mod_2n_refs(
        circ,
        &zero_extended_target,
        &boundary,
        carry,
        compare_borrow,
    );
    circ.ccx(control, compare_borrow, target);
    cuccaro_add_mod_2n_refs(
        circ,
        &zero_extended_target,
        &boundary,
        carry,
        compare_borrow,
    );

    circ.x(carry);
    cuccaro_sub_mod_2n(
        circ,
        l_s,
        &boundary[..length_width],
        carry,
        &boundary[length_width],
    );
    circ.x(carry);
    for (source, destination) in l_t.iter().zip(&boundary) {
        circ.cx(source, destination);
    }
    free_clean(circ, boundary);
}

fn increment_mod_2n_refs(circ: &mut Circuit, register: &[&QReg], carries: &[&QReg]) {
    let width = register.len();
    if width == 0 {
        return;
    }
    if width == 1 {
        circ.x(register[0]);
        return;
    }
    assert!(carries.len() >= width - 1);
    circ.cx(register[0], carries[0]);
    for index in 1..width - 1 {
        circ.ccx(register[index], carries[index - 1], carries[index]);
    }
    circ.cx(carries[width - 2], register[width - 1]);
    for index in (1..width - 1).rev() {
        circ.ccx(register[index], carries[index - 1], carries[index]);
        circ.cx(carries[index - 1], register[index]);
    }
    circ.cx(register[0], carries[0]);
    circ.x(register[0]);
}

/// Replace `register` by `constant - register (mod 2^n)` without allocating.
/// The transform is its own inverse. Addition of the classical constant is
/// decomposed into suffix increments, so the supplied clean scratch is reused
/// and restored after every term.
fn affine_complement_constant_refs(
    circ: &mut Circuit,
    register: &[&QReg],
    constant: usize,
    scratch: &[&QReg],
) {
    assert!(!register.is_empty());
    assert!(register.len() < usize::BITS as usize);
    assert!(scratch.len() >= register.len().saturating_sub(1));
    for lane in register {
        circ.x(lane);
    }
    let mask = (1usize << register.len()) - 1;
    let addend = constant.wrapping_add(1) & mask;
    for bit in 0..register.len() {
        if ((addend >> bit) & 1) != 0 {
            increment_mod_2n_refs(circ, &register[bit..], scratch);
        }
    }
}

fn increment_mod_2n_dirty_ladder_refs(
    circ: &mut Circuit,
    register: &[&QReg],
    dirty: &[&QReg],
) {
    assert!(!register.is_empty());
    for bit in (1..register.len()).rev() {
        mcx_dirty_ladder(circ, &register[..bit], register[bit], dirty);
    }
    circ.x(register[0]);
}

fn decrement_mod_2n_dirty_ladder_refs(
    circ: &mut Circuit,
    register: &[&QReg],
    dirty: &[&QReg],
) {
    assert!(!register.is_empty());
    circ.x(register[0]);
    for bit in 1..register.len() {
        mcx_dirty_ladder(circ, &register[..bit], register[bit], dirty);
    }
}

fn add_constant_dirty_refs(
    circ: &mut Circuit,
    register: &[&QReg],
    constant: usize,
    dirty: &[&QReg],
) {
    for bit in 0..register.len() {
        if ((constant >> bit) & 1) != 0 {
            increment_mod_2n_dirty_ladder_refs(circ, &register[bit..], dirty);
        }
    }
}

fn sub_constant_dirty_refs(
    circ: &mut Circuit,
    register: &[&QReg],
    constant: usize,
    dirty: &[&QReg],
) {
    for bit in (0..register.len()).rev() {
        if ((constant >> bit) & 1) != 0 {
            decrement_mod_2n_dirty_ladder_refs(circ, &register[bit..], dirty);
        }
    }
}

fn affine_complement_constant_dirty_refs(
    circ: &mut Circuit,
    register: &[&QReg],
    constant: usize,
    dirty: &[&QReg],
) {
    assert!(!register.is_empty());
    assert!(register.len() < usize::BITS as usize);
    assert!(dirty.len() >= register.len().saturating_sub(3));
    for lane in register {
        circ.x(lane);
    }
    let mask = (1usize << register.len()) - 1;
    let addend = constant.wrapping_add(1) & mask;
    add_constant_dirty_refs(circ, register, addend, dirty);
}

fn assert_controlled_short_carry_preconditions(
    control: &QReg,
    register: &[QReg],
    carries: &[QReg],
) {
    assert!(
        register.len() >= 2,
        "short-carry increment requires at least two bits"
    );
    assert_eq!(
        carries.len(),
        register.len() - 2,
        "short-carry increment requires exactly n-2 clean carry lanes"
    );
    for (index, lane) in register.iter().chain(carries).enumerate() {
        assert_ne!(lane.id(), control.id(), "short-carry lane {index} aliases control");
    }
    for (index, lane) in register.iter().enumerate() {
        assert!(register[..index]
            .iter()
            .all(|other| other.id() != lane.id()));
        assert!(carries.iter().all(|other| other.id() != lane.id()));
    }
    for (index, lane) in carries.iter().enumerate() {
        assert!(carries[..index]
            .iter()
            .all(|other| other.id() != lane.id()));
    }
}

/// Add `control` modulo `2^n` with `n-2` clean carries. The final ripple
/// carry is used directly as the second control of the top-bit Toffoli instead
/// of being materialized in a redundant `n-1`st lane.
fn controlled_increment_mod_2n_short_carry(
    circ: &mut Circuit,
    control: &QReg,
    register: &[QReg],
    carries: &[QReg],
) {
    assert_controlled_short_carry_preconditions(control, register, carries);
    if register.len() == 2 {
        circ.ccx(control, &register[0], &register[1]);
        circ.cx(control, &register[0]);
        return;
    }

    let last_carry = carries.len() - 1;
    circ.ccx(control, &register[0], &carries[0]);
    circ.cx(control, &register[0]);
    for index in 1..register.len() - 2 {
        circ.ccx(&register[index], &carries[index - 1], &carries[index]);
    }
    circ.ccx(
        &register[register.len() - 2],
        &carries[last_carry],
        &register[register.len() - 1],
    );
    circ.cx(&carries[last_carry], &register[register.len() - 2]);
    for index in (1..register.len() - 2).rev() {
        circ.ccx(&register[index], &carries[index - 1], &carries[index]);
        circ.cx(&carries[index - 1], &register[index]);
    }
    circ.cx(control, &register[0]);
    circ.ccx(control, &register[0], &carries[0]);
    circ.cx(control, &register[0]);
}

/// Exact gate-stream inverse of [`controlled_increment_mod_2n_short_carry`].
fn controlled_decrement_mod_2n_short_carry(
    circ: &mut Circuit,
    control: &QReg,
    register: &[QReg],
    carries: &[QReg],
) {
    assert_controlled_short_carry_preconditions(control, register, carries);
    if register.len() == 2 {
        circ.cx(control, &register[0]);
        circ.ccx(control, &register[0], &register[1]);
        return;
    }

    let last_carry = carries.len() - 1;
    circ.cx(control, &register[0]);
    circ.ccx(control, &register[0], &carries[0]);
    circ.cx(control, &register[0]);
    for index in 1..register.len() - 2 {
        circ.cx(&carries[index - 1], &register[index]);
        circ.ccx(&register[index], &carries[index - 1], &carries[index]);
    }
    circ.cx(&carries[last_carry], &register[register.len() - 2]);
    circ.ccx(
        &register[register.len() - 2],
        &carries[last_carry],
        &register[register.len() - 1],
    );
    for index in (1..register.len() - 2).rev() {
        circ.ccx(&register[index], &carries[index - 1], &carries[index]);
    }
    circ.cx(control, &register[0]);
    circ.ccx(control, &register[0], &carries[0]);
}

fn assert_q830_hybrid_block_counter_preconditions(
    control: &QReg,
    register: &[QReg],
    dirty: &[QReg],
    clean: &[QReg],
    boundary: &QReg,
) {
    assert_eq!(register.len(), REFERENCE_LENGTH_WIDTH);
    assert!(dirty.len() >= 4);
    assert_eq!(clean.len(), 2);
    let mut ids = Vec::new();
    for lane in std::iter::once(control)
        .chain(register)
        .chain(dirty)
        .chain(clean)
        .chain(std::iter::once(boundary))
    {
        assert!(!ids.contains(&lane.id()), "q830 hybrid counter lane alias");
        ids.push(lane.id());
    }
}

fn controlled_increment_mod_2n_hybrid_block(
    circ: &mut Circuit,
    control: &QReg,
    register: &[QReg],
    dirty: &[QReg],
    clean: &[QReg],
    boundary: &QReg,
) {
    const LOW_WIDTH: usize = 5;
    assert_q830_hybrid_block_counter_preconditions(control, register, dirty, clean, boundary);
    let dirty_refs = dirty.iter().collect::<Vec<_>>();
    let boundary_controls = std::iter::once(control)
        .chain(register[..LOW_WIDTH].iter())
        .collect::<Vec<_>>();
    mcx_dirty_ladder(circ, &boundary_controls, boundary, &dirty_refs);
    let high_carries = clean.iter().map(QReg::borrowed_alias).collect::<Vec<_>>();
    controlled_increment_mod_2n_short_carry(
        circ,
        boundary,
        &register[LOW_WIDTH..],
        &high_carries,
    );
    drop(high_carries);
    mcx_dirty_ladder(circ, &boundary_controls, boundary, &dirty_refs);
    let low_carries = std::iter::once(boundary.borrowed_alias())
        .chain(clean.iter().map(QReg::borrowed_alias))
        .collect::<Vec<_>>();
    controlled_increment_mod_2n_short_carry(
        circ,
        control,
        &register[..LOW_WIDTH],
        &low_carries,
    );
}

fn controlled_decrement_mod_2n_hybrid_block(
    circ: &mut Circuit,
    control: &QReg,
    register: &[QReg],
    dirty: &[QReg],
    clean: &[QReg],
    boundary: &QReg,
) {
    const LOW_WIDTH: usize = 5;
    assert_q830_hybrid_block_counter_preconditions(control, register, dirty, clean, boundary);
    let low_carries = std::iter::once(boundary.borrowed_alias())
        .chain(clean.iter().map(QReg::borrowed_alias))
        .collect::<Vec<_>>();
    controlled_decrement_mod_2n_short_carry(
        circ,
        control,
        &register[..LOW_WIDTH],
        &low_carries,
    );
    drop(low_carries);
    let dirty_refs = dirty.iter().collect::<Vec<_>>();
    let boundary_controls = std::iter::once(control)
        .chain(register[..LOW_WIDTH].iter())
        .collect::<Vec<_>>();
    mcx_dirty_ladder(circ, &boundary_controls, boundary, &dirty_refs);
    let high_carries = clean.iter().map(QReg::borrowed_alias).collect::<Vec<_>>();
    controlled_decrement_mod_2n_short_carry(
        circ,
        boundary,
        &register[LOW_WIDTH..],
        &high_carries,
    );
    mcx_dirty_ladder(circ, &boundary_controls, boundary, &dirty_refs);
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum Q839UlsCallbackDirection {
    Forward,
    Reverse,
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
struct Q839UlsCallbackAuditTrace {
    forward_calls: usize,
    reverse_calls: usize,
    forward_callback_slices: usize,
    reverse_callback_slices: usize,
    lender_ids_checked: usize,
    callback_ops_checked: usize,
    lender_role_checks: usize,
}

thread_local! {
    static Q839_ULS_CALLBACK_AUDIT_TRACE: std::cell::Cell<(
        bool,
        Q839UlsCallbackAuditTrace,
    )> = std::cell::Cell::new((false, Q839UlsCallbackAuditTrace {
        forward_calls: 0,
        reverse_calls: 0,
        forward_callback_slices: 0,
        reverse_callback_slices: 0,
        lender_ids_checked: 0,
        callback_ops_checked: 0,
        lender_role_checks: 0,
    }));
}

fn begin_q839_uls_callback_audit() {
    Q839_ULS_CALLBACK_AUDIT_TRACE.with(|trace| {
        trace.set((true, Q839UlsCallbackAuditTrace::default()));
    });
}

fn finish_q839_uls_callback_audit() -> Q839UlsCallbackAuditTrace {
    Q839_ULS_CALLBACK_AUDIT_TRACE.with(|trace| {
        let (_, snapshot) = trace.get();
        trace.set((false, snapshot));
        snapshot
    })
}

fn record_q839_uls_callback_call(
    direction: Q839UlsCallbackDirection,
    lenders: &[&QReg],
) -> bool {
    Q839_ULS_CALLBACK_AUDIT_TRACE.with(|trace| {
        let (enabled, mut snapshot) = trace.get();
        if !enabled {
            return false;
        }
        assert_eq!(lenders.len(), 2, "q837 ULS requires exactly two lenders");
        assert_ne!(lenders[0].id(), lenders[1].id());
        match direction {
            Q839UlsCallbackDirection::Forward => snapshot.forward_calls += 1,
            Q839UlsCallbackDirection::Reverse => snapshot.reverse_calls += 1,
        }
        snapshot.lender_ids_checked += lenders.len();
        trace.set((enabled, snapshot));
        true
    })
}

fn begin_q839_uls_callback_slice(circ: &Circuit) -> usize {
    assert!(!circ.b.count_only, "q837 callback audit requires emitted ops");
    circ.b.ops.len()
}

fn finish_q839_uls_callback_slice(
    circ: &Circuit,
    direction: Q839UlsCallbackDirection,
    lenders: &[&QReg],
    start_ops_idx: usize,
) {
    assert_eq!(lenders.len(), 2, "q837 ULS requires exactly two lenders");
    let ops = &circ.b.ops[start_ops_idx..];
    assert!(!ops.is_empty(), "q837 production callback emitted no operations");
    for op in ops {
        for lender in lenders {
            let lender_id = u64::from(lender.id());
            assert_ne!(op.q_control1.0, lender_id, "q837 ULS lender used as control1");
            assert_ne!(op.q_control2.0, lender_id, "q837 ULS lender used as control2");
            assert_ne!(op.q_target.0, lender_id, "q837 ULS lender used as target");
        }
    }
    Q839_ULS_CALLBACK_AUDIT_TRACE.with(|trace| {
        let (enabled, mut snapshot) = trace.get();
        assert!(enabled, "q837 callback audit stopped inside a callback");
        match direction {
            Q839UlsCallbackDirection::Forward => snapshot.forward_callback_slices += 1,
            Q839UlsCallbackDirection::Reverse => snapshot.reverse_callback_slices += 1,
        }
        snapshot.callback_ops_checked += ops.len();
        snapshot.lender_role_checks += ops.len() * lenders.len() * 3;
        trace.set((enabled, snapshot));
    });
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct Q839SevenPlateauLenderProofReport {
    pub short_carry_widths_checked: usize,
    pub short_carry_basis_states_checked: usize,
    pub short_carry_increment_checks: usize,
    pub short_carry_decrement_checks: usize,
    pub short_carry_inverse_checks: usize,
    pub short_carry_reversed_op_stream_checks: usize,
    pub uls_basis_states_checked: usize,
    pub uls_peak_lanes_removed: usize,
    pub uls_production_layouts_checked: usize,
    pub uls_forward_calls: usize,
    pub uls_reverse_calls: usize,
    pub uls_forward_callback_slices: usize,
    pub uls_reverse_callback_slices: usize,
    pub uls_callback_lender_ids_checked: usize,
    pub uls_callback_ops_checked: usize,
    pub uls_callback_lender_role_checks: usize,
    pub support_lender_basis_states_checked: usize,
    pub support_production_directions_checked: usize,
    pub support_borrow_windows_checked: usize,
    pub support_borrowed_owned_shots_checked: usize,
    pub support_lender_clean_entry_checks: usize,
    pub support_lender_restore_checks: usize,
    pub support_phase_clean_checks: usize,
    pub support_ancilla_clean_checks: usize,
    pub support_peak_lanes_removed: usize,
}

struct Q839SevenPlateauProofEnvironment {
    saved: Vec<(&'static str, Option<std::ffi::OsString>)>,
}

impl Q839SevenPlateauProofEnvironment {
    fn capture() -> Self {
        const NAMES: &[&str] = &[
            "POINT_ADD_COUNT_ONLY",
            "LOWQ_DIRECT_PREFIX_BITLEN",
            "LOWQ_REUSE_ZERO_CARRIES_FOR_PREFIX",
            "LOWQ_REUSE_ZERO_CARRIES_FOR_FULL_PREFIX_SCRATCH",
            "LOWQ_REUSE_ROTATED_BITLEN_SCRATCH",
            "LOWQ_REUSE_COEFFICIENT_COMPARATOR_SCRATCH",
            "LOWQ_FUSED_ZERO_PREFIX_BITLEN",
            "LOWQ_DIRECT_PREFIX_DIRTY_UPDATE",
            "LOWQ_DIRECT_PREFIX_NO_FLAG",
            "LOWQ_RS_DIRTY_ZERO_CORRECTION",
            "LOWQ_REGISTER_SHARED_REVERSE_DECREMENT_STREAM",
            "LOWQ_CALLER_SCRATCH_KG_REVERSE_DECREMENT",
            "LOWQ_SUB800_ULS_DIRTY_MCX3",
            COEFFICIENT_RAW_BITLEN_LOAN_FLAG,
            INPLACE_ROTATED_BITLEN_BOUNDARY_FLAG,
            FUSED_PREFIX_SCRATCH_LOAN_FLAG,
            PROMISED_LQ_SWAP_BORROW_FLAG,
            SPLIT_COEFFICIENT_ROTATION_LIFETIME_FLAG,
            COEFFICIENT_LESS_THAN_LANE_REUSE_FLAG,
            CLEAN_CHAIN_COEFFICIENT_ADD_LENDER_FLAG,
            PRESERVED_DY_TOP_PREFIX_LOAN_FLAG,
            MIXED_WIDTH_L_R_PRIME_FLAG,
            PAIRED_BITLEN_SOURCE_COMPLEMENT_FLAG,
            COEFFICIENT_NONNEGATIVE_X_CANCEL_FLAG,
            Q845_LIFETIME_COEFFICIENT_FUSION_FLAG,
            PROMISED_SWAP_SUPPORT_LIFETIME_FUSION_FLAG,
            Q845_SWAP_ONLY_T_PRIME_LENGTH_FLAG,
            Q851_TRUNCATED_SWAP_ONLY_GUARD_FLAG,
            Q851_FIXED_SIGN_EVENT_FLAG,
            SUB800_INPLACE_GUARD_ADDRESS_FLAG,
            SUB800_RAW_PREFIX_PRESERVED_LENDER_FLAG,
            SUB800_RAW_PREFIX_PREDICATE_LENDER_FLAG,
            SUB800_MIXED_BOUNDARY_SCRATCH_EXTENSION_FLAG,
            SUB800_BORROWED_ROTATED_UNDERFLOW_FLAG,
            SUB800_SPLIT_MIXED_ROTATED_LENGTH_FLAG,
            SUB800_SPLIT_SAME_ROTATED_LENGTH_FLAG,
            SUB800_SPLIT_TWO_HIGH_ROTATED_LENGTH_FLAG,
            SUB800_SPLIT_THREE_HIGH_ROTATED_LENGTH_FLAG,
            SUB800_SPLIT_FOUR_HIGH_ROTATED_LENGTH_FLAG,
            SUB800_ULS_CLEAN_LENDER_FLAG,
            SUB800_ULS_FUSED_TARGET_FLAG,
            SUB800_ULS_DIRECT_SELECTOR_FLAG,
            Q839_SEVEN_PLATEAU_LENDERS_FLAG,
        ];
        Self {
            saved: NAMES
                .iter()
                .copied()
                .map(|name| (name, std::env::var_os(name)))
                .collect(),
        }
    }
}

impl Drop for Q839SevenPlateauProofEnvironment {
    fn drop(&mut self) {
        for (name, value) in self.saved.drain(..) {
            match value {
                Some(value) => std::env::set_var(name, value),
                None => std::env::remove_var(name),
            }
        }
    }
}

fn configure_q839_seven_plateau_proof_environment() {
    std::env::remove_var("POINT_ADD_COUNT_ONLY");
    for name in [
        "LOWQ_DIRECT_PREFIX_BITLEN",
        "LOWQ_REUSE_ZERO_CARRIES_FOR_PREFIX",
        "LOWQ_REUSE_ZERO_CARRIES_FOR_FULL_PREFIX_SCRATCH",
        "LOWQ_REUSE_ROTATED_BITLEN_SCRATCH",
        "LOWQ_REUSE_COEFFICIENT_COMPARATOR_SCRATCH",
        "LOWQ_FUSED_ZERO_PREFIX_BITLEN",
        "LOWQ_REGISTER_SHARED_REVERSE_DECREMENT_STREAM",
        "LOWQ_CALLER_SCRATCH_KG_REVERSE_DECREMENT",
        "LOWQ_SUB800_ULS_DIRTY_MCX3",
        COEFFICIENT_RAW_BITLEN_LOAN_FLAG,
        INPLACE_ROTATED_BITLEN_BOUNDARY_FLAG,
        FUSED_PREFIX_SCRATCH_LOAN_FLAG,
        PROMISED_LQ_SWAP_BORROW_FLAG,
        SPLIT_COEFFICIENT_ROTATION_LIFETIME_FLAG,
        COEFFICIENT_LESS_THAN_LANE_REUSE_FLAG,
        CLEAN_CHAIN_COEFFICIENT_ADD_LENDER_FLAG,
        PRESERVED_DY_TOP_PREFIX_LOAN_FLAG,
        MIXED_WIDTH_L_R_PRIME_FLAG,
        PAIRED_BITLEN_SOURCE_COMPLEMENT_FLAG,
        COEFFICIENT_NONNEGATIVE_X_CANCEL_FLAG,
        Q845_LIFETIME_COEFFICIENT_FUSION_FLAG,
        PROMISED_SWAP_SUPPORT_LIFETIME_FUSION_FLAG,
        Q845_SWAP_ONLY_T_PRIME_LENGTH_FLAG,
        Q851_TRUNCATED_SWAP_ONLY_GUARD_FLAG,
        Q851_FIXED_SIGN_EVENT_FLAG,
        SUB800_INPLACE_GUARD_ADDRESS_FLAG,
        SUB800_RAW_PREFIX_PRESERVED_LENDER_FLAG,
        SUB800_RAW_PREFIX_PREDICATE_LENDER_FLAG,
        SUB800_MIXED_BOUNDARY_SCRATCH_EXTENSION_FLAG,
        SUB800_BORROWED_ROTATED_UNDERFLOW_FLAG,
        SUB800_SPLIT_MIXED_ROTATED_LENGTH_FLAG,
        SUB800_SPLIT_SAME_ROTATED_LENGTH_FLAG,
        SUB800_SPLIT_TWO_HIGH_ROTATED_LENGTH_FLAG,
        SUB800_SPLIT_THREE_HIGH_ROTATED_LENGTH_FLAG,
        SUB800_SPLIT_FOUR_HIGH_ROTATED_LENGTH_FLAG,
        SUB800_ULS_CLEAN_LENDER_FLAG,
        SUB800_ULS_FUSED_TARGET_FLAG,
        Q839_SEVEN_PLATEAU_LENDERS_FLAG,
    ] {
        std::env::set_var(name, "1");
    }
    for name in [
        "LOWQ_DIRECT_PREFIX_DIRTY_UPDATE",
        "LOWQ_DIRECT_PREFIX_NO_FLAG",
        "LOWQ_RS_DIRTY_ZERO_CORRECTION",
        SUB800_ULS_DIRECT_SELECTOR_FLAG,
    ] {
        std::env::remove_var(name);
    }
}

fn build_q839_uls_production_callback_harness(inverse: bool) -> B {
    const PHYSICAL_WORK_WIDTH: usize = 259;
    const SCAN_WIDTH: usize = 257;

    let mut circ = q845_fusion_proof_circuit();
    let phase1 = circ.alloc_qreg("q837.uls-audit.phase1");
    let phase2 = circ.alloc_qreg("q837.uls-audit.phase2");
    let sign = circ.alloc_qreg("q837.uls-audit.sign");
    let work1 = circ.alloc_qreg_bits("q837.uls-audit.work1", PHYSICAL_WORK_WIDTH);
    let work2 = circ.alloc_qreg_bits("q837.uls-audit.work2", PHYSICAL_WORK_WIDTH);
    let l_t = circ.alloc_qreg_bits("q837.uls-audit.l-t", REFERENCE_LENGTH_WIDTH);
    let l_t_prime =
        circ.alloc_qreg_bits("q837.uls-audit.l-t-prime", REFERENCE_LENGTH_WIDTH);
    let l_s = circ.alloc_qreg_bits("q837.uls-audit.l-s", REFERENCE_LENGTH_WIDTH);
    let l_r_prime = circ.alloc_qreg_bits("q837.uls-audit.l-r-prime", REFERENCE_R_LENGTH_WIDTH);
    let enable = circ.alloc_qreg("q837.uls-audit.enable");
    let above_guard = circ.alloc_qreg("q837.uls-audit.above-guard");
    let add_only = circ.alloc_qreg("q837.uls-audit.add-only");
    let chain = circ.alloc_qreg_bits("q837.uls-audit.chain", 2);
    let support_lender = circ.alloc_qreg("q837.uls-audit.support-lender");

    coefficient_fused_data_and_sign_q845_swap_only(
        &mut circ,
        &phase1,
        &phase2,
        &sign,
        &work1[..SCAN_WIDTH],
        &work2[..SCAN_WIDTH],
        PHYSICAL_WORK_WIDTH,
        &l_t,
        &l_t_prime,
        &l_s,
        &l_r_prime,
        &enable,
        &above_guard,
        &add_only,
        &chain,
        inverse,
        Some(&support_lender),
        None,
    );
    circ.into_builder()
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
struct Q839SupportProductionAudit {
    borrowed_owned_shots_checked: usize,
    lender_clean_entry_checks: usize,
    lender_restore_checks: usize,
    phase_clean_checks: usize,
    ancilla_clean_checks: usize,
}

impl Q839SupportProductionAudit {
    fn add(&mut self, other: Self) {
        self.borrowed_owned_shots_checked += other.borrowed_owned_shots_checked;
        self.lender_clean_entry_checks += other.lender_clean_entry_checks;
        self.lender_restore_checks += other.lender_restore_checks;
        self.phase_clean_checks += other.phase_clean_checks;
        self.ancilla_clean_checks += other.ancilla_clean_checks;
    }
}

fn q839_support_sample_mask(id: u32) -> u64 {
    let mut value = u64::from(id).wrapping_add(0x9e37_79b9_7f4a_7c15);
    value = (value ^ (value >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    value = (value ^ (value >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    (value ^ (value >> 31)) & !1
}

fn verify_q839_support_production_harness(
    label: &[u8],
    owned: &PromisedLqSwapProofHarness,
    borrowed: &PromisedLqSwapProofHarness,
) -> Q839SupportProductionAudit {
    use crate::circuit::QubitId;
    use crate::sim::Simulator;
    use sha3::{
        digest::{ExtendableOutput, Update},
        Shake128,
    };

    const SHOTS: usize = 64;
    const FIRST_ZERO_PREDICATE_SHOTS: u64 = 0xff;
    const ZERO_Q_NONZERO_S_SHOT: u64 = 1 << 8;
    const NONZERO_Q_ZERO_S_SHOT: u64 = 1 << 9;
    const ZERO_Q_ZERO_S_SHOT: u64 = 1 << 10;
    const NONZERO_Q_NONZERO_S_SHOT: u64 = 1 << 11;
    const FOUR_ZERO_PREDICATE_CLASSES: u64 = ZERO_Q_NONZERO_S_SHOT
        | NONZERO_Q_ZERO_S_SHOT
        | ZERO_Q_ZERO_S_SHOT
        | NONZERO_Q_NONZERO_S_SHOT;

    assert_eq!(owned.data_ids, borrowed.data_ids);
    assert_eq!(owned.l_t_prime_ids, borrowed.l_t_prime_ids);
    assert_eq!(owned.l_q_ids, borrowed.l_q_ids);
    assert_eq!(owned.l_s_ids, borrowed.l_s_ids);
    assert_eq!(owned.preserved_dy_top_id, borrowed.preserved_dy_top_id);
    assert!(owned.q839_support_lender_windows.is_empty());
    assert_eq!(borrowed.q839_support_lender_windows.len(), 1);
    let window = borrowed.q839_support_lender_windows[0];
    assert_eq!(window.lender_id, borrowed.preserved_dy_top_id);
    assert_eq!(owned.builder.peak_qubits, borrowed.builder.peak_qubits + 1);

    let mut owned_seed = Shake128::default();
    owned_seed.update(label);
    owned_seed.update(b"-owned");
    let mut owned_xof = owned_seed.finalize_xof();
    let mut borrowed_seed = Shake128::default();
    borrowed_seed.update(label);
    borrowed_seed.update(b"-borrowed");
    let mut borrowed_xof = borrowed_seed.finalize_xof();
    let mut owned_simulator = Simulator::new(
        owned.builder.next_qubit as usize,
        owned.builder.next_bit as usize,
        &mut owned_xof,
    );
    let mut borrowed_simulator = Simulator::new(
        borrowed.builder.next_qubit as usize,
        borrowed.builder.next_bit as usize,
        &mut borrowed_xof,
    );
    for &id in &owned.data_ids {
        let mut mask = q839_support_sample_mask(id);
        if owned.l_t_prime_ids.contains(&id) {
            mask = 0;
        } else if owned.l_q_ids.contains(&id) {
            mask &= !(FIRST_ZERO_PREDICATE_SHOTS
                | ZERO_Q_NONZERO_S_SHOT
                | ZERO_Q_ZERO_S_SHOT);
            mask |= NONZERO_Q_ZERO_S_SHOT | NONZERO_Q_NONZERO_S_SHOT;
        } else if owned.l_s_ids.contains(&id) {
            mask &= !(FIRST_ZERO_PREDICATE_SHOTS
                | NONZERO_Q_ZERO_S_SHOT
                | ZERO_Q_ZERO_S_SHOT);
            mask |= ZERO_Q_NONZERO_S_SHOT | NONZERO_Q_NONZERO_S_SHOT;
        }
        *owned_simulator.qubit_mut(QubitId(u64::from(id))) = mask;
        *borrowed_simulator.qubit_mut(QubitId(u64::from(id))) = mask;
    }
    let zero_predicate_mask = |ids: &[u32]| {
        ids.iter().fold(u64::MAX, |mask, &id| {
            mask & !owned_simulator.qubit(QubitId(u64::from(id)))
        })
    };
    let zero_q_shots = zero_predicate_mask(&owned.l_q_ids);
    let zero_s_shots = zero_predicate_mask(&owned.l_s_ids);
    assert_eq!(
        zero_q_shots,
        FIRST_ZERO_PREDICATE_SHOTS | ZERO_Q_NONZERO_S_SHOT | ZERO_Q_ZERO_S_SHOT,
        "{label:?} q837 support Zq shots"
    );
    assert_eq!(
        zero_s_shots,
        FIRST_ZERO_PREDICATE_SHOTS | NONZERO_Q_ZERO_S_SHOT | ZERO_Q_ZERO_S_SHOT,
        "{label:?} q837 support Zs shots"
    );
    assert_eq!(
        zero_q_shots & !zero_s_shots,
        ZERO_Q_NONZERO_S_SHOT,
        "{label:?} q837 support omitted Zq=1,Zs=0 shot"
    );
    assert_eq!(
        (!zero_q_shots & zero_s_shots) & FOUR_ZERO_PREDICATE_CLASSES,
        NONZERO_Q_ZERO_S_SHOT,
        "{label:?} q837 support omitted Zq=0,Zs=1 shot"
    );
    assert_eq!(
        (zero_q_shots & zero_s_shots) & FOUR_ZERO_PREDICATE_CLASSES,
        ZERO_Q_ZERO_S_SHOT,
        "{label:?} q837 support omitted Zq=1,Zs=1 shot"
    );
    assert_eq!(
        (!zero_q_shots & !zero_s_shots) & FOUR_ZERO_PREDICATE_CLASSES,
        NONZERO_Q_NONZERO_S_SHOT,
        "{label:?} q837 support omitted Zq=0,Zs=0 shot"
    );

    borrowed_simulator.apply_iter(borrowed.builder.ops[..window.entry_ops_idx].iter());
    assert_eq!(
        borrowed_simulator.qubit(QubitId(u64::from(window.lender_id))),
        0,
        "{label:?} q837 support lender dirty at entry"
    );
    assert_eq!(
        borrowed_simulator.phase, 0,
        "{label:?} q837 support phase dirty at entry"
    );
    borrowed_simulator.apply_iter(
        borrowed.builder.ops[window.entry_ops_idx..window.restore_ops_idx].iter(),
    );
    assert_eq!(
        borrowed_simulator.qubit(QubitId(u64::from(window.lender_id))),
        0,
        "{label:?} q837 support lender dirty at restore"
    );
    assert_eq!(
        borrowed_simulator.phase, 0,
        "{label:?} q837 support phase dirty at restore"
    );
    borrowed_simulator.apply_iter(borrowed.builder.ops[window.restore_ops_idx..].iter());
    owned_simulator.apply_iter(owned.builder.ops.iter());
    assert_eq!(owned_simulator.phase, 0, "{label:?} owned support phase");
    assert_eq!(borrowed_simulator.phase, 0, "{label:?} borrowed support phase");

    for &id in &owned.data_ids {
        assert_eq!(
            owned_simulator.qubit(QubitId(u64::from(id))),
            borrowed_simulator.qubit(QubitId(u64::from(id))),
            "{label:?} borrowed support changed external q{id}"
        );
    }
    assert_eq!(
        owned_simulator.qubit(QubitId(u64::from(owned.preserved_dy_top_id))),
        0
    );
    assert_eq!(
        borrowed_simulator.qubit(QubitId(u64::from(borrowed.preserved_dy_top_id))),
        0
    );

    let mut owned_external = vec![false; owned.builder.next_qubit as usize];
    let mut borrowed_external = vec![false; borrowed.builder.next_qubit as usize];
    for &id in &owned.data_ids {
        owned_external[id as usize] = true;
        borrowed_external[id as usize] = true;
    }
    owned_external[owned.preserved_dy_top_id as usize] = true;
    borrowed_external[borrowed.preserved_dy_top_id as usize] = true;
    for id in 0..owned.builder.next_qubit {
        if !owned_external[id as usize] {
            assert_eq!(
                owned_simulator.qubit(QubitId(u64::from(id))),
                0,
                "{label:?} owned support left q{id} dirty"
            );
        }
    }
    for id in 0..borrowed.builder.next_qubit {
        if !borrowed_external[id as usize] {
            assert_eq!(
                borrowed_simulator.qubit(QubitId(u64::from(id))),
                0,
                "{label:?} borrowed support left q{id} dirty"
            );
        }
    }

    Q839SupportProductionAudit {
        borrowed_owned_shots_checked: SHOTS,
        lender_clean_entry_checks: SHOTS,
        lender_restore_checks: SHOTS,
        phase_clean_checks: 2 * SHOTS,
        ancilla_clean_checks: 2 * SHOTS,
    }
}

/// Diagnostic proof for the lender primitives used by the gated seven-plateau
/// cut. It exhausts short carries through production width 9, audits the actual
/// production coefficient callback slices against both ULS lenders, and checks
/// the production mixed-width swap with its preserved-top support lender.
#[doc(hidden)]
#[must_use]
pub fn exhaustive_q839_seven_plateau_lender_check() -> Q839SevenPlateauLenderProofReport {
    use crate::point_add::trailmix_port::arith::khattar_gidney::{
        unary_iterate_log_star_toggle_target_with_clean_lenders,
    };

    let _environment = Q839SevenPlateauProofEnvironment::capture();
    configure_q839_seven_plateau_proof_environment();

    let mut short_carry_basis_states_checked = 0usize;
    let mut short_carry_increment_checks = 0usize;
    let mut short_carry_decrement_checks = 0usize;
    let mut short_carry_inverse_checks = 0usize;
    let mut short_carry_reversed_op_stream_checks = 0usize;
    for width in 2usize..=REFERENCE_LENGTH_WIDTH {
        let build = |decrement: bool| {
            let mut circ = Circuit::new();
            let control = circ.alloc_input_qreg_bits("q839.short-carry.control", 1);
            let register = circ.alloc_input_qreg_bits("q839.short-carry.register", width);
            let carries = circ.alloc_qreg_bits("q839.short-carry.carries", width - 2);
            let lender = circ.alloc_input_qreg_bits("q839.short-carry.uls-lender", 1);
            if decrement {
                controlled_decrement_mod_2n_short_carry(
                    &mut circ,
                    &control[0],
                    &register,
                    &carries,
                );
            } else {
                controlled_increment_mod_2n_short_carry(
                    &mut circ,
                    &control[0],
                    &register,
                    &carries,
                );
            }
            (
                circ.into_builder(),
                control[0].id(),
                register.iter().map(QReg::id).collect::<Vec<_>>(),
                carries.iter().map(QReg::id).collect::<Vec<_>>(),
                lender[0].id(),
            )
        };
        let (increment, control_id, register_ids, carry_ids, lender_id) = build(false);
        let (
            decrement,
            decrement_control_id,
            decrement_register_ids,
            decrement_carry_ids,
            decrement_lender_id,
        ) = build(true);
        assert_eq!(control_id, decrement_control_id);
        assert_eq!(register_ids, decrement_register_ids);
        assert_eq!(carry_ids, decrement_carry_ids);
        assert_eq!(lender_id, decrement_lender_id);
        assert_eq!(
            decrement.ops,
            increment.ops.iter().rev().copied().collect::<Vec<_>>()
        );
        short_carry_reversed_op_stream_checks += 1;
        let mask = (1usize << width) - 1;
        for control in 0usize..=1 {
            for value in 0usize..=mask {
                for lender in 0usize..=1 {
                    let mut input = (control as u64) << control_id;
                    input |= (lender as u64) << lender_id;
                    for (bit, id) in register_ids.iter().enumerate() {
                        input |= (((value >> bit) & 1) as u64) << id;
                    }
                    let incremented = apply_scalar(&increment.ops, input);
                    let actual = register_ids.iter().enumerate().fold(
                        0usize,
                        |word, (bit, id)| {
                            word | ((((incremented >> id) & 1) as usize) << bit)
                        },
                    );
                    assert_eq!(
                        actual,
                        value.wrapping_add(control) & mask,
                        "short-carry width={width} control={control} value={value} lender={lender}"
                    );
                    assert_eq!((incremented >> lender_id) & 1, lender as u64);
                    assert!(carry_ids
                        .iter()
                        .all(|id| ((incremented >> id) & 1) == 0));
                    let decremented = apply_scalar(&decrement.ops, input);
                    let actual_decrement = register_ids.iter().enumerate().fold(
                        0usize,
                        |word, (bit, id)| {
                            word | ((((decremented >> id) & 1) as usize) << bit)
                        },
                    );
                    assert_eq!(
                        actual_decrement,
                        value.wrapping_sub(control) & mask,
                        "short-carry decrement width={width} control={control} value={value} lender={lender}"
                    );
                    assert_eq!((decremented >> lender_id) & 1, lender as u64);
                    assert!(carry_ids
                        .iter()
                        .all(|id| ((decremented >> id) & 1) == 0));
                    assert_eq!(apply_scalar(&decrement.ops, incremented), input);
                    assert_eq!(apply_scalar(&increment.ops, decremented), input);
                    short_carry_basis_states_checked += 1;
                    short_carry_increment_checks += 1;
                    short_carry_decrement_checks += 1;
                    short_carry_inverse_checks += 2;
                }
            }
        }
    }

    let build_uls = |lender_count: usize| {
        let mut circ = Circuit::new();
        let lenders = circ.alloc_qreg_bits("q839.uls-proof.lenders", 2);
        let counter = circ.alloc_input_qreg_bits("q839.uls-proof.counter", 9);
        let target = circ.alloc_input_qreg_bits("q839.uls-proof.target", 1);
        let counter_refs: Vec<&QReg> = counter.iter().collect();
        let lender_refs: Vec<&QReg> = lenders.iter().take(lender_count).collect();
        unary_iterate_log_star_toggle_target_with_clean_lenders(
            &mut circ,
            &counter_refs,
            16,
            &lender_refs,
            &target[0],
            true,
            |_, _| {},
        );
        let external_ids: Vec<u32> = lenders
            .iter()
            .chain(&counter)
            .chain(std::iter::once(&target[0]))
            .map(QReg::id)
            .collect();
        let external_mask = external_ids
            .iter()
            .fold(0u64, |mask, id| mask | (1u64 << id));
        (circ.into_builder(), external_ids, external_mask)
    };
    let (uls_baseline, baseline_ids, baseline_mask) = build_uls(0);
    let (uls_candidate, candidate_ids, candidate_mask) = build_uls(2);
    assert_eq!(baseline_ids, candidate_ids);
    assert_eq!(baseline_mask, candidate_mask);
    assert_eq!(uls_candidate.peak_qubits + 2, uls_baseline.peak_qubits);
    let mut uls_basis_states_checked = 0usize;
    for value in 0u64..(1u64 << 10) {
        let mut input = 0u64;
        for (bit, id) in baseline_ids[2..].iter().enumerate() {
            input |= ((value >> bit) & 1) << id;
        }
        let baseline_output = apply_scalar(&uls_baseline.ops, input);
        let candidate_output = apply_scalar(&uls_candidate.ops, input);
        assert_eq!(baseline_output & baseline_mask, candidate_output & candidate_mask);
        assert_eq!(baseline_output & !baseline_mask, 0);
        assert_eq!(candidate_output & !candidate_mask, 0);
        uls_basis_states_checked += 1;
    }

    reset_sub800_q839_route_coverage();
    begin_q839_uls_callback_audit();
    let production_forward = build_q839_uls_production_callback_harness(false);
    let production_inverse = build_q839_uls_production_callback_harness(true);
    let uls_callback_audit = finish_q839_uls_callback_audit();
    assert!(!production_forward.ops.is_empty());
    assert!(!production_inverse.ops.is_empty());
    let production_callback_slices_per_direction = 2 * 257;
    assert_eq!(uls_callback_audit.forward_calls, 2);
    assert_eq!(uls_callback_audit.reverse_calls, 2);
    assert_eq!(
        uls_callback_audit.forward_callback_slices,
        production_callback_slices_per_direction
    );
    assert_eq!(
        uls_callback_audit.reverse_callback_slices,
        production_callback_slices_per_direction
    );
    assert_eq!(uls_callback_audit.lender_ids_checked, 8);
    assert!(uls_callback_audit.callback_ops_checked > 0);
    assert_eq!(
        uls_callback_audit.lender_role_checks,
        uls_callback_audit.callback_ops_checked * 2 * 3
    );
    let callback_coverage = sub800_q839_route_coverage();
    assert_eq!(callback_coverage.seven_plateau_uls_forward_loans, 2);
    assert_eq!(callback_coverage.seven_plateau_uls_reverse_loans, 2);
    assert_eq!(
        callback_coverage.seven_plateau_short_increments,
        production_callback_slices_per_direction
    );
    assert_eq!(
        callback_coverage.seven_plateau_short_decrements,
        production_callback_slices_per_direction
    );

    let mut support_lender_basis_states_checked = 0usize;
    let mut support = Circuit::new();
    let zero_q = support.alloc_input_qreg_bits("q839.support-proof.zero-q", 1);
    let zero_s = support.alloc_input_qreg_bits("q839.support-proof.zero-s", 1);
    let target = support.alloc_input_qreg_bits("q839.support-proof.target", 1);
    let lender = support.alloc_qreg("q839.support-proof.lender");
    support.ccx(&zero_q[0], &zero_s[0], &lender);
    support.cx(&lender, &target[0]);
    support.ccx(&zero_q[0], &zero_s[0], &lender);
    let support_ops = support.into_builder().ops;
    for input in 0u64..8 {
        let encoded = (input & 1) << zero_q[0].id()
            | ((input >> 1) & 1) << zero_s[0].id()
            | ((input >> 2) & 1) << target[0].id();
        let output = apply_scalar(&support_ops, encoded);
        assert_eq!((output >> lender.id()) & 1, 0);
        assert_eq!(
            (output >> target[0].id()) & 1,
            ((input >> 2) & 1) ^ ((input & 1) & ((input >> 1) & 1))
        );
        support_lender_basis_states_checked += 1;
    }

    let mut support_production_audit = Q839SupportProductionAudit::default();
    for (inverse, label) in [
        (false, b"q837-support-forward".as_slice()),
        (true, b"q837-support-inverse".as_slice()),
    ] {
        std::env::remove_var(Q839_SEVEN_PLATEAU_LENDERS_FLAG);
        let owned = build_promised_l_q_swap_proof_harness_with_widths(
            259,
            REFERENCE_LENGTH_WIDTH,
            REFERENCE_R_LENGTH_WIDTH,
            inverse,
            PromisedLqSwapRoute::BorrowLq,
            true,
        );
        std::env::set_var(Q839_SEVEN_PLATEAU_LENDERS_FLAG, "1");
        let borrowed = build_promised_l_q_swap_proof_harness_with_widths(
            259,
            REFERENCE_LENGTH_WIDTH,
            REFERENCE_R_LENGTH_WIDTH,
            inverse,
            PromisedLqSwapRoute::BorrowLq,
            true,
        );
        support_production_audit.add(verify_q839_support_production_harness(
            label, &owned, &borrowed,
        ));
    }
    let support_coverage = sub800_q839_route_coverage();
    assert_eq!(support_coverage.seven_plateau_support_loans, 2);
    assert_eq!(support_coverage.seven_plateau_support_fallbacks, 0);

    Q839SevenPlateauLenderProofReport {
        short_carry_widths_checked: REFERENCE_LENGTH_WIDTH - 1,
        short_carry_basis_states_checked,
        short_carry_increment_checks,
        short_carry_decrement_checks,
        short_carry_inverse_checks,
        short_carry_reversed_op_stream_checks,
        uls_basis_states_checked,
        uls_peak_lanes_removed: 2,
        uls_production_layouts_checked: 2,
        uls_forward_calls: uls_callback_audit.forward_calls,
        uls_reverse_calls: uls_callback_audit.reverse_calls,
        uls_forward_callback_slices: uls_callback_audit.forward_callback_slices,
        uls_reverse_callback_slices: uls_callback_audit.reverse_callback_slices,
        uls_callback_lender_ids_checked: uls_callback_audit.lender_ids_checked,
        uls_callback_ops_checked: uls_callback_audit.callback_ops_checked,
        uls_callback_lender_role_checks: uls_callback_audit.lender_role_checks,
        support_lender_basis_states_checked,
        support_production_directions_checked: 2,
        support_borrow_windows_checked: 2,
        support_borrowed_owned_shots_checked: support_production_audit
            .borrowed_owned_shots_checked,
        support_lender_clean_entry_checks: support_production_audit.lender_clean_entry_checks,
        support_lender_restore_checks: support_production_audit.lender_restore_checks,
        support_phase_clean_checks: support_production_audit.phase_clean_checks,
        support_ancilla_clean_checks: support_production_audit.ancilla_clean_checks,
        support_peak_lanes_removed: 1,
    }
}

#[cfg(test)]
mod q839_seven_plateau_lender_tests {
    #[test]
    fn production_lender_proof_hardening() {
        let report = super::exhaustive_q839_seven_plateau_lender_check();
        assert_eq!(report.short_carry_widths_checked, 8);
        assert_eq!(report.short_carry_basis_states_checked, 4_080);
        assert_eq!(report.short_carry_reversed_op_stream_checks, 8);
        assert_eq!(report.uls_forward_calls, report.uls_reverse_calls);
        assert_eq!(
            report.uls_forward_callback_slices,
            report.uls_reverse_callback_slices
        );
        assert_eq!(report.uls_peak_lanes_removed, 2);
        assert_eq!(report.support_borrow_windows_checked, 2);
        assert_eq!(report.support_peak_lanes_removed, 1);
    }
}

enum Q845SwapOnlyAddress<'a> {
    Allocated(Vec<QReg>),
    InPlaceLS(&'a [QReg]),
}

impl Q845SwapOnlyAddress<'_> {
    fn lanes(&self) -> &[QReg] {
        match self {
            Self::Allocated(address) => address,
            Self::InPlaceLS(address) => address,
        }
    }

    fn is_in_place(&self) -> bool {
        matches!(self, Self::InPlaceLS(_))
    }
}

/// Dynamic coefficient guard used when `l_t_prime` is intentionally zero
/// between swap sites. The global cache becomes a reversible population
/// counter. The baseline allocates a nine-lane address; the sub-800 candidate
/// stores the phase-local affine coordinates directly in `l_s` and restores
/// `l_s` before returning.
struct Q845SwapOnlyCoefficientGuard<'a> {
    address: Q845SwapOnlyAddress<'a>,
    physical_work_width: usize,
    uls_clean_lender: Option<&'a QReg>,
    relocation_l_q: Option<&'a [QReg]>,
}

fn toggle_register_geq_constant_vchain(
    circ: &mut Circuit,
    register: &[QReg],
    value: usize,
    target: &QReg,
    scratch: &[QReg],
) {
    let width = register.len();
    assert!(width < usize::BITS as usize);
    let modulus = 1usize << width;
    if value == 0 {
        circ.x(target);
        return;
    }
    if value >= modulus {
        return;
    }
    assert!(scratch.len() >= width.saturating_sub(2));
    let scratch_refs = scratch.iter().collect::<Vec<_>>();

    let mut toggle_pattern = |required: &[(usize, bool)]| {
        for &(index, expected) in required {
            if !expected {
                circ.x(&register[index]);
            }
        }
        let controls = required
            .iter()
            .map(|&(index, _)| &register[index])
            .collect::<Vec<_>>();
        multi_controlled_x_vchain_borrowed(circ, &controls, target, &scratch_refs);
        for &(index, expected) in required.iter().rev() {
            if !expected {
                circ.x(&register[index]);
            }
        }
    };

    // These terms are disjoint: either the register equals the constant, or
    // its highest differing bit is one where the constant has zero.
    for differing_bit in (0..width).rev() {
        if ((value >> differing_bit) & 1) != 0 {
            continue;
        }
        let mut required = Vec::with_capacity(width - differing_bit);
        required.push((differing_bit, true));
        for high_bit in differing_bit + 1..width {
            required.push((high_bit, ((value >> high_bit) & 1) != 0));
        }
        toggle_pattern(&required);
    }
    let equality = (0..width)
        .map(|index| (index, ((value >> index) & 1) != 0))
        .collect::<Vec<_>>();
    toggle_pattern(&equality);
}

fn toggle_register_geq_constant_dirty(
    circ: &mut Circuit,
    register: &[QReg],
    value: usize,
    target: &QReg,
    dirty: &[QReg],
) {
    let width = register.len();
    assert!(width < usize::BITS as usize);
    let modulus = 1usize << width;
    if value == 0 {
        circ.x(target);
        return;
    }
    if value >= modulus {
        return;
    }
    assert!(dirty.len() >= width.saturating_sub(2));
    let dirty_refs = dirty.iter().collect::<Vec<_>>();
    let mut toggle_pattern = |required: &[(usize, bool)]| {
        for &(index, expected) in required {
            if !expected {
                circ.x(&register[index]);
            }
        }
        let controls = required
            .iter()
            .map(|&(index, _)| &register[index])
            .collect::<Vec<_>>();
        mcx_dirty_ladder(circ, &controls, target, &dirty_refs);
        for &(index, expected) in required.iter().rev() {
            if !expected {
                circ.x(&register[index]);
            }
        }
    };
    for differing_bit in (0..width).rev() {
        if ((value >> differing_bit) & 1) != 0 {
            continue;
        }
        let mut required = Vec::with_capacity(width - differing_bit);
        required.push((differing_bit, true));
        for high_bit in differing_bit + 1..width {
            required.push((high_bit, ((value >> high_bit) & 1) != 0));
        }
        toggle_pattern(&required);
    }
    let equality = (0..width)
        .map(|index| (index, ((value >> index) & 1) != 0))
        .collect::<Vec<_>>();
    toggle_pattern(&equality);
}

fn toggle_nonzero_dirty(
    circ: &mut Circuit,
    control: &QReg,
    count: &[QReg],
    target: &QReg,
    dirty: &[QReg],
) {
    assert!(dirty.len() >= count.len().saturating_sub(1));
    circ.cx(control, target);
    for bit in count {
        circ.x(bit);
    }
    let controls = std::iter::once(control).chain(count).collect::<Vec<_>>();
    let dirty_refs = dirty.iter().collect::<Vec<_>>();
    mcx_dirty_ladder(circ, &controls, target, &dirty_refs);
    for bit in count.iter().rev() {
        circ.x(bit);
    }
}

impl<'a> Q845SwapOnlyCoefficientGuard<'a> {
    fn prepare(
        circ: &mut Circuit,
        physical_work_width: usize,
        count: &[QReg],
        l_s: &'a [QReg],
        l_r_prime: &[QReg],
        carry: &QReg,
        overflow: &QReg,
        uls_clean_lender: Option<&'a QReg>,
        relocation_l_q: Option<&'a [QReg]>,
    ) -> Self {
        assert_eq!(count.len(), l_s.len());
        assert_l_r_prime_metadata_width(count.len(), l_r_prime.len());
        assert!(physical_work_width < (1usize << count.len()));
        if q830_coefficient_counter_relocation_requested() {
            assert_eq!(
                relocation_l_q.map(|lanes| lanes.len()),
                Some(SUB820_L_Q_WIDTH),
                "relocated coefficient count requires the eight l_q lenders"
            );
        }
        if sub800_inplace_guard_address_requested() {
            let address = l_s.iter().collect::<Vec<_>>();
            let scratch = count.iter().collect::<Vec<_>>();
            affine_complement_constant_refs(
                circ,
                &address,
                physical_work_width,
                &scratch,
            );
            let source = l_r_prime
                .iter()
                .chain(count.iter().skip(l_r_prime.len()))
                .collect::<Vec<_>>();
            assert_eq!(source.len(), l_s.len());
            cuccaro_sub_mod_2n_no_overflow_refs(circ, &source, l_s, carry);
            return Self {
                address: Q845SwapOnlyAddress::InPlaceLS(l_s),
                physical_work_width,
                uls_clean_lender,
                relocation_l_q,
            };
        }

        let address = circ.alloc_qreg_bits("rs.q845-swap-only.address", count.len());
        toggle_constant(circ, &address, physical_work_width);

        for (source, destination) in l_s.iter().zip(count) {
            circ.cx(source, destination);
        }
        cuccaro_sub_mod_2n(circ, count, &address, carry, overflow);
        for (source, destination) in l_s.iter().zip(count) {
            circ.cx(source, destination);
        }
        for (source, destination) in l_r_prime.iter().zip(count) {
            circ.cx(source, destination);
        }
        cuccaro_sub_mod_2n(circ, count, &address, carry, overflow);
        for (source, destination) in l_r_prime.iter().zip(count) {
            circ.cx(source, destination);
        }

        Self {
            address: Q845SwapOnlyAddress::Allocated(address),
            physical_work_width,
            uls_clean_lender,
            relocation_l_q,
        }
    }

    fn accumulate(
        &self,
        circ: &mut Circuit,
        lower_negative: &QReg,
        active: &QReg,
        work_bit: &QReg,
        count: &[QReg],
        count_scratch: &[QReg],
        condition: &QReg,
        condition_scratch: &QReg,
    ) {
        multi_controlled_x_vchain_borrowed(
            circ,
            &[lower_negative, active, work_bit],
            condition,
            &[condition_scratch],
        );
        if q830_coefficient_counter_relocation_requested() {
            assert_eq!(count.len(), REFERENCE_LENGTH_WIDTH);
            assert_eq!(count_scratch.len(), 7);
            controlled_increment_mod_2n_hybrid_block(
                circ,
                condition,
                count,
                &count_scratch[..5],
                &count_scratch[5..],
                condition_scratch,
            );
        } else if q839_seven_plateau_lenders_requested() {
            assert_eq!(count.len(), REFERENCE_LENGTH_WIDTH);
            assert_eq!(count_scratch.len(), count.len() - 2);
            Q839_SEVEN_PLATEAU_SHORT_INCREMENTS.fetch_add(1, Ordering::Relaxed);
            controlled_increment_mod_2n_short_carry(circ, condition, count, count_scratch);
        } else {
            controlled_increment_mod_2n(circ, condition, count, count_scratch);
        }
        multi_controlled_x_vchain_borrowed(
            circ,
            &[lower_negative, active, work_bit],
            condition,
            &[condition_scratch],
        );
    }

    fn unaccumulate(
        &self,
        circ: &mut Circuit,
        lower_negative: &QReg,
        active: &QReg,
        work_bit: &QReg,
        count: &[QReg],
        count_scratch: &[QReg],
        condition: &QReg,
        condition_scratch: &QReg,
    ) {
        multi_controlled_x_vchain_borrowed(
            circ,
            &[lower_negative, active, work_bit],
            condition,
            &[condition_scratch],
        );
        if q830_coefficient_counter_relocation_requested() {
            assert_eq!(count.len(), REFERENCE_LENGTH_WIDTH);
            assert_eq!(count_scratch.len(), 7);
            controlled_decrement_mod_2n_hybrid_block(
                circ,
                condition,
                count,
                &count_scratch[..5],
                &count_scratch[5..],
                condition_scratch,
            );
        } else if q839_seven_plateau_lenders_requested() {
            assert_eq!(count.len(), REFERENCE_LENGTH_WIDTH);
            assert_eq!(count_scratch.len(), count.len() - 2);
            Q839_SEVEN_PLATEAU_SHORT_DECREMENTS.fetch_add(1, Ordering::Relaxed);
            controlled_decrement_mod_2n_short_carry(circ, condition, count, count_scratch);
        } else {
            controlled_decrement_mod_2n(circ, condition, count, count_scratch);
        }
        multi_controlled_x_vchain_borrowed(
            circ,
            &[lower_negative, active, work_bit],
            condition,
            &[condition_scratch],
        );
    }

    fn for_each_forward<F>(
        &self,
        circ: &mut Circuit,
        scan_width: usize,
        active: &QReg,
        range_scratch: &[QReg],
        mut body: F,
    )
    where
        F: FnMut(&mut Circuit, usize, &QReg),
    {
        use crate::point_add::trailmix_port::arith::khattar_gidney::{
            sub800_uls_production_callback_index,
            unary_iterate_direct_toggle_target_with_clean_scratch,
            unary_iterate_log_star_toggle_target_with_clean_lender,
            unary_iterate_log_star_toggle_target_with_clean_lenders,
            unary_iterate_log_star_with_clean_lender,
        };

        assert!(scan_width <= self.physical_work_width);
        circ.x(active);
        let address = self.address.lanes().iter().collect::<Vec<_>>();
        let truncated = q851_truncated_swap_only_guard_requested()
            && scan_width < self.physical_work_width;
        let n_iters = if truncated {
            scan_width + 1
        } else {
            self.physical_work_width + 1
        };
        if sub800_uls_direct_selector_requested() {
            SUB800_Q838_ULS_DIRECT_FORWARD_CALLS.fetch_add(1, Ordering::Relaxed);
            let selector_scratch = range_scratch.iter().collect::<Vec<_>>();
            unary_iterate_direct_toggle_target_with_clean_scratch(
                circ,
                &address,
                n_iters,
                &selector_scratch,
                active,
                true,
                |circ, callback_index| {
                    if let Some(index) = sub800_uls_production_callback_index(
                        self.physical_work_width,
                        scan_width,
                        truncated,
                        false,
                        callback_index,
                    ) {
                        body(circ, index, active);
                    }
                },
            );
        } else if sub800_uls_fused_target_requested() {
            SUB800_Q839_ULS_FUSED_FORWARD_CALLS.fetch_add(1, Ordering::Relaxed);
            if q839_seven_plateau_lenders_requested() {
                assert_eq!(range_scratch.len(), REFERENCE_LENGTH_WIDTH - 1);
                let lenders: Vec<&QReg> = if let Some(l_q) = self.relocation_l_q {
                    vec![&l_q[6], &l_q[7]]
                } else {
                    let cursor_lender = range_scratch
                        .last()
                        .expect("seven-plateau ULS cursor lender");
                    self.uls_clean_lender
                        .into_iter()
                        .chain(std::iter::once(cursor_lender))
                        .collect()
                };
                Q839_SEVEN_PLATEAU_ULS_FORWARD_LOANS.fetch_add(1, Ordering::Relaxed);
                let audit_callbacks =
                    record_q839_uls_callback_call(Q839UlsCallbackDirection::Forward, &lenders);
                unary_iterate_log_star_toggle_target_with_clean_lenders(
                    circ,
                    &address,
                    n_iters,
                    &lenders,
                    active,
                    true,
                    |circ, index| {
                        if index < scan_width {
                            if audit_callbacks {
                                let callback_start = begin_q839_uls_callback_slice(circ);
                                body(circ, index, active);
                                finish_q839_uls_callback_slice(
                                    circ,
                                    Q839UlsCallbackDirection::Forward,
                                    &lenders,
                                    callback_start,
                                );
                            } else {
                                body(circ, index, active);
                            }
                        }
                    },
                );
            } else {
                unary_iterate_log_star_toggle_target_with_clean_lender(
                    circ,
                    &address,
                    n_iters,
                    self.uls_clean_lender,
                    active,
                    true,
                    |circ, index| {
                        if index < scan_width {
                            body(circ, index, active);
                        }
                    },
                );
            }
        } else {
            unary_iterate_log_star_with_clean_lender(
                circ,
                &address,
                n_iters,
                self.uls_clean_lender,
                |circ, index, gate| {
                    circ.cx(gate, active);
                    if index < scan_width {
                        body(circ, index, active);
                    }
                },
            );
        }
        if truncated {
            if let Some(l_q) = self.relocation_l_q {
                toggle_register_geq_constant_dirty(
                    circ,
                    self.address.lanes(),
                    scan_width + 1,
                    active,
                    l_q,
                );
            } else {
                toggle_register_geq_constant_vchain(
                    circ,
                    self.address.lanes(),
                    scan_width + 1,
                    active,
                    range_scratch,
                );
            }
        }
    }

    fn for_each_reverse<F>(
        &self,
        circ: &mut Circuit,
        scan_width: usize,
        active: &QReg,
        range_scratch: &[QReg],
        constant_scratch: &[&QReg],
        constant_carry: &QReg,
        mut body: F,
    )
    where
        F: FnMut(&mut Circuit, usize, &QReg),
    {
        use crate::point_add::trailmix_port::arith::khattar_gidney::{
            sub800_uls_production_callback_index,
            unary_iterate_direct_toggle_target_with_clean_scratch,
            unary_iterate_log_star_toggle_target_with_clean_lender,
            unary_iterate_log_star_toggle_target_with_clean_lenders,
            unary_iterate_log_star_with_clean_lender,
        };

        assert!(scan_width <= self.physical_work_width);
        assert_eq!(constant_scratch.len(), self.address.lanes().len());
        let truncated = q851_truncated_swap_only_guard_requested()
            && scan_width < self.physical_work_width;
        let delta = self.physical_work_width - scan_width;
        if truncated {
            if let Some(l_q) = self.relocation_l_q {
                let address = self.address.lanes().iter().collect::<Vec<_>>();
                let dirty = l_q.iter().collect::<Vec<_>>();
                sub_constant_dirty_refs(circ, &address, delta, &dirty);
                toggle_register_geq_constant_dirty(
                    circ,
                    self.address.lanes(),
                    scan_width + 1,
                    active,
                    l_q,
                );
            } else {
                for (index, bit) in constant_scratch.iter().enumerate() {
                    if ((delta >> index) & 1) != 0 {
                        circ.x(bit);
                    }
                }
                cuccaro_sub_mod_2n_no_overflow_refs(
                    circ,
                    constant_scratch,
                    self.address.lanes(),
                    constant_carry,
                );
                for (index, bit) in constant_scratch.iter().enumerate() {
                    if ((delta >> index) & 1) != 0 {
                        circ.x(bit);
                    }
                }
                toggle_register_geq_constant_vchain(
                    circ,
                    self.address.lanes(),
                    scan_width + 1,
                    active,
                    range_scratch,
                );
            }
        }
        let address = self.address.lanes().iter().collect::<Vec<_>>();
        let n_iters = if truncated {
            scan_width + 1
        } else {
            self.physical_work_width + 1
        };
        if sub800_uls_direct_selector_requested() {
            SUB800_Q838_ULS_DIRECT_REVERSE_CALLS.fetch_add(1, Ordering::Relaxed);
            let selector_scratch = range_scratch.iter().collect::<Vec<_>>();
            unary_iterate_direct_toggle_target_with_clean_scratch(
                circ,
                &address,
                n_iters,
                &selector_scratch,
                active,
                false,
                |circ, reverse_index| {
                    let index = if truncated {
                        scan_width - reverse_index
                    } else {
                        self.physical_work_width - reverse_index
                    };
                    if index < scan_width {
                        body(circ, index, active);
                    }
                },
            );
        } else if sub800_uls_fused_target_requested() {
            SUB800_Q839_ULS_FUSED_REVERSE_CALLS.fetch_add(1, Ordering::Relaxed);
            if q839_seven_plateau_lenders_requested() {
                assert_eq!(range_scratch.len(), REFERENCE_LENGTH_WIDTH - 1);
                let lenders: Vec<&QReg> = if let Some(l_q) = self.relocation_l_q {
                    vec![&l_q[6], &l_q[7]]
                } else {
                    let cursor_lender = range_scratch
                        .last()
                        .expect("seven-plateau ULS cursor lender");
                    self.uls_clean_lender
                        .into_iter()
                        .chain(std::iter::once(cursor_lender))
                        .collect()
                };
                Q839_SEVEN_PLATEAU_ULS_REVERSE_LOANS.fetch_add(1, Ordering::Relaxed);
                let audit_callbacks =
                    record_q839_uls_callback_call(Q839UlsCallbackDirection::Reverse, &lenders);
                unary_iterate_log_star_toggle_target_with_clean_lenders(
                    circ,
                    &address,
                    n_iters,
                    &lenders,
                    active,
                    false,
                    |circ, reverse_index| {
                        if let Some(index) = sub800_uls_production_callback_index(
                            self.physical_work_width,
                            scan_width,
                            truncated,
                            true,
                            reverse_index,
                        ) {
                            if audit_callbacks {
                                let callback_start = begin_q839_uls_callback_slice(circ);
                                body(circ, index, active);
                                finish_q839_uls_callback_slice(
                                    circ,
                                    Q839UlsCallbackDirection::Reverse,
                                    &lenders,
                                    callback_start,
                                );
                            } else {
                                body(circ, index, active);
                            }
                        }
                    },
                );
            } else {
                unary_iterate_log_star_toggle_target_with_clean_lender(
                    circ,
                    &address,
                    n_iters,
                    self.uls_clean_lender,
                    active,
                    false,
                    |circ, reverse_index| {
                        if let Some(index) = sub800_uls_production_callback_index(
                            self.physical_work_width,
                            scan_width,
                            truncated,
                            true,
                            reverse_index,
                        ) {
                            body(circ, index, active);
                        }
                    },
                );
            }
        } else {
            unary_iterate_log_star_with_clean_lender(
                circ,
                &address,
                n_iters,
                self.uls_clean_lender,
                |circ, reverse_index, gate| {
                    let index = if truncated {
                        scan_width - reverse_index
                    } else {
                        self.physical_work_width - reverse_index
                    };
                    if index < scan_width {
                        body(circ, index, active);
                    }
                    circ.cx(gate, active);
                },
            );
        }
        circ.x(active);
        if truncated {
            if let Some(l_q) = self.relocation_l_q {
                let address = self.address.lanes().iter().collect::<Vec<_>>();
                let dirty = l_q.iter().collect::<Vec<_>>();
                add_constant_dirty_refs(circ, &address, delta, &dirty);
            } else {
                for (index, bit) in constant_scratch.iter().enumerate() {
                    if ((delta >> index) & 1) != 0 {
                        circ.x(bit);
                    }
                }
                cuccaro_add_mod_2n_no_overflow_refs(
                    circ,
                    constant_scratch,
                    self.address.lanes(),
                    constant_carry,
                );
                for (index, bit) in constant_scratch.iter().enumerate() {
                    if ((delta >> index) & 1) != 0 {
                        circ.x(bit);
                    }
                }
            }
        }
    }

    fn toggle_nonzero(
        &self,
        circ: &mut Circuit,
        control: &QReg,
        count: &[QReg],
        target: &QReg,
        scratch: &[QReg],
    ) {
        if let Some(l_q) = self.relocation_l_q {
            toggle_nonzero_dirty(circ, control, count, target, l_q);
            return;
        }
        assert!(scratch.len() >= count.len().saturating_sub(1));
        circ.cx(control, target);
        for bit in count {
            circ.x(bit);
        }
        let controls: Vec<&QReg> = std::iter::once(control).chain(count).collect();
        let scratch_refs: Vec<&QReg> = scratch.iter().collect();
        multi_controlled_x_vchain_borrowed(circ, &controls, target, &scratch_refs);
        for bit in count {
            circ.x(bit);
        }
    }

    fn prepare_reverse_boundary(
        &self,
        circ: &mut Circuit,
        source: &[&QReg],
        l_s: &[QReg],
        l_r_prime: &[QReg],
        carry: &QReg,
        overflow: &QReg,
    ) {
        assert_eq!(source.len(), self.address.lanes().len());
        if self.address.is_in_place() {
            let address = self.address.lanes().iter().collect::<Vec<_>>();
            if let Some(l_q) = self.relocation_l_q {
                let dirty = l_q.iter().take(6).collect::<Vec<_>>();
                affine_complement_constant_dirty_refs(
                    circ,
                    &address,
                    self.physical_work_width,
                    &dirty,
                );
            } else {
                affine_complement_constant_refs(
                    circ,
                    &address,
                    self.physical_work_width,
                    source,
                );
            }
            return;
        }
        for (bit, destination) in l_r_prime.iter().zip(source) {
            circ.cx(bit, destination);
        }
        cuccaro_add_mod_2n_refs(circ, source, self.address.lanes(), carry, overflow);
        for (bit, destination) in l_r_prime.iter().zip(source) {
            circ.cx(bit, destination);
        }
        for (bit, destination) in l_s.iter().zip(source) {
            circ.cx(bit, destination);
        }
        cuccaro_add_mod_2n_refs(circ, source, self.address.lanes(), carry, overflow);
        for (bit, destination) in l_s.iter().zip(source) {
            circ.cx(bit, destination);
        }
        toggle_constant(circ, self.address.lanes(), self.physical_work_width);

        for (bit, destination) in l_s.iter().zip(source) {
            circ.cx(bit, destination);
        }
        cuccaro_add_mod_2n_refs(circ, source, self.address.lanes(), carry, overflow);
        for (bit, destination) in l_s.iter().zip(source) {
            circ.cx(bit, destination);
        }
        for (bit, destination) in l_r_prime.iter().zip(source) {
            circ.cx(bit, destination);
        }
        cuccaro_add_mod_2n_refs(circ, source, self.address.lanes(), carry, overflow);
        for (bit, destination) in l_r_prime.iter().zip(source) {
            circ.cx(bit, destination);
        }
    }

    fn finish(
        self,
        circ: &mut Circuit,
        count: &[QReg],
        l_s: &[QReg],
        l_r_prime: &[QReg],
        carry: &QReg,
        overflow: &QReg,
    ) {
        match self.address {
            Q845SwapOnlyAddress::Allocated(address) => {
                for (bit, destination) in l_r_prime.iter().zip(count) {
                    circ.cx(bit, destination);
                }
                cuccaro_sub_mod_2n(circ, count, &address, carry, overflow);
                for (bit, destination) in l_r_prime.iter().zip(count) {
                    circ.cx(bit, destination);
                }
                for (bit, destination) in l_s.iter().zip(count) {
                    circ.cx(bit, destination);
                }
                cuccaro_sub_mod_2n(circ, count, &address, carry, overflow);
                for (bit, destination) in l_s.iter().zip(count) {
                    circ.cx(bit, destination);
                }
                free_clean(circ, address);
            }
            Q845SwapOnlyAddress::InPlaceLS(address) => {
                assert!(address
                    .iter()
                    .zip(l_s)
                    .all(|(left, right)| left.id() == right.id()));
                let source = l_r_prime
                    .iter()
                    .chain(count.iter().skip(l_r_prime.len()))
                    .collect::<Vec<_>>();
                assert_eq!(source.len(), address.len());
                cuccaro_sub_mod_2n_no_overflow_refs(circ, &source, address, carry);
            }
        }
    }
}

#[derive(Clone, Copy)]
enum Sub800GuardProofStage {
    Prepare,
    ReverseBoundary,
    FullRoundtrip,
}

struct Sub800GuardProofHarness {
    builder: B,
    l_s_ids: Vec<u32>,
    l_r_prime_ids: Vec<u32>,
    scratch_ids: Vec<u32>,
}

fn build_sub800_guard_proof_harness(stage: Sub800GuardProofStage) -> Sub800GuardProofHarness {
    let mut circ = Circuit::new();
    let l_s = circ.alloc_qreg_bits("sub800.guard-proof.l-s", REFERENCE_LENGTH_WIDTH);
    let l_r_prime =
        circ.alloc_qreg_bits("sub800.guard-proof.l-r-prime", REFERENCE_R_LENGTH_WIDTH);
    let count = circ.alloc_qreg_bits("sub800.guard-proof.count", REFERENCE_LENGTH_WIDTH);
    let carry = circ.alloc_qreg("sub800.guard-proof.carry");
    let overflow = circ.alloc_qreg("sub800.guard-proof.overflow");
    let reverse_source =
        circ.alloc_qreg_bits("sub800.guard-proof.reverse-source", REFERENCE_LENGTH_WIDTH);

    {
        let guard = Q845SwapOnlyCoefficientGuard::prepare(
            &mut circ,
            259,
            &count,
            &l_s,
            &l_r_prime,
            &carry,
            &overflow,
            None,
            None,
        );
        if matches!(
            stage,
            Sub800GuardProofStage::ReverseBoundary | Sub800GuardProofStage::FullRoundtrip
        ) {
            let source = reverse_source.iter().collect::<Vec<_>>();
            guard.prepare_reverse_boundary(
                &mut circ,
                &source,
                &l_s,
                &l_r_prime,
                &carry,
                &overflow,
            );
        }
        if matches!(stage, Sub800GuardProofStage::FullRoundtrip) {
            guard.finish(
                &mut circ,
                &count,
                &l_s,
                &l_r_prime,
                &carry,
                &overflow,
            );
        }
    }

    let l_s_ids = l_s.iter().map(QReg::id).collect::<Vec<_>>();
    let l_r_prime_ids = l_r_prime.iter().map(QReg::id).collect::<Vec<_>>();
    let scratch_ids = count
        .iter()
        .chain(std::iter::once(&carry))
        .chain(std::iter::once(&overflow))
        .chain(&reverse_source)
        .map(QReg::id)
        .collect::<Vec<_>>();
    Sub800GuardProofHarness {
        builder: circ.into_builder(),
        l_s_ids,
        l_r_prime_ids,
        scratch_ids,
    }
}

fn sub800_guard_expected(stage: Sub800GuardProofStage, l_s: usize, l_r_prime: usize) -> usize {
    match stage {
        Sub800GuardProofStage::Prepare => 259usize.wrapping_sub(l_s).wrapping_sub(l_r_prime) & 511,
        Sub800GuardProofStage::ReverseBoundary => l_s.wrapping_add(l_r_prime) & 511,
        Sub800GuardProofStage::FullRoundtrip => l_s,
    }
}

fn verify_sub800_guard_harness(
    stage: Sub800GuardProofStage,
    harness: &Sub800GuardProofHarness,
) -> (usize, usize, usize, usize, usize, usize) {
    use crate::circuit::QubitId;
    use crate::sim::Simulator;
    use sha3::{
        digest::{ExtendableOutput, Update},
        Shake128,
    };

    let mut basis_states_checked = 0usize;
    let mut coordinate_checks = 0usize;
    let mut inverse_pair_checks = 0usize;
    let mut scratch_clean_checks = 0usize;
    let mut phase_clean_checks = 0usize;
    let mut ancilla_clean_checks = 0usize;
    let state_count = 512usize * 256usize;
    for batch_start in (0..state_count).step_by(64) {
        let mut seed = Shake128::default();
        seed.update(b"sub800-inplace-guard-address-v1");
        seed.update(&[stage as u8]);
        seed.update(&(batch_start as u64).to_le_bytes());
        let mut xof = seed.finalize_xof();
        let mut simulator = Simulator::new(
            harness.builder.next_qubit as usize,
            harness.builder.next_bit as usize,
            &mut xof,
        );
        for shot in 0..64 {
            let state = batch_start + shot;
            let l_s = state & 511;
            let l_r_prime = (state >> 9) & 255;
            for (bit, id) in harness.l_s_ids.iter().enumerate() {
                if ((l_s >> bit) & 1) != 0 {
                    *simulator.qubit_mut(QubitId(u64::from(*id))) |= 1u64 << shot;
                }
            }
            for (bit, id) in harness.l_r_prime_ids.iter().enumerate() {
                if ((l_r_prime >> bit) & 1) != 0 {
                    *simulator.qubit_mut(QubitId(u64::from(*id))) |= 1u64 << shot;
                }
            }
        }

        simulator.apply_iter(harness.builder.ops.iter());
        assert_eq!(simulator.phase, 0, "sub800 guard forward phase garbage");
        for (bit, id) in harness.l_s_ids.iter().enumerate() {
            let expected = (0..64).fold(0u64, |plane, shot| {
                let state = batch_start + shot;
                let l_s = state & 511;
                let l_r_prime = (state >> 9) & 255;
                plane
                    | ((((sub800_guard_expected(stage, l_s, l_r_prime) >> bit) & 1) as u64)
                        << shot)
            });
            assert_eq!(simulator.qubit(QubitId(u64::from(*id))), expected);
        }
        for (bit, id) in harness.l_r_prime_ids.iter().enumerate() {
            let expected = (0..64).fold(0u64, |plane, shot| {
                let l_r_prime = ((batch_start + shot) >> 9) & 255;
                plane | ((((l_r_prime >> bit) & 1) as u64) << shot)
            });
            assert_eq!(simulator.qubit(QubitId(u64::from(*id))), expected);
        }
        for id in &harness.scratch_ids {
            assert_eq!(simulator.qubit(QubitId(u64::from(*id))), 0);
        }

        simulator.apply_iter(harness.builder.ops.iter().rev());
        assert_eq!(simulator.phase, 0, "sub800 guard inverse phase garbage");
        for (bit, id) in harness.l_s_ids.iter().enumerate() {
            let expected = (0..64).fold(0u64, |plane, shot| {
                let l_s = (batch_start + shot) & 511;
                plane | ((((l_s >> bit) & 1) as u64) << shot)
            });
            assert_eq!(simulator.qubit(QubitId(u64::from(*id))), expected);
        }
        for (bit, id) in harness.l_r_prime_ids.iter().enumerate() {
            let expected = (0..64).fold(0u64, |plane, shot| {
                let l_r_prime = ((batch_start + shot) >> 9) & 255;
                plane | ((((l_r_prime >> bit) & 1) as u64) << shot)
            });
            assert_eq!(simulator.qubit(QubitId(u64::from(*id))), expected);
        }
        for id in &harness.scratch_ids {
            assert_eq!(simulator.qubit(QubitId(u64::from(*id))), 0);
        }

        basis_states_checked += 64;
        coordinate_checks += 64;
        inverse_pair_checks += 64;
        scratch_clean_checks += 2 * 64;
        phase_clean_checks += 2 * 64;
        ancilla_clean_checks += 2 * 64;
    }
    (
        basis_states_checked,
        coordinate_checks,
        inverse_pair_checks,
        scratch_clean_checks,
        phase_clean_checks,
        ancilla_clean_checks,
    )
}

/// Exhaustively verify the three affine-coordinate stages used by the
/// default-off sub-800 guard-address candidate.
#[doc(hidden)]
#[must_use]
pub fn exhaustive_sub800_inplace_guard_address_check() -> Sub800InplaceGuardAddressProofReport {
    let saved = std::env::var_os(SUB800_INPLACE_GUARD_ADDRESS_FLAG);
    std::env::set_var(SUB800_INPLACE_GUARD_ADDRESS_FLAG, "1");
    let stages = [
        Sub800GuardProofStage::Prepare,
        Sub800GuardProofStage::ReverseBoundary,
        Sub800GuardProofStage::FullRoundtrip,
    ];
    let harnesses = stages.map(build_sub800_guard_proof_harness);
    let mut totals = [0usize; 6];
    for (stage, harness) in stages.into_iter().zip(&harnesses) {
        let report = verify_sub800_guard_harness(stage, harness);
        for (total, value) in totals.iter_mut().zip([
            report.0, report.1, report.2, report.3, report.4, report.5,
        ]) {
            *total += value;
        }
    }
    match saved {
        Some(value) => std::env::set_var(SUB800_INPLACE_GUARD_ADDRESS_FLAG, value),
        None => std::env::remove_var(SUB800_INPLACE_GUARD_ADDRESS_FLAG),
    }

    Sub800InplaceGuardAddressProofReport {
        stages_checked: stages.len(),
        basis_states_checked: totals[0],
        coordinate_checks: totals[1],
        inverse_pair_checks: totals[2],
        scratch_clean_checks: totals[3],
        phase_clean_checks: totals[4],
        ancilla_clean_checks: totals[5],
        allocated_address_lanes: 0,
        prepare: gate_counts(&harnesses[0].builder.ops),
        reverse_boundary: gate_counts(&harnesses[1].builder.ops),
        full_roundtrip: gate_counts(&harnesses[2].builder.ops),
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum Q851CoefficientCursorDirection {
    Decrement,
    Increment,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum Q851CoefficientCursorTraversal {
    InverseForward,
    InverseReverse,
    ForwardForward,
    ForwardReverse,
}

fn q851_coefficient_cursor_transition(
    traversal: Q851CoefficientCursorTraversal,
    index: usize,
    scan_width: usize,
) -> Option<(usize, Q851CoefficientCursorDirection)> {
    assert!(index < scan_width);
    match traversal {
        Q851CoefficientCursorTraversal::InverseForward => index
            .checked_sub(1)
            .map(|event| (event, Q851CoefficientCursorDirection::Decrement)),
        Q851CoefficientCursorTraversal::InverseReverse => (index + 1 != scan_width)
            .then_some((index, Q851CoefficientCursorDirection::Increment)),
        Q851CoefficientCursorTraversal::ForwardForward => (index + 1 != scan_width)
            .then_some((index, Q851CoefficientCursorDirection::Decrement)),
        Q851CoefficientCursorTraversal::ForwardReverse => index
            .checked_sub(1)
            .map(|event| (event, Q851CoefficientCursorDirection::Increment)),
    }
}

fn toggle_q851_fixed_sign_event(
    circ: &mut Circuit,
    signed_length: &[QReg],
    transition_index: usize,
    scratch: &[QReg],
) {
    assert_eq!(signed_length.len(), REFERENCE_LENGTH_WIDTH);
    assert!(
        scratch.len() >= REFERENCE_LENGTH_WIDTH - 3,
        "fixed-sign cursor event requires six clean v-chain lanes"
    );
    assert!(transition_index < 256);
    for (index, lane) in signed_length.iter().enumerate() {
        for other in &signed_length[index + 1..] {
            assert_ne!(lane.id(), other.id(), "coefficient cursor aliases itself");
        }
        for other in scratch {
            assert_ne!(lane.id(), other.id(), "coefficient cursor aliases scratch");
        }
    }
    for (index, lane) in scratch.iter().enumerate() {
        for other in &scratch[index + 1..] {
            assert_ne!(lane.id(), other.id(), "coefficient scratch aliases itself");
        }
    }

    if q830_dirty_fixed_sign_event_requested()
        || q830_coefficient_counter_relocation_requested()
    {
        let low = &signed_length[..REFERENCE_LENGTH_WIDTH - 1];
        for (bit, lane) in low.iter().enumerate() {
            if transition_index & (1usize << bit) == 0 {
                circ.x(lane);
            }
        }
        let controls = low.iter().collect::<Vec<_>>();
        let dirty = scratch[..REFERENCE_LENGTH_WIDTH - 3]
            .iter()
            .collect::<Vec<_>>();
        mcx_dirty_ladder(
            circ,
            &controls,
            signed_length.last().expect("nonempty signed length"),
            &dirty,
        );
        for (bit, lane) in low.iter().enumerate().rev() {
            if transition_index & (1usize << bit) == 0 {
                circ.x(lane);
            }
        }
        return;
    }

    // For 0 <= i <= 256,
    // MSB((L - i) mod 512) = MSB(L) XOR [i > (L mod 256)].
    // The transition after body i therefore toggles only when L mod 256 = i.
    let low = &signed_length[..REFERENCE_LENGTH_WIDTH - 1];
    for (bit, lane) in low.iter().enumerate() {
        if transition_index & (1usize << bit) == 0 {
            circ.x(lane);
        }
    }
    let controls = low.iter().collect::<Vec<_>>();
    multi_controlled_x_vchain(
        circ,
        &controls,
        signed_length.last().expect("nonempty signed length"),
        scratch,
    );
    for (bit, lane) in low.iter().enumerate().rev() {
        if transition_index & (1usize << bit) == 0 {
            circ.x(lane);
        }
    }
}

fn transition_q851_coefficient_cursor(
    circ: &mut Circuit,
    signed_length: &[QReg],
    scan_width: usize,
    transition_index: usize,
    scratch: &[QReg],
    direction: Q851CoefficientCursorDirection,
) {
    let fixed_sign_domain = signed_length.len() == REFERENCE_LENGTH_WIDTH && scan_width <= 257;
    if q851_fixed_sign_event_requested() && fixed_sign_domain {
        toggle_q851_fixed_sign_event(circ, signed_length, transition_index, scratch);
        return;
    }
    match direction {
        Q851CoefficientCursorDirection::Decrement => {
            decrement_mod_2n(circ, signed_length, scratch)
        }
        Q851CoefficientCursorDirection::Increment => {
            increment_mod_2n(circ, signed_length, scratch)
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn transition_q851_coefficient_cursor_at_body(
    circ: &mut Circuit,
    signed_length: &[QReg],
    scan_width: usize,
    index: usize,
    scratch: &[QReg],
    traversal: Q851CoefficientCursorTraversal,
) {
    if let Some((transition_index, direction)) =
        q851_coefficient_cursor_transition(traversal, index, scan_width)
    {
        transition_q851_coefficient_cursor(
            circ,
            signed_length,
            scan_width,
            transition_index,
            scratch,
            direction,
        );
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum Q845FusedCoefficientBody {
    AddForward,
    AddReverse,
    SubForward,
    SubReverse,
}

fn emit_q845_fused_coefficient_body(
    circ: &mut Circuit,
    body: Q845FusedCoefficientBody,
    active: &QReg,
    carry: &QReg,
    work1: &QReg,
    work2: &QReg,
    tmp: &QReg,
) {
    match body {
        Q845FusedCoefficientBody::AddForward => {
            circ.ccx(active, carry, work2);
            circ.ccx(active, carry, work1);
            multi_controlled_x_vchain(
                circ,
                &[active, work1, work2],
                carry,
                std::slice::from_ref(tmp),
            );
        }
        Q845FusedCoefficientBody::AddReverse => {
            multi_controlled_x_vchain(
                circ,
                &[active, work1, work2],
                carry,
                std::slice::from_ref(tmp),
            );
            circ.ccx(active, carry, work1);
            circ.ccx(active, work1, work2);
        }
        Q845FusedCoefficientBody::SubForward => {
            circ.ccx(active, work1, work2);
            circ.ccx(active, carry, work1);
            multi_controlled_x_vchain(
                circ,
                &[active, work1, work2],
                carry,
                std::slice::from_ref(tmp),
            );
        }
        Q845FusedCoefficientBody::SubReverse => {
            multi_controlled_x_vchain(
                circ,
                &[active, work1, work2],
                carry,
                std::slice::from_ref(tmp),
            );
            circ.ccx(active, carry, work1);
            circ.ccx(active, carry, work2);
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn toggle_output_coefficient_enable_q845_inline_underflow(
    circ: &mut Circuit,
    phase1: &QReg,
    phase2: &QReg,
    sign: &QReg,
    carry: &QReg,
    above_guard: &QReg,
    enable: &QReg,
    scratch: &[&QReg],
) {
    assert_eq!(scratch.len(), 3);
    // phase1 AND (phase2 OR sign_out OR (carry AND !above_guard)), split
    // into disjoint terms. The final five-control term borrows add_only as its
    // third v-chain lane only at cuts where add_only is proved zero.
    circ.ccx(phase1, phase2, enable);
    circ.x(phase2);
    multi_controlled_x_vchain_borrowed(circ, &[phase1, phase2, sign], enable, scratch);
    circ.x(sign);
    circ.x(above_guard);
    multi_controlled_x_vchain_borrowed(
        circ,
        &[phase1, phase2, sign, carry, above_guard],
        enable,
        scratch,
    );
    circ.x(above_guard);
    circ.x(sign);
    circ.x(phase2);
}

#[allow(clippy::too_many_arguments)]
fn coefficient_fused_data_and_sign_q845_swap_only(
    circ: &mut Circuit,
    phase1: &QReg,
    phase2: &QReg,
    sign: &QReg,
    work1: &[QReg],
    work2: &[QReg],
    physical_work_width: usize,
    l_t: &[QReg],
    l_t_prime: &[QReg],
    l_s: &[QReg],
    l_r_prime: &[QReg],
    enable: &QReg,
    above_guard: &QReg,
    add_only: &QReg,
    chain: &[QReg],
    inverse: bool,
    uls_clean_lender: Option<&QReg>,
    relocation_l_q: Option<&[QReg]>,
) {
    assert_eq!(work1.len(), work2.len());
    assert!(!work1.is_empty());
    assert_eq!(l_t.len(), l_s.len());
    assert_l_r_prime_metadata_width(l_t.len(), l_r_prime.len());
    assert_eq!(chain.len(), 2);
    let relocate_count = q830_coefficient_counter_relocation_requested();
    if relocate_count {
        assert_eq!(l_t.len(), REFERENCE_LENGTH_WIDTH);
        assert_eq!(l_t_prime.len(), 1);
        assert_eq!(relocation_l_q.map(|lanes| lanes.len()), Some(SUB820_L_Q_WIDTH));
        assert!(uls_clean_lender.is_some());
    } else {
        assert_eq!(l_t.len(), l_t_prime.len());
        assert!(relocation_l_q.is_none());
    }

    let relocation_cursor_scratch = relocate_count.then(|| {
        circ.alloc_qreg_bits(
            "rs.q845-swap-only.cursor-scratch",
            l_t.len().saturating_sub(1),
        )
    });
    let count_storage = if let Some(cursor_scratch) = relocation_cursor_scratch.as_ref() {
        cursor_scratch
            .iter()
            .map(QReg::borrowed_alias)
            .chain(l_t_prime.iter().map(QReg::borrowed_alias))
            .collect::<Vec<_>>()
    } else {
        l_t_prime.iter().map(QReg::borrowed_alias).collect::<Vec<_>>()
    };

    let guard = Q845SwapOnlyCoefficientGuard::prepare(
        circ,
        physical_work_width,
        &count_storage,
        l_s,
        l_r_prime,
        &chain[0],
        &chain[1],
        uls_clean_lender,
        relocation_l_q,
    );
    let cursor_scratch = relocation_cursor_scratch.unwrap_or_else(|| {
        circ.alloc_qreg_bits(
            "rs.q845-swap-only.cursor-scratch",
            l_t.len().saturating_sub(1),
        )
    });
    let callback_cursor_scratch = if let Some(l_q) = relocation_l_q {
        l_q[..6]
            .iter()
            .map(QReg::borrowed_alias)
            .collect::<Vec<_>>()
    } else if q839_seven_plateau_lenders_requested() {
        assert_eq!(cursor_scratch.len(), REFERENCE_LENGTH_WIDTH - 1);
        cursor_scratch[..cursor_scratch.len() - 1]
            .iter()
            .map(QReg::borrowed_alias)
            .collect::<Vec<_>>()
    } else {
        cursor_scratch.iter().map(QReg::borrowed_alias).collect::<Vec<_>>()
    };
    let callback_counter_scratch = if let Some(l_q) = relocation_l_q {
        l_q[..6]
            .iter()
            .map(QReg::borrowed_alias)
            .chain(uls_clean_lender.into_iter().map(QReg::borrowed_alias))
            .collect::<Vec<_>>()
    } else {
        callback_cursor_scratch
            .iter()
            .map(QReg::borrowed_alias)
            .collect::<Vec<_>>()
    };
    assert_eq!(count_storage.len(), REFERENCE_LENGTH_WIDTH);
    let carry = circ.alloc_qreg("rs.q845-swap-only.carry");
    let coefficient_active = &chain[0];
    let tmp = &chain[1];
    let output_scratch = [&chain[0], &chain[1], add_only];
    let lower_negative = l_t.last().expect("nonempty coefficient cursor");
    let bracket = production_coefficient_nonnegative_bracket();

    if inverse {
        guard.for_each_forward(
            circ,
            work1.len(),
            above_guard,
            &cursor_scratch,
            |circ, index, guard_active| {
            transition_q851_coefficient_cursor_at_body(
                circ,
                l_t,
                work1.len(),
                index,
                &callback_cursor_scratch,
                Q851CoefficientCursorTraversal::InverseForward,
            );
            guard.accumulate(
                circ,
                lower_negative,
                guard_active,
                &work2[index],
                &count_storage,
                &callback_counter_scratch,
                coefficient_active,
                tmp,
            );
            begin_coefficient_nonnegative_sign_bracket(
                circ,
                phase1,
                lower_negative,
                coefficient_active,
                bracket,
            );
            emit_q845_fused_coefficient_body(
                circ,
                Q845FusedCoefficientBody::SubForward,
                coefficient_active,
                &carry,
                &work1[index],
                &work2[index],
                tmp,
            );
            end_coefficient_nonnegative_sign_bracket(
                circ,
                phase1,
                lower_negative,
                coefficient_active,
                bracket,
            );
            },
        );
        guard.toggle_nonzero(circ, phase1, &count_storage, above_guard, &cursor_scratch);

        toggle_output_coefficient_enable_q845_inline_underflow(
            circ,
            phase1,
            phase2,
            sign,
            &carry,
            above_guard,
            enable,
            &output_scratch,
        );
        circ.x(enable);
        circ.ccx(phase1, enable, add_only);
        circ.x(enable);
        circ.x(above_guard);
        circ.ccx(&carry, above_guard, sign);
        circ.x(above_guard);
        circ.cx(phase1, sign);
        guard.toggle_nonzero(circ, enable, &count_storage, above_guard, &cursor_scratch);

        let reverse_source: Vec<&QReg> = if relocate_count {
            count_storage.iter().collect()
        } else {
            cursor_scratch
                .iter()
                .chain(std::iter::once(coefficient_active))
                .collect()
        };
        guard.prepare_reverse_boundary(
            circ,
            &reverse_source,
            l_s,
            l_r_prime,
            tmp,
            above_guard,
        );
        guard.for_each_reverse(
            circ,
            work1.len(),
            above_guard,
            &cursor_scratch,
            &reverse_source,
            tmp,
            |circ, index, guard_active| {
            transition_q851_coefficient_cursor_at_body(
                circ,
                l_t,
                work1.len(),
                index,
                &callback_cursor_scratch,
                Q851CoefficientCursorTraversal::InverseReverse,
            );
            begin_coefficient_nonnegative_sign_bracket(
                circ,
                add_only,
                lower_negative,
                coefficient_active,
                bracket,
            );
            emit_q845_fused_coefficient_body(
                circ,
                Q845FusedCoefficientBody::SubReverse,
                coefficient_active,
                &carry,
                &work1[index],
                &work2[index],
                tmp,
            );
            end_coefficient_nonnegative_sign_bracket(
                circ,
                add_only,
                lower_negative,
                coefficient_active,
                bracket,
            );
            begin_coefficient_nonnegative_sign_bracket(
                circ,
                enable,
                lower_negative,
                coefficient_active,
                bracket,
            );
            emit_q845_fused_coefficient_body(
                circ,
                Q845FusedCoefficientBody::AddReverse,
                coefficient_active,
                &carry,
                &work1[index],
                &work2[index],
                tmp,
            );
            end_coefficient_nonnegative_sign_bracket(
                circ,
                enable,
                lower_negative,
                coefficient_active,
                bracket,
            );
            guard.unaccumulate(
                circ,
                lower_negative,
                guard_active,
                &work2[index],
                &count_storage,
                &callback_counter_scratch,
                coefficient_active,
                tmp,
            );
            },
        );
        circ.x(enable);
        circ.ccx(phase1, enable, add_only);
        circ.x(enable);
        toggle_initial_coefficient_enable(circ, phase1, phase2, sign, enable, chain);
    } else {
        toggle_initial_coefficient_enable(circ, phase1, phase2, sign, enable, chain);
        circ.x(enable);
        circ.ccx(phase1, enable, add_only);
        circ.x(enable);
        guard.for_each_forward(
            circ,
            work1.len(),
            above_guard,
            &cursor_scratch,
            |circ, index, guard_active| {
            guard.accumulate(
                circ,
                lower_negative,
                guard_active,
                &work2[index],
                &count_storage,
                &callback_counter_scratch,
                coefficient_active,
                tmp,
            );
            begin_coefficient_nonnegative_sign_bracket(
                circ,
                enable,
                lower_negative,
                coefficient_active,
                bracket,
            );
            emit_q845_fused_coefficient_body(
                circ,
                Q845FusedCoefficientBody::SubForward,
                coefficient_active,
                &carry,
                &work1[index],
                &work2[index],
                tmp,
            );
            end_coefficient_nonnegative_sign_bracket(
                circ,
                enable,
                lower_negative,
                coefficient_active,
                bracket,
            );
            begin_coefficient_nonnegative_sign_bracket(
                circ,
                add_only,
                lower_negative,
                coefficient_active,
                bracket,
            );
            emit_q845_fused_coefficient_body(
                circ,
                Q845FusedCoefficientBody::AddForward,
                coefficient_active,
                &carry,
                &work1[index],
                &work2[index],
                tmp,
            );
            end_coefficient_nonnegative_sign_bracket(
                circ,
                add_only,
                lower_negative,
                coefficient_active,
                bracket,
            );
            transition_q851_coefficient_cursor_at_body(
                circ,
                l_t,
                work1.len(),
                index,
                &callback_cursor_scratch,
                Q851CoefficientCursorTraversal::ForwardForward,
            );
            },
        );
        guard.toggle_nonzero(circ, enable, &count_storage, above_guard, &cursor_scratch);

        circ.cx(phase1, sign);
        circ.x(above_guard);
        circ.ccx(&carry, above_guard, sign);
        circ.x(above_guard);
        circ.x(enable);
        circ.ccx(phase1, enable, add_only);
        circ.x(enable);
        toggle_output_coefficient_enable_q845_inline_underflow(
            circ,
            phase1,
            phase2,
            sign,
            &carry,
            above_guard,
            enable,
            &output_scratch,
        );
        guard.toggle_nonzero(circ, phase1, &count_storage, above_guard, &cursor_scratch);

        let reverse_source: Vec<&QReg> = if relocate_count {
            count_storage.iter().collect()
        } else {
            cursor_scratch
                .iter()
                .chain(std::iter::once(coefficient_active))
                .collect()
        };
        guard.prepare_reverse_boundary(
            circ,
            &reverse_source,
            l_s,
            l_r_prime,
            tmp,
            above_guard,
        );
        guard.for_each_reverse(
            circ,
            work1.len(),
            above_guard,
            &cursor_scratch,
            &reverse_source,
            tmp,
            |circ, index, guard_active| {
            begin_coefficient_nonnegative_sign_bracket(
                circ,
                phase1,
                lower_negative,
                coefficient_active,
                bracket,
            );
            emit_q845_fused_coefficient_body(
                circ,
                Q845FusedCoefficientBody::AddReverse,
                coefficient_active,
                &carry,
                &work1[index],
                &work2[index],
                tmp,
            );
            end_coefficient_nonnegative_sign_bracket(
                circ,
                phase1,
                lower_negative,
                coefficient_active,
                bracket,
            );
            guard.unaccumulate(
                circ,
                lower_negative,
                guard_active,
                &work2[index],
                &count_storage,
                &callback_counter_scratch,
                coefficient_active,
                tmp,
            );
            transition_q851_coefficient_cursor_at_body(
                circ,
                l_t,
                work1.len(),
                index,
                &callback_cursor_scratch,
                Q851CoefficientCursorTraversal::ForwardReverse,
            );
            },
        );
    }

    guard.finish(
        circ,
        &count_storage,
        l_s,
        l_r_prime,
        coefficient_active,
        tmp,
    );
    circ.zero_and_free(carry);
    free_clean(circ, cursor_scratch);
}

#[allow(clippy::too_many_arguments)]
fn coefficient_fused_data_and_sign_q845(
    circ: &mut Circuit,
    phase1: &QReg,
    phase2: &QReg,
    sign: &QReg,
    work1: &[QReg],
    work2: &[QReg],
    l_t: &[QReg],
    enable: &QReg,
    above_guard: &QReg,
    add_only: &QReg,
    chain: &[QReg],
    inverse: bool,
) {
    assert_eq!(work1.len(), work2.len());
    assert!(!work1.is_empty());
    assert!(!l_t.is_empty());
    assert_eq!(chain.len(), 2);
    assert_ne!(chain[0].id(), chain[1].id());
    let cursor_scratch = circ.alloc_qreg_bits(
        "rs.q845-coeff-fused.cursor-scratch",
        l_t.len().saturating_sub(1),
    );
    let carry = circ.alloc_qreg("rs.q845-coeff-fused.carry");
    let active = &chain[0];
    let tmp = &chain[1];
    let output_scratch = [&chain[0], &chain[1], add_only];
    let bracket = production_coefficient_nonnegative_bracket();

    if inverse {
        for index in 0..work1.len() {
            if index != 0 {
                decrement_mod_2n(circ, l_t, &cursor_scratch);
            }
            begin_coefficient_nonnegative_bracket(circ, phase1, l_t, active, bracket);
            emit_q845_fused_coefficient_body(
                circ,
                Q845FusedCoefficientBody::SubForward,
                active,
                &carry,
                &work1[index],
                &work2[index],
                tmp,
            );
            end_coefficient_nonnegative_bracket(circ, phase1, l_t, active, bracket);
        }

        toggle_output_coefficient_enable_q845_inline_underflow(
            circ,
            phase1,
            phase2,
            sign,
            &carry,
            above_guard,
            enable,
            &output_scratch,
        );
        circ.x(enable);
        circ.ccx(phase1, enable, add_only);
        circ.x(enable);
        circ.x(above_guard);
        circ.ccx(&carry, above_guard, sign);
        circ.x(above_guard);
        circ.cx(phase1, sign);

        for index in (0..work1.len()).rev() {
            if index + 1 != work1.len() {
                increment_mod_2n(circ, l_t, &cursor_scratch);
            }
            begin_coefficient_nonnegative_bracket(circ, add_only, l_t, active, bracket);
            emit_q845_fused_coefficient_body(
                circ,
                Q845FusedCoefficientBody::SubReverse,
                active,
                &carry,
                &work1[index],
                &work2[index],
                tmp,
            );
            end_coefficient_nonnegative_bracket(circ, add_only, l_t, active, bracket);
            begin_coefficient_nonnegative_bracket(circ, enable, l_t, active, bracket);
            emit_q845_fused_coefficient_body(
                circ,
                Q845FusedCoefficientBody::AddReverse,
                active,
                &carry,
                &work1[index],
                &work2[index],
                tmp,
            );
            end_coefficient_nonnegative_bracket(circ, enable, l_t, active, bracket);
        }
        circ.x(enable);
        circ.ccx(phase1, enable, add_only);
        circ.x(enable);
        toggle_initial_coefficient_enable(circ, phase1, phase2, sign, enable, chain);
    } else {
        circ.x(enable);
        circ.ccx(phase1, enable, add_only);
        circ.x(enable);
        for index in 0..work1.len() {
            begin_coefficient_nonnegative_bracket(circ, enable, l_t, active, bracket);
            emit_q845_fused_coefficient_body(
                circ,
                Q845FusedCoefficientBody::SubForward,
                active,
                &carry,
                &work1[index],
                &work2[index],
                tmp,
            );
            end_coefficient_nonnegative_bracket(circ, enable, l_t, active, bracket);
            begin_coefficient_nonnegative_bracket(circ, add_only, l_t, active, bracket);
            emit_q845_fused_coefficient_body(
                circ,
                Q845FusedCoefficientBody::AddForward,
                active,
                &carry,
                &work1[index],
                &work2[index],
                tmp,
            );
            end_coefficient_nonnegative_bracket(circ, add_only, l_t, active, bracket);
            if index + 1 != work1.len() {
                decrement_mod_2n(circ, l_t, &cursor_scratch);
            }
        }

        circ.cx(phase1, sign);
        circ.x(above_guard);
        circ.ccx(&carry, above_guard, sign);
        circ.x(above_guard);
        circ.x(enable);
        circ.ccx(phase1, enable, add_only);
        circ.x(enable);
        toggle_output_coefficient_enable_q845_inline_underflow(
            circ,
            phase1,
            phase2,
            sign,
            &carry,
            above_guard,
            enable,
            &output_scratch,
        );

        for index in (0..work1.len()).rev() {
            begin_coefficient_nonnegative_bracket(circ, phase1, l_t, active, bracket);
            emit_q845_fused_coefficient_body(
                circ,
                Q845FusedCoefficientBody::AddReverse,
                active,
                &carry,
                &work1[index],
                &work2[index],
                tmp,
            );
            end_coefficient_nonnegative_bracket(circ, phase1, l_t, active, bracket);
            if index != 0 {
                increment_mod_2n(circ, l_t, &cursor_scratch);
            }
        }
    }

    circ.zero_and_free(carry);
    free_clean(circ, cursor_scratch);
}

#[allow(clippy::too_many_arguments)]
fn coefficient_phase_block_fused_q845(
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
    uls_clean_lender: Option<&QReg>,
    relocation_l_q: Option<&[QReg]>,
) {
    let enable = circ.alloc_qreg("rs.q845-coeff-fused.enable");
    let above_guard = circ.alloc_qreg("rs.q845-coeff-fused.above-guard");
    let add_only = circ.alloc_qreg("rs.q845-coeff-fused.add-only");
    let chain = circ.alloc_qreg_bits("rs.q845-coeff-fused.chain", 2);
    let guard_scratch = [&chain[0], &chain[1], &add_only];

    if q845_swap_only_t_prime_length_requested() {
        coefficient_fused_data_and_sign_q845_swap_only(
            circ,
            phase1,
            phase2,
            sign,
            work1,
            work2,
            full_work2.len(),
            l_t,
            l_t_prime,
            l_s,
            l_r_prime,
            &enable,
            &above_guard,
            &add_only,
            &chain,
            inverse,
            uls_clean_lender,
            relocation_l_q,
        );
        free_clean(circ, chain);
        circ.zero_and_free(add_only);
        circ.zero_and_free(above_guard);
        circ.zero_and_free(enable);
        return;
    }

    if inverse {
        toggle_q845_coefficient_length_above_guarded_boundary(
            circ,
            phase1,
            l_t_prime,
            l_t,
            l_s,
            &above_guard,
            &guard_scratch,
        );
        controlled_xor_raw_t_prime_bit_length(
            circ,
            phase1,
            l_s,
            l_r_prime,
            full_work2,
            l_t_prime,
            &chain,
        );
        coefficient_fused_data_and_sign_q845(
            circ,
            phase1,
            phase2,
            sign,
            work1,
            work2,
            l_t,
            &enable,
            &above_guard,
            &add_only,
            &chain,
            true,
        );
        controlled_xor_raw_t_prime_bit_length(
            circ,
            phase1,
            l_s,
            l_r_prime,
            full_work2,
            l_t_prime,
            &chain,
        );
        toggle_initial_coefficient_enable(circ, phase1, phase2, sign, &enable, &chain);
        toggle_q845_coefficient_length_above_guarded_boundary(
            circ,
            &enable,
            l_t_prime,
            l_t,
            l_s,
            &above_guard,
            &guard_scratch,
        );
        toggle_initial_coefficient_enable(circ, phase1, phase2, sign, &enable, &chain);
    } else {
        toggle_initial_coefficient_enable(circ, phase1, phase2, sign, &enable, &chain);
        toggle_q845_coefficient_length_above_guarded_boundary(
            circ,
            &enable,
            l_t_prime,
            l_t,
            l_s,
            &above_guard,
            &guard_scratch,
        );
        controlled_xor_raw_t_prime_bit_length(
            circ,
            phase1,
            l_s,
            l_r_prime,
            full_work2,
            l_t_prime,
            &chain,
        );
        coefficient_fused_data_and_sign_q845(
            circ,
            phase1,
            phase2,
            sign,
            work1,
            work2,
            l_t,
            &enable,
            &above_guard,
            &add_only,
            &chain,
            false,
        );
        controlled_xor_raw_t_prime_bit_length(
            circ,
            phase1,
            l_s,
            l_r_prime,
            full_work2,
            l_t_prime,
            &chain,
        );
        toggle_q845_coefficient_length_above_guarded_boundary(
            circ,
            phase1,
            l_t_prime,
            l_t,
            l_s,
            &above_guard,
            &guard_scratch,
        );
    }

    free_clean(circ, chain);
    circ.zero_and_free(add_only);
    circ.zero_and_free(above_guard);
    circ.zero_and_free(enable);
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
    coefficient_phase_block_with_uls_clean_lender(
        circ,
        phase1,
        phase2,
        sign,
        work1,
        work2,
        full_work2,
        l_t,
        l_t_prime,
        l_s,
        l_r_prime,
        inverse,
        None,
        None,
    );
}

#[allow(clippy::too_many_arguments)]
fn coefficient_phase_block_with_uls_clean_lender(
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
    uls_clean_lender: Option<&QReg>,
    relocation_l_q: Option<&[QReg]>,
) {
    assert!(
        q845_swap_only_coefficient_dependencies_satisfied(),
        "Q845 swap-only t-prime lifecycle requires Q845 coefficient fusion"
    );
    if q845_lifetime_coefficient_fusion_requested() {
        coefficient_phase_block_fused_q845(
            circ, phase1, phase2, sign, work1, work2, full_work2, l_t, l_t_prime, l_s, l_r_prime,
            inverse, uls_clean_lender, relocation_l_q,
        );
        return;
    }

    assert!(
        uls_clean_lender.is_none(),
        "ULS clean lending requires the Q845 fused coefficient route"
    );
    assert!(
        relocation_l_q.is_none(),
        "counter relocation requires the Q845 fused coefficient route"
    );

    let enable = circ.alloc_qreg("rs.coeff-block.enable");
    let less_than = circ.alloc_qreg("rs.coeff-block.less-than");
    let add_only = circ.alloc_qreg("rs.coeff-block.add-only");
    let borrowed_comparator = borrowed_coefficient_comparator_requested();
    let borrowed_less_than = coefficient_less_than_lane_reuse_requested();
    let chain = circ.alloc_qreg_bits("rs.coeff-block.chain", 2);
    // The v-chain restores both lanes after every enable update. `add_only` is
    // zero before the first comparison and is exactly uncomputed before the
    // second comparison in either direction. Each boundary comparator restores
    // these lanes before the less-than carry loop borrows them.
    let compare_scratch = if borrowed_comparator || borrowed_less_than {
        vec![&chain[0], &chain[1], &add_only]
    } else {
        Vec::new()
    };
    let coefficient_add_lenders = [&chain[0], &chain[1]];
    // The raw-bitlength loan is deliberately narrower: add_only is live as
    // its control, so only the two restored chain lanes are passed onward.
    assert_eq!(l_t.len(), l_t_prime.len());

    if inverse {
        toggle_coefficient_less_than(
            circ,
            phase1,
            work1,
            work2,
            l_t,
            l_s,
            l_t_prime,
            &less_than,
            &compare_scratch,
        );
        toggle_output_coefficient_enable(circ, phase1, phase2, sign, &less_than, &enable, &chain);

        circ.x(&enable);
        circ.ccx(phase1, &enable, &add_only);
        circ.x(&enable);
        with_bit_length_callsite(
            "p.bitlen.rs.coeff-inverse.pre-add.deposit",
            "p.bitlen.rs.coeff-inverse.pre-add.erase",
            || {
                controlled_xor_raw_t_prime_bit_length(
                    circ, &add_only, l_s, l_r_prime, full_work2, l_t_prime, &chain,
                );
            },
        );
        coefficient_add_data_only(
            circ,
            &add_only,
            work1,
            work2,
            l_t,
            true,
            &coefficient_add_lenders,
        );
        with_bit_length_callsite(
            "p.bitlen.rs.coeff-inverse.post-add.deposit",
            "p.bitlen.rs.coeff-inverse.post-add.erase",
            || {
                controlled_xor_raw_t_prime_bit_length(
                    circ, &add_only, l_s, l_r_prime, full_work2, l_t_prime, &chain,
                );
            },
        );
        circ.x(&enable);
        circ.ccx(phase1, &enable, &add_only);
        circ.x(&enable);

        circ.cx(&less_than, sign);
        circ.cx(phase1, sign);
        toggle_coefficient_less_than(
            circ,
            &enable,
            work1,
            work2,
            l_t,
            l_s,
            l_t_prime,
            &less_than,
            &compare_scratch,
        );
        toggle_initial_coefficient_enable(circ, phase1, phase2, sign, &enable, &chain);
    } else {
        toggle_initial_coefficient_enable(circ, phase1, phase2, sign, &enable, &chain);
        toggle_coefficient_less_than(
            circ,
            &enable,
            work1,
            work2,
            l_t,
            l_s,
            l_t_prime,
            &less_than,
            &compare_scratch,
        );
        circ.cx(phase1, sign);
        circ.cx(&less_than, sign);

        circ.x(&enable);
        circ.ccx(phase1, &enable, &add_only);
        circ.x(&enable);
        with_bit_length_callsite(
            "p.bitlen.rs.coeff-forward.pre-add.deposit",
            "p.bitlen.rs.coeff-forward.pre-add.erase",
            || {
                controlled_xor_raw_t_prime_bit_length(
                    circ, &add_only, l_s, l_r_prime, full_work2, l_t_prime, &chain,
                );
            },
        );
        coefficient_add_data_only(
            circ,
            &add_only,
            work1,
            work2,
            l_t,
            false,
            &coefficient_add_lenders,
        );
        with_bit_length_callsite(
            "p.bitlen.rs.coeff-forward.post-add.deposit",
            "p.bitlen.rs.coeff-forward.post-add.erase",
            || {
                controlled_xor_raw_t_prime_bit_length(
                    circ, &add_only, l_s, l_r_prime, full_work2, l_t_prime, &chain,
                );
            },
        );
        circ.x(&enable);
        circ.ccx(phase1, &enable, &add_only);
        circ.x(&enable);

        toggle_output_coefficient_enable(circ, phase1, phase2, sign, &less_than, &enable, &chain);
        toggle_coefficient_less_than(
            circ,
            phase1,
            work1,
            work2,
            l_t,
            l_s,
            l_t_prime,
            &less_than,
            &compare_scratch,
        );
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
    assert_l_r_prime_metadata_width(l_t.len(), l_r_prime.len());
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
        remainder_scratch_width(length_width, l_r_prime.len()),
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
        conditional_work_and_length_swap_under_zero_predicate(
            circ,
            zero_q,
            zero_s,
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
            if rotated_bitlen_scratch_reuse_requested() {
                chain
            } else {
                &[]
            },
            &[],
            false,
            PromisedLqSwapRoute::Configured,
        );
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
    register_shared_scheduled_step_with_preserved_dy_top(
        circ,
        step,
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
        None,
    );
}

#[allow(clippy::too_many_arguments)]
fn register_shared_scheduled_step_with_preserved_dy_top(
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
    preserved_dy_top: Option<&QReg>,
) {
    assert_eq!(work1.len(), 259);
    assert_eq!(work2.len(), 259);
    assert_eq!(l_t.len(), REFERENCE_LENGTH_WIDTH);
    assert_eq!(
        l_t_prime.len(),
        if q830_coefficient_counter_relocation_requested() {
            1
        } else {
            REFERENCE_LENGTH_WIDTH
        }
    );
    assert_eq!(l_q.len(), SUB820_L_Q_WIDTH);
    assert_eq!(l_s.len(), REFERENCE_LENGTH_WIDTH);
    assert_eq!(l_r_prime.len(), production_l_r_prime_width());
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
        remainder_scratch_width(REFERENCE_LENGTH_WIDTH, l_r_prime.len()),
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
    coefficient_phase_block_with_uls_clean_lender(
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
        preserved_dy_top.filter(|_| sub800_uls_clean_lender_requested()),
        q830_coefficient_counter_relocation_requested().then_some(l_q),
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
        let preserved_dy_top_scratch: Vec<&QReg> = preserved_dy_top
            .filter(|_| preserved_dy_top_prefix_loan_requested())
            .into_iter()
            .collect();
        conditional_work_and_length_swap_under_zero_predicate(
            circ,
            zero_q,
            zero_s,
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
            if rotated_bitlen_scratch_reuse_requested() {
                chain
            } else {
                &[]
            },
            &preserved_dy_top_scratch,
            false,
            PromisedLqSwapRoute::Configured,
        );
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
    register_shared_scheduled_step_inverse_with_preserved_dy_top(
        circ,
        step,
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
        None,
    );
}

#[allow(clippy::too_many_arguments)]
fn register_shared_scheduled_step_inverse_with_preserved_dy_top(
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
    preserved_dy_top: Option<&QReg>,
) {
    assert_eq!(work1.len(), 259);
    assert_eq!(work2.len(), 259);
    assert_eq!(l_t.len(), REFERENCE_LENGTH_WIDTH);
    assert_eq!(
        l_t_prime.len(),
        if q830_coefficient_counter_relocation_requested() {
            1
        } else {
            REFERENCE_LENGTH_WIDTH
        }
    );
    assert_eq!(l_q.len(), SUB820_L_Q_WIDTH);
    assert_eq!(l_s.len(), REFERENCE_LENGTH_WIDTH);
    assert_eq!(l_r_prime.len(), production_l_r_prime_width());
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
        let preserved_dy_top_scratch: Vec<&QReg> = preserved_dy_top
            .filter(|_| preserved_dy_top_prefix_loan_requested())
            .into_iter()
            .collect();
        conditional_work_and_length_swap_under_zero_predicate(
            circ,
            zero_q,
            zero_s,
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
            if rotated_bitlen_scratch_reuse_requested() {
                chain
            } else {
                &[]
            },
            &preserved_dy_top_scratch,
            true,
            PromisedLqSwapRoute::Configured,
        );
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
    coefficient_phase_block_with_uls_clean_lender(
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
        preserved_dy_top.filter(|_| sub800_uls_clean_lender_requested()),
        q830_coefficient_counter_relocation_requested().then_some(l_q),
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
        remainder_scratch_width(REFERENCE_LENGTH_WIDTH, l_r_prime.len()),
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

    let l_r_prime = circ.alloc_qreg_bits("rs.divider.l-r-prime", production_l_r_prime_width());
    let reflected_source_width = if mixed_width_l_r_prime_requested() {
        255
    } else {
        256
    };
    let source: Vec<&QReg> = dx.iter().take(reflected_source_width).collect();
    bit_length_lean(circ, &source, &l_r_prime, false);

    dx.push(circ.alloc_qreg("rs.divider.work2-pad0"));
    dx.push(circ.alloc_qreg("rs.divider.work2-pad1"));
    dx.reverse();
    let work2 = dx;

    let work1 = circ.alloc_qreg_bits("rs.divider.work1", REGISTER_SHARED_WORK_WIDTH);
    toggle_initial_work1(circ, &work1);
    let l_t = circ.alloc_qreg_bits("rs.divider.l-t", REFERENCE_LENGTH_WIDTH);
    let l_t_prime_width = if q830_coefficient_counter_relocation_requested() {
        assert!(
            q830_direct_swap_metadata_requested(),
            "counter relocation requires direct swap metadata"
        );
        1
    } else {
        REFERENCE_LENGTH_WIDTH
    };
    let l_t_prime = circ.alloc_qreg_bits("rs.divider.l-t-prime", l_t_prime_width);
    let l_q = circ.alloc_qreg_bits("rs.divider.l-q", SUB820_L_Q_WIDTH);
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

fn register_shared_forward(
    circ: &mut Circuit,
    core: &RegisterSharedCore,
    preserved_dy_top: Option<&QReg>,
) {
    for step in 1..=REFERENCE_STEPS {
        register_shared_scheduled_step_with_preserved_dy_top(
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
            preserved_dy_top,
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
    let l_q = circ.alloc_qreg_bits("rs.divider.l-q.rebuilt", SUB820_L_Q_WIDTH);
    let l_r_prime =
        circ.alloc_qreg_bits("rs.divider.l-r-prime.rebuilt", production_l_r_prime_width());
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

fn register_shared_reverse(
    circ: &mut Circuit,
    core: &RegisterSharedCore,
    preserved_dy_top: Option<&QReg>,
) {
    for step in (1..=REFERENCE_STEPS).rev() {
        register_shared_scheduled_step_inverse_with_preserved_dy_top(
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
            preserved_dy_top,
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

    let reflected_source_width = if mixed_width_l_r_prime_requested() {
        255
    } else {
        256
    };
    let source: Vec<&QReg> = core.work2.iter().take(reflected_source_width).collect();
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
    // The canonical top lane is zero throughout the EEA schedule and is not
    // otherwise touched until terminal multiplication. Each local borrower
    // restores it before returning.
    let preserved_dy_top = &dy[REGISTER_SHARED_FIELD_WIDTH - 1];
    register_shared_forward(circ, &core, Some(preserved_dy_top));
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
    let dy_ghosts: Vec<_> = dy.iter().map(|lane| circ.hmr_ghost(lane)).collect();
    free_clean(circ, dy);

    let core = register_shared_rebuild_terminal(circ, terminal);
    // Keep the already-clean canonical top lane through the rebuilt reverse
    // schedule. It hosts the omitted l_q high bit at swap scans and replaces
    // the third ULS ancilla, so this lifetime extension does not recreate the
    // removed persistent lane at the global peak.
    register_shared_reverse(circ, &core, Some(&lambda_top));
    let dx = register_shared_finish(circ, core);

    lambda.push(lambda_top);
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
    // `dy` is the canonical product restored by divide-forward. Its clean top
    // lane is restored by every borrower before this schedule returns to the
    // terminal multiplication below, and remains available to the rebuilt
    // inverse schedule after that multiplication is undone.
    let preserved_dy_top = &dy[REGISTER_SHARED_FIELD_WIDTH - 1];
    register_shared_forward(circ, &core, Some(preserved_dy_top));
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
    register_shared_reverse(circ, &core, Some(preserved_dy_top));
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
        conditional_work_and_length_swap_under_zero_predicate(
            circ,
            zero_q,
            zero_s,
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
            if rotated_bitlen_scratch_reuse_requested() {
                chain
            } else {
                &[]
            },
            &[],
            true,
            PromisedLqSwapRoute::Configured,
        );
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
        remainder_scratch_width(length_width, l_r_prime.len()),
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

fn build_coefficient_boundary_comparator(
    length_width: usize,
    borrowed: bool,
    controlled_roundtrip: bool,
) -> B {
    assert!(length_width > 0);
    let previous_flag = std::env::var_os("LOWQ_REUSE_COEFFICIENT_COMPARATOR_SCRATCH");
    if borrowed {
        std::env::set_var("LOWQ_REUSE_COEFFICIENT_COMPARATOR_SCRATCH", "1");
    } else {
        std::env::remove_var("LOWQ_REUSE_COEFFICIENT_COMPARATOR_SCRATCH");
    }

    let mut circ = Circuit::new();
    let control = controlled_roundtrip.then(|| circ.alloc_qreg("rs.coeff-proof.control"));
    let target_length = circ.alloc_qreg_bits("rs.coeff-proof.target-length", length_width);
    let l_t = circ.alloc_qreg_bits("rs.coeff-proof.l-t", length_width);
    let l_s = circ.alloc_qreg_bits("rs.coeff-proof.l-s", length_width);
    let target = circ.alloc_qreg("rs.coeff-proof.target");
    let scratch = if borrowed {
        circ.alloc_qreg_bits("rs.coeff-proof.borrowed-scratch", 3)
    } else {
        Vec::new()
    };
    let scratch_refs: Vec<&QReg> = scratch.iter().collect();

    if let Some(control) = control.as_ref() {
        let predicate = circ.alloc_qreg("rs.coeff-proof.predicate");
        toggle_coefficient_length_above_boundary(
            &mut circ,
            &target_length,
            &l_t,
            &l_s,
            &predicate,
            &scratch_refs,
        );
        circ.ccx(control, &predicate, &target);
        toggle_coefficient_length_above_boundary(
            &mut circ,
            &target_length,
            &l_t,
            &l_s,
            &predicate,
            &scratch_refs,
        );
        circ.zero_and_free(predicate);
    } else {
        toggle_coefficient_length_above_boundary(
            &mut circ,
            &target_length,
            &l_t,
            &l_s,
            &target,
            &scratch_refs,
        );
    }
    free_clean(&mut circ, scratch);

    match previous_flag {
        Some(value) => std::env::set_var("LOWQ_REUSE_COEFFICIENT_COMPARATOR_SCRATCH", value),
        None => std::env::remove_var("LOWQ_REUSE_COEFFICIENT_COMPARATOR_SCRATCH"),
    }
    circ.into_builder()
}

#[must_use]
pub fn exhaustive_borrowed_coefficient_comparator_check(
) -> ReferenceBorrowedCoefficientComparatorProofReport {
    let mut basis_states_checked = 0usize;
    let mut control_off_identity_checks = 0usize;
    let mut equality_boundary_checks = 0usize;
    let mut subtraction_underflow_checks = 0usize;
    let mut addition_overflow_checks = 0usize;
    let mut oracle_equivalence_checks = 0usize;
    let mut inverse_pair_checks = 0usize;
    let mut operand_restore_checks = 0usize;
    let mut ancilla_clean_checks = 0usize;

    for length_width in 1..=3 {
        let baseline = build_coefficient_boundary_comparator(length_width, false, true);
        let borrowed = build_coefficient_boundary_comparator(length_width, true, true);
        let mask = (1u64 << length_width) - 1;
        let target_offset = 1 + 3 * length_width;
        let data_width = target_offset + 1;
        let preserved_mask = (1u64 << target_offset) - 1;

        for input in 0..(1u64 << data_width) {
            let control = input & 1 != 0;
            let target_length = (input >> 1) & mask;
            let l_t = (input >> (1 + length_width)) & mask;
            let l_s = (input >> (1 + 2 * length_width)) & mask;
            let full_sum = l_t + l_s;
            let boundary = full_sum & mask;
            let predicate = target_length > boundary;
            let expected = input ^ (u64::from(control && predicate) << target_offset);

            let baseline_output = apply_scalar(&baseline.ops, input);
            let borrowed_output = apply_scalar(&borrowed.ops, input);
            assert_eq!(baseline_output, expected);
            assert_eq!(borrowed_output, expected);
            assert_eq!(borrowed_output, baseline_output);
            assert_eq!(baseline_output >> data_width, 0);
            assert_eq!(borrowed_output >> data_width, 0);
            assert_eq!(baseline_output & preserved_mask, input & preserved_mask);
            assert_eq!(borrowed_output & preserved_mask, input & preserved_mask);
            assert_eq!(apply_scalar(&baseline.ops, baseline_output), input);
            assert_eq!(apply_scalar(&borrowed.ops, borrowed_output), input);

            basis_states_checked += 1;
            oracle_equivalence_checks += 1;
            inverse_pair_checks += 2;
            operand_restore_checks += 2;
            ancilla_clean_checks += 2;
            if !control {
                assert_eq!(borrowed_output, input);
                control_off_identity_checks += 1;
            }
            if target_length == boundary {
                equality_boundary_checks += 1;
            }
            if predicate {
                subtraction_underflow_checks += 1;
            }
            if full_sum > mask {
                addition_overflow_checks += 1;
            }
        }
    }

    let baseline_reference9_builder =
        build_coefficient_boundary_comparator(REFERENCE_LENGTH_WIDTH, false, false);
    let borrowed_reference9_builder =
        build_coefficient_boundary_comparator(REFERENCE_LENGTH_WIDTH, true, false);
    let baseline_reference9 = measurement_classical_gate_counts(&baseline_reference9_builder.ops);
    let borrowed_reference9 = measurement_classical_gate_counts(&borrowed_reference9_builder.ops);
    let baseline_reference9_active_qubits = baseline_reference9_builder.active_qubits as usize;
    let baseline_reference9_peak_qubits = baseline_reference9_builder.peak_qubits as usize;
    let borrowed_reference9_active_qubits = borrowed_reference9_builder.active_qubits as usize;
    let borrowed_reference9_peak_qubits = borrowed_reference9_builder.peak_qubits as usize;
    let baseline_reference9_temporary_qubits =
        baseline_reference9_peak_qubits - baseline_reference9_active_qubits;
    let borrowed_reference9_temporary_qubits =
        borrowed_reference9_peak_qubits - borrowed_reference9_active_qubits;

    assert_eq!(
        baseline_reference9_active_qubits,
        3 * REFERENCE_LENGTH_WIDTH + 1
    );
    assert_eq!(
        borrowed_reference9_active_qubits,
        baseline_reference9_active_qubits
    );
    assert_eq!(baseline_reference9_temporary_qubits, 33);
    assert_eq!(borrowed_reference9_temporary_qubits, 13);
    assert!(borrowed_reference9.ccx < baseline_reference9.ccx);
    assert!(borrowed_reference9_peak_qubits < baseline_reference9_peak_qubits);

    ReferenceBorrowedCoefficientComparatorProofReport {
        widths_checked: 3,
        basis_states_checked,
        control_off_identity_checks,
        equality_boundary_checks,
        subtraction_underflow_checks,
        addition_overflow_checks,
        oracle_equivalence_checks,
        inverse_pair_checks,
        operand_restore_checks,
        ancilla_clean_checks,
        baseline_reference9,
        borrowed_reference9,
        baseline_reference9_active_qubits,
        baseline_reference9_peak_qubits,
        baseline_reference9_temporary_qubits,
        borrowed_reference9_active_qubits,
        borrowed_reference9_peak_qubits,
        borrowed_reference9_temporary_qubits,
        borrowed_caller_lanes: 3,
        borrowed_reference9_fresh_qubits: REFERENCE_LENGTH_WIDTH + 1,
        reference9_toffoli_reduction: baseline_reference9.ccx - borrowed_reference9.ccx,
        reference9_standalone_peak_reduction: baseline_reference9_peak_qubits
            - borrowed_reference9_peak_qubits,
        production_incremental_caller_qubits: 0,
        reference9_production_local_peak_reduction: baseline_reference9_temporary_qubits
            - (REFERENCE_LENGTH_WIDTH + 1),
    }
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct CoefficientRawBitLengthLoanLocalResources {
    pub active_qubits: usize,
    pub peak_qubits: usize,
    pub temporary_qubits: usize,
    pub emitted_ops: usize,
    pub emitted_toffoli: usize,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct CoefficientRawBitLengthLoanProofReport {
    pub configurations_checked: usize,
    pub zero_modes_checked: usize,
    pub directions_checked: usize,
    pub basis_states_checked: usize,
    pub baseline_equivalence_checks: usize,
    pub simulator_equivalence_checks: usize,
    pub phase_clean_checks: usize,
    pub inverse_pair_checks: usize,
    pub control_off_checks: usize,
    pub add_only_one_checks: usize,
    pub comparator_boundary_checks: usize,
    pub borrowed_lane_clean_checks: usize,
    pub ancilla_clean_checks: usize,
    pub no_overflow_widths_checked: usize,
    pub no_overflow_basis_states_checked: usize,
    pub no_overflow_inverse_pair_checks: usize,
    pub default_stream_identity_checks: usize,
    pub precondition_rejections: usize,
    pub legacy_baseline_zero_carry_allocations: usize,
    pub legacy_loan_zero_carry_allocations: usize,
    pub fused_loan_zero_flag_allocations: usize,
    pub fused_loan_zero_carry_allocations: usize,
    pub fused_loan_zero_prefix_allocations: usize,
    pub baseline_raw_rotation_carry_allocations: usize,
    pub baseline_raw_rotation_overflow_allocations: usize,
    pub baseline_rotated_carry_allocations: usize,
    pub baseline_rotated_overflow_allocations: usize,
    pub baseline_rotated_enabled_allocations: usize,
    pub loan_raw_rotation_carry_allocations: usize,
    pub loan_raw_rotation_overflow_allocations: usize,
    pub loan_rotated_carry_allocations: usize,
    pub loan_rotated_overflow_allocations: usize,
    pub loan_rotated_enabled_allocations: usize,
    pub baseline_local: CoefficientRawBitLengthLoanLocalResources,
    pub loan_local: CoefficientRawBitLengthLoanLocalResources,
    pub local_qubit_delta: i64,
    pub local_toffoli_delta: i64,
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct InplaceRotatedBoundaryLocalResources {
    pub active_qubits: usize,
    pub peak_qubits: usize,
    pub temporary_qubits: usize,
    pub emitted_ops: usize,
    pub emitted_toffoli: usize,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct InplaceRotatedBoundaryProofReport {
    pub configurations_checked: usize,
    pub zero_modes_checked: usize,
    pub scratch_modes_checked: usize,
    pub basis_states_checked: usize,
    pub scalar_equivalence_checks: usize,
    pub simulator_equivalence_checks: usize,
    pub phase_clean_checks: usize,
    pub inverse_pair_checks: usize,
    pub control_off_checks: usize,
    pub saturation_underflow_checks: usize,
    pub nonnegative_difference_checks: usize,
    pub preserved_input_checks: usize,
    pub borrowed_scratch_clean_checks: usize,
    pub ancilla_clean_checks: usize,
    pub arithmetic_widths_checked: usize,
    pub arithmetic_basis_states_checked: usize,
    pub arithmetic_inverse_pair_checks: usize,
    pub default_stream_identity_checks: usize,
    pub configured_stream_selection_checks: usize,
    pub baseline_boundary_qubits_allocated: usize,
    pub candidate_boundary_qubits_allocated: usize,
    pub candidate_inplace_boundary_uses: usize,
    pub baseline_local: InplaceRotatedBoundaryLocalResources,
    pub candidate_local: InplaceRotatedBoundaryLocalResources,
    pub local_qubit_delta: i64,
    pub local_toffoli_delta: i64,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct BorrowedRotatedUnderflowProofReport {
    pub configurations_checked: usize,
    pub boundary_forms_checked: usize,
    pub paired_source_modes_checked: usize,
    pub basis_states_checked: usize,
    pub scalar_equivalence_checks: usize,
    pub simulator_equivalence_checks: usize,
    pub phase_clean_checks: usize,
    pub inverse_pair_checks: usize,
    pub control_off_checks: usize,
    pub saturation_underflow_checks: usize,
    pub nonnegative_difference_checks: usize,
    pub preserved_input_checks: usize,
    pub borrowed_scratch_clean_checks: usize,
    pub ancilla_clean_checks: usize,
    pub same_boundary_baseline: InplaceRotatedBoundaryLocalResources,
    pub same_boundary_candidate: InplaceRotatedBoundaryLocalResources,
    pub mixed_boundary_baseline: InplaceRotatedBoundaryLocalResources,
    pub mixed_boundary_candidate: InplaceRotatedBoundaryLocalResources,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct SplitMixedRotatedLengthProofReport {
    pub configurations_checked: usize,
    pub paired_source_modes_checked: usize,
    pub basis_states_checked: usize,
    pub scalar_equivalence_checks: usize,
    pub simulator_equivalence_checks: usize,
    pub phase_clean_checks: usize,
    pub inverse_pair_checks: usize,
    pub control_off_checks: usize,
    pub high_indicator_checks: usize,
    pub saturation_underflow_checks: usize,
    pub nonnegative_difference_checks: usize,
    pub preserved_input_checks: usize,
    pub borrowed_scratch_clean_checks: usize,
    pub ancilla_clean_checks: usize,
    pub baseline_local: InplaceRotatedBoundaryLocalResources,
    pub candidate_local: InplaceRotatedBoundaryLocalResources,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct SplitSameRotatedLengthProofReport {
    pub configurations_checked: usize,
    pub paired_source_modes_checked: usize,
    pub basis_states_checked: usize,
    pub scalar_equivalence_checks: usize,
    pub simulator_equivalence_checks: usize,
    pub phase_clean_checks: usize,
    pub inverse_pair_checks: usize,
    pub control_off_checks: usize,
    pub high_indicator_checks: usize,
    pub saturation_underflow_checks: usize,
    pub nonnegative_difference_checks: usize,
    pub preserved_input_checks: usize,
    pub borrowed_scratch_clean_checks: usize,
    pub ancilla_clean_checks: usize,
    pub baseline_local: InplaceRotatedBoundaryLocalResources,
    pub candidate_local: InplaceRotatedBoundaryLocalResources,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct SplitTwoHighRotatedLengthProofReport {
    pub borrow_truth_table_cases: usize,
    pub bit_lengths_checked: usize,
    pub high_indicator_basis_states: usize,
    pub high_indicator_inverse_checks: usize,
    pub high_indicator_source_restore_checks: usize,
    pub high_indicator_scratch_clean_checks: usize,
    pub same_arithmetic_cases: usize,
    pub mixed_arithmetic_cases: usize,
    pub same_circuit_basis_states: usize,
    pub mixed_circuit_basis_states: usize,
    pub scalar_equivalence_checks: usize,
    pub inverse_pair_checks: usize,
    pub control_off_checks: usize,
    pub phase_clean_checks: usize,
    pub scratch_restore_checks: usize,
    pub same_boundary_baseline: InplaceRotatedBoundaryLocalResources,
    pub same_boundary_candidate: InplaceRotatedBoundaryLocalResources,
    pub mixed_boundary_baseline: InplaceRotatedBoundaryLocalResources,
    pub mixed_boundary_candidate: InplaceRotatedBoundaryLocalResources,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct SplitThreeHighRotatedLengthProofReport {
    pub borrow_truth_table_cases: usize,
    pub bit_lengths_checked: usize,
    pub regression_192_255_cases: usize,
    pub complement_modes_checked: usize,
    pub initial_high_states_checked: usize,
    pub high_indicator_basis_states: usize,
    pub high_indicator_inverse_checks: usize,
    pub high_indicator_source_restore_checks: usize,
    pub high_indicator_scratch_clean_checks: usize,
    pub same_arithmetic_cases: usize,
    pub mixed_arithmetic_cases: usize,
    pub same_circuit_basis_states: usize,
    pub mixed_circuit_basis_states: usize,
    pub scalar_equivalence_checks: usize,
    pub inverse_pair_checks: usize,
    pub control_off_checks: usize,
    pub phase_clean_checks: usize,
    pub scratch_restore_checks: usize,
    pub dirty_lender_restore_checks: usize,
    pub route_precedence_checks: usize,
    pub same_boundary_two_high: InplaceRotatedBoundaryLocalResources,
    pub same_boundary_candidate: InplaceRotatedBoundaryLocalResources,
    pub mixed_boundary_two_high: InplaceRotatedBoundaryLocalResources,
    pub mixed_boundary_candidate: InplaceRotatedBoundaryLocalResources,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct SplitFourHighRotatedLengthProofReport {
    pub borrow_truth_table_cases: usize,
    pub bit_lengths_checked: usize,
    pub regression_224_255_cases: usize,
    pub complement_modes_checked: usize,
    pub source_classes_checked: usize,
    pub initial_high_states_checked: usize,
    pub high_indicator_basis_states: usize,
    pub high_indicator_inverse_checks: usize,
    pub high_indicator_source_restore_checks: usize,
    pub high_indicator_scratch_clean_checks: usize,
    pub same_arithmetic_cases: usize,
    pub mixed_arithmetic_cases: usize,
    pub same_circuit_basis_states: usize,
    pub mixed_circuit_basis_states: usize,
    pub production_same_circuit_cases: usize,
    pub production_mixed_circuit_cases: usize,
    pub production_inverse_checks: usize,
    pub production_scratch_restore_checks: usize,
    pub production_dirty_lender_restore_checks: usize,
    pub scalar_equivalence_checks: usize,
    pub inverse_pair_checks: usize,
    pub control_off_checks: usize,
    pub phase_clean_checks: usize,
    pub scratch_restore_checks: usize,
    pub dirty_lender_restore_checks: usize,
    pub route_precedence_checks: usize,
    pub same_boundary_three_high: InplaceRotatedBoundaryLocalResources,
    pub same_boundary_candidate: InplaceRotatedBoundaryLocalResources,
    pub mixed_boundary_three_high: InplaceRotatedBoundaryLocalResources,
    pub mixed_boundary_candidate: InplaceRotatedBoundaryLocalResources,
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct PairedBitLengthSourceComplementLocalResources {
    pub active_qubits: usize,
    pub peak_qubits: usize,
    pub temporary_qubits: usize,
    pub emitted_ops: usize,
    pub emitted_x: usize,
    pub emitted_toffoli: usize,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct PairedBitLengthSourceComplementProofReport {
    pub source_widths_checked: usize,
    pub maximum_source_width: usize,
    pub boundary_forms_checked: usize,
    pub scratch_modes_checked: usize,
    pub boundary_routes_checked: usize,
    pub configurations_checked: usize,
    pub basis_states_checked: usize,
    pub oracle_checks: usize,
    pub scalar_equivalence_checks: usize,
    pub simulator_equivalence_checks: usize,
    pub phase_clean_checks: usize,
    pub inverse_pair_checks: usize,
    pub control_off_checks: usize,
    pub source_restore_checks: usize,
    pub boundary_restore_checks: usize,
    pub borrowed_scratch_clean_checks: usize,
    pub ancilla_clean_checks: usize,
    pub default_stream_identity_checks: usize,
    pub non_x_kind_identity_checks: usize,
    pub local_source_width: usize,
    pub local_output_width: usize,
    pub local_boundary_width: usize,
    pub local_baseline: PairedBitLengthSourceComplementLocalResources,
    pub local_optimized: PairedBitLengthSourceComplementLocalResources,
    pub local_qubit_delta: i64,
    pub local_ops_delta: i64,
    pub local_x_delta: i64,
    pub local_toffoli_delta: i64,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct CoefficientNonnegativeXCancelProofReport {
    pub cursor_widths_checked: usize,
    pub body_kinds_checked: usize,
    pub configurations_checked: usize,
    pub basis_states_checked: usize,
    pub scalar_equivalence_checks: usize,
    pub simulator_equivalence_checks: usize,
    pub phase_clean_checks: usize,
    pub inverse_pair_checks: usize,
    pub cursor_restore_checks: usize,
    pub active_restore_checks: usize,
    pub scratch_clean_checks: usize,
    pub control_off_identity_checks: usize,
    pub default_stream_identity_checks: usize,
    pub non_x_kind_identity_checks: usize,
    pub local_ops_delta: i64,
    pub local_x_delta: i64,
    pub local_toffoli_delta: i64,
    pub scheduled_steps: usize,
    pub active_coefficient_positions: usize,
    pub bracket_pairs_per_phase_block_position: usize,
    pub coefficient_phase_blocks_per_point_add: usize,
    pub removed_x_per_position: usize,
    pub total_removed_x: usize,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct Q845LifetimeCoefficientFusionProofReport {
    pub guard_widths_checked: usize,
    pub guard_basis_states_checked: usize,
    pub guard_scratch_clean_checks: usize,
    pub guard_carry_out_cases_checked: usize,
    pub packed_widths_checked: usize,
    pub packed_basis_states_checked: usize,
    pub packed_rotation_mapping_checks: usize,
    pub packed_wrapped_residue_checks: usize,
    pub packed_source_guard_clean_checks: usize,
    pub packed_underflow_equivalence_checks: usize,
    pub packed_add_headroom_checks: usize,
    pub minimum_spill_width: usize,
    pub work_widths_checked: usize,
    pub active_width_cases_checked: usize,
    pub promised_basis_states_checked: usize,
    pub oracle_transition_checks: usize,
    pub inverse_pair_checks: usize,
    pub control_off_identity_checks: usize,
    pub underflow_checks: usize,
    pub above_guard_checks: usize,
    pub add_headroom_checks: usize,
    pub scratch_clean_checks: usize,
    pub cursor_restore_checks: usize,
    pub default_off_stream_identity_checks: usize,
    pub dispatch_stream_identity_checks: usize,
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct FusedPrefixScratchLoanLocalResources {
    pub active_qubits: usize,
    pub peak_qubits: usize,
    pub temporary_qubits: usize,
    pub emitted_ops: usize,
    pub emitted_toffoli: usize,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct FusedPrefixScratchLoanProofReport {
    pub source_widths_checked: usize,
    pub maximum_source_width: usize,
    pub accumulator_width: usize,
    pub lender_modes_checked: usize,
    pub directions_checked: usize,
    pub basis_states_checked: usize,
    pub scalar_equivalence_checks: usize,
    pub simulator_equivalence_checks: usize,
    pub phase_clean_checks: usize,
    pub inverse_pair_checks: usize,
    pub source_restore_checks: usize,
    pub lender_clean_checks: usize,
    pub ancilla_clean_checks: usize,
    pub default_stream_identity_checks: usize,
    pub kg_reverse_composition_checks: usize,
    pub kg_reverse_simulator_checks: usize,
    pub kg_reverse_phase_clean_checks: usize,
    pub kg_reverse_inverse_pair_checks: usize,
    pub kg_reverse_scratch_clean_checks: usize,
    pub kg_reverse_changed_streams: usize,
    pub alias_rejections: usize,
    pub three_lender_owned_lanes: usize,
    pub three_lender_borrowed_lanes: usize,
    pub seven_lender_owned_lanes: usize,
    pub seven_lender_borrowed_lanes: usize,
    pub eight_lender_owned_lanes: usize,
    pub eight_lender_borrowed_lanes: usize,
    pub nine_lender_owned_lanes: usize,
    pub nine_lender_borrowed_lanes: usize,
    pub baseline_local_trace: FusedPrefixScratchLoanAllocationTrace,
    pub candidate_local_trace: FusedPrefixScratchLoanAllocationTrace,
    pub baseline_local: FusedPrefixScratchLoanLocalResources,
    pub candidate_local: FusedPrefixScratchLoanLocalResources,
    pub local_qubit_delta: i64,
    pub local_ops_delta: i64,
    pub local_toffoli_delta: i64,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum RawBitLengthZeroMode {
    Legacy,
    Fused,
}

struct RawBitLengthProofEnvironment {
    saved: Vec<(&'static str, Option<std::ffi::OsString>)>,
}

impl RawBitLengthProofEnvironment {
    fn capture() -> Self {
        const NAMES: [&str; 30] = [
            COEFFICIENT_RAW_BITLEN_LOAN_FLAG,
            INPLACE_ROTATED_BITLEN_BOUNDARY_FLAG,
            FUSED_PREFIX_SCRATCH_LOAN_FLAG,
            PRESERVED_DY_TOP_PREFIX_LOAN_FLAG,
            MIXED_WIDTH_L_R_PRIME_FLAG,
            PROMISED_LQ_SWAP_BORROW_FLAG,
            SPLIT_COEFFICIENT_ROTATION_LIFETIME_FLAG,
            COEFFICIENT_LESS_THAN_LANE_REUSE_FLAG,
            super::shrunken_pz_state_machine::CALLER_SCRATCH_KG_REVERSE_DECREMENT_FLAG,
            CLEAN_CHAIN_COEFFICIENT_ADD_LENDER_FLAG,
            "LOWQ_REUSE_COEFFICIENT_COMPARATOR_SCRATCH",
            "LOWQ_DIRECT_PREFIX_BITLEN",
            "LOWQ_REUSE_ZERO_CARRIES_FOR_PREFIX",
            "LOWQ_REUSE_ZERO_CARRIES_FOR_FULL_PREFIX_SCRATCH",
            "LOWQ_FUSED_ZERO_PREFIX_BITLEN",
            "LOWQ_DIRECT_PREFIX_DIRTY_UPDATE",
            "LOWQ_DIRECT_PREFIX_NO_FLAG",
            "LOWQ_RS_DIRTY_ZERO_CORRECTION",
            PAIRED_BITLEN_SOURCE_COMPLEMENT_FLAG,
            PROMISED_SWAP_SUPPORT_LIFETIME_FUSION_FLAG,
            Q845_SWAP_ONLY_T_PRIME_LENGTH_FLAG,
            SUB800_MIXED_BOUNDARY_SCRATCH_EXTENSION_FLAG,
            SUB800_BORROWED_ROTATED_UNDERFLOW_FLAG,
            SUB800_SPLIT_MIXED_ROTATED_LENGTH_FLAG,
            SUB800_SPLIT_SAME_ROTATED_LENGTH_FLAG,
            SUB800_SPLIT_TWO_HIGH_ROTATED_LENGTH_FLAG,
            SUB800_SPLIT_THREE_HIGH_ROTATED_LENGTH_FLAG,
            SUB800_SPLIT_FOUR_HIGH_ROTATED_LENGTH_FLAG,
            SUB800_ULS_FUSED_TARGET_FLAG,
            SUB800_ULS_DIRECT_SELECTOR_FLAG,
        ];
        Self {
            saved: NAMES
                .into_iter()
                .map(|name| (name, std::env::var_os(name)))
                .collect(),
        }
    }
}

impl Drop for RawBitLengthProofEnvironment {
    fn drop(&mut self) {
        for (name, value) in self.saved.drain(..) {
            match value {
                Some(value) => std::env::set_var(name, value),
                None => std::env::remove_var(name),
            }
        }
    }
}

fn configure_raw_bit_length_loan_proof(mode: RawBitLengthZeroMode, loan: bool) {
    std::env::set_var("LOWQ_REUSE_COEFFICIENT_COMPARATOR_SCRATCH", "1");
    std::env::set_var("LOWQ_DIRECT_PREFIX_BITLEN", "1");
    std::env::set_var("LOWQ_REUSE_ZERO_CARRIES_FOR_PREFIX", "1");
    for name in [
        "LOWQ_DIRECT_PREFIX_DIRTY_UPDATE",
        "LOWQ_DIRECT_PREFIX_NO_FLAG",
        "LOWQ_RS_DIRTY_ZERO_CORRECTION",
    ] {
        std::env::remove_var(name);
    }
    match mode {
        RawBitLengthZeroMode::Legacy => {
            std::env::remove_var("LOWQ_REUSE_ZERO_CARRIES_FOR_FULL_PREFIX_SCRATCH");
            std::env::remove_var("LOWQ_FUSED_ZERO_PREFIX_BITLEN");
        }
        RawBitLengthZeroMode::Fused => {
            std::env::set_var("LOWQ_REUSE_ZERO_CARRIES_FOR_FULL_PREFIX_SCRATCH", "1");
            std::env::set_var("LOWQ_FUSED_ZERO_PREFIX_BITLEN", "1");
        }
    }
    if loan {
        std::env::set_var(COEFFICIENT_RAW_BITLEN_LOAN_FLAG, "1");
    } else {
        std::env::remove_var(COEFFICIENT_RAW_BITLEN_LOAN_FLAG);
    }
    std::env::remove_var(INPLACE_ROTATED_BITLEN_BOUNDARY_FLAG);
    std::env::remove_var(FUSED_PREFIX_SCRATCH_LOAN_FLAG);
    std::env::remove_var(PROMISED_LQ_SWAP_BORROW_FLAG);
    std::env::remove_var(PROMISED_SWAP_SUPPORT_LIFETIME_FUSION_FLAG);
    std::env::remove_var(Q845_SWAP_ONLY_T_PRIME_LENGTH_FLAG);
    std::env::remove_var(SPLIT_COEFFICIENT_ROTATION_LIFETIME_FLAG);
    std::env::remove_var(COEFFICIENT_LESS_THAN_LANE_REUSE_FLAG);
    std::env::remove_var(
        super::shrunken_pz_state_machine::CALLER_SCRATCH_KG_REVERSE_DECREMENT_FLAG,
    );
    std::env::remove_var(CLEAN_CHAIN_COEFFICIENT_ADD_LENDER_FLAG);
}

struct RawBitLengthProofHarness {
    builder: B,
    data_ids: Vec<u32>,
    external_mask: u64,
    control_mask: u64,
    output_mask: u64,
    preserved_mask: u64,
    chain_mask: u64,
    add_only_mask: u64,
    borrowed_mask: u64,
    after_first_comparator: usize,
    raw_start: usize,
    raw_end: usize,
    before_second_comparator: usize,
    raw_trace: RawBitLengthAllocationTrace,
    zero_trace: DynamicBitLengthZeroAllocationTrace,
}

fn qreg_mask<'a>(lanes: impl IntoIterator<Item = &'a QReg>) -> u64 {
    lanes.into_iter().fold(0u64, |mask, lane| {
        mask | 1u64.checked_shl(lane.id()).unwrap_or(0)
    })
}

fn build_raw_bit_length_proof_harness(
    length_width: usize,
    work_width: usize,
    inverse: bool,
    direct_baseline: bool,
    with_comparators: bool,
) -> RawBitLengthProofHarness {
    assert!(length_width > 0);
    assert!(work_width <= (1usize << length_width) - 1);
    let mut circ = Circuit::new();
    let control = circ.alloc_qreg("rs.raw-loan-proof.control");
    let l_s = circ.alloc_qreg_bits("rs.raw-loan-proof.l-s", length_width);
    let l_r_prime = circ.alloc_qreg_bits("rs.raw-loan-proof.l-r-prime", length_width);
    let work2 = circ.alloc_qreg_bits("rs.raw-loan-proof.work2", work_width);
    let output = circ.alloc_qreg_bits("rs.raw-loan-proof.output", length_width);
    let comparator_before = circ.alloc_qreg("rs.raw-loan-proof.comparator-before");
    let comparator_after = circ.alloc_qreg("rs.raw-loan-proof.comparator-after");
    let chain = circ.alloc_qreg_bits("rs.raw-loan-proof.chain", 2);
    let add_only = circ.alloc_qreg("rs.raw-loan-proof.add-only");

    let data_ids: Vec<u32> = std::iter::once(&control)
        .chain(&l_s)
        .chain(&l_r_prime)
        .chain(&work2)
        .chain(&output)
        .map(QReg::id)
        .collect();
    let external_mask = qreg_mask(
        std::iter::once(&control)
            .chain(&l_s)
            .chain(&l_r_prime)
            .chain(&work2)
            .chain(&output)
            .chain(std::iter::once(&comparator_before))
            .chain(std::iter::once(&comparator_after))
            .chain(&chain)
            .chain(std::iter::once(&add_only)),
    );
    let output_mask = qreg_mask(&output);
    let preserved_mask = qreg_mask(
        std::iter::once(&control)
            .chain(&l_s)
            .chain(&l_r_prime)
            .chain(&work2),
    );
    let chain_mask = qreg_mask(&chain);
    let add_only_mask = 1u64 << add_only.id();
    let borrowed_mask = chain_mask | add_only_mask;
    let compare_scratch = vec![&chain[0], &chain[1], &add_only];
    let (first_target, second_target) = if inverse {
        (&comparator_after, &comparator_before)
    } else {
        (&comparator_before, &comparator_after)
    };

    begin_raw_bit_length_allocation_trace();
    begin_dynamic_bit_length_zero_allocation_trace();
    if with_comparators {
        toggle_coefficient_length_above_boundary(
            &mut circ,
            &output,
            &l_r_prime,
            &l_s,
            first_target,
            &compare_scratch,
        );
    }
    let after_first_comparator = circ.total_ops() as usize;
    circ.cx(&control, &add_only);
    let raw_start = circ.total_ops() as usize;
    if direct_baseline {
        controlled_xor_raw_t_prime_bit_length_allocated(
            &mut circ, &add_only, &l_s, &l_r_prime, &work2, &output,
        );
    } else {
        controlled_xor_raw_t_prime_bit_length(
            &mut circ, &add_only, &l_s, &l_r_prime, &work2, &output, &chain,
        );
    }
    let raw_end = circ.total_ops() as usize;
    circ.cx(&control, &add_only);
    let before_second_comparator = circ.total_ops() as usize;
    if with_comparators {
        toggle_coefficient_length_above_boundary(
            &mut circ,
            &output,
            &l_r_prime,
            &l_s,
            second_target,
            &compare_scratch,
        );
    }
    let zero_trace = finish_dynamic_bit_length_zero_allocation_trace();
    let raw_trace = finish_raw_bit_length_allocation_trace();
    drop(compare_scratch);

    RawBitLengthProofHarness {
        builder: circ.into_builder(),
        data_ids,
        external_mask,
        control_mask: 1u64 << control.id(),
        output_mask,
        preserved_mask,
        chain_mask,
        add_only_mask,
        borrowed_mask,
        after_first_comparator,
        raw_start,
        raw_end,
        before_second_comparator,
        raw_trace,
        zero_trace,
    }
}

fn raw_bit_length_input(harness: &RawBitLengthProofHarness, value: u64) -> u64 {
    harness
        .data_ids
        .iter()
        .enumerate()
        .fold(0u64, |state, (bit, id)| {
            state | (((value >> bit) & 1) << id)
        })
}

fn assert_raw_bit_length_harness_clean(
    harness: &RawBitLengthProofHarness,
    state: u64,
    context: &str,
) {
    assert_eq!(
        state & harness.borrowed_mask,
        0,
        "{context}: borrowed coefficient lanes were not restored"
    );
    assert_eq!(
        state & !harness.external_mask,
        0,
        "{context}: raw bit-length helper left an internal lane dirty"
    );
}

fn verify_raw_bit_length_simulator_equivalence(
    baseline: &RawBitLengthProofHarness,
    loan: &RawBitLengthProofHarness,
) -> (usize, usize) {
    use crate::circuit::QubitId;
    use crate::sim::Simulator;
    use sha3::{
        digest::{ExtendableOutput, Update},
        Shake128,
    };

    assert_eq!(baseline.data_ids, loan.data_ids);
    assert_eq!(baseline.external_mask, loan.external_mask);
    let states = 1usize << baseline.data_ids.len();
    let mut cases_checked = 0usize;
    let mut phase_clean_checks = 0usize;
    for batch_start in (0..states).step_by(64) {
        let shots = (states - batch_start).min(64);
        let live = if shots == 64 {
            u64::MAX
        } else {
            (1u64 << shots) - 1
        };
        let mut baseline_seed = Shake128::default();
        baseline_seed.update(b"coefficient-raw-bitlen-loan-baseline");
        baseline_seed.update(&(batch_start as u64).to_le_bytes());
        let mut baseline_xof = baseline_seed.finalize_xof();
        let mut baseline_simulator = Simulator::new(
            baseline.builder.next_qubit as usize,
            baseline.builder.next_bit as usize,
            &mut baseline_xof,
        );
        let mut loan_seed = Shake128::default();
        loan_seed.update(b"coefficient-raw-bitlen-loan-baseline");
        loan_seed.update(&(batch_start as u64).to_le_bytes());
        let mut loan_xof = loan_seed.finalize_xof();
        let mut loan_simulator = Simulator::new(
            loan.builder.next_qubit as usize,
            loan.builder.next_bit as usize,
            &mut loan_xof,
        );

        for shot in 0..shots {
            let value = (batch_start + shot) as u64;
            for (bit, (&baseline_id, &loan_id)) in
                baseline.data_ids.iter().zip(&loan.data_ids).enumerate()
            {
                if (value >> bit) & 1 != 0 {
                    *baseline_simulator.qubit_mut(QubitId(u64::from(baseline_id))) |= 1u64 << shot;
                    *loan_simulator.qubit_mut(QubitId(u64::from(loan_id))) |= 1u64 << shot;
                }
            }
        }
        baseline_simulator.apply_iter(baseline.builder.ops.iter());
        loan_simulator.apply_iter(loan.builder.ops.iter());
        assert_eq!(baseline_simulator.phase & live, 0);
        assert_eq!(loan_simulator.phase & live, 0);
        phase_clean_checks += 2 * shots;

        for id in 0..baseline.builder.next_qubit {
            let value = baseline_simulator.qubit(QubitId(u64::from(id))) & live;
            if baseline.external_mask & (1u64 << id) != 0 {
                assert_eq!(value, loan_simulator.qubit(QubitId(u64::from(id))) & live);
            } else {
                assert_eq!(value, 0, "baseline simulator left q{id} dirty");
            }
        }
        for id in 0..loan.builder.next_qubit {
            if loan.external_mask & (1u64 << id) == 0 {
                assert_eq!(
                    loan_simulator.qubit(QubitId(u64::from(id))) & live,
                    0,
                    "loan simulator left q{id} dirty"
                );
            }
        }
        cases_checked += shots;
    }
    (cases_checked, phase_clean_checks)
}

fn verify_raw_bit_length_loan_boundaries(harness: &RawBitLengthProofHarness, input: u64) -> usize {
    let after_first = apply_scalar(
        &harness.builder.ops[..harness.after_first_comparator],
        input,
    );
    assert_eq!(after_first & harness.borrowed_mask, 0);

    let at_raw_start = apply_scalar(&harness.builder.ops[..harness.raw_start], input);
    assert_eq!(at_raw_start & harness.chain_mask, 0);
    assert_eq!(
        at_raw_start & harness.add_only_mask,
        if input & harness.control_mask == 0 {
            0
        } else {
            harness.add_only_mask
        }
    );

    let at_raw_end = apply_scalar(&harness.builder.ops[..harness.raw_end], input);
    assert_eq!(at_raw_end & harness.chain_mask, 0);
    assert_eq!(
        at_raw_end & harness.add_only_mask,
        at_raw_start & harness.add_only_mask
    );
    assert_eq!(at_raw_end & !harness.external_mask, 0);

    let before_second = apply_scalar(
        &harness.builder.ops[..harness.before_second_comparator],
        input,
    );
    assert_eq!(before_second & harness.borrowed_mask, 0);
    4
}

fn raw_bit_length_local_resources(
    harness: &RawBitLengthProofHarness,
) -> CoefficientRawBitLengthLoanLocalResources {
    let counts = measurement_classical_gate_counts(&harness.builder.ops);
    let active_qubits = harness.builder.active_qubits as usize;
    let peak_qubits = harness.builder.peak_qubits as usize;
    CoefficientRawBitLengthLoanLocalResources {
        active_qubits,
        peak_qubits,
        temporary_qubits: peak_qubits - active_qubits,
        emitted_ops: harness.builder.ops.len(),
        emitted_toffoli: counts.ccx,
    }
}

fn assert_raw_bit_length_default_stream(
    default: &RawBitLengthProofHarness,
    direct: &RawBitLengthProofHarness,
) {
    assert_eq!(default.builder.ops, direct.builder.ops);
    assert_eq!(default.builder.next_qubit, direct.builder.next_qubit);
    assert_eq!(default.builder.next_bit, direct.builder.next_bit);
    assert_eq!(default.builder.active_qubits, direct.builder.active_qubits);
    assert_eq!(default.builder.peak_qubits, direct.builder.peak_qubits);
    assert_eq!(default.builder.free_qubits, direct.builder.free_qubits);
    assert_eq!(
        default.builder.allocation_serial,
        direct.builder.allocation_serial
    );
}

fn build_no_overflow_cuccaro(width: usize, subtract: bool) -> B {
    let mut circ = Circuit::new();
    let a = circ.alloc_qreg_bits("rs.raw-loan-proof.no-overflow-a", width);
    let b = circ.alloc_qreg_bits("rs.raw-loan-proof.no-overflow-b", width);
    let carry = circ.alloc_qreg("rs.raw-loan-proof.no-overflow-carry");
    if subtract {
        cuccaro_sub_mod_2n_no_overflow(&mut circ, &a, &b, &carry);
    } else {
        cuccaro_add_mod_2n_no_overflow(&mut circ, &a, &b, &carry);
    }
    circ.into_builder()
}

fn exhaustive_no_overflow_cuccaro_check() -> (usize, usize, usize) {
    let mut basis_states_checked = 0usize;
    let mut inverse_pair_checks = 0usize;
    for width in 1..=4 {
        let add = build_no_overflow_cuccaro(width, false);
        let sub = build_no_overflow_cuccaro(width, true);
        let mask = (1u64 << width) - 1;
        for a in 0..=mask {
            for b in 0..=mask {
                let input = a | (b << width);
                let added = apply_scalar(&add.ops, input);
                let subtracted = apply_scalar(&sub.ops, input);
                assert_eq!(added, a | (((a + b) & mask) << width));
                assert_eq!(subtracted, a | ((b.wrapping_sub(a) & mask) << width));
                assert_eq!(apply_scalar(&sub.ops, added), input);
                assert_eq!(apply_scalar(&add.ops, subtracted), input);
                basis_states_checked += 1;
                inverse_pair_checks += 2;
            }
        }
    }
    (4, basis_states_checked, inverse_pair_checks)
}

fn coefficient_raw_bitlen_loan_precondition_rejections() -> usize {
    use std::panic::{catch_unwind, AssertUnwindSafe};

    let previous_hook = std::panic::take_hook();
    std::panic::set_hook(Box::new(|_| {}));
    let short_chain = catch_unwind(AssertUnwindSafe(|| {
        let mut circ = Circuit::new();
        let control = circ.alloc_qreg("rs.raw-loan-reject.control");
        let l_s = circ.alloc_qreg_bits("rs.raw-loan-reject.l-s", 2);
        let l_r = circ.alloc_qreg_bits("rs.raw-loan-reject.l-r", 2);
        let work = circ.alloc_qreg_bits("rs.raw-loan-reject.work", 2);
        let output = circ.alloc_qreg_bits("rs.raw-loan-reject.output", 2);
        let chain = circ.alloc_qreg_bits("rs.raw-loan-reject.chain", 1);
        controlled_xor_raw_t_prime_bit_length_loaned(
            &mut circ, &control, &l_s, &l_r, &work, &output, &chain,
        );
    }));
    let aliased_chain = catch_unwind(AssertUnwindSafe(|| {
        let mut circ = Circuit::new();
        let control = circ.alloc_qreg("rs.raw-loan-reject.control");
        let l_s = circ.alloc_qreg_bits("rs.raw-loan-reject.l-s", 2);
        let l_r = circ.alloc_qreg_bits("rs.raw-loan-reject.l-r", 2);
        let work = circ.alloc_qreg_bits("rs.raw-loan-reject.work", 2);
        let output = circ.alloc_qreg_bits("rs.raw-loan-reject.output", 2);
        controlled_xor_raw_t_prime_bit_length_loaned(
            &mut circ, &control, &l_s, &l_r, &work, &output, &work,
        );
    }));
    std::panic::set_hook(previous_hook);
    assert!(short_chain.is_err());
    assert!(aliased_chain.is_err());
    2
}

/// Exhaustively verify the coefficient raw-bitlength loan at reduced widths.
/// The harness places the borrowed comparator on both sides of the raw call,
/// checks the live `add_only` boundary, and separately proves the no-overflow
/// Cuccaro pair used by the outer rotation.
#[doc(hidden)]
pub fn exhaustive_coefficient_raw_bitlength_loan_check() -> CoefficientRawBitLengthLoanProofReport {
    assert!(
        std::env::var_os(COEFFICIENT_RAW_BITLEN_LOAN_FLAG).is_none(),
        "the coefficient raw bit-length loan must default off"
    );
    let _environment = RawBitLengthProofEnvironment::capture();
    let configurations = [(1usize, 1usize), (2, 1), (2, 2), (2, 3)];
    let modes = [RawBitLengthZeroMode::Legacy, RawBitLengthZeroMode::Fused];

    let mut basis_states_checked = 0usize;
    let mut baseline_equivalence_checks = 0usize;
    let mut simulator_equivalence_checks = 0usize;
    let mut phase_clean_checks = 0usize;
    let mut inverse_pair_checks = 0usize;
    let mut control_off_checks = 0usize;
    let mut add_only_one_checks = 0usize;
    let mut comparator_boundary_checks = 0usize;
    let mut borrowed_lane_clean_checks = 0usize;
    let mut ancilla_clean_checks = 0usize;
    let mut legacy_baseline_trace = None;
    let mut legacy_loan_trace = None;
    let mut fused_loan_zero_trace = None;
    let mut baseline_raw_trace = None;
    let mut loan_raw_trace = None;

    for mode in modes {
        for &(length_width, work_width) in &configurations {
            configure_raw_bit_length_loan_proof(mode, false);
            let baseline_forward =
                build_raw_bit_length_proof_harness(length_width, work_width, false, false, true);
            let baseline_inverse =
                build_raw_bit_length_proof_harness(length_width, work_width, true, false, true);
            configure_raw_bit_length_loan_proof(mode, true);
            let loan_forward =
                build_raw_bit_length_proof_harness(length_width, work_width, false, false, true);
            let loan_inverse =
                build_raw_bit_length_proof_harness(length_width, work_width, true, false, true);

            for harness in [&loan_forward, &loan_inverse] {
                assert_eq!(harness.raw_trace.raw_rotation_carry_allocations, 1);
                assert_eq!(harness.raw_trace.raw_rotation_overflow_allocations, 0);
                assert_eq!(harness.raw_trace.rotated_carry_allocations, 0);
                assert_eq!(harness.raw_trace.rotated_overflow_allocations, 0);
                assert_eq!(harness.raw_trace.rotated_enabled_allocations, 0);
                if mode == RawBitLengthZeroMode::Fused {
                    assert_eq!(
                        harness.zero_trace,
                        DynamicBitLengthZeroAllocationTrace::default()
                    );
                }
            }
            for harness in [&baseline_forward, &baseline_inverse] {
                assert_eq!(harness.raw_trace.raw_rotation_carry_allocations, 1);
                assert_eq!(harness.raw_trace.raw_rotation_overflow_allocations, 1);
                assert_eq!(harness.raw_trace.rotated_carry_allocations, 1);
                assert_eq!(harness.raw_trace.rotated_overflow_allocations, 1);
                assert_eq!(harness.raw_trace.rotated_enabled_allocations, 1);
            }
            assert_eq!(
                measurement_classical_gate_counts(&loan_forward.builder.ops).ccx,
                measurement_classical_gate_counts(&baseline_forward.builder.ops).ccx
            );
            assert_eq!(
                measurement_classical_gate_counts(&loan_inverse.builder.ops).ccx,
                measurement_classical_gate_counts(&baseline_inverse.builder.ops).ccx
            );
            for (baseline, loan) in [
                (&baseline_forward, &loan_forward),
                (&baseline_inverse, &loan_inverse),
            ] {
                let (cases, phase_checks) =
                    verify_raw_bit_length_simulator_equivalence(baseline, loan);
                simulator_equivalence_checks += cases;
                phase_clean_checks += phase_checks;
            }

            if (length_width, work_width) == (2, 3) {
                match mode {
                    RawBitLengthZeroMode::Legacy => {
                        legacy_baseline_trace = Some(baseline_forward.zero_trace);
                        legacy_loan_trace = Some(loan_forward.zero_trace);
                    }
                    RawBitLengthZeroMode::Fused => {
                        fused_loan_zero_trace = Some(loan_forward.zero_trace);
                        baseline_raw_trace = Some(baseline_forward.raw_trace);
                        loan_raw_trace = Some(loan_forward.raw_trace);
                    }
                }
            }

            assert_eq!(baseline_forward.data_ids, loan_forward.data_ids);
            assert_eq!(baseline_inverse.data_ids, loan_inverse.data_ids);
            let data_states = 1u64 << baseline_forward.data_ids.len();
            for value in 0..data_states {
                let input = raw_bit_length_input(&baseline_forward, value);
                let baseline_output = apply_scalar(&baseline_forward.builder.ops, input);
                let loan_output = apply_scalar(&loan_forward.builder.ops, input);
                assert_eq!(loan_output, baseline_output);
                assert_raw_bit_length_harness_clean(
                    &baseline_forward,
                    baseline_output,
                    "baseline forward",
                );
                assert_raw_bit_length_harness_clean(&loan_forward, loan_output, "loan forward");

                let baseline_inverse_output = apply_scalar(&baseline_inverse.builder.ops, input);
                let loan_inverse_output = apply_scalar(&loan_inverse.builder.ops, input);
                assert_eq!(loan_inverse_output, baseline_inverse_output);
                assert_raw_bit_length_harness_clean(
                    &baseline_inverse,
                    baseline_inverse_output,
                    "baseline inverse",
                );
                assert_raw_bit_length_harness_clean(
                    &loan_inverse,
                    loan_inverse_output,
                    "loan inverse",
                );

                assert_eq!(
                    apply_scalar(&baseline_inverse.builder.ops, baseline_output),
                    input
                );
                assert_eq!(apply_scalar(&loan_inverse.builder.ops, loan_output), input);
                assert_eq!(
                    loan_output & loan_forward.preserved_mask,
                    input & loan_forward.preserved_mask
                );
                if input & loan_forward.control_mask == 0 {
                    assert_eq!(
                        loan_output & loan_forward.output_mask,
                        input & loan_forward.output_mask
                    );
                    control_off_checks += 2;
                } else {
                    add_only_one_checks += 2;
                }

                comparator_boundary_checks +=
                    verify_raw_bit_length_loan_boundaries(&loan_forward, input);
                comparator_boundary_checks +=
                    verify_raw_bit_length_loan_boundaries(&loan_inverse, input);
                borrowed_lane_clean_checks += 4;
                ancilla_clean_checks += 4;
                basis_states_checked += 2;
                baseline_equivalence_checks += 2;
                inverse_pair_checks += 2;
            }
        }
    }

    configure_raw_bit_length_loan_proof(RawBitLengthZeroMode::Fused, false);
    let default_forward = build_raw_bit_length_proof_harness(2, 3, false, false, true);
    let direct_forward = build_raw_bit_length_proof_harness(2, 3, false, true, true);
    assert_raw_bit_length_default_stream(&default_forward, &direct_forward);
    let default_inverse = build_raw_bit_length_proof_harness(2, 3, true, false, true);
    let direct_inverse = build_raw_bit_length_proof_harness(2, 3, true, true, true);
    assert_raw_bit_length_default_stream(&default_inverse, &direct_inverse);

    configure_raw_bit_length_loan_proof(RawBitLengthZeroMode::Fused, false);
    let baseline_local = raw_bit_length_local_resources(&build_raw_bit_length_proof_harness(
        REFERENCE_LENGTH_WIDTH,
        259,
        false,
        false,
        false,
    ));
    configure_raw_bit_length_loan_proof(RawBitLengthZeroMode::Fused, true);
    let loan_local = raw_bit_length_local_resources(&build_raw_bit_length_proof_harness(
        REFERENCE_LENGTH_WIDTH,
        259,
        false,
        false,
        false,
    ));
    assert_eq!(loan_local.active_qubits, baseline_local.active_qubits);
    assert_eq!(loan_local.emitted_toffoli, baseline_local.emitted_toffoli);
    assert!(loan_local.peak_qubits < baseline_local.peak_qubits);

    let legacy_baseline_trace = legacy_baseline_trace.expect("legacy baseline trace");
    let legacy_loan_trace = legacy_loan_trace.expect("legacy loan trace");
    assert_eq!(legacy_baseline_trace.flag_allocations, 2);
    assert_eq!(legacy_loan_trace.flag_allocations, 2);
    assert_eq!(legacy_baseline_trace.carry_allocations, 4);
    assert_eq!(legacy_loan_trace.carry_allocations, 2);
    assert!(legacy_baseline_trace.prefix_allocations > 0);
    assert_eq!(
        legacy_loan_trace.prefix_allocations,
        legacy_baseline_trace.prefix_allocations
    );
    let fused_loan_zero_trace = fused_loan_zero_trace.expect("fused loan trace");
    assert_eq!(
        fused_loan_zero_trace,
        DynamicBitLengthZeroAllocationTrace::default()
    );
    let baseline_raw_trace = baseline_raw_trace.expect("baseline raw trace");
    let loan_raw_trace = loan_raw_trace.expect("loan raw trace");
    let (no_overflow_widths, no_overflow_states, no_overflow_inverse_pairs) =
        exhaustive_no_overflow_cuccaro_check();
    let precondition_rejections = coefficient_raw_bitlen_loan_precondition_rejections();

    CoefficientRawBitLengthLoanProofReport {
        configurations_checked: configurations.len(),
        zero_modes_checked: modes.len(),
        directions_checked: 2,
        basis_states_checked,
        baseline_equivalence_checks,
        simulator_equivalence_checks,
        phase_clean_checks,
        inverse_pair_checks,
        control_off_checks,
        add_only_one_checks,
        comparator_boundary_checks,
        borrowed_lane_clean_checks,
        ancilla_clean_checks,
        no_overflow_widths_checked: no_overflow_widths,
        no_overflow_basis_states_checked: no_overflow_states,
        no_overflow_inverse_pair_checks: no_overflow_inverse_pairs,
        default_stream_identity_checks: 2,
        precondition_rejections,
        legacy_baseline_zero_carry_allocations: legacy_baseline_trace.carry_allocations,
        legacy_loan_zero_carry_allocations: legacy_loan_trace.carry_allocations,
        fused_loan_zero_flag_allocations: fused_loan_zero_trace.flag_allocations,
        fused_loan_zero_carry_allocations: fused_loan_zero_trace.carry_allocations,
        fused_loan_zero_prefix_allocations: fused_loan_zero_trace.prefix_allocations,
        baseline_raw_rotation_carry_allocations: baseline_raw_trace.raw_rotation_carry_allocations,
        baseline_raw_rotation_overflow_allocations: baseline_raw_trace
            .raw_rotation_overflow_allocations,
        baseline_rotated_carry_allocations: baseline_raw_trace.rotated_carry_allocations,
        baseline_rotated_overflow_allocations: baseline_raw_trace.rotated_overflow_allocations,
        baseline_rotated_enabled_allocations: baseline_raw_trace.rotated_enabled_allocations,
        loan_raw_rotation_carry_allocations: loan_raw_trace.raw_rotation_carry_allocations,
        loan_raw_rotation_overflow_allocations: loan_raw_trace.raw_rotation_overflow_allocations,
        loan_rotated_carry_allocations: loan_raw_trace.rotated_carry_allocations,
        loan_rotated_overflow_allocations: loan_raw_trace.rotated_overflow_allocations,
        loan_rotated_enabled_allocations: loan_raw_trace.rotated_enabled_allocations,
        baseline_local,
        loan_local,
        local_qubit_delta: loan_local.peak_qubits as i64 - baseline_local.peak_qubits as i64,
        local_toffoli_delta: loan_local.emitted_toffoli as i64
            - baseline_local.emitted_toffoli as i64,
    }
}

struct InplaceRotatedBoundaryHarness {
    builder: B,
    data_ids: Vec<u32>,
    external_mask: u64,
    control_mask: u64,
    output_mask: u64,
    preserved_mask: u64,
    scratch_mask: u64,
    raw_trace: RawBitLengthAllocationTrace,
}

fn build_inplace_rotated_boundary_harness(
    length_width: usize,
    source_width: usize,
    scratch_lanes: usize,
    route: SaturatingDifferenceBoundaryRoute,
) -> InplaceRotatedBoundaryHarness {
    assert!(length_width > 0);
    assert!(source_width <= (1usize << length_width) - 1);
    assert!(scratch_lanes == 0 || scratch_lanes >= 3);
    let mut circ = Circuit::new();
    let control = circ.alloc_qreg("rs.inplace-boundary-proof.control");
    let boundary = circ.alloc_qreg_bits("rs.inplace-boundary-proof.boundary", length_width);
    let source = circ.alloc_qreg_bits("rs.inplace-boundary-proof.source", source_width);
    let output = circ.alloc_qreg_bits("rs.inplace-boundary-proof.output", length_width);
    let scratch = circ.alloc_qreg_bits("rs.inplace-boundary-proof.scratch", scratch_lanes);
    let source_refs: Vec<&QReg> = source.iter().collect();
    let scratch_refs: Vec<&QReg> = scratch.iter().collect();
    let data_ids: Vec<u32> = std::iter::once(&control)
        .chain(&boundary)
        .chain(&source)
        .chain(&output)
        .map(QReg::id)
        .collect();
    let external_mask = qreg_mask(
        std::iter::once(&control)
            .chain(&boundary)
            .chain(&source)
            .chain(&output)
            .chain(&scratch),
    );
    let preserved_mask = qreg_mask(std::iter::once(&control).chain(&boundary).chain(&source));
    let output_mask = qreg_mask(&output);
    let scratch_mask = qreg_mask(&scratch);

    begin_raw_bit_length_allocation_trace();
    controlled_xor_saturating_bit_length_difference_with_route(
        &mut circ,
        &control,
        &boundary,
        &source_refs,
        &output,
        &scratch_refs,
        None,
        route,
    );
    let raw_trace = finish_raw_bit_length_allocation_trace();
    InplaceRotatedBoundaryHarness {
        builder: circ.into_builder(),
        data_ids,
        external_mask,
        control_mask: 1u64 << control.id(),
        output_mask,
        preserved_mask,
        scratch_mask,
        raw_trace,
    }
}

fn inplace_rotated_boundary_input(harness: &InplaceRotatedBoundaryHarness, value: u64) -> u64 {
    harness
        .data_ids
        .iter()
        .enumerate()
        .fold(0u64, |state, (bit, id)| {
            state | (((value >> bit) & 1) << id)
        })
}

fn assert_inplace_rotated_boundary_clean(
    harness: &InplaceRotatedBoundaryHarness,
    state: u64,
    context: &str,
) {
    assert_eq!(
        state & harness.scratch_mask,
        0,
        "{context}: borrowed scratch dirty"
    );
    assert_eq!(
        state & !harness.external_mask,
        0,
        "{context}: internal scratch dirty"
    );
}

fn assert_inplace_rotated_boundary_stream_identity(
    left: &InplaceRotatedBoundaryHarness,
    right: &InplaceRotatedBoundaryHarness,
) {
    assert_eq!(left.builder.ops, right.builder.ops);
    assert_eq!(left.builder.next_qubit, right.builder.next_qubit);
    assert_eq!(left.builder.next_bit, right.builder.next_bit);
    assert_eq!(left.builder.active_qubits, right.builder.active_qubits);
    assert_eq!(left.builder.peak_qubits, right.builder.peak_qubits);
    assert_eq!(left.builder.free_qubits, right.builder.free_qubits);
    assert_eq!(
        left.builder.allocation_serial,
        right.builder.allocation_serial
    );
}

fn verify_inplace_rotated_boundary_simulator_equivalence(
    baseline: &InplaceRotatedBoundaryHarness,
    candidate: &InplaceRotatedBoundaryHarness,
) -> (usize, usize) {
    use crate::circuit::QubitId;
    use crate::sim::Simulator;
    use sha3::{
        digest::{ExtendableOutput, Update},
        Shake128,
    };

    assert_eq!(baseline.data_ids, candidate.data_ids);
    assert_eq!(baseline.external_mask, candidate.external_mask);
    let states = 1usize << baseline.data_ids.len();
    let mut cases_checked = 0usize;
    let mut phase_clean_checks = 0usize;
    for batch_start in (0..states).step_by(64) {
        let shots = (states - batch_start).min(64);
        let live = if shots == 64 {
            u64::MAX
        } else {
            (1u64 << shots) - 1
        };
        let mut baseline_seed = Shake128::default();
        baseline_seed.update(b"inplace-rotated-boundary-proof");
        baseline_seed.update(&(batch_start as u64).to_le_bytes());
        let mut baseline_xof = baseline_seed.finalize_xof();
        let mut baseline_simulator = Simulator::new(
            baseline.builder.next_qubit as usize,
            baseline.builder.next_bit as usize,
            &mut baseline_xof,
        );
        let mut candidate_seed = Shake128::default();
        candidate_seed.update(b"inplace-rotated-boundary-proof");
        candidate_seed.update(&(batch_start as u64).to_le_bytes());
        let mut candidate_xof = candidate_seed.finalize_xof();
        let mut candidate_simulator = Simulator::new(
            candidate.builder.next_qubit as usize,
            candidate.builder.next_bit as usize,
            &mut candidate_xof,
        );

        for shot in 0..shots {
            let value = (batch_start + shot) as u64;
            for (bit, (&baseline_id, &candidate_id)) in baseline
                .data_ids
                .iter()
                .zip(&candidate.data_ids)
                .enumerate()
            {
                if (value >> bit) & 1 != 0 {
                    *baseline_simulator.qubit_mut(QubitId(u64::from(baseline_id))) |= 1u64 << shot;
                    *candidate_simulator.qubit_mut(QubitId(u64::from(candidate_id))) |=
                        1u64 << shot;
                }
            }
        }
        baseline_simulator.apply_iter(baseline.builder.ops.iter());
        candidate_simulator.apply_iter(candidate.builder.ops.iter());
        assert_eq!(baseline_simulator.phase & live, 0);
        assert_eq!(candidate_simulator.phase & live, 0);
        phase_clean_checks += 2 * shots;

        for id in 0..baseline.builder.next_qubit {
            let baseline_value = baseline_simulator.qubit(QubitId(u64::from(id))) & live;
            if baseline.external_mask & (1u64 << id) != 0 {
                assert_eq!(
                    baseline_value,
                    candidate_simulator.qubit(QubitId(u64::from(id))) & live
                );
            } else {
                assert_eq!(baseline_value, 0, "baseline simulator left q{id} dirty");
            }
        }
        for id in 0..candidate.builder.next_qubit {
            if candidate.external_mask & (1u64 << id) == 0 {
                assert_eq!(
                    candidate_simulator.qubit(QubitId(u64::from(id))) & live,
                    0,
                    "candidate simulator left q{id} dirty"
                );
            }
        }
        cases_checked += shots;
    }
    (cases_checked, phase_clean_checks)
}

fn inplace_rotated_boundary_local_resources(
    harness: &InplaceRotatedBoundaryHarness,
) -> InplaceRotatedBoundaryLocalResources {
    let counts = measurement_classical_gate_counts(&harness.builder.ops);
    let active_qubits = harness.builder.active_qubits as usize;
    let peak_qubits = harness.builder.peak_qubits as usize;
    InplaceRotatedBoundaryLocalResources {
        active_qubits,
        peak_qubits,
        temporary_qubits: peak_qubits - active_qubits,
        emitted_ops: harness.builder.ops.len(),
        emitted_toffoli: counts.ccx,
    }
}

fn configure_inplace_rotated_boundary_proof(mode: RawBitLengthZeroMode) {
    configure_raw_bit_length_loan_proof(mode, false);
    std::env::remove_var(INPLACE_ROTATED_BITLEN_BOUNDARY_FLAG);
    std::env::remove_var(SUB800_BORROWED_ROTATED_UNDERFLOW_FLAG);
    std::env::remove_var(SUB800_SPLIT_MIXED_ROTATED_LENGTH_FLAG);
    std::env::remove_var(SUB800_SPLIT_SAME_ROTATED_LENGTH_FLAG);
    std::env::remove_var(SUB800_SPLIT_TWO_HIGH_ROTATED_LENGTH_FLAG);
    std::env::remove_var(SUB800_SPLIT_THREE_HIGH_ROTATED_LENGTH_FLAG);
    std::env::remove_var(SUB800_SPLIT_FOUR_HIGH_ROTATED_LENGTH_FLAG);
    std::env::remove_var(SUB800_ULS_FUSED_TARGET_FLAG);
    std::env::remove_var(SUB800_ULS_DIRECT_SELECTOR_FLAG);
}

fn bit_length_usize(value: usize) -> usize {
    if value == 0 {
        0
    } else {
        usize::BITS as usize - value.leading_zeros() as usize
    }
}

fn build_signed_boundary_arithmetic_kernel(width: usize, inplace: bool, inverse: bool) -> B {
    let mut circ = Circuit::new();
    if inplace {
        let boundary = circ.alloc_qreg_bits("rs.inplace-arithmetic.boundary", width);
        let length = circ.alloc_qreg_bits("rs.inplace-arithmetic.length", width + 1);
        let carry = circ.alloc_qreg("rs.inplace-arithmetic.carry");
        let (low_length, sign_lane) = length.split_at(width);
        if inverse {
            cuccaro_add_mod_2n(&mut circ, &boundary, low_length, &carry, &sign_lane[0]);
        } else {
            cuccaro_sub_mod_2n(&mut circ, &boundary, low_length, &carry, &sign_lane[0]);
        }
        circ.into_builder()
    } else {
        let boundary = circ.alloc_qreg_bits("rs.materialized-arithmetic.boundary", width + 1);
        let length = circ.alloc_qreg_bits("rs.materialized-arithmetic.length", width + 1);
        let carry = circ.alloc_qreg("rs.materialized-arithmetic.carry");
        let overflow = circ.alloc_qreg("rs.materialized-arithmetic.overflow");
        if inverse {
            cuccaro_add_mod_2n(&mut circ, &boundary, &length, &carry, &overflow);
        } else {
            cuccaro_sub_mod_2n(&mut circ, &boundary, &length, &carry, &overflow);
        }
        circ.into_builder()
    }
}

fn exhaustive_inplace_signed_boundary_arithmetic_check() -> (usize, usize, usize) {
    let mut states_checked = 0usize;
    let mut inverse_pair_checks = 0usize;
    for width in 1..=4 {
        let materialized_sub = build_signed_boundary_arithmetic_kernel(width, false, false);
        let materialized_add = build_signed_boundary_arithmetic_kernel(width, false, true);
        let inplace_sub = build_signed_boundary_arithmetic_kernel(width, true, false);
        let inplace_add = build_signed_boundary_arithmetic_kernel(width, true, true);
        let value_mask = (1u64 << width) - 1;
        let signed_mask = (1u64 << (width + 1)) - 1;
        for boundary in 0..=value_mask {
            for length in 0..=value_mask {
                let materialized_input = boundary | (length << (width + 1));
                let inplace_input = boundary | (length << width);
                let materialized_output = apply_scalar(&materialized_sub.ops, materialized_input);
                let inplace_output = apply_scalar(&inplace_sub.ops, inplace_input);
                let materialized_length = (materialized_output >> (width + 1)) & signed_mask;
                let inplace_length = (inplace_output >> width) & signed_mask;
                assert_eq!(
                    inplace_length, materialized_length,
                    "signed-boundary mismatch width={width} boundary={boundary} length={length}"
                );
                let materialized_added = apply_scalar(&materialized_add.ops, materialized_input);
                let inplace_added = apply_scalar(&inplace_add.ops, inplace_input);
                let materialized_added_length = (materialized_added >> (width + 1)) & signed_mask;
                let inplace_added_length = (inplace_added >> width) & signed_mask;
                assert_eq!(inplace_added_length, materialized_added_length);
                assert_eq!(
                    apply_scalar(&materialized_sub.ops, materialized_added),
                    materialized_input,
                    "materialized inverse mismatch width={width} boundary={boundary} length={length}"
                );
                assert_eq!(
                    apply_scalar(&inplace_sub.ops, inplace_added),
                    inplace_input,
                    "in-place inverse mismatch width={width} boundary={boundary} length={length}"
                );
                states_checked += 1;
                inverse_pair_checks += 2;
            }
        }
    }
    (4, states_checked, inverse_pair_checks)
}

/// Prove the allocation-free signed-boundary arithmetic against the previous
/// materialized-boundary helper, including every underflow case at reduced
/// widths and both production bit-length-zero implementations.
#[doc(hidden)]
pub fn exhaustive_inplace_rotated_bitlen_boundary_check() -> InplaceRotatedBoundaryProofReport {
    assert!(
        std::env::var_os(INPLACE_ROTATED_BITLEN_BOUNDARY_FLAG).is_none(),
        "the in-place rotated boundary route must default off"
    );
    let _environment = RawBitLengthProofEnvironment::capture();
    let configurations = [(1usize, 1usize), (2, 1), (2, 3), (3, 4)];
    let modes = [RawBitLengthZeroMode::Legacy, RawBitLengthZeroMode::Fused];
    let scratch_modes = [0usize, 3usize];

    let mut basis_states_checked = 0usize;
    let mut scalar_equivalence_checks = 0usize;
    let mut simulator_equivalence_checks = 0usize;
    let mut phase_clean_checks = 0usize;
    let mut inverse_pair_checks = 0usize;
    let mut control_off_checks = 0usize;
    let mut saturation_underflow_checks = 0usize;
    let mut nonnegative_difference_checks = 0usize;
    let mut preserved_input_checks = 0usize;
    let mut borrowed_scratch_clean_checks = 0usize;
    let mut ancilla_clean_checks = 0usize;

    for mode in modes {
        configure_inplace_rotated_boundary_proof(mode);
        for scratch_lanes in scratch_modes {
            for &(length_width, source_width) in &configurations {
                let baseline = build_inplace_rotated_boundary_harness(
                    length_width,
                    source_width,
                    scratch_lanes,
                    SaturatingDifferenceBoundaryRoute::Materialized,
                );
                let candidate = build_inplace_rotated_boundary_harness(
                    length_width,
                    source_width,
                    scratch_lanes,
                    SaturatingDifferenceBoundaryRoute::Inplace,
                );
                assert_eq!(baseline.data_ids, candidate.data_ids);
                let (simulator_cases, simulator_phases) =
                    verify_inplace_rotated_boundary_simulator_equivalence(&baseline, &candidate);
                simulator_equivalence_checks += simulator_cases;
                phase_clean_checks += simulator_phases;

                let data_states = 1u64 << baseline.data_ids.len();
                let length_mask = (1u64 << length_width) - 1;
                let source_mask = (1u64 << source_width) - 1;
                for value in 0..data_states {
                    let input = inplace_rotated_boundary_input(&baseline, value);
                    let baseline_output = apply_scalar(&baseline.builder.ops, input);
                    let candidate_output = apply_scalar(&candidate.builder.ops, input);
                    assert_eq!(candidate_output, baseline_output);
                    assert_inplace_rotated_boundary_clean(
                        &baseline,
                        baseline_output,
                        "materialized boundary",
                    );
                    assert_inplace_rotated_boundary_clean(
                        &candidate,
                        candidate_output,
                        "in-place boundary",
                    );
                    assert_eq!(
                        candidate_output & candidate.preserved_mask,
                        input & candidate.preserved_mask
                    );
                    assert_eq!(apply_scalar(&baseline.builder.ops, baseline_output), input);
                    assert_eq!(
                        apply_scalar(&candidate.builder.ops, candidate_output),
                        input
                    );
                    if input & candidate.control_mask == 0 {
                        assert_eq!(
                            candidate_output & candidate.output_mask,
                            input & candidate.output_mask
                        );
                        control_off_checks += 1;
                    }
                    let boundary = (value >> 1) & length_mask;
                    let source = (value >> (1 + length_width)) & source_mask;
                    if bit_length_usize(source as usize) < boundary as usize {
                        saturation_underflow_checks += 1;
                    } else {
                        nonnegative_difference_checks += 1;
                    }
                    if scratch_lanes != 0 {
                        borrowed_scratch_clean_checks += 2;
                    }
                    ancilla_clean_checks += 2;
                    preserved_input_checks += 1;
                    inverse_pair_checks += 2;
                    scalar_equivalence_checks += 1;
                    basis_states_checked += 1;
                }
            }
        }
    }

    configure_inplace_rotated_boundary_proof(RawBitLengthZeroMode::Fused);
    let mut default_stream_identity_checks = 0usize;
    let mut configured_stream_selection_checks = 0usize;
    for scratch_lanes in scratch_modes {
        let configured_default = build_inplace_rotated_boundary_harness(
            2,
            3,
            scratch_lanes,
            SaturatingDifferenceBoundaryRoute::Configured,
        );
        let direct_materialized = build_inplace_rotated_boundary_harness(
            2,
            3,
            scratch_lanes,
            SaturatingDifferenceBoundaryRoute::Materialized,
        );
        assert_inplace_rotated_boundary_stream_identity(&configured_default, &direct_materialized);
        default_stream_identity_checks += 1;

        std::env::set_var(INPLACE_ROTATED_BITLEN_BOUNDARY_FLAG, "1");
        let configured_candidate = build_inplace_rotated_boundary_harness(
            2,
            3,
            scratch_lanes,
            SaturatingDifferenceBoundaryRoute::Configured,
        );
        let direct_candidate = build_inplace_rotated_boundary_harness(
            2,
            3,
            scratch_lanes,
            SaturatingDifferenceBoundaryRoute::Inplace,
        );
        assert_inplace_rotated_boundary_stream_identity(&configured_candidate, &direct_candidate);
        configured_stream_selection_checks += 1;
        std::env::remove_var(INPLACE_ROTATED_BITLEN_BOUNDARY_FLAG);
    }

    let baseline_local_harness = build_inplace_rotated_boundary_harness(
        REFERENCE_LENGTH_WIDTH,
        259,
        7,
        SaturatingDifferenceBoundaryRoute::Materialized,
    );
    let candidate_local_harness = build_inplace_rotated_boundary_harness(
        REFERENCE_LENGTH_WIDTH,
        259,
        7,
        SaturatingDifferenceBoundaryRoute::Inplace,
    );
    let baseline_local = inplace_rotated_boundary_local_resources(&baseline_local_harness);
    let candidate_local = inplace_rotated_boundary_local_resources(&candidate_local_harness);
    assert_eq!(candidate_local.active_qubits, baseline_local.active_qubits);
    assert_eq!(candidate_local.peak_qubits + 1, baseline_local.peak_qubits);
    assert!(candidate_local.emitted_toffoli < baseline_local.emitted_toffoli);
    assert_eq!(
        baseline_local_harness
            .raw_trace
            .rotated_boundary_qubits_allocated,
        REFERENCE_LENGTH_WIDTH + 1
    );
    assert_eq!(
        candidate_local_harness
            .raw_trace
            .rotated_boundary_qubits_allocated,
        0
    );
    assert_eq!(
        candidate_local_harness
            .raw_trace
            .rotated_inplace_boundary_uses,
        1
    );

    let (arithmetic_widths, arithmetic_states, arithmetic_inverse_pairs) =
        exhaustive_inplace_signed_boundary_arithmetic_check();
    InplaceRotatedBoundaryProofReport {
        configurations_checked: configurations.len(),
        zero_modes_checked: modes.len(),
        scratch_modes_checked: scratch_modes.len(),
        basis_states_checked,
        scalar_equivalence_checks,
        simulator_equivalence_checks,
        phase_clean_checks,
        inverse_pair_checks,
        control_off_checks,
        saturation_underflow_checks,
        nonnegative_difference_checks,
        preserved_input_checks,
        borrowed_scratch_clean_checks,
        ancilla_clean_checks,
        arithmetic_widths_checked: arithmetic_widths,
        arithmetic_basis_states_checked: arithmetic_states,
        arithmetic_inverse_pair_checks: arithmetic_inverse_pairs,
        default_stream_identity_checks,
        configured_stream_selection_checks,
        baseline_boundary_qubits_allocated: baseline_local_harness
            .raw_trace
            .rotated_boundary_qubits_allocated,
        candidate_boundary_qubits_allocated: candidate_local_harness
            .raw_trace
            .rotated_boundary_qubits_allocated,
        candidate_inplace_boundary_uses: candidate_local_harness
            .raw_trace
            .rotated_inplace_boundary_uses,
        baseline_local,
        candidate_local,
        local_qubit_delta: candidate_local.peak_qubits as i64 - baseline_local.peak_qubits as i64,
        local_toffoli_delta: candidate_local.emitted_toffoli as i64
            - baseline_local.emitted_toffoli as i64,
    }
}

fn configure_paired_bitlength_source_complement_proof(
    prefix_scratch_loan: bool,
    configured_inplace: bool,
) {
    configure_raw_bit_length_loan_proof(RawBitLengthZeroMode::Fused, false);
    if prefix_scratch_loan {
        std::env::set_var(FUSED_PREFIX_SCRATCH_LOAN_FLAG, "1");
    } else {
        std::env::remove_var(FUSED_PREFIX_SCRATCH_LOAN_FLAG);
    }
    if configured_inplace {
        std::env::set_var(INPLACE_ROTATED_BITLEN_BOUNDARY_FLAG, "1");
    } else {
        std::env::remove_var(INPLACE_ROTATED_BITLEN_BOUNDARY_FLAG);
    }
    std::env::remove_var(PAIRED_BITLEN_SOURCE_COMPLEMENT_FLAG);
    std::env::remove_var(SUB800_BORROWED_ROTATED_UNDERFLOW_FLAG);
    std::env::remove_var(SUB800_SPLIT_MIXED_ROTATED_LENGTH_FLAG);
    std::env::remove_var(SUB800_SPLIT_SAME_ROTATED_LENGTH_FLAG);
    std::env::remove_var(SUB800_SPLIT_TWO_HIGH_ROTATED_LENGTH_FLAG);
    std::env::remove_var(SUB800_SPLIT_THREE_HIGH_ROTATED_LENGTH_FLAG);
    std::env::remove_var(SUB800_SPLIT_FOUR_HIGH_ROTATED_LENGTH_FLAG);
    std::env::remove_var(SUB800_ULS_FUSED_TARGET_FLAG);
    std::env::remove_var(SUB800_ULS_DIRECT_SELECTOR_FLAG);
}

fn set_paired_bitlength_source_complement_proof_mode(enabled: Option<bool>) {
    match enabled {
        Some(true) => std::env::set_var(PAIRED_BITLEN_SOURCE_COMPLEMENT_FLAG, "1"),
        Some(false) => std::env::set_var(PAIRED_BITLEN_SOURCE_COMPLEMENT_FLAG, "0"),
        None => std::env::remove_var(PAIRED_BITLEN_SOURCE_COMPLEMENT_FLAG),
    }
}

fn build_paired_bitlength_source_complement_harness(
    output_width: usize,
    source_width: usize,
    scratch_lanes: usize,
    route: SaturatingDifferenceBoundaryRoute,
    mixed_boundary: bool,
    enabled: Option<bool>,
) -> InplaceRotatedBoundaryHarness {
    assert!(output_width > usize::from(mixed_boundary));
    assert!(source_width > 1);
    assert!(source_width <= (1usize << output_width) - 1);
    assert!(scratch_lanes == 0 || scratch_lanes >= 3);
    set_paired_bitlength_source_complement_proof_mode(enabled);

    let boundary_width = output_width - usize::from(mixed_boundary);
    let mut circ = Circuit::new();
    let control = circ.alloc_qreg("rs.paired-bitlen-proof.control");
    let boundary =
        circ.alloc_qreg_bits("rs.paired-bitlen-proof.boundary", boundary_width);
    let source = circ.alloc_qreg_bits("rs.paired-bitlen-proof.source", source_width);
    let output = circ.alloc_qreg_bits("rs.paired-bitlen-proof.output", output_width);
    let scratch = circ.alloc_qreg_bits("rs.paired-bitlen-proof.scratch", scratch_lanes);
    let source_refs: Vec<&QReg> = source.iter().collect();
    let scratch_refs: Vec<&QReg> = scratch.iter().collect();
    let data_ids: Vec<u32> = std::iter::once(&control)
        .chain(&boundary)
        .chain(&source)
        .chain(&output)
        .map(QReg::id)
        .collect();
    let external_mask = qreg_mask(
        std::iter::once(&control)
            .chain(&boundary)
            .chain(&source)
            .chain(&output)
            .chain(&scratch),
    );
    let preserved_mask = qreg_mask(std::iter::once(&control).chain(&boundary).chain(&source));
    let output_mask = qreg_mask(&output);
    let scratch_mask = qreg_mask(&scratch);

    begin_raw_bit_length_allocation_trace();
    controlled_xor_saturating_bit_length_difference_with_route(
        &mut circ,
        &control,
        &boundary,
        &source_refs,
        &output,
        &scratch_refs,
        None,
        route,
    );
    let raw_trace = finish_raw_bit_length_allocation_trace();
    InplaceRotatedBoundaryHarness {
        builder: circ.into_builder(),
        data_ids,
        external_mask,
        control_mask: 1u64 << control.id(),
        output_mask,
        preserved_mask,
        scratch_mask,
        raw_trace,
    }
}

fn paired_bitlength_source_complement_local_resources(
    harness: &InplaceRotatedBoundaryHarness,
) -> PairedBitLengthSourceComplementLocalResources {
    use crate::circuit::OperationType;

    let active_qubits = harness.builder.active_qubits as usize;
    let peak_qubits = harness.builder.peak_qubits as usize;
    PairedBitLengthSourceComplementLocalResources {
        active_qubits,
        peak_qubits,
        temporary_qubits: peak_qubits - active_qubits,
        emitted_ops: harness.builder.ops.len(),
        emitted_x: harness.builder.counted_kind_ops[OperationType::X as usize],
        emitted_toffoli: harness.builder.counted_kind_ops[OperationType::CCX as usize]
            + harness.builder.counted_kind_ops[OperationType::CCZ as usize],
    }
}

/// Prove that the signed high lane can be replaced by a clean borrowed
/// underflow witness. The Cuccaro carry is restored after subtraction, so the
/// same lane can hold the enable predicate before the inverse addition.
#[doc(hidden)]
#[must_use]
pub fn exhaustive_borrowed_rotated_underflow_check(
) -> BorrowedRotatedUnderflowProofReport {
    assert!(
        std::env::var_os(SUB800_BORROWED_ROTATED_UNDERFLOW_FLAG).is_none(),
        "the borrowed rotated-underflow feature must default off"
    );
    let _environment = RawBitLengthProofEnvironment::capture();
    configure_paired_bitlength_source_complement_proof(true, true);
    std::env::set_var(SUB800_MIXED_BOUNDARY_SCRATCH_EXTENSION_FLAG, "1");

    let mut configurations_checked = 0usize;
    let mut basis_states_checked = 0usize;
    let mut scalar_equivalence_checks = 0usize;
    let mut simulator_equivalence_checks = 0usize;
    let mut phase_clean_checks = 0usize;
    let mut inverse_pair_checks = 0usize;
    let mut control_off_checks = 0usize;
    let mut saturation_underflow_checks = 0usize;
    let mut nonnegative_difference_checks = 0usize;
    let mut preserved_input_checks = 0usize;
    let mut borrowed_scratch_clean_checks = 0usize;
    let mut ancilla_clean_checks = 0usize;

    for output_width in [2usize, 3] {
        let maximum_source_width = ((1usize << output_width) - 1).min(4);
        for source_width in 2..=maximum_source_width {
            for mixed_boundary in [false, true] {
                let boundary_width = output_width - usize::from(mixed_boundary);
                for paired_source in [false, true] {
                    std::env::remove_var(SUB800_BORROWED_ROTATED_UNDERFLOW_FLAG);
                    let baseline = build_paired_bitlength_source_complement_harness(
                        output_width,
                        source_width,
                        3,
                        SaturatingDifferenceBoundaryRoute::Inplace,
                        mixed_boundary,
                        Some(paired_source),
                    );
                    std::env::set_var(SUB800_BORROWED_ROTATED_UNDERFLOW_FLAG, "1");
                    let candidate = build_paired_bitlength_source_complement_harness(
                        output_width,
                        source_width,
                        3,
                        SaturatingDifferenceBoundaryRoute::Inplace,
                        mixed_boundary,
                        Some(paired_source),
                    );
                    assert_eq!(baseline.data_ids, candidate.data_ids);
                    assert_eq!(baseline.external_mask, candidate.external_mask);
                    let (simulator_cases, simulator_phases) =
                        verify_inplace_rotated_boundary_simulator_equivalence(
                            &baseline,
                            &candidate,
                        );
                    simulator_equivalence_checks += simulator_cases;
                    phase_clean_checks += simulator_phases;

                    let data_states = 1u64 << baseline.data_ids.len();
                    let boundary_mask = (1u64 << boundary_width) - 1;
                    let source_mask = (1u64 << source_width) - 1;
                    for value in 0..data_states {
                        let input = inplace_rotated_boundary_input(&baseline, value);
                        let baseline_output = apply_scalar(&baseline.builder.ops, input);
                        let candidate_output = apply_scalar(&candidate.builder.ops, input);
                        assert_eq!(
                            candidate_output, baseline_output,
                            "borrowed underflow mismatch output_width={output_width} \
                             source_width={source_width} mixed={mixed_boundary} \
                             paired={paired_source} value={value}"
                        );
                        assert_inplace_rotated_boundary_clean(
                            &baseline,
                            baseline_output,
                            "owned signed lane",
                        );
                        assert_inplace_rotated_boundary_clean(
                            &candidate,
                            candidate_output,
                            "borrowed underflow lane",
                        );
                        assert_eq!(
                            baseline_output & baseline.preserved_mask,
                            input & baseline.preserved_mask
                        );
                        assert_eq!(
                            candidate_output & candidate.preserved_mask,
                            input & candidate.preserved_mask
                        );
                        assert_eq!(
                            apply_scalar(&baseline.builder.ops, baseline_output),
                            input
                        );
                        assert_eq!(
                            apply_scalar(&candidate.builder.ops, candidate_output),
                            input
                        );
                        if input & candidate.control_mask == 0 {
                            assert_eq!(
                                candidate_output & candidate.output_mask,
                                input & candidate.output_mask
                            );
                            control_off_checks += 1;
                        }

                        let boundary = (value >> 1) & boundary_mask;
                        let source =
                            (value >> (1 + boundary_width)) & source_mask;
                        if bit_length_usize(source as usize) < boundary as usize {
                            saturation_underflow_checks += 1;
                        } else {
                            nonnegative_difference_checks += 1;
                        }
                        preserved_input_checks += 2;
                        borrowed_scratch_clean_checks += 2;
                        ancilla_clean_checks += 2;
                        inverse_pair_checks += 2;
                        scalar_equivalence_checks += 1;
                        basis_states_checked += 1;
                    }
                    configurations_checked += 1;
                }
            }
        }
    }

    configure_paired_bitlength_source_complement_proof(true, true);
    std::env::set_var(SUB800_MIXED_BOUNDARY_SCRATCH_EXTENSION_FLAG, "1");
    std::env::remove_var(SUB800_BORROWED_ROTATED_UNDERFLOW_FLAG);
    let same_boundary_baseline = build_paired_bitlength_source_complement_harness(
        REFERENCE_LENGTH_WIDTH,
        259,
        7,
        SaturatingDifferenceBoundaryRoute::Inplace,
        false,
        Some(true),
    );
    let mixed_boundary_baseline = build_paired_bitlength_source_complement_harness(
        REFERENCE_LENGTH_WIDTH,
        259,
        7,
        SaturatingDifferenceBoundaryRoute::Inplace,
        true,
        Some(true),
    );
    std::env::set_var(SUB800_BORROWED_ROTATED_UNDERFLOW_FLAG, "1");
    let same_boundary_candidate = build_paired_bitlength_source_complement_harness(
        REFERENCE_LENGTH_WIDTH,
        259,
        7,
        SaturatingDifferenceBoundaryRoute::Inplace,
        false,
        Some(true),
    );
    let mixed_boundary_candidate = build_paired_bitlength_source_complement_harness(
        REFERENCE_LENGTH_WIDTH,
        259,
        7,
        SaturatingDifferenceBoundaryRoute::Inplace,
        true,
        Some(true),
    );

    let same_boundary_baseline =
        inplace_rotated_boundary_local_resources(&same_boundary_baseline);
    let same_boundary_candidate =
        inplace_rotated_boundary_local_resources(&same_boundary_candidate);
    let mixed_boundary_baseline =
        inplace_rotated_boundary_local_resources(&mixed_boundary_baseline);
    let mixed_boundary_candidate =
        inplace_rotated_boundary_local_resources(&mixed_boundary_candidate);
    assert_eq!(
        same_boundary_candidate.peak_qubits + 1,
        same_boundary_baseline.peak_qubits
    );
    assert_eq!(
        mixed_boundary_candidate.peak_qubits + 1,
        mixed_boundary_baseline.peak_qubits
    );
    assert!(
        same_boundary_candidate.emitted_toffoli
            <= same_boundary_baseline.emitted_toffoli
    );
    assert!(
        mixed_boundary_candidate.emitted_toffoli
            <= mixed_boundary_baseline.emitted_toffoli
    );

    BorrowedRotatedUnderflowProofReport {
        configurations_checked,
        boundary_forms_checked: 2,
        paired_source_modes_checked: 2,
        basis_states_checked,
        scalar_equivalence_checks,
        simulator_equivalence_checks,
        phase_clean_checks,
        inverse_pair_checks,
        control_off_checks,
        saturation_underflow_checks,
        nonnegative_difference_checks,
        preserved_input_checks,
        borrowed_scratch_clean_checks,
        ancilla_clean_checks,
        same_boundary_baseline,
        same_boundary_candidate,
        mixed_boundary_baseline,
        mixed_boundary_candidate,
    }
}

/// Prove the mixed-width identity
///
/// `(256h+l)-b = 256(h XOR borrow)+(l-b mod 256)`
///
/// together with saturation on `borrow AND NOT h`. The high bit `h` is
/// reconstructed from the top source lanes and held in restored scratch.
#[doc(hidden)]
#[must_use]
pub fn exhaustive_split_mixed_rotated_length_check(
) -> SplitMixedRotatedLengthProofReport {
    assert!(
        std::env::var_os(SUB800_SPLIT_MIXED_ROTATED_LENGTH_FLAG).is_none(),
        "the split mixed rotated-length feature must default off"
    );
    let _environment = RawBitLengthProofEnvironment::capture();
    configure_paired_bitlength_source_complement_proof(true, true);
    std::env::set_var(SUB800_MIXED_BOUNDARY_SCRATCH_EXTENSION_FLAG, "1");
    std::env::set_var(SUB800_BORROWED_ROTATED_UNDERFLOW_FLAG, "1");

    let mut configurations_checked = 0usize;
    let mut basis_states_checked = 0usize;
    let mut scalar_equivalence_checks = 0usize;
    let mut simulator_equivalence_checks = 0usize;
    let mut phase_clean_checks = 0usize;
    let mut inverse_pair_checks = 0usize;
    let mut control_off_checks = 0usize;
    let mut high_indicator_checks = 0usize;
    let mut saturation_underflow_checks = 0usize;
    let mut nonnegative_difference_checks = 0usize;
    let mut preserved_input_checks = 0usize;
    let mut borrowed_scratch_clean_checks = 0usize;
    let mut ancilla_clean_checks = 0usize;

    for output_width in [2usize, 3] {
        let maximum_source_width = ((1usize << output_width) - 1).min(4);
        let boundary_width = output_width - 1;
        for source_width in 2..=maximum_source_width {
            for paired_source in [false, true] {
                std::env::remove_var(SUB800_SPLIT_MIXED_ROTATED_LENGTH_FLAG);
                let baseline = build_paired_bitlength_source_complement_harness(
                    output_width,
                    source_width,
                    3,
                    SaturatingDifferenceBoundaryRoute::Inplace,
                    true,
                    Some(paired_source),
                );
                std::env::set_var(SUB800_SPLIT_MIXED_ROTATED_LENGTH_FLAG, "1");
                let candidate = build_paired_bitlength_source_complement_harness(
                    output_width,
                    source_width,
                    3,
                    SaturatingDifferenceBoundaryRoute::Inplace,
                    true,
                    Some(paired_source),
                );
                assert_eq!(baseline.data_ids, candidate.data_ids);
                assert_eq!(baseline.external_mask, candidate.external_mask);
                let (simulator_cases, simulator_phases) =
                    verify_inplace_rotated_boundary_simulator_equivalence(
                        &baseline,
                        &candidate,
                    );
                simulator_equivalence_checks += simulator_cases;
                phase_clean_checks += simulator_phases;

                let data_states = 1u64 << baseline.data_ids.len();
                let boundary_mask = (1u64 << boundary_width) - 1;
                let source_mask = (1u64 << source_width) - 1;
                for value in 0..data_states {
                    let input = inplace_rotated_boundary_input(&baseline, value);
                    let baseline_output = apply_scalar(&baseline.builder.ops, input);
                    let candidate_output = apply_scalar(&candidate.builder.ops, input);
                    assert_eq!(
                        candidate_output, baseline_output,
                        "split mixed length mismatch output_width={output_width} \
                         source_width={source_width} paired={paired_source} value={value}"
                    );
                    assert_inplace_rotated_boundary_clean(
                        &baseline,
                        baseline_output,
                        "nine-lane mixed length",
                    );
                    assert_inplace_rotated_boundary_clean(
                        &candidate,
                        candidate_output,
                        "split mixed length",
                    );
                    assert_eq!(
                        baseline_output & baseline.preserved_mask,
                        input & baseline.preserved_mask
                    );
                    assert_eq!(
                        candidate_output & candidate.preserved_mask,
                        input & candidate.preserved_mask
                    );
                    assert_eq!(
                        apply_scalar(&baseline.builder.ops, baseline_output),
                        input
                    );
                    assert_eq!(
                        apply_scalar(&candidate.builder.ops, candidate_output),
                        input
                    );
                    if input & candidate.control_mask == 0 {
                        assert_eq!(
                            candidate_output & candidate.output_mask,
                            input & candidate.output_mask
                        );
                        control_off_checks += 1;
                    }

                    let boundary = (value >> 1) & boundary_mask;
                    let source = (value >> (1 + boundary_width)) & source_mask;
                    let bit_length = bit_length_usize(source as usize);
                    if bit_length >= (1usize << (output_width - 1)) {
                        high_indicator_checks += 1;
                    }
                    if bit_length < boundary as usize {
                        saturation_underflow_checks += 1;
                    } else {
                        nonnegative_difference_checks += 1;
                    }
                    preserved_input_checks += 2;
                    borrowed_scratch_clean_checks += 2;
                    ancilla_clean_checks += 2;
                    inverse_pair_checks += 2;
                    scalar_equivalence_checks += 1;
                    basis_states_checked += 1;
                }
                configurations_checked += 1;
            }
        }
    }

    configure_paired_bitlength_source_complement_proof(true, true);
    std::env::set_var(SUB800_MIXED_BOUNDARY_SCRATCH_EXTENSION_FLAG, "1");
    std::env::set_var(SUB800_BORROWED_ROTATED_UNDERFLOW_FLAG, "1");
    std::env::remove_var(SUB800_SPLIT_MIXED_ROTATED_LENGTH_FLAG);
    let baseline_local = build_paired_bitlength_source_complement_harness(
        REFERENCE_LENGTH_WIDTH,
        259,
        7,
        SaturatingDifferenceBoundaryRoute::Inplace,
        true,
        Some(true),
    );
    std::env::set_var(SUB800_SPLIT_MIXED_ROTATED_LENGTH_FLAG, "1");
    let candidate_local = build_paired_bitlength_source_complement_harness(
        REFERENCE_LENGTH_WIDTH,
        259,
        7,
        SaturatingDifferenceBoundaryRoute::Inplace,
        true,
        Some(true),
    );
    let baseline_local = inplace_rotated_boundary_local_resources(&baseline_local);
    let candidate_local = inplace_rotated_boundary_local_resources(&candidate_local);
    assert_eq!(candidate_local.peak_qubits + 1, baseline_local.peak_qubits);

    SplitMixedRotatedLengthProofReport {
        configurations_checked,
        paired_source_modes_checked: 2,
        basis_states_checked,
        scalar_equivalence_checks,
        simulator_equivalence_checks,
        phase_clean_checks,
        inverse_pair_checks,
        control_off_checks,
        high_indicator_checks,
        saturation_underflow_checks,
        nonnegative_difference_checks,
        preserved_input_checks,
        borrowed_scratch_clean_checks,
        ancilla_clean_checks,
        baseline_local,
        candidate_local,
    }
}

/// Prove the full one-bit high-stage subtraction used when the boundary and
/// output have equal width. The low Cuccaro borrow `u`, source high bit `h`,
/// and boundary high bit `g` determine the final underflow through
/// `g XOR u XOR gu XOR hg XOR hu`.
#[doc(hidden)]
#[must_use]
pub fn exhaustive_split_same_rotated_length_check(
) -> SplitSameRotatedLengthProofReport {
    assert!(
        std::env::var_os(SUB800_SPLIT_SAME_ROTATED_LENGTH_FLAG).is_none(),
        "the split same-width rotated-length feature must default off"
    );
    let _environment = RawBitLengthProofEnvironment::capture();
    configure_paired_bitlength_source_complement_proof(true, true);
    std::env::set_var(SUB800_MIXED_BOUNDARY_SCRATCH_EXTENSION_FLAG, "1");
    std::env::set_var(SUB800_BORROWED_ROTATED_UNDERFLOW_FLAG, "1");

    let mut configurations_checked = 0usize;
    let mut basis_states_checked = 0usize;
    let mut scalar_equivalence_checks = 0usize;
    let mut simulator_equivalence_checks = 0usize;
    let mut phase_clean_checks = 0usize;
    let mut inverse_pair_checks = 0usize;
    let mut control_off_checks = 0usize;
    let mut high_indicator_checks = 0usize;
    let mut saturation_underflow_checks = 0usize;
    let mut nonnegative_difference_checks = 0usize;
    let mut preserved_input_checks = 0usize;
    let mut borrowed_scratch_clean_checks = 0usize;
    let mut ancilla_clean_checks = 0usize;

    for output_width in [2usize, 3] {
        let maximum_source_width = ((1usize << output_width) - 1).min(4);
        let boundary_width = output_width;
        for source_width in 2..=maximum_source_width {
            for paired_source in [false, true] {
                std::env::remove_var(SUB800_SPLIT_SAME_ROTATED_LENGTH_FLAG);
                let baseline = build_paired_bitlength_source_complement_harness(
                    output_width,
                    source_width,
                    3,
                    SaturatingDifferenceBoundaryRoute::Inplace,
                    false,
                    Some(paired_source),
                );
                std::env::set_var(SUB800_SPLIT_SAME_ROTATED_LENGTH_FLAG, "1");
                let candidate = build_paired_bitlength_source_complement_harness(
                    output_width,
                    source_width,
                    3,
                    SaturatingDifferenceBoundaryRoute::Inplace,
                    false,
                    Some(paired_source),
                );
                assert_eq!(baseline.data_ids, candidate.data_ids);
                assert_eq!(baseline.external_mask, candidate.external_mask);
                let (simulator_cases, simulator_phases) =
                    verify_inplace_rotated_boundary_simulator_equivalence(
                        &baseline,
                        &candidate,
                    );
                simulator_equivalence_checks += simulator_cases;
                phase_clean_checks += simulator_phases;

                let data_states = 1u64 << baseline.data_ids.len();
                let boundary_mask = (1u64 << boundary_width) - 1;
                let source_mask = (1u64 << source_width) - 1;
                for value in 0..data_states {
                    let input = inplace_rotated_boundary_input(&baseline, value);
                    let baseline_output = apply_scalar(&baseline.builder.ops, input);
                    let candidate_output = apply_scalar(&candidate.builder.ops, input);
                    assert_eq!(
                        candidate_output, baseline_output,
                        "split same-width length mismatch output_width={output_width} \
                         source_width={source_width} paired={paired_source} value={value}"
                    );
                    assert_inplace_rotated_boundary_clean(
                        &baseline,
                        baseline_output,
                        "nine-lane same-width length",
                    );
                    assert_inplace_rotated_boundary_clean(
                        &candidate,
                        candidate_output,
                        "split same-width length",
                    );
                    assert_eq!(
                        baseline_output & baseline.preserved_mask,
                        input & baseline.preserved_mask
                    );
                    assert_eq!(
                        candidate_output & candidate.preserved_mask,
                        input & candidate.preserved_mask
                    );
                    assert_eq!(
                        apply_scalar(&baseline.builder.ops, baseline_output),
                        input
                    );
                    assert_eq!(
                        apply_scalar(&candidate.builder.ops, candidate_output),
                        input
                    );
                    if input & candidate.control_mask == 0 {
                        assert_eq!(
                            candidate_output & candidate.output_mask,
                            input & candidate.output_mask
                        );
                        control_off_checks += 1;
                    }

                    let boundary = (value >> 1) & boundary_mask;
                    let source = (value >> (1 + boundary_width)) & source_mask;
                    let bit_length = bit_length_usize(source as usize);
                    if bit_length >= (1usize << (output_width - 1)) {
                        high_indicator_checks += 1;
                    }
                    if bit_length < boundary as usize {
                        saturation_underflow_checks += 1;
                    } else {
                        nonnegative_difference_checks += 1;
                    }
                    preserved_input_checks += 2;
                    borrowed_scratch_clean_checks += 2;
                    ancilla_clean_checks += 2;
                    inverse_pair_checks += 2;
                    scalar_equivalence_checks += 1;
                    basis_states_checked += 1;
                }
                configurations_checked += 1;
            }
        }
    }

    configure_paired_bitlength_source_complement_proof(true, true);
    std::env::set_var(SUB800_MIXED_BOUNDARY_SCRATCH_EXTENSION_FLAG, "1");
    std::env::set_var(SUB800_BORROWED_ROTATED_UNDERFLOW_FLAG, "1");
    std::env::remove_var(SUB800_SPLIT_SAME_ROTATED_LENGTH_FLAG);
    let baseline_local = build_paired_bitlength_source_complement_harness(
        REFERENCE_LENGTH_WIDTH,
        259,
        7,
        SaturatingDifferenceBoundaryRoute::Inplace,
        false,
        Some(true),
    );
    std::env::set_var(SUB800_SPLIT_SAME_ROTATED_LENGTH_FLAG, "1");
    let candidate_local = build_paired_bitlength_source_complement_harness(
        REFERENCE_LENGTH_WIDTH,
        259,
        7,
        SaturatingDifferenceBoundaryRoute::Inplace,
        false,
        Some(true),
    );
    let baseline_local = inplace_rotated_boundary_local_resources(&baseline_local);
    let candidate_local = inplace_rotated_boundary_local_resources(&candidate_local);
    assert_eq!(candidate_local.peak_qubits + 1, baseline_local.peak_qubits);

    SplitSameRotatedLengthProofReport {
        configurations_checked,
        paired_source_modes_checked: 2,
        basis_states_checked,
        scalar_equivalence_checks,
        simulator_equivalence_checks,
        phase_clean_checks,
        inverse_pair_checks,
        control_off_checks,
        high_indicator_checks,
        saturation_underflow_checks,
        nonnegative_difference_checks,
        preserved_input_checks,
        borrowed_scratch_clean_checks,
        ancilla_clean_checks,
        baseline_local,
        candidate_local,
    }
}

fn read_paired_bitlength_register(state: u64, ids: &[u32]) -> u64 {
    ids.iter().enumerate().fold(0u64, |value, (bit, id)| {
        value | (((state >> id) & 1) << bit)
    })
}

struct SplitTwoHighStageHarness {
    builder: B,
    data_ids: Vec<u32>,
    control_id: u32,
    boundary_ids: Vec<u32>,
    length_ids: Vec<u32>,
    high7_id: u32,
    high8_id: u32,
    output_ids: Vec<u32>,
    preserved_mask: u64,
    output_mask: u64,
    clean_scratch_mask: u64,
}

fn build_split_two_high_stage_harness(mixed_boundary: bool) -> SplitTwoHighStageHarness {
    const LOW_WIDTH: usize = 2;
    const OUTPUT_WIDTH: usize = LOW_WIDTH + 2;

    let mut circ = Circuit::new();
    let control = circ.alloc_qreg("sub800.two-high-proof.control");
    let boundary = circ.alloc_qreg_bits(
        "sub800.two-high-proof.boundary",
        OUTPUT_WIDTH - usize::from(mixed_boundary),
    );
    let length = circ.alloc_qreg_bits("sub800.two-high-proof.length", LOW_WIDTH);
    let output = circ.alloc_qreg_bits("sub800.two-high-proof.output", OUTPUT_WIDTH);
    let scratch = circ.alloc_qreg_bits("sub800.two-high-proof.scratch", 7);
    let scratch_refs: Vec<&QReg> = scratch.iter().collect();
    if mixed_boundary {
        controlled_xor_saturating_difference_inplace_mixed_boundary(
            &mut circ,
            &control,
            &boundary,
            &[],
            &length,
            &output,
            &scratch_refs,
        );
    } else {
        controlled_xor_saturating_difference_inplace_boundary(
            &mut circ,
            &control,
            &boundary,
            &[],
            &length,
            &output,
            &scratch_refs,
        );
    }
    drop(scratch_refs);

    let data_ids: Vec<u32> = std::iter::once(&control)
        .chain(&boundary)
        .chain(&length)
        .chain(std::iter::once(&scratch[5]))
        .chain(std::iter::once(&scratch[6]))
        .chain(&output)
        .map(QReg::id)
        .collect();
    let preserved_mask = qreg_mask(
        std::iter::once(&control)
            .chain(&boundary)
            .chain(&length)
            .chain(std::iter::once(&scratch[5]))
            .chain(std::iter::once(&scratch[6])),
    );
    let output_mask = qreg_mask(&output);
    let clean_scratch_mask = qreg_mask(&scratch[..5]);
    SplitTwoHighStageHarness {
        builder: circ.into_builder(),
        data_ids,
        control_id: control.id(),
        boundary_ids: boundary.iter().map(QReg::id).collect(),
        length_ids: length.iter().map(QReg::id).collect(),
        high7_id: scratch[5].id(),
        high8_id: scratch[6].id(),
        output_ids: output.iter().map(QReg::id).collect(),
        preserved_mask,
        output_mask,
        clean_scratch_mask,
    }
}

fn split_two_high_stage_input(harness: &SplitTwoHighStageHarness, value: u64) -> u64 {
    harness
        .data_ids
        .iter()
        .enumerate()
        .fold(0u64, |state, (bit, id)| {
            state | (((value >> bit) & 1) << id)
        })
}

fn full_subtraction_borrow_anf(x: usize, y: usize, borrow: usize) -> usize {
    y ^ borrow ^ (y & borrow) ^ (x & y) ^ (x & borrow)
}

struct SplitTwoHighIndicatorHarness {
    builder: B,
    source_ids: Vec<u32>,
    high7_id: u32,
    high8_id: u32,
    clean_scratch_ids: Vec<u32>,
}

fn build_split_two_high_indicator_harness(
    source_is_complemented: bool,
) -> SplitTwoHighIndicatorHarness {
    const SOURCE_WIDTH: usize = 259;

    let mut circ = Circuit::new();
    let source = circ.alloc_input_qreg_bits("sub800.two-high-indicator.source", SOURCE_WIDTH);
    let scratch = circ.alloc_qreg_bits("sub800.two-high-indicator.scratch", 7);
    let source_refs: Vec<&QReg> = source.iter().collect();
    let scratch_refs: Vec<&QReg> = scratch.iter().collect();
    toggle_split_bit_length_two_high(
        &mut circ,
        &source_refs,
        &scratch[5],
        &scratch[6],
        &scratch_refs,
        source_is_complemented,
    );
    drop(source_refs);
    drop(scratch_refs);

    SplitTwoHighIndicatorHarness {
        builder: circ.into_builder(),
        source_ids: source.iter().map(QReg::id).collect(),
        high7_id: scratch[5].id(),
        high8_id: scratch[6].id(),
        clean_scratch_ids: scratch[..5].iter().map(QReg::id).collect(),
    }
}

fn verify_split_two_high_indicator_harness(
    harness: &SplitTwoHighIndicatorHarness,
    source_is_complemented: bool,
) -> (usize, usize, usize, usize) {
    use crate::circuit::OperationType;

    assert!(harness.builder.ops.iter().all(|operation| matches!(
        operation.kind,
        OperationType::X | OperationType::CX | OperationType::CCX
    )));
    let total_qubits = harness.builder.next_qubit as usize;
    let mut basis_states = 0usize;
    let mut inverse_checks = 0usize;
    let mut source_restore_checks = 0usize;
    let mut scratch_clean_checks = 0usize;
    for bit_length in 0usize..=259 {
        for initial_highs in 0usize..4 {
            let mut input = vec![false; total_qubits];
            for (bit, id) in harness.source_ids.iter().copied().enumerate() {
                let logical = bit_length != 0 && bit == bit_length - 1;
                input[id as usize] = logical ^ source_is_complemented;
            }
            input[harness.high7_id as usize] = (initial_highs & 1) != 0;
            input[harness.high8_id as usize] = (initial_highs & 2) != 0;

            let output = apply_basis_vector(&harness.builder.ops, input.clone());
            let expected_high7 = ((bit_length >> 7) & 1) != 0;
            let expected_high8 = ((bit_length >> 8) & 1) != 0;
            assert_eq!(
                output[harness.high7_id as usize],
                input[harness.high7_id as usize] ^ expected_high7
            );
            assert_eq!(
                output[harness.high8_id as usize],
                input[harness.high8_id as usize] ^ expected_high8
            );
            assert!(harness
                .source_ids
                .iter()
                .all(|id| output[*id as usize] == input[*id as usize]));
            assert!(harness
                .clean_scratch_ids
                .iter()
                .all(|id| !output[*id as usize]));
            assert_eq!(apply_basis_vector(&harness.builder.ops, output), input);
            basis_states += 1;
            inverse_checks += 1;
            source_restore_checks += 1;
            scratch_clean_checks += 1;
        }
    }
    (
        basis_states,
        inverse_checks,
        source_restore_checks,
        scratch_clean_checks,
    )
}

fn verify_split_two_high_stage(
    harness: &SplitTwoHighStageHarness,
) -> (usize, usize, usize, usize, usize) {
    use crate::circuit::OperationType;

    assert!(harness.builder.ops.iter().all(|op| matches!(
        op.kind,
        OperationType::X | OperationType::CX | OperationType::CCX
    )));
    let states = 1u64 << harness.data_ids.len();
    let mut scalar_equivalence_checks = 0usize;
    let mut inverse_pair_checks = 0usize;
    let mut control_off_checks = 0usize;
    let mut phase_clean_checks = 0usize;
    let mut scratch_restore_checks = 0usize;
    for value in 0..states {
        let input = split_two_high_stage_input(harness, value);
        let output = apply_scalar(&harness.builder.ops, input);
        let control = ((input >> harness.control_id) & 1) as usize;
        let boundary = read_paired_bitlength_register(input, &harness.boundary_ids) as usize;
        let low = read_paired_bitlength_register(input, &harness.length_ids) as usize;
        let high7 = ((input >> harness.high7_id) & 1) as usize;
        let high8 = ((input >> harness.high8_id) & 1) as usize;
        let length = low | (high7 << harness.length_ids.len())
            | (high8 << (harness.length_ids.len() + 1));
        let initial_output = read_paired_bitlength_register(input, &harness.output_ids) as usize;
        let expected_xor = if control == 1 && length >= boundary {
            length - boundary
        } else {
            0
        };
        let expected_output = initial_output ^ expected_xor;
        assert_eq!(
            read_paired_bitlength_register(output, &harness.output_ids) as usize,
            expected_output
        );
        assert_eq!(output & harness.preserved_mask, input & harness.preserved_mask);
        assert_eq!(output & harness.clean_scratch_mask, 0);
        assert_eq!(apply_scalar(&harness.builder.ops, output), input);
        if control == 0 {
            assert_eq!(output & harness.output_mask, input & harness.output_mask);
            control_off_checks += 1;
        }
        scalar_equivalence_checks += 1;
        inverse_pair_checks += 1;
        phase_clean_checks += 1;
        scratch_restore_checks += 1;
    }
    (
        scalar_equivalence_checks,
        inverse_pair_checks,
        control_off_checks,
        phase_clean_checks,
        scratch_restore_checks,
    )
}

/// Prove the two-bit high-stage decomposition used by the Q840 candidate.
/// The production source support has every attainable bit length 0..259.
#[doc(hidden)]
#[must_use]
pub fn exhaustive_split_two_high_rotated_length_check(
) -> SplitTwoHighRotatedLengthProofReport {
    assert!(
        std::env::var_os(SUB800_SPLIT_TWO_HIGH_ROTATED_LENGTH_FLAG).is_none(),
        "the split-two-high rotated-length feature must default off"
    );
    let _environment = RawBitLengthProofEnvironment::capture();

    let mut borrow_truth_table_cases = 0usize;
    for x in 0usize..=1 {
        for y in 0usize..=1 {
            for borrow in 0usize..=1 {
                let expected = usize::from((2 * x) < (2 * y + borrow));
                assert_eq!(full_subtraction_borrow_anf(x, y, borrow), expected);
                borrow_truth_table_cases += 1;
            }
        }
    }

    let mut bit_lengths_checked = 0usize;
    let mut same_arithmetic_cases = 0usize;
    let mut mixed_arithmetic_cases = 0usize;
    for bit_length in 0usize..=259 {
        let low = bit_length & 0x7f;
        let high7 = (bit_length >> 7) & 1;
        let high8 = (bit_length >> 8) & 1;
        let z128 = usize::from(bit_length < 128);
        let z256 = usize::from(bit_length < 256);
        assert_eq!(high7, z128 ^ z256);
        assert_eq!(high8, 1 ^ z256);
        bit_lengths_checked += 1;

        for boundary in 0usize..512 {
            let boundary_low = boundary & 0x7f;
            let boundary_high7 = (boundary >> 7) & 1;
            let boundary_high8 = (boundary >> 8) & 1;
            let low_borrow = usize::from(low < boundary_low);
            let high_borrow =
                full_subtraction_borrow_anf(high7, boundary_high7, low_borrow);
            let underflow =
                full_subtraction_borrow_anf(high8, boundary_high8, high_borrow);
            let difference = ((low + 128 - boundary_low) & 0x7f)
                | ((high7 ^ boundary_high7 ^ low_borrow) << 7)
                | ((high8 ^ boundary_high8 ^ high_borrow) << 8);
            assert_eq!(underflow, usize::from(bit_length < boundary));
            assert_eq!(difference, bit_length.wrapping_sub(boundary) & 0x1ff);
            same_arithmetic_cases += 1;
        }

        for boundary in 0usize..256 {
            let boundary_low = boundary & 0x7f;
            let boundary_high7 = (boundary >> 7) & 1;
            let low_borrow = usize::from(low < boundary_low);
            let high_borrow =
                full_subtraction_borrow_anf(high7, boundary_high7, low_borrow);
            let underflow = full_subtraction_borrow_anf(high8, 0, high_borrow);
            let difference = ((low + 128 - boundary_low) & 0x7f)
                | ((high7 ^ boundary_high7 ^ low_borrow) << 7)
                | ((high8 ^ high_borrow) << 8);
            assert_eq!(underflow, usize::from(bit_length < boundary));
            assert_eq!(difference, bit_length.wrapping_sub(boundary) & 0x1ff);
            mixed_arithmetic_cases += 1;
        }
    }

    let same_stage = build_split_two_high_stage_harness(false);
    let mixed_stage = build_split_two_high_stage_harness(true);
    let high_indicator = build_split_two_high_indicator_harness(false);
    let complemented_high_indicator = build_split_two_high_indicator_harness(true);
    let same_circuit_basis_states = 1usize << same_stage.data_ids.len();
    let mixed_circuit_basis_states = 1usize << mixed_stage.data_ids.len();
    let same_checks = verify_split_two_high_stage(&same_stage);
    let mixed_checks = verify_split_two_high_stage(&mixed_stage);
    let high_checks = verify_split_two_high_indicator_harness(&high_indicator, false);
    let complemented_high_checks =
        verify_split_two_high_indicator_harness(&complemented_high_indicator, true);

    configure_paired_bitlength_source_complement_proof(true, true);
    std::env::set_var(SUB800_MIXED_BOUNDARY_SCRATCH_EXTENSION_FLAG, "1");
    std::env::set_var(SUB800_BORROWED_ROTATED_UNDERFLOW_FLAG, "1");
    std::env::set_var(SUB800_SPLIT_MIXED_ROTATED_LENGTH_FLAG, "1");
    std::env::set_var(SUB800_SPLIT_SAME_ROTATED_LENGTH_FLAG, "1");
    std::env::remove_var(SUB800_SPLIT_TWO_HIGH_ROTATED_LENGTH_FLAG);
    let same_boundary_baseline = build_paired_bitlength_source_complement_harness(
        REFERENCE_LENGTH_WIDTH,
        259,
        7,
        SaturatingDifferenceBoundaryRoute::Inplace,
        false,
        Some(true),
    );
    let mixed_boundary_baseline = build_paired_bitlength_source_complement_harness(
        REFERENCE_LENGTH_WIDTH,
        259,
        7,
        SaturatingDifferenceBoundaryRoute::Inplace,
        true,
        Some(true),
    );
    std::env::set_var(SUB800_SPLIT_TWO_HIGH_ROTATED_LENGTH_FLAG, "1");
    let same_boundary_candidate = build_paired_bitlength_source_complement_harness(
        REFERENCE_LENGTH_WIDTH,
        259,
        7,
        SaturatingDifferenceBoundaryRoute::Inplace,
        false,
        Some(true),
    );
    let mixed_boundary_candidate = build_paired_bitlength_source_complement_harness(
        REFERENCE_LENGTH_WIDTH,
        259,
        7,
        SaturatingDifferenceBoundaryRoute::Inplace,
        true,
        Some(true),
    );
    let same_boundary_baseline =
        inplace_rotated_boundary_local_resources(&same_boundary_baseline);
    let same_boundary_candidate =
        inplace_rotated_boundary_local_resources(&same_boundary_candidate);
    let mixed_boundary_baseline =
        inplace_rotated_boundary_local_resources(&mixed_boundary_baseline);
    let mixed_boundary_candidate =
        inplace_rotated_boundary_local_resources(&mixed_boundary_candidate);
    assert_eq!(same_boundary_candidate.peak_qubits + 1, same_boundary_baseline.peak_qubits);
    assert_eq!(
        mixed_boundary_candidate.peak_qubits + 1,
        mixed_boundary_baseline.peak_qubits
    );

    SplitTwoHighRotatedLengthProofReport {
        borrow_truth_table_cases,
        bit_lengths_checked,
        high_indicator_basis_states: high_checks.0 + complemented_high_checks.0,
        high_indicator_inverse_checks: high_checks.1 + complemented_high_checks.1,
        high_indicator_source_restore_checks: high_checks.2 + complemented_high_checks.2,
        high_indicator_scratch_clean_checks: high_checks.3 + complemented_high_checks.3,
        same_arithmetic_cases,
        mixed_arithmetic_cases,
        same_circuit_basis_states,
        mixed_circuit_basis_states,
        scalar_equivalence_checks: same_checks.0 + mixed_checks.0,
        inverse_pair_checks: same_checks.1 + mixed_checks.1,
        control_off_checks: same_checks.2 + mixed_checks.2,
        phase_clean_checks: same_checks.3 + mixed_checks.3,
        scratch_restore_checks: same_checks.4 + mixed_checks.4,
        same_boundary_baseline,
        same_boundary_candidate,
        mixed_boundary_baseline,
        mixed_boundary_candidate,
    }
}

struct SplitThreeHighStageHarness {
    builder: B,
    data_ids: Vec<u32>,
    control_id: u32,
    boundary_ids: Vec<u32>,
    length_ids: Vec<u32>,
    high6_id: u32,
    high7_id: u32,
    high8_id: u32,
    output_ids: Vec<u32>,
    dirty_id: u32,
    preserved_mask: u64,
    output_mask: u64,
    clean_scratch_mask: u64,
}

fn build_split_three_high_stage_harness(mixed_boundary: bool) -> SplitThreeHighStageHarness {
    const LOW_WIDTH: usize = 2;
    const OUTPUT_WIDTH: usize = LOW_WIDTH + 3;

    let mut circ = Circuit::new();
    let control = circ.alloc_qreg("sub800.three-high-proof.control");
    let boundary = circ.alloc_qreg_bits(
        "sub800.three-high-proof.boundary",
        OUTPUT_WIDTH - usize::from(mixed_boundary),
    );
    let length = circ.alloc_qreg_bits("sub800.three-high-proof.length", LOW_WIDTH);
    let output = circ.alloc_qreg_bits("sub800.three-high-proof.output", OUTPUT_WIDTH);
    let scratch = circ.alloc_qreg_bits("sub800.three-high-proof.scratch", 8);
    let dirty = circ.alloc_qreg("sub800.three-high-proof.dirty");
    let source = [&dirty];
    let scratch_refs: Vec<&QReg> = scratch.iter().collect();
    if mixed_boundary {
        controlled_xor_saturating_difference_inplace_mixed_boundary(
            &mut circ,
            &control,
            &boundary,
            &source,
            &length,
            &output,
            &scratch_refs,
        );
    } else {
        controlled_xor_saturating_difference_inplace_boundary(
            &mut circ,
            &control,
            &boundary,
            &source,
            &length,
            &output,
            &scratch_refs,
        );
    }
    drop(scratch_refs);

    let data_ids: Vec<u32> = std::iter::once(&control)
        .chain(&boundary)
        .chain(&length)
        .chain(std::iter::once(&scratch[5]))
        .chain(std::iter::once(&scratch[6]))
        .chain(std::iter::once(&scratch[7]))
        .chain(&output)
        .chain(std::iter::once(&dirty))
        .map(QReg::id)
        .collect();
    let preserved_mask = qreg_mask(
        std::iter::once(&control)
            .chain(&boundary)
            .chain(&length)
            .chain(std::iter::once(&scratch[5]))
            .chain(std::iter::once(&scratch[6]))
            .chain(std::iter::once(&scratch[7]))
            .chain(std::iter::once(&dirty)),
    );
    let output_mask = qreg_mask(&output);
    let clean_scratch_mask = qreg_mask(&scratch[..5]);
    SplitThreeHighStageHarness {
        builder: circ.into_builder(),
        data_ids,
        control_id: control.id(),
        boundary_ids: boundary.iter().map(QReg::id).collect(),
        length_ids: length.iter().map(QReg::id).collect(),
        high6_id: scratch[5].id(),
        high7_id: scratch[6].id(),
        high8_id: scratch[7].id(),
        output_ids: output.iter().map(QReg::id).collect(),
        dirty_id: dirty.id(),
        preserved_mask,
        output_mask,
        clean_scratch_mask,
    }
}

fn split_three_high_stage_input(harness: &SplitThreeHighStageHarness, value: u64) -> u64 {
    harness
        .data_ids
        .iter()
        .enumerate()
        .fold(0u64, |state, (bit, id)| {
            state | (((value >> bit) & 1) << id)
        })
}

struct SplitThreeHighIndicatorHarness {
    builder: B,
    source_ids: Vec<u32>,
    high6_id: u32,
    high7_id: u32,
    high8_id: u32,
    clean_scratch_ids: Vec<u32>,
}

fn build_split_three_high_indicator_harness(
    source_is_complemented: bool,
) -> SplitThreeHighIndicatorHarness {
    const SOURCE_WIDTH: usize = 259;

    let mut circ = Circuit::new();
    let source = circ.alloc_input_qreg_bits("sub800.three-high-indicator.source", SOURCE_WIDTH);
    let scratch = circ.alloc_qreg_bits("sub800.three-high-indicator.scratch", 8);
    let source_refs: Vec<&QReg> = source.iter().collect();
    let scratch_refs: Vec<&QReg> = scratch.iter().collect();
    toggle_split_bit_length_three_high(
        &mut circ,
        &source_refs,
        &scratch[5],
        &scratch[6],
        &scratch[7],
        &scratch_refs,
        source_is_complemented,
    );
    drop(source_refs);
    drop(scratch_refs);

    SplitThreeHighIndicatorHarness {
        builder: circ.into_builder(),
        source_ids: source.iter().map(QReg::id).collect(),
        high6_id: scratch[5].id(),
        high7_id: scratch[6].id(),
        high8_id: scratch[7].id(),
        clean_scratch_ids: scratch[..5].iter().map(QReg::id).collect(),
    }
}

fn verify_split_three_high_indicator_harness(
    harness: &SplitThreeHighIndicatorHarness,
    source_is_complemented: bool,
) -> (usize, usize, usize, usize) {
    use crate::circuit::OperationType;

    assert!(harness.builder.ops.iter().all(|operation| matches!(
        operation.kind,
        OperationType::X | OperationType::CX | OperationType::CCX
    )));
    let total_qubits = harness.builder.next_qubit as usize;
    let mut basis_states = 0usize;
    let mut inverse_checks = 0usize;
    let mut source_restore_checks = 0usize;
    let mut scratch_clean_checks = 0usize;
    for bit_length in 0usize..=259 {
        for initial_highs in 0usize..8 {
            let mut input = vec![false; total_qubits];
            for (bit, id) in harness.source_ids.iter().copied().enumerate() {
                let logical = bit_length != 0 && bit == bit_length - 1;
                input[id as usize] = logical ^ source_is_complemented;
            }
            input[harness.high6_id as usize] = (initial_highs & 1) != 0;
            input[harness.high7_id as usize] = (initial_highs & 2) != 0;
            input[harness.high8_id as usize] = (initial_highs & 4) != 0;

            let output = apply_basis_vector(&harness.builder.ops, input.clone());
            let expected_high6 = ((bit_length >> 6) & 1) != 0;
            let expected_high7 = ((bit_length >> 7) & 1) != 0;
            let expected_high8 = ((bit_length >> 8) & 1) != 0;
            assert_eq!(
                output[harness.high6_id as usize],
                input[harness.high6_id as usize] ^ expected_high6
            );
            assert_eq!(
                output[harness.high7_id as usize],
                input[harness.high7_id as usize] ^ expected_high7
            );
            assert_eq!(
                output[harness.high8_id as usize],
                input[harness.high8_id as usize] ^ expected_high8
            );
            assert!(harness
                .source_ids
                .iter()
                .all(|id| output[*id as usize] == input[*id as usize]));
            assert!(harness
                .clean_scratch_ids
                .iter()
                .all(|id| !output[*id as usize]));
            assert_eq!(apply_basis_vector(&harness.builder.ops, output), input);
            basis_states += 1;
            inverse_checks += 1;
            source_restore_checks += 1;
            scratch_clean_checks += 1;
        }
    }
    (
        basis_states,
        inverse_checks,
        source_restore_checks,
        scratch_clean_checks,
    )
}

fn verify_split_three_high_stage(
    harness: &SplitThreeHighStageHarness,
) -> (usize, usize, usize, usize, usize, usize) {
    use crate::circuit::OperationType;

    assert!(harness.builder.ops.iter().all(|operation| matches!(
        operation.kind,
        OperationType::X | OperationType::CX | OperationType::CCX
    )));
    let states = 1u64 << harness.data_ids.len();
    let mut scalar_equivalence_checks = 0usize;
    let mut inverse_pair_checks = 0usize;
    let mut control_off_checks = 0usize;
    let mut phase_clean_checks = 0usize;
    let mut scratch_restore_checks = 0usize;
    let mut dirty_lender_restore_checks = 0usize;
    for value in 0..states {
        let input = split_three_high_stage_input(harness, value);
        let output = apply_scalar(&harness.builder.ops, input);
        let control = ((input >> harness.control_id) & 1) as usize;
        let boundary = read_paired_bitlength_register(input, &harness.boundary_ids) as usize;
        let low = read_paired_bitlength_register(input, &harness.length_ids) as usize;
        let high6 = ((input >> harness.high6_id) & 1) as usize;
        let high7 = ((input >> harness.high7_id) & 1) as usize;
        let high8 = ((input >> harness.high8_id) & 1) as usize;
        let length = low
            | (high6 << harness.length_ids.len())
            | (high7 << (harness.length_ids.len() + 1))
            | (high8 << (harness.length_ids.len() + 2));
        let initial_output = read_paired_bitlength_register(input, &harness.output_ids) as usize;
        let expected_xor = if control == 1 && length >= boundary {
            length - boundary
        } else {
            0
        };
        assert_eq!(
            read_paired_bitlength_register(output, &harness.output_ids) as usize,
            initial_output ^ expected_xor
        );
        assert_eq!(output & harness.preserved_mask, input & harness.preserved_mask);
        assert_eq!(output & harness.clean_scratch_mask, 0);
        assert_eq!(
            (output >> harness.dirty_id) & 1,
            (input >> harness.dirty_id) & 1
        );
        assert_eq!(apply_scalar(&harness.builder.ops, output), input);
        if control == 0 {
            assert_eq!(output & harness.output_mask, input & harness.output_mask);
            control_off_checks += 1;
        }
        scalar_equivalence_checks += 1;
        inverse_pair_checks += 1;
        phase_clean_checks += 1;
        scratch_restore_checks += 1;
        dirty_lender_restore_checks += 1;
    }
    (
        scalar_equivalence_checks,
        inverse_pair_checks,
        control_off_checks,
        phase_clean_checks,
        scratch_restore_checks,
        dirty_lender_restore_checks,
    )
}

/// Prove the six-low/three-high decomposition for the audited 259-bit source.
#[doc(hidden)]
#[must_use]
pub fn exhaustive_split_three_high_rotated_length_check(
) -> SplitThreeHighRotatedLengthProofReport {
    assert!(
        std::env::var_os(SUB800_SPLIT_THREE_HIGH_ROTATED_LENGTH_FLAG).is_none(),
        "the split-three-high rotated-length feature must default off"
    );
    let _environment = RawBitLengthProofEnvironment::capture();

    let mut borrow_truth_table_cases = 0usize;
    for x in 0usize..=1 {
        for y in 0usize..=1 {
            for borrow in 0usize..=1 {
                let expected = usize::from((2 * x) < (2 * y + borrow));
                assert_eq!(full_subtraction_borrow_anf(x, y, borrow), expected);
                borrow_truth_table_cases += 1;
            }
        }
    }

    let mut bit_lengths_checked = 0usize;
    let mut regression_192_255_cases = 0usize;
    let mut same_arithmetic_cases = 0usize;
    let mut mixed_arithmetic_cases = 0usize;
    for bit_length in 0usize..=259 {
        let low = bit_length & 0x3f;
        let high6 = (bit_length >> 6) & 1;
        let high7 = (bit_length >> 7) & 1;
        let high8 = (bit_length >> 8) & 1;
        let z64 = usize::from(bit_length < 64);
        let z128 = usize::from(bit_length < 128);
        let z192 = usize::from(bit_length < 192);
        let z256 = usize::from(bit_length < 256);
        assert_eq!(high6, z64 ^ z128 ^ z192 ^ z256);
        assert_eq!(high7, z128 ^ z256);
        assert_eq!(high8, 1 ^ z256);
        assert_eq!(low | (high6 << 6) | (high7 << 7) | (high8 << 8), bit_length);
        if (192..=255).contains(&bit_length) {
            assert_eq!((high6, high7, high8), (1, 1, 0));
            regression_192_255_cases += 1;
        }
        bit_lengths_checked += 1;

        for boundary in 0usize..512 {
            let boundary_low = boundary & 0x3f;
            let boundary_high6 = (boundary >> 6) & 1;
            let boundary_high7 = (boundary >> 7) & 1;
            let boundary_high8 = (boundary >> 8) & 1;
            let low_borrow = usize::from(low < boundary_low);
            let high6_borrow =
                full_subtraction_borrow_anf(high6, boundary_high6, low_borrow);
            let high7_borrow =
                full_subtraction_borrow_anf(high7, boundary_high7, high6_borrow);
            let underflow =
                full_subtraction_borrow_anf(high8, boundary_high8, high7_borrow);
            let difference = ((low + 64 - boundary_low) & 0x3f)
                | ((high6 ^ boundary_high6 ^ low_borrow) << 6)
                | ((high7 ^ boundary_high7 ^ high6_borrow) << 7)
                | ((high8 ^ boundary_high8 ^ high7_borrow) << 8);
            assert_eq!(underflow, usize::from(bit_length < boundary));
            assert_eq!(difference, bit_length.wrapping_sub(boundary) & 0x1ff);
            same_arithmetic_cases += 1;
        }

        for boundary in 0usize..256 {
            let boundary_low = boundary & 0x3f;
            let boundary_high6 = (boundary >> 6) & 1;
            let boundary_high7 = (boundary >> 7) & 1;
            let low_borrow = usize::from(low < boundary_low);
            let high6_borrow =
                full_subtraction_borrow_anf(high6, boundary_high6, low_borrow);
            let high7_borrow =
                full_subtraction_borrow_anf(high7, boundary_high7, high6_borrow);
            let underflow = full_subtraction_borrow_anf(high8, 0, high7_borrow);
            let difference = ((low + 64 - boundary_low) & 0x3f)
                | ((high6 ^ boundary_high6 ^ low_borrow) << 6)
                | ((high7 ^ boundary_high7 ^ high6_borrow) << 7)
                | ((high8 ^ high7_borrow) << 8);
            assert_eq!(underflow, usize::from(bit_length < boundary));
            assert_eq!(difference, bit_length.wrapping_sub(boundary) & 0x1ff);
            mixed_arithmetic_cases += 1;
        }
    }

    let same_stage = build_split_three_high_stage_harness(false);
    let mixed_stage = build_split_three_high_stage_harness(true);
    let high_indicator = build_split_three_high_indicator_harness(false);
    let complemented_high_indicator = build_split_three_high_indicator_harness(true);
    let same_circuit_basis_states = 1usize << same_stage.data_ids.len();
    let mixed_circuit_basis_states = 1usize << mixed_stage.data_ids.len();
    let same_checks = verify_split_three_high_stage(&same_stage);
    let mixed_checks = verify_split_three_high_stage(&mixed_stage);
    let high_checks = verify_split_three_high_indicator_harness(&high_indicator, false);
    let complemented_high_checks =
        verify_split_three_high_indicator_harness(&complemented_high_indicator, true);

    configure_paired_bitlength_source_complement_proof(true, true);
    std::env::set_var(SUB800_MIXED_BOUNDARY_SCRATCH_EXTENSION_FLAG, "1");
    std::env::set_var(SUB800_BORROWED_ROTATED_UNDERFLOW_FLAG, "1");
    std::env::set_var(SUB800_SPLIT_MIXED_ROTATED_LENGTH_FLAG, "1");
    std::env::set_var(SUB800_SPLIT_SAME_ROTATED_LENGTH_FLAG, "1");
    std::env::set_var(SUB800_SPLIT_TWO_HIGH_ROTATED_LENGTH_FLAG, "1");
    std::env::remove_var(SUB800_SPLIT_THREE_HIGH_ROTATED_LENGTH_FLAG);
    let same_boundary_two_high = build_paired_bitlength_source_complement_harness(
        REFERENCE_LENGTH_WIDTH,
        259,
        8,
        SaturatingDifferenceBoundaryRoute::Inplace,
        false,
        Some(true),
    );
    let mixed_boundary_two_high = build_paired_bitlength_source_complement_harness(
        REFERENCE_LENGTH_WIDTH,
        259,
        8,
        SaturatingDifferenceBoundaryRoute::Inplace,
        true,
        Some(true),
    );
    // Keep two-high enabled: the candidate must take precedence over it.
    std::env::set_var(SUB800_SPLIT_THREE_HIGH_ROTATED_LENGTH_FLAG, "1");
    let same_boundary_candidate = build_paired_bitlength_source_complement_harness(
        REFERENCE_LENGTH_WIDTH,
        259,
        8,
        SaturatingDifferenceBoundaryRoute::Inplace,
        false,
        Some(true),
    );
    let mixed_boundary_candidate = build_paired_bitlength_source_complement_harness(
        REFERENCE_LENGTH_WIDTH,
        259,
        8,
        SaturatingDifferenceBoundaryRoute::Inplace,
        true,
        Some(true),
    );
    let same_boundary_two_high =
        inplace_rotated_boundary_local_resources(&same_boundary_two_high);
    let same_boundary_candidate =
        inplace_rotated_boundary_local_resources(&same_boundary_candidate);
    let mixed_boundary_two_high =
        inplace_rotated_boundary_local_resources(&mixed_boundary_two_high);
    let mixed_boundary_candidate =
        inplace_rotated_boundary_local_resources(&mixed_boundary_candidate);
    assert_eq!(
        same_boundary_candidate.peak_qubits + 1,
        same_boundary_two_high.peak_qubits
    );
    assert_eq!(
        mixed_boundary_candidate.peak_qubits + 1,
        mixed_boundary_two_high.peak_qubits
    );

    SplitThreeHighRotatedLengthProofReport {
        borrow_truth_table_cases,
        bit_lengths_checked,
        regression_192_255_cases,
        complement_modes_checked: 2,
        initial_high_states_checked: 8,
        high_indicator_basis_states: high_checks.0 + complemented_high_checks.0,
        high_indicator_inverse_checks: high_checks.1 + complemented_high_checks.1,
        high_indicator_source_restore_checks: high_checks.2 + complemented_high_checks.2,
        high_indicator_scratch_clean_checks: high_checks.3 + complemented_high_checks.3,
        same_arithmetic_cases,
        mixed_arithmetic_cases,
        same_circuit_basis_states,
        mixed_circuit_basis_states,
        scalar_equivalence_checks: same_checks.0 + mixed_checks.0,
        inverse_pair_checks: same_checks.1 + mixed_checks.1,
        control_off_checks: same_checks.2 + mixed_checks.2,
        phase_clean_checks: same_checks.3 + mixed_checks.3,
        scratch_restore_checks: same_checks.4 + mixed_checks.4,
        dirty_lender_restore_checks: same_checks.5 + mixed_checks.5,
        route_precedence_checks: 2,
        same_boundary_two_high,
        same_boundary_candidate,
        mixed_boundary_two_high,
        mixed_boundary_candidate,
    }
}

struct SplitFourHighStageHarness {
    builder: B,
    data_ids: Vec<u32>,
    control_id: u32,
    boundary_ids: Vec<u32>,
    length_ids: Vec<u32>,
    high5_id: u32,
    high6_id: u32,
    high7_id: u32,
    high8_id: u32,
    output_ids: Vec<u32>,
    dirty_id: u32,
    preserved_mask: u64,
    output_mask: u64,
    clean_scratch_mask: u64,
}

fn build_split_four_high_stage_harness(
    mixed_boundary: bool,
    low_width: usize,
) -> SplitFourHighStageHarness {
    assert!(low_width >= 1);
    let output_width = low_width + 4;

    let mut circ = Circuit::new();
    let control = circ.alloc_qreg("sub800.four-high-proof.control");
    let boundary = circ.alloc_qreg_bits(
        "sub800.four-high-proof.boundary",
        output_width - usize::from(mixed_boundary),
    );
    let length = circ.alloc_qreg_bits("sub800.four-high-proof.length", low_width);
    let output = circ.alloc_qreg_bits("sub800.four-high-proof.output", output_width);
    let scratch = circ.alloc_qreg_bits("sub800.four-high-proof.scratch", 9);
    let dirty = circ.alloc_qreg("sub800.four-high-proof.dirty");
    let source = [&dirty];
    let scratch_refs: Vec<&QReg> = scratch.iter().collect();
    if mixed_boundary {
        controlled_xor_saturating_difference_inplace_mixed_boundary(
            &mut circ,
            &control,
            &boundary,
            &source,
            &length,
            &output,
            &scratch_refs,
        );
    } else {
        controlled_xor_saturating_difference_inplace_boundary(
            &mut circ,
            &control,
            &boundary,
            &source,
            &length,
            &output,
            &scratch_refs,
        );
    }
    drop(scratch_refs);

    let data_ids: Vec<u32> = std::iter::once(&control)
        .chain(&boundary)
        .chain(&length)
        .chain(&scratch[5..9])
        .chain(&output)
        .chain(std::iter::once(&dirty))
        .map(QReg::id)
        .collect();
    let preserved_mask = qreg_mask(
        std::iter::once(&control)
            .chain(&boundary)
            .chain(&length)
            .chain(&scratch[5..9])
            .chain(std::iter::once(&dirty)),
    );
    let output_mask = qreg_mask(&output);
    let clean_scratch_mask = qreg_mask(&scratch[..5]);
    SplitFourHighStageHarness {
        builder: circ.into_builder(),
        data_ids,
        control_id: control.id(),
        boundary_ids: boundary.iter().map(QReg::id).collect(),
        length_ids: length.iter().map(QReg::id).collect(),
        high5_id: scratch[5].id(),
        high6_id: scratch[6].id(),
        high7_id: scratch[7].id(),
        high8_id: scratch[8].id(),
        output_ids: output.iter().map(QReg::id).collect(),
        dirty_id: dirty.id(),
        preserved_mask,
        output_mask,
        clean_scratch_mask,
    }
}

fn split_four_high_stage_input(harness: &SplitFourHighStageHarness, value: u64) -> u64 {
    harness
        .data_ids
        .iter()
        .enumerate()
        .fold(0u64, |state, (bit, id)| {
            state | (((value >> bit) & 1) << id)
        })
}

struct SplitFourHighIndicatorHarness {
    builder: B,
    source_ids: Vec<u32>,
    high5_id: u32,
    high6_id: u32,
    high7_id: u32,
    high8_id: u32,
    clean_scratch_ids: Vec<u32>,
}

fn build_split_four_high_indicator_harness(
    source_is_complemented: bool,
) -> SplitFourHighIndicatorHarness {
    const SOURCE_WIDTH: usize = 259;

    let mut circ = Circuit::new();
    let source = circ.alloc_input_qreg_bits("sub800.four-high-indicator.source", SOURCE_WIDTH);
    let scratch = circ.alloc_qreg_bits("sub800.four-high-indicator.scratch", 9);
    let source_refs: Vec<&QReg> = source.iter().collect();
    let scratch_refs: Vec<&QReg> = scratch.iter().collect();
    toggle_split_bit_length_four_high(
        &mut circ,
        &source_refs,
        &scratch[5],
        &scratch[6],
        &scratch[7],
        &scratch[8],
        &scratch_refs,
        source_is_complemented,
    );
    drop(source_refs);
    drop(scratch_refs);

    SplitFourHighIndicatorHarness {
        builder: circ.into_builder(),
        source_ids: source.iter().map(QReg::id).collect(),
        high5_id: scratch[5].id(),
        high6_id: scratch[6].id(),
        high7_id: scratch[7].id(),
        high8_id: scratch[8].id(),
        clean_scratch_ids: scratch[..5].iter().map(QReg::id).collect(),
    }
}

fn verify_split_four_high_indicator_harness(
    harness: &SplitFourHighIndicatorHarness,
    source_is_complemented: bool,
) -> (usize, usize, usize, usize) {
    use crate::circuit::OperationType;

    assert!(harness.builder.ops.iter().all(|operation| matches!(
        operation.kind,
        OperationType::X | OperationType::CX | OperationType::CCX
    )));
    let total_qubits = harness.builder.next_qubit as usize;
    let mut basis_states = 0usize;
    let mut inverse_checks = 0usize;
    let mut source_restore_checks = 0usize;
    let mut scratch_clean_checks = 0usize;
    for bit_length in 0usize..=259 {
        for source_class in 0usize..3 {
            for initial_highs in 0usize..16 {
                let mut input = vec![false; total_qubits];
                for (bit, id) in harness.source_ids.iter().copied().enumerate() {
                    let logical = if bit_length == 0 {
                        false
                    } else {
                        let top = bit_length - 1;
                        match source_class {
                            0 => bit == top,
                            1 => bit <= top,
                            2 => bit == top || (bit < top && bit % 2 == 0),
                            _ => unreachable!(),
                        }
                    };
                    input[id as usize] = logical ^ source_is_complemented;
                }
                input[harness.high5_id as usize] = (initial_highs & 1) != 0;
                input[harness.high6_id as usize] = (initial_highs & 2) != 0;
                input[harness.high7_id as usize] = (initial_highs & 4) != 0;
                input[harness.high8_id as usize] = (initial_highs & 8) != 0;

                let output = apply_basis_vector(&harness.builder.ops, input.clone());
                for (id, expected) in [
                    (harness.high5_id, ((bit_length >> 5) & 1) != 0),
                    (harness.high6_id, ((bit_length >> 6) & 1) != 0),
                    (harness.high7_id, ((bit_length >> 7) & 1) != 0),
                    (harness.high8_id, ((bit_length >> 8) & 1) != 0),
                ] {
                    assert_eq!(output[id as usize], input[id as usize] ^ expected);
                }
                assert!(harness
                    .source_ids
                    .iter()
                    .all(|id| output[*id as usize] == input[*id as usize]));
                assert!(harness
                    .clean_scratch_ids
                    .iter()
                    .all(|id| !output[*id as usize]));
                assert_eq!(apply_basis_vector(&harness.builder.ops, output), input);
                basis_states += 1;
                inverse_checks += 1;
                source_restore_checks += 1;
                scratch_clean_checks += 1;
            }
        }
    }
    (
        basis_states,
        inverse_checks,
        source_restore_checks,
        scratch_clean_checks,
    )
}

fn verify_split_four_high_stage(
    harness: &SplitFourHighStageHarness,
) -> (usize, usize, usize, usize, usize, usize) {
    use crate::circuit::OperationType;

    assert!(harness.builder.ops.iter().all(|operation| matches!(
        operation.kind,
        OperationType::X | OperationType::CX | OperationType::CCX
    )));
    let states = 1u64 << harness.data_ids.len();
    let mut scalar_equivalence_checks = 0usize;
    let mut inverse_pair_checks = 0usize;
    let mut control_off_checks = 0usize;
    let mut phase_clean_checks = 0usize;
    let mut scratch_restore_checks = 0usize;
    let mut dirty_lender_restore_checks = 0usize;
    for value in 0..states {
        let input = split_four_high_stage_input(harness, value);
        let output = apply_scalar(&harness.builder.ops, input);
        let control = ((input >> harness.control_id) & 1) as usize;
        let boundary = read_paired_bitlength_register(input, &harness.boundary_ids) as usize;
        let low = read_paired_bitlength_register(input, &harness.length_ids) as usize;
        let high5 = ((input >> harness.high5_id) & 1) as usize;
        let high6 = ((input >> harness.high6_id) & 1) as usize;
        let high7 = ((input >> harness.high7_id) & 1) as usize;
        let high8 = ((input >> harness.high8_id) & 1) as usize;
        let length = low
            | (high5 << harness.length_ids.len())
            | (high6 << (harness.length_ids.len() + 1))
            | (high7 << (harness.length_ids.len() + 2))
            | (high8 << (harness.length_ids.len() + 3));
        let initial_output = read_paired_bitlength_register(input, &harness.output_ids) as usize;
        let expected_xor = if control == 1 && length >= boundary {
            length - boundary
        } else {
            0
        };
        assert_eq!(
            read_paired_bitlength_register(output, &harness.output_ids) as usize,
            initial_output ^ expected_xor
        );
        assert_eq!(output & harness.preserved_mask, input & harness.preserved_mask);
        assert_eq!(output & harness.clean_scratch_mask, 0);
        assert_eq!(
            (output >> harness.dirty_id) & 1,
            (input >> harness.dirty_id) & 1
        );
        assert_eq!(apply_scalar(&harness.builder.ops, output), input);
        if control == 0 {
            assert_eq!(output & harness.output_mask, input & harness.output_mask);
            control_off_checks += 1;
        }
        scalar_equivalence_checks += 1;
        inverse_pair_checks += 1;
        phase_clean_checks += 1;
        scratch_restore_checks += 1;
        dirty_lender_restore_checks += 1;
    }
    (
        scalar_equivalence_checks,
        inverse_pair_checks,
        control_off_checks,
        phase_clean_checks,
        scratch_restore_checks,
        dirty_lender_restore_checks,
    )
}

fn split_four_high_structured_input(
    harness: &SplitFourHighStageHarness,
    bit_length: usize,
    boundary: usize,
    initial_output: usize,
    dirty: bool,
) -> u64 {
    let mut input = 1u64 << harness.control_id;
    for (bit, id) in harness.boundary_ids.iter().copied().enumerate() {
        input |= (((boundary >> bit) & 1) as u64) << id;
    }
    for (bit, id) in harness.length_ids.iter().copied().enumerate() {
        input |= (((bit_length >> bit) & 1) as u64) << id;
    }
    for (bit, id) in [
        harness.high5_id,
        harness.high6_id,
        harness.high7_id,
        harness.high8_id,
    ]
    .into_iter()
    .enumerate()
    {
        input |= (((bit_length >> (harness.length_ids.len() + bit)) & 1) as u64) << id;
    }
    for (bit, id) in harness.output_ids.iter().copied().enumerate() {
        input |= (((initial_output >> bit) & 1) as u64) << id;
    }
    input | (u64::from(dirty) << harness.dirty_id)
}

fn verify_split_four_high_production_stage(
    harness: &SplitFourHighStageHarness,
    mixed_boundary: bool,
) -> (usize, usize, usize, usize) {
    use crate::circuit::OperationType;

    assert_eq!(harness.length_ids.len(), 5);
    assert_eq!(harness.output_ids.len(), 9);
    assert_eq!(harness.boundary_ids.len(), 9 - usize::from(mixed_boundary));
    assert!(harness.builder.ops.iter().all(|operation| matches!(
        operation.kind,
        OperationType::X | OperationType::CX | OperationType::CCX
    )));
    let boundary_limit = 1usize << harness.boundary_ids.len();
    let output_mask = (1usize << harness.output_ids.len()) - 1;
    let mut cases_checked = 0usize;
    let mut inverse_checks = 0usize;
    let mut scratch_restore_checks = 0usize;
    let mut dirty_lender_restore_checks = 0usize;
    for bit_length in 0usize..=259 {
        for boundary in 0usize..boundary_limit {
            for output_class in 0usize..2 {
                let initial_output = if output_class == 0 {
                    0
                } else {
                    (bit_length.wrapping_mul(257)
                        ^ boundary.wrapping_mul(17)
                        ^ 0x155)
                        & output_mask
                };
                let input = split_four_high_structured_input(
                    harness,
                    bit_length,
                    boundary,
                    initial_output,
                    output_class != 0,
                );
                let output = apply_scalar(&harness.builder.ops, input);
                let expected_xor = if bit_length >= boundary {
                    bit_length - boundary
                } else {
                    0
                };
                assert_eq!(
                    read_paired_bitlength_register(output, &harness.output_ids) as usize,
                    initial_output ^ expected_xor
                );
                assert_eq!(output & harness.preserved_mask, input & harness.preserved_mask);
                assert_eq!(output & harness.clean_scratch_mask, 0);
                assert_eq!(
                    (output >> harness.dirty_id) & 1,
                    (input >> harness.dirty_id) & 1
                );
                assert_eq!(apply_scalar(&harness.builder.ops, output), input);
                cases_checked += 1;
                inverse_checks += 1;
                scratch_restore_checks += 1;
                dirty_lender_restore_checks += 1;
            }
        }
    }
    (
        cases_checked,
        inverse_checks,
        scratch_restore_checks,
        dirty_lender_restore_checks,
    )
}

/// Prove the five-low/four-high decomposition for the audited 259-bit source.
#[doc(hidden)]
#[must_use]
pub fn exhaustive_split_four_high_rotated_length_check(
) -> SplitFourHighRotatedLengthProofReport {
    assert!(
        std::env::var_os(SUB800_SPLIT_FOUR_HIGH_ROTATED_LENGTH_FLAG).is_none(),
        "the split-four-high proof process must not inherit a route override"
    );
    let _environment = RawBitLengthProofEnvironment::capture();

    let mut borrow_truth_table_cases = 0usize;
    for x in 0usize..=1 {
        for y in 0usize..=1 {
            for borrow in 0usize..=1 {
                let expected = usize::from((2 * x) < (2 * y + borrow));
                assert_eq!(full_subtraction_borrow_anf(x, y, borrow), expected);
                borrow_truth_table_cases += 1;
            }
        }
    }

    let mut bit_lengths_checked = 0usize;
    let mut regression_224_255_cases = 0usize;
    let mut same_arithmetic_cases = 0usize;
    let mut mixed_arithmetic_cases = 0usize;
    for bit_length in 0usize..=259 {
        let low = bit_length & 0x1f;
        let high5 = (bit_length >> 5) & 1;
        let high6 = (bit_length >> 6) & 1;
        let high7 = (bit_length >> 7) & 1;
        let high8 = (bit_length >> 8) & 1;
        let z32 = usize::from(bit_length < 32);
        let z64 = usize::from(bit_length < 64);
        let z96 = usize::from(bit_length < 96);
        let z128 = usize::from(bit_length < 128);
        let z160 = usize::from(bit_length < 160);
        let z192 = usize::from(bit_length < 192);
        let z224 = usize::from(bit_length < 224);
        let z256 = usize::from(bit_length < 256);
        assert_eq!(high5, z32 ^ z64 ^ z96 ^ z128 ^ z160 ^ z192 ^ z224 ^ z256);
        assert_eq!(high6, z64 ^ z128 ^ z192 ^ z256);
        assert_eq!(high7, z128 ^ z256);
        assert_eq!(high8, 1 ^ z256);
        assert_eq!(
            low | (high5 << 5) | (high6 << 6) | (high7 << 7) | (high8 << 8),
            bit_length
        );
        if (224..=255).contains(&bit_length) {
            assert_eq!((high5, high6, high7, high8), (1, 1, 1, 0));
            regression_224_255_cases += 1;
        }
        bit_lengths_checked += 1;

        for boundary in 0usize..512 {
            let boundary_low = boundary & 0x1f;
            let boundary_high5 = (boundary >> 5) & 1;
            let boundary_high6 = (boundary >> 6) & 1;
            let boundary_high7 = (boundary >> 7) & 1;
            let boundary_high8 = (boundary >> 8) & 1;
            let low_borrow = usize::from(low < boundary_low);
            let high5_borrow =
                full_subtraction_borrow_anf(high5, boundary_high5, low_borrow);
            let high6_borrow =
                full_subtraction_borrow_anf(high6, boundary_high6, high5_borrow);
            let high7_borrow =
                full_subtraction_borrow_anf(high7, boundary_high7, high6_borrow);
            let underflow =
                full_subtraction_borrow_anf(high8, boundary_high8, high7_borrow);
            let difference = ((low + 32 - boundary_low) & 0x1f)
                | ((high5 ^ boundary_high5 ^ low_borrow) << 5)
                | ((high6 ^ boundary_high6 ^ high5_borrow) << 6)
                | ((high7 ^ boundary_high7 ^ high6_borrow) << 7)
                | ((high8 ^ boundary_high8 ^ high7_borrow) << 8);
            assert_eq!(underflow, usize::from(bit_length < boundary));
            assert_eq!(difference, bit_length.wrapping_sub(boundary) & 0x1ff);
            same_arithmetic_cases += 1;
        }

        for boundary in 0usize..256 {
            let boundary_low = boundary & 0x1f;
            let boundary_high5 = (boundary >> 5) & 1;
            let boundary_high6 = (boundary >> 6) & 1;
            let boundary_high7 = (boundary >> 7) & 1;
            let low_borrow = usize::from(low < boundary_low);
            let high5_borrow =
                full_subtraction_borrow_anf(high5, boundary_high5, low_borrow);
            let high6_borrow =
                full_subtraction_borrow_anf(high6, boundary_high6, high5_borrow);
            let high7_borrow =
                full_subtraction_borrow_anf(high7, boundary_high7, high6_borrow);
            let underflow = full_subtraction_borrow_anf(high8, 0, high7_borrow);
            let difference = ((low + 32 - boundary_low) & 0x1f)
                | ((high5 ^ boundary_high5 ^ low_borrow) << 5)
                | ((high6 ^ boundary_high6 ^ high5_borrow) << 6)
                | ((high7 ^ boundary_high7 ^ high6_borrow) << 7)
                | ((high8 ^ high7_borrow) << 8);
            assert_eq!(underflow, usize::from(bit_length < boundary));
            assert_eq!(difference, bit_length.wrapping_sub(boundary) & 0x1ff);
            mixed_arithmetic_cases += 1;
        }
    }

    let same_stage = build_split_four_high_stage_harness(false, 1);
    let mixed_stage = build_split_four_high_stage_harness(true, 1);
    let production_same_stage = build_split_four_high_stage_harness(false, 5);
    let production_mixed_stage = build_split_four_high_stage_harness(true, 5);
    let high_indicator = build_split_four_high_indicator_harness(false);
    let complemented_high_indicator = build_split_four_high_indicator_harness(true);
    let same_circuit_basis_states = 1usize << same_stage.data_ids.len();
    let mixed_circuit_basis_states = 1usize << mixed_stage.data_ids.len();
    let same_checks = verify_split_four_high_stage(&same_stage);
    let mixed_checks = verify_split_four_high_stage(&mixed_stage);
    let production_same_checks =
        verify_split_four_high_production_stage(&production_same_stage, false);
    let production_mixed_checks =
        verify_split_four_high_production_stage(&production_mixed_stage, true);
    let high_checks = verify_split_four_high_indicator_harness(&high_indicator, false);
    let complemented_high_checks =
        verify_split_four_high_indicator_harness(&complemented_high_indicator, true);

    reset_sub800_q839_route_coverage();
    configure_paired_bitlength_source_complement_proof(true, true);
    std::env::set_var(SUB800_MIXED_BOUNDARY_SCRATCH_EXTENSION_FLAG, "1");
    std::env::set_var(SUB800_BORROWED_ROTATED_UNDERFLOW_FLAG, "1");
    std::env::set_var(SUB800_SPLIT_MIXED_ROTATED_LENGTH_FLAG, "1");
    std::env::set_var(SUB800_SPLIT_SAME_ROTATED_LENGTH_FLAG, "1");
    std::env::set_var(SUB800_SPLIT_TWO_HIGH_ROTATED_LENGTH_FLAG, "1");
    std::env::set_var(SUB800_SPLIT_THREE_HIGH_ROTATED_LENGTH_FLAG, "1");
    std::env::remove_var(SUB800_SPLIT_FOUR_HIGH_ROTATED_LENGTH_FLAG);
    let same_boundary_three_high = build_paired_bitlength_source_complement_harness(
        REFERENCE_LENGTH_WIDTH,
        259,
        9,
        SaturatingDifferenceBoundaryRoute::Inplace,
        false,
        Some(true),
    );
    let mixed_boundary_three_high = build_paired_bitlength_source_complement_harness(
        REFERENCE_LENGTH_WIDTH,
        259,
        9,
        SaturatingDifferenceBoundaryRoute::Inplace,
        true,
        Some(true),
    );
    // Keep all earlier split routes enabled: four-high must take precedence.
    std::env::set_var(SUB800_SPLIT_FOUR_HIGH_ROTATED_LENGTH_FLAG, "1");
    let same_boundary_candidate = build_paired_bitlength_source_complement_harness(
        REFERENCE_LENGTH_WIDTH,
        259,
        9,
        SaturatingDifferenceBoundaryRoute::Inplace,
        false,
        Some(true),
    );
    let mixed_boundary_candidate = build_paired_bitlength_source_complement_harness(
        REFERENCE_LENGTH_WIDTH,
        259,
        9,
        SaturatingDifferenceBoundaryRoute::Inplace,
        true,
        Some(true),
    );
    let same_boundary_three_high =
        inplace_rotated_boundary_local_resources(&same_boundary_three_high);
    let same_boundary_candidate =
        inplace_rotated_boundary_local_resources(&same_boundary_candidate);
    let mixed_boundary_three_high =
        inplace_rotated_boundary_local_resources(&mixed_boundary_three_high);
    let mixed_boundary_candidate =
        inplace_rotated_boundary_local_resources(&mixed_boundary_candidate);
    let route_coverage = sub800_q839_route_coverage();
    assert_eq!(route_coverage.split_three_same_calls, 1);
    assert_eq!(route_coverage.split_three_mixed_calls, 1);
    assert_eq!(route_coverage.split_four_same_calls, 1);
    assert_eq!(route_coverage.split_four_mixed_calls, 1);
    assert_eq!(
        same_boundary_candidate.peak_qubits + 1,
        same_boundary_three_high.peak_qubits
    );
    assert_eq!(
        mixed_boundary_candidate.peak_qubits + 1,
        mixed_boundary_three_high.peak_qubits
    );

    SplitFourHighRotatedLengthProofReport {
        borrow_truth_table_cases,
        bit_lengths_checked,
        regression_224_255_cases,
        complement_modes_checked: 2,
        source_classes_checked: 3,
        initial_high_states_checked: 16,
        high_indicator_basis_states: high_checks.0 + complemented_high_checks.0,
        high_indicator_inverse_checks: high_checks.1 + complemented_high_checks.1,
        high_indicator_source_restore_checks: high_checks.2 + complemented_high_checks.2,
        high_indicator_scratch_clean_checks: high_checks.3 + complemented_high_checks.3,
        same_arithmetic_cases,
        mixed_arithmetic_cases,
        same_circuit_basis_states,
        mixed_circuit_basis_states,
        production_same_circuit_cases: production_same_checks.0,
        production_mixed_circuit_cases: production_mixed_checks.0,
        production_inverse_checks: production_same_checks.1 + production_mixed_checks.1,
        production_scratch_restore_checks: production_same_checks.2
            + production_mixed_checks.2,
        production_dirty_lender_restore_checks: production_same_checks.3
            + production_mixed_checks.3,
        scalar_equivalence_checks: same_checks.0 + mixed_checks.0,
        inverse_pair_checks: same_checks.1 + mixed_checks.1,
        control_off_checks: same_checks.2 + mixed_checks.2,
        phase_clean_checks: same_checks.3 + mixed_checks.3,
        scratch_restore_checks: same_checks.4 + mixed_checks.4,
        dirty_lender_restore_checks: same_checks.5 + mixed_checks.5,
        route_precedence_checks: route_coverage.split_four_same_calls
            + route_coverage.split_four_mixed_calls,
        same_boundary_three_high,
        same_boundary_candidate,
        mixed_boundary_three_high,
        mixed_boundary_candidate,
    }
}

/// Prove the exact cancellation
///
/// `X_S K_add X_S V X_S K_sub X_S = X_S K_add V K_sub X_S`
///
/// for the fused direct-prefix bit-length route. The reduced-width sweep
/// covers both boundary representations, owned and seven-lane borrowed
/// scratch, both arithmetic routes, every external basis state, phase
/// cleanliness, inverse pairing, and complete scratch restoration.
#[doc(hidden)]
#[must_use]
pub fn exhaustive_paired_bitlength_source_complement_check(
) -> PairedBitLengthSourceComplementProofReport {
    use crate::circuit::OperationType;

    const OUTPUT_WIDTH: usize = 3;
    const MAX_SOURCE_WIDTH: usize = 5;
    const LOCAL_SOURCE_WIDTH: usize = 259;
    const LOCAL_OUTPUT_WIDTH: usize = REFERENCE_LENGTH_WIDTH;
    const LOCAL_BOUNDARY_WIDTH: usize = REFERENCE_R_LENGTH_WIDTH;

    assert_ne!(
        std::env::var(PAIRED_BITLEN_SOURCE_COMPLEMENT_FLAG)
            .ok()
            .as_deref(),
        Some("1"),
        "paired source complement must default off"
    );
    let _environment = RawBitLengthProofEnvironment::capture();
    let source_widths = 2usize..=MAX_SOURCE_WIDTH;
    let boundary_forms = [false, true];
    let scratch_modes = [(0usize, false), (7usize, true)];
    let routes = [
        SaturatingDifferenceBoundaryRoute::Materialized,
        SaturatingDifferenceBoundaryRoute::Inplace,
    ];

    let mut configurations_checked = 0usize;
    let mut basis_states_checked = 0usize;
    let mut oracle_checks = 0usize;
    let mut scalar_equivalence_checks = 0usize;
    let mut simulator_equivalence_checks = 0usize;
    let mut phase_clean_checks = 0usize;
    let mut inverse_pair_checks = 0usize;
    let mut control_off_checks = 0usize;
    let mut source_restore_checks = 0usize;
    let mut boundary_restore_checks = 0usize;
    let mut borrowed_scratch_clean_checks = 0usize;
    let mut ancilla_clean_checks = 0usize;
    let mut default_stream_identity_checks = 0usize;
    let mut non_x_kind_identity_checks = 0usize;

    for source_width in source_widths.clone() {
        for mixed_boundary in boundary_forms {
            let boundary_width = OUTPUT_WIDTH - usize::from(mixed_boundary);
            for (scratch_lanes, prefix_scratch_loan) in scratch_modes {
                configure_paired_bitlength_source_complement_proof(
                    prefix_scratch_loan,
                    false,
                );
                for route in routes {
                    let default = build_paired_bitlength_source_complement_harness(
                        OUTPUT_WIDTH,
                        source_width,
                        scratch_lanes,
                        route,
                        mixed_boundary,
                        None,
                    );
                    let baseline = build_paired_bitlength_source_complement_harness(
                        OUTPUT_WIDTH,
                        source_width,
                        scratch_lanes,
                        route,
                        mixed_boundary,
                        Some(false),
                    );
                    let optimized = build_paired_bitlength_source_complement_harness(
                        OUTPUT_WIDTH,
                        source_width,
                        scratch_lanes,
                        route,
                        mixed_boundary,
                        Some(true),
                    );
                    assert_inplace_rotated_boundary_stream_identity(&default, &baseline);
                    default_stream_identity_checks += 1;
                    assert_eq!(baseline.data_ids, optimized.data_ids);
                    assert_eq!(baseline.builder.active_qubits, optimized.builder.active_qubits);
                    assert_eq!(baseline.builder.peak_qubits, optimized.builder.peak_qubits);
                    assert_eq!(
                        baseline.builder.ops.len() - optimized.builder.ops.len(),
                        2 * source_width
                    );
                    assert_eq!(
                        baseline.builder.counted_kind_ops[OperationType::X as usize]
                            - optimized.builder.counted_kind_ops[OperationType::X as usize],
                        2 * source_width
                    );
                    for kind in 0..baseline.builder.counted_kind_ops.len() {
                        if kind != OperationType::X as usize {
                            assert_eq!(
                                baseline.builder.counted_kind_ops[kind],
                                optimized.builder.counted_kind_ops[kind]
                            );
                            non_x_kind_identity_checks += 1;
                        }
                    }

                    let (simulator_cases, simulator_phases) =
                        verify_inplace_rotated_boundary_simulator_equivalence(
                            &baseline,
                            &optimized,
                        );
                    simulator_equivalence_checks += simulator_cases;
                    phase_clean_checks += simulator_phases;

                    let data_states = 1u64 << baseline.data_ids.len();
                    let boundary_mask = (1u64 << boundary_width) - 1;
                    let source_mask = (1u64 << source_width) - 1;
                    let output_mask = (1u64 << OUTPUT_WIDTH) - 1;
                    let boundary_ids = &baseline.data_ids[1..1 + boundary_width];
                    let source_start = 1 + boundary_width;
                    let source_ids = &baseline.data_ids[source_start..source_start + source_width];
                    let output_ids = &baseline.data_ids[source_start + source_width..];
                    for value in 0..data_states {
                        let input = inplace_rotated_boundary_input(&baseline, value);
                        let baseline_output = apply_scalar(&baseline.builder.ops, input);
                        let optimized_output = apply_scalar(&optimized.builder.ops, input);
                        assert_eq!(optimized_output, baseline_output);
                        assert_inplace_rotated_boundary_clean(
                            &baseline,
                            baseline_output,
                            "paired-complement baseline",
                        );
                        assert_inplace_rotated_boundary_clean(
                            &optimized,
                            optimized_output,
                            "paired-complement optimized",
                        );

                        let control = value & 1;
                        let boundary = (value >> 1) & boundary_mask;
                        let source = (value >> (1 + boundary_width)) & source_mask;
                        let output =
                            (value >> (1 + boundary_width + source_width)) & output_mask;
                        let difference = (bit_length_usize(source as usize) as u64)
                            .saturating_sub(boundary);
                        let expected_output = output ^ if control == 1 { difference } else { 0 };
                        for state in [baseline_output, optimized_output] {
                            assert_eq!(read_paired_bitlength_register(state, source_ids), source);
                            assert_eq!(
                                read_paired_bitlength_register(state, boundary_ids),
                                boundary
                            );
                            assert_eq!(
                                read_paired_bitlength_register(state, output_ids),
                                expected_output
                            );
                        }
                        oracle_checks += 2;
                        source_restore_checks += 2;
                        boundary_restore_checks += 2;
                        if control == 0 {
                            assert_eq!(expected_output, output);
                            control_off_checks += 1;
                        }
                        assert_eq!(
                            apply_scalar(&baseline.builder.ops, baseline_output),
                            input
                        );
                        assert_eq!(
                            apply_scalar(&optimized.builder.ops, optimized_output),
                            input
                        );
                        inverse_pair_checks += 2;
                        if scratch_lanes != 0 {
                            borrowed_scratch_clean_checks += 2;
                        }
                        ancilla_clean_checks += 2;
                        scalar_equivalence_checks += 1;
                        basis_states_checked += 1;
                    }
                    configurations_checked += 1;
                }
            }
        }
    }

    configure_paired_bitlength_source_complement_proof(true, true);
    let local_baseline_harness = build_paired_bitlength_source_complement_harness(
        LOCAL_OUTPUT_WIDTH,
        LOCAL_SOURCE_WIDTH,
        7,
        SaturatingDifferenceBoundaryRoute::Configured,
        true,
        Some(false),
    );
    let local_optimized_harness = build_paired_bitlength_source_complement_harness(
        LOCAL_OUTPUT_WIDTH,
        LOCAL_SOURCE_WIDTH,
        7,
        SaturatingDifferenceBoundaryRoute::Configured,
        true,
        Some(true),
    );
    let local_baseline =
        paired_bitlength_source_complement_local_resources(&local_baseline_harness);
    let local_optimized =
        paired_bitlength_source_complement_local_resources(&local_optimized_harness);
    assert_eq!(local_baseline.active_qubits, local_optimized.active_qubits);
    assert_eq!(local_baseline.peak_qubits, local_optimized.peak_qubits);
    assert_eq!(
        local_baseline.emitted_ops - local_optimized.emitted_ops,
        2 * LOCAL_SOURCE_WIDTH
    );
    assert_eq!(
        local_baseline.emitted_x - local_optimized.emitted_x,
        2 * LOCAL_SOURCE_WIDTH
    );
    assert_eq!(
        local_baseline.emitted_toffoli,
        local_optimized.emitted_toffoli
    );

    PairedBitLengthSourceComplementProofReport {
        source_widths_checked: source_widths.count(),
        maximum_source_width: MAX_SOURCE_WIDTH,
        boundary_forms_checked: boundary_forms.len(),
        scratch_modes_checked: scratch_modes.len(),
        boundary_routes_checked: routes.len(),
        configurations_checked,
        basis_states_checked,
        oracle_checks,
        scalar_equivalence_checks,
        simulator_equivalence_checks,
        phase_clean_checks,
        inverse_pair_checks,
        control_off_checks,
        source_restore_checks,
        boundary_restore_checks,
        borrowed_scratch_clean_checks,
        ancilla_clean_checks,
        default_stream_identity_checks,
        non_x_kind_identity_checks,
        local_source_width: LOCAL_SOURCE_WIDTH,
        local_output_width: LOCAL_OUTPUT_WIDTH,
        local_boundary_width: LOCAL_BOUNDARY_WIDTH,
        local_baseline,
        local_optimized,
        local_qubit_delta: local_optimized.peak_qubits as i64
            - local_baseline.peak_qubits as i64,
        local_ops_delta: local_optimized.emitted_ops as i64 - local_baseline.emitted_ops as i64,
        local_x_delta: local_optimized.emitted_x as i64 - local_baseline.emitted_x as i64,
        local_toffoli_delta: local_optimized.emitted_toffoli as i64
            - local_baseline.emitted_toffoli as i64,
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum CoefficientBracketBody {
    AddForward,
    AddReverse,
    SubForward,
    SubReverse,
}

impl CoefficientBracketBody {
    fn inverse(self) -> Self {
        match self {
            Self::AddForward => Self::SubReverse,
            Self::AddReverse => Self::SubForward,
            Self::SubForward => Self::AddReverse,
            Self::SubReverse => Self::AddForward,
        }
    }
}

const COEFFICIENT_BRACKET_BODIES: [CoefficientBracketBody; 4] = [
    CoefficientBracketBody::AddForward,
    CoefficientBracketBody::AddReverse,
    CoefficientBracketBody::SubForward,
    CoefficientBracketBody::SubReverse,
];

#[allow(clippy::too_many_arguments)]
fn emit_coefficient_bracket_body(
    circ: &mut Circuit,
    body: CoefficientBracketBody,
    active: &QReg,
    carry: &QReg,
    work1: &QReg,
    work2: &QReg,
    tmp: &QReg,
) {
    match body {
        CoefficientBracketBody::AddForward => {
            circ.ccx(active, carry, work2);
            circ.ccx(active, carry, work1);
            multi_controlled_x_vchain(
                circ,
                &[active, work1, work2],
                carry,
                std::slice::from_ref(tmp),
            );
        }
        CoefficientBracketBody::AddReverse => {
            multi_controlled_x_vchain(
                circ,
                &[active, work1, work2],
                carry,
                std::slice::from_ref(tmp),
            );
            circ.ccx(active, carry, work1);
            circ.ccx(active, work1, work2);
        }
        CoefficientBracketBody::SubForward => {
            circ.ccx(active, work1, work2);
            circ.ccx(active, carry, work1);
            multi_controlled_x_vchain(
                circ,
                &[active, work1, work2],
                carry,
                std::slice::from_ref(tmp),
            );
        }
        CoefficientBracketBody::SubReverse => {
            multi_controlled_x_vchain(
                circ,
                &[active, work1, work2],
                carry,
                std::slice::from_ref(tmp),
            );
            circ.ccx(active, carry, work1);
            circ.ccx(active, carry, work2);
        }
    }
}

struct CoefficientBracketHarness {
    builder: B,
    data_ids: Vec<u32>,
    cursor_mask: u64,
    active_mask: u64,
    tmp_id: u32,
}

fn set_coefficient_nonnegative_x_cancel_proof_mode(enabled: Option<bool>) {
    match enabled {
        Some(true) => std::env::set_var(COEFFICIENT_NONNEGATIVE_X_CANCEL_FLAG, "1"),
        Some(false) => std::env::set_var(COEFFICIENT_NONNEGATIVE_X_CANCEL_FLAG, "0"),
        None => std::env::remove_var(COEFFICIENT_NONNEGATIVE_X_CANCEL_FLAG),
    }
}

fn build_coefficient_bracket_harness(
    cursor_width: usize,
    body: CoefficientBracketBody,
    enabled: Option<bool>,
) -> CoefficientBracketHarness {
    assert!(cursor_width > 0);
    set_coefficient_nonnegative_x_cancel_proof_mode(enabled);
    let mut circ = Circuit::new();
    let control = circ.alloc_qreg("rs.coeff-xcancel-proof.control");
    let cursor = circ.alloc_qreg_bits("rs.coeff-xcancel-proof.cursor", cursor_width);
    let active = circ.alloc_qreg("rs.coeff-xcancel-proof.active");
    let carry = circ.alloc_qreg("rs.coeff-xcancel-proof.carry");
    let work1 = circ.alloc_qreg("rs.coeff-xcancel-proof.work1");
    let work2 = circ.alloc_qreg("rs.coeff-xcancel-proof.work2");
    let tmp = circ.alloc_qreg("rs.coeff-xcancel-proof.tmp");
    let bracket = production_coefficient_nonnegative_bracket();
    begin_coefficient_nonnegative_bracket(&mut circ, &control, &cursor, &active, bracket);
    emit_coefficient_bracket_body(
        &mut circ,
        body,
        &active,
        &carry,
        &work1,
        &work2,
        &tmp,
    );
    end_coefficient_nonnegative_bracket(&mut circ, &control, &cursor, &active, bracket);

    let data_ids: Vec<u32> = std::iter::once(&control)
        .chain(&cursor)
        .chain([&active, &carry, &work1, &work2])
        .map(QReg::id)
        .collect();
    CoefficientBracketHarness {
        builder: circ.into_builder(),
        data_ids,
        cursor_mask: qreg_mask(&cursor),
        active_mask: 1u64 << active.id(),
        tmp_id: tmp.id(),
    }
}

fn assert_coefficient_bracket_stream_identity(
    left: &CoefficientBracketHarness,
    right: &CoefficientBracketHarness,
) {
    assert_eq!(left.builder.ops, right.builder.ops);
    assert_eq!(left.builder.next_qubit, right.builder.next_qubit);
    assert_eq!(left.builder.next_bit, right.builder.next_bit);
    assert_eq!(left.builder.active_qubits, right.builder.active_qubits);
    assert_eq!(left.builder.peak_qubits, right.builder.peak_qubits);
    assert_eq!(left.builder.free_qubits, right.builder.free_qubits);
    assert_eq!(left.builder.allocation_serial, right.builder.allocation_serial);
}

fn coefficient_bracket_input(harness: &CoefficientBracketHarness, value: u64) -> u64 {
    harness
        .data_ids
        .iter()
        .enumerate()
        .fold(0u64, |state, (bit, id)| {
            state | (((value >> bit) & 1) << id)
        })
}

fn verify_coefficient_bracket_simulator_equivalence(
    baseline: &CoefficientBracketHarness,
    optimized: &CoefficientBracketHarness,
) -> (usize, usize, usize) {
    use crate::circuit::QubitId;
    use crate::sim::Simulator;
    use sha3::{
        digest::{ExtendableOutput, Update},
        Shake128,
    };

    assert_eq!(baseline.data_ids, optimized.data_ids);
    let states = 1usize << baseline.data_ids.len();
    let mut cases_checked = 0usize;
    let mut phase_clean_checks = 0usize;
    let mut scratch_clean_checks = 0usize;
    for batch_start in (0..states).step_by(64) {
        let shots = (states - batch_start).min(64);
        let live = if shots == 64 {
            u64::MAX
        } else {
            (1u64 << shots) - 1
        };
        let mut baseline_seed = Shake128::default();
        baseline_seed.update(b"coefficient-nonnegative-x-cancel");
        baseline_seed.update(&(batch_start as u64).to_le_bytes());
        let mut baseline_xof = baseline_seed.clone().finalize_xof();
        let mut optimized_xof = baseline_seed.finalize_xof();
        let mut baseline_simulator = Simulator::new(
            baseline.builder.next_qubit as usize,
            baseline.builder.next_bit as usize,
            &mut baseline_xof,
        );
        let mut optimized_simulator = Simulator::new(
            optimized.builder.next_qubit as usize,
            optimized.builder.next_bit as usize,
            &mut optimized_xof,
        );
        for shot in 0..shots {
            let value = (batch_start + shot) as u64;
            for (bit, (&baseline_id, &optimized_id)) in baseline
                .data_ids
                .iter()
                .zip(&optimized.data_ids)
                .enumerate()
            {
                if (value >> bit) & 1 != 0 {
                    *baseline_simulator.qubit_mut(QubitId(u64::from(baseline_id))) |=
                        1u64 << shot;
                    *optimized_simulator.qubit_mut(QubitId(u64::from(optimized_id))) |=
                        1u64 << shot;
                }
            }
        }
        baseline_simulator.apply_iter(baseline.builder.ops.iter());
        optimized_simulator.apply_iter(optimized.builder.ops.iter());
        assert_eq!(baseline_simulator.phase & live, 0);
        assert_eq!(optimized_simulator.phase & live, 0);
        phase_clean_checks += 2 * shots;
        for id in 0..baseline.builder.next_qubit {
            assert_eq!(
                baseline_simulator.qubit(QubitId(u64::from(id))) & live,
                optimized_simulator.qubit(QubitId(u64::from(id))) & live
            );
        }
        assert_eq!(
            baseline_simulator.qubit(QubitId(u64::from(baseline.tmp_id))) & live,
            0
        );
        assert_eq!(
            optimized_simulator.qubit(QubitId(u64::from(optimized.tmp_id))) & live,
            0
        );
        scratch_clean_checks += 2 * shots;
        cases_checked += shots;
    }
    (cases_checked, phase_clean_checks, scratch_clean_checks)
}

/// Exhaustively prove that the two sign-bit X gates surrounding a
/// source-disjoint coefficient body cancel, then derive the exact production
/// saving from the sealed active-window schedule.
#[doc(hidden)]
#[must_use]
pub fn exhaustive_coefficient_nonnegative_x_cancel_check(
) -> CoefficientNonnegativeXCancelProofReport {
    use crate::circuit::OperationType;

    assert_ne!(
        std::env::var(COEFFICIENT_NONNEGATIVE_X_CANCEL_FLAG)
            .ok()
            .as_deref(),
        Some("1"),
        "coefficient nonnegative X cancellation must default off"
    );
    let saved = std::env::var_os(COEFFICIENT_NONNEGATIVE_X_CANCEL_FLAG);
    let cursor_widths = 1usize..=4;
    let mut configurations_checked = 0usize;
    let mut basis_states_checked = 0usize;
    let mut scalar_equivalence_checks = 0usize;
    let mut simulator_equivalence_checks = 0usize;
    let mut phase_clean_checks = 0usize;
    let mut inverse_pair_checks = 0usize;
    let mut cursor_restore_checks = 0usize;
    let mut active_restore_checks = 0usize;
    let mut scratch_clean_checks = 0usize;
    let mut control_off_identity_checks = 0usize;
    let mut default_stream_identity_checks = 0usize;
    let mut non_x_kind_identity_checks = 0usize;

    for cursor_width in cursor_widths.clone() {
        for body in COEFFICIENT_BRACKET_BODIES {
            let default = build_coefficient_bracket_harness(cursor_width, body, None);
            let baseline = build_coefficient_bracket_harness(cursor_width, body, Some(false));
            let optimized = build_coefficient_bracket_harness(cursor_width, body, Some(true));
            let inverse =
                build_coefficient_bracket_harness(cursor_width, body.inverse(), Some(true));
            assert_coefficient_bracket_stream_identity(&default, &baseline);
            default_stream_identity_checks += 1;
            assert_eq!(baseline.data_ids, optimized.data_ids);
            assert_eq!(baseline.builder.active_qubits, optimized.builder.active_qubits);
            assert_eq!(baseline.builder.peak_qubits, optimized.builder.peak_qubits);
            assert_eq!(baseline.builder.ops.len() - optimized.builder.ops.len(), 2);
            assert_eq!(
                baseline.builder.counted_kind_ops[OperationType::X as usize]
                    - optimized.builder.counted_kind_ops[OperationType::X as usize],
                2
            );
            for kind in 0..baseline.builder.counted_kind_ops.len() {
                if kind != OperationType::X as usize {
                    assert_eq!(
                        baseline.builder.counted_kind_ops[kind],
                        optimized.builder.counted_kind_ops[kind]
                    );
                    non_x_kind_identity_checks += 1;
                }
            }
            let (simulator_cases, phases, scratch) =
                verify_coefficient_bracket_simulator_equivalence(&baseline, &optimized);
            simulator_equivalence_checks += simulator_cases;
            phase_clean_checks += phases;
            scratch_clean_checks += scratch;

            let states = 1u64 << baseline.data_ids.len();
            for value in 0..states {
                let input = coefficient_bracket_input(&baseline, value);
                let baseline_output = apply_scalar(&baseline.builder.ops, input);
                let optimized_output = apply_scalar(&optimized.builder.ops, input);
                assert_eq!(optimized_output, baseline_output);
                assert_eq!(optimized_output & (1u64 << optimized.tmp_id), 0);
                assert_eq!(
                    optimized_output & optimized.cursor_mask,
                    input & optimized.cursor_mask
                );
                assert_eq!(
                    optimized_output & optimized.active_mask,
                    input & optimized.active_mask
                );
                assert_eq!(apply_scalar(&inverse.builder.ops, optimized_output), input);
                let control = value & 1;
                let active_bit = (value >> (cursor_width + 1)) & 1;
                if control == 0 && active_bit == 0 {
                    assert_eq!(optimized_output, input);
                    control_off_identity_checks += 1;
                }
                basis_states_checked += 1;
                scalar_equivalence_checks += 1;
                inverse_pair_checks += 1;
                cursor_restore_checks += 1;
                active_restore_checks += 1;
                scratch_clean_checks += 1;
            }
            configurations_checked += 1;
        }
    }

    match saved {
        Some(value) => std::env::set_var(COEFFICIENT_NONNEGATIVE_X_CANCEL_FLAG, value),
        None => std::env::remove_var(COEFFICIENT_NONNEGATIVE_X_CANCEL_FLAG),
    }
    let schedule = exhaustive_reference_schedule_check();
    let active_coefficient_positions = schedule.t_window_sum;
    let bracket_pairs_per_phase_block_position = 6usize;
    let coefficient_phase_blocks_per_point_add = 4usize;
    let removed_x_per_position =
        2 * bracket_pairs_per_phase_block_position * coefficient_phase_blocks_per_point_add;
    let total_removed_x = active_coefficient_positions * removed_x_per_position;
    assert_eq!(active_coefficient_positions, 249_543);
    assert_eq!(removed_x_per_position, 48);
    assert_eq!(total_removed_x, 11_978_064);

    CoefficientNonnegativeXCancelProofReport {
        cursor_widths_checked: cursor_widths.count(),
        body_kinds_checked: COEFFICIENT_BRACKET_BODIES.len(),
        configurations_checked,
        basis_states_checked,
        scalar_equivalence_checks,
        simulator_equivalence_checks,
        phase_clean_checks,
        inverse_pair_checks,
        cursor_restore_checks,
        active_restore_checks,
        scratch_clean_checks,
        control_off_identity_checks,
        default_stream_identity_checks,
        non_x_kind_identity_checks,
        local_ops_delta: -2,
        local_x_delta: -2,
        local_toffoli_delta: 0,
        scheduled_steps: schedule.steps_checked,
        active_coefficient_positions,
        bracket_pairs_per_phase_block_position,
        coefficient_phase_blocks_per_point_add,
        removed_x_per_position,
        total_removed_x,
    }
}

struct Q845FusionGuardHarness {
    builder: B,
    control_id: u32,
    target_length_ids: Vec<u32>,
    source_length_ids: Vec<u32>,
    shift_ids: Vec<u32>,
    target_id: u32,
    data_mask: u64,
}

struct Q845FusionCoreHarness {
    builder: B,
    phase1_id: u32,
    phase2_id: u32,
    sign_id: u32,
    work1_ids: Vec<u32>,
    work2_ids: Vec<u32>,
    length_ids: Vec<u32>,
    above_guard_id: u32,
    data_mask: u64,
}

fn q845_fusion_id_mask(ids: &[u32]) -> u64 {
    ids.iter().fold(0u64, |mask, &id| mask | (1u64 << id))
}

fn q845_fusion_proof_circuit() -> Circuit {
    let mut circ = Circuit::new();
    circ.b.count_only = false;
    circ.b.fiat_hash = None;
    circ
}

fn apply_q845_fusion_classical_with_clean_resets(
    ops: &[crate::circuit::Op],
    mut state: u64,
) -> u64 {
    use crate::circuit::OperationType;

    let bit = |word: u64, id: u64| ((word >> id) & 1) != 0;
    for operation in ops {
        match operation.kind {
            OperationType::X => state ^= 1u64 << operation.q_target.0,
            OperationType::CX => {
                if bit(state, operation.q_control1.0) {
                    state ^= 1u64 << operation.q_target.0;
                }
            }
            OperationType::CCX => {
                if bit(state, operation.q_control1.0)
                    && bit(state, operation.q_control2.0)
                {
                    state ^= 1u64 << operation.q_target.0;
                }
            }
            OperationType::R => assert!(
                !bit(state, operation.q_target.0),
                "Q845 lifetime fusion reset dirty q{}",
                operation.q_target.0
            ),
            other => panic!("Q845 lifetime fusion emitted nonclassical operation {other:?}"),
        }
    }
    state
}

fn q845_fusion_set_register(mut value: u64, ids: &[u32], register: usize) -> u64 {
    for (index, &id) in ids.iter().enumerate() {
        let mask = 1u64 << id;
        value &= !mask;
        if register & (1usize << index) != 0 {
            value |= mask;
        }
    }
    value
}

fn q845_fusion_read_register(value: u64, ids: &[u32]) -> usize {
    ids.iter().enumerate().fold(0usize, |result, (index, id)| {
        result | ((((value >> id) & 1) as usize) << index)
    })
}

fn build_q845_fusion_guard_harness(length_width: usize) -> Q845FusionGuardHarness {
    let mut circ = q845_fusion_proof_circuit();
    let control = circ.alloc_qreg("q845.coeff-fusion-proof.guard.control");
    let target_length =
        circ.alloc_qreg_bits("q845.coeff-fusion-proof.guard.target-length", length_width);
    let source_length =
        circ.alloc_qreg_bits("q845.coeff-fusion-proof.guard.source-length", length_width);
    let shift = circ.alloc_qreg_bits("q845.coeff-fusion-proof.guard.shift", length_width);
    let target = circ.alloc_qreg("q845.coeff-fusion-proof.guard.target");
    let scratch = circ.alloc_qreg_bits("q845.coeff-fusion-proof.guard.scratch", 3);
    let scratch_refs = scratch.iter().collect::<Vec<_>>();
    toggle_q845_coefficient_length_above_guarded_boundary(
        &mut circ,
        &control,
        &target_length,
        &source_length,
        &shift,
        &target,
        &scratch_refs,
    );
    free_clean(&mut circ, scratch);
    let control_id = control.id();
    let target_length_ids = target_length.iter().map(QReg::id).collect::<Vec<_>>();
    let source_length_ids = source_length.iter().map(QReg::id).collect::<Vec<_>>();
    let shift_ids = shift.iter().map(QReg::id).collect::<Vec<_>>();
    let target_id = target.id();
    let data_mask = (1u64 << control_id)
        | q845_fusion_id_mask(&target_length_ids)
        | q845_fusion_id_mask(&source_length_ids)
        | q845_fusion_id_mask(&shift_ids)
        | (1u64 << target_id);
    Q845FusionGuardHarness {
        builder: circ.into_builder(),
        control_id,
        target_length_ids,
        source_length_ids,
        shift_ids,
        target_id,
        data_mask,
    }
}

fn build_q845_fusion_core_harness(
    work_width: usize,
    length_width: usize,
    inverse: bool,
) -> Q845FusionCoreHarness {
    let mut circ = q845_fusion_proof_circuit();
    let phase1 = circ.alloc_qreg("q845.coeff-fusion-proof.core.phase1");
    let phase2 = circ.alloc_qreg("q845.coeff-fusion-proof.core.phase2");
    let sign = circ.alloc_qreg("q845.coeff-fusion-proof.core.sign");
    let work1 = circ.alloc_qreg_bits("q845.coeff-fusion-proof.core.work1", work_width);
    let work2 = circ.alloc_qreg_bits("q845.coeff-fusion-proof.core.work2", work_width);
    let l_t = circ.alloc_qreg_bits("q845.coeff-fusion-proof.core.l-t", length_width);
    let above_guard = circ.alloc_qreg("q845.coeff-fusion-proof.core.above-guard");
    let enable = circ.alloc_qreg("q845.coeff-fusion-proof.core.enable");
    let add_only = circ.alloc_qreg("q845.coeff-fusion-proof.core.add-only");
    let chain = circ.alloc_qreg_bits("q845.coeff-fusion-proof.core.chain", 2);

    if inverse {
        coefficient_fused_data_and_sign_q845(
            &mut circ,
            &phase1,
            &phase2,
            &sign,
            &work1,
            &work2,
            &l_t,
            &enable,
            &above_guard,
            &add_only,
            &chain,
            true,
        );
    } else {
        toggle_initial_coefficient_enable(
            &mut circ,
            &phase1,
            &phase2,
            &sign,
            &enable,
            &chain,
        );
        coefficient_fused_data_and_sign_q845(
            &mut circ,
            &phase1,
            &phase2,
            &sign,
            &work1,
            &work2,
            &l_t,
            &enable,
            &above_guard,
            &add_only,
            &chain,
            false,
        );
    }

    free_clean(&mut circ, chain);
    circ.zero_and_free(add_only);
    circ.zero_and_free(enable);
    let phase1_id = phase1.id();
    let phase2_id = phase2.id();
    let sign_id = sign.id();
    let work1_ids = work1.iter().map(QReg::id).collect::<Vec<_>>();
    let work2_ids = work2.iter().map(QReg::id).collect::<Vec<_>>();
    let length_ids = l_t.iter().map(QReg::id).collect::<Vec<_>>();
    let above_guard_id = above_guard.id();
    let data_mask = (1u64 << phase1_id)
        | (1u64 << phase2_id)
        | (1u64 << sign_id)
        | q845_fusion_id_mask(&work1_ids)
        | q845_fusion_id_mask(&work2_ids)
        | q845_fusion_id_mask(&length_ids)
        | (1u64 << above_guard_id);
    Q845FusionCoreHarness {
        builder: circ.into_builder(),
        phase1_id,
        phase2_id,
        sign_id,
        work1_ids,
        work2_ids,
        length_ids,
        above_guard_id,
        data_mask,
    }
}

fn set_q845_lifetime_fusion_proof_mode(enabled: Option<bool>) {
    match enabled {
        Some(true) => std::env::set_var(Q845_LIFETIME_COEFFICIENT_FUSION_FLAG, "1"),
        Some(false) => std::env::set_var(Q845_LIFETIME_COEFFICIENT_FUSION_FLAG, "0"),
        None => std::env::remove_var(Q845_LIFETIME_COEFFICIENT_FUSION_FLAG),
    }
}

fn build_q845_fusion_dispatch_harness(enabled: Option<bool>, direct: bool) -> B {
    set_q845_lifetime_fusion_proof_mode(enabled);
    let mut circ = q845_fusion_proof_circuit();
    let phase1 = circ.alloc_qreg("q845.coeff-fusion-proof.dispatch.phase1");
    let phase2 = circ.alloc_qreg("q845.coeff-fusion-proof.dispatch.phase2");
    let sign = circ.alloc_qreg("q845.coeff-fusion-proof.dispatch.sign");
    let work1 = circ.alloc_qreg_bits("q845.coeff-fusion-proof.dispatch.work1", 4);
    let work2 = circ.alloc_qreg_bits("q845.coeff-fusion-proof.dispatch.work2", 4);
    let l_t = circ.alloc_qreg_bits("q845.coeff-fusion-proof.dispatch.l-t", 3);
    let l_t_prime = circ.alloc_qreg_bits("q845.coeff-fusion-proof.dispatch.l-t-prime", 3);
    let l_s = circ.alloc_qreg_bits("q845.coeff-fusion-proof.dispatch.l-s", 3);
    let l_r_prime = circ.alloc_qreg_bits("q845.coeff-fusion-proof.dispatch.l-r-prime", 3);
    if direct {
        coefficient_phase_block_fused_q845(
            &mut circ,
            &phase1,
            &phase2,
            &sign,
            &work1,
            &work2,
            &work2,
            &l_t,
            &l_t_prime,
            &l_s,
            &l_r_prime,
            false,
            None,
            None,
        );
    } else {
        coefficient_phase_block(
            &mut circ,
            &phase1,
            &phase2,
            &sign,
            &work1,
            &work2,
            &work2,
            &l_t,
            &l_t_prime,
            &l_s,
            &l_r_prime,
            false,
        );
    }
    circ.into_builder()
}

fn assert_q845_fusion_builder_identity(left: &B, right: &B) {
    assert_eq!(left.ops, right.ops);
    assert_eq!(left.next_qubit, right.next_qubit);
    assert_eq!(left.next_bit, right.next_bit);
    assert_eq!(left.active_qubits, right.active_qubits);
    assert_eq!(left.peak_qubits, right.peak_qubits);
    assert_eq!(left.free_qubits, right.free_qubits);
    assert_eq!(left.allocation_serial, right.allocation_serial);
}

struct Q845SwapOnlyCoreHarness {
    builder: B,
    phase1_id: u32,
    phase2_id: u32,
    sign_id: u32,
    work1_ids: Vec<u32>,
    work2_ids: Vec<u32>,
    l_t_ids: Vec<u32>,
    l_t_prime_ids: Vec<u32>,
    l_s_ids: Vec<u32>,
    l_r_prime_ids: Vec<u32>,
    above_guard_id: u32,
    data_mask: u64,
}

fn build_q845_swap_only_core_harness(
    physical_work_width: usize,
    scan_width: usize,
    coefficient_width: usize,
    length_width: usize,
    r_length_width: usize,
    inverse: bool,
    swap_only: bool,
    truncated_guard: bool,
) -> Q845SwapOnlyCoreHarness {
    assert!(scan_width <= physical_work_width);
    assert!(coefficient_width <= scan_width);
    assert_l_r_prime_metadata_width(length_width, r_length_width);
    assert!(!truncated_guard || swap_only);
    if truncated_guard {
        std::env::set_var(Q851_TRUNCATED_SWAP_ONLY_GUARD_FLAG, "1");
    } else {
        std::env::remove_var(Q851_TRUNCATED_SWAP_ONLY_GUARD_FLAG);
    }
    let mut circ = q845_fusion_proof_circuit();
    let phase1 = circ.alloc_qreg("q845.swap-only-proof.phase1");
    let phase2 = circ.alloc_qreg("q845.swap-only-proof.phase2");
    let sign = circ.alloc_qreg("q845.swap-only-proof.sign");
    let work1 = circ.alloc_qreg_bits("q845.swap-only-proof.work1", scan_width);
    let work2 = circ.alloc_qreg_bits("q845.swap-only-proof.work2", physical_work_width);
    let l_t = circ.alloc_qreg_bits("q845.swap-only-proof.l-t", length_width);
    let l_t_prime = circ.alloc_qreg_bits("q845.swap-only-proof.l-t-prime", length_width);
    let l_s = circ.alloc_qreg_bits("q845.swap-only-proof.l-s", length_width);
    let l_r_prime = circ.alloc_qreg_bits("q845.swap-only-proof.l-r-prime", r_length_width);
    let above_guard = circ.alloc_qreg("q845.swap-only-proof.above-guard");
    let enable = circ.alloc_qreg("q845.swap-only-proof.enable");
    let add_only = circ.alloc_qreg("q845.swap-only-proof.add-only");
    let chain = circ.alloc_qreg_bits("q845.swap-only-proof.chain", 2);

    if swap_only {
        coefficient_fused_data_and_sign_q845_swap_only(
            &mut circ,
            &phase1,
            &phase2,
            &sign,
            &work1[..coefficient_width],
            &work2[..coefficient_width],
            physical_work_width,
            &l_t,
            &l_t_prime,
            &l_s,
            &l_r_prime,
            &enable,
            &above_guard,
            &add_only,
            &chain,
            inverse,
            None,
            None,
        );
    } else if inverse {
        coefficient_fused_data_and_sign_q845(
            &mut circ,
            &phase1,
            &phase2,
            &sign,
            &work1[..coefficient_width],
            &work2[..coefficient_width],
            &l_t,
            &enable,
            &above_guard,
            &add_only,
            &chain,
            true,
        );
    } else {
        toggle_initial_coefficient_enable(
            &mut circ,
            &phase1,
            &phase2,
            &sign,
            &enable,
            &chain,
        );
        coefficient_fused_data_and_sign_q845(
            &mut circ,
            &phase1,
            &phase2,
            &sign,
            &work1[..coefficient_width],
            &work2[..coefficient_width],
            &l_t,
            &enable,
            &above_guard,
            &add_only,
            &chain,
            false,
        );
    }

    free_clean(&mut circ, chain);
    circ.zero_and_free(add_only);
    circ.zero_and_free(enable);
    let phase1_id = phase1.id();
    let phase2_id = phase2.id();
    let sign_id = sign.id();
    let work1_ids = work1.iter().map(QReg::id).collect::<Vec<_>>();
    let work2_ids = work2.iter().map(QReg::id).collect::<Vec<_>>();
    let l_t_ids = l_t.iter().map(QReg::id).collect::<Vec<_>>();
    let l_t_prime_ids = l_t_prime.iter().map(QReg::id).collect::<Vec<_>>();
    let l_s_ids = l_s.iter().map(QReg::id).collect::<Vec<_>>();
    let l_r_prime_ids = l_r_prime.iter().map(QReg::id).collect::<Vec<_>>();
    let above_guard_id = above_guard.id();
    let data_mask = (1u64 << phase1_id)
        | (1u64 << phase2_id)
        | (1u64 << sign_id)
        | q845_fusion_id_mask(&work1_ids)
        | q845_fusion_id_mask(&work2_ids)
        | q845_fusion_id_mask(&l_t_ids)
        | q845_fusion_id_mask(&l_t_prime_ids)
        | q845_fusion_id_mask(&l_s_ids)
        | q845_fusion_id_mask(&l_r_prime_ids)
        | (1u64 << above_guard_id);
    Q845SwapOnlyCoreHarness {
        builder: circ.into_builder(),
        phase1_id,
        phase2_id,
        sign_id,
        work1_ids,
        work2_ids,
        l_t_ids,
        l_t_prime_ids,
        l_s_ids,
        l_r_prime_ids,
        above_guard_id,
        data_mask,
    }
}

fn build_q845_swap_only_dispatch_harness(enabled: Option<bool>) -> B {
    match enabled {
        Some(true) => std::env::set_var(Q845_SWAP_ONLY_T_PRIME_LENGTH_FLAG, "1"),
        Some(false) => std::env::set_var(Q845_SWAP_ONLY_T_PRIME_LENGTH_FLAG, "0"),
        None => std::env::remove_var(Q845_SWAP_ONLY_T_PRIME_LENGTH_FLAG),
    }
    std::env::set_var(Q845_LIFETIME_COEFFICIENT_FUSION_FLAG, "1");
    let mut circ = q845_fusion_proof_circuit();
    let phase1 = circ.alloc_qreg("q845.swap-only-dispatch.phase1");
    let phase2 = circ.alloc_qreg("q845.swap-only-dispatch.phase2");
    let sign = circ.alloc_qreg("q845.swap-only-dispatch.sign");
    let work1 = circ.alloc_qreg_bits("q845.swap-only-dispatch.work1", 4);
    let work2 = circ.alloc_qreg_bits("q845.swap-only-dispatch.work2", 5);
    let l_t = circ.alloc_qreg_bits("q845.swap-only-dispatch.l-t", 3);
    let l_t_prime = circ.alloc_qreg_bits("q845.swap-only-dispatch.l-t-prime", 3);
    let l_s = circ.alloc_qreg_bits("q845.swap-only-dispatch.l-s", 3);
    let l_r_prime = circ.alloc_qreg_bits("q845.swap-only-dispatch.l-r-prime", 2);
    coefficient_phase_block(
        &mut circ,
        &phase1,
        &phase2,
        &sign,
        &work1,
        &work2[..4],
        &work2,
        &l_t,
        &l_t_prime,
        &l_s,
        &l_r_prime,
        false,
    );
    circ.into_builder()
}

struct Q845EphemeralSwapHarness {
    builder: B,
    data_ids: Vec<u32>,
    external_mask: u64,
    iteration_id: u32,
    work1_ids: Vec<u32>,
    work2_ids: Vec<u32>,
    l_t_ids: Vec<u32>,
    l_t_prime_ids: Vec<u32>,
    l_q_ids: Vec<u32>,
    l_s_ids: Vec<u32>,
    l_r_prime_ids: Vec<u32>,
}

fn build_q845_ephemeral_swap_harness(
    work_width: usize,
    length_width: usize,
    r_length_width: usize,
    inverse: bool,
    swap_only: bool,
    fuse_support_lifetime: bool,
) -> Q845EphemeralSwapHarness {
    if swap_only {
        std::env::set_var(Q845_SWAP_ONLY_T_PRIME_LENGTH_FLAG, "1");
    } else {
        std::env::remove_var(Q845_SWAP_ONLY_T_PRIME_LENGTH_FLAG);
    }
    std::env::set_var(PROMISED_LQ_SWAP_BORROW_FLAG, "1");
    if fuse_support_lifetime {
        std::env::set_var(PROMISED_SWAP_SUPPORT_LIFETIME_FUSION_FLAG, "1");
    } else {
        std::env::remove_var(PROMISED_SWAP_SUPPORT_LIFETIME_FUSION_FLAG);
    }
    let mut circ = q845_fusion_proof_circuit();
    let iteration = circ.alloc_qreg("q845.ephemeral-swap.iteration");
    let work1 = circ.alloc_qreg_bits("q845.ephemeral-swap.work1", work_width);
    let work2 = circ.alloc_qreg_bits("q845.ephemeral-swap.work2", work_width);
    let l_t = circ.alloc_qreg_bits("q845.ephemeral-swap.l-t", length_width);
    let l_t_prime = circ.alloc_qreg_bits("q845.ephemeral-swap.l-t-prime", length_width);
    let l_q = circ.alloc_qreg_bits("q845.ephemeral-swap.l-q", length_width);
    let l_s = circ.alloc_qreg_bits("q845.ephemeral-swap.l-s", length_width);
    let l_r_prime = circ.alloc_qreg_bits("q845.ephemeral-swap.l-r-prime", r_length_width);
    let predicate_chain_width = 3usize.max(length_width.saturating_sub(2));
    let condition_scratch =
        circ.alloc_qreg_bits("q845.ephemeral-swap.condition", 3 + predicate_chain_width);
    let zero_q = &condition_scratch[0];
    let zero_s = &condition_scratch[1];
    let control = &condition_scratch[2];
    let chain = &condition_scratch[3..];
    compute_zero(&mut circ, &l_q, zero_q, chain);
    compute_zero(&mut circ, &l_s, zero_s, chain);
    conditional_work_and_length_swap_under_zero_predicate(
        &mut circ,
        zero_q,
        zero_s,
        control,
        &iteration,
        &work1,
        &work2,
        &l_t,
        &l_t_prime,
        &l_q,
        &l_s,
        &l_r_prime,
        (1, work_width),
        (1, work_width),
        chain,
        &[],
        inverse,
        PromisedLqSwapRoute::Configured,
    );
    uncompute_zero(&mut circ, &l_s, zero_s, chain);
    uncompute_zero(&mut circ, &l_q, zero_q, chain);
    free_clean(&mut circ, condition_scratch);

    let iteration_id = iteration.id();
    let work1_ids = work1.iter().map(QReg::id).collect::<Vec<_>>();
    let work2_ids = work2.iter().map(QReg::id).collect::<Vec<_>>();
    let l_t_ids = l_t.iter().map(QReg::id).collect::<Vec<_>>();
    let l_t_prime_ids = l_t_prime.iter().map(QReg::id).collect::<Vec<_>>();
    let l_q_ids = l_q.iter().map(QReg::id).collect::<Vec<_>>();
    let l_s_ids = l_s.iter().map(QReg::id).collect::<Vec<_>>();
    let l_r_prime_ids = l_r_prime.iter().map(QReg::id).collect::<Vec<_>>();
    let data_ids = std::iter::once(iteration_id)
        .chain(work1_ids.iter().copied())
        .chain(work2_ids.iter().copied())
        .chain(l_t_ids.iter().copied())
        .chain(l_t_prime_ids.iter().copied())
        .chain(l_q_ids.iter().copied())
        .chain(l_s_ids.iter().copied())
        .chain(l_r_prime_ids.iter().copied())
        .collect::<Vec<_>>();
    let external_mask = data_ids
        .iter()
        .fold(0u64, |mask, id| mask | (1u64 << id));
    Q845EphemeralSwapHarness {
        builder: circ.into_builder(),
        data_ids,
        external_mask,
        iteration_id,
        work1_ids,
        work2_ids,
        l_t_ids,
        l_t_prime_ids,
        l_q_ids,
        l_s_ids,
        l_r_prime_ids,
    }
}

fn q845_test_bit_length(value: usize) -> usize {
    if value == 0 {
        0
    } else {
        usize::BITS as usize - value.leading_zeros() as usize
    }
}

fn q845_test_pack_work1(width: usize, t: usize, r: usize) -> usize {
    let l_t = q845_test_bit_length(t);
    let l_r = q845_test_bit_length(r);
    assert!(l_t + 1 + l_r <= width);
    let mut packed = t;
    for bit in 0..l_r {
        if (r >> bit) & 1 != 0 {
            packed |= 1usize << (width - 1 - bit);
        }
    }
    packed
}

fn q845_test_pack_work2(width: usize, t_prime: usize, r_prime: usize) -> usize {
    let l_t_prime = q845_test_bit_length(t_prime);
    let l_r_prime = q845_test_bit_length(r_prime);
    assert!(l_t_prime + l_r_prime <= width);
    let mut packed = t_prime;
    for bit in 0..l_r_prime {
        if (r_prime >> bit) & 1 != 0 {
            packed |= 1usize << (width - 1 - bit);
        }
    }
    packed
}

fn q845_pack_data_value(data_ids: &[u32], state: u64) -> usize {
    data_ids
        .iter()
        .enumerate()
        .fold(0usize, |packed, (bit, id)| {
            packed | ((((state >> id) & 1) as usize) << bit)
        })
}

struct Q851RangeComparatorHarness {
    builder: B,
    register_ids: Vec<u32>,
    target_id: u32,
    external_mask: u64,
}

fn build_q851_range_comparator_harness(
    width: usize,
    value: usize,
) -> Q851RangeComparatorHarness {
    let mut circ = q845_fusion_proof_circuit();
    let register = circ.alloc_qreg_bits("q851.range-proof.register", width);
    let target = circ.alloc_qreg("q851.range-proof.target");
    let scratch = circ.alloc_qreg_bits("q851.range-proof.scratch", width.saturating_sub(2));
    toggle_register_geq_constant_vchain(&mut circ, &register, value, &target, &scratch);
    free_clean(&mut circ, scratch);
    let register_ids = register.iter().map(QReg::id).collect::<Vec<_>>();
    let target_id = target.id();
    let external_mask = q845_fusion_id_mask(&register_ids) | (1u64 << target_id);
    Q851RangeComparatorHarness {
        builder: circ.into_builder(),
        register_ids,
        target_id,
        external_mask,
    }
}

struct Q851TruncatedGuardTraceHarness {
    builder: B,
    l_s_ids: Vec<u32>,
    l_r_prime_ids: Vec<u32>,
    forward_trace_ids: Vec<u32>,
    reverse_trace_ids: Vec<u32>,
    external_mask: u64,
}

fn build_q851_truncated_guard_trace_harness(
    physical_work_width: usize,
    scan_width: usize,
    length_width: usize,
    truncated: bool,
) -> Q851TruncatedGuardTraceHarness {
    if truncated {
        std::env::set_var(Q851_TRUNCATED_SWAP_ONLY_GUARD_FLAG, "1");
    } else {
        std::env::remove_var(Q851_TRUNCATED_SWAP_ONLY_GUARD_FLAG);
    }
    let mut circ = q845_fusion_proof_circuit();
    let count = circ.alloc_qreg_bits("q851.guard-trace.count", length_width);
    let l_s = circ.alloc_qreg_bits("q851.guard-trace.l-s", length_width);
    let l_r_prime = circ.alloc_qreg_bits("q851.guard-trace.l-r-prime", length_width);
    let carry = circ.alloc_qreg("q851.guard-trace.carry");
    let overflow = circ.alloc_qreg("q851.guard-trace.overflow");
    let active = circ.alloc_qreg("q851.guard-trace.active");
    let forward_trace = circ.alloc_qreg_bits("q851.guard-trace.forward", scan_width);
    let reverse_trace = circ.alloc_qreg_bits("q851.guard-trace.reverse", scan_width);

    let guard = Q845SwapOnlyCoefficientGuard::prepare(
        &mut circ,
        physical_work_width,
        &count,
        &l_s,
        &l_r_prime,
        &carry,
        &overflow,
        None,
        None,
    );
    guard.for_each_forward(
        &mut circ,
        scan_width,
        &active,
        &count,
        |circ, index, guard_active| circ.cx(guard_active, &forward_trace[index]),
    );
    let constant_scratch = count.iter().collect::<Vec<_>>();
    guard.prepare_reverse_boundary(
        &mut circ,
        &constant_scratch,
        &l_s,
        &l_r_prime,
        &carry,
        &active,
    );
    guard.for_each_reverse(
        &mut circ,
        scan_width,
        &active,
        &count,
        &constant_scratch,
        &carry,
        |circ, index, guard_active| circ.cx(guard_active, &reverse_trace[index]),
    );
    guard.finish(
        &mut circ,
        &count,
        &l_s,
        &l_r_prime,
        &carry,
        &overflow,
    );
    free_clean(&mut circ, count);
    circ.zero_and_free(active);
    circ.zero_and_free(overflow);
    circ.zero_and_free(carry);

    let l_s_ids = l_s.iter().map(QReg::id).collect::<Vec<_>>();
    let l_r_prime_ids = l_r_prime.iter().map(QReg::id).collect::<Vec<_>>();
    let forward_trace_ids = forward_trace.iter().map(QReg::id).collect::<Vec<_>>();
    let reverse_trace_ids = reverse_trace.iter().map(QReg::id).collect::<Vec<_>>();
    let external_mask = q845_fusion_id_mask(&l_s_ids)
        | q845_fusion_id_mask(&l_r_prime_ids)
        | q845_fusion_id_mask(&forward_trace_ids)
        | q845_fusion_id_mask(&reverse_trace_ids);
    Q851TruncatedGuardTraceHarness {
        builder: circ.into_builder(),
        l_s_ids,
        l_r_prime_ids,
        forward_trace_ids,
        reverse_trace_ids,
        external_mask,
    }
}

struct Q851FixedSignEventHarness {
    builder: B,
    cursor_ids: Vec<u32>,
    scratch_ids: Vec<u32>,
    allocation_free_transitions: usize,
}

fn build_q851_fixed_sign_event_harness(
    transition_index: usize,
    direction: Q851CoefficientCursorDirection,
) -> Q851FixedSignEventHarness {
    let mut circ = q845_fusion_proof_circuit();
    let cursor = circ.alloc_qreg_bits("q851.fixed-sign-proof.cursor", REFERENCE_LENGTH_WIDTH);
    let scratch = circ.alloc_qreg_bits(
        "q851.fixed-sign-proof.scratch",
        REFERENCE_LENGTH_WIDTH - 1,
    );
    let before = circ.b.next_qubit;
    transition_q851_coefficient_cursor(
        &mut circ,
        &cursor,
        257,
        transition_index,
        &scratch,
        direction,
    );
    assert_eq!(circ.b.next_qubit, before, "fixed sign event allocated a qubit");
    Q851FixedSignEventHarness {
        builder: circ.into_builder(),
        cursor_ids: cursor.iter().map(QReg::id).collect(),
        scratch_ids: scratch.iter().map(QReg::id).collect(),
        allocation_free_transitions: 1,
    }
}

fn build_q851_fixed_sign_sequence_harness(
    width: usize,
    direction: Q851CoefficientCursorDirection,
) -> Q851FixedSignEventHarness {
    assert!((1..=257).contains(&width));
    let mut circ = q845_fusion_proof_circuit();
    let cursor = circ.alloc_qreg_bits("q851.fixed-sign-sequence.cursor", REFERENCE_LENGTH_WIDTH);
    let scratch = circ.alloc_qreg_bits(
        "q851.fixed-sign-sequence.scratch",
        REFERENCE_LENGTH_WIDTH - 1,
    );
    let mut allocation_free_transitions = 0usize;
    match direction {
        Q851CoefficientCursorDirection::Decrement => {
            for transition_index in 0..width - 1 {
                let before = circ.b.next_qubit;
                transition_q851_coefficient_cursor(
                    &mut circ,
                    &cursor,
                    width,
                    transition_index,
                    &scratch,
                    direction,
                );
                assert_eq!(circ.b.next_qubit, before, "forward sign event allocated a qubit");
                allocation_free_transitions += 1;
            }
        }
        Q851CoefficientCursorDirection::Increment => {
            for transition_index in (0..width - 1).rev() {
                let before = circ.b.next_qubit;
                transition_q851_coefficient_cursor(
                    &mut circ,
                    &cursor,
                    width,
                    transition_index,
                    &scratch,
                    direction,
                );
                assert_eq!(circ.b.next_qubit, before, "reverse sign event allocated a qubit");
                allocation_free_transitions += 1;
            }
        }
    }
    Q851FixedSignEventHarness {
        builder: circ.into_builder(),
        cursor_ids: cursor.iter().map(QReg::id).collect(),
        scratch_ids: scratch.iter().map(QReg::id).collect(),
        allocation_free_transitions,
    }
}

fn q851_fixed_sign_cursor_value(initial: usize, body_index: usize) -> usize {
    assert!(initial < 512);
    assert!(body_index <= 256);
    let low = initial & 255;
    let sign = ((initial >> 8) & 1) ^ usize::from(body_index > low);
    low | (sign << 8)
}

fn build_q851_fixed_sign_domain_harness(
    cursor_width: usize,
    scan_width: usize,
    direction: Q851CoefficientCursorDirection,
) -> B {
    assert!(cursor_width > 1);
    let mut circ = q845_fusion_proof_circuit();
    let cursor = circ.alloc_qreg_bits("q851.fixed-sign-domain.cursor", cursor_width);
    let scratch = circ.alloc_qreg_bits("q851.fixed-sign-domain.scratch", cursor_width - 1);
    let before = circ.b.next_qubit;
    transition_q851_coefficient_cursor(
        &mut circ,
        &cursor,
        scan_width,
        0,
        &scratch,
        direction,
    );
    assert_eq!(circ.b.next_qubit, before, "domain fallback allocated a qubit");
    circ.into_builder()
}

fn q851_apply_baseline_cursor_transition(
    value: usize,
    direction: Q851CoefficientCursorDirection,
) -> usize {
    match direction {
        Q851CoefficientCursorDirection::Decrement => value.wrapping_sub(1) & 511,
        Q851CoefficientCursorDirection::Increment => value.wrapping_add(1) & 511,
    }
}

fn q851_apply_fixed_sign_event(value: usize, transition_index: usize) -> usize {
    assert!(value < 512);
    assert!(transition_index < 256);
    if value & 255 == transition_index {
        value ^ 256
    } else {
        value
    }
}

fn q851_apply_cursor_transition_pair(
    baseline: &mut usize,
    candidate: &mut usize,
    traversal: Q851CoefficientCursorTraversal,
    index: usize,
    scan_width: usize,
) -> bool {
    let Some((transition_index, direction)) =
        q851_coefficient_cursor_transition(traversal, index, scan_width)
    else {
        return false;
    };
    assert!(transition_index < 256);
    *baseline = q851_apply_baseline_cursor_transition(*baseline, direction);
    *candidate = q851_apply_fixed_sign_event(*candidate, transition_index);
    true
}

fn verify_q851_fixed_sign_harness<Input, Expected>(
    label: &[u8],
    harness: &Q851FixedSignEventHarness,
    input_value: Input,
    expected_value: Expected,
) -> (usize, usize, usize, usize)
where
    Input: Fn(usize) -> usize,
    Expected: Fn(usize) -> usize,
{
    use crate::circuit::QubitId;
    use crate::sim::Simulator;
    use sha3::{
        digest::{ExtendableOutput, Update},
        Shake128,
    };

    assert_eq!(harness.cursor_ids.len(), REFERENCE_LENGTH_WIDTH);
    assert_eq!(harness.scratch_ids.len(), REFERENCE_LENGTH_WIDTH - 1);
    let mut cases_checked = 0usize;
    let mut scratch_clean_checks = 0usize;
    let mut phase_clean_checks = 0usize;
    let mut ancilla_clean_checks = 0usize;
    for batch_start in (0usize..512).step_by(64) {
        let mut seed = Shake128::default();
        seed.update(label);
        seed.update(&(batch_start as u64).to_le_bytes());
        let mut xof = seed.finalize_xof();
        let mut simulator = Simulator::new(
            harness.builder.next_qubit as usize,
            harness.builder.next_bit as usize,
            &mut xof,
        );
        for shot in 0..64 {
            let input = input_value(batch_start + shot);
            for (bit, id) in harness.cursor_ids.iter().enumerate() {
                if input & (1usize << bit) != 0 {
                    *simulator.qubit_mut(QubitId(u64::from(*id))) |= 1u64 << shot;
                }
            }
        }
        simulator.apply_iter(harness.builder.ops.iter());
        assert_eq!(simulator.phase, 0, "{label:?} left phase garbage");
        for (bit, id) in harness.cursor_ids.iter().enumerate() {
            let expected = (0..64).fold(0u64, |plane, shot| {
                let value = expected_value(batch_start + shot);
                plane | ((((value >> bit) & 1) as u64) << shot)
            });
            assert_eq!(
                simulator.qubit(QubitId(u64::from(*id))),
                expected,
                "{label:?} cursor bit {bit} mismatch at batch {batch_start}"
            );
        }
        for id in &harness.scratch_ids {
            assert_eq!(
                simulator.qubit(QubitId(u64::from(*id))),
                0,
                "{label:?} left scratch q{id} dirty"
            );
        }
        cases_checked += 64;
        scratch_clean_checks += 64;
        phase_clean_checks += 64;
        ancilla_clean_checks += 64;
    }
    (
        cases_checked,
        scratch_clean_checks,
        phase_clean_checks,
        ancilla_clean_checks,
    )
}

fn q851_baseline_transition9_counts() -> RegisterSharedGateCounts {
    let mut circ = q845_fusion_proof_circuit();
    let cursor = circ.alloc_qreg_bits("q851.fixed-sign-proof.baseline", REFERENCE_LENGTH_WIDTH);
    let scratch = circ.alloc_qreg_bits(
        "q851.fixed-sign-proof.baseline-scratch",
        REFERENCE_LENGTH_WIDTH - 1,
    );
    let before = circ.b.next_qubit;
    increment_mod_2n(&mut circ, &cursor, &scratch);
    assert_eq!(circ.b.next_qubit, before);
    gate_counts(&circ.into_builder().ops)
}

/// Exhaustively prove the fixed Q851 coefficient-sign event before enabling it
/// in a whole-route count. The proof keeps the low cursor byte fixed and checks
/// every body index against the exact modulo-512 identity.
#[doc(hidden)]
#[must_use]
pub fn exhaustive_q851_fixed_sign_event_check() -> Q851FixedSignEventProofReport {
    let saved = std::env::var_os(Q851_FIXED_SIGN_EVENT_FLAG);
    std::env::remove_var(Q851_FIXED_SIGN_EVENT_FLAG);
    let mut default_off_stream_identity_checks = 0usize;
    for direction in [
        Q851CoefficientCursorDirection::Decrement,
        Q851CoefficientCursorDirection::Increment,
    ] {
        let default = build_q851_fixed_sign_event_harness(0, direction);
        std::env::set_var(Q851_FIXED_SIGN_EVENT_FLAG, "0");
        let explicit_off = build_q851_fixed_sign_event_harness(0, direction);
        assert_q845_fusion_builder_identity(&default.builder, &explicit_off.builder);
        default_off_stream_identity_checks += 1;
        std::env::remove_var(Q851_FIXED_SIGN_EVENT_FLAG);
    }

    let mut identity_pairs_checked = 0usize;
    let mut body_256_checks = 0usize;
    for initial in 0usize..512 {
        for body_index in 0usize..=256 {
            let direct = ((initial + 512 - body_index) & 511) >> 8;
            let identity = q851_fixed_sign_cursor_value(initial, body_index) >> 8;
            assert_eq!(direct, identity);
            identity_pairs_checked += 1;
            body_256_checks += usize::from(body_index == 256);
        }
    }

    let mut production_schedule_widths_checked = 0usize;
    for step in 1..=REFERENCE_STEPS {
        let width = reference_active_windows(256, step).t_add_sub.1;
        assert!((1..=257).contains(&width));
        production_schedule_widths_checked += 1;
    }

    let mut domain_fallback_stream_identity_checks = 0usize;
    for direction in [
        Q851CoefficientCursorDirection::Decrement,
        Q851CoefficientCursorDirection::Increment,
    ] {
        for (cursor_width, scan_width) in [(5usize, 5usize), (9usize, 259usize)] {
            std::env::remove_var(Q851_FIXED_SIGN_EVENT_FLAG);
            let baseline =
                build_q851_fixed_sign_domain_harness(cursor_width, scan_width, direction);
            std::env::set_var(Q851_FIXED_SIGN_EVENT_FLAG, "1");
            let candidate =
                build_q851_fixed_sign_domain_harness(cursor_width, scan_width, direction);
            assert_q845_fusion_builder_identity(&baseline, &candidate);
            domain_fallback_stream_identity_checks += 1;
        }
    }

    std::env::set_var(Q851_FIXED_SIGN_EVENT_FLAG, "1");
    let mut transition_events_checked = 0usize;
    let mut transition_basis_states_checked = 0usize;
    let mut direction_stream_identity_checks = 0usize;
    let mut allocation_free_microkernels_checked = 0usize;
    let mut scratch_clean_checks = 0usize;
    let mut phase_clean_checks = 0usize;
    let mut ancilla_clean_checks = 0usize;
    let mut transition_x_min = usize::MAX;
    let mut transition_x_max = 0usize;
    let mut transition_ops_min = usize::MAX;
    let mut transition_ops_max = 0usize;
    let mut transition_toffoli = None;
    for transition_index in 0usize..256 {
        let forward = build_q851_fixed_sign_event_harness(
            transition_index,
            Q851CoefficientCursorDirection::Decrement,
        );
        let reverse = build_q851_fixed_sign_event_harness(
            transition_index,
            Q851CoefficientCursorDirection::Increment,
        );
        assert_q845_fusion_builder_identity(&forward.builder, &reverse.builder);
        direction_stream_identity_checks += 1;
        allocation_free_microkernels_checked +=
            forward.allocation_free_transitions + reverse.allocation_free_transitions;

        let counts = gate_counts(&forward.builder.ops);
        assert_eq!(counts.ccx, 13);
        assert_eq!(counts.cx, 0);
        assert_eq!(
            counts.x,
            2 * (REFERENCE_LENGTH_WIDTH - 1 - transition_index.count_ones() as usize)
        );
        match transition_toffoli {
            Some(expected) => assert_eq!(counts.ccx, expected),
            None => transition_toffoli = Some(counts.ccx),
        }
        transition_x_min = transition_x_min.min(counts.x);
        transition_x_max = transition_x_max.max(counts.x);
        transition_ops_min = transition_ops_min.min(counts.total);
        transition_ops_max = transition_ops_max.max(counts.total);

        let (cases, scratch, phase, ancilla) = verify_q851_fixed_sign_harness(
            b"q851-fixed-sign-event",
            &forward,
            |initial| initial,
            |initial| {
                if initial & 255 == transition_index {
                    initial ^ 256
                } else {
                    initial
                }
            },
        );
        transition_basis_states_checked += cases;
        scratch_clean_checks += scratch;
        phase_clean_checks += phase;
        ancilla_clean_checks += ancilla;
        transition_events_checked += 1;
    }

    let mut sequence_widths_checked = 0usize;
    let mut forward_sequence_cases_checked = 0usize;
    let mut reverse_sequence_cases_checked = 0usize;
    let mut exact_reverse_stream_checks = 0usize;
    let mut cursor_restore_checks = 0usize;
    let mut allocation_free_sequence_transitions_checked = 0usize;
    for width in 1usize..=257 {
        let forward = build_q851_fixed_sign_sequence_harness(
            width,
            Q851CoefficientCursorDirection::Decrement,
        );
        let reverse = build_q851_fixed_sign_sequence_harness(
            width,
            Q851CoefficientCursorDirection::Increment,
        );
        assert_eq!(forward.cursor_ids, reverse.cursor_ids);
        assert_eq!(forward.scratch_ids, reverse.scratch_ids);
        assert_eq!(
            reverse.builder.ops,
            forward.builder.ops.iter().rev().copied().collect::<Vec<_>>()
        );
        exact_reverse_stream_checks += 1;
        allocation_free_sequence_transitions_checked +=
            forward.allocation_free_transitions + reverse.allocation_free_transitions;
        let forward_counts = gate_counts(&forward.builder.ops);
        let reverse_counts = gate_counts(&reverse.builder.ops);
        assert_eq!(forward_counts, reverse_counts);
        assert_eq!(forward_counts.ccx, 13 * (width - 1));

        let body_index = width - 1;
        let (cases, scratch, phase, ancilla) = verify_q851_fixed_sign_harness(
            b"q851-fixed-sign-forward-sequence",
            &forward,
            |initial| initial,
            |initial| q851_fixed_sign_cursor_value(initial, body_index),
        );
        forward_sequence_cases_checked += cases;
        scratch_clean_checks += scratch;
        phase_clean_checks += phase;
        ancilla_clean_checks += ancilla;

        let (cases, scratch, phase, ancilla) = verify_q851_fixed_sign_harness(
            b"q851-fixed-sign-reverse-sequence",
            &reverse,
            |initial| q851_fixed_sign_cursor_value(initial, body_index),
            |initial| initial,
        );
        reverse_sequence_cases_checked += cases;
        cursor_restore_checks += cases;
        scratch_clean_checks += scratch;
        phase_clean_checks += phase;
        ancilla_clean_checks += ancilla;
        sequence_widths_checked += 1;
    }

    // The production coefficient body observes the cursor only through
    // `lower_negative`; the shared traversal function below also drives all
    // four production transition sites. Therefore equality of every observed
    // sign plane plus exact route-exit restoration is a compositional miter for
    // the unchanged guard, data, phase, and cleanup operations around it.
    let mut route_branch_cases_checked = 0usize;
    let mut route_transition_index_checks = 0usize;
    let mut route_body_sign_observation_checks = 0usize;
    let mut route_cursor_restore_checks = 0usize;
    for width in 1usize..=257 {
        for initial in 0usize..512 {
            let mut baseline = initial;
            let mut candidate = initial;
            for index in 0usize..width {
                route_transition_index_checks += usize::from(
                    q851_apply_cursor_transition_pair(
                        &mut baseline,
                        &mut candidate,
                        Q851CoefficientCursorTraversal::InverseForward,
                        index,
                        width,
                    ),
                );
                assert_eq!(baseline, (initial + 512 - index) & 511);
                assert_eq!(candidate >> 8, baseline >> 8);
                route_body_sign_observation_checks += 1;
            }
            for index in (0usize..width).rev() {
                route_transition_index_checks += usize::from(
                    q851_apply_cursor_transition_pair(
                        &mut baseline,
                        &mut candidate,
                        Q851CoefficientCursorTraversal::InverseReverse,
                        index,
                        width,
                    ),
                );
                assert_eq!(baseline, (initial + 512 - index) & 511);
                assert_eq!(candidate >> 8, baseline >> 8);
                route_body_sign_observation_checks += 1;
            }
            assert_eq!(baseline, initial);
            assert_eq!(candidate, initial);
            route_branch_cases_checked += 1;
            route_cursor_restore_checks += 1;

            let mut baseline = initial;
            let mut candidate = initial;
            for index in 0usize..width {
                assert_eq!(baseline, (initial + 512 - index) & 511);
                assert_eq!(candidate >> 8, baseline >> 8);
                route_body_sign_observation_checks += 1;
                route_transition_index_checks += usize::from(
                    q851_apply_cursor_transition_pair(
                        &mut baseline,
                        &mut candidate,
                        Q851CoefficientCursorTraversal::ForwardForward,
                        index,
                        width,
                    ),
                );
            }
            for index in (0usize..width).rev() {
                assert_eq!(baseline, (initial + 512 - index) & 511);
                assert_eq!(candidate >> 8, baseline >> 8);
                route_body_sign_observation_checks += 1;
                route_transition_index_checks += usize::from(
                    q851_apply_cursor_transition_pair(
                        &mut baseline,
                        &mut candidate,
                        Q851CoefficientCursorTraversal::ForwardReverse,
                        index,
                        width,
                    ),
                );
            }
            assert_eq!(baseline, initial);
            assert_eq!(candidate, initial);
            route_branch_cases_checked += 1;
            route_cursor_restore_checks += 1;
        }
    }

    let baseline_transition9 = q851_baseline_transition9_counts();
    assert_eq!(baseline_transition9.x, 1);
    assert_eq!(baseline_transition9.cx, 10);
    assert_eq!(baseline_transition9.ccx, 14);
    assert_eq!(baseline_transition9.total, 25);
    let transition_toffoli = transition_toffoli.expect("at least one transition event");
    assert_eq!(transition_x_min, 0);
    assert_eq!(transition_x_max, 16);
    assert_eq!(transition_ops_min, 13);
    assert_eq!(transition_ops_max, 29);

    match saved {
        Some(value) => std::env::set_var(Q851_FIXED_SIGN_EVENT_FLAG, value),
        None => std::env::remove_var(Q851_FIXED_SIGN_EVENT_FLAG),
    }

    Q851FixedSignEventProofReport {
        identity_pairs_checked,
        body_256_checks,
        production_schedule_widths_checked,
        domain_fallback_stream_identity_checks,
        route_branch_cases_checked,
        route_transition_index_checks,
        route_body_sign_observation_checks,
        route_cursor_restore_checks,
        transition_events_checked,
        transition_basis_states_checked,
        direction_stream_identity_checks,
        sequence_widths_checked,
        forward_sequence_cases_checked,
        reverse_sequence_cases_checked,
        exact_reverse_stream_checks,
        cursor_restore_checks,
        scratch_clean_checks,
        phase_clean_checks,
        ancilla_clean_checks,
        default_off_stream_identity_checks,
        allocation_free_microkernels_checked,
        allocation_free_sequence_transitions_checked,
        transition_toffoli,
        transition_x_min,
        transition_x_max,
        transition_ops_min,
        transition_ops_max,
        baseline_transition9,
    }
}

#[doc(hidden)]
#[must_use]
pub fn exhaustive_q845_swap_only_coefficient_check() -> Q845SwapOnlyCoefficientProofReport {
    let saved_swap_only = std::env::var_os(Q845_SWAP_ONLY_T_PRIME_LENGTH_FLAG);
    let saved_truncated_guard = std::env::var_os(Q851_TRUNCATED_SWAP_ONLY_GUARD_FLAG);
    let saved_inplace_guard = std::env::var_os(SUB800_INPLACE_GUARD_ADDRESS_FLAG);
    let saved_fusion = std::env::var_os(Q845_LIFETIME_COEFFICIENT_FUSION_FLAG);
    let saved_promised = std::env::var_os(PROMISED_LQ_SWAP_BORROW_FLAG);
    let saved_support_fusion =
        std::env::var_os(PROMISED_SWAP_SUPPORT_LIFETIME_FUSION_FLAG);
    let saved_paired_source = std::env::var_os(PAIRED_BITLEN_SOURCE_COMPLEMENT_FLAG);
    let saved_coefficient_x = std::env::var_os(COEFFICIENT_NONNEGATIVE_X_CANCEL_FLAG);
    let saved_direct_prefix = std::env::var_os("LOWQ_DIRECT_PREFIX_BITLEN");
    let saved_fused_zero_prefix = std::env::var_os("LOWQ_FUSED_ZERO_PREFIX_BITLEN");
    std::env::remove_var(PAIRED_BITLEN_SOURCE_COMPLEMENT_FLAG);
    std::env::remove_var(COEFFICIENT_NONNEGATIVE_X_CANCEL_FLAG);
    std::env::remove_var(Q851_TRUNCATED_SWAP_ONLY_GUARD_FLAG);
    std::env::remove_var(SUB800_INPLACE_GUARD_ADDRESS_FLAG);
    let default_stream = build_q845_swap_only_dispatch_harness(None);
    let explicit_off_stream = build_q845_swap_only_dispatch_harness(Some(false));
    assert_q845_fusion_builder_identity(&default_stream, &explicit_off_stream);

    std::env::remove_var(Q851_TRUNCATED_SWAP_ONLY_GUARD_FLAG);
    let default_guard_stream = build_q845_swap_only_dispatch_harness(Some(true));
    std::env::set_var(Q851_TRUNCATED_SWAP_ONLY_GUARD_FLAG, "0");
    let explicit_off_guard_stream = build_q845_swap_only_dispatch_harness(Some(true));
    assert_q845_fusion_builder_identity(&default_guard_stream, &explicit_off_guard_stream);
    let truncated_guard_default_off_stream_identity_checks = 1usize;

    std::env::set_var(Q845_SWAP_ONLY_T_PRIME_LENGTH_FLAG, "1");
    std::env::remove_var(Q845_LIFETIME_COEFFICIENT_FUSION_FLAG);
    assert!(!q845_swap_only_coefficient_dependencies_satisfied());
    std::env::set_var(Q845_LIFETIME_COEFFICIENT_FUSION_FLAG, "1");
    std::env::remove_var(PROMISED_LQ_SWAP_BORROW_FLAG);
    assert!(!q845_swap_only_swap_dependencies_satisfied());
    std::env::set_var(PROMISED_LQ_SWAP_BORROW_FLAG, "1");
    assert!(q845_swap_only_coefficient_dependencies_satisfied());
    assert!(q845_swap_only_swap_dependencies_satisfied());
    let mut dependency_checks = 4usize;
    std::env::set_var(Q851_TRUNCATED_SWAP_ONLY_GUARD_FLAG, "1");
    std::env::remove_var(Q845_SWAP_ONLY_T_PRIME_LENGTH_FLAG);
    assert!(!q845_swap_only_coefficient_dependencies_satisfied());
    std::env::set_var(Q845_SWAP_ONLY_T_PRIME_LENGTH_FLAG, "1");
    assert!(q845_swap_only_coefficient_dependencies_satisfied());
    dependency_checks += 2;

    std::env::set_var("LOWQ_DIRECT_PREFIX_BITLEN", "1");
    std::env::set_var("LOWQ_FUSED_ZERO_PREFIX_BITLEN", "1");
    assert_eq!(
        std::env::var("LOWQ_DIRECT_PREFIX_BITLEN").ok().as_deref(),
        Some("1")
    );
    assert_eq!(
        std::env::var("LOWQ_FUSED_ZERO_PREFIX_BITLEN")
            .ok()
            .as_deref(),
        Some("1")
    );
    dependency_checks += 2;
    std::env::set_var(PAIRED_BITLEN_SOURCE_COMPLEMENT_FLAG, "1");
    std::env::set_var(COEFFICIENT_NONNEGATIVE_X_CANCEL_FLAG, "1");
    std::env::set_var(PROMISED_SWAP_SUPPORT_LIFETIME_FUSION_FLAG, "1");
    let full_feature_environment_checks = [
        PAIRED_BITLEN_SOURCE_COMPLEMENT_FLAG,
        COEFFICIENT_NONNEGATIVE_X_CANCEL_FLAG,
        Q845_LIFETIME_COEFFICIENT_FUSION_FLAG,
        PROMISED_LQ_SWAP_BORROW_FLAG,
        PROMISED_SWAP_SUPPORT_LIFETIME_FUSION_FLAG,
        Q845_SWAP_ONLY_T_PRIME_LENGTH_FLAG,
        Q851_TRUNCATED_SWAP_ONLY_GUARD_FLAG,
    ]
    .into_iter()
    .filter(|flag| std::env::var(flag).ok().as_deref() == Some("1"))
    .count();
    assert_eq!(full_feature_environment_checks, 7);

    let mut truncated_range_comparator_cases_checked = 0usize;
    for width in 2usize..=4 {
        let modulus = 1usize << width;
        for value in 0..=modulus {
            let harness = build_q851_range_comparator_harness(width, value);
            assert!(harness.builder.next_qubit < 64);
            for register in 0..modulus {
                for target in [false, true] {
                    let mut input = q845_fusion_set_register(
                        0,
                        &harness.register_ids,
                        register,
                    );
                    input |= u64::from(target) << harness.target_id;
                    let output = apply_q845_fusion_classical_with_clean_resets(
                        &harness.builder.ops,
                        input,
                    );
                    let expected = target ^ (register >= value);
                    assert_eq!(((output >> harness.target_id) & 1) != 0, expected);
                    assert_eq!(
                        q845_fusion_read_register(output, &harness.register_ids),
                        register
                    );
                    assert_eq!(output & !harness.external_mask, 0);
                    truncated_range_comparator_cases_checked += 1;
                }
            }
        }
    }

    let trace_physical_width = 7usize;
    let trace_length_width = 3usize;
    let mut truncated_guard_layout_cases_checked = 0usize;
    let mut truncated_guard_trace_checks = 0usize;
    for trace_scan_width in 1..=trace_physical_width {
        let baseline = build_q851_truncated_guard_trace_harness(
            trace_physical_width,
            trace_scan_width,
            trace_length_width,
            false,
        );
        let candidate = build_q851_truncated_guard_trace_harness(
            trace_physical_width,
            trace_scan_width,
            trace_length_width,
            true,
        );
        assert_eq!(baseline.l_s_ids, candidate.l_s_ids);
        assert_eq!(baseline.l_r_prime_ids, candidate.l_r_prime_ids);
        assert_eq!(baseline.forward_trace_ids, candidate.forward_trace_ids);
        assert_eq!(baseline.reverse_trace_ids, candidate.reverse_trace_ids);
        assert_eq!(baseline.external_mask, candidate.external_mask);
        truncated_guard_layout_cases_checked += 1;
        for shift in 0..=trace_physical_width {
            for remainder_length in 0..=trace_physical_width - shift {
                let mut input = q845_fusion_set_register(0, &baseline.l_s_ids, shift);
                input = q845_fusion_set_register(
                    input,
                    &baseline.l_r_prime_ids,
                    remainder_length,
                );
                let baseline_output = apply_q845_fusion_classical_with_clean_resets(
                    &baseline.builder.ops,
                    input,
                );
                let candidate_output = apply_q845_fusion_classical_with_clean_resets(
                    &candidate.builder.ops,
                    input,
                );
                assert_eq!(candidate_output, baseline_output);
                assert_eq!(candidate_output & !candidate.external_mask, 0);
                let coefficient_width = trace_physical_width - shift - remainder_length;
                let active_width = coefficient_width.min(trace_scan_width);
                let expected_trace = (1usize << active_width) - 1;
                assert_eq!(
                    q845_fusion_read_register(candidate_output, &candidate.forward_trace_ids),
                    expected_trace
                );
                assert_eq!(
                    q845_fusion_read_register(candidate_output, &candidate.reverse_trace_ids),
                    expected_trace
                );
                truncated_guard_trace_checks += 1;
            }
        }
    }

    let mut production_truncated_address_cases_checked = 0usize;
    let production_modulus = 1usize << REFERENCE_LENGTH_WIDTH;
    for step in 1..=REFERENCE_STEPS {
        let scan = reference_active_windows(256, step).t_add_sub.1;
        let delta = REGISTER_SHARED_WORK_WIDTH - scan;
        let mut excluded_values = vec![0usize, delta, REGISTER_SHARED_WORK_WIDTH];
        if delta > 0 {
            excluded_values.push(delta - 1);
        }
        if delta < REGISTER_SHARED_WORK_WIDTH {
            excluded_values.push(delta + 1);
        }
        excluded_values.push(REGISTER_SHARED_WORK_WIDTH - 1);
        excluded_values.sort_unstable();
        excluded_values.dedup();
        for excluded in excluded_values {
            let coefficient_width = REGISTER_SHARED_WORK_WIDTH - excluded;
            let reverse_address =
                (excluded + production_modulus - delta) % production_modulus;
            let high = coefficient_width > scan;
            assert_eq!(reverse_address > scan, high);

            let mut forward_active = true;
            for index in 0..=scan {
                if coefficient_width == index {
                    forward_active = !forward_active;
                }
                if index < scan {
                    assert_eq!(forward_active, index < coefficient_width);
                }
            }
            forward_active ^= high;
            assert!(!forward_active);

            let mut reverse_active = high;
            for reverse_index in 0..=scan {
                let index = scan - reverse_index;
                if index < scan {
                    assert_eq!(reverse_active, index < coefficient_width);
                }
                if reverse_address == reverse_index {
                    reverse_active = !reverse_active;
                }
            }
            reverse_active = !reverse_active;
            assert!(!reverse_active);
            production_truncated_address_cases_checked += 1;
        }
    }

    let physical_work_width = 5usize;
    let scan_width = 4usize;
    let length_width = 3usize;
    let r_length_width = 2usize;
    let layouts = [(0usize, 0usize), (1, 0), (0, 1), (1, 1), (0, 2)];
    let mut layout_cases_checked = 0usize;
    let mut internal_boundary_layout_cases_checked = 0usize;
    let mut promised_basis_states_checked = 0usize;
    let mut oracle_transition_checks = 0usize;
    let mut inverse_pair_checks = 0usize;
    let mut scratch_clean_checks = 0usize;
    let mut cursor_restore_checks = 0usize;
    let mut count_restore_checks = 0usize;
    let mut residue_preservation_checks = 0usize;
    let mut excluded_suffix_preservation_checks = 0usize;

    for &(shift, r_length) in &layouts {
        let coefficient_width = physical_work_width - shift - r_length;
        assert!(coefficient_width > 0);
        std::env::remove_var(SUB800_INPLACE_GUARD_ADDRESS_FLAG);
        let baseline_forward = build_q845_swap_only_core_harness(
            physical_work_width,
            scan_width,
            coefficient_width.min(scan_width),
            length_width,
            r_length_width,
            false,
            false,
            false,
        );
        let baseline_inverse = build_q845_swap_only_core_harness(
            physical_work_width,
            scan_width,
            coefficient_width.min(scan_width),
            length_width,
            r_length_width,
            true,
            false,
            false,
        );
        std::env::set_var(SUB800_INPLACE_GUARD_ADDRESS_FLAG, "1");
        let candidate_forward = build_q845_swap_only_core_harness(
            physical_work_width,
            scan_width,
            scan_width,
            length_width,
            r_length_width,
            false,
            true,
            true,
        );
        let candidate_inverse = build_q845_swap_only_core_harness(
            physical_work_width,
            scan_width,
            scan_width,
            length_width,
            r_length_width,
            true,
            true,
            true,
        );
        assert!(candidate_forward.builder.next_qubit < 64);
        assert_eq!(baseline_forward.data_mask, candidate_forward.data_mask);
        assert_eq!(baseline_forward.work1_ids, candidate_forward.work1_ids);
        assert_eq!(baseline_forward.work2_ids, candidate_forward.work2_ids);
        let all_mask = (1u64 << candidate_forward.builder.next_qubit) - 1;
        let scratch_mask = all_mask & !candidate_forward.data_mask;
        let comparison_mask = candidate_forward.data_mask
            & !q845_fusion_id_mask(&candidate_forward.l_t_prime_ids)
            & !(1u64 << candidate_forward.above_guard_id);
        let excluded_suffix_mask = candidate_forward
            .work1_ids
            .iter()
            .skip(coefficient_width.min(scan_width))
            .chain(
                candidate_forward
                    .work2_ids
                    .iter()
                    .skip(coefficient_width.min(scan_width)),
            )
            .fold(0u64, |mask, id| mask | (1u64 << id));
        layout_cases_checked += 1;
        internal_boundary_layout_cases_checked += usize::from(coefficient_width < scan_width);
        for active_width in 1..=coefficient_width.min(scan_width) {
            let active_mask = (1usize << active_width) - 1;
            let source_mask = (1usize << active_width.saturating_sub(1)) - 1;
            for phase1 in [false, true] {
                for phase2 in [false, true] {
                    for sign in [false, true] {
                        let enable = phase1 && (phase2 || !sign);
                        let add_only = phase1 && !enable;
                        for work1 in 0..(1usize << scan_width) {
                            if work1 & (1usize << (active_width - 1)) != 0 {
                                continue;
                            }
                            let source = work1 & source_mask;
                            for work2 in 0..(1usize << physical_work_width) {
                                let target =
                                    work2 & ((1usize << coefficient_width.min(scan_width)) - 1);
                                let target_low = target & active_mask;
                                let above_guard = target >> active_width != 0;
                                if add_only
                                    && (above_guard
                                        || target_low
                                            >= (1usize << active_width.saturating_sub(1)))
                                {
                                    continue;
                                }
                                if add_only && target_low + source > active_mask {
                                    continue;
                                }

                                let mut baseline_input = 0u64;
                                baseline_input |= u64::from(phase1) << baseline_forward.phase1_id;
                                baseline_input |= u64::from(phase2) << baseline_forward.phase2_id;
                                baseline_input |= u64::from(sign) << baseline_forward.sign_id;
                                baseline_input = q845_fusion_set_register(
                                    baseline_input,
                                    &baseline_forward.work1_ids,
                                    work1,
                                );
                                baseline_input = q845_fusion_set_register(
                                    baseline_input,
                                    &baseline_forward.work2_ids,
                                    work2,
                                );
                                baseline_input = q845_fusion_set_register(
                                    baseline_input,
                                    &baseline_forward.l_t_ids,
                                    active_width - 1,
                                );
                                baseline_input = q845_fusion_set_register(
                                    baseline_input,
                                    &baseline_forward.l_s_ids,
                                    shift,
                                );
                                baseline_input = q845_fusion_set_register(
                                    baseline_input,
                                    &baseline_forward.l_r_prime_ids,
                                    r_length,
                                );
                                baseline_input |=
                                    u64::from(above_guard) << baseline_forward.above_guard_id;
                                let candidate_input =
                                    baseline_input & !(1u64 << baseline_forward.above_guard_id);

                                let baseline_output =
                                    apply_q845_fusion_classical_with_clean_resets(
                                        &baseline_forward.builder.ops,
                                        baseline_input,
                                    );
                                let candidate_output =
                                    apply_q845_fusion_classical_with_clean_resets(
                                        &candidate_forward.builder.ops,
                                        candidate_input,
                                    );
                                assert_eq!(
                                    candidate_output & comparison_mask,
                                    baseline_output & comparison_mask
                                );
                                assert_eq!(candidate_output & scratch_mask, 0);
                                assert_eq!(
                                    candidate_output
                                        & q845_fusion_id_mask(
                                            &candidate_forward.l_t_prime_ids,
                                        ),
                                    0
                                );
                                assert_eq!(
                                    candidate_output
                                        & (1u64 << candidate_forward.above_guard_id),
                                    0
                                );
                                assert_eq!(
                                    apply_q845_fusion_classical_with_clean_resets(
                                        &candidate_inverse.builder.ops,
                                        candidate_output,
                                    ),
                                    candidate_input
                                );
                                assert_eq!(
                                    apply_q845_fusion_classical_with_clean_resets(
                                        &baseline_inverse.builder.ops,
                                        baseline_output,
                                    ),
                                    baseline_input
                                );
                                assert_eq!(
                                    candidate_output & excluded_suffix_mask,
                                    candidate_input & excluded_suffix_mask
                                );
                                promised_basis_states_checked += 1;
                                oracle_transition_checks += 1;
                                inverse_pair_checks += 1;
                                scratch_clean_checks += 1;
                                cursor_restore_checks += 1;
                                count_restore_checks += 1;
                                residue_preservation_checks += 1;
                                excluded_suffix_preservation_checks += 1;
                            }
                        }
                    }
                }
            }
        }
    }

    std::env::remove_var(SUB800_INPLACE_GUARD_ADDRESS_FLAG);

    let swap_work_width = 3usize;
    let baseline_swap_forward = build_q845_ephemeral_swap_harness(
        swap_work_width,
        length_width,
        r_length_width,
        false,
        true,
        false,
    );
    let baseline_swap_inverse = build_q845_ephemeral_swap_harness(
        swap_work_width,
        length_width,
        r_length_width,
        true,
        true,
        false,
    );
    let candidate_swap_forward = build_q845_ephemeral_swap_harness(
        swap_work_width,
        length_width,
        r_length_width,
        false,
        true,
        true,
    );
    let candidate_swap_inverse = build_q845_ephemeral_swap_harness(
        swap_work_width,
        length_width,
        r_length_width,
        true,
        true,
        true,
    );
    let persistent_swap_forward = build_q845_ephemeral_swap_harness(
        swap_work_width,
        length_width,
        r_length_width,
        false,
        false,
        true,
    );
    let persistent_swap_inverse = build_q845_ephemeral_swap_harness(
        swap_work_width,
        length_width,
        r_length_width,
        true,
        false,
        true,
    );
    assert_eq!(
        baseline_swap_forward.data_ids,
        candidate_swap_forward.data_ids
    );
    assert_eq!(
        baseline_swap_forward.external_mask,
        candidate_swap_forward.external_mask
    );
    assert_eq!(
        persistent_swap_forward.data_ids,
        candidate_swap_forward.data_ids
    );
    assert_eq!(
        persistent_swap_forward.external_mask,
        candidate_swap_forward.external_mask
    );
    let lifecycle_comparison_mask = candidate_swap_forward.external_mask
        & !q845_fusion_id_mask(&candidate_swap_forward.l_t_prime_ids);
    let mut swap_inputs = Vec::new();
    let mut swap_outputs = Vec::new();
    let mut ephemeral_swap_cases_checked = 0usize;
    let mut ephemeral_control_on_checks = 0usize;
    let mut ephemeral_control_off_checks = 0usize;
    let mut persistent_lifecycle_equivalence_checks = 0usize;
    let mut ephemeral_inverse_pair_checks = 0usize;
    let mut ephemeral_l_t_prime_zero_checks = 0usize;
    let swap_mask = (1usize << swap_work_width) - 1;
    for t in 1..=swap_mask {
        for r in 0..=swap_mask {
            let l_t = q845_test_bit_length(t);
            let l_r = q845_test_bit_length(r);
            if l_t + 1 + l_r > swap_work_width {
                continue;
            }
            let work1 = q845_test_pack_work1(swap_work_width, t, r);
            for t_prime in 0..=swap_mask {
                for r_prime in 0..=swap_mask {
                    let l_t_prime = q845_test_bit_length(t_prime);
                    let l_r_prime = q845_test_bit_length(r_prime);
                    if l_t_prime + l_r_prime > swap_work_width {
                        continue;
                    }
                    let work2 =
                        q845_test_pack_work2(swap_work_width, t_prime, r_prime);
                    for iteration in [false, true] {
                        let mut input = u64::from(iteration)
                            << candidate_swap_forward.iteration_id;
                        input = q845_fusion_set_register(
                            input,
                            &candidate_swap_forward.work1_ids,
                            work1,
                        );
                        input = q845_fusion_set_register(
                            input,
                            &candidate_swap_forward.work2_ids,
                            work2,
                        );
                        input = q845_fusion_set_register(
                            input,
                            &candidate_swap_forward.l_t_ids,
                            l_t,
                        );
                        input = q845_fusion_set_register(
                            input,
                            &candidate_swap_forward.l_r_prime_ids,
                            l_r_prime,
                        );
                        let baseline_output = apply_scalar(
                            &baseline_swap_forward.builder.ops,
                            input,
                        );
                        let candidate_output = apply_scalar(
                            &candidate_swap_forward.builder.ops,
                            input,
                        );
                        let persistent_input = q845_fusion_set_register(
                            input,
                            &persistent_swap_forward.l_t_prime_ids,
                            l_t_prime,
                        );
                        let persistent_output = apply_scalar(
                            &persistent_swap_forward.builder.ops,
                            persistent_input,
                        );
                        assert_eq!(candidate_output, baseline_output);
                        assert_eq!(
                            candidate_output & lifecycle_comparison_mask,
                            persistent_output & lifecycle_comparison_mask
                        );
                        assert_eq!(
                            candidate_output & !candidate_swap_forward.external_mask,
                            0
                        );
                        assert_eq!(
                            q845_fusion_set_register(
                                candidate_output,
                                &candidate_swap_forward.work1_ids,
                                work2,
                            ) & q845_fusion_id_mask(&candidate_swap_forward.work1_ids),
                            candidate_output & q845_fusion_id_mask(&candidate_swap_forward.work1_ids)
                        );
                        assert_eq!(
                            q845_fusion_set_register(
                                candidate_output,
                                &candidate_swap_forward.work2_ids,
                                work1,
                            ) & q845_fusion_id_mask(&candidate_swap_forward.work2_ids),
                            candidate_output & q845_fusion_id_mask(&candidate_swap_forward.work2_ids)
                        );
                        assert_eq!(
                            q845_fusion_set_register(
                                candidate_output,
                                &candidate_swap_forward.l_t_ids,
                                l_t_prime,
                            ) & q845_fusion_id_mask(&candidate_swap_forward.l_t_ids),
                            candidate_output & q845_fusion_id_mask(&candidate_swap_forward.l_t_ids)
                        );
                        assert_eq!(
                            q845_fusion_set_register(
                                candidate_output,
                                &candidate_swap_forward.l_r_prime_ids,
                                l_r,
                            ) & q845_fusion_id_mask(&candidate_swap_forward.l_r_prime_ids),
                            candidate_output
                                & q845_fusion_id_mask(&candidate_swap_forward.l_r_prime_ids)
                        );
                        assert_eq!(
                            candidate_output
                                & q845_fusion_id_mask(
                                    &candidate_swap_forward.l_t_prime_ids,
                                ),
                            0
                        );
                        assert_eq!(
                            candidate_output
                                & q845_fusion_id_mask(&candidate_swap_forward.l_q_ids),
                            0
                        );
                        assert_eq!(
                            candidate_output
                                & q845_fusion_id_mask(&candidate_swap_forward.l_s_ids),
                            0
                        );
                        assert_eq!(
                            ((candidate_output >> candidate_swap_forward.iteration_id) & 1) != 0,
                            !iteration
                        );
                        assert_eq!(
                            apply_scalar(&candidate_swap_inverse.builder.ops, candidate_output),
                            input
                        );
                        assert_eq!(
                            apply_scalar(&baseline_swap_inverse.builder.ops, baseline_output),
                            input
                        );
                        assert_eq!(
                            apply_scalar(
                                &persistent_swap_inverse.builder.ops,
                                persistent_output,
                            ),
                            persistent_input
                        );
                        swap_inputs.push(q845_pack_data_value(
                            &candidate_swap_forward.data_ids,
                            input,
                        ));
                        swap_outputs.push(q845_pack_data_value(
                            &candidate_swap_forward.data_ids,
                            candidate_output,
                        ));
                        ephemeral_swap_cases_checked += 1;
                        ephemeral_control_on_checks += 1;
                        persistent_lifecycle_equivalence_checks += 1;
                        ephemeral_inverse_pair_checks += 1;
                        ephemeral_l_t_prime_zero_checks += 1;

                        for (blocked_l_q, blocked_l_s) in [(1usize, 0usize), (0, 1)] {
                            let mut blocked_input = q845_fusion_set_register(
                                input,
                                &candidate_swap_forward.l_q_ids,
                                blocked_l_q,
                            );
                            blocked_input = q845_fusion_set_register(
                                blocked_input,
                                &candidate_swap_forward.l_s_ids,
                                blocked_l_s,
                            );
                            let blocked_persistent_input = q845_fusion_set_register(
                                blocked_input,
                                &persistent_swap_forward.l_t_prime_ids,
                                l_t_prime,
                            );
                            let blocked_baseline_output = apply_scalar(
                                &baseline_swap_forward.builder.ops,
                                blocked_input,
                            );
                            let blocked_candidate_output = apply_scalar(
                                &candidate_swap_forward.builder.ops,
                                blocked_input,
                            );
                            let blocked_persistent_output = apply_scalar(
                                &persistent_swap_forward.builder.ops,
                                blocked_persistent_input,
                            );
                            assert_eq!(blocked_baseline_output, blocked_input);
                            assert_eq!(blocked_candidate_output, blocked_input);
                            assert_eq!(blocked_persistent_output, blocked_persistent_input);
                            assert_eq!(
                                blocked_candidate_output & lifecycle_comparison_mask,
                                blocked_persistent_output & lifecycle_comparison_mask
                            );
                            assert_eq!(
                                apply_scalar(
                                    &candidate_swap_inverse.builder.ops,
                                    blocked_candidate_output,
                                ),
                                blocked_input
                            );
                            assert_eq!(
                                apply_scalar(
                                    &baseline_swap_inverse.builder.ops,
                                    blocked_baseline_output,
                                ),
                                blocked_input
                            );
                            assert_eq!(
                                apply_scalar(
                                    &persistent_swap_inverse.builder.ops,
                                    blocked_persistent_output,
                                ),
                                blocked_persistent_input
                            );
                            swap_inputs.push(q845_pack_data_value(
                                &candidate_swap_forward.data_ids,
                                blocked_input,
                            ));
                            swap_outputs.push(q845_pack_data_value(
                                &candidate_swap_forward.data_ids,
                                blocked_candidate_output,
                            ));
                            ephemeral_swap_cases_checked += 1;
                            ephemeral_control_off_checks += 1;
                            persistent_lifecycle_equivalence_checks += 1;
                            ephemeral_inverse_pair_checks += 1;
                            ephemeral_l_t_prime_zero_checks += 1;
                        }
                    }
                }
            }
        }
    }
    let (_, forward_phase_checks, forward_ancilla_checks) =
        verify_q847_selected_simulator_equivalence(
            b"q845-ephemeral-swap-forward",
            &baseline_swap_forward.builder,
            &candidate_swap_forward.builder,
            &baseline_swap_forward.data_ids,
            baseline_swap_forward.external_mask,
            &swap_inputs,
        );
    let (_, inverse_phase_checks, inverse_ancilla_checks) =
        verify_q847_selected_simulator_equivalence(
            b"q845-ephemeral-swap-inverse",
            &baseline_swap_inverse.builder,
            &candidate_swap_inverse.builder,
            &baseline_swap_inverse.data_ids,
            baseline_swap_inverse.external_mask,
            &swap_outputs,
        );
    let ephemeral_phase_clean_checks = forward_phase_checks + inverse_phase_checks;
    let ephemeral_ancilla_clean_checks = forward_ancilla_checks + inverse_ancilla_checks;

    match saved_swap_only {
        Some(value) => std::env::set_var(Q845_SWAP_ONLY_T_PRIME_LENGTH_FLAG, value),
        None => std::env::remove_var(Q845_SWAP_ONLY_T_PRIME_LENGTH_FLAG),
    }
    match saved_truncated_guard {
        Some(value) => std::env::set_var(Q851_TRUNCATED_SWAP_ONLY_GUARD_FLAG, value),
        None => std::env::remove_var(Q851_TRUNCATED_SWAP_ONLY_GUARD_FLAG),
    }
    match saved_inplace_guard {
        Some(value) => std::env::set_var(SUB800_INPLACE_GUARD_ADDRESS_FLAG, value),
        None => std::env::remove_var(SUB800_INPLACE_GUARD_ADDRESS_FLAG),
    }
    match saved_fusion {
        Some(value) => std::env::set_var(Q845_LIFETIME_COEFFICIENT_FUSION_FLAG, value),
        None => std::env::remove_var(Q845_LIFETIME_COEFFICIENT_FUSION_FLAG),
    }
    match saved_promised {
        Some(value) => std::env::set_var(PROMISED_LQ_SWAP_BORROW_FLAG, value),
        None => std::env::remove_var(PROMISED_LQ_SWAP_BORROW_FLAG),
    }
    match saved_support_fusion {
        Some(value) => {
            std::env::set_var(PROMISED_SWAP_SUPPORT_LIFETIME_FUSION_FLAG, value)
        }
        None => std::env::remove_var(PROMISED_SWAP_SUPPORT_LIFETIME_FUSION_FLAG),
    }
    match saved_paired_source {
        Some(value) => std::env::set_var(PAIRED_BITLEN_SOURCE_COMPLEMENT_FLAG, value),
        None => std::env::remove_var(PAIRED_BITLEN_SOURCE_COMPLEMENT_FLAG),
    }
    match saved_coefficient_x {
        Some(value) => std::env::set_var(COEFFICIENT_NONNEGATIVE_X_CANCEL_FLAG, value),
        None => std::env::remove_var(COEFFICIENT_NONNEGATIVE_X_CANCEL_FLAG),
    }
    match saved_direct_prefix {
        Some(value) => std::env::set_var("LOWQ_DIRECT_PREFIX_BITLEN", value),
        None => std::env::remove_var("LOWQ_DIRECT_PREFIX_BITLEN"),
    }
    match saved_fused_zero_prefix {
        Some(value) => std::env::set_var("LOWQ_FUSED_ZERO_PREFIX_BITLEN", value),
        None => std::env::remove_var("LOWQ_FUSED_ZERO_PREFIX_BITLEN"),
    }

    Q845SwapOnlyCoefficientProofReport {
        dependency_checks,
        full_feature_environment_checks,
        truncated_guard_default_off_stream_identity_checks,
        truncated_range_comparator_cases_checked,
        truncated_guard_layout_cases_checked,
        truncated_guard_trace_checks,
        production_truncated_address_cases_checked,
        layout_cases_checked,
        internal_boundary_layout_cases_checked,
        promised_basis_states_checked,
        oracle_transition_checks,
        inverse_pair_checks,
        scratch_clean_checks,
        cursor_restore_checks,
        count_restore_checks,
        residue_preservation_checks,
        excluded_suffix_preservation_checks,
        default_off_stream_identity_checks: 1,
        ephemeral_swap_cases_checked,
        ephemeral_control_on_checks,
        ephemeral_control_off_checks,
        persistent_lifecycle_equivalence_checks,
        ephemeral_inverse_pair_checks,
        ephemeral_l_t_prime_zero_checks,
        ephemeral_phase_clean_checks,
        ephemeral_ancilla_clean_checks,
    }
}

/// Prove the carry-in guarded comparison and the zero-spill fused coefficient
/// core independently of the full point-add count.
#[doc(hidden)]
#[must_use]
pub fn exhaustive_q845_lifetime_coefficient_fusion_check(
) -> Q845LifetimeCoefficientFusionProofReport {
    const LENGTH_WIDTH: usize = 4;

    let saved = std::env::var_os(Q845_LIFETIME_COEFFICIENT_FUSION_FLAG);
    let default_stream = build_q845_fusion_dispatch_harness(None, false);
    let explicit_off_stream = build_q845_fusion_dispatch_harness(Some(false), false);
    assert_q845_fusion_builder_identity(&default_stream, &explicit_off_stream);
    let default_off_stream_identity_checks = 1;
    let candidate_stream = build_q845_fusion_dispatch_harness(Some(true), false);
    let direct_stream = build_q845_fusion_dispatch_harness(Some(false), true);
    assert_q845_fusion_builder_identity(&candidate_stream, &direct_stream);
    let dispatch_stream_identity_checks = 1;
    set_q845_lifetime_fusion_proof_mode(Some(true));

    let mut guard_basis_states_checked = 0usize;
    let mut guard_scratch_clean_checks = 0usize;
    let mut guard_carry_out_cases_checked = 0usize;
    for length_width in 1..=4 {
        let guard = build_q845_fusion_guard_harness(length_width);
        assert!(guard.builder.next_qubit < 64);
        let all_mask = (1u64 << guard.builder.next_qubit) - 1;
        let scratch_mask = all_mask & !guard.data_mask;
        let modulus = 1usize << length_width;
        for control in [false, true] {
            for target_bit in [false, true] {
                for target_length in 0..modulus {
                    for source_length in 0..modulus {
                        for shift in 0..modulus {
                            let mut input = 0u64;
                            input |= u64::from(control) << guard.control_id;
                            input |= u64::from(target_bit) << guard.target_id;
                            input = q845_fusion_set_register(
                                input,
                                &guard.target_length_ids,
                                target_length,
                            );
                            input = q845_fusion_set_register(
                                input,
                                &guard.source_length_ids,
                                source_length,
                            );
                            input = q845_fusion_set_register(input, &guard.shift_ids, shift);
                            let output = apply_q845_fusion_classical_with_clean_resets(
                                &guard.builder.ops,
                                input,
                            );
                            let toggle =
                                control && target_length > source_length + shift + 1;
                            let expected = if toggle {
                                input ^ (1u64 << guard.target_id)
                            } else {
                                input
                            };
                            assert_eq!(output & guard.data_mask, expected);
                            assert_eq!(output & scratch_mask, 0);
                            guard_basis_states_checked += 1;
                            guard_scratch_clean_checks += 1;
                            guard_carry_out_cases_checked +=
                                usize::from(source_length + shift + 1 >= modulus);
                        }
                    }
                }
            }
        }
    }

    let mut packed_basis_states_checked = 0usize;
    let mut packed_rotation_mapping_checks = 0usize;
    let mut packed_wrapped_residue_checks = 0usize;
    let mut packed_source_guard_clean_checks = 0usize;
    let mut packed_underflow_equivalence_checks = 0usize;
    let mut packed_add_headroom_checks = 0usize;
    for packed_width in 2..=8 {
        let packed_mask = (1usize << packed_width) - 1;
        for shift in 0..packed_width {
            for active_width in 1..=packed_width - shift {
                let active_mask = (1usize << active_width) - 1;
                let source_mask = (1usize << (active_width - 1)) - 1;
                for target in 0..=packed_mask {
                    let rotated = if shift == 0 {
                        target
                    } else {
                        ((target >> shift) | (target << (packed_width - shift))) & packed_mask
                    };
                    let quotient = target >> shift;
                    assert_eq!(rotated & active_mask, quotient & active_mask);
                    packed_rotation_mapping_checks += 1;
                    if shift > 0 {
                        assert_eq!(
                            rotated >> (packed_width - shift),
                            target & ((1usize << shift) - 1)
                        );
                        packed_wrapped_residue_checks += 1;
                    }
                    let above_guard = target >> (shift + active_width) != 0;
                    assert_eq!(above_guard, quotient > active_mask);
                    for source in 0..=source_mask {
                        let shifted_source = source << shift;
                        let rotated_source = if shift == 0 {
                            shifted_source
                        } else {
                            ((shifted_source >> shift)
                                | (shifted_source << (packed_width - shift)))
                                & packed_mask
                        };
                        assert_eq!(rotated_source, source);
                        assert_eq!(rotated_source & (1usize << (active_width - 1)), 0);
                        packed_source_guard_clean_checks += 1;
                        let fused_underflow =
                            !above_guard && (quotient & active_mask) < source;
                        assert_eq!(fused_underflow, target < shifted_source);
                        packed_underflow_equivalence_checks += 1;

                        if target < (1usize << (shift + active_width - 1)) {
                            let sum = target + shifted_source;
                            assert!(sum < (1usize << (shift + active_width)));
                            let rotated_sum = if shift == 0 {
                                sum
                            } else {
                                ((sum >> shift) | (sum << (packed_width - shift))) & packed_mask
                            };
                            assert_eq!(
                                rotated_sum & active_mask,
                                (quotient & active_mask) + source
                            );
                            if shift > 0 {
                                assert_eq!(
                                    rotated_sum >> (packed_width - shift),
                                    rotated >> (packed_width - shift)
                                );
                            }
                            packed_add_headroom_checks += 1;
                        }
                    }
                    packed_basis_states_checked += 1;
                }
            }
        }
    }

    let mut active_width_cases_checked = 0usize;
    let mut promised_basis_states_checked = 0usize;
    let mut oracle_transition_checks = 0usize;
    let mut inverse_pair_checks = 0usize;
    let mut control_off_identity_checks = 0usize;
    let mut underflow_checks = 0usize;
    let mut above_guard_checks = 0usize;
    let mut add_headroom_checks = 0usize;
    let mut scratch_clean_checks = 0usize;
    let mut cursor_restore_checks = 0usize;
    for work_width in 2..=5 {
        let forward = build_q845_fusion_core_harness(work_width, LENGTH_WIDTH, false);
        let inverse = build_q845_fusion_core_harness(work_width, LENGTH_WIDTH, true);
        assert_eq!(forward.data_mask, inverse.data_mask);
        assert!(forward.builder.next_qubit < 64);
        let all_mask = (1u64 << forward.builder.next_qubit) - 1;
        let scratch_mask = all_mask & !forward.data_mask;
        for active_width in 1..=work_width {
            active_width_cases_checked += 1;
            let active_mask = (1usize << active_width) - 1;
            let source_mask = (1usize << active_width.saturating_sub(1)) - 1;
            for phase1 in [false, true] {
                for phase2 in [false, true] {
                    for sign in [false, true] {
                        let enable = phase1 && (phase2 || !sign);
                        let add_only = phase1 && !enable;
                        for work1 in 0..(1usize << work_width) {
                            if work1 & (1usize << (active_width - 1)) != 0 {
                                continue;
                            }
                            let source = work1 & source_mask;
                            for target in 0..(1usize << work_width) {
                                let target_low = target & active_mask;
                                for above_guard in [false, true] {
                                    if add_only
                                        && (above_guard
                                            || target_low
                                                >= (1usize
                                                    << active_width.saturating_sub(1)))
                                    {
                                        continue;
                                    }
                                    let sum = target_low + source;
                                    if add_only && sum > active_mask {
                                        continue;
                                    }
                                    let underflow =
                                        enable && !above_guard && target_low < source;
                                    let mut input = 0u64;
                                    input |= u64::from(phase1) << forward.phase1_id;
                                    input |= u64::from(phase2) << forward.phase2_id;
                                    input |= u64::from(sign) << forward.sign_id;
                                    input = q845_fusion_set_register(
                                        input,
                                        &forward.work1_ids,
                                        work1,
                                    );
                                    input = q845_fusion_set_register(
                                        input,
                                        &forward.work2_ids,
                                        target,
                                    );
                                    input = q845_fusion_set_register(
                                        input,
                                        &forward.length_ids,
                                        active_width - 1,
                                    );
                                    input |=
                                        u64::from(above_guard) << forward.above_guard_id;
                                    let output = apply_q845_fusion_classical_with_clean_resets(
                                        &forward.builder.ops,
                                        input,
                                    );
                                    assert_eq!(output & scratch_mask, 0);
                                    scratch_clean_checks += 1;
                                    let mut expected = input;
                                    if add_only {
                                        let expected_target = (target & !active_mask) | sum;
                                        expected = q845_fusion_set_register(
                                            expected,
                                            &forward.work2_ids,
                                            expected_target,
                                        );
                                        add_headroom_checks += 1;
                                    }
                                    if phase1 ^ underflow {
                                        expected ^= 1u64 << forward.sign_id;
                                    }
                                    assert_eq!(output & forward.data_mask, expected);
                                    oracle_transition_checks += 1;
                                    assert_eq!(
                                        apply_q845_fusion_classical_with_clean_resets(
                                            &inverse.builder.ops,
                                            output,
                                        ),
                                        input
                                    );
                                    inverse_pair_checks += 1;
                                    if !phase1 {
                                        assert_eq!(output, input);
                                        control_off_identity_checks += 1;
                                    }
                                    underflow_checks += usize::from(underflow);
                                    above_guard_checks += usize::from(enable && above_guard);
                                    cursor_restore_checks += 1;
                                    promised_basis_states_checked += 1;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    match saved {
        Some(value) => std::env::set_var(Q845_LIFETIME_COEFFICIENT_FUSION_FLAG, value),
        None => std::env::remove_var(Q845_LIFETIME_COEFFICIENT_FUSION_FLAG),
    }
    Q845LifetimeCoefficientFusionProofReport {
        guard_widths_checked: 4,
        guard_basis_states_checked,
        guard_scratch_clean_checks,
        guard_carry_out_cases_checked,
        packed_widths_checked: 7,
        packed_basis_states_checked,
        packed_rotation_mapping_checks,
        packed_wrapped_residue_checks,
        packed_source_guard_clean_checks,
        packed_underflow_equivalence_checks,
        packed_add_headroom_checks,
        minimum_spill_width: 0,
        work_widths_checked: 4,
        active_width_cases_checked,
        promised_basis_states_checked,
        oracle_transition_checks,
        inverse_pair_checks,
        control_off_identity_checks,
        underflow_checks,
        above_guard_checks,
        add_headroom_checks,
        scratch_clean_checks,
        cursor_restore_checks,
        default_off_stream_identity_checks,
        dispatch_stream_identity_checks,
    }
}

struct FusedPrefixScratchLoanHarness {
    builder: B,
    data_ids: Vec<u32>,
    external_mask: u64,
    source_mask: u64,
    lender_mask: u64,
    trace: FusedPrefixScratchLoanAllocationTrace,
}

fn build_fused_prefix_scratch_loan_harness(
    source_width: usize,
    accumulator_width: usize,
    lender_count: usize,
    decrement: bool,
    direct_unloaned_entry: bool,
) -> FusedPrefixScratchLoanHarness {
    let mut circ = Circuit::new();
    let source = circ.alloc_qreg_bits("rs.fused-prefix-loan-proof.source", source_width);
    let accumulator =
        circ.alloc_qreg_bits("rs.fused-prefix-loan-proof.accumulator", accumulator_width);
    let lenders = circ.alloc_qreg_bits("rs.fused-prefix-loan-proof.lender", lender_count);
    let source_refs: Vec<&QReg> = source.iter().collect();
    let lender_refs: Vec<&QReg> = lenders.iter().collect();
    let data_ids: Vec<u32> = source.iter().chain(&accumulator).map(QReg::id).collect();
    let external_mask = qreg_mask(source.iter().chain(&accumulator).chain(&lenders));
    let source_mask = qreg_mask(&source);
    let lender_mask = qreg_mask(&lenders);

    begin_fused_prefix_scratch_loan_allocation_trace();
    if direct_unloaned_entry {
        bit_length_lean_allow_zero(&mut circ, &source_refs, &accumulator, decrement);
    } else {
        bit_length_lean_allow_zero_with_borrowed_scratch(
            &mut circ,
            &source_refs,
            &accumulator,
            decrement,
            None,
            &lender_refs,
        );
    }
    let trace = finish_fused_prefix_scratch_loan_allocation_trace();
    FusedPrefixScratchLoanHarness {
        builder: circ.into_builder(),
        data_ids,
        external_mask,
        source_mask,
        lender_mask,
        trace,
    }
}

fn fused_prefix_scratch_loan_input(harness: &FusedPrefixScratchLoanHarness, value: u64) -> u64 {
    harness
        .data_ids
        .iter()
        .enumerate()
        .fold(0u64, |state, (bit, id)| {
            state | (((value >> bit) & 1) << id)
        })
}

fn assert_fused_prefix_scratch_loan_clean(
    harness: &FusedPrefixScratchLoanHarness,
    state: u64,
    context: &str,
) {
    assert_eq!(state & harness.lender_mask, 0, "{context}: lender dirty");
    assert_eq!(
        state & !harness.external_mask,
        0,
        "{context}: internal ancilla dirty"
    );
}

fn assert_fused_prefix_scratch_loan_stream_identity(
    left: &FusedPrefixScratchLoanHarness,
    right: &FusedPrefixScratchLoanHarness,
) {
    assert_eq!(left.builder.ops, right.builder.ops);
    assert_eq!(left.builder.next_qubit, right.builder.next_qubit);
    assert_eq!(left.builder.next_bit, right.builder.next_bit);
    assert_eq!(left.builder.active_qubits, right.builder.active_qubits);
    assert_eq!(left.builder.peak_qubits, right.builder.peak_qubits);
    assert_eq!(left.builder.free_qubits, right.builder.free_qubits);
    assert_eq!(
        left.builder.allocation_serial,
        right.builder.allocation_serial
    );
}

fn verify_fused_prefix_scratch_loan_simulator_equivalence(
    baseline: &FusedPrefixScratchLoanHarness,
    candidate: &FusedPrefixScratchLoanHarness,
) -> (usize, usize) {
    use crate::circuit::QubitId;
    use crate::sim::Simulator;
    use sha3::{
        digest::{ExtendableOutput, Update},
        Shake128,
    };

    assert_eq!(baseline.data_ids, candidate.data_ids);
    assert_eq!(baseline.external_mask, candidate.external_mask);
    let states = 1usize << baseline.data_ids.len();
    let mut cases_checked = 0usize;
    let mut phase_clean_checks = 0usize;
    for batch_start in (0..states).step_by(64) {
        let shots = (states - batch_start).min(64);
        let live = if shots == 64 {
            u64::MAX
        } else {
            (1u64 << shots) - 1
        };
        let mut baseline_seed = Shake128::default();
        baseline_seed.update(b"fused-prefix-scratch-loan-proof");
        baseline_seed.update(&(batch_start as u64).to_le_bytes());
        let mut baseline_xof = baseline_seed.finalize_xof();
        let mut baseline_simulator = Simulator::new(
            baseline.builder.next_qubit as usize,
            baseline.builder.next_bit as usize,
            &mut baseline_xof,
        );
        let mut candidate_seed = Shake128::default();
        candidate_seed.update(b"fused-prefix-scratch-loan-proof");
        candidate_seed.update(&(batch_start as u64).to_le_bytes());
        let mut candidate_xof = candidate_seed.finalize_xof();
        let mut candidate_simulator = Simulator::new(
            candidate.builder.next_qubit as usize,
            candidate.builder.next_bit as usize,
            &mut candidate_xof,
        );

        for shot in 0..shots {
            let value = (batch_start + shot) as u64;
            for (bit, (&baseline_id, &candidate_id)) in baseline
                .data_ids
                .iter()
                .zip(&candidate.data_ids)
                .enumerate()
            {
                if (value >> bit) & 1 != 0 {
                    *baseline_simulator.qubit_mut(QubitId(u64::from(baseline_id))) |= 1u64 << shot;
                    *candidate_simulator.qubit_mut(QubitId(u64::from(candidate_id))) |=
                        1u64 << shot;
                }
            }
        }
        baseline_simulator.apply_iter(baseline.builder.ops.iter());
        candidate_simulator.apply_iter(candidate.builder.ops.iter());
        assert_eq!(baseline_simulator.phase & live, 0);
        assert_eq!(candidate_simulator.phase & live, 0);
        phase_clean_checks += 2 * shots;

        for id in 0..baseline.builder.next_qubit {
            let baseline_value = baseline_simulator.qubit(QubitId(u64::from(id))) & live;
            if baseline.external_mask & (1u64 << id) != 0 {
                assert_eq!(
                    baseline_value,
                    candidate_simulator.qubit(QubitId(u64::from(id))) & live
                );
            } else {
                assert_eq!(baseline_value, 0, "baseline simulator left q{id} dirty");
            }
        }
        for id in 0..candidate.builder.next_qubit {
            if candidate.external_mask & (1u64 << id) == 0 {
                assert_eq!(
                    candidate_simulator.qubit(QubitId(u64::from(id))) & live,
                    0,
                    "candidate simulator left q{id} dirty"
                );
            }
        }
        cases_checked += shots;
    }
    (cases_checked, phase_clean_checks)
}

fn build_fused_prefix_scratch_loan_local_builder(
    source_width: usize,
    accumulator_width: usize,
    lender_count: usize,
    decrement: bool,
) -> (B, FusedPrefixScratchLoanAllocationTrace) {
    let mut circ = Circuit::new();
    let source = circ.alloc_qreg_bits("rs.fused-prefix-loan-proof.source", source_width);
    let accumulator =
        circ.alloc_qreg_bits("rs.fused-prefix-loan-proof.accumulator", accumulator_width);
    let lenders = circ.alloc_qreg_bits("rs.fused-prefix-loan-proof.lender", lender_count);
    let source_refs: Vec<&QReg> = source.iter().collect();
    let lender_refs: Vec<&QReg> = lenders.iter().collect();

    begin_fused_prefix_scratch_loan_allocation_trace();
    bit_length_lean_allow_zero_with_borrowed_scratch(
        &mut circ,
        &source_refs,
        &accumulator,
        decrement,
        None,
        &lender_refs,
    );
    let trace = finish_fused_prefix_scratch_loan_allocation_trace();
    (circ.into_builder(), trace)
}

fn fused_prefix_scratch_loan_local_resources(builder: &B) -> FusedPrefixScratchLoanLocalResources {
    let counts = measurement_classical_gate_counts(&builder.ops);
    let active_qubits = builder.active_qubits as usize;
    let peak_qubits = builder.peak_qubits as usize;
    FusedPrefixScratchLoanLocalResources {
        active_qubits,
        peak_qubits,
        temporary_qubits: peak_qubits - active_qubits,
        emitted_ops: builder.ops.len(),
        emitted_toffoli: counts.ccx,
    }
}

fn configure_fused_prefix_scratch_loan_proof(enabled: bool) {
    configure_raw_bit_length_loan_proof(RawBitLengthZeroMode::Fused, false);
    if enabled {
        std::env::set_var(FUSED_PREFIX_SCRATCH_LOAN_FLAG, "1");
    } else {
        std::env::remove_var(FUSED_PREFIX_SCRATCH_LOAN_FLAG);
    }
}

fn configure_fused_prefix_scratch_loan_kg_reverse_proof() {
    configure_fused_prefix_scratch_loan_proof(true);
    std::env::set_var(
        super::shrunken_pz_state_machine::CALLER_SCRATCH_KG_REVERSE_DECREMENT_FLAG,
        "1",
    );
}

fn fused_prefix_scratch_loan_rejections() -> usize {
    use std::panic::{catch_unwind, AssertUnwindSafe};

    fn build_rejection(kind: usize) {
        configure_fused_prefix_scratch_loan_proof(true);
        let mut circ = Circuit::new();
        let source = circ.alloc_qreg_bits("rs.fused-prefix-loan-reject.source", 2);
        let output = circ.alloc_qreg_bits("rs.fused-prefix-loan-reject.output", 3);
        let lenders = circ.alloc_qreg_bits("rs.fused-prefix-loan-reject.lender", 10);
        let source_refs: Vec<&QReg> = source.iter().collect();
        let lender_refs: Vec<&QReg> = match kind {
            0 => vec![&lenders[0], &lenders[0]],
            1 => vec![&source[0]],
            2 => vec![&output[0]],
            3 => lenders.iter().collect(),
            4 | 5 => vec![&lenders[0]],
            _ => unreachable!(),
        };
        if kind == 4 {
            std::env::remove_var("LOWQ_FUSED_ZERO_PREFIX_BITLEN");
        }
        if kind == 5 {
            std::env::remove_var("LOWQ_REUSE_ZERO_CARRIES_FOR_FULL_PREFIX_SCRATCH");
        }
        bit_length_lean_allow_zero_with_borrowed_scratch(
            &mut circ,
            &source_refs,
            &output,
            false,
            None,
            &lender_refs,
        );
    }

    let previous_hook = std::panic::take_hook();
    std::panic::set_hook(Box::new(|_| {}));
    let mut rejected = 0usize;
    for kind in 0..6 {
        let result = catch_unwind(AssertUnwindSafe(|| build_rejection(kind)));
        assert!(
            result.is_err(),
            "fused-prefix loan rejection {kind} unexpectedly passed"
        );
        rejected += 1;
    }
    std::panic::set_hook(previous_hook);
    configure_fused_prefix_scratch_loan_proof(false);
    rejected
}

/// Exhaustively verify partial borrowing for the nine-lane fused-prefix
/// layout. The three-, seven-, eight-, and nine-lender shapes are exact
/// production layouts in the coefficient and sub-800 swap routes.
#[doc(hidden)]
pub fn exhaustive_fused_prefix_scratch_loan_check() -> FusedPrefixScratchLoanProofReport {
    assert!(
        std::env::var_os(FUSED_PREFIX_SCRATCH_LOAN_FLAG).is_none(),
        "the fused-prefix scratch loan must default off"
    );
    assert!(
        std::env::var_os(
            super::shrunken_pz_state_machine::CALLER_SCRATCH_KG_REVERSE_DECREMENT_FLAG,
        )
        .is_none(),
        "the caller-scratch KG reverse decrement must default off"
    );
    let _environment = RawBitLengthProofEnvironment::capture();
    const MAX_SOURCE_WIDTH: usize = 8;
    const ACCUMULATOR_WIDTH: usize = 5;
    let lender_modes = [3usize, 7usize, 8usize, 9usize];
    let mut basis_states_checked = 0usize;
    let mut scalar_equivalence_checks = 0usize;
    let mut simulator_equivalence_checks = 0usize;
    let mut phase_clean_checks = 0usize;
    let mut inverse_pair_checks = 0usize;
    let mut source_restore_checks = 0usize;
    let mut lender_clean_checks = 0usize;
    let mut ancilla_clean_checks = 0usize;
    let mut kg_reverse_composition_checks = 0usize;
    let mut kg_reverse_simulator_checks = 0usize;
    let mut kg_reverse_phase_clean_checks = 0usize;
    let mut kg_reverse_inverse_pair_checks = 0usize;
    let mut kg_reverse_scratch_clean_checks = 0usize;
    let mut kg_reverse_changed_streams = 0usize;
    let mut three_lender_trace = None;
    let mut seven_lender_trace = None;
    let mut eight_lender_trace = None;
    let mut nine_lender_trace = None;

    for source_width in 0..=MAX_SOURCE_WIDTH {
        for lender_count in lender_modes {
            configure_fused_prefix_scratch_loan_proof(false);
            let baseline_add = build_fused_prefix_scratch_loan_harness(
                source_width,
                ACCUMULATOR_WIDTH,
                lender_count,
                false,
                false,
            );
            let baseline_sub = build_fused_prefix_scratch_loan_harness(
                source_width,
                ACCUMULATOR_WIDTH,
                lender_count,
                true,
                false,
            );
            configure_fused_prefix_scratch_loan_proof(true);
            let candidate_add = build_fused_prefix_scratch_loan_harness(
                source_width,
                ACCUMULATOR_WIDTH,
                lender_count,
                false,
                false,
            );
            let candidate_sub = build_fused_prefix_scratch_loan_harness(
                source_width,
                ACCUMULATOR_WIDTH,
                lender_count,
                true,
                false,
            );
            configure_fused_prefix_scratch_loan_kg_reverse_proof();
            let composed_add = build_fused_prefix_scratch_loan_harness(
                source_width,
                ACCUMULATOR_WIDTH,
                lender_count,
                false,
                false,
            );
            let composed_sub = build_fused_prefix_scratch_loan_harness(
                source_width,
                ACCUMULATOR_WIDTH,
                lender_count,
                true,
                false,
            );
            if source_width != 0 {
                for baseline in [&baseline_add, &baseline_sub] {
                    assert_eq!(baseline.trace.calls, 1);
                    assert_eq!(baseline.trace.owned_lanes, 9);
                    assert_eq!(baseline.trace.borrowed_lanes, 0);
                }
                for candidate in [&candidate_add, &candidate_sub] {
                    assert_eq!(candidate.trace.calls, 1);
                    assert_eq!(candidate.trace.owned_lanes, 9 - lender_count);
                    assert_eq!(candidate.trace.borrowed_lanes, lender_count);
                }
                assert_eq!(composed_add.trace, candidate_add.trace);
                assert_eq!(composed_sub.trace, candidate_sub.trace);
                if source_width == MAX_SOURCE_WIDTH {
                    match lender_count {
                        3 => three_lender_trace = Some(candidate_add.trace),
                        7 => seven_lender_trace = Some(candidate_add.trace),
                        8 => eight_lender_trace = Some(candidate_add.trace),
                        9 => nine_lender_trace = Some(candidate_add.trace),
                        _ => unreachable!(),
                    }
                }
            }

            for (baseline, candidate) in [
                (&baseline_add, &candidate_add),
                (&baseline_sub, &candidate_sub),
            ] {
                let (cases, phases) =
                    verify_fused_prefix_scratch_loan_simulator_equivalence(baseline, candidate);
                simulator_equivalence_checks += cases;
                phase_clean_checks += phases;
            }
            for (candidate, composed) in [
                (&candidate_add, &composed_add),
                (&candidate_sub, &composed_sub),
            ] {
                if candidate.builder.ops != composed.builder.ops {
                    kg_reverse_changed_streams += 1;
                }
                let (cases, phases) =
                    verify_fused_prefix_scratch_loan_simulator_equivalence(candidate, composed);
                kg_reverse_simulator_checks += cases;
                kg_reverse_phase_clean_checks += phases;
            }

            let data_states = 1u64 << baseline_add.data_ids.len();
            for value in 0..data_states {
                let input = fused_prefix_scratch_loan_input(&baseline_add, value);
                let baseline_added = apply_scalar(&baseline_add.builder.ops, input);
                let candidate_added = apply_scalar(&candidate_add.builder.ops, input);
                let baseline_subtracted = apply_scalar(&baseline_sub.builder.ops, input);
                let candidate_subtracted = apply_scalar(&candidate_sub.builder.ops, input);
                let composed_added = apply_scalar(&composed_add.builder.ops, input);
                let composed_subtracted = apply_scalar(&composed_sub.builder.ops, input);
                assert_eq!(candidate_added, baseline_added);
                assert_eq!(candidate_subtracted, baseline_subtracted);
                assert_eq!(composed_added, candidate_added);
                assert_eq!(composed_subtracted, candidate_subtracted);
                for (harness, state, context) in [
                    (&baseline_add, baseline_added, "baseline add"),
                    (&candidate_add, candidate_added, "candidate add"),
                    (&baseline_sub, baseline_subtracted, "baseline subtract"),
                    (&candidate_sub, candidate_subtracted, "candidate subtract"),
                ] {
                    assert_fused_prefix_scratch_loan_clean(harness, state, context);
                }
                for (harness, state, context) in [
                    (&composed_add, composed_added, "KG reverse composed add"),
                    (
                        &composed_sub,
                        composed_subtracted,
                        "KG reverse composed subtract",
                    ),
                ] {
                    assert_fused_prefix_scratch_loan_clean(harness, state, context);
                }
                assert_eq!(
                    candidate_added & candidate_add.source_mask,
                    input & candidate_add.source_mask
                );
                assert_eq!(
                    candidate_subtracted & candidate_sub.source_mask,
                    input & candidate_sub.source_mask
                );
                assert_eq!(
                    apply_scalar(&candidate_sub.builder.ops, candidate_added),
                    input
                );
                assert_eq!(
                    apply_scalar(&candidate_add.builder.ops, candidate_subtracted),
                    input
                );
                assert_eq!(
                    apply_scalar(&composed_sub.builder.ops, composed_added),
                    input
                );
                assert_eq!(
                    apply_scalar(&composed_add.builder.ops, composed_subtracted),
                    input
                );
                basis_states_checked += 2;
                scalar_equivalence_checks += 2;
                inverse_pair_checks += 2;
                source_restore_checks += 2;
                lender_clean_checks += 2;
                ancilla_clean_checks += 4;
                kg_reverse_composition_checks += 2;
                kg_reverse_inverse_pair_checks += 2;
                kg_reverse_scratch_clean_checks += 2;
            }
        }
    }

    let mut default_stream_identity_checks = 0usize;
    configure_fused_prefix_scratch_loan_proof(false);
    for source_width in 0..=MAX_SOURCE_WIDTH {
        for lender_count in lender_modes {
            for decrement in [false, true] {
                let configured = build_fused_prefix_scratch_loan_harness(
                    source_width,
                    ACCUMULATOR_WIDTH,
                    lender_count,
                    decrement,
                    false,
                );
                let direct = build_fused_prefix_scratch_loan_harness(
                    source_width,
                    ACCUMULATOR_WIDTH,
                    lender_count,
                    decrement,
                    true,
                );
                assert_fused_prefix_scratch_loan_stream_identity(&configured, &direct);
                default_stream_identity_checks += 1;
            }
        }
    }

    configure_fused_prefix_scratch_loan_proof(false);
    let (baseline_local_builder, baseline_local_trace) =
        build_fused_prefix_scratch_loan_local_builder(259, 10, 9, false);
    configure_fused_prefix_scratch_loan_proof(true);
    let (candidate_local_builder, candidate_local_trace) =
        build_fused_prefix_scratch_loan_local_builder(259, 10, 9, false);
    let baseline_local = fused_prefix_scratch_loan_local_resources(&baseline_local_builder);
    let candidate_local = fused_prefix_scratch_loan_local_resources(&candidate_local_builder);
    assert_eq!(baseline_local.active_qubits, candidate_local.active_qubits);
    assert_eq!(candidate_local.peak_qubits + 9, baseline_local.peak_qubits);
    assert_eq!(
        candidate_local.emitted_toffoli,
        baseline_local.emitted_toffoli
    );
    assert_eq!(candidate_local.emitted_ops + 9, baseline_local.emitted_ops);

    let alias_rejections = fused_prefix_scratch_loan_rejections();
    let three_lender_trace = three_lender_trace.expect("three-lender allocation trace");
    let seven_lender_trace = seven_lender_trace.expect("seven-lender allocation trace");
    let eight_lender_trace = eight_lender_trace.expect("eight-lender allocation trace");
    let nine_lender_trace = nine_lender_trace.expect("nine-lender allocation trace");
    FusedPrefixScratchLoanProofReport {
        source_widths_checked: MAX_SOURCE_WIDTH + 1,
        maximum_source_width: MAX_SOURCE_WIDTH,
        accumulator_width: ACCUMULATOR_WIDTH,
        lender_modes_checked: lender_modes.len(),
        directions_checked: 2,
        basis_states_checked,
        scalar_equivalence_checks,
        simulator_equivalence_checks,
        phase_clean_checks,
        inverse_pair_checks,
        source_restore_checks,
        lender_clean_checks,
        ancilla_clean_checks,
        default_stream_identity_checks,
        kg_reverse_composition_checks,
        kg_reverse_simulator_checks,
        kg_reverse_phase_clean_checks,
        kg_reverse_inverse_pair_checks,
        kg_reverse_scratch_clean_checks,
        kg_reverse_changed_streams,
        alias_rejections,
        three_lender_owned_lanes: three_lender_trace.owned_lanes,
        three_lender_borrowed_lanes: three_lender_trace.borrowed_lanes,
        seven_lender_owned_lanes: seven_lender_trace.owned_lanes,
        seven_lender_borrowed_lanes: seven_lender_trace.borrowed_lanes,
        eight_lender_owned_lanes: eight_lender_trace.owned_lanes,
        eight_lender_borrowed_lanes: eight_lender_trace.borrowed_lanes,
        nine_lender_owned_lanes: nine_lender_trace.owned_lanes,
        nine_lender_borrowed_lanes: nine_lender_trace.borrowed_lanes,
        baseline_local_trace,
        candidate_local_trace,
        baseline_local,
        candidate_local,
        local_qubit_delta: candidate_local.peak_qubits as i64 - baseline_local.peak_qubits as i64,
        local_ops_delta: candidate_local.emitted_ops as i64 - baseline_local.emitted_ops as i64,
        local_toffoli_delta: candidate_local.emitted_toffoli as i64
            - baseline_local.emitted_toffoli as i64,
    }
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

fn build_remainder_kernel_for_layout(
    work_width: usize,
    total_work_width: usize,
    window_upper: usize,
    length_width: usize,
    kernel: RemainderKernel,
    layout: RemainderScratchLayout,
) -> B {
    assert!(
        std::env::var_os(PHASE_OVERLAID_REMAINDER_SCRATCH_FLAG).is_none(),
        "remainder layout proof requires the feature flag to be initially absent"
    );
    if layout == RemainderScratchLayout::PhaseOverlaid {
        std::env::set_var(PHASE_OVERLAID_REMAINDER_SCRATCH_FLAG, "1");
    }
    let builder = build_remainder_kernel(
        work_width,
        total_work_width,
        window_upper,
        length_width,
        kernel,
    );
    std::env::remove_var(PHASE_OVERLAID_REMAINDER_SCRATCH_FLAG);
    builder
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
struct RemainderPhaseLengths {
    compute_operation: usize,
    toggle_enable: usize,
    prepare_range: usize,
    unprepare_range: usize,
    uncompute_operation: usize,
}

fn remainder_phase_lengths(
    total_work_width: usize,
    window_upper: usize,
    length_width: usize,
    add: bool,
) -> RemainderPhaseLengths {
    let mut circ = Circuit::new();
    let phase1 = circ.alloc_qreg("rs.remainder-phase-lengths.phase1");
    let phase2 = circ.alloc_qreg("rs.remainder-phase-lengths.phase2");
    let sign = circ.alloc_qreg("rs.remainder-phase-lengths.sign");
    let l_t = circ.alloc_qreg_bits("rs.remainder-phase-lengths.l-t", length_width);
    let l_q = circ.alloc_qreg_bits("rs.remainder-phase-lengths.l-q", length_width);
    let l_s = circ.alloc_qreg_bits("rs.remainder-phase-lengths.l-s", length_width);
    let l_r_prime = circ.alloc_qreg_bits("rs.remainder-phase-lengths.l-r-prime", length_width);
    let raw_scratch = circ.alloc_qreg_bits(
        "rs.remainder-phase-lengths.scratch",
        phase_overlaid_remainder_scratch_width(length_width, length_width),
    );
    let scratch = split_phase_overlaid_remainder_scratch(&raw_scratch, length_width, length_width);

    let start = circ.total_ops() as usize;
    compute_remainder_operation(&mut circ, &phase1, &l_r_prime, &scratch);
    let after_compute = circ.total_ops() as usize;
    if add {
        toggle_remainder_add_enable(
            &mut circ,
            scratch.operation,
            &phase2,
            &sign,
            scratch.phase_sign,
            scratch.enable,
        );
    }
    let after_enable = circ.total_ops() as usize;
    prepare_remainder_range(
        &mut circ,
        total_work_width,
        window_upper,
        &l_t,
        &l_q,
        &l_s,
        &scratch,
    );
    let after_prepare = circ.total_ops() as usize;
    unprepare_remainder_range(
        &mut circ,
        total_work_width,
        window_upper,
        &l_t,
        &l_q,
        &l_s,
        &scratch,
    );
    let after_unprepare = circ.total_ops() as usize;
    if add {
        toggle_remainder_add_enable(
            &mut circ,
            scratch.operation,
            &phase2,
            &sign,
            scratch.phase_sign,
            scratch.enable,
        );
    }
    let after_disable = circ.total_ops() as usize;
    uncompute_remainder_operation(&mut circ, &phase1, &l_r_prime, &scratch);
    let after_uncompute = circ.total_ops() as usize;

    RemainderPhaseLengths {
        compute_operation: after_compute - start,
        toggle_enable: after_enable - after_compute,
        prepare_range: after_prepare - after_enable,
        unprepare_range: after_unprepare - after_prepare,
        uncompute_operation: after_uncompute - after_disable,
    }
}

fn remainder_phase_boundaries(
    emitted_ops: usize,
    total_work_width: usize,
    window_upper: usize,
    length_width: usize,
    kernel: RemainderKernel,
) -> Vec<(usize, u64)> {
    let add = matches!(kernel, RemainderKernel::Add | RemainderKernel::AddInverse);
    let lengths = remainder_phase_lengths(total_work_width, window_upper, length_width, add);
    let predicate_mask = 0b11u64;
    let enabled_mask = predicate_mask | (1 << 2);
    let constant_overflow_mask = 1 << (length_width + 5);
    let prepared_mask = if add {
        enabled_mask | (1 << 3) | constant_overflow_mask
    } else {
        predicate_mask | (1 << 3) | constant_overflow_mask
    };
    let enable_prefix = lengths.compute_operation + lengths.toggle_enable;
    let prepare_end = enable_prefix + lengths.prepare_range;
    let suffix = lengths.unprepare_range + lengths.toggle_enable + lengths.uncompute_operation;
    let window_end = emitted_ops
        .checked_sub(suffix)
        .expect("remainder suffix exceeds complete kernel");
    assert!(prepare_end <= window_end);
    let unprepare_end = window_end + lengths.unprepare_range;
    let disable_end = unprepare_end + lengths.toggle_enable;
    assert_eq!(
        disable_end + lengths.uncompute_operation,
        emitted_ops,
        "remainder phase accounting must cover the complete kernel"
    );

    if add {
        vec![
            (lengths.compute_operation, predicate_mask),
            (enable_prefix, enabled_mask),
            (prepare_end, prepared_mask),
            (window_end, prepared_mask),
            (unprepare_end, enabled_mask),
            (disable_end, predicate_mask),
            (emitted_ops, 0),
        ]
    } else {
        vec![
            (lengths.compute_operation, predicate_mask),
            (prepare_end, prepared_mask),
            (window_end, prepared_mask),
            (unprepare_end, predicate_mask),
            (emitted_ops, 0),
        ]
    }
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

#[must_use]
pub fn exhaustive_phase_overlaid_remainder_scratch_check(
) -> ReferencePhaseOverlaidRemainderScratchProofReport {
    assert!(
        std::env::var_os(PHASE_OVERLAID_REMAINDER_SCRATCH_FLAG).is_none(),
        "phase-overlaid remainder proof must start with the feature flag absent"
    );
    assert_eq!(
        remainder_scratch_width(REFERENCE_LENGTH_WIDTH, REFERENCE_LENGTH_WIDTH),
        baseline_remainder_scratch_width(REFERENCE_LENGTH_WIDTH, REFERENCE_LENGTH_WIDTH),
        "the default-off route must retain the baseline scratch layout"
    );

    // Includes a one-bit edge case, shifted windows, and both lower and upper
    // window boundaries. Every data/control basis state is checked per row.
    let configurations = [
        (1usize, 1usize, 1usize, 1usize),
        (2, 1, 2, 1),
        (2, 2, 3, 2),
        (2, 3, 3, 3),
    ];
    let kernels = [
        RemainderKernel::Add,
        RemainderKernel::AddInverse,
        RemainderKernel::Sub,
        RemainderKernel::SubInverse,
    ];

    let mut basis_states_checked = 0usize;
    let mut baseline_equivalence_checks = 0usize;
    let mut inverse_pair_checks = 0usize;
    let mut phase_boundary_clean_checks = 0usize;
    let mut ancilla_clean_checks = 0usize;
    let mut control_combinations_checked = 0usize;
    let mut zero_remainder_checks = 0usize;
    let mut range_equality_checks = 0usize;
    let mut range_boundary_checks = 0usize;

    for &(length_width, work_width, total_work_width, window_upper) in &configurations {
        let baseline: Vec<B> = kernels
            .iter()
            .map(|&kernel| {
                build_remainder_kernel_for_layout(
                    work_width,
                    total_work_width,
                    window_upper,
                    length_width,
                    kernel,
                    RemainderScratchLayout::Baseline,
                )
            })
            .collect();
        let overlaid: Vec<B> = kernels
            .iter()
            .map(|&kernel| {
                build_remainder_kernel_for_layout(
                    work_width,
                    total_work_width,
                    window_upper,
                    length_width,
                    kernel,
                    RemainderScratchLayout::PhaseOverlaid,
                )
            })
            .collect();
        let boundaries: Vec<Vec<(usize, u64)>> = kernels
            .iter()
            .zip(&overlaid)
            .map(|(&kernel, builder)| {
                remainder_phase_boundaries(
                    builder.ops.len(),
                    total_work_width,
                    window_upper,
                    length_width,
                    kernel,
                )
            })
            .collect();

        let data_width = 3 + 2 * work_width + 4 * length_width;
        let scratch_width = phase_overlaid_remainder_scratch_width(length_width, length_width);
        assert!(data_width + scratch_width < u64::BITS as usize);
        let data_mask = (1u64 << data_width) - 1;
        let length_mask = (1u64 << length_width) - 1;
        let lengths_offset = 3 + 2 * work_width;
        let mut controls_seen = [false; 8];

        for input in 0..(1u64 << data_width) {
            controls_seen[(input & 0b111) as usize] = true;
            let l_t = (input >> lengths_offset) & length_mask;
            let l_q = (input >> (lengths_offset + length_width)) & length_mask;
            let l_s = (input >> (lengths_offset + 2 * length_width)) & length_mask;
            let l_r_prime = (input >> (lengths_offset + 3 * length_width)) & length_mask;
            let q_boundary = ((window_upper - 1) as u64) & length_mask;
            let s_boundary = ((total_work_width - window_upper + 1) as u64) & length_mask;
            let q_sum = l_t.wrapping_add(l_q) & length_mask;

            if l_r_prime == 0 {
                zero_remainder_checks += 1;
            }
            if q_sum == q_boundary || l_s == s_boundary {
                range_equality_checks += 1;
            }
            let q_before = q_boundary.wrapping_sub(1) & length_mask;
            let q_after = q_boundary.wrapping_add(1) & length_mask;
            let s_before = s_boundary.wrapping_sub(1) & length_mask;
            let s_after = s_boundary.wrapping_add(1) & length_mask;
            if [q_before, q_boundary, q_after].contains(&q_sum)
                || [s_before, s_boundary, s_after].contains(&l_s)
            {
                range_boundary_checks += 1;
            }

            let mut overlaid_outputs = [0u64; 4];
            for index in 0..kernels.len() {
                let baseline_output = apply_scalar(&baseline[index].ops, input);
                assert_eq!(
                    baseline_output >> data_width,
                    0,
                    "baseline remainder scratch dirty"
                );
                ancilla_clean_checks += 1;

                let mut overlaid_output = input;
                let mut previous_cut = 0usize;
                for &(cut, allowed_dirty) in &boundaries[index] {
                    assert!(previous_cut <= cut && cut <= overlaid[index].ops.len());
                    overlaid_output =
                        apply_scalar(&overlaid[index].ops[previous_cut..cut], overlaid_output);
                    let scratch_state = overlaid_output >> data_width;
                    assert_eq!(
                        scratch_state & !allowed_dirty,
                        0,
                        "phase boundary exposed a dirty lane outside its live set"
                    );
                    phase_boundary_clean_checks += 1;
                    previous_cut = cut;
                }
                assert_eq!(previous_cut, overlaid[index].ops.len());
                assert_eq!(
                    overlaid_output >> data_width,
                    0,
                    "phase-overlaid remainder scratch dirty"
                );
                assert_eq!(
                    overlaid_output & data_mask,
                    baseline_output & data_mask,
                    "phase-overlaid remainder kernel differs from baseline"
                );
                ancilla_clean_checks += 1;
                baseline_equivalence_checks += 1;
                basis_states_checked += 1;
                overlaid_outputs[index] = overlaid_output;
            }

            for (forward_index, inverse_index) in [(0usize, 1usize), (2, 3)] {
                assert_eq!(
                    apply_scalar(
                        &overlaid[inverse_index].ops,
                        overlaid_outputs[forward_index]
                    ),
                    input,
                    "phase-overlaid forward/inverse pair failed"
                );
                assert_eq!(
                    apply_scalar(
                        &overlaid[forward_index].ops,
                        overlaid_outputs[inverse_index]
                    ),
                    input,
                    "phase-overlaid inverse/forward pair failed"
                );
                inverse_pair_checks += 2;
            }
        }

        assert!(controls_seen.iter().all(|&seen| seen));
        control_combinations_checked += controls_seen.len();
    }

    assert!(zero_remainder_checks > 0);
    assert!(range_equality_checks > 0);
    assert!(range_boundary_checks > 0);

    let baseline_add257 = build_remainder_kernel_for_layout(
        257,
        259,
        259,
        REFERENCE_LENGTH_WIDTH,
        RemainderKernel::Add,
        RemainderScratchLayout::Baseline,
    );
    let overlaid_add257 = build_remainder_kernel_for_layout(
        257,
        259,
        259,
        REFERENCE_LENGTH_WIDTH,
        RemainderKernel::Add,
        RemainderScratchLayout::PhaseOverlaid,
    );
    let baseline_sub257 = build_remainder_kernel_for_layout(
        257,
        259,
        259,
        REFERENCE_LENGTH_WIDTH,
        RemainderKernel::Sub,
        RemainderScratchLayout::Baseline,
    );
    let overlaid_sub257 = build_remainder_kernel_for_layout(
        257,
        259,
        259,
        REFERENCE_LENGTH_WIDTH,
        RemainderKernel::Sub,
        RemainderScratchLayout::PhaseOverlaid,
    );
    for kernel in [RemainderKernel::AddInverse, RemainderKernel::SubInverse] {
        let baseline = build_remainder_kernel_for_layout(
            257,
            259,
            259,
            REFERENCE_LENGTH_WIDTH,
            kernel,
            RemainderScratchLayout::Baseline,
        );
        let overlaid = build_remainder_kernel_for_layout(
            257,
            259,
            259,
            REFERENCE_LENGTH_WIDTH,
            kernel,
            RemainderScratchLayout::PhaseOverlaid,
        );
        assert_eq!(gate_counts(&baseline.ops), gate_counts(&overlaid.ops));
    }

    let baseline_add_counts = gate_counts(&baseline_add257.ops);
    let overlaid_add_counts = gate_counts(&overlaid_add257.ops);
    let baseline_sub_counts = gate_counts(&baseline_sub257.ops);
    let overlaid_sub_counts = gate_counts(&overlaid_sub257.ops);
    assert_eq!(baseline_add_counts, overlaid_add_counts);
    assert_eq!(baseline_sub_counts, overlaid_sub_counts);

    let baseline_reference9_scratch_lanes =
        baseline_remainder_scratch_width(REFERENCE_LENGTH_WIDTH, REFERENCE_LENGTH_WIDTH);
    let overlaid_reference9_scratch_lanes =
        phase_overlaid_remainder_scratch_width(REFERENCE_LENGTH_WIDTH, REFERENCE_LENGTH_WIDTH);
    assert_eq!(baseline_reference9_scratch_lanes, 35);
    assert_eq!(overlaid_reference9_scratch_lanes, 15);
    assert_eq!(
        baseline_add257.peak_qubits as usize - overlaid_add257.peak_qubits as usize,
        baseline_reference9_scratch_lanes - overlaid_reference9_scratch_lanes
    );
    assert!(
        std::env::var_os(PHASE_OVERLAID_REMAINDER_SCRATCH_FLAG).is_none(),
        "phase-overlaid remainder proof must restore the default-off environment"
    );

    ReferencePhaseOverlaidRemainderScratchProofReport {
        length_widths_checked: 2,
        window_configurations_checked: configurations.len(),
        basis_states_checked,
        baseline_equivalence_checks,
        inverse_pair_checks,
        phase_boundary_clean_checks,
        ancilla_clean_checks,
        control_combinations_checked,
        zero_remainder_checks,
        range_equality_checks,
        range_boundary_checks,
        baseline_reference9_scratch_lanes,
        overlaid_reference9_scratch_lanes,
        scratch_lanes_saved: baseline_reference9_scratch_lanes - overlaid_reference9_scratch_lanes,
        baseline_reference9_peak_qubits: baseline_add257.peak_qubits as usize,
        overlaid_reference9_peak_qubits: overlaid_add257.peak_qubits as usize,
        baseline_add257: baseline_add_counts,
        overlaid_add257: overlaid_add_counts,
        baseline_sub257: baseline_sub_counts,
        overlaid_sub257: overlaid_sub_counts,
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
    let scratch = if rotated_bitlen_scratch_reuse_requested() {
        circ.alloc_qreg_bits("rs.swap-length.borrowed-rotated-bitlen", 3)
    } else {
        Vec::new()
    };
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

fn q847_lifetime_local_resources(builder: &B) -> Q847LifetimeLocalResources {
    let counts = measurement_classical_gate_counts(&builder.ops);
    let active_qubits = builder.active_qubits as usize;
    let peak_qubits = builder.peak_qubits as usize;
    Q847LifetimeLocalResources {
        active_qubits,
        peak_qubits,
        temporary_qubits: peak_qubits - active_qubits,
        emitted_ops: builder.ops.len(),
        emitted_toffoli: counts.ccx,
    }
}

fn assert_q847_default_stream_identity(configured: &B, direct: &B) {
    assert_eq!(configured.ops, direct.ops);
    assert_eq!(configured.next_qubit, direct.next_qubit);
    assert_eq!(configured.next_bit, direct.next_bit);
    assert_eq!(configured.active_qubits, direct.active_qubits);
    assert_eq!(configured.peak_qubits, direct.peak_qubits);
    assert_eq!(configured.free_qubits, direct.free_qubits);
    assert_eq!(configured.allocation_serial, direct.allocation_serial);
}

fn q847_basis_input(data_ids: &[u32], value: usize) -> u64 {
    data_ids.iter().enumerate().fold(0u64, |state, (bit, id)| {
        state | ((((value >> bit) & 1) as u64) << id)
    })
}

fn verify_q847_simulator_equivalence_impl(
    label: &[u8],
    baseline: &B,
    candidate: &B,
    data_ids: &[u32],
    external_mask: u64,
    require_zero_phase: bool,
) -> (usize, usize, usize) {
    use crate::circuit::QubitId;
    use crate::sim::Simulator;
    use sha3::{
        digest::{ExtendableOutput, Update},
        Shake128,
    };

    assert!(baseline.next_qubit <= 64);
    assert!(candidate.next_qubit <= 64);
    let states = 1usize << data_ids.len();
    let mut cases_checked = 0usize;
    let mut phase_clean_checks = 0usize;
    let mut ancilla_clean_checks = 0usize;
    for batch_start in (0..states).step_by(64) {
        let shots = (states - batch_start).min(64);
        let live = if shots == 64 {
            u64::MAX
        } else {
            (1u64 << shots) - 1
        };
        let mut baseline_seed = Shake128::default();
        baseline_seed.update(label);
        baseline_seed.update(&(batch_start as u64).to_le_bytes());
        let mut baseline_xof = baseline_seed.finalize_xof();
        let mut baseline_simulator = Simulator::new(
            baseline.next_qubit as usize,
            baseline.next_bit as usize,
            &mut baseline_xof,
        );
        let mut candidate_seed = Shake128::default();
        candidate_seed.update(label);
        candidate_seed.update(&(batch_start as u64).to_le_bytes());
        let mut candidate_xof = candidate_seed.finalize_xof();
        let mut candidate_simulator = Simulator::new(
            candidate.next_qubit as usize,
            candidate.next_bit as usize,
            &mut candidate_xof,
        );

        for shot in 0..shots {
            let value = batch_start + shot;
            for (bit, &id) in data_ids.iter().enumerate() {
                if (value >> bit) & 1 != 0 {
                    *baseline_simulator.qubit_mut(QubitId(u64::from(id))) |= 1u64 << shot;
                    *candidate_simulator.qubit_mut(QubitId(u64::from(id))) |= 1u64 << shot;
                }
            }
        }
        baseline_simulator.apply_iter(baseline.ops.iter());
        candidate_simulator.apply_iter(candidate.ops.iter());
        if require_zero_phase {
            assert_eq!(
                baseline_simulator.phase & live,
                0,
                "baseline {label:?} left phase garbage"
            );
            assert_eq!(
                candidate_simulator.phase & live,
                0,
                "candidate {label:?} left phase garbage"
            );
            phase_clean_checks += 2 * shots;
        }

        for id in 0..baseline.next_qubit {
            let baseline_value = baseline_simulator.qubit(QubitId(u64::from(id))) & live;
            if external_mask & (1u64 << id) != 0 {
                assert_eq!(
                    baseline_value,
                    candidate_simulator.qubit(QubitId(u64::from(id))) & live
                );
            } else {
                assert_eq!(baseline_value, 0, "baseline {label:?} left q{id} dirty");
            }
        }
        for id in 0..candidate.next_qubit {
            if external_mask & (1u64 << id) == 0 {
                assert_eq!(
                    candidate_simulator.qubit(QubitId(u64::from(id))) & live,
                    0,
                    "candidate {label:?} left q{id} dirty"
                );
            }
        }
        cases_checked += shots;
        ancilla_clean_checks += 2 * shots;
    }
    (cases_checked, phase_clean_checks, ancilla_clean_checks)
}

fn verify_q847_simulator_equivalence(
    label: &[u8],
    baseline: &B,
    candidate: &B,
    data_ids: &[u32],
    external_mask: u64,
) -> (usize, usize, usize) {
    verify_q847_simulator_equivalence_impl(
        label,
        baseline,
        candidate,
        data_ids,
        external_mask,
        true,
    )
}

fn verify_q847_selected_simulator_equivalence(
    label: &[u8],
    baseline: &B,
    candidate: &B,
    data_ids: &[u32],
    external_mask: u64,
    values: &[usize],
) -> (usize, usize, usize) {
    use crate::circuit::QubitId;
    use crate::sim::Simulator;
    use sha3::{
        digest::{ExtendableOutput, Update},
        Shake128,
    };

    assert!(baseline.next_qubit <= 64 && candidate.next_qubit <= 64);
    let mut cases_checked = 0usize;
    let mut phase_clean_checks = 0usize;
    let mut ancilla_clean_checks = 0usize;
    for (batch_index, batch) in values.chunks(64).enumerate() {
        let shots = batch.len();
        let live = if shots == 64 {
            u64::MAX
        } else {
            (1u64 << shots) - 1
        };
        let mut baseline_seed = Shake128::default();
        baseline_seed.update(label);
        baseline_seed.update(&(batch_index as u64).to_le_bytes());
        let mut baseline_xof = baseline_seed.finalize_xof();
        let mut baseline_simulator = Simulator::new(
            baseline.next_qubit as usize,
            baseline.next_bit as usize,
            &mut baseline_xof,
        );
        let mut candidate_seed = Shake128::default();
        candidate_seed.update(label);
        candidate_seed.update(&(batch_index as u64).to_le_bytes());
        let mut candidate_xof = candidate_seed.finalize_xof();
        let mut candidate_simulator = Simulator::new(
            candidate.next_qubit as usize,
            candidate.next_bit as usize,
            &mut candidate_xof,
        );
        for (shot, &value) in batch.iter().enumerate() {
            for (bit, &id) in data_ids.iter().enumerate() {
                if (value >> bit) & 1 != 0 {
                    *baseline_simulator.qubit_mut(QubitId(u64::from(id))) |= 1u64 << shot;
                    *candidate_simulator.qubit_mut(QubitId(u64::from(id))) |= 1u64 << shot;
                }
            }
        }
        baseline_simulator.apply_iter(baseline.ops.iter());
        candidate_simulator.apply_iter(candidate.ops.iter());
        assert_eq!(baseline_simulator.phase & live, 0, "baseline phase");
        assert_eq!(candidate_simulator.phase & live, 0, "candidate phase");
        for id in 0..baseline.next_qubit {
            let baseline_value = baseline_simulator.qubit(QubitId(u64::from(id))) & live;
            if external_mask & (1u64 << id) != 0 {
                assert_eq!(
                    baseline_value,
                    candidate_simulator.qubit(QubitId(u64::from(id))) & live
                );
            } else {
                assert_eq!(baseline_value, 0, "baseline {label:?} left q{id} dirty");
            }
        }
        for id in 0..candidate.next_qubit {
            if external_mask & (1u64 << id) == 0 {
                assert_eq!(
                    candidate_simulator.qubit(QubitId(u64::from(id))) & live,
                    0,
                    "candidate {label:?} left q{id} dirty"
                );
            }
        }
        cases_checked += shots;
        phase_clean_checks += 2 * shots;
        ancilla_clean_checks += 2 * shots;
    }
    (cases_checked, phase_clean_checks, ancilla_clean_checks)
}

fn verify_preserved_dy_top_borrow_windows(
    label: &[u8],
    candidate: &PromisedLqSwapProofHarness,
) -> (usize, usize, usize) {
    use crate::circuit::QubitId;
    use crate::sim::Simulator;
    use sha3::{
        digest::{ExtendableOutput, Update},
        Shake128,
    };

    assert!(candidate.builder.next_qubit <= 64);
    assert!(!candidate.preserved_dy_top_windows.is_empty());
    let states = 1usize << candidate.data_ids.len();
    let mut entry_checks = 0usize;
    let mut restore_checks = 0usize;
    let mut phase_checks = 0usize;
    for (window_index, window) in candidate.preserved_dy_top_windows.iter().enumerate() {
        assert_eq!(candidate.preserved_dy_top_mask, 1u64 << window.lender_id);
        for batch_start in (0..states).step_by(64) {
            let shots = (states - batch_start).min(64);
            let live = if shots == 64 {
                u64::MAX
            } else {
                (1u64 << shots) - 1
            };
            let mut seed = Shake128::default();
            seed.update(label);
            seed.update(&(window_index as u64).to_le_bytes());
            seed.update(&(batch_start as u64).to_le_bytes());
            let mut xof = seed.finalize_xof();
            let mut simulator = Simulator::new(
                candidate.builder.next_qubit as usize,
                candidate.builder.next_bit as usize,
                &mut xof,
            );
            for shot in 0..shots {
                let value = batch_start + shot;
                for (bit, &id) in candidate.data_ids.iter().enumerate() {
                    if (value >> bit) & 1 != 0 {
                        *simulator.qubit_mut(QubitId(u64::from(id))) |= 1u64 << shot;
                    }
                }
            }
            simulator.apply_iter(candidate.builder.ops[..window.entry_ops_idx].iter());
            assert_eq!(
                simulator.qubit(QubitId(u64::from(window.lender_id))) & live,
                0,
                "{label:?} lender dirty at borrow entry"
            );
            assert_eq!(
                simulator.phase & live,
                0,
                "{label:?} phase dirty at borrow entry"
            );
            simulator.apply_iter(
                candidate.builder.ops[window.entry_ops_idx..window.restore_ops_idx].iter(),
            );
            assert_eq!(
                simulator.qubit(QubitId(u64::from(window.lender_id))) & live,
                0,
                "{label:?} lender dirty after restore"
            );
            assert_eq!(
                simulator.phase & live,
                0,
                "{label:?} phase dirty after restore"
            );
            entry_checks += shots;
            restore_checks += shots;
            phase_checks += 2 * shots;
        }
    }
    (entry_checks, restore_checks, phase_checks)
}

fn preserved_dy_top_alias_rejections() -> usize {
    use std::panic::{catch_unwind, AssertUnwindSafe};

    fn build_rejection(kind: usize) {
        let mut circ = Circuit::new();
        let control = circ.alloc_qreg("q846.reject.control");
        let left_length = circ.alloc_qreg_bits("q846.reject.left-length", 2);
        let work = circ.alloc_qreg_bits("q846.reject.work", 3);
        let output = circ.alloc_qreg_bits("q846.reject.output", 2);
        let scratch = circ.alloc_qreg_bits("q846.reject.scratch", 3);
        let lender = match kind {
            0 => &control,
            1 => &left_length[0],
            2 => &work[0],
            3 => &output[0],
            4 => &scratch[0],
            _ => unreachable!(),
        };
        controlled_xor_rotated_suffix_bit_length(
            &mut circ,
            &control,
            &left_length,
            &work,
            &output,
            &scratch,
            &[lender],
        );
    }

    let previous_hook = std::panic::take_hook();
    std::panic::set_hook(Box::new(|_| {}));
    let mut rejected = 0usize;
    for kind in 0..5 {
        assert!(catch_unwind(AssertUnwindSafe(|| build_rejection(kind))).is_err());
        rejected += 1;
    }
    std::panic::set_hook(previous_hook);
    rejected
}

fn verify_q847_simulator_data_and_ancilla_equivalence(
    label: &[u8],
    baseline: &B,
    candidate: &B,
    data_ids: &[u32],
    external_mask: u64,
) -> (usize, usize) {
    let (cases, phases, ancillas) = verify_q847_simulator_equivalence_impl(
        label,
        baseline,
        candidate,
        data_ids,
        external_mask,
        false,
    );
    assert_eq!(phases, 0);
    (cases, ancillas)
}

fn verify_coefficient_add_lender_window_phase_clean(
    label: &[u8],
    candidate: &B,
    data_ids: &[u32],
    trace: CoefficientAddLenderTrace,
) -> usize {
    use crate::circuit::QubitId;
    use crate::sim::Simulator;
    use sha3::{
        digest::{ExtendableOutput, Update},
        Shake128,
    };

    assert!(candidate.next_qubit <= 64);
    assert_eq!(trace.calls, 1);
    let states = 1usize << data_ids.len();
    let mut phase_clean_checks = 0usize;
    for batch_start in (0..states).step_by(64) {
        let shots = (states - batch_start).min(64);
        let live = if shots == 64 {
            u64::MAX
        } else {
            (1u64 << shots) - 1
        };
        let mut seed = Shake128::default();
        seed.update(label);
        seed.update(&(batch_start as u64).to_le_bytes());
        let mut xof = seed.finalize_xof();
        let mut simulator = Simulator::new(
            candidate.next_qubit as usize,
            candidate.next_bit as usize,
            &mut xof,
        );
        for shot in 0..shots {
            let value = batch_start + shot;
            for (bit, &id) in data_ids.iter().enumerate() {
                if (value >> bit) & 1 != 0 {
                    *simulator.qubit_mut(QubitId(u64::from(id))) |= 1u64 << shot;
                }
            }
        }
        simulator.apply_iter(candidate.ops[..trace.entry_ops_idx].iter());
        assert_eq!(
            simulator.phase & live,
            0,
            "{label:?} dirty lender entry phase"
        );
        for id in 0..candidate.next_qubit {
            if trace.lender_mask & (1u64 << id) != 0 {
                assert_eq!(
                    simulator.qubit(QubitId(u64::from(id))) & live,
                    0,
                    "{label:?} dirty lender q{id} at entry"
                );
            }
        }
        simulator.apply_iter(candidate.ops[trace.entry_ops_idx..trace.restore_ops_idx].iter());
        assert_eq!(
            simulator.phase & live,
            0,
            "{label:?} dirty lender restore phase"
        );
        for id in 0..candidate.next_qubit {
            if trace.lender_mask & (1u64 << id) != 0 {
                assert_eq!(
                    simulator.qubit(QubitId(u64::from(id))) & live,
                    0,
                    "{label:?} dirty lender q{id} at restore"
                );
            }
        }
        phase_clean_checks += 2 * shots;
    }
    phase_clean_checks
}

fn configure_q847_lifetime_prerequisites(coefficient_loan: bool) {
    configure_raw_bit_length_loan_proof(RawBitLengthZeroMode::Fused, coefficient_loan);
    std::env::set_var(INPLACE_ROTATED_BITLEN_BOUNDARY_FLAG, "1");
    std::env::set_var(FUSED_PREFIX_SCRATCH_LOAN_FLAG, "1");
    std::env::remove_var(PRESERVED_DY_TOP_PREFIX_LOAN_FLAG);
    std::env::remove_var(MIXED_WIDTH_L_R_PRIME_FLAG);
    std::env::remove_var(PROMISED_LQ_SWAP_BORROW_FLAG);
    std::env::remove_var(PROMISED_SWAP_SUPPORT_LIFETIME_FUSION_FLAG);
    std::env::remove_var(Q845_SWAP_ONLY_T_PRIME_LENGTH_FLAG);
    std::env::remove_var(SPLIT_COEFFICIENT_ROTATION_LIFETIME_FLAG);
    std::env::remove_var(COEFFICIENT_LESS_THAN_LANE_REUSE_FLAG);
    std::env::remove_var(CLEAN_CHAIN_COEFFICIENT_ADD_LENDER_FLAG);
}

struct PromisedLqSwapProofHarness {
    builder: B,
    data_ids: Vec<u32>,
    l_t_prime_ids: Vec<u32>,
    l_q_ids: Vec<u32>,
    l_s_ids: Vec<u32>,
    preserved_dy_top_id: u32,
    external_mask: u64,
    iteration_mask: u64,
    l_q_mask: u64,
    l_s_mask: u64,
    preserved_dy_top_mask: u64,
    preserved_dy_top_windows: Vec<PreservedDyTopBorrowWindow>,
    q839_support_lender_windows: Vec<Q839SupportLenderBorrowWindow>,
    prefix_trace: FusedPrefixScratchLoanAllocationTrace,
}

fn build_promised_l_q_swap_proof_harness(
    work_width: usize,
    length_width: usize,
    inverse: bool,
    route: PromisedLqSwapRoute,
) -> PromisedLqSwapProofHarness {
    build_promised_l_q_swap_proof_harness_with_preserved_dy_top(
        work_width,
        length_width,
        inverse,
        route,
        false,
    )
}

fn build_promised_l_q_swap_proof_harness_with_preserved_dy_top(
    work_width: usize,
    length_width: usize,
    inverse: bool,
    route: PromisedLqSwapRoute,
    borrow_preserved_dy_top: bool,
) -> PromisedLqSwapProofHarness {
    build_promised_l_q_swap_proof_harness_with_widths(
        work_width,
        length_width,
        length_width,
        inverse,
        route,
        borrow_preserved_dy_top,
    )
}

fn build_promised_l_q_swap_proof_harness_with_widths(
    work_width: usize,
    length_width: usize,
    r_length_width: usize,
    inverse: bool,
    route: PromisedLqSwapRoute,
    borrow_preserved_dy_top: bool,
) -> PromisedLqSwapProofHarness {
    assert!(length_width > 0);
    assert_l_r_prime_metadata_width(length_width, r_length_width);
    let mut circ = Circuit::new();
    let iteration = circ.alloc_qreg("q847.swap.iteration");
    let work1 = circ.alloc_qreg_bits("q847.swap.work1", work_width);
    let work2 = circ.alloc_qreg_bits("q847.swap.work2", work_width);
    let l_t = circ.alloc_qreg_bits("q847.swap.l-t", length_width);
    let l_t_prime = circ.alloc_qreg_bits("q847.swap.l-t-prime", length_width);
    let l_q_width = if q845_swap_only_t_prime_length_requested() {
        r_length_width
    } else {
        length_width
    };
    let l_q = circ.alloc_qreg_bits("q847.swap.l-q", l_q_width);
    let l_s = circ.alloc_qreg_bits("q847.swap.l-s", length_width);
    let l_r_prime = circ.alloc_qreg_bits("q847.swap.l-r-prime", r_length_width);
    // This proof-local lane models either ec3.ty_ov or
    // rs.divider.dy-restored[256]. It is intentionally excluded from the
    // basis inputs, so every tested schedule enters with canonical dy clean.
    let preserved_dy_top = circ.alloc_qreg("q846.swap.preserved-dy-top");
    let data_ids: Vec<u32> = std::iter::once(&iteration)
        .chain(&work1)
        .chain(&work2)
        .chain(&l_t)
        .chain(&l_t_prime)
        .chain(&l_q)
        .chain(&l_s)
        .chain(&l_r_prime)
        .map(QReg::id)
        .collect();
    let external_mask = qreg_mask(
        std::iter::once(&iteration)
            .chain(&work1)
            .chain(&work2)
            .chain(&l_t)
            .chain(&l_t_prime)
            .chain(&l_q)
            .chain(&l_s)
            .chain(&l_r_prime),
    ) | 1u64.checked_shl(preserved_dy_top.id()).unwrap_or(0);
    let iteration_mask = 1u64 << iteration.id();
    let l_q_mask = qreg_mask(&l_q);
    let l_s_mask = qreg_mask(&l_s);
    let preserved_dy_top_mask = 1u64.checked_shl(preserved_dy_top.id()).unwrap_or(0);
    let l_t_prime_ids = l_t_prime.iter().map(QReg::id).collect::<Vec<_>>();
    let l_q_ids = l_q.iter().map(QReg::id).collect::<Vec<_>>();
    let l_s_ids = l_s.iter().map(QReg::id).collect::<Vec<_>>();
    let preserved_dy_top_id = preserved_dy_top.id();

    let predicate_chain_width = 3usize.max(length_width.saturating_sub(2));
    let condition_scratch =
        circ.alloc_qreg_bits("q847.swap.zero-predicate", 3 + predicate_chain_width);
    let zero_q = &condition_scratch[0];
    let zero_s = &condition_scratch[1];
    let control = &condition_scratch[2];
    let chain = &condition_scratch[3..];
    compute_zero(&mut circ, &l_q, zero_q, chain);
    compute_zero(&mut circ, &l_s, zero_s, chain);
    let preserved_dy_top_scratch = borrow_preserved_dy_top
        .then_some(&preserved_dy_top)
        .into_iter()
        .collect::<Vec<_>>();
    begin_preserved_dy_top_borrow_trace();
    begin_q839_support_lender_borrow_trace();
    begin_fused_prefix_scratch_loan_allocation_trace();
    conditional_work_and_length_swap_under_zero_predicate(
        &mut circ,
        zero_q,
        zero_s,
        control,
        &iteration,
        &work1,
        &work2,
        &l_t,
        &l_t_prime,
        &l_q,
        &l_s,
        &l_r_prime,
        (1, work_width),
        (1, work_width),
        chain,
        &preserved_dy_top_scratch,
        inverse,
        route,
    );
    let prefix_trace = finish_fused_prefix_scratch_loan_allocation_trace();
    let q839_support_lender_windows = finish_q839_support_lender_borrow_trace();
    let preserved_dy_top_windows = finish_preserved_dy_top_borrow_trace();
    uncompute_zero(&mut circ, &l_s, zero_s, chain);
    uncompute_zero(&mut circ, &l_q, zero_q, chain);
    free_clean(&mut circ, condition_scratch);

    PromisedLqSwapProofHarness {
        builder: circ.into_builder(),
        data_ids,
        l_t_prime_ids,
        l_q_ids,
        l_s_ids,
        preserved_dy_top_id,
        external_mask,
        iteration_mask,
        l_q_mask,
        l_s_mask,
        preserved_dy_top_mask,
        preserved_dy_top_windows,
        q839_support_lender_windows,
        prefix_trace,
    }
}

/// Exhaustively check the production zero predicate around the promised
/// `l_q` swap lender, including arbitrary `l_q` values whenever control is off.
#[doc(hidden)]
#[must_use]
pub fn exhaustive_promised_l_q_swap_borrow_check() -> PromisedLqSwapBorrowProofReport {
    assert!(
        std::env::var_os(PROMISED_LQ_SWAP_BORROW_FLAG).is_none(),
        "the promised l_q swap borrow must default off"
    );
    let _environment = RawBitLengthProofEnvironment::capture();
    configure_q847_lifetime_prerequisites(true);
    let configurations = [(1usize, 1usize), (2, 2)];
    let mut basis_states_checked = 0usize;
    let mut scalar_equivalence_checks = 0usize;
    let mut simulator_equivalence_checks = 0usize;
    let mut control_off_checks = 0usize;
    let mut control_off_nonzero_l_q_checks = 0usize;
    let mut control_off_invalid_support_checks = 0usize;
    let mut control_on_promised_support_checks = 0usize;
    let mut lender_restore_checks = 0usize;
    let mut inverse_pair_checks = 0usize;
    let mut phase_clean_checks = 0usize;
    let mut ancilla_clean_checks = 0usize;

    for &(work_width, length_width) in &configurations {
        std::env::remove_var(PROMISED_LQ_SWAP_BORROW_FLAG);
        let legacy = build_promised_l_q_swap_proof_harness(
            work_width,
            length_width,
            false,
            PromisedLqSwapRoute::Allocated,
        );
        let baseline = build_promised_l_q_swap_proof_harness(
            work_width,
            length_width,
            false,
            PromisedLqSwapRoute::PromisedAllocated,
        );
        std::env::set_var(PROMISED_LQ_SWAP_BORROW_FLAG, "1");
        let candidate = build_promised_l_q_swap_proof_harness(
            work_width,
            length_width,
            false,
            PromisedLqSwapRoute::Configured,
        );
        let candidate_inverse = build_promised_l_q_swap_proof_harness(
            work_width,
            length_width,
            true,
            PromisedLqSwapRoute::Configured,
        );
        assert_eq!(baseline.data_ids, candidate.data_ids);
        assert_eq!(baseline.external_mask, candidate.external_mask);
        let (simulator_cases, phases, ancillas) = verify_q847_simulator_equivalence(
            b"q847-promised-lq-swap",
            &baseline.builder,
            &candidate.builder,
            &baseline.data_ids,
            baseline.external_mask,
        );
        simulator_equivalence_checks += simulator_cases;
        phase_clean_checks += phases;
        ancilla_clean_checks += ancillas;

        for value in 0..(1usize << baseline.data_ids.len()) {
            let input = q847_basis_input(&baseline.data_ids, value);
            let baseline_output = apply_scalar(&baseline.builder.ops, input);
            let candidate_output = apply_scalar(&candidate.builder.ops, input);
            assert_eq!(candidate_output, baseline_output);
            assert_eq!(candidate_output & !candidate.external_mask, 0);
            assert_eq!(
                candidate_output & candidate.l_q_mask,
                input & candidate.l_q_mask
            );
            let control_on = (candidate_output ^ input) & candidate.iteration_mask != 0;
            if control_on {
                assert_eq!(
                    candidate_output,
                    apply_scalar(&legacy.builder.ops, input),
                    "support-qualified borrow diverged from the legacy swap"
                );
                control_on_promised_support_checks += 1;
            } else {
                assert_eq!(candidate_output, input, "control-off swap must be identity");
                control_off_checks += 1;
                control_off_nonzero_l_q_checks += usize::from(input & candidate.l_q_mask != 0);
                control_off_invalid_support_checks +=
                    usize::from(input & candidate.l_q_mask == 0 && input & candidate.l_s_mask == 0);
            }
            assert_eq!(
                apply_scalar(&candidate_inverse.builder.ops, candidate_output),
                input
            );
            basis_states_checked += 1;
            scalar_equivalence_checks += 1;
            lender_restore_checks += 1;
            inverse_pair_checks += 1;
        }
    }

    std::env::remove_var(PROMISED_LQ_SWAP_BORROW_FLAG);
    let configured_default =
        build_promised_l_q_swap_proof_harness(2, 2, false, PromisedLqSwapRoute::Configured);
    let direct_default =
        build_promised_l_q_swap_proof_harness(2, 2, false, PromisedLqSwapRoute::Allocated);
    assert_q847_default_stream_identity(&configured_default.builder, &direct_default.builder);
    let configured_default_inverse =
        build_promised_l_q_swap_proof_harness(2, 2, true, PromisedLqSwapRoute::Configured);
    let direct_default_inverse =
        build_promised_l_q_swap_proof_harness(2, 2, true, PromisedLqSwapRoute::Allocated);
    assert_q847_default_stream_identity(
        &configured_default_inverse.builder,
        &direct_default_inverse.builder,
    );

    let baseline_local = q847_lifetime_local_resources(
        &build_promised_l_q_swap_proof_harness(
            259,
            REFERENCE_LENGTH_WIDTH,
            false,
            PromisedLqSwapRoute::Allocated,
        )
        .builder,
    );
    let candidate_local = q847_lifetime_local_resources(
        &build_promised_l_q_swap_proof_harness(
            259,
            REFERENCE_LENGTH_WIDTH,
            false,
            PromisedLqSwapRoute::BorrowLq,
        )
        .builder,
    );
    assert_eq!(baseline_local.active_qubits, candidate_local.active_qubits);
    assert_eq!(baseline_local.peak_qubits - candidate_local.peak_qubits, 8);
    assert!(candidate_local.emitted_ops > baseline_local.emitted_ops);
    assert!(candidate_local.emitted_toffoli > baseline_local.emitted_toffoli);
    let whole_point_add_invocations = 4 * REFERENCE_STEPS.div_ceil(4);
    let local_ops_delta = candidate_local.emitted_ops as i64 - baseline_local.emitted_ops as i64;
    let local_toffoli_delta =
        candidate_local.emitted_toffoli as i64 - baseline_local.emitted_toffoli as i64;

    PromisedLqSwapBorrowProofReport {
        configurations_checked: configurations.len(),
        basis_states_checked,
        scalar_equivalence_checks,
        simulator_equivalence_checks,
        control_off_checks,
        control_off_nonzero_l_q_checks,
        control_off_invalid_support_checks,
        control_on_promised_support_checks,
        lender_restore_checks,
        inverse_pair_checks,
        phase_clean_checks,
        ancilla_clean_checks,
        default_stream_identity_checks: 2,
        reference_lender_lanes: REFERENCE_LENGTH_WIDTH,
        reset_ops_removed_per_invocation: REFERENCE_LENGTH_WIDTH - 1,
        whole_point_add_invocations,
        whole_point_add_ops_delta: whole_point_add_invocations as i64 * local_ops_delta,
        whole_point_add_toffoli_delta: whole_point_add_invocations as i64 * local_toffoli_delta,
        baseline_local,
        candidate_local,
        local_qubit_delta: candidate_local.peak_qubits as i64 - baseline_local.peak_qubits as i64,
        local_ops_delta,
        local_toffoli_delta,
    }
}

/// Prove that the support discrepancy can remain in `l_q` across the promised
/// swap. The support-qualified branch has discrepancy zero; every gate that
/// temporarily uses `l_q` is controlled by that branch, so arbitrary
/// discrepancy values on the control-off branch are preserved exactly.
#[doc(hidden)]
#[must_use]
pub fn exhaustive_promised_swap_support_lifetime_fusion_check(
) -> PromisedSwapSupportLifetimeFusionProofReport {
    assert!(
        std::env::var_os(PROMISED_SWAP_SUPPORT_LIFETIME_FUSION_FLAG).is_none(),
        "the promised-swap support lifetime fusion must default off"
    );
    let _environment = RawBitLengthProofEnvironment::capture();
    configure_q847_lifetime_prerequisites(true);
    std::env::set_var(PROMISED_LQ_SWAP_BORROW_FLAG, "1");

    let configurations = [(1usize, 1usize, 1usize), (2, 2, 2), (2, 2, 1)];
    let mut basis_states_checked = 0usize;
    let mut scalar_equivalence_checks = 0usize;
    let mut simulator_equivalence_checks = 0usize;
    let mut control_off_checks = 0usize;
    let mut control_off_nonzero_l_q_checks = 0usize;
    let mut control_on_promised_support_checks = 0usize;
    let mut excluded_mixed_width_overflow_states = 0usize;
    let mut lender_restore_checks = 0usize;
    let mut inverse_pair_checks = 0usize;
    let mut phase_clean_checks = 0usize;
    let mut ancilla_clean_checks = 0usize;

    for &(work_width, length_width, r_length_width) in &configurations {
        std::env::remove_var(PROMISED_SWAP_SUPPORT_LIFETIME_FUSION_FLAG);
        let baseline = build_promised_l_q_swap_proof_harness_with_widths(
            work_width,
            length_width,
            r_length_width,
            false,
            PromisedLqSwapRoute::BorrowLq,
            false,
        );
        let baseline_inverse = build_promised_l_q_swap_proof_harness_with_widths(
            work_width,
            length_width,
            r_length_width,
            true,
            PromisedLqSwapRoute::BorrowLq,
            false,
        );
        std::env::set_var(PROMISED_SWAP_SUPPORT_LIFETIME_FUSION_FLAG, "1");
        let candidate = build_promised_l_q_swap_proof_harness_with_widths(
            work_width,
            length_width,
            r_length_width,
            false,
            PromisedLqSwapRoute::BorrowLq,
            false,
        );
        let candidate_inverse = build_promised_l_q_swap_proof_harness_with_widths(
            work_width,
            length_width,
            r_length_width,
            true,
            PromisedLqSwapRoute::BorrowLq,
            false,
        );
        assert_eq!(baseline.data_ids, candidate.data_ids);
        assert_eq!(baseline.external_mask, candidate.external_mask);

        let mut supported_values = Vec::new();
        let mut supported_inverse_values = Vec::new();
        for value in 0..(1usize << baseline.data_ids.len()) {
            let input = q847_basis_input(&baseline.data_ids, value);
            let baseline_output = apply_scalar(&baseline.builder.ops, input);
            let candidate_output = apply_scalar(&candidate.builder.ops, input);
            let baseline_control_on =
                (baseline_output ^ input) & baseline.iteration_mask != 0;
            let baseline_lender_restored =
                baseline_output & baseline.l_q_mask == input & baseline.l_q_mask;
            if baseline_control_on && !baseline_lender_restored {
                // With r_length_width = length_width - 1, arbitrary reduced
                // states can request a suffix length outside the representable
                // metadata range. The production mixed-width route separately
                // proves that this overflow is unreachable.
                excluded_mixed_width_overflow_states += 1;
                basis_states_checked += 1;
                continue;
            }
            assert_eq!(
                candidate_output,
                baseline_output,
                "support-lifetime forward mismatch work={work_width} len={length_width} rlen={r_length_width} value={value} input={input:#x} baseline={baseline_output:#x} candidate={candidate_output:#x} lq={:#x} ls={:#x}",
                input & candidate.l_q_mask,
                input & candidate.l_s_mask,
            );
            assert_eq!(candidate_output & !candidate.external_mask, 0);
            assert_eq!(candidate_output & candidate.l_q_mask, input & candidate.l_q_mask);

            let baseline_inverse_output = apply_scalar(&baseline_inverse.builder.ops, input);
            let candidate_inverse_output = apply_scalar(&candidate_inverse.builder.ops, input);
            assert_eq!(candidate_inverse_output, baseline_inverse_output);
            assert_eq!(
                apply_scalar(&candidate_inverse.builder.ops, candidate_output),
                input
            );

            if baseline_control_on {
                supported_values.push(value);
                supported_inverse_values.push(
                    candidate
                        .data_ids
                        .iter()
                        .enumerate()
                        .fold(0usize, |packed, (bit, id)| {
                            packed | ((((candidate_output >> id) & 1) as usize) << bit)
                        }),
                );
                control_on_promised_support_checks += 1;
            } else {
                control_off_checks += 1;
                control_off_nonzero_l_q_checks +=
                    usize::from(input & candidate.l_q_mask != 0);
            }
            basis_states_checked += 1;
            scalar_equivalence_checks += 2;
            lender_restore_checks += 1;
            inverse_pair_checks += 1;
        }
        let (_, supported_phases, supported_ancillas) =
            verify_q847_selected_simulator_equivalence(
                b"support-phase-clean",
                &baseline.builder,
                &candidate.builder,
                &baseline.data_ids,
                baseline.external_mask,
                &supported_values,
            );
        simulator_equivalence_checks += supported_values.len();
        phase_clean_checks += supported_phases;
        ancilla_clean_checks += supported_ancillas;
        let (_, inverse_phases, inverse_ancillas) = verify_q847_selected_simulator_equivalence(
            b"support-inverse-phase-clean",
            &baseline_inverse.builder,
            &candidate_inverse.builder,
            &baseline_inverse.data_ids,
            baseline_inverse.external_mask,
            &supported_inverse_values,
        );
        simulator_equivalence_checks += supported_inverse_values.len();
        phase_clean_checks += inverse_phases;
        ancilla_clean_checks += inverse_ancillas;
    }

    std::env::remove_var(PROMISED_SWAP_SUPPORT_LIFETIME_FUSION_FLAG);
    let default_forward = build_promised_l_q_swap_proof_harness_with_widths(
        2,
        2,
        1,
        false,
        PromisedLqSwapRoute::BorrowLq,
        false,
    );
    let default_inverse = build_promised_l_q_swap_proof_harness_with_widths(
        2,
        2,
        1,
        true,
        PromisedLqSwapRoute::BorrowLq,
        false,
    );
    std::env::set_var(PROMISED_SWAP_SUPPORT_LIFETIME_FUSION_FLAG, "0");
    let explicit_off_forward = build_promised_l_q_swap_proof_harness_with_widths(
        2,
        2,
        1,
        false,
        PromisedLqSwapRoute::BorrowLq,
        false,
    );
    let explicit_off_inverse = build_promised_l_q_swap_proof_harness_with_widths(
        2,
        2,
        1,
        true,
        PromisedLqSwapRoute::BorrowLq,
        false,
    );
    assert_q847_default_stream_identity(&default_forward.builder, &explicit_off_forward.builder);
    assert_q847_default_stream_identity(&default_inverse.builder, &explicit_off_inverse.builder);

    std::env::remove_var(PROMISED_SWAP_SUPPORT_LIFETIME_FUSION_FLAG);
    let baseline_local = q847_lifetime_local_resources(
        &build_promised_l_q_swap_proof_harness_with_widths(
            259,
            REFERENCE_LENGTH_WIDTH,
            REFERENCE_R_LENGTH_WIDTH,
            false,
            PromisedLqSwapRoute::BorrowLq,
            false,
        )
        .builder,
    );
    std::env::set_var(PROMISED_SWAP_SUPPORT_LIFETIME_FUSION_FLAG, "1");
    let candidate_local = q847_lifetime_local_resources(
        &build_promised_l_q_swap_proof_harness_with_widths(
            259,
            REFERENCE_LENGTH_WIDTH,
            REFERENCE_R_LENGTH_WIDTH,
            false,
            PromisedLqSwapRoute::BorrowLq,
            false,
        )
        .builder,
    );
    let local_qubit_delta = candidate_local.peak_qubits as i64 - baseline_local.peak_qubits as i64;
    let local_ops_delta = candidate_local.emitted_ops as i64 - baseline_local.emitted_ops as i64;
    let local_toffoli_delta =
        candidate_local.emitted_toffoli as i64 - baseline_local.emitted_toffoli as i64;
    assert_eq!(local_qubit_delta, 0);
    assert!(local_ops_delta < 0);
    assert!(local_toffoli_delta < 0);

    PromisedSwapSupportLifetimeFusionProofReport {
        configurations_checked: configurations.len(),
        basis_states_checked,
        scalar_equivalence_checks,
        simulator_equivalence_checks,
        control_off_checks,
        control_off_nonzero_l_q_checks,
        control_on_promised_support_checks,
        excluded_mixed_width_overflow_states,
        lender_restore_checks,
        inverse_pair_checks,
        phase_clean_checks,
        ancilla_clean_checks,
        default_stream_identity_checks: 2,
        baseline_local,
        candidate_local,
        local_qubit_delta,
        local_ops_delta,
        local_toffoli_delta,
    }
}

/// Exhaustively compare the preserved canonical `dy[256]` lender against the
/// independent owned fused-prefix scratch route at reduced widths.
#[doc(hidden)]
#[must_use]
pub fn exhaustive_preserved_dy_top_prefix_loan_check() -> PreservedDyTopPrefixLoanProofReport {
    assert!(
        std::env::var_os(PRESERVED_DY_TOP_PREFIX_LOAN_FLAG).is_none(),
        "the preserved dy top prefix loan must default off"
    );
    let _environment = RawBitLengthProofEnvironment::capture();
    configure_q847_lifetime_prerequisites(true);
    std::env::set_var(PROMISED_LQ_SWAP_BORROW_FLAG, "1");

    let configurations = [(1usize, 1usize), (2, 2)];
    let mut basis_states_checked = 0usize;
    let mut scalar_equivalence_checks = 0usize;
    let mut simulator_equivalence_checks = 0usize;
    let mut control_off_checks = 0usize;
    let mut lender_clean_entry_checks = 0usize;
    let mut lender_restore_checks = 0usize;
    let mut inverse_pair_checks = 0usize;
    let mut phase_clean_checks = 0usize;
    let mut ancilla_clean_checks = 0usize;
    let mut borrow_windows_checked = 0usize;

    for &(work_width, length_width) in &configurations {
        let mut candidates = Vec::with_capacity(2);
        let mut baselines = Vec::with_capacity(2);
        for &inverse in &[false, true] {
            let baseline = build_promised_l_q_swap_proof_harness_with_preserved_dy_top(
                work_width,
                length_width,
                inverse,
                PromisedLqSwapRoute::Configured,
                false,
            );
            let candidate = build_promised_l_q_swap_proof_harness_with_preserved_dy_top(
                work_width,
                length_width,
                inverse,
                PromisedLqSwapRoute::Configured,
                true,
            );
            assert_eq!(baseline.data_ids, candidate.data_ids);
            assert_eq!(baseline.external_mask, candidate.external_mask);
            assert!(baseline.preserved_dy_top_windows.is_empty());
            assert_eq!(baseline.prefix_trace.calls, candidate.prefix_trace.calls);
            assert!(baseline.prefix_trace.calls > 0);
            assert_eq!(
                baseline.prefix_trace.maximum_owned_lanes,
                candidate.prefix_trace.maximum_owned_lanes + 1
            );
            assert_eq!(
                candidate.prefix_trace.maximum_borrowed_lanes,
                baseline.prefix_trace.maximum_borrowed_lanes + 1
            );
            borrow_windows_checked += candidate.preserved_dy_top_windows.len();

            let label = if inverse {
                b"q846-preserved-dy-top-inverse".as_slice()
            } else {
                b"q846-preserved-dy-top-forward".as_slice()
            };
            let (simulator_cases, phases, ancillas) = verify_q847_simulator_equivalence(
                label,
                &baseline.builder,
                &candidate.builder,
                &baseline.data_ids,
                baseline.external_mask,
            );
            simulator_equivalence_checks += simulator_cases;
            phase_clean_checks += phases;
            ancilla_clean_checks += ancillas;
            let (entries, restores, window_phases) =
                verify_preserved_dy_top_borrow_windows(label, &candidate);
            lender_clean_entry_checks += entries;
            lender_restore_checks += restores;
            phase_clean_checks += window_phases;

            for value in 0..(1usize << baseline.data_ids.len()) {
                let input = q847_basis_input(&baseline.data_ids, value);
                let baseline_output = apply_scalar(&baseline.builder.ops, input);
                let candidate_output = apply_scalar(&candidate.builder.ops, input);
                assert_eq!(candidate_output, baseline_output);
                assert_eq!(candidate_output & !candidate.external_mask, 0);
                assert_eq!(candidate_output & candidate.preserved_dy_top_mask, 0);
                if (candidate_output ^ input) & candidate.iteration_mask == 0 {
                    assert_eq!(
                        candidate_output, input,
                        "control-off route must be identity"
                    );
                    control_off_checks += 1;
                }
                basis_states_checked += 1;
                scalar_equivalence_checks += 1;
                lender_restore_checks += 1;
            }
            baselines.push(baseline);
            candidates.push(candidate);
        }

        let candidate_forward = &candidates[0];
        let candidate_inverse = &candidates[1];
        for value in 0..(1usize << candidate_forward.data_ids.len()) {
            let input = q847_basis_input(&candidate_forward.data_ids, value);
            let output = apply_scalar(&candidate_forward.builder.ops, input);
            assert_eq!(
                apply_scalar(&candidate_inverse.builder.ops, output),
                input,
                "preserved dy top forward/inverse pair failed"
            );
            inverse_pair_checks += 1;
        }
        drop(baselines);
    }

    let baseline_local = build_promised_l_q_swap_proof_harness_with_preserved_dy_top(
        259,
        REFERENCE_LENGTH_WIDTH,
        false,
        PromisedLqSwapRoute::Configured,
        false,
    );
    let candidate_local = build_promised_l_q_swap_proof_harness_with_preserved_dy_top(
        259,
        REFERENCE_LENGTH_WIDTH,
        false,
        PromisedLqSwapRoute::Configured,
        true,
    );
    let baseline_local_resources = q847_lifetime_local_resources(&baseline_local.builder);
    let candidate_local_resources = q847_lifetime_local_resources(&candidate_local.builder);
    assert_eq!(baseline_local.prefix_trace.maximum_owned_lanes, 2);
    assert_eq!(baseline_local.prefix_trace.maximum_borrowed_lanes, 7);
    assert_eq!(candidate_local.prefix_trace.maximum_owned_lanes, 1);
    assert_eq!(candidate_local.prefix_trace.maximum_borrowed_lanes, 8);
    let baseline_owned_prefix_lanes = baseline_local.prefix_trace.maximum_owned_lanes;
    let candidate_owned_prefix_lanes = candidate_local.prefix_trace.maximum_owned_lanes;
    assert_eq!(
        baseline_local_resources.active_qubits,
        candidate_local_resources.active_qubits
    );
    assert_eq!(
        baseline_local_resources.peak_qubits - candidate_local_resources.peak_qubits,
        1
    );
    assert_eq!(
        baseline_local_resources.emitted_toffoli,
        candidate_local_resources.emitted_toffoli
    );
    assert!(candidate_local_resources.emitted_ops < baseline_local_resources.emitted_ops);

    let alias_rejections = preserved_dy_top_alias_rejections();
    PreservedDyTopPrefixLoanProofReport {
        configurations_checked: configurations.len(),
        directions_checked: 2,
        basis_states_checked,
        scalar_equivalence_checks,
        simulator_equivalence_checks,
        control_off_checks,
        lender_clean_entry_checks,
        lender_restore_checks,
        inverse_pair_checks,
        phase_clean_checks,
        ancilla_clean_checks,
        borrow_windows_checked,
        alias_rejections,
        baseline_owned_prefix_lanes,
        candidate_owned_prefix_lanes,
        local_qubit_delta: candidate_local_resources.peak_qubits as i64
            - baseline_local_resources.peak_qubits as i64,
        local_ops_delta: candidate_local_resources.emitted_ops as i64
            - baseline_local_resources.emitted_ops as i64,
        local_toffoli_delta: candidate_local_resources.emitted_toffoli as i64
            - baseline_local_resources.emitted_toffoli as i64,
    }
}

struct MixedLrPrimeHarness {
    builder: B,
    data_ids: Vec<u32>,
    external_mask: u64,
    control_mask: u64,
    omitted_high_mask: u64,
}

fn build_mixed_l_r_prime_boundary_harness(
    source_width: usize,
    output_width: usize,
    mixed: bool,
    applications: usize,
) -> MixedLrPrimeHarness {
    assert!(output_width >= 2);
    let mut circ = Circuit::new();
    let control = circ.alloc_qreg("q845.boundary.control");
    let source = circ.alloc_qreg_bits("q845.boundary.source", source_width);
    let boundary = circ.alloc_qreg_bits("q845.boundary.l-r-prime-storage", output_width);
    let output = circ.alloc_qreg_bits("q845.boundary.output", output_width);
    let scratch = circ.alloc_qreg_bits("q845.boundary.scratch", 3);
    let logical_boundary = &boundary[..output_width - 1];
    let data_ids: Vec<u32> = std::iter::once(&control)
        .chain(&source)
        .chain(logical_boundary)
        .chain(&output)
        .map(QReg::id)
        .collect();
    let external_mask = qreg_mask(
        std::iter::once(&control)
            .chain(&source)
            .chain(&boundary)
            .chain(&output),
    );
    let source_refs: Vec<&QReg> = source.iter().collect();
    let scratch_refs: Vec<&QReg> = scratch.iter().collect();
    let boundary_view = if mixed { logical_boundary } else { &boundary };
    for _ in 0..applications {
        controlled_xor_saturating_bit_length_difference_with_route(
            &mut circ,
            &control,
            boundary_view,
            &source_refs,
            &output,
            &scratch_refs,
            None,
            SaturatingDifferenceBoundaryRoute::Inplace,
        );
    }
    free_clean(&mut circ, scratch);
    MixedLrPrimeHarness {
        builder: circ.into_builder(),
        data_ids,
        external_mask,
        control_mask: 1u64 << control.id(),
        omitted_high_mask: 1u64 << boundary[output_width - 1].id(),
    }
}

fn build_mixed_l_r_prime_swap_harness(
    work_width: usize,
    length_width: usize,
    mixed: bool,
    inverse: bool,
) -> MixedLrPrimeHarness {
    assert!(length_width >= 2);
    let mut circ = Circuit::new();
    let iteration = circ.alloc_qreg("q845.swap.iteration");
    let work1 = circ.alloc_qreg_bits("q845.swap.work1", work_width);
    let work2 = circ.alloc_qreg_bits("q845.swap.work2", work_width);
    let l_t = circ.alloc_qreg_bits("q845.swap.l-t", length_width);
    let l_t_prime = circ.alloc_qreg_bits("q845.swap.l-t-prime", length_width);
    let l_q = circ.alloc_qreg_bits("q845.swap.l-q", length_width);
    let l_s = circ.alloc_qreg_bits("q845.swap.l-s", length_width);
    let l_r_prime = circ.alloc_qreg_bits("q845.swap.l-r-prime-storage", length_width);
    let logical_l_r_prime = &l_r_prime[..length_width - 1];
    let data_ids: Vec<u32> = std::iter::once(&iteration)
        .chain(&work1)
        .chain(&work2)
        .chain(&l_t)
        .chain(&l_t_prime)
        .chain(&l_q)
        .chain(&l_s)
        .chain(logical_l_r_prime)
        .map(QReg::id)
        .collect();
    let external_mask = qreg_mask(
        std::iter::once(&iteration)
            .chain(&work1)
            .chain(&work2)
            .chain(&l_t)
            .chain(&l_t_prime)
            .chain(&l_q)
            .chain(&l_s)
            .chain(&l_r_prime),
    );
    let l_r_prime_view = if mixed { logical_l_r_prime } else { &l_r_prime };
    let predicate_chain_width = 3usize.max(length_width.saturating_sub(2));
    let condition_scratch =
        circ.alloc_qreg_bits("q845.swap.zero-predicate", 3 + predicate_chain_width);
    let zero_q = &condition_scratch[0];
    let zero_s = &condition_scratch[1];
    let control = &condition_scratch[2];
    let chain = &condition_scratch[3..];
    compute_zero(&mut circ, &l_q, zero_q, chain);
    compute_zero(&mut circ, &l_s, zero_s, chain);
    conditional_work_and_length_swap_under_zero_predicate(
        &mut circ,
        zero_q,
        zero_s,
        control,
        &iteration,
        &work1,
        &work2,
        &l_t,
        &l_t_prime,
        &l_q,
        &l_s,
        l_r_prime_view,
        (1, work_width),
        (1, work_width),
        chain,
        &[],
        inverse,
        PromisedLqSwapRoute::BorrowLq,
    );
    uncompute_zero(&mut circ, &l_s, zero_s, chain);
    uncompute_zero(&mut circ, &l_q, zero_q, chain);
    free_clean(&mut circ, condition_scratch);
    MixedLrPrimeHarness {
        builder: circ.into_builder(),
        data_ids,
        external_mask,
        control_mask: 1u64 << iteration.id(),
        omitted_high_mask: 1u64 << l_r_prime[length_width - 1].id(),
    }
}

/// Compare the eight-lane representation against an independent nine-lane
/// zero-high reference for signed decoding and forward/inverse metadata swaps.
#[doc(hidden)]
#[must_use]
pub fn exhaustive_mixed_l_r_prime_check() -> MixedLrPrimeProofReport {
    let _environment = RawBitLengthProofEnvironment::capture();
    configure_q847_lifetime_prerequisites(true);
    let boundary_configurations = [(2usize, 2usize), (4, 3), (5, 4)];
    let swap_configurations = [(1usize, 2usize), (2, 2)];
    let mut boundary_basis_states_checked = 0usize;
    let mut swap_basis_states_checked = 0usize;
    let mut scalar_equivalence_checks = 0usize;
    let mut simulator_equivalence_checks = 0usize;
    let mut control_off_checks = 0usize;
    let mut unsupported_control_on_states = 0usize;
    let mut omitted_high_lane_clean_checks = 0usize;
    let mut inverse_pair_checks = 0usize;
    let mut phase_clean_checks = 0usize;
    let mut ancilla_clean_checks = 0usize;

    for &(source_width, output_width) in &boundary_configurations {
        let baseline = build_mixed_l_r_prime_boundary_harness(source_width, output_width, false, 1);
        let candidate = build_mixed_l_r_prime_boundary_harness(source_width, output_width, true, 1);
        let roundtrip = build_mixed_l_r_prime_boundary_harness(source_width, output_width, true, 2);
        assert_eq!(baseline.data_ids, candidate.data_ids);
        let (cases, phases, ancillas) = verify_q847_simulator_equivalence(
            b"q845-mixed-lrp-boundary",
            &baseline.builder,
            &candidate.builder,
            &baseline.data_ids,
            baseline.external_mask,
        );
        simulator_equivalence_checks += cases;
        phase_clean_checks += phases;
        ancilla_clean_checks += ancillas;
        for value in 0..(1usize << baseline.data_ids.len()) {
            let input = q847_basis_input(&baseline.data_ids, value);
            let baseline_output = apply_scalar(&baseline.builder.ops, input);
            let candidate_output = apply_scalar(&candidate.builder.ops, input);
            assert_eq!(candidate_output, baseline_output);
            assert_eq!(candidate_output & candidate.omitted_high_mask, 0);
            assert_eq!(apply_scalar(&roundtrip.builder.ops, input), input);
            if input & candidate.control_mask == 0 {
                assert_eq!(candidate_output, input);
                control_off_checks += 1;
            }
            boundary_basis_states_checked += 1;
            scalar_equivalence_checks += 1;
            omitted_high_lane_clean_checks += 1;
            inverse_pair_checks += 1;
        }
    }

    for &(work_width, length_width) in &swap_configurations {
        let mut candidate_directions = Vec::with_capacity(2);
        let mut supported_directions = Vec::with_capacity(2);
        for &inverse in &[false, true] {
            let baseline =
                build_mixed_l_r_prime_swap_harness(work_width, length_width, false, inverse);
            let candidate =
                build_mixed_l_r_prime_swap_harness(work_width, length_width, true, inverse);
            assert_eq!(baseline.data_ids, candidate.data_ids);
            let label = if inverse {
                b"q845-mixed-lrp-swap-inverse".as_slice()
            } else {
                b"q845-mixed-lrp-swap-forward".as_slice()
            };
            // The legacy nine-lane circuit is the independent support oracle:
            // the mixed representation is valid exactly when its ninth output
            // lane remains zero. States rejected this way must have fired the
            // support-qualified control; every control-off state is retained.
            let mut supported_values = Vec::new();
            for value in 0..(1usize << baseline.data_ids.len()) {
                let input = q847_basis_input(&baseline.data_ids, value);
                let baseline_output = apply_scalar(&baseline.builder.ops, input);
                let control_on = (baseline_output ^ input) & baseline.control_mask != 0;
                if baseline_output & baseline.omitted_high_mask != 0 {
                    assert!(control_on);
                    unsupported_control_on_states += 1;
                    continue;
                }
                supported_values.push(value);
            }
            let (cases, phases, ancillas) = verify_q847_selected_simulator_equivalence(
                label,
                &baseline.builder,
                &candidate.builder,
                &baseline.data_ids,
                baseline.external_mask,
                &supported_values,
            );
            simulator_equivalence_checks += cases;
            phase_clean_checks += phases;
            ancilla_clean_checks += ancillas;
            for &value in &supported_values {
                let input = q847_basis_input(&baseline.data_ids, value);
                let baseline_output = apply_scalar(&baseline.builder.ops, input);
                let candidate_output = apply_scalar(&candidate.builder.ops, input);
                assert_eq!(candidate_output, baseline_output);
                assert_eq!(candidate_output & candidate.omitted_high_mask, 0);
                if (candidate_output ^ input) & candidate.control_mask == 0 {
                    assert_eq!(candidate_output, input);
                    control_off_checks += 1;
                }
                swap_basis_states_checked += 1;
                scalar_equivalence_checks += 1;
                omitted_high_lane_clean_checks += 1;
            }
            candidate_directions.push(candidate);
            supported_directions.push(supported_values);
        }
        for &value in &supported_directions[0] {
            let input = q847_basis_input(&candidate_directions[0].data_ids, value);
            let forward = apply_scalar(&candidate_directions[0].builder.ops, input);
            assert_eq!(
                apply_scalar(&candidate_directions[1].builder.ops, forward),
                input
            );
            inverse_pair_checks += 1;
        }
    }

    MixedLrPrimeProofReport {
        boundary_configurations_checked: boundary_configurations.len(),
        swap_configurations_checked: swap_configurations.len(),
        directions_checked: 2,
        boundary_basis_states_checked,
        swap_basis_states_checked,
        scalar_equivalence_checks,
        simulator_equivalence_checks,
        control_off_checks,
        unsupported_control_on_states,
        omitted_high_lane_clean_checks,
        inverse_pair_checks,
        phase_clean_checks,
        ancilla_clean_checks,
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum SplitRotationProofRoute {
    Configured,
    Continuous,
    Split,
}

struct SplitRotationProofHarness {
    builder: B,
    data_ids: Vec<u32>,
    external_mask: u64,
    control_mask: u64,
    lender_mask: u64,
    trace: RawBitLengthAllocationTrace,
}

fn build_split_rotation_proof_harness(
    length_width: usize,
    work_width: usize,
    loaned: bool,
    route: SplitRotationProofRoute,
) -> SplitRotationProofHarness {
    assert!(length_width > 0);
    assert!(work_width <= (1usize << length_width) - 1);
    let mut circ = Circuit::new();
    let control = circ.alloc_qreg("q847.rotation.control");
    let l_s = circ.alloc_qreg_bits("q847.rotation.l-s", length_width);
    let l_r_prime = circ.alloc_qreg_bits("q847.rotation.l-r-prime", length_width);
    let work2 = circ.alloc_qreg_bits("q847.rotation.work2", work_width);
    let output = circ.alloc_qreg_bits("q847.rotation.output", length_width);
    let chain = circ.alloc_qreg_bits("q847.rotation.chain", 2);
    let data_ids: Vec<u32> = std::iter::once(&control)
        .chain(&l_s)
        .chain(&l_r_prime)
        .chain(&work2)
        .chain(&output)
        .map(QReg::id)
        .collect();
    let external_mask = qreg_mask(
        std::iter::once(&control)
            .chain(&l_s)
            .chain(&l_r_prime)
            .chain(&work2)
            .chain(&output)
            .chain(&chain),
    );
    let control_mask = 1u64 << control.id();
    let lender_mask = qreg_mask(&chain);

    begin_raw_bit_length_allocation_trace();
    match route {
        SplitRotationProofRoute::Configured => controlled_xor_raw_t_prime_bit_length(
            &mut circ, &control, &l_s, &l_r_prime, &work2, &output, &chain,
        ),
        SplitRotationProofRoute::Continuous if loaned => {
            controlled_xor_raw_t_prime_bit_length_loaned_with_lifetime(
                &mut circ,
                &control,
                &l_s,
                &l_r_prime,
                &work2,
                &output,
                &chain,
                RawTPrimeRotationLifetime::Continuous,
            )
        }
        SplitRotationProofRoute::Split if loaned => {
            controlled_xor_raw_t_prime_bit_length_loaned_with_lifetime(
                &mut circ,
                &control,
                &l_s,
                &l_r_prime,
                &work2,
                &output,
                &chain,
                RawTPrimeRotationLifetime::Split,
            )
        }
        SplitRotationProofRoute::Continuous => {
            controlled_xor_raw_t_prime_bit_length_allocated_with_lifetime(
                &mut circ,
                &control,
                &l_s,
                &l_r_prime,
                &work2,
                &output,
                RawTPrimeRotationLifetime::Continuous,
            )
        }
        SplitRotationProofRoute::Split => {
            controlled_xor_raw_t_prime_bit_length_allocated_with_lifetime(
                &mut circ,
                &control,
                &l_s,
                &l_r_prime,
                &work2,
                &output,
                RawTPrimeRotationLifetime::Split,
            )
        }
    }
    let trace = finish_raw_bit_length_allocation_trace();
    SplitRotationProofHarness {
        builder: circ.into_builder(),
        data_ids,
        external_mask,
        control_mask,
        lender_mask,
        trace,
    }
}

/// Exhaustively prove that releasing and recomputing the physical rotation
/// around raw `t'` bit-length extraction preserves the reversible map.
#[doc(hidden)]
#[must_use]
pub fn exhaustive_split_coefficient_rotation_lifetime_check(
) -> SplitCoefficientRotationLifetimeProofReport {
    assert!(
        std::env::var_os(SPLIT_COEFFICIENT_ROTATION_LIFETIME_FLAG).is_none(),
        "the split coefficient rotation lifetime must default off"
    );
    let _environment = RawBitLengthProofEnvironment::capture();
    let configurations = [(1usize, 1usize), (2, 1), (2, 2), (2, 3)];
    let lender_modes = [false, true];
    let mut basis_states_checked = 0usize;
    let mut scalar_equivalence_checks = 0usize;
    let mut simulator_equivalence_checks = 0usize;
    let mut control_off_checks = 0usize;
    let mut lender_restore_checks = 0usize;
    let mut inverse_pair_checks = 0usize;
    let mut phase_clean_checks = 0usize;
    let mut ancilla_clean_checks = 0usize;
    let mut split_release_checks = 0usize;
    let mut split_recompute_checks = 0usize;

    for &loaned in &lender_modes {
        configure_q847_lifetime_prerequisites(loaned);
        for &(length_width, work_width) in &configurations {
            let baseline = build_split_rotation_proof_harness(
                length_width,
                work_width,
                loaned,
                SplitRotationProofRoute::Continuous,
            );
            std::env::set_var(SPLIT_COEFFICIENT_ROTATION_LIFETIME_FLAG, "1");
            let candidate = build_split_rotation_proof_harness(
                length_width,
                work_width,
                loaned,
                SplitRotationProofRoute::Configured,
            );
            assert_eq!(candidate.trace.raw_rotation_split_releases, 1);
            assert_eq!(candidate.trace.raw_rotation_split_recomputes, 1);
            assert_eq!(candidate.trace.raw_rotation_lanes_released, length_width);
            split_release_checks += 1;
            split_recompute_checks += 1;
            assert_eq!(baseline.data_ids, candidate.data_ids);
            assert_eq!(baseline.external_mask, candidate.external_mask);
            let (simulator_cases, phases, ancillas) = verify_q847_simulator_equivalence(
                b"q847-split-coefficient-rotation",
                &baseline.builder,
                &candidate.builder,
                &baseline.data_ids,
                baseline.external_mask,
            );
            simulator_equivalence_checks += simulator_cases;
            phase_clean_checks += phases;
            ancilla_clean_checks += ancillas;

            for value in 0..(1usize << baseline.data_ids.len()) {
                let input = q847_basis_input(&baseline.data_ids, value);
                let baseline_output = apply_scalar(&baseline.builder.ops, input);
                let candidate_output = apply_scalar(&candidate.builder.ops, input);
                assert_eq!(candidate_output, baseline_output);
                assert_eq!(candidate_output & !candidate.external_mask, 0);
                assert_eq!(candidate_output & candidate.lender_mask, 0);
                assert_eq!(
                    apply_scalar(&candidate.builder.ops, candidate_output),
                    input
                );
                if input & candidate.control_mask == 0 {
                    assert_eq!(candidate_output, input);
                    control_off_checks += 1;
                }
                lender_restore_checks += usize::from(loaned);
                basis_states_checked += 1;
                scalar_equivalence_checks += 1;
                inverse_pair_checks += 1;
            }
            std::env::remove_var(SPLIT_COEFFICIENT_ROTATION_LIFETIME_FLAG);
        }
    }

    for &loaned in &lender_modes {
        configure_q847_lifetime_prerequisites(loaned);
        let configured_default =
            build_split_rotation_proof_harness(2, 3, loaned, SplitRotationProofRoute::Configured);
        let direct_default =
            build_split_rotation_proof_harness(2, 3, loaned, SplitRotationProofRoute::Continuous);
        assert_q847_default_stream_identity(&configured_default.builder, &direct_default.builder);
    }

    configure_q847_lifetime_prerequisites(true);
    let baseline_local = q847_lifetime_local_resources(
        &build_split_rotation_proof_harness(
            REFERENCE_LENGTH_WIDTH,
            259,
            true,
            SplitRotationProofRoute::Continuous,
        )
        .builder,
    );
    let candidate_local = q847_lifetime_local_resources(
        &build_split_rotation_proof_harness(
            REFERENCE_LENGTH_WIDTH,
            259,
            true,
            SplitRotationProofRoute::Split,
        )
        .builder,
    );
    assert_eq!(baseline_local.active_qubits, candidate_local.active_qubits);
    assert_eq!(baseline_local.peak_qubits - candidate_local.peak_qubits, 9);
    assert_eq!(
        candidate_local.emitted_toffoli - baseline_local.emitted_toffoli,
        36
    );
    assert_eq!(
        candidate_local.emitted_ops - baseline_local.emitted_ops,
        135
    );
    let whole_point_add_invocations = 8 * REFERENCE_STEPS;
    let whole_point_add_ops_delta =
        whole_point_add_invocations * (candidate_local.emitted_ops - baseline_local.emitted_ops);
    let whole_point_add_toffoli_delta = whole_point_add_invocations
        * (candidate_local.emitted_toffoli - baseline_local.emitted_toffoli);
    assert_eq!(whole_point_add_toffoli_delta, 425_952);

    SplitCoefficientRotationLifetimeProofReport {
        configurations_checked: configurations.len() * lender_modes.len(),
        lender_modes_checked: lender_modes.len(),
        basis_states_checked,
        scalar_equivalence_checks,
        simulator_equivalence_checks,
        control_off_checks,
        lender_restore_checks,
        inverse_pair_checks,
        phase_clean_checks,
        ancilla_clean_checks,
        default_stream_identity_checks: lender_modes.len(),
        split_release_checks,
        split_recompute_checks,
        reference_rotation_lanes_released: REFERENCE_LENGTH_WIDTH,
        whole_point_add_invocations,
        whole_point_add_ops_delta: whole_point_add_ops_delta as i64,
        whole_point_add_toffoli_delta: whole_point_add_toffoli_delta as i64,
        baseline_local,
        candidate_local,
        local_qubit_delta: candidate_local.peak_qubits as i64 - baseline_local.peak_qubits as i64,
        local_ops_delta: candidate_local.emitted_ops as i64 - baseline_local.emitted_ops as i64,
        local_toffoli_delta: candidate_local.emitted_toffoli as i64
            - baseline_local.emitted_toffoli as i64,
    }
}

/// Proof-local copy of the pre-in-place less-than route from parent c562c6c7.
/// It deliberately materializes `l_t` so the differential proof does not send
/// both sides through the production in-place cursor implementation.
#[allow(clippy::too_many_arguments)]
fn legacy_toggle_coefficient_less_than_copied_cursor_for_proof(
    circ: &mut Circuit,
    control: &QReg,
    work1: &[QReg],
    work2: &[QReg],
    l_t: &[QReg],
    l_s: &[QReg],
    target_length: &[QReg],
    target: &QReg,
    compare_scratch: &[&QReg],
) {
    assert_coefficient_less_than_lane_reuse_preconditions(
        control,
        work1,
        work2,
        l_t,
        l_s,
        target_length,
        target,
        compare_scratch,
    );
    let above_t = circ.alloc_qreg("proof.legacy-coeff-compare.above-t");
    toggle_coefficient_length_above_boundary(
        circ,
        target_length,
        l_t,
        l_s,
        &above_t,
        compare_scratch,
    );
    record_coefficient_less_than_lifetime_boundary(circ, true);

    let cursor = circ.alloc_qreg_bits("proof.legacy-coeff-compare.cursor", l_t.len());
    for (source, destination) in l_t.iter().zip(&cursor) {
        circ.cx(source, destination);
    }
    let cursor_scratch = circ.alloc_qreg_bits(
        "proof.legacy-coeff-compare.cursor-scratch",
        l_t.len().saturating_sub(1),
    );
    decrement_mod_2n(circ, &cursor, &cursor_scratch);
    let carry = compare_scratch[0];
    let active = compare_scratch[1];
    let tmp = compare_scratch[2];

    for index in 0..work1.len() {
        toggle_control_and_nonnegative(circ, control, &cursor, active);
        circ.ccx(active, &work1[index], &work2[index]);
        circ.ccx(active, carry, &work1[index]);
        multi_controlled_x_vchain(
            circ,
            &[active, &work1[index], &work2[index]],
            carry,
            std::slice::from_ref(tmp),
        );
        toggle_control_and_nonnegative(circ, control, &cursor, active);
        if index + 1 != work1.len() {
            decrement_mod_2n(circ, &cursor, &cursor_scratch);
        }
    }

    circ.x(&above_t);
    circ.ccx(carry, &above_t, target);
    circ.x(&above_t);

    for index in (0..work1.len()).rev() {
        if index + 1 != work1.len() {
            increment_mod_2n(circ, &cursor, &cursor_scratch);
        }
        toggle_control_and_nonnegative(circ, control, &cursor, active);
        multi_controlled_x_vchain(
            circ,
            &[active, &work1[index], &work2[index]],
            carry,
            std::slice::from_ref(tmp),
        );
        circ.ccx(active, carry, &work1[index]);
        circ.ccx(active, &work1[index], &work2[index]);
        toggle_control_and_nonnegative(circ, control, &cursor, active);
    }
    increment_mod_2n(circ, &cursor, &cursor_scratch);
    free_clean(circ, cursor_scratch);
    for (source, destination) in l_t.iter().zip(&cursor) {
        circ.cx(source, destination);
    }
    free_clean(circ, cursor);

    record_coefficient_less_than_lifetime_boundary(circ, false);
    toggle_coefficient_length_above_boundary(
        circ,
        target_length,
        l_t,
        l_s,
        &above_t,
        compare_scratch,
    );
    circ.zero_and_free(above_t);
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum CoefficientLessThanProofRoute {
    Configured,
    Owned,
    Borrowed,
    LegacyCopiedBorrowed,
}

struct CoefficientLessThanProofHarness {
    builder: B,
    data_ids: Vec<u32>,
    external_mask: u64,
    control_mask: u64,
    length_mask: u64,
    target_mask: u64,
    lender_mask: u64,
    trace: CoefficientLessThanLifetimeTrace,
}

fn build_coefficient_less_than_proof_harness(
    length_width: usize,
    work_width: usize,
    route: CoefficientLessThanProofRoute,
) -> CoefficientLessThanProofHarness {
    assert!(length_width > 0);
    assert!(work_width > 0);
    let mut circ = Circuit::new();
    let control = circ.alloc_qreg("q847.less-than.control");
    let work1 = circ.alloc_qreg_bits("q847.less-than.work1", work_width);
    let work2 = circ.alloc_qreg_bits("q847.less-than.work2", work_width);
    let l_t = circ.alloc_qreg_bits("q847.less-than.l-t", length_width);
    let l_s = circ.alloc_qreg_bits("q847.less-than.l-s", length_width);
    let target_length = circ.alloc_qreg_bits("q847.less-than.target-length", length_width);
    let target = circ.alloc_qreg("q847.less-than.target");
    let chain = circ.alloc_qreg_bits("q847.less-than.chain", 2);
    let add_only = circ.alloc_qreg("q847.less-than.add-only");
    let compare_scratch = vec![&chain[0], &chain[1], &add_only];
    let data_ids: Vec<u32> = std::iter::once(&control)
        .chain(&work1)
        .chain(&work2)
        .chain(&l_t)
        .chain(&l_s)
        .chain(&target_length)
        .chain(std::iter::once(&target))
        .map(QReg::id)
        .collect();
    let external_mask = qreg_mask(
        std::iter::once(&control)
            .chain(&work1)
            .chain(&work2)
            .chain(&l_t)
            .chain(&l_s)
            .chain(&target_length)
            .chain(std::iter::once(&target))
            .chain(&chain)
            .chain(std::iter::once(&add_only)),
    );
    let control_mask = 1u64 << control.id();
    let length_mask = qreg_mask(&l_t);
    let target_mask = qreg_mask(std::iter::once(&target));
    let lender_mask = qreg_mask(&chain) | qreg_mask(std::iter::once(&add_only));

    begin_coefficient_less_than_lifetime_trace();
    match route {
        CoefficientLessThanProofRoute::Configured => toggle_coefficient_less_than(
            &mut circ,
            &control,
            &work1,
            &work2,
            &l_t,
            &l_s,
            &target_length,
            &target,
            &compare_scratch,
        ),
        CoefficientLessThanProofRoute::Owned => toggle_coefficient_less_than_with_lane_route(
            &mut circ,
            &control,
            &work1,
            &work2,
            &l_t,
            &l_s,
            &target_length,
            &target,
            &compare_scratch,
            false,
        ),
        CoefficientLessThanProofRoute::Borrowed => toggle_coefficient_less_than_with_lane_route(
            &mut circ,
            &control,
            &work1,
            &work2,
            &l_t,
            &l_s,
            &target_length,
            &target,
            &compare_scratch,
            true,
        ),
        CoefficientLessThanProofRoute::LegacyCopiedBorrowed => {
            legacy_toggle_coefficient_less_than_copied_cursor_for_proof(
                &mut circ,
                &control,
                &work1,
                &work2,
                &l_t,
                &l_s,
                &target_length,
                &target,
                &compare_scratch,
            )
        }
    }
    let trace = finish_coefficient_less_than_lifetime_trace();
    drop(compare_scratch);
    CoefficientLessThanProofHarness {
        builder: circ.into_builder(),
        data_ids,
        external_mask,
        control_mask,
        length_mask,
        target_mask,
        lender_mask,
        trace,
    }
}

fn coefficient_less_than_lane_alias_rejections() -> usize {
    use std::panic::{catch_unwind, AssertUnwindSafe};

    fn build_rejection(kind: usize) {
        let mut circ = Circuit::new();
        let control = circ.alloc_qreg("q847.alias.control");
        let work1 = circ.alloc_qreg_bits("q847.alias.work1", 1);
        let work2 = circ.alloc_qreg_bits("q847.alias.work2", 1);
        let l_t = circ.alloc_qreg_bits("q847.alias.l-t", 2);
        let l_s = circ.alloc_qreg_bits("q847.alias.l-s", 2);
        let target_length = circ.alloc_qreg_bits("q847.alias.target-length", 2);
        let target = circ.alloc_qreg("q847.alias.target");
        let lenders = circ.alloc_qreg_bits("q847.alias.lenders", 3);
        let scratch: Vec<&QReg> = match kind {
            0 => vec![&lenders[0], &lenders[1]],
            1 => vec![&lenders[0], &lenders[0], &lenders[2]],
            2 => vec![&control, &lenders[1], &lenders[2]],
            3 => vec![&work1[0], &lenders[1], &lenders[2]],
            4 => vec![&work2[0], &lenders[1], &lenders[2]],
            5 => vec![&l_t[0], &lenders[1], &lenders[2]],
            6 => vec![&l_s[0], &lenders[1], &lenders[2]],
            7 => vec![&target_length[0], &lenders[1], &lenders[2]],
            8 => vec![&target, &lenders[1], &lenders[2]],
            _ => unreachable!(),
        };
        toggle_coefficient_less_than_with_lane_route(
            &mut circ,
            &control,
            &work1,
            &work2,
            &l_t,
            &l_s,
            &target_length,
            &target,
            &scratch,
            true,
        );
    }

    let previous_hook = std::panic::take_hook();
    std::panic::set_hook(Box::new(|_| {}));
    let mut rejected = 0usize;
    for kind in 0..9 {
        let result = catch_unwind(AssertUnwindSafe(|| build_rejection(kind)));
        assert!(
            result.is_err(),
            "coefficient less-than alias rejection {kind} unexpectedly passed"
        );
        rejected += 1;
    }
    std::panic::set_hook(previous_hook);
    rejected
}

/// Exhaustively verify the sequential alias of `(chain[0], chain[1],
/// add_only)` from the restored boundary comparator into the less-than carry
/// loop, including both cut points and explicit disjointness failures.
#[doc(hidden)]
#[must_use]
pub fn exhaustive_coefficient_less_than_lane_reuse_check() -> CoefficientLessThanLaneReuseProofReport
{
    assert!(
        std::env::var_os(COEFFICIENT_LESS_THAN_LANE_REUSE_FLAG).is_none(),
        "the coefficient less-than lane reuse must default off"
    );
    let _environment = RawBitLengthProofEnvironment::capture();
    configure_q847_lifetime_prerequisites(true);
    let configurations = [(1usize, 1usize), (2, 1), (2, 2), (2, 3)];
    let mut basis_states_checked = 0usize;
    let mut scalar_equivalence_checks = 0usize;
    let mut simulator_equivalence_checks = 0usize;
    let mut control_off_checks = 0usize;
    let mut boundary_restore_checks = 0usize;
    let mut lender_restore_checks = 0usize;
    let mut inverse_pair_checks = 0usize;
    let mut phase_clean_checks = 0usize;
    let mut ancilla_clean_checks = 0usize;

    for &(length_width, work_width) in &configurations {
        std::env::remove_var(COEFFICIENT_LESS_THAN_LANE_REUSE_FLAG);
        let baseline = build_coefficient_less_than_proof_harness(
            length_width,
            work_width,
            CoefficientLessThanProofRoute::Owned,
        );
        std::env::set_var(COEFFICIENT_LESS_THAN_LANE_REUSE_FLAG, "1");
        let candidate = build_coefficient_less_than_proof_harness(
            length_width,
            work_width,
            CoefficientLessThanProofRoute::Configured,
        );
        assert!(candidate.trace.after_first_boundary > 0);
        assert!(candidate.trace.before_second_boundary > candidate.trace.after_first_boundary);
        assert_eq!(baseline.data_ids, candidate.data_ids);
        assert_eq!(baseline.external_mask, candidate.external_mask);
        let (simulator_cases, phases, ancillas) = verify_q847_simulator_equivalence(
            b"q847-coefficient-less-than-lanes",
            &baseline.builder,
            &candidate.builder,
            &baseline.data_ids,
            baseline.external_mask,
        );
        simulator_equivalence_checks += simulator_cases;
        phase_clean_checks += phases;
        ancilla_clean_checks += ancillas;

        for value in 0..(1usize << baseline.data_ids.len()) {
            let input = q847_basis_input(&baseline.data_ids, value);
            let after_boundary = apply_scalar(
                &candidate.builder.ops[..candidate.trace.after_first_boundary],
                input,
            );
            let before_second_boundary = apply_scalar(
                &candidate.builder.ops[..candidate.trace.before_second_boundary],
                input,
            );
            assert_eq!(after_boundary & candidate.lender_mask, 0);
            assert_eq!(before_second_boundary & candidate.lender_mask, 0);
            assert_eq!(
                after_boundary & candidate.length_mask,
                input & candidate.length_mask
            );
            assert_eq!(
                before_second_boundary & candidate.length_mask,
                input & candidate.length_mask
            );
            boundary_restore_checks += 2;

            let baseline_output = apply_scalar(&baseline.builder.ops, input);
            let candidate_output = apply_scalar(&candidate.builder.ops, input);
            assert_eq!(candidate_output, baseline_output);
            assert_eq!(candidate_output & !candidate.external_mask, 0);
            assert_eq!(candidate_output & candidate.lender_mask, 0);
            assert_eq!((candidate_output ^ input) & !candidate.target_mask, 0);
            assert_eq!(
                apply_scalar(&candidate.builder.ops, candidate_output),
                input
            );
            if input & candidate.control_mask == 0 {
                assert_eq!(candidate_output, input);
                control_off_checks += 1;
            }
            basis_states_checked += 1;
            scalar_equivalence_checks += 1;
            lender_restore_checks += 1;
            inverse_pair_checks += 1;
        }
        std::env::remove_var(COEFFICIENT_LESS_THAN_LANE_REUSE_FLAG);
    }

    let configured_default =
        build_coefficient_less_than_proof_harness(2, 3, CoefficientLessThanProofRoute::Configured);
    let direct_default =
        build_coefficient_less_than_proof_harness(2, 3, CoefficientLessThanProofRoute::Owned);
    assert_q847_default_stream_identity(&configured_default.builder, &direct_default.builder);
    let alias_rejections = coefficient_less_than_lane_alias_rejections();

    let baseline_local = q847_lifetime_local_resources(
        &build_coefficient_less_than_proof_harness(
            REFERENCE_LENGTH_WIDTH,
            259,
            CoefficientLessThanProofRoute::Owned,
        )
        .builder,
    );
    let candidate_local = q847_lifetime_local_resources(
        &build_coefficient_less_than_proof_harness(
            REFERENCE_LENGTH_WIDTH,
            259,
            CoefficientLessThanProofRoute::Borrowed,
        )
        .builder,
    );
    assert_eq!(baseline_local.active_qubits, candidate_local.active_qubits);
    assert_eq!(baseline_local.peak_qubits - candidate_local.peak_qubits, 1);
    assert_eq!(baseline_local.emitted_ops - candidate_local.emitted_ops, 3);
    assert_eq!(
        baseline_local.emitted_toffoli,
        candidate_local.emitted_toffoli
    );
    let whole_point_add_invocations = 8 * REFERENCE_STEPS;

    CoefficientLessThanLaneReuseProofReport {
        configurations_checked: configurations.len(),
        basis_states_checked,
        scalar_equivalence_checks,
        simulator_equivalence_checks,
        control_off_checks,
        boundary_restore_checks,
        lender_restore_checks,
        inverse_pair_checks,
        phase_clean_checks,
        ancilla_clean_checks,
        default_stream_identity_checks: 1,
        alias_rejections,
        caller_lanes_reused: 3,
        reset_ops_removed_per_invocation: 3,
        whole_point_add_invocations,
        whole_point_add_ops_delta: -((whole_point_add_invocations * 3) as i64),
        baseline_local,
        candidate_local,
        local_qubit_delta: candidate_local.peak_qubits as i64 - baseline_local.peak_qubits as i64,
        local_ops_delta: candidate_local.emitted_ops as i64 - baseline_local.emitted_ops as i64,
        local_toffoli_delta: candidate_local.emitted_toffoli as i64
            - baseline_local.emitted_toffoli as i64,
    }
}

/// Proof-local copy of the pre-in-place coefficient-add route from parent
/// c562c6c7. The copied cursor is retained solely as an independent baseline.
#[allow(clippy::too_many_arguments)]
fn legacy_coefficient_add_data_only_copied_cursor_for_proof(
    circ: &mut Circuit,
    control: &QReg,
    work1: &[QReg],
    work2: &[QReg],
    l_t: &[QReg],
    inverse: bool,
    lenders: &[&QReg],
) {
    assert_clean_chain_coefficient_add_lender_preconditions(control, work1, work2, l_t, lenders);
    record_coefficient_add_lender_entry(circ, lenders);
    let cursor = circ.alloc_qreg_bits("proof.legacy-coeff-add.cursor", l_t.len());
    for (source, destination) in l_t.iter().zip(&cursor) {
        circ.cx(source, destination);
    }
    let cursor_scratch = circ.alloc_qreg_bits(
        "proof.legacy-coeff-add.cursor-scratch",
        l_t.len().saturating_sub(1),
    );
    let carry = circ.alloc_qreg("proof.legacy-coeff-add.carry");
    let active = lenders[0];
    let tmp = lenders[1];

    if inverse {
        for index in 0..work1.len() {
            if index != 0 {
                decrement_mod_2n(circ, &cursor, &cursor_scratch);
            }
            toggle_control_and_nonnegative(circ, control, &cursor, active);
            circ.ccx(active, &work1[index], &work2[index]);
            circ.ccx(active, &carry, &work1[index]);
            multi_controlled_x_vchain(
                circ,
                &[active, &work1[index], &work2[index]],
                &carry,
                std::slice::from_ref(tmp),
            );
            toggle_control_and_nonnegative(circ, control, &cursor, active);
        }
        for index in (0..work1.len()).rev() {
            if index + 1 != work1.len() {
                increment_mod_2n(circ, &cursor, &cursor_scratch);
            }
            toggle_control_and_nonnegative(circ, control, &cursor, active);
            multi_controlled_x_vchain(
                circ,
                &[active, &work1[index], &work2[index]],
                &carry,
                std::slice::from_ref(tmp),
            );
            circ.ccx(active, &carry, &work1[index]);
            circ.ccx(active, &carry, &work2[index]);
            toggle_control_and_nonnegative(circ, control, &cursor, active);
        }
    } else {
        for index in 0..work1.len() {
            toggle_control_and_nonnegative(circ, control, &cursor, active);
            circ.ccx(active, &carry, &work2[index]);
            circ.ccx(active, &carry, &work1[index]);
            multi_controlled_x_vchain(
                circ,
                &[active, &work1[index], &work2[index]],
                &carry,
                std::slice::from_ref(tmp),
            );
            toggle_control_and_nonnegative(circ, control, &cursor, active);
            if index + 1 != work1.len() {
                decrement_mod_2n(circ, &cursor, &cursor_scratch);
            }
        }
        for index in (0..work1.len()).rev() {
            toggle_control_and_nonnegative(circ, control, &cursor, active);
            multi_controlled_x_vchain(
                circ,
                &[active, &work1[index], &work2[index]],
                &carry,
                std::slice::from_ref(tmp),
            );
            circ.ccx(active, &carry, &work1[index]);
            circ.ccx(active, &work1[index], &work2[index]);
            toggle_control_and_nonnegative(circ, control, &cursor, active);
            if index != 0 {
                increment_mod_2n(circ, &cursor, &cursor_scratch);
            }
        }
    }

    record_coefficient_add_lender_restore(circ);
    circ.zero_and_free(carry);
    free_clean(circ, cursor_scratch);
    for (source, destination) in l_t.iter().zip(&cursor) {
        circ.cx(source, destination);
    }
    free_clean(circ, cursor);
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum CoefficientAddLenderProofRoute {
    Configured,
    Owned,
    Borrowed,
    LegacyCopiedBorrowed,
}

struct CoefficientAddLenderProofHarness {
    builder: B,
    data_ids: Vec<u32>,
    external_mask: u64,
    control_mask: u64,
    length_mask: u64,
    lender_mask: u64,
    trace: CoefficientAddLenderTrace,
}

fn build_coefficient_add_lender_proof_harness(
    length_width: usize,
    work_width: usize,
    inverse: bool,
    route: CoefficientAddLenderProofRoute,
) -> CoefficientAddLenderProofHarness {
    assert!(length_width > 0);
    assert!(work_width > 0);
    let mut circ = Circuit::new();
    let control = circ.alloc_qreg("q847.coeff-add.control");
    let work1 = circ.alloc_qreg_bits("q847.coeff-add.work1", work_width);
    let work2 = circ.alloc_qreg_bits("q847.coeff-add.work2", work_width);
    let l_t = circ.alloc_qreg_bits("q847.coeff-add.l-t", length_width);
    let lenders = circ.alloc_qreg_bits("q847.coeff-add.clean-chain", 2);
    let lender_refs = [&lenders[0], &lenders[1]];
    let data_ids: Vec<u32> = std::iter::once(&control)
        .chain(&work1)
        .chain(&work2)
        .chain(&l_t)
        .map(QReg::id)
        .collect();
    let external_mask = qreg_mask(
        std::iter::once(&control)
            .chain(&work1)
            .chain(&work2)
            .chain(&l_t)
            .chain(&lenders),
    );
    let control_mask = 1u64 << control.id();
    let length_mask = qreg_mask(&l_t);
    let lender_mask = qreg_mask(&lenders);

    begin_coefficient_add_lender_trace();
    match route {
        CoefficientAddLenderProofRoute::Configured => coefficient_add_data_only(
            &mut circ,
            &control,
            &work1,
            &work2,
            &l_t,
            inverse,
            &lender_refs,
        ),
        CoefficientAddLenderProofRoute::Owned => coefficient_add_data_only_with_lane_route(
            &mut circ,
            &control,
            &work1,
            &work2,
            &l_t,
            inverse,
            &lender_refs,
            false,
        ),
        CoefficientAddLenderProofRoute::Borrowed => coefficient_add_data_only_with_lane_route(
            &mut circ,
            &control,
            &work1,
            &work2,
            &l_t,
            inverse,
            &lender_refs,
            true,
        ),
        CoefficientAddLenderProofRoute::LegacyCopiedBorrowed => {
            legacy_coefficient_add_data_only_copied_cursor_for_proof(
                &mut circ,
                &control,
                &work1,
                &work2,
                &l_t,
                inverse,
                &lender_refs,
            )
        }
    }
    let trace = finish_coefficient_add_lender_trace();

    CoefficientAddLenderProofHarness {
        builder: circ.into_builder(),
        data_ids,
        external_mask,
        control_mask,
        length_mask,
        lender_mask,
        trace,
    }
}

fn clean_chain_coefficient_add_distinctness_rejections() -> usize {
    use std::panic::{catch_unwind, AssertUnwindSafe};

    fn build_rejection(kind: usize) {
        let mut circ = Circuit::new();
        let control = circ.alloc_qreg("q847.coeff-add-alias.control");
        let work1 = circ.alloc_qreg_bits("q847.coeff-add-alias.work1", 1);
        let work2 = circ.alloc_qreg_bits("q847.coeff-add-alias.work2", 1);
        let l_t = circ.alloc_qreg_bits("q847.coeff-add-alias.l-t", 2);
        let lenders = circ.alloc_qreg_bits("q847.coeff-add-alias.lenders", 3);
        let scratch: Vec<&QReg> = match kind {
            0 => vec![&lenders[0]],
            1 => vec![&lenders[0], &lenders[0]],
            2 => vec![&control, &lenders[1]],
            3 => vec![&work1[0], &lenders[1]],
            4 => vec![&work2[0], &lenders[1]],
            5 => vec![&l_t[0], &lenders[1]],
            6 => vec![&lenders[0], &lenders[1], &lenders[2]],
            _ => unreachable!(),
        };
        coefficient_add_data_only_with_lane_route(
            &mut circ, &control, &work1, &work2, &l_t, false, &scratch, true,
        );
    }

    let previous_hook = std::panic::take_hook();
    std::panic::set_hook(Box::new(|_| {}));
    let mut rejected = 0usize;
    for kind in 0..7 {
        let result = catch_unwind(AssertUnwindSafe(|| build_rejection(kind)));
        assert!(
            result.is_err(),
            "clean-chain coefficient-add rejection {kind} unexpectedly passed"
        );
        rejected += 1;
    }
    std::panic::set_hook(previous_hook);
    rejected
}

struct CoefficientPhaseBlockProofHarness {
    builder: B,
    data_ids: Vec<u32>,
    external_mask: u64,
    length_mask: u64,
    trace: CoefficientAddLenderTrace,
}

fn build_coefficient_phase_block_lender_proof_harness(
    length_width: usize,
    work_width: usize,
    inverse: bool,
) -> CoefficientPhaseBlockProofHarness {
    assert!(length_width > 0);
    assert!(work_width > 0);
    let mut circ = Circuit::new();
    let phase1 = circ.alloc_qreg("q847.coeff-block.phase1");
    let phase2 = circ.alloc_qreg("q847.coeff-block.phase2");
    let sign = circ.alloc_qreg("q847.coeff-block.sign");
    let work1 = circ.alloc_qreg_bits("q847.coeff-block.work1", work_width);
    let work2 = circ.alloc_qreg_bits("q847.coeff-block.work2", work_width);
    let l_t = circ.alloc_qreg_bits("q847.coeff-block.l-t", length_width);
    let l_t_prime = circ.alloc_qreg_bits("q847.coeff-block.l-t-prime", length_width);
    let l_s = circ.alloc_qreg_bits("q847.coeff-block.l-s", length_width);
    let l_r_prime = circ.alloc_qreg_bits("q847.coeff-block.l-r-prime", length_width);
    let data_ids: Vec<u32> = [&phase1, &phase2, &sign]
        .into_iter()
        .chain(&work1)
        .chain(&work2)
        .chain(&l_t)
        .chain(&l_t_prime)
        .chain(&l_s)
        .chain(&l_r_prime)
        .map(QReg::id)
        .collect();
    let external_mask = qreg_mask(
        [&phase1, &phase2, &sign]
            .into_iter()
            .chain(&work1)
            .chain(&work2)
            .chain(&l_t)
            .chain(&l_t_prime)
            .chain(&l_s)
            .chain(&l_r_prime),
    );
    let length_mask = qreg_mask(&l_t);

    begin_coefficient_add_lender_trace();
    coefficient_phase_block(
        &mut circ, &phase1, &phase2, &sign, &work1, &work2, &work2, &l_t, &l_t_prime, &l_s,
        &l_r_prime, inverse,
    );
    let trace = finish_coefficient_add_lender_trace();

    CoefficientPhaseBlockProofHarness {
        builder: circ.into_builder(),
        data_ids,
        external_mask,
        length_mask,
        trace,
    }
}

fn configure_q849_lifetime_cuts_for_proof() {
    configure_q847_lifetime_prerequisites(true);
    std::env::set_var(PROMISED_LQ_SWAP_BORROW_FLAG, "1");
    std::env::set_var(SPLIT_COEFFICIENT_ROTATION_LIFETIME_FLAG, "1");
    std::env::set_var(COEFFICIENT_LESS_THAN_LANE_REUSE_FLAG, "1");
}

/// Exhaustively prove the default-off clean-chain lender route for the
/// coefficient data update and its reduced-width `coefficient_phase_block`
/// composition. Phase checks cover the exact lender entry/restore window; the
/// pre-existing reset stream outside that window is compared for data and
/// ancilla equivalence but is not promoted into a trusted phase claim. No
/// source-baked route is exercised by this proof.
#[doc(hidden)]
#[must_use]
pub fn exhaustive_clean_chain_coefficient_add_lender_check(
) -> CleanChainCoefficientAddLenderProofReport {
    assert!(
        std::env::var_os(CLEAN_CHAIN_COEFFICIENT_ADD_LENDER_FLAG).is_none(),
        "the clean-chain coefficient-add lender must default off"
    );
    let _environment = RawBitLengthProofEnvironment::capture();
    configure_q847_lifetime_prerequisites(true);
    let configurations = [(1usize, 1usize), (2, 1), (2, 2), (2, 3), (3, 4)];
    let mut basis_states_checked = 0usize;
    let mut scalar_equivalence_checks = 0usize;
    let mut simulator_equivalence_checks = 0usize;
    let mut control_off_identity_checks = 0usize;
    let mut length_preservation_checks = 0usize;
    let mut lender_clean_entry_checks = 0usize;
    let mut lender_restore_checks = 0usize;
    let mut roundtrip_checks = 0usize;
    let mut phase_clean_checks = 0usize;
    let mut ancilla_clean_checks = 0usize;

    for &(length_width, work_width) in &configurations {
        std::env::remove_var(CLEAN_CHAIN_COEFFICIENT_ADD_LENDER_FLAG);
        let baseline_forward = build_coefficient_add_lender_proof_harness(
            length_width,
            work_width,
            false,
            CoefficientAddLenderProofRoute::Owned,
        );
        let baseline_inverse = build_coefficient_add_lender_proof_harness(
            length_width,
            work_width,
            true,
            CoefficientAddLenderProofRoute::Owned,
        );
        std::env::set_var(CLEAN_CHAIN_COEFFICIENT_ADD_LENDER_FLAG, "1");
        let candidate_forward = build_coefficient_add_lender_proof_harness(
            length_width,
            work_width,
            false,
            CoefficientAddLenderProofRoute::Configured,
        );
        let candidate_inverse = build_coefficient_add_lender_proof_harness(
            length_width,
            work_width,
            true,
            CoefficientAddLenderProofRoute::Configured,
        );

        for (inverse, baseline, candidate, opposite) in [
            (
                false,
                &baseline_forward,
                &candidate_forward,
                &candidate_inverse,
            ),
            (
                true,
                &baseline_inverse,
                &candidate_inverse,
                &candidate_forward,
            ),
        ] {
            assert_eq!(baseline.data_ids, candidate.data_ids);
            assert_eq!(baseline.external_mask, candidate.external_mask);
            assert_eq!(candidate.trace.calls, 1);
            assert_eq!(candidate.trace.lender_mask, candidate.lender_mask);
            assert!(candidate.trace.restore_ops_idx > candidate.trace.entry_ops_idx);
            let label: &[u8] = if inverse {
                b"q847-clean-chain-coefficient-add-inverse"
            } else {
                b"q847-clean-chain-coefficient-add-forward"
            };
            let (simulator_cases, ancillas) = verify_q847_simulator_data_and_ancilla_equivalence(
                label,
                &baseline.builder,
                &candidate.builder,
                &baseline.data_ids,
                baseline.external_mask,
            );
            simulator_equivalence_checks += simulator_cases;
            phase_clean_checks += verify_coefficient_add_lender_window_phase_clean(
                label,
                &candidate.builder,
                &candidate.data_ids,
                candidate.trace,
            );
            ancilla_clean_checks += ancillas;

            for value in 0..(1usize << baseline.data_ids.len()) {
                let input = q847_basis_input(&baseline.data_ids, value);
                let at_entry = apply_scalar(
                    &candidate.builder.ops[..candidate.trace.entry_ops_idx],
                    input,
                );
                let at_restore = apply_scalar(
                    &candidate.builder.ops[..candidate.trace.restore_ops_idx],
                    input,
                );
                assert_eq!(at_entry & candidate.lender_mask, 0);
                assert_eq!(at_restore & candidate.lender_mask, 0);
                assert_eq!(
                    at_entry & candidate.length_mask,
                    input & candidate.length_mask
                );
                assert_eq!(
                    at_restore & candidate.length_mask,
                    input & candidate.length_mask
                );

                let baseline_output = apply_scalar(&baseline.builder.ops, input);
                let candidate_output = apply_scalar(&candidate.builder.ops, input);
                assert_eq!(candidate_output, baseline_output);
                assert_eq!(candidate_output & !candidate.external_mask, 0);
                assert_eq!(
                    candidate_output & candidate.length_mask,
                    input & candidate.length_mask
                );
                assert_eq!(candidate_output & candidate.lender_mask, 0);
                assert_eq!(apply_scalar(&opposite.builder.ops, candidate_output), input);
                if input & candidate.control_mask == 0 {
                    assert_eq!(candidate_output, input);
                    control_off_identity_checks += 1;
                }

                basis_states_checked += 1;
                scalar_equivalence_checks += 1;
                length_preservation_checks += 1;
                lender_clean_entry_checks += 1;
                lender_restore_checks += 1;
                roundtrip_checks += 1;
            }
        }
    }

    std::env::remove_var(CLEAN_CHAIN_COEFFICIENT_ADD_LENDER_FLAG);
    let mut default_stream_identity_checks = 0usize;
    for inverse in [false, true] {
        let configured = build_coefficient_add_lender_proof_harness(
            2,
            3,
            inverse,
            CoefficientAddLenderProofRoute::Configured,
        );
        let direct = build_coefficient_add_lender_proof_harness(
            2,
            3,
            inverse,
            CoefficientAddLenderProofRoute::Owned,
        );
        assert_q847_default_stream_identity(&configured.builder, &direct.builder);
        default_stream_identity_checks += 1;
    }
    let distinctness_rejections = clean_chain_coefficient_add_distinctness_rejections();

    configure_q849_lifetime_cuts_for_proof();
    let composition_length_width = 1usize;
    let composition_work_width = 1usize;
    std::env::remove_var(CLEAN_CHAIN_COEFFICIENT_ADD_LENDER_FLAG);
    let composition_baseline_forward = build_coefficient_phase_block_lender_proof_harness(
        composition_length_width,
        composition_work_width,
        false,
    );
    let composition_baseline_inverse = build_coefficient_phase_block_lender_proof_harness(
        composition_length_width,
        composition_work_width,
        true,
    );
    std::env::set_var(CLEAN_CHAIN_COEFFICIENT_ADD_LENDER_FLAG, "1");
    let composition_candidate_forward = build_coefficient_phase_block_lender_proof_harness(
        composition_length_width,
        composition_work_width,
        false,
    );
    let composition_candidate_inverse = build_coefficient_phase_block_lender_proof_harness(
        composition_length_width,
        composition_work_width,
        true,
    );
    let mut composition_basis_states_checked = 0usize;
    let mut composition_equivalence_checks = 0usize;
    let mut composition_lender_clean_entry_checks = 0usize;
    let mut composition_lender_restore_checks = 0usize;
    let mut composition_phase_clean_checks = 0usize;
    let mut composition_ancilla_clean_checks = 0usize;
    for (inverse, baseline, candidate) in [
        (
            false,
            &composition_baseline_forward,
            &composition_candidate_forward,
        ),
        (
            true,
            &composition_baseline_inverse,
            &composition_candidate_inverse,
        ),
    ] {
        assert_eq!(baseline.data_ids, candidate.data_ids);
        assert_eq!(baseline.external_mask, candidate.external_mask);
        assert_eq!(candidate.trace.calls, 1);
        assert!(candidate.trace.restore_ops_idx > candidate.trace.entry_ops_idx);
        let label: &[u8] = if inverse {
            b"q847-clean-chain-coefficient-block-inverse"
        } else {
            b"q847-clean-chain-coefficient-block-forward"
        };
        let (simulator_cases, ancillas) = verify_q847_simulator_data_and_ancilla_equivalence(
            label,
            &baseline.builder,
            &candidate.builder,
            &baseline.data_ids,
            baseline.external_mask,
        );
        assert_eq!(simulator_cases, 1usize << baseline.data_ids.len());
        composition_phase_clean_checks += verify_coefficient_add_lender_window_phase_clean(
            label,
            &candidate.builder,
            &candidate.data_ids,
            candidate.trace,
        );
        composition_ancilla_clean_checks += ancillas;

        for value in 0..(1usize << baseline.data_ids.len()) {
            let input = q847_basis_input(&baseline.data_ids, value);
            let at_entry = apply_scalar(
                &candidate.builder.ops[..candidate.trace.entry_ops_idx],
                input,
            );
            let at_restore = apply_scalar(
                &candidate.builder.ops[..candidate.trace.restore_ops_idx],
                input,
            );
            assert_eq!(at_entry & candidate.trace.lender_mask, 0);
            assert_eq!(at_restore & candidate.trace.lender_mask, 0);
            let baseline_output = apply_scalar(&baseline.builder.ops, input);
            let candidate_output = apply_scalar(&candidate.builder.ops, input);
            assert_eq!(candidate_output, baseline_output);
            assert_eq!(candidate_output & !candidate.external_mask, 0);
            assert_eq!(
                at_entry & candidate.length_mask,
                input & candidate.length_mask
            );
            assert_eq!(
                at_restore & candidate.length_mask,
                input & candidate.length_mask
            );

            composition_basis_states_checked += 1;
            composition_equivalence_checks += 1;
            composition_lender_clean_entry_checks += 1;
            composition_lender_restore_checks += 1;
        }
    }

    configure_q847_lifetime_prerequisites(true);
    let baseline_local = q847_lifetime_local_resources(
        &build_coefficient_add_lender_proof_harness(
            REFERENCE_LENGTH_WIDTH,
            259,
            false,
            CoefficientAddLenderProofRoute::Owned,
        )
        .builder,
    );
    let candidate_local = q847_lifetime_local_resources(
        &build_coefficient_add_lender_proof_harness(
            REFERENCE_LENGTH_WIDTH,
            259,
            false,
            CoefficientAddLenderProofRoute::Borrowed,
        )
        .builder,
    );
    assert_eq!(baseline_local.active_qubits, candidate_local.active_qubits);
    assert_eq!(
        baseline_local.emitted_toffoli,
        candidate_local.emitted_toffoli
    );
    let whole_point_add_invocations = 4 * REFERENCE_STEPS;
    let local_ops_delta = candidate_local.emitted_ops as i64 - baseline_local.emitted_ops as i64;

    CleanChainCoefficientAddLenderProofReport {
        configurations_checked: configurations.len(),
        directions_checked: 2,
        basis_states_checked,
        scalar_equivalence_checks,
        simulator_equivalence_checks,
        control_off_identity_checks,
        length_preservation_checks,
        lender_clean_entry_checks,
        lender_restore_checks,
        roundtrip_checks,
        lender_window_phase_clean_checks: phase_clean_checks,
        ancilla_clean_checks,
        default_stream_identity_checks,
        distinctness_rejections,
        composition_basis_states_checked,
        composition_equivalence_checks,
        composition_lender_clean_entry_checks,
        composition_lender_restore_checks,
        composition_lender_window_phase_clean_checks: composition_phase_clean_checks,
        composition_ancilla_clean_checks,
        caller_lanes_reused: 2,
        reset_ops_removed_per_invocation: baseline_local.emitted_ops - candidate_local.emitted_ops,
        whole_point_add_invocations,
        whole_point_add_ops_delta: whole_point_add_invocations as i64 * local_ops_delta,
        baseline_local,
        candidate_local,
        local_qubit_delta: candidate_local.peak_qubits as i64 - baseline_local.peak_qubits as i64,
        local_ops_delta,
        local_toffoli_delta: candidate_local.emitted_toffoli as i64
            - baseline_local.emitted_toffoli as i64,
    }
}

fn assert_legacy_copied_cursor_stream_delta(legacy: &B, in_place: &B, width: usize) {
    use crate::circuit::OperationType;

    assert_eq!(legacy.counted_ops, legacy.ops.len());
    assert_eq!(in_place.counted_ops, in_place.ops.len());
    assert_eq!(legacy.counted_ops - in_place.counted_ops, 3 * width);
    for kind in 0..legacy.counted_kind_ops.len() {
        let expected = if kind == OperationType::CX as usize {
            2 * width
        } else if kind == OperationType::R as usize {
            width
        } else {
            0
        };
        assert_eq!(
            legacy.counted_kind_ops[kind] - in_place.counted_kind_ops[kind],
            expected,
            "legacy copied-cursor operation delta drift for kind {kind}"
        );
    }
}

/// Exhaustively compare the production in-place `l_t` cursor routes against
/// proof-local copies of the parent c562c6c7 materialized-cursor schedules.
/// Both directions, all reduced-width basis states, phase, ancillas, control
/// off, length restoration, and inverse pairing are checked independently of
/// the route-to-route lender proofs.
#[doc(hidden)]
#[must_use]
pub fn exhaustive_inplace_lt_cursor_legacy_differential_check(
) -> InPlaceLtCursorLegacyDifferentialProofReport {
    let _environment = RawBitLengthProofEnvironment::capture();
    configure_q847_lifetime_prerequisites(true);
    let less_than_configurations = [(1usize, 1usize), (2, 1), (2, 2), (2, 3)];
    let coefficient_add_configurations = [(1usize, 1usize), (2, 1), (2, 2), (2, 3), (3, 4)];
    let mut basis_states_checked = 0usize;
    let mut scalar_equivalence_checks = 0usize;
    let mut simulator_equivalence_checks = 0usize;
    let mut control_off_checks = 0usize;
    let mut length_restore_checks = 0usize;
    let mut inverse_pair_checks = 0usize;
    let mut phase_clean_checks = 0usize;
    let mut ancilla_clean_checks = 0usize;
    let mut legacy_streams_checked = 0usize;
    let mut copied_cursor_lanes_removed = 0usize;

    for &(length_width, work_width) in &less_than_configurations {
        let legacy = build_coefficient_less_than_proof_harness(
            length_width,
            work_width,
            CoefficientLessThanProofRoute::LegacyCopiedBorrowed,
        );
        let in_place = build_coefficient_less_than_proof_harness(
            length_width,
            work_width,
            CoefficientLessThanProofRoute::Borrowed,
        );
        assert_eq!(legacy.data_ids, in_place.data_ids);
        assert_eq!(legacy.external_mask, in_place.external_mask);
        assert_legacy_copied_cursor_stream_delta(&legacy.builder, &in_place.builder, length_width);
        legacy_streams_checked += 1;
        copied_cursor_lanes_removed += length_width;

        let (simulator_cases, phases, ancillas) = verify_q847_simulator_equivalence(
            b"q847-inplace-lt-cursor-legacy-less-than",
            &legacy.builder,
            &in_place.builder,
            &legacy.data_ids,
            legacy.external_mask,
        );
        simulator_equivalence_checks += simulator_cases;
        phase_clean_checks += phases;
        ancilla_clean_checks += ancillas;

        for value in 0..(1usize << legacy.data_ids.len()) {
            let input = q847_basis_input(&legacy.data_ids, value);
            let legacy_output = apply_scalar(&legacy.builder.ops, input);
            let in_place_output = apply_scalar(&in_place.builder.ops, input);
            assert_eq!(in_place_output, legacy_output);
            assert_eq!(legacy_output & !legacy.external_mask, 0);
            assert_eq!(in_place_output & !in_place.external_mask, 0);
            assert_eq!(
                legacy_output & legacy.length_mask,
                input & legacy.length_mask
            );
            assert_eq!(
                in_place_output & in_place.length_mask,
                input & in_place.length_mask
            );
            assert_eq!((in_place_output ^ input) & !in_place.target_mask, 0);
            assert_eq!(apply_scalar(&legacy.builder.ops, legacy_output), input);
            assert_eq!(apply_scalar(&in_place.builder.ops, in_place_output), input);
            if input & in_place.control_mask == 0 {
                assert_eq!(legacy_output, input);
                assert_eq!(in_place_output, input);
                control_off_checks += 1;
            }
            basis_states_checked += 1;
            scalar_equivalence_checks += 1;
            length_restore_checks += 1;
            inverse_pair_checks += 2;
        }
    }

    for &(length_width, work_width) in &coefficient_add_configurations {
        let legacy_forward = build_coefficient_add_lender_proof_harness(
            length_width,
            work_width,
            false,
            CoefficientAddLenderProofRoute::LegacyCopiedBorrowed,
        );
        let legacy_inverse = build_coefficient_add_lender_proof_harness(
            length_width,
            work_width,
            true,
            CoefficientAddLenderProofRoute::LegacyCopiedBorrowed,
        );
        let in_place_forward = build_coefficient_add_lender_proof_harness(
            length_width,
            work_width,
            false,
            CoefficientAddLenderProofRoute::Borrowed,
        );
        let in_place_inverse = build_coefficient_add_lender_proof_harness(
            length_width,
            work_width,
            true,
            CoefficientAddLenderProofRoute::Borrowed,
        );

        for (inverse, legacy, in_place, legacy_opposite, in_place_opposite) in [
            (
                false,
                &legacy_forward,
                &in_place_forward,
                &legacy_inverse,
                &in_place_inverse,
            ),
            (
                true,
                &legacy_inverse,
                &in_place_inverse,
                &legacy_forward,
                &in_place_forward,
            ),
        ] {
            assert_eq!(legacy.data_ids, in_place.data_ids);
            assert_eq!(legacy.external_mask, in_place.external_mask);
            assert_legacy_copied_cursor_stream_delta(
                &legacy.builder,
                &in_place.builder,
                length_width,
            );
            legacy_streams_checked += 1;
            copied_cursor_lanes_removed += length_width;
            let label: &[u8] = if inverse {
                b"q847-inplace-lt-cursor-legacy-add-inverse"
            } else {
                b"q847-inplace-lt-cursor-legacy-add-forward"
            };
            let (simulator_cases, phases, ancillas) = verify_q847_simulator_equivalence(
                label,
                &legacy.builder,
                &in_place.builder,
                &legacy.data_ids,
                legacy.external_mask,
            );
            simulator_equivalence_checks += simulator_cases;
            phase_clean_checks += phases;
            ancilla_clean_checks += ancillas;

            for value in 0..(1usize << legacy.data_ids.len()) {
                let input = q847_basis_input(&legacy.data_ids, value);
                let legacy_output = apply_scalar(&legacy.builder.ops, input);
                let in_place_output = apply_scalar(&in_place.builder.ops, input);
                assert_eq!(in_place_output, legacy_output);
                assert_eq!(legacy_output & !legacy.external_mask, 0);
                assert_eq!(in_place_output & !in_place.external_mask, 0);
                assert_eq!(
                    legacy_output & legacy.length_mask,
                    input & legacy.length_mask
                );
                assert_eq!(
                    in_place_output & in_place.length_mask,
                    input & in_place.length_mask
                );
                assert_eq!(
                    apply_scalar(&legacy_opposite.builder.ops, legacy_output),
                    input
                );
                assert_eq!(
                    apply_scalar(&in_place_opposite.builder.ops, in_place_output),
                    input
                );
                if input & in_place.control_mask == 0 {
                    assert_eq!(legacy_output, input);
                    assert_eq!(in_place_output, input);
                    control_off_checks += 1;
                }
                basis_states_checked += 1;
                scalar_equivalence_checks += 1;
                length_restore_checks += 1;
                inverse_pair_checks += 2;
            }
        }
    }

    assert_eq!(legacy_streams_checked, 14);
    assert_eq!(copied_cursor_lanes_removed, 27);
    InPlaceLtCursorLegacyDifferentialProofReport {
        less_than_configurations_checked: less_than_configurations.len(),
        coefficient_add_configurations_checked: coefficient_add_configurations.len(),
        coefficient_add_directions_checked: 2,
        basis_states_checked,
        scalar_equivalence_checks,
        simulator_equivalence_checks,
        control_off_checks,
        length_restore_checks,
        inverse_pair_checks,
        phase_clean_checks,
        ancilla_clean_checks,
        legacy_streams_checked,
        copied_cursor_lanes_removed,
        local_ops_removed: 3 * copied_cursor_lanes_removed,
        local_cx_removed: 2 * copied_cursor_lanes_removed,
        local_resets_removed: copied_cursor_lanes_removed,
    }
}
