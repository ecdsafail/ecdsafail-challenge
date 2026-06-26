use super::exhaustive_gated_compare_check;

#[test]
fn gated_compare_and_gate_hold_exhaustive_widths_1_through_5() {
    let report = exhaustive_gated_compare_check();
    assert_eq!(report.widths_checked, 5);
    assert_eq!(report.comparator_states_checked, 5_456);
    assert_eq!(report.gate_hold_states_checked, 10_912);
    assert_eq!(report.max_comparator_extra_qubits, 1);
    assert_eq!(report.max_gate_hold_extra_qubits, 1);
}
