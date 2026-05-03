//! Scratch-600 architecture frontier tests.
//!
//! Executable accounting for candidate architectures that could plausibly live
//! in the Google-low-qubit regime: tx,ty plus <=600--663 live quantum scratch.
//! This keeps selector/parser/cleanup costs visible before any full hook-up.

#![cfg(test)]

#[derive(Clone, Copy, Debug)]
struct Candidate {
    name: &'static str,
    scratch_bits: usize,
    charged_toffoli: Option<usize>,
    blocker: &'static str,
}

#[test]
fn scratch600_frontier_requires_selector_or_parser_breakthrough() {
    const STRICT_SCRATCH: usize = 600;
    const GOOGLE_LOW_QUBIT_SCRATCH: usize = 663; // 1175 total - tx,ty=512.
    const GOOGLE_LOW_QUBIT_TOFFOLI: usize = 2_700_000;

    let candidates = [
        Candidate {
            name: "streamed_mask_qoffset_plus_lowword_selector",
            scratch_bits: 510,
            charged_toffoli: Some(2_765_676),
            blocker: "lowword selector is 120480 CCX over the 87840 selector margin",
        },
        Candidate {
            name: "by_consumed_high_state_selector",
            scratch_bits: 3_892,
            charged_toffoli: Some(3_917_624),
            blocker: "consumed lowword q/high-state update projects 1217624 CCX over target before matrix selection and q-history cleanup",
        },
        Candidate {
            name: "by_tiny_consumed_high_state_selector",
            scratch_bits: 3_916,
            charged_toffoli: Some(4_370_270),
            blocker: "smaller w=1/2/4 lowword windows do not rescue the consumed-high update; best w4 projection is 1670270 CCX over target before matrix selection and q-history cleanup",
        },
        Candidate {
            name: "partial_prefix32_qoffset_lowword_model",
            scratch_bits: 542,
            charged_toffoli: None,
            blocker: "one-DIV local pieces project 2697524, but adversarial two-denominator ledger misses by 1368262",
        },
        Candidate {
            name: "partial_prefix48_qoffset_lowword_model",
            scratch_bits: 558,
            charged_toffoli: None,
            blocker: "one-DIV local pieces project 2652404, but no charged algebra deletes the second denominator/replay",
        },
        Candidate {
            name: "partial_prefix80_qoffset_lowword_model",
            scratch_bits: 590,
            charged_toffoli: None,
            blocker: "one-DIV local pieces project 2562164, but only 10 scratch bits remain and two-denominator point-add is not viable",
        },
        Candidate {
            name: "partial_prefix90_qoffset_lowword_model",
            scratch_bits: 600,
            charged_toffoli: None,
            blocker: "one-DIV local pieces project 2533964 at strict scratch cap, but two-denominator ledger projects 4068262",
        },
        Candidate {
            name: "streamed_mask_qoffset_replay_body_only",
            scratch_bits: 510,
            charged_toffoli: None,
            blocker: "replay body projects 2645196 but selector is deliberately uncharged",
        },
        Candidate {
            name: "tiny_lowword_selector_without_den_update",
            scratch_bits: 510,
            charged_toffoli: None,
            blocker: "w1 selector-only model projects 2664876, but the best tiny-window fixed-matrix update is still 304132 CCX over selector slack",
        },
        Candidate {
            name: "full_ratio_by_selector_state",
            scratch_bits: 560,
            charged_toffoli: Some(9_952_686),
            blocker: "state fits, but A-step ratio inverse proxy projects to 9952686 total",
        },
        Candidate {
            name: "compact_by_denpair_plus_sidecar",
            scratch_bits: 564,
            charged_toffoli: Some(3_793_920),
            blocker: "state fits, direct denominator compute+uncompute is too costly",
        },
        Candidate {
            name: "plusminus_raw_k_stream_without_parser",
            scratch_bits: 564,
            charged_toffoli: None,
            blocker: "raw stream fits only before boundary/rank/live-parser cost is charged",
        },
        Candidate {
            name: "plusminus_scaled_konly_slack_denominator_blocked",
            scratch_bits: 517,
            charged_toffoli: None,
            blocker: "sampled active-chain/Solinas model treats quantum k-history as an executed-gate filter; emitted 257-bit active step is 138771 CCX, so two-DIV step-only is 56063484",
        },
        Candidate {
            name: "centered_euclid_raw_q_stream_without_parser",
            scratch_bits: 592,
            charged_toffoli: None,
            blocker: "raw stream fits only before parser/rank/live-recompute cost is charged",
        },
        Candidate {
            name: "direct_centered_signnorm_raw_digits_only",
            scratch_bits: 653,
            charged_toffoli: None,
            blocker: "raw sign-normalized digits fit, but exact cneg p99 is 2792914; norm signs have dense MBU parity and exact toy reverse collisions",
        },
        Candidate {
            name: "direct_centered_signnorm_logical_coeff_signs",
            scratch_bits: 765,
            charged_toffoli: Some(2_746_960),
            blocker: "logical coefficient signs reduce exact cneg work, but split p99 is still 46960 over target and normalization-sign scratch remains 765 p99",
        },
        Candidate {
            name: "direct_centered_restoring_final_stored_alignment",
            scratch_bits: 602,
            charged_toffoli: None,
            blocker: "restoring-final select1 now has phase-clean toy cleanup; stored alignment+branch decoder averages 2645270 with raw variable metadata p99 602, but delimited/gamma metadata is 719/809 and exact packed parser/cleanup is unbuilt",
        },
        Candidate {
            name: "direct_centered_signnorm_rank_compressed_signs",
            scratch_bits: 765,
            charged_toffoli: None,
            blocker: "even combinatorial/rank-compressed normalization signs need 765 p99 scratch bits, 102 over Google",
        },
        Candidate {
            name: "halfgcd_first_matrix_checkpoint_only",
            scratch_bits: 524,
            charged_toffoli: None,
            blocker: "matrix alone fits, but matrix+residual/tail exceeds scratch",
        },
        Candidate {
            name: "halfgcd_det_compressed_matrix_tail_payload",
            scratch_bits: 564,
            charged_toffoli: None,
            blocker: "compressed payload/replay fits, but straight-line prefix generation needs 769 bits and optimistic in-loop determinant recovery projects 4491940 Toffoli",
        },
        Candidate {
            name: "halfgcd_second_column_tail_stream",
            scratch_bits: 514,
            charged_toffoli: None,
            blocker: "second-column exact decoder average model fits at 2606688 and fixed-bound active prefix toy cleans, but active-control Bennett cleanup averages 3462517 and p99 remains 3705990",
        },
        Candidate {
            name: "folded_kaliski_one_pair_plus_required_sidecar",
            scratch_bits: 512 + 255,
            charged_toffoli: Some(4_089_274),
            blocker: "branch-recovery sidecar pushes folded Kaliski over scratch",
        },
    ];

    let best_state = candidates.iter().map(|c| c.scratch_bits).min().unwrap();
    let best_charged_sota_shaped = candidates
        .iter()
        .filter(|c| c.scratch_bits <= STRICT_SCRATCH)
        .filter_map(|c| c.charged_toffoli.map(|t| (c.name, c.scratch_bits, t)))
        .min_by_key(|(_, _, t)| *t)
        .unwrap();

    let streamed_selector_budget = 87_840usize;
    let streamed_lowword_selector = 208_320usize;
    let streamed_selector_shortfall = streamed_lowword_selector - streamed_selector_budget;
    let streamed_gap_to_google = best_charged_sota_shaped.2 as isize - GOOGLE_LOW_QUBIT_TOFFOLI as isize;

    let streamed_replay_body_projection = 2_645_196usize;
    let streamed_replay_unfunded_selector_budget =
        GOOGLE_LOW_QUBIT_TOFFOLI - streamed_replay_body_projection;
    let tiny_lowword_w1_selector_projection = 2_664_876usize;
    let tiny_lowword_w1_selector_slack =
        GOOGLE_LOW_QUBIT_TOFFOLI - tiny_lowword_w1_selector_projection;
    let tiny_lowword_best_fixed_update_excess = 304_132usize;
    let partial_prefix32_projection = 2_697_524usize;
    let partial_prefix48_projection = 2_652_404usize;
    let partial_prefix80_projection = 2_562_164usize;
    let partial_prefix90_projection = 2_533_964usize;
    let partial_prefix_two_den_projection = 4_068_262usize;
    let partial_prefix32_gap = partial_prefix32_projection as isize - GOOGLE_LOW_QUBIT_TOFFOLI as isize;
    let partial_prefix48_gap = partial_prefix48_projection as isize - GOOGLE_LOW_QUBIT_TOFFOLI as isize;
    let partial_prefix80_gap = partial_prefix80_projection as isize - GOOGLE_LOW_QUBIT_TOFFOLI as isize;
    let partial_prefix90_gap = partial_prefix90_projection as isize - GOOGLE_LOW_QUBIT_TOFFOLI as isize;
    let partial_prefix_two_den_gap = partial_prefix_two_den_projection as isize - GOOGLE_LOW_QUBIT_TOFFOLI as isize;
    let by_consumed_high_update_mean_compute_ccx = 515_494usize;
    let by_consumed_high_update_compute_uncompute_ccx = 1_030_988usize;
    let by_consumed_high_q_oracle_total_ccx = 329_280usize;
    let by_consumed_high_optimistic_pointadd = 3_917_624usize;
    let by_consumed_high_gap_to_2700k =
        by_consumed_high_optimistic_pointadd as isize - GOOGLE_LOW_QUBIT_TOFFOLI as isize;
    let by_consumed_high_max_peak_q = 3_892usize;
    let by_tiny_consumed_high_best_w = 4usize;
    let by_tiny_consumed_high_q_oracle_total_ccx = 168_000usize;
    let by_tiny_consumed_high_update_compute_ccx = 822_457usize;
    let by_tiny_consumed_high_update_compute_uncompute_ccx = 1_644_914usize;
    let by_tiny_consumed_high_optimistic_pointadd = 4_370_270usize;
    let by_tiny_consumed_high_gap_to_2700k =
        by_tiny_consumed_high_optimistic_pointadd as isize - GOOGLE_LOW_QUBIT_TOFFOLI as isize;
    let by_tiny_consumed_high_max_peak_q = 3_916usize;
    let centered_raw_scratch = 592usize;
    let centered_boundary_scratch_p99 = 710usize;
    let centered_parser_over_strict = centered_boundary_scratch_p99 - STRICT_SCRATCH;
    let direct_signnorm_raw_digit_scratch_p99 = 653usize;
    let direct_signnorm_rank_scratch_p99 = 765usize;
    let direct_signnorm_ambiguous_rank_scratch_p99 = 764usize;
    let direct_signnorm_rank_over_google =
        direct_signnorm_rank_scratch_p99 - GOOGLE_LOW_QUBIT_SCRATCH;
    let direct_signnorm_ambiguous_rank_over_google =
        direct_signnorm_ambiguous_rank_scratch_p99 - GOOGLE_LOW_QUBIT_SCRATCH;
    let direct_signnorm_exact_split_p99 = 2_792_914usize;
    let direct_signnorm_exact_split_gap =
        direct_signnorm_exact_split_p99 as isize - GOOGLE_LOW_QUBIT_TOFFOLI as isize;
    let direct_signnorm_logsign_once_p99 = 2_723_992usize;
    let direct_signnorm_logsign_split_p99 = 2_746_960usize;
    let direct_signnorm_logsign_once_gap =
        direct_signnorm_logsign_once_p99 as isize - GOOGLE_LOW_QUBIT_TOFFOLI as isize;
    let direct_signnorm_logsign_split_gap =
        direct_signnorm_logsign_split_p99 as isize - GOOGLE_LOW_QUBIT_TOFFOLI as isize;
    let direct_signnorm_mbu_degree_n14 = 13usize;
    let direct_signnorm_mbu_density_n14 = 8_208usize;
    let direct_signnorm_mbu_max_count_n14 = 8usize;
    let direct_signnorm_reverse_collisions_n14 = 2_658usize;
    let direct_signnorm_reverse_states_n14 = 64_178usize;
    let direct_signnorm_reverse_total_steps_n14 = 89_008usize;
    let direct_restoring_final_coeff_width_p99 = 47_654usize;
    let direct_restoring_final_digit_payload_p99 = 362usize;
    let direct_restoring_final_raw_digit_scratch_p99 = 256usize + direct_restoring_final_digit_payload_p99;
    let direct_restoring_final_raw_digit_over_strict =
        direct_restoring_final_raw_digit_scratch_p99 as isize - STRICT_SCRATCH as isize;
    let direct_restoring_final_raw_digit_gap_google =
        direct_restoring_final_raw_digit_scratch_p99 as isize - GOOGLE_LOW_QUBIT_SCRATCH as isize;
    let direct_restoring_final_no_unit_digits_p99 = 69usize;
    let direct_restoring_final_count_p99 = 118usize;
    let direct_restoring_final_select1x_p99 = 2_423_946usize;
    let direct_restoring_final_select2x_p99 = 2_531_042usize;
    let direct_restoring_final_select3x_p99 = 2_638_012usize;
    let direct_restoring_final_select1x_gap =
        direct_restoring_final_select1x_p99 as isize - GOOGLE_LOW_QUBIT_TOFFOLI as isize;
    let direct_restoring_final_select2x_gap =
        direct_restoring_final_select2x_p99 as isize - GOOGLE_LOW_QUBIT_TOFFOLI as isize;
    let direct_restoring_final_select3x_gap =
        direct_restoring_final_select3x_p99 as isize - GOOGLE_LOW_QUBIT_TOFFOLI as isize;
    let direct_restoring_final_toy_ccx = 306usize;
    let direct_restoring_final_toy_peak_q = 85usize;
    let direct_restoring_final_toy_neg2_cases = 10_120usize;
    let direct_restoring_final_toy_zero_final_cases = 10_010usize;
    let direct_restoring_final_bennett_fast_inverse_toy_ccx = 632usize;
    let direct_restoring_final_bennett_fast_inverse_toy_peak_q = 104usize;
    let direct_restoring_final_single_selector_toy_ccx = 298usize;
    let direct_restoring_final_single_selector_toy_peak_q = 88usize;
    let direct_restoring_final_single_selector_bennett_toy_ccx = 596usize;
    let direct_restoring_final_single_selector_bennett_toy_peak_q = 108usize;
    let direct_restoring_final_payload_mbu_degree_n14 = 13usize;
    let direct_restoring_final_payload_mbu_density_n14 = 8_284usize;
    let direct_restoring_final_payload_max_n14 = 26usize;
    let direct_restoring_final_reverse_q_collisions_n14 = 0usize;
    let direct_restoring_final_reverse_q_states_n14 = 89_008usize;
    let direct_restoring_final_reverse_q_total_steps_n14 = 89_008usize;
    let direct_restoring_final_reverse_q_max_mult_n14 = 1usize;
    let direct_restoring_final_residual_q_collisions_n14 = 4_248usize;
    let direct_restoring_final_residual_q_states_n14 = 61_008usize;
    let direct_restoring_final_residual_q_total_steps_n14 = 89_008usize;
    let direct_restoring_final_residual_q_max_mult_n14 = 64usize;
    let direct_restoring_final_reverse_coeff_candidates_transitions_n14 = 105_388usize;
    let direct_restoring_final_reverse_coeff_candidates_endpoints_n14 = 16_380usize;
    let direct_restoring_final_reverse_coeff_candidates_low_n14 = 48_896usize;
    let direct_restoring_final_reverse_coeff_candidates_high_n14 = 40_112usize;
    let direct_restoring_final_reverse_coeff_candidates_exact_n14 = 16_380usize;
    let direct_restoring_final_reverse_coeff_high_branch_degree_n14 = 13usize;
    let direct_restoring_final_reverse_coeff_high_branch_density_n14 = 8_208usize;
    let direct_restoring_final_reverse_coeff_high_branch_max_count_n14 = 8usize;
    let direct_restoring_final_reverse_coeff_high_branch_total_n14 = 40_112usize;
    let direct_restoring_final_reverse_coeff_candidates_max_q_bits_n14 = 14usize;
    let direct_restoring_final_reverse_coeff_candidates_max_coeff_abs_bits_n14 = 14usize;
    let direct_restoring_final_coeff_decoder_exact_p99 = 185_694usize;
    let direct_restoring_final_coeff_decoder_digit_width_p99 = 46_950usize;
    let direct_restoring_final_coeff_decoder_scan_p99 = 138_744usize;
    let direct_restoring_final_coeff_decoder_steps_p99 = 119usize;
    let direct_restoring_final_coeff_decoder_digits_p99 = 358usize;
    let direct_restoring_final_coeff_decoder_oneway_margin = 15_497usize;
    let direct_restoring_final_coeff_decoder_margin = -170_197isize;
    let direct_restoring_final_coeff_decoder_augmented_pointadd_p99 = 3_380_788usize;
    let direct_restoring_final_coeff_decoder_augmented_gap = 680_788isize;
    let direct_restoring_final_avg_select3_mean = 2_480_906usize;
    let direct_restoring_final_avg_select3_first64 = 2_486_059usize;
    let direct_restoring_final_avg_select3_p99 = 2_638_276usize;
    let direct_restoring_final_avg_decoder_exact_mean = 166_144usize;
    let direct_restoring_final_avg_decoder_exact_p99 = 185_875usize;
    let direct_restoring_final_avg_decoder_noscan_mean = 44_342usize;
    let direct_restoring_final_avg_decoder_noscan_p99 = 47_094usize;
    let direct_restoring_final_avg_exact_select3_mean = 3_145_482usize;
    let direct_restoring_final_avg_exact_select3_first64 = 3_152_638usize;
    let direct_restoring_final_avg_exact_select3_p99 = 3_375_572usize;
    let direct_restoring_final_avg_exact_select3_gap =
        direct_restoring_final_avg_exact_select3_mean as isize - GOOGLE_LOW_QUBIT_TOFFOLI as isize;
    let direct_restoring_final_avg_noscan_select1_mean = 2_465_688usize;
    let direct_restoring_final_avg_noscan_select1_first64 = 2_470_688usize;
    let direct_restoring_final_avg_noscan_select1_p99 = 2_610_296usize;
    let direct_restoring_final_avg_noscan_select1_gap =
        direct_restoring_final_avg_noscan_select1_mean as isize - GOOGLE_LOW_QUBIT_TOFFOLI as isize;
    let direct_restoring_final_avg_noscan_select2_mean = 2_561_982usize;
    let direct_restoring_final_avg_noscan_select2_first64 = 2_567_281usize;
    let direct_restoring_final_avg_noscan_select2_p99 = 2_716_722usize;
    let direct_restoring_final_avg_noscan_select2_gap =
        direct_restoring_final_avg_noscan_select2_mean as isize - GOOGLE_LOW_QUBIT_TOFFOLI as isize;
    let direct_restoring_final_avg_noscan_select3_mean = 2_658_276usize;
    let direct_restoring_final_avg_noscan_select3_first64 = 2_663_875usize;
    let direct_restoring_final_avg_noscan_select3_gap =
        direct_restoring_final_avg_noscan_select3_mean as isize - GOOGLE_LOW_QUBIT_TOFFOLI as isize;
    let direct_restoring_final_avg_noscan_select3_p99 = 2_823_264usize;
    let direct_restoring_final_stored_align_select1_mean = 2_537_430usize;
    let direct_restoring_final_stored_align_select1_first64 = 2_532_052usize;
    let direct_restoring_final_stored_align_select1_p99 = 2_689_752usize;
    let direct_restoring_final_stored_align_branch_select1_mean = 2_645_270usize;
    let direct_restoring_final_stored_align_branch_select1_first64 = 2_639_465usize;
    let direct_restoring_final_stored_align_branch_select1_p99 = 2_812_592usize;
    let direct_restoring_final_stored_align_branch_select1_gap =
        direct_restoring_final_stored_align_branch_select1_mean as isize
            - GOOGLE_LOW_QUBIT_TOFFOLI as isize;
    let direct_restoring_final_stored_align_fixed_scratch_p99 = 1_318usize;
    let direct_restoring_final_stored_align_variable_scratch_p99 = 602usize;
    let direct_restoring_final_stored_align_variable_scratch_max = 615usize;
    let direct_restoring_final_stored_align_delimited_scratch_p99 = 719usize;
    let direct_restoring_final_stored_align_gamma_scratch_p99 = 809usize;
    let direct_restoring_final_stored_align_pop_barrel_p99 = 19_928usize;
    let direct_restoring_final_stored_align_branch_select_p99 = 31_033usize;
    let direct_restoring_final_stored_align_branch_count_p99 = 117usize;
    let direct_restoring_final_coeff_decoder_alignment_degree_n14 = 13usize;
    let direct_restoring_final_coeff_decoder_alignment_density_n14 = 8_278usize;
    let direct_restoring_final_coeff_decoder_alignment_max_n14 = 13usize;
    let plusminus_raw_scratch = 564usize;
    let plusminus_unary_scratch_p99 = 640usize;
    let plusminus_parser_over_strict = plusminus_unary_scratch_p99 - STRICT_SCRATCH;
    let plusminus_scaled_slack_scratch_max = 517usize;
    let plusminus_scaled_solinas_projected_max = 2_027_038usize;
    let plusminus_scaled_solinas_gap_max = plusminus_scaled_solinas_projected_max as isize - GOOGLE_LOW_QUBIT_TOFFOLI as isize;
    let plusminus_active_quantum_forward_ccx = 138_771usize;
    let plusminus_active_quantum_two_div_step_only = 56_063_484usize;
    let plusminus_active_quantum_gap_to_2700k =
        plusminus_active_quantum_two_div_step_only as isize - GOOGLE_LOW_QUBIT_TOFFOLI as isize;
    let halfgcd_matrix_only = 524usize;
    let halfgcd_matrix_tail_raw = 689usize;
    let halfgcd_tail_over_google = halfgcd_matrix_tail_raw - GOOGLE_LOW_QUBIT_SCRATCH;
    let halfgcd_det_compressed_tail = 564usize;
    let halfgcd_det_compressed_tail_gap =
        halfgcd_det_compressed_tail as isize - GOOGLE_LOW_QUBIT_SCRATCH as isize;
    let halfgcd_det_recovery_num_bits_p99 = 262usize;
    let halfgcd_det_recovery_den_bits_p99 = 128usize;
    let halfgcd_tail_raw_rank_max_mult_n14 = 1usize;
    let halfgcd_tail_raw_rank_degree_n14 = 0usize;
    let halfgcd_tail_raw_rank_density_n14 = 0usize;
    let halfgcd_tail_raw_compressed_rank_max_mult_n14 = 1usize;
    let halfgcd_tail_raw_compressed_rank_degree_n14 = 0usize;
    let halfgcd_tail_raw_compressed_rank_density_n14 = 0usize;
    let halfgcd_matrix_apply_p99_ccx = 236_313usize;
    let halfgcd_tail_replay_p99_ccx = 102_725usize;
    let halfgcd_det_recovery_floor_p99_ccx = 52_757usize;
    let halfgcd_replay_with_recovery_floor_pointadd_p99 = 1_410_512usize;
    let halfgcd_replay_with_recovery_floor_gap_to_2700k =
        halfgcd_replay_with_recovery_floor_pointadd_p99 as isize - GOOGLE_LOW_QUBIT_TOFFOLI as isize;
    let halfgcd_full_prefix_live_p99_bits = 769usize;
    let halfgcd_full_prefix_live_gap_google =
        halfgcd_full_prefix_live_p99_bits as isize - GOOGLE_LOW_QUBIT_SCRATCH as isize;
    let halfgcd_compressed_residual_live_p99_bits = 646usize;
    let halfgcd_compressed_tail_stream_peak_p99_bits = 646usize;
    let halfgcd_compressed_tail_stream_peak_gap_google =
        halfgcd_compressed_tail_stream_peak_p99_bits as isize - GOOGLE_LOW_QUBIT_SCRATCH as isize;
    let halfgcd_inloop_prefix_steps_p99 = 92usize;
    let halfgcd_inloop_recovery_floor_p99_ccx = 1_540_714usize;
    let halfgcd_inloop_recovery_pointadd_p99 = 4_491_940usize;
    let halfgcd_inloop_recovery_gap_to_2700k =
        halfgcd_inloop_recovery_pointadd_p99 as isize - GOOGLE_LOW_QUBIT_TOFFOLI as isize;
    let halfgcd_second_col_bits_p99 = 265usize;
    let halfgcd_second_col_residual_bits_p99 = 256usize;
    let halfgcd_second_col_residual_live_p99_bits = 514usize;
    let halfgcd_second_col_tail_raw_bits_p99 = 434usize;
    let halfgcd_second_col_tail_stream_peak_p99_bits = 514usize;
    let halfgcd_second_col_tail_stream_peak_gap_google =
        halfgcd_second_col_tail_stream_peak_p99_bits as isize - GOOGLE_LOW_QUBIT_SCRATCH as isize;
    let halfgcd_second_col_tail_raw_rank_max_mult_n14 = 1usize;
    let halfgcd_second_col_tail_raw_rank_degree_n14 = 0usize;
    let halfgcd_second_col_tail_raw_rank_density_n14 = 0usize;
    let halfgcd_second_col_prefix_final_bd_max_mult_n14 = 1usize;
    let halfgcd_second_col_prefix_local_reverse_max_mult_n14 = 1usize;
    let halfgcd_second_col_prefix_local_reverse_collisions_n14 = 0usize;
    let halfgcd_second_col_prefix_transitions_n14 = 82_028usize;
    let halfgcd_second_col_prefix_residual_q_collisions_n14 = 2_184usize;
    let halfgcd_second_col_prefix_residual_q_states_n14 = 77_893usize;
    let halfgcd_second_col_prefix_residual_q_total_steps_n14 = 82_028usize;
    let halfgcd_second_col_prefix_residual_q_max_mult_n14 = 53usize;
    let halfgcd_second_col_prefix_reverse_formula_transitions_n14 = 82_028usize;
    let halfgcd_second_col_prefix_reverse_formula_endpoints_n14 = 16_380usize;
    let halfgcd_second_col_prefix_reverse_formula_coeff_steps_n14 = 65_648usize;
    let halfgcd_second_col_prefix_reverse_formula_max_q_bits_n14 = 14usize;
    let halfgcd_second_col_prefix_reverse_formula_max_coeff_abs_bits_n14 = 14usize;
    let halfgcd_second_col_prefix_coeff_decoder_exact_p99 = 81_879usize;
    let halfgcd_second_col_prefix_coeff_decoder_digit_width_p99 = 15_046usize;
    let halfgcd_second_col_prefix_coeff_decoder_final_fix_p99 = 12_077usize;
    let halfgcd_second_col_prefix_coeff_decoder_scan_p99 = 54_756usize;
    let halfgcd_second_col_prefix_coeff_decoder_steps_p99 = 91usize;
    let halfgcd_second_col_prefix_coeff_decoder_digits_p99 = 224usize;
    let halfgcd_second_col_prefix_coeff_decoder_final_negative_p99 = 25usize;
    let halfgcd_second_col_prefix_augmented_extraction_p99 = 397_001usize;
    let halfgcd_second_col_prefix_augmented_pointadd_p99 = 2_866_082usize;
    let halfgcd_second_col_prefix_augmented_gap_to_2700k =
        halfgcd_second_col_prefix_augmented_pointadd_p99 as isize - GOOGLE_LOW_QUBIT_TOFFOLI as isize;
    let halfgcd_second_col_prefix_steps_p99 = 91usize;
    let halfgcd_second_col_prefix_digits_p99 = 221usize;
    let halfgcd_second_col_prefix_final_negative_p99 = 37usize;
    let halfgcd_second_col_prefix_bounded_barrel_bits = 5usize;
    let halfgcd_second_col_prefix_residual_digit_width_p99 = 42_514usize;
    let halfgcd_second_col_prefix_coeff_digit_width_p99 = 14_949usize;
    let halfgcd_second_col_prefix_final_fix_width_p99 = 46_732usize;
    let halfgcd_second_col_prefix_oneway_budget_ccx = 345_059usize;
    let halfgcd_second_col_prefix_bounded_extraction_p99 = 244_769usize;
    let halfgcd_second_col_prefix_exact_extraction_p99 = 315_122usize;
    let halfgcd_second_col_prefix_exact_pointadd_p99 = 2_538_566usize;
    let halfgcd_second_col_prefix_exact_gap_to_2700k =
        halfgcd_second_col_prefix_exact_pointadd_p99 as isize - GOOGLE_LOW_QUBIT_TOFFOLI as isize;
    let halfgcd_second_col_prefix_coeff_decoder_no_scan_p99 = 26_796usize;
    let halfgcd_second_col_prefix_coeff_decoder_scan_budget = 3_141usize;
    let halfgcd_second_col_prefix_coeff_decoder_scan_over_budget = 51_975isize;
    let halfgcd_second_col_prefix_avg_exact_base_mean = 2_336_737usize;
    let halfgcd_second_col_prefix_avg_exact_base_first64 = 2_347_685usize;
    let halfgcd_second_col_prefix_avg_exact_base_p99 = 2_539_226usize;
    let halfgcd_second_col_prefix_avg_decoder_exact_mean = 67_488usize;
    let halfgcd_second_col_prefix_avg_decoder_exact_p99 = 81_425usize;
    let halfgcd_second_col_prefix_avg_decoder_noscan_mean = 23_150usize;
    let halfgcd_second_col_prefix_avg_decoder_noscan_p99 = 26_641usize;
    let halfgcd_second_col_prefix_avg_aug_exact_mean = 2_606_688usize;
    let halfgcd_second_col_prefix_avg_aug_exact_first64 = 2_624_897usize;
    let halfgcd_second_col_prefix_avg_aug_exact_p99 = 2_856_574usize;
    let halfgcd_second_col_prefix_avg_aug_exact_gap =
        halfgcd_second_col_prefix_avg_aug_exact_mean as isize - GOOGLE_LOW_QUBIT_TOFFOLI as isize;
    let halfgcd_second_col_prefix_avg_aug_noscan_mean = 2_429_337usize;
    let halfgcd_second_col_prefix_avg_aug_noscan_first64 = 2_442_020usize;
    let halfgcd_second_col_prefix_avg_aug_noscan_p99 = 2_643_668usize;
    let halfgcd_second_col_prefix_avg_aug_noscan_gap =
        halfgcd_second_col_prefix_avg_aug_noscan_mean as isize - GOOGLE_LOW_QUBIT_TOFFOLI as isize;
    let halfgcd_second_col_prefix_step_toy_ccx = 308usize;
    let halfgcd_second_col_prefix_step_toy_peak_q = 106usize;
    let halfgcd_second_col_prefix_step_toy_final_negative_cases = 39_270usize;
    let halfgcd_second_col_prefix_fixed_bound_active_toy_ccx = 11_464usize;
    let halfgcd_second_col_prefix_fixed_bound_active_toy_peak_q = 173usize;
    let halfgcd_second_col_prefix_fixed_bound_active_toy_active_slots = 123_900usize;
    let halfgcd_second_col_prefix_fixed_bound_active_toy_inactive_slots = 77_700usize;
    let halfgcd_second_col_prefix_fixed_bound_active_toy_halted_inputs = 40_180usize;
    let halfgcd_second_col_prefix_fixed_bound_active_toy_full_bound_inputs = 10_220usize;
    let halfgcd_second_col_prefix_fixed_bound_active_toy_dirty_restore_cases = 0usize;
    let halfgcd_second_col_prefix_fixed_bound_active_toy_dirty_history_cases = 0usize;
    let halfgcd_second_col_prefix_fixed_bound_active_toy_dirty_phase_cases = 0usize;
    let halfgcd_second_col_prefix_active_model_base_mean = 2_337_422usize;
    let halfgcd_second_col_prefix_active_model_base_first64 = 2_343_957usize;
    let halfgcd_second_col_prefix_active_model_oneway_mean = 548_634usize;
    let halfgcd_second_col_prefix_active_model_oneway_p99 = 607_653usize;
    let halfgcd_second_col_prefix_active_model_pointadd_mean = 3_462_517usize;
    let halfgcd_second_col_prefix_active_model_pointadd_first64 = 3_470_918usize;
    let halfgcd_second_col_prefix_active_model_pointadd_p99 = 3_705_990usize;
    let halfgcd_second_col_prefix_active_model_gap_to_2700k = 762_517isize;
    let halfgcd_second_col_prefix_active_model_over_exact_mean = 1_125_095usize;
    let halfgcd_second_col_prefix_active_model_over_exact_p99 = 1_177_392usize;

    eprintln!("\nScratch-600 architecture frontier:");
    for c in candidates {
        eprintln!(
            "  {:45} scratch={:4} charged_toffoli={:?} blocker={}",
            c.name, c.scratch_bits, c.charged_toffoli, c.blocker
        );
    }
    eprintln!(
        "best charged <=600-scratch row: {} scratch={} toffoli={} gap_to_2.7M={streamed_gap_to_google}",
        best_charged_sota_shaped.0, best_charged_sota_shaped.1, best_charged_sota_shaped.2,
    );

    println!("METRIC scratch600_frontier_best_scratch_bits={best_state}");
    println!("METRIC scratch600_frontier_best_charged_scratch_bits={}", best_charged_sota_shaped.1);
    println!("METRIC scratch600_frontier_best_charged_toffoli={}", best_charged_sota_shaped.2);
    println!("METRIC scratch600_frontier_best_charged_gap_to_2700k={streamed_gap_to_google}");
    println!("METRIC scratch600_streamed_replay_body_projected_toffoli={streamed_replay_body_projection}");
    println!("METRIC scratch600_streamed_unfunded_selector_budget_ccx={streamed_replay_unfunded_selector_budget}");
    println!("METRIC scratch600_streamed_selector_budget_ccx={streamed_selector_budget}");
    println!("METRIC scratch600_streamed_lowword_selector_ccx={streamed_lowword_selector}");
    println!("METRIC scratch600_streamed_selector_shortfall_ccx={streamed_selector_shortfall}");
    println!("METRIC scratch600_tiny_lowword_w1_selector_projection={tiny_lowword_w1_selector_projection}");
    println!("METRIC scratch600_tiny_lowword_w1_selector_slack={tiny_lowword_w1_selector_slack}");
    println!("METRIC scratch600_tiny_lowword_best_fixed_update_excess={tiny_lowword_best_fixed_update_excess}");
    println!("METRIC scratch600_partial_prefix32_projected_toffoli={partial_prefix32_projection}");
    println!("METRIC scratch600_partial_prefix32_gap_to_2700k={partial_prefix32_gap}");
    println!("METRIC scratch600_partial_prefix48_projected_toffoli={partial_prefix48_projection}");
    println!("METRIC scratch600_partial_prefix48_gap_to_2700k={partial_prefix48_gap}");
    println!("METRIC scratch600_partial_prefix80_projected_toffoli={partial_prefix80_projection}");
    println!("METRIC scratch600_partial_prefix80_gap_to_2700k={partial_prefix80_gap}");
    println!("METRIC scratch600_partial_prefix90_projected_toffoli={partial_prefix90_projection}");
    println!("METRIC scratch600_partial_prefix90_gap_to_2700k={partial_prefix90_gap}");
    println!("METRIC scratch600_partial_prefix_two_den_projected_toffoli={partial_prefix_two_den_projection}");
    println!("METRIC scratch600_partial_prefix_two_den_gap_to_2700k={partial_prefix_two_den_gap}");
    println!("METRIC scratch600_by_consumed_high_update_mean_compute_ccx={by_consumed_high_update_mean_compute_ccx}");
    println!("METRIC scratch600_by_consumed_high_update_compute_uncompute_ccx={by_consumed_high_update_compute_uncompute_ccx}");
    println!("METRIC scratch600_by_consumed_high_q_oracle_total_ccx={by_consumed_high_q_oracle_total_ccx}");
    println!("METRIC scratch600_by_consumed_high_optimistic_pointadd={by_consumed_high_optimistic_pointadd}");
    println!("METRIC scratch600_by_consumed_high_gap_to_2700k={by_consumed_high_gap_to_2700k}");
    println!("METRIC scratch600_by_consumed_high_max_peak_q={by_consumed_high_max_peak_q}");
    println!("METRIC scratch600_by_tiny_consumed_high_best_w={by_tiny_consumed_high_best_w}");
    println!("METRIC scratch600_by_tiny_consumed_high_q_oracle_total_ccx={by_tiny_consumed_high_q_oracle_total_ccx}");
    println!("METRIC scratch600_by_tiny_consumed_high_update_compute_ccx={by_tiny_consumed_high_update_compute_ccx}");
    println!("METRIC scratch600_by_tiny_consumed_high_update_compute_uncompute_ccx={by_tiny_consumed_high_update_compute_uncompute_ccx}");
    println!("METRIC scratch600_by_tiny_consumed_high_optimistic_pointadd={by_tiny_consumed_high_optimistic_pointadd}");
    println!("METRIC scratch600_by_tiny_consumed_high_gap_to_2700k={by_tiny_consumed_high_gap_to_2700k}");
    println!("METRIC scratch600_by_tiny_consumed_high_max_peak_q={by_tiny_consumed_high_max_peak_q}");
    println!("METRIC scratch600_centered_raw_scratch_bits={centered_raw_scratch}");
    println!("METRIC scratch600_centered_boundary_scratch_p99={centered_boundary_scratch_p99}");
    println!("METRIC scratch600_centered_parser_over_strict_bits={centered_parser_over_strict}");
    println!("METRIC scratch600_direct_signnorm_raw_digit_scratch_p99={direct_signnorm_raw_digit_scratch_p99}");
    println!("METRIC scratch600_direct_signnorm_rank_scratch_p99={direct_signnorm_rank_scratch_p99}");
    println!("METRIC scratch600_direct_signnorm_rank_over_google_bits={direct_signnorm_rank_over_google}");
    println!("METRIC scratch600_direct_signnorm_ambiguous_rank_scratch_p99={direct_signnorm_ambiguous_rank_scratch_p99}");
    println!("METRIC scratch600_direct_signnorm_ambiguous_rank_over_google_bits={direct_signnorm_ambiguous_rank_over_google}");
    println!("METRIC scratch600_direct_signnorm_exact_split_p99={direct_signnorm_exact_split_p99}");
    println!("METRIC scratch600_direct_signnorm_exact_split_gap_to_2700k={direct_signnorm_exact_split_gap}");
    println!("METRIC scratch600_direct_signnorm_logsign_once_p99={direct_signnorm_logsign_once_p99}");
    println!("METRIC scratch600_direct_signnorm_logsign_split_p99={direct_signnorm_logsign_split_p99}");
    println!("METRIC scratch600_direct_signnorm_logsign_once_gap_to_2700k={direct_signnorm_logsign_once_gap}");
    println!("METRIC scratch600_direct_signnorm_logsign_split_gap_to_2700k={direct_signnorm_logsign_split_gap}");
    println!("METRIC scratch600_direct_signnorm_mbu_degree_n14={direct_signnorm_mbu_degree_n14}");
    println!("METRIC scratch600_direct_signnorm_mbu_density_n14={direct_signnorm_mbu_density_n14}");
    println!("METRIC scratch600_direct_signnorm_mbu_max_count_n14={direct_signnorm_mbu_max_count_n14}");
    println!("METRIC scratch600_direct_signnorm_reverse_collisions_n14={direct_signnorm_reverse_collisions_n14}");
    println!("METRIC scratch600_direct_signnorm_reverse_states_n14={direct_signnorm_reverse_states_n14}");
    println!("METRIC scratch600_direct_signnorm_reverse_total_steps_n14={direct_signnorm_reverse_total_steps_n14}");
    println!("METRIC scratch600_direct_restoring_final_coeff_width_p99={direct_restoring_final_coeff_width_p99}");
    println!("METRIC scratch600_direct_restoring_final_digit_payload_p99={direct_restoring_final_digit_payload_p99}");
    println!("METRIC scratch600_direct_restoring_final_raw_digit_scratch_p99={direct_restoring_final_raw_digit_scratch_p99}");
    println!("METRIC scratch600_direct_restoring_final_raw_digit_over_strict_bits={direct_restoring_final_raw_digit_over_strict}");
    println!("METRIC scratch600_direct_restoring_final_raw_digit_gap_google_bits={direct_restoring_final_raw_digit_gap_google}");
    println!("METRIC scratch600_direct_restoring_final_no_unit_digits_p99={direct_restoring_final_no_unit_digits_p99}");
    println!("METRIC scratch600_direct_restoring_final_count_p99={direct_restoring_final_count_p99}");
    println!("METRIC scratch600_direct_restoring_final_select1x_p99={direct_restoring_final_select1x_p99}");
    println!("METRIC scratch600_direct_restoring_final_select2x_p99={direct_restoring_final_select2x_p99}");
    println!("METRIC scratch600_direct_restoring_final_select3x_p99={direct_restoring_final_select3x_p99}");
    println!("METRIC scratch600_direct_restoring_final_select1x_gap_to_2700k={direct_restoring_final_select1x_gap}");
    println!("METRIC scratch600_direct_restoring_final_select2x_gap_to_2700k={direct_restoring_final_select2x_gap}");
    println!("METRIC scratch600_direct_restoring_final_select3x_gap_to_2700k={direct_restoring_final_select3x_gap}");
    println!("METRIC scratch600_direct_restoring_final_toy_ccx={direct_restoring_final_toy_ccx}");
    println!("METRIC scratch600_direct_restoring_final_toy_peak_q={direct_restoring_final_toy_peak_q}");
    println!("METRIC scratch600_direct_restoring_final_toy_neg2_cases={direct_restoring_final_toy_neg2_cases}");
    println!("METRIC scratch600_direct_restoring_final_toy_zero_final_cases={direct_restoring_final_toy_zero_final_cases}");
    println!("METRIC scratch600_direct_restoring_final_bennett_fast_inverse_toy_ccx={direct_restoring_final_bennett_fast_inverse_toy_ccx}");
    println!("METRIC scratch600_direct_restoring_final_bennett_fast_inverse_toy_peak_q={direct_restoring_final_bennett_fast_inverse_toy_peak_q}");
    println!("METRIC scratch600_direct_restoring_final_single_selector_toy_ccx={direct_restoring_final_single_selector_toy_ccx}");
    println!("METRIC scratch600_direct_restoring_final_single_selector_toy_peak_q={direct_restoring_final_single_selector_toy_peak_q}");
    println!("METRIC scratch600_direct_restoring_final_single_selector_bennett_toy_ccx={direct_restoring_final_single_selector_bennett_toy_ccx}");
    println!("METRIC scratch600_direct_restoring_final_single_selector_bennett_toy_peak_q={direct_restoring_final_single_selector_bennett_toy_peak_q}");
    println!("METRIC scratch600_direct_restoring_final_payload_mbu_degree_n14={direct_restoring_final_payload_mbu_degree_n14}");
    println!("METRIC scratch600_direct_restoring_final_payload_mbu_density_n14={direct_restoring_final_payload_mbu_density_n14}");
    println!("METRIC scratch600_direct_restoring_final_payload_max_n14={direct_restoring_final_payload_max_n14}");
    println!("METRIC scratch600_direct_restoring_final_reverse_q_collisions_n14={direct_restoring_final_reverse_q_collisions_n14}");
    println!("METRIC scratch600_direct_restoring_final_reverse_q_states_n14={direct_restoring_final_reverse_q_states_n14}");
    println!("METRIC scratch600_direct_restoring_final_reverse_q_total_steps_n14={direct_restoring_final_reverse_q_total_steps_n14}");
    println!("METRIC scratch600_direct_restoring_final_reverse_q_max_mult_n14={direct_restoring_final_reverse_q_max_mult_n14}");
    println!("METRIC scratch600_direct_restoring_final_residual_q_collisions_n14={direct_restoring_final_residual_q_collisions_n14}");
    println!("METRIC scratch600_direct_restoring_final_residual_q_states_n14={direct_restoring_final_residual_q_states_n14}");
    println!("METRIC scratch600_direct_restoring_final_residual_q_total_steps_n14={direct_restoring_final_residual_q_total_steps_n14}");
    println!("METRIC scratch600_direct_restoring_final_residual_q_max_mult_n14={direct_restoring_final_residual_q_max_mult_n14}");
    println!("METRIC scratch600_direct_restoring_final_reverse_coeff_candidates_transitions_n14={direct_restoring_final_reverse_coeff_candidates_transitions_n14}");
    println!("METRIC scratch600_direct_restoring_final_reverse_coeff_candidates_endpoints_n14={direct_restoring_final_reverse_coeff_candidates_endpoints_n14}");
    println!("METRIC scratch600_direct_restoring_final_reverse_coeff_candidates_low_n14={direct_restoring_final_reverse_coeff_candidates_low_n14}");
    println!("METRIC scratch600_direct_restoring_final_reverse_coeff_candidates_high_n14={direct_restoring_final_reverse_coeff_candidates_high_n14}");
    println!("METRIC scratch600_direct_restoring_final_reverse_coeff_candidates_exact_n14={direct_restoring_final_reverse_coeff_candidates_exact_n14}");
    println!("METRIC scratch600_direct_restoring_final_reverse_coeff_high_branch_degree_n14={direct_restoring_final_reverse_coeff_high_branch_degree_n14}");
    println!("METRIC scratch600_direct_restoring_final_reverse_coeff_high_branch_density_n14={direct_restoring_final_reverse_coeff_high_branch_density_n14}");
    println!("METRIC scratch600_direct_restoring_final_reverse_coeff_high_branch_max_count_n14={direct_restoring_final_reverse_coeff_high_branch_max_count_n14}");
    println!("METRIC scratch600_direct_restoring_final_reverse_coeff_high_branch_total_n14={direct_restoring_final_reverse_coeff_high_branch_total_n14}");
    println!("METRIC scratch600_direct_restoring_final_reverse_coeff_candidates_max_q_bits_n14={direct_restoring_final_reverse_coeff_candidates_max_q_bits_n14}");
    println!("METRIC scratch600_direct_restoring_final_reverse_coeff_candidates_max_coeff_abs_bits_n14={direct_restoring_final_reverse_coeff_candidates_max_coeff_abs_bits_n14}");
    println!("METRIC scratch600_direct_restoring_final_coeff_decoder_exact_p99={direct_restoring_final_coeff_decoder_exact_p99}");
    println!("METRIC scratch600_direct_restoring_final_coeff_decoder_digit_width_p99={direct_restoring_final_coeff_decoder_digit_width_p99}");
    println!("METRIC scratch600_direct_restoring_final_coeff_decoder_scan_p99={direct_restoring_final_coeff_decoder_scan_p99}");
    println!("METRIC scratch600_direct_restoring_final_coeff_decoder_steps_p99={direct_restoring_final_coeff_decoder_steps_p99}");
    println!("METRIC scratch600_direct_restoring_final_coeff_decoder_digits_p99={direct_restoring_final_coeff_decoder_digits_p99}");
    println!("METRIC scratch600_direct_restoring_final_coeff_decoder_oneway_margin={direct_restoring_final_coeff_decoder_oneway_margin}");
    println!("METRIC scratch600_direct_restoring_final_coeff_decoder_margin={direct_restoring_final_coeff_decoder_margin}");
    println!("METRIC scratch600_direct_restoring_final_coeff_decoder_augmented_pointadd_p99={direct_restoring_final_coeff_decoder_augmented_pointadd_p99}");
    println!("METRIC scratch600_direct_restoring_final_coeff_decoder_augmented_gap_to_2700k={direct_restoring_final_coeff_decoder_augmented_gap}");
    println!("METRIC scratch600_direct_restoring_final_avg_select3_mean={direct_restoring_final_avg_select3_mean}");
    println!("METRIC scratch600_direct_restoring_final_avg_select3_first64={direct_restoring_final_avg_select3_first64}");
    println!("METRIC scratch600_direct_restoring_final_avg_select3_p99={direct_restoring_final_avg_select3_p99}");
    println!("METRIC scratch600_direct_restoring_final_avg_decoder_exact_mean={direct_restoring_final_avg_decoder_exact_mean}");
    println!("METRIC scratch600_direct_restoring_final_avg_decoder_exact_p99={direct_restoring_final_avg_decoder_exact_p99}");
    println!("METRIC scratch600_direct_restoring_final_avg_decoder_noscan_mean={direct_restoring_final_avg_decoder_noscan_mean}");
    println!("METRIC scratch600_direct_restoring_final_avg_decoder_noscan_p99={direct_restoring_final_avg_decoder_noscan_p99}");
    println!("METRIC scratch600_direct_restoring_final_avg_exact_select3_mean={direct_restoring_final_avg_exact_select3_mean}");
    println!("METRIC scratch600_direct_restoring_final_avg_exact_select3_first64={direct_restoring_final_avg_exact_select3_first64}");
    println!("METRIC scratch600_direct_restoring_final_avg_exact_select3_p99={direct_restoring_final_avg_exact_select3_p99}");
    println!("METRIC scratch600_direct_restoring_final_avg_exact_select3_gap_to_2700k={direct_restoring_final_avg_exact_select3_gap}");
    println!("METRIC scratch600_direct_restoring_final_avg_noscan_select1_mean={direct_restoring_final_avg_noscan_select1_mean}");
    println!("METRIC scratch600_direct_restoring_final_avg_noscan_select1_first64={direct_restoring_final_avg_noscan_select1_first64}");
    println!("METRIC scratch600_direct_restoring_final_avg_noscan_select1_p99={direct_restoring_final_avg_noscan_select1_p99}");
    println!("METRIC scratch600_direct_restoring_final_avg_noscan_select1_gap_to_2700k={direct_restoring_final_avg_noscan_select1_gap}");
    println!("METRIC scratch600_direct_restoring_final_avg_noscan_select2_mean={direct_restoring_final_avg_noscan_select2_mean}");
    println!("METRIC scratch600_direct_restoring_final_avg_noscan_select2_first64={direct_restoring_final_avg_noscan_select2_first64}");
    println!("METRIC scratch600_direct_restoring_final_avg_noscan_select2_p99={direct_restoring_final_avg_noscan_select2_p99}");
    println!("METRIC scratch600_direct_restoring_final_avg_noscan_select2_gap_to_2700k={direct_restoring_final_avg_noscan_select2_gap}");
    println!("METRIC scratch600_direct_restoring_final_avg_noscan_select3_mean={direct_restoring_final_avg_noscan_select3_mean}");
    println!("METRIC scratch600_direct_restoring_final_avg_noscan_select3_first64={direct_restoring_final_avg_noscan_select3_first64}");
    println!("METRIC scratch600_direct_restoring_final_avg_noscan_select3_gap_to_2700k={direct_restoring_final_avg_noscan_select3_gap}");
    println!("METRIC scratch600_direct_restoring_final_avg_noscan_select3_p99={direct_restoring_final_avg_noscan_select3_p99}");
    println!("METRIC scratch600_direct_restoring_final_stored_align_select1_mean={direct_restoring_final_stored_align_select1_mean}");
    println!("METRIC scratch600_direct_restoring_final_stored_align_select1_first64={direct_restoring_final_stored_align_select1_first64}");
    println!("METRIC scratch600_direct_restoring_final_stored_align_select1_p99={direct_restoring_final_stored_align_select1_p99}");
    println!("METRIC scratch600_direct_restoring_final_stored_align_branch_select1_mean={direct_restoring_final_stored_align_branch_select1_mean}");
    println!("METRIC scratch600_direct_restoring_final_stored_align_branch_select1_first64={direct_restoring_final_stored_align_branch_select1_first64}");
    println!("METRIC scratch600_direct_restoring_final_stored_align_branch_select1_p99={direct_restoring_final_stored_align_branch_select1_p99}");
    println!("METRIC scratch600_direct_restoring_final_stored_align_branch_select1_gap_to_2700k={direct_restoring_final_stored_align_branch_select1_gap}");
    println!("METRIC scratch600_direct_restoring_final_stored_align_fixed_scratch_p99={direct_restoring_final_stored_align_fixed_scratch_p99}");
    println!("METRIC scratch600_direct_restoring_final_stored_align_variable_scratch_p99={direct_restoring_final_stored_align_variable_scratch_p99}");
    println!("METRIC scratch600_direct_restoring_final_stored_align_variable_scratch_max={direct_restoring_final_stored_align_variable_scratch_max}");
    println!("METRIC scratch600_direct_restoring_final_stored_align_delimited_scratch_p99={direct_restoring_final_stored_align_delimited_scratch_p99}");
    println!("METRIC scratch600_direct_restoring_final_stored_align_gamma_scratch_p99={direct_restoring_final_stored_align_gamma_scratch_p99}");
    println!("METRIC scratch600_direct_restoring_final_stored_align_pop_barrel_p99={direct_restoring_final_stored_align_pop_barrel_p99}");
    println!("METRIC scratch600_direct_restoring_final_stored_align_branch_select_p99={direct_restoring_final_stored_align_branch_select_p99}");
    println!("METRIC scratch600_direct_restoring_final_stored_align_branch_count_p99={direct_restoring_final_stored_align_branch_count_p99}");
    println!("METRIC scratch600_direct_restoring_final_coeff_decoder_alignment_degree_n14={direct_restoring_final_coeff_decoder_alignment_degree_n14}");
    println!("METRIC scratch600_direct_restoring_final_coeff_decoder_alignment_density_n14={direct_restoring_final_coeff_decoder_alignment_density_n14}");
    println!("METRIC scratch600_direct_restoring_final_coeff_decoder_alignment_max_n14={direct_restoring_final_coeff_decoder_alignment_max_n14}");
    println!("METRIC scratch600_plusminus_raw_scratch_bits={plusminus_raw_scratch}");
    println!("METRIC scratch600_plusminus_unary_scratch_p99={plusminus_unary_scratch_p99}");
    println!("METRIC scratch600_plusminus_parser_over_strict_bits={plusminus_parser_over_strict}");
    println!("METRIC scratch600_plusminus_scaled_slack_scratch_max={plusminus_scaled_slack_scratch_max}");
    println!("METRIC scratch600_plusminus_scaled_solinas_projected_max={plusminus_scaled_solinas_projected_max}");
    println!("METRIC scratch600_plusminus_scaled_solinas_gap_max_to_2700k={plusminus_scaled_solinas_gap_max}");
    println!("METRIC scratch600_plusminus_active_quantum_forward_ccx={plusminus_active_quantum_forward_ccx}");
    println!("METRIC scratch600_plusminus_active_quantum_two_div_step_only={plusminus_active_quantum_two_div_step_only}");
    println!("METRIC scratch600_plusminus_active_quantum_gap_to_2700k={plusminus_active_quantum_gap_to_2700k}");
    println!("METRIC scratch600_halfgcd_matrix_only_bits={halfgcd_matrix_only}");
    println!("METRIC scratch600_halfgcd_matrix_tail_raw_bits={halfgcd_matrix_tail_raw}");
    println!("METRIC scratch600_halfgcd_tail_over_google_bits={halfgcd_tail_over_google}");
    println!("METRIC scratch600_halfgcd_det_compressed_tail_bits={halfgcd_det_compressed_tail}");
    println!("METRIC scratch600_halfgcd_det_compressed_tail_gap_google={halfgcd_det_compressed_tail_gap}");
    println!("METRIC scratch600_halfgcd_det_recovery_num_bits_p99={halfgcd_det_recovery_num_bits_p99}");
    println!("METRIC scratch600_halfgcd_det_recovery_den_bits_p99={halfgcd_det_recovery_den_bits_p99}");
    println!("METRIC scratch600_halfgcd_tail_raw_rank_max_mult_n14={halfgcd_tail_raw_rank_max_mult_n14}");
    println!("METRIC scratch600_halfgcd_tail_raw_rank_degree_n14={halfgcd_tail_raw_rank_degree_n14}");
    println!("METRIC scratch600_halfgcd_tail_raw_rank_density_n14={halfgcd_tail_raw_rank_density_n14}");
    println!("METRIC scratch600_halfgcd_tail_raw_compressed_rank_max_mult_n14={halfgcd_tail_raw_compressed_rank_max_mult_n14}");
    println!("METRIC scratch600_halfgcd_tail_raw_compressed_rank_degree_n14={halfgcd_tail_raw_compressed_rank_degree_n14}");
    println!("METRIC scratch600_halfgcd_tail_raw_compressed_rank_density_n14={halfgcd_tail_raw_compressed_rank_density_n14}");
    println!("METRIC scratch600_halfgcd_matrix_apply_p99_ccx={halfgcd_matrix_apply_p99_ccx}");
    println!("METRIC scratch600_halfgcd_tail_replay_p99_ccx={halfgcd_tail_replay_p99_ccx}");
    println!("METRIC scratch600_halfgcd_det_recovery_floor_p99_ccx={halfgcd_det_recovery_floor_p99_ccx}");
    println!("METRIC scratch600_halfgcd_replay_with_recovery_floor_pointadd_p99={halfgcd_replay_with_recovery_floor_pointadd_p99}");
    println!("METRIC scratch600_halfgcd_replay_with_recovery_floor_gap_to_2700k={halfgcd_replay_with_recovery_floor_gap_to_2700k}");
    println!("METRIC scratch600_halfgcd_full_prefix_live_p99_bits={halfgcd_full_prefix_live_p99_bits}");
    println!("METRIC scratch600_halfgcd_full_prefix_live_gap_google={halfgcd_full_prefix_live_gap_google}");
    println!("METRIC scratch600_halfgcd_compressed_residual_live_p99_bits={halfgcd_compressed_residual_live_p99_bits}");
    println!("METRIC scratch600_halfgcd_compressed_tail_stream_peak_p99_bits={halfgcd_compressed_tail_stream_peak_p99_bits}");
    println!("METRIC scratch600_halfgcd_compressed_tail_stream_peak_gap_google={halfgcd_compressed_tail_stream_peak_gap_google}");
    println!("METRIC scratch600_halfgcd_inloop_prefix_steps_p99={halfgcd_inloop_prefix_steps_p99}");
    println!("METRIC scratch600_halfgcd_inloop_recovery_floor_p99_ccx={halfgcd_inloop_recovery_floor_p99_ccx}");
    println!("METRIC scratch600_halfgcd_inloop_recovery_pointadd_p99={halfgcd_inloop_recovery_pointadd_p99}");
    println!("METRIC scratch600_halfgcd_inloop_recovery_gap_to_2700k={halfgcd_inloop_recovery_gap_to_2700k}");
    println!("METRIC scratch600_halfgcd_second_col_bits_p99={halfgcd_second_col_bits_p99}");
    println!("METRIC scratch600_halfgcd_second_col_residual_bits_p99={halfgcd_second_col_residual_bits_p99}");
    println!("METRIC scratch600_halfgcd_second_col_residual_live_p99_bits={halfgcd_second_col_residual_live_p99_bits}");
    println!("METRIC scratch600_halfgcd_second_col_tail_raw_bits_p99={halfgcd_second_col_tail_raw_bits_p99}");
    println!("METRIC scratch600_halfgcd_second_col_tail_stream_peak_p99_bits={halfgcd_second_col_tail_stream_peak_p99_bits}");
    println!("METRIC scratch600_halfgcd_second_col_tail_stream_peak_gap_google={halfgcd_second_col_tail_stream_peak_gap_google}");
    println!("METRIC scratch600_halfgcd_second_col_tail_raw_rank_max_mult_n14={halfgcd_second_col_tail_raw_rank_max_mult_n14}");
    println!("METRIC scratch600_halfgcd_second_col_tail_raw_rank_degree_n14={halfgcd_second_col_tail_raw_rank_degree_n14}");
    println!("METRIC scratch600_halfgcd_second_col_tail_raw_rank_density_n14={halfgcd_second_col_tail_raw_rank_density_n14}");
    println!("METRIC scratch600_halfgcd_second_col_prefix_final_bd_max_mult_n14={halfgcd_second_col_prefix_final_bd_max_mult_n14}");
    println!("METRIC scratch600_halfgcd_second_col_prefix_local_reverse_max_mult_n14={halfgcd_second_col_prefix_local_reverse_max_mult_n14}");
    println!("METRIC scratch600_halfgcd_second_col_prefix_local_reverse_collisions_n14={halfgcd_second_col_prefix_local_reverse_collisions_n14}");
    println!("METRIC scratch600_halfgcd_second_col_prefix_transitions_n14={halfgcd_second_col_prefix_transitions_n14}");
    println!("METRIC scratch600_halfgcd_second_col_prefix_residual_q_collisions_n14={halfgcd_second_col_prefix_residual_q_collisions_n14}");
    println!("METRIC scratch600_halfgcd_second_col_prefix_residual_q_states_n14={halfgcd_second_col_prefix_residual_q_states_n14}");
    println!("METRIC scratch600_halfgcd_second_col_prefix_residual_q_total_steps_n14={halfgcd_second_col_prefix_residual_q_total_steps_n14}");
    println!("METRIC scratch600_halfgcd_second_col_prefix_residual_q_max_mult_n14={halfgcd_second_col_prefix_residual_q_max_mult_n14}");
    println!("METRIC scratch600_halfgcd_second_col_prefix_reverse_formula_transitions_n14={halfgcd_second_col_prefix_reverse_formula_transitions_n14}");
    println!("METRIC scratch600_halfgcd_second_col_prefix_reverse_formula_endpoints_n14={halfgcd_second_col_prefix_reverse_formula_endpoints_n14}");
    println!("METRIC scratch600_halfgcd_second_col_prefix_reverse_formula_coeff_steps_n14={halfgcd_second_col_prefix_reverse_formula_coeff_steps_n14}");
    println!("METRIC scratch600_halfgcd_second_col_prefix_reverse_formula_max_q_bits_n14={halfgcd_second_col_prefix_reverse_formula_max_q_bits_n14}");
    println!("METRIC scratch600_halfgcd_second_col_prefix_reverse_formula_max_coeff_abs_bits_n14={halfgcd_second_col_prefix_reverse_formula_max_coeff_abs_bits_n14}");
    println!("METRIC scratch600_halfgcd_second_col_prefix_coeff_decoder_exact_p99={halfgcd_second_col_prefix_coeff_decoder_exact_p99}");
    println!("METRIC scratch600_halfgcd_second_col_prefix_coeff_decoder_digit_width_p99={halfgcd_second_col_prefix_coeff_decoder_digit_width_p99}");
    println!("METRIC scratch600_halfgcd_second_col_prefix_coeff_decoder_final_fix_p99={halfgcd_second_col_prefix_coeff_decoder_final_fix_p99}");
    println!("METRIC scratch600_halfgcd_second_col_prefix_coeff_decoder_scan_p99={halfgcd_second_col_prefix_coeff_decoder_scan_p99}");
    println!("METRIC scratch600_halfgcd_second_col_prefix_coeff_decoder_steps_p99={halfgcd_second_col_prefix_coeff_decoder_steps_p99}");
    println!("METRIC scratch600_halfgcd_second_col_prefix_coeff_decoder_digits_p99={halfgcd_second_col_prefix_coeff_decoder_digits_p99}");
    println!("METRIC scratch600_halfgcd_second_col_prefix_coeff_decoder_final_negative_p99={halfgcd_second_col_prefix_coeff_decoder_final_negative_p99}");
    println!("METRIC scratch600_halfgcd_second_col_prefix_coeff_decoder_no_scan_p99={halfgcd_second_col_prefix_coeff_decoder_no_scan_p99}");
    println!("METRIC scratch600_halfgcd_second_col_prefix_coeff_decoder_scan_budget={halfgcd_second_col_prefix_coeff_decoder_scan_budget}");
    println!("METRIC scratch600_halfgcd_second_col_prefix_coeff_decoder_scan_over_budget={halfgcd_second_col_prefix_coeff_decoder_scan_over_budget}");
    println!("METRIC scratch600_halfgcd_second_col_prefix_augmented_extraction_p99={halfgcd_second_col_prefix_augmented_extraction_p99}");
    println!("METRIC scratch600_halfgcd_second_col_prefix_augmented_pointadd_p99={halfgcd_second_col_prefix_augmented_pointadd_p99}");
    println!("METRIC scratch600_halfgcd_second_col_prefix_augmented_gap_to_2700k={halfgcd_second_col_prefix_augmented_gap_to_2700k}");
    println!("METRIC scratch600_halfgcd_second_col_prefix_steps_p99={halfgcd_second_col_prefix_steps_p99}");
    println!("METRIC scratch600_halfgcd_second_col_prefix_digits_p99={halfgcd_second_col_prefix_digits_p99}");
    println!("METRIC scratch600_halfgcd_second_col_prefix_final_negative_p99={halfgcd_second_col_prefix_final_negative_p99}");
    println!("METRIC scratch600_halfgcd_second_col_prefix_bounded_barrel_bits={halfgcd_second_col_prefix_bounded_barrel_bits}");
    println!("METRIC scratch600_halfgcd_second_col_prefix_residual_digit_width_p99={halfgcd_second_col_prefix_residual_digit_width_p99}");
    println!("METRIC scratch600_halfgcd_second_col_prefix_coeff_digit_width_p99={halfgcd_second_col_prefix_coeff_digit_width_p99}");
    println!("METRIC scratch600_halfgcd_second_col_prefix_final_fix_width_p99={halfgcd_second_col_prefix_final_fix_width_p99}");
    println!("METRIC scratch600_halfgcd_second_col_prefix_oneway_budget_ccx={halfgcd_second_col_prefix_oneway_budget_ccx}");
    println!("METRIC scratch600_halfgcd_second_col_prefix_bounded_extraction_p99={halfgcd_second_col_prefix_bounded_extraction_p99}");
    println!("METRIC scratch600_halfgcd_second_col_prefix_exact_extraction_p99={halfgcd_second_col_prefix_exact_extraction_p99}");
    println!("METRIC scratch600_halfgcd_second_col_prefix_exact_pointadd_p99={halfgcd_second_col_prefix_exact_pointadd_p99}");
    println!("METRIC scratch600_halfgcd_second_col_prefix_exact_gap_to_2700k={halfgcd_second_col_prefix_exact_gap_to_2700k}");
    println!("METRIC scratch600_halfgcd_second_col_prefix_avg_exact_base_mean={halfgcd_second_col_prefix_avg_exact_base_mean}");
    println!("METRIC scratch600_halfgcd_second_col_prefix_avg_exact_base_first64={halfgcd_second_col_prefix_avg_exact_base_first64}");
    println!("METRIC scratch600_halfgcd_second_col_prefix_avg_exact_base_p99={halfgcd_second_col_prefix_avg_exact_base_p99}");
    println!("METRIC scratch600_halfgcd_second_col_prefix_avg_decoder_exact_mean={halfgcd_second_col_prefix_avg_decoder_exact_mean}");
    println!("METRIC scratch600_halfgcd_second_col_prefix_avg_decoder_exact_p99={halfgcd_second_col_prefix_avg_decoder_exact_p99}");
    println!("METRIC scratch600_halfgcd_second_col_prefix_avg_decoder_noscan_mean={halfgcd_second_col_prefix_avg_decoder_noscan_mean}");
    println!("METRIC scratch600_halfgcd_second_col_prefix_avg_decoder_noscan_p99={halfgcd_second_col_prefix_avg_decoder_noscan_p99}");
    println!("METRIC scratch600_halfgcd_second_col_prefix_avg_aug_exact_mean={halfgcd_second_col_prefix_avg_aug_exact_mean}");
    println!("METRIC scratch600_halfgcd_second_col_prefix_avg_aug_exact_first64={halfgcd_second_col_prefix_avg_aug_exact_first64}");
    println!("METRIC scratch600_halfgcd_second_col_prefix_avg_aug_exact_p99={halfgcd_second_col_prefix_avg_aug_exact_p99}");
    println!("METRIC scratch600_halfgcd_second_col_prefix_avg_aug_exact_gap_to_2700k={halfgcd_second_col_prefix_avg_aug_exact_gap}");
    println!("METRIC scratch600_halfgcd_second_col_prefix_avg_aug_noscan_mean={halfgcd_second_col_prefix_avg_aug_noscan_mean}");
    println!("METRIC scratch600_halfgcd_second_col_prefix_avg_aug_noscan_first64={halfgcd_second_col_prefix_avg_aug_noscan_first64}");
    println!("METRIC scratch600_halfgcd_second_col_prefix_avg_aug_noscan_p99={halfgcd_second_col_prefix_avg_aug_noscan_p99}");
    println!("METRIC scratch600_halfgcd_second_col_prefix_avg_aug_noscan_gap_to_2700k={halfgcd_second_col_prefix_avg_aug_noscan_gap}");
    println!("METRIC scratch600_halfgcd_second_col_prefix_step_toy_ccx={halfgcd_second_col_prefix_step_toy_ccx}");
    println!("METRIC scratch600_halfgcd_second_col_prefix_step_toy_peak_q={halfgcd_second_col_prefix_step_toy_peak_q}");
    println!("METRIC scratch600_halfgcd_second_col_prefix_step_toy_final_negative_cases={halfgcd_second_col_prefix_step_toy_final_negative_cases}");
    println!("METRIC scratch600_halfgcd_second_col_prefix_fixed_bound_active_toy_ccx={halfgcd_second_col_prefix_fixed_bound_active_toy_ccx}");
    println!("METRIC scratch600_halfgcd_second_col_prefix_fixed_bound_active_toy_peak_q={halfgcd_second_col_prefix_fixed_bound_active_toy_peak_q}");
    println!("METRIC scratch600_halfgcd_second_col_prefix_fixed_bound_active_toy_active_slots={halfgcd_second_col_prefix_fixed_bound_active_toy_active_slots}");
    println!("METRIC scratch600_halfgcd_second_col_prefix_fixed_bound_active_toy_inactive_slots={halfgcd_second_col_prefix_fixed_bound_active_toy_inactive_slots}");
    println!("METRIC scratch600_halfgcd_second_col_prefix_fixed_bound_active_toy_halted_inputs={halfgcd_second_col_prefix_fixed_bound_active_toy_halted_inputs}");
    println!("METRIC scratch600_halfgcd_second_col_prefix_fixed_bound_active_toy_full_bound_inputs={halfgcd_second_col_prefix_fixed_bound_active_toy_full_bound_inputs}");
    println!("METRIC scratch600_halfgcd_second_col_prefix_fixed_bound_active_toy_dirty_restore_cases={halfgcd_second_col_prefix_fixed_bound_active_toy_dirty_restore_cases}");
    println!("METRIC scratch600_halfgcd_second_col_prefix_fixed_bound_active_toy_dirty_history_cases={halfgcd_second_col_prefix_fixed_bound_active_toy_dirty_history_cases}");
    println!("METRIC scratch600_halfgcd_second_col_prefix_fixed_bound_active_toy_dirty_phase_cases={halfgcd_second_col_prefix_fixed_bound_active_toy_dirty_phase_cases}");
    println!("METRIC scratch600_halfgcd_second_col_prefix_active_model_base_mean={halfgcd_second_col_prefix_active_model_base_mean}");
    println!("METRIC scratch600_halfgcd_second_col_prefix_active_model_base_first64={halfgcd_second_col_prefix_active_model_base_first64}");
    println!("METRIC scratch600_halfgcd_second_col_prefix_active_model_oneway_mean={halfgcd_second_col_prefix_active_model_oneway_mean}");
    println!("METRIC scratch600_halfgcd_second_col_prefix_active_model_oneway_p99={halfgcd_second_col_prefix_active_model_oneway_p99}");
    println!("METRIC scratch600_halfgcd_second_col_prefix_active_model_pointadd_mean={halfgcd_second_col_prefix_active_model_pointadd_mean}");
    println!("METRIC scratch600_halfgcd_second_col_prefix_active_model_pointadd_first64={halfgcd_second_col_prefix_active_model_pointadd_first64}");
    println!("METRIC scratch600_halfgcd_second_col_prefix_active_model_pointadd_p99={halfgcd_second_col_prefix_active_model_pointadd_p99}");
    println!("METRIC scratch600_halfgcd_second_col_prefix_active_model_gap_to_2700k={halfgcd_second_col_prefix_active_model_gap_to_2700k}");
    println!("METRIC scratch600_halfgcd_second_col_prefix_active_model_over_exact_mean={halfgcd_second_col_prefix_active_model_over_exact_mean}");
    println!("METRIC scratch600_halfgcd_second_col_prefix_active_model_over_exact_p99={halfgcd_second_col_prefix_active_model_over_exact_p99}");

    assert!(best_state <= STRICT_SCRATCH, "at least some state shapes fit");
    assert!(streamed_gap_to_google > 0, "no fully charged <=600-scratch row should be counted as solved yet");
    assert!(streamed_selector_shortfall > 0, "streamed-mask route still needs a selector breakthrough");
    assert!(
        tiny_lowword_w1_selector_slack > 0 && tiny_lowword_best_fixed_update_excess > 250_000,
        "tiny lowword selector/update tradeoff changed; revisit streamed BY route"
    );
    assert!(
        by_consumed_high_gap_to_2700k > 1_000_000 && by_consumed_high_max_peak_q > GOOGLE_LOW_QUBIT_SCRATCH,
        "consumed high-state BY selector should stay demoted until a fused low-peak update exists"
    );
    assert!(
        by_tiny_consumed_high_best_w == 4
            && by_tiny_consumed_high_gap_to_2700k > 1_000_000
            && by_tiny_consumed_high_max_peak_q > GOOGLE_LOW_QUBIT_SCRATCH,
        "tiny-window consumed high-state BY selector should stay demoted until a fused low-peak update exists"
    );
    assert!(centered_parser_over_strict > 0 && plusminus_parser_over_strict > 0, "raw streams must not be counted before parser cost");
    assert!(
        plusminus_active_quantum_gap_to_2700k > 50_000_000,
        "plus-minus active-chain quantum-control blocker changed; revisit physical integration"
    );
    assert!(
        direct_signnorm_rank_over_google > 0 && direct_signnorm_ambiguous_rank_over_google > 0,
        "sign-normalized direct route should stay blocked until normalization signs fit Google scratch"
    );
    assert!(
        direct_signnorm_exact_split_gap > 0,
        "phase-clean exact sign normalization should not be counted as p99 low-qubit solved"
    );
    assert!(
        direct_signnorm_logsign_once_gap > 0 && direct_signnorm_logsign_split_gap > 0,
        "logical coefficient signs now reach p99 low-qubit target; promote direct sign-normalized route"
    );
    assert!(
        direct_signnorm_mbu_degree_n14 + 1 >= 14
            && direct_signnorm_mbu_density_n14 > (1usize << 14) / 4
            && direct_signnorm_mbu_max_count_n14 > 4,
        "normalization-sign MBU parity changed; revisit sign-normalized direct route"
    );
    assert!(
        direct_signnorm_reverse_collisions_n14 > 2_000
            && direct_signnorm_reverse_states_n14 > 60_000
            && direct_signnorm_reverse_total_steps_n14 > 80_000,
        "normalization signs may be reverse-recoverable now; revisit sign-normalized direct route"
    );
    assert!(
        direct_restoring_final_raw_digit_over_strict > 0
            && direct_restoring_final_raw_digit_gap_google <= 0,
        "restoring-final direct route scratch state changed; update the Google-scratch ledger"
    );
    assert!(
        direct_restoring_final_select3x_gap < 0
            && direct_restoring_final_toy_neg2_cases > 0
            && direct_restoring_final_toy_zero_final_cases > 0,
        "restoring-final direct route lost its modeled low-qubit margin or toy coverage"
    );
    assert!(
        direct_restoring_final_bennett_fast_inverse_toy_ccx > 2 * direct_restoring_final_toy_ccx
            && direct_restoring_final_bennett_fast_inverse_toy_ccx < 700
            && direct_restoring_final_bennett_fast_inverse_toy_peak_q > direct_restoring_final_toy_peak_q,
        "restoring-final fast-inverse cleanup changed; revisit production packing budget"
    );
    assert!(
        direct_restoring_final_single_selector_toy_ccx < direct_restoring_final_toy_ccx
            && direct_restoring_final_single_selector_toy_peak_q
                > direct_restoring_final_toy_peak_q
            && direct_restoring_final_single_selector_bennett_toy_ccx
                < direct_restoring_final_bennett_fast_inverse_toy_ccx
            && direct_restoring_final_single_selector_bennett_toy_peak_q
                > direct_restoring_final_bennett_fast_inverse_toy_peak_q,
        "restoring-final select1 toy changed; revisit selector-factor ledger"
    );
    assert!(
        direct_restoring_final_payload_mbu_degree_n14 + 1 >= 14
            && direct_restoring_final_payload_mbu_density_n14 > (1usize << 14) / 4
            && direct_restoring_final_payload_max_n14 > 14,
        "restoring-final payload MBU toy result changed; revisit parser shortcut"
    );
    assert!(
        direct_restoring_final_reverse_q_collisions_n14 == 0
            && direct_restoring_final_reverse_q_states_n14
                == direct_restoring_final_reverse_q_total_steps_n14
            && direct_restoring_final_reverse_q_max_mult_n14 == 1,
        "restoring-final reverse q recovery changed; revisit no-payload decoder route"
    );
    assert!(
        direct_restoring_final_residual_q_collisions_n14 > 0
            && direct_restoring_final_residual_q_states_n14
                < direct_restoring_final_residual_q_total_steps_n14
            && direct_restoring_final_residual_q_max_mult_n14 > 1,
        "residual-only restoring-final reverse q is no longer ambiguous"
    );
    assert!(
        direct_restoring_final_reverse_coeff_candidates_transitions_n14
            == direct_restoring_final_reverse_coeff_candidates_exact_n14
                + direct_restoring_final_reverse_coeff_candidates_low_n14
                + direct_restoring_final_reverse_coeff_candidates_high_n14
            && direct_restoring_final_reverse_coeff_candidates_endpoints_n14
                == direct_restoring_final_reverse_coeff_candidates_exact_n14
            && direct_restoring_final_reverse_coeff_candidates_low_n14 > 0
            && direct_restoring_final_reverse_coeff_candidates_high_n14 > 0,
        "restoring-final coefficient candidate accounting changed; revisit branch selector"
    );
    assert!(
        direct_restoring_final_reverse_coeff_high_branch_degree_n14 + 1 >= 14
            && direct_restoring_final_reverse_coeff_high_branch_density_n14
                > (1usize << 14) / 4
            && direct_restoring_final_reverse_coeff_high_branch_max_count_n14 > 0
            && direct_restoring_final_reverse_coeff_high_branch_total_n14
                == direct_restoring_final_reverse_coeff_candidates_high_n14,
        "restoring-final coefficient high-branch selector stopped looking dense"
    );
    assert!(
        direct_restoring_final_coeff_decoder_exact_p99
            > direct_restoring_final_coeff_decoder_oneway_margin
            && direct_restoring_final_coeff_decoder_margin < 0
            && direct_restoring_final_coeff_decoder_augmented_pointadd_p99
                > GOOGLE_LOW_QUBIT_TOFFOLI
            && direct_restoring_final_coeff_decoder_augmented_gap > 0,
        "restoring-final coefficient decoder now fits; promote no-payload cleanup"
    );
    assert!(
        direct_restoring_final_avg_exact_select3_gap > 0
            && direct_restoring_final_avg_exact_select3_first64 > GOOGLE_LOW_QUBIT_TOFFOLI
            && direct_restoring_final_avg_exact_select3_p99 > GOOGLE_LOW_QUBIT_TOFFOLI
            && direct_restoring_final_avg_noscan_select1_gap < 0
            && direct_restoring_final_avg_noscan_select1_first64 < GOOGLE_LOW_QUBIT_TOFFOLI
            && direct_restoring_final_avg_noscan_select1_p99 < GOOGLE_LOW_QUBIT_TOFFOLI
            && direct_restoring_final_avg_noscan_select2_gap < 0
            && direct_restoring_final_avg_noscan_select2_first64 < GOOGLE_LOW_QUBIT_TOFFOLI
            && direct_restoring_final_avg_noscan_select2_p99 > GOOGLE_LOW_QUBIT_TOFFOLI
            && direct_restoring_final_avg_noscan_select3_gap < 0
            && direct_restoring_final_avg_noscan_select3_first64 < GOOGLE_LOW_QUBIT_TOFFOLI
            && direct_restoring_final_avg_noscan_select3_p99 > GOOGLE_LOW_QUBIT_TOFFOLI,
        "restoring-final average gate changed; revisit selector factor and scan-free decoder"
    );
    assert!(
        direct_restoring_final_stored_align_branch_select1_mean < GOOGLE_LOW_QUBIT_TOFFOLI
            && direct_restoring_final_stored_align_branch_select1_first64
                < GOOGLE_LOW_QUBIT_TOFFOLI
            && direct_restoring_final_stored_align_fixed_scratch_p99 > GOOGLE_LOW_QUBIT_SCRATCH
            && direct_restoring_final_stored_align_variable_scratch_p99 <= STRICT_SCRATCH + 2
            && direct_restoring_final_stored_align_delimited_scratch_p99
                > GOOGLE_LOW_QUBIT_SCRATCH
            && direct_restoring_final_stored_align_gamma_scratch_p99 > GOOGLE_LOW_QUBIT_SCRATCH,
        "restoring-final stored-alignment metadata changed; revisit exact parser/packing blocker"
    );
    assert!(
        direct_restoring_final_coeff_decoder_alignment_degree_n14 + 1 >= 14
            && direct_restoring_final_coeff_decoder_alignment_density_n14 > (1usize << 14) / 4
            && direct_restoring_final_coeff_decoder_alignment_max_n14 + 1 >= 14,
        "restoring-final coefficient decoder alignment metadata stopped looking dense"
    );
    assert!(halfgcd_tail_over_google > 0, "half-GCD checkpoint must be fused before it fits");
    assert!(
        halfgcd_det_compressed_tail_gap < 0 && halfgcd_det_recovery_num_bits_p99 > 256,
        "half-GCD determinant compression state changed; update recovery blocker"
    );
    assert!(
        halfgcd_tail_raw_rank_max_mult_n14 == 1
            && halfgcd_tail_raw_rank_degree_n14 == 0
            && halfgcd_tail_raw_rank_density_n14 == 0,
        "half-GCD raw-tail parser toy result changed; update frontier blocker"
    );
    assert!(
        halfgcd_tail_raw_compressed_rank_max_mult_n14 == 1
            && halfgcd_tail_raw_compressed_rank_degree_n14 == 0
            && halfgcd_tail_raw_compressed_rank_density_n14 == 0,
        "half-GCD compressed raw-tail parser toy result changed; update frontier blocker"
    );
    assert!(
        halfgcd_replay_with_recovery_floor_gap_to_2700k < 0,
        "half-GCD arithmetic replay floor changed; update matrix-extraction blocker"
    );
    assert!(
        halfgcd_full_prefix_live_gap_google > 0 && halfgcd_compressed_tail_stream_peak_gap_google <= 0,
        "half-GCD checkpoint extraction schedule changed; update prefix-compression blocker"
    );
    assert!(
        halfgcd_inloop_recovery_gap_to_2700k > 0,
        "half-GCD in-loop determinant recovery floor no longer blocks; revisit compressed checkpoint route"
    );
    assert!(
        halfgcd_second_col_tail_stream_peak_gap_google < 0
            && halfgcd_second_col_tail_raw_rank_max_mult_n14 == 1
            && halfgcd_second_col_tail_raw_rank_degree_n14 == 0
            && halfgcd_second_col_tail_raw_rank_density_n14 == 0
            && halfgcd_second_col_prefix_final_bd_max_mult_n14 == 1
            && halfgcd_second_col_prefix_local_reverse_max_mult_n14 == 1
            && halfgcd_second_col_prefix_local_reverse_collisions_n14 == 0
            && halfgcd_second_col_prefix_reverse_formula_transitions_n14 == halfgcd_second_col_prefix_transitions_n14
            && halfgcd_second_col_prefix_reverse_formula_endpoints_n14 == 16_380
            && halfgcd_second_col_prefix_reverse_formula_coeff_steps_n14 == 65_648,
        "half-GCD second-column tail-stream scratch/parser changed; revisit prefix extraction route"
    );
    assert!(
        halfgcd_second_col_prefix_residual_q_collisions_n14 > 0
            && halfgcd_second_col_prefix_residual_q_states_n14
                < halfgcd_second_col_prefix_residual_q_total_steps_n14
            && halfgcd_second_col_prefix_residual_q_max_mult_n14 > 1,
        "half-GCD second-column residual-only reverse q is no longer ambiguous"
    );
    assert!(
        halfgcd_second_col_prefix_exact_extraction_p99 < halfgcd_second_col_prefix_oneway_budget_ccx
            && halfgcd_second_col_prefix_exact_gap_to_2700k < 0,
        "half-GCD second-column optimistic prefix ledger no longer has low-qubit margin"
    );
    assert!(
        halfgcd_second_col_prefix_augmented_extraction_p99
            > halfgcd_second_col_prefix_oneway_budget_ccx
            && halfgcd_second_col_prefix_augmented_gap_to_2700k > 0,
        "half-GCD coefficient q-cleanup no longer blocks the exact-prefix route"
    );
    assert!(
        halfgcd_second_col_prefix_coeff_decoder_no_scan_p99
            < halfgcd_second_col_prefix_oneway_budget_ccx
                - halfgcd_second_col_prefix_exact_extraction_p99
            && halfgcd_second_col_prefix_coeff_decoder_scan_budget < 4_000
            && halfgcd_second_col_prefix_coeff_decoder_scan_over_budget > 40_000,
        "half-GCD coefficient decoder alignment budget changed; revisit q-cleanup route"
    );
    assert!(
        halfgcd_second_col_prefix_avg_aug_exact_gap < 0
            && halfgcd_second_col_prefix_avg_aug_exact_first64 < GOOGLE_LOW_QUBIT_TOFFOLI
            && halfgcd_second_col_prefix_avg_aug_exact_p99 > GOOGLE_LOW_QUBIT_TOFFOLI,
        "half-GCD exact-decoder average route changed; revisit implementation priority"
    );
    assert!(
        halfgcd_second_col_prefix_avg_aug_noscan_gap < 0
            && halfgcd_second_col_prefix_avg_aug_noscan_p99 < GOOGLE_LOW_QUBIT_TOFFOLI,
        "half-GCD scan-free lower bound no longer has full sampled margin"
    );
    assert!(
        halfgcd_second_col_prefix_step_toy_ccx == 308
            && halfgcd_second_col_prefix_step_toy_peak_q > 100
            && halfgcd_second_col_prefix_step_toy_final_negative_cases > 10_000,
        "half-GCD prefix-step toy evidence changed; revisit extractor implementation risk"
    );
    assert!(
        halfgcd_second_col_prefix_fixed_bound_active_toy_ccx == 11_464
            && halfgcd_second_col_prefix_fixed_bound_active_toy_peak_q == 173
            && halfgcd_second_col_prefix_fixed_bound_active_toy_active_slots > 0
            && halfgcd_second_col_prefix_fixed_bound_active_toy_inactive_slots > 0
            && halfgcd_second_col_prefix_fixed_bound_active_toy_halted_inputs > 0
            && halfgcd_second_col_prefix_fixed_bound_active_toy_full_bound_inputs > 0
            && halfgcd_second_col_prefix_fixed_bound_active_toy_dirty_restore_cases == 0
            && halfgcd_second_col_prefix_fixed_bound_active_toy_dirty_history_cases == 0
            && halfgcd_second_col_prefix_fixed_bound_active_toy_dirty_phase_cases == 0,
        "half-GCD fixed-bound active prefix toy changed; revisit variable-length prefix implementation risk"
    );
    assert!(
        halfgcd_second_col_prefix_active_model_base_mean < GOOGLE_LOW_QUBIT_TOFFOLI
            && halfgcd_second_col_prefix_active_model_pointadd_mean > GOOGLE_LOW_QUBIT_TOFFOLI
            && halfgcd_second_col_prefix_active_model_pointadd_first64 > GOOGLE_LOW_QUBIT_TOFFOLI
            && halfgcd_second_col_prefix_active_model_pointadd_p99 > GOOGLE_LOW_QUBIT_TOFFOLI
            && halfgcd_second_col_prefix_active_model_gap_to_2700k > 700_000
            && halfgcd_second_col_prefix_active_model_over_exact_mean > 1_000_000
            && halfgcd_second_col_prefix_active_model_over_exact_p99 > 1_000_000,
        "half-GCD active-control prefix model no longer blocks the Bennett cleanup route"
    );
}
