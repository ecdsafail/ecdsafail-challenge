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
            blocker: "second-column live/tail state fits; optimistic exact-barrel prefix ledger projects 2385096, but phase-clean extractor implementation is still unbuilt",
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
    let direct_signnorm_mbu_degree_n14 = 13usize;
    let direct_signnorm_mbu_density_n14 = 8_208usize;
    let direct_signnorm_mbu_max_count_n14 = 8usize;
    let direct_signnorm_reverse_collisions_n14 = 2_658usize;
    let direct_signnorm_reverse_states_n14 = 64_178usize;
    let direct_signnorm_reverse_total_steps_n14 = 89_008usize;
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
    let halfgcd_second_col_prefix_steps_p99 = 91usize;
    let halfgcd_second_col_prefix_digits_p99 = 253usize;
    let halfgcd_second_col_prefix_bounded_barrel_bits = 5usize;
    let halfgcd_second_col_prefix_residual_digit_width_p99 = 48_776usize;
    let halfgcd_second_col_prefix_coeff_digit_width_p99 = 17_316usize;
    let halfgcd_second_col_prefix_oneway_budget_ccx = 345_059usize;
    let halfgcd_second_col_prefix_bounded_extraction_p99 = 206_323usize;
    let halfgcd_second_col_prefix_exact_extraction_p99 = 276_840usize;
    let halfgcd_second_col_prefix_exact_pointadd_p99 = 2_385_096usize;
    let halfgcd_second_col_prefix_exact_gap_to_2700k =
        halfgcd_second_col_prefix_exact_pointadd_p99 as isize - GOOGLE_LOW_QUBIT_TOFFOLI as isize;

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
    println!("METRIC scratch600_direct_signnorm_mbu_degree_n14={direct_signnorm_mbu_degree_n14}");
    println!("METRIC scratch600_direct_signnorm_mbu_density_n14={direct_signnorm_mbu_density_n14}");
    println!("METRIC scratch600_direct_signnorm_mbu_max_count_n14={direct_signnorm_mbu_max_count_n14}");
    println!("METRIC scratch600_direct_signnorm_reverse_collisions_n14={direct_signnorm_reverse_collisions_n14}");
    println!("METRIC scratch600_direct_signnorm_reverse_states_n14={direct_signnorm_reverse_states_n14}");
    println!("METRIC scratch600_direct_signnorm_reverse_total_steps_n14={direct_signnorm_reverse_total_steps_n14}");
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
    println!("METRIC scratch600_halfgcd_second_col_prefix_steps_p99={halfgcd_second_col_prefix_steps_p99}");
    println!("METRIC scratch600_halfgcd_second_col_prefix_digits_p99={halfgcd_second_col_prefix_digits_p99}");
    println!("METRIC scratch600_halfgcd_second_col_prefix_bounded_barrel_bits={halfgcd_second_col_prefix_bounded_barrel_bits}");
    println!("METRIC scratch600_halfgcd_second_col_prefix_residual_digit_width_p99={halfgcd_second_col_prefix_residual_digit_width_p99}");
    println!("METRIC scratch600_halfgcd_second_col_prefix_coeff_digit_width_p99={halfgcd_second_col_prefix_coeff_digit_width_p99}");
    println!("METRIC scratch600_halfgcd_second_col_prefix_oneway_budget_ccx={halfgcd_second_col_prefix_oneway_budget_ccx}");
    println!("METRIC scratch600_halfgcd_second_col_prefix_bounded_extraction_p99={halfgcd_second_col_prefix_bounded_extraction_p99}");
    println!("METRIC scratch600_halfgcd_second_col_prefix_exact_extraction_p99={halfgcd_second_col_prefix_exact_extraction_p99}");
    println!("METRIC scratch600_halfgcd_second_col_prefix_exact_pointadd_p99={halfgcd_second_col_prefix_exact_pointadd_p99}");
    println!("METRIC scratch600_halfgcd_second_col_prefix_exact_gap_to_2700k={halfgcd_second_col_prefix_exact_gap_to_2700k}");

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
            && halfgcd_second_col_tail_raw_rank_density_n14 == 0,
        "half-GCD second-column tail-stream scratch/parser changed; revisit prefix extraction route"
    );
    assert!(
        halfgcd_second_col_prefix_exact_extraction_p99 < halfgcd_second_col_prefix_oneway_budget_ccx
            && halfgcd_second_col_prefix_exact_gap_to_2700k < 0,
        "half-GCD second-column optimistic prefix ledger no longer has low-qubit margin"
    );
}
