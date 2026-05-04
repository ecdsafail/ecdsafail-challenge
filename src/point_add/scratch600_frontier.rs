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
            name: "by_centered_exactparity_fast_signed_clean_replay",
            scratch_bits: 2_200,
            charged_toffoli: Some(5_878_716),
            blocker: "one clean 560-step replay body is 2618000 CCX at 2720q; exact parity cleanup alone exceeds the per-DIV budget and two clean DIVs project 3178716 over target",
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
            name: "scaled_by_compressed_pattern_fixed_id_decode",
            scratch_bits: 600,
            charged_toffoli: Some(2_944_889),
            blocker: "two fast scaled-BY replays project 2637286 before compressed-ID expansion; sampled fixed-ID row-touch floor adds 307603 and misses by 244889, so this needs a structured pattern parser",
        },
        Candidate {
            name: "scaled_by_raw_pattern_streaming_parser",
            scratch_bits: 671,
            charged_toffoli: Some(2_701_606),
            blocker: "raw 560-bit pattern plus single A only fits 663 scratch if the delta parser is non-reversible; sampled reversible delta checkpoint needs 5 bits and 666 scratch, retained A history is p99 218 bits, and two exact clean pattern decoders miss by 1606 before compressed expansion",
        },
        Candidate {
            name: "scaled_by_h_only_compressed_history_budget",
            scratch_bits: 2_224,
            charged_toffoli: None,
            blocker: "h-only compressed-history model is sub-MToffoli and under 2800q on sampled modular arithmetic, but exact W=4 toys show next-h update payload saturates a full rank per window (n14 rank_p99=28 over 7 windows, max_next_h_choices=16), so it still needs charged high-ratio/rank history before integration",
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
            name: "plusminus_scaled_solinas_history_scale_packed",
            scratch_bits: 822,
            charged_toffoli: Some(2_230_850),
            blocker: "optimistic Solinas history-carry scale model is gate-shaped, but the actual k22 split multihalve chunk peaks at 822q; even granting one reusable 256-bit lane remains 159 scratch bits over Google, and naive overlap is 1078 scratch",
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
            blocker: "raw sign-normalized digits fit, but exact cneg p99 is 2792914; norm signs have dense MBU parity and magnitude-only exact toy reverse collisions. Post-step rows disambiguate signs, and det-low2 xor coeff_v_sign recovers the norm sign on exact toys (n14 formula mismatches 0/89008), but raw physical cneg remains too expensive",
        },
        Candidate {
            name: "direct_centered_signnorm_logical_coeff_signs",
            scratch_bits: 657,
            charged_toffoli: Some(2_723_992),
            blocker: "det-low2 xor coeff_v_sign removes the normalization-sign sidecar in exact toys, and the local predicate toy is phase-clean at 14 CCX. Deleting physical rem cneg would clear p99 by 2720, but the tested signed-remainder recurrence jumps to 3136080 from 180 p99 steps, and the signed-domain non-restoring body needs a relative-negative zero predicate per digit (toy 416 CCX, scaled predicate 1025 CCX), so there is no simple two's-complement signed-domain wire-in",
        },
        Candidate {
            name: "direct_centered_restoring_final_stored_alignment",
            scratch_bits: 662,
            charged_toffoli: Some(2_709_483),
            blocker: "restoring-final select1 has phase-clean toy cleanup; 7-symbol branch-conditioned blocks fit 662 scratch and lower-bound to 2655117 average, but a 2x binary compare/subtract parser floor still pushes 9483 over with 52 non-contiguous alignment-support steps before rank mapping, renormalization, or cleanup",
        },
        Candidate {
            name: "direct_centered_restoring_final_mixed67_parser",
            scratch_bits: 663,
            charged_toffoli: Some(2_708_727),
            blocker: "period-5 mask-9 mixed 6/7 branch-conditioned blocks fit 663 scratch and reduce the 2x binary lookup miss to 8727 Toffoli, but still require sub-1.679x lookup cleanup to reach target",
        },
        Candidate {
            name: "direct_centered_restoring_final_mixed4to8_parser",
            scratch_bits: 663,
            charged_toffoli: Some(2_708_680),
            blocker: "period-4 code-8656 mixed 4..8 branch-conditioned blocks fit 663 scratch and reduce the 2x binary lookup miss to 8680 Toffoli, but still require sub-1.681x lookup cleanup; selective adjacent-pair grouping saves only 26.9 of 1084.9 needed, and coherent full-tree lookup still misses",
        },
        Candidate {
            name: "direct_centered_restoring_final_low_branch_digit_mixed4to8_floor",
            scratch_bits: 663,
            charged_toffoli: Some(2_643_614),
            blocker: "low-candidate branch-as-final-digit lower bound clears binary lookup by 56386, high_q=low_q+1 on the sample set, and a 23-CCX branch digit toy is Bennett-clean; superseded by the alignment-only parser floor, but still needs full reverse-decoder integration",
        },
        Candidate {
            name: "direct_centered_restoring_final_low_branch_align_only_prefix_tree_floor",
            scratch_bits: 580,
            charged_toffoli: None,
            blocker: "branch-as-final-digit removes branch symbols from the parser stream; low-alignment block2 fits 580 scratch and prefix-tree node floor projects 2593870, but this raw row excludes support-weighted selected add/sub and variable-support decoder integration",
        },
        Candidate {
            name: "direct_centered_restoring_final_low_branch_weighted_prefix_span_floor",
            scratch_bits: 580,
            charged_toffoli: None,
            blocker: "support-weighted selected add/sub raises the low-branch prefix projection to 2666583 with 33417 margin, but a Shannon-style variable block2 decoder misses by 6303 and offset-1 pairing still misses by 5674; superseded by selective length-flattening at 663 scratch",
        },
        Candidate {
            name: "direct_centered_restoring_final_low_branch_selective_prefix_flatten_floor",
            scratch_bits: 663,
            charged_toffoli: None,
            blocker: "p99-only selective length-flattening has 394 sampled max prefix bits and would need 676 scratch; trimming 9 balanced steps gives a 381-bit sampled max, fits 663 scratch, and projects 2661534 with 38466 margin; support-2..18 generated balanced block2 selected-add/sub roundtrip family is phase-clean across 289 pairs with max 804 CCX, and peak-fit mixed schedule codebooks decode 856854 sampled symbols with no collisions or mismatches; disjoint 8192 secp holdout already has 182 missing symbols, 170 missing traces, 7 over-budget rows, and 391 seen bits; a parity scaling probe at 65536 training / 32768 holdout still leaves 105 missing symbols and only about 631 Toffoli margin; toy exact-domain train/exhaust probes miss symbols in all 4 cases, support-only toys find modest but real misses (n16: 26 symbols over 11 steps, exact contiguous span 16), and charged interval probes do not rescue it: guard4 covers 0/4 toy domains while a full 0..n per-step envelope covers 4/4 but fits 0/4 and reaches 4268 over-budget traces / 37 bits vs n16 budget 24; promotion now needs a non-sampled support proof or a different decoder, not the simple interval fallback",
        },
        Candidate {
            name: "direct_centered_restoring_final_mixed4to8_joint_binary_floor",
            scratch_bits: 663,
            charged_toffoli: Some(2_693_369),
            blocker: "joint block-pattern binary-depth floor would clear 2.7M by 6631 at 663 scratch, but assumes a phase-clean block-rank decoder; exact n14 rank parity is degree 14 and 8098/16384 dense, selective adjacent-pair grouping saves only 26.9 of 1084.9 needed, and arbitrary full-scan support is 68058 rows and misses by 498777",
        },
        Candidate {
            name: "direct_centered_restoring_final_mixed67_huffman_floor",
            scratch_bits: 663,
            charged_toffoli: Some(2_690_447),
            blocker: "distribution-aware Huffman path floor would clear 2.7M by 9553 at 663 scratch, but coherent tree execution reverts to the full scan and misses by 105208 unless a phase-clean classical-path decoder is found",
        },
        Candidate {
            name: "direct_centered_restoring_final_cond_block6_parser",
            scratch_bits: 665,
            charged_toffoli: Some(2_708_135),
            blocker: "6-symbol branch-conditioned blocks cut the 2x binary lookup miss to 8135 Toffoli, but require 665 scratch p99, 2 over the Google scratch model",
        },
        Candidate {
            name: "direct_centered_restoring_final_block4_parser",
            scratch_bits: 675,
            charged_toffoli: Some(2_705_475),
            blocker: "4-symbol parser blocks cut the 2x binary lookup miss to 5475 Toffoli but require 675 scratch p99, 12 over the Google scratch model",
        },
        Candidate {
            name: "direct_centered_signnorm_rank_compressed_signs",
            scratch_bits: 765,
            charged_toffoli: None,
            blocker: "superseded by det-low2 coefficient-sign recovery; rank-compressed normalization signs still document that generic sign sidecars need 765 p99 scratch bits",
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
            name: "halfgcd_second_column_fixed_depth64_dynamic_barrel_model",
            scratch_bits: 515,
            charged_toffoli: None,
            blocker: "if alignment layers are BitId conditions, depth64 dynamic barrels average 1986713 with p99 2047416, but simulator stats do not discount quantum controls; HMR controls are random rather than alignment values and generic alignment-control MBU phase is dense at n14",
        },
        Candidate {
            name: "halfgcd_second_column_fixed_depth64_tail_stream",
            scratch_bits: 515,
            charged_toffoli: Some(2_934_322),
            blocker: "fixed-depth64 popcount-priced coefficient application averages 2740052, but coefficient bits are quantum data; a generous static binary application floor averages 2934322, so this route needs classical coefficient controls or a static window recode",
        },
        Candidate {
            name: "halfgcd_second_column_fixed_depth64_static_window_floor",
            scratch_bits: 515,
            charged_toffoli: Some(2_748_271),
            blocker: "joint static-window scan improves to w6 average 2749506 (+49506) under the exact bit-product floor; sparse signed wNAF recoding lowers the source-product floor to 2748271 (+48271), but still needs selector/recoder cost below 86824 one-way instead of 99575; free-active compact NAF w2 would clear at 2691392, but the omitted active/zero predicate is 38097 one-way against 4304 slack, with active-only toy parity dense at n14 (degree 14, 8322/16384) and active support 29/30 slots; table-only w4 would be 2559198 before data application, but row-controlled source products make the best table-source floor w2 average 3956644, generic cleanup is dense at n14 (plain 8194/16384, wNAF 8162/16384), and exact toy support leaves 27/28 coefficient bit positions live",
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
    let scaled_by_pattern_fixed_id_bits = 481usize;
    let scaled_by_pattern_fixed_id_distinct_rows = 307_603usize;
    let scaled_by_pattern_fixed_id_max_window_rows = 9_339usize;
    let scaled_by_pattern_fixed_id_nonzero_table_bits = 2_457_030usize;
    let scaled_by_pattern_fixed_id_two_replay_before_decode = 2_637_286usize;
    let scaled_by_pattern_fixed_id_remaining_to_2700k =
        GOOGLE_LOW_QUBIT_TOFFOLI - scaled_by_pattern_fixed_id_two_replay_before_decode;
    let scaled_by_pattern_fixed_id_row_floor_gap =
        (scaled_by_pattern_fixed_id_two_replay_before_decode
            + scaled_by_pattern_fixed_id_distinct_rows) as isize
            - GOOGLE_LOW_QUBIT_TOFFOLI as isize;
    let scaled_by_pattern_fixed_id_bit_floor_gap =
        (scaled_by_pattern_fixed_id_two_replay_before_decode
            + scaled_by_pattern_fixed_id_nonzero_table_bits) as isize
            - GOOGLE_LOW_QUBIT_TOFFOLI as isize;
    let scaled_by_raw_pattern_bits = 560usize;
    let scaled_by_raw_pattern_delta_bits = 10usize;
    let scaled_by_raw_pattern_single_a_scratch = 661usize;
    let scaled_by_raw_pattern_one_checkpoint_scratch = 671usize;
    let scaled_by_raw_pattern_window_a_scratch = 676usize;
    let scaled_by_raw_pattern_delta_checkpoint_max_rows = 25usize;
    let scaled_by_raw_pattern_delta_checkpoint_bits = 5usize;
    let scaled_by_raw_pattern_delta_checkpoint_scratch = 666usize;
    let scaled_by_raw_pattern_delta_checkpoint_scratch_slack = 2usize;
    let scaled_by_raw_pattern_ambiguous_a_bits_mean_milli = 199_653usize;
    let scaled_by_raw_pattern_ambiguous_a_bits_p99 = 218usize;
    let scaled_by_raw_pattern_ambiguous_a_bits_max = 232usize;
    let scaled_by_raw_pattern_two_replay_before_branch_decode = 2_577_286usize;
    let scaled_by_raw_pattern_exact_decoder_per_replay = 62_160usize;
    let scaled_by_raw_pattern_exact_two_decoder_projection =
        scaled_by_raw_pattern_two_replay_before_branch_decode
            + 2 * scaled_by_raw_pattern_exact_decoder_per_replay;
    let scaled_by_raw_pattern_exact_two_decoder_gap =
        scaled_by_raw_pattern_exact_two_decoder_projection as isize
            - GOOGLE_LOW_QUBIT_TOFFOLI as isize;
    let scaled_by_h_only_model_modular_windows = 35usize;
    let scaled_by_h_only_model_modular_toffoli = 672_650usize;
    let scaled_by_h_only_model_peak = 2_736usize;
    let scaled_by_h_only_model_history_bits = 480usize;
    let scaled_by_h_only_next_ratio_toy_n14_windows = 7usize;
    let scaled_by_h_only_next_ratio_toy_n14_keys = 802usize;
    let scaled_by_h_only_next_ratio_toy_n14_ambiguous_keys = 742usize;
    let scaled_by_h_only_next_ratio_toy_n14_max_next_h_choices = 16usize;
    let scaled_by_h_only_next_ratio_toy_n14_rank_p99 = 28usize;
    let scaled_by_h_only_next_ratio_toy_n14_rank_max = 28usize;
    let scaled_by_h_only_next_ratio_toy_n14_rank_mean_milli = 27_718usize;
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
    let by_centered_exactparity_clean_replay_ccx = 2_618_000usize;
    let by_centered_exactparity_clean_peak_q = 2_720usize;
    let by_centered_exactparity_clean_scratch_bits = 2_200usize;
    let by_centered_exactparity_clean_per_div_budget =
        (GOOGLE_LOW_QUBIT_TOFFOLI - 642_716usize) / 2;
    let by_centered_exactparity_two_clean_div_projection = 5_878_716usize;
    let by_centered_exactparity_two_clean_div_gap =
        by_centered_exactparity_two_clean_div_projection as isize
            - GOOGLE_LOW_QUBIT_TOFFOLI as isize;
    let centered_raw_scratch = 592usize;
    let centered_boundary_scratch_p99 = 710usize;
    let centered_parser_over_strict = centered_boundary_scratch_p99 - STRICT_SCRATCH;
    let direct_signnorm_raw_digit_scratch_p99 = 653usize;
    let direct_signnorm_det_coeffsign_scratch_p99 = direct_signnorm_raw_digit_scratch_p99 + 4usize;
    let direct_signnorm_det_coeffsign_scratch_gap_google =
        direct_signnorm_det_coeffsign_scratch_p99 as isize - GOOGLE_LOW_QUBIT_SCRATCH as isize;
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
    let direct_signnorm_logsign_no_rem_cneg_projection_p99 = 2_697_280usize;
    let direct_signnorm_logsign_no_rem_cneg_gap =
        direct_signnorm_logsign_no_rem_cneg_projection_p99 as isize
            - GOOGLE_LOW_QUBIT_TOFFOLI as isize;
    let direct_signnorm_prefinal_signed_remainder_p99 = 3_136_080usize;
    let direct_signnorm_prefinal_signed_remainder_gap =
        direct_signnorm_prefinal_signed_remainder_p99 as isize
            - GOOGLE_LOW_QUBIT_TOFFOLI as isize;
    let direct_signnorm_prefinal_signed_remainder_count_p99 = 180usize;
    let direct_signnorm_prefinal_signed_remainder_digit_payload_p99 = 498usize;
    let direct_signnorm_prefinal_signed_remainder_width_extra_max = 1usize;
    let direct_signnorm_logsign_direct_rem_toy_ccx = 148usize;
    let direct_signnorm_logsign_direct_rem_toy_peak_q = 80usize;
    let direct_signnorm_logsign_direct_rem_toy_phase_dirty_cases = 0usize;
    let direct_signnorm_logsign_exact_cneg257 = 512usize;
    let direct_signnorm_logsign_exact_rem_p99 = 26_712usize;
    let direct_signnorm_logsign_exact_once_p99 = 2_746_960usize;
    let direct_signnorm_logsign_exact_split_p99 = 2_794_228usize;
    let direct_signnorm_logsign_once_gap =
        direct_signnorm_logsign_once_p99 as isize - GOOGLE_LOW_QUBIT_TOFFOLI as isize;
    let direct_signnorm_logsign_split_gap =
        direct_signnorm_logsign_split_p99 as isize - GOOGLE_LOW_QUBIT_TOFFOLI as isize;
    let direct_signnorm_logsign_exact_once_gap =
        direct_signnorm_logsign_exact_once_p99 as isize - GOOGLE_LOW_QUBIT_TOFFOLI as isize;
    let direct_signnorm_logsign_exact_split_gap =
        direct_signnorm_logsign_exact_split_p99 as isize - GOOGLE_LOW_QUBIT_TOFFOLI as isize;
    let direct_signnorm_mbu_degree_n14 = 13usize;
    let direct_signnorm_mbu_density_n14 = 8_208usize;
    let direct_signnorm_mbu_max_count_n14 = 8usize;
    let direct_signnorm_reverse_collisions_n14 = 2_658usize;
    let direct_signnorm_reverse_states_n14 = 64_178usize;
    let direct_signnorm_reverse_total_steps_n14 = 89_008usize;
    let direct_signnorm_coeff_reverse_collisions_n14 = 0usize;
    let direct_signnorm_coeff_reverse_states_n14 = 89_008usize;
    let direct_signnorm_coeff_reverse_total_steps_n14 = 89_008usize;
    let direct_signnorm_coeff_reverse_max_mult_n14 = 1usize;
    let direct_signnorm_coeff_reverse_zero_coeff_cases_n14 = 0usize;
    let direct_signnorm_det_sign_reverse_collisions_n14 = 2_654usize;
    let direct_signnorm_det_sign_reverse_states_n14 = 70_742usize;
    let direct_signnorm_det_sign_reverse_max_mult_n14 = 2usize;
    let direct_signnorm_det_coeffsign_reverse_collisions_n14 = 0usize;
    let direct_signnorm_det_coeffsign_reverse_states_n14 = 73_396usize;
    let direct_signnorm_det_coeffsign_reverse_total_steps_n14 = 89_008usize;
    let direct_signnorm_det_coeffsign_reverse_max_mult_n14 = 1usize;
    let direct_signnorm_det_coeffsign_bad_det_cases_n14 = 0usize;
    let direct_signnorm_det_coeffsign_low2_mismatches_n14 = 0usize;
    let direct_signnorm_det_coeffsign_formula_mismatches_n14 = 0usize;
    let direct_signnorm_det_coeffsign_predicate_p1_ccx = 14usize;
    let direct_signnorm_det_coeffsign_predicate_p1_peak_q = 18usize;
    let direct_signnorm_det_coeffsign_predicate_p1_valid_odd_det_cases = 3_072usize;
    let direct_signnorm_det_coeffsign_predicate_p3_ccx = 14usize;
    let direct_signnorm_det_coeffsign_predicate_p3_peak_q = 18usize;
    let direct_signnorm_det_coeffsign_predicate_p3_valid_odd_det_cases = 3_072usize;
    let direct_signnorm_signed_domain_relative_negative_toy_ccx = 45usize;
    let direct_signnorm_signed_domain_relative_negative_257_ccx = 1_025usize;
    let direct_signnorm_signed_domain_floor_toy_ccx = 416usize;
    let direct_signnorm_signed_domain_floor_toy_peak_q = 62usize;
    let direct_signnorm_signed_domain_floor_toy_final_negative_cases = 1_984usize;
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
    let direct_restoring_final_branch_digit_toy_branch_ccx = 23usize;
    let direct_restoring_final_branch_digit_toy_forward_ccx = 321usize;
    let direct_restoring_final_branch_digit_toy_roundtrip_ccx = 642usize;
    let direct_restoring_final_branch_digit_toy_peak_q = 109usize;
    let direct_restoring_final_branch_digit_toy_branch_one_cases = 40_425usize;
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
    let direct_restoring_final_low_branch_adjacent_transitions_n14 = 105_388usize;
    let direct_restoring_final_low_branch_adjacent_ambiguous_n14 = 89_008usize;
    let direct_restoring_final_low_branch_adjacent_high_n14 = 40_112usize;
    let direct_restoring_final_low_branch_adjacent_violations_n14 = 0usize;
    let direct_restoring_final_low_branch_adjacent_max_delta_n14 = 1usize;
    let direct_restoring_final_low_branch_adjacent_max_alignment_n14 = 13usize;
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
    let direct_restoring_final_stored_align_length_rank_scratch_p99 = 748usize;
    let direct_restoring_final_stored_align_length_rank_scratch_max = 769usize;
    let direct_restoring_final_stored_align_public_len_mismatch_p99 = 50usize;
    let direct_restoring_final_stored_align_public_len_rank_p99 = 109usize;
    let direct_restoring_final_stored_align_public_len_position_only_p99 = 706usize;
    let direct_restoring_final_stored_align_public_len_position_plus3_p99 = 840usize;
    let direct_restoring_final_stored_align_q_len_position_only_p99 = 695usize;
    let direct_restoring_final_stored_align_digit_len_position_only_p99 = 686usize;
    let direct_restoring_final_stored_align_joint_len_position_only_p99 = 681usize;
    let direct_restoring_final_stored_align_joint_len_position_plus3_p99 = 755usize;
    let direct_restoring_final_stored_align_pop_barrel_p99 = 19_928usize;
    let direct_restoring_final_stored_align_branch_select_p99 = 31_033usize;
    let direct_restoring_final_stored_align_branch_count_p99 = 117usize;
    let direct_restoring_final_branch_final_current_branch_select_mean = 26_960.072f64;
    let direct_restoring_final_branch_final_current_branch_select_p99 = 31_033usize;
    let direct_restoring_final_branch_final_width_mean = 13_531.836f64;
    let direct_restoring_final_branch_final_width_p99 = 15_573usize;
    let direct_restoring_final_branch_final_width_minus1_mean = 13_428.236f64;
    let direct_restoring_final_branch_final_selected_width_saving_mean = 13_428.236f64;
    let direct_restoring_final_branch_final_low_path_width_saving_mean = 16_266.246f64;
    let direct_restoring_final_branch_final_low_path_width_minus1_saving_mean = 16_369.846f64;
    let direct_restoring_final_branch_final_low_extra_touch_mean = -2_838.010f64;
    let direct_restoring_final_branch_final_branch_count_mean = 103.600f64;
    let direct_restoring_final_branch_final_branch_count_p99 = 117usize;
    let direct_restoring_final_branch_final_high_branch_mean = 45.614f64;
    let direct_restoring_final_branch_final_high_branch_p99 = 59usize;
    let direct_restoring_final_branch_final_alignment_diff_p99 = 26usize;
    let direct_restoring_final_branch_final_digit_len_diff_p99 = 38usize;
    let direct_restoring_final_branch_final_high_adjacent_violations = 0usize;
    let direct_restoring_final_branch_final_current_mixed4to8_gap = 8_679.472f64;
    let direct_restoring_final_branch_final_current_scan_mixed4to8_gap = 105_160.608f64;
    let direct_restoring_final_branch_final_selected_width_lookup_target_mean = 12_424.944f64;
    let direct_restoring_final_branch_final_low_path_width_lookup_target_mean = 13_843.949f64;
    let direct_restoring_final_branch_final_low_path_width_minus1_lookup_target_mean =
        13_895.749f64;
    let direct_restoring_final_branch_final_selected_width_lookup_multiplier_budget =
        1.828_338f64;
    let direct_restoring_final_branch_final_low_path_width_lookup_multiplier_budget =
        2.037_145f64;
    let direct_restoring_final_branch_final_low_path_width_minus1_lookup_multiplier_budget =
        2.044_768f64;
    let direct_restoring_final_branch_final_selected_width_mixed4to8_gap = -45_033.471f64;
    let direct_restoring_final_branch_final_low_path_width_mixed4to8_gap = -56_385.513f64;
    let direct_restoring_final_branch_final_low_path_width_minus1_mixed4to8_gap = -56_799.914f64;
    let direct_restoring_final_branch_final_selected_width_scan_mixed4to8_gap = 51_447.665f64;
    let direct_restoring_final_branch_final_low_path_width_scan_mixed4to8_gap = 40_095.623f64;
    let direct_restoring_final_branch_final_low_path_width_minus1_scan_mixed4to8_gap =
        39_681.222f64;
    let direct_restoring_final_low_branch_align_only_low_path_width_saving_mean =
        16_267.701f64;
    let direct_restoring_final_low_branch_align_only_model_precision_bits = 13usize;
    let direct_restoring_final_low_branch_align_only_raw_scratch_p99 = 471usize;
    let direct_restoring_final_low_branch_align_only_raw_scratch_max = 478usize;
    let direct_restoring_final_low_branch_align_only_step_entropy_scratch_p99 = 530usize;
    let direct_restoring_final_low_branch_align_only_step_entropy_scratch_max = 547usize;
    let direct_restoring_final_low_branch_align_only_step_prefix_scratch_p99 = 578usize;
    let direct_restoring_final_low_branch_align_only_step_prefix_scratch_max = 593usize;
    let direct_restoring_final_low_branch_align_only_best_block_symbols = 2usize;
    let direct_restoring_final_low_branch_align_only_best_touch_floor_mean = 542.607f64;
    let direct_restoring_final_low_branch_align_only_best_touch_floor_p99 = 594usize;
    let direct_restoring_final_low_branch_align_only_best_compressed_bits_p99 = 298usize;
    let direct_restoring_final_low_branch_align_only_best_live_scratch_p99 = 580usize;
    let direct_restoring_final_low_branch_align_only_best_symbol_count_p99 = 118usize;
    let direct_restoring_final_low_branch_align_only_best_augmented_gap = -117_630.377f64;
    let direct_restoring_final_low_branch_align_only_mixed4to8_schedule_code = 44usize;
    let direct_restoring_final_low_branch_align_only_mixed4to8_touch_floor_mean = 1_044.569f64;
    let direct_restoring_final_low_branch_align_only_mixed4to8_live_scratch_p99 = 568usize;
    let direct_restoring_final_low_branch_align_only_scan_lookup_floor_mean = 18_687.901f64;
    let direct_restoring_final_low_branch_align_only_scan_lookup_floor_p99 = 20_384usize;
    let direct_restoring_final_low_branch_align_only_binary_lookup_floor_mean = 5_541.966f64;
    let direct_restoring_final_low_branch_align_only_binary_lookup_floor_p99 = 6_214usize;
    let direct_restoring_final_low_branch_align_only_huffman_lookup_floor_mean = 3_336.980f64;
    let direct_restoring_final_low_branch_align_only_huffman_lookup_floor_p99 = 3_666usize;
    let direct_restoring_final_low_branch_align_only_prefix_tree_node_floor_mean = 1_437.531f64;
    let direct_restoring_final_low_branch_align_only_prefix_tree_node_floor_p99 = 1_568usize;
    let direct_restoring_final_low_branch_align_only_best_scan_gap = 31_872.834f64;
    let direct_restoring_final_low_branch_align_only_best_binary_gap = -73_294.652f64;
    let direct_restoring_final_low_branch_align_only_best_huffman_gap = -90_934.535f64;
    let direct_restoring_final_low_branch_align_only_best_prefix_tree_gap = -106_130.130f64;
    let direct_restoring_final_low_branch_align_only_mixed4to8_scan_gap = 33_880.682f64;
    let direct_restoring_final_low_branch_align_only_mixed4to8_prefix_tree_gap = -104_122.283f64;
    let direct_restoring_final_low_branch_align_only_best_lookup_target_mean = 14_703.797f64;
    let direct_restoring_final_low_branch_align_only_best_lookup_multiplier_budget =
        2.653_174f64;
    let direct_restoring_final_low_branch_align_only_scan_over_binary_multiplier = 3.372_071f64;
    let direct_restoring_final_low_branch_align_only_huffman_over_binary_multiplier =
        0.602_129f64;
    let direct_restoring_final_low_branch_align_only_prefix_tree_over_binary_multiplier =
        0.259_390f64;
    let direct_restoring_final_low_branch_align_only_support_noncontig_steps = 61usize;
    let direct_restoring_final_low_branch_align_only_support_max_span = 24usize;
    let direct_restoring_final_prefix_bit_reader_toy_eq_ccx = 4usize;
    let direct_restoring_final_prefix_bit_reader_toy_dynamic_read_ccx = 16usize;
    let direct_restoring_final_prefix_bit_reader_toy_reader_forward_ccx = 20usize;
    let direct_restoring_final_prefix_bit_reader_toy_tree_forward_ccx = 6usize;
    let direct_restoring_final_prefix_bit_reader_toy_full_forward_ccx = 26usize;
    let direct_restoring_final_prefix_bit_reader_toy_roundtrip_ccx = 52usize;
    let direct_restoring_final_prefix_bit_reader_toy_peak_q = 31usize;
    let direct_restoring_final_prefix_bit_reader_toy_cursor_states = 4usize;
    let direct_restoring_final_prefix_bit_reader_toy_internal_nodes = 4usize;
    let direct_restoring_final_prefix_bit_reader_toy_reader_over_tree = 3.333_333f64;
    let direct_restoring_final_prefix_bit_reader_toy_tree_over_node_roundtrip = 1.500_000f64;
    let direct_restoring_final_prefix_bit_reader_toy_full_over_node_roundtrip = 6.500_000f64;
    let direct_restoring_final_prefix_bit_reader_toy_roundtrip_ratio_budget = 10.228_508f64;
    let direct_restoring_final_prefix_bit_reader_toy_tree_only_scaled_gap = -100_380.006f64;
    let direct_restoring_final_prefix_bit_reader_toy_cursor_scaled_gap = -42_878.766f64;
    let direct_restoring_final_prefix_bit_reader_toy_dirty_restore_cases = 0usize;
    let direct_restoring_final_prefix_bit_reader_toy_dirty_history_cases = 0usize;
    let direct_restoring_final_prefix_bit_reader_toy_dirty_phase_cases = 0usize;
    let direct_restoring_final_prefix_cursor_advance_toy_ccx = 6usize;
    let direct_restoring_final_prefix_cursor_advance_toy_peak_q = 11usize;
    let direct_restoring_final_prefix_cursor_advance_toy_combined_roundtrip_ccx = 58usize;
    let direct_restoring_final_prefix_cursor_advance_toy_combined_over_node_roundtrip =
        7.250_000f64;
    let direct_restoring_final_prefix_cursor_advance_toy_roundtrip_ratio_budget = 10.228_508f64;
    let direct_restoring_final_prefix_cursor_advance_toy_combined_scaled_gap = -34_253.580f64;
    let direct_restoring_final_prefix_cursor_advance_toy_dirty_restore_cases = 0usize;
    let direct_restoring_final_prefix_cursor_advance_toy_dirty_history_cases = 0usize;
    let direct_restoring_final_prefix_cursor_advance_toy_dirty_phase_cases = 0usize;
    let direct_restoring_final_prefix_block2_toy_tree_ccx = 12usize;
    let direct_restoring_final_prefix_block2_toy_read2_ccx = 16usize;
    let direct_restoring_final_prefix_block2_toy_cursor_add_ccx = 10usize;
    let direct_restoring_final_prefix_block2_toy_decode_forward_ccx = 28usize;
    let direct_restoring_final_prefix_block2_toy_total_ccx = 66usize;
    let direct_restoring_final_prefix_block2_toy_peak_q = 52usize;
    let direct_restoring_final_prefix_block2_toy_over_node_roundtrip = 4.125_000f64;
    let direct_restoring_final_prefix_block2_toy_roundtrip_ratio_budget = 10.228_508f64;
    let direct_restoring_final_prefix_block2_toy_scaled_gap = -70_191.855f64;
    let direct_restoring_final_prefix_block2_toy_dirty_restore_cases = 0usize;
    let direct_restoring_final_prefix_block2_toy_dirty_history_cases = 0usize;
    let direct_restoring_final_prefix_block2_toy_dirty_phase_cases = 0usize;
    let direct_restoring_final_prefix_block2_consume_toy_decode_forward_ccx = 28usize;
    let direct_restoring_final_prefix_block2_consume_toy_cursor_add_ccx = 10usize;
    let direct_restoring_final_prefix_block2_consume_toy_consume_ccx = 4usize;
    let direct_restoring_final_prefix_block2_consume_toy_parser_transient_ccx = 76usize;
    let direct_restoring_final_prefix_block2_consume_toy_total_ccx = 80usize;
    let direct_restoring_final_prefix_block2_consume_toy_peak_q = 56usize;
    let direct_restoring_final_prefix_block2_consume_toy_parser_over_node_roundtrip =
        4.750_000f64;
    let direct_restoring_final_prefix_block2_consume_toy_roundtrip_ratio_budget = 10.228_508f64;
    let direct_restoring_final_prefix_block2_consume_toy_parser_scaled_gap = -63_004.200f64;
    let direct_restoring_final_prefix_block2_consume_toy_dirty_restore_cases = 0usize;
    let direct_restoring_final_prefix_block2_consume_toy_dirty_history_cases = 0usize;
    let direct_restoring_final_prefix_block2_consume_toy_dirty_phase_cases = 0usize;
    let direct_restoring_final_prefix_block2_leaf_touch_toy_decode_forward_ccx = 28usize;
    let direct_restoring_final_prefix_block2_leaf_touch_toy_leaf_touch_ccx = 40usize;
    let direct_restoring_final_prefix_block2_leaf_touch_toy_parser_transient_ccx = 56usize;
    let direct_restoring_final_prefix_block2_leaf_touch_toy_total_ccx = 96usize;
    let direct_restoring_final_prefix_block2_leaf_touch_toy_peak_q = 43usize;
    let direct_restoring_final_prefix_block2_leaf_touch_toy_parser_over_node_roundtrip =
        3.500_000f64;
    let direct_restoring_final_prefix_block2_leaf_touch_toy_total_over_node_roundtrip =
        6.000_000f64;
    let direct_restoring_final_prefix_block2_leaf_touch_toy_roundtrip_ratio_budget =
        10.228_508f64;
    let direct_restoring_final_prefix_block2_leaf_touch_toy_parser_scaled_gap = -77_379.510f64;
    let direct_restoring_final_prefix_block2_leaf_touch_toy_total_scaled_gap = -48_628.890f64;
    let direct_restoring_final_prefix_block2_leaf_touch_toy_dirty_restore_cases = 0usize;
    let direct_restoring_final_prefix_block2_leaf_touch_toy_dirty_history_cases = 0usize;
    let direct_restoring_final_prefix_block2_leaf_touch_toy_dirty_phase_cases = 0usize;
    let direct_restoring_final_prefix_block2_selected_addsub_toy_decode_forward_ccx = 28usize;
    let direct_restoring_final_prefix_block2_selected_addsub_toy_select_shift_ccx = 60usize;
    let direct_restoring_final_prefix_block2_selected_addsub_toy_addsub_ccx = 12usize;
    let direct_restoring_final_prefix_block2_selected_addsub_toy_arithmetic_ccx = 72usize;
    let direct_restoring_final_prefix_block2_selected_addsub_toy_parser_transient_ccx = 56usize;
    let direct_restoring_final_prefix_block2_selected_addsub_toy_total_ccx = 128usize;
    let direct_restoring_final_prefix_block2_selected_addsub_toy_peak_q = 56usize;
    let direct_restoring_final_prefix_block2_selected_addsub_toy_parser_over_node_roundtrip =
        3.500_000f64;
    let direct_restoring_final_prefix_block2_selected_addsub_toy_arithmetic_over_node_roundtrip =
        4.500_000f64;
    let direct_restoring_final_prefix_block2_selected_addsub_toy_total_over_node_roundtrip =
        8.000_000f64;
    let direct_restoring_final_prefix_block2_selected_addsub_toy_roundtrip_ratio_budget =
        10.228_508f64;
    let direct_restoring_final_prefix_block2_selected_addsub_toy_parser_scaled_gap =
        -77_379.510f64;
    let direct_restoring_final_prefix_block2_selected_addsub_toy_total_scaled_gap =
        -25_628.394f64;
    let direct_restoring_final_prefix_block2_selected_addsub_toy_dirty_restore_cases = 0usize;
    let direct_restoring_final_prefix_block2_selected_addsub_toy_dirty_history_cases = 0usize;
    let direct_restoring_final_prefix_block2_selected_addsub_toy_dirty_phase_cases = 0usize;
    let direct_restoring_final_prefix_block2_balanced_selected_addsub_toy_tree_ccx = 12usize;
    let direct_restoring_final_prefix_block2_balanced_selected_addsub_toy_read2_ccx = 6usize;
    let direct_restoring_final_prefix_block2_balanced_selected_addsub_toy_decode_forward_ccx =
        18usize;
    let direct_restoring_final_prefix_block2_balanced_selected_addsub_toy_select_shift_ccx =
        60usize;
    let direct_restoring_final_prefix_block2_balanced_selected_addsub_toy_addsub_ccx = 12usize;
    let direct_restoring_final_prefix_block2_balanced_selected_addsub_toy_arithmetic_ccx =
        72usize;
    let direct_restoring_final_prefix_block2_balanced_selected_addsub_toy_parser_transient_ccx =
        36usize;
    let direct_restoring_final_prefix_block2_balanced_selected_addsub_toy_total_ccx = 108usize;
    let direct_restoring_final_prefix_block2_balanced_selected_addsub_toy_peak_q = 51usize;
    let direct_restoring_final_prefix_block2_balanced_selected_addsub_toy_parser_over_node_roundtrip =
        2.250_000f64;
    let direct_restoring_final_prefix_block2_balanced_selected_addsub_toy_arithmetic_over_node_roundtrip =
        4.500_000f64;
    let direct_restoring_final_prefix_block2_balanced_selected_addsub_toy_total_over_node_roundtrip =
        6.750_000f64;
    let direct_restoring_final_prefix_block2_balanced_selected_addsub_toy_roundtrip_ratio_budget =
        10.228_508f64;
    let direct_restoring_final_prefix_block2_balanced_selected_addsub_toy_total_scaled_gap =
        -40_003.704f64;
    let direct_restoring_final_prefix_block2_balanced_selected_addsub_toy_dirty_restore_cases =
        0usize;
    let direct_restoring_final_prefix_block2_balanced_selected_addsub_toy_dirty_history_cases =
        0usize;
    let direct_restoring_final_prefix_block2_balanced_selected_addsub_toy_dirty_phase_cases =
        0usize;
    let direct_restoring_final_prefix_block2_selected_addsub_roundtrip_toy_forward_decode_ccx =
        28usize;
    let direct_restoring_final_prefix_block2_selected_addsub_roundtrip_toy_reverse_decode_ccx =
        28usize;
    let direct_restoring_final_prefix_block2_selected_addsub_roundtrip_toy_forward_select_shift_ccx =
        60usize;
    let direct_restoring_final_prefix_block2_selected_addsub_roundtrip_toy_reverse_select_shift_ccx =
        60usize;
    let direct_restoring_final_prefix_block2_selected_addsub_roundtrip_toy_forward_addsub_ccx =
        12usize;
    let direct_restoring_final_prefix_block2_selected_addsub_roundtrip_toy_reverse_addsub_ccx =
        12usize;
    let direct_restoring_final_prefix_block2_selected_addsub_roundtrip_toy_forward_ccx = 128usize;
    let direct_restoring_final_prefix_block2_selected_addsub_roundtrip_toy_reverse_ccx = 128usize;
    let direct_restoring_final_prefix_block2_selected_addsub_roundtrip_toy_arithmetic_ccx =
        144usize;
    let direct_restoring_final_prefix_block2_selected_addsub_roundtrip_toy_parser_transient_ccx =
        112usize;
    let direct_restoring_final_prefix_block2_selected_addsub_roundtrip_toy_total_ccx = 256usize;
    let direct_restoring_final_prefix_block2_selected_addsub_roundtrip_toy_peak_q = 56usize;
    let direct_restoring_final_prefix_block2_selected_addsub_roundtrip_toy_parser_over_node_roundtrip =
        3.500_000f64;
    let direct_restoring_final_prefix_block2_selected_addsub_roundtrip_toy_arithmetic_over_node_roundtrip =
        4.500_000f64;
    let direct_restoring_final_prefix_block2_selected_addsub_roundtrip_toy_total_over_node_roundtrip =
        8.000_000f64;
    let direct_restoring_final_prefix_block2_selected_addsub_roundtrip_toy_roundtrip_ratio_budget =
        10.228_508f64;
    let direct_restoring_final_prefix_block2_selected_addsub_roundtrip_toy_parser_scaled_gap =
        -77_379.510f64;
    let direct_restoring_final_prefix_block2_selected_addsub_roundtrip_toy_total_scaled_gap =
        -25_628.394f64;
    let direct_restoring_final_prefix_block2_selected_addsub_roundtrip_toy_dirty_restore_cases =
        0usize;
    let direct_restoring_final_prefix_block2_selected_addsub_roundtrip_toy_dirty_history_cases =
        0usize;
    let direct_restoring_final_prefix_block2_selected_addsub_roundtrip_toy_dirty_phase_cases =
        0usize;
    let direct_restoring_final_prefix_block2_span24_roundtrip_toy_forward_decode_ccx = 28usize;
    let direct_restoring_final_prefix_block2_span24_roundtrip_toy_reverse_decode_ccx = 28usize;
    let direct_restoring_final_prefix_block2_span24_roundtrip_toy_forward_select_shift_ccx =
        60usize;
    let direct_restoring_final_prefix_block2_span24_roundtrip_toy_reverse_select_shift_ccx =
        60usize;
    let direct_restoring_final_prefix_block2_span24_roundtrip_toy_forward_addsub_ccx = 50usize;
    let direct_restoring_final_prefix_block2_span24_roundtrip_toy_reverse_addsub_ccx = 50usize;
    let direct_restoring_final_prefix_block2_span24_roundtrip_toy_total_ccx = 332usize;
    let direct_restoring_final_prefix_block2_span24_roundtrip_toy_peak_q = 113usize;
    let direct_restoring_final_prefix_block2_span24_roundtrip_toy_parser_over_node_roundtrip =
        3.500_000f64;
    let direct_restoring_final_prefix_block2_span24_roundtrip_toy_arithmetic_over_node_roundtrip =
        6.875_000f64;
    let direct_restoring_final_prefix_block2_span24_roundtrip_toy_total_over_node_roundtrip =
        10.375_000f64;
    let direct_restoring_final_prefix_block2_span24_roundtrip_toy_roundtrip_ratio_budget =
        10.228_508f64;
    let direct_restoring_final_prefix_block2_span24_roundtrip_toy_total_scaled_gap =
        1_684.695f64;
    let direct_restoring_final_prefix_block2_span24_taper_materialized_full_add_per_digit =
        55usize;
    let direct_restoring_final_prefix_block2_span24_taper_add_per_digit_floor = 63usize;
    let direct_restoring_final_prefix_block2_span24_taper_arithmetic_floor = 252usize;
    let direct_restoring_final_prefix_block2_span24_taper_total_floor = 364usize;
    let direct_restoring_final_prefix_block2_span24_taper_total_over_node_roundtrip =
        11.375_000f64;
    let direct_restoring_final_prefix_block2_span24_taper_scaled_gap = 13_184.943f64;
    let direct_restoring_final_prefix_block2_span24_roundtrip_toy_dirty_restore_cases = 0usize;
    let direct_restoring_final_prefix_block2_span24_roundtrip_toy_dirty_history_cases = 0usize;
    let direct_restoring_final_prefix_block2_span24_roundtrip_toy_dirty_phase_cases = 0usize;
    let direct_restoring_final_low_branch_prefix_support_weighted_prefix_node_mean =
        1_437.531f64;
    let direct_restoring_final_low_branch_prefix_support_weighted_prefix_node_p99 = 1_568usize;
    let direct_restoring_final_low_branch_prefix_support_weighted_materialized_digit_mean =
        10_990.740f64;
    let direct_restoring_final_low_branch_prefix_support_weighted_materialized_digit_p99 =
        12_023usize;
    let direct_restoring_final_low_branch_prefix_support_weighted_tree_decode_mean =
        2_665.870f64;
    let direct_restoring_final_low_branch_prefix_support_weighted_dynamic_even_mean =
        7_330.417f64;
    let direct_restoring_final_low_branch_prefix_support_weighted_variable_decode_mean =
        9_996.286f64;
    let direct_restoring_final_low_branch_prefix_support_weighted_variable_offset1_decode_mean =
        9_917.655f64;
    let direct_restoring_final_low_branch_prefix_support_weighted_arithmetic_over_node_roundtrip =
        3.822_784f64;
    let direct_restoring_final_low_branch_prefix_support_weighted_total_over_node_roundtrip =
        7.322_784f64;
    let direct_restoring_final_low_branch_prefix_support_weighted_variable_total_over_node_roundtrip =
        10.776_573f64;
    let direct_restoring_final_low_branch_prefix_support_weighted_variable_offset1_total_over_node_roundtrip =
        10.721_874f64;
    let direct_restoring_final_low_branch_prefix_support_weighted_balanced_total_over_node_roundtrip =
        5.916_745f64;
    let direct_restoring_final_low_branch_prefix_support_weighted_selective_total_over_node_roundtrip =
        6.491_114f64;
    let direct_restoring_final_low_branch_prefix_support_weighted_ratio_budget = 10.228_508f64;
    let direct_restoring_final_low_branch_prefix_support_weighted_gap = -33_416.553f64;
    let direct_restoring_final_low_branch_prefix_support_weighted_variable_gap = 6_302.872f64;
    let direct_restoring_final_low_branch_prefix_support_weighted_variable_offset1_gap =
        5_673.821f64;
    let direct_restoring_final_low_branch_prefix_support_weighted_balanced_gap = -49_586.352f64;
    let direct_restoring_final_low_branch_prefix_support_weighted_selective_gap = -42_980.963f64;
    let direct_restoring_final_low_branch_prefix_support_weighted_projected_toffoli =
        2_666_583.447f64;
    let direct_restoring_final_low_branch_prefix_support_weighted_variable_projected_toffoli =
        2_706_302.872f64;
    let direct_restoring_final_low_branch_prefix_support_weighted_variable_offset1_projected_toffoli =
        2_705_673.821f64;
    let direct_restoring_final_low_branch_prefix_support_weighted_balanced_projected_toffoli =
        2_650_413.648f64;
    let direct_restoring_final_low_branch_prefix_support_weighted_selective_projected_toffoli =
        2_657_019.037f64;
    let direct_restoring_final_low_branch_prefix_support_weighted_shannon_prefix_bit_p99 =
        322usize;
    let direct_restoring_final_low_branch_prefix_support_weighted_balanced_prefix_bit_p99 =
        415usize;
    let direct_restoring_final_low_branch_prefix_support_weighted_selective_prefix_bit_p99 =
        381usize;
    let direct_restoring_final_low_branch_prefix_support_weighted_selective_prefix_bit_max =
        394usize;
    let direct_restoring_final_low_branch_prefix_support_weighted_selective_prefix_scratch_p99 =
        663usize;
    let direct_restoring_final_low_branch_prefix_support_weighted_selective_flatten_steps = 92usize;
    let direct_restoring_final_low_branch_prefix_support_weighted_maxconstrained_prefix_bit_mean =
        338.549f64;
    let direct_restoring_final_low_branch_prefix_support_weighted_maxconstrained_prefix_bit_p99 =
        372usize;
    let direct_restoring_final_low_branch_prefix_support_weighted_maxconstrained_prefix_bit_max =
        381usize;
    let direct_restoring_final_low_branch_prefix_support_weighted_maxconstrained_prefix_scratch_max =
        663usize;
    let direct_restoring_final_low_branch_prefix_support_weighted_maxconstrained_flatten_steps =
        83usize;
    let direct_restoring_final_low_branch_prefix_support_weighted_maxconstrained_trimmed_steps =
        9usize;
    let direct_restoring_final_low_branch_prefix_support_weighted_maxconstrained_over_budget_rows =
        0usize;
    let direct_restoring_final_low_branch_prefix_support_weighted_maxconstrained_over_budget_mass =
        0isize;
    let direct_restoring_final_low_branch_prefix_support_weighted_selective_schedule_codebook_steps =
        126usize;
    let direct_restoring_final_low_branch_prefix_support_weighted_selective_schedule_max_code_len =
        13usize;
    let direct_restoring_final_low_branch_prefix_support_weighted_selective_schedule_max_len_classes =
        12usize;
    let direct_restoring_final_low_branch_prefix_support_weighted_selective_schedule_bit_mean =
        348.210f64;
    let direct_restoring_final_low_branch_prefix_support_weighted_selective_schedule_bit_p99 =
        381usize;
    let direct_restoring_final_low_branch_prefix_support_weighted_selective_schedule_max_bits =
        394usize;
    let direct_restoring_final_low_branch_prefix_support_weighted_selective_schedule_decoded_symbols =
        856_854usize;
    let direct_restoring_final_low_branch_prefix_support_weighted_selective_schedule_prefix_collisions =
        0usize;
    let direct_restoring_final_low_branch_prefix_support_weighted_selective_schedule_decode_mismatches =
        0usize;
    let direct_restoring_final_low_branch_prefix_support_weighted_selective_schedule_cursor_mismatches =
        0usize;
    let direct_restoring_final_low_branch_prefix_support_weighted_maxconstrained_schedule_codebook_steps =
        126usize;
    let direct_restoring_final_low_branch_prefix_support_weighted_maxconstrained_schedule_max_code_len =
        13usize;
    let direct_restoring_final_low_branch_prefix_support_weighted_maxconstrained_schedule_max_len_classes =
        12usize;
    let direct_restoring_final_low_branch_prefix_support_weighted_maxconstrained_schedule_bit_mean =
        338.549f64;
    let direct_restoring_final_low_branch_prefix_support_weighted_maxconstrained_schedule_bit_p99 =
        372usize;
    let direct_restoring_final_low_branch_prefix_support_weighted_maxconstrained_schedule_max_bits =
        381usize;
    let direct_restoring_final_low_branch_prefix_support_weighted_maxconstrained_schedule_decoded_symbols =
        856_854usize;
    let direct_restoring_final_low_branch_prefix_support_weighted_maxconstrained_schedule_prefix_collisions =
        0usize;
    let direct_restoring_final_low_branch_prefix_support_weighted_maxconstrained_schedule_decode_mismatches =
        0usize;
    let direct_restoring_final_low_branch_prefix_support_weighted_maxconstrained_schedule_cursor_mismatches =
        0usize;
    let direct_restoring_final_low_branch_prefix_support_weighted_selective_dynamic_even_mean =
        1_169.937f64;
    let direct_restoring_final_low_branch_prefix_support_weighted_selective_variable_decode_mean =
        3_835.807f64;
    let direct_restoring_final_low_branch_prefix_support_weighted_maxconstrained_dynamic_even_mean =
        1_734.331f64;
    let direct_restoring_final_low_branch_prefix_support_weighted_maxconstrained_dynamic_even_p99 =
        1_807usize;
    let direct_restoring_final_low_branch_prefix_support_weighted_maxconstrained_variable_decode_mean =
        4_400.200f64;
    let direct_restoring_final_low_branch_prefix_support_weighted_maxconstrained_variable_decode_p99 =
        4_707usize;
    let direct_restoring_final_low_branch_prefix_support_weighted_maxconstrained_total_over_node_roundtrip =
        6.883_727f64;
    let direct_restoring_final_low_branch_prefix_support_weighted_maxconstrained_gap =
        -38_465.817f64;
    let direct_restoring_final_low_branch_prefix_support_weighted_maxconstrained_projected_toffoli =
        2_661_534.183f64;
    let direct_restoring_final_peakfit_toy_cases_with_sample_gap = 4usize;
    let direct_restoring_final_peakfit_toy_largest_missing_symbols = 598usize;
    let direct_restoring_final_peakfit_toy_largest_sample_over_budget_traces = 3_050usize;
    let direct_restoring_final_peakfit_toy_largest_exact_over_budget_traces = 3_132usize;
    let direct_restoring_final_low_branch_support_toy_cases_with_missing = 4usize;
    let direct_restoring_final_low_branch_support_toy_largest_missing_symbols = 26usize;
    let direct_restoring_final_low_branch_support_toy_largest_missing_steps = 11usize;
    let direct_restoring_final_low_branch_support_toy_largest_span_gap = 4usize;
    let direct_restoring_final_low_branch_support_toy_largest_exact_span = 16usize;
    let direct_restoring_final_low_branch_interval_toy_guard4_cover_cases = 0usize;
    let direct_restoring_final_low_branch_interval_toy_guard4_largest_missing_symbols =
        25usize;
    let direct_restoring_final_low_branch_interval_toy_guard4_largest_over_budget_traces =
        4_090usize;
    let direct_restoring_final_low_branch_interval_toy_guard4_largest_max_bits = 37usize;
    let direct_restoring_final_low_branch_interval_toy_full_cover_cases = 4usize;
    let direct_restoring_final_low_branch_interval_toy_full_fit_cases = 0usize;
    let direct_restoring_final_low_branch_interval_toy_full_largest_over_budget_traces =
        4_268usize;
    let direct_restoring_final_low_branch_interval_toy_full_largest_max_bits = 37usize;
    let direct_restoring_final_peakfit_holdout_missing_symbols = 182usize;
    let direct_restoring_final_peakfit_holdout_missing_traces = 170usize;
    let direct_restoring_final_peakfit_holdout_over_budget_rows = 7usize;
    let direct_restoring_final_peakfit_holdout_max_seen_bits = 391usize;
    let direct_restoring_final_peakfit_scaled_probe_train_samples = 65_536usize;
    let direct_restoring_final_peakfit_scaled_probe_holdout_samples = 32_768usize;
    let direct_restoring_final_peakfit_scaled_probe_flatten_steps = 45usize;
    let direct_restoring_final_peakfit_scaled_probe_missing_symbols = 105usize;
    let direct_restoring_final_peakfit_scaled_probe_missing_traces = 104usize;
    let direct_restoring_final_peakfit_scaled_probe_over_budget_rows = 0usize;
    let direct_restoring_final_peakfit_scaled_probe_max_seen_bits = 380usize;
    let direct_restoring_final_peakfit_scaled_probe_gap = -631.0f64;
    let direct_restoring_final_low_branch_prefix_support_weighted_span24_uniform_gap =
        1_684.686f64;
    let direct_restoring_final_low_branch_prefix_support_weighted_span24_symbol_mean = 1.000f64;
    let direct_restoring_final_low_branch_prefix_support_weighted_span24_symbol_p99 = 1usize;
    let direct_restoring_final_low_branch_prefix_support_weighted_support_noncontig_steps = 61usize;
    let direct_restoring_final_low_branch_prefix_support_weighted_support_max_span = 24usize;
    let direct_restoring_final_low_branch_prefix_support_weighted_support_max_symbols = 18usize;
    let direct_restoring_final_prefix_block2_balanced_family_toy_checked_circuits = 289usize;
    let direct_restoring_final_prefix_block2_balanced_family_toy_simulated_circuits = 49usize;
    let direct_restoring_final_prefix_block2_balanced_family_toy_simulated_cases = 26_728usize;
    let direct_restoring_final_prefix_block2_balanced_family_toy_max_support = 18usize;
    let direct_restoring_final_prefix_block2_balanced_family_toy_max_tree_ccx = 64usize;
    let direct_restoring_final_prefix_block2_balanced_family_toy_max_read2_ccx = 10usize;
    let direct_restoring_final_prefix_block2_balanced_family_toy_max_decode_forward_ccx = 74usize;
    let direct_restoring_final_prefix_block2_balanced_family_toy_max_select_shift_ccx = 216usize;
    let direct_restoring_final_prefix_block2_balanced_family_toy_max_addsub_ccx = 38usize;
    let direct_restoring_final_prefix_block2_balanced_family_toy_max_total_ccx = 804usize;
    let direct_restoring_final_prefix_block2_balanced_family_toy_max_peak_q = 186usize;
    let direct_restoring_final_prefix_block2_balanced_family_toy_max_total_over_node_roundtrip =
        7.666_667f64;
    let direct_restoring_final_prefix_block2_balanced_family_toy_max_ratio_support0 = 3usize;
    let direct_restoring_final_prefix_block2_balanced_family_toy_max_ratio_support1 = 2usize;
    let direct_restoring_final_prefix_block2_balanced_family_toy_max_total_scaled_gap =
        -29_461.810f64;
    let direct_restoring_final_prefix_block2_balanced_family_toy_dirty_restore_cases = 0usize;
    let direct_restoring_final_prefix_block2_balanced_family_toy_dirty_history_cases = 0usize;
    let direct_restoring_final_prefix_block2_balanced_family_toy_dirty_phase_cases = 0usize;
    let direct_restoring_final_coeff_decoder_alignment_degree_n14 = 13usize;
    let direct_restoring_final_coeff_decoder_alignment_density_n14 = 8_278usize;
    let direct_restoring_final_coeff_decoder_alignment_max_n14 = 13usize;
    let direct_restoring_final_align_entropy_variable_scratch_p99 = 602usize;
    let direct_restoring_final_align_entropy_variable_scratch_max = 620usize;
    let direct_restoring_final_align_entropy_global_scratch_p99 = 627usize;
    let direct_restoring_final_align_entropy_global_scratch_max = 650usize;
    let direct_restoring_final_align_entropy_step_scratch_p99 = 622usize;
    let direct_restoring_final_align_entropy_step_scratch_max = 640usize;
    let direct_restoring_final_align_prefix_global_scratch_p99 = 756usize;
    let direct_restoring_final_align_prefix_global_scratch_max = 781usize;
    let direct_restoring_final_align_prefix_step_scratch_p99 = 752usize;
    let direct_restoring_final_align_prefix_step_scratch_max = 771usize;
    let direct_restoring_final_align_entropy_branch_count_p99 = 117usize;
    let direct_restoring_final_align_entropy_branch_count_max = 125usize;
    let direct_restoring_final_range_parser_model_precision_bits = 13usize;
    let direct_restoring_final_range_parser_state_bits_p99 = 366usize;
    let direct_restoring_final_range_parser_live_scratch_p99 = 648usize;
    let direct_restoring_final_range_parser_symbol_count_p99 = 235usize;
    let direct_restoring_final_range_parser_state_touch_floor_mean = 71_167usize;
    let direct_restoring_final_range_parser_state_touch_floor_p99 = 84_835usize;
    let direct_restoring_final_range_parser_oneway_budget = 13_682usize;
    let direct_restoring_final_range_parser_augmented_mean_gap = 229_938isize;
    let direct_restoring_final_block_parser_model_precision_bits = 13usize;
    let direct_restoring_final_block_parser_oneway_budget = 13_682.500f64;
    let direct_restoring_final_block_parser_best_block_symbols = 8usize;
    let direct_restoring_final_block_parser_best_touch_floor_mean = 2_815.547f64;
    let direct_restoring_final_block_parser_best_touch_floor_p99 = 3_028usize;
    let direct_restoring_final_block_parser_best_compressed_bits_p99 = 381usize;
    let direct_restoring_final_block_parser_best_live_scratch_p99 = 663usize;
    let direct_restoring_final_block_parser_best_symbol_count_p99 = 235usize;
    let direct_restoring_final_block_parser_best_augmented_gap = -43_467.810f64;
    let direct_restoring_final_block32_touch_floor_mean = 10_719.669f64;
    let direct_restoring_final_block32_touch_floor_p99 = 11_618usize;
    let direct_restoring_final_block32_compressed_bits_p99 = 369usize;
    let direct_restoring_final_block32_live_scratch_p99 = 651usize;
    let direct_restoring_final_block32_symbol_count_p99 = 235usize;
    let direct_restoring_final_block32_augmented_gap = -11_851.323f64;
    let direct_restoring_final_block4_touch_floor_mean = 1_458.384f64;
    let direct_restoring_final_block4_touch_floor_p99 = 1_568usize;
    let direct_restoring_final_block4_compressed_bits_p99 = 393usize;
    let direct_restoring_final_block4_live_scratch_p99 = 675usize;
    let direct_restoring_final_block4_symbol_count_p99 = 235usize;
    let direct_restoring_final_block4_augmented_gap = -48_896.464f64;
    let direct_restoring_final_block5_touch_floor_mean = 1_797.940f64;
    let direct_restoring_final_block5_touch_floor_p99 = 1_936usize;
    let direct_restoring_final_block5_compressed_bits_p99 = 388usize;
    let direct_restoring_final_block5_live_scratch_p99 = 670usize;
    let direct_restoring_final_block5_symbol_count_p99 = 235usize;
    let direct_restoring_final_block5_augmented_gap = -47_538.241f64;
    let direct_restoring_final_block6_touch_floor_mean = 2_131.752f64;
    let direct_restoring_final_block6_touch_floor_p99 = 2_293usize;
    let direct_restoring_final_block6_compressed_bits_p99 = 384usize;
    let direct_restoring_final_block6_live_scratch_p99 = 666usize;
    let direct_restoring_final_block6_symbol_count_p99 = 235usize;
    let direct_restoring_final_block6_augmented_gap = -46_202.990f64;
    let direct_restoring_final_block7_touch_floor_mean = 2_477.987f64;
    let direct_restoring_final_block7_touch_floor_p99 = 2_664usize;
    let direct_restoring_final_block7_compressed_bits_p99 = 382usize;
    let direct_restoring_final_block7_live_scratch_p99 = 664usize;
    let direct_restoring_final_block7_symbol_count_p99 = 235usize;
    let direct_restoring_final_block7_augmented_gap = -44_818.050f64;
    let direct_restoring_final_block_parser_best_qrom_row_floor = 115_056usize;
    let direct_restoring_final_block_parser_best_qrom_max_rows_in_block = 4_934usize;
    let direct_restoring_final_block_parser_best_qrom_block_count_p99 = 30usize;
    let direct_restoring_final_block_parser_best_qrom_gap = 405_494.000f64;
    let direct_restoring_final_block32_qrom_row_floor = 56_169usize;
    let direct_restoring_final_block32_qrom_max_rows_in_block = 8_192usize;
    let direct_restoring_final_block32_qrom_block_count_p99 = 8usize;
    let direct_restoring_final_block_parser_lookup_scan_floor_mean = 18_856.559f64;
    let direct_restoring_final_block_parser_lookup_scan_floor_p99 = 20_579usize;
    let direct_restoring_final_block_parser_cond_branch_lookup_scan_floor_mean = 18_855.902f64;
    let direct_restoring_final_block_parser_cond_branch_lookup_scan_floor_p99 = 20_579usize;
    let direct_restoring_final_block_parser_best_with_lookup_mean = 21_672.106f64;
    let direct_restoring_final_block_parser_best_with_lookup_gap = 31_958.425f64;
    let direct_restoring_final_block_parser_binary_lookup_floor_mean = 6_796.417f64;
    let direct_restoring_final_block_parser_binary_lookup_floor_p99 = 7_605usize;
    let direct_restoring_final_block_parser_huffman_lookup_floor_mean = 4_511.413f64;
    let direct_restoring_final_block_parser_huffman_lookup_floor_p99 = 4_888usize;
    let direct_restoring_final_block_parser_best_with_binary_lookup_mean = 9_611.964f64;
    let direct_restoring_final_block_parser_best_with_binary_lookup_gap = -16_282.144f64;
    let direct_restoring_final_block_parser_best_with_binary_lookup_2x_mean = 16_408.380f64;
    let direct_restoring_final_block_parser_best_with_binary_lookup_2x_gap = 10_903.522f64;
    let direct_restoring_final_block4_with_binary_lookup_2x_mean = 15_051.217f64;
    let direct_restoring_final_block4_with_binary_lookup_2x_gap = 5_474.868f64;
    let direct_restoring_final_block4_lookup_multiplier_budget = 1.798_612f64;
    let direct_restoring_final_block5_with_binary_lookup_2x_mean = 15_390.773f64;
    let direct_restoring_final_block5_with_binary_lookup_2x_gap = 6_833.091f64;
    let direct_restoring_final_block5_lookup_multiplier_budget = 1.748_651f64;
    let direct_restoring_final_block7_with_binary_lookup_2x_mean = 16_070.820f64;
    let direct_restoring_final_block7_with_binary_lookup_2x_gap = 9_553.282f64;
    let direct_restoring_final_block7_lookup_multiplier_budget = 1.648_591f64;
    let direct_restoring_final_block_parser_cond_branch_best_block_symbols = 7usize;
    let direct_restoring_final_block_parser_cond_branch_touch_floor_mean = 2_461.660f64;
    let direct_restoring_final_block_parser_cond_branch_touch_floor_p99 = 2_647usize;
    let direct_restoring_final_block_parser_cond_branch_compressed_bits_p99 = 380usize;
    let direct_restoring_final_block_parser_cond_branch_live_scratch_p99 = 662usize;
    let direct_restoring_final_block_parser_cond_branch_augmented_gap = -44_883.358f64;
    let direct_restoring_final_block_parser_cond_branch_binary_lookup_floor_mean = 6_795.760f64;
    let direct_restoring_final_block_parser_cond_branch_binary_lookup_floor_p99 = 7_605usize;
    let direct_restoring_final_block_parser_cond_branch_huffman_lookup_floor_mean = 4_510.756f64;
    let direct_restoring_final_block_parser_cond_branch_huffman_lookup_floor_p99 = 4_888usize;
    let direct_restoring_final_block_parser_cond_branch_best_with_binary_lookup_mean =
        9_257.420f64;
    let direct_restoring_final_block_parser_cond_branch_best_with_binary_lookup_gap =
        -17_700.320f64;
    let direct_restoring_final_block_parser_cond_branch_best_with_binary_lookup_2x_mean =
        16_053.179f64;
    let direct_restoring_final_block_parser_cond_branch_best_with_binary_lookup_2x_gap =
        9_482.718f64;
    let direct_restoring_final_cond_block4_touch_floor_mean = 1_460.321f64;
    let direct_restoring_final_cond_block4_touch_floor_p99 = 1_572usize;
    let direct_restoring_final_cond_block4_compressed_bits_p99 = 394usize;
    let direct_restoring_final_cond_block4_live_scratch_p99 = 676usize;
    let direct_restoring_final_cond_block4_symbol_count_p99 = 235usize;
    let direct_restoring_final_cond_block4_augmented_gap = -48_888.715f64;
    let direct_restoring_final_cond_block5_touch_floor_mean = 1_787.362f64;
    let direct_restoring_final_cond_block5_touch_floor_p99 = 1_924usize;
    let direct_restoring_final_cond_block5_compressed_bits_p99 = 386usize;
    let direct_restoring_final_cond_block5_live_scratch_p99 = 668usize;
    let direct_restoring_final_cond_block5_symbol_count_p99 = 235usize;
    let direct_restoring_final_cond_block5_augmented_gap = -47_580.554f64;
    let direct_restoring_final_cond_block6_touch_floor_mean = 2_124.620f64;
    let direct_restoring_final_cond_block6_touch_floor_p99 = 2_287usize;
    let direct_restoring_final_cond_block6_compressed_bits_p99 = 383usize;
    let direct_restoring_final_cond_block6_live_scratch_p99 = 665usize;
    let direct_restoring_final_cond_block6_symbol_count_p99 = 235usize;
    let direct_restoring_final_cond_block6_augmented_gap = -46_231.521f64;
    let direct_restoring_final_cond_block7_touch_floor_mean = 2_461.660f64;
    let direct_restoring_final_cond_block7_touch_floor_p99 = 2_647usize;
    let direct_restoring_final_cond_block7_compressed_bits_p99 = 380usize;
    let direct_restoring_final_cond_block7_live_scratch_p99 = 662usize;
    let direct_restoring_final_cond_block7_symbol_count_p99 = 235usize;
    let direct_restoring_final_cond_block7_augmented_gap = -44_883.358f64;
    let direct_restoring_final_cond_block4_with_binary_lookup_2x_mean = 15_051.840f64;
    let direct_restoring_final_cond_block4_with_binary_lookup_2x_gap = 5_477.361f64;
    let direct_restoring_final_cond_block4_lookup_multiplier_budget = 1.798_501f64;
    let direct_restoring_final_cond_block5_with_binary_lookup_2x_mean = 15_378.881f64;
    let direct_restoring_final_cond_block5_with_binary_lookup_2x_gap = 6_785.522f64;
    let direct_restoring_final_cond_block5_lookup_multiplier_budget = 1.750_377f64;
    let direct_restoring_final_cond_block6_with_binary_lookup_2x_mean = 15_716.139f64;
    let direct_restoring_final_cond_block6_with_binary_lookup_2x_gap = 8_134.555f64;
    let direct_restoring_final_cond_block6_lookup_multiplier_budget = 1.700_749f64;
    let direct_restoring_final_cond_block7_with_binary_lookup_2x_mean = 16_053.179f64;
    let direct_restoring_final_cond_block7_with_binary_lookup_2x_gap = 9_482.718f64;
    let direct_restoring_final_cond_block7_lookup_multiplier_budget = 1.651_153f64;
    let direct_restoring_final_cond_mixed67_best_period = 5usize;
    let direct_restoring_final_cond_mixed67_best_mask = 9usize;
    let direct_restoring_final_cond_mixed67_best_seven_count = 2usize;
    let direct_restoring_final_cond_mixed67_touch_floor_mean = 2_272.649f64;
    let direct_restoring_final_cond_mixed67_touch_floor_p99 = 2_450usize;
    let direct_restoring_final_cond_mixed67_compressed_bits_p99 = 381usize;
    let direct_restoring_final_cond_mixed67_live_scratch_p99 = 663usize;
    let direct_restoring_final_cond_mixed67_symbol_count_p99 = 235usize;
    let direct_restoring_final_cond_mixed67_augmented_gap = -45_639.404f64;
    let direct_restoring_final_cond_mixed67_with_binary_lookup_2x_mean = 15_864.168f64;
    let direct_restoring_final_cond_mixed67_with_binary_lookup_2x_gap = 8_726.672f64;
    let direct_restoring_final_cond_mixed67_lookup_multiplier_budget = 1.678_966f64;
    let direct_restoring_final_cond_mixed67_with_cond_scan_lookup_2x_mean = 39_984.452f64;
    let direct_restoring_final_cond_mixed67_with_cond_scan_lookup_2x_gap = 105_207.810f64;
    let direct_restoring_final_cond_mixed67_with_huffman_lookup_2x_mean = 11_294.160f64;
    let direct_restoring_final_cond_mixed67_with_huffman_lookup_2x_gap = -9_553.359f64;
    let direct_restoring_final_cond_mixed67_huffman_lookup_multiplier_budget = 2.529_477f64;
    let direct_restoring_final_cond_mixed4to8_best_period = 4usize;
    let direct_restoring_final_cond_mixed4to8_schedule_code = 8_656usize;
    let direct_restoring_final_cond_mixed4to8_touch_floor_mean = 2_260.848f64;
    let direct_restoring_final_cond_mixed4to8_touch_floor_p99 = 2_439usize;
    let direct_restoring_final_cond_mixed4to8_compressed_bits_p99 = 381usize;
    let direct_restoring_final_cond_mixed4to8_live_scratch_p99 = 663usize;
    let direct_restoring_final_cond_mixed4to8_symbol_count_p99 = 235usize;
    let direct_restoring_final_cond_mixed4to8_augmented_gap = -45_686.608f64;
    let direct_restoring_final_cond_mixed4to8_with_binary_lookup_2x_mean =
        15_852.367f64;
    let direct_restoring_final_cond_mixed4to8_with_binary_lookup_2x_gap = 8_679.468f64;
    let direct_restoring_final_cond_mixed4to8_lookup_multiplier_budget = 1.680_703f64;
    let direct_restoring_final_cond_mixed4to8_block_joint_binary_lookup_mean =
        4_881.906f64;
    let direct_restoring_final_cond_mixed4to8_block_joint_binary_lookup_p99 =
        5_356usize;
    let direct_restoring_final_cond_mixed4to8_block_joint_support_row_floor =
        68_058usize;
    let direct_restoring_final_cond_mixed4to8_block_joint_max_patterns = 4_368usize;
    let direct_restoring_final_cond_mixed4to8_block_joint_block_count_p99 = 38usize;
    let direct_restoring_final_cond_mixed4to8_with_block_joint_binary_lookup_2x_mean =
        12_024.660f64;
    let direct_restoring_final_cond_mixed4to8_with_block_joint_binary_lookup_2x_gap =
        -6_631.358f64;
    let direct_restoring_final_cond_mixed4to8_block_joint_lookup_multiplier_budget =
        2.339_589f64;
    let direct_restoring_final_cond_mixed4to8_with_block_joint_scan_lookup_2x_mean =
        138_376.848f64;
    let direct_restoring_final_cond_mixed4to8_with_block_joint_scan_lookup_2x_gap =
        498_777.392f64;
    let direct_restoring_final_selective_pair_lookup_baseline_mean = 6_795.760f64;
    let direct_restoring_final_selective_pair_lookup_selected_saving_mean = 26.854f64;
    let direct_restoring_final_selective_pair_lookup_required_saving_mean = 1_084.933f64;
    let direct_restoring_final_selective_pair_lookup_mean = 6_768.906f64;
    let direct_restoring_final_selective_pair_lookup_target_mean = 5_710.826f64;
    let direct_restoring_final_selective_pair_lookup_gap = 8_464.638f64;
    let direct_restoring_final_selective_pair_lookup_selected_positions = 3usize;
    let direct_restoring_final_selective_pair_lookup_support_rows = 132usize;
    let direct_restoring_final_selective_pair_lookup_max_patterns = 90usize;
    let direct_restoring_final_block_joint_rank_degree_n14 = 14usize;
    let direct_restoring_final_block_joint_rank_density_n14 = 8_098usize;
    let direct_restoring_final_block_joint_rank_max_rank_n14 = 2_938usize;
    let direct_restoring_final_block_joint_rank_max_patterns_n14 = 2_939usize;
    let direct_restoring_final_block_joint_rank_support_rows_n14 = 3_474usize;
    let direct_restoring_final_block_joint_rank_max_blocks_n14 = 4usize;
    let direct_restoring_final_huffman_tree_toy_compare_ccx = 30usize;
    let direct_restoring_final_huffman_tree_toy_forward_ccx = 34usize;
    let direct_restoring_final_huffman_tree_toy_roundtrip_ccx = 68usize;
    let direct_restoring_final_huffman_tree_toy_peak_q = 21usize;
    let direct_restoring_final_huffman_tree_toy_weighted_path_depth = 1.888_889f64;
    let direct_restoring_final_huffman_tree_toy_full_tree_nodes = 5usize;
    let direct_restoring_final_huffman_tree_toy_path_compare_ccx_mean = 11.333f64;
    let direct_restoring_final_huffman_tree_toy_full_over_path_ratio = 2.647_059f64;
    let direct_restoring_final_huffman_tree_toy_dirty_restore_cases = 0usize;
    let direct_restoring_final_huffman_tree_toy_dirty_history_cases = 0usize;
    let direct_restoring_final_huffman_tree_toy_dirty_phase_cases = 0usize;
    let direct_restoring_final_block_parser_align_support_noncontig_steps = 52usize;
    let direct_restoring_final_block_parser_align_support_offset_steps = 127usize;
    let direct_restoring_final_block_parser_align_support_max_span = 20usize;
    let plusminus_raw_scratch = 564usize;
    let plusminus_unary_scratch_p99 = 640usize;
    let plusminus_parser_over_strict = plusminus_unary_scratch_p99 - STRICT_SCRATCH;
    let plusminus_scaled_slack_scratch_max = 517usize;
    let plusminus_scaled_solinas_projected_max = 2_230_850usize;
    let plusminus_scaled_solinas_gap_max = plusminus_scaled_solinas_projected_max as isize - GOOGLE_LOW_QUBIT_TOFFOLI as isize;
    let plusminus_solinas_scale_chunk_no_threshold_ccx = 3_390usize;
    let plusminus_solinas_scale_chunk_no_threshold_peak = 822usize;
    let plusminus_solinas_scale_chunk_exact_ccx = 7_564usize;
    let plusminus_solinas_scale_chunk_exact_peak = 822usize;
    let plusminus_solinas_scale_chunk_primitive_extra = 566usize;
    let plusminus_solinas_scale_chunk_naive_overlap_scratch = 1_078usize;
    let plusminus_solinas_scale_chunk_naive_over_google =
        plusminus_solinas_scale_chunk_naive_overlap_scratch as isize
            - GOOGLE_LOW_QUBIT_SCRATCH as isize;
    let plusminus_solinas_scale_chunk_one_lane_reuse_scratch = 822usize;
    let plusminus_solinas_scale_chunk_one_lane_reuse_over_google =
        plusminus_solinas_scale_chunk_one_lane_reuse_scratch as isize
            - GOOGLE_LOW_QUBIT_SCRATCH as isize;
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
    let halfgcd_second_col_fixed_depth64_scratch_p99 = 515usize;
    let halfgcd_second_col_fixed_depth64_prefix_extract_width_sum_p99 = 16_504usize;
    let halfgcd_second_col_fixed_depth64_prefix_max_digits_p99 = 14usize;
    let halfgcd_second_col_fixed_depth64_prefix_bounded_barrel_bits = 5usize;
    let halfgcd_second_col_fixed_depth64_decoder_width_sum_p99 = 4_592usize;
    let halfgcd_second_col_fixed_depth64_decoder_max_digits_p99 = 14usize;
    let halfgcd_second_col_fixed_depth64_decoder_bounded_barrel_bits = 5usize;
    let halfgcd_second_col_fixed_depth64_prefix_adversarial_prefix_max_digits = 256usize;
    let halfgcd_second_col_fixed_depth64_prefix_adversarial_decoder_max_digits = 256usize;
    let halfgcd_second_col_fixed_depth64_prefix_adversarial_required_barrel_bits = 8usize;
    let halfgcd_second_col_fixed_depth64_prefix_adversarial_missing_layers = 3usize;
    let halfgcd_second_col_fixed_depth64_prefix_adversarial_prefix_width_sum = 516usize;
    let halfgcd_second_col_fixed_depth64_prefix_adversarial_decoder_width_sum = 257usize;
    let halfgcd_second_col_fixed_depth64_prefix_full_domain_avg_gap_floor = 40_866isize;
    let halfgcd_second_col_fixed_depth64_tail_bits_p99 = 225usize;
    let halfgcd_second_col_fixed_depth64_tail_count_p99 = 108usize;
    let halfgcd_second_col_fixed_depth64_tail_width_sum_p99 = 9_136usize;
    let halfgcd_second_col_fixed_depth64_tail_max_q_bits_p99 = 14usize;
    let halfgcd_second_col_fixed_depth64_tail_bounded_barrel_bits = 5usize;
    let halfgcd_second_col_fixed_depth64_tail_adversarial_q_bits = 169usize;
    let halfgcd_second_col_fixed_depth64_tail_adversarial_required_barrel_bits = 8usize;
    let halfgcd_second_col_fixed_depth64_tail_adversarial_missing_layers = 3usize;
    let halfgcd_second_col_fixed_depth64_tail_adversarial_width_sum = 944usize;
    let halfgcd_second_col_fixed_depth64_tail_adversarial_count = 33usize;
    let halfgcd_second_col_fixed_depth64_tail_full_domain_avg_gap_floor = 40_860isize;
    let halfgcd_second_col_fixed_depth64_tail_extract_floor_p99 = 66_389usize;
    let halfgcd_second_col_fixed_depth64_tail_bounded_barrel_floor_p99 = 45_680usize;
    let halfgcd_second_col_fixed_depth64_tail_logbarrel_floor_p99 = 73_088usize;
    let halfgcd_second_col_fixed_depth64_exact_mean = 2_332_242usize;
    let halfgcd_second_col_fixed_depth64_exact_p99 = 2_474_014usize;
    let halfgcd_second_col_fixed_depth64_exact_tail_floor_mean = 2_533_612usize;
    let halfgcd_second_col_fixed_depth64_exact_tail_floor_p99 = 2_616_556usize;
    let halfgcd_second_col_fixed_depth64_exact_tail_floor_gap =
        halfgcd_second_col_fixed_depth64_exact_tail_floor_mean as isize
            - GOOGLE_LOW_QUBIT_TOFFOLI as isize;
    let halfgcd_second_col_fixed_depth64_exact_tail_bounded_barrel_mean = 2_663_148usize;
    let halfgcd_second_col_fixed_depth64_exact_tail_bounded_barrel_p99 = 2_729_614usize;
    let halfgcd_second_col_fixed_depth64_exact_tail_bounded_barrel_gap =
        halfgcd_second_col_fixed_depth64_exact_tail_bounded_barrel_mean as isize
            - GOOGLE_LOW_QUBIT_TOFFOLI as isize;
    let halfgcd_second_col_fixed_depth64_exact_tail_bounded_plus_one_width_mean =
        2_689_056usize;
    let halfgcd_second_col_fixed_depth64_exact_tail_bounded_plus_one_width_p99 =
        2_759_506usize;
    let halfgcd_second_col_fixed_depth64_exact_tail_bounded_plus_one_width_gap =
        halfgcd_second_col_fixed_depth64_exact_tail_bounded_plus_one_width_mean as isize
            - GOOGLE_LOW_QUBIT_TOFFOLI as isize;
    let halfgcd_second_col_fixed_depth64_exact_tail_bounded_plus_two_width_mean =
        2_714_963usize;
    let halfgcd_second_col_fixed_depth64_exact_tail_bounded_plus_two_width_p99 =
        2_791_246usize;
    let halfgcd_second_col_fixed_depth64_exact_tail_bounded_plus_two_width_gap =
        halfgcd_second_col_fixed_depth64_exact_tail_bounded_plus_two_width_mean as isize
            - GOOGLE_LOW_QUBIT_TOFFOLI as isize;
    let halfgcd_second_col_fixed_depth64_noscan_tail_bounded_barrel_mean = 2_535_006usize;
    let halfgcd_second_col_fixed_depth64_noscan_tail_bounded_barrel_p99 = 2_602_866usize;
    let halfgcd_second_col_fixed_depth64_exact_tail_logbarrel_mean = 2_740_870usize;
    let halfgcd_second_col_fixed_depth64_exact_tail_logbarrel_p99 = 2_822_826usize;
    let halfgcd_second_col_fixed_depth64_exact_tail_logbarrel_gap =
        halfgcd_second_col_fixed_depth64_exact_tail_logbarrel_mean as isize
            - GOOGLE_LOW_QUBIT_TOFFOLI as isize;
    let halfgcd_second_col_fixed_depth64_noscan_tail_logbarrel_mean = 2_612_728usize;
    let halfgcd_second_col_fixed_depth64_noscan_tail_logbarrel_p99 = 2_706_032usize;
    let halfgcd_second_col_fixed_depth64_exact_prefix_bounded_tail_logbarrel_mean =
        2_500_182usize;
    let halfgcd_second_col_fixed_depth64_exact_prefix_bounded_tail_logbarrel_p99 =
        2_585_888usize;
    let halfgcd_second_col_fixed_depth64_exact_prefix_bounded_tail_logbarrel_gap =
        halfgcd_second_col_fixed_depth64_exact_prefix_bounded_tail_logbarrel_mean as isize
            - GOOGLE_LOW_QUBIT_TOFFOLI as isize;
    let halfgcd_second_col_fixed_depth64_exact_prefix_bounded_tail_logbarrel_plus_one_width_mean =
        2_580_412usize;
    let halfgcd_second_col_fixed_depth64_exact_prefix_bounded_tail_logbarrel_plus_one_width_p99 =
        2_664_928usize;
    let halfgcd_second_col_fixed_depth64_exact_prefix_bounded_tail_logbarrel_plus_one_width_gap =
        halfgcd_second_col_fixed_depth64_exact_prefix_bounded_tail_logbarrel_plus_one_width_mean
            as isize
            - GOOGLE_LOW_QUBIT_TOFFOLI as isize;
    let halfgcd_second_col_fixed_depth64_exact_prefix_bounded_tail_logbarrel_plus_two_width_mean =
        2_660_641usize;
    let halfgcd_second_col_fixed_depth64_exact_prefix_bounded_tail_logbarrel_plus_two_width_p99 =
        2_743_340usize;
    let halfgcd_second_col_fixed_depth64_exact_prefix_bounded_tail_logbarrel_plus_two_width_gap =
        halfgcd_second_col_fixed_depth64_exact_prefix_bounded_tail_logbarrel_plus_two_width_mean
            as isize
            - GOOGLE_LOW_QUBIT_TOFFOLI as isize;
    let halfgcd_second_col_fixed_depth64_dynamic_barrel_static_mean = 2_740_052usize;
    let halfgcd_second_col_fixed_depth64_dynamic_barrel_static_p99 = 2_824_674usize;
    let halfgcd_second_col_fixed_depth64_dynamic_barrel_mean = 1_986_713usize;
    let halfgcd_second_col_fixed_depth64_dynamic_barrel_p99 = 2_047_416usize;
    let halfgcd_second_col_fixed_depth64_dynamic_barrel_gap =
        halfgcd_second_col_fixed_depth64_dynamic_barrel_mean as isize
            - GOOGLE_LOW_QUBIT_TOFFOLI as isize;
    let halfgcd_second_col_fixed_depth64_dynamic_barrel_savings_mean = 753_339usize;
    let halfgcd_second_col_fixed_depth64_dynamic_barrel_scratch_p99 = 515usize;
    let halfgcd_second_col_fixed_depth64_dynamic_prefix_decoder_static_mean = 160_541usize;
    let halfgcd_second_col_fixed_depth64_dynamic_prefix_decoder_mean = 19_250usize;
    let halfgcd_second_col_fixed_depth64_dynamic_tail_static_mean = 51_671usize;
    let halfgcd_second_col_fixed_depth64_dynamic_tail_mean = 4_628usize;
    let halfgcd_second_col_fixed_depth64_dynamic_high_layer_hits_p99 = 0usize;
    let halfgcd_second_col_fixed_depth64_static_app_mean = 2_934_322usize;
    let halfgcd_second_col_fixed_depth64_static_app_p99 = 3_010_096usize;
    let halfgcd_second_col_fixed_depth64_static_app_gap =
        halfgcd_second_col_fixed_depth64_static_app_mean as isize
            - GOOGLE_LOW_QUBIT_TOFFOLI as isize;
    let halfgcd_second_col_fixed_depth64_app_popcount_mean = 181_518usize;
    let halfgcd_second_col_fixed_depth64_app_static_floor_mean = 278_653usize;
    let halfgcd_second_col_fixed_depth64_app_static_over_popcount_mean = 97_135usize;
    let halfgcd_second_col_fixed_depth64_static_sep4_app_mean = 2_601_124usize;
    let halfgcd_second_col_fixed_depth64_static_sep4_app_p99 = 2_693_054usize;
    let halfgcd_second_col_fixed_depth64_static_sep4_app_gap =
        halfgcd_second_col_fixed_depth64_static_sep4_app_mean as isize
            - GOOGLE_LOW_QUBIT_TOFFOLI as isize;
    let halfgcd_second_col_fixed_depth64_static_joint4_app_mean = 2_545_006usize;
    let halfgcd_second_col_fixed_depth64_static_joint4_app_p99 = 2_644_284usize;
    let halfgcd_second_col_fixed_depth64_static_joint4_app_gap =
        halfgcd_second_col_fixed_depth64_static_joint4_app_mean as isize
            - GOOGLE_LOW_QUBIT_TOFFOLI as isize;
    let halfgcd_second_col_fixed_depth64_static_sep4_selector_budget_oneway = 49_438usize;
    let halfgcd_second_col_fixed_depth64_static_joint4_selector_budget_oneway = 77_497usize;
    let halfgcd_second_col_fixed_depth64_app_static_sep4_floor_mean = 112_054usize;
    let halfgcd_second_col_fixed_depth64_app_static_joint4_floor_mean = 83_995usize;
    let halfgcd_second_col_fixed_depth64_static_sep4_with_selector_floor_mean =
        2_824_277usize;
    let halfgcd_second_col_fixed_depth64_static_sep4_with_selector_floor_p99 =
        2_896_612usize;
    let halfgcd_second_col_fixed_depth64_static_joint4_with_selector_floor_mean =
        2_768_159usize;
    let halfgcd_second_col_fixed_depth64_static_joint4_with_selector_floor_p99 =
        2_845_292usize;
    let halfgcd_second_col_fixed_depth64_static_joint4_with_selector_floor_gap =
        halfgcd_second_col_fixed_depth64_static_joint4_with_selector_floor_mean as isize
            - GOOGLE_LOW_QUBIT_TOFFOLI as isize;
    let halfgcd_second_col_fixed_depth64_app_static_selector_floor_mean = 111_577usize;
    let halfgcd_second_col_fixed_depth64_app_static_selector_floor_over_joint4_budget =
        34_080usize;
    let halfgcd_second_col_fixed_depth64_static_window_scan_best_w = 6usize;
    let halfgcd_second_col_fixed_depth64_static_window_scan_best_mean = 2_749_506usize;
    let halfgcd_second_col_fixed_depth64_static_window_scan_best_p99 = 2_827_898usize;
    let halfgcd_second_col_fixed_depth64_static_window_scan_best_gap =
        halfgcd_second_col_fixed_depth64_static_window_scan_best_mean as isize
            - GOOGLE_LOW_QUBIT_TOFFOLI as isize;
    let halfgcd_second_col_fixed_depth64_static_window_scan_best_app_mean = 74_668usize;
    let halfgcd_second_col_fixed_depth64_static_window_scan_best_selector_mean =
        111_577usize;
    let halfgcd_second_col_fixed_depth64_static_window_scan_best_table_row_mean =
        76_657usize;
    let halfgcd_second_col_fixed_depth64_static_window_table_only_best_w = 4usize;
    let halfgcd_second_col_fixed_depth64_static_window_table_only_best_mean =
        2_559_198usize;
    let halfgcd_second_col_fixed_depth64_static_window_table_only_best_p99 =
        2_657_270usize;
    let halfgcd_second_col_fixed_depth64_static_window_table_only_best_gap =
        halfgcd_second_col_fixed_depth64_static_window_table_only_best_mean as isize
            - GOOGLE_LOW_QUBIT_TOFFOLI as isize;
    let halfgcd_second_col_fixed_depth64_static_window_table_source_best_w = 2usize;
    let halfgcd_second_col_fixed_depth64_static_window_table_source_mean =
        3_956_644usize;
    let halfgcd_second_col_fixed_depth64_static_window_table_source_p99 =
        4_257_224usize;
    let halfgcd_second_col_fixed_depth64_static_window_table_source_gap =
        halfgcd_second_col_fixed_depth64_static_window_table_source_mean as isize
            - GOOGLE_LOW_QUBIT_TOFFOLI as isize;
    let halfgcd_second_col_fixed_depth64_static_window_table_source_product_floor_mean =
        677_829usize;
    let halfgcd_second_col_fixed_depth64_static_window_required_selector_mean =
        86_824usize;
    let halfgcd_second_col_fixed_depth64_static_window_selector_cut_needed =
        24_753usize;
    let halfgcd_second_col_fixed_depth64_static_window_table_margin = 10_166usize;
    let halfgcd_second_col_fixed_depth64_static_window_source_product_best_w = 6usize;
    let halfgcd_second_col_fixed_depth64_static_window_source_product_mean =
        2_756_381usize;
    let halfgcd_second_col_fixed_depth64_static_window_source_product_p99 =
        2_834_890usize;
    let halfgcd_second_col_fixed_depth64_static_window_source_product_gap =
        halfgcd_second_col_fixed_depth64_static_window_source_product_mean as isize
            - GOOGLE_LOW_QUBIT_TOFFOLI as isize;
    let halfgcd_second_col_fixed_depth64_static_window_source_product_floor_mean =
        115_014usize;
    let halfgcd_second_col_fixed_depth64_static_window_source_product_table_row_mean =
        76_657usize;
    let halfgcd_second_col_fixed_depth64_static_window_wnaf_best_w = 6usize;
    let halfgcd_second_col_fixed_depth64_static_window_wnaf_mean = 2_748_271usize;
    let halfgcd_second_col_fixed_depth64_static_window_wnaf_p99 = 2_826_444usize;
    let halfgcd_second_col_fixed_depth64_static_window_wnaf_gap =
        halfgcd_second_col_fixed_depth64_static_window_wnaf_mean as isize
            - GOOGLE_LOW_QUBIT_TOFFOLI as isize;
    let halfgcd_second_col_fixed_depth64_static_window_wnaf_app_mean = 86_052usize;
    let halfgcd_second_col_fixed_depth64_static_window_wnaf_selector_floor_mean =
        99_575usize;
    let halfgcd_second_col_fixed_depth64_static_window_wnaf_source_product_floor_mean =
        99_575usize;
    let halfgcd_second_col_fixed_depth64_static_window_wnaf_table_row_floor_mean =
        32_463usize;
    let halfgcd_second_col_fixed_depth64_static_window_wnaf_positions_mean = 29.837f64;
    let halfgcd_second_col_fixed_depth64_static_window_wnaf_compact_best_w = 2usize;
    let halfgcd_second_col_fixed_depth64_static_window_wnaf_compact_mean =
        2_691_392usize;
    let halfgcd_second_col_fixed_depth64_static_window_wnaf_compact_p99 =
        2_775_864usize;
    let halfgcd_second_col_fixed_depth64_static_window_wnaf_compact_gap =
        halfgcd_second_col_fixed_depth64_static_window_wnaf_compact_mean as isize
            - GOOGLE_LOW_QUBIT_TOFFOLI as isize;
    let halfgcd_second_col_fixed_depth64_static_window_wnaf_compact_app_mean =
        119_092usize;
    let halfgcd_second_col_fixed_depth64_static_window_wnaf_compact_source_product_floor_mean =
        38_097usize;
    let halfgcd_second_col_fixed_depth64_static_window_wnaf_compact_missing_active_floor_mean =
        38_097usize;
    let halfgcd_second_col_fixed_depth64_static_window_wnaf_compact_active_slack_oneway =
        4_304usize;
    let halfgcd_second_col_fixed_depth64_static_window_wnaf_compact_table_row_floor_mean =
        497usize;
    let halfgcd_second_col_compact_wnaf_active_degree_n14 = 14usize;
    let halfgcd_second_col_compact_wnaf_active_density_n14 = 8_322usize;
    let halfgcd_second_col_compact_wnaf_active_positions_n14 = 15usize;
    let halfgcd_second_col_compact_wnaf_active_pair_positions_n14 = 15usize;
    let halfgcd_second_col_compact_wnaf_active_slots_n14 = 29usize;
    let halfgcd_second_col_compact_wnaf_active_full_slots_n14 = 30usize;
    let halfgcd_second_col_compact_wnaf_active_max_pair_n14 = 3usize;
    let halfgcd_second_col_alignment_mbu_degree_n14 = 14usize;
    let halfgcd_second_col_alignment_mbu_density_n14 = 8_142usize;
    let halfgcd_second_col_alignment_mbu_max_alignment_n14 = 13usize;
    let halfgcd_second_col_static_window_mbu_degree_n14 = 13usize;
    let halfgcd_second_col_static_window_mbu_density_n14 = 8_194usize;
    let halfgcd_second_col_static_window_mbu_max_coeff_bits_n14 = 14usize;
    let halfgcd_second_col_static_window_mbu_max_pair_n14 = 63usize;
    let halfgcd_second_col_static_window_wnaf_mbu_degree_n14 = 13usize;
    let halfgcd_second_col_static_window_wnaf_mbu_density_n14 = 8_162usize;
    let halfgcd_second_col_static_window_wnaf_mbu_max_positions_n14 = 15usize;
    let halfgcd_second_col_static_window_wnaf_mbu_max_pair_n14 = 63usize;
    let halfgcd_second_col_static_window_support_rows_n14 = 213usize;
    let halfgcd_second_col_static_window_support_full_rows_n14 = 315usize;
    let halfgcd_second_col_static_window_support_ppm_n14 = 676_190usize;
    let halfgcd_second_col_static_window_support_saturated_windows_n14 = 1usize;
    let halfgcd_second_col_static_window_support_windows_n14 = 5usize;
    let halfgcd_second_col_static_window_bit_support_n14 = 27usize;
    let halfgcd_second_col_static_window_full_bits_n14 = 28usize;
    let halfgcd_second_col_static_window_bit_support_ppm_n14 = 964_285usize;

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
    println!("METRIC scratch600_scaled_by_pattern_fixed_id_bits={scaled_by_pattern_fixed_id_bits}");
    println!("METRIC scratch600_scaled_by_pattern_fixed_id_distinct_rows={scaled_by_pattern_fixed_id_distinct_rows}");
    println!("METRIC scratch600_scaled_by_pattern_fixed_id_max_window_rows={scaled_by_pattern_fixed_id_max_window_rows}");
    println!("METRIC scratch600_scaled_by_pattern_fixed_id_nonzero_table_bits={scaled_by_pattern_fixed_id_nonzero_table_bits}");
    println!("METRIC scratch600_scaled_by_pattern_fixed_id_two_replay_before_decode={scaled_by_pattern_fixed_id_two_replay_before_decode}");
    println!("METRIC scratch600_scaled_by_pattern_fixed_id_remaining_to_2700k={scaled_by_pattern_fixed_id_remaining_to_2700k}");
    println!("METRIC scratch600_scaled_by_pattern_fixed_id_row_floor_gap={scaled_by_pattern_fixed_id_row_floor_gap}");
    println!("METRIC scratch600_scaled_by_pattern_fixed_id_bit_floor_gap={scaled_by_pattern_fixed_id_bit_floor_gap}");
    println!("METRIC scratch600_scaled_by_raw_pattern_bits={scaled_by_raw_pattern_bits}");
    println!("METRIC scratch600_scaled_by_raw_pattern_delta_bits={scaled_by_raw_pattern_delta_bits}");
    println!("METRIC scratch600_scaled_by_raw_pattern_single_a_scratch={scaled_by_raw_pattern_single_a_scratch}");
    println!("METRIC scratch600_scaled_by_raw_pattern_one_checkpoint_scratch={scaled_by_raw_pattern_one_checkpoint_scratch}");
    println!("METRIC scratch600_scaled_by_raw_pattern_window_a_scratch={scaled_by_raw_pattern_window_a_scratch}");
    println!("METRIC scratch600_scaled_by_raw_pattern_delta_checkpoint_max_rows={scaled_by_raw_pattern_delta_checkpoint_max_rows}");
    println!("METRIC scratch600_scaled_by_raw_pattern_delta_checkpoint_bits={scaled_by_raw_pattern_delta_checkpoint_bits}");
    println!("METRIC scratch600_scaled_by_raw_pattern_delta_checkpoint_scratch={scaled_by_raw_pattern_delta_checkpoint_scratch}");
    println!("METRIC scratch600_scaled_by_raw_pattern_delta_checkpoint_scratch_slack={scaled_by_raw_pattern_delta_checkpoint_scratch_slack}");
    println!("METRIC scratch600_scaled_by_raw_pattern_ambiguous_a_bits_mean_milli={scaled_by_raw_pattern_ambiguous_a_bits_mean_milli}");
    println!("METRIC scratch600_scaled_by_raw_pattern_ambiguous_a_bits_p99={scaled_by_raw_pattern_ambiguous_a_bits_p99}");
    println!("METRIC scratch600_scaled_by_raw_pattern_ambiguous_a_bits_max={scaled_by_raw_pattern_ambiguous_a_bits_max}");
    println!("METRIC scratch600_scaled_by_raw_pattern_two_replay_before_branch_decode={scaled_by_raw_pattern_two_replay_before_branch_decode}");
    println!("METRIC scratch600_scaled_by_raw_pattern_exact_decoder_per_replay={scaled_by_raw_pattern_exact_decoder_per_replay}");
    println!("METRIC scratch600_scaled_by_raw_pattern_exact_two_decoder_projection={scaled_by_raw_pattern_exact_two_decoder_projection}");
    println!("METRIC scratch600_scaled_by_raw_pattern_exact_two_decoder_gap={scaled_by_raw_pattern_exact_two_decoder_gap}");
    println!("METRIC scratch600_scaled_by_h_only_model_modular_windows={scaled_by_h_only_model_modular_windows}");
    println!("METRIC scratch600_scaled_by_h_only_model_modular_toffoli={scaled_by_h_only_model_modular_toffoli}");
    println!("METRIC scratch600_scaled_by_h_only_model_peak={scaled_by_h_only_model_peak}");
    println!("METRIC scratch600_scaled_by_h_only_model_history_bits={scaled_by_h_only_model_history_bits}");
    println!("METRIC scratch600_scaled_by_h_only_next_ratio_toy_n14_windows={scaled_by_h_only_next_ratio_toy_n14_windows}");
    println!("METRIC scratch600_scaled_by_h_only_next_ratio_toy_n14_keys={scaled_by_h_only_next_ratio_toy_n14_keys}");
    println!("METRIC scratch600_scaled_by_h_only_next_ratio_toy_n14_ambiguous_keys={scaled_by_h_only_next_ratio_toy_n14_ambiguous_keys}");
    println!("METRIC scratch600_scaled_by_h_only_next_ratio_toy_n14_max_next_h_choices={scaled_by_h_only_next_ratio_toy_n14_max_next_h_choices}");
    println!("METRIC scratch600_scaled_by_h_only_next_ratio_toy_n14_rank_p99={scaled_by_h_only_next_ratio_toy_n14_rank_p99}");
    println!("METRIC scratch600_scaled_by_h_only_next_ratio_toy_n14_rank_max={scaled_by_h_only_next_ratio_toy_n14_rank_max}");
    println!("METRIC scratch600_scaled_by_h_only_next_ratio_toy_n14_rank_mean_milli={scaled_by_h_only_next_ratio_toy_n14_rank_mean_milli}");
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
    println!("METRIC scratch600_by_centered_exactparity_clean_replay_ccx={by_centered_exactparity_clean_replay_ccx}");
    println!("METRIC scratch600_by_centered_exactparity_clean_peak_q={by_centered_exactparity_clean_peak_q}");
    println!("METRIC scratch600_by_centered_exactparity_clean_scratch_bits={by_centered_exactparity_clean_scratch_bits}");
    println!("METRIC scratch600_by_centered_exactparity_clean_per_div_budget={by_centered_exactparity_clean_per_div_budget}");
    println!("METRIC scratch600_by_centered_exactparity_two_clean_div_projection={by_centered_exactparity_two_clean_div_projection}");
    println!("METRIC scratch600_by_centered_exactparity_two_clean_div_gap_to_2700k={by_centered_exactparity_two_clean_div_gap}");
    println!("METRIC scratch600_centered_raw_scratch_bits={centered_raw_scratch}");
    println!("METRIC scratch600_centered_boundary_scratch_p99={centered_boundary_scratch_p99}");
    println!("METRIC scratch600_centered_parser_over_strict_bits={centered_parser_over_strict}");
    println!("METRIC scratch600_direct_signnorm_raw_digit_scratch_p99={direct_signnorm_raw_digit_scratch_p99}");
    println!("METRIC scratch600_direct_signnorm_det_coeffsign_scratch_p99={direct_signnorm_det_coeffsign_scratch_p99}");
    println!("METRIC scratch600_direct_signnorm_det_coeffsign_scratch_gap_google={direct_signnorm_det_coeffsign_scratch_gap_google}");
    println!("METRIC scratch600_direct_signnorm_rank_scratch_p99={direct_signnorm_rank_scratch_p99}");
    println!("METRIC scratch600_direct_signnorm_rank_over_google_bits={direct_signnorm_rank_over_google}");
    println!("METRIC scratch600_direct_signnorm_ambiguous_rank_scratch_p99={direct_signnorm_ambiguous_rank_scratch_p99}");
    println!("METRIC scratch600_direct_signnorm_ambiguous_rank_over_google_bits={direct_signnorm_ambiguous_rank_over_google}");
    println!("METRIC scratch600_direct_signnorm_exact_split_p99={direct_signnorm_exact_split_p99}");
    println!("METRIC scratch600_direct_signnorm_exact_split_gap_to_2700k={direct_signnorm_exact_split_gap}");
    println!("METRIC scratch600_direct_signnorm_logsign_once_p99={direct_signnorm_logsign_once_p99}");
    println!("METRIC scratch600_direct_signnorm_logsign_split_p99={direct_signnorm_logsign_split_p99}");
    println!("METRIC scratch600_direct_signnorm_logsign_direct_rem_toy_ccx={direct_signnorm_logsign_direct_rem_toy_ccx}");
    println!("METRIC scratch600_direct_signnorm_logsign_direct_rem_toy_peak_q={direct_signnorm_logsign_direct_rem_toy_peak_q}");
    println!("METRIC scratch600_direct_signnorm_logsign_direct_rem_toy_phase_dirty_cases={direct_signnorm_logsign_direct_rem_toy_phase_dirty_cases}");
    println!("METRIC scratch600_direct_signnorm_logsign_exact_cneg257={direct_signnorm_logsign_exact_cneg257}");
    println!("METRIC scratch600_direct_signnorm_logsign_exact_rem_p99={direct_signnorm_logsign_exact_rem_p99}");
    println!("METRIC scratch600_direct_signnorm_logsign_exact_once_p99={direct_signnorm_logsign_exact_once_p99}");
    println!("METRIC scratch600_direct_signnorm_logsign_exact_split_p99={direct_signnorm_logsign_exact_split_p99}");
    println!("METRIC scratch600_direct_signnorm_logsign_no_rem_cneg_projection_p99={direct_signnorm_logsign_no_rem_cneg_projection_p99}");
    println!("METRIC scratch600_direct_signnorm_logsign_once_gap_to_2700k={direct_signnorm_logsign_once_gap}");
    println!("METRIC scratch600_direct_signnorm_logsign_split_gap_to_2700k={direct_signnorm_logsign_split_gap}");
    println!("METRIC scratch600_direct_signnorm_logsign_exact_once_gap_to_2700k={direct_signnorm_logsign_exact_once_gap}");
    println!("METRIC scratch600_direct_signnorm_logsign_exact_split_gap_to_2700k={direct_signnorm_logsign_exact_split_gap}");
    println!("METRIC scratch600_direct_signnorm_logsign_no_rem_cneg_gap_to_2700k={direct_signnorm_logsign_no_rem_cneg_gap}");
    println!("METRIC scratch600_direct_signnorm_prefinal_signed_remainder_p99={direct_signnorm_prefinal_signed_remainder_p99}");
    println!("METRIC scratch600_direct_signnorm_prefinal_signed_remainder_gap_to_2700k={direct_signnorm_prefinal_signed_remainder_gap}");
    println!("METRIC scratch600_direct_signnorm_prefinal_signed_remainder_count_p99={direct_signnorm_prefinal_signed_remainder_count_p99}");
    println!("METRIC scratch600_direct_signnorm_prefinal_signed_remainder_digit_payload_p99={direct_signnorm_prefinal_signed_remainder_digit_payload_p99}");
    println!("METRIC scratch600_direct_signnorm_prefinal_signed_remainder_width_extra_max={direct_signnorm_prefinal_signed_remainder_width_extra_max}");
    println!("METRIC scratch600_direct_signnorm_mbu_degree_n14={direct_signnorm_mbu_degree_n14}");
    println!("METRIC scratch600_direct_signnorm_mbu_density_n14={direct_signnorm_mbu_density_n14}");
    println!("METRIC scratch600_direct_signnorm_mbu_max_count_n14={direct_signnorm_mbu_max_count_n14}");
    println!("METRIC scratch600_direct_signnorm_reverse_collisions_n14={direct_signnorm_reverse_collisions_n14}");
    println!("METRIC scratch600_direct_signnorm_reverse_states_n14={direct_signnorm_reverse_states_n14}");
    println!("METRIC scratch600_direct_signnorm_reverse_total_steps_n14={direct_signnorm_reverse_total_steps_n14}");
    println!("METRIC scratch600_direct_signnorm_coeff_reverse_collisions_n14={direct_signnorm_coeff_reverse_collisions_n14}");
    println!("METRIC scratch600_direct_signnorm_coeff_reverse_states_n14={direct_signnorm_coeff_reverse_states_n14}");
    println!("METRIC scratch600_direct_signnorm_coeff_reverse_total_steps_n14={direct_signnorm_coeff_reverse_total_steps_n14}");
    println!("METRIC scratch600_direct_signnorm_coeff_reverse_max_mult_n14={direct_signnorm_coeff_reverse_max_mult_n14}");
    println!("METRIC scratch600_direct_signnorm_coeff_reverse_zero_coeff_cases_n14={direct_signnorm_coeff_reverse_zero_coeff_cases_n14}");
    println!("METRIC scratch600_direct_signnorm_det_sign_reverse_collisions_n14={direct_signnorm_det_sign_reverse_collisions_n14}");
    println!("METRIC scratch600_direct_signnorm_det_sign_reverse_states_n14={direct_signnorm_det_sign_reverse_states_n14}");
    println!("METRIC scratch600_direct_signnorm_det_sign_reverse_max_mult_n14={direct_signnorm_det_sign_reverse_max_mult_n14}");
    println!("METRIC scratch600_direct_signnorm_det_coeffsign_reverse_collisions_n14={direct_signnorm_det_coeffsign_reverse_collisions_n14}");
    println!("METRIC scratch600_direct_signnorm_det_coeffsign_reverse_states_n14={direct_signnorm_det_coeffsign_reverse_states_n14}");
    println!("METRIC scratch600_direct_signnorm_det_coeffsign_reverse_total_steps_n14={direct_signnorm_det_coeffsign_reverse_total_steps_n14}");
    println!("METRIC scratch600_direct_signnorm_det_coeffsign_reverse_max_mult_n14={direct_signnorm_det_coeffsign_reverse_max_mult_n14}");
    println!("METRIC scratch600_direct_signnorm_det_coeffsign_bad_det_cases_n14={direct_signnorm_det_coeffsign_bad_det_cases_n14}");
    println!("METRIC scratch600_direct_signnorm_det_coeffsign_low2_mismatches_n14={direct_signnorm_det_coeffsign_low2_mismatches_n14}");
    println!("METRIC scratch600_direct_signnorm_det_coeffsign_formula_mismatches_n14={direct_signnorm_det_coeffsign_formula_mismatches_n14}");
    println!("METRIC scratch600_direct_signnorm_det_coeffsign_predicate_p1_ccx={direct_signnorm_det_coeffsign_predicate_p1_ccx}");
    println!("METRIC scratch600_direct_signnorm_det_coeffsign_predicate_p1_peak_q={direct_signnorm_det_coeffsign_predicate_p1_peak_q}");
    println!("METRIC scratch600_direct_signnorm_det_coeffsign_predicate_p1_valid_odd_det_cases={direct_signnorm_det_coeffsign_predicate_p1_valid_odd_det_cases}");
    println!("METRIC scratch600_direct_signnorm_det_coeffsign_predicate_p3_ccx={direct_signnorm_det_coeffsign_predicate_p3_ccx}");
    println!("METRIC scratch600_direct_signnorm_det_coeffsign_predicate_p3_peak_q={direct_signnorm_det_coeffsign_predicate_p3_peak_q}");
    println!("METRIC scratch600_direct_signnorm_det_coeffsign_predicate_p3_valid_odd_det_cases={direct_signnorm_det_coeffsign_predicate_p3_valid_odd_det_cases}");
    println!("METRIC scratch600_direct_signnorm_signed_domain_relative_negative_toy_ccx={direct_signnorm_signed_domain_relative_negative_toy_ccx}");
    println!("METRIC scratch600_direct_signnorm_signed_domain_relative_negative_257_ccx={direct_signnorm_signed_domain_relative_negative_257_ccx}");
    println!("METRIC scratch600_direct_signnorm_signed_domain_floor_toy_ccx={direct_signnorm_signed_domain_floor_toy_ccx}");
    println!("METRIC scratch600_direct_signnorm_signed_domain_floor_toy_peak_q={direct_signnorm_signed_domain_floor_toy_peak_q}");
    println!("METRIC scratch600_direct_signnorm_signed_domain_floor_toy_final_negative_cases={direct_signnorm_signed_domain_floor_toy_final_negative_cases}");
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
    println!("METRIC scratch600_direct_restoring_final_branch_digit_toy_branch_ccx={direct_restoring_final_branch_digit_toy_branch_ccx}");
    println!("METRIC scratch600_direct_restoring_final_branch_digit_toy_forward_ccx={direct_restoring_final_branch_digit_toy_forward_ccx}");
    println!("METRIC scratch600_direct_restoring_final_branch_digit_toy_roundtrip_ccx={direct_restoring_final_branch_digit_toy_roundtrip_ccx}");
    println!("METRIC scratch600_direct_restoring_final_branch_digit_toy_peak_q={direct_restoring_final_branch_digit_toy_peak_q}");
    println!("METRIC scratch600_direct_restoring_final_branch_digit_toy_branch_one_cases={direct_restoring_final_branch_digit_toy_branch_one_cases}");
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
    println!("METRIC scratch600_direct_restoring_final_low_branch_adjacent_transitions_n14={direct_restoring_final_low_branch_adjacent_transitions_n14}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_adjacent_ambiguous_n14={direct_restoring_final_low_branch_adjacent_ambiguous_n14}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_adjacent_high_n14={direct_restoring_final_low_branch_adjacent_high_n14}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_adjacent_violations_n14={direct_restoring_final_low_branch_adjacent_violations_n14}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_adjacent_max_delta_n14={direct_restoring_final_low_branch_adjacent_max_delta_n14}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_adjacent_max_alignment_n14={direct_restoring_final_low_branch_adjacent_max_alignment_n14}");
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
    println!("METRIC scratch600_direct_restoring_final_stored_align_length_rank_scratch_p99={direct_restoring_final_stored_align_length_rank_scratch_p99}");
    println!("METRIC scratch600_direct_restoring_final_stored_align_length_rank_scratch_max={direct_restoring_final_stored_align_length_rank_scratch_max}");
    println!("METRIC scratch600_direct_restoring_final_stored_align_public_len_mismatch_p99={direct_restoring_final_stored_align_public_len_mismatch_p99}");
    println!("METRIC scratch600_direct_restoring_final_stored_align_public_len_rank_p99={direct_restoring_final_stored_align_public_len_rank_p99}");
    println!("METRIC scratch600_direct_restoring_final_stored_align_public_len_position_only_p99={direct_restoring_final_stored_align_public_len_position_only_p99}");
    println!("METRIC scratch600_direct_restoring_final_stored_align_public_len_position_plus3_p99={direct_restoring_final_stored_align_public_len_position_plus3_p99}");
    println!("METRIC scratch600_direct_restoring_final_stored_align_q_len_position_only_p99={direct_restoring_final_stored_align_q_len_position_only_p99}");
    println!("METRIC scratch600_direct_restoring_final_stored_align_digit_len_position_only_p99={direct_restoring_final_stored_align_digit_len_position_only_p99}");
    println!("METRIC scratch600_direct_restoring_final_stored_align_joint_len_position_only_p99={direct_restoring_final_stored_align_joint_len_position_only_p99}");
    println!("METRIC scratch600_direct_restoring_final_stored_align_joint_len_position_plus3_p99={direct_restoring_final_stored_align_joint_len_position_plus3_p99}");
    println!("METRIC scratch600_direct_restoring_final_stored_align_pop_barrel_p99={direct_restoring_final_stored_align_pop_barrel_p99}");
    println!("METRIC scratch600_direct_restoring_final_stored_align_branch_select_p99={direct_restoring_final_stored_align_branch_select_p99}");
    println!("METRIC scratch600_direct_restoring_final_stored_align_branch_count_p99={direct_restoring_final_stored_align_branch_count_p99}");
    println!("METRIC scratch600_direct_restoring_final_branch_final_current_branch_select_mean={direct_restoring_final_branch_final_current_branch_select_mean:.3}");
    println!("METRIC scratch600_direct_restoring_final_branch_final_current_branch_select_p99={direct_restoring_final_branch_final_current_branch_select_p99}");
    println!("METRIC scratch600_direct_restoring_final_branch_final_width_mean={direct_restoring_final_branch_final_width_mean:.3}");
    println!("METRIC scratch600_direct_restoring_final_branch_final_width_p99={direct_restoring_final_branch_final_width_p99}");
    println!("METRIC scratch600_direct_restoring_final_branch_final_width_minus1_mean={direct_restoring_final_branch_final_width_minus1_mean:.3}");
    println!("METRIC scratch600_direct_restoring_final_branch_final_selected_width_saving_mean={direct_restoring_final_branch_final_selected_width_saving_mean:.3}");
    println!("METRIC scratch600_direct_restoring_final_branch_final_low_path_width_saving_mean={direct_restoring_final_branch_final_low_path_width_saving_mean:.3}");
    println!("METRIC scratch600_direct_restoring_final_branch_final_low_path_width_minus1_saving_mean={direct_restoring_final_branch_final_low_path_width_minus1_saving_mean:.3}");
    println!("METRIC scratch600_direct_restoring_final_branch_final_low_extra_touch_mean={direct_restoring_final_branch_final_low_extra_touch_mean:.3}");
    println!("METRIC scratch600_direct_restoring_final_branch_final_branch_count_mean={direct_restoring_final_branch_final_branch_count_mean:.3}");
    println!("METRIC scratch600_direct_restoring_final_branch_final_branch_count_p99={direct_restoring_final_branch_final_branch_count_p99}");
    println!("METRIC scratch600_direct_restoring_final_branch_final_high_branch_mean={direct_restoring_final_branch_final_high_branch_mean:.3}");
    println!("METRIC scratch600_direct_restoring_final_branch_final_high_branch_p99={direct_restoring_final_branch_final_high_branch_p99}");
    println!("METRIC scratch600_direct_restoring_final_branch_final_alignment_diff_p99={direct_restoring_final_branch_final_alignment_diff_p99}");
    println!("METRIC scratch600_direct_restoring_final_branch_final_digit_len_diff_p99={direct_restoring_final_branch_final_digit_len_diff_p99}");
    println!("METRIC scratch600_direct_restoring_final_branch_final_high_adjacent_violations={direct_restoring_final_branch_final_high_adjacent_violations}");
    println!("METRIC scratch600_direct_restoring_final_branch_final_current_mixed4to8_gap_to_2700k={direct_restoring_final_branch_final_current_mixed4to8_gap:.3}");
    println!("METRIC scratch600_direct_restoring_final_branch_final_current_scan_mixed4to8_gap_to_2700k={direct_restoring_final_branch_final_current_scan_mixed4to8_gap:.3}");
    println!("METRIC scratch600_direct_restoring_final_branch_final_selected_width_lookup_target_mean={direct_restoring_final_branch_final_selected_width_lookup_target_mean:.3}");
    println!("METRIC scratch600_direct_restoring_final_branch_final_low_path_width_lookup_target_mean={direct_restoring_final_branch_final_low_path_width_lookup_target_mean:.3}");
    println!("METRIC scratch600_direct_restoring_final_branch_final_low_path_width_minus1_lookup_target_mean={direct_restoring_final_branch_final_low_path_width_minus1_lookup_target_mean:.3}");
    println!("METRIC scratch600_direct_restoring_final_branch_final_selected_width_lookup_multiplier_budget={direct_restoring_final_branch_final_selected_width_lookup_multiplier_budget:.6}");
    println!("METRIC scratch600_direct_restoring_final_branch_final_low_path_width_lookup_multiplier_budget={direct_restoring_final_branch_final_low_path_width_lookup_multiplier_budget:.6}");
    println!("METRIC scratch600_direct_restoring_final_branch_final_low_path_width_minus1_lookup_multiplier_budget={direct_restoring_final_branch_final_low_path_width_minus1_lookup_multiplier_budget:.6}");
    println!("METRIC scratch600_direct_restoring_final_branch_final_selected_width_mixed4to8_gap_to_2700k={direct_restoring_final_branch_final_selected_width_mixed4to8_gap:.3}");
    println!("METRIC scratch600_direct_restoring_final_branch_final_low_path_width_mixed4to8_gap_to_2700k={direct_restoring_final_branch_final_low_path_width_mixed4to8_gap:.3}");
    println!("METRIC scratch600_direct_restoring_final_branch_final_low_path_width_minus1_mixed4to8_gap_to_2700k={direct_restoring_final_branch_final_low_path_width_minus1_mixed4to8_gap:.3}");
    println!("METRIC scratch600_direct_restoring_final_branch_final_selected_width_scan_mixed4to8_gap_to_2700k={direct_restoring_final_branch_final_selected_width_scan_mixed4to8_gap:.3}");
    println!("METRIC scratch600_direct_restoring_final_branch_final_low_path_width_scan_mixed4to8_gap_to_2700k={direct_restoring_final_branch_final_low_path_width_scan_mixed4to8_gap:.3}");
    println!("METRIC scratch600_direct_restoring_final_branch_final_low_path_width_minus1_scan_mixed4to8_gap_to_2700k={direct_restoring_final_branch_final_low_path_width_minus1_scan_mixed4to8_gap:.3}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_align_only_low_path_width_saving_mean={direct_restoring_final_low_branch_align_only_low_path_width_saving_mean:.3}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_align_only_model_precision_bits={direct_restoring_final_low_branch_align_only_model_precision_bits}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_align_only_raw_scratch_p99={direct_restoring_final_low_branch_align_only_raw_scratch_p99}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_align_only_raw_scratch_max={direct_restoring_final_low_branch_align_only_raw_scratch_max}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_align_only_step_entropy_scratch_p99={direct_restoring_final_low_branch_align_only_step_entropy_scratch_p99}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_align_only_step_entropy_scratch_max={direct_restoring_final_low_branch_align_only_step_entropy_scratch_max}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_align_only_step_prefix_scratch_p99={direct_restoring_final_low_branch_align_only_step_prefix_scratch_p99}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_align_only_step_prefix_scratch_max={direct_restoring_final_low_branch_align_only_step_prefix_scratch_max}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_align_only_best_block_symbols={direct_restoring_final_low_branch_align_only_best_block_symbols}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_align_only_best_touch_floor_mean={direct_restoring_final_low_branch_align_only_best_touch_floor_mean:.3}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_align_only_best_touch_floor_p99={direct_restoring_final_low_branch_align_only_best_touch_floor_p99}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_align_only_best_compressed_bits_p99={direct_restoring_final_low_branch_align_only_best_compressed_bits_p99}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_align_only_best_live_scratch_p99={direct_restoring_final_low_branch_align_only_best_live_scratch_p99}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_align_only_best_symbol_count_p99={direct_restoring_final_low_branch_align_only_best_symbol_count_p99}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_align_only_best_augmented_gap_to_2700k={direct_restoring_final_low_branch_align_only_best_augmented_gap:.3}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_align_only_mixed4to8_schedule_code={direct_restoring_final_low_branch_align_only_mixed4to8_schedule_code}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_align_only_mixed4to8_touch_floor_mean={direct_restoring_final_low_branch_align_only_mixed4to8_touch_floor_mean:.3}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_align_only_mixed4to8_live_scratch_p99={direct_restoring_final_low_branch_align_only_mixed4to8_live_scratch_p99}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_align_only_scan_lookup_floor_mean={direct_restoring_final_low_branch_align_only_scan_lookup_floor_mean:.3}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_align_only_scan_lookup_floor_p99={direct_restoring_final_low_branch_align_only_scan_lookup_floor_p99}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_align_only_binary_lookup_floor_mean={direct_restoring_final_low_branch_align_only_binary_lookup_floor_mean:.3}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_align_only_binary_lookup_floor_p99={direct_restoring_final_low_branch_align_only_binary_lookup_floor_p99}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_align_only_huffman_lookup_floor_mean={direct_restoring_final_low_branch_align_only_huffman_lookup_floor_mean:.3}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_align_only_huffman_lookup_floor_p99={direct_restoring_final_low_branch_align_only_huffman_lookup_floor_p99}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_align_only_prefix_tree_node_floor_mean={direct_restoring_final_low_branch_align_only_prefix_tree_node_floor_mean:.3}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_align_only_prefix_tree_node_floor_p99={direct_restoring_final_low_branch_align_only_prefix_tree_node_floor_p99}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_align_only_best_scan_gap_to_2700k={direct_restoring_final_low_branch_align_only_best_scan_gap:.3}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_align_only_best_binary_gap_to_2700k={direct_restoring_final_low_branch_align_only_best_binary_gap:.3}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_align_only_best_huffman_gap_to_2700k={direct_restoring_final_low_branch_align_only_best_huffman_gap:.3}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_align_only_best_prefix_tree_gap_to_2700k={direct_restoring_final_low_branch_align_only_best_prefix_tree_gap:.3}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_align_only_mixed4to8_scan_gap_to_2700k={direct_restoring_final_low_branch_align_only_mixed4to8_scan_gap:.3}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_align_only_mixed4to8_prefix_tree_gap_to_2700k={direct_restoring_final_low_branch_align_only_mixed4to8_prefix_tree_gap:.3}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_align_only_best_lookup_target_mean={direct_restoring_final_low_branch_align_only_best_lookup_target_mean:.3}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_align_only_best_lookup_multiplier_budget={direct_restoring_final_low_branch_align_only_best_lookup_multiplier_budget:.6}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_align_only_scan_over_binary_multiplier={direct_restoring_final_low_branch_align_only_scan_over_binary_multiplier:.6}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_align_only_huffman_over_binary_multiplier={direct_restoring_final_low_branch_align_only_huffman_over_binary_multiplier:.6}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_align_only_prefix_tree_over_binary_multiplier={direct_restoring_final_low_branch_align_only_prefix_tree_over_binary_multiplier:.6}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_align_only_support_noncontig_steps={direct_restoring_final_low_branch_align_only_support_noncontig_steps}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_align_only_support_max_span={direct_restoring_final_low_branch_align_only_support_max_span}");
    println!("METRIC scratch600_direct_restoring_final_prefix_bit_reader_toy_eq_ccx={direct_restoring_final_prefix_bit_reader_toy_eq_ccx}");
    println!("METRIC scratch600_direct_restoring_final_prefix_bit_reader_toy_dynamic_read_ccx={direct_restoring_final_prefix_bit_reader_toy_dynamic_read_ccx}");
    println!("METRIC scratch600_direct_restoring_final_prefix_bit_reader_toy_reader_forward_ccx={direct_restoring_final_prefix_bit_reader_toy_reader_forward_ccx}");
    println!("METRIC scratch600_direct_restoring_final_prefix_bit_reader_toy_tree_forward_ccx={direct_restoring_final_prefix_bit_reader_toy_tree_forward_ccx}");
    println!("METRIC scratch600_direct_restoring_final_prefix_bit_reader_toy_full_forward_ccx={direct_restoring_final_prefix_bit_reader_toy_full_forward_ccx}");
    println!("METRIC scratch600_direct_restoring_final_prefix_bit_reader_toy_roundtrip_ccx={direct_restoring_final_prefix_bit_reader_toy_roundtrip_ccx}");
    println!("METRIC scratch600_direct_restoring_final_prefix_bit_reader_toy_peak_q={direct_restoring_final_prefix_bit_reader_toy_peak_q}");
    println!("METRIC scratch600_direct_restoring_final_prefix_bit_reader_toy_cursor_states={direct_restoring_final_prefix_bit_reader_toy_cursor_states}");
    println!("METRIC scratch600_direct_restoring_final_prefix_bit_reader_toy_internal_nodes={direct_restoring_final_prefix_bit_reader_toy_internal_nodes}");
    println!("METRIC scratch600_direct_restoring_final_prefix_bit_reader_toy_reader_over_tree={direct_restoring_final_prefix_bit_reader_toy_reader_over_tree:.6}");
    println!("METRIC scratch600_direct_restoring_final_prefix_bit_reader_toy_tree_over_node_roundtrip={direct_restoring_final_prefix_bit_reader_toy_tree_over_node_roundtrip:.6}");
    println!("METRIC scratch600_direct_restoring_final_prefix_bit_reader_toy_full_over_node_roundtrip={direct_restoring_final_prefix_bit_reader_toy_full_over_node_roundtrip:.6}");
    println!("METRIC scratch600_direct_restoring_final_prefix_bit_reader_toy_roundtrip_ratio_budget={direct_restoring_final_prefix_bit_reader_toy_roundtrip_ratio_budget:.6}");
    println!("METRIC scratch600_direct_restoring_final_prefix_bit_reader_toy_tree_only_scaled_gap_to_2700k={direct_restoring_final_prefix_bit_reader_toy_tree_only_scaled_gap:.3}");
    println!("METRIC scratch600_direct_restoring_final_prefix_bit_reader_toy_cursor_scaled_gap_to_2700k={direct_restoring_final_prefix_bit_reader_toy_cursor_scaled_gap:.3}");
    println!("METRIC scratch600_direct_restoring_final_prefix_bit_reader_toy_dirty_restore_cases={direct_restoring_final_prefix_bit_reader_toy_dirty_restore_cases}");
    println!("METRIC scratch600_direct_restoring_final_prefix_bit_reader_toy_dirty_history_cases={direct_restoring_final_prefix_bit_reader_toy_dirty_history_cases}");
    println!("METRIC scratch600_direct_restoring_final_prefix_bit_reader_toy_dirty_phase_cases={direct_restoring_final_prefix_bit_reader_toy_dirty_phase_cases}");
    println!("METRIC scratch600_direct_restoring_final_prefix_cursor_advance_toy_ccx={direct_restoring_final_prefix_cursor_advance_toy_ccx}");
    println!("METRIC scratch600_direct_restoring_final_prefix_cursor_advance_toy_peak_q={direct_restoring_final_prefix_cursor_advance_toy_peak_q}");
    println!("METRIC scratch600_direct_restoring_final_prefix_cursor_advance_toy_combined_roundtrip_ccx={direct_restoring_final_prefix_cursor_advance_toy_combined_roundtrip_ccx}");
    println!("METRIC scratch600_direct_restoring_final_prefix_cursor_advance_toy_combined_over_node_roundtrip={direct_restoring_final_prefix_cursor_advance_toy_combined_over_node_roundtrip:.6}");
    println!("METRIC scratch600_direct_restoring_final_prefix_cursor_advance_toy_roundtrip_ratio_budget={direct_restoring_final_prefix_cursor_advance_toy_roundtrip_ratio_budget:.6}");
    println!("METRIC scratch600_direct_restoring_final_prefix_cursor_advance_toy_combined_scaled_gap_to_2700k={direct_restoring_final_prefix_cursor_advance_toy_combined_scaled_gap:.3}");
    println!("METRIC scratch600_direct_restoring_final_prefix_cursor_advance_toy_dirty_restore_cases={direct_restoring_final_prefix_cursor_advance_toy_dirty_restore_cases}");
    println!("METRIC scratch600_direct_restoring_final_prefix_cursor_advance_toy_dirty_history_cases={direct_restoring_final_prefix_cursor_advance_toy_dirty_history_cases}");
    println!("METRIC scratch600_direct_restoring_final_prefix_cursor_advance_toy_dirty_phase_cases={direct_restoring_final_prefix_cursor_advance_toy_dirty_phase_cases}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_toy_tree_ccx={direct_restoring_final_prefix_block2_toy_tree_ccx}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_toy_read2_ccx={direct_restoring_final_prefix_block2_toy_read2_ccx}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_toy_cursor_add_ccx={direct_restoring_final_prefix_block2_toy_cursor_add_ccx}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_toy_decode_forward_ccx={direct_restoring_final_prefix_block2_toy_decode_forward_ccx}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_toy_total_ccx={direct_restoring_final_prefix_block2_toy_total_ccx}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_toy_peak_q={direct_restoring_final_prefix_block2_toy_peak_q}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_toy_over_node_roundtrip={direct_restoring_final_prefix_block2_toy_over_node_roundtrip:.6}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_toy_roundtrip_ratio_budget={direct_restoring_final_prefix_block2_toy_roundtrip_ratio_budget:.6}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_toy_scaled_gap_to_2700k={direct_restoring_final_prefix_block2_toy_scaled_gap:.3}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_toy_dirty_restore_cases={direct_restoring_final_prefix_block2_toy_dirty_restore_cases}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_toy_dirty_history_cases={direct_restoring_final_prefix_block2_toy_dirty_history_cases}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_toy_dirty_phase_cases={direct_restoring_final_prefix_block2_toy_dirty_phase_cases}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_consume_toy_decode_forward_ccx={direct_restoring_final_prefix_block2_consume_toy_decode_forward_ccx}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_consume_toy_cursor_add_ccx={direct_restoring_final_prefix_block2_consume_toy_cursor_add_ccx}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_consume_toy_consume_ccx={direct_restoring_final_prefix_block2_consume_toy_consume_ccx}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_consume_toy_parser_transient_ccx={direct_restoring_final_prefix_block2_consume_toy_parser_transient_ccx}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_consume_toy_total_ccx={direct_restoring_final_prefix_block2_consume_toy_total_ccx}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_consume_toy_peak_q={direct_restoring_final_prefix_block2_consume_toy_peak_q}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_consume_toy_parser_over_node_roundtrip={direct_restoring_final_prefix_block2_consume_toy_parser_over_node_roundtrip:.6}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_consume_toy_roundtrip_ratio_budget={direct_restoring_final_prefix_block2_consume_toy_roundtrip_ratio_budget:.6}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_consume_toy_parser_scaled_gap_to_2700k={direct_restoring_final_prefix_block2_consume_toy_parser_scaled_gap:.3}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_consume_toy_dirty_restore_cases={direct_restoring_final_prefix_block2_consume_toy_dirty_restore_cases}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_consume_toy_dirty_history_cases={direct_restoring_final_prefix_block2_consume_toy_dirty_history_cases}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_consume_toy_dirty_phase_cases={direct_restoring_final_prefix_block2_consume_toy_dirty_phase_cases}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_leaf_touch_toy_decode_forward_ccx={direct_restoring_final_prefix_block2_leaf_touch_toy_decode_forward_ccx}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_leaf_touch_toy_leaf_touch_ccx={direct_restoring_final_prefix_block2_leaf_touch_toy_leaf_touch_ccx}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_leaf_touch_toy_parser_transient_ccx={direct_restoring_final_prefix_block2_leaf_touch_toy_parser_transient_ccx}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_leaf_touch_toy_total_ccx={direct_restoring_final_prefix_block2_leaf_touch_toy_total_ccx}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_leaf_touch_toy_peak_q={direct_restoring_final_prefix_block2_leaf_touch_toy_peak_q}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_leaf_touch_toy_parser_over_node_roundtrip={direct_restoring_final_prefix_block2_leaf_touch_toy_parser_over_node_roundtrip:.6}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_leaf_touch_toy_total_over_node_roundtrip={direct_restoring_final_prefix_block2_leaf_touch_toy_total_over_node_roundtrip:.6}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_leaf_touch_toy_roundtrip_ratio_budget={direct_restoring_final_prefix_block2_leaf_touch_toy_roundtrip_ratio_budget:.6}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_leaf_touch_toy_parser_scaled_gap_to_2700k={direct_restoring_final_prefix_block2_leaf_touch_toy_parser_scaled_gap:.3}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_leaf_touch_toy_total_scaled_gap_to_2700k={direct_restoring_final_prefix_block2_leaf_touch_toy_total_scaled_gap:.3}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_leaf_touch_toy_dirty_restore_cases={direct_restoring_final_prefix_block2_leaf_touch_toy_dirty_restore_cases}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_leaf_touch_toy_dirty_history_cases={direct_restoring_final_prefix_block2_leaf_touch_toy_dirty_history_cases}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_leaf_touch_toy_dirty_phase_cases={direct_restoring_final_prefix_block2_leaf_touch_toy_dirty_phase_cases}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_selected_addsub_toy_decode_forward_ccx={direct_restoring_final_prefix_block2_selected_addsub_toy_decode_forward_ccx}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_selected_addsub_toy_select_shift_ccx={direct_restoring_final_prefix_block2_selected_addsub_toy_select_shift_ccx}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_selected_addsub_toy_addsub_ccx={direct_restoring_final_prefix_block2_selected_addsub_toy_addsub_ccx}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_selected_addsub_toy_arithmetic_ccx={direct_restoring_final_prefix_block2_selected_addsub_toy_arithmetic_ccx}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_selected_addsub_toy_parser_transient_ccx={direct_restoring_final_prefix_block2_selected_addsub_toy_parser_transient_ccx}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_selected_addsub_toy_total_ccx={direct_restoring_final_prefix_block2_selected_addsub_toy_total_ccx}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_selected_addsub_toy_peak_q={direct_restoring_final_prefix_block2_selected_addsub_toy_peak_q}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_selected_addsub_toy_parser_over_node_roundtrip={direct_restoring_final_prefix_block2_selected_addsub_toy_parser_over_node_roundtrip:.6}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_selected_addsub_toy_arithmetic_over_node_roundtrip={direct_restoring_final_prefix_block2_selected_addsub_toy_arithmetic_over_node_roundtrip:.6}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_selected_addsub_toy_total_over_node_roundtrip={direct_restoring_final_prefix_block2_selected_addsub_toy_total_over_node_roundtrip:.6}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_selected_addsub_toy_roundtrip_ratio_budget={direct_restoring_final_prefix_block2_selected_addsub_toy_roundtrip_ratio_budget:.6}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_selected_addsub_toy_parser_scaled_gap_to_2700k={direct_restoring_final_prefix_block2_selected_addsub_toy_parser_scaled_gap:.3}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_selected_addsub_toy_total_scaled_gap_to_2700k={direct_restoring_final_prefix_block2_selected_addsub_toy_total_scaled_gap:.3}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_selected_addsub_toy_dirty_restore_cases={direct_restoring_final_prefix_block2_selected_addsub_toy_dirty_restore_cases}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_selected_addsub_toy_dirty_history_cases={direct_restoring_final_prefix_block2_selected_addsub_toy_dirty_history_cases}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_selected_addsub_toy_dirty_phase_cases={direct_restoring_final_prefix_block2_selected_addsub_toy_dirty_phase_cases}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_balanced_selected_addsub_toy_tree_ccx={direct_restoring_final_prefix_block2_balanced_selected_addsub_toy_tree_ccx}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_balanced_selected_addsub_toy_read2_ccx={direct_restoring_final_prefix_block2_balanced_selected_addsub_toy_read2_ccx}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_balanced_selected_addsub_toy_decode_forward_ccx={direct_restoring_final_prefix_block2_balanced_selected_addsub_toy_decode_forward_ccx}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_balanced_selected_addsub_toy_select_shift_ccx={direct_restoring_final_prefix_block2_balanced_selected_addsub_toy_select_shift_ccx}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_balanced_selected_addsub_toy_addsub_ccx={direct_restoring_final_prefix_block2_balanced_selected_addsub_toy_addsub_ccx}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_balanced_selected_addsub_toy_arithmetic_ccx={direct_restoring_final_prefix_block2_balanced_selected_addsub_toy_arithmetic_ccx}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_balanced_selected_addsub_toy_parser_transient_ccx={direct_restoring_final_prefix_block2_balanced_selected_addsub_toy_parser_transient_ccx}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_balanced_selected_addsub_toy_total_ccx={direct_restoring_final_prefix_block2_balanced_selected_addsub_toy_total_ccx}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_balanced_selected_addsub_toy_peak_q={direct_restoring_final_prefix_block2_balanced_selected_addsub_toy_peak_q}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_balanced_selected_addsub_toy_parser_over_node_roundtrip={direct_restoring_final_prefix_block2_balanced_selected_addsub_toy_parser_over_node_roundtrip:.6}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_balanced_selected_addsub_toy_arithmetic_over_node_roundtrip={direct_restoring_final_prefix_block2_balanced_selected_addsub_toy_arithmetic_over_node_roundtrip:.6}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_balanced_selected_addsub_toy_total_over_node_roundtrip={direct_restoring_final_prefix_block2_balanced_selected_addsub_toy_total_over_node_roundtrip:.6}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_balanced_selected_addsub_toy_roundtrip_ratio_budget={direct_restoring_final_prefix_block2_balanced_selected_addsub_toy_roundtrip_ratio_budget:.6}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_balanced_selected_addsub_toy_total_scaled_gap_to_2700k={direct_restoring_final_prefix_block2_balanced_selected_addsub_toy_total_scaled_gap:.3}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_balanced_selected_addsub_toy_dirty_restore_cases={direct_restoring_final_prefix_block2_balanced_selected_addsub_toy_dirty_restore_cases}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_balanced_selected_addsub_toy_dirty_history_cases={direct_restoring_final_prefix_block2_balanced_selected_addsub_toy_dirty_history_cases}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_balanced_selected_addsub_toy_dirty_phase_cases={direct_restoring_final_prefix_block2_balanced_selected_addsub_toy_dirty_phase_cases}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_selected_addsub_roundtrip_toy_forward_decode_ccx={direct_restoring_final_prefix_block2_selected_addsub_roundtrip_toy_forward_decode_ccx}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_selected_addsub_roundtrip_toy_reverse_decode_ccx={direct_restoring_final_prefix_block2_selected_addsub_roundtrip_toy_reverse_decode_ccx}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_selected_addsub_roundtrip_toy_forward_select_shift_ccx={direct_restoring_final_prefix_block2_selected_addsub_roundtrip_toy_forward_select_shift_ccx}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_selected_addsub_roundtrip_toy_reverse_select_shift_ccx={direct_restoring_final_prefix_block2_selected_addsub_roundtrip_toy_reverse_select_shift_ccx}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_selected_addsub_roundtrip_toy_forward_addsub_ccx={direct_restoring_final_prefix_block2_selected_addsub_roundtrip_toy_forward_addsub_ccx}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_selected_addsub_roundtrip_toy_reverse_addsub_ccx={direct_restoring_final_prefix_block2_selected_addsub_roundtrip_toy_reverse_addsub_ccx}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_selected_addsub_roundtrip_toy_forward_ccx={direct_restoring_final_prefix_block2_selected_addsub_roundtrip_toy_forward_ccx}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_selected_addsub_roundtrip_toy_reverse_ccx={direct_restoring_final_prefix_block2_selected_addsub_roundtrip_toy_reverse_ccx}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_selected_addsub_roundtrip_toy_arithmetic_ccx={direct_restoring_final_prefix_block2_selected_addsub_roundtrip_toy_arithmetic_ccx}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_selected_addsub_roundtrip_toy_parser_transient_ccx={direct_restoring_final_prefix_block2_selected_addsub_roundtrip_toy_parser_transient_ccx}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_selected_addsub_roundtrip_toy_total_ccx={direct_restoring_final_prefix_block2_selected_addsub_roundtrip_toy_total_ccx}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_selected_addsub_roundtrip_toy_peak_q={direct_restoring_final_prefix_block2_selected_addsub_roundtrip_toy_peak_q}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_selected_addsub_roundtrip_toy_parser_over_node_roundtrip={direct_restoring_final_prefix_block2_selected_addsub_roundtrip_toy_parser_over_node_roundtrip:.6}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_selected_addsub_roundtrip_toy_arithmetic_over_node_roundtrip={direct_restoring_final_prefix_block2_selected_addsub_roundtrip_toy_arithmetic_over_node_roundtrip:.6}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_selected_addsub_roundtrip_toy_total_over_node_roundtrip={direct_restoring_final_prefix_block2_selected_addsub_roundtrip_toy_total_over_node_roundtrip:.6}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_selected_addsub_roundtrip_toy_roundtrip_ratio_budget={direct_restoring_final_prefix_block2_selected_addsub_roundtrip_toy_roundtrip_ratio_budget:.6}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_selected_addsub_roundtrip_toy_parser_scaled_gap_to_2700k={direct_restoring_final_prefix_block2_selected_addsub_roundtrip_toy_parser_scaled_gap:.3}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_selected_addsub_roundtrip_toy_total_scaled_gap_to_2700k={direct_restoring_final_prefix_block2_selected_addsub_roundtrip_toy_total_scaled_gap:.3}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_selected_addsub_roundtrip_toy_dirty_restore_cases={direct_restoring_final_prefix_block2_selected_addsub_roundtrip_toy_dirty_restore_cases}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_selected_addsub_roundtrip_toy_dirty_history_cases={direct_restoring_final_prefix_block2_selected_addsub_roundtrip_toy_dirty_history_cases}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_selected_addsub_roundtrip_toy_dirty_phase_cases={direct_restoring_final_prefix_block2_selected_addsub_roundtrip_toy_dirty_phase_cases}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_span24_roundtrip_toy_forward_decode_ccx={direct_restoring_final_prefix_block2_span24_roundtrip_toy_forward_decode_ccx}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_span24_roundtrip_toy_reverse_decode_ccx={direct_restoring_final_prefix_block2_span24_roundtrip_toy_reverse_decode_ccx}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_span24_roundtrip_toy_forward_select_shift_ccx={direct_restoring_final_prefix_block2_span24_roundtrip_toy_forward_select_shift_ccx}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_span24_roundtrip_toy_reverse_select_shift_ccx={direct_restoring_final_prefix_block2_span24_roundtrip_toy_reverse_select_shift_ccx}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_span24_roundtrip_toy_forward_addsub_ccx={direct_restoring_final_prefix_block2_span24_roundtrip_toy_forward_addsub_ccx}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_span24_roundtrip_toy_reverse_addsub_ccx={direct_restoring_final_prefix_block2_span24_roundtrip_toy_reverse_addsub_ccx}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_span24_roundtrip_toy_total_ccx={direct_restoring_final_prefix_block2_span24_roundtrip_toy_total_ccx}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_span24_roundtrip_toy_peak_q={direct_restoring_final_prefix_block2_span24_roundtrip_toy_peak_q}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_span24_roundtrip_toy_parser_over_node_roundtrip={direct_restoring_final_prefix_block2_span24_roundtrip_toy_parser_over_node_roundtrip:.6}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_span24_roundtrip_toy_arithmetic_over_node_roundtrip={direct_restoring_final_prefix_block2_span24_roundtrip_toy_arithmetic_over_node_roundtrip:.6}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_span24_roundtrip_toy_total_over_node_roundtrip={direct_restoring_final_prefix_block2_span24_roundtrip_toy_total_over_node_roundtrip:.6}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_span24_roundtrip_toy_roundtrip_ratio_budget={direct_restoring_final_prefix_block2_span24_roundtrip_toy_roundtrip_ratio_budget:.6}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_span24_roundtrip_toy_total_scaled_gap_to_2700k={direct_restoring_final_prefix_block2_span24_roundtrip_toy_total_scaled_gap:.3}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_span24_taper_materialized_full_add_per_digit={direct_restoring_final_prefix_block2_span24_taper_materialized_full_add_per_digit}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_span24_taper_add_per_digit_floor={direct_restoring_final_prefix_block2_span24_taper_add_per_digit_floor}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_span24_taper_arithmetic_floor={direct_restoring_final_prefix_block2_span24_taper_arithmetic_floor}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_span24_taper_total_floor={direct_restoring_final_prefix_block2_span24_taper_total_floor}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_span24_taper_total_over_node_roundtrip={direct_restoring_final_prefix_block2_span24_taper_total_over_node_roundtrip:.6}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_span24_taper_scaled_gap_to_2700k={direct_restoring_final_prefix_block2_span24_taper_scaled_gap:.3}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_span24_roundtrip_toy_dirty_restore_cases={direct_restoring_final_prefix_block2_span24_roundtrip_toy_dirty_restore_cases}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_span24_roundtrip_toy_dirty_history_cases={direct_restoring_final_prefix_block2_span24_roundtrip_toy_dirty_history_cases}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_span24_roundtrip_toy_dirty_phase_cases={direct_restoring_final_prefix_block2_span24_roundtrip_toy_dirty_phase_cases}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_prefix_support_weighted_prefix_node_mean={direct_restoring_final_low_branch_prefix_support_weighted_prefix_node_mean:.3}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_prefix_support_weighted_prefix_node_p99={direct_restoring_final_low_branch_prefix_support_weighted_prefix_node_p99}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_prefix_support_weighted_materialized_digit_mean={direct_restoring_final_low_branch_prefix_support_weighted_materialized_digit_mean:.3}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_prefix_support_weighted_materialized_digit_p99={direct_restoring_final_low_branch_prefix_support_weighted_materialized_digit_p99}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_prefix_support_weighted_tree_decode_mean={direct_restoring_final_low_branch_prefix_support_weighted_tree_decode_mean:.3}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_prefix_support_weighted_dynamic_even_mean={direct_restoring_final_low_branch_prefix_support_weighted_dynamic_even_mean:.3}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_prefix_support_weighted_variable_decode_mean={direct_restoring_final_low_branch_prefix_support_weighted_variable_decode_mean:.3}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_prefix_support_weighted_variable_offset1_decode_mean={direct_restoring_final_low_branch_prefix_support_weighted_variable_offset1_decode_mean:.3}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_prefix_support_weighted_arithmetic_over_node_roundtrip={direct_restoring_final_low_branch_prefix_support_weighted_arithmetic_over_node_roundtrip:.6}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_prefix_support_weighted_total_over_node_roundtrip={direct_restoring_final_low_branch_prefix_support_weighted_total_over_node_roundtrip:.6}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_prefix_support_weighted_variable_total_over_node_roundtrip={direct_restoring_final_low_branch_prefix_support_weighted_variable_total_over_node_roundtrip:.6}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_prefix_support_weighted_variable_offset1_total_over_node_roundtrip={direct_restoring_final_low_branch_prefix_support_weighted_variable_offset1_total_over_node_roundtrip:.6}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_prefix_support_weighted_balanced_total_over_node_roundtrip={direct_restoring_final_low_branch_prefix_support_weighted_balanced_total_over_node_roundtrip:.6}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_prefix_support_weighted_selective_total_over_node_roundtrip={direct_restoring_final_low_branch_prefix_support_weighted_selective_total_over_node_roundtrip:.6}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_prefix_support_weighted_ratio_budget={direct_restoring_final_low_branch_prefix_support_weighted_ratio_budget:.6}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_prefix_support_weighted_gap_to_2700k={direct_restoring_final_low_branch_prefix_support_weighted_gap:.3}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_prefix_support_weighted_variable_gap_to_2700k={direct_restoring_final_low_branch_prefix_support_weighted_variable_gap:.3}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_prefix_support_weighted_variable_offset1_gap_to_2700k={direct_restoring_final_low_branch_prefix_support_weighted_variable_offset1_gap:.3}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_prefix_support_weighted_balanced_gap_to_2700k={direct_restoring_final_low_branch_prefix_support_weighted_balanced_gap:.3}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_prefix_support_weighted_selective_gap_to_2700k={direct_restoring_final_low_branch_prefix_support_weighted_selective_gap:.3}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_prefix_support_weighted_projected_toffoli={direct_restoring_final_low_branch_prefix_support_weighted_projected_toffoli:.3}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_prefix_support_weighted_variable_projected_toffoli={direct_restoring_final_low_branch_prefix_support_weighted_variable_projected_toffoli:.3}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_prefix_support_weighted_variable_offset1_projected_toffoli={direct_restoring_final_low_branch_prefix_support_weighted_variable_offset1_projected_toffoli:.3}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_prefix_support_weighted_balanced_projected_toffoli={direct_restoring_final_low_branch_prefix_support_weighted_balanced_projected_toffoli:.3}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_prefix_support_weighted_selective_projected_toffoli={direct_restoring_final_low_branch_prefix_support_weighted_selective_projected_toffoli:.3}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_prefix_support_weighted_shannon_prefix_bit_p99={direct_restoring_final_low_branch_prefix_support_weighted_shannon_prefix_bit_p99}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_prefix_support_weighted_balanced_prefix_bit_p99={direct_restoring_final_low_branch_prefix_support_weighted_balanced_prefix_bit_p99}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_prefix_support_weighted_selective_prefix_bit_p99={direct_restoring_final_low_branch_prefix_support_weighted_selective_prefix_bit_p99}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_prefix_support_weighted_selective_prefix_bit_max={direct_restoring_final_low_branch_prefix_support_weighted_selective_prefix_bit_max}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_prefix_support_weighted_selective_prefix_scratch_p99={direct_restoring_final_low_branch_prefix_support_weighted_selective_prefix_scratch_p99}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_prefix_support_weighted_selective_flatten_steps={direct_restoring_final_low_branch_prefix_support_weighted_selective_flatten_steps}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_prefix_support_weighted_maxconstrained_prefix_bit_mean={direct_restoring_final_low_branch_prefix_support_weighted_maxconstrained_prefix_bit_mean:.3}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_prefix_support_weighted_maxconstrained_prefix_bit_p99={direct_restoring_final_low_branch_prefix_support_weighted_maxconstrained_prefix_bit_p99}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_prefix_support_weighted_maxconstrained_prefix_bit_max={direct_restoring_final_low_branch_prefix_support_weighted_maxconstrained_prefix_bit_max}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_prefix_support_weighted_maxconstrained_prefix_scratch_max={direct_restoring_final_low_branch_prefix_support_weighted_maxconstrained_prefix_scratch_max}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_prefix_support_weighted_maxconstrained_flatten_steps={direct_restoring_final_low_branch_prefix_support_weighted_maxconstrained_flatten_steps}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_prefix_support_weighted_maxconstrained_trimmed_steps={direct_restoring_final_low_branch_prefix_support_weighted_maxconstrained_trimmed_steps}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_prefix_support_weighted_maxconstrained_over_budget_rows={direct_restoring_final_low_branch_prefix_support_weighted_maxconstrained_over_budget_rows}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_prefix_support_weighted_maxconstrained_over_budget_mass={direct_restoring_final_low_branch_prefix_support_weighted_maxconstrained_over_budget_mass}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_prefix_support_weighted_selective_schedule_codebook_steps={direct_restoring_final_low_branch_prefix_support_weighted_selective_schedule_codebook_steps}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_prefix_support_weighted_selective_schedule_max_code_len={direct_restoring_final_low_branch_prefix_support_weighted_selective_schedule_max_code_len}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_prefix_support_weighted_selective_schedule_max_len_classes={direct_restoring_final_low_branch_prefix_support_weighted_selective_schedule_max_len_classes}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_prefix_support_weighted_selective_schedule_bit_mean={direct_restoring_final_low_branch_prefix_support_weighted_selective_schedule_bit_mean:.3}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_prefix_support_weighted_selective_schedule_bit_p99={direct_restoring_final_low_branch_prefix_support_weighted_selective_schedule_bit_p99}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_prefix_support_weighted_selective_schedule_max_bits={direct_restoring_final_low_branch_prefix_support_weighted_selective_schedule_max_bits}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_prefix_support_weighted_selective_schedule_decoded_symbols={direct_restoring_final_low_branch_prefix_support_weighted_selective_schedule_decoded_symbols}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_prefix_support_weighted_selective_schedule_prefix_collisions={direct_restoring_final_low_branch_prefix_support_weighted_selective_schedule_prefix_collisions}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_prefix_support_weighted_selective_schedule_decode_mismatches={direct_restoring_final_low_branch_prefix_support_weighted_selective_schedule_decode_mismatches}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_prefix_support_weighted_selective_schedule_cursor_mismatches={direct_restoring_final_low_branch_prefix_support_weighted_selective_schedule_cursor_mismatches}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_prefix_support_weighted_maxconstrained_schedule_codebook_steps={direct_restoring_final_low_branch_prefix_support_weighted_maxconstrained_schedule_codebook_steps}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_prefix_support_weighted_maxconstrained_schedule_max_code_len={direct_restoring_final_low_branch_prefix_support_weighted_maxconstrained_schedule_max_code_len}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_prefix_support_weighted_maxconstrained_schedule_max_len_classes={direct_restoring_final_low_branch_prefix_support_weighted_maxconstrained_schedule_max_len_classes}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_prefix_support_weighted_maxconstrained_schedule_bit_mean={direct_restoring_final_low_branch_prefix_support_weighted_maxconstrained_schedule_bit_mean:.3}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_prefix_support_weighted_maxconstrained_schedule_bit_p99={direct_restoring_final_low_branch_prefix_support_weighted_maxconstrained_schedule_bit_p99}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_prefix_support_weighted_maxconstrained_schedule_max_bits={direct_restoring_final_low_branch_prefix_support_weighted_maxconstrained_schedule_max_bits}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_prefix_support_weighted_maxconstrained_schedule_decoded_symbols={direct_restoring_final_low_branch_prefix_support_weighted_maxconstrained_schedule_decoded_symbols}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_prefix_support_weighted_maxconstrained_schedule_prefix_collisions={direct_restoring_final_low_branch_prefix_support_weighted_maxconstrained_schedule_prefix_collisions}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_prefix_support_weighted_maxconstrained_schedule_decode_mismatches={direct_restoring_final_low_branch_prefix_support_weighted_maxconstrained_schedule_decode_mismatches}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_prefix_support_weighted_maxconstrained_schedule_cursor_mismatches={direct_restoring_final_low_branch_prefix_support_weighted_maxconstrained_schedule_cursor_mismatches}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_prefix_support_weighted_selective_dynamic_even_mean={direct_restoring_final_low_branch_prefix_support_weighted_selective_dynamic_even_mean:.3}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_prefix_support_weighted_selective_variable_decode_mean={direct_restoring_final_low_branch_prefix_support_weighted_selective_variable_decode_mean:.3}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_prefix_support_weighted_maxconstrained_dynamic_even_mean={direct_restoring_final_low_branch_prefix_support_weighted_maxconstrained_dynamic_even_mean:.3}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_prefix_support_weighted_maxconstrained_dynamic_even_p99={direct_restoring_final_low_branch_prefix_support_weighted_maxconstrained_dynamic_even_p99}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_prefix_support_weighted_maxconstrained_variable_decode_mean={direct_restoring_final_low_branch_prefix_support_weighted_maxconstrained_variable_decode_mean:.3}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_prefix_support_weighted_maxconstrained_variable_decode_p99={direct_restoring_final_low_branch_prefix_support_weighted_maxconstrained_variable_decode_p99}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_prefix_support_weighted_maxconstrained_total_over_node_roundtrip={direct_restoring_final_low_branch_prefix_support_weighted_maxconstrained_total_over_node_roundtrip:.6}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_prefix_support_weighted_maxconstrained_gap_to_2700k={direct_restoring_final_low_branch_prefix_support_weighted_maxconstrained_gap:.3}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_prefix_support_weighted_maxconstrained_projected_toffoli={direct_restoring_final_low_branch_prefix_support_weighted_maxconstrained_projected_toffoli:.3}");
    println!("METRIC scratch600_direct_restoring_final_peakfit_toy_cases_with_sample_gap={direct_restoring_final_peakfit_toy_cases_with_sample_gap}");
    println!("METRIC scratch600_direct_restoring_final_peakfit_toy_largest_missing_symbols={direct_restoring_final_peakfit_toy_largest_missing_symbols}");
    println!("METRIC scratch600_direct_restoring_final_peakfit_toy_largest_sample_over_budget_traces={direct_restoring_final_peakfit_toy_largest_sample_over_budget_traces}");
    println!("METRIC scratch600_direct_restoring_final_peakfit_toy_largest_exact_over_budget_traces={direct_restoring_final_peakfit_toy_largest_exact_over_budget_traces}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_support_toy_cases_with_missing={direct_restoring_final_low_branch_support_toy_cases_with_missing}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_support_toy_largest_missing_symbols={direct_restoring_final_low_branch_support_toy_largest_missing_symbols}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_support_toy_largest_missing_steps={direct_restoring_final_low_branch_support_toy_largest_missing_steps}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_support_toy_largest_span_gap={direct_restoring_final_low_branch_support_toy_largest_span_gap}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_support_toy_largest_exact_span={direct_restoring_final_low_branch_support_toy_largest_exact_span}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_interval_toy_guard4_cover_cases={direct_restoring_final_low_branch_interval_toy_guard4_cover_cases}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_interval_toy_guard4_largest_missing_symbols={direct_restoring_final_low_branch_interval_toy_guard4_largest_missing_symbols}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_interval_toy_guard4_largest_over_budget_traces={direct_restoring_final_low_branch_interval_toy_guard4_largest_over_budget_traces}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_interval_toy_guard4_largest_max_bits={direct_restoring_final_low_branch_interval_toy_guard4_largest_max_bits}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_interval_toy_full_cover_cases={direct_restoring_final_low_branch_interval_toy_full_cover_cases}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_interval_toy_full_fit_cases={direct_restoring_final_low_branch_interval_toy_full_fit_cases}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_interval_toy_full_largest_over_budget_traces={direct_restoring_final_low_branch_interval_toy_full_largest_over_budget_traces}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_interval_toy_full_largest_max_bits={direct_restoring_final_low_branch_interval_toy_full_largest_max_bits}");
    println!("METRIC scratch600_direct_restoring_final_peakfit_holdout_missing_symbols={direct_restoring_final_peakfit_holdout_missing_symbols}");
    println!("METRIC scratch600_direct_restoring_final_peakfit_holdout_missing_traces={direct_restoring_final_peakfit_holdout_missing_traces}");
    println!("METRIC scratch600_direct_restoring_final_peakfit_holdout_over_budget_rows={direct_restoring_final_peakfit_holdout_over_budget_rows}");
    println!("METRIC scratch600_direct_restoring_final_peakfit_holdout_max_seen_bits={direct_restoring_final_peakfit_holdout_max_seen_bits}");
    println!("METRIC scratch600_direct_restoring_final_peakfit_scaled_probe_train_samples={direct_restoring_final_peakfit_scaled_probe_train_samples}");
    println!("METRIC scratch600_direct_restoring_final_peakfit_scaled_probe_holdout_samples={direct_restoring_final_peakfit_scaled_probe_holdout_samples}");
    println!("METRIC scratch600_direct_restoring_final_peakfit_scaled_probe_flatten_steps={direct_restoring_final_peakfit_scaled_probe_flatten_steps}");
    println!("METRIC scratch600_direct_restoring_final_peakfit_scaled_probe_missing_symbols={direct_restoring_final_peakfit_scaled_probe_missing_symbols}");
    println!("METRIC scratch600_direct_restoring_final_peakfit_scaled_probe_missing_traces={direct_restoring_final_peakfit_scaled_probe_missing_traces}");
    println!("METRIC scratch600_direct_restoring_final_peakfit_scaled_probe_over_budget_rows={direct_restoring_final_peakfit_scaled_probe_over_budget_rows}");
    println!("METRIC scratch600_direct_restoring_final_peakfit_scaled_probe_max_seen_bits={direct_restoring_final_peakfit_scaled_probe_max_seen_bits}");
    println!("METRIC scratch600_direct_restoring_final_peakfit_scaled_probe_gap_to_2700k={direct_restoring_final_peakfit_scaled_probe_gap:.3}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_prefix_support_weighted_span24_uniform_gap_to_2700k={direct_restoring_final_low_branch_prefix_support_weighted_span24_uniform_gap:.3}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_prefix_support_weighted_span24_symbol_mean={direct_restoring_final_low_branch_prefix_support_weighted_span24_symbol_mean:.3}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_prefix_support_weighted_span24_symbol_p99={direct_restoring_final_low_branch_prefix_support_weighted_span24_symbol_p99}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_prefix_support_weighted_support_noncontig_steps={direct_restoring_final_low_branch_prefix_support_weighted_support_noncontig_steps}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_prefix_support_weighted_support_max_span={direct_restoring_final_low_branch_prefix_support_weighted_support_max_span}");
    println!("METRIC scratch600_direct_restoring_final_low_branch_prefix_support_weighted_support_max_symbols={direct_restoring_final_low_branch_prefix_support_weighted_support_max_symbols}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_balanced_family_toy_checked_circuits={direct_restoring_final_prefix_block2_balanced_family_toy_checked_circuits}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_balanced_family_toy_simulated_circuits={direct_restoring_final_prefix_block2_balanced_family_toy_simulated_circuits}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_balanced_family_toy_simulated_cases={direct_restoring_final_prefix_block2_balanced_family_toy_simulated_cases}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_balanced_family_toy_max_support={direct_restoring_final_prefix_block2_balanced_family_toy_max_support}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_balanced_family_toy_max_tree_ccx={direct_restoring_final_prefix_block2_balanced_family_toy_max_tree_ccx}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_balanced_family_toy_max_read2_ccx={direct_restoring_final_prefix_block2_balanced_family_toy_max_read2_ccx}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_balanced_family_toy_max_decode_forward_ccx={direct_restoring_final_prefix_block2_balanced_family_toy_max_decode_forward_ccx}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_balanced_family_toy_max_select_shift_ccx={direct_restoring_final_prefix_block2_balanced_family_toy_max_select_shift_ccx}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_balanced_family_toy_max_addsub_ccx={direct_restoring_final_prefix_block2_balanced_family_toy_max_addsub_ccx}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_balanced_family_toy_max_total_ccx={direct_restoring_final_prefix_block2_balanced_family_toy_max_total_ccx}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_balanced_family_toy_max_peak_q={direct_restoring_final_prefix_block2_balanced_family_toy_max_peak_q}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_balanced_family_toy_max_total_over_node_roundtrip={direct_restoring_final_prefix_block2_balanced_family_toy_max_total_over_node_roundtrip:.6}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_balanced_family_toy_max_ratio_support0={direct_restoring_final_prefix_block2_balanced_family_toy_max_ratio_support0}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_balanced_family_toy_max_ratio_support1={direct_restoring_final_prefix_block2_balanced_family_toy_max_ratio_support1}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_balanced_family_toy_max_total_scaled_gap_to_2700k={direct_restoring_final_prefix_block2_balanced_family_toy_max_total_scaled_gap:.3}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_balanced_family_toy_dirty_restore_cases={direct_restoring_final_prefix_block2_balanced_family_toy_dirty_restore_cases}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_balanced_family_toy_dirty_history_cases={direct_restoring_final_prefix_block2_balanced_family_toy_dirty_history_cases}");
    println!("METRIC scratch600_direct_restoring_final_prefix_block2_balanced_family_toy_dirty_phase_cases={direct_restoring_final_prefix_block2_balanced_family_toy_dirty_phase_cases}");
    println!("METRIC scratch600_direct_restoring_final_coeff_decoder_alignment_degree_n14={direct_restoring_final_coeff_decoder_alignment_degree_n14}");
    println!("METRIC scratch600_direct_restoring_final_coeff_decoder_alignment_density_n14={direct_restoring_final_coeff_decoder_alignment_density_n14}");
    println!("METRIC scratch600_direct_restoring_final_coeff_decoder_alignment_max_n14={direct_restoring_final_coeff_decoder_alignment_max_n14}");
    println!("METRIC scratch600_direct_restoring_final_align_entropy_variable_scratch_p99={direct_restoring_final_align_entropy_variable_scratch_p99}");
    println!("METRIC scratch600_direct_restoring_final_align_entropy_variable_scratch_max={direct_restoring_final_align_entropy_variable_scratch_max}");
    println!("METRIC scratch600_direct_restoring_final_align_entropy_global_scratch_p99={direct_restoring_final_align_entropy_global_scratch_p99}");
    println!("METRIC scratch600_direct_restoring_final_align_entropy_global_scratch_max={direct_restoring_final_align_entropy_global_scratch_max}");
    println!("METRIC scratch600_direct_restoring_final_align_entropy_step_scratch_p99={direct_restoring_final_align_entropy_step_scratch_p99}");
    println!("METRIC scratch600_direct_restoring_final_align_entropy_step_scratch_max={direct_restoring_final_align_entropy_step_scratch_max}");
    println!("METRIC scratch600_direct_restoring_final_align_prefix_global_scratch_p99={direct_restoring_final_align_prefix_global_scratch_p99}");
    println!("METRIC scratch600_direct_restoring_final_align_prefix_global_scratch_max={direct_restoring_final_align_prefix_global_scratch_max}");
    println!("METRIC scratch600_direct_restoring_final_align_prefix_step_scratch_p99={direct_restoring_final_align_prefix_step_scratch_p99}");
    println!("METRIC scratch600_direct_restoring_final_align_prefix_step_scratch_max={direct_restoring_final_align_prefix_step_scratch_max}");
    println!("METRIC scratch600_direct_restoring_final_align_entropy_branch_count_p99={direct_restoring_final_align_entropy_branch_count_p99}");
    println!("METRIC scratch600_direct_restoring_final_align_entropy_branch_count_max={direct_restoring_final_align_entropy_branch_count_max}");
    println!("METRIC scratch600_direct_restoring_final_range_parser_model_precision_bits={direct_restoring_final_range_parser_model_precision_bits}");
    println!("METRIC scratch600_direct_restoring_final_range_parser_state_bits_p99={direct_restoring_final_range_parser_state_bits_p99}");
    println!("METRIC scratch600_direct_restoring_final_range_parser_live_scratch_p99={direct_restoring_final_range_parser_live_scratch_p99}");
    println!("METRIC scratch600_direct_restoring_final_range_parser_symbol_count_p99={direct_restoring_final_range_parser_symbol_count_p99}");
    println!("METRIC scratch600_direct_restoring_final_range_parser_state_touch_floor_mean={direct_restoring_final_range_parser_state_touch_floor_mean}");
    println!("METRIC scratch600_direct_restoring_final_range_parser_state_touch_floor_p99={direct_restoring_final_range_parser_state_touch_floor_p99}");
    println!("METRIC scratch600_direct_restoring_final_range_parser_oneway_budget={direct_restoring_final_range_parser_oneway_budget}");
    println!("METRIC scratch600_direct_restoring_final_range_parser_augmented_mean_gap_to_2700k={direct_restoring_final_range_parser_augmented_mean_gap}");
    println!("METRIC scratch600_direct_restoring_final_block_parser_model_precision_bits={direct_restoring_final_block_parser_model_precision_bits}");
    println!("METRIC scratch600_direct_restoring_final_block_parser_oneway_budget={direct_restoring_final_block_parser_oneway_budget:.3}");
    println!("METRIC scratch600_direct_restoring_final_block_parser_best_block_symbols={direct_restoring_final_block_parser_best_block_symbols}");
    println!("METRIC scratch600_direct_restoring_final_block_parser_best_touch_floor_mean={direct_restoring_final_block_parser_best_touch_floor_mean:.3}");
    println!("METRIC scratch600_direct_restoring_final_block_parser_best_touch_floor_p99={direct_restoring_final_block_parser_best_touch_floor_p99}");
    println!("METRIC scratch600_direct_restoring_final_block_parser_best_compressed_bits_p99={direct_restoring_final_block_parser_best_compressed_bits_p99}");
    println!("METRIC scratch600_direct_restoring_final_block_parser_best_live_scratch_p99={direct_restoring_final_block_parser_best_live_scratch_p99}");
    println!("METRIC scratch600_direct_restoring_final_block_parser_best_symbol_count_p99={direct_restoring_final_block_parser_best_symbol_count_p99}");
    println!("METRIC scratch600_direct_restoring_final_block_parser_best_augmented_gap_to_2700k={direct_restoring_final_block_parser_best_augmented_gap:.3}");
    println!("METRIC scratch600_direct_restoring_final_block32_touch_floor_mean={direct_restoring_final_block32_touch_floor_mean:.3}");
    println!("METRIC scratch600_direct_restoring_final_block32_touch_floor_p99={direct_restoring_final_block32_touch_floor_p99}");
    println!("METRIC scratch600_direct_restoring_final_block32_compressed_bits_p99={direct_restoring_final_block32_compressed_bits_p99}");
    println!("METRIC scratch600_direct_restoring_final_block32_live_scratch_p99={direct_restoring_final_block32_live_scratch_p99}");
    println!("METRIC scratch600_direct_restoring_final_block32_symbol_count_p99={direct_restoring_final_block32_symbol_count_p99}");
    println!("METRIC scratch600_direct_restoring_final_block32_augmented_gap_to_2700k={direct_restoring_final_block32_augmented_gap:.3}");
    println!("METRIC scratch600_direct_restoring_final_block4_touch_floor_mean={direct_restoring_final_block4_touch_floor_mean:.3}");
    println!("METRIC scratch600_direct_restoring_final_block4_touch_floor_p99={direct_restoring_final_block4_touch_floor_p99}");
    println!("METRIC scratch600_direct_restoring_final_block4_compressed_bits_p99={direct_restoring_final_block4_compressed_bits_p99}");
    println!("METRIC scratch600_direct_restoring_final_block4_live_scratch_p99={direct_restoring_final_block4_live_scratch_p99}");
    println!("METRIC scratch600_direct_restoring_final_block4_symbol_count_p99={direct_restoring_final_block4_symbol_count_p99}");
    println!("METRIC scratch600_direct_restoring_final_block4_augmented_gap_to_2700k={direct_restoring_final_block4_augmented_gap:.3}");
    println!("METRIC scratch600_direct_restoring_final_block5_touch_floor_mean={direct_restoring_final_block5_touch_floor_mean:.3}");
    println!("METRIC scratch600_direct_restoring_final_block5_touch_floor_p99={direct_restoring_final_block5_touch_floor_p99}");
    println!("METRIC scratch600_direct_restoring_final_block5_compressed_bits_p99={direct_restoring_final_block5_compressed_bits_p99}");
    println!("METRIC scratch600_direct_restoring_final_block5_live_scratch_p99={direct_restoring_final_block5_live_scratch_p99}");
    println!("METRIC scratch600_direct_restoring_final_block5_symbol_count_p99={direct_restoring_final_block5_symbol_count_p99}");
    println!("METRIC scratch600_direct_restoring_final_block5_augmented_gap_to_2700k={direct_restoring_final_block5_augmented_gap:.3}");
    println!("METRIC scratch600_direct_restoring_final_block6_touch_floor_mean={direct_restoring_final_block6_touch_floor_mean:.3}");
    println!("METRIC scratch600_direct_restoring_final_block6_touch_floor_p99={direct_restoring_final_block6_touch_floor_p99}");
    println!("METRIC scratch600_direct_restoring_final_block6_compressed_bits_p99={direct_restoring_final_block6_compressed_bits_p99}");
    println!("METRIC scratch600_direct_restoring_final_block6_live_scratch_p99={direct_restoring_final_block6_live_scratch_p99}");
    println!("METRIC scratch600_direct_restoring_final_block6_symbol_count_p99={direct_restoring_final_block6_symbol_count_p99}");
    println!("METRIC scratch600_direct_restoring_final_block6_augmented_gap_to_2700k={direct_restoring_final_block6_augmented_gap:.3}");
    println!("METRIC scratch600_direct_restoring_final_block7_touch_floor_mean={direct_restoring_final_block7_touch_floor_mean:.3}");
    println!("METRIC scratch600_direct_restoring_final_block7_touch_floor_p99={direct_restoring_final_block7_touch_floor_p99}");
    println!("METRIC scratch600_direct_restoring_final_block7_compressed_bits_p99={direct_restoring_final_block7_compressed_bits_p99}");
    println!("METRIC scratch600_direct_restoring_final_block7_live_scratch_p99={direct_restoring_final_block7_live_scratch_p99}");
    println!("METRIC scratch600_direct_restoring_final_block7_symbol_count_p99={direct_restoring_final_block7_symbol_count_p99}");
    println!("METRIC scratch600_direct_restoring_final_block7_augmented_gap_to_2700k={direct_restoring_final_block7_augmented_gap:.3}");
    println!("METRIC scratch600_direct_restoring_final_block_parser_best_qrom_row_floor={direct_restoring_final_block_parser_best_qrom_row_floor}");
    println!("METRIC scratch600_direct_restoring_final_block_parser_best_qrom_max_rows_in_block={direct_restoring_final_block_parser_best_qrom_max_rows_in_block}");
    println!("METRIC scratch600_direct_restoring_final_block_parser_best_qrom_block_count_p99={direct_restoring_final_block_parser_best_qrom_block_count_p99}");
    println!("METRIC scratch600_direct_restoring_final_block_parser_best_qrom_gap_to_2700k={direct_restoring_final_block_parser_best_qrom_gap:.3}");
    println!("METRIC scratch600_direct_restoring_final_block32_qrom_row_floor={direct_restoring_final_block32_qrom_row_floor}");
    println!("METRIC scratch600_direct_restoring_final_block32_qrom_max_rows_in_block={direct_restoring_final_block32_qrom_max_rows_in_block}");
    println!("METRIC scratch600_direct_restoring_final_block32_qrom_block_count_p99={direct_restoring_final_block32_qrom_block_count_p99}");
    println!("METRIC scratch600_direct_restoring_final_block_parser_lookup_scan_floor_mean={direct_restoring_final_block_parser_lookup_scan_floor_mean:.3}");
    println!("METRIC scratch600_direct_restoring_final_block_parser_lookup_scan_floor_p99={direct_restoring_final_block_parser_lookup_scan_floor_p99}");
    println!("METRIC scratch600_direct_restoring_final_block_parser_cond_branch_lookup_scan_floor_mean={direct_restoring_final_block_parser_cond_branch_lookup_scan_floor_mean:.3}");
    println!("METRIC scratch600_direct_restoring_final_block_parser_cond_branch_lookup_scan_floor_p99={direct_restoring_final_block_parser_cond_branch_lookup_scan_floor_p99}");
    println!("METRIC scratch600_direct_restoring_final_block_parser_best_with_lookup_mean={direct_restoring_final_block_parser_best_with_lookup_mean:.3}");
    println!("METRIC scratch600_direct_restoring_final_block_parser_best_with_lookup_gap_to_2700k={direct_restoring_final_block_parser_best_with_lookup_gap:.3}");
    println!("METRIC scratch600_direct_restoring_final_block_parser_binary_lookup_floor_mean={direct_restoring_final_block_parser_binary_lookup_floor_mean:.3}");
    println!("METRIC scratch600_direct_restoring_final_block_parser_binary_lookup_floor_p99={direct_restoring_final_block_parser_binary_lookup_floor_p99}");
    println!("METRIC scratch600_direct_restoring_final_block_parser_huffman_lookup_floor_mean={direct_restoring_final_block_parser_huffman_lookup_floor_mean:.3}");
    println!("METRIC scratch600_direct_restoring_final_block_parser_huffman_lookup_floor_p99={direct_restoring_final_block_parser_huffman_lookup_floor_p99}");
    println!("METRIC scratch600_direct_restoring_final_block_parser_best_with_binary_lookup_mean={direct_restoring_final_block_parser_best_with_binary_lookup_mean:.3}");
    println!("METRIC scratch600_direct_restoring_final_block_parser_best_with_binary_lookup_gap_to_2700k={direct_restoring_final_block_parser_best_with_binary_lookup_gap:.3}");
    println!("METRIC scratch600_direct_restoring_final_block_parser_best_with_binary_lookup_2x_mean={direct_restoring_final_block_parser_best_with_binary_lookup_2x_mean:.3}");
    println!("METRIC scratch600_direct_restoring_final_block_parser_best_with_binary_lookup_2x_gap_to_2700k={direct_restoring_final_block_parser_best_with_binary_lookup_2x_gap:.3}");
    println!("METRIC scratch600_direct_restoring_final_block4_with_binary_lookup_2x_mean={direct_restoring_final_block4_with_binary_lookup_2x_mean:.3}");
    println!("METRIC scratch600_direct_restoring_final_block4_with_binary_lookup_2x_gap_to_2700k={direct_restoring_final_block4_with_binary_lookup_2x_gap:.3}");
    println!("METRIC scratch600_direct_restoring_final_block4_lookup_multiplier_budget={direct_restoring_final_block4_lookup_multiplier_budget:.6}");
    println!("METRIC scratch600_direct_restoring_final_block5_with_binary_lookup_2x_mean={direct_restoring_final_block5_with_binary_lookup_2x_mean:.3}");
    println!("METRIC scratch600_direct_restoring_final_block5_with_binary_lookup_2x_gap_to_2700k={direct_restoring_final_block5_with_binary_lookup_2x_gap:.3}");
    println!("METRIC scratch600_direct_restoring_final_block5_lookup_multiplier_budget={direct_restoring_final_block5_lookup_multiplier_budget:.6}");
    println!("METRIC scratch600_direct_restoring_final_block7_with_binary_lookup_2x_mean={direct_restoring_final_block7_with_binary_lookup_2x_mean:.3}");
    println!("METRIC scratch600_direct_restoring_final_block7_with_binary_lookup_2x_gap_to_2700k={direct_restoring_final_block7_with_binary_lookup_2x_gap:.3}");
    println!("METRIC scratch600_direct_restoring_final_block7_lookup_multiplier_budget={direct_restoring_final_block7_lookup_multiplier_budget:.6}");
    println!("METRIC scratch600_direct_restoring_final_block_parser_cond_branch_best_block_symbols={direct_restoring_final_block_parser_cond_branch_best_block_symbols}");
    println!("METRIC scratch600_direct_restoring_final_block_parser_cond_branch_touch_floor_mean={direct_restoring_final_block_parser_cond_branch_touch_floor_mean:.3}");
    println!("METRIC scratch600_direct_restoring_final_block_parser_cond_branch_touch_floor_p99={direct_restoring_final_block_parser_cond_branch_touch_floor_p99}");
    println!("METRIC scratch600_direct_restoring_final_block_parser_cond_branch_compressed_bits_p99={direct_restoring_final_block_parser_cond_branch_compressed_bits_p99}");
    println!("METRIC scratch600_direct_restoring_final_block_parser_cond_branch_live_scratch_p99={direct_restoring_final_block_parser_cond_branch_live_scratch_p99}");
    println!("METRIC scratch600_direct_restoring_final_block_parser_cond_branch_augmented_gap_to_2700k={direct_restoring_final_block_parser_cond_branch_augmented_gap:.3}");
    println!("METRIC scratch600_direct_restoring_final_block_parser_cond_branch_binary_lookup_floor_mean={direct_restoring_final_block_parser_cond_branch_binary_lookup_floor_mean:.3}");
    println!("METRIC scratch600_direct_restoring_final_block_parser_cond_branch_binary_lookup_floor_p99={direct_restoring_final_block_parser_cond_branch_binary_lookup_floor_p99}");
    println!("METRIC scratch600_direct_restoring_final_block_parser_cond_branch_huffman_lookup_floor_mean={direct_restoring_final_block_parser_cond_branch_huffman_lookup_floor_mean:.3}");
    println!("METRIC scratch600_direct_restoring_final_block_parser_cond_branch_huffman_lookup_floor_p99={direct_restoring_final_block_parser_cond_branch_huffman_lookup_floor_p99}");
    println!("METRIC scratch600_direct_restoring_final_block_parser_cond_branch_best_with_binary_lookup_mean={direct_restoring_final_block_parser_cond_branch_best_with_binary_lookup_mean:.3}");
    println!("METRIC scratch600_direct_restoring_final_block_parser_cond_branch_best_with_binary_lookup_gap_to_2700k={direct_restoring_final_block_parser_cond_branch_best_with_binary_lookup_gap:.3}");
    println!("METRIC scratch600_direct_restoring_final_block_parser_cond_branch_best_with_binary_lookup_2x_mean={direct_restoring_final_block_parser_cond_branch_best_with_binary_lookup_2x_mean:.3}");
    println!("METRIC scratch600_direct_restoring_final_block_parser_cond_branch_best_with_binary_lookup_2x_gap_to_2700k={direct_restoring_final_block_parser_cond_branch_best_with_binary_lookup_2x_gap:.3}");
    println!("METRIC scratch600_direct_restoring_final_cond_block4_touch_floor_mean={direct_restoring_final_cond_block4_touch_floor_mean:.3}");
    println!("METRIC scratch600_direct_restoring_final_cond_block4_touch_floor_p99={direct_restoring_final_cond_block4_touch_floor_p99}");
    println!("METRIC scratch600_direct_restoring_final_cond_block4_compressed_bits_p99={direct_restoring_final_cond_block4_compressed_bits_p99}");
    println!("METRIC scratch600_direct_restoring_final_cond_block4_live_scratch_p99={direct_restoring_final_cond_block4_live_scratch_p99}");
    println!("METRIC scratch600_direct_restoring_final_cond_block4_symbol_count_p99={direct_restoring_final_cond_block4_symbol_count_p99}");
    println!("METRIC scratch600_direct_restoring_final_cond_block4_augmented_gap_to_2700k={direct_restoring_final_cond_block4_augmented_gap:.3}");
    println!("METRIC scratch600_direct_restoring_final_cond_block5_touch_floor_mean={direct_restoring_final_cond_block5_touch_floor_mean:.3}");
    println!("METRIC scratch600_direct_restoring_final_cond_block5_touch_floor_p99={direct_restoring_final_cond_block5_touch_floor_p99}");
    println!("METRIC scratch600_direct_restoring_final_cond_block5_compressed_bits_p99={direct_restoring_final_cond_block5_compressed_bits_p99}");
    println!("METRIC scratch600_direct_restoring_final_cond_block5_live_scratch_p99={direct_restoring_final_cond_block5_live_scratch_p99}");
    println!("METRIC scratch600_direct_restoring_final_cond_block5_symbol_count_p99={direct_restoring_final_cond_block5_symbol_count_p99}");
    println!("METRIC scratch600_direct_restoring_final_cond_block5_augmented_gap_to_2700k={direct_restoring_final_cond_block5_augmented_gap:.3}");
    println!("METRIC scratch600_direct_restoring_final_cond_block6_touch_floor_mean={direct_restoring_final_cond_block6_touch_floor_mean:.3}");
    println!("METRIC scratch600_direct_restoring_final_cond_block6_touch_floor_p99={direct_restoring_final_cond_block6_touch_floor_p99}");
    println!("METRIC scratch600_direct_restoring_final_cond_block6_compressed_bits_p99={direct_restoring_final_cond_block6_compressed_bits_p99}");
    println!("METRIC scratch600_direct_restoring_final_cond_block6_live_scratch_p99={direct_restoring_final_cond_block6_live_scratch_p99}");
    println!("METRIC scratch600_direct_restoring_final_cond_block6_symbol_count_p99={direct_restoring_final_cond_block6_symbol_count_p99}");
    println!("METRIC scratch600_direct_restoring_final_cond_block6_augmented_gap_to_2700k={direct_restoring_final_cond_block6_augmented_gap:.3}");
    println!("METRIC scratch600_direct_restoring_final_cond_block7_touch_floor_mean={direct_restoring_final_cond_block7_touch_floor_mean:.3}");
    println!("METRIC scratch600_direct_restoring_final_cond_block7_touch_floor_p99={direct_restoring_final_cond_block7_touch_floor_p99}");
    println!("METRIC scratch600_direct_restoring_final_cond_block7_compressed_bits_p99={direct_restoring_final_cond_block7_compressed_bits_p99}");
    println!("METRIC scratch600_direct_restoring_final_cond_block7_live_scratch_p99={direct_restoring_final_cond_block7_live_scratch_p99}");
    println!("METRIC scratch600_direct_restoring_final_cond_block7_symbol_count_p99={direct_restoring_final_cond_block7_symbol_count_p99}");
    println!("METRIC scratch600_direct_restoring_final_cond_block7_augmented_gap_to_2700k={direct_restoring_final_cond_block7_augmented_gap:.3}");
    println!("METRIC scratch600_direct_restoring_final_cond_block4_with_binary_lookup_2x_mean={direct_restoring_final_cond_block4_with_binary_lookup_2x_mean:.3}");
    println!("METRIC scratch600_direct_restoring_final_cond_block4_with_binary_lookup_2x_gap_to_2700k={direct_restoring_final_cond_block4_with_binary_lookup_2x_gap:.3}");
    println!("METRIC scratch600_direct_restoring_final_cond_block4_lookup_multiplier_budget={direct_restoring_final_cond_block4_lookup_multiplier_budget:.6}");
    println!("METRIC scratch600_direct_restoring_final_cond_block5_with_binary_lookup_2x_mean={direct_restoring_final_cond_block5_with_binary_lookup_2x_mean:.3}");
    println!("METRIC scratch600_direct_restoring_final_cond_block5_with_binary_lookup_2x_gap_to_2700k={direct_restoring_final_cond_block5_with_binary_lookup_2x_gap:.3}");
    println!("METRIC scratch600_direct_restoring_final_cond_block5_lookup_multiplier_budget={direct_restoring_final_cond_block5_lookup_multiplier_budget:.6}");
    println!("METRIC scratch600_direct_restoring_final_cond_block6_with_binary_lookup_2x_mean={direct_restoring_final_cond_block6_with_binary_lookup_2x_mean:.3}");
    println!("METRIC scratch600_direct_restoring_final_cond_block6_with_binary_lookup_2x_gap_to_2700k={direct_restoring_final_cond_block6_with_binary_lookup_2x_gap:.3}");
    println!("METRIC scratch600_direct_restoring_final_cond_block6_lookup_multiplier_budget={direct_restoring_final_cond_block6_lookup_multiplier_budget:.6}");
    println!("METRIC scratch600_direct_restoring_final_cond_block7_with_binary_lookup_2x_mean={direct_restoring_final_cond_block7_with_binary_lookup_2x_mean:.3}");
    println!("METRIC scratch600_direct_restoring_final_cond_block7_with_binary_lookup_2x_gap_to_2700k={direct_restoring_final_cond_block7_with_binary_lookup_2x_gap:.3}");
    println!("METRIC scratch600_direct_restoring_final_cond_block7_lookup_multiplier_budget={direct_restoring_final_cond_block7_lookup_multiplier_budget:.6}");
    println!("METRIC scratch600_direct_restoring_final_cond_mixed67_best_period={direct_restoring_final_cond_mixed67_best_period}");
    println!("METRIC scratch600_direct_restoring_final_cond_mixed67_best_mask={direct_restoring_final_cond_mixed67_best_mask}");
    println!("METRIC scratch600_direct_restoring_final_cond_mixed67_best_seven_count={direct_restoring_final_cond_mixed67_best_seven_count}");
    println!("METRIC scratch600_direct_restoring_final_cond_mixed67_touch_floor_mean={direct_restoring_final_cond_mixed67_touch_floor_mean:.3}");
    println!("METRIC scratch600_direct_restoring_final_cond_mixed67_touch_floor_p99={direct_restoring_final_cond_mixed67_touch_floor_p99}");
    println!("METRIC scratch600_direct_restoring_final_cond_mixed67_compressed_bits_p99={direct_restoring_final_cond_mixed67_compressed_bits_p99}");
    println!("METRIC scratch600_direct_restoring_final_cond_mixed67_live_scratch_p99={direct_restoring_final_cond_mixed67_live_scratch_p99}");
    println!("METRIC scratch600_direct_restoring_final_cond_mixed67_symbol_count_p99={direct_restoring_final_cond_mixed67_symbol_count_p99}");
    println!("METRIC scratch600_direct_restoring_final_cond_mixed67_augmented_gap_to_2700k={direct_restoring_final_cond_mixed67_augmented_gap:.3}");
    println!("METRIC scratch600_direct_restoring_final_cond_mixed67_with_binary_lookup_2x_mean={direct_restoring_final_cond_mixed67_with_binary_lookup_2x_mean:.3}");
    println!("METRIC scratch600_direct_restoring_final_cond_mixed67_with_binary_lookup_2x_gap_to_2700k={direct_restoring_final_cond_mixed67_with_binary_lookup_2x_gap:.3}");
    println!("METRIC scratch600_direct_restoring_final_cond_mixed67_lookup_multiplier_budget={direct_restoring_final_cond_mixed67_lookup_multiplier_budget:.6}");
    println!("METRIC scratch600_direct_restoring_final_cond_mixed4to8_best_period={direct_restoring_final_cond_mixed4to8_best_period}");
    println!("METRIC scratch600_direct_restoring_final_cond_mixed4to8_schedule_code={direct_restoring_final_cond_mixed4to8_schedule_code}");
    println!("METRIC scratch600_direct_restoring_final_cond_mixed4to8_touch_floor_mean={direct_restoring_final_cond_mixed4to8_touch_floor_mean:.3}");
    println!("METRIC scratch600_direct_restoring_final_cond_mixed4to8_touch_floor_p99={direct_restoring_final_cond_mixed4to8_touch_floor_p99}");
    println!("METRIC scratch600_direct_restoring_final_cond_mixed4to8_compressed_bits_p99={direct_restoring_final_cond_mixed4to8_compressed_bits_p99}");
    println!("METRIC scratch600_direct_restoring_final_cond_mixed4to8_live_scratch_p99={direct_restoring_final_cond_mixed4to8_live_scratch_p99}");
    println!("METRIC scratch600_direct_restoring_final_cond_mixed4to8_symbol_count_p99={direct_restoring_final_cond_mixed4to8_symbol_count_p99}");
    println!("METRIC scratch600_direct_restoring_final_cond_mixed4to8_augmented_gap_to_2700k={direct_restoring_final_cond_mixed4to8_augmented_gap:.3}");
    println!("METRIC scratch600_direct_restoring_final_cond_mixed4to8_with_binary_lookup_2x_mean={direct_restoring_final_cond_mixed4to8_with_binary_lookup_2x_mean:.3}");
    println!("METRIC scratch600_direct_restoring_final_cond_mixed4to8_with_binary_lookup_2x_gap_to_2700k={direct_restoring_final_cond_mixed4to8_with_binary_lookup_2x_gap:.3}");
    println!("METRIC scratch600_direct_restoring_final_cond_mixed4to8_lookup_multiplier_budget={direct_restoring_final_cond_mixed4to8_lookup_multiplier_budget:.6}");
    println!("METRIC scratch600_direct_restoring_final_cond_mixed4to8_block_joint_binary_lookup_mean={direct_restoring_final_cond_mixed4to8_block_joint_binary_lookup_mean:.3}");
    println!("METRIC scratch600_direct_restoring_final_cond_mixed4to8_block_joint_binary_lookup_p99={direct_restoring_final_cond_mixed4to8_block_joint_binary_lookup_p99}");
    println!("METRIC scratch600_direct_restoring_final_cond_mixed4to8_block_joint_support_row_floor={direct_restoring_final_cond_mixed4to8_block_joint_support_row_floor}");
    println!("METRIC scratch600_direct_restoring_final_cond_mixed4to8_block_joint_max_patterns={direct_restoring_final_cond_mixed4to8_block_joint_max_patterns}");
    println!("METRIC scratch600_direct_restoring_final_cond_mixed4to8_block_joint_block_count_p99={direct_restoring_final_cond_mixed4to8_block_joint_block_count_p99}");
    println!("METRIC scratch600_direct_restoring_final_cond_mixed4to8_with_block_joint_binary_lookup_2x_mean={direct_restoring_final_cond_mixed4to8_with_block_joint_binary_lookup_2x_mean:.3}");
    println!("METRIC scratch600_direct_restoring_final_cond_mixed4to8_with_block_joint_binary_lookup_2x_gap_to_2700k={direct_restoring_final_cond_mixed4to8_with_block_joint_binary_lookup_2x_gap:.3}");
    println!("METRIC scratch600_direct_restoring_final_cond_mixed4to8_block_joint_lookup_multiplier_budget={direct_restoring_final_cond_mixed4to8_block_joint_lookup_multiplier_budget:.6}");
    println!("METRIC scratch600_direct_restoring_final_cond_mixed4to8_with_block_joint_scan_lookup_2x_mean={direct_restoring_final_cond_mixed4to8_with_block_joint_scan_lookup_2x_mean:.3}");
    println!("METRIC scratch600_direct_restoring_final_cond_mixed4to8_with_block_joint_scan_lookup_2x_gap_to_2700k={direct_restoring_final_cond_mixed4to8_with_block_joint_scan_lookup_2x_gap:.3}");
    println!("METRIC scratch600_direct_restoring_final_selective_pair_lookup_baseline_mean={direct_restoring_final_selective_pair_lookup_baseline_mean:.3}");
    println!("METRIC scratch600_direct_restoring_final_selective_pair_lookup_selected_saving_mean={direct_restoring_final_selective_pair_lookup_selected_saving_mean:.3}");
    println!("METRIC scratch600_direct_restoring_final_selective_pair_lookup_required_saving_mean={direct_restoring_final_selective_pair_lookup_required_saving_mean:.3}");
    println!("METRIC scratch600_direct_restoring_final_selective_pair_lookup_mean={direct_restoring_final_selective_pair_lookup_mean:.3}");
    println!("METRIC scratch600_direct_restoring_final_selective_pair_lookup_target_mean={direct_restoring_final_selective_pair_lookup_target_mean:.3}");
    println!("METRIC scratch600_direct_restoring_final_selective_pair_lookup_gap_to_2700k={direct_restoring_final_selective_pair_lookup_gap:.3}");
    println!("METRIC scratch600_direct_restoring_final_selective_pair_lookup_selected_positions={direct_restoring_final_selective_pair_lookup_selected_positions}");
    println!("METRIC scratch600_direct_restoring_final_selective_pair_lookup_support_rows={direct_restoring_final_selective_pair_lookup_support_rows}");
    println!("METRIC scratch600_direct_restoring_final_selective_pair_lookup_max_patterns={direct_restoring_final_selective_pair_lookup_max_patterns}");
    println!("METRIC scratch600_direct_restoring_final_block_joint_rank_degree_n14={direct_restoring_final_block_joint_rank_degree_n14}");
    println!("METRIC scratch600_direct_restoring_final_block_joint_rank_density_n14={direct_restoring_final_block_joint_rank_density_n14}");
    println!("METRIC scratch600_direct_restoring_final_block_joint_rank_max_rank_n14={direct_restoring_final_block_joint_rank_max_rank_n14}");
    println!("METRIC scratch600_direct_restoring_final_block_joint_rank_max_patterns_n14={direct_restoring_final_block_joint_rank_max_patterns_n14}");
    println!("METRIC scratch600_direct_restoring_final_block_joint_rank_support_rows_n14={direct_restoring_final_block_joint_rank_support_rows_n14}");
    println!("METRIC scratch600_direct_restoring_final_block_joint_rank_max_blocks_n14={direct_restoring_final_block_joint_rank_max_blocks_n14}");
    println!("METRIC scratch600_direct_restoring_final_cond_mixed67_with_cond_scan_lookup_2x_mean={direct_restoring_final_cond_mixed67_with_cond_scan_lookup_2x_mean:.3}");
    println!("METRIC scratch600_direct_restoring_final_cond_mixed67_with_cond_scan_lookup_2x_gap_to_2700k={direct_restoring_final_cond_mixed67_with_cond_scan_lookup_2x_gap:.3}");
    println!("METRIC scratch600_direct_restoring_final_cond_mixed67_with_huffman_lookup_2x_mean={direct_restoring_final_cond_mixed67_with_huffman_lookup_2x_mean:.3}");
    println!("METRIC scratch600_direct_restoring_final_cond_mixed67_with_huffman_lookup_2x_gap_to_2700k={direct_restoring_final_cond_mixed67_with_huffman_lookup_2x_gap:.3}");
    println!("METRIC scratch600_direct_restoring_final_cond_mixed67_huffman_lookup_multiplier_budget={direct_restoring_final_cond_mixed67_huffman_lookup_multiplier_budget:.6}");
    println!("METRIC scratch600_direct_restoring_final_huffman_tree_toy_compare_ccx={direct_restoring_final_huffman_tree_toy_compare_ccx}");
    println!("METRIC scratch600_direct_restoring_final_huffman_tree_toy_forward_ccx={direct_restoring_final_huffman_tree_toy_forward_ccx}");
    println!("METRIC scratch600_direct_restoring_final_huffman_tree_toy_roundtrip_ccx={direct_restoring_final_huffman_tree_toy_roundtrip_ccx}");
    println!("METRIC scratch600_direct_restoring_final_huffman_tree_toy_peak_q={direct_restoring_final_huffman_tree_toy_peak_q}");
    println!("METRIC scratch600_direct_restoring_final_huffman_tree_toy_weighted_path_depth={direct_restoring_final_huffman_tree_toy_weighted_path_depth:.6}");
    println!("METRIC scratch600_direct_restoring_final_huffman_tree_toy_full_tree_nodes={direct_restoring_final_huffman_tree_toy_full_tree_nodes}");
    println!("METRIC scratch600_direct_restoring_final_huffman_tree_toy_path_compare_ccx_mean={direct_restoring_final_huffman_tree_toy_path_compare_ccx_mean:.3}");
    println!("METRIC scratch600_direct_restoring_final_huffman_tree_toy_full_over_path_ratio={direct_restoring_final_huffman_tree_toy_full_over_path_ratio:.6}");
    println!("METRIC scratch600_direct_restoring_final_huffman_tree_toy_dirty_restore_cases={direct_restoring_final_huffman_tree_toy_dirty_restore_cases}");
    println!("METRIC scratch600_direct_restoring_final_huffman_tree_toy_dirty_history_cases={direct_restoring_final_huffman_tree_toy_dirty_history_cases}");
    println!("METRIC scratch600_direct_restoring_final_huffman_tree_toy_dirty_phase_cases={direct_restoring_final_huffman_tree_toy_dirty_phase_cases}");
    println!("METRIC scratch600_direct_restoring_final_block_parser_align_support_noncontig_steps={direct_restoring_final_block_parser_align_support_noncontig_steps}");
    println!("METRIC scratch600_direct_restoring_final_block_parser_align_support_offset_steps={direct_restoring_final_block_parser_align_support_offset_steps}");
    println!("METRIC scratch600_direct_restoring_final_block_parser_align_support_max_span={direct_restoring_final_block_parser_align_support_max_span}");
    println!("METRIC scratch600_plusminus_raw_scratch_bits={plusminus_raw_scratch}");
    println!("METRIC scratch600_plusminus_unary_scratch_p99={plusminus_unary_scratch_p99}");
    println!("METRIC scratch600_plusminus_parser_over_strict_bits={plusminus_parser_over_strict}");
    println!("METRIC scratch600_plusminus_scaled_slack_scratch_max={plusminus_scaled_slack_scratch_max}");
    println!("METRIC scratch600_plusminus_scaled_solinas_projected_max={plusminus_scaled_solinas_projected_max}");
    println!("METRIC scratch600_plusminus_scaled_solinas_gap_max_to_2700k={plusminus_scaled_solinas_gap_max}");
    println!("METRIC scratch600_plusminus_solinas_scale_chunk_no_threshold_ccx={plusminus_solinas_scale_chunk_no_threshold_ccx}");
    println!("METRIC scratch600_plusminus_solinas_scale_chunk_no_threshold_peak_q={plusminus_solinas_scale_chunk_no_threshold_peak}");
    println!("METRIC scratch600_plusminus_solinas_scale_chunk_exact_ccx={plusminus_solinas_scale_chunk_exact_ccx}");
    println!("METRIC scratch600_plusminus_solinas_scale_chunk_exact_peak_q={plusminus_solinas_scale_chunk_exact_peak}");
    println!("METRIC scratch600_plusminus_solinas_scale_chunk_primitive_extra_q={plusminus_solinas_scale_chunk_primitive_extra}");
    println!("METRIC scratch600_plusminus_solinas_scale_chunk_naive_overlap_scratch={plusminus_solinas_scale_chunk_naive_overlap_scratch}");
    println!("METRIC scratch600_plusminus_solinas_scale_chunk_naive_over_google_bits={plusminus_solinas_scale_chunk_naive_over_google}");
    println!("METRIC scratch600_plusminus_solinas_scale_chunk_one_lane_reuse_scratch={plusminus_solinas_scale_chunk_one_lane_reuse_scratch}");
    println!("METRIC scratch600_plusminus_solinas_scale_chunk_one_lane_reuse_over_google_bits={plusminus_solinas_scale_chunk_one_lane_reuse_over_google}");
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
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_scratch_p99={halfgcd_second_col_fixed_depth64_scratch_p99}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_prefix_extract_width_sum_p99={halfgcd_second_col_fixed_depth64_prefix_extract_width_sum_p99}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_prefix_max_digits_p99={halfgcd_second_col_fixed_depth64_prefix_max_digits_p99}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_prefix_bounded_barrel_bits={halfgcd_second_col_fixed_depth64_prefix_bounded_barrel_bits}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_decoder_width_sum_p99={halfgcd_second_col_fixed_depth64_decoder_width_sum_p99}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_decoder_max_digits_p99={halfgcd_second_col_fixed_depth64_decoder_max_digits_p99}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_decoder_bounded_barrel_bits={halfgcd_second_col_fixed_depth64_decoder_bounded_barrel_bits}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_prefix_adversarial_prefix_max_digits={halfgcd_second_col_fixed_depth64_prefix_adversarial_prefix_max_digits}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_prefix_adversarial_decoder_max_digits={halfgcd_second_col_fixed_depth64_prefix_adversarial_decoder_max_digits}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_prefix_adversarial_required_barrel_bits={halfgcd_second_col_fixed_depth64_prefix_adversarial_required_barrel_bits}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_prefix_adversarial_missing_layers={halfgcd_second_col_fixed_depth64_prefix_adversarial_missing_layers}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_prefix_adversarial_prefix_width_sum={halfgcd_second_col_fixed_depth64_prefix_adversarial_prefix_width_sum}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_prefix_adversarial_decoder_width_sum={halfgcd_second_col_fixed_depth64_prefix_adversarial_decoder_width_sum}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_prefix_full_domain_avg_gap_floor={halfgcd_second_col_fixed_depth64_prefix_full_domain_avg_gap_floor}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_tail_bits_p99={halfgcd_second_col_fixed_depth64_tail_bits_p99}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_tail_count_p99={halfgcd_second_col_fixed_depth64_tail_count_p99}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_tail_width_sum_p99={halfgcd_second_col_fixed_depth64_tail_width_sum_p99}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_tail_max_q_bits_p99={halfgcd_second_col_fixed_depth64_tail_max_q_bits_p99}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_tail_bounded_barrel_bits={halfgcd_second_col_fixed_depth64_tail_bounded_barrel_bits}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_tail_adversarial_q_bits={halfgcd_second_col_fixed_depth64_tail_adversarial_q_bits}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_tail_adversarial_required_barrel_bits={halfgcd_second_col_fixed_depth64_tail_adversarial_required_barrel_bits}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_tail_adversarial_missing_layers={halfgcd_second_col_fixed_depth64_tail_adversarial_missing_layers}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_tail_adversarial_width_sum={halfgcd_second_col_fixed_depth64_tail_adversarial_width_sum}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_tail_adversarial_count={halfgcd_second_col_fixed_depth64_tail_adversarial_count}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_tail_full_domain_avg_gap_floor={halfgcd_second_col_fixed_depth64_tail_full_domain_avg_gap_floor}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_tail_extract_floor_p99={halfgcd_second_col_fixed_depth64_tail_extract_floor_p99}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_tail_bounded_barrel_floor_p99={halfgcd_second_col_fixed_depth64_tail_bounded_barrel_floor_p99}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_tail_logbarrel_floor_p99={halfgcd_second_col_fixed_depth64_tail_logbarrel_floor_p99}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_exact_mean={halfgcd_second_col_fixed_depth64_exact_mean}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_exact_p99={halfgcd_second_col_fixed_depth64_exact_p99}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_exact_tail_floor_mean={halfgcd_second_col_fixed_depth64_exact_tail_floor_mean}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_exact_tail_floor_p99={halfgcd_second_col_fixed_depth64_exact_tail_floor_p99}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_exact_tail_floor_gap_to_2700k={halfgcd_second_col_fixed_depth64_exact_tail_floor_gap}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_exact_tail_bounded_barrel_mean={halfgcd_second_col_fixed_depth64_exact_tail_bounded_barrel_mean}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_exact_tail_bounded_barrel_p99={halfgcd_second_col_fixed_depth64_exact_tail_bounded_barrel_p99}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_exact_tail_bounded_barrel_gap_to_2700k={halfgcd_second_col_fixed_depth64_exact_tail_bounded_barrel_gap}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_exact_tail_bounded_plus_one_width_mean={halfgcd_second_col_fixed_depth64_exact_tail_bounded_plus_one_width_mean}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_exact_tail_bounded_plus_one_width_p99={halfgcd_second_col_fixed_depth64_exact_tail_bounded_plus_one_width_p99}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_exact_tail_bounded_plus_one_width_gap_to_2700k={halfgcd_second_col_fixed_depth64_exact_tail_bounded_plus_one_width_gap}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_exact_tail_bounded_plus_two_width_mean={halfgcd_second_col_fixed_depth64_exact_tail_bounded_plus_two_width_mean}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_exact_tail_bounded_plus_two_width_p99={halfgcd_second_col_fixed_depth64_exact_tail_bounded_plus_two_width_p99}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_exact_tail_bounded_plus_two_width_gap_to_2700k={halfgcd_second_col_fixed_depth64_exact_tail_bounded_plus_two_width_gap}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_noscan_tail_bounded_barrel_mean={halfgcd_second_col_fixed_depth64_noscan_tail_bounded_barrel_mean}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_noscan_tail_bounded_barrel_p99={halfgcd_second_col_fixed_depth64_noscan_tail_bounded_barrel_p99}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_exact_tail_logbarrel_mean={halfgcd_second_col_fixed_depth64_exact_tail_logbarrel_mean}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_exact_tail_logbarrel_p99={halfgcd_second_col_fixed_depth64_exact_tail_logbarrel_p99}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_exact_tail_logbarrel_gap_to_2700k={halfgcd_second_col_fixed_depth64_exact_tail_logbarrel_gap}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_noscan_tail_logbarrel_mean={halfgcd_second_col_fixed_depth64_noscan_tail_logbarrel_mean}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_noscan_tail_logbarrel_p99={halfgcd_second_col_fixed_depth64_noscan_tail_logbarrel_p99}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_exact_prefix_bounded_tail_logbarrel_mean={halfgcd_second_col_fixed_depth64_exact_prefix_bounded_tail_logbarrel_mean}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_exact_prefix_bounded_tail_logbarrel_p99={halfgcd_second_col_fixed_depth64_exact_prefix_bounded_tail_logbarrel_p99}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_exact_prefix_bounded_tail_logbarrel_gap_to_2700k={halfgcd_second_col_fixed_depth64_exact_prefix_bounded_tail_logbarrel_gap}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_exact_prefix_bounded_tail_logbarrel_plus_one_width_mean={halfgcd_second_col_fixed_depth64_exact_prefix_bounded_tail_logbarrel_plus_one_width_mean}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_exact_prefix_bounded_tail_logbarrel_plus_one_width_p99={halfgcd_second_col_fixed_depth64_exact_prefix_bounded_tail_logbarrel_plus_one_width_p99}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_exact_prefix_bounded_tail_logbarrel_plus_one_width_gap_to_2700k={halfgcd_second_col_fixed_depth64_exact_prefix_bounded_tail_logbarrel_plus_one_width_gap}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_exact_prefix_bounded_tail_logbarrel_plus_two_width_mean={halfgcd_second_col_fixed_depth64_exact_prefix_bounded_tail_logbarrel_plus_two_width_mean}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_exact_prefix_bounded_tail_logbarrel_plus_two_width_p99={halfgcd_second_col_fixed_depth64_exact_prefix_bounded_tail_logbarrel_plus_two_width_p99}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_exact_prefix_bounded_tail_logbarrel_plus_two_width_gap_to_2700k={halfgcd_second_col_fixed_depth64_exact_prefix_bounded_tail_logbarrel_plus_two_width_gap}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_dynamic_barrel_static_mean={halfgcd_second_col_fixed_depth64_dynamic_barrel_static_mean}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_dynamic_barrel_static_p99={halfgcd_second_col_fixed_depth64_dynamic_barrel_static_p99}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_dynamic_barrel_mean={halfgcd_second_col_fixed_depth64_dynamic_barrel_mean}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_dynamic_barrel_p99={halfgcd_second_col_fixed_depth64_dynamic_barrel_p99}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_dynamic_barrel_gap_to_2700k={halfgcd_second_col_fixed_depth64_dynamic_barrel_gap}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_dynamic_barrel_savings_mean={halfgcd_second_col_fixed_depth64_dynamic_barrel_savings_mean}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_dynamic_barrel_scratch_p99={halfgcd_second_col_fixed_depth64_dynamic_barrel_scratch_p99}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_dynamic_prefix_decoder_static_mean={halfgcd_second_col_fixed_depth64_dynamic_prefix_decoder_static_mean}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_dynamic_prefix_decoder_mean={halfgcd_second_col_fixed_depth64_dynamic_prefix_decoder_mean}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_dynamic_tail_static_mean={halfgcd_second_col_fixed_depth64_dynamic_tail_static_mean}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_dynamic_tail_mean={halfgcd_second_col_fixed_depth64_dynamic_tail_mean}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_dynamic_high_layer_hits_p99={halfgcd_second_col_fixed_depth64_dynamic_high_layer_hits_p99}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_static_app_mean={halfgcd_second_col_fixed_depth64_static_app_mean}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_static_app_p99={halfgcd_second_col_fixed_depth64_static_app_p99}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_static_app_gap_to_2700k={halfgcd_second_col_fixed_depth64_static_app_gap}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_app_popcount_mean={halfgcd_second_col_fixed_depth64_app_popcount_mean}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_app_static_floor_mean={halfgcd_second_col_fixed_depth64_app_static_floor_mean}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_app_static_over_popcount_mean={halfgcd_second_col_fixed_depth64_app_static_over_popcount_mean}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_static_sep4_app_mean={halfgcd_second_col_fixed_depth64_static_sep4_app_mean}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_static_sep4_app_p99={halfgcd_second_col_fixed_depth64_static_sep4_app_p99}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_static_sep4_app_gap_to_2700k={halfgcd_second_col_fixed_depth64_static_sep4_app_gap}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_static_joint4_app_mean={halfgcd_second_col_fixed_depth64_static_joint4_app_mean}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_static_joint4_app_p99={halfgcd_second_col_fixed_depth64_static_joint4_app_p99}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_static_joint4_app_gap_to_2700k={halfgcd_second_col_fixed_depth64_static_joint4_app_gap}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_static_sep4_selector_budget_oneway={halfgcd_second_col_fixed_depth64_static_sep4_selector_budget_oneway}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_static_joint4_selector_budget_oneway={halfgcd_second_col_fixed_depth64_static_joint4_selector_budget_oneway}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_app_static_sep4_floor_mean={halfgcd_second_col_fixed_depth64_app_static_sep4_floor_mean}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_app_static_joint4_floor_mean={halfgcd_second_col_fixed_depth64_app_static_joint4_floor_mean}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_static_sep4_with_selector_floor_mean={halfgcd_second_col_fixed_depth64_static_sep4_with_selector_floor_mean}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_static_sep4_with_selector_floor_p99={halfgcd_second_col_fixed_depth64_static_sep4_with_selector_floor_p99}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_static_joint4_with_selector_floor_mean={halfgcd_second_col_fixed_depth64_static_joint4_with_selector_floor_mean}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_static_joint4_with_selector_floor_p99={halfgcd_second_col_fixed_depth64_static_joint4_with_selector_floor_p99}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_static_joint4_with_selector_floor_gap_to_2700k={halfgcd_second_col_fixed_depth64_static_joint4_with_selector_floor_gap}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_app_static_selector_floor_mean={halfgcd_second_col_fixed_depth64_app_static_selector_floor_mean}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_app_static_selector_floor_over_joint4_budget={halfgcd_second_col_fixed_depth64_app_static_selector_floor_over_joint4_budget}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_static_window_scan_best_w={halfgcd_second_col_fixed_depth64_static_window_scan_best_w}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_static_window_scan_best_mean={halfgcd_second_col_fixed_depth64_static_window_scan_best_mean}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_static_window_scan_best_p99={halfgcd_second_col_fixed_depth64_static_window_scan_best_p99}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_static_window_scan_best_gap_to_2700k={halfgcd_second_col_fixed_depth64_static_window_scan_best_gap}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_static_window_scan_best_app_mean={halfgcd_second_col_fixed_depth64_static_window_scan_best_app_mean}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_static_window_scan_best_selector_mean={halfgcd_second_col_fixed_depth64_static_window_scan_best_selector_mean}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_static_window_scan_best_table_row_mean={halfgcd_second_col_fixed_depth64_static_window_scan_best_table_row_mean}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_static_window_table_only_best_w={halfgcd_second_col_fixed_depth64_static_window_table_only_best_w}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_static_window_table_only_best_mean={halfgcd_second_col_fixed_depth64_static_window_table_only_best_mean}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_static_window_table_only_best_p99={halfgcd_second_col_fixed_depth64_static_window_table_only_best_p99}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_static_window_table_only_best_gap_to_2700k={halfgcd_second_col_fixed_depth64_static_window_table_only_best_gap}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_static_window_table_source_best_w={halfgcd_second_col_fixed_depth64_static_window_table_source_best_w}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_static_window_table_source_mean={halfgcd_second_col_fixed_depth64_static_window_table_source_mean}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_static_window_table_source_p99={halfgcd_second_col_fixed_depth64_static_window_table_source_p99}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_static_window_table_source_gap_to_2700k={halfgcd_second_col_fixed_depth64_static_window_table_source_gap}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_static_window_table_source_product_floor_mean={halfgcd_second_col_fixed_depth64_static_window_table_source_product_floor_mean}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_static_window_required_selector_mean={halfgcd_second_col_fixed_depth64_static_window_required_selector_mean}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_static_window_selector_cut_needed={halfgcd_second_col_fixed_depth64_static_window_selector_cut_needed}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_static_window_table_margin={halfgcd_second_col_fixed_depth64_static_window_table_margin}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_static_window_source_product_best_w={halfgcd_second_col_fixed_depth64_static_window_source_product_best_w}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_static_window_source_product_mean={halfgcd_second_col_fixed_depth64_static_window_source_product_mean}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_static_window_source_product_p99={halfgcd_second_col_fixed_depth64_static_window_source_product_p99}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_static_window_source_product_gap_to_2700k={halfgcd_second_col_fixed_depth64_static_window_source_product_gap}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_static_window_source_product_floor_mean={halfgcd_second_col_fixed_depth64_static_window_source_product_floor_mean}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_static_window_source_product_table_row_mean={halfgcd_second_col_fixed_depth64_static_window_source_product_table_row_mean}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_static_window_wnaf_best_w={halfgcd_second_col_fixed_depth64_static_window_wnaf_best_w}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_static_window_wnaf_mean={halfgcd_second_col_fixed_depth64_static_window_wnaf_mean}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_static_window_wnaf_p99={halfgcd_second_col_fixed_depth64_static_window_wnaf_p99}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_static_window_wnaf_gap_to_2700k={halfgcd_second_col_fixed_depth64_static_window_wnaf_gap}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_static_window_wnaf_app_mean={halfgcd_second_col_fixed_depth64_static_window_wnaf_app_mean}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_static_window_wnaf_selector_floor_mean={halfgcd_second_col_fixed_depth64_static_window_wnaf_selector_floor_mean}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_static_window_wnaf_source_product_floor_mean={halfgcd_second_col_fixed_depth64_static_window_wnaf_source_product_floor_mean}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_static_window_wnaf_table_row_floor_mean={halfgcd_second_col_fixed_depth64_static_window_wnaf_table_row_floor_mean}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_static_window_wnaf_positions_mean={halfgcd_second_col_fixed_depth64_static_window_wnaf_positions_mean:.3}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_static_window_wnaf_compact_best_w={halfgcd_second_col_fixed_depth64_static_window_wnaf_compact_best_w}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_static_window_wnaf_compact_mean={halfgcd_second_col_fixed_depth64_static_window_wnaf_compact_mean}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_static_window_wnaf_compact_p99={halfgcd_second_col_fixed_depth64_static_window_wnaf_compact_p99}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_static_window_wnaf_compact_gap_to_2700k={halfgcd_second_col_fixed_depth64_static_window_wnaf_compact_gap}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_static_window_wnaf_compact_app_mean={halfgcd_second_col_fixed_depth64_static_window_wnaf_compact_app_mean}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_static_window_wnaf_compact_source_product_floor_mean={halfgcd_second_col_fixed_depth64_static_window_wnaf_compact_source_product_floor_mean}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_static_window_wnaf_compact_missing_active_floor_mean={halfgcd_second_col_fixed_depth64_static_window_wnaf_compact_missing_active_floor_mean}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_static_window_wnaf_compact_active_slack_oneway={halfgcd_second_col_fixed_depth64_static_window_wnaf_compact_active_slack_oneway}");
    println!("METRIC scratch600_halfgcd_second_col_fixed_depth64_static_window_wnaf_compact_table_row_floor_mean={halfgcd_second_col_fixed_depth64_static_window_wnaf_compact_table_row_floor_mean}");
    println!("METRIC scratch600_halfgcd_second_col_compact_wnaf_active_degree_n14={halfgcd_second_col_compact_wnaf_active_degree_n14}");
    println!("METRIC scratch600_halfgcd_second_col_compact_wnaf_active_density_n14={halfgcd_second_col_compact_wnaf_active_density_n14}");
    println!("METRIC scratch600_halfgcd_second_col_compact_wnaf_active_positions_n14={halfgcd_second_col_compact_wnaf_active_positions_n14}");
    println!("METRIC scratch600_halfgcd_second_col_compact_wnaf_active_pair_positions_n14={halfgcd_second_col_compact_wnaf_active_pair_positions_n14}");
    println!("METRIC scratch600_halfgcd_second_col_compact_wnaf_active_slots_n14={halfgcd_second_col_compact_wnaf_active_slots_n14}");
    println!("METRIC scratch600_halfgcd_second_col_compact_wnaf_active_full_slots_n14={halfgcd_second_col_compact_wnaf_active_full_slots_n14}");
    println!("METRIC scratch600_halfgcd_second_col_compact_wnaf_active_max_pair_n14={halfgcd_second_col_compact_wnaf_active_max_pair_n14}");
    println!("METRIC scratch600_halfgcd_second_col_alignment_mbu_degree_n14={halfgcd_second_col_alignment_mbu_degree_n14}");
    println!("METRIC scratch600_halfgcd_second_col_alignment_mbu_density_n14={halfgcd_second_col_alignment_mbu_density_n14}");
    println!("METRIC scratch600_halfgcd_second_col_alignment_mbu_max_alignment_n14={halfgcd_second_col_alignment_mbu_max_alignment_n14}");
    println!("METRIC scratch600_halfgcd_second_col_static_window_mbu_degree_n14={halfgcd_second_col_static_window_mbu_degree_n14}");
    println!("METRIC scratch600_halfgcd_second_col_static_window_mbu_density_n14={halfgcd_second_col_static_window_mbu_density_n14}");
    println!("METRIC scratch600_halfgcd_second_col_static_window_mbu_max_coeff_bits_n14={halfgcd_second_col_static_window_mbu_max_coeff_bits_n14}");
    println!("METRIC scratch600_halfgcd_second_col_static_window_mbu_max_pair_n14={halfgcd_second_col_static_window_mbu_max_pair_n14}");
    println!("METRIC scratch600_halfgcd_second_col_static_window_wnaf_mbu_degree_n14={halfgcd_second_col_static_window_wnaf_mbu_degree_n14}");
    println!("METRIC scratch600_halfgcd_second_col_static_window_wnaf_mbu_density_n14={halfgcd_second_col_static_window_wnaf_mbu_density_n14}");
    println!("METRIC scratch600_halfgcd_second_col_static_window_wnaf_mbu_max_positions_n14={halfgcd_second_col_static_window_wnaf_mbu_max_positions_n14}");
    println!("METRIC scratch600_halfgcd_second_col_static_window_wnaf_mbu_max_pair_n14={halfgcd_second_col_static_window_wnaf_mbu_max_pair_n14}");
    println!("METRIC scratch600_halfgcd_second_col_static_window_support_rows_n14={halfgcd_second_col_static_window_support_rows_n14}");
    println!("METRIC scratch600_halfgcd_second_col_static_window_support_full_rows_n14={halfgcd_second_col_static_window_support_full_rows_n14}");
    println!("METRIC scratch600_halfgcd_second_col_static_window_support_ppm_n14={halfgcd_second_col_static_window_support_ppm_n14}");
    println!("METRIC scratch600_halfgcd_second_col_static_window_support_saturated_windows_n14={halfgcd_second_col_static_window_support_saturated_windows_n14}");
    println!("METRIC scratch600_halfgcd_second_col_static_window_support_windows_n14={halfgcd_second_col_static_window_support_windows_n14}");
    println!("METRIC scratch600_halfgcd_second_col_static_window_bit_support_n14={halfgcd_second_col_static_window_bit_support_n14}");
    println!("METRIC scratch600_halfgcd_second_col_static_window_full_bits_n14={halfgcd_second_col_static_window_full_bits_n14}");
    println!("METRIC scratch600_halfgcd_second_col_static_window_bit_support_ppm_n14={halfgcd_second_col_static_window_bit_support_ppm_n14}");

    assert!(best_state <= STRICT_SCRATCH, "at least some state shapes fit");
    assert!(streamed_gap_to_google > 0, "no fully charged <=600-scratch row should be counted as solved yet");
    assert!(streamed_selector_shortfall > 0, "streamed-mask route still needs a selector breakthrough");
    assert!(
        tiny_lowword_w1_selector_slack > 0 && tiny_lowword_best_fixed_update_excess > 250_000,
        "tiny lowword selector/update tradeoff changed; revisit streamed BY route"
    );
    assert!(
        scaled_by_pattern_fixed_id_bits < 560
            && scaled_by_pattern_fixed_id_distinct_rows
                > scaled_by_pattern_fixed_id_remaining_to_2700k
            && scaled_by_pattern_fixed_id_row_floor_gap > 200_000
            && scaled_by_pattern_fixed_id_bit_floor_gap > 2_000_000,
        "compressed scaled-BY pattern IDs now have a generic decoder budget; revisit fixed-ID history"
    );
    assert!(
        scaled_by_raw_pattern_single_a_scratch <= GOOGLE_LOW_QUBIT_SCRATCH
            && scaled_by_raw_pattern_one_checkpoint_scratch > GOOGLE_LOW_QUBIT_SCRATCH
            && scaled_by_raw_pattern_window_a_scratch > GOOGLE_LOW_QUBIT_SCRATCH
            && scaled_by_raw_pattern_delta_checkpoint_bits
                > scaled_by_raw_pattern_delta_checkpoint_scratch_slack
            && scaled_by_raw_pattern_delta_checkpoint_scratch > GOOGLE_LOW_QUBIT_SCRATCH
            && scaled_by_raw_pattern_ambiguous_a_bits_p99 > 200
            && scaled_by_raw_pattern_exact_two_decoder_gap > 0,
        "raw-pattern scaled-BY streaming now has reversible scratch/decode margin; revisit raw history"
    );
    assert!(
        scaled_by_h_only_model_modular_windows == 35
            && scaled_by_h_only_model_modular_toffoli < 900_000
            && scaled_by_h_only_model_peak < 2_800
            && scaled_by_h_only_model_history_bits == 480
            && scaled_by_h_only_next_ratio_toy_n14_max_next_h_choices == 16
            && scaled_by_h_only_next_ratio_toy_n14_rank_p99 >= 28
            && scaled_by_h_only_next_ratio_toy_n14_rank_mean_milli > 27_000
            && scaled_by_h_only_next_ratio_toy_n14_ambiguous_keys > 700,
        "h-only BY next-ratio payload became cheap; revisit compressed history update"
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
    assert!(
        by_centered_exactparity_clean_replay_ccx
            > by_centered_exactparity_clean_per_div_budget
            && by_centered_exactparity_clean_peak_q < 2_800
            && by_centered_exactparity_clean_scratch_bits > GOOGLE_LOW_QUBIT_SCRATCH
            && by_centered_exactparity_two_clean_div_gap > 3_000_000,
        "centered exact-parity clean BY replay changed; revisit whether exact parity cleanup can be promoted"
    );
    assert!(centered_parser_over_strict > 0 && plusminus_parser_over_strict > 0, "raw streams must not be counted before parser cost");
    assert!(
        plusminus_active_quantum_gap_to_2700k > 50_000_000,
        "plus-minus active-chain quantum-control blocker changed; revisit physical integration"
    );
    assert!(
        plusminus_scaled_solinas_projected_max < GOOGLE_LOW_QUBIT_TOFFOLI
            && plusminus_solinas_scale_chunk_no_threshold_ccx == 3_390
            && plusminus_solinas_scale_chunk_no_threshold_peak > GOOGLE_LOW_QUBIT_SCRATCH
            && plusminus_solinas_scale_chunk_one_lane_reuse_over_google > 0
            && plusminus_solinas_scale_chunk_naive_over_google > 400,
        "plus-minus Solinas scale chunk now fits packed scratch; revisit scale-history route"
    );
    assert!(
        direct_signnorm_rank_over_google > 0 && direct_signnorm_ambiguous_rank_over_google > 0,
        "archival rank-compressed normalization-sign sidecar changed; update signnorm ledger"
    );
    assert!(
        direct_signnorm_det_coeffsign_scratch_gap_google <= 0,
        "det-low2 coefficient-sign recovery no longer fits Google scratch"
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
        direct_signnorm_logsign_direct_rem_toy_ccx == 148
            && direct_signnorm_logsign_direct_rem_toy_peak_q == 80
            && direct_signnorm_logsign_direct_rem_toy_phase_dirty_cases == 0
            && direct_signnorm_logsign_exact_cneg257 == 512
            && direct_signnorm_logsign_exact_rem_p99 > 20_000
            && direct_signnorm_logsign_exact_once_gap == direct_signnorm_logsign_split_gap
            && direct_signnorm_logsign_exact_split_gap > direct_signnorm_logsign_split_gap,
        "logical-sign rem-cneg phase/cost changed; revisit direct sign-normalized route"
    );
    assert!(
        direct_signnorm_logsign_no_rem_cneg_gap < 0
            && direct_signnorm_prefinal_signed_remainder_gap > 400_000
            && direct_signnorm_prefinal_signed_remainder_count_p99 >= 180
            && direct_signnorm_prefinal_signed_remainder_digit_payload_p99 >= 498
            && direct_signnorm_prefinal_signed_remainder_width_extra_max == 1,
        "sign-normalized rem-cneg escape changed; revisit logical-remainder route"
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
        direct_signnorm_coeff_reverse_collisions_n14 == 0
            && direct_signnorm_coeff_reverse_states_n14 == direct_signnorm_reverse_total_steps_n14
            && direct_signnorm_coeff_reverse_total_steps_n14
                == direct_signnorm_reverse_total_steps_n14
            && direct_signnorm_coeff_reverse_max_mult_n14 == 1
            && direct_signnorm_coeff_reverse_zero_coeff_cases_n14 == 0,
        "coefficient rows stopped disambiguating sign-normalized reverse signs; require explicit sign history again"
    );
    assert!(
        direct_signnorm_det_sign_reverse_collisions_n14 > 2_000
            && direct_signnorm_det_sign_reverse_states_n14 > 70_000
            && direct_signnorm_det_sign_reverse_max_mult_n14 == 2
            && direct_signnorm_det_coeffsign_reverse_collisions_n14 == 0
            && direct_signnorm_det_coeffsign_reverse_states_n14 > 73_000
            && direct_signnorm_det_coeffsign_reverse_total_steps_n14
                == direct_signnorm_reverse_total_steps_n14
            && direct_signnorm_det_coeffsign_reverse_max_mult_n14 == 1
            && direct_signnorm_det_coeffsign_bad_det_cases_n14 == 0
            && direct_signnorm_det_coeffsign_low2_mismatches_n14 == 0
            && direct_signnorm_det_coeffsign_formula_mismatches_n14 == 0,
        "det-low2 xor coeff_v_sign stopped recovering sign-normalized norm signs"
    );
    assert!(
        direct_signnorm_det_coeffsign_predicate_p1_ccx == 14
            && direct_signnorm_det_coeffsign_predicate_p1_peak_q <= 18
            && direct_signnorm_det_coeffsign_predicate_p1_valid_odd_det_cases > 3_000
            && direct_signnorm_det_coeffsign_predicate_p3_ccx == 14
            && direct_signnorm_det_coeffsign_predicate_p3_peak_q <= 18
            && direct_signnorm_det_coeffsign_predicate_p3_valid_odd_det_cases > 3_000,
        "det-low2 coefficient-sign recovery predicate toy changed"
    );
    assert!(
        direct_signnorm_signed_domain_relative_negative_toy_ccx == 45
            && direct_signnorm_signed_domain_relative_negative_257_ccx
                > direct_signnorm_logsign_exact_cneg257
            && direct_signnorm_signed_domain_floor_toy_ccx == 416
            && direct_signnorm_signed_domain_floor_toy_peak_q == 62
            && direct_signnorm_signed_domain_floor_toy_final_negative_cases > 1_900,
        "signed-domain floor body needs an expensive zero-guarded relative sign predicate"
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
        direct_restoring_final_branch_digit_toy_branch_ccx
            < direct_restoring_final_single_selector_toy_ccx / 10
            && direct_restoring_final_branch_digit_toy_forward_ccx
                == direct_restoring_final_single_selector_toy_ccx
                    + direct_restoring_final_branch_digit_toy_branch_ccx
            && direct_restoring_final_branch_digit_toy_roundtrip_ccx
                == 2 * direct_restoring_final_branch_digit_toy_forward_ccx
            && direct_restoring_final_branch_digit_toy_peak_q
                == direct_restoring_final_single_selector_bennett_toy_peak_q + 1
            && direct_restoring_final_branch_digit_toy_branch_one_cases > 0,
        "restoring-final branch-digit toy changed; revisit low-branch fold integration"
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
        direct_restoring_final_low_branch_adjacent_transitions_n14
            == direct_restoring_final_reverse_coeff_candidates_transitions_n14
            && direct_restoring_final_low_branch_adjacent_ambiguous_n14
                == direct_restoring_final_reverse_coeff_candidates_low_n14
                    + direct_restoring_final_reverse_coeff_candidates_high_n14
            && direct_restoring_final_low_branch_adjacent_ambiguous_n14
                == direct_restoring_final_reverse_q_total_steps_n14
            && direct_restoring_final_low_branch_adjacent_high_n14
                == direct_restoring_final_reverse_coeff_candidates_high_n14
            && direct_restoring_final_low_branch_adjacent_violations_n14 == 0
            && direct_restoring_final_low_branch_adjacent_max_delta_n14 == 1
            && direct_restoring_final_low_branch_adjacent_max_alignment_n14 + 1 >= 14,
        "low-branch high candidate stopped being an exact final digit"
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
            && direct_restoring_final_stored_align_gamma_scratch_p99 > GOOGLE_LOW_QUBIT_SCRATCH
            && direct_restoring_final_stored_align_length_rank_scratch_p99
                > GOOGLE_LOW_QUBIT_SCRATCH
            && direct_restoring_final_stored_align_public_len_position_only_p99
                > GOOGLE_LOW_QUBIT_SCRATCH
            && direct_restoring_final_stored_align_public_len_position_plus3_p99
                > GOOGLE_LOW_QUBIT_SCRATCH
            && direct_restoring_final_stored_align_q_len_position_only_p99
                > GOOGLE_LOW_QUBIT_SCRATCH
            && direct_restoring_final_stored_align_digit_len_position_only_p99
                > GOOGLE_LOW_QUBIT_SCRATCH
            && direct_restoring_final_stored_align_joint_len_position_only_p99
                > GOOGLE_LOW_QUBIT_SCRATCH,
        "restoring-final stored-alignment metadata changed; revisit exact parser/packing blocker"
    );
    assert!(
        direct_restoring_final_branch_final_high_adjacent_violations == 0
            && direct_restoring_final_branch_final_current_branch_select_p99
                == direct_restoring_final_stored_align_branch_select_p99
            && direct_restoring_final_branch_final_branch_count_p99
                == direct_restoring_final_stored_align_branch_count_p99
            && direct_restoring_final_branch_final_selected_width_saving_mean > 13_000.0
            && direct_restoring_final_branch_final_low_path_width_saving_mean > 16_000.0
            && direct_restoring_final_branch_final_low_extra_touch_mean < 0.0
            && direct_restoring_final_branch_final_selected_width_mixed4to8_gap < 0.0
            && direct_restoring_final_branch_final_low_path_width_mixed4to8_gap < -50_000.0
            && direct_restoring_final_branch_final_current_mixed4to8_gap > 8_000.0
            && direct_restoring_final_branch_final_low_path_width_scan_mixed4to8_gap > 39_000.0
            && direct_restoring_final_branch_final_low_path_width_lookup_multiplier_budget > 2.0
            && direct_restoring_final_branch_final_low_path_width_lookup_multiplier_budget < 2.1
            && direct_restoring_final_branch_final_low_path_width_lookup_target_mean
                < direct_restoring_final_block_parser_cond_branch_lookup_scan_floor_mean
            && direct_restoring_final_branch_final_alignment_diff_p99 > 0
            && direct_restoring_final_branch_final_digit_len_diff_p99 > 0,
        "branch-as-final-digit restoring-final lower bound changed; revisit branch-fold circuit target"
    );
    assert!(
        direct_restoring_final_low_branch_align_only_model_precision_bits == 13
            && direct_restoring_final_low_branch_align_only_raw_scratch_p99
                <= GOOGLE_LOW_QUBIT_SCRATCH
            && direct_restoring_final_low_branch_align_only_step_prefix_scratch_p99
                < GOOGLE_LOW_QUBIT_SCRATCH
            && direct_restoring_final_low_branch_align_only_best_block_symbols == 2
            && direct_restoring_final_low_branch_align_only_best_live_scratch_p99
                < GOOGLE_LOW_QUBIT_SCRATCH
            && direct_restoring_final_low_branch_align_only_best_augmented_gap < -100_000.0
            && direct_restoring_final_low_branch_align_only_best_binary_gap < -70_000.0
            && direct_restoring_final_low_branch_align_only_best_prefix_tree_gap < -100_000.0
            && direct_restoring_final_low_branch_align_only_best_scan_gap > 30_000.0
            && direct_restoring_final_low_branch_align_only_best_lookup_multiplier_budget > 2.6
            && direct_restoring_final_low_branch_align_only_scan_over_binary_multiplier
                > direct_restoring_final_low_branch_align_only_best_lookup_multiplier_budget
            && direct_restoring_final_low_branch_align_only_huffman_over_binary_multiplier < 1.0
            && direct_restoring_final_low_branch_align_only_prefix_tree_over_binary_multiplier < 0.3
            && direct_restoring_final_low_branch_align_only_support_noncontig_steps > 50,
        "low-branch alignment-only parser budget changed; revisit prefix/lookup decoder target"
    );
    assert!(
        direct_restoring_final_prefix_bit_reader_toy_eq_ccx == 4
            && direct_restoring_final_prefix_bit_reader_toy_dynamic_read_ccx == 16
            && direct_restoring_final_prefix_bit_reader_toy_tree_forward_ccx == 6
            && direct_restoring_final_prefix_bit_reader_toy_roundtrip_ccx == 52
            && direct_restoring_final_prefix_bit_reader_toy_peak_q == 31
            && direct_restoring_final_prefix_bit_reader_toy_full_over_node_roundtrip
                < direct_restoring_final_prefix_bit_reader_toy_roundtrip_ratio_budget
            && direct_restoring_final_prefix_bit_reader_toy_cursor_scaled_gap < 0.0
            && direct_restoring_final_prefix_bit_reader_toy_dirty_restore_cases == 0
            && direct_restoring_final_prefix_bit_reader_toy_dirty_history_cases == 0
            && direct_restoring_final_prefix_bit_reader_toy_dirty_phase_cases == 0,
        "prefix bit-reader toy changed; revisit low-branch cursor/parser integration"
    );
    assert!(
        direct_restoring_final_prefix_cursor_advance_toy_ccx == 6
            && direct_restoring_final_prefix_cursor_advance_toy_peak_q == 11
            && direct_restoring_final_prefix_cursor_advance_toy_combined_roundtrip_ccx == 58
            && direct_restoring_final_prefix_cursor_advance_toy_combined_over_node_roundtrip
                < direct_restoring_final_prefix_cursor_advance_toy_roundtrip_ratio_budget
            && direct_restoring_final_prefix_cursor_advance_toy_combined_scaled_gap < -30_000.0
            && direct_restoring_final_prefix_cursor_advance_toy_dirty_restore_cases == 0
            && direct_restoring_final_prefix_cursor_advance_toy_dirty_history_cases == 0
            && direct_restoring_final_prefix_cursor_advance_toy_dirty_phase_cases == 0,
        "prefix cursor advance toy changed; revisit low-branch parser cursor budget"
    );
    assert!(
        direct_restoring_final_prefix_block2_toy_tree_ccx == 12
            && direct_restoring_final_prefix_block2_toy_read2_ccx == 16
            && direct_restoring_final_prefix_block2_toy_cursor_add_ccx == 10
            && direct_restoring_final_prefix_block2_toy_decode_forward_ccx == 28
            && direct_restoring_final_prefix_block2_toy_total_ccx == 66
            && direct_restoring_final_prefix_block2_toy_peak_q == 52
            && direct_restoring_final_prefix_block2_toy_over_node_roundtrip
                < direct_restoring_final_prefix_block2_toy_roundtrip_ratio_budget
            && direct_restoring_final_prefix_block2_toy_scaled_gap < -60_000.0
            && direct_restoring_final_prefix_block2_toy_dirty_restore_cases == 0
            && direct_restoring_final_prefix_block2_toy_dirty_history_cases == 0
            && direct_restoring_final_prefix_block2_toy_dirty_phase_cases == 0,
        "prefix block2 toy changed; revisit low-branch parser block budget"
    );
    assert!(
        direct_restoring_final_prefix_block2_consume_toy_decode_forward_ccx == 28
            && direct_restoring_final_prefix_block2_consume_toy_cursor_add_ccx == 10
            && direct_restoring_final_prefix_block2_consume_toy_consume_ccx == 4
            && direct_restoring_final_prefix_block2_consume_toy_parser_transient_ccx == 76
            && direct_restoring_final_prefix_block2_consume_toy_total_ccx == 80
            && direct_restoring_final_prefix_block2_consume_toy_peak_q == 56
            && direct_restoring_final_prefix_block2_consume_toy_parser_over_node_roundtrip
                < direct_restoring_final_prefix_block2_consume_toy_roundtrip_ratio_budget
            && direct_restoring_final_prefix_block2_consume_toy_parser_scaled_gap < -60_000.0
            && direct_restoring_final_prefix_block2_consume_toy_dirty_restore_cases == 0
            && direct_restoring_final_prefix_block2_consume_toy_dirty_history_cases == 0
            && direct_restoring_final_prefix_block2_consume_toy_dirty_phase_cases == 0,
        "prefix block2 consume/uncompute toy changed; revisit transient parser output cleanup"
    );
    assert!(
        direct_restoring_final_prefix_block2_leaf_touch_toy_decode_forward_ccx == 28
            && direct_restoring_final_prefix_block2_leaf_touch_toy_leaf_touch_ccx == 40
            && direct_restoring_final_prefix_block2_leaf_touch_toy_parser_transient_ccx == 56
            && direct_restoring_final_prefix_block2_leaf_touch_toy_total_ccx == 96
            && direct_restoring_final_prefix_block2_leaf_touch_toy_peak_q == 43
            && direct_restoring_final_prefix_block2_leaf_touch_toy_parser_over_node_roundtrip
                < direct_restoring_final_prefix_block2_leaf_touch_toy_roundtrip_ratio_budget
            && direct_restoring_final_prefix_block2_leaf_touch_toy_total_over_node_roundtrip
                < direct_restoring_final_prefix_block2_leaf_touch_toy_roundtrip_ratio_budget
            && direct_restoring_final_prefix_block2_leaf_touch_toy_parser_scaled_gap < -70_000.0
            && direct_restoring_final_prefix_block2_leaf_touch_toy_total_scaled_gap < -35_000.0
            && direct_restoring_final_prefix_block2_leaf_touch_toy_dirty_restore_cases == 0
            && direct_restoring_final_prefix_block2_leaf_touch_toy_dirty_history_cases == 0
            && direct_restoring_final_prefix_block2_leaf_touch_toy_dirty_phase_cases == 0,
        "prefix block2 leaf-touch toy changed; revisit parser-to-state integration"
    );
    assert!(
        direct_restoring_final_prefix_block2_selected_addsub_toy_decode_forward_ccx == 28
            && direct_restoring_final_prefix_block2_selected_addsub_toy_select_shift_ccx == 60
            && direct_restoring_final_prefix_block2_selected_addsub_toy_addsub_ccx == 12
            && direct_restoring_final_prefix_block2_selected_addsub_toy_arithmetic_ccx == 72
            && direct_restoring_final_prefix_block2_selected_addsub_toy_parser_transient_ccx == 56
            && direct_restoring_final_prefix_block2_selected_addsub_toy_total_ccx == 128
            && direct_restoring_final_prefix_block2_selected_addsub_toy_peak_q == 56
            && direct_restoring_final_prefix_block2_selected_addsub_toy_parser_over_node_roundtrip
                < direct_restoring_final_prefix_block2_selected_addsub_toy_roundtrip_ratio_budget
            && direct_restoring_final_prefix_block2_selected_addsub_toy_arithmetic_over_node_roundtrip
                < direct_restoring_final_prefix_block2_selected_addsub_toy_roundtrip_ratio_budget
            && direct_restoring_final_prefix_block2_selected_addsub_toy_total_over_node_roundtrip
                < direct_restoring_final_prefix_block2_selected_addsub_toy_roundtrip_ratio_budget
            && direct_restoring_final_prefix_block2_selected_addsub_toy_parser_scaled_gap < -70_000.0
            && direct_restoring_final_prefix_block2_selected_addsub_toy_total_scaled_gap < -20_000.0
            && direct_restoring_final_prefix_block2_selected_addsub_toy_dirty_restore_cases == 0
            && direct_restoring_final_prefix_block2_selected_addsub_toy_dirty_history_cases == 0
            && direct_restoring_final_prefix_block2_selected_addsub_toy_dirty_phase_cases == 0,
        "prefix block2 selected add/sub toy changed; revisit parser-to-arithmetic integration"
    );
    assert!(
        direct_restoring_final_prefix_block2_balanced_selected_addsub_toy_tree_ccx == 12
            && direct_restoring_final_prefix_block2_balanced_selected_addsub_toy_read2_ccx == 6
            && direct_restoring_final_prefix_block2_balanced_selected_addsub_toy_decode_forward_ccx
                == 18
            && direct_restoring_final_prefix_block2_balanced_selected_addsub_toy_select_shift_ccx
                == 60
            && direct_restoring_final_prefix_block2_balanced_selected_addsub_toy_addsub_ccx == 12
            && direct_restoring_final_prefix_block2_balanced_selected_addsub_toy_arithmetic_ccx
                == 72
            && direct_restoring_final_prefix_block2_balanced_selected_addsub_toy_parser_transient_ccx
                == 36
            && direct_restoring_final_prefix_block2_balanced_selected_addsub_toy_total_ccx == 108
            && direct_restoring_final_prefix_block2_balanced_selected_addsub_toy_peak_q == 51
            && direct_restoring_final_prefix_block2_balanced_selected_addsub_toy_total_over_node_roundtrip
                < direct_restoring_final_prefix_block2_balanced_selected_addsub_toy_roundtrip_ratio_budget
            && direct_restoring_final_prefix_block2_balanced_selected_addsub_toy_total_scaled_gap
                < -35_000.0
            && direct_restoring_final_prefix_block2_balanced_selected_addsub_toy_dirty_restore_cases
                == 0
            && direct_restoring_final_prefix_block2_balanced_selected_addsub_toy_dirty_history_cases
                == 0
            && direct_restoring_final_prefix_block2_balanced_selected_addsub_toy_dirty_phase_cases
                == 0,
        "balanced prefix block2 selected add/sub toy changed; revisit selective flattening"
    );
    assert!(
        direct_restoring_final_prefix_block2_selected_addsub_roundtrip_toy_forward_decode_ccx == 28
            && direct_restoring_final_prefix_block2_selected_addsub_roundtrip_toy_reverse_decode_ccx
                == 28
            && direct_restoring_final_prefix_block2_selected_addsub_roundtrip_toy_forward_select_shift_ccx
                == 60
            && direct_restoring_final_prefix_block2_selected_addsub_roundtrip_toy_reverse_select_shift_ccx
                == 60
            && direct_restoring_final_prefix_block2_selected_addsub_roundtrip_toy_forward_addsub_ccx
                == 12
            && direct_restoring_final_prefix_block2_selected_addsub_roundtrip_toy_reverse_addsub_ccx
                == 12
            && direct_restoring_final_prefix_block2_selected_addsub_roundtrip_toy_forward_ccx
                == 128
            && direct_restoring_final_prefix_block2_selected_addsub_roundtrip_toy_reverse_ccx
                == 128
            && direct_restoring_final_prefix_block2_selected_addsub_roundtrip_toy_arithmetic_ccx
                == 144
            && direct_restoring_final_prefix_block2_selected_addsub_roundtrip_toy_parser_transient_ccx
                == 112
            && direct_restoring_final_prefix_block2_selected_addsub_roundtrip_toy_total_ccx == 256
            && direct_restoring_final_prefix_block2_selected_addsub_roundtrip_toy_peak_q == 56
            && direct_restoring_final_prefix_block2_selected_addsub_roundtrip_toy_parser_over_node_roundtrip
                < direct_restoring_final_prefix_block2_selected_addsub_roundtrip_toy_roundtrip_ratio_budget
            && direct_restoring_final_prefix_block2_selected_addsub_roundtrip_toy_arithmetic_over_node_roundtrip
                < direct_restoring_final_prefix_block2_selected_addsub_roundtrip_toy_roundtrip_ratio_budget
            && direct_restoring_final_prefix_block2_selected_addsub_roundtrip_toy_total_over_node_roundtrip
                < direct_restoring_final_prefix_block2_selected_addsub_roundtrip_toy_roundtrip_ratio_budget
            && direct_restoring_final_prefix_block2_selected_addsub_roundtrip_toy_parser_scaled_gap < -70_000.0
            && direct_restoring_final_prefix_block2_selected_addsub_roundtrip_toy_total_scaled_gap < -20_000.0
            && direct_restoring_final_prefix_block2_selected_addsub_roundtrip_toy_dirty_restore_cases == 0
            && direct_restoring_final_prefix_block2_selected_addsub_roundtrip_toy_dirty_history_cases == 0
            && direct_restoring_final_prefix_block2_selected_addsub_roundtrip_toy_dirty_phase_cases == 0,
        "prefix block2 selected add/sub roundtrip toy changed; revisit parser-to-arithmetic cleanup"
    );
    assert!(
        direct_restoring_final_prefix_block2_span24_roundtrip_toy_forward_decode_ccx == 28
            && direct_restoring_final_prefix_block2_span24_roundtrip_toy_reverse_decode_ccx == 28
            && direct_restoring_final_prefix_block2_span24_roundtrip_toy_forward_select_shift_ccx
                == 60
            && direct_restoring_final_prefix_block2_span24_roundtrip_toy_reverse_select_shift_ccx
                == 60
            && direct_restoring_final_prefix_block2_span24_roundtrip_toy_forward_addsub_ccx == 50
            && direct_restoring_final_prefix_block2_span24_roundtrip_toy_reverse_addsub_ccx == 50
            && direct_restoring_final_prefix_block2_span24_roundtrip_toy_total_ccx == 332
            && direct_restoring_final_prefix_block2_span24_roundtrip_toy_peak_q == 113
            && direct_restoring_final_prefix_block2_span24_roundtrip_toy_parser_over_node_roundtrip
                < direct_restoring_final_prefix_block2_span24_roundtrip_toy_roundtrip_ratio_budget
            && direct_restoring_final_prefix_block2_span24_roundtrip_toy_total_over_node_roundtrip
                > direct_restoring_final_prefix_block2_span24_roundtrip_toy_roundtrip_ratio_budget
            && direct_restoring_final_prefix_block2_span24_roundtrip_toy_total_scaled_gap > 0.0
            && direct_restoring_final_prefix_block2_span24_taper_add_per_digit_floor
                > direct_restoring_final_prefix_block2_span24_taper_materialized_full_add_per_digit
            && direct_restoring_final_prefix_block2_span24_taper_scaled_gap
                > direct_restoring_final_prefix_block2_span24_roundtrip_toy_total_scaled_gap
            && direct_restoring_final_prefix_block2_span24_roundtrip_toy_dirty_restore_cases == 0
            && direct_restoring_final_prefix_block2_span24_roundtrip_toy_dirty_history_cases == 0
            && direct_restoring_final_prefix_block2_span24_roundtrip_toy_dirty_phase_cases == 0,
        "span24 prefix selected-addsub accounting changed; revisit low-branch range integration"
    );
    assert!(
        (direct_restoring_final_low_branch_prefix_support_weighted_prefix_node_mean
            - direct_restoring_final_low_branch_align_only_prefix_tree_node_floor_mean)
            .abs()
            < 0.01
            && direct_restoring_final_low_branch_prefix_support_weighted_prefix_node_p99
                == direct_restoring_final_low_branch_align_only_prefix_tree_node_floor_p99
            && direct_restoring_final_low_branch_prefix_support_weighted_total_over_node_roundtrip
                < direct_restoring_final_low_branch_prefix_support_weighted_ratio_budget
            && direct_restoring_final_low_branch_prefix_support_weighted_gap < -30_000.0
            && direct_restoring_final_low_branch_prefix_support_weighted_projected_toffoli
                < 2_670_000.0
            && direct_restoring_final_low_branch_prefix_support_weighted_variable_gap > 0.0
            && direct_restoring_final_low_branch_prefix_support_weighted_variable_offset1_gap
                > 0.0
            && direct_restoring_final_low_branch_prefix_support_weighted_balanced_gap < -40_000.0
            && direct_restoring_final_low_branch_prefix_support_weighted_balanced_prefix_bit_p99
                > 381
            && direct_restoring_final_low_branch_prefix_support_weighted_selective_prefix_bit_p99
                == 381
            && direct_restoring_final_low_branch_prefix_support_weighted_selective_prefix_bit_max
                > direct_restoring_final_low_branch_prefix_support_weighted_selective_prefix_bit_p99
            && direct_restoring_final_low_branch_prefix_support_weighted_selective_prefix_scratch_p99
                == GOOGLE_LOW_QUBIT_SCRATCH
            && direct_restoring_final_low_branch_prefix_support_weighted_selective_flatten_steps
                > 80
            && direct_restoring_final_low_branch_prefix_support_weighted_selective_schedule_codebook_steps
                == 126
            && direct_restoring_final_low_branch_prefix_support_weighted_selective_schedule_max_code_len
                == 13
            && direct_restoring_final_low_branch_prefix_support_weighted_selective_schedule_max_len_classes
                == 12
            && direct_restoring_final_low_branch_prefix_support_weighted_selective_schedule_bit_p99
                == direct_restoring_final_low_branch_prefix_support_weighted_selective_prefix_bit_p99
            && direct_restoring_final_low_branch_prefix_support_weighted_selective_schedule_max_bits
                < 400
            && direct_restoring_final_low_branch_prefix_support_weighted_selective_schedule_decoded_symbols
                > 800_000
            && direct_restoring_final_low_branch_prefix_support_weighted_selective_schedule_prefix_collisions
                == 0
            && direct_restoring_final_low_branch_prefix_support_weighted_selective_schedule_decode_mismatches
                == 0
            && direct_restoring_final_low_branch_prefix_support_weighted_selective_schedule_cursor_mismatches
                == 0
            && direct_restoring_final_low_branch_prefix_support_weighted_selective_total_over_node_roundtrip
                < direct_restoring_final_low_branch_prefix_support_weighted_ratio_budget
            && direct_restoring_final_low_branch_prefix_support_weighted_selective_gap < -40_000.0
            && direct_restoring_final_low_branch_prefix_support_weighted_selective_projected_toffoli
                < 2_660_000.0
            && direct_restoring_final_low_branch_prefix_support_weighted_maxconstrained_prefix_bit_p99
                <= 381
            && direct_restoring_final_low_branch_prefix_support_weighted_maxconstrained_prefix_bit_max
                == 381
            && direct_restoring_final_low_branch_prefix_support_weighted_maxconstrained_prefix_scratch_max
                == GOOGLE_LOW_QUBIT_SCRATCH
            && direct_restoring_final_low_branch_prefix_support_weighted_maxconstrained_flatten_steps
                == 83
            && direct_restoring_final_low_branch_prefix_support_weighted_maxconstrained_trimmed_steps
                == 9
            && direct_restoring_final_low_branch_prefix_support_weighted_maxconstrained_over_budget_rows
                == 0
            && direct_restoring_final_low_branch_prefix_support_weighted_maxconstrained_over_budget_mass
                == 0
            && direct_restoring_final_low_branch_prefix_support_weighted_maxconstrained_schedule_codebook_steps
                == 126
            && direct_restoring_final_low_branch_prefix_support_weighted_maxconstrained_schedule_max_code_len
                == 13
            && direct_restoring_final_low_branch_prefix_support_weighted_maxconstrained_schedule_max_len_classes
                == 12
            && direct_restoring_final_low_branch_prefix_support_weighted_maxconstrained_schedule_bit_p99
                == direct_restoring_final_low_branch_prefix_support_weighted_maxconstrained_prefix_bit_p99
            && direct_restoring_final_low_branch_prefix_support_weighted_maxconstrained_schedule_max_bits
                == direct_restoring_final_low_branch_prefix_support_weighted_maxconstrained_prefix_bit_max
            && direct_restoring_final_low_branch_prefix_support_weighted_maxconstrained_schedule_decoded_symbols
                > 800_000
            && direct_restoring_final_low_branch_prefix_support_weighted_maxconstrained_schedule_prefix_collisions
                == 0
            && direct_restoring_final_low_branch_prefix_support_weighted_maxconstrained_schedule_decode_mismatches
                == 0
            && direct_restoring_final_low_branch_prefix_support_weighted_maxconstrained_schedule_cursor_mismatches
                == 0
            && direct_restoring_final_low_branch_prefix_support_weighted_maxconstrained_dynamic_even_p99
                == 1_807
            && direct_restoring_final_low_branch_prefix_support_weighted_maxconstrained_variable_decode_p99
                == 4_707
            && direct_restoring_final_low_branch_prefix_support_weighted_maxconstrained_total_over_node_roundtrip
                < direct_restoring_final_low_branch_prefix_support_weighted_ratio_budget
            && direct_restoring_final_low_branch_prefix_support_weighted_maxconstrained_gap
                < -30_000.0
            && direct_restoring_final_low_branch_prefix_support_weighted_maxconstrained_projected_toffoli
                < 2_670_000.0
            && direct_restoring_final_peakfit_toy_cases_with_sample_gap == 4
            && direct_restoring_final_peakfit_toy_largest_missing_symbols > 500
            && direct_restoring_final_peakfit_toy_largest_sample_over_budget_traces > 3_000
            && direct_restoring_final_peakfit_toy_largest_exact_over_budget_traces > 3_000
            && direct_restoring_final_low_branch_support_toy_cases_with_missing == 4
            && direct_restoring_final_low_branch_support_toy_largest_missing_symbols == 26
            && direct_restoring_final_low_branch_support_toy_largest_missing_steps == 11
            && direct_restoring_final_low_branch_support_toy_largest_span_gap == 4
            && direct_restoring_final_low_branch_support_toy_largest_exact_span == 16
            && direct_restoring_final_low_branch_interval_toy_guard4_cover_cases == 0
            && direct_restoring_final_low_branch_interval_toy_guard4_largest_missing_symbols
                == 25
            && direct_restoring_final_low_branch_interval_toy_guard4_largest_over_budget_traces
                > 4_000
            && direct_restoring_final_low_branch_interval_toy_guard4_largest_max_bits == 37
            && direct_restoring_final_low_branch_interval_toy_full_cover_cases == 4
            && direct_restoring_final_low_branch_interval_toy_full_fit_cases == 0
            && direct_restoring_final_low_branch_interval_toy_full_largest_over_budget_traces
                > 4_000
            && direct_restoring_final_low_branch_interval_toy_full_largest_max_bits == 37
            && direct_restoring_final_peakfit_holdout_missing_symbols == 182
            && direct_restoring_final_peakfit_holdout_missing_traces == 170
            && direct_restoring_final_peakfit_holdout_over_budget_rows == 7
            && direct_restoring_final_peakfit_holdout_max_seen_bits > 381
            && direct_restoring_final_peakfit_scaled_probe_train_samples == 65_536
            && direct_restoring_final_peakfit_scaled_probe_holdout_samples == 32_768
            && direct_restoring_final_peakfit_scaled_probe_flatten_steps == 45
            && direct_restoring_final_peakfit_scaled_probe_missing_symbols > 100
            && direct_restoring_final_peakfit_scaled_probe_missing_traces > 100
            && direct_restoring_final_peakfit_scaled_probe_over_budget_rows == 0
            && direct_restoring_final_peakfit_scaled_probe_max_seen_bits <= 381
            && direct_restoring_final_peakfit_scaled_probe_gap > -1_000.0
            && direct_restoring_final_peakfit_scaled_probe_gap < 0.0
            && direct_restoring_final_low_branch_prefix_support_weighted_span24_uniform_gap > 0.0
            && direct_restoring_final_low_branch_prefix_support_weighted_span24_symbol_p99 == 1
            && direct_restoring_final_low_branch_prefix_support_weighted_support_noncontig_steps
                == direct_restoring_final_low_branch_align_only_support_noncontig_steps
            && direct_restoring_final_low_branch_prefix_support_weighted_support_max_span
                == direct_restoring_final_low_branch_align_only_support_max_span,
        "support-weighted low-branch span accounting changed; revisit prefix decoder promotion"
    );
    assert!(
        direct_restoring_final_prefix_block2_balanced_family_toy_checked_circuits == 289
            && direct_restoring_final_prefix_block2_balanced_family_toy_simulated_circuits == 49
            && direct_restoring_final_prefix_block2_balanced_family_toy_simulated_cases > 25_000
            && direct_restoring_final_prefix_block2_balanced_family_toy_max_support
                == direct_restoring_final_low_branch_prefix_support_weighted_support_max_symbols
            && direct_restoring_final_prefix_block2_balanced_family_toy_max_tree_ccx == 64
            && direct_restoring_final_prefix_block2_balanced_family_toy_max_read2_ccx == 10
            && direct_restoring_final_prefix_block2_balanced_family_toy_max_decode_forward_ccx == 74
            && direct_restoring_final_prefix_block2_balanced_family_toy_max_select_shift_ccx == 216
            && direct_restoring_final_prefix_block2_balanced_family_toy_max_addsub_ccx == 38
            && direct_restoring_final_prefix_block2_balanced_family_toy_max_total_ccx == 804
            && direct_restoring_final_prefix_block2_balanced_family_toy_max_peak_q <= 190
            && direct_restoring_final_prefix_block2_balanced_family_toy_max_total_over_node_roundtrip
                < 8.0
            && direct_restoring_final_prefix_block2_balanced_family_toy_max_total_scaled_gap
                < -25_000.0
            && direct_restoring_final_prefix_block2_balanced_family_toy_dirty_restore_cases == 0
            && direct_restoring_final_prefix_block2_balanced_family_toy_dirty_history_cases == 0
            && direct_restoring_final_prefix_block2_balanced_family_toy_dirty_phase_cases == 0,
        "balanced support-family prefix decoder evidence drifted"
    );
    assert!(
        direct_restoring_final_coeff_decoder_alignment_degree_n14 + 1 >= 14
            && direct_restoring_final_coeff_decoder_alignment_density_n14 > (1usize << 14) / 4
            && direct_restoring_final_coeff_decoder_alignment_max_n14 + 1 >= 14,
        "restoring-final coefficient decoder alignment metadata stopped looking dense"
    );
    assert!(
        direct_restoring_final_align_entropy_variable_scratch_p99 <= GOOGLE_LOW_QUBIT_SCRATCH
            && direct_restoring_final_align_entropy_variable_scratch_max
                <= GOOGLE_LOW_QUBIT_SCRATCH
            && direct_restoring_final_align_entropy_global_scratch_p99
                <= GOOGLE_LOW_QUBIT_SCRATCH
            && direct_restoring_final_align_entropy_global_scratch_max
                <= GOOGLE_LOW_QUBIT_SCRATCH
            && direct_restoring_final_align_entropy_step_scratch_p99
                <= GOOGLE_LOW_QUBIT_SCRATCH
            && direct_restoring_final_align_entropy_step_scratch_max
                <= GOOGLE_LOW_QUBIT_SCRATCH
            && direct_restoring_final_align_prefix_global_scratch_p99
                > GOOGLE_LOW_QUBIT_SCRATCH
            && direct_restoring_final_align_prefix_step_scratch_p99
                > GOOGLE_LOW_QUBIT_SCRATCH
            && direct_restoring_final_align_entropy_branch_count_p99
                == direct_restoring_final_stored_align_branch_count_p99,
        "restoring-final metadata coding frontier changed; revisit parser route"
    );
    assert!(
        direct_restoring_final_range_parser_model_precision_bits == 13
            && direct_restoring_final_range_parser_live_scratch_p99 <= GOOGLE_LOW_QUBIT_SCRATCH
            && direct_restoring_final_range_parser_state_touch_floor_mean
                > direct_restoring_final_range_parser_oneway_budget
            && direct_restoring_final_range_parser_augmented_mean_gap > 0,
        "restoring-final range-parser hard-piece accounting changed; revisit parser route"
    );
    assert!(
        direct_restoring_final_block_parser_model_precision_bits == 13
            && direct_restoring_final_block_parser_best_block_symbols == 8
            && direct_restoring_final_block_parser_best_live_scratch_p99
                <= GOOGLE_LOW_QUBIT_SCRATCH
            && direct_restoring_final_block_parser_best_touch_floor_mean
                < direct_restoring_final_block_parser_oneway_budget
            && direct_restoring_final_block_parser_best_augmented_gap < 0.0
            && direct_restoring_final_block32_live_scratch_p99
                <= GOOGLE_LOW_QUBIT_SCRATCH
            && direct_restoring_final_block32_touch_floor_mean
                < direct_restoring_final_block_parser_oneway_budget,
        "restoring-final block range-parser lower bound changed; revisit toy parser work"
    );
    assert!(
        direct_restoring_final_block4_touch_floor_mean
            < direct_restoring_final_block_parser_best_touch_floor_mean
            && direct_restoring_final_block4_live_scratch_p99 > GOOGLE_LOW_QUBIT_SCRATCH
            && direct_restoring_final_block5_live_scratch_p99 > GOOGLE_LOW_QUBIT_SCRATCH
            && direct_restoring_final_block6_live_scratch_p99 > GOOGLE_LOW_QUBIT_SCRATCH
            && direct_restoring_final_block7_live_scratch_p99 > GOOGLE_LOW_QUBIT_SCRATCH
            && direct_restoring_final_block4_with_binary_lookup_2x_gap > 0.0
            && direct_restoring_final_block4_with_binary_lookup_2x_gap
                < direct_restoring_final_block_parser_best_with_binary_lookup_2x_gap
            && direct_restoring_final_block7_with_binary_lookup_2x_gap > 0.0
            && direct_restoring_final_block7_with_binary_lookup_2x_gap
                < direct_restoring_final_block_parser_best_with_binary_lookup_2x_gap
            && direct_restoring_final_block4_lookup_multiplier_budget > 1.7
            && direct_restoring_final_block4_lookup_multiplier_budget < 2.0,
        "small restoring-final parser blocks now fit scratch or close the 2x lookup gap; revisit parser packing"
    );
    assert!(
        direct_restoring_final_block_parser_lookup_scan_floor_mean
            > direct_restoring_final_block_parser_oneway_budget
            && direct_restoring_final_block_parser_best_with_lookup_mean
                > direct_restoring_final_block_parser_oneway_budget
            && direct_restoring_final_block_parser_best_with_lookup_gap > 0.0,
        "restoring-final threshold-scan lookup floor now fits; build a block parser toy"
    );
    assert!(
        direct_restoring_final_block_parser_best_with_binary_lookup_mean
            < direct_restoring_final_block_parser_oneway_budget
            && direct_restoring_final_block_parser_best_with_binary_lookup_gap < 0.0
            && direct_restoring_final_block_parser_best_with_binary_lookup_2x_mean
                > direct_restoring_final_block_parser_oneway_budget
            && direct_restoring_final_block_parser_best_with_binary_lookup_2x_gap > 0.0
            && direct_restoring_final_block_parser_align_support_noncontig_steps > 0,
        "restoring-final binary-threshold parser accounting changed; revisit rank-map/parser route"
    );
    assert!(
        direct_restoring_final_block_parser_cond_branch_best_block_symbols == 7
            && direct_restoring_final_block_parser_cond_branch_live_scratch_p99
                < direct_restoring_final_block_parser_best_live_scratch_p99
            && direct_restoring_final_block_parser_cond_branch_live_scratch_p99
                > STRICT_SCRATCH
            && direct_restoring_final_block_parser_cond_branch_binary_lookup_floor_mean
                <= direct_restoring_final_block_parser_binary_lookup_floor_mean
            && direct_restoring_final_block_parser_cond_branch_best_with_binary_lookup_2x_gap
                > 9_000.0,
        "conditional branch model now fits restoring-final strict scratch/binary budget; build parser"
    );
    assert!(
        direct_restoring_final_cond_block6_live_scratch_p99 == GOOGLE_LOW_QUBIT_SCRATCH + 2
            && direct_restoring_final_cond_block6_with_binary_lookup_2x_gap > 0.0
            && direct_restoring_final_cond_block6_with_binary_lookup_2x_gap
                < direct_restoring_final_block_parser_cond_branch_best_with_binary_lookup_2x_gap
            && direct_restoring_final_cond_block6_lookup_multiplier_budget > 1.69
            && direct_restoring_final_cond_block6_lookup_multiplier_budget < 1.71
            && direct_restoring_final_cond_block5_live_scratch_p99
                > direct_restoring_final_cond_block6_live_scratch_p99
            && direct_restoring_final_cond_block4_live_scratch_p99
                > direct_restoring_final_cond_block5_live_scratch_p99,
        "branch-conditioned block6 parser now fits scratch or closes the binary lookup gap; update direct-centered frontier"
    );
    assert!(
        direct_restoring_final_cond_mixed67_best_period == 5
            && direct_restoring_final_cond_mixed67_best_mask == 9
            && direct_restoring_final_cond_mixed67_best_seven_count == 2
            && direct_restoring_final_cond_mixed67_live_scratch_p99
                == GOOGLE_LOW_QUBIT_SCRATCH
            && direct_restoring_final_cond_mixed67_with_binary_lookup_2x_gap > 0.0
            && direct_restoring_final_cond_mixed67_with_binary_lookup_2x_gap
                < direct_restoring_final_block_parser_cond_branch_best_with_binary_lookup_2x_gap
            && direct_restoring_final_cond_mixed67_lookup_multiplier_budget > 1.67
            && direct_restoring_final_cond_mixed67_lookup_multiplier_budget < 1.69,
        "mixed 6/7 branch-conditioned parser now reaches target or changed schedule; update direct-centered frontier"
    );
    assert!(
        direct_restoring_final_cond_mixed4to8_best_period == 4
            && direct_restoring_final_cond_mixed4to8_schedule_code == 8_656
            && direct_restoring_final_cond_mixed4to8_live_scratch_p99
                == GOOGLE_LOW_QUBIT_SCRATCH
            && direct_restoring_final_cond_mixed4to8_with_binary_lookup_2x_gap > 0.0
            && direct_restoring_final_cond_mixed4to8_with_binary_lookup_2x_gap
                < direct_restoring_final_cond_mixed67_with_binary_lookup_2x_gap
            && direct_restoring_final_cond_mixed4to8_lookup_multiplier_budget > 1.68
            && direct_restoring_final_cond_mixed4to8_lookup_multiplier_budget < 1.69,
        "mixed 4..8 branch-conditioned parser now reaches target or changed schedule; update direct-centered frontier"
    );
    assert!(
        direct_restoring_final_cond_mixed4to8_with_block_joint_binary_lookup_2x_gap < 0.0
            && direct_restoring_final_cond_mixed4to8_with_block_joint_scan_lookup_2x_gap
                > 400_000.0
            && direct_restoring_final_cond_mixed4to8_block_joint_support_row_floor
                > direct_restoring_final_block32_qrom_row_floor
            && direct_restoring_final_cond_mixed4to8_block_joint_lookup_multiplier_budget
                > 2.3,
        "block-joint lookup floor no longer has the expected binary-depth opening/full-scan blocker"
    );
    assert!(
        direct_restoring_final_selective_pair_lookup_selected_saving_mean * 20.0
            < direct_restoring_final_selective_pair_lookup_required_saving_mean
            && direct_restoring_final_selective_pair_lookup_gap > 8_000.0
            && direct_restoring_final_selective_pair_lookup_selected_positions > 0
            && direct_restoring_final_selective_pair_lookup_support_rows < 1_000
            && direct_restoring_final_selective_pair_lookup_mean
                > direct_restoring_final_selective_pair_lookup_target_mean,
        "selective adjacent-pair lookup now closes the restoring-final parser gap; revisit rank decoder"
    );
    assert!(
        direct_restoring_final_block_joint_rank_degree_n14 >= 14
            && direct_restoring_final_block_joint_rank_density_n14 > (1usize << 14) / 4
            && direct_restoring_final_block_joint_rank_max_patterns_n14 > 2_000
            && direct_restoring_final_block_joint_rank_support_rows_n14 > 3_000,
        "block-joint rank cleanup became sparse enough to revisit the direct parser decoder"
    );
    assert!(
        direct_restoring_final_block_parser_huffman_lookup_floor_mean
            < direct_restoring_final_block_parser_binary_lookup_floor_mean
            && direct_restoring_final_block_parser_cond_branch_huffman_lookup_floor_mean
                < direct_restoring_final_block_parser_cond_branch_binary_lookup_floor_mean
            && direct_restoring_final_cond_mixed67_with_huffman_lookup_2x_gap < 0.0
            && direct_restoring_final_cond_mixed67_huffman_lookup_multiplier_budget > 2.5,
        "distribution-aware restoring parser lookup floor no longer clears target; update parser target"
    );
    assert!(
        direct_restoring_final_block_parser_cond_branch_lookup_scan_floor_mean
            > direct_restoring_final_block_parser_cond_branch_huffman_lookup_floor_mean
            && direct_restoring_final_cond_mixed67_with_cond_scan_lookup_2x_gap > 100_000.0
            && direct_restoring_final_huffman_tree_toy_compare_ccx == 30
            && direct_restoring_final_huffman_tree_toy_roundtrip_ccx
                == 2 * direct_restoring_final_huffman_tree_toy_forward_ccx
            && direct_restoring_final_huffman_tree_toy_full_over_path_ratio
                > direct_restoring_final_cond_mixed67_huffman_lookup_multiplier_budget
            && direct_restoring_final_huffman_tree_toy_dirty_restore_cases == 0
            && direct_restoring_final_huffman_tree_toy_dirty_history_cases == 0
            && direct_restoring_final_huffman_tree_toy_dirty_phase_cases == 0,
        "coherent restoring parser decision-tree blocker changed; revisit Huffman path floor"
    );
    assert!(
        direct_restoring_final_block_parser_best_qrom_row_floor as f64
            > direct_restoring_final_block_parser_oneway_budget
            && direct_restoring_final_block32_qrom_row_floor as f64
                > direct_restoring_final_block_parser_oneway_budget
            && direct_restoring_final_block_parser_best_qrom_gap > 0.0,
        "restoring-final block-QROM row floor now fits; table decoder may revive parser route"
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
    assert!(
        halfgcd_second_col_fixed_depth64_scratch_p99 < GOOGLE_LOW_QUBIT_SCRATCH
            && halfgcd_second_col_fixed_depth64_prefix_bounded_barrel_bits < 8
            && halfgcd_second_col_fixed_depth64_decoder_bounded_barrel_bits < 8
            && halfgcd_second_col_fixed_depth64_prefix_adversarial_prefix_max_digits > 32
            && halfgcd_second_col_fixed_depth64_prefix_adversarial_decoder_max_digits > 32
            && halfgcd_second_col_fixed_depth64_prefix_adversarial_required_barrel_bits == 8
            && halfgcd_second_col_fixed_depth64_prefix_adversarial_missing_layers == 3
            && halfgcd_second_col_fixed_depth64_prefix_full_domain_avg_gap_floor > 0
            && halfgcd_second_col_fixed_depth64_exact_mean < GOOGLE_LOW_QUBIT_TOFFOLI
            && halfgcd_second_col_fixed_depth64_exact_p99 < GOOGLE_LOW_QUBIT_TOFFOLI
            && halfgcd_second_col_fixed_depth64_exact_tail_floor_gap < 0
            && halfgcd_second_col_fixed_depth64_exact_tail_floor_p99 < GOOGLE_LOW_QUBIT_TOFFOLI
            && halfgcd_second_col_fixed_depth64_tail_bounded_barrel_bits < 8
            && halfgcd_second_col_fixed_depth64_tail_adversarial_q_bits > 32
            && halfgcd_second_col_fixed_depth64_tail_adversarial_required_barrel_bits == 8
            && halfgcd_second_col_fixed_depth64_tail_adversarial_missing_layers == 3
            && halfgcd_second_col_fixed_depth64_tail_full_domain_avg_gap_floor > 0
            && halfgcd_second_col_fixed_depth64_exact_tail_bounded_barrel_gap < 0
            && halfgcd_second_col_fixed_depth64_exact_tail_bounded_barrel_p99
                > GOOGLE_LOW_QUBIT_TOFFOLI
            && halfgcd_second_col_fixed_depth64_exact_tail_bounded_plus_one_width_gap < 0
            && halfgcd_second_col_fixed_depth64_exact_tail_bounded_plus_two_width_gap > 0
            && halfgcd_second_col_fixed_depth64_noscan_tail_bounded_barrel_p99
                < GOOGLE_LOW_QUBIT_TOFFOLI
            && halfgcd_second_col_fixed_depth64_exact_tail_logbarrel_gap > 0
            && halfgcd_second_col_fixed_depth64_exact_tail_logbarrel_p99 > GOOGLE_LOW_QUBIT_TOFFOLI
            && halfgcd_second_col_fixed_depth64_noscan_tail_logbarrel_mean
                < GOOGLE_LOW_QUBIT_TOFFOLI
            && halfgcd_second_col_fixed_depth64_noscan_tail_logbarrel_p99
                > GOOGLE_LOW_QUBIT_TOFFOLI
            && halfgcd_second_col_fixed_depth64_exact_prefix_bounded_tail_logbarrel_gap < 0
            && halfgcd_second_col_fixed_depth64_exact_prefix_bounded_tail_logbarrel_p99
                < GOOGLE_LOW_QUBIT_TOFFOLI
            && halfgcd_second_col_fixed_depth64_exact_prefix_bounded_tail_logbarrel_plus_one_width_gap
                < 0
            && halfgcd_second_col_fixed_depth64_exact_prefix_bounded_tail_logbarrel_plus_one_width_p99
                < GOOGLE_LOW_QUBIT_TOFFOLI
            && halfgcd_second_col_fixed_depth64_exact_prefix_bounded_tail_logbarrel_plus_two_width_gap
                < 0
            && halfgcd_second_col_fixed_depth64_exact_prefix_bounded_tail_logbarrel_plus_two_width_p99
                > GOOGLE_LOW_QUBIT_TOFFOLI,
        "half-GCD fixed-depth64 tail-alignment frontier changed; revisit tail parser priority"
    );
    assert!(
        halfgcd_second_col_fixed_depth64_static_app_mean > GOOGLE_LOW_QUBIT_TOFFOLI
            && halfgcd_second_col_fixed_depth64_static_app_p99 > GOOGLE_LOW_QUBIT_TOFFOLI
            && halfgcd_second_col_fixed_depth64_static_app_gap > 200_000
            && halfgcd_second_col_fixed_depth64_app_static_floor_mean
                > halfgcd_second_col_fixed_depth64_app_popcount_mean
            && halfgcd_second_col_fixed_depth64_app_static_over_popcount_mean
                > halfgcd_second_col_fixed_depth64_exact_tail_logbarrel_gap as usize,
        "half-GCD coefficient application control accounting changed; revisit whether popcount application can be promoted"
    );
    assert!(
        halfgcd_second_col_fixed_depth64_static_sep4_app_gap < 0
            && halfgcd_second_col_fixed_depth64_static_sep4_app_p99 < GOOGLE_LOW_QUBIT_TOFFOLI
            && halfgcd_second_col_fixed_depth64_static_joint4_app_gap < 0
            && halfgcd_second_col_fixed_depth64_static_joint4_app_p99 < GOOGLE_LOW_QUBIT_TOFFOLI
            && halfgcd_second_col_fixed_depth64_static_sep4_selector_budget_oneway < 50_000
            && halfgcd_second_col_fixed_depth64_static_joint4_selector_budget_oneway < 80_000
            && halfgcd_second_col_fixed_depth64_app_static_joint4_floor_mean
                < halfgcd_second_col_fixed_depth64_app_static_sep4_floor_mean
            && halfgcd_second_col_fixed_depth64_app_static_sep4_floor_mean
                < halfgcd_second_col_fixed_depth64_app_popcount_mean,
        "half-GCD static-window coefficient floor changed; revisit selector/precompute budget"
    );
    assert!(
        halfgcd_second_col_fixed_depth64_static_joint4_with_selector_floor_gap > 0
            && halfgcd_second_col_fixed_depth64_static_joint4_with_selector_floor_p99
                > GOOGLE_LOW_QUBIT_TOFFOLI
            && halfgcd_second_col_fixed_depth64_static_sep4_with_selector_floor_mean
                > halfgcd_second_col_fixed_depth64_static_joint4_with_selector_floor_mean
            && halfgcd_second_col_fixed_depth64_app_static_selector_floor_mean
                > halfgcd_second_col_fixed_depth64_static_joint4_selector_budget_oneway
            && halfgcd_second_col_fixed_depth64_app_static_selector_floor_over_joint4_budget
                > 30_000,
        "half-GCD static-window bit-product selector floor changed; revisit whether a real selector can fit"
    );
    assert!(
        halfgcd_second_col_fixed_depth64_static_window_scan_best_w == 6
            && halfgcd_second_col_fixed_depth64_static_window_scan_best_gap > 0
            && halfgcd_second_col_fixed_depth64_static_window_scan_best_mean
                < halfgcd_second_col_fixed_depth64_static_joint4_with_selector_floor_mean
            && halfgcd_second_col_fixed_depth64_static_window_scan_best_selector_mean
                == halfgcd_second_col_fixed_depth64_app_static_selector_floor_mean
            && halfgcd_second_col_fixed_depth64_static_window_scan_best_table_row_mean
                > 70_000,
        "half-GCD static-window scan now fits or changed selector/table tradeoff; revisit route"
    );
    assert!(
        halfgcd_second_col_fixed_depth64_static_window_table_only_best_w == 4
            && halfgcd_second_col_fixed_depth64_static_window_table_only_best_gap < 0
            && halfgcd_second_col_fixed_depth64_static_window_required_selector_mean
                < halfgcd_second_col_fixed_depth64_static_window_scan_best_selector_mean
            && halfgcd_second_col_fixed_depth64_static_window_required_selector_mean
                > halfgcd_second_col_fixed_depth64_static_window_scan_best_table_row_mean
            && halfgcd_second_col_fixed_depth64_static_window_selector_cut_needed > 20_000
            && halfgcd_second_col_fixed_depth64_static_window_table_margin > 5_000,
        "half-GCD static-window selector breakthrough budget changed; revisit selector design"
    );
    assert!(
        halfgcd_second_col_fixed_depth64_static_window_table_source_best_w == 2
            && halfgcd_second_col_fixed_depth64_static_window_table_source_gap > 1_000_000
            && halfgcd_second_col_fixed_depth64_static_window_table_source_product_floor_mean
                > halfgcd_second_col_fixed_depth64_static_window_source_product_floor_mean,
        "half-GCD table row source-products now fit; revisit static table selector design"
    );
    assert!(
        halfgcd_second_col_fixed_depth64_static_window_source_product_best_w == 6
            && halfgcd_second_col_fixed_depth64_static_window_source_product_gap
                > halfgcd_second_col_fixed_depth64_static_window_scan_best_gap
            && halfgcd_second_col_fixed_depth64_static_window_source_product_floor_mean
                > halfgcd_second_col_fixed_depth64_static_window_required_selector_mean
            && halfgcd_second_col_fixed_depth64_static_window_source_product_floor_mean
                > halfgcd_second_col_fixed_depth64_static_window_source_product_table_row_mean,
        "half-GCD static-window source-product selector floor now fits; build the structural selector"
    );
    assert!(
        halfgcd_second_col_fixed_depth64_static_window_wnaf_best_w == 6
            && halfgcd_second_col_fixed_depth64_static_window_wnaf_gap > 0
            && halfgcd_second_col_fixed_depth64_static_window_wnaf_gap
                < halfgcd_second_col_fixed_depth64_static_window_source_product_gap
            && halfgcd_second_col_fixed_depth64_static_window_wnaf_selector_floor_mean
                > halfgcd_second_col_fixed_depth64_static_window_required_selector_mean
            && halfgcd_second_col_fixed_depth64_static_window_wnaf_source_product_floor_mean
                > halfgcd_second_col_fixed_depth64_static_window_wnaf_table_row_floor_mean,
        "half-GCD sparse signed-window selector floor now fits; revisit wNAF coefficient application"
    );
    assert!(
        halfgcd_second_col_fixed_depth64_static_window_wnaf_compact_best_w == 2
            && halfgcd_second_col_fixed_depth64_static_window_wnaf_compact_gap < 0
            && halfgcd_second_col_fixed_depth64_static_window_wnaf_compact_missing_active_floor_mean
                > halfgcd_second_col_fixed_depth64_static_window_wnaf_compact_active_slack_oneway
            && halfgcd_second_col_fixed_depth64_static_window_wnaf_compact_source_product_floor_mean
                > halfgcd_second_col_fixed_depth64_static_window_wnaf_compact_table_row_floor_mean,
        "half-GCD compact NAF active predicate now fits; build the recoder-cleanup toy"
    );
    assert!(
        halfgcd_second_col_compact_wnaf_active_degree_n14 == 14
            && halfgcd_second_col_compact_wnaf_active_density_n14 > 8_000
            && halfgcd_second_col_compact_wnaf_active_pair_positions_n14
                == halfgcd_second_col_compact_wnaf_active_positions_n14
            && halfgcd_second_col_compact_wnaf_active_slots_n14 + 1
                >= halfgcd_second_col_compact_wnaf_active_full_slots_n14
            && halfgcd_second_col_compact_wnaf_active_max_pair_n14 == 3,
        "half-GCD compact NAF active predicate became structurally cheap; revisit wNAF recoding"
    );
    assert!(
        halfgcd_second_col_fixed_depth64_dynamic_barrel_static_mean
            > GOOGLE_LOW_QUBIT_TOFFOLI
            && halfgcd_second_col_fixed_depth64_dynamic_barrel_mean
                < GOOGLE_LOW_QUBIT_TOFFOLI
            && halfgcd_second_col_fixed_depth64_dynamic_barrel_p99
                < GOOGLE_LOW_QUBIT_TOFFOLI
            && halfgcd_second_col_fixed_depth64_dynamic_barrel_savings_mean
                > halfgcd_second_col_fixed_depth64_exact_tail_logbarrel_gap as usize
            && halfgcd_second_col_fixed_depth64_dynamic_barrel_scratch_p99
                == halfgcd_second_col_fixed_depth64_scratch_p99
            && halfgcd_second_col_alignment_mbu_degree_n14 >= 14
            && halfgcd_second_col_alignment_mbu_density_n14 > (1usize << 14) / 4,
        "half-GCD dynamic barrel/classical-control frontier changed; revisit measurement-clean alignment priority"
    );
    assert!(
        halfgcd_second_col_static_window_mbu_max_coeff_bits_n14 >= 14
            && halfgcd_second_col_static_window_mbu_max_pair_n14 == 63
            && halfgcd_second_col_static_window_mbu_degree_n14 + 2 >= 14
            && halfgcd_second_col_static_window_mbu_density_n14 > (1usize << 14) / 4,
        "half-GCD static-window selector bits may be generic-MBU clean; revisit table-only route"
    );
    assert!(
        halfgcd_second_col_static_window_wnaf_mbu_max_positions_n14 >= 14
            && halfgcd_second_col_static_window_wnaf_mbu_max_pair_n14 == 63
            && halfgcd_second_col_static_window_wnaf_mbu_degree_n14 + 2 >= 14
            && halfgcd_second_col_static_window_wnaf_mbu_density_n14 > (1usize << 14) / 4,
        "half-GCD wNAF static-window selector bits may be generic-MBU clean; revisit table-only route"
    );
    assert!(
        halfgcd_second_col_static_window_support_rows_n14 * 5
            > halfgcd_second_col_static_window_support_full_rows_n14 * 3
            && halfgcd_second_col_static_window_support_saturated_windows_n14 > 0
            && halfgcd_second_col_static_window_support_saturated_windows_n14
                < halfgcd_second_col_static_window_support_windows_n14
            && halfgcd_second_col_static_window_bit_support_n14 + 2
                >= halfgcd_second_col_static_window_full_bits_n14
            && halfgcd_second_col_static_window_bit_support_ppm_n14 > 900_000,
        "half-GCD static-window support changed enough to revisit source-product selector route"
    );
}
