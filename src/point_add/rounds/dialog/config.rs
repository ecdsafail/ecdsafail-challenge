
use super::*;

pub const DIALOG_GCD_ACTIVE_ITERATIONS_ENV: &str = "DIALOG_GCD_ACTIVE_ITERATIONS";
pub const DIALOG_GCD_COMPARE_BITS_ENV: &str = "DIALOG_GCD_COMPARE_BITS";
pub const DIALOG_GCD_PA9024_COMPARE_SCHEDULE_ENV: &str = "DIALOG_GCD_PA9024_COMPARE_SCHEDULE";
pub const DIALOG_GCD_PA9024_COMPARE_SCHEDULE_FLOOR_ENV: &str =
    "DIALOG_GCD_PA9024_COMPARE_SCHEDULE_FLOOR";
pub const DIALOG_GCD_APPLY_CLEAN_COMPARE_BITS_ENV: &str = "DIALOG_GCD_APPLY_CLEAN_COMPARE_BITS";
pub const DIALOG_GCD_COMPRESSED_SIDECAR_LOG_ENV: &str = "DIALOG_GCD_COMPRESSED_SIDECAR_LOG";
pub const DIALOG_GCD_COMPRESSED_BLOCK_LIFECYCLE_ENV: &str = "DIALOG_GCD_COMPRESSED_BLOCK_LIFECYCLE";
pub const DIALOG_GCD_RAW_APPLY_DIRECT_SPECIAL_ADD_ENV: &str =
    "DIALOG_GCD_RAW_APPLY_DIRECT_SPECIAL_ADD";
pub const DIALOG_GCD_RAW_APPLY_MATERIALIZED_SPECIAL_ADD_ENV: &str =
    "DIALOG_GCD_RAW_APPLY_MATERIALIZED_SPECIAL_ADD";
pub const DIALOG_GCD_RAW_APPLY_REVERSE_FAST_SUB_ENV: &str = "DIALOG_GCD_RAW_APPLY_REVERSE_FAST_SUB";
pub const DIALOG_GCD_RAW_APPLY_REVERSE_MATERIALIZED_SPECIAL_SUB_ENV: &str =
    "DIALOG_GCD_RAW_APPLY_REVERSE_MATERIALIZED_SPECIAL_SUB";
pub const DIALOG_GCD_RAW_TOBITVECTOR_MATERIALIZED_SUB_ENV: &str =
    "DIALOG_GCD_RAW_TOBITVECTOR_MATERIALIZED_SUB";
pub const DIALOG_GCD_RAW_TOBITVECTOR_VARIABLE_WIDTH_ENV: &str =
    "DIALOG_GCD_RAW_TOBITVECTOR_VARIABLE_WIDTH";
pub const DIALOG_GCD_RAW_TOBITVECTOR_BORROW_FUTURE_LOG_CARRIES_ENV: &str =
    "DIALOG_GCD_RAW_TOBITVECTOR_BORROW_FUTURE_LOG_CARRIES";
pub const DIALOG_GCD_RAW_IPMUL_TERMINAL_REUSE_ENV: &str = "DIALOG_GCD_RAW_IPMUL_TERMINAL_REUSE";
pub const DIALOG_GCD_RAW_IPMUL_CLEAR_P_RESIDUAL_ENV: &str = "DIALOG_GCD_RAW_IPMUL_CLEAR_P_RESIDUAL";
pub const DIALOG_GCD_RAW_QUOTIENT_TERMINAL_REUSE_ENV: &str =
    "DIALOG_GCD_RAW_QUOTIENT_TERMINAL_REUSE";
pub const DIALOG_GCD_RAW_QUOTIENT_KEEP_TERMINAL_U_ENV: &str =
    "DIALOG_GCD_RAW_QUOTIENT_KEEP_TERMINAL_U";
pub const DIALOG_GCD_RAW_APPLY_TRUNCATED_CLEAN_ENV: &str = "DIALOG_GCD_RAW_APPLY_TRUNCATED_CLEAN";
pub const DIALOG_GCD_RAW_PA_ENV: &str = "DIALOG_GCD_RAW_PA";
pub const DIALOG_GCD_RAW_PA_STOP_AFTER_QUOTIENT_ENV: &str = "DIALOG_GCD_RAW_PA_STOP_AFTER_QUOTIENT";
pub const DIALOG_GCD_RAW_PA_STOP_AFTER_XTAIL_ENV: &str = "DIALOG_GCD_RAW_PA_STOP_AFTER_XTAIL";
pub const DIALOG_GCD_RAW_PA_STOP_AFTER_C_ENV: &str = "DIALOG_GCD_RAW_PA_STOP_AFTER_C";
pub const DIALOG_GCD_RAW_PA_STOP_AFTER_PAIR2_ENV: &str = "DIALOG_GCD_RAW_PA_STOP_AFTER_PAIR2";

pub(crate) fn dialog_gcd_raw_apply_direct_special_add_enabled() -> bool {
    std::env::var(DIALOG_GCD_RAW_APPLY_DIRECT_SPECIAL_ADD_ENV)
        .ok()
        .as_deref()
        == Some("1")
}

pub(crate) fn dialog_gcd_raw_apply_materialized_special_add_enabled() -> bool {
    std::env::var(DIALOG_GCD_RAW_APPLY_MATERIALIZED_SPECIAL_ADD_ENV)
        .ok()
        .as_deref()
        == Some("1")
}

pub(crate) fn dialog_gcd_raw_apply_reverse_fast_sub_enabled() -> bool {
    std::env::var(DIALOG_GCD_RAW_APPLY_REVERSE_FAST_SUB_ENV)
        .ok()
        .as_deref()
        == Some("1")
}

pub(crate) fn dialog_gcd_raw_apply_reverse_materialized_special_sub_enabled() -> bool {
    std::env::var(DIALOG_GCD_RAW_APPLY_REVERSE_MATERIALIZED_SPECIAL_SUB_ENV)
        .ok()
        .as_deref()
        == Some("1")
}

pub(crate) fn dialog_gcd_apply_chunked_f_blocks() -> Option<usize> {
    std::env::var("DIALOG_GCD_APPLY_CHUNKED_F_BLOCKS")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .filter(|&blocks| blocks >= 2)
}

pub(crate) fn dialog_gcd_apply_chunked_f_cuts() -> Option<Vec<usize>> {
    std::env::var("DIALOG_GCD_APPLY_CHUNKED_F_CUTS")
        .ok()
        .filter(|value| !value.trim().is_empty())
        .map(|value| {
            value
                .split(',')
                .map(|item| {
                    item.trim()
                        .parse::<usize>()
                        .expect("DIALOG_GCD_APPLY_CHUNKED_F_CUTS")
                })
                .collect()
        })
}

pub(crate) fn dialog_gcd_apply_chunked_f_cut() -> Option<usize> {
    std::env::var("DIALOG_GCD_APPLY_CHUNKED_F_CUT")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .filter(|&cut| (1..N).contains(&cut))
}

pub(crate) fn dialog_gcd_apply_chunked_f_cut2() -> Option<usize> {
    std::env::var("DIALOG_GCD_APPLY_CHUNKED_F_CUT2")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .filter(|&cut| (1..N).contains(&cut))
}

pub(crate) fn dialog_gcd_apply_chunked_f_cut3() -> Option<usize> {
    std::env::var("DIALOG_GCD_APPLY_CHUNKED_F_CUT3")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .filter(|&cut| (1..N).contains(&cut))
}

pub(crate) fn dialog_gcd_apply_chunked_f_cut4() -> Option<usize> {
    std::env::var("DIALOG_GCD_APPLY_CHUNKED_F_CUT4")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .filter(|&cut| (1..N).contains(&cut))
}

pub(crate) fn dialog_gcd_apply_chunked_f_custom4_enabled() -> bool {
    std::env::var("DIALOG_GCD_APPLY_CHUNKED_F_CUSTOM4")
        .ok()
        .as_deref()
        == Some("1")
}

pub(crate) fn dialog_gcd_apply_chunked_f_custom5_enabled() -> bool {
    std::env::var("DIALOG_GCD_APPLY_CHUNKED_F_CUSTOM5")
        .ok()
        .as_deref()
        == Some("1")
}

pub(crate) fn dialog_gcd_apply_chunked_f_reuse_cin_zero_enabled() -> bool {
    std::env::var("DIALOG_GCD_APPLY_CHUNKED_F_REUSE_CIN_ZERO")
        .ok()
        .as_deref()
        != Some("0")
}

pub(crate) fn dialog_gcd_apply_chunked_f_fuse_boundary_clears_enabled() -> bool {
    std::env::var("DIALOG_GCD_APPLY_CHUNKED_F_FUSE_BOUNDARY_CLEARS")
        .ok()
        .as_deref()
        != Some("0")
}

pub(crate) fn dialog_gcd_apply_borrow_future_boundary_carries_enabled() -> bool {
    std::env::var("DIALOG_GCD_APPLY_BORROW_FUTURE_BOUNDARY_CARRIES")
        .ok()
        .as_deref()
        == Some("1")
}

pub(crate) fn dialog_gcd_apply_boundary_free_owned_during_replay_enabled() -> bool {
    std::env::var("DIALOG_GCD_APPLY_BOUNDARY_FREE_OWNED_DURING_REPLAY")
        .ok()
        .as_deref()
        == Some("1")
}

pub(crate) fn dialog_gcd_apply_implicit_high_zero_enabled() -> bool {
    std::env::var("DIALOG_GCD_APPLY_IMPLICIT_HIGH_ZERO")
        .ok()
        .as_deref()
        == Some("1")
}

pub(crate) fn dialog_gcd_apply_chunked_f_auto_topclean_target() -> Option<u32> {
    std::env::var("DIALOG_GCD_APPLY_CHUNKED_F_AUTO_TOPCLEAN_TARGET")
        .ok()
        .and_then(|value| value.parse::<u32>().ok())
        .filter(|&target| target > 0)
}

pub(crate) fn dialog_gcd_apply_chunked_f_auto_topclean_max_bits() -> usize {
    std::env::var("DIALOG_GCD_APPLY_CHUNKED_F_AUTO_TOPCLEAN_MAX_BITS")
        .ok()
        .and_then(|value| value.parse::<usize>().ok())
        .unwrap_or(1)
}

pub(crate) fn dialog_gcd_apply_final_lowq_enabled() -> bool {
    std::env::var("DIALOG_GCD_APPLY_FINAL_LOWQ").ok().as_deref() == Some("1")
}

pub(crate) fn dialog_gcd_fused_hclear_measured_enabled() -> bool {
    std::env::var("DIALOG_GCD_FUSED_HCLEAR_MEASURED").ok().as_deref() == Some("1")
}

pub(crate) fn dialog_gcd_fused_dclear_measured_enabled() -> bool {
    std::env::var("DIALOG_GCD_FUSED_DCLEAR_MEASURED").ok().as_deref() == Some("1")
}

pub(crate) fn dialog_gcd_fused_ovfclear_measured_enabled() -> bool {
    std::env::var("DIALOG_GCD_FUSED_OVFCLEAR_MEASURED").ok().as_deref() == Some("1")
}

pub(crate) fn dialog_gcd_fused_halve_edclear_measured_enabled() -> bool {
    std::env::var("DIALOG_GCD_FUSED_HALVE_EDCLEAR_MEASURED").ok().as_deref() == Some("1")
}

pub(crate) fn dialog_gcd_apply_final_windowed_fast_blocks() -> Option<usize> {
    std::env::var("DIALOG_GCD_APPLY_FINAL_WINDOWED_FAST_BLOCKS")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .filter(|&blocks| blocks >= 2)
}

pub(crate) fn dialog_gcd_apply_final_topclean_bits() -> usize {
    std::env::var("DIALOG_GCD_APPLY_FINAL_TOPCLEAN")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(0)
}

pub(crate) fn dialog_gcd_apply_boundary_split() -> Option<usize> {
    std::env::var("DIALOG_GCD_APPLY_BOUNDARY_SPLIT")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .filter(|&split| split > 0)
}

pub(crate) fn dialog_gcd_apply_boundary_conditional_replay_enabled() -> bool {
    std::env::var("DIALOG_GCD_APPLY_BOUNDARY_CONDITIONAL_REPLAY")
        .ok()
        .as_deref()
        == Some("1")
}

pub(crate) fn dialog_gcd_reverse_branch_conditional_replay_enabled() -> bool {
    std::env::var("DIALOG_GCD_REVERSE_BRANCH_CONDITIONAL_REPLAY")
        .ok()
        .as_deref()
        == Some("1")
}

pub(crate) fn dialog_gcd_special_clean_conditional_replay_enabled() -> bool {
    std::env::var("DIALOG_GCD_SPECIAL_CLEAN_CONDITIONAL_REPLAY")
        .ok()
        .as_deref()
        == Some("1")
}

pub(crate) fn dialog_gcd_apply_replay_swap_host_enabled() -> bool {

    std::env::var("DIALOG_GCD_APPLY_REPLAY_SWAP_HOST")
        .ok()
        .as_deref()
        == Some("1")
}

pub(crate) fn dialog_gcd_raw_tobitvector_materialized_sub_enabled() -> bool {
    std::env::var(DIALOG_GCD_RAW_TOBITVECTOR_MATERIALIZED_SUB_ENV)
        .ok()
        .as_deref()
        == Some("1")
}

pub(crate) fn dialog_gcd_ctrl_body_vented_enabled() -> bool {
    std::env::var("DIALOG_GCD_CTRL_BODY_VENTED")
        .ok()
        .as_deref()
        != Some("0")
}

pub(crate) fn dialog_gcd_raw_tobitvector_variable_width_enabled() -> bool {
    std::env::var(DIALOG_GCD_RAW_TOBITVECTOR_VARIABLE_WIDTH_ENV)
        .ok()
        .as_deref()
        == Some("1")
}

pub(crate) fn dialog_gcd_raw_tobitvector_borrow_future_log_carries_enabled() -> bool {
    std::env::var(DIALOG_GCD_RAW_TOBITVECTOR_BORROW_FUTURE_LOG_CARRIES_ENV)
        .ok()
        .as_deref()
        == Some("1")
}

pub(crate) fn dialog_gcd_raw_ipmul_terminal_reuse_enabled() -> bool {
    std::env::var(DIALOG_GCD_RAW_IPMUL_TERMINAL_REUSE_ENV)
        .ok()
        .as_deref()
        == Some("1")
}

pub(crate) fn dialog_gcd_raw_ipmul_clear_p_residual_enabled() -> bool {
    std::env::var(DIALOG_GCD_RAW_IPMUL_CLEAR_P_RESIDUAL_ENV)
        .ok()
        .as_deref()
        == Some("1")
}

pub(crate) fn dialog_gcd_raw_quotient_terminal_reuse_enabled() -> bool {
    if let Ok(value) = std::env::var(DIALOG_GCD_RAW_QUOTIENT_TERMINAL_REUSE_ENV) {
        return value == "1";
    }
    dialog_gcd_raw_ipmul_terminal_reuse_enabled()
}

pub(crate) fn dialog_gcd_raw_quotient_keep_terminal_u_enabled() -> bool {
    std::env::var(DIALOG_GCD_RAW_QUOTIENT_KEEP_TERMINAL_U_ENV)
        .ok()
        .as_deref()
        == Some("1")
}

pub(crate) fn dialog_gcd_raw_apply_truncated_clean_enabled() -> bool {
    std::env::var(DIALOG_GCD_RAW_APPLY_TRUNCATED_CLEAN_ENV)
        .ok()
        .as_deref()
        == Some("1")
}

pub(crate) fn dialog_gcd_raw_pa_stop_after_quotient_enabled() -> bool {
    std::env::var(DIALOG_GCD_RAW_PA_STOP_AFTER_QUOTIENT_ENV)
        .ok()
        .as_deref()
        == Some("1")
}

pub(crate) fn dialog_gcd_raw_pa_stop_after_xtail_enabled() -> bool {
    std::env::var(DIALOG_GCD_RAW_PA_STOP_AFTER_XTAIL_ENV)
        .ok()
        .as_deref()
        == Some("1")
}

pub(crate) fn dialog_gcd_raw_pa_stop_after_c_enabled() -> bool {
    std::env::var(DIALOG_GCD_RAW_PA_STOP_AFTER_C_ENV)
        .ok()
        .as_deref()
        == Some("1")
}

pub(crate) fn dialog_gcd_raw_pa_stop_after_pair2_enabled() -> bool {
    std::env::var(DIALOG_GCD_RAW_PA_STOP_AFTER_PAIR2_ENV)
        .ok()
        .as_deref()
        == Some("1")
}

pub(crate) const DIALOG_GCD_MAX_ITERATIONS: usize = 402;
pub(crate) const DIALOG_GCD_RAW_LOG_BITS: usize = 2 * DIALOG_GCD_MAX_ITERATIONS;
pub(crate) const DIALOG_GCD_SPECIAL_ADD_LSBS: usize = 73;
pub(crate) const DIALOG_GCD_DEFAULT_COMPARE_BITS: usize = 77;
pub(crate) const DIALOG_GCD_HIGH_TAIL_ALIAS_GROUP_SIZE: usize = 3;
pub(crate) const DIALOG_GCD_HIGH_TAIL_ALIAS_BLOCK_BITS: usize = 5;

pub(crate) fn dialog_gcd_compressed_sidecar_log_enabled() -> bool {
    std::env::var(DIALOG_GCD_COMPRESSED_SIDECAR_LOG_ENV)
        .ok()
        .as_deref()
        == Some("1")
}

pub(crate) fn dialog_gcd_compressed_block_lifecycle_enabled() -> bool {
    if dialog_gcd_k5_clean_block_enabled() {
        return true;
    }
    std::env::var(DIALOG_GCD_COMPRESSED_BLOCK_LIFECYCLE_ENV)
        .ok()
        .as_deref()
        == Some("1")
}

pub(crate) fn dialog_gcd_k2_enabled() -> bool {
    std::env::var("DIALOG_GCD_K2").ok().as_deref() == Some("1")
}

pub(crate) fn dialog_gcd_k5_clean_block_enabled() -> bool {
    dialog_gcd_k2_enabled()
        && std::env::var("DIALOG_GCD_K5_CLEAN_BLOCK")
            .ok()
            .as_deref()
            == Some("1")
}

pub(crate) fn dialog_gcd_block_bits() -> usize {
    if dialog_gcd_k5_clean_block_enabled() {
        12
    } else if dialog_gcd_k2_pair_compress_enabled() {

        DIALOG_GCD_HIGH_TAIL_ALIAS_BLOCK_BITS
    } else if dialog_gcd_k2_enabled() {
        DIALOG_GCD_HIGH_TAIL_ALIAS_BLOCK_BITS + DIALOG_GCD_HIGH_TAIL_ALIAS_GROUP_SIZE
    } else {
        DIALOG_GCD_HIGH_TAIL_ALIAS_BLOCK_BITS
    }
}

pub(crate) fn dialog_gcd_raw_block_len() -> usize {
    if dialog_gcd_k2_enabled() {
        3 * dialog_gcd_sidecar_group_size()
    } else {
        2 * dialog_gcd_sidecar_group_size()
    }
}

pub const DIALOG_GCD_PA9024_COMPARE_SCHEDULE: [usize; 258] = [
    22, 21, 24, 24, 28, 25, 29, 26, 29, 30, 33, 35, 31, 32, 31, 33, 33, 34, 30, 32, 33, 35, 33, 35,
    34, 33, 35, 35, 35, 34, 33, 33, 33, 34, 34, 38, 35, 35, 33, 36, 34, 36, 37, 36, 38, 36, 38, 36,
    42, 36, 37, 37, 39, 36, 36, 39, 37, 39, 35, 38, 38, 38, 36, 44, 37, 38, 36, 39, 41, 38, 37, 41,
    40, 35, 36, 37, 41, 38, 38, 38, 37, 37, 39, 37, 37, 37, 38, 39, 38, 37, 42, 40, 38, 38, 39, 43,
    41, 39, 40, 42, 40, 39, 44, 39, 44, 40, 43, 40, 40, 41, 42, 41, 43, 42, 45, 41, 43, 42, 43, 42,
    43, 42, 44, 42, 46, 41, 44, 42, 44, 42, 43, 46, 44, 43, 48, 50, 48, 44, 44, 44, 55, 46, 46, 44,
    43, 49, 44, 45, 44, 48, 44, 46, 45, 46, 45, 44, 45, 46, 48, 46, 45, 46, 50, 48, 44, 47, 47, 46,
    45, 46, 45, 48, 47, 49, 47, 47, 46, 49, 48, 49, 46, 48, 50, 51, 47, 54, 49, 48, 47, 48, 51, 50,
    53, 54, 50, 52, 50, 51, 53, 52, 49, 52, 50, 52, 49, 52, 49, 53, 51, 55, 52, 51, 51, 51, 49, 47,
    47, 45, 45, 43, 43, 41, 41, 39, 39, 37, 37, 35, 35, 33, 33, 31, 31, 29, 29, 27, 27, 25, 25, 23,
    23, 21, 21, 19, 19, 17, 17, 15, 15, 13, 13, 11, 11, 9, 9, 7, 7, 5,
];

pub(crate) fn dialog_gcd_active_iterations() -> usize {
    std::env::var(DIALOG_GCD_ACTIVE_ITERATIONS_ENV)
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .filter(|&iters| (1..=DIALOG_GCD_MAX_ITERATIONS).contains(&iters))
        .unwrap_or(DIALOG_GCD_MAX_ITERATIONS)
}

pub(crate) fn dialog_gcd_compare_bits() -> usize {
    std::env::var(DIALOG_GCD_COMPARE_BITS_ENV)
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .filter(|&bits| (1..=N).contains(&bits))
        .unwrap_or(DIALOG_GCD_DEFAULT_COMPARE_BITS)
}

pub(crate) fn dialog_gcd_apply_clean_compare_bits() -> usize {
    std::env::var(DIALOG_GCD_APPLY_CLEAN_COMPARE_BITS_ENV)
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .filter(|&bits| (1..=N).contains(&bits))
        .unwrap_or_else(dialog_gcd_compare_bits)
}

pub(crate) fn dialog_gcd_pa9024_compare_schedule_enabled() -> bool {
    std::env::var(DIALOG_GCD_PA9024_COMPARE_SCHEDULE_ENV)
        .ok()
        .as_deref()
        == Some("1")
}

pub(crate) fn dialog_gcd_pa9024_compare_schedule_floor() -> usize {
    std::env::var(DIALOG_GCD_PA9024_COMPARE_SCHEDULE_FLOOR_ENV)
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .filter(|&bits| bits <= N)
        .unwrap_or(1)
        .max(1)
}

pub(crate) fn dialog_gcd_pa9024_compare_schedule_margin() -> usize {
    std::env::var("DIALOG_GCD_PA9024_COMPARE_SCHEDULE_MARGIN")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(0)
}

fn dialog_gcd_compare_step_bits(step: usize) -> Option<usize> {
    let map = std::env::var("DIALOG_GCD_COMPARE_STEP_BITS").ok()?;
    map.split(',').rev().find_map(|entry| {
        let (raw_step, raw_bits) = entry.trim().split_once(':')?;
        if raw_step.trim().parse::<usize>().ok()? != step {
            return None;
        }
        raw_bits.trim().parse::<usize>().ok()
    })
}

pub(crate) fn dialog_gcd_compare_bits_for_step(step: usize, active_width: usize) -> usize {
    if let Some(bits) = dialog_gcd_compare_step_bits(step) {
        return bits.clamp(1, active_width);
    }
    let global = dialog_gcd_compare_bits().min(active_width);
    if dialog_gcd_pa9024_compare_schedule_enabled() {
        let scheduled = (DIALOG_GCD_PA9024_COMPARE_SCHEDULE
            .get(step)
            .copied()
            .unwrap_or(global)
            + dialog_gcd_pa9024_compare_schedule_margin())
        .max(dialog_gcd_pa9024_compare_schedule_floor())
        .min(active_width);
        return scheduled.min(global).max(1);
    }
    global.max(1)
}

pub(crate) fn dialog_gcd_fused_branch_bits_enabled() -> bool {
    std::env::var("DIALOG_GCD_FUSED_BRANCH_BITS")
        .ok()
        .as_deref()
        == Some("1")
}

pub(crate) fn dialog_gcd_odd_u_lowbit_fastpath_enabled() -> bool {
    std::env::var("DIALOG_GCD_ODD_U_LOWBIT_FASTPATH")
        .ok()
        .as_deref()
        == Some("1")
}

pub(crate) fn dialog_gcd_k2_pair_compress_enabled() -> bool {
    dialog_gcd_k2_enabled()
        && !dialog_gcd_k5_clean_block_enabled()
        && std::env::var("DIALOG_GCD_K2_PAIR_COMPRESS")
            .ok()
            .as_deref()
            == Some("1")
}

pub(crate) fn dialog_gcd_sidecar_group_size() -> usize {
    if dialog_gcd_k5_clean_block_enabled() {
        5
    } else if dialog_gcd_k2_pair_compress_enabled() {
        2
    } else {
        DIALOG_GCD_HIGH_TAIL_ALIAS_GROUP_SIZE
    }
}

pub(crate) fn dialog_gcd_apply_fused_fold_enabled() -> bool {
    std::env::var("DIALOG_GCD_APPLY_FUSED_FOLD").ok().as_deref() == Some("1")
}

pub const DIALOG_FUSE_C_FORM_ENV: &str = "DIALOG_FUSE_C_FORM";
pub(crate) fn dialog_fuse_c_form_enabled() -> bool {
    std::env::var(DIALOG_FUSE_C_FORM_ENV).ok().as_deref() == Some("1")
}

pub const DIALOG_FUSE_X_RESTORE_ENV: &str = "DIALOG_FUSE_X_RESTORE";
pub(crate) fn dialog_fuse_x_restore_enabled() -> bool {
    std::env::var(DIALOG_FUSE_X_RESTORE_ENV).ok().as_deref() == Some("1")
}
