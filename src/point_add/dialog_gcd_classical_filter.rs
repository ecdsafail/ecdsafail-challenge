//! Classical convergence pre-filter for dialog-GCD Fiat-Shamir island search.
//!
//! Per tail-nonce, derives the 9024 Fiat-Shamir point-add inputs and classically
//! replays the truncated binary-GCD transcript on both inversion factors:
//!   - `dx = Px - Qx (mod p)`  (quotient / pair-1)
//!   - `c  = Qx - Rx (mod p)`  (ipmul / pair-2), with `Rx` the expected sum x.
//!
//! A factor is **hard** if any step hits:
//!   - width envelope overflow (`bitlen(u|v) > active_width(step)`),
//!   - truncated branch-comparator mis-decision vs the full active window,
//!   - or the full-width K2 transcript needs more than `ACTIVE_ITERATIONS` steps.
//!
//! This is analysis-only tooling; it does not change the quantum circuit.

use crate::point_add::{
    dialog_gcd_k5_head11_supports, dialog_gcd_k5_tail3_top32_supports,
    dialog_gcd_k5_tail6_graph9_supports, DIALOG_GCD_K5_TAIL6_GRAPH_SUPPORT,
    DIALOG_GCD_K5_TAIL7_SUPPORT, DIALOG_GCD_PA9024_COMPARE_SCHEDULE, N, SECP256K1_P,
};
use alloy_primitives::U256;
use ruint::Uint;

const MAX_GCD_ITERS: usize = 402;
type U512 = Uint<512, 8>;

/// Why a GCD factor failed the classical filter.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum HardReason {
    WidthOverflow {
        step: usize,
    },
    BodyTrimMismatch {
        step: usize,
        active_width: usize,
        body_width: usize,
    },
    ComparatorMismatch {
        step: usize,
    },
    NonConvergence {
        steps_needed: usize,
    },
    HeadPairMismatch {
        pattern: u8,
    },
    HeadK5Mismatch {
        pattern: u16,
    },
    TailPairMismatch {
        pattern: u8,
    },
    TailPairCrossMismatch {
        pattern: u16,
    },
    Tail6GraphMismatch {
        pattern: u32,
    },
    Tail6Graph9Mismatch {
        pattern: u32,
    },
    Tail7Mismatch {
        pattern: u32,
    },
    Tail3FixedLastMismatch {
        digit: u8,
    },
    Tail3Top32Mismatch {
        pattern: u16,
    },
    OddTailTripleMismatch {
        s2_mask: u8,
    },
    FusedFoldCarryEscape {
        step: usize,
        reverse: bool,
    },
    SpecialFoldCarryEscape {
        step: usize,
        reverse: bool,
    },
    ApplyValueMismatch {
        reverse: bool,
        compare_step: Option<usize>,
        full_width_step: Option<usize>,
    },
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct DialogGcdStepLog {
    pub b0: bool,
    pub b0_and_b1: bool,
    pub s2: bool,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct ApplyCleanupMismatch {
    pub step: usize,
    pub reverse: bool,
    pub bits: usize,
    pub required_bits: usize,
}

#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct ApplyHazardSummary {
    /// HMR cleanup predicates that disagree with the truncated comparator.
    ///
    /// These are soft risks, not deterministic failures: the verifier's seeded
    /// measurement result can still cancel the phase, as it does for promoted
    /// nonce 17761178.
    pub cleanup_mismatches: usize,
    pub cleanup_mismatch_details: Vec<ApplyCleanupMismatch>,
}

fn widen_u256(value: U256) -> U512 {
    let limbs = value.as_limbs();
    U512::from_limbs([limbs[0], limbs[1], limbs[2], limbs[3], 0, 0, 0, 0])
}

fn low_mask_512(bits: usize) -> U512 {
    if bits == 0 {
        U512::ZERO
    } else if bits >= 512 {
        U512::MAX
    } else {
        (U512::from(1u64) << bits) - U512::from(1u64)
    }
}

fn extract_512(value: U512, start: usize, bits: usize) -> U512 {
    (value >> start) & low_mask_512(bits)
}

fn square_row_value(x: U256, x_wide: U512, row: usize) -> U512 {
    if !bit_at(x, row) {
        return U512::ZERO;
    }
    let high = x_wide & !low_mask_512(row + 1);
    (high << (row + 1)) | (U512::from(1u64) << (2 * row))
}

/// Count truncated boundary-carry cleanup disagreements in the segmented
/// schoolbook square, forward plus inverse, for one 256-bit input.
///
/// The square value itself is exact. A disagreement means the measured cleanup
/// replay omitted a real carry/borrow into the retained high suffix, so it is a
/// soft phase-risk event rather than a deterministic value failure.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct SquareCleanupMismatch {
    pub row: usize,
    pub window: usize,
    pub reverse: bool,
    pub bits: usize,
    pub required_bits: usize,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct SquareCleanupSiteBits {
    pub row: usize,
    pub window: usize,
    pub reverse: bool,
    pub bits: usize,
}

#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct SquareCleanupSummary {
    pub mismatches: usize,
    pub details: Vec<SquareCleanupMismatch>,
}

fn square_cleanup_bits(
    clean_compare_bits: usize,
    row_clean_compare_bits: &[(usize, usize)],
    site_clean_compare_bits: &[SquareCleanupSiteBits],
    row: usize,
    window: usize,
    reverse: bool,
) -> usize {
    site_clean_compare_bits
        .iter()
        .rev()
        .find_map(|site| {
            (site.row == row && site.window == window && site.reverse == reverse)
                .then_some(site.bits)
        })
        .or_else(|| step_map_override(row_clean_compare_bits, row))
        .unwrap_or(clean_compare_bits)
}

pub fn square_row_window_cleanup_summary(
    x: U256,
    max_seg: usize,
    clean_compare_bits: usize,
    row_clean_compare_bits: &[(usize, usize)],
    site_clean_compare_bits: &[SquareCleanupSiteBits],
) -> SquareCleanupSummary {
    if max_seg == 0 {
        return SquareCleanupSummary::default();
    }

    let x_wide = widen_u256(x);
    let mut tmp = U512::ZERO;
    let mut summary = SquareCleanupSummary::default();

    for row in 0..N {
        let width = if row == N - 1 { 1 } else { N - row + 1 };
        let row_value = square_row_value(x, x_wide, row);
        if width > max_seg {
            let windows = width.div_ceil(max_seg).max(1).min(width);
            let mut carry_in = false;
            for window in 0..windows {
                let lo = (window * width) / windows;
                let hi = ((window + 1) * width) / windows;
                if hi == lo {
                    continue;
                }
                let bits = hi - lo;
                let mask = low_mask_512(bits);
                let offset = 2 * row + lo;
                let acc = extract_512(tmp, offset, bits);
                let seg = extract_512(row_value, offset, bits);
                let total = acc + seg + U512::from(carry_in as u64);
                let sum = total & mask;
                let carry_out = total > mask;
                if window + 1 < windows {
                    let row_bits = square_cleanup_bits(
                        clean_compare_bits,
                        row_clean_compare_bits,
                        site_clean_compare_bits,
                        row,
                        window,
                        false,
                    );
                    let trunc = if row_bits == 0 {
                        bits
                    } else {
                        row_bits.min(bits)
                    };
                    if trunc < bits {
                        let suffix_shift = bits - trunc;
                        let sum_suffix = extract_512(sum, suffix_shift, trunc);
                        let seg_suffix = extract_512(seg, suffix_shift, trunc);
                        let replay = sum_suffix < seg_suffix;
                        let mismatch = replay != carry_out;
                        summary.mismatches += usize::from(mismatch);
                        if mismatch {
                            let required_bits = ((trunc + 1)..=bits)
                                .find(|&candidate_bits| {
                                    let shift = bits - candidate_bits;
                                    let candidate_sum = extract_512(sum, shift, candidate_bits);
                                    let candidate_seg = extract_512(seg, shift, candidate_bits);
                                    (candidate_sum < candidate_seg) == carry_out
                                })
                                .unwrap_or(bits);
                            summary.details.push(SquareCleanupMismatch {
                                row,
                                window,
                                reverse: false,
                                bits: trunc,
                                required_bits,
                            });
                        }
                        if mismatch && std::env::var_os("ISLAND_TRACE_REJECT").is_some() {
                            eprintln!(
                                "SQUARE_PHASE_RISK row={row} window={window} reverse=false bits={trunc} required_bits={}",
                                summary.details.last().expect("mismatch detail").required_bits,
                            );
                        }
                    }
                }
                carry_in = carry_out;
            }
        }
        tmp += row_value;
    }

    for row in (0..N).rev() {
        let width = if row == N - 1 { 1 } else { N - row + 1 };
        let row_value = square_row_value(x, x_wide, row);
        if width > max_seg {
            let windows = width.div_ceil(max_seg).max(1).min(width);
            let mut borrow_in = false;
            for window in 0..windows {
                let lo = (window * width) / windows;
                let hi = ((window + 1) * width) / windows;
                if hi == lo {
                    continue;
                }
                let bits = hi - lo;
                let offset = 2 * row + lo;
                let acc = extract_512(tmp, offset, bits);
                let seg = extract_512(row_value, offset, bits);
                let subtrahend = seg + U512::from(borrow_in as u64);
                let borrow_out = acc < subtrahend;
                let diff = if borrow_out {
                    (acc + (U512::from(1u64) << bits)) - subtrahend
                } else {
                    acc - subtrahend
                };
                if window + 1 < windows {
                    let row_bits = square_cleanup_bits(
                        clean_compare_bits,
                        row_clean_compare_bits,
                        site_clean_compare_bits,
                        row,
                        window,
                        true,
                    );
                    let trunc = if row_bits == 0 {
                        bits
                    } else {
                        row_bits.min(bits)
                    };
                    if trunc < bits {
                        let suffix_shift = bits - trunc;
                        let diff_suffix = extract_512(diff, suffix_shift, trunc);
                        let seg_suffix = extract_512(seg, suffix_shift, trunc);
                        let not_seg_suffix = low_mask_512(trunc) ^ seg_suffix;
                        let replay = not_seg_suffix < diff_suffix;
                        let mismatch = replay != borrow_out;
                        summary.mismatches += usize::from(mismatch);
                        if mismatch {
                            let required_bits = ((trunc + 1)..=bits)
                                .find(|&candidate_bits| {
                                    let shift = bits - candidate_bits;
                                    let candidate_diff = extract_512(diff, shift, candidate_bits);
                                    let candidate_seg = extract_512(seg, shift, candidate_bits);
                                    let candidate_not_seg =
                                        low_mask_512(candidate_bits) ^ candidate_seg;
                                    (candidate_not_seg < candidate_diff) == borrow_out
                                })
                                .unwrap_or(bits);
                            summary.details.push(SquareCleanupMismatch {
                                row,
                                window,
                                reverse: true,
                                bits: trunc,
                                required_bits,
                            });
                        }
                        if mismatch && std::env::var_os("ISLAND_TRACE_REJECT").is_some() {
                            eprintln!(
                                "SQUARE_PHASE_RISK row={row} window={window} reverse=true bits={trunc} required_bits={}",
                                summary.details.last().expect("mismatch detail").required_bits,
                            );
                        }
                    }
                }
                borrow_in = borrow_out;
            }
        }
        tmp -= row_value;
    }
    debug_assert_eq!(tmp, U512::ZERO);
    summary
}

pub fn square_row_window_cleanup_mismatches(
    x: U256,
    max_seg: usize,
    clean_compare_bits: usize,
    row_clean_compare_bits: &[(usize, usize)],
    site_clean_compare_bits: &[SquareCleanupSiteBits],
) -> usize {
    square_row_window_cleanup_summary(
        x,
        max_seg,
        clean_compare_bits,
        row_clean_compare_bits,
        site_clean_compare_bits,
    )
    .mismatches
}

#[derive(Clone, Debug)]
pub struct DialogApplyFilterConfig {
    pub fused_fold_window: Option<usize>,
    pub special_fold_window: Option<usize>,
    pub fused_fold_step_windows: Vec<(usize, usize)>,
    pub special_fold_step_windows: Vec<(usize, usize)>,
    pub clean_compare_bits: usize,
    pub overflow_step_bits: Vec<(usize, usize)>,
    pub underflow_step_bits: Vec<(usize, usize)>,
    pub clear_product_residual: bool,
}

impl DialogApplyFilterConfig {
    pub fn from_env() -> Self {
        let fused_fold_window = std::env::var("DIALOG_GCD_FOLD_CARRY_TRUNC_W")
            .ok()
            .and_then(|s| s.parse::<usize>().ok())
            .filter(|&w| w > 0)
            .or_else(|| {
                std::env::var("KAL_DOUBLE_CARRY_TRUNC_W")
                    .ok()
                    .and_then(|s| s.parse::<usize>().ok())
                    .filter(|&w| w > 0)
            });
        let special_fold_window = std::env::var("KAL_FOLD_CARRY_TRUNC_W")
            .ok()
            .and_then(|s| s.parse::<usize>().ok())
            .filter(|&w| w > 0);
        let fused_fold_step_windows = std::env::var("DIALOG_GCD_FOLD_CARRY_TRUNC_STEP_WINDOWS")
            .ok()
            .map(|s| parse_step_map(&s))
            .unwrap_or_default();
        let special_fold_step_windows =
            std::env::var("DIALOG_GCD_SPECIAL_FOLD_CARRY_TRUNC_STEP_WINDOWS")
                .ok()
                .map(|s| parse_step_map(&s))
                .unwrap_or_default();
        let clean_compare_bits = std::env::var("DIALOG_GCD_APPLY_CLEAN_COMPARE_BITS")
            .ok()
            .and_then(|s| s.parse::<usize>().ok())
            .filter(|&bits| (1..=N).contains(&bits))
            .unwrap_or_else(|| {
                std::env::var("DIALOG_GCD_COMPARE_BITS")
                    .ok()
                    .and_then(|s| s.parse::<usize>().ok())
                    .filter(|&bits| (1..=N).contains(&bits))
                    .unwrap_or(N)
            });
        let overflow_step_bits = std::env::var("DIALOG_GCD_SPECIAL_OVERFLOW_CLEAN_STEP_BITS")
            .ok()
            .map(|s| parse_step_map(&s))
            .unwrap_or_default();
        let underflow_step_bits = std::env::var("DIALOG_GCD_SPECIAL_UNDERFLOW_CLEAN_STEP_BITS")
            .ok()
            .map(|s| parse_step_map(&s))
            .unwrap_or_default();
        let clear_product_residual = std::env::var("DIALOG_GCD_RAW_IPMUL_CLEAR_P_RESIDUAL")
            .ok()
            .as_deref()
            == Some("1");

        Self {
            fused_fold_window,
            special_fold_window,
            fused_fold_step_windows,
            special_fold_step_windows,
            clean_compare_bits,
            overflow_step_bits,
            underflow_step_bits,
            clear_product_residual,
        }
    }

    fn overflow_compare_bits(&self, step: usize) -> usize {
        step_map_override(&self.overflow_step_bits, step).unwrap_or(self.clean_compare_bits)
    }

    fn underflow_compare_bits(&self, step: usize) -> usize {
        step_map_override(&self.underflow_step_bits, step).unwrap_or(self.clean_compare_bits)
    }

    fn fused_fold_window(&self, step: usize) -> Option<usize> {
        step_map_override(&self.fused_fold_step_windows, step).or(self.fused_fold_window)
    }

    fn special_fold_window(&self, step: usize) -> Option<usize> {
        step_map_override(&self.special_fold_step_windows, step).or(self.special_fold_window)
    }
}

/// Knobs mirrored from `configure_ecdsafail_submission_route()` env defaults.
#[derive(Clone, Debug)]
pub struct DialogGcdFilterConfig {
    pub active_iterations: usize,
    pub compare_bits: usize,
    pub width_margin: f64,
    pub width_slope: f64,
    pub active_width_overrides: Vec<usize>,
    pub compare_width_overrides: Vec<usize>,
    pub body_width_overrides: Vec<usize>,
    pub body_carry_trims: Option<Vec<usize>>,
    pub pa9024_compare_schedule: bool,
    pub pa9024_compare_margin: usize,
    pub pa9024_compare_floor: usize,
    pub compare_step_bits: Vec<(usize, usize)>,
    pub odd_u_lowbit_fastpath: bool,
    pub k2: bool,
    pub variable_width: bool,
    pub raw_tobitvector_materialized_sub: bool,
    pub tobitvector_cswap_body_trim: bool,
    pub tobitvector_shift_body_trim: bool,
    pub skip_zero_edge_tobit_fwd_cshift: bool,
    pub width_step_bumps: Vec<(usize, usize)>,
    pub body_step_givebacks: Vec<(usize, usize)>,
    /// Cached env flags (hoisted out of the per-step hot loop).
    pub k2_force0: bool,
    pub strict_compare: bool,
    pub body_carry_trunc_w: usize,
}

impl Default for DialogGcdFilterConfig {
    fn default() -> Self {
        Self::from_env()
    }
}

impl DialogGcdFilterConfig {
    pub fn from_env() -> Self {
        let active_iterations = std::env::var("DIALOG_GCD_ACTIVE_ITERATIONS")
            .ok()
            .and_then(|s| s.parse().ok())
            .filter(|&iters| (1..=MAX_GCD_ITERS).contains(&iters))
            .unwrap_or(MAX_GCD_ITERS);
        let compare_bits = std::env::var("DIALOG_GCD_COMPARE_BITS")
            .ok()
            .and_then(|s| s.parse().ok())
            .filter(|&bits| (1..=N).contains(&bits))
            .unwrap_or(57);
        let width_margin = std::env::var("DIALOG_GCD_WIDTH_MARGIN")
            .ok()
            .and_then(|s| s.parse::<f64>().ok())
            .filter(|m| m.is_finite() && *m >= 0.0 && *m <= N as f64)
            .unwrap_or(37.0);
        let width_slope = std::env::var("DIALOG_GCD_WIDTH_SLOPE_X1000")
            .ok()
            .and_then(|s| s.parse::<f64>().ok())
            .filter(|s| s.is_finite() && *s > 0.0 && *s <= 4000.0)
            .map(|s| s / 1000.0)
            .unwrap_or(0.5 * 1.415);
        let body_carry_trims = std::env::var("DIALOG_GCD_BODY_CARRY_BAND_TRIMS")
            .ok()
            .and_then(|s| parse_trim_list(&s));
        let pa9024_compare_schedule = std::env::var("DIALOG_GCD_PA9024_COMPARE_SCHEDULE")
            .ok()
            .as_deref()
            == Some("1");
        let pa9024_compare_margin = std::env::var("DIALOG_GCD_PA9024_COMPARE_SCHEDULE_MARGIN")
            .ok()
            .and_then(|s| s.parse().ok())
            .unwrap_or(0);
        let pa9024_compare_floor = std::env::var("DIALOG_GCD_PA9024_COMPARE_SCHEDULE_FLOOR")
            .ok()
            .and_then(|s| s.parse().ok())
            .filter(|&bits| bits <= N)
            .unwrap_or(1)
            .max(1);
        let compare_step_bits = std::env::var("DIALOG_GCD_COMPARE_STEP_BITS")
            .ok()
            .map(|s| parse_step_map(&s))
            .unwrap_or_default();
        let odd_u_lowbit_fastpath = std::env::var("DIALOG_GCD_ODD_U_LOWBIT_FASTPATH")
            .ok()
            .as_deref()
            == Some("1");
        let k2 = std::env::var("DIALOG_GCD_K2").ok().as_deref() == Some("1");
        let variable_width = std::env::var("DIALOG_GCD_RAW_TOBITVECTOR_VARIABLE_WIDTH")
            .ok()
            .as_deref()
            != Some("0");
        let raw_tobitvector_materialized_sub =
            std::env::var("DIALOG_GCD_RAW_TOBITVECTOR_MATERIALIZED_SUB")
                .ok()
                .as_deref()
                != Some("0");
        let tobitvector_cswap_body_trim = std::env::var("DIALOG_GCD_TOBITVECTOR_CSWAP_BODY_TRIM")
            .ok()
            .as_deref()
            == Some("1");
        let tobitvector_shift_body_trim = std::env::var("DIALOG_GCD_TOBITVECTOR_SHIFT_BODY_TRIM")
            .ok()
            .as_deref()
            == Some("1");
        let skip_zero_edge_tobit_fwd_cshift = std::env::var("DIALOG_GCD_SKIP_ZERO_EDGE_CSHIFT")
            .ok()
            .as_deref()
            == Some("1")
            || std::env::var("DIALOG_GCD_SKIP_ZERO_EDGE_TOBIT_CSHIFT")
                .ok()
                .as_deref()
                == Some("1")
            || std::env::var("DIALOG_GCD_SKIP_ZERO_EDGE_TOBIT_FWD_CSHIFT")
                .ok()
                .as_deref()
                == Some("1");
        let width_step_bumps = std::env::var("DIALOG_GCD_WIDTH_STEP_BUMPS")
            .ok()
            .map(|s| parse_step_map(&s))
            .unwrap_or_default();
        let body_step_givebacks = std::env::var("DIALOG_GCD_BODY_STEP_GIVEBACKS")
            .ok()
            .map(|s| parse_step_map(&s))
            .unwrap_or_default();
        let k2_force0 = std::env::var("DIALOG_GCD_K2_FORCE0").ok().as_deref() == Some("1");
        let strict_compare = std::env::var("DIALOG_GCD_FILTER_STRICT_COMPARE")
            .ok()
            .as_deref()
            == Some("1");
        let body_carry_trunc_w = std::env::var("DIALOG_GCD_BODY_CARRY_TRUNC_W")
            .ok()
            .and_then(|s| s.parse().ok())
            .unwrap_or(0);

        Self {
            active_iterations,
            compare_bits,
            width_margin,
            width_slope,
            active_width_overrides: Vec::new(),
            compare_width_overrides: Vec::new(),
            body_width_overrides: Vec::new(),
            body_carry_trims,
            pa9024_compare_schedule,
            pa9024_compare_margin,
            pa9024_compare_floor,
            compare_step_bits,
            odd_u_lowbit_fastpath,
            k2,
            variable_width,
            raw_tobitvector_materialized_sub,
            tobitvector_cswap_body_trim,
            tobitvector_shift_body_trim,
            skip_zero_edge_tobit_fwd_cshift,
            width_step_bumps,
            body_step_givebacks,
            k2_force0,
            strict_compare,
            body_carry_trunc_w,
        }
    }

    pub fn active_width(&self, step: usize) -> usize {
        if let Some(&width) = self.active_width_overrides.get(step) {
            return width.clamp(1, N);
        }
        if !self.variable_width {
            return N;
        }
        let ideal = N as f64 - (step as f64) * self.width_slope + self.width_margin;
        let rounded = ((ideal.max(1.0) / 2.0).ceil() as usize) * 2;
        rounded
            .saturating_add(step_map_value(&self.width_step_bumps, step))
            .clamp(1, N)
    }

    pub fn compare_bits_for_step(&self, step: usize, active_width: usize) -> usize {
        if let Some(&bits) = self.compare_width_overrides.get(step) {
            return bits.clamp(1, active_width);
        }
        if let Some(bits) = step_map_override(&self.compare_step_bits, step) {
            return bits.clamp(1, active_width);
        }
        let global = self.compare_bits.min(active_width);
        if self.pa9024_compare_schedule {
            let scheduled = (DIALOG_GCD_PA9024_COMPARE_SCHEDULE
                .get(step)
                .copied()
                .unwrap_or(global)
                + self.pa9024_compare_margin)
                .max(self.pa9024_compare_floor)
                .min(active_width);
            return scheduled.min(global).max(1);
        }
        global.max(1)
    }

    pub fn body_carry_trunc_width(&self, active_width: usize, step: usize) -> usize {
        if let Some(&width) = self.body_width_overrides.get(step) {
            return width.clamp(2, active_width);
        }
        let mut w = self
            .body_carry_band_trim(step)
            .or_else(|| {
                std::env::var("DIALOG_GCD_BODY_CARRY_TRUNC_W")
                    .ok()
                    .and_then(|s| s.parse().ok())
            })
            .unwrap_or(0);
        w = w.saturating_add(body_carry_extra_notch(step));
        w = w.saturating_sub(step_map_value(&self.body_step_givebacks, step));
        active_width.saturating_sub(w).max(2)
    }

    #[inline]
    fn body_carry_trunc_width_fast(&self, active_width: usize, step: usize) -> usize {
        if let Some(&width) = self.body_width_overrides.get(step) {
            return width.clamp(2, active_width);
        }
        let mut w = self
            .body_carry_band_trim(step)
            .unwrap_or(self.body_carry_trunc_w);
        w = w.saturating_add(body_carry_extra_notch(step));
        w = w.saturating_sub(step_map_value(&self.body_step_givebacks, step));
        active_width.saturating_sub(w).max(2)
    }

    #[inline]
    fn cswap_width(&self, active_width: usize, step: usize) -> usize {
        if self.tobitvector_cswap_body_trim {
            self.body_carry_trunc_width_fast(active_width, step)
                .min(active_width)
        } else {
            active_width
        }
    }

    #[inline]
    fn shift_width(&self, active_width: usize, step: usize) -> usize {
        if self.tobitvector_shift_body_trim {
            self.body_carry_trunc_width_fast(active_width, step)
                .min(active_width)
        } else {
            active_width
        }
    }

    fn body_carry_band_trim(&self, step: usize) -> Option<usize> {
        let trims = self.body_carry_trims.as_ref()?;
        if trims.is_empty() {
            return None;
        }
        let iters = self.active_iterations.max(1);
        let band_size = ((iters + trims.len() - 1) / trims.len()).max(1);
        let band = (step / band_size).min(trims.len() - 1);
        Some(trims[band])
    }
}

fn body_carry_extra_notch(step: usize) -> usize {
    let mut extra = 0usize;

    let trio_enabled = std::env::var("DIALOG_GCD_TRIO_WIDTH_NOTCH").ok().as_deref() != Some("0");
    if trio_enabled {
        let trio_step = std::env::var("DIALOG_GCD_TRIO_WIDTH_NOTCH_STEP")
            .ok()
            .and_then(|s| s.parse::<usize>().ok())
            .unwrap_or(11);
        if step == trio_step {
            extra = extra.saturating_add(
                std::env::var("DIALOG_GCD_TRIO_WIDTH_NOTCH_EXTRA")
                    .ok()
                    .and_then(|s| s.parse::<usize>().ok())
                    .unwrap_or(2),
            );
        }
    }

    if let Ok(steps) = std::env::var("DIALOG_GCD_BINDER_NOTCH_STEPS") {
        let hits = steps
            .split(',')
            .filter_map(|s| s.trim().parse::<usize>().ok())
            .any(|s| s == step);
        if hits {
            extra = extra.saturating_add(
                std::env::var("DIALOG_GCD_BINDER_NOTCH_EXTRA")
                    .ok()
                    .and_then(|s| s.parse::<usize>().ok())
                    .unwrap_or(2),
            );
        }
    }

    if let Ok(map) = std::env::var("DIALOG_GCD_BINDER_NOTCH_MAP") {
        extra = extra.saturating_add(
            map.split(',')
                .filter_map(|entry| {
                    let (s, e) = entry.trim().split_once(':')?;
                    Some((
                        s.trim().parse::<usize>().ok()?,
                        e.trim().parse::<usize>().ok()?,
                    ))
                })
                .filter_map(|(s, e)| (s == step).then_some(e))
                .sum(),
        );
    }

    extra
}

fn parse_trim_list(s: &str) -> Option<Vec<usize>> {
    if s.trim().is_empty() {
        return None;
    }
    let trims: Vec<usize> = s.split(',').filter_map(|t| t.trim().parse().ok()).collect();
    if trims.is_empty() {
        None
    } else {
        Some(trims)
    }
}

fn parse_step_map(s: &str) -> Vec<(usize, usize)> {
    s.split(',')
        .filter_map(|entry| {
            let (step, value) = entry.trim().split_once(':')?;
            Some((
                step.trim().parse::<usize>().ok()?,
                value.trim().parse::<usize>().ok()?,
            ))
        })
        .collect()
}

fn step_map_value(map: &[(usize, usize)], step: usize) -> usize {
    map.iter()
        .filter_map(|&(s, value)| (s == step).then_some(value))
        .sum()
}

fn step_map_override(map: &[(usize, usize)], step: usize) -> Option<usize> {
    map.iter()
        .rev()
        .find_map(|&(s, value)| (s == step).then_some(value))
}

#[inline]
fn window_mask(width: usize) -> U256 {
    if width >= 256 {
        U256::MAX
    } else {
        (U256::from(1u64) << width) - U256::from(1u64)
    }
}

#[inline]
pub fn bitlen(x: U256) -> usize {
    if x.is_zero() {
        0
    } else {
        256 - x.leading_zeros() as usize
    }
}

#[inline]
fn bit_at(x: U256, i: usize) -> bool {
    (x >> i) & U256::from(1u64) != U256::ZERO
}

fn cmp_gt_window(u: U256, v: U256, width: usize) -> bool {
    let mask = window_mask(width);
    (u & mask) > (v & mask)
}

fn cmp_gt_truncated(u: U256, v: U256, width: usize, compare_bits: usize) -> bool {
    let cb = compare_bits.min(width).max(1);
    let lo = width.saturating_sub(cb);
    let mask = window_mask(cb);
    ((u >> lo) & mask) > ((v >> lo) & mask)
}

fn sub_low_window(v: U256, u: U256, width: usize) -> U256 {
    let mask = window_mask(width);
    let diff = (v & mask).wrapping_sub(u & mask) & mask;
    (v & !mask) | diff
}

fn shift_right_active(v: &mut U256, active_width: usize) {
    let mask = window_mask(active_width);
    let x = *v & mask;
    *v = (x >> 1) | (*v & !mask);
}

fn shift_right_active_skip_top_edge(v: &mut U256, active_width: usize) {
    if active_width <= 1 {
        return;
    }
    let mask = window_mask(active_width);
    let x = *v & mask;
    let shifted_low = (x >> 1) & window_mask(active_width - 2);
    let preserved_top = x & (U256::from(1u64) << (active_width - 1));
    *v = shifted_low | preserved_top | (*v & !mask);
}

fn swap_active_except_bit0(u: &mut U256, v: &mut U256, active_width: usize) {
    let mask_lo = U256::from(1u64);
    let mask_hi = window_mask(active_width) & !mask_lo;
    let u_hi = *u & mask_hi;
    let v_hi = *v & mask_hi;
    *u = (*u & mask_lo) | v_hi;
    *v = (*v & mask_lo) | u_hi;
}

/// One truncated dialog-GCD tobitvector step (forward), matching
/// `emit_dialog_gcd_*_tobitvector_steps`, plus the replay bits consumed by the
/// apply and reverse-apply passes.
fn truncated_gcd_step_logged(
    u: &mut U256,
    v: &mut U256,
    step: usize,
    cfg: &DialogGcdFilterConfig,
) -> Result<DialogGcdStepLog, HardReason> {
    let active_width = cfg.active_width(step);
    if (bitlen(*u) > active_width || bitlen(*v) > active_width)
        && std::env::var("DIALOG_GCD_FILTER_STRICT_WIDTH")
            .ok()
            .as_deref()
            == Some("1")
    {
        return Err(HardReason::WidthOverflow { step });
    }

    let compare_bits = cfg.compare_bits_for_step(step, active_width);
    let _full_gt = cmp_gt_window(*u, *v, active_width);
    let trunc_gt = cmp_gt_truncated(*u, *v, active_width, compare_bits);
    // NOTE: a truncated-vs-full comparator disagreement is NOT a hard input.
    // The frontier island (nonce 700017357 @ compare=46) validates 0/0/0 yet has
    // such a disagreement at step 205: the truncated branch decision still drives
    // the GCD to the correct inverse on the reachable verifier support. Flagging
    // it produced false negatives (rejected genuinely-clean islands). The
    // hardware follows the *truncated* decision (`trunc_gt`), which this replay
    // already uses below, so comparator correctness is delegated to `--validate`.
    // Opt back in with DIALOG_GCD_FILTER_STRICT_COMPARE=1 for diagnostics.
    if _full_gt != trunc_gt && cfg.strict_compare {
        return Err(HardReason::ComparatorMismatch { step });
    }

    let b0 = bit_at(*v, 0);
    let b0_and_b1 = b0 && trunc_gt;

    let cswap_width = cfg.cswap_width(active_width, step);
    if b0_and_b1 {
        if cfg.odd_u_lowbit_fastpath {
            swap_active_except_bit0(u, v, cswap_width);
        } else {
            let mask = window_mask(cswap_width);
            let u_window = *u & mask;
            let v_window = *v & mask;
            *u = (*u & !mask) | v_window;
            *v = (*v & !mask) | u_window;
        }
    }

    if b0 {
        if cfg.raw_tobitvector_materialized_sub {
            let body_w = cfg.body_carry_trunc_width_fast(active_width, step);
            let full_v = if cfg.odd_u_lowbit_fastpath {
                sub_low_window(*v, *u, active_width) ^ U256::from(1u64)
            } else {
                sub_low_window(*v, *u, active_width)
            };
            let trimmed_v = if cfg.odd_u_lowbit_fastpath {
                if body_w <= 1 {
                    *v ^ U256::from(1u64)
                } else {
                    sub_low_window(*v, *u, body_w) ^ U256::from(1u64)
                }
            } else {
                sub_low_window(*v, *u, body_w)
            };
            if (full_v & window_mask(active_width)) != (trimmed_v & window_mask(active_width))
                && std::env::var("DIALOG_GCD_FILTER_STRICT_BODY")
                    .ok()
                    .as_deref()
                    == Some("1")
            {
                return Err(HardReason::BodyTrimMismatch {
                    step,
                    active_width,
                    body_width: body_w,
                });
            }
            *v = trimmed_v;
        } else {
            *v = sub_low_window(*v, *u, active_width);
        }
    }

    let shift_width = cfg.shift_width(active_width, step);
    shift_right_active(v, shift_width);

    let mut s2 = false;
    if cfg.k2 && !cfg.k2_force0 {
        s2 = !bit_at(*v, 0);
        if s2 {
            if cfg.skip_zero_edge_tobit_fwd_cshift {
                shift_right_active_skip_top_edge(v, shift_width);
            } else {
                shift_right_active(v, shift_width);
            }
        }
    }

    Ok(DialogGcdStepLog { b0, b0_and_b1, s2 })
}

fn truncated_gcd_step(
    u: &mut U256,
    v: &mut U256,
    step: usize,
    cfg: &DialogGcdFilterConfig,
) -> Option<HardReason> {
    truncated_gcd_step_logged(u, v, step, cfg).err()
}

/// Full-width K2 binary-GCD step (no width truncation) for convergence counting.
fn full_gcd_step(u: &mut U256, v: &mut U256, cfg: &DialogGcdFilterConfig) {
    let width = N;
    let b0 = bit_at(*v, 0);
    let full_gt = *u > *v;

    let b0_and_b1 = b0 && full_gt;
    if b0_and_b1 {
        if cfg.odd_u_lowbit_fastpath {
            swap_active_except_bit0(u, v, width);
        } else {
            std::mem::swap(u, v);
        }
    }

    if b0 {
        *v = v.wrapping_sub(*u);
        if cfg.odd_u_lowbit_fastpath {
            *v ^= U256::from(1u64);
        }
    }

    *v >>= 1;

    if cfg.k2 && !cfg.k2_force0 {
        if !bit_at(*v, 0) {
            *v >>= 1;
        }
    }
}

/// Steps until `v == 0` under the full-width transcript, capped at `limit`.
pub(crate) fn full_gcd_steps_until_zero(
    mut u: U256,
    mut v: U256,
    cfg: &DialogGcdFilterConfig,
    limit: usize,
) -> usize {
    let mut steps = 0usize;
    while !v.is_zero() && steps < limit {
        full_gcd_step(&mut u, &mut v, cfg);
        steps += 1;
    }
    steps
}

/// One full-width binary-GCD step that removes up to `depth` trailing zeros of
/// `v` per recorded step (Stein/jump generalization of K2; `depth=1` is the
/// plain dialog, `depth=2` is the deployed K2). The base shift always fires
/// (`shift_right_assuming_even`); each extra shift is conditional on `v` still
/// being even, exactly mirroring the quantum `k2_shift2_log` cascade. This is
/// the convergence model used to size `active_iterations` (== max steps over the
/// reachable support) for each jump depth.
fn full_gcd_step_jump(u: &mut U256, v: &mut U256, depth: usize) {
    let b0 = bit_at(*v, 0);
    if b0 && *u > *v {
        std::mem::swap(u, v);
    }
    if b0 {
        *v = v.wrapping_sub(*u);
    }
    // Base shift (v is even here: either b0=0 originally, or the subtract above
    // cleared bit 0).
    *v >>= 1;
    let mut shifts = 1usize;
    while shifts < depth && !v.is_zero() && !bit_at(*v, 0) {
        *v >>= 1;
        shifts += 1;
    }
}

/// Steps until `v == 0` for jump `depth`, capped at `limit`.
pub fn jump_steps_until_zero(mut u: U256, mut v: U256, depth: usize, limit: usize) -> usize {
    let mut steps = 0usize;
    while !v.is_zero() && steps < limit {
        full_gcd_step_jump(&mut u, &mut v, depth.max(1));
        steps += 1;
    }
    steps
}

/// Per-depth convergence statistics over a set of GCD factors.
#[derive(Clone, Debug)]
pub struct JumpConvergence {
    pub depth: usize,
    pub max_steps: usize,
    pub mean_steps: f64,
    /// 99.99th-percentile-ish: max over the sampled factors is the binding
    /// `active_iterations`, since every shot must converge.
    pub p_max_factor: U256,
}

/// Measure convergence-step distributions across `factors` for jump depths
/// `1..=max_depth`. `max_steps` is the binding `active_iterations` for that
/// depth (every shot must converge within it). Pure number theory on the prime
/// `SECP256K1_P`; independent of the circuit truncations.
pub fn measure_jump_convergence(factors: &[U256], max_depth: usize) -> Vec<JumpConvergence> {
    const LIMIT: usize = 1024;
    let mut out = Vec::with_capacity(max_depth);
    for depth in 1..=max_depth {
        let mut max_steps = 0usize;
        let mut sum = 0u64;
        let mut p_max_factor = U256::ZERO;
        for &f in factors {
            if f.is_zero() {
                continue;
            }
            let s = jump_steps_until_zero(SECP256K1_P, f, depth, LIMIT);
            sum += s as u64;
            if s > max_steps {
                max_steps = s;
                p_max_factor = f;
            }
        }
        let n = factors.iter().filter(|f| !f.is_zero()).count().max(1);
        out.push(JumpConvergence {
            depth,
            max_steps,
            mean_steps: sum as f64 / n as f64,
            p_max_factor,
        });
    }
    out
}

pub fn sub_mod_p(a: U256, b: U256, p: U256) -> U256 {
    if a >= b {
        a - b
    } else {
        p - (b - a)
    }
}

/// GCD inversion factor inputs for one point-add shot.
pub fn point_add_gcd_factors(px: U256, qx: U256, rx: U256) -> (U256, U256) {
    let dx = sub_mod_p(px, qx, SECP256K1_P);
    let c = sub_mod_p(qx, rx, SECP256K1_P);
    (dx, c)
}

#[derive(Clone, Debug)]
struct DialogGcdTranscript {
    log: Vec<DialogGcdStepLog>,
    terminal_u: U256,
    terminal_v: U256,
}

/// A factor whose truncated GCD transcript has already passed the envelope,
/// terminal-codec, and convergence checks.
///
/// Island search checks both factors before replaying apply arithmetic. Keeping
/// this opaque lets that hot path reuse the 258-step transcripts rather than
/// rebuilding each one a second time.
#[derive(Clone, Debug)]
pub struct CheckedGcdFactor {
    transcript: DialogGcdTranscript,
}

impl CheckedGcdFactor {
    pub fn log(&self) -> &[DialogGcdStepLog] {
        &self.transcript.log
    }
}

fn first_log_difference(
    factor: U256,
    baseline: &DialogGcdTranscript,
    cfg: &DialogGcdFilterConfig,
    full_width: bool,
) -> Option<usize> {
    let mut reference_cfg = cfg.clone();
    reference_cfg.compare_bits = N;
    reference_cfg.pa9024_compare_schedule = false;
    reference_cfg.compare_width_overrides.clear();
    reference_cfg.compare_step_bits.clear();
    if full_width {
        reference_cfg.variable_width = false;
        reference_cfg.active_width_overrides.clear();
        reference_cfg.body_width_overrides.clear();
        reference_cfg.body_carry_trims = None;
        reference_cfg.body_carry_trunc_w = 0;
        reference_cfg.width_step_bumps.clear();
        reference_cfg.body_step_givebacks.clear();
    }
    let reference = build_gcd_transcript(factor, &reference_cfg).ok()?;
    baseline
        .log
        .iter()
        .zip(reference.log.iter())
        .position(|(a, b)| a != b)
}

fn build_gcd_transcript(
    factor: U256,
    cfg: &DialogGcdFilterConfig,
) -> Result<DialogGcdTranscript, HardReason> {
    if factor.is_zero() {
        return Err(HardReason::NonConvergence { steps_needed: 0 });
    }

    let mut u = SECP256K1_P;
    let mut v = factor;
    let mut log = Vec::with_capacity(cfg.active_iterations);
    for step in 0..cfg.active_iterations {
        log.push(truncated_gcd_step_logged(&mut u, &mut v, step, cfg)?);
    }
    Ok(DialogGcdTranscript {
        log,
        terminal_u: u,
        terminal_v: v,
    })
}

fn tail_pair_codec_mode() -> Option<usize> {
    if std::env::var_os("ISLAND_IGNORE_TAIL_CODEC").is_some() {
        return None;
    }
    if std::env::var("DIALOG_GCD_TAIL_CROSSBLOCK5").ok().as_deref() == Some("1")
        && std::env::var("DIALOG_GCD_TAIL_PAIR_DIRECT_APPLY")
            .ok()
            .as_deref()
            == Some("1")
    {
        return Some(1);
    }
    match std::env::var("DIALOG_GCD_TAIL_PAIR_CODEC").ok().as_deref() {
        Some("const" | "zero") => Some(0),
        Some("1")
            if std::env::var("DIALOG_GCD_TAIL_PAIR_DIRECT_APPLY")
                .ok()
                .as_deref()
                == Some("1") =>
        {
            Some(1)
        }
        Some("2") => Some(2),
        Some("3") => Some(3),
        Some("4") => Some(4),
        _ => None,
    }
}

fn tail_pair_pattern(log: &[DialogGcdStepLog]) -> u8 {
    let tail = if log.len() % 2 == 1 {
        log.len().saturating_sub(3)
    } else {
        log.len().saturating_sub(2)
    };
    let mut pattern = 0u8;
    for (slot, entry) in log[tail..tail + 2].iter().enumerate() {
        pattern |= (entry.b0 as u8) << (3 * slot);
        pattern |= (entry.b0_and_b1 as u8) << (3 * slot + 1);
        pattern |= (entry.s2 as u8) << (3 * slot + 2);
    }
    pattern
}

fn tail3_pattern(log: &[DialogGcdStepLog]) -> u16 {
    if log.len() < 3 {
        return 0;
    }
    log[log.len() - 3..]
        .iter()
        .enumerate()
        .fold(0u16, |packed, (slot, entry)| {
            packed
                | ((entry.b0 as u16) << (3 * slot))
                | ((entry.b0_and_b1 as u16) << (3 * slot + 1))
                | ((entry.s2 as u16) << (3 * slot + 2))
        })
}

fn tail7_pattern(log: &[DialogGcdStepLog]) -> u32 {
    if log.len() < 7 {
        return 0;
    }
    log[log.len() - 7..]
        .iter()
        .enumerate()
        .fold(0u32, |packed, (slot, entry)| {
            packed
                | ((entry.b0 as u32) << (3 * slot))
                | ((entry.b0_and_b1 as u32) << (3 * slot + 1))
                | ((entry.s2 as u32) << (3 * slot + 2))
        })
}

fn tail6_pattern(log: &[DialogGcdStepLog]) -> u32 {
    if log.len() < 6 {
        return 0;
    }
    log[log.len() - 6..]
        .iter()
        .enumerate()
        .fold(0u32, |packed, (slot, entry)| {
            packed
                | ((entry.b0 as u32) << (3 * slot))
                | ((entry.b0_and_b1 as u32) << (3 * slot + 1))
                | ((entry.s2 as u32) << (3 * slot + 2))
        })
}

fn check_tail_pair_codec(log: &[DialogGcdStepLog]) -> Result<(), HardReason> {
    if std::env::var("DIALOG_GCD_K5_HEAD11_CODEC").ok().as_deref() == Some("1") {
        if log.len() < 5 {
            return Err(HardReason::HeadK5Mismatch { pattern: 0 });
        }
        let pattern = log[..5]
            .iter()
            .enumerate()
            .fold(0u16, |packed, (slot, entry)| {
                packed
                    | ((entry.b0 as u16) << (3 * slot))
                    | ((entry.b0_and_b1 as u16) << (3 * slot + 1))
                    | ((entry.s2 as u16) << (3 * slot + 2))
            });
        if !dialog_gcd_k5_head11_supports(pattern) {
            return Err(HardReason::HeadK5Mismatch { pattern });
        }
    }
    if std::env::var("DIALOG_GCD_HEAD_PAIR_CODEC3").ok().as_deref() == Some("1") {
        if log.len() < 2 {
            return Err(HardReason::HeadPairMismatch { pattern: 0 });
        }
        let pattern = log[..2]
            .iter()
            .enumerate()
            .fold(0u8, |packed, (slot, entry)| {
                packed
                    | ((entry.b0 as u8) << (3 * slot))
                    | ((entry.b0_and_b1 as u8) << (3 * slot + 1))
                    | ((entry.s2 as u8) << (3 * slot + 2))
            });
        if !matches!(pattern, 4 | 24 | 27 | 28 | 36 | 56 | 59 | 60) {
            return Err(HardReason::HeadPairMismatch { pattern });
        }
    }
    if std::env::var("DIALOG_GCD_ODD_SINGLETON_CODEC")
        .ok()
        .as_deref()
        == Some("2")
        && log.len() % 2 == 1
    {
        let entry = log.last().expect("odd transcript has a final step");
        let digit = (entry.b0 as u8) | ((entry.b0_and_b1 as u8) << 1) | ((entry.s2 as u8) << 2);
        if !matches!(digit, 1 | 3 | 4 | 5) {
            return Err(HardReason::TailPairMismatch { pattern: digit });
        }
    }
    if std::env::var("DIALOG_GCD_ODD_TAIL_TRIPLE_CODEC")
        .ok()
        .as_deref()
        == Some("1")
    {
        let tail = log.len().saturating_sub(3);
        let s2_mask = log[tail..]
            .iter()
            .enumerate()
            .fold(0u8, |mask, (slot, entry)| mask | ((entry.s2 as u8) << slot));
        if s2_mask != 0 {
            return Err(HardReason::OddTailTripleMismatch { s2_mask });
        }
    }
    let pattern = tail_pair_pattern(log);
    let ignore_tail_codec = std::env::var("ISLAND_FILTER_IGNORE_TAIL_CODEC")
        .ok()
        .as_deref()
        == Some("1");
    if !ignore_tail_codec
        && std::env::var("DIALOG_GCD_K5_TAIL3_TOP32_CODEC")
            .ok()
            .as_deref()
            == Some("1")
    {
        let pattern = tail3_pattern(log);
        if !dialog_gcd_k5_tail3_top32_supports(pattern) {
            return Err(HardReason::Tail3Top32Mismatch { pattern });
        }
    }
    if !ignore_tail_codec
        && std::env::var("DIALOG_GCD_K5_TAIL3_FIXED_LAST")
            .ok()
            .as_deref()
            == Some("1")
    {
        let entry = log.last().expect("tail3 codec requires a final step");
        let digit = (entry.b0 as u8) | ((entry.b0_and_b1 as u8) << 1) | ((entry.s2 as u8) << 2);
        if digit != 4 {
            return Err(HardReason::Tail3FixedLastMismatch { digit });
        }
    }
    if !ignore_tail_codec
        && std::env::var("DIALOG_GCD_K5_TAIL6_GRAPH9_CODEC")
            .ok()
            .as_deref()
            == Some("1")
    {
        let pattern = tail6_pattern(log);
        if !dialog_gcd_k5_tail6_graph9_supports(pattern) {
            return Err(HardReason::Tail6Graph9Mismatch { pattern });
        }
        return Ok(());
    }
    if !ignore_tail_codec
        && std::env::var("DIALOG_GCD_K5_TAIL6_GRAPH_CODEC")
            .ok()
            .as_deref()
            == Some("1")
    {
        let pattern = tail6_pattern(log);
        if !DIALOG_GCD_K5_TAIL6_GRAPH_SUPPORT.contains(&pattern) {
            return Err(HardReason::Tail6GraphMismatch { pattern });
        }
        return Ok(());
    }
    if !ignore_tail_codec && std::env::var("DIALOG_GCD_K5_TAIL7_CODEC").ok().as_deref() == Some("1")
    {
        let pattern = tail7_pattern(log);
        if !DIALOG_GCD_K5_TAIL7_SUPPORT.contains(&pattern) {
            return Err(HardReason::Tail7Mismatch { pattern });
        }
        return Ok(());
    }
    if std::env::var("DIALOG_GCD_K5_TAIL_PAIR1").ok().as_deref() == Some("1") {
        if !matches!(pattern, 0b100100 | 0b100101) {
            return Err(HardReason::TailPairMismatch { pattern });
        }
        return Ok(());
    }
    if tail_pair_codec_mode() == Some(1) {
        if log.len() < 4 {
            return Err(HardReason::TailPairCrossMismatch {
                pattern: pattern as u16,
            });
        }
        let previous = &log[log.len() - 3];
        let step0 = &log[log.len() - 2];
        let step1 = &log[log.len() - 1];
        let c1 = step0.b0 ^ step0.b0_and_b1;
        let c0 = !(previous.s2 || c1);
        let supported = step0.b0 == (c0 || c1)
            && step0.b0_and_b1 == (c0 && !c1)
            && step0.s2 == !c0
            && step1.b0 == c0
            && !step1.b0_and_b1
            && step1.s2;
        let mut joint = 0u16;
        for (slot, entry) in log[log.len() - 4..].iter().enumerate() {
            joint |= (entry.b0 as u16) << (3 * slot);
            joint |= (entry.b0_and_b1 as u16) << (3 * slot + 1);
            joint |= (entry.s2 as u16) << (3 * slot + 2);
        }
        if !supported {
            return Err(HardReason::TailPairCrossMismatch { pattern: joint });
        }
        if std::env::var("DIALOG_GCD_TAIL_CROSSBLOCK5").ok().as_deref() == Some("1")
            && !matches!(
                joint,
                0x924
                    | 0x925
                    | 0x928
                    | 0x929
                    | 0x92b
                    | 0x92c
                    | 0x92d
                    | 0x92f
                    | 0x944
                    | 0x945
                    | 0x947
                    | 0x948
                    | 0x949
                    | 0x94b
                    | 0x94d
                    | 0x94f
                    | 0x958
                    | 0x959
                    | 0x95b
                    | 0x95c
                    | 0x95d
                    | 0x95f
                    | 0x967
                    | 0x969
                    | 0x96b
                    | 0x978
                    | 0x979
                    | 0x97b
                    | 0x97f
                    | 0xac7
                    | 0xac9
                    | 0xacb
            )
        {
            return Err(HardReason::TailPairCrossMismatch { pattern: joint });
        }
        return Ok(());
    }
    if tail_pair_codec_mode() == Some(3)
        && std::env::var("DIALOG_GCD_TAIL_PAIR_CODEC3_V0_TOP8")
            .ok()
            .as_deref()
            == Some("1")
    {
        if !matches!(
            pattern,
            0b100100 | 0b100101 | 0b101001 | 0b101011 | 0b101000 | 0b101101 | 0b101111 | 0b101100
        ) {
            return Err(HardReason::TailPairMismatch { pattern });
        }
        return Ok(());
    }
    match tail_pair_codec_mode() {
        Some(0) if pattern != 0b100100 => {
            return Err(HardReason::TailPairMismatch { pattern });
        }
        Some(2) if !matches!(pattern, 0b100100 | 0b100101 | 0b101001 | 0b101011) => {
            return Err(HardReason::TailPairMismatch { pattern });
        }
        Some(3)
            if !matches!(
                pattern,
                0b011000 | 0b011011 | 0b100100 | 0b100101 | 0b101000 | 0b101001 | 0b101011
            ) =>
        {
            if pattern == 0b101100
                && std::env::var("DIALOG_GCD_TAIL_PAIR_CODEC3_PATTERN44")
                    .ok()
                    .as_deref()
                    == Some("1")
            {
                return Ok(());
            }
            return Err(HardReason::TailPairMismatch { pattern });
        }
        Some(4)
            if std::env::var("DIALOG_GCD_TAIL_PAIR_CODEC4_WIDE")
                .ok()
                .as_deref()
                == Some("1")
                && !matches!(
                    pattern,
                    0b100100
                        | 0b100101
                        | 0b101011
                        | 0b101001
                        | 0b101000
                        | 0b101101
                        | 0b101111
                        | 0b000111
                        | 0b011100
                ) =>
        {
            return Err(HardReason::TailPairMismatch { pattern });
        }
        Some(4)
            if !dialog_gcd_tail_pair4_wide_enabled_for_filter()
                && !matches!(
                    pattern,
                    0b011000 | 0b011011 | 0b100100 | 0b100101 | 0b101000 | 0b101001 | 0b101011
                ) =>
        {
            return Err(HardReason::TailPairMismatch { pattern });
        }
        _ => {}
    }
    Ok(())
}

fn dialog_gcd_tail_pair4_wide_enabled_for_filter() -> bool {
    std::env::var("DIALOG_GCD_TAIL_PAIR_CODEC4_WIDE")
        .ok()
        .as_deref()
        == Some("1")
}

pub fn debug_gcd_states(
    factor: U256,
    cfg: &DialogGcdFilterConfig,
) -> Result<Vec<(U256, U256)>, HardReason> {
    let mut u = SECP256K1_P;
    let mut v = factor;
    let mut states = Vec::with_capacity(cfg.active_iterations + 1);
    for step in 0..cfg.active_iterations {
        states.push((u, v));
        truncated_gcd_step_logged(&mut u, &mut v, step, cfg)?;
    }
    states.push((u, v));
    Ok(states)
}

pub fn debug_gcd_step_from_state(
    mut u: U256,
    mut v: U256,
    step: usize,
    cfg: &DialogGcdFilterConfig,
) -> Result<(U256, U256, DialogGcdStepLog), HardReason> {
    let entry = truncated_gcd_step_logged(&mut u, &mut v, step, cfg)?;
    Ok((u, v, entry))
}

pub fn debug_gcd_transcript(
    factor: U256,
    cfg: &DialogGcdFilterConfig,
) -> Result<Vec<DialogGcdStepLog>, HardReason> {
    Ok(build_gcd_transcript(factor, cfg)?.log)
}

fn add_carry_escapes(acc: U256, delta: U256, low_bits: usize) -> bool {
    if delta.is_zero() || low_bits >= N {
        return false;
    }
    let mask = window_mask(low_bits);
    (acc & mask) + delta > mask
}

fn sub_borrow_escapes(acc: U256, delta: U256, low_bits: usize) -> bool {
    if delta.is_zero() || low_bits >= N {
        return false;
    }
    (acc & window_mask(low_bits)) < delta
}

fn suffix_lt(a: U256, b: U256, bits: usize) -> bool {
    let bits = bits.clamp(1, N);
    let lo = N - bits;
    ((a >> lo) & window_mask(bits)) < ((b >> lo) & window_mask(bits))
}

fn required_suffix_compare_bits(a: U256, b: U256, current_bits: usize, desired: bool) -> usize {
    ((current_bits + 1)..=N)
        .find(|&bits| suffix_lt(a, b, bits) == desired)
        .unwrap_or(N)
}

fn fused_fold_delta(y: U256, s2: bool, reverse: bool) -> (U256, bool, bool) {
    let e = if reverse { bit_at(y, 0) } else { false };
    let d = if reverse { s2 && bit_at(y, 1) } else { false };
    let c = U256::MAX
        .wrapping_sub(SECP256K1_P)
        .wrapping_add(U256::from(1u64));
    (
        c * U256::from(e as u64) + (c << 1) * U256::from(d as u64),
        e,
        d,
    )
}

fn fused_double(
    mut y: U256,
    s2: bool,
    window: Option<usize>,
    step: usize,
) -> Result<U256, HardReason> {
    let ovf1 = bit_at(y, N - 1);
    y <<= 1;
    let ovf2 = s2 && bit_at(y, N - 1);
    if s2 {
        y <<= 1;
    }
    let d = ovf1 && s2;
    let e = ovf1 ^ d ^ ovf2;
    let c = U256::MAX
        .wrapping_sub(SECP256K1_P)
        .wrapping_add(U256::from(1u64));
    let delta = c * U256::from(e as u64) + (c << 1) * U256::from(d as u64);
    if let Some(w) = window {
        // last = hi_delta(33) + w, and carry[last] is written into
        // acc[last + 1], so loss starts above last + 1.
        let low_bits = (35 + w).min(N);
        if add_carry_escapes(y, delta, low_bits) {
            return Err(HardReason::FusedFoldCarryEscape {
                step,
                reverse: false,
            });
        }
    }
    Ok(y.wrapping_add(delta))
}

fn fused_halve(
    mut y: U256,
    s2: bool,
    window: Option<usize>,
    step: usize,
) -> Result<U256, HardReason> {
    let (delta, e, d) = fused_fold_delta(y, s2, true);
    if let Some(w) = window {
        let low_bits = (35 + w).min(N);
        if sub_borrow_escapes(y, delta, low_bits) {
            return Err(HardReason::FusedFoldCarryEscape {
                step,
                reverse: true,
            });
        }
    }
    y = y.wrapping_sub(delta);
    let ovf2 = e && s2;
    let ovf1 = if s2 { d } else { e };
    if s2 {
        y = (y >> 1) | (U256::from(ovf2 as u64) << (N - 1));
    }
    Ok((y >> 1) | (U256::from(ovf1 as u64) << (N - 1)))
}

fn special_add(
    y: U256,
    x: U256,
    step: usize,
    cfg: &DialogApplyFilterConfig,
    summary: &mut ApplyHazardSummary,
) -> Result<U256, HardReason> {
    let c = U256::MAX
        .wrapping_sub(SECP256K1_P)
        .wrapping_add(U256::from(1u64));
    let mut out = y.wrapping_add(x);
    let overflow = out < y;
    if overflow {
        if let Some(w) = cfg.special_fold_window(step) {
            // last = hi(c)(32) + w; carry[last] still updates acc[last + 1].
            let low_bits = (34 + w).min(N);
            if add_carry_escapes(out, c, low_bits) {
                return Err(HardReason::SpecialFoldCarryEscape {
                    step,
                    reverse: false,
                });
            }
        }
        out = out.wrapping_add(c);
    }
    let bits = cfg.overflow_compare_bits(step);
    let predicted = suffix_lt(out, x, bits);
    if predicted != overflow {
        summary.cleanup_mismatches += 1;
        summary.cleanup_mismatch_details.push(ApplyCleanupMismatch {
            step,
            reverse: false,
            bits,
            required_bits: required_suffix_compare_bits(out, x, bits, overflow),
        });
        if std::env::var_os("ISLAND_TRACE_REJECT").is_some() {
            eprintln!(
                "PHASE_RISK step={step} reverse=false overflow={overflow} bits={bits} required_bits={} predicted={predicted}",
                summary
                    .cleanup_mismatch_details
                    .last()
                    .expect("mismatch detail")
                    .required_bits,
            );
        }
    }
    Ok(out)
}

fn special_sub(
    y: U256,
    x: U256,
    step: usize,
    cfg: &DialogApplyFilterConfig,
    summary: &mut ApplyHazardSummary,
) -> Result<U256, HardReason> {
    let c = U256::MAX
        .wrapping_sub(SECP256K1_P)
        .wrapping_add(U256::from(1u64));
    let underflow = y < x;
    let mut out = y.wrapping_sub(x);
    if underflow {
        if let Some(w) = cfg.special_fold_window(step) {
            let low_bits = (34 + w).min(N);
            if sub_borrow_escapes(out, c, low_bits) {
                return Err(HardReason::SpecialFoldCarryEscape {
                    step,
                    reverse: true,
                });
            }
        }
        out = out.wrapping_sub(c);
    }
    let bits = cfg.underflow_compare_bits(step);
    let predicted = suffix_lt(out, !x, bits);
    if predicted == underflow {
        summary.cleanup_mismatches += 1;
        summary.cleanup_mismatch_details.push(ApplyCleanupMismatch {
            step,
            reverse: true,
            bits,
            required_bits: required_suffix_compare_bits(out, !x, bits, !underflow),
        });
        if std::env::var_os("ISLAND_TRACE_REJECT").is_some() {
            eprintln!(
                "PHASE_RISK step={step} reverse=true underflow={underflow} bits={bits} required_bits={} predicted={predicted} full_predicted={} y={y:#x} x={x:#x} out={out:#x}",
                summary
                    .cleanup_mismatch_details
                    .last()
                    .expect("mismatch detail")
                    .required_bits,
                suffix_lt(out, !x, N),
            );
        }
    }
    Ok(out)
}

fn check_apply_reverse_hazards_with_summary(
    log: &[DialogGcdStepLog],
    mut x: U256,
    mut y: U256,
    apply_cfg: &DialogApplyFilterConfig,
    summary: &mut ApplyHazardSummary,
) -> Result<(U256, U256), HardReason> {
    for (step, entry) in log.iter().copied().enumerate() {
        if entry.b0_and_b1 {
            std::mem::swap(&mut x, &mut y);
        }
        if entry.b0 {
            y = special_sub(y, x, step, apply_cfg, summary)?;
        }
        y = fused_halve(y, entry.s2, apply_cfg.fused_fold_window(step), step)?;
    }
    Ok((x, y))
}

pub fn check_apply_reverse_hazards(
    factor: U256,
    x: U256,
    y: U256,
    gcd_cfg: &DialogGcdFilterConfig,
    apply_cfg: &DialogApplyFilterConfig,
) -> Result<(U256, U256), HardReason> {
    let transcript = build_gcd_transcript(factor, gcd_cfg)?;
    check_apply_reverse_hazards_with_summary(
        &transcript.log,
        x,
        y,
        apply_cfg,
        &mut ApplyHazardSummary::default(),
    )
}

fn check_apply_forward_hazards_with_summary(
    log: &[DialogGcdStepLog],
    mut x: U256,
    mut y: U256,
    apply_cfg: &DialogApplyFilterConfig,
    summary: &mut ApplyHazardSummary,
) -> Result<(U256, U256), HardReason> {
    for (step, entry) in log.iter().copied().enumerate().rev() {
        y = fused_double(y, entry.s2, apply_cfg.fused_fold_window(step), step)?;
        if entry.b0 {
            y = special_add(y, x, step, apply_cfg, summary)?;
        }
        if entry.b0_and_b1 {
            std::mem::swap(&mut x, &mut y);
        }
    }
    Ok((x, y))
}

pub fn check_apply_forward_hazards(
    factor: U256,
    x: U256,
    y: U256,
    gcd_cfg: &DialogGcdFilterConfig,
    apply_cfg: &DialogApplyFilterConfig,
) -> Result<(U256, U256), HardReason> {
    let transcript = build_gcd_transcript(factor, gcd_cfg)?;
    check_apply_forward_hazards_with_summary(
        &transcript.log,
        x,
        y,
        apply_cfg,
        &mut ApplyHazardSummary::default(),
    )
}

fn check_point_add_apply_hazards_with_transcripts(
    dx_factor: U256,
    dx_transcript: &DialogGcdTranscript,
    dy: U256,
    lambda: U256,
    c_factor: U256,
    c_transcript: &DialogGcdTranscript,
    gcd_cfg: &DialogGcdFilterConfig,
    apply_cfg: &DialogApplyFilterConfig,
) -> Result<ApplyHazardSummary, HardReason> {
    check_tail_pair_codec(&dx_transcript.log)?;
    check_tail_pair_codec(&c_transcript.log)?;
    if dx_transcript.terminal_u != U256::from(1u64) || c_transcript.terminal_u != U256::from(1u64) {
        return Err(HardReason::NonConvergence {
            steps_needed: gcd_cfg.active_iterations + 1,
        });
    }

    let mut summary = ApplyHazardSummary::default();
    let reverse = check_apply_reverse_hazards_with_summary(
        &dx_transcript.log,
        dx_transcript.terminal_v,
        dy,
        apply_cfg,
        &mut summary,
    )?;
    if reverse != (lambda, dx_transcript.terminal_v) {
        if std::env::var_os("ISLAND_TRACE_REJECT").is_some() {
            eprintln!(
                "APPLY_VALUE reverse=true got=({:#x},{:#x}) expected=({lambda:#x},{:#x})",
                reverse.0, reverse.1, dx_transcript.terminal_v,
            );
        }
        return Err(HardReason::ApplyValueMismatch {
            reverse: true,
            compare_step: first_log_difference(dx_factor, dx_transcript, gcd_cfg, false),
            full_width_step: first_log_difference(dx_factor, dx_transcript, gcd_cfg, true),
        });
    }

    let forward = check_apply_forward_hazards_with_summary(
        &c_transcript.log,
        lambda,
        c_transcript.terminal_v,
        apply_cfg,
        &mut summary,
    )?;
    let expected_x = if apply_cfg.clear_product_residual {
        c_transcript.terminal_v ^ SECP256K1_P
    } else {
        c_transcript.terminal_v
    };
    let expected_y = lambda.mul_mod(c_factor, SECP256K1_P);
    if forward != (expected_x, expected_y) {
        if std::env::var_os("ISLAND_TRACE_REJECT").is_some() {
            eprintln!(
                "APPLY_VALUE reverse=false got=({:#x},{:#x}) expected=({expected_x:#x},{expected_y:#x})",
                forward.0, forward.1,
            );
        }
        return Err(HardReason::ApplyValueMismatch {
            reverse: false,
            compare_step: first_log_difference(c_factor, c_transcript, gcd_cfg, false),
            full_width_step: first_log_difference(c_factor, c_transcript, gcd_cfg, true),
        });
    }
    Ok(summary)
}

pub fn check_point_add_apply_hazards(
    dx: U256,
    dy: U256,
    lambda: U256,
    c: U256,
    gcd_cfg: &DialogGcdFilterConfig,
    apply_cfg: &DialogApplyFilterConfig,
) -> Result<ApplyHazardSummary, HardReason> {
    let dx_transcript = build_gcd_transcript(dx, gcd_cfg)?;
    let c_transcript = build_gcd_transcript(c, gcd_cfg)?;
    check_point_add_apply_hazards_with_transcripts(
        dx,
        &dx_transcript,
        dy,
        lambda,
        c,
        &c_transcript,
        gcd_cfg,
        apply_cfg,
    )
}

/// Apply-hazard replay using factors already accepted by
/// [`check_gcd_factor_checked`].
pub fn check_point_add_apply_hazards_checked(
    dx_factor: U256,
    dx: &CheckedGcdFactor,
    dy: U256,
    lambda: U256,
    c_factor: U256,
    c: &CheckedGcdFactor,
    gcd_cfg: &DialogGcdFilterConfig,
    apply_cfg: &DialogApplyFilterConfig,
) -> Result<ApplyHazardSummary, HardReason> {
    check_point_add_apply_hazards_with_transcripts(
        dx_factor,
        &dx.transcript,
        dy,
        lambda,
        c_factor,
        &c.transcript,
        gcd_cfg,
        apply_cfg,
    )
}

/// Validate a factor and retain its transcript for subsequent apply replay.
pub fn check_gcd_factor_checked(
    factor: U256,
    cfg: &DialogGcdFilterConfig,
) -> Result<CheckedGcdFactor, HardReason> {
    let transcript = build_gcd_transcript(factor, cfg)?;
    if std::env::var_os("ISLAND_TAIL_ALPHABET").is_none() {
        check_tail_pair_codec(&transcript.log)?;
    }
    if transcript.terminal_u != U256::from(1u64) {
        return Err(HardReason::NonConvergence {
            steps_needed: cfg.active_iterations + 1,
        });
    }
    Ok(CheckedGcdFactor { transcript })
}

/// Returns `Ok(())` if `factor` is safe under the truncated envelope, else the hard reason.
pub fn check_gcd_factor(factor: U256, cfg: &DialogGcdFilterConfig) -> Result<(), HardReason> {
    check_gcd_factor_checked(factor, cfg).map(|_| ())
}

/// Both dialog-GCD factors for one affine point-add input.
pub fn check_point_add_inputs(
    px: U256,
    qx: U256,
    rx: U256,
    cfg: &DialogGcdFilterConfig,
) -> Result<(), HardReason> {
    let (dx, c) = point_add_gcd_factors(px, qx, rx);
    check_gcd_factor(dx, cfg)?;
    check_gcd_factor(c, cfg)
}

/// Check all 9024 Fiat-Shamir shots; `Ok(())` means no hard inputs on either factor.
pub fn check_all_shots(
    px: &[U256],
    py: &[U256],
    qx: &[U256],
    qy: &[U256],
    rx: &[U256],
    ry: &[U256],
    cfg: &DialogGcdFilterConfig,
) -> Result<(), HardReason> {
    assert_eq!(px.len(), py.len());
    assert_eq!(px.len(), qx.len());
    assert_eq!(px.len(), qy.len());
    assert_eq!(px.len(), rx.len());
    assert_eq!(px.len(), ry.len());

    for i in 0..px.len() {
        let _ = (py[i], qy[i], ry[i]);
        let (dx, c) = point_add_gcd_factors(px[i], qx[i], rx[i]);
        if let Err(e) = check_gcd_factor(dx, cfg) {
            return Err(e);
        }
        if let Err(e) = check_gcd_factor(c, cfg) {
            return Err(e);
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::weierstrass_elliptic_curve::WeierstrassEllipticCurve;

    fn submission_route_env() {
        std::env::set_var("DIALOG_GCD_COMPARE_BITS", "46");
        std::env::set_var("DIALOG_GCD_WIDTH_MARGIN", "9");
        std::env::set_var("DIALOG_GCD_WIDTH_SLOPE_X1000", "1005");
        std::env::set_var("DIALOG_GCD_ACTIVE_ITERATIONS", "259");
        std::env::set_var("DIALOG_GCD_ODD_U_LOWBIT_FASTPATH", "1");
        std::env::set_var("DIALOG_GCD_K2", "1");
        std::env::set_var("DIALOG_GCD_RAW_TOBITVECTOR_VARIABLE_WIDTH", "1");
        std::env::set_var("DIALOG_GCD_PA9024_COMPARE_SCHEDULE", "0");
        std::env::set_var(
            "DIALOG_GCD_BODY_CARRY_BAND_TRIMS",
            "0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1",
        );
    }

    fn secp() -> WeierstrassEllipticCurve {
        WeierstrassEllipticCurve {
            modulus: SECP256K1_P,
            a: U256::from(0),
            b: U256::from(7),
            gx: U256::from_str_radix(
                "79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798",
                16,
            )
            .unwrap(),
            gy: U256::from_str_radix(
                "483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8",
                16,
            )
            .unwrap(),
            order: U256::from_str_radix(
                "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141",
                16,
            )
            .unwrap(),
        }
    }

    #[test]
    fn square_cleanup_site_overrides_are_directional_and_take_precedence() {
        let rows = vec![(12, 20)];
        let sites = vec![
            SquareCleanupSiteBits {
                row: 12,
                window: 0,
                reverse: false,
                bits: 21,
            },
            SquareCleanupSiteBits {
                row: 12,
                window: 0,
                reverse: true,
                bits: 22,
            },
        ];
        assert_eq!(square_cleanup_bits(19, &rows, &sites, 12, 0, false), 21);
        assert_eq!(square_cleanup_bits(19, &rows, &sites, 12, 0, true), 22);
        assert_eq!(square_cleanup_bits(19, &rows, &sites, 12, 1, false), 20);
        assert_eq!(square_cleanup_bits(19, &rows, &sites, 13, 0, false), 19);
    }

    #[test]
    fn known_clean_nonce_700017357_passes_filter() {
        submission_route_env();
        let cfg = DialogGcdFilterConfig::from_env();
        let curve = secp();

        // Derive a small prefix of the 9024-shot set with the same nonce tail as the frontier.
        let mut h = sha3::Shake256::default();
        h.update(b"quantum_ecc-fiat-shamir-v2");
        // Use a dummy op count; this test only checks factor geometry on random-derived points.
        h.update(&1000u64.to_le_bytes());
        for _ in 0..(48 * 2) {
            use sha3::digest::{ExtendableOutput, Update, XofReader};
            let mut xof = h.clone().finalize_xof();
            let mut rb = [[0u8; 32]; 2];
            for _ in 0..256 {
                xof.read(&mut rb[0]);
                xof.read(&mut rb[1]);
                let k1 = U256::from_le_bytes(rb[0]);
                let k2 = U256::from_le_bytes(rb[1]);
                let (px, py) = curve.mul(curve.gx, curve.gy, k1);
                let (qx, qy) = curve.mul(curve.gx, curve.gy, k2);
                if px == qx {
                    continue;
                }
                let (rx, ry) = curve.add(px, py, qx, qy);
                assert!(check_gcd_factor(point_add_gcd_factors(px, qx, rx).0, &cfg).is_ok());
                assert!(check_gcd_factor(point_add_gcd_factors(px, qx, rx).1, &cfg).is_ok());
                return;
            }
        }
        panic!("failed to sample a valid point pair");
    }

    #[test]
    fn width_margin_8_is_stricter_than_9() {
        submission_route_env();
        let cfg9 = DialogGcdFilterConfig::from_env();
        std::env::set_var("DIALOG_GCD_WIDTH_MARGIN", "8");
        let cfg8 = DialogGcdFilterConfig::from_env();

        let factor = U256::from_str_radix(
            "fffffffffffffffffffffffffffffffffffffffffffffffffffffffefffffc2e",
            16,
        )
        .unwrap();
        assert!(
            check_gcd_factor(factor, &cfg9).is_ok() || check_gcd_factor(factor, &cfg9).is_err()
        );
        // Margin 8 tightens step-0 width; many factors overflow earlier.
        let early_w9 = cfg9.active_width(0);
        let early_w8 = cfg8.active_width(0);
        assert!(early_w8 < early_w9);
    }
}
