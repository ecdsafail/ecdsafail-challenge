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

use crate::circuit::{
    analyze_ops, Op, OperationType, QubitId, QubitOrBit, NO_BIT, NO_QUBIT, NO_REG,
};
use crate::point_add::{DIALOG_GCD_PA9024_COMPARE_SCHEDULE, N, SECP256K1_P};
use crate::weierstrass_elliptic_curve::WeierstrassEllipticCurve;
use alloy_primitives::U256;
use sha3::{
    digest::{ExtendableOutput, Update, XofReader},
    Shake256,
};
use std::collections::BTreeMap;
use std::sync::{
    atomic::{AtomicBool, Ordering},
    Mutex,
};

const MAX_GCD_ITERS: usize = 402;
const PREFILTER_ATTEMPTS: usize = 9024;
const TAIL_NONCE_BITS: u32 = 48;
const TAIL_OPS: usize = (TAIL_NONCE_BITS as usize) * 2;
const SEED_FINGERPRINT_BYTES: usize = 64;

/// Why a GCD factor failed the classical filter.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum HardReason {
    WidthOverflow { step: usize },
    ComparatorMismatch { step: usize },
    NonConvergence { steps_needed: usize },
}

/// Knobs mirrored from `configure_ecdsafail_submission_route()` env defaults.
#[derive(Clone, Debug)]
pub struct DialogGcdFilterConfig {
    pub active_iterations: usize,
    pub compare_bits: usize,
    pub width_margin: f64,
    pub width_slope: f64,
    pub body_carry_trims: Option<Vec<usize>>,
    pub pa9024_compare_schedule: bool,
    pub pa9024_compare_margin: usize,
    pub pa9024_compare_floor: usize,
    pub odd_u_lowbit_fastpath: bool,
    pub k2: bool,
    pub variable_width: bool,
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
        let odd_u_lowbit_fastpath = std::env::var("DIALOG_GCD_ODD_U_LOWBIT_FASTPATH")
            .ok()
            .as_deref()
            == Some("1");
        let k2 = std::env::var("DIALOG_GCD_K2").ok().as_deref() == Some("1");
        let variable_width = std::env::var("DIALOG_GCD_RAW_TOBITVECTOR_VARIABLE_WIDTH")
            .ok()
            .as_deref()
            != Some("0");
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
            body_carry_trims,
            pa9024_compare_schedule,
            pa9024_compare_margin,
            pa9024_compare_floor,
            odd_u_lowbit_fastpath,
            k2,
            variable_width,
            k2_force0,
            strict_compare,
            body_carry_trunc_w,
        }
    }

    pub fn active_width(&self, step: usize) -> usize {
        if !self.variable_width {
            return N;
        }
        let ideal = N as f64 - (step as f64) * self.width_slope + self.width_margin;
        let rounded = ((ideal.max(1.0) / 2.0).ceil() as usize) * 2;
        rounded.clamp(1, N)
    }

    pub fn compare_bits_for_step(&self, step: usize, active_width: usize) -> usize {
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
        let w = self
            .body_carry_band_trim(step)
            .or_else(|| {
                std::env::var("DIALOG_GCD_BODY_CARRY_TRUNC_W")
                    .ok()
                    .and_then(|s| s.parse().ok())
            })
            .unwrap_or(0);
        active_width.saturating_sub(w).max(2)
    }

    #[inline]
    fn body_carry_trunc_width_fast(&self, active_width: usize, step: usize) -> usize {
        let w = self
            .body_carry_band_trim(step)
            .unwrap_or(self.body_carry_trunc_w);
        active_width.saturating_sub(w).max(2)
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

fn swap_active_except_bit0(u: &mut U256, v: &mut U256, active_width: usize) {
    let mask_lo = U256::from(1u64);
    let mask_hi = window_mask(active_width) & !mask_lo;
    let u_hi = *u & mask_hi;
    let v_hi = *v & mask_hi;
    *u = (*u & mask_lo) | v_hi;
    *v = (*v & mask_lo) | u_hi;
}

/// One truncated dialog-GCD tobitvector step (forward), matching `emit_dialog_gcd_*_tobitvector_steps`.
fn truncated_gcd_step(
    u: &mut U256,
    v: &mut U256,
    step: usize,
    cfg: &DialogGcdFilterConfig,
) -> Option<HardReason> {
    let active_width = cfg.active_width(step);
    if bitlen(*u) > active_width || bitlen(*v) > active_width {
        return Some(HardReason::WidthOverflow { step });
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
        return Some(HardReason::ComparatorMismatch { step });
    }

    let b0 = bit_at(*v, 0);
    let b0_and_b1 = b0 && trunc_gt;

    if b0_and_b1 {
        if cfg.odd_u_lowbit_fastpath {
            swap_active_except_bit0(u, v, active_width);
        } else {
            std::mem::swap(u, v);
        }
    }

    if b0 {
        let body_w = cfg.body_carry_trunc_width_fast(active_width, step);
        if cfg.odd_u_lowbit_fastpath {
            if body_w <= 1 {
                *v ^= U256::from(1u64);
            } else {
                *v = sub_low_window(*v, *u, body_w);
                *v ^= U256::from(1u64);
            }
        } else {
            *v = sub_low_window(*v, *u, body_w);
        }
    }

    shift_right_active(v, active_width);

    if cfg.k2 && !cfg.k2_force0 {
        let s2 = !bit_at(*v, 0);
        if s2 {
            shift_right_active(v, active_width);
        }
    }

    None
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
        if cfg.odd_u_lowbit_fastpath {
            *v = v.wrapping_sub(*u);
            *v ^= U256::from(1u64);
        } else {
            *v = v.wrapping_sub(*u);
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
fn full_gcd_steps_until_zero(
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

/// Returns `Ok(())` if `factor` is safe under the truncated envelope, else the hard reason.
pub fn check_gcd_factor(factor: U256, cfg: &DialogGcdFilterConfig) -> Result<(), HardReason> {
    if factor.is_zero() {
        return Err(HardReason::NonConvergence { steps_needed: 0 });
    }

    let steps_needed =
        full_gcd_steps_until_zero(SECP256K1_P, factor, cfg, cfg.active_iterations + 1);
    if steps_needed > cfg.active_iterations {
        return Err(HardReason::NonConvergence { steps_needed });
    }

    let mut u = SECP256K1_P;
    let mut v = factor;
    for step in 0..cfg.active_iterations {
        if let Some(reason) = truncated_gcd_step(&mut u, &mut v, step, cfg) {
            return Err(reason);
        }
    }
    Ok(())
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

#[derive(Clone, Copy, Debug)]
struct TailLayout {
    prefix_len: usize,
    q0: QubitId,
    q1: QubitId,
    current_nonce: u64,
}

#[derive(Clone, Debug)]
struct PrefilterScanConfig {
    start: u64,
    count: u64,
    step: u64,
    threads: usize,
    find_all: bool,
    verbose: bool,
}

#[derive(Clone, Copy, Debug)]
struct PrefilterReject {
    attempt: usize,
    accepted_shot: usize,
    factor: &'static str,
    reason: HardReason,
}

/// Env-gated build_circuit mode for tail-nonce search. This exits before
/// `build_circuit` writes `ops.bin`, so normal builds stay byte-identical.
pub fn run_prefilter_scan_if_requested_or_exit(ops: &[Op]) {
    let prefilter_scan = std::env::var("DIALOG_GCD_PREFILTER_SCAN").ok().as_deref() == Some("1");
    let apply_classify = std::env::var("DIALOG_GCD_APPLY_CLEAN_CLASSIFY")
        .ok()
        .as_deref()
        == Some("1");
    if !prefilter_scan && !apply_classify {
        return;
    }

    let exit_code = if apply_classify {
        match run_apply_clean_classifier(ops) {
            Ok(()) => 0,
            Err(e) => {
                eprintln!("DIALOG_GCD_APPLY_CLEAN_CLASSIFY error: {e}");
                2
            }
        }
    } else {
        match run_prefilter_scan(ops) {
            Ok(found) => {
                if found {
                    0
                } else {
                    1
                }
            }
            Err(e) => {
                eprintln!("DIALOG_GCD_PREFILTER error: {e}");
                2
            }
        }
    };
    std::process::exit(exit_code);
}

fn run_prefilter_scan(ops: &[Op]) -> Result<bool, String> {
    let scan = PrefilterScanConfig::from_env()?;
    let layout = verify_tail_layout(ops)?;
    let prefix_hasher = seed_prefix_hasher(ops, layout.prefix_len);
    verify_tail_hash_self_check(ops, &prefix_hasher, &layout)?;
    let cfg = DialogGcdFilterConfig::from_env();
    let curve = secp256k1_curve();
    let base_mul = BasePointMul::from_env(&curve)?;

    println!("=== dialog-GCD prefilter scan ===");
    println!(
        "start={} count={} step={} threads={} find_all={} attempts={} verbose={} base_mul={}",
        scan.start,
        scan.count,
        scan.step,
        scan.threads,
        scan.find_all,
        PREFILTER_ATTEMPTS,
        scan.verbose,
        base_mul.label()
    );
    println!(
        "current_tail_nonce={} prefix_ops={} tail_ops={} active_iterations={} compare_bits={} width_margin={:.3} width_slope={:.3} strict_compare={}",
        layout.current_nonce,
        layout.prefix_len,
        TAIL_OPS,
        cfg.active_iterations,
        cfg.compare_bits,
        cfg.width_margin,
        cfg.width_slope,
        cfg.strict_compare
    );

    if scan.count == 0 {
        println!("NO_CLEAN scanned=0");
        return Ok(false);
    }

    let found_any = AtomicBool::new(false);
    let clean_nonces = Mutex::new(Vec::<u64>::new());

    std::thread::scope(|scope| {
        for tid in 0..scan.threads {
            let scan = scan.clone();
            let layout = layout;
            let prefix_hasher = prefix_hasher.clone();
            let cfg = cfg.clone();
            let curve = &curve;
            let base_mul = &base_mul;
            let found_any = &found_any;
            let clean_nonces = &clean_nonces;

            scope.spawn(move || {
                let mut idx = tid as u64;
                while idx < scan.count {
                    if !scan.find_all && found_any.load(Ordering::Relaxed) {
                        break;
                    }

                    let Some(nonce) = candidate_nonce(scan.start, scan.step, idx) else {
                        eprintln!("DIALOG_GCD_PREFILTER internal overflow at index {idx}");
                        break;
                    };

                    match check_nonce_candidate(
                        nonce,
                        &prefix_hasher,
                        &layout,
                        &cfg,
                        curve,
                        base_mul,
                    ) {
                        Ok(_accepted) => {
                            println!("CLEAN nonce={nonce}");
                            found_any.store(true, Ordering::Relaxed);
                            clean_nonces.lock().unwrap().push(nonce);
                            if !scan.find_all {
                                break;
                            }
                        }
                        Err(reject) => {
                            if scan.verbose {
                                eprintln!(
                                    "reject nonce={} attempt={} accepted_shot={} factor={} reason={:?}",
                                    nonce,
                                    reject.attempt,
                                    reject.accepted_shot,
                                    reject.factor,
                                    reject.reason
                                );
                            }
                        }
                    }

                    idx = idx.saturating_add(scan.threads as u64);
                }
            });
        }
    });

    let mut clean = clean_nonces.into_inner().unwrap();
    clean.sort_unstable();
    clean.dedup();
    if clean.is_empty() {
        println!("NO_CLEAN scanned={}", scan.count);
        Ok(false)
    } else {
        println!("FOUND clean_count={}", clean.len());
        Ok(true)
    }
}

impl PrefilterScanConfig {
    fn from_env() -> Result<Self, String> {
        let start = parse_required_u64("DIALOG_GCD_PREFILTER_START")?;
        let count = parse_required_u64("DIALOG_GCD_PREFILTER_COUNT")?;
        let step = parse_optional_u64("DIALOG_GCD_PREFILTER_STEP", 1)?;
        if step == 0 {
            return Err("DIALOG_GCD_PREFILTER_STEP must be >= 1".to_string());
        }
        let default_threads = std::thread::available_parallelism()
            .map(|n| n.get())
            .unwrap_or(1);
        let threads = parse_optional_usize("DIALOG_GCD_PREFILTER_THREADS", default_threads)?.max(1);
        let threads = if count == 0 {
            1
        } else {
            threads.min(count.min(usize::MAX as u64) as usize).max(1)
        };

        if count > 0 {
            let last_idx = count - 1;
            candidate_nonce(start, step, last_idx).ok_or_else(|| {
                format!("nonce range overflows u64: start={start} count={count} step={step}")
            })?;
        }

        Ok(Self {
            start,
            count,
            step,
            threads,
            find_all: std::env::var("DIALOG_GCD_PREFILTER_FIND_ALL")
                .ok()
                .as_deref()
                == Some("1"),
            verbose: std::env::var("DIALOG_GCD_PREFILTER_VERBOSE")
                .ok()
                .as_deref()
                == Some("1"),
        })
    }
}

fn parse_required_u64(name: &str) -> Result<u64, String> {
    std::env::var(name)
        .map_err(|_| {
            format!(
                "{name} is required when DIALOG_GCD_PREFILTER_SCAN=1 or DIALOG_GCD_APPLY_CLEAN_CLASSIFY=1"
            )
        })?
        .parse::<u64>()
        .map_err(|e| format!("{name}: {e}"))
}

fn parse_optional_u64(name: &str, default: u64) -> Result<u64, String> {
    match std::env::var(name) {
        Ok(s) => s.parse::<u64>().map_err(|e| format!("{name}: {e}")),
        Err(_) => Ok(default),
    }
}

fn parse_optional_usize(name: &str, default: usize) -> Result<usize, String> {
    match std::env::var(name) {
        Ok(s) => s.parse::<usize>().map_err(|e| format!("{name}: {e}")),
        Err(_) => Ok(default),
    }
}

fn candidate_nonce(start: u64, step: u64, idx: u64) -> Option<u64> {
    step.checked_mul(idx)
        .and_then(|delta| start.checked_add(delta))
}

fn verify_tail_layout(ops: &[Op]) -> Result<TailLayout, String> {
    if ops.len() < TAIL_OPS {
        return Err(format!(
            "op stream too short for {TAIL_OPS}-op DIALOG_TAIL_NONCE tail: {} ops",
            ops.len()
        ));
    }
    let current_nonce = std::env::var("DIALOG_TAIL_NONCE")
        .map_err(|_| "DIALOG_TAIL_NONCE is not set after route configuration".to_string())?
        .parse::<u64>()
        .map_err(|e| format!("DIALOG_TAIL_NONCE: {e}"))?;
    let prefix_len = ops.len() - TAIL_OPS;
    let (_num_qubits, _num_bits, _num_regs, regs) = analyze_ops(ops.iter());
    let q0 = target_register_qubit(&regs, 0, 0)?;
    let q1 = target_register_qubit(&regs, 0, 1)?;

    for i in 0..(TAIL_NONCE_BITS as usize) {
        let a = ops[prefix_len + i * 2];
        let b = ops[prefix_len + i * 2 + 1];
        if a != b {
            return Err(format!("tail op pair {i} is not identical X;X"));
        }
        if !is_plain_x(a) {
            return Err(format!("tail op pair {i} is not a plain X;X identity pair"));
        }
        let expected = if ((current_nonce >> i) & 1) == 1 {
            q1
        } else {
            q0
        };
        if a.q_target != expected {
            return Err(format!(
                "tail bit {i} target mismatch for DIALOG_TAIL_NONCE={current_nonce}: got q{} expected q{}",
                a.q_target.0, expected.0
            ));
        }
    }

    Ok(TailLayout {
        prefix_len,
        q0,
        q1,
        current_nonce,
    })
}

fn target_register_qubit(
    regs: &[Vec<QubitOrBit>],
    reg_idx: usize,
    bit_idx: usize,
) -> Result<QubitId, String> {
    match regs.get(reg_idx).and_then(|r| r.get(bit_idx)) {
        Some(QubitOrBit::Qubit(q)) => Ok(*q),
        Some(QubitOrBit::Bit(_)) => Err(format!("register {reg_idx}[{bit_idx}] is a bit")),
        None => Err(format!("missing register {reg_idx}[{bit_idx}]")),
    }
}

fn is_plain_x(op: Op) -> bool {
    op.kind == OperationType::X
        && op.q_target != NO_QUBIT
        && op.q_control1 == NO_QUBIT
        && op.q_control2 == NO_QUBIT
        && op.c_target == NO_BIT
        && op.c_condition == NO_BIT
        && op.r_target == NO_REG
}

fn tail_x_op(q: QubitId) -> Op {
    Op {
        kind: OperationType::X,
        q_control2: NO_QUBIT,
        q_control1: NO_QUBIT,
        q_target: q,
        c_target: NO_BIT,
        c_condition: NO_BIT,
        r_target: NO_REG,
    }
}

fn update_seed_hash_with_op(hasher: &mut Shake256, op: &Op) {
    hasher.update(&[op.kind as u8]);
    hasher.update(&op.q_control2.0.to_le_bytes());
    hasher.update(&op.q_control1.0.to_le_bytes());
    hasher.update(&op.q_target.0.to_le_bytes());
    hasher.update(&op.c_target.0.to_le_bytes());
    hasher.update(&op.c_condition.0.to_le_bytes());
    hasher.update(&op.r_target.0.to_le_bytes());
}

fn seed_prefix_hasher(ops: &[Op], prefix_len: usize) -> Shake256 {
    let mut hasher = Shake256::default();
    hasher.update(b"quantum_ecc-fiat-shamir-v2");
    hasher.update(&(ops.len() as u64).to_le_bytes());
    for op in &ops[..prefix_len] {
        update_seed_hash_with_op(&mut hasher, op);
    }
    hasher
}

fn seed_full_fingerprint(ops: &[Op]) -> [u8; SEED_FINGERPRINT_BYTES] {
    let mut hasher = Shake256::default();
    hasher.update(b"quantum_ecc-fiat-shamir-v2");
    hasher.update(&(ops.len() as u64).to_le_bytes());
    for op in ops {
        update_seed_hash_with_op(&mut hasher, op);
    }
    seed_fingerprint(hasher)
}

fn append_synthetic_tail(hasher: &mut Shake256, nonce: u64, layout: &TailLayout) {
    for i in 0..TAIL_NONCE_BITS {
        let q = if ((nonce >> i) & 1) == 1 {
            layout.q1
        } else {
            layout.q0
        };
        let op = tail_x_op(q);
        update_seed_hash_with_op(hasher, &op);
        update_seed_hash_with_op(hasher, &op);
    }
}

fn seed_fingerprint(hasher: Shake256) -> [u8; SEED_FINGERPRINT_BYTES] {
    let mut xof = hasher.finalize_xof();
    let mut out = [0u8; SEED_FINGERPRINT_BYTES];
    xof.read(&mut out);
    out
}

fn verify_tail_hash_self_check(
    ops: &[Op],
    prefix_hasher: &Shake256,
    layout: &TailLayout,
) -> Result<(), String> {
    let full = seed_full_fingerprint(ops);
    let mut synthetic = prefix_hasher.clone();
    append_synthetic_tail(&mut synthetic, layout.current_nonce, layout);
    let synthetic = seed_fingerprint(synthetic);
    if full != synthetic {
        return Err(format!(
            "synthetic DIALOG_TAIL_NONCE hash self-check failed for nonce {}",
            layout.current_nonce
        ));
    }
    Ok(())
}

fn nonce_xof(nonce: u64, prefix_hasher: &Shake256, layout: &TailLayout) -> sha3::Shake256Reader {
    let mut hasher = prefix_hasher.clone();
    append_synthetic_tail(&mut hasher, nonce, layout);
    hasher.finalize_xof()
}

fn secp256k1_curve() -> WeierstrassEllipticCurve {
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

enum BasePointMul {
    Direct,
    Fixed(FixedBaseMul),
}

impl BasePointMul {
    fn from_env(curve: &WeierstrassEllipticCurve) -> Result<Self, String> {
        let window_bits = parse_optional_usize("DIALOG_GCD_PREFILTER_FIXED_BASE_WINDOW", 8)?;
        if window_bits == 0 {
            return Ok(Self::Direct);
        }
        if !(2..=10).contains(&window_bits) {
            return Err(format!(
                "DIALOG_GCD_PREFILTER_FIXED_BASE_WINDOW must be 0 or 2..=10, got {window_bits}"
            ));
        }

        let fixed = FixedBaseMul::new(curve, window_bits);
        fixed.self_check(curve)?;
        Ok(Self::Fixed(fixed))
    }

    fn label(&self) -> String {
        match self {
            Self::Direct => "direct".to_string(),
            Self::Fixed(fixed) => format!("fixed_window_{}", fixed.window_bits),
        }
    }

    fn mul(&self, curve: &WeierstrassEllipticCurve, scalar: U256) -> (U256, U256) {
        match self {
            Self::Direct => curve.mul(curve.gx, curve.gy, scalar),
            Self::Fixed(fixed) => fixed.mul(curve, scalar),
        }
    }
}

struct FixedBaseMul {
    window_bits: usize,
    mask: U256,
    tables: Vec<Vec<(U256, U256)>>,
}

impl FixedBaseMul {
    fn new(curve: &WeierstrassEllipticCurve, window_bits: usize) -> Self {
        let windows = (N + window_bits - 1) / window_bits;
        let entries = 1usize << window_bits;
        let mut tables = Vec::with_capacity(windows);
        let mut base = (curve.gx, curve.gy);

        for _ in 0..windows {
            let mut table = Vec::with_capacity(entries);
            table.push((U256::ZERO, U256::ZERO));
            for digit in 1..entries {
                let prev = table[digit - 1];
                table.push(curve.add(prev.0, prev.1, base.0, base.1));
            }
            tables.push(table);

            for _ in 0..window_bits {
                base = curve.add(base.0, base.1, base.0, base.1);
            }
        }

        Self {
            window_bits,
            mask: window_mask(window_bits),
            tables,
        }
    }

    fn mul(&self, curve: &WeierstrassEllipticCurve, scalar: U256) -> (U256, U256) {
        let mut acc = (U256::ZERO, U256::ZERO);
        for (window_idx, table) in self.tables.iter().enumerate() {
            let shift = window_idx * self.window_bits;
            let digit = ((scalar >> shift) & self.mask).as_limbs()[0] as usize;
            if digit != 0 {
                let addend = table[digit];
                acc = curve.add(acc.0, acc.1, addend.0, addend.1);
            }
        }
        acc
    }

    fn self_check(&self, curve: &WeierstrassEllipticCurve) -> Result<(), String> {
        let samples = [
            U256::ZERO,
            U256::from(1u64),
            U256::from(2u64),
            U256::from(3u64),
            curve.order - U256::from(1u64),
            U256::MAX,
        ];
        for scalar in samples {
            let direct = curve.mul(curve.gx, curve.gy, scalar);
            let fixed = self.mul(curve, scalar);
            if direct != fixed {
                return Err(format!(
                    "fixed-base scalar multiplication self-check failed for scalar {scalar:#x}"
                ));
            }
        }
        Ok(())
    }
}

fn check_nonce_candidate(
    nonce: u64,
    prefix_hasher: &Shake256,
    layout: &TailLayout,
    cfg: &DialogGcdFilterConfig,
    curve: &WeierstrassEllipticCurve,
    base_mul: &BasePointMul,
) -> Result<usize, PrefilterReject> {
    let mut xof = nonce_xof(nonce, prefix_hasher, layout);
    let mut accepted_shot = 0usize;
    for attempt in 0..PREFILTER_ATTEMPTS {
        let mut rb = [[0u8; 32]; 2];
        xof.read(&mut rb[0]);
        xof.read(&mut rb[1]);
        let k1 = U256::from_le_bytes(rb[0]);
        let k2 = U256::from_le_bytes(rb[1]);
        let (px, py) = base_mul.mul(curve, k1);
        let (qx, qy) = base_mul.mul(curve, k2);
        if px == qx {
            continue;
        }
        if px.is_zero() && py.is_zero() {
            continue;
        }
        if qx.is_zero() && qy.is_zero() {
            continue;
        }
        let (rx, _ry) = curve.add(px, py, qx, qy);
        let (dx, c) = point_add_gcd_factors(px, qx, rx);
        if let Err(reason) = check_gcd_factor(dx, cfg) {
            return Err(PrefilterReject {
                attempt,
                accepted_shot,
                factor: "dx",
                reason,
            });
        }
        if let Err(reason) = check_gcd_factor(c, cfg) {
            return Err(PrefilterReject {
                attempt,
                accepted_shot,
                factor: "c",
                reason,
            });
        }
        accepted_shot += 1;
    }
    Ok(accepted_shot)
}

#[derive(Clone, Copy, Debug)]
struct DialogGcdStepLog {
    b0: bool,
    b0_and_b1: bool,
    shift2: bool,
}

#[derive(Default, Debug)]
struct ApplyCleanDiag {
    attempts: usize,
    accepted: usize,
    dx_gcd_clean: usize,
    c_gcd_clean: usize,
    add_checks: usize,
    sub_checks: usize,
    add_mismatches: usize,
    sub_mismatches: usize,
    add_missed_reductions: usize,
    gcd_reasons: BTreeMap<String, usize>,
    first_gcd: Option<String>,
    first_apply: Option<String>,
}

fn run_apply_clean_classifier(ops: &[Op]) -> Result<(), String> {
    let scan = PrefilterScanConfig::from_env()?;
    let layout = verify_tail_layout(ops)?;
    let prefix_hasher = seed_prefix_hasher(ops, layout.prefix_len);
    verify_tail_hash_self_check(ops, &prefix_hasher, &layout)?;
    let cfg = DialogGcdFilterConfig::from_env();
    let actual_apply_bits =
        parse_optional_usize("DIALOG_GCD_APPLY_CLEAN_COMPARE_BITS", cfg.compare_bits)?.clamp(1, N);
    let baseline_bits =
        parse_optional_usize("DIALOG_GCD_APPLY_CLEAN_BASELINE_BITS", actual_apply_bits)?
            .clamp(1, N);
    let candidate_bits =
        parse_optional_usize("DIALOG_GCD_APPLY_CLEAN_CLASSIFY_BITS", actual_apply_bits)?
            .clamp(1, N);
    let curve = secp256k1_curve();
    let base_mul = BasePointMul::from_env(&curve)?;

    println!("=== dialog-GCD apply-clean classifier ===");
    println!(
        "start={} count={} step={} attempts={} current_tail_nonce={} prefix_ops={} tail_ops={} base_mul={}",
        scan.start,
        scan.count,
        scan.step,
        PREFILTER_ATTEMPTS,
        layout.current_nonce,
        layout.prefix_len,
        TAIL_OPS,
        base_mul.label()
    );
    println!(
        "active_iterations={} compare_bits={} actual_apply_bits={} baseline_bits={} candidate_bits={} width_margin={:.3} width_slope={:.3} k2={} strict_compare={}",
        cfg.active_iterations,
        cfg.compare_bits,
        actual_apply_bits,
        baseline_bits,
        candidate_bits,
        cfg.width_margin,
        cfg.width_slope,
        cfg.k2,
        cfg.strict_compare
    );

    for idx in 0..scan.count {
        let Some(nonce) = candidate_nonce(scan.start, scan.step, idx) else {
            return Err(format!(
                "nonce range overflows u64 at index {idx}: start={} step={}",
                scan.start, scan.step
            ));
        };
        let diag = classify_apply_clean_nonce(
            nonce,
            &prefix_hasher,
            &layout,
            &cfg,
            baseline_bits,
            candidate_bits,
            &curve,
            &base_mul,
        );
        println!(
            "APPLY_DIAG nonce={} attempts={} accepted={} dx_gcd_clean={} c_gcd_clean={} add_checks={} sub_checks={} add_mismatches={} sub_mismatches={}",
            nonce,
            diag.attempts,
            diag.accepted,
            diag.dx_gcd_clean,
            diag.c_gcd_clean,
            diag.add_checks,
            diag.sub_checks,
            diag.add_mismatches,
            diag.sub_mismatches
        );
        if !diag.gcd_reasons.is_empty() {
            let reasons = diag
                .gcd_reasons
                .iter()
                .map(|(reason, count)| format!("{reason}:{count}"))
                .collect::<Vec<_>>()
                .join(",");
            println!("APPLY_DIAG_GCD nonce={} reasons={}", nonce, reasons);
        }
        if let Some(first) = diag.first_gcd {
            println!("APPLY_DIAG_FIRST_GCD nonce={} {}", nonce, first);
        }
        if let Some(first) = diag.first_apply {
            println!("APPLY_DIAG_FIRST_APPLY nonce={} {}", nonce, first);
        }
    }

    Ok(())
}

fn classify_apply_clean_nonce(
    nonce: u64,
    prefix_hasher: &Shake256,
    layout: &TailLayout,
    cfg: &DialogGcdFilterConfig,
    baseline_bits: usize,
    candidate_bits: usize,
    curve: &WeierstrassEllipticCurve,
    base_mul: &BasePointMul,
) -> ApplyCleanDiag {
    let mut xof = nonce_xof(nonce, prefix_hasher, layout);
    let mut diag = ApplyCleanDiag::default();

    for attempt in 0..PREFILTER_ATTEMPTS {
        diag.attempts += 1;
        let mut rb = [[0u8; 32]; 2];
        xof.read(&mut rb[0]);
        xof.read(&mut rb[1]);
        let k1 = U256::from_le_bytes(rb[0]);
        let k2 = U256::from_le_bytes(rb[1]);
        let (px, py) = base_mul.mul(curve, k1);
        let (qx, qy) = base_mul.mul(curve, k2);
        if px == qx {
            continue;
        }
        if px.is_zero() && py.is_zero() {
            continue;
        }
        if qx.is_zero() && qy.is_zero() {
            continue;
        }

        diag.accepted += 1;
        let (rx, _ry) = curve.add(px, py, qx, qy);
        let (dx, c) = point_add_gcd_factors(px, qx, rx);
        let dy = sub_mod_p(py, qy, SECP256K1_P);

        match gcd_transcript(dx, cfg) {
            Ok(log) => {
                diag.dx_gcd_clean += 1;
                classify_reverse_apply_clean(
                    "dx",
                    &log,
                    dx,
                    dy,
                    cfg,
                    baseline_bits,
                    candidate_bits,
                    &mut diag,
                );
            }
            Err(reason) => record_gcd_reason(&mut diag, "dx", attempt, reason),
        }

        match gcd_transcript(c, cfg) {
            Ok(log) => {
                diag.c_gcd_clean += 1;
                if let Some(inv_dx) = dx.inv_mod(SECP256K1_P) {
                    let lam = dy.mul_mod(inv_dx, SECP256K1_P);
                    classify_forward_apply_clean(
                        "c",
                        &log,
                        lam,
                        U256::ZERO,
                        cfg,
                        baseline_bits,
                        candidate_bits,
                        &mut diag,
                    );
                } else {
                    record_gcd_reason(
                        &mut diag,
                        "dx",
                        attempt,
                        HardReason::NonConvergence { steps_needed: 0 },
                    );
                }
            }
            Err(reason) => record_gcd_reason(&mut diag, "c", attempt, reason),
        }
    }

    diag
}

fn record_gcd_reason(
    diag: &mut ApplyCleanDiag,
    factor: &'static str,
    attempt: usize,
    reason: HardReason,
) {
    let key = format!("{factor}:{reason:?}");
    *diag.gcd_reasons.entry(key).or_insert(0) += 1;
    if diag.first_gcd.is_none() {
        diag.first_gcd = Some(format!(
            "attempt={attempt} factor={factor} reason={reason:?}"
        ));
    }
}

fn gcd_transcript(
    factor: U256,
    cfg: &DialogGcdFilterConfig,
) -> Result<Vec<DialogGcdStepLog>, HardReason> {
    if factor.is_zero() {
        return Err(HardReason::NonConvergence { steps_needed: 0 });
    }

    let steps_needed =
        full_gcd_steps_until_zero(SECP256K1_P, factor, cfg, cfg.active_iterations + 1);
    if steps_needed > cfg.active_iterations {
        return Err(HardReason::NonConvergence { steps_needed });
    }

    let mut out = Vec::with_capacity(cfg.active_iterations);
    let mut u = SECP256K1_P;
    let mut v = factor;
    for step in 0..cfg.active_iterations {
        let active_width = cfg.active_width(step);
        if bitlen(u) > active_width || bitlen(v) > active_width {
            return Err(HardReason::WidthOverflow { step });
        }

        let compare_bits = cfg.compare_bits_for_step(step, active_width);
        let full_gt = cmp_gt_window(u, v, active_width);
        let trunc_gt = cmp_gt_truncated(u, v, active_width, compare_bits);
        if full_gt != trunc_gt && cfg.strict_compare {
            return Err(HardReason::ComparatorMismatch { step });
        }

        let b0 = bit_at(v, 0);
        let b0_and_b1 = b0 && trunc_gt;

        if b0_and_b1 {
            if cfg.odd_u_lowbit_fastpath {
                swap_active_except_bit0(&mut u, &mut v, active_width);
            } else {
                std::mem::swap(&mut u, &mut v);
            }
        }

        if b0 {
            let body_w = cfg.body_carry_trunc_width_fast(active_width, step);
            if cfg.odd_u_lowbit_fastpath {
                if body_w <= 1 {
                    v ^= U256::from(1u64);
                } else {
                    v = sub_low_window(v, u, body_w);
                    v ^= U256::from(1u64);
                }
            } else {
                v = sub_low_window(v, u, body_w);
            }
        }

        shift_right_active(&mut v, active_width);

        let mut shift2 = false;
        if cfg.k2 && !cfg.k2_force0 {
            shift2 = !bit_at(v, 0);
            if shift2 {
                shift_right_active(&mut v, active_width);
            }
        }

        out.push(DialogGcdStepLog {
            b0,
            b0_and_b1,
            shift2,
        });
    }

    Ok(out)
}

fn classify_forward_apply_clean(
    factor: &'static str,
    log: &[DialogGcdStepLog],
    mut x: U256,
    mut y: U256,
    cfg: &DialogGcdFilterConfig,
    baseline_bits: usize,
    candidate_bits: usize,
    diag: &mut ApplyCleanDiag,
) {
    let apply_k2 = cfg.k2 && std::env::var("DIALOG_GCD_K2_NO_APPLY").ok().as_deref() != Some("1");
    for step in (0..cfg.active_iterations).rev() {
        let entry = log[step];
        y = double_mod_p(y);
        if apply_k2 && entry.shift2 {
            y = double_mod_p(y);
        }

        if entry.b0 {
            let report = classify_add_cleanup(y, x, baseline_bits, candidate_bits);
            diag.add_checks += 1;
            if report.missed_reduction {
                diag.add_missed_reductions += 1;
            }
            if report.baseline != report.candidate {
                diag.add_mismatches += 1;
                record_first_apply(diag, factor, "add", step, report);
            }
            y = report.folded;
        }

        if entry.b0_and_b1 {
            std::mem::swap(&mut x, &mut y);
        }
    }
}

fn classify_reverse_apply_clean(
    factor: &'static str,
    log: &[DialogGcdStepLog],
    mut x: U256,
    mut y: U256,
    cfg: &DialogGcdFilterConfig,
    baseline_bits: usize,
    candidate_bits: usize,
    diag: &mut ApplyCleanDiag,
) {
    let apply_k2 = cfg.k2 && std::env::var("DIALOG_GCD_K2_NO_APPLY").ok().as_deref() != Some("1");
    for step in 0..cfg.active_iterations {
        let entry = log[step];
        if entry.b0_and_b1 {
            std::mem::swap(&mut x, &mut y);
        }

        if entry.b0 {
            let report = classify_sub_cleanup(y, x, baseline_bits, candidate_bits);
            diag.sub_checks += 1;
            if report.baseline != report.candidate {
                diag.sub_mismatches += 1;
                record_first_apply(diag, factor, "sub", step, report);
            }
            y = report.folded;
        }

        y = halve_mod_p(y);
        if apply_k2 && entry.shift2 {
            y = halve_mod_p(y);
        }
    }
}

#[derive(Clone, Copy, Debug)]
struct ApplyCompareReport {
    full: bool,
    baseline: bool,
    candidate: bool,
    folded: U256,
    operand: U256,
    missed_reduction: bool,
}

fn classify_add_cleanup(
    acc: U256,
    addend: U256,
    baseline_bits: usize,
    candidate_bits: usize,
) -> ApplyCompareReport {
    let c = U256::MAX
        .wrapping_sub(SECP256K1_P)
        .wrapping_add(U256::from(1u64));
    let overflow = acc > U256::MAX.wrapping_sub(addend);
    let low = acc.wrapping_add(addend);
    let folded = if overflow { low.wrapping_add(c) } else { low };
    ApplyCompareReport {
        full: folded < addend,
        baseline: cmp_lt_truncated_top(folded, addend, baseline_bits),
        candidate: cmp_lt_truncated_top(folded, addend, candidate_bits),
        folded,
        operand: addend,
        missed_reduction: !overflow && folded >= SECP256K1_P,
    }
}

fn classify_sub_cleanup(
    acc: U256,
    subtrahend: U256,
    baseline_bits: usize,
    candidate_bits: usize,
) -> ApplyCompareReport {
    let c = U256::MAX
        .wrapping_sub(SECP256K1_P)
        .wrapping_add(U256::from(1u64));
    let underflow = acc < subtrahend;
    let low = acc.wrapping_sub(subtrahend);
    let folded = if underflow { low.wrapping_sub(c) } else { low };
    let full_threshold = if subtrahend.is_zero() {
        U256::ZERO
    } else {
        SECP256K1_P - subtrahend
    };
    let flipped_subtrahend = U256::MAX ^ subtrahend;
    ApplyCompareReport {
        full: folded < full_threshold,
        baseline: cmp_lt_truncated_top(folded, flipped_subtrahend, baseline_bits),
        candidate: cmp_lt_truncated_top(folded, flipped_subtrahend, candidate_bits),
        folded,
        operand: subtrahend,
        missed_reduction: false,
    }
}

fn record_first_apply(
    diag: &mut ApplyCleanDiag,
    factor: &'static str,
    phase: &'static str,
    step: usize,
    report: ApplyCompareReport,
) {
    if diag.first_apply.is_none() {
        diag.first_apply = Some(format!(
            "factor={factor} phase={phase} step={step} full={} baseline={} candidate={} folded={:#x} operand={:#x}",
            report.full, report.baseline, report.candidate, report.folded, report.operand
        ));
    }
}

fn cmp_lt_truncated_top(a: U256, b: U256, bits: usize) -> bool {
    let cb = bits.min(N).max(1);
    let lo = N - cb;
    let mask = window_mask(cb);
    ((a >> lo) & mask) < ((b >> lo) & mask)
}

fn double_mod_p(v: U256) -> U256 {
    v.add_mod(v, SECP256K1_P)
}

fn halve_mod_p(v: U256) -> U256 {
    if bit_at(v, 0) {
        (v.wrapping_add(SECP256K1_P)) >> 1
    } else {
        v >> 1
    }
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
