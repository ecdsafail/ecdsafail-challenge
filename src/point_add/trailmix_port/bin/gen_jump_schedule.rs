//! Register-shrink schedule generator for the jump-GCD (Schrottenloher),
//! modeled on `gen_shrunken_pz_schedule.rs`: sample tens of millions of binary
//! jump-GCD trajectories, record the per-step working-width EXTREME (no
//! additive margin -- the bound IS the observed extreme), set the iters budget
//! from a deep convergence quantile, validate on a held-out batch, and emit
//! `src/arith/schrottenloher/jump_schedule.rs` (committed; forward and reverse
//! GCD read the same arrays so they stay in lockstep).
//!
//! For each step i, current_n[i] = max over samples of max(bitlen(u),bitlen(v))
//! at the START of step i. The quantum register at step i holds current_n[i]
//! bits; a held-out shot "misses" if its width exceeds current_n[i] at some
//! step (truncation) or it has not converged by the budget.
//!
//! Env: JS_TRAIN (sample count/jump, default 30_000_000),
//!      JS_VALID (held-out count/jump, default 5_000_000),
//!      JS_JUMPS (csv, default "1,2,3,4").
//!
//! Run: cargo run --release --bin gen_jump_schedule

use alloy_primitives::U512;
use rand::RngCore;
use std::io::Write;

const N_STEPS: usize = 460;
const N_W: usize = 257; // width bins 0..=256
const CONV_BUDGET: f64 = 1e-7; // iters = 1 - 1e-7 convergence quantile + buffer

#[inline]
fn bl(x: &U512) -> usize {
    512 - x.leading_zeros()
}

fn secp_q() -> U512 {
    // q = 2^256 - 2^32 - 977 = F_SECP256K1 modulus the GCD inverts against.
    (U512::from(1u64) << 256) - (U512::from(1u64) << 32) - U512::from(977u64)
}

fn rand_x(rng: &mut impl RngCore, q: &U512) -> U512 {
    loop {
        let mut limbs = [0u64; 4];
        for l in limbs.iter_mut() {
            *l = rng.next_u64();
        }
        let r = U512::from(limbs[0])
            ^ (U512::from(limbs[1]) << 64)
            ^ (U512::from(limbs[2]) << 128)
            ^ (U512::from(limbs[3]) << 192);
        if r != U512::ZERO && r < *q {
            return r;
        }
    }
}

/// One jump-GCD trajectory. Updates per-step width extremes; returns the
/// convergence step (or N_STEPS if it did not converge). Mirrors the step in
/// `gcd_jump.rs` / `to_bitvector_classical_jump` exactly.
fn run_sample(x: U512, q: &U512, jump: usize, hist: &mut [u32]) -> usize {
    let one = U512::from(1u64);
    let (mut u, mut v) = (*q, x);
    let mut conv = N_STEPS;
    for step in 0..N_STEPS {
        let w = bl(&u).max(bl(&v));
        hist[step * N_W + w] += 1;
        // jump-before-swap: shift v to odd (up to jump) FIRST, then subtract.
        let mut j = 0;
        while j < jump && v != U512::ZERO && !v.bit(0) {
            v >>= 1;
            j += 1;
        }
        if v.bit(0) {
            if u > v {
                std::mem::swap(&mut u, &mut v);
            }
            v -= u;
        }
        if u == one && v == U512::ZERO {
            conv = step + 1;
            break;
        }
    }
    conv
}

/// Held-out check: does this trajectory FIT the schedule? (width never exceeds
/// sched[step] within `iters`, and it converges by `iters`).
fn fits(x: U512, q: &U512, jump: usize, sched: &[u16], iters: usize) -> bool {
    let one = U512::from(1u64);
    let (mut u, mut v) = (*q, x);
    for step in 0..iters {
        if (bl(&u).max(bl(&v)) as u16) > sched[step] {
            return false;
        }
        let mut j = 0;
        while j < jump && v != U512::ZERO && !v.bit(0) {
            v >>= 1;
            j += 1;
        }
        if v.bit(0) {
            if u > v {
                std::mem::swap(&mut u, &mut v);
            }
            v -= u;
        }
        if u == one && v == U512::ZERO {
            return true;
        }
    }
    false
}

/// Max comparator gap over a forward GCD pass: at each swap-decision step the
/// narrow comparator scans the top `trunc` bits of width `sched[step]`; the gap
/// is `sched[step] - bitlen(u^v)` (the differing bit's depth below the schedule
/// top). A top-`trunc` comparator mis-decides the pass iff this max gap >= trunc,
/// so one pass yields the miss rate for EVERY trunc.
fn cmp_max_gap(x: U512, q: &U512, jump: usize, sched: &[u16], iters: usize) -> usize {
    let one = U512::from(1u64);
    let (mut u, mut v) = (*q, x);
    let mut maxgap = 0usize;
    for step in 0..iters {
        let cn = sched[step] as usize;
        let mut j = 0;
        while j < jump && v != U512::ZERO && !v.bit(0) {
            v >>= 1;
            j += 1;
        }
        if v.bit(0) {
            if u != v {
                let gap = cn.saturating_sub(bl(&(u ^ v)));
                if gap > maxgap {
                    maxgap = gap;
                }
            }
            if u > v {
                std::mem::swap(&mut u, &mut v);
            }
            v -= u;
        }
        if u == one && v == U512::ZERO {
            break;
        }
    }
    maxgap
}

struct JumpResult {
    jump: usize,
    sched: Vec<u16>, // length == iters
    iters: usize,
    sum_cn: u64,
    peak_cn: u16,
    misses: u64,
    valid: u64,
}

fn gen_jump(jump: usize, train: usize, valid: usize, q: &U512) -> JumpResult {
    let nthreads = std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(4);
    let per = train / nthreads;
    let width_tail: f64 = std::env::var("JS_WTAIL")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(1e-5);
    // --- training: per-step width HISTOGRAM + convergence histogram.
    let parts: Vec<(Vec<u32>, Vec<u64>)> = std::thread::scope(|sc| {
        let handles: Vec<_> = (0..nthreads)
            .map(|t| {
                let cnt = if t == nthreads - 1 {
                    train - per * (nthreads - 1)
                } else {
                    per
                };
                sc.spawn(move || {
                    let mut rng = rand::thread_rng();
                    let mut h = vec![0u32; N_STEPS * N_W];
                    let mut ch = vec![0u64; N_STEPS + 1];
                    for _ in 0..cnt {
                        let x = rand_x(&mut rng, q);
                        let c = run_sample(x, q, jump, &mut h);
                        ch[c] += 1;
                    }
                    (h, ch)
                })
            })
            .collect();
        handles.into_iter().map(|h| h.join().unwrap()).collect()
    });
    let mut hist = vec![0u64; N_STEPS * N_W];
    let mut conv_hist = vec![0u64; N_STEPS + 1];
    for (h, ch) in &parts {
        for i in 0..N_STEPS * N_W {
            hist[i] += h[i] as u64;
        }
        for c in 0..=N_STEPS {
            conv_hist[c] += ch[c];
        }
    }
    // iters = deepest convergence step with tail prob <= CONV_BUDGET, + buffer.
    let total: u64 = conv_hist.iter().sum();
    let thr = (CONV_BUDGET * total as f64).ceil() as u64;
    let mut tail = 0u64;
    let mut n_used = N_STEPS;
    for nu in (1..N_STEPS).rev() {
        tail += conv_hist[nu + 1];
        if tail > thr {
            n_used = nu + 1;
            break;
        }
        n_used = nu;
    }
    // buffer above the quantile, rounded UP to a multiple of 3 so the base-5
    // packer (3 symbols/window) tiles the dialog with no partial window.
    let iters = ((n_used + 4).div_ceil(3) * 3).min(N_STEPS);
    // per-step width = the (1 - width_tail) quantile of max(bitlen(u),bitlen(v))
    // over the shots that REACH that step (truncating the slowest width_tail).
    let sched: Vec<u16> = (0..iters)
        .map(|step| {
            let row = &hist[step * N_W..step * N_W + N_W];
            let tot: u64 = row.iter().sum();
            let keep = ((1.0 - width_tail) * tot as f64).ceil() as u64;
            let mut cum = 0u64;
            let mut w = N_W - 1;
            for (width, &c) in row.iter().enumerate() {
                cum += c;
                if cum >= keep {
                    w = width;
                    break;
                }
            }
            w as u16
        })
        .collect();
    let sum_cn: u64 = sched.iter().map(|&c| c as u64).sum();
    let peak_cn = *sched.iter().max().unwrap();

    // --- held-out validation: whole-pass miss rate against this schedule.
    let vper = valid / nthreads;
    let miss_parts: Vec<u64> = std::thread::scope(|sc| {
        let sref = &sched;
        let handles: Vec<_> = (0..nthreads)
            .map(|t| {
                let cnt = if t == nthreads - 1 {
                    valid - vper * (nthreads - 1)
                } else {
                    vper
                };
                sc.spawn(move || {
                    let mut rng = rand::thread_rng();
                    let mut miss = 0u64;
                    for _ in 0..cnt {
                        let x = rand_x(&mut rng, q);
                        if !fits(x, q, jump, sref, iters) {
                            miss += 1;
                        }
                    }
                    miss
                })
            })
            .collect();
        handles.into_iter().map(|h| h.join().unwrap()).collect()
    });
    let misses: u64 = miss_parts.iter().sum();
    eprintln!(
        "[gen] jump={jump}: iters={iters} sum(cn)={sum_cn} peak={peak_cn} \
         held-out miss={misses}/{valid} (rel9k~{:.5})",
        (1.0 - misses as f64 / valid as f64).powi(9000)
    );

    // Comparator-gap reliability: histogram each pass's MAX comparator gap, then
    // report P(max-gap >= trunc) for a trunc sweep. A top-`trunc` swap comparator
    // mis-decides a GCD pass iff its max gap >= trunc; a point-add runs ~4 such
    // passes and the fuzz bar is 9000 Fiat-Shamir points (=> ~36000 exponent).
    const GMAX: usize = 200;
    let gap_parts: Vec<Vec<u64>> = std::thread::scope(|sc| {
        let sref = &sched;
        let handles: Vec<_> = (0..nthreads)
            .map(|t| {
                let cnt = if t == nthreads - 1 {
                    valid - vper * (nthreads - 1)
                } else {
                    vper
                };
                sc.spawn(move || {
                    let mut rng = rand::thread_rng();
                    let mut h = vec![0u64; GMAX + 1];
                    for _ in 0..cnt {
                        let x = rand_x(&mut rng, q);
                        h[cmp_max_gap(x, q, jump, sref, iters).min(GMAX)] += 1;
                    }
                    h
                })
            })
            .collect();
        handles.into_iter().map(|h| h.join().unwrap()).collect()
    });
    let mut ghist = vec![0u64; GMAX + 1];
    for h in &gap_parts {
        for g in 0..=GMAX {
            ghist[g] += h[g];
        }
    }
    for &trunc in &[40usize, 48, 56, 64, 72, 80, 88, 96] {
        let tail: u64 = ghist[trunc..=GMAX].iter().sum();
        let p = tail as f64 / valid as f64;
        eprintln!(
            "[cmp] jump={jump} trunc={trunc}: pass-miss={tail}/{valid} (p={p:.3e}) \
             rel(4pass*9k)~{:.5}",
            (1.0 - p).powi(9000 * 4)
        );
    }
    let maxg = (0..=GMAX).rev().find(|&g| ghist[g] > 0).unwrap_or(0);
    eprintln!("[cmp] jump={jump}: observed max comparator gap = {maxg}");

    JumpResult {
        jump,
        sched,
        iters,
        sum_cn,
        peak_cn,
        misses,
        valid: valid as u64,
    }
}

fn main() {
    let train: usize = std::env::var("JS_TRAIN")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(30_000_000);
    let valid: usize = std::env::var("JS_VALID")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(5_000_000);
    let jumps: Vec<usize> = std::env::var("JS_JUMPS")
        .unwrap_or_else(|_| "1,2,3,4".into())
        .split(',')
        .filter_map(|s| s.trim().parse().ok())
        .collect();
    let q = secp_q();
    eprintln!("[gen] train={train}/jump valid={valid}/jump jumps={jumps:?}");

    let results: Vec<JumpResult> = jumps
        .iter()
        .map(|&j| gen_jump(j, train, valid, &q))
        .collect();

    let mut out = String::new();
    out.push_str("//! GENERATED by src/bin/gen_jump_schedule.rs -- do not edit by hand.\n");
    out.push_str("//!\n");
    out.push_str("//! Per-step register-shrink schedule for the jump-GCD: current_n[i] is the\n");
    out.push_str("//! observed max(bitlen(u),bitlen(v)) extreme at step i (no additive margin).\n");
    out.push_str(
        "//! Forward and reverse GCD shrink u,v to current_n[i] bits at step i and read\n",
    );
    out.push_str(&format!(
        "//! THESE arrays in lockstep. train={train}/jump, valid={valid}/jump.\n//!\n"
    ));
    for r in &results {
        out.push_str(&format!(
            "//! jump={}: iters={} sum(cn)={} peak={} held-out miss={}/{} (rel9k~{:.5})\n",
            r.jump,
            r.iters,
            r.sum_cn,
            r.peak_cn,
            r.misses,
            r.valid,
            (1.0 - r.misses as f64 / r.valid as f64).powi(9000)
        ));
    }
    out.push('\n');
    out.push_str("/// (schedule, iters_budget) for `jump`; iters_budget == schedule.len().\n");
    out.push_str("pub fn jump_schedule(jump: usize) -> (&'static [u16], usize) {\n");
    out.push_str("    let s: &'static [u16] = match jump {\n");
    for r in &results {
        out.push_str(&format!("        {} => SCHED_J{},\n", r.jump, r.jump));
    }
    out.push_str("        _ => panic!(\"jump_schedule: no schedule for jump={jump}\"),\n");
    out.push_str("    };\n    (s, s.len())\n}\n\n");
    for r in &results {
        let body: Vec<String> = r.sched.iter().map(|c| c.to_string()).collect();
        out.push_str(&format!(
            "static SCHED_J{}: &[u16] = &[{}];\n\n",
            r.jump,
            body.join(", ")
        ));
    }

    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("src/arith/schrottenloher/jump_schedule.rs");
    let mut f = std::fs::File::create(&path).expect("create jump_schedule.rs");
    f.write_all(out.as_bytes()).expect("write");
    eprintln!("[gen] wrote {}", path.display());
}
