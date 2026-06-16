//! FAST (Rust) generator for the shrunken-PZ per-step schedule. Ports the Python
//! scripts/shrunken_pz_schedule_gen.py run_widths sampler to U512 so it can sample
//! enough trajectories (tens of millions) to estimate the deep per-step extremes
//! WITHOUT an additive margin -- the bounds ARE the observed extremes, and the
//! whole-pass tail is driven down by SAMPLE COUNT (a better quantile), not margin.
//!
//! Schedule bounds (margin=0):
//!   _W  (width)      = max bitlen each register HOLDS (incl. B<<s / cb<<s2 transients)
//!   _LO (window low) = (min value-bitlen) - 1  (clz/comparator scan src[LO..W])
//!   _SD (shift)      = rb_round(max shift)      (rotator distance ceiling, 2^k-1)
//!   NSTEPS           = max convergence + buffer (step count is tof-only/free)
//!
//! Env: SPZ_TRAIN (sample count, default 50_000_000), SPZ_VALID (held-out whole-pass
//! validation, default 5_000_000), SPZ_NOWRITE=1 (measure only). Writes
//! src/inversion/shrunken_pz_schedule.rs.

use alloy_primitives::U512;
use rand::RngCore;
use std::io::Write;

const N_STEPS: usize = 640;

#[inline]
fn bl(x: &U512) -> usize {
    512 - x.leading_zeros()
}
#[inline]
fn blq(q: u128) -> usize {
    128 - q.leading_zeros() as usize
}

fn secp_p() -> U512 {
    // 2^256 - 2^32 - 977
    (U512::from(1u64) << 256) - (U512::from(1u64) << 32) - U512::from(977u64)
}

fn rand_x(rng: &mut impl RngCore, p: &U512) -> U512 {
    loop {
        let mut limbs = [0u64; 4];
        for l in limbs.iter_mut() {
            *l = rng.next_u64();
        }
        let r = U512::from(limbs[0])
            ^ (U512::from(limbs[1]) << 64)
            ^ (U512::from(limbs[2]) << 128)
            ^ (U512::from(limbs[3]) << 192);
        if r != U512::ZERO && r < *p {
            return r;
        }
    }
}

/// One trajectory. Updates the per-step accumulators (maxw / minvbl / maxsd). Returns
/// the convergence step (or N_STEPS if it didn't converge in `n_steps`).
fn run_sample(
    x_orig: &U512,
    p: &U512,
    half: &U512,
    n_steps: usize,
    maxw: &mut [[u16; 5]],
    minvbl: &mut [[u16; 5]],
    maxsd: &mut [[i32; 2]],
) -> usize {
    let one = U512::from(1u64);
    let sgn = *x_orig > *half;
    let xx = if sgn { *p - *x_orig } else { *x_orig };
    let mut a = *p;
    let mut b = xx;
    let mut ca = U512::ZERO;
    let mut cb = one;
    let mut q: u128 = 0;
    let mut conv = n_steps;
    let mut converged = false;
    for step in 0..n_steps {
        if a.is_zero() && b == one && q == 0 {
            if !converged {
                conv = step;
                converged = true;
            }
            let bca = bl(&ca) as u16;
            let mw = &mut maxw[step];
            mw[0] = mw[0].max(0);
            mw[1] = mw[1].max(1);
            mw[2] = mw[2].max(bca);
            let mv = &mut minvbl[step];
            mv[1] = mv[1].min(1);
            mv[2] = if mv[2] == 0 { bca } else { mv[2].min(bca) };
            continue;
        }
        let (mut wa, mut wb, mut wca, mut wcb, mut wq) = (
            bl(&a) as u16,
            bl(&b) as u16,
            bl(&ca) as u16,
            bl(&cb) as u16,
            blq(q) as u16,
        );
        // value-bitlens at step start (clz/comparators read these) -> min window lo.
        {
            let mv = &mut minvbl[step];
            let vb = [
                bl(&a) as u16,
                bl(&b) as u16,
                bl(&ca) as u16,
                bl(&cb) as u16,
                blq(q) as u16,
            ];
            for r in 0..5 {
                if vb[r] != 0 {
                    mv[r] = if mv[r] == 0 { vb[r] } else { mv[r].min(vb[r]) };
                }
            }
        }
        let (mut s_div, mut s2_mul): (i32, i32) = (-1, -1);
        // MULTIPLY (A<B): cb<<s2 transient, ca grows. q!=0 here (bits to consume).
        if a < b && q != 0 {
            let s2 = q.trailing_zeros() as usize;
            s2_mul = s2 as i32;
            let cbs = cb << s2;
            wcb = wcb.max(bl(&cbs) as u16);
            q ^= 1u128 << s2;
            ca += cbs;
            wca = wca.max(bl(&ca) as u16);
            wq = wq.max(blq(q) as u16);
        }
        // DIVISION (ca<cb): B<<s transient, A shrinks.
        if ca < cb {
            let mut s: i64 = bl(&a) as i64 - bl(&b) as i64;
            let offset = if s >= 0 {
                a < (b << (s as usize))
            } else {
                false
            };
            if offset {
                s -= 1;
            }
            if s >= 0 {
                s_div = s as i32;
                let bsh = b << (s as usize);
                wb = wb.max(bl(&bsh) as u16);
                if a >= bsh {
                    a -= bsh;
                    q ^= 1u128 << s;
                    wq = wq.max(blq(q) as u16);
                }
                wa = wa.max(bl(&a) as u16);
            }
        }
        // SWAP
        if q == 0 && !a.is_zero() {
            std::mem::swap(&mut a, &mut b);
            std::mem::swap(&mut ca, &mut cb);
        }
        let mw = &mut maxw[step];
        let w = [wa, wb, wca, wcb, wq];
        for r in 0..5 {
            mw[r] = mw[r].max(w[r]);
        }
        let ms = &mut maxsd[step];
        ms[0] = ms[0].max(s_div);
        ms[1] = ms[1].max(s2_mul);
    }
    conv
}

fn rb_round(v: i32) -> u16 {
    let v = v.max(1) as u32;
    ((1u32 << (32 - v.leading_zeros())) - 1) as u16
}

/// Whole-pass fit check against a built schedule (sw,lo,sd) over n_used steps.
/// Returns None if fits, else the FIRST failing cause.
fn fits_why(
    x_orig: &U512,
    p: &U512,
    half: &U512,
    n_used: usize,
    sw: &[[u16; 5]],
    lo: &[[u16; 5]],
    sd: &[[u16; 2]],
) -> Option<&'static str> {
    let one = U512::from(1u64);
    let sgn = *x_orig > *half;
    let xx = if sgn { *p - *x_orig } else { *x_orig };
    let mut a = *p;
    let mut b = xx;
    let mut ca = U512::ZERO;
    let mut cb = one;
    let mut q: u128 = 0;
    for step in 0..n_used {
        if a.is_zero() && b == one && q == 0 {
            return None;
        }
        let vb = [
            bl(&a) as u16,
            bl(&b) as u16,
            bl(&ca) as u16,
            bl(&cb) as u16,
            blq(q) as u16,
        ];
        let (mut wa, mut wb, mut wca, mut wcb, mut wq) = (vb[0], vb[1], vb[2], vb[3], vb[4]);
        for r in 0..5 {
            if vb[r] != 0 && vb[r] <= lo[step][r] {
                return Some("window");
            }
        }
        let (mut s_div, mut s2_mul): (i32, i32) = (-1, -1);
        if a < b && q != 0 {
            let s2 = q.trailing_zeros() as usize;
            s2_mul = s2 as i32;
            let cbs = cb << s2;
            wcb = wcb.max(bl(&cbs) as u16);
            q ^= 1u128 << s2;
            ca += cbs;
            wca = wca.max(bl(&ca) as u16);
            wq = wq.max(blq(q) as u16);
        }
        if ca < cb {
            let mut s: i64 = bl(&a) as i64 - bl(&b) as i64;
            let offset = if s >= 0 {
                a < (b << (s as usize))
            } else {
                false
            };
            if offset {
                s -= 1;
            }
            if s >= 0 {
                s_div = s as i32;
                let bsh = b << (s as usize);
                wb = wb.max(bl(&bsh) as u16);
                if a >= bsh {
                    a -= bsh;
                    q ^= 1u128 << s;
                    wq = wq.max(blq(q) as u16);
                }
                wa = wa.max(bl(&a) as u16);
            }
        }
        if q == 0 && !a.is_zero() {
            std::mem::swap(&mut a, &mut b);
            std::mem::swap(&mut ca, &mut cb);
        }
        let w = [wa, wb, wca, wcb, wq];
        for r in 0..5 {
            if w[r] > sw[step][r] {
                return Some("width");
            }
        }
        if s_div > sd[step][0] as i32 || s2_mul > sd[step][1] as i32 {
            return Some("shift");
        }
    }
    Some("noconv")
}

/// Max comparator differing-bit gap over a pass, per compare type:
///   [0] a<b  (multiply branch selector), anchored at max(sw[A],sw[B])
///   [1] ca<cb (division branch selector), anchored at max(sw[ca],sw[cb])
///   [2] a<(b<<s) (the aligned divstep decision), anchored at sw[A]
/// gap = W - bitlen(x^y); a top-`trunc` comparator window mis-decides that
/// compare iff its gap >= trunc. Returns the per-type max over the pass.
/// Per compare type, returns [sched-anchored gap; msb-anchored gap]. sched gap =
/// W_schedule - diffbit (fixed top window from the register width); msb gap =
/// max(bitlen) - diffbit (window anchored at the operands' runtime leading bit,
/// i.e. how far below the MSB they differ -- the GCD-closeness term).
fn cmp_max_gaps(
    x_orig: &U512,
    p: &U512,
    half: &U512,
    n_used: usize,
    sw: &[[u16; 5]],
) -> ([usize; 3], [usize; 3]) {
    let one = U512::from(1u64);
    let sgn = *x_orig > *half;
    let xx = if sgn { *p - *x_orig } else { *x_orig };
    let (mut a, mut b, mut ca, mut cb) = (*p, xx, U512::ZERO, one);
    let mut q: u128 = 0;
    let mut gs = [0usize; 3]; // schedule-anchored
    let mut gm = [0usize; 3]; // msb-anchored
    let mut rec = |i: usize, x: &U512, y: &U512, wsched: usize| {
        let d = bl(&(*x ^ *y));
        if d != 0 {
            gs[i] = gs[i].max(wsched.saturating_sub(d));
            let msb = bl(x).max(bl(y));
            gm[i] = gm[i].max(msb.saturating_sub(d));
        }
    };
    for step in 0..n_used {
        if a.is_zero() && b == one && q == 0 {
            break;
        }
        if q != 0 {
            rec(0, &a, &b, sw[step][0].max(sw[step][1]) as usize);
        }
        if a < b && q != 0 {
            let s2 = q.trailing_zeros() as usize;
            ca += cb << s2;
            q ^= 1u128 << s2;
        }
        rec(1, &ca, &cb, sw[step][2].max(sw[step][3]) as usize);
        if ca < cb {
            let mut s: i64 = bl(&a) as i64 - bl(&b) as i64;
            if s >= 0 {
                let bsh = b << (s as usize);
                rec(2, &a, &bsh, sw[step][0] as usize);
                if a < bsh {
                    s -= 1;
                }
            }
            if s >= 0 {
                let bsh = b << (s as usize);
                if a >= bsh {
                    a -= bsh;
                    q ^= 1u128 << s;
                }
            }
        }
        if q == 0 && !a.is_zero() {
            std::mem::swap(&mut a, &mut b);
            std::mem::swap(&mut ca, &mut cb);
        }
    }
    (gs, gm)
}

fn main() {
    let train: usize = std::env::var("SPZ_TRAIN")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(50_000_000);
    let valid: usize = std::env::var("SPZ_VALID")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(5_000_000);
    let p = secp_p();
    let half = p >> 1;
    let nthreads = std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(4);
    let per = train / nthreads;
    eprintln!("[gen] sampling {train} trajectories across {nthreads} threads...");
    let parts: Vec<(Vec<[u16; 5]>, Vec<[u16; 5]>, Vec<[i32; 2]>, Vec<u64>)> =
        std::thread::scope(|sc| {
            let handles: Vec<_> = (0..nthreads)
                .map(|t| {
                    let cnt = if t == nthreads - 1 {
                        train - per * (nthreads - 1)
                    } else {
                        per
                    };
                    sc.spawn(move || {
                        let mut rng = rand::thread_rng();
                        let mut mw = vec![[0u16; 5]; N_STEPS];
                        let mut mv = vec![[0u16; 5]; N_STEPS];
                        let mut ms = vec![[-1i32; 2]; N_STEPS];
                        let mut ch = vec![0u64; N_STEPS + 1];
                        for _ in 0..cnt {
                            let x = rand_x(&mut rng, &p);
                            let c = run_sample(&x, &p, &half, N_STEPS, &mut mw, &mut mv, &mut ms);
                            ch[c] += 1;
                        }
                        (mw, mv, ms, ch)
                    })
                })
                .collect();
            handles.into_iter().map(|h| h.join().unwrap()).collect()
        });
    // reduce across threads.
    let mut maxw = vec![[0u16; 5]; N_STEPS];
    let mut minvbl = vec![[0u16; 5]; N_STEPS];
    let mut maxsd = vec![[-1i32; 2]; N_STEPS];
    let mut conv_hist = vec![0u64; N_STEPS + 1];
    for (mw, mv, ms, ch) in &parts {
        for i in 0..N_STEPS {
            for r in 0..5 {
                maxw[i][r] = maxw[i][r].max(mw[i][r]);
                if mv[i][r] != 0 {
                    minvbl[i][r] = if minvbl[i][r] == 0 {
                        mv[i][r]
                    } else {
                        minvbl[i][r].min(mv[i][r])
                    };
                }
            }
            maxsd[i][0] = maxsd[i][0].max(ms[i][0]);
            maxsd[i][1] = maxsd[i][1].max(ms[i][1]);
        }
        for c in 0..=N_STEPS {
            conv_hist[c] += ch[c];
        }
    }
    // n_used = smallest k with P(conv > k) <= CONV_BUDGET (step count is tof-only/free,
    // so a deep convergence quantile; NOT the noisy single-sample max).
    let conv_budget = 1e-7_f64;
    let total: u64 = conv_hist.iter().sum();
    let thr = (conv_budget * total as f64).ceil() as u64;
    let max_conv = (0..=N_STEPS).rev().find(|&c| conv_hist[c] > 0).unwrap_or(0);
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
    if conv_hist[N_STEPS] > thr {
        eprintln!(
            "[gen] WARNING: {} trajectories did not converge in N_STEPS={N_STEPS} (> budget {thr}); raise N_STEPS",
            conv_hist[N_STEPS]
        );
    }
    eprintln!("[gen] max_conv={max_conv}, conv 1-{conv_budget:e} quantile -> n_used={n_used}");

    // Build schedule (NO additive margin): bounds ARE the observed extremes.
    let sw: Vec<[u16; 5]> = maxw.clone();
    let lo: Vec<[u16; 5]> = (0..N_STEPS)
        .map(|i| {
            let mut r = [0u16; 5];
            for j in 0..5 {
                r[j] = minvbl[i][j].saturating_sub(1);
            }
            r
        })
        .collect();
    let sd: Vec<[u16; 2]> = (0..N_STEPS)
        .map(|i| [rb_round(maxsd[i][0]), rb_round(maxsd[i][1])])
        .collect();

    let peak_step = (0..n_used)
        .max_by_key(|&i| sw[i].iter().map(|&w| w as u32).sum::<u32>())
        .unwrap();
    let peak: u32 = sw[peak_step].iter().map(|&w| w as u32).sum();
    eprintln!(
        "[gen] PEAK A+B+ca+cb+q={peak} at step {peak_step}: A={} B={} ca={} cb={} q={}",
        sw[peak_step][0], sw[peak_step][1], sw[peak_step][2], sw[peak_step][3], sw[peak_step][4]
    );

    // Comparator differing-bit gap sweep. Per compare type, SCHEDULE-anchored
    // (fixed top window from the register width) and MSB-anchored (window from
    // the operands' runtime leading bit) max gaps + the smallest top-`trunc`
    // keeping the whole-pass mis-decide tail < 1e-7.
    {
        const GMAX: usize = 260;
        let gper = valid / nthreads;
        let gparts: Vec<[Vec<u64>; 6]> = std::thread::scope(|sc| {
            let swref = &sw;
            let handles: Vec<_> = (0..nthreads)
                .map(|t| {
                    let cnt = if t == nthreads - 1 {
                        valid - gper * (nthreads - 1)
                    } else {
                        gper
                    };
                    sc.spawn(move || {
                        let mut rng = rand::thread_rng();
                        let mut h: [Vec<u64>; 6] = std::array::from_fn(|_| vec![0u64; GMAX + 1]);
                        for _ in 0..cnt {
                            let x = rand_x(&mut rng, &p);
                            let (gs, gm) = cmp_max_gaps(&x, &p, &half, n_used, swref);
                            for ti in 0..3 {
                                h[ti][gs[ti].min(GMAX)] += 1;
                                h[3 + ti][gm[ti].min(GMAX)] += 1;
                            }
                        }
                        h
                    })
                })
                .collect();
            handles.into_iter().map(|h| h.join().unwrap()).collect()
        });
        let mut gh: [Vec<u64>; 6] = std::array::from_fn(|_| vec![0u64; GMAX + 1]);
        for h in &gparts {
            for ti in 0..6 {
                for k in 0..=GMAX {
                    gh[ti][k] += h[ti][k];
                }
            }
        }
        let names = ["a<b", "ca<cb", "a<(b<<s)"];
        for ti in 0..6 {
            let total: u64 = gh[ti].iter().sum();
            let maxg = (0..=GMAX).rev().find(|&k| gh[ti][k] > 0).unwrap_or(0);
            let mut tail = 0u64;
            let mut safe = 0usize;
            for tr in (0..=GMAX).rev() {
                tail += gh[ti][tr];
                if (tail as f64) > 1e-7 * total as f64 {
                    safe = tr + 1;
                    break;
                }
            }
            let anchor = if ti < 3 { "sched" } else { "msb " };
            eprintln!(
                "[cmp] {:>9} ({anchor}): maxgap={maxg} safe_trunc={safe}",
                names[ti % 3]
            );
        }
    }

    // Validate whole-pass on a held-out set (threaded).
    eprintln!("[gen] validating {valid} held-out trajectories across {nthreads} threads...");
    let vper = valid / nthreads;
    let vparts: Vec<(usize, std::collections::BTreeMap<&'static str, usize>)> =
        std::thread::scope(|sc| {
            let handles: Vec<_> = (0..nthreads)
                .map(|t| {
                    let cnt = if t == nthreads - 1 {
                        valid - vper * (nthreads - 1)
                    } else {
                        vper
                    };
                    let (swr, lor, sdr) = (&sw, &lo, &sd);
                    sc.spawn(move || {
                        let mut rng = rand::thread_rng();
                        let mut m = 0usize;
                        let mut w: std::collections::BTreeMap<&'static str, usize> =
                            Default::default();
                        for _ in 0..cnt {
                            let x = rand_x(&mut rng, &p);
                            if let Some(c) = fits_why(&x, &p, &half, n_used, swr, lor, sdr) {
                                m += 1;
                                *w.entry(c).or_default() += 1;
                            }
                        }
                        (m, w)
                    })
                })
                .collect();
            handles.into_iter().map(|h| h.join().unwrap()).collect()
        });
    let mut miss = 0usize;
    let mut why: std::collections::BTreeMap<&'static str, usize> = Default::default();
    for (m, w) in vparts {
        miss += m;
        for (k, v) in w {
            *why.entry(k).or_default() += v;
        }
    }
    let rate = 1.0 - (miss as f64) / (valid as f64);
    eprintln!(
        "[gen] WHOLE-PASS (margin=0): {:.6}%  ({miss}/{valid} miss, causes={why:?})",
        rate * 100.0
    );
    let f = (miss as f64) / (valid as f64);
    // 95% over 9000 EC-adds (2 inversions each) => per-inversion f <= ~2.85e-6.
    let p9000 = (1.0 - 2.0 * f).max(0.0).powi(9000);
    eprintln!(
        "[gen] est P(all 9000 EC-adds ok) = {:.4} (target >= 0.95)",
        p9000
    );

    if std::env::var("SPZ_NOWRITE").ok().as_deref() == Some("1") {
        eprintln!("[gen] SPZ_NOWRITE=1: schedule NOT written");
        return;
    }

    // Write src/inversion/shrunken_pz_schedule.rs.
    let path = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/src/inversion/shrunken_pz_schedule.rs"
    );
    let mut s = String::new();
    s.push_str("// GENERATED by src/bin/gen_shrunken_pz_schedule.rs -- cross-gated shrunken-PZ\n");
    s.push_str(&format!(
        "// single-q inversion per-step schedule (margin=0; bounds = per-step extremes\n\
         // over {train} samples; held-out whole-pass {:.5}% = {miss}/{valid} miss).\n",
        rate * 100.0
    ));
    s.push_str(&format!(
        "// peak A+B+ca+cb+q={peak} at step {peak_step}.\n"
    ));
    s.push_str("// _W = register width (= max bitlen incl transient); _LO = clz window\n");
    s.push_str("// low bound (scan src[LO..W], MSB guaranteed >= LO); _SD = shift bound.\n\n");
    s.push_str(&format!(
        "#[allow(dead_code)]\npub const SHRUNKEN_PZ_NSTEPS: usize = {n_used};\n\n"
    ));
    let arr = |name: &str, vals: &dyn Fn(usize) -> u16| -> String {
        let body: Vec<String> = (0..n_used).map(|i| vals(i).to_string()).collect();
        format!(
            "#[allow(dead_code)]\npub const {name}: [u16; SHRUNKEN_PZ_NSTEPS] = [{}];\n",
            body.join(", ")
        )
    };
    for (nm, idx) in [
        ("SHRUNKEN_PZ_A", 0),
        ("SHRUNKEN_PZ_B", 1),
        ("SHRUNKEN_PZ_CA", 2),
        ("SHRUNKEN_PZ_CB", 3),
        ("SHRUNKEN_PZ_Q", 4),
    ] {
        s.push_str(&arr(nm, &|i| sw[i][idx]));
    }
    for (nm, idx) in [
        ("SHRUNKEN_PZ_A_LO", 0),
        ("SHRUNKEN_PZ_B_LO", 1),
        ("SHRUNKEN_PZ_CA_LO", 2),
        ("SHRUNKEN_PZ_CB_LO", 3),
        ("SHRUNKEN_PZ_Q_LO", 4),
    ] {
        s.push_str(&arr(nm, &|i| lo[i][idx]));
    }
    s.push_str(&arr("SHRUNKEN_PZ_SDIV", &|i| sd[i][0]));
    s.push_str(&arr("SHRUNKEN_PZ_S2", &|i| sd[i][1]));
    s.push_str(
        r#"
/// Per-step register widths (A, B, ca, cb, q). Out-of-range clamps to last.
#[allow(dead_code)]
pub fn reg_widths(i: usize) -> (usize, usize, usize, usize, usize) {
    let j = i.min(SHRUNKEN_PZ_NSTEPS - 1);
    (SHRUNKEN_PZ_A[j] as usize, SHRUNKEN_PZ_B[j] as usize, SHRUNKEN_PZ_CA[j] as usize,
     SHRUNKEN_PZ_CB[j] as usize, SHRUNKEN_PZ_Q[j] as usize)
}

/// Per-step clz-window low bounds (A, B, ca, cb, q): scan src[LO..W], the MSB
/// is guaranteed in [LO, W) for whole-pass-fitting inputs.
#[allow(dead_code)]
pub fn reg_los(i: usize) -> (usize, usize, usize, usize, usize) {
    let j = i.min(SHRUNKEN_PZ_NSTEPS - 1);
    (SHRUNKEN_PZ_A_LO[j] as usize, SHRUNKEN_PZ_B_LO[j] as usize, SHRUNKEN_PZ_CA_LO[j] as usize,
     SHRUNKEN_PZ_CB_LO[j] as usize, SHRUNKEN_PZ_Q_LO[j] as usize)
}

/// Per-step shift bounds (division s, multiply s2) -> rotator distance ceiling.
#[allow(dead_code)]
pub fn shift_bounds(i: usize) -> (usize, usize) {
    let j = i.min(SHRUNKEN_PZ_NSTEPS - 1);
    (SHRUNKEN_PZ_SDIV[j] as usize, SHRUNKEN_PZ_S2[j] as usize)
}
"#,
    );
    std::fs::File::create(path)
        .unwrap()
        .write_all(s.as_bytes())
        .unwrap();
    eprintln!("[gen] wrote {path}");
}
