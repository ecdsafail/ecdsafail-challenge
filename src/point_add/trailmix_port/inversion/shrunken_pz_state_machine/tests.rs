//! Tests for the shrunken-PZ reversible inversion state machine
//! (extracted from the former single-file module).

use super::*;
use num_bigint::BigUint;
use num_traits::{One, Zero};
use rand::Rng;

fn secp_p() -> BigUint {
    (BigUint::from(1u32) << 256u32) - (BigUint::from(1u32) << 32u32) - BigUint::from(977u32)
}
fn bl(x: &BigUint) -> i64 {
    x.bits() as i64
}

/// Classical predicate: does the divstep trajectory for this `(x, p)` input
/// stay within every per-step register width budget of the schedule? Used to
/// filter test shots for the windowed `shrunken_pz_divide_*` drivers.
fn shrunken_pz_shot_fits(x: &BigUint, p: &BigUint, n: usize) -> bool {
    use crate::point_add::trailmix_port::inversion::shrunken_pz_schedule::{reg_los, reg_widths, shift_bounds};
    let one = BigUint::one();
    let half = p >> 1u32;
    let sgn = x > &half;
    let mut a = p.clone();
    let mut b = if sgn { p - x } else { x.clone() };
    let (mut ca, mut cb) = (BigUint::zero(), one.clone());
    let mut q: u128 = 0;
    for i in 0..n {
        let (wa, wb, wca, wcb, wq) = reg_widths(i);
        let (lo_a, lo_b, ca_window, cb_window, _) = reg_los(i);
        let (sdb, s2b) = shift_bounds(i);
        let (wg, wc) = (wa.max(wb), wca.max(wcb));
        if a.is_zero() && b == one && q == 0 {
            return true;
        }
        let mut hg = a.bits().max(b.bits()) as usize;
        let mut hc = ca.bits().max(cb.bits()) as usize;
        let mut hq = (128 - (q | 1).leading_zeros()) as usize;
        if a < b {
            let s2 = q.trailing_zeros() as usize;
            if s2 > s2b || (cb.bits() as usize) <= cb_window {
                return false;
            }
            let cbs = &cb << s2;
            hc = hc.max(cbs.bits() as usize);
            q ^= 1u128 << s2;
            ca += &cbs;
            hc = hc.max(ca.bits() as usize);
            hq = hq.max((128 - (q | 1).leading_zeros()) as usize);
            if (ca.bits() as usize) <= ca_window || (cbs.bits() as usize) <= ca_window {
                return false;
            }
        }
        if ca < cb {
            if (a.bits() as usize) <= lo_a || (b.bits() as usize) <= lo_b {
                return false;
            }
            let mut s: i64 = a.bits() as i64 - b.bits() as i64;
            if s >= 0 && (s as usize) > sdb {
                return false;
            }
            let mut offset = false;
            if s >= 0 {
                offset = a < (&b << (s as u32));
            }
            if offset {
                s -= 1;
            }
            if s >= 0 {
                let bsh = &b << (s as u32);
                hg = hg.max(bsh.bits() as usize);
                if a >= bsh {
                    a -= &bsh;
                    q ^= 1u128 << (s as u32);
                    hq = hq.max((128 - (q | 1).leading_zeros()) as usize);
                }
                hg = hg.max(a.bits() as usize);
            }
        }
        if q == 0 && !a.is_zero() {
            std::mem::swap(&mut a, &mut b);
            std::mem::swap(&mut ca, &mut cb);
        }
        if hg > wg || hc > wc || hq > wq {
            return false;
        }
    }
    false
}

/// Faithful integer port of scripts/kaliski_test.py `pz_big_step` -- the
/// circuit's reference. Returns (inverse a, iteration count). The iteration
/// count drives the fixed-circuit length (n_iters = its p99.99 + the schedule).
/// Returns (inverse, iteration count, peak bitlen-sum bl(A)+bl(B)+bl(ca)+bl(cb)
/// over the iterations -- the unpacked SM register peak before pads/ancillae).
/// `wmax[k][i]` = running max bitlen of register k in {A,B,ca,cb} at iter i,
/// across the caller's shots. The UNPACKED register peak is
/// max_i (wmax[0][i]+wmax[1][i]+wmax[2][i]+wmax[3][i]) -- independent per
/// register, NOT the per-shot sum (that's the cursor/packed peak).
fn pz_big_step_ref(
    x_in: &BigUint,
    s_hist: &mut [u64; 16],
    wmax: &mut [[u16; 512]; 4],
) -> (BigUint, usize, i64) {
    let p = secp_p();
    let half = &p >> 1;
    let sgn = x_in > &half;
    let x = if sgn { &p - x_in } else { x_in.clone() };

    let (mut a_gcd, mut b_gcd) = (p.clone(), x.clone());
    let (mut ca, mut cb) = (BigUint::zero(), BigUint::one());
    let (mut q_div, mut q_mul): (u128, u128) = (0, 0);
    let mut mul_active = false;
    let mut div_active = true;
    let mut parity = true;
    let mut iters = 0usize;
    let mut peak: i64 = 0;

    while div_active || mul_active {
        peak = peak.max(bl(&a_gcd) + bl(&b_gcd) + bl(&ca) + bl(&cb));
        if iters < 512 {
            let row = [bl(&a_gcd), bl(&b_gcd), bl(&ca), bl(&cb)];
            for k in 0..4 {
                wmax[k][iters] = wmax[k][iters].max(row[k] as u16);
            }
        }
        iters += 1;
        // ---- division substep (runs every iter; self-gates via A<B) ----
        {
            let mut s = bl(&a_gcd) - bl(&b_gcd);
            if div_active && s >= 0 {
                s_hist[(s as usize).min(15)] += 1; // rotation amount distribution
            }
            let mut offset = false;
            if s >= 0 {
                b_gcd <<= s as usize; // align B<<s
                offset = a_gcd < b_gcd;
            }
            if offset {
                s -= 1;
                if s >= 0 {
                    b_gcd >>= 1u32;
                }
            }
            if s >= 0 && a_gcd >= b_gcd {
                a_gcd = &a_gcd - &b_gcd; // subtract
                b_gcd >>= s as usize; // restore B>>s
                q_div ^= 1u128 << s; // set quotient bit s
            } else {
                // can't subtract (A<B) -> division of this quotient is done
                div_active = false;
            }
        }
        // ---- transition (division finished, hand the quotient to multiply) ----
        if !div_active && !mul_active && q_div != 0 {
            std::mem::swap(&mut q_div, &mut q_mul);
            if !a_gcd.is_zero() {
                std::mem::swap(&mut a_gcd, &mut b_gcd);
            }
            mul_active = true;
            div_active = true;
        }
        // ---- multiply substep (pipelined) ----
        if mul_active {
            if q_mul != 0 {
                let mut s = q_mul.trailing_zeros();
                q_mul ^= 1u128 << s;
                cb <<= s; // b<<s
                ca = &ca + &cb; // a += b<<s
                if bl(&ca) != bl(&cb) {
                    cb <<= 1u32;
                    s += 1;
                }
                cb >>= s; // restore b (s incremented above if bitlens differed)
            } else {
                std::mem::swap(&mut ca, &mut cb);
                parity = !parity;
                mul_active = false;
            }
        }
        // ---- transition again (covers div-done landing after a mul step) ----
        if !div_active && !mul_active && q_div != 0 {
            std::mem::swap(&mut q_div, &mut q_mul);
            if !a_gcd.is_zero() {
                std::mem::swap(&mut a_gcd, &mut b_gcd);
            }
            mul_active = true;
            div_active = true;
        }
        assert!(iters < 100_000, "pz_big_step_ref did not terminate");
    }
    // end: A==0, B==1, b==P. sign flips are mod p (ca is a magnitude).
    let mut a_out = ca % &p;
    if parity {
        a_out = (&p - &a_out) % &p;
    }
    if sgn {
        a_out = (&p - &a_out) % &p;
    }
    (a_out, iters, peak)
}

/// Measure the codebase KG incrementer (consume-once) tof at n=256 to check
/// the paper's 3n claim and compare to bit_length_lean's per-sweep cost.
#[test]
fn prof_kg_incrementer() {
    let n = 256usize;
    let mut rng = rand::thread_rng();
    let mut c = Circuit::new();
    c.set_max_qubit_peak(400);
    let a = c.alloc_qreg_bits("a", n);
    for shot in 0..64 {
        load_big(
            &mut c,
            &a,
            &BigUint::from_bytes_le(&rng.gen::<[u8; 32]>()),
            shot,
        );
    }
    let aref: Vec<&QReg> = a.iter().collect();
    crate::point_add::trailmix_port::arith::khattar_gidney::inc_khattar_gidney_refs(&mut c, &aref);
    eprintln!(
        "inc_khattar_gidney n=256: peak {} q, tof {} (3n would be {})",
        c.peak_qubits,
        c.executed_toffoli_shots / 64,
        3 * n
    );
    let mut outs = vec![];
    outs.extend(a);
    let _ = c.destroy_sim(outs);
}

/// Lean prefix-AND bit_length: s += bitlen(src) on random values spanning
/// the full bitlen range (random value >> random shift). Checks value,
/// reversibility (dec restores s to 0 + src clean), and tof.
#[test]
fn bit_length_lean_correct() {
    let n = 256usize;
    let mut rng = rand::thread_rng();
    let mut c = Circuit::new();
    c.set_max_qubit_peak(700);
    let src = c.alloc_qreg_bits("src", n);
    let s = c.alloc_qreg_bits("s", 9);
    let mut svals = Vec::new();
    for shot in 0..64 {
        // random value, random right-shift -> bitlen uniformly spans 1..=256.
        let v = loop {
            let sh = (rng.gen::<u32>() % 256) as usize;
            let x = BigUint::from_bytes_le(&rng.gen::<[u8; 32]>()) >> sh;
            if !x.is_zero() {
                break x;
            }
        };
        load_big(&mut c, &src, &v, shot);
        svals.push(v);
    }
    let srcref: Vec<&QReg> = src.iter().collect();
    bit_length_lean(&mut c, &srcref, &s, false); // s += bitlen(src)
    {
        let (sr, sv) = (&s, svals.clone());
        c.contract_check("bll", move |view, shot| {
            let got = rd_big(&view, sr, shot);
            let exp = BigUint::from(sv[shot].bits());
            if got != exp {
                return Err(format!("bitlen got {got} exp {exp}"));
            }
            Ok(())
        });
    }
    eprintln!(
        "bit_length_lean ok n=256: peak {} q, tof {} (per scan)",
        c.peak_qubits,
        c.executed_toffoli_shots / 64
    );
    for (sec, t) in &c.executed_toffoli_by_section {
        if sec.starts_with("bll.") {
            eprintln!("    {}: {} tof", sec, t / 64);
        }
    }
    bit_length_lean(&mut c, &srcref, &s, true); // s -= bitlen(src) -> 0
    for q in s {
        c.zero_and_free(q);
    }
    c.assert_phase_clean(); // HMR (clear_and) corrections must net phase-clean
    let mut outs = vec![];
    outs.extend(src);
    let _ = c.destroy_sim(outs);
}

fn rd_big(view: &crate::point_add::trailmix_port::circuit::ContractSimView, reg: &[QReg], shot: usize) -> BigUint {
    let mut x = BigUint::zero();
    for (j, q) in reg.iter().enumerate() {
        if view.contract_read_bit_shot(q, shot) {
            x |= BigUint::one() << j;
        }
    }
    x
}

/// Pipeline reference -- forward carries the pack (A,B,ca,cb), in-flight pads
/// (q_div,q_mul), and the in-flight flag mul_active (+ done counter in the real
/// driver). div_active is data (A>=B). Transition on (A<B) AND !mul_active AND
/// (q_div!=0); multiply drains q_mul bit-by-bit, finish (swap ca<->cb, parity,
/// mul_active=0) the iter AFTER the last bit (mul_active lag). State after `k`.
#[allow(clippy::type_complexity)]
fn pz_ref_at(
    x: &BigUint,
    k: usize,
) -> (BigUint, BigUint, BigUint, BigUint, u128, u128, bool, bool) {
    let p = secp_p();
    let (mut a, mut b) = (p.clone(), x.clone());
    let (mut ca, mut cb) = (BigUint::zero(), BigUint::one());
    let (mut q_div, mut q_mul): (u128, u128) = (0, 0);
    let mut mul_active = false;
    let mut parity = true;
    let bl = |v: &BigUint| -> i64 { v.bits() as i64 };
    let transition =
        |a: &mut BigUint, b: &mut BigUint, qd: &mut u128, qm: &mut u128, ma: &mut bool| {
            if *a < *b && !*ma && *qd != 0 {
                std::mem::swap(qd, qm);
                if !a.is_zero() {
                    std::mem::swap(a, b);
                }
                *ma = true;
            }
        };
    for _ in 0..k {
        {
            let mut s = bl(&a) - bl(&b);
            let mut offset = false;
            if s >= 0 {
                b <<= s as usize;
                offset = a < b;
            }
            if offset {
                s -= 1;
                if s >= 0 {
                    b >>= 1u32;
                }
            }
            if s >= 0 && a >= b {
                a = &a - &b;
                b >>= s as usize;
                q_div ^= 1u128 << s;
            } else if s >= 0 {
                b >>= s as usize;
            }
        }
        transition(&mut a, &mut b, &mut q_div, &mut q_mul, &mut mul_active);
        if mul_active {
            if q_mul != 0 {
                let mut s2 = q_mul.trailing_zeros();
                q_mul ^= 1u128 << s2;
                cb <<= s2;
                ca = &ca + &cb;
                if bl(&ca) != bl(&cb) {
                    cb <<= 1u32;
                    s2 += 1;
                }
                cb >>= s2;
            } else {
                std::mem::swap(&mut ca, &mut cb);
                parity = !parity;
                mul_active = false;
            }
        }
        transition(&mut a, &mut b, &mut q_div, &mut q_mul, &mut mul_active);
    }
    (a, b, ca, cb, q_div, q_mul, mul_active, parity)
}

/// Like pz_ref_at but runs `full` complete iters then `sub` sub-phases of the
/// next iter: 1=after division, 2=+transition1, 3=+multiply, 4=+transition2.
/// For phase-level localization of the quantum pipeline.
#[allow(clippy::type_complexity)]
fn pz_ref_sub(
    x: &BigUint,
    full: usize,
    sub: usize,
) -> (BigUint, BigUint, BigUint, BigUint, u128, u128, bool, bool) {
    let p = secp_p();
    let (mut a, mut b) = (p.clone(), x.clone());
    let (mut ca, mut cb) = (BigUint::zero(), BigUint::one());
    let (mut q_div, mut q_mul): (u128, u128) = (0, 0);
    let mut mul_active = false;
    let mut parity = true;
    let bl = |v: &BigUint| -> i64 { v.bits() as i64 };
    macro_rules! division {
        () => {{
            let mut s = bl(&a) - bl(&b);
            let mut offset = false;
            if s >= 0 {
                b <<= s as usize;
                offset = a < b;
            }
            if offset {
                s -= 1;
                if s >= 0 {
                    b >>= 1u32;
                }
            }
            if s >= 0 && a >= b {
                a = &a - &b;
                b >>= s as usize;
                q_div ^= 1u128 << s;
            } else if s >= 0 {
                b >>= s as usize;
            }
        }};
    }
    macro_rules! trans {
        () => {{
            if a < b && !mul_active && q_div != 0 {
                std::mem::swap(&mut q_div, &mut q_mul);
                if !a.is_zero() {
                    std::mem::swap(&mut a, &mut b);
                }
                mul_active = true;
            }
        }};
    }
    macro_rules! multiply {
        () => {{
            if mul_active {
                if q_mul != 0 {
                    let mut s2 = q_mul.trailing_zeros();
                    q_mul ^= 1u128 << s2;
                    cb <<= s2;
                    ca = &ca + &cb;
                    if bl(&ca) != bl(&cb) {
                        cb <<= 1u32;
                        s2 += 1;
                    }
                    cb >>= s2;
                } else {
                    std::mem::swap(&mut ca, &mut cb);
                    parity = !parity;
                    mul_active = false;
                }
            }
        }};
    }
    for _ in 0..full {
        division!();
        trans!();
        multiply!();
        trans!();
    }
    if sub >= 1 {
        division!();
    }
    if sub >= 2 {
        trans!();
    }
    if sub >= 3 {
        multiply!();
    }
    if sub >= 4 {
        trans!();
    }
    (a, b, ca, cb, q_div, q_mul, mul_active, parity)
}

/// Pure-integer self-check: the data-derived pipeline reference computes 1/x.
#[test]
fn pz_ref_at_computes_inverse() {
    let p = secp_p();
    let mut rng = rand::thread_rng();
    for _ in 0..2000 {
        let x = (BigUint::from_bytes_le(&rng.gen::<[u8; 32]>()) % (&p - BigUint::one()))
            + BigUint::one();
        // run to termination (A==0,B==1).
        let (a, b, ca, cb, qd, qm, ma, parity) = pz_ref_at(&x, 480);
        assert!(
            a.is_zero() && b == BigUint::one(),
            "not terminated: A={a} B={b}"
        );
        assert!(qd == 0 && qm == 0 && !ma, "pads/flag not clean");
        // inverse = parity ? P - ca : ca.
        let inv = if parity { (&p - &ca) % &p } else { ca.clone() };
        assert_eq!((&inv * &x) % &p, BigUint::one(), "wrong inverse for x={x}");
        let _ = cb;
    }
}

fn load_big(c: &mut Circuit, reg: &[QReg], v: &BigUint, shot: usize) {
    let mut bytes = v.to_bytes_le();
    bytes.resize((reg.len() + 7) / 8, 0);
    c.sim_load_reg_bytes_shot(reg, &bytes, shot);
}

/// WINDOWED division substep vs the integer model + tof vs full-width.
/// A,B in a known bitlen window [W-D, W] (A>=B), so lo/rot_bits are exact.
/// Validates A->A-B<<s_eff, q_div bit, scratch clean, and the tof saving.
#[test]
fn division_substep_windowed_correct() {
    use crate::point_add::trailmix_port::inversion::shrunken_pz_schedule::{
        reg_los, reg_widths, shift_bounds, SHRUNKEN_PZ_NSTEPS,
    };
    // Use a REAL schedule step's clz window (the schedule was calibrated on the
    // actual EEA-trajectory value distribution). Pick a mid-trajectory step
    // whose A/B windows are comfortably wide, so the OFFSET comparison
    // `A < B<<diff` -- a SAME-bitlen compare (after B<<diff, bitlen == bitlen(A))
    // whose discriminating bit lies below the shared MSB -- reliably falls inside
    // [lo, n). (A tight window like step 0's 8-bit one only works for the real
    // EEA values the schedule was trained on, not independent random A,B.)
    let step = (0..SHRUNKEN_PZ_NSTEPS)
        .find(|&i| {
            let (wa, wb, ..) = reg_widths(i);
            let (alo, blo, ..) = reg_los(i);
            wa.min(wb).saturating_sub(alo.max(blo)) >= 45
        })
        .expect("a wide-window schedule step exists");
    let (wa, wb, _, _, _) = reg_widths(step);
    let (alo, blo, _, _, _) = reg_los(step);
    let (sdb, _) = shift_bounds(step);
    let n = wa.max(wb);
    let rot_bits = if sdb == 0 {
        1
    } else {
        64 - (sdb as u64).leading_zeros() as usize
    };
    let hi = wa.min(wb); // keep both MSBs in their windows
    let lobnd = alo.max(blo);
    let mut rng = rand::thread_rng();
    let mut c = Circuit::new();
    c.set_max_qubit_peak(1900);
    let a = c.alloc_qreg_bits("A", n);
    let b = c.alloc_qreg_bits("B", n);
    let q_div = c.alloc_qreg_bits("qd", n);
    let s_rot = c.alloc_qreg_bits("srot", 9);
    let offset = c.alloc_qreg("off");
    let active = c.alloc_qreg("act");
    let mut avs = Vec::new();
    let mut bvs = Vec::new();
    for shot in 0..64 {
        // bitlen(A) near the top of the window (ample room below alo for the
        // offset's discriminating bit); bitlen(B) = bitlen(A) - shift, shift in
        // [1, sdb] with B's MSB kept >= blo. bitlen(A) > bitlen(B) => A > B.
        let bl_a = rng.gen_range((hi - 3)..=hi);
        let max_sh = sdb.min(bl_a - lobnd - 2).max(1);
        let sh = rng.gen_range(1..=max_sh);
        let bl_b = bl_a - sh;
        let mut mk = |bl: usize| -> BigUint {
            (BigUint::one() << (bl - 1))
                | (BigUint::from_bytes_le(&rng.gen::<[u8; 32]>()) % (BigUint::one() << (bl - 1)))
        };
        let av = mk(bl_a);
        let bv = mk(bl_b);
        for (reg, val) in [(&a, &av), (&b, &bv)] {
            load_big(&mut c, reg, val, shot);
        }
        avs.push(av);
        bvs.push(bv);
    }
    c.x(&active);

    let t0 = c.executed_toffoli_shots;
    division_substep_windowed(
        &mut c, &a, &b, &q_div, &s_rot, &offset, &active, alo, blo, rot_bits,
    );
    let win_tof = (c.executed_toffoli_shots - t0) / 64;

    {
        let (ar, qr, av2, bv2) = (&a, &q_div, avs.clone(), bvs.clone());
        c.contract_check("divw", move |view, shot| {
            let (av, bv) = (&av2[shot], &bv2[shot]);
            let mut sh = av.bits() as i64 - bv.bits() as i64;
            if sh >= 0 && *av < (bv << sh as usize) {
                sh -= 1;
            }
            let exp_a = av - (bv << sh as usize);
            let got_a = rd_big(&view, ar, shot);
            if got_a != exp_a {
                return Err(format!("A: got {got_a:x} exp {exp_a:x} (s={sh})"));
            }
            let got_q = rd_big(&view, qr, shot);
            if got_q != (BigUint::one() << sh as usize) {
                return Err(format!("q_div wrong (s={sh}): got {got_q:x}"));
            }
            Ok(())
        });
    }
    for q in s_rot.into_iter() {
        c.zero_and_free(q);
    }
    c.x(&active); // uncompute the active=1 setup
    for q in [offset, active] {
        c.zero_and_free(q);
    }
    eprintln!(
        "division_substep_windowed (sched step {step}, A-window {} bits, rot {rot_bits}): \
         {win_tof} tof",
        wa - alo,
    );
    c.assert_phase_clean();
    let mut outs = vec![];
    for r in [a, b, q_div] {
        outs.extend(r);
    }
    let _ = c.destroy_sim(outs);
}

/// ROUND-TRIP: division_substep_windowed then its inverse must restore the
/// initial state (A,B,q_div back to start; s/s_rot/offset/active clean).
/// Validates the backward-pass substep gate-for-gate.
#[test]
fn division_substep_windowed_roundtrip() {
    let n = 256usize;
    let (wlo, whi) = (248usize, 256usize);
    let lo = wlo - 1;
    let rot_bits = 4usize;
    let mut rng = rand::thread_rng();
    let mut c = Circuit::new();
    c.set_max_qubit_peak(1200);
    let a = c.alloc_qreg_bits("A", n);
    let b = c.alloc_qreg_bits("B", n);
    let q_div = c.alloc_qreg_bits("qd", n);
    let s = c.alloc_qreg_bits("s", 9);
    let s_rot = c.alloc_qreg_bits("srot", 9);
    let offset = c.alloc_qreg("off");
    let active = c.alloc_qreg("act");
    let mut avs = Vec::new();
    for shot in 0..64 {
        let mut mk = || -> BigUint {
            let bl = rng.gen_range(wlo..=whi);
            (BigUint::one() << (bl - 1))
                | (BigUint::from_bytes_le(&rng.gen::<[u8; 32]>()) % (BigUint::one() << (bl - 1)))
        };
        let (v1, v2) = (mk(), mk());
        let (av, bv) = if v1 >= v2 { (v1, v2) } else { (v2, v1) };
        load_big(&mut c, &a, &av, shot);
        load_big(&mut c, &b, &bv, shot);
        avs.push((av, bv));
    }
    c.x(&active);
    division_substep_windowed(
        &mut c, &a, &b, &q_div, &s_rot, &offset, &active, lo, lo, rot_bits,
    );
    division_substep_windowed_inv(
        &mut c, &a, &b, &q_div, &s_rot, &offset, &active, lo, lo, rot_bits,
    );
    // A,B restored, q_div back to 0.
    {
        let (ar, br, qr, vals) = (&a, &b, &q_div, avs.clone());
        c.contract_check("divrt", move |view, shot| {
            let (av, bv) = &vals[shot];
            let ga = rd_big(&view, ar, shot);
            let gb = rd_big(&view, br, shot);
            let gq = rd_big(&view, qr, shot);
            if &ga != av || &gb != bv || !gq.is_zero() {
                return Err(format!(
                    "roundtrip: A {ga:x}/{av:x} B {gb:x}/{bv:x} q {gq:x}"
                ));
            }
            Ok(())
        });
    }
    c.x(&active);
    for q in s.into_iter().chain(s_rot) {
        c.zero_and_free(q);
    }
    for q in [offset, active] {
        c.zero_and_free(q);
    }
    for q in q_div {
        c.zero_and_free(q);
    }
    c.assert_phase_clean();
    eprintln!("division_substep_windowed_roundtrip ok (forward.inverse = id)");
    let mut outs = vec![];
    outs.extend(a);
    outs.extend(b);
    let _ = c.destroy_sim(outs);
}

/// ROUND-TRIP: multiply_substep_windowed then its inverse = identity
/// (ca,cb,q_mul restored; s/s_rot/off/active clean).
#[test]
fn multiply_substep_windowed_roundtrip() {
    let n = 260usize;
    let lo = 247usize;
    let rot_bits = 3usize;
    let mut rng = rand::thread_rng();
    let mut c = Circuit::new();
    c.set_max_qubit_peak(1200);
    let ca = c.alloc_qreg_bits("ca", n);
    let cb = c.alloc_qreg_bits("cb", n);
    let q = c.alloc_qreg_bits("q", 16);
    let s = c.alloc_qreg_bits("s", 9);
    let s_rot = c.alloc_qreg_bits("srot", 9);
    let off = c.alloc_qreg("off");
    let active = c.alloc_qreg("act");
    let mut vals = Vec::new();
    for shot in 0..64 {
        let mut mk = || -> BigUint {
            let bl = rng.gen_range(248usize..=252);
            (BigUint::one() << (bl - 1))
                | (BigUint::from_bytes_le(&rng.gen::<[u8; 32]>()) % (BigUint::one() << (bl - 1)))
        };
        let (cav, cbv) = (mk(), mk());
        let s2 = rng.gen_range(1u32..=6);
        let qv = BigUint::one() << s2; // single bit -> ctz = s2
        load_big(&mut c, &ca, &cav, shot);
        load_big(&mut c, &cb, &cbv, shot);
        load_big(&mut c, &q, &qv, shot);
        vals.push((cav, cbv, qv));
    }
    c.x(&active);
    multiply_substep_windowed(
        &mut c, &ca, &cb, &q, &s_rot, &off, &active, lo, lo, rot_bits,
    );
    multiply_substep_windowed_inv(
        &mut c, &ca, &cb, &q, &s_rot, &off, &active, lo, lo, rot_bits,
    );
    {
        let (car, cbr, qr, vv) = (&ca, &cb, &q, vals.clone());
        c.contract_check("mulrt", move |view, shot| {
            let (cav, cbv, qv) = &vv[shot];
            let gca = rd_big(&view, car, shot);
            let gcb = rd_big(&view, cbr, shot);
            let gq = rd_big(&view, qr, shot);
            if &gca != cav || &gcb != cbv || &gq != qv {
                return Err(format!(
                    "mul roundtrip: ca {gca:x}/{cav:x} cb {gcb:x}/{cbv:x} q {gq:x}/{qv:x}"
                ));
            }
            Ok(())
        });
    }
    c.x(&active);
    for qreg in s.into_iter().chain(s_rot) {
        c.zero_and_free(qreg);
    }
    for qreg in [off, active] {
        c.zero_and_free(qreg);
    }
    c.assert_phase_clean();
    eprintln!("multiply_substep_windowed_roundtrip ok (forward.inverse = id)");
    let mut outs = vec![];
    outs.extend(ca);
    outs.extend(cb);
    outs.extend(q); // restored to its input -> keep, don't free
    let _ = c.destroy_sim(outs);
}

/// FULL FORWARD-BACK: run the windowed shrunken_pz inversion forward (x -> 1/x in cb),
/// copy 1/x out, run the gate-for-gate backward pass, and verify x is restored
/// and all state (counter, parity, ancillas) returns to its initial value --
/// i.e. the inversion is properly reversible. cb_out holds the correct 1/x.
#[test]
fn shrunken_pz_forward_back() {
    use crate::point_add::trailmix_port::inversion::shrunken_pz_schedule::{
        reg_los, reg_widths, shift_bounds, SHRUNKEN_PZ_NSTEPS,
    };
    fn resize(c: &mut Circuit, reg: &mut Vec<QReg>, target: usize, name: &'static str) {
        while reg.len() > target {
            let q = reg.pop().unwrap();
            c.zero_and_free(q);
        }
        while reg.len() < target {
            let k = reg.len();
            reg.push(c.alloc_qreg(&format!("{name}[{k}]")));
        }
    }
    // schedule-fit (widths + clz windows + shift bounds), same as the fwd test.
    fn shot_fits(x_orig: &BigUint, p: &BigUint, n: usize) -> bool {
        let one = BigUint::one();
        let half = p >> 1u32;
        let sgn = x_orig > &half;
        let mut a = p.clone();
        let mut b = if sgn { p - x_orig } else { x_orig.clone() };
        let (mut ca, mut cb) = (BigUint::zero(), one.clone());
        let mut q: u128 = 0;
        for i in 0..n {
            let (wa, wb, wca, wcb, wq) = reg_widths(i);
            let (alo, blo, calo, cblo, _) = reg_los(i);
            let (sdb, s2b) = shift_bounds(i);
            let (wg, wc) = (wa.max(wb), wca.max(wcb));
            if a.is_zero() && b == one && q == 0 {
                return true;
            }
            let mut hg = a.bits().max(b.bits()) as usize;
            let mut hc = ca.bits().max(cb.bits()) as usize;
            let mut hq = (128 - (q | 1).leading_zeros()) as usize;
            if a < b {
                let s2 = q.trailing_zeros() as usize;
                if s2 > s2b || (cb.bits() as usize) <= cblo {
                    return false;
                }
                let cbs = &cb << s2;
                hc = hc.max(cbs.bits() as usize);
                q ^= 1u128 << s2;
                ca += &cbs;
                hc = hc.max(ca.bits() as usize);
                hq = hq.max((128 - (q | 1).leading_zeros()) as usize);
                if (ca.bits() as usize) <= calo || (cbs.bits() as usize) <= calo {
                    return false;
                }
            }
            if ca < cb {
                if (a.bits() as usize) <= alo || (b.bits() as usize) <= blo {
                    return false;
                }
                let mut s: i64 = a.bits() as i64 - b.bits() as i64;
                if s >= 0 && (s as usize) > sdb {
                    return false;
                }
                let mut offset = false;
                if s >= 0 {
                    offset = a < (&b << (s as u32));
                }
                if offset {
                    s -= 1;
                }
                if s >= 0 {
                    let bsh = &b << (s as u32);
                    hg = hg.max(bsh.bits() as usize);
                    if a >= bsh {
                        a -= &bsh;
                        q ^= 1u128 << (s as u32);
                        hq = hq.max((128 - (q | 1).leading_zeros()) as usize);
                    }
                    hg = hg.max(a.bits() as usize);
                }
            }
            if q == 0 && !a.is_zero() {
                std::mem::swap(&mut a, &mut b);
                std::mem::swap(&mut ca, &mut cb);
            }
            if hg > wg || hc > wc || hq > wq {
                return false;
            }
        }
        false
    }

    let p = secp_p();
    let half = &p >> 1u32;
    let one = BigUint::one();
    let n = SHRUNKEN_PZ_NSTEPS;
    let mut rng = rand::thread_rng();
    let mut shots: Vec<BigUint> = Vec::new();
    while shots.len() < 64 {
        let x = (BigUint::from_bytes_le(&rng.gen::<[u8; 32]>()) % (&p - &one)) + &one;
        if shot_fits(&x, &p, n) {
            shots.push(x);
        }
    }

    let mut c = Circuit::new();
    c.set_max_qubit_peak(1400);
    let (a0, b0, ca0, cb0, q0) = reg_widths(0);
    let (wg0, wc0) = (a0.max(b0), ca0.max(cb0));
    let mut aa: Vec<QReg> = c.alloc_qreg_bits("A", wg0);
    let mut bb: Vec<QReg> = c.alloc_qreg_bits("B", wg0);
    let mut cca: Vec<QReg> = c.alloc_qreg_bits("ca", wc0);
    let mut ccb: Vec<QReg> = c.alloc_qreg_bits("cb", wc0);
    let mut qq: Vec<QReg> = c.alloc_qreg_bits("q", q0.max(1));
    let s = c.alloc_qreg_bits("s", 9);
    let s_rot = c.alloc_qreg_bits("srot", 9);
    let off = c.alloc_qreg("off");
    let parity = c.alloc_qreg("par");
    let counter = c.alloc_qreg_bits("done_ctr", 10);

    for (shot, x) in shots.iter().enumerate() {
        let xp = if x > &half { &p - x } else { x.clone() };
        load_big(&mut c, &aa, &p, shot);
        load_big(&mut c, &bb, &xp, shot);
        load_big(&mut c, &ccb, &one, shot);
    }
    c.x(&parity);

    // FORWARD pass (via the extracted driver).
    shrunken_pz_invert_forward(
        &mut c, &mut aa, &mut bb, &mut cca, &mut ccb, &mut qq, &counter, &parity, &s_rot, &off,
    );
    let fwd_peak = c.peak_qubits;
    let fwd_tof = c.executed_toffoli_shots / 64;

    // COPY the cofactor cb (the inverse, up to parity/sign) and the forward
    // parity out into persistent registers, so the backward can restore the
    // working state while we keep the result.
    let cb_out = c.alloc_qreg_bits("cb_out", ccb.len());
    for j in 0..ccb.len() {
        c.cx(&ccb[j], &cb_out[j]);
    }
    let par_out = c.alloc_qreg("par_out");
    c.cx(&parity, &par_out);

    // BACKWARD pass (via the extracted driver).
    shrunken_pz_invert_backward(
        &mut c, &mut aa, &mut bb, &mut cca, &mut ccb, &mut qq, &counter, &parity, &s_rot, &off,
    );
    let tot_tof = c.executed_toffoli_shots / 64;
    eprintln!(
        "shrunken_pz_forward_back: fwd peak {fwd_peak} q / {fwd_tof} tof, round-trip {tot_tof} tof"
    );

    // SHPZ_DEBUG=1 drops into the REPL with the full fwd+back circuit built,
    // for `prof whole tof top N` / `prof whole tof contains p.` exploration.
    if std::env::var("SHPZ_DEBUG").is_ok() {
        let mut d = crate::debugger::Debugger::attach(&mut c);
        d.repl();
    }

    // After backward: registers restored to initial (A=p, B=xp, ca=0, cb=1,
    // q=0, counter=0, parity=1). cb_out = 1/x. Verify, then free clean state.
    c.x(&parity); // parity should be back to 1 -> uncompute to 0 for free
    for q in s.into_iter().chain(s_rot).chain(counter) {
        c.zero_and_free(q);
    }
    for q in [off, parity] {
        c.zero_and_free(q);
    }
    // cca should be 0 -> free clean.
    for q in cca {
        c.zero_and_free(q);
    }
    for q in qq {
        c.zero_and_free(q);
    }

    let aa_len = aa.len();
    let bb_len = bb.len();
    let mut outs: Vec<QReg> = Vec::new();
    outs.extend(aa);
    let bb_off = outs.len();
    outs.extend(bb);
    let ccb_off = outs.len();
    outs.extend(ccb);
    let cbo_off = outs.len();
    let cbo_len = cb_out.len();
    outs.extend(cb_out);
    let par_idx = outs.len();
    outs.push(par_out);
    let (sim, det) = c.destroy_sim(outs);
    let rd = |det: &[QReg], off: usize, len: usize, shot: usize| -> BigUint {
        let mut v = BigUint::zero();
        for j in 0..len {
            if sim.read_bit_shot(&det[off + j], shot) == 1 {
                v |= BigUint::one() << j;
            }
        }
        v
    };
    for (shot, x) in shots.iter().enumerate() {
        let sgn = x > &half;
        let xp = if sgn { &p - x } else { x.clone() };
        let got_a = rd(&det, 0, aa_len, shot);
        let got_b = rd(&det, bb_off, bb_len, shot);
        let got_cb = rd(&det, ccb_off, cbo_off - ccb_off, shot);
        let got_cbo = rd(&det, cbo_off, cbo_len, shot);
        assert_eq!(got_a, p, "shot {shot}: A not restored to p");
        assert_eq!(got_b, xp, "shot {shot}: B not restored to xp");
        assert_eq!(got_cb, one, "shot {shot}: cb not restored to 1");
        // apply the shrunken_pz parity/sign adjustment to the raw cofactor (mirrors
        // shrunken_pz_crossgate_sched's inverse extraction).
        let par = sim.read_bit_shot(&det[par_idx], shot) == 1;
        let mut inv = if par {
            got_cbo.clone()
        } else {
            (&p - &got_cbo % &p) % &p
        };
        inv %= &p;
        if sgn {
            inv = (&p - inv) % &p;
        }
        let want = x.modpow(&(&p - BigUint::from(2u32)), &p);
        assert_eq!(inv, want, "shot {shot}: cb_out (adjusted) != x^-1");
    }
}

/// shrunken_pz_divide_forward: lambda = dy/dx mod p with dx & dy PRESERVED (dy via the
/// HMR-ghost trick). The EC-add slope primitive. Checks value, reversibility
/// (dx, dy restored), and phase-cleanliness over 64 random schedule-fit shots.
#[test]
fn shrunken_pz_divide_forward_test() {
    use crate::point_add::trailmix_port::inversion::shrunken_pz_schedule::SHRUNKEN_PZ_NSTEPS;

    let p = secp_p();
    let half = &p >> 1u32;
    let one = BigUint::one();
    let n = SHRUNKEN_PZ_NSTEPS;
    let mut rng = rand::thread_rng();
    let mut shots: Vec<(BigUint, BigUint)> = Vec::new();
    while shots.len() < 64 {
        let dx = (BigUint::from_bytes_le(&rng.gen::<[u8; 32]>()) % (&p - &one)) + &one;
        let dy = BigUint::from_bytes_le(&rng.gen::<[u8; 32]>()) % &p;
        let absdx = if dx > half { &p - &dx } else { dx.clone() };
        if shrunken_pz_shot_fits(&absdx, &p, n) {
            shots.push((dx, dy));
        }
    }

    let mut c = Circuit::new();
    c.set_max_qubit_peak(1300);
    let dx_reg: Vec<QReg> = c.alloc_qreg_bits("dx", 257);
    let dy_reg: Vec<QReg> = c.alloc_qreg_bits("dy", 257);
    for (shot, (dx, dy)) in shots.iter().enumerate() {
        load_big(&mut c, &dx_reg, dx, shot);
        load_big(&mut c, &dy_reg, dy, shot);
    }
    let (dx_out, dy_out, lambda) = shrunken_pz_divide_forward(&mut c, dx_reg, dy_reg);
    let peak = c.peak_qubits;
    eprintln!(
        "shrunken_pz_divide_forward_test: peak {peak}q, {} tof",
        c.executed_toffoli_shots / 64
    );
    c.assert_phase_clean();

    let (dxl, dyl, lal) = (dx_out.len(), dy_out.len(), lambda.len());
    let mut outs: Vec<QReg> = Vec::new();
    outs.extend(dx_out);
    let dyo = outs.len();
    outs.extend(dy_out);
    let lao = outs.len();
    outs.extend(lambda);
    let (sim, det) = c.destroy_sim(outs);
    let rd = |off: usize, len: usize, shot: usize| -> BigUint {
        let mut v = BigUint::zero();
        for j in 0..len {
            if sim.read_bit_shot(&det[off + j], shot) == 1 {
                v |= BigUint::one() << j;
            }
        }
        v
    };
    for (shot, (dx, dy)) in shots.iter().enumerate() {
        let got_dx = rd(0, dxl, shot);
        let got_dy = rd(dyo, dyl, shot);
        let got_lambda = rd(lao, lal, shot);
        assert_eq!(&got_dx, dx, "dx not preserved (shot {shot})");
        assert_eq!(&got_dy, dy, "dy not preserved (shot {shot})");
        let inv = dx.modpow(&(&p - 2u32), &p);
        let want = (dy * &inv) % &p;
        assert_eq!(got_lambda, want, "lambda != dy/dx (shot {shot})");
    }
}

/// DYNAMIC-W cross-gated shrunken_pz inversion: registers resized per step from the
/// shrunken_pz_schedule schedule (shrink gcd pair, grow cofactor pair), so every
/// barrel-rotate / comparator / clz runs at the live width, not flat 290.
/// Convergence is self-freezing: g_mul=(A<B AND A!=0), g_div=(ca<cb) both
/// evaluate to 0 once (A,B,q)=(0,1,0) and ca=P>cb. Full inversion; checks
/// cb == x^-1 mod p. Shots pre-filtered to fit the schedule (the ~0.1%
/// whole-pass tail is the accepted Shor failure rate).
#[test]
fn shrunken_pz_crossgate_sched() {
    use crate::point_add::trailmix_port::inversion::shrunken_pz_schedule::{
        reg_los, reg_widths, shift_bounds, SHRUNKEN_PZ_NSTEPS,
    };
    // rotator shift-bits to cover a max shift `b`: ceil-log2 so 2^bits > b.
    fn rot_bits_for(b: usize) -> usize {
        if b == 0 {
            1
        } else {
            64 - (b as u64).leading_zeros() as usize
        }
    }

    fn resize(c: &mut Circuit, reg: &mut Vec<QReg>, target: usize, name: &'static str) {
        while reg.len() > target {
            let q = reg.pop().unwrap();
            c.zero_and_free(q);
        }
        while reg.len() < target {
            // tag grown bits `name[k]` (consistent with alloc_qreg_bits) so
            // the dump hook / read_tagged_reg follows the full register.
            let k = reg.len();
            reg.push(c.alloc_qreg(&format!("{name}[{k}]")));
        }
    }

    // Classical fit check: run shrunken_pz and verify every register (incl B<<s /
    // cb<<s2 transients) stays within the PAIRED scheduled width at every
    // step, and converges within n_steps.
    fn shot_fits(x_orig: &BigUint, p: &BigUint, n: usize) -> bool {
        let one = BigUint::one();
        let half = p >> 1u32;
        let sgn = x_orig > &half;
        let mut a = p.clone();
        let mut b = if sgn { p - x_orig } else { x_orig.clone() };
        let (mut ca, mut cb) = (BigUint::zero(), one.clone());
        let mut q: u128 = 0;
        for i in 0..n {
            let (wa, wb, wca, wcb, wq) = reg_widths(i);
            let (alo, blo, calo, cblo, _) = reg_los(i);
            let (sdb, s2b) = shift_bounds(i);
            let wg = wa.max(wb);
            let wc = wca.max(wcb);
            if a.is_zero() && b == one && q == 0 {
                return true;
            }
            let mut hg = a.bits().max(b.bits()) as usize; // max gcd-pair bitlen this step
            let mut hc = ca.bits().max(cb.bits()) as usize;
            let mut hq = (128 - (q | 1).leading_zeros()) as usize;
            if a < b {
                // MULTIPLY active: windowed clz on (ca,cb<<s2) [ca window] and
                // (cb,ca) [cb/ca windows]; rotate distance s2 <= s2b.
                let s2 = q.trailing_zeros() as usize;
                if s2 > s2b {
                    return false;
                }
                let cbs = &cb << s2;
                if (cb.bits() as usize) <= cblo {
                    return false;
                }
                hc = hc.max(cbs.bits() as usize); // cb<<s2 transient
                q ^= 1u128 << s2;
                ca += &cbs;
                hc = hc.max(ca.bits() as usize);
                hq = hq.max((128 - (q | 1).leading_zeros()) as usize);
                if (ca.bits() as usize) <= calo || (cbs.bits() as usize) <= calo {
                    return false; // grown ca / cb<<s2 must be in the ca window
                }
            }
            if ca < cb {
                // DIVISION active: windowed clz must see A,B MSB in the window,
                // and the rotate distance (s_initial = bl(A)-bl(B)) must fit.
                if (a.bits() as usize) <= alo || (b.bits() as usize) <= blo {
                    return false;
                }
                let mut s: i64 = a.bits() as i64 - b.bits() as i64;
                if s >= 0 && (s as usize) > sdb {
                    return false; // rotate distance exceeds rot_bits bound
                }
                let mut offset = false;
                if s >= 0 {
                    offset = a < (&b << (s as u32));
                }
                if offset {
                    s -= 1;
                }
                if s >= 0 {
                    let bsh = &b << (s as u32);
                    hg = hg.max(bsh.bits() as usize); // B<<s transient
                    if a >= bsh {
                        a -= &bsh;
                        q ^= 1u128 << (s as u32);
                        hq = hq.max((128 - (q | 1).leading_zeros()) as usize);
                    }
                    hg = hg.max(a.bits() as usize);
                }
            }
            if q == 0 && !a.is_zero() {
                std::mem::swap(&mut a, &mut b);
                std::mem::swap(&mut ca, &mut cb);
            }
            if hg > wg || hc > wc || hq > wq {
                return false;
            }
        }
        false
    }

    // shrunken_pz model trajectory: per step, the POST-MULTIPLY state (after the
    // multiply substep, before division) and the FINAL state (after
    // division+swap). Lets us localize which substep first diverges.
    type St = (BigUint, BigUint, BigUint, BigUint, u128);
    fn shrunken_pz_traj(x_orig: &BigUint, p: &BigUint, n: usize) -> Vec<(St, St)> {
        let one = BigUint::one();
        let half = p >> 1u32;
        let sgn = x_orig > &half;
        let mut a = p.clone();
        let mut b = if sgn { p - x_orig } else { x_orig.clone() };
        let (mut ca, mut cb) = (BigUint::zero(), one.clone());
        let mut q: u128 = 0;
        let mut out = Vec::with_capacity(n);
        for _ in 0..n {
            let converged = a.is_zero() && b == one && q == 0;
            if !converged && a < b {
                let s2 = q.trailing_zeros();
                q ^= 1u128 << s2;
                ca += &cb << s2;
            }
            let post_mul = (a.clone(), b.clone(), ca.clone(), cb.clone(), q);
            if !converged && ca < cb {
                let mut s: i64 = a.bits() as i64 - b.bits() as i64;
                let mut offset = false;
                if s >= 0 {
                    offset = a < (&b << (s as u32));
                }
                if offset {
                    s -= 1;
                }
                if s >= 0 {
                    let bsh = &b << (s as u32);
                    if a >= bsh {
                        a -= &bsh;
                        q ^= 1u128 << (s as u32);
                    }
                }
            }
            if !converged && q == 0 && !a.is_zero() {
                std::mem::swap(&mut a, &mut b);
                std::mem::swap(&mut ca, &mut cb);
            }
            let fin = (a.clone(), b.clone(), ca.clone(), cb.clone(), q);
            out.push((post_mul, fin));
        }
        out
    }

    let p = secp_p();
    let half = &p >> 1u32;
    let one = BigUint::one();
    let n = SHRUNKEN_PZ_NSTEPS;
    let mut rng = rand::thread_rng();

    // collect 64 schedule-fitting shots + their model trajectories
    let mut shots: Vec<BigUint> = Vec::new();
    let mut trajs: Vec<Vec<(St, St)>> = Vec::new();
    while shots.len() < 64 {
        let x = (BigUint::from_bytes_le(&rng.gen::<[u8; 32]>()) % (&p - &one)) + &one;
        if shot_fits(&x, &p, n) {
            trajs.push(shrunken_pz_traj(&x, &p, n));
            shots.push(x);
        }
    }

    let mut c = Circuit::new();
    c.set_max_qubit_peak(900);
    let (a0, b0, ca0, cb0, q0) = reg_widths(0);
    let wg0 = a0.max(b0);
    let wc0 = ca0.max(cb0);
    let mut aa: Vec<QReg> = c.alloc_qreg_bits("A", wg0);
    let mut bb: Vec<QReg> = c.alloc_qreg_bits("B", wg0);
    let mut cca: Vec<QReg> = c.alloc_qreg_bits("ca", wc0);
    let mut ccb: Vec<QReg> = c.alloc_qreg_bits("cb", wc0);
    let mut qq: Vec<QReg> = c.alloc_qreg_bits("q", q0.max(1));
    let s = c.alloc_qreg_bits("s", 9);
    let s_rot = c.alloc_qreg_bits("srot", 9);
    let off = c.alloc_qreg("off");
    let parity = c.alloc_qreg("par");
    // DONE COUNTER: counts steps spent in the converged fixed point. active
    // = (counter==0) gates the substeps; once converged, counter starts and
    // freezes everything. The per-step `done` flag cleans via (counter!=0)
    // (user's recipe) so there is no leak. ~10 bits covers the over-run tail.
    let counter = c.alloc_qreg_bits("done_ctr", 10);

    // Debugger dump hook: `dump shrunken_pz <shot>` prints the algebraic state
    // (A,B,ca,cb,q,counter as integers + the EEA residual ca*x,cb*x mod p)
    // at the cursor, following the dynamic-W registers by tag.
    {
        let xs = shots.clone();
        let p_h = p.clone();
        let half_h = half.clone();
        c.register_debug_dump("shrunken_pz", move |dbg, args| {
            let shot: usize = args.first().and_then(|s| s.parse().ok()).unwrap_or(0);
            if shot >= xs.len() {
                return format!("shot {shot} out of range");
            }
            let a = dbg.read_tagged_reg("A", shot, 300);
            let b = dbg.read_tagged_reg("B", shot, 300);
            let ca = dbg.read_tagged_reg("ca", shot, 300);
            let cb = dbg.read_tagged_reg("cb", shot, 300);
            let q = dbg.read_tagged_reg("q", shot, 64);
            let ctr = dbg.read_tagged_reg("done_ctr", shot, 16);
            let x = &xs[shot];
            let xi = if x > &half_h { &p_h - x } else { x.clone() };
            let cax = (&ca * &xi) % &p_h;
            let cbx = (&cb * &xi) % &p_h;
            format!(
                "[shrunken_pz] shot {shot} cursor={} sec={}\n  A ={a:x}\n  B ={b:x}\n  \
                 ca={ca:x}\n  cb={cb:x}\n  q ={q:x}\n  done_ctr={ctr}\n  \
                 ca*x%p={cax:x}\n  cb*x%p={cbx:x}  (x_int={xi:x})",
                dbg.cursor(),
                dbg.current_section()
            )
        });
    }

    for (shot, x) in shots.iter().enumerate() {
        let xp = if x > &half { &p - x } else { x.clone() };
        load_big(&mut c, &aa, &p, shot);
        load_big(&mut c, &bb, &xp, shot);
        load_big(&mut c, &ccb, &one, shot);
    }
    c.x(&parity);

    for i in 0..n {
        let (wa, wb, wca, wcb, wq) = reg_widths(i);
        let wg = wa.max(wb);
        let wc = wca.max(wcb);
        resize(&mut c, &mut aa, wg, "A");
        resize(&mut c, &mut bb, wg, "B");
        resize(&mut c, &mut cca, wc, "ca");
        resize(&mut c, &mut ccb, wc, "cb");
        resize(&mut c, &mut qq, wq.max(1), "q");

        let ps_mul = c.push_section("shrunken_pz.mul");
        // active = (counter == 0): gates the substeps. Once converged the
        // counter starts and active drops, freezing everything. counter is
        // not touched by the substeps, so active recompute-cleans.
        let active = c.alloc_qreg("active");
        or_is_zero(&mut c, &counter, &active);

        // MULTIPLY, gated g_mul=(A<B AND active). Leaves A,B fixed -> clean.
        // (No A!=0 -- the final q-drain steps have A=0 but must still run;
        // convergence is handled by the done-counter, not this gate.)
        let g_mul = c.alloc_qreg("g_mul");
        {
            let ar: Vec<&QReg> = aa.iter().collect();
            let br: Vec<&QReg> = bb.iter().collect();
            let lt = c.alloc_qreg("gm.lt");
            borrow_compare_refs(&mut c, &ar, &br, &lt); // lt=(A<B)
            c.ccx(&lt, &active, &g_mul); // g_mul = (A<B) AND active
            borrow_compare_refs(&mut c, &ar, &br, &lt);
            c.zero_and_free(lt);
        }
        {
            let (_, _, calo, cblo, _) = reg_los(i);
            let (_, s2b) = shift_bounds(i);
            multiply_substep_windowed(
                &mut c,
                &cca,
                &ccb,
                &qq,
                &s_rot,
                &off,
                &g_mul,
                calo,
                cblo,
                rot_bits_for(s2b),
            );
        }
        {
            let ar: Vec<&QReg> = aa.iter().collect();
            let br: Vec<&QReg> = bb.iter().collect();
            let lt = c.alloc_qreg("gm.lt2");
            borrow_compare_refs(&mut c, &ar, &br, &lt);
            c.ccx(&lt, &active, &g_mul);
            borrow_compare_refs(&mut c, &ar, &br, &lt);
            c.zero_and_free(lt);
        }
        c.zero_and_free(g_mul);
        c.pop_section(&ps_mul);

        // POST-MULTIPLY invariant: catches if the MULTIPLY substep alone
        // diverges (before division runs).
        {
            let (ar, br, car, cbr, qr) = (&aa, &bb, &cca, &ccb, &qq);
            let tr = &trajs;
            let (wgi, wci, wqi) = (aa.len(), cca.len(), qq.len());
            c.contract_check("shrunken_pz_postmul", |view, shot| {
                let (ea, eb, eca, ecb, eq) = &tr[shot][i].0;
                let ga = rd_big(&view, ar, shot);
                let gb = rd_big(&view, br, shot);
                let gca = rd_big(&view, car, shot);
                let gcb = rd_big(&view, cbr, shot);
                let gq = rd_big(&view, qr, shot);
                let eqb = BigUint::from(*eq);
                if &ga == ea && &gb == eb && &gca == eca && &gcb == ecb && gq == eqb {
                    return Ok(());
                }
                let d = |tag: &str, c: &BigUint, m: &BigUint| {
                    let mk = if c == m { "  " } else { "!!" };
                    format!("\n  {mk} {tag:3} c={c:x}\n        m={m:x}")
                };
                Err(format!(
                    "POSTMUL step {i} shot {shot} (gcd={wgi} cof={wci} q={wqi}){}{}{}{}{}",
                    d("A", &ga, ea),
                    d("B", &gb, eb),
                    d("ca", &gca, eca),
                    d("cb", &gcb, ecb),
                    d("q", &gq, &eqb),
                ))
            });
        }

        let ps_div = c.push_section("shrunken_pz.div");
        // DIVISION, gated g_div=(ca<cb AND active). Leaves ca,cb fixed.
        let g_div = c.alloc_qreg("g_div");
        {
            let car: Vec<&QReg> = cca.iter().collect();
            let cbr: Vec<&QReg> = ccb.iter().collect();
            let lt = c.alloc_qreg("gd.lt");
            borrow_compare_refs(&mut c, &car, &cbr, &lt); // lt=(ca<cb)
            c.ccx(&lt, &active, &g_div);
            borrow_compare_refs(&mut c, &car, &cbr, &lt);
            c.zero_and_free(lt);
        }
        {
            let (alo, blo, _, _, _) = reg_los(i);
            let (sdb, _) = shift_bounds(i);
            division_substep_windowed(
                &mut c,
                &aa,
                &bb,
                &qq,
                &s_rot,
                &off,
                &g_div,
                alo,
                blo,
                rot_bits_for(sdb),
            );
        }
        {
            let car: Vec<&QReg> = cca.iter().collect();
            let cbr: Vec<&QReg> = ccb.iter().collect();
            let lt = c.alloc_qreg("gd.lt2");
            borrow_compare_refs(&mut c, &car, &cbr, &lt);
            c.ccx(&lt, &active, &g_div);
            borrow_compare_refs(&mut c, &car, &cbr, &lt);
            c.zero_and_free(lt);
        }
        c.zero_and_free(g_div);
        c.pop_section(&ps_div);

        let ps_swap = c.push_section("shrunken_pz.swap");
        // SWAP, gated g_swap=(q==0 & A!=0 & active).
        let g_swap = c.alloc_qreg("g_swap");
        let mk_gswap = |c: &mut Circuit, qq: &[QReg], aa: &[QReg], active: &QReg, g: &QReg| {
            let qz = c.alloc_qreg("sw.qz");
            let anz = c.alloc_qreg("sw.anz");
            let t = c.alloc_qreg("sw.t");
            or_is_zero(c, qq, &qz);
            or_nonzero(c, aa, &anz);
            c.ccx(&qz, &anz, &t); // t = (q==0 & A!=0)
            c.ccx(&t, active, g); // g_swap = t AND active
            c.ccx(&qz, &anz, &t);
            or_nonzero(c, aa, &anz);
            or_is_zero(c, qq, &qz);
            c.zero_and_free(t);
            c.zero_and_free(anz);
            c.zero_and_free(qz);
        };
        mk_gswap(&mut c, &qq, &aa, &active, &g_swap);
        for j in 0..wg {
            c.cswap(&g_swap, &aa[j], &bb[j]);
        }
        for j in 0..wc {
            c.cswap(&g_swap, &cca[j], &ccb[j]);
        }
        c.cx(&g_swap, &parity);
        mk_gswap(&mut c, &qq, &aa, &active, &g_swap); // uncompute (preserved by swap)
        c.zero_and_free(g_swap);
        c.pop_section(&ps_swap);

        let ps_ctr = c.push_section("shrunken_pz.ctr");
        // uncompute active (counter untouched by substeps -> clean).
        or_is_zero(&mut c, &counter, &active);
        c.zero_and_free(active);

        // DONE-COUNTER update (user's recipe): done = (A==0 & q==0); if done
        // counter += 1; clean done via done ^= (counter != 0). No leak.
        // (A==0 & q==0) <=> converged: A=0 only after the gcd completes, and
        // q=0 with A=0 is the fixed point -- q=0 with A!=0 triggers a swap.
        {
            let done = c.alloc_qreg("done");
            let az = c.alloc_qreg("d.az");
            let qz = c.alloc_qreg("d.qz");
            or_is_zero(&mut c, &aa, &az);
            or_is_zero(&mut c, &qq, &qz);
            c.ccx(&az, &qz, &done); // done = (A==0 & q==0)
            or_is_zero(&mut c, &qq, &qz);
            or_is_zero(&mut c, &aa, &az);
            c.zero_and_free(qz);
            c.zero_and_free(az);
            ctrl_inc(&mut c, &done, &counter); // counter += done
            let cnz = c.alloc_qreg("d.cnz");
            or_nonzero(&mut c, &counter, &cnz);
            c.cx(&cnz, &done); // done ^= (counter != 0) -> clears done
            or_nonzero(&mut c, &counter, &cnz);
            c.zero_and_free(cnz);
            c.zero_and_free(done);
        }
        c.pop_section(&ps_ctr);

        // per-step INVARIANT: circuit (A,B,ca,cb,q) must match the shrunken_pz model.
        // Fires at the FIRST divergent step (DEBUG_ON_FAIL attaches there).
        {
            let (ar, br, car, cbr, qr) = (&aa, &bb, &cca, &ccb, &qq);
            let tr = &trajs;
            let (wgi, wci, wqi) = (aa.len(), cca.len(), qq.len());
            c.contract_check("shrunken_pz_inv", |view, shot| {
                let (ea, eb, eca, ecb, eq) = &tr[shot][i].1;
                let ga = rd_big(&view, ar, shot);
                let gb = rd_big(&view, br, shot);
                let gca = rd_big(&view, car, shot);
                let gcb = rd_big(&view, cbr, shot);
                let gq = rd_big(&view, qr, shot);
                let eqb = BigUint::from(*eq);
                let ok = &ga == ea && &gb == eb && &gca == eca && &gcb == ecb && gq == eqb;
                if !ok {
                    // ALGEBRAIC STATE DUMP: circuit (c=) vs model (m=).
                    let d = |tag: &str, c: &BigUint, m: &BigUint| {
                        let mark = if c == m { "  " } else { "!!" };
                        format!("\n  {mark} {tag:3} c={c:x}\n        m={m:x}")
                    };
                    return Err(format!(
                        "step {i} shot {shot} (widths gcd={wgi} cof={wci} q={wqi}){}{}{}{}{}",
                        d("A", &ga, ea),
                        d("B", &gb, eb),
                        d("ca", &gca, eca),
                        d("cb", &gcb, ecb),
                        d("q", &gq, &eqb),
                    ));
                }
                Ok(())
            });
        }
    }

    let peak = c.peak_qubits;
    let tof = c.executed_toffoli_shots / 64;
    eprintln!("shrunken_pz_crossgate_sched: {n} steps, peak {peak} q, tof {tof}");

    for q in s {
        c.zero_and_free(q);
    }
    for q in s_rot {
        c.zero_and_free(q);
    }
    c.zero_and_free(off);

    // SHPZ_DEBUG=1 drops into the REPL here (forward pass built, scratch freed)
    // for interactive `prof whole tof top N` / `dump shrunken_pz <shot>` exploration.
    if std::env::var("SHPZ_DEBUG").is_ok() {
        let mut d = crate::debugger::Debugger::attach(&mut c);
        d.repl();
    }

    // read final cb + parity, check inverse. Move registers into destroy_sim.
    let cb_len = ccb.len();
    let cb_off = 0usize;
    let mut outs: Vec<QReg> = Vec::new();
    outs.extend(ccb); // cb at [0..cb_len)
    let par_idx = outs.len();
    outs.push(parity);
    outs.extend(aa);
    outs.extend(bb);
    outs.extend(cca);
    outs.extend(qq);
    outs.extend(counter); // persistent over-run count (not garbage)
    let (sim, detached) = c.destroy_sim(outs);

    for (shot, x) in shots.iter().enumerate() {
        let sgn = x > &half;
        let mut cb_val = BigUint::zero();
        for j in 0..cb_len {
            if sim.read_bit_shot(&detached[cb_off + j], shot) == 1 {
                cb_val |= BigUint::one() << j;
            }
        }
        let par = sim.read_bit_shot(&detached[par_idx], shot) == 1;
        let mut inv = if par {
            cb_val.clone()
        } else {
            &p - (&cb_val % &p)
        };
        inv %= &p;
        if sgn {
            inv = (&p - inv) % &p;
        }
        let want = x.modpow(&(&p - BigUint::from(2u32)), &p);
        assert_eq!(inv, want, "shot {shot}: cb != x^-1 (x={x:x})");
    }
}

#[test]
fn pz_big_step_ref_correct() {
    let p = secp_p();
    let mut rng = rand::thread_rng();
    let mut next = || {
        let mut bytes = [0u8; 32];
        rng.fill(&mut bytes);
        BigUint::from_bytes_le(&bytes) % &p
    };
    let mut max_iters = 0usize;
    let mut max_peak: i64 = 0;
    let mut iter_hist = vec![0usize; 0];
    let mut s_hist = [0u64; 16];
    let mut wmax = Box::new([[0u16; 512]; 4]);
    for _ in 0..3000 {
        let x = {
            let v = next();
            if v.is_zero() {
                BigUint::one()
            } else {
                v
            }
        };
        let (inv, iters, peak) = pz_big_step_ref(&x, &mut s_hist, &mut wmax);
        assert_eq!((&inv * &x) % &p, BigUint::one(), "inverse wrong for x={x}");
        max_iters = max_iters.max(iters);
        max_peak = max_peak.max(peak);
        iter_hist.push(iters);
    }
    iter_hist.sort_unstable();
    let p99 = iter_hist[iter_hist.len() * 99 / 100];
    let tot: u64 = s_hist.iter().sum();
    let le2: u64 = s_hist[0] + s_hist[1] + s_hist[2];
    eprintln!(
        "pz_big_step_ref OK on 3000 secp inputs: max_iters={max_iters} p99={p99} \
         per-shot/cursor peak(bl A+B+a+b)={max_peak}"
    );
    eprintln!(
        "division rotation s: s<=2 = {:.4}% ; hist[0..8]={:?} (>=8 in last bins)",
        100.0 * le2 as f64 / tot as f64,
        &s_hist[0..8]
    );
    let unpacked_peak = (0..512)
        .map(|i| wmax[0][i] as u32 + wmax[1][i] as u32 + wmax[2][i] as u32 + wmax[3][i] as u32)
        .max()
        .unwrap();
    eprintln!(
        "UNPACKED register peak (sum of per-iter per-reg max over 3000 shots) = {unpacked_peak} \
         [vs per-shot/cursor peak {max_peak}]"
    );
}
