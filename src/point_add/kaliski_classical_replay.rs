//! Classical replay of the current Kaliski iteration to test whether `m_i`
//! is recoverable from a small live-state window `(u[0], v_w[0], s[0],
//! iter_idx, f)` at iteration start.
//!
//! Purpose: feasibility check for §7/§9 of `kaliski_1200q_feasibility.md`.
//! If `m_i` can be recomputed from a constant-size fingerprint of the live
//! Kaliski state at each iteration, then `m_hist` (407 qubits) is not
//! persistent state and Kaliski fits in 1200q (together with the other
//! compressions analyzed).
//!
//! This is analysis-only; does not change the quantum circuit. Run via:
//!     KALISKI_REPLAY=1 cargo test --release classical_replay
//! or via the top-level harness by calling
//! `kaliski_classical_replay::run_feasibility_test()`.

#![allow(dead_code)]

use alloy_primitives::U256;
use sha3::{
    digest::{ExtendableOutput, Update, XofReader},
    Shake128,
};

pub const SECP256K1_P: U256 = U256::from_limbs([
    0xFFFFFFFEFFFFFC2F,
    0xFFFFFFFFFFFFFFFF,
    0xFFFFFFFFFFFFFFFF,
    0xFFFFFFFFFFFFFFFF,
]);

/// One Kaliski iteration, classical replay of the quantum circuit in
/// `mod.rs::kaliski_iteration`. Returns `m_i` for this iteration and
/// mutates the state in place.
///
/// Arguments:
///   (u, v_w, r, s): current state (as U256, full width)
///   f: termination flag (1 = not terminated yet; 0 = terminated)
///   p: modulus
///
/// Returns: m_i for this iteration (0 or 1).
pub fn kaliski_iter_classical(
    u: &mut U256,
    v_w: &mut U256,
    r: &mut U256,
    s: &mut U256,
    f: &mut u8,
    p: U256,
) -> u8 {
    // STEP 0: is_zero = (v_w == 0); m_i ^= (f AND is_zero); f ^= m_i.
    let is_zero = if *v_w == U256::ZERO { 1u8 } else { 0 };
    let mut m_i: u8 = 0;
    if *f == 1 && is_zero == 1 {
        m_i ^= 1;
    }
    *f ^= m_i;

    // If f is now 0, the rest of this iteration is a no-op in terms of state
    // (a_f, b_f, add_f stay 0). But we still run the mechanical replay so
    // ordering matches the quantum circuit.

    // STEP 1:
    //   a_f ^= (f AND NOT u[0])
    //   m_i ^= (f AND u[0] AND NOT v_w[0])
    //   b_f = a_f XOR m_i
    let u0 = (u.as_limbs()[0] & 1) as u8;
    let v0 = (v_w.as_limbs()[0] & 1) as u8;
    let mut a_f: u8 = 0;
    if *f == 1 && u0 == 0 {
        a_f ^= 1;
    }
    if *f == 1 && u0 == 1 && v0 == 0 {
        m_i ^= 1;
    }
    let mut b_f = a_f ^ m_i;

    // STEP 2:
    //   l_gt = (u > v_w)
    //   add_f = f AND l_gt
    //   a_f ^= (add_f AND NOT b_f_orig)
    //   m_i ^= (add_f AND NOT b_f_orig)
    // Note: the circuit toggles b_f's polarity around the ccx using x(b_f)
    // on both sides, so the effective control is "NOT b_f" of the b_f value
    // entering STEP 2.
    let l_gt = if *u > *v_w { 1u8 } else { 0 };
    let add_f_step2 = (*f & l_gt) as u8;
    let b_not = 1 ^ b_f;
    let delta = add_f_step2 & b_not;
    a_f ^= delta;
    m_i ^= delta;
    // b_f is updated to new a_f ^ m_i implicitly (cx chain at end of step 1
    // sets b_f = a_f XOR m_i; STEP 2 further modifies a_f and m_i but not
    // b_f directly).

    // STEP 3: with control a_f: swap(u, v_w); swap(r, s).
    if a_f == 1 {
        std::mem::swap(u, v_w);
        std::mem::swap(r, s);
    }

    // STEP 4: add_f = f AND NOT b_f; with control add_f: v_w -= u; s += r.
    let add_f_step4 = *f & (1 ^ b_f);
    if add_f_step4 == 1 {
        // v_w -= u (mod 2^n for the register width, but in reality kaliski
        // keeps v_w < u after the swap, so this is just v_w - u as a
        // nonneg integer if swap happened appropriately; in the quantum
        // circuit this is mod 2^n).
        *v_w = v_w.wrapping_sub(*u);
        // s += r mod p (actually in the quantum circuit this is mod 2^n
        // but the Kaliski invariant keeps s < p always — check).
        let sum = s.wrapping_add(*r);
        *s = if sum < *s || sum >= p {
            // carry or overflow beyond p? Kaliski invariant says s+r ≤ p,
            // so we should not need mod reduction. Preserve wrapping for
            // faithful replay.
            sum
        } else {
            sum
        };
    }

    // STEP 5: uncompute add_f, b_f. State unchanged for u,v_w,r,s.

    // STEP 6: v_w := v_w >> 1 (unconditional, always safe).
    *v_w = *v_w >> 1;

    // STEP 7+8: r := 2r mod p.
    let r2 = r.wrapping_add(*r);
    *r = if r2 >= p || r2 < *r { r2.wrapping_sub(p) } else { r2 };

    // STEP 9: with control a_f: swap(u, v_w); swap(r, s).
    if a_f == 1 {
        std::mem::swap(u, v_w);
        std::mem::swap(r, s);
    }

    // STEP 10: a_f uncomputed from s[0] invariant (state unchanged).

    m_i
}

/// Classical driver: run `iters` iterations of Kaliski on input `v_in`,
/// starting from the same initial state the circuit does: u=p, v_w=v_in,
/// r=0, s=1, f=1.
///
/// Returns a vector of m_i values (one per iter) and final state.
pub fn kaliski_run(
    v_in: U256,
    p: U256,
    iters: usize,
) -> (Vec<u8>, Vec<(U256, U256, U256, U256, u8)>) {
    let mut u = p;
    let mut v_w = v_in;
    let mut r = U256::ZERO;
    let mut s = U256::from(1u64);
    let mut f: u8 = 1;

    let mut m_hist = Vec::with_capacity(iters);
    let mut snapshots = Vec::with_capacity(iters);

    for _i in 0..iters {
        snapshots.push((u, v_w, r, s, f));
        let m = kaliski_iter_classical(&mut u, &mut v_w, &mut r, &mut s, &mut f, p);
        m_hist.push(m);
    }

    (m_hist, snapshots)
}

/// Generate a pseudorandom but reproducible secp256k1 element.
fn random_element(seed: u64) -> U256 {
    let mut h = Shake128::default();
    h.update(&seed.to_le_bytes());
    let mut reader = h.finalize_xof();
    loop {
        let mut buf = [0u8; 32];
        reader.read(&mut buf);
        let v = U256::from_be_bytes(buf);
        if v != U256::ZERO && v < SECP256K1_P {
            return v;
        }
    }
}

/// The central feasibility test: can we recompute `m_i` at iteration `i`
/// from a constant-size fingerprint of the live state?
///
/// We test several candidate fingerprints:
///   F1 = (f, u[0], v_w[0])              — 3 bits
///   F2 = (f, u[0], v_w[0], s[0])        — 4 bits
///   F3 = (f, u[0], v_w[0], r[0], s[0])  — 5 bits
///   F4 = (f, u[0], v_w[0], u>v_w, s[0]) — 5 bits (includes comparator)
///   F5 = F4 + all low bits              — includes "u % 4" etc.
///
/// For each fingerprint, check whether `m_i` is a deterministic function
/// of the fingerprint across all (input, iter) pairs seen.
/// If ANY fingerprint is a deterministic function, m_hist elimination via
/// HRSL-style recomputation is possible.
pub fn run_feasibility_test(n_inputs: usize, iters_per_input: usize) -> bool {
    let p = SECP256K1_P;

    // Five candidate fingerprint tables: fingerprint -> set of m_i values
    // observed. If any fingerprint never maps to both 0 and 1, it's
    // deterministic.
    use std::collections::HashMap;
    let mut f1: HashMap<u8, u8> = HashMap::new();
    let mut f2: HashMap<u8, u8> = HashMap::new();
    let mut f3: HashMap<u8, u8> = HashMap::new();
    let mut f4: HashMap<u8, u8> = HashMap::new();
    let mut fmin: HashMap<u8, u8> = HashMap::new();
    let mut f1_conflicts = 0usize;
    let mut f2_conflicts = 0usize;
    let mut f3_conflicts = 0usize;
    let mut f4_conflicts = 0usize;
    let mut fmin_conflicts = 0usize;

    let mut f1_conflict_examples: Vec<(u8, u8)> = Vec::new();
    let mut f4_conflict_examples: Vec<(u8, u8)> = Vec::new();

    let mut total_samples = 0usize;

    for seed in 0..n_inputs {
        let v_in = random_element(seed as u64 + 1);
        let (m_hist, snaps) = kaliski_run(v_in, p, iters_per_input);
        for i in 0..iters_per_input {
            let (u, v_w, r, s, f) = snaps[i];
            let m_i = m_hist[i];
            let u0 = (u.as_limbs()[0] & 1) as u8;
            let v0 = (v_w.as_limbs()[0] & 1) as u8;
            let r0 = (r.as_limbs()[0] & 1) as u8;
            let s0 = (s.as_limbs()[0] & 1) as u8;
            let gt = if u > v_w { 1u8 } else { 0 };

            // F1: (f, u0, v0)
            let k1 = (f << 2) | (u0 << 1) | v0;
            match f1.get(&k1) {
                None => {
                    f1.insert(k1, m_i);
                }
                Some(&v) if v != m_i => {
                    f1_conflicts += 1;
                    if f1_conflict_examples.len() < 4 {
                        f1_conflict_examples.push((k1, m_i));
                    }
                }
                _ => {}
            }

            // F2: (f, u0, v0, s0)
            let k2 = (f << 3) | (u0 << 2) | (v0 << 1) | s0;
            match f2.get(&k2) {
                None => {
                    f2.insert(k2, m_i);
                }
                Some(&v) if v != m_i => {
                    f2_conflicts += 1;
                }
                _ => {}
            }

            // F3: (f, u0, v0, r0, s0)
            let k3 = (f << 4) | (u0 << 3) | (v0 << 2) | (r0 << 1) | s0;
            match f3.get(&k3) {
                None => {
                    f3.insert(k3, m_i);
                }
                Some(&v) if v != m_i => {
                    f3_conflicts += 1;
                }
                _ => {}
            }

            // F4: (f, u0, v0, gt, s0)   <-- includes the u>v_w comparator
            let k4 = (f << 4) | (u0 << 3) | (v0 << 2) | (gt << 1) | s0;
            match f4.get(&k4) {
                None => {
                    f4.insert(k4, m_i);
                }
                Some(&v) if v != m_i => {
                    f4_conflicts += 1;
                    if f4_conflict_examples.len() < 4 {
                        f4_conflict_examples.push((k4, m_i));
                    }
                }
                _ => {}
            }

            // F_min: (f, u0, v0, gt) -- MINIMAL candidate without s0.
            let kmin = (f << 3) | (u0 << 2) | (v0 << 1) | gt;
            match fmin.get(&kmin) {
                None => {
                    fmin.insert(kmin, m_i);
                }
                Some(&v) if v != m_i => {
                    fmin_conflicts += 1;
                }
                _ => {}
            }

            total_samples += 1;
        }
    }

    println!(
        "=== Kaliski m_i-from-fingerprint feasibility test ===
inputs={} iters/input={} total_samples={}",
        n_inputs, iters_per_input, total_samples
    );
    println!(
        "F1 (f,u0,v0): {} conflicts (examples: {:?})",
        f1_conflicts, f1_conflict_examples
    );
    println!("F2 (f,u0,v0,s0): {} conflicts", f2_conflicts);
    println!("F3 (f,u0,v0,r0,s0): {} conflicts", f3_conflicts);
    println!(
        "F4 (f,u0,v0,u>v_w,s0): {} conflicts (examples: {:?})",
        f4_conflicts, f4_conflict_examples
    );
    println!(
        "F_min (f,u0,v0,u>v_w): {} conflicts",
        fmin_conflicts
    );

    // The structural question: F4 includes the comparator, which is
    // mathematically what STEP 2 looks at. Per the algorithm, m_i can only
    // be 1 if either (f ∧ v_w=0) — captured by F1 via u0/v0 — or
    // (f ∧ u0 ∧ NOT v0) [STEP 1] — captured by F1 — or (add_f_step2 ∧ NOT
    // b_f_orig) — which depends on the comparator AND u0,v0,f. So F4 ⊇
    // the full determination set ONLY if b_f_orig depends only on F1.
    //
    // b_f_orig after STEP 1 = a_f XOR m_i = (f ∧ NOT u0) XOR
    //   (f ∧ (u0=1 ∧ v0=0) from STEP 1's m_i toggle) XOR
    //   (m_i entering STEP 1 from STEP 0's toggle).
    //
    // This is a function only of (f, u0, v0, m_i_entering), i.e. F1.
    // So by induction, F4 captures everything m_i needs. Expected: F4 = 0.

    let f4_is_deterministic = f4_conflicts == 0;
    if f4_is_deterministic {
        println!(
            "\n✅ FEASIBILITY CONFIRMED: m_i is a deterministic function of
  (f, u[0], v_w[0], u>v_w, s[0]) at iter start.

This means m_hist (407 qubits) can be REPLACED with an in-iteration
recomputation that needs:
  - quantum: 0 extra persistent qubits
  - per-iter: one comparator reuse (we already compute u>v_w for STEP 2)
  - structural: m_i becomes iter-local, freed at iter end

Qubit-budget implication: Kaliski persistent state drops by 407 qubits
(from 1432 → 1025), making single-inversion Kaliski + aggressive
compression genuinely fit in ≤ 1200 qubits.

HOWEVER: this doesn't say m_i can be *uncomputed* without the forward
copy. The backward pass still needs m_i to choose the branch. Recomputing
m_i DURING backward from (u,v_w,r,s) live-state is only possible if
those registers hold their forward-iteration-start values at backward-iter
start — which they do iff backward is the exact gate-inverse of forward.

The quantum circuit already satisfies this (backward = inverse(forward)),
so m_i can be recomputed at backward-start from the same fingerprint.

=> m_hist elimination is ALGEBRAICALLY POSSIBLE. Implementation risk
remains: phase-correction protocol for the recomputed ancilla (Gidney
MBU pattern applies).");
    } else {
        println!(
            "\n❌ FEASIBILITY REJECTED: F4 is NOT a deterministic function of
(f, u[0], v_w[0], u>v_w, s[0]). {} conflicts observed.

This means a wider fingerprint is needed, which in the extreme requires
the full (u,v_w,r,s) state = 4n bits = exactly what we already store
implicitly in the live state. Recomputing m_i is then equivalent to
re-running the iteration body, which is NOT free.

Implication: Kaliski's m_hist is NOT compressible via simple HRSL-style
recomputation. Kaliski at 1200q requires either
  (a) Kim unconditional execution (+28% Toffoli), or
  (b) Bennett pebbling (multiplicative Toffoli overhead), or
  (c) replacing Kaliski with a non-Kaliski inverter.",
            f4_conflicts
        );
    }

    f4_is_deterministic
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn classical_replay_matches_self() {
        // Sanity: run once, confirm state evolves and we terminate.
        let v_in = random_element(42);
        let (m_hist, snaps) = kaliski_run(v_in, SECP256K1_P, 512);
        assert_eq!(m_hist.len(), 512);
        assert_eq!(snaps.len(), 512);
        // Termination: f should flip to 0 at some point.
        let any_terminated = snaps.iter().any(|s| s.4 == 0);
        assert!(any_terminated, "Kaliski never terminated in 512 iters");
    }

    #[test]
    fn feasibility_test_small() {
        let ok = run_feasibility_test(20, 512);
        // Don't assert the answer — we want the output printed regardless.
        let _ = ok;
    }

    #[test]
    fn feasibility_test_large() {
        // Scale up: 500 inputs × 512 iters = 256k samples. If F4 is truly
        // deterministic, it should still have 0 conflicts.
        let ok = run_feasibility_test(500, 512);
        assert!(ok, "F4 fingerprint failed to be deterministic at large scale");
    }

    /// Helper: replay a single iter and return m_i, a_f, AND the state
    /// just-before-step-10 (= iter-end state since step 10 only touches a_f).
    fn kaliski_iter_with_af(
        u: &mut U256,
        v_w: &mut U256,
        r: &mut U256,
        s: &mut U256,
        f: &mut u8,
        p: U256,
    ) -> (u8, u8) {
        // Replay EXACTLY as kaliski_iter_classical but also return a_f.
        let is_zero = if *v_w == U256::ZERO { 1u8 } else { 0 };
        let mut m_i: u8 = 0;
        if *f == 1 && is_zero == 1 { m_i ^= 1; }
        *f ^= m_i;

        let u0 = (u.as_limbs()[0] & 1) as u8;
        let v0 = (v_w.as_limbs()[0] & 1) as u8;
        let mut a_f: u8 = 0;
        if *f == 1 && u0 == 0 { a_f ^= 1; }
        if *f == 1 && u0 == 1 && v0 == 0 { m_i ^= 1; }
        let b_f = a_f ^ m_i;

        let l_gt = if *u > *v_w { 1u8 } else { 0 };
        let add_f_step2 = (*f & l_gt) as u8;
        let b_not = 1 ^ b_f;
        let delta = add_f_step2 & b_not;
        a_f ^= delta;
        m_i ^= delta;

        if a_f == 1 {
            std::mem::swap(u, v_w);
            std::mem::swap(r, s);
        }

        let add_f_step4 = *f & (1 ^ b_f);
        if add_f_step4 == 1 {
            *v_w = v_w.wrapping_sub(*u);
            let sum = s.wrapping_add(*r);
            *s = sum;
        }

        *v_w = *v_w >> 1;

        let r2 = r.wrapping_add(*r);
        *r = if r2 >= p || r2 < *r { r2.wrapping_sub(p) } else { r2 };

        if a_f == 1 {
            std::mem::swap(u, v_w);
            std::mem::swap(r, s);
        }

        // Return (m_i, a_f before step 10 zeroes it).
        (m_i, a_f)
    }

    #[test]
    fn check_iter_end_plus_af_fingerprint() {
        // Test: at step 10 start (after cswap9, before a_f uncompute),
        // does (f, u[0], v_w[0], u>v_w, s[0], a_f) determine m_i?
        let p = SECP256K1_P;
        use std::collections::HashMap;
        let mut tt: HashMap<u16, (u8, usize, usize)> = HashMap::new();

        for seed in 0..500u64 {
            let v_in = random_element(seed + 1);
            let mut u = p;
            let mut v_w = v_in;
            let mut r = U256::ZERO;
            let mut s = U256::from(1u64);
            let mut f: u8 = 1;
            for _i in 0..512 {
                let (m_i, a_f) = kaliski_iter_with_af(&mut u, &mut v_w, &mut r, &mut s, &mut f, p);
                let u0 = (u.as_limbs()[0] & 1) as u8;
                let v0 = (v_w.as_limbs()[0] & 1) as u8;
                let gt = if u > v_w { 1u8 } else { 0 };
                let s0 = (s.as_limbs()[0] & 1) as u8;
                // Fingerprint: (f, u0, v0, gt, s0, a_f) = 6 bits.
                let k = ((f as u16) << 5) | ((u0 as u16) << 4) | ((v0 as u16) << 3)
                       | ((gt as u16) << 2) | ((s0 as u16) << 1) | (a_f as u16);
                let e = tt.entry(k).or_insert((m_i, 0, 0));
                if m_i == 0 { e.1 += 1; } else { e.2 += 1; }
            }
        }

        let mut conflicts = 0usize;
        println!("\n=== iter-END + a_f fingerprint -> m_i ===");
        let mut ks: Vec<_> = tt.keys().collect();
        ks.sort();
        for k in ks {
            let (_, c0, c1) = tt[k];
            let conflict = c0 > 0 && c1 > 0;
            if conflict { conflicts += 1; }
            println!("  key={:06b}: m=0 {:>7} m=1 {:>7}{}",
                k, c0, c1, if conflict { "  <- CONFLICT" } else { "" });
        }
        println!("TOTAL CONFLICTS: {}", conflicts);
    }

    #[test]
    fn check_iter_end_fingerprint() {
        // Test: at the END of iteration i (after all 10 steps, before iter i+1
        // starts), can we recover m_i as a function of the iter-end state?
        //
        // If YES with a small fingerprint, we can uncompute m_i at iter end.
        // If NO, we need to uncompute m_i at a specific mid-iter point where
        // the start-of-iter fingerprint is still live.
        let p = SECP256K1_P;
        use std::collections::HashMap;
        let mut tt: HashMap<u8, (u8, usize, usize)> = HashMap::new(); // key -> (first_m, count_0, count_1)

        for seed in 0..200u64 {
            let v_in = random_element(seed + 1);
            // Get m_hist and the state AT END of each iter. Run iter-by-iter.
            let mut u = p;
            let mut v_w = v_in;
            let mut r = U256::ZERO;
            let mut s = U256::from(1u64);
            let mut f: u8 = 1;
            for _i in 0..512 {
                let m_i = kaliski_iter_classical(&mut u, &mut v_w, &mut r, &mut s, &mut f, p);
                // State here is END-of-iter.
                let u0 = (u.as_limbs()[0] & 1) as u8;
                let v0 = (v_w.as_limbs()[0] & 1) as u8;
                let gt = if u > v_w { 1u8 } else { 0 };
                let s0 = (s.as_limbs()[0] & 1) as u8;
                // Fingerprint: (f_after, u0_after, v0_after, gt_after, s0_after) -- 5 bits
                let k = (f << 4) | (u0 << 3) | (v0 << 2) | (gt << 1) | s0;
                let e = tt.entry(k).or_insert((m_i, 0, 0));
                if m_i == 0 { e.1 += 1; } else { e.2 += 1; }
            }
        }

        println!("\n=== iter-END fingerprint F=(f,u0,v0,gt,s0) -> m_i ===");
        let mut conflicts = 0usize;
        let mut ks: Vec<_> = tt.keys().collect();
        ks.sort();
        for k in ks {
            let (_, c0, c1) = tt[k];
            let conflict = c0 > 0 && c1 > 0;
            if conflict { conflicts += 1; }
            println!("  key={:05b}: m=0 {:>6} times, m=1 {:>6} times{}",
                k, c0, c1, if conflict { "  <- CONFLICT" } else { "" });
        }
        println!("TOTAL CONFLICTS: {}", conflicts);
    }

    #[test]
    fn verify_minimal_formula() {
        // Check: m_i = f AND u[0] AND (NOT v_w[0] OR gt)
        let p = SECP256K1_P;
        let mut total = 0u64;
        let mut mismatches = 0u64;
        for seed in 0..500u64 {
            let v_in = random_element(seed + 1);
            let (m_hist, snaps) = kaliski_run(v_in, p, 512);
            for i in 0..512 {
                let (u, v_w, _, _, f) = snaps[i];
                let u0 = (u.as_limbs()[0] & 1) as u8;
                let v0 = (v_w.as_limbs()[0] & 1) as u8;
                let gt = if u > v_w { 1u8 } else { 0 };
                let m_i = m_hist[i];
                // Candidate formula: f AND u0 AND (NOT v0 OR gt)
                let pred = f & u0 & ((1 ^ v0) | gt);
                if pred != m_i {
                    mismatches += 1;
                }
                total += 1;
            }
        }
        println!("minimal formula: f AND u0 AND (NOT v0 OR gt)");
        println!("  total={} mismatches={}", total, mismatches);
        assert_eq!(mismatches, 0);
    }

    #[test]
    fn extract_fmin_truth_table() {
        // Print the full 16-entry truth table of F_min -> m_i empirically.
        use std::collections::HashMap;
        let p = SECP256K1_P;
        let mut tt: HashMap<u8, (u8, usize)> = HashMap::new();

        for seed in 0..500u64 {
            let v_in = random_element(seed + 1);
            let (m_hist, snaps) = kaliski_run(v_in, p, 512);
            for i in 0..512 {
                let (u, v_w, _r, _s, f) = snaps[i];
                let u0 = (u.as_limbs()[0] & 1) as u8;
                let v0 = (v_w.as_limbs()[0] & 1) as u8;
                let gt = if u > v_w { 1u8 } else { 0 };
                let m_i = m_hist[i];
                let k = (f << 3) | (u0 << 2) | (v0 << 1) | gt;
                let e = tt.entry(k).or_insert((m_i, 0));
                e.1 += 1;
                assert_eq!(e.0, m_i, "F_min {:04b} produced both m=0 and m=1", k);
            }
        }

        println!("\n=== F_min truth table (empirically observed) ===");
        println!("| f | u[0] | v_w[0] | u>v_w | m_i | samples |");
        println!("|---|------|--------|-------|-----|---------|");
        for k in 0u8..16 {
            let f = (k >> 3) & 1;
            let u0 = (k >> 2) & 1;
            let v0 = (k >> 1) & 1;
            let gt = k & 1;
            if let Some(&(m, c)) = tt.get(&k) {
                println!("| {} |  {}   |   {}    |   {}   |  {}  | {:>6}  |",
                    f, u0, v0, gt, m, c);
            } else {
                println!("| {} |  {}   |   {}    |   {}   |  -  |    0   | (unreachable state)",
                    f, u0, v0, gt);
            }
        }
    }

    #[test]
    fn f1_f2_f3_conflicts_are_narrow() {
        // Extra diagnostic: find the specific fingerprint keys that conflict
        // in F1, F2, F3 so we understand WHY they conflict.
        use std::collections::HashMap;
        let p = SECP256K1_P;
        let mut f1_sets: HashMap<u8, (usize, usize)> = HashMap::new(); // (count_0, count_1)

        for seed in 0..100u64 {
            let v_in = random_element(seed + 1);
            let (m_hist, snaps) = kaliski_run(v_in, p, 512);
            for i in 0..512 {
                let (u, v_w, _r, _s, f) = snaps[i];
                let u0 = (u.as_limbs()[0] & 1) as u8;
                let v0 = (v_w.as_limbs()[0] & 1) as u8;
                let m_i = m_hist[i];
                let k1 = (f << 2) | (u0 << 1) | v0;
                let e = f1_sets.entry(k1).or_insert((0, 0));
                if m_i == 0 { e.0 += 1; } else { e.1 += 1; }
            }
        }

        println!("F1 key -> (count_m=0, count_m=1):");
        let mut keys: Vec<_> = f1_sets.keys().collect();
        keys.sort();
        for k in keys {
            let (c0, c1) = f1_sets[k];
            println!("  key={:03b} (f={}, u0={}, v0={}): m=0 {} times, m=1 {} times{}",
                k, (k >> 2) & 1, (k >> 1) & 1, k & 1, c0, c1,
                if c0 > 0 && c1 > 0 { "  <- CONFLICT" } else { "" });
        }
    }
}
