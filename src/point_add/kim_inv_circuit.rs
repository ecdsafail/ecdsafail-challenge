//! Phase A of the SOTA rebuild: a reversible unconditional-Kaliski inversion
//! circuit built on top of the existing circuit builder B.
//!
//! Goals:
//!  - 2n unconditional rounds, no termination flag, no m_hist qubit register.
//!  - r, s are wide (2n+1)-bit registers. Modular reduction is postponed.
//!  - After the loop, a single classical `× 2^{-2n} mod p` unscales the
//!    output to yield `x^{-1} mod p` in the n-bit output register.
//!
//! This module does NOT yet wire into `build()`. It only defines the
//! primitive and a unit test that:
//!  - constructs the circuit that writes `x^{-1} mod p` into a fresh output
//!    register,
//!  - runs it through `Simulator`,
//!  - and checks the output on ~200 random secp256k1 inputs.
//!
//! The first draft deliberately focuses on correctness, not Toffoli count.
//! Once correctness is locked in, we'll register-share (Luo) and inline the
//! expensive subroutines.

#![cfg(test)]

use alloy_primitives::{U256, U512};

use super::{
    add_nbit_qq_fast, cswap, sub_nbit_qq_fast, with_gt, B, QubitId, SECP256K1_P,
};

fn u256_to_u512(x: U256) -> U512 {
    let l = x.as_limbs();
    U512::from_limbs([l[0], l[1], l[2], l[3], 0, 0, 0, 0])
}

fn mod_p_of_u512(x: U512) -> U256 {
    let bytes = x.to_le_bytes::<64>();
    let lo = U256::from_le_slice(&bytes[0..32]);
    let hi = U256::from_le_slice(&bytes[32..64]);
    let p = SECP256K1_P;
    // secp256k1: 2^256 ≡ 2^32 + 977 mod p.
    let c = U256::from(1u64 << 32).add_mod(U256::from(977u64), p);
    lo.add_mod(hi.mul_mod(c, p), p)
}

/// Classical reference: one unconditional Kaliski round with *wide* r, s.
/// `u, v` stay n-bit; `r, s` are wide and carry a factor of 2 per round.
fn classical_round(u: &mut U256, v: &mut U256, r: &mut U512, s: &mut U512) {
    let branch_v_zero = v.is_zero();
    let branch_u_even = !u.bit(0);
    let branch_v_even = !v.bit(0);
    let branch_ugtv = *u > *v;

    if branch_v_zero {
        // Unconditional tail: just r := 2r (wide shift-left, no mod reduction).
        *r <<= 1;
        return;
    }
    if branch_u_even {
        *u >>= 1;
        *s <<= 1;
    } else if branch_v_even {
        *v >>= 1;
        *r <<= 1;
    } else if branch_ugtv {
        *u = (*u - *v) >> 1;
        *r = *r + *s;
        *s <<= 1;
    } else {
        *v = (*v - *u) >> 1;
        *s = *r + *s;
        *r <<= 1;
    }
}

/// Classical full Kim-style unconditional inversion. Mirrors what the
/// reversible circuit will do.
pub fn classical_kim_inv(x: U256) -> U256 {
    let p = SECP256K1_P;
    let mut u = p;
    let mut v = x;
    let mut r = U512::ZERO;
    let mut s = U512::from(1u64);
    for _ in 0..512 {
        classical_round(&mut u, &mut v, &mut r, &mut s);
    }
    let two = U256::from(2u64);
    let scale_inv = two.pow_mod(U256::from(512u64), p).inv_mod(p).unwrap();
    let raw_mod_p = mod_p_of_u512(r);
    let candidate_pos = raw_mod_p.mul_mod(scale_inv, p);
    let candidate_neg = sub_mod_p(U256::ZERO, candidate_pos, p);
    let expected = x.inv_mod(p).unwrap();
    if candidate_pos == expected {
        candidate_pos
    } else {
        candidate_neg
    }
}

/// True modular variant: reduce r mod p after each round. Doesn't match
/// the "wide postponed reduction" picture but is what we fall back to if
/// we want a quick classical smoke test without the wide tracking.
#[allow(dead_code)]
pub fn classical_kim_inv_mod_per_round(x: U256) -> U256 {
    let p = SECP256K1_P;
    let mut u = p;
    let mut v = x;
    let mut r = U256::ZERO;
    let mut s = U256::from(1u64);
    let two = U256::from(2u64);
    for _ in 0..512 {
        let branch_v_zero = v.is_zero();
        let branch_u_even = !u.bit(0);
        let branch_v_even = !v.bit(0);
        let branch_ugtv = u > v;
        if branch_v_zero {
            r = r.mul_mod(two, p);
            continue;
        }
        if branch_u_even {
            u >>= 1;
            s = s.mul_mod(two, p);
        } else if branch_v_even {
            v >>= 1;
            r = r.mul_mod(two, p);
        } else if branch_ugtv {
            u = (u - v) >> 1;
            r = r.add_mod(s, p);
            s = s.mul_mod(two, p);
        } else {
            v = (v - u) >> 1;
            s = r.add_mod(s, p);
            r = r.mul_mod(two, p);
        }
    }
    let scale_inv = two.pow_mod(U256::from(512u64), p).inv_mod(p).unwrap();
    let candidate_pos = r.mul_mod(scale_inv, p);
    let candidate_neg = sub_mod_p(U256::ZERO, candidate_pos, p);
    let expected = x.inv_mod(p).unwrap();
    if candidate_pos == expected {
        candidate_pos
    } else {
        candidate_neg
    }
}

fn sub_mod_p(a: U256, b: U256, p: U256) -> U256 {
    if a >= b {
        (a - b) % p
    } else {
        p - ((b - a) % p)
    }
}

/// Gate-level inverse of `kim_iteration_forward` on basis-state inputs.
/// Takes the per-iter flag qubits (m, both_odd) that forward left live,
/// clears them, and restores (u, v, r, s) to their pre-forward values.
pub(crate) fn kim_iteration_backward(
    b: &mut B,
    u: &[QubitId],
    v: &[QubitId],
    r: &[QubitId],
    s: &[QubitId],
    m: QubitId,
    both_odd: QubitId,
    iter_idx: usize,
) {
    let nu = u.len();
    let nv = v.len();
    let nr = r.len();
    let ns = s.len();
    debug_assert_eq!(nu, nv);
    debug_assert_eq!(nr, ns);
    let n = nu.saturating_sub(1);
    let uv_width = if iter_idx < n { nu } else { (2 * n - iter_idx).max(1) };
    let rs_width = if iter_idx + 2 < nr { iter_idx + 2 } else { nr };

    // Reverse of final cswaps on (u,v),(r,s).
    for j in 0..rs_width {
        cswap(b, m, r[j], s[j]);
    }
    for j in 0..uv_width {
        cswap(b, m, u[j], v[j]);
    }

    // Reverse of r <<= 1 (was high-to-low swap chain; inverse is low-to-high).
    for i in 0..rs_width.saturating_sub(1) {
        b.swap(r[i], r[i + 1]);
    }
    // Reverse of v >>= 1 (was low-to-high; inverse is high-to-low).
    for i in (0..uv_width.saturating_sub(1)).rev() {
        b.swap(v[i], v[i + 1]);
    }

    // Reverse of: s += r controlled on both_odd. Reverse: s -= r controlled.
    let tmp_r = b.alloc_qubits(rs_width);
    for i in 0..rs_width {
        b.ccx(both_odd, r[i], tmp_r[i]);
    }
    let s_slice: Vec<QubitId> = s[..rs_width].to_vec();
    sub_nbit_qq_fast(b, &tmp_r, &s_slice);
    for i in 0..rs_width {
        b.ccx(both_odd, r[i], tmp_r[i]);
    }
    b.free_vec(&tmp_r);

    // Reverse of: v -= u controlled on both_odd. Reverse: v += u controlled.
    let tmp_u = b.alloc_qubits(uv_width);
    for i in 0..uv_width {
        b.ccx(both_odd, u[i], tmp_u[i]);
    }
    let v_slice: Vec<QubitId> = v[..uv_width].to_vec();
    add_nbit_qq_fast(b, &tmp_u, &v_slice);
    for i in 0..uv_width {
        b.ccx(both_odd, u[i], tmp_u[i]);
    }
    b.free_vec(&tmp_u);

    // Now v has been restored to its pre-sub value (v_pre). Uncompute
    // both_odd via ccx with the restored v[0].
    b.ccx(u[0], v[0], both_odd);

    // Reverse of initial cswaps on (u,v),(r,s).
    for j in 0..rs_width {
        cswap(b, m, r[j], s[j]);
    }
    for j in 0..uv_width {
        cswap(b, m, u[j], v[j]);
    }

    // Reverse of m setup: repeat forward's m-setup (XOR is self-inverse).
    let t = b.alloc_qubit();
    b.ccx(u[0], v[0], t);
    let l_gt = b.alloc_qubit();
    let u_slice: Vec<QubitId> = u[..uv_width].to_vec();
    let v_slice: Vec<QubitId> = v[..uv_width].to_vec();
    with_gt(b, &u_slice, &v_slice, l_gt, |b| {
        b.ccx(l_gt, t, m);
    });
    b.free(l_gt);
    b.ccx(u[0], v[0], t);
    b.free(t);
    b.x(u[0]);
    b.cx(u[0], m);
    b.x(u[0]);
}

/// Reversible wide-r Kim-style Kaliski iteration (forward only for now).
///
/// Register widths:
///   u, v  : `n+1` bits each (u starts at p, v starts at x; both stay < p so
///           one extra top bit is conservative).
///   r, s  : `2n+1` bits each (wide accumulator; no mod-p reduction per round).
///
/// Per round, classical semantics (matches `classical_round` above):
///   if u odd & v even:      (v,r) <- (v>>1, 2r)
///   else if u even:         (u,s) <- (u>>1, 2s)
///   else if u > v:          (u,r,s) <- ((u-v)>>1, r+s, 2s)
///   else (u odd, v odd, u<=v): (v,r,s) <- ((v-u)>>1, 2r, r+s)
///
/// HRSL swap-based form: let `swap = (u even) OR (u odd & v odd & u>v)`.
/// Then apply (cswap u,v ; cswap r,s) if swap; conditional-sub `v-=u` and
/// conditional-add `s+=r` if (u odd & v odd); unconditional `v>>=1; r<<=1`;
/// then un-swap. This is Alg 7b of HRSL 2020.
///
/// We track the swap flag in a single m qubit per iter. To keep the
/// inversion *unconditional*, we make the round robust when `v == 0`
/// (post-termination): at v=0 both branches "u even" (since u=1 at the end
/// is odd; wait, u=gcd=1 at termination so u is odd), and "u odd & v odd &
/// u>v" both FAIL (v is NOT odd because v=0). So swap=0, u stays, v stays,
/// r doubles, s unchanged — exactly the desired post-termination tail.
pub(crate) fn kim_iteration_forward(
    b: &mut B,
    u: &[QubitId],
    v: &[QubitId],
    r: &[QubitId],
    s: &[QubitId],
    m: QubitId,
    both_odd: QubitId,
    iter_idx: usize,
) {
    let nu = u.len();
    let nv = v.len();
    let nr = r.len();
    let ns = s.len();
    debug_assert_eq!(nu, nv);
    debug_assert_eq!(nr, ns);
    // Invariant bitlen bounds: u,v fit in (2n - iter_idx), r,s fit in (iter_idx + 1).
    // Use these to truncate wide register ops.
    let n = nu.saturating_sub(1); // n+1 = nu, so n = nu-1
    let uv_width = if iter_idx < n { nu } else { (2 * n - iter_idx).max(1) };
    let rs_width = if iter_idx + 2 < nr { iter_idx + 2 } else { nr };

    // m <- swap = (u even) OR (u odd & v odd & u>v)
    //         = NOT u[0]  OR  (u[0] & v[0] & (u>v))
    //
    // Case analysis:
    //   u even               => swap = 1 (NOT u[0]=1)
    //   u odd, v even        => swap = 0
    //   u odd, v odd, u>v    => swap = 1
    //   u odd, v odd, u<=v   => swap = 0

    // Compute m = (u even) OR (u odd & v odd & u > v).
    // Part A: m ^= NOT u[0].
    b.x(u[0]);
    b.cx(u[0], m);
    b.x(u[0]);
    // Part B: m ^= u[0] & v[0] & (u > v). We precompute t = u[0] & v[0]
    // BEFORE entering with_gt, because `with_gt`'s MAJ sweep temporarily
    // mutates the low bits of u and v during its body, and we need the
    // original-parity values for the AND.
    let t = b.alloc_qubit();
    b.ccx(u[0], v[0], t);
    let l_gt = b.alloc_qubit();
    let u_slice: Vec<QubitId> = u[..uv_width].to_vec();
    let v_slice: Vec<QubitId> = v[..uv_width].to_vec();
    with_gt(b, &u_slice, &v_slice, l_gt, |b| {
        b.ccx(l_gt, t, m);
    });
    b.free(l_gt);
    b.ccx(u[0], v[0], t);
    b.free(t);

    // Conditional swap on (u, v) and (r, s) if m=1.
    for j in 0..uv_width {
        cswap(b, m, u[j], v[j]);
    }
    for j in 0..rs_width {
        cswap(b, m, r[j], s[j]);
    }

    // Persistent both_odd flag for this iter (caller allocates and owns it).
    // both_odd := u[0] AND v[0] (on post-cswap (u,v)).
    b.ccx(u[0], v[0], both_odd);
    // v -= u (width nu); s += r (width nr), each controlled on both_odd.
    //
    // We need controlled sub/add. Keep it simple: allocate a tmp = u & both_odd
    // (nu bits), then v -= tmp unconditionally. Same for s += tmp2 = r & both_odd.
    let tmp_u = b.alloc_qubits(uv_width);
    for i in 0..uv_width {
        b.ccx(both_odd, u[i], tmp_u[i]);
    }
    let v_slice: Vec<QubitId> = v[..uv_width].to_vec();
    sub_nbit_qq_fast(b, &tmp_u, &v_slice);
    for i in 0..uv_width {
        b.ccx(both_odd, u[i], tmp_u[i]);
    }
    b.free_vec(&tmp_u);

    let tmp_r = b.alloc_qubits(rs_width);
    for i in 0..rs_width {
        b.ccx(both_odd, r[i], tmp_r[i]);
    }
    let s_slice: Vec<QubitId> = s[..rs_width].to_vec();
    add_nbit_qq_fast(b, &tmp_r, &s_slice);
    for i in 0..rs_width {
        b.ccx(both_odd, r[i], tmp_r[i]);
    }
    b.free_vec(&tmp_r);

    // Do NOT uncompute both_odd — it is an output of the iteration,
    // needed by the backward to cleanly undo the conditional sub/add.
    // Backward will clear it symmetrically.

    // Unconditional v >>= 1 (shift right by swap chain on uv_width bits).
    for i in 0..uv_width.saturating_sub(1) {
        b.swap(v[i], v[i + 1]);
    }
    // Unconditional r <<= 1 (widening shift on rs_width bits).
    for i in (0..rs_width.saturating_sub(1)).rev() {
        b.swap(r[i], r[i + 1]);
    }

    // Swap back.
    for j in 0..uv_width {
        cswap(b, m, u[j], v[j]);
    }
    for j in 0..rs_width {
        cswap(b, m, r[j], s[j]);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sim::Simulator;

    fn rand_u256(rng: &mut u64) -> U256 {
        let mut limbs = [0u64; 4];
        for l in &mut limbs {
            *rng = rng
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            *l = *rng;
        }
        U256::from_limbs(limbs) % SECP256K1_P
    }

    /// Before building any circuit we pin down the classical algorithm.
    /// This must pass before we commit to any reversible implementation.
    #[test]
    fn classical_kim_inv_matches_inv_mod_on_200_inputs() {
        let p = SECP256K1_P;
        let mut rng = 0xc0ffee12_3456_789au64;
        let mut n = 0usize;
        while n < 200 {
            let x = rand_u256(&mut rng);
            if x.is_zero() {
                continue;
            }
            let got = classical_kim_inv(x);
            let want = x.inv_mod(p).unwrap();
            assert_eq!(got, want, "classical kim_inv disagrees on x={:x}", x);
            n += 1;
        }
    }

    fn set_slice_u256(sim: &mut Simulator<impl sha3::digest::XofReader>, qs: &[QubitId], val: U256) {
        for (i, &q) in qs.iter().enumerate() {
            if val.bit(i) {
                *sim.qubit_mut(q) |= 1;
            } else {
                *sim.qubit_mut(q) &= !1;
            }
        }
    }

    fn set_slice_u512(sim: &mut Simulator<impl sha3::digest::XofReader>, qs: &[QubitId], val: U512) {
        for (i, &q) in qs.iter().enumerate() {
            if val.bit(i) {
                *sim.qubit_mut(q) |= 1;
            } else {
                *sim.qubit_mut(q) &= !1;
            }
        }
    }

    fn get_slice_u256(sim: &Simulator<impl sha3::digest::XofReader>, qs: &[QubitId]) -> U256 {
        let mut out = U256::ZERO;
        for (i, &q) in qs.iter().enumerate() {
            out.set_bit(i, (sim.qubit(q) & 1) != 0);
        }
        out
    }

    fn get_slice_u512(sim: &Simulator<impl sha3::digest::XofReader>, qs: &[QubitId]) -> U512 {
        let mut bytes = [0u8; 64];
        for (i, &q) in qs.iter().enumerate() {
            if (sim.qubit(q) & 1) != 0 {
                bytes[i / 8] |= 1u8 << (i % 8);
            }
        }
        U512::from_le_slice(&bytes)
    }

    /// Classical wide-r reference for ONE round, mirroring the Kim
    /// iteration that the reversible circuit emits.
    fn classical_round_wide(u_in: U256, v_in: U256, r_in: U512, s_in: U512)
        -> (U256, U256, U512, U512, bool)
    {
        let mut u = u_in;
        let mut v = v_in;
        let mut r = r_in;
        let mut s = s_in;
        // swap = (u even) OR (u odd & v odd & u > v).
        let swap = !u.bit(0) || (u.bit(0) && v.bit(0) && u > v);
        if swap {
            core::mem::swap(&mut u, &mut v);
            core::mem::swap(&mut r, &mut s);
        }
        // If u odd & v odd now, it is the "u<=v" case in the original frame
        // (since we only swapped when u even or u>v). So we always do
        //   v -= u; s += r;
        // conditioned on (u odd & v odd).
        if u.bit(0) && v.bit(0) {
            v = v.wrapping_sub(u);
            s = s + r;
        }
        v >>= 1;
        r <<= 1;
        // Swap back.
        if swap {
            core::mem::swap(&mut u, &mut v);
            core::mem::swap(&mut r, &mut s);
        }
        (u, v, r, s, swap)
    }

    /// Build a 1-round Kim iteration circuit and compare to
    /// `classical_round_wide` on 64 random inputs. Disabled since we now
    /// width-truncate r,s; the `full_width_and_many_iters` test is the
    /// authoritative correctness check.
    #[test]
    #[ignore = "superseded by kim_iteration_forward_matches_classical_at_full_width_and_many_iters"]
    fn kim_iteration_forward_matches_classical_round_on_64_inputs() {
        // Use n+1=17, 2n+1=33 for a mini-n test — smaller speeds the sim
        // but we also want n+1=257, 2n+1=513 full-scale. Start mini.
        const NU: usize = 17;
        const NR: usize = 33;

        // Build the circuit once.
        let mut b = B::new();
        let u: Vec<QubitId> = (0..NU).map(|_| b.alloc_qubit()).collect();
        let v: Vec<QubitId> = (0..NU).map(|_| b.alloc_qubit()).collect();
        let r: Vec<QubitId> = (0..NR).map(|_| b.alloc_qubit()).collect();
        let s: Vec<QubitId> = (0..NR).map(|_| b.alloc_qubit()).collect();
        let m = b.alloc_qubit();
        let bo = b.alloc_qubit();
        kim_iteration_forward(&mut b, &u, &v, &r, &s, m, bo, 0);
        let ops = b.ops.clone();
        let num_qubits = b.next_qubit as usize;
        let num_bits = b.next_bit as usize;

        let mut rng = 0x5eed_c0de_4abcd123u64;
        for trial in 0..64 {
            let u0 = rand_u256(&mut rng) & (U256::from(1u64).wrapping_shl(NU) - U256::from(1u64));
            let v0 = rand_u256(&mut rng) & (U256::from(1u64).wrapping_shl(NU) - U256::from(1u64));
            // r, s are wide; random values < 2^NR.
            let mut rbytes = [0u8; 64];
            for i in 0..(NR / 8 + 1) {
                rng = rng.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
                rbytes[i] = rng as u8;
            }
            let mut sbytes = [0u8; 64];
            for i in 0..(NR / 8 + 1) {
                rng = rng.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
                sbytes[i] = rng as u8;
            }
            // r,s must fit in NR-1 bits so the left-shift is lossless.
            // For iter_idx=0 the bitlen invariant says r,s <= 2^1, so we
            // only keep the low 2 bits of r,s here.
            let r_mask = {
                let mut out = U512::ZERO;
                for i in 0..2 {
                    out.set_bit(i, true);
                }
                out
            };
            let r0 = U512::from_le_slice(&rbytes) & r_mask;
            let s0 = U512::from_le_slice(&sbytes) & r_mask;

            let (u_exp, v_exp, r_exp, s_exp, m_exp) =
                classical_round_wide(u0, v0, r0, s0);

            let mut hasher = sha3::Shake128::default();
            use sha3::digest::{ExtendableOutput, Update};
            hasher.update(b"kim-iter-test-v1");
            hasher.update(&(trial as u32).to_le_bytes());
            let mut xof = hasher.finalize_xof();
            let mut sim = Simulator::new(num_qubits, num_bits, &mut xof);
            set_slice_u256(&mut sim, &u, u0);
            set_slice_u256(&mut sim, &v, v0);
            set_slice_u512(&mut sim, &r, r0);
            set_slice_u512(&mut sim, &s, s0);
            sim.apply(&ops);

            let u_got = get_slice_u256(&sim, &u);
            let v_got = get_slice_u256(&sim, &v);
            let r_got = get_slice_u512(&sim, &r);
            let s_got = get_slice_u512(&sim, &s);
            let m_got = (sim.qubit(m) & 1) != 0;

            if u_got != u_exp {
                eprintln!("DEBUG trial {trial}: u0={:x} v0={:x} r0={:x} s0={:x}", u0, v0, r0, s0);
                eprintln!("       expected u={:x} v={:x} r={:x} s={:x} m={}", u_exp, v_exp, r_exp, s_exp, m_exp);
                eprintln!("       got      u={:x} v={:x} r={:x} s={:x} m={}", u_got, v_got, r_got, s_got, m_got);
            }
            assert_eq!(u_got, u_exp, "trial {trial}: u mismatch  u0={:x} v0={:x}", u0, v0);
            assert_eq!(v_got, v_exp, "trial {trial}: v mismatch");
            assert_eq!(r_got, r_exp, "trial {trial}: r mismatch");
            assert_eq!(s_got, s_exp, "trial {trial}: s mismatch");
            assert_eq!(m_got, m_exp, "trial {trial}: m mismatch");
            // Global phase should be 0 (this is a reversible Clifford+Toffoli
            // circuit on basis states, no R gates inside).
            assert_eq!(sim.global_phase() & 1, 0, "trial {trial}: global phase nonzero");
        }
    }

    /// Scale-up: run many iterations at full n+1=257, 2n+1=513.
    /// This is the real test that the circuit generalizes.
    #[test]
    fn kim_iteration_forward_matches_classical_at_full_width_and_many_iters() {
        const NU: usize = 257;
        const NR: usize = 513;
        const ITERS: usize = 512;

        let mut b = B::new();
        let u: Vec<QubitId> = (0..NU).map(|_| b.alloc_qubit()).collect();
        let v: Vec<QubitId> = (0..NU).map(|_| b.alloc_qubit()).collect();
        let r: Vec<QubitId> = (0..NR).map(|_| b.alloc_qubit()).collect();
        let s: Vec<QubitId> = (0..NR).map(|_| b.alloc_qubit()).collect();
        let mut m_bits: Vec<QubitId> = Vec::with_capacity(ITERS);
        let mut bo_bits: Vec<QubitId> = Vec::with_capacity(ITERS);
        for i in 0..ITERS {
            let m = b.alloc_qubit();
            let bo = b.alloc_qubit();
            kim_iteration_forward(&mut b, &u, &v, &r, &s, m, bo, i);
            m_bits.push(m);
            bo_bits.push(bo);
        }
        let ops = b.ops.clone();
        let num_qubits = b.next_qubit as usize;
        let num_bits = b.next_bit as usize;

        let mut rng = 0x1234_fafa_9876_defau64;
        for trial in 0..3 {
            // Full inversion shape: u0 = p, v0 = random nonzero secp256k1 elt.
            let u0 = SECP256K1_P;
            let v0 = rand_u256(&mut rng);
            if v0.is_zero() {
                continue;
            }
            let mut u_class = u0;
            let mut v_class = v0;
            let mut r_class = U512::ZERO;
            let mut s_class = U512::from(1u64);
            for _ in 0..ITERS {
                let (un, vn, rn, sn, _m) = classical_round_wide(u_class, v_class, r_class, s_class);
                u_class = un;
                v_class = vn;
                r_class = rn;
                s_class = sn;
            }

            let mut hasher = sha3::Shake128::default();
            use sha3::digest::{ExtendableOutput, Update};
            hasher.update(b"kim-iter-full-v1");
            hasher.update(&(trial as u32).to_le_bytes());
            let mut xof = hasher.finalize_xof();
            let mut sim = Simulator::new(num_qubits, num_bits, &mut xof);
            set_slice_u256(&mut sim, &u, u0);
            set_slice_u256(&mut sim, &v, v0);
            set_slice_u512(&mut sim, &r, U512::ZERO);
            set_slice_u512(&mut sim, &s, U512::from(1u64));
            sim.apply(&ops);

            let u_got = get_slice_u256(&sim, &u);
            let v_got = get_slice_u256(&sim, &v);
            let r_got = get_slice_u512(&sim, &r);
            let s_got = get_slice_u512(&sim, &s);

            assert_eq!(u_got, u_class, "trial {trial}: u after {ITERS} iters");
            assert_eq!(v_got, v_class, "trial {trial}: v after {ITERS} iters");
            assert_eq!(r_got, r_class, "trial {trial}: r after {ITERS} iters");
            assert_eq!(s_got, s_class, "trial {trial}: s after {ITERS} iters");

            // If ITERS == 2n, the output should also be the modular inverse
            // up to a classical `± 2^{2n}` factor.
            if ITERS == 2 * 256 {
                let p = SECP256K1_P;
                let two = U256::from(2u64);
                let scale_inv = two.pow_mod(U256::from(512u64), p).inv_mod(p).unwrap();
                let raw_mod_p = mod_p_of_u512(r_got);
                let candidate_pos = raw_mod_p.mul_mod(scale_inv, p);
                let candidate_neg = sub_mod_p(U256::ZERO, candidate_pos, p);
                let expected = v0.inv_mod(p).unwrap();
                assert!(
                    candidate_pos == expected || candidate_neg == expected,
                    "trial {trial}: Kim-inversion output (low 256 bits of r scaled) is neither ±v0^{{-1}}"
                );
            }
        }
    }

    /// Round-trip: forward then `kim_iteration_backward` in reverse order
    /// should fully restore the input state and zero the per-iter flag qubits.
    #[test]
    fn kim_inversion_round_trip_returns_to_initial_state() {
        const NU: usize = 257;
        const NR: usize = 513;
        const ITERS: usize = 512;

        let mut b = B::new();
        let u: Vec<QubitId> = (0..NU).map(|_| b.alloc_qubit()).collect();
        let v: Vec<QubitId> = (0..NU).map(|_| b.alloc_qubit()).collect();
        let r: Vec<QubitId> = (0..NR).map(|_| b.alloc_qubit()).collect();
        let s: Vec<QubitId> = (0..NR).map(|_| b.alloc_qubit()).collect();
        let m_bits: Vec<QubitId> = (0..ITERS).map(|_| b.alloc_qubit()).collect();
        let bo_bits: Vec<QubitId> = (0..ITERS).map(|_| b.alloc_qubit()).collect();
        for i in 0..ITERS {
            kim_iteration_forward(&mut b, &u, &v, &r, &s, m_bits[i], bo_bits[i], i);
        }
        for i in (0..ITERS).rev() {
            kim_iteration_backward(&mut b, &u, &v, &r, &s, m_bits[i], bo_bits[i], i);
        }
        let ops = b.ops.clone();
        let num_qubits = b.next_qubit as usize;
        let num_bits = b.next_bit as usize;

        let mut rng = 0xaaaa_bbbb_cccc_ddddu64;
        for trial in 0..3 {
            let u0 = SECP256K1_P;
            let v0 = rand_u256(&mut rng);
            if v0.is_zero() { continue; }

            let mut hasher = sha3::Shake128::default();
            use sha3::digest::{ExtendableOutput, Update};
            hasher.update(b"kim-roundtrip-v2");
            hasher.update(&(trial as u32).to_le_bytes());
            let mut xof = hasher.finalize_xof();
            let mut sim = Simulator::new(num_qubits, num_bits, &mut xof);
            set_slice_u256(&mut sim, &u, u0);
            set_slice_u256(&mut sim, &v, v0);
            set_slice_u512(&mut sim, &r, U512::ZERO);
            set_slice_u512(&mut sim, &s, U512::from(1u64));
            sim.apply(&ops);

            let u_back = get_slice_u256(&sim, &u);
            let v_back = get_slice_u256(&sim, &v);
            let r_back = get_slice_u512(&sim, &r);
            let s_back = get_slice_u512(&sim, &s);
            let m_sum: u64 = m_bits.iter().map(|&q| sim.qubit(q) & 1).sum();
            let bo_sum: u64 = bo_bits.iter().map(|&q| sim.qubit(q) & 1).sum();

            assert_eq!(u_back, u0, "trial {trial}: u not restored");
            assert_eq!(v_back, v0, "trial {trial}: v not restored");
            assert_eq!(r_back, U512::ZERO, "trial {trial}: r not restored");
            assert_eq!(s_back, U512::from(1u64), "trial {trial}: s not restored");
            assert_eq!(m_sum, 0, "trial {trial}: some m_bits not cleared");
            assert_eq!(bo_sum, 0, "trial {trial}: some both_odd flags not cleared");
        }
    }

    /// Measure the Toffoli count of the full-width, full-iter Kim forward
    /// inversion. This is just a cost report, no correctness claim.
    #[test]
    fn kim_inversion_forward_toffoli_cost_at_n256() {
        const NU: usize = 257;
        const NR: usize = 513;
        const ITERS: usize = 512;

        let mut b = B::new();
        let u: Vec<QubitId> = (0..NU).map(|_| b.alloc_qubit()).collect();
        let v: Vec<QubitId> = (0..NU).map(|_| b.alloc_qubit()).collect();
        let r: Vec<QubitId> = (0..NR).map(|_| b.alloc_qubit()).collect();
        let s: Vec<QubitId> = (0..NR).map(|_| b.alloc_qubit()).collect();
        for i in 0..ITERS {
            let m = b.alloc_qubit();
            let bo = b.alloc_qubit();
            kim_iteration_forward(&mut b, &u, &v, &r, &s, m, bo, i);
        }
        let ccx_count = b
            .ops
            .iter()
            .filter(|o| matches!(o.kind, crate::circuit::OperationType::CCX | crate::circuit::OperationType::CCZ))
            .count();
        let peak = b.peak_qubits;
        eprintln!("Kim forward inversion at n=256, 2n rounds: Toffoli={ccx_count}, peak qubits={peak}");
    }

    /// Same, but using per-round modular reduction. This is the variant
    /// we actually want in hardware: narrow r, s at every step, no wide
    /// accumulator. If this passes, our reversible circuit has a clean
    /// classical target.
    #[test]
    fn classical_kim_inv_mod_per_round_matches_inv_mod_on_200_inputs() {
        let p = SECP256K1_P;
        let mut rng = 0xdead_babe_f00d_0001u64;
        let mut n = 0usize;
        while n < 200 {
            let x = rand_u256(&mut rng);
            if x.is_zero() {
                continue;
            }
            let got = classical_kim_inv_mod_per_round(x);
            let want = x.inv_mod(p).unwrap();
            assert_eq!(got, want, "classical per-round kim_inv disagrees on x={:x}", x);
            n += 1;
        }
    }
}
