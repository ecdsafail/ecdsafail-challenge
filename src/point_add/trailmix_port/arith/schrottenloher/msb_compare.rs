//! Approximate top-k-bit comparators per Schrottenloher 2026
//! Sec. 4 / Schrottenloher's `qc_clt_uint` (controlled less-than).
//!
//! `a < b` is decided on the top `k` MSBs only. Used to drive the
//! reduction branch in approximate mod-add. Failure mode: when the
//! top-k MSBs of `a` and `b` are equal, the comparator returns 0
//! (i.e. "false"), regardless of the full-precision result. On uniform
//! inputs this happens with probability ~2^-k per call.
//!
//! For secp256k1 with padding = 30, k = 30 gives ~10^-9 per call.
use crate::point_add::trailmix_port::circuit::{Circuit, QReg};

/// STRICT controlled-less-than on top-k MSBs.
///
/// `target ^= ctrl AND (a_msbs < b_msbs)`. STRICT: at top-k equality
/// `a_top == b_top`, the predicate is FALSE (no XOR). This matches
/// Schrottenloher Alg 10's `qc_clt_uint` semantic; using `>=` would
/// spuriously fire in the no-overflow case where `a == b` at top-k,
/// leaving `anc_y` dirty for the next iteration.
///
/// Implementation: uses `compare_geq_physical_middle` so the ccx
/// touching `target` is the LAST gate on `target` before the
/// uncompute UMA pass — which contains no further allocations. This
/// matters: an interleaved alloc after `target`'s last touch trips
/// the strict-dealloc retention check when `target` (a reused
/// per-iter ancilla) drops at end-of-function.
pub fn controlled_lt_msbs(
    circ: &mut Circuit,
    ctrl: &QReg,
    a: &[QReg],
    b: &[QReg],
    k: usize,
    target: &QReg,
) {
    assert!(k <= a.len() && k <= b.len(), "k must fit in both registers");
    let a_top = &a[a.len() - k..];
    let b_top = &b[b.len() - k..];
    let lt_flag = circ.alloc_qreg("lt_flag");
    crate::point_add::trailmix_port::arith::compare::compare_geq_physical_middle(circ, a_top, b_top, &lt_flag, |c, flag| {
        // Inside: flag = (a >= b). target ^= ctrl AND (a < b) =
        // ctrl AND NOT flag.
        c.x(flag);
        c.ccx(ctrl, flag, target);
        c.x(flag);
    });
    // lt_flag returns to 0 after middle's UMA; drop and drain.
    drop(lt_flag);
}

/// UNCONTROLLED top-k less-than: `target ^= (a_msbs < b_msbs)`.
///
/// The ctrl-free form of [`controlled_lt_msbs`], used by the
/// unconditional pseudo-Mersenne mod-add cleanup. STRICT at top-k
/// equality (predicate FALSE), same as the controlled version. The
/// inner `cx(flag, target)` (flag = a >= b) gives `target ^= NOT flag
/// = target ^= (a < b)`, with `flag` restored to 0 by the comparator's
/// UMA pass.
pub fn lt_msbs(circ: &mut Circuit, a: &[QReg], b: &[QReg], k: usize, target: &QReg) {
    assert!(k <= a.len() && k <= b.len(), "k must fit in both registers");
    let a_top = &a[a.len() - k..];
    let b_top = &b[b.len() - k..];
    let lt_flag = circ.alloc_qreg("lt_flag");
    crate::point_add::trailmix_port::arith::compare::compare_geq_physical_middle(circ, a_top, b_top, &lt_flag, |c, flag| {
        c.x(flag);
        c.cx(flag, target);
        c.x(flag);
    });
    drop(lt_flag);
}

/// Gidney measure-uncompute variant of [`controlled_lt_msbs`]: identical
/// semantics, but the top-`k` comparator uses
/// [`crate::point_add::trailmix_port::arith::compare::compare_geq_gidney_middle`] (n Toffoli + measurement-erased
/// carries, peak +(k+1)) instead of the 2n Cuccaro comparator. Use where the
/// ancilla headroom exists (the GCD loop): halves the per-call comparator
/// Toffoli at the cost of +(k+1) peak qubits.
pub fn controlled_lt_msbs_gidney(
    circ: &mut Circuit,
    ctrl: &QReg,
    a: &[QReg],
    b: &[QReg],
    k: usize,
    target: &QReg,
) {
    assert!(k <= a.len() && k <= b.len(), "k must fit in both registers");
    let a_top = &a[a.len() - k..];
    let b_top = &b[b.len() - k..];
    let lt_flag = circ.alloc_qreg("lt_flag");
    crate::point_add::trailmix_port::arith::compare::compare_geq_gidney_middle(circ, a_top, b_top, &lt_flag, |c, flag| {
        c.x(flag);
        c.ccx(ctrl, flag, target);
        c.x(flag);
    });
    drop(lt_flag);
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::point_add::trailmix_port::circuit::Circuit;

    /// Verify that controlled_lt_msbs agrees with the classical
    /// `(a >> shift) < (b >> shift)` on a large random set, with
    /// disagreements bounded by the equality-tail probability.
    #[test]
    fn controlled_lt_msbs_random_agreement() {
        use rand::rngs::StdRng;
        use rand::{Rng, SeedableRng};

        let n = 64usize;
        let k = 16usize; // top 16 of 64 bits
        let shift = n - k;

        let mut circ = Circuit::new();
        let a = circ.alloc_qreg_bits("a", n);
        let b = circ.alloc_qreg_bits("b", n);
        let ctrl = circ.alloc_qreg("ctrl");
        let flag = circ.alloc_qreg("flag");

        let mut rng = StdRng::seed_from_u64(0xc31d0_a017_b1ab);
        let mut a_vals = [0u64; 64];
        let mut b_vals = [0u64; 64];
        let mut ctrl_vals = [false; 64];
        for shot in 0..64 {
            a_vals[shot] = rng.gen();
            b_vals[shot] = rng.gen();
            ctrl_vals[shot] = (shot & 1) == 0;
            let abytes = a_vals[shot].to_le_bytes();
            let bbytes = b_vals[shot].to_le_bytes();
            circ.sim_load_reg_bytes_shot(&a, &abytes, shot);
            circ.sim_load_reg_bytes_shot(&b, &bbytes, shot);
            if ctrl_vals[shot] {
                circ.sim_load_reg_bytes_shot(std::slice::from_ref(&ctrl), &[1u8], shot);
            }
        }

        controlled_lt_msbs(&mut circ, &ctrl, &a, &b, k, &flag);

        let mut outs: Vec<QReg> = Vec::new();
        outs.extend(a);
        outs.extend(b);
        outs.push(ctrl);
        outs.push(flag);
        let (sim, detached) = circ.destroy_sim(outs);
        let flag_d = &detached[2 * n + 1];

        let mismatches = 0usize;
        for shot in 0..64 {
            let want_lt = (a_vals[shot] >> shift) < (b_vals[shot] >> shift);
            let want_eq = (a_vals[shot] >> shift) == (b_vals[shot] >> shift);
            let want_flag = ctrl_vals[shot] && want_lt;
            let got_flag = sim.read_bit_shot(flag_d, shot) == 1;
            if got_flag != want_flag {
                let _ = want_eq;
                panic!(
                    "shot {}: a={:#x} b={:#x} k={} ctrl={} want_lt={} got={}",
                    shot, a_vals[shot], b_vals[shot], k, ctrl_vals[shot], want_lt, got_flag
                );
            }
        }
        let _ = mismatches;
    }
}
