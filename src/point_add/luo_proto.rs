//! Classical prototype / qubit-budget scratchpad for Luo-style register sharing.
//!
//! Goal of this file: validate *cheaply* whether a Luo/PZ-style inversion track
//! is even compatible with the user's budget:
//!   - only ~600 qubits over the 512 input-point-coordinate qubits,
//!   - i.e. total target around 1100–1200 for the full point-add, or at least
//!     meaningfully below our current 2716q.
//!
//! We do NOT attempt a reversible circuit here. We only:
//!   1. model the qubit budget implied by Luo's register sharing,
//!   2. compare it to our current Kaliski scaffold,
//!   3. record the minimum structural consequences.
//!
//! Literature anchor: `/tmp/luo_ec_clean.txt`, especially Table 1 and Algorithm 3.

#![cfg(test)]

use alloy_primitives::U256;

use super::SECP256K1_P;

/// Very coarse qubit budget for the current live affine scaffold.
#[derive(Debug, Clone, Copy)]
struct Budget {
    tx_ty: usize,
    inversion_state: usize,
    lambda_and_mul_state: usize,
    classical_bits: usize,
    total: usize,
}

/// Current live build (best stable before the 511 detour):
/// - tx,ty = 2n = 512
/// - Kaliski persistent state = u,v_w,r,s,m_hist,f ≈ 4n + iters + 1
///   with iters ≈ 407/404 → ~1432
/// - live body state + mul transients explain the observed 2716 peak.
fn current_budget_estimate(n: usize, iters: usize) -> Budget {
    let tx_ty = 2 * n;
    let inversion_state = 4 * n + iters + 1;
    // Remaining gap to the observed peak (2716) is dominated by lam,
    // tmp_ext, carries, and a few flags.
    let total = 2716;
    let lambda_and_mul_state = total - tx_ty - inversion_state;
    Budget {
        tx_ty,
        inversion_state,
        lambda_and_mul_state,
        classical_bits: 2 * n, // ox, oy
        total,
    }
}

/// Luo-style inversion state from `/tmp/luo_ec_clean.txt`:
/// Table 1 says inversion can be done in roughly `3n + 4 log2 n + O(1)`
/// qubits *total* for the inversion component.
///
/// For n=256 this is about 3*256 + 4*8 = 800 qubits total, INCLUDING the
/// input/output pair of the inversion itself.
///
/// In our point-add context the inversion input is one n-bit value (dx or
/// similar), and we still need the 2n point coordinates live. So the key
/// number is the non-tx/ty overhead: about n + O(log n), not 4n+iters.
fn luo_inversion_total_qubits(n: usize) -> usize {
    3 * n + 4 * (n.ilog2() as usize)
}

/// Conservative point-add budget if we swapped our Kaliski block for a Luo/PZ
/// block *without* changing anything else in the affine scaffold.
fn naive_luo_point_add_budget_conservative(n: usize, current_other_peak: usize) -> usize {
    // Keep tx,ty live. Add the full Luo inversion block. Keep the rest of the
    // current non-inversion transients as-is.
    let tx_ty = 2 * n;
    let inversion_total = luo_inversion_total_qubits(n);
    tx_ty + inversion_total + current_other_peak
}

/// Optimistic overlap model for the same swap.
///
/// Luo's `3n + 4 log n` inversion count already includes the n-bit inversion
/// input register. In our point-add scaffold that register is already part of
/// tx/ty, so only `luo_total - n` is really "extra".
fn naive_luo_point_add_budget_optimistic(n: usize, current_other_peak: usize) -> usize {
    let tx_ty = 2 * n;
    let inversion_extra = luo_inversion_total_qubits(n) - n;
    tx_ty + inversion_extra + current_other_peak
}

/// Clean arithmetic helper for a tiny classical sanity check.
fn sub_mod(a: U256, b: U256, p: U256) -> U256 {
    if a >= b {
        (a - b) % p
    } else {
        p - ((b - a) % p)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn luo_budget_is_qubit_relevant() {
        let n = 256usize;
        let cur = current_budget_estimate(n, 407);
        eprintln!("current budget estimate: {:?}", cur);

        let current_other_peak = cur.lambda_and_mul_state;
        let luo_cons = naive_luo_point_add_budget_conservative(n, current_other_peak);
        let luo_opt = naive_luo_point_add_budget_optimistic(n, current_other_peak);
        eprintln!("naive Luo swap-in peak estimate: conservative={luo_cons}, optimistic={luo_opt}");

        assert!(luo_cons < cur.total, "Luo-style inversion must reduce peak (conservative)");
        assert!(luo_opt < cur.total, "Luo-style inversion must reduce peak (optimistic)");

        // User budget is ~600 qubits over the 512 input coords.
        assert!(luo_opt > 1112, "If this flips, Luo-alone got us into the user budget");
    }

    #[test]
    fn luo_alone_is_not_sota_but_is_structural() {
        let n = 256usize;
        let cur = current_budget_estimate(n, 407);
        let current_other_peak = cur.lambda_and_mul_state;
        let luo_cons = naive_luo_point_add_budget_conservative(n, current_other_peak);
        let luo_opt = naive_luo_point_add_budget_optimistic(n, current_other_peak);

        eprintln!(
            "luo_cons={}, luo_opt={}, current_total={}, saved_cons={}, saved_opt={}",
            luo_cons,
            luo_opt,
            cur.total,
            cur.total - luo_cons,
            cur.total - luo_opt
        );
        assert!(cur.total - luo_cons >= 500, "Luo should save ~500+ qubits conservatively");
        assert!(cur.total - luo_opt >= 800, "Luo should save ~800+ qubits optimistically");
    }

    #[test]
    fn even_free_inversion_needs_scaffold_collapse_to_hit_user_budget() {
        let n = 256usize;
        let cur = current_budget_estimate(n, 407);
        let target_total = 512 + 600;
        let inversion_free_total = cur.tx_ty + cur.lambda_and_mul_state;
        eprintln!(
            "inversion-free affine scaffold total={}, target_total={}, excess={}",
            inversion_free_total,
            target_total,
            inversion_free_total - target_total
        );
        assert_eq!(inversion_free_total, 1284);
        assert_eq!(inversion_free_total - target_total, 172);
        assert!(inversion_free_total > target_total);
    }

    #[test]
    fn optimistic_luo_still_needs_hundreds_more_qubits_cut() {
        let n = 256usize;
        let cur = current_budget_estimate(n, 407);
        let target_total = 512 + 600;
        let luo_opt = naive_luo_point_add_budget_optimistic(n, cur.lambda_and_mul_state);
        eprintln!(
            "optimistic Luo total={}, target_total={}, remaining_gap={}",
            luo_opt,
            target_total,
            luo_opt - target_total
        );
        assert_eq!(luo_opt, 1828);
        assert_eq!(luo_opt - target_total, 716);
    }

    #[test]
    fn dy_py_relation_sanity() {
        // Tiny guard against the kind of algebra drift we had earlier.
        let p = SECP256K1_P;
        let px = U256::from(123u64);
        let py = U256::from(456u64);
        let qx = U256::from(17u64);
        let qy = U256::from(31u64);
        let dx = sub_mod(px, qx, p);
        let dy = sub_mod(py, qy, p);
        assert_eq!(dx, U256::from(106u64));
        assert_eq!(dy, U256::from(425u64));
    }
}
