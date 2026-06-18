//! MBU (measurement-based uncomputation) primitives.
//!
//! Given a qubit q with val(q) = f(alive witnesses), free q cleanly:
//!
//!   1. HMR(q) → classical bit b. Kickback: `(-1)^(val(q)·b)`.
//!   2. Apply a phase gate sequence realising `(-1)^(f·b)`.
//!   3. Free q.
//!
//! Decompositions:
//!
//!   | f                   | phase gates                                 |
//!   |---------------------|---------------------------------------------|
//!   | 0 (constant)        | (nothing; HMR of |0> is trivial)            |
//!   | x                   | `z_if_bit(x`, b)                              |
//!   | x AND y             | `cz_if_bit(x`, y, b)                          |
//!   | x AND y AND z       | `ccz_if_bit(x`, y, z, b)                      |
//!   | x XOR y             | `z_if_bit(x`, b); `z_if_bit(y`, b)              |
//!   | x XOR y XOR z       | composed via one intermediate XOR qubit     |
//!   | x OR y = xy+x+y     | composed (xor of AND + xor of copies)       |
//!   | MAJ(x,y,z) = xy+yz+xz | composed (xor of three AND qubits)        |
//!
//! Why the (x XOR y) decomposition is z-of-x-THEN-z-of-y (not ccz):
//!
//!   (-1)^((x XOR y)·b) = (-1)^((x+y)·b) = (-1)^(xb) · (-1)^(yb)
//!
//! since in F2 `x XOR y = x+y`, and the product `(-1)^a · (-1)^b`
//! equals `(-1)^(a+b)` with `+` being F2 sum. Each `z_if_bit` applies
//! `(-1)^(q·b)` phase.
//!
//! The tracker must "see" q's identity to match the obligation.
//! For AND/XOR/COPY/AndOf3, the tracker's native transfer functions
//! track through. For OR/MAJ/3XOR we compose via intermediate ancillae
//! so q's `AbsVal` stays representable.

use crate::point_add::trailmix_port::arith::mcx::{mcx_clean_k, mcx_dirty_any_k, mcx_dirty_any_k_consume};
use crate::point_add::trailmix_port::circuit::{Circuit, QReg};

#[cfg(test)]
use crate::point_add::trailmix_port::arith::{cuccaro::*, mcx::*};

// =========================================================================
// Simple primitives — tracker tracks natively through forward ops.
// =========================================================================

// =========================================================================
// Composed primitives — allocate intermediates so tracker can follow.
// These leave the circuit with no net ancillae (all intermediates freed).
// =========================================================================

// and_tree_compute (Bennett-style O(log n)-anc AND reduction) was
// deleted: it's semantically a C^n X operation, and the unified
// primitive `mcx_dirty_any_k` (Theorem 3 recursion) serves the same
// purpose with 1 dirty ancilla instead of O(log n) clean ones.
// Callers pass a `dirty_bank` so the AND compute borrows its dirty
// from alive registers (paper's Theorem 4 pattern).

// =========================================================================
// The earlier Expr / compare_geq_const_witness path allocated one qubit
// per expression-tree node (O(depth × log run) peak ancillae, O(n) for
// adversarial constants), so it was removed in favor of
// compare_geq_theorem3 (polylog ancillae; ops are currently O(n^1.58),
// pending the V_2-based Theorem-3 construction).
// =========================================================================

#[cfg(test)]
mod tests {
    use super::{
        cinc_khattar_gidney, inc_khattar_gidney, kg_prefix_ancilla_count, mcx_clean_k,
        xor_and_of_khattar_gidney, KgPrefixAnd,
    };
    use crate::point_add::trailmix_port::circuit::{Circuit, QReg};

    // === Negative test: declare_and_of catches mismatched identity ===
    //
    // Compute q = x AND y (value = 1 when x=y=1).
    // Then call declare_and_of(q, x, y_wrong) where y_wrong = NOT y.
    // sim_mask check should fail and panic.
    #[test]
    #[should_panic(expected = "declare_and_of")]
    fn test_declare_and_of_catches_mismatch() {
        let mut circ = Circuit::new();
        let x = circ.alloc_qreg("x");
        let y = circ.alloc_qreg("y");
        circ.x(&x);
        circ.x(&y); // x=1, y=1
        let q = circ.alloc_qreg("q");
        circ.ccx(&x, &y, &q); // q = 1
                              // y_wrong = NOT y = 0. x AND y_wrong = 0 ≠ q = 1.
        let y_wrong = circ.alloc_qreg("y_wrong");
        // Leave y_wrong as |0> (= NOT 1 conceptually).
        circ.declare_and_of(&q, &x, &y_wrong); // should panic.
    }

    // =====================================================================
    // 1-bit adder / rolling compare tests — DELETED along with the
    // underlying primitives (misleading "2q peak" claim on compound
    // constant patterns; see git history for the Expr-extended version).
    // Theorem 3 (Vandaele 2026, Θ(n) gates + 1 dirty ancilla classical
    // comparator) is the real replacement.
    // =====================================================================

    // Expr witness evaluator tests deleted — primitive removed.

    // === Negative test: declare_and3_of catches mismatch ===
    #[test]
    #[should_panic(expected = "declare_and3_of")]
    fn test_declare_and3_of_catches_mismatch() {
        let mut circ = Circuit::new();
        let x = circ.alloc_qreg("x");
        let y = circ.alloc_qreg("y");
        let z = circ.alloc_qreg("z");
        circ.x(&x);
        circ.x(&y);
        circ.x(&z); // all 1
        let q = circ.alloc_qreg("q");
        circ.cx(&x, &q);
        circ.cx(&y, &q); // q = x XOR y = 0, not x AND y AND z = 1.
        circ.declare_and3_of(&q, &x, &y, &z); // should panic.
    }

    fn run_inc_khattar_gidney_case(n: usize, a_init: u64) {
        let mut circ = Circuit::new();
        let a = (0..n)
            .map(|i| circ.alloc_qreg(&format!("a{i}")))
            .collect::<Vec<_>>();
        // Use sim_load_reg_bytes_shot to set initial state without emitting X gates.
        // Direct X-gate setup would leave X(a[i]) immediately before inc_khattar_gidney's
        // own X(a[0]) for n=1 / a_init=1, triggering the redundant-op detector.
        {
            let mut bytes = vec![0u8; n.div_ceil(8)];
            for i in 0..n {
                if (a_init >> i) & 1 == 1 {
                    bytes[i / 8] |= 1u8 << (i % 8);
                }
            }
            circ.sim_load_reg_bytes_shot(&a, &bytes, 0);
        }
        inc_khattar_gidney(&mut circ, &a);
        let (sim, detached) = circ.destroy_sim(a);
        let got: u64 = (0..n)
            .map(|i| (sim.qubit_mask(&detached[i]) & 1) << i)
            .sum();
        let exp = (a_init + 1) & ((1u64 << n) - 1);
        assert_eq!(
            got,
            exp,
            "inc_khattar_gidney n={} a={:0w$b}",
            n,
            a_init,
            w = n
        );
        assert_eq!(sim.phase_mask(), 0, "inc_khattar_gidney phase n={}", n);
    }

    #[test]
    fn inc_khattar_gidney_n1_all() {
        for a in 0..(1u64 << 1) {
            run_inc_khattar_gidney_case(1, a);
        }
    }

    #[test]
    fn inc_khattar_gidney_n2_all() {
        for a in 0..(1u64 << 2) {
            run_inc_khattar_gidney_case(2, a);
        }
    }

    #[test]
    fn inc_khattar_gidney_n3_all() {
        for a in 0..(1u64 << 3) {
            run_inc_khattar_gidney_case(3, a);
        }
    }

    #[test]
    fn inc_khattar_gidney_n4_all() {
        for a in 0..(1u64 << 4) {
            run_inc_khattar_gidney_case(4, a);
        }
    }

    #[test]
    fn inc_khattar_gidney_n5_all() {
        for a in 0..(1u64 << 5) {
            run_inc_khattar_gidney_case(5, a);
        }
    }

    fn run_cinc_khattar_gidney_case(n: usize, ctrl_init: u64, a_init: u64) {
        let mut circ = Circuit::new();
        let ctrl = circ.alloc_qreg("ctrl");
        let a = (0..n)
            .map(|i| circ.alloc_qreg(&format!("a{i}")))
            .collect::<Vec<_>>();
        if ctrl_init == 1 {
            circ.x(&ctrl);
        }
        for i in 0..n {
            if (a_init >> i) & 1 == 1 {
                circ.x(&a[i]);
            }
        }
        cinc_khattar_gidney(&mut circ, &a, &ctrl);
        let mut outputs: Vec<crate::point_add::trailmix_port::circuit::QReg> = Vec::new();
        outputs.push(ctrl);
        outputs.extend(a);
        let (sim, detached) = circ.destroy_sim(outputs);
        let ctrl_d = &detached[0];
        let a_d = &detached[1..1 + n];
        let got_ctrl = sim.qubit_mask(ctrl_d) & 1;
        let got_a: u64 = (0..n).map(|i| (sim.qubit_mask(&a_d[i]) & 1) << i).sum();
        let exp_a = if ctrl_init == 1 {
            (a_init + 1) & ((1u64 << n) - 1)
        } else {
            a_init
        };
        assert_eq!(got_ctrl, ctrl_init, "cinc_khattar_gidney ctrl drift n={n}");
        assert_eq!(
            got_a,
            exp_a,
            "cinc_khattar_gidney n={} ctrl={} a={:0w$b}",
            n,
            ctrl_init,
            a_init,
            w = n,
        );
        assert_eq!(sim.phase_mask(), 0, "cinc_khattar_gidney phase n={}", n);
    }

    #[test]
    fn cinc_khattar_gidney_n4_all() {
        for ctrl in 0..2 {
            for a in 0..(1u64 << 4) {
                run_cinc_khattar_gidney_case(4, ctrl, a);
            }
        }
    }

    fn run_xor_and_of_khattar_gidney_case(n: usize, bits_init: u64, target_init: u64) {
        let mut circ = Circuit::new();
        let bits = (0..n)
            .map(|i| circ.alloc_qreg(&format!("b{i}")))
            .collect::<Vec<_>>();
        let target = circ.alloc_qreg("target");
        if target_init == 1 {
            circ.x(&target);
        }
        for i in 0..n {
            if (bits_init >> i) & 1 == 1 {
                circ.x(&bits[i]);
            }
        }
        xor_and_of_khattar_gidney(&mut circ, &bits, &target);
        let mut outputs: Vec<crate::point_add::trailmix_port::circuit::QReg> = Vec::new();
        outputs.push(target);
        outputs.extend(bits);
        let (sim, detached) = circ.destroy_sim(outputs);
        let target_d = &detached[0];
        let bits_d = &detached[1..1 + n];
        let got_target = sim.qubit_mask(target_d) & 1;
        let exp_and = if bits_init == (1u64 << n) - 1 { 1 } else { 0 };
        assert_eq!(
            got_target,
            target_init ^ exp_and,
            "xor_and_of_khattar_gidney target n={} bits={:0w$b} target={}",
            n,
            bits_init,
            target_init,
            w = n,
        );
        for i in 0..n {
            let got = sim.qubit_mask(&bits_d[i]) & 1;
            let exp = (bits_init >> i) & 1;
            assert_eq!(
                got, exp,
                "xor_and_of_khattar_gidney changed bit {} for n={}",
                i, n
            );
        }
        assert_eq!(
            sim.phase_mask(),
            0,
            "xor_and_of_khattar_gidney phase n={}",
            n
        );
    }

    #[test]
    fn xor_and_of_khattar_gidney_n6_all() {
        for n in 1usize..=6 {
            for bits in 0..(1u64 << n) {
                for target in 0..2 {
                    run_xor_and_of_khattar_gidney_case(n, bits, target);
                }
            }
        }
    }

    /// Verify `KgPrefixAnd::consume_with_body` yields the correct
    /// prefix-AND at every i ∈ [0, n], over random input patterns.
    /// Body for the test: `target[i] ^= AND(ctrls)` via `mcx_clean_k`
    /// (which handles the 0/1/2-ctrl small cases correctly).
    ///
    /// Includes a contract that captures bits_init pre-prefix and
    /// checks all target[i] post-consume against the classical
    /// prefix-AND.
    fn run_kg_prefix_and_streaming_case(n: usize, bits_init: u64) {
        let mut circ = Circuit::new();
        let q: Vec<QReg> = (0..n).map(|i| circ.alloc_qreg(&format!("q{i}"))).collect();
        for i in 0..n {
            if (bits_init >> i) & 1 == 1 {
                circ.x(&q[i]);
            }
        }
        let targets: Vec<QReg> = (0..=n).map(|i| circ.alloc_qreg(&format!("t{i}"))).collect();

        // CONTRACT [pre]: capture bits_init.
        {
            let q_refs_capture: Vec<&QReg> = q.iter().collect();
            circ.contract_capture(
                "kg_prefix_and_streaming",
                move |view, shot| -> Result<u64, String> {
                    let mut v = 0u64;
                    for (i, q) in q_refs_capture.iter().enumerate() {
                        if view.contract_read_bit_shot(q, shot) {
                            v |= 1u64 << i;
                        }
                    }
                    Ok(v)
                },
            );
        }

        // Separate target sets for forward + reverse so we can verify
        // each direction independently. targets_fwd is written by the
        // forward body, targets_rev by the reverse body.
        let targets_fwd: Vec<QReg> = targets;
        let targets_rev: Vec<QReg> = (0..=n)
            .map(|i| circ.alloc_qreg(&format!("rt{i}")))
            .collect();

        let anc_owned = circ.alloc_qreg_bits("kg_pa_anc", kg_prefix_ancilla_count(n));
        let anc_refs: Vec<&QReg> = anc_owned.iter().collect();
        let q_refs: Vec<&QReg> = q.iter().collect();

        let targets_fwd_refs: Vec<&QReg> = targets_fwd.iter().collect();
        let targets_rev_refs: Vec<&QReg> = targets_rev.iter().collect();
        KgPrefixAnd::new(&q_refs, &anc_refs)
            .forward(&mut circ, |c, i, ctrls| {
                let ctrl_owned: Vec<&QReg> = ctrls.to_vec();
                mcx_clean_k(c, &ctrl_owned, targets_fwd_refs[i]);
            })
            .reverse(&mut circ, |c, i, ctrls| {
                let ctrl_owned: Vec<&QReg> = ctrls.to_vec();
                mcx_clean_k(c, &ctrl_owned, targets_rev_refs[i]);
            });
        for q in anc_owned {
            circ.zero_and_free(q);
        }

        // CONTRACT [post]: both targets_fwd[i] and targets_rev[i] must
        // equal AND(q[0..i]) — the forward body and reverse body see
        // the SAME conditionally-clean ctrls per position.
        {
            let fwd_refs: Vec<&QReg> = targets_fwd.iter().collect();
            let rev_refs: Vec<&QReg> = targets_rev.iter().collect();
            let n_cap = n;
            circ.contract_pop_and_check::<u64, _>(
                "kg_prefix_and_streaming",
                move |captured, view, shot| -> Result<(), String> {
                    let bits = *captured;
                    for i in 0..=n_cap {
                        let mask_below = if i == 0 { 0 } else { (1u64 << i) - 1 };
                        let exp = if (bits & mask_below) == mask_below {
                            1u8
                        } else {
                            0
                        };
                        let got_fwd = if view.contract_read_bit_shot(fwd_refs[i], shot) {
                            1u8
                        } else {
                            0
                        };
                        let got_rev = if view.contract_read_bit_shot(rev_refs[i], shot) {
                            1u8
                        } else {
                            0
                        };
                        if got_fwd != exp {
                            return Err(format!(
                                "shot {}: targets_fwd[{}] = {} expected {} (bits={:#x}, n={})",
                                shot, i, got_fwd, exp, bits, n_cap,
                            ));
                        }
                        if got_rev != exp {
                            return Err(format!(
                                "shot {}: targets_rev[{}] = {} expected {} (bits={:#x}, n={})",
                                shot, i, got_rev, exp, bits, n_cap,
                            ));
                        }
                    }
                    Ok(())
                },
            );
        }

        let mut outs: Vec<QReg> = Vec::new();
        outs.extend(targets_fwd);
        outs.extend(targets_rev);
        outs.extend(q);
        let _ = circ.destroy_sim(outs);
    }

    #[test]
    fn kg_prefix_and_streaming_n6_all() {
        for n in 1usize..=6 {
            for bits in 0..(1u64 << n) {
                run_kg_prefix_and_streaming_case(n, bits);
            }
        }
    }

    #[test]
    fn xor_and_of_khattar_gidney_large_samples() {
        for &n in &[22usize, 33, 223] {
            let samples = [
                (vec![true; n], 0u64, 1u64),
                (
                    {
                        let mut v = vec![true; n];
                        v[0] = false;
                        v
                    },
                    0u64,
                    0u64,
                ),
                (
                    {
                        let mut v = vec![true; n];
                        v[n - 1] = false;
                        v
                    },
                    1u64,
                    1u64,
                ),
            ];
            for (bits_init, target_init, exp_target) in samples {
                let mut circ = Circuit::new();
                let bits = (0..n)
                    .map(|i| circ.alloc_qreg(&format!("b{i}")))
                    .collect::<Vec<_>>();
                let target = circ.alloc_qreg("target");
                if target_init == 1 {
                    circ.x(&target);
                }
                for i in 0..n {
                    if bits_init[i] {
                        circ.x(&bits[i]);
                    }
                }
                xor_and_of_khattar_gidney(&mut circ, &bits, &target);
                let mut outputs: Vec<crate::point_add::trailmix_port::circuit::QReg> = Vec::new();
                outputs.push(target);
                outputs.extend(bits);
                let (sim, detached) = circ.destroy_sim(outputs);
                let target_d = &detached[0];
                let bits_d = &detached[1..1 + n];
                let got_target = sim.qubit_mask(target_d) & 1;
                assert_eq!(
                    got_target, exp_target,
                    "xor_and_of_khattar_gidney large sample failed for n={n}"
                );
                for i in 0..n {
                    let got = sim.qubit_mask(&bits_d[i]) & 1;
                    let exp = bits_init[i] as u64;
                    assert_eq!(
                        got, exp,
                        "xor_and_of_khattar_gidney changed bit {} for n={}",
                        i, n
                    );
                }
                assert_eq!(
                    sim.phase_mask(),
                    0,
                    "xor_and_of_khattar_gidney phase n={}",
                    n
                );
            }
        }
    }
}

/// Gidney's n-bit incrementer (ZEROED ancilla variant), as drawn in
/// Vandaele 2026 Fig 8(b) (arXiv:2603.12917 ref [10]). n data bits + n-2
/// clean ancillae; ancillae end in |0⟩. 2(n-2) CCX + (n-1) CX + (2n-3) X.
///
/// Convention: a[0] = LSB, a[n-1] = MSB. Requires n ≥ 2.
///
/// Four slices (per paper page 21 caption):
///   Slice 1: forward CCX ladder anc[0] = a[0]·a[1], anc[k] = anc[k-1]·a[k+1]
///   Slice 2: CX(anc[k-1], a[k+1]) + X(a[k+1]) for k=1..n-2, plus CX(a[0],a[1])+X(a[1])
///            at the top and CX(anc[n-3], a[n-1]) at the bottom (no X on MSB)
///   Slice 3: reverse CCX ladder (bottom-up) zeroes the ancs using the
///            Bennett identity (anc[k] · `a_pre_slice2_relation` works out)
///   Slice 4: X on every data bit EXCEPT MSB a[n-1]
///
/// Verified exhaustively for n=6 (all 64 inputs, ancs zeroed) via Python
/// trace. Base case for Vandaele Theorem 4 (Lemma 7 substitutes the ancs
/// with a promise register).

/// ## Theorem 4 construction (Vandaele 2026, 1-dirty-ancilla INC)
///
/// Current:    `inc_gidney_fig8b` — n-2 CLEAN ancs, Θ(n) gates.  DONE.
/// Current:    `inc_gidney_fig8b_ctrl` — Lemma 7 (k=1), n-2 CLEAN
///             promise, Θ(n) gates.  DONE.
/// Next:       Lemma 8 — controlled strong-promise INC with 2⌈√n⌉
///             promise qubits. Construction via Fig 10 (paper p. 25):
///             split data into k=⌈√n⌉ blocks, 2 rounds alternating odd/
///             even block-INCs using Lemma 7, with X-flip shuffles on
///             promise markers between rounds. Each block-INC with ♢
///             promise expands via Eq. 43 triple (+1; +1; -1 on dirty).
/// Next:       Theorem 4 — INC with 1 dirty ancilla via Eq. 44 recursion:
///             let α = 2⌈√n⌉, β = n-α, ψ = dirty ancilla
///               1. X(α)
///               2. `Lemma8_promise_INC(β`, α as promise)
///               3. X(α)
///               4. fan-out CX(α → β, α → ψ) via Eq. 37
///               5. X(α)
///               6. `Lemma8_promise_DEC(β`, α as promise)
///               7. X(α)
///               8. fan-out CX(α → β, α → ψ) again
///               9. INC(α) [recurse; base case: direct n ≤ 3]
/// Final:      Corollary 7 wrapper — `conditional_increment(ctrl`, a, p) =
///             INC_{n-p+1}([ctrl, a[p..]]); X(ctrl). Uses Theorem 4 for
///             the inner INC, giving Θ(n) gates + 1 dirty ancilla.

/// Compute the full prefix-AND ladder for `bits` into `ladder`.
///
/// Semantics for `bits.len() >= 2`:
/// - `ladder[0] ^= bits[0] & bits[1]`
/// - `ladder[k] ^= ladder[k-1]_pre & bits[k+1]` for `k >= 1`
///
/// The ladder is a conditionally-clean workspace pattern: callers are
/// expected to run the exact inverse with [`prefix_and_ladder_rev_refs`] to
/// restore it. This is the persistent prefix substrate used by Gidney's
/// Fig. 8(b) incrementer and is also the right shape for a future
/// Khattar-Gidney Section 6.1 producer/consumer ladder.
pub(crate) fn prefix_and_ladder_fwd_refs(circ: &mut Circuit, bits: &[&QReg], ladder: &[&QReg]) {
    let n = bits.len();
    assert!(n >= 2, "prefix_and_ladder_fwd: n >= 2");
    assert_eq!(
        ladder.len(),
        n - 2,
        "prefix_and_ladder_fwd: expected {} ladder qubits, got {}",
        n - 2,
        ladder.len(),
    );
    if n == 2 {
        return;
    }
    circ.ccx(bits[0], bits[1], ladder[0]);
    for k in 1..(n - 2) {
        circ.ccx(ladder[k - 1], bits[k + 1], ladder[k]);
    }
}

/// Exact inverse of [`prefix_and_ladder_fwd_refs`].
pub(crate) fn prefix_and_ladder_rev_refs(circ: &mut Circuit, bits: &[&QReg], ladder: &[&QReg]) {
    let n = bits.len();
    assert!(n >= 2, "prefix_and_ladder_rev: n >= 2");
    assert_eq!(
        ladder.len(),
        n - 2,
        "prefix_and_ladder_rev: expected {} ladder qubits, got {}",
        n - 2,
        ladder.len(),
    );
    if n == 2 {
        return;
    }
    for k in (1..(n - 2)).rev() {
        circ.ccx(ladder[k - 1], bits[k + 1], ladder[k]);
    }
    circ.ccx(bits[0], bits[1], ladder[0]);
}

/// Reverse only the TOP `r` links of the prefix-AND ladder (highest-index
/// `r` links, top-down). Pairs with [`prefix_and_ladder_partial_fwd_refs`].
/// The untouched lower links stay live. Used by [`unary_iterate`] to update a
/// big-endian prefix-AND when only the low (= last-entering) bits change.
pub(crate) fn prefix_and_ladder_partial_rev_refs(
    circ: &mut Circuit,
    bits: &[&QReg],
    ladder: &[&QReg],
    r: usize,
) {
    let nlinks = ladder.len();
    let r = r.min(nlinks);
    for k in (nlinks - r..nlinks).rev() {
        if k == 0 {
            circ.ccx(bits[0], bits[1], ladder[0]);
        } else {
            circ.ccx(ladder[k - 1], bits[k + 1], ladder[k]);
        }
    }
}

/// Forward (recompute) only the TOP `r` links of the prefix-AND ladder,
/// bottom-up. Inverse of [`prefix_and_ladder_partial_rev_refs`].
pub(crate) fn prefix_and_ladder_partial_fwd_refs(
    circ: &mut Circuit,
    bits: &[&QReg],
    ladder: &[&QReg],
    r: usize,
) {
    let nlinks = ladder.len();
    let r = r.min(nlinks);
    for k in nlinks - r..nlinks {
        if k == 0 {
            circ.ccx(bits[0], bits[1], ladder[0]);
        } else {
            circ.ccx(ladder[k - 1], bits[k + 1], ladder[k]);
        }
    }
}

/// Unary iteration: for `i in 0..n_iters`, run `body(circ, i, gate)` where
/// `gate = (c == i)` is a freshly computed 1-qubit control, live during the
/// body and uncomputed right after. The n-bit little-endian counter `c` is
/// restored to its input value `v` on exit.
///
/// Cost: a single big-endian linear prefix-AND is built once; between steps
/// only the low bits of `c` change (gray-code update `c ^= i^(i+1)`), so the
/// prefix-AND is patched by partial reverse/forward of just the affected top
/// links. Amortized ~2 CCX/step for the patch + 2 CCX/step for the `gate`
/// detect, vs ~2n CCX/step for a full equality-MCX per step. See
/// `notes/live_intermediate_and_unary.md`.
///
/// PRECONDITION: `n_iters >= 1` and `n_iters <= 2^n`. For `v >= n_iters` the
/// gate never fires (v out of range); for `v < n_iters` it fires once, at `i = v`.
pub fn unary_iterate<F>(circ: &mut Circuit, c: &[&QReg], n_iters: usize, mut body: F)
where
    F: FnMut(&mut Circuit, usize, &QReg),
{
    let n = c.len();
    assert!(n >= 2, "unary_iterate: need n >= 2 counter bits");
    assert!(n_iters >= 1, "unary_iterate: n_iters >= 1");
    assert!(
        n >= 63 || n_iters <= (1usize << n),
        "unary_iterate: n_iters {n_iters} exceeds 2^{n}"
    );

    // Big-endian bit order: bits_be[j] = c[n-1-j]. The LSB c[0] is the separate
    // final control (= bits_be[n-1]); flipping it touches NO ladder link.
    let bits_be: Vec<&QReg> = c.iter().rev().copied().collect();

    // setup: c ^= all_ones  =>  c = ~v.
    for q in c {
        circ.x(q);
    }
    let ladder_owned = circ.alloc_qreg_bits("unary_ladder", n - 2);
    let ladder: Vec<&QReg> = ladder_owned.iter().collect();
    prefix_and_ladder_fwd_refs(circ, &bits_be, &ladder);

    // `top` = AND(c[1..n]); full all-ones = top & c[0].
    let top: &QReg = if ladder.is_empty() {
        bits_be[0]
    } else {
        ladder[ladder.len() - 1]
    };
    let c0 = c[0];
    let gate = circ.alloc_qreg("unary_gate");

    for i in 0..n_iters {
        // gate = (c == all_ones) = (~v ^ i == ~0) = (v == i).
        circ.ccx(top, c0, &gate);
        body(circ, i, &gate);
        circ.ccx(top, c0, &gate); // uncompute gate

        if i + 1 < n_iters {
            let m = i ^ (i + 1); // low-contiguous run of b bits (= 2^b - 1)
            let b = m.count_ones() as usize;
            // c[0] is the separate final control; c[1..b-1] live in the top
            // (b-1) ladder links. Patch them.
            let r = b.saturating_sub(1);
            prefix_and_ladder_partial_rev_refs(circ, &bits_be, &ladder, r);
            for j in 0..b.min(n) {
                circ.x(c[j]); // c[0..b-1] ^= m
            }
            prefix_and_ladder_partial_fwd_refs(circ, &bits_be, &ladder, r);
        }
    }

    circ.zero_and_free(gate);
    prefix_and_ladder_rev_refs(circ, &bits_be, &ladder); // uncompute ladder
    drop(ladder);
    for q in ladder_owned {
        circ.zero_and_free(q);
    }

    // restore c = v: current c = ~v ^ (n_iters-1); XOR ~(n_iters-1).
    let last = n_iters - 1;
    for (j, q) in c.iter().enumerate() {
        if (last >> j) & 1 == 0 {
            circ.x(q);
        }
    }
}

/// Lowest layer index whose ops reference ANY of `qubits` (by pointer).
/// `usize::MAX` if none touch it. Used to bound the gray-code partial rewind:
/// reverse layers `[k..]` (which undoes, top-down, every CCX that consumed the
/// changed bits + the conditionally-clean ancillae built on them) before
/// flipping, then re-run `[k..]` forward.
fn lowest_layer_touching(layers: &[KgPrefixLayer], qubits: &[&QReg]) -> usize {
    let mut k = usize::MAX;
    for (i, layer) in layers.iter().enumerate() {
        let hit = layer.ops.iter().any(|op| match op {
            KgPrefixOp::X(q) => qubits.iter().any(|cq| std::ptr::eq(*cq, *q)),
            KgPrefixOp::Ccx(a, b, t) => qubits
                .iter()
                .any(|cq| std::ptr::eq(*cq, *a) || std::ptr::eq(*cq, *b) || std::ptr::eq(*cq, *t)),
        });
        if hit {
            k = k.min(i);
            break; // layers only touch a bit at/after its entry; first hit is the lowest
        }
    }
    k
}

/// Same contract as [`unary_iterate`] but on the Khattar-Gidney LOG\* prefix-AND
/// (`kg_prefix_ancilla_count(n-1)` ≈ log\*(n) ancillae) instead of the linear
/// `prefix_and_ladder` (n-2). The all-ones detector `(c == i)` is `AND(top-prefix)
/// AND c[0]`; between steps the gray-code `i^(i+1)` (a contiguous LSB run) is
/// applied by PARTIAL-rewinding only the KG layer suffix that touches the changed
/// bits: reverse `[k..]`, flip, forward `[k..]`. `n-1` counter bits feed the
/// prefix-AND (the top), `c[0]` is the separate final control.
pub fn unary_iterate_log_star<F>(circ: &mut Circuit, c: &[&QReg], n_iters: usize, mut body: F)
where
    F: FnMut(&mut Circuit, usize, &QReg),
{
    let n = c.len();
    assert!(n >= 2, "unary_iterate_log_star: need n >= 2 counter bits");
    assert!(n_iters >= 1, "unary_iterate_log_star: n_iters >= 1");
    assert!(
        n >= 63 || n_iters <= (1usize << n),
        "unary_iterate_log_star: n_iters {n_iters} exceeds 2^{n}"
    );

    // c ^= all_ones => c = ~v, so (c == all_ones) == (v == i).
    for q in c {
        circ.x(q);
    }
    let bits_be: Vec<&QReg> = c.iter().rev().copied().collect(); // [c[n-1] .. c[0]]
    let c0 = c[0]; // = bits_be[n-1], separate final control
    let nb = n - 1; // prefix-AND over the top nb bits = AND(c[1..n])
    let gate = circ.alloc_qreg("uls_gate");

    if nb == 1 {
        // top = bits_be[0] = c[n-1]; gate = top AND c0.
        let top = bits_be[0];
        for i in 0..n_iters {
            circ.ccx(top, c0, &gate);
            body(circ, i, &gate);
            circ.ccx(top, c0, &gate);
            if i + 1 < n_iters {
                let b = (i ^ (i + 1)).count_ones() as usize;
                for j in 0..b.min(n) {
                    circ.x(c[j]);
                }
            }
        }
        circ.zero_and_free(gate);
    } else {
        let pa_bits: Vec<&QReg> = bits_be[0..nb].to_vec();
        let anc_owned = circ.alloc_qreg_bits("uls_anc", kg_prefix_ancilla_count(nb));
        let anc: Vec<&QReg> = anc_owned.iter().collect();
        let layers = kg_get_layers_for_prefix_and(&pa_bits, &anc);
        // forward: compute the prefix-AND (held in the conditionally-clean ancs).
        for layer in &layers {
            for &op in &layer.ops {
                op.emit(circ);
            }
        }
        // all-ones detector = AND(layers[nb].ctrls) (the full prefix at position nb).
        let base: Vec<&QReg> = layers[nb].ctrls.clone();
        for i in 0..n_iters {
            // gate = AND(base) AND c0 = (c == all_ones) = (v == i).
            let mut gc: Vec<&QReg> = base.clone();
            gc.push(c0);
            mcx_clean_k(circ, &gc, &gate);
            body(circ, i, &gate);
            mcx_clean_k(circ, &gc, &gate); // uncompute
            if i + 1 < n_iters {
                let b = (i ^ (i + 1)).count_ones() as usize;
                // changed counter bits c[0..b-1]: c[0]=c0 (separate), c[1..b-1] = bits_be[n-1-j].
                let changed_pa: Vec<&QReg> = (1..b.min(n)).map(|j| bits_be[n - 1 - j]).collect();
                let k = lowest_layer_touching(&layers, &changed_pa);
                if k != usize::MAX {
                    for layer in layers[k..].iter().rev() {
                        for &op in layer.ops.iter().rev() {
                            op.emit(circ);
                        }
                    }
                }
                for j in 0..b.min(n) {
                    circ.x(c[j]); // flip changed bits (c0 + the prefix bits)
                }
                if k != usize::MAX {
                    for layer in &layers[k..] {
                        for &op in &layer.ops {
                            op.emit(circ);
                        }
                    }
                }
            }
        }
        circ.zero_and_free(gate);
        // reverse all layers (uncompute the prefix-AND).
        for layer in layers.iter().rev() {
            for &op in layer.ops.iter().rev() {
                op.emit(circ);
            }
        }
        for q in anc_owned {
            circ.zero_and_free(q);
        }
    }

    // restore c = v: current c = ~v ^ (n_iters-1); XOR ~(n_iters-1).
    let last = n_iters - 1;
    for (j, q) in c.iter().enumerate() {
        if (last >> j) & 1 == 0 {
            circ.x(q);
        }
    }
}

#[derive(Clone, Copy, Debug)]
enum KgPrefixOp<'a> {
    X(&'a QReg),
    Ccx(&'a QReg, &'a QReg, &'a QReg),
}

impl KgPrefixOp<'_> {
    #[inline]
    fn emit(self, circ: &mut Circuit) {
        match self {
            KgPrefixOp::X(q) => circ.x(q),
            KgPrefixOp::Ccx(a, b, t) => circ.ccx(a, b, t),
        }
    }
}

#[derive(Clone, Debug)]
struct KgPrefixLayer<'a> {
    ctrls: Vec<&'a QReg>,
    ops: Vec<KgPrefixOp<'a>>,
}

fn kg_get_layer_id(x: usize) -> usize {
    let mut layer_id = 0usize;
    let mut s = 0usize;
    while s <= x {
        s += (1usize << layer_id) + 1;
        layer_id += 1;
    }
    layer_id - 1
}

fn kg_start_layer(layer_id: usize) -> usize {
    let mut s = 0usize;
    for i in 0..layer_id {
        s += (1usize << i) + 1;
    }
    s
}

/// Upper-bound ancilla budget for the Khattar-Gidney prefix layer
/// decomposition. Some `n` cause the recursive layer builder to
/// reference fewer ancillae than this bound — call
/// [`kg_prefix_ancilla_count_exact`] for the precise count.
#[must_use]
pub fn kg_prefix_ancilla_count(n: usize) -> usize {
    if n <= 1 {
        return 0;
    }
    let targets_len = kg_get_layer_id(n - 1) + 1;
    if targets_len <= 2 {
        1
    } else {
        2 + kg_prefix_ancilla_count(targets_len)
    }
}

fn kg_apply_prefix_controlled_x(circ: &mut Circuit, ctrls: &[&QReg], target: &QReg) {
    match ctrls {
        [] => circ.x(target),
        [c] => circ.cx(c, target),
        [a, b] => circ.ccx(a, b, target),
        _ => panic!(
            "kg_apply_prefix_controlled_x: expected <=2 ctrls, got {}",
            ctrls.len()
        ),
    }
}

fn kg_anc_index(len: usize, idx: isize) -> usize {
    if idx >= 0 {
        idx as usize
    } else {
        (len as isize + idx) as usize
    }
}

fn kg_get_layers_for_prefix_and<'a>(
    q: &[&'a QReg],
    inp_anc: &[&'a QReg],
) -> Vec<KgPrefixLayer<'a>> {
    assert!(
        !q.is_empty(),
        "kg_get_layers_for_prefix_and: q must be non-empty"
    );
    if q.len() == 1 {
        return vec![
            KgPrefixLayer {
                ctrls: Vec::new(),
                ops: Vec::new(),
            },
            KgPrefixLayer {
                ctrls: vec![q[0]],
                ops: Vec::new(),
            },
        ];
    }
    assert!(
        inp_anc.len() >= kg_prefix_ancilla_count(q.len()),
        "kg_get_layers_for_prefix_and: expected at least {} ancillae for n={}, got {}",
        kg_prefix_ancilla_count(q.len()),
        q.len(),
        inp_anc.len(),
    );

    let n = q.len();
    let n_layers = kg_get_layer_id(q.len() - 1);
    let mut ret = vec![KgPrefixLayer {
        ctrls: Vec::new(),
        ops: Vec::new(),
    }];
    let mut targets: Vec<&'a QReg> = Vec::new();
    let mut anc: Vec<&'a QReg> = vec![inp_anc[0]];

    for layer_id in 0..=n_layers {
        let st = kg_start_layer(layer_id);
        let en = n.min(kg_start_layer(layer_id + 1));

        let mut layer_ctrls = targets.clone();
        layer_ctrls.push(q[st]);
        ret.push(KgPrefixLayer {
            ctrls: layer_ctrls,
            ops: Vec::new(),
        });

        for i in (st + 1)..en {
            let offset = i - st;
            let anc_len = anc.len();
            let q0 = q[i];
            let (q1, t) = if offset == 1 {
                (q[i - 1], anc[kg_anc_index(anc_len, -1)])
            } else {
                (
                    anc[kg_anc_index(anc_len, -(offset as isize - 1))],
                    anc[kg_anc_index(anc_len, -(offset as isize))],
                )
            };
            let mut ops = Vec::new();
            if std::ptr::eq(t, inp_anc[0]) {
                ops.push(KgPrefixOp::Ccx(q0, q1, t));
            } else {
                ops.push(KgPrefixOp::X(t));
                ops.push(KgPrefixOp::Ccx(q0, q1, t));
            }
            let mut ctrls = targets.clone();
            ctrls.push(t);
            ret.push(KgPrefixLayer { ctrls, ops });
        }

        let layer_len = en - st;
        let push_idx = kg_anc_index(anc.len(), 1 - layer_len as isize);
        targets.push(anc[push_idx]);

        let slice_start = kg_anc_index(anc.len(), 2 - layer_len as isize);
        let mut next_anc = anc[slice_start..].to_vec();
        next_anc.extend(q[st..en].iter());
        anc = next_anc;
    }

    if targets.len() <= 2 {
        return ret;
    }

    ret.push(KgPrefixLayer {
        ctrls: Vec::new(),
        ops: Vec::new(),
    });
    let target_prefix_layers = kg_get_layers_for_prefix_and(&targets, &inp_anc[2..]);
    for layer_id in 1..=n_layers {
        let st = kg_start_layer(layer_id);
        let en = n.min(kg_start_layer(layer_id + 1));
        let target_prefix_targets = target_prefix_layers[layer_id].ctrls.clone();
        ret[st + 1]
            .ops
            .extend_from_slice(&target_prefix_layers[layer_id].ops);

        let temp_target = if target_prefix_targets.len() == 1 {
            target_prefix_targets[0]
        } else {
            assert_eq!(target_prefix_targets.len(), 2);
            ret[st + 1].ops.push(KgPrefixOp::Ccx(
                target_prefix_targets[0],
                target_prefix_targets[1],
                inp_anc[1],
            ));
            inp_anc[1]
        };

        for i in st..en {
            let local = *ret[i + 1]
                .ctrls
                .last()
                .expect("kg_get_layers_for_prefix_and: empty local ctrl");
            ret[i + 1].ctrls = vec![temp_target, local];
        }

        if target_prefix_targets.len() == 2 {
            ret[en + 1].ops.push(KgPrefixOp::Ccx(
                target_prefix_targets[0],
                target_prefix_targets[1],
                temp_target,
            ));
        }
    }

    ret
}

/// Streaming Khattar-Gidney prefix-AND (Sec 4 / Fig 4 of KG 2025).
///
/// Builds the prefix-AND ladder with `log*(n)` clean ancillae and
/// exposes a per-position control-set callback so callers can run
/// arbitrary bodies (e.g. a strided-XOR demux for bitlen-via-popcount)
/// against `AND(q[0..i])` for each `i` without paying the full
/// w-ancilla prefix-OR scratch of a naive thermometer.
///
/// USAGE (matches the conditionally-clean construction — body is
/// invoked layer-by-layer in DESCENDING order, interleaved with the
/// per-layer reverse-ops as in `inc_khattar_gidney_refs_inner`):
/// ```text
/// let anc_owned = circ.alloc_qreg_bits("kg_pa", kg_prefix_ancilla_count(n));
/// let anc_refs: Vec<&QReg> = anc_owned.iter().collect();
/// let q_refs: Vec<&QReg> = q.iter().collect();
/// // new() emits the forward sweep; consume_with_body emits the
/// // reverse sweep, calling `body(i, ctrls)` once per layer i in
/// // descending order (i = n .. 0). The body sees the layer's
/// // ctrls in their conditionally-clean state, where the AND of the
/// // ctrls equals AND(q[0..i]) at the moment of the call.
/// KgPrefixAnd::new(circ, &q_refs, &anc_refs)
///     .consume_with_body(circ, |c, i, ctrls| {
///         // body example: for each k such that 2^k | i,
///         // emit `mcx_clean_k(ctrls, clz[k])`.
///         for k in 0..clz.len() {
///             if i > 0 && (i & ((1 << k) - 1)) == 0 {
///                 mcx_clean_k(c, ctrls, &clz[k]);
///             }
///         }
///     });
/// for q in anc_owned { circ.zero_and_free(q); }
/// ```
///
/// INVARIANTS for the body:
/// - MUST treat ctrls as read-only.
/// - MUST be reversible across the call (XOR-style writes into
///   caller-owned output bits are fine).
///
/// Cost: forward+reverse sweep = ~2(2n-3) Toffoli + linear X's, plus
/// whatever the body emits per layer.
/// Ancillae: `kg_prefix_ancilla_count(n)` ≈ `log*(n)`.
///
/// NOTE on the conditionally-clean trick — during the forward sweep,
/// some input q-bits get temporarily X-bracketed (used as "borrowed"
/// ancillae). The ctrls' raw bits do NOT equal the prefix-AND in
/// isolation; the AND of the layer's ctrls equals the prefix-AND
/// only at the specific reverse-iteration point for THAT layer.
/// This is why the API interleaves body+reverse-op per layer rather
/// than exposing a static `ctrls_at(i)` lookup.
/// Phase-1 of the streaming prefix-AND. `KgPrefixAnd::new()` returns
/// this; the caller must then `.forward(circ, body)` to emit the
/// forward sweep, which yields a [`KgPrefixAndForwardDone`] that the
/// caller `.reverse(circ, body)`s. Rust's type system enforces the
/// ordering at compile time — you cannot call reverse before forward,
/// and you cannot skip forward.
pub struct KgPrefixAnd<'a> {
    layers: Vec<KgPrefixLayer<'a>>,
    /// = `q.len()` at construction. `layers.len()` may exceed `n+1`
    /// because the recursion appends sync-placeholder layers.
    n: usize,
}

/// Phase-2 of the streaming prefix-AND, after the forward sweep has
/// been emitted. The only thing you can do with this is `.reverse(...)`.
pub struct KgPrefixAndForwardDone<'a> {
    layers: Vec<KgPrefixLayer<'a>>,
    n: usize,
}

impl<'a> KgPrefixAnd<'a> {
    /// Allocate the prefix-AND plan. Emits NO quantum ops; just builds
    /// the layer schedule. Call `.forward(circ, body)` to actually
    /// emit the forward sweep.
    ///
    /// `q`: the input bits (length n; assumed in some pure state).
    /// `anc_refs`: at least `kg_prefix_ancilla_count(q.len())` qubits,
    ///             each in |0>. The caller owns the underlying `QRegs`.
    #[track_caller]
    #[must_use]
    pub fn new(q: &[&'a QReg], anc_refs: &[&'a QReg]) -> Self {
        assert!(!q.is_empty(), "KgPrefixAnd::new: q must be non-empty");
        let needed = kg_prefix_ancilla_count(q.len());
        assert!(
            anc_refs.len() >= needed,
            "KgPrefixAnd::new: needed {} ancillae for n={}, got {}",
            needed,
            q.len(),
            anc_refs.len()
        );
        let n = q.len();
        let layers = kg_get_layers_for_prefix_and(q, anc_refs);
        Self { layers, n }
    }

    /// Number of input bits `n` (= `q.len()`).
    #[must_use]
    pub fn n(&self) -> usize {
        self.n
    }

    /// Emit the forward sweep with an ASCENDING body. For each
    /// position layer i ∈ [0, n], emits layer i's forward ops and
    /// THEN calls `body(circ, i, &layer.ctrls)`. After all forward
    /// ops are emitted, returns [`KgPrefixAndForwardDone`] which the
    /// caller can `.reverse(...)`.
    ///
    /// At the moment of the body call, the AND of `ctrls` equals
    /// the prefix-AND `AND(q[0..i])` (the conditionally-clean
    /// identity holds both immediately after layer i's forward ops
    /// AND at the corresponding reverse-iter moment).
    ///
    /// Layers at i > n are recursion-sync placeholders — their forward
    /// ops are emitted but body is skipped.
    ///
    /// Pass `|_, _, _| {}` as the body if you only need the reverse
    /// pass to do work.
    pub fn forward(
        self,
        circ: &mut Circuit,
        mut body: impl FnMut(&mut Circuit, usize, &[&'a QReg]),
    ) -> KgPrefixAndForwardDone<'a> {
        for (i, layer) in self.layers.iter().enumerate() {
            for &op in &layer.ops {
                op.emit(circ);
            }
            if i <= self.n {
                body(circ, i, &layer.ctrls);
            }
        }
        KgPrefixAndForwardDone {
            layers: self.layers,
            n: self.n,
        }
    }
}

impl<'a> KgPrefixAndForwardDone<'a> {
    /// Number of input bits `n` (= `q.len()`).
    #[must_use]
    pub fn n(&self) -> usize {
        self.n
    }

    /// Emit the reverse sweep with a DESCENDING body. For each
    /// position layer i ∈ [n, 0], calls `body(circ, i, &layer.ctrls)`
    /// FIRST and then emits layer i's reverse ops. Consumes self;
    /// after return, all ancillae are restored to |0> (caller still
    /// owns the `QRegs` and must `zero_and_free` them).
    ///
    /// Pass `|_, _, _| {}` as the body if you only needed the forward
    /// pass to do work.
    pub fn reverse(
        self,
        circ: &mut Circuit,
        mut body: impl FnMut(&mut Circuit, usize, &[&'a QReg]),
    ) {
        for (i, layer) in self.layers.iter().enumerate().rev() {
            if i <= self.n {
                body(circ, i, &layer.ctrls);
            }
            for &op in layer.ops.iter().rev() {
                op.emit(circ);
            }
        }
    }
}

/// Khattar-Gidney 2025 incrementer: recursively produce/consume the
/// prefix-AND ladder with `log*_2(n)` clean ancillae.
///
/// This is a direct port of the Zenodo Qualtran reference artifact's
/// `get_layers_for_prefix_and` + incrementer wrapper into this repo's
/// gate set. The internal decomposition ensures every target flip is
/// controlled by at most 2 qubits; the recursion lives in the prefix
/// layer producer, not in the final increment consumer.
pub fn inc_khattar_gidney(circ: &mut Circuit, a: &[QReg]) {
    let a_refs: Vec<&QReg> = a.iter().collect();
    inc_khattar_gidney_refs(circ, &a_refs);
}

/// Reference-slice variant of [`inc_khattar_gidney`] for callers that
/// have a `Vec<QReg>` (e.g. when prepending a ctrl qubit). Avoids the
/// owned-slice constraint of the public API.
pub fn inc_khattar_gidney_refs(circ: &mut Circuit, a: &[&QReg]) {
    inc_khattar_gidney_refs_inner(circ, a, /*skip_lsb_x=*/ false);
}

/// Khattar-Gidney increment with an optional skip of the i=0
/// reverse-layer X(a[0]). `cinc_khattar_gidney` uses skip=true to
/// fold its trailing X(ctrl) into this routine — the two X's would
/// otherwise be a redundant pair.
fn inc_khattar_gidney_refs_inner(circ: &mut Circuit, a: &[&QReg], skip_lsb_x: bool) {
    let n = a.len();
    if n == 0 {
        return;
    }
    if n == 1 {
        if !skip_lsb_x {
            circ.x(a[0]);
        }
        return;
    }

    // Use the over-bound for safety: kg_prefix_ancilla_count_exact's
    // dry-run only counts ancs the layers fn writes/reads in ITS dry-run
    // pass, which uses the over-bound's ancs as input. The actual layer
    // builder, when given fewer ancs, may underflow on the recursion's
    // inp_anc[2..] slice. Over-bound + zero_and_free per-anc satisfies
    // strict-dealloc (each free emits an R touch).
    let anc_owned = circ.alloc_qreg_bits("kg_inc_anc", kg_prefix_ancilla_count(n - 1));
    let anc_refs: Vec<&QReg> = anc_owned.iter().collect();
    let a_top: &[&QReg] = &a[..n - 1];
    let layers = kg_get_layers_for_prefix_and(a_top, &anc_refs);

    for layer in &layers {
        for &op in &layer.ops {
            op.emit(circ);
        }
    }
    for (i, layer) in layers.iter().enumerate().rev() {
        if i < n && !(i == 0 && skip_lsb_x) {
            kg_apply_prefix_controlled_x(circ, &layer.ctrls, a[i]);
        }
        // Emit inverse ops.
        for &op in layer.ops.iter().rev() {
            op.emit(circ);
        }
    }
    // Free anc qubits via zero_and_free. The kg_prefix_ancilla_count_exact
    // upper bound (max_used+1) overcounts when the index range is sparse;
    // some ancillae are allocated but never touched by any layer op. Calling
    // zero_and_free emits an R gate per ancilla, satisfying the strict
    // "must be touched before free" check (R is the canonical end-of-life
    // marker, not a dummy gate). All ancillae are |0> here because the
    // forward+reverse layer pass is reversible.
    drop(layers);
    drop(anc_refs);
    for q in anc_owned {
        circ.zero_and_free(q);
    }
}

/// `target ^= AND(bits)` using the Khattar-Gidney prefix decomposition.
///
/// This reuses the same `log*_2(n)`-clean prefix producer as
/// [`inc_khattar_gidney`], but consumes only the full-prefix control.
pub fn xor_and_of_khattar_gidney(circ: &mut Circuit, bits: &[QReg], target: &QReg) {
    let bits_refs: Vec<&QReg> = bits.iter().collect();
    xor_and_of_khattar_gidney_refs(circ, &bits_refs, target);
}

/// Same operation as [`xor_and_of_khattar_gidney_refs`], but uses caller-owned
/// clean ancillae instead of allocating new lanes. The caller must provide at
/// least [`kg_prefix_ancilla_count(bits.len())`] qubits, all initialized to
/// |0>; they are restored to |0> on return.
pub fn xor_and_of_khattar_gidney_refs_with_anc<'a>(
    circ: &mut Circuit,
    bits: &[&'a QReg],
    target: &QReg,
    anc_refs: &[&'a QReg],
) {
    match bits.len() {
        0 => {
            circ.x(target);
            return;
        }
        1 => {
            circ.cx(bits[0], target);
            return;
        }
        2 => {
            circ.ccx(bits[0], bits[1], target);
            return;
        }
        _ => {}
    }
    assert!(
        anc_refs.len() >= kg_prefix_ancilla_count(bits.len()),
        "xor_and_of_khattar_gidney_refs_with_anc: need {} clean ancillae for n={}, got {}",
        kg_prefix_ancilla_count(bits.len()),
        bits.len(),
        anc_refs.len(),
    );

    // PRE: capture (AND(bits)_pre, target_pre).
    {
        let bits_for_capture: Vec<&QReg> = bits.to_vec();
        let target_ref = target;
        circ.contract_capture(
            "mbu.xor_and_kg_refs_with_anc.pre",
            move |view, shot| -> Result<(bool, bool), String> {
                let mut and_v = true;
                for q in &bits_for_capture {
                    and_v &= view.contract_read_bit_shot(q, shot);
                }
                let t = view.contract_read_bit_shot(target_ref, shot);
                Ok((and_v, t))
            },
        );
    }

    let layers = kg_get_layers_for_prefix_and(bits, anc_refs);

    for (i, layer) in layers.iter().enumerate() {
        if i > bits.len() {
            break;
        }
        for &op in &layer.ops {
            op.emit(circ);
        }
    }

    for (i, layer) in layers.iter().enumerate().rev() {
        if i > bits.len() {
            continue;
        }
        if i == bits.len() {
            kg_apply_prefix_controlled_x(circ, &layer.ctrls, target);
        }
        for &op in layer.ops.iter().rev() {
            op.emit(circ);
        }
    }
    drop(layers);

    // POST: target ^= AND(bits); bits unchanged. The provided ancillae are
    // restored structurally by the prefix ladder; their eventual zero/free is
    // the caller's responsibility.
    {
        let bits_for_check: Vec<&QReg> = bits.to_vec();
        let target_ref = target;
        circ.contract_pop_and_check::<(bool, bool), _>(
            "mbu.xor_and_kg_refs_with_anc.pre",
            move |cap, view, shot| -> Result<(), String> {
                let (and_pre, t_pre) = *cap;
                let mut and_post = true;
                for q in &bits_for_check {
                    and_post &= view.contract_read_bit_shot(q, shot);
                }
                if and_post != and_pre {
                    return Err(format!(
                        "xor_and_kg_with_anc: bits AND changed {} -> {}",
                        u8::from(and_pre),
                        u8::from(and_post)
                    ));
                }
                let t_post = view.contract_read_bit_shot(target_ref, shot);
                let expected = t_pre ^ and_pre;
                if t_post != expected {
                    return Err(format!(
                        "xor_and_kg_with_anc: target {}->{} expected {} (t_pre={}, AND={})",
                        u8::from(t_pre),
                        u8::from(t_post),
                        u8::from(expected),
                        u8::from(t_pre),
                        u8::from(and_pre),
                    ));
                }
                Ok(())
            },
        );
    }
}

/// Variant of [`xor_and_of_khattar_gidney_refs`] that ALSO frees
/// `target` at its last gate-touch (the prefix-controlled-X). The
/// ancilla cleanup pass that follows does not touch `target`, so the
/// strict-dealloc gap is zero.
pub fn xor_and_of_khattar_gidney_refs_consume(circ: &mut Circuit, bits: &[&QReg], target: QReg) {
    match bits.len() {
        0 => {
            circ.x(&target);
            drop(target);
            return;
        }
        1 => {
            circ.cx(bits[0], &target);
            drop(target);
            return;
        }
        2 => {
            circ.ccx(bits[0], bits[1], &target);
            drop(target);
            return;
        }
        _ => {}
    }
    let anc_owned = circ.alloc_qreg_bits("kg_and_anc", kg_prefix_ancilla_count(bits.len()));
    let anc_refs: Vec<&QReg> = anc_owned.iter().collect();
    let layers = kg_get_layers_for_prefix_and(bits, &anc_refs);

    for (i, layer) in layers.iter().enumerate() {
        if i > bits.len() {
            break;
        }
        for &op in &layer.ops {
            op.emit(circ);
        }
    }

    let mut target_slot = Some(target);
    for (i, layer) in layers.iter().enumerate().rev() {
        if i > bits.len() {
            continue;
        }
        if i == bits.len() {
            let t = target_slot.as_ref().expect("target only consumed once");
            kg_apply_prefix_controlled_x(circ, &layer.ctrls, t);
            // Last gate-touch on target. Free it now so the strict-
            // dealloc check sees gap=0.
            drop(target_slot.take());
        }
        for &op in layer.ops.iter().rev() {
            op.emit(circ);
        }
    }
    drop(layers);
    drop(anc_refs);
    for q in anc_owned {
        circ.zero_and_free(q);
    }
}

/// Reference-slice variant of [`xor_and_of_khattar_gidney`].
///
/// `target ^= AND(bits[0..])` via the Khattar–Gidney Sec 5.3 / Sec 6.1
/// prefix-AND ladder (Fig 4 in the paper). 2n-3 Toffolis, log*_2(n)
/// clean ancillae, O(log n) depth.
pub fn xor_and_of_khattar_gidney_refs(circ: &mut Circuit, bits: &[&QReg], target: &QReg) {
    // PRE: capture (AND(bits)_pre, target_pre).
    {
        let bits_for_capture: Vec<&QReg> = bits.to_vec();
        let target_ref = target;
        circ.contract_capture(
            "mbu.xor_and_kg_refs.pre",
            move |view, shot| -> Result<(bool, bool), String> {
                let mut and_v = true;
                for q in &bits_for_capture {
                    and_v &= view.contract_read_bit_shot(q, shot);
                }
                let t = view.contract_read_bit_shot(target_ref, shot);
                Ok((and_v, t))
            },
        );
    }

    xor_and_of_khattar_gidney_refs_inner(circ, bits, target);

    // POST: target ^= AND(bits); bits unchanged.
    {
        let bits_for_check: Vec<&QReg> = bits.to_vec();
        let target_ref = target;
        circ.contract_pop_and_check::<(bool, bool), _>(
            "mbu.xor_and_kg_refs.pre",
            move |cap, view, shot| -> Result<(), String> {
                let (and_pre, t_pre) = *cap;
                let mut and_post = true;
                for q in &bits_for_check {
                    and_post &= view.contract_read_bit_shot(q, shot);
                }
                if and_post != and_pre {
                    return Err(format!(
                        "xor_and_kg: bits AND changed {} -> {}",
                        u8::from(and_pre),
                        u8::from(and_post)
                    ));
                }
                let t_post = view.contract_read_bit_shot(target_ref, shot);
                let expected = t_pre ^ and_pre;
                if t_post != expected {
                    return Err(format!(
                        "xor_and_kg: target {}->{} expected {} (t_pre={}, AND={})",
                        u8::from(t_pre),
                        u8::from(t_post),
                        u8::from(expected),
                        u8::from(t_pre),
                        u8::from(and_pre),
                    ));
                }
                Ok(())
            },
        );
    }
}

fn xor_and_of_khattar_gidney_refs_inner(circ: &mut Circuit, bits: &[&QReg], target: &QReg) {
    match bits.len() {
        0 => {
            circ.x(target);
            return;
        }
        1 => {
            circ.cx(bits[0], target);
            return;
        }
        2 => {
            circ.ccx(bits[0], bits[1], target);
            return;
        }
        _ => {}
    }

    // Use the over-bound here (matches what kg_get_layers_for_prefix_and
    // asserts). The "_exact" count is sometimes lower than the layer
    // builder's actual recursion needs (specifically when an outer call's
    // inp_anc[2..] needs to satisfy an inner call's _exact requirement
    // — the outer _exact can under-count what the inner needs). Until
    // _exact is rewritten to recurse via kg_prefix_ancilla_count_exact
    // on the inner targets, use the over-bound to avoid panics.
    let anc_owned = circ.alloc_qreg_bits("kg_and_anc", kg_prefix_ancilla_count(bits.len()));
    let anc_refs: Vec<&QReg> = anc_owned.iter().collect();
    let layers = kg_get_layers_for_prefix_and(bits, &anc_refs);

    // The target XOR fires exactly at layer index bits.len(). Layers with
    // index > bits.len() hold no target injection and their computed ancillae
    // serve as controls only for those higher-index layers — never for the
    // actual target XOR. They are dead computations (their ops cancel exactly
    // in the forward+reverse pair) and must be omitted to avoid the
    // redundant-op detector firing on the seam.
    for (i, layer) in layers.iter().enumerate() {
        if i > bits.len() {
            break;
        }
        for &op in &layer.ops {
            op.emit(circ);
        }
    }

    for (i, layer) in layers.iter().enumerate().rev() {
        if i > bits.len() {
            continue;
        }
        if i == bits.len() {
            kg_apply_prefix_controlled_x(circ, &layer.ctrls, target);
        }
        // Emit inverse ops.
        for &op in layer.ops.iter().rev() {
            op.emit(circ);
        }
    }
    // Free via zero_and_free — see inc_khattar_gidney_refs note about
    // sparse ancilla indices needing R to satisfy strict-dealloc.
    drop(layers);
    drop(anc_refs);
    for q in anc_owned {
        circ.zero_and_free(q);
    }
}

/// Controlled increment via the standard `[ctrl] ++ a` wrapper:
/// `a += ctrl (mod 2^n)`.
pub fn cinc_khattar_gidney(circ: &mut Circuit, a: &[QReg], ctrl: &QReg) {
    let a_refs: Vec<&QReg> = a.iter().collect();
    cinc_khattar_gidney_refs(circ, &a_refs, ctrl);
}

/// Reference-slice variant of [`cinc_khattar_gidney`].
pub fn cinc_khattar_gidney_refs(circ: &mut Circuit, a: &[&QReg], ctrl: &QReg) {
    if a.is_empty() {
        return;
    }

    // PRE: capture (a_pre, ctrl_pre).
    let n = a.len();
    {
        let a_for_capture: Vec<&QReg> = a.to_vec();
        let ctrl_ref = ctrl;
        circ.contract_capture(
            "mbu.cinc_kg_refs.pre",
            move |view, shot| -> Result<(u128, bool), String> {
                let cap = if n >= 128 { 128 } else { n };
                let mut av: u128 = 0;
                for b in 0..cap {
                    if view.contract_read_bit_shot(a_for_capture[b], shot) {
                        av |= 1u128 << b;
                    }
                }
                let cv = view.contract_read_bit_shot(ctrl_ref, shot);
                Ok((av, cv))
            },
        );
    }

    // cinc(a, ctrl) = inc(combined=[ctrl, a]) followed by X(ctrl) to
    // undo the LSB flip — but inc_khattar_gidney's i=0 reverse-layer
    // op IS that X, so the published "inc-then-X" pair cancels. Use
    // the skip-LSB-X variant to emit the optimized sequence directly.
    let mut combined: Vec<&QReg> = Vec::with_capacity(1 + a.len());
    combined.push(ctrl);
    combined.extend(a.iter().copied());
    inc_khattar_gidney_refs_inner(circ, &combined, /*skip_lsb_x=*/ true);

    // POST: a == (a_pre + ctrl_pre) mod 2^n; ctrl unchanged.
    {
        let a_for_check: Vec<&QReg> = a.to_vec();
        let ctrl_ref = ctrl;
        circ.contract_pop_and_check::<(u128, bool), _>(
            "mbu.cinc_kg_refs.pre",
            move |cap, view, shot| -> Result<(), String> {
                let (a_pre, c_pre) = *cap;
                let cap_n = if n >= 128 { 128 } else { n };
                let mut a_post: u128 = 0;
                for b in 0..cap_n {
                    if view.contract_read_bit_shot(a_for_check[b], shot) {
                        a_post |= 1u128 << b;
                    }
                }
                let mask = if cap_n >= 128 {
                    !0u128
                } else {
                    (1u128 << cap_n) - 1
                };
                let expected = (a_pre.wrapping_add(u128::from(c_pre))) & mask;
                if a_post != expected {
                    return Err(format!(
                        "cinc_kg: a {:#x}->{:#x}, expected {:#x} (a_pre={:#x}, ctrl={})",
                        a_pre,
                        a_post,
                        expected,
                        a_pre,
                        u8::from(c_pre),
                    ));
                }
                let c_post = view.contract_read_bit_shot(ctrl_ref, shot);
                if c_post != c_pre {
                    return Err(format!(
                        "cinc_kg: ctrl changed {} -> {}",
                        u8::from(c_pre),
                        u8::from(c_post)
                    ));
                }
                Ok(())
            },
        );
    }
}

// [DELETED 2026-05-30] `controlled_add_cuccaro` (10n CCX via cccx
// shared-scratch ancilla) and its post-check have been removed. All
// callers route through [`controlled_add_cuccaro_3n`] (3n CCX). See
// git log for the deleted body.

// [DELETED 2026-05-30] `controlled_add_cuccaro_mbu` and
// `controlled_add_cuccaro_mbu_refs` (8n CCX via streaming MBU-AND on
// each cccx) have been removed. All callers route through
// [`controlled_add_cuccaro_3n`] / [`controlled_add_cuccaro_3n_refs`]
// (3n CCX). See git log for the deleted body.

/// Classical-quantum compare. CURRENTLY IMPLEMENTED INCORRECTLY
/// per HARD RULE — ops are O(n^1.58) (should be Θ(n)); polylog ancs
/// are honest. Bridge to Theorem 3 via `V_2` stack is under construction
/// (`compare_lt_qq` done, `compare_geq_cq` via temp register is REJECTED
/// because it's n transient ancs).
///
/// Do not call on performance-critical paths until `V_2` Theorem 3 lands.
///
/// Recursive halving:
///   x >= c  iff  (`x_hi` > `c_hi`) OR (`x_hi` == `c_hi` AND `x_lo` >= `c_lo`)
///
/// Sub-calls:
///   1. compute s1 = (`x_lo` >= `c_lo`)     [recurse on lo]
///   2. out ^= (`x_hi` > `c_hi`)            [recurse on hi with c+1]
///   3. out ^= s1 AND (`x_hi` == `c_hi`)    [AND-tree eq + conjunction]
///   4. uncompute s1                    [recurse self-inverse]
///
/// Base cases: n ≤ 2 direct.
///
/// Ops: Θ(n log n) (master theorem). Ancs: O(log² n) (one scratch
/// per recursion level + log-depth AND trees).
pub fn compare_geq_theorem3(circ: &mut Circuit, x: &[QReg], c: &[u8], out: &QReg) {
    let n = x.len();
    if n == 0 {
        // Empty x vs empty (or larger) c: vacuously x=0 >= c=0, so
        // out ^= 1 iff c is numerically 0.
        if bytes_ge_pow2(c, 0) {
            // c >= 2^0 = 1, so c > 0, x < c, out unchanged.
        } else {
            circ.x(out);
        }
        return;
    }
    // If c >= 2^n, x < c always, out unchanged.
    if bytes_ge_pow2(c, n) {
        return;
    }
    // If c == 0, x >= 0 always, out ^= 1.
    if bytes_is_zero(c) {
        circ.x(out);
        return;
    }
    if n == 1 {
        // c > 0 and c < 2^1 = 2, so c == 1.
        // x >= 1 iff x == 1. out ^= x[0].
        circ.cx(&x[0], out);
        return;
    }
    if n == 2 {
        // c in {1, 2, 3}.
        let c0 = bit_of(c, 0);
        let c1 = bit_of(c, 1);
        match (c1, c0) {
            (false, true) => {
                // c=1: x >= 1 iff x != 0. out ^= (x[0] OR x[1]).
                //   x[0] OR x[1] = NOT(NOT x[0] AND NOT x[1]).
                //   Simpler: out ^= x[0]; out ^= x[1]; out ^= x[0]·x[1].
                circ.cx(&x[0], out);
                circ.cx(&x[1], out);
                circ.ccx(&x[0], &x[1], out);
            }
            (true, false) => {
                // c=2: x >= 2 iff x[1] = 1. out ^= x[1].
                circ.cx(&x[1], out);
            }
            (true, true) => {
                // c=3: x >= 3 iff x = 3, i.e. x[0]·x[1].
                circ.ccx(&x[0], &x[1], out);
            }
            _ => unreachable!(),
        }
        return;
    }

    // General n >= 3: dispatch to compare_lt_cq_paper (Vandaele 2026
    // Theorem 3 / Fig 7 / Eq 32 with Fig 2(a) dirty upgrade).
    //
    // compare_lt_cq_paper gives z ^= 1[x < c] with O(n log n) gates and
    // 1 dirty ancilla. We want out ^= 1[x >= c] = out ^= 1 ^ 1[x < c],
    // so we X(out) to pick up the constant-1 contribution and then call
    // compare_lt_cq_paper to XOR in 1[x < c].
    //
    // c passed to compare_lt_cq_paper must be exactly n bits, so
    // construct a Vec<u8> from the low n bits of c.
    let c_bits: Vec<u8> = (0..n).map(|i| u8::from(bit_of(c, i))).collect();
    circ.x(out);
    compare_lt_cq_paper(circ, x, &c_bits, out);
}

/// Reference-slice variant of [`compare_geq_theorem3`].
pub fn compare_geq_theorem3_refs(circ: &mut Circuit, x: &[&QReg], c: &[u8], out: &QReg) {
    let n = x.len();
    if n == 0 {
        if bytes_ge_pow2(c, 0) {
            // c >= 1, x = 0 < c, out unchanged.
        } else {
            circ.x(out);
        }
        return;
    }
    if bytes_ge_pow2(c, n) {
        return;
    }
    if bytes_is_zero(c) {
        circ.x(out);
        return;
    }
    if n == 1 {
        circ.cx(x[0], out);
        return;
    }
    if n == 2 {
        let c0 = bit_of(c, 0);
        let c1 = bit_of(c, 1);
        match (c1, c0) {
            (false, true) => {
                circ.cx(x[0], out);
                circ.cx(x[1], out);
                circ.ccx(x[0], x[1], out);
            }
            (true, false) => {
                circ.cx(x[1], out);
            }
            (true, true) => {
                circ.ccx(x[0], x[1], out);
            }
            _ => unreachable!(),
        }
        return;
    }

    let c_bits: Vec<u8> = (0..n).map(|i| u8::from(bit_of(c, i))).collect();
    circ.x(out);
    compare_lt_cq_paper_refs(circ, x, &c_bits, out);
}

/// Reference-slice variant of [`compare_geq_theorem3_free_out`].
pub fn compare_geq_theorem3_free_out_refs(circ: &mut Circuit, x: &[&QReg], c: &[u8], out: QReg) {
    let n = x.len();
    if n == 0 {
        if bytes_ge_pow2(c, 0) {
            drop(out);
        } else {
            circ.x(&out);
            drop(out);
        }
        return;
    }
    if bytes_ge_pow2(c, n) {
        return;
    }
    if bytes_is_zero(c) {
        circ.x(&out);
        drop(out);
        return;
    }
    if n == 1 {
        circ.cx(x[0], &out);
        drop(out);
        return;
    }
    if n == 2 {
        let c0 = bit_of(c, 0);
        let c1 = bit_of(c, 1);
        match (c1, c0) {
            (false, true) => {
                circ.cx(x[0], &out);
                circ.cx(x[1], &out);
                circ.ccx(x[0], x[1], &out);
                drop(out);
            }
            (true, false) => {
                circ.cx(x[1], &out);
                drop(out);
            }
            (true, true) => {
                circ.ccx(x[0], x[1], &out);
                drop(out);
            }
            _ => unreachable!(),
        }
        return;
    }
    let c_bits: Vec<u8> = (0..n).map(|i| u8::from(bit_of(c, i))).collect();
    circ.x(&out);
    compare_lt_cq_paper_free_z_refs(circ, x, &c_bits, out);
}

/// Extract bit `i` from a little-endian byte vec.
fn bit_of(bytes: &[u8], i: usize) -> bool {
    let byte_idx = i / 8;
    if byte_idx >= bytes.len() {
        return false;
    }
    (bytes[byte_idx] >> (i % 8)) & 1 == 1
}

fn bytes_is_zero(bytes: &[u8]) -> bool {
    bytes.iter().all(|&b| b == 0)
}

/// True iff the numeric value of `bytes` is >= 2^k.
fn bytes_ge_pow2(bytes: &[u8], k: usize) -> bool {
    // Any bit at position >= k being 1 ⇒ value >= 2^k.
    for i in 0..bytes.len() * 8 {
        if i >= k && (bytes[i / 8] >> (i % 8)) & 1 == 1 {
            return true;
        }
    }
    false
}

pub fn cinc_gidney_halving(circ: &mut Circuit, a: &[QReg], ctrl: &QReg) {
    let n = a.len();
    if n == 0 {
        return;
    }
    if n == 1 {
        circ.cx(ctrl, &a[0]);
        return;
    }
    if n == 2 {
        circ.ccx(ctrl, &a[0], &a[1]);
        circ.cx(ctrl, &a[0]);
        return;
    }

    let m = n / 2;

    // g := AND(a[0..m], ctrl). Borrow dirty from a[m] (high half).
    let g = circ.alloc_qreg("ghalv_g");
    let mut ctrls: Vec<&QReg> = Vec::with_capacity(m + 1);
    ctrls.extend(a[..m].iter());
    ctrls.push(ctrl);
    let dirty = &a[m];
    mcx_dirty_any_k(circ, &ctrls, &g, dirty);

    // Propagate carry into high half.
    cinc_gidney_halving(circ, &a[m..], &g);

    // Uncompute g (self-inverse: a[0..m] and ctrl are unchanged above).
    mcx_dirty_any_k(circ, &ctrls, &g, dirty);

    drop(g);

    // Low-half increment.
    cinc_gidney_halving(circ, &a[..m], ctrl);
}

/// Theorem 5 (Vandaele 2026): classical-quantum adder `x += c mod 2^n`.
/// Θ(n log n) gates; uses the caller-supplied dirty ancilla `g`, which
/// is returned to its original (unknown) value on exit.
///
/// Recursive Häner-et-al. structure: split x,c at m=⌈n/2⌉. The high
/// half xH takes the carry from xL+cL via a controlled INC, which is
/// surrounded by two self-inverse CARRY compares against the same
/// threshold — xL is untouched between them so the second call cleanly
/// zeroes g before the recursive sub-adds run.
///
///   CARRY(xL ≥ 2^m − cL → g)   [g ^= carry]
///   cinc(xH, ctrl=g, p=0)       [xH += g]
///   CARRY(xL ≥ 2^m − cL → g)   [g restored]
///   `add_classical(xL`, cL, g)    [recurse on low half]
///   `add_classical(xH`, cH, g)    [recurse on high half]
///
/// CARRY: `compare_geq_theorem3` (polylog ancs; ops O(n^1.58) pending
/// V_2-based rewrite to Θ(n)). cinc: `cinc_gidney_halving`.
///
/// Correctness note: xL+cL ≥ 2^m iff xL ≥ 2^m − cL, so the forward
/// comparator gives the carry. Because xL is not modified between
/// the two CARRY calls (cinc only touches xH, and the first recursive
/// call on xL happens AFTER the uncompute), the second CARRY XORs the
/// same value back into g.
pub fn classical_quantum_add(circ: &mut Circuit, x: &[QReg], c: &[u8], g: &QReg) {
    let n = x.len();
    if n == 0 {
        return;
    }
    // All-zero c: adding 0 is a no-op. Prevents recursion from emitting
    // wasted cinc_gidney_halving / recursive calls through all-zero subtrees, which
    // is the bulk of cost when c is sparse (e.g. c = R = 2^32 + 977 on
    // a 256-bit x: top-half c is all zero).
    if c.iter().all(|&b| b == 0) {
        return;
    }
    if n == 1 {
        if (c[0] & 1) == 1 {
            circ.x(&x[0]);
        }
        return;
    }
    let m = n.div_ceil(2);

    // Split c into low m bits (c_lo) and high n-m bits (c_hi).
    let c_lo = extract_low_bits(c, m);
    let c_hi = extract_bit_range(c, m, n);
    let x_lo = &x[..m];
    let x_hi = &x[m..];

    let c_lo_zero = c_lo.iter().all(|&b| b == 0);
    let c_hi_zero = c_hi.iter().all(|&b| b == 0);

    // Threshold T = 2^m - c_lo for the CARRY compare x_lo >= T.
    // Byte-level subtraction (u128 overflows at m=129).
    let (t_is_zero, t_bytes) = two_pow_m_minus(&c_lo, m);

    // If c_lo=0: no carry from low half, skip the CARRY compare AND
    // the controlled INC (which would fire with g=0 = no-op, but still
    // emits gates). Only the high-half recursion has work to do.
    if !c_lo_zero {
        // Forward CARRY: g ^= 1[xL >= T]. T=0 ↔ c_lo=0 ↔ carry impossible.
        if !t_is_zero {
            compare_geq_theorem3(circ, x_lo, &t_bytes, g);
        }

        // Controlled INC: x_hi += g.
        cinc_gidney_halving(circ, x_hi, g);

        // Reverse CARRY: self-inverse — restores g.
        if !t_is_zero {
            compare_geq_theorem3(circ, x_lo, &t_bytes, g);
        }
    }

    // Recurse on halves, skipping all-zero subtrees.
    if !c_lo_zero {
        classical_quantum_add(circ, x_lo, &c_lo, g);
    }
    if !c_hi_zero {
        classical_quantum_add(circ, x_hi, &c_hi, g);
    }
}

/// Compute `(2^m − val)` over m-bit unsigned integers, returning
/// `(is_zero, bytes_le)`. When val == 0 the result is 2^m which is
/// representationally zero in m-bit arithmetic — flag it so callers
/// can skip the compare (no m-bit x can satisfy x ≥ 2^m).
/// Byte-level subtraction so it handles m > 128 without u128 overflow.
fn two_pow_m_minus(val: &[u8], m: usize) -> (bool, Vec<u8>) {
    if m == 0 {
        return (true, vec![]);
    }
    // Check zero-ness cheaply.
    let is_zero = val.iter().all(|&b| b == 0);
    if is_zero {
        return (true, vec![0u8; m.div_ceil(8)]);
    }
    // Compute t = 2^m − val. Equivalent to (~val (over m bits)) + 1.
    let mut t = vec![0u8; m.div_ceil(8)];
    for i in 0..m {
        let byte_idx = i / 8;
        let bit = if byte_idx < val.len() {
            (val[byte_idx] >> (i % 8)) & 1
        } else {
            0
        };
        if bit == 0 {
            t[byte_idx] |= 1u8 << (i % 8);
        }
    }
    // Now t = ~val (over m bits). Add 1.
    let mut carry: u16 = 1;
    for byte in &mut t {
        let sum = u16::from(*byte) + carry;
        *byte = (sum & 0xff) as u8;
        carry = sum >> 8;
        if carry == 0 {
            break;
        }
    }
    // Mask off bits beyond m in the top byte.
    let top_bits = m % 8;
    if top_bits != 0 {
        let mask = (1u8 << top_bits) - 1;
        let top = t.len() - 1;
        t[top] &= mask;
    }
    (false, t)
}

/// Corollary 8 (Vandaele 2026): 1-controlled classical-quantum adder.
/// If ctrl=1, a += val; else a unchanged. Θ(n log n) gates, 1 dirty
/// ancilla (internally allocated).
///
/// Implementation: a ctrl-gated Theorem 5. Every gate emitted by the
/// adder has ctrl added to its control list — CX becomes CCX, X becomes
/// CX, CCX becomes C³X (expanded via a scratch AND(ctrl, other)). The
/// recursive structure, CARRY compares, and controlled INC all inherit
/// the outer ctrl.
///
/// Here we achieve that by computing `g_eff = ctrl AND g_internal`
/// once per recursion level, so the controlled INC sees the combined
/// (ctrl AND carry) and the CARRY compares are guarded by ctrl via
/// an extra CCX layer on the final XOR.
/// Returns true iff `ctrl_cq_add_impl` will ever execute a CARRY step
/// (i.e., will ever touch the g ancilla).
fn ctrl_cq_add_uses_g(n: usize, c: &[u8]) -> bool {
    if n <= 1 {
        return false;
    }
    if c.iter().all(|&b| b == 0) {
        return false;
    }
    let m = n.div_ceil(2);
    let c_lo = extract_low_bits(c, m);
    let c_hi = extract_bit_range(c, m, n);
    let c_lo_zero = c_lo.iter().all(|&b| b == 0);
    if !c_lo_zero {
        // CARRY step fires here, touching g.
        return true;
    }
    // c_lo is zero; g only used if recursive x_hi call uses it.
    let c_hi_zero = c_hi.iter().all(|&b| b == 0);
    if c_hi_zero {
        return false;
    }
    ctrl_cq_add_uses_g(n - m, &c_hi)
}

pub fn controlled_classical_quantum_add(circ: &mut Circuit, ctrl: &QReg, a: &[QReg], val: &[u8]) {
    let a_refs: Vec<&QReg> = a.iter().collect();
    controlled_classical_quantum_add_refs(circ, ctrl, &a_refs, val);
}

/// Reference-slice variant of [`controlled_classical_quantum_add`].
pub fn controlled_classical_quantum_add_refs(
    circ: &mut Circuit,
    ctrl: &QReg,
    a: &[&QReg],
    val: &[u8],
) {
    if !ctrl_cq_add_uses_g(a.len(), val) {
        ctrl_cq_add_impl_refs(circ, ctrl, a, val, ctrl);
        return;
    }
    let g = circ.alloc_qreg("cadd_dirty_g");
    ctrl_cq_add_impl_consume_g_refs(circ, ctrl, a, val, g);
}

fn ctrl_cq_add_impl_refs(circ: &mut Circuit, ctrl: &QReg, x: &[&QReg], c: &[u8], g: &QReg) {
    let n = x.len();
    if n == 0 {
        return;
    }
    if c.iter().all(|&b| b == 0) {
        return;
    }
    if n == 1 {
        if (c[0] & 1) == 1 {
            circ.cx(ctrl, x[0]);
        }
        return;
    }
    let m = n.div_ceil(2);
    let c_lo = extract_low_bits(c, m);
    let c_hi = extract_bit_range(c, m, n);
    let x_lo = &x[..m];
    let x_hi = &x[m..];

    let c_lo_zero = c_lo.iter().all(|&b| b == 0);
    let c_hi_zero = c_hi.iter().all(|&b| b == 0);

    let (t_is_zero, t_bytes) = two_pow_m_minus(&c_lo, m);

    if !c_lo_zero {
        if !t_is_zero {
            let s = circ.alloc_qreg("cqadd_cmp_s");
            compare_geq_theorem3_refs(circ, x_lo, &t_bytes, &s);
            circ.ccx(ctrl, &s, g);
            compare_geq_theorem3_free_out_refs(circ, x_lo, &t_bytes, s);
        }

        let cg = circ.alloc_qreg("cqadd_cg");
        circ.ccx(ctrl, g, &cg);
        cinc_khattar_gidney_refs(circ, x_hi, &cg);
        circ.ccx(ctrl, g, &cg);
        drop(cg);

        if !t_is_zero {
            let s = circ.alloc_qreg("cqadd_cmp_s2");
            compare_geq_theorem3_refs(circ, x_lo, &t_bytes, &s);
            circ.ccx(ctrl, &s, g);
            compare_geq_theorem3_free_out_refs(circ, x_lo, &t_bytes, s);
        }
    }

    if !c_lo_zero {
        ctrl_cq_add_impl_refs(circ, ctrl, x_lo, &c_lo, g);
    }
    if !c_hi_zero {
        ctrl_cq_add_impl_refs(circ, ctrl, x_hi, &c_hi, g);
    }
}

/// Like `ctrl_cq_add_impl_refs` but frees `g` immediately after its last
/// gate-touch, deep inside the recursion. Avoids the strict-dealloc
/// gap that occurs when the caller frees `g` after trailing CX ops
/// that don't touch `g`.
///
/// Invariant: `ctrl_cq_add_uses_g(x.len(), c)` must be true.
fn ctrl_cq_add_impl_consume_g_refs(
    circ: &mut Circuit,
    ctrl: &QReg,
    x: &[&QReg],
    c: &[u8],
    g: QReg,
) {
    let n = x.len();
    debug_assert!(
        n >= 2 && ctrl_cq_add_uses_g(n, c),
        "ctrl_cq_add_impl_consume_g_refs: g is not used (n={n})"
    );

    let m = n.div_ceil(2);
    let c_lo = extract_low_bits(c, m);
    let c_hi = extract_bit_range(c, m, n);
    let x_lo = &x[..m];
    let x_hi = &x[m..];
    let c_lo_zero = c_lo.iter().all(|&b| b == 0);
    let c_hi_zero = c_hi.iter().all(|&b| b == 0);
    let (t_is_zero, t_bytes) = two_pow_m_minus(&c_lo, m);

    // Determine which sub-recursion is the LAST to touch g.
    let hi_uses_g = !c_hi_zero && ctrl_cq_add_uses_g(x_hi.len(), &c_hi);
    let lo_uses_g = !c_lo_zero && ctrl_cq_add_uses_g(x_lo.len(), &c_lo);
    // If neither sub-recursion uses g, the CARRY section at this level is
    // the last user. We restructure the CARRY to free g at the exact
    // last gate-touch, before the compare uncomputation trailing ops.
    let carry_is_last_g_user = !hi_uses_g && !lo_uses_g;

    let mut g_holder = Some(g);

    if !c_lo_zero {
        if !t_is_zero {
            let s = circ.alloc_qreg("cqadd_cmp_s");
            compare_geq_theorem3_refs(circ, x_lo, &t_bytes, &s);
            circ.ccx(ctrl, &s, g_holder.as_ref().expect("g alive"));
            compare_geq_theorem3_free_out_refs(circ, x_lo, &t_bytes, s);
        }

        let cg = circ.alloc_qreg("cqadd_cg");
        circ.ccx(ctrl, g_holder.as_ref().expect("g alive"), &cg);
        cinc_khattar_gidney_refs(circ, x_hi, &cg);
        circ.ccx(ctrl, g_holder.as_ref().expect("g alive"), &cg);
        if carry_is_last_g_user && t_is_zero {
            let _ = g_holder.take();
        }
        drop(cg);

        if !t_is_zero {
            let s = circ.alloc_qreg("cqadd_cmp_s2");
            compare_geq_theorem3_refs(circ, x_lo, &t_bytes, &s);
            circ.ccx(ctrl, &s, g_holder.as_ref().expect("g alive"));
            if carry_is_last_g_user {
                let _ = g_holder.take();
            }
            compare_geq_theorem3_free_out_refs(circ, x_lo, &t_bytes, s);
        }
    }

    if hi_uses_g {
        let g = g_holder.take().expect("g alive for hi_uses_g recursion");
        if !c_lo_zero {
            ctrl_cq_add_impl_refs(circ, ctrl, x_lo, &c_lo, &g);
        }
        ctrl_cq_add_impl_consume_g_refs(circ, ctrl, x_hi, &c_hi, g);
    } else if lo_uses_g {
        let g = g_holder.take().expect("g alive for lo_uses_g recursion");
        ctrl_cq_add_impl_consume_g_refs(circ, ctrl, x_lo, &c_lo, g);
        debug_assert!(
            c_hi_zero || !ctrl_cq_add_uses_g(x_hi.len(), &c_hi),
            "ctrl_cq_add_impl_consume_g_refs: lo_uses_g branch needs hi-side g"
        );
        if !c_hi_zero {
            ctrl_cq_add_impl_refs(circ, ctrl, x_hi, &c_hi, ctrl);
        }
    } else {
        if !c_lo_zero {
            ctrl_cq_add_impl_refs(circ, ctrl, x_lo, &c_lo, ctrl);
        }
        if !c_hi_zero {
            ctrl_cq_add_impl_refs(circ, ctrl, x_hi, &c_hi, ctrl);
        }
    }
}

/// Extract the low `bits` bits of `src` (byte-packed, LSB first) into
/// a byte vector sized to hold `bits` bits.
fn extract_low_bits(src: &[u8], bits: usize) -> Vec<u8> {
    let n_bytes = bits.div_ceil(8);
    let mut out = vec![0u8; n_bytes];
    for i in 0..bits {
        let byte_idx = i / 8;
        if byte_idx < src.len() && (src[byte_idx] >> (i % 8)) & 1 == 1 {
            out[i / 8] |= 1 << (i % 8);
        }
    }
    out
}

/// Extract bits [lo..hi) of `src` into a byte vector aligned to the new LSB.
/// Bits beyond `src.len()`*8 are treated as 0 (zero-extension).
fn extract_bit_range(src: &[u8], lo: usize, hi: usize) -> Vec<u8> {
    let bits = hi - lo;
    let n_bytes = bits.div_ceil(8);
    let mut out = vec![0u8; n_bytes];
    for i in 0..bits {
        let sidx = lo + i;
        let byte_idx = sidx / 8;
        if byte_idx < src.len() && (src[byte_idx] >> (sidx % 8)) & 1 == 1 {
            out[i / 8] |= 1 << (i % 8);
        }
    }
    out
}

// =========================================================================
// Vandaele V_2 stack (Theorem 2 / Theorem 3 machinery).
//
// Layered bottom-up:
//   - l2_naive        Definition 2.3 (Eq. 5) for k=2: CCX ladder on 2n+1 qubits.
//                     Ancilla-free, Θ(n) gates, O(n) depth.
//                     (Paper's Lemma 4 gives log-depth via n ancs; we skip
//                     that optimization — we care about gate count + ancs,
//                     not depth. All paper ops bounds still hold.)
//   - v2_naive        Definition 2.4 (Eq. 6) for k=2: V-shape of two L_2
//                     ladders. Ancilla-free.
//   - compare_geq_v2  Theorem 3 via V_2 per Eq. 30-32: X-mask + slice-2
//                     structure with dirty-anc wiring. 1 dirty ancilla.
//
// All operations use only {CCX, CX, X} and are classically reversible.
// =========================================================================

/// `L_2^(n)` operator (Vandaele Def 2.3, Eq. 5 for k=2).
///
/// Acts on 2n+1 qubits `wire[0..=2n]` as a CCX ladder:
///   for i=1..n:  CCX(wire[2i-2], wire[2i-1], wire[2i])
///
/// Classically computes: wire[2i] ^= prefix-AND-pattern. Ancilla-free.
/// Gates: n CCX. Depth: O(n) (log-depth via Lemma 4 deferred).
/// Self-inverse since CCX is self-inverse and gates target non-overlapping
/// positions' targets (each wire[2i] is touched by one CCX).

fn v2_naive_refs(circ: &mut Circuit, wire: &[&QReg]) {
    let len = wire.len();
    assert!(
        len >= 3 && len % 2 == 1,
        "V_2 needs 2n+1 qubits (n≥1), got {len}"
    );
    let n = (len - 1) / 2;
    // L_2^(n-1) forward ladder on wire[0..2n-1]:
    //   CCX(wire[0], wire[1], wire[2]); CCX(wire[2], wire[3], wire[4]); ...
    //   ...; CCX(wire[2n-4], wire[2n-3], wire[2n-2]).
    for i in 1..n {
        circ.ccx(wire[2 * i - 2], wire[2 * i - 1], wire[2 * i]);
    }
    // Middle CCX on the last triple:
    //   CCX(wire[2n-2], wire[2n-1], wire[2n]).
    circ.ccx(wire[2 * n - 2], wire[2 * n - 1], wire[2 * n]);
    // L_2^(n-1) reverse ladder (same gates in reverse; CCX is self-inverse):
    //   CCX(wire[2n-4], wire[2n-3], wire[2n-2]); ...; CCX(wire[0], wire[1], wire[2]).
    for i in (1..n).rev() {
        circ.ccx(wire[2 * i - 2], wire[2 * i - 1], wire[2 * i]);
    }
}

/// Variant of [`v2_naive`] that frees `wire[last]` right after the
/// middle CCX (its last gate-touch). The reverse ladder only touches
/// `wire[2..2n-2]`, so `wire[2n]` can be freed before it. Caller passes
/// the trailing wire (z) by value so we can drop it in place.
#[allow(dead_code)]
fn v2_naive_free_last(circ: &mut Circuit, wire_prefix: &[&QReg], z: QReg) {
    let len = wire_prefix.len() + 1;
    assert!(
        len >= 3 && len % 2 == 1,
        "v2_naive_free_last: needs 2n+1 qubits, got {len}"
    );
    let n = (len - 1) / 2;
    for i in 1..n {
        circ.ccx(
            wire_prefix[2 * i - 2],
            wire_prefix[2 * i - 1],
            wire_prefix[2 * i],
        );
    }
    // Middle CCX — last gate touching z = wire[2n].
    circ.ccx(wire_prefix[2 * n - 2], wire_prefix[2 * n - 1], &z);
    // Drop z right after its last touch; the reverse ladder doesn't touch z.
    drop(z);
    for i in (1..n).rev() {
        circ.ccx(
            wire_prefix[2 * i - 2],
            wire_prefix[2 * i - 1],
            wire_prefix[2 * i],
        );
    }
}

/// Quantum-quantum comparator per Vandaele 2026 Fig 5.
///
/// Empirical semantic (traced for n=2 exhaustive): `z ^= 1[a > b]`
/// with register convention `a[0], b[0]` at top, `a[n-1], b[n-1]`
/// at bottom. Paper labels output as `z ⊕ (a < b)`, but trace
/// shows `1[a > b]` — the sign discrepancy is likely paper
/// labeling convention; either interpretation is trivially
/// reversible via `z ^= 1` post-compare.
///
/// Ancilla: **0** (paper's Theorem 2). Gates: O(n) total.
///
/// Fig 5 structure (n=5 example, generalizes):
///   Slice 1:
///     (a) X on every `b_i` (n X)
///     (b) `CX(a_i`, `b_i`) for i=1..n-1 (n-1 CX)
///     (c) CX(a_{n-1}, z) — captures MSB carry
///     (d) CX ladder on a: CX(a_{i-1}, `a_i`) for i = n-1 down to 2 (n-2 CX)
///   Slice 2 (the `V_2` operator):
///     Palindromic CCX chain on interleaved wire [`a_0`, `b_0`, ..., a_{n-1}, b_{n-1}, z]
///     = `CCX(a_0,b_0,a_1)`; `CCX(a_1,b_1,a_2)`; ...; CCX(a_{n-1},b_{n-1},z); reverse
///     = 2n-1 CCX (via `v2_naive`)
///   Slice 3: inverse of slice 1 EXCEPT col 4 (CX(a_{n-1},z) stays — it's the output).
///
/// Total: (2n-1) CCX + (4n-3) CX + 2n X. 0 ancillae. O(n) gates.

/// Quantum-quantum comparator via Cuccaro MAJ/reverse-MAJ with 1 dirty
/// ancilla. Retained for comparison; the paper's `compare_lt_qq_paper`
/// is ancilla-free.
///
/// Builds `b + (~a)` via Cuccaro MAJ cascade, extracts carry into z
/// (which = 1 iff b > a iff a < b since we skip the +1), then reverses
/// the MAJ cascade to restore a and b.
///
/// Not the paper's V_2-based Θ(log n)-depth Theorem 2 — we trade depth
/// for simplicity. Ops budget (the constraint that matters here) is
/// Θ(n) either way. a and b preserved; z XOR-ed.

/// Build the wire sequence for the top-half `V_2^(h)` call in Eq 32.
///
/// Slots (length 2h+1):
///   [`g_0`, `a_0_data`, `g_1`, `a_1_data`, ..., g_{h-1}, `last_data`, target]
/// with `g_i` <- a[n-h+i] (bottom-half a's play dirty g-slots), `a_i_data` <- a[i]
/// for i=0..h-2, `last_data` = anc0 (holds AND(a[h-1..n])), target = anc1 = z.
fn build_top_wires_refs<'a>(
    a: &[&'a QReg],
    anc0: &'a QReg,
    anc1: &'a QReg,
    h: usize,
    n: usize,
) -> Vec<&'a QReg> {
    debug_assert!(h >= 1 && n >= h);
    let mut w = Vec::with_capacity(2 * h + 1);
    for i in 0..h - 1 {
        w.push(a[n - h + i]);
        w.push(a[i]);
    }
    w.push(a[n - 1]);
    w.push(anc0);
    w.push(anc1);
    w
}

/// Build the wire sequence for the bottom-half `V_2^(l)` call in Eq 32.
///
/// Slots (length 2l+1):
///   [`g_0`, `a_h`, `g_1`, a_{h+1}, ..., g_{l-1}, a_{n-1}, target]
/// with `g_i` <- a[i] (top-half a's play dirty g-slots), data <- a[h+i]
/// for i=0..l-1, target = anc1 = z.
fn build_bot_wires_refs<'a>(a: &[&'a QReg], anc1: &'a QReg, h: usize, l: usize) -> Vec<&'a QReg> {
    debug_assert!(l >= 1);
    let mut w = Vec::with_capacity(2 * l + 1);
    for i in 0..l {
        w.push(a[i]);
        w.push(a[h + i]);
    }
    w.push(anc1);
    w
}

/// Emit the Eq 32 `V_2` decomposition of Fig 7's slice 2 multi-ctrl X cascade,
/// with Fig 2(a) clean→dirty upgrade so anc0 is dirty.
///
/// Structure (see `notes/theorem3_eq32_gates.md)`:
///   glue C^(l+1)X(a[h-1..n]; anc0);
///   (`V_2^(h)` on `top_wires`; `top_cXOR_wall)^2`;     // ctrl-U #1 per Fig 2(a)
///   glue C^(l+1)X(a[h-1..n]; anc0);               // uncompute/re-toggle
///   (`V_2^(h)` on `top_wires`; `top_cXOR_wall)^2`;     // ctrl-U #2 per Fig 2(a)
///   (`V_2^(l)` on `bot_wires`; `bot_cXOR_wall)^2`;     // bottom, no anc0 involvement
///
/// Preconditions: n >= 2. Caller has emitted slice 1 (X-mask + c-CX), col 4
/// X(z)-iff-c_{n-1}=1, and computed `c_eff` (the ladder-updated classical
/// values) before calling this.
///
/// Glue gate: `C^(l+1)X` — for k <= 5 uses existing `mcx_dirty` with a[0] as
/// psi (a[0] is always outside `glue_ctrls` = a[h-1..n] when n >= 4; for
/// n=2,3 k<=2 and no psi is needed).
/// For k >= 6 (n >= 10), awaits a separate multi-dirty extension.
/// Clean-ancilla variant of `slice2_eq32`: `anc0` must be |0⟩ on entry
/// (and will be returned to |0⟩). Skips Fig 2(a)'s dirty-ancilla
/// doubling — the top half runs `glue · ctrl-U · glue` (one pair)
/// instead of `glue · ctrl-U · glue · ctrl-U` (two pairs). Halves
/// the top-block cost.
fn slice2_eq32_clean_refs(circ: &mut Circuit, a: &[&QReg], c_eff: &[u8], z: &QReg) {
    let n = a.len();
    debug_assert_eq!(c_eff.len(), n);
    debug_assert!(n >= 2, "slice2_eq32_clean requires n >= 2");

    let h = n.div_ceil(2);
    let l = n / 2;

    // When c_eff[0..h] are all zero, the top block's cXOR wall emits no
    // X gates, making each emit_top_block call identical. Two consecutive
    // v2_naive calls (a palindromic CCX sequence) would produce the same
    // gate at the seam — triggering the redundancy detector. Since
    // V_2 · V_2 = identity (V_2 is self-inverse), skipping both top-block
    // calls (and the glue gates that exist only to enable them) is correct.
    let top_c_zero = c_eff[..h].iter().all(|&x| x == 0);

    let glue_ctrls: Vec<&QReg> = a[h - 1..n].to_vec();

    if !top_c_zero {
        // Allocate anc0 and free it inside mcx_dirty_any_k_consume RIGHT
        // AFTER the second glue restores it to |0>. The bot_block half
        // doesn't touch anc0, so leaving it live there burns the gap.
        let anc0 = circ.alloc_qreg("t3_anc0");
        {
            let top_wires: Vec<&QReg> = build_top_wires_refs(a, &anc0, z, h, n);
            let emit_top_block = |circ: &mut Circuit, top_wires: &[&QReg]| {
                v2_naive_refs(circ, top_wires);
                for i in 0..h {
                    if c_eff[i] == 1 {
                        circ.x(a[n - h + i]);
                    }
                }
            };

            // Clean-anc form: glue · ctrl-U · glue. ctrl-U = (V_2 · cXOR)^2
            // = 2 top_blocks. Total top half: 2 glue + 2 top_blocks (half the
            // dirty version's 2 glue + 4 top_blocks).
            mcx_dirty_any_k(circ, &glue_ctrls, &anc0, a[0]);
            emit_top_block(circ, &top_wires);
            emit_top_block(circ, &top_wires);
        }
        // top_wires dropped — anc0 movable.
        // Second glue using the consume variant: anc0 is freed right after
        // the last gate that touches it (inside mcx_dirty_any_k), before
        // any trailing X-wrap ops that restore glue_ctrls.
        mcx_dirty_any_k_consume(circ, &glue_ctrls, anc0, a[0]);
    }

    // Bottom half — no ancilla involvement.
    // Same guard: when c_eff[h..n] are all zero, both bot-block calls are
    // identical (empty cXOR wall), V_2 · V_2 = identity, skip both.
    //
    // The second bot-block's cXOR wall is NOT emitted here. It is merged
    // with the caller's slice3 inversions: for each i in 0..l, the second
    // cXOR would emit X(a[i]) iff c_eff[h+i]=1, and the standard slice3
    // emits X(a[i]) iff (i=0 OR c[i]=0). When both would fire they cancel
    // (net 0); when only one fires the caller emits the single net X.
    // This avoids the redundancy-detector panic: v2_naive call 2 ends with
    // a CCX on a[i], so after V_2 the last_op_for(a[i]) is CCX (not X),
    // and the merged slice3 X(a[i]) is never seen as adjacent to an X.
    // See compare_lt_cq_paper_refs for the merged slice3 emission.
    if l >= 1 && c_eff[h..].iter().any(|&x| x != 0) {
        let bot_wires: Vec<&QReg> = build_bot_wires_refs(a, z, h, l);
        // First bot_block: V_2 then cXOR (full).
        v2_naive_refs(circ, &bot_wires);
        for i in 0..l {
            if c_eff[h + i] == 1 {
                circ.x(a[i]);
            }
        }
        // Second bot_block: V_2 only; caller emits the merged cXOR+slice3.
        v2_naive_refs(circ, &bot_wires);
    }
}

/// Classical-quantum comparator per Vandaele 2026 Theorem 3 / Fig 7.
///
/// Computes: `z ^= 1[a < c]` under the convention that a[0] and c[0] are
/// the LSBs. Paper labels the Fig 7 output as `z ⊕ (c < a)`, but
/// exhaustive tracing shows the actual semantic is `z ⊕ (a < c)` — same
/// convention flip observed in Fig 5 (where paper labels `z ⊕ (a < b)`
/// but the circuit computes `z ^= 1[a > b]`). Either interpretation is
/// a correct comparator; callers can XOR z with 1 to invert.
///
/// Ancilla budget: **1 dirty** (supplied as `dirty` param). Caller must
/// pass a qubit whose state is allowed to be arbitrary before the call;
/// it is restored to its entry state at the end.
///
/// Gate count: O(n) for n <= 9 (given current `mcx_dirty` k <= 5 support).
/// For n >= 10 the glue C^(l+1)X needs an extension of `mcx_dirty` —
/// not yet wired.
///
/// Implementation:
///   Slice 1: X-mask on a (unconditional + c_i-guarded for i>=1).
///   Slice 2: col 4 X(z) iff c_{n-1}=1; compile-time c-ladder; Eq 32
///            `V_2` decomposition with Fig 2(a) clean→dirty upgrade.
///   Slice 3: inverse of slice 1.
///
/// See `notes/theorem3_progress.md` and `notes/theorem3_eq32_gates.md` for
/// the full gate-by-gate derivation from the Vandaele `TikZ` source.
pub fn compare_lt_cq_paper(circ: &mut Circuit, a: &[QReg], c: &[u8], z: &QReg) {
    let a_refs: Vec<&QReg> = a.iter().collect();
    compare_lt_cq_paper_refs(circ, &a_refs, c, z);
}

/// Reference-slice variant of [`compare_lt_cq_paper`].
pub(crate) fn compare_lt_cq_paper_refs(circ: &mut Circuit, a: &[&QReg], c: &[u8], z: &QReg) {
    let n = a.len();
    assert_eq!(c.len(), n, "compare_lt_cq_paper: a/c length mismatch");

    if n == 0 {
        return;
    }
    // c=0 short-circuit: 1[a < 0] = 0 always, so z unchanged. The
    // general path emits slice1 X's (on a) and slice3 X's that exactly
    // invert them when slice2 is also empty (which it is when
    // c_eff is fully zero). Without this short-circuit the slice1 and
    // slice3 X's appear adjacent on a-registers and trigger the
    // redundancy detector.
    if c.iter().all(|&x| x == 0) {
        return;
    }
    if n == 1 {
        // z ^= 1[a < c] = ~a[0] · c[0].
        // Fig 7's n=5 structure degenerates incorrectly at n=1 (trace
        // shows it computes c·a, not ~a·c), so we special-case.
        if c[0] == 1 {
            circ.x(a[0]);
            circ.cx(a[0], z);
            circ.x(a[0]);
        }
        return;
    }

    // Compile-time c-ladder (computed first, needed for top_c_zero/bot_c_zero
    // which determine the semi-isolated index skip before slice1).
    //   For j = n-1 down to 2: c_j ^= c_{j-1}.
    //   c_0 and c_1 are unchanged.
    let mut c_eff: Vec<u8> = c.to_vec();
    for j in (2..n).rev() {
        c_eff[j] ^= c_eff[j - 1];
    }

    let h = n.div_ceil(2);
    let l = n / 2;
    let top_c_zero = c_eff[..h].iter().all(|&x| x == 0);

    // Semi-isolated index (odd n only): the index n-h = (n-1)/2 appears in
    // top_wires but NOT in bot_wires. When top_c_zero is true (top block
    // AND its glue are skipped entirely) and the bot block runs, a[n-h] is
    // never touched by any V_2 or glue gate between slice1 and slice3.
    // The slice1 X(a[n-h]) and slice3 X(a[n-h]) form a pure no-op pair
    // (a[n-h] serves no role in the active V_2 computation). Removing both
    // preserves semantics and eliminates the adjacent X-X redundancy.
    //
    // For even n: the corresponding "bot-only" index is a[h-1], but for
    // even n when bot_c_zero=true the top block always runs, and the glue
    // (which uses a[h-1..n] as controls, including a[h-1]) touches a[h-1]
    // between slice1 and slice3. So no adjacency, no skip needed.
    let semi_idx = n - h; // = (n-1)/2 for odd n (unused for even n)
    let skip_semi = n % 2 == 1 && top_c_zero && c[n - h] == 0;

    // Slice 1 cells 2+3 merged: net effect a[i] ^= (1 XOR c[i]).
    // Emitting cell 2 (X all a_i) then cell 3 (X(a_i) iff c[i]=1, i>=1)
    // produces adjacent X-X on a[i] when c[i]=1, rejected by the detector.
    // Merged: emit X(a[i]) only when net flip is odd: always for i=0
    // (cell 3 skips i=0), and for i>=1 only when c[i]=0.
    // Semi-isolated index (odd n, top_c_zero case) skipped entirely.
    circ.x(a[0]);
    for i in 1..n {
        if i == semi_idx && skip_semi {
            continue;
        }
        if c[i] == 0 {
            circ.x(a[i]);
        }
    }

    // Slice 2 cells (n+3)..(2n+2): the multi-ctrl X cascade, emitted
    // via Eq 32 V_2 decomposition. Alloc a CLEAN anc0 internally
    // (slice2_eq32_clean handles it) — we trade +1 peak ancilla per
    // recursion level (polylog aggregated) for halving all ops vs
    // the dirty-ancilla doubling Fig 2(a) requires.
    //
    // NOTE: slice2_eq32_clean_refs does NOT emit the second bot-block cXOR
    // wall. It emits (V_2·cXOR)·V_2 for the bot half. The deferred cXOR is
    // merged with slice3 below to avoid adjacent X-X at the seam.
    slice2_eq32_clean_refs(circ, a, &c_eff, z);

    // Slice 2 cell 4: X(z) iff c_{n-1}=1 (uses ORIGINAL c_{n-1}, pre-ladder).
    // Emitted AFTER slice2_eq32_clean_refs rather than before, because
    // compare_geq emits X(out=z) just before calling compare_lt, and X(z)
    // col4 would produce adjacent X-X on z (separated only by slice1's X ops
    // on a-registers). Since X(z) commutes with V_2 (z is only a target in
    // V_2, never a control), reordering to post-V_2 is semantics-preserving.
    // After slice2_eq32_clean_refs, last_op_for(z) is a CCX from v2_naive,
    // so X(z) here sees no adjacency with the prior X(z) from compare_geq.
    if c[n - 1] == 1 {
        circ.x(z);
    }

    // Merged (bot_cXOR_call2 · slice3_standard), emitted high-to-low.
    //
    // slice2_eq32_clean_refs deferred the second bot-block cXOR wall
    // (indices 0..l, fires for c_eff[h+i]=1). slice3_standard flips a[i]
    // iff (i=0 OR c[i]=0). Net flip for each i = XOR of both; if net=1 we
    // emit X(a[i]), otherwise the two cancel and neither is emitted.
    //
    // After v2_naive_call2 (the last gate of the deferred bot-block),
    // last_op_for(a[i]) is a CCX — not X — so the first X here never
    // triggers the redundancy detector. High-to-low order ensures no
    // internal adjacency: each a[i] is touched at most once in this section.
    //
    // Semi-isolated indices are also skipped here (matching slice1 skip).
    let bot_ran = l >= 1 && c_eff[h..].iter().any(|&x| x != 0);
    // High indices (l..n): only slice3_standard contributes (no bot cXOR).
    for i in (l..n).rev() {
        // i is always >= 1 here (l >= 1 when n >= 2, and i >= l >= 1).
        if i == semi_idx && skip_semi {
            continue;
        }
        if c[i] == 0 {
            circ.x(a[i]);
        }
    }
    // Low indices (0..l): merge deferred bot_cXOR_call2 with slice3_standard.
    for i in (0..l).rev() {
        if i == semi_idx && skip_semi {
            continue;
        }
        let cxor2 = bot_ran && c_eff[h + i] == 1;
        let s3 = i == 0 || c[i] == 0;
        if cxor2 ^ s3 {
            circ.x(a[i]);
        }
    }
}

/// Reference-slice variant of the free-z compare (frees `z` right after
/// its last gate-touch, before any trailing X-wrap ops on `a`).
pub(crate) fn compare_lt_cq_paper_free_z_refs(circ: &mut Circuit, a: &[&QReg], c: &[u8], z: QReg) {
    let n = a.len();
    debug_assert_eq!(c.len(), n);
    if n == 0 {
        return;
    }
    if n == 1 {
        if c[0] == 1 {
            circ.x(a[0]);
            circ.cx(a[0], &z);
            // Last touch on z is the cx above; drop now.
            drop(z);
            circ.x(a[0]);
        }
        return;
    }

    // Compile-time c-ladder and isolation analysis (same as compare_lt_cq_paper_refs).
    let mut c_eff: Vec<u8> = c.to_vec();
    for j in (2..n).rev() {
        c_eff[j] ^= c_eff[j - 1];
    }

    let h = n.div_ceil(2);
    let l = n / 2;
    let top_c_zero = c_eff[..h].iter().all(|&x| x == 0);

    // Semi-isolated index (odd n only, same logic as compare_lt_cq_paper_refs).
    let semi_idx = n - h; // = (n-1)/2 for odd n
    let skip_semi = n % 2 == 1 && top_c_zero && c[n - h] == 0;

    // Slice 1 cells 2+3 merged (same merge as compare_lt_cq_paper_refs,
    // plus skip_semi for the semi-isolated index).
    circ.x(a[0]);
    for i in 1..n {
        if i == semi_idx && skip_semi {
            continue;
        }
        if c[i] == 0 {
            circ.x(a[i]);
        }
    }
    // Note: col4 X(z) is NOT emitted here. See the col4 parameter below.

    // Slice 2: use the free-z variant so z is freed inside.
    // NOTE: does NOT emit the second bot-block cXOR wall (same as
    // slice2_eq32_clean_refs); it is merged with slice3 below.
    // Pass col4 = c[n-1]==1 so the X(z) flip is deferred to inside
    // slice2 (after the first V_2 gate on z), avoiding the adjacent
    // X-X that would occur if emitted here after compare_geq's X(out).
    let col4 = c[n - 1] == 1;
    slice2_eq32_clean_free_z_refs(circ, a, &c_eff, z, col4);

    // Merged (bot_cXOR_call2 · slice3_standard), same logic as in
    // compare_lt_cq_paper_refs. Semi-isolated indices skipped here too.
    let bot_ran = l >= 1 && c_eff[h..].iter().any(|&x| x != 0);
    for i in (l..n).rev() {
        if i == semi_idx && skip_semi {
            continue;
        }
        if c[i] == 0 {
            circ.x(a[i]);
        }
    }
    for i in (0..l).rev() {
        if i == semi_idx && skip_semi {
            continue;
        }
        let cxor2 = bot_ran && c_eff[h + i] == 1;
        let s3 = i == 0 || c[i] == 0;
        if cxor2 ^ s3 {
            circ.x(a[i]);
        }
    }
}

/// Variant of `slice2_eq32_clean` that frees `z` right after its last
/// gate-touch (the middle CCX of the second bot-block `v2_naive` call).
///
/// `col4`: if true, emit X(z) after the first `v2_naive` gate on z and
/// before `v2_naive_free_last`. This is the Vandaele "cell 4" flip, which
/// the caller (`compare_lt_cq_paper_free_z_refs`) cannot safely emit before
/// calling this function when `compare_geq`'s X(out=z) immediately precedes
/// and would produce an adjacent X-X pair on z (triggering the redundancy
/// detector). Since X(z) commutes with `V_2` (z is only a target, never a
/// control), deferring it to this interior position is semantics-preserving.
fn slice2_eq32_clean_free_z_refs(
    circ: &mut Circuit,
    a: &[&QReg],
    c_eff: &[u8],
    z: QReg,
    col4: bool,
) {
    let n = a.len();
    let h = n.div_ceil(2);
    let l = n / 2;

    // Guard: same logic as slice2_eq32_clean_refs. When c_eff[0..h] are all
    // zero the top-block cXOR wall is empty, making two consecutive v2_naive
    // calls produce adjacent identical CCXs at the seam. V_2·V_2=identity,
    // so both calls (plus the glue pair) can be skipped entirely.
    let top_c_zero = c_eff[..h].iter().all(|&x| x == 0);
    let bot_c_zero = l >= 1 && c_eff[h..].iter().all(|&x| x == 0);

    let glue_ctrls: Vec<&QReg> = a[h - 1..n].to_vec();
    // z stays alive only through the top block's v2_naive calls when
    // bot doesn't use it. The second glue (mcx_dirty_any_k_consume)
    // allocates internal mcxk_t/t3 ancillae which would advance
    // last_alloc_op_idx past z's last touch (v2_naive_call2's middle
    // CCX). Wrap z in Option so we can drop it early when bot is
    // skipped, before the offending allocs.
    let z_used_in_bot = l >= 1 && !bot_c_zero;
    let mut z_holder: Option<QReg> = Some(z);

    if !top_c_zero {
        let anc0 = circ.alloc_qreg("t3_anc0");
        {
            let z_ref = z_holder.as_ref().expect("z alive entering top block");
            let top_wires: Vec<&QReg> = build_top_wires_refs(a, &anc0, z_ref, h, n);
            let emit_top_block = |circ: &mut Circuit, top_wires: &[&QReg]| {
                v2_naive_refs(circ, top_wires);
                for i in 0..h {
                    if c_eff[i] == 1 {
                        circ.x(a[n - h + i]);
                    }
                }
            };

            mcx_dirty_any_k(circ, &glue_ctrls, &anc0, a[0]);
            emit_top_block(circ, &top_wires);
            emit_top_block(circ, &top_wires);
        }
        // top_wires borrow released. If bot won't use z, drop it now —
        // BEFORE the second glue's mcx_dirty_any_k_consume internal
        // allocs advance last_alloc_op_idx past z's last touch.
        if !z_used_in_bot {
            let z_owned = z_holder.take().expect("z still alive after top");
            if col4 {
                circ.x(&z_owned);
            }
            drop(z_owned);
        }
        mcx_dirty_any_k_consume(circ, &glue_ctrls, anc0, a[0]);
    }

    if z_used_in_bot {
        let z_owned = z_holder.take().expect("z used in bot");
        // Build bot_prefix (without z) up-front from `a` alone — z
        // appears only as the trailing entry and we need to free it
        // early via v2_naive_free_last.
        let bot_prefix: Vec<&QReg> = {
            let mut w: Vec<&QReg> = Vec::with_capacity(2 * l);
            for i in 0..l {
                w.push(a[i]);
                w.push(a[h + i]);
            }
            w
        };
        // First bot_block: emit full v2_naive (with z), then cXOR wall.
        {
            let mut bot_wires: Vec<&QReg> = bot_prefix.clone();
            bot_wires.push(&z_owned);
            v2_naive_refs(circ, &bot_wires);
        }
        for i in 0..l {
            if c_eff[h + i] == 1 {
                circ.x(a[i]);
            }
        }
        if col4 {
            circ.x(&z_owned);
        }
        // Second bot_block: V_2 only (frees z at last gate-touch);
        // the second cXOR is deferred to the caller's merged slice3.
        v2_naive_free_last(circ, &bot_prefix, z_owned);
    } else if let Some(z_owned) = z_holder.take() {
        // top_c_zero AND bot_c_zero (c_eff all-zero) — caller should
        // have short-circuited, but handle defensively.
        if col4 {
            circ.x(&z_owned);
        }
        drop(z_owned);
    }
}

#[cfg(test)]
mod v2_tests {
    use super::*;

    fn run_compare_lt_cq_paper(n: usize, a_val: u64, c_val: u64) {
        let mut circ = Circuit::new();
        let a: Vec<QReg> = (0..n)
            .map(|i| circ.alloc_qreg(&format!("a{}", i)))
            .collect();
        let z = circ.alloc_qreg("z");
        let dirty = circ.alloc_qreg("dirty");
        {
            let mut bytes = vec![0u8; n.div_ceil(8)];
            for i in 0..n {
                if (a_val >> i) & 1 == 1 {
                    bytes[i / 8] |= 1u8 << (i % 8);
                }
            }
            circ.sim_load_reg_bytes_shot(&a, &bytes, 0);
        }
        // Capture dirty's initial value (random per input bit 0 ^ 1 here
        // we just use alloc_input_qubit which gives |0>; include a flip
        // sometimes to exercise nonzero psi).
        let dirty_init_bit = (a_val ^ c_val).wrapping_mul(0x9E37_79B9_u64) & 1;
        circ.sim_load_reg_bytes_shot(std::slice::from_ref(&dirty), &[dirty_init_bit as u8], 0);

        let c: Vec<u8> = (0..n).map(|i| ((c_val >> i) & 1) as u8).collect();
        compare_lt_cq_paper(&mut circ, &a, &c, &z);
        let mut outputs: Vec<crate::point_add::trailmix_port::circuit::QReg> = Vec::new();
        outputs.extend(a);
        outputs.push(z);
        outputs.push(dirty);
        let (sim, detached) = circ.destroy_sim(outputs);
        let a_d = &detached[..n];
        let z_d = &detached[n];
        let dirty_d = &detached[n + 1];
        let got_a: u64 = (0..n).map(|i| (sim.qubit_mask(&a_d[i]) & 1) << i).sum();
        let got_z = sim.qubit_mask(z_d) & 1;
        let got_dirty = sim.qubit_mask(dirty_d) & 1;
        // Empirical semantic (confirmed by Fig 7 TikZ trace for n=2):
        // z ^= 1[a < c]. Paper labels z ⊕ (c<a) but trace shows a<c.
        let expected = if a_val < c_val { 1 } else { 0 };
        assert_eq!(
            got_z,
            expected,
            "compare_lt_cq_paper n={} a={} c={}: got_z={} exp={} (a<c={})",
            n,
            a_val,
            c_val,
            got_z,
            expected,
            a_val < c_val
        );
        assert_eq!(
            got_a, a_val,
            "a drift n={} a={} c={}: got_a={}",
            n, a_val, c_val, got_a
        );
        assert_eq!(
            got_dirty, dirty_init_bit,
            "dirty drift n={} a={} c={}: got={} init={}",
            n, a_val, c_val, got_dirty, dirty_init_bit
        );
        assert_eq!(
            sim.phase_mask(),
            0,
            "phase n={} a={} c={} dirty_init={}: sim_phase={:#x}",
            n,
            a_val,
            c_val,
            dirty_init_bit,
            sim.phase_mask()
        );
    }

    #[test]
    fn compare_lt_cq_paper_n1_all() {
        for a in 0..2 {
            for c in 0..2 {
                run_compare_lt_cq_paper(1, a, c);
            }
        }
    }
    #[test]
    fn compare_lt_cq_paper_n2_all() {
        for a in 0..4 {
            for c in 0..4 {
                run_compare_lt_cq_paper(2, a, c);
            }
        }
    }
    #[test]
    fn compare_lt_cq_paper_n3_all() {
        for a in 0..8 {
            for c in 0..8 {
                run_compare_lt_cq_paper(3, a, c);
            }
        }
    }
    #[test]
    fn compare_lt_cq_paper_n4_all() {
        for a in 0..16 {
            for c in 0..16 {
                run_compare_lt_cq_paper(4, a, c);
            }
        }
    }
    #[test]
    fn compare_lt_cq_paper_n5_all() {
        for a in 0..32 {
            for c in 0..32 {
                run_compare_lt_cq_paper(5, a, c);
            }
        }
    }
    #[test]
    fn compare_lt_cq_paper_n8_sample() {
        let mix = |s: u64| {
            s.wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407)
        };
        for seed in 0..256u64 {
            let r = mix(seed);
            run_compare_lt_cq_paper(8, r & 0xFF, (r >> 8) & 0xFF);
        }
    }
    #[test]
    fn compare_lt_cq_paper_n9_sample() {
        let mix = |s: u64| {
            s.wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407)
        };
        for seed in 0..256u64 {
            let r = mix(seed);
            run_compare_lt_cq_paper(9, r & 0x1FF, (r >> 9) & 0x1FF);
        }
    }
    #[test]
    fn compare_lt_cq_paper_n10_sample() {
        // n=10 triggers glue k=6, which recurses into compare_lt_cq_paper(6).
        let mix = |s: u64| {
            s.wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407)
        };
        for seed in 0..256u64 {
            let r = mix(seed);
            run_compare_lt_cq_paper(10, r & 0x3FF, (r >> 10) & 0x3FF);
        }
    }
    #[test]
    fn compare_lt_cq_paper_n12_sample() {
        let mix = |s: u64| {
            s.wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407)
        };
        for seed in 0..256u64 {
            let r = mix(seed);
            run_compare_lt_cq_paper(12, r & 0xFFF, (r >> 12) & 0xFFF);
        }
    }
    #[test]
    fn compare_lt_cq_paper_n16_sample() {
        let mix = |s: u64| {
            s.wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407)
        };
        for seed in 0..128u64 {
            let r = mix(seed);
            run_compare_lt_cq_paper(16, r & 0xFFFF, (r >> 16) & 0xFFFF);
        }
    }
    #[test]
    fn compare_lt_cq_paper_n32_sample() {
        let mix = |s: u64| {
            s.wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407)
        };
        for seed in 0..64u64 {
            let r = mix(seed);
            run_compare_lt_cq_paper(32, r & 0xFFFFFFFF, (r >> 32) ^ (r & 0xDEADBEEF));
        }
    }

    fn run_compare_lt_cq_paper_wide(n: usize, a_bits: &[bool], c_bits: &[bool]) {
        assert_eq!(a_bits.len(), n);
        assert_eq!(c_bits.len(), n);
        let mut circ = Circuit::new();
        let a: Vec<QReg> = (0..n)
            .map(|i| circ.alloc_qreg(&format!("a{}", i)))
            .collect();
        let z = circ.alloc_qreg("z");
        let dirty = circ.alloc_qreg("dirty");
        {
            let mut bytes = vec![0u8; n.div_ceil(8)];
            for i in 0..n {
                if a_bits[i] {
                    bytes[i / 8] |= 1u8 << (i % 8);
                }
            }
            circ.sim_load_reg_bytes_shot(&a, &bytes, 0);
        }
        let dirty_init: u8 = (a_bits[0] ^ c_bits[0]) as u8;
        circ.sim_load_reg_bytes_shot(std::slice::from_ref(&dirty), &[dirty_init], 0);
        let c: Vec<u8> = c_bits.iter().map(|b| *b as u8).collect();

        let ops_before = circ.ops.len();
        compare_lt_cq_paper(&mut circ, &a, &c, &z);
        let ops_count = circ.ops.len() - ops_before;
        let mut outputs: Vec<crate::point_add::trailmix_port::circuit::QReg> = Vec::new();
        outputs.extend(a);
        outputs.push(z);
        outputs.push(dirty);
        let (sim, detached) = circ.destroy_sim(outputs);
        let a_d = &detached[..n];
        let z_d = &detached[n];
        let dirty_d = &detached[n + 1];
        let got_z = sim.qubit_mask(z_d) & 1;

        // Compare a < c as arbitrary-precision ints (MSB-first comparison
        // from index n-1 down to 0).
        let mut lt: u64 = 0;
        for i in (0..n).rev() {
            let a_bit = a_bits[i] as u8;
            let c_bit = c_bits[i] as u8;
            if a_bit != c_bit {
                lt = if a_bit < c_bit { 1 } else { 0 };
                break;
            }
        }
        assert_eq!(
            got_z, lt,
            "compare_lt_cq_paper n={}: got_z={} exp={}",
            n, got_z, lt
        );
        for i in 0..n {
            let got = (sim.qubit_mask(&a_d[i]) & 1) as u8;
            assert_eq!(got, a_bits[i] as u8, "a drift n={} bit {}", n, i);
        }
        let got_dirty = (sim.qubit_mask(dirty_d) & 1) as u8;
        assert_eq!(got_dirty, dirty_init, "dirty drift n={}", n);
        assert_eq!(sim.phase_mask(), 0, "phase n={}", n);
        println!("compare_lt_cq_paper n={} ops={}", n, ops_count);
    }

    #[test]
    fn compare_lt_cq_paper_n64_sample() {
        let mix = |s: u64| {
            s.wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407)
        };
        for seed in 0..4u64 {
            let r1 = mix(seed);
            let r2 = mix(seed ^ 0xAAAA);
            let a_bits: Vec<bool> = (0..64).map(|i| ((r1 >> (i & 63)) & 1) == 1).collect();
            let c_bits: Vec<bool> = (0..64).map(|i| ((r2 >> (i & 63)) & 1) == 1).collect();
            run_compare_lt_cq_paper_wide(64, &a_bits, &c_bits);
        }
    }
    #[test]
    fn compare_lt_cq_paper_n128_sample() {
        let mix = |s: u64| {
            s.wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407)
        };
        for seed in 0..2u64 {
            let r1 = mix(seed);
            let r2 = mix(seed ^ 0xAAAA);
            let a_bits: Vec<bool> = (0..128).map(|i| (mix(r1 ^ (i as u64)) & 1) == 1).collect();
            let c_bits: Vec<bool> = (0..128).map(|i| (mix(r2 ^ (i as u64)) & 1) == 1).collect();
            run_compare_lt_cq_paper_wide(128, &a_bits, &c_bits);
        }
    }
    #[test]
    fn compare_lt_cq_paper_n257_sample() {
        let mix = |s: u64| {
            s.wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407)
        };
        for seed in 0..1u64 {
            let r1 = mix(seed);
            let r2 = mix(seed ^ 0xAAAA);
            let a_bits: Vec<bool> = (0..257).map(|i| (mix(r1 ^ (i as u64)) & 1) == 1).collect();
            let c_bits: Vec<bool> = (0..257).map(|i| (mix(r2 ^ (i as u64)) & 1) == 1).collect();
            run_compare_lt_cq_paper_wide(257, &a_bits, &c_bits);
        }
    }
}

#[cfg(test)]
mod cond_inc_tests {
    use super::*;
    use crate::point_add::trailmix_port::circuit::Circuit;

    #[test]
    fn mcx_clean_k_ops_table() {
        for &k in &[5usize, 7, 9, 11, 13] {
            let mut circ = Circuit::new();
            let ctrls: Vec<QReg> = (0..k).map(|_| circ.alloc_qreg("c")).collect();
            for q in &ctrls {
                circ.x(q);
            }
            let t = circ.alloc_qreg("t");
            let t0 = circ.ops.len();
            let ctrl_refs: Vec<&QReg> = ctrls.iter().collect();
            mcx_clean_k(&mut circ, &ctrl_refs, &t);
            let ops = circ.ops.len() - t0;
            eprintln!("mcx_clean_k k={:>2} ops={:>5}", k, ops);
            drop(ctrl_refs);
            let mut outs = ctrls;
            outs.push(t);
            let _ = circ.destroy_sim(outs);
        }
    }

    fn run_cqadd_case(n: usize, x_init: u64, c_val: u64) {
        let mut circ = Circuit::new();
        let x: Vec<QReg> = (0..n)
            .map(|i| circ.alloc_qreg(&format!("x{}", i)))
            .collect();
        // Dirty ancilla g: start in |0>, preserved on exit.
        let g = circ.alloc_qreg("g_dirty");
        {
            let mut bytes = vec![0u8; n.div_ceil(8)];
            for i in 0..n {
                if (x_init >> i) & 1 == 1 {
                    bytes[i / 8] |= 1u8 << (i % 8);
                }
            }
            circ.sim_load_reg_bytes_shot(&x, &bytes, 0);
        }
        // Classical constant c as little-endian bytes covering n bits.
        let n_bytes = n.div_ceil(8);
        let mut c_bytes = vec![0u8; n_bytes];
        for i in 0..n {
            if (c_val >> i) & 1 == 1 {
                c_bytes[i / 8] |= 1 << (i % 8);
            }
        }
        classical_quantum_add(&mut circ, &x, &c_bytes, &g);

        let mut outputs: Vec<crate::point_add::trailmix_port::circuit::QReg> = Vec::new();
        outputs.extend(x);
        outputs.push(g);
        let (sim, detached) = circ.destroy_sim(outputs);
        let x_d = &detached[..n];
        let g_d = &detached[n];
        let got: u64 = (0..n).map(|i| (sim.qubit_mask(&x_d[i]) & 1) << i).sum();
        let mask = (1u64 << n) - 1;
        let expected = (x_init.wrapping_add(c_val)) & mask;
        assert_eq!(
            got,
            expected,
            "cqadd n={} x={:0w$b} c={:0w$b}: got={:0w$b} exp={:0w$b}",
            n,
            x_init,
            c_val,
            got,
            expected,
            w = n
        );
        assert_eq!(
            sim.qubit_mask(g_d) & 1,
            0,
            "cqadd n={} x={:0w$b} c={:0w$b}: g leaked ({})",
            n,
            x_init,
            c_val,
            sim.qubit_mask(g_d) & 1,
            w = n
        );
        assert_eq!(
            sim.phase_mask(),
            0,
            "cqadd phase n={} x={} c={}: {:#x}",
            n,
            x_init,
            c_val,
            sim.phase_mask()
        );
    }

    fn run_mcx_dirty_case(k: usize, ctrls_val: u64, psi_val: u64, t_val: u64) {
        let mut circ = Circuit::new();
        let ctrls: Vec<QReg> = (0..k)
            .map(|i| circ.alloc_qreg(&format!("c{}", i)))
            .collect();
        let psi = circ.alloc_qreg("psi");
        let t = circ.alloc_qreg("t");
        {
            let mut bytes = vec![0u8; k.div_ceil(8)];
            for i in 0..k {
                if (ctrls_val >> i) & 1 == 1 {
                    bytes[i / 8] |= 1u8 << (i % 8);
                }
            }
            circ.sim_load_reg_bytes_shot(&ctrls, &bytes, 0);
        }
        circ.sim_load_reg_bytes_shot(std::slice::from_ref(&psi), &[psi_val as u8], 0);
        circ.sim_load_reg_bytes_shot(std::slice::from_ref(&t), &[t_val as u8], 0);
        let ctrl_refs: Vec<&QReg> = ctrls.iter().collect();
        mcx_dirty(&mut circ, &ctrl_refs, &t, &psi);
        drop(ctrl_refs);
        let mut outputs: Vec<crate::point_add::trailmix_port::circuit::QReg> = Vec::new();
        outputs.extend(ctrls);
        outputs.push(psi);
        outputs.push(t);
        let (sim, detached) = circ.destroy_sim(outputs);
        let ctrls_d = &detached[..k];
        let psi_d = &detached[k];
        let t_d = &detached[k + 1];
        let got_psi = sim.qubit_mask(psi_d) & 1;
        let got_t = sim.qubit_mask(t_d) & 1;
        let and_ctrls = (0..k).fold(1u64, |acc, i| acc & ((ctrls_val >> i) & 1));
        let expected_t = t_val ^ and_ctrls;
        assert_eq!(
            got_psi, psi_val,
            "mcx_dirty k={} ctrls={:b} psi={} t={}: psi corrupted (got {})",
            k, ctrls_val, psi_val, t_val, got_psi
        );
        assert_eq!(
            got_t, expected_t,
            "mcx_dirty k={} ctrls={:b} psi={} t={}: target got {} expected {}",
            k, ctrls_val, psi_val, t_val, got_t, expected_t
        );
        for (i, q) in ctrls_d.iter().enumerate() {
            let v = sim.qubit_mask(q) & 1;
            let exp = (ctrls_val >> i) & 1;
            assert_eq!(
                v, exp,
                "mcx_dirty k={} ctrl c{} changed: {} -> {}",
                k, i, exp, v
            );
        }
        assert_eq!(sim.phase_mask(), 0, "mcx_dirty k={} phase", k);
    }

    #[test]
    fn mcx_dirty_k3_all() {
        for bits in 0..(1u64 << 5) {
            let cv = bits & 7;
            let psi = (bits >> 3) & 1;
            let t = (bits >> 4) & 1;
            run_mcx_dirty_case(3, cv, psi, t);
        }
    }

    #[test]
    fn mcx_dirty_k4_all() {
        for bits in 0..(1u64 << 6) {
            let cv = bits & 0xF;
            let psi = (bits >> 4) & 1;
            let t = (bits >> 5) & 1;
            run_mcx_dirty_case(4, cv, psi, t);
        }
    }

    #[test]
    fn mcx_dirty_k5_all() {
        for bits in 0..(1u64 << 7) {
            let cv = bits & 0x1F;
            let psi = (bits >> 5) & 1;
            let t = (bits >> 6) & 1;
            run_mcx_dirty_case(5, cv, psi, t);
        }
    }

    fn run_mcx_clean_case(k: usize, ctrls_val: u64, t_val: u64) {
        let mut circ = Circuit::new();
        let ctrls: Vec<QReg> = (0..k)
            .map(|i| circ.alloc_qreg(&format!("c{}", i)))
            .collect();
        let t = circ.alloc_qreg("t");
        {
            let mut bytes = vec![0u8; k.div_ceil(8)];
            for i in 0..k {
                if (ctrls_val >> i) & 1 == 1 {
                    bytes[i / 8] |= 1u8 << (i % 8);
                }
            }
            circ.sim_load_reg_bytes_shot(&ctrls, &bytes, 0);
        }
        circ.sim_load_reg_bytes_shot(std::slice::from_ref(&t), &[t_val as u8], 0);
        let ops_before = circ.ops.len();
        let peak_before = circ.peak_qubits;
        let ctrl_refs: Vec<&QReg> = ctrls.iter().collect();
        mcx_clean_k(&mut circ, &ctrl_refs, &t);
        drop(ctrl_refs);
        let ops = circ.ops.len() - ops_before;
        let peak_delta = circ.peak_qubits.saturating_sub(peak_before);
        let mut outputs: Vec<crate::point_add::trailmix_port::circuit::QReg> = Vec::new();
        outputs.extend(ctrls);
        outputs.push(t);
        let (sim, detached) = circ.destroy_sim(outputs);
        let ctrls_d = &detached[..k];
        let t_d = &detached[k];
        let got_t = sim.qubit_mask(t_d) & 1;
        let and_ctrls = if k == 0 {
            1
        } else {
            (0..k).fold(1u64, |acc, i| acc & ((ctrls_val >> i) & 1))
        };
        let expected_t = t_val ^ and_ctrls;
        assert_eq!(
            got_t, expected_t,
            "mcx_clean_k k={} ctrls={:b} t={}: target got {} expected {} (ops={}, peak_delta={})",
            k, ctrls_val, t_val, got_t, expected_t, ops, peak_delta
        );
        for (i, q) in ctrls_d.iter().enumerate() {
            let v = sim.qubit_mask(q) & 1;
            let exp = (ctrls_val >> i) & 1;
            assert_eq!(
                v, exp,
                "mcx_clean_k k={} ctrl c{} changed: {} -> {}",
                k, i, exp, v
            );
        }
        assert_eq!(sim.phase_mask(), 0, "mcx_clean_k k={} phase", k);
    }

    /// Same validation as run_unary_case but for the KG log* variant.
    fn run_unary_ls_case(n: usize, v: u64, n_iters: usize) {
        use super::unary_iterate_log_star;
        let mut circ = Circuit::new();
        let c: Vec<QReg> = (0..n)
            .map(|i| circ.alloc_qreg(&format!("c{}", i)))
            .collect();
        let res: Vec<QReg> = (0..n)
            .map(|i| circ.alloc_qreg(&format!("r{}", i)))
            .collect();
        let fired = circ.alloc_qreg("fired");
        let mut bytes = vec![0u8; n.div_ceil(8)];
        for i in 0..n {
            if (v >> i) & 1 == 1 {
                bytes[i / 8] |= 1u8 << (i % 8);
            }
        }
        circ.sim_load_reg_bytes_shot(&c, &bytes, 0);
        let c_refs: Vec<&QReg> = c.iter().collect();
        let res_refs: Vec<&QReg> = res.iter().collect();
        let fired_ref = &fired;
        unary_iterate_log_star(&mut circ, &c_refs, n_iters, |circ, i, gate| {
            circ.cx(gate, fired_ref);
            for bit in 0..n {
                if (i >> bit) & 1 == 1 {
                    circ.cx(gate, res_refs[bit]);
                }
            }
        });
        drop(c_refs);
        drop(res_refs);
        let mut outs: Vec<QReg> = Vec::new();
        outs.extend(c);
        outs.extend(res);
        outs.push(fired);
        let (sim, det) = circ.destroy_sim(outs);
        let c_d = &det[..n];
        let res_d = &det[n..2 * n];
        let fired_d = &det[2 * n];
        let in_range = (v as usize) < n_iters;
        let mut res_v: u64 = 0;
        for (b, q) in res_d.iter().enumerate() {
            res_v |= (sim.qubit_mask(q) & 1) << b;
        }
        let expect_res = if in_range { v } else { 0 };
        assert_eq!(res_v, expect_res, "uls n={} v={} res", n, v);
        assert_eq!(
            sim.qubit_mask(fired_d) & 1,
            in_range as u64,
            "uls n={} v={} fired",
            n,
            v
        );
        let mut c_out: u64 = 0;
        for (b, q) in c_d.iter().enumerate() {
            c_out |= (sim.qubit_mask(q) & 1) << b;
        }
        assert_eq!(c_out, v, "uls n={} v={} counter not restored", n, v);
        assert_eq!(sim.phase_mask(), 0, "uls n={} v={} phase", n, v);
    }

    #[test]
    fn unary_iterate_log_star_full_range() {
        for n in 2..=6 {
            let l = 1usize << n;
            for v in 0..(1u64 << n) {
                run_unary_ls_case(n, v, l);
            }
        }
    }

    #[test]
    fn unary_iterate_log_star_partial_range() {
        for &(n, l) in &[(4usize, 11usize), (5, 20), (6, 50), (7, 100)] {
            for v in 0..(1u64 << n) {
                run_unary_ls_case(n, v, l);
            }
        }
    }

    #[test]
    fn mcx_clean_k3_all() {
        for bits in 0..(1u64 << 4) {
            let cv = bits & 7;
            let tv = (bits >> 3) & 1;
            run_mcx_clean_case(3, cv, tv);
        }
    }

    #[test]
    fn mcx_clean_k4_all() {
        for bits in 0..(1u64 << 5) {
            let cv = bits & 0xF;
            let tv = (bits >> 4) & 1;
            run_mcx_clean_case(4, cv, tv);
        }
    }

    #[test]
    fn mcx_clean_k5_all() {
        for bits in 0..(1u64 << 6) {
            let cv = bits & 0x1F;
            let tv = (bits >> 5) & 1;
            run_mcx_clean_case(5, cv, tv);
        }
    }

    #[test]
    fn mcx_clean_k6_all() {
        for bits in 0..(1u64 << 7) {
            let cv = bits & 0x3F;
            let tv = (bits >> 6) & 1;
            run_mcx_clean_case(6, cv, tv);
        }
    }

    #[test]
    fn mcx_clean_k7_all() {
        for bits in 0..(1u64 << 8) {
            let cv = bits & 0x7F;
            let tv = (bits >> 7) & 1;
            run_mcx_clean_case(7, cv, tv);
        }
    }

    #[test]
    fn mcx_clean_k8_all() {
        for bits in 0..(1u64 << 9) {
            let cv = bits & 0xFF;
            let tv = (bits >> 8) & 1;
            run_mcx_clean_case(8, cv, tv);
        }
    }

    #[test]
    fn mcx_clean_k10_sample() {
        // Full exhaustive would be 2^11 = 2048 cases; manageable but
        // sampling the boundary + random interior saves test time.
        for &cv in &[0u64, 0x3FFu64, 0x3FEu64, 0x1FFu64, 0x2AAu64, 0x155u64] {
            for tv in 0..2 {
                run_mcx_clean_case(10, cv, tv);
            }
        }
    }

    #[test]
    fn mcx_clean_k16_sample() {
        for &cv in &[0u64, 0xFFFFu64, 0xFFFEu64, 0x7FFFu64, 0xAAAAu64, 0x5555u64] {
            for tv in 0..2 {
                run_mcx_clean_case(16, cv, tv);
            }
        }
    }

    #[test]
    fn mcx_clean_k32_sample() {
        for &cv in &[0u64, 0xFFFFFFFFu64, 0xFFFFFFFEu64, 0xAAAAAAAAu64] {
            for tv in 0..2 {
                run_mcx_clean_case(32, cv, tv);
            }
        }
    }

    #[test]
    fn mcx_clean_k64_sample() {
        for &cv in &[0u64, u64::MAX, u64::MAX - 1, 0xAAAA_AAAA_AAAA_AAAAu64] {
            for tv in 0..2 {
                run_mcx_clean_case(64, cv, tv);
            }
        }
    }

    #[test]
    fn mcx_clean_k_bench_ops() {
        // Report ops counts for documentation / cost tracking.
        for &k in &[4usize, 6, 8, 16, 32, 64, 128] {
            let mut circ = Circuit::new();
            let ctrls: Vec<QReg> = (0..k)
                .map(|i| circ.alloc_qreg(&format!("c{}", i)))
                .collect();
            let t = circ.alloc_qreg("t");
            for q in &ctrls {
                circ.x(q);
            }
            let ops_before = circ.ops.len();
            let peak_before = circ.peak_qubits;
            let ctrl_refs: Vec<&QReg> = ctrls.iter().collect();
            mcx_clean_k(&mut circ, &ctrl_refs, &t);
            drop(ctrl_refs);
            let ops = circ.ops.len() - ops_before;
            let peak_delta = circ.peak_qubits.saturating_sub(peak_before);
            eprintln!(
                "mcx_clean_k k={:>4} ops={:>6} peak_delta={:>3}",
                k, ops, peak_delta
            );
            let mut outs = ctrls;
            outs.push(t);
            let _ = circ.destroy_sim(outs);
        }
    }

    #[test]
    fn cqadd_n4_all() {
        for x in 0..16 {
            for c in 0..16 {
                run_cqadd_case(4, x, c);
            }
        }
    }

    #[test]
    fn cqadd_n8_all() {
        for x in 0..256 {
            for c in 0..256 {
                run_cqadd_case(8, x, c);
            }
        }
    }

    fn run_ctrl_cqadd_case(n: usize, ctrl_bit: u64, x_init: u64, c_val: u64) {
        let mut circ = Circuit::new();
        let ctrl = circ.alloc_qreg("ctrl");
        let x: Vec<QReg> = (0..n)
            .map(|i| circ.alloc_qreg(&format!("x{}", i)))
            .collect();
        circ.sim_load_reg_bytes_shot(std::slice::from_ref(&ctrl), &[ctrl_bit as u8], 0);
        {
            let mut bytes = vec![0u8; n.div_ceil(8)];
            for i in 0..n {
                if (x_init >> i) & 1 == 1 {
                    bytes[i / 8] |= 1u8 << (i % 8);
                }
            }
            circ.sim_load_reg_bytes_shot(&x, &bytes, 0);
        }
        let n_bytes = n.div_ceil(8);
        let mut c_bytes = vec![0u8; n_bytes];
        for i in 0..n {
            if (c_val >> i) & 1 == 1 {
                c_bytes[i / 8] |= 1 << (i % 8);
            }
        }
        controlled_classical_quantum_add(&mut circ, &ctrl, &x, &c_bytes);
        let mut outputs: Vec<crate::point_add::trailmix_port::circuit::QReg> = Vec::new();
        outputs.push(ctrl);
        outputs.extend(x);
        let (sim, detached) = circ.destroy_sim(outputs);
        let ctrl_d = &detached[0];
        let x_d = &detached[1..1 + n];
        let got: u64 = (0..n).map(|i| (sim.qubit_mask(&x_d[i]) & 1) << i).sum();
        let mask = (1u64 << n) - 1;
        let expected = if ctrl_bit == 1 {
            x_init.wrapping_add(c_val) & mask
        } else {
            x_init & mask
        };
        assert_eq!(
            got, expected,
            "ctrl_cqadd n={} ctrl={} x={:b} c={:b}: got={:b} exp={:b}",
            n, ctrl_bit, x_init, c_val, got, expected
        );
        assert_eq!(
            sim.qubit_mask(ctrl_d) & 1,
            ctrl_bit,
            "ctrl mutated n={} ctrl={}",
            n,
            ctrl_bit
        );
        assert_eq!(
            sim.phase_mask(),
            0,
            "ctrl_cqadd phase n={} ctrl={} x={} c={}: {:#x}",
            n,
            ctrl_bit,
            x_init,
            c_val,
            sim.phase_mask()
        );
    }

    #[test]
    fn ctrl_cqadd_n4_all() {
        for ctrl in 0..2 {
            for x in 0..16 {
                for c in 0..16 {
                    run_ctrl_cqadd_case(4, ctrl, x, c);
                }
            }
        }
    }

    #[test]
    fn ctrl_cqadd_n8_sample() {
        let mix = |s: u64| {
            s.wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407)
        };
        for seed in 0..200u64 {
            let r = mix(seed);
            run_ctrl_cqadd_case(8, r & 1, (r >> 1) & 0xFF, (r >> 9) & 0xFF);
        }
    }

    #[test]
    fn cqadd_n16_random() {
        // 2^32 exhaustive is too slow; sample 4096 pairs.
        let mix = |s: u64| {
            s.wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407)
        };
        for seed in 0..4096u64 {
            let r = mix(seed);
            let x = r & 0xFFFF;
            let c = (r >> 32) & 0xFFFF;
            run_cqadd_case(16, x, c);
        }
    }

    // Sanity: dec_lemma8_ctrl is the exact inverse of inc_lemma8_ctrl.
    // For every (c, data, prom) classical state, running INC then DEC
    // must return data, c, prom to their starting values.

    // Boundary-heavy cases: exercise α = all 1s (triggers β carry, the
    // hard case of the Eq. 44 derivation), data near wrap, random psi.

    // u128 variant for n > 64.

    // Big-integer variant for n > 128. data represented as a bit-vec.

    fn run_add_cuccaro_case(n: usize, a_init: u64, b_init: u64) {
        let mut circ = Circuit::new();
        let a: Vec<QReg> = (0..n)
            .map(|i| circ.alloc_qreg(&format!("a{}", i)))
            .collect();
        let b: Vec<QReg> = (0..n)
            .map(|i| circ.alloc_qreg(&format!("b{}", i)))
            .collect();
        {
            let mut a_bytes = vec![0u8; n.div_ceil(8)];
            let mut b_bytes = vec![0u8; n.div_ceil(8)];
            for i in 0..n {
                if (a_init >> i) & 1 == 1 {
                    a_bytes[i / 8] |= 1u8 << (i % 8);
                }
                if (b_init >> i) & 1 == 1 {
                    b_bytes[i / 8] |= 1u8 << (i % 8);
                }
            }
            circ.sim_load_reg_bytes_shot(&a, &a_bytes, 0);
            circ.sim_load_reg_bytes_shot(&b, &b_bytes, 0);
        }
        add_cuccaro(&mut circ, &a, &b);
        let mut outputs: Vec<crate::point_add::trailmix_port::circuit::QReg> = Vec::new();
        outputs.extend(a);
        outputs.extend(b);
        let (sim, detached) = circ.destroy_sim(outputs);
        let a_d = &detached[..n];
        let b_d = &detached[n..2 * n];
        let got_a: u64 = (0..n).map(|i| (sim.qubit_mask(&a_d[i]) & 1) << i).sum();
        let got_b: u64 = (0..n).map(|i| (sim.qubit_mask(&b_d[i]) & 1) << i).sum();
        let mask = if n == 64 { u64::MAX } else { (1u64 << n) - 1 };
        let expected_a = a_init.wrapping_add(b_init) & mask;
        assert_eq!(
            got_a,
            expected_a,
            "add_cuccaro n={} a={:0w$b} b={:0w$b}: got_a={:0w$b} exp={:0w$b}",
            n,
            a_init,
            b_init,
            got_a,
            expected_a,
            w = n
        );
        assert_eq!(
            got_b,
            b_init & mask,
            "add_cuccaro b drift n={}: got_b={} exp={}",
            n,
            got_b,
            b_init & mask
        );
        assert_eq!(sim.phase_mask(), 0, "add_cuccaro phase n={}", n);
    }

    #[test]
    fn add_cuccaro_n2_all() {
        for a in 0..4 {
            for b in 0..4 {
                run_add_cuccaro_case(2, a, b);
            }
        }
    }
    #[test]
    fn add_cuccaro_n4_all() {
        for a in 0..16 {
            for b in 0..16 {
                run_add_cuccaro_case(4, a, b);
            }
        }
    }
    #[test]
    fn add_cuccaro_n8_all() {
        for a in 0..256 {
            for b in 0..256 {
                run_add_cuccaro_case(8, a, b);
            }
        }
    }
    #[test]
    fn add_cuccaro_n16_sample() {
        let mix = |s: u64| {
            s.wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407)
        };
        for seed in 0..1024u64 {
            let r = mix(seed);
            run_add_cuccaro_case(16, r & 0xFFFF, (r >> 16) & 0xFFFF);
        }
    }
    #[test]
    fn add_cuccaro_n32_sample() {
        let mix = |s: u64| {
            s.wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407)
        };
        for seed in 0..256u64 {
            let r = mix(seed);
            run_add_cuccaro_case(32, r & 0xFFFFFFFF, (r >> 32) & 0xFFFFFFFF);
        }
    }

    fn run_ctrl_add_cuccaro_ovf_case(n: usize, ctrl_val: u64, a_init: u64, b_init: u64) {
        let mut circ = Circuit::new();
        let ctrl = circ.alloc_qreg("ctrl");
        let a_ext: Vec<QReg> = (0..=n)
            .map(|i| circ.alloc_qreg(&format!("a{}", i)))
            .collect();
        let b: Vec<QReg> = (0..n)
            .map(|i| circ.alloc_qreg(&format!("b{}", i)))
            .collect();
        circ.sim_load_reg_bytes_shot(std::slice::from_ref(&ctrl), &[ctrl_val as u8], 0);
        {
            // Load only the lower n qubits of a_ext; a_ext[n] starts at 0.
            let mut a_bytes = vec![0u8; n.div_ceil(8)];
            let mut b_bytes = vec![0u8; n.div_ceil(8)];
            for i in 0..n {
                if (a_init >> i) & 1 == 1 {
                    a_bytes[i / 8] |= 1u8 << (i % 8);
                }
                if (b_init >> i) & 1 == 1 {
                    b_bytes[i / 8] |= 1u8 << (i % 8);
                }
            }
            circ.sim_load_reg_bytes_shot(&a_ext[..n], &a_bytes, 0);
            circ.sim_load_reg_bytes_shot(&b, &b_bytes, 0);
        }
        controlled_add_cuccaro_with_overflow(&mut circ, &ctrl, &a_ext, &b);
        let mut outputs: Vec<crate::point_add::trailmix_port::circuit::QReg> = Vec::new();
        outputs.push(ctrl);
        outputs.extend(a_ext);
        outputs.extend(b);
        let (sim, detached) = circ.destroy_sim(outputs);
        let ctrl_d = &detached[0];
        // a_ext has n+1 elements (indices 0..=n)
        let a_ext_d = &detached[1..1 + n + 1];
        let b_d = &detached[1 + n + 1..1 + n + 1 + n];
        let got_sum: u64 = (0..n).map(|i| (sim.qubit_mask(&a_ext_d[i]) & 1) << i).sum();
        let got_ovf = sim.qubit_mask(&a_ext_d[n]) & 1;
        let got_b: u64 = (0..n).map(|i| (sim.qubit_mask(&b_d[i]) & 1) << i).sum();
        let got_ctrl = sim.qubit_mask(ctrl_d) & 1;
        let mask = if n == 64 { u64::MAX } else { (1u64 << n) - 1 };
        let (expected_sum, expected_ovf) = if ctrl_val == 1 {
            let full = a_init.wrapping_add(b_init);
            (full & mask, if n == 64 { 0 } else { (full >> n) & 1 })
        } else {
            (a_init & mask, 0)
        };
        assert_eq!(
            got_sum, expected_sum,
            "n={} ctrl={} a={} b={}",
            n, ctrl_val, a_init, b_init
        );
        assert_eq!(got_ovf, expected_ovf, "ovf n={} ctrl={}", n, ctrl_val);
        assert_eq!(got_b, b_init & mask, "b drift");
        assert_eq!(got_ctrl, ctrl_val, "ctrl drift");
        assert_eq!(sim.phase_mask(), 0, "phase");
    }

    #[test]
    fn ctrl_add_cuccaro_ovf_n3_all() {
        for c in 0..2 {
            for a in 0..8 {
                for b in 0..8 {
                    run_ctrl_add_cuccaro_ovf_case(3, c, a, b);
                }
            }
        }
    }

    #[test]
    fn ctrl_add_cuccaro_ovf_n4_all() {
        for c in 0..2 {
            for a in 0..16 {
                for b in 0..16 {
                    run_ctrl_add_cuccaro_ovf_case(4, c, a, b);
                }
            }
        }
    }

    #[test]
    fn ctrl_add_cuccaro_ovf_n8_sample() {
        let mix = |s: u64| {
            s.wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407)
        };
        for seed in 0..256u64 {
            let r = mix(seed);
            run_ctrl_add_cuccaro_ovf_case(8, r & 1, (r >> 1) & 0xFF, (r >> 9) & 0xFF);
        }
    }

    fn run_add_cuccaro_overflow_case(n: usize, a_init: u64, b_init: u64) {
        let mut circ = Circuit::new();
        let a_ext: Vec<QReg> = (0..=n)
            .map(|i| circ.alloc_qreg(&format!("a{}", i)))
            .collect();
        let b: Vec<QReg> = (0..n)
            .map(|i| circ.alloc_qreg(&format!("b{}", i)))
            .collect();
        {
            // Load only the lower n qubits of a_ext; a_ext[n] starts at 0.
            let mut a_bytes = vec![0u8; n.div_ceil(8)];
            let mut b_bytes = vec![0u8; n.div_ceil(8)];
            for i in 0..n {
                if (a_init >> i) & 1 == 1 {
                    a_bytes[i / 8] |= 1u8 << (i % 8);
                }
                if (b_init >> i) & 1 == 1 {
                    b_bytes[i / 8] |= 1u8 << (i % 8);
                }
            }
            circ.sim_load_reg_bytes_shot(&a_ext[..n], &a_bytes, 0);
            circ.sim_load_reg_bytes_shot(&b, &b_bytes, 0);
        }
        // a_ext[n] starts at 0.
        add_cuccaro_with_overflow(&mut circ, &a_ext, &b);
        let mut outputs: Vec<crate::point_add::trailmix_port::circuit::QReg> = Vec::new();
        outputs.extend(a_ext);
        outputs.extend(b);
        let (sim, detached) = circ.destroy_sim(outputs);
        // a_ext has n+1 elements
        let a_ext_d = &detached[..n + 1];
        let b_d = &detached[n + 1..n + 1 + n];
        let got_sum: u64 = (0..n).map(|i| (sim.qubit_mask(&a_ext_d[i]) & 1) << i).sum();
        let got_ovf = sim.qubit_mask(&a_ext_d[n]) & 1;
        let got_b: u64 = (0..n).map(|i| (sim.qubit_mask(&b_d[i]) & 1) << i).sum();
        let full = a_init.wrapping_add(b_init);
        let mask = if n == 64 { u64::MAX } else { (1u64 << n) - 1 };
        let expected_sum = full & mask;
        let expected_ovf = if n == 64 { 0 } else { (full >> n) & 1 };
        assert_eq!(
            got_sum, expected_sum,
            "add_cuccaro_ovf sum n={} a={:b} b={:b}: got={:b} exp={:b}",
            n, a_init, b_init, got_sum, expected_sum
        );
        assert_eq!(
            got_ovf, expected_ovf,
            "add_cuccaro_ovf n={} a={} b={}: ovf got={} exp={}",
            n, a_init, b_init, got_ovf, expected_ovf
        );
        assert_eq!(got_b, b_init & mask, "b drift n={}", n);
        assert_eq!(sim.phase_mask(), 0, "phase n={}", n);
    }

    #[test]
    fn add_cuccaro_ovf_n4_all() {
        for a in 0..16 {
            for b in 0..16 {
                run_add_cuccaro_overflow_case(4, a, b);
            }
        }
    }
    #[test]
    fn add_cuccaro_ovf_n8_all() {
        for a in 0..256 {
            for b in 0..256 {
                run_add_cuccaro_overflow_case(8, a, b);
            }
        }
    }
    // [DELETED 2026-05-30] tests for `controlled_add_cuccaro` (10n CCX
    // variant) removed alongside the primitive itself. Coverage is now
    // provided by the `ctrl_add_cuccaro_3n_*` test family below.

    fn run_ctrl_add_cuccaro_3n_case(n: usize, ctrl_val: u64, a_init: u64, b_init: u64) {
        let mut circ = Circuit::new();
        let ctrl = circ.alloc_qreg("ctrl");
        let a: Vec<QReg> = (0..n)
            .map(|i| circ.alloc_qreg(&format!("a{}", i)))
            .collect();
        let b: Vec<QReg> = (0..n)
            .map(|i| circ.alloc_qreg(&format!("b{}", i)))
            .collect();
        circ.sim_load_reg_bytes_shot(std::slice::from_ref(&ctrl), &[ctrl_val as u8], 0);
        {
            let mut a_bytes = vec![0u8; n.div_ceil(8)];
            let mut b_bytes = vec![0u8; n.div_ceil(8)];
            for i in 0..n {
                if (a_init >> i) & 1 == 1 {
                    a_bytes[i / 8] |= 1u8 << (i % 8);
                }
                if (b_init >> i) & 1 == 1 {
                    b_bytes[i / 8] |= 1u8 << (i % 8);
                }
            }
            circ.sim_load_reg_bytes_shot(&a, &a_bytes, 0);
            circ.sim_load_reg_bytes_shot(&b, &b_bytes, 0);
        }
        crate::point_add::trailmix_port::arith::cuccaro::controlled_add_cuccaro_3n(&mut circ, &ctrl, &a, &b);
        let mut outputs: Vec<crate::point_add::trailmix_port::circuit::QReg> = Vec::new();
        outputs.push(ctrl);
        outputs.extend(a);
        outputs.extend(b);
        let (sim, detached) = circ.destroy_sim(outputs);
        let ctrl_d = &detached[0];
        let a_d = &detached[1..1 + n];
        let b_d = &detached[1 + n..1 + 2 * n];
        let got_a: u64 = (0..n).map(|i| (sim.qubit_mask(&a_d[i]) & 1) << i).sum();
        let got_b: u64 = (0..n).map(|i| (sim.qubit_mask(&b_d[i]) & 1) << i).sum();
        let got_ctrl = sim.qubit_mask(ctrl_d) & 1;
        let mask = if n == 64 { u64::MAX } else { (1u64 << n) - 1 };
        let expected_a = if ctrl_val == 1 {
            a_init.wrapping_add(b_init) & mask
        } else {
            a_init & mask
        };
        assert_eq!(
            got_a, expected_a,
            "ctrl_add_cuccaro_3n n={n} ctrl={ctrl_val} a={a_init:x} b={b_init:x}: got {got_a:x}, exp {expected_a:x}"
        );
        assert_eq!(got_b, b_init & mask, "ctrl_add_cuccaro_3n: b drift n={n}");
        assert_eq!(got_ctrl, ctrl_val, "ctrl_add_cuccaro_3n: ctrl drift n={n}");
        assert_eq!(sim.phase_mask(), 0, "ctrl_add_cuccaro_3n: phase n={n}");
    }

    #[test]
    fn ctrl_add_cuccaro_3n_n3_all() {
        for c in 0..2 {
            for a in 0..8 {
                for b in 0..8 {
                    run_ctrl_add_cuccaro_3n_case(3, c, a, b);
                }
            }
        }
    }

    #[test]
    fn ctrl_add_cuccaro_3n_n4_all() {
        for c in 0..2 {
            for a in 0..16 {
                for b in 0..16 {
                    run_ctrl_add_cuccaro_3n_case(4, c, a, b);
                }
            }
        }
    }

    #[test]
    fn ctrl_add_cuccaro_3n_n8_sample() {
        let mix = |s: u64| {
            s.wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407)
        };
        for seed in 0..512u64 {
            let r = mix(seed);
            run_ctrl_add_cuccaro_3n_case(8, r & 1, (r >> 1) & 0xFF, (r >> 9) & 0xFF);
        }
    }

    #[test]
    fn ctrl_add_cuccaro_3n_n16_sample() {
        let mix = |s: u64| {
            s.wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407)
        };
        for seed in 0..256u64 {
            let r = mix(seed);
            run_ctrl_add_cuccaro_3n_case(16, r & 1, (r >> 1) & 0xFFFF, (r >> 17) & 0xFFFF);
        }
    }

    /// Gate-count check: controlled_add_cuccaro_3n should use 3n - 2 CCX
    /// for n >= 2 (n CCX forward MAJ + 2n CCX reverse, less 2 for the
    /// truncated boundary).
    ///
    /// Vs `controlled_add_cuccaro_mbu` (8n CCX) this is a ~2.6x reduction
    /// per call -- the dominant savings for mod_mul's controlled-add path.
    #[test]
    fn ctrl_add_cuccaro_3n_tof_count_n32() {
        let n = 32;
        let mut circ = Circuit::new();
        let ctrl = circ.alloc_qreg("ctrl");
        let a: Vec<QReg> = (0..n).map(|i| circ.alloc_qreg(&format!("a{i}"))).collect();
        let b: Vec<QReg> = (0..n).map(|i| circ.alloc_qreg(&format!("b{i}"))).collect();
        let tof_before = circ.ccx_emitted;
        crate::point_add::trailmix_port::arith::cuccaro::controlled_add_cuccaro_3n(&mut circ, &ctrl, &a, &b);
        let tof = (circ.ccx_emitted - tof_before) as usize;
        let expected = 3 * n - 2;
        assert_eq!(
            tof, expected,
            "controlled_add_cuccaro_3n n={n}: expected {expected} CCX, got {tof}"
        );
    }

    #[test]
    fn add_cuccaro_ovf_n16_sample() {
        let mix = |s: u64| {
            s.wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407)
        };
        for seed in 0..256u64 {
            let r = mix(seed);
            run_add_cuccaro_overflow_case(16, r & 0xFFFF, (r >> 16) & 0xFFFF);
        }
    }

    // [DELETED 2026-05-30] Tests for `controlled_add_cuccaro_mbu` and
    // `controlled_add_cuccaro_mbu_refs` (8n CCX streaming-MBU variant)
    // have been removed alongside the primitives themselves. The
    // `ctrl_add_cuccaro_3n_*` family below covers all behaviours,
    // including the reversed-physical refs-view case (which the 3n
    // primitive also supports via `controlled_add_cuccaro_3n_refs`).

    fn run_compare_geq_t3_case(n: usize, x_init: u64, c_val: u64) {
        let mut circ = Circuit::new();
        let x: Vec<QReg> = (0..n)
            .map(|i| circ.alloc_qreg(&format!("x{}", i)))
            .collect();
        let out = circ.alloc_qreg("out");
        {
            let mut bytes = vec![0u8; n.div_ceil(8)];
            for i in 0..n {
                if (x_init >> i) & 1 == 1 {
                    bytes[i / 8] |= 1u8 << (i % 8);
                }
            }
            circ.sim_load_reg_bytes_shot(&x, &bytes, 0);
        }
        let n_bytes = n.div_ceil(8).max(1);
        let mut c_bytes = vec![0u8; n_bytes];
        for i in 0..64 {
            if (c_val >> i) & 1 == 1 && i / 8 < n_bytes {
                c_bytes[i / 8] |= 1 << (i % 8);
            }
        }

        compare_geq_theorem3(&mut circ, &x, &c_bytes, &out);

        let mut outputs: Vec<crate::point_add::trailmix_port::circuit::QReg> = Vec::new();
        outputs.extend(x);
        outputs.push(out);
        let (sim, detached) = circ.destroy_sim(outputs);
        let x_d = &detached[..n];
        let out_d = &detached[n];
        let got = sim.qubit_mask(out_d) & 1;
        let mask = if n == 64 { u64::MAX } else { (1u64 << n) - 1 };
        let x_masked = x_init & mask;
        let c_masked = c_val & mask;
        let expected = if x_masked >= c_masked { 1 } else { 0 };
        assert_eq!(
            got,
            expected,
            "compare_t3 n={} x={:0w$b} c={:0w$b}: got={} exp={}",
            n,
            x_masked,
            c_masked,
            got,
            expected,
            w = n
        );
        // x must be preserved.
        for i in 0..n {
            let got_xi = (sim.qubit_mask(&x_d[i]) & 1) == 1;
            let exp_xi = (x_init >> i) & 1 == 1;
            assert_eq!(
                got_xi, exp_xi,
                "compare_t3 x drift n={} bit {}: got={} exp={}",
                n, i, got_xi, exp_xi
            );
        }
        assert_eq!(
            sim.phase_mask(),
            0,
            "compare_t3 phase n={} x={} c={}",
            n,
            x_init,
            c_val
        );
    }

    #[test]
    fn compare_t3_n1_all() {
        for x in 0..2 {
            for c in 0..2 {
                run_compare_geq_t3_case(1, x, c);
            }
        }
    }
    #[test]
    fn compare_t3_n2_all() {
        for x in 0..4 {
            for c in 0..4 {
                run_compare_geq_t3_case(2, x, c);
            }
        }
    }
    #[test]
    fn compare_t3_n3_all() {
        for x in 0..8 {
            for c in 0..8 {
                run_compare_geq_t3_case(3, x, c);
            }
        }
    }
    #[test]
    fn compare_t3_n4_all() {
        for x in 0..16 {
            for c in 0..16 {
                run_compare_geq_t3_case(4, x, c);
            }
        }
    }
    #[test]
    fn compare_t3_n8_all() {
        for x in 0..256 {
            for c in 0..256 {
                run_compare_geq_t3_case(8, x, c);
            }
        }
    }

    #[test]
    fn compare_t3_n16_sample() {
        let mix = |s: u64| {
            s.wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407)
        };
        for seed in 0..1024u64 {
            let r = mix(seed);
            let x = r & 0xFFFF;
            let c = (r >> 16) & 0xFFFF;
            run_compare_geq_t3_case(16, x, c);
        }
    }

    #[test]
    fn compare_t3_n16_boundary() {
        // Hardest: adversarial 0xAA pattern (alternation).
        for x in 0..16u64 {
            run_compare_geq_t3_case(16, x * 0x1111, 0xAAAA);
            run_compare_geq_t3_case(16, x * 0x1111, 0x5555);
        }
        // All bits
        run_compare_geq_t3_case(16, 0, 0);
        run_compare_geq_t3_case(16, 0xFFFF, 0xFFFF);
        run_compare_geq_t3_case(16, 0xFFFF, 0);
        run_compare_geq_t3_case(16, 0, 0xFFFF);
    }

    #[test]
    fn compare_t3_n32_sample() {
        let mix = |s: u64| {
            s.wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407)
        };
        for seed in 0..256u64 {
            let r = mix(seed);
            let x = r & 0xFFFFFFFF;
            let c = (r >> 32) & 0xFFFFFFFF;
            run_compare_geq_t3_case(32, x, c);
        }
    }
}
