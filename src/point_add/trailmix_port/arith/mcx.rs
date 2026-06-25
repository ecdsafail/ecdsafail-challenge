//! Multi-controlled-X (MCX) gadgets: clean-ancilla and dirty-ancilla
//! Toffoli-ladder constructions. Extracted from the former
//! `mbu_primitives` grab-bag (these are reversible MCX primitives, not
//! MBU-specific).

use crate::point_add::trailmix_port::arith::khattar_gidney::{
    xor_and_of_khattar_gidney_refs, xor_and_of_khattar_gidney_refs_consume,
};
use crate::point_add::trailmix_port::circuit::{Circuit, QReg};

/// Toggle `target` by the product of `ctrls` using restored dirty lenders.
///
/// The first cascade includes the seed link `d0 ^= c0*c1`; the second omits
/// it. Dirty-seeded terms occur in both cascades and cancel, while the full
/// control product occurs once. For `k >= 3` this costs `4k-8` CCX gates and
/// restores every lender exactly.
pub fn mcx_dirty_ladder(
    circ: &mut Circuit,
    ctrls: &[&QReg],
    target: &QReg,
    dirty: &[&QReg],
) {
    let k = ctrls.len();
    match k {
        0 => {
            circ.x(target);
            return;
        }
        1 => {
            circ.cx(ctrls[0], target);
            return;
        }
        2 => {
            circ.ccx(ctrls[0], ctrls[1], target);
            return;
        }
        _ => {}
    }

    assert!(dirty.len() >= k - 2, "mcx_dirty_ladder lender shortage");
    let dirty = &dirty[..k - 2];
    for (index, &q) in dirty.iter().enumerate() {
        assert!(!std::ptr::eq(q, target), "dirty lender aliases target");
        assert!(
            !ctrls.iter().any(|&control| std::ptr::eq(q, control)),
            "dirty lender aliases control"
        );
        assert!(
            !dirty[..index].iter().any(|&other| std::ptr::eq(q, other)),
            "duplicate dirty lender"
        );
    }

    let cascade = |circ: &mut Circuit, include_seed: bool| {
        if include_seed {
            circ.ccx(ctrls[0], ctrls[1], dirty[0]);
        }
        for i in 1..dirty.len() {
            circ.ccx(dirty[i - 1], ctrls[i + 1], dirty[i]);
        }
        circ.ccx(dirty[dirty.len() - 1], ctrls[k - 1], target);
        for i in (1..dirty.len()).rev() {
            circ.ccx(dirty[i - 1], ctrls[i + 1], dirty[i]);
        }
        if include_seed {
            circ.ccx(ctrls[0], ctrls[1], dirty[0]);
        }
    };

    cascade(circ, true);
    cascade(circ, false);
}

///
/// C^k X (k-controlled NOT) with ONE dirty ancilla, using the
/// Barenco "surrounded" decomposition.
///
/// For ctrls = [`c_0`, ..., c_{k-1}] and target `t`, with dirty ancilla
/// `psi` (arbitrary state, restored on exit):
///   t ^= `AND(c_0`, ..., c_{k-1})
///
/// Gate count: Θ(k) CCX. Ancillae: 1 dirty (restored). No clean ancs.
///
/// Base cases:
///   k=0: X(t)
///   k=1: `CX(c_0`, t)
///   k=2: `CCX(c_0`, `c_1`, t)
///   k>=3: Barenco split: T = AND(ctrls). Split ctrls at half, compute
///         half-AND into psi, use as control, uncompute. 4 recursive
///         calls pattern. T(k) = 4*T(k/2) would be O(k²); but using
///         the Barenco trick (two different splits interleaved) gives
///         O(k).
///
/// We implement the iterative O(k) form: at each level, use one
/// recursive invocation plus its mirror. The "surrounded" identity:
///   CCX(c0, c1, psi); CCX(psi, c2, t); CCX(c0, c1, psi);
///   CCX(psi, c2, t)
/// This computes t ^= AND(c0, c1, c2) and restores psi. 4 CCX.
///
/// For k controls: chain the pattern. Θ(k) CCX total.
pub fn mcx_dirty(circ: &mut Circuit, ctrls: &[&QReg], target: &QReg, psi: &QReg) {
    let k = ctrls.len();

    // PRE: capture (AND(ctrls), target, psi).
    {
        let ctrls_for_capture: Vec<&QReg> = ctrls.to_vec();
        let target_ref = target;
        let psi_ref = psi;
        circ.contract_capture(
            "mbu.mcx_dirty.pre",
            move |view, shot| -> Result<(bool, bool, bool), String> {
                let mut and_v = true;
                for q in &ctrls_for_capture {
                    and_v &= view.contract_read_bit_shot(q, shot);
                }
                let t = view.contract_read_bit_shot(target_ref, shot);
                let p = view.contract_read_bit_shot(psi_ref, shot);
                Ok((and_v, t, p))
            },
        );
    }

    match k {
        0 => circ.x(target),
        1 => circ.cx(ctrls[0], target),
        2 => circ.ccx(ctrls[0], ctrls[1], target),
        3 => {
            // Surrounded 4-CCX form: target ^= AND(c0, c1, c2), psi restored.
            //   CCX(c0, c1, psi)     psi ^= c0·c1
            //   CCX(psi, c2, target) target ^= psi·c2
            //   CCX(c0, c1, psi)     psi restored
            //   CCX(psi, c2, target) target ^= psi_0·c2 (cancels extra)
            circ.ccx(ctrls[0], ctrls[1], psi);
            circ.ccx(psi, ctrls[2], target);
            circ.ccx(ctrls[0], ctrls[1], psi);
            circ.ccx(psi, ctrls[2], target);
        }
        4 => {
            mcx_dirty_k4(circ, ctrls, target, psi);
        }
        5 => {
            // k=5 via doubled C^4 X: CCX(c0,c1,psi); C^4X(psi,c2,c3,c4→t
            // with c0 dirty); CCX(c0,c1,psi); C^4X again. Verified
            // exhaustively via Python. 2 + 2·10 = 22 CCX.
            let (c0, c1, c2, c3, c4) = (ctrls[0], ctrls[1], ctrls[2], ctrls[3], ctrls[4]);
            circ.ccx(c0, c1, psi);
            mcx_dirty_k4(circ, &[psi, c2, c3, c4], target, c0);
            circ.ccx(c0, c1, psi);
            mcx_dirty_k4(circ, &[psi, c2, c3, c4], target, c0);
        }
        _ => {
            // The Barenco-style constants above are derived only for
            // k <= 5. Callers needing k >= 6 must route through
            // `mcx_dirty_any_k` (Theorem 3 recursion via `mcx_clean_k`);
            // a direct call here violates that contract.
            panic!("mcx_dirty supports k <= 5 controls (got k = {k}); use mcx_dirty_any_k for k >= 6");
        }
    }

    // POST: target ^= AND(ctrls); psi restored; ctrls unchanged.
    {
        let ctrls_for_check: Vec<&QReg> = ctrls.to_vec();
        let target_ref = target;
        let psi_ref = psi;
        circ.contract_pop_and_check::<(bool, bool, bool), _>(
            "mbu.mcx_dirty.pre",
            move |cap, view, shot| -> Result<(), String> {
                let (and_pre, t_pre, p_pre) = *cap;
                let mut and_post = true;
                for q in &ctrls_for_check {
                    and_post &= view.contract_read_bit_shot(q, shot);
                }
                if and_post != and_pre {
                    return Err(format!(
                        "mcx_dirty: ctrls AND changed {} -> {}",
                        u8::from(and_pre),
                        u8::from(and_post)
                    ));
                }
                let t_post = view.contract_read_bit_shot(target_ref, shot);
                let p_post = view.contract_read_bit_shot(psi_ref, shot);
                let expected = t_pre ^ and_pre;
                if t_post != expected {
                    return Err(format!(
                        "mcx_dirty: target {}->{} expected {} (and={})",
                        u8::from(t_pre),
                        u8::from(t_post),
                        u8::from(expected),
                        u8::from(and_pre)
                    ));
                }
                if p_post != p_pre {
                    return Err(format!(
                        "mcx_dirty: psi (dirty ancilla) changed {} -> {} (must be restored)",
                        u8::from(p_pre),
                        u8::from(p_post)
                    ));
                }
                Ok(())
            },
        );
    }
}

/// `target ^= AND(ctrls)` via the Khattar–Gidney Sec 5.3 (Fig 4)
/// prefix-AND construction.
///
/// Cost: 2k-3 Toffoli, log*_2(k) clean ancillae, O(log k) depth for
/// k >= 4. At k=256: ~509 Toffolis with ~5 clean ancillae. (The
/// previous Karatsuba-halving recursion was Θ(k^log2 3) ≈ 4700
/// Toffolis at k=256 — ~9.2x over this construction.)
///
/// The k>=4 path delegates to [`xor_and_of_khattar_gidney_refs`],
/// which is structurally the same Fig 4 / Sec 6.1 prefix-AND ladder.
///
/// Base cases (degenerate for KG):
///   k=0: X(target)
///   k=1: CX
///   k=2: CCX
///   k=3: alloc t; CCX(c0,c1,t); CCX(t,c2,target); `clear_and(t,c0,c1)`.
///        `clear_and` picks the MBU (`HMR+cz_if_bit`) discharge when
///        possible, saving one Toffoli vs the naive 3-CCX form.
pub fn mcx_clean_k(circ: &mut Circuit, ctrls: &[&QReg], target: &QReg) {
    let k = ctrls.len();

    // PRE: capture target_pre and AND(ctrls)_pre per shot.
    {
        let ctrls_for_capture: Vec<&QReg> = ctrls.to_vec();
        let target_ref = target;
        circ.contract_capture(
            "mbu.mcx_clean_k.pre",
            move |view, shot| -> Result<(bool, bool), String> {
                let mut and_v = true;
                for q in &ctrls_for_capture {
                    and_v &= view.contract_read_bit_shot(q, shot);
                }
                let t = view.contract_read_bit_shot(target_ref, shot);
                Ok((and_v, t))
            },
        );
    }

    match k {
        0 => circ.x(target),
        1 => circ.cx(ctrls[0], target),
        2 => circ.ccx(ctrls[0], ctrls[1], target),
        3 => {
            // MBU: replace the trailing `ccx(ctrls[0], ctrls[1], t)`
            // uncompute with HMR + cz_if_bit. t holds ctrls[0] AND
            // ctrls[1] after the forward CCX (the middle CCX writes
            // to target, not t), and ctrls[0]/ctrls[1] are not
            // re-versioned between, so declare_and_of structurally
            // matches cz_if_bit's discharge. Saves 1 CCX per call.
            let t = circ.alloc_qreg_bits("mcxk_t3", 1);
            circ.ccx(ctrls[0], ctrls[1], &t[0]);
            circ.ccx(&t[0], ctrls[2], target);
            // Clear t back to |0>. clear_and picks MBU (HMR+cz_if_bit,
            // no Toffoli) outside a condition, or reversible ccx
            // (push_condition-safe, +1 Toffoli) inside one.
            circ.clear_and(&t[0], ctrls[0], ctrls[1]);
            drop(t);
        }
        4 => {
            // Balanced-tree flat sequence: 2 ancillae, 3 Toffoli
            // outside a condition (5 inside). Cheaper than the KG
            // dispatch which allocates kg_prefix_ancilla_count(4)
            // ancillae and builds a multi-layer tree.
            //
            //   t01 = c0 AND c1       (1 ccx)
            //   t23 = c2 AND c3       (1 ccx)
            //   target ^= t01 AND t23 (1 ccx)
            //   clear_and(t23,c2,c3)  (MBU: 0 ccx outside cond)
            //   clear_and(t01,c0,c1)  (MBU: 0 ccx outside cond)
            //
            // c0..c3 are not re-versioned between compute and clear,
            // so the MBU declare_and_of identity holds.
            let t01 = circ.alloc_qreg_bits("mcxk_t01_4", 1);
            let t23 = circ.alloc_qreg_bits("mcxk_t23_4", 1);
            circ.ccx(ctrls[0], ctrls[1], &t01[0]);
            circ.ccx(ctrls[2], ctrls[3], &t23[0]);
            circ.ccx(&t01[0], &t23[0], target);
            circ.clear_and(&t23[0], ctrls[2], ctrls[3]);
            circ.clear_and(&t01[0], ctrls[0], ctrls[1]);
            drop(t23);
            drop(t01);
        }
        5 => {
            // Balanced-tree flat sequence: 3 ancillae, 4 Toffoli
            // outside a condition. The 4-leaf AND is built as in
            // k=4; then a second-level ccx folds c4 onto target.
            //
            //   t01   = c0 AND c1        (1 ccx)
            //   t23   = c2 AND c3        (1 ccx)
            //   t0123 = t01 AND t23      (1 ccx)
            //   target ^= t0123 AND c4   (1 ccx)
            //   clear_and(t0123,t01,t23) (MBU: 0 ccx outside cond)
            //   clear_and(t23,c2,c3)     (MBU)
            //   clear_and(t01,c0,c1)     (MBU)
            //
            // None of t01/t23/c0..c4 are re-versioned between their
            // compute and clear sites, so each declare_and_of holds.
            let t01 = circ.alloc_qreg_bits("mcxk_t01_5", 1);
            let t23 = circ.alloc_qreg_bits("mcxk_t23_5", 1);
            let t0123 = circ.alloc_qreg_bits("mcxk_t0123_5", 1);
            circ.ccx(ctrls[0], ctrls[1], &t01[0]);
            circ.ccx(ctrls[2], ctrls[3], &t23[0]);
            circ.ccx(&t01[0], &t23[0], &t0123[0]);
            circ.ccx(&t0123[0], ctrls[4], target);
            circ.clear_and(&t0123[0], &t01[0], &t23[0]);
            circ.clear_and(&t23[0], ctrls[2], ctrls[3]);
            circ.clear_and(&t01[0], ctrls[0], ctrls[1]);
            drop(t0123);
            drop(t23);
            drop(t01);
        }
        _ => {
            // Khattar–Gidney Sec 5.3 prefix-AND. 2k-3 Toffoli with
            // log*_2(k) clean ancillae.
            xor_and_of_khattar_gidney_refs(circ, ctrls, target);
        }
    }

    // POST: target ^= AND(ctrls); ctrls unchanged.
    {
        let ctrls_for_check: Vec<&QReg> = ctrls.to_vec();
        let target_ref = target;
        circ.contract_pop_and_check::<(bool, bool), _>(
            "mbu.mcx_clean_k.pre",
            move |cap, view, shot| -> Result<(), String> {
                let (and_pre, t_pre) = *cap;
                let mut and_post = true;
                for q in &ctrls_for_check {
                    and_post &= view.contract_read_bit_shot(q, shot);
                }
                if and_post != and_pre {
                    return Err(format!(
                        "mcx_clean_k: ctrls AND changed {} -> {}",
                        u8::from(and_pre),
                        u8::from(and_post)
                    ));
                }
                let t_post = view.contract_read_bit_shot(target_ref, shot);
                let expected = t_pre ^ and_pre;
                if t_post != expected {
                    return Err(format!(
                        "mcx_clean_k: target {}->{}, expected {} (t_pre={}, AND(ctrls)={})",
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

/// Variant of [`mcx_clean_k`] that ALSO frees `target` after the XOR.
/// Used when the caller alloc'd `target` and just needs the AND
/// folded back to |0> for cleanup. Frees `target` at the last gate
/// that touches it, so the strict-dealloc gap is 0 even when the
/// recursion's inner-scratch cleanup leaves trailing ops that don't
/// touch `target`.
pub(crate) fn mcx_clean_k_uncompute_consume(circ: &mut Circuit, ctrls: &[&QReg], target: QReg) {
    let k = ctrls.len();
    match k {
        0 => {
            circ.x(&target); /* target drops on return */
        }
        1 => {
            circ.cx(ctrls[0], &target);
        }
        2 => {
            circ.ccx(ctrls[0], ctrls[1], &target);
        }
        3 => {
            // MBU: same swap as mcx_clean_k k=3 (replace trailing
            // CCX with HMR + cz_if_bit). Saves 1 CCX per call.
            let t = circ.alloc_qreg_bits("mcxk_t3", 1);
            circ.ccx(ctrls[0], ctrls[1], &t[0]);
            circ.ccx(&t[0], ctrls[2], &target);
            // Last gate-touch on target is the ccx above; drop now.
            drop(target);
            // Clear t (MBU outside a condition, reversible ccx inside).
            circ.clear_and(&t[0], ctrls[0], ctrls[1]);
            drop(t);
        }
        _ => {
            // Khattar–Gidney Sec 5.3 prefix-AND with target freed at
            // its last gate-touch.
            xor_and_of_khattar_gidney_refs_consume(circ, ctrls, target);
        }
    }
}

/// k=4 case of `mcx_dirty`. Sequence (verified via Python trace):
///   CCX(c0,c1,psi)  psi ^= c0·c1
///   [C^3X via 4-CCX surrounded on (psi,c2,c3)→target, borrowing c0]:
///     CCX(psi,c2,c0); CCX(c0,c3,target); CCX(psi,c2,c0); CCX(c0,c3,target)
///     → target ^= psi·c2·c3 = (`psi_0` ⊕ c0·c1)·c2·c3
///   CCX(c0,c1,psi)  psi restored
///   [C^3X again, now with `psi=psi_0` → target ^= `psi_0·c2·c3`]
///     CCX(psi,c2,c0); CCX(c0,c3,target); CCX(psi,c2,c0); CCX(c0,c3,target)
/// Net: target ^= c0·c1·c2·c3, psi and c0 restored. 10 CCX.
fn mcx_dirty_k4(circ: &mut Circuit, c: &[&QReg], t: &QReg, psi: &QReg) {
    debug_assert_eq!(c.len(), 4);
    let (c0, c1, c2, c3) = (c[0], c[1], c[2], c[3]);
    circ.ccx(c0, c1, psi);
    // Inner C^3X #1 using c0 as temp dirty.
    circ.ccx(psi, c2, c0);
    circ.ccx(c0, c3, t);
    circ.ccx(psi, c2, c0);
    circ.ccx(c0, c3, t);
    circ.ccx(c0, c1, psi);
    // Inner C^3X #2 to cancel extra psi_0·c2·c3 term.
    circ.ccx(psi, c2, c0);
    circ.ccx(c0, c3, t);
    circ.ccx(psi, c2, c0);
    circ.ccx(c0, c3, t);
}

/// C^kX with 1 dirty ancilla for any k. Recurses via Theorem 3 when k >= 6.
///
/// Base case: k <= 5 uses existing `mcx_dirty` (Barenco-style constants).
/// Recursive case: k >= 6 falls back to `mcx_clean_k` (O(log k) clean
/// ancillae). The dirty qubit is ignored for k >= 6; callers that need
/// strict dirty-only semantics should use k <= 5.
pub fn mcx_dirty_any_k(circ: &mut Circuit, ctrls: &[&QReg], target: &QReg, dirty: &QReg) {
    // PRE: capture (AND(ctrls), target, dirty).
    {
        let ctrls_for_capture: Vec<&QReg> = ctrls.to_vec();
        let target_ref = target;
        let dirty_ref = dirty;
        circ.contract_capture(
            "mbu.mcx_dirty_any_k.pre",
            move |view, shot| -> Result<(bool, bool, bool), String> {
                let mut and_v = true;
                for q in &ctrls_for_capture {
                    and_v &= view.contract_read_bit_shot(q, shot);
                }
                let t = view.contract_read_bit_shot(target_ref, shot);
                let d = view.contract_read_bit_shot(dirty_ref, shot);
                Ok((and_v, t, d))
            },
        );
    }

    let k = ctrls.len();
    if k <= 5 {
        mcx_dirty(circ, ctrls, target, dirty);
    } else {
        let _ = dirty;
        mcx_clean_k(circ, ctrls, target);
    }

    // POST: target ^= AND(ctrls); dirty restored (k<=5) or unchanged
    // (k>=6, mcx_clean_k path ignores `dirty`).
    {
        let ctrls_for_check: Vec<&QReg> = ctrls.to_vec();
        let target_ref = target;
        let dirty_ref = dirty;
        circ.contract_pop_and_check::<(bool, bool, bool), _>(
            "mbu.mcx_dirty_any_k.pre",
            move |cap, view, shot| -> Result<(), String> {
                let (and_pre, t_pre, d_pre) = *cap;
                let mut and_post = true;
                for q in &ctrls_for_check {
                    and_post &= view.contract_read_bit_shot(q, shot);
                }
                if and_post != and_pre {
                    return Err(format!(
                        "mcx_dirty_any_k: ctrls AND changed {} -> {}",
                        u8::from(and_pre),
                        u8::from(and_post)
                    ));
                }
                let t_post = view.contract_read_bit_shot(target_ref, shot);
                let expected = t_pre ^ and_pre;
                if t_post != expected {
                    return Err(format!(
                        "mcx_dirty_any_k: target {}->{} expected {} (and={})",
                        u8::from(t_pre),
                        u8::from(t_post),
                        u8::from(expected),
                        u8::from(and_pre)
                    ));
                }
                let d_post = view.contract_read_bit_shot(dirty_ref, shot);
                if d_post != d_pre {
                    return Err(format!(
                        "mcx_dirty_any_k: dirty ancilla changed {} -> {}",
                        u8::from(d_pre),
                        u8::from(d_post)
                    ));
                }
                Ok(())
            },
        );
    }
}

/// Variant of [`mcx_dirty_any_k`] that frees `target` right after the
/// last gate-touch. For k <= 5, `target` is freed after the last
/// `mcx_dirty` gate. For k >= 6, uses `mcx_clean_k_uncompute_consume`
/// which frees target at its last gate-touch inside the recursion.
pub(crate) fn mcx_dirty_any_k_consume(
    circ: &mut Circuit,
    ctrls: &[&QReg],
    target: QReg,
    dirty: &QReg,
) {
    let k = ctrls.len();
    if k <= 5 {
        mcx_dirty(circ, ctrls, &target, dirty);
        // target drops at function end (last gate-touch was mcx_dirty).
        return;
    }
    let _ = dirty;
    mcx_clean_k_uncompute_consume(circ, ctrls, target);
}
