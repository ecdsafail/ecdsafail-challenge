//! The GCD swap-decision comparator for the product-min secp256k1 EC-add,
//! built on this crate's `B` builder.
//!
//! ## What the schedule needs
//! At GCD jump step `i` the swap decision compares the GCD cofactors `u`,`v`
//! over a narrow window, the per-step `GAP_J2[i]`. The comparator scans the
//! top `k` bits of the operands; it mis-decides the `u <-> v` swap iff the
//! highest differing bit of `u`,`v` sits below the window, i.e. it returns
//! "equal -> no swap" when the top-`k` MSBs of the two operands agree. On
//! uniform operands this happens with probability ~2^-k per call.
//!
//! ## Carry handling: chunked Cuccaro + Gidney held carries
//! The backend is `compare_geq_chunked_middle`. The bottom `[0, n-k)` bits run an
//! in-place Cuccaro `a >= b` MAJ chain (one live carry, uncomputed exactly by the
//! self-inverse UMA). The top `[n-k, n)` bits hold `k` Gidney carries that are
//! measure-erased (MBU) on the reverse, so only `k+1` carries are live when the
//! caller body runs. The held-carry count `k` is supplied per call from the
//! schedule (`next_cmp_k`): `k = 0` is pure in-place Cuccaro (peak-safe), `k >= n`
//! is full Gidney.

use super::{B, BExt};
use crate::circuit::{QubitId};

/// Chunked `a >= b` comparator with a middle callback, the backend the GCD swap
/// comparator uses. `k` is the held-carry count (clamped to `n`).
/// `body(circ, flag)` sees `flag = (a >= b)`; `a`/`b` restored, `flag` cleaned.
pub fn compare_geq_chunked_middle<F: FnOnce(&mut B, &QubitId)>(
    circ: &mut B,
    a: &[QubitId],
    b: &[QubitId],
    flag: &QubitId,
    body: F,
    k: usize,
) {
    let n = a.len();
    assert_eq!(b.len(), n, "compare_geq_chunked_middle: a,b equal width");
    if n == 0 {
        circ.x(*flag);
        body(circ, flag);
        circ.x(*flag);
        return;
    }
    let k = k.min(n);
    let split = n - k; // bottom [0, split) in-place; top [split, n) held.
    let mut cy: Vec<Option<QubitId>> = (0..=n).map(|_| None).collect();
    let c = circ.alloc_qubit();
    circ.x(c); // c_0 = 1
    // Forward bottom: in-place Cuccaro MAJ (only `c` live).
    for i in 0..split {
        circ.x(b[i]);
        circ.cx(c, b[i]);
        circ.cx(c, a[i]);
        circ.ccx(a[i], b[i], c); // c = c_{i+1}
    }
    cy[split] = Some(c);
    // Forward top: Gidney held carries.
    for i in split..n {
        let next = circ.alloc_qubit();
        {
            let ci = cy[i].as_ref().unwrap();
            circ.x(b[i]);
            circ.cx(*ci, b[i]);
            circ.cx(*ci, a[i]);
            circ.ccx(a[i], b[i], next);
            circ.cx(*ci, next); // next = c_{i+1}
        }
        cy[i + 1] = Some(next);
    }
    let cn = cy[n].as_ref().unwrap();
    circ.cx(*cn, *flag); // flag = c_n = (a >= b)
    body(circ, flag);
    let cn = cy[n].as_ref().unwrap();
    circ.cx(*cn, *flag); // clean flag
    // Reverse top: measure-erase the held carries.
    for i in (split..n).rev() {
        let next = cy[i + 1].take().unwrap();
        circ.cx(*cy[i].as_ref().unwrap(), next); // c_{i+1} -> ta&tb
        let bit = circ.alloc_bit();
        circ.hmr(next, bit);
        circ.zero_and_free(next);
        circ.cz_if_bit(a[i], b[i], bit);
        circ.cx(*cy[i].as_ref().unwrap(), a[i]);
        circ.cx(*cy[i].as_ref().unwrap(), b[i]);
        circ.x(b[i]);
    }
    // Reverse bottom: in-place UMA.
    let c = cy[split].take().unwrap();
    for i in (0..split).rev() {
        circ.ccx(a[i], b[i], c);
        circ.cx(c, a[i]);
        circ.cx(c, b[i]);
        circ.x(b[i]);
    }
    circ.x(c); // c_0 = 1 -> 0
    circ.zero_and_free(c);
}

/// Controlled form: `target ^= ctrl AND (u_top < v_top)` on the top `k` MSBs.
/// Used where the swap decision is itself gated.
pub fn controlled_swap_decision_lt_truncated(
    circ: &mut B,
    ctrl: &QubitId,
    u: &[QubitId],
    v: &[QubitId],
    k: usize,
    target: &QubitId,
) {
    assert!(k <= u.len() && k <= v.len(), "k must fit in both operands");
    let u_top: Vec<QubitId> = u[u.len() - k..].to_vec();
    let v_top: Vec<QubitId> = v[v.len() - k..].to_vec();
    // Held-carry count for the chunked comparator, supplied per call by the
    // schedule (`ck` Gidney carries held, `ck+1` live).
    let ck = super::next_cmp_k();
    let lt_flag = circ.alloc_qubit();
    compare_geq_chunked_middle(
        circ,
        &u_top,
        &v_top,
        &lt_flag,
        |c, flag| {
            c.x(*flag);
            c.ccx(*ctrl, *flag, *target);
            c.x(*flag);
        },
        ck,
    );
    circ.zero_and_free(lt_flag);
}

/// Measurement-vented `a + b + cin` carry chain with a middle callback. Computes
/// the chain `cy[i+1] = carry of (a + ~b + ~cin)` bit-by-bit, then at the top bit
/// hands `(ta = a_top ^ c, tb = ~b_top ^ c, c_{n-1})` to `body` -- the final carry
/// `cy[n] = (ta AND tb) XOR c_{n-1}` is NOT built, so `body` deposits its phase via
/// a bare Z/CZ/CCZ on those three (no value flip), riding through the reverse
/// measure-uncompute. Reverse vents each internal carry by `hmr` + `cz_if_bit`.
/// `a`,`b`,`cin` restored. Equal-width `a`,`b` (the chunked-erase caller).
///
/// `carry-out(a + b + cin) = NOT carry-out(a + ~b + ~cin)`, so a caller wanting
/// to test `[a + b + cin >= 2^n]` reads the complement of the built predicate.
pub fn compare_geq_cin_middle<F: FnOnce(&mut B, &QubitId, &QubitId, &QubitId)>(
    circ: &mut B,
    a: &[QubitId],
    b: &[QubitId],
    cin: &QubitId,
    body: F,
) {
    let n = a.len();
    assert_eq!(b.len(), n, "compare_geq_cin_middle: a,b equal width");
    assert!(n >= 1, "needs >= 1 bit");
    let mut cy: Vec<Option<QubitId>> = Vec::with_capacity(n);
    let c0 = circ.alloc_qubit();
    circ.x(c0);
    circ.cx(*cin, c0); // cy[0] = 1 ^ cin = ~cin (carry-in of a + ~b + ~cin)
    cy.push(Some(c0));
    for i in 0..n - 1 {
        let next = circ.alloc_qubit();
        let ci = cy[i].as_ref().unwrap();
        circ.x(b[i]);
        circ.cx(*ci, b[i]);
        circ.cx(*ci, a[i]);
        circ.ccx(a[i], b[i], next);
        circ.cx(*ci, next);
        cy.push(Some(next));
    }
    // Top bit: fold only, hand (ta, tb, c_{n-1}) to body.
    {
        let i = n - 1;
        let ci = cy[i].as_ref().unwrap();
        circ.x(b[i]);
        circ.cx(*ci, b[i]);
        circ.cx(*ci, a[i]);
        body(circ, &a[i], &b[i], ci);
        circ.cx(*ci, a[i]);
        circ.cx(*ci, b[i]);
        circ.x(b[i]);
    }
    // Reverse: vent cy[1..n-1] via hmr, restore a/b.
    for i in (0..n - 1).rev() {
        let next = cy[i + 1].take().unwrap();
        let ci_raw = cy[i].as_ref().unwrap();
        circ.cx(*ci_raw, next); // next = ta_i & tb_i
        let bit = circ.alloc_bit();
        circ.hmr(next, bit);
        circ.zero_and_free(next);
        circ.cz_if_bit(a[i], b[i], bit);
        circ.cx(*cy[i].as_ref().unwrap(), a[i]);
        circ.cx(*cy[i].as_ref().unwrap(), b[i]);
        circ.x(b[i]);
    }
    let c0 = cy[0].take().unwrap();
    circ.cx(*cin, c0); // ~cin -> 1
    circ.x(c0); // 1 -> 0
    circ.zero_and_free(c0);
}

/// Vented uncompute of a GCD swap-decision flag that holds `ctrl AND (v_top <
/// u_top)` (the forward decision). HMR the flag to |0> (0 Toffoli), then
/// under the HMR `push_condition` recompute the predicate as a deferred Z, using
/// `compare_geq_chunked_middle(v_top, u_top)`. `v`,`u` restored. The forward
/// computes the flag normally (a value); only this reverse clear vents.
///
/// `[v >= u] = carryout(v + ~u + 1)` over the top-`k` window; the predicate is
/// `ctrl AND (v < u) = ctrl AND NOT[v>=u]`. With `lt_flag = [v_top >= u_top]`,
/// deposit `Z^(ctrl AND NOT lt_flag)` via `X(lt_flag); CZ(ctrl, lt_flag);
/// X(lt_flag)`.
pub fn swap_decision_uncompute_vented(
    circ: &mut B,
    ctrl: &QubitId,
    v: &[QubitId],
    u: &[QubitId],
    k: usize,
    flag: &QubitId,
) {
    assert!(k <= v.len() && k <= u.len(), "k must fit in both operands");
    let v_top: Vec<QubitId> = v[v.len() - k..].to_vec();
    let u_top: Vec<QubitId> = u[u.len() - k..].to_vec();
    // Held-carry count for the chunked comparator, supplied per call by the
    // schedule (matches the forward decision).
    let ck = super::next_cmp_k();
    let bit = circ.alloc_bit();
    circ.hmr(*flag, bit);
    circ.push_condition(bit);
    let lt_flag = circ.alloc_qubit(); // = [v >= u]
    compare_geq_chunked_middle(
        circ,
        &v_top,
        &u_top,
        &lt_flag,
        |c, fl| {
            // deposit Z^(ctrl AND NOT fl) = Z^(ctrl AND [v < u]), gated by the HMR
            // condition (push_condition). Same phase as the cin (ta,tb,c_prev) form.
            c.x(*fl);
            c.cz(*ctrl, *fl);
            c.x(*fl);
        },
        ck,
    );
    circ.zero_and_free(lt_flag);
    circ.pop_condition();
}

