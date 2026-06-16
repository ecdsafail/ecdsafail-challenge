//! Physical arithmetic primitives for the secp256k1 circuit.
//!
//! Uses MBUC (HMR + CZ phase correction) for carry cleanup.
//! All circuits are physical-only: no selfwire, no overlap CCX.
//!
//! A handful of helper fns (`controlled_mod_halve_secp256k1`,
//! `alloc_work_reg`, `varwidth_add`, etc.) are kept around as
//! reference-implementations for future work.

use crate::point_add::trailmix_port::arith::const_add::get_const_bit;
use crate::point_add::trailmix_port::circuit::{Circuit, QReg};

#[cfg(test)]
use crate::point_add::trailmix_port::arith::compare::*;

// === Addition / Subtraction ===

/// Addition: a += b via canonical Cuccaro MAJ/UMA (arXiv:
/// quant-ph/0410184). 1 clean ancilla, 2n Toffoli, 4n CX. Polylog
/// peak (vs `add_physical`'s n-1 AND ancs).
pub fn add(circ: &mut Circuit, a: &[QReg], b: &[QReg]) {
    crate::point_add::trailmix_port::arith::cuccaro::add_cuccaro(circ, a, b);
}

/// Subtraction: a -= b via bit-complement wrap around `add_physical`.
pub fn sub(circ: &mut Circuit, a: &[QReg], b: &[QReg]) {
    for q in a {
        circ.x(q);
    }
    add(circ, a, b);
    for q in a {
        circ.x(q);
    }
}

// === Comparison ===

// =====================================================================
// Inline-phase comparators (notes/MBUC_GADGETS.md §7).
//
// These do `forward MAJ + push_condition(bit); <phase block>;
// pop_condition() + backward UMA` in a single call. They fold the
// `temp + compare-twice` pattern in the mod_*_mbu phase corrections
// into a single pass -- about half the comparison cost per call.
//
// Each helper internally allocates one carry qubit and any
// register-extension qubits, and frees them via R (zeroed by
// the backward UMA).
// =====================================================================

// === Modular arithmetic ===

/// Controlled sub-constant: if ctrl=1, a -= val. XORs `a` into
/// two's-complement form, adds, XORs back -- but XORs are gated
/// by ctrl via CX(ctrl, a[i]) would permanently flip; instead,
/// wrap via `a := ~a; a += ctrl*val; a := ~a` only if ctrl,
/// which is wasteful. Simpler approach: just delegate by `XORing`
/// val (classical NOT) and adding ctrl and `ctrl` itself
/// (two's complement +1). Since val is constant, ~val is also
/// constant, so we can call `controlled_add_const` with ~val and
/// additionally add ctrl at position 0.
/// `a += c (mod 2^a.len())` where `c` is a classical-bit register.
///
/// Mirrors `controlled_add_const` but the per-bit decision is a runtime
/// classical Cbit instead of a compile-time constant: for each i in
/// `0..c.len()` we load `c[i]` into a fresh `QReg` `ctrl`, run an
/// unconditional Häner-style halving `cinc_gidney_halving(a[i..], ctrl)`,
/// and uncompute `ctrl` back to |0> with a second `x_if_bit`. Cost is
/// ~n inc calls (`O(n^2)` CCX + CX).
///
/// Per-iteration overhead vs. the old `with_condition` form: 2 single-qubit
/// `x_if_bit` gates (load + uncompute) and 1 ancilla alloc/free, in exchange
/// for being able to call primitives that internally use `R` / `zero_and_free`
/// (e.g. Khattar-Gidney mcx ancilla cleanup) — those are forbidden inside
/// `push_condition` blocks.
pub fn add_creg(circ: &mut Circuit, a: &[QReg], c: &[crate::point_add::trailmix_port::circuit::Cbit]) {
    let n = a.len();
    if n == 0 || c.is_empty() {
        return;
    }
    {
        let a_for_capture: Vec<&QReg> = a.iter().collect();
        let c_ids: Vec<u32> = c.iter().map(|b| b.raw()).collect();
        let n_cap = n;
        circ.contract_capture(
            "poc_arith.add_creg",
            move |view, shot| -> Result<(crate::point_add::trailmix_port::num_bigint::BigUint, crate::point_add::trailmix_port::num_bigint::BigUint), String> {
                let mut a_pre = crate::point_add::trailmix_port::num_bigint::BigUint::from(0u32);
                for (i, q) in a_for_capture.iter().enumerate() {
                    if view.contract_read_bit_shot(q, shot) {
                        a_pre |= crate::point_add::trailmix_port::num_bigint::BigUint::from(1u32) << i;
                    }
                }
                let mut c_val = crate::point_add::trailmix_port::num_bigint::BigUint::from(0u32);
                for (i, id) in c_ids.iter().enumerate() {
                    if (view.bit_mask(*id) >> shot) & 1 == 1 {
                        c_val |= crate::point_add::trailmix_port::num_bigint::BigUint::from(1u32) << i;
                    }
                }
                let modulus = crate::point_add::trailmix_port::num_bigint::BigUint::from(1u32) << n_cap;
                Ok((a_pre, c_val % modulus))
            },
        );
    }
    let lim = c.len().min(n);
    for i in 0..lim {
        let bit = c[i];
        let sub = &a[i..];
        // Load classical bit into a fresh quantum ctrl, run the
        // unconditional halving cinc, then uncompute ctrl back to |0>.
        let ctrl = circ.alloc_qreg("creg.add.ctrl");
        circ.x_if_bit(&ctrl, bit);
        crate::point_add::trailmix_port::arith::khattar_gidney::cinc_gidney_halving(circ, sub, &ctrl);
        circ.x_if_bit(&ctrl, bit);
        circ.zero_and_free(ctrl);
    }
    {
        let a_for_check: Vec<&QReg> = a.iter().collect();
        let n_cap = n;
        circ.contract_pop_and_check::<(crate::point_add::trailmix_port::num_bigint::BigUint, crate::point_add::trailmix_port::num_bigint::BigUint), _>(
            "poc_arith.add_creg",
            move |captured, view, shot| -> Result<(), String> {
                let (a_pre, c_val) = captured;
                let mut a_post = crate::point_add::trailmix_port::num_bigint::BigUint::from(0u32);
                for (i, q) in a_for_check.iter().enumerate() {
                    if view.contract_read_bit_shot(q, shot) {
                        a_post |= crate::point_add::trailmix_port::num_bigint::BigUint::from(1u32) << i;
                    }
                }
                let modulus = crate::point_add::trailmix_port::num_bigint::BigUint::from(1u32) << n_cap;
                let expected = (a_pre + c_val) % &modulus;
                if a_post != expected {
                    return Err(format!(
                        "shot {shot}: a_post = {a_post:#x}, expected a_pre + c mod 2^{n_cap} = {expected:#x}"
                    ));
                }
                Ok(())
            },
        );
    }
}

/// `a -= c (mod 2^a.len())` where `c` is a classical-bit register.
///
/// Identity: `a - c = ~(~a + c)` (two's complement). Implemented as
/// `X-flip a; add_creg(a, c); X-flip a`. The bracketing X's are cheap
/// (2n gates) compared to the inner `add_creg`, and avoid duplicating
/// the per-bit-decrement logic.
pub fn sub_creg(circ: &mut Circuit, a: &[QReg], c: &[crate::point_add::trailmix_port::circuit::Cbit]) {
    if a.is_empty() {
        return;
    }
    let n = a.len();
    {
        let a_for_capture: Vec<&QReg> = a.iter().collect();
        let c_ids: Vec<u32> = c.iter().map(|b| b.raw()).collect();
        let n_cap = n;
        circ.contract_capture(
            "poc_arith.sub_creg",
            move |view, shot| -> Result<(crate::point_add::trailmix_port::num_bigint::BigUint, crate::point_add::trailmix_port::num_bigint::BigUint), String> {
                let mut a_pre = crate::point_add::trailmix_port::num_bigint::BigUint::from(0u32);
                for (i, q) in a_for_capture.iter().enumerate() {
                    if view.contract_read_bit_shot(q, shot) {
                        a_pre |= crate::point_add::trailmix_port::num_bigint::BigUint::from(1u32) << i;
                    }
                }
                let mut c_val = crate::point_add::trailmix_port::num_bigint::BigUint::from(0u32);
                for (i, id) in c_ids.iter().enumerate() {
                    if (view.bit_mask(*id) >> shot) & 1 == 1 {
                        c_val |= crate::point_add::trailmix_port::num_bigint::BigUint::from(1u32) << i;
                    }
                }
                let modulus = crate::point_add::trailmix_port::num_bigint::BigUint::from(1u32) << n_cap;
                Ok((a_pre, c_val % modulus))
            },
        );
    }
    for q in a {
        circ.x(q);
    }
    add_creg(circ, a, c);
    for q in a {
        circ.x(q);
    }
    {
        let a_for_check: Vec<&QReg> = a.iter().collect();
        let n_cap = n;
        circ.contract_pop_and_check::<(crate::point_add::trailmix_port::num_bigint::BigUint, crate::point_add::trailmix_port::num_bigint::BigUint), _>(
            "poc_arith.sub_creg",
            move |captured, view, shot| -> Result<(), String> {
                let (a_pre, c_val) = captured;
                let mut a_post = crate::point_add::trailmix_port::num_bigint::BigUint::from(0u32);
                for (i, q) in a_for_check.iter().enumerate() {
                    if view.contract_read_bit_shot(q, shot) {
                        a_post |= crate::point_add::trailmix_port::num_bigint::BigUint::from(1u32) << i;
                    }
                }
                let modulus = crate::point_add::trailmix_port::num_bigint::BigUint::from(1u32) << n_cap;
                let expected = if a_pre >= c_val {
                    a_pre - c_val
                } else {
                    &modulus + a_pre - c_val
                };
                if a_post != expected {
                    return Err(format!(
                        "shot {shot}: a_post = {a_post:#x}, expected a_pre - c mod 2^{n_cap} = {expected:#x}"
                    ));
                }
                Ok(())
            },
        );
    }
}

/// Unconditional add-constant: a += val (mod 2^|a|).
///
/// Same dispatch as `controlled_add_const` but omits the ctrl qubit
/// entirely from inner gate structures. Saves O(log n) Toffoli per
/// cinc and the per-AND ctrl input in the Vandaele CQ-add path.
pub fn add_const(circ: &mut Circuit, a: &[QReg], val: &[u8]) {
    let n = a.len();
    if n == 0 {
        return;
    }
    let mut lo_bit = usize::MAX;
    let mut pop = 0usize;
    for i in 0..n {
        if get_const_bit(val, i) {
            if lo_bit == usize::MAX {
                lo_bit = i;
            }
            pop += 1;
        }
    }
    if pop == 0 {
        return;
    }
    if pop == 1 {
        // Single bit: just one inc from that position.
        crate::point_add::trailmix_port::arith::khattar_gidney::inc_khattar_gidney(circ, &a[lo_bit..]);
        return;
    }
    // General case: unconditional Vandaele CQ-add. Allocate the carry
    // ancilla locally; classical_quantum_add zeros it within each
    // recursion level.
    let g = circ.alloc_qreg("add_const_g");
    crate::point_add::trailmix_port::arith::khattar_gidney::classical_quantum_add(circ, a, val, &g);
    circ.zero_and_free(g);
}

/// Unconditional subtract-constant: a -= val (mod 2^|a|).
pub fn sub_const(circ: &mut Circuit, a: &[QReg], val: &[u8]) {
    let n = a.len();
    if n == 0 {
        return;
    }
    // Compute -val mod 2^n = ~val + 1 (n-bit two's complement) so we
    // can subtract via a single add_const.
    let mut neg_bits = vec![false; n];
    let mut carry = true;
    for i in 0..n {
        let inv = !get_const_bit(val, i);
        neg_bits[i] = inv ^ carry;
        carry = inv && carry;
    }
    let mut neg_val = vec![0u8; n.div_ceil(8)];
    for i in 0..n {
        if neg_bits[i] {
            neg_val[i / 8] |= 1u8 << (i % 8);
        }
    }
    add_const(circ, a, &neg_val);
}

/// Controlled add (quantum b): if ctrl=1, a += b.
/// Cuccaro (polylog peak) when ctrl does not alias a or b. Falls
/// back to `controlled_add_physical` when aliasing is detected
/// (Cuccaro's CCX(ctrl, b[i], b[i-1]) self-wires at i=|b|-1 when
/// ctrl=b[i]; the physical form has per-bit aliasing guards).
/// The aliasing case happens on squaring (result = x^2) where the
/// multiplicand and multiplier are the same register.
pub fn controlled_add(circ: &mut Circuit, ctrl: &QReg, a: &[QReg], b: &[QReg]) {
    let aliases_a = qreg_slice_contains(a, ctrl);
    let aliases_b = qreg_slice_contains(b, ctrl);
    // Pre/post contract: only on the non-aliasing path (the only one
    // the EC inversion exercises; aliasing path is squaring-specific
    // and out of scope for the cap_a / multi_add cascade debug.)
    let do_contract = !aliases_a && !aliases_b;
    if do_contract {
        let a_for_capture: Vec<&QReg> = a.iter().collect();
        let b_for_capture: Vec<&QReg> = b.iter().collect();
        let ctrl_cap: &QReg = ctrl;
        circ.contract_capture(
            "poc_arith.controlled_add",
            |view, shot| -> Result<(crate::point_add::trailmix_port::num_bigint::BigUint, crate::point_add::trailmix_port::num_bigint::BigUint, bool), String> {
                let read = |regs: &[&QReg]| -> crate::point_add::trailmix_port::num_bigint::BigUint {
                    let mut v = crate::point_add::trailmix_port::num_bigint::BigUint::from(0u32);
                    for (i, q) in regs.iter().enumerate() {
                        if view.contract_read_bit_shot(q, shot) {
                            v |= crate::point_add::trailmix_port::num_bigint::BigUint::from(1u32) << i;
                        }
                    }
                    v
                };
                let a_v = read(&a_for_capture);
                let b_v = read(&b_for_capture);
                let c_v = view.contract_read_bit_shot(ctrl_cap, shot);
                Ok((a_v, b_v, c_v))
            },
        );
    }
    if !aliases_a && !aliases_b {
        // 3n CCX controlled Cuccaro. Forward MAJ on a single carry
        // register (1 CCX per bit = n CCX); reverse pass per bit
        // (CCX(a,b,c) + CCX(ctrl,b,a) = 2 CCX = 2n CCX). Total 3n CCX,
        // vs controlled_add_cuccaro_mbu's 8n CCX. Same semantics: a += b
        // when ctrl=1, unchanged when ctrl=0; b and ctrl preserved.
        crate::point_add::trailmix_port::arith::cuccaro::controlled_add_cuccaro_3n(circ, ctrl, a, b);
    } else {
        // ctrl aliases a would be problematic: Cuccaro modifies a, so ctrl's
        // value would drift mid-add. The controlled_mod_add_rfold_mbu
        // contract only calls this with ctrl aliasing b (never a), so we
        // reject the a-alias case here to catch misuse early.
        assert!(
            !aliases_a,
            "controlled_add: ctrl aliases a register -- unsupported"
        );

        // ctrl aliases b -- copy to fresh scratch and use the 3n variant.
        // The 3n form preserves b across the add, so ctrl (= b[i]) is
        // restored at the end and the final cx(ctrl, scratch) zeros
        // scratch cleanly.
        // Peak: +2 ancillae (scratch + 3n's carry register). Polylog.
        let scratch = circ.alloc_qreg("cadd_alias_scratch");
        circ.cx(ctrl, &scratch);
        circ.declare_copy_of(&scratch, ctrl);
        crate::point_add::trailmix_port::arith::cuccaro::controlled_add_cuccaro_3n(circ, &scratch, a, b);
        circ.cx(ctrl, &scratch);
        // scratch drops here; drain fires at next gate (gap=0).
    }
    if do_contract {
        let a_for_check: Vec<&QReg> = a.iter().collect();
        let n_cap = a.len();
        circ.contract_pop_and_check::<(crate::point_add::trailmix_port::num_bigint::BigUint, crate::point_add::trailmix_port::num_bigint::BigUint, bool), _>(
            "poc_arith.controlled_add",
            move |captured, view, shot| -> Result<(), String> {
                let (a_pre, b_pre, c_v) = captured;
                let mut a_post = crate::point_add::trailmix_port::num_bigint::BigUint::from(0u32);
                for (i, q) in a_for_check.iter().enumerate() {
                    if view.contract_read_bit_shot(q, shot) {
                        a_post |= crate::point_add::trailmix_port::num_bigint::BigUint::from(1u32) << i;
                    }
                }
                let modulus = crate::point_add::trailmix_port::num_bigint::BigUint::from(1u32) << n_cap;
                let expected = if *c_v {
                    (a_pre + b_pre) % &modulus
                } else {
                    a_pre.clone()
                };
                if a_post != expected {
                    return Err(format!(
                        "shot {}: a_post = {:#x}, expected (ctrl={}: a_pre {} b_pre) mod 2^{} = {:#x}",
                        shot, a_post, c_v, if *c_v { "+" } else { "unchanged" }, n_cap, expected
                    ));
                }
                Ok(())
            },
        );
    }
}

/// Pointer-equality test for `QReg` slice membership. Since `QReg` is non-Copy
/// and identity is by qubit-id (which is module-private), we compare by
/// reference identity: a slice contains `q` iff one of its elements is the
/// same `QReg` instance. (For the alias-detection use case, callers pass the
/// same `QReg` references, so reference identity matches qubit identity.)
fn qreg_slice_contains(slice: &[QReg], q: &QReg) -> bool {
    slice.iter().any(|s| std::ptr::eq(s, q))
}

/// Top-K add-overflow phase-correction MBU. HMRs `q_to_hmr` with the
/// identity `q_to_hmr ≡ ctrl AND 1[a_top_k + b_top_k overflows]`.
///
/// Builds the K-bit ripple-carry MAJ chain over the top K bits of
/// (a + b) (without materializing the sum) to read the carry-out,
/// then runs the matching Cuccaro UMA chain to restore a, b, c.
///
/// Cost: ~2K Toffoli (MAJ chain + UMA chain) + 2 ancilla qubits.
///
/// This is the structural counterpart of `controlled_lt_msbs` for the
/// FORWARD mod-sub Alg-11 cleanup: forward sub leaves y[n] = borrow,
/// and `borrow ≡ ctrl AND 1[y_top + x_top overflows top K]` modulo a
/// 2^-K approximation tail (matches forward mod-add's 2^-K tail).
pub fn controlled_add_overflow_msbs_phase_correction_mbu(
    circ: &mut Circuit,
    ctrl: &QReg,
    a: &[QReg],
    b: &[QReg],
    q_to_hmr: &QReg,
    k: usize,
) {
    let n = a.len();
    assert_eq!(n, b.len(), "topk add-overflow requires equal a/b lengths");
    if n == 0 || k == 0 {
        let bit = circ.alloc_bit();
        circ.hmr(q_to_hmr, bit);
        circ.free_bit(bit);
        return;
    }
    let k = k.min(n);
    let lo = n - k;

    let ctrl_copy = circ.alloc_qreg("addovf_ctrl_copy");
    circ.cx(ctrl, &ctrl_copy);

    let carry = circ.alloc_qreg("addovf_carry");

    // carry init = 1, so the MAJ chain effectively computes a + b + 1
    // and returns carry-out = 1[a_top + b_top >= 2^K - 1]. This handles
    // the boundary case where a_top + b_top = 2^K - 1 (which corresponds
    // to a borrow case in the full mod-sub when low bits propagate a
    // carry up). The forward mod-add's `controlled_lt_msbs` cleanup
    // has the analogous +1 hidden inside the borrow-chain init.
    circ.x(&carry);

    // Cuccaro MAJ chain over top K bits with carry-in = 1.
    for i in lo..n {
        circ.cx(&carry, &b[i]);
        circ.cx(&carry, &a[i]);
        circ.ccx(&a[i], &b[i], &carry);
    }

    // Capture: q_to_hmr identity = ctrl AND carry.
    circ.declare_and_of(q_to_hmr, &ctrl_copy, &carry);
    let bit = circ.alloc_bit();
    circ.hmr(q_to_hmr, bit);
    circ.cz_if_bit(&ctrl_copy, &carry, bit);
    circ.free_bit(bit);

    // UMA chain (Cuccaro un-MAJ) to restore a, b, carry. NO inner
    // controlled add — we only want the carry-out, not the sum.
    for i in (lo..n).rev() {
        circ.ccx(&a[i], &b[i], &carry);
        circ.cx(&carry, &a[i]);
        circ.cx(&carry, &b[i]);
    }
    // Restore carry to |0> by undoing the initial X.
    circ.x(&carry);
    drop(carry);

    circ.cx(ctrl, &ctrl_copy);
    drop(ctrl_copy);
}

/// UNCONTROLLED top-k add-overflow flag clear: clears `q_to_clean`
/// knowing it equals the top-`k` add-overflow of (a, b).
///
/// The ctrl-free form of [`controlled_add_overflow_msbs_phase_correction_mbu`],
/// used by the unconditional pseudo-Mersenne mod-sub cleanup. After a
/// mod-sub, `q_to_clean` holds the borrow flag, and the identity
/// `borrow == 1[a_top + b_top + 1 overflows K bits]` (= the carry-out of
/// the carry-in-1 MAJ chain) lets us clear it with a single reversible
/// `cx(carry, q_to_clean)` — no Toffoli, no HMR — once the tracker is
/// told `q_to_clean` is a copy of `carry`.
pub fn add_overflow_msbs_phase_correction(
    circ: &mut Circuit,
    a: &[QReg],
    b: &[QReg],
    q_to_clean: &QReg,
    k: usize,
) {
    let n = a.len();
    assert_eq!(n, b.len(), "topk add-overflow requires equal a/b lengths");
    assert!(n > 0 && k > 0, "uncontrolled add-overflow needs k >= 1");
    let k = k.min(n);
    let lo = n - k;

    let carry = circ.alloc_qreg("addovf_carry");
    // carry-in = 1: the MAJ chain returns carry-out = 1[a_top + b_top >= 2^K - 1].
    circ.x(&carry);
    for i in lo..n {
        circ.cx(&carry, &b[i]);
        circ.cx(&carry, &a[i]);
        circ.ccx(&a[i], &b[i], &carry);
    }
    // q_to_clean == carry. Tell the tracker, then clear reversibly (carry
    // is read-only here, so the MAJ window's restore stays valid).
    circ.declare_copy_of(q_to_clean, &carry);
    circ.cx(&carry, q_to_clean);
    // UMA chain (un-MAJ) to restore a, b, carry.
    for i in (lo..n).rev() {
        circ.ccx(&a[i], &b[i], &carry);
        circ.cx(&carry, &a[i]);
        circ.cx(&carry, &b[i]);
    }
    circ.x(&carry);
    drop(carry);
}

/// Controlled sub (quantum b): if ctrl=1, a -= b.
/// X-sandwich around `controlled_add` (Cuccaro). `ctrl=0` case: the
/// NOT/NOT wraps cancel and `controlled_add` adds 0 -> a unchanged.
/// `ctrl=1` case: a <- ~(~a + b) = a - b (mod 2^n). Correct.
///
/// Peak 1 anc (the 3n controlled-add's carry register). Gates: 3n
/// Toffoli + (CX from the inner add) + 2n X (the outer sandwich).
pub fn controlled_sub(circ: &mut Circuit, ctrl: &QReg, a: &[QReg], b: &[QReg]) {
    let n = a.len();
    if n == 0 {
        return;
    }
    // Callers must provide b at least as wide as a (padding into a
    // synthetic Vec<QReg> would require either Clone or fresh-anc
    // copies; the original Qubit-typed code padded with fresh zero
    // qubits and freed them after — under the Qubit-private regime
    // we require the caller to pass the padded slice in directly).
    assert!(
        b.len() >= n,
        "controlled_sub: b ({} bits) shorter than a ({} bits); \
         caller must provide a same-width or wider b slice",
        b.len(),
        n
    );
    let b_eq = &b[..n];
    for q in a {
        circ.x(q);
    }
    controlled_add(circ, ctrl, a, b_eq);
    for q in a {
        circ.x(q);
    }
}

// === Variable-width helpers ===

#[cfg(test)]
mod tests {
    use super::compare_geq_gidney_middle;
    use crate::point_add::trailmix_port::circuit::QReg;
    use crate::point_add::trailmix_port::circuit::Circuit;

    #[test]
    fn compare_geq_gidney_middle_random() {
        use rand::Rng;
        let nbits = 16usize;
        let mut circ = Circuit::new();
        let a = circ.alloc_qreg_bits("a", nbits);
        let b = circ.alloc_qreg_bits("b", nbits);
        let flag = circ.alloc_qreg("flag");
        let target = circ.alloc_qreg("target");

        let mut rng = rand::thread_rng();
        let mut a_pre = [0u32; 64];
        let mut b_pre = [0u32; 64];
        for shot in 0..64 {
            a_pre[shot] = rng.gen::<u32>() & 0xffff;
            b_pre[shot] = rng.gen::<u32>() & 0xffff;
            circ.sim_load_reg_bytes_shot(&a, &a_pre[shot].to_le_bytes()[..2], shot);
            circ.sim_load_reg_bytes_shot(&b, &b_pre[shot].to_le_bytes()[..2], shot);
        }

        compare_geq_gidney_middle(&mut circ, &a, &b, &flag, |c, fl| {
            c.cx(fl, &target); // capture (a >= b) into target
        });
        circ.assert_phase_clean();

        let mut outputs: Vec<QReg> = Vec::new();
        outputs.extend(a);
        outputs.extend(b);
        outputs.push(flag);
        outputs.push(target);
        let (sim, det) = circ.destroy_sim(outputs);
        for shot in 0..64 {
            let got_t = sim.read_bit_shot(&det[2 * nbits + 1], shot);
            let exp = if a_pre[shot] >= b_pre[shot] { 1 } else { 0 };
            assert_eq!(got_t, exp, "shot {shot}: (a>=b) mismatch");
            assert_eq!(
                sim.read_bit_shot(&det[2 * nbits], shot),
                0,
                "shot {shot}: flag not 0"
            );
            let mut got_a = 0u32;
            let mut got_b = 0u32;
            for i in 0..nbits {
                if sim.read_bit_shot(&det[i], shot) == 1 {
                    got_a |= 1 << i;
                }
                if sim.read_bit_shot(&det[nbits + i], shot) == 1 {
                    got_b |= 1 << i;
                }
            }
            assert_eq!(got_a, a_pre[shot], "shot {shot}: a not restored");
            assert_eq!(got_b, b_pre[shot], "shot {shot}: b not restored");
        }
    }
}
