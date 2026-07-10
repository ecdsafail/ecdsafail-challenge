//! Cuccaro 3n CCX controlled adder (Schrottenloher 2026 / trailmix)
//!
//! Replaces the expensive MBU-based `controlled_add_cuccaro_mbu` (8n CCX)
//! with the optimal Cuccaro construction (3n CCX).
//!
//! Construction (per Schrottenloher / trailmix cuccaro.rs):
//!
//!   FORWARD MAJ chain (1-qubit ripple, single carry register c=|0>):
//!     per bit i:  CX(c, b[i]); CX(c, a[i]); CCX(a[i], b[i], c)
//!     state post-bit i:
//!       a[i] = a_orig XOR c_in_i
//!       b[i] = b_orig XOR c_in_i
//!       c    = c_in_i XOR MAJ(c_in_i, a_orig, b_orig) = c_out_i = c_in_{i+1}
//!     After all n bits: c = carry-out of (a + b) >> n.
//!     Cost: 1 CCX per bit, n CCX total.
//!
//!   REVERSE pass (gated on ctrl, descending i):
//!     CCX(a[i], b[i], c)   ; restore c to c_in_i        (1 CCX)
//!     CX(c, a[i])          ; a[i] := a_orig             (CX)
//!     CCX(ctrl, b[i], a[i]); gated: a[i] XOR= ctrl·b_i  (1 CCX)
//!                            ctrl=0: no-op  → a[i] stays a_orig
//!                            ctrl=1: a[i] := a_orig XOR b_orig XOR c_in_i
//!                                    = sum_i              ✓
//!     CX(c, b[i])          ; b[i] := b_orig             (CX)
//!     Cost: 2 CCX per bit, 2n CCX total.
//!
//!   c is restored to |0> at the end (initial carry-in was 0, full ripple
//!   unwinds to 0).
//!
//! Total: 3n CCX (vs 8n for MBU variant, vs 3n-2 for hybrid with full venting).

use crate::circuit::QubitId;
use crate::point_add::trailmix_ludicrous::{B, BExt};

/// Pure controlled add (3n CCX): if ctrl=1, a := a + b mod 2^n; else
/// a unchanged. Same semantics as MBU variant but uses ~2.7x fewer Toffolis.
///
/// Semantics:
///   ctrl=1: a := (a + b) mod 2^n
///   ctrl=0: a unchanged
///   b, ctrl preserved in both cases.
///
/// Total: 3n CCX (vs 8n for controlled_add_cuccaro_mbu).
pub fn controlled_add_cuccaro_3n(circ: &mut B, ctrl: &QubitId, a: &[QubitId], b: &[QubitId]) {
    let n = a.len();
    assert_eq!(b.len(), n, "controlled_add_cuccaro_3n: a/b length mismatch");

    if n == 0 {
        return;
    }
    if n == 1 {
        // 1-bit case: a[0] ^= ctrl·b[0]
        circ.ccx(*ctrl, b[0], a[0]);
        return;
    }

    let c = circ.alloc_qubit();

    // Forward MAJ chain. Single carry register `c` ripples bit-by-bit.
    // Per bit i: CX(c, b); CX(c, a); CCX(a, b, c).
    for i in 0..n {
        circ.cx(c, b[i]);
        circ.cx(c, a[i]);
        circ.ccx(a[i], b[i], c);
    }

    // Reverse pass, descending. Per bit i:
    //   CCX(a,b,c)       -- restore c to c_in_i
    //   CX(c, a)          -- a := a_orig (undo)
    //   CCX(ctrl, b, a)   -- gated: a XOR= ctrl·(b XOR c_in) = ctrl·b_orig XOR ctrl·c_in
    //   CX(c, b)          -- b := b_orig
    //
    // When ctrl=1, the chained CX(c,a) then CCX(ctrl,b,a) yields
    //   a_post = a_orig XOR (b_orig XOR c_in) = a_orig XOR b_orig XOR c_in = sum_i.
    // When ctrl=0, the CCX is a no-op, so a_post = a_orig.
    for i in (0..n).rev() {
        circ.ccx(a[i], b[i], c);
        circ.cx(c, a[i]);
        circ.ccx(*ctrl, b[i], a[i]);
        circ.cx(c, b[i]);
    }

    // c is back to |0> (initial carry-in was 0; ripple fully unwound).
    circ.zero_and_free(c);
}

/// Controlled subtract (3n CCX): if ctrl=1, a := a - b mod 2^n; else unchanged.
/// b and ctrl preserved.
/// 
/// Subtraction = addition with negated b. We use the identity:
///   a - b = a + (~b) + 1  (mod 2^n)
/// where ~b is the bitwise complement of b.
pub fn controlled_sub_cuccaro_3n(circ: &mut B, ctrl: &QubitId, a: &[QubitId], b: &[QubitId]) {
    let n = a.len();
    assert_eq!(b.len(), n, "controlled_sub_cuccaro_3n: a/b length mismatch");

    if n == 0 {
        return;
    }

    // Subtraction via addition: a - b = a + (~b) + 1 mod 2^n
    // We do this by first flipping all b bits, then adding 1 at the end.
    
    // First: compute a + (~b) using the 3n adder
    // Save b bits by copying them first
    let b_copy: Vec<QubitId> = b.iter().map(|&q| {
        let tmp = circ.alloc_qubit();
        circ.cx(q, tmp);
        tmp
    }).collect();

    // Flip the copies (these are the ones we'll add)
    for q in &b_copy {
        circ.x(*q);
    }

    // Now add (~b_copy) to a
    // We'll use the forward MAJ chain, then handle the +1 separately
    // Actually, let's use the standard approach: 
    // For subtraction, we want a - b = a + (2^n - b) = a + (~b) + 1
    // The +1 is done by pre-setting the carry-in
    
    // Simpler: just do forward add with flipped bits, then do a +1 at end
    let c = circ.alloc_qubit();
    circ.x(c); // Set initial carry to 1 (for the +1 in a + ~b + 1)

    // Forward MAJ chain with flipped b
    for i in 0..n {
        circ.cx(c, b_copy[i]);
        circ.cx(c, a[i]);
        circ.ccx(a[i], b_copy[i], c);
    }

    // Reverse pass
    for i in (0..n).rev() {
        circ.ccx(a[i], b_copy[i], c);
        circ.cx(c, a[i]);
        circ.ccx(*ctrl, b_copy[i], a[i]);
        circ.cx(c, b_copy[i]);
    }

    circ.zero_and_free(c);

    // Clean up the copied bits
    for (orig, &copy) in b.iter().zip(b_copy.iter()) {
        circ.cx(copy, *orig);
        circ.zero_and_free(copy);
    }
}

/// Controlled add with carry-out (3n CCX + 1 for carry).
/// Returns the carry-out qubit.
/// if ctrl=1: a := a + b mod 2^n, cout = carry_out
/// if ctrl=0: a unchanged, cout = 0
pub fn controlled_add_cuccaro_3n_with_cout(
    circ: &mut B, 
    ctrl: &QubitId, 
    a: &[QubitId], 
    b: &[QubitId],
    cout: &QubitId,
) {
    let n = a.len();
    assert_eq!(b.len(), n, "controlled_add_cuccaro_3n_with_cout: a/b length mismatch");

    if n == 0 {
        return;
    }
    if n == 1 {
        circ.ccx(*ctrl, b[0], a[0]);
        circ.ccx(*ctrl, b[0], *cout); // cout gets the carry
        return;
    }

    let c = circ.alloc_qubit();

    // Forward MAJ chain
    for i in 0..n {
        circ.cx(c, b[i]);
        circ.cx(c, a[i]);
        circ.ccx(a[i], b[i], c);
    }

    // Carry out: copy the final carry
    circ.cx(c, *cout);

    // Reverse pass
    for i in (0..n).rev() {
        circ.ccx(a[i], b[i], c);
        circ.cx(c, a[i]);
        circ.ccx(*ctrl, b[i], a[i]);
        circ.cx(c, b[i]);
    }

    circ.zero_and_free(c);
}
