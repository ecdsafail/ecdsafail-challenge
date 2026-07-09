//! Gidney 2025 venting adder primitives (arxiv 2507.23079).
//!
//! These primitives implement classical-quantum addition with O(1) clean
//! ancilla qubits, by "venting" carry qubits (measuring them in X basis
//! and deferring the corresponding phase-flip tasks to the end via
//! Häner-Roetteler-Soeken's carry-xor construction).
//!
//! Python reference: https://zenodo.org/doi/10.5281/zenodo.15866587
//!
//! The key primitives:
//! - [`xor_right_shifted_carries_into`]: Häner carry-xor.
//!   Performs `Q_dst ^= carry(Q_src, offset, carry_in) >> 1` in ~2n CCX
//!   using 0 clean ancilla.
//! - [`add_vented_2clean`]: streaming vented add. 2 clean ancilla, ~n CCX,
//!   leaves n-2 phase-flip tasks behind.
//! - [`iadd_3clean`]: full const-quantum add. 3 clean ancilla, 4n CCX.
//!
//! Status: initial port, API subject to change. Tests in the unit-test
//! module at the bottom.

use super::{BitId, QubitId, B};
use crate::circuit::{Op, OperationType};

/// Performs `Q_dst ^= carry(Q_src, offset, carry_in) >> 1` in-place.
///
/// Here `carry(x, d, c0)` returns an n-bit value where bit k is the carry
/// into bit k of the addition `x + d + c0` (with c0 being the bit-0
/// carry-in). The `>> 1` means we skip the LSB of the carry (which equals
/// the carry-in and is trivially accessible).
///
/// `offset` may be classical or quantum. When classical, `offset[k]` is
/// a `BitId` whose value is the k-th bit of the constant offset. When
/// quantum, `offset[k]` is a `QubitId`.
///
/// Cost: ~2n CCX, 0 clean ancilla.
///
/// # Arguments
/// - `q_src`: n+1 qubits (or n) representing the "target" of the
///   reference addition.
/// - `offset`: n classical bits (the constant to add).
/// - `q_dst`: n qubits to XOR the right-shifted carries into.
/// - `carry_in`: classical bit (0 or 1) for the LSB carry-in.
#[allow(dead_code)]
pub(crate) fn xor_right_shifted_carries_into_classical(
    b: &mut B,
    q_src: &[QubitId],
    offset_bits: u64,
    q_dst: &[QubitId],
    carry_in: bool,
) {
    let n = q_dst.len();
    assert!(n <= q_src.len() && q_src.len() <= n + 1, "len mismatch");
    if n == 0 {
        return;
    }

    // Helper: bit k of the classical offset.
    let bit = |k: usize| -> bool {
        if k >= 64 {
            false
        } else {
            (offset_bits >> k) & 1 != 0
        }
    };

    // Helper: apply CCX(ctrl_a, ctrl_b, target) with each control
    // possibly classically-inverted. The original `a ^ offset[k]` means:
    // if offset[k] = 0, use `a` directly; if offset[k] = 1, use `NOT a`.
    // We implement this via `X(a)` before and after the CCX.
    let ccx_inv =
        |b: &mut B, ctrl_a: QubitId, inv_a: bool, ctrl_b: QubitId, inv_b: bool, target: QubitId| {
            if inv_a {
                b.x(ctrl_a);
            }
            if inv_b {
                b.x(ctrl_b);
            }
            b.ccx(ctrl_a, ctrl_b, target);
            if inv_b {
                b.x(ctrl_b);
            }
            if inv_a {
                b.x(ctrl_a);
            }
        };

    // First loop (reversed over k=1..n):
    //   ccx(Q_src[k] ^ offset[k], Q_dst[k-1], Q_dst[k])
    for k in (1..n).rev() {
        ccx_inv(b, q_src[k], bit(k), q_dst[k - 1], false, q_dst[k]);
    }

    // broadcast_cx(offset, Q_dst): for each k, if offset[k]: X(Q_dst[k]).
    // (This is equivalent to XORing the classical offset into Q_dst.)
    for k in 0..n {
        if bit(k) {
            b.x(q_dst[k]);
        }
    }

    // ccx(Q_src[0] ^ offset[0], carry_in ^ offset[0], Q_dst[0])
    // carry_in is CLASSICAL here. If (carry_in XOR offset[0]) = 0, the
    // CCX has a classical-0 control and does nothing. If it's 1, the CCX
    // reduces to CX(q_src[0] with inv, q_dst[0]).
    let carry_in_xor_offset0 = carry_in ^ bit(0);
    if carry_in_xor_offset0 {
        // CX(q_src[0] ^ offset[0], q_dst[0]).
        if bit(0) {
            b.x(q_src[0]);
        }
        b.cx(q_src[0], q_dst[0]);
        if bit(0) {
            b.x(q_src[0]);
        }
    }

    // Second loop (k=1..n):
    //   ccx(Q_src[k] ^ offset[k], Q_dst[k-1] ^ offset[k], Q_dst[k])
    for k in 1..n {
        ccx_inv(b, q_src[k], bit(k), q_dst[k - 1], bit(k), q_dst[k]);
    }
}

/// Gidney 2025 streaming vented adder (Figure 2, arxiv 2507.23079).
///
/// Performs `Q_target += offset + carry_in` (mod 2^n) while using only
/// 2 clean ancilla qubits. Leaves behind n-2 "vent" phase-flip tasks in
/// classical bits `vent_keys[1..n-1]`; these must be corrected by a
/// subsequent `xor_right_shifted_carries_into` + classical-CZ sandwich
/// (see Figure 4's second half).
///
/// Uses the X-basis demolition measurement (HMR) to "vent" carries
/// eagerly as they're computed, freeing each carry qubit for reuse
/// immediately after it stops being needed by the ripple.
///
/// Cost: n ± O(1) CCX, 2 clean ancilla, n-2 classical bits for vent_keys.
///
/// # Arguments
/// - `q_target`: n qubits. On exit: target + offset + carry_in mod 2^n.
///   PLUS residual phase-flip tasks indexed by `vent_keys`.
/// - `q_clean2`: 2 clean ancilla qubits.
/// - `offset_bits`: classical n-bit offset (bit k is `(offset_bits >> k) & 1`).
/// - `carry_in`: classical carry-in bit.
/// - `vent_keys`: n classical bits. On exit: `vent_keys[k]` for k in 1..n-1
///   holds the random measurement outcome that needs phase correction later.
///   `vent_keys[0]` and `vent_keys[n-1]` are unused.
pub(crate) fn add_vented_2clean_classical(
    b: &mut B,
    q_target: &[QubitId],
    q_clean2: &[QubitId; 2],
    offset_bits: u64,
    carry_in: bool,
    vent_keys: &[BitId],
) {
    add_vented_2clean_classical_cxt(
        b,
        q_target,
        q_clean2,
        offset_bits,
        carry_in,
        vent_keys,
        None,
    );
}

/// Extended vented adder supporting optional `carry_xor_target`: during the
/// k-th ripple step, if `carry_xor_target[k]` is Some(q), emit
/// `cx(carries[k], q)` — XORing the computed carry into a target qubit
/// before it gets vented. This is used by `iadd_dirty_2clean` to merge
/// the Gidney Figure 4 carry-xor pass into the vented add itself.
pub(crate) fn add_vented_2clean_classical_cxt(
    b: &mut B,
    q_target: &[QubitId],
    q_clean2: &[QubitId; 2],
    offset_bits: u64,
    carry_in: bool,
    vent_keys: &[BitId],
    carry_xor_target: Option<&[Option<QubitId>]>,
) {
    let n = q_target.len();
    if n == 0 {
        return;
    }
    let bit = |k: usize| -> bool {
        if k >= 64 {
            false
        } else {
            (offset_bits >> k) & 1 != 0
        }
    };

    if n == 1 {
        if carry_in {
            b.x(q_target[0]);
        }
        if bit(0) {
            b.x(q_target[0]);
        }
        return;
    }

    // carries[0] = carry_in (classical).
    // carries[k] = q_clean2[k % 2] for k in 1..n-1.
    // carries[n-1] = q_target[n-1].
    // We represent carry_in as classical via branching on its value.

    // broadcast_cx(offset, q_target): for each k, if offset[k]: X(q_target[k]).
    for k in 0..n {
        if bit(k) {
            b.x(q_target[k]);
        }
    }

    // Helper to apply the CCX with classical-inverted control, and when
    // the control source is carry_in (classical), simplify.
    // carries[k] for k=0 is classical carry_in; for k=n-1 is q_target[n-1]; else ancilla.
    let get_carry_qubit = |k: usize| -> Option<QubitId> {
        if k == 0 {
            None // classical carry_in
        } else if k == n - 1 {
            Some(q_target[n - 1])
        } else {
            Some(q_clean2[k % 2])
        }
    };

    for k in 0..n - 1 {
        // if k < n-2: rz(carries[k+1]) (reset the NEXT carry qubit to |0>).
        // Since q_clean2 qubits are reused in alternation, the qubit
        // q_clean2[(k+1) % 2] needs to be at |0> before we write into it.
        // The `rz` op = R (reset to |0>).
        if k < n - 2 {
            if let Some(q) = get_carry_qubit(k + 1) {
                // Reset via R op.
                let mut op = Op::empty();
                op.kind = OperationType::R;
                op.q_target = q;
                b.ops.push(op);
            }
        }

        // ccx(q_target[k], carries[k] XOR offset[k], carries[k+1])
        // Cases based on carries[k]'s source:
        //   k==0: carries[0] = carry_in (classical bit).
        //     carries[k] XOR offset[k] = carry_in XOR bit(0), which is a classical bit.
        //     If false: CCX becomes no-op (classical-0 control).
        //     If true: CCX becomes CX(q_target[k], carries[k+1]).
        //   k>=1: carries[k] is a qubit. offset[k] inverts it.
        if k == 0 {
            let eff_carry = carry_in ^ bit(0);
            if eff_carry {
                // CX(q_target[0], carries[1])
                if let Some(q) = get_carry_qubit(1) {
                    b.cx(q_target[0], q);
                }
            }
        } else {
            let carry_q = get_carry_qubit(k).expect("non-boundary carry");
            let carry_next = get_carry_qubit(k + 1).expect("non-boundary next carry");
            if bit(k) {
                b.x(carry_q);
                b.ccx(q_target[k], carry_q, carry_next);
                b.x(carry_q);
            } else {
                b.ccx(q_target[k], carry_q, carry_next);
            }
        }

        // cx(carries[k], q_target[k])
        if k == 0 {
            if carry_in {
                b.x(q_target[0]);
            }
        } else {
            let carry_q = get_carry_qubit(k).expect("non-boundary carry");
            b.cx(carry_q, q_target[k]);
        }

        // Optional: cx(carries[k], carry_xor_target[k]) if provided.
        // (Python reference: `out.cx(carries[k], carry_xor_target[k])`)
        if let Some(cxt) = carry_xor_target {
            if k < cxt.len() {
                if let Some(dst) = cxt[k] {
                    if k == 0 {
                        if carry_in {
                            b.x(dst);
                        }
                    } else {
                        let carry_q = get_carry_qubit(k).expect("non-boundary carry");
                        b.cx(carry_q, dst);
                    }
                }
            }
        }

        // mx(carries[k], out=vent_keys[k]) for k > 0
        if k > 0 {
            let carry_q = get_carry_qubit(k).expect("non-boundary carry");
            b.hmr(carry_q, vent_keys[k]);
        }

        // cx(offset[k], carries[k+1]): if offset[k] classical: if set, X(carries[k+1]).
        if bit(k) {
            if let Some(q) = get_carry_qubit(k + 1) {
                b.x(q);
            }
        }
    }
}

/// HRS 2017 adder (arxiv 1709.06648): `Q_target += offset + carry_in`
/// using n-2 clean ancilla qubits as carry storage.
///
/// Cost: n ± O(1) CCX.
///
/// # Arguments
/// - `q_target`: n qubits (the destination register).
/// - `q_clean`: at least n-2 clean ancilla qubits.
/// - `offset_bits`: classical n-bit offset.
/// - `carry_in`: classical carry-in.
pub(crate) fn iadd_linear_clean_classical(
    b: &mut B,
    q_target: &[QubitId],
    q_clean: &[QubitId],
    offset_bits: u64,
    carry_in: bool,
) {
    let n = q_target.len();
    if n == 0 {
        return;
    }
    assert!(q_clean.len() >= n.saturating_sub(2), "need n-2 clean");
    let q_clean = &q_clean[..n.saturating_sub(2)];

    let bit = |k: usize| -> bool {
        if k >= 64 {
            false
        } else {
            (offset_bits >> k) & 1 != 0
        }
    };

    // Special case n==1:
    if n == 1 {
        if bit(0) {
            b.x(q_target[0]);
        }
        if carry_in {
            b.x(q_target[0]);
        }
        return;
    }
    // Special case n==2:
    if n == 2 {
        // carries = [carry_in, q_target[1]].
        // broadcast_cx(offset[:1], carries[1:]): if offset[0]: X(q_target[1]).
        if bit(0) {
            b.x(q_target[1]);
        }
        // broadcast_cx(offset, q_target): if offset[k]: X(q_target[k]).
        for k in 0..2 {
            if bit(k) {
                b.x(q_target[k]);
            }
        }
        // ccx loop: k=0. carries[0]=cin, carries[1]=q_target[1].
        // ccx(q_target[0], carries[0] XOR offset[0], carries[1]).
        let eff0 = carry_in ^ bit(0);
        if eff0 {
            b.cx(q_target[0], q_target[1]);
        }
        // uncompute loop: empty for n==2.
        // cx(carries[0], q_target[0]): if carry_in: X(q_target[0]).
        if carry_in {
            b.x(q_target[0]);
        }
        return;
    }

    // Reset clean ancilla (they may be dirty).
    // Python did `out.rz(q)` which is our `R` op.
    for &q in q_clean.iter() {
        let mut op = Op::empty();
        op.kind = OperationType::R;
        op.q_target = q;
        b.ops.push(op);
    }

    // carries[0] = cin (classical); carries[1..n-1] = q_clean[0..n-2]; carries[n-1] = q_target[n-1].
    let get_carry = |k: usize| -> Option<QubitId> {
        if k == 0 {
            None
        } else if k == n - 1 {
            Some(q_target[n - 1])
        } else {
            Some(q_clean[k - 1])
        }
    };

    // broadcast_cx(offset[:n-1], carries[1:]).
    // i.e. for k in 0..n-1: if offset[k]: X(carries[k+1]).
    for k in 0..n - 1 {
        if bit(k) {
            if let Some(q) = get_carry(k + 1) {
                b.x(q);
            }
        }
    }
    // broadcast_cx(offset, q_target): for k in 0..n: if offset[k]: X(q_target[k]).
    for k in 0..n {
        if bit(k) {
            b.x(q_target[k]);
        }
    }

    // Forward compute loop.
    for k in 0..n - 1 {
        // ccx(q_target[k], carries[k] XOR offset[k], carries[k+1]).
        let next = get_carry(k + 1).expect("k+1 in bounds");
        if k == 0 {
            // carries[0] = cin. cin XOR offset[0]: classical.
            let eff = carry_in ^ bit(0);
            if eff {
                b.cx(q_target[0], next);
            }
        } else {
            let cur = get_carry(k).expect("k in bounds");
            if bit(k) {
                b.x(cur);
                b.ccx(q_target[k], cur, next);
                b.x(cur);
            } else {
                b.ccx(q_target[k], cur, next);
            }
        }
    }

    // Uncompute loop (reversed, with HMR + CZ + CCZ).
    for k in (0..n - 2).rev() {
        // cx(carries[k+1], q_target[k+1]).
        let next = get_carry(k + 1).expect("k+1 in bounds");
        b.cx(next, q_target[k + 1]);
        // mx(carries[k+1], out=m). This measures next.
        let m = b.alloc_bit();
        b.hmr(next, m);
        // cz(m, offset[k]): classically conditional CZ, but offset[k] is
        // classical. So this is a phase flip if both m=1 and offset[k]=1.
        // We implement as: if bit(k): Z_if(???, m) - but CZ on a classical value is...
        // Actually, `cz(m, offset[k])` means CZ conditional on classical m AND classical offset[k].
        // If either is 0 classically, no-op. If both 1, apply Z to... nothing?
        // Wait - `cz` in the CircuitBuilder takes two args. When one is a classical bit,
        // it's a phase flip conditional on that bit. Here `m` is a Bit and offset[k] is a Bit.
        // If both are classical bits, cz(m, bk) = apply neg if both are 1.
        // In our framework: neg_if(m) if bit(k) is 1 (classical).
        if bit(k) {
            let mut op = Op::empty();
            op.kind = OperationType::Neg;
            op.c_condition = m;
            b.ops.push(op);
        }
        // ccz(m, q_target[k], carries[k] XOR offset[k]).
        // This is CZ(q_target[k], carries[k] with inv based on offset[k])
        // classically conditioned on m.
        if k == 0 {
            // carries[0] = cin. Classical. cin XOR offset[0] = bool.
            let eff = carry_in ^ bit(0);
            if eff {
                // ccz(m, q_target[k], 1) = cz(m, q_target[k]) = z_if(q_target[k], m)?
                // Actually ccz(m, q, 1) applies negative phase iff m=1 AND q=1 AND 1=1.
                // That's just z_if(q, m).
                let mut op = Op::empty();
                op.kind = OperationType::Z;
                op.q_target = q_target[k];
                op.c_condition = m;
                b.ops.push(op);
            }
        } else {
            let cur = get_carry(k).expect("k in bounds");
            // CCZ(q_target[k], cur, ???, m). We need a third qubit; but
            // Gidney's ccz was a 2-qubit Z (CZ with classical cond). Our
            // ccz_if takes 3 qubits. Since we only want CZ on (q_target, cur)
            // conditioned on m, and Neg op is global phase flip on m, we use
            // `cz_if(q_target[k], cur, m)` instead.
            if bit(k) {
                b.x(cur);
                b.cz_if(q_target[k], cur, m);
                b.x(cur);
            } else {
                b.cz_if(q_target[k], cur, m);
            }
        }
    }
    // cx(carries[0], q_target[0]): if cin: X(q_target[0]).
    if carry_in {
        b.x(q_target[0]);
    }
}

/// Gidney 2025 adder with 2 clean + (n-2) dirty ancilla (Figure 4).
/// Performs `Q_target += offset + carry_in` using 3n ± O(1) CCX.
///
/// Uses the vented 2-clean adder then corrects via a pair of carry-xors
/// sandwiching classically-controlled Z gates (to convert vent bits into
/// actual phase flips).
///
/// **STATUS**: initial port but correctness is INCOMPLETE. The Python
/// reference merges the carry-xor into the vented add via
/// `carry_xor_target=[None]+Q_dirty`; our port does them separately,
/// which produces correct sum in q_target but LEAKS PHASE and perturbs
/// q_dirty. Needs: (a) extend add_vented_2clean_classical with a
/// `carry_xor_target` parameter, OR (b) figure out the correct
/// sequencing of carry-xor + vent-key phase-fix.
///
/// # Arguments
/// - `q_target`: n qubits (destination).
/// - `q_dirty`: at least n-2 dirty ancilla qubits (value preserved).
/// - `q_clean2`: at least 2 clean ancilla.
/// - `offset_bits`: classical offset.
/// - `carry_in`: classical carry-in.
#[allow(dead_code)]
pub(crate) fn iadd_dirty_2clean_classical(
    b: &mut B,
    q_target: &[QubitId],
    q_dirty: &[QubitId],
    q_clean2: &[QubitId; 2],
    offset_bits: u64,
    carry_in: bool,
) {
    let n = q_target.len();
    if n == 0 {
        return;
    }
    // Fall back to HRS linear-clean if we have enough clean qubits.
    // (Here we only have 2 clean. HRS needs n-2. If n<=4, q_clean2 suffices.)
    if n <= 4 {
        iadd_linear_clean_classical(b, q_target, q_clean2, offset_bits, carry_in);
        return;
    }
    assert!(q_dirty.len() >= n - 2, "need n-2 dirty qubits");
    let q_dirty = &q_dirty[..n - 2];

    // Vent_keys: n classical bits.
    let vent_keys: Vec<BitId> = (0..n).map(|_| b.alloc_bit()).collect();

    // carry_xor_target matches Python's [None] + Q_dirty (length n). At step
    // k (for k >= 1), XOR carries[k] into q_dirty[k-1].
    let cxt: Vec<Option<QubitId>> = (0..n)
        .map(|k| {
            if k == 0 {
                None
            } else {
                q_dirty.get(k - 1).copied()
            }
        })
        .collect();

    // Run the vented 2-clean adder WITH carry_xor_target merged.
    add_vented_2clean_classical_cxt(
        b,
        q_target,
        q_clean2,
        offset_bits,
        carry_in,
        &vent_keys,
        Some(&cxt),
    );

    // Broadcast_x on q_target (NOT each bit).
    for k in 0..n {
        b.x(q_target[k]);
    }
    // Broadcast_cz(q_dirty, vent_keys[1:]): for k in 0..n-2, z_if(q_dirty[k], vent_keys[k+1]).
    for k in 0..n - 2 {
        let mut op = Op::empty();
        op.kind = OperationType::Z;
        op.q_target = q_dirty[k];
        op.c_condition = vent_keys[k + 1];
        b.ops.push(op);
    }
    // carry_xor into q_dirty (src is now the bit-inverted q_target, which by
    // Gidney eq. 8 produces the same carries as the original pre-add target).
    // Python: Q_src=Q_target[:-1] (n-1 bits), after broadcast_x. We're
    // already in the broadcast-x sandwich, so use q_target directly.
    xor_right_shifted_carries_into_classical(b, &q_target[..n - 1], offset_bits, q_dirty, carry_in);
    for k in 0..n - 2 {
        let mut op = Op::empty();
        op.kind = OperationType::Z;
        op.q_target = q_dirty[k];
        op.c_condition = vent_keys[k + 1];
        b.ops.push(op);
    }
    for k in 0..n {
        b.x(q_target[k]);
    }
}

/// Controlled variant of `iadd_dirty_2clean_classical`: performs
/// `if ctrl: q_target += offset + carry_in` using the Gidney replacement
/// rule "replace every offset bit that's 1 with the control qubit".
///
/// # Note
/// carry_in is assumed classical (not controlled). If you need the
/// carry_in to be conditional on ctrl too, pre-process it.
pub(crate) fn ciadd_dirty_2clean_classical(
    b: &mut B,
    q_target: &[QubitId],
    q_dirty: &[QubitId],
    q_clean2: &[QubitId; 2],
    offset_bits: u64,
    ctrl: QubitId,
    carry_in: bool,
) {
    // When ctrl=0, we want NO add at all. Classical carry_in is only
    // actually applied when ctrl=1. Effective carry_in = ctrl AND
    // classical_carry_in. Since classical_carry_in is a compile-time bool,
    // when it's true we need carry_in = ctrl (quantum); when false, 0.
    // The rest of ciadd_dirty_2clean passes `carry_in: bool` = classical.
    // Work around by transforming: if carry_in=true, we effectively want
    // the adder to add (offset + 1) when ctrl=1. But offset+1 might change
    // many bits of offset (carry chain). Simpler: if carry_in=true, we
    // temporarily set q_target[0] ^= ctrl, then run the add with cin=false,
    // then... hmm this changes the add's trajectory.
    //
    // Cleanest fix: support `carry_in_q: Option<QubitId>` where Some(ctrl)
    // means the carry-in is a qubit. For now, require caller to pass
    // carry_in=false when using the controlled variant.
    assert!(
        !carry_in,
        "ciadd_dirty_2clean_classical requires carry_in=false; pre-process if needed"
    );
    let n = q_target.len();
    if n == 0 {
        return;
    }
    if n <= 4 {
        // Fallback: use HRS variant. For simplicity, apply X gates controlled
        // on ctrl to simulate controlled-add by CX-loading `offset` into a
        // temp n-bit register (this defeats the ancilla-saving purpose for
        // small n but is correct).
        let a: Vec<QubitId> = (0..n).map(|_| b.alloc_qubit()).collect();
        for i in 0..n {
            if (offset_bits >> i) & 1 != 0 {
                b.cx(ctrl, a[i]);
            }
        }
        // Use HRS linear-clean with a and q_target; treat a as the offset via
        // a CX-loaded classical constant. But HRS takes CLASSICAL offset.
        // So we'd need to do a quantum-quantum add here. Simpler: just do it
        // as we already do (via ccx to load f).
        // Actually our q_clean2 has 2 clean. For n<=4 we need n-2<=2 clean
        // which HRS supports. But HRS needs CLASSICAL offset; here offset is
        // quantum (a). Different primitive needed.
        //
        // For now: just bail out and use the caller's existing code path.
        // We'll skip this branch by asserting n>4.
        for i in 0..n {
            if (offset_bits >> i) & 1 != 0 {
                b.cx(ctrl, a[i]);
            }
        }
        for q in a {
            b.free(q);
        }
        panic!("ciadd_dirty_2clean: n<=4 fallback not implemented; use uncontrolled path");
    }
    assert!(q_dirty.len() >= n - 2, "need n-2 dirty qubits");
    let q_dirty = &q_dirty[..n - 2];

    // Vent_keys: n classical bits.
    let vent_keys: Vec<BitId> = (0..n).map(|_| b.alloc_bit()).collect();

    // carry_xor_target (Python's [None] + Q_dirty).
    let cxt: Vec<Option<QubitId>> = (0..n)
        .map(|k| {
            if k == 0 {
                None
            } else {
                q_dirty.get(k - 1).copied()
            }
        })
        .collect();

    // Controlled vented add. When offset_bits[k] = 1, the operations that
    // would have unconditionally used `1` now use ctrl.
    // The add_vented_2clean_classical_cxt takes a classical offset, so we
    // can't directly use it here. Write an inline controlled version.
    c_add_vented_2clean_inline(
        b,
        q_target,
        q_clean2,
        offset_bits,
        ctrl,
        carry_in,
        &vent_keys,
        &cxt,
    );

    for k in 0..n {
        // Replace broadcast_x with controlled X.
        b.cx(ctrl, q_target[k]);
    }
    for k in 0..n - 2 {
        // Z on q_dirty[k] conditional on vent_keys[k+1].
        // But Gidney's controlled variant: Z should also be controlled by ctrl.
        // Actually no: the phase fix is wrt the ACTUAL vent measurements,
        // which already include ctrl via the vented add. So Z is just
        // applied iff vent_keys[k+1]=1 (classical).
        let mut op = Op::empty();
        op.kind = OperationType::Z;
        op.q_target = q_dirty[k];
        op.c_condition = vent_keys[k + 1];
        b.ops.push(op);
    }
    // The carry_xor should also be controlled. For simplicity, fall back:
    // use a controlled version of xor_right_shifted_carries_into.
    c_xor_right_shifted_carries_into_classical(
        b,
        &q_target[..n - 1],
        offset_bits,
        ctrl,
        q_dirty,
        carry_in,
    );
    for k in 0..n - 2 {
        let mut op = Op::empty();
        op.kind = OperationType::Z;
        op.q_target = q_dirty[k];
        op.c_condition = vent_keys[k + 1];
        b.ops.push(op);
    }
    for k in 0..n {
        b.cx(ctrl, q_target[k]);
    }
}

/// Controlled vented add (inline). Matches `add_vented_2clean_classical_cxt`
/// but with each offset_bits[k]=1 behaving as if controlled by `ctrl`.
fn c_add_vented_2clean_inline(
    b: &mut B,
    q_target: &[QubitId],
    q_clean2: &[QubitId; 2],
    offset_bits: u64,
    ctrl: QubitId,
    carry_in: bool,
    vent_keys: &[BitId],
    carry_xor_target: &[Option<QubitId>],
) {
    let n = q_target.len();
    if n < 2 {
        // Degenerate case: for n=1, just do CCX(ctrl, offset[0] == 1, q_target[0]).
        if n == 1 {
            if carry_in {
                b.cx(ctrl, q_target[0]);
            }
            if (offset_bits & 1) != 0 {
                b.cx(ctrl, q_target[0]);
            }
        }
        return;
    }

    let bit = |k: usize| -> bool {
        if k >= 64 {
            false
        } else {
            (offset_bits >> k) & 1 != 0
        }
    };
    // broadcast_cx(offset, q_target) becomes: for k, if offset[k]=1: CX(ctrl, q_target[k]).
    for k in 0..n {
        if bit(k) {
            b.cx(ctrl, q_target[k]);
        }
    }
    // Helpers
    let get_carry_qubit = |k: usize| -> Option<QubitId> {
        if k == 0 {
            None
        } else if k == n - 1 {
            Some(q_target[n - 1])
        } else {
            Some(q_clean2[k % 2])
        }
    };

    for k in 0..n - 1 {
        // Reset next carry (if it's a clean ancilla).
        if k < n - 2 {
            if let Some(q) = get_carry_qubit(k + 1) {
                let mut op = Op::empty();
                op.kind = OperationType::R;
                op.q_target = q;
                b.ops.push(op);
            }
        }

        // CCX(q_target[k], carries[k] XOR (ctrl * offset[k]), carries[k+1])
        // For k=0: carries[0] = cin (classical).
        //   carries[0] XOR (ctrl * offset[0]) = cin XOR (ctrl AND bit(0)).
        //   If bit(0)=1: = cin XOR ctrl (= ~cin if ctrl=1, cin if ctrl=0).
        //   If bit(0)=0: = cin (classical).
        // The CCX's three inputs: q_target[k], the above, carries[k+1].
        // Use classical carry_in -> either trivial or becomes a CCX with ctrl.
        if k == 0 {
            let next = get_carry_qubit(1);
            if let Some(next_q) = next {
                if bit(0) {
                    // CCX(q_target[0], cin XOR ctrl, next_q). Use CX if cin=1
                    // (inverts ctrl control) and CCX otherwise.
                    if carry_in {
                        // cin XOR ctrl = NOT ctrl.
                        b.x(ctrl);
                        b.ccx(q_target[0], ctrl, next_q);
                        b.x(ctrl);
                    } else {
                        b.ccx(q_target[0], ctrl, next_q);
                    }
                } else if carry_in {
                    // CCX(q_target[0], 1, next_q) = CX(q_target[0], next_q).
                    b.cx(q_target[0], next_q);
                }
                // else: both inputs 0, no op.
            }
        } else {
            let cur = get_carry_qubit(k).expect("non-boundary carry");
            let next = get_carry_qubit(k + 1).expect("non-boundary next carry");
            if bit(k) {
                // carries[k] XOR ctrl. Use CCCX-style decomp: flip cur via CX(ctrl, cur),
                // then CCX(q_target[k], cur, next), then flip back.
                b.cx(ctrl, cur);
                b.ccx(q_target[k], cur, next);
                b.cx(ctrl, cur);
            } else {
                b.ccx(q_target[k], cur, next);
            }
        }

        // CX(carries[k], q_target[k]).
        if k == 0 {
            if carry_in {
                b.x(q_target[0]);
            }
        } else {
            let cur = get_carry_qubit(k).expect("non-boundary carry");
            b.cx(cur, q_target[k]);
        }

        // Optional carry_xor_target.
        if k < carry_xor_target.len() {
            if let Some(dst) = carry_xor_target[k] {
                if k == 0 {
                    if carry_in {
                        b.x(dst);
                    }
                } else {
                    let cur = get_carry_qubit(k).expect("non-boundary carry");
                    b.cx(cur, dst);
                }
            }
        }

        // Measure vent.
        if k > 0 {
            let cur = get_carry_qubit(k).expect("non-boundary carry");
            b.hmr(cur, vent_keys[k]);
        }

        // CX(offset[k], carries[k+1]) becomes CX(ctrl, carries[k+1]) if offset[k]=1.
        if bit(k) {
            if let Some(q) = get_carry_qubit(k + 1) {
                b.cx(ctrl, q);
            }
        }
    }
}

// ============================================================================
// Quantum-offset variants (for use when the offset is a quantum register,
// not a classical constant). The Gidney replacement rule: where classical
// offset[k]=1 triggered an operation, quantum offset[k] now CONTROLS that
// operation.
// ============================================================================

/// Quantum-offset variant of `add_vented_2clean_classical_cxt`: performs
/// `q_target += q_offset + carry_in` (mod 2^n) where q_offset is quantum.
///
/// Cost: 2n±O(1) CCX, 2 clean ancilla, n classical vent_keys.
/// (vs. our Cuccaro-based add_nbit_qq_fast at n-1 CCX + n-1 carry ancilla.)
///
/// # Peak win
/// Peak transient during this add: 2 clean + 1 c_in = 3 extra qubits.
/// vs Cuccaro fast which needs n-1 carry ancilla = n+O(1) extra qubits.
/// Saves ~n qubits at peak.
pub(crate) fn add_vented_2clean_qoffset(
    b: &mut B,
    q_target: &[QubitId],
    q_clean2: &[QubitId; 2],
    q_offset: &[QubitId],
    carry_in: bool,
    vent_keys: &[BitId],
    carry_xor_target: Option<&[Option<QubitId>]>,
) {
    let n = q_target.len();
    assert_eq!(q_offset.len(), n, "q_offset length must match q_target");
    if n == 0 {
        return;
    }
    if n == 1 {
        if carry_in {
            b.x(q_target[0]);
        }
        b.cx(q_offset[0], q_target[0]);
        return;
    }

    // broadcast_cx(q_offset, q_target): CX(q_offset[k], q_target[k]).
    for k in 0..n {
        b.cx(q_offset[k], q_target[k]);
    }

    let get_carry_qubit = |k: usize| -> Option<QubitId> {
        if k == 0 {
            None
        } else if k == n - 1 {
            Some(q_target[n - 1])
        } else {
            Some(q_clean2[k % 2])
        }
    };

    for k in 0..n - 1 {
        if k < n - 2 {
            if let Some(q) = get_carry_qubit(k + 1) {
                let mut op = Op::empty();
                op.kind = OperationType::R;
                op.q_target = q;
                b.ops.push(op);
            }
        }

        // CCX(q_target[k], carries[k] XOR q_offset[k], carries[k+1])
        // For k=0: carries[0] = cin (classical). CCX(q_target[k], cin XOR q_offset[0], next).
        // If cin=0: CCX(q_target[0], q_offset[0], next).
        // If cin=1: CCX(q_target[0], NOT q_offset[0], next).
        if k == 0 {
            let next = get_carry_qubit(1);
            if let Some(next_q) = next {
                if carry_in {
                    b.x(q_offset[0]);
                    b.ccx(q_target[0], q_offset[0], next_q);
                    b.x(q_offset[0]);
                } else {
                    b.ccx(q_target[0], q_offset[0], next_q);
                }
            }
        } else {
            let cur = get_carry_qubit(k).expect("non-boundary carry");
            let next = get_carry_qubit(k + 1).expect("non-boundary next carry");
            // CCX(q_target[k], cur XOR q_offset[k], next).
            // Do: CX(q_offset[k], cur); CCX(q_target[k], cur, next); CX(q_offset[k], cur).
            b.cx(q_offset[k], cur);
            b.ccx(q_target[k], cur, next);
            b.cx(q_offset[k], cur);
        }

        // CX(carries[k], q_target[k])
        if k == 0 {
            if carry_in {
                b.x(q_target[0]);
            }
        } else {
            let cur = get_carry_qubit(k).expect("non-boundary carry");
            b.cx(cur, q_target[k]);
        }

        // Optional carry_xor_target
        if let Some(cxt) = carry_xor_target {
            if k < cxt.len() {
                if let Some(dst) = cxt[k] {
                    if k == 0 {
                        if carry_in {
                            b.x(dst);
                        }
                    } else {
                        let cur = get_carry_qubit(k).expect("non-boundary carry");
                        b.cx(cur, dst);
                    }
                }
            }
        }

        // Vent: mx(carries[k], vent_keys[k])
        if k > 0 {
            let cur = get_carry_qubit(k).expect("non-boundary carry");
            b.hmr(cur, vent_keys[k]);
        }

        // CX(q_offset[k], carries[k+1])
        if let Some(q) = get_carry_qubit(k + 1) {
            b.cx(q_offset[k], q);
        }
    }
}

/// Quantum-offset version of xor_right_shifted_carries_into.
/// `Q_dst ^= carry(Q_src, q_offset, carry_in) >> 1`.
pub(crate) fn xor_right_shifted_carries_into_qoffset(
    b: &mut B,
    q_src: &[QubitId],
    q_offset: &[QubitId],
    q_dst: &[QubitId],
    carry_in: bool,
) {
    let n = q_dst.len();
    assert!(n <= q_src.len() && q_src.len() <= n + 1, "len mismatch");
    if n == 0 {
        return;
    }
    // Helper to apply CCX(src[k] XOR q_offset[k], dst_prev XOR q_offset[k], dst[k]).
    // We do this by CX(q_offset[k], src[k]); CX(q_offset[k], dst_prev); CCX; CX; CX.
    let ccx_with_qxor = |b: &mut B,
                         ctrl_a: QubitId,
                         xor_a: Option<QubitId>,
                         ctrl_b: QubitId,
                         xor_b: Option<QubitId>,
                         target: QubitId| {
        if let Some(x) = xor_a {
            b.cx(x, ctrl_a);
        }
        if let Some(x) = xor_b {
            b.cx(x, ctrl_b);
        }
        b.ccx(ctrl_a, ctrl_b, target);
        if let Some(x) = xor_b {
            b.cx(x, ctrl_b);
        }
        if let Some(x) = xor_a {
            b.cx(x, ctrl_a);
        }
    };

    for k in (1..n).rev() {
        ccx_with_qxor(b, q_src[k], Some(q_offset[k]), q_dst[k - 1], None, q_dst[k]);
    }
    // broadcast_cx(q_offset, q_dst): CX(q_offset[k], q_dst[k]).
    for k in 0..n {
        b.cx(q_offset[k], q_dst[k]);
    }
    // ccx(q_src[0] XOR q_offset[0], cin XOR q_offset[0], q_dst[0]).
    // For classical cin: if cin=1, the second control is NOT q_offset[0].
    b.cx(q_offset[0], q_src[0]);
    if carry_in {
        b.x(q_offset[0]);
    }
    b.ccx(q_src[0], q_offset[0], q_dst[0]);
    if carry_in {
        b.x(q_offset[0]);
    }
    b.cx(q_offset[0], q_src[0]);

    for k in 1..n {
        ccx_with_qxor(
            b,
            q_src[k],
            Some(q_offset[k]),
            q_dst[k - 1],
            Some(q_offset[k]),
            q_dst[k],
        );
    }
}

/// Quantum-offset version of iadd_dirty_2clean: `q_target += q_offset + cin`
/// using 2 clean + n-2 dirty ancilla. Cost ~3n CCX.
pub(crate) fn iadd_dirty_2clean_qoffset(
    b: &mut B,
    q_target: &[QubitId],
    q_dirty: &[QubitId],
    q_clean2: &[QubitId; 2],
    q_offset: &[QubitId],
    carry_in: bool,
) {
    let n = q_target.len();
    assert_eq!(q_offset.len(), n);
    if n == 0 {
        return;
    }
    if n <= 4 {
        panic!("iadd_dirty_2clean_qoffset: n<=4 not supported yet, use cuccaro_add");
    }
    assert!(q_dirty.len() >= n - 2, "need n-2 dirty qubits");
    let q_dirty = &q_dirty[..n - 2];

    let vent_keys: Vec<BitId> = (0..n).map(|_| b.alloc_bit()).collect();
    let cxt: Vec<Option<QubitId>> = (0..n)
        .map(|k| {
            if k == 0 {
                None
            } else {
                q_dirty.get(k - 1).copied()
            }
        })
        .collect();

    add_vented_2clean_qoffset(
        b,
        q_target,
        q_clean2,
        q_offset,
        carry_in,
        &vent_keys,
        Some(&cxt),
    );

    for k in 0..n {
        b.x(q_target[k]);
    }
    for k in 0..n - 2 {
        let mut op = Op::empty();
        op.kind = OperationType::Z;
        op.q_target = q_dirty[k];
        op.c_condition = vent_keys[k + 1];
        b.ops.push(op);
    }
    xor_right_shifted_carries_into_qoffset(b, &q_target[..n - 1], q_offset, q_dirty, carry_in);
    for k in 0..n - 2 {
        let mut op = Op::empty();
        op.kind = OperationType::Z;
        op.q_target = q_dirty[k];
        op.c_condition = vent_keys[k + 1];
        b.ops.push(op);
    }
    for k in 0..n {
        b.x(q_target[k]);
    }
}














pub(crate) fn isub_dirty_2clean_qoffset(
    b: &mut B,
    q_target: &[QubitId],
    q_dirty: &[QubitId],
    q_clean2: &[QubitId; 2],
    q_offset: &[QubitId],
) {
    let n = q_target.len();
    for k in 0..n {
        b.x(q_target[k]);
    }
    iadd_dirty_2clean_qoffset(b, q_target, q_dirty, q_clean2, q_offset, false);
    for k in 0..n {
        b.x(q_target[k]);
    }
}

/// Controlled variant of xor_right_shifted_carries_into.
fn c_xor_right_shifted_carries_into_classical(
    b: &mut B,
    q_src: &[QubitId],
    offset_bits: u64,
    ctrl: QubitId,
    q_dst: &[QubitId],
    carry_in: bool,
) {
    let n = q_dst.len();
    assert!(n <= q_src.len() && q_src.len() <= n + 1);
    if n == 0 {
        return;
    }
    let bit = |k: usize| -> bool {
        if k >= 64 {
            false
        } else {
            (offset_bits >> k) & 1 != 0
        }
    };

    // Helper for CCX where both controls may be "inverted" by XOR with ctrl.
    // The original has `Q_src[k] ^ offset[k]`; controlled version: if offset[k]=1,
    // the effective control is (Q_src[k] XOR ctrl); if offset[k]=0, it's just Q_src[k].
    let ccx_ctrl_mix = |b: &mut B,
                        ctrl_a: QubitId,
                        a_xor_ctrl: bool,
                        ctrl_b: QubitId,
                        b_xor_ctrl: bool,
                        target: QubitId| {
        if a_xor_ctrl {
            b.cx(ctrl, ctrl_a);
        }
        if b_xor_ctrl {
            b.cx(ctrl, ctrl_b);
        }
        b.ccx(ctrl_a, ctrl_b, target);
        if b_xor_ctrl {
            b.cx(ctrl, ctrl_b);
        }
        if a_xor_ctrl {
            b.cx(ctrl, ctrl_a);
        }
    };

    for k in (1..n).rev() {
        ccx_ctrl_mix(b, q_src[k], bit(k), q_dst[k - 1], false, q_dst[k]);
    }
    // broadcast_cx(offset, q_dst): for k, if offset[k]: CX(ctrl, q_dst[k]).
    for k in 0..n {
        if bit(k) {
            b.cx(ctrl, q_dst[k]);
        }
    }
    // ccx(q_src[0] XOR offset[0], carry_in XOR offset[0], q_dst[0])
    // carry_in XOR ctrl*offset[0]: if offset[0]=0 then just cin; if offset[0]=1 then cin XOR ctrl.
    let cin_eff_uses_ctrl = bit(0);
    let cin_classical_part = carry_in ^ false; // base carry_in, ctrl XOR handled separately
    if cin_eff_uses_ctrl {
        // Effective second control = ctrl XOR carry_in.
        // CCX(q_src[0] XOR (ctrl*offset[0]=ctrl), ctrl XOR cin, q_dst[0]).
        // We do this by: first adjusting q_src[0] based on ctrl (if bit(0)=1),
        // then the effective control is q_src[0]_adj AND (ctrl_XOR_cin).
        // Simpler: handle as CCX with cur=ctrl (since bit(0)=1) and
        // effective 2nd = ctrl XOR cin = ~ctrl if cin=1, else ctrl.
        // If cin=1: CCX(q_src[0] XOR ctrl, ~ctrl, q_dst[0]) = ...
        //   = ccx with both controls on ctrl in some form.
        // This is getting complex. Let's just compute the effective controls inline.
        if carry_in {
            // CCX(q_src[0] XOR ctrl, ~ctrl, q_dst[0]):
            //   flip q_src[0] via CX(ctrl, q_src[0]); flip ctrl via X; CCX; flip back
            b.cx(ctrl, q_src[0]);
            b.x(ctrl);
            b.ccx(q_src[0], ctrl, q_dst[0]);
            b.x(ctrl);
            b.cx(ctrl, q_src[0]);
        } else {
            // CCX(q_src[0] XOR ctrl, ctrl, q_dst[0]):
            b.cx(ctrl, q_src[0]);
            b.ccx(q_src[0], ctrl, q_dst[0]);
            b.cx(ctrl, q_src[0]);
        }
    } else {
        // offset[0]=0. CCX(q_src[0], cin, q_dst[0]).
        if cin_classical_part {
            // CCX(q_src[0], 1, q_dst[0]) = CX(q_src[0], q_dst[0]).
            b.cx(q_src[0], q_dst[0]);
        }
        // else both classical 0, no-op.
    }
    for k in 1..n {
        ccx_ctrl_mix(b, q_src[k], bit(k), q_dst[k - 1], bit(k), q_dst[k]);
    }
}

/// Controlled sub by classical constant: `if ctrl: q_target -= c` using
/// the identity `x - c = ~(~x + c)` and the venting `ciadd_dirty_2clean`.
///
/// Requires 2 clean + n-2 dirty ancilla. Cost: ~3n CCX + 2n CX.
pub(crate) fn cisub_dirty_2clean_classical(
    b: &mut B,
    q_target: &[QubitId],
    q_dirty: &[QubitId],
    q_clean2: &[QubitId; 2],
    c_bits: u64,
    ctrl: QubitId,
) {
    let n = q_target.len();
    // if ctrl: x = ~x
    for k in 0..n {
        b.cx(ctrl, q_target[k]);
    }
    ciadd_dirty_2clean_classical(
        b, q_target, q_dirty, q_clean2, c_bits, ctrl,
        false, // carry_in=false (controlled variant requires this)
    );
    for k in 0..n {
        b.cx(ctrl, q_target[k]);
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::sim::Simulator;
    use sha3::{
        digest::{ExtendableOutput, Update},
        Shake256,
    };

    fn anf_degree_density_from_truth_table(mut table: Vec<u8>, vars: usize) -> (usize, usize) {
        let states = 1usize << vars;
        // Möbius transform from truth table to ANF coefficients.
        for bit in 0..vars {
            let step = 1usize << bit;
            for mask in 0..states {
                if (mask & step) != 0 {
                    table[mask] ^= table[mask ^ step];
                }
            }
        }

        let mut degree = 0usize;
        let mut density = 0usize;
        for (mask, &coeff) in table.iter().enumerate() {
            if coeff != 0 {
                density += 1;
                degree = degree.max(mask.count_ones() as usize);
            }
        }
        (degree, density)
    }

    fn product_phase_anf_degree_density(n: usize, phase_mask: u64) -> (usize, usize) {
        assert!(n > 0 && n <= 10, "test keeps exhaustive table small");
        let vars = 2 * n;
        let states = 1usize << vars;
        let x_mask = (1u64 << n) - 1;
        let mut table = vec![0u8; states];
        for state in 0..states {
            let x = (state as u64) & x_mask;
            let y = ((state as u64) >> n) & x_mask;
            let prod = x * y;
            table[state] = ((prod & phase_mask).count_ones() & 1) as u8;
        }
        anf_degree_density_from_truth_table(table, vars)
    }

    fn carry_save_product_bits_for_phase_test(n: usize, x: u64, y: u64) -> Vec<u8> {
        // Deterministic carry-save compression of the n×n partial products down
        // to at most two wires per weight column.  This models the most tempting
        // redundant-product MBUC rescue: avoid a final carry-propagate product,
        // then X-measure the two carry-save rows instead of the binary product.
        let mut cols = vec![Vec::<u8>::new(); 2 * n + 8];
        for i in 0..n {
            for j in 0..n {
                let bit = (((x >> i) & 1) & ((y >> j) & 1)) as u8;
                cols[i + j].push(bit);
            }
        }
        for k in 0..cols.len() - 1 {
            while cols[k].len() > 2 {
                let a = cols[k].pop().unwrap();
                let b = cols[k].pop().unwrap();
                let c = cols[k].pop().unwrap();
                let sum = a ^ b ^ c;
                let carry = (a & b) ^ (a & c) ^ (b & c);
                cols[k].push(sum);
                cols[k + 1].push(carry);
            }
        }
        let mut out = Vec::with_capacity(4 * n + 4);
        for col in cols.iter().take(2 * n + 2) {
            out.push(*col.get(0).unwrap_or(&0));
            out.push(*col.get(1).unwrap_or(&0));
        }
        out
    }

    fn carry_save_product_phase_anf_degree_density(
        n: usize,
        top_column_only: bool,
    ) -> (usize, usize) {
        assert!(
            n > 0 && n <= 8,
            "test keeps exhaustive carry-save table small"
        );
        let vars = 2 * n;
        let states = 1usize << vars;
        let x_mask = (1u64 << n) - 1;
        let mut table = vec![0u8; states];
        for state in 0..states {
            let x = (state as u64) & x_mask;
            let y = ((state as u64) >> n) & x_mask;
            let bits = carry_save_product_bits_for_phase_test(n, x, y);
            table[state] = if top_column_only {
                let k = 2 * (2 * n - 2);
                bits[k] ^ bits[k + 1]
            } else {
                bits.iter().fold(0u8, |acc, &b| acc ^ b)
            };
        }
        anf_degree_density_from_truth_table(table, vars)
    }

    #[test]
    fn raw_product_measurement_phase_is_dense_not_free_kickmix() {
        // If a 2n-bit schoolbook product scratch `t=x*y` were simply X-measured,
        // the random measurement outcomes request phases of the form
        //     (-1)^(mask · (x*y))
        // on the preserved x/y registers.  The low product bit is quadratic, but
        // typical masks also touch carry-dependent high bits.  Exhaustive ANF on
        // toy widths shows these phase functions are already high-degree and
        // dense, so raw product-scratch MBUC is not the missing cheap IMUL.
        for &n in &[4usize, 6, 8, 10] {
            let full_mask = if 2 * n == 64 {
                u64::MAX
            } else {
                (1u64 << (2 * n)) - 1
            };
            let high_mask = 1u64 << (2 * n - 2);
            let (deg_full, dens_full) = product_phase_anf_degree_density(n, full_mask);
            let (deg_high, dens_high) = product_phase_anf_degree_density(n, high_mask);
            eprintln!(
                "raw_product_phase n={n} full_deg={deg_full} full_density={dens_full} high_deg={deg_high} high_density={dens_high}"
            );
            if n == 10 {
                println!("METRIC raw_product_mbu_fullmask_degree_n10={deg_full}");
                println!("METRIC raw_product_mbu_fullmask_density_n10={dens_full}");
                println!("METRIC raw_product_mbu_highbit_degree_n10={deg_high}");
                println!("METRIC raw_product_mbu_highbit_density_n10={dens_high}");
            }
        }

        let (deg_full, dens_full) = product_phase_anf_degree_density(10, (1u64 << 20) - 1);
        let (deg_high, dens_high) = product_phase_anf_degree_density(10, 1u64 << 18);
        assert_eq!(deg_full, 19);
        assert_eq!(dens_full, 427_812);
        assert_eq!(deg_high, 19);
        assert_eq!(dens_high, 120_581);
    }

    #[test]
    fn carry_save_product_scratch_mbu_still_has_dense_phases() {
        // Maybe the raw binary product was the wrong representation: a
        // carry-save product avoids the final carry-propagation chain.  But the
        // carry-save compressor still contains majority carries, and measuring
        // the final redundant rows asks for phases of those carry functions.
        // Exhaustive toy ANFs are already full-degree at n=8.
        for &n in &[4usize, 6, 8] {
            let (deg_all, dens_all) = carry_save_product_phase_anf_degree_density(n, false);
            let (deg_top, dens_top) = carry_save_product_phase_anf_degree_density(n, true);
            eprintln!(
                "carry_save_product_phase n={n} all_deg={deg_all} all_density={dens_all} top_deg={deg_top} top_density={dens_top}"
            );
            if n == 8 {
                println!("METRIC carry_save_product_mbu_all_degree_n8={deg_all}");
                println!("METRIC carry_save_product_mbu_all_density_n8={dens_all}");
                println!("METRIC carry_save_product_mbu_top_degree_n8={deg_top}");
                println!("METRIC carry_save_product_mbu_top_density_n8={dens_top}");
            }
        }
        let (deg_all, dens_all) = carry_save_product_phase_anf_degree_density(8, false);
        let (deg_top, dens_top) = carry_save_product_phase_anf_degree_density(8, true);
        assert_eq!(deg_all, 16);
        assert_eq!(dens_all, 20_440);
        assert_eq!(deg_top, 15);
        assert_eq!(dens_top, 3_602);
    }

    /// Classical reference: compute bit-k of carry(x, d, cin).
    /// The carry bit into position k (c_k) is defined by:
    ///   c_0 = cin
    ///   c_{k+1} = MAJ(c_k, x_k, d_k)
    fn classical_carry(x: u64, d: u64, cin: bool, n: usize) -> u64 {
        // Compute bit-by-bit.
        let mut c: u64 = 0;
        let mut prev = cin;
        for k in 0..n {
            let xk = (x >> k) & 1 != 0;
            let dk = (d >> k) & 1 != 0;
            // new carry = MAJ(prev, xk, dk)
            let new_carry = (prev && xk) || (prev && dk) || (xk && dk);
            if new_carry {
                c |= 1 << (k + 1);
            }
            prev = new_carry;
        }
        // Also set bit 0 to cin (the "carry into bit 0")
        if cin {
            c |= 1;
        }
        c
    }

    fn run_xor_rsh_carries(n: usize, trials: usize) -> bool {
        let mut hasher = Shake256::default();
        hasher.update(&[n as u8, trials as u8, 42]);
        use sha3::digest::XofReader;
        let mut xof = <sha3::Shake256 as sha3::digest::ExtendableOutput>::finalize_xof(hasher);
        for _trial in 0..trials {
            let mut buf = [0u8; 32];
            xof.read(&mut buf);
            let src_raw = u64::from_le_bytes(buf[0..8].try_into().unwrap());
            let dst_raw = u64::from_le_bytes(buf[8..16].try_into().unwrap());
            let offset_raw = u64::from_le_bytes(buf[16..24].try_into().unwrap());
            let cin_raw = buf[24];
            let src = if n < 64 {
                src_raw & ((1u64 << n) - 1)
            } else {
                src_raw
            };
            let dst = if n < 64 {
                dst_raw & ((1u64 << n) - 1)
            } else {
                dst_raw
            };
            let offset = if n < 64 {
                offset_raw & ((1u64 << n) - 1)
            } else {
                offset_raw
            };
            let cin = (cin_raw & 1) != 0;

            // Build circuit with src, dst qubits.
            let mut bb = B::new();
            let q_src: Vec<QubitId> = bb.alloc_qubits(n);
            let q_dst: Vec<QubitId> = bb.alloc_qubits(n);

            xor_right_shifted_carries_into_classical(&mut bb, &q_src, offset, &q_dst, cin);

            let ops = bb.ops.clone();
            let num_qubits = bb.next_qubit as usize;
            let num_bits = 0usize;
            let mut inner_hasher = Shake256::default();
            inner_hasher.update(&[77u8]);
            let mut inner_xof =
                <sha3::Shake256 as sha3::digest::ExtendableOutput>::finalize_xof(inner_hasher);
            let mut sim = Simulator::new(num_qubits, num_bits, &mut inner_xof);
            sim.clear_for_shot();
            // Set src[k] = (src >> k) & 1 for shot 0.
            for k in 0..n {
                if (src >> k) & 1 != 0 {
                    *sim.qubit_mut(q_src[k]) = 1; // set bit for shot 0
                }
                if (dst >> k) & 1 != 0 {
                    *sim.qubit_mut(q_dst[k]) = 1;
                }
            }
            sim.apply(&ops);

            let expected_carries = classical_carry(src, offset, cin, n + 1);
            let expected_rsh = expected_carries >> 1; // carries shifted right by 1
            let expected_dst = (dst ^ expected_rsh) & ((1u64 << n) - 1);

            let mut got_dst: u64 = 0;
            for k in 0..n {
                if sim.qubit(q_dst[k]) & 1 != 0 {
                    got_dst |= 1 << k;
                }
            }
            if got_dst != expected_dst {
                eprintln!(
                    "n={} src={:#x} dst={:#x} offset={:#x} cin={} got={:#x} exp={:#x}",
                    n, src, dst, offset, cin, got_dst, expected_dst
                );
                return false;
            }
        }
        true
    }

    #[test]
    fn test_xor_rsh_carries_small() {
        for n in 1..=8 {
            assert!(run_xor_rsh_carries(n, 20), "failed at n={n}");
        }
    }

    /// Test the vented 2-clean adder followed by phase-correction.
    /// Full protocol (Figure 4 in Gidney paper):
    /// 1. Run vented add on q_target with 2 clean ancilla, collecting
    ///    vent_keys.
    /// 2. Apply correction: broadcast_x(q_dst_xor_target); broadcast_cz(workspace, vent_keys);
    ///    xor_right_shifted_carries_into(...); broadcast_cz; xor_right_shifted_carries_into;
    ///    broadcast_x.
    ///
    /// For this test we use a DIRECT approach: add completes, then we
    /// simulate and verify:
    ///   (a) q_target holds correct sum.
    ///   (b) With vent_keys' phase contributions, global_phase is consistent.
    fn run_vented_add_2clean(n: usize, trials: usize) -> (usize, usize) {
        let mut hasher = Shake256::default();
        hasher.update(&[n as u8, trials as u8, 51]);
        use sha3::digest::XofReader;
        let mut xof = <sha3::Shake256 as sha3::digest::ExtendableOutput>::finalize_xof(hasher);
        let mut ok = 0;
        let mut bad = 0;
        for _trial in 0..trials {
            let mut buf = [0u8; 24];
            xof.read(&mut buf);
            let target_raw = u64::from_le_bytes(buf[0..8].try_into().unwrap());
            let offset_raw = u64::from_le_bytes(buf[8..16].try_into().unwrap());
            let cin_raw = buf[16];
            let mask = if n < 64 { (1u64 << n) - 1 } else { u64::MAX };
            let target = target_raw & mask;
            let offset = offset_raw & mask;
            let cin = (cin_raw & 1) != 0;

            let mut bb = B::new();
            let q_target: Vec<QubitId> = bb.alloc_qubits(n);
            let q_clean2: [QubitId; 2] = [bb.alloc_qubit(), bb.alloc_qubit()];
            let vent_keys: Vec<BitId> = (0..n).map(|_| bb.alloc_bit()).collect();

            add_vented_2clean_classical(&mut bb, &q_target, &q_clean2, offset, cin, &vent_keys);

            let ops = bb.ops.clone();
            let num_qubits = bb.next_qubit as usize;
            let num_bits = bb.next_bit as usize;
            let mut inner_hasher = Shake256::default();
            inner_hasher.update(&[101u8]);
            let mut inner_xof =
                <sha3::Shake256 as sha3::digest::ExtendableOutput>::finalize_xof(inner_hasher);
            let mut sim = Simulator::new(num_qubits, num_bits, &mut inner_xof);
            sim.clear_for_shot();
            for k in 0..n {
                if (target >> k) & 1 != 0 {
                    *sim.qubit_mut(q_target[k]) = 1;
                }
            }
            sim.apply(&ops);

            let expected_sum = (target.wrapping_add(offset).wrapping_add(cin as u64)) & mask;
            let mut got: u64 = 0;
            for k in 0..n {
                if sim.qubit(q_target[k]) & 1 != 0 {
                    got |= 1 << k;
                }
            }
            if got == expected_sum {
                ok += 1;
            } else {
                bad += 1;
                if bad < 3 {
                    eprintln!(
                        "vented add FAIL n={} t={:#x} o={:#x} cin={} got={:#x} exp={:#x}",
                        n, target, offset, cin, got, expected_sum
                    );
                }
            }
        }
        (ok, bad)
    }

    #[test]
    fn test_vented_add_2clean_small() {
        for n in 2..=8 {
            let (ok, bad) = run_vented_add_2clean(n, 20);
            assert_eq!(bad, 0, "n={n}: {ok}/{} passed", ok + bad);
        }
    }

    fn run_linear_clean_add(n: usize, trials: usize) -> (usize, usize) {
        let mut hasher = Shake256::default();
        hasher.update(&[n as u8, trials as u8, 73]);
        use sha3::digest::XofReader;
        let mut xof = <sha3::Shake256 as sha3::digest::ExtendableOutput>::finalize_xof(hasher);
        let mut ok = 0;
        let mut bad = 0;
        for _trial in 0..trials {
            let mut buf = [0u8; 24];
            xof.read(&mut buf);
            let target_raw = u64::from_le_bytes(buf[0..8].try_into().unwrap());
            let offset_raw = u64::from_le_bytes(buf[8..16].try_into().unwrap());
            let cin_raw = buf[16];
            let mask = if n < 64 { (1u64 << n) - 1 } else { u64::MAX };
            let target = target_raw & mask;
            let offset = offset_raw & mask;
            let cin = (cin_raw & 1) != 0;

            let mut bb = B::new();
            let q_target: Vec<QubitId> = bb.alloc_qubits(n);
            let n_clean = n.saturating_sub(2).max(2);
            let q_clean: Vec<QubitId> = bb.alloc_qubits(n_clean);

            iadd_linear_clean_classical(&mut bb, &q_target, &q_clean, offset, cin);

            let ops = bb.ops.clone();
            let num_qubits = bb.next_qubit as usize;
            let num_bits = bb.next_bit as usize;
            let mut inner_hasher = Shake256::default();
            inner_hasher.update(&[151u8]);
            let mut inner_xof =
                <sha3::Shake256 as sha3::digest::ExtendableOutput>::finalize_xof(inner_hasher);
            let mut sim = Simulator::new(num_qubits, num_bits, &mut inner_xof);
            sim.clear_for_shot();
            for k in 0..n {
                if (target >> k) & 1 != 0 {
                    *sim.qubit_mut(q_target[k]) = 1;
                }
            }
            sim.apply(&ops);

            let expected_sum = (target.wrapping_add(offset).wrapping_add(cin as u64)) & mask;
            let mut got: u64 = 0;
            for k in 0..n {
                if sim.qubit(q_target[k]) & 1 != 0 {
                    got |= 1 << k;
                }
            }
            if got == expected_sum {
                ok += 1;
            } else {
                bad += 1;
                if bad < 3 {
                    eprintln!(
                        "HRS FAIL n={} t={:#x} o={:#x} cin={} got={:#x} exp={:#x}",
                        n, target, offset, cin, got, expected_sum
                    );
                }
            }
        }
        (ok, bad)
    }

    #[test]
    fn test_iadd_linear_clean_small() {
        for n in 1..=8 {
            let (ok, bad) = run_linear_clean_add(n, 20);
            assert_eq!(bad, 0, "n={n}: {ok}/{} passed", ok + bad);
        }
    }

    fn run_iadd_dirty_2clean(n: usize, trials: usize) -> (usize, usize) {
        let mut hasher = Shake256::default();
        hasher.update(&[n as u8, trials as u8, 97]);
        use sha3::digest::XofReader;
        let mut xof = <sha3::Shake256 as sha3::digest::ExtendableOutput>::finalize_xof(hasher);
        let mut ok = 0;
        let mut bad = 0;
        for _trial in 0..trials {
            let mut buf = [0u8; 32];
            xof.read(&mut buf);
            let target_raw = u64::from_le_bytes(buf[0..8].try_into().unwrap());
            let offset_raw = u64::from_le_bytes(buf[8..16].try_into().unwrap());
            let dirty_raw = u64::from_le_bytes(buf[16..24].try_into().unwrap());
            let cin_raw = buf[24];
            let mask = if n < 64 { (1u64 << n) - 1 } else { u64::MAX };
            let target = target_raw & mask;
            let offset = offset_raw & mask;
            let dirty_init = dirty_raw & mask;
            let cin = (cin_raw & 1) != 0;

            let mut bb = B::new();
            let q_target: Vec<QubitId> = bb.alloc_qubits(n);
            let q_dirty: Vec<QubitId> = bb.alloc_qubits(n.saturating_sub(2).max(1));
            let q_clean2: [QubitId; 2] = [bb.alloc_qubit(), bb.alloc_qubit()];

            iadd_dirty_2clean_classical(&mut bb, &q_target, &q_dirty, &q_clean2, offset, cin);

            let ops = bb.ops.clone();
            let num_qubits = bb.next_qubit as usize;
            let num_bits = bb.next_bit as usize;
            let mut inner_hasher = Shake256::default();
            inner_hasher.update(&[201u8]);
            let mut inner_xof =
                <sha3::Shake256 as sha3::digest::ExtendableOutput>::finalize_xof(inner_hasher);
            let mut sim = Simulator::new(num_qubits, num_bits, &mut inner_xof);
            sim.clear_for_shot();
            for k in 0..n {
                if (target >> k) & 1 != 0 {
                    *sim.qubit_mut(q_target[k]) = 1;
                }
            }
            // Dirty init
            for (k, &q) in q_dirty.iter().enumerate() {
                if (dirty_init >> k) & 1 != 0 {
                    *sim.qubit_mut(q) = 1;
                }
            }
            sim.apply(&ops);

            let expected_sum = (target.wrapping_add(offset).wrapping_add(cin as u64)) & mask;
            let mut got: u64 = 0;
            for k in 0..n {
                if sim.qubit(q_target[k]) & 1 != 0 {
                    got |= 1 << k;
                }
            }
            // Check dirty is preserved (when n > 4, the dirty path is used).
            let mut got_dirty: u64 = 0;
            for (k, &q) in q_dirty.iter().enumerate() {
                if sim.qubit(q) & 1 != 0 {
                    got_dirty |= 1 << k;
                }
            }
            let dirty_ok = if n > 4 {
                got_dirty == (dirty_init & ((1u64 << q_dirty.len()) - 1).min(mask))
            } else {
                true
            };
            // Check phase is 0
            let phase = sim.global_phase() & 1;

            if got == expected_sum && dirty_ok && phase == 0 {
                ok += 1;
            } else {
                bad += 1;
                if bad < 3 {
                    eprintln!(
                        "iadd_dirty_2clean FAIL n={} t={:#x} o={:#x} d={:#x} cin={} got={:#x} exp={:#x} dirty_ok={} phase={}",
                        n, target, offset, dirty_init, cin, got, expected_sum, dirty_ok, phase
                    );
                }
            }
        }
        (ok, bad)
    }

    #[test]
    fn test_iadd_dirty_2clean_small() {
        for n in 2..=8 {
            let (ok, bad) = run_iadd_dirty_2clean(n, 10);
            assert_eq!(bad, 0, "n={n}: {ok}/{} passed", ok + bad);
        }
    }

    fn run_ciadd_dirty_2clean(n: usize, trials: usize) -> (usize, usize) {
        let mut hasher = Shake256::default();
        hasher.update(&[n as u8, trials as u8, 113]);
        use sha3::digest::XofReader;
        let mut xof = <sha3::Shake256 as sha3::digest::ExtendableOutput>::finalize_xof(hasher);
        let mut ok = 0;
        let mut bad = 0;
        for _trial in 0..trials {
            let mut buf = [0u8; 40];
            xof.read(&mut buf);
            let target_raw = u64::from_le_bytes(buf[0..8].try_into().unwrap());
            let offset_raw = u64::from_le_bytes(buf[8..16].try_into().unwrap());
            let dirty_raw = u64::from_le_bytes(buf[16..24].try_into().unwrap());
            let cin_raw = buf[24];
            let ctrl_raw = buf[25];
            let mask = if n < 64 { (1u64 << n) - 1 } else { u64::MAX };
            let target = target_raw & mask;
            let offset = offset_raw & mask;
            let dirty_init = dirty_raw & mask;
            let cin = false; // controlled variant requires classical cin=false
            let _ = cin_raw;
            let ctrl_val = (ctrl_raw & 1) != 0;

            let mut bb = B::new();
            let q_target: Vec<QubitId> = bb.alloc_qubits(n);
            let q_dirty: Vec<QubitId> = bb.alloc_qubits(n.saturating_sub(2).max(1));
            let q_clean2: [QubitId; 2] = [bb.alloc_qubit(), bb.alloc_qubit()];
            let q_ctrl = bb.alloc_qubit();

            ciadd_dirty_2clean_classical(
                &mut bb, &q_target, &q_dirty, &q_clean2, offset, q_ctrl, cin,
            );

            let ops = bb.ops.clone();
            let num_qubits = bb.next_qubit as usize;
            let num_bits = bb.next_bit as usize;
            let mut inner_hasher = Shake256::default();
            inner_hasher.update(&[211u8]);
            let mut inner_xof =
                <sha3::Shake256 as sha3::digest::ExtendableOutput>::finalize_xof(inner_hasher);
            let mut sim = Simulator::new(num_qubits, num_bits, &mut inner_xof);
            sim.clear_for_shot();
            for k in 0..n {
                if (target >> k) & 1 != 0 {
                    *sim.qubit_mut(q_target[k]) = 1;
                }
            }
            for (k, &q) in q_dirty.iter().enumerate() {
                if (dirty_init >> k) & 1 != 0 {
                    *sim.qubit_mut(q) = 1;
                }
            }
            if ctrl_val {
                *sim.qubit_mut(q_ctrl) = 1;
            }
            sim.apply(&ops);

            let expected_sum = if ctrl_val {
                (target.wrapping_add(offset).wrapping_add(cin as u64)) & mask
            } else {
                target
            };
            let mut got: u64 = 0;
            for k in 0..n {
                if sim.qubit(q_target[k]) & 1 != 0 {
                    got |= 1 << k;
                }
            }
            let mut got_dirty: u64 = 0;
            for (k, &q) in q_dirty.iter().enumerate() {
                if sim.qubit(q) & 1 != 0 {
                    got_dirty |= 1 << k;
                }
            }
            let dirty_ok = got_dirty == (dirty_init & ((1u64 << q_dirty.len()) - 1).min(mask));
            let phase = sim.global_phase() & 1;
            let ctrl_preserved = sim.qubit(q_ctrl) & 1 == (ctrl_val as u64);

            if got == expected_sum && dirty_ok && phase == 0 && ctrl_preserved {
                ok += 1;
            } else {
                bad += 1;
                if bad < 3 {
                    eprintln!(
                        "ciadd_dirty FAIL n={} t={:#x} o={:#x} d={:#x} cin={} ctrl={} got={:#x} exp={:#x} d_ok={} phase={} ctrl_preserved={}",
                        n, target, offset, dirty_init, cin, ctrl_val, got, expected_sum, dirty_ok, phase, ctrl_preserved
                    );
                }
            }
        }
        (ok, bad)
    }

    #[test]
    fn test_ciadd_dirty_2clean_small() {
        for n in 5..=10 {
            let (ok, bad) = run_ciadd_dirty_2clean(n, 10);
            assert_eq!(bad, 0, "n={n}: {ok}/{} passed", ok + bad);
        }
    }

    fn run_cisub_dirty(n: usize, trials: usize) -> (usize, usize) {
        let mut hasher = Shake256::default();
        hasher.update(&[n as u8, trials as u8, 179]);
        use sha3::digest::XofReader;
        let mut xof = <sha3::Shake256 as sha3::digest::ExtendableOutput>::finalize_xof(hasher);
        let mut ok = 0;
        let mut bad = 0;
        for _trial in 0..trials {
            let mut buf = [0u8; 40];
            xof.read(&mut buf);
            let target_raw = u64::from_le_bytes(buf[0..8].try_into().unwrap());
            let c_raw = u64::from_le_bytes(buf[8..16].try_into().unwrap());
            let dirty_raw = u64::from_le_bytes(buf[16..24].try_into().unwrap());
            let ctrl_raw = buf[25];
            let mask = if n < 64 { (1u64 << n) - 1 } else { u64::MAX };
            let target = target_raw & mask;
            let c = c_raw & mask;
            let dirty_init = dirty_raw & mask;
            let ctrl_val = (ctrl_raw & 1) != 0;

            let mut bb = B::new();
            let q_target: Vec<QubitId> = bb.alloc_qubits(n);
            let q_dirty: Vec<QubitId> = bb.alloc_qubits(n.saturating_sub(2).max(1));
            let q_clean2: [QubitId; 2] = [bb.alloc_qubit(), bb.alloc_qubit()];
            let q_ctrl = bb.alloc_qubit();

            cisub_dirty_2clean_classical(&mut bb, &q_target, &q_dirty, &q_clean2, c, q_ctrl);

            let ops = bb.ops.clone();
            let num_qubits = bb.next_qubit as usize;
            let num_bits = bb.next_bit as usize;
            let mut inner_hasher = Shake256::default();
            inner_hasher.update(&[221u8]);
            let mut inner_xof =
                <sha3::Shake256 as sha3::digest::ExtendableOutput>::finalize_xof(inner_hasher);
            let mut sim = Simulator::new(num_qubits, num_bits, &mut inner_xof);
            sim.clear_for_shot();
            for k in 0..n {
                if (target >> k) & 1 != 0 {
                    *sim.qubit_mut(q_target[k]) = 1;
                }
            }
            for (k, &q) in q_dirty.iter().enumerate() {
                if (dirty_init >> k) & 1 != 0 {
                    *sim.qubit_mut(q) = 1;
                }
            }
            if ctrl_val {
                *sim.qubit_mut(q_ctrl) = 1;
            }
            sim.apply(&ops);

            let expected = if ctrl_val {
                target.wrapping_sub(c) & mask
            } else {
                target
            };
            let mut got: u64 = 0;
            for k in 0..n {
                if sim.qubit(q_target[k]) & 1 != 0 {
                    got |= 1 << k;
                }
            }
            let mut got_dirty: u64 = 0;
            for (k, &q) in q_dirty.iter().enumerate() {
                if sim.qubit(q) & 1 != 0 {
                    got_dirty |= 1 << k;
                }
            }
            let dirty_ok = got_dirty == (dirty_init & ((1u64 << q_dirty.len()) - 1).min(mask));
            let phase = sim.global_phase() & 1;

            if got == expected && dirty_ok && phase == 0 {
                ok += 1;
            } else {
                bad += 1;
                if bad < 3 {
                    eprintln!(
                        "cisub FAIL n={} t={:#x} c={:#x} ctrl={} got={:#x} exp={:#x} d_ok={} phase={}",
                        n, target, c, ctrl_val, got, expected, dirty_ok, phase
                    );
                }
            }
        }
        (ok, bad)
    }

    #[test]
    fn test_cisub_dirty_small() {
        for n in 5..=10 {
            let (ok, bad) = run_cisub_dirty(n, 10);
            assert_eq!(bad, 0, "n={n}: {ok}/{} passed", ok + bad);
        }
    }

    #[test]
    fn test_cisub_dirty_large() {
        let n = 256;
        let mut hasher = Shake256::default();
        hasher.update(&[n as u8, 50u8, 17]);
        use sha3::digest::XofReader;
        let mut xof = <sha3::Shake256 as sha3::digest::ExtendableOutput>::finalize_xof(hasher);
        let trials = 50;
        let c_low = 0x1_0000_03D1u64;
        let mut ok = 0;
        let mut bad = 0;
        for _trial in 0..trials {
            let mut buf = [0u8; 40];
            xof.read(&mut buf);
            let target = u64::from_le_bytes(buf[0..8].try_into().unwrap());
            let dirty_init = u64::from_le_bytes(buf[8..16].try_into().unwrap());
            let ctrl_val = (buf[16] & 1) != 0;

            let mut bb = B::new();
            let q_target: Vec<QubitId> = bb.alloc_qubits(n);
            let q_dirty: Vec<QubitId> = bb.alloc_qubits(n - 2);
            let q_clean2: [QubitId; 2] = [bb.alloc_qubit(), bb.alloc_qubit()];
            let q_ctrl = bb.alloc_qubit();

            cisub_dirty_2clean_classical(&mut bb, &q_target, &q_dirty, &q_clean2, c_low, q_ctrl);

            let ops = bb.ops.clone();
            let num_qubits = bb.next_qubit as usize;
            let num_bits = bb.next_bit as usize;
            let mut inner_hasher = Shake256::default();
            inner_hasher.update(&[19u8]);
            let mut inner_xof =
                <sha3::Shake256 as sha3::digest::ExtendableOutput>::finalize_xof(inner_hasher);
            let mut sim = Simulator::new(num_qubits, num_bits, &mut inner_xof);
            sim.clear_for_shot();
            for k in 0..64 {
                if (target >> k) & 1 != 0 {
                    *sim.qubit_mut(q_target[k]) = 1;
                }
            }
            for (k, &q) in q_dirty.iter().enumerate().take(64) {
                if (dirty_init >> k) & 1 != 0 {
                    *sim.qubit_mut(q) = 1;
                }
            }
            if ctrl_val {
                *sim.qubit_mut(q_ctrl) = 1;
            }
            sim.apply(&ops);

            let expected = if ctrl_val {
                target.wrapping_sub(c_low)
            } else {
                target
            };
            let mut got: u64 = 0;
            for k in 0..64 {
                if sim.qubit(q_target[k]) & 1 != 0 {
                    got |= 1 << k;
                }
            }
            let mut got_dirty: u64 = 0;
            for (k, &q) in q_dirty.iter().enumerate().take(64) {
                if sim.qubit(q) & 1 != 0 {
                    got_dirty |= 1 << k;
                }
            }
            let dirty_ok = got_dirty == dirty_init;
            let phase = sim.global_phase() & 1;
            let ctrl_preserved = sim.qubit(q_ctrl) & 1 == (ctrl_val as u64);
            if got == expected && dirty_ok && phase == 0 && ctrl_preserved {
                ok += 1;
            } else {
                bad += 1;
                if bad < 3 {
                    eprintln!(
                        "cisub n=256 FAIL t={:#x} d={:#x} ctrl={} got={:#x} exp={:#x} d_ok={} phase={} ctrl_ok={}",
                        target, dirty_init, ctrl_val, got, expected, dirty_ok, phase, ctrl_preserved
                    );
                }
            }
        }
        assert_eq!(bad, 0, "n=256 cisub: {ok}/{trials} passed");
    }

    #[test]
    fn test_cisub_dirty_kaliski_pattern() {
        // Test with dirty qubits in Kaliski-specific patterns.
        let n = 256;
        let c_low = 0x1_0000_03D1u64;
        let trials = 50;
        let mut hasher = Shake256::default();
        hasher.update(&[99u8]);
        use sha3::digest::XofReader;
        let mut xof = <sha3::Shake256 as sha3::digest::ExtendableOutput>::finalize_xof(hasher);
        let mut ok = 0;
        let mut bad = 0;
        for _trial in 0..trials {
            let mut buf = [0u8; 16];
            xof.read(&mut buf);
            let target = u64::from_le_bytes(buf[0..8].try_into().unwrap());
            let dirty_u_lsb = (buf[8] & 1) != 0; // u[0] simulator
            let ctrl_val = (buf[9] & 1) != 0;

            let mut bb = B::new();
            let q_target: Vec<QubitId> = bb.alloc_qubits(n);
            let q_dirty: Vec<QubitId> = bb.alloc_qubits(n - 2);
            let q_clean2: [QubitId; 2] = [bb.alloc_qubit(), bb.alloc_qubit()];
            let q_ctrl = bb.alloc_qubit();

            cisub_dirty_2clean_classical(&mut bb, &q_target, &q_dirty, &q_clean2, c_low, q_ctrl);

            let ops = bb.ops.clone();
            let num_qubits = bb.next_qubit as usize;
            let num_bits = bb.next_bit as usize;
            let mut inner_hasher = Shake256::default();
            inner_hasher.update(&[77u8]);
            let mut inner_xof =
                <sha3::Shake256 as sha3::digest::ExtendableOutput>::finalize_xof(inner_hasher);
            let mut sim = Simulator::new(num_qubits, num_bits, &mut inner_xof);
            sim.clear_for_shot();
            for k in 0..64 {
                if (target >> k) & 1 != 0 {
                    *sim.qubit_mut(q_target[k]) = 1;
                }
            }
            // Kaliski pattern: dirty[0] = u[0]=1 (at termination). Rest = 0.
            if dirty_u_lsb {
                *sim.qubit_mut(q_dirty[0]) = 1;
            }
            if ctrl_val {
                *sim.qubit_mut(q_ctrl) = 1;
            }
            sim.apply(&ops);

            let expected = if ctrl_val {
                target.wrapping_sub(c_low)
            } else {
                target
            };
            let mut got: u64 = 0;
            for k in 0..64 {
                if sim.qubit(q_target[k]) & 1 != 0 {
                    got |= 1 << k;
                }
            }
            let got_dirty0 = sim.qubit(q_dirty[0]) & 1 != 0;
            let dirty_ok = got_dirty0 == dirty_u_lsb;
            let phase = sim.global_phase() & 1;
            let ctrl_preserved = sim.qubit(q_ctrl) & 1 == (ctrl_val as u64);
            if got == expected && dirty_ok && phase == 0 && ctrl_preserved {
                ok += 1;
            } else {
                bad += 1;
                if bad < 3 {
                    eprintln!(
                        "cisub kaliski FAIL t={:#x} d0={} ctrl={} got={:#x} exp={:#x} d_ok={} phase={} ctrl_ok={}",
                        target, dirty_u_lsb, ctrl_val, got, expected, dirty_ok, phase, ctrl_preserved
                    );
                }
            }
        }
        assert_eq!(bad, 0, "kaliski pattern cisub: {ok}/{trials} passed");
    }

    fn run_iadd_qoffset_dirty(n: usize, trials: usize) -> (usize, usize) {
        let mut hasher = Shake256::default();
        hasher.update(&[n as u8, trials as u8, 199]);
        use sha3::digest::XofReader;
        let mut xof = <sha3::Shake256 as sha3::digest::ExtendableOutput>::finalize_xof(hasher);
        let mut ok = 0;
        let mut bad = 0;
        for _trial in 0..trials {
            let mut buf = [0u8; 40];
            xof.read(&mut buf);
            let target_raw = u64::from_le_bytes(buf[0..8].try_into().unwrap());
            let offset_raw = u64::from_le_bytes(buf[8..16].try_into().unwrap());
            let dirty_raw = u64::from_le_bytes(buf[16..24].try_into().unwrap());
            let cin_raw = buf[24];
            let mask = if n < 64 { (1u64 << n) - 1 } else { u64::MAX };
            let target = target_raw & mask;
            let offset = offset_raw & mask;
            let dirty_init = dirty_raw & mask;
            let cin = (cin_raw & 1) != 0;

            let mut bb = B::new();
            let q_target: Vec<QubitId> = bb.alloc_qubits(n);
            let q_offset: Vec<QubitId> = bb.alloc_qubits(n);
            let q_dirty: Vec<QubitId> = bb.alloc_qubits(n.saturating_sub(2).max(1));
            let q_clean2: [QubitId; 2] = [bb.alloc_qubit(), bb.alloc_qubit()];

            iadd_dirty_2clean_qoffset(&mut bb, &q_target, &q_dirty, &q_clean2, &q_offset, cin);

            let ops = bb.ops.clone();
            let num_qubits = bb.next_qubit as usize;
            let num_bits = bb.next_bit as usize;
            let mut inner_hasher = Shake256::default();
            inner_hasher.update(&[231u8]);
            let mut inner_xof =
                <sha3::Shake256 as sha3::digest::ExtendableOutput>::finalize_xof(inner_hasher);
            let mut sim = Simulator::new(num_qubits, num_bits, &mut inner_xof);
            sim.clear_for_shot();
            for k in 0..n {
                if (target >> k) & 1 != 0 {
                    *sim.qubit_mut(q_target[k]) = 1;
                }
                if (offset >> k) & 1 != 0 {
                    *sim.qubit_mut(q_offset[k]) = 1;
                }
            }
            for (k, &q) in q_dirty.iter().enumerate() {
                if (dirty_init >> k) & 1 != 0 {
                    *sim.qubit_mut(q) = 1;
                }
            }
            sim.apply(&ops);

            let expected = (target.wrapping_add(offset).wrapping_add(cin as u64)) & mask;
            let mut got: u64 = 0;
            for k in 0..n {
                if sim.qubit(q_target[k]) & 1 != 0 {
                    got |= 1 << k;
                }
            }
            // Check q_offset preserved.
            let mut got_offset: u64 = 0;
            for k in 0..n {
                if sim.qubit(q_offset[k]) & 1 != 0 {
                    got_offset |= 1 << k;
                }
            }
            let mut got_dirty: u64 = 0;
            for (k, &q) in q_dirty.iter().enumerate() {
                if sim.qubit(q) & 1 != 0 {
                    got_dirty |= 1 << k;
                }
            }
            let dirty_ok = got_dirty == (dirty_init & ((1u64 << q_dirty.len()) - 1).min(mask));
            let offset_ok = got_offset == offset;
            let phase = sim.global_phase() & 1;

            if got == expected && dirty_ok && offset_ok && phase == 0 {
                ok += 1;
            } else {
                bad += 1;
                if bad < 3 {
                    eprintln!(
                        "iadd_qoffset FAIL n={} t={:#x} o={:#x} d={:#x} cin={} got={:#x} exp={:#x} d_ok={} o_ok={} phase={}",
                        n, target, offset, dirty_init, cin, got, expected, dirty_ok, offset_ok, phase
                    );
                }
            }
        }
        (ok, bad)
    }

    #[test]
    fn test_iadd_qoffset_dirty_small() {
        for n in 5..=10 {
            let (ok, bad) = run_iadd_qoffset_dirty(n, 10);
            assert_eq!(bad, 0, "n={n}: {ok}/{} passed", ok + bad);
        }
    }

    fn run_iadd_qoffset_narrow(n: usize, m: usize, trials: usize) -> (usize, usize) {
        let mut hasher = Shake256::default();
        hasher.update(&[n as u8, m as u8, trials as u8, 211]);
        use sha3::digest::XofReader;
        let mut xof = <sha3::Shake256 as sha3::digest::ExtendableOutput>::finalize_xof(hasher);
        let mut ok = 0;
        let mut bad = 0;
        let nmask = if n < 64 { (1u64 << n) - 1 } else { u64::MAX };
        let mmask = if m < 64 { (1u64 << m) - 1 } else { u64::MAX };
        for _trial in 0..trials {
            let mut buf = [0u8; 40];
            xof.read(&mut buf);
            let target = u64::from_le_bytes(buf[0..8].try_into().unwrap()) & nmask;
            let offset = u64::from_le_bytes(buf[8..16].try_into().unwrap()) & mmask;
            let dirty_init = u64::from_le_bytes(buf[16..24].try_into().unwrap()) & nmask;
            let cin = (buf[24] & 1) != 0;

            let mut bb = B::new();
            let q_target: Vec<QubitId> = bb.alloc_qubits(n);
            let q_offset: Vec<QubitId> = bb.alloc_qubits(m);
            let q_dirty: Vec<QubitId> = bb.alloc_qubits(n.saturating_sub(2).max(1));
            let q_clean2: [QubitId; 2] = [bb.alloc_qubit(), bb.alloc_qubit()];

            iadd_dirty_2clean_qoffset_narrow(
                &mut bb, &q_target, &q_dirty, &q_clean2, &q_offset, cin,
            );

            let ops = bb.ops.clone();
            let num_qubits = bb.next_qubit as usize;
            let num_bits = bb.next_bit as usize;
            let mut inner_hasher = Shake256::default();
            inner_hasher.update(&[233u8]);
            let mut inner_xof =
                <sha3::Shake256 as sha3::digest::ExtendableOutput>::finalize_xof(inner_hasher);
            let mut sim = Simulator::new(num_qubits, num_bits, &mut inner_xof);
            sim.clear_for_shot();
            for k in 0..n {
                if k < 64 && (target >> k) & 1 != 0 {
                    *sim.qubit_mut(q_target[k]) = 1;
                }
            }
            for k in 0..m {
                if k < 64 && (offset >> k) & 1 != 0 {
                    *sim.qubit_mut(q_offset[k]) = 1;
                }
            }
            for (k, &q) in q_dirty.iter().enumerate() {
                if k < 64 && (dirty_init >> k) & 1 != 0 {
                    *sim.qubit_mut(q) = 1;
                }
            }
            sim.apply(&ops);

            let expected: u128 = (target as u128) + (offset as u128) + (cin as u128);
            let readbits = n.min(127);
            let expected = expected & ((1u128 << readbits) - 1);
            let mut got: u128 = 0;
            for k in 0..readbits {
                if sim.qubit(q_target[k]) & 1 != 0 {
                    got |= 1u128 << k;
                }
            }
            let mut got_offset: u64 = 0;
            for k in 0..m {
                if sim.qubit(q_offset[k]) & 1 != 0 {
                    got_offset |= 1 << k;
                }
            }
            let mut got_dirty: u64 = 0;
            for (k, &q) in q_dirty.iter().enumerate() {
                if sim.qubit(q) & 1 != 0 {
                    got_dirty |= 1 << k;
                }
            }
            let dmask = if q_dirty.len() < 64 {
                (1u64 << q_dirty.len()) - 1
            } else {
                u64::MAX
            };
            let dirty_ok = got_dirty == (dirty_init & dmask & nmask);
            let offset_ok = got_offset == offset;
            let phase = sim.global_phase() & 1;

            if got == expected && dirty_ok && offset_ok && phase == 0 {
                ok += 1;
            } else {
                bad += 1;
                if bad < 3 {
                    eprintln!(
                        "narrow FAIL n={} m={} t={:#x} o={:#x} d={:#x} cin={} got={:#x} exp={:#x} d_ok={} o_ok={} phase={}",
                        n, m, target, offset, dirty_init, cin, got, expected, dirty_ok, offset_ok, phase
                    );
                }
            }
        }
        (ok, bad)
    }

    #[test]
    fn test_iadd_qoffset_narrow_small() {
        // m < n covers the shift22 use (short spill into wide register).
        for n in 5..=12 {
            for m in 1..n {
                let (ok, bad) = run_iadd_qoffset_narrow(n, m, 12);
                assert_eq!(bad, 0, "n={n} m={m}: {ok}/{} passed", ok + bad);
            }
        }
    }

    #[test]
    fn test_iadd_qoffset_narrow_wide() {
        // n=256, m=22: the exact shift22 shape.
        let (ok, bad) = run_iadd_qoffset_narrow(256, 22, 40);
        assert_eq!(bad, 0, "n=256 m=22: {ok}/{} passed", ok + bad);
    }
}
