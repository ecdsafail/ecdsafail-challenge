//! Gidney 2025 classical-quantum constant adder (arXiv:2507.23079),
//! ported via the multi-term ghost discharge API.
//!
//! Adds a classical constant `c` into an `n`-bit register `a` in place
//! using `n-1` BORROWED dirty bits (arbitrary values, restored on exit)
//! and only O(1) clean ancillae. ~3n Toffoli. The dirty bits can be any
//! already-live qubits the add does not touch (e.g. the high bits of the
//! register being added to), so the carry scratch costs ~0 peak qubits.
//!
//! Mechanism. A clean-ancilla ripple would hold all `n-1` carries at
//! once (peak +n). Instead each carry is *measurement-vented* (`hmr_ghost`)
//! as soon as the next is computed, and its value is `XORed` into a dirty
//! bit. The vented carry's deferred phase is corrected by two
//! `Z(dirty[i])` deposits — one before and one after `XORCarries`
//! restores the dirty bit:
//!
//!   dirty[i]@before = `dirty_orig`[i] XOR carry_{i+1}
//!   dirty[i]@after  = `dirty_orig`[i]
//!   term1 XOR term2 = carry_{i+1}  == the vented value  ✓
//!
//! `ghost_xor_z` accumulates each term's 64-shot sim mask; `close_ghost`
//! requires the accumulated XOR to equal the vented value's mask, so the
//! tracker verifies the cancellation on every shot before clearing the
//! obligation.

use crate::point_add::trailmix_port::circuit::{Circuit, QReg};

fn cbit(c: &[u8], i: usize) -> bool {
    let byte = i / 8;
    byte < c.len() && (c[byte] >> (i % 8)) & 1 == 1
}

/// `a += c (mod 2^n)` using `dirty` (>= n-1 borrowed bits, restored).
pub fn add_const_gidney(circ: &mut Circuit, a: &[QReg], c: &[u8], dirty: &[QReg]) {
    let n = a.len();
    if n == 0 {
        return;
    }
    if n == 1 {
        if cbit(c, 0) {
            circ.x(&a[0]);
        }
        return;
    }
    assert!(dirty.len() >= n - 1, "need n-1 borrowed dirty bits");
    let prev = circ.push_section("gidney_const");

    // ---- Forward vent pass: ripple carries, store carry_{i+1} into
    // dirty[i], form the sum in place, vent each clean carry as a ghost.
    let mut ghosts: Vec<crate::point_add::trailmix_port::tracker::ghost::Ghost> = Vec::with_capacity(n - 1);
    let mut cy = circ.alloc_qreg("gc_cy"); // carry_0 = 0
    for i in 0..(n - 1) {
        let new = circ.alloc_qreg("gc_carry");
        let anc = circ.alloc_qreg("gc_anc");
        if cbit(c, i) {
            circ.x(&anc);
        }
        circ.cx(&cy, &anc); // anc = c_i XOR carry_i
        circ.cx(&cy, &a[i]); // a[i] = a_i XOR carry_i  (= t_i)
        circ.ccx(&a[i], &anc, &new); // new = t_i AND (c_i XOR carry_i)
        circ.cx(&cy, &new); // new = MAJ = carry_{i+1}
        circ.cx(&new, &dirty[i]); // dirty[i] ^= carry_{i+1}
        circ.cx(&cy, &anc); // restore anc = c_i
        if cbit(c, i) {
            circ.x(&anc); // anc = 0
            circ.x(&a[i]); // a[i] = sum_i
        }
        circ.zero_and_free(anc);

        if i > 0 {
            ghosts.push(circ.hmr_ghost(&cy)); // vent carry_i
            circ.zero_and_free(cy);
        } else {
            circ.zero_and_free(cy); // carry_0 = 0
        }
        cy = new;
    }
    if cbit(c, n - 1) {
        circ.x(&a[n - 1]);
    }
    circ.cx(&cy, &a[n - 1]); // a[n-1] = sum
    ghosts.push(circ.hmr_ghost(&cy)); // vent carry_{n-1}
    circ.zero_and_free(cy);
    debug_assert_eq!(ghosts.len(), n - 1);
    // ghosts[i] vents carry_{i+1}; dirty[i] = dirty_orig[i] XOR carry_{i+1}.

    // ---- Correction term 1: Z(dirty[i]) (= dirty_orig XOR carry_{i+1}).
    for i in 0..(n - 1) {
        circ.ghost_xor_z(&mut ghosts[i], &dirty[i]);
    }

    // ---- Restore the dirty bits: XOR the carries back out.
    for q in a {
        circ.x(q);
    }
    xor_carries(circ, a, c, dirty);
    for q in a {
        circ.x(q);
    }

    // ---- Correction term 2: Z(dirty[i]) (= dirty_orig) + close.
    for (i, mut g) in ghosts.into_iter().enumerate() {
        circ.ghost_xor_z(&mut g, &dirty[i]);
        circ.close_ghost(g); // verifies term1 XOR term2 == carry_{i+1}
    }

    circ.pop_section(&prev);
}

/// Controlled `a += ctrl * c (mod 2^n)` using `dirty` (>= n-1 borrowed
/// bits, restored). c-loads are gated on `ctrl`; for `ctrl=0`, `a` is
/// unchanged and all vented carries are 0. Same ghost machinery.
pub fn controlled_add_const_gidney(
    circ: &mut Circuit,
    ctrl: &QReg,
    a: &[QReg],
    c: &[u8],
    dirty: &[QReg],
) {
    let ar: Vec<&QReg> = a.iter().collect();
    let dr: Vec<&QReg> = dirty.iter().collect();
    controlled_add_const_gidney_refs(circ, ctrl, &ar, c, &dr);
}

/// Reference-slice variant of [`controlled_add_const_gidney`].
pub fn controlled_add_const_gidney_refs(
    circ: &mut Circuit,
    ctrl: &QReg,
    a: &[&QReg],
    c: &[u8],
    dirty: &[&QReg],
) {
    let n = a.len();
    if n == 0 {
        return;
    }
    if n == 1 {
        if cbit(c, 0) {
            circ.cx(ctrl, a[0]);
        }
        return;
    }
    assert!(dirty.len() >= n - 1, "need n-1 borrowed dirty bits");
    let prev = circ.push_section("gidney_cadd");

    let mut ghosts: Vec<crate::point_add::trailmix_port::tracker::ghost::Ghost> = Vec::with_capacity(n - 1);
    let mut cy = circ.alloc_qreg("gcc_cy");
    for i in 0..(n - 1) {
        let new = circ.alloc_qreg("gcc_carry");
        let anc = circ.alloc_qreg("gcc_anc");
        if cbit(c, i) {
            circ.cx(ctrl, &anc); // anc = ctrl*c_i
        }
        circ.cx(&cy, &anc);
        circ.cx(&cy, a[i]);
        circ.ccx(a[i], &anc, &new);
        circ.cx(&cy, &new); // new = carry_{i+1}
        circ.cx(&new, dirty[i]);
        circ.cx(&cy, &anc);
        if cbit(c, i) {
            circ.cx(ctrl, &anc); // anc = 0
            circ.cx(ctrl, a[i]); // a[i] = sum_i
        }
        circ.zero_and_free(anc);

        if i > 0 {
            ghosts.push(circ.hmr_ghost(&cy));
            circ.zero_and_free(cy);
        } else {
            circ.zero_and_free(cy);
        }
        cy = new;
    }
    if cbit(c, n - 1) {
        circ.cx(ctrl, a[n - 1]);
    }
    circ.cx(&cy, a[n - 1]);
    ghosts.push(circ.hmr_ghost(&cy));
    circ.zero_and_free(cy);

    for i in 0..(n - 1) {
        circ.ghost_xor_z(&mut ghosts[i], dirty[i]);
    }
    for q in a {
        circ.x(q);
    }
    xor_carries_ctrl_refs(circ, ctrl, a, c, dirty);
    for q in a {
        circ.x(q);
    }
    for (i, mut g) in ghosts.into_iter().enumerate() {
        circ.ghost_xor_z(&mut g, dirty[i]);
        circ.close_ghost(g);
    }
    circ.pop_section(&prev);
}

/// Controlled variant of [`xor_carries`]: condition flips gated on `ctrl`.
fn xor_carries_ctrl_refs(circ: &mut Circuit, ctrl: &QReg, a: &[&QReg], c: &[u8], out: &[&QReg]) {
    let n = a.len();
    let ccx_cond = |circ: &mut Circuit, c1: &QReg, c2: &QReg, t: &QReg, b0: bool, b1: bool| {
        if b0 {
            circ.cx(ctrl, c1);
        }
        if b1 {
            circ.cx(ctrl, c2);
        }
        circ.ccx(c1, c2, t);
        if b0 {
            circ.cx(ctrl, c1);
        }
        if b1 {
            circ.cx(ctrl, c2);
        }
    };
    for i in (1..(n - 1)).rev() {
        ccx_cond(circ, a[i], out[i - 1], out[i], cbit(c, i), false);
    }
    for i in 0..(n - 1) {
        if cbit(c, i) {
            circ.cx(ctrl, out[i]);
        }
    }
    let cin = circ.alloc_qreg("xcc_cin");
    ccx_cond(circ, &cin, a[0], out[0], cbit(c, 0), cbit(c, 0));
    circ.zero_and_free(cin);
    for i in 1..(n - 1) {
        ccx_cond(circ, a[i], out[i - 1], out[i], cbit(c, i), cbit(c, i));
    }
}

/// Involutory `XORCarries`: recompute the `n-1` carries of `a + c` (with
/// `a` the complemented sum, per the caller) and XOR them into `out`
/// (= dirty). Composed with the forward `dirty ^= carry`, restores `out`.
fn xor_carries(circ: &mut Circuit, a: &[QReg], c: &[u8], out: &[QReg]) {
    let n = a.len();
    let ccx_cond = |circ: &mut Circuit, c1: &QReg, c2: &QReg, t: &QReg, b0: bool, b1: bool| {
        if b0 {
            circ.x(c1);
        }
        if b1 {
            circ.x(c2);
        }
        circ.ccx(c1, c2, t);
        if b0 {
            circ.x(c1);
        }
        if b1 {
            circ.x(c2);
        }
    };
    for i in (1..(n - 1)).rev() {
        ccx_cond(circ, &a[i], &out[i - 1], &out[i], cbit(c, i), false);
    }
    for i in 0..(n - 1) {
        if cbit(c, i) {
            circ.x(&out[i]);
        }
    }
    let cin = circ.alloc_qreg("xc_cin");
    ccx_cond(circ, &cin, &a[0], &out[0], cbit(c, 0), cbit(c, 0));
    circ.zero_and_free(cin);
    for i in 1..(n - 1) {
        ccx_cond(circ, &a[i], &out[i - 1], &out[i], cbit(c, i), cbit(c, i));
    }
}

/// Restore variant of [`xor_carries`] for the comparator: `out` holds ALL
/// `n` carries `carry_1..carry_n` of `a + c` (XOR'd into `out[0..n]`), and
/// this XORs them back out (re-deriving from `a`). Mirror of `xor_carries`
/// with the carry index extended to `n` (the overflow `carry_n` included).
fn xor_carries_all_refs(circ: &mut Circuit, a: &[&QReg], c: &[u8], out: &[&QReg]) {
    let n = a.len();
    let ccx_cond = |circ: &mut Circuit, c1: &QReg, c2: &QReg, t: &QReg, b0: bool, b1: bool| {
        if b0 {
            circ.x(c1);
        }
        if b1 {
            circ.x(c2);
        }
        circ.ccx(c1, c2, t);
        if b0 {
            circ.x(c1);
        }
        if b1 {
            circ.x(c2);
        }
    };
    for i in (1..n).rev() {
        ccx_cond(circ, a[i], out[i - 1], out[i], cbit(c, i), false);
    }
    for i in 0..n {
        if cbit(c, i) {
            circ.x(out[i]);
        }
    }
    let cin = circ.alloc_qreg("xca_cin");
    ccx_cond(circ, &cin, a[0], out[0], cbit(c, 0), cbit(c, 0));
    circ.zero_and_free(cin);
    for i in 1..n {
        ccx_cond(circ, a[i], out[i - 1], out[i], cbit(c, i), cbit(c, i));
    }
}

/// `out ^= (a >= k)` for a classical constant `k` (n = `a.len()` bits), `a`
/// preserved. The cheap (~3n Toffoli) low-clean-peak Gidney constant
/// comparator: ripple the carry of `a + (2^n - k)` (whose overflow `carry_n`
/// is `1 iff a >= k`) using `n` BORROWED dirty bits (restored) and O(1)
/// clean ancilla, venting each carry via X-basis measurement (Clifford, no
/// Toffoli). Unlike the adder it does NOT form the sum (a is restored each
/// column) and it grabs the overflow carry into `out`. Replaces the
/// O(n log n) `compare_geq_theorem3` for hot-path constant compares.
pub fn compare_geq_const_gidney(
    circ: &mut Circuit,
    a: &[QReg],
    k: &[u8],
    out: &QReg,
    dirty: &[QReg],
) {
    let ar: Vec<&QReg> = a.iter().collect();
    let dr: Vec<&QReg> = dirty.iter().collect();
    compare_geq_const_gidney_refs(circ, &ar, k, out, &dr);
}

/// Reference-slice variant of [`compare_geq_const_gidney`].
pub fn compare_geq_const_gidney_refs(
    circ: &mut Circuit,
    a: &[&QReg],
    k: &[u8],
    out: &QReg,
    dirty: &[&QReg],
) {
    let n = a.len();
    // k as integer (n fits in u128 for our small comparator widths).
    let kv: u128 = (0..n.min(128))
        .filter(|&i| cbit(k, i))
        .map(|i| 1u128 << i)
        .sum();
    if n == 0 {
        if kv == 0 {
            circ.x(out);
        }
        return;
    }
    if kv == 0 {
        circ.x(out); // a >= 0 always
        return;
    }
    if kv >= (1u128 << n) {
        return; // a < k always; out unchanged
    }
    // c = 2^n - k  (n-bit constant, in [1, 2^n))
    let cv: u128 = (1u128 << n) - kv;
    let c: Vec<u8> = (0..n.div_ceil(8)).map(|b| (cv >> (8 * b)) as u8).collect();
    assert!(
        dirty.len() >= n,
        "compare_geq_const_gidney needs n borrowed dirty bits"
    );

    let prev = circ.push_section("cmp_geq_gidney");
    let mut ghosts: Vec<crate::point_add::trailmix_port::tracker::ghost::Ghost> = Vec::with_capacity(n);
    let mut cy = circ.alloc_qreg("cmpg_cy"); // carry_0 = 0
    for i in 0..n {
        let new = circ.alloc_qreg("cmpg_carry");
        let anc = circ.alloc_qreg("cmpg_anc");
        if cbit(&c, i) {
            circ.x(&anc);
        }
        circ.cx(&cy, &anc); // anc = c_i XOR carry_i
        circ.cx(&cy, a[i]); // a[i] = t_i = a_i XOR carry_i
        circ.ccx(a[i], &anc, &new); // new = t_i AND (c_i XOR carry_i)
        circ.cx(&cy, &new); // new = MAJ = carry_{i+1}
        circ.cx(&new, dirty[i]); // dirty[i] ^= carry_{i+1}
        if i == n - 1 {
            circ.cx(&new, out); // grab carry_n = (a >= k)
        }
        circ.cx(&cy, a[i]); // RESTORE a[i] = a_i (vs adder: forms sum)
        circ.cx(&cy, &anc); // restore anc = c_i
        if cbit(&c, i) {
            circ.x(&anc);
        }
        circ.zero_and_free(anc);
        if i > 0 {
            ghosts.push(circ.hmr_ghost(&cy)); // vent carry_i
        }
        circ.zero_and_free(cy);
        cy = new;
    }
    ghosts.push(circ.hmr_ghost(&cy)); // vent carry_n
    circ.zero_and_free(cy);
    debug_assert_eq!(ghosts.len(), n);
    // ghosts[i] vents carry_{i+1}; dirty[i] = dirty_orig[i] XOR carry_{i+1}.

    for i in 0..n {
        circ.ghost_xor_z(&mut ghosts[i], dirty[i]); // term1 = dirty_orig ^ carry_{i+1}
    }
    xor_carries_all_refs(circ, a, &c, dirty); // restore dirty -> dirty_orig
    for (i, mut g) in ghosts.into_iter().enumerate() {
        circ.ghost_xor_z(&mut g, dirty[i]); // term2 = dirty_orig
        circ.close_ghost(g); // verifies term1 ^ term2 == carry_{i+1}
    }
    circ.pop_section(&prev);
}

/// Controlled hybrid TTK-Gidney register adder (Schrottenloher's
/// `ControlledHybridAdder`, itself Gidney 2018 arXiv:1709.06648 Fig.4a fused
/// with the TTK in-place carry trick). Computes `a += ctrl * b (mod 2^n)`;
/// `b` and `ctrl` are preserved.
///
/// Carries are threaded in place through `b` (TTK), so the only extra qubits
/// are the `vents` measurement-vent ancillae. Each vent ancilla replaces one
/// carry-*uncompute* Toffoli with a measurement (Gidney's measure-and-fixup
/// AND erasure, Fig.3 bottom): the AND `a[i] & b[i]` is computed into the
/// ancilla (1 Toffoli), then erased by an X-basis measurement (`hmr_ghost`)
/// whose phase kickback is cancelled by `CZ(a[i], b[i])` gated on the measured
/// bit (`ghost_xor_cz`). The AND inputs `a[i]`, `b[i]` are untouched between
/// compute (forward step i) and erase (reverse step i), so the CZ targets are
/// alive and `close_ghost` sim-verifies the cancellation on all 64 shots.
///
/// Total Toffoli = `3n - 2 - vents`, where `vents = min(vents_budget, n-1)`:
///   * `vents = 0`     -> TTK/Cuccaro `3n` controlled adder,
///   * `vents = n - 1` -> Gidney `2n` controlled adder.
/// Peak qubits: `+vents` (held between the forward and reverse carry chains).
pub fn controlled_hybrid_add(
    circ: &mut Circuit,
    ctrl: &QReg,
    a: &[QReg],
    b: &[QReg],
    vents_budget: usize,
) {
    let aref: Vec<&QReg> = a.iter().collect();
    let bref: Vec<&QReg> = b.iter().collect();
    controlled_hybrid_add_refs(circ, ctrl, &aref, &bref, vents_budget);
}

/// Refs variant of [`controlled_hybrid_add`] for non-contiguous operand windows
/// (e.g. the shifted/scattered registers in the shrunken-PZ divstep). Identical
/// gate sequence; the vents are freshly-allocated measurement ancillae, so there
/// is no contiguity/borrowed-dirty assumption on `a`/`b`.
pub fn controlled_hybrid_add_refs(
    circ: &mut Circuit,
    ctrl: &QReg,
    a: &[&QReg],
    b: &[&QReg],
    vents_budget: usize,
) {
    let n = a.len();
    assert_eq!(b.len(), n, "controlled_hybrid_add: a, b must match width");
    if n == 0 {
        return;
    }
    if n == 1 {
        circ.ccx(ctrl, b[0], a[0]);
        return;
    }
    let vents = vents_budget.min(n - 1);
    let prev = circ.push_section("hybrid_cadd");

    // qarton qr_x = b (addend, carry-threaded), qr_y = a (target).
    for i in 1..n {
        circ.cx(b[i], a[i]);
    }
    for i in (1..n - 1).rev() {
        circ.cx(b[i], b[i + 1]);
    }

    // Forward carry chain. The first `vents` carries land in measurement-vent
    // ancillae; the rest are Toffoli'd straight into b.
    let mut vent_ancs: Vec<Option<QReg>> = (0..n - 1).map(|_| None).collect();
    for i in 0..n - 1 {
        if i < vents {
            let anc = circ.alloc_qreg("hyb_vent");
            circ.ccx(a[i], b[i], &anc); // anc = a[i] & b[i]
            circ.cx(&anc, b[i + 1]);
            vent_ancs[i] = Some(anc);
        } else {
            circ.ccx(a[i], b[i], b[i + 1]);
        }
    }

    // Reverse: write the controlled sum bit, then uncompute each carry.
    for i in (0..n - 1).rev() {
        circ.ccx(ctrl, b[i + 1], a[i + 1]); // controlled sum bit i+1
        if i < vents {
            let anc = vent_ancs[i].take().unwrap();
            circ.cx(&anc, b[i + 1]); // undo the forward cx
                                     // Measure-and-fixup AND erasure: HMR(anc), then CZ(a[i], b[i]).
            let mut g = circ.hmr_ghost(&anc);
            circ.zero_and_free(anc);
            circ.ghost_xor_cz(&mut g, a[i], b[i]);
            circ.close_ghost(g);
        } else {
            circ.ccx(a[i], b[i], b[i + 1]);
        }
    }

    for i in 1..n - 1 {
        circ.cx(b[i], b[i + 1]);
    }
    circ.ccx(ctrl, b[0], a[0]);
    for i in 1..n {
        circ.cx(b[i], a[i]);
    }
    circ.pop_section(&prev);
}

/// UNCONDITIONAL measurement-vented adder `a += b mod 2^n` (b restored). Same
/// vented carry chain as [`controlled_hybrid_add_refs`] but the sum bits are
/// plain `cx` (no control) -- so it costs ONLY the carry chain: ~n Toffoli at
/// `vents = n-1` (vs Cuccaro's ~2n), using `vents` clean measurement ancillae.
/// Carry-out beyond `a.len()` is dropped (mod 2^n). Internally uses HMR for the
/// vent erasure (self-contained); the call is NOT gate-reversible via
/// `emit_reverse_since` -- hand-reverse if you need its inverse.
pub fn hybrid_add_refs(circ: &mut Circuit, a: &[&QReg], b: &[&QReg], vents_budget: usize) {
    let n = a.len();
    assert_eq!(b.len(), n, "hybrid_add: a, b must match width");
    if n == 0 {
        return;
    }
    if n == 1 {
        circ.cx(b[0], a[0]);
        return;
    }
    let vents = vents_budget.min(n - 1);
    let prev = circ.push_section("hybrid_add");
    for i in 1..n {
        circ.cx(b[i], a[i]);
    }
    for i in (1..n - 1).rev() {
        circ.cx(b[i], b[i + 1]);
    }
    let mut vent_ancs: Vec<Option<QReg>> = (0..n - 1).map(|_| None).collect();
    for i in 0..n - 1 {
        if i < vents {
            let anc = circ.alloc_qreg("hyb_vent");
            circ.ccx(a[i], b[i], &anc);
            circ.cx(&anc, b[i + 1]);
            vent_ancs[i] = Some(anc);
        } else {
            circ.ccx(a[i], b[i], b[i + 1]);
        }
    }
    for i in (0..n - 1).rev() {
        circ.cx(b[i + 1], a[i + 1]); // UNCONDITIONAL sum bit i+1
        if i < vents {
            let anc = vent_ancs[i].take().unwrap();
            circ.cx(&anc, b[i + 1]);
            let mut g = circ.hmr_ghost(&anc);
            circ.zero_and_free(anc);
            circ.ghost_xor_cz(&mut g, a[i], b[i]);
            circ.close_ghost(g);
        } else {
            circ.ccx(a[i], b[i], b[i + 1]);
        }
    }
    for i in 1..n - 1 {
        circ.cx(b[i], b[i + 1]);
    }
    circ.cx(b[0], a[0]); // UNCONDITIONAL sum bit 0
    for i in 1..n {
        circ.cx(b[i], a[i]);
    }
    circ.pop_section(&prev);
}

/// Slice wrapper for [`hybrid_add_refs`].
pub fn hybrid_add(circ: &mut Circuit, a: &[QReg], b: &[QReg], vents_budget: usize) {
    let aref: Vec<&QReg> = a.iter().collect();
    let bref: Vec<&QReg> = b.iter().collect();
    hybrid_add_refs(circ, &aref, &bref, vents_budget);
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::{thread_rng, Rng};

    #[test]
    fn gidney_hybrid_cadd_value_and_phase_clean_n8() {
        let n = 8usize;
        let mut rng = thread_rng();
        // Exercise pure-TTK (0), mixed, and full-Gidney (n-1) vent budgets.
        for &vents in &[0usize, 1, 4, 7, 100] {
            let mut circ = Circuit::new();
            let ctrl = circ.alloc_qreg("ctrl");
            let a = circ.alloc_qreg_bits("a", n);
            let b = circ.alloc_qreg_bits("b", n);

            let mut a_in = [0u64; 64];
            let mut b_in = [0u64; 64];
            let mut ctrl_in = [0u8; 64];
            for shot in 0..64 {
                let cv: u8 = rng.gen::<u8>() & 1;
                ctrl_in[shot] = cv;
                if cv == 1 {
                    circ.sim_load_reg_bytes_shot(std::slice::from_ref(&ctrl), &[1u8], shot);
                }
                let av: u64 = rng.gen::<u64>() & ((1 << n) - 1);
                let bv: u64 = rng.gen::<u64>() & ((1 << n) - 1);
                a_in[shot] = av;
                b_in[shot] = bv;
                circ.sim_load_reg_bytes_shot(&a, &av.to_le_bytes(), shot);
                circ.sim_load_reg_bytes_shot(&b, &bv.to_le_bytes(), shot);
            }

            controlled_hybrid_add(&mut circ, &ctrl, &a, &b, vents);
            circ.assert_phase_clean();

            let mut outs: Vec<QReg> = vec![ctrl];
            outs.extend(a);
            outs.extend(b);
            let (sim, det) = circ.destroy_sim(outs);
            for shot in 0..64 {
                let mut got_a: u64 = 0;
                for i in 0..n {
                    if sim.read_bit_shot(&det[1 + i], shot) == 1 {
                        got_a |= 1 << i;
                    }
                }
                let want = if ctrl_in[shot] == 1 {
                    (a_in[shot] + b_in[shot]) & ((1 << n) - 1)
                } else {
                    a_in[shot]
                };
                assert_eq!(got_a, want, "vents={vents} shot {shot}: a+=ctrl*b");
                let mut got_b: u64 = 0;
                for i in 0..n {
                    if sim.read_bit_shot(&det[1 + n + i], shot) == 1 {
                        got_b |= 1 << i;
                    }
                }
                assert_eq!(
                    got_b, b_in[shot],
                    "vents={vents} shot {shot}: b not restored"
                );
            }
        }
    }

    #[test]
    fn gidney_const_value_and_phase_clean_n8() {
        let n = 8usize;
        let mut rng = thread_rng();
        let c: u64 = 0b10110101;
        let c_bytes = c.to_le_bytes();

        let mut circ = Circuit::new();
        let a = circ.alloc_qreg_bits("a", n);
        let dirty = circ.alloc_qreg_bits("dirty", n - 1);

        let mut a_in = [0u64; 64];
        let mut d_in = [0u64; 64];
        for shot in 0..64 {
            let av: u64 = rng.gen::<u64>() & ((1 << n) - 1);
            let dv: u64 = rng.gen::<u64>() & ((1 << (n - 1)) - 1);
            a_in[shot] = av;
            d_in[shot] = dv;
            circ.sim_load_reg_bytes_shot(&a, &av.to_le_bytes(), shot);
            circ.sim_load_reg_bytes_shot(&dirty, &dv.to_le_bytes(), shot);
        }

        add_const_gidney(&mut circ, &a, &c_bytes, &dirty);
        circ.assert_phase_clean();

        let mut outs: Vec<QReg> = Vec::new();
        outs.extend(a);
        outs.extend(dirty);
        let (sim, det) = circ.destroy_sim(outs);
        for shot in 0..64 {
            let mut got_a: u64 = 0;
            for i in 0..n {
                if sim.read_bit_shot(&det[i], shot) == 1 {
                    got_a |= 1 << i;
                }
            }
            assert_eq!(got_a, (a_in[shot] + c) & ((1 << n) - 1), "shot {shot}: a+c");
            let mut got_d: u64 = 0;
            for i in 0..(n - 1) {
                if sim.read_bit_shot(&det[n + i], shot) == 1 {
                    got_d |= 1 << i;
                }
            }
            assert_eq!(got_d, d_in[shot], "shot {shot}: dirty not restored");
        }
    }

    #[test]
    fn compare_geq_const_gidney_value_and_phase_clean_n8() {
        let n = 8usize;
        let mut rng = thread_rng();
        for &k in &[1u64, 3, 47, 81, 128, 162, 200, 255] {
            let k_bytes = k.to_le_bytes();
            let mut circ = Circuit::new();
            let a = circ.alloc_qreg_bits("a", n);
            let dirty = circ.alloc_qreg_bits("dirty", n); // n borrowed dirty
            let out = circ.alloc_qreg("out");
            let mut a_in = [0u64; 64];
            let mut d_in = [0u64; 64];
            for shot in 0..64 {
                let av: u64 = rng.gen::<u64>() & ((1 << n) - 1);
                let dv: u64 = rng.gen::<u64>() & ((1 << n) - 1);
                a_in[shot] = av;
                d_in[shot] = dv;
                circ.sim_load_reg_bytes_shot(&a, &av.to_le_bytes(), shot);
                circ.sim_load_reg_bytes_shot(&dirty, &dv.to_le_bytes(), shot);
            }
            let ccx0 = circ.ccx_emitted;
            let ccz0 = circ.ccz_emitted;
            compare_geq_const_gidney(&mut circ, &a, &k_bytes, &out, &dirty);
            if k == 81 {
                eprintln!(
                    "  compare_geq_const_gidney(n=8) tof={}",
                    (circ.ccx_emitted - ccx0) + (circ.ccz_emitted - ccz0)
                );
            }
            circ.assert_phase_clean();
            let mut outs: Vec<QReg> = Vec::new();
            outs.extend(a);
            outs.extend(dirty);
            outs.push(out);
            let (sim, det) = circ.destroy_sim(outs);
            for shot in 0..64 {
                let mut got_a: u64 = 0;
                for i in 0..n {
                    if sim.read_bit_shot(&det[i], shot) == 1 {
                        got_a |= 1 << i;
                    }
                }
                assert_eq!(got_a, a_in[shot], "k={k} shot {shot}: a not preserved");
                let mut got_d: u64 = 0;
                for i in 0..n {
                    if sim.read_bit_shot(&det[n + i], shot) == 1 {
                        got_d |= 1 << i;
                    }
                }
                assert_eq!(got_d, d_in[shot], "k={k} shot {shot}: dirty not restored");
                let got_out = sim.read_bit_shot(&det[2 * n], shot);
                let want = u8::from(a_in[shot] >= k);
                assert_eq!(
                    got_out, want,
                    "k={k} shot {shot}: a={} >= {k}? want {want} got {got_out}",
                    a_in[shot]
                );
            }
        }
    }

    #[test]
    fn gidney_controlled_const_value_and_phase_clean_n8() {
        let n = 8usize;
        let mut rng = thread_rng();
        let c: u64 = 0b01101011;
        let c_bytes = c.to_le_bytes();

        let mut circ = Circuit::new();
        let ctrl = circ.alloc_qreg("ctrl");
        let a = circ.alloc_qreg_bits("a", n);
        let dirty = circ.alloc_qreg_bits("dirty", n - 1);

        let mut a_in = [0u64; 64];
        let mut d_in = [0u64; 64];
        let mut ctrl_in = [0u8; 64];
        for shot in 0..64 {
            let cv: u8 = rng.gen::<u8>() & 1;
            ctrl_in[shot] = cv;
            if cv == 1 {
                circ.sim_load_reg_bytes_shot(std::slice::from_ref(&ctrl), &[1u8], shot);
            }
            let av: u64 = rng.gen::<u64>() & ((1 << n) - 1);
            let dv: u64 = rng.gen::<u64>() & ((1 << (n - 1)) - 1);
            a_in[shot] = av;
            d_in[shot] = dv;
            circ.sim_load_reg_bytes_shot(&a, &av.to_le_bytes(), shot);
            circ.sim_load_reg_bytes_shot(&dirty, &dv.to_le_bytes(), shot);
        }

        controlled_add_const_gidney(&mut circ, &ctrl, &a, &c_bytes, &dirty);
        circ.assert_phase_clean();

        let mut outs: Vec<QReg> = vec![ctrl];
        outs.extend(a);
        outs.extend(dirty);
        let (sim, det) = circ.destroy_sim(outs);
        for shot in 0..64 {
            let mut got_a: u64 = 0;
            for i in 0..n {
                if sim.read_bit_shot(&det[1 + i], shot) == 1 {
                    got_a |= 1 << i;
                }
            }
            let want = if ctrl_in[shot] == 1 {
                (a_in[shot] + c) & ((1 << n) - 1)
            } else {
                a_in[shot]
            };
            assert_eq!(got_a, want, "shot {shot}: ctrl={} a+c", ctrl_in[shot]);
            let mut got_d: u64 = 0;
            for i in 0..(n - 1) {
                if sim.read_bit_shot(&det[1 + n + i], shot) == 1 {
                    got_d |= 1 << i;
                }
            }
            assert_eq!(got_d, d_in[shot], "shot {shot}: dirty not restored");
        }
    }
}
