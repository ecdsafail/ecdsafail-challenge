//! Controlled qubit-qubit (q-q) vented-adder family (the GCD-subtract + apply cofactor-add
//! path) for the reversible secp256k1 EC point-add circuit, built on this
//! crate's `B` builder. Carries are vented two ways:
//!   * AND-carry erase: `let b = alloc_bit(); hmr(q,b); cz_if_bit(..,b)`
//!     -- an X-basis HMR measurement followed by a conditional-Z phase correction.
//!   * conditional erase: `let b = alloc_bit(); hmr(carry,b); push_condition(b);
//!     middle(..|..|{z/cz/ccz/neg}..); pop_condition()` -- the
//!     boundary carry is measured, then the comparator middle applies the gated
//!     phase correction.
//! `k` is the available qubit headroom, passed in from the schedule.

use super::comparator::compare_geq_cin_middle;
use super::{B, BExt};
use crate::circuit::{QubitId};

// ============================================================================
// controlled_hybrid_add_refs  (plain Gidney AND-carry adder)
// ============================================================================
pub fn controlled_hybrid_add_refs(circ: &mut B, ctrl: &QubitId, a: &[&QubitId], b: &[&QubitId]) {
    controlled_hybrid_add_refs_impl(circ, ctrl, a, b, false);
}

/// As [`controlled_hybrid_add_refs`] but skips the bit-0 controlled sum
/// `a[0] ^= ctrl AND b[0]`. Used by [`controlled_vented_chunk_add`], whose `|1>`
/// low pad's bit-0 sum is immediately undone by the caller's `ccx(ctrl,cin,one)`
/// restore -- the two cancel. Skipping the sum here and dropping that restore
/// reproduces the elided circuit directly (a[0] is a discarded pad, so its sum
/// value is irrelevant; the bit-0 carry into the chunk is still formed). Saves
/// 2 CCX per chunk.
fn controlled_hybrid_add_refs_skiplow(circ: &mut B, ctrl: &QubitId, a: &[&QubitId], b: &[&QubitId]) {
    controlled_hybrid_add_refs_impl(circ, ctrl, a, b, true);
}

fn controlled_hybrid_add_refs_impl(circ: &mut B, ctrl: &QubitId, a: &[&QubitId], b: &[&QubitId], skip_low_ctrl_sum: bool) {
    let n = a.len();
    assert_eq!(b.len(), n, "controlled_hybrid_add: a, b must match width");
    if n == 0 {
        return;
    }
    if n == 1 {
        circ.ccx(*ctrl, *b[0], *a[0]);
        return;
    }
    // Effective vent count, read from the baked schedule; no decision logic here.
    // When no schedule is loaded (standalone primitive unit tests) `next_hyb_v`
    // returns usize::MAX, which the `i < vents` loop reads as "vent every carry"
    // -- still value-correct.
    let vents = super::next_hyb_v();

    // b is the addend (carry-threaded); a is the target.
    for i in 1..n {
        circ.cx(*b[i], *a[i]);
    }
    for i in (1..n - 1).rev() {
        circ.cx(*b[i], *b[i + 1]);
    }

    // Forward carry chain. First `vents` carries land in measurement-vent ancillae;
    // the rest are Toffoli'd straight into b.
    let mut vent_ancs: Vec<Option<QubitId>> = (0..n - 1).map(|_| None).collect();
    for i in 0..n - 1 {
        if i < vents {
            let anc = circ.alloc_qubit();
            circ.ccx(*a[i], *b[i], anc); // anc = a[i] & b[i]
            circ.cx(anc, *b[i + 1]);
            vent_ancs[i] = Some(anc);
        } else {
            circ.ccx(*a[i], *b[i], *b[i + 1]);
        }
    }

    // Reverse: write the controlled sum bit, then uncompute each carry. Vented
    // lanes measure-erase; the rest re-Toffoli.
    for i in (0..n - 1).rev() {
        circ.ccx(*ctrl, *b[i + 1], *a[i + 1]); // controlled sum bit i+1
        if let Some(anc) = vent_ancs[i].take() {
            circ.cx(anc, *b[i + 1]); // undo the forward cx
            // AND-carry erase.
            let bit = circ.alloc_bit();
            circ.hmr(anc, bit);
            circ.zero_and_free(anc);
            circ.cz_if_bit(*a[i], *b[i], bit);
        } else {
            circ.ccx(*a[i], *b[i], *b[i + 1]);
        }
    }

    for i in 1..n - 1 {
        circ.cx(*b[i], *b[i + 1]);
    }
    if !skip_low_ctrl_sum {
        circ.ccx(*ctrl, *b[0], *a[0]);
    }
    for i in 1..n {
        circ.cx(*b[i], *a[i]);
    }
}

// ============================================================================
// controlled_clean_add_threaded  (size-s chunk add, threaded cin/cout, `vents`)
// ============================================================================
fn controlled_clean_add_threaded(
    circ: &mut B,
    ctrl: &QubitId,
    a: &[&QubitId],
    b: &[&QubitId],
    cin: Option<&QubitId>,
    cout: Option<&QubitId>,
    vents: usize,
) {
    let s = a.len();
    if s == 0 {
        if let (Some(ci), Some(co)) = (cin, cout) {
            circ.ccx(*ctrl, *ci, *co);
        }
        return;
    }
    let n_inner = if cout.is_some() { s } else { s - 1 };
    let mut inner: Vec<Option<QubitId>> = (0..n_inner).map(|_| Some(circ.alloc_qubit())).collect();
    let produces = |i: usize| cout.is_some() || i + 1 < s;
    // Forward MAJ (unconditional -- not gated by ctrl).
    for i in 0..s {
        if !produces(i) {
            continue;
        }
        let co = inner[i].as_ref().unwrap();
        let ci: Option<&QubitId> = if i == 0 { cin } else { inner[i - 1].as_ref() };
        if let Some(ci) = ci {
            circ.cx(*ci, *a[i]);
            circ.cx(*ci, *b[i]);
            circ.ccx(*a[i], *b[i], *co);
            circ.cx(*ci, *co);
        } else {
            circ.ccx(*a[i], *b[i], *co);
        }
    }
    // Extract the gated boundary carry: cout = ctrl AND top_carry (c_s).
    if let Some(cout) = cout {
        circ.ccx(*ctrl, *inner[s - 1].as_ref().unwrap(), *cout);
    }
    // Reverse: gated sums, clear every internal carry.
    for i in (0..s).rev() {
        if !produces(i) {
            let ci: Option<&QubitId> = if i == 0 { cin } else { inner[i - 1].as_ref() };
            if let Some(ci) = ci {
                circ.cx(*ci, *b[i]);
            }
            circ.ccx(*ctrl, *b[i], *a[i]);
            if let Some(ci) = ci {
                circ.cx(*ci, *b[i]);
            }
            continue;
        }
        let co = inner[i].take().unwrap();
        let ci: Option<&QubitId> = if i == 0 { cin } else { inner[i - 1].as_ref() };
        if let Some(ci) = ci {
            circ.cx(*ci, co); // co = a[i] & b[i] (folded)
        }
        if i < vents {
            let bit = circ.alloc_bit();
            circ.hmr(co, bit);
            circ.zero_and_free(co);
            circ.cz_if_bit(*a[i], *b[i], bit);
        } else {
            circ.ccx(*a[i], *b[i], co);
            circ.zero_and_free(co);
        }
        if let Some(ci) = ci {
            circ.cx(*ci, *a[i]); // a = a_i (remove the forward fold)
        }
        circ.ccx(*ctrl, *b[i], *a[i]); // a ^= ctrl & (b_i ^ ci) -> sum
        if let Some(ci) = ci {
            circ.cx(*ci, *b[i]); // restore b
        }
    }
}

// ============================================================================
// controlled_erase_carry_gated[_capped]  (the boundary-carry gated vented erase)
// ============================================================================
fn deref(s: &[&QubitId]) -> Vec<QubitId> {
    s.iter().map(|q| **q).collect()
}

fn controlled_erase_carry_gated(circ: &mut B, ctrl: &QubitId, a: &[&QubitId], b: &[&QubitId], cin: &QubitId, carry: QubitId) {
    let bit = circ.alloc_bit();
    circ.hmr(carry, bit);
    circ.push_condition(bit);
    let (av, bv) = (deref(a), deref(b));
    let ctrl = *ctrl;
    compare_geq_cin_middle(circ, &av, &bv, cin, |c, ta, tb, c_prev| {
        // ctrl . NOT c_n = ctrl.1 ^ ctrl.(ta&tb) ^ ctrl.c_prev:
        c.z(ctrl);
        c.ccz(ctrl, *ta, *tb);
        c.cz(ctrl, *c_prev);
    });
    circ.pop_condition();
    circ.zero_and_free(carry);
}

fn controlled_erase_carry_gated_capped(
    circ: &mut B,
    ctrl: &QubitId,
    a: &[&QubitId],
    b: &[&QubitId],
    cin: &QubitId,
    carry: QubitId,
    cap: usize,
) {
    let s = a.len();
    if s <= cap {
        controlled_erase_carry_gated(circ, ctrl, a, b, cin, carry);
        return;
    }
    let lo = s - cap;
    let bit = circ.alloc_bit();
    circ.hmr(carry, bit);
    circ.push_condition(bit);
    let zcin = circ.alloc_qubit();
    let (av, bv) = (deref(&a[lo..]), deref(&b[lo..]));
    let ctrl = *ctrl;
    compare_geq_cin_middle(circ, &av, &bv, &zcin, |c, ta, tb, c_prev| {
        c.z(ctrl);
        c.ccz(ctrl, *ta, *tb);
        c.cz(ctrl, *c_prev);
    });
    circ.zero_and_free(zcin);
    circ.pop_condition();
    circ.zero_and_free(carry);
}

// ============================================================================
// controlled_vented_chunk_add  (padded plain-Gidney add with cin/cout)
// ============================================================================
fn controlled_vented_chunk_add(circ: &mut B, ctrl: &QubitId, a_chunk: &[&QubitId], b_chunk: &[&QubitId], cin: &QubitId, cout: &QubitId) {
    let one = circ.alloc_qubit();
    circ.x(one);
    let zero = circ.alloc_qubit();
    let mut aext: Vec<&QubitId> = Vec::with_capacity(a_chunk.len() + 2);
    aext.push(&one);
    aext.extend_from_slice(a_chunk);
    aext.push(cout);
    let mut bext: Vec<&QubitId> = Vec::with_capacity(b_chunk.len() + 2);
    bext.push(cin);
    bext.extend_from_slice(b_chunk);
    bext.push(&zero);
    // Skip the adder's bit-0 sum `one ^= ctrl AND cin` and its restore (they
    // cancel). `one` is left at |1> by the adder, so x() returns it to |0>
    // directly. Saves 2 CCX/chunk.
    controlled_hybrid_add_refs_skiplow(circ, ctrl, &aext, &bext);
    circ.x(one);
    circ.zero_and_free(one);
    circ.zero_and_free(zero);
}

// ============================================================================
// varchunk_schedule / varchunk_cost
// ============================================================================
fn varchunk_schedule(n: usize, k: usize) -> Vec<usize> {
    const RESERVE: usize = 4;
    let mut sizes = Vec::new();
    let (mut covered, mut held) = (0usize, 0usize);
    while covered < n {
        let room = k.saturating_sub(held + RESERVE);
        if room == 0 {
            return Vec::new();
        }
        let s = room.min(n - covered);
        sizes.push(s);
        covered += s;
        held += 1;
    }
    sizes
}

fn varchunk_cost(n: usize, k: usize, cap: usize) -> usize {
    let sizes = varchunk_schedule(n, k);
    if sizes.is_empty() {
        return usize::MAX;
    }
    let erase: usize = sizes.iter().map(|&s| s.min(cap) / 2).sum();
    n + erase
}

// ============================================================================
// adaptive_layout / adaptive_add_cost_tof
// ============================================================================
pub(crate) struct AdaptiveLayout {
    pub(crate) c: usize,
    pub(crate) chunked_len: usize,
    pub(crate) plain_len: usize,
}
pub(crate) const ADAPTIVE_RES: usize = 5;

pub(crate) fn adaptive_layout(n: usize, k: usize) -> AdaptiveLayout {
    let c = ((n as f64).sqrt() as usize).clamp(1, n);
    let mut plain = 0usize;
    while plain < n {
        let l = n - (plain + 1);
        let nch = l.div_ceil(c);
        if nch + (plain + 1) <= k {
            plain += 1;
        } else {
            break;
        }
    }
    AdaptiveLayout { c, chunked_len: n - plain, plain_len: plain }
}

fn adaptive_add_cost_tof(n: usize, k: usize, controlled: bool) -> u64 {
    if n == 0 {
        return 0;
    }
    let base = if controlled { 3 * n } else { 2 * n };
    let s2 = 2 * (n as f64).sqrt() as usize;
    let saved = if k >= n {
        n
    } else if k < s2 {
        (k * k) / 8
    } else {
        n / 2 + (k - s2) / 2
    };
    (base.saturating_sub(saved)) as u64
}

// ============================================================================
// controlled_chunked_then_cuccaro
// ============================================================================

fn controlled_chunked_then_cuccaro(circ: &mut B, ctrl: &QubitId, a: &[&QubitId], b: &[&QubitId], cout: Option<&QubitId>, k: usize) {
    let n = a.len();
    if n == 0 {
        return;
    }
    let cin0 = circ.alloc_qubit();
    let mut bounds: Vec<(usize, usize)> = Vec::new();
    let (mut lo, mut i) = (0usize, 0usize);
    while lo < n && k > i + 2 {
        let cc = (k - 2 - i).min(n - lo);
        bounds.push((lo, lo + cc));
        lo += cc;
        i += 1;
    }
    let chunked_len = lo;
    let mut carries: Vec<QubitId> = Vec::with_capacity(bounds.len());
    for (j, &(clo, chi)) in bounds.iter().enumerate() {
        let cy = circ.alloc_qubit();
        let cin: &QubitId = if j == 0 { &cin0 } else { &carries[j - 1] };
        controlled_clean_add_threaded(circ, ctrl, &a[clo..chi], &b[clo..chi], Some(cin), Some(&cy), chi - clo);
        carries.push(cy);
    }
    if chunked_len < n {
        let cin: &QubitId = carries.last().unwrap_or(&cin0);
        // cuccaro_carry(ctrl, x, y, ..) accumulates `y += ctrl*x` (sum into the
        // second reg, first restored). This adder computes `a += b`, so `a` must be
        // the accumulator: pass it as `y` (second), the addend `b` as `x` (first).
        let at = deref(&a[chunked_len..n]);
        let bt = deref(&b[chunked_len..n]);
        super::arith::cuccaro_carry(circ, Some(ctrl), &bt, &at, Some(cin), cout);
    } else if let Some(co) = cout {
        circ.cx(*carries.last().unwrap_or(&cin0), *co);
    }
    for j in (0..bounds.len()).rev() {
        let (clo, chi) = bounds[j];
        let carry = carries.pop().expect("carry present");
        let cin: &QubitId = if j == 0 { &cin0 } else { &carries[j - 1] };
        controlled_erase_carry_gated(circ, ctrl, &a[clo..chi], &b[clo..chi], cin, carry);
    }
    circ.zero_and_free(cin0);
}

// ============================================================================
// controlled_hybrid_add_adaptive_refs
// ============================================================================
fn controlled_hybrid_add_adaptive_refs(circ: &mut B, ctrl: &QubitId, a: &[&QubitId], b: &[&QubitId], k: usize) {
    let n = a.len();
    assert_eq!(b.len(), n, "controlled adaptive add: a,b width mismatch");
    if n == 0 {
        return;
    }
    let c = ((n as f64).sqrt() as usize).clamp(1, n);
    if n <= 4 || k.saturating_add(2 * c) >= n {
        controlled_hybrid_add_refs(circ, ctrl, a, b);
        return;
    }
    let tight = k < n.div_ceil(c) + c + ADAPTIVE_RES;
    let cov = (k.saturating_sub(2).saturating_mul(k.saturating_sub(1)) / 2).min(n);
    if tight && cov < n {
        if cov > 2 * k {
            controlled_chunked_then_cuccaro(circ, ctrl, a, b, None, k);
        } else {
            controlled_hybrid_add_refs(circ, ctrl, a, b);
        }
        return;
    }
    let cin0 = circ.alloc_qubit();
    let mut bounds: Vec<(usize, usize)> = Vec::new();
    let (l, plain_len) = if tight {
        let (mut lo, mut i) = (0usize, 0usize);
        while lo < n && k > i + 2 {
            let cc = (k - 2 - i).min(n - lo);
            bounds.push((lo, lo + cc));
            lo += cc;
            i += 1;
        }
        (n, 0)
    } else {
        let lay = adaptive_layout(n, k);
        let mut lo = 0;
        while lo < lay.chunked_len {
            let hi = (lo + lay.c).min(lay.chunked_len);
            bounds.push((lo, hi));
            lo = hi;
        }
        (lay.chunked_len, lay.plain_len)
    };
    let mut carries: Vec<QubitId> = Vec::with_capacity(bounds.len());
    for (j, &(lo, hi)) in bounds.iter().enumerate() {
        let cout = circ.alloc_qubit();
        let cin: &QubitId = if j == 0 { &cin0 } else { &carries[j - 1] };
        controlled_clean_add_threaded(circ, ctrl, &a[lo..hi], &b[lo..hi], Some(cin), Some(&cout), hi - lo);
        carries.push(cout);
    }
    if plain_len > 0 {
        let cin: &QubitId = carries.last().unwrap_or(&cin0);
        controlled_clean_add_threaded(circ, ctrl, &a[l..n], &b[l..n], Some(cin), None, plain_len);
    }
    for j in (0..bounds.len()).rev() {
        let (lo, hi) = bounds[j];
        let carry = carries.pop().expect("carry present");
        let cin: &QubitId = if j == 0 { &cin0 } else { &carries[j - 1] };
        controlled_erase_carry_gated(circ, ctrl, &a[lo..hi], &b[lo..hi], cin, carry);
    }
    circ.zero_and_free(cin0);
}

// ============================================================================
// controlled_hybrid_add_varchunk_gated_refs
// ============================================================================
fn controlled_hybrid_add_varchunk_gated_refs(circ: &mut B, ctrl: &QubitId, a: &[&QubitId], b: &[&QubitId], k: usize, cap: usize) {
    let n = a.len();
    assert_eq!(b.len(), n, "varchunk add: a,b width mismatch");
    if n == 0 {
        return;
    }
    let sizes = varchunk_schedule(n, k);
    assert!(!sizes.is_empty(), "varchunk infeasible at k={k} for n={n}");
    let cin0 = circ.alloc_qubit();
    let mut carries: Vec<QubitId> = Vec::with_capacity(sizes.len());
    let mut bounds: Vec<(usize, usize)> = Vec::with_capacity(sizes.len());
    let mut lo = 0usize;
    for (j, &s) in sizes.iter().enumerate() {
        let hi = lo + s;
        let cout = circ.alloc_qubit();
        let cin: &QubitId = if j == 0 { &cin0 } else { &carries[j - 1] };
        controlled_vented_chunk_add(circ, ctrl, &a[lo..hi], &b[lo..hi], cin, &cout);
        carries.push(cout);
        bounds.push((lo, hi));
        lo = hi;
    }
    for j in (0..sizes.len()).rev() {
        let (lo, hi) = bounds[j];
        let carry = carries.pop().expect("carry present");
        let cin: &QubitId = if j == 0 { &cin0 } else { &carries[j - 1] };
        controlled_erase_carry_gated_capped(circ, ctrl, &a[lo..hi], &b[lo..hi], cin, carry, cap);
    }
    circ.zero_and_free(cin0);
}

// ============================================================================
// controlled_hybrid_add_knob_capped_refs / controlled_hybrid_add_capped
// `k` comes from the caller (schedule).
// ============================================================================
fn controlled_hybrid_add_knob_capped_refs(circ: &mut B, ctrl: &QubitId, a: &[&QubitId], b: &[&QubitId], k: usize, cap: usize) {
    let n = a.len();
    if cap < n
        && !varchunk_schedule(n, k).is_empty()
        && (varchunk_cost(n, k, cap) as u64 + n as u64) < adaptive_add_cost_tof(n, k, true)
    {
        controlled_hybrid_add_varchunk_gated_refs(circ, ctrl, a, b, k, cap);
    } else {
        controlled_hybrid_add_adaptive_refs(circ, ctrl, a, b, k);
    }
}

/// Branch-dispatch variant: `branch` is the baked dispatch decision from the
/// schedule (0=plain, 1=varchunk, 2=adaptive; chunked-then-Cuccaro is k-selected
/// inside the adaptive path). The varchunk-vs-rest
/// choice is the cost-model-dependent knob decision -- we take it from the schedule
/// (instead of re-deciding via the cost model); the plain/adaptive/cuccaro
/// sub-choice is k-determined, so `controlled_hybrid_add_adaptive_refs(k)`
/// reproduces it.
pub fn controlled_hybrid_add_capped_branch(circ: &mut B, ctrl: &QubitId, a: &[&QubitId], b: &[&QubitId], k: usize, cap: usize, branch: u8) {
    let n = a.len();
    if branch == 1 && n > 0 && !varchunk_schedule(n, k).is_empty() {
        controlled_hybrid_add_varchunk_gated_refs(circ, ctrl, a, b, k, cap);
    } else if branch == 0 {
        // Plain (full Gidney, ~1 Toffoli/bit) -- the high-headroom branch
        // (k+2c>=n).
        controlled_hybrid_add_refs(circ, ctrl, a, b);
    } else if branch == 255 {
        // no schedule: fall back to the cost-model knob.
        controlled_hybrid_add_knob_capped_refs(circ, ctrl, a, b, k, cap);
    } else {
        controlled_hybrid_add_adaptive_refs(circ, ctrl, a, b, k);
    }
}

// ============================================================================
// controlled_hybrid_add_cout_refs  -- `ctrl ? (a += b mod 2^n)` depositing
// `ctrl AND carry_out(bit n-1)` into the caller-owned transient `cout` (|0> on
// entry). The apply cofactor-add's q-q core. `k` from the schedule.
// ============================================================================
pub fn controlled_hybrid_add_cout_refs(circ: &mut B, ctrl: &QubitId, a: &[&QubitId], b: &[&QubitId], cout: &QubitId, k: usize) {
    let n = a.len();
    assert_eq!(b.len(), n, "controlled cout add: a,b width mismatch");
    assert!(n >= 1, "controlled cout add: empty operands");
    let c = ((n as f64).sqrt() as usize).clamp(1, n);
    let lay = adaptive_layout(n, k);
    let tight = k < n.div_ceil(c) + c + ADAPTIVE_RES;
    let cov = (k.saturating_sub(2).saturating_mul(k.saturating_sub(1)) / 2).min(n);
    if n > 4 && k.saturating_add(2 * c) < n && tight && cov > 2 * k {
        controlled_chunked_then_cuccaro(circ, ctrl, a, b, Some(cout), k);
        return;
    }
    if n <= 4 || k < n.div_ceil(c) + c + ADAPTIVE_RES || k.saturating_add(2 * c) >= n || lay.plain_len == 0 {
        let zpad = circ.alloc_qubit();
        let mut aref: Vec<&QubitId> = a.to_vec();
        aref.push(cout);
        let mut bref: Vec<&QubitId> = b.to_vec();
        bref.push(&zpad);
        controlled_hybrid_add_refs(circ, ctrl, &aref, &bref);
        circ.zero_and_free(zpad);
        return;
    }
    let l = lay.chunked_len;
    let mut bounds: Vec<(usize, usize)> = Vec::new();
    let mut lo = 0;
    while lo < l {
        let hi = (lo + lay.c).min(l);
        bounds.push((lo, hi));
        lo = hi;
    }
    let mut carries: Vec<QubitId> = Vec::with_capacity(bounds.len());
    for (j, &(lo, hi)) in bounds.iter().enumerate() {
        let cy = circ.alloc_qubit();
        let cin: Option<&QubitId> = if j == 0 { None } else { Some(&carries[j - 1]) };
        controlled_clean_add_threaded(circ, ctrl, &a[lo..hi], &b[lo..hi], cin, Some(&cy), hi - lo);
        carries.push(cy);
    }
    controlled_clean_add_threaded(circ, ctrl, &a[l..n], &b[l..n], carries.last(), Some(cout), lay.plain_len);
    for j in (0..bounds.len()).rev() {
        let (lo, hi) = bounds[j];
        let carry = carries.pop().expect("carry present");
        if j == 0 {
            let cin0 = circ.alloc_qubit();
            controlled_erase_carry_gated(circ, ctrl, &a[lo..hi], &b[lo..hi], &cin0, carry);
            circ.zero_and_free(cin0);
        } else {
            controlled_erase_carry_gated(circ, ctrl, &a[lo..hi], &b[lo..hi], &carries[j - 1], carry);
        }
    }
}
