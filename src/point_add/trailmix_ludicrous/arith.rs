//! Modular-arithmetic primitives the product-min secp256k1 EC-add uses,
//! built on this crate's `B` builder: explicit qubit alloc/free, raw
//! measurement-vented carries (hmr + cz_if_bit), no higher-level ancilla-management
//! or phase-tracker machinery.
//!
//! ## secp256k1 product-min constants (hardcoded)
//!   q   = 2^256 - f,   f = 2^32 + 977   (bitlen(f) = 33)
//!   PAD = 21
//!   +f window width  lsbs = PAD + bitlen(f) = 54   (carry beyond dropped)
//!   less-than comparator   msbs = PAD = 21          (top-k less-than)
//!
//! ## Normal uncompute vs MBU (measurement-based uncomputation) choice (per primitive)
//! - register add (the n-bit `y += ctrl*x`): Cuccaro ripple, uncomputed normally. The
//!   uncompute is an exact gate-inverse (no measurement), so it is already
//!   ancilla-clean; the +1 carry ancilla is freed. Choosing the MBU
//!   Gidney vent here would only trade Toffolis for measurement on the live
//!   qubits.
//! - `+f` fold (`add_f_window`): MBU Gidney clean
//!   constant adder. The carry chain feeds each carry into the next carry's
//!   AND, so the carries CANNOT be uncomputed normally in place without a
//!   second forward pass; the clean Gidney form instead measures (`hmr`)
//!   each carry and corrects its deferred phase with a `CZ` gated on the hmr
//!   bit. Cost: +(lsbs-1) clean carry qubits, all freed.
//! - `mod_double`: normal shift + the (MBU) `+f` fold. The shift and
//!   the ancilla cleanup are normal CX; only the shared `+f` fold vents.

use super::{B, BExt};
use crate::circuit::{BitId, QubitId};

/// secp256k1 reduction constant f = 2^256 - q.
pub const F_SECP256K1: u64 = (1u64 << 32) + 977;
/// bitlen(f).
pub const F_BITLEN: usize = 33;
/// Profile padding.
pub const PAD: usize = 21;
/// `+f` fold window width: carry beyond bit `LSBS-1` is dropped (~2^-PAD miss).
pub const LSBS: usize = PAD + F_BITLEN; // 54
/// Top-k less-than comparator width for the mod-add/sub overflow cleanup.
pub const MSBS: usize = PAD; // 21
/// Chunk width for the measurement-vented chunked-gated register adder used on the
/// peak-bound apply path. Sized so the per-call working set (one chunk's `W`
/// carries + the `n/W` boundary carries + the erase comparator's `W` carries) fits
/// the apply's ~1170 ceiling headroom over the 514+603 resident base.
pub const APPLY_CHUNK: usize = 40;

#[inline]
fn cbit(c: &[u8], i: usize) -> bool {
    let byte = i / 8;
    byte < c.len() && (c[byte] >> (i % 8)) & 1 == 1
}

// =====================================================================
// Single-ancilla controlled Cuccaro ripple-add (the LOW-PEAK adder)
// =====================================================================

/// `y += ctrl * x  (mod 2^s)` over `s = y.len()` bits, riding the carry through
/// one ancilla (the addend `x` is used only as a control and is restored). When
/// `cout` is `Some`, deposits the overall carry-out (`ctrl AND carry(bit s-1)`)
/// into it. `cin`, if `Some`, is a caller-owned carry-in (restored); else a fresh
/// |0> ancilla ripples and is freed. `ctrl = None` makes it an UNCONTROLLED add.
///
/// The controlled adder the low-product apply / GCD-subtract use; 1 live
/// ancilla. ~3 Toffoli/bit controlled (forward MAJ + reverse MAJ-restore + gated
/// sum); ~2/bit uncontrolled (the gated sum degenerates to a CX).
pub fn cuccaro_carry(
    circ: &mut B,
    ctrl: Option<&QubitId>,
    x: &[QubitId],
    y: &[QubitId],
    cin: Option<&QubitId>,
    cout: Option<&QubitId>,
) {
    let s = y.len();
    assert_eq!(x.len(), s, "cuccaro_carry: x,y width mismatch");
    let fresh = if cin.is_none() { Some(circ.alloc_qubit()) } else { None };
    let c: &QubitId = cin.unwrap_or_else(|| fresh.as_ref().unwrap());
    let sum = |circ: &mut B, xi: &QubitId, yi: &QubitId| match ctrl {
        Some(ct) => circ.ccx(*ct, *xi, *yi), // y ^= ctrl*(x ^ c_in) = gated sum
        None => circ.cx(*xi, *yi),
    };
    let gated_carry = |circ: &mut B, co: &QubitId| match ctrl {
        Some(ct) => circ.ccx(*ct, *c, *co),
        None => circ.cx(*c, *co),
    };
    if s == 0 {
        if let Some(co) = cout {
            gated_carry(circ, co);
        }
    } else {
        // Forward MAJ ripple: `c` is the running carry of x + y + cin.
        for i in 0..s {
            circ.cx(*c, y[i]);
            circ.cx(*c, x[i]);
            circ.ccx(x[i], y[i], *c);
        }
        if let Some(co) = cout {
            gated_carry(circ, co); // c holds carry-out(bit s-1)
        }
        // Reverse: restore each carry, write the (gated) sum into y, restore x.
        for i in (0..s).rev() {
            circ.ccx(x[i], y[i], *c);
            circ.cx(*c, y[i]);
            sum(circ, &x[i], &y[i]);
            circ.cx(*c, x[i]);
        }
    }
    if let Some(f) = fresh {
        circ.zero_and_free(f);
    }
}

/// Shared body of [`controlled_clean_add_threaded`] (`ctrl = Some`) and
/// [`add_threaded`] (`ctrl = None`, sums via CX, ~1 Toffoli/bit).
fn clean_add_threaded_opt(
    circ: &mut B,
    ctrl: Option<&QubitId>,
    x: &[QubitId],
    y: &[QubitId],
    cin: Option<&QubitId>,
    cout: Option<&QubitId>,
) {
    let s = y.len();
    assert_eq!(x.len(), s, "vented add: x,y width mismatch");
    // gated sum: y ^= ctrl?(x_i) -- ccx if controlled, cx if not.
    let gated_sum = |circ: &mut B, xi: &QubitId, yi: &QubitId| match ctrl {
        Some(ct) => circ.ccx(*ct, *xi, *yi),
        None => circ.cx(*xi, *yi),
    };
    if s == 0 {
        if let (Some(ci), Some(co)) = (cin, cout) {
            match ctrl {
                Some(ct) => circ.ccx(*ct, *ci, *co),
                None => circ.cx(*ci, *co),
            }
        }
        return;
    }
    let n_inner = if cout.is_some() { s } else { s - 1 };
    let mut inner: Vec<Option<QubitId>> = (0..n_inner).map(|_| Some(circ.alloc_qubit())).collect();
    let produces = |i: usize| cout.is_some() || i + 1 < s;
    // Forward MAJ (UNCONDITIONAL): carry of x + y (+ cin).
    for i in 0..s {
        if !produces(i) {
            continue;
        }
        let co = inner[i].as_ref().unwrap();
        let ci: Option<&QubitId> = if i == 0 { cin } else { inner[i - 1].as_ref() };
        if let Some(ci) = ci {
            circ.cx(*ci, x[i]);
            circ.cx(*ci, y[i]);
            circ.ccx(x[i], y[i], *co);
            circ.cx(*ci, *co);
        } else {
            circ.ccx(x[i], y[i], *co);
        }
    }
    if let Some(cout) = cout {
        let top = inner[s - 1].as_ref().unwrap();
        match ctrl {
            Some(ct) => circ.ccx(*ct, *top, *cout),
            None => circ.cx(*top, *cout),
        }
    }
    // Reverse: gated sums, vent every internal carry. The forward folded both
    // x[i] and y[i] by ci; here y (the accumulator) is unfolded before the gated
    // sum and x (the addend) restored after.
    for i in (0..s).rev() {
        if !produces(i) {
            // Top mod bit (no carry-out): sum y ^= [ctrl &] (x_i ^ ci). x folded by
            // ci just for the sum, then unfolded.
            let ci: Option<&QubitId> = if i == 0 { cin } else { inner[i - 1].as_ref() };
            if let Some(ci) = ci {
                circ.cx(*ci, x[i]);
            }
            gated_sum(circ, &x[i], &y[i]);
            if let Some(ci) = ci {
                circ.cx(*ci, x[i]);
            }
            continue;
        }
        let co = inner[i].take().unwrap();
        let ci: Option<&QubitId> = if i == 0 { cin } else { inner[i - 1].as_ref() };
        if let Some(ci) = ci {
            circ.cx(*ci, co); // co = x[i] & y[i] (re-fold for the AND identity)
        }
        // Vent the carry AND: hmr(co) + cz_if_bit(x[i], y[i]) gated on the hmr bit.
        let bit = circ.alloc_bit();
        circ.hmr(co, bit);
        circ.zero_and_free(co);
        circ.cz_if_bit(x[i], y[i], bit);
        if let Some(ci) = ci {
            circ.cx(*ci, y[i]); // unfold the accumulator y (remove forward fold)
        }
        gated_sum(circ, &x[i], &y[i]); // y ^= [ctrl &] (x_i ^ ci) -> sum
        if let Some(ci) = ci {
            circ.cx(*ci, x[i]); // restore the addend x
        }
    }
}

/// Gated VENTED erase of an inter-chunk boundary carry: `carry` holds
/// `[ctrl AND] carryout(a + b + cin)`. HMR it (0 Toffoli, resets to |0>), then on
/// the ~50% kickback shots recompute the predicate as a deferred Z. Net: carry ->
/// |0>, phase-clean. `carryout(a+b+cin) = 1 ^ (ta&tb) ^ c_prev` (the complement of
/// the [`compare_geq_cin_middle`] built carry), so deposit (with `ctrl = Some`)
/// `Z^(ctrl AND (1 ^ (ta&tb) ^ c_prev))` = `Z(ctrl) ^ CCZ(ctrl,ta,tb) ^ CZ(ctrl,c_prev)`,
/// or (uncontrolled) `Z^(1 ^ (ta&tb) ^ c_prev)` = `neg ^ CZ(ta,tb) ^ Z(c_prev)`.
///
/// PAD-capped: when the chunk width exceeds `cap`, the boundary predicate
/// is recomputed from only the top `cap` bits with a fresh |0> carry-in -- a
/// phase truncation that mis-clears with probability ~2^-cap (the schedule's
/// PAD term accounts for it). `cap = None` = exact full-width erase.
pub(crate) fn erase_carry_gated_opt(
    circ: &mut B,
    ctrl: Option<&QubitId>,
    a: &[QubitId],
    b: &[QubitId],
    cin: &QubitId,
    carry: &QubitId,
    cap: Option<usize>,
) {
    let s = a.len();
    let bit = circ.alloc_bit();
    circ.hmr(*carry, bit);
    circ.push_condition(bit);
    let deposit = |c: &mut B, ta: &QubitId, tb: &QubitId, c_prev: &QubitId| match ctrl {
        Some(ct) => {
            c.z(*ct);
            c.ccz(*ct, *ta, *tb);
            c.cz(*ct, *c_prev);
        }
        None => {
            c.neg(); // the constant 1
            c.cz(*ta, *tb);
            c.z(*c_prev);
        }
    };
    match cap {
        Some(k) if k < s => {
            // PAD-truncated: recompute on the top `k` bits with a |0> carry-in.
            let lo = s - k;
            let zcin = circ.alloc_qubit();
            super::comparator::compare_geq_cin_middle(circ, &a[lo..], &b[lo..], &zcin, deposit);
            circ.zero_and_free(zcin);
        }
        _ => {
            super::comparator::compare_geq_cin_middle(circ, a, b, cin, deposit);
        }
    }
    circ.pop_condition();
}

/// A vented chunked add that optionally captures the overall
/// carry-out (`cout ^= ctrl AND carryout(x+y)`) into `cout` (|0> on entry). Exact
/// boundary erases (`cap = None`): the apply mod-add/sub operands are generic field
/// elements whose carries propagate the full chunk width, so the boundary predicate
/// must be recomputed in full (the PAD cap is only valid where the carry RUN is
/// schedule-bounded -- the GCD subtract, which here uses the erase-free
/// full-threaded adder instead).
pub fn controlled_add_vented_chunked_cout(
    circ: &mut B,
    ctrl: &QubitId,
    x: &[QubitId],
    y: &[QubitId],
    chunk: usize,
    cout: Option<&QubitId>,
) {
    add_vented_chunked_opt(circ, Some(ctrl), x, y, chunk, cout, None);
}

/// The product-min operating ceiling. The adder's vent budget =
/// `CEILING - active_qubits` -- the available qubit headroom the adaptive adder
/// fills.
pub const CEILING: usize = 1167;

/// Shared chunk-emit body: clean-threaded low chunks (each holding its boundary
/// carry) + a plain comparator-free top region + reverse `cap`-capped boundary
/// erases. Factored out of [`add_vented_chunked_opt`] so the GCD adaptive layout
/// and the apply chunked layout share one gate-emitter.
fn emit_chunked_capped(
    circ: &mut B,
    ctrl: Option<&QubitId>,
    x: &[QubitId],
    y: &[QubitId],
    bounds: &[(usize, usize)],
    plain_len: usize,
    cout: Option<&QubitId>,
    cap: Option<usize>,
) {
    let n = y.len();
    let l = n - plain_len; // chunked low-region length
    let cin0 = circ.alloc_qubit();
    let mut carries: Vec<QubitId> = Vec::with_capacity(bounds.len());
    for (j, &(lo, hi)) in bounds.iter().enumerate() {
        let cy = circ.alloc_qubit();
        let cin: &QubitId = if j == 0 { &cin0 } else { &carries[j - 1] };
        clean_add_threaded_opt(circ, ctrl, &x[lo..hi], &y[lo..hi], Some(cin), Some(&cy));
        carries.push(cy);
    }
    if l < n {
        let top_cin: &QubitId = carries.last().unwrap_or(&cin0);
        clean_add_threaded_opt(circ, ctrl, &x[l..n], &y[l..n], Some(top_cin), cout);
    } else if let Some(co) = cout {
        circ.cx(*carries.last().unwrap(), *co);
    }
    for j in (0..bounds.len()).rev() {
        let (lo, hi) = bounds[j];
        let carry = carries.pop().expect("carry present");
        let cin: &QubitId = if j == 0 { &cin0 } else { &carries[j - 1] };
        erase_carry_gated_opt(circ, ctrl, &y[lo..hi], &x[lo..hi], cin, &carry, cap);
        circ.zero_and_free(carry);
    }
    circ.zero_and_free(cin0);
}

/// UNCONTROLLED plain Gidney measurement-vented add `a += b mod 2^n` (b restored),
/// the first `vents = min(vents_budget, n-1)` carries measurement-vented (HMR + a
/// gated `cz`) and the rest normally uncomputed. The vented AND-carries use bare
/// `hmr`/`cz_if_bit`; sums are UNCONDITIONAL `cx`. `vents` is fixed -- so this is
/// fully determined by `(n, vents_budget)`, no schedule needed.
fn hybrid_add_plain(circ: &mut B, a: &[QubitId], b: &[QubitId], vents_budget: usize) {
    let n = a.len();
    assert_eq!(b.len(), n, "hybrid_add: a,b width mismatch");
    if n == 0 {
        return;
    }
    if n == 1 {
        circ.cx(b[0], a[0]);
        return;
    }
    let vents = vents_budget.min(n - 1);
    for i in 1..n {
        circ.cx(b[i], a[i]);
    }
    for i in (1..n - 1).rev() {
        circ.cx(b[i], b[i + 1]);
    }
    let mut vent_ancs: Vec<Option<QubitId>> = (0..n - 1).map(|_| None).collect();
    for i in 0..n - 1 {
        if i < vents {
            let anc = circ.alloc_qubit();
            circ.ccx(a[i], b[i], anc);
            circ.cx(anc, b[i + 1]);
            vent_ancs[i] = Some(anc);
        } else {
            circ.ccx(a[i], b[i], b[i + 1]);
        }
    }
    for i in (0..n - 1).rev() {
        circ.cx(b[i + 1], a[i + 1]); // UNCONDITIONAL sum bit i+1
        if i < vents {
            let anc = vent_ancs[i].take().unwrap();
            circ.cx(anc, b[i + 1]);
            let bit = circ.alloc_bit();
            circ.hmr(anc, bit);
            circ.zero_and_free(anc);
            circ.cz_if_bit(a[i], b[i], bit);
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
}

/// UNCONTROLLED exact adaptive add `a += b mod 2^n` (b restored) at headroom `k`.
/// Degenerates to the plain Gidney add when headroom is ample (`k + 2c >= n`) or
/// `n <= 4`; otherwise a sqrt(n)-chunked low region (boundary carries gated-erased)
/// + a comparator-free plain top region sized by [`super::gidney::adaptive_layout`].
/// `k` is the exact qubit headroom from the baked square row-add schedule -- NO
/// `active_qubits` read, NO cap. Reuses [`emit_chunked_capped`] (the same chunked
/// gate-emitter the apply uses) for the adaptive_layout branch.
pub(crate) fn hybrid_add_adaptive(circ: &mut B, a: &[QubitId], b: &[QubitId], k: usize) {
    let n = a.len();
    assert_eq!(b.len(), n, "adaptive add: a,b width mismatch");
    if n == 0 {
        return;
    }
    let c = ((n as f64).sqrt() as usize).clamp(1, n);
    if n <= 4 || k.saturating_add(2 * c) >= n {
        hybrid_add_plain(circ, a, b, k);
        return;
    }
    if k < n.div_ceil(c) + c + super::gidney::ADAPTIVE_RES {
        let cov = (k.saturating_mul(k.saturating_sub(1)) / 2).min(n);
        if cov > 2 * k {
            // The tight chunked-then-cuccaro branch. Unreachable for the
            // square: k is the fixed square headroom (~130) and n <= 258, so the
            // tight gate `k < ~2*sqrt(n)` never fires. Panic if a schedule
            // change ever violates this assumption.
            unreachable!("square adaptive add hit the tight chunked_then_cuccaro branch (n={n}, k={k})");
        }
        hybrid_add_plain(circ, a, b, k);
        return;
    }
    let lay = super::gidney::adaptive_layout(n, k);
    let l = lay.chunked_len;
    let mut bounds: Vec<(usize, usize)> = Vec::new();
    let mut lo = 0;
    while lo < l {
        let hi = (lo + lay.c).min(l);
        bounds.push((lo, hi));
        lo = hi;
    }
    // a += b: accumulator = a (y), addend = b (x). emit_chunked_capped(x=addend,
    // y=accumulator); cout=None (mod 2^n), cap=None (exact erases, integer square).
    emit_chunked_capped(circ, None, b, a, &bounds, lay.plain_len, None, None);
}

fn add_vented_chunked_opt(
    circ: &mut B,
    ctrl: Option<&QubitId>,
    x: &[QubitId],
    y: &[QubitId],
    chunk: usize,
    cout: Option<&QubitId>,
    cap: Option<usize>,
) {
    add_vented_chunked_opt_capped(circ, ctrl, x, y, chunk, cout, cap, usize::MAX);
}

#[allow(clippy::too_many_arguments)]
fn add_vented_chunked_opt_capped(
    circ: &mut B,
    ctrl: Option<&QubitId>,
    x: &[QubitId],
    y: &[QubitId],
    chunk: usize,
    cout: Option<&QubitId>,
    cap: Option<usize>,
    max_vents: usize,
) {
    let n = y.len();
    assert_eq!(x.len(), n, "chunked add: x,y width mismatch");
    if n == 0 {
        return;
    }
    // Headroom-adaptive layout -- the decision hardcoded as `CEILING - live`.
    // The top `plain_len` bits are a single plain
    // threaded add whose internal carries are measurement-vented and whose overall
    // carry goes to `cout` -- NO normal boundary erase. The low bits are chunked
    // into width-`c` blocks whose boundary carries are normally erased. Held carries
    // = `ceil(chunked_len/c) + plain_len`, balanced to the vent budget `k`, so
    // `plain_len = (k*c - n)/(c-1)`. More headroom => bigger plain region => fewer
    // erases; an early GCD subtract (large headroom) is fully plain (zero erases).
    let c = chunk.clamp(1, n); // low-region block width
    let live = circ.active_qubits as usize;
    // vent budget (held-carry cap); capped at `max_vents` (e.g. ROW_ADD_VENTS for the
    // square row-adds).
    let k = CEILING.saturating_sub(live).clamp(1, n).min(max_vents);
    let plain_len = if k >= n {
        n
    } else if c <= 1 {
        0
    } else {
        ((k * c).saturating_sub(n) / (c - 1)).min(n)
    };
    let l = n - plain_len; // chunked low-region length
    let mut bounds: Vec<(usize, usize)> = Vec::new();
    let mut lo = 0;
    while lo < l {
        let hi = (lo + c).min(l);
        bounds.push((lo, hi));
        lo = hi;
    }
    emit_chunked_capped(circ, ctrl, x, y, &bounds, plain_len, cout, cap);
}

// =====================================================================
// Clean (measurement-vented) Gidney +f constant fold
// =====================================================================

// Conditional CCX helper: optionally `cx(ctrl, c1)` / `cx(ctrl, c2)` around a
// `ccx(c1, c2, t)`, gated by `b0`/`b1` (constant-bit flags). The wrapping CX
// pairs cancel, so each is applied iff its flag is set.
fn ccx_cond(circ: &mut B, ctrl: &QubitId, c1: &QubitId, c2: &QubitId, t: &QubitId, b0: bool, b1: bool) {
    if b0 { circ.cx(*ctrl, *c1); }
    if b1 { circ.cx(*ctrl, *c2); }
    circ.ccx(*c1, *c2, *t);
    if b0 { circ.cx(*ctrl, *c1); }
    if b1 { circ.cx(*ctrl, *c2); }
}

// Recompute the n-1 carries of `a + ctrl*c[off..]` (a = complemented sum) and XOR them
// back out of the borrowed dirty bits `out`, restoring them. carry-in = `cin` (read-only).
fn xor_carries_off_cin(circ: &mut B, ctrl: &QubitId, a: &[QubitId], c: &[u8], off: usize, out: &[QubitId], cin: &QubitId) {
    let n = a.len();
    for i in (1..n - 1).rev() {
        ccx_cond(circ, ctrl, &a[i], &out[i - 1], &out[i], cbit(c, off + i), false);
    }
    for i in 0..n - 1 {
        if cbit(c, off + i) { circ.cx(*ctrl, out[i]); }
    }
    ccx_cond(circ, ctrl, cin, &a[0], &out[0], cbit(c, off), cbit(c, off));
    for i in 1..n - 1 {
        ccx_cond(circ, ctrl, &a[i], &out[i - 1], &out[i], cbit(c, off + i), cbit(c, off + i));
    }
}

// Borrowed-dirty controlled const add `a += ctrl*c[off..off+n] (mod 2^n)` with carry-IN
// `cin` (read as cy_0, NOT freed -- owned by the caller). Carries vented via hmr->bit and
// Z-discharged (z_if_bit(dirty[i]) before+after the xor_carries restore = Z^(bit&carry)).
fn dirty_carryin(circ: &mut B, ctrl: &QubitId, a: &[QubitId], c: &[u8], off: usize, dirty: &[QubitId], cin: &QubitId) {
    let n = a.len();
    debug_assert!(n >= 2 && dirty.len() >= n - 1);
    let mut bits: Vec<BitId> = Vec::with_capacity(n - 1);
    let mut cy_owned: Option<QubitId> = None;
    for i in 0..(n - 1) {
        let new = circ.alloc_qubit();
        let anc = circ.alloc_qubit();
        let on = cbit(c, off + i);
        let cyref: QubitId = match cy_owned { Some(q) => q, None => *cin };
        if on { circ.cx(*ctrl, anc); }
        circ.cx(cyref, anc);
        circ.cx(cyref, a[i]);
        circ.ccx(a[i], anc, new);
        circ.cx(cyref, new);
        circ.cx(new, dirty[i]);
        circ.cx(cyref, anc);
        if on { circ.cx(*ctrl, anc); circ.cx(*ctrl, a[i]); }
        circ.zero_and_free(anc);
        if let Some(old) = cy_owned.take() {
            let b = circ.alloc_bit();
            circ.hmr(old, b);
            bits.push(b);
            circ.zero_and_free(old);
        }
        cy_owned = Some(new);
    }
    let cy_top = cy_owned.take().unwrap();
    if cbit(c, off + n - 1) { circ.cx(*ctrl, a[n - 1]); }
    circ.cx(cy_top, a[n - 1]);
    {
        let b = circ.alloc_bit();
        circ.hmr(cy_top, b);
        bits.push(b);
    }
    circ.zero_and_free(cy_top);
    for i in 0..(n - 1) { circ.z_if_bit(dirty[i], bits[i]); }
    for q in a { circ.x(*q); }
    xor_carries_off_cin(circ, ctrl, a, c, off, dirty, cin);
    for q in a { circ.x(*q); }
    for i in 0..(n - 1) { circ.z_if_bit(dirty[i], bits[i]); }
}

// Hybrid: clean prefix carries [0,k) + borrowed-dirty suffix [k,n). cy_k carries into
// the suffix; the clean prefix carries are vented (hmr) with a CZ discharge (clean carry
// = ta & tb, an AND of two live qubits). k = clean carry count (1 <= k <= n-2).
// ===================================================================
// Graduated chunked-gated const add (the +f hybrid suffix path). Used for large
// suffixes where the borrowed-dirty path would otherwise be hit.
// ===================================================================

fn graduated_const_fits(n: usize, k: usize) -> bool {
    k >= 4 && (k - 3) * (k - 2) / 2 >= n
}
fn graduated_const_kmin(n: usize) -> usize {
    (4..).find(|&k| graduated_const_fits(n, k)).unwrap()
}

/// Native clean const carry chain over a chunk: `a += ctrl*c[coff..] mod 2^s`, carry-OUT
/// kept in `cout`, internal carries measurement-vented.
fn const_chunk_add_clean(circ: &mut B, ctrl: &QubitId, a: &[QubitId], c: &[u8], coff: usize, cin: &QubitId, cout: &QubitId) {
    let s = a.len();
    if s == 0 {
        return;
    }
    let mut int: Vec<Option<QubitId>> = (0..s - 1).map(|_| Some(circ.alloc_qubit())).collect();
    for i in 0..s {
        let on = cbit(c, coff + i);
        let cin_ref: QubitId = if i == 0 { *cin } else { *int[i - 1].as_ref().unwrap() };
        let cout_ref: QubitId = if i == s - 1 { *cout } else { *int[i].as_ref().unwrap() };
        circ.cx(cin_ref, a[i]);
        if on {
            circ.cx(*ctrl, cin_ref);
        }
        circ.ccx(a[i], cin_ref, cout_ref);
        if on {
            circ.cx(*ctrl, cin_ref);
        }
        circ.cx(cin_ref, cout_ref);
    }
    for i in 0..s {
        if cbit(c, coff + i) {
            circ.cx(*ctrl, a[i]);
        }
    }
    for i in (0..s - 1).rev() {
        let on = cbit(c, coff + i);
        let int_i = int[i].take().unwrap();
        let cin_ref: QubitId = if i == 0 { *cin } else { *int[i - 1].as_ref().unwrap() };
        if on {
            circ.cx(*ctrl, a[i]);
        }
        circ.cx(cin_ref, int_i);
        if on {
            circ.cx(*ctrl, cin_ref);
        }
        let b = circ.alloc_bit();
        circ.hmr(int_i, b);
        circ.zero_and_free(int_i);
        circ.cz_if_bit(a[i], cin_ref, b);
        if on {
            circ.cx(*ctrl, cin_ref);
            circ.cx(*ctrl, a[i]);
        }
    }
}

/// No-temp const carry comparator with a middle callback handing `(a_top, cy_top,
/// const_top)`.
fn compare_geq_const_cin_middle<F: FnOnce(&mut B, &QubitId, &QubitId, bool)>(circ: &mut B, a: &[QubitId], c: &[u8], coff: usize, cin: &QubitId, body: F) {
    let s = a.len();
    let mut cy: Vec<Option<QubitId>> = Vec::with_capacity(s);
    let c0 = circ.alloc_qubit();
    circ.x(c0);
    circ.cx(*cin, c0);
    cy.push(Some(c0));
    for i in 0..s - 1 {
        let on = cbit(c, coff + i);
        let next = circ.alloc_qubit();
        let ci = *cy[i].as_ref().unwrap();
        circ.ccx(a[i], ci, next);
        if !on {
            circ.cx(a[i], next);
            circ.cx(ci, next);
        }
        cy.push(Some(next));
    }
    {
        let i = s - 1;
        let on = cbit(c, coff + i);
        let ci = *cy[i].as_ref().unwrap();
        body(circ, &a[i], &ci, on);
    }
    for i in (0..s - 1).rev() {
        let on = cbit(c, coff + i);
        let next = cy[i + 1].take().unwrap();
        let ci = *cy[i].as_ref().unwrap();
        if !on {
            circ.cx(ci, next);
            circ.cx(a[i], next);
        }
        let b = circ.alloc_bit();
        circ.hmr(next, b);
        circ.zero_and_free(next);
        circ.cz_if_bit(a[i], ci, b);
    }
    let c0 = cy[0].take().unwrap();
    circ.cx(*cin, c0);
    circ.x(c0);
    circ.zero_and_free(c0);
}

/// Gated-erase the inter-chunk carry against the classical chunk constant: hmr the
/// carry, then deposit the predicate phase under a push_condition on the hmr bit.
fn controlled_erase_carry_gated_const(circ: &mut B, ctrl: &QubitId, a: &[QubitId], c: &[u8], coff: usize, cin: &QubitId, carry: QubitId) {
    let bit = circ.alloc_bit();
    circ.hmr(carry, bit);
    circ.push_condition(bit);
    compare_geq_const_cin_middle(circ, a, c, coff, cin, |cc, a_top, cy_top, ctop| {
        // Z^(ctrl . NOT cy_s); const_top=1 (AND): cy_s=a&cy; =0 (OR): cy_s=a|cy.
        cc.z(*ctrl);
        cc.ccz(*ctrl, *a_top, *cy_top);
        if !ctop {
            cc.cz(*ctrl, *a_top);
            cc.cz(*ctrl, *cy_top);
        }
    });
    circ.pop_condition();
    circ.zero_and_free(carry);
}

/// GRADUATED staircase const add on a suffix: chunk `i` width `k-3-i`.
fn controlled_add_const_chunked_graduated_off(circ: &mut B, ctrl: &QubitId, a: &[QubitId], c: &[u8], coff: usize, cin: &QubitId, k: usize) {
    let n = a.len();
    if n == 0 {
        return;
    }
    let mut bounds: Vec<(usize, usize)> = Vec::new();
    let (mut lo, mut i) = (0usize, 0usize);
    while lo < n && k > i + 3 {
        let cc = (k - 3 - i).min(n - lo);
        bounds.push((lo, lo + cc));
        lo += cc;
        i += 1;
    }
    assert_eq!(lo, n, "graduated staircase (k={k}) covers {lo} < n={n}");
    let mut carries: Vec<QubitId> = Vec::with_capacity(bounds.len());
    for (j, &(clo, chi)) in bounds.iter().enumerate() {
        let cout = circ.alloc_qubit();
        let cin_ref: QubitId = if j == 0 { *cin } else { carries[j - 1] };
        const_chunk_add_clean(circ, ctrl, &a[clo..chi], c, coff + clo, &cin_ref, &cout);
        carries.push(cout);
    }
    for j in (0..bounds.len()).rev() {
        let (clo, chi) = bounds[j];
        let carry = carries.pop().expect("carry present");
        let cin_ref: QubitId = if j == 0 { *cin } else { carries[j - 1] };
        controlled_erase_carry_gated_const(circ, ctrl, &a[clo..chi], c, coff + clo, &cin_ref, carry);
    }
}

/// `reg[..lsbs] += ctrl * c (mod 2^lsbs)` via the CLEAN measurement-vented
/// Gidney constant adder. Carry-out of bit `lsbs-1` is dropped (matches the
/// ludicrous `+f` window's ~2^-PAD approximation). Allocates `lsbs-1` clean
/// carry qubits, measurement-vents (`hmr`) each on the reverse pass, and
/// corrects its deferred phase with a `CZ`/`CCZ` gated on the hmr bit.
///
/// MBU by necessity: the nested carry chain (each carry feeds the next
/// AND) has no in-place normal uncompute.
///
/// Borrowed-bit hybrid: clean measurement-vented carries for the prefix `[0,k)`
/// up to the live headroom, borrowed-dirty register bits for the overflow suffix.
#[allow(clippy::needless_range_loop)] // i indexes a[] while reading the constant c at bit i
fn add_f_window_hybrid(circ: &mut B, ctrl: &QubitId, reg: &[QubitId], lsbs: usize, c: &[u8], k: usize) {
    let n = lsbs;
    let a: Vec<QubitId> = reg[..n].to_vec();
    let suf_dirty = n - k - 1; // suffix has n-k bits -> n-k-1 dirty
    assert!(reg.len() >= lsbs + suf_dirty, "+f hybrid: not enough high bits to borrow");
    let dirty: Vec<QubitId> = (lsbs..lsbs + suf_dirty).map(|i| reg[i]).collect();
    let mut cy: Vec<Option<QubitId>> = (0..k).map(|_| Some(circ.alloc_qubit())).collect();
    // prefix forward bits 0..k (clean ta/tb chain).
    if cbit(c, 0) { circ.ccx(*ctrl, a[0], *cy[0].as_ref().unwrap()); }
    for i in 1..k {
        let ci = *cy[i - 1].as_ref().unwrap();
        let next = *cy[i].as_ref().unwrap();
        circ.cx(ci, a[i]);
        if cbit(c, i) { circ.cx(*ctrl, ci); }
        circ.ccx(a[i], ci, next);
        if cbit(c, i) { circ.cx(*ctrl, ci); }
        circ.cx(ci, next);
    }
    // prefix sums bits 0..k.
    for i in 0..k { if cbit(c, i) { circ.cx(*ctrl, a[i]); } }
    // suffix [k,n): carry-in cy_k = cy[k-1]. Graduated-chunked when the remaining
    // headroom covers it (avoids the large-suffix borrowed-dirty path), else
    // borrowed-dirty.
    {
        let a_hi: Vec<QubitId> = a[k..].to_vec();
        let cin = *cy[k - 1].as_ref().unwrap();
        let sn = n - k;
        // Graduated-chunked suffix. With `g` (the prefix size) schedule-driven
        // (deterministic), graduated runs with fixed params -> phase-clean. kmin
        // always fits. sn < 2 falls to the borrowed-dirty path.
        if sn >= 2 {
            controlled_add_const_chunked_graduated_off(circ, ctrl, &a_hi, c, k, &cin, graduated_const_kmin(sn));
        } else {
            dirty_carryin(circ, ctrl, &a_hi, c, k, &dirty, &cin);
        }
    }
    // prefix reverse bits k-1..1: vent cy[i] (=carry_{i+1}) via CZ discharge.
    for i in (1..k).rev() {
        if cbit(c, i) { circ.cx(*ctrl, a[i]); } // sum_i -> ta_i
        let ci = *cy[i - 1].as_ref().unwrap();
        let next = *cy[i].as_ref().unwrap();
        circ.cx(ci, next); // next = ta & tb (undo cy XOR)
        if cbit(c, i) { circ.cx(*ctrl, ci); } // ci -> tb_i
        let nq = cy[i].take().unwrap();
        let b = circ.alloc_bit();
        circ.hmr(nq, b);
        circ.zero_and_free(nq);
        circ.cz_if_bit(a[i], ci, b);
        if cbit(c, i) { circ.cx(*ctrl, ci); circ.cx(*ctrl, a[i]); }
    }
    // reverse bit 0: vent cy[0] = carry_1 = a_0 & ctrl.
    let cy0 = cy[0].take().unwrap();
    if cbit(c, 0) {
        circ.cx(*ctrl, a[0]);
        let b = circ.alloc_bit();
        circ.hmr(cy0, b);
        circ.zero_and_free(cy0);
        circ.cz_if_bit(a[0], *ctrl, b);
        circ.cx(*ctrl, a[0]);
    } else {
        circ.zero_and_free(cy0);
    }
}

/// `reg[..lsbs] += ctrl * c (mod 2^lsbs)` -- the +f overflow fold, headroom-adaptive:
/// `g = min(CEILING-live, n-1)` clean measurement-vented carries + borrowed-dirty
/// register bits for the overflow, so the fold never inflates the peak past the
/// ceiling.
fn add_f_window(circ: &mut B, ctrl: &QubitId, reg: &[QubitId], lsbs: usize, c: &[u8], g_sched: Option<usize>) {
    let n = lsbs;
    assert!(n <= reg.len(), "register too short for +f window");
    if n == 0 { return; }
    if n == 1 {
        if cbit(c, 0) { circ.cx(*ctrl, reg[0]); }
        return;
    }
    // `g` = clean-vent count. Schedule-driven for the apply cofactor folds (g_sched =
    // the baked `g`, deterministic -> phase-clean); live-headroom read (CEILING - live)
    // for the doublings (always clean -- high headroom).
    let g = g_sched.map_or_else(|| CEILING.saturating_sub(circ.active_qubits as usize), |g| g).min(n - 1);
    if g >= n - 1 {
        add_f_window_clean(circ, ctrl, reg, lsbs, c); // all-clean path
    } else if g == 0 {
        let cin = circ.alloc_qubit(); // carry_0 = 0
        let a_full: Vec<QubitId> = reg[..n].to_vec();
        let dirty: Vec<QubitId> = (lsbs..lsbs + (n - 1)).map(|i| reg[i]).collect();
        dirty_carryin(circ, ctrl, &a_full, c, 0, &dirty, &cin);
        circ.zero_and_free(cin);
    } else {
        add_f_window_hybrid(circ, ctrl, reg, lsbs, c, g);
    }
}

fn add_f_window_clean(circ: &mut B, ctrl: &QubitId, reg: &[QubitId], lsbs: usize, c: &[u8]) {
    let n = lsbs;
    assert!(n <= reg.len(), "register too short for +f window");
    if n == 0 {
        return;
    }
    if n == 1 {
        if cbit(c, 0) {
            circ.cx(*ctrl, reg[0]);
        }
        return;
    }
    let a: Vec<QubitId> = reg[..n].to_vec();
    // carries[i] = cy_{i+1} (carry OUT of bit i); cy_0 = 0 elided. Need n-1.
    let mut cy: Vec<Option<QubitId>> = (0..n - 1).map(|_| Some(circ.alloc_qubit())).collect();

    // Forward bit 0: cy_1 = ctrl & c_0 & a_0.
    if cbit(c, 0) {
        circ.ccx(*ctrl, a[0], *cy[0].as_ref().unwrap());
    }
    // Forward bits 1..n-1: carry_{i+1} = cy_i ^ ((a_i ^ cy_i) & (eff_i ^ cy_i)).
    for i in 1..n - 1 {
        let ci = cy[i - 1].take().unwrap();
        let next = cy[i].take().unwrap();
        circ.cx(ci, a[i]); // a[i] = ta = a_i ^ cy_i
        if cbit(c, i) {
            circ.cx(*ctrl, ci); // ci = cy_i ^ ctrl = tb
        }
        circ.ccx(a[i], ci, next); // next = ta & tb
        if cbit(c, i) {
            circ.cx(*ctrl, ci); // ci -> cy_i
        }
        circ.cx(ci, next); // next = carry_{i+1}
        cy[i - 1] = Some(ci);
        cy[i] = Some(next);
    }

    // Sums: a[i] ^= eff_i for i<n-1; top bit a[n-1] ^= eff_{n-1} ^ cy_{n-1}.
    for i in 0..n - 1 {
        if cbit(c, i) {
            circ.cx(*ctrl, a[i]);
        }
    }
    if cbit(c, n - 1) {
        circ.cx(*ctrl, a[n - 1]);
    }
    circ.cx(*cy[n - 2].as_ref().unwrap(), a[n - 1]); // cy_{n-1} into top bit (dropped beyond)

    // Reverse bits n-2..1: undo sum -> ta, vent the carry AND, restore.
    for i in (1..n - 1).rev() {
        if cbit(c, i) {
            circ.cx(*ctrl, a[i]); // a[i] = ta
        }
        let next = cy[i].take().unwrap(); // cy_{i+1}
        let ci = cy[i - 1].take().unwrap(); // cy_i
        circ.cx(ci, next); // next = ta & tb (undo carry XOR)
        if cbit(c, i) {
            circ.cx(*ctrl, ci); // ci = tb
        }
        // Vent carry AND: hmr(next), then CZ(ta=a[i], tb=ci) gated on the hmr
        // bit cancels the deferred phase (carry == a[i] & ci on the vented shots).
        let mbit = circ.alloc_bit();
        circ.hmr(next, mbit);
        circ.zero_and_free(next);
        circ.cz_if_bit(a[i], ci, mbit);
        if cbit(c, i) {
            circ.cx(*ctrl, ci); // ci -> cy_i
            circ.cx(*ctrl, a[i]); // ta -> sum_i
        }
        cy[i - 1] = Some(ci);
    }
    // Reverse bit 0: vent cy_1 = a_0 & (ctrl & c_0).
    let cy1 = cy[0].take().unwrap();
    if cbit(c, 0) {
        circ.cx(*ctrl, a[0]); // sum_0 -> ta_0 = a_0
        let mbit = circ.alloc_bit();
        circ.hmr(cy1, mbit);
        circ.zero_and_free(cy1);
        // carry_1 = a_0 & ctrl ; CZ(a[0], ctrl) gated on the hmr bit.
        circ.cz_if_bit(a[0], *ctrl, mbit);
        circ.cx(*ctrl, a[0]); // ta_0 -> sum_0
    } else {
        // c_0 = 0 -> cy_1 never set, still |0>.
        circ.zero_and_free(cy1);
    }
}


/// `reg[..lsbs] -= ctrl * c` (X-sandwich of [`add_f_window`]).
fn sub_f_window(circ: &mut B, ctrl: &QubitId, reg: &[QubitId], lsbs: usize, c: &[u8]) {
    for q in &reg[..lsbs] {
        circ.x(*q);
    }
    add_f_window(circ, ctrl, reg, lsbs, c, None);
    for q in &reg[..lsbs] {
        circ.x(*q);
    }
}

// =====================================================================
// Truncated top-k less-than overflow cleanup (normal)
// =====================================================================

/// Conditional (measurement-vented) clear of a flag that already holds
/// `ctrl AND (a_top < b_top)` on entry (the post-add+fold overflow anc). HMR the
/// flag (0 Toffoli, resets it to |0>), then on the ~50% of shots where the measured
/// bit fired, recompute the truncated top-`k` predicate and deposit it as a Z-phase
/// correction gated on the HMR bit (a `push_condition` on the HMR bit + a
/// `cz_if_bit` deposit). Net: flag -> |0>, phase-clean, with the comparator running
/// on only `k` bits and only ~half the shots. The HMR (vs the normal
/// `zero_and_free`) also tolerates the rare top-`k` mis-clear (it never asserts
/// |0>), which is what lets the comparator be truncated to `k = MSBS` (~2^-MSBS
/// per-call mis-decide -- the schedule-sim's CMP/PAD term accounts for it).
///
/// PRECONDITION: `target` holds `ctrl AND (a_top < b_top)` on entry. NOT valid for
/// a flag-recreating reverse (which must use the normal `controlled_lt_msbs`).
fn controlled_lt_msbs_conditional(circ: &mut B, ctrl: Option<&QubitId>, a: &[QubitId], b: &[QubitId], k: usize, target: QubitId) {
    let a_top: Vec<QubitId> = a[a.len() - k..].to_vec();
    let b_top: Vec<QubitId> = b[b.len() - k..].to_vec();
    let bit = circ.alloc_bit();
    circ.hmr(target, bit); // measure the vented flag p = ctrl AND (a<b); reset to |0>.
    // Free the flag now -- before the recompute comparator -- so it is not held live
    // across the (chunked, window-carry) comparator. The predicate is recomputed from
    // the live operands a,b, not from `target`, so this is the same free-then-recompute
    // ordering the mod-sub borrow clean uses (`mod_sub_vented`). Holding it through the
    // comparator would cost +1 qubit at the forward-multiply cofactor-add peak.
    circ.zero_and_free(target);
    let ctrl = ctrl.copied();
    circ.push_condition(bit);
    // On the gated (HMR-bit) shots, recompute the predicate as a deferred Z through the
    // headroom-adaptive (chunked) comparator backend
    // (`compare_geq_chunked_middle`): a flag-based `[a_top >= b_top]`, with
    // the held-carry count = the full window `k` (effk == k under the ample +f-fold
    // headroom). Deposit Z^(ctrl AND NOT flag) = Z^(ctrl AND (a < b)), gated by the
    // active HMR condition.
    let lt_flag = circ.alloc_qubit();
    super::comparator::compare_geq_chunked_middle(
        circ,
        &a_top,
        &b_top,
        &lt_flag,
        |c, flag| {
            c.x(*flag);
            match &ctrl {
                Some(ct) => c.cz(*ct, *flag),
                None => c.z(*flag),
            }
            c.x(*flag);
        },
        k,
    );
    circ.zero_and_free(lt_flag);
    circ.pop_condition();
}

/// CONDITIONAL (measurement-vented) clear of a flag holding
/// `ctrl AND carryout(a_top + b_top)` on entry (the mod-SUB borrow anc). HMR the
/// flag, then on the gated shots recompute the top-`k` add-carry predicate as a
/// Z-phase. Same truncation/tolerance rationale as
/// [`controlled_lt_msbs_conditional`].
///
/// Identity (as in the normal form): carryout(a_top + b_top) over `k` bits =
/// `NOT(~b_top >= a_top)`. Flip `b_top -> ~b_top`, recompute `flag = (~b_top >=
/// a_top)`, deposit `Z^(ctrl AND NOT flag)`, un-flip.
///
/// PRECONDITION: `target` holds `ctrl AND carryout(...)` on entry.
fn controlled_add_carry_msbs_conditional(circ: &mut B, ctrl: Option<&QubitId>, a: &[QubitId], b: &[QubitId], k: usize, target: &QubitId) {
    let a_top: Vec<QubitId> = a[a.len() - k..].to_vec();
    let b_top: Vec<QubitId> = b[b.len() - k..].to_vec();
    let bit = circ.alloc_bit();
    circ.hmr(*target, bit);
    circ.push_condition(bit);
    for q in &b_top {
        circ.x(*q); // b_top -> ~b_top
    }
    // Vented chunked comparator backend (`compare_geq_chunked_middle`); held-carry
    // count = the full window `k`. Predicate unchanged: flag = (~b_top >= a_top),
    // carry = NOT flag, deposit Z^(ctrl AND carry). UNCONTROLLED (`ctrl = None`, the
    // unconditional coord/square subs): bare `Z^carry`.
    let ctrl = ctrl.copied();
    let lt_flag = circ.alloc_qubit();
    super::comparator::compare_geq_chunked_middle(circ, &b_top, &a_top, &lt_flag, |c, flag| {
        c.x(*flag);
        match &ctrl {
            Some(ct) => c.cz(*ct, *flag),
            None => c.z(*flag),
        }
        c.x(*flag);
    }, k);
    circ.zero_and_free(lt_flag);
    for q in &b_top {
        circ.x(*q); // restore b_top
    }
    circ.pop_condition();
}

// =====================================================================
// Public primitives
// =====================================================================

/// Controlled mod-add with an optional SCHEDULE carry-cap `k`. `Some(k)` (the
/// apply cofactor add) emits the adaptive cout decomposition at the baked `k`.
/// `None` (the off-peak coordinate steps) uses the local exact cout adder.
pub fn controlled_mod_add_k(circ: &mut B, ctrl: &QubitId, x: &[QubitId], y: &[QubitId], sched_k: Option<usize>, ffg_g: Option<usize>) {
    let n = x.len();
    assert_eq!(y.len(), n, "x,y must both be n=256 bits");
    assert_eq!(n, 256, "secp256k1 controlled_mod_add expects n=256");
    let f_bytes = F_SECP256K1.to_le_bytes();
    let anc = circ.alloc_qubit();
    // 1) y += ctrl*x; carry-out (overflow) captured into the transient anc.
    match sched_k {
        Some(k) => {
            let yr: Vec<&QubitId> = y.iter().collect();
            let xr: Vec<&QubitId> = x.iter().collect();
            super::gidney::controlled_hybrid_add_cout_refs(circ, ctrl, &yr, &xr, &anc, k);
        }
        None => controlled_add_vented_chunked_cout(circ, ctrl, x, y, APPLY_CHUNK, Some(&anc)),
    }
    // 2) gated +f reduction (anc holds ctrl AND overflow); carry beyond LSBS dropped.
    add_f_window(circ, &anc, y, LSBS, &f_bytes, ffg_g);
    // 3) less-than comparator erases anc back to |0>: anc holds `ctrl AND (y_final < x)`. The
    //    ludicrous profile truncates this comparator to the top `MSBS = PAD = 21`
    //    bits and clears the flag by measurement-vent (HMR + gated top-k Z), not a
    //    normal full-width comparator + `zero_and_free`. The HMR never asserts
    //    |0>, so the ~2^-MSBS top-k mis-clear is tolerated (the schedule-sim's
    //    CMP/PAD term accounts for it). ~21 vs ~256 CCX per call (the full-width
    //    k=256 normal form), on ~half the shots.
    debug_assert_eq!(MSBS, PAD); // the +f comparator is truncated to the top PAD bits
    // The less-than erase consumes `anc`: it HMRs it then frees it BEFORE the recompute
    // comparator (so the overflow flag is not held live across the chunked comparator).
    controlled_lt_msbs_conditional(circ, Some(ctrl), &y[..n], &x[..n], MSBS, anc);
}

/// In-place pseudo-Mersenne modular subtraction `y := y - x (mod q)` over
/// n=256-bit `x`,`y`. The register sub is normal UNCONTROLLED Cuccaro
/// (`cuccaro_carry(None)` -- `cx` sums, NO |1>-gated `ccx`): X-sandwich +
/// Cuccaro add + gated `-f` fold + top-MSBS add-carry borrow clean. The
/// borrow-clean deposit is `Z^carry` -- reused via the |1>-ctrl form (the
/// `cz(|1>,flag)` degenerates to a free `Z`; the comparator CCX is
/// ctrl-independent), so no register-sub CCX waste.
pub fn mod_sub(circ: &mut B, x: &[QubitId], y: &[QubitId]) {
    let n = x.len();
    assert_eq!(y.len(), n, "x,y must both be n=256 bits");
    assert_eq!(n, 256, "secp256k1 mod_sub expects n=256");
    let f_bytes = F_SECP256K1.to_le_bytes();
    let anc = circ.alloc_qubit();
    // X-sandwich: ~y += x => y -= x; cout = borrow (UNCONTROLLED, cx sums).
    for q in y {
        circ.x(*q);
    }
    cuccaro_carry(circ, None, x, y, None, Some(&anc));
    for q in y {
        circ.x(*q);
    }
    // gated -f fold on the borrow.
    sub_f_window(circ, &anc, y, LSBS, &f_bytes);
    // clean anc: top-MSBS add-carry predicate (UNCONTROLLED, bare Z^carry).
    controlled_add_carry_msbs_conditional(circ, None, &y[..n], &x[..n], MSBS, &anc);
    circ.zero_and_free(anc);
}

/// UNCONTROLLED VENTED add `y += x (mod 2^n)`, carry-out of bit `n-1` -> `cout`
/// (caller-owned |0>). Full-clean path: zero-pad to `n+1` and run the uncontrolled
/// vented `hybrid_add_plain` at full-clean headroom (vents = n). ~1 Toffoli/bit, NO
/// |1>-gated `eff` (the coord adds are unconditional). Temp-clean regardless of
/// input: `zpad` is algebraically restored, the vent ancillae HMR-reset.
fn add_cout_vented_unctrl(circ: &mut B, x: &[QubitId], y: &[QubitId], cout: &QubitId) {
    let n = y.len();
    assert_eq!(x.len(), n, "add_cout_vented_unctrl: x,y width mismatch");
    let zpad = circ.alloc_qubit();
    let mut a: Vec<QubitId> = y.to_vec();
    a.push(*cout);
    let mut b: Vec<QubitId> = x.to_vec();
    b.push(zpad);
    hybrid_add_plain(circ, &a, &b, n); // vents = n => full clean over the (n+1)-bit add
    circ.zero_and_free(zpad);
}

/// UNCONTROLLED `y := y + x (mod q)` over n=256-bit x,y. VENTED register add
/// (overflow -> anc) + CLEAN gated `+f` fold + top-MSBS `(y<x)` overflow clean.
/// Temp-clean (anc consumed by the clt). The overflow clean deposits a bare `Z`
/// (UNCONTROLLED -- no |1> qubit).
pub fn mod_add(circ: &mut B, x: &[QubitId], y: &[QubitId]) {
    let n = x.len();
    assert_eq!(y.len(), n, "mod_add: x,y must both be n=256 bits");
    assert_eq!(n, 256, "secp256k1 mod_add expects n=256");
    let f_bytes = F_SECP256K1.to_le_bytes();
    let anc = circ.alloc_qubit();
    add_cout_vented_unctrl(circ, x, y, &anc);
    // gated +f fold on the overflow (clean carries, windowed -- same tail as the
    // controlled coord add, already modeled in the schedule sim).
    add_f_window(circ, &anc, y, LSBS, &f_bytes, Some(LSBS - 1));
    // clean anc: anc ^= (y_top < x_top) over the top MSBS bits (consumes anc).
    controlled_lt_msbs_conditional(circ, None, &y[..n], &x[..n], MSBS, anc);
}

/// UNCONTROLLED VENTED `y := y - x (mod q)`. X-sandwich + VENTED register sub
/// (`add_cout_vented_unctrl`) + CLEAN gated `-f` fold + top-MSBS borrow clean.
/// Distinct from [`mod_sub`] (which uses a normal Cuccaro register sub).
/// Temp-clean (anc HMR-reset before free).
pub fn mod_sub_vented(circ: &mut B, x: &[QubitId], y: &[QubitId]) {
    let n = x.len();
    assert_eq!(y.len(), n, "mod_sub_vented: x,y must both be n=256 bits");
    assert_eq!(n, 256, "secp256k1 mod_sub_vented expects n=256");
    let f_bytes = F_SECP256K1.to_le_bytes();
    let anc = circ.alloc_qubit();
    for q in y {
        circ.x(*q);
    }
    add_cout_vented_unctrl(circ, x, y, &anc);
    for q in y {
        circ.x(*q);
    }
    // gated -f fold on the borrow = X-sandwich of a forced-CLEAN +f window.
    for q in &y[..LSBS] {
        circ.x(*q);
    }
    add_f_window(circ, &anc, y, LSBS, &f_bytes, Some(LSBS - 1));
    for q in &y[..LSBS] {
        circ.x(*q);
    }
    controlled_add_carry_msbs_conditional(circ, None, &y[..n], &x[..n], MSBS, &anc);
    circ.zero_and_free(anc);
}

/// In-place modular negate `x := q - x (mod q)` for x in (0,q). Identity (since
/// `q = 2^256 - f`): `q - x = ~(x + (f-1))`. So one exact full-width const-add of
/// `(f-1)` (no carry escapes 2^256 since `x + (f-1) < q + f - 1 = 2^256 - 1`) then
/// flip all 256 bits. Folding the `+1` of a flip+inc+sub_const(f) form into the
/// constant drops the increment and the 257th carry bit. The full-width
/// `add_f_window_clean` is exact (no fold tail), so this adds nothing to the
/// schedule sim. Boundary: x=0 -> q (out of range); the EC tail only negates a
/// generic post-subtract value.
pub fn mod_neg(circ: &mut B, x: &[QubitId]) {
    let n = x.len();
    assert_eq!(n, 256, "secp256k1 mod_neg expects n=256");
    let f_minus_1 = (F_SECP256K1 - 1).to_le_bytes();
    add_const_window_clean(circ, x, n, &f_minus_1); // x += (f-1), full-width exact
    for q in x {
        circ.x(*q); // flip all -> ~(x + f - 1) = q - x
    }
}

/// UNCONDITIONAL clean const-add `reg[..lsbs] += c (mod 2^lsbs)` (carry beyond
/// `lsbs-1` dropped; at `lsbs = reg.len()` it is exact). Specialization of
/// [`add_f_window_clean`] with the control hardwired to |1>: gated `ccx(|1>,·)` ->
/// `cx`, `cx(|1>,·)` -> `x`, the bit-0 vent `cz_if_bit(·,|1>,·)` -> `z_if_bit`.
/// NO constant-control qubit. Measurement-vented carries (~lsbs-1 Toffoli).
fn add_const_window_clean(circ: &mut B, reg: &[QubitId], lsbs: usize, c: &[u8]) {
    let n = lsbs;
    assert!(n <= reg.len(), "register too short for const window");
    if n == 0 {
        return;
    }
    if n == 1 {
        if cbit(c, 0) {
            circ.x(reg[0]);
        }
        return;
    }
    let a: Vec<QubitId> = reg[..n].to_vec();
    let mut cy: Vec<Option<QubitId>> = (0..n - 1).map(|_| Some(circ.alloc_qubit())).collect();
    // Forward bit 0: cy_1 = a_0 & c_0.
    if cbit(c, 0) {
        circ.cx(a[0], *cy[0].as_ref().unwrap());
    }
    // Forward bits 1..n-1: carry_{i+1} = cy_i ^ ((a_i ^ cy_i) & (c_i ^ cy_i)).
    for i in 1..n - 1 {
        let ci = cy[i - 1].take().unwrap();
        let next = cy[i].take().unwrap();
        circ.cx(ci, a[i]); // a[i] = ta = a_i ^ cy_i
        if cbit(c, i) {
            circ.x(ci); // ci = cy_i ^ c_i = tb
        }
        circ.ccx(a[i], ci, next); // next = ta & tb
        if cbit(c, i) {
            circ.x(ci); // ci -> cy_i
        }
        circ.cx(ci, next); // next = carry_{i+1}
        cy[i - 1] = Some(ci);
        cy[i] = Some(next);
    }
    // Sums: a[i] ^= c_i for i<n-1; top bit also gets cy_{n-1}.
    for i in 0..n - 1 {
        if cbit(c, i) {
            circ.x(a[i]);
        }
    }
    if cbit(c, n - 1) {
        circ.x(a[n - 1]);
    }
    circ.cx(*cy[n - 2].as_ref().unwrap(), a[n - 1]);
    // Reverse bits n-2..1: undo sum -> ta, measurement-vent the carry AND, restore.
    for i in (1..n - 1).rev() {
        if cbit(c, i) {
            circ.x(a[i]); // a[i] = ta
        }
        let next = cy[i].take().unwrap();
        let ci = cy[i - 1].take().unwrap();
        circ.cx(ci, next); // next = ta & tb (undo carry XOR)
        if cbit(c, i) {
            circ.x(ci); // ci = tb
        }
        let mbit = circ.alloc_bit();
        circ.hmr(next, mbit);
        circ.zero_and_free(next);
        circ.cz_if_bit(a[i], ci, mbit); // CZ(ta, tb) cancels the deferred carry phase
        if cbit(c, i) {
            circ.x(ci); // ci -> cy_i
            circ.x(a[i]); // ta -> sum_i
        }
        cy[i - 1] = Some(ci);
    }
    // Reverse bit 0: vent cy_1 = a_0 & c_0.
    let cy1 = cy[0].take().unwrap();
    if cbit(c, 0) {
        circ.x(a[0]); // sum_0 -> ta_0 = a_0
        let mbit = circ.alloc_bit();
        circ.hmr(cy1, mbit);
        circ.zero_and_free(cy1);
        circ.z_if_bit(a[0], mbit); // carry_1 = a_0 (c_0 = 1); Z^(a_0) gated on the hmr bit
        circ.x(a[0]); // ta_0 -> sum_0
    } else {
        circ.zero_and_free(cy1);
    }
}

/// In-place pseudo-Mersenne modular doubling: `a := 2*a (mod q)` (Alg 7).
/// `a` is `n+1 = 257` bits: `a[0..n]` holds x in [0,q); `a[n] = |0>` (overflow
/// slot, restored).
///
/// normal shift + the (MBU) `+f` fold + a normal CX ancilla cleanup.
pub fn mod_double(circ: &mut B, a: &[QubitId]) {
    let n = a.len() - 1;
    assert_eq!(n, 256, "secp256k1 mod_double expects 257-bit a");
    let f_bytes = F_SECP256K1.to_le_bytes();
    // 1) a := 2*a: value bit i -> slot i+1; old MSB lands in a[n] (overflow).
    //    A reversible left-shift by 1 = swap chain from the top down.
    for i in (0..n).rev() {
        circ.swap(a[i], a[i + 1]);
    }
    // 2) if a[n] (overflow), add f to the bottom LSBS bits (carry beyond dropped).
    add_f_window(circ, &a[n], a, LSBS, &f_bytes, None);
    // 3) clean a[n]: after step 1, a[0] = 0; after step 2 with f odd (bit 0 = 1),
    //    a[0] == a[n]. So CX(a[0], a[n]) clears a[n].
    circ.cx(a[0], a[n]);
}

/// Inverse of [`mod_double`]: `a := a / 2 (mod q)` (maps `2x mod q -> x`).
/// Gate-for-gate reverse of [`mod_double`] -- undo the CX cleanup, X-sandwich the
/// `+f` window into a `-f`, then reverse the left-shift. `a` is `n+1 = 257` bits
/// (a[256] = |0> overflow slot, restored). Used to restore an operand after a
/// doubling ramp (e.g. the square's `f*hi` reduce).
pub fn mod_double_reverse(circ: &mut B, a: &[QubitId]) {
    let n = a.len() - 1;
    assert_eq!(n, 256, "secp256k1 mod_double_reverse expects 257-bit a");
    let f_bytes = F_SECP256K1.to_le_bytes();
    // 3') undo the CX(a[0], a[n]) cleanup.
    circ.cx(a[0], a[n]);
    // 2') -f fold gated on a[n] (X-sandwich the +f window).
    sub_f_window(circ, &a[n], a, LSBS, &f_bytes);
    // 1') reverse the left-shift (swap chain, bottom-up).
    for i in 0..n {
        circ.swap(a[i], a[i + 1]);
    }
}

/// Public re-export of [`add_f_window`] for the gcd apply's controlled doubling
/// (which needs the gated `+f` fold on an arbitrary control). Same contract.
pub fn add_f_window_pub(circ: &mut B, ctrl: &QubitId, reg: &[QubitId], lsbs: usize, c: &[u8], g_sched: Option<usize>) {
    add_f_window(circ, ctrl, reg, lsbs, c, g_sched);
}

