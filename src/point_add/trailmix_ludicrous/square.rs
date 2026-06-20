//! Modular square-subtract `output -= lambda^2 mod q` (secp256k1) used by the
//! EC point-add, built on the sibling `super::arith` mod-sub / mod-double
//! primitives.
//!
//! - [`symmetric_square_into_prod`]: now uses 1-level Karatsuba square for n=256
//!   (3 symmetric schoolbook subs on h=128 + combines) to cut cross-product
//!   CCX from n(n-1)/2=32640 down to ~3*h(h-1)/2 + combine overhead. Falls
//!   back to schoolbook for other sizes. Returns aux state for reverse.
//!   (See Karatsuba warning on qubit/Toffoli trade-off below.)
//! - The old pure schoolbook is retained as `schoolbook_*` helpers.
//! - [`mod_square_sub_pm_secp256k1_symmetric`]: the unconditional Stage-2 reduce
//!   `output -= lo + f*hi mod q`, built from `super::arith::{mod_double,
//!   mod_sub}`.
//!
//! ## secp256k1 constants
//!   q   = 2^256 - f,   f = 2^32 + 977   (bits {0,4,6,7,8,9,32})
//!   PAD = 21  (the +f window carry-drop -> ~2^-PAD per-fire approximation,
//!              inherited from `super::arith`'s mod-sub / mod-double folds).

use super::arith::{mod_double, mod_double_reverse, mod_sub};
use super::{B, BExt};
use crate::circuit::{QubitId};

const N: usize = 256;

/// Toffoli-free AND-uncompute (HMR + conditional-Z): `t` holds `a AND b` (here a
/// square cross-product `x_i AND x_j`); the HMR measures it to |0> and the
/// `cz_if_bit` cancels the deferred phase. Replaces the explicit reverse `ccx`
/// (1 Toffoli) with a measurement (0 Toffoli).
fn clear_and(circ: &mut B, t: &QubitId, a: &QubitId, b: &QubitId) {
    let bit = circ.alloc_bit();
    circ.hmr(*t, bit);
    circ.cz_if_bit(*a, *b, bit);
}

/// Set bits of f = 2^32 + 977 = bits {0,4,6,7,8,9,32}. `lambda^2 == lo + f*hi`,
/// so `f*hi = sum_{j in F_BITS} hi*2^j` -- one mod-double ramp + mod-sub per bit.
const F_BITS: [usize; 7] = [0, 4, 6, 7, 8, 9, 32];

/// `slice += row` (mod 2^slice.len) via `arith::hybrid_add_adaptive`. `slice` is
/// exactly one bit wider than `row` (one carry slot); the row carry rides into that top
/// slot (or, when this slice is an interior window of a wider accumulator, into
/// the already-populated high bits of `prod` -- the caller sizes the slice so the
/// final carry lands in a real |0> or populated slot, never dropped).
///
/// One clean zero-pad qubit, freed.
fn add_into(circ: &mut B, slice: &[QubitId], row: &[QubitId]) {
    let m = row.len();
    assert_eq!(slice.len(), m + 1, "slice must be one wider than row");
    if m == 0 {
        return;
    }
    // Zero-pad `row` to the slice width and run the UNCONTROLLED exact adaptive add
    // `slice += row_padded` (mod 2^(m+1)); the row carry rides into slice[m] (the
    // pad keeps the addend's top bit |0>). The adder's headroom `k` is the value
    // baked into the row-add schedule (SQ_ROW_K), read via next_sqrow_k().
    let pad = circ.alloc_qubit();
    let mut b: Vec<QubitId> = row.to_vec();
    b.push(pad);
    let k = super::next_sqrow_k();
    super::arith::hybrid_add_adaptive(circ, slice, &b, k);
    circ.zero_and_free(pad);
}

/// Karatsuba half-sum helpers (for z1 = (lo + hi)^2).
/// lo + hi (with possible carry out into extra bit).
fn karatsuba_half_sum_compute(circ: &mut B, lo: &[QubitId], hi: &[QubitId], sum: &[QubitId]) {
    let n = lo.len();
    assert_eq!(hi.len(), n);
    assert_eq!(sum.len(), n + 1);
    for i in 0..n {
        circ.cx(lo[i], sum[i]);
        circ.cx(hi[i], sum[i]);
    }
    // simple carry into the extra bit (can be improved with hybrid later)
    let mut carry = circ.alloc_qubit();
    circ.ccx(lo[0], hi[0], carry); // placeholder for proper carry chain
    circ.cx(carry, sum[n]);
    circ.zero_and_free(carry);
}

fn karatsuba_half_sum_uncompute(circ: &mut B, lo: &[QubitId], hi: &[QubitId], sum: &[QubitId]) {
    let n = lo.len();
    let mut carry = circ.alloc_qubit();
    circ.cx(carry, sum[n]);
    circ.ccx(lo[0], hi[0], carry);
    for i in 0..n {
        circ.cx(hi[i], sum[i]);
        circ.cx(lo[i], sum[i]);
    }
    circ.zero_and_free(carry);
}

// ─────────────────────────────────────────────────────────────────────────────
// Karatsuba helpers for symmetric square (1-level, n=256 power-of-2 split).
// Schoolbook cross-products: n(n-1)/2  →  3*(h(h-1)/2) with h=n/2.
// Cross CCX saving ~ (n(n-1)/2 - 3 h (h-1)/2) = 8256 for n=256.
// Trade-off: +~258 qubits peak during reduction (z1_reg lives across Stage 2
// to enable exact reverse-combine + z1 uncompute). Net T win, Q penalty.
// The combine steps (z1 -=z0, z1-=z2, mid+=z1) use hybrid_add (consume sqrow
// schedule or fallback); for first version we accept schedule overread (MAX
// vents for extra) and some cursor shift for sub-row k's.
// ─────────────────────────────────────────────────────────────────────────────

/// Uncontrolled half-sum `acc[0..h+1] := lo + hi` (acc initially |0|, sized h+1).
/// acc high bit receives the carry. Uses CX for lo + hybrid for hi (matches
/// karatsuba in arith but using trailmix hybrid + next_sqrow_k for vent budget).
fn karatsuba_half_sum_compute(circ: &mut B, lo: &[QubitId], hi: &[QubitId], acc: &[QubitId]) {
    let h = lo.len();
    debug_assert_eq!(hi.len(), h);
    debug_assert_eq!(acc.len(), h + 1);
    for i in 0..h {
        circ.cx(lo[i], acc[i]);
    }
    let hi_pad = circ.alloc_qubit();
    let mut hi_ext = hi.to_vec();
    hi_ext.push(hi_pad);
    let k = super::next_sqrow_k();
    super::arith::hybrid_add_adaptive(circ, acc, &hi_ext, k);
    circ.zero_and_free(hi_pad);
}

/// Gate-reverse of half-sum: acc -= hi then CX-restore lo. Uses X-sandwich +
/// hybrid_add (mirrors how reverse square subtracts rows).
fn karatsuba_half_sum_uncompute(circ: &mut B, lo: &[QubitId], hi: &[QubitId], acc: &[QubitId]) {
    let h = lo.len();
    let hi_pad = circ.alloc_qubit();
    let mut hi_ext = hi.to_vec();
    hi_ext.push(hi_pad);
    // acc -= hi_ext via X-sandwich + add (see symmetric reverse for the trick)
    for q in acc {
        circ.x(*q);
    }
    let k = super::next_sqrow_k();
    super::arith::hybrid_add_adaptive(circ, acc, &hi_ext, k);
    for q in acc {
        circ.x(*q);
    }
    circ.zero_and_free(hi_pad);
    for i in 0..h {
        circ.cx(lo[i], acc[i]);
    }
}

/// Symmetric schoolbook square of x[m] written directly into pre-allocated
/// `target` (exactly 2m |0> qubits). Mirrors the row logic of the original
/// but targets fixed slices (for z0/z2 placement inside a big prod).
fn schoolbook_square_symmetric_into(circ: &mut B, x: &[QubitId], target: &[QubitId]) {
    let m = x.len();
    debug_assert_eq!(target.len(), 2 * m);
    for i in 0..m {
        let num_cross = m.saturating_sub(i + 1);
        let width = if i == m - 1 { 1 } else { m - i + 1 };
        let row: Vec<QubitId> = (0..width).map(|_| circ.alloc_qubit()).collect();
        circ.cx(x[i], row[0]);
        for k in 0..num_cross {
            circ.ccx(x[i], x[i + 1 + k], row[k + 2]);
        }
        let end = (2 * i + width + 1).min(2 * m);
        add_into(circ, &target[2 * i..end], &row);
        for k in 0..num_cross {
            clear_and(circ, &row[k + 2], &x[i], &x[i + 1 + k]);
        }
        circ.cx(x[i], row[0]);
        for q in row {
            circ.zero_and_free(q);
        }
    }
}

/// Reverse of [`schoolbook_square_symmetric_into`]: rebuilds rows and subtracts
/// (via X-sandwich add) from target slices; uncomputes crosses. Target bits
/// are driven back to |0>.
fn schoolbook_square_symmetric_into_reverse(circ: &mut B, x: &[QubitId], target: &[QubitId]) {
    let m = x.len();
    debug_assert_eq!(target.len(), 2 * m);
    let mut sim_len = 2 * m; // simulate original shrinking prod.len() during rev pops for matching windows
    for i in (0..m).rev() {
        let num_cross = m.saturating_sub(i + 1);
        let width = if i == m - 1 { 1 } else { m - i + 1 };
        let hi = (2 * i + width + 1).min(sim_len);
        let row: Vec<QubitId> = (0..width).map(|_| circ.alloc_qubit()).collect();
        circ.cx(x[i], row[0]);
        for k in 0..num_cross {
            circ.ccx(x[i], x[i + 1 + k], row[k + 2]);
        }
        let sl = &target[2 * i..hi];
        for q in sl {
            circ.x(*q);
        }
        add_into(circ, sl, &row);
        for q in sl {
            circ.x(*q);
        }
        for k in 0..num_cross {
            clear_and(circ, &row[k + 2], &x[i], &x[i + 1 + k]);
        }
        circ.cx(x[i], row[0]);
        for q in row {
            circ.zero_and_free(q);
        }
        // simulate keep shrink (mirror original)
        let keep = (m + i + 1).min(2 * m);
        if sim_len > keep {
            sim_len = keep;
        }
    }
}

/// Schoolbook wrapper -- exact original lazy logic (to preserve proven behavior for non-kara and z1).
fn schoolbook_symmetric_square_into_prod(circ: &mut B, x: &[QubitId], prod: &mut Vec<QubitId>) {
    let n = x.len();
    assert!(prod.is_empty(), "prod is grown lazily; pass an empty Vec");
    for i in 0..n {
        let num_cross = n.saturating_sub(i + 1);
        let width = if i == n - 1 { 1 } else { n - i + 1 };
        let hi = (2 * i + width + 1).min(2 * n);
        while prod.len() < hi {
            prod.push(circ.alloc_qubit());
        }
        let row: Vec<QubitId> = (0..width).map(|_| circ.alloc_qubit()).collect();
        circ.cx(x[i], row[0]);
        for k in 0..num_cross {
            circ.ccx(x[i], x[i + 1 + k], row[k + 2]);
        }
        add_into(circ, &prod[2 * i..hi], &row);
        for k in 0..num_cross {
            clear_and(circ, &row[k + 2], &x[i], &x[i + 1 + k]);
        }
        circ.cx(x[i], row[0]);
        for q in row {
            circ.zero_and_free(q);
        }
    }
    debug_assert_eq!(prod.len(), 2 * n, "prod must reach 2n after the build");
}

/// Schoolbook reverse wrapper -- exact original (with lazy free during).
fn schoolbook_symmetric_square_into_prod_reverse(circ: &mut B, x: &[QubitId], mut prod: Vec<QubitId>) {
    let n = x.len();
    assert_eq!(prod.len(), 2 * n);
    for i in (0..n).rev() {
        let num_cross = n.saturating_sub(i + 1);
        let width = if i == n - 1 { 1 } else { n - i + 1 };
        let row: Vec<QubitId> = (0..width).map(|_| circ.alloc_qubit()).collect();
        circ.cx(x[i], row[0]);
        for k in 0..num_cross {
            circ.ccx(x[i], x[i + 1 + k], row[k + 2]);
        }
        let hi = (2 * i + width + 1).min(prod.len());
        for q in &prod[2 * i..hi] {
            circ.x(*q);
        }
        add_into(circ, &prod[2 * i..hi], &row);
        for q in &prod[2 * i..hi] {
            circ.x(*q);
        }
        for k in 0..num_cross {
            clear_and(circ, &row[k + 2], &x[i], &x[i + 1 + k]);
        }
        circ.cx(x[i], row[0]);
        for q in row {
            circ.zero_and_free(q);
        }
        let keep = (n + i + 1).min(2 * n);
        while prod.len() > keep {
            circ.zero_and_free(prod.pop().unwrap());
        }
    }
    for q in prod {
        circ.zero_and_free(q);
    }
}

/// 1-level Karatsuba decomposition for the prod build (only for n==256).
/// z0=lo^2, z2=hi^2 written directly into their final positions in prod via
/// schoolbook subs; z1=(lo+hi)^2 computed early (before prod alloc) then
/// reduced to 2·lo·hi and folded at offset h. Returns the post-combine z1
/// (holding 2·lo·hi) so caller can keep it live across reduction for reverse.
fn karatsuba_symmetric_square_into_prod(circ: &mut B, x: &[QubitId], prod: &mut Vec<QubitId>) -> Vec<QubitId> {
    let n = x.len();
    debug_assert_eq!(n, 256);
    assert!(prod.is_empty());
    let h = n / 2;
    let x_lo = &x[0..h];
    let x_hi = &x[h..n];

    // z1 first (low peak): (lo+hi)^2 before allocating the 512-bit prod.
    let mut z1: Vec<QubitId> = Vec::new();
    {
        let x_sum: Vec<QubitId> = (0..=h).map(|_| circ.alloc_qubit()).collect();
        karatsuba_half_sum_compute(circ, x_lo, x_hi, &x_sum);
        schoolbook_symmetric_square_into_prod(circ, &x_sum, &mut z1);
        karatsuba_half_sum_uncompute(circ, x_lo, x_hi, &x_sum);
        for q in x_sum {
            circ.zero_and_free(q);
        }
    }
    debug_assert_eq!(z1.len(), 2 * (h + 1));

    // Pre-allocate full 2n prod (all |0>) now that z1 operand is freed.
    for _ in 0..(2 * n) {
        prod.push(circ.alloc_qubit());
    }

    // Build z0 and z2 in separate temps (using proven schoolbook), add their
    // values into the main prod at final locations, then use the temps for
    // z1 subs (pure values). Uncompute+free temps promptly.
    let mut z0p: Vec<QubitId> = Vec::new();
    schoolbook_symmetric_square_into_prod(circ, x_lo, &mut z0p);
    // add z0p into prod low (exact width 256, fits)
    {
        let low: Vec<QubitId> = prod[0..2 * h].to_vec();
        let k = 1000usize;
        super::arith::hybrid_add_adaptive(circ, &low, &z0p, k);
    }
    // z2
    let mut z2p: Vec<QubitId> = Vec::new();
    schoolbook_symmetric_square_into_prod(circ, x_hi, &mut z2p);
    {
        let hi: Vec<QubitId> = prod[2 * h..4 * h].to_vec();
        let k = 1000usize;
        super::arith::hybrid_add_adaptive(circ, &hi, &z2p, k);
    }

    // Combine on z1 using the pure z*p temps: z1 -= z0p; z1 -= z2p
    {
        let mut z0_ext: Vec<QubitId> = z0p.clone();
        let p0 = circ.alloc_qubit();
        let p1 = circ.alloc_qubit();
        z0_ext.push(p0);
        z0_ext.push(p1);
        for q in &z1 {
            circ.x(*q);
        }
        let k = 1000usize;
        super::arith::hybrid_add_adaptive(circ, &z1, &z0_ext, k);
        for q in &z1 {
            circ.x(*q);
        }
        circ.zero_and_free(p0);
        circ.zero_and_free(p1);
    }
    {
        let mut z2_ext: Vec<QubitId> = z2p.clone();
        let p0 = circ.alloc_qubit();
        let p1 = circ.alloc_qubit();
        z2_ext.push(p0);
        z2_ext.push(p1);
        for q in &z1 {
            circ.x(*q);
        }
        let k = 1000usize;
        super::arith::hybrid_add_adaptive(circ, &z1, &z2_ext, k);
        for q in &z1 {
            circ.x(*q);
        }
        circ.zero_and_free(p0);
        circ.zero_and_free(p1);
    }

    // mid add (from z1 reduced)
    {
        let mid_start = h;
        let mid_w = 3 * h; // 384
        let mut z1_ext: Vec<QubitId> = z1.clone();
        let np = mid_w.saturating_sub(z1_ext.len());
        let pads: Vec<QubitId> = (0..np).map(|_| circ.alloc_qubit()).collect();
        z1_ext.extend(pads);
        let mid_t: Vec<QubitId> = prod[mid_start..mid_start + mid_w].to_vec();
        let k = 1000usize;
        super::arith::hybrid_add_adaptive(circ, &mid_t, &z1_ext, k);
        for q in z1_ext.iter().skip(z1.len()) {
            circ.zero_and_free(*q);
        }
    }

    // Uncompute z0p and z2p (reverse their builds), freeing temps. Main prod
    // retains the placed values.
    schoolbook_symmetric_square_into_prod_reverse(circ, x_lo, z0p);
    schoolbook_symmetric_square_into_prod_reverse(circ, x_hi, z2p);

    // Leave z1 live (holds 2·lo·hi) for reduction + reverse. Caller frees later.
    z1
}

/// Reverse for the karatsuba prod build. Mirrors forward ordering exactly
/// (using the passed-in z1 that was kept live by caller).
fn karatsuba_symmetric_square_into_prod_reverse(circ: &mut B, x: &[QubitId], prod: Vec<QubitId>, z1: Vec<QubitId>) {
    let n = x.len();
    debug_assert_eq!(n, 256);
    debug_assert_eq!(prod.len(), 2 * n);
    debug_assert_eq!(z1.len(), 2 * (n / 2 + 1));
    let h = n / 2;
    let x_lo = &x[0..h];
    let x_hi = &x[h..n];

    // Undo mid add first: mid -= z1 (z1 currently 2lohi)
    {
        let mid_start = h;
        let mid_w = 3 * h;
        let mut z1_ext: Vec<QubitId> = z1.clone();
        let np = mid_w.saturating_sub(z1_ext.len());
        let pads: Vec<QubitId> = (0..np).map(|_| circ.alloc_qubit()).collect();
        z1_ext.extend(pads);
        let mid_t: Vec<QubitId> = prod[mid_start..mid_start + mid_w].to_vec();
        for q in &mid_t {
            circ.x(*q);
        }
        let k = 1000usize;
        super::arith::hybrid_add_adaptive(circ, &mid_t, &z1_ext, k);
        for q in &mid_t {
            circ.x(*q);
        }
        for q in z1_ext.iter().skip(z1.len()) {
            circ.zero_and_free(*q);
        }
    }

    // Rebuild pure z0p / z2p via school (from live x), use to restore z1 and
    // to undo the adds into main prod (subtract from low/high), then uncompute temps.
    let mut z0p: Vec<QubitId> = Vec::new();
    schoolbook_symmetric_square_into_prod(circ, x_lo, &mut z0p);
    let mut z2p: Vec<QubitId> = Vec::new();
    schoolbook_symmetric_square_into_prod(circ, x_hi, &mut z2p);

    // Restore z1 : += z0p ; += z2p
    {
        let mut z0_ext: Vec<QubitId> = z0p.clone();
        let p0 = circ.alloc_qubit();
        let p1 = circ.alloc_qubit();
        z0_ext.push(p0);
        z0_ext.push(p1);
        let k = 1000usize;
        super::arith::hybrid_add_adaptive(circ, &z1, &z0_ext, k);
        circ.zero_and_free(p0);
        circ.zero_and_free(p1);
    }
    {
        let mut z2_ext: Vec<QubitId> = z2p.clone();
        let p0 = circ.alloc_qubit();
        let p1 = circ.alloc_qubit();
        z2_ext.push(p0);
        z2_ext.push(p1);
        let k = 1000usize;
        super::arith::hybrid_add_adaptive(circ, &z1, &z2_ext, k);
        circ.zero_and_free(p0);
        circ.zero_and_free(p1);
    }

    // Undo the z adds to main prod (subtract recomputed z0p/z2p from slices)
    {
        let low: Vec<QubitId> = prod[0..2 * h].to_vec();
        let k = 1000usize;
        for q in &low {
            circ.x(*q);
        }
        super::arith::hybrid_add_adaptive(circ, &low, &z0p, k);
        for q in &low {
            circ.x(*q);
        }
    }
    {
        let hi: Vec<QubitId> = prod[2 * h..4 * h].to_vec();
        let k = 1000usize;
        for q in &hi {
            circ.x(*q);
        }
        super::arith::hybrid_add_adaptive(circ, &hi, &z2p, k);
        for q in &hi {
            circ.x(*q);
        }
    }

    // Uncompute the temps (this zeros z*p).
    schoolbook_symmetric_square_into_prod_reverse(circ, x_lo, z0p);
    schoolbook_symmetric_square_into_prod_reverse(circ, x_hi, z2p);

    // All prod bits now |0>; free them.
    for q in prod {
        circ.zero_and_free(q);
    }

    // Uncompute z1 last (rebuild operand, reverse the z1 square, un-sum).
    {
        let x_sum: Vec<QubitId> = (0..=h).map(|_| circ.alloc_qubit()).collect();
        karatsuba_half_sum_compute(circ, x_lo, x_hi, &x_sum);
        // pass z1 (the vec) to school reverse which drains/frees it
        schoolbook_symmetric_square_into_prod_reverse(circ, &x_sum, z1);
        karatsuba_half_sum_uncompute(circ, x_lo, x_hi, &x_sum);
        for q in x_sum {
            circ.zero_and_free(q);
        }
    }
}

/// Build `prod[0..2n] = value(x[0..n])^2` (integer, no reduction) via
/// Karatsuba (n=256) or fall back to symmetric schoolbook.
///
/// The public entry now returns auxiliary state (z1 for Karatsuba case) that
/// must be passed through to the matching reverse to enable exact uncompute.
fn symmetric_square_into_prod(circ: &mut B, x: &[QubitId], prod: &mut Vec<QubitId>) -> Vec<QubitId> {
    let n = x.len();
    assert!(prod.is_empty(), "prod is grown lazily; pass an empty Vec");
    if n == 256 && std::env::var("LUD_KARATSUBA_SQUARE").ok().as_deref() == Some("1") {
        karatsuba_symmetric_square_into_prod(circ, x, prod)
    } else {
        schoolbook_symmetric_square_into_prod(circ, x, prod);
        vec![]
    }
}

/// Gate-reverse of [`symmetric_square_into_prod`]: rebuilds each row and
/// SUBTRACTS it from `prod`, draining `prod` back to |0>. Now takes optional
/// auxiliary z1 state returned from the forward (non-empty only for Karatsuba
/// n=256 path). For schoolbook, pass empty vec.
fn symmetric_square_into_prod_reverse(circ: &mut B, x: &[QubitId], prod: Vec<QubitId>, z1: Vec<QubitId>) {
    let n = x.len();
    assert_eq!(prod.len(), 2 * n);
    if n == 256 && std::env::var("LUD_KARATSUBA_SQUARE").ok().as_deref() == Some("1") {
        karatsuba_symmetric_square_into_prod_reverse(circ, x, prod, z1);
    } else {
        assert!(z1.is_empty(), "z1 state only for Karatsuba path");
        schoolbook_symmetric_square_into_prod_reverse(circ, x, prod);
    }
}

/// Unconditional `output_reg -= lambda^2 mod q` (secp256k1), normal throughout.
///
/// `lambda` is `n = 256` bits (lambda < q); `output_reg` is `n = 256` bits and
/// holds a value < q on entry (the EC-add keeps output reduced).
///
/// Stage 1: build the 2n-bit integer product `prod = lambda^2`
/// with [`symmetric_square_into_prod`] (Karatsuba for n=256: cross CCX
/// ~24384 vs 32640, net win after combine costs).
///   WARNING: Karatsuba square increases peak qubit count by ~258 during
///   Stage 2 reduction (z1_reg kept live across mod_doubles/subs to enable
///   exact reverse). This trades Q for lower T/Toffoli. May impact score
///   (product-min is sensitive to 1166q floor) even if T drops.
/// Stage 2 (reduce): `lambda < q < 2^256 => lambda^2 < q^2 < 2^512`,
/// so `hi = prod>>256 < q`. With `2^256 == f (mod q)`, `lambda^2 == lo + f*hi`.
/// Subtract `lo` from `output`, then for each set bit j of f walk `hi` in place
/// by [`mod_double`] and subtract `hi*2^j mod q`; restore `hi` with the matched
/// reverse doublings. Uses `arith::mod_sub` (uncontrolled normal Cuccaro
/// register sub).
/// Stage 3: uncompute `prod` (gate-reverse of Stage 1).
///
/// Value note (carried-over miss probability): each `mod_double` / `mod_sub`
/// inherits `super::arith`'s `+f`-window carry drop -- a documented ~2^-PAD
/// (PAD=21) per-fire approximation. The common path is exact; the only legal
/// divergence is that rare large-input +f-window miss.
pub fn mod_square_sub_pm_secp256k1_symmetric(circ: &mut B, lambda: &[QubitId], output_reg: &[QubitId]) {
    let n = N;
    assert_eq!(lambda.len(), n, "lambda must be n=256 bits (< q)");
    assert_eq!(output_reg.len(), n, "output must be n=256 bits (< q)");

    // Stage 1: prod = lambda^2 (integer, 2n bits).
    // For n=256 the Karatsuba path returns auxiliary z1 state (lives across
    // reduction to support reverse combine/uncompute; see tradeoff warning).
    let mut prod: Vec<QubitId> = Vec::with_capacity(2 * n);
    let z1_state = symmetric_square_into_prod(circ, lambda, &mut prod);

    // Stage 2: output -= (lo + f*hi) mod q, operating on prod's own halves.
    //   lo = prod[0..n]                                  (n-bit, lo can be >= q)
    //   hi = prod[n..2n]                                 (n-bit, hi < q)
    // mod_double needs a 257-bit operand whose top bit is |0>; one inserted pad
    // above hi supplies that overflow slot (restored, removed at the end).
    {

        // --- lo term: output -= lo mod q ---
        // lo = prod[0..n] is a full integer < 2^256 (not pre-reduced), but
        // mod_sub subtracts mod q, which is the value we want:
        // lambda^2 mod q == (lo + f*hi) mod q == (lo mod q + ...).
        // UNCONTROLLED: no |1>-gated register-sub CCX.
        mod_sub(circ, &prod[0..n], output_reg);

        let pad_hi = circ.alloc_qubit();
        let mut hi_ext: Vec<QubitId> = prod[n..2 * n].to_vec();
        hi_ext.push(pad_hi); // index 2n = overflow slot, |0>

        // --- hi terms: output -= (hi*2^j mod q) for each f-bit j ---
        // Explicit unrolled schedule over F_BITS=[0,4,6,7,8,9,32] (tweak of
        // the doubling schedule for better-square research; equivalent to the
        // prior loop but makes runs of mod_double visible for potential
        // replacement by wider shifts or Litinski-style in follow-ups).
        mod_sub(circ, &hi_ext[0..n], output_reg); // j=0: 0 doubles
        for _ in 0..4 { mod_double(circ, &hi_ext); }
        mod_sub(circ, &hi_ext[0..n], output_reg); // j=4
        for _ in 0..2 { mod_double(circ, &hi_ext); }
        mod_sub(circ, &hi_ext[0..n], output_reg); // j=6
        mod_double(circ, &hi_ext);
        mod_sub(circ, &hi_ext[0..n], output_reg); // j=7
        mod_double(circ, &hi_ext);
        mod_sub(circ, &hi_ext[0..n], output_reg); // j=8
        mod_double(circ, &hi_ext);
        mod_sub(circ, &hi_ext[0..n], output_reg); // j=9
        for _ in 0..23 { mod_double(circ, &hi_ext); }
        mod_sub(circ, &hi_ext[0..n], output_reg); // j=32
        // restore: exactly 32 reverse halvings (4+2+1+1+1+23)
        for _ in 0..32 {
            mod_double_reverse(circ, &hi_ext);
        }

        circ.zero_and_free(pad_hi);
    }

    // Stage 3: uncompute prod (gate-reverse of Stage 1).
    symmetric_square_into_prod_reverse(circ, lambda, prod, z1_state);
}

