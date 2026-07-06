//! Modular square-subtract `output -= lambda^2 mod q` (secp256k1) used by the
//! EC point-add, built on the sibling `super::arith` mod-sub / mod-double
//! primitives.
//!
//! - [`symmetric_square_into_prod`]: the symmetric schoolbook square -- each
//!   cross-product x_i*x_j once, ~n^2/2 CCX. The row-add is `arith::
//!   hybrid_add_adaptive`; the cross ANDs are uncomputed by `clear_and` (HMR +
//!   conditional-Z), the diagonal by `cx`.
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

/// Build `prod[0..2n] += value(x[0..n])^2` (integer, no reduction) via the
/// symmetric schoolbook square: each off-diagonal cross-product x[i]*x[j] (i<j)
/// is computed once, halving the AND/Toffoli count vs the full schoolbook.
///
///   x^2 = sum_i x[i]*2^(2i)  +  sum_{i<j} 2*x[i]*x[j]*2^(i+j)
///
/// Row `i` (added at product position 2i):
///   bit 0      = diagonal x[i]               (pos 2i)         via CX
///   bit 1      = 0 (gap)
///   bit k+2    = cross x[i] AND x[i+1+k]      (pos 2i+2+k)     via CCX
///
/// `prod` is grown lazily (only up to the highest bit written so far) so the
/// per-row register recycles the not-yet-allocated high slots. Pass an empty Vec.
fn symmetric_square_into_prod(circ: &mut B, x: &[QubitId], prod: &mut Vec<QubitId>) {
    let n = x.len();
    assert!(prod.is_empty(), "prod is grown lazily; pass an empty Vec");
    for i in 0..n {
        // Row i has (n-1-i) crosses; the top cross lands at row-bit (n-1-i)+1 =
        // n-i, so width = n-i+1 (i == n-1: only the diagonal, width 1).
        let num_cross = n.saturating_sub(i + 1);
        let width = if i == n - 1 { 1 } else { n - i + 1 };
        // Row-add writes prod[2i .. 2i+width+1] (one carry slot). Grow prod up to
        // the highest bit written so far.
        let hi = (2 * i + width + 1).min(2 * n);
        while prod.len() < hi {
            prod.push(circ.alloc_qubit());
        }
        let row: Vec<QubitId> = (0..width).map(|_| circ.alloc_qubit()).collect();
        circ.cx(x[i], row[0]); // diagonal
        for k in 0..num_cross {
            circ.ccx(x[i], x[i + 1 + k], row[k + 2]); // cross x[i] & x[i+1+k]
        }
        add_into(circ, &prod[2 * i..hi], &row);
        // Uncompute the row: each cross `row[k+2] = x[i] AND x[i+1+k]` is a clean
        // AND (add_into restored `row`), so measurement-vent it (clear_and: HMR +
        // cz, 0 Toffoli) instead of a reverse ccx. The diagonal is a CX.
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

/// Gate-reverse of [`symmetric_square_into_prod`]: rebuilds each row and
/// SUBTRACTS it from `prod`, draining `prod` back to |0>. Rows run in reverse
/// order; `prod` is freed lazily (mirror of the forward lazy growth).
fn symmetric_square_into_prod_reverse(circ: &mut B, x: &[QubitId], mut prod: Vec<QubitId>) {
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
        // subtract the row (X-sandwiched add).
        for q in &prod[2 * i..hi] {
            circ.x(*q);
        }
        add_into(circ, &prod[2 * i..hi], &row);
        for q in &prod[2 * i..hi] {
            circ.x(*q);
        }
        // Vent the cross AND-uncompute (clean ANDs; see the forward build).
        for k in 0..num_cross {
            clear_and(circ, &row[k + 2], &x[i], &x[i + 1 + k]);
        }
        circ.cx(x[i], row[0]);
        for q in row {
            circ.zero_and_free(q);
        }
        // Rows below i reach at most prod index n+i, so all indices > n+i are now
        // |0> and can be freed (mirror of the forward lazy growth).
        let keep = (n + i + 1).min(2 * n);
        while prod.len() > keep {
            circ.zero_and_free(prod.pop().unwrap());
        }
    }
    for q in prod {
        circ.zero_and_free(q);
    }
}

/// Unconditional `output_reg -= lambda^2 mod q` (secp256k1), normal throughout.
///
/// `lambda` is `n = 256` bits (lambda < q); `output_reg` is `n = 256` bits and
/// holds a value < q on entry (the EC-add keeps output reduced).
///
/// Stage 1: build the 2n-bit integer product `prod = lambda^2`
/// with [`symmetric_square_into_prod`] (~n(n-1)/2 CCX).
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
    let mut prod: Vec<QubitId> = Vec::with_capacity(2 * n);
    symmetric_square_into_prod(circ, lambda, &mut prod);

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
        let mut shifted = 0usize;
        for &j in &F_BITS {
            while shifted < j {
                mod_double(circ, &hi_ext); // hi *= 2 mod q (257-bit, pad restored)
                shifted += 1;
            }
            mod_sub(circ, &hi_ext[0..n], output_reg);
        }
        // restore hi to its raw value (matched-reverse the doublings). mod_double
        // is built from a swap chain (self-inverse) + a +f window (NOT self-
        // inverse), so a literal re-run does NOT undo it -- use the value-inverse
        // mod_double_reverse (one halving per forward doubling).
        while shifted > 0 {
            mod_double_reverse(circ, &hi_ext);
            shifted -= 1;
        }

        circ.zero_and_free(pad_hi);
    }

    // Stage 3: uncompute prod (gate-reverse of Stage 1).
    symmetric_square_into_prod_reverse(circ, lambda, prod);
}

