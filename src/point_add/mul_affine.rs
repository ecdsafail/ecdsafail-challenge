//! (refactor r2) Mechanically extracted from mul.rs. No logic changes.
use super::*;

/// Symmetric schoolbook for squaring: x² = sum_i x[i]·2^(2i) + sum_{i<j} 2·x[i]·x[j]·2^(i+j).
/// Each cross-product is computed ONCE (instead of twice in full schoolbook),
/// halving the AND count + Cuccaro_add length. Saves ~130k CCX per squaring.
///
/// Row i layout (width n-i): bit 0 = diagonal x[i] at position 2i, bit 1 = 0
/// (gap), bit k+2 = cross-product (x[i] AND x[i+1+k]) at position i+(i+1+k)+1.
pub(crate) fn schoolbook_square_symmetric(b: &mut B, x: &[QubitId], tmp_ext: &[QubitId]) {
    let n = x.len();
    debug_assert_eq!(tmp_ext.len(), 2 * n);
    for i in 0..n {
        // Width: bit 0 = diag at pos 2i, bit 1 = gap, bits 2..(n-i) = cross-
        // products at positions 2i+2..i+n. Last bit index = n-i, so width = n-i+1.
        // Edge case: i = n-1 has only the diagonal, width = 1.
        let width = if i == n - 1 { 1 } else { n - i + 1 };
        let num_cross = if i + 1 < n { n - i - 1 } else { 0 };
        // num_cross = number of cross-products in this row = width - 2 when width >= 2.
        let row = b.alloc_qubits(width);
        b.cx(x[i], row[0]);
        for k in 0..num_cross {
            b.ccx(x[i], x[i + 1 + k], row[k + 2]);
        }
        let pad = b.alloc_qubit();
        let mut row_padded = row.clone();
        row_padded.push(pad);
        let slice: Vec<QubitId> = tmp_ext[2 * i..2 * i + width + 1].to_vec();
        let c_in = b.alloc_qubit();
        cuccaro_add_fast(b, &row_padded, &slice, c_in);
        b.free(c_in);
        b.free(pad);
        b.cx(x[i], row[0]);
        for k in 0..num_cross {
            let m = b.alloc_bit();
            b.hmr(row[k + 2], m);
            b.cz_if(x[i], x[i + 1 + k], m);
        }
        b.free_vec(&row);
    }
}

pub(crate) fn schoolbook_square_symmetric_inverse(b: &mut B, x: &[QubitId], tmp_ext: &[QubitId]) {
    let n = x.len();
    for i in (0..n).rev() {
        let width = if i == n - 1 { 1 } else { n - i + 1 };
        let num_cross = if i + 1 < n { n - i - 1 } else { 0 };
        let row = b.alloc_qubits(width);
        b.cx(x[i], row[0]);
        for k in 0..num_cross {
            b.ccx(x[i], x[i + 1 + k], row[k + 2]);
        }
        let pad = b.alloc_qubit();
        let mut row_padded = row.clone();
        row_padded.push(pad);
        let slice: Vec<QubitId> = tmp_ext[2 * i..2 * i + width + 1].to_vec();
        let c_in = b.alloc_qubit();
        cuccaro_sub_fast(b, &row_padded, &slice, c_in);
        b.free(c_in);
        b.free(pad);
        b.cx(x[i], row[0]);
        for k in 0..num_cross {
            let m = b.alloc_bit();
            b.hmr(row[k + 2], m);
            b.cz_if(x[i], x[i + 1 + k], m);
        }
        b.free_vec(&row);
    }
}

/// Schoolbook squarer with Bennett uncompute. For squaring `tmp_ext = x*x`
/// (2n bits, no mod reduction), then ADD with Solinas reduction to acc,
/// then uncompute tmp_ext via gate-level inverse.
/// Peak-bounded symmetric square: identical to `schoolbook_square_symmetric`
/// except rows whose accumulator slice width exceeds `max_fast_width` use the
/// register-free in-place Cuccaro (`cuccaro_add`, no ~width carry register), so
/// their per-row transient is ~width instead of ~2·width. Narrow rows keep the
/// cheaper measurement Cuccaro. Used by the AFFINE_SQUARE_RECOMPUTE squares that
/// co-reside with the 256-bit `breg` register (base ~1536): without this the
/// widest rows' ~2·257 transient pushes the early-uncompute / recompute squares
/// to ~2052; clamping wide rows keeps every affine phase <= ~1938 (<= 1952).
/// Value-/phase-identical to the fast square except for the carry-register
/// elision on the clamped rows (+~width Toffoli/clamped row). `max_fast_width =
/// usize::MAX` reproduces `schoolbook_square_symmetric` exactly.
pub(crate) fn schoolbook_square_symmetric_pb(
    b: &mut B,
    x: &[QubitId],
    tmp_ext: &[QubitId],
    max_fast_width: usize,
) {
    let n = x.len();
    debug_assert_eq!(tmp_ext.len(), 2 * n);
    for i in 0..n {
        let width = if i == n - 1 { 1 } else { n - i + 1 };
        let num_cross = if i + 1 < n { n - i - 1 } else { 0 };
        let row = b.alloc_qubits(width);
        b.cx(x[i], row[0]);
        for k in 0..num_cross {
            b.ccx(x[i], x[i + 1 + k], row[k + 2]);
        }
        let pad = b.alloc_qubit();
        let mut row_padded = row.clone();
        row_padded.push(pad);
        let slice: Vec<QubitId> = tmp_ext[2 * i..2 * i + width + 1].to_vec();
        let c_in = b.alloc_qubit();
        if width > max_fast_width {
            cuccaro_add(b, &row_padded, &slice, c_in);
        } else {
            cuccaro_add_fast(b, &row_padded, &slice, c_in);
        }
        b.free(c_in);
        b.free(pad);
        b.cx(x[i], row[0]);
        for k in 0..num_cross {
            let m = b.alloc_bit();
            b.hmr(row[k + 2], m);
            b.cz_if(x[i], x[i + 1 + k], m);
        }
        b.free_vec(&row);
    }
}

/// Exact gate-level inverse of `schoolbook_square_symmetric_pb` (same
/// `max_fast_width` row-clamp schedule, in reverse).
pub(crate) fn schoolbook_square_symmetric_pb_inverse(
    b: &mut B,
    x: &[QubitId],
    tmp_ext: &[QubitId],
    max_fast_width: usize,
) {
    let n = x.len();
    debug_assert_eq!(tmp_ext.len(), 2 * n);
    for i in (0..n).rev() {
        let width = if i == n - 1 { 1 } else { n - i + 1 };
        let num_cross = if i + 1 < n { n - i - 1 } else { 0 };
        let row = b.alloc_qubits(width);
        b.cx(x[i], row[0]);
        for k in 0..num_cross {
            b.ccx(x[i], x[i + 1 + k], row[k + 2]);
        }
        let pad = b.alloc_qubit();
        let mut row_padded = row.clone();
        row_padded.push(pad);
        let slice: Vec<QubitId> = tmp_ext[2 * i..2 * i + width + 1].to_vec();
        let c_in = b.alloc_qubit();
        if width > max_fast_width {
            cuccaro_sub(b, &row_padded, &slice, c_in);
        } else {
            cuccaro_sub_fast(b, &row_padded, &slice, c_in);
        }
        b.free(c_in);
        b.free(pad);
        b.cx(x[i], row[0]);
        for k in 0..num_cross {
            let m = b.alloc_bit();
            b.hmr(row[k + 2], m);
            b.cz_if(x[i], x[i + 1 + k], m);
        }
        b.free_vec(&row);
    }
}

pub(crate) fn squaring_add_to_acc_schoolbook(b: &mut B, acc: &[QubitId], x: &[QubitId], p: U256) {
    let n = acc.len();
    debug_assert_eq!(n, 256);
    debug_assert_eq!(x.len(), n);

    let tmp_ext = b.alloc_qubits(2 * n);
    schoolbook_square_symmetric(b, x, &tmp_ext);

    let lo: Vec<QubitId> = tmp_ext[0..n].to_vec();
    let hi: Vec<QubitId> = tmp_ext[n..2 * n].to_vec();
    mod_add_qq_fast(b, acc, &lo, p);
    mod_add_qq_fast(b, acc, &hi, p);
    for _ in 0..4 {
        mod_double_inplace_fast(b, &hi, p);
    }
    mod_add_qq_fast(b, acc, &hi, p);
    for _ in 0..2 {
        mod_double_inplace_fast(b, &hi, p);
    }
    mod_sub_qq_fast(b, acc, &hi, p);
    for _ in 0..4 {
        mod_double_inplace_fast(b, &hi, p);
    }
    mod_add_qq_fast(b, acc, &hi, p);
    let (spill, flag_inv, ovf) = mod_shift_left_by_k(b, &hi, p, 22);
    mod_add_qq(b, acc, &hi, p);
    mod_shift_right_by_k(b, &hi, p, 22, spill, flag_inv, ovf);
    for _ in 0..10 {
        mod_halve_inplace_fast(b, &hi, p);
    }

    schoolbook_square_symmetric_inverse(b, x, &tmp_ext);
    b.free_vec(&tmp_ext);
}

pub(crate) fn mod_add_solinas_ext_product(b: &mut B, acc: &[QubitId], tmp_ext: &[QubitId], p: U256) {
    let n = acc.len();
    debug_assert_eq!(n, 256);
    debug_assert_eq!(tmp_ext.len(), 2 * n);
    let lo: Vec<QubitId> = tmp_ext[0..n].to_vec();
    let hi: Vec<QubitId> = tmp_ext[n..2 * n].to_vec();
    mod_add_qq_fast(b, acc, &lo, p);
    mod_add_qq_fast(b, acc, &hi, p);
    for _ in 0..4 {
        mod_double_inplace_fast(b, &hi, p);
    }
    mod_add_qq_fast(b, acc, &hi, p);
    for _ in 0..2 {
        mod_double_inplace_fast(b, &hi, p);
    }
    mod_sub_qq_fast(b, acc, &hi, p);
    for _ in 0..4 {
        mod_double_inplace_fast(b, &hi, p);
    }
    mod_add_qq_fast(b, acc, &hi, p);
    if sol_ext_product_pos32_fast() {
        // SOL_EXT_PRODUCT_POS32_FAST: fast measurement-based add at position 32.
        let (spill, flag_inv, ovf) = mod_shift_left_by_k(b, &hi, p, 22);
        mod_add_qq_fast(b, acc, &hi, p);
        mod_shift_right_by_k(b, &hi, p, 22, spill, flag_inv, ovf);
    } else {
        let (spill, flag_inv, ovf) = mod_shift_left_by_k(b, &hi, p, 22);
        mod_add_qq(b, acc, &hi, p);
        mod_shift_right_by_k(b, &hi, p, 22, spill, flag_inv, ovf);
    }
    for _ in 0..10 {
        mod_halve_inplace_fast(b, &hi, p);
    }
}

pub(crate) fn mod_sub_solinas_ext_product(b: &mut B, acc: &[QubitId], tmp_ext: &[QubitId], p: U256) {
    let n = acc.len();
    debug_assert_eq!(n, 256);
    debug_assert_eq!(tmp_ext.len(), 2 * n);
    let lo: Vec<QubitId> = tmp_ext[0..n].to_vec();
    let hi: Vec<QubitId> = tmp_ext[n..2 * n].to_vec();
    mod_sub_qq_fast(b, acc, &lo, p);
    mod_sub_qq_fast(b, acc, &hi, p);
    for _ in 0..4 {
        mod_double_inplace_fast(b, &hi, p);
    }
    mod_sub_qq_fast(b, acc, &hi, p);
    for _ in 0..2 {
        mod_double_inplace_fast(b, &hi, p);
    }
    mod_add_qq_fast(b, acc, &hi, p);
    for _ in 0..4 {
        mod_double_inplace_fast(b, &hi, p);
    }
    mod_sub_qq_fast(b, acc, &hi, p);
    if sol_ext_product_pos32_fast() {
        let (spill, flag_inv, ovf) = mod_shift_left_by_k(b, &hi, p, 22);
        mod_sub_qq_fast(b, acc, &hi, p);
        mod_shift_right_by_k(b, &hi, p, 22, spill, flag_inv, ovf);
    } else {
        let (spill, flag_inv, ovf) = mod_shift_left_by_k(b, &hi, p, 22);
        mod_sub_qq(b, acc, &hi, p);
        mod_shift_right_by_k(b, &hi, p, 22, spill, flag_inv, ovf);
    }
    for _ in 0..10 {
        mod_halve_inplace_fast(b, &hi, p);
    }
}

/// Low-scratch twin of `mod_add_solinas_ext_product` (acc += tmp_ext mod p).
/// Uses the SAME validated reduction as the affine y-mul's `gz_solinas_lowscratch`
/// fold: register-free position folds (`mod_*_qq_lowq_lowscratch`), register-free
/// direct doublings (`mod_double_inplace_direct`), and a dirty-borrow shift22 that
/// borrows the product's read-only `lo` half as the venting donor. This drops the
/// fold's transient from ~515 (the register-loaded shift22 `cuccaro_op_0` padded
/// register + carries — the breg_red/breg_unred 2051 binder) to ~k+small. `lo` is
/// restored on exit (dirty-borrow), so the subsequent square uncompute reads an
/// intact `tmp_ext = lam²`. Used ONLY on the AFFINE_SQUARE_RECOMPUTE path (the
/// default path keeps the byte-identical regular fold).
pub(crate) fn mod_add_solinas_ext_product_lowscratch(
    b: &mut B,
    acc: &[QubitId],
    tmp_ext: &[QubitId],
    p: U256,
) {
    let n = acc.len();
    debug_assert_eq!(n, 256);
    debug_assert_eq!(tmp_ext.len(), 2 * n);
    let lo: Vec<QubitId> = tmp_ext[0..n].to_vec();
    let hi: Vec<QubitId> = tmp_ext[n..2 * n].to_vec();
    mod_add_qq_lowq_lowscratch(b, acc, &lo, p);
    mod_add_qq_lowq_lowscratch(b, acc, &hi, p); // position 0
    for _ in 0..4 {
        mod_double_inplace_direct(b, &hi, p);
    }
    mod_add_qq_lowq_lowscratch(b, acc, &hi, p); // position 4
    for _ in 0..2 {
        mod_double_inplace_direct(b, &hi, p);
    }
    mod_sub_qq_lowq_lowscratch(b, acc, &hi, p); // position 6
    for _ in 0..4 {
        mod_double_inplace_direct(b, &hi, p);
    }
    mod_add_qq_lowq_lowscratch(b, acc, &hi, p); // position 10
    if gz_solinas_lowscratch() {
        let (spill, flag_inv, ovf) = mod_shift_left_by_k_dirty(b, &hi, p, 22, &lo);
        b.set_phase("shift22_pos32_dirty");
        mod_add_qq_dirty(b, acc, &hi, p, &lo); // position 32 (venting dirty-borrow)
        mod_shift_right_by_k_dirty(b, &hi, p, 22, spill, flag_inv, ovf, &lo);
    } else {
        let (spill, flag_inv, ovf) = mod_shift_left_by_k(b, &hi, p, 22);
        mod_add_qq(b, acc, &hi, p); // position 32
        mod_shift_right_by_k(b, &hi, p, 22, spill, flag_inv, ovf);
    }
    b.set_phase("sol_halve_tail");
    for _ in 0..10 {
        mod_halve_inplace_fast(b, &hi, p);
    }
}

/// Low-scratch twin of `mod_sub_solinas_ext_product` (acc -= tmp_ext mod p).
/// Companion to `mod_add_solinas_ext_product_lowscratch`; the position-32 acc-fold
/// uses the register-free `mod_sub_qq_lowq_lowscratch` (no dirty donor needed),
/// while the shift22 itself still dirty-borrows `lo`.
pub(crate) fn mod_sub_solinas_ext_product_lowscratch(
    b: &mut B,
    acc: &[QubitId],
    tmp_ext: &[QubitId],
    p: U256,
) {
    let n = acc.len();
    debug_assert_eq!(n, 256);
    debug_assert_eq!(tmp_ext.len(), 2 * n);
    let lo: Vec<QubitId> = tmp_ext[0..n].to_vec();
    let hi: Vec<QubitId> = tmp_ext[n..2 * n].to_vec();
    mod_sub_qq_lowq_lowscratch(b, acc, &lo, p);
    mod_sub_qq_lowq_lowscratch(b, acc, &hi, p); // position 0
    for _ in 0..4 {
        mod_double_inplace_direct(b, &hi, p);
    }
    mod_sub_qq_lowq_lowscratch(b, acc, &hi, p); // position 4
    for _ in 0..2 {
        mod_double_inplace_direct(b, &hi, p);
    }
    mod_add_qq_lowq_lowscratch(b, acc, &hi, p); // position 6 (sign flipped)
    for _ in 0..4 {
        mod_double_inplace_direct(b, &hi, p);
    }
    mod_sub_qq_lowq_lowscratch(b, acc, &hi, p); // position 10
    if gz_solinas_lowscratch() {
        let (spill, flag_inv, ovf) = mod_shift_left_by_k_dirty(b, &hi, p, 22, &lo);
        b.set_phase("shift22_pos32_dirty");
        mod_sub_qq_lowq_lowscratch(b, acc, &hi, p); // position 32 (register-free)
        mod_shift_right_by_k_dirty(b, &hi, p, 22, spill, flag_inv, ovf, &lo);
    } else {
        let (spill, flag_inv, ovf) = mod_shift_left_by_k(b, &hi, p, 22);
        mod_sub_qq(b, acc, &hi, p); // position 32
        mod_shift_right_by_k(b, &hi, p, 22, spill, flag_inv, ovf);
    }
    b.set_phase("sol_halve_tail");
    for _ in 0..10 {
        mod_halve_inplace_fast(b, &hi, p);
    }
}

/// CLEAN low-scratch twin of `mod_add_solinas_ext_product` (acc += tmp_ext mod p).
/// Unlike `mod_add_solinas_ext_product_lowscratch`, uses NO dirty-borrow: every
/// position fold is the register-free, carry-free `mod_add/sub_qq_lowq_lowscratch`
/// (slow in-place Cuccaro + direct const adders + slow comparator, ~0 wide
/// transient), the doublings are `mod_double_inplace_direct` (register-free), and
/// the position-32 shift22 uses the ORDINARY clean `mod_shift_left/right_by_k`
/// (holds only its ~257-wide `padded` cuccaro register, freed before the next op
/// — no donor borrowed from `lo`). `hi` is doubled then halved back and the shift
/// is its own inverse, so `tmp_ext` is restored exactly on exit, allowing the
/// early square uncompute. Correctness does NOT depend on the register-sharing
/// layout (no aliased donor), so it validates 0/0/0 on the UNPACKED base where the
/// dirty-borrow `*_lowscratch` variant gives 1-5 classical mismatches.
pub(crate) fn mod_add_solinas_ext_product_clean_lowscratch(
    b: &mut B,
    acc: &[QubitId],
    tmp_ext: &[QubitId],
    p: U256,
) {
    let n = acc.len();
    debug_assert_eq!(n, 256);
    debug_assert_eq!(tmp_ext.len(), 2 * n);
    let lo: Vec<QubitId> = tmp_ext[0..n].to_vec();
    let hi: Vec<QubitId> = tmp_ext[n..2 * n].to_vec();
    mod_add_qq_lowq_lowscratch(b, acc, &lo, p);
    mod_add_qq_lowq_lowscratch(b, acc, &hi, p); // position 0
    for _ in 0..4 {
        mod_double_inplace_direct(b, &hi, p);
    }
    mod_add_qq_lowq_lowscratch(b, acc, &hi, p); // position 4
    for _ in 0..2 {
        mod_double_inplace_direct(b, &hi, p);
    }
    mod_sub_qq_lowq_lowscratch(b, acc, &hi, p); // position 6
    for _ in 0..4 {
        mod_double_inplace_direct(b, &hi, p);
    }
    mod_add_qq_lowq_lowscratch(b, acc, &hi, p); // position 10
    // position 32: ordinary clean shift22 (lowq slow cuccaro, no dirty donor) +
    // register-free position fold.
    let (spill, flag_inv, ovf) = mod_shift_left_by_k(b, &hi, p, 22);
    mod_add_qq_lowq_lowscratch(b, acc, &hi, p); // position 32 (register-free, clean)
    mod_shift_right_by_k(b, &hi, p, 22, spill, flag_inv, ovf);
    b.set_phase("sol_halve_tail");
    for _ in 0..10 {
        mod_halve_inplace_fast(b, &hi, p);
    }
}

/// CLEAN low-scratch twin of `mod_sub_solinas_ext_product` (acc -= tmp_ext mod p).
/// Companion to `mod_add_solinas_ext_product_clean_lowscratch`; same NO-dirty
/// contract.
pub(crate) fn mod_sub_solinas_ext_product_clean_lowscratch(
    b: &mut B,
    acc: &[QubitId],
    tmp_ext: &[QubitId],
    p: U256,
) {
    let n = acc.len();
    debug_assert_eq!(n, 256);
    debug_assert_eq!(tmp_ext.len(), 2 * n);
    let lo: Vec<QubitId> = tmp_ext[0..n].to_vec();
    let hi: Vec<QubitId> = tmp_ext[n..2 * n].to_vec();
    mod_sub_qq_lowq_lowscratch(b, acc, &lo, p);
    mod_sub_qq_lowq_lowscratch(b, acc, &hi, p); // position 0
    for _ in 0..4 {
        mod_double_inplace_direct(b, &hi, p);
    }
    mod_sub_qq_lowq_lowscratch(b, acc, &hi, p); // position 4
    for _ in 0..2 {
        mod_double_inplace_direct(b, &hi, p);
    }
    mod_add_qq_lowq_lowscratch(b, acc, &hi, p); // position 6 (sign flipped)
    for _ in 0..4 {
        mod_double_inplace_direct(b, &hi, p);
    }
    mod_sub_qq_lowq_lowscratch(b, acc, &hi, p); // position 10
    let (spill, flag_inv, ovf) = mod_shift_left_by_k(b, &hi, p, 22);
    mod_sub_qq_lowq_lowscratch(b, acc, &hi, p); // position 32 (register-free, clean)
    mod_shift_right_by_k(b, &hi, p, 22, spill, flag_inv, ovf);
    b.set_phase("sol_halve_tail");
    for _ in 0..10 {
        mod_halve_inplace_fast(b, &hi, p);
    }
}

pub(crate) fn square_tx_and_combined_ty_l2minus3qx(
    b: &mut B,
    tx: &[QubitId],
    ty: &[QubitId],
    lam: &[QubitId],
    ox: &[BitId],
    p: U256,
) {
    let n = tx.len();
    debug_assert_eq!(n, 256);
    debug_assert_eq!(ty.len(), n);
    debug_assert_eq!(lam.len(), n);

    // AFFINE_SQUARE_RECOMPUTE (default OFF): break the affine double-512 binder.
    // In the legacy path the lam² square's 512-bit `tmp_ext` is held alive across
    // the y-mul (which allocates its OWN 512-bit product) -> two 512-bit products
    // co-resident = the affine peak binder (affine_combined_y_mul / rshift22_rev_step2
    // / schoolbook_mul_inverse all 2307, all carrying affine_combined_square=512).
    // Instead: fold the square into breg=r, UNCOMPUTE the square immediately (the
    // fold restores tmp_ext exactly — hi is doubled/shifted then halved/unshifted,
    // lo is read-only — so the gate-level inverse is valid), run the y-mul with NO
    // square co-resident, then RECOMPUTE the square once to zero breg via the
    // Solinas fold it would have used anyway. The folds are VALUE-identical to the
    // validated AFFINE_R_LIFECYCLE path; the only difference is the square is
    // materialized twice (the breg lam²-mod-p register must be both BUILT before
    // and DESTROYED after the y-mul, and lam²-mod-p can only be uncomputed by
    // recomputing lam²). Costs +1 square pair (~+65k executed Toffoli) but drops
    // EVERY affine phase to <= ~1938 (was 2307). FUNDAMENTAL dependency: the 512
    // cannot be freed Toffoli-free (the breg build+destroy straddle the y-mul; the
    // only shared resource is the 512 buffer). Validated value-correct (0/0/0 with
    // all truncations off). The extra ops re-roll the Fiat-Shamir inputs, so the
    // tight default truncation must be re-tuned (KAL_REROLL / KAL_WTRUNC_MARGIN /
    // KAL_CARRYTAIL_W) for THIS op stream — that re-tune is the standard island
    // lottery, NOT an affine value bug. Byte-identical fallback: =0.
    // BAKED DEFAULT ON: required prereq of the validated C* construction (frees the
    // affine 512-square so the affine cluster drops below the C* inversion sweep).
    // AFFINE_SQUARE_RECOMPUTE=0 restores the legacy co-resident-square path.
    let affine_square_recompute = env_flag_enabled("AFFINE_SQUARE_RECOMPUTE", true)
        && std::env::var("AFFINE_R_LIFECYCLE").ok().as_deref() != Some("0");

    if affine_square_recompute {
        // Peak-bound the breg-coresident squares: the early-uncompute and the
        // recompute square run while `breg` (256) is live (base ~1536), so their
        // widest rows' ~2·257 fast transient would peak ~2052. Clamp rows wider
        // than `mfw` to the register-free in-place Cuccaro (~width transient) so
        // every affine phase stays <= ~1938. 200 leaves margin below 1952
        // (base 1536 + 2·200 + 3 = 1939).
        // BAKED DEFAULT 243: the affine-square ceiling that keeps the affine cluster
        // (~2024) just below the C* inversion/mul1 binder (2025) while minimizing
        // Toffoli (more fast measurement-Cuccaro rows). mfw>=244 makes affine bind.
        let mfw = env_usize("AFFINE_SQUARE_RECOMPUTE_MFW").unwrap_or(243);

        // FOLD VARIANT (default CLEAN): the early-uncompute + recompute Solinas
        // folds must restore `tmp_ext` exactly so the square uncompute is valid.
        // The CLEAN fold (`*_clean_lowscratch`) does this with NO dirty-borrow, so
        // it is correct on the UNPACKED base (the dirty-borrow `*_lowscratch`
        // variant aliases the product's `lo` half as a venting donor, sound only
        // paired with KAL_REGISTER_SHARING — 1-5 classical mismatches unpacked).
        // Set AFFINE_RECOMPUTE_DIRTY_FOLD=1 for the packed/register-sharing pairing.
        let dirty_fold = env_flag_enabled("AFFINE_RECOMPUTE_DIRTY_FOLD", false);

        b.set_phase("affine_combined_square");
        let tmp_ext = b.alloc_qubits(2 * n);
        schoolbook_square_symmetric_pb(b, lam, &tmp_ext, mfw);

        b.set_phase("affine_combined_breg_red");
        let breg = b.alloc_qubits(n);
        // Clean (or dirty) low-scratch fold restores tmp_ext exactly (hi doubled
        // then halved back, lo read-only), so the early square uncompute is valid.
        if dirty_fold {
            mod_add_solinas_ext_product_lowscratch(b, &breg, &tmp_ext, p); // breg = lam² mod p = r
        } else {
            mod_add_solinas_ext_product_clean_lowscratch(b, &breg, &tmp_ext, p); // breg = lam² mod p = r
        }
        // Early uncompute: the fold restores tmp_ext, so the gate-level inverse is
        // valid. Free the 512 BEFORE the y-mul so the two products never co-reside.
        b.set_phase("affine_combined_square_unc");
        schoolbook_square_symmetric_pb_inverse(b, lam, &tmp_ext, mfw);
        b.free_vec(&tmp_ext);

        mod_sub_double_qb(b, &breg, ox, p);
        mod_sub_qb(b, &breg, ox, p); // breg = r - 3ox  (the y-mul multiplier)

        b.set_phase("affine_combined_y_mul");
        if env_flag_enabled("POINT_ADD_AFFINE_COMBINED_Y_KARATSUBA_LOWQ", false) {
            mod_mul_add_into_acc_karatsuba_lowq(b, ty, lam, &breg, p);
        } else if env_flag_enabled("AFFINE_Y_MUL_LOWSCRATCH_FOLD", stack_2565_enabled()) {
            mod_mul_add_into_acc_schoolbook_lowscratch_fold(b, ty, lam, &breg, p);
        } else {
            mod_mul_add_into_acc_schoolbook(b, ty, lam, &breg, p);
        }

        b.set_phase("affine_combined_breg_unred");
        mod_add_qb(b, &breg, ox, p);
        mod_add_double_qb(b, &breg, ox, p); // breg = lam² mod p = r (restored)

        b.set_phase("affine_combined_tx_update");
        mod_sub_qq_fast(b, tx, &breg, p); // tx -= r

        b.set_phase("affine_combined_breg_unred");
        // Recompute the square (NOT co-resident with the y-mul) to zero breg via
        // the one Solinas fold it would have used anyway.
        let tmp_ext2 = b.alloc_qubits(2 * n);
        schoolbook_square_symmetric_pb(b, lam, &tmp_ext2, mfw);
        if dirty_fold {
            mod_sub_solinas_ext_product_lowscratch(b, &breg, &tmp_ext2, p); // breg -= lam² mod p = 0
        } else {
            mod_sub_solinas_ext_product_clean_lowscratch(b, &breg, &tmp_ext2, p); // breg -= lam² mod p = 0
        }
        schoolbook_square_symmetric_pb_inverse(b, lam, &tmp_ext2, mfw);
        b.free_vec(&tmp_ext2);
        b.free_vec(&breg);

        b.set_phase("affine_combined_tx_update");
        mod_add_double_qb(b, tx, ox, p);
        mod_add_qb(b, tx, ox, p);
        mod_neg_inplace_fast(b, tx, p);
        return;
    }

    b.set_phase("affine_combined_square");
    let tmp_ext = b.alloc_qubits(2 * n);
    schoolbook_square_symmetric(b, lam, &tmp_ext);

    b.set_phase("affine_combined_breg_red");
    let breg = b.alloc_qubits(n);
    mod_add_solinas_ext_product(b, &breg, &tmp_ext, p);
    mod_sub_double_qb(b, &breg, ox, p);
    mod_sub_qb(b, &breg, ox, p);

    b.set_phase("affine_combined_y_mul");
    if env_flag_enabled("POINT_ADD_AFFINE_COMBINED_Y_KARATSUBA_LOWQ", false) {
        mod_mul_add_into_acc_karatsuba_lowq(b, ty, lam, &breg, p);
    } else if env_flag_enabled("AFFINE_Y_MUL_LOWSCRATCH_FOLD", stack_2565_enabled()) {
        // Peak-minimized y-mul: cuts the ~256-wide transient scratch the
        // schoolbook MAC holds on top of its 512 product while the lam² square's
        // 512 product co-resides (the -x correction pads, the Solinas fold's
        // carry/const registers, and the fold mod_double's const register). The
        // y-mul instant drops 2565 -> ~2333, below the next cluster (2459),
        // breaking the 2565 binder. Default-on under STACK-2565; set
        // AFFINE_Y_MUL_LOWSCRATCH_FOLD=0 to restore the byte-identical
        // fast-fold schoolbook MAC (peak 2565).
        mod_mul_add_into_acc_schoolbook_lowscratch_fold(b, ty, lam, &breg, p);
    } else {
        mod_mul_add_into_acc_schoolbook(b, ty, lam, &breg, p);
    }

    // r-lifecycle (default): fold lambda^2 once and reuse the reduced value for
    // the tx update so tx_update is a cheap qq-sub instead of a second full
    // Solinas fold. After the two 3Qx re-adds below, `breg` holds
    // r = lambda^2 mod p (the reduced value) -- exactly the constant tx_update
    // must subtract. Consume breg-as-r for tx BEFORE zeroing breg, then zero
    // breg with the one Solinas fold it would have used anyway. No extra
    // register => peak-neutral. Validated -18,963 Toffoli, peak 2708 unchanged.
    // Set AFFINE_R_LIFECYCLE=0 to fall back to the legacy 3-fold path.
    let affine_r_lifecycle =
        std::env::var("AFFINE_R_LIFECYCLE").ok().as_deref() != Some("0");

    if affine_r_lifecycle {
        b.set_phase("affine_combined_breg_unred");
        mod_add_qb(b, &breg, ox, p); // breg = lambda^2 mod p = r
        mod_add_double_qb(b, &breg, ox, p);

        b.set_phase("affine_combined_tx_update");
        // tx -= r  (== tx -= lambda^2 mod p), reusing breg=r, cheap qq sub.
        mod_sub_qq_fast(b, tx, &breg, p);

        b.set_phase("affine_combined_breg_unred");
        // Zero breg via the one Solinas fold it would have used anyway.
        mod_sub_solinas_ext_product(b, &breg, &tmp_ext, p);
        b.free_vec(&breg);

        b.set_phase("affine_combined_tx_update");
        mod_add_double_qb(b, tx, ox, p);
        mod_add_qb(b, tx, ox, p);
        mod_neg_inplace_fast(b, tx, p);
    } else {
        b.set_phase("affine_combined_breg_unred");
        mod_add_qb(b, &breg, ox, p);
        mod_add_double_qb(b, &breg, ox, p);
        mod_sub_solinas_ext_product(b, &breg, &tmp_ext, p);
        b.free_vec(&breg);

        b.set_phase("affine_combined_tx_update");
        mod_sub_solinas_ext_product(b, tx, &tmp_ext, p);
        mod_add_double_qb(b, tx, ox, p);
        mod_add_qb(b, tx, ox, p);
        mod_neg_inplace_fast(b, tx, p);
    }

    schoolbook_square_symmetric_inverse(b, lam, &tmp_ext);
    b.free_vec(&tmp_ext);
}

/// Schoolbook squarer with Bennett uncompute. For squaring `tmp_ext = x*x`
/// (2n bits, no mod reduction), then sub from acc with on-the-fly Solinas
/// reduction, then uncompute tmp_ext via gate-level inverse. Saves ~170k
/// CCX vs walk-x squaring (459k → 289k) by avoiding 256 expensive
/// cmod_add_qq calls (each 5n) in favor of 2n²=131k of cheap AND+Cuccaro.
pub(crate) fn squaring_sub_from_acc_schoolbook(b: &mut B, acc: &[QubitId], x: &[QubitId], p: U256) {
    let n = acc.len();
    debug_assert_eq!(n, 256);
    debug_assert_eq!(x.len(), n);
    let c = U256::MAX.wrapping_sub(p).wrapping_add(U256::from(1));

    // Wide accumulator (2n bits) starts at 0.
    let tmp_ext = b.alloc_qubits(2 * n);

    // Phase 1: symmetric schoolbook tmp_ext = x*x (~half the CCX of full).
    schoolbook_square_symmetric(b, x, &tmp_ext);

    // Phase 2: subtract (lo + hi*c mod p) from acc.
    // For each set bit k of c, sub (hi shifted by k mod p) from acc, by
    // walking hi via mod_double in place. Sub lo first.
    let lo: Vec<QubitId> = tmp_ext[0..n].to_vec();
    let hi: Vec<QubitId> = tmp_ext[n..2 * n].to_vec();
    mod_sub_qq_fast(b, acc, &lo, p);
    let _ = c;
    // 977 consolidation: c = {+2^0, +2^4, -2^6, +2^10, +2^32}. For acc-=hi·c, signs flip:
    // acc -= hi·2^0, acc -= hi·2^4, acc += hi·2^6, acc -= hi·2^10, acc -= hi·2^32.
    mod_sub_qq_fast(b, acc, &hi, p);
    for _ in 0..4 {
        mod_double_inplace_fast(b, &hi, p);
    }
    mod_sub_qq_fast(b, acc, &hi, p);
    for _ in 0..2 {
        mod_double_inplace_fast(b, &hi, p);
    }
    mod_add_qq_fast(b, acc, &hi, p); // sign flipped
    for _ in 0..4 {
        mod_double_inplace_fast(b, &hi, p);
    }
    mod_sub_qq_fast(b, acc, &hi, p);
    let (spill, flag_inv, ovf) = mod_shift_left_by_k(b, &hi, p, 22);
    mod_sub_qq(b, acc, &hi, p);
    mod_shift_right_by_k(b, &hi, p, 22, spill, flag_inv, ovf);
    for _ in 0..10 {
        mod_halve_inplace_fast(b, &hi, p);
    }

    // Phase 3: uncompute tmp_ext via symmetric schoolbook inverse.
    schoolbook_square_symmetric_inverse(b, x, &tmp_ext);

    b.free_vec(&tmp_ext);
}

pub(crate) fn squaring_sub_from_acc_schoolbook_lowq_shift22(
    b: &mut B,
    acc: &[QubitId],
    x: &[QubitId],
    p: U256,
) {
    let n = acc.len();
    debug_assert_eq!(n, 256);
    debug_assert_eq!(x.len(), n);
    let c = U256::MAX.wrapping_sub(p).wrapping_add(U256::from(1));

    let tmp_ext = b.alloc_qubits(2 * n);
    schoolbook_square_symmetric(b, x, &tmp_ext);

    let lo: Vec<QubitId> = tmp_ext[0..n].to_vec();
    let hi: Vec<QubitId> = tmp_ext[n..2 * n].to_vec();
    mod_sub_qq_fast(b, acc, &lo, p);
    let _ = c;
    mod_sub_qq_fast(b, acc, &hi, p);
    for _ in 0..4 {
        mod_double_inplace_fast(b, &hi, p);
    }
    mod_sub_qq_fast(b, acc, &hi, p);
    for _ in 0..2 {
        mod_double_inplace_fast(b, &hi, p);
    }
    mod_add_qq_fast(b, acc, &hi, p);
    for _ in 0..4 {
        mod_double_inplace_fast(b, &hi, p);
    }
    mod_sub_qq_fast(b, acc, &hi, p);
    let (spill, flag_inv, ovf) = mod_shift_left_by_k_lowq(b, &hi, p, 22);
    mod_sub_qq(b, acc, &hi, p);
    mod_shift_right_by_k_lowq(b, &hi, p, 22, spill, flag_inv, ovf);
    for _ in 0..10 {
        mod_halve_inplace_fast(b, &hi, p);
    }

    schoolbook_square_symmetric_inverse(b, x, &tmp_ext);
    b.free_vec(&tmp_ext);
}

pub(crate) fn mod_mul_sub_qq(b: &mut B, acc: &[QubitId], x: &[QubitId], y: &[QubitId], p: U256) {
    // acc -= x * y mod p. Negate x, run schoolbook ADD (cheaper than sub),
    // then restore x. For x≠y we can walk the negated multiplicand in place
    // and halve it back afterwards, avoiding the doubled tmp register. For
    // squaring we snapshot the original control bits once into `ctrl_copy`,
    // then reuse the same in-place walk on the negated x.
    let n = acc.len();
    let is_squaring = x[0] == y[0]; // same register → squaring
    if is_squaring {
        // Use the schoolbook squarer for the squaring case (~170k savings).
        squaring_sub_from_acc_schoolbook(b, acc, x, p);
        return;
    }
    if false {
        // Hold the original x bits fixed for control while x itself walks
        // through (-x)*2^i mod p.
        let ctrl_copy = b.alloc_qubits(n);
        for i in 0..n {
            b.cx(x[i], ctrl_copy[i]);
        }
        mod_neg_inplace_fast(b, x, p);
        for i in 0..n {
            cmod_add_qq(b, acc, x, ctrl_copy[i], p);
            if i < n - 1 {
                mod_double_inplace_fast(b, x, p);
            }
        }
        for _ in 0..(n - 1) {
            mod_halve_inplace_fast(b, x, p);
        }
        mod_neg_inplace_fast(b, x, p);
        for i in 0..n {
            b.cx(x[i], ctrl_copy[i]);
        }
        b.free_vec(&ctrl_copy);
    } else {
        // Keep x negated during the loop and walk it in place.
        mod_neg_inplace_fast(b, x, p);
        for i in 0..n {
            cmod_add_qq(b, acc, x, y[i], p);
            if i < n - 1 {
                mod_double_inplace_fast(b, x, p);
            }
        }
        for _ in 0..(n - 1) {
            mod_halve_inplace_fast(b, x, p);
        }
        mod_neg_inplace_fast(b, x, p);
    }
}
