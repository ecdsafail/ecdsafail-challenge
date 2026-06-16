//! Curve-generic "dialog" in-place modular multiplication: the GCD pack
//! (`gcd_pack.rs`) + Bezout unpack (`bezout_unpack.rs`) + `IPModMul`
//! (`mod_mul_eea.rs`), parameterized by curve.
//!
//! The GCD pack is PURELY integer arithmetic on `(u, v)` — its only
//! curve dependency is the constant `q` loaded into `u_full`, so it takes
//! `(n, q_bytes)` directly. The apply / mul direction invokes the curve's
//! field-modular primitives through the `ModRed` trait.
//!
//! Direct curve-generic port of the validated secp256k1 functions:
//!   - `gcd_pack::forward_gcd_pack_quantum_secp256k1` (and `_reverse`)
//!   - `bezout_unpack::apply_bitvector_quantum_secp256k1` (and `_inv`)
//!   - `mod_mul_eea::{mod_mul_in_place_eea_secp256k1, _reverse,
//!     clear_zeroed_drift_reg_secp256k1}`

use crate::point_add::trailmix_port::arith::schrottenloher::gcd_pack::{expected_iterations, u_padding};
use crate::point_add::trailmix_port::circuit::{Circuit, QReg};

use super::mod_arith::ModRed;

/// Dialog window geometry, MUST match `gcd_pack::{DIALOG_M, DIALOG_PACK}`
/// (verified: `gcd_pack.rs:19-20`).
const DIALOG_M: usize = 5;
const DIALOG_PACK: usize = 8;

/// Measurement-vent budget for the GCD csub's controlled hybrid adder.
/// `usize::MAX` clamps to the full Gidney 2n adder (the GCD has full
/// ancilla headroom; the peak is set by `apply_bv`, not the GCD).
const GCD_CSUB_VENTS: usize = usize::MAX;

// ===========================================================================
// Forward / reverse GCD pack — parameterized by (n, q_bytes).
// Curve-generic copy of gcd_pack::forward_gcd_pack_quantum_secp256k1{,_reverse}.
// ===========================================================================

/// Quantum forward GCD packing. Drives `(u, v)` from `(q, x)` to `(1, 0)`,
/// packing per-iter `(b0, b0_and_b1)` choices into `garbage`.
///
/// Copy of `gcd_pack::forward_gcd_pack_quantum_secp256k1` with `n` and the
/// `u_full = q` init taken as params (`q_bytes` little-endian). Everything
/// else (`current_n` schedule, cswap, csub, `right_shift`, comparator, the
/// `gcd_compress5` window absorb, the strict-dealloc truncation) is
/// curve-agnostic.
pub fn forward_gcd_pack_quantum(
    circ: &mut Circuit,
    v_full: &mut Vec<QReg>,
    n: usize,
    q_bytes: &[u8],
) -> Vec<QReg> {
    use crate::point_add::trailmix_port::arith::schrottenloher::gcd_compress5::compress_5iter_refs;
    use crate::point_add::trailmix_port::arith::schrottenloher::msb_compare;

    let pad = u_padding(n);
    assert!(
        v_full.len() >= n,
        "v_full must be at least n bits (pad above n is never active in the current_n schedule)"
    );
    let iters = expected_iterations(n);
    let garbage_len = iters / DIALOG_M * DIALOG_PACK;

    let prev = circ.push_section("fwd_gcd");

    // Load u_full = q (n bits) || 0 (pad bits). Mutable so we can
    // truncate as `current_n` shrinks per iter (strict-dealloc
    // requires inactive tail bits to be freed before the comparator
    // allocates fresh ancillae in subsequent iters).
    let mut u_full = circ.alloc_qreg_bits("u_full", n);
    for i in 0..n {
        let byte_idx = i / 8;
        if byte_idx < q_bytes.len() && ((q_bytes[byte_idx] >> (i % 8)) & 1) == 1 {
            circ.x(&u_full[i]);
        }
    }

    // Per-iter ancillae (reused).
    let b0 = circ.alloc_qreg("gcd.b0");
    let b0_and_b1 = circ.alloc_qreg("gcd.b0_and_b1");

    // Paper Sec. 3.1 TRUNCATE constant. Combined with u_padding the
    // approximate compare looks at top (TRUNCATE + u_padding) bits of the
    // CURRENT (shrunken) effective register width.
    const TRUNCATE: usize = 40;
    let trunc = TRUNCATE + pad;

    // Incremental garbage register. The freed bits of u_full enter the
    // allocator pool as iters progress; per-pack allocs reuse them.
    let mut garbage_vec: Vec<QReg> = Vec::with_capacity(garbage_len);
    // The current DIALOG_M-iter dialog window, held DECOMPRESSED (2*M bits)
    // and compressed ONCE at the window end.
    let mut cur_win: Vec<QReg> = Vec::new();

    const STEP: usize = 1;
    let current_n_at = |i: usize| -> usize {
        let raw = (n as f64) - (i as f64) * 0.5 * 1.415_f64 + pad as f64;
        let raw_int = raw.ceil().max(0.0) as usize;
        let stepped = raw_int.div_ceil(STEP) * STEP;
        stepped.min(n)
    };

    for i in 0..iters {
        let iter_section = format!("fwd_iter_{i:04}");
        let iter_prev = circ.push_section(&iter_section);

        let current_n = current_n_at(i);

        // Free u_full AND v_full tail bits that just became inactive.
        while u_full.len() > current_n.max(1) {
            let q = u_full.pop().expect("u_full nonempty");
            circ.zero_and_free(q);
        }
        while v_full.len() > current_n.max(1) {
            let q = v_full.pop().expect("v_full nonempty");
            circ.zero_and_free(q);
        }

        // Top `trunc` bits inside the shrunken view.
        let cmp_eff = trunc.min(current_n);
        let cmp_lo = current_n.saturating_sub(cmp_eff);

        // 1) b0 = v_full[0] (parity of v).
        circ.cx(&v_full[0], &b0);

        // 2) b0_and_b1 ^= b0 AND (v_top < u_top) on top-`trunc` bits.
        msb_compare::controlled_lt_msbs_gidney(
            circ,
            &b0,
            &v_full[cmp_lo..current_n],
            &u_full[cmp_lo..current_n],
            cmp_eff,
            &b0_and_b1,
        );

        // 3) cswap(b0_and_b1, u_full, v_full) over `current_n` bits.
        let s_cs = circ.push_section("cswap");
        for j in 0..current_n {
            circ.cswap(&b0_and_b1, &u_full[j], &v_full[j]);
        }
        circ.pop_section(&s_cs);

        // 4) v_full -= b0 * u_full over `current_n` bits.
        let s_sub = circ.push_section("csub");
        let v_slice = &v_full[..current_n];
        let u_slice = &u_full[..current_n];
        for q in v_slice {
            circ.x(q);
        }
        crate::point_add::trailmix_port::arith::gidney_const_adder::controlled_hybrid_add(
            circ,
            &b0,
            v_slice,
            u_slice,
            GCD_CSUB_VENTS,
        );
        for q in v_slice {
            circ.x(q);
        }
        circ.pop_section(&s_sub);

        // 5) v_full[..current_n] >>= 1.
        crate::point_add::trailmix_port::arith::shift::right_shift(circ, v_slice);

        // 6-7) Per-window dialog packing.
        let slot = i % DIALOG_M;
        if slot == 0 {
            cur_win = (0..2 * DIALOG_M)
                .map(|_| circ.alloc_qreg("eea_win"))
                .collect();
        }
        circ.swap(&b0, &cur_win[2 * slot]);
        circ.swap(&b0_and_b1, &cur_win[2 * slot + 1]);
        if slot == DIALOG_M - 1 {
            let ds: Vec<QReg> = (0..8).map(|_| circ.alloc_qreg("c5_dirty")).collect();
            {
                let win_refs: [&QReg; 10] = std::array::from_fn(|k| &cur_win[k]);
                let dirty: Vec<&QReg> = ds.iter().collect();
                compress_5iter_refs(circ, &win_refs, &dirty);
            }
            for q in ds.into_iter().rev() {
                circ.zero_and_free(q);
            }
            while cur_win.len() > DIALOG_PACK {
                let freed = cur_win.pop().expect("high window bit |0> post-compress");
                circ.zero_and_free(freed);
            }
            garbage_vec.append(&mut cur_win);
        }

        // 8) Strict-dealloc touch.
        for k in 0..current_n {
            circ.cx(&b0, &u_full[k]);
            circ.cx(&b0_and_b1, &u_full[k]);
        }

        circ.pop_section(&iter_prev);
    }

    // End: drive u_full[0] toward 0.
    circ.x(&u_full[0]);

    while v_full.len() > 1 {
        let q = v_full.pop().expect("v_full nonempty");
        circ.zero_and_free(q);
    }

    drop(b0);
    drop(b0_and_b1);
    drop(u_full);
    assert_eq!(garbage_vec.len(), garbage_len);
    circ.pop_section(&prev);
    garbage_vec
}

/// Reverse of `forward_gcd_pack_quantum`. Restores `v_full[..n]` to its
/// original value and drives `garbage` back to all-|0>.
///
/// Copy of `gcd_pack::forward_gcd_pack_quantum_secp256k1_reverse` with `n`
/// and the `q_bytes` de-init taken as params.
pub fn forward_gcd_pack_quantum_reverse(
    circ: &mut Circuit,
    v_full: &mut Vec<QReg>,
    garbage: &mut Vec<QReg>,
    n: usize,
    q_bytes: &[u8],
) {
    use crate::point_add::trailmix_port::arith::schrottenloher::gcd_compress5::compress_5iter_reverse_refs;
    use crate::point_add::trailmix_port::arith::schrottenloher::msb_compare;

    let pad = u_padding(n);
    let iters = expected_iterations(n);
    let garbage_len = iters / DIALOG_M * DIALOG_PACK;
    assert!(garbage.len() >= garbage_len);

    let prev = circ.push_section("fwd_gcd_rev");

    let mut u_full: Vec<QReg> = Vec::with_capacity(n);
    u_full.push(circ.alloc_qreg("u_full_re[0]"));
    circ.x(&u_full[0]);

    let b0 = circ.alloc_qreg("gcd_re.b0");
    let b0_and_b1 = circ.alloc_qreg("gcd_re.b0_and_b1");

    const TRUNCATE: usize = 40;
    let trunc = TRUNCATE + pad;
    const STEP: usize = 1;
    let current_n_at = |i: usize| -> usize {
        let raw = (n as f64) - (i as f64) * 0.5 * 1.415_f64 + pad as f64;
        let raw_int = raw.ceil().max(0.0) as usize;
        let stepped = raw_int.div_ceil(STEP) * STEP;
        stepped.min(n)
    };
    let mut cur_win: Vec<QReg> = Vec::new();

    for i in (0..iters).rev() {
        let iter_section = format!("rev_iter_{i:04}");
        let iter_prev = circ.push_section(&iter_section);

        let current_n = current_n_at(i).max(1);
        while u_full.len() < current_n {
            u_full.push(circ.alloc_qreg("u_full_re"));
        }
        while v_full.len() < current_n {
            v_full.push(circ.alloc_qreg("v_full_re"));
        }

        let slot = i % DIALOG_M;
        if slot == DIALOG_M - 1 {
            let glen = garbage.len();
            cur_win = garbage.split_off(glen - DIALOG_PACK);
            for _ in 0..(2 * DIALOG_M - DIALOG_PACK) {
                cur_win.push(circ.alloc_qreg("eea_win_re"));
            }
            let ds: Vec<QReg> = (0..8).map(|_| circ.alloc_qreg("c5_dirty")).collect();
            {
                let win_refs: [&QReg; 10] = std::array::from_fn(|k| &cur_win[k]);
                let dirty: Vec<&QReg> = ds.iter().collect();
                compress_5iter_reverse_refs(circ, &win_refs, &dirty);
            }
            for q in ds.into_iter().rev() {
                circ.zero_and_free(q);
            }
        }
        circ.swap(&b0, &cur_win[2 * slot]);
        circ.swap(&b0_and_b1, &cur_win[2 * slot + 1]);

        let v_slice = &v_full[..current_n];
        crate::point_add::trailmix_port::arith::shift::left_shift(circ, v_slice);

        let s_add = circ.push_section("rev_cadd");
        let u_slice = &u_full[..current_n];
        crate::point_add::trailmix_port::arith::gidney_const_adder::controlled_hybrid_add(
            circ,
            &b0,
            v_slice,
            u_slice,
            GCD_CSUB_VENTS,
        );
        circ.pop_section(&s_add);

        let s_cs = circ.push_section("rev_cswap");
        for j in 0..current_n {
            circ.cswap(&b0_and_b1, &u_full[j], &v_full[j]);
        }
        circ.pop_section(&s_cs);

        let cmp_eff = trunc.min(current_n);
        let cmp_lo = current_n.saturating_sub(cmp_eff);
        msb_compare::controlled_lt_msbs_gidney(
            circ,
            &b0,
            &v_full[cmp_lo..current_n],
            &u_full[cmp_lo..current_n],
            cmp_eff,
            &b0_and_b1,
        );

        circ.cx(&v_full[0], &b0);

        for k in 0..current_n {
            circ.cx(&b0, &u_full[k]);
            circ.cx(&b0_and_b1, &u_full[k]);
        }

        circ.pop_section(&iter_prev);

        if i % DIALOG_M == 0 {
            for q in cur_win.drain(..) {
                circ.zero_and_free(q);
            }
        }
    }
    debug_assert!(garbage.is_empty(), "garbage tape not fully drained");

    // After all reverse iters, u_full should be q; de-init by X-flipping
    // the set bits, then free as |0>.
    for i in 0..n.min(u_full.len()) {
        let byte_idx = i / 8;
        if byte_idx < q_bytes.len() && ((q_bytes[byte_idx] >> (i % 8)) & 1) == 1 {
            circ.x(&u_full[i]);
        }
    }
    while let Some(q) = u_full.pop() {
        circ.zero_and_free(q);
    }

    drop(b0);
    drop(b0_and_b1);
    circ.pop_section(&prev);
}

// ===========================================================================
// Apply / Bezout unpack — generic over R: ModRed.
// Curve-generic copy of bezout_unpack::apply_bitvector_quantum_secp256k1{,_inv}.
// ===========================================================================

/// Quantum Bezout reconstruction. On input `(x_reg, y_reg) = (z, 0)` and
/// garbage from `to_bitvector(q, x_orig)`, produces `(0, z * x_orig mod q)`.
///
/// Copy of `bezout_unpack::apply_bitvector_quantum_secp256k1`: `n = mr.n()`,
/// the per-iter `mod_double` + `controlled_mod_add` routed through `mr`.
pub fn apply_bitvector_quantum<R: ModRed>(
    circ: &mut Circuit,
    garbage: &[QReg],
    x_reg: &[QReg],
    y_reg: &[QReg],
    mr: &R,
) {
    use crate::point_add::trailmix_port::arith::schrottenloher::gcd_compress5::{
        compress_5iter_refs, compress_5iter_reverse_refs,
    };

    let n = mr.n();
    assert_eq!(x_reg.len(), n + 1, "x_reg must be n+1 bits");
    assert_eq!(y_reg.len(), n + 1, "y_reg must be n+1 bits");
    let iters = expected_iterations(n);
    let garbage_len = iters / DIALOG_M * DIALOG_PACK;
    assert!(
        garbage.len() >= garbage_len,
        "garbage must have at least {garbage_len} bits"
    );

    let prev = circ.push_section("apply_bv");
    let b0 = circ.alloc_qreg("apply.b0");
    let b0_and_b1 = circ.alloc_qreg("apply.b0_and_b1");

    let mut win_anc: Option<(QReg, QReg)> = None;

    for i in (0..iters).rev() {
        let pack_off = DIALOG_PACK * (i / DIALOG_M);
        let slot = i % DIALOG_M;

        if slot == DIALOG_M - 1 {
            let a0 = circ.alloc_qreg("apply_win0");
            let a1 = circ.alloc_qreg("apply_win1");
            {
                let w: [&QReg; 10] = std::array::from_fn(|k| match k {
                    8 => &a0,
                    9 => &a1,
                    _ => &garbage[pack_off + k],
                });
                let dirty: [&QReg; 8] = std::array::from_fn(|j| &x_reg[j]);
                compress_5iter_reverse_refs(circ, &w, &dirty);
            }
            win_anc = Some((a0, a1));
        }

        {
            let pair = win_anc.as_ref().expect("window open during apply_bv");
            let (a0, a1) = (&pair.0, &pair.1);
            let w: [&QReg; 10] = std::array::from_fn(|k| match k {
                8 => a0,
                9 => a1,
                _ => &garbage[pack_off + k],
            });
            circ.swap(&b0, w[2 * slot]);
            circ.swap(&b0_and_b1, w[2 * slot + 1]);
        }

        // y := 2y mod q ; y += b0 * x mod q ; cswap(b0_and_b1, x, y).
        mr.mod_double(circ, y_reg);
        mr.controlled_mod_add(circ, &b0, x_reg, y_reg);
        for j in 0..=n {
            circ.cswap(&b0_and_b1, &x_reg[j], &y_reg[j]);
        }

        {
            let pair = win_anc.as_ref().expect("window open during apply_bv");
            let (a0, a1) = (&pair.0, &pair.1);
            let w: [&QReg; 10] = std::array::from_fn(|k| match k {
                8 => a0,
                9 => a1,
                _ => &garbage[pack_off + k],
            });
            circ.swap(&b0, w[2 * slot]);
            circ.swap(&b0_and_b1, w[2 * slot + 1]);
        }

        if slot == 0 {
            let (a0, a1) = win_anc.take().expect("window open during apply_bv");
            {
                let w: [&QReg; 10] = std::array::from_fn(|k| match k {
                    8 => &a0,
                    9 => &a1,
                    _ => &garbage[pack_off + k],
                });
                let dirty: [&QReg; 8] = std::array::from_fn(|j| &x_reg[j]);
                compress_5iter_refs(circ, &w, &dirty);
            }
            circ.zero_and_free(a1);
            circ.zero_and_free(a0);
        }
    }
    debug_assert!(win_anc.is_none(), "all apply_bv windows closed");

    // Strict-dealloc settle.
    circ.cx(&b0, &x_reg[0]);
    circ.cx(&b0_and_b1, &x_reg[0]);
    drop(b0);
    drop(b0_and_b1);
    circ.pop_section(&prev);
}

/// Forward-reading Algorithm 3 in the INVERSE-OPS direction (the division
/// direction). For input `(u, v) = (0, y)` and garbage from
/// `to_bitvector(q, x_orig)`, produces `(y * x_orig^-1, 0)`.
///
/// Copy of `bezout_unpack::apply_bitvector_quantum_secp256k1_inv`:
/// `controlled_mod_sub` and `mod_halve` routed through `mr`.
pub fn apply_bitvector_quantum_inv<R: ModRed>(
    circ: &mut Circuit,
    garbage: &[QReg],
    x_reg: &[QReg],
    y_reg: &[QReg],
    mr: &R,
) {
    use crate::point_add::trailmix_port::arith::schrottenloher::gcd_compress5::{
        compress_5iter_refs, compress_5iter_reverse_refs,
    };

    let n = mr.n();
    assert_eq!(x_reg.len(), n + 1, "x_reg must be n+1 bits");
    assert_eq!(y_reg.len(), n + 1, "y_reg must be n+1 bits");
    let iters = expected_iterations(n);
    let garbage_len = iters / DIALOG_M * DIALOG_PACK;
    assert!(garbage.len() >= garbage_len);

    let prev = circ.push_section("apply_bv_inv");
    let b0 = circ.alloc_qreg("apply_inv.b0");
    let b0_and_b1 = circ.alloc_qreg("apply_inv.b0_and_b1");

    let mut win_anc: Option<(QReg, QReg)> = None;

    for i in 0..iters {
        let pack_off = DIALOG_PACK * (i / DIALOG_M);
        let slot = i % DIALOG_M;

        if slot == 0 {
            let a0 = circ.alloc_qreg("apply_inv_win0");
            let a1 = circ.alloc_qreg("apply_inv_win1");
            {
                let w: [&QReg; 10] = std::array::from_fn(|k| match k {
                    8 => &a0,
                    9 => &a1,
                    _ => &garbage[pack_off + k],
                });
                let dirty: [&QReg; 8] = std::array::from_fn(|j| &x_reg[j]);
                compress_5iter_reverse_refs(circ, &w, &dirty);
            }
            win_anc = Some((a0, a1));
        }

        {
            let pair = win_anc.as_ref().expect("window open during apply_bv_inv");
            let (a0, a1) = (&pair.0, &pair.1);
            let w: [&QReg; 10] = std::array::from_fn(|k| match k {
                8 => a0,
                9 => a1,
                _ => &garbage[pack_off + k],
            });
            circ.swap(&b0, w[2 * slot]);
            circ.swap(&b0_and_b1, w[2 * slot + 1]);
        }

        // if b0_and_b1: cswap x, y.
        for j in 0..=n {
            circ.cswap(&b0_and_b1, &x_reg[j], &y_reg[j]);
        }

        // y -= b0 * x mod q ; y *= 2^-1 mod q.
        mr.controlled_mod_sub(circ, &b0, x_reg, y_reg);
        mr.mod_halve(circ, y_reg);

        {
            let pair = win_anc.as_ref().expect("window open during apply_bv_inv");
            let (a0, a1) = (&pair.0, &pair.1);
            let w: [&QReg; 10] = std::array::from_fn(|k| match k {
                8 => a0,
                9 => a1,
                _ => &garbage[pack_off + k],
            });
            circ.swap(&b0, w[2 * slot]);
            circ.swap(&b0_and_b1, w[2 * slot + 1]);
        }

        if slot == DIALOG_M - 1 {
            let (a0, a1) = win_anc.take().expect("window open during apply_bv_inv");
            {
                let w: [&QReg; 10] = std::array::from_fn(|k| match k {
                    8 => &a0,
                    9 => &a1,
                    _ => &garbage[pack_off + k],
                });
                let dirty: [&QReg; 8] = std::array::from_fn(|j| &x_reg[j]);
                compress_5iter_refs(circ, &w, &dirty);
            }
            circ.zero_and_free(a1);
            circ.zero_and_free(a0);
        }
    }

    // Strict-dealloc settle.
    circ.cx(&b0, &x_reg[0]);
    circ.cx(&b0_and_b1, &x_reg[0]);
    drop(b0);
    drop(b0_and_b1);
    circ.pop_section(&prev);
}

// ===========================================================================
// IPModMul — generic over R: ModRed.
// Curve-generic copy of mod_mul_eea::{mod_mul_in_place_eea_secp256k1, _reverse,
// clear_zeroed_drift_reg_secp256k1}.
// ===========================================================================

/// Clear a register holding the zeroed apply-bitvector output (value `≡ 0
/// (mod q)` whose representation is exactly `0` or `q`). X each bit where
/// `q` has a 1 (0 Toffoli).
///
/// Copy of `mod_mul_eea::clear_zeroed_drift_reg_secp256k1`, parameterized.
pub fn clear_zeroed_drift_reg(circ: &mut Circuit, reg: &[QReg], n: usize, q_bytes: &[u8]) {
    for i in 0..n {
        let byte = i / 8;
        if byte < q_bytes.len() && (q_bytes[byte] >> (i % 8)) & 1 == 1 {
            circ.x(&reg[i]);
        }
    }
}

/// In-place modular multiplication: `(x_full, y_full) -> (x_full, x*y mod q)`.
///
/// Copy of `mod_mul_eea::mod_mul_in_place_eea_secp256k1`, calling the
/// parameterized `forward_gcd_pack` / apply / clear with `mr.n()`,
/// `mr.q_bytes()`, and `mr`.
pub fn mod_mul_in_place_eea<R: ModRed>(
    circ: &mut Circuit,
    x_full: &mut Vec<QReg>,
    y_full: &[QReg],
    mr: &R,
) {
    let n = mr.n();
    assert!(x_full.len() >= n, "x_full must be at least n = {n} bits");
    assert_eq!(y_full.len(), n + 1, "y_full must be n+1 bits");

    let prev = circ.push_section("mod_mul_eea");
    let original_len = x_full.len();

    // Step 1: x_full -> garbage (allocated incrementally inside).
    let mut garbage = forward_gcd_pack_quantum(circ, x_full, n, mr.q_bytes());

    // Allocate the apply-bv scratchpad AFTER forward_gcd.
    let tmp_full: Vec<QReg> = (0..=n)
        .map(|i| circ.alloc_qreg(&format!("mmul_tmp[{i}]")))
        .collect();

    // Step 2: apply garbage. After: y_full = 0 (mod q), tmp_full = y*x_orig.
    apply_bitvector_quantum(circ, &garbage, y_full, &tmp_full, mr);

    // Step 3: swap y_full <-> tmp_full.
    let s_sw = circ.push_section("eea_swap");
    for j in 0..=n {
        circ.swap(&y_full[j], &tmp_full[j]);
    }
    circ.pop_section(&s_sw);

    // tmp_full now holds the zeroed register (≡ 0 mod q). The ModRed impl
    // clears it to |0> per its representation (PseudoMersenne: X-clear q;
    // GenericPrime: already canonical 0, no-op).
    mr.clear_zeroed(circ, &tmp_full);
    for q in tmp_full {
        circ.zero_and_free(q);
    }

    // Step 4: reverse the GCD pack, restoring x_full from garbage.
    forward_gcd_pack_quantum_reverse(circ, x_full, &mut garbage, n, mr.q_bytes());
    assert!(garbage.is_empty(), "reverse should drain garbage tape");

    while x_full.len() < original_len {
        x_full.push(circ.alloc_qreg("x_pad_restore"));
    }

    circ.pop_section(&prev);
}

/// Inverse direction: `(x_full, y_full) -> (x_full, y * x^-1 mod q)`.
///
/// Copy of `mod_mul_eea::mod_mul_in_place_eea_secp256k1_reverse`, routed
/// through `mr`.
pub fn mod_mul_in_place_eea_reverse<R: ModRed>(
    circ: &mut Circuit,
    x_full: &mut Vec<QReg>,
    y_full: &[QReg],
    mr: &R,
) {
    let n = mr.n();
    assert!(x_full.len() >= n);
    assert_eq!(y_full.len(), n + 1);

    let prev = circ.push_section("mod_mul_eea_rev");
    let original_len = x_full.len();

    // Step 1 of forward inverted = forward GCD on x.
    let mut garbage = forward_gcd_pack_quantum(circ, x_full, n, mr.q_bytes());

    // Allocate scratchpad AFTER forward_gcd.
    let tmp_full: Vec<QReg> = (0..=n)
        .map(|i| circ.alloc_qreg(&format!("mmul_rev_tmp[{i}]")))
        .collect();

    // Step 3 of forward inverted = swap y_full <-> tmp_full.
    let s_sw = circ.push_section("eea_rev_swap");
    for j in 0..=n {
        circ.swap(&y_full[j], &tmp_full[j]);
    }
    circ.pop_section(&s_sw);

    // Step 2: forward-reading Algorithm 3 (inverse ops). Produces
    // y * x_orig^-1; leaves tmp_full at canonical 0.
    apply_bitvector_quantum_inv(circ, &garbage, y_full, &tmp_full, mr);

    for q in tmp_full {
        circ.zero_and_free(q);
    }

    // Step 4 of forward inverted = reverse forward GCD = restore x.
    forward_gcd_pack_quantum_reverse(circ, x_full, &mut garbage, n, mr.q_bytes());
    assert!(garbage.is_empty(), "reverse should drain garbage tape");

    while x_full.len() < original_len {
        x_full.push(circ.alloc_qreg("x_pad_restore"));
    }

    circ.pop_section(&prev);
}
