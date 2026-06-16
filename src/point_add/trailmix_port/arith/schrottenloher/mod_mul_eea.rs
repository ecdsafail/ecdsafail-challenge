//! Alg 4: in-place modular multiplication via the EEA dialog.
//!
//!   Convert x -> garbage g via `gcd_pack` (consumes `x_reg`).
//!   Map y -> y*x mod q via `bezout_unpack` on (y, tmp=0) controlled by g.
//!   Convert g -> x via inverse `gcd_pack` (restores `x_reg`).
//!
//! Net: |x, y> -> |x, y*x mod q>.
//!
//! Direct port of Schrottenloher's `IPModMul` (gcd.py:438).

use crate::point_add::trailmix_port::arith::schrottenloher::bezout_unpack::{
    apply_bitvector_quantum_secp256k1_inv_m, apply_bitvector_quantum_secp256k1_m,
};
use crate::point_add::trailmix_port::arith::schrottenloher::gcd_jump::{
    apply_bitvector_jump_packed_inv_secp256k1, apply_bitvector_jump_packed_secp256k1,
    apply_bitvector_jump_quantum_secp256k1, apply_bitvector_jump_quantum_secp256k1_inv,
    compress_dialog_jump, decompress_dialog_jump, forward_gcd_jump_quantum_secp256k1,
    forward_gcd_jump_quantum_secp256k1_reverse,
};
use crate::point_add::trailmix_port::arith::schrottenloher::gcd_pack::{
    expected_iterations, forward_gcd_pack_quantum_secp256k1_m,
    forward_gcd_pack_quantum_secp256k1_reverse_m, u_padding, DIALOG_M,
};
use crate::point_add::trailmix_port::circuit::{Circuit, QReg};

/// In-place modular multiplication for secp256k1: `(x_full, y_full) ->
/// (x_full, x*y mod q)`.
///
/// Convention:
///   - `x_full`: caller-allocated `n + u_padding(n) = 293` bits. Pre:
///     low `n=256` bits hold `x`; high bits 0. Post: same (restored).
///   - `y_full`: caller-allocated `n + 1 = 257` bits. Pre: low `n`
///     bits hold `y`; bit 256 is 0. Post: low 256 bits hold `x*y
///     mod q` (representation drift up to a few multiples of q);
///     bit 256 = 0.
///   - `tmp_full`: caller-allocated `n + 1 = 257` bits, pre 0. Post:
///     contents are `≡ 0 (mod q)` but possibly = k*q for small k
///     (pseudo-Mersenne drift). Caller is responsible for any final
///     reduction or for accepting the drift through subsequent
///     operations. In the IPModMul-style pointadd driver, `tmp_full`
///     gets used as the next mul's accumulator slot and drift cancels
///     structurally.
pub fn mod_mul_in_place_eea_secp256k1(circ: &mut Circuit, x_full: &mut Vec<QReg>, y_full: &[QReg]) {
    mod_mul_in_place_eea_secp256k1_m(circ, x_full, y_full, DIALOG_M, 0);
}

/// Window-size- and vent-parameterized in-place modular multiply.
/// `dialog_m = 5` is the base-3-packed (SP1-validated, 1173q) config;
/// `dialog_m = 1` is the RAW higher-qubit config (wider tape, no
/// `compress_5iter`). `vents` is the apply_bv-peak measurement-vent budget
/// (each -1 Toffoli / +1 peak; reused across the sequential Horner ops).
pub fn mod_mul_in_place_eea_secp256k1_m(
    circ: &mut Circuit,
    x_full: &mut Vec<QReg>,
    y_full: &[QReg],
    dialog_m: usize,
    vents: usize,
) {
    let n = 256usize;
    let pad = u_padding(n);
    let iters = expected_iterations(n);
    let garbage_len = crate::point_add::trailmix_port::arith::schrottenloher::gcd_pack::expected_garbage_m(n, dialog_m);

    assert!(x_full.len() >= n, "x_full must be at least n = {n} bits");
    assert_eq!(y_full.len(), n + 1, "y_full must be n+1 = 257 bits");

    let prev = circ.push_section("mod_mul_eea");
    let original_len = x_full.len();

    // Step 1: x_full -> garbage (allocated incrementally inside).
    // x_full's Vec is shrunk to len 1 by the forward pass. tmp is NOT
    // allocated yet: during forward_gcd (the peak section) we do not
    // want 257 idle scratch qubits. This matches the author's IPModMul,
    // where ToBitVector runs before the tmp register exists.
    let mut garbage = forward_gcd_pack_quantum_secp256k1_m(circ, x_full, dialog_m);

    // Allocate the apply-bv scratchpad AFTER forward_gcd. Peak now occurs
    // in apply_bv (garbage + y + tmp), not forward_gcd.
    let tmp_full: Vec<QReg> = (0..=n)
        .map(|i| circ.alloc_qreg(&format!("mmul_tmp[{i}]")))
        .collect();

    // Step 2: apply garbage. After: y_full = 0 (mod q), tmp_full =
    // y * x_orig (mod q).
    apply_bitvector_quantum_secp256k1_m(circ, &garbage, y_full, &tmp_full, dialog_m, vents);

    // Step 3: swap y_full <-> tmp_full so the product ends up in
    // y_full and tmp_full returns to (drift-)0.
    let s_sw = circ.push_section("eea_swap");
    for j in 0..=n {
        circ.swap(&y_full[j], &tmp_full[j]);
    }
    circ.pop_section(&s_sw);

    // tmp_full now holds the zeroed register: value is ≡ 0 (mod q) with
    // representation in [0, 2^256), i.e. exactly 0 or q. Clear it to |0>
    // by XOR-ing the constant q (X gates, 0 Toffoli) and free, BEFORE
    // reverse_gcd, so tmp is not live during the GCD passes.
    clear_zeroed_drift_reg_secp256k1(circ, &tmp_full);
    for q in tmp_full {
        circ.zero_and_free(q);
    }

    // Step 4: reverse the GCD pack, restoring x_full from garbage.
    // Both x_full (regrows) and garbage (drained pack-by-pack as each
    // becomes |0>) are mutated.
    forward_gcd_pack_quantum_secp256k1_reverse_m(circ, x_full, &mut garbage, dialog_m);
    assert!(garbage.is_empty(), "reverse should drain garbage tape");

    // The reverse grew x_full back to n=256; restore to original_len
    // (the caller's slot for the anc bit at position n).
    while x_full.len() < original_len {
        x_full.push(circ.alloc_qreg("x_pad_restore"));
    }

    let _ = (pad, iters, garbage_len);
    circ.pop_section(&prev);
}

/// Clear a register that holds the zeroed apply-bitvector output: a value
/// `≡ 0 (mod q)` whose representation in `[0, 2^256)` is exactly `0` or
/// `q`. Apply X to each bit where `q` has a 1; if the register held `q`
/// this maps it to 0, and if it held 0 this would map it to q (caught by
/// the subsequent `zero_and_free` if the determinism assumption breaks).
fn clear_zeroed_drift_reg_secp256k1(circ: &mut Circuit, reg: &[QReg]) {
    let n = 256usize;
    let q: crate::point_add::trailmix_port::compat::num_bigint::BigUint =
        (crate::point_add::trailmix_port::compat::num_bigint::BigUint::from(1u32) << 256u32) - crate::point_add::trailmix_port::compat::num_bigint::BigUint::from(super::F_SECP256K1);
    let q_bytes = q.to_bytes_le();
    for i in 0..n {
        let byte = i / 8;
        if byte < q_bytes.len() && (q_bytes[byte] >> (i % 8)) & 1 == 1 {
            circ.x(&reg[i]);
        }
    }
}

/// Inverse direction: `(x_full, y_full) -> (x_full, y * x^-1 mod q)`.
///
/// Per Schrottenloher 2026 Sec 3: "Notice that the inverse circuit
/// will multiply by x^-1 instead. So, using the same circuit, we can
/// implement both in-place multiplication steps found in Algorithm 1."
/// Alg 1 line 6 (`x2, y2 ← x2, y2 * x2^-1`) uses this reverse
/// direction; line 11 (`x2, y2 ← x2, y2 * x2`) uses the forward.
///
/// Implementation: gate-for-gate inverse of `mod_mul_in_place_eea`.
/// Same caller conventions on `x_full` / `y_full` / `tmp_full`.
pub fn mod_mul_in_place_eea_secp256k1_reverse(
    circ: &mut Circuit,
    x_full: &mut Vec<QReg>,
    y_full: &[QReg],
) {
    mod_mul_in_place_eea_secp256k1_reverse_m(circ, x_full, y_full, DIALOG_M, 0);
}

/// Window-size- and vent-parameterized in-place modular division
/// (`y -> y * x^-1`). `dialog_m = 1` = RAW higher-qubit; `dialog_m = 5` =
/// base-3-packed (SP1-validated). `vents` = apply_bv-peak vent budget.
pub fn mod_mul_in_place_eea_secp256k1_reverse_m(
    circ: &mut Circuit,
    x_full: &mut Vec<QReg>,
    y_full: &[QReg],
    dialog_m: usize,
    vents: usize,
) {
    let n = 256usize;
    let pad = u_padding(n);
    let iters = expected_iterations(n);

    assert!(x_full.len() >= n);
    assert_eq!(y_full.len(), n + 1);

    let prev = circ.push_section("mod_mul_eea_rev");
    let _ = iters;
    let _ = pad;
    let original_len = x_full.len();

    // Step 1 of forward inverted = forward GCD on x (driving x -> 0).
    let mut garbage = forward_gcd_pack_quantum_secp256k1_m(circ, x_full, dialog_m);

    // Allocate scratchpad AFTER forward_gcd (peak avoidance, see forward).
    let tmp_full: Vec<QReg> = (0..=n)
        .map(|i| circ.alloc_qreg(&format!("mmul_rev_tmp[{i}]")))
        .collect();

    // Step 3 of forward inverted = swap y_full <-> tmp_full.
    let s_sw = circ.push_section("eea_rev_swap");
    for j in 0..=n {
        circ.swap(&y_full[j], &tmp_full[j]);
    }
    circ.pop_section(&s_sw);

    // Step 2: forward-reading Algorithm 3 with approximate mod-halve and
    // controlled mod-sub. Produces y * x_orig^-1; leaves tmp_full at the
    // zeroed (≡0 mod q) representation.
    apply_bitvector_quantum_secp256k1_inv_m(circ, &garbage, y_full, &tmp_full, dialog_m, vents);

    // tmp_full holds the zeroed register. Unlike the forward direction
    // (which leaves the non-canonical representative q), the inverse
    // apply leaves it at canonical 0, so free directly.
    for q in tmp_full {
        circ.zero_and_free(q);
    }

    // Step 4 of forward inverted = reverse forward GCD = restore x.
    // Reverse drains garbage tape pack-by-pack as each completes.
    forward_gcd_pack_quantum_secp256k1_reverse_m(circ, x_full, &mut garbage, dialog_m);
    assert!(garbage.is_empty(), "reverse should drain garbage tape");

    // Restore caller's high pad bits (alloc'd 0, untouched by gcd_pack).
    while x_full.len() < original_len {
        x_full.push(circ.alloc_qreg("x_pad_restore"));
    }

    let _ = (pad, iters);
    circ.pop_section(&prev);
}

/// Jump-GCD in-place modular multiply: `(x_full, y_full) -> (x_full, y*x mod q)`.
/// Mirror of `mod_mul_in_place_eea_secp256k1_m` using the jump-before-swap
/// GCD/apply (fewer steps -> fewer per-step adders).
pub fn mod_mul_in_place_jump_secp256k1(
    circ: &mut Circuit,
    x_full: &mut Vec<QReg>,
    y_full: &[QReg],
    jump: usize,
    vents: usize,
    coupled: bool,
    packed: bool,
) {
    let n = 256usize;
    assert!(x_full.len() >= n, "x_full must be at least n = {n} bits");
    assert_eq!(y_full.len(), n + 1, "y_full must be n+1 = 257 bits");
    let prev = circ.push_section("mod_mul_jump");
    let original_len = x_full.len();

    let mut garbage = forward_gcd_jump_quantum_secp256k1(circ, x_full, jump);
    // base-5 PACK the dialog (off-peak) so the apply peak drops ~180q.
    if packed {
        compress_dialog_jump(circ, &mut garbage);
    }
    let tmp_full: Vec<QReg> = (0..=n)
        .map(|i| circ.alloc_qreg(&format!("jmul_tmp[{i}]")))
        .collect();
    // y -> 0, tmp -> y*x mod q.
    if packed {
        apply_bitvector_jump_packed_secp256k1(circ, &garbage, y_full, &tmp_full, vents, coupled);
    } else {
        apply_bitvector_jump_quantum_secp256k1(circ, &garbage, y_full, &tmp_full, jump);
    }
    // product into y_full; tmp returns to the (drift-)0 representation.
    let s_sw = circ.push_section("jmul_swap");
    for j in 0..=n {
        circ.swap(&y_full[j], &tmp_full[j]);
    }
    circ.pop_section(&s_sw);
    clear_zeroed_drift_reg_secp256k1(circ, &tmp_full);
    for q in tmp_full {
        circ.zero_and_free(q);
    }
    if packed {
        decompress_dialog_jump(circ, &mut garbage);
    }
    forward_gcd_jump_quantum_secp256k1_reverse(circ, x_full, &mut garbage, jump);
    assert!(garbage.is_empty(), "reverse should drain garbage tape");
    while x_full.len() < original_len {
        x_full.push(circ.alloc_qreg("x_pad_restore"));
    }
    circ.pop_section(&prev);
}

/// Jump-GCD in-place modular division: `(x_full, y_full) -> (x_full, y*x^-1 mod q)`.
/// Mirror of `mod_mul_in_place_eea_secp256k1_reverse_m` using the jump `apply_inv`.
pub fn mod_mul_in_place_jump_reverse_secp256k1(
    circ: &mut Circuit,
    x_full: &mut Vec<QReg>,
    y_full: &[QReg],
    jump: usize,
    vents: usize,
    coupled: bool,
    packed: bool,
) {
    let n = 256usize;
    assert!(x_full.len() >= n);
    assert_eq!(y_full.len(), n + 1);
    let prev = circ.push_section("mod_mul_jump_rev");
    let original_len = x_full.len();

    let mut garbage = forward_gcd_jump_quantum_secp256k1(circ, x_full, jump);
    if packed {
        compress_dialog_jump(circ, &mut garbage);
    }
    let tmp_full: Vec<QReg> = (0..=n)
        .map(|i| circ.alloc_qreg(&format!("jmul_rev_tmp[{i}]")))
        .collect();
    let s_sw = circ.push_section("jmul_rev_swap");
    for j in 0..=n {
        circ.swap(&y_full[j], &tmp_full[j]);
    }
    circ.pop_section(&s_sw);
    // forward-reading inverse apply -> y*x^-1; tmp left at canonical 0.
    if packed {
        apply_bitvector_jump_packed_inv_secp256k1(
            circ, &garbage, y_full, &tmp_full, vents, coupled,
        );
    } else {
        apply_bitvector_jump_quantum_secp256k1_inv(circ, &garbage, y_full, &tmp_full, jump);
    }
    for q in tmp_full {
        circ.zero_and_free(q);
    }
    if packed {
        decompress_dialog_jump(circ, &mut garbage);
    }
    forward_gcd_jump_quantum_secp256k1_reverse(circ, x_full, &mut garbage, jump);
    assert!(garbage.is_empty(), "reverse should drain garbage tape");
    while x_full.len() < original_len {
        x_full.push(circ.alloc_qreg("x_pad_restore"));
    }
    circ.pop_section(&prev);
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::point_add::trailmix_port::compat::num_bigint::BigUint;
    use crate::point_add::trailmix_port::compat::num_traits::Zero;

    // The previous mul_then_inverse_restores test required forward
    // mod_mul to be the GATE-INVERSE of reverse mod_mul. With the new
    // reverse direction using exact apply_inv (mod_halve_pm_general +
    // controlled_mod_sub_mbu) and forward still using approximate
    // pseudo-Mersenne primitives, the two don't compose as
    // gate-inverses — forward output can have drift in [0, 2^256)
    // which the exact reverse rejects via tracker MBU identity checks.
    //
    // Both directions ARE correct mod q on their respective domain
    // (forward: any input → product mod q; reverse: any input →
    // quotient mod q). They just don't compose gate-by-gate. The
    // EC-add Alg 1 calls them on disjoint logical slots, which works.

    #[test]
    fn mod_mul_in_place_eea_secp256k1_value_correct() {
        use rand::{thread_rng, Rng};

        let n = 256usize;
        let pad = u_padding(n);
        let total = n + pad;
        let q: BigUint = (BigUint::from(1u32) << 256u32) - BigUint::from(super::super::F_SECP256K1);

        let mut rng = thread_rng();
        let mut shot_data: Vec<(BigUint, BigUint, BigUint)> = Vec::new(); // (x, y, want)
        while shot_data.len() < 64 {
            let raw_x: [u8; 32] = rng.gen();
            let raw_y: [u8; 32] = rng.gen();
            let mut x = BigUint::from_bytes_le(&raw_x) % &q;
            if !x.bit(0) {
                x += BigUint::from(1u32);
                x %= &q;
                if !x.bit(0) {
                    continue;
                }
            }
            // Skip GCD budget tail (rare).
            if super::super::gcd_pack::to_bitvector_classical(q.clone(), x.clone(), n).is_none() {
                continue;
            }
            let y = BigUint::from_bytes_le(&raw_y) % &q;
            let want = (&x * &y) % &q;
            shot_data.push((x, y, want));
        }

        let mut circ = crate::point_add::trailmix_port::circuit::Circuit::new();
        circ.set_max_qubit_peak(1180);
        let mut x_full = circ.alloc_qreg_bits("x_full", total);
        let y_full = circ.alloc_qreg_bits("y_full", n + 1);

        for (shot, (x, y, _want)) in shot_data.iter().enumerate() {
            let mut xbuf = [0u8; 32];
            for (i, b) in x.to_bytes_le().iter().take(32).enumerate() {
                xbuf[i] = *b;
            }
            circ.sim_load_reg_bytes_shot(&x_full[..n], &xbuf, shot);

            let mut ybuf = [0u8; 32];
            for (i, b) in y.to_bytes_le().iter().take(32).enumerate() {
                ybuf[i] = *b;
            }
            circ.sim_load_reg_bytes_shot(&y_full[..n], &ybuf, shot);
        }

        let ops0 = circ.ops.len();
        let ccx0 = circ.ccx_emitted;
        let ccz0 = circ.ccz_emitted;
        mod_mul_in_place_eea_secp256k1(&mut circ, &mut x_full, &y_full);
        let mul_ops = circ.ops.len() - ops0;
        let mul_tof = (circ.ccx_emitted - ccx0) + (circ.ccz_emitted - ccz0);
        let peak = circ.peak_qubits;
        eprintln!(
            "  cost(mod_mul_in_place_eea_secp256k1, n=256): ops={mul_ops} tof={mul_tof} peak+={peak}"
        );

        let mut outs: Vec<crate::point_add::trailmix_port::circuit::QReg> = Vec::new();
        outs.extend(x_full);
        outs.extend(y_full);
        let (sim, detached) = circ.destroy_sim(outs);

        for (shot, (x_in, _y_in, want)) in shot_data.iter().enumerate() {
            // x_full restored: low n bits = x_in (mod q congruence ok).
            let mut got_x = BigUint::zero();
            for i in 0..n {
                if sim.read_bit_shot(&detached[i], shot) == 1 {
                    got_x |= BigUint::from(1u32) << i;
                }
            }
            let got_x_mod = &got_x % &q;
            assert_eq!(got_x_mod, *x_in, "shot {shot}: x not restored");

            // y_full holds x*y mod q.
            let y_off = total;
            let mut got_y = BigUint::zero();
            for i in 0..n {
                if sim.read_bit_shot(&detached[y_off + i], shot) == 1 {
                    got_y |= BigUint::from(1u32) << i;
                }
            }
            let got_y_mod = &got_y % &q;
            assert_eq!(got_y_mod, *want, "shot {shot}: y != x*y mod q");
        }
    }
}
