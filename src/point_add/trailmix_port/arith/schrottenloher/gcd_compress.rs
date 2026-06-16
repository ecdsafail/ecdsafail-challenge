//! Fig.1: in-place 5-Toffoli compression mapping three successive
//! `(b0, b0&b1)` pairs (6 bits, 27 valid values out of 64) into 5
//! useful bits + 1 reusable zero.
//!
//! Direct port of Schrottenloher's `Compressor` (`point_add/compressor.py:138`).
//! Synthesized via SAT in the paper; we replicate the exact gate
//! sequence here.
//!
//! Layout: input register `r` of width 6 (little-endian). Pairs are
//! `(r[0], r[1])`, `(r[2], r[3])`, `(r[4], r[5])`. Each pair is one of
//! `(0,0)`, `(1,0)`, `(1,1)` per the GCD invariant. After the circuit,
//! the 5-bit compressed value lives in `r[0..5]` and `r[5]` is
//! guaranteed to be 0 on valid inputs.

use crate::point_add::trailmix_port::circuit::{Circuit, QReg};

/// Refs variant of the Fig.1 compressor: takes the 6 qubits as `&QReg`
/// references in any storage layout. Useful when the 6th qubit is
/// allocated separately from the first 5 (as in the Swapper pattern,
/// where the 5 packed bits live in a shared garbage register and the
/// 6th is a per-call ancilla).
pub fn compress_3iter_refs(circ: &mut Circuit, r: &[&QReg; 6]) {
    // Direct port from Schrottenloher's Compressor circuit, point_add/compressor.py:151.
    circ.cx(r[1], r[0]);
    circ.cx(r[3], r[2]);
    circ.cx(r[5], r[4]);

    circ.cx(r[0], r[2]);
    circ.cx(r[5], r[3]);
    circ.x(r[4]);
    circ.ccx(r[1], r[3], r[5]);
    circ.cx(r[1], r[4]);
    circ.x(r[2]);
    circ.ccx(r[3], r[4], r[5]);
    circ.ccx(r[4], r[5], r[1]);
    circ.ccx(r[2], r[5], r[0]);
    circ.ccx(r[0], r[1], r[5]);
}

/// Refs variant of the Fig.1 inverse compressor: restores the 6-bit
/// pair-encoded form from the 5-bit compressed bit-string in `r[0..5]`,
/// with `r[5] = 0` pre-call.
pub fn compress_3iter_reverse_refs(circ: &mut Circuit, r: &[&QReg; 6]) {
    circ.ccx(r[0], r[1], r[5]);
    circ.ccx(r[2], r[5], r[0]);
    circ.ccx(r[4], r[5], r[1]);
    circ.ccx(r[3], r[4], r[5]);
    circ.x(r[2]);
    circ.cx(r[1], r[4]);
    circ.ccx(r[1], r[3], r[5]);
    circ.x(r[4]);
    circ.cx(r[5], r[3]);
    circ.cx(r[0], r[2]);

    circ.cx(r[5], r[4]);
    circ.cx(r[3], r[2]);
    circ.cx(r[1], r[0]);
}

/// Schrottenloher's "Swapper" (compressor.py:181): in-place
/// swap of the 2-bit `bb` register with the `(2*slot, 2*slot+1)` slot
/// inside the compressed 5-bit pack `pack5`.
///
/// Allocates its own 1-bit ancilla internally and frees it before
/// returning. (Earlier API took `pack_anc` as a parameter, but
/// callers had to keep it alive across multi-iter loops where the
/// subsequent compute had no natural touch on it — triggering the
/// strict-dealloc retention check. Internal management makes the
/// ancilla per-call.)
///
/// Procedure:
///   1) Decompress `pack5+pack_anc` to a 6-bit pair-encoded layout.
///   2) Swap bb[0] with the lo bit of slot, bb[1] with the hi bit.
///   3) Recompress; `pack_anc` returns to |0> and is freed.
///
/// Self-inverse: a second call un-swaps. On valid inputs (the slot
/// currently holds `(0, 0)`), this lifts the (b0, b0&b1) pair into
/// the pack; calling again zeroes (b0, b0&b1) back out.
pub fn swapper(circ: &mut Circuit, slot: usize, bb: &[&QReg; 2], pack5: &[QReg]) {
    assert!(slot < 3, "slot must be 0, 1, or 2");
    assert_eq!(pack5.len(), 5);
    let pack_anc = circ.alloc_qreg("swp_anc");
    let pack6: [&QReg; 6] = [
        &pack5[0], &pack5[1], &pack5[2], &pack5[3], &pack5[4], &pack_anc,
    ];
    compress_3iter_reverse_refs(circ, &pack6);
    circ.swap(bb[0], pack6[2 * slot]);
    circ.swap(bb[1], pack6[2 * slot + 1]);
    compress_3iter_refs(circ, &pack6);
    circ.zero_and_free(pack_anc);
}
