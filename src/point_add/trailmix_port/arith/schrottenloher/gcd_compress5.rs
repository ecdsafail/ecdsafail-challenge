//! M=5 dialog compressor: in-place radix-3 -> binary encoder mapping
//! five successive `(b0, b0&b1)` pairs (10 bits, 243 valid of 1024)
//! into an 8-bit value + 2 reusable zeros.
//!
//! Each pair is one of `(0,0)`, `(1,0)`, `(1,1)` per the GCD invariant
//! (the `(0,1)` combo is impossible since `b0&b1=1` forces `b0=1`), so
//! the window is base-3: digit `d_i = b0 + (b0&b1) in {0,1,2}`. The
//! packed value is the literal base-3 number `sum_i d_i * 3^i in [0,243)`.
//!
//! Density 8/5 = 1.600 bits/pair vs the M=3 compressor's 5/3 = 1.667,
//! so the GCD/Bezout garbage tape shrinks ~3.3% (670 -> 648 bits at
//! n=256, 402 iters), dropping the `apply_bv` peak ~22 qubits.
//!
//! Construction (low-to-high radix merge, fully in-place):
//!   1. pair -> digit: `cx(w[2i+1], w[2i])` makes `w[2i]=(d_i==1)`,
//!      `w[2i+1]=(d_i==2)`, whose 2-bit value is exactly `d_i`.
//!   2. for i = 1..=4: accumulate `e += 3^i * d_i` (two controlled
//!      constant-adds, on the mutually-exclusive lo/hi digit bits),
//!      then clear the digit bits by `d_i = e // 3^i` recovered from
//!      two constant comparisons (e is preserved, `e_prev` < 3^i so the
//!      quotient is exactly `d_i`). For i<4 the digit is first staged
//!      into 2 scratch bits so the add (which writes e's range,
//!      overlapping w[2i],w[2i+1]) cannot clobber its own controls;
//!      for i=4 the digit lives at w[8],w[9], above e's 8-bit range,
//!      so it controls the add directly.
//!
//! After the call: `w[0..8]` holds the base-3 value, `w[8]=w[9]=0`,
//! all scratch returned to |0>. The decompressor is the exact inverse.

use crate::point_add::trailmix_port::arith::gidney_const_adder::{
    compare_geq_const_gidney_refs, controlled_add_const_gidney_refs,
};
use crate::point_add::trailmix_port::circuit::{Circuit, QReg};

/// 3^i for i in 0..=5.
const POW3: [u16; 6] = [1, 3, 9, 27, 81, 243];

/// Bit-width of the accumulator after merging digit i (value < 3^{i+1}).
/// i: 1->4 (<9), 2->5 (<27), 3->7 (<81), 4->8 (<243).
const ACC_WIDTH: [usize; 5] = [2, 4, 5, 7, 8];

/// Compress 5 successive `(b0, b0&b1)` pairs held in `w[0..10]` (little
/// endian pairs `(w[0],w[1]) .. (w[8],w[9])`) into the 8-bit base-3
/// value in `w[0..8]`. `w[8]` and `w[9]` are guaranteed 0 afterward on
/// valid input. Allocates and frees 2 scratch qubits internally.
/// `dirty`: >= 8 BORROWED qubits, disjoint from `w`, restored on exit
/// (the constant comparator vents its carries into them). Callers with
/// peak headroom (GCD pack) pass freshly-allocated |0> bits; the peak-
/// critical `apply_bv` borrows idle data-register bits so the compressor
/// adds ~0 clean peak.
pub fn compress_5iter_refs(circ: &mut Circuit, w: &[&QReg; 10], dirty: &[&QReg]) {
    let c5prev = circ.push_section("c5pack");
    // 1. pair -> digit encoding
    for i in 0..5 {
        circ.cx(w[2 * i + 1], w[2 * i]);
    }

    // staging scratch reused by merges i=1,2,3 (i=4 needs none)
    let s0 = circ.alloc_qreg("c5_stage0");
    let s1 = circ.alloc_qreg("c5_stage1");

    for i in 1..=3 {
        merge_digit(circ, w, i, Some((&s0, &s1)), dirty);
    }
    // s0,s1 are |0> after merge 3; free before merge 4 (which allocs).
    circ.zero_and_free(s1);
    circ.zero_and_free(s0);

    merge_digit(circ, w, 4, None, dirty);
    circ.pop_section(&c5prev);
}

/// Inverse of `compress_5iter_refs`. `dirty` per `compress_5iter_refs`.
pub fn compress_5iter_reverse_refs(circ: &mut Circuit, w: &[&QReg; 10], dirty: &[&QReg]) {
    let c5prev = circ.push_section("c5pack");
    unmerge_digit(circ, w, 4, None, dirty);

    let s0 = circ.alloc_qreg("c5_stage0");
    let s1 = circ.alloc_qreg("c5_stage1");
    for i in (1..=3).rev() {
        unmerge_digit(circ, w, i, Some((&s0, &s1)), dirty);
    }
    circ.zero_and_free(s1);
    circ.zero_and_free(s0);

    for i in 0..5 {
        circ.cx(w[2 * i + 1], w[2 * i]);
    }
    circ.pop_section(&c5prev);
}

/// Forward merge of digit `i` into the accumulator. `stage` is the
/// scratch pair for i<4; for i==4 it's None and the digit controls at
/// `w[8],w[9]` are used directly. `dirty` (>= wi borrowed bits, disjoint
/// from `w`) is the constant comparator's vented-carry scratch.
fn merge_digit(
    circ: &mut Circuit,
    w: &[&QReg; 10],
    i: usize,
    stage: Option<(&QReg, &QReg)>,
    dirty: &[&QReg],
) {
    let r = POW3[i] as u8; // 3^i
    let r2 = (2 * POW3[i]) as u8; // 2*3^i
    let wi = ACC_WIDTH[i];
    let e: Vec<&QReg> = (0..wi).map(|k| w[k]).collect();

    let (clo, chi): (&QReg, &QReg) = match stage {
        Some((s0, s1)) => {
            // move digit out of e's range into scratch
            circ.swap(w[2 * i], s0);
            circ.swap(w[2 * i + 1], s1);
            (s0, s1)
        }
        None => (w[8], w[9]),
    };

    // e += 3^i * d_i  (controls are mutually exclusive) via the cheap ~3n
    // borrowed-dirty Gidney controlled const-add (vs the dense CQ add_const).
    controlled_add_const_gidney_refs(circ, clo, &e, &[r], dirty);
    controlled_add_const_gidney_refs(circ, chi, &e, &[r2], dirty);

    // clear the digit: clo ^= (e>=r)^(e>=2r) = (d==1); chi ^= (e>=2r) = (d==2)
    // via the cheap ~3n vented-carry constant comparator (borrowed dirty).
    compare_geq_const_gidney_refs(circ, &e, &[r], clo, dirty);
    compare_geq_const_gidney_refs(circ, &e, &[r2], clo, dirty);
    compare_geq_const_gidney_refs(circ, &e, &[r2], chi, dirty);
}

/// Inverse of `merge_digit`.
fn unmerge_digit(
    circ: &mut Circuit,
    w: &[&QReg; 10],
    i: usize,
    stage: Option<(&QReg, &QReg)>,
    dirty: &[&QReg],
) {
    let r = POW3[i] as u8;
    let r2 = (2 * POW3[i]) as u8;
    let wi = ACC_WIDTH[i];
    let e: Vec<&QReg> = (0..wi).map(|k| w[k]).collect();

    let (clo, chi): (&QReg, &QReg) = match stage {
        Some((s0, s1)) => (s0, s1),
        None => (w[8], w[9]),
    };

    // un-clear: reverse the 3 XORs to restore clo=(d==1), chi=(d==2)
    compare_geq_const_gidney_refs(circ, &e, &[r2], chi, dirty);
    compare_geq_const_gidney_refs(circ, &e, &[r2], clo, dirty);
    compare_geq_const_gidney_refs(circ, &e, &[r], clo, dirty);

    // un-add: e -= chi*r2 then e -= clo*r, each via Gidney add of the
    // wi-bit two's complement 2^wi - k (popcount-independent, ~3n).
    let cr2 = ((1u16 << wi) - u16::from(r2)).to_le_bytes();
    let cr = ((1u16 << wi) - u16::from(r)).to_le_bytes();
    controlled_add_const_gidney_refs(circ, chi, &e, &cr2, dirty);
    controlled_add_const_gidney_refs(circ, clo, &e, &cr, dirty);

    if let Some((s0, s1)) = stage {
        // swap digit back from scratch into w[2i],w[2i+1]
        circ.swap(w[2 * i], s0);
        circ.swap(w[2 * i + 1], s1);
    }
}

/// M=5 analogue of `gcd_compress::swapper`: in-place swap of the 2-bit
/// `bb` register with the `(2*slot, 2*slot+1)` pair inside the compressed
/// 8-bit pack `pack8`. Allocates its own 2 scratch bits (the decompressed
/// window's high pair) and frees them before returning.
///
/// Procedure: decompress pack8 + 2 scratch to the 10-bit pair layout,
/// swap `bb` with the slot, recompress (scratch returns to |0>).
/// Self-inverse: a second call un-swaps. On valid input (the slot holds
/// `(0,0)`), this lifts the `(b0, b0&b1)` pair into the pack; calling
/// again zeroes it back out.
pub fn swapper5(circ: &mut Circuit, slot: usize, bb: &[&QReg; 2], pack8: &[QReg]) {
    assert!(slot < 5, "slot must be in 0..5");
    assert_eq!(pack8.len(), 8);
    let a0 = circ.alloc_qreg("swp5_a0");
    let a1 = circ.alloc_qreg("swp5_a1");
    // Fresh |0> comparator dirty scratch (8 bits, disjoint from the pack).
    // swapper5 is the per-iter path (debug/tests) where peak headroom exists;
    // the per-window apply_bv calls compress_5iter_*_refs directly with
    // BORROWED dirty to avoid the +8 peak.
    let dscratch: Vec<QReg> = (0..8).map(|_| circ.alloc_qreg("swp5_dirty")).collect();
    let dirty: Vec<&QReg> = dscratch.iter().collect();
    let w: [&QReg; 10] = [
        &pack8[0], &pack8[1], &pack8[2], &pack8[3], &pack8[4], &pack8[5], &pack8[6], &pack8[7],
        &a0, &a1,
    ];
    compress_5iter_reverse_refs(circ, &w, &dirty);
    circ.swap(bb[0], w[2 * slot]);
    circ.swap(bb[1], w[2 * slot + 1]);
    compress_5iter_refs(circ, &w, &dirty);
    drop(dirty);
    for q in dscratch.into_iter().rev() {
        circ.zero_and_free(q);
    }
    circ.zero_and_free(a1);
    circ.zero_and_free(a0);
}

/// Classical decompress: inverse of `compress_classical_5`. Returns the
/// 10-bit pair-encoded input, or `None` if `packed` is not in the image
/// (i.e. >= 243).
#[must_use]
pub fn uncompress_classical_5(packed: u8) -> Option<u16> {
    if packed >= 243 {
        return None;
    }
    // base-3 digits low-to-high, re-encoded as (b0, b0&b1) pairs:
    // d=0 -> (0,0), d=1 -> (1,0), d=2 -> (1,1).
    let mut v = u16::from(packed);
    let mut out = 0u16;
    for i in 0..5u32 {
        let d = v % 3;
        v /= 3;
        let (b0, b0b1) = match d {
            0 => (0u16, 0u16),
            1 => (1, 0),
            _ => (1, 1),
        };
        out |= b0 << (2 * i);
        out |= b0b1 << (2 * i + 1);
    }
    Some(out)
}

/// Classical reference: map a 10-bit input (5 pairs) to its base-3
/// value, or `None` if any pair is the invalid `(0,1)` encoding.
///
/// Input bit layout: `(b0_0, b0&b1_0, b0_1, b0&b1_1, ..., b0_4, b0&b1_4)`
/// = bits 0..10. Each pair must be in `{(0,0),(1,0),(1,1)}`.
#[must_use]
pub fn compress_classical_5(input: u16) -> Option<u8> {
    let bit = |k: u32| (input >> k) & 1;
    let mut acc: u16 = 0;
    for i in 0..5u32 {
        let b0 = bit(2 * i);
        let b0b1 = bit(2 * i + 1);
        if b0 == 0 && b0b1 == 1 {
            return None; // invalid: b0&b1=1 requires b0=1
        }
        let d = b0 + b0b1; // (0,0)->0, (1,0)->1, (1,1)->2
        acc += d * POW3[i as usize];
    }
    debug_assert!(acc < 243);
    Some(acc as u8)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::point_add::trailmix_port::circuit::Circuit;

    fn valid_inputs() -> Vec<u16> {
        (0u16..1024)
            .filter(|&x| compress_classical_5(x).is_some())
            .collect()
    }

    #[test]
    fn compress_classical_5_is_bijection_onto_0_242() {
        let valid = valid_inputs();
        assert_eq!(valid.len(), 243, "must be 3^5 = 243 valid inputs");
        let mut seen = std::collections::HashSet::new();
        for &x in &valid {
            let v = compress_classical_5(x).unwrap();
            assert!(v < 243, "value {v} out of range for input {x:#012b}");
            assert!(seen.insert(v), "collision at value {v}");
        }
        assert_eq!(seen.len(), 243);
        // invalid encodings rejected
        assert!(compress_classical_5(0b0000000010).is_none()); // pair0 = (0,1)
        assert!(compress_classical_5(0b1000000000).is_none()); // pair4 = (0,1)
    }

    /// The gate sequence agrees with the classical reference on all 243
    /// valid inputs (run in shot-parallel batches of 64), with w[8],w[9]
    /// cleared and the circuit phase-clean.
    #[test]
    fn compress_5iter_quantum_matches_classical() {
        let valid = valid_inputs();
        for chunk in valid.chunks(64) {
            let mut circ = Circuit::new();
            let w = circ.alloc_qreg_bits("c5_in", 10);
            for (shot, &x) in chunk.iter().enumerate() {
                circ.sim_load_reg_bytes_shot(&w, &[(x & 0xff) as u8, (x >> 8) as u8], shot);
            }
            let wr: [&QReg; 10] = std::array::from_fn(|k| &w[k]);
            let ds: Vec<QReg> = (0..8).map(|_| circ.alloc_qreg("c5_d")).collect();
            let dirty: Vec<&QReg> = ds.iter().collect();
            compress_5iter_refs(&mut circ, &wr, &dirty);
            drop(dirty);
            for q in ds.into_iter().rev() {
                circ.zero_and_free(q);
            }
            circ.assert_phase_clean();

            let (sim, detached) = circ.destroy_sim(w);
            for (shot, &x) in chunk.iter().enumerate() {
                let want = compress_classical_5(x).unwrap();
                let mut got = 0u16;
                for k in 0..8 {
                    if sim.read_bit_shot(&detached[k], shot) == 1 {
                        got |= 1 << k;
                    }
                }
                assert_eq!(
                    got as u8, want,
                    "shot {shot} input={x:#012b}: got={got} want={want}"
                );
                assert_eq!(
                    sim.read_bit_shot(&detached[8], shot),
                    0,
                    "shot {shot}: w[8] dirty"
                );
                assert_eq!(
                    sim.read_bit_shot(&detached[9], shot),
                    0,
                    "shot {shot}: w[9] dirty"
                );
            }
        }
    }

    #[test]
    fn uncompress_classical_5_is_inverse() {
        for x in valid_inputs() {
            let packed = compress_classical_5(x).unwrap();
            assert_eq!(uncompress_classical_5(packed), Some(x), "input {x:#012b}");
        }
        assert!(uncompress_classical_5(243).is_none());
        assert!(uncompress_classical_5(255).is_none());
    }

    /// swapper5 loads a `(b0,b0&b1)` pair into a slot of the pack and
    /// unloads it, matching the classical pack of the updated window.
    /// Checks self-inverse (load+unload == identity), ancilla/phase clean.
    /// `slot` is a classical parameter, so one circuit per slot.
    #[test]
    fn swapper5_load_unload_matches_classical() {
        let pair_bits = [(0u16, 0u16), (1, 0), (1, 1)];
        for slot in 0..5usize {
            let mut scases: Vec<(u16, (u16, u16))> = Vec::new();
            'fill: for base in valid_inputs() {
                if (base >> (2 * slot)) & 0b11 != 0 {
                    continue;
                }
                for &p in &pair_bits {
                    scases.push((base, p));
                    if scases.len() == 64 {
                        break 'fill;
                    }
                }
            }

            let mut circ = Circuit::new();
            let pack = circ.alloc_qreg_bits("swp5_pack", 8);
            let bb0 = circ.alloc_qreg("swp5_bb0");
            let bb1 = circ.alloc_qreg("swp5_bb1");
            for (shot, (base, (b0, b0b1))) in scases.iter().enumerate() {
                let packed = compress_classical_5(*base).unwrap();
                circ.sim_load_reg_bytes_shot(&pack, &[packed], shot);
                circ.sim_load_reg_bytes_shot(std::slice::from_ref(&bb0), &[*b0 as u8], shot);
                circ.sim_load_reg_bytes_shot(std::slice::from_ref(&bb1), &[*b0b1 as u8], shot);
            }
            // load
            swapper5(&mut circ, slot, &[&bb0, &bb1], &pack);
            // unload
            swapper5(&mut circ, slot, &[&bb0, &bb1], &pack);
            circ.assert_phase_clean();

            let mut outs = pack;
            outs.push(bb0);
            outs.push(bb1);
            let (sim, det) = circ.destroy_sim(outs);
            for (shot, (base, (b0, b0b1))) in scases.iter().enumerate() {
                let want_pack = compress_classical_5(*base).unwrap();
                let mut got_pack = 0u16;
                for k in 0..8 {
                    if sim.read_bit_shot(&det[k], shot) == 1 {
                        got_pack |= 1 << k;
                    }
                }
                assert_eq!(
                    got_pack as u8, want_pack,
                    "slot {slot} shot {shot}: pack not restored after load+unload"
                );
                assert_eq!(sim.read_bit_shot(&det[8], shot) as u16, *b0, "bb0 restored");
                assert_eq!(
                    sim.read_bit_shot(&det[9], shot) as u16,
                    *b0b1,
                    "bb1 restored"
                );
            }
        }
    }

    /// Cost probe: Toffoli per compress_5iter (drives the apply_bv
    /// per-iter vs per-window decision -- the swapper runs it ~3560x in
    /// the EC-add, so each Toffoli here is ~3560 in the total).
    #[test]
    fn compress_5iter_cost_probe() {
        let mut circ = Circuit::new();
        let w = circ.alloc_qreg_bits("cost", 10);
        let wr: [&QReg; 10] = std::array::from_fn(|k| &w[k]);
        let ds: Vec<QReg> = (0..8).map(|_| circ.alloc_qreg("c5_d")).collect();
        let dirty: Vec<&QReg> = ds.iter().collect();
        let ccx0 = circ.ccx_emitted;
        let ccz0 = circ.ccz_emitted;
        let ops0 = circ.ops.len();
        compress_5iter_refs(&mut circ, &wr, &dirty);
        let tof = (circ.ccx_emitted - ccx0) + (circ.ccz_emitted - ccz0);
        let ops = circ.ops.len() - ops0;
        eprintln!("  compress_5iter cost: tof={tof} ops={ops}");
        drop(dirty);
        for q in ds.into_iter().rev() {
            circ.zero_and_free(q);
        }
        let _ = circ.destroy_sim(w);
    }

    /// Forward then reverse restores the original 10-bit input.
    #[test]
    fn compress_5iter_roundtrip() {
        let valid = valid_inputs();
        for chunk in valid.chunks(64) {
            let mut circ = Circuit::new();
            let w = circ.alloc_qreg_bits("c5_rt", 10);
            for (shot, &x) in chunk.iter().enumerate() {
                circ.sim_load_reg_bytes_shot(&w, &[(x & 0xff) as u8, (x >> 8) as u8], shot);
            }
            let wr: [&QReg; 10] = std::array::from_fn(|k| &w[k]);
            let ds: Vec<QReg> = (0..8).map(|_| circ.alloc_qreg("c5_d")).collect();
            let dirty: Vec<&QReg> = ds.iter().collect();
            compress_5iter_refs(&mut circ, &wr, &dirty);
            compress_5iter_reverse_refs(&mut circ, &wr, &dirty);
            drop(dirty);
            for q in ds.into_iter().rev() {
                circ.zero_and_free(q);
            }
            circ.assert_phase_clean();

            let (sim, detached) = circ.destroy_sim(w);
            for (shot, &x) in chunk.iter().enumerate() {
                let mut got = 0u16;
                for k in 0..10 {
                    if sim.read_bit_shot(&detached[k], shot) == 1 {
                        got |= 1 << k;
                    }
                }
                assert_eq!(got, x, "shot {shot}: roundtrip failed");
            }
        }
    }
}
