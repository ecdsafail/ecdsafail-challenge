//! Alg 2: Euclidean algorithm packing per-iter `(b0, b0&b1)` choices
//! into a garbage bit-vector. 1.413n + `c_iter`*sqrt(n) iters,
//! register-shared `u, v` with `c_pad`*sqrt(n) padding, total
//! 2.355n + O(sqrt(n)) garbage bits via the Fig.1 compression.
//!
//! Direct port of Schrottenloher's `to_bitvector` (`gcd_functions.py:51`) and
//! `ToBitVector` (gcd.py:150).

use crate::point_add::trailmix_port::arith::schrottenloher::gcd_compress5::{compress_classical_5, uncompress_classical_5};
use crate::point_add::trailmix_port::circuit::{Circuit, QReg};
use num_bigint::BigUint;
use num_traits::{One, Zero};

/// Dialog window geometry. M=5 base-3 pairs pack into PACK=8 bits
/// (3^5 = 243 < 256), density 8/5 = 1.600 vs the older M=3's 5/3 =
/// 1.667. Iters round UP to a multiple of M so every window is full
/// (no partial-window handling). See `gcd_compress5` for the in-place
/// radix-3 encoder.
pub const DIALOG_M: usize = 5;
pub const DIALOG_PACK: usize = 8;

/// Measurement-vent ancilla budget for the GCD `csub`'s controlled hybrid
/// adder. Each vent ancilla turns one carry-uncompute Toffoli into a
/// measurement (`3n - 2 - vents` Toffoli per csub) at a cost of +1 peak qubit.
/// The inversion's peak is set by the `apply_bv` reconstruction phase, NOT the
/// GCD, so the GCD csub has full headroom: `usize::MAX` clamps to `current_n - 1`
/// per iteration (full Gidney 2n adder). The standalone inversion peaks at 1191q
/// (apply_bv); within the low-qubit EC-add config the peak is 1173q.
/// (Schrottenloher's space-opt config caps this at `ancilla_budget=183` because
/// their peak IS during the GCD; ours is not, so we go all the way.)
const GCD_CSUB_VENTS: usize = usize::MAX;

/// Empirical iteration-count slack (4 standard deviations). Matches
/// Schrottenloher's `ITERATIONS_VAR = 2.4`.
pub const ITERATIONS_VAR: f64 = 2.4;

/// Empirical padding for register-sharing. Matches Schrottenloher's
/// `U_PAD_VAR = 2.3`.
pub const U_PAD_VAR: f64 = 2.3;

/// Number of GCD iterations for an n-bit modulus. Rounded UP to the
/// nearest multiple of `DIALOG_M` so the garbage register is exactly
/// `DIALOG_PACK` bits per `DIALOG_M` iters (one full compression block).
#[must_use]
pub fn expected_iterations(n: usize) -> usize {
    let raw = (1.413_f64 * n as f64) + ITERATIONS_VAR * (n as f64).sqrt();
    ((raw / DIALOG_M as f64).ceil() as usize) * DIALOG_M
}

/// Padding qubits added to each of `u`, `v` so register-sharing with
/// the garbage tape stays valid.
#[must_use]
pub fn u_padding(n: usize) -> usize {
    (U_PAD_VAR * (n as f64).sqrt()).ceil() as usize
}

/// Total garbage bit-vector width for an n-bit modulus. Exactly
/// `DIALOG_PACK` bits per `DIALOG_M` GCD iterations.
#[must_use]
pub fn expected_garbage(n: usize) -> usize {
    expected_iterations(n) / DIALOG_M * DIALOG_PACK
}

/// Packed bits per dialog window for a given window size. Only the two
/// supported modes: `dialog_m = 5` (base-3 radix compression, 8 bits per
/// 5 pairs — the SP1-validated low-qubit config) and `dialog_m = 1` (RAW,
/// 2 bits per pair, no compression — the higher-qubit / lower-Toffoli
/// config: it drops every `compress_5iter` call at the cost of a wider
/// garbage tape).
#[must_use]
pub fn dialog_pack(dialog_m: usize) -> usize {
    match dialog_m {
        1 => 2,           // RAW: 2 bits/pair, 0 compression Toffoli
        3 => 5,           // M=3 SAT compressor: 3 pairs -> 5 bits, ~5 tof/window
        5 => DIALOG_PACK, // M=5 base-3 arithmetic: 5 pairs -> 8 bits, ~316 tof/window
        _ => panic!("dialog_m must be 1 (raw), 3 (SAT), or 5 (arith), got {dialog_m}"),
    }
}

/// Garbage tape width for a given window size. The iteration count is the
/// same as `expected_iterations` (convergence requirement, a multiple of
/// 5); only the per-window packing changes.
#[must_use]
pub fn expected_garbage_m(n: usize, dialog_m: usize) -> usize {
    expected_iterations(n) / dialog_m * dialog_pack(dialog_m)
}

/// Classical reference: forward GCD of `(u, v)` packing
/// `(b0, b0&b1)` choices into a compressed garbage bit-vector. After
/// `expected_iterations(n)` iters, the algorithm should reach
/// `(u, v) = (1, 0)`. Returns the garbage bit-vector on success, or
/// `None` if the iteration budget didn't suffice (statistical tail
/// per `ITERATIONS_VAR`).
///
/// Mirrors Schrottenloher's `to_bitvector(u, v)` (`gcd_functions.py:51`).
///
/// Pre: `u` is odd, `u, v < 2^n`.
#[must_use]
pub fn to_bitvector_classical(mut u: BigUint, mut v: BigUint, n: usize) -> Option<Vec<u8>> {
    assert!(u.bit(0), "u must be odd");

    let iters = expected_iterations(n);
    let n_garbage_bits = iters / DIALOG_M * DIALOG_PACK;
    let n_garbage_bytes = n_garbage_bits.div_ceil(8);
    let mut garbage = vec![0u8; n_garbage_bytes];

    // compress_classical_5(all-(0,0)) = 0, so the all-zero garbage already
    // encodes every window as five (0,0) pairs -- no init pass needed.

    for i in 0..iters {
        let b0 = u8::from(v.bit(0));
        let b1 = u8::from(u > v);
        let b0_and_b1 = b0 & b1;

        // Read current 8-bit pack at offset DIALOG_PACK * (i / DIALOG_M).
        let pack_off_bit = DIALOG_PACK * (i / DIALOG_M);
        let mut pack = 0u8;
        for k in 0..DIALOG_PACK {
            let g_bit = (garbage[(pack_off_bit + k) / 8] >> ((pack_off_bit + k) % 8)) & 1;
            pack |= g_bit << k;
        }
        // Decompress and write the (b0, b0&b1) pair at slot (i % DIALOG_M).
        let decompressed =
            uncompress_classical_5(pack).expect("decompress must succeed on prior compress output");
        let slot = (i % DIALOG_M) * 2;
        // Sanity: the slot must currently be (0,0).
        debug_assert_eq!(
            (decompressed >> slot) & 0b11,
            0,
            "absorber slot {slot} must be (0,0) pre-fill"
        );
        let new_decomp =
            decompressed | (u16::from(b0) << slot) | (u16::from(b0_and_b1) << (slot + 1));
        let new_pack =
            compress_classical_5(new_decomp).expect("recompress must succeed on valid pair");
        // Write back.
        for k in 0..DIALOG_PACK {
            let want_bit = (new_pack >> k) & 1;
            let idx = pack_off_bit + k;
            let mask = 1u8 << (idx % 8);
            if want_bit == 1 {
                garbage[idx / 8] |= mask;
            } else {
                garbage[idx / 8] &= !mask;
            }
        }

        if b0_and_b1 == 1 {
            std::mem::swap(&mut u, &mut v);
        }
        if b0 == 1 {
            v = if v >= u {
                v - &u
            } else {
                (&u + &v) - 2u32 * &u
            };
            // After (u, v) potentially swapped, v >= u when b0_and_b1=1. When
            // b0=1 and b0_and_b1=0, b1=0 means u <= v, so v >= u. Safe.
        }
        v >>= 1;
    }
    if !v.is_zero() || u != BigUint::one() {
        return None; // budget exhausted before convergence
    }
    Some(garbage)
}

/// Inverse of `to_bitvector_classical`: read the garbage in reverse,
/// rebuild `(u, v)` from the dialog. Mirrors `from_bitvector`
/// (`gcd_functions.py:123`).
#[must_use]
pub fn from_bitvector_classical(garbage: &[u8], n: usize) -> (BigUint, BigUint) {
    let iters = expected_iterations(n);
    let mut u = BigUint::one();
    let mut v = BigUint::zero();
    for i in (0..iters).rev() {
        let pack_off_bit = DIALOG_PACK * (i / DIALOG_M);
        let mut pack = 0u8;
        for k in 0..DIALOG_PACK {
            let g_bit = (garbage[(pack_off_bit + k) / 8] >> ((pack_off_bit + k) % 8)) & 1;
            pack |= g_bit << k;
        }
        let decompressed = uncompress_classical_5(pack).expect("valid pack expected");
        let slot = (i % DIALOG_M) * 2;
        let b0 = (decompressed >> slot) & 1;
        let b0_and_b1 = (decompressed >> (slot + 1)) & 1;

        v <<= 1;
        if b0 == 1 {
            v += &u;
        }
        if b0_and_b1 == 1 {
            std::mem::swap(&mut u, &mut v);
        }
    }
    (u, v)
}

/// Quantum forward GCD packing for secp256k1. Drives `(u, v)` from
/// `(q, x)` to `(1, 0)`, packing per-iter `(b0, b0_and_b1)` choices
/// into `garbage`. After the call, the caller's `v_full` register is
/// zeroed (the value has been consumed into the garbage); `u_full`
/// and all ancillae are internally cleaned to |0>.
///
/// Inputs:
///   - `v_full`: caller-allocated `n + u_padding(n)` bits. Pre: low
///     `n` bits hold the value `x` to encode; high `u_padding(n)`
///     bits are 0. Post: all bits 0.
///   - `garbage`: caller-allocated `expected_garbage(n)` bits, pre 0.
///     Post: GCD-step pack matching `to_bitvector_classical(q, x)`.
///
/// Internally allocates `u_full = n + u_padding(n)` bits initialized
/// to `q = 2^n - F_SECP256K1` (via X gates on the appropriate bits),
/// plus per-iter `b0, b0_and_b1` and `pack_anc` ancillae. All
/// internal ancillae are cleaned to |0> on exit.
///
/// Matches Schrottenloher's `ToBitVector` (gcd.py:150) at secp256k1 n=256.
/// Forward GCD packing with full per-iter register shrinking.
/// Both `u_full` (internal) AND `v_full` (caller-owned, passed as
/// &mut Vec) are popped as the `current_n` schedule shrinks; the
/// freed bits enter the allocator's free pool and feed the garbage
/// register's incremental allocation. Returns the
/// `expected_garbage(n) = 648` bit garbage register (M=5 tape).
///
/// After the call, `v_full`'s Vec has length 1 (just bit 0, holding
/// the residual value before the final X-flip). The reverse pass
/// regrows it.
pub fn forward_gcd_pack_quantum_secp256k1(circ: &mut Circuit, v_full: &mut Vec<QReg>) -> Vec<QReg> {
    forward_gcd_pack_quantum_secp256k1_m(circ, v_full, DIALOG_M)
}

/// Window-size-parameterized forward GCD pack. `dialog_m = 5` is the
/// base-3-compressed (SP1-validated) low-qubit config; `dialog_m = 1` is
/// the RAW (2 bits/pair, no `compress_5iter`) higher-qubit / lower-Toffoli
/// config. The iteration count is identical (convergence); only the
/// per-window packing differs.
pub fn forward_gcd_pack_quantum_secp256k1_m(
    circ: &mut Circuit,
    v_full: &mut Vec<QReg>,
    dialog_m: usize,
) -> Vec<QReg> {
    use crate::point_add::trailmix_port::arith::schrottenloher::gcd_compress5::compress_5iter_refs;
    use crate::point_add::trailmix_port::arith::schrottenloher::msb_compare;

    let n = 256usize;
    let pad = u_padding(n);
    assert!(
        v_full.len() >= n,
        "v_full must be at least n = 256 bits (pad above n is never active in the current_n schedule)"
    );
    let pack = dialog_pack(dialog_m);
    let iters = expected_iterations(n);
    let garbage_len = iters / dialog_m * pack;

    let prev = circ.push_section("fwd_gcd");

    // Load u_full = q (n bits) || 0 (pad bits). Mutable so we can
    // truncate as `current_n` shrinks per iter (strict-dealloc
    // requires inactive tail bits to be freed before the comparator
    // allocates fresh ancillae in subsequent iters).
    let mut u_full = circ.alloc_qreg_bits("u_full", n);
    let q_bigint: num_bigint::BigUint = (num_bigint::BigUint::from(1u32) << n as u32)
        - num_bigint::BigUint::from(super::F_SECP256K1);
    let q_bytes = q_bigint.to_bytes_le();
    for i in 0..n {
        let byte_idx = i / 8;
        if byte_idx < q_bytes.len() && ((q_bytes[byte_idx] >> (i % 8)) & 1) == 1 {
            circ.x(&u_full[i]);
        }
    }

    // Per-iter ancillae (reused).
    let b0 = circ.alloc_qreg("gcd.b0");
    let b0_and_b1 = circ.alloc_qreg("gcd.b0_and_b1");

    // Paper Sec. 3.1 TRUNCATE constant. Schrottenloher sets TRUNCATE = 40
    // (gcd.py header constant). Combined with u_padding the
    // approximate compare looks at top (TRUNCATE + u_padding) bits
    // of the CURRENT (shrunken) effective register width.
    const TRUNCATE: usize = 40;
    let trunc = TRUNCATE + pad;

    // Incremental garbage register. The freed bits of u_full enter the
    // allocator pool as iters progress; per-pack allocs reuse them.
    // De-facto register sharing without bit-reversed interleave.
    let mut garbage_vec: Vec<QReg> = Vec::with_capacity(garbage_len);
    // The current DIALOG_M-iter dialog window, held DECOMPRESSED (2*M bits =
    // M (0,0)/(1,0)/(1,1) pairs) and compressed ONCE at the window end. The
    // absorb is two cheap swaps per iter; the radix-3 encoder fires once per
    // M iters.
    let mut cur_win: Vec<QReg> = Vec::new();

    // STEP constant in gcd.py is used to round current_n; we follow
    // verbatim — register width shrinks per iter (paper's
    // register-sharing trick), letting the approximate compare on
    // top-`trunc` bits stay tight.
    const STEP: usize = 1;
    let current_n_at = |i: usize| -> usize {
        // current_n = ceil(n - i * 0.5 * 1.415 + u_padding), clamped
        // to [0, n]. Rounded up to STEP. (gcd.py:228-230.)
        let raw = (n as f64) - (i as f64) * 0.5 * 1.415_f64 + pad as f64;
        let raw_int = raw.ceil().max(0.0) as usize;
        let stepped = raw_int.div_ceil(STEP) * STEP;
        stepped.min(n)
    };

    for i in 0..iters {
        let iter_section = format!("fwd_iter_{i:04}");
        let iter_prev = circ.push_section(&iter_section);

        // Shrunken active width: u_reg = u_full[..current_n], same
        // for v_reg.
        let current_n = current_n_at(i);

        // Free u_full AND v_full tail bits that just became inactive.
        // They were never touched in this or any subsequent iter and
        // are |0> (cswap only operates on the active prefix; the
        // ctrl-sub and right_shift work in-slice). zero_and_free
        // returns the qubits to the allocator pool for garbage reuse.
        while u_full.len() > current_n.max(1) {
            let q = u_full.pop().expect("u_full nonempty");
            circ.zero_and_free(q);
        }
        while v_full.len() > current_n.max(1) {
            let q = v_full.pop().expect("v_full nonempty");
            circ.zero_and_free(q);
        }

        // Top `trunc` bits inside the shrunken view are
        // u_reg[current_n - trunc..current_n].
        let cmp_eff = trunc.min(current_n);
        let cmp_lo = current_n.saturating_sub(cmp_eff);

        // 1) b0 = v_full[0] (parity of v).
        circ.cx(&v_full[0], &b0);

        // 2) b0_and_b1 ^= b0 AND (v_top < u_top) on top-`trunc` bits
        //    of the shrunken effective width. Gidney measure-uncompute
        //    comparator (n Toffoli, +cmp_eff peak) -- the GCD has the
        //    ancilla headroom (peak is set by apply_bv, not here).
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

        // 6-7) Per-window dialog packing. At a fresh window boundary,
        //    allocate the 2*M-bit DECOMPRESSED window (all |0> = M (0,0)
        //    pairs; freed u_full bits feed the allocator pool). Absorb each
        //    (b0, b0_and_b1) into its slot with two swaps (b0/b0_and_b1
        //    return to |0> = the old (0,0) slot value). At the window end
        //    compress 2*M -> PACK + (2*M - PACK) |0>, keep the PACK packed
        //    bits, and free the zeroed high bits. Identical garbage content
        //    to the per-iter swapper, but the encoder fires once per M iters.
        let slot = i % dialog_m;
        if slot == 0 {
            cur_win = (0..2 * dialog_m)
                .map(|_| circ.alloc_qreg("eea_win"))
                .collect();
        }
        circ.swap(&b0, &cur_win[2 * slot]);
        circ.swap(&b0_and_b1, &cur_win[2 * slot + 1]);
        if slot == dialog_m - 1 {
            if dialog_m == DIALOG_M {
                // Fresh |0> comparator dirty (8 bits). The compress runs below the
                // csub peak (the GCD's binding op), so +8 here doesn't raise the
                // GCD-phase peak; apply_bv borrows instead.
                let ds: Vec<QReg> = (0..8).map(|_| circ.alloc_qreg("c5_dirty")).collect();
                {
                    let win_refs: [&QReg; 10] = std::array::from_fn(|k| &cur_win[k]);
                    let dirty: Vec<&QReg> = ds.iter().collect();
                    compress_5iter_refs(circ, &win_refs, &dirty);
                }
                for q in ds.into_iter().rev() {
                    circ.zero_and_free(q);
                }
                // free the (2*M - PACK) high zeros, keep PACK packed bits
                while cur_win.len() > pack {
                    let freed = cur_win.pop().expect("high window bit |0> post-compress");
                    circ.zero_and_free(freed);
                }
            } else if dialog_m == 3 {
                // M=3 SAT permutation compressor: 6 pair-bits -> 5 packed (+1 |0>),
                // ~5 Toffoli, no scratch (gcd_compress.rs).
                {
                    let win_refs: [&QReg; 6] = std::array::from_fn(|k| &cur_win[k]);
                    crate::point_add::trailmix_port::arith::schrottenloher::gcd_compress::compress_3iter_refs(
                        circ, &win_refs,
                    );
                }
                while cur_win.len() > pack {
                    let freed = cur_win.pop().expect("high window bit |0> post-compress");
                    circ.zero_and_free(freed);
                }
            }
            // RAW (dialog_m == 1): cur_win is the 2 raw (b0, b0&b1) bits;
            // no compression, append directly (pack == 2 == 2*dialog_m).
            garbage_vec.append(&mut cur_win); // packed bits, in window order
        }

        // 8) Strict-dealloc touch: bump last_touched_op on every alive
        //    u_full bit past swp_anc's allocation. Without this, when a
        //    later iter (or function-end) frees u_full[k], the check
        //    sees swp_anc's alloc happened after u_full[k]'s last data
        //    touch in cswap/csub, and fires. b0 and b0_and_b1 are both
        //    |0> here, so two CX gates with different controls have no
        //    semantic effect and aren't redundancy-elided (different
        //    operands).
        for k in 0..current_n {
            circ.cx(&b0, &u_full[k]);
            circ.cx(&b0_and_b1, &u_full[k]);
        }

        circ.pop_section(&iter_prev);
    }

    // End: u_full ≈ 1, v_full ≈ 0 (mod q drift); X u_full[0] to drive
    // u_full toward 0. Approximate primitives may have left residual
    // multiples of q in u_full / v_full; we do NOT assert exact 0
    // here. Callers / S9 will compose with the reverse pass that
    // restores the value cleanly.
    circ.x(&u_full[0]);

    // The Euclidean (u, v) recurrence is EXACT integer arithmetic
    // (Cuccaro adds + shifts), so v converges to exactly 0. The schedule
    // floor leaves v_full at current_n_min (~10) bits, all holding 0.
    // Free them down to a single bit so the caller's x register is not a
    // ~10-qubit dead weight during the apply_bitvector peak. The reverse
    // GCD regrows v_full from length 1 symmetrically.
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

/// Reverse of `forward_gcd_pack_quantum_secp256k1`. Restores
/// `v_full[..n]` to its original value `x` (consumed by the forward
/// pass) and drives `garbage` back to all-|0>. Uses a fresh `u_full`
/// internally — the forward's `u_full` was freed, so the reverse must
/// regrow it.
///
/// Inputs:
///   - `v_full`: caller-allocated `n + u_padding(n)` bits, currently 0
///     (per the forward's exit state).
///   - `garbage`: caller-allocated `expected_garbage(n)` bits, currently
///     filled with the forward's per-iter pack.
///
/// Outputs:
///   - `v_full[..n]` = original `x`, `v_full[n..total]` = 0.
///   - `garbage` = all 0.
pub fn forward_gcd_pack_quantum_secp256k1_reverse(
    circ: &mut Circuit,
    v_full: &mut Vec<QReg>,
    garbage: &mut Vec<QReg>,
) {
    forward_gcd_pack_quantum_secp256k1_reverse_m(circ, v_full, garbage, DIALOG_M);
}

/// Window-size-parameterized reverse GCD pack (see
/// `forward_gcd_pack_quantum_secp256k1_m`). `dialog_m = 1` drains a RAW
/// 2-bit-per-pair tape (no decompression); `dialog_m = 5` decompresses
/// the base-3 packs.
pub fn forward_gcd_pack_quantum_secp256k1_reverse_m(
    circ: &mut Circuit,
    v_full: &mut Vec<QReg>,
    garbage: &mut Vec<QReg>,
    dialog_m: usize,
) {
    use crate::point_add::trailmix_port::arith::schrottenloher::gcd_compress5::compress_5iter_reverse_refs;
    use crate::point_add::trailmix_port::arith::schrottenloher::msb_compare;

    let n = 256usize;
    let pad = u_padding(n);
    let pack = dialog_pack(dialog_m);
    // Reverse GROWS v_full back from whatever size it has (usually 1
    // after forward shrunk it). The shrink-grow API is symmetric.
    let iters = expected_iterations(n);
    let garbage_len = iters / dialog_m * pack;
    assert!(garbage.len() >= garbage_len);

    let prev = circ.push_section("fwd_gcd_rev");

    // Fresh u_full, starting empty: bit 0 first, then grow as the
    // reverse iter sequence asks for more bits.
    let mut u_full: Vec<QReg> = Vec::with_capacity(n);
    u_full.push(circ.alloc_qreg("u_full_re[0]"));

    // The forward terminated with: u_full[0] = 1 (algorithm's u_final),
    // then X(u_full[0]) flipped it to 0, then dropped. To reverse:
    // start u_full[0] = 0 (fresh-alloc), apply X(u_full[0]) (= reverse
    // of the forward's X), giving u_full[0] = 1.
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
    // Current window, decompressed while its pairs are extracted; see the
    // forward pass for the per-window scheme (DIALOG_M pairs per window).
    let mut cur_win: Vec<QReg> = Vec::new();

    for i in (0..iters).rev() {
        let iter_section = format!("rev_iter_{i:04}");
        let iter_prev = circ.push_section(&iter_section);

        let current_n = current_n_at(i).max(1);
        // Grow u_full back to current_n bits. New bits are fresh |0>.
        while u_full.len() < current_n {
            u_full.push(circ.alloc_qreg("u_full_re"));
        }
        // Symmetrically regrow v_full as iter sequence asks.
        while v_full.len() < current_n {
            v_full.push(circ.alloc_qreg("v_full_re"));
        }

        // Inverse of the forward per-window packing. At a window start
        // (reverse order => slot M-1), pull this window's PACK packed bits
        // off the end of the garbage tape and DECOMPRESS them (+ 2*M-PACK
        // fresh |0>) into the 2*M-bit window. Extract the slot's pair into
        // (b0, b0_and_b1); the reverse body below re-derives and cancels
        // them back to |0>, zeroing the slot. At the window end (slot 0)
        // the window is all |0> and is freed (draining the tape).
        let slot = i % dialog_m;
        if slot == dialog_m - 1 {
            let glen = garbage.len();
            cur_win = garbage.split_off(glen - pack);
            if dialog_m == DIALOG_M {
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
            } else if dialog_m == 3 {
                // M=3 SAT decompress: 5 packed -> pad to 6 (+1 |0>), inverse permute.
                for _ in 0..(2 * 3 - 5) {
                    cur_win.push(circ.alloc_qreg("eea_win_re"));
                }
                {
                    let win_refs: [&QReg; 6] = std::array::from_fn(|k| &cur_win[k]);
                    crate::point_add::trailmix_port::arith::schrottenloher::gcd_compress::compress_3iter_reverse_refs(
                        circ, &win_refs,
                    );
                }
            }
            // RAW (dialog_m == 1): cur_win is exactly the 2 raw pair bits
            // split off the tape; no decompression.
        }
        circ.swap(&b0, &cur_win[2 * slot]);
        circ.swap(&b0_and_b1, &cur_win[2 * slot + 1]);

        // Inverse of right_shift = left_shift (rotates the other way).
        let v_slice = &v_full[..current_n];
        crate::point_add::trailmix_port::arith::shift::left_shift(circ, v_slice);

        // Inverse of csub (X-sandwich + cadd). Sandwich Xs cancel; net
        // inverse is plain controlled_add on the SAME slice.
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

        // cswap is involutory.
        let s_cs = circ.push_section("rev_cswap");
        for j in 0..current_n {
            circ.cswap(&b0_and_b1, &u_full[j], &v_full[j]);
        }
        circ.pop_section(&s_cs);

        // compare is involutory (controlled_lt_msbs is compute/use/uncompute).
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

        // cx(v[0], b0) is involutory.
        circ.cx(&v_full[0], &b0);

        // End-of-iter strict-dealloc touch on u_full[..current_n] (see
        // matching block in forward GCD). b0, b0_and_b1 are both |0>
        // here; the two CX gates with different controls are
        // semantically no-op and not redundancy-elided.
        for k in 0..current_n {
            circ.cx(&b0, &u_full[k]);
            circ.cx(&b0_and_b1, &u_full[k]);
        }

        circ.pop_section(&iter_prev);

        // Window end (reverse, slot 0): every pair has been extracted, so
        // the decompressed window is all |0>; free its 2*M qubits, draining
        // this window off the garbage tape (PACK packed bits were split off
        // at slot M-1; the rest were decompress ancillae). RAW: just the 2
        // pair bits, now |0>.
        if i % dialog_m == 0 {
            for q in cur_win.drain(..) {
                circ.zero_and_free(q);
            }
        }
    }
    debug_assert!(garbage.is_empty(), "garbage tape not fully drained");

    // After all reverse iters, u_full should be q (the algorithm's
    // u_init in forward). Deinit q by X-flipping the set bits, then
    // free u_full as |0>.
    let q_bigint: num_bigint::BigUint = (num_bigint::BigUint::from(1u32) << n as u32)
        - num_bigint::BigUint::from(super::F_SECP256K1);
    let q_bytes = q_bigint.to_bytes_le();
    for i in 0..n.min(u_full.len()) {
        let byte_idx = i / 8;
        if byte_idx < q_bytes.len() && ((q_bytes[byte_idx] >> (i % 8)) & 1) == 1 {
            circ.x(&u_full[i]);
        }
    }
    // u_full is now |0>; free each bit in order.
    while let Some(q) = u_full.pop() {
        circ.zero_and_free(q);
    }

    drop(b0);
    drop(b0_and_b1);
    circ.pop_section(&prev);
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Round-trip: to_bitvector then from_bitvector must recover (u, v).
    /// Tested at small n (= 16) where the algorithm converges
    /// deterministically and the iteration budget is comfortably large.
    #[test]
    fn classical_gcd_roundtrip_n16() {
        let n = 16usize;
        let test_pairs: Vec<(u32, u32)> = vec![
            (1, 0),
            (3, 5),
            (7, 12),
            (15, 8),
            (0xabcd, 0x1234),
            (0xffff, 0x8000),
            (0x8001, 0x4001),
        ];
        for (uu, vv) in test_pairs {
            let u = BigUint::from(uu);
            let v = BigUint::from(vv);
            let g = to_bitvector_classical(u.clone(), v.clone(), n)
                .unwrap_or_else(|| panic!("budget exhausted on (u={uu}, v={vv})"));
            let (ru, rv) = from_bitvector_classical(&g, n);
            assert_eq!(ru, u, "u mismatch on (u={uu}, v={vv})");
            assert_eq!(rv, v, "v mismatch on (u={uu}, v={vv})");
        }
    }

    /// Quantum forward GCD: load v with a random odd value, run, then
    /// compare the QUANTUM garbage register against the classical
    /// to_bitvector(q, x) output. 64 random shots, no RNG seeding.
    #[test]
    fn forward_gcd_pack_quantum_secp256k1_garbage_correct() {
        use rand::{thread_rng, Rng};
        let n = 256usize;
        let pad = u_padding(n);
        let total = n + pad;
        let iters = expected_iterations(n);
        let garbage_len = iters / DIALOG_M * DIALOG_PACK;
        let q_le: BigUint =
            (BigUint::from(1u32) << 256u32) - BigUint::from(super::super::F_SECP256K1);

        let mut rng = thread_rng();
        let mut shot_data: Vec<(BigUint, Vec<u8>)> = Vec::new(); // (x, expected_garbage)
        while shot_data.len() < 64 {
            let raw: [u8; 32] = rng.gen();
            let mut x = BigUint::from_bytes_le(&raw) % &q_le;
            if x.is_zero() {
                x = BigUint::from(1u32);
            }
            let g = match to_bitvector_classical(q_le.clone(), x.clone(), n) {
                Some(g) => g,
                None => continue,
            };
            shot_data.push((x, g));
        }

        let mut circ = crate::point_add::trailmix_port::circuit::Circuit::new();
        // Forward-GCD peak occurs mid-run (~iter 52): u/v still full width
        // (current_n=256) while the garbage tape has accumulated as NEW
        // qubits (register-sharing only reclaims freed u/v bits once they
        // start shrinking, i>~52). Measured: 861 forward / 864 fwd+rev,
        // ~equal for M=3 and M=5. This phase legitimately needs ~864 (u/v
        // stay full-width while the garbage tape accrues), so we set a tight
        // per-test cap of 870 here, well under the EC-add apply_bv peak ~1169.
        circ.set_max_qubit_peak(870);
        let mut v_full = circ.alloc_qreg_bits("v_full", total);

        for (shot, (x, _g)) in shot_data.iter().enumerate() {
            let mut buf = [0u8; 32];
            for (i, b) in x.to_bytes_le().iter().take(32).enumerate() {
                buf[i] = *b;
            }
            circ.sim_load_reg_bytes_shot(&v_full[..n], &buf, shot);
        }

        let ops0 = circ.ops.len();
        let ccx0 = circ.ccx_emitted;
        let ccz0 = circ.ccz_emitted;
        let garbage = forward_gcd_pack_quantum_secp256k1(&mut circ, &mut v_full);
        let f_ops = circ.ops.len() - ops0;
        let f_tof = (circ.ccx_emitted - ccx0) + (circ.ccz_emitted - ccz0);
        let peak = circ.peak_qubits;
        eprintln!(
            "  cost(forward_gcd_pack_quantum_secp256k1, n=256): ops={f_ops} tof={f_tof} peak+={peak} (iters={iters})"
        );

        // forward_gcd frees the inactive u/v tail per iteration and the
        // converged v_full remnant at the end, so v_full is now ~1 bit, not
        // `total`. The garbage register follows it in `outs`, so index it at
        // the actual (shrunk) v_full length, not the original `total`.
        let vlen = v_full.len();
        let mut outs: Vec<crate::point_add::trailmix_port::circuit::QReg> = Vec::new();
        outs.extend(v_full);
        outs.extend(garbage);
        let (sim, detached) = circ.destroy_sim(outs);

        for (shot, (_x, expected_g)) in shot_data.iter().enumerate() {
            for i in 0..garbage_len {
                let want_bit = (expected_g[i / 8] >> (i % 8)) & 1;
                let got_bit = sim.read_bit_shot(&detached[vlen + i], shot);
                assert_eq!(
                    got_bit, want_bit,
                    "shot {} bit {}: garbage mismatch want={} got={}",
                    shot, i, want_bit, got_bit
                );
            }
        }
    }

    /// Forward GCD + reverse GCD should restore v_full to its initial
    /// value and zero the garbage register. 64 random shots.
    #[test]
    fn forward_gcd_pack_quantum_roundtrip_secp256k1() {
        use rand::{thread_rng, Rng};
        let n = 256usize;
        let pad = u_padding(n);
        let total = n + pad;
        let iters = expected_iterations(n);
        let q_le: BigUint =
            (BigUint::from(1u32) << 256u32) - BigUint::from(super::super::F_SECP256K1);

        let mut rng = thread_rng();
        let mut shot_data: Vec<BigUint> = Vec::new();
        while shot_data.len() < 64 {
            let raw: [u8; 32] = rng.gen();
            let mut x = BigUint::from_bytes_le(&raw) % &q_le;
            if x.is_zero() {
                x = BigUint::from(1u32);
            }
            // Make sure the budget will succeed on this x; otherwise
            // skip.
            if to_bitvector_classical(q_le.clone(), x.clone(), n).is_none() {
                continue;
            }
            shot_data.push(x);
        }

        let mut circ = crate::point_add::trailmix_port::circuit::Circuit::new();
        // Forward-GCD peak occurs mid-run (~iter 52): u/v still full width
        // (current_n=256) while the garbage tape has accumulated as NEW
        // qubits (register-sharing only reclaims freed u/v bits once they
        // start shrinking, i>~52). Measured: 861 forward / 864 fwd+rev,
        // ~equal for M=3 and M=5. This phase legitimately needs ~864 (u/v
        // stay full-width while the garbage tape accrues), so we set a tight
        // per-test cap of 870 here, well under the EC-add apply_bv peak ~1169.
        circ.set_max_qubit_peak(870);
        let mut v_full = circ.alloc_qreg_bits("v_full", total);

        for (shot, x) in shot_data.iter().enumerate() {
            let mut buf = [0u8; 32];
            for (i, b) in x.to_bytes_le().iter().take(32).enumerate() {
                buf[i] = *b;
            }
            circ.sim_load_reg_bytes_shot(&v_full[..n], &buf, shot);
        }

        let ops0 = circ.ops.len();
        let ccx0 = circ.ccx_emitted;
        let ccz0 = circ.ccz_emitted;
        let mut garbage = forward_gcd_pack_quantum_secp256k1(&mut circ, &mut v_full);
        forward_gcd_pack_quantum_secp256k1_reverse(&mut circ, &mut v_full, &mut garbage);
        assert!(garbage.is_empty(), "roundtrip should drain garbage tape");
        let rt_ops = circ.ops.len() - ops0;
        let rt_tof = (circ.ccx_emitted - ccx0) + (circ.ccz_emitted - ccz0);
        let peak = circ.peak_qubits;
        eprintln!(
            "  cost(forward+reverse GCD secp256k1): ops={rt_ops} tof={rt_tof} peak+={peak} (iters={iters})"
        );

        // After fwd+rev the register-sharing schedule has freed the pad bits
        // and drained the garbage tape, so v_full now holds the original
        // value in its low n bits with no pad/garbage qubits remaining
        // (`garbage.is_empty()` above is the drain check). Phase must be
        // clean: every measurement-vent ghost in the hybrid-adder csubs was
        // resolved by its own close_ghost.
        circ.assert_phase_clean();
        let mut outs: Vec<crate::point_add::trailmix_port::circuit::QReg> = Vec::new();
        outs.extend(v_full);
        let (sim, detached) = circ.destroy_sim(outs);

        for (shot, x) in shot_data.iter().enumerate() {
            // v_full[..n] must equal x
            let mut got_v = BigUint::zero();
            for i in 0..n {
                if sim.read_bit_shot(&detached[i], shot) == 1 {
                    got_v |= BigUint::from(1u32) << i;
                }
            }
            assert_eq!(got_v, *x, "shot {} v not restored", shot);
        }
    }

    /// At n=256, run on a handful of random (q, x) pairs (matching the
    /// invert-x-mod-q use case). u_init = q, v = x.
    #[test]
    fn classical_gcd_random_n256() {
        use rand::rngs::StdRng;
        use rand::{Rng, SeedableRng};

        let n = 256usize;
        // secp256k1 q
        let q_le: BigUint =
            (BigUint::from(1u32) << 256u32) - BigUint::from(super::super::F_SECP256K1);
        let mut rng = StdRng::seed_from_u64(0xc4d_ec_011a_b1cd);
        let n_trials = 16usize;
        let mut budget_exhausted = 0usize;
        for trial in 0..n_trials {
            let raw: [u8; 32] = rng.gen();
            let mut x = BigUint::from_bytes_le(&raw) % &q_le;
            if x.is_zero() {
                x = BigUint::one();
            }
            let g = to_bitvector_classical(q_le.clone(), x.clone(), n);
            match g {
                None => {
                    budget_exhausted += 1;
                }
                Some(g) => {
                    let (ru, rv) = from_bitvector_classical(&g, n);
                    assert_eq!(ru, q_le, "trial {trial}: u recovery");
                    assert_eq!(rv, x, "trial {trial}: v recovery");
                }
            }
        }
        eprintln!("n=256 classical GCD: {budget_exhausted}/{n_trials} budget-exhausted");
        // Per the paper, 1 in 10000-20000 inputs hit the iteration tail;
        // with 16 trials we expect ~0 exhaustions. >2 is suspicious.
        assert!(
            budget_exhausted <= 2,
            "{budget_exhausted} budget exhaustions out of {n_trials} is high"
        );
    }
}
