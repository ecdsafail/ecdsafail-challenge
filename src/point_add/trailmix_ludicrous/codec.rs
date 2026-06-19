//! Dialog (de)compression codec for the product-min jump-GCD EC-add,
//! built on this crate's `B` builder.
//!
//! The jump-GCD inversion records a per-step "dialog": a few qubits
//! per GCD step `(subtracted, swap, s_2)`. To save qubits, windows of symbols are
//! packed by this reversible codec before the reverse pass consumes them. The
//! product-min config uses the all-triple codec: each window of 3 base-5
//! symbols (9 raw bits) compresses to a 7-bit code and back. Only 5 of the 8
//! `(subtracted, swap, s_2)` patterns are reachable, so a window of 3 carries
//! log2(5^3) ~ 6.97 bits and fits in 7. The Step-0 boundary symbol and a small
//! Pair/Raw tail finish the tiling at the schedule length (`iters = 258`).
//!
//! Every ancilla is allocated via `alloc_qubit` and freed via `zero_and_free`, and every
//! AND-uncompute is a measurement-vented `clear_and` (HMR + conditional-Z,
//! zero Toffoli).

use super::{B, BExt};
use crate::circuit::{QubitId};

// ===================================================================
// clear_and: measurement-vented AND-uncompute (HMR + conditional-Z).
//
// `t` currently holds `a AND b` (e.g. from a prior `ccx(a, b, t)`). The bare
// circuit's `hmr` measures `t` out to |0> with a random outcome bit that XORs
// `a&b` into the global phase; a `cz_if_bit(a, b, outcome)` gated on that bit
// cancels exactly that phase. No Toffoli is spent (the measurement replaces the
// reverse CCX).
// ===================================================================

fn clear_and(circ: &mut B, t: &QubitId, a: &QubitId, b: &QubitId) {
    let bit = circ.alloc_bit();
    circ.hmr(*t, bit);
    circ.cz_if_bit(*a, *b, bit);
}

// ===================================================================
// Pairs core (6->5). SAT-synthesized in-place reversible circuit mapping the 25
// valid 2-symbol dialog inputs to 25 distinct 5-bit codes with wire 5 freeable.
// ===================================================================

/// Compress 2 base-5 symbols on `w[0..6]` -> 5-bit code on `w[0..5]`, `w[5]`->|0>.
fn compress_2sym_fast(circ: &mut B, w: &[&QubitId; 6]) {
    circ.x(*w[3]);
    circ.cx(*w[5], *w[1]);
    circ.cx(*w[4], *w[0]);
    circ.x(*w[2]);
    circ.ccx(*w[1], *w[3], *w[5]);
    circ.cx(*w[3], *w[5]);
    circ.cx(*w[3], *w[0]);
    circ.cx(*w[1], *w[5]);
    circ.cx(*w[5], *w[3]);
    circ.ccx(*w[5], *w[0], *w[4]);
    // Terminal AND-uncompute: w[5] holds w[3] AND w[4] here (the freed wire).
    // Vent it via clear_and; the reverse rebuilds it with a CCX.
    clear_and(circ, w[5], w[3], w[4]);
}

/// Inverse of [`compress_2sym_fast`]: the same gates reversed (each X/CX/CCX is
/// self-inverse; the terminal AND-uncompute is rebuilt with a CCX).
fn compress_2sym_fast_reverse(circ: &mut B, w: &[&QubitId; 6]) {
    circ.ccx(*w[3], *w[4], *w[5]);
    circ.ccx(*w[5], *w[0], *w[4]);
    circ.cx(*w[5], *w[3]);
    circ.cx(*w[1], *w[5]);
    circ.cx(*w[3], *w[0]);
    circ.cx(*w[3], *w[5]);
    circ.ccx(*w[1], *w[3], *w[5]);
    circ.x(*w[2]);
    circ.cx(*w[4], *w[0]);
    circ.cx(*w[5], *w[1]);
    circ.x(*w[3]);
}

// ===================================================================
// Triple codec (9->7). pair(s0,s1) -> normalize -> merge25(s2). 18 executed
// Toffoli (forward) after the clear_and vents.
// ===================================================================

/// Data wires holding the 7-bit code after [`compress_3sym`].
pub const TRIPLE_DATA_WIRES: [usize; 7] = [0, 1, 2, 3, 4, 7, 8];
/// Wires freed to |0> after compress (pair's wire 5 + merge's wire 6).
pub const TRIPLE_FREED_WIRES: [usize; 2] = [5, 6];

/// Affine pair-code normalizer (s2,s3 pair-code set -> {0..24}). `(kind,a,b,c)`:
/// 0 = `x a`, 1 = `cx a->b`, 2 = `ccx a,b->c`. Wires addressed in the 6..11
/// reference layout; the 9->7 triple puts the same structure 6 wires lower.
#[rustfmt::skip]
const NORMALIZER_OPS: &[(u8, u8, u8, u8)] = &[
    (1,10,9,0), (1,9,6,0), (1,10,6,0), (1,6,10,0), (1,10,6,0), (0,8,0,0), (0,9,0,0), (2,7,9,10),
    (1,10,9,0), (1,10,7,0), (1,10,9,0), (1,9,10,0), (1,10,9,0), (1,8,10,0), (1,9,8,0), (1,8,9,0),
    (1,9,8,0), (1,8,7,0), (1,7,8,0), (1,8,7,0), (1,8,6,0), (1,6,8,0), (1,8,6,0), (2,6,8,10),
    (1,9,7,0), (1,10,9,0), (1,9,10,0), (1,10,9,0), (1,9,7,0), (1,7,9,0), (1,9,7,0), (1,7,6,0),
    (1,6,7,0), (1,7,6,0), (2,6,10,9), (1,10,9,0), (1,10,8,0), (1,10,7,0), (1,9,10,0), (1,9,8,0),
    (1,9,7,0), (1,10,8,0), (1,8,10,0), (1,10,8,0), (1,7,6,0), (1,6,10,0), (0,10,0,0), (2,6,7,8),
    (1,10,9,0), (1,10,8,0), (1,10,7,0), (1,10,6,0), (1,8,10,0), (1,8,9,0), (1,8,7,0), (1,8,6,0),
    (1,7,8,0), (1,6,9,0), (1,6,8,0), (0,8,0,0), (2,6,7,8), (1,10,9,0), (1,10,8,0), (1,10,7,0),
    (1,10,9,0), (1,9,10,0), (1,10,9,0), (1,8,10,0), (1,9,8,0), (1,8,9,0), (1,9,8,0), (1,7,6,0),
    (1,6,10,0), (1,6,9,0), (0,9,0,0), (0,10,0,0), (2,6,10,8), (1,10,9,0), (1,9,8,0), (1,9,7,0),
    (1,9,6,0), (1,8,7,0), (1,8,6,0), (1,10,8,0), (1,8,10,0), (1,10,8,0), (1,7,6,0), (1,6,9,0),
    (1,7,6,0), (1,6,7,0), (1,7,6,0), (0,6,0,0), (0,8,0,0), (2,7,8,9), (1,6,8,0), (1,7,8,0),
    (1,7,6,0), (1,6,8,0), (1,7,6,0), (1,6,7,0), (1,7,6,0), (0,6,0,0), (0,9,0,0), (0,10,0,0),
];

/// merge25: fold normalized pair23 (w[6..11]) + raw s4 (w[12..15]) into the high
/// code bits using clean ancilla w[15..17]. Same `(kind,a,b,c)` encoding.
#[rustfmt::skip]
const MERGE25_OPS: &[(u8, u8, u8, u8)] = &[
    (1,12,9,0), (1,14,10,0), (2,10,12,14), (1,13,9,0), (2,9,12,13), (2,13,14,12), (1,12,6,0), (1,7,10,0),
    (2,10,12,7), (1,6,9,0), (2,9,12,6), (1,12,8,0), (0,12,0,0), (1,14,12,0), (1,7,10,0), (0,10,0,0),
    (1,6,9,0), (0,13,0,0), (2,8,13,15), (2,14,15,16), (2,8,13,15), (2,10,9,15), (2,16,15,12), (2,10,9,15),
    (2,8,13,15), (2,14,15,16), (2,8,13,15), (0,13,0,0), (1,6,9,0), (0,10,0,0), (1,7,10,0),
];

// merge25 op indices whose CCX target is driven back to |0> -- the target holds
// (op.1 AND op.2) and the XOR clears it. These become `clear_and` vents. The
// set differs by direction.
const MERGE25_CLEAR_FWD: [usize; 5] = [20, 22, 23, 25, 26];
const MERGE25_CLEAR_REV: [usize; 4] = [18, 19, 21, 24];

/// Apply one normalizer/merge25 op against a wire slice with operand offset `off`
/// (0 = the 17-wire reference layout, 6 = the 9->7 triple's 11-wire layout).
/// Placeholder operands on x/cx ops are never read, so they need no remap.
#[inline]
fn apply_op_off(circ: &mut B, w: &[&QubitId], op: (u8, u8, u8, u8), off: u8) {
    let m = |i: u8| w[(i - off) as usize];
    match op.0 {
        0 => circ.x(*m(op.1)),
        1 => circ.cx(*m(op.1), *m(op.2)),
        2 => circ.ccx(*m(op.1), *m(op.2), *m(op.3)),
        _ => unreachable!("bad codec op kind"),
    }
}

/// Apply [`MERGE25_OPS`] (forward or reversed) with the AND-uncompute targets
/// vented via [`clear_and`]. `off` selects the wire layout.
fn apply_merge25(circ: &mut B, w: &[&QubitId], off: u8, reverse: bool) {
    let clear: &[usize] = if reverse { &MERGE25_CLEAR_REV } else { &MERGE25_CLEAR_FWD };
    let n = MERGE25_OPS.len();
    for step in 0..n {
        let i = if reverse { n - 1 - step } else { step };
        let op = MERGE25_OPS[i];
        if clear.contains(&i) {
            clear_and(
                circ,
                w[(op.3 - off) as usize],
                w[(op.1 - off) as usize],
                w[(op.2 - off) as usize],
            );
        } else {
            apply_op_off(circ, w, op, off);
        }
    }
}

/// Compress 3 base-5 symbols (s0 on `w[0..3]`, s1 on `w[3..6]`, s2 on `w[6..9]`)
/// plus 2 clean ancilla (`w[9]`, `w[10]`) into a 7-bit code on
/// [`TRIPLE_DATA_WIRES`]. [`TRIPLE_FREED_WIRES`] -> |0>, ancilla restored.
fn compress_3sym(circ: &mut B, w: &[&QubitId; 11]) {
    compress_2sym_fast(circ, &[w[0], w[1], w[2], w[3], w[4], w[5]]);
    for &op in NORMALIZER_OPS {
        apply_op_off(circ, &w[..], op, 6);
    }
    apply_merge25(circ, &w[..], 6, false);
}

/// Decompress: inverse of [`compress_3sym`] (ops reversed; X/CX/CCX self-inverse).
fn compress_3sym_reverse(circ: &mut B, w: &[&QubitId; 11]) {
    apply_merge25(circ, &w[..], 6, true);
    for &op in NORMALIZER_OPS.iter().rev() {
        apply_op_off(circ, &w[..], op, 6);
    }
    compress_2sym_fast_reverse(circ, &[w[0], w[1], w[2], w[3], w[4], w[5]]);
}

// ===================================================================
// Dialog codec abstraction + tiling. The product-min all-triple tiling at len 258
// emits Step0, Triple, Pair (Raw and the 5-symbol codec are defined but unreached
// for this config).
// ===================================================================

/// The per-region codecs the product-min tiling emits.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum DialogCodec {
    /// 2 symbols (6 raw bits) -> 5-bit code.
    Pair,
    /// 3 symbols (9 raw bits) -> 7-bit code. The product-min codec.
    Triple,
    /// 1 symbol kept raw (3 bits, no codec).
    Raw,
    /// Step-0 only: at the first GCD step `swap == subtracted` is known, so this
    /// 1-symbol codec keeps `subtracted` + `s_2` (2 bits) and drops the redundant
    /// `swap`: compress `CX(sub->swap)` clears swap to |0>; decompress restores it.
    Step0,
}

impl DialogCodec {
    /// Symbols per window.
    pub fn syms(self) -> usize {
        match self {
            Self::Pair => 2,
            Self::Triple => 3,
            Self::Raw | Self::Step0 => 1,
        }
    }
    /// Compressed code-bit width per window.
    pub fn code_bits(self) -> usize {
        match self {
            Self::Pair => 5,
            Self::Triple => 7,
            Self::Raw => 3,
            Self::Step0 => 2,
        }
    }
    /// Transient clean ancilla (|0> in and out), separate from the symbol wires.
    fn clean_anc(self) -> usize {
        match self {
            Self::Pair | Self::Raw | Self::Step0 => 0,
            Self::Triple => 2,
        }
    }
    /// Symbol-window wires holding the code after compress (kept by compaction).
    fn data_wires(self) -> &'static [usize] {
        match self {
            Self::Pair => &[0, 1, 2, 3, 4],
            Self::Triple => &TRIPLE_DATA_WIRES,
            Self::Raw => &[0, 1, 2],
            Self::Step0 => &[0, 2], // keep subtracted(0) + s_2(2); swap(1) is freed
        }
    }
    /// Symbol-window wires cleared to |0> by compress (freed by compaction).
    fn freed_wires(self) -> &'static [usize] {
        match self {
            Self::Pair => &[5],
            Self::Triple => &TRIPLE_FREED_WIRES,
            Self::Raw => &[],
            Self::Step0 => &[1], // swap (= subtracted; CX'd to |0>)
        }
    }
    /// Forward codec on a prebuilt window (compress: symbols -> code).
    fn compress(self, circ: &mut B, win: &[&QubitId]) {
        match self {
            Self::Pair => compress_2sym_fast(circ, win.try_into().unwrap()),
            Self::Triple => compress_3sym(circ, win.try_into().unwrap()),
            Self::Raw => {}
            // swap == subtracted -> CX(sub, swap) clears swap to |0>.
            Self::Step0 => circ.cx(*win[0], *win[1]),
        }
    }
    /// Reverse codec (decompress: code -> symbols).
    fn decompress(self, circ: &mut B, win: &[&QubitId]) {
        match self {
            Self::Pair => compress_2sym_fast_reverse(circ, win.try_into().unwrap()),
            Self::Triple => compress_3sym_reverse(circ, win.try_into().unwrap()),
            Self::Raw => {}
            // restore swap = subtracted on the freshly expanded |0> swap wire.
            Self::Step0 => circ.cx(*win[0], *win[1]),
        }
    }

    /// Incremental decompress of one window: given the window's `code_bits`
    /// compressed data qubits (`data`, in [`data_wires`] order), expand it back
    /// to the `syms*3` raw symbol slots (re-inserting freed |0> wires) and run
    /// the reverse codec. Returns the raw symbol slots in contiguous
    /// `[sym0.sub, sym0.swap, sym0.s2, sym1..]` order. The transient clean
    /// ancilla are allocated and freed inside.
    #[must_use]
    pub fn decompress_window(self, circ: &mut B, data: &[QubitId]) -> Vec<QubitId> {
        assert_eq!(data.len(), self.code_bits(), "data len != code_bits");
        // Build the syms*3 raw slot vector: data qubits at data_wires positions,
        // fresh |0> at freed_wires positions.
        let mut slots: Vec<Option<QubitId>> = (0..self.syms() * 3).map(|_| None).collect();
        let mut it = data.iter();
        for &d in self.data_wires() {
            slots[d] = Some(*it.next().expect("data bit"));
        }
        for &f in self.freed_wires() {
            slots[f] = Some(circ.alloc_qubit());
        }
        let raw: Vec<QubitId> = slots.into_iter().map(|s| s.expect("slot")).collect();
        let clean: Vec<QubitId> = (0..self.clean_anc()).map(|_| circ.alloc_qubit()).collect();
        let win: Vec<&QubitId> = raw.iter().chain(clean.iter()).collect();
        self.decompress(circ, &win);
        for q in clean {
            circ.zero_and_free(q);
        }
        raw
    }

    /// Incremental recompress of one window: inverse of [`decompress_window`].
    /// Given the `syms*3` raw symbol slots, run the forward codec, free the
    /// freed-wire slots (now |0>), and return the `code_bits` compressed data
    /// qubits (in [`data_wires`] order).
    #[must_use]
    pub fn compress_window(self, circ: &mut B, raw: &[QubitId]) -> Vec<QubitId> {
        assert_eq!(raw.len(), self.syms() * 3, "raw len != syms*3");
        let clean: Vec<QubitId> = (0..self.clean_anc()).map(|_| circ.alloc_qubit()).collect();
        let win: Vec<&QubitId> = raw.iter().chain(clean.iter()).collect();
        self.compress(circ, &win);
        for q in clean {
            circ.zero_and_free(q);
        }
        // Keep data_wires, free freed_wires (now |0>).
        let mut data: Vec<QubitId> = Vec::with_capacity(self.code_bits());
        let dset = self.data_wires();
        for (k, &q) in raw.iter().enumerate() {
            if dset.contains(&k) {
                data.push(q);
            } else {
                circ.zero_and_free(q);
            }
        }
        data
    }
}

/// Ordered region schedule for the all-triple product-min config: Step0 first,
/// then triples, then a small Pair/Raw tail. `(n3, n5)` are sized as the
/// product-min selector does (`n3 = iters/3`, `n5 = 0`); `iters` is the schedule
/// length (258 for ludicrous). Restricted to the codecs this config emits: the
/// 5-symbol codec and a dedicated raw tail are omitted (n5 == 0 and the raw tail
/// is fixed at 0 for product-min).
#[must_use]
pub fn jump_dialog_regions(n3: usize, iters: usize) -> Vec<(DialogCodec, usize)> {
    // Symbol 0 is the 1-symbol Step0 codec. The remaining iters-1 symbols tile.
    let codec_syms = iters - 1;
    let mut n3 = n3;
    while 3 * n3 > codec_syms {
        n3 -= 1;
    }
    let rem = codec_syms - 3 * n3;
    let mut r = vec![(DialogCodec::Step0, 1)];
    if n3 > 0 {
        r.push((DialogCodec::Triple, n3));
    }
    // Tail: a 3-symbol leftover packs denser as one Triple (tight regime: n3>0).
    let tight = n3 > 0;
    match rem {
        3 if tight => r.push((DialogCodec::Triple, 1)),
        _ => {
            if rem / 2 > 0 {
                r.push((DialogCodec::Pair, rem / 2));
            }
            if rem % 2 == 1 {
                r.push((DialogCodec::Raw, 1));
            }
        }
    }
    r
}

/// Persistent compressed-tape qubits for the `(n3)` all-triple dialog codec at
/// `iters` symbols: `1` (t1) + sum over regions of `code_bits * count`.
#[must_use]
pub fn dialog_tape_qubits(n3: usize, iters: usize) -> usize {
    1 + jump_dialog_regions(n3, iters)
        .into_iter()
        .map(|(codec, count)| codec.code_bits() * count)
        .sum::<usize>()
}

