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

use super::{BExt, B};
use crate::circuit::QubitId;

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

const TAIL4_TOP32_CODE_BITS: usize = 5;
const TAIL4_TOP32_CODE_CONSTANT: u8 = 22;
const TAIL4_TOP32_ENCODER_ANF: [&[u16]; TAIL4_TOP32_CODE_BITS] = [
    &[2, 4, 8, 16, 32, 1024, 65, 528],
    &[1, 2, 4, 16, 32, 520],
    &[1, 2, 4, 8, 1024, 72],
    &[4, 16, 32, 256, 520],
    &[1, 4, 1024, 10, 66],
];
const TAIL4_TOP32_DECODER_ANF: [&[u16]; 12] = [
    &[0, 3, 7, 18, 20, 21, 25, 27, 29, 31],
    &[5, 6, 7, 9, 11, 13, 16, 18, 19, 20, 21, 31],
    &[0, 11, 15, 16, 18, 20, 23, 24, 26, 27, 28, 29, 31],
    &[3, 4, 5, 6, 11, 18, 19, 20, 21, 24, 26, 27, 28, 30],
    &[0, 2, 3, 5, 6, 9, 13, 18, 19, 20, 23, 25, 31],
    &[30, 31],
    &[3, 7, 9, 13, 19, 21, 25, 30],
    &[],
    &[2, 3, 5, 6, 8, 9, 11, 13, 16, 19, 25, 27, 29],
    &[0, 1, 4, 9, 13, 18, 19, 20, 27, 29],
    &[0, 3, 7, 9, 13, 19, 21, 25, 30],
    &[0],
];

fn tail4_top32_enabled() -> bool {
    std::env::var("TLM_TAIL4_TOP32").ok().as_deref() == Some("1")
}

fn toggle_mcx_with_dirty(circ: &mut B, controls: &[QubitId], dirty: &[QubitId], target: QubitId) {
    debug_assert!(!controls.contains(&target));
    match controls.len() {
        0 => circ.x(target),
        1 => circ.cx(controls[0], target),
        2 => circ.ccx(controls[0], controls[1], target),
        count => {
            let bridge = dirty
                .iter()
                .copied()
                .find(|q| *q != target && !controls.contains(q))
                .expect("tail4 codec needs a disjoint dirty bridge");
            let rest: Vec<QubitId> = dirty.iter().copied().filter(|q| *q != bridge).collect();
            toggle_mcx_with_dirty(circ, &controls[..count - 1], &rest, bridge);
            circ.ccx(bridge, controls[count - 1], target);
            toggle_mcx_with_dirty(circ, &controls[..count - 1], &rest, bridge);
            circ.ccx(bridge, controls[count - 1], target);
        }
    }
}

fn toggle_anf_with_dirty(
    circ: &mut B,
    controls: &[QubitId],
    target: QubitId,
    dirty: &[QubitId],
    terms: &[u16],
) {
    for &mask in terms {
        let term_controls: Vec<QubitId> = controls
            .iter()
            .enumerate()
            .filter_map(|(i, q)| ((mask >> i) & 1 != 0).then_some(*q))
            .collect();
        toggle_mcx_with_dirty(circ, &term_controls, dirty, target);
    }
}

// The historical K5 codec stores each symbol's first two flags contiguously,
// followed by all four s2 flags. The current tape stores contiguous triples.
fn tail4_reordered_raw(raw: &[QubitId]) -> [QubitId; 12] {
    assert_eq!(raw.len(), 12, "tail4 raw window must contain four symbols");
    [
        raw[0], raw[1], raw[3], raw[4], raw[6], raw[7], raw[9], raw[10], raw[2], raw[5], raw[8],
        raw[11],
    ]
}

fn tail4_toggle_code_from_raw(circ: &mut B, code: &[QubitId], raw: &[QubitId]) {
    assert_eq!(code.len(), TAIL4_TOP32_CODE_BITS);
    let wires = tail4_reordered_raw(raw);
    for (i, terms) in TAIL4_TOP32_ENCODER_ANF.iter().enumerate() {
        if (TAIL4_TOP32_CODE_CONSTANT >> i) & 1 != 0 {
            circ.x(code[i]);
        }
        for &mask in *terms {
            let controls: Vec<QubitId> = wires
                .iter()
                .enumerate()
                .filter_map(|(j, q)| ((mask >> j) & 1 != 0).then_some(*q))
                .collect();
            toggle_mcx_with_dirty(circ, &controls, &wires, code[i]);
        }
    }
}

fn tail4_toggle_raw_from_code(circ: &mut B, code: &[QubitId], raw: &[QubitId]) {
    assert_eq!(code.len(), TAIL4_TOP32_CODE_BITS);
    let wires = tail4_reordered_raw(raw);
    for (i, terms) in TAIL4_TOP32_DECODER_ANF.iter().enumerate() {
        toggle_anf_with_dirty(circ, code, wires[i], &wires, terms);
    }
}

fn compress_tail4_top32_payload(circ: &mut B, raw: &[QubitId]) -> Vec<QubitId> {
    let code: Vec<QubitId> = (0..TAIL4_TOP32_CODE_BITS)
        .map(|_| circ.alloc_qubit())
        .collect();
    tail4_toggle_code_from_raw(circ, &code, raw);
    tail4_toggle_raw_from_code(circ, &code, raw);
    for &q in raw {
        circ.zero_and_free(q);
    }
    code
}

fn decompress_tail4_top32_payload(circ: &mut B, code: &[QubitId]) -> Vec<QubitId> {
    let raw: Vec<QubitId> = (0..12).map(|_| circ.alloc_qubit()).collect();
    tail4_toggle_raw_from_code(circ, code, &raw);
    tail4_toggle_code_from_raw(circ, code, &raw);
    for &q in code {
        circ.zero_and_free(q);
    }
    raw
}

fn compress_tail4_top32(circ: &mut B, raw: &[QubitId]) -> Vec<QubitId> {
    assert_eq!(
        raw.len(),
        15,
        "tail4 hybrid window must contain five symbols"
    );
    let mut data = raw[..3].to_vec();
    data.extend(compress_tail4_top32_payload(circ, &raw[3..]));
    data
}

fn decompress_tail4_top32(circ: &mut B, data: &[QubitId]) -> Vec<QubitId> {
    assert_eq!(data.len(), 3 + TAIL4_TOP32_CODE_BITS);
    let mut raw = data[..3].to_vec();
    raw.extend(decompress_tail4_top32_payload(circ, &data[3..]));
    raw
}

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
    let clear: &[usize] = if reverse {
        &MERGE25_CLEAR_REV
    } else {
        &MERGE25_CLEAR_FWD
    };
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
    /// One generic symbol plus four terminal symbols on the frozen support (15 -> 8).
    Tail4Top32,
}

impl DialogCodec {
    /// Symbols per window.
    pub fn syms(self) -> usize {
        match self {
            Self::Pair => 2,
            Self::Triple => 3,
            Self::Tail4Top32 => 5,
            Self::Raw | Self::Step0 => 1,
        }
    }
    /// Compressed code-bit width per window.
    pub fn code_bits(self) -> usize {
        match self {
            Self::Pair => 5,
            Self::Triple => 7,
            Self::Tail4Top32 => 3 + TAIL4_TOP32_CODE_BITS,
            Self::Raw => 3,
            Self::Step0 => 2,
        }
    }
    /// Transient clean ancilla (|0> in and out), separate from the symbol wires.
    fn clean_anc(self) -> usize {
        match self {
            Self::Pair | Self::Raw | Self::Step0 | Self::Tail4Top32 => 0,
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
            Self::Tail4Top32 => &[], // separate code wires; handled by the special path
        }
    }
    /// Symbol-window wires cleared to |0> by compress (freed by compaction).
    fn freed_wires(self) -> &'static [usize] {
        match self {
            Self::Pair => &[5],
            Self::Triple => &TRIPLE_FREED_WIRES,
            Self::Raw => &[],
            Self::Step0 => &[1], // swap (= subtracted; CX'd to |0>)
            Self::Tail4Top32 => &[],
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
            Self::Tail4Top32 => unreachable!("tail4 uses separate code wires"),
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
            Self::Tail4Top32 => unreachable!("tail4 uses separate code wires"),
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
        if self == Self::Tail4Top32 {
            return decompress_tail4_top32(circ, data);
        }
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
        if self == Self::Tail4Top32 {
            return compress_tail4_top32(circ, raw);
        }
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

/// Step-0-only codec that folds the separate `t1` prefix into the first dialog
/// symbol. The reachable triples `(t1, subtracted, s2)` are only four cases:
/// `(0,1,0)`, `(1,1,0)`, `(1,1,1)`, `(1,0,1)`. They are not an affine subspace,
/// so one nonlinear gate is needed to clear a wire while keeping a reversible map.
#[must_use]
pub fn compress_step0_with_t1(circ: &mut B, t1: QubitId, raw: &[QubitId]) -> Vec<QubitId> {
    assert_eq!(raw.len(), 3, "step0 raw symbol is [sub, swap, s2]");
    let sub = raw[0];
    let swap = raw[1];
    let s2 = raw[2];
    circ.cx(sub, swap);
    circ.cx(sub, t1);
    circ.x(sub);
    circ.ccx(t1, s2, sub);
    circ.zero_and_free(sub);
    circ.zero_and_free(swap);
    vec![t1, s2]
}

/// Inverse of [`compress_step0_with_t1`]. Returns the restored `t1` qubit plus
/// the raw `[sub, swap, s2]` symbol slots.
#[must_use]
pub fn decompress_step0_with_t1(circ: &mut B, data: &[QubitId]) -> (QubitId, Vec<QubitId>) {
    assert_eq!(data.len(), 2, "step0+t1 code is two bits");
    let t1 = data[0];
    let s2 = data[1];
    let sub = circ.alloc_qubit();
    let swap = circ.alloc_qubit();
    circ.ccx(t1, s2, sub);
    circ.x(sub);
    circ.cx(sub, t1);
    circ.cx(sub, swap);
    (t1, vec![sub, swap, s2])
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
    let tail4 = usize::from(tail4_top32_enabled() && iters >= 6) * 5;
    let codec_syms = iters - 1 - tail4;
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
    if tail4 != 0 {
        r.push((DialogCodec::Tail4Top32, 1));
    }
    r
}

/// Persistent compressed-tape qubits for the `(n3)` all-triple dialog codec at
/// `iters` symbols: sum over regions of `code_bits * count`. Step0's code bits
/// include the first-shift `t1` bit, so no separate prefix is resident.
#[must_use]
pub fn dialog_tape_qubits(n3: usize, iters: usize) -> usize {
    jump_dialog_regions(n3, iters)
        .into_iter()
        .map(|(codec, count)| codec.code_bits() * count)
        .sum::<usize>()
}
