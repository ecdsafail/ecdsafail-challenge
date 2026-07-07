
use super::{B, BExt};
use crate::circuit::{QubitId};

fn clear_and(circ: &mut B, t: &QubitId, a: &QubitId, b: &QubitId) {
    let bit = circ.alloc_bit();
    circ.hmr(*t, bit);
    circ.cz_if_bit(*a, *b, bit);
}

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

    clear_and(circ, w[5], w[3], w[4]);
}

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

pub const TRIPLE_DATA_WIRES: [usize; 7] = [0, 1, 2, 3, 4, 7, 8];

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

fn toggle_mcx_with_dirty(
    circ: &mut B,
    controls: &[QubitId],
    dirty: &[QubitId],
    target: QubitId,
) {
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

fn tail4_reordered_raw(raw: &[QubitId]) -> [QubitId; 12] {
    assert_eq!(raw.len(), 12, "tail4 raw window must contain four symbols");
    [
        raw[0], raw[1], raw[3], raw[4], raw[6], raw[7], raw[9], raw[10],
        raw[2], raw[5], raw[8], raw[11],
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
    assert_eq!(raw.len(), 15, "tail4 hybrid window must contain five symbols");
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

#[rustfmt::skip]
const MERGE25_OPS: &[(u8, u8, u8, u8)] = &[
    (1,12,9,0), (1,14,10,0), (2,10,12,14), (1,13,9,0), (2,9,12,13), (2,13,14,12), (1,12,6,0), (1,7,10,0),
    (2,10,12,7), (1,6,9,0), (2,9,12,6), (1,12,8,0), (0,12,0,0), (1,14,12,0), (1,7,10,0), (0,10,0,0),
    (1,6,9,0), (0,13,0,0), (2,8,13,15), (2,14,15,16), (2,8,13,15), (2,10,9,15), (2,16,15,12), (2,10,9,15),
    (2,8,13,15), (2,14,15,16), (2,8,13,15), (0,13,0,0), (1,6,9,0), (0,10,0,0), (1,7,10,0),
];

const MERGE25_CLEAR_FWD: [usize; 5] = [20, 22, 23, 25, 26];
const MERGE25_CLEAR_REV: [usize; 4] = [18, 19, 21, 24];

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

fn compress_3sym(circ: &mut B, w: &[&QubitId; 11]) {
    compress_2sym_fast(circ, &[w[0], w[1], w[2], w[3], w[4], w[5]]);
    for &op in NORMALIZER_OPS {
        apply_op_off(circ, &w[..], op, 6);
    }
    apply_merge25(circ, &w[..], 6, false);
}

fn compress_3sym_reverse(circ: &mut B, w: &[&QubitId; 11]) {
    apply_merge25(circ, &w[..], 6, true);
    for &op in NORMALIZER_OPS.iter().rev() {
        apply_op_off(circ, &w[..], op, 6);
    }
    compress_2sym_fast_reverse(circ, &[w[0], w[1], w[2], w[3], w[4], w[5]]);
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum DialogCodec {

    Pair,

    Triple,

    Raw,

    Step0,

    Tail4Top32,
}

impl DialogCodec {

    pub fn syms(self) -> usize {
        match self {
            Self::Pair => 2,
            Self::Triple => 3,
            Self::Tail4Top32 => 5,
            Self::Raw | Self::Step0 => 1,
        }
    }

    pub fn code_bits(self) -> usize {
        match self {
            Self::Pair => 5,
            Self::Triple => 7,
            Self::Tail4Top32 => 3 + TAIL4_TOP32_CODE_BITS,
            Self::Raw => 3,
            Self::Step0 => 2,
        }
    }

    fn clean_anc(self) -> usize {
        match self {
            Self::Pair | Self::Raw | Self::Step0 | Self::Tail4Top32 => 0,
            Self::Triple => 2,
        }
    }

    fn data_wires(self) -> &'static [usize] {
        match self {
            Self::Pair => &[0, 1, 2, 3, 4],
            Self::Triple => &TRIPLE_DATA_WIRES,
            Self::Raw => &[0, 1, 2],
            Self::Step0 => &[0, 2],
            Self::Tail4Top32 => &[],
        }
    }

    fn freed_wires(self) -> &'static [usize] {
        match self {
            Self::Pair => &[5],
            Self::Triple => &TRIPLE_FREED_WIRES,
            Self::Raw => &[],
            Self::Step0 => &[1],
            Self::Tail4Top32 => &[],
        }
    }

    fn compress(self, circ: &mut B, win: &[&QubitId]) {
        match self {
            Self::Pair => compress_2sym_fast(circ, win.try_into().unwrap()),
            Self::Triple => compress_3sym(circ, win.try_into().unwrap()),
            Self::Raw => {}

            Self::Step0 => circ.cx(*win[0], *win[1]),
            Self::Tail4Top32 => unreachable!("tail4 uses separate code wires"),
        }
    }

    fn decompress(self, circ: &mut B, win: &[&QubitId]) {
        match self {
            Self::Pair => compress_2sym_fast_reverse(circ, win.try_into().unwrap()),
            Self::Triple => compress_3sym_reverse(circ, win.try_into().unwrap()),
            Self::Raw => {}

            Self::Step0 => circ.cx(*win[0], *win[1]),
            Self::Tail4Top32 => unreachable!("tail4 uses separate code wires"),
        }
    }

    #[must_use]
    pub fn decompress_window(self, circ: &mut B, data: &[QubitId]) -> Vec<QubitId> {
        assert_eq!(data.len(), self.code_bits(), "data len != code_bits");
        if self == Self::Tail4Top32 {
            return decompress_tail4_top32(circ, data);
        }

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

#[must_use]
pub fn jump_dialog_regions(n3: usize, iters: usize) -> Vec<(DialogCodec, usize)> {

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

#[must_use]
pub fn dialog_tape_qubits(n3: usize, iters: usize) -> usize {
    jump_dialog_regions(n3, iters)
        .into_iter()
        .map(|(codec, count)| codec.code_bits() * count)
        .sum::<usize>()
}
