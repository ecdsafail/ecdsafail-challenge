//! Base-5 dialog compressor for the jump-before-swap GCD (jump=2): packs 3
//! successive 5-symbol windows into a 7-bit radix-5 value (5^3=125 < 2^7=128),
//! ~2.33 bits/step vs the raw 3. In-place radix merge into a SEPARATE 7-bit
//! accumulator (no digit-0 overlap, no per-digit move-to-scratch).
//!
//! Each symbol is (sub, swap, s2) = w[3i], w[3i+1], w[3i+2] with digit
//!   `d_i` = `sym_to_digit(sub,swap,s2)`: (1,0,0)=0,(1,1,0)=1,(1,0,1)=2,(1,1,1)=3,
//!   (0,0,1)=4 (overflow).
//! Merge adds `d_i`*5^i to e via `e += swap*r + s2*2r + (!sub)*2r` (3 controlled
//! const-adds, no Toffoli on the controls). Recovery clears the symbol bits
//! from e: s2 = (e>=2r); sub = !(e>=4r); swap = parity (e>=r)^(e>=2r)^(e>=3r)^
//! (e>=4r). After all 3 merges w[0..9] are |0> and e holds the code; e is
//! swapped into w[0..7], leaving w[7],w[8] = |0> (the 2 freed bits).

use crate::point_add::trailmix_port::arith::gidney_const_adder::{
    compare_geq_const_gidney_refs, controlled_add_const_gidney_refs,
};
use crate::point_add::trailmix_port::arith::schrottenloher::gcd_jump::jump_bs_digit_to_sym;
use crate::point_add::trailmix_port::circuit::{Circuit, QReg};
use crate::point_add::trailmix_port::tracker::ghost::Ghost;

const POW5: [u8; 3] = [1, 5, 25];
/// accumulator width after merging digit i (value < 5^{i+1}): 5,25,125.
const ACC5: [usize; 3] = [3, 5, 7];

/// number of valid base-5 codes (5^3); codes 125,126,127 never occur.
const N_CODES: usize = 125;

/// The 9-bit symbol pattern produced by code `c` (0..125): bit 3i+0 = `sub_i`,
/// 3i+1 = `swap_i`, 3i+2 = `s2_i`, where (`sub_i,swap_i,s2_i`) = `digit_to_sym(c`'s
/// i-th base-5 digit). This is the QROM table value: `v == code_pattern(code)`
/// for a freshly-built code, so XOR-ing it clears v (and writing it expands).
fn code_pattern(c: usize) -> u16 {
    let mut p = 0u16;
    let mut c = c;
    for i in 0..3 {
        let d = (c % 5) as u8;
        c /= 5;
        let (sub, swap, s2) = jump_bs_digit_to_sym(d);
        p |= u16::from(sub) << (3 * i);
        p |= u16::from(swap) << (3 * i + 1);
        p |= u16::from(s2) << (3 * i + 2);
    }
    p
}

/// Build the dense 7-bit base-5 code into accumulator `e` from the 3 symbols in
/// `w[0..9]` (cheap controlled-const-adds, NO recovery comparators): for each
/// digit i, e += `swap_i`*5^i + `s2_i`*2*5^i + (!`sub_i`)*2*5^i = `digit_i` * 5^i.
/// Leaves w unchanged (the sub-flip cancels).
fn build_code(circ: &mut Circuit, w: &[&QReg; 9], e: &[QReg], dirty: &[&QReg]) {
    // digit 0: e is |0>, so d0 = swap + 2(s2 & sub) + 4(!sub) (the bit form,
    // equal to the value form swap + 2*s2 + 2*!sub on every valid symbol)
    // writes DIRECTLY into e[0..3] with no carry -- one CCX for bit 1, CX for
    // the rest. The full const-adders into a zero accumulator are wasted work.
    {
        let (sub, swap, s2) = (w[0], w[1], w[2]);
        circ.cx(swap, &e[0]); // e[0] = swap
        circ.ccx(s2, sub, &e[1]); // e[1] = s2 & sub
        circ.x(&e[2]);
        circ.cx(sub, &e[2]); // e[2] = !sub
    }
    for i in 1..3 {
        let r = POW5[i];
        let ei: Vec<&QReg> = (0..ACC5[i]).map(|k| &e[k]).collect();
        let (sub, swap, s2) = (w[3 * i], w[3 * i + 1], w[3 * i + 2]);
        controlled_add_const_gidney_refs(circ, swap, &ei, &[r], dirty);
        controlled_add_const_gidney_refs(circ, s2, &ei, &[2 * r], dirty);
        circ.x(sub);
        controlled_add_const_gidney_refs(circ, sub, &ei, &[2 * r], dirty);
        circ.x(sub);
    }
}

/// Inverse of `build_code`: drain `e` back to |0> using the symbols in `w`
/// (two's-complement controlled-subtracts).
fn unbuild_code(circ: &mut Circuit, w: &[&QReg; 9], e: &[QReg], dirty: &[&QReg]) {
    for i in (0..3).rev() {
        let r = POW5[i];
        let wi = ACC5[i];
        let ei: Vec<&QReg> = (0..wi).map(|k| &e[k]).collect();
        let (sub, swap, s2) = (w[3 * i], w[3 * i + 1], w[3 * i + 2]);
        let m = 1u16 << wi;
        let c_r = (m - u16::from(r)).to_le_bytes();
        let c_2r = (m - u16::from(2 * r)).to_le_bytes();
        circ.x(sub);
        controlled_add_const_gidney_refs(circ, sub, &ei, &c_2r, dirty);
        circ.x(sub);
        controlled_add_const_gidney_refs(circ, s2, &ei, &c_2r, dirty);
        controlled_add_const_gidney_refs(circ, swap, &ei, &c_r, dirty);
    }
}

/// XOR the symbol pattern for the code held in `e` into `w[0..9]` via ONE log*
/// unary-iterate keyed on the dense code: at step `code`, the one-hot `gate`
/// fires iff e==code, and we CX `code_pattern(code)` onto w (Clifford, free).
/// Self-inverse; used to clear w (pack) or to write w (unpack).
/// Babbush/Gidney unary iteration over codes `[base, base+2^bits.len())` keyed
/// on the address bits `bits` (LSB-first within the code) with the parent
/// selector `ctrl`: 1 Toffoli per internal node (the CNOT-sibling-reuse trick
/// reflects `ctrl&top` to `ctrl&!top` for free) and a measurement-free
/// `clear_and` uncompute, so a 7-bit address costs ~2^7 = 127 executed Toffoli.
/// `body(circ, code, ctrl)` runs at each in-range leaf with `ctrl = (e==code)`.
/// Top-level entry: the MSB of the address selects the two halves DIRECTLY
/// (no constant-1 root ctrl -- `one & top` would declare as `CopyOf(top)` but
/// discharge as `AndOf(one, top)`, which the structural matcher can't reconcile
/// since it doesn't know `one == 1`). High half uses `top` as ctrl; low half
/// flips `top` to act as `!top`.
fn unary_iter_codes<F: FnMut(&mut Circuit, usize, &QReg)>(
    circ: &mut Circuit,
    addr: &[&QReg],
    body: &mut F,
) {
    debug_assert!(!addr.is_empty(), "unary_iter_codes: need >= 1 address bit");
    let top = addr[addr.len() - 1];
    let rest = &addr[..addr.len() - 1];
    let half = 1usize << rest.len();
    unary_iter_tree(circ, top, rest, half, body); // MSB = 1
    circ.x(top);
    unary_iter_tree(circ, top, rest, 0, body); // MSB = 0
    circ.x(top);
}

fn unary_iter_tree<F: FnMut(&mut Circuit, usize, &QReg)>(
    circ: &mut Circuit,
    ctrl: &QReg,
    bits: &[&QReg],
    base: usize,
    body: &mut F,
) {
    if base >= N_CODES {
        return; // whole subtree out of range; never selected (e < N_CODES)
    }
    if bits.is_empty() {
        body(circ, base, ctrl);
        return;
    }
    let top = bits[bits.len() - 1];
    let rest = &bits[..bits.len() - 1];
    let half = 1usize << rest.len();
    let anc = circ.alloc_qreg("qrom_anc");
    circ.ccx(ctrl, top, &anc); // anc = ctrl AND top
    unary_iter_tree(circ, &anc, rest, base + half, body); // top=1 (high half)
    circ.cx(ctrl, &anc); // anc = ctrl AND NOT top (free reflect)
    unary_iter_tree(circ, &anc, rest, base, body); // top=0 (low half)
    circ.cx(ctrl, &anc); // restore anc = ctrl AND top
    circ.clear_and(&anc, ctrl, top); // measurement-free uncompute
    circ.zero_and_free(anc);
}

/// Discharge variant of the unary tree for the ghost-CLEAR (phase-only): the
/// leaf selector `AND(ctrl, last_bit)` is never materialized -- since the ghost
/// bit is CLASSICAL, `Z(AND(ctrl,last_bit))` gated on it IS `cz_if_bit(ctrl,
/// last_bit, ghost)`, a free 2-qubit CZ. So the entire bottom level of ANDs
/// (64 of the ~126 ccx) vanishes. `deposit(circ, code, a, b)` does the per-code
/// discharge as CZ on the two AND-inputs (a, b). Used ONLY for the clear (Z's);
/// the decompress CX-write still needs the materialized ctrl (a CX, not a phase,
/// and multi-target so the shared selector wins).
fn discharge_tree<F: FnMut(&mut Circuit, usize, &QReg, &QReg)>(
    circ: &mut Circuit,
    ctrl: &QReg,
    bits: &[&QReg],
    base: usize,
    deposit: &mut F,
) {
    if base >= N_CODES {
        return;
    }
    if bits.len() == 1 {
        // deepest level: discharge both leaves via CZ on (ctrl, bit0) -- no AND.
        let bit0 = bits[0];
        if base + 1 < N_CODES {
            deposit(circ, base + 1, ctrl, bit0); // bit0 = 1
        }
        circ.x(bit0);
        deposit(circ, base, ctrl, bit0); // bit0 = 0 (CZ on !bit0)
        circ.x(bit0);
        return;
    }
    let top = bits[bits.len() - 1];
    let rest = &bits[..bits.len() - 1];
    let half = 1usize << rest.len();
    let anc = circ.alloc_qreg("qrom_anc");
    circ.ccx(ctrl, top, &anc);
    discharge_tree(circ, &anc, rest, base + half, deposit);
    circ.cx(ctrl, &anc);
    discharge_tree(circ, &anc, rest, base, deposit);
    circ.cx(ctrl, &anc);
    circ.clear_and(&anc, ctrl, top);
    circ.zero_and_free(anc);
}

/// Root wrapper for [`discharge_tree`] (MSB selects the halves directly).
fn discharge_codes<F: FnMut(&mut Circuit, usize, &QReg, &QReg)>(
    circ: &mut Circuit,
    addr: &[&QReg],
    deposit: &mut F,
) {
    let top = addr[addr.len() - 1];
    let rest = &addr[..addr.len() - 1];
    let half = 1usize << rest.len();
    discharge_tree(circ, top, rest, half, deposit);
    circ.x(top);
    discharge_tree(circ, top, rest, 0, deposit);
    circ.x(top);
}

/// QROM pack (9->7): build the dense base-5 code (the 9->7 arithmetic), swap it
/// into `w[0..7]` (so the output code lands there and the 9 symbol bits are
/// displaced into `e[0..7] ++ w[7..9]`), then GHOST-clear those 9 displaced
/// bits keyed on the code: each is HMR'd out (the tracker can't prove the
/// would-be CX clear is |0>, but a measured-out bit is freed), and its phase is
/// discharged by depositing `Z(leaf_ctrl)` into the ghost at every code whose
/// symbol pattern has that bit set (the unary tree). `close_ghost` asserts the
/// deposited terms XOR to the measured value -- the "`resolved_ghost`" check.
/// Post: `w[0..7]` = code, `w[7]`,`w[8]` = |0> (the caller's freed pair).
pub fn compress_3sym_qrom_refs(circ: &mut Circuit, w: &[&QReg; 9], dirty: &[&QReg]) {
    let prev = circ.push_section("q5enc");
    let e: Vec<QReg> = (0..7).map(|_| circ.alloc_qreg("q5_acc")).collect();
    build_code(circ, w, &e, dirty); // e = code, w = symbols
    for k in 0..7 {
        circ.swap(w[k], &e[k]); // w[0..7] = code; symbols now in e[0..7] ++ w[7..9]
    }
    // ghost-clear the 9 displaced symbol bits, keyed on the code in w[0..7].
    let mut ghosts: Vec<Ghost> = Vec::with_capacity(9);
    for j in 0..9 {
        let slot: &QReg = if j < 7 { &e[j] } else { w[j] };
        ghosts.push(circ.hmr_ghost(slot));
    }
    let addr: Vec<&QReg> = (0..7).map(|k| w[k]).collect();
    discharge_codes(circ, &addr, &mut |circ, code, a, b| {
        let pat = code_pattern(code);
        for bit in 0..9 {
            if (pat >> bit) & 1 == 1 {
                // Z(AND(a,b)) gated on the classical ghost bit = CZ(a,b) -- free.
                circ.ghost_xor_cz(&mut ghosts[bit], a, b);
            }
        }
    });
    for g in ghosts {
        circ.close_ghost(g); // assert the deposited CZ terms == measured value
    }
    for q in e.into_iter().rev() {
        circ.zero_and_free(q); // e[0..7] HMR'd to |0>
    }
    circ.pop_section(&prev);
}

/// Inverse of `compress_3sym_qrom_refs`: code in `w[0..7]` (`w[7..9]`=|0>) ->
/// 3 symbols. The clean mirror -- CX-write the symbol pattern into the |0>
/// slots `e[0..7] ++ w[7..9]` keyed on the code (writing into |0> IS provable,
/// so no ghosts), swap the code back out of `w[0..7]`, then drain `e`.
pub fn compress_3sym_qrom_reverse_refs(circ: &mut Circuit, w: &[&QReg; 9], dirty: &[&QReg]) {
    let prev = circ.push_section("q5dec");
    let e: Vec<QReg> = (0..7).map(|_| circ.alloc_qreg("q5_acc")).collect();
    let addr: Vec<&QReg> = (0..7).map(|k| w[k]).collect();
    unary_iter_codes(circ, &addr, &mut |circ, code, ctrl| {
        let pat = code_pattern(code);
        for b in 0..9 {
            if (pat >> b) & 1 == 1 {
                let slot: &QReg = if b < 7 { &e[b] } else { w[b] };
                circ.cx(ctrl, slot); // |0> -> symbol bit (provable write)
            }
        }
    });
    // now e[0..7] ++ w[7..9] = symbols, w[0..7] = code.
    for k in 0..7 {
        circ.swap(w[k], &e[k]); // w[0..7] = symbols[0..7]; e = code
    }
    unbuild_code(circ, w, &e, dirty); // e -= base-5 code -> |0>
    for q in e.into_iter().rev() {
        circ.zero_and_free(q);
    }
    circ.pop_section(&prev);
}

/// Forward merge of digit `i` (symbol bits w[3i..3i+3]) into accumulator `e`,
/// then recover/clear the symbol bits from `e`.
fn merge_digit5(circ: &mut Circuit, w: &[&QReg; 9], e: &[QReg], i: usize, dirty: &[&QReg]) {
    let r = POW5[i];
    let wi = ACC5[i];
    let ei: Vec<&QReg> = (0..wi).map(|k| &e[k]).collect();
    let sub = w[3 * i];
    let swap = w[3 * i + 1];
    let s2 = w[3 * i + 2];
    // e += swap*r + s2*2r + (!sub)*2r
    controlled_add_const_gidney_refs(circ, swap, &ei, &[r], dirty);
    controlled_add_const_gidney_refs(circ, s2, &ei, &[2 * r], dirty);
    circ.x(sub);
    controlled_add_const_gidney_refs(circ, sub, &ei, &[2 * r], dirty);
    circ.x(sub);
    // recover/clear: s2 ^= (e>=2r); sub: x then ^= (e>=4r); swap ^= parity.
    compare_geq_const_gidney_refs(circ, &ei, &[2 * r], s2, dirty);
    circ.x(sub);
    compare_geq_const_gidney_refs(circ, &ei, &[4 * r], sub, dirty);
    compare_geq_const_gidney_refs(circ, &ei, &[r], swap, dirty);
    compare_geq_const_gidney_refs(circ, &ei, &[2 * r], swap, dirty);
    compare_geq_const_gidney_refs(circ, &ei, &[3 * r], swap, dirty);
    compare_geq_const_gidney_refs(circ, &ei, &[4 * r], swap, dirty);
}

/// Inverse of `merge_digit5`: re-derive the symbol bits from `e`, then subtract
/// `d_i`*5^i (two's-complement controlled add of 2^wi - k).
fn unmerge_digit5(circ: &mut Circuit, w: &[&QReg; 9], e: &[QReg], i: usize, dirty: &[&QReg]) {
    let r = POW5[i];
    let wi = ACC5[i];
    let ei: Vec<&QReg> = (0..wi).map(|k| &e[k]).collect();
    let sub = w[3 * i];
    let swap = w[3 * i + 1];
    let s2 = w[3 * i + 2];
    // un-recover (reverse order of the 6 comparators).
    compare_geq_const_gidney_refs(circ, &ei, &[4 * r], swap, dirty);
    compare_geq_const_gidney_refs(circ, &ei, &[3 * r], swap, dirty);
    compare_geq_const_gidney_refs(circ, &ei, &[2 * r], swap, dirty);
    compare_geq_const_gidney_refs(circ, &ei, &[r], swap, dirty);
    compare_geq_const_gidney_refs(circ, &ei, &[4 * r], sub, dirty);
    circ.x(sub);
    compare_geq_const_gidney_refs(circ, &ei, &[2 * r], s2, dirty);
    // un-merge: e -= swap*r + s2*2r + (!sub)*2r via two's-complement adds.
    let m = 1u16 << wi;
    let c_r = (m - u16::from(r)).to_le_bytes();
    let c_2r = (m - u16::from(2 * r)).to_le_bytes();
    circ.x(sub);
    controlled_add_const_gidney_refs(circ, sub, &ei, &c_2r, dirty);
    circ.x(sub);
    controlled_add_const_gidney_refs(circ, s2, &ei, &c_2r, dirty);
    controlled_add_const_gidney_refs(circ, swap, &ei, &c_r, dirty);
}

/// Compress 3 successive jump-before-swap symbols (`w[0..9]`, each
/// (sub,swap,s2)) into the 7-bit base-5 code in `w[0..7]`; `w[7]`,`w[8]` are
/// |0> afterward. `dirty`: >= 13 BORROWED bits disjoint from `w`, restored.
pub fn compress_3sym_base5_refs(circ: &mut Circuit, w: &[&QReg; 9], dirty: &[&QReg]) {
    let prev = circ.push_section("b5pack");
    let e: Vec<QReg> = (0..7).map(|_| circ.alloc_qreg("b5_acc")).collect();
    for i in 0..3 {
        merge_digit5(circ, w, &e, i, dirty);
    }
    // w[0..9] are |0>; e holds the code. Move e into w[0..7].
    for k in 0..7 {
        circ.swap(w[k], &e[k]);
    }
    for q in e.into_iter().rev() {
        circ.zero_and_free(q);
    }
    circ.pop_section(&prev);
}

/// Inverse of `compress_3sym_base5_refs`: 7-bit code in `w[0..7]` (and
/// `w[7]`,`w[8]` = |0>) -> 3 symbols in `w[0..9]`.
pub fn compress_3sym_base5_reverse_refs(circ: &mut Circuit, w: &[&QReg; 9], dirty: &[&QReg]) {
    let prev = circ.push_section("b5pack");
    let e: Vec<QReg> = (0..7).map(|_| circ.alloc_qreg("b5_acc")).collect();
    for k in 0..7 {
        circ.swap(w[k], &e[k]);
    }
    for i in (0..3).rev() {
        unmerge_digit5(circ, w, &e, i, dirty);
    }
    for q in e.into_iter().rev() {
        circ.zero_and_free(q);
    }
    circ.pop_section(&prev);
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::point_add::trailmix_port::arith::schrottenloher::gcd_jump::{jump_bs_digit_to_sym, jump_bs_pack3};

    /// Compress 3 random valid symbols -> 7-bit code; check the code equals
    /// pack3 and w[7..9]=0; then decompress and check the symbols return.
    #[test]
    fn base5_compress_roundtrip() {
        use rand::{thread_rng, Rng};
        let mut rng = thread_rng();
        let mut circ = Circuit::new();
        circ.set_max_qubit_peak(64);
        // 9 window bits + 13 dirty (borrowed, loaded with junk to prove restore).
        let w: Vec<QReg> = (0..9).map(|i| circ.alloc_qreg(&format!("w{i}"))).collect();
        let dirty: Vec<QReg> = (0..13).map(|i| circ.alloc_qreg(&format!("d{i}"))).collect();

        let mut digits: Vec<[u8; 3]> = Vec::new();
        for shot in 0..64 {
            let d = [
                rng.gen_range(0..5u8),
                rng.gen_range(0..5u8),
                rng.gen_range(0..5u8),
            ];
            digits.push(d);
            for (i, &di) in d.iter().enumerate() {
                let (sub, sw, s2) = jump_bs_digit_to_sym(di);
                if sub == 1 {
                    circ.sim_load_reg_bytes_shot(std::slice::from_ref(&w[3 * i]), &[1], shot);
                }
                if sw == 1 {
                    circ.sim_load_reg_bytes_shot(std::slice::from_ref(&w[3 * i + 1]), &[1], shot);
                }
                if s2 == 1 {
                    circ.sim_load_reg_bytes_shot(std::slice::from_ref(&w[3 * i + 2]), &[1], shot);
                }
            }
            // dirty borrowed: random junk, must be restored.
            for d in &dirty {
                if rng.gen::<bool>() {
                    circ.sim_load_reg_bytes_shot(std::slice::from_ref(d), &[1], shot);
                }
            }
        }

        let wref: [&QReg; 9] = std::array::from_fn(|k| &w[k]);
        let dref: Vec<&QReg> = dirty.iter().collect();
        let ccx0 = circ.ccx_emitted;
        let ccz0 = circ.ccz_emitted;
        compress_3sym_base5_refs(&mut circ, &wref, &dref);
        // emitted Toffoli for one window (NOT per-shot -- ccx_emitted is a total).
        let comp_tof = circ.ccx_emitted - ccx0 + circ.ccz_emitted - ccz0;
        eprintln!("  b5 compress emitted-tof/window = {comp_tof}");
        // check the code on each shot via a contract (mid-circuit read).
        circ.contract_check("b5_code", |view, shot| {
            let mut code = 0u8;
            for k in 0..7 {
                if view.read_bit_shot(&w[k], shot) {
                    code |= 1 << k;
                }
            }
            let want = jump_bs_pack3(digits[shot]);
            if code != want {
                return Err(format!("code {code} != pack3 {want}"));
            }
            if view.read_bit_shot(&w[7], shot) || view.read_bit_shot(&w[8], shot) {
                return Err("w[7]/w[8] not cleared".into());
            }
            Ok(())
        });
        compress_3sym_base5_reverse_refs(&mut circ, &wref, &dref);

        circ.assert_phase_clean();
        let mut outs: Vec<QReg> = Vec::new();
        outs.extend(w);
        outs.extend(dirty);
        let (sim, detached) = circ.destroy_sim(outs);
        for shot in 0..64 {
            for (i, &di) in digits[shot].iter().enumerate() {
                let (sub, sw, s2) = jump_bs_digit_to_sym(di);
                assert_eq!(
                    sim.read_bit_shot(&detached[3 * i], shot),
                    sub,
                    "shot {shot} sym{i} sub"
                );
                assert_eq!(
                    sim.read_bit_shot(&detached[3 * i + 1], shot),
                    sw,
                    "shot {shot} sym{i} swap"
                );
                assert_eq!(
                    sim.read_bit_shot(&detached[3 * i + 2], shot),
                    s2,
                    "shot {shot} sym{i} s2"
                );
            }
        }
    }

    /// QROM pack/unpack: roundtrip correctness + emitted-Toffoli/window vs the
    /// radix packer (both are structural, so one window = the per-op cost).
    #[test]
    fn qrom_compress_roundtrip() {
        use crate::point_add::trailmix_port::arith::schrottenloher::gcd_jump::jump_bs_pack3;
        use rand::{thread_rng, Rng};
        let mut rng = thread_rng();
        let mut circ = Circuit::new();
        circ.set_max_qubit_peak(64);
        let w: Vec<QReg> = (0..9).map(|i| circ.alloc_qreg(&format!("w{i}"))).collect();
        let dirty: Vec<QReg> = (0..13).map(|i| circ.alloc_qreg(&format!("d{i}"))).collect();

        let mut digits: Vec<[u8; 3]> = Vec::new();
        for shot in 0..64 {
            let d = [
                rng.gen_range(0..5u8),
                rng.gen_range(0..5u8),
                rng.gen_range(0..5u8),
            ];
            digits.push(d);
            for (i, &di) in d.iter().enumerate() {
                let (sub, sw, s2) = jump_bs_digit_to_sym(di);
                for (off, bit) in [(0, sub), (1, sw), (2, s2)] {
                    if bit == 1 {
                        circ.sim_load_reg_bytes_shot(
                            std::slice::from_ref(&w[3 * i + off]),
                            &[1],
                            shot,
                        );
                    }
                }
            }
            for d in &dirty {
                if rng.gen::<bool>() {
                    circ.sim_load_reg_bytes_shot(std::slice::from_ref(d), &[1], shot);
                }
            }
        }

        let wref: [&QReg; 9] = std::array::from_fn(|k| &w[k]);
        let dref: Vec<&QReg> = dirty.iter().collect();
        let c0 = circ.ccx_emitted + circ.ccz_emitted;
        compress_3sym_qrom_refs(&mut circ, &wref, &dref);
        let pack_tof = circ.ccx_emitted + circ.ccz_emitted - c0;
        eprintln!("  QROM pack emitted-tof/window = {pack_tof} (radix was 18 comparators+9 adds)");

        circ.contract_check("qrom_code", |view, shot| {
            let mut code = 0u8;
            for k in 0..7 {
                if view.read_bit_shot(&w[k], shot) {
                    code |= 1 << k;
                }
            }
            let want = jump_bs_pack3(digits[shot]);
            if code != want {
                return Err(format!("code {code} != pack3 {want}"));
            }
            if view.read_bit_shot(&w[7], shot) || view.read_bit_shot(&w[8], shot) {
                return Err("w[7]/w[8] not cleared".into());
            }
            Ok(())
        });
        compress_3sym_qrom_reverse_refs(&mut circ, &wref, &dref);
        circ.assert_phase_clean();

        let mut outs: Vec<QReg> = Vec::new();
        outs.extend(w);
        outs.extend(dirty);
        let (sim, detached) = circ.destroy_sim(outs);
        for shot in 0..64 {
            for (i, &di) in digits[shot].iter().enumerate() {
                let (sub, sw, s2) = jump_bs_digit_to_sym(di);
                assert_eq!(
                    sim.read_bit_shot(&detached[3 * i], shot),
                    sub,
                    "shot {shot} sym{i} sub"
                );
                assert_eq!(
                    sim.read_bit_shot(&detached[3 * i + 1], shot),
                    sw,
                    "shot {shot} sym{i} swap"
                );
                assert_eq!(
                    sim.read_bit_shot(&detached[3 * i + 2], shot),
                    s2,
                    "shot {shot} sym{i} s2"
                );
            }
        }
    }
}
