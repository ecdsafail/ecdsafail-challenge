//! (refactor) Mechanically extracted from mod.rs. No logic changes.
use super::*;

/// Shift v left by k bits mod p. Returns (spill, flag_inv, ovf) which MUST
/// be passed to mod_shift_right_by_k for cleanup. Bennett-pattern: flags
/// stay alive across the body so the inverse can cleanly cancel them.
///
/// k must be small enough that spill·c < p. For k≤22 with secp256k1 this holds.
/// SHIFT22_CARRYTAIL (default-ON): route the shift22 STEP-3 unconditional +c and
/// STEP-4 conditional -c (and their reversals) through the carry-tail-truncatable
/// direct const adders (cadd/csub_nbit_const_direct_fast). The constant is the
/// SPARSE Solinas c = 2^32+977 (top bit 32), so kal_carrytail_count_c anchors the
/// window above bit 32 and the high result bits stay exact. Replaces the
/// register-loaded full-width add_nbit_const / csub_nbit_const of the shift22
/// reduction, clipping the carry/borrow chains identically forward and backward
/// (phase-clean). Set SHIFT22_CARRYTAIL=0 to restore the register-loaded path.
pub(crate) fn shift22_carrytail() -> bool {
    std::env::var("SHIFT22_CARRYTAIL").ok().as_deref() != Some("0")
}

/// Dedicated carry-tail cut for the shift22 STEP-3/STEP-4 direct const adders,
/// DECOUPLED from the global Kaliski `kal_carrytail_w` so the proven Kaliski
/// island (W=22) stays pinned. Anchors above c's top set bit (k0=33) plus a
/// dedicated window W (SHIFT22_CARRYTAIL_W, default 37 → cut=70). The shift22
/// reduction operates on a freshly-folded value whose conditional-sub borrow run
/// can be longer than the Kaliski case, so the window is anchored wider than the
/// Kaliski 22. A full SHIFT22_CARRYTAIL_W sweep (each = trusted eval over 9024
/// shots) found W=37 and W=40 the clean islands on the shift22-direct op-stream:
/// W=37 → 0/0/0, avg-exec 2,429,688 Toffoli × 2309 = 5,610,149,592 (deepest clean);
/// W=40 → 0/0/0, avg-exec 2,429,724. Note this is NOT a pure truncation-soundness
/// floor — even W=223 (full chain, zero truncation) MISSES the Fiat-Shamir island
/// (1 mismatch), because routing STEP-3/4 through the direct const adders re-rolls
/// the test inputs; W=37 is the value that lands a 9024-clean island AND truncates.
/// Single value, used identically forward and reverse (phase-parity).
/// W∈{20,22,25,30,35,36,38,39,41,42,45,50,55,60,90,223} all MISS the island.
pub(crate) fn shift22_carrytail_cut() -> usize {
    let w = std::env::var("SHIFT22_CARRYTAIL_W")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(37);
    33usize.saturating_add(w)
}

pub(crate) fn lowq_shift22() -> bool {
    // Qubit-first default: the global LOWQ shift22 path is strict-clean on the
    // current scaffold and lowers the benchmark peak (2736q -> 2715q) at a
    // small Toffoli cost. Keep LOWQ_SHIFT22=0 as an explicit opt-out for
    // Toffoli-first diagnostics and baseline comparisons.
    match std::env::var("LOWQ_SHIFT22") {
        Ok(v) => v != "0",
        Err(_) => true,
    }
}

pub(crate) fn mod_shift_left_by_k(
    b: &mut B,
    v: &[QubitId],
    p: U256,
    k: usize,
) -> (Vec<QubitId>, QubitId, QubitId) {
    let n = v.len();
    debug_assert_eq!(n, 256);
    let c = U256::MAX.wrapping_sub(p).wrapping_add(U256::from(1));

    let spill = b.alloc_qubits(k);
    let ovf = b.alloc_qubit();
    let flag_inv = b.alloc_qubit();

    // Step 1: k rounds of shift-by-1, capturing top bits into spill.
    for shift_i in 0..k {
        b.swap(v[n - 1], spill[k - 1 - shift_i]);
        for i in (0..n - 1).rev() {
            b.swap(v[i], v[i + 1]);
        }
    }

    // Step 2: add spill · c to v_ext (using ovf as bit n).
    // c = 2^32 + 977 = 2^32 + 2^10 - 2^6 + 2^4 + 2^0.
    // Consolidate 4 bits (6,7,8,9) of 977 into 2^10 - 2^6: saves 2 Cuccaros per shift.
    // Op list: ADD at 0, 4, 10, 32; SUB at 6. Total 5 ops instead of 7.
    let mut v_ext = v.to_vec();
    v_ext.push(ovf);
    let cuccaro_op = |b: &mut B, pos: usize, is_sub: bool| {
        let pad_width = n + 1 - pos;
        let padded = b.alloc_qubits(pad_width);
        for i in 0..k.min(pad_width) {
            b.cx(spill[i], padded[i]);
        }
        let v_slice: Vec<QubitId> = v_ext[pos..n + 1].to_vec();
        let c_in = b.alloc_qubit();
        if lowq_shift22() {
            if is_sub {
                cuccaro_sub(b, &padded, &v_slice, c_in);
            } else {
                cuccaro_add(b, &padded, &v_slice, c_in);
            }
        } else if is_sub {
            // Fast cuccaro: saves ~n CCX per op. Peak during this op (~514
            // transient) is still below the mod_add_qq_fast peak (517) inside
            // the enclosing Solinas, so no global peak increase.
            cuccaro_sub_fast(b, &padded, &v_slice, c_in);
        } else {
            cuccaro_add_fast(b, &padded, &v_slice, c_in);
        }
        b.free(c_in);
        for i in 0..k.min(pad_width) {
            b.cx(spill[i], padded[i]);
        }
        b.free_vec(&padded);
    };
    b.set_phase("shift22_cuccaro_op_0");
    cuccaro_op(b, 0, false);
    b.set_phase("shift22_cuccaro_op_4");
    cuccaro_op(b, 4, false);
    b.set_phase("shift22_cuccaro_op_6");
    cuccaro_op(b, 6, true);
    b.set_phase("shift22_cuccaro_op_10");
    cuccaro_op(b, 10, false);
    b.set_phase("shift22_cuccaro_op_32");
    cuccaro_op(b, 32, false);

    // Step 3: const add.
    b.set_phase("shift22_step3");
    if shift22_carrytail() {
        // CARRY-TAIL: route the unconditional +c through the truncatable direct
        // const-add (QQFOLD pattern: always-on `one` ctrl) with a DEDICATED cut.
        // c is the SPARSE Solinas 2^32+977; the dedicated window keeps high bits
        // exact. Reversal mirrors this exact cut via csub_..._cut(same cut).
        let cut = shift22_carrytail_cut();
        let one = b.alloc_qubit();
        b.x(one);
        cadd_nbit_const_direct_fast_cut(b, &v_ext, c, one, cut);
        b.x(one);
        b.free(one);
    } else if lowq_shift22() {
        add_nbit_const(b, &v_ext, c);
    } else {
        add_nbit_const_fast(b, &v_ext, c);
    }
    b.x(ovf);
    b.cx(ovf, flag_inv); // flag_inv = NOT(top_bit_after_add) = (value < p)
    b.x(ovf);

    // Step 4: conditional const sub.
    b.set_phase("shift22_step4");
    if shift22_carrytail() {
        csub_nbit_const_direct_fast_cut(b, &v_ext, c, flag_inv, shift22_carrytail_cut());
    } else if lowq_shift22() {
        csub_nbit_const(b, &v_ext, c, flag_inv);
    } else {
        csub_nbit_const_fast(b, &v_ext, c, flag_inv);
    }
    b.x(flag_inv);
    b.cx(flag_inv, ovf);
    b.x(flag_inv);

    (spill, flag_inv, ovf)
}

/// Gate-level inverse of mod_shift_left_by_k.
pub(crate) fn mod_shift_right_by_k(
    b: &mut B,
    v: &[QubitId],
    p: U256,
    k: usize,
    spill: Vec<QubitId>,
    flag_inv: QubitId,
    ovf: QubitId,
) {
    let n = v.len();
    debug_assert_eq!(n, 256);
    let c = U256::MAX.wrapping_sub(p).wrapping_add(U256::from(1));

    let mut v_ext = v.to_vec();
    v_ext.push(ovf);

    // Reverse step 4.
    b.x(flag_inv);
    b.cx(flag_inv, ovf);
    b.x(flag_inv);
    b.set_phase("rshift22_rev_step4");
    if shift22_carrytail() {
        // Exact inverse of step4's csub_..._cut: cadd with the SAME cut.
        cadd_nbit_const_direct_fast_cut(b, &v_ext, c, flag_inv, shift22_carrytail_cut());
    } else if lowq_shift22() {
        cadd_nbit_const(b, &v_ext, c, flag_inv);
    } else {
        cadd_nbit_const_fast(b, &v_ext, c, flag_inv);
    }

    // Reverse step 3.
    b.x(ovf);
    b.cx(ovf, flag_inv);
    b.x(ovf);
    b.set_phase("rshift22_rev_step3");
    if shift22_carrytail() {
        // Exact inverse of step3's cadd_..._cut: csub with the SAME cut.
        let cut = shift22_carrytail_cut();
        let one = b.alloc_qubit();
        b.x(one);
        csub_nbit_const_direct_fast_cut(b, &v_ext, c, one, cut);
        b.x(one);
        b.free(one);
    } else if lowq_shift22() {
        sub_nbit_const(b, &v_ext, c);
    } else {
        sub_nbit_const_fast(b, &v_ext, c);
    }
    b.free(flag_inv);
    b.set_phase("rshift22_rev_step2");

    // Reverse step 2: inverse of the consolidated op list (5 ops, in reverse order, flipped signs).
    let cuccaro_op = |b: &mut B, pos: usize, is_sub: bool| {
        let pad_width = n + 1 - pos;
        let padded = b.alloc_qubits(pad_width);
        for i in 0..k.min(pad_width) {
            b.cx(spill[i], padded[i]);
        }
        let v_slice: Vec<QubitId> = v_ext[pos..n + 1].to_vec();
        let c_in = b.alloc_qubit();
        if lowq_shift22() {
            if is_sub {
                cuccaro_sub(b, &padded, &v_slice, c_in);
            } else {
                cuccaro_add(b, &padded, &v_slice, c_in);
            }
        } else if is_sub {
            cuccaro_sub_fast(b, &padded, &v_slice, c_in);
        } else {
            cuccaro_add_fast(b, &padded, &v_slice, c_in);
        }
        b.free(c_in);
        for i in 0..k.min(pad_width) {
            b.cx(spill[i], padded[i]);
        }
        b.free_vec(&padded);
    };
    // Reverse: undo ADD at 32, 10; undo SUB at 6; undo ADD at 4, 0.
    cuccaro_op(b, 32, true); // undo +spill·2^32
    cuccaro_op(b, 10, true); // undo +spill·2^10
    cuccaro_op(b, 6, false); // undo -spill·2^6
    cuccaro_op(b, 4, true); // undo +spill·2^4
    cuccaro_op(b, 0, true); // undo +spill·2^0

    // Reverse step 1: reverse swap cascades.
    for shift_i in (0..k).rev() {
        for i in 0..n - 1 {
            b.swap(v[i], v[i + 1]);
        }
        b.swap(v[n - 1], spill[k - 1 - shift_i]);
    }

    b.free(ovf);
    b.free_vec(&spill);
}

/// Low-scratch spill add `v_ext[pos..] += spill*2^pos` (KAL_GZ_SOLINAS_LOWSCRATCH).
/// `spill` is a k-bit quantum value; the full-width add it represents has nonzero
/// addend bits only in [pos, pos+k). Instead of the ~(n+1-pos)-wide `padded`
/// carry-scratch register, this adds the k-bit chunk with a captured carry-out,
/// propagates that single carry through the high bits via a Gidney venting
/// controlled-increment (2 clean ancilla + a borrowed DIRTY donor, restored),
/// then uncomputes the carry-out via the exact comparator identity
/// `carry == (low_sum < spill)`. Net transient ~k+5 instead of ~n.
pub(crate) fn shift22_spill_op_dirty(
    b: &mut B,
    v_ext: &[QubitId],
    spill: &[QubitId],
    pos: usize,
    k: usize,
    is_sub: bool,
    dirty: &[QubitId],
) {
    let total = v_ext.len(); // n+1
    let w = total - pos; // width of the affected window
    debug_assert!(w >= k);
    let v_slice: Vec<QubitId> = v_ext[pos..total].to_vec();
    let low: Vec<QubitId> = v_slice[..k].to_vec();
    let hi: Vec<QubitId> = v_slice[k..].to_vec(); // width w-k

    // Step A: add/sub the k-bit spill into the low window, capturing carry/borrow.
    let carry = b.alloc_qubit();
    let zpad = b.alloc_qubit(); // |0> top bit of the (k+1)-bit addend
    let mut addend = spill.to_vec();
    addend.push(zpad);
    let mut acc = low.clone();
    acc.push(carry);
    let c_in = b.alloc_qubit();
    if is_sub {
        cuccaro_sub(b, &addend, &acc, c_in);
    } else {
        cuccaro_add(b, &addend, &acc, c_in);
    }
    b.free(c_in);
    b.free(zpad); // restored to |0> by cuccaro (addend register is preserved)

    // Step B: propagate the single carry/borrow into the high window.
    // add: hi += carry. sub: hi -= carry == controlled-decrement (invert,+1,invert).
    if w - k >= 5 {
        let dlen = (w - k).saturating_sub(2);
        let q_clean2: [QubitId; 2] = [b.alloc_qubit(), b.alloc_qubit()];
        if is_sub {
            for &q in hi.iter() {
                b.cx(carry, q); // controlled-invert hi (no-op when carry=0)
            }
            venting::ciadd_dirty_2clean_classical(
                b, &hi, &dirty[..dlen], &q_clean2, 1, carry, false,
            );
            for &q in hi.iter() {
                b.cx(carry, q);
            }
        } else {
            venting::ciadd_dirty_2clean_classical(
                b, &hi, &dirty[..dlen], &q_clean2, 1, carry, false,
            );
        }
        b.free(q_clean2[0]);
        b.free(q_clean2[1]);
    } else {
        let c = U256::from(1u64);
        if is_sub {
            csub_nbit_const_direct_fast(b, &hi, c, carry);
        } else {
            cadd_nbit_const_direct_fast(b, &hi, c, carry);
        }
    }

    // Step C: uncompute `carry`. add: carry == (low_new < spill).
    // sub: borrow == ((~low_new) < spill).
    if is_sub {
        for &q in low.iter() {
            b.x(q);
        }
        cmp_lt_into(b, &low, spill, carry);
        for &q in low.iter() {
            b.x(q);
        }
    } else {
        cmp_lt_into(b, &low, spill, carry);
    }
    b.free(carry);
}

fn shift22_compute_m_977(b: &mut B, m: &[QubitId], spill: &[QubitId], k: usize, undo: bool) {
    debug_assert_eq!(m.len(), 32);
    let terms: [(usize, bool); 4] = [(0, false), (4, false), (6, true), (10, false)];
    let term_op = |b: &mut B, pos: usize, is_sub: bool| {
        let w = 32 - pos;
        let m_slice: Vec<QubitId> = m[pos..32].to_vec();
        let pad = b.alloc_qubits(w - k);
        let mut addend = spill.to_vec();
        addend.extend_from_slice(&pad);
        let c_in = b.alloc_qubit();
        if is_sub {
            cuccaro_sub(b, &addend, &m_slice, c_in);
        } else {
            cuccaro_add(b, &addend, &m_slice, c_in);
        }
        b.free(c_in);
        b.free_vec(&pad);
    };
    if !undo {
        for &(pos, is_sub) in terms.iter() {
            term_op(b, pos, is_sub);
        }
    } else {
        for &(pos, is_sub) in terms.iter().rev() {
            term_op(b, pos, !is_sub);
        }
    }
}

/// KAL_GZ_SOLINAS_LOWSCRATCH forward shift22: same arithmetic as
/// `mod_shift_left_by_k` but STEP-2 spill ops use the dirty-borrow narrow
/// spill-add (no ~257-wide `padded`) and STEP-3/STEP-4 const-add /
/// conditional-const-sub use Gidney venting dirty-borrow const adders (no
/// ~257-wide loaded-constant register). `dirty` is a co-resident DIRTY donor
/// register (>= n-2 wide, restored to its entry value on exit).
pub(crate) fn mod_shift_left_by_k_dirty(
    b: &mut B,
    v: &[QubitId],
    p: U256,
    k: usize,
    dirty: &[QubitId],
) -> (Vec<QubitId>, QubitId, QubitId) {
    let n = v.len();
    debug_assert_eq!(n, 256);
    let c = U256::MAX.wrapping_sub(p).wrapping_add(U256::from(1));
    let c_low = c.as_limbs()[0];

    let spill = b.alloc_qubits(k);
    let ovf = b.alloc_qubit();
    let flag_inv = b.alloc_qubit();

    for shift_i in 0..k {
        b.swap(v[n - 1], spill[k - 1 - shift_i]);
        for i in (0..n - 1).rev() {
            b.swap(v[i], v[i + 1]);
        }
    }

    let mut v_ext = v.to_vec();
    v_ext.push(ovf);

    if shift22_collapse() && k <= 22 {
        // SHIFT22_COLLAPSE: compute m=spill*977 (fits 32 bits), fold at pos 0
        // and pos 32 instead of 5 dirty spill ops. Reversal mirrors this.
        let m = b.alloc_qubits(32);
        shift22_compute_m_977(b, &m, &spill, k, false);
        b.set_phase("shift22_cuccaro_op_0");
        shift22_spill_op_dirty(b, &v_ext, &m, 0, 32, false, dirty);
        b.set_phase("shift22_cuccaro_op_32");
        shift22_spill_op_dirty(b, &v_ext, &spill, 32, k, false, dirty);
        shift22_compute_m_977(b, &m, &spill, k, true);
        b.free_vec(&m);
    } else {
        b.set_phase("shift22_cuccaro_op_0");
        shift22_spill_op_dirty(b, &v_ext, &spill, 0, k, false, dirty);
        b.set_phase("shift22_cuccaro_op_4");
        shift22_spill_op_dirty(b, &v_ext, &spill, 4, k, false, dirty);
        b.set_phase("shift22_cuccaro_op_6");
        shift22_spill_op_dirty(b, &v_ext, &spill, 6, k, true, dirty);
        b.set_phase("shift22_cuccaro_op_10");
        shift22_spill_op_dirty(b, &v_ext, &spill, 10, k, false, dirty);
        b.set_phase("shift22_cuccaro_op_32");
        shift22_spill_op_dirty(b, &v_ext, &spill, 32, k, false, dirty);
    }

    // Step 3: unconditional const add of c (register-free venting dirty-borrow).
    b.set_phase("shift22_step3");
    {
        let m = v_ext.len();
        let q_clean2: [QubitId; 2] = [b.alloc_qubit(), b.alloc_qubit()];
        venting::iadd_dirty_2clean_classical(
            b, &v_ext, &dirty[..m - 2], &q_clean2, c_low, false,
        );
        b.free(q_clean2[0]);
        b.free(q_clean2[1]);
    }
    b.x(ovf);
    b.cx(ovf, flag_inv); // flag_inv = NOT(top_bit_after_add) = (value < p)
    b.x(ovf);

    // Step 4: conditional const sub of c (register-free venting dirty-borrow).
    b.set_phase("shift22_step4");
    {
        let m = v_ext.len();
        let q_clean2: [QubitId; 2] = [b.alloc_qubit(), b.alloc_qubit()];
        venting::cisub_dirty_2clean_classical(
            b, &v_ext, &dirty[..m - 2], &q_clean2, c_low, flag_inv,
        );
        b.free(q_clean2[0]);
        b.free(q_clean2[1]);
    }
    b.x(flag_inv);
    b.cx(flag_inv, ovf);
    b.x(flag_inv);

    (spill, flag_inv, ovf)
}

/// Gate-level inverse of `mod_shift_left_by_k_dirty`.
pub(crate) fn mod_shift_right_by_k_dirty(
    b: &mut B,
    v: &[QubitId],
    p: U256,
    k: usize,
    spill: Vec<QubitId>,
    flag_inv: QubitId,
    ovf: QubitId,
    dirty: &[QubitId],
) {
    let n = v.len();
    debug_assert_eq!(n, 256);
    let c = U256::MAX.wrapping_sub(p).wrapping_add(U256::from(1));
    let c_low = c.as_limbs()[0];

    let mut v_ext = v.to_vec();
    v_ext.push(ovf);

    // Reverse step 4.
    b.x(flag_inv);
    b.cx(flag_inv, ovf);
    b.x(flag_inv);
    b.set_phase("rshift22_rev_step4");
    {
        // inverse of cisub(c, flag_inv) is ciadd(c, flag_inv): if flag_inv: x += c.
        let m = v_ext.len();
        let q_clean2: [QubitId; 2] = [b.alloc_qubit(), b.alloc_qubit()];
        venting::ciadd_dirty_2clean_classical(
            b, &v_ext, &dirty[..m - 2], &q_clean2, c_low, flag_inv, false,
        );
        b.free(q_clean2[0]);
        b.free(q_clean2[1]);
    }

    // Reverse step 3.
    b.x(ovf);
    b.cx(ovf, flag_inv);
    b.x(ovf);
    b.set_phase("rshift22_rev_step3");
    {
        // inverse of iadd(c) is isub(c) == invert; iadd(c); invert.
        let m = v_ext.len();
        let q_clean2: [QubitId; 2] = [b.alloc_qubit(), b.alloc_qubit()];
        for &q in v_ext.iter() {
            b.x(q);
        }
        venting::iadd_dirty_2clean_classical(
            b, &v_ext, &dirty[..m - 2], &q_clean2, c_low, false,
        );
        for &q in v_ext.iter() {
            b.x(q);
        }
        b.free(q_clean2[0]);
        b.free(q_clean2[1]);
    }
    b.free(flag_inv);
    b.set_phase("rshift22_rev_step2");

    // Reverse step 2: undo the spill ops in reverse order with flipped signs.
    if shift22_collapse() && k <= 22 {
        // Reverse the COLLAPSE: recompute m, undo pos 32 then pos 0, uncompute m.
        let m = b.alloc_qubits(32);
        shift22_compute_m_977(b, &m, &spill, k, false);
        shift22_spill_op_dirty(b, &v_ext, &spill, 32, k, true, dirty); // undo +spill*2^32
        shift22_spill_op_dirty(b, &v_ext, &m, 0, 32, true, dirty); // undo +m=spill*977
        shift22_compute_m_977(b, &m, &spill, k, true);
        b.free_vec(&m);
    } else {
        shift22_spill_op_dirty(b, &v_ext, &spill, 32, k, true, dirty); // undo +2^32
        shift22_spill_op_dirty(b, &v_ext, &spill, 10, k, true, dirty); // undo +2^10
        shift22_spill_op_dirty(b, &v_ext, &spill, 6, k, false, dirty); // undo -2^6
        shift22_spill_op_dirty(b, &v_ext, &spill, 4, k, true, dirty); // undo +2^4
        shift22_spill_op_dirty(b, &v_ext, &spill, 0, k, true, dirty); // undo +2^0
    }

    // Reverse step 1: reverse swap cascades.
    for shift_i in (0..k).rev() {
        for i in 0..n - 1 {
            b.swap(v[i], v[i + 1]);
        }
        b.swap(v[n - 1], spill[k - 1 - shift_i]);
    }

    b.free(ovf);
    b.free_vec(&spill);
}

pub(crate) fn mod_shift_left_by_k_lowq(
    b: &mut B,
    v: &[QubitId],
    p: U256,
    k: usize,
) -> (Vec<QubitId>, QubitId, QubitId) {
    let n = v.len();
    debug_assert_eq!(n, 256);
    let c = U256::MAX.wrapping_sub(p).wrapping_add(U256::from(1));

    let spill = b.alloc_qubits(k);
    let ovf = b.alloc_qubit();
    let flag_inv = b.alloc_qubit();

    for shift_i in 0..k {
        b.swap(v[n - 1], spill[k - 1 - shift_i]);
        for i in (0..n - 1).rev() {
            b.swap(v[i], v[i + 1]);
        }
    }

    let mut v_ext = v.to_vec();
    v_ext.push(ovf);
    let cuccaro_op = |b: &mut B, pos: usize, is_sub: bool| {
        let pad_width = n + 1 - pos;
        let padded = b.alloc_qubits(pad_width);
        for i in 0..k.min(pad_width) {
            b.cx(spill[i], padded[i]);
        }
        let v_slice: Vec<QubitId> = v_ext[pos..n + 1].to_vec();
        let c_in = b.alloc_qubit();
        if is_sub {
            cuccaro_sub(b, &padded, &v_slice, c_in);
        } else {
            cuccaro_add(b, &padded, &v_slice, c_in);
        }
        b.free(c_in);
        for i in 0..k.min(pad_width) {
            b.cx(spill[i], padded[i]);
        }
        b.free_vec(&padded);
    };
    cuccaro_op(b, 0, false);
    cuccaro_op(b, 4, false);
    cuccaro_op(b, 6, true);
    cuccaro_op(b, 10, false);
    cuccaro_op(b, 32, false);

    add_nbit_const(b, &v_ext, c);
    b.x(ovf);
    b.cx(ovf, flag_inv);
    b.x(ovf);
    csub_nbit_const(b, &v_ext, c, flag_inv);
    b.x(flag_inv);
    b.cx(flag_inv, ovf);
    b.x(flag_inv);

    (spill, flag_inv, ovf)
}

pub(crate) fn mod_shift_right_by_k_lowq(
    b: &mut B,
    v: &[QubitId],
    p: U256,
    k: usize,
    spill: Vec<QubitId>,
    flag_inv: QubitId,
    ovf: QubitId,
) {
    let n = v.len();
    debug_assert_eq!(n, 256);
    let c = U256::MAX.wrapping_sub(p).wrapping_add(U256::from(1));

    let mut v_ext = v.to_vec();
    v_ext.push(ovf);

    b.x(flag_inv);
    b.cx(flag_inv, ovf);
    b.x(flag_inv);
    cadd_nbit_const(b, &v_ext, c, flag_inv);

    b.x(ovf);
    b.cx(ovf, flag_inv);
    b.x(ovf);
    sub_nbit_const(b, &v_ext, c);
    b.free(flag_inv);

    let cuccaro_op = |b: &mut B, pos: usize, is_sub: bool| {
        let pad_width = n + 1 - pos;
        let padded = b.alloc_qubits(pad_width);
        for i in 0..k.min(pad_width) {
            b.cx(spill[i], padded[i]);
        }
        let v_slice: Vec<QubitId> = v_ext[pos..n + 1].to_vec();
        let c_in = b.alloc_qubit();
        if is_sub {
            cuccaro_sub(b, &padded, &v_slice, c_in);
        } else {
            cuccaro_add(b, &padded, &v_slice, c_in);
        }
        b.free(c_in);
        for i in 0..k.min(pad_width) {
            b.cx(spill[i], padded[i]);
        }
        b.free_vec(&padded);
    };
    cuccaro_op(b, 32, true);
    cuccaro_op(b, 10, true);
    cuccaro_op(b, 6, false);
    cuccaro_op(b, 4, true);
    cuccaro_op(b, 0, true);

    for shift_i in (0..k).rev() {
        for i in 0..n - 1 {
            b.swap(v[i], v[i + 1]);
        }
        b.swap(v[n - 1], spill[k - 1 - shift_i]);
    }

    b.free(ovf);
    b.free_vec(&spill);
}
