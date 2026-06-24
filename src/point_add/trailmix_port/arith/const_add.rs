//! Add/subtract of a compile-time constant into a quantum register, with
//! clustered / windowed / sparse / runs-forced encodings that exploit the
//! constant's bit structure. Extracted from `poc_arith`.

use crate::point_add::trailmix_port::circuit::{Circuit, QReg};

/// Step descriptor for the run-based carry automaton. We hold a Vec<QReg>
/// of all the "carry" `QRegs` alive concurrently (each step has its own); the
/// step records the indexes into that Vec rather than the `QReg` itself
/// (which is non-Copy / non-Clone).
#[derive(Clone, Copy)]
enum RunCarryStep {
    FirstOne {
        lo: usize,
        hi: usize,
        carry_idx: usize,
    },
    Zero {
        lo: usize,
        hi: usize,
        prev_carry_idx: usize,
        carry_idx: usize,
    },
    One {
        lo: usize,
        hi: usize,
        prev_carry_idx: usize,
        carry_idx: usize,
    },
}

#[must_use]
pub fn get_const_bit(bytes: &[u8], i: usize) -> bool {
    let byte_idx = i / 8;
    let bit_idx = i % 8;
    byte_idx < bytes.len() && (bytes[byte_idx] >> bit_idx) & 1 == 1
}

/// Controlled add-constant: if ctrl=1, a += val.
///
/// Dispatches between two backends:
/// - **Sparse** (popcount x 5n threshold): iterate over the set bits
///   of `val` and, for each, call `cinc_khattar_gidney`
///   from that position. Cost ~= popcount x cinc(n-pos). Best when
///   popcount is small (e.g. rfold R = 2^32+977 at popcount 7).
/// - **Dense** (Theorem 5 via `controlled_classical_quantum_add`):
///   Theta(n log^2 n) single pass. Best when popcount ~= n/2.
///
/// Threshold is popcount <= log2(n) (a conservative value; Theorem 5
/// beats sparse at roughly popcount > log2(n) * const).
pub fn controlled_add_const(circ: &mut Circuit, ctrl: &QReg, a: &[QReg], val: &[u8]) {
    let n = a.len();
    if n == 0 {
        return;
    }

    // Find low/high set bits and popcount of val within first n bits.
    let mut lo_bit = usize::MAX;
    let mut hi_bit = 0usize;
    let mut pop = 0usize;
    for i in 0..n {
        if get_const_bit(val, i) {
            if lo_bit == usize::MAX {
                lo_bit = i;
            }
            hi_bit = i;
            pop += 1;
        }
    }
    if pop == 0 {
        return;
    }
    if pop == 1 {
        // Single bit: just one cinc from that position.
        crate::point_add::trailmix_port::arith::khattar_gidney::cinc_khattar_gidney(circ, &a[lo_bit..], ctrl);
        return;
    }
    let runs = one_runs(val, n);

    // Windowed path: set bits all lie in [lo_bit, hi_bit]. If the
    // window width is small relative to n, adding via
    //   (1) controlled classq_add on a[lo..=hi] with c = val restricted
    //       to window bits, and
    //   (2) a single controlled cinc_khattar_gidney on a[hi+1..] for the carry-out,
    // is much cheaper than a cinc_khattar_gidney per set bit. Concretely for
    // val = R = 2^32 + 977 (popcount=7, window=[0,32]): classq_add(33)
    // + witness + cinc(223) ~= 15K ops vs 79K for 7 separate cincs.
    //
    // Heuristic: clustered constants benefit from the run-based
    // decomposition below; narrow windows benefit from the generic
    // classq-add path; very sparse wide constants still prefer the
    // per-bit suffix increments.
    let window = hi_bit - lo_bit + 1;
    let segment_count = runs.len().saturating_mul(2).saturating_sub(1);
    // Empirically (profile_add_const_f for f=2^32+977 into 63-bit reg):
    //   runs_forced  359 Tof   <- always cheapest for sparse multi-bit constants
    //   sparse      1109 Tof
    //   windowed    1458 Tof
    //   Theorem 5   2116 Tof
    // The old heuristic gated `runs_forced` on `window <= n/3` which is FALSE
    // for f (window=33 vs lsbs=63/3=21), routing to Theorem 5. Empirically,
    // runs_forced wins as long as #runs is modest — extend the threshold so
    // sparse-but-wide constants like f use it.
    if runs.len() <= 8 && segment_count <= 15 {
        controlled_add_const_runs_forced(circ, ctrl, a, val);
    } else if window <= n / 3 {
        controlled_add_const_windowed(circ, ctrl, a, val, lo_bit, hi_bit);
    } else if pop <= (n.trailing_zeros() as usize).max(1) + 4 {
        controlled_add_const_sparse(circ, ctrl, a, val);
    } else {
        crate::point_add::trailmix_port::arith::khattar_gidney::controlled_classical_quantum_add(circ, ctrl, a, val);
    }
}

fn one_runs(val: &[u8], n: usize) -> Vec<(usize, usize)> {
    let mut runs = Vec::new();
    let mut i = 0usize;
    while i < n {
        if !get_const_bit(val, i) {
            i += 1;
            continue;
        }
        let lo = i;
        while i + 1 < n && get_const_bit(val, i + 1) {
            i += 1;
        }
        runs.push((lo, i));
        i += 1;
    }
    runs
}

fn xor_all_ones(circ: &mut Circuit, block: &[QReg], target: &QReg) {
    let block_refs: Vec<&QReg> = block.iter().collect();
    crate::point_add::trailmix_port::arith::mcx::mcx_clean_k(circ, &block_refs, target);
}

/// Like `xor_all_ones` but frees `target` at its last gate-touch inside
/// `mcx_clean_k_uncompute_consume`, before the uncompute step allocs ancillae.
fn xor_all_ones_consume_free(circ: &mut Circuit, block: &[QReg], target: QReg) {
    let block_refs: Vec<&QReg> = block.iter().collect();
    crate::point_add::trailmix_port::arith::mcx::mcx_clean_k_uncompute_consume(circ, &block_refs, target);
}

fn apply_conditional_decrement(circ: &mut Circuit, block: &[QReg], ctrl: &QReg) {
    // Decrement by ctrl via X-sandwich + cinc_khattar_gidney (O(n log* n)
    // CCX/CX). The leading/trailing X-loops on `block` produce adjacent
    // X-X pairs only on bits the inner cinc never touches; auto-elide
    // cancels those at push time, leaving the optimal sequence.
    for q in block {
        circ.x(q);
    }
    crate::point_add::trailmix_port::arith::khattar_gidney::cinc_khattar_gidney(circ, block, ctrl);
    for q in block {
        circ.x(q);
    }
}

fn carry_after_first_one_run(circ: &mut Circuit, ctrl: &QReg, block: &[QReg], carry: &QReg) {
    let all_ones = circ.alloc_qreg("run_all_ones");
    xor_all_ones(circ, block, &all_ones);
    circ.cx(ctrl, carry);
    circ.ccx(ctrl, &all_ones, carry);
    // Free all_ones at its last gate-touch before any subsequent allocs.
    xor_all_ones_consume_free(circ, block, all_ones);
}

/// Uncompute variant: frees carry before `xor_all_ones_consume_free` allocs ancillae.
fn carry_after_first_one_run_uncompute_free_carry(
    circ: &mut Circuit,
    ctrl: &QReg,
    block: &[QReg],
    carry: QReg,
) {
    let all_ones = circ.alloc_qreg("run_all_ones");
    xor_all_ones(circ, block, &all_ones);
    circ.cx(ctrl, &carry);
    circ.ccx(ctrl, &all_ones, &carry);
    drop(carry);
    xor_all_ones_consume_free(circ, block, all_ones);
}

fn carry_after_zero_run(circ: &mut Circuit, prev_carry: &QReg, block: &[QReg], carry: &QReg) {
    // Inlined to avoid cancelling X(block[i]) at the boundary between
    // xor_all_zeros and xor_all_zeros_consume_free: the trailing ~block
    // of the first and the leading ~block of the second are adjacent
    // X-X pairs on each block bit (separated only by the ccx on
    // prev_carry/all_zeros/carry, which doesn't touch block).
    let all_zeros = circ.alloc_qreg("run_all_zeros");
    let block_refs: Vec<&QReg> = block.iter().collect();
    for q in block {
        circ.x(q);
    }
    crate::point_add::trailmix_port::arith::mcx::mcx_clean_k(circ, &block_refs, &all_zeros);
    circ.ccx(prev_carry, &all_zeros, carry);
    crate::point_add::trailmix_port::arith::mcx::mcx_clean_k_uncompute_consume(circ, &block_refs, all_zeros);
    for q in block {
        circ.x(q);
    }
}

/// Uncompute variant of `carry_after_zero_run`: frees `carry` at its last
/// gate-touch (the ccx), BEFORE `xor_all_zeros_consume_free` allocs ancillae
/// that would push `last_alloc_op_idx` past carry's last touch.
fn carry_after_zero_run_uncompute_free_carry(
    circ: &mut Circuit,
    prev_carry: &QReg,
    block: &[QReg],
    carry: QReg,
) {
    let all_zeros = circ.alloc_qreg("run_all_zeros");
    let block_refs: Vec<&QReg> = block.iter().collect();
    for q in block {
        circ.x(q);
    }
    crate::point_add::trailmix_port::arith::mcx::mcx_clean_k(circ, &block_refs, &all_zeros);
    circ.ccx(prev_carry, &all_zeros, &carry);
    drop(carry);
    crate::point_add::trailmix_port::arith::mcx::mcx_clean_k_uncompute_consume(circ, &block_refs, all_zeros);
    for q in block {
        circ.x(q);
    }
}

fn carry_after_one_run(
    circ: &mut Circuit,
    ctrl: &QReg,
    prev_carry: &QReg,
    block: &[QReg],
    carry: &QReg,
) {
    let dec_ctrl = circ.alloc_qreg("run_dec_ctrl");
    circ.cx(ctrl, &dec_ctrl);
    circ.cx(prev_carry, &dec_ctrl);

    let all_ones = circ.alloc_qreg("run_all_ones");
    xor_all_ones(circ, block, &all_ones);

    circ.cx(ctrl, carry);
    circ.ccx(&dec_ctrl, &all_ones, carry);

    // Free all_ones at its last gate-touch before the uncompute allocs ancillae.
    xor_all_ones_consume_free(circ, block, all_ones);

    circ.cx(prev_carry, &dec_ctrl);
    circ.cx(ctrl, &dec_ctrl);
    drop(dec_ctrl);
}

/// Uncompute variant: frees carry before `xor_all_ones_consume_free` allocs ancillae.
fn carry_after_one_run_uncompute_free_carry(
    circ: &mut Circuit,
    ctrl: &QReg,
    prev_carry: &QReg,
    block: &[QReg],
    carry: QReg,
) {
    let dec_ctrl = circ.alloc_qreg("run_dec_ctrl");
    circ.cx(ctrl, &dec_ctrl);
    circ.cx(prev_carry, &dec_ctrl);

    let all_ones = circ.alloc_qreg("run_all_ones");
    xor_all_ones(circ, block, &all_ones);

    circ.cx(ctrl, &carry);
    circ.ccx(&dec_ctrl, &all_ones, &carry);
    // carry's last touch is the ccx above. Free before consume allocs.
    drop(carry);

    xor_all_ones_consume_free(circ, block, all_ones);

    circ.cx(prev_carry, &dec_ctrl);
    circ.cx(ctrl, &dec_ctrl);
    drop(dec_ctrl);
}

/// Run-based controlled add for constants with clustered `1` bits.
///
/// This is a true chained carry automaton over alternating 1-runs and
/// 0-runs inside the active window `[runs[0].0, runs.last().1]`.
/// It pays one final suffix increment, not one suffix increment per
/// run. That is the whole point of this backend.
#[doc(hidden)]
pub fn controlled_add_const_runs_forced(circ: &mut Circuit, ctrl: &QReg, a: &[QReg], val: &[u8]) {
    use crate::point_add::trailmix_port::arith::khattar_gidney::cinc_khattar_gidney;

    let runs = one_runs(val, a.len());
    if runs.is_empty() {
        return;
    }

    let mut steps: Vec<RunCarryStep> = Vec::new();
    // All carry QRegs live in this Vec until they're consumed in the
    // uncompute pass. Steps reference into it by index.
    let mut carries: Vec<Option<QReg>> = Vec::new();

    let (first_lo, first_hi) = runs[0];
    let first_block = &a[first_lo..=first_hi];
    apply_conditional_decrement(circ, first_block, ctrl);

    let need_first_carry = runs.len() > 1 || first_hi + 1 < a.len();
    let mut prev_carry_idx: Option<usize> = if need_first_carry {
        let carry = circ.alloc_qreg("run_carry");
        carry_after_first_one_run(circ, ctrl, first_block, &carry);
        let idx = carries.len();
        carries.push(Some(carry));
        steps.push(RunCarryStep::FirstOne {
            lo: first_lo,
            hi: first_hi,
            carry_idx: idx,
        });
        Some(idx)
    } else {
        None
    };

    for pair in runs.windows(2) {
        let (prev_lo, prev_hi) = pair[0];
        let (lo, hi) = pair[1];
        let zero_lo = prev_hi + 1;
        let zero_hi = lo - 1;
        debug_assert!(zero_lo <= zero_hi);
        let carry_in_idx = prev_carry_idx.expect("carry missing before zero-run");

        let zero_block = &a[zero_lo..=zero_hi];
        {
            let carry_in = carries[carry_in_idx].as_ref().expect("carry alive");
            cinc_khattar_gidney(circ, zero_block, carry_in);
        }
        let zero_carry = circ.alloc_qreg("run_carry");
        {
            let carry_in = carries[carry_in_idx].as_ref().expect("carry alive");
            carry_after_zero_run(circ, carry_in, zero_block, &zero_carry);
        }
        let zero_carry_idx = carries.len();
        carries.push(Some(zero_carry));
        steps.push(RunCarryStep::Zero {
            lo: zero_lo,
            hi: zero_hi,
            prev_carry_idx: carry_in_idx,
            carry_idx: zero_carry_idx,
        });

        let one_block = &a[lo..=hi];
        let dec_ctrl = circ.alloc_qreg("run_dec_ctrl");
        circ.cx(ctrl, &dec_ctrl);
        {
            let zero_carry_ref = carries[zero_carry_idx].as_ref().expect("carry alive");
            circ.cx(zero_carry_ref, &dec_ctrl);
        }
        apply_conditional_decrement(circ, one_block, &dec_ctrl);
        {
            let zero_carry_ref = carries[zero_carry_idx].as_ref().expect("carry alive");
            circ.cx(zero_carry_ref, &dec_ctrl);
        }
        circ.cx(ctrl, &dec_ctrl);
        drop(dec_ctrl);

        let need_one_carry = hi + 1 < a.len();
        prev_carry_idx = if need_one_carry {
            let one_carry = circ.alloc_qreg("run_carry");
            {
                let zero_carry_ref = carries[zero_carry_idx].as_ref().expect("carry alive");
                carry_after_one_run(circ, ctrl, zero_carry_ref, one_block, &one_carry);
            }
            let idx = carries.len();
            carries.push(Some(one_carry));
            steps.push(RunCarryStep::One {
                lo,
                hi,
                prev_carry_idx: zero_carry_idx,
                carry_idx: idx,
            });
            Some(idx)
        } else {
            None
        };

        let _ = prev_lo;
    }

    let last_hi = runs.last().unwrap().1;
    if let Some(carry_idx) = prev_carry_idx {
        let carry_ref = carries[carry_idx].as_ref().expect("carry alive");
        cinc_khattar_gidney(circ, &a[last_hi + 1..], carry_ref);
    }

    for step in steps.into_iter().rev() {
        match step {
            RunCarryStep::FirstOne { lo, hi, carry_idx } => {
                let carry = carries[carry_idx].take().expect("carry alive");
                carry_after_first_one_run_uncompute_free_carry(circ, ctrl, &a[lo..=hi], carry);
            }
            RunCarryStep::Zero {
                lo,
                hi,
                prev_carry_idx,
                carry_idx,
            } => {
                let carry = carries[carry_idx].take().expect("carry alive");
                // Use the uncompute variant that frees carry before
                // xor_all_zeros_consume_free allocs intermediate ancillae.
                let prev_ref = carries[prev_carry_idx].as_ref().expect("prev carry alive");
                // We can't pass `prev_ref` directly because we need the
                // uncompute callee to take ownership of `carry`; pass a
                // borrow of prev (still alive at this point in the chain).
                carry_after_zero_run_uncompute_free_carry(circ, prev_ref, &a[lo..=hi], carry);
            }
            RunCarryStep::One {
                lo,
                hi,
                prev_carry_idx,
                carry_idx,
            } => {
                let carry = carries[carry_idx].take().expect("carry alive");
                let prev_ref = carries[prev_carry_idx].as_ref().expect("prev carry alive");
                carry_after_one_run_uncompute_free_carry(circ, ctrl, prev_ref, &a[lo..=hi], carry);
            }
        }
    }
    // Any remaining carries in the Vec (none should remain after the
    // reverse pass consumed them all) drop here.
    drop(carries);
}

/// Windowed controlled add-constant. Requires val's set bits to all
/// lie in [lo, hi]. Does one `classq_add` on the window, computes the
/// carry-out via a classical compare witness, then propagates via a
/// single controlled increment on a[hi+1..].
pub fn controlled_add_const_windowed(
    circ: &mut Circuit,
    ctrl: &QReg,
    a: &[QReg],
    val: &[u8],
    lo: usize,
    hi: usize,
) {
    use crate::point_add::trailmix_port::arith::khattar_gidney::{
        cinc_khattar_gidney, compare_geq_theorem3, controlled_classical_quantum_add,
    };
    let n = a.len();
    assert!(hi < n);

    let window = hi - lo + 1;

    // Build c_window: val shifted down by lo, masked to window bits.
    let mut c_window = vec![0u8; window.div_ceil(8)];
    for i in 0..window {
        if get_const_bit(val, lo + i) {
            c_window[i / 8] |= 1u8 << (i % 8);
        }
    }

    // Step 1: a[lo..=hi] += ctrl * c_window (mod 2^window).
    controlled_classical_quantum_add(circ, ctrl, &a[lo..=hi], &c_window);

    // Step 2: carry_out = ctrl AND (a[lo..=hi]_new < c_window).
    //   (Because a_new = a_old + ctrl*c mod 2^window, and overflow
    //    happened iff a_old + ctrl*c >= 2^window iff a_new < ctrl*c;
    //    for ctrl=0 this is a_new < 0 = false.)
    let v = circ.alloc_qreg("rfold_win_v");
    compare_geq_theorem3(circ, &a[lo..=hi], &c_window, &v);
    circ.x(&v); // v = (a[lo..=hi] < c_window).

    let carry = circ.alloc_qreg("rfold_win_c");
    circ.ccx(&v, ctrl, &carry);

    // Step 3: propagate the carry into a[hi+1..].
    if hi + 1 < n {
        cinc_khattar_gidney(circ, &a[hi + 1..], &carry);
    }

    // Uncompute carry and v.
    circ.ccx(&v, ctrl, &carry);
    drop(carry); // last touch was ccx above; drain at next gate (gap=0).

    circ.x(&v);
    compare_geq_theorem3(circ, &a[lo..=hi], &c_window, &v);
    // v drops here; drain fires at next gate (gap=0).
}

/// Sparse-constant controlled add: iterates over set bits of `val`,
/// emits a `cinc_khattar_gidney` from each bit position.
///
/// Semantic: a += ctrl * val (mod 2^n).
///
/// Cost: sum over set bits i of cinc(n-i), still O(popcount*n) in the
/// worst case but with a much smaller constant than Theorem 4.
pub fn controlled_add_const_sparse(circ: &mut Circuit, ctrl: &QReg, a: &[QReg], val: &[u8]) {
    let n = a.len();
    if n == 0 {
        return;
    }
    for i in 0..n {
        if get_const_bit(val, i) {
            // Add 2^i to a, controlled by ctrl: increment a[i..] by 1
            // conditioned on ctrl.
            crate::point_add::trailmix_port::arith::khattar_gidney::cinc_khattar_gidney(circ, &a[i..], ctrl);
        }
    }
}

pub fn controlled_sub_const(circ: &mut Circuit, ctrl: &QReg, a: &[QReg], val: &[u8]) {
    let n = a.len();
    if n == 0 {
        return;
    }
    // Compute -val mod 2^n = ~val + 1 (n-bit two's complement) so we
    // can subtract via a single controlled_add_const. The X-sandwich
    // form (~a; add_const; ~a) leaves bits of `a` untouched by the
    // inner add_const with cancelling X's at start and end — the
    // redundant-op detector flags those as wasted gates.
    let mut neg_bits = vec![false; n];
    let mut carry = true;
    for i in 0..n {
        let inv = !get_const_bit(val, i);
        neg_bits[i] = inv ^ carry;
        carry = inv && carry;
    }
    let mut neg_val = vec![0u8; n.div_ceil(8)];
    for i in 0..n {
        if neg_bits[i] {
            neg_val[i / 8] |= 1u8 << (i % 8);
        }
    }
    controlled_add_const(circ, ctrl, a, &neg_val);
}

/// Reference-slice variant of [`controlled_add_const`].
///
/// Routes directly to the dense `controlled_classical_quantum_add_refs`
/// (Theorem 5) implementation. The dispatch heuristics in the
/// `&[QReg]`-shaped variant (sparse, runs-forced, windowed) are
/// performance optimizations for known-shape constants; for the
/// view-shaped path we just use the always-correct cqadd path.
/// Cost: O(n log^2 n), polylog peak ancs.
pub fn controlled_add_const_refs(circ: &mut Circuit, ctrl: &QReg, a: &[&QReg], val: &[u8]) {
    let n = a.len();
    if n == 0 {
        return;
    }
    crate::point_add::trailmix_port::arith::khattar_gidney::controlled_classical_quantum_add_refs(circ, ctrl, a, val);
}

/// Reference-slice variant of [`controlled_sub_const`].
///
/// X-sandwich form: a := a + (-val mod 2^n) = a - val.
pub fn controlled_sub_const_refs(circ: &mut Circuit, ctrl: &QReg, a: &[&QReg], val: &[u8]) {
    let n = a.len();
    if n == 0 {
        return;
    }
    let mut neg_bits = vec![false; n];
    let mut carry = true;
    for i in 0..n {
        let inv = !get_const_bit(val, i);
        neg_bits[i] = inv ^ carry;
        carry = inv && carry;
    }
    let mut neg_val = vec![0u8; n.div_ceil(8)];
    for i in 0..n {
        if neg_bits[i] {
            neg_val[i / 8] |= 1u8 << (i % 8);
        }
    }
    controlled_add_const_refs(circ, ctrl, a, &neg_val);
}
