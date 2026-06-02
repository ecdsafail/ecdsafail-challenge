//! (refactor) Mechanically extracted from kaliski.rs. No logic changes.
use super::*;
/// Phase-clean variant of [`mul_by_const_acc`].  It uses exact Cuccaro based
/// add/double/halve blocks rather than the measurement-based fast variants.
/// This is too costly for production, but useful as an algebra-validating
/// fallback when the fast constant multiplier introduces alt-seed phase.
pub(crate) fn mul_by_const_acc_phase_clean(
    b: &mut B,
    x: &[QubitId],
    c: U256,
    acc: &[QubitId],
    p: U256,
    subtract: bool,
) {
    mul_by_const_acc_impl(b, x, c, acc, p, subtract, false, false);
}

/// Mixed variant for diagnosing the prescaler phase: exact q-q add/sub at the
/// sparse constant bits, but fast modular double/halve to walk between bit
/// positions.  If this is phase-clean, the culprit is the fast q-q add/sub, not
/// the scale-walk itself.
pub(crate) fn mul_by_const_acc_exact_adds_fast_shifts(
    b: &mut B,
    x: &[QubitId],
    c: U256,
    acc: &[QubitId],
    p: U256,
    subtract: bool,
) {
    mul_by_const_acc_impl(b, x, c, acc, p, subtract, false, true);
}

pub(crate) fn shift_tmp_up_for_sparse_const(
    b: &mut B,
    tmp: &[QubitId],
    p: U256,
    mut delta: usize,
    undo: &mut Vec<SparseConstShiftUndo>,
) {
    while delta >= 22 {
        let (spill, flag_inv, ovf) = mod_shift_left_by_k(b, tmp, p, 22);
        undo.push(SparseConstShiftUndo::Chunk(22, spill, flag_inv, ovf));
        delta -= 22;
    }
    if delta >= 12 {
        let (spill, flag_inv, ovf) = mod_shift_left_by_k(b, tmp, p, delta);
        undo.push(SparseConstShiftUndo::Chunk(delta, spill, flag_inv, ovf));
    } else if delta > 0 {
        for _ in 0..delta {
            mod_double_inplace_fast(b, tmp, p);
        }
        undo.push(SparseConstShiftUndo::Doubles(delta));
    }
}

pub(crate) fn undo_sparse_const_shifts(b: &mut B, tmp: &[QubitId], p: U256, undo: Vec<SparseConstShiftUndo>) {
    for item in undo.into_iter().rev() {
        match item {
            SparseConstShiftUndo::Doubles(k) => {
                for _ in 0..k {
                    mod_halve_inplace_fast(b, tmp, p);
                }
            }
            SparseConstShiftUndo::Chunk(k, spill, flag_inv, ovf) => {
                mod_shift_right_by_k(b, tmp, p, k, spill, flag_inv, ovf);
            }
        }
    }
}

/// `acc ±= x * c mod p` using exact q-q add/sub at sparse constant bits, but
/// jumping between distant bit positions with the Solinas k-bit shifter instead
/// of one modular double per zero bit.  This borrows `x` itself as the moving
/// 2^i*x lane and restores it before returning, removing the field-sized tmp
/// register from prescaled Kaliski initialization.
pub(crate) fn mul_by_const_acc_chunked_shifts_inplace_src(
    b: &mut B,
    x: &[QubitId],
    c: U256,
    acc: &[QubitId],
    p: U256,
    subtract: bool,
) {
    if c == U256::ZERO {
        return;
    }

    let mut positions = Vec::new();
    for i in 0..256 {
        if bit(c, i) {
            positions.push(i);
        }
    }

    let mut undo = Vec::new();
    let mut cur = 0usize;
    for pos in positions {
        shift_tmp_up_for_sparse_const(b, x, p, pos - cur, &mut undo);
        cur = pos;
        if subtract {
            mod_sub_qq(b, acc, x, p);
        } else {
            mod_add_qq(b, acc, x, p);
        }
    }

    undo_sparse_const_shifts(b, x, p, undo);
}

pub(crate) fn mul_by_const_acc_impl(
    b: &mut B,
    x: &[QubitId],
    c: U256,
    acc: &[QubitId],
    p: U256,
    subtract: bool,
    fast_adds: bool,
    fast_shifts: bool,
) {
    let n = x.len();
    if c == U256::ZERO {
        return;
    }

    // tmp := x  (via CX copy)
    let tmp = b.alloc_qubits(n);
    for i in 0..n {
        b.cx(x[i], tmp[i]);
    }

    // Iterate bits of c from LSB to MSB. At step i, tmp holds x * 2^i mod p.
    // Add tmp to acc if bit i of c is set. Then double tmp for the next step.
    //
    // We iterate up through the highest set bit of c, plus any trailing zero
    // bits (we must double enough times to make uncomputation clean).
    let mut top = 0usize;
    for i in 0..256 {
        if bit(c, i) {
            top = i;
        }
    }

    for i in 0..=top {
        if bit(c, i) {
            if fast_adds {
                if subtract {
                    mod_sub_qq_fast(b, acc, &tmp, p);
                } else {
                    mod_add_qq_fast(b, acc, &tmp, p);
                }
            } else if subtract {
                mod_sub_qq(b, acc, &tmp, p);
            } else {
                mod_add_qq(b, acc, &tmp, p);
            }
        }
        if i < top {
            if fast_shifts {
                mod_double_inplace_fast(b, &tmp, p);
            } else {
                mod_double_inplace(b, &tmp, p);
            }
        }
    }

    // At this point tmp = x * 2^top mod p. Halve it back `top` times to
    // recover x, then uncompute via cx.
    for _ in 0..top {
        if fast_shifts {
            mod_halve_inplace_fast(b, &tmp, p);
        } else {
            mod_halve_inplace(b, &tmp, p);
        }
    }
    for i in 0..n {
        b.cx(x[i], tmp[i]);
    }
    b.free_vec(&tmp);
}

pub(crate) fn kaliski_forward_with_coeff_caps(
    b: &mut B,
    v_in: &[QubitId],
    st: &KaliskiState,
    p: U256,
    iters: usize,
    coeff: Option<(&[QubitId], &[QubitId])>,
    bulk_caps: BulkPrefixCaps,
) {
    let n = v_in.len();
    debug_assert!(iters <= st.m_hist.len());
    if let Some((cr, cs)) = coeff {
        assert_eq!(cr.len(), n);
        assert_eq!(cs.len(), n);
    }

    // ─── Init ───
    // u := p (classical load)
    for i in 0..n {
        if bit(p, i) {
            b.x(st.u[i]);
        }
    }
    // v_w := v_in  (CX-copy; v_in unchanged)
    for i in 0..n {
        b.cx(v_in[i], st.v_w[i]);
    }
    // s := 1
    b.x(st.s[0]);
    // f := 1
    b.x(st.f_flag);

    // ─── Iterations ───
    let use_bulk_prefix3 = bulk_prefix_enabled();
    let mut frame: Option<QubitId> = None;
    for i in 0..iters {
        let is_last = i + 1 == iters;
        if use_bulk_prefix3 && i < bulk_caps.forward {
            kaliski_iteration_bulk_prefix3(
                b,
                p,
                &st.u,
                &st.v_w,
                &st.r,
                &st.s,
                st.m_hist[i],
                &st.m_hist[i + 1..iters],
                i,
                coeff,
                &mut frame,
                is_last,
            );
        } else {
            kaliski_iteration(
                b,
                p,
                &st.u,
                &st.v_w,
                &st.r,
                &st.s,
                st.m_hist[i],
                &st.m_hist[i + 1..iters],
                st.f_flag,
                i,
                coeff,
                &mut frame,
                is_last,
            );
        }
    }
    debug_assert!(frame.is_none(), "kaliski forward frame not consumed");

    // After the loop for nonzero v_in, classical invariants give:
    //   u = 1, v_w = 0, f = 0, a = b = add = 0
    //   r = raw coefficient (the NEGATIVE form: r = -v^{-1} * 2^{2n} mod p)
    //   s = some coefficient
    // We skip the `x(r); add_nbit_const(r, p+1)` negation (~2n CCX per call,
    // 4 calls total ≈ 8n Toffoli saved). Callers compensate by using the
    // negated inv: body multiplications that would normally `mul_add` with
    // +inv become `mul_sub` with -inv, and vice versa.
}

pub(crate) fn kaliski_backward_caps(
    b: &mut B,
    v_in: &[QubitId],
    st: &KaliskiState,
    p: U256,
    iters: usize,
    bulk_caps: BulkPrefixCaps,
) {
    let n = v_in.len();
    debug_assert!(iters <= st.m_hist.len());

    let use_bulk_prefix3 = bulk_prefix_enabled();
    // ─── Reverse iterations (in reverse order) ───
    let mut frame: Option<QubitId> = None;
    for i in (0..iters).rev() {
        let is_last = i + 1 == iters;
        if use_bulk_prefix3 && i < bulk_caps.backward {
            kaliski_iteration_bulk_prefix3_backward(
                b,
                p,
                &st.u,
                &st.v_w,
                &st.r,
                &st.s,
                st.m_hist[i],
                &st.m_hist[i + 1..iters],
                i,
                &mut frame,
                is_last,
            );
        } else {
            kaliski_iteration_backward(
                b,
                p,
                &st.u,
                &st.v_w,
                &st.r,
                &st.s,
                st.m_hist[i],
                &st.m_hist[i + 1..iters],
                st.f_flag,
                i,
                &mut frame,
                is_last,
            );
        }
    }
    debug_assert!(frame.is_none(), "kaliski backward frame not consumed");

    // ─── Reverse Init ───
    b.x(st.f_flag);
    b.x(st.s[0]);
    for i in 0..n {
        b.cx(v_in[i], st.v_w[i]);
    }
    for i in 0..n {
        if bit(p, i) {
            b.x(st.u[i]);
        }
    }
}

/// Run `body` with `inv` holding `v_in^{-1} mod p`, leaving `v_in`
/// unchanged. Allocates the kaliski state and `inv` register itself, then
/// frees them at the end. The body must NOT touch `st` or `v_in`.
///
/// Implementation keeps `st` live across the body, so we only run
/// `kaliski_forward` ONCE (and its emit_inverse once), instead of the
/// 4-call structure of the previous Bennett-cleaned `kal_compute_into`.
/// Halves the dominant kaliski cost.
pub(crate) fn emit_inverse_hmr_safe<F: FnOnce(&mut B)>(b: &mut B, f: F) {
    let start = b.ops.len();
    f(b);
    let end = b.ops.len();
    let fwd: Vec<_> = b.ops[start..end].to_vec();
    b.ops.truncate(start);
    for op in fwd.into_iter().rev() {
        match op.kind {
            OperationType::X
            | OperationType::Z
            | OperationType::CX
            | OperationType::CZ
            | OperationType::CCX
            | OperationType::CCZ
            | OperationType::Swap => b.ops.push(op),
            OperationType::R
            | OperationType::Hmr
            | OperationType::Register
            | OperationType::AppendToRegister
            | OperationType::DebugPrint => {}
            _ => panic!(
                "emit_inverse_hmr_safe: non-invertible op kind {:?} inside forward block",
                op.kind
            ),
        }
    }
}

pub(crate) fn with_kal_inv_raw<F: FnOnce(&mut B, &[QubitId])>(
    b: &mut B,
    v_in: &[QubitId],
    p: U256,
    iters: usize,
    body: F,
) {
    with_kal_inv_raw_coeff_caps(b, v_in, p, iters, None, bulk_prefix_caps(KalPair::Default), body);
}

pub(crate) fn with_kal_inv_raw_pair<F: FnOnce(&mut B, &[QubitId])>(
    b: &mut B,
    v_in: &[QubitId],
    p: U256,
    iters: usize,
    pair: KalPair,
    body: F,
) {
    with_kal_inv_raw_coeff_caps(b, v_in, p, iters, None, bulk_prefix_caps(pair), body);
}

pub(crate) fn kaliski_forward_alias_v_w_caps(
    b: &mut B,
    st: &KaliskiState,
    p: U256,
    iters: usize,
    bulk_caps: BulkPrefixCaps,
) {
    let n = st.v_w.len();
    debug_assert!(iters <= st.m_hist.len());

    for i in 0..n {
        if bit(p, i) {
            b.x(st.u[i]);
        }
    }
    b.x(st.s[0]);
    b.x(st.f_flag);

    let use_bulk_prefix3 = bulk_prefix_enabled();
    let mut frame: Option<QubitId> = None;
    for i in 0..iters {
        let is_last = i + 1 == iters;
        if use_bulk_prefix3 && i < bulk_caps.forward {
            kaliski_iteration_bulk_prefix3(
                b,
                p,
                &st.u,
                &st.v_w,
                &st.r,
                &st.s,
                st.m_hist[i],
                &st.m_hist[i + 1..iters],
                i,
                None,
                &mut frame,
                is_last,
            );
        } else {
            kaliski_iteration(
                b,
                p,
                &st.u,
                &st.v_w,
                &st.r,
                &st.s,
                st.m_hist[i],
                &st.m_hist[i + 1..iters],
                st.f_flag,
                i,
                None,
                &mut frame,
                is_last,
            );
        }
    }
    debug_assert!(frame.is_none(), "kaliski alias forward frame not consumed");
}

pub(crate) fn kaliski_backward_alias_v_w_caps(
    b: &mut B,
    st: &KaliskiState,
    p: U256,
    iters: usize,
    bulk_caps: BulkPrefixCaps,
) {
    debug_assert!(iters <= st.m_hist.len());

    let use_bulk_prefix3 = bulk_prefix_enabled();
    let mut frame: Option<QubitId> = None;
    for i in (0..iters).rev() {
        let is_last = i + 1 == iters;
        if use_bulk_prefix3 && i < bulk_caps.backward {
            kaliski_iteration_bulk_prefix3_backward(
                b,
                p,
                &st.u,
                &st.v_w,
                &st.r,
                &st.s,
                st.m_hist[i],
                &st.m_hist[i + 1..iters],
                i,
                &mut frame,
                is_last,
            );
        } else {
            kaliski_iteration_backward(
                b,
                p,
                &st.u,
                &st.v_w,
                &st.r,
                &st.s,
                st.m_hist[i],
                &st.m_hist[i + 1..iters],
                st.f_flag,
                i,
                &mut frame,
                is_last,
            );
        }
    }
    debug_assert!(frame.is_none(), "kaliski alias backward frame not consumed");

    b.x(st.f_flag);
    b.x(st.s[0]);
    for i in 0..st.u.len() {
        if bit(p, i) {
            b.x(st.u[i]);
        }
    }
}

pub(crate) fn with_kal_inv_raw_borrow_v_w_pair<F: FnOnce(&mut B, &[QubitId])>(
    b: &mut B,
    alias_v_w: &[QubitId],
    p: U256,
    iters: usize,
    pair: KalPair,
    body: F,
) {
    let n = alias_v_w.len();
    // Borrow the live denominator register as Kaliski's v_w. The callback must
    // not read or write alias_v_w: it is consumed to zero until backward restores it.
    //
    // DIALOG HISTORY-FOLD (KAL_DIALOG_FOLD=1): route foldable m_hist[j] into idle
    // HIGH bits of the borrowed v_w (alias_v_w[c]) instead of fresh qubits; only
    // the irreducible step-4 carry-pool slots stay dedicated.  Allocation ORDER
    // is preserved (u, r, s, m_hist, f_flag) so fold-OFF is byte-identical.
    let fold = kal_dialog_fold_enabled();
    let u = b.alloc_qubits(n);
    let r = b.alloc_qubits(n);
    let s = b.alloc_qubits(n);
    let mut dedicated_m_hist: Vec<QubitId> = Vec::new();
    let m_hist: Vec<QubitId> = if fold {
        let map = dialog_fold_vw_slotmap(n, iters);
        (0..iters)
            .map(|j| match map[j] {
                Some(c) => alias_v_w[c], // ride in the vacated v_w high bit
                None => {
                    let q = b.alloc_qubit();
                    dedicated_m_hist.push(q);
                    q
                }
            })
            .collect()
    } else {
        let v = b.alloc_qubits(iters);
        dedicated_m_hist = v.clone();
        v
    };
    let f_flag = b.alloc_qubit();
    let mut st = KaliskiState {
        u,
        v_w: alias_v_w.to_vec(),
        r,
        s,
        m_hist,
        f_flag,
    };
    let bulk_caps = bulk_prefix_caps(pair);
    let keep_full_state = std::env::var("KAL_KEEP_FULL_STATE").ok().as_deref() == Some("1");
    let keep_u = keep_full_state || std::env::var("KAL_KEEP_U").ok().as_deref() == Some("1");
    let free_s = !keep_full_state && std::env::var("KAL_FREE_S").ok().as_deref() != Some("0");

    kaliski_forward_alias_v_w_caps(b, &st, p, iters, bulk_caps);

    // Keep f_flag live across the body. Free/realloc of the terminal sentinel is
    // phase-fragile in alias envelopes.
    if !keep_u {
        b.x(st.u[0]);
        b.free_vec(&st.u);
    }
    if free_s {
        for i in 0..n {
            if bit(p, i) {
                b.x(st.s[i]);
            }
        }
        b.free_vec(&st.s);
    }

    let r_low: Vec<QubitId> = st.r[..n].to_vec();
    body(b, &r_low);

    if !keep_u {
        st.u = b.alloc_qubits(n);
        b.x(st.u[0]);
    }
    if free_s {
        st.s = b.alloc_qubits(n);
        for i in 0..n {
            if bit(p, i) {
                b.x(st.s[i]);
            }
        }
    }

    kaliski_backward_alias_v_w_caps(b, &st, p, iters, bulk_caps);

    b.free(st.f_flag);
    // Free only the DEDICATED m_hist qubits; folded slots are alias_v_w high
    // bits owned by the caller (restored to |0> by the backward m_i peel and to
    // the original tx value in their low bits). For fold-OFF this == st.m_hist.
    b.free_vec(&dedicated_m_hist);
    b.free_vec(&st.s);
    b.free_vec(&st.r);
    b.free_vec(&st.u);
}

pub(crate) fn kaliski_forward_prescaled_mixed(
    b: &mut B,
    v_in: &[QubitId],
    st: &KaliskiState,
    p: U256,
    iters: usize,
    scale: U256,
) {
    kaliski_forward_prescaled_kind(b, v_in, st, p, iters, scale, false);
}

pub(crate) fn kaliski_forward_prescaled_chunked(
    b: &mut B,
    v_in: &[QubitId],
    st: &KaliskiState,
    p: U256,
    iters: usize,
    scale: U256,
) {
    kaliski_forward_prescaled_kind(b, v_in, st, p, iters, scale, true);
}

pub(crate) fn kaliski_forward_prescaled_kind(
    b: &mut B,
    v_in: &[QubitId],
    st: &KaliskiState,
    p: U256,
    iters: usize,
    scale: U256,
    chunked: bool,
) {
    let n = v_in.len();
    debug_assert!(iters <= st.m_hist.len());

    for i in 0..n {
        if bit(p, i) {
            b.x(st.u[i]);
        }
    }
    if chunked {
        mul_by_const_acc_chunked_shifts_inplace_src(b, v_in, scale, &st.v_w, p, false);
    } else {
        mul_by_const_acc_exact_adds_fast_shifts(b, v_in, scale, &st.v_w, p, false);
    }
    b.x(st.s[0]);
    b.x(st.f_flag);

    let use_bulk_prefix3 = bulk_prefix_enabled();
    let bulk_prefix_iters = bulk_prefix_safe_iters();
    let mut frame: Option<QubitId> = None;
    for i in 0..iters {
        let is_last = i + 1 == iters;
        if use_bulk_prefix3 && i < bulk_prefix_iters {
            kaliski_iteration_bulk_prefix3(
                b,
                p,
                &st.u,
                &st.v_w,
                &st.r,
                &st.s,
                st.m_hist[i],
                &st.m_hist[i + 1..iters],
                i,
                None,
                &mut frame,
                is_last,
            );
        } else {
            kaliski_iteration(
                b,
                p,
                &st.u,
                &st.v_w,
                &st.r,
                &st.s,
                st.m_hist[i],
                &st.m_hist[i + 1..iters],
                st.f_flag,
                i,
                None,
                &mut frame,
                is_last,
            );
        }
    }
    debug_assert!(frame.is_none(), "kaliski prescaled forward frame not consumed");
}

pub(crate) fn kaliski_backward_prescaled_mixed(
    b: &mut B,
    v_in: &[QubitId],
    st: &KaliskiState,
    p: U256,
    iters: usize,
    scale: U256,
) {
    kaliski_backward_prescaled_kind(b, v_in, st, p, iters, scale, false);
}

pub(crate) fn kaliski_backward_prescaled_chunked(
    b: &mut B,
    v_in: &[QubitId],
    st: &KaliskiState,
    p: U256,
    iters: usize,
    scale: U256,
) {
    kaliski_backward_prescaled_kind(b, v_in, st, p, iters, scale, true);
}

pub(crate) fn kaliski_backward_prescaled_kind(
    b: &mut B,
    v_in: &[QubitId],
    st: &KaliskiState,
    p: U256,
    iters: usize,
    scale: U256,
    chunked: bool,
) {
    let n = v_in.len();
    debug_assert!(iters <= st.m_hist.len());

    let use_bulk_prefix3 = bulk_prefix_enabled();
    let bulk_prefix_iters = bulk_prefix_safe_iters();
    let mut frame: Option<QubitId> = None;
    for i in (0..iters).rev() {
        let is_last = i + 1 == iters;
        if use_bulk_prefix3 && i < bulk_prefix_iters {
            kaliski_iteration_bulk_prefix3_backward(
                b,
                p,
                &st.u,
                &st.v_w,
                &st.r,
                &st.s,
                st.m_hist[i],
                &st.m_hist[i + 1..iters],
                i,
                &mut frame,
                is_last,
            );
        } else {
            kaliski_iteration_backward(
                b,
                p,
                &st.u,
                &st.v_w,
                &st.r,
                &st.s,
                st.m_hist[i],
                &st.m_hist[i + 1..iters],
                st.f_flag,
                i,
                &mut frame,
                is_last,
            );
        }
    }
    debug_assert!(frame.is_none(), "kaliski prescaled backward frame not consumed");

    b.x(st.f_flag);
    b.x(st.s[0]);
    if chunked {
        mul_by_const_acc_chunked_shifts_inplace_src(b, v_in, scale, &st.v_w, p, true);
    } else {
        mul_by_const_acc_exact_adds_fast_shifts(b, v_in, scale, &st.v_w, p, true);
    }
    for i in 0..n {
        if bit(p, i) {
            b.x(st.u[i]);
        }
    }
}

pub(crate) fn with_kal_inv_raw_prescaled_mixed<F: FnOnce(&mut B, &[QubitId])>(
    b: &mut B,
    v_in: &[QubitId],
    p: U256,
    iters: usize,
    body: F,
) {
    with_kal_inv_raw_prescaled_kind(b, v_in, p, iters, false, body);
}

pub(crate) fn with_kal_inv_raw_prescaled_chunked<F: FnOnce(&mut B, &[QubitId])>(
    b: &mut B,
    v_in: &[QubitId],
    p: U256,
    iters: usize,
    body: F,
) {
    with_kal_inv_raw_prescaled_kind(b, v_in, p, iters, true, body);
}

pub(crate) fn with_kal_inv_raw_prescaled_kind<F: FnOnce(&mut B, &[QubitId])>(
    b: &mut B,
    v_in: &[QubitId],
    p: U256,
    iters: usize,
    chunked: bool,
    body: F,
) {
    let n = v_in.len();
    let mut st = alloc_kaliski_state(b, n, iters);
    let scale = pow_mod_2_k(p, iters);
    let keep_full_state = std::env::var("KAL_KEEP_FULL_STATE").ok().as_deref() == Some("1");
    let keep_u = keep_full_state || std::env::var("KAL_KEEP_U").ok().as_deref() == Some("1");
    let keep_v = keep_full_state || std::env::var("KAL_KEEP_V").ok().as_deref() == Some("1");
    let keep_f = keep_full_state || std::env::var("KAL_KEEP_F").ok().as_deref() == Some("1");
    let free_s = !keep_full_state && std::env::var("KAL_FREE_S").ok().as_deref() != Some("0");

    if chunked {
        kaliski_forward_prescaled_chunked(b, v_in, &st, p, iters, scale);
    } else {
        kaliski_forward_prescaled_mixed(b, v_in, &st, p, iters, scale);
    }

    if !keep_v {
        b.free_vec(&st.v_w);
    }
    if !keep_f {
        b.free(st.f_flag);
    }
    if !keep_u {
        b.x(st.u[0]);
        b.free_vec(&st.u);
    }
    if free_s {
        for i in 0..n {
            if bit(p, i) {
                b.x(st.s[i]);
            }
        }
        b.free_vec(&st.s);
    }

    let r_low: Vec<QubitId> = st.r[..n].to_vec();
    body(b, &r_low);

    if !keep_u {
        st.u = b.alloc_qubits(n);
        b.x(st.u[0]);
    }
    if !keep_f {
        st.f_flag = b.alloc_qubit();
    }
    if !keep_v {
        st.v_w = b.alloc_qubits(n);
    }
    if free_s {
        st.s = b.alloc_qubits(n);
        for i in 0..n {
            if bit(p, i) {
                b.x(st.s[i]);
            }
        }
    }

    if chunked {
        kaliski_backward_prescaled_chunked(b, v_in, &st, p, iters, scale);
    } else {
        kaliski_backward_prescaled_mixed(b, v_in, &st, p, iters, scale);
    }
    free_kaliski_state(b, st);
}

pub(crate) fn kaliski_xor_inv_raw_into_keep_alias_vw(
    b: &mut B,
    v_in: &[QubitId],
    alias_v_w: &[QubitId],
    p: U256,
    iters: usize,
    pair: KalPair,
    inv_keep: &[QubitId],
    caller_owns_v_w: bool,
) {
    let n = v_in.len();
    assert_eq!(alias_v_w.len(), n);
    assert_eq!(inv_keep.len(), n);
    let mut st = KaliskiState {
        u: b.alloc_qubits(n),
        v_w: alias_v_w.to_vec(),
        r: b.alloc_qubits(n),
        s: b.alloc_qubits(n),
        m_hist: b.alloc_qubits(iters),
        f_flag: b.alloc_qubit(),
    };
    let bulk_caps = cleanup_bulk_prefix_caps(pair);

    // H194/H199: mirror with_kal_inv_raw_coeff_caps's keep_u/keep_v/keep_f/free_s
    // envelope inside the cleanup helper so the forward Kaliski round-trip is
    // structurally identical to the production primary-helper round-trip.
    //
    // H199 bisect (attempt-198, this branch's 8-cell sweep) located the unique
    // envelope axis that closes the cleanup phase batches at both iters=0
    // (locator) and iters=374 (strict bulk-prefix3): `keep_u=false,
    // keep_f=true, free_s=false`.  Truth table (altseed_phase_batches_total):
    //
    //   (U,F,S)   iters=0   iters=374
    //   (0,0,0)     0          2
    //   (0,0,1)     0          1
    //   (0,1,0)     0          0   ← LOCKED DEFAULT
    //   (0,1,1)     0          0
    //   (1,0,0)     1          0
    //   (1,0,1)     0          1
    //   (1,1,0)     1          0
    //   (1,1,1)     0          2
    //
    // (0,1,0) and (0,1,1) are the only cells altseed-clean at BOTH iters=0
    // and iters=374; we pick (0,1,0) as the minimal-axis change (only
    // keep_f flips from the production-mirror default).  free_s is left
    // false (no `s` mutation in cleanup) and keep_u false (free `u` like
    // production).  caller_owns_v_w forces keep_v=true.
    //
    // env_keep_v always true because v_w aliases the caller-provided `ty`.
    let env_keep_u = std::env::var("KAL_PAIR1_INVKEEP_CLEANUP_ENV_KEEP_U")
        .ok()
        .as_deref()
        == Some("1");
    let env_keep_v = std::env::var("KAL_PAIR1_INVKEEP_CLEANUP_ENV_KEEP_V")
        .ok()
        .as_deref()
        != Some("0");
    // H199: default keep_f=true (the unique iters=374 closer); env override
    // wins so the bisect harness can still flip this.
    let env_keep_f = std::env::var("KAL_PAIR1_INVKEEP_CLEANUP_ENV_KEEP_F")
        .ok()
        .as_deref()
        .map(|s| s == "1")
        .unwrap_or(true);
    // H199: default free_s=false (no `s` mutation in cleanup); env override
    // wins.  (free_s=true is equivalent at iters=374 but adds 2n X-gates
    // around an alloc/realloc on `s`, so the minimal lock is false.)
    let env_free_s = std::env::var("KAL_PAIR1_INVKEEP_CLEANUP_ENV_FREE_S")
        .ok()
        .as_deref()
        .map(|s| s == "1")
        .unwrap_or(false);
    // When the helper uses emit_inverse_hmr_safe(forward) for the reverse
    // pass, forward and backward must see the SAME qubit ids; an envelope
    // that frees+reallocates would break this.  Disable when the user
    // requested generalized-reverse mode.
    let envelope_active = std::env::var("KAL_BULK3_GENERALIZED_REVERSE").is_err();
    // Honor alias contract: never free the caller-owned v_w.
    let keep_v_effective = env_keep_v || caller_owns_v_w;

    if std::env::var("TRACE_PHASE_LOCAL_PEAK")
        .ok()
        .map(|v| v.starts_with("pair1_invkeep") || v.starts_with("pair1_outside"))
        .unwrap_or(false)
    {
        eprintln!(
            "INVKEEP_CLEANUP_BULK_CAPS forward={} backward={}",
            bulk_caps.forward, bulk_caps.backward
        );
        eprintln!(
            "INVKEEP_CLEANUP_ENV keep_u={} keep_v={} keep_f={} free_s={} env_active={} caller_owns_v_w={}",
            env_keep_u, keep_v_effective, env_keep_f, env_free_s, envelope_active, caller_owns_v_w
        );
    }

    kaliski_forward_with_coeff_caps(b, v_in, &st, p, iters, None, bulk_caps);

    // Free envelope components between forward and backward, mirroring
    // with_kal_inv_raw_coeff_caps.  v_w is never freed here because it aliases
    // the caller's register (caller_owns_v_w guard).
    if envelope_active && !env_keep_u {
        // Forward end-state invariant: u[0] = 1, u[1..] = 0.  X-clear u[0]
        // then free.
        b.x(st.u[0]);
        b.free_vec(&st.u);
    }
    if envelope_active && !env_keep_f {
        b.free(st.f_flag);
    }
    if envelope_active && env_free_s {
        // Forward end-state invariant: s == p.  X-clear bits of p then free.
        for i in 0..n {
            if bit(p, i) {
                b.x(st.s[i]);
            }
        }
        b.free_vec(&st.s);
    }

    // Body: copy r_low into inv_keep via CNOTs (n-bit fan-out).  r is a
    // deterministic classical state at this point so the body is phase-free.
    for i in 0..n {
        b.cx(st.r[i], inv_keep[i]);
    }

    // Re-allocate envelope components before backward, exactly mirroring
    // production.  Note: st.v_w retains the alias; we never touch it.
    if envelope_active && !env_keep_u {
        st.u = b.alloc_qubits(n);
        b.x(st.u[0]);
    }
    if envelope_active && !env_keep_f {
        st.f_flag = b.alloc_qubit();
    }
    if envelope_active && env_free_s {
        st.s = b.alloc_qubits(n);
        for i in 0..n {
            if bit(p, i) {
                b.x(st.s[i]);
            }
        }
    }

    if std::env::var("KAL_BULK3_GENERALIZED_REVERSE").is_ok() {
        emit_inverse_hmr_safe(b, |b| {
            kaliski_forward_with_coeff_caps(b, v_in, &st, p, iters, None, bulk_caps)
        });
    } else {
        kaliski_backward_caps(b, v_in, &st, p, iters, bulk_caps);
    }
    b.free(st.f_flag);
    b.free_vec(&st.m_hist);
    b.free_vec(&st.s);
    b.free_vec(&st.r);
    if !caller_owns_v_w {
        b.free_vec(&st.v_w);
    }
    b.free_vec(&st.u);
}

pub(crate) fn with_kal_inv_raw_coeff<F: FnOnce(&mut B, &[QubitId])>(
    b: &mut B,
    v_in: &[QubitId],
    p: U256,
    iters: usize,
    coeff: Option<(&[QubitId], &[QubitId])>,
    body: F,
) {
    with_kal_inv_raw_coeff_caps(
        b,
        v_in,
        p,
        iters,
        coeff,
        bulk_prefix_caps(KalPair::Default),
        body,
    );
}


pub(crate) fn with_kal_inv_raw_coeff_caps<F: FnOnce(&mut B, &[QubitId])>(
    b: &mut B,
    v_in: &[QubitId],
    p: U256,
    iters: usize,
    coeff: Option<(&[QubitId], &[QubitId])>,
    bulk_caps: BulkPrefixCaps,
    body: F,
) {
    let n = v_in.len();
    let mut st = alloc_kaliski_state(b, n, iters);
    let keep_full_state = std::env::var("KAL_KEEP_FULL_STATE").ok().as_deref() == Some("1");
    let keep_u = keep_full_state || std::env::var("KAL_KEEP_U").ok().as_deref() == Some("1");
    let keep_v = keep_full_state || std::env::var("KAL_KEEP_V").ok().as_deref() == Some("1");
    let keep_f = keep_full_state || std::env::var("KAL_KEEP_F").ok().as_deref() == Some("1");
    // KAL_FREE_S=1 (default ON in this branch): at end of forward Kaliski,
    // the s register provably equals p (the modulus) when iters >= ~407
    // (verified classically for our specific Kaliski variant). Free s by
    // X-ing the bits of p, then re-load before backward.
    let free_s = !keep_full_state && std::env::var("KAL_FREE_S").ok().as_deref() != Some("0");

    // Forward kaliski. st.r[..n] holds raw = v_in^{-1} * 2^(2n) mod p.
    // If coeff is supplied, the same branch controls also transform that
    // external coefficient pair, but the ordinary qrisp sentinel state remains
    // available for clean branch-flag uncomputation.
    kaliski_forward_with_coeff_caps(b, v_in, &st, p, iters, coeff, bulk_caps);

    if !keep_v {
        b.free_vec(&st.v_w);
    }
    if !keep_f {
        b.free(st.f_flag);
    }
    if !keep_u {
        b.x(st.u[0]);
        b.free_vec(&st.u);
    }
    if free_s {
        // s = p at this point. X each bit of p to zero it.
        for i in 0..n {
            if bit(p, i) {
                b.x(st.s[i]);
            }
        }
        b.free_vec(&st.s);
    }

    let r_low: Vec<QubitId> = st.r[..n].to_vec();
    body(b, &r_low);

    if !keep_u {
        // Re-alloc at |0> for the backward pass; restore u[0] = 1.
        st.u = b.alloc_qubits(n);
        b.x(st.u[0]);
    }
    if !keep_f {
        st.f_flag = b.alloc_qubit();
    }
    if !keep_v {
        st.v_w = b.alloc_qubits(n);
    }
    if free_s {
        // Re-allocate s and load p back.
        st.s = b.alloc_qubits(n);
        for i in 0..n {
            if bit(p, i) {
                b.x(st.s[i]);
            }
        }
    }

    // Experimental mode: use the exact reversed forward block shape, but skip
    // HMR/R in the reverse replay. This is heavier than the explicit backward,
    // but it keeps the specialized prefix and its matching global reverse in a
    // single contract. The hope is to eliminate the residual phase mismatch.
    if std::env::var("KAL_BULK3_GENERALIZED_REVERSE").is_ok() {
        emit_inverse_hmr_safe(b, |b| {
            kaliski_forward_with_coeff_caps(b, v_in, &st, p, iters, None, bulk_caps)
        });
    } else {
        // Explicit backward pass (uses measurement-based uncompute, saves
        // ~511 CCX per iteration vs the emit_inverse version).  Use the same
        // promoted/pair-specific cap family selected for the forward pass so
        // a 378th bulk step can be enabled only where it is phase-clean.
        kaliski_backward_caps(b, v_in, &st, p, iters, bulk_caps);
    }

    free_kaliski_state(b, st);
}
