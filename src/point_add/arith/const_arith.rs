use super::*;

#[inline]
fn maj1_inputs_distinct(a: QubitId, k: QubitId, carry: QubitId, target: QubitId) -> bool {
    a != k && a != carry && a != target && k != carry && k != target && carry != target
}

#[inline]
fn fold_maj1_enabled() -> bool {
    std::env::var("DIALOG_GCD_FOLD_MAJ1").ok().as_deref() == Some("1")
}

fn emit_fold_maj1(b: &mut B, a: QubitId, k: QubitId, carry: QubitId, target: QubitId) {
    debug_assert!(maj1_inputs_distinct(a, k, carry, target));
    b.cx(carry, target);
    b.cx(carry, a);
    b.cx(carry, k);
    b.ccx(a, k, target);
    b.cx(carry, k);
    b.cx(carry, a);
}

fn emit_fold_majority(
    b: &mut B,
    a: QubitId,
    k: QubitId,
    carry: QubitId,
    target: QubitId,
    maj2: bool,
) {
    if fold_maj1_enabled() && maj1_inputs_distinct(a, k, carry, target) {
        emit_fold_maj1(b, a, k, carry, target);
    } else if maj2 {
        b.ccx(a, carry, target);
        b.cx(a, carry);
        b.ccx(k, carry, target);
        b.cx(a, carry);
    } else {
        b.ccx(a, carry, target);
        b.ccx(k, a, target);
        b.ccx(k, carry, target);
    }
}

pub(crate) fn csub_nbit_const(b: &mut B, acc: &[QubitId], c: U256, ctrl: QubitId) {
    let n = acc.len();
    let a = b.alloc_qubits(n);
    for i in 0..n {
        if bit(c, i) {
            b.cx(ctrl, a[i]);
        }
    }
    sub_nbit_qq(b, &a, acc);
    for i in 0..n {
        if bit(c, i) {
            b.cx(ctrl, a[i]);
        }
    }
    b.free_vec(&a);
}

pub(crate) fn cadd_nbit_const(b: &mut B, acc: &[QubitId], c: U256, ctrl: QubitId) {
    let n = acc.len();
    let a = b.alloc_qubits(n);
    for i in 0..n {
        if bit(c, i) {
            b.cx(ctrl, a[i]);
        }
    }
    add_nbit_qq(b, &a, acc);
    for i in 0..n {
        if bit(c, i) {
            b.cx(ctrl, a[i]);
        }
    }
    b.free_vec(&a);
}

pub(crate) fn csub_nbit_const_fast(b: &mut B, acc: &[QubitId], c: U256, ctrl: QubitId) {
    let n = acc.len();
    let a = b.alloc_qubits(n);
    for i in 0..n {
        if bit(c, i) {
            b.cx(ctrl, a[i]);
        }
    }
    sub_nbit_qq_fast(b, &a, acc);
    for i in 0..n {
        if bit(c, i) {
            b.cx(ctrl, a[i]);
        }
    }
    b.free_vec(&a);
}

pub(crate) fn csub_nbit_const_direct_fast(b: &mut B, acc: &[QubitId], c: U256, ctrl: QubitId) {
    let n = acc.len();
    if n == 0 {
        return;
    }
    if n == 1 {
        if bit(c, 0) {
            b.cx(ctrl, acc[0]);
        }
        return;
    }

    let borrows = b.alloc_qubits(n - 1);

    for i in 0..n - 1 {
        let target = borrows[i];
        let borrow_in = if i == 0 { None } else { Some(borrows[i - 1]) };
        if bit(c, i) {
            b.x(acc[i]);
            if let Some(bi) = borrow_in {
                emit_fold_majority(b, acc[i], ctrl, bi, target, false);
            } else {
                b.ccx(acc[i], ctrl, target);
            }
            b.x(acc[i]);
        } else if let Some(bi) = borrow_in {
            b.x(acc[i]);
            b.ccx(acc[i], bi, target);
            b.x(acc[i]);
        }
    }

    for i in 0..n {
        if bit(c, i) {
            b.cx(ctrl, acc[i]);
        }
        if i > 0 {
            b.cx(borrows[i - 1], acc[i]);
        }
    }

    for i in (0..n - 1).rev() {
        let m = b.alloc_bit();
        b.hmr(borrows[i], m);
        let borrow_in = if i == 0 { None } else { Some(borrows[i - 1]) };
        if bit(c, i) {
            if let Some(bi) = borrow_in {
                b.cz_if(acc[i], ctrl, m);
                b.cz_if(acc[i], bi, m);
                b.cz_if(ctrl, bi, m);
            } else {
                b.cz_if(acc[i], ctrl, m);
            }
        } else if let Some(bi) = borrow_in {
            b.cz_if(acc[i], bi, m);
        }
    }

    b.free_vec(&borrows);
}

pub(crate) fn cadd_nbit_const_fast(b: &mut B, acc: &[QubitId], c: U256, ctrl: QubitId) {
    let n = acc.len();
    let a = b.alloc_qubits(n);
    for i in 0..n {
        if bit(c, i) {
            b.cx(ctrl, a[i]);
        }
    }
    add_nbit_qq_fast(b, &a, acc);
    for i in 0..n {
        if bit(c, i) {
            b.cx(ctrl, a[i]);
        }
    }
    b.free_vec(&a);
}

pub(crate) fn cadd_nbit_const_direct_fast(b: &mut B, acc: &[QubitId], c: U256, ctrl: QubitId) {
    let n = acc.len();
    if n == 0 {
        return;
    }
    if n == 1 {
        if bit(c, 0) {
            b.cx(ctrl, acc[0]);
        }
        return;
    }

    let carries = b.alloc_qubits(n - 1);

    for i in 0..n - 1 {
        let target = carries[i];
        let carry_in = if i == 0 { None } else { Some(carries[i - 1]) };
        if bit(c, i) {
            if let Some(ci) = carry_in {
                emit_fold_majority(b, acc[i], ctrl, ci, target, false);
            } else {
                b.ccx(acc[i], ctrl, target);
            }
        } else if let Some(ci) = carry_in {
            b.ccx(acc[i], ci, target);
        }
    }

    for i in 0..n {
        if bit(c, i) {
            b.cx(ctrl, acc[i]);
        }
        if i > 0 {
            b.cx(carries[i - 1], acc[i]);
        }
    }

    for i in (0..n - 1).rev() {
        let m = b.alloc_bit();
        b.hmr(carries[i], m);
        let carry_in = if i == 0 { None } else { Some(carries[i - 1]) };
        if bit(c, i) {
            b.x(acc[i]);
            if let Some(ci) = carry_in {
                b.cz_if(acc[i], ctrl, m);
                b.cz_if(acc[i], ci, m);
                b.x(acc[i]);
                b.cz_if(ctrl, ci, m);
            } else {
                b.cz_if(acc[i], ctrl, m);
                b.x(acc[i]);
            }
        } else if let Some(ci) = carry_in {
            b.x(acc[i]);
            b.cz_if(acc[i], ci, m);
            b.x(acc[i]);
        }
    }

    b.free_vec(&carries);
}

pub(crate) fn add_nbit_const_extcarry_clean(b: &mut B, acc_ext: &[QubitId], c: U256) {
    add_nbit_const_extcarry_clean_with_cin(b, acc_ext, c, None);
}

pub(crate) fn add_nbit_const_extcarry_clean_with_cin(
    b: &mut B,
    acc_ext: &[QubitId],
    c: U256,
    borrow_cin: Option<QubitId>,
) {
    let ext = acc_ext.len();
    debug_assert!(ext >= 1);
    let n = ext - 1;
    let ca = load_const(b, n, c);
    let (c_in, fresh) = match borrow_cin {
        Some(q) => (q, false),
        None => (b.alloc_qubit(), true),
    };
    cuccaro_add_low_to_ext_clean(b, &ca, acc_ext, c_in);
    if fresh {
        b.free(c_in);
    }
    unload_const(b, &ca, c);
}

pub(crate) fn sub_nbit_const_extcarry_clean(b: &mut B, acc_ext: &[QubitId], c: U256) {
    let ext = acc_ext.len();
    debug_assert!(ext >= 1);
    let n = ext - 1;
    let ca = load_const(b, n, c);
    let c_in = b.alloc_qubit();
    cuccaro_sub_low_to_ext_clean(b, &ca, acc_ext, c_in);
    b.free(c_in);
    unload_const(b, &ca, c);
}

pub(crate) fn cadd_nbit_const_extcarry_clean(
    b: &mut B,
    acc_ext: &[QubitId],
    c: U256,
    ctrl: QubitId,
) {
    let ext = acc_ext.len();
    debug_assert!(ext >= 1);
    let n = ext - 1;
    let ca = b.alloc_qubits(n);
    for i in 0..n {
        if bit(c, i) {
            b.cx(ctrl, ca[i]);
        }
    }
    let c_in = b.alloc_qubit();
    cuccaro_add_low_to_ext_clean(b, &ca, acc_ext, c_in);
    b.free(c_in);
    for i in 0..n {
        if bit(c, i) {
            b.cx(ctrl, ca[i]);
        }
    }
    b.free_vec(&ca);
}

pub(crate) fn csub_nbit_const_extcarry_clean(
    b: &mut B,
    acc_ext: &[QubitId],
    c: U256,
    ctrl: QubitId,
) {
    csub_nbit_const_extcarry_clean_with_cin(b, acc_ext, c, ctrl, None);
}

pub(crate) fn csub_nbit_const_extcarry_clean_with_cin(
    b: &mut B,
    acc_ext: &[QubitId],
    c: U256,
    ctrl: QubitId,
    borrow_cin: Option<QubitId>,
) {
    let ext = acc_ext.len();
    debug_assert!(ext >= 1);
    let n = ext - 1;
    let ca = b.alloc_qubits(n);
    for i in 0..n {
        if bit(c, i) {
            b.cx(ctrl, ca[i]);
        }
    }
    let (c_in, fresh) = match borrow_cin {
        Some(q) => (q, false),
        None => (b.alloc_qubit(), true),
    };
    cuccaro_sub_low_to_ext_clean(b, &ca, acc_ext, c_in);
    if fresh {
        b.free(c_in);
    }
    for i in 0..n {
        if bit(c, i) {
            b.cx(ctrl, ca[i]);
        }
    }
    b.free_vec(&ca);
}

pub(crate) fn add_nbit_const_direct_uncontrolled_fast(b: &mut B, acc: &[QubitId], c: U256) {
    let ctrl = b.alloc_qubit();
    b.x(ctrl);
    cadd_nbit_const_direct_fast(b, acc, c, ctrl);
    b.x(ctrl);
    b.free(ctrl);
}

pub(crate) fn sub_nbit_const_direct_uncontrolled_fast(b: &mut B, acc: &[QubitId], c: U256) {
    let ctrl = b.alloc_qubit();
    b.x(ctrl);
    csub_nbit_const_direct_fast(b, acc, c, ctrl);
    b.x(ctrl);
    b.free(ctrl);
}

pub(crate) fn add_nbit_const_fast(b: &mut B, acc: &[QubitId], c: U256) {
    if secp_direct_const_arith_enabled() {
        add_nbit_const_direct_uncontrolled_fast(b, acc, c);
        return;
    }
    let n = acc.len();
    let a = load_const(b, n, c);
    add_nbit_qq_fast(b, &a, acc);
    unload_const(b, &a, c);
}

pub(crate) fn sub_nbit_const_fast(b: &mut B, acc: &[QubitId], c: U256) {
    if secp_direct_const_arith_enabled() {
        sub_nbit_const_direct_uncontrolled_fast(b, acc, c);
        return;
    }
    let n = acc.len();
    let a = load_const(b, n, c);
    sub_nbit_qq_fast(b, &a, acc);
    unload_const(b, &a, c);
}

pub(crate) fn highest_set_bit(c: U256) -> usize {
    let mut hi = 0usize;
    for i in 0..256 {
        if bit(c, i) {
            hi = i;
        }
    }
    hi
}

pub(crate) fn double_carry_trunc_window() -> Option<usize> {
    std::env::var("KAL_DOUBLE_CARRY_TRUNC_W")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .filter(|&w| w > 0)
}

pub(crate) fn fold_carry_trunc_window() -> Option<usize> {
    std::env::var("KAL_FOLD_CARRY_TRUNC_W")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .filter(|&w| w > 0)
}

pub(crate) fn perpos_maj2_enabled() -> bool {
    std::env::var("DIALOG_GCD_PERPOS_MAJ2").ok().as_deref() == Some("1")
}

pub(crate) fn fold_maj2_enabled() -> bool {
    std::env::var("DIALOG_GCD_FOLD_MAJ2").ok().as_deref() == Some("1")
}

fn borrowed_const_fold_carries(
    b: &mut B,
    need: usize,
    borrowed: &[QubitId],
) -> (Vec<QubitId>, Vec<QubitId>) {
    let borrowed_len = borrowed.len().min(need);
    let owned = b.alloc_qubits(need - borrowed_len);
    let mut carries = Vec::with_capacity(need);
    carries.extend_from_slice(&borrowed[..borrowed_len]);
    carries.extend_from_slice(&owned);
    (carries, owned)
}

pub(crate) fn cadd_nbit_const_direct_trunc_fast(
    b: &mut B,
    acc: &[QubitId],
    c: U256,
    ctrl: QubitId,
    window: usize,
) {
    cadd_nbit_const_direct_trunc_fast_borrowed_carries(b, acc, c, ctrl, window, &[]);
}

pub(crate) fn cadd_nbit_const_direct_trunc_fast_borrowed_carries(
    b: &mut B,
    acc: &[QubitId],
    c: U256,
    ctrl: QubitId,
    window: usize,
    borrowed_carries: &[QubitId],
) {
    let n = acc.len();
    if n == 0 {
        return;
    }
    if n == 1 {
        if bit(c, 0) {
            b.cx(ctrl, acc[0]);
        }
        return;
    }

    let hi = highest_set_bit(c);
    let last = core::cmp::min(n - 2, hi.saturating_add(window));
    let maj2 = fold_maj2_enabled();
    let (carries, owned_carries) = borrowed_const_fold_carries(b, last + 1, borrowed_carries);

    for i in 0..=last {
        let target = carries[i];
        let carry_in = if i == 0 { None } else { Some(carries[i - 1]) };
        if bit(c, i) {
            if let Some(ci) = carry_in {
                emit_fold_majority(b, acc[i], ctrl, ci, target, maj2);
            } else {
                b.ccx(acc[i], ctrl, target);
            }
        } else if let Some(ci) = carry_in {
            b.ccx(acc[i], ci, target);
        }
    }

    for i in 0..n {
        if bit(c, i) {
            b.cx(ctrl, acc[i]);
        }
        if i > 0 && i - 1 <= last {
            b.cx(carries[i - 1], acc[i]);
        }
    }

    for i in (0..=last).rev() {
        let m = b.alloc_bit();
        b.hmr(carries[i], m);
        let carry_in = if i == 0 { None } else { Some(carries[i - 1]) };
        if bit(c, i) {
            b.x(acc[i]);
            if let Some(ci) = carry_in {
                b.cz_if(acc[i], ctrl, m);
                b.cz_if(acc[i], ci, m);
                b.x(acc[i]);
                b.cz_if(ctrl, ci, m);
            } else {
                b.cz_if(acc[i], ctrl, m);
                b.x(acc[i]);
            }
        } else if let Some(ci) = carry_in {
            b.x(acc[i]);
            b.cz_if(acc[i], ci, m);
            b.x(acc[i]);
        }
    }

    b.free_vec(&owned_carries);
}

pub(crate) fn csub_nbit_const_direct_trunc_fast(
    b: &mut B,
    acc: &[QubitId],
    c: U256,
    ctrl: QubitId,
    window: usize,
) {
    csub_nbit_const_direct_trunc_fast_borrowed_carries(b, acc, c, ctrl, window, &[]);
}

pub(crate) fn csub_nbit_const_direct_trunc_fast_borrowed_carries(
    b: &mut B,
    acc: &[QubitId],
    c: U256,
    ctrl: QubitId,
    window: usize,
    borrowed_carries: &[QubitId],
) {
    let n = acc.len();
    if n == 0 {
        return;
    }
    if n == 1 {
        if bit(c, 0) {
            b.cx(ctrl, acc[0]);
        }
        return;
    }

    let hi = highest_set_bit(c);
    let last = core::cmp::min(n - 2, hi.saturating_add(window));
    let maj2 = fold_maj2_enabled();
    let (borrows, owned_borrows) = borrowed_const_fold_carries(b, last + 1, borrowed_carries);

    for i in 0..=last {
        let target = borrows[i];
        let borrow_in = if i == 0 { None } else { Some(borrows[i - 1]) };
        if bit(c, i) {
            b.x(acc[i]);
            if let Some(bi) = borrow_in {
                emit_fold_majority(b, acc[i], ctrl, bi, target, maj2);
            } else {
                b.ccx(acc[i], ctrl, target);
            }
            b.x(acc[i]);
        } else if let Some(bi) = borrow_in {
            b.x(acc[i]);
            b.ccx(acc[i], bi, target);
            b.x(acc[i]);
        }
    }

    for i in 0..n {
        if bit(c, i) {
            b.cx(ctrl, acc[i]);
        }
        if i > 0 && i - 1 <= last {
            b.cx(borrows[i - 1], acc[i]);
        }
    }

    for i in (0..=last).rev() {
        let m = b.alloc_bit();
        b.hmr(borrows[i], m);
        let borrow_in = if i == 0 { None } else { Some(borrows[i - 1]) };
        if bit(c, i) {
            if let Some(bi) = borrow_in {
                b.cz_if(acc[i], ctrl, m);
                b.cz_if(acc[i], bi, m);
                b.cz_if(ctrl, bi, m);
            } else {
                b.cz_if(acc[i], ctrl, m);
            }
        } else if let Some(bi) = borrow_in {
            b.cz_if(acc[i], bi, m);
        }
    }

    b.free_vec(&owned_borrows);
}

fn special_fold_park_low_carries() -> usize {
    std::env::var("DIALOG_GCD_SPECIAL_FOLD_PARK_LOW_CARRIES")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(0)
}

fn special_fold_park_low_carries_at_step(step: Option<usize>) -> usize {
    let mapped = step.and_then(|step| {
        let map = std::env::var("DIALOG_GCD_SPECIAL_FOLD_PARK_LOW_CARRIES_STEP_MAP").ok()?;
        map.split(',').rev().find_map(|entry| {
            let (raw_step, raw_value) = entry.trim().split_once(':')?;
            if raw_step.trim().parse::<usize>().ok()? != step {
                return None;
            }
            raw_value.trim().parse::<usize>().ok()
        })
    });
    mapped.unwrap_or_else(special_fold_park_low_carries)
}

fn cconst_nbit_direct_trunc_fast_parked(
    b: &mut B,
    acc: &[QubitId],
    c: U256,
    ctrl: QubitId,
    window: usize,
    park_low: usize,
    is_add: bool,
) {
    let n = acc.len();
    if n <= 1 {
        if n == 1 && bit(c, 0) {
            b.cx(ctrl, acc[0]);
        }
        return;
    }

    let hi = highest_set_bit(c);
    let last = core::cmp::min(n - 2, hi.saturating_add(window));
    let park_low = core::cmp::min(park_low, last.saturating_sub(hi));
    if park_low == 0 {
        if is_add {
            cadd_nbit_const_direct_trunc_fast(b, acc, c, ctrl, window);
        } else {
            csub_nbit_const_direct_trunc_fast(b, acc, c, ctrl, window);
        }
        return;
    }

    let split = last - park_low;
    let maj2 = fold_maj2_enabled();
    let prefix = b.alloc_qubits(split + 1);
    let kctrl = |i: usize| bit(c, i).then_some(ctrl);

    for i in 0..=split {
        let target = prefix[i];
        let carry_in = if i == 0 { None } else { Some(prefix[i - 1]) };
        if is_add {
            if let Some(kc) = kctrl(i) {
                if let Some(ci) = carry_in {
                    emit_fold_majority(b, acc[i], kc, ci, target, maj2);
                } else {
                    b.ccx(acc[i], kc, target);
                }
            } else if let Some(ci) = carry_in {
                b.ccx(acc[i], ci, target);
            }
        } else if let Some(kc) = kctrl(i) {
            b.x(acc[i]);
            if let Some(ci) = carry_in {
                emit_fold_majority(b, acc[i], kc, ci, target, maj2);
            } else {
                b.ccx(acc[i], kc, target);
            }
            b.x(acc[i]);
        } else if let Some(ci) = carry_in {
            b.x(acc[i]);
            b.ccx(acc[i], ci, target);
            b.x(acc[i]);
        }
    }

    for i in 0..=split {
        if let Some(kc) = kctrl(i) {
            b.cx(kc, acc[i]);
        }
        if i > 0 {
            b.cx(prefix[i - 1], acc[i]);
        }
    }

    for i in (0..park_low).rev() {
        let measured = b.alloc_bit();
        b.hmr(prefix[i], measured);
        let carry_in = if i == 0 { None } else { Some(prefix[i - 1]) };
        fold_postsum_carry_phase_uncompute(b, acc, kctrl(i), carry_in, measured, i, is_add);
        b.free(prefix[i]);
    }

    let tail = b.alloc_qubits(park_low);
    let carry = |i: usize| {
        if i <= split {
            prefix[i]
        } else {
            tail[i - split - 1]
        }
    };
    for i in split + 1..=last {
        let target = carry(i);
        let carry_in = carry(i - 1);
        if is_add {
            if let Some(kc) = kctrl(i) {
                emit_fold_majority(b, acc[i], kc, carry_in, target, maj2);
            } else {
                b.ccx(acc[i], carry_in, target);
            }
        } else if let Some(kc) = kctrl(i) {
            b.x(acc[i]);
            emit_fold_majority(b, acc[i], kc, carry_in, target, maj2);
            b.x(acc[i]);
        } else {
            b.x(acc[i]);
            b.ccx(acc[i], carry_in, target);
            b.x(acc[i]);
        }
    }

    for i in split + 1..n {
        if let Some(kc) = kctrl(i) {
            b.cx(kc, acc[i]);
        }
        if i - 1 <= last {
            b.cx(carry(i - 1), acc[i]);
        }
    }

    for i in (split + 1..=last).rev() {
        let measured = b.alloc_bit();
        b.hmr(carry(i), measured);
        fold_postsum_carry_phase_uncompute(
            b,
            acc,
            kctrl(i),
            Some(carry(i - 1)),
            measured,
            i,
            is_add,
        );
        b.free(carry(i));
    }
    drop(tail);

    for i in 0..park_low {
        b.reacquire(prefix[i]);
        let carry_in = if i == 0 { None } else { Some(prefix[i - 1]) };
        fold_postsum_carry_compute(b, acc, kctrl(i), carry_in, prefix[i], i, is_add);
    }

    for i in (0..=split).rev() {
        let measured = b.alloc_bit();
        b.hmr(prefix[i], measured);
        let carry_in = if i == 0 { None } else { Some(prefix[i - 1]) };
        fold_postsum_carry_phase_uncompute(b, acc, kctrl(i), carry_in, measured, i, is_add);
        b.free(prefix[i]);
    }
}

pub(crate) fn cadd_nbit_const_direct_trunc_fast_releasing_scratch(
    b: &mut B,
    acc: &[QubitId],
    c: U256,
    ctrl: QubitId,
    window: usize,
    releasable_scratch: &[QubitId],
) {
    cadd_nbit_const_direct_trunc_fast_releasing_scratch_at_step(
        b,
        acc,
        c,
        ctrl,
        window,
        releasable_scratch,
        None,
    );
}

pub(crate) fn cadd_nbit_const_direct_trunc_fast_releasing_scratch_at_step(
    b: &mut B,
    acc: &[QubitId],
    c: U256,
    ctrl: QubitId,
    window: usize,
    releasable_scratch: &[QubitId],
    step: Option<usize>,
) {
    let park_low = special_fold_park_low_carries_at_step(step);
    if park_low == 0 || releasable_scratch.is_empty() {
        cadd_nbit_const_direct_trunc_fast_borrowed_carries(
            b,
            acc,
            c,
            ctrl,
            window,
            releasable_scratch,
        );
        return;
    }
    b.free_vec(releasable_scratch);
    cconst_nbit_direct_trunc_fast_parked(b, acc, c, ctrl, window, park_low, true);
    b.reacquire_vec(releasable_scratch);
}

pub(crate) fn csub_nbit_const_direct_trunc_fast_releasing_scratch(
    b: &mut B,
    acc: &[QubitId],
    c: U256,
    ctrl: QubitId,
    window: usize,
    releasable_scratch: &[QubitId],
) {
    csub_nbit_const_direct_trunc_fast_releasing_scratch_at_step(
        b,
        acc,
        c,
        ctrl,
        window,
        releasable_scratch,
        None,
    );
}

pub(crate) fn csub_nbit_const_direct_trunc_fast_releasing_scratch_at_step(
    b: &mut B,
    acc: &[QubitId],
    c: U256,
    ctrl: QubitId,
    window: usize,
    releasable_scratch: &[QubitId],
    step: Option<usize>,
) {
    let park_low = special_fold_park_low_carries_at_step(step);
    if park_low == 0 || releasable_scratch.is_empty() {
        csub_nbit_const_direct_trunc_fast_borrowed_carries(
            b,
            acc,
            c,
            ctrl,
            window,
            releasable_scratch,
        );
        return;
    }
    b.free_vec(releasable_scratch);
    cconst_nbit_direct_trunc_fast_parked(b, acc, c, ctrl, window, park_low, false);
    b.reacquire_vec(releasable_scratch);
}

pub(crate) fn cadd_per_position_controls_trunc(
    b: &mut B,
    acc: &[QubitId],
    controls: &[Option<QubitId>],
    last: usize,
) {
    let n = acc.len();
    debug_assert!(last < n);
    debug_assert!(controls.len() <= n);
    let kctrl = |i: usize| -> Option<QubitId> {
        if i < controls.len() {
            controls[i]
        } else {
            None
        }
    };
    let maj2 = perpos_maj2_enabled();
    let carries = b.alloc_qubits(last + 1);

    for i in 0..=last {
        let target = carries[i];
        let carry_in = if i == 0 { None } else { Some(carries[i - 1]) };
        if let Some(kc) = kctrl(i) {
            if let Some(ci) = carry_in {
                emit_fold_majority(b, acc[i], kc, ci, target, maj2);
            } else {
                b.ccx(acc[i], kc, target);
            }
        } else if let Some(ci) = carry_in {
            b.ccx(acc[i], ci, target);
        }
    }

    for i in 0..n {
        if let Some(kc) = kctrl(i) {
            b.cx(kc, acc[i]);
        }
        if i > 0 && i - 1 <= last {
            b.cx(carries[i - 1], acc[i]);
        }
    }

    for i in (0..=last).rev() {
        let m = b.alloc_bit();
        b.hmr(carries[i], m);
        let carry_in = if i == 0 { None } else { Some(carries[i - 1]) };
        if let Some(kc) = kctrl(i) {
            b.x(acc[i]);
            if let Some(ci) = carry_in {
                b.cz_if(acc[i], kc, m);
                b.cz_if(acc[i], ci, m);
                b.x(acc[i]);
                b.cz_if(kc, ci, m);
            } else {
                b.cz_if(acc[i], kc, m);
                b.x(acc[i]);
            }
        } else if let Some(ci) = carry_in {
            b.x(acc[i]);
            b.cz_if(acc[i], ci, m);
            b.x(acc[i]);
        }
    }

    b.free_vec(&carries);
}

pub(crate) fn csub_per_position_controls_trunc(
    b: &mut B,
    acc: &[QubitId],
    controls: &[Option<QubitId>],
    last: usize,
) {
    let n = acc.len();
    debug_assert!(last < n);
    debug_assert!(controls.len() <= n);
    let kctrl = |i: usize| -> Option<QubitId> {
        if i < controls.len() {
            controls[i]
        } else {
            None
        }
    };
    let maj2 = perpos_maj2_enabled();
    let borrows = b.alloc_qubits(last + 1);

    for i in 0..=last {
        let target = borrows[i];
        let borrow_in = if i == 0 { None } else { Some(borrows[i - 1]) };
        if let Some(kc) = kctrl(i) {
            b.x(acc[i]);
            if let Some(bi) = borrow_in {
                emit_fold_majority(b, acc[i], kc, bi, target, maj2);
            } else {
                b.ccx(acc[i], kc, target);
            }
            b.x(acc[i]);
        } else if let Some(bi) = borrow_in {
            b.x(acc[i]);
            b.ccx(acc[i], bi, target);
            b.x(acc[i]);
        }
    }

    for i in 0..n {
        if let Some(kc) = kctrl(i) {
            b.cx(kc, acc[i]);
        }
        if i > 0 && i - 1 <= last {
            b.cx(borrows[i - 1], acc[i]);
        }
    }

    for i in (0..=last).rev() {
        let m = b.alloc_bit();
        b.hmr(borrows[i], m);
        let borrow_in = if i == 0 { None } else { Some(borrows[i - 1]) };
        if let Some(kc) = kctrl(i) {
            if let Some(bi) = borrow_in {
                b.cz_if(acc[i], kc, m);
                b.cz_if(acc[i], bi, m);
                b.cz_if(kc, bi, m);
            } else {
                b.cz_if(acc[i], kc, m);
            }
        } else if let Some(bi) = borrow_in {
            b.cz_if(acc[i], bi, m);
        }
    }

    b.free_vec(&borrows);
}

pub(crate) fn fold_freed_tail_enabled() -> bool {
    std::env::var("DIALOG_GCD_FOLD_FREED_TAIL").ok().as_deref() == Some("1")
}

pub(crate) fn fold_freed_tail_ed_enabled() -> bool {
    std::env::var("DIALOG_GCD_FOLD_FREED_TAIL_ED")
        .ok()
        .as_deref()
        == Some("1")
}

pub(crate) fn fold_host_derived_controls_enabled() -> bool {
    std::env::var("DIALOG_GCD_FOLD_HOST_DERIVED_CONTROLS")
        .ok()
        .as_deref()
        == Some("1")
}

fn fold_host_n10_enabled() -> bool {
    std::env::var("DIALOG_GCD_FOLD_HOST_N10").ok().as_deref() == Some("1")
}

fn fold_host_h_n10_enabled() -> bool {
    std::env::var("DIALOG_GCD_FOLD_HOST_H_N10").ok().as_deref() == Some("1")
}

fn fold_host_h_xed_n10_enabled() -> bool {
    std::env::var("DIALOG_GCD_FOLD_HOST_H_XED_N10")
        .ok()
        .as_deref()
        == Some("1")
}

fn fold_host_e_enabled() -> bool {
    std::env::var("DIALOG_GCD_FOLD_HOST_E").ok().as_deref() == Some("1")
}

fn fold_host_d_enabled() -> bool {
    std::env::var("DIALOG_GCD_FOLD_HOST_D").ok().as_deref() == Some("1")
}

pub(crate) fn fold_only_carry_trunc_window() -> Option<usize> {
    std::env::var("DIALOG_GCD_FOLD_CARRY_TRUNC_W")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .filter(|&w| w > 0)
}

fn fold_park_low_carries() -> usize {
    std::env::var("DIALOG_GCD_FOLD_PARK_LOW_CARRIES")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(0)
}

pub(crate) fn fold_park_low_carries_at_step(step: Option<usize>) -> usize {
    let mapped = step.and_then(|step| {
        let map = std::env::var("DIALOG_GCD_FOLD_PARK_LOW_CARRIES_STEP_MAP").ok()?;
        map.split(',').rev().find_map(|entry| {
            let (raw_step, raw_value) = entry.trim().split_once(':')?;
            if raw_step.trim().parse::<usize>().ok()? != step {
                return None;
            }
            raw_value.trim().parse::<usize>().ok()
        })
    });
    mapped.unwrap_or_else(fold_park_low_carries)
}

pub(crate) fn fold_stream_controls_enabled() -> bool {
    std::env::var("DIALOG_GCD_FOLD_STREAM_CONTROLS")
        .ok()
        .as_deref()
        == Some("1")
        && fold_park_low_carries() >= 12
}

fn fold_host_streamed_control_enabled() -> bool {
    std::env::var("DIALOG_GCD_FOLD_HOST_STREAMED_CONTROL")
        .ok()
        .as_deref()
        == Some("1")
        && fold_park_low_carries() >= 13
}

fn fold_host_e_top_carry_enabled() -> bool {
    std::env::var("DIALOG_GCD_FOLD_HOST_E_TOP_CARRY")
        .ok()
        .as_deref()
        == Some("1")
}

fn fold_host_d_carry12_enabled() -> bool {
    std::env::var("DIALOG_GCD_FOLD_HOST_D_CARRY12")
        .ok()
        .as_deref()
        == Some("1")
        && fold_park_low_carries() >= 14
}

fn fold_host_ovf2_carry13_enabled() -> bool {
    std::env::var("DIALOG_GCD_FOLD_HOST_OVF2_CARRY13")
        .ok()
        .as_deref()
        == Some("1")
        && fold_park_low_carries() >= 15
}

fn fold_free_first_high_carry_enabled() -> bool {
    std::env::var("DIALOG_GCD_FOLD_FREE_FIRST_HIGH_CARRY")
        .ok()
        .as_deref()
        == Some("1")
}

fn fold_stream_profile_phase(
    b: &mut B,
    add_phase: &'static str,
    sub_phase: &'static str,
    is_add: bool,
) {
    if std::env::var("DIALOG_GCD_FOLD_PROFILE_PHASES")
        .ok()
        .as_deref()
        == Some("1")
    {
        b.set_phase(if is_add { add_phase } else { sub_phase });
    }
}

fn fold_postsum_carry_phase_uncompute(
    b: &mut B,
    acc: &[QubitId],
    kctrl: Option<QubitId>,
    carry_in: Option<QubitId>,
    measured: BitId,
    i: usize,
    is_add: bool,
) {
    if is_add {
        if let Some(kc) = kctrl {
            b.x(acc[i]);
            if let Some(ci) = carry_in {
                b.cz_if(acc[i], kc, measured);
                b.cz_if(acc[i], ci, measured);
                b.x(acc[i]);
                b.cz_if(kc, ci, measured);
            } else {
                b.cz_if(acc[i], kc, measured);
                b.x(acc[i]);
            }
        } else if let Some(ci) = carry_in {
            b.x(acc[i]);
            b.cz_if(acc[i], ci, measured);
            b.x(acc[i]);
        }
    } else if let Some(kc) = kctrl {
        if let Some(ci) = carry_in {
            b.cz_if(acc[i], kc, measured);
            b.cz_if(acc[i], ci, measured);
            b.cz_if(kc, ci, measured);
        } else {
            b.cz_if(acc[i], kc, measured);
        }
    } else if let Some(ci) = carry_in {
        b.cz_if(acc[i], ci, measured);
    }
}

fn fold_postsum_carry_compute(
    b: &mut B,
    acc: &[QubitId],
    kctrl: Option<QubitId>,
    carry_in: Option<QubitId>,
    target: QubitId,
    i: usize,
    is_add: bool,
) {
    if is_add {
        if let Some(kc) = kctrl {
            b.x(acc[i]);
            if let Some(ci) = carry_in {
                emit_fold_majority(b, acc[i], kc, ci, target, perpos_maj2_enabled());
                b.x(acc[i]);
            } else {
                b.ccx(acc[i], kc, target);
                b.x(acc[i]);
            }
        } else if let Some(ci) = carry_in {
            b.x(acc[i]);
            b.ccx(acc[i], ci, target);
            b.x(acc[i]);
        }
    } else if let Some(kc) = kctrl {
        if let Some(ci) = carry_in {
            emit_fold_majority(b, acc[i], kc, ci, target, perpos_maj2_enabled());
        } else {
            b.ccx(acc[i], kc, target);
        }
    } else if let Some(ci) = carry_in {
        b.ccx(acc[i], ci, target);
    }
}

fn fold_presum_carry_compute_and_sum(
    b: &mut B,
    acc: &[QubitId],
    kctrl: Option<QubitId>,
    carry_in: Option<QubitId>,
    target: QubitId,
    i: usize,
    is_add: bool,
    maj2: bool,
) {
    if is_add {
        if let Some(kc) = kctrl {
            if let Some(ci) = carry_in {
                emit_fold_majority(b, acc[i], kc, ci, target, maj2);
            } else {
                b.ccx(acc[i], kc, target);
            }
        } else if let Some(ci) = carry_in {
            b.ccx(acc[i], ci, target);
        }
    } else if let Some(kc) = kctrl {
        b.x(acc[i]);
        if let Some(ci) = carry_in {
            emit_fold_majority(b, acc[i], kc, ci, target, maj2);
        } else {
            b.ccx(acc[i], kc, target);
        }
        b.x(acc[i]);
    } else if let Some(ci) = carry_in {
        b.x(acc[i]);
        b.ccx(acc[i], ci, target);
        b.x(acc[i]);
    }
    if let Some(kc) = kctrl {
        b.cx(kc, acc[i]);
    }
    if let Some(ci) = carry_in {
        b.cx(ci, acc[i]);
    }
}

pub(crate) fn secp_fold_controls(
    e: QubitId,
    d: QubitId,
    h: QubitId,
    xed: QubitId,
    eord: QubitId,
    n10: QubitId,
    hi_delta: usize,
    hi_c: usize,
) -> Vec<Option<QubitId>> {
    let mut controls: Vec<Option<QubitId>> = vec![None; hi_delta + 1];
    controls[0] = Some(e);
    controls[1] = Some(d);
    controls[4] = Some(e);
    controls[5] = Some(d);
    controls[6] = Some(e);
    controls[7] = Some(xed);
    controls[8] = Some(eord);
    controls[9] = Some(eord);
    controls[10] = Some(n10);
    controls[11] = Some(h);
    controls[hi_c] = Some(e);
    controls[hi_delta] = Some(d);
    controls
}

pub(crate) fn fold_ripple_freed_tail(
    b: &mut B,
    acc: &[QubitId],
    e: QubitId,
    d: QubitId,
    h: QubitId,
    xed: QubitId,
    eord: QubitId,
    n10: QubitId,
    last: usize,
    is_add: bool,
) {
    fold_ripple_freed_tail_ed(b, acc, e, d, h, xed, eord, n10, None, None, last, is_add);
}

pub(crate) fn fold_ripple_freed_tail_ed_streamed(
    b: &mut B,
    acc: &[QubitId],
    e: QubitId,
    d: QubitId,
    ed: Option<(QubitId, QubitId, QubitId)>,
    park_low: usize,
    last: usize,
    is_add: bool,
) {
    let free_ed = ed.is_some() && fold_freed_tail_ed_enabled();
    let n = acc.len();
    let hi_delta = 33usize;
    debug_assert!(last < n);
    debug_assert!(last > hi_delta, "freed-tail requires a nonempty high tail");
    let park_low = core::cmp::min(park_low, hi_delta);
    assert!(
        park_low >= 12,
        "streamed fold controls require at least 12 parked carries"
    );
    let host_streamed = fold_host_streamed_control_enabled();
    let host_e_top = free_ed && fold_host_e_top_carry_enabled();
    let host_d_carry12 = host_e_top && fold_host_d_carry12_enabled();
    let host_ovf2_carry13 = host_d_carry12 && fold_host_ovf2_carry13_enabled();
    let maj2 = perpos_maj2_enabled();
    let kctrl = |i: usize| match i {
        0 | 4 | 6 | 32 => Some(e),
        1 | 5 | 33 => Some(d),
        _ => None,
    };
    fold_stream_profile_phase(
        b,
        "dialog_gcd_streamed_double_active",
        "dialog_gcd_streamed_halve_active",
        is_add,
    );
    let low_chain_last = if host_e_top { hi_delta - 1 } else { hi_delta };
    let mut low = b.alloc_qubits(
        low_chain_last + 1 - usize::from(host_d_carry12) - usize::from(host_ovf2_carry13),
    );
    if host_d_carry12 {
        low.insert(12, d);
    }
    if host_ovf2_carry13 {
        let (_, ovf2, _) = ed.expect("host_ovf2_carry13 implies ed is Some");
        low.insert(13, ovf2);
        let (ovf1, _, _) = ed.expect("host_ovf2_carry13 implies ed is Some");
        b.cx(ovf1, ovf2);
        b.cx(d, ovf2);
        b.cx(e, ovf2);
    }
    let streamed_slot = low[park_low - 1];

    for i in 0..7 {
        let carry_in = if i == 0 { None } else { Some(low[i - 1]) };
        fold_presum_carry_compute_and_sum(b, acc, kctrl(i), carry_in, low[i], i, is_add, maj2);
    }
    let streamed = if host_streamed {
        streamed_slot
    } else {
        b.alloc_qubit()
    };
    b.cx(e, streamed);
    b.cx(d, streamed);
    fold_presum_carry_compute_and_sum(
        b,
        acc,
        Some(streamed),
        Some(low[6]),
        low[7],
        7,
        is_add,
        maj2,
    );
    b.ccx(e, d, streamed);
    for i in 8..10 {
        fold_presum_carry_compute_and_sum(
            b,
            acc,
            Some(streamed),
            Some(low[i - 1]),
            low[i],
            i,
            is_add,
            maj2,
        );
    }
    b.cx(e, streamed);
    fold_presum_carry_compute_and_sum(
        b,
        acc,
        Some(streamed),
        Some(low[9]),
        low[10],
        10,
        is_add,
        maj2,
    );
    b.cx(d, streamed);
    fold_presum_carry_compute_and_sum(
        b,
        acc,
        Some(streamed),
        Some(low[10]),
        low[11],
        11,
        is_add,
        maj2,
    );
    if host_streamed {
        b.ccx(e, d, streamed);
    }
    if host_d_carry12 {
        let (ovf1, _, s2) = ed.expect("host_d_carry12 implies ed is Some");
        b.ccx(ovf1, s2, d);
    }
    for i in 12..=low_chain_last {
        fold_presum_carry_compute_and_sum(
            b,
            acc,
            kctrl(i),
            Some(low[i - 1]),
            low[i],
            i,
            is_add,
            maj2,
        );
    }

    let free_first_high_carry = fold_free_first_high_carry_enabled()
        && park_low < low_chain_last
        && !(host_d_carry12 && park_low == 12)
        && !(host_ovf2_carry13 && park_low == 13);
    if free_first_high_carry {
        let m = b.alloc_bit();
        b.hmr(low[park_low], m);
        fold_postsum_carry_phase_uncompute(
            b,
            acc,
            kctrl(park_low),
            Some(low[park_low - 1]),
            m,
            park_low,
            is_add,
        );
        b.free(low[park_low]);
    }

    for i in (12..park_low).rev() {
        let m = b.alloc_bit();
        b.hmr(low[i], m);
        fold_postsum_carry_phase_uncompute(b, acc, kctrl(i), Some(low[i - 1]), m, i, is_add);
        b.free(low[i]);
    }
    if host_d_carry12 {
        let (ovf1, _, s2) = ed.expect("host_d_carry12 implies ed is Some");
        b.reacquire(d);
        b.ccx(ovf1, s2, d);
    }
    if host_streamed {
        b.reacquire(streamed);
        b.ccx(e, d, streamed);
    }
    let m11 = b.alloc_bit();
    b.hmr(low[11], m11);
    fold_postsum_carry_phase_uncompute(b, acc, Some(streamed), Some(low[10]), m11, 11, is_add);
    b.free(low[11]);
    b.cx(d, streamed);
    let m10 = b.alloc_bit();
    b.hmr(low[10], m10);
    fold_postsum_carry_phase_uncompute(b, acc, Some(streamed), Some(low[9]), m10, 10, is_add);
    b.free(low[10]);
    b.cx(e, streamed);
    for i in (8..10).rev() {
        let m = b.alloc_bit();
        b.hmr(low[i], m);
        fold_postsum_carry_phase_uncompute(b, acc, Some(streamed), Some(low[i - 1]), m, i, is_add);
        b.free(low[i]);
    }
    b.ccx(e, d, streamed);
    let m7 = b.alloc_bit();
    b.hmr(low[7], m7);
    fold_postsum_carry_phase_uncompute(b, acc, Some(streamed), Some(low[6]), m7, 7, is_add);
    b.free(low[7]);
    b.cx(d, streamed);
    b.cx(e, streamed);
    b.free(streamed);
    for i in (0..7).rev() {
        let m = b.alloc_bit();
        b.hmr(low[i], m);
        let carry_in = if i == 0 { None } else { Some(low[i - 1]) };
        fold_postsum_carry_phase_uncompute(b, acc, kctrl(i), carry_in, m, i, is_add);
        b.free(low[i]);
    }

    if host_ovf2_carry13 {
        let (ovf1, ovf2, _) = ed.expect("host_ovf2_carry13 implies ed is Some");
        b.reacquire(ovf2);
        b.cx(ovf1, ovf2);
        b.cx(d, ovf2);
        b.cx(e, ovf2);
    }

    if host_e_top {
        let (ovf1, ovf2, _) = ed.expect("host_e_top implies ed is Some");
        b.cx(ovf1, e);
        b.cx(d, e);
        b.cx(ovf2, e);
        fold_presum_carry_compute_and_sum(
            b,
            acc,
            Some(d),
            Some(low[hi_delta - 1]),
            e,
            hi_delta,
            is_add,
            maj2,
        );
    }

    if free_ed {
        let (ovf1, ovf2, s2) = ed.expect("free_ed implies ed is Some");
        if !host_e_top {
            b.cx(ovf1, e);
            b.cx(d, e);
            b.cx(ovf2, e);
            b.free(e);
        }
        let md = b.alloc_bit();
        b.hmr(d, md);
        b.cz_if(ovf1, s2, md);
        b.free(d);
    }

    fold_stream_profile_phase(
        b,
        "dialog_gcd_streamed_double_tail",
        "dialog_gcd_streamed_halve_tail",
        is_add,
    );
    let tail_len = last - hi_delta;
    let tail = b.alloc_qubits(tail_len);
    let cw = |i: usize| -> QubitId {
        if i < hi_delta {
            low[i]
        } else if i == hi_delta {
            if host_e_top {
                e
            } else {
                low[i]
            }
        } else {
            tail[i - hi_delta - 1]
        }
    };
    for i in hi_delta + 1..=last {
        if is_add {
            b.ccx(acc[i], cw(i - 1), cw(i));
        } else {
            b.x(acc[i]);
            b.ccx(acc[i], cw(i - 1), cw(i));
            b.x(acc[i]);
        }
    }
    for i in hi_delta + 1..n {
        if i - 1 <= last {
            b.cx(cw(i - 1), acc[i]);
        }
    }
    for i in (hi_delta + 1..=last).rev() {
        let m = b.alloc_bit();
        b.hmr(cw(i), m);
        let carry_in = cw(i - 1);
        if is_add {
            b.x(acc[i]);
            b.cz_if(acc[i], carry_in, m);
            b.x(acc[i]);
        } else {
            b.cz_if(acc[i], carry_in, m);
        }
        b.free(cw(i));
    }
    drop(tail);

    fold_stream_profile_phase(
        b,
        "dialog_gcd_streamed_double_reverse",
        "dialog_gcd_streamed_halve_reverse",
        is_add,
    );
    if free_ed {
        let (ovf1, ovf2, s2) = ed.expect("free_ed implies ed is Some");
        b.reacquire(d);
        b.ccx(ovf1, s2, d);
        if host_e_top {
            let m_top = b.alloc_bit();
            b.hmr(e, m_top);
            fold_postsum_carry_phase_uncompute(
                b,
                acc,
                Some(d),
                Some(low[hi_delta - 1]),
                m_top,
                hi_delta,
                is_add,
            );
            b.free(e);
        }
        b.reacquire(e);
        b.cx(ovf1, e);
        b.cx(d, e);
        b.cx(ovf2, e);
    }
    if host_ovf2_carry13 {
        let (ovf1, ovf2, _) = ed.expect("host_ovf2_carry13 implies ed is Some");
        b.cx(ovf1, ovf2);
        b.cx(d, ovf2);
        b.cx(e, ovf2);
    }

    for i in 0..park_low {
        if !(host_d_carry12 && i == 12) && !(host_ovf2_carry13 && i == 13) {
            b.reacquire(low[i]);
        }
    }
    for i in 0..7 {
        let carry_in = if i == 0 { None } else { Some(low[i - 1]) };
        fold_postsum_carry_compute(b, acc, kctrl(i), carry_in, low[i], i, is_add);
    }
    let streamed = if host_streamed {
        streamed_slot
    } else {
        b.alloc_qubit()
    };
    b.cx(e, streamed);
    b.cx(d, streamed);
    fold_postsum_carry_compute(b, acc, Some(streamed), Some(low[6]), low[7], 7, is_add);
    b.ccx(e, d, streamed);
    for i in 8..10 {
        fold_postsum_carry_compute(b, acc, Some(streamed), Some(low[i - 1]), low[i], i, is_add);
    }
    b.cx(e, streamed);
    fold_postsum_carry_compute(b, acc, Some(streamed), Some(low[9]), low[10], 10, is_add);
    b.cx(d, streamed);
    fold_postsum_carry_compute(b, acc, Some(streamed), Some(low[10]), low[11], 11, is_add);
    if host_streamed {
        b.ccx(e, d, streamed);
    }
    if host_d_carry12 {
        let (ovf1, _, s2) = ed.expect("host_d_carry12 implies ed is Some");
        b.ccx(ovf1, s2, d);
    }
    for i in 12..park_low {
        fold_postsum_carry_compute(b, acc, kctrl(i), Some(low[i - 1]), low[i], i, is_add);
    }
    if free_first_high_carry {
        b.reacquire(low[park_low]);
        fold_postsum_carry_compute(
            b,
            acc,
            kctrl(park_low),
            Some(low[park_low - 1]),
            low[park_low],
            park_low,
            is_add,
        );
    }

    for i in (12..=low_chain_last).rev() {
        let m = b.alloc_bit();
        b.hmr(low[i], m);
        fold_postsum_carry_phase_uncompute(b, acc, kctrl(i), Some(low[i - 1]), m, i, is_add);
        b.free(low[i]);
    }
    if host_d_carry12 {
        let (ovf1, _, s2) = ed.expect("host_d_carry12 implies ed is Some");
        b.reacquire(d);
        b.ccx(ovf1, s2, d);
    }
    if host_ovf2_carry13 {
        let (ovf1, ovf2, _) = ed.expect("host_ovf2_carry13 implies ed is Some");
        b.reacquire(ovf2);
        b.cx(ovf1, ovf2);
        b.cx(d, ovf2);
        b.cx(e, ovf2);
    }
    if host_streamed {
        b.reacquire(streamed);
        b.ccx(e, d, streamed);
    }
    let m11 = b.alloc_bit();
    b.hmr(low[11], m11);
    fold_postsum_carry_phase_uncompute(b, acc, Some(streamed), Some(low[10]), m11, 11, is_add);
    b.free(low[11]);
    b.cx(d, streamed);
    let m10 = b.alloc_bit();
    b.hmr(low[10], m10);
    fold_postsum_carry_phase_uncompute(b, acc, Some(streamed), Some(low[9]), m10, 10, is_add);
    b.free(low[10]);
    b.cx(e, streamed);
    for i in (8..10).rev() {
        let m = b.alloc_bit();
        b.hmr(low[i], m);
        fold_postsum_carry_phase_uncompute(b, acc, Some(streamed), Some(low[i - 1]), m, i, is_add);
        b.free(low[i]);
    }
    b.ccx(e, d, streamed);
    let m7 = b.alloc_bit();
    b.hmr(low[7], m7);
    fold_postsum_carry_phase_uncompute(b, acc, Some(streamed), Some(low[6]), m7, 7, is_add);
    b.free(low[7]);
    b.cx(d, streamed);
    b.cx(e, streamed);
    b.free(streamed);
    for i in (0..7).rev() {
        let m = b.alloc_bit();
        b.hmr(low[i], m);
        let carry_in = if i == 0 { None } else { Some(low[i - 1]) };
        fold_postsum_carry_phase_uncompute(b, acc, kctrl(i), carry_in, m, i, is_add);
        b.free(low[i]);
    }
    drop(low);
}

pub(crate) fn fold_ripple_freed_tail_ed(
    b: &mut B,
    acc: &[QubitId],
    e: QubitId,
    d: QubitId,
    h: QubitId,
    xed: QubitId,
    eord: QubitId,
    n10: QubitId,
    ed: Option<(QubitId, QubitId, QubitId)>,
    step: Option<usize>,
    last: usize,
    is_add: bool,
) {
    let configured_park_low = fold_park_low_carries_at_step(step);
    if fold_host_derived_controls_enabled() && configured_park_low <= 7 {
        fold_ripple_freed_tail_ed_hosted(
            b,
            acc,
            e,
            d,
            h,
            xed,
            eord,
            n10,
            ed,
            configured_park_low,
            last,
            is_add,
        );
        return;
    }
    if fold_stream_controls_enabled() && configured_park_low >= 12 {
        b.cx(h, n10);
        b.cx(d, n10);
        b.cx(h, eord);
        b.cx(xed, eord);
        b.cx(d, xed);
        b.cx(e, xed);
        b.free(n10);
        b.free(eord);
        b.free(xed);
        let mh = b.alloc_bit();
        b.hmr(h, mh);
        b.cz_if(e, d, mh);
        b.free(h);
        fold_ripple_freed_tail_ed_streamed(b, acc, e, d, ed, configured_park_low, last, is_add);
        b.reacquire(h);
        b.ccx(e, d, h);
        b.reacquire(xed);
        b.cx(e, xed);
        b.cx(d, xed);
        b.reacquire(eord);
        b.cx(xed, eord);
        b.cx(h, eord);
        b.reacquire(n10);
        b.cx(d, n10);
        b.cx(h, n10);
        return;
    }

    let free_ed = ed.is_some() && fold_freed_tail_ed_enabled();
    let n = acc.len();
    let hi_delta = 33usize;
    let hi_c = 32usize;
    debug_assert!(last < n);
    debug_assert!(last > hi_delta, "freed-tail requires a nonempty high tail");
    let controls = secp_fold_controls(e, d, h, xed, eord, n10, hi_delta, hi_c);
    let kctrl = |i: usize| controls.get(i).copied().flatten();
    let maj2 = perpos_maj2_enabled();
    let park_low = core::cmp::min(configured_park_low, hi_delta);
    let host_all_derived = fold_host_derived_controls_enabled() && park_low >= 15;
    let host_h_xed_n10 = (fold_host_h_xed_n10_enabled() || host_all_derived) && park_low >= 14;
    let host_h_n10 = (fold_host_h_n10_enabled() || host_h_xed_n10) && park_low >= 13;
    let host_xed = host_h_xed_n10;
    let host_eord = host_all_derived;
    let host_e = fold_host_e_enabled() && host_all_derived && free_ed && park_low >= 17;
    let host_d = fold_host_d_enabled() && host_e && park_low >= 18;
    let host_n10 = (fold_host_n10_enabled() || host_h_n10) && park_low >= 12;
    let stream_controls = fold_stream_controls_enabled() && park_low >= 12 && !host_n10;

    if stream_controls {
        b.cx(h, n10);
        b.cx(d, n10);
        b.cx(h, eord);
        b.cx(xed, eord);
        b.cx(d, xed);
        b.cx(e, xed);
        b.free(n10);
        b.free(eord);
        b.free(xed);
        let mh = b.alloc_bit();
        b.hmr(h, mh);
        b.cz_if(e, d, mh);
        b.free(h);
    }

    let low = if host_h_n10 {
        b.cx(h, n10);
        b.cx(d, n10);
        b.free(n10);
        if host_eord {
            b.cx(h, eord);
            b.cx(xed, eord);
            b.free(eord);
        }
        if host_xed {
            b.cx(d, xed);
            b.cx(e, xed);
            b.free(xed);
        }
        b.ccx(e, d, h);
        b.free(h);
        if host_e {
            let (ovf1, ovf2, _) = ed.expect("host_e requires overflow controls");
            b.cx(ovf1, e);
            b.cx(d, e);
            b.cx(ovf2, e);
            b.free(e);
        }

        let d_slot = host_d.then_some(d);
        let e_slot = if host_e {
            let slot = b.alloc_qubit();
            debug_assert_eq!(slot, e);
            Some(slot)
        } else {
            None
        };
        let h_slot = b.alloc_qubit();
        debug_assert_eq!(h_slot, h);
        let xed_slot = if host_xed {
            let slot = b.alloc_qubit();
            debug_assert_eq!(slot, xed);
            Some(slot)
        } else {
            None
        };
        let eord_slot = if host_eord {
            let slot = b.alloc_qubit();
            debug_assert_eq!(slot, eord);
            Some(slot)
        } else {
            None
        };
        let n10_slot = b.alloc_qubit();
        debug_assert_eq!(n10_slot, n10);
        let regular = b.alloc_qubits(
            hi_delta
                - 1
                - usize::from(host_xed)
                - usize::from(host_eord)
                - usize::from(host_e)
                - usize::from(host_d),
        );
        let mut low = Vec::with_capacity(hi_delta + 1);
        let mut next_regular = 0usize;
        for i in 0..=hi_delta {
            if i == 28 && host_d {
                low.push(n10_slot);
            } else if i == 29 && host_d {
                low.push(d_slot.expect("hosted d slot"));
            } else if i == 29 && host_e {
                low.push(n10_slot);
            } else if i == 30 {
                low.push(h_slot);
            } else if i == 31 && host_xed {
                low.push(xed_slot.expect("hosted xed slot"));
            } else if i == 32 && host_eord {
                low.push(eord_slot.expect("hosted eord slot"));
            } else if i == hi_delta {
                low.push(e_slot.unwrap_or(n10_slot));
            } else {
                low.push(regular[next_regular]);
                next_regular += 1;
            }
        }
        debug_assert_eq!(next_regular, regular.len());
        low
    } else if host_n10 {
        b.cx(h, n10);
        b.cx(d, n10);
        b.free(n10);
        let n10_slot = b.alloc_qubit();
        debug_assert_eq!(n10_slot, n10);
        let mut low = b.alloc_qubits(hi_delta);
        low.push(n10_slot);
        low
    } else {
        b.alloc_qubits(hi_delta + 1)
    };

    let mut tail_d = None;
    let mut streamed_forward = None;
    if host_h_n10 {
        let e_host = host_e.then_some(low[hi_delta]);
        let h_host = low[30];
        let xed_host = host_xed.then_some(low[31]);
        let eord_host = host_eord.then_some(low[32]);
        let d_host = host_d.then_some(low[29]);
        let n10_host = if host_d {
            low[28]
        } else if host_e {
            low[29]
        } else {
            low[hi_delta]
        };
        let e_ctrl = e_host.unwrap_or(e);
        let d_ctrl = d_host.unwrap_or(d);
        if let Some(e_host) = e_host {
            let (ovf1, ovf2, _) = ed.expect("host_e requires overflow controls");
            b.cx(ovf1, e_host);
            b.cx(d_ctrl, e_host);
            b.cx(ovf2, e_host);
        }
        b.ccx(e_ctrl, d_ctrl, h_host);
        if let Some(xed_host) = xed_host {
            b.cx(e_ctrl, xed_host);
            b.cx(d_ctrl, xed_host);
        }
        if let Some(eord_host) = eord_host {
            b.cx(xed_host.expect("hosted xed for eord"), eord_host);
            b.cx(h_host, eord_host);
        }
        b.cx(d_ctrl, n10_host);
        b.cx(h_host, n10_host);
        for i in 0..=hi_delta {
            if (host_d && i == 28) || (!host_d && host_e && i == 29) || (!host_e && i == 30) {
                b.cx(h_host, n10_host);
                b.cx(d_ctrl, n10_host);
                if let Some(eord_host) = eord_host {
                    b.cx(h_host, eord_host);
                    b.cx(xed_host.expect("hosted xed for eord"), eord_host);
                }
                if let Some(xed_host) = xed_host {
                    b.cx(d_ctrl, xed_host);
                    b.cx(e_ctrl, xed_host);
                }
                b.ccx(e_ctrl, d_ctrl, h_host);
            }
            if host_d && i == 29 {
                let (ovf1, _, s2) = ed.expect("host_d requires overflow controls");
                b.ccx(ovf1, s2, d_ctrl);
            }
            if host_e && i == hi_delta {
                let (ovf1, ovf2, _) = ed.expect("host_e requires overflow controls");
                b.cx(ovf2, e_ctrl);
                if host_d {
                    b.cx(tail_d.expect("tail d is live at bit 33"), e_ctrl);
                } else {
                    b.cx(d_ctrl, e_ctrl);
                }
                b.cx(ovf1, e_ctrl);
            }
            let target = low[i];
            let carry_in = if i == 0 { None } else { Some(low[i - 1]) };
            let kc = match i {
                0 | 4 | 6 | 32 if host_e => Some(e_ctrl),
                1 | 5 if host_d => d_host,
                7 if host_xed => xed_host,
                8 | 9 if host_eord => eord_host,
                10 => Some(n10_host),
                11 => Some(h_host),
                33 if host_d => tail_d,
                _ => kctrl(i),
            };
            if is_add {
                if let Some(kc) = kc {
                    if let Some(ci) = carry_in {
                        emit_fold_majority(b, acc[i], kc, ci, target, maj2);
                    } else {
                        b.ccx(acc[i], kc, target);
                    }
                } else if let Some(ci) = carry_in {
                    b.ccx(acc[i], ci, target);
                }
            } else if let Some(kc) = kc {
                b.x(acc[i]);
                if let Some(ci) = carry_in {
                    emit_fold_majority(b, acc[i], kc, ci, target, maj2);
                } else {
                    b.ccx(acc[i], kc, target);
                }
                b.x(acc[i]);
            } else if let Some(ci) = carry_in {
                b.x(acc[i]);
                b.ccx(acc[i], ci, target);
                b.x(acc[i]);
            }
            if let Some(kc) = kc {
                b.cx(kc, acc[i]);
            }
            if i > 0 {
                b.cx(low[i - 1], acc[i]);
            }
            if host_d && i == hi_delta - 1 {
                let measured = b.alloc_bit();
                b.hmr(low[31], measured);
                fold_postsum_carry_phase_uncompute(
                    b,
                    acc,
                    None,
                    Some(low[30]),
                    measured,
                    31,
                    is_add,
                );
                b.free(low[31]);
                let slot = b.alloc_qubit();
                debug_assert_eq!(slot, low[31]);
                let (ovf1, _, s2) = ed.expect("host_d requires overflow controls");
                b.ccx(ovf1, s2, slot);
                tail_d = Some(slot);
            }
        }
    } else if host_n10 {
        let n10_host = low[hi_delta];
        b.cx(d, n10_host);
        b.cx(h, n10_host);
        for i in 0..=hi_delta {
            if i == hi_delta {
                b.cx(h, n10_host);
                b.cx(d, n10_host);
            }
            let target = low[i];
            let carry_in = if i == 0 { None } else { Some(low[i - 1]) };
            let kc = if i == 10 { Some(n10_host) } else { kctrl(i) };
            if is_add {
                if let Some(kc) = kc {
                    if let Some(ci) = carry_in {
                        emit_fold_majority(b, acc[i], kc, ci, target, maj2);
                    } else {
                        b.ccx(acc[i], kc, target);
                    }
                } else if let Some(ci) = carry_in {
                    b.ccx(acc[i], ci, target);
                }
            } else if let Some(kc) = kc {
                b.x(acc[i]);
                if let Some(ci) = carry_in {
                    emit_fold_majority(b, acc[i], kc, ci, target, maj2);
                } else {
                    b.ccx(acc[i], kc, target);
                }
                b.x(acc[i]);
            } else if let Some(ci) = carry_in {
                b.x(acc[i]);
                b.ccx(acc[i], ci, target);
                b.x(acc[i]);
            }
            if let Some(kc) = kc {
                b.cx(kc, acc[i]);
            }
            if i > 0 {
                b.cx(low[i - 1], acc[i]);
            }
        }
    } else if stream_controls {
        for i in 0..7 {
            let carry_in = if i == 0 { None } else { Some(low[i - 1]) };
            fold_presum_carry_compute_and_sum(b, acc, kctrl(i), carry_in, low[i], i, is_add, maj2);
        }
        let streamed = b.alloc_qubit();
        b.cx(e, streamed);
        b.cx(d, streamed);
        fold_presum_carry_compute_and_sum(
            b,
            acc,
            Some(streamed),
            Some(low[6]),
            low[7],
            7,
            is_add,
            maj2,
        );
        b.ccx(e, d, streamed);
        for i in 8..10 {
            fold_presum_carry_compute_and_sum(
                b,
                acc,
                Some(streamed),
                Some(low[i - 1]),
                low[i],
                i,
                is_add,
                maj2,
            );
        }
        b.cx(e, streamed);
        fold_presum_carry_compute_and_sum(
            b,
            acc,
            Some(streamed),
            Some(low[9]),
            low[10],
            10,
            is_add,
            maj2,
        );
        b.cx(d, streamed);
        fold_presum_carry_compute_and_sum(
            b,
            acc,
            Some(streamed),
            Some(low[10]),
            low[11],
            11,
            is_add,
            maj2,
        );
        for i in 12..=hi_delta {
            fold_presum_carry_compute_and_sum(
                b,
                acc,
                kctrl(i),
                Some(low[i - 1]),
                low[i],
                i,
                is_add,
                maj2,
            );
        }
        streamed_forward = Some(streamed);
    } else {
        for i in 0..=hi_delta {
            let target = low[i];
            let carry_in = if i == 0 { None } else { Some(low[i - 1]) };
            if is_add {
                if let Some(kc) = kctrl(i) {
                    if let Some(ci) = carry_in {
                        emit_fold_majority(b, acc[i], kc, ci, target, maj2);
                    } else {
                        b.ccx(acc[i], kc, target);
                    }
                } else if let Some(ci) = carry_in {
                    b.ccx(acc[i], ci, target);
                }
            } else if let Some(kc) = kctrl(i) {
                b.x(acc[i]);
                if let Some(ci) = carry_in {
                    emit_fold_majority(b, acc[i], kc, ci, target, maj2);
                } else {
                    b.ccx(acc[i], kc, target);
                }
                b.x(acc[i]);
            } else if let Some(ci) = carry_in {
                b.x(acc[i]);
                b.ccx(acc[i], ci, target);
                b.x(acc[i]);
            }
        }

        for i in 0..=hi_delta {
            if let Some(kc) = kctrl(i) {
                b.cx(kc, acc[i]);
            }
            if i > 0 {
                b.cx(low[i - 1], acc[i]);
            }
        }
    }

    if park_low > 0 {
        if host_h_n10 {
            for i in (12..park_low).rev() {
                let m = b.alloc_bit();
                b.hmr(low[i], m);
                fold_postsum_carry_phase_uncompute(
                    b,
                    acc,
                    kctrl(i),
                    Some(low[i - 1]),
                    m,
                    i,
                    is_add,
                );
                b.free(low[i]);
            }

            let d_ctrl = if host_d {
                tail_d.expect("tail d is live during parked-carry cleanup")
            } else {
                d
            };
            let rev_e = if host_e {
                let (ovf1, ovf2, _) = ed.expect("host_e requires overflow controls");
                let rev_e = b.alloc_qubit();
                b.cx(ovf1, rev_e);
                b.cx(d_ctrl, rev_e);
                b.cx(ovf2, rev_e);
                Some(rev_e)
            } else {
                None
            };
            let e_ctrl = rev_e.unwrap_or(e);
            let rev_h = b.alloc_qubit();
            b.ccx(e_ctrl, d_ctrl, rev_h);
            let measured_h = b.alloc_bit();
            b.hmr(low[11], measured_h);
            fold_postsum_carry_phase_uncompute(
                b,
                acc,
                Some(rev_h),
                Some(low[10]),
                measured_h,
                11,
                is_add,
            );
            b.free(low[11]);

            let rev_n10 = b.alloc_qubit();
            debug_assert_eq!(rev_n10, low[11]);
            b.cx(d_ctrl, rev_n10);
            b.cx(rev_h, rev_n10);
            let measured_n10 = b.alloc_bit();
            b.hmr(low[10], measured_n10);
            fold_postsum_carry_phase_uncompute(
                b,
                acc,
                Some(rev_n10),
                Some(low[9]),
                measured_n10,
                10,
                is_add,
            );
            b.free(low[10]);
            b.cx(rev_h, rev_n10);
            b.cx(d_ctrl, rev_n10);
            b.free(rev_n10);

            let rev_xed = if host_xed {
                let rev_xed = b.alloc_qubit();
                b.cx(e_ctrl, rev_xed);
                b.cx(d_ctrl, rev_xed);
                Some(rev_xed)
            } else {
                None
            };
            let rev_eord = if host_eord {
                let rev_eord = b.alloc_qubit();
                b.cx(rev_xed.expect("hosted xed for eord"), rev_eord);
                b.cx(rev_h, rev_eord);
                Some(rev_eord)
            } else {
                None
            };
            if !host_xed {
                b.ccx(e, d, rev_h);
                b.free(rev_h);
            }
            for i in (0..10).rev() {
                let m = b.alloc_bit();
                b.hmr(low[i], m);
                let carry_in = if i == 0 { None } else { Some(low[i - 1]) };
                let kc = match i {
                    0 | 4 | 6 if host_e => rev_e,
                    1 | 5 if host_d => Some(d_ctrl),
                    7 if host_xed => rev_xed,
                    8 | 9 if host_eord => rev_eord,
                    _ => kctrl(i),
                };
                fold_postsum_carry_phase_uncompute(b, acc, kc, carry_in, m, i, is_add);
                b.free(low[i]);
            }
            if let Some(rev_eord) = rev_eord {
                b.cx(rev_h, rev_eord);
                b.cx(rev_xed.expect("hosted xed for eord"), rev_eord);
                b.free(rev_eord);
            }
            if let Some(rev_xed) = rev_xed {
                b.cx(d_ctrl, rev_xed);
                b.cx(e_ctrl, rev_xed);
                b.free(rev_xed);
            }
            if host_xed {
                b.ccx(e_ctrl, d_ctrl, rev_h);
                b.free(rev_h);
            }
            if let Some(rev_e) = rev_e {
                let (ovf1, ovf2, _) = ed.expect("host_e requires overflow controls");
                b.cx(ovf2, rev_e);
                b.cx(d_ctrl, rev_e);
                b.cx(ovf1, rev_e);
                b.free(rev_e);
            }
        } else if host_n10 {
            for i in (11..park_low).rev() {
                let m = b.alloc_bit();
                b.hmr(low[i], m);
                fold_postsum_carry_phase_uncompute(
                    b,
                    acc,
                    kctrl(i),
                    Some(low[i - 1]),
                    m,
                    i,
                    is_add,
                );
                b.free(low[i]);
            }
            let rev_n10 = b.alloc_qubit();
            b.cx(d, rev_n10);
            b.cx(h, rev_n10);
            for i in (0..11).rev() {
                let m = b.alloc_bit();
                b.hmr(low[i], m);
                let carry_in = if i == 0 { None } else { Some(low[i - 1]) };
                let kc = if i == 10 { Some(rev_n10) } else { kctrl(i) };
                fold_postsum_carry_phase_uncompute(b, acc, kc, carry_in, m, i, is_add);
                b.free(low[i]);
            }
            b.cx(h, rev_n10);
            b.cx(d, rev_n10);
            b.free(rev_n10);
        } else if stream_controls {
            let streamed = streamed_forward
                .take()
                .expect("streamed forward control is live");
            for i in (12..park_low).rev() {
                let m = b.alloc_bit();
                b.hmr(low[i], m);
                fold_postsum_carry_phase_uncompute(
                    b,
                    acc,
                    kctrl(i),
                    Some(low[i - 1]),
                    m,
                    i,
                    is_add,
                );
                b.free(low[i]);
            }
            let m11 = b.alloc_bit();
            b.hmr(low[11], m11);
            fold_postsum_carry_phase_uncompute(
                b,
                acc,
                Some(streamed),
                Some(low[10]),
                m11,
                11,
                is_add,
            );
            b.free(low[11]);
            b.cx(d, streamed);
            let m10 = b.alloc_bit();
            b.hmr(low[10], m10);
            fold_postsum_carry_phase_uncompute(
                b,
                acc,
                Some(streamed),
                Some(low[9]),
                m10,
                10,
                is_add,
            );
            b.free(low[10]);
            b.cx(e, streamed);
            for i in (8..10).rev() {
                let m = b.alloc_bit();
                b.hmr(low[i], m);
                fold_postsum_carry_phase_uncompute(
                    b,
                    acc,
                    Some(streamed),
                    Some(low[i - 1]),
                    m,
                    i,
                    is_add,
                );
                b.free(low[i]);
            }
            b.ccx(e, d, streamed);
            let m7 = b.alloc_bit();
            b.hmr(low[7], m7);
            fold_postsum_carry_phase_uncompute(b, acc, Some(streamed), Some(low[6]), m7, 7, is_add);
            b.free(low[7]);
            b.cx(d, streamed);
            b.cx(e, streamed);
            b.free(streamed);
            for i in (0..7).rev() {
                let m = b.alloc_bit();
                b.hmr(low[i], m);
                let carry_in = if i == 0 { None } else { Some(low[i - 1]) };
                fold_postsum_carry_phase_uncompute(b, acc, kctrl(i), carry_in, m, i, is_add);
                b.free(low[i]);
            }
        } else {
            for i in (0..park_low).rev() {
                let m = b.alloc_bit();
                b.hmr(low[i], m);
                let carry_in = if i == 0 { None } else { Some(low[i - 1]) };
                fold_postsum_carry_phase_uncompute(b, acc, kctrl(i), carry_in, m, i, is_add);
                b.free(low[i]);
            }
        }
    }
    debug_assert!(streamed_forward.is_none());

    if !stream_controls {
        if !host_n10 {
            b.cx(h, n10);
            b.cx(d, n10);
            b.free(n10);
        }
        if !host_eord {
            if host_h_n10 {
                let rev_h = b.alloc_qubit();
                b.ccx(e, d, rev_h);
                b.cx(rev_h, eord);
                b.ccx(e, d, rev_h);
                b.free(rev_h);
            } else {
                b.cx(h, eord);
            }
            if host_xed {
                let rev_xed = b.alloc_qubit();
                b.cx(e, rev_xed);
                b.cx(d, rev_xed);
                b.cx(rev_xed, eord);
                b.cx(d, rev_xed);
                b.cx(e, rev_xed);
                b.free(rev_xed);
            } else {
                b.cx(xed, eord);
            }
            b.free(eord);
        }
        if !host_xed {
            b.cx(d, xed);
            b.cx(e, xed);
            b.free(xed);
        }
        if !host_h_n10 {
            let mh = b.alloc_bit();
            b.hmr(h, mh);
            b.cz_if(e, d, mh);
            b.free(h);
        }
    }

    if free_ed {
        let (ovf1, ovf2, s2) = ed.expect("free_ed implies ed is Some");
        if !host_e {
            b.cx(ovf1, e);
            b.cx(d, e);
            b.cx(ovf2, e);
            b.free(e);
        }
        if !host_d {
            let md = b.alloc_bit();
            b.hmr(d, md);
            b.cz_if(ovf1, s2, md);
            b.free(d);
        }
    }

    let tail_len = last - hi_delta;
    let tail = b.alloc_qubits(tail_len);
    let cw = |i: usize| -> QubitId {
        if i <= hi_delta {
            low[i]
        } else {
            tail[i - hi_delta - 1]
        }
    };

    for i in (hi_delta + 1)..=last {
        if is_add {
            b.ccx(acc[i], cw(i - 1), cw(i));
        } else {
            b.x(acc[i]);
            b.ccx(acc[i], cw(i - 1), cw(i));
            b.x(acc[i]);
        }
    }

    for i in (hi_delta + 1)..n {
        if i - 1 <= last {
            b.cx(cw(i - 1), acc[i]);
        }
    }

    for i in (hi_delta + 1..=last).rev() {
        let m = b.alloc_bit();
        b.hmr(cw(i), m);
        let carry_in = cw(i - 1);
        if is_add {
            b.x(acc[i]);
            b.cz_if(acc[i], carry_in, m);
            b.x(acc[i]);
        } else {
            b.cz_if(acc[i], carry_in, m);
        }
        b.free(cw(i));
    }
    drop(tail);

    if free_ed && !host_d {
        let (ovf1, ovf2, s2) = ed.expect("free_ed implies ed is Some");
        b.reacquire(d);
        b.ccx(ovf1, s2, d);
        if !host_e {
            b.reacquire(e);
            b.cx(ovf1, e);
            b.cx(d, e);
            b.cx(ovf2, e);
        }
    }
    if host_h_n10 {
        if host_e {
            let measured = b.alloc_bit();
            b.hmr(low[hi_delta], measured);
            if host_d {
                fold_postsum_carry_phase_uncompute(
                    b,
                    acc,
                    tail_d,
                    Some(low[hi_delta - 1]),
                    measured,
                    hi_delta,
                    is_add,
                );
            } else {
                fold_postsum_carry_phase_uncompute(
                    b,
                    acc,
                    Some(d),
                    Some(low[hi_delta - 1]),
                    measured,
                    hi_delta,
                    is_add,
                );
            }
            b.free(low[hi_delta]);
            let (ovf1, ovf2, _) = ed.expect("host_e requires overflow controls");
            b.reacquire(e);
            b.cx(ovf1, e);
            if host_d {
                b.cx(tail_d.expect("tail d is live while restoring e"), e);
            } else {
                b.cx(d, e);
            }
            b.cx(ovf2, e);
        }
        if host_d {
            let carry31 = b.alloc_qubit();
            fold_postsum_carry_compute(b, acc, None, Some(low[30]), carry31, 31, is_add);
            let measured32 = b.alloc_bit();
            b.hmr(low[32], measured32);
            fold_postsum_carry_phase_uncompute(
                b,
                acc,
                Some(e),
                Some(carry31),
                measured32,
                32,
                is_add,
            );
            b.free(low[32]);
            let measured31 = b.alloc_bit();
            b.hmr(carry31, measured31);
            fold_postsum_carry_phase_uncompute(b, acc, None, Some(low[30]), measured31, 31, is_add);
            b.free(carry31);

            for i in (28..=30).rev() {
                let measured = b.alloc_bit();
                b.hmr(low[i], measured);
                fold_postsum_carry_phase_uncompute(
                    b,
                    acc,
                    None,
                    Some(low[i - 1]),
                    measured,
                    i,
                    is_add,
                );
                b.free(low[i]);
            }

            let d_tail = tail_d.expect("tail d is live during d restoration");
            b.reacquire(d);
            b.cx(d_tail, d);
            b.cx(d, d_tail);
            b.free(d_tail);
        } else {
            let high_start = if host_e { 29 } else { 30 };
            let high_end = if host_e { hi_delta - 1 } else { hi_delta };
            for i in (high_start..=high_end).rev() {
                let measured = b.alloc_bit();
                b.hmr(low[i], measured);
                let kc = match i {
                    32 => Some(e),
                    33 => Some(d),
                    _ => None,
                };
                fold_postsum_carry_phase_uncompute(
                    b,
                    acc,
                    kc,
                    Some(low[i - 1]),
                    measured,
                    i,
                    is_add,
                );
                b.free(low[i]);
            }
        }
    } else if host_n10 {
        let measured = b.alloc_bit();
        b.hmr(low[hi_delta], measured);
        fold_postsum_carry_phase_uncompute(
            b,
            acc,
            Some(d),
            Some(low[hi_delta - 1]),
            measured,
            hi_delta,
            is_add,
        );
        b.free(low[hi_delta]);
    }
    if stream_controls {
        for i in 0..park_low {
            b.reacquire(low[i]);
        }

        for i in 0..7 {
            let carry_in = if i == 0 { None } else { Some(low[i - 1]) };
            fold_postsum_carry_compute(b, acc, kctrl(i), carry_in, low[i], i, is_add);
        }

        let streamed = b.alloc_qubit();
        b.cx(e, streamed);
        b.cx(d, streamed);
        fold_postsum_carry_compute(b, acc, Some(streamed), Some(low[6]), low[7], 7, is_add);
        b.ccx(e, d, streamed);
        for i in 8..10 {
            fold_postsum_carry_compute(b, acc, Some(streamed), Some(low[i - 1]), low[i], i, is_add);
        }
        b.cx(e, streamed);
        fold_postsum_carry_compute(b, acc, Some(streamed), Some(low[9]), low[10], 10, is_add);
        b.cx(d, streamed);
        fold_postsum_carry_compute(b, acc, Some(streamed), Some(low[10]), low[11], 11, is_add);
        for i in 12..park_low {
            fold_postsum_carry_compute(b, acc, kctrl(i), Some(low[i - 1]), low[i], i, is_add);
        }

        for i in (12..=hi_delta).rev() {
            let m = b.alloc_bit();
            b.hmr(low[i], m);
            fold_postsum_carry_phase_uncompute(b, acc, kctrl(i), Some(low[i - 1]), m, i, is_add);
            b.free(low[i]);
        }

        let m11 = b.alloc_bit();
        b.hmr(low[11], m11);
        fold_postsum_carry_phase_uncompute(b, acc, Some(streamed), Some(low[10]), m11, 11, is_add);
        b.free(low[11]);
        b.cx(d, streamed);

        let m10 = b.alloc_bit();
        b.hmr(low[10], m10);
        fold_postsum_carry_phase_uncompute(b, acc, Some(streamed), Some(low[9]), m10, 10, is_add);
        b.free(low[10]);
        b.cx(e, streamed);

        for i in (8..10).rev() {
            let m = b.alloc_bit();
            b.hmr(low[i], m);
            fold_postsum_carry_phase_uncompute(
                b,
                acc,
                Some(streamed),
                Some(low[i - 1]),
                m,
                i,
                is_add,
            );
            b.free(low[i]);
        }
        b.ccx(e, d, streamed);

        let m7 = b.alloc_bit();
        b.hmr(low[7], m7);
        fold_postsum_carry_phase_uncompute(b, acc, Some(streamed), Some(low[6]), m7, 7, is_add);
        b.free(low[7]);
        b.cx(d, streamed);
        b.cx(e, streamed);
        b.free(streamed);

        for i in (0..7).rev() {
            let m = b.alloc_bit();
            b.hmr(low[i], m);
            let carry_in = if i == 0 { None } else { Some(low[i - 1]) };
            fold_postsum_carry_phase_uncompute(b, acc, kctrl(i), carry_in, m, i, is_add);
            b.free(low[i]);
        }
        drop(low);

        b.reacquire(h);
        b.ccx(e, d, h);
        b.reacquire(xed);
        b.cx(e, xed);
        b.cx(d, xed);
        b.reacquire(eord);
        b.cx(xed, eord);
        b.cx(h, eord);
        b.reacquire(n10);
        b.cx(d, n10);
        b.cx(h, n10);
    } else {
        b.reacquire(h);
        b.ccx(e, d, h);
        b.reacquire(xed);
        b.cx(e, xed);
        b.cx(d, xed);
        b.reacquire(eord);
        b.cx(xed, eord);
        b.cx(h, eord);
        b.reacquire(n10);
        b.cx(d, n10);
        b.cx(h, n10);

        if park_low > 0 {
            for i in 0..park_low {
                b.reacquire(low[i]);
                let carry_in = if i == 0 { None } else { Some(low[i - 1]) };
                fold_postsum_carry_compute(b, acc, kctrl(i), carry_in, low[i], i, is_add);
            }
        }

        let low_top = if host_d {
            27
        } else if host_e {
            28
        } else if host_h_n10 {
            29
        } else if host_n10 {
            hi_delta - 1
        } else {
            hi_delta
        };
        for i in (0..=low_top).rev() {
            let m = b.alloc_bit();
            b.hmr(low[i], m);
            let carry_in = if i == 0 { None } else { Some(low[i - 1]) };
            fold_postsum_carry_phase_uncompute(b, acc, kctrl(i), carry_in, m, i, is_add);
            b.free(low[i]);
        }
        drop(low);
    }
}

fn fold_ripple_freed_tail_ed_hosted(
    b: &mut B,
    acc: &[QubitId],
    e: QubitId,
    d: QubitId,
    h: QubitId,
    xed: QubitId,
    eord: QubitId,
    n10: QubitId,
    ed: Option<(QubitId, QubitId, QubitId)>,
    park_low: usize,
    last: usize,
    is_add: bool,
) {
    let free_ed = ed.is_some() && fold_freed_tail_ed_enabled();
    let n = acc.len();
    let hi_delta = 33usize;
    let split = 11usize;
    debug_assert!(last < n);
    debug_assert!(last > hi_delta, "freed-tail requires a nonempty high tail");
    let maj2 = perpos_maj2_enabled();
    let park_low = core::cmp::min(park_low, split + 1);

    b.cx(h, n10);
    b.cx(d, n10);
    b.free(n10);
    b.cx(h, eord);
    b.cx(xed, eord);
    b.cx(d, xed);
    b.cx(e, xed);
    b.free(eord);
    b.free(xed);
    let mh = b.alloc_bit();
    b.hmr(h, mh);
    b.cz_if(e, d, mh);
    b.free(h);

    let low = b.alloc_qubits(hi_delta + 1);

    let h_host = low[30];
    let xed_host = low[31];
    let eord_host = low[32];
    let n10_host = low[33];
    b.ccx(e, d, h_host);
    b.cx(e, xed_host);
    b.cx(d, xed_host);
    b.cx(xed_host, eord_host);
    b.cx(h_host, eord_host);
    b.cx(d, n10_host);
    b.cx(h_host, n10_host);

    let hosted_kctrl = |i: usize| -> Option<QubitId> {
        match i {
            0 | 4 | 6 | 32 => Some(e),
            1 | 5 | 33 => Some(d),
            7 => Some(xed_host),
            8 | 9 => Some(eord_host),
            10 => Some(n10_host),
            11 => Some(h_host),
            _ => None,
        }
    };

    for i in 0..=split {
        let target = low[i];
        let carry_in = if i == 0 { None } else { Some(low[i - 1]) };
        if is_add {
            if let Some(kc) = hosted_kctrl(i) {
                if let Some(ci) = carry_in {
                    emit_fold_majority(b, acc[i], kc, ci, target, maj2);
                } else {
                    b.ccx(acc[i], kc, target);
                }
            } else if let Some(ci) = carry_in {
                b.ccx(acc[i], ci, target);
            }
        } else if let Some(kc) = hosted_kctrl(i) {
            b.x(acc[i]);
            if let Some(ci) = carry_in {
                emit_fold_majority(b, acc[i], kc, ci, target, maj2);
            } else {
                b.ccx(acc[i], kc, target);
            }
            b.x(acc[i]);
        } else if let Some(ci) = carry_in {
            b.x(acc[i]);
            b.ccx(acc[i], ci, target);
            b.x(acc[i]);
        }
    }
    for i in 0..=split {
        if let Some(kc) = hosted_kctrl(i) {
            b.cx(kc, acc[i]);
        }
        if i > 0 {
            b.cx(low[i - 1], acc[i]);
        }
    }

    b.cx(h_host, n10_host);
    b.cx(d, n10_host);
    b.cx(h_host, eord_host);
    b.cx(xed_host, eord_host);
    b.cx(d, xed_host);
    b.cx(e, xed_host);
    let mh_host = b.alloc_bit();
    b.hmr(h_host, mh_host);
    b.cz_if(e, d, mh_host);

    for i in split + 1..=hi_delta {
        let target = low[i];
        let carry_in = Some(low[i - 1]);
        let kctrl = match i {
            32 => Some(e),
            33 => Some(d),
            _ => None,
        };
        if is_add {
            if let Some(kc) = kctrl {
                emit_fold_majority(
                    b,
                    acc[i],
                    kc,
                    carry_in.expect("high carry-in"),
                    target,
                    maj2,
                );
            } else {
                b.ccx(acc[i], carry_in.expect("high carry-in"), target);
            }
        } else if let Some(kc) = kctrl {
            b.x(acc[i]);
            emit_fold_majority(
                b,
                acc[i],
                kc,
                carry_in.expect("high carry-in"),
                target,
                maj2,
            );
            b.x(acc[i]);
        } else {
            b.x(acc[i]);
            b.ccx(acc[i], carry_in.expect("high carry-in"), target);
            b.x(acc[i]);
        }
    }
    for i in split + 1..=hi_delta {
        let kctrl = match i {
            32 => Some(e),
            33 => Some(d),
            _ => None,
        };
        if let Some(kc) = kctrl {
            b.cx(kc, acc[i]);
        }
        b.cx(low[i - 1], acc[i]);
    }

    if park_low > 0 {
        let controls =
            secp_fold_controls(e, d, h_host, xed_host, eord_host, n10_host, hi_delta, 32);
        for i in (0..park_low).rev() {
            let m = b.alloc_bit();
            b.hmr(low[i], m);
            let carry_in = if i == 0 { None } else { Some(low[i - 1]) };
            fold_postsum_carry_phase_uncompute(
                b,
                acc,
                controls.get(i).copied().flatten(),
                carry_in,
                m,
                i,
                is_add,
            );
            b.free(low[i]);
        }
    }

    if free_ed {
        let (ovf1, ovf2, s2) = ed.expect("free_ed implies ed is Some");
        b.cx(ovf1, e);
        b.cx(d, e);
        b.cx(ovf2, e);
        b.free(e);
        let md = b.alloc_bit();
        b.hmr(d, md);
        b.cz_if(ovf1, s2, md);
        b.free(d);
    }

    let tail_len = last - hi_delta;
    let tail = b.alloc_qubits(tail_len);
    let cw = |i: usize| -> QubitId {
        if i <= hi_delta {
            low[i]
        } else {
            tail[i - hi_delta - 1]
        }
    };

    for i in (hi_delta + 1)..=last {
        if is_add {
            b.ccx(acc[i], cw(i - 1), cw(i));
        } else {
            b.x(acc[i]);
            b.ccx(acc[i], cw(i - 1), cw(i));
            b.x(acc[i]);
        }
    }
    for i in (hi_delta + 1)..n {
        if i - 1 <= last {
            b.cx(cw(i - 1), acc[i]);
        }
    }
    for i in (hi_delta + 1..=last).rev() {
        let m = b.alloc_bit();
        b.hmr(cw(i), m);
        let carry_in = cw(i - 1);
        if is_add {
            b.x(acc[i]);
            b.cz_if(acc[i], carry_in, m);
            b.x(acc[i]);
        } else {
            b.cz_if(acc[i], carry_in, m);
        }
        b.free(cw(i));
    }
    drop(tail);

    if free_ed {
        let (ovf1, ovf2, s2) = ed.expect("free_ed implies ed is Some");
        b.reacquire(d);
        b.ccx(ovf1, s2, d);
        b.reacquire(e);
        b.cx(ovf1, e);
        b.cx(d, e);
        b.cx(ovf2, e);
    }

    for i in (split + 1..=hi_delta).rev() {
        let m = b.alloc_bit();
        b.hmr(low[i], m);
        let kctrl = match i {
            32 => Some(e),
            33 => Some(d),
            _ => None,
        };
        fold_postsum_carry_phase_uncompute(b, acc, kctrl, Some(low[i - 1]), m, i, is_add);
        b.free(low[i]);
    }

    let rev_h = b.alloc_qubit();
    let rev_xed = b.alloc_qubit();
    let rev_eord = b.alloc_qubit();
    let rev_n10 = b.alloc_qubit();
    b.ccx(e, d, rev_h);
    b.cx(e, rev_xed);
    b.cx(d, rev_xed);
    b.cx(rev_xed, rev_eord);
    b.cx(rev_h, rev_eord);
    b.cx(d, rev_n10);
    b.cx(rev_h, rev_n10);
    let controls = secp_fold_controls(e, d, rev_h, rev_xed, rev_eord, rev_n10, hi_delta, 32);

    if park_low > 0 {
        for i in 0..park_low {
            b.reacquire(low[i]);
            let carry_in = if i == 0 { None } else { Some(low[i - 1]) };
            fold_postsum_carry_compute(
                b,
                acc,
                controls.get(i).copied().flatten(),
                carry_in,
                low[i],
                i,
                is_add,
            );
        }
    }

    for i in (0..=split).rev() {
        let m = b.alloc_bit();
        b.hmr(low[i], m);
        let carry_in = if i == 0 { None } else { Some(low[i - 1]) };
        fold_postsum_carry_phase_uncompute(
            b,
            acc,
            controls.get(i).copied().flatten(),
            carry_in,
            m,
            i,
            is_add,
        );
        b.free(low[i]);
    }
    drop(low);

    b.reacquire(h);
    b.cx(rev_h, h);
    b.reacquire(xed);
    b.cx(rev_xed, xed);
    b.reacquire(eord);
    b.cx(rev_eord, eord);
    b.reacquire(n10);
    b.cx(rev_n10, n10);

    b.cx(rev_h, rev_n10);
    b.cx(d, rev_n10);
    b.cx(rev_h, rev_eord);
    b.cx(rev_xed, rev_eord);
    b.cx(d, rev_xed);
    b.cx(e, rev_xed);
    b.free(rev_n10);
    b.free(rev_eord);
    b.free(rev_xed);
    let mh_rev = b.alloc_bit();
    b.hmr(rev_h, mh_rev);
    b.cz_if(e, d, mh_rev);
    b.free(rev_h);
}
