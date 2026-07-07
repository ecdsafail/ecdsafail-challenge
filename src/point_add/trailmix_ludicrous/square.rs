
use super::arith::{self, cuccaro_carry, mod_add_lowpeak, mod_add_shifted_low, mod_sub, mod_sub_shifted_low, F_SECP256K1, LSBS};
use super::{B, BExt};
use crate::circuit::{QubitId};

const N: usize = 256;

fn clear_and(circ: &mut B, t: &QubitId, a: &QubitId, b: &QubitId) {
    let bit = circ.alloc_bit();
    circ.hmr(*t, bit);
    circ.cz_if_bit(*a, *b, bit);
}

const F_NAF_TERMS: [(usize, ShiftOp); 5] = [
    (0, ShiftOp::Sub),
    (4, ShiftOp::Sub),
    (6, ShiftOp::Add),
    (10, ShiftOp::Sub),
    (32, ShiftOp::Sub),
];

#[derive(Copy, Clone)]
enum ShiftOp {
    Add,
    Sub,
}

fn add_f_window_shifted(circ: &mut B, ctrl: &QubitId, reg: &[QubitId], offset: usize) {
    let f_bytes = F_SECP256K1.to_le_bytes();
    arith::add_f_window_pub(circ, ctrl, &reg[offset..], LSBS, &f_bytes, None);
}

fn sub_f_window_shifted(circ: &mut B, ctrl: &QubitId, reg: &[QubitId], offset: usize) {
    for q in &reg[offset..offset + LSBS] {
        circ.x(*q);
    }
    add_f_window_shifted(circ, ctrl, reg, offset);
    for q in &reg[offset..offset + LSBS] {
        circ.x(*q);
    }
}

fn apply_shifted_hi_term(
    circ: &mut B,
    hi: &[QubitId],
    output_reg: &[QubitId],
    shift: usize,
    op: ShiftOp,
) {
    let n = hi.len();
    assert_eq!(n, 256, "hi must be 256 bits");
    assert!(shift < n, "shift must be less than 256");

    match op {
        ShiftOp::Add => mod_add_shifted_low(circ, &hi[..n - shift], output_reg, shift),
        ShiftOp::Sub => {
            if shift == 0 {
                mod_sub(circ, hi, output_reg);
            } else {
                mod_sub_shifted_low(circ, &hi[..n - shift], output_reg, shift);
            }
        }
    }

    for t in 0..shift {
        let ctrl = &hi[n - shift + t];
        match op {
            ShiftOp::Add => add_f_window_shifted(circ, ctrl, output_reg, t),
            ShiftOp::Sub => sub_f_window_shifted(circ, ctrl, output_reg, t),
        }
    }
}

fn add_into(circ: &mut B, slice: &[QubitId], row: &[QubitId]) {
    let m = row.len();
    assert_eq!(slice.len(), m + 1, "slice must be one wider than row");
    if m == 0 {
        return;
    }

    let pad = circ.alloc_qubit();
    let mut b: Vec<QubitId> = row.to_vec();
    b.push(pad);
    let k = super::next_sqrow_k();
    super::arith::hybrid_add_adaptive(circ, slice, &b, k);
    circ.zero_and_free(pad);
}

fn symmetric_square_into_prod(circ: &mut B, x: &[QubitId], prod: &mut Vec<QubitId>) {
    let n = x.len();
    if std::env::var("TLM_SQ_TRACE").ok().as_deref() == Some("1") {
        eprintln!("SQ_CALL fwd n={n} crosses={}", n * (n - 1) / 2);
    }
    assert!(prod.is_empty(), "prod is grown lazily; pass an empty Vec");

    if square_addsub_enabled() && !(square_addsub_skip_c() && n == 129) {
        for _ in 0..(2 * n) {
            prod.push(circ.alloc_qubit());
        }
        if square_addsub_local_diag() {
            crate::point_add::arith::square_addsub_local(circ, x, prod);
        } else {
            crate::point_add::arith::square_addsub_vented(circ, x, prod);
        }
        return;
    }
    for i in 0..n {

        let num_cross = n.saturating_sub(i + 1);
        let width = if i == n - 1 { 1 } else { n - i + 1 };

        let hi = (2 * i + width + 1).min(2 * n);
        while prod.len() < hi {
            prod.push(circ.alloc_qubit());
        }
        let row: Vec<QubitId> = (0..width).map(|_| circ.alloc_qubit()).collect();
        circ.cx(x[i], row[0]);

        let skip_and = square_addsub_probe();
        if !skip_and {
            for k in 0..num_cross {
                circ.ccx(x[i], x[i + 1 + k], row[k + 2]);
            }
        }
        add_into(circ, &prod[2 * i..hi], &row);

        if !skip_and {
            for k in 0..num_cross {
                clear_and(circ, &row[k + 2], &x[i], &x[i + 1 + k]);
            }
        }
        circ.cx(x[i], row[0]);
        for q in row {
            circ.zero_and_free(q);
        }
    }
    debug_assert_eq!(prod.len(), 2 * n, "prod must reach 2n after the build");
}

fn square_addsub_enabled() -> bool {
    true
}

fn square_addsub_local_diag() -> bool {
    std::env::var("TLM_SQUARE_ADDSUB_LOCAL").ok().as_deref() == Some("1")
}

fn square_addsub_skip_c() -> bool {
    std::env::var("TLM_SQUARE_ADDSUB_SKIP_C").ok().as_deref() == Some("1")
}

fn square_addsub_probe() -> bool {
    std::env::var("TLM_SQUARE_ADDSUB_PROBE").ok().as_deref() == Some("1")
}

fn symmetric_square_into_prod_reverse(circ: &mut B, x: &[QubitId], mut prod: Vec<QubitId>) {
    let n = x.len();
    if std::env::var("TLM_SQ_TRACE").ok().as_deref() == Some("1") {
        eprintln!("SQ_CALL rev n={n} crosses={}", n * (n - 1) / 2);
    }
    assert_eq!(prod.len(), 2 * n);

    if square_addsub_enabled() && !(square_addsub_skip_c() && n == 129) {
        if square_addsub_local_diag() {
            crate::point_add::arith::square_addsub_local_inverse(circ, x, &prod);
        } else {
            crate::point_add::arith::square_addsub_vented_inverse(circ, x, &prod);
        }
        for q in prod {
            circ.zero_and_free(q);
        }
        return;
    }
    for i in (0..n).rev() {
        let num_cross = n.saturating_sub(i + 1);
        let width = if i == n - 1 { 1 } else { n - i + 1 };
        let row: Vec<QubitId> = (0..width).map(|_| circ.alloc_qubit()).collect();
        circ.cx(x[i], row[0]);
        let skip_and = square_addsub_probe();
        if !skip_and {
            for k in 0..num_cross {
                circ.ccx(x[i], x[i + 1 + k], row[k + 2]);
            }
        }
        let hi = (2 * i + width + 1).min(prod.len());

        for q in &prod[2 * i..hi] {
            circ.x(*q);
        }
        add_into(circ, &prod[2 * i..hi], &row);
        for q in &prod[2 * i..hi] {
            circ.x(*q);
        }

        if !skip_and {
            for k in 0..num_cross {
                clear_and(circ, &row[k + 2], &x[i], &x[i + 1 + k]);
            }
        }
        circ.cx(x[i], row[0]);
        for q in row {
            circ.zero_and_free(q);
        }

        let keep = (n + i + 1).min(2 * n);
        while prod.len() > keep {
            circ.zero_and_free(prod.pop().unwrap());
        }
    }
    for q in prod {
        circ.zero_and_free(q);
    }
}

fn alloc_zeroes(circ: &mut B, n: usize) -> Vec<QubitId> {
    (0..n).map(|_| circ.alloc_qubit()).collect()
}

fn free_zeroes(circ: &mut B, qs: Vec<QubitId>) {
    for q in qs {
        circ.zero_and_free(q);
    }
}

fn flipped(op: ShiftOp) -> ShiftOp {
    match op {
        ShiftOp::Add => ShiftOp::Sub,
        ShiftOp::Sub => ShiftOp::Add,
    }
}

fn apply_full_width(circ: &mut B, operand: &[QubitId], output_reg: &[QubitId], op: ShiftOp) {
    assert_eq!(operand.len(), N, "full-width modular operand must be 256 bits");
    match op {
        ShiftOp::Add => mod_add_lowpeak(circ, operand, output_reg),
        ShiftOp::Sub => mod_sub(circ, operand, output_reg),
    }
}

fn apply_unshifted_value(circ: &mut B, value: &[QubitId], output_reg: &[QubitId], op: ShiftOp) {
    assert!(value.len() <= N, "unshifted value must fit in 256 bits");
    let pads = alloc_zeroes(circ, N - value.len());
    let mut operand = Vec::with_capacity(N);
    operand.extend_from_slice(value);
    operand.extend_from_slice(&pads);
    apply_full_width(circ, &operand, output_reg, op);
    free_zeroes(circ, pads);
}

fn apply_shifted_value_direct(
    circ: &mut B,
    value: &[QubitId],
    output_reg: &[QubitId],
    shift: usize,
    op: ShiftOp,
) {
    assert!(value.len() + shift <= N, "shifted value must fit in 256 bits");
    let low_pads = alloc_zeroes(circ, shift);
    let high_pads = alloc_zeroes(circ, N - shift - value.len());
    let mut operand = Vec::with_capacity(N);
    operand.extend_from_slice(&low_pads);
    operand.extend_from_slice(value);
    operand.extend_from_slice(&high_pads);
    apply_full_width(circ, &operand, output_reg, op);
    free_zeroes(circ, high_pads);
    free_zeroes(circ, low_pads);
}

fn apply_shifted_value_low(
    circ: &mut B,
    value: &[QubitId],
    output_reg: &[QubitId],
    shift: usize,
    op: ShiftOp,
) {
    assert!(value.len() + shift <= N, "shifted value must fit in 256 bits");
    if shift == 0 {
        apply_unshifted_value(circ, value, output_reg, op);
        return;
    }

    let high_pads = alloc_zeroes(circ, N - shift - value.len());
    let mut operand = Vec::with_capacity(N - shift);
    operand.extend_from_slice(value);
    operand.extend_from_slice(&high_pads);
    match op {
        ShiftOp::Add => mod_add_shifted_low(circ, &operand, output_reg, shift),
        ShiftOp::Sub => mod_sub_shifted_low(circ, &operand, output_reg, shift),
    }
    free_zeroes(circ, high_pads);
}

fn env_tag_enabled(var: &str, tag: &str) -> bool {
    std::env::var(var)
        .ok()
        .map(|tags| tags.split(',').any(|t| t.trim() == tag))
        .unwrap_or(false)
}

fn apply_f_times_value_tagged(circ: &mut B, value: &[QubitId], output_reg: &[QubitId], op: ShiftOp, tag: &str) {
    assert!(value.len() <= N, "f-fold value must fit in 256 bits");
    if value.len() + 32 <= N
        && (std::env::var("TLM_SQUARE_F_RAMP10_DIRECT32").ok().as_deref() == Some("1")
            || env_tag_enabled("TLM_SQUARE_F_RAMP10_DIRECT32_TAGS", tag))
    {
        let pads = alloc_zeroes(circ, N + 1 - value.len());
        let mut ext = Vec::with_capacity(N + 1);
        ext.extend_from_slice(value);
        ext.extend_from_slice(&pads);

        let mut shifted = 0usize;
        for &(shift, sub_f_op) in &F_NAF_TERMS {
            let term_op = match op {
                ShiftOp::Sub => sub_f_op,
                ShiftOp::Add => flipped(sub_f_op),
            };
            if shift == 32 {
                continue;
            }
            while shifted < shift {
                arith::mod_double(circ, &ext);
                shifted += 1;
            }
            apply_full_width(circ, &ext[..N], output_reg, term_op);
        }
        while shifted > 0 {
            arith::mod_double_reverse(circ, &ext);
            shifted -= 1;
        }
        free_zeroes(circ, pads);

        let term_op = match op {
            ShiftOp::Sub => ShiftOp::Sub,
            ShiftOp::Add => ShiftOp::Add,
        };
        apply_shifted_value_direct(circ, value, output_reg, 32, term_op);
        return;
    }

    if env_tag_enabled("TLM_SQUARE_F_DIRECT_TAGS", tag) && value.len() + 32 <= N {
        for &(shift, sub_f_op) in &F_NAF_TERMS {
            let term_op = match op {
                ShiftOp::Sub => sub_f_op,
                ShiftOp::Add => flipped(sub_f_op),
            };
            apply_shifted_value_direct(circ, value, output_reg, shift, term_op);
        }
        return;
    }

    if std::env::var("TLM_SQUARE_F_SHIFTED_LOW").ok().as_deref() == Some("1")
        && value.len() + 32 <= N
    {
        for &(shift, sub_f_op) in &F_NAF_TERMS {
            let term_op = match op {
                ShiftOp::Sub => sub_f_op,
                ShiftOp::Add => flipped(sub_f_op),
            };
            apply_shifted_value_low(circ, value, output_reg, shift, term_op);
        }
        return;
    }

    if std::env::var("TLM_SQUARE_F_DIRECT_SHIFT").ok().as_deref() == Some("1")
        && value.len() + 32 <= N
    {
        for &(shift, sub_f_op) in &F_NAF_TERMS {
            let term_op = match op {
                ShiftOp::Sub => sub_f_op,
                ShiftOp::Add => flipped(sub_f_op),
            };
            apply_shifted_value_direct(circ, value, output_reg, shift, term_op);
        }
        return;
    }

    if value.len() == N {
        for &(shift, sub_f_op) in &F_NAF_TERMS {
            let term_op = match op {
                ShiftOp::Sub => sub_f_op,
                ShiftOp::Add => flipped(sub_f_op),
            };
            apply_shifted_hi_term(circ, value, output_reg, shift, term_op);
        }
        return;
    }

    let pads = alloc_zeroes(circ, N + 1 - value.len());
    let mut ext = Vec::with_capacity(N + 1);
    ext.extend_from_slice(value);
    ext.extend_from_slice(&pads);

    let mut shifted = 0usize;
    for &(shift, sub_f_op) in &F_NAF_TERMS {
        while shifted < shift {
            arith::mod_double(circ, &ext);
            shifted += 1;
        }
        let term_op = match op {
            ShiftOp::Sub => sub_f_op,
            ShiftOp::Add => flipped(sub_f_op),
        };
        apply_full_width(circ, &ext[..N], output_reg, term_op);
    }
    while shifted > 0 {
        arith::mod_double_reverse(circ, &ext);
        shifted -= 1;
    }

    free_zeroes(circ, pads);
}

fn apply_f_times_value(circ: &mut B, value: &[QubitId], output_reg: &[QubitId], op: ShiftOp) {
    apply_f_times_value_tagged(circ, value, output_reg, op, "generic");
}

fn apply_shifted_128_tagged(circ: &mut B, value: &[QubitId], output_reg: &[QubitId], op: ShiftOp, tag: &str) {
    assert!(value.len() <= N + 2, "128-shifted half product must be at most 258 bits");
    let low_len = value.len().min(128);
    let low_pads = alloc_zeroes(circ, 128);
    let high_pads = alloc_zeroes(circ, 128 - low_len);
    let mut operand = Vec::with_capacity(N);
    operand.extend_from_slice(&low_pads);
    operand.extend_from_slice(&value[..low_len]);
    operand.extend_from_slice(&high_pads);
    apply_full_width(circ, &operand, output_reg, op);
    free_zeroes(circ, high_pads);
    free_zeroes(circ, low_pads);

    if value.len() > 128 {
        if matches!(tag, "a" | "b" | "c") {
            arith::with_shifted_square_ffg_prefix_scope(|| {
                apply_f_times_value_tagged(circ, &value[128..], output_reg, op, tag);
            });
        } else {
            apply_f_times_value_tagged(circ, &value[128..], output_reg, op, tag);
        }
    }
}

fn build_sum_hi_lo(circ: &mut B, lambda: &[QubitId]) -> Vec<QubitId> {
    let sum = alloc_zeroes(circ, 129);
    for i in 0..128 {
        circ.cx(lambda[i], sum[i]);
    }
    cuccaro_carry(circ, None, &lambda[128..N], &sum[..128], None, Some(&sum[128]));
    sum
}

fn unbuild_sum_hi_lo(circ: &mut B, lambda: &[QubitId], sum: Vec<QubitId>) {
    let hi_pad = circ.alloc_qubit();
    let mut hi_ext = Vec::with_capacity(129);
    hi_ext.extend_from_slice(&lambda[128..N]);
    hi_ext.push(hi_pad);

    for q in &sum {
        circ.x(*q);
    }
    cuccaro_carry(circ, None, &hi_ext, &sum, None, None);
    for q in &sum {
        circ.x(*q);
    }

    circ.zero_and_free(hi_pad);
    for i in 0..128 {
        circ.cx(lambda[i], sum[i]);
    }
    free_zeroes(circ, sum);
}

pub fn mod_square_sub_pm_secp256k1_symmetric(circ: &mut B, lambda: &[QubitId], output_reg: &[QubitId]) {
    let n = N;
    assert_eq!(lambda.len(), n, "lambda must be n=256 bits (< q)");
    assert_eq!(output_reg.len(), n, "output must be n=256 bits (< q)");

    circ.set_phase("square_sum_hi_lo");
    let sum = build_sum_hi_lo(circ, lambda);

    circ.set_phase("square_c_sum_build");
    let mut c_prod: Vec<QubitId> = Vec::with_capacity(2 * sum.len());
    symmetric_square_into_prod(circ, &sum, &mut c_prod);
    circ.set_phase("square_c_sum_apply_shifted_128_sub");
    apply_shifted_128_tagged(circ, &c_prod, output_reg, ShiftOp::Sub, "c");
    circ.set_phase("square_c_sum_unbuild");
    symmetric_square_into_prod_reverse(circ, &sum, c_prod);

    circ.set_phase("square_a_lo_build");
    let lo = &lambda[..128];
    let mut a_prod: Vec<QubitId> = Vec::with_capacity(2 * lo.len());
    symmetric_square_into_prod(circ, lo, &mut a_prod);
    circ.set_phase("square_a_lo_apply_unshifted_sub");
    apply_unshifted_value(circ, &a_prod, output_reg, ShiftOp::Sub);
    circ.set_phase("square_a_lo_apply_shifted_128_add");
    apply_shifted_128_tagged(circ, &a_prod, output_reg, ShiftOp::Add, "a");
    circ.set_phase("square_a_lo_unbuild");
    symmetric_square_into_prod_reverse(circ, lo, a_prod);

    circ.set_phase("square_b_hi_build");
    let hi = &lambda[128..N];
    let mut b_prod: Vec<QubitId> = Vec::with_capacity(2 * hi.len());
    symmetric_square_into_prod(circ, hi, &mut b_prod);
    circ.set_phase("square_b_hi_apply_shifted_128_add");
    apply_shifted_128_tagged(circ, &b_prod, output_reg, ShiftOp::Add, "b");
    circ.set_phase("square_b_hi_apply_f_times_sub");
    apply_f_times_value(circ, &b_prod, output_reg, ShiftOp::Sub);
    circ.set_phase("square_b_hi_unbuild");
    symmetric_square_into_prod_reverse(circ, hi, b_prod);

    circ.set_phase("square_sum_hi_lo_unbuild");
    unbuild_sum_hi_lo(circ, lambda, sum);
}
