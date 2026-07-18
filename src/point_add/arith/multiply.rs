
use super::*;

#[allow(dead_code)]
pub(crate) fn mod_mul_write_into_zero_acc_schoolbook_lowq(
    b: &mut B,
    acc: &[QubitId],
    x: &[QubitId],
    y: &[QubitId],
    p: U256,
) {
    let n = acc.len();
    debug_assert_eq!(n, 256);

    let tmp_ext = b.alloc_qubits(2 * n);
    schoolbook_mul_into_addsub_lowq(b, x, y, &tmp_ext);

    let lo: Vec<QubitId> = tmp_ext[0..n].to_vec();
    let hi: Vec<QubitId> = tmp_ext[n..2 * n].to_vec();
    mod_add_qq_fast_from_zero(b, acc, &lo, p);
    mod_add_qq_fast(b, acc, &hi, p);
    for _ in 0..4 {
        mod_double_inplace_fast(b, &hi, p);
    }
    mod_add_qq_fast(b, acc, &hi, p);
    for _ in 0..2 {
        mod_double_inplace_fast(b, &hi, p);
    }
    mod_sub_qq_fast(b, acc, &hi, p);
    for _ in 0..4 {
        mod_double_inplace_fast(b, &hi, p);
    }
    mod_add_qq_fast(b, acc, &hi, p);
    let (spill, flag_inv, ovf) = mod_shift_left_by_k(b, &hi, p, 22);
    mod_add_qq(b, acc, &hi, p);
    mod_shift_right_by_k(b, &hi, p, 22, spill, flag_inv, ovf);
    for _ in 0..10 {
        mod_halve_inplace_fast(b, &hi, p);
    }

    schoolbook_mul_into_addsub_lowq_inverse(b, x, y, &tmp_ext);
    b.free_vec(&tmp_ext);
}

pub(crate) fn controlled_add_subtract_lowq(b: &mut B, x: &[QubitId], acc: &[QubitId], ctrl: QubitId) {
    let n = x.len();
    debug_assert_eq!(acc.len(), n + 1);

    let pad = b.alloc_qubit();
    let mut x_ext = x.to_vec();
    x_ext.push(pad);

    let c_in = b.alloc_qubit();

    b.x(ctrl);
    for i in 0..n {
        b.cx(ctrl, x_ext[i]);
    }
    b.cx(ctrl, c_in);

    cuccaro_add(b, &x_ext, acc, c_in);

    b.cx(ctrl, c_in);
    for i in 0..n {
        b.cx(ctrl, x_ext[i]);
    }
    b.x(ctrl);

    b.free(c_in);
    b.free(pad);
}

pub(crate) fn controlled_add_subtract_lowq_inverse(b: &mut B, x: &[QubitId], acc: &[QubitId], ctrl: QubitId) {
    let n = x.len();
    debug_assert_eq!(acc.len(), n + 1);

    let pad = b.alloc_qubit();
    let mut x_ext = x.to_vec();
    x_ext.push(pad);

    let c_in = b.alloc_qubit();

    b.x(ctrl);
    for i in 0..n {
        b.cx(ctrl, x_ext[i]);
    }
    b.cx(ctrl, c_in);

    cuccaro_sub(b, &x_ext, acc, c_in);

    b.cx(ctrl, c_in);
    for i in 0..n {
        b.cx(ctrl, x_ext[i]);
    }
    b.x(ctrl);

    b.free(c_in);
    b.free(pad);
}

pub(crate) fn schoolbook_mul_into_addsub_lowq(b: &mut B, x: &[QubitId], y: &[QubitId], tmp_ext: &[QubitId]) {
    let n = x.len();
    debug_assert_eq!(y.len(), n);
    debug_assert_eq!(tmp_ext.len(), 2 * n);

    let low = b.alloc_qubit();
    let mut wide: Vec<QubitId> = Vec::with_capacity(2 * n + 1);
    wide.push(low);
    wide.extend_from_slice(tmp_ext);

    for k in 0..n {
        let slice: Vec<QubitId> = wide[k..k + n + 1].to_vec();
        controlled_add_subtract_lowq(b, x, &slice, y[k]);
    }

    {
        let pad = b.alloc_qubit();
        let mut y_ext = y.to_vec();
        y_ext.push(pad);
        let slice: Vec<QubitId> = wide[n..2 * n + 1].to_vec();
        let c_in = b.alloc_qubit();
        b.x(c_in);
        cuccaro_add(b, &y_ext, &slice, c_in);
        b.x(c_in);
        b.free(c_in);
        b.free(pad);
    }

    b.x(wide[2 * n]);

    {
        let mut x_ext: Vec<QubitId> = x.to_vec();
        while x_ext.len() < 2 * n + 1 {
            x_ext.push(b.alloc_qubit());
        }
        let c_in = b.alloc_qubit();
        cuccaro_sub(b, &x_ext, &wide, c_in);
        b.free(c_in);
        for _ in n..2 * n + 1 {
            let q = x_ext.pop().unwrap();
            b.free(q);
        }
    }

    {
        let pad = b.alloc_qubit();
        let mut x_ext = x.to_vec();
        x_ext.push(pad);
        let slice: Vec<QubitId> = wide[n..2 * n + 1].to_vec();
        let c_in = b.alloc_qubit();
        cuccaro_add(b, &x_ext, &slice, c_in);
        b.free(c_in);
        b.free(pad);
    }

    b.free(low);
}

pub(crate) fn schoolbook_mul_into_addsub_lowq_inverse(
    b: &mut B,
    x: &[QubitId],
    y: &[QubitId],
    tmp_ext: &[QubitId],
) {
    let n = x.len();
    debug_assert_eq!(y.len(), n);
    debug_assert_eq!(tmp_ext.len(), 2 * n);

    let low = b.alloc_qubit();
    let mut wide: Vec<QubitId> = Vec::with_capacity(2 * n + 1);
    wide.push(low);
    wide.extend_from_slice(tmp_ext);

    {
        let pad = b.alloc_qubit();
        let mut x_ext = x.to_vec();
        x_ext.push(pad);
        let slice: Vec<QubitId> = wide[n..2 * n + 1].to_vec();
        let c_in = b.alloc_qubit();
        cuccaro_sub(b, &x_ext, &slice, c_in);
        b.free(c_in);
        b.free(pad);
    }

    {
        let mut x_ext: Vec<QubitId> = x.to_vec();
        while x_ext.len() < 2 * n + 1 {
            x_ext.push(b.alloc_qubit());
        }
        let c_in = b.alloc_qubit();
        cuccaro_add(b, &x_ext, &wide, c_in);
        b.free(c_in);
        for _ in n..2 * n + 1 {
            let q = x_ext.pop().unwrap();
            b.free(q);
        }
    }

    b.x(wide[2 * n]);

    {
        let pad = b.alloc_qubit();
        let mut y_ext = y.to_vec();
        y_ext.push(pad);
        let slice: Vec<QubitId> = wide[n..2 * n + 1].to_vec();
        let c_in = b.alloc_qubit();
        b.x(c_in);
        cuccaro_sub(b, &y_ext, &slice, c_in);
        b.x(c_in);
        b.free(c_in);
        b.free(pad);
    }
    for k in (0..n).rev() {
        let slice: Vec<QubitId> = wide[k..k + n + 1].to_vec();
        controlled_add_subtract_lowq_inverse(b, x, &slice, y[k]);
    }

    b.free(low);
}

pub(crate) fn karatsuba_half_sum_compute(b: &mut B, lo: &[QubitId], hi: &[QubitId], acc: &[QubitId]) {
    let h = lo.len();
    debug_assert_eq!(h, hi.len());
    debug_assert_eq!(acc.len(), h + 1);
    for i in 0..h {
        b.cx(lo[i], acc[i]);
    }
    let hi_pad = b.alloc_qubit();
    let mut hi_ext = hi.to_vec();
    hi_ext.push(hi_pad);
    add_nbit_qq_fast(b, &hi_ext, acc);
    b.free(hi_pad);
}

pub(crate) fn karatsuba_half_sum_uncompute(b: &mut B, lo: &[QubitId], hi: &[QubitId], acc: &[QubitId]) {
    let h = lo.len();
    let hi_pad = b.alloc_qubit();
    let mut hi_ext = hi.to_vec();
    hi_ext.push(hi_pad);
    sub_nbit_qq_fast(b, &hi_ext, acc);
    b.free(hi_pad);
    for i in 0..h {
        b.cx(lo[i], acc[i]);
    }
}

pub(crate) fn schoolbook_square_symmetric(b: &mut B, x: &[QubitId], tmp_ext: &[QubitId]) {
    let n = x.len();
    debug_assert_eq!(tmp_ext.len(), 2 * n);
    for i in 0..n {

        let width = if i == n - 1 { 1 } else { n - i + 1 };
        let num_cross = if i + 1 < n { n - i - 1 } else { 0 };

        let row = b.alloc_qubits(width);
        b.cx(x[i], row[0]);
        for k in 0..num_cross {
            b.ccx(x[i], x[i + 1 + k], row[k + 2]);
        }
        let pad = b.alloc_qubit();
        let mut row_padded = row.clone();
        row_padded.push(pad);
        let slice: Vec<QubitId> = tmp_ext[2 * i..2 * i + width + 1].to_vec();
        let c_in = b.alloc_qubit();
        cuccaro_add_fast(b, &row_padded, &slice, c_in);
        b.free(c_in);
        b.free(pad);
        b.cx(x[i], row[0]);
        for k in 0..num_cross {
            let m = b.alloc_bit();
            b.hmr(row[k + 2], m);
            b.cz_if(x[i], x[i + 1 + k], m);
        }
        b.free_vec(&row);
    }
}

pub(crate) fn schoolbook_square_symmetric_inverse(b: &mut B, x: &[QubitId], tmp_ext: &[QubitId]) {
    let n = x.len();
    for i in (0..n).rev() {
        let width = if i == n - 1 { 1 } else { n - i + 1 };
        let num_cross = if i + 1 < n { n - i - 1 } else { 0 };
        let row = b.alloc_qubits(width);
        b.cx(x[i], row[0]);
        for k in 0..num_cross {
            b.ccx(x[i], x[i + 1 + k], row[k + 2]);
        }
        let pad = b.alloc_qubit();
        let mut row_padded = row.clone();
        row_padded.push(pad);
        let slice: Vec<QubitId> = tmp_ext[2 * i..2 * i + width + 1].to_vec();
        let c_in = b.alloc_qubit();
        cuccaro_sub_fast(b, &row_padded, &slice, c_in);
        b.free(c_in);
        b.free(pad);
        b.cx(x[i], row[0]);
        for k in 0..num_cross {
            let m = b.alloc_bit();
            b.hmr(row[k + 2], m);
            b.cz_if(x[i], x[i + 1 + k], m);
        }
        b.free_vec(&row);
    }
}

pub(crate) fn schoolbook_square_symmetric_lowq(b: &mut B, x: &[QubitId], tmp_ext: &[QubitId]) {
    let n = x.len();
    debug_assert_eq!(tmp_ext.len(), 2 * n);
    for i in 0..n {
        let width = if i == n - 1 { 1 } else { n - i + 1 };
        let num_cross = if i + 1 < n { n - i - 1 } else { 0 };
        let row = b.alloc_qubits(width);
        b.cx(x[i], row[0]);
        for k in 0..num_cross {
            b.ccx(x[i], x[i + 1 + k], row[k + 2]);
        }
        let pad = b.alloc_qubit();
        let mut row_padded = row.clone();
        row_padded.push(pad);
        let slice: Vec<QubitId> = tmp_ext[2 * i..2 * i + width + 1].to_vec();
        let c_in = b.alloc_qubit();
        cuccaro_add(b, &row_padded, &slice, c_in);
        b.free(c_in);
        b.free(pad);
        b.cx(x[i], row[0]);
        for k in 0..num_cross {
            let m = b.alloc_bit();
            b.hmr(row[k + 2], m);
            b.cz_if(x[i], x[i + 1 + k], m);
        }
        b.free_vec(&row);
    }
}

pub(crate) fn schoolbook_square_symmetric_lowq_inverse(b: &mut B, x: &[QubitId], tmp_ext: &[QubitId]) {
    let n = x.len();
    for i in (0..n).rev() {
        let width = if i == n - 1 { 1 } else { n - i + 1 };
        let num_cross = if i + 1 < n { n - i - 1 } else { 0 };
        let row = b.alloc_qubits(width);
        b.cx(x[i], row[0]);
        for k in 0..num_cross {
            b.ccx(x[i], x[i + 1 + k], row[k + 2]);
        }
        let pad = b.alloc_qubit();
        let mut row_padded = row.clone();
        row_padded.push(pad);
        let slice: Vec<QubitId> = tmp_ext[2 * i..2 * i + width + 1].to_vec();
        let c_in = b.alloc_qubit();
        cuccaro_sub(b, &row_padded, &slice, c_in);
        b.free(c_in);
        b.free(pad);
        b.cx(x[i], row[0]);
        for k in 0..num_cross {
            let m = b.alloc_bit();
            b.hmr(row[k + 2], m);
            b.cz_if(x[i], x[i + 1 + k], m);
        }
        b.free_vec(&row);
    }
}

pub(crate) fn schoolbook_square_symmetric_hosted(
    b: &mut B,
    x: &[QubitId],
    tmp_ext: &[QubitId],
    host: &[QubitId],
) {
    let n = x.len();
    debug_assert_eq!(tmp_ext.len(), 2 * n);
    if square_selfhost_safe_lane_reuse_enabled() {
        assert_qubit_slices_disjoint(&[x, tmp_ext, host]);
    }
    for i in 0..n {
        let width = if i == n - 1 { 1 } else { n - i + 1 };
        let num_cross = if i + 1 < n { n - i - 1 } else { 0 };
        let row = b.alloc_qubits(width);
        b.cx(x[i], row[0]);
        for k in 0..num_cross {
            b.ccx(x[i], x[i + 1 + k], row[k + 2]);
        }
        let slice: Vec<QubitId> = tmp_ext[2 * i..2 * i + width + 1].to_vec();
        if square_selfhost_safe_lane_reuse_enabled() {

            assert!(host.len() > width);
            cuccaro_add_fast_low_to_ext_borrowed_carries(
                b,
                &row,
                &slice,
                host[width],
                &host[..width],
            );
        } else {
            let pad = b.alloc_qubit();
            let mut row_padded = row.clone();
            row_padded.push(pad);
            let c_in = b.alloc_qubit();
            cuccaro_add_fast_borrowed_carries(
                b,
                &row_padded,
                &slice,
                c_in,
                &host[..row_padded.len() - 1],
            );
            b.free(c_in);
            b.free(pad);
        }
        b.cx(x[i], row[0]);
        for k in 0..num_cross {
            let m = b.alloc_bit();
            b.hmr(row[k + 2], m);
            b.cz_if(x[i], x[i + 1 + k], m);
        }
        b.free_vec(&row);
    }
}

pub(crate) fn schoolbook_square_symmetric_hosted_inverse(
    b: &mut B,
    x: &[QubitId],
    tmp_ext: &[QubitId],
    host: &[QubitId],
) {
    let n = x.len();
    if square_selfhost_safe_lane_reuse_enabled() {
        assert_qubit_slices_disjoint(&[x, tmp_ext, host]);
    }
    for i in (0..n).rev() {
        let width = if i == n - 1 { 1 } else { n - i + 1 };
        let num_cross = if i + 1 < n { n - i - 1 } else { 0 };
        let row = b.alloc_qubits(width);
        b.cx(x[i], row[0]);
        for k in 0..num_cross {
            b.ccx(x[i], x[i + 1 + k], row[k + 2]);
        }
        let slice: Vec<QubitId> = tmp_ext[2 * i..2 * i + width + 1].to_vec();
        if square_selfhost_safe_lane_reuse_enabled() {
            assert!(host.len() > width);
            cuccaro_sub_fast_low_to_ext_borrowed_carries(
                b,
                &row,
                &slice,
                host[width],
                &host[..width],
            );
        } else {
            let pad = b.alloc_qubit();
            let mut row_padded = row.clone();
            row_padded.push(pad);
            let c_in = b.alloc_qubit();
            cuccaro_sub_fast_borrowed_carries(
                b,
                &row_padded,
                &slice,
                c_in,
                &host[..row_padded.len() - 1],
            );
            b.free(c_in);
            b.free(pad);
        }
        b.cx(x[i], row[0]);
        for k in 0..num_cross {
            let m = b.alloc_bit();
            b.hmr(row[k + 2], m);
            b.cz_if(x[i], x[i + 1 + k], m);
        }
        b.free_vec(&row);
    }
}

pub(crate) fn square_selfhost_safe_lane_reuse_enabled() -> bool {
    std::env::var("SQUARE_SELFHOST_SAFE_LANE_REUSE")
        .ok()
        .as_deref()
        == Some("1")
}

pub(crate) fn assert_qubit_slices_disjoint(slices: &[&[QubitId]]) {
    let mut seen = std::collections::BTreeSet::new();
    for slice in slices {
        for &q in *slice {
            assert!(seen.insert(q), "scratch lane q{} aliases an operand", q.0);
        }
    }
}

pub(crate) fn square_selfhost_gate_suffix_carries(n: usize) -> usize {
    std::env::var("SQUARE_SELFHOST_GATE_SUFFIX_CARRIES")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(0)
        .min(n.saturating_sub(1))
}

pub(crate) fn square_row_windows() -> usize {
    std::env::var("SQUARE_ROW_WINDOWS")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(0)
}

fn square_row_window_min_width() -> usize {
    std::env::var("SQUARE_ROW_WINDOW_MIN_WIDTH")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(96)
}

fn square_row_max_seg() -> usize {
    std::env::var("SQUARE_ROW_MAX_SEG")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(0)
}

fn square_cleanup_direction(raw: &str) -> Option<bool> {
    match raw.trim().to_ascii_lowercase().as_str() {
        "f" | "forward" | "false" | "0" => Some(false),
        "r" | "reverse" | "true" | "1" => Some(true),
        _ => None,
    }
}

fn square_row_window_clean_compare_bits(
    row: usize,
    window: usize,
    reverse: bool,
) -> usize {
    let default_bits = std::env::var("SQUARE_ROW_WINDOW_CLEAN_COMPARE_BITS")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(0);
    let row_bits = std::env::var("SQUARE_ROW_WINDOW_CLEAN_ROW_BITS")
        .ok()
        .and_then(|spec| {
            spec.split(',').rev().find_map(|item| {
                let (raw_row, raw_bits) = item.trim().split_once(':')?;
                if raw_row.trim().parse::<usize>().ok()? != row {
                    return None;
                }
                raw_bits
                    .trim()
                    .parse::<usize>()
                    .ok()
                    .filter(|bits| (1..=N).contains(bits))
            })
        })
        .unwrap_or(default_bits);
    let Ok(spec) = std::env::var("SQUARE_ROW_WINDOW_CLEAN_SITE_BITS") else {
        return row_bits;
    };
    for item in spec.split(',').rev() {
        let fields: Vec<_> = item.trim().split(':').map(str::trim).collect();
        if fields.len() != 4 {
            continue;
        }
        if fields[0].parse::<usize>().ok() != Some(row)
            || fields[1].parse::<usize>().ok() != Some(window)
            || square_cleanup_direction(fields[2]) != Some(reverse)
        {
            continue;
        }
        if let Ok(bits) = fields[3].parse::<usize>() {
            if (1..=N).contains(&bits) {
                return bits;
            }
        }
    }
    row_bits
}

fn square_row_window_measured_carry_clear_enabled() -> bool {
    std::env::var("SQUARE_ROW_WINDOW_MEASURED_CARRY_CLEAR")
        .ok()
        .as_deref()
        == Some("1")
}

fn square_row_bit_set(b: &mut B, x: &[QubitId], i: usize, j: usize, t: QubitId) {
    if j == 0 {
        b.cx(x[i], t);
    } else if j == 1 {

    } else {
        b.ccx(x[i], x[i + 1 + (j - 2)], t);
    }
}

fn square_row_bit_clear_hmr(b: &mut B, x: &[QubitId], i: usize, j: usize, t: QubitId) {
    if j == 0 {
        b.cx(x[i], t);
    } else if j == 1 {

    } else {
        let m = b.alloc_bit();
        b.hmr(t, m);
        b.cz_if(x[i], x[i + 1 + (j - 2)], m);
    }
}

fn square_row_windowed_apply(
    b: &mut B,
    x: &[QubitId],
    tmp_ext: &[QubitId],
    i: usize,
    width: usize,
    windows: usize,
    forward: bool,
) {
    let base = 2 * i;
    let windows = windows.max(1).min(width);

    let bounds: Vec<(usize, usize)> = (0..windows)
        .map(|w| {
            let lo = (w * width) / windows;
            let hi = ((w + 1) * width) / windows;
            (lo, hi)
        })
        .filter(|&(lo, hi)| hi > lo)
        .collect();
    let nwin = bounds.len();

    let build_seg = |b: &mut B, lo: usize, hi: usize| -> Vec<QubitId> {
        let seg = b.alloc_qubits(hi - lo);
        for (k, &q) in seg.iter().enumerate() {
            square_row_bit_set(b, x, i, lo + k, q);
        }
        seg
    };
    let clear_seg = |b: &mut B, lo: usize, seg: &[QubitId]| {
        for (k, &q) in seg.iter().enumerate() {
            square_row_bit_clear_hmr(b, x, i, lo + k, q);
        }
        b.free_vec(seg);
    };

    let row_top = base + width + 1;
    let borrow_lane = |b: &mut B, _need: usize| -> Vec<QubitId> {

        tmp_ext[row_top..row_top + _need].to_vec()
    };

    let mut carry_in = b.alloc_qubit();
    let first_carry = carry_in;
    let mut couts: Vec<(QubitId, usize, usize, QubitId, usize)> = Vec::new();
    for (wi, &(lo, hi)) in bounds.iter().enumerate() {
        let last = wi == nwin - 1;
        let seg = build_seg(b, lo, hi);
        let seg_w = hi - lo;

        let pad = b.alloc_qubit();
        let mut a_block = seg.clone();
        a_block.push(pad);
        let high = if last {

            tmp_ext[base + hi]
        } else {
            b.alloc_qubit()
        };
        let mut acc_block: Vec<QubitId> = tmp_ext[base + lo..base + hi].to_vec();
        acc_block.push(high);
        let nblk = a_block.len();
        let carries = borrow_lane(b, nblk - 1);
        if forward {
            cuccaro_add_fast_borrowed_carries(b, &a_block, &acc_block, carry_in, &carries);
        } else {
            cuccaro_sub_fast_borrowed_carries(b, &a_block, &acc_block, carry_in, &carries);
        }
        b.free(pad);
        if last {

        } else {
            couts.push((high, lo, hi, carry_in, wi));
            carry_in = high;
        }
        clear_seg(b, lo, &seg);
    }

    let slow_cmp = std::env::var("SQUARE_ROW_WINDOW_SLOW_CMP").ok().as_deref() == Some("1");
    let measured_clear = square_row_window_measured_carry_clear_enabled();
    for &(cout, lo, hi, cin, window) in couts.iter().rev() {
        let clean_cmp_bits =
            square_row_window_clean_compare_bits(i, window, !forward);
        let seg_w = hi - lo;
        let trunc_w = if clean_cmp_bits == 0 {
            seg_w
        } else {
            clean_cmp_bits.min(seg_w)
        };
        if trunc_w < seg_w {
            let suffix_lo = hi - trunc_w;
            let seg = build_seg(b, suffix_lo, hi);
            let carries = tmp_ext[row_top..row_top + trunc_w].to_vec();
            let cmp_cin = b.alloc_qubit();
            if forward {
                if measured_clear {
                    let phase = b.alloc_bit();
                    b.hmr(cout, phase);
                    cmp_lt_phase_conditioned_with_cin_borrowed_carries(
                        b,
                        &tmp_ext[base + suffix_lo..base + hi],
                        &seg,
                        cmp_cin,
                        &carries,
                        phase,
                    );
                } else {
                    cmp_lt_into_fast_with_cin_borrowed_carries(
                        b,
                        &tmp_ext[base + suffix_lo..base + hi],
                        &seg,
                        cmp_cin,
                        cout,
                        &carries,
                    );
                }
            } else {
                for &q in &seg {
                    b.x(q);
                }
                if measured_clear {
                    let phase = b.alloc_bit();
                    b.hmr(cout, phase);
                    cmp_lt_phase_conditioned_with_cin_borrowed_carries(
                        b,
                        &seg,
                        &tmp_ext[base + suffix_lo..base + hi],
                        cmp_cin,
                        &carries,
                        phase,
                    );
                } else {
                    cmp_lt_into_fast_with_cin_borrowed_carries(
                        b,
                        &seg,
                        &tmp_ext[base + suffix_lo..base + hi],
                        cmp_cin,
                        cout,
                        &carries,
                    );
                }
                for &q in &seg {
                    b.x(q);
                }
            }
            b.free(cmp_cin);
            clear_seg(b, suffix_lo, &seg);
        } else {
            let seg = build_seg(b, lo, hi);
            let carries = tmp_ext[row_top..row_top + seg_w].to_vec();
            if forward {

                if measured_clear {
                    let phase = b.alloc_bit();
                    b.hmr(cout, phase);
                    cmp_lt_phase_conditioned_with_cin_borrowed_carries(
                        b,
                        &tmp_ext[base + lo..base + hi],
                        &seg,
                        cin,
                        &carries,
                        phase,
                    );
                } else if slow_cmp {
                    cmp_lt_into_with_cin_slow(b, &tmp_ext[base + lo..base + hi], &seg, cin, cout);
                } else {
                    cmp_lt_into_fast_with_cin_borrowed_carries(
                        b, &tmp_ext[base + lo..base + hi], &seg, cin, cout, &carries,
                    );
                }
            } else {

                for k in 0..seg_w {
                    b.x(seg[k]);
                }
                if measured_clear {
                    let phase = b.alloc_bit();
                    b.hmr(cout, phase);
                    cmp_lt_phase_conditioned_with_cin_borrowed_carries(
                        b,
                        &seg,
                        &tmp_ext[base + lo..base + hi],
                        cin,
                        &carries,
                        phase,
                    );
                } else if slow_cmp {
                    cmp_lt_into_with_cin_slow(b, &seg, &tmp_ext[base + lo..base + hi], cin, cout);
                } else {
                    cmp_lt_into_fast_with_cin_borrowed_carries(
                        b, &seg, &tmp_ext[base + lo..base + hi], cin, cout, &carries,
                    );
                }
                for k in 0..seg_w {
                    b.x(seg[k]);
                }
            }
            clear_seg(b, lo, &seg);
        }
        b.free(cout);
    }
    b.free(first_carry);
}

pub(crate) fn schoolbook_square_symmetric_lowq_selfhosted(b: &mut B, x: &[QubitId], tmp_ext: &[QubitId]) {
    schoolbook_square_symmetric_lowq_selfhosted_with_clean_supplement(b, x, tmp_ext, &[]);
}

pub(crate) fn schoolbook_square_symmetric_lowq_selfhosted_with_clean_supplement(
    b: &mut B,
    x: &[QubitId],
    tmp_ext: &[QubitId],
    clean_supplement: &[QubitId],
) {
    let n = x.len();
    debug_assert_eq!(tmp_ext.len(), 2 * n);
    let safe_reuse = square_selfhost_safe_lane_reuse_enabled();
    if safe_reuse {
        assert_qubit_slices_disjoint(&[x, tmp_ext, clean_supplement]);
    }
    let gate_prefix_rows = std::env::var("SQUARE_SELFHOST_GATE_PREFIX_ROWS")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(0);
    let row_windows = square_row_windows();
    let row_window_min = square_row_window_min_width();
    let max_seg = square_row_max_seg();
    for i in 0..n {
        let width = if i == n - 1 { 1 } else { n - i + 1 };
        let num_cross = if i + 1 < n { n - i - 1 } else { 0 };
        if max_seg > 0 && i >= gate_prefix_rows && width > max_seg {
            let w = width.div_ceil(max_seg);
            square_row_windowed_apply(b, x, tmp_ext, i, width, w, true);
            continue;
        }
        if max_seg == 0 && row_windows >= 1 && i >= gate_prefix_rows && width >= row_window_min {
            square_row_windowed_apply(b, x, tmp_ext, i, width, row_windows, true);
            continue;
        }
        let row = b.alloc_qubits(width);
        b.cx(x[i], row[0]);
        for k in 0..num_cross {
            b.ccx(x[i], x[i + 1 + k], row[k + 2]);
        }
        let hi = 2 * i + width + 1;
        let slice: Vec<QubitId> = tmp_ext[2 * i..hi].to_vec();
        if i < gate_prefix_rows {
            let pad = b.alloc_qubit();
            let mut row_padded = row.clone();
            row_padded.push(pad);
            let c_in = b.alloc_qubit();
            cuccaro_add(b, &row_padded, &slice, c_in);
            b.free(c_in);
            b.free(pad);
        } else if safe_reuse {
            let need = row.len() - square_selfhost_gate_suffix_carries(row.len());
            let avail = tmp_ext.len() - hi;
            let from_tmp = need.min(avail);
            let from_supplement = (need - from_tmp).min(clean_supplement.len());
            let from_global = need - from_tmp - from_supplement;
            let gpool = b.alloc_qubits(from_global);
            let mut carries: Vec<QubitId> = tmp_ext[hi..hi + from_tmp].to_vec();
            carries.extend_from_slice(&clean_supplement[..from_supplement]);
            carries.extend_from_slice(&gpool);
            cuccaro_add_fast_low_to_ext_borrowed_carries_no_cin(b, &row, &slice, &carries);
            b.free_vec(&gpool);
        } else {
            let pad = b.alloc_qubit();
            let mut row_padded = row.clone();
            row_padded.push(pad);
            let c_in = b.alloc_qubit();
            let need = row_padded.len() - 1;
            let avail = tmp_ext.len() - hi;
            let from_tmp = need.min(avail);
            let from_global = need - from_tmp;
            let gpool = b.alloc_qubits(from_global);
            let mut carries: Vec<QubitId> = tmp_ext[hi..hi + from_tmp].to_vec();
            carries.extend_from_slice(&gpool);
            cuccaro_add_fast_borrowed_carries(b, &row_padded, &slice, c_in, &carries);
            b.free(c_in);
            b.free_vec(&gpool);
            b.free(pad);
        }
        b.cx(x[i], row[0]);
        for k in 0..num_cross {
            let m = b.alloc_bit();
            b.hmr(row[k + 2], m);
            b.cz_if(x[i], x[i + 1 + k], m);
        }
        b.free_vec(&row);
    }
}

pub(crate) fn schoolbook_square_symmetric_lowq_selfhosted_inverse(
    b: &mut B,
    x: &[QubitId],
    tmp_ext: &[QubitId],
) {
    schoolbook_square_symmetric_lowq_selfhosted_inverse_with_clean_supplement(b, x, tmp_ext, &[]);
}

pub(crate) fn schoolbook_square_symmetric_lowq_selfhosted_inverse_with_clean_supplement(
    b: &mut B,
    x: &[QubitId],
    tmp_ext: &[QubitId],
    clean_supplement: &[QubitId],
) {
    let n = x.len();
    debug_assert_eq!(tmp_ext.len(), 2 * n);
    let safe_reuse = square_selfhost_safe_lane_reuse_enabled();
    if safe_reuse {
        assert_qubit_slices_disjoint(&[x, tmp_ext, clean_supplement]);
    }
    let gate_prefix_rows = std::env::var("SQUARE_SELFHOST_GATE_PREFIX_ROWS")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(0);
    let row_windows = square_row_windows();
    let row_window_min = square_row_window_min_width();
    let max_seg = square_row_max_seg();
    for i in (0..n).rev() {
        let width = if i == n - 1 { 1 } else { n - i + 1 };
        let num_cross = if i + 1 < n { n - i - 1 } else { 0 };
        if max_seg > 0 && i >= gate_prefix_rows && width > max_seg {
            let w = width.div_ceil(max_seg);
            square_row_windowed_apply(b, x, tmp_ext, i, width, w, false);
            continue;
        }
        if max_seg == 0 && row_windows >= 1 && i >= gate_prefix_rows && width >= row_window_min {
            square_row_windowed_apply(b, x, tmp_ext, i, width, row_windows, false);
            continue;
        }
        let row = b.alloc_qubits(width);
        b.cx(x[i], row[0]);
        for k in 0..num_cross {
            b.ccx(x[i], x[i + 1 + k], row[k + 2]);
        }
        let hi = 2 * i + width + 1;
        let slice: Vec<QubitId> = tmp_ext[2 * i..hi].to_vec();
        if i < gate_prefix_rows {
            let pad = b.alloc_qubit();
            let mut row_padded = row.clone();
            row_padded.push(pad);
            let c_in = b.alloc_qubit();
            cuccaro_sub(b, &row_padded, &slice, c_in);
            b.free(c_in);
            b.free(pad);
        } else if safe_reuse {
            let need = row.len() - square_selfhost_gate_suffix_carries(row.len());
            let avail = tmp_ext.len() - hi;
            let from_tmp = need.min(avail);
            let from_supplement = (need - from_tmp).min(clean_supplement.len());
            let from_global = need - from_tmp - from_supplement;
            let gpool = b.alloc_qubits(from_global);
            let mut carries: Vec<QubitId> = tmp_ext[hi..hi + from_tmp].to_vec();
            carries.extend_from_slice(&clean_supplement[..from_supplement]);
            carries.extend_from_slice(&gpool);
            cuccaro_sub_fast_low_to_ext_borrowed_carries_no_cin(b, &row, &slice, &carries);
            b.free_vec(&gpool);
        } else {
            let pad = b.alloc_qubit();
            let mut row_padded = row.clone();
            row_padded.push(pad);
            let c_in = b.alloc_qubit();
            let need = row_padded.len() - 1;
            let avail = tmp_ext.len() - hi;
            let from_tmp = need.min(avail);
            let from_global = need - from_tmp;
            let gpool = b.alloc_qubits(from_global);
            let mut carries: Vec<QubitId> = tmp_ext[hi..hi + from_tmp].to_vec();
            carries.extend_from_slice(&gpool);
            cuccaro_sub_fast_borrowed_carries(b, &row_padded, &slice, c_in, &carries);
            b.free(c_in);
            b.free_vec(&gpool);
            b.free(pad);
        }
        b.cx(x[i], row[0]);
        for k in 0..num_cross {
            let m = b.alloc_bit();
            b.hmr(row[k + 2], m);
            b.cz_if(x[i], x[i + 1 + k], m);
        }
        b.free_vec(&row);
    }
}

pub(crate) fn kara_z2_selfhost_enabled() -> bool {
    std::env::var("KARA_Z2_SELFHOST").ok().as_deref() != Some("0")
}

pub(crate) fn xtail_sq_selfhost_enabled() -> bool {
    std::env::var("XTAIL_SQ_SELFHOST").ok().as_deref() != Some("0")
}

fn round84_inplace_solinas_fold_enabled() -> bool {
    std::env::var("ROUND84_INPLACE_SOLINAS_FOLD")
        .ok()
        .as_deref()
        == Some("1")
}

fn round84_fold_fast_add_enabled() -> bool {
    std::env::var("ROUND84_FOLD_FAST_ADD").ok().as_deref() == Some("1")
}
#[inline]
fn round84_add_small(b: &mut B, a: &[QubitId], acc: &[QubitId]) {
    if round84_fold_fast_add_enabled() {
        add_nbit_qq_fast(b, a, acc);
    } else {
        add_nbit_qq(b, a, acc);
    }
}
#[inline]
fn round84_sub_small(b: &mut B, a: &[QubitId], acc: &[QubitId]) {
    if round84_fold_fast_add_enabled() {
        sub_nbit_qq_fast(b, a, acc);
    } else {
        sub_nbit_qq(b, a, acc);
    }
}

fn round84_inplace_quotient_carry_trunc_window() -> usize {
    std::env::var("ROUND84_INPLACE_QUOTIENT_CARRY_TRUNC_W")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(21)
        .max(1)
}

fn round84_inplace_vent_carry_enabled() -> bool {
    std::env::var("ROUND84_INPLACE_VENT_CARRY")
        .ok()
        .as_deref()
        == Some("1")
}

fn round84_correction_wrap_borrow_quotient_top_enabled() -> bool {
    std::env::var("ROUND84_CORRECTION_WRAP_BORROW_QUOTIENT_TOP")
        .ok()
        .as_deref()
        == Some("1")
}

fn round84_keep_quotient_product_enabled() -> bool {
    std::env::var("ROUND84_KEEP_QUOTIENT_PRODUCT")
        .ok()
        .as_deref()
        == Some("1")
}

fn round84_qprod_naf_enabled() -> bool {
    std::env::var("R84_QPROD_NAF").ok().as_deref() == Some("1")
}

fn round84_qprod_short_enabled() -> bool {
    std::env::var("ROUND84_QPROD_SHORT").ok().as_deref() == Some("1")
}

struct Round84FoldStep {
    shift: usize,
    add: bool,
    wrap: QubitId,
}

struct Round84AggregateFold {
    steps: Vec<Round84FoldStep>,
    quotient: Vec<QubitId>,
    correction_wrap: QubitId,
    correction_wrap_owned: bool,
    product: Option<Vec<QubitId>>,
}

fn round84_update_fold_quotient(
    b: &mut B,
    quotient: &[QubitId],
    hi: &[QubitId],
    step: &Round84FoldStep,
    inverse: bool,
) {
    let add = step.add != inverse;
    let update_wrap = |b: &mut B| {
        if add {
            cadd_nbit_const_direct_fast(b, quotient, U256::from(1), step.wrap);
        } else {
            csub_nbit_const_direct_fast(b, quotient, U256::from(1), step.wrap);
        }
    };
    let update_spill = |b: &mut B| {
        if step.shift == 0 {
            return;
        }
        let pad = b.alloc_qubits(quotient.len() - step.shift);
        let mut spill = hi[hi.len() - step.shift..].to_vec();
        spill.extend_from_slice(&pad);
        if add {
            add_nbit_qq(b, &spill, quotient);
        } else {
            sub_nbit_qq(b, &spill, quotient);
        }
        b.free_vec(&pad);
    };

    if inverse {
        update_wrap(b);
        update_spill(b);
    } else {
        update_spill(b);
        update_wrap(b);
    }
}

fn round84_compute_quotient_c_product(b: &mut B, quotient: &[QubitId], dirty: &[QubitId]) -> Vec<QubitId> {

    let q = &quotient[..33];
    let product = b.alloc_qubits(66);
    for i in 0..q.len() {
        b.cx(q[i], product[i]);
    }
    if round84_qprod_naf_enabled() {
        for (shift, add) in [(10usize, true), (32, true), (5, false), (4, false)] {
            if round84_qprod_vent_pad_enabled()
                && (product.len() - shift - q.len()) >= round84_qprod_vent_pad_min_width()
            {
                round84_qprod_shifted_addsub_vented(b, q, &product, shift, add, dirty);
                continue;
            }
            let target = &product[shift..];
            if round84_qprod_short_enabled() {
                if add {
                    add_short_to_long_qq_fast_no_cin(b, q, target);
                } else {
                    sub_short_to_long_qq_fast_no_cin(b, q, target);
                }
            } else {
                let pad = b.alloc_qubits(target.len() - q.len());
                let mut source = q.to_vec();
                source.extend_from_slice(&pad);
                if add {
                    round84_add_small(b, &source, target);
                } else {
                    round84_sub_small(b, &source, target);
                }
                b.free_vec(&pad);
            }
        }
    } else {
        for shift in [4usize, 6, 7, 8, 9, 32] {
            let target = &product[shift..];
            if round84_qprod_short_enabled() {
                add_short_to_long_qq_fast_no_cin(b, q, target);
            } else {
                let pad = b.alloc_qubits(target.len() - q.len());
                let mut source = q.to_vec();
                source.extend_from_slice(&pad);
                round84_add_small(b, &source, target);
                b.free_vec(&pad);
            }
        }
    }
    product
}

fn round84_uncompute_quotient_c_product(b: &mut B, quotient: &[QubitId], product: &[QubitId], dirty: &[QubitId]) {
    let q = &quotient[..33];
    if round84_qprod_naf_enabled() {
        for (shift, add) in [(10usize, true), (32, true), (5, false), (4, false)]
            .into_iter()
            .rev()
        {
            if round84_qprod_vent_pad_enabled()
                && (product.len() - shift - q.len()) >= round84_qprod_vent_pad_min_width()
            {

                round84_qprod_shifted_addsub_vented(b, q, product, shift, !add, dirty);
                continue;
            }
            let target = &product[shift..];
            if round84_qprod_short_enabled() {
                if add {
                    sub_short_to_long_qq_fast_no_cin(b, q, target);
                } else {
                    add_short_to_long_qq_fast_no_cin(b, q, target);
                }
            } else {
                let pad = b.alloc_qubits(target.len() - q.len());
                let mut source = q.to_vec();
                source.extend_from_slice(&pad);
                if add {
                    round84_sub_small(b, &source, target);
                } else {
                    round84_add_small(b, &source, target);
                }
                b.free_vec(&pad);
            }
        }
    } else {
        for shift in [4usize, 6, 7, 8, 9, 32].into_iter().rev() {
            let target = &product[shift..];
            if round84_qprod_short_enabled() {
                sub_short_to_long_qq_fast_no_cin(b, q, target);
            } else {
                let pad = b.alloc_qubits(target.len() - q.len());
                let mut source = q.to_vec();
                source.extend_from_slice(&pad);
                round84_sub_small(b, &source, target);
                b.free_vec(&pad);
            }
        }
    }
    for i in 0..q.len() {
        b.cx(q[i], product[i]);
    }
    b.free_vec(product);
}

fn round84_add_narrow_correction(
    b: &mut B,
    lo: &[QubitId],
    product: &[QubitId],
    dirty: &[QubitId],
    borrowed_wrap: Option<QubitId>,
) -> (QubitId, bool) {
    let (wrap, owned_wrap) = borrowed_wrap.map_or_else(|| (b.alloc_qubit(), true), |q| (q, false));
    let source_top = b.alloc_qubit();
    let mut target_ext = lo[..product.len()].to_vec();
    target_ext.push(wrap);
    let mut source_ext = product.to_vec();
    source_ext.push(source_top);
    round84_add_small(b, &source_ext, &target_ext);
    b.free(source_top);
    let high = &lo[product.len()..];
    if round84_inplace_vent_carry_enabled() {
        let clean2 = [b.alloc_qubit(), b.alloc_qubit()];
        venting::ciadd_dirty_2clean_classical(
            b,
            high,
            &dirty[..high.len() - 2],
            &clean2,
            1,
            wrap,
            false,
        );
        b.free(clean2[0]);
        b.free(clean2[1]);
    } else {
        cadd_nbit_const_direct_trunc_fast(
            b,
            high,
            U256::from(1),
            wrap,
            round84_inplace_quotient_carry_trunc_window(),
        );
    }
    (wrap, owned_wrap)
}

fn round84_sub_narrow_correction(
    b: &mut B,
    lo: &[QubitId],
    product: &[QubitId],
    wrap: QubitId,
    dirty: &[QubitId],
    owned_wrap: bool,
) {
    let high = &lo[product.len()..];
    if round84_inplace_vent_carry_enabled() {
        let clean2 = [b.alloc_qubit(), b.alloc_qubit()];
        venting::cisub_dirty_2clean_classical(b, high, &dirty[..high.len() - 2], &clean2, 1, wrap);
        b.free(clean2[0]);
        b.free(clean2[1]);
    } else {
        csub_nbit_const_direct_trunc_fast(
            b,
            high,
            U256::from(1),
            wrap,
            round84_inplace_quotient_carry_trunc_window(),
        );
    }
    let source_top = b.alloc_qubit();
    let mut target_ext = lo[..product.len()].to_vec();
    target_ext.push(wrap);
    let mut source_ext = product.to_vec();
    source_ext.push(source_top);
    round84_sub_small(b, &source_ext, &target_ext);
    b.free(source_top);
    if owned_wrap {
        b.free(wrap);
    }
}

fn round84_fold_hi_into_lo_aggregate(
    b: &mut B,
    lo: &[QubitId],
    hi: &[QubitId],
    dirty: &[QubitId],
) -> Round84AggregateFold {
    let n = lo.len();
    let quotient = b.alloc_qubits(34);

    let terms = [
        (0usize, true),
        (4, false),
        (5, false),
        (10, true),
        (32, true),
    ];
    let mut steps = Vec::with_capacity(terms.len());

    for (shift, add) in terms {
        let width = n - shift;
        let wrap = b.alloc_qubit();
        let source_top = b.alloc_qubit();
        let mut target_ext = lo[shift..].to_vec();
        target_ext.push(wrap);
        let mut source_ext = hi[..width].to_vec();
        source_ext.push(source_top);
        if add {
            add_nbit_qq(b, &source_ext, &target_ext);
        } else {
            sub_nbit_qq(b, &source_ext, &target_ext);
        }
        b.free(source_top);

        let step = Round84FoldStep { shift, add, wrap };
        round84_update_fold_quotient(b, &quotient, hi, &step, false);
        steps.push(step);
    }

    let product = round84_compute_quotient_c_product(b, &quotient, dirty);
    let borrowed_correction_wrap = round84_correction_wrap_borrow_quotient_top_enabled()
        .then_some(quotient[33]);
    let (correction_wrap, correction_wrap_owned) =
        round84_add_narrow_correction(b, lo, &product, dirty, borrowed_correction_wrap);
    let product = if round84_keep_quotient_product_enabled() {
        Some(product)
    } else {
        round84_uncompute_quotient_c_product(b, &quotient, &product, dirty);
        None
    };
    Round84AggregateFold {
        steps,
        quotient,
        correction_wrap,
        correction_wrap_owned,
        product,
    }
}

fn round84_unfold_hi_from_lo_aggregate(
    b: &mut B,
    lo: &[QubitId],
    hi: &[QubitId],
    dirty: &[QubitId],
    state: Round84AggregateFold,
) {
    let product = state
        .product
        .unwrap_or_else(|| round84_compute_quotient_c_product(b, &state.quotient, dirty));
    round84_sub_narrow_correction(
        b,
        lo,
        &product,
        state.correction_wrap,
        dirty,
        state.correction_wrap_owned,
    );
    round84_uncompute_quotient_c_product(b, &state.quotient, &product, dirty);

    for step in state.steps.into_iter().rev() {
        round84_update_fold_quotient(b, &state.quotient, hi, &step, true);
        let width = lo.len() - step.shift;
        let source_top = b.alloc_qubit();
        let mut target_ext = lo[step.shift..].to_vec();
        target_ext.push(step.wrap);
        let mut source_ext = hi[..width].to_vec();
        source_ext.push(source_top);
        if step.add {
            sub_nbit_qq(b, &source_ext, &target_ext);
        } else {
            add_nbit_qq(b, &source_ext, &target_ext);
        }
        b.free(source_top);
        b.free(step.wrap);
    }
    b.free_vec(&state.quotient);
}

pub(crate) fn squaring_sub_from_acc_schoolbook(b: &mut B, acc: &[QubitId], x: &[QubitId], p: U256) {
    let n = acc.len();
    debug_assert_eq!(n, 256);
    debug_assert_eq!(x.len(), n);
    let c = U256::MAX.wrapping_sub(p).wrapping_add(U256::from(1));

    let tmp_ext = b.alloc_qubits(2 * n);

    schoolbook_square_symmetric(b, x, &tmp_ext);

    let lo: Vec<QubitId> = tmp_ext[0..n].to_vec();
    let hi: Vec<QubitId> = tmp_ext[n..2 * n].to_vec();
    mod_sub_qq_fast(b, acc, &lo, p);
    let _ = c;

    mod_sub_qq_fast(b, acc, &hi, p);
    for _ in 0..4 {
        mod_double_inplace_fast(b, &hi, p);
    }
    mod_sub_qq_fast(b, acc, &hi, p);
    for _ in 0..2 {
        mod_double_inplace_fast(b, &hi, p);
    }
    mod_add_qq_fast(b, acc, &hi, p);
    for _ in 0..4 {
        mod_double_inplace_fast(b, &hi, p);
    }
    mod_sub_qq_fast(b, acc, &hi, p);
    let (spill, flag_inv, ovf) = mod_shift_left_by_k(b, &hi, p, 22);
    mod_sub_qq(b, acc, &hi, p);
    mod_shift_right_by_k(b, &hi, p, 22, spill, flag_inv, ovf);
    for _ in 0..10 {
        mod_halve_inplace_fast(b, &hi, p);
    }

    schoolbook_square_symmetric_inverse(b, x, &tmp_ext);

    b.free_vec(&tmp_ext);
}

pub(crate) fn squaring_sub_from_acc_karatsuba(b: &mut B, acc: &[QubitId], x: &[QubitId], p: U256) {
    let n = acc.len();
    debug_assert_eq!(n, 256);
    debug_assert_eq!(x.len(), n);
    let h = n / 2;
    let x_lo: Vec<QubitId> = x[0..h].to_vec();
    let x_hi: Vec<QubitId> = x[h..n].to_vec();

    let mut z1_reg = b.alloc_qubits(2 * (h + 1));

    let free_z1_top = std::env::var("KARA_FREE_Z1_TOPBIT").ok().as_deref() == Some("1");

    let z02_lowq = std::env::var("KARA_Z02_LOWQ").ok().as_deref() == Some("1");

    {
        let x_sum = b.alloc_qubits(h + 1);
        karatsuba_half_sum_compute(b, &x_lo, &x_hi, &x_sum);
        schoolbook_square_symmetric(b, &x_sum, &z1_reg);
        karatsuba_half_sum_uncompute(b, &x_lo, &x_hi, &x_sum);
        b.free_vec(&x_sum);
    }

    let tmp_ext = b.alloc_qubits(2 * n);

    {
        let slice: Vec<QubitId> = tmp_ext[0..2 * h].to_vec();
        if z02_lowq {

            let host: Vec<QubitId> = tmp_ext[2 * h..4 * h].to_vec();
            schoolbook_square_symmetric_hosted(b, &x_lo, &slice, &host);
        } else {
            schoolbook_square_symmetric(b, &x_lo, &slice);
        }
    }
    {
        let slice: Vec<QubitId> = tmp_ext[2 * h..4 * h].to_vec();
        if z02_lowq {
            if kara_z2_selfhost_enabled() {
                if square_selfhost_safe_lane_reuse_enabled() {

                    let clean_square_bits = [z1_reg[1], tmp_ext[1]];
                    schoolbook_square_symmetric_lowq_selfhosted_with_clean_supplement(
                        b,
                        &x_hi,
                        &slice,
                        &clean_square_bits,
                    );
                } else {
                    schoolbook_square_symmetric_lowq_selfhosted(b, &x_hi, &slice);
                }
            } else {
                schoolbook_square_symmetric_lowq(b, &x_hi, &slice);
            }
        } else {
            schoolbook_square_symmetric(b, &x_hi, &slice);
        }
    }

    {
        let pad = b.alloc_qubits(2);
        let mut z0_ext: Vec<QubitId> = tmp_ext[0..2 * h].to_vec();
        z0_ext.extend_from_slice(&pad);
        sub_nbit_qq(b, &z0_ext, &z1_reg);
        b.free_vec(&pad);
    }
    {
        let pad = b.alloc_qubits(2);
        let mut z2_ext: Vec<QubitId> = tmp_ext[2 * h..4 * h].to_vec();
        z2_ext.extend_from_slice(&pad);
        sub_nbit_qq(b, &z2_ext, &z1_reg);
        b.free_vec(&pad);
    }

    if free_z1_top {
        let top = z1_reg.pop().expect("z1_reg width 2*(h+1) >= 2");
        b.free(top);
    }
    {
        let pad = b.alloc_qubits(3 * h - z1_reg.len());
        let mut z1_ext: Vec<QubitId> = z1_reg.to_vec();
        z1_ext.extend_from_slice(&pad);
        let acc_slice: Vec<QubitId> = tmp_ext[h..4 * h].to_vec();
        add_nbit_qq(b, &z1_ext, &acc_slice);
        b.free_vec(&pad);
    }

    let mod_fast = std::env::var("KARA_SOL_MOD_FAST").ok().as_deref() == Some("1");
    let dbl_fast = std::env::var("KARA_SOL_DBL_FAST").ok().as_deref() == Some("1");
    let shift_fast = std::env::var("KARA_SOL_SHIFT_FAST").ok().as_deref() == Some("1");
    let lo: Vec<QubitId> = tmp_ext[0..n].to_vec();
    let hi: Vec<QubitId> = tmp_ext[n..2 * n].to_vec();

    let mod_vent = std::env::var("KARA_SOL_MOD_VENT").ok().as_deref() == Some("1");
    let mod_sub = |b: &mut B, acc: &[QubitId], a: &[QubitId]| {
        if mod_vent {
            mod_sub_qq_vent(b, acc, a, p);
        } else if mod_fast {
            mod_sub_qq_fast(b, acc, a, p);
        } else {
            mod_sub_qq(b, acc, a, p);
        }
    };
    let mod_add = |b: &mut B, acc: &[QubitId], a: &[QubitId]| {
        if mod_vent {
            mod_add_qq_vent(b, acc, a, p);
        } else if mod_fast {
            mod_add_qq_fast(b, acc, a, p);
        } else {
            mod_add_qq(b, acc, a, p);
        }
    };
    let mod_dbl = |b: &mut B, v: &[QubitId]| {
        if dbl_fast {
            mod_double_inplace_fast(b, v, p);
        } else {
            mod_double_inplace_direct_const_fast(b, v, p);
        }
    };
    let mod_hlv = |b: &mut B, v: &[QubitId]| {
        if dbl_fast {
            mod_halve_inplace_fast(b, v, p);
        } else {
            mod_halve_inplace_direct_const_fast(b, v, p);
        }
    };
    b.set_phase("r84k_sol_subadd");
    mod_sub(b, acc, &lo);
    mod_sub(b, acc, &hi);
    for _ in 0..4 {
        mod_dbl(b, &hi);
    }
    mod_sub(b, acc, &hi);
    for _ in 0..2 {
        mod_dbl(b, &hi);
    }
    mod_add(b, acc, &hi);
    for _ in 0..4 {
        mod_dbl(b, &hi);
    }
    mod_sub(b, acc, &hi);
    b.set_phase("r84k_sol_shift");

    let shift_dirty = std::env::var("ROUND84_XTAIL_BORROW_CARRIES")
        .ok()
        .as_deref()
        == Some("1");
    if shift_dirty {

        b.set_phase("r84k_sol_dbl22");
        for _ in 0..22 {
            mod_dbl(b, &hi);
        }
        b.set_phase("r84k_sol_midsub");
        mod_sub(b, acc, &hi);
        b.set_phase("r84k_sol_hlv22");
        for _ in 0..22 {
            mod_hlv(b, &hi);
        }
    } else {
        b.set_phase("r84k_sol_shiftL");
        let (spill, flag_inv, ovf) = if shift_fast {
            mod_shift_left_by_k(b, &hi, p, 22)
        } else {
            mod_shift_left_by_k_lowq(b, &hi, p, 22)
        };
        b.set_phase("r84k_sol_midsub");
        mod_sub(b, acc, &hi);
        b.set_phase("r84k_sol_shiftR");
        if shift_fast {
            mod_shift_right_by_k(b, &hi, p, 22, spill, flag_inv, ovf);
        } else {
            mod_shift_right_by_k_lowq(b, &hi, p, 22, spill, flag_inv, ovf);
        }
    }
    b.set_phase("r84k_sol_halve");
    for _ in 0..10 {
        mod_hlv(b, &hi);
    }

    b.set_phase("r84k_inv_combine");
    {
        let pad = b.alloc_qubits(3 * h - z1_reg.len());
        let mut z1_ext: Vec<QubitId> = z1_reg.to_vec();
        z1_ext.extend_from_slice(&pad);
        let acc_slice: Vec<QubitId> = tmp_ext[h..4 * h].to_vec();
        sub_nbit_qq(b, &z1_ext, &acc_slice);
        b.free_vec(&pad);
    }

    if free_z1_top {
        let top = b.alloc_qubit();
        z1_reg.push(top);
    }
    {
        let pad = b.alloc_qubits(2);
        let mut z2_ext: Vec<QubitId> = tmp_ext[2 * h..4 * h].to_vec();
        z2_ext.extend_from_slice(&pad);
        add_nbit_qq(b, &z2_ext, &z1_reg);
        b.free_vec(&pad);
    }
    {
        let pad = b.alloc_qubits(2);
        let mut z0_ext: Vec<QubitId> = tmp_ext[0..2 * h].to_vec();
        z0_ext.extend_from_slice(&pad);
        add_nbit_qq(b, &z0_ext, &z1_reg);
        b.free_vec(&pad);
    }

    b.set_phase("r84k_z_inv_squares");
    {
        let slice: Vec<QubitId> = tmp_ext[2 * h..4 * h].to_vec();
        if z02_lowq {
            if kara_z2_selfhost_enabled() {
                if square_selfhost_safe_lane_reuse_enabled() {

                    let clean_square_bits = [z1_reg[1], tmp_ext[1]];
                    schoolbook_square_symmetric_lowq_selfhosted_inverse_with_clean_supplement(
                        b,
                        &x_hi,
                        &slice,
                        &clean_square_bits,
                    );
                } else {
                    schoolbook_square_symmetric_lowq_selfhosted_inverse(b, &x_hi, &slice);
                }
            } else {
                schoolbook_square_symmetric_lowq_inverse(b, &x_hi, &slice);
            }
        } else {
            schoolbook_square_symmetric_inverse(b, &x_hi, &slice);
        }
    }
    {
        let slice: Vec<QubitId> = tmp_ext[0..2 * h].to_vec();
        if z02_lowq {

            let host: Vec<QubitId> = tmp_ext[2 * h..4 * h].to_vec();
            schoolbook_square_symmetric_hosted_inverse(b, &x_lo, &slice, &host);
        } else {
            schoolbook_square_symmetric_inverse(b, &x_lo, &slice);
        }
    }
    b.free_vec(&tmp_ext);

    {
        let x_sum = b.alloc_qubits(h + 1);
        karatsuba_half_sum_compute(b, &x_lo, &x_hi, &x_sum);
        schoolbook_square_symmetric_inverse(b, &x_sum, &z1_reg);
        karatsuba_half_sum_uncompute(b, &x_lo, &x_hi, &x_sum);
        b.free_vec(&x_sum);
    }

    b.free_vec(&z1_reg);
}

pub(crate) fn squaring_sub_from_acc_schoolbook_lowq_shift22(
    b: &mut B,
    acc: &[QubitId],
    x: &[QubitId],
    p: U256,
) {
    let n = acc.len();
    debug_assert_eq!(n, 256);
    debug_assert_eq!(x.len(), n);

    let tmp_ext = b.alloc_qubits(2 * n);
    b.set_phase("round84_inplace_solinas_square_forward");
    if xtail_sq_selfhost_enabled() {
        schoolbook_square_symmetric_lowq_selfhosted(b, x, &tmp_ext);
    } else {
        schoolbook_square_symmetric_lowq(b, x, &tmp_ext);
    }

    let lo: Vec<QubitId> = tmp_ext[0..n].to_vec();
    let hi: Vec<QubitId> = tmp_ext[n..2 * n].to_vec();
    if round84_inplace_solinas_fold_enabled() {
        b.set_phase("round84_inplace_solinas_fold");
        let state = round84_fold_hi_into_lo_aggregate(b, &lo, &hi, acc);
        b.set_phase("round84_inplace_solinas_sub");
        mod_sub_qq_vent(b, acc, &lo, p);
        b.set_phase("round84_inplace_solinas_unfold");
        round84_unfold_hi_from_lo_aggregate(b, &lo, &hi, acc, state);
    } else {
        mod_sub_qq(b, acc, &lo, p);
        mod_sub_qq(b, acc, &hi, p);
        for _ in 0..4 {
            mod_double_inplace_direct_const_fast(b, &hi, p);
        }
        mod_sub_qq(b, acc, &hi, p);
        for _ in 0..2 {
            mod_double_inplace_direct_const_fast(b, &hi, p);
        }
        mod_add_qq(b, acc, &hi, p);
        for _ in 0..4 {
            mod_double_inplace_direct_const_fast(b, &hi, p);
        }
        mod_sub_qq(b, acc, &hi, p);
        let (spill, flag_inv, ovf) = mod_shift_left_by_k_lowq(b, &hi, p, 22);
        if r84_lowq_enabled() {
            mod_sub_qq_lowq(b, acc, &hi, p);
        } else {
            mod_sub_qq(b, acc, &hi, p);
        }
        mod_shift_right_by_k_lowq(b, &hi, p, 22, spill, flag_inv, ovf);
        for _ in 0..10 {
            mod_halve_inplace_direct_const_fast(b, &hi, p);
        }
    }

    b.set_phase("round84_inplace_solinas_square_inverse");
    if xtail_sq_selfhost_enabled() {
        schoolbook_square_symmetric_lowq_selfhosted_inverse(b, x, &tmp_ext);
    } else {
        schoolbook_square_symmetric_lowq_inverse(b, x, &tmp_ext);
    }
    b.free_vec(&tmp_ext);
}

pub(crate) fn squaring_sub_from_acc_walk_controls_lowq(b: &mut B, acc: &[QubitId], x: &[QubitId], p: U256) {
    let n = acc.len();
    debug_assert_eq!(n, 256);
    debug_assert_eq!(x.len(), n);

    let ctrl_copy = b.alloc_qubits(n);
    for i in 0..n {
        b.cx(x[i], ctrl_copy[i]);
    }

    mod_neg_inplace_fast(b, x, p);
    for i in 0..n {
        cmod_add_qq(b, acc, x, ctrl_copy[i], p);
        if i < n - 1 {
            mod_double_inplace_fast(b, x, p);
        }
    }
    for _ in 0..(n - 1) {
        mod_halve_inplace_fast(b, x, p);
    }
    mod_neg_inplace_fast(b, x, p);

    for i in 0..n {
        b.cx(x[i], ctrl_copy[i]);
    }
    b.free_vec(&ctrl_copy);
}

fn round84_qprod_vent_pad_enabled() -> bool {
    std::env::var("ROUND84_QPROD_VENT_PAD").ok().as_deref() == Some("1")
}

fn round84_qprod_vent_pad_min_width() -> usize {
    std::env::var("ROUND84_QPROD_VENT_PAD_MINW")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(29)
}

fn round84_qprod_shifted_addsub_vented(
    b: &mut B,
    q: &[QubitId],
    product: &[QubitId],
    shift: usize,
    add: bool,
    dirty: &[QubitId],
) {
    let m = q.len();
    let total = product.len() - shift;
    debug_assert!(total >= m);
    let high_w = total - m;

    if high_w < 5 || dirty.len() < high_w.saturating_sub(2) {
        let target = &product[shift..];
        let pad = b.alloc_qubits(high_w);
        let mut source = q.to_vec();
        source.extend_from_slice(&pad);
        if add {
            round84_add_small(b, &source, target);
        } else {
            round84_sub_small(b, &source, target);
        }
        b.free_vec(&pad);
        return;
    }

    let wrap = b.alloc_qubit();
    let mut low_ext = product[shift..shift + m].to_vec();
    low_ext.push(wrap);
    let high = &product[shift + m..];

    if add {

        let c_in = b.alloc_qubit();
        cuccaro_add_low_to_ext_clean(b, q, &low_ext, c_in);
        b.free(c_in);

        let clean2 = [b.alloc_qubit(), b.alloc_qubit()];
        venting::ciadd_dirty_2clean_classical(
            b,
            high,
            &dirty[..high_w - 2],
            &clean2,
            1,
            wrap,
            false,
        );
        b.free(clean2[1]);
        b.free(clean2[0]);

        cmp_lt_into(b, &product[shift..shift + m], q, wrap);
    } else {

        let c_in = b.alloc_qubit();
        cuccaro_sub_low_to_ext_clean(b, q, &low_ext, c_in);
        b.free(c_in);

        let clean2 = [b.alloc_qubit(), b.alloc_qubit()];
        venting::cisub_dirty_2clean_classical(b, high, &dirty[..high_w - 2], &clean2, 1, wrap);
        b.free(clean2[1]);
        b.free(clean2[0]);

        for &qb in q {
            b.x(qb);
        }
        cmp_lt_into(b, q, &product[shift..shift + m], wrap);
        for &qb in q {
            b.x(qb);
        }
    }
    b.free(wrap);
}

pub(crate) fn cross_addsub_stage1(
    b: &mut B,
    prod: &[QubitId],
    off: usize,
    addend: &[QubitId],
    ctrl: QubitId,
) {
    let w = prod.len();
    let m = addend.len();
    if m == 0 {
        return;
    }
    debug_assert!(off + m <= w);
    let t = b.alloc_qubits(w);
    for k in 0..m {
        b.cx(addend[k], t[off + k]);
    }
    b.x(ctrl);
    for k in 0..w {
        b.cx(ctrl, t[k]);
    }
    let cin = b.alloc_qubit();
    b.cx(ctrl, cin);
    cuccaro_add(b, &t, prod, cin);
    b.cx(ctrl, cin);
    for k in 0..w {
        b.cx(ctrl, t[k]);
    }
    b.x(ctrl);
    for k in 0..m {
        b.cx(addend[k], t[off + k]);
    }
    b.free(cin);
    b.free_vec(&t);
}

pub(crate) fn square_addsub_stage1(b: &mut B, x: &[QubitId], prod: &[QubitId]) {
    let n = x.len();
    debug_assert_eq!(prod.len(), 2 * n);

    for i in 0..n {
        let m = n - 1 - i;
        if m == 0 {
            continue;
        }
        let off = 2 * i + 1;
        cross_addsub_stage1(b, prod, off, &x[i + 1..n], x[i]);
    }

    let t = b.alloc_qubits(2 * n);
    for i in 0..n {
        let p = 2 * i + 1;
        if p < 2 * n {
            b.cx(x[i], t[p]);
        }
    }
    let cin = b.alloc_qubit();
    cuccaro_add(b, &t, prod, cin);
    b.free(cin);
    for i in 0..n {
        let p = 2 * i + 1;
        if p < 2 * n {
            b.cx(x[i], t[p]);
        }
    }
    b.free_vec(&t);

    let zero_ctrl = b.alloc_qubit();
    cross_addsub_stage1(b, prod, 0, x, zero_ctrl);
    b.free(zero_ctrl);
}

fn square_corr_forward(b: &mut B, x: &[QubitId], prod: &[QubitId]) {
    let n = x.len();

    let zeros = b.alloc_qubits(n);
    let mut a2d: Vec<QubitId> = Vec::with_capacity(2 * n);
    for i in 0..n {
        a2d.push(zeros[i]);
        a2d.push(x[i]);
    }
    let cin = b.alloc_qubit();
    cuccaro_add(b, &a2d, prod, cin);
    b.free(cin);
    b.free_vec(&zeros);

    let pad = b.alloc_qubits(n);
    let mut xext = x.to_vec();
    xext.extend_from_slice(&pad);
    let cinx = b.alloc_qubit();
    cuccaro_sub(b, &xext, prod, cinx);
    b.free(cinx);
    b.free_vec(&pad);

    let p1 = b.alloc_qubit();
    let mut a = x[0..n - 1].to_vec();
    a.push(p1);
    let high: Vec<QubitId> = prod[n..2 * n].to_vec();
    let cinl = b.alloc_qubit();
    b.x(cinl);
    cuccaro_add(b, &a, &high, cinl);
    b.x(cinl);
    b.free(cinl);
    b.free(p1);

    b.x(prod[2 * n - 1]);
}

fn square_corr_inverse(b: &mut B, x: &[QubitId], prod: &[QubitId]) {
    let n = x.len();

    b.x(prod[2 * n - 1]);

    let p1 = b.alloc_qubit();
    let mut a = x[0..n - 1].to_vec();
    a.push(p1);
    let high: Vec<QubitId> = prod[n..2 * n].to_vec();
    let cinl = b.alloc_qubit();
    b.x(cinl);
    cuccaro_sub(b, &a, &high, cinl);
    b.x(cinl);
    b.free(cinl);
    b.free(p1);

    let pad = b.alloc_qubits(n);
    let mut xext = x.to_vec();
    xext.extend_from_slice(&pad);
    let cinx = b.alloc_qubit();
    cuccaro_add(b, &xext, prod, cinx);
    b.free(cinx);
    b.free_vec(&pad);

    let zeros = b.alloc_qubits(n);
    let mut a2d: Vec<QubitId> = Vec::with_capacity(2 * n);
    for i in 0..n {
        a2d.push(zeros[i]);
        a2d.push(x[i]);
    }
    let cin = b.alloc_qubit();
    cuccaro_sub(b, &a2d, prod, cin);
    b.free(cin);
    b.free_vec(&zeros);
}

pub(crate) fn square_addsub_local(b: &mut B, x: &[QubitId], prod: &[QubitId]) {
    let n = x.len();
    debug_assert_eq!(prod.len(), 2 * n);

    for i in 0..n {
        let m = n - 1 - i;
        if m == 0 {
            continue;
        }
        let off = 2 * i + 1;
        let slice: Vec<QubitId> = prod[off..off + m + 1].to_vec();
        controlled_add_subtract_lowq(b, &x[i + 1..n], &slice, x[i]);
    }

    let t = b.alloc_qubits(2 * n);
    for i in 0..n {
        let p = 2 * i + 1;
        if p < 2 * n {
            b.cx(x[i], t[p]);
        }
    }
    let cin = b.alloc_qubit();
    cuccaro_add(b, &t, prod, cin);
    b.free(cin);
    for i in 0..n {
        let p = 2 * i + 1;
        if p < 2 * n {
            b.cx(x[i], t[p]);
        }
    }
    b.free_vec(&t);

    let zc = b.alloc_qubit();
    cross_addsub_stage1(b, prod, 0, x, zc);
    b.free(zc);

    if n >= 2 {
        let oc = b.alloc_qubit();
        b.x(oc);
        cross_addsub_stage1(b, prod, n, &x[0..n - 1], oc);
        b.x(oc);
        b.free(oc);
    }

    let t2 = b.alloc_qubits(2 * n);
    b.x(t2[2 * n - 1]);
    b.x(t2[n]);
    let cin2 = b.alloc_qubit();
    cuccaro_add(b, &t2, prod, cin2);
    b.free(cin2);
    b.x(t2[2 * n - 1]);
    b.x(t2[n]);
    b.free_vec(&t2);
}

pub(crate) fn controlled_add_subtract_vented_borrowed(
    b: &mut B,
    x: &[QubitId],
    acc: &[QubitId],
    ctrl: QubitId,
    carries: &[QubitId],
) {
    let n = x.len();
    debug_assert_eq!(acc.len(), n + 1);
    let pad = b.alloc_qubit();
    let mut x_ext = x.to_vec();
    x_ext.push(pad);
    let c_in = b.alloc_qubit();
    b.x(ctrl);
    for i in 0..n {
        b.cx(ctrl, x_ext[i]);
    }
    b.cx(ctrl, c_in);
    cuccaro_add_fast_borrowed_carries(b, &x_ext, acc, c_in, carries);
    b.cx(ctrl, c_in);
    for i in 0..n {
        b.cx(ctrl, x_ext[i]);
    }
    b.x(ctrl);
    b.free(c_in);
    b.free(pad);
}

pub(crate) fn square_addsub_vented(b: &mut B, x: &[QubitId], prod: &[QubitId]) {
    let n = x.len();
    debug_assert_eq!(prod.len(), 2 * n);
    for i in 0..n {
        let m = n - 1 - i;
        if m == 0 {
            continue;
        }
        let off = 2 * i + 1;
        let slice: Vec<QubitId> = prod[off..off + m + 1].to_vec();

        let hi = off + m + 1;
        let need = m;
        let carries: Vec<QubitId> = prod[hi..hi + need].to_vec();
        controlled_add_subtract_vented_borrowed(b, &x[i + 1..n], &slice, x[i], &carries);
    }

    square_corr_forward(b, x, prod);
}

pub(crate) fn square_addsub_local_inverse(b: &mut B, x: &[QubitId], prod: &[QubitId]) {
    let n = x.len();
    debug_assert_eq!(prod.len(), 2 * n);
    let t2 = b.alloc_qubits(2 * n);
    b.x(t2[2 * n - 1]);
    b.x(t2[n]);
    let cin2 = b.alloc_qubit();
    cuccaro_sub(b, &t2, prod, cin2);
    b.free(cin2);
    b.x(t2[2 * n - 1]);
    b.x(t2[n]);
    b.free_vec(&t2);
    if n >= 2 {
        let zc = b.alloc_qubit();
        cross_addsub_stage1(b, prod, n, &x[0..n - 1], zc);
        b.free(zc);
    }
    let oc = b.alloc_qubit();
    b.x(oc);
    cross_addsub_stage1(b, prod, 0, x, oc);
    b.x(oc);
    b.free(oc);
    let t = b.alloc_qubits(2 * n);
    for i in 0..n {
        let p = 2 * i + 1;
        if p < 2 * n {
            b.cx(x[i], t[p]);
        }
    }
    let cin = b.alloc_qubit();
    cuccaro_sub(b, &t, prod, cin);
    b.free(cin);
    for i in 0..n {
        let p = 2 * i + 1;
        if p < 2 * n {
            b.cx(x[i], t[p]);
        }
    }
    b.free_vec(&t);
    for i in (0..n).rev() {
        let m = n - 1 - i;
        if m == 0 {
            continue;
        }
        let off = 2 * i + 1;
        let slice: Vec<QubitId> = prod[off..off + m + 1].to_vec();
        controlled_add_subtract_lowq_inverse(b, &x[i + 1..n], &slice, x[i]);
    }
}

pub(crate) fn controlled_add_subtract_vented_borrowed_inverse(
    b: &mut B,
    x: &[QubitId],
    acc: &[QubitId],
    ctrl: QubitId,
    carries: &[QubitId],
) {
    let n = x.len();
    debug_assert_eq!(acc.len(), n + 1);
    let pad = b.alloc_qubit();
    let mut x_ext = x.to_vec();
    x_ext.push(pad);
    let c_in = b.alloc_qubit();
    b.x(ctrl);
    for i in 0..n {
        b.cx(ctrl, x_ext[i]);
    }
    b.cx(ctrl, c_in);
    cuccaro_sub_fast_borrowed_carries(b, &x_ext, acc, c_in, carries);
    b.cx(ctrl, c_in);
    for i in 0..n {
        b.cx(ctrl, x_ext[i]);
    }
    b.x(ctrl);
    b.free(c_in);
    b.free(pad);
}

pub(crate) fn square_addsub_vented_inverse(b: &mut B, x: &[QubitId], prod: &[QubitId]) {
    let n = x.len();
    debug_assert_eq!(prod.len(), 2 * n);

    square_corr_inverse(b, x, prod);

    for i in (0..n).rev() {
        let m = n - 1 - i;
        if m == 0 {
            continue;
        }
        let off = 2 * i + 1;
        let slice: Vec<QubitId> = prod[off..off + m + 1].to_vec();
        let hi = off + m + 1;
        let carries: Vec<QubitId> = prod[hi..hi + m].to_vec();
        controlled_add_subtract_vented_borrowed_inverse(b, &x[i + 1..n], &slice, x[i], &carries);
    }
}

pub(crate) mod square_addsub_selftest {
    use super::*;
    use crate::sim::Simulator;
    use crate::circuit::OperationType;
    use sha3::digest::{ExtendableOutput, Update, XofReader};

    fn count_tof(ops: &[crate::circuit::Op]) -> usize {
        ops.iter()
            .filter(|o| matches!(o.kind, OperationType::CCX | OperationType::CCZ))
            .count()
    }

    pub(crate) fn toffoli_compare() {
        for &n in &[128usize, 129] {
            let mut b1 = B::new();
            let x1 = b1.alloc_qubits(n);
            let p1 = b1.alloc_qubits(2 * n);
            schoolbook_square_symmetric(&mut b1, &x1, &p1);
            let cur = count_tof(&b1.ops);
            let peak_cur = b1.peak_qubits;

            let mut b2 = B::new();
            let x2 = b2.alloc_qubits(n);
            let p2 = b2.alloc_qubits(2 * n);
            square_addsub_vented(&mut b2, &x2, &p2);
            let new = count_tof(&b2.ops);
            let peak_new = b2.peak_qubits;

            println!(
                "  SQ_TOF n={n}: current(AND) CCX={cur} peakQ={peak_cur} | addsub_vented CCX={new} peakQ={peak_new} | delta={}",
                cur as i64 - new as i64
            );
        }
    }

    pub(crate) fn run() {
        let big = std::env::var("TLM_SQ_SELFTEST_BIG").ok().as_deref() == Some("1");
        exhaustive_small_n();
        random_large_n(if big { 4096 } else { 64 });
        println!("  SQ_SELFTEST (vented): bit-exact vs classical x^2 — 0 divergence");
        inverse_drains_to_zero(if big { 4096 } else { 256 });
        println!("  SQ_SELFTEST (inverse): forward+inverse drains prod to 0, clean");
        toffoli_compare();
    }

    fn inverse_drains_to_zero(count: usize) {
        for &n in &[1usize, 2, 3, 8, 32, 127, 128, 129] {
            let mut seed = sha3::Shake256::default();
            seed.update(b"missed3-inv");
            seed.update(&[n as u8]);
            let mut xof = seed.finalize_xof();
            let mut buf = [0u8; 32];
            let batches = (count / 64).max(1);
            for batch in 0..batches {
                let mut xs = Vec::with_capacity(64);
                for _ in 0..64 {
                    xof.read(&mut buf);
                    let mut v = U256::from_le_bytes(buf);
                    if n < 256 {
                        v &= (U256::from(1u64) << n) - U256::from(1u64);
                    }
                    xs.push(v);
                }
                let mut b = B::new();
                let x = b.alloc_qubits(n);
                let prod = b.alloc_qubits(2 * n);
                square_addsub_vented(&mut b, &x, &prod);
                square_addsub_vented_inverse(&mut b, &x, &prod);
                let nq = b.next_qubit as usize;
                let nb = b.next_bit as usize;
                let mut s2 = sha3::Shake256::default();
                s2.update(b"missed3-inv-sim");
                s2.update(&[n as u8, batch as u8]);
                let mut xof2 = s2.finalize_xof();
                let mut sim = Simulator::new(nq, nb, &mut xof2);
                sim.clear_for_shot();
                for (shot, xv) in xs.iter().enumerate() {
                    for i in 0..n {
                        if xv.bit(i) {
                            *sim.qubit_mut(x[i]) |= 1u64 << shot;
                        }
                    }
                }
                sim.apply_iter(b.ops.iter());
                assert_eq!(sim.phase, 0, "inv n={n} b{batch}: phase garbage");

                let mut is_x = vec![false; nq];
                for &q in x.iter() {
                    is_x[q.0 as usize] = true;
                }
                for q in 0..nq {
                    if !is_x[q] {
                        assert_eq!(
                            sim.qubit(QubitId(q as u64)),
                            0,
                            "inv n={n} b{batch}: nonzero q{q} (prod/ancilla not drained)"
                        );
                    }
                }
                for (shot, xv) in xs.iter().enumerate() {
                    for i in 0..n {
                        let got = (sim.qubit(x[i]) >> shot) & 1 == 1;
                        assert_eq!(got, xv.bit(i), "inv n={n} b{batch}: x[{i}] corrupted");
                    }
                }
            }
        }
    }

    fn check_square(n: usize, xs: &[U256], label: &str) {
        assert!(xs.len() <= 64);
        let mut b = B::new();
        let x = b.alloc_qubits(n);
        let prod = b.alloc_qubits(2 * n);
        square_addsub_vented(&mut b, &x, &prod);
        let nq = b.next_qubit as usize;
        let nb = b.next_bit as usize;

        let mut seed = sha3::Shake256::default();
        seed.update(b"missed3-square-addsub-stage1");
        seed.update(label.as_bytes());
        let mut xof = seed.finalize_xof();

        let mut sim = Simulator::new(nq, nb, &mut xof);
        sim.clear_for_shot();
        for (shot, xv) in xs.iter().enumerate() {
            for i in 0..n {
                if xv.bit(i) {
                    *sim.qubit_mut(x[i]) |= 1u64 << shot;
                }
            }
        }
        sim.apply_iter(b.ops.iter());

        let cond_mask: u64 = if xs.len() == 64 {
            u64::MAX
        } else {
            (1u64 << xs.len()) - 1
        };
        if sim.phase & cond_mask != 0 {
            let mut is_reg = vec![false; nq];
            for &q in x.iter().chain(prod.iter()) {
                is_reg[q.0 as usize] = true;
            }
            for q in 0..nq {
                let v = sim.qubit(QubitId(q as u64)) & cond_mask;
                if !is_reg[q] && v != 0 {
                    eprintln!("  DIRTY q{q} = {v:#018x}");
                }
            }
            panic!("{label}: phase garbage {:#018x}", sim.phase & cond_mask);
        }

        for (shot, xv) in xs.iter().enumerate() {
            let mut out = U256::ZERO;
            for i in 0..(2 * n) {
                if (sim.qubit(prod[i]) >> shot) & 1 == 1 {
                    out |= U256::from(1u64) << i;
                }
            }
            let expect = xv.wrapping_mul(*xv);
            assert_eq!(out, expect, "{label}: shot {shot} x={xv:#x} got {out:#x}");
        }

        let mut is_reg = vec![false; nq];
        for &q in x.iter().chain(prod.iter()) {
            is_reg[q.0 as usize] = true;
        }
        for q in 0..nq {
            if !is_reg[q] {
                assert_eq!(
                    sim.qubit(QubitId(q as u64)) & cond_mask,
                    0,
                    "{label}: dirty ancilla q{q}"
                );
            }
        }
    }

    fn exhaustive_small_n() {
        for n in 1..=6usize {
            let limit = 1usize << n;
            let xs: Vec<U256> = (0..limit).map(|v| U256::from(v as u64)).collect();

            for chunk in xs.chunks(64) {
                check_square(n, chunk, &format!("exhaustive-n{n}"));
            }
        }
    }

    fn random_large_n(batches: usize) {
        for &n in &[8usize, 16, 32, 64, 127, 128, 129] {
            let mut seed = sha3::Shake256::default();
            seed.update(b"missed3-rand");
            seed.update(&[n as u8]);
            let mut xof = seed.finalize_xof();
            let mut buf = [0u8; 32];
            for batch in 0..batches {
                let mut xs = Vec::with_capacity(64);
                for _ in 0..64 {
                    xof.read(&mut buf);
                    let mut v = U256::from_le_bytes(buf);

                    if n < 256 {
                        v &= (U256::from(1u64) << n) - U256::from(1u64);
                    }
                    xs.push(v);
                }
                check_square(n, &xs, &format!("rand-n{n}-b{batch}"));
            }
        }
    }
}
