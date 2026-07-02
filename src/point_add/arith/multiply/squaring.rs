//! Symmetric squaring: schoolbook symmetric square (+ low-peak and
//! self-hosted variants) and the windowed square-row machinery.
use super::*;

// ─── 2-level Karatsuba variants (recursive on inner half-mults) ───
// Costs 2 extra z1_inner registers of ~2*(n/4+1) qubits each (~260 total for n=256).
// Higher peak qubits; use only at low-peak mul sites.

/// Symmetric schoolbook for squaring: x² = sum_i x[i]·2^(2i) + sum_{i<j} 2·x[i]·x[j]·2^(i+j).
/// Each cross-product is computed ONCE (instead of twice in full schoolbook),
/// halving the AND count + Cuccaro_add length. Saves ~130k CCX per squaring.
///
/// Row i layout (width n-i): bit 0 = diagonal x[i] at position 2i, bit 1 = 0
/// (gap), bit k+2 = cross-product (x[i] AND x[i+1+k]) at position i+(i+1+k)+1.
#[allow(dead_code)] // retained reference/alternative impl; not on active build path
pub(crate) fn schoolbook_square_symmetric(b: &mut B, x: &[QubitId], tmp_ext: &[QubitId]) {
    let n = x.len();
    debug_assert_eq!(tmp_ext.len(), 2 * n);
    for i in 0..n {
        // Width: bit 0 = diag at pos 2i, bit 1 = gap, bits 2..(n-i) = cross-
        // products at positions 2i+2..i+n. Last bit index = n-i, so width = n-i+1.
        // Edge case: i = n-1 has only the diagonal, width = 1.
        let width = if i == n - 1 { 1 } else { n - i + 1 };
        let num_cross = if i + 1 < n { n - i - 1 } else { 0 };
        // num_cross = number of cross-products in this row = width - 2 when width >= 2.
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

#[allow(dead_code)] // retained reference/alternative impl; not on active build path
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

/// Like `schoolbook_square_symmetric` (fast, measurement UMA) but the per-row
/// Cuccaro carry lane is hosted on a caller-supplied clean register `host`
/// (returned clean) instead of a fresh allocation. Toffoli-identical to the
/// fast square, peak-identical to the lowq square — used for the z0 lobe of the
/// round84 Karatsuba square, where the not-yet-written z2 slice is clean scratch.
#[allow(dead_code)] // retained reference/alternative impl; not on active build path
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
            // The z2 sibling host is clean and disjoint from x and z0.  It has
            // ample room for both the width carry lanes and one clean c_in.
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

/// Experimental square-only reclaim.  This is deliberately opt-in: every lane
/// borrowed by the prototype is either an untouched high tail of the square
/// accumulator, a caller-proved square bit that is exactly zero, or a clean
/// sibling square destination.  Dirty-but-idle data and operand aliases are not
/// eligible.
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

/// Like `schoolbook_square_symmetric_lowq` but converts the per-row Cuccaro
/// UMA-uncompute (CCX, executed every shot) into measurement-based (fast)
/// uncompute, WITHOUT a separate clean host register. The fast carry lane is
/// hosted on the slice's OWN not-yet-written high zeros
/// (`tmp_ext[2i+width+1 ..]`, which rows 0..=i never touch) topped up with a
/// small global remainder (<=3 qubits, since the lane width exceeds the clean
/// tail by exactly the 3-bit diagonal/gap/pad overhead). Unlike
/// `schoolbook_square_symmetric_hosted` this needs no sibling clean register,
/// so it applies where the sibling slice is occupied (the Karatsuba z2 square).
/// Peak rises only by the global remainder (<=3); Toffoli drops by the whole
/// UMA-uncompute. Under `SQUARE_SELFHOST_SAFE_LANE_REUSE=1`, the source-high
/// zero is represented structurally (no allocated `pad`) and an optional
/// caller-proved clean supplement is consumed before the global remainder. The
/// borrowed carries are returned clean by the HMR uncompute.
/// Peak-bounded row window for the selfhosted square. When set (>=2), each
/// schoolbook square row's add into `tmp_ext` is sliced into this many windows;
/// the transient row register holds only one window's worth of cross-term
/// qubits at a time (peak ~= 1024 + width/windows + boundary carries) instead of
/// the full row (peak ~= 1024 + 257). Value-exact: the same product lands in
/// `tmp_ext`. Cost: a per-boundary carry-clean comparator that rebuilds the row
/// prefix (extra CCX), traded for the dropped peak qubits.
pub(crate) fn square_row_windows() -> usize {
    std::env::var("SQUARE_ROW_WINDOWS")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(0)
}

/// Minimum row width below which a row is built monolithically (windowing a
/// narrow row buys no peak but still pays the comparator tax).
fn square_row_window_min_width() -> usize {
    std::env::var("SQUARE_ROW_WINDOW_MIN_WIDTH")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(96)
}

/// When >0, each row is windowed into the *minimum* number of windows that
/// keeps every window's source segment <= this width. Rows narrow enough to fit
/// in one segment are built monolithically (no comparator tax). This minimizes
/// the carry-recovery comparator overhead: only the rows wide enough to break
/// the peak budget get windowed, and only into as many windows as needed. When
/// set, it overrides the fixed SQUARE_ROW_WINDOWS count.
fn square_row_max_seg() -> usize {
    std::env::var("SQUARE_ROW_MAX_SEG")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(0)
}

/// Optional truncation for the row-window boundary-carry cleanup comparator.
/// Default 0 means exact/full-width.  When set below the segment width, cleanup
/// compares only the high suffix of the segment and final partial sum.  This is
/// a deliberate island-hunt knob: it keeps the same low peak and saves Toffoli,
/// but wrong suffix ties leave the boundary carry dirty.
fn square_cleanup_direction(raw: &str) -> Option<bool> {
    match raw.trim().to_ascii_lowercase().as_str() {
        "f" | "forward" | "false" | "0" => Some(false),
        "r" | "reverse" | "true" | "1" => Some(true),
        _ => None,
    }
}

fn square_row_window_clean_compare_bits(row: usize, window: usize, reverse: bool) -> usize {
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

/// Set row bit `j` of square row `i` into `t`. Bit 0 = x_i (diagonal low),
/// bit 1 = 0 (gap), bit 2+k = x_i & x_{i+1+k} (doubled cross term).
fn square_row_bit_set(b: &mut B, x: &[QubitId], i: usize, j: usize, t: QubitId) {
    if j == 0 {
        b.cx(x[i], t);
    } else if j == 1 {
        // gap bit: zero, nothing to do
    } else {
        b.ccx(x[i], x[i + 1 + (j - 2)], t);
    }
}

/// Measurement-based clear of a row bit set by `square_row_bit_set` (the bit is
/// known to equal its set expression at clear time).
fn square_row_bit_clear_hmr(b: &mut B, x: &[QubitId], i: usize, j: usize, t: QubitId) {
    if j == 0 {
        b.cx(x[i], t);
    } else if j == 1 {
        // gap bit: nothing
    } else {
        let m = b.alloc_bit();
        b.hmr(t, m);
        b.cz_if(x[i], x[i + 1 + (j - 2)], m);
    }
}

/// Windowed selfhosted square row add: `tmp_ext[2i ..] += row_i` where
/// `row_i` has `width` bits, built one window at a time. `forward=true` adds,
/// `forward=false` subtracts (the inverse). Value-identical to a single
/// `cuccaro_{add,sub}` of the full row into `tmp_ext[2i..2i+width+1]`.
///
/// The full-width add is split into a chain of low-to-ext adds. Window `w`
/// covers row bits `[lo..hi)` and writes `tmp_ext[base+lo .. base+hi+1]` (the
/// extra high cell absorbs the window carry). Because windows are contiguous in
/// `tmp_ext`, window `w`'s carry lands in `tmp_ext[base+hi]`, which is the low
/// cell of window `w+1` — so the carry chains *through* `tmp_ext` with no
/// separate carry-out ancilla and no boundary comparators. The per-window
/// Cuccaro carry lane is borrowed from `tmp_ext`'s not-yet-written high zeros
/// (rows `0..=i` never touch `tmp_ext[2i+width+1 ..]`), topped up by a small
/// global remainder, so the transient overhead is only the `seg_w`-wide source
/// window. Forward order low→high; inverse must mirror it high→low.
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

    // Window boundaries over the row bit range [0, width).
    let bounds: Vec<(usize, usize)> = (0..windows)
        .map(|w| {
            let lo = (w * width) / windows;
            let hi = ((w + 1) * width) / windows;
            (lo, hi)
        })
        .filter(|&(lo, hi)| hi > lo)
        .collect();
    let nwin = bounds.len();

    // Each interior window's carry-out (forward) / borrow-out (inverse) is
    // captured in a fresh clean ancilla `cout`, fed as the carry/borrow-IN of
    // the next window so the carry ripples across the boundary into the
    // accumulator. The final window captures its carry into tmp_ext[base+width].
    // Interior couts are NOT clean after being consumed as the next c_in
    // (Cuccaro restores c_in to the carry value), so they are uncomputed by a
    // *local* width-bounded comparator that recovers the carry from the final
    // partial sum and the rebuilt source window — peak stays ~1024 + 2*seg_w.
    //
    // The inverse (sub) is built as the structural mirror of the forward (add):
    // same window order and carry chaining, add->sub, with the borrow-recovery
    // comparator X-wrapped per the Cuccaro sub convention. It SUBTRACTS the same
    // row value the forward ADDED, so tmp_ext returns to its pre-square state.

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

    // The Cuccaro carry lane for each window add/sub is borrowed from tmp_ext's
    // clean high tail (positions beyond this row's footprint base+width+1, which
    // no row 0..=i touches), so the per-window transient overhead is only the
    // seg_w source bits + a 0-pad + the cout ancilla (~seg_w+2), never an
    // allocated carry array. The interior carry-out cleanup uses the *slow*
    // (carry-array-free) comparator, so cleanup is peak-flat (+0 beyond seg).
    let row_top = base + width + 1; // first clean tmp_ext cell above the row.
    let borrow_lane = |_b: &mut B, _need: usize| -> Vec<QubitId> {
        // Always available: tmp_ext beyond row_top is clean and >= seg_w wide
        // for every window (seg_w <= width and the high tail is wide enough).
        tmp_ext[row_top..row_top + _need].to_vec()
    };

    // carry/borrow-in for window 0 is a clean zero.
    let mut carry_in = b.alloc_qubit();
    let first_carry = carry_in;
    let mut couts: Vec<(QubitId, usize, usize, QubitId, usize)> = Vec::new();
    for (wi, &(lo, hi)) in bounds.iter().enumerate() {
        let last = wi == nwin - 1;
        let seg = build_seg(b, lo, hi);
        // Build a_block = seg ++ 0pad, acc_block = tmp[lo..hi] ++ high, n = seg_w+1.
        let pad = b.alloc_qubit();
        let mut a_block = seg.clone();
        a_block.push(pad);
        let high = if last {
            // Final window: high carry lands in the (clean) tmp_ext[base+width].
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
            // nothing extra: carry already in tmp_ext[base+width].
        } else {
            couts.push((high, lo, hi, carry_in, wi));
            carry_in = high;
        }
        clear_seg(b, lo, &seg);
    }
    // Reverse sweep: clean each interior cout with a local comparator. The
    // measured-uncompute fast comparator (~n CCX) borrows its n-wide carry lane
    // from tmp_ext's clean high tail, so cleanup adds no peak qubits. Setting
    // SQUARE_ROW_WINDOW_SLOW_CMP=1 falls back to the carry-array-free slow
    // comparator (~2n CCX, also peak-flat) for cross-checking.
    let slow_cmp = std::env::var("SQUARE_ROW_WINDOW_SLOW_CMP").ok().as_deref() == Some("1");
    let measured_clear = square_row_window_measured_carry_clear_enabled();
    for &(cout, lo, hi, cin, window) in couts.iter().rev() {
        let clean_cmp_bits = square_row_window_clean_compare_bits(i, window, !forward);
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
                // carry_out = (partial_sum < seg + cin)
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
                        b,
                        &tmp_ext[base + lo..base + hi],
                        &seg,
                        cin,
                        cout,
                        &carries,
                    );
                }
            } else {
                // borrow_out = (seg + cin > partial_diff)
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
                        b,
                        &seg,
                        &tmp_ext[base + lo..base + hi],
                        cin,
                        cout,
                        &carries,
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

pub(crate) fn schoolbook_square_symmetric_lowq_selfhosted(
    b: &mut B,
    x: &[QubitId],
    tmp_ext: &[QubitId],
) {
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

#[allow(dead_code)] // retained reference/alternative impl; not on active build path
struct Round84FoldStep {
    shift: usize,
    add: bool,
    wrap: QubitId,
}

#[allow(dead_code)] // retained reference/alternative impl; not on active build path
struct Round84AggregateFold {
    steps: Vec<Round84FoldStep>,
    quotient: Vec<QubitId>,
    correction_wrap: QubitId,
    correction_wrap_owned: bool,
    product: Option<Vec<QubitId>>,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sim::Simulator;
    use sha3::{
        digest::{ExtendableOutput, Update, XofReader},
        Shake128,
    };

    fn get<R: XofReader>(sim: &Simulator<'_, R>, qs: &[QubitId], shot: usize) -> u64 {
        let mut v = 0u64;
        for (i, &q) in qs.iter().enumerate() {
            v |= ((sim.qubit(q) >> shot) & 1) << i;
        }
        v
    }

    /// `schoolbook_square_symmetric` writes `x²` (2n-bit) into a zero `tmp_ext`,
    /// leaving `x` unchanged and every internal ancilla back at |0>. Exhaustive
    /// over all n=4 inputs (x ∈ [0,16)).
    #[test]
    fn schoolbook_square_symmetric_computes_x_squared() {
        const N: usize = 4;
        let mut b = B::new();
        let x = b.alloc_qubits(N);
        let tmp = b.alloc_qubits(2 * N);
        schoolbook_square_symmetric(&mut b, &x, &tmp);
        let nq = b.next_qubit as usize;
        let nb = b.next_bit as usize;
        let inputs: std::collections::HashSet<u64> =
            x.iter().chain(tmp.iter()).map(|q| q.0).collect();

        let mut seed = Shake128::default();
        seed.update(b"sqsym-n4");
        let mut xof = seed.finalize_xof();
        let mut sim = Simulator::new(nq, nb, &mut xof);
        for shot in 0..(1 << N) {
            for (i, &q) in x.iter().enumerate() {
                if (shot >> i) & 1 != 0 {
                    *sim.qubit_mut(q) |= 1u64 << shot;
                }
            }
        }
        sim.apply_iter(b.ops.iter());
        assert_eq!(sim.phase, 0, "phase garbage");
        for shot in 0..(1 << N) {
            let xv = shot as u64;
            assert_eq!(get(&sim, &tmp, shot), xv * xv, "x={xv}: tmp != x^2");
            assert_eq!(get(&sim, &x, shot), xv, "x={xv} changed");
        }
        for q in 0..nq as u64 {
            if !inputs.contains(&q) {
                assert_eq!(sim.qubit(QubitId(q)), 0, "ancilla q{q} not clean");
            }
        }
    }
}
