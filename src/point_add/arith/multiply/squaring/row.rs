use super::*;

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
pub(crate) fn square_row_windowed_apply(
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
