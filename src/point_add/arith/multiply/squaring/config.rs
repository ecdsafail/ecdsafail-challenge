use super::*;

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
pub(crate) fn square_row_window_min_width() -> usize {
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
pub(crate) fn square_row_max_seg() -> usize {
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

pub(crate) fn square_row_window_clean_compare_bits(
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

pub(crate) fn square_row_window_measured_carry_clear_enabled() -> bool {
    std::env::var("SQUARE_ROW_WINDOW_MEASURED_CARRY_CLEAR")
        .ok()
        .as_deref()
        == Some("1")
}
