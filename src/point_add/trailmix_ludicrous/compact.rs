//! Optional reset-bounded qubit-id compaction.
//!
//! The builder reuses freed qubits, but the scored qubit count is the maximum
//! id appearing in `ops.bin`. This pass recolors non-IO temp lifetimes after
//! unconditional resets/HMRs, preserving the four benchmark registers.

use crate::circuit::{analyze_ops, Op, OperationType, QubitId, NO_BIT, NO_QUBIT};
use std::cmp::Reverse;
use std::collections::BinaryHeap;

#[derive(Clone, Debug)]
struct Segment {
    q: usize,
    ordinal: usize,
    start: usize,
    end: usize,
    assigned: usize,
}

fn q_operands(op: &Op) -> [QubitId; 3] {
    [op.q_control2, op.q_control1, op.q_target]
}

fn max_qubit_id(ops: &[Op], fixed_qubits: &[QubitId]) -> usize {
    let mut max_q = 0usize;
    for op in ops {
        for q in q_operands(op) {
            if q != NO_QUBIT {
                max_q = max_q.max(q.0 as usize + 1);
            }
        }
    }
    for &q in fixed_qubits {
        max_q = max_q.max(q.0 as usize + 1);
    }
    max_q
}

fn is_unconditional_reset(op: &Op, condition_depth: usize) -> bool {
    matches!(op.kind, OperationType::R | OperationType::Hmr)
        && op.q_target != NO_QUBIT
        && op.c_condition == NO_BIT
        && condition_depth == 0
}

fn start_segment(
    segments: &mut Vec<Segment>,
    lookup: &mut [Vec<usize>],
    current: &mut [Option<usize>],
    ordinals: &[usize],
    q: usize,
    op_idx: usize,
) -> usize {
    let idx = segments.len();
    segments.push(Segment {
        q,
        ordinal: ordinals[q],
        start: op_idx,
        end: op_idx,
        assigned: q,
    });
    lookup[q].push(idx);
    current[q] = Some(idx);
    idx
}

fn release_finished(
    active: &mut BinaryHeap<Reverse<(usize, usize)>>,
    free_colors: &mut BinaryHeap<Reverse<usize>>,
    next_start: usize,
) {
    while let Some(&Reverse((end, color))) = active.peek() {
        if end >= next_start {
            break;
        }
        active.pop();
        free_colors.push(Reverse(color));
    }
}

fn next_available_color(next_color: &mut usize, fixed: &[bool]) -> usize {
    while *next_color < fixed.len() && fixed[*next_color] {
        *next_color += 1;
    }
    let color = *next_color;
    *next_color += 1;
    color
}

fn assign_colors(segments: &mut [Segment], fixed: &[bool]) {
    let mut order: Vec<usize> = (0..segments.len()).collect();
    order.sort_by_key(|&idx| (segments[idx].start, segments[idx].end, segments[idx].q));

    let mut active: BinaryHeap<Reverse<(usize, usize)>> = BinaryHeap::new();
    let mut free_colors: BinaryHeap<Reverse<usize>> = BinaryHeap::new();
    let mut next_color = 0usize;

    for idx in order {
        release_finished(&mut active, &mut free_colors, segments[idx].start);
        let color = if let Some(Reverse(color)) = free_colors.pop() {
            color
        } else {
            next_available_color(&mut next_color, fixed)
        };
        segments[idx].assigned = color;
        active.push(Reverse((segments[idx].end, color)));
    }
}

fn remap_one(
    q: QubitId,
    fixed: &[bool],
    ordinals: &[usize],
    lookup: &[Vec<usize>],
    segments: &[Segment],
) -> QubitId {
    if q == NO_QUBIT {
        return q;
    }
    let qi = q.0 as usize;
    if fixed.get(qi).copied().unwrap_or(false) {
        return q;
    }
    let seg_idx = lookup
        .get(qi)
        .and_then(|rows| rows.get(ordinals[qi]))
        .copied()
        .expect("temp qubit operand without a compact segment");
    QubitId(segments[seg_idx].assigned as u64)
}

pub fn run(mut ops: Vec<Op>, fixed_qubits: &[QubitId]) -> Vec<Op> {
    let max_q = max_qubit_id(&ops, fixed_qubits);
    let mut fixed = vec![false; max_q];
    for &q in fixed_qubits {
        if q != NO_QUBIT {
            fixed[q.0 as usize] = true;
        }
    }

    let mut segments: Vec<Segment> = Vec::new();
    let mut lookup: Vec<Vec<usize>> = vec![Vec::new(); max_q];
    let mut current: Vec<Option<usize>> = vec![None; max_q];
    let mut ordinals: Vec<usize> = vec![0; max_q];
    let mut condition_depth = 0usize;

    for (op_idx, op) in ops.iter().enumerate() {
        if op.kind == OperationType::PushCondition {
            condition_depth += 1;
            continue;
        }
        if op.kind == OperationType::PopCondition {
            condition_depth = condition_depth.saturating_sub(1);
            continue;
        }

        for q in q_operands(op) {
            if q == NO_QUBIT {
                continue;
            }
            let qi = q.0 as usize;
            if fixed.get(qi).copied().unwrap_or(false) {
                continue;
            }
            let seg_idx = current[qi].unwrap_or_else(|| {
                start_segment(
                    &mut segments,
                    &mut lookup,
                    &mut current,
                    &ordinals,
                    qi,
                    op_idx,
                )
            });
            segments[seg_idx].end = op_idx;
        }

        if is_unconditional_reset(op, condition_depth) {
            let qi = op.q_target.0 as usize;
            if !fixed.get(qi).copied().unwrap_or(false) {
                current[qi] = None;
                ordinals[qi] += 1;
            }
        }
    }

    assign_colors(&mut segments, &fixed);

    let mut rewrite_ordinals: Vec<usize> = vec![0; max_q];
    condition_depth = 0;
    for op in ops.iter_mut() {
        if op.kind == OperationType::PushCondition {
            condition_depth += 1;
            continue;
        }
        if op.kind == OperationType::PopCondition {
            condition_depth = condition_depth.saturating_sub(1);
            continue;
        }

        op.q_control2 = remap_one(op.q_control2, &fixed, &rewrite_ordinals, &lookup, &segments);
        op.q_control1 = remap_one(op.q_control1, &fixed, &rewrite_ordinals, &lookup, &segments);
        let original_target = op.q_target;
        op.q_target = remap_one(op.q_target, &fixed, &rewrite_ordinals, &lookup, &segments);

        if is_unconditional_reset(op, condition_depth) && original_target != NO_QUBIT {
            let qi = original_target.0 as usize;
            if !fixed.get(qi).copied().unwrap_or(false) {
                rewrite_ordinals[qi] += 1;
            }
        }
    }

    if std::env::var("TRACE_RESET_BOUNDED_COMPACT").is_ok() {
        let before = max_q;
        let (after, _, _, _) = analyze_ops(ops.iter());
        eprintln!(
            "RESET_BOUNDED_COMPACT before_q={} after_q={} temp_segments={}",
            before,
            after,
            segments.len()
        );
    }

    ops
}
