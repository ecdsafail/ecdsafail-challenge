//! Dialog-GCD compressed-sidecar path: the round763 block compressor, the
//! runway / composite scratch layout helpers, and the
//! `emit_dialog_gcd_compressed_sidecar_*` block-lifecycle emitters
//! (tobitvector / apply / ipmul / quotient). An alternate, lower-peak encoding
//! of the GCD transcript log; shares the raw-path config levers and comparators
//! from the parent `dialog` module.
use super::*;

pub(crate) fn round763_dedup_enabled() -> bool {
    // EXACT rewrite: the pair ccx(1,3->4) ... ccx(1,3->4) bracketing cx(1->0)
    // cancels (nothing between them touches 1/3/4), so it reduces to bare cx(1->0).
    // 2 CCX -> 0 per direction x ~1064 sites. Default OFF (op-stream reseed).
    std::env::var("DIALOG_GCD_ROUND763_DEDUP").ok().as_deref() == Some("1")
}

pub(crate) fn round763_compress_lever_enabled() -> bool {
    // Reachable-support rewrite of the round763 6->5 sidecar packer. Each raw
    // slot is (b0, b0_and_b1), with b0_and_b1 = b0 & (v<u), so state (0,1) is
    // unreachable on the verifier support. On that support, three CCX collapse
    // to CX and the compressor drops from 9 CCX to 4 CCX per direction.
    std::env::var("DIALOG_GCD_ROUND763_COMPRESS_LEVER")
        .ok()
        .as_deref()
        == Some("1")
}

pub(crate) fn emit_dialog_gcd_round763_compressor(b: &mut B, block: &[QubitId]) {
    assert_eq!(block.len(), 6);
    if round763_compress_lever_enabled() {
        b.cx(block[5], block[3]);
        b.ccx(block[3], block[4], block[5]);
        b.cx(block[1], block[4]);
        b.cx(block[1], block[0]);
        b.ccx(block[4], block[5], block[1]);
        b.cx(block[0], block[2]);
        b.ccx(block[2], block[5], block[0]);
        b.ccx(block[0], block[1], block[5]);
        return;
    }
    b.ccx(block[4], block[5], block[3]);
    b.ccx(block[3], block[4], block[5]);
    b.ccx(block[1], block[2], block[4]);
    if round763_dedup_enabled() {
        b.cx(block[1], block[0]);
    } else {
        b.ccx(block[1], block[3], block[4]);
        b.cx(block[1], block[0]);
        b.ccx(block[1], block[3], block[4]);
    }
    b.ccx(block[4], block[5], block[1]);
    b.ccx(block[0], block[5], block[2]);
    b.ccx(block[2], block[5], block[0]);
    b.ccx(block[0], block[1], block[5]);
}

pub(crate) fn emit_dialog_gcd_round763_compressor_inverse(b: &mut B, block: &[QubitId]) {
    assert_eq!(block.len(), 6);
    if round763_compress_lever_enabled() {
        b.ccx(block[0], block[1], block[5]);
        b.ccx(block[2], block[5], block[0]);
        b.cx(block[0], block[2]);
        b.ccx(block[4], block[5], block[1]);
        b.cx(block[1], block[0]);
        b.cx(block[1], block[4]);
        b.ccx(block[3], block[4], block[5]);
        b.cx(block[5], block[3]);
        return;
    }
    b.ccx(block[0], block[1], block[5]);
    b.ccx(block[2], block[5], block[0]);
    b.ccx(block[0], block[5], block[2]);
    b.ccx(block[4], block[5], block[1]);
    if round763_dedup_enabled() {
        b.cx(block[1], block[0]);
    } else {
        b.ccx(block[1], block[3], block[4]);
        b.cx(block[1], block[0]);
        b.ccx(block[1], block[3], block[4]);
    }
    b.ccx(block[1], block[2], block[4]);
    b.ccx(block[3], block[4], block[5]);
    b.ccx(block[4], block[5], block[3]);
}

pub(crate) fn emit_dialog_gcd_round763_compressed_block_swapper(
    b: &mut B,
    pair: &[QubitId],
    compressed_block: &[QubitId],
    scratch: QubitId,
    slot: usize,
) {
    assert_eq!(pair.len(), 2);
    assert_eq!(compressed_block.len(), 5);
    assert!(slot < 3);
    let mut block = compressed_block.to_vec();
    block.push(scratch);
    emit_dialog_gcd_round763_compressor_inverse(b, &block);
    b.swap(pair[0], block[2 * slot]);
    b.swap(pair[1], block[2 * slot + 1]);
    emit_dialog_gcd_round763_compressor(b, &block);
}

pub(crate) fn dialog_gcd_compressed_sidecar_blocks() -> usize {
    let group_size = dialog_gcd_sidecar_group_size();
    (dialog_gcd_active_iterations() + group_size - 1) / group_size
}

pub(crate) fn dialog_gcd_compressed_sidecar_bits() -> usize {
    dialog_gcd_compressed_sidecar_blocks() * dialog_gcd_block_bits()
}

pub(crate) fn dialog_gcd_compressed_sidecar_block(compressed_log: &[QubitId], step: usize) -> &[QubitId] {
    let block = step / dialog_gcd_sidecar_group_size();
    let bb = dialog_gcd_block_bits();
    let start = block * bb;
    &compressed_log[start..start + bb]
}

pub(crate) fn dialog_gcd_compressed_log_u_high_runway_enabled() -> bool {
    // Prototype, deliberately NOT enabled by configure_ecdsafail_submission_route.
    //
    // The wrapper used to allocate all of u and the complete compressed
    // transcript at once.  Instead, a late transcript suffix can use high u
    // lanes: those cells are not touched until forward replay has shrunk u below
    // their hosts, stay live across terminal-reuse apply, and are consumed by
    // reverse replay before u grows back into them.
    //
    // This is an experimental support-envelope optimization: it relies on the
    // same terminal convergence and width envelope as terminal reuse and
    // variable-width tobitvector.  Default OFF keeps the accepted route
    // byte-identical.
    // K=2: runway layout is now block_bits()-aware (8-bit stride), so it is safe
    // to host the wider K2 transcript blocks on u-high — this is the peak lever.
    std::env::var("DIALOG_GCD_COMPRESSED_LOG_U_HIGH_RUNWAY")
        .ok()
        .as_deref()
        == Some("1")
}

pub(crate) fn dialog_gcd_compressed_log_u_high_runway_blocks() -> usize {
    // Optional tuning cap for the prototype.  The uncapped layout parks the
    // longest suffix; lowering the cap is useful when balancing wrapper savings
    // against reverse-replay scratch pressure.  On the accepted a8d8d5a route,
    // 16 whole blocks is the largest prefix-independent tail runway before the
    // reverse add loses its cheap scratch host.  Keep larger schedules available
    // as an explicit experiment, but default the opt-in prototype to that safe
    // subset.
    std::env::var("DIALOG_GCD_COMPRESSED_LOG_U_HIGH_RUNWAY_BLOCKS")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(16)
}

#[derive(Clone, Debug)]
pub(crate) struct DialogGcdCompressedLogUHighRunway {
    remapped_log: Vec<QubitId>,
    parked_u_indices: Vec<usize>,
}

pub(crate) fn dialog_gcd_slice_intersects(a: &[QubitId], b: &[QubitId]) -> bool {
    a.iter().any(|q| b.contains(q))
}

pub(crate) fn dialog_gcd_runway_layout() -> Vec<(usize, usize)> {
    // Leave the top six u lanes unparked.  The accepted a8d8d5a route hosts a
    // raw 3-step block there whenever the tail is wide enough; reserving those
    // lanes keeps that scratch host disjoint from parked transcript cells.
    let raw_block_bits = 2 * DIALOG_GCD_HIGH_TAIL_ALIAS_GROUP_SIZE;
    let Some(highest_host) = N.checked_sub(raw_block_bits + 1) else {
        return Vec::new();
    };
    let blocks = dialog_gcd_compressed_sidecar_blocks();

    // Find the longest whole-block suffix that fits.  Blocks are assigned in
    // forward order to descending u positions: the earliest parked block gets
    // the highest hosts because it is replayed last and therefore needs the
    // widest inactive-u threshold.
    let first_allowed = blocks.saturating_sub(dialog_gcd_compressed_log_u_high_runway_blocks());
    for first_block in first_allowed..blocks {
        let mut next_host = highest_host;
        let bb = dialog_gcd_block_bits();
        let mut layout = Vec::with_capacity((blocks - first_block) * bb);
        let mut fits = true;
        for block in first_block..blocks {
            let (start, end) = dialog_gcd_compressed_sidecar_block_step_range(block);
            let active_threshold = (start..end)
                .map(dialog_gcd_tobitvector_active_width)
                .max()
                .unwrap_or(1);
            for slot in 0..bb {
                if next_host < active_threshold {
                    fits = false;
                    break;
                }
                layout.push((block * bb + slot, next_host));
                let Some(next) = next_host.checked_sub(1) else {
                    fits = false;
                    break;
                };
                next_host = next;
            }
            if !fits {
                break;
            }
        }
        if fits {
            return layout;
        }
    }
    Vec::new()
}

pub(crate) fn dialog_gcd_allocated_compressed_sidecar_bits() -> usize {
    if dialog_gcd_compressed_log_u_high_runway_enabled() {
        dialog_gcd_compressed_sidecar_bits() - dialog_gcd_runway_layout().len()
    } else {
        dialog_gcd_compressed_sidecar_bits()
    }
}

pub(crate) fn dialog_gcd_build_compressed_log_u_high_runway(
    u: &[QubitId],
    allocated_log: &[QubitId],
) -> Option<DialogGcdCompressedLogUHighRunway> {
    if !dialog_gcd_compressed_log_u_high_runway_enabled() {
        return None;
    }
    assert_eq!(u.len(), N);
    let layout = dialog_gcd_runway_layout();
    if layout.is_empty() {
        return None;
    }

    let expected_allocated = dialog_gcd_compressed_sidecar_bits() - layout.len();
    assert_eq!(allocated_log.len(), expected_allocated);
    let first_relocated = layout[0].0;
    assert_eq!(first_relocated, allocated_log.len());
    let mut remapped_log = allocated_log.to_vec();
    let mut parked_u_indices = Vec::with_capacity(layout.len());
    for (log_index, u_index) in layout {
        // These logical transcript cells are not needed until their late
        // forward blocks, when the width envelope guarantees that u[u_index] is
        // inactive and |0>.  Reverse consumes them before u grows back into the
        // same hosts.
        assert_eq!(log_index, remapped_log.len());
        remapped_log.push(u[u_index]);
        parked_u_indices.push(u_index);
    }
    assert_eq!(remapped_log.len(), dialog_gcd_compressed_sidecar_bits());
    Some(DialogGcdCompressedLogUHighRunway {
        remapped_log,
        parked_u_indices,
    })
}

pub(crate) fn dialog_gcd_release_terminal_u(
    b: &mut B,
    u: &[QubitId],
    runway: Option<&DialogGcdCompressedLogUHighRunway>,
) {
    for (index, &q) in u.iter().enumerate() {
        if runway.is_none_or(|r| !r.parked_u_indices.contains(&index)) {
            b.free(q);
        }
    }
}

pub(crate) fn dialog_gcd_reacquire_terminal_u(
    b: &mut B,
    u: &[QubitId],
    runway: Option<&DialogGcdCompressedLogUHighRunway>,
) {
    for (index, &q) in u.iter().enumerate() {
        if runway.is_none_or(|r| !r.parked_u_indices.contains(&index)) {
            b.reacquire(q);
        }
    }
}

pub(crate) fn dialog_gcd_runway_safe_future_prefix<'a>(
    future: Option<&'a [QubitId]>,
    u: &[QubitId],
    active_width: usize,
) -> Option<&'a [QubitId]> {
    let active_u = &u[..active_width];
    future
        .map(|slice| {
            let safe = slice
                .iter()
                .position(|q| active_u.contains(q))
                .unwrap_or(slice.len());
            &slice[..safe]
        })
        .filter(|slice| !slice.is_empty())
}

pub(crate) fn dialog_gcd_composite_scratch_enabled() -> bool {
    std::env::var("DIALOG_GCD_COMPOSITE_SCRATCH")
        .ok()
        .as_deref()
        == Some("1")
}

pub(crate) fn dialog_gcd_borrow_current_block_enabled() -> bool {
    // The GCD-walk peak (compress_block / shift / reverse_add, all at the same
    // height) is pinned by the composite body-scratch DEFICIT: at the widest
    // (early) steps the materialized sub/add wants ~2*active_width-1 clean lanes
    // for gated+carries, but the only |0> borrow there is the unwritten
    // future-log (block k+1..), leaving a fresh-allocated deficit on top of the
    // resident tx+ty+u+log.
    //
    // Novel observation: the CURRENT block's own compressed cells are also |0>
    // for the entire duration of that block's steps -- forward they are written
    // only by compress_block AFTER every step, reverse they are decompressed
    // into raw_block BEFORE every step -- yet the future-carry slice deliberately
    // starts at block k+1 and never offers them. Folding block k's own cells into
    // the body-scratch borrow shrinks the deficit (a pure qubit relabel, 0 added
    // Toffoli) and is value-exact: the body's measured uncompute restores them to
    // |0> before compress_block/decompress consumes them.
    std::env::var("DIALOG_GCD_BORROW_CURRENT_BLOCK")
        .ok()
        .as_deref()
        == Some("1")
}

pub(crate) fn dialog_gcd_borrow_current_s2_enabled() -> bool {
    // Successor lever to BORROW_CURRENT_BLOCK for the K2 path. The current step's
    // own shift2 (`s2`) cell is provably |0> across its sub/add body window
    // (forward: written only by the later shift phase; reverse: already
    // uncomputed by reverse_unshift) and is restored to |0> by the body's
    // measured uncompute before the shift/unshift consumer. Folding it into the
    // composite-scratch borrow removes one fresh-allocated deficit lane at the
    // width-clamped GCD-walk binder steps (where active_width is pinned at N and
    // the future-log borrow has already shrunk a block), dropping the three
    // compressed-block tobitvector near-binders one qubit. Pure relabel, 0 added
    // Toffoli, value-exact on the reachable GCD support. Default off keeps the
    // accepted op stream byte-identical.
    std::env::var("DIALOG_GCD_BORROW_CURRENT_S2")
        .ok()
        .as_deref()
        == Some("1")
}

pub(crate) struct DialogGcdCompositeScratch {
    lanes: Vec<QubitId>,
    owned: Vec<QubitId>,
}

pub(crate) fn dialog_gcd_build_composite_scratch(
    b: &mut B,
    future: Option<&[QubitId]>,
    u: &[QubitId],
    v: &[QubitId],
    compressed_log: &[QubitId],
    raw_block: &[QubitId],
    active_width: usize,
    step: usize,
) -> DialogGcdCompositeScratch {
    // The selected add/sub body is the dominant consumer of this composite
    // scratch (gated host + borrowed carries). Under the no-physical-c_in body
    // it needs only 2*body_len-1 == 2*body_w-3 lanes (vs 2*active_width-1), and
    // for the untrimmed fastpath body_w == active_width, so the demand drops by
    // exactly 2 lanes — the -1 peak qubit after the gap lane is also reclaimed.
    let body_start = if dialog_gcd_odd_u_lowbit_fastpath_enabled() {
        1
    } else {
        0
    };
    let body_w = dialog_gcd_body_carry_trunc_width(active_width, step);
    let body_len = body_w.saturating_sub(body_start);
    let nocin = dialog_gcd_selected_body_nocin_enabled()
        && !dialog_gcd_selected_body_nocin_keep_pool()
        && body_start >= 1
        && body_len >= 1;
    let want = if nocin {
        // Match the body's exact host demand; never exceed the legacy ask.
        (2 * body_len - 1).min(2 * active_width - 1)
    } else {
        2 * active_width - 1
    };
    let mut lanes = Vec::with_capacity(want);
    let mut push = |q: QubitId| {
        if lanes.len() < want
            && !lanes.contains(&q)
            && !raw_block.contains(&q)
            && !u[..active_width].contains(&q)
            && !v[..active_width].contains(&q)
        {
            lanes.push(q);
        }
    };
    if let Some(future) = dialog_gcd_runway_safe_future_prefix(future, u, active_width) {
        for &q in future {
            push(q);
        }
    }
    if dialog_gcd_borrow_current_block_enabled() {
        // Current block's own compressed cells: |0> across this block's steps
        // (forward written only at compress_block, reverse decompressed before
        // steps). They sit just BELOW the future-carry slice's start (k+1) and
        // are otherwise idle scratch. Restored to |0> by the body's measured
        // uncompute. Skip any that the runway parked onto active u (excluded by
        // push's active-u guard anyway, but kept explicit for clarity).
        let block_cells = dialog_gcd_compressed_sidecar_block(compressed_log, step);
        for &q in block_cells {
            push(q);
        }
    }
    for &q in &v[active_width..] {
        push(q);
    }
    for &q in &u[active_width..] {
        if !compressed_log.contains(&q) {
            push(q);
        }
    }
    if dialog_gcd_borrow_current_s2_enabled() && !raw_block.is_empty() {
        // The CURRENT step's own K2 shift2 (`s2`) cell is |0> across this step's
        // body window: forward it is written only by the later SHIFT phase
        // (after the sub body), reverse it has just been uncomputed by
        // reverse_unshift (before the add body). It is restored to |0> by the
        // body's measured uncompute before either consumer runs. Folding it into
        // the body-scratch borrow shrinks the fresh deficit by one lane at the
        // width-clamped binder steps (the same retiming trick as the current-block
        // compressed cells; pure relabel, 0 added Toffoli). The `push` closure
        // excludes all raw_block cells, so add it explicitly with the same
        // operand/duplicate guards. Disjoint from b0/b0_and_b1 (different slot
        // offset) and from u/v (raw_block is its own register).
        let group_size = dialog_gcd_sidecar_group_size();
        let slot = step % group_size;
        let s2 = raw_block[2 * group_size + slot];
        if lanes.len() < want
            && !lanes.contains(&s2)
            && !u[..active_width].contains(&s2)
            && !v[..active_width].contains(&s2)
        {
            lanes.push(s2);
        }
    }
    let owned = b.alloc_qubits(want - lanes.len());
    lanes.extend_from_slice(&owned);
    DialogGcdCompositeScratch { lanes, owned }
}

pub(crate) fn dialog_gcd_pick_runway_safe_borrow_slice<'a>(
    future: Option<&'a [QubitId]>,
    u: &'a [QubitId],
    compressed_log: &[QubitId],
    active_width: usize,
) -> Option<&'a [QubitId]> {
    if !dialog_gcd_compressed_log_u_high_runway_enabled() {
        return dialog_gcd_pick_borrow_slice(future, u, active_width);
    }

    let safe_future = dialog_gcd_runway_safe_future_prefix(future, u, active_width);
    if dialog_gcd_late_borrow_uv_high_enabled() && active_width >= 1 {
        let want = 2 * active_width - 1;
        let short = safe_future.map_or(true, |slice| slice.len() < want);
        if short && u.len() >= active_width + want {
            let candidate = &u[active_width..active_width + want];
            // Parked cells can still carry unread transcript data.  Be
            // conservative: only use an in-place high-u fallback when it is
            // disjoint from every logical transcript cell, including clean
            // parked cells already consumed by reverse replay.
            if !dialog_gcd_slice_intersects(candidate, compressed_log) {
                return Some(candidate);
            }
        }
    }
    safe_future
}

pub(crate) fn dialog_gcd_host_reverse_raw_block_enabled() -> bool {
    // K=2 packer bring-up: disable raw-block hosting (the raw_block is wider, 9 vs
    // 6, and the hosts assume 6). Allocate raw_block fresh; re-enable K2-aware later.
    if dialog_gcd_k2_enabled() {
        return false;
    }
    std::env::var("DIALOG_GCD_HOST_REVERSE_RAW_BLOCK")
        .ok()
        .as_deref()
        == Some("1")
}

pub(crate) fn dialog_gcd_reverse_raw_block_host<'a>(
    u: &'a [QubitId],
    compressed_log: &'a [QubitId],
    block: usize,
) -> Option<&'a [QubitId]> {
    if !dialog_gcd_host_reverse_raw_block_enabled() {
        return None;
    }
    let (start, _) = dialog_gcd_compressed_sidecar_block_step_range(block);
    let active_width = dialog_gcd_tobitvector_active_width(start);
    let want = 2 * active_width - 1;
    if u.len().saturating_sub(active_width) >= want + 2 * dialog_gcd_sidecar_group_size() {
        let candidate = &u[u.len() - 2 * dialog_gcd_sidecar_group_size()..];
        if !dialog_gcd_compressed_log_u_high_runway_enabled()
            || !dialog_gcd_slice_intersects(candidate, compressed_log)
        {
            return Some(candidate);
        }
    }
    let future_start = (block + 1) * dialog_gcd_block_bits();
    let future = compressed_log.get(future_start..)?;
    let raw_bits = 2 * dialog_gcd_sidecar_group_size();
    if future.len() < want + raw_bits {
        return None;
    }
    if !dialog_gcd_compressed_log_u_high_runway_enabled() {
        return Some(&future[future.len() - raw_bits..]);
    }
    // Keep the raw host after the largest possible carry+gated prefix and away
    // from active u.  With remapped runway cells the old final-six shortcut can
    // alias the growing reverse u prefix.
    future[want..]
        .windows(raw_bits)
        .rev()
        .find(|candidate| !dialog_gcd_slice_intersects(candidate, &u[..active_width]))
}

pub(crate) fn dialog_gcd_forward_raw_block_host<'a>(
    u: &'a [QubitId],
    compressed_log: &'a [QubitId],
    block: usize,
) -> Option<&'a [QubitId]> {
    if !dialog_gcd_host_reverse_raw_block_enabled() {
        return None;
    }
    let (start, _) = dialog_gcd_compressed_sidecar_block_step_range(block);
    let active_width = dialog_gcd_tobitvector_active_width(start);
    let want = 2 * active_width - 1;
    let future_start = (block + 1) * dialog_gcd_block_bits();
    if let Some(future) = compressed_log.get(future_start..) {
        if future.len() >= want + 2 * dialog_gcd_sidecar_group_size() {
            let raw_bits = 2 * dialog_gcd_sidecar_group_size();
            if !dialog_gcd_compressed_log_u_high_runway_enabled() {
                return Some(&future[future.len() - raw_bits..]);
            }
            if let Some(candidate) = future[want..]
                .windows(raw_bits)
                .rev()
                .find(|candidate| !dialog_gcd_slice_intersects(candidate, &u[..active_width]))
            {
                return Some(candidate);
            }
        }
    }
    if u.len().saturating_sub(active_width) >= want + 2 * dialog_gcd_sidecar_group_size() {
        let candidate = &u[u.len() - 2 * dialog_gcd_sidecar_group_size()..];
        if !dialog_gcd_compressed_log_u_high_runway_enabled()
            || !dialog_gcd_slice_intersects(candidate, compressed_log)
        {
            Some(candidate)
        } else {
            None
        }
    } else {
        None
    }
}

pub(crate) fn dialog_gcd_compressed_sidecar_future_carry_slice(
    compressed_log: &[QubitId],
    step: usize,
    active_width: usize,
) -> Option<&[QubitId]> {
    if !dialog_gcd_raw_tobitvector_borrow_future_log_carries_enabled() {
        return None;
    }
    let carry_need = active_width.saturating_sub(1);
    // When hosting the gated register too, request up to carry(n-1)+gated(n)=2n-1
    // clean slots; the consumer splits the returned slice. Graceful: never return
    // fewer than carry_need (so carry borrowing is preserved), never more than
    // what the future region holds.
    let want = if dialog_gcd_host_gated_enabled() {
        2 * active_width - 1
    } else {
        carry_need
    };
    let next_block = step / dialog_gcd_sidecar_group_size() + 1;
    let start = next_block * dialog_gcd_block_bits();
    compressed_log
        .get(start..)
        .filter(|future| future.len() >= carry_need)
        .map(|future| &future[..future.len().min(want)])
}

pub(crate) fn dialog_gcd_compressed_sidecar_block_step_range(block: usize) -> (usize, usize) {
    let group_size = dialog_gcd_sidecar_group_size();
    let start = block * group_size;
    let end = (start + group_size).min(dialog_gcd_active_iterations());
    (start, end)
}

pub(crate) fn dialog_gcd_copy_compressed_block_to_raw(
    b: &mut B,
    compressed_block: &[QubitId],
    raw_block: &[QubitId],
    steps: usize,
) {
    if dialog_gcd_k2_pair_compress_enabled() {
        dialog_gcd_k2_pair_copy_compressed_block_to_raw(b, compressed_block, raw_block, steps);
        return;
    }
    let base_bits = DIALOG_GCD_HIGH_TAIL_ALIAS_BLOCK_BITS; // 5
    let raw_base = 2 * dialog_gcd_sidecar_group_size(); // 6
    assert_eq!(compressed_block.len(), dialog_gcd_block_bits());
    assert_eq!(raw_block.len(), dialog_gcd_raw_block_len());
    let swap_host = dialog_gcd_apply_replay_swap_host_enabled();
    for i in 0..base_bits {
        if swap_host {
            b.swap(compressed_block[i], raw_block[i]);
        } else {
            b.cx(compressed_block[i], raw_block[i]);
        }
    }
    emit_dialog_gcd_round763_compressor_inverse(b, &raw_block[0..raw_base]);
    // K=2 shift2 tail: compressed[5..] -> raw[6..] (raw, no compression).
    for j in base_bits..dialog_gcd_block_bits() {
        let r = raw_base + (j - base_bits);
        if swap_host {
            b.swap(compressed_block[j], raw_block[r]);
        } else {
            b.cx(compressed_block[j], raw_block[r]);
        }
    }
}

pub(crate) fn dialog_gcd_clear_raw_block_copy(
    b: &mut B,
    compressed_block: &[QubitId],
    raw_block: &[QubitId],
    steps: usize,
) {
    if dialog_gcd_k2_pair_compress_enabled() {
        dialog_gcd_k2_pair_clear_raw_block_copy(b, compressed_block, raw_block, steps);
        return;
    }
    let base_bits = DIALOG_GCD_HIGH_TAIL_ALIAS_BLOCK_BITS;
    let raw_base = 2 * dialog_gcd_sidecar_group_size();
    assert_eq!(compressed_block.len(), dialog_gcd_block_bits());
    assert_eq!(raw_block.len(), dialog_gcd_raw_block_len());
    let swap_host = dialog_gcd_apply_replay_swap_host_enabled();
    // Inverse of copy: clear the shift2 tail first, then recompress the base.
    for j in base_bits..dialog_gcd_block_bits() {
        let r = raw_base + (j - base_bits);
        if swap_host {
            b.swap(compressed_block[j], raw_block[r]);
        } else {
            b.cx(compressed_block[j], raw_block[r]);
        }
    }
    emit_dialog_gcd_round763_compressor(b, &raw_block[0..raw_base]);
    for i in 0..base_bits {
        if swap_host {
            b.swap(compressed_block[i], raw_block[i]);
        } else {
            b.cx(compressed_block[i], raw_block[i]);
        }
    }
}

pub(crate) fn emit_dialog_gcd_compressed_sidecar_tobitvector_steps_block_lifecycle(
    b: &mut B,
    u: &[QubitId],
    v: &[QubitId],
    compressed_log: &[QubitId],
    raw_block: &[QubitId],
) {
    assert_eq!(u.len(), N);
    assert_eq!(v.len(), N);
    assert!(raw_block.is_empty() || raw_block.len() == dialog_gcd_raw_block_len());
    assert!(compressed_log.len() >= dialog_gcd_compressed_sidecar_bits());

    for block in 0..dialog_gcd_compressed_sidecar_blocks() {
        let (start, end) = dialog_gcd_compressed_sidecar_block_step_range(block);
        let hosted_raw_block = dialog_gcd_forward_raw_block_host(u, compressed_log, block);
        let owned_raw_block =
            if dialog_gcd_host_reverse_raw_block_enabled() && hosted_raw_block.is_none() {
                b.alloc_qubits(dialog_gcd_raw_block_len())
            } else {
                Vec::new()
            };
        let raw_block = hosted_raw_block.unwrap_or_else(|| {
            if owned_raw_block.is_empty() {
                raw_block
            } else {
                &owned_raw_block
            }
        });
        for step in start..end {
            let slot = step - start;
            let b0 = raw_block[2 * slot];
            let b0_and_b1 = raw_block[2 * slot + 1];
            let active_width = dialog_gcd_tobitvector_active_width(step);
            let u_active = &u[..active_width];
            let v_active = &v[..active_width];
            let compare_bits = dialog_gcd_compare_bits_for_step(step, active_width);

            let future = dialog_gcd_compressed_sidecar_future_carry_slice(
                compressed_log,
                step,
                active_width,
            );
            let composite_scratch = dialog_gcd_composite_scratch_enabled().then(|| {
                dialog_gcd_build_composite_scratch(
                    b,
                    future,
                    u,
                    v,
                    compressed_log,
                    raw_block,
                    active_width,
                    step,
                )
            });
            let borrowed_carries = composite_scratch.as_ref().map_or_else(
                || {
                    dialog_gcd_pick_runway_safe_borrow_slice(
                        future,
                        u,
                        compressed_log,
                        active_width,
                    )
                },
                |scratch| Some(scratch.lanes.as_slice()),
            );

            b.set_phase("dialog_gcd_compressed_block_tobitvector_branch_bits");
            b.cx(v[0], b0);
            if dialog_gcd_fused_branch_bits_enabled() {
                // Fused path derives b0_and_b1 from the in-flight comparator carry
                // and never materializes a separate `cmp` ancilla. Allocating it
                // here would add a dead live-qubit at the branch_bits peak instant
                // (peak is measured by simultaneously-live count, not qubit-id reuse),
                // so it is allocated only on the non-fused branch below.
                if dialog_gcd_branch_bits_host_comparator_enabled() {
                    // Host the comparator's c_in+carries transient on the idle
                    // future-log slice (the same slice the subtract borrows below;
                    // it is unwritten at the comparator instant) so branch_bits no
                    // longer allocates its own peak qubit. Value-exact; the slice is
                    // returned clean by the measured uncompute sweep.
                    dialog_gcd_ccx_cmp_gt_truncated_into_width_hosted(
                        b,
                        u_active,
                        v_active,
                        b0,
                        b0_and_b1,
                        compare_bits,
                        borrowed_carries,
                    );
                } else {
                    dialog_gcd_ccx_cmp_gt_truncated_into_width(
                        b,
                        u_active,
                        v_active,
                        b0,
                        b0_and_b1,
                        compare_bits,
                    );
                }
            } else {
                let cmp = b.alloc_qubit();
                dialog_gcd_cmp_gt_truncated_into_width(b, u_active, v_active, cmp, compare_bits);
                b.ccx(b0, cmp, b0_and_b1);
                dialog_gcd_cmp_gt_truncated_into_width(b, u_active, v_active, cmp, compare_bits);
                b.free(cmp);
            }

            b.set_phase("dialog_gcd_compressed_block_tobitvector_cswap");
            for (i, (&ui, &vi)) in u_active.iter().zip(v_active.iter()).enumerate() {
                if i == 0 && dialog_gcd_odd_u_lowbit_fastpath_enabled() {
                    continue;
                }
                cswap(b, b0_and_b1, ui, vi);
            }

            b.set_phase("dialog_gcd_compressed_block_tobitvector_subtract");
            dialog_gcd_controlled_sub_selected(b, u_active, v_active, b0, borrowed_carries, step);

            b.set_phase("dialog_gcd_compressed_block_tobitvector_shift");
            dialog_gcd_shift_right_assuming_even(b, v_active);
            if dialog_gcd_k2_enabled() {
                // K=2: record shift2 = NOT v_active[0] (v still even after the
                // first shift) into the sidecar, then conditionally shift v_active
                // right once more. Free 1-bit shift is a relabel; this 2nd shift is
                // data-dependent (cswap cascade), ~aw CCX.
                let s2 = raw_block[2 * dialog_gcd_sidecar_group_size() + slot];
                let v0 = v_active[0];
                if std::env::var("DIALOG_GCD_K2_FORCE0").ok().as_deref() != Some("1") {
                    b.cx(v0, s2);
                    b.x(s2);
                }
                for i in 0..v_active.len().saturating_sub(1) {
                    let (lo, hi) = (v_active[i], v_active[i + 1]);
                    cswap(b, s2, lo, hi);
                }
            }
            if let Some(scratch) = composite_scratch {
                b.free_vec(&scratch.owned);
            }
        }

        b.set_phase("dialog_gcd_compressed_block_tobitvector_compress_block");
        let base_bits = DIALOG_GCD_HIGH_TAIL_ALIAS_BLOCK_BITS; // 5
        let compressed_block = dialog_gcd_compressed_sidecar_block(compressed_log, start);
        if dialog_gcd_compressed_log_u_high_runway_enabled() {
            // A parked forward block is first written only after its high-u
            // hosts have left the active prefix.
            assert!(
                !dialog_gcd_slice_intersects(
                    compressed_block,
                    &u[..dialog_gcd_tobitvector_active_width(start)]
                ),
                "compressed-log runway overlaps active forward u prefix at block {block}"
            );
        }
        if dialog_gcd_k2_pair_compress_enabled() {
            dialog_gcd_k2_pair_clear_raw_block_copy(b, compressed_block, raw_block, end - start);
        } else {
            let raw_base = 2 * dialog_gcd_sidecar_group_size(); // 6
            emit_dialog_gcd_round763_compressor(b, &raw_block[0..raw_base]);
            for i in 0..base_bits {
                b.swap(raw_block[i], compressed_block[i]);
            }
            // K=2: stash the shift2 bits raw[raw_base..] into compressed_block[5..].
            for j in base_bits..dialog_gcd_block_bits() {
                b.swap(raw_block[raw_base + (j - base_bits)], compressed_block[j]);
            }
        }
        if !owned_raw_block.is_empty() {
            b.free_vec(&owned_raw_block);
        }
    }
}

pub(crate) fn emit_dialog_gcd_compressed_sidecar_tobitvector_steps_reverse_block_lifecycle(
    b: &mut B,
    u: &[QubitId],
    v: &[QubitId],
    compressed_log: &[QubitId],
    raw_block: &[QubitId],
) {
    assert_eq!(u.len(), N);
    assert_eq!(v.len(), N);
    assert!(raw_block.is_empty() || raw_block.len() == dialog_gcd_raw_block_len());
    assert!(compressed_log.len() >= dialog_gcd_compressed_sidecar_bits());

    for block in (0..dialog_gcd_compressed_sidecar_blocks()).rev() {
        let (start, end) = dialog_gcd_compressed_sidecar_block_step_range(block);
        let compressed_block = dialog_gcd_compressed_sidecar_block(compressed_log, start);
        let hosted_raw_block = dialog_gcd_reverse_raw_block_host(u, compressed_log, block);
        let owned_raw_block =
            if dialog_gcd_host_reverse_raw_block_enabled() && hosted_raw_block.is_none() {
                b.alloc_qubits(dialog_gcd_raw_block_len())
            } else {
                Vec::new()
            };
        let raw_block = hosted_raw_block.unwrap_or_else(|| {
            if owned_raw_block.is_empty() {
                raw_block
            } else {
                &owned_raw_block
            }
        });

        b.set_phase("dialog_gcd_compressed_block_tobitvector_reverse_decompress_block");
        if dialog_gcd_compressed_log_u_high_runway_enabled() {
            // A parked block must be consumed while all of its high-u hosts are
            // outside this block's active prefix.
            assert!(
                !dialog_gcd_slice_intersects(
                    compressed_block,
                    &u[..dialog_gcd_tobitvector_active_width(start)]
                ),
                "compressed-log runway overlaps active reverse u prefix at block {block}"
            );
        }
        {
            let base_bits = DIALOG_GCD_HIGH_TAIL_ALIAS_BLOCK_BITS; // 5
            if dialog_gcd_k2_pair_compress_enabled() {
                dialog_gcd_k2_pair_copy_compressed_block_to_raw(
                    b,
                    compressed_block,
                    raw_block,
                    end - start,
                );
            } else {
                let raw_base = 2 * dialog_gcd_sidecar_group_size(); // 6
                for i in 0..base_bits {
                    b.swap(compressed_block[i], raw_block[i]);
                }
                emit_dialog_gcd_round763_compressor_inverse(b, &raw_block[0..raw_base]);
                // K=2: bring the shift2 bits compressed[5..] -> raw[raw_base..].
                for j in base_bits..dialog_gcd_block_bits() {
                    b.swap(compressed_block[j], raw_block[raw_base + (j - base_bits)]);
                }
            }
        }

        for step in (start..end).rev() {
            let slot = step - start;
            let b0 = raw_block[2 * slot];
            let b0_and_b1 = raw_block[2 * slot + 1];
            let active_width = dialog_gcd_tobitvector_active_width(step);
            let u_active = &u[..active_width];
            let v_active = &v[..active_width];
            let compare_bits = dialog_gcd_compare_bits_for_step(step, active_width);

            b.set_phase("dialog_gcd_compressed_block_tobitvector_reverse_unshift");
            if dialog_gcd_k2_enabled() {
                // mirror of forward K=2: conditional un-shift (reverse cswap order),
                // then uncompute s2 back to |0> (v_active[0] is restored after the
                // un-shift to the value s2 was derived from).
                let s2 = raw_block[2 * dialog_gcd_sidecar_group_size() + slot];
                for i in (0..v_active.len().saturating_sub(1)).rev() {
                    let (lo, hi) = (v_active[i], v_active[i + 1]);
                    cswap(b, s2, lo, hi);
                }
                let v0 = v_active[0];
                if std::env::var("DIALOG_GCD_K2_FORCE0").ok().as_deref() != Some("1") {
                    b.x(s2);
                    b.cx(v0, s2);
                }
            }
            dialog_gcd_unshift_right_assuming_even(b, v_active);

            b.set_phase("dialog_gcd_compressed_block_tobitvector_reverse_add");
            let future = dialog_gcd_compressed_sidecar_future_carry_slice(
                compressed_log,
                step,
                active_width,
            );
            let composite_scratch = dialog_gcd_composite_scratch_enabled().then(|| {
                dialog_gcd_build_composite_scratch(
                    b,
                    future,
                    u,
                    v,
                    compressed_log,
                    raw_block,
                    active_width,
                    step,
                )
            });
            let borrowed_carries = composite_scratch.as_ref().map_or_else(
                || {
                    dialog_gcd_pick_runway_safe_borrow_slice(
                        future,
                        u,
                        compressed_log,
                        active_width,
                    )
                },
                |scratch| Some(scratch.lanes.as_slice()),
            );
            dialog_gcd_controlled_add_selected(b, u_active, v_active, b0, borrowed_carries, step);

            b.set_phase("dialog_gcd_compressed_block_tobitvector_reverse_cswap");
            for (i, (&ui, &vi)) in u_active.iter().zip(v_active.iter()).enumerate() {
                if i == 0 && dialog_gcd_odd_u_lowbit_fastpath_enabled() {
                    continue;
                }
                cswap(b, b0_and_b1, ui, vi);
            }

            b.set_phase("dialog_gcd_compressed_block_tobitvector_reverse_branch_bits");
            if dialog_gcd_fused_branch_bits_enabled() {
                // Fused path: no separate `cmp` ancilla (derives b0_and_b1 from the
                // comparator carry). Allocating it would add a dead live-qubit at the
                // reverse_branch_bits peak instant, so allocate only on the non-fused
                // branch below. See forward lifecycle for the rationale.
                if dialog_gcd_branch_bits_host_comparator_enabled() {
                    // Mirror of the forward path: host the comparator transient on
                    // the idle future-log slice (same slice the add borrowed above).
                    dialog_gcd_ccx_cmp_gt_truncated_into_width_hosted(
                        b,
                        u_active,
                        v_active,
                        b0,
                        b0_and_b1,
                        compare_bits,
                        borrowed_carries,
                    );
                } else {
                    dialog_gcd_ccx_cmp_gt_truncated_into_width(
                        b,
                        u_active,
                        v_active,
                        b0,
                        b0_and_b1,
                        compare_bits,
                    );
                }
            } else {
                let cmp = b.alloc_qubit();
                dialog_gcd_cmp_gt_truncated_into_width(b, u_active, v_active, cmp, compare_bits);
                b.ccx(b0, cmp, b0_and_b1);
                dialog_gcd_cmp_gt_truncated_into_width(b, u_active, v_active, cmp, compare_bits);
                b.free(cmp);
            }
            b.cx(v[0], b0);
            if let Some(scratch) = composite_scratch {
                b.free_vec(&scratch.owned);
            }
        }
        if !owned_raw_block.is_empty() {
            b.free_vec(&owned_raw_block);
        }
    }
}

pub(crate) fn emit_dialog_gcd_compressed_sidecar_apply_bitvector_block_lifecycle(
    b: &mut B,
    compressed_log: &[QubitId],
    x: &[QubitId],
    y: &[QubitId],
    p: U256,
    raw_block: &[QubitId],
) {
    assert_eq!(x.len(), N);
    assert_eq!(y.len(), N);
    assert_eq!(raw_block.len(), dialog_gcd_raw_block_len());

    for block in (0..dialog_gcd_compressed_sidecar_blocks()).rev() {
        let (start, end) = dialog_gcd_compressed_sidecar_block_step_range(block);
        let compressed_block = dialog_gcd_compressed_sidecar_block(compressed_log, start);

        b.set_phase("dialog_gcd_compressed_block_apply_decompress_block");
        dialog_gcd_copy_compressed_block_to_raw(b, compressed_block, raw_block, end - start);
        let clean_scratch = if dialog_gcd_apply_replay_swap_host_enabled() {
            compressed_block
        } else {
            &[]
        };

        for step in (start..end).rev() {
            let slot = step - start;
            let b0 = raw_block[2 * slot];
            let b0_and_b1 = raw_block[2 * slot + 1];

            b.set_phase("dialog_gcd_compressed_block_apply_double_y");
            let apply_k2 = dialog_gcd_k2_enabled()
                && std::env::var("DIALOG_GCD_K2_NO_APPLY").ok().as_deref() != Some("1");
            if apply_k2 && dialog_gcd_apply_fused_fold_enabled() {
                // Fuse mod_double_inplace_fast + cmod_double_inplace_lazy into a
                // single shared carry chain (value-identical; see fn doc).
                let s2 = raw_block[2 * dialog_gcd_sidecar_group_size() + slot];
                dialog_gcd_fused_double_y(b, y, p, s2);
            } else {
                mod_double_inplace_fast(b, y, p);
                if apply_k2 {
                    // mirror the forward K=2 second shift: conditional 2nd double of y.
                    // MUST use the lazy (Solinas, truncated) controlled double so it
                    // composes with the uncontrolled mod_double_inplace_fast above.
                    let s2 = raw_block[2 * dialog_gcd_sidecar_group_size() + slot];
                    cmod_double_inplace_lazy(b, y, p, s2);
                }
            }

            b.set_phase("dialog_gcd_compressed_block_apply_cadd");
            if dialog_gcd_raw_apply_materialized_special_add_enabled() {
                dialog_gcd_cmod_add_materialized_pseudomersenne_with_clean_scratch(
                    b,
                    y,
                    x,
                    b0,
                    p,
                    clean_scratch,
                );
            } else if dialog_gcd_raw_apply_direct_special_add_enabled() {
                dialog_gcd_cmod_add_pseudomersenne_lowq(b, y, x, b0, p);
            } else {
                cmod_add_qq_lowq(b, y, x, b0, p);
            }

            b.set_phase("dialog_gcd_compressed_block_apply_cswap");
            for (&xi, &yi) in x.iter().zip(y.iter()) {
                cswap(b, b0_and_b1, xi, yi);
            }
        }

        b.set_phase("dialog_gcd_compressed_block_apply_clear_block_copy");
        dialog_gcd_clear_raw_block_copy(b, compressed_block, raw_block, end - start);
    }
}

pub(crate) fn emit_dialog_gcd_compressed_sidecar_apply_bitvector_reverse_exact_block_lifecycle(
    b: &mut B,
    compressed_log: &[QubitId],
    x: &[QubitId],
    y: &[QubitId],
    p: U256,
    raw_block: &[QubitId],
) {
    assert_eq!(x.len(), N);
    assert_eq!(y.len(), N);
    assert_eq!(raw_block.len(), dialog_gcd_raw_block_len());

    for block in 0..dialog_gcd_compressed_sidecar_blocks() {
        let (start, end) = dialog_gcd_compressed_sidecar_block_step_range(block);
        let compressed_block = dialog_gcd_compressed_sidecar_block(compressed_log, start);

        b.set_phase("dialog_gcd_compressed_block_apply_reverse_decompress_block");
        dialog_gcd_copy_compressed_block_to_raw(b, compressed_block, raw_block, end - start);
        let clean_scratch = if dialog_gcd_apply_replay_swap_host_enabled() {
            compressed_block
        } else {
            &[]
        };

        for step in start..end {
            let slot = step - start;
            let b0 = raw_block[2 * slot];
            let b0_and_b1 = raw_block[2 * slot + 1];

            b.set_phase("dialog_gcd_compressed_block_apply_reverse_cswap");
            for (&xi, &yi) in x.iter().zip(y.iter()) {
                cswap(b, b0_and_b1, xi, yi);
            }

            b.set_phase("dialog_gcd_compressed_block_apply_reverse_csub");
            if dialog_gcd_raw_apply_reverse_materialized_special_sub_enabled() {
                dialog_gcd_cmod_sub_materialized_pseudomersenne_with_clean_scratch(
                    b,
                    y,
                    x,
                    b0,
                    p,
                    clean_scratch,
                );
            } else if dialog_gcd_raw_apply_reverse_fast_sub_enabled() {
                cmod_sub_qq(b, y, x, b0, p);
            } else {
                cmod_sub_qq_lowq(b, y, x, b0, p);
            }

            b.set_phase("dialog_gcd_compressed_block_apply_reverse_halve_y");
            let apply_k2 = dialog_gcd_k2_enabled()
                && std::env::var("DIALOG_GCD_K2_NO_APPLY").ok().as_deref() != Some("1");
            if apply_k2
                && dialog_gcd_apply_fused_fold_enabled()
                && std::env::var("DIALOG_GCD_FUSE_HALVE_OFF").ok().as_deref() != Some("1")
            {
                // Fuse mod_halve_inplace_fast + cmod_halve_inplace_lazy into a
                // single shared borrow chain (exact inverse of the fused double;
                // see fn doc on dialog_gcd_fused_halve_y).
                let s2 = raw_block[2 * dialog_gcd_sidecar_group_size() + slot];
                dialog_gcd_fused_halve_y(b, y, p, s2);
            } else {
                mod_halve_inplace_fast(b, y, p);
                if apply_k2 {
                    // mirror the forward K=2 second shift: conditional 2nd halve of y.
                    // MUST use the lazy (Solinas, truncated) controlled halve to match.
                    let s2 = raw_block[2 * dialog_gcd_sidecar_group_size() + slot];
                    cmod_halve_inplace_lazy(b, y, p, s2);
                }
            }
        }

        b.set_phase("dialog_gcd_compressed_block_apply_reverse_clear_block_copy");
        dialog_gcd_clear_raw_block_copy(b, compressed_block, raw_block, end - start);
    }
}

pub(crate) fn emit_dialog_gcd_compressed_sidecar_tobitvector_steps(
    b: &mut B,
    u: &[QubitId],
    v: &[QubitId],
    compressed_log: &[QubitId],
    pair: &[QubitId],
    scratch: QubitId,
) {
    assert_eq!(u.len(), N);
    assert_eq!(v.len(), N);
    assert_eq!(pair.len(), 2);
    assert!(compressed_log.len() >= dialog_gcd_compressed_sidecar_bits());

    for step in 0..dialog_gcd_active_iterations() {
        let b0 = pair[0];
        let b0_and_b1 = pair[1];
        let cmp = b.alloc_qubit();
        let active_width = dialog_gcd_tobitvector_active_width(step);
        let u_active = &u[..active_width];
        let v_active = &v[..active_width];
        let compare_bits = dialog_gcd_compare_bits_for_step(step, active_width);

        b.set_phase("dialog_gcd_compressed_sidecar_tobitvector_branch_bits");
        b.cx(v[0], b0);
        if dialog_gcd_fused_branch_bits_enabled() {
            dialog_gcd_ccx_cmp_gt_truncated_into_width(
                b,
                u_active,
                v_active,
                b0,
                b0_and_b1,
                compare_bits,
            );
        } else {
            dialog_gcd_cmp_gt_truncated_into_width(b, u_active, v_active, cmp, compare_bits);
            b.ccx(b0, cmp, b0_and_b1);
            dialog_gcd_cmp_gt_truncated_into_width(b, u_active, v_active, cmp, compare_bits);
        }
        b.free(cmp);

        b.set_phase("dialog_gcd_compressed_sidecar_tobitvector_cswap");
        for (i, (&ui, &vi)) in u_active.iter().zip(v_active.iter()).enumerate() {
            if i == 0 && dialog_gcd_odd_u_lowbit_fastpath_enabled() {
                continue;
            }
            cswap(b, b0_and_b1, ui, vi);
        }

        b.set_phase("dialog_gcd_compressed_sidecar_tobitvector_subtract");
        let borrowed_carries =
            dialog_gcd_compressed_sidecar_future_carry_slice(compressed_log, step, active_width);
        dialog_gcd_controlled_sub_selected(b, u_active, v_active, b0, borrowed_carries, step);

        b.set_phase("dialog_gcd_compressed_sidecar_tobitvector_shift");
        dialog_gcd_shift_right_assuming_even(b, v_active);

        b.set_phase("dialog_gcd_compressed_sidecar_tobitvector_absorb_pair");
        let block = dialog_gcd_compressed_sidecar_block(compressed_log, step);
        emit_dialog_gcd_round763_compressed_block_swapper(
            b,
            pair,
            block,
            scratch,
            step % DIALOG_GCD_HIGH_TAIL_ALIAS_GROUP_SIZE,
        );
    }
}

pub(crate) fn emit_dialog_gcd_compressed_sidecar_tobitvector_steps_reverse(
    b: &mut B,
    u: &[QubitId],
    v: &[QubitId],
    compressed_log: &[QubitId],
    pair: &[QubitId],
    scratch: QubitId,
) {
    assert_eq!(u.len(), N);
    assert_eq!(v.len(), N);
    assert_eq!(pair.len(), 2);
    assert!(compressed_log.len() >= dialog_gcd_compressed_sidecar_bits());

    for step in (0..dialog_gcd_active_iterations()).rev() {
        b.set_phase("dialog_gcd_compressed_sidecar_tobitvector_reverse_load_pair");
        let block = dialog_gcd_compressed_sidecar_block(compressed_log, step);
        emit_dialog_gcd_round763_compressed_block_swapper(
            b,
            pair,
            block,
            scratch,
            step % DIALOG_GCD_HIGH_TAIL_ALIAS_GROUP_SIZE,
        );

        let b0 = pair[0];
        let b0_and_b1 = pair[1];
        let cmp = b.alloc_qubit();
        let active_width = dialog_gcd_tobitvector_active_width(step);
        let u_active = &u[..active_width];
        let v_active = &v[..active_width];
        let compare_bits = dialog_gcd_compare_bits_for_step(step, active_width);

        b.set_phase("dialog_gcd_compressed_sidecar_tobitvector_reverse_unshift");
        dialog_gcd_unshift_right_assuming_even(b, v_active);

        b.set_phase("dialog_gcd_compressed_sidecar_tobitvector_reverse_add");
        let borrowed_carries =
            dialog_gcd_compressed_sidecar_future_carry_slice(compressed_log, step, active_width);
        dialog_gcd_controlled_add_selected(b, u_active, v_active, b0, borrowed_carries, step);

        b.set_phase("dialog_gcd_compressed_sidecar_tobitvector_reverse_cswap");
        for (i, (&ui, &vi)) in u_active.iter().zip(v_active.iter()).enumerate() {
            if i == 0 && dialog_gcd_odd_u_lowbit_fastpath_enabled() {
                continue;
            }
            cswap(b, b0_and_b1, ui, vi);
        }

        b.set_phase("dialog_gcd_compressed_sidecar_tobitvector_reverse_branch_bits");
        if dialog_gcd_fused_branch_bits_enabled() {
            dialog_gcd_ccx_cmp_gt_truncated_into_width(
                b,
                u_active,
                v_active,
                b0,
                b0_and_b1,
                compare_bits,
            );
        } else {
            dialog_gcd_cmp_gt_truncated_into_width(b, u_active, v_active, cmp, compare_bits);
            b.ccx(b0, cmp, b0_and_b1);
            dialog_gcd_cmp_gt_truncated_into_width(b, u_active, v_active, cmp, compare_bits);
        }
        b.free(cmp);
        b.cx(v[0], b0);
    }
}

pub(crate) fn emit_dialog_gcd_compressed_sidecar_apply_bitvector(
    b: &mut B,
    compressed_log: &[QubitId],
    x: &[QubitId],
    y: &[QubitId],
    p: U256,
    pair: &[QubitId],
    scratch: QubitId,
) {
    assert_eq!(x.len(), N);
    assert_eq!(y.len(), N);
    assert_eq!(pair.len(), 2);

    for step in (0..dialog_gcd_active_iterations()).rev() {
        b.set_phase("dialog_gcd_compressed_sidecar_apply_load_pair");
        let block = dialog_gcd_compressed_sidecar_block(compressed_log, step);
        emit_dialog_gcd_round763_compressed_block_swapper(
            b,
            pair,
            block,
            scratch,
            step % DIALOG_GCD_HIGH_TAIL_ALIAS_GROUP_SIZE,
        );

        let b0 = pair[0];
        let b0_and_b1 = pair[1];

        b.set_phase("dialog_gcd_compressed_sidecar_apply_double_y");
        mod_double_inplace_fast(b, y, p);

        b.set_phase("dialog_gcd_compressed_sidecar_apply_cadd");
        if dialog_gcd_raw_apply_materialized_special_add_enabled() {
            dialog_gcd_cmod_add_materialized_pseudomersenne(b, y, x, b0, p);
        } else if dialog_gcd_raw_apply_direct_special_add_enabled() {
            dialog_gcd_cmod_add_pseudomersenne_lowq(b, y, x, b0, p);
        } else {
            cmod_add_qq_lowq(b, y, x, b0, p);
        }

        b.set_phase("dialog_gcd_compressed_sidecar_apply_cswap");
        for (&xi, &yi) in x.iter().zip(y.iter()) {
            cswap(b, b0_and_b1, xi, yi);
        }

        b.set_phase("dialog_gcd_compressed_sidecar_apply_unload_pair");
        let block = dialog_gcd_compressed_sidecar_block(compressed_log, step);
        emit_dialog_gcd_round763_compressed_block_swapper(
            b,
            pair,
            block,
            scratch,
            step % DIALOG_GCD_HIGH_TAIL_ALIAS_GROUP_SIZE,
        );
    }
}

pub(crate) fn emit_dialog_gcd_compressed_sidecar_apply_bitvector_reverse_exact(
    b: &mut B,
    compressed_log: &[QubitId],
    x: &[QubitId],
    y: &[QubitId],
    p: U256,
    pair: &[QubitId],
    scratch: QubitId,
) {
    assert_eq!(x.len(), N);
    assert_eq!(y.len(), N);
    assert_eq!(pair.len(), 2);

    for step in 0..dialog_gcd_active_iterations() {
        b.set_phase("dialog_gcd_compressed_sidecar_apply_reverse_load_pair");
        let block = dialog_gcd_compressed_sidecar_block(compressed_log, step);
        emit_dialog_gcd_round763_compressed_block_swapper(
            b,
            pair,
            block,
            scratch,
            step % DIALOG_GCD_HIGH_TAIL_ALIAS_GROUP_SIZE,
        );

        let b0 = pair[0];
        let b0_and_b1 = pair[1];

        b.set_phase("dialog_gcd_compressed_sidecar_apply_reverse_cswap");
        for (&xi, &yi) in x.iter().zip(y.iter()) {
            cswap(b, b0_and_b1, xi, yi);
        }

        b.set_phase("dialog_gcd_compressed_sidecar_apply_reverse_csub");
        if dialog_gcd_raw_apply_reverse_materialized_special_sub_enabled() {
            dialog_gcd_cmod_sub_materialized_pseudomersenne(b, y, x, b0, p);
        } else if dialog_gcd_raw_apply_reverse_fast_sub_enabled() {
            cmod_sub_qq(b, y, x, b0, p);
        } else {
            cmod_sub_qq_lowq(b, y, x, b0, p);
        }

        b.set_phase("dialog_gcd_compressed_sidecar_apply_reverse_halve_y");
        mod_halve_inplace_fast(b, y, p);

        b.set_phase("dialog_gcd_compressed_sidecar_apply_reverse_unload_pair");
        let block = dialog_gcd_compressed_sidecar_block(compressed_log, step);
        emit_dialog_gcd_round763_compressed_block_swapper(
            b,
            pair,
            block,
            scratch,
            step % DIALOG_GCD_HIGH_TAIL_ALIAS_GROUP_SIZE,
        );
    }
}

pub(crate) fn emit_dialog_gcd_compressed_sidecar_ipmul_block_lifecycle(
    b: &mut B,
    factor: &[QubitId],
    target: &[QubitId],
    p: U256,
) {
    assert_eq!(factor.len(), N);
    assert_eq!(target.len(), N);

    let compressed_log = b.alloc_qubits(dialog_gcd_allocated_compressed_sidecar_bits());
    let raw_block = if dialog_gcd_host_reverse_raw_block_enabled() {
        Vec::new()
    } else {
        b.alloc_qubits(dialog_gcd_raw_block_len())
    };
    let u = b.alloc_qubits(N);
    let runway = dialog_gcd_build_compressed_log_u_high_runway(&u, &compressed_log);
    let replay_log = runway
        .as_ref()
        .map_or(compressed_log.as_slice(), |r| r.remapped_log.as_slice());
    b.set_phase("dialog_gcd_compressed_block_ipmul_load_p");
    for i in 0..N {
        if bit(p, i) {
            b.x(u[i]);
        }
    }

    b.set_phase("dialog_gcd_compressed_block_ipmul_tobitvector");
    emit_dialog_gcd_compressed_sidecar_tobitvector_steps_block_lifecycle(
        b, &u, factor, replay_log, &raw_block,
    );

    if dialog_gcd_raw_ipmul_terminal_reuse_enabled() {
        b.set_phase("dialog_gcd_compressed_block_ipmul_release_terminal_u");
        b.x(u[0]);
        dialog_gcd_release_terminal_u(b, &u, runway.as_ref());

        b.set_phase("dialog_gcd_compressed_block_ipmul_apply_bitvector_reuse_factor_zero");
        let apply_raw_block = if dialog_gcd_host_reverse_raw_block_enabled() {
            b.alloc_qubits(dialog_gcd_raw_block_len())
        } else {
            Vec::new()
        };
        emit_dialog_gcd_compressed_sidecar_apply_bitvector_block_lifecycle(
            b,
            replay_log,
            target,
            factor,
            p,
            if apply_raw_block.is_empty() {
                &raw_block
            } else {
                &apply_raw_block
            },
        );
        if !apply_raw_block.is_empty() {
            b.free_vec(&apply_raw_block);
        }

        if dialog_gcd_raw_ipmul_clear_p_residual_enabled() {
            b.set_phase("dialog_gcd_compressed_block_ipmul_clear_p_residual_source_lane");
            for i in 0..N {
                if bit(p, i) {
                    b.x(target[i]);
                }
            }
        }

        b.set_phase("dialog_gcd_compressed_block_ipmul_swap_product_into_target");
        for i in 0..N {
            b.swap(target[i], factor[i]);
        }

        b.set_phase("dialog_gcd_compressed_block_ipmul_reacquire_terminal_u");
        dialog_gcd_reacquire_terminal_u(b, &u, runway.as_ref());
        b.set_phase("dialog_gcd_compressed_block_ipmul_seed_terminal_u");
        b.x(u[0]);

        b.set_phase("dialog_gcd_compressed_block_ipmul_uncompute_tobitvector");
        emit_dialog_gcd_compressed_sidecar_tobitvector_steps_reverse_block_lifecycle(
            b, &u, factor, replay_log, &raw_block,
        );

        b.set_phase("dialog_gcd_compressed_block_ipmul_unload_p");
        for i in 0..N {
            if bit(p, i) {
                b.x(u[i]);
            }
        }
        if !b.k2_shift2_log.is_empty() {
            let log = std::mem::take(&mut b.k2_shift2_log);
            b.free_vec(&log);
        }
        b.free_vec(&u);
        if !raw_block.is_empty() {
            b.free_vec(&raw_block);
        }
        b.free_vec(&compressed_log);
        return;
    }

    let tmp = b.alloc_qubits(N);
    b.set_phase("dialog_gcd_compressed_block_ipmul_apply_bitvector");
    emit_dialog_gcd_compressed_sidecar_apply_bitvector_block_lifecycle(
        b, replay_log, target, &tmp, p, &raw_block,
    );

    b.set_phase("dialog_gcd_compressed_block_ipmul_swap_product_into_target");
    for i in 0..N {
        b.swap(target[i], tmp[i]);
    }

    b.set_phase("dialog_gcd_compressed_block_ipmul_free_zero_tmp");
    b.free_vec(&tmp);

    b.set_phase("dialog_gcd_compressed_block_ipmul_uncompute_tobitvector");
    emit_dialog_gcd_compressed_sidecar_tobitvector_steps_reverse_block_lifecycle(
        b, &u, factor, replay_log, &raw_block,
    );

    b.set_phase("dialog_gcd_compressed_block_ipmul_unload_p");
    for i in 0..N {
        if bit(p, i) {
            b.x(u[i]);
        }
    }
    if !b.k2_shift2_log.is_empty() {
        let log = std::mem::take(&mut b.k2_shift2_log);
        b.free_vec(&log);
    }
    b.free_vec(&u);
    b.free_vec(&raw_block);
    b.free_vec(&compressed_log);
}

pub(crate) fn emit_dialog_gcd_compressed_sidecar_ipmul(
    b: &mut B,
    factor: &[QubitId],
    target: &[QubitId],
    p: U256,
) {
    assert_eq!(factor.len(), N);
    assert_eq!(target.len(), N);

    if dialog_gcd_compressed_block_lifecycle_enabled() {
        emit_dialog_gcd_compressed_sidecar_ipmul_block_lifecycle(b, factor, target, p);
        return;
    }

    let compressed_log = b.alloc_qubits(dialog_gcd_compressed_sidecar_bits());
    let pair = b.alloc_qubits(2);
    let compressor_scratch = b.alloc_qubit();
    let u = b.alloc_qubits(N);
    b.set_phase("dialog_gcd_compressed_sidecar_ipmul_load_p");
    for i in 0..N {
        if bit(p, i) {
            b.x(u[i]);
        }
    }

    b.set_phase("dialog_gcd_compressed_sidecar_ipmul_tobitvector");
    emit_dialog_gcd_compressed_sidecar_tobitvector_steps(
        b,
        &u,
        factor,
        &compressed_log,
        &pair,
        compressor_scratch,
    );

    if dialog_gcd_raw_ipmul_terminal_reuse_enabled() {
        b.set_phase("dialog_gcd_compressed_sidecar_ipmul_release_terminal_u");
        b.x(u[0]);
        b.free_vec(&u);

        b.set_phase("dialog_gcd_compressed_sidecar_ipmul_apply_bitvector_reuse_factor_zero");
        emit_dialog_gcd_compressed_sidecar_apply_bitvector(
            b,
            &compressed_log,
            target,
            factor,
            p,
            &pair,
            compressor_scratch,
        );

        if dialog_gcd_raw_ipmul_clear_p_residual_enabled() {
            b.set_phase("dialog_gcd_compressed_sidecar_ipmul_clear_p_residual_source_lane");
            for i in 0..N {
                if bit(p, i) {
                    b.x(target[i]);
                }
            }
        }

        b.set_phase("dialog_gcd_compressed_sidecar_ipmul_swap_product_into_target");
        for i in 0..N {
            b.swap(target[i], factor[i]);
        }

        b.set_phase("dialog_gcd_compressed_sidecar_ipmul_reacquire_terminal_u");
        b.reacquire_vec(&u);
        b.set_phase("dialog_gcd_compressed_sidecar_ipmul_seed_terminal_u");
        b.x(u[0]);

        b.set_phase("dialog_gcd_compressed_sidecar_ipmul_uncompute_tobitvector");
        emit_dialog_gcd_compressed_sidecar_tobitvector_steps_reverse(
            b,
            &u,
            factor,
            &compressed_log,
            &pair,
            compressor_scratch,
        );

        b.set_phase("dialog_gcd_compressed_sidecar_ipmul_unload_p");
        for i in 0..N {
            if bit(p, i) {
                b.x(u[i]);
            }
        }
        b.free_vec(&u);
        b.free(compressor_scratch);
        b.free_vec(&pair);
        b.free_vec(&compressed_log);
        return;
    }

    let tmp = b.alloc_qubits(N);
    b.set_phase("dialog_gcd_compressed_sidecar_ipmul_apply_bitvector");
    emit_dialog_gcd_compressed_sidecar_apply_bitvector(
        b,
        &compressed_log,
        target,
        &tmp,
        p,
        &pair,
        compressor_scratch,
    );

    b.set_phase("dialog_gcd_compressed_sidecar_ipmul_swap_product_into_target");
    for i in 0..N {
        b.swap(target[i], tmp[i]);
    }

    b.set_phase("dialog_gcd_compressed_sidecar_ipmul_free_zero_tmp");
    b.free_vec(&tmp);

    b.set_phase("dialog_gcd_compressed_sidecar_ipmul_uncompute_tobitvector");
    emit_dialog_gcd_compressed_sidecar_tobitvector_steps_reverse(
        b,
        &u,
        factor,
        &compressed_log,
        &pair,
        compressor_scratch,
    );

    b.set_phase("dialog_gcd_compressed_sidecar_ipmul_unload_p");
    for i in 0..N {
        if bit(p, i) {
            b.x(u[i]);
        }
    }
    b.free_vec(&u);
    b.free(compressor_scratch);
    b.free_vec(&pair);
    b.free_vec(&compressed_log);
}

pub(crate) fn emit_dialog_gcd_compressed_sidecar_quotient_block_lifecycle(
    b: &mut B,
    factor: &[QubitId],
    target: &[QubitId],
    p: U256,
) {
    assert_eq!(factor.len(), N);
    assert_eq!(target.len(), N);

    let compressed_log = b.alloc_qubits(dialog_gcd_allocated_compressed_sidecar_bits());
    let raw_block = if dialog_gcd_host_reverse_raw_block_enabled() {
        Vec::new()
    } else {
        b.alloc_qubits(dialog_gcd_raw_block_len())
    };
    let u = b.alloc_qubits(N);
    let runway = dialog_gcd_build_compressed_log_u_high_runway(&u, &compressed_log);
    let replay_log = runway
        .as_ref()
        .map_or(compressed_log.as_slice(), |r| r.remapped_log.as_slice());
    b.set_phase("dialog_gcd_compressed_block_quotient_load_p");
    for i in 0..N {
        if bit(p, i) {
            b.x(u[i]);
        }
    }

    b.set_phase("dialog_gcd_compressed_block_quotient_tobitvector");
    emit_dialog_gcd_compressed_sidecar_tobitvector_steps_block_lifecycle(
        b, &u, factor, replay_log, &raw_block,
    );

    if dialog_gcd_raw_quotient_terminal_reuse_enabled() {
        b.set_phase("dialog_gcd_compressed_block_quotient_release_terminal_u");
        b.x(u[0]);
        dialog_gcd_release_terminal_u(b, &u, runway.as_ref());

        b.set_phase("dialog_gcd_compressed_block_quotient_apply_reverse_reuse_factor_zero");
        let apply_raw_block = if dialog_gcd_host_reverse_raw_block_enabled() {
            b.alloc_qubits(dialog_gcd_raw_block_len())
        } else {
            Vec::new()
        };
        emit_dialog_gcd_compressed_sidecar_apply_bitvector_reverse_exact_block_lifecycle(
            b,
            replay_log,
            factor,
            target,
            p,
            if apply_raw_block.is_empty() {
                &raw_block
            } else {
                &apply_raw_block
            },
        );
        if !apply_raw_block.is_empty() {
            b.free_vec(&apply_raw_block);
        }

        b.set_phase("dialog_gcd_compressed_block_quotient_swap_quotient_into_target");
        for i in 0..N {
            b.swap(target[i], factor[i]);
        }

        b.set_phase("dialog_gcd_compressed_block_quotient_reacquire_terminal_u");
        dialog_gcd_reacquire_terminal_u(b, &u, runway.as_ref());
        b.set_phase("dialog_gcd_compressed_block_quotient_seed_terminal_u");
        b.x(u[0]);

        b.set_phase("dialog_gcd_compressed_block_quotient_uncompute_tobitvector");
        emit_dialog_gcd_compressed_sidecar_tobitvector_steps_reverse_block_lifecycle(
            b, &u, factor, replay_log, &raw_block,
        );

        b.set_phase("dialog_gcd_compressed_block_quotient_unload_p");
        for i in 0..N {
            if bit(p, i) {
                b.x(u[i]);
            }
        }
        if !b.k2_shift2_log.is_empty() {
            let log = std::mem::take(&mut b.k2_shift2_log);
            b.free_vec(&log);
        }
        b.free_vec(&u);
        if !raw_block.is_empty() {
            b.free_vec(&raw_block);
        }
        b.free_vec(&compressed_log);
        return;
    }

    b.set_phase("dialog_gcd_compressed_block_quotient_apply_reverse");
    emit_dialog_gcd_compressed_sidecar_apply_bitvector_reverse_exact_block_lifecycle(
        b, replay_log, factor, target, p, &raw_block,
    );

    b.set_phase("dialog_gcd_compressed_block_quotient_uncompute_tobitvector");
    emit_dialog_gcd_compressed_sidecar_tobitvector_steps_reverse_block_lifecycle(
        b, &u, factor, replay_log, &raw_block,
    );

    b.set_phase("dialog_gcd_compressed_block_quotient_unload_p");
    for i in 0..N {
        if bit(p, i) {
            b.x(u[i]);
        }
    }
    if !b.k2_shift2_log.is_empty() {
        let log = std::mem::take(&mut b.k2_shift2_log);
        b.free_vec(&log);
    }
    b.free_vec(&u);
    b.free_vec(&raw_block);
    b.free_vec(&compressed_log);
}

pub(crate) fn emit_dialog_gcd_compressed_sidecar_quotient(
    b: &mut B,
    factor: &[QubitId],
    target: &[QubitId],
    p: U256,
) {
    assert_eq!(factor.len(), N);
    assert_eq!(target.len(), N);

    if dialog_gcd_compressed_block_lifecycle_enabled() {
        emit_dialog_gcd_compressed_sidecar_quotient_block_lifecycle(b, factor, target, p);
        return;
    }

    let compressed_log = b.alloc_qubits(dialog_gcd_compressed_sidecar_bits());
    let pair = b.alloc_qubits(2);
    let compressor_scratch = b.alloc_qubit();
    let u = b.alloc_qubits(N);
    b.set_phase("dialog_gcd_compressed_sidecar_quotient_load_p");
    for i in 0..N {
        if bit(p, i) {
            b.x(u[i]);
        }
    }

    b.set_phase("dialog_gcd_compressed_sidecar_quotient_tobitvector");
    emit_dialog_gcd_compressed_sidecar_tobitvector_steps(
        b,
        &u,
        factor,
        &compressed_log,
        &pair,
        compressor_scratch,
    );

    if dialog_gcd_raw_quotient_terminal_reuse_enabled() {
        b.set_phase("dialog_gcd_compressed_sidecar_quotient_release_terminal_u");
        b.x(u[0]);
        b.free_vec(&u);

        b.set_phase("dialog_gcd_compressed_sidecar_quotient_apply_reverse_reuse_factor_zero");
        emit_dialog_gcd_compressed_sidecar_apply_bitvector_reverse_exact(
            b,
            &compressed_log,
            factor,
            target,
            p,
            &pair,
            compressor_scratch,
        );

        b.set_phase("dialog_gcd_compressed_sidecar_quotient_swap_quotient_into_target");
        for i in 0..N {
            b.swap(target[i], factor[i]);
        }

        b.set_phase("dialog_gcd_compressed_sidecar_quotient_reacquire_terminal_u");
        b.reacquire_vec(&u);
        b.set_phase("dialog_gcd_compressed_sidecar_quotient_seed_terminal_u");
        b.x(u[0]);

        b.set_phase("dialog_gcd_compressed_sidecar_quotient_uncompute_tobitvector");
        emit_dialog_gcd_compressed_sidecar_tobitvector_steps_reverse(
            b,
            &u,
            factor,
            &compressed_log,
            &pair,
            compressor_scratch,
        );

        b.set_phase("dialog_gcd_compressed_sidecar_quotient_unload_p");
        for i in 0..N {
            if bit(p, i) {
                b.x(u[i]);
            }
        }
        b.free_vec(&u);
        b.free(compressor_scratch);
        b.free_vec(&pair);
        b.free_vec(&compressed_log);
        return;
    }

    b.set_phase("dialog_gcd_compressed_sidecar_quotient_apply_reverse");
    emit_dialog_gcd_compressed_sidecar_apply_bitvector_reverse_exact(
        b,
        &compressed_log,
        factor,
        target,
        p,
        &pair,
        compressor_scratch,
    );

    b.set_phase("dialog_gcd_compressed_sidecar_quotient_uncompute_tobitvector");
    emit_dialog_gcd_compressed_sidecar_tobitvector_steps_reverse(
        b,
        &u,
        factor,
        &compressed_log,
        &pair,
        compressor_scratch,
    );

    b.set_phase("dialog_gcd_compressed_sidecar_quotient_unload_p");
    for i in 0..N {
        if bit(p, i) {
            b.x(u[i]);
        }
    }
    b.free_vec(&u);
    b.free(compressor_scratch);
    b.free_vec(&pair);
    b.free_vec(&compressed_log);
}


pub(crate) fn emit_dialog_gcd_k2_pair_core_encoder(b: &mut B, core: &[QubitId]) {
    assert_eq!(core.len(), 5);
    b.x(core[1]);
    b.cx(core[0], core[3]);
    b.ccx(core[1], core[3], core[0]);
    b.cx(core[0], core[1]);
    b.cx(core[2], core[3]);
    b.ccx(core[0], core[3], core[2]);
    b.ccx(core[1], core[2], core[0]);
    b.x(core[0]);
    b.x(core[3]);
    b.ccx(core[1], core[4], core[0]);
    b.ccx(core[0], core[2], core[4]);
    b.ccx(core[1], core[4], core[0]);
}

pub(crate) fn emit_dialog_gcd_k2_pair_core_encoder_inverse(b: &mut B, core: &[QubitId]) {
    assert_eq!(core.len(), 5);
    b.ccx(core[1], core[4], core[0]);
    b.ccx(core[0], core[2], core[4]);
    b.ccx(core[1], core[4], core[0]);
    b.x(core[3]);
    b.x(core[0]);
    b.ccx(core[1], core[2], core[0]);
    b.ccx(core[0], core[3], core[2]);
    b.cx(core[2], core[3]);
    b.cx(core[0], core[1]);
    b.ccx(core[1], core[3], core[0]);
    b.cx(core[0], core[3]);
    b.x(core[1]);
}

pub(crate) fn dialog_gcd_k2_pair_core(raw_block: &[QubitId]) -> [QubitId; 5] {
    assert_eq!(raw_block.len(), 6);
    [
        raw_block[0], // first step b0
        raw_block[1], // first step b0_and_b1
        raw_block[4], // first step shift2
        raw_block[2], // second step b0
        raw_block[3], // second step b0_and_b1
    ]
}

pub(crate) fn dialog_gcd_k2_pair_copy_compressed_block_to_raw(
    b: &mut B,
    compressed_block: &[QubitId],
    raw_block: &[QubitId],
    steps: usize,
) {
    assert_eq!(compressed_block.len(), DIALOG_GCD_HIGH_TAIL_ALIAS_BLOCK_BITS);
    assert_eq!(raw_block.len(), 6);
    assert!((1..=2).contains(&steps));
    let swap_host = dialog_gcd_apply_replay_swap_host_enabled();
    if steps == 1 {
        let raw_encoded = [raw_block[0], raw_block[1], raw_block[4]];
        for (&c, &r) in compressed_block.iter().take(3).zip(raw_encoded.iter()) {
            if swap_host {
                b.swap(c, r);
            } else {
                b.cx(c, r);
            }
        }
        return;
    }
    let raw_encoded = [raw_block[1], raw_block[4], raw_block[2], raw_block[3], raw_block[5]];
    for (&c, &r) in compressed_block.iter().zip(raw_encoded.iter()) {
        if swap_host {
            b.swap(c, r);
        } else {
            b.cx(c, r);
        }
    }
    let core = dialog_gcd_k2_pair_core(raw_block);
    emit_dialog_gcd_k2_pair_core_encoder_inverse(b, &core);
}

pub(crate) fn dialog_gcd_k2_pair_clear_raw_block_copy(
    b: &mut B,
    compressed_block: &[QubitId],
    raw_block: &[QubitId],
    steps: usize,
) {
    assert_eq!(compressed_block.len(), DIALOG_GCD_HIGH_TAIL_ALIAS_BLOCK_BITS);
    assert_eq!(raw_block.len(), 6);
    assert!((1..=2).contains(&steps));
    let swap_host = dialog_gcd_apply_replay_swap_host_enabled();
    if steps == 1 {
        let raw_encoded = [raw_block[0], raw_block[1], raw_block[4]];
        for (&c, &r) in compressed_block.iter().take(3).zip(raw_encoded.iter()) {
            if swap_host {
                b.swap(c, r);
            } else {
                b.cx(c, r);
            }
        }
        return;
    }
    let core = dialog_gcd_k2_pair_core(raw_block);
    emit_dialog_gcd_k2_pair_core_encoder(b, &core);
    let raw_encoded = [raw_block[1], raw_block[4], raw_block[2], raw_block[3], raw_block[5]];
    for (&c, &r) in compressed_block.iter().zip(raw_encoded.iter()) {
        if swap_host {
            b.swap(c, r);
        } else {
            b.cx(c, r);
        }
    }
}

pub(crate) fn dialog_gcd_fused_double_y(b: &mut B, y: &[QubitId], p: U256, s2: QubitId) {
    let n = y.len();
    debug_assert_eq!(n, 256);
    let c = U256::MAX.wrapping_sub(p).wrapping_add(U256::from(1));

    // ── shift1 (unconditional left shift): ovf1 = old y[255]; y[0] = 0 ──
    let ovf1 = b.alloc_qubit();
    b.swap(y[n - 1], ovf1);
    for i in (0..n - 1).rev() {
        b.swap(y[i], y[i + 1]);
    }

    // ── cond-shift2 (left shift gated by s2) on the UNFOLDED register ──
    // ovf2 = s2 & top(Y0); y[0] = 0 (and y[1] = 0 iff s2, used by cleanup).
    let ovf2 = b.alloc_qubit();
    cswap(b, s2, y[n - 1], ovf2);
    for i in (0..n - 1).rev() {
        cswap(b, s2, y[i], y[i + 1]);
    }

    // ── derive the fold controls e, d, h (3 CCX) + free-CX combinations ──
    let e = b.alloc_qubit();
    let d = b.alloc_qubit();
    let h = b.alloc_qubit();
    // d = ovf1 & s2
    b.ccx(ovf1, s2, d);
    // e = (ovf1 & ¬s2) ^ ovf2   (the two terms are mutually exclusive)
    b.x(s2);
    b.ccx(ovf1, s2, e);
    b.x(s2);
    b.cx(ovf2, e);
    // h = e & d  (= ovf2 & d, since e == ovf2 whenever d == 1)
    b.ccx(ovf2, d, h);
    // free-CX derived controls
    let xed = b.alloc_qubit(); // e ^ d
    b.cx(e, xed);
    b.cx(d, xed);
    let eord = b.alloc_qubit(); // e | d  = (e^d) ^ (e&d)
    b.cx(xed, eord);
    b.cx(h, eord);
    let n10 = b.alloc_qubit(); // ¬e & d = d ^ (e&d)
    b.cx(d, n10);
    b.cx(h, n10);

    // ── combined fold: y += δ = c·e + 2c·d, one truncated ripple ──
    // Per-position controls of the constant whose bits are
    //   c·e bits {0,4,6,7,8,9,32} and 2c·d bits {1,5,7,8,9,10,33}:
    //   0:e 1:d 4:e 5:d 6:e 7:(e^d) 8:(e|d) 9:(e|d) 10:(¬e&d) 11:(e&d) 32:e 33:d
    let hi_delta = highest_set_bit(c) + 1; // = 33 for secp256k1
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
    controls[highest_set_bit(c)] = Some(e); // bit 32
    controls[hi_delta] = Some(d); // bit 33
    let last = match double_carry_trunc_window() {
        Some(w) => core::cmp::min(n - 2, hi_delta.saturating_add(w)),
        None => n - 2,
    };
    cadd_per_position_controls_trunc(b, y, &controls, last);

    // ── cleanup: return all 8 ancilla to |0⟩ (deterministic, phase-free) ──
    // After the fold y[0] = e and (iff s2) y[1] = ovf1 (the second clean low bit).
    // Uncompute derived controls first (reverse free CX, while e,d,h still hold).
    b.cx(h, n10);
    b.cx(d, n10);
    b.cx(h, eord);
    b.cx(xed, eord);
    b.cx(d, xed);
    b.cx(e, xed);
    b.free(n10);
    b.free(eord);
    b.free(xed);
    // Clear h = ovf2 & d (ovf2, d still live).
    b.ccx(ovf2, d, h);
    b.free(h);
    // Clear e via parity: y[0] == e.
    b.cx(y[0], e);
    // Clear d via the register: d == s2 & y[1].
    b.ccx(s2, y[1], d);
    b.free(d);
    b.free(e);
    // Clear ovf1 == (s2 ? y[1] : y[0]).
    b.ccx(s2, y[1], ovf1);
    b.x(s2);
    b.ccx(s2, y[0], ovf1);
    b.x(s2);
    b.free(ovf1);
    // Clear ovf2 == s2 & y[0].
    b.ccx(s2, y[0], ovf2);
    b.free(ovf2);
}

pub(crate) fn dialog_gcd_fused_halve_y(b: &mut B, y: &[QubitId], p: U256, s2: QubitId) {
    let n = y.len();
    debug_assert_eq!(n, 256);
    let c = U256::MAX.wrapping_sub(p).wrapping_add(U256::from(1));

    // ── recover the fold controls e, d, h directly from y_new ──
    let e = b.alloc_qubit();
    let d = b.alloc_qubit();
    let h = b.alloc_qubit();
    // e = y_new[0]
    b.cx(y[0], e);
    // d = s2 & y_new[1]
    b.ccx(s2, y[1], d);
    // h = e & d
    b.ccx(e, d, h);
    // free-CX derived controls (identical to the forward fold)
    let xed = b.alloc_qubit(); // e ^ d
    b.cx(e, xed);
    b.cx(d, xed);
    let eord = b.alloc_qubit(); // e | d  = (e^d) ^ (e&d)
    b.cx(xed, eord);
    b.cx(h, eord);
    let n10 = b.alloc_qubit(); // ¬e & d = d ^ (e&d)
    b.cx(d, n10);
    b.cx(h, n10);

    // ── reconstruct the overflow bits for the reverse shifts ──
    let ovf2 = b.alloc_qubit(); // ovf2 = e & s2
    let ovf1 = b.alloc_qubit(); // ovf1 = (s2 ? d : e)
    b.ccx(e, s2, ovf2);
    b.ccx(s2, d, ovf1);
    b.x(s2);
    b.ccx(s2, e, ovf1);
    b.x(s2);

    // ── combined fold inverse: y −= δ = c·e + 2c·d, one truncated ripple ──
    // Same per-position controls as the forward fused double.
    let hi_delta = highest_set_bit(c) + 1; // = 33 for secp256k1
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
    controls[highest_set_bit(c)] = Some(e); // bit 32
    controls[hi_delta] = Some(d); // bit 33
    let last = match double_carry_trunc_window() {
        Some(w) => core::cmp::min(n - 2, hi_delta.saturating_add(w)),
        None => n - 2,
    };
    csub_per_position_controls_trunc(b, y, &controls, last);

    // ── uncompute derived controls (reverse free CX, while e,d,h still hold) ──
    b.cx(h, n10);
    b.cx(d, n10);
    b.cx(h, eord);
    b.cx(xed, eord);
    b.cx(d, xed);
    b.cx(e, xed);
    b.free(n10);
    b.free(eord);
    b.free(xed);
    // Clear h = e & d (e, d still live).
    b.ccx(e, d, h);
    b.free(h);
    // Clear e and d via the live overflow qubits (the register low bits are now
    // cleared by the csub, so we cannot read them off y any more):
    //   e == (s2 ? ovf2 : ovf1);   d == (s2 ? ovf1 : 0).
    b.x(s2);
    b.ccx(s2, ovf1, e); // s2=0: e ^= ovf1
    b.x(s2);
    b.ccx(s2, ovf2, e); // s2=1: e ^= ovf2
    b.ccx(s2, ovf1, d); // s2=1: d ^= ovf1   (s2=0: d already 0)
    b.free(e);
    b.free(d);

    // ── un-cond-shift2 (right shift gated by s2), re-inserting ovf2 at top ──
    for i in 0..n - 1 {
        cswap(b, s2, y[i], y[i + 1]);
    }
    cswap(b, s2, y[n - 1], ovf2);
    // The boundary cswap already pulled the vacated top bit (0) into ovf2, so
    // ovf2 is |0> here. (A `ccx(s2, y[n-1], ovf2)` would WRONGLY re-set it to
    // s2&y[n-1] = e, dirtying the ancilla — the free's reset then masks the
    // value error but leaks global phase. So: no extra clear.)
    b.free(ovf2);

    // ── un-shift1 (unconditional right shift), re-inserting ovf1 at top ──
    for i in 0..n - 1 {
        b.swap(y[i], y[i + 1]);
    }
    b.swap(y[n - 1], ovf1);
    // The swap already pulled the vacated top bit (0) into ovf1, so ovf1 is |0>
    // here. (A `cx(y[n-1], ovf1)` would re-dirty it — see ovf2 note above.)
    b.free(ovf1);
}
