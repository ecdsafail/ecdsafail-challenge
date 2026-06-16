//! Time-travel debugger for the 64-parallel kmx simulator.
//!
//! Attach to a fully-constructed Circuit with `Debugger::attach`,
//! then drive an interactive REPL (or call the programmatic API
//! directly). The debugger parses the Circuit's `ops` log and
//! replays gate-by-gate using the same transfer functions as our
//! in-process sim, so its reported state matches what `prove_zero`
//! and `phase_assert_zero` observed at build time.
//!
//! ## Replay model
//!
//! State at cursor k = sim state AFTER ops[0..k] have been applied,
//! starting from `circ.initial_sim_state` (captured at first
//! `set_section`, i.e. post-`sim_load_*`, pre-gate).
//!
//! ## Checkpoints and deltas
//!
//! Two-tier undo:
//!   1. Full state snapshots every `checkpoint_interval` ops
//!      (default 10M) — coarse anchors for long jumps.
//!   2. A bounded ring-buffer of per-op `Delta`s spanning
//!      [cursor - `delta_log.len()`, cursor) — O(1) single-step
//!      reverse. Capacity `DELTA_LOG_CAP` (default 500K).
//!
//! `goto(k)`:
//!   - `k == cursor`: no-op.
//!   - `k > cursor`: forward-apply, pushing deltas (ring-drops
//!     oldest past capacity).
//!   - `k < cursor` and `(cursor - k) <= delta_log.len()`:
//!     pop deltas and reverse — O(cursor - k).
//!   - otherwise: restore nearest snapshot ≤ k, clear log,
//!     forward-replay to k (log repopulates on the way if the
//!     remaining forward distance ≤ cap).
//!
//! ## Breakpoints
//!
//! - `Op(usize)` — break when cursor reaches that op index.
//! - `SectionStart(String)` — break on entering that section.
//! - `PhaseNonzero` — break the first op where `sim_phase` != 0.
//! - `ROnNonzero` — break the first R-on-non-zero event.
//!
//! ## REPL commands
//!
//! ```text
//! s / step            advance 1 op
//! b / back            rewind 1 op
//! n / next            advance to end of current section
//! c / continue        run until next breakpoint or end
//! g <idx>             jump to op index
//! g @<name>           jump to start of named section
//! p q<N>              print qubit N's 64-bit mask
//! p q<N>.s<S>         print qubit N's bit in shot S
//! p b<N>              print classical bit N's 64-bit mask
//! p phase             print sim_phase (one bit per shot)
//! p reg <start> <n>   print qubits [start..start+n] as U256 shot 0
//! p section           print current section
//! where / w           print cursor, section, progress
//! ghosts / gh         pending spooky-pebble ghosts at cursor
//!                       (id, anchor, bit, 64-shot tape mask, create site).
//!                       Reconstructed by replaying the ghost-event log up
//!                       to the cursor; use it to chase ghost-discharge
//!                       desync in windowed-tape sweeps.
//! list / l            show ~10 ops around cursor
//! break op <idx>      add op-index breakpoint
//! break section <n>   add section-start breakpoint
//! break phase         break on first phase != 0
//! break r-nonzero     break on first R-on-non-zero event
//! breakpoints / bp    list breakpoints
//! clear <id>          remove breakpoint by list-index
//! q / quit            exit
//! h / help            show this help
//! ```

use std::fmt::Write as _;
use std::io::{BufRead, Write};

/// The full sim state at some cursor position.
#[derive(Clone)]
struct SimState {
    qubits: Vec<u64>,
    bits: Vec<u64>,
    phase: u64,
    cond_stack: Vec<u32>,
    hmr_counter: u64,
    r_on_nonzero_events: u64, // cumulative count
}

#[derive(Clone)]
struct Checkpoint {
    op_idx: usize,
    state: SimState,
    profile: CursorProfileState,
}

use crate::point_add::trailmix_port::circuit::Op;

#[derive(Debug, Clone)]
pub enum Breakpoint {
    Op(usize),
    SectionStart(String),
    PhaseNonzero,
    ROnNonzero,
    /// Fire when qubit `q`'s 64-shot `sim_mask` equals `expected`.
    QubitValue {
        q: u32,
        expected: u64,
    },
    /// Fire when qubit `q`'s mask differs from the last-seen value.
    /// Useful for "when did q change?" walks, forward or backward.
    QubitChange {
        q: u32,
        last: u64,
    },
}

#[derive(Debug, Clone)]
struct SectionProfile {
    name: String,
    ops: u64,
    tof: u64,
}

#[derive(Debug, Clone)]
struct ProfileCache {
    rows: Vec<SectionProfile>,
    total_ops: u64,
    total_tof: u64,
}

#[derive(Debug, Clone, Default)]
struct CursorProfileState {
    by_section: std::collections::BTreeMap<String, (u64, u64)>,
    total_ops: u64,
    total_tof: u64,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum ProfileMode {
    Whole,
    Cursor,
    Peak,
}

/// Custom dump hook callable from the debugger REPL via `dump <name> [args]`.
/// Hook receives `&Debugger` (for cursor + qubit reads + `section_marks`) and
/// arg strings, returns a string to print.
pub type DumpHook = Box<dyn Fn(&Debugger, &[&str]) -> String + Send + Sync>;

pub struct Debugger {
    ops: Vec<Option<Op>>,
    dump_hooks: std::collections::HashMap<String, DumpHook>,
    /// Global op-index of the first op in `ops` (equal to
    /// `ops_truncated` at the time of attach, or 0 if no streaming
    /// truncation). `cursor` is stored LOCAL to `ops` (`0..=ops.len()`),
    /// but user-facing APIs (`goto`, `cursor`, `num_ops`) translate
    /// to/from global indices by adding this offset so op indices
    /// printed by panic messages, `section_marks`, and `qubit_alloc_log`
    /// line up with the debugger.
    ops_start_idx: usize,
    section_marks: Vec<(usize, String)>,
    /// Per-op scope frame index (see circuit.rs `ScopeFrameLog`).
    op_scope: Vec<u32>,
    scope_frames_log: Vec<crate::point_add::trailmix_port::circuit::ScopeFrameLog>,
    /// Full alloc log from the circuit — lets `qubit_name_at(q, op)`
    /// return the tag that was current at that `op_idx`.
    qubit_alloc_log: Vec<crate::point_add::trailmix_port::circuit::QubitAllocEvent>,
    /// Parallel log for classical-bit allocs. Anonymous `alloc_bit`
    /// records a `b{N}` tag; `alloc_bit_named` / `alloc_bits_named`
    /// record user labels (with `name[i]` suffix for batches).
    bit_alloc_log: Vec<crate::point_add::trailmix_port::circuit::BitAllocEvent>,
    /// Spooky-pebble ghost create/resolve events. Replayed up to the
    /// cursor to reconstruct the pending-ghost set (`ghosts` command).
    ghost_event_log: Vec<crate::point_add::trailmix_port::circuit::GhostEvent>,
    hmr_seed: u64,

    /// Sim-state checkpoints. Populated during emission via Circuit's
    /// `live_checkpoints` (every `checkpoint_interval` ops) and copied
    /// in by `attach`. `restore_snapshot_and_replay` snaps to the
    /// nearest one and forward-steps from there. Profile fields are
    /// initialised empty — cursor-mode profile data after a
    /// snapshot-restore reflects only the steps from the checkpoint
    /// forward.
    checkpoints: Vec<Checkpoint>,

    /// Local cursor into `ops` (`0..=ops.len()`). Use `cursor()` for
    /// the global op index that matches `circuit.ops_truncated` + local.
    cursor: usize,
    state: SimState,
    cursor_profile: CursorProfileState,

    /// Ring-buffer of deltas covering ops
    /// [cursor - `delta_log.len()`, cursor). Oldest at index 0.
    delta_log: std::collections::VecDeque<Delta>,
    delta_log_cap: usize,

    breakpoints: Vec<Breakpoint>,
    whole_profile_cache: ProfileCache,
    peak_profile_cache: ProfileCache,
}

impl Debugger {
    /// Attach a debugger by MOVING the heavy Vecs out of `circ`. The
    /// caller is expected to either panic afterwards (`DEBUG_ON_FAIL` site)
    /// or treat the circuit as gutted (`ops`, `op_scope`, `qubit_alloc_log`,
    /// etc. are now empty). Avoids the multi-GB clones that the previous
    /// `&Circuit` signature forced on 100M-op traces.
    pub fn attach(circ: &mut crate::point_add::trailmix_port::circuit::Circuit) -> Self {
        // Memory accounting at attach: log Vec sizes/capacities and current
        // RSS so OOMs are debuggable. Cheap (~10 field reads); no env gate.
        {
            let stats = circ.debug_stats();
            let rss_mib = crate::point_add::trailmix_port::circuit::Circuit::current_rss_kib().unwrap_or(0) / 1024;
            eprintln!(
                "[debugger-attach] RSS={}MiB ops.len={} cap={} \
                 op_scope.len={} qubit_alloc_log.len={} \
                 scope_frames_log.len={} live_checkpoints={} \
                 peak_live_tags={} elide_deltas={} deferred_contracts={}",
                rss_mib,
                circ.ops.len(),
                circ.ops.capacity(),
                stats.op_scope,
                stats.qubit_alloc_log,
                stats.scope_frames_log,
                circ.live_checkpoints.len(),
                circ.peak_live_tags.len(),
                circ.elide_deltas_len(),
                circ.deferred_contracts_len(),
            );
        }
        // Preserve None slots so absolute op_idx → position remains
        // 1:1 inside the debugger. None is skipped at apply time
        // (no sim/phase effect from the elide).
        // MOVE rather than clone: the panic path drops `circ` after this,
        // and a clone of a 100M-op Vec doubles peak memory transiently
        // (~5 GiB), tripping the 8G systemd cap.
        let ops: Vec<Option<Op>> = std::mem::take(&mut circ.ops);
        // Choose the starting sim state for replay:
        //   * Streaming on + prior truncation → use the snapshot
        //     taken at the most recent truncation. `ops_start_idx`
        //     becomes that truncation's global op index, and the
        //     debugger's valid range is [ops_start_idx,
        //     ops_start_idx + ops.len()].
        //   * No truncation → fall back to initial_sim_state at op 0.
        let (initial, ops_start_idx) = if let Some(snap) = &circ.truncation_snapshot {
            (
                SimState {
                    qubits: snap.qubits.clone(),
                    bits: snap.bits.clone(),
                    phase: snap.phase,
                    cond_stack: snap.cond_stack.clone(),
                    hmr_counter: snap.hmr_counter,
                    r_on_nonzero_events: snap.r_on_nonzero_events,
                },
                snap.ops_truncated as usize,
            )
        } else if let Some((q, b)) = &circ.initial_sim_state {
            (
                SimState {
                    qubits: q.clone(),
                    bits: b.clone(),
                    phase: 0,
                    cond_stack: Vec::new(),
                    hmr_counter: 0,
                    r_on_nonzero_events: 0,
                },
                0,
            )
        } else {
            eprintln!(
                "[debugger] warning: no initial_sim_state \
                captured; starting from all-zero state"
            );
            (
                SimState {
                    qubits: vec![0u64; 4096],
                    bits: vec![0u64; 1024],
                    phase: 0,
                    cond_stack: Vec::new(),
                    hmr_counter: 0,
                    r_on_nonzero_events: 0,
                },
                0,
            )
        };

        let delta_log_cap = 500_000;
        let empty_profile = CursorProfileState::default();
        // Live checkpoints captured during emission. Each carries the
        // sim state AFTER ops[..op_idx] were applied. Drop any whose
        // op_idx falls outside the debugger's replayable range
        // [ops_start_idx, ops_start_idx + ops.len()].
        let max_op_idx = ops_start_idx + ops.len();
        let mut checkpoints: Vec<Checkpoint> = Vec::new();
        // Anchor at op_idx = ops_start_idx so restore_snapshot_and_replay
        // always finds at least one checkpoint <= target.
        checkpoints.push(Checkpoint {
            op_idx: 0,
            state: initial.clone(),
            profile: empty_profile.clone(),
        });
        for c in &circ.live_checkpoints {
            let op_idx = c.ops_truncated as usize;
            if op_idx <= ops_start_idx || op_idx > max_op_idx {
                continue;
            }
            checkpoints.push(Checkpoint {
                // Stored as LOCAL index (relative to ops_start_idx) to
                // match the rest of the debugger's cursor convention.
                op_idx: op_idx - ops_start_idx,
                state: SimState {
                    qubits: c.qubits.clone(),
                    bits: c.bits.clone(),
                    phase: c.phase,
                    cond_stack: c.cond_stack.clone(),
                    hmr_counter: c.hmr_counter,
                    r_on_nonzero_events: c.r_on_nonzero_events,
                },
                profile: empty_profile.clone(),
            });
        }
        // Final state at attach time: use Circuit's CURRENT sim state
        // directly. No replay required — Circuit already simulated every
        // op during emission.
        let end_state = if let Some(sim) = &circ.sim {
            SimState {
                qubits: sim.clone(),
                bits: circ.sim_bits.clone().unwrap_or_default(),
                phase: circ.sim_phase,
                cond_stack: circ.sim_condition_stack.clone(),
                hmr_counter: circ.sim_hmr_counter,
                r_on_nonzero_events: u64::from(circ.sim_phase_errors),
            }
        } else {
            initial.clone()
        };
        let cursor_local = ops.len();
        eprintln!(
            "[debugger] sizes pre-clone: ops.len={} section_marks.len={} op_scope.len={} \
             scope_frames_log.len={} qubit_alloc_log.len={} checkpoints={}",
            ops.len(),
            circ.section_marks.len(),
            circ.op_scope.len(),
            circ.scope_frames_log.len(),
            circ.qubit_alloc_log.len(),
            checkpoints.len(),
        );
        Self {
            ops,
            dump_hooks: std::mem::take(&mut circ.dump_hooks),
            ops_start_idx,
            section_marks: std::mem::take(&mut circ.section_marks),
            op_scope: std::mem::take(&mut circ.op_scope),
            scope_frames_log: std::mem::take(&mut circ.scope_frames_log),
            qubit_alloc_log: std::mem::take(&mut circ.qubit_alloc_log),
            bit_alloc_log: std::mem::take(&mut circ.bit_alloc_log),
            ghost_event_log: std::mem::take(&mut circ.ghost_event_log),
            hmr_seed: circ.initial_hmr_seed,
            checkpoints,
            cursor: cursor_local,
            state: end_state,
            cursor_profile: empty_profile,
            delta_log: std::collections::VecDeque::with_capacity(delta_log_cap.min(1 << 16)),
            delta_log_cap,
            breakpoints: Vec::new(),
            whole_profile_cache: Self::build_whole_profile_cache(circ),
            peak_profile_cache: Self::build_peak_profile_cache(circ),
        }
    }

    fn build_whole_profile_cache(circ: &crate::point_add::trailmix_port::circuit::Circuit) -> ProfileCache {
        let mut keys = std::collections::BTreeSet::<String>::new();
        keys.extend(circ.executed_ops_by_section.keys().cloned());
        keys.extend(circ.executed_toffoli_by_section.keys().cloned());
        let mut rows = Vec::new();
        let mut total_ops = 0u64;
        let mut total_tof = 0u64;
        for name in keys {
            let ops = circ
                .executed_ops_by_section
                .get(&name)
                .copied()
                .unwrap_or(0);
            let tof = circ
                .executed_toffoli_by_section
                .get(&name)
                .copied()
                .unwrap_or(0);
            total_ops += ops;
            total_tof += tof;
            rows.push(SectionProfile { name, ops, tof });
        }
        rows.sort_by(|a, b| {
            b.ops
                .cmp(&a.ops)
                .then_with(|| b.tof.cmp(&a.tof))
                .then_with(|| a.name.cmp(&b.name))
        });
        ProfileCache {
            rows,
            total_ops,
            total_tof,
        }
    }

    fn normalize_peak_tag_class(tag_class: &str) -> String {
        let mut out = String::with_capacity(tag_class.len());
        let mut chars = tag_class.chars().peekable();
        while let Some(ch) = chars.next() {
            if ch == '[' {
                out.push('[');
                for next in chars.by_ref() {
                    if next == ']' {
                        out.push(']');
                        break;
                    }
                }
            } else {
                out.push(ch);
            }
        }
        out
    }

    fn build_peak_profile_cache(circ: &crate::point_add::trailmix_port::circuit::Circuit) -> ProfileCache {
        let mut grouped = std::collections::BTreeMap::<String, u64>::new();
        for tag in &circ.peak_live_tags {
            let name = match tag.rsplit_once('/') {
                Some((prefix, tag_class)) => {
                    format!("{}/{}", prefix, Self::normalize_peak_tag_class(tag_class))
                }
                None => Self::normalize_peak_tag_class(tag),
            };
            *grouped.entry(name).or_insert(0) += 1;
        }
        let mut rows = grouped
            .into_iter()
            .map(|(name, ops)| SectionProfile { name, ops, tof: 0 })
            .collect::<Vec<_>>();
        rows.sort_by(|a, b| b.ops.cmp(&a.ops).then_with(|| a.name.cmp(&b.name)));
        let total_ops = rows.iter().map(|r| r.ops).sum();
        ProfileCache {
            rows,
            total_ops,
            total_tof: 0,
        }
    }

    /// Global op index of the debugger's current position (adds
    /// `ops_start_idx` to the local cursor so the value matches
    /// `Circuit::total_ops()` semantics even when streaming has
    /// truncated earlier ops).
    #[must_use]
    pub fn cursor(&self) -> usize {
        self.cursor + self.ops_start_idx
    }
    pub fn section_marks_iter(&self) -> impl Iterator<Item = (usize, &str)> {
        self.section_marks.iter().map(|(i, n)| (*i, n.as_str()))
    }
    /// Global op index of the last op in `ops` (== `circuit.total_ops`
    /// at attach time).
    #[must_use]
    pub fn num_ops(&self) -> usize {
        self.ops.len() + self.ops_start_idx
    }
    /// Smallest op index the debugger can reach (= `ops_start_idx`).
    /// Earlier ops were dropped by `CIRC_OPS_CAP` streaming.
    #[must_use]
    pub fn min_op(&self) -> usize {
        self.ops_start_idx
    }

    #[must_use]
    pub fn current_section(&self) -> &str {
        let mut last = "";
        let gc = self.cursor();
        for (idx, name) in &self.section_marks {
            if *idx <= gc {
                last = name.as_str();
            } else {
                break;
            }
        }
        last
    }

    /// The spooky-pebble ghosts that are PENDING (created but not yet
    /// resolved) as of the current cursor position. Reconstructed by
    /// replaying `ghost_event_log` up to `cursor()`: a create adds the
    /// ghost, a resolve removes it. Returns the create-events still
    /// pending, in creation order.
    #[must_use]
    pub fn pending_ghosts_at_cursor(&self) -> Vec<&crate::point_add::trailmix_port::circuit::GhostEvent> {
        let gc = self.cursor();
        let mut pending: Vec<&crate::point_add::trailmix_port::circuit::GhostEvent> = Vec::new();
        for ev in &self.ghost_event_log {
            if ev.op_idx > gc {
                break; // log is in op order; nothing later applies yet
            }
            if ev.create {
                pending.push(ev);
            } else {
                pending.retain(|c| c.id != ev.id);
            }
        }
        pending
    }

    /// Human-readable summary of the pending-ghost set at the cursor,
    /// for the `ghosts` REPL command.
    #[must_use]
    pub fn ghosts_line(&self) -> String {
        let pending = self.pending_ghosts_at_cursor();
        if self.ghost_event_log.is_empty() {
            return "ghosts: (none — circuit emitted no hmr_ghost)".to_string();
        }
        if pending.is_empty() {
            return format!(
                "ghosts: 0 pending at op {} ({} create/resolve events total)",
                self.cursor(),
                self.ghost_event_log.len(),
            );
        }
        let mut s = format!(
            "ghosts: {} pending at op {} —\n",
            pending.len(),
            self.cursor(),
        );
        for g in &pending {
            let _ = writeln!(
                s,
                "  #{:<4} anchor={:<4} bit=b{:<5} tape_mask={:#018x} \
                 created@op {} '{}'",
                g.id, g.anchor_id, g.bit_raw, g.mask_at_hmr, g.op_idx, g.section,
            );
        }
        s.pop(); // trailing newline
        s
    }

    #[must_use]
    pub fn qubit(&self, q: u32) -> u64 {
        self.state.qubits.get(q as usize).copied().unwrap_or(0)
    }
    #[must_use]
    pub fn bit(&self, b: u32) -> u64 {
        self.state.bits.get(b as usize).copied().unwrap_or(0)
    }
    #[must_use]
    pub fn phase(&self) -> u64 {
        self.state.phase
    }
    #[must_use]
    pub fn r_events(&self) -> u64 {
        self.state.r_on_nonzero_events
    }

    /// Forward-step one op, appending its delta to the log (and
    /// ring-dropping the oldest entry if at capacity).
    fn forward_one(&mut self) {
        debug_assert!(self.cursor < self.ops.len());
        let global_idx = self.ops_start_idx + self.cursor;
        let section = self.section_name_at_global(global_idx);
        if let Some(op) = self.ops[self.cursor] {
            // tof counts STATIC Toffoli emissions (see whole-pass
            // profiler): 1 per CCX/CCZ that fires on at least one
            // shot, 0 otherwise.
            let tof = match &op {
                Op::Ccx(_, _, _) | Op::Ccz(_, _, _) if cond_mask(&self.state) != 0 => 1,
                _ => 0,
            };
            let entry = self
                .cursor_profile
                .by_section
                .entry(section)
                .or_insert((0, 0));
            entry.0 += 1;
            entry.1 += tof;
            self.cursor_profile.total_ops += 1;
            self.cursor_profile.total_tof += tof;
            let d = apply_op_recorded(&op, &mut self.state, self.hmr_seed);
            if self.delta_log.len() == self.delta_log_cap {
                self.delta_log.pop_front();
            }
            self.delta_log.push_back(d);
        } else {
            // Elided slot (None): no sim/phase work; record a no-op
            // delta so back_one_via_log balances correctly.
            if self.delta_log.len() == self.delta_log_cap {
                self.delta_log.pop_front();
            }
            self.delta_log.push_back(Delta::None);
        }
        self.cursor += 1;
    }

    /// Backward-step one op via the delta log (O(1)).
    /// Caller must ensure the log is non-empty.
    fn back_one_via_log(&mut self) {
        debug_assert!(!self.delta_log.is_empty());
        debug_assert!(self.cursor > 0);
        let d = self.delta_log.pop_back().expect("log non-empty");
        if let Some(op) = self.ops[self.cursor - 1].as_ref() {
            reverse_op(op, &d, &mut self.state, self.hmr_seed);
            let global_idx = self.ops_start_idx + self.cursor - 1;
            let section = self.section_name_at_global(global_idx);
            // Match the static-tof convention used in forward_one /
            // bulk_replay: 1 per CCX/CCZ that fires on at least one
            // shot, 0 otherwise.
            let tof = match op {
                Op::Ccx(_, _, _) | Op::Ccz(_, _, _) if cond_mask(&self.state) != 0 => 1,
                _ => 0,
            };
            let remove_entry = if let Some(entry) = self.cursor_profile.by_section.get_mut(&section)
            {
                entry.0 = entry.0.saturating_sub(1);
                entry.1 = entry.1.saturating_sub(tof);
                entry.0 == 0 && entry.1 == 0
            } else {
                false
            };
            if remove_entry {
                self.cursor_profile.by_section.remove(&section);
            }
            self.cursor_profile.total_ops = self.cursor_profile.total_ops.saturating_sub(1);
            self.cursor_profile.total_tof = self.cursor_profile.total_tof.saturating_sub(tof);
        }
        // Elided slot: nothing to reverse, profile already wasn't bumped.
        self.cursor -= 1;
    }

    /// Rebuild state at `target` from the nearest snapshot, clearing
    /// the delta log. Used when the target is outside the current
    /// delta window.
    fn restore_snapshot_and_replay(&mut self, target: usize) {
        let target = target.min(self.ops.len());
        let cp = self
            .checkpoints
            .iter()
            .rposition(|c| c.op_idx <= target)
            .map(|i| &self.checkpoints[i])
            .expect("checkpoint chain always has op_idx=0");
        self.state = cp.state.clone();
        self.cursor_profile = cp.profile.clone();
        self.cursor = cp.op_idx;
        self.delta_log.clear();
        while self.cursor < target {
            self.forward_one();
        }
    }

    /// Set cursor to global op index `target_global`. Clamps to
    /// `[ops_start_idx, ops_start_idx + ops.len()]` (the valid
    /// replayable range). Uses the delta log for cheap local travel
    /// and snapshot+replay for long jumps.
    pub fn goto(&mut self, target_global: usize) {
        let max_global = self.ops_start_idx + self.ops.len();
        let tg = target_global.min(max_global).max(self.ops_start_idx);
        let target = tg - self.ops_start_idx; // local
        if target == self.cursor {
            return;
        }
        if target > self.cursor {
            // Forward. If the gap is enormous, hopping to a snapshot
            // nearer the target is cheaper than serial forward-step.
            let cp_ahead = self
                .checkpoints
                .iter()
                .rposition(|c| c.op_idx <= target && c.op_idx > self.cursor);
            if let Some(i) = cp_ahead {
                let cp = &self.checkpoints[i];
                if cp.op_idx - self.cursor > self.delta_log_cap {
                    self.state = cp.state.clone();
                    self.cursor = cp.op_idx;
                    self.delta_log.clear();
                }
            }
            while self.cursor < target {
                self.forward_one();
            }
            return;
        }
        // Backward.
        let distance = self.cursor - target;
        if distance <= self.delta_log.len() {
            for _ in 0..distance {
                self.back_one_via_log();
            }
        } else {
            self.restore_snapshot_and_replay(target);
        }
    }

    pub fn step(&mut self, n: usize) {
        self.goto(self.cursor().saturating_add(n));
    }
    pub fn back(&mut self, n: usize) {
        self.goto(self.cursor().saturating_sub(n));
    }

    /// Advance until next section boundary or end.
    pub fn next_section(&mut self) {
        let current = self.current_section().to_string();
        let gc = self.cursor();
        let end = self.num_ops();
        let target = self
            .section_marks
            .iter()
            .find(|(idx, _)| *idx > gc)
            .map_or(end, |(idx, _)| *idx);
        let _ = current;
        self.goto(target);
    }

    /// Advance until a breakpoint condition fires or end reached.
    /// Returns the triggering breakpoint if any.
    pub fn run_until_break(&mut self) -> Option<Breakpoint> {
        while self.cursor < self.ops.len() {
            let next_global = self.cursor + 1 + self.ops_start_idx;
            let op_bp = self
                .breakpoints
                .iter()
                .find(|bp| matches!(bp, Breakpoint::Op(i) if *i == next_global))
                .cloned();
            if let Some(bp) = op_bp {
                self.goto(next_global);
                return Some(bp);
            }
            let sec_bp = self
                .breakpoints
                .iter()
                .find(|bp| match bp {
                    Breakpoint::SectionStart(name) => self
                        .section_marks
                        .iter()
                        .any(|(idx, n)| *idx == next_global && n == name),
                    _ => false,
                })
                .cloned();
            if let Some(bp) = sec_bp {
                self.goto(next_global);
                return Some(bp);
            }
            let had_phase = self.state.phase != 0;
            let had_r = self.state.r_on_nonzero_events;
            // Snapshot qubit-value watches BEFORE stepping.
            let qv_before: Vec<(usize, u64)> = self
                .breakpoints
                .iter()
                .enumerate()
                .filter_map(|(i, bp)| match bp {
                    Breakpoint::QubitValue { q, .. } | Breakpoint::QubitChange { q, .. } => {
                        Some((i, self.qubit(*q)))
                    }
                    _ => None,
                })
                .collect();
            self.forward_one();
            if !had_phase && self.state.phase != 0 {
                for bp in &self.breakpoints {
                    if matches!(bp, Breakpoint::PhaseNonzero) {
                        return Some(bp.clone());
                    }
                }
            }
            if self.state.r_on_nonzero_events > had_r {
                for bp in &self.breakpoints {
                    if matches!(bp, Breakpoint::ROnNonzero) {
                        return Some(bp.clone());
                    }
                }
            }
            // Check qubit-value / qubit-change breakpoints.
            for (i, before) in &qv_before {
                match &self.breakpoints[*i] {
                    Breakpoint::QubitValue { q, expected } => {
                        let now = self.qubit(*q);
                        if now == *expected && *before != *expected {
                            return Some(self.breakpoints[*i].clone());
                        }
                    }
                    Breakpoint::QubitChange { q, .. } => {
                        let now = self.qubit(*q);
                        if now != *before {
                            return Some(Breakpoint::QubitChange {
                                q: *q,
                                last: *before,
                            });
                        }
                    }
                    _ => {}
                }
            }
        }
        None
    }

    /// Walk backward through ops until any active breakpoint fires
    /// (`QubitValue` / `QubitChange` / Op). Useful for "when did q change
    /// from 0 to 1?" style investigations. Returns the firing bp or
    /// None if we hit op 0 without matching.
    pub fn run_backward_until_break(&mut self) -> Option<Breakpoint> {
        while self.cursor > 0 {
            // Check Op breakpoints first (cursor-matching).
            let gc = self.cursor();
            let op_bp = self
                .breakpoints
                .iter()
                .find(|bp| matches!(bp, Breakpoint::Op(i) if *i == gc))
                .cloned();
            if let Some(bp) = op_bp {
                return Some(bp);
            }
            // Snapshot qubit state BEFORE stepping back.
            let qv_before: Vec<(usize, u64)> = self
                .breakpoints
                .iter()
                .enumerate()
                .filter_map(|(i, bp)| match bp {
                    Breakpoint::QubitValue { q, .. } | Breakpoint::QubitChange { q, .. } => {
                        Some((i, self.qubit(*q)))
                    }
                    _ => None,
                })
                .collect();
            self.back(1);
            for (i, before) in &qv_before {
                match &self.breakpoints[*i] {
                    Breakpoint::QubitValue { q, expected } => {
                        let now = self.qubit(*q);
                        if now == *expected && *before != *expected {
                            return Some(self.breakpoints[*i].clone());
                        }
                    }
                    Breakpoint::QubitChange { q, .. } => {
                        let now = self.qubit(*q);
                        if now != *before {
                            return Some(Breakpoint::QubitChange {
                                q: *q,
                                last: *before,
                            });
                        }
                    }
                    _ => {}
                }
            }
        }
        None
    }

    pub fn add_breakpoint(&mut self, bp: Breakpoint) {
        self.breakpoints.push(bp);
    }
    pub fn clear_breakpoint(&mut self, idx: usize) {
        if idx < self.breakpoints.len() {
            self.breakpoints.remove(idx);
        }
    }
    #[must_use]
    pub fn breakpoints(&self) -> &[Breakpoint] {
        &self.breakpoints
    }

    #[must_use]
    pub fn op_text(&self, idx: usize) -> Option<String> {
        self.ops
            .get(idx)
            .and_then(|s| s.as_ref())
            .map(|op| self.format_op_named(idx, op))
    }

    /// Name-resolving formatter for ops. Replaces raw qubit ids with
    /// the tag they had at `op_idx`. Falls back to `q{N}` only if
    /// nothing was ever recorded for the id.
    #[must_use]
    pub fn format_op_named(&self, op_idx: usize, op: &Op) -> String {
        let q = |n: u32| self.qubit_name_at(n, op_idx);
        match *op {
            Op::Register(r) => format!("REGISTER {r}"),
            Op::AppendQubit(r, qb) => format!("APPEND_Q {} {}", r, q(qb)),
            Op::AppendBit(r, b) => format!("APPEND_B {r} b{b}"),
            Op::X(qb) => format!("X {}", q(qb)),
            Op::Z(qb) => format!("Z {}", q(qb)),
            Op::Cx(c, t) => format!("CX {} {}", q(c), q(t)),
            Op::Cz(a, b) => format!("CZ {} {}", q(a), q(b)),
            Op::Ccx(a, b, c) => format!("CCX {} {} {}", q(a), q(b), q(c)),
            Op::Ccz(a, b, c) => format!("CCZ {} {} {}", q(a), q(b), q(c)),
            Op::Swap(a, b) => format!("SWAP {} {}", q(a), q(b)),
            Op::Hmr(qb, b) => format!("HMR {} b{}", q(qb), b),
            Op::R(qb) => format!("R {}", q(qb)),
            Op::Neg => "NEG".into(),
            Op::PushCondition(b) => format!("PUSH_COND b{b}"),
            Op::PopCondition => "POP_COND".into(),
            Op::BitInvert(b) => format!("BIT_INV b{b}"),
            Op::BitStore0(b) => format!("BIT_STORE0 b{b}"),
            Op::BitStore1(b) => format!("BIT_STORE1 b{b}"),
        }
    }

    /// Return the tag for qubit `qid` as of `op_idx`. Walks the alloc
    /// log (sorted by `op_idx`) for the most recent alloc at or before
    /// `op_idx`. Falls back to `q{N}` if no entry exists.
    #[must_use]
    pub fn qubit_name_at(&self, qid: u32, op_idx: usize) -> String {
        let mut best: Option<&str> = None;
        let mut best_op: usize = 0;
        for ev in &self.qubit_alloc_log {
            if ev.qubit == qid && ev.op_idx <= op_idx && (best.is_none() || ev.op_idx >= best_op) {
                best = Some(&*ev.tag);
                best_op = ev.op_idx;
            }
        }
        best.map_or_else(|| format!("q{qid}"), std::string::ToString::to_string)
    }

    /// Resolve a tag (full or suffix match) to the most recent qubit
    /// id allocated with that tag at or before the cursor. Matches on
    /// exact string first; falls back to any tag ENDING with `tag`
    /// (handy for dropping scope prefixes, e.g. `det_sign` matches
    /// `pornin.outer_12/det_sign`).
    #[must_use]
    pub fn qubit_id_from_tag(&self, tag: &str) -> Option<u32> {
        let mut best: Option<u32> = None;
        let mut best_op: usize = 0;
        for ev in &self.qubit_alloc_log {
            if ev.op_idx > self.cursor() {
                continue;
            }
            if &*ev.tag == tag && (best.is_none() || ev.op_idx >= best_op) {
                best = Some(ev.qubit);
                best_op = ev.op_idx;
            }
        }
        if best.is_some() {
            return best;
        }
        // Suffix fallback.
        for ev in &self.qubit_alloc_log {
            if ev.op_idx > self.cursor() {
                continue;
            }
            if ev.tag.ends_with(tag) && (best.is_none() || ev.op_idx >= best_op) {
                best = Some(ev.qubit);
                best_op = ev.op_idx;
            }
        }
        best
    }

    /// Cbit analogue of `qubit_id_from_tag`. Walks `bit_alloc_log`
    /// for an exact tag match (most recent at or before the cursor),
    /// then falls back to a tag-ends-with suffix match.
    #[must_use]
    pub fn cbit_id_from_tag(&self, tag: &str) -> Option<u32> {
        let mut best: Option<u32> = None;
        let mut best_op: usize = 0;
        for ev in &self.bit_alloc_log {
            if ev.op_idx > self.cursor() {
                continue;
            }
            if &*ev.tag == tag && (best.is_none() || ev.op_idx >= best_op) {
                best = Some(ev.bit);
                best_op = ev.op_idx;
            }
        }
        if best.is_some() {
            return best;
        }
        for ev in &self.bit_alloc_log {
            if ev.op_idx > self.cursor() {
                continue;
            }
            if ev.tag.ends_with(tag) && (best.is_none() || ev.op_idx >= best_op) {
                best = Some(ev.bit);
                best_op = ev.op_idx;
            }
        }
        best
    }

    /// Cbit analogue of `qubit_name_at`. Returns the most recent tag
    /// for cbit `bid` at or before `op_idx`, falling back to `b{N}`.
    #[must_use]
    pub fn cbit_name_at(&self, bid: u32, op_idx: usize) -> String {
        let mut best: Option<&str> = None;
        let mut best_op: usize = 0;
        for ev in &self.bit_alloc_log {
            if ev.bit == bid && ev.op_idx <= op_idx && (best.is_none() || ev.op_idx >= best_op) {
                best = Some(&*ev.tag);
                best_op = ev.op_idx;
            }
        }
        best.map_or_else(|| format!("b{bid}"), std::string::ToString::to_string)
    }

    /// Dump a window of ops around the cursor, annotating section
    /// marks and the cursor position.
    #[must_use]
    pub fn list(&self, ctx: usize) -> String {
        // Print ops in the GLOBAL index space. Local `i` iterates
        // over retained ops; `gi = i + ops_start_idx` is what shows
        // in the output (matching panic-site op indices).
        let lo = self.cursor.saturating_sub(ctx);
        let hi = (self.cursor + ctx).min(self.ops.len());
        let gc = self.cursor();
        let mut out = String::new();
        for i in lo..hi {
            let gi = i + self.ops_start_idx;
            for (idx, name) in &self.section_marks {
                if *idx == gi {
                    let _ = writeln!(out, "                  ┌── [{name}]");
                }
            }
            let marker = if gi == gc { ">>" } else { "  " };
            let op_str = self.ops[i]
                .as_ref()
                .map_or_else(|| "(elided)".into(), |op| self.format_op_named(gi, op));
            let _ = writeln!(out, "  {marker} {gi:>9}  {op_str}");
        }
        if self.cursor == self.ops.len() {
            let _ = writeln!(out, "  >> {gc:>9}  (end)");
        }
        out
    }

    /// Read a register of qubits as a little-endian U256 for shot 0.
    /// Takes raw qubit ids (u32) — `Qubit` is module-private to
    /// `circuit.rs` and cannot be reconstructed here.
    #[must_use]
    pub fn read_reg_shot(&self, qs: &[u32], shot: usize) -> Vec<u8> {
        let n = qs.len();
        let mut bytes = vec![0u8; n.div_ceil(8)];
        for (i, &q) in qs.iter().enumerate() {
            if (self.qubit(q) >> shot) & 1 == 1 {
                bytes[i / 8] |= 1 << (i % 8);
            }
        }
        bytes
    }

    /// Read a register addressed by per-bit TAGS `label[0]`, `label[1]`, ...
    /// at the current cursor, as an integer for `shot`. Follows registers
    /// that grow/shrink (dynamic-W): each `label[k]` resolves to whatever
    /// qubit currently carries that tag. Stops at the first `k` whose tag is
    /// no longer live (resolved qubit's current name doesn't match), capped
    /// at `max_bits`. Freed-but-unreused high bits read as 0, so the value is
    /// still correct. Lets a dump hook print resizing registers algebraically.
    #[must_use]
    pub fn read_tagged_reg(
        &self,
        label: &str,
        shot: usize,
        max_bits: usize,
    ) -> num_bigint::BigUint {
        let cursor = self.cursor();
        let mut val = num_bigint::BigUint::from(0u32);
        for k in 0..max_bits {
            let tag = format!("{label}[{k}]");
            let Some(qid) = self.qubit_id_from_tag(&tag) else {
                break;
            };
            if !self.qubit_name_at(qid, cursor).ends_with(&tag) {
                break;
            }
            if (self.qubit(qid) >> shot) & 1 == 1 {
                val |= num_bigint::BigUint::from(1u32) << k;
            }
        }
        val
    }

    /// Enter an interactive REPL over stdin/stdout.
    pub fn repl(&mut self) {
        let stdin = std::io::stdin();
        let stdout = std::io::stdout();
        let mut out = stdout.lock();
        writeln!(
            out,
            "[debugger] attached — {} ops in {} sections. \
            Type `h` for help.",
            self.ops.len(),
            self.section_marks.len()
        )
        .ok();
        writeln!(out, "{}", self.where_line()).ok();

        for line in stdin.lock().lines() {
            let Ok(line) = line else {
                break;
            };
            let line = line.trim();
            if line.is_empty() {
                write!(out, "(dbg) ").ok();
                out.flush().ok();
                continue;
            }
            match self.run_command(line) {
                ReplAction::Quit => break,
                ReplAction::Print(s) => {
                    writeln!(out, "{s}").ok();
                }
            }
            write!(out, "(dbg) ").ok();
            out.flush().ok();
        }
    }

    fn where_line(&self) -> String {
        let mut s = format!(
            "  op={}/{}  section={}  phase={:#018x}  r-nonzero={}  range=[{}..{})",
            self.cursor(),
            self.num_ops(),
            self.current_section(),
            self.state.phase,
            self.state.r_on_nonzero_events,
            self.ops_start_idx,
            self.num_ops(),
        );
        // Scope chain from the cursor's current op (= last retained if
        // streaming truncated). Cheaper than `src` since we only need the
        // names not the file:line.
        if !self.op_scope.is_empty() {
            let local = self.cursor.min(self.op_scope.len() - 1);
            let mut cur = self.op_scope[local];
            if cur != u32::MAX {
                s.push_str("\n  scope:");
                let mut names: Vec<&str> = Vec::new();
                while cur != u32::MAX {
                    let f = &self.scope_frames_log[cur as usize];
                    names.push(&f.name);
                    cur = f.parent.unwrap_or(u32::MAX);
                }
                // Print root → innermost (= reverse of walk order).
                for (i, name) in names.iter().rev().enumerate() {
                    if i > 0 {
                        s.push_str(" / ");
                    } else {
                        s.push(' ');
                    }
                    s.push_str(name);
                }
            }
        }
        s
    }

    /// Render the scope chain (innermost → root) for a given GLOBAL op
    /// index, each frame on its own line. Used by the `src` REPL command.
    /// `op_scope` is LOCAL-indexed (length = `self.ops.len()`); convert global
    /// → local. If the requested op is below the retained tail (= dropped
    /// by streaming truncation), fall back to the FIRST retained op's
    /// scope (the scope is unlikely to have changed for a streaming-context
    /// query) and tag the output so the caller knows.
    fn source_trace(&self, op_idx_global: usize) -> String {
        if self.op_scope.is_empty() {
            return format!(
                "op {op_idx_global}: no retained ops (streaming truncation cleared op_scope; \
                 enable a smaller CIRC_OPS_CAP or rerun with debug-on-fail \
                 earlier in the trace)"
            );
        }
        let (local, truncated_note) = if op_idx_global < self.ops_start_idx {
            (
                0usize,
                format!(
                    " [note: op {} is below retained tail (start={}); showing scope \
                 at first retained op {}]",
                    op_idx_global, self.ops_start_idx, self.ops_start_idx,
                ),
            )
        } else {
            let off = op_idx_global - self.ops_start_idx;
            if off >= self.op_scope.len() {
                let last = self.op_scope.len() - 1;
                (
                    last,
                    format!(
                        " [note: op {} is past last retained op {}; showing scope \
                     at last retained op]",
                        op_idx_global,
                        self.ops_start_idx + last,
                    ),
                )
            } else {
                (off, String::new())
            }
        };
        let mut cur = self.op_scope[local];
        if cur == u32::MAX {
            return format!(
                "op {op_idx_global}: no scope frame at local idx {local} (use enter_scope! to tag){truncated_note}",
            );
        }
        let mut out = format!(
            "op {op_idx_global} (local {local}): scope chain (innermost → root):{truncated_note}\n",
        );
        while cur != u32::MAX {
            let f = &self.scope_frames_log[cur as usize];
            let _ = writeln!(out, "  {} at {}:{}", f.name, f.file, f.line);
            cur = f.parent.unwrap_or(u32::MAX);
        }
        out
    }

    fn section_name_at_global(&self, global_idx: usize) -> String {
        // section_marks is sorted by op_idx (push/pop append in emission
        // order). Binary search for the largest entry with idx <=
        // global_idx via partition_point — O(log M) vs the prior linear
        // scan that turned build_checkpoints into O(N × M).
        let pos = self
            .section_marks
            .partition_point(|(idx, _)| *idx <= global_idx);
        if pos == 0 {
            return String::new();
        }
        self.section_marks[pos - 1].1.clone()
    }

    fn cursor_profile_rows(&self) -> Vec<SectionProfile> {
        let mut rows: Vec<SectionProfile> = self
            .cursor_profile
            .by_section
            .iter()
            .map(|(name, (ops, tof))| SectionProfile {
                name: name.clone(),
                ops: *ops,
                tof: *tof,
            })
            .collect();
        rows.sort_by(|a, b| {
            b.ops
                .cmp(&a.ops)
                .then_with(|| b.tof.cmp(&a.tof))
                .then_with(|| a.name.cmp(&b.name))
        });
        rows
    }

    fn profile_source(&self, mode: ProfileMode) -> (Vec<SectionProfile>, u64, u64) {
        match mode {
            ProfileMode::Whole => (
                self.whole_profile_cache.rows.clone(),
                self.whole_profile_cache.total_ops,
                self.whole_profile_cache.total_tof,
            ),
            ProfileMode::Cursor => (
                self.cursor_profile_rows(),
                self.cursor_profile.total_ops,
                self.cursor_profile.total_tof,
            ),
            ProfileMode::Peak => (
                self.peak_profile_cache.rows.clone(),
                self.peak_profile_cache.total_ops,
                0,
            ),
        }
    }

    fn format_profile_rows<'a, I>(
        rows: I,
        total_ops: u64,
        total_tof: u64,
        limit: usize,
        sort_by_tof: bool,
        peak_mode: bool,
    ) -> String
    where
        I: IntoIterator<Item = &'a SectionProfile>,
    {
        let mut rows: Vec<&SectionProfile> = rows.into_iter().collect();
        if sort_by_tof {
            rows.sort_by(|a, b| {
                b.tof
                    .cmp(&a.tof)
                    .then_with(|| b.ops.cmp(&a.ops))
                    .then_with(|| a.name.cmp(&b.name))
            });
        } else {
            rows.sort_by(|a, b| {
                b.ops
                    .cmp(&a.ops)
                    .then_with(|| b.tof.cmp(&a.tof))
                    .then_with(|| a.name.cmp(&b.name))
            });
        }

        let matched_ops: u64 = rows.iter().map(|r| r.ops).sum();
        let matched_tof: u64 = rows.iter().map(|r| r.tof).sum();
        let shown = rows.len().min(limit);

        let mut out = String::new();
        if peak_mode {
            let _ = writeln!(
                out,
                "matched={} shown={} total_qubits={} matched_qubits={}",
                rows.len(),
                shown,
                total_ops,
                matched_ops,
            );
            out.push_str("   qubits     q%  section\n");
            for row in rows.into_iter().take(limit) {
                let q_pct = if total_ops == 0 {
                    0.0
                } else {
                    100.0 * (row.ops as f64) / (total_ops as f64)
                };
                let _ = writeln!(out, "{:>9} {:>6.1}  {}", row.ops, q_pct, row.name);
            }
        } else {
            // tof is summed across 64 shots (= total Toffoli fires); display
            // per-shot average by dividing by 64.
            const SHOTS: u64 = 64;
            let _ = writeln!(
                out,
                "matched={} shown={} total_ops={} total_tof={} matched_ops={} matched_tof={}",
                rows.len(),
                shown,
                total_ops,
                total_tof / SHOTS,
                matched_ops,
                matched_tof / SHOTS,
            );
            out.push_str("      ops    tof/shot   op%   tof%  section\n");
            for row in rows.into_iter().take(limit) {
                let op_pct = if total_ops == 0 {
                    0.0
                } else {
                    100.0 * (row.ops as f64) / (total_ops as f64)
                };
                let tof_pct = if total_tof == 0 {
                    0.0
                } else {
                    100.0 * (row.tof as f64) / (total_tof as f64)
                };
                let _ = writeln!(
                    out,
                    "{:>9} {:>10} {:>5.1} {:>6.1}  {}",
                    row.ops,
                    row.tof / SHOTS,
                    op_pct,
                    tof_pct,
                    row.name
                );
            }
        }
        out
    }

    fn print_profile(&mut self, args: &[&str]) -> String {
        const USAGE: &str =
            "usage: prof [whole|cursor|peak] [ops|tof] {top [n] | exact <name> [n] | prefix <prefix> [n] | contains <substr> [n] | current [n] | split <prefix|current> [n]}";

        let mut mode = ProfileMode::Whole;
        let mut sort_by_tof = false;
        let mut idx = 0usize;
        if let Some(&scope) = args.get(idx) {
            match scope {
                "whole" => {
                    mode = ProfileMode::Whole;
                    idx += 1;
                }
                "cursor" => {
                    mode = ProfileMode::Cursor;
                    idx += 1;
                }
                "peak" => {
                    mode = ProfileMode::Peak;
                    idx += 1;
                }
                _ => {}
            }
        }
        if let Some(&metric) = args.get(idx) {
            match metric {
                "ops" => idx += 1,
                "tof" => {
                    sort_by_tof = true;
                    idx += 1;
                }
                "q" | "qubits" if mode == ProfileMode::Peak => idx += 1,
                _ => {}
            }
        }

        let (cache_rows, total_ops, total_tof) = self.profile_source(mode);
        let mut limit = 20usize;

        let rows: Vec<&SectionProfile> = match args.get(idx).copied() {
            None | Some("top") => {
                if matches!(args.get(idx), Some(&"top")) {
                    idx += 1;
                }
                if let Some(n) = args.get(idx).and_then(|s| s.parse::<usize>().ok()) {
                    limit = n;
                }
                cache_rows.iter().collect()
            }
            Some("exact") => {
                idx += 1;
                let Some(name) = args.get(idx) else {
                    return USAGE.into();
                };
                if let Some(n) = args.get(idx + 1).and_then(|s| s.parse::<usize>().ok()) {
                    limit = n;
                }
                cache_rows.iter().filter(|r| r.name == *name).collect()
            }
            Some("prefix") => {
                idx += 1;
                let Some(prefix) = args.get(idx) else {
                    return USAGE.into();
                };
                if let Some(n) = args.get(idx + 1).and_then(|s| s.parse::<usize>().ok()) {
                    limit = n;
                }
                cache_rows
                    .iter()
                    .filter(|r| r.name.starts_with(*prefix))
                    .collect()
            }
            Some("contains") => {
                idx += 1;
                let Some(substr) = args.get(idx) else {
                    return USAGE.into();
                };
                if let Some(n) = args.get(idx + 1).and_then(|s| s.parse::<usize>().ok()) {
                    limit = n;
                }
                cache_rows
                    .iter()
                    .filter(|r| r.name.contains(*substr))
                    .collect()
            }
            Some("current") => {
                idx += 1;
                if let Some(n) = args.get(idx).and_then(|s| s.parse::<usize>().ok()) {
                    limit = n;
                }
                let current = self.current_section().to_string();
                let child_prefix = format!("{current}/");
                cache_rows
                    .iter()
                    .filter(|r| r.name == current || r.name.starts_with(&child_prefix))
                    .collect()
            }
            Some("split") => {
                idx += 1;
                let Some(prefix_raw) = args.get(idx) else {
                    return USAGE.into();
                };
                if let Some(n) = args.get(idx + 1).and_then(|s| s.parse::<usize>().ok()) {
                    limit = n;
                }
                let prefix = if *prefix_raw == "current" {
                    self.current_section().to_string()
                } else {
                    (*prefix_raw).to_string()
                };
                let prefix = prefix.trim_end_matches('/').to_string();
                let child_prefix = format!("{prefix}/");
                let mut grouped = std::collections::BTreeMap::<String, (u64, u64)>::new();
                for row in &cache_rows {
                    if row.name == prefix {
                        let entry = grouped.entry(prefix.clone()).or_insert((0, 0));
                        entry.0 += row.ops;
                        entry.1 += row.tof;
                        continue;
                    }
                    if let Some(rest) = row.name.strip_prefix(&child_prefix) {
                        let child = rest.split('/').next().unwrap_or(rest);
                        let name = format!("{prefix}/{child}");
                        let entry = grouped.entry(name).or_insert((0, 0));
                        entry.0 += row.ops;
                        entry.1 += row.tof;
                    }
                }
                let rows: Vec<SectionProfile> = grouped
                    .into_iter()
                    .map(|(name, (ops, tof))| SectionProfile { name, ops, tof })
                    .collect();
                return Self::format_profile_rows(
                    rows.iter(),
                    total_ops,
                    total_tof,
                    limit,
                    sort_by_tof,
                    mode == ProfileMode::Peak,
                );
            }
            Some(_) => {
                return USAGE.into();
            }
        };

        Self::format_profile_rows(
            rows.iter().copied(),
            total_ops,
            total_tof,
            limit,
            sort_by_tof,
            mode == ProfileMode::Peak,
        )
    }

    pub fn profile(&mut self, args: &[&str]) -> String {
        self.print_profile(args)
    }

    fn run_command(&mut self, line: &str) -> ReplAction {
        let mut parts = line.split_whitespace();
        let cmd = parts.next().unwrap_or("");
        match cmd {
            "s" | "step" => {
                let n: usize = parts.next().and_then(|s| s.parse().ok()).unwrap_or(1);
                self.step(n);
                ReplAction::Print(self.where_line())
            }
            "b" | "back" => {
                let n: usize = parts.next().and_then(|s| s.parse().ok()).unwrap_or(1);
                self.back(n);
                ReplAction::Print(self.where_line())
            }
            "n" | "next" => {
                self.next_section();
                ReplAction::Print(self.where_line())
            }
            "c" | "continue" => match self.run_until_break() {
                Some(bp) => ReplAction::Print(format!("break: {:?}\n{}", bp, self.where_line())),
                None => ReplAction::Print(format!("(reached end)\n{}", self.where_line())),
            },
            "rc" | "rcontinue" | "rb" | "run-back" => match self.run_backward_until_break() {
                Some(bp) => {
                    ReplAction::Print(format!("rev-break: {:?}\n{}", bp, self.where_line()))
                }
                None => {
                    ReplAction::Print(format!("(hit op 0 without match)\n{}", self.where_line()))
                }
            },
            "g" | "goto" => {
                if let Some(tok) = parts.next() {
                    if let Some(name) = tok.strip_prefix('@') {
                        let found = self
                            .section_marks
                            .iter()
                            .find(|(_, n)| n == name)
                            .map(|(i, _)| *i);
                        match found {
                            Some(idx) => {
                                self.goto(idx);
                                ReplAction::Print(self.where_line())
                            }
                            None => ReplAction::Print(format!("no section named '{name}'")),
                        }
                    } else if let Ok(idx) = tok.parse::<usize>() {
                        self.goto(idx);
                        ReplAction::Print(self.where_line())
                    } else {
                        ReplAction::Print("usage: g <idx> | g @<section>".into())
                    }
                } else {
                    ReplAction::Print("usage: g <idx> | g @<section>".into())
                }
            }
            "p" | "print" | "show" => {
                let what = parts.next().unwrap_or("");
                ReplAction::Print(self.print_expr(what, &parts.collect::<Vec<_>>()))
            }
            "dump" => {
                let name = parts.next().unwrap_or("");
                let args: Vec<&str> = parts.collect();
                if name.is_empty() {
                    let names: Vec<String> = self.dump_hooks.keys().cloned().collect();
                    ReplAction::Print(format!("registered dump hooks: {names:?}"))
                } else if let Some(hook) = self.dump_hooks.get(name) {
                    ReplAction::Print(hook(self, &args))
                } else {
                    ReplAction::Print(format!("no such dump hook '{name}'"))
                }
            }
            "w" | "where" => ReplAction::Print(self.where_line()),
            "ghosts" | "gh" => ReplAction::Print(self.ghosts_line()),
            "l" | "list" => {
                let ctx: usize = parts.next().and_then(|s| s.parse().ok()).unwrap_or(10);
                ReplAction::Print(self.list(ctx))
            }
            "break" | "watch" => {
                let sub = parts.next().unwrap_or("");
                // Resolve "q<N>", bare "<N>", or a tag label to a
                // qubit id. Tags let callers use descriptive names.
                let parse_q_or_tag = |d: &Self, s: &str| -> Option<u32> {
                    if let Some(t) = s.strip_prefix('q') {
                        if let Ok(n) = t.parse() {
                            return Some(n);
                        }
                    }
                    if let Ok(n) = s.parse::<u32>() {
                        return Some(n);
                    }
                    d.qubit_id_from_tag(s)
                };
                let parse_mask = |s: &str| -> Option<u64> {
                    if let Some(hex) = s.strip_prefix("0x") {
                        u64::from_str_radix(hex, 16).ok()
                    } else if s == "all" {
                        Some(u64::MAX)
                    } else {
                        s.parse().ok()
                    }
                };
                let bp = match sub {
                    "op" => parts
                        .next()
                        .and_then(|s| s.parse().ok())
                        .map(Breakpoint::Op),
                    "section" => parts
                        .next()
                        .map(|s| Breakpoint::SectionStart(s.to_string())),
                    "phase" => Some(Breakpoint::PhaseNonzero),
                    "r-nonzero" => Some(Breakpoint::ROnNonzero),
                    // `watch <label-or-q<N>> [= <mask>]`
                    s if parse_q_or_tag(self, s).is_some() => {
                        let q = parse_q_or_tag(self, s).unwrap();
                        let rest: Vec<_> = parts.collect();
                        if rest.len() >= 2 && rest[0] == "=" {
                            parse_mask(rest[1]).map(|m| Breakpoint::QubitValue { q, expected: m })
                        } else {
                            let last = self.qubit(q);
                            Some(Breakpoint::QubitChange { q, last })
                        }
                    }
                    _ => None,
                };
                match bp {
                    Some(bp) => {
                        self.breakpoints.push(bp.clone());
                        ReplAction::Print(format!("added: {bp:?}"))
                    }
                    None => ReplAction::Print(
                        "usage: break {op <idx> | section <n> | phase \
                         | r-nonzero | <label-or-q<N>> [= <mask>]}"
                            .into(),
                    ),
                }
            }
            "breakpoints" | "bp" => {
                let mut s = String::new();
                for (i, bp) in self.breakpoints.iter().enumerate() {
                    let _ = writeln!(s, "  {i}: {bp:?}");
                }
                if s.is_empty() {
                    s = "(none)".into();
                }
                ReplAction::Print(s)
            }
            "clear" => {
                if let Some(idx) = parts.next().and_then(|s| s.parse().ok()) {
                    self.clear_breakpoint(idx);
                    ReplAction::Print("ok".into())
                } else {
                    ReplAction::Print("usage: clear <idx>".into())
                }
            }
            "src" | "source" => {
                let idx: usize = parts
                    .next()
                    .and_then(|s| s.parse().ok())
                    .unwrap_or_else(|| self.cursor());
                ReplAction::Print(self.source_trace(idx))
            }
            "prof" | "profile" => {
                let rest: Vec<_> = parts.collect();
                ReplAction::Print(self.print_profile(&rest))
            }
            "tag" => {
                // Accept q<N>, <N>, or a tag label. For labels,
                // resolve to the qubit id via qubit_id_from_tag.
                let arg = parts.next().unwrap_or("");
                let qid: Option<u32> = arg
                    .strip_prefix('q')
                    .and_then(|s| s.parse().ok())
                    .or_else(|| arg.parse().ok())
                    .or_else(|| self.qubit_id_from_tag(arg));
                let gc = self.cursor();
                match qid {
                    Some(q) => ReplAction::Print(format!(
                        "q{} @ op {} -> {}",
                        q,
                        gc,
                        self.qubit_name_at(q, gc)
                    )),
                    None => ReplAction::Print("usage: tag q<N> | tag <N> | tag <label>".into()),
                }
            }
            "q" | "quit" | "exit" => ReplAction::Quit,
            "h" | "help" | "?" => ReplAction::Print(HELP.into()),
            _ => ReplAction::Print(format!("unknown command '{cmd}'; type `h` for help")),
        }
    }

    fn print_expr(&self, what: &str, rest: &[&str]) -> String {
        if what == "phase" {
            return format!(
                "phase = {:#018x} (shots with bit=1: {:?})",
                self.state.phase,
                (0..64)
                    .filter(|i| (self.state.phase >> i) & 1 != 0)
                    .collect::<Vec<_>>()
            );
        }
        if what == "section" {
            return self.current_section().to_string();
        }
        if what == "reg" {
            let start: u32 = rest.first().and_then(|s| s.parse().ok()).unwrap_or(0);
            let n: usize = rest.get(1).and_then(|s| s.parse().ok()).unwrap_or(256);
            let shot: usize = rest
                .get(2)
                .and_then(|s| {
                    s.strip_prefix("shot")
                        .and_then(|x| x.parse().ok())
                        .or_else(|| s.parse().ok())
                })
                .unwrap_or(0);
            let qs: Vec<u32> = (start..start + n as u32).collect();
            let b = self.read_reg_shot(&qs, shot);
            let mut hex = String::from("0x");
            for &x in b.iter().rev() {
                let _ = write!(hex, "{x:02x}");
            }
            return format!(
                "reg q{}..q{} shot{} = {}",
                start,
                start + n as u32 - 1,
                shot,
                hex
            );
        }
        if what == "qreg" {
            // p qreg <label> — read scattered qubits of a multi-bit
            // register allocated as `<label>[0]`, `<label>[1]`, etc.
            // Resolves each index via alloc_log, then assembles a U256.
            // `<label>:be` reverses bit order (label[0] = MSB).
            let label_raw = match rest.first() {
                Some(l) => *l,
                None => return "usage: p qreg <label>[:be]".into(),
            };
            let (label, be_order): (&str, bool) = match label_raw.strip_suffix(":be") {
                Some(s) => (s, true),
                None => (label_raw, false),
            };
            // Match alloc events whose tag ends in "<label>[N]". Tags
            // may include a section path prefix (e.g. "ec/.../label[0]"),
            // so suffix-match the label part.
            let mut by_idx: std::collections::HashMap<usize, (u32, usize)> =
                std::collections::HashMap::new();
            for ev in &self.qubit_alloc_log {
                if ev.op_idx > self.cursor() {
                    continue;
                }
                let tag: &str = &ev.tag;
                let Some(bracket_close) = tag.strip_suffix(']') else {
                    continue;
                };
                let Some(bracket_open_pos) = bracket_close.rfind('[') else {
                    continue;
                };
                let idx_str = &bracket_close[bracket_open_pos + 1..];
                let label_part = &bracket_close[..bracket_open_pos];
                let matches = label_part == label || label_part.ends_with(&format!("/{label}"));
                if !matches {
                    continue;
                }
                if let Ok(idx) = idx_str.parse::<usize>() {
                    by_idx
                        .entry(idx)
                        .and_modify(|(q, op)| {
                            if ev.op_idx >= *op {
                                *q = ev.qubit;
                                *op = ev.op_idx;
                            }
                        })
                        .or_insert((ev.qubit, ev.op_idx));
                }
            }
            if by_idx.is_empty() {
                return format!("no alloc events for label '{label}'");
            }
            let max_idx = *by_idx.keys().max().unwrap();
            let qs: Vec<u32> = (0..=max_idx)
                .map(|i| by_idx.get(&i).map_or(u32::MAX, |(q, _)| *q))
                .collect();
            // Read each qubit's bit for the requested shot. With be_order,
            // the qubit at index `i` is treated as bit (max_idx - i).
            // Optional rest[1] = "shotN" or plain integer selects shot.
            let shot: usize = rest
                .get(1)
                .and_then(|s| {
                    s.strip_prefix("shot")
                        .and_then(|x| x.parse().ok())
                        .or_else(|| s.parse().ok())
                })
                .unwrap_or(0);
            let n_bits = max_idx + 1;
            let mut bits: Vec<u8> = vec![0; n_bits.div_ceil(8)];
            for (i, &q) in qs.iter().enumerate() {
                if q == u32::MAX {
                    continue;
                }
                let bit = (self.qubit(q) >> shot) & 1;
                let pos = if be_order { max_idx - i } else { i };
                if bit != 0 {
                    bits[pos / 8] |= 1 << (pos % 8);
                }
            }
            let mut hex = String::from("0x");
            for &x in bits.iter().rev() {
                let _ = write!(hex, "{x:02x}");
            }
            let order_tag = if be_order { " (BE: label[0]=MSB)" } else { "" };
            return format!("{label}[0..{n_bits}] shot{shot} = {hex}{order_tag}  (qubits: {qs:?})");
        }
        if what == "creg" {
            // p creg <label> — read scattered cbits of a multi-bit
            // classical register allocated as `<label>[0]`, `<label>[1]`,
            // etc., mirroring `p qreg`. `<label>:be` reverses bit order.
            let label_raw = match rest.first() {
                Some(l) => *l,
                None => return "usage: p creg <label>[:be]".into(),
            };
            let (label, be_order): (&str, bool) = match label_raw.strip_suffix(":be") {
                Some(s) => (s, true),
                None => (label_raw, false),
            };
            let mut by_idx: std::collections::HashMap<usize, (u32, usize)> =
                std::collections::HashMap::new();
            for ev in &self.bit_alloc_log {
                if ev.op_idx > self.cursor() {
                    continue;
                }
                let tag: &str = &ev.tag;
                let Some(bracket_close) = tag.strip_suffix(']') else {
                    continue;
                };
                let Some(bracket_open_pos) = bracket_close.rfind('[') else {
                    continue;
                };
                let idx_str = &bracket_close[bracket_open_pos + 1..];
                let label_part = &bracket_close[..bracket_open_pos];
                let matches = label_part == label || label_part.ends_with(&format!("/{label}"));
                if !matches {
                    continue;
                }
                if let Ok(idx) = idx_str.parse::<usize>() {
                    by_idx
                        .entry(idx)
                        .and_modify(|(b, op)| {
                            if ev.op_idx >= *op {
                                *b = ev.bit;
                                *op = ev.op_idx;
                            }
                        })
                        .or_insert((ev.bit, ev.op_idx));
                }
            }
            if by_idx.is_empty() {
                return format!("no cbit alloc events for label '{label}'");
            }
            let max_idx = *by_idx.keys().max().unwrap();
            let bs: Vec<u32> = (0..=max_idx)
                .map(|i| by_idx.get(&i).map_or(u32::MAX, |(b, _)| *b))
                .collect();
            let n_bits = max_idx + 1;
            let mut bits: Vec<u8> = vec![0; n_bits.div_ceil(8)];
            for (i, &b) in bs.iter().enumerate() {
                if b == u32::MAX {
                    continue;
                }
                let bit = self.bit(b) & 1; // shot 0
                let pos = if be_order { max_idx - i } else { i };
                if bit != 0 {
                    bits[pos / 8] |= 1 << (pos % 8);
                }
            }
            let mut hex = String::from("0x");
            for &x in bits.iter().rev() {
                let _ = write!(hex, "{x:02x}");
            }
            let order_tag = if be_order { " (BE: label[0]=MSB)" } else { "" };
            return format!("{label}[0..{n_bits}] shot0 = {hex}{order_tag}  (cbits: {bs:?})");
        }
        // qN, qN.sS, bN, bN.sS, or a tag label (with optional .sS).
        let (head, shot_str) = match what.find(".s") {
            Some(i) => (&what[..i], Some(&what[i + 2..])),
            None => (what, None),
        };
        let qid: Option<u32> = head
            .strip_prefix('q')
            .and_then(|s| s.parse().ok())
            .or_else(|| self.qubit_id_from_tag(head));
        if let Some(id) = qid {
            let mask = self.qubit(id);
            let name = self.qubit_name_at(id, self.cursor());
            if let Some(shot) = shot_str.and_then(|s| s.parse::<usize>().ok()) {
                return format!("{} (q{}).s{} = {}", name, id, shot, (mask >> shot) & 1);
            }
            return format!(
                "{} (q{}) = {:#018x} ({} shots=1)",
                name,
                id,
                mask,
                mask.count_ones()
            );
        }
        if let Some(stripped) = what.strip_prefix('b') {
            let (id_str, shot_str) = match stripped.find(".s") {
                Some(i) => (&stripped[..i], Some(&stripped[i + 2..])),
                None => (stripped, None),
            };
            if let Ok(id) = id_str.parse::<u32>() {
                let mask = self.bit(id);
                if let Some(shot) = shot_str.and_then(|s| s.parse::<usize>().ok()) {
                    return format!("b{}.s{} = {}", id, shot, (mask >> shot) & 1);
                }
                return format!("b{id} = {mask:#018x}");
            }
        }
        // Fallback: try as a cbit tag (named alloc_bit_named label).
        // Lets `p <cbit-tag>` resolve without an explicit `creg` prefix,
        // mirroring `p <qubit-tag>`.
        if let Some(id) = self.cbit_id_from_tag(head) {
            let mask = self.bit(id);
            let name = self.cbit_name_at(id, self.cursor());
            if let Some(shot) = shot_str.and_then(|s| s.parse::<usize>().ok()) {
                return format!("{} (b{}).s{} = {}", name, id, shot, (mask >> shot) & 1);
            }
            return format!("{name} (b{id}) = {mask:#018x}");
        }
        format!("don't know how to print '{what}'")
    }
}

enum ReplAction {
    Quit,
    Print(String),
}

#[cfg(test)]
mod tests;

const HELP: &str = "Commands:\n\
    s / step [n]       advance n ops (default 1)\n\
    b / back [n]       rewind n ops\n\
    n / next           advance to next section\n\
    c / continue       run until breakpoint or end\n\
    g / goto <idx>     jump to op index\n\
    g @<section>       jump to start of named section\n\
    p / show q<N>[.s<S>]      print qubit mask or specific shot bit\n\
    p / show b<N>[.s<S>]      print classical bit mask\n\
    p / show phase            print sim_phase mask + failing shots\n\
    p / show reg <start> <n>  print n-qubit register as U256 hex (shot 0)\n\
    p / show qreg <label>     print scattered multi-bit register by label (e.g. 'p qreg ct_scratch_pa')\n\
    p / show creg <label>     print scattered multi-bit classical register by label\n\
    p / show section          current section name\n\
    w / where          cursor, section, phase, r-nonzero count\n\
    ghosts / gh        pending spooky-pebble ghosts at cursor (id, anchor, bit, tape mask, create site)\n\
    l / list [n]       window of n ops around cursor\n\
    break op <idx>     break at op index\n\
    break section <n>  break entering section\n\
    break phase        break when sim_phase first goes non-zero\n\
    break r-nonzero    break on first R-on-non-zero\n\
    breakpoints / bp   list breakpoints\n\
    clear <idx>        remove breakpoint by list index\n\
    src / source [idx] print scope chain (file:line) for op idx\n\
    prof / profile     section / peak-qubit profiler\n\
    tag <q-or-label>   print tag for qubit (accepts q<N>, <N>, or label)\n\
    rb / run-back      walk backward until a breakpoint fires\n\
    break / watch <spec>  watch qubit changes / values (see below)\n\
    q / quit           exit debugger\n\
    \n\
  break/watch specs: op <idx> | section <name> | phase | r-nonzero\n\
    | q<N> | <label> | <q-or-label> = <mask>  (mask: 0x<hex>|<dec>|all)\n\
  profiler:\n\
    prof [whole|cursor|peak] [ops|tof|q] top [n]                  biggest exact sections / peak tags\n\
    prof [whole|cursor|peak] [ops|tof|q] exact <name> [n]         one exact section / peak tag\n\
    prof [whole|cursor|peak] [ops|tof|q] prefix <prefix> [n]      all exact sections / peak tags under prefix\n\
    prof [whole|cursor|peak] [ops|tof|q] contains <substr> [n]    substring match\n\
    prof [whole|cursor|peak] [ops|tof|q] current [n]              current section and its children (whole/cursor only)\n\
    prof [whole|cursor|peak] [ops|tof|q] split <prefix|current> [n] immediate children rolled up";

// ── Sim transfer functions (mirror circuit.rs's gate methods) ──

fn cond_mask(state: &SimState) -> u64 {
    let mut m = u64::MAX;
    for &b in &state.cond_stack {
        m &= *state.bits.get(b as usize).unwrap_or(&0);
    }
    m
}

fn ensure_qubit(state: &mut SimState, q: u32) {
    let qi = q as usize;
    if qi >= state.qubits.len() {
        state.qubits.resize(qi + 256, 0);
    }
}
fn ensure_bit(state: &mut SimState, b: u32) {
    let bi = b as usize;
    if bi >= state.bits.len() {
        state.bits.resize(bi + 256, 0);
    }
}

fn rng_u64(counter: u64, seed: u64) -> u64 {
    let c = counter;
    let mut z = c.wrapping_mul(0x9E37_79B9_7F4A_7C15).wrapping_add(seed);
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

/// Per-op undo data. Most ops are self-inverse (re-applying them
/// reverses the state change and doesn't touch `hmr_counter` or
/// `r_on_nonzero_events`); for those we store `SelfInverse`. The
/// rest record minimal pre-op state.
#[derive(Clone)]
enum Delta {
    SelfInverse,
    Hmr { old_q: u64, old_b: u64 },
    R { old_q: u64, was_nonzero: bool },
    PushCond,
    PopCond { restored_bit: Option<u32> },
    BitStore { old: u64 },
    None,
}

/// Apply an op forward while capturing undo info.
fn apply_op_recorded(op: &Op, state: &mut SimState, hmr_seed: u64) -> Delta {
    match *op {
        Op::Hmr(q, b) => {
            ensure_qubit(state, q);
            ensure_bit(state, b);
            let old_q = state.qubits[q as usize];
            let old_b = state.bits[b as usize];
            state.hmr_counter = state.hmr_counter.wrapping_add(1);
            let rng = rng_u64(state.hmr_counter, hmr_seed);
            state.bits[b as usize] = rng;
            state.phase ^= old_q & rng;
            state.qubits[q as usize] = 0;
            Delta::Hmr { old_q, old_b }
        }
        Op::R(q) => {
            ensure_qubit(state, q);
            let cm = cond_mask(state);
            let old_q = state.qubits[q as usize];
            state.hmr_counter = state.hmr_counter.wrapping_add(1);
            let rng = rng_u64(state.hmr_counter, hmr_seed);
            state.phase ^= old_q & rng & cm;
            let was_nonzero = old_q != 0;
            if was_nonzero {
                state.r_on_nonzero_events += 1;
            }
            state.qubits[q as usize] = 0;
            Delta::R { old_q, was_nonzero }
        }
        Op::PushCondition(b) => {
            ensure_bit(state, b);
            state.cond_stack.push(b);
            Delta::PushCond
        }
        Op::PopCondition => {
            let restored_bit = state.cond_stack.pop();
            Delta::PopCond { restored_bit }
        }
        Op::BitStore0(b) => {
            ensure_bit(state, b);
            let old = state.bits[b as usize];
            let cm = cond_mask(state);
            state.bits[b as usize] &= !cm;
            Delta::BitStore { old }
        }
        Op::BitStore1(b) => {
            ensure_bit(state, b);
            let old = state.bits[b as usize];
            let cm = cond_mask(state);
            state.bits[b as usize] |= cm;
            Delta::BitStore { old }
        }
        // BIT_INVERT is self-inverse only if b is not in the
        // cond_stack (directly or transitively). If b is conditioning
        // itself, applying twice yields a different cm on the second
        // pass. Always record the old value for safety; it is also
        // only 8 bytes per entry.
        Op::BitInvert(b) => {
            ensure_bit(state, b);
            let old = state.bits[b as usize];
            let cm = cond_mask(state);
            state.bits[b as usize] ^= cm;
            Delta::BitStore { old }
        }
        // Register/Append ops have no sim effect.
        Op::Register(_) | Op::AppendQubit(_, _) | Op::AppendBit(_, _) => Delta::None,
        // All remaining variants are self-inverse: applying the op
        // a second time reverts its state effect and never touches
        // the hmr counter or r-on-nonzero counter.
        _ => {
            apply_op(op, state, hmr_seed);
            Delta::SelfInverse
        }
    }
}

/// Undo an op using its recorded delta. For self-inverse ops we
/// re-apply the op; for destructive ops we restore captured state.
fn reverse_op(op: &Op, delta: &Delta, state: &mut SimState, hmr_seed: u64) {
    match (op, delta) {
        (_, Delta::SelfInverse) => apply_op(op, state, hmr_seed),
        (_, Delta::None) => {}
        (Op::Hmr(q, b), Delta::Hmr { old_q, old_b }) => {
            let rng = rng_u64(state.hmr_counter, hmr_seed);
            state.phase ^= *old_q & rng;
            state.bits[*b as usize] = *old_b;
            state.qubits[*q as usize] = *old_q;
            state.hmr_counter = state.hmr_counter.wrapping_sub(1);
        }
        (Op::R(q), Delta::R { old_q, was_nonzero }) => {
            let rng = rng_u64(state.hmr_counter, hmr_seed);
            // cond_mask at forward time equals cond_mask now (the
            // condition stack is unaffected between forward and
            // reverse of R).
            let cm = cond_mask(state);
            state.phase ^= *old_q & rng & cm;
            state.qubits[*q as usize] = *old_q;
            if *was_nonzero {
                state.r_on_nonzero_events -= 1;
            }
            state.hmr_counter = state.hmr_counter.wrapping_sub(1);
        }
        (Op::PushCondition(_), Delta::PushCond) => {
            state.cond_stack.pop();
        }
        (Op::PopCondition, Delta::PopCond { restored_bit }) => {
            if let Some(b) = restored_bit {
                state.cond_stack.push(*b);
            }
        }
        (Op::BitStore0(b) | Op::BitStore1(b) | Op::BitInvert(b), Delta::BitStore { old }) => {
            state.bits[*b as usize] = *old;
        }
        _ => panic!("reverse_op: delta/op mismatch {op:?}"),
    }
}

fn apply_op(op: &Op, state: &mut SimState, hmr_seed: u64) {
    match *op {
        Op::X(q) => {
            ensure_qubit(state, q);
            let c = cond_mask(state);
            state.qubits[q as usize] ^= c;
        }
        Op::Z(q) => {
            ensure_qubit(state, q);
            let c = cond_mask(state);
            state.phase ^= c & state.qubits[q as usize];
        }
        Op::Cx(c1, t) => {
            ensure_qubit(state, c1);
            ensure_qubit(state, t);
            let cm = cond_mask(state);
            let v = cm & state.qubits[c1 as usize];
            state.qubits[t as usize] ^= v;
        }
        Op::Cz(a, b) => {
            ensure_qubit(state, a);
            ensure_qubit(state, b);
            let cm = cond_mask(state);
            state.phase ^= cm & state.qubits[a as usize] & state.qubits[b as usize];
        }
        Op::Ccx(c1, c2, t) => {
            ensure_qubit(state, c1);
            ensure_qubit(state, c2);
            ensure_qubit(state, t);
            let cm = cond_mask(state);
            let v = cm & state.qubits[c1 as usize] & state.qubits[c2 as usize];
            state.qubits[t as usize] ^= v;
        }
        Op::Ccz(a, b, c) => {
            ensure_qubit(state, a);
            ensure_qubit(state, b);
            ensure_qubit(state, c);
            let cm = cond_mask(state);
            state.phase ^=
                cm & state.qubits[a as usize] & state.qubits[b as usize] & state.qubits[c as usize];
        }
        Op::Swap(a, b) => {
            ensure_qubit(state, a);
            ensure_qubit(state, b);
            let cm = cond_mask(state);
            let mut qa = state.qubits[a as usize];
            let mut qb = state.qubits[b as usize];
            qa ^= qb;
            qb ^= cm & qa;
            qa ^= qb;
            state.qubits[a as usize] = qa;
            state.qubits[b as usize] = qb;
        }
        Op::Hmr(q, b) => {
            ensure_qubit(state, q);
            ensure_bit(state, b);
            state.hmr_counter = state.hmr_counter.wrapping_add(1);
            let rng = rng_u64(state.hmr_counter, hmr_seed);
            let qval = state.qubits[q as usize];
            state.bits[b as usize] = rng;
            state.phase ^= qval & rng;
            state.qubits[q as usize] = 0;
        }
        Op::R(q) => {
            ensure_qubit(state, q);
            let cm = cond_mask(state);
            state.hmr_counter = state.hmr_counter.wrapping_add(1);
            let rng = rng_u64(state.hmr_counter, hmr_seed);
            let qval = state.qubits[q as usize];
            state.phase ^= qval & rng & cm;
            if qval != 0 {
                state.r_on_nonzero_events += 1;
            }
            state.qubits[q as usize] = 0;
        }
        Op::Neg => {
            let cm = cond_mask(state);
            state.phase ^= cm;
        }
        Op::PushCondition(b) => {
            ensure_bit(state, b);
            state.cond_stack.push(b);
        }
        Op::PopCondition => {
            state.cond_stack.pop();
        }
        Op::BitInvert(b) => {
            ensure_bit(state, b);
            let cm = cond_mask(state);
            state.bits[b as usize] ^= cm;
        }
        Op::BitStore0(b) => {
            ensure_bit(state, b);
            let cm = cond_mask(state);
            state.bits[b as usize] &= !cm;
        }
        Op::BitStore1(b) => {
            ensure_bit(state, b);
            let cm = cond_mask(state);
            state.bits[b as usize] |= cm;
        }
        // Register/Append ops have no sim effect.
        Op::Register(_) | Op::AppendQubit(_, _) | Op::AppendBit(_, _) => {}
    }
}
