//! Spooky-pebble / deferred-discharge HMR API.
//!
//! A `Ghost` captures the 64-shot mask and anchor of an HMR'd qubit; the
//! `ghost_xor_*` calls accumulate the deferred phase obligation, and
//! `resolve_ghost` / `close_ghost` discharge it (sim-verified).
//!
//! Quick mental model:
//!
//!   forward pass:  q := `f(alive_regs)`
//!                  let g = `circ.hmr_ghost(&q)`;  // q is now |0>
//!                  `circ.zero_and_free(q)`;
//!                  // ... many ops, possibly mutating `alive_regs` ...
//!
//!   reverse pass:  let r = `circ.alloc_qreg`(...);
//!                  // recompute r := f(restored `alive_regs`)
//!                  `circ.resolve_ghost(g`, &r);  // discharges phase
//!                  // uncompute r; `zero_and_free`
//!
//! The Ghost type is `!Clone` / `!Copy` / `#[must_use]`. Dropping
//! without resolve is a panic (phase leak).

use std::panic::Location;

use crate::point_add::trailmix_port::circuit::{Cbit, Circuit, QReg};

// ── The opaque receipt ─────────────────────────────────────────

/// Opaque handle to a deferred HMR phase obligation (the "spooky
/// pebble"). Created by [`Circuit::hmr_ghost`]; resolved exactly
/// once via [`Circuit::resolve_ghost`].
///
/// Dropping a Ghost without resolve is a **phase leak** and panics
/// (with a debugger hook under `DEBUG_ON_FAIL=1`).
///
/// `!Clone`, `!Copy`, `#[must_use]`. The internal fields are
/// crate-private; outside callers can only thread the value
/// through the resolve API.
#[must_use = "Ghost must be resolved via Circuit::resolve_ghost; \
              dropping it is a phase leak"]
#[derive(Debug)]
pub struct Ghost {
    /// Sequential id for diagnostics — also used by
    /// `Circuit::pending_ghosts` as the bookkeeping key.
    pub(crate) id: u64,

    /// The random classical bit allocated at HMR time. Stays
    /// alive until `resolve_ghost` discharges + frees it.
    pub(crate) bit: Cbit,

    /// Anchor id installed in the tracker at HMR time. The
    /// obligation's `AbsVal` is `Anchor(anchor_id)`; the resolve
    /// target is promoted to the same `Anchor(anchor_id)` before
    /// `z_if_bit` fires, so structural cancellation succeeds.
    pub(crate) anchor_id: u64,

    /// 64-shot bitmask of `q`'s value at HMR time. Bit i = shot
    /// i's value. Used by `resolve_ghost` to sim-verify the
    /// discharge target matches.
    pub(crate) mask_at_hmr: u64,

    /// Diagnostics: section path at HMR time.
    pub(crate) hmr_section: String,

    /// Diagnostics: total op index at HMR time
    /// (`ops_truncated + ops.len()`).
    pub(crate) hmr_op_idx: usize,

    /// Diagnostics: source location of the `hmr_ghost` call, via
    /// `#[track_caller]` + `Location::caller()`.
    pub(crate) hmr_caller: &'static Location<'static>,

    /// Set by `resolve_ghost` immediately before forgetting the
    /// Ghost via consumption. Drop checks this; if `false` and
    /// not already unwinding, Drop panics.
    pub(crate) consumed: bool,

    /// Running XOR of the 64-shot sim masks of the discharge terms
    /// deposited via `ghost_xor_z` / `ghost_xor_cz` / `ghost_xor_ccz`.
    /// `close_ghost` asserts this equals `mask_at_hmr` — i.e. the
    /// multi-term discharge reproduces the vented value on every shot,
    /// so the deposited Z/CZ/CCZ phases cancel the HMR kickback for any
    /// random bit value. Stays 0 for the single-term `resolve_ghost`
    /// path (which uses anchor matching instead).
    pub(crate) acc_xor: u64,
}

impl Ghost {
    /// Internal constructor. Not exposed; outside callers go
    /// through `Circuit::hmr_ghost`.
    pub(crate) fn new(
        id: u64,
        bit: Cbit,
        anchor_id: u64,
        mask_at_hmr: u64,
        hmr_section: String,
        hmr_op_idx: usize,
        hmr_caller: &'static Location<'static>,
    ) -> Self {
        Self {
            id,
            bit,
            anchor_id,
            mask_at_hmr,
            hmr_section,
            hmr_op_idx,
            hmr_caller,
            consumed: false,
            acc_xor: 0,
        }
    }
}

impl Drop for Ghost {
    fn drop(&mut self) {
        if self.consumed || std::thread::panicking() {
            return;
        }
        // DEBUG_ON_FAIL: print the HMR origin so the user can goto
        // it in a debugger session. We can't attach the debugger
        // here because we don't hold &mut Circuit; printing the
        // op idx is the next best thing.
        if std::env::var("DEBUG_ON_FAIL").is_ok() {
            eprintln!(
                "[ghost] Ghost #{} dropped without resolve.\n\
                 [ghost]   HMR at {} (op #{}, section '{}')\n\
                 [ghost]   To inspect: rerun with DEBUG_ON_FAIL=1 \
                 against a panic earlier in the build, or attach \
                 the debugger and `goto {}`.",
                self.id, self.hmr_caller, self.hmr_op_idx, self.hmr_section, self.hmr_op_idx,
            );
        }
        panic!(
            "Ghost #{} from {} (op #{}, section '{}') dropped \
             without resolve_ghost — phase leak.",
            self.id, self.hmr_caller, self.hmr_op_idx, self.hmr_section,
        );
    }
}

// ── Per-Circuit pending-ghost bookkeeping ──────────────────────

/// One row of `Circuit::pending_ghosts`. Mirrors the diagnostic
/// fields of `Ghost` so the Circuit-level Drop check can produce
/// a useful error message even if the Ghost itself was leaked
/// via `std::mem::forget`.
#[derive(Debug, Clone)]
pub struct GhostRecord {
    pub id: u64,
    pub anchor_id: u64,
    pub bit_raw: u32,
    /// 64-shot bitmask of the HMR'd qubit's value at ghost-creation
    /// time (bit i = shot i). This IS the tape bit the ghost stands
    /// for; contracts use it to assert the pending-ghost set matches a
    /// reference decision tape, and that a discharge witness reproduces
    /// it. Mirrors `Ghost::mask_at_hmr`.
    pub mask_at_hmr: u64,
    pub hmr_section: String,
    pub hmr_op_idx: usize,
    pub hmr_caller: &'static Location<'static>,
}

impl From<&Ghost> for GhostRecord {
    fn from(g: &Ghost) -> Self {
        GhostRecord {
            id: g.id,
            anchor_id: g.anchor_id,
            bit_raw: g.bit.raw(),
            mask_at_hmr: g.mask_at_hmr,
            hmr_section: g.hmr_section.clone(),
            hmr_op_idx: g.hmr_op_idx,
            hmr_caller: g.hmr_caller,
        }
    }
}

// Circuit-side ghost API: `hmr_ghost` measures a qubit into a deferred-discharge
// receipt (promoting its tracker value to an anchor); `resolve_ghost` sim-verifies
// the discharge register and emits the cancelling phase fix-up.

impl Circuit {
    /// HMR `q` and return a [`Ghost`] receipt for deferred
    /// discharge.
    ///
    /// Equivalent to `hmr(&q, fresh_bit)` plus a tracker anchor
    /// promotion that detaches the HMR obligation from any specific
    /// `(qubit, version)` pair. The returned `Ghost` carries the
    /// Cbit, the 64-shot sim mask of `q` at HMR time, and the anchor
    /// id; pass it to [`Circuit::resolve_ghost`] together with a
    /// register whose per-shot value matches.
    ///
    /// Forbidden inside `push_condition` (same as `hmr`).
    #[track_caller]
    pub fn hmr_ghost(&mut self, q: &QReg) -> Ghost {
        let caller = Location::caller();
        // 1. HMR is forbidden under push_condition (same rule as
        //    `hmr()`). We capture an EMPTY
        //    condition mask here; the strict-semantic check at
        //    resolve_ghost re-asserts emptiness.
        assert!(
            self.sim_condition_stack.is_empty(),
            "hmr_ghost q{} cannot be inside a push_condition block \
             — HMR under a condition is forbidden (same rule as \
             `hmr`). Called at {}.",
            q.id(),
            caller,
        );
        // 2. Snapshot the per-shot value of q BEFORE on_hmr zeros
        //    it. resolve_ghost will verify the discharge target
        //    matches this mask. With sim disabled, sim_get_mask
        //    returns 0; the resolve-time check is then trivially
        //    skipped because r_mask will also be 0 (both sim-off).
        let mask = self.sim_get_mask(q.id());
        // 3. Anchor-promote q in the tracker BEFORE emitting the
        //    HMR. The HMR's `on_hmr` reads val(q), so the
        //    obligation gets `val = Anchor(anchor_id)` (which is
        //    always version-current).
        let anchor_id = self.phase.anchor_qubit_and_get_id(q.id());
        // 4. Allocate the random classical bit.
        let bit = self.alloc_bit();
        // 5. Standard HMR: emits Op::Hmr, zeros q's sim, records
        //    the obligation. We don't follow with a discharge here
        //    — that's the whole point of Ghost.
        self.hmr(q, bit);
        // 5b. Record the vented value's 64-shot mask on the obligation
        //     so a multi-term (`ghost_xor_*` + `close_ghost`) discharge
        //     can be sim-verified by the tracker.
        self.phase.set_obligation_hmr_mask(bit.raw(), mask);
        // 6. Build the receipt + bookkeeping row.
        let id = self.next_ghost_id;
        self.next_ghost_id += 1;
        let hmr_op_idx = self.ops_truncated as usize + self.ops.len();
        let section = self.current_section.clone();
        // Debugger ghost-event log (create).
        self.ghost_event_log
            .push(crate::point_add::trailmix_port::circuit::GhostEvent {
                op_idx: hmr_op_idx,
                create: true,
                id,
                anchor_id,
                mask_at_hmr: mask,
                bit_raw: bit.raw(),
                section: section.clone(),
            });
        let g = Ghost::new(id, bit, anchor_id, mask, section, hmr_op_idx, caller);
        self.pending_ghosts.push(GhostRecord::from(&g));
        g
    }

    /// Discharge a [`Ghost`] by providing a register `r` whose
    /// per-shot value matches `q`'s value at the HMR call site.
    /// Emits `z_if_bit(r, g.bit)`, frees `g.bit`, and removes the
    /// ghost from the pending-ghost table.
    ///
    /// Forbidden inside `push_condition` (strict semantic). The HMR side also
    /// requires an empty stack, so both ends of a Ghost agree on
    /// the trivial-condition case.
    ///
    /// Panics (with `DEBUG_ON_FAIL=1` attach) on:
    ///   - mismatched 64-shot sim mask (`r != q` on some shot)
    ///   - non-empty condition stack
    #[track_caller]
    pub fn resolve_ghost(&mut self, mut g: Ghost, r: &QReg) {
        let caller = Location::caller();
        // 1. Strict semantic: empty condition stack required at
        //    resolve time too.
        assert!(
            self.sim_condition_stack.is_empty(),
            "resolve_ghost: condition stack must be empty at resolve \
             time (depth = {}). Ghost #{} HMR at {} (op #{}, \
             section '{}'); resolve attempted at {}.",
            self.sim_condition_stack.len(),
            g.id,
            g.hmr_caller,
            g.hmr_op_idx,
            g.hmr_section,
            caller,
        );
        // 2. Sim-verify that r matches the captured q-mask. Mirror
        //    declare_identity_raw's DEBUG_ON_FAIL hook so the
        //    debugger drops in at the failing op for inspection.
        let r_mask = self.sim_get_mask(r.id());
        if r_mask != g.mask_at_hmr {
            let last_op = self.ops_truncated as usize + self.ops.len();
            let section = self.current_section.clone();
            let diff = r_mask ^ g.mask_at_hmr;
            if std::env::var("DEBUG_ON_FAIL").is_ok() {
                eprintln!(
                    "[resolve_ghost] r_mask mismatch: \
                     r_mask=q{}={:#x} vs g.mask_at_hmr={:#x} (diff={:#x}) \
                     Ghost #{} (anchor_id={}, bit=b{}) \
                     HMR at {} (op #{}, section '{}') \
                     resolve at {} (op #{}, section '{}')",
                    r.id(),
                    r_mask,
                    g.mask_at_hmr,
                    diff,
                    g.id,
                    g.anchor_id,
                    g.bit.raw(),
                    g.hmr_caller,
                    g.hmr_op_idx,
                    g.hmr_section,
                    caller,
                    last_op,
                    section,
                );
                eprintln!("[debugger] DEBUG_ON_FAIL set — attaching");
                let mut dbg =
                    crate::point_add::trailmix_port::tracker::debugger::Debugger::attach(self);
                dbg.goto(last_op);
                dbg.repl();
            }
            // Mark consumed BEFORE we panic so the Ghost's own Drop
            // doesn't double-panic during unwinding (the
            // panicking() guard already handles unwind, but being
            // explicit is harmless and matches the success-path
            // mark below).
            g.consumed = true;
            panic!(
                "resolve_ghost: r_mask mismatch (Ghost #{} from {} \
                 op #{}): r=q{} mask {:#x} != captured {:#x} (diff \
                 {:#x}); section '{}' op #{}.",
                g.id,
                g.hmr_caller,
                g.hmr_op_idx,
                r.id(),
                r_mask,
                g.mask_at_hmr,
                diff,
                section,
                last_op,
            );
        }
        // 3. Re-anchor r to the same id. The next on_z(r) under the
        //    condition stack [bit] will deposit Anchor(id) into the
        //    obligation's discharges (via discharge_to_cond).
        self.phase.re_anchor_for_resolve(r.id(), g.anchor_id);
        // 4. Discharge: z_if_bit(r, bit) := push_cond(bit), z(r),
        //    pop_cond. (circuit.rs:3519-3524.) Emits the actual
        //    phase-cancellation gate.
        self.z_if_bit(r, g.bit);
        // 5. Free the Cbit. Inside on_free_bit, the tracker
        //    notices obligation_is_clean(ob) (one Anchor(id)
        //    discharge matches the Anchor(id) obligation) and
        //    removes the entry — so this Cbit can be reused
        //    without a stale-obligation footgun.
        self.free_bit(g.bit);
        // Debugger ghost-event log (resolve).
        let resolve_op_idx = self.ops_truncated as usize + self.ops.len();
        self.ghost_event_log
            .push(crate::point_add::trailmix_port::circuit::GhostEvent {
                op_idx: resolve_op_idx,
                create: false,
                id: g.id,
                anchor_id: g.anchor_id,
                mask_at_hmr: 0,
                bit_raw: g.bit.raw(),
                section: self.current_section.clone(),
            });
        // 6. Drop the bookkeeping row + mark the Ghost consumed so
        //    its Drop is a no-op.
        self.pending_ghosts.retain(|rec| rec.id != g.id);
        g.consumed = true;
        // Ghost is moved-in by value; it drops at end of scope but
        // `consumed = true` short-circuits the panic branch.
    }

    // ── Multi-term (XOR-of-Z/CZ/CCZ) ghost discharge ───────────────
    //
    // For deferred phases that cancel only as a SUM of terms — e.g.
    // Gidney's borrowed-dirty constant adder, where a vented carry's
    // phase is corrected by `Z(dirty@v1) + Z(dirty@v2)` with
    // `dirty@v1 XOR dirty@v2 == carry`. Each `ghost_xor_*` emits one
    // phase gate gated on the ghost's bit and accumulates that term's
    // 64-shot sim mask; `close_ghost` requires the accumulated XOR to
    // equal the vented value's mask before clearing the obligation.

    /// Deposit a `Z(r)` discharge term toward ghost `g` (gated on
    /// `g.bit`), accumulating `r`'s 64-shot mask.
    #[track_caller]
    pub fn ghost_xor_z(&mut self, g: &mut Ghost, r: &QReg) {
        assert!(
            self.sim_condition_stack.is_empty(),
            "ghost_xor_z forbidden inside push_condition"
        );
        let m = self.sim_get_mask(r.id());
        g.acc_xor ^= m;
        self.z_if_bit(r, g.bit);
        self.phase.accumulate_obligation_discharge(g.bit.raw(), m);
    }

    /// Deposit a `CZ(r1, r2)` discharge term toward ghost `g` — the
    /// vented value's AND component. Accumulates `mask(r1) & mask(r2)`.
    #[track_caller]
    pub fn ghost_xor_cz(&mut self, g: &mut Ghost, r1: &QReg, r2: &QReg) {
        assert!(
            self.sim_condition_stack.is_empty(),
            "ghost_xor_cz forbidden inside push_condition"
        );
        let m = self.sim_get_mask(r1.id()) & self.sim_get_mask(r2.id());
        g.acc_xor ^= m;
        self.cz_if_bit(r1, r2, g.bit);
        self.phase.accumulate_obligation_discharge(g.bit.raw(), m);
    }

    /// Deposit a `CCZ(r1, r2, r3)` discharge term toward ghost `g`.
    /// Accumulates `mask(r1) & mask(r2) & mask(r3)`.
    #[track_caller]
    pub fn ghost_xor_ccz(&mut self, g: &mut Ghost, r1: &QReg, r2: &QReg, r3: &QReg) {
        assert!(
            self.sim_condition_stack.is_empty(),
            "ghost_xor_ccz forbidden inside push_condition"
        );
        let m =
            self.sim_get_mask(r1.id()) & self.sim_get_mask(r2.id()) & self.sim_get_mask(r3.id());
        g.acc_xor ^= m;
        self.ccz_if_bit(r1, r2, r3, g.bit);
        self.phase.accumulate_obligation_discharge(g.bit.raw(), m);
    }

    /// Finalize a multi-term ghost: require that the accumulated
    /// discharge XOR equals the vented value on every shot (the tracker
    /// checks `discharge_xor == hmr_mask`), then free the bit and drop
    /// the bookkeeping. Panics (with `DEBUG_ON_FAIL=1` attach) on
    /// mismatch — the phases would NOT cancel for a random bit.
    #[track_caller]
    pub fn close_ghost(&mut self, mut g: Ghost) {
        let caller = Location::caller();
        match self.phase.resolve_masked_obligation(g.bit.raw()) {
            Ok(()) => {}
            Err((hmr_mask, discharge_xor)) => {
                let last_op = self.ops_truncated as usize + self.ops.len();
                if std::env::var("DEBUG_ON_FAIL").is_ok() {
                    eprintln!(
                        "[close_ghost] discharge mismatch: discharge_xor={:#x} \
                         != hmr_mask={:#x} (diff={:#x}) Ghost #{} HMR at {} \
                         (op #{}); close at {} (op #{}). The Z/CZ/CCZ terms do \
                         not reproduce the vented value, so the phase would not \
                         cancel for a random HMR bit.",
                        discharge_xor,
                        hmr_mask,
                        discharge_xor ^ hmr_mask,
                        g.id,
                        g.hmr_caller,
                        g.hmr_op_idx,
                        caller,
                        last_op,
                    );
                    let mut dbg =
                        crate::point_add::trailmix_port::tracker::debugger::Debugger::attach(self);
                    dbg.goto(last_op);
                    dbg.repl();
                }
                g.consumed = true;
                panic!(
                    "close_ghost: discharge_xor={:#x} != hmr_mask={:#x} \
                     (Ghost #{} HMR at {} op #{}); acc_xor={:#x}",
                    discharge_xor, hmr_mask, g.id, g.hmr_caller, g.hmr_op_idx, g.acc_xor,
                );
            }
        }
        // Defense-in-depth: the Ghost's own accumulator must agree with
        // the tracker's (both XOR the same per-term masks).
        debug_assert_eq!(
            g.acc_xor, g.mask_at_hmr,
            "ghost acc_xor disagrees with hmr mask"
        );
        self.free_bit(g.bit);
        self.pending_ghosts.retain(|rec| rec.id != g.id);
        g.consumed = true;
    }

    /// Number of unresolved [`Ghost`]s currently held by this
    /// Circuit. Useful for asserting in test harnesses
    /// ("forward pass produced N ghosts; reverse pass should
    /// resolve all N before `assert_phase_clean`").
    #[must_use]
    pub fn pending_ghost_count(&self) -> usize {
        self.pending_ghosts.len()
    }

    /// Tear down a circuit, handing over the pending ghosts EXACTLY as
    /// `destroy_sim` hands over the live `QRegs` as `outputs`. Every pending ghost
    /// MUST be in `ghosts` (the wrapped `destroy_sim` asserts none remain), so
    /// nothing leaks silently. The handed-over ghosts are CONSUMED here: each
    /// bookkeeping row is dropped and the `Ghost` is marked consumed so its Drop
    /// is a no-op.
    ///
    /// For a real circuit you `resolve_ghost` every ghost first (emitting the
    /// phase cancellation), leaving none pending — then plain `destroy_sim`
    /// suffices. For an intentionally-incomplete COST-MEASUREMENT fragment (a
    /// forward pass whose spooky tape has no discharge yet) you pass the
    /// undischarged tape here; it is ABANDONED (NOT phase-discharged — the
    /// circuit is not phase-clean), which is the explicit, accounted-for way to
    /// measure peak/Toffoli without leaking the obligation.
    pub fn destroy_sim_ghosts(
        &mut self,
        outputs: Vec<crate::point_add::trailmix_port::circuit::QReg>,
        ghosts: Vec<Ghost>,
    ) -> (
        crate::point_add::trailmix_port::circuit::DestroyedSimState,
        Vec<crate::point_add::trailmix_port::circuit::QReg>,
    ) {
        // The handed-over ghosts must be the COMPLETE set of phase obligations:
        // assert the lattice is clean except for exactly these ghosts' HMR
        // obligations. If it is, then resolving all of them WOULD make the phase
        // clean -- so abandoning them is honest, and there is no un-ghosted phase
        // obligation hiding a real bug.
        let ghost_bits: std::collections::HashSet<u32> =
            ghosts.iter().map(|g| g.bit.raw()).collect();
        self.phase.assert_clean_except(&ghost_bits);
        for mut g in ghosts {
            self.pending_ghosts.retain(|rec| rec.id != g.id);
            g.consumed = true;
        }
        // destroy_sim asserts pending_ghosts is now empty (every ghost handed
        // over) and detaches the QRegs.
        self.destroy_sim(outputs)
    }
}

// ── Tests ──────────────────────────────────────────────────────
//
// Tracker-level smoke tests live in
// phase_lattice.rs (anchor lifecycle); the tests here exercise the
// full Circuit-side dance with a live 64-shot sim.

#[cfg(test)]
mod tests {
    use crate::point_add::trailmix_port::circuit::Circuit;

    /// Happy path (design §9.1, scoped to a single-bit input):
    ///   - allocate input bits a, b with random per-shot values
    ///   - compute q = a AND b (CCX), HMR-ghost q, zero+free q
    ///   - recompute r = a AND b (CCX), resolve_ghost(g, &r)
    ///   - uncompute r, free
    ///   - destroy_sim with [a, b] as outputs (avoids the
    ///     classical-load cleanup dance for inputs)
    ///   - assert_phase_clean passes; phase_mask() == 0.
    #[test]
    fn ghost_resolves_against_recomputed_witness() {
        let mut c = Circuit::new();
        c.set_section("ghost_happy");
        // Inputs: mark as F2 atoms so tracker sees them as Top
        // (not Zero from alloc). Without `alloc_input_qreg`,
        // q = ccx(a, b, q) tracker-collapses to Zero and the test
        // is trivially clean for the wrong reason.
        let a = c.alloc_input_qreg("a");
        let b = c.alloc_input_qreg("b");
        // Per-shot random bits — 64 independent (a, b) patterns
        // in one run, which is the spec's PBT-equivalent.
        use rand::Rng;
        let mut rng = rand::thread_rng();
        for shot in 0..64 {
            let av: u8 = rng.gen::<u8>() & 1;
            let bv: u8 = rng.gen::<u8>() & 1;
            c.sim_load_reg_bytes_shot(std::slice::from_ref(&a), &[av], shot);
            c.sim_load_reg_bytes_shot(std::slice::from_ref(&b), &[bv], shot);
        }

        let q = c.alloc_qreg("q_witness");
        c.ccx(&a, &b, &q); // q := a AND b

        let g = c.hmr_ghost(&q);
        assert_eq!(c.pending_ghost_count(), 1);
        c.zero_and_free(q); // q is post-HMR Zero on every shot

        // Some intervening ops just to prove the anchor survives
        // tracker version churn (a stand-in for the "forward pass
        // → reverse pass" gap in real use sites).
        let scratch = c.alloc_qreg("scratch");
        c.cx(&a, &scratch);
        c.cx(&a, &scratch);
        c.zero_and_free(scratch);

        // Compute a fresh witness register with the SAME identity.
        let r = c.alloc_qreg("r_recompute");
        c.ccx(&a, &b, &r); // r := a AND b matches q at HMR time

        c.resolve_ghost(g, &r);
        assert_eq!(c.pending_ghost_count(), 0);

        // Uncompute r so it reads |0> on every shot, then free.
        c.ccx(&a, &b, &r);
        c.zero_and_free(r);

        c.assert_phase_clean();

        // Drop a, b via destroy_sim: bypasses the classical
        // sim-cleanup which is init-only (sim_load_reg_bytes_shot
        // refuses to fire after the first gate emission). The
        // phase_mask is the snapshotted sim_phase.
        let (snap, _outs) = c.destroy_sim(vec![a, b]);
        assert_eq!(
            snap.phase_mask(),
            0,
            "sim_phase should be 0 after a clean ghost dance",
        );
    }

    /// Dropping a Ghost without `resolve_ghost` is a phase leak;
    /// the Ghost's own Drop impl panics with the HMR call site.
    #[test]
    #[should_panic(expected = "phase leak")]
    fn ghost_panics_on_unresolved_drop() {
        let mut c = Circuit::new();
        c.set_section("ghost_drop_leak");
        let q = c.alloc_qreg("q_drop");
        let g = c.hmr_ghost(&q);
        c.zero_and_free(q);
        drop(g); // panics: "Ghost ... dropped without resolve_ghost — phase leak."
    }

    /// resolve_ghost with a mismatched per-shot value panics with
    /// an `r_mask` diagnostic. Force the mismatch by loading `a`
    /// to all-ones, computing q := a via cx, then resolving
    /// against a fresh |0> register.
    #[test]
    #[should_panic(expected = "r_mask")]
    fn ghost_panics_on_wrong_resolve_target() {
        let mut c = Circuit::new();
        c.set_section("ghost_mismatch");
        let a = c.alloc_input_qreg("a");
        // All-ones across all 64 shots: a's mask = u64::MAX.
        for shot in 0..64 {
            c.sim_load_reg_bytes_shot(std::slice::from_ref(&a), &[1], shot);
        }
        let q = c.alloc_qreg("q");
        c.cx(&a, &q); // q := a
        let g = c.hmr_ghost(&q);
        // q's sim is now |0>; tracker also has q = Zero.
        c.zero_and_free(q);
        let r = c.alloc_qreg("r_wrong"); // |0> on every shot
        c.resolve_ghost(g, &r); // PANIC: r_mask=0 != captured u64::MAX
    }

    /// resolve_ghost under an active push_condition panics with
    /// the strict-semantic message from design §3.2.
    #[test]
    #[should_panic(expected = "condition stack")]
    fn ghost_resolve_inside_push_condition_panics() {
        let mut c = Circuit::new();
        c.set_section("ghost_under_cond");
        let q = c.alloc_qreg("q");
        let g = c.hmr_ghost(&q);
        c.zero_and_free(q);
        let r = c.alloc_qreg("r"); // |0>, matches q's |0> mask
        let cond = c.alloc_bit();
        c.push_condition(cond);
        c.resolve_ghost(g, &r); // PANIC: cond stack non-empty
    }

    /// hmr_ghost under an active push_condition panics with the
    /// strict-conditionality message (mirrors `hmr`'s check).
    #[test]
    #[should_panic(expected = "push_condition")]
    fn ghost_hmr_inside_push_condition_panics() {
        let mut c = Circuit::new();
        c.set_section("ghost_hmr_under_cond");
        let q = c.alloc_qreg("q");
        let cond = c.alloc_bit();
        c.push_condition(cond);
        let _g = c.hmr_ghost(&q); // PANIC
    }

    /// Circuit dropped while a forgotten Ghost is in the pending
    /// table panics with the Circuit-side diagnostic. Uses
    /// `mem::forget` to bypass Ghost::drop so the Circuit-level
    /// check is exclusively what fires.
    #[test]
    #[should_panic(expected = "unresolved Ghost")]
    fn circuit_panics_on_drop_with_pending_ghosts() {
        let mut c = Circuit::new();
        c.set_section("ghost_circ_drop");
        let q = c.alloc_qreg("q");
        let g = c.hmr_ghost(&q);
        c.zero_and_free(q);
        // SAFETY: deliberately bypass Ghost::drop to exercise the
        // Circuit-side fallback diagnostic. Don't do this in
        // production code.
        std::mem::forget(g);
        drop(c); // panics: "Circuit dropped with 1 unresolved Ghost(s) ..."
    }

    /// assert_phase_clean fires the Ghost-specific pre-check when
    /// a Ghost is still pending. Message includes the HMR location.
    #[test]
    #[should_panic(expected = "unresolved Ghost")]
    fn assert_phase_clean_reports_pending_ghosts() {
        let mut c = Circuit::new();
        c.set_section("ghost_apc_pending");
        let q = c.alloc_qreg("q");
        let g = c.hmr_ghost(&q);
        c.zero_and_free(q);
        std::mem::forget(g);
        c.assert_phase_clean(); // PANIC: pending_ghosts non-empty
                                // Don't let Circuit::drop fire its own panic during unwind:
                                // unreachable, but a guard for completeness if behaviour changes.
    }

    /// A contract can read the pending-ghost set: count, and each
    /// ghost's tape-bit mask per shot. This is the load-bearing
    /// assertion shape for windowed-tape spooky sweeps ("after the
    /// forward pass, exactly K ghosts pend and ghost i's tape bit ==
    /// the reference decision bit for step i").
    #[test]
    fn contract_reads_pending_ghosts() {
        let mut c = Circuit::new();
        c.set_section("ghost_contract");
        let a = c.alloc_input_qreg("a");
        let b = c.alloc_input_qreg("b");
        use rand::Rng;
        let mut rng = rand::thread_rng();
        let mut av = [0u8; 64];
        let mut bv = [0u8; 64];
        for shot in 0..64 {
            av[shot] = rng.gen::<u8>() & 1;
            bv[shot] = rng.gen::<u8>() & 1;
            c.sim_load_reg_bytes_shot(std::slice::from_ref(&a), &[av[shot]], shot);
            c.sim_load_reg_bytes_shot(std::slice::from_ref(&b), &[bv[shot]], shot);
        }
        let q = c.alloc_qreg("q");
        c.ccx(&a, &b, &q); // q := a AND b
        let g = c.hmr_ghost(&q);
        let gid = g.id;
        c.zero_and_free(q);

        // Contract: exactly one ghost pends, and its tape bit equals
        // a AND b on every shot.
        c.contract_check("one_ghost_tape", |view, shot| {
            if view.pending_ghost_count() != 1 {
                return Err(format!(
                    "expected 1 pending ghost, got {}",
                    view.pending_ghost_count()
                ));
            }
            let expect = (av[shot] & bv[shot]) == 1;
            match view.ghost_value_shot(gid, shot) {
                Some(v) if v == expect => Ok(()),
                other => Err(format!(
                    "ghost {} tape bit {:?} != expected {}",
                    gid, other, expect
                )),
            }
        });

        let r = c.alloc_qreg("r_recompute");
        c.ccx(&a, &b, &r); // r := a AND b matches q at HMR time
        c.resolve_ghost(g, &r);

        // After resolve, the ghost set is empty.
        c.contract_check("no_ghost", |view, _shot| {
            if view.pending_ghost_count() == 0 {
                Ok(())
            } else {
                Err(format!(
                    "{} ghosts still pending",
                    view.pending_ghost_count()
                ))
            }
        });

        c.ccx(&a, &b, &r);
        c.zero_and_free(r);
        c.assert_phase_clean();
        let (snap, _outs) = c.destroy_sim(vec![a, b]);
        assert_eq!(snap.phase_mask(), 0);
    }
}
