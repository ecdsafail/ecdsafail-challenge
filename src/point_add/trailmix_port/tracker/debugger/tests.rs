//! Tests for the time-travel debugger.

use super::*;
use crate::point_add::trailmix_port::circuit::{Cbit, Circuit};
use rand::Rng;

/// Build a random circuit, attach the debugger, replay to the end,
/// and assert the debugger's sim state matches circ.sim qubit-by-
/// qubit, bit-by-bit. Runs with many seeds / shapes / pre-gate
/// loads to catch replay divergences (initial_sim_state capture
/// timing, op-semantics drift, HMR counter misalignment, etc.).
///
/// Regression: this would have caught the "set_section snapshots
/// initial_sim_state" bug where user code emitted gates before
/// the first set_section, causing 8705-qubit divergence in a
/// pornin forward build.
fn differential_one_seed() {
    let mut rng = rand::thread_rng();
    let mut c = Circuit::new();
    // Allocate input + internal qubits and bits.
    let n_in = rng.gen_range(1..=8);
    let n_int = rng.gen_range(1..=8);
    let n_bits = rng.gen_range(0..=4);
    let inputs = c.alloc_qreg_bits("inputs", n_in);
    let internals = c.alloc_qreg_bits("internal", n_int);
    let bits: Vec<Cbit> = (0..n_bits).map(|_| c.alloc_input_bit()).collect();
    // Randomize input values.
    for q in &inputs {
        if rng.gen::<bool>() {
            let mut bytes = [0u8; 32];
            bytes[0] = 1;
            c.sim_load_reg_bytes_shot(std::slice::from_ref(q), &bytes, 0);
        }
    }
    for &b in &bits {
        if rng.gen::<bool>() {
            let mut bytes = [0u8; 32];
            bytes[0] = 1;
            c.sim_load_bits_bytes(&[b], &bytes);
        }
    }

    // Track all qubits via index. inputs come first (0..n_in),
    // then internals (n_in..n_in+n_int).
    let n_total: usize = n_in + n_int;
    let all_b: Vec<Cbit> = bits.clone();

    // Sometimes call set_section BEFORE any gates; sometimes not.
    let use_early_section = rng.gen::<bool>();
    if use_early_section {
        c.set_section("start");
    }

    // Helper: pick a qubit index, then resolve to &QReg. Returning
    // &'a QReg with explicit lifetime ties the borrow to whichever
    // input vector owns the picked qubit, so the borrow checker
    // accepts it without needing a raw-pointer detour.
    fn qref<'a>(
        idx: usize,
        inputs: &'a [crate::point_add::trailmix_port::circuit::QReg],
        internals: &'a [crate::point_add::trailmix_port::circuit::QReg],
    ) -> &'a crate::point_add::trailmix_port::circuit::QReg {
        if idx < inputs.len() {
            &inputs[idx]
        } else {
            &internals[idx - inputs.len()]
        }
    }

    // Emit a random sequence of gates.
    let n_ops = rng.gen_range(5..100);
    for _ in 0..n_ops {
        let kind = rng.gen_range(0..10);
        let pick = |r: &mut rand::rngs::ThreadRng| r.gen_range(0..n_total);
        match kind {
            0 => {
                let q = pick(&mut rng);
                c.x(qref(q, &inputs, &internals));
            }
            1 => {
                let q = pick(&mut rng);
                c.z(qref(q, &inputs, &internals));
            }
            2 => {
                let a = pick(&mut rng);
                let b_cands: Vec<usize> = (0..n_total).filter(|&x| x != a).collect();
                if !b_cands.is_empty() {
                    let b = b_cands[rng.gen_range(0..b_cands.len())];
                    c.cx(qref(a, &inputs, &internals), qref(b, &inputs, &internals));
                }
            }
            3 => {
                let a = pick(&mut rng);
                let b_cands: Vec<usize> = (0..n_total).filter(|&x| x != a).collect();
                if !b_cands.is_empty() {
                    let b = b_cands[rng.gen_range(0..b_cands.len())];
                    c.cz(qref(a, &inputs, &internals), qref(b, &inputs, &internals));
                }
            }
            4 => {
                if n_total >= 3 {
                    let i = rng.gen_range(0..n_total);
                    let mut j = rng.gen_range(0..n_total);
                    while j == i {
                        j = rng.gen_range(0..n_total);
                    }
                    let mut k = rng.gen_range(0..n_total);
                    while k == i || k == j {
                        k = rng.gen_range(0..n_total);
                    }
                    c.ccx(
                        qref(i, &inputs, &internals),
                        qref(j, &inputs, &internals),
                        qref(k, &inputs, &internals),
                    );
                }
            }
            5 => {
                let a = pick(&mut rng);
                let b_cands: Vec<usize> = (0..n_total).filter(|&x| x != a).collect();
                if !b_cands.is_empty() {
                    let b = b_cands[rng.gen_range(0..b_cands.len())];
                    c.swap(qref(a, &inputs, &internals), qref(b, &inputs, &internals));
                }
            }
            6 if !all_b.is_empty() => {
                let q = pick(&mut rng);
                let b = all_b[rng.gen_range(0..all_b.len())];
                c.hmr(qref(q, &inputs, &internals), b);
            }
            7 => {
                // X twice (no-op). Tests consecutive same ops.
                let q = pick(&mut rng);
                c.x(qref(q, &inputs, &internals));
                c.x(qref(q, &inputs, &internals));
            }
            8 if rng.gen::<bool>() => c.neg(),
            _ => {
                let q = pick(&mut rng);
                c.x(qref(q, &inputs, &internals));
            }
        }
    }

    // Attach debugger, goto end, diff. After Debugger::attach, drop
    // the QReg vectors so the Circuit's pending-frees queue is empty
    // by the time destroy_sim runs (alloc_qreg_bits qregs would
    // otherwise hold strong Rcs into the circuit).
    let mut d = Debugger::attach(&mut c);
    d.goto(d.num_ops()); // c.ops was moved into the debugger by attach; use d.num_ops()
    let total_q = c.total_qubits() as usize;
    let total_b_count = c.total_bits() as usize;
    // Convert to u64 phase before destroy.
    let pre_phase = d.phase();
    // Destroy sim consuming the Vec of inputs+internals so the QRegs
    // in the test bound to the same circuit can be dropped after.
    // Concatenate into one Vec<QReg> for destroy_sim.
    let mut all_qregs: Vec<crate::point_add::trailmix_port::circuit::QReg> = Vec::with_capacity(n_total);
    all_qregs.extend(inputs);
    all_qregs.extend(internals);
    let (sim, _surviving) = c.destroy_sim(all_qregs);
    for q in 0..total_q {
        let real = sim[q];
        let dbg_v = d.qubit(q as u32);
        assert_eq!(
            real, dbg_v,
            "q{}: real={:#x} dbg={:#x} (early_section={})",
            q, real, dbg_v, use_early_section
        );
    }
    for b in 0..total_b_count {
        let real = sim.bit_mask(b as u32);
        let dbg_v = d.bit(b as u32);
        assert_eq!(real, dbg_v, "b{}: real={:#x} dbg={:#x}", b, real, dbg_v);
    }
    assert_eq!(
        sim.phase_mask(),
        pre_phase,
        "phase mismatch (real={:#x} dbg={:#x})",
        sim.phase_mask(),
        pre_phase
    );
}

#[test]
fn differential_debugger_matches_sim() {
    // Run 200 random circuits. Half use early set_section, half
    // don't — both should diff-match after the circuit.rs fix that
    // moved initial_sim_state capture to first gate emission.
    for _ in 0..200 {
        differential_one_seed();
    }
}

/// Regression test for the specific "gates-before-set_section"
/// bug: lazily captured initial_sim_state must be taken BEFORE
/// any gate's sim modification. Emit gates + sim_load_reg_bytes
/// BEFORE any set_section, then attach debugger and verify the
/// replay matches.
#[test]
fn debugger_captures_initial_state_before_any_gate() {
    let mut c = Circuit::new();
    let y = c.alloc_qreg("y");
    let a = c.alloc_qreg("a");
    // Load y = 1.
    let mut bytes = [0u8; 32];
    bytes[0] = 1;
    c.sim_load_reg_bytes_shot(std::slice::from_ref(&y), &bytes, 0);
    // Emit a CX BEFORE any set_section. Historically this lost
    // the initial y=1 load from the debugger's replay.
    c.cx(&y, &a);
    // Now add a set_section (lazy snapshot would fire here if
    // we hadn't fixed it).
    c.set_section("post_cx");
    c.cx(&a, &y); // any further gate

    let mut d = Debugger::attach(&mut c);
    d.goto(d.num_ops()); // c.ops was moved into the debugger by attach; use d.num_ops()
    let total_q = c.total_qubits() as usize;
    let (sim, _) = c.destroy_sim(vec![y, a]);
    for q in 0..total_q {
        let real = sim[q];
        let dbg_v = d.qubit(q as u32);
        assert_eq!(real, dbg_v, "q{}: real={:#x} dbg={:#x}", q, real, dbg_v);
    }
}

/// Build a small reversible dance and confirm the debugger's
/// forward/backward reach identical intermediate states.
#[test]
fn round_trip_step_back() {
    let mut c = Circuit::new();
    // QRegs allocated in order: a → id 0, b → id 1, anc → id 2.
    let a = c.alloc_qreg("a");
    let b = c.alloc_qreg("b");
    let anc = c.alloc_qreg("anc");
    let (a_id, b_id, anc_id) = (0u32, 1u32, 2u32);
    c.sim_load_reg_bytes_shot(std::slice::from_ref(&a), &[1], 0); // a = 1 in shot 0
    c.sim_load_reg_bytes_shot(std::slice::from_ref(&b), &[1], 0); // b = 1
    c.set_section("dance");
    c.ccx(&a, &b, &anc);
    c.hmr(&anc, Cbit(c.peak_bits)); // nonsense bit, but tests op
    c.x(&a);
    c.x(&a);

    // Compute end-state values via debugger BEFORE destroy_sim
    // (debugger holds a snapshot independent of the live sim).
    let mut d = Debugger::attach(&mut c);
    d.goto(d.num_ops()); // c.ops was moved into the debugger by attach; use d.num_ops()
    let end_phase = d.phase();
    let end_q_a = (d.qubit(a_id) & 1) as u8;
    // Cursor at end by default.
    assert_eq!(d.cursor(), d.num_ops());

    // Rewind to start, step forward, back, and compare checkpoints.
    d.goto(0);
    assert_eq!(d.cursor(), 0);
    // Walk forward to midway, record state, continue, rewind to mid,
    // and verify state matches.
    d.step(1);
    let mid_cursor = d.cursor();
    let mid_q_a = d.qubit(a_id);
    let mid_q_b = d.qubit(b_id);
    let mid_anc = d.qubit(anc_id);
    let mid_phase = d.phase();
    d.goto(d.num_ops());
    d.goto(mid_cursor);
    assert_eq!(d.qubit(a_id), mid_q_a);
    assert_eq!(d.qubit(b_id), mid_q_b);
    assert_eq!(d.qubit(anc_id), mid_anc);
    assert_eq!(d.phase(), mid_phase);
    // Sanity-check end values match what the live sim produced.
    let (sim, _) = c.destroy_sim(vec![a, b, anc]);
    assert_eq!(sim.phase_mask(), end_phase);
    assert_eq!((sim[a_id as usize] & 1) as u8, end_q_a);
}

#[test]
fn phase_breakpoint_fires_at_neg() {
    let mut c = Circuit::new();
    c.set_section("main");
    c.neg(); // flips all 64 shots
    let mut d = Debugger::attach(&mut c);
    d.goto(0);
    d.add_breakpoint(Breakpoint::PhaseNonzero);
    let hit = d.run_until_break();
    assert!(matches!(hit, Some(Breakpoint::PhaseNonzero)));
    assert_ne!(d.phase(), 0);
}

/// Exercises every op variant forward-then-back and asserts the
/// final state matches a fresh snapshot-only replay.
#[test]
fn delta_back_matches_snapshot_replay() {
    let mut c = Circuit::new();
    let a = c.alloc_qreg("a");
    let b = c.alloc_qreg("b");
    let t = c.alloc_qreg("t");
    c.sim_load_reg_bytes_shot(std::slice::from_ref(&a), &[0x5a], 0);
    c.sim_load_reg_bytes_shot(std::slice::from_ref(&b), &[0xa5], 0);
    let flag = c.alloc_bit();
    c.set_section("mix");
    // HMR + Push/Pop + BitStore + self-inverse ops. (R was
    // formerly emitted via c.r(t); zero_and_free now is the only
    // public path and it consumes the QReg, so the trailing R is
    // emitted by dropping/freeing t cleanly; the test still walks
    // every other op variant.)
    c.ccx(&a, &b, &t);
    c.cx(&a, &t);
    c.z(&a);
    c.cz(&a, &b);
    c.swap(&a, &b);
    c.neg();
    c.bit_store1(flag);
    c.push_condition(flag);
    c.x(&t);
    c.bit_invert(flag);
    c.pop_condition();
    c.bit_store0(flag);
    // Clean t back to |0>: undo the c.x(&t) that fired with cond=1
    // (bit_store1 set flag pre-push). After pop, t state can be
    // computed via sim. Easiest: x(t) once more so the conditional
    // x's parity is even when bit_store1 was applied.
    // Simpler: emit c.x(&t) under a no-op sequence to put t back to
    // its original computed state, then zero_and_free which checks |0>.
    // To avoid arithmetic gymnastics, just emit an X if needed.
    // The combined sequence (ccx + cx + conditional-x + ...) leaves
    // t in some state; the test's role is only to exercise the
    // op variants, not enforce |0>. Use r_internal-equivalent path:
    // emit a fresh ancilla zero-and-free is the only public R.
    // Pad t to |0> manually:
    // After ccx(a,b,t): t = (a&b). a=0x5a, b=0xa5; a&b=0; t=0.
    // After cx(a,t): t ^= a. t = 0x5a.
    // swap(a,b): a<->b; a=0xa5, b=0x5a; t unchanged 0x5a.
    // x(t) under cond=1 (flag=1, then bit_invert flips to 0; the x
    // executes under flag=1 only): t ^= 1 in shot 0 (bit 0 was 0):
    // 0x5a ^ 1 = 0x5b.
    // After pop_condition + bit_store0: flag=0.
    // To zero t: x_if_bit isn't available; we just XOR t by 0x5b
    // unconditionally to bring it to |0>.
    // 0x5b = 0101_1011. Apply X to bits 0,1,3,4,6 of t (single qubit).
    // But t is a single QReg so its low bit is 0x5b & 1 = 1.
    // For a single qubit t, value across shots is 0x5b. We x-ate
    // many times; the cleanest way is to load via sim and just R.
    // Use the R via the destroy_sim path: append x(&t) iff sim
    // reports t == 1, then zero_and_free.
    // Simplest path: just give up cleaning t and detach it before
    // destroy_sim by passing it as an output to destroy_sim.
    let hmr_bit = c.alloc_bit();
    c.x(&a);
    c.hmr(&a, hmr_bit);

    let mut d = Debugger::attach(&mut c);
    // Snapshot-only reference: rebuild each state via restore+replay.
    let reference_states: Vec<SimState> = (0..=d.num_ops())
        .map(|k| {
            d.restore_snapshot_and_replay(k);
            d.state.clone()
        })
        .collect();

    // Now walk forward op-by-op, then back op-by-op; compare to ref.
    d.restore_snapshot_and_replay(0);
    assert_eq!(d.cursor(), 0);
    assert!(d.delta_log.is_empty());
    // Forward walk.
    for k in 1..=d.num_ops() {
        d.forward_one();
        assert_eq!(d.cursor(), k);
        assert_eq!(
            d.state.qubits, reference_states[k].qubits,
            "qubit mismatch forward k={}",
            k
        );
        assert_eq!(
            d.state.bits, reference_states[k].bits,
            "bit mismatch forward k={}",
            k
        );
        assert_eq!(
            d.state.phase, reference_states[k].phase,
            "phase mismatch forward k={}",
            k
        );
        assert_eq!(d.state.hmr_counter, reference_states[k].hmr_counter);
        assert_eq!(
            d.state.r_on_nonzero_events,
            reference_states[k].r_on_nonzero_events
        );
        assert_eq!(d.state.cond_stack, reference_states[k].cond_stack);
    }
    // Delta-log backward walk to 0.
    while d.cursor() > 0 {
        let k = d.cursor();
        d.back_one_via_log();
        assert_eq!(d.cursor(), k - 1);
        assert_eq!(
            d.state.qubits,
            reference_states[k - 1].qubits,
            "qubit mismatch reverse after k={}->{}",
            k,
            k - 1
        );
        assert_eq!(
            d.state.bits,
            reference_states[k - 1].bits,
            "bit mismatch reverse after k={}->{}",
            k,
            k - 1
        );
        assert_eq!(
            d.state.phase,
            reference_states[k - 1].phase,
            "phase mismatch reverse after k={}->{}",
            k,
            k - 1
        );
        assert_eq!(d.state.hmr_counter, reference_states[k - 1].hmr_counter);
        assert_eq!(
            d.state.r_on_nonzero_events,
            reference_states[k - 1].r_on_nonzero_events
        );
        assert_eq!(d.state.cond_stack, reference_states[k - 1].cond_stack);
    }
    // Detach the live qregs so the strict-dealloc check inside
    // Drop doesn't fire on residue values left in t.
    let _ = c.destroy_sim(vec![a, b, t]);
}

/// Small ring capacity: forward past cap, then back should fall
/// through to snapshot+replay without losing correctness.
#[test]
fn delta_log_ring_drop_recovers_via_snapshot() {
    let mut c = Circuit::new();
    let q = c.alloc_qreg("q");
    c.sim_load_reg_bytes_shot(std::slice::from_ref(&q), &[1], 0);
    c.set_section("spin");
    // 50 X(q) all on the same qubit: the redundant-op
    // auto-eliminator cancels them in pairs, leaving 25 elided
    // (None) slots in self.ops and a net q value unchanged from
    // the sim_load_reg_bytes_shot above. We verify the cancel
    // happened and exercise the snapshot+replay machinery.
    for _ in 0..50 {
        c.x(&q);
    }
    let emitted_x: usize = c
        .ops
        .iter()
        .filter(|op| matches!(op, Some(Op::X(_))))
        .count();
    assert_eq!(
        emitted_x, 0,
        "50 X(q) on the same qubit must auto-elide to 0 emitted X ops"
    );

    let mut d = Debugger::attach(&mut c);
    d.delta_log_cap = 8; // force many ring drops
    d.restore_snapshot_and_replay(0);
    // Walk forward to the end; log can only hold 8.
    let total = d.num_ops();
    while d.cursor() < total {
        d.forward_one();
    }
    assert_eq!(d.cursor(), total);
    assert!(d.delta_log.len() <= 8);
    // Request a back to the start — forces snapshot replay.
    d.goto(0);
    assert_eq!(d.cursor(), 0);
    // Walk forward again; sim value of q at the final cursor must
    // match the loaded value (50 X's cancelled to identity).
    while d.cursor() < total {
        d.forward_one();
    }
    // q is the first (and only) qubit allocated, so its index in
    // state.qubits is 0.
    assert_eq!(
        d.state.qubits[0] & 1,
        1,
        "q's shot 0 was loaded with 1 and 50 X's must cancel to identity"
    );
    let _ = c.destroy_sim(vec![q]);
}

#[test]
fn section_breakpoint_fires() {
    let mut c = Circuit::new();
    c.set_section("a");
    let qa = c.alloc_qreg("qa");
    c.x(&qa);
    c.set_section("b");
    let qb = c.alloc_qreg("qb");
    c.x(&qb);
    let mut d = Debugger::attach(&mut c);
    d.goto(0);
    d.add_breakpoint(Breakpoint::SectionStart("b".into()));
    let hit = d.run_until_break();
    assert!(matches!(hit, Some(Breakpoint::SectionStart(_))));
    assert_eq!(d.current_section(), "b");
    let _ = c.destroy_sim(vec![qa, qb]);
}

#[test]
fn profiler_supports_exact_prefix_current_and_split_queries() {
    let mut c = Circuit::new();
    let a = c.alloc_qreg("a");
    let b = c.alloc_qreg("b");
    let t = c.alloc_qreg("t");

    c.set_section("root/a/x");
    c.x(&a);
    c.ccx(&a, &b, &t);

    c.set_section("root/a/y");
    c.cx(&a, &b);
    c.ccx(&a, &b, &t);

    c.set_section("root/b/z");
    c.x(&b);

    let mut d = Debugger::attach(&mut c);
    d.goto(d.num_ops()); // c.ops was moved into the debugger by attach; use d.num_ops()

    let exact = d.print_profile(&["whole", "exact", "root/a/x"]);
    assert!(exact.contains("root/a/x"), "{}", exact);
    assert!(exact.contains("matched=1"), "{}", exact);

    let prefix = d.print_profile(&["whole", "prefix", "root/a"]);
    assert!(prefix.contains("root/a/x"), "{}", prefix);
    assert!(prefix.contains("root/a/y"), "{}", prefix);

    let current = d.print_profile(&["whole", "current"]);
    assert!(current.contains("root/b/z"), "{}", current);

    let split = d.print_profile(&["whole", "split", "root"]);
    assert!(split.contains("root/a"), "{}", split);
    assert!(split.contains("root/b"), "{}", split);

    let tof = d.print_profile(&["whole", "tof", "top", "5"]);
    // 2 CCX in the circuit; Toffoli count is a per-circuit property, not ×64 shots.
    assert!(tof.contains("total_tof=2"), "{}", tof);
    let _ = c.destroy_sim(vec![a, b, t]);
}

#[test]
fn cursor_profiler_tracks_step_back_and_goto() {
    let mut c = Circuit::new();
    let a = c.alloc_qreg("a");
    let b = c.alloc_qreg("b");
    let t = c.alloc_qreg("t");

    c.set_section("root/a/x");
    c.x(&a);
    c.ccx(&a, &b, &t);

    c.set_section("root/a/y");
    c.cx(&a, &b);
    c.ccx(&a, &b, &t);

    c.set_section("root/b/z");
    c.x(&b);

    let mut d = Debugger::attach(&mut c);
    d.goto(0);
    let empty = d.print_profile(&["cursor", "top", "5"]);
    assert!(empty.contains("total_ops=0"), "{}", empty);

    d.step(2);
    let mid = d.print_profile(&["cursor", "prefix", "root/a"]);
    assert!(mid.contains("root/a/x"), "{}", mid);
    assert!(mid.contains("matched_ops=2"), "{}", mid);

    d.back(1);
    let back = d.print_profile(&["cursor", "prefix", "root/a"]);
    assert!(back.contains("matched_ops=1"), "{}", back);
    assert!(!back.contains("matched_ops=2"), "{}", back);

    d.goto(d.num_ops()); // c.ops was moved into the debugger by attach; use d.num_ops()
    let end = d.print_profile(&["cursor", "split", "root"]);
    assert!(end.contains("root/a"), "{}", end);
    assert!(end.contains("root/b"), "{}", end);
    assert!(end.contains("total_ops=5"), "{}", end);
    let _ = c.destroy_sim(vec![a, b, t]);
}

#[test]
fn peak_profiler_supports_top_exact_prefix_and_contains() {
    let mut c = Circuit::new();
    c.set_section("shell/inv");
    let a0 = c.alloc_qreg("row_add_spill");
    let a1 = c.alloc_qreg("row_add_spill");
    let a2 = c.alloc_qreg("carry_push_many");
    c.set_section("shell/mul");
    let b0 = c.alloc_qreg("carry_push_many");

    let mut d = Debugger::attach(&mut c);
    let top = d.print_profile(&["peak", "top", "10"]);
    assert!(top.contains("total_qubits=4"), "{}", top);

    let exact = d.print_profile(&["peak", "exact", "shell/inv/row_add_spill"]);
    assert!(exact.contains("shell/inv/row_add_spill"), "{}", exact);
    assert!(exact.contains("matched_qubits=2"), "{}", exact);

    let exact_mul = d.print_profile(&["peak", "exact", "shell/mul/carry_push_many"]);
    assert!(
        exact_mul.contains("shell/mul/carry_push_many"),
        "{}",
        exact_mul
    );
    assert!(exact_mul.contains("matched_qubits=1"), "{}", exact_mul);

    let prefix = d.print_profile(&["peak", "prefix", "shell"]);
    assert!(prefix.contains("shell/inv/row_add_spill"), "{}", prefix);
    assert!(prefix.contains("shell/inv/carry_push_many"), "{}", prefix);
    assert!(prefix.contains("shell/mul/carry_push_many"), "{}", prefix);

    let contains = d.print_profile(&["peak", "contains", "mul"]);
    assert!(
        contains.contains("shell/mul/carry_push_many"),
        "{}",
        contains
    );
    let _ = c.destroy_sim(vec![a0, a1, a2, b0]);
}

#[test]
fn source_trace_records_scope_chain() {
    // Exercise enter_scope!/exit_scope and verify the debugger
    // can render a nested scope chain for an op inside. Use
    // distinct qubits per X so the redundant-op auto-elide
    // doesn't cancel them away (which would shift op_idx values
    // and break the captured inner_op_idx_mid).
    let mut c = Circuit::new();
    let q1 = c.alloc_qreg("q1");
    let q2 = c.alloc_qreg("q2");
    let q3 = c.alloc_qreg("q3");

    let outer = crate::enter_scope!(&mut c, "outer_block");
    let outer_file = file!();
    c.x(&q1);
    let inner_op_idx_mid;
    {
        let inner = crate::enter_scope!(&mut c, "inner_block");
        inner_op_idx_mid = c.ops.len(); // next op's idx
        c.x(&q2);
        c.exit_scope(inner);
    }
    c.x(&q3);
    c.exit_scope(outer);

    let d = Debugger::attach(&mut c);
    let trace = d.source_trace(inner_op_idx_mid);
    assert!(trace.contains("inner_block"), "trace: {}", trace);
    assert!(trace.contains("outer_block"), "trace: {}", trace);
    assert!(trace.contains(outer_file), "trace: {}", trace);
    let _ = c.destroy_sim(vec![q1, q2, q3]);
}

/// Regression test for the broken `src`/`source` debugger command.
///
/// Before the fix, `scope_frames_log` was only populated by manual
/// `enter_scope!` calls, which NO emission code ever made — so in a
/// real `--release` run `op_scope` was entirely `u32::MAX` and
/// `src <op>` always reported "no scope frame at local idx ... (use
/// enter_scope! to tag)", making the single most important debugging
/// command useless.
///
/// The fix makes the gate-emission boundary `#[track_caller]` and
/// auto-interns each distinct caller `(file, line)` into
/// `scope_frames_log`, so EVERY emitted op resolves to the source
/// site that emitted it with zero manual tagging. This test emits
/// gates from two distinct call sites (helpers `emit_x` / `emit_cx`)
/// with NO `enter_scope!` anywhere and asserts:
///   * `op_scope[i] != u32::MAX` for every op (auto frame present),
///   * `scope_frames_log` is non-empty after attach,
///   * `src <op>` returns a real `file:line`, NOT the old
///     "no scope frame" message,
///   * ops from different call sites get different frames / lines.
#[test]
fn src_resolves_file_line_without_manual_enter_scope() {
    #[track_caller]
    fn emit_x(c: &mut Circuit, q: &crate::point_add::trailmix_port::circuit::QReg) {
        c.x(q);
    }
    #[track_caller]
    fn emit_cx(c: &mut Circuit, a: &crate::point_add::trailmix_port::circuit::QReg, b: &crate::point_add::trailmix_port::circuit::QReg) {
        c.cx(a, b);
    }

    let mut c = Circuit::new();
    let q1 = c.alloc_qreg("q1");
    let q2 = c.alloc_qreg("q2");

    // No enter_scope! anywhere — exactly the real-circuit situation.
    let x_op_idx = c.ops.len();
    emit_x(&mut c, &q1); // call site A (this file, this line)
    let cx_op_idx = c.ops.len();
    emit_cx(&mut c, &q1, &q2); // call site B (this file, a different line)

    // Every op got an auto-captured frame (the old bug left these
    // all u32::MAX).
    assert_eq!(c.op_scope.len(), c.ops.len());
    for (i, &s) in c.op_scope.iter().enumerate() {
        assert_ne!(s, u32::MAX, "op {} has no scope frame (regression)", i);
    }
    // Distinct call sites get distinct interned frames.
    assert_ne!(
        c.op_scope[x_op_idx], c.op_scope[cx_op_idx],
        "distinct call sites should intern to distinct frames"
    );

    let d = Debugger::attach(&mut c);
    assert!(
        !d.scope_frames_log.is_empty(),
        "scope_frames_log empty after attach — `src` would be dead"
    );

    let this_file = file!();
    let x_trace = d.source_trace(x_op_idx);
    assert!(
        !x_trace.contains("no scope frame"),
        "src still reports no scope frame (regression): {}",
        x_trace
    );
    assert!(
        x_trace.contains(this_file),
        "src did not resolve to this file: {}",
        x_trace
    );
    assert!(
        x_trace.contains("scope chain"),
        "src output missing scope chain header: {}",
        x_trace
    );

    let cx_trace = d.source_trace(cx_op_idx);
    assert!(
        cx_trace.contains(this_file) && !cx_trace.contains("no scope frame"),
        "src did not resolve cx op to file:line: {}",
        cx_trace
    );
    // The two call sites are on different source lines, so the
    // rendered file:line scope chains must differ.
    assert_ne!(
        x_trace, cx_trace,
        "distinct call sites rendered identical scope chains"
    );

    let _ = c.destroy_sim(vec![q1, q2]);
}
