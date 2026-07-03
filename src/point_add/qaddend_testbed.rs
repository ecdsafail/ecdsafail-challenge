//! Tier B (issue #27/#28), quantum-addend point-add TESTBED — first increment.
//!
//! ADRs 0011/0012/0013 all converge on one deferred build: a *functionally
//! correct* QROM-fed **quantum** addend. ADR 0012 showed the scored PA already
//! runs an uncontrolled q-q add over an addend it loads into qubits; ADR 0013
//! measured the width cost of holding that addend resident. What was still
//! missing is the **composition itself** — read `P[k]` from a quantum table INTO
//! an addend register, have an adder consume it, and unread — proven to compute
//! the right sum with every ancilla returned to |0>. ADR 0014 decides to build it
//! self-contained, width-parametric, and verify by simulation.
//!
//! This harness composes, on fresh registers:
//!   1. **QROM read** — the same unary-iteration selector (Gidney 2018 §III.C) as
//!      `ladder_lookup_cost.py` (ADR 0010) and `ladder_full.rs`, but now WITH the
//!      leaf **data-writes** those cost-only harnesses omit. Loads the addend
//!      register `addend := P[k]` from a classical table, addressed by the quantum
//!      window register `k`. Selector Toffoli is `2^(w+1)−4` (asserted — ties the
//!      construction back to ADR 0010).
//!   2. **q-q add** — an uncontrolled Cuccaro ripple-carry `acc += addend`
//!      (addend preserved), the exact shape `coord_addsub` uses (ADR 0012).
//!   3. **QROM unread** — a second application of the reversible read, returning
//!      `addend` to |0>.
//!
//! Verified by masked multi-shot simulation over ALL `2^w` window values and
//! several accumulator inputs: `acc' == (acc + P[k]) mod 2^n`, and the addend,
//! selector ancilla, carry, and window register all clean/preserved. The peak
//! width and its register breakdown are reported — the small-scale, executable
//! analogue of ADR 0013's full-width +256..512.
//!
//! v1 scope (ADR 0014): integer add mod `2^n` (the field-modular reduction tail is
//! a fixed, addend-independent extension, as `coord_addsub`'s `mod_add` shows) and
//! a single read→add→unread (not the full ladder). Issue #28's EC exceptional
//! cases (`P==Q`, `dx=0`, ∞) need the group law on top — the next increments on
//! this testbed.
//!
//! `#[cfg(test)]` only; never compiled into the scored circuit.

use crate::circuit::{analyze_ops, Op, OperationType, QubitId};
use crate::point_add::B;
use crate::sim::Simulator;

// ── circuit fragments (built on the real B op-emitter) ──────────────────────

/// Cuccaro MAJ: propagate the carry chain. `c` = carry-in line, `b` = accumulator
/// bit, `a` = addend bit (temporarily holds the outgoing carry).
fn maj(circ: &mut B, c: QubitId, b: QubitId, a: QubitId) {
    circ.cx(a, b);
    circ.cx(a, c);
    circ.ccx(c, b, a);
}

/// Cuccaro UMA: unwind the carry chain, restoring `a` and finishing `b`.
fn uma(circ: &mut B, c: QubitId, b: QubitId, a: QubitId) {
    circ.ccx(c, b, a);
    circ.cx(a, c);
    circ.cx(c, b);
}

/// `acc += addend (mod 2^n)`, **addend preserved**, `carry` a single |0> ancilla
/// (returned to |0>). Uncontrolled q-q ripple-carry — the `coord_addsub` shape.
fn cuccaro_add(circ: &mut B, addend: &[QubitId], acc: &[QubitId], carry: QubitId) {
    let n = addend.len();
    assert_eq!(acc.len(), n);
    maj(circ, carry, acc[0], addend[0]);
    for i in 1..n {
        maj(circ, addend[i - 1], acc[i], addend[i]);
    }
    // mod 2^n: the carry-out (would be cx(addend[n-1], cout)) is intentionally dropped.
    for i in (1..n).rev() {
        uma(circ, addend[i - 1], acc[i], addend[i]);
    }
    uma(circ, carry, acc[0], addend[0]);
}

/// Controlled-X with an optional control (top level of the unary iteration has
/// none — the always-on root).
fn cx_opt(circ: &mut B, ctrl: Option<QubitId>, tgt: QubitId) {
    match ctrl {
        None => circ.x(tgt),
        Some(c) => circ.cx(c, tgt),
    }
}
fn and_into(circ: &mut B, ctrl: Option<QubitId>, addr_bit: QubitId, anc: QubitId) {
    match ctrl {
        None => circ.cx(addr_bit, anc),
        Some(c) => circ.ccx(c, addr_bit, anc),
    }
}

/// Emit a full unary-iteration QROM read: `addend ^= table[k]`, where `k` is the
/// value of the quantum `win` register. `anc[level]` is the one-ancilla-per-level
/// selector spine; at the leaf the current spine line is the one-hot signal for
/// the resolved address, under which the classical constant's set bits are `CX`ed
/// into `addend`. Self-inverse on `addend` (applying it twice clears the addend),
/// and self-uncomputing on `anc`.
#[allow(clippy::too_many_arguments)]
fn qrom_read(circ: &mut B, win: &[QubitId], anc: &[QubitId], addend: &[QubitId], table: &[u64]) {
    fn rec(
        circ: &mut B,
        level: usize,
        ctrl: Option<QubitId>,
        addr: usize,
        win: &[QubitId],
        anc: &[QubitId],
        addend: &[QubitId],
        table: &[u64],
    ) {
        if level == win.len() {
            let val = table[addr];
            for (bit, &q) in addend.iter().enumerate() {
                if (val >> bit) & 1 == 1 {
                    // leaf ctrl is always Some for w >= 1
                    cx_opt(circ, ctrl, q);
                }
            }
            return;
        }
        let a = anc[level];
        let wq = win[level];
        // a = ctrl AND win[level]
        and_into(circ, ctrl, wq, a);
        rec(
            circ,
            level + 1,
            Some(a),
            addr | (1 << level),
            win,
            anc,
            addend,
            table,
        );
        // a = ctrl AND NOT win[level]
        cx_opt(circ, ctrl, a);
        rec(circ, level + 1, Some(a), addr, win, anc, addend, table);
        // restore a = ctrl AND win[level], then uncompute a -> 0
        cx_opt(circ, ctrl, a);
        and_into(circ, ctrl, wq, a);
    }
    // Shape invariants: turn silent OOB/shift panics into clear failures if this
    // reusable helper is ever called with inconsistent registers.
    assert!(
        win.len() <= 32,
        "qrom_read: window width {} unreasonably large (2^w table)",
        win.len()
    );
    assert_eq!(
        anc.len(),
        win.len(),
        "qrom_read: need one selector ancilla per window bit ({} anc vs {} win)",
        anc.len(),
        win.len()
    );
    assert_eq!(
        table.len(),
        1usize << win.len(),
        "qrom_read: table must have 2^w = {} entries, got {}",
        1usize << win.len(),
        table.len()
    );
    assert!(
        addend.len() <= 64,
        "qrom_read: u64 table entries support <= 64 addend bits, got {}",
        addend.len()
    );
    rec(circ, 0, None, 0, win, anc, addend, table);
}

fn toffoli_count<'a>(ops: impl Iterator<Item = &'a Op>) -> u64 {
    ops.filter(|o| matches!(o.kind, OperationType::CCX | OperationType::CCZ))
        .count() as u64
}

/// splitmix64 — deterministic table constants (no `rand`, no `Math::random`).
fn splitmix(mut z: u64) -> u64 {
    z = z.wrapping_add(0x9E37_79B9_7F4A_7C15);
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

/// Build read→add→unread for one (n, w) instance and check it by simulation over
/// all `2^w` windows × several accumulators. Returns (peak_qubits, per-read
/// selector Toffoli) for the width report.
fn run_instance(n: usize, w: usize) -> (u64, u64) {
    assert!((1..=6).contains(&w) && (2..=60).contains(&n));
    let table: Vec<u64> = (0..(1usize << w))
        .map(|a| splitmix(0x0ADD_0000 ^ (a as u64)) & ((1u64 << n) - 1))
        .collect();

    // Registers: window | addend | accumulator | selector spine | carry.
    let mut circ = B::new_for_test();
    let win = circ.alloc_qubits(w);
    let addend = circ.alloc_qubits(n);
    let acc = circ.alloc_qubits(n);
    let anc = circ.alloc_qubits(w);
    let carry = circ.alloc_qubits(1);

    qrom_read(&mut circ, &win, &anc, &addend, &table); // addend := P[k]
    cuccaro_add(&mut circ, &addend, &acc, carry[0]); // acc += P[k]  (addend preserved)
    qrom_read(&mut circ, &win, &anc, &addend, &table); // addend := 0  (unread)

    let ops = circ.take_ops();
    let (peak_qubits, nbits, _r, _regs) = analyze_ops(ops.iter());

    // One read's selector Toffoli must equal the unary-iteration cost (ADR 0010).
    let mut probe = B::new_for_test();
    let pw = probe.alloc_qubits(w);
    let pa = probe.alloc_qubits(n);
    let pn = probe.alloc_qubits(w);
    qrom_read(&mut probe, &pw, &pn, &pa, &table);
    let read_tof = toffoli_count(probe.take_ops().iter());
    assert_eq!(
        read_tof,
        (1u64 << (w + 1)) - 4,
        "unary-iteration read selector Toffoli != 2^(w+1)-4 (n={n}, w={w})"
    );

    // Masked multi-shot: lane s carries window k_s and accumulator y_s. Cover all
    // 2^w windows across the low bits, several accumulators across the high bits.
    let n_win = 1usize << w;
    let shots = 64usize;
    let mask_n = if n >= 64 { u64::MAX } else { (1u64 << n) - 1 };
    let k_of = |s: usize| s % n_win;
    let y_of = |s: usize| splitmix(0x0ACC_1100 ^ (s / n_win) as u64) & mask_n;

    let mut seed = sha3::Shake128::default();
    sha3::digest::Update::update(&mut seed, b"qaddend-testbed");
    let mut xof = sha3::digest::ExtendableOutput::finalize_xof(seed);
    let mut sim = Simulator::new(peak_qubits as usize, nbits as usize, &mut xof);
    sim.clear_for_shot();
    // Load window + accumulator masks (addend/anc/carry stay |0>).
    for (j, &q) in win.iter().enumerate() {
        let mut m = 0u64;
        for s in 0..shots {
            m |= (((k_of(s) >> j) & 1) as u64) << s;
        }
        *sim.qubit_mut(q) = m;
    }
    for (i, &q) in acc.iter().enumerate() {
        let mut m = 0u64;
        for s in 0..shots {
            m |= ((y_of(s) >> i) & 1) << s;
        }
        *sim.qubit_mut(q) = m;
    }
    sim.apply_iter(ops.iter());

    // Reconstruct per-shot results and check acc == (y + P[k]) mod 2^n.
    let read_reg = |sim: &Simulator<_>, reg: &[QubitId], s: usize| -> u64 {
        let mut v = 0u64;
        for (i, &q) in reg.iter().enumerate() {
            v |= ((sim.qubit(q) >> s) & 1) << i;
        }
        v
    };
    for s in 0..shots {
        let k = k_of(s);
        let y = y_of(s);
        let expect = y.wrapping_add(table[k]) & mask_n;
        let got = read_reg(&sim, &acc, s);
        assert_eq!(got, expect, "acc mismatch (n={n}, w={w}, shot={s}, k={k})");
        // window preserved; addend, spine, carry all clean.
        assert_eq!(
            read_reg(&sim, &win, s),
            k as u64,
            "window register perturbed"
        );
        assert_eq!(read_reg(&sim, &addend, s), 0, "addend not returned to |0>");
        assert_eq!(read_reg(&sim, &anc, s), 0, "selector ancilla dirty");
    }
    assert_eq!(sim.qubit(carry[0]), 0, "carry ancilla dirty");
    assert_eq!(sim.phase, 0, "unexpected phase (no phase gates emitted)");

    (peak_qubits, read_tof)
}

#[test]
fn qrom_fed_quantum_addend_add() {
    // A few (n, w) instances: exercises width-parametricity of the composition.
    let cases = [(8usize, 3usize), (6, 2), (10, 4), (16, 3)];
    eprintln!("\n=== issue #27/#28 quantum-addend testbed: QROM read -> q-q add -> unread ===");
    eprintln!("  (functional correctness by simulation; ADR 0014)");
    for (n, w) in cases {
        let (peak, read_tof) = run_instance(n, w);
        // Register breakdown: window(w) + addend(n) + acc(n) + spine(w) + carry(1).
        let expect_peak = (w + n + n + w + 1) as u64;
        assert_eq!(peak, expect_peak, "unexpected peak width for n={n}, w={w}");
        eprintln!(
            "  n={n:<2} w={w}: PASS  peak={peak} qubits [win {w} + addend {n} + acc {n} + spine {w} + carry 1], \
             read selector Toffoli={read_tof} (=2^(w+1)-4)"
        );
        eprintln!(
            "        -> QROM overhead over the bare adder (acc {n} + carry 1): addend {n} + spine {w} = {} qubits",
            n + w
        );
    }
    eprintln!("  => the QROM-fed quantum addend computes acc+P[k] correctly with all ancilla");
    eprintln!("     clean; its addend+spine ride ON TOP of the adder (small-scale ADR 0013).");
}
