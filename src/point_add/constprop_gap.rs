//! Tier B (issue #27), first increment: MEASURE the classical-vs-quantum-addend
//! gap — the one assumption behind `ecdlp_estimate.py`'s ladder and behind
//! `ladder_full.rs`'s `28·PA` (ADR 0011). Finding: the gap is **negligible**.
//!
//! `ecdlp_estimate.py` warns that the scored point-add folds a *classical,
//! compile-time* addend, so its cost "may not fully transfer to the [quantum-
//! addend] windowed setting." This harness shows that concern is essentially
//! moot: the scored PA does **not** exploit the constant addend to cheapen its
//! arithmetic.
//!
//! Structural reason (`coord_addsub`, `trailmix_ludicrous/ec_add.rs`): to consume
//! the classical addend `coord: &[BitId]`, the circuit allocates a fresh **qubit**
//! register, loads the addend into it with `x_if_bit` (Clifford `X`, not Toffoli),
//! runs an **uncontrolled quantum-quantum** vented mod-add/sub (a full Cuccaro-
//! class adder — the same Toffoli cost a *quantum* addend would pay), then unloads.
//! So the Toffoli count is addend-*value*-independent; the classical value only
//! drives the (Clifford) load/unload.
//!
//! Measured confirmation (build twice, in-process, count static CCX/CCZ):
//!   - the only addend-value-dependent optimization is the peephole
//!     constant-propagation pass (`CONSTPROP_DISABLE`), worth ~770 Toffoli
//!     (~0.05% of PA);
//!   - the direct-constant-arithmetic knobs (`SECP_DIRECT_CONST_ARITH`,
//!     `KAL_DIRECT_CONST_WALKS`) are **inert** for the trailmix build path
//!     (identical Toffoli), i.e. the scored circuit already takes the load-into-
//!     qubits path.
//!
//! => the classical-vs-quantum-addend gap is ≤ ~0.05% of PA. The scored PA — and
//! therefore `ladder_full.rs`'s `28·PA` and `ecdlp_estimate.py`'s headline —
//! already reflect the quantum-addend *arithmetic* cost. What #27 still needs is
//! the register-overlap / width question (A2: does the QROM-provided addend
//! register reuse PA's workspace) and functional correctness of a QROM-fed add;
//! the Toffoli gap, the part `ecdlp_estimate.py` flagged, is resolved here.
//!
//! `#[cfg(test)]` only; never compiled into the scored circuit. This test toggles
//! process-global env vars, so run it ALONE, e.g.
//!   `cargo test constprop_addend_gap -- --ignored --exact \
//!      point_add::constprop_gap::constprop_addend_gap`
//! (running it in parallel with the other build()-calling #[ignore] tests could
//! perturb their measurements).

use super::build;
use crate::circuit::{analyze_ops, OperationType};

fn toffoli_count<'a>(ops: impl Iterator<Item = &'a crate::circuit::Op>) -> u64 {
    ops.filter(|o| matches!(o.kind, OperationType::CCX | OperationType::CCZ))
        .count() as u64
}

fn build_and_measure() -> (u64, u64) {
    let ops = build();
    let tof = toffoli_count(ops.iter());
    let (qubits, _b, _r, _regs) = analyze_ops(ops.iter());
    (tof, qubits)
}

/// Sets an env var for its lifetime and restores the *prior* value on drop, so a
/// panic mid-test cannot leak a mutated process environment into other tests
/// (`build()` reads these vars to pick circuit-generation paths).
struct EnvGuard {
    key: &'static str,
    prior: Option<std::ffi::OsString>,
}

impl EnvGuard {
    fn set(key: &'static str, val: &str) -> Self {
        let prior = std::env::var_os(key);
        std::env::set_var(key, val);
        Self { key, prior }
    }
}

impl Drop for EnvGuard {
    fn drop(&mut self) {
        match &self.prior {
            Some(v) => std::env::set_var(self.key, v),
            None => std::env::remove_var(self.key),
        }
    }
}

#[test]
#[ignore = "heavy; toggles env vars, run ALONE with --ignored --exact"]
fn constprop_addend_gap() {
    // These toggles are process-global; serialize so a manual multi-test run can
    // never observe another test's env (poison-tolerant: we only need exclusion).
    static SERIAL: std::sync::Mutex<()> = std::sync::Mutex::new(());
    let _serial = SERIAL.lock().unwrap_or_else(|e| e.into_inner());

    // Direct-constant-arithmetic ON: would exploit the classical constant in the
    // adders instead of loading it into a qubit register. Guards restore the prior
    // env when the scope ends (or on unwind), so each phase builds under exactly
    // the vars it sets and nothing leaks afterward.
    let (dc_tof, _dc_q) = {
        let _g1 = EnvGuard::set("SECP_DIRECT_CONST_ARITH", "1");
        let _g2 = EnvGuard::set("KAL_DIRECT_CONST_WALKS", "1");
        build_and_measure()
    };

    // Peephole constprop OFF, then restored so the default (on) build below is
    // measured with the pre-existing environment.
    let (off_tof, _off_q) = {
        let _g = EnvGuard::set("CONSTPROP_DISABLE", "1");
        build_and_measure()
    };
    let (on_tof, on_q) = build_and_measure();

    assert!(on_tof > 0, "degenerate build");
    // Constprop can only remove / downgrade Toffoli, never add them.
    assert!(
        off_tof >= on_tof,
        "constprop-off Toffoli ({off_tof}) < constprop-on ({on_tof}) — unexpected"
    );

    let constprop_saving = off_tof - on_tof;
    let pct = 100.0 * constprop_saving as f64 / on_tof as f64;

    eprintln!("\n=== issue #27 classical-vs-quantum-addend gap ===");
    eprintln!("  scored PA (default)              : toffoli={on_tof}  qubits={on_q}");
    eprintln!("  peephole constprop OFF           : toffoli={off_tof}  (+{constprop_saving} = the");
    eprintln!("      only addend-value-dependent Toffoli, {pct:.2}% of PA)");
    eprintln!("  direct-const-arith ON            : toffoli={dc_tof}  (inert: trailmix already");
    eprintln!("      loads the addend into qubits + uncontrolled q-q Cuccaro add)");
    eprintln!("  => classical-vs-quantum-addend Toffoli gap <= {pct:.2}% of PA — negligible.");
    eprintln!("     ladder_full.rs's 28·PA already reflects the quantum-addend arithmetic.");

    // Locked findings:
    // (1) the classical-addend advantage in Toffoli is < 1% of PA (the arithmetic
    //     is addend-value-independent: loaded into qubits, not const-exploited).
    assert!(
        constprop_saving * 100 < on_tof, // saving < 1% of PA
        "classical-addend Toffoli advantage {constprop_saving} exceeds 1% of PA {on_tof}"
    );
    // (2) the direct-const-arith knobs do not change the trailmix Toffoli, i.e.
    //     the scored path is already the load-into-qubits (quantum-equivalent) one.
    assert_eq!(
        dc_tof, on_tof,
        "direct-const-arith changed the Toffoli — the gap analysis needs revisiting"
    );
}
