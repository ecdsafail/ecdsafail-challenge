//! Per-gate **hotness** census: for every CCX/CCZ in the scored stream, the number
//! of the 9,024 official shots on which its classical condition stack is satisfied.
//!
//! ## Why this is the quantity that matters
//!
//! The scorer charges a Toffoli by the popcount of its condition mask, and by
//! nothing else (`src/sim.rs:77-86`, verbatim):
//!
//! ```text
//! let mut cond = current_base_condition;                       // PushCondition stack
//! if op.c_condition != NO_BIT { cond &= self.bit(op.c_condition); }
//! let executed_shots = cond.count_ones() as u64;
//! match op.kind {
//!     OperationType::CCZ | OperationType::CCX => {
//!         self.stats.toffoli_gates += executed_shots;
//!     }
//! ```
//!
//! The charge does **not** depend on the control qubits. So a gate is charged in
//! full even when its controls make it a no-op, and a gate whose condition stack
//! is false on a shot is free on that shot. `score = round(total/9024) x Q`.
//!
//! That splits the cost surface into two independent levers, only one of which
//! the deep strip addresses:
//!
//!   - **delete** a gate: needs a proof it never *fires* (effect mask always zero).
//!     That is what a fire-census certifies, and it is expensive to establish and
//!     fragile under stream edits. See `../07-syntactic-certification.md` for why
//!     the sampled version of that proof cannot be made stream-agnostic.
//!   - **cool** a gate: tighten its condition stack so it is charged on fewer
//!     shots. Semantics are preserved for free whenever the added condition is
//!     implied by the existing one on every shot where the gate can fire.
//!
//! A gate at hotness h costs exactly h Toffoli-shots. Cooling a gate from h to h'
//! saves (h - h')/9024 of average Toffoli, whether or not it ever fires.
//!
//! ## Method
//!
//! This replays the **official** measurement: `point_add::build()`, `analyze_ops`
//! for the register layout, the Fiat-Shamir XOF over the whole op stream, the same
//! 9,024-shot draw, the same 64-shot batching, and the condition/dispatch loop
//! copied verbatim from `src/sim.rs::apply_iter` with one line added to attribute
//! `executed_shots` to the op that incurred it.
//!
//! **Gate on the instrument.** Nothing here is trustworthy unless the attributed
//! charges sum to the scorer's own total. The tool asserts
//! `sum(charge) == sim.stats.toffoli_gates` and prints avgT so it can be checked
//! against `eval_circuit` by eye. Run it on a circuit you have just scored.
//!
//! ## Usage
//!
//! There is no `[[example]]` entry in `Cargo.toml`, so this is built the same way
//! `hyperplane_mitm.cpp` is: copied out, compiled, and the scratch directory removed.
//! Run from the repository root.
//!
//!     mkdir -p examples && cp src/point_add/memory/repro/hotness.rs examples/
//!     cargo build --release --offline --example hotness
//!     ./target/release/examples/hotness out        # writes out.hot.tsv
//!     ./target/release/examples/hotness -          # summary only, no files
//!     rm -rf examples
//!
//! Honours `SUB4_APPLY_STRIP` / `TLM_SCHED_J2_DELTA` like every other instrument;
//! run it with the same environment as the build you are pricing. Default (no env)
//! is the scored circuit.
//!
//! Output `<prefix>.hot.tsv`, one line per CCX/CCZ in stream order:
//!     opidx  kind  c2  c1  t  cond  charge  fire  head  hotness
//! where `charge` is shots-charged out of `n_shots`, `fire` is shots on which the
//! effect mask is non-zero, `head = charge - fire` is the per-gate cooling ceiling,
//! and `hotness = charge/n_shots`.

use quantum_ecc::circuit::{analyze_ops, BitId, Op, OperationType, QubitId, QubitOrBit, NO_BIT};
use quantum_ecc::point_add;
use quantum_ecc::sim::Simulator;
use quantum_ecc::weierstrass_elliptic_curve::WeierstrassEllipticCurve;
use ruint::aliases::U256;
use sha3::digest::{ExtendableOutput, Update, XofReader};
use sha3::Shake256;
use std::collections::HashMap;
use std::io::Write;

const NUM_TESTS: usize = 9024;
const BATCH: usize = 64;

// Verbatim from src/bin/eval_circuit.rs.
fn secp256k1() -> WeierstrassEllipticCurve {
    WeierstrassEllipticCurve {
        a: U256::from(0u64),
        b: U256::from(7u64),
        modulus: U256::from_str_radix(
            "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F",
            16,
        )
        .unwrap(),
        gx: U256::from_str_radix(
            "79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798",
            16,
        )
        .unwrap(),
        gy: U256::from_str_radix(
            "483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8",
            16,
        )
        .unwrap(),
        order: U256::from_str_radix(
            "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141",
            16,
        )
        .unwrap(),
    }
}

// Verbatim from src/bin/eval_circuit.rs.
fn fiat_shamir_seed(ops: &[Op]) -> sha3::Shake256Reader {
    let mut hasher = Shake256::default();
    hasher.update(b"quantum_ecc-fiat-shamir-v2");
    hasher.update(&(ops.len() as u64).to_le_bytes());
    for op in ops {
        hasher.update(&[op.kind as u8]);
        hasher.update(&op.q_control2.0.to_le_bytes());
        hasher.update(&op.q_control1.0.to_le_bytes());
        hasher.update(&op.q_target.0.to_le_bytes());
        hasher.update(&op.c_target.0.to_le_bytes());
        hasher.update(&op.c_condition.0.to_le_bytes());
        hasher.update(&op.r_target.0.to_le_bytes());
    }
    hasher.finalize_xof()
}

/// `src/sim.rs::apply_iter`, copied dispatch-for-dispatch, plus per-op attribution
/// of `executed_shots`. Any divergence from the trusted loop shows up as a failed
/// total-charge assertion, which is why that assertion is not optional.
fn apply_iter_charged<R: XofReader>(
    sim: &mut Simulator<'_, R>,
    ops: &[Op],
    charge: &mut [u64],
    fire: &mut [u64],
) {
    let mut condition_stack: Vec<u64> = Vec::new();
    let mut current_base_condition = u64::MAX;

    for (i, op) in ops.iter().enumerate() {
        let mut cond = current_base_condition;
        if op.c_condition != NO_BIT {
            cond &= sim.bit(op.c_condition);
        }

        let executed_shots = cond.count_ones() as u64;

        match op.kind {
            OperationType::CCZ | OperationType::CCX => {
                sim.stats.toffoli_gates += executed_shots;
                charge[i] += executed_shots;
            }
            OperationType::CX
            | OperationType::CZ
            | OperationType::Swap
            | OperationType::R
            | OperationType::Hmr => {
                sim.stats.clifford_gates += executed_shots;
            }
            _ => {}
        }

        match op.kind {
            OperationType::CCX => {
                let v = cond & sim.qubit(op.q_control1) & sim.qubit(op.q_control2);
                fire[i] += v.count_ones() as u64;
                *sim.qubit_mut(op.q_target) ^= v;
            }
            OperationType::CX => {
                let v = cond & sim.qubit(op.q_control1);
                *sim.qubit_mut(op.q_target) ^= v;
            }
            OperationType::Swap => {
                let mut q_c1 = sim.qubit(op.q_control1);
                let mut q_t = sim.qubit(op.q_target);
                q_c1 ^= q_t;
                q_t ^= cond & q_c1;
                q_c1 ^= q_t;
                *sim.qubit_mut(op.q_control1) = q_c1;
                *sim.qubit_mut(op.q_target) = q_t;
            }
            OperationType::X => {
                *sim.qubit_mut(op.q_target) ^= cond;
            }
            OperationType::CCZ => {
                let v = cond
                    & sim.qubit(op.q_target)
                    & sim.qubit(op.q_control1)
                    & sim.qubit(op.q_control2);
                fire[i] += v.count_ones() as u64;
                sim.phase ^= v;
            }
            OperationType::CZ => {
                let v = cond & sim.qubit(op.q_target) & sim.qubit(op.q_control1);
                sim.phase ^= v;
            }
            OperationType::Z => {
                let v = cond & sim.qubit(op.q_target);
                sim.phase ^= v;
            }
            OperationType::Neg => {
                sim.phase ^= cond;
            }
            OperationType::Hmr => {
                let mut buf = [0u8; 8];
                sim.xof.read(&mut buf);
                let rng_val = u64::from_le_bytes(buf);
                *sim.bit_mut(op.c_target) &= !cond;
                *sim.bit_mut(op.c_target) ^= rng_val & cond;
                sim.phase ^= sim.qubit(op.q_target) & rng_val & cond;
                *sim.qubit_mut(op.q_target) &= !cond;
            }
            OperationType::R => {
                let mut buf = [0u8; 8];
                sim.xof.read(&mut buf);
                let rng_val = u64::from_le_bytes(buf);
                sim.phase ^= sim.qubit(op.q_target) & rng_val & cond;
                *sim.qubit_mut(op.q_target) &= !cond;
            }
            OperationType::BitInvert => {
                *sim.bit_mut(op.c_target) ^= cond;
            }
            OperationType::BitStore0 => {
                *sim.bit_mut(op.c_target) &= !cond;
            }
            OperationType::BitStore1 => {
                *sim.bit_mut(op.c_target) |= cond;
            }
            OperationType::AppendToRegister
            | OperationType::Register
            | OperationType::DebugPrint => {}
            OperationType::PushCondition => {
                condition_stack.push(current_base_condition);
                current_base_condition &= sim.bit(op.c_condition);
            }
            OperationType::PopCondition => {
                if let Some(val) = condition_stack.pop() {
                    current_base_condition = val;
                }
            }
        }
    }
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let prefix = args.get(1).cloned().unwrap_or_else(|| "hot".to_string());
    let write_files = prefix != "-";

    let ops: Vec<Op> = point_add::build();
    let (total_qubits, num_bits, _nregs, regs) = analyze_ops(ops.iter());
    println!("ops={} qubits={} bits={}", ops.len(), total_qubits, num_bits);

    let curve = secp256k1();
    let mut xof = fiat_shamir_seed(&ops);

    // Same draw and same rejection rules as eval_circuit::run_tests.
    let mut targets = Vec::with_capacity(NUM_TESTS);
    let mut offsets = Vec::with_capacity(NUM_TESTS);
    for _ in 0..NUM_TESTS {
        let mut rb = [[0u8; 32]; 2];
        XofReader::read(&mut xof, &mut rb[0]);
        XofReader::read(&mut xof, &mut rb[1]);
        let k1 = U256::from_le_bytes(rb[0]);
        let k2 = U256::from_le_bytes(rb[1]);
        let t = curve.mul(curve.gx, curve.gy, k1);
        let o = curve.mul(curve.gx, curve.gy, k2);
        if t.0 == o.0 {
            continue;
        }
        if t.0.is_zero() && t.1.is_zero() {
            continue;
        }
        if o.0.is_zero() && o.1.is_zero() {
            continue;
        }
        targets.push(t);
        offsets.push(o);
    }
    let n = targets.len();
    println!("shots={n}");

    let mut charge = vec![0u64; ops.len()];
    let mut fire = vec![0u64; ops.len()];
    let mut sim = Simulator::new(total_qubits as usize, num_bits as usize, &mut xof);

    assert_eq!(n % BATCH, 0, "partial final batch would need lane masking for fire counts");
    let num_batches = (n + BATCH - 1) / BATCH;
    for batch in 0..num_batches {
        let bs = BATCH.min(n - batch * BATCH);
        sim.clear_for_shot();
        for shot in 0..bs {
            let i = batch * BATCH + shot;
            sim.set_register(&regs[0], targets[i].0, shot);
            sim.set_register(&regs[1], targets[i].1, shot);
            sim.set_register(&regs[2], offsets[i].0, shot);
            sim.set_register(&regs[3], offsets[i].1, shot);
        }
        apply_iter_charged(&mut sim, &ops, &mut charge, &mut fire);
    }

    // ---- the gate: attribution must reconstruct the scorer's own total ----
    let attributed: u64 = ops
        .iter()
        .enumerate()
        .filter(|(_, op)| matches!(op.kind, OperationType::CCX | OperationType::CCZ))
        .map(|(i, _)| charge[i])
        .sum();
    let total = sim.stats.toffoli_gates;
    assert_eq!(
        attributed, total,
        "attribution does not reconstruct the scorer total"
    );
    let avg_t = total as f64 / n as f64;
    let gates_n = ops
        .iter()
        .filter(|op| matches!(op.kind, OperationType::CCX | OperationType::CCZ))
        .count();
    println!("GATE ok: attributed={attributed} == toffoli_gates={total}");
    println!("avgT={avg_t:.3}  (compare eval_circuit's 'avg executed Toffoli')");

    // ---- the cooling bound ----
    //
    // A semantics-preserving condition must be TRUE on every shot where the gate
    // fires, so the most any condition can save on gate g is
    //     charge(g) - fire(g)
    // shot-charges. Summed over the stream this is an absolute ceiling on the
    // "make gates conditionally cold" lever, and it needs no candidate bit to be
    // found -- it grants a perfect oracle condition to every gate at once.
    let mut head_total: u64 = 0;
    let mut fired_total: u64 = 0;
    let mut dead_gates = 0usize;
    let mut dead_charge = 0u64;
    let mut always_fire = 0usize;
    for (i, op) in ops.iter().enumerate() {
        if !matches!(op.kind, OperationType::CCX | OperationType::CCZ) {
            continue;
        }
        head_total += charge[i] - fire[i];
        fired_total += fire[i];
        if fire[i] == 0 {
            dead_gates += 1;
            dead_charge += charge[i];
        }
        if fire[i] == charge[i] {
            always_fire += 1;
        }
    }
    println!("\n=== COOLING BOUND (perfect-oracle condition on every gate) ===");
    println!("total charge         : {total}");
    println!("total fired          : {fired_total}  ({:.4}% of charge)", 100.0*fired_total as f64/total as f64);
    println!("COOLING HEADROOM     : {head_total} shot-charges = {:.3} avgT ({:.4}% of score)",
             head_total as f64 / n as f64, 100.0*head_total as f64/total as f64);
    println!("gates that NEVER fire: {dead_gates} carrying {dead_charge} charge = {:.3} avgT",
             dead_charge as f64 / n as f64);
    println!("gates fired on every charged shot (uncoolable): {always_fire} ({:.3}%)",
             100.0*always_fire as f64 / gates_n as f64);

    // ---- distribution ----
    let gates: Vec<(usize, &Op)> = ops
        .iter()
        .enumerate()
        .filter(|(_, op)| matches!(op.kind, OperationType::CCX | OperationType::CCZ))
        .collect();
    println!("gates={}", gates.len());

    let mut cold = 0usize;
    let mut full = 0usize;
    let mut partial = 0usize;
    let mut charge_from_partial = 0u64;
    for (i, _) in &gates {
        let c = charge[*i];
        if c == 0 {
            cold += 1;
        } else if c == n as u64 {
            full += 1;
        } else {
            partial += 1;
            charge_from_partial += c;
        }
    }
    println!(
        "COLD (charge=0)      : {cold} gates ({:.3}% of gates), 0 charge",
        100.0 * cold as f64 / gates.len() as f64
    );
    println!(
        "FULL (charge={n})   : {full} gates ({:.3}%), {} charge ({:.3}% of total)",
        100.0 * full as f64 / gates.len() as f64,
        full as u64 * n as u64,
        100.0 * (full as u64 * n as u64) as f64 / total as f64
    );
    println!(
        "PARTIAL              : {partial} gates ({:.3}%), {charge_from_partial} charge ({:.3}% of total)",
        100.0 * partial as f64 / gates.len() as f64,
        100.0 * charge_from_partial as f64 / total as f64
    );

    // Decile histogram over hotness, with the charge each decile carries.
    let mut hist = [0usize; 11];
    let mut hist_charge = [0u64; 11];
    for (i, _) in &gates {
        let h = charge[*i] as f64 / n as f64;
        let b = ((h * 10.0).floor() as usize).min(10);
        hist[b] += 1;
        hist_charge[b] += charge[*i];
    }
    println!("\nHOTNESS DECILES (bucket = floor(hotness*10), 10 = exactly 1.0)");
    println!("bucket\tgates\tshare_gates\tcharge\tshare_charge");
    for b in 0..=10 {
        println!(
            "{}\t{}\t{:.4}%\t{}\t{:.4}%",
            b,
            hist[b],
            100.0 * hist[b] as f64 / gates.len() as f64,
            hist_charge[b],
            100.0 * hist_charge[b] as f64 / total as f64
        );
    }

    // Charge grouped by the op's own condition bit (NO_BIT = conditioned only by
    // the enclosing PushCondition stack, if any).
    let mut by_cond: HashMap<u64, (usize, u64)> = HashMap::new();
    for (i, op) in &gates {
        let e = by_cond.entry(op.c_condition.0).or_insert((0, 0));
        e.0 += 1;
        e.1 += charge[*i];
    }
    let mut bc: Vec<(u64, (usize, u64))> = by_cond.into_iter().collect();
    bc.sort_by_key(|(_, (_, c))| std::cmp::Reverse(*c));
    println!("\nTOP CONDITION BITS BY CHARGE (cond bit, gates, charge, share, mean hotness)");
    for (cb, (ng, c)) in bc.iter().take(25) {
        let label = if *cb == u64::MAX {
            "NO_BIT".to_string()
        } else {
            cb.to_string()
        };
        println!(
            "{}\t{}\t{}\t{:.4}%\t{:.4}",
            label,
            ng,
            c,
            100.0 * *c as f64 / total as f64,
            *c as f64 / (*ng as f64 * n as f64)
        );
    }
    println!("distinct_cond_bits={}", bc.len());

    // Charge grouped by operand tuple -- the unit the strip keys on, so this is
    // directly comparable to the census tables.
    let mut by_tuple: HashMap<(u8, u64, u64, u64, u64), (usize, u64)> = HashMap::new();
    for (i, op) in &gates {
        let k = (
            op.kind as u8,
            op.q_control2.0,
            op.q_control1.0,
            op.q_target.0,
            op.c_condition.0,
        );
        let e = by_tuple.entry(k).or_insert((0, 0));
        e.0 += 1;
        e.1 += charge[*i];
    }
    let mut bt: Vec<_> = by_tuple.into_iter().collect();
    bt.sort_by_key(|(_, (_, c))| std::cmp::Reverse(*c));
    println!("\nTOP OPERAND TUPLES BY CHARGE (kind,c2,c1,t,cond, occurrences, charge, share)");
    for (k, (ng, c)) in bt.iter().take(25) {
        println!(
            "{},{},{},{},{}\t{}\t{}\t{:.4}%",
            k.0,
            k.1,
            k.2,
            k.3,
            if k.4 == u64::MAX {
                "NO_BIT".to_string()
            } else {
                k.4.to_string()
            },
            ng,
            c,
            100.0 * *c as f64 / total as f64
        );
    }
    println!("distinct_tuples={}", bt.len());

    if write_files {
        let f = std::fs::File::create(format!("{prefix}.hot.tsv")).unwrap();
        let mut w = std::io::BufWriter::new(f);
        writeln!(w, "opidx\tkind\tc2\tc1\tt\tcond\tcharge\tfire\thead\thotness").unwrap();
        for (i, op) in &gates {
            writeln!(
                w,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.6}",
                i,
                op.kind as u8,
                op.q_control2.0,
                op.q_control1.0,
                op.q_target.0,
                op.c_condition.0,
                charge[*i],
                fire[*i],
                charge[*i] - fire[*i],
                charge[*i] as f64 / n as f64
            )
            .unwrap();
        }
        println!("\nwrote {prefix}.hot.tsv");
    }

    let _ = (BitId(0), QubitId(0), QubitOrBit::Qubit(QubitId(0)));
}
