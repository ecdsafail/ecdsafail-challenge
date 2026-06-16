pub(crate) mod primitives;
pub(crate) mod rfold;
pub(crate) mod schedule;
pub(crate) mod state_machine;

use crate::sim::Simulator;

use super::B;

pub(crate) fn schedule_core_peak() -> (usize, usize) {
    let mut best_step = 0;
    let mut best_peak = 0;
    for step in 0..schedule::SHRUNKEN_PZ_NSTEPS {
        let (a, b, ca, cb, q) = schedule::reg_widths(step);
        let peak = a + b + ca + cb + q;
        if peak > best_peak {
            best_peak = peak;
            best_step = step;
        }
    }
    (best_step, best_peak)
}

pub(crate) fn schedule_selftest() -> Result<(), String> {
    let (step, peak) = schedule_core_peak();
    if peak != 741 {
        return Err(format!(
            "unexpected shrunken-PZ core peak {peak} at step {step}; expected 741"
        ));
    }
    let (a, b, ca, cb, q) = schedule::reg_widths(step);
    if a + b + ca + cb + q != peak {
        return Err(format!("core peak components do not sum at step {step}"));
    }
    Ok(())
}

pub(crate) fn primitives_selftest() -> Result<(), String> {
    use sha3::{
        digest::{ExtendableOutput, Update},
        Shake256,
    };

    const A_BITS: usize = 6;
    const B_BITS: usize = 3;
    const Q_BITS: usize = A_BITS - B_BITS + 1;

    let mut b = B::new();
    let a = b.alloc_qubits(A_BITS);
    let d = b.alloc_qubits(B_BITS);
    let q = b.alloc_qubits(Q_BITS);
    primitives::long_division(&mut b, &a, &d, &q);
    primitives::long_division_reverse(&mut b, &a, &d, &q);

    let mut seed = Shake256::default();
    seed.update(b"shrunken-pz-primitives");
    let mut xof = seed.finalize_xof();
    let mut sim = Simulator::new(b.next_qubit as usize, b.next_bit as usize, &mut xof);
    let mut expected_a = [0u64; 64];
    let mut expected_d = [0u64; 64];
    for shot in 0..64usize {
        let d_val = 1 + ((shot as u64 * 5 + 1) % ((1u64 << B_BITS) - 1));
        let a_val = (shot as u64 * 11 + 7) % (1u64 << A_BITS);
        expected_a[shot] = a_val;
        expected_d[shot] = d_val;
        set_reg(&mut sim, &a, a_val, shot);
        set_reg(&mut sim, &d, d_val, shot);
    }
    sim.apply_iter(b.ops.iter());
    if sim.phase != 0 {
        return Err(format!(
            "phase garbage after long-division roundtrip: 0x{:016x}",
            sim.phase
        ));
    }
    for shot in 0..64usize {
        let got_a = get_reg(&sim, &a, shot);
        let got_d = get_reg(&sim, &d, shot);
        let got_q = get_reg(&sim, &q, shot);
        if got_a != expected_a[shot] || got_d != expected_d[shot] || got_q != 0 {
            return Err(format!(
                "shot {shot}: got a={got_a}, d={got_d}, q={got_q}; expected a={}, d={}, q=0",
                expected_a[shot], expected_d[shot]
            ));
        }
    }
    Ok(())
}

fn set_reg<R: sha3::digest::XofReader>(
    sim: &mut Simulator<'_, R>,
    qs: &[crate::circuit::QubitId],
    value: u64,
    shot: usize,
) {
    for (i, &q) in qs.iter().enumerate() {
        if (value >> i) & 1 != 0 {
            *sim.qubit_mut(q) |= 1u64 << shot;
        }
    }
}

fn get_reg<R: sha3::digest::XofReader>(
    sim: &Simulator<'_, R>,
    qs: &[crate::circuit::QubitId],
    shot: usize,
) -> u64 {
    let mut out = 0u64;
    for (i, &q) in qs.iter().enumerate() {
        out |= ((sim.qubit(q) >> shot) & 1) << i;
    }
    out
}
