//! Proof harness for omitting the physical low bit of `l_s`.
//!
//! At every completed scheduled-step boundary, `l_s mod 2 = step mod 2`.
//! This module proves the gate-level pre/post-shift consequences before the
//! representation is admitted into the production divider.

use super::register_shared_eea_microkernels::{
    post_shift, post_shift_inverse, pre_shift, pre_shift_inverse,
};
use crate::circuit::{OperationType, QubitId};
use crate::point_add::trailmix_port::circuit::{Circuit, QReg};
use crate::point_add::B;
use crate::sim::Simulator;
use sha3::{
    digest::{ExtendableOutput, Update},
    Shake128,
};

const PRODUCTION_HIGH_WIDTH: usize = 8;
const PRODUCTION_WORK_WIDTH: usize = 259;
const REDUCED_WORK_WIDTH: usize = 4;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum ProofMode {
    Pre,
    Post,
    RoundTrip,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
struct InputState {
    phase1: bool,
    phase2: bool,
    high: usize,
    work: u64,
    work_seed: Option<u64>,
}

struct Harness {
    builder: B,
    phase1: u32,
    phase2: u32,
    work: Vec<u32>,
    high: Vec<u32>,
    low_or_host: u32,
    external: Vec<bool>,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct LsParityProofReport {
    pub high_widths_checked: Vec<usize>,
    pub modes_checked: usize,
    pub step_parities_checked: usize,
    pub reduced_basis_states_checked: usize,
    pub production_sample_states_checked: usize,
    pub logical_register_checks: usize,
    pub host_restoration_checks: usize,
    pub phase_checks: usize,
    pub ancilla_cleanup_checks: usize,
    pub production_baseline_ops: [usize; 3],
    pub production_candidate_ops: [usize; 3],
    pub production_baseline_toffoli: [usize; 3],
    pub production_candidate_toffoli: [usize; 3],
}

fn free_clean(circ: &mut Circuit, registers: Vec<QReg>) {
    for register in registers {
        circ.zero_and_free(register);
    }
}

fn materialize_constant(circ: &mut Circuit, host: &QReg, value: bool) {
    if value {
        circ.x(host);
    }
}

/// Toggle `host` by `p xor 1 xor phase1`, where `p` is the entry parity.
fn toggle_mid_parity(circ: &mut Circuit, host: &QReg, phase1: &QReg, p: bool) {
    if !p {
        circ.x(host);
    }
    circ.cx(phase1, host);
}

fn full_view<'a>(host: &'a QReg, high: &'a [QReg]) -> Vec<&'a QReg> {
    std::iter::once(host).chain(high).collect()
}

fn run_candidate_mode(
    circ: &mut Circuit,
    mode: ProofMode,
    p: bool,
    phase1: &QReg,
    phase2: &QReg,
    work: &[QReg],
    high: &[QReg],
    host: &QReg,
) {
    match mode {
        ProofMode::Pre => {
            materialize_constant(circ, host, p);
            let view = full_view(host, high);
            let scratch = circ.alloc_qreg_bits("ls-parity.pre.scratch", view.len() + 4);
            let owned_view = view
                .iter()
                .map(|lane| lane.borrowed_alias())
                .collect::<Vec<_>>();
            pre_shift(circ, phase1, phase2, work, &owned_view, &scratch);
            free_clean(circ, scratch);
            toggle_mid_parity(circ, host, phase1, p);
        }
        ProofMode::Post => {
            toggle_mid_parity(circ, host, phase1, p);
            let view = full_view(host, high);
            let scratch = circ.alloc_qreg_bits("ls-parity.post.scratch", view.len() + 4);
            let owned_view = view
                .iter()
                .map(|lane| lane.borrowed_alias())
                .collect::<Vec<_>>();
            post_shift(circ, phase1, phase2, work, &owned_view, &scratch);
            free_clean(circ, scratch);
            materialize_constant(circ, host, !p);
        }
        ProofMode::RoundTrip => {
            materialize_constant(circ, host, p);
            let view = full_view(host, high);
            let owned_view = view
                .iter()
                .map(|lane| lane.borrowed_alias())
                .collect::<Vec<_>>();
            let pre_scratch =
                circ.alloc_qreg_bits("ls-parity.roundtrip.pre.scratch", view.len() + 4);
            pre_shift(circ, phase1, phase2, work, &owned_view, &pre_scratch);
            free_clean(circ, pre_scratch);
            toggle_mid_parity(circ, host, phase1, p);

            toggle_mid_parity(circ, host, phase1, p);
            let post_scratch =
                circ.alloc_qreg_bits("ls-parity.roundtrip.post.scratch", view.len() + 4);
            post_shift(circ, phase1, phase2, work, &owned_view, &post_scratch);
            free_clean(circ, post_scratch);
            materialize_constant(circ, host, !p);

            materialize_constant(circ, host, !p);
            let post_inverse_scratch =
                circ.alloc_qreg_bits("ls-parity.roundtrip.post-inverse.scratch", view.len() + 4);
            post_shift_inverse(
                circ,
                phase1,
                phase2,
                work,
                &owned_view,
                &post_inverse_scratch,
            );
            free_clean(circ, post_inverse_scratch);
            toggle_mid_parity(circ, host, phase1, p);

            toggle_mid_parity(circ, host, phase1, p);
            let pre_inverse_scratch =
                circ.alloc_qreg_bits("ls-parity.roundtrip.pre-inverse.scratch", view.len() + 4);
            pre_shift_inverse(
                circ,
                phase1,
                phase2,
                work,
                &owned_view,
                &pre_inverse_scratch,
            );
            free_clean(circ, pre_inverse_scratch);
            materialize_constant(circ, host, p);
        }
    }
}

fn build_harness(
    candidate: bool,
    mode: ProofMode,
    p: bool,
    high_width: usize,
    work_width: usize,
) -> Harness {
    assert!(high_width > 0);
    assert!(work_width >= 3);
    let mut circ = Circuit::new();
    let phase1 = circ.alloc_qreg("ls-parity.phase1");
    let phase2 = circ.alloc_qreg("ls-parity.phase2");
    let work = circ.alloc_qreg_bits("ls-parity.work", work_width);
    let high = circ.alloc_qreg_bits("ls-parity.high", high_width);
    let low_or_host = circ.alloc_qreg("ls-parity.low-or-host");

    if candidate {
        run_candidate_mode(
            &mut circ,
            mode,
            p,
            &phase1,
            &phase2,
            &work,
            &high,
            &low_or_host,
        );
    } else {
        let full = std::iter::once(low_or_host.borrowed_alias())
            .chain(high.iter().map(QReg::borrowed_alias))
            .collect::<Vec<_>>();
        match mode {
            ProofMode::Pre => {
                let scratch = circ.alloc_qreg_bits("ls-parity.baseline.pre", full.len() + 4);
                pre_shift(&mut circ, &phase1, &phase2, &work, &full, &scratch);
                free_clean(&mut circ, scratch);
            }
            ProofMode::Post => {
                let scratch = circ.alloc_qreg_bits("ls-parity.baseline.post", full.len() + 4);
                post_shift(&mut circ, &phase1, &phase2, &work, &full, &scratch);
                free_clean(&mut circ, scratch);
            }
            ProofMode::RoundTrip => {
                let pre_scratch = circ.alloc_qreg_bits("ls-parity.baseline.pre", full.len() + 4);
                pre_shift(&mut circ, &phase1, &phase2, &work, &full, &pre_scratch);
                free_clean(&mut circ, pre_scratch);
                let post_scratch = circ.alloc_qreg_bits("ls-parity.baseline.post", full.len() + 4);
                post_shift(&mut circ, &phase1, &phase2, &work, &full, &post_scratch);
                free_clean(&mut circ, post_scratch);
                let post_inverse_scratch =
                    circ.alloc_qreg_bits("ls-parity.baseline.post-inverse", full.len() + 4);
                post_shift_inverse(
                    &mut circ,
                    &phase1,
                    &phase2,
                    &work,
                    &full,
                    &post_inverse_scratch,
                );
                free_clean(&mut circ, post_inverse_scratch);
                let pre_inverse_scratch =
                    circ.alloc_qreg_bits("ls-parity.baseline.pre-inverse", full.len() + 4);
                pre_shift_inverse(
                    &mut circ,
                    &phase1,
                    &phase2,
                    &work,
                    &full,
                    &pre_inverse_scratch,
                );
                free_clean(&mut circ, pre_inverse_scratch);
            }
        }
    }

    let phase1_id = phase1.id();
    let phase2_id = phase2.id();
    let work_ids = work.iter().map(QReg::id).collect::<Vec<_>>();
    let high_ids = high.iter().map(QReg::id).collect::<Vec<_>>();
    let low_or_host_id = low_or_host.id();
    let builder = circ.into_builder();
    let mut external = vec![false; builder.next_qubit as usize];
    for id in std::iter::once(phase1_id)
        .chain(std::iter::once(phase2_id))
        .chain(work_ids.iter().copied())
        .chain(high_ids.iter().copied())
        .chain(std::iter::once(low_or_host_id))
    {
        external[id as usize] = true;
    }
    Harness {
        builder,
        phase1: phase1_id,
        phase2: phase2_id,
        work: work_ids,
        high: high_ids,
        low_or_host: low_or_host_id,
        external,
    }
}

fn splitmix64(mut value: u64) -> u64 {
    value = value.wrapping_add(0x9e37_79b9_7f4a_7c15);
    value = (value ^ (value >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    value = (value ^ (value >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    value ^ (value >> 31)
}

fn state_work_bit(state: InputState, bit: usize) -> bool {
    if let Some(seed) = state.work_seed {
        ((splitmix64(seed ^ ((bit / 64) as u64)) >> (bit % 64)) & 1) != 0
    } else {
        ((state.work >> bit) & 1) != 0
    }
}

fn initial_low(mode: ProofMode, p: bool, phase1: bool) -> bool {
    match mode {
        ProofMode::Pre | ProofMode::RoundTrip => p,
        ProofMode::Post => p ^ true ^ phase1,
    }
}

fn expected_output_low(mode: ProofMode, p: bool, phase1: bool) -> bool {
    match mode {
        ProofMode::Pre => p ^ true ^ phase1,
        ProofMode::Post => !p,
        ProofMode::RoundTrip => p,
    }
}

fn load_harness<R: sha3::digest::XofReader>(
    simulator: &mut Simulator<'_, R>,
    harness: &Harness,
    states: &[InputState],
    mode: ProofMode,
    p: bool,
    baseline: bool,
) -> u64 {
    assert!(states.len() <= 64);
    let active = if states.len() == 64 {
        u64::MAX
    } else {
        (1u64 << states.len()) - 1
    };
    let mask = |predicate: fn(InputState) -> bool| -> u64 {
        states.iter().enumerate().fold(0u64, |bits, (shot, state)| {
            bits | ((predicate(*state) as u64) << shot)
        })
    };
    *simulator.qubit_mut(QubitId(u64::from(harness.phase1))) = mask(|state| state.phase1);
    *simulator.qubit_mut(QubitId(u64::from(harness.phase2))) = mask(|state| state.phase2);
    for (bit, &id) in harness.work.iter().enumerate() {
        let value = states.iter().enumerate().fold(0u64, |bits, (shot, state)| {
            bits | ((state_work_bit(*state, bit) as u64) << shot)
        });
        *simulator.qubit_mut(QubitId(u64::from(id))) = value;
    }
    for (bit, &id) in harness.high.iter().enumerate() {
        let value = states.iter().enumerate().fold(0u64, |bits, (shot, state)| {
            bits | ((((state.high >> bit) & 1) as u64) << shot)
        });
        *simulator.qubit_mut(QubitId(u64::from(id))) = value;
    }
    if baseline {
        let value = states.iter().enumerate().fold(0u64, |bits, (shot, state)| {
            bits | ((initial_low(mode, p, state.phase1) as u64) << shot)
        });
        *simulator.qubit_mut(QubitId(u64::from(harness.low_or_host))) = value;
    }
    active
}

fn new_simulator<'a>(
    harness: &Harness,
    seed: &'static [u8],
    xof: &'a mut sha3::Shake128Reader,
) -> Simulator<'a, sha3::Shake128Reader> {
    let mut hasher = Shake128::default();
    hasher.update(seed);
    *xof = hasher.finalize_xof();
    Simulator::new(
        harness.builder.next_qubit as usize,
        harness.builder.next_bit as usize,
        xof,
    )
}

fn assert_equal_mask(label: &str, left: u64, right: u64, active: u64) {
    let difference = (left ^ right) & active;
    assert_eq!(
        difference,
        0,
        "{label}: first differing shot {}",
        difference.trailing_zeros()
    );
}

fn prove_batch(
    baseline: &Harness,
    candidate: &Harness,
    states: &[InputState],
    mode: ProofMode,
    p: bool,
    report: &mut LsParityProofReport,
) {
    let mut baseline_xof = Shake128::default().finalize_xof();
    let mut candidate_xof = Shake128::default().finalize_xof();
    let mut baseline_sim = new_simulator(baseline, b"ls-parity-baseline", &mut baseline_xof);
    let mut candidate_sim = new_simulator(candidate, b"ls-parity-candidate", &mut candidate_xof);
    let active = load_harness(&mut baseline_sim, baseline, states, mode, p, true);
    assert_eq!(
        load_harness(&mut candidate_sim, candidate, states, mode, p, false),
        active
    );
    baseline_sim.apply_iter(baseline.builder.ops.iter());
    candidate_sim.apply_iter(candidate.builder.ops.iter());

    for (&baseline_id, &candidate_id) in baseline.work.iter().zip(&candidate.work) {
        assert_equal_mask(
            "work",
            baseline_sim.qubit(QubitId(u64::from(baseline_id))),
            candidate_sim.qubit(QubitId(u64::from(candidate_id))),
            active,
        );
        report.logical_register_checks += states.len();
    }
    for (&baseline_id, &candidate_id) in baseline.high.iter().zip(&candidate.high) {
        assert_equal_mask(
            "high l_s",
            baseline_sim.qubit(QubitId(u64::from(baseline_id))),
            candidate_sim.qubit(QubitId(u64::from(candidate_id))),
            active,
        );
        report.logical_register_checks += states.len();
    }
    for (&baseline_id, &candidate_id) in [baseline.phase1, baseline.phase2]
        .iter()
        .zip([candidate.phase1, candidate.phase2].iter())
    {
        assert_equal_mask(
            "controls",
            baseline_sim.qubit(QubitId(u64::from(baseline_id))),
            candidate_sim.qubit(QubitId(u64::from(candidate_id))),
            active,
        );
        report.logical_register_checks += states.len();
    }

    let expected_low = states.iter().enumerate().fold(0u64, |bits, (shot, state)| {
        bits | ((expected_output_low(mode, p, state.phase1) as u64) << shot)
    });
    assert_equal_mask(
        "baseline low l_s",
        baseline_sim.qubit(QubitId(u64::from(baseline.low_or_host))),
        expected_low,
        active,
    );
    assert_equal_mask(
        "candidate host restoration",
        candidate_sim.qubit(QubitId(u64::from(candidate.low_or_host))),
        0,
        active,
    );
    report.host_restoration_checks += states.len();

    assert_equal_mask("phase", baseline_sim.phase, candidate_sim.phase, active);
    assert_equal_mask("baseline phase", baseline_sim.phase, 0, active);
    report.phase_checks += states.len();

    for (id, &external) in baseline.external.iter().enumerate() {
        if !external {
            assert_equal_mask(
                "baseline ancilla",
                baseline_sim.qubit(QubitId(id as u64)),
                0,
                active,
            );
        }
    }
    for (id, &external) in candidate.external.iter().enumerate() {
        if !external {
            assert_equal_mask(
                "candidate ancilla",
                candidate_sim.qubit(QubitId(id as u64)),
                0,
                active,
            );
        }
    }
    report.ancilla_cleanup_checks += states.len();
}

fn emitted_toffoli(builder: &B) -> usize {
    builder
        .ops
        .iter()
        .filter(|op| matches!(op.kind, OperationType::CCX | OperationType::CCZ))
        .count()
}

fn mode_index(mode: ProofMode) -> usize {
    match mode {
        ProofMode::Pre => 0,
        ProofMode::Post => 1,
        ProofMode::RoundTrip => 2,
    }
}

#[must_use]
pub fn prove_ls_parity_shift_representation() -> LsParityProofReport {
    let modes = [ProofMode::Pre, ProofMode::Post, ProofMode::RoundTrip];
    let mut report = LsParityProofReport {
        high_widths_checked: (1..=PRODUCTION_HIGH_WIDTH).collect(),
        modes_checked: modes.len(),
        step_parities_checked: 2,
        reduced_basis_states_checked: 0,
        production_sample_states_checked: 0,
        logical_register_checks: 0,
        host_restoration_checks: 0,
        phase_checks: 0,
        ancilla_cleanup_checks: 0,
        production_baseline_ops: [0; 3],
        production_candidate_ops: [0; 3],
        production_baseline_toffoli: [0; 3],
        production_candidate_toffoli: [0; 3],
    };

    for high_width in 1..=PRODUCTION_HIGH_WIDTH {
        for mode in modes {
            for p in [false, true] {
                let baseline = build_harness(false, mode, p, high_width, REDUCED_WORK_WIDTH);
                let candidate = build_harness(true, mode, p, high_width, REDUCED_WORK_WIDTH);
                let mut states = Vec::new();
                for phase1 in [false, true] {
                    for phase2 in [false, true] {
                        for high in 0..(1usize << high_width) {
                            for work in 0..(1u64 << REDUCED_WORK_WIDTH) {
                                states.push(InputState {
                                    phase1,
                                    phase2,
                                    high,
                                    work,
                                    work_seed: None,
                                });
                            }
                        }
                    }
                }
                for batch in states.chunks(64) {
                    prove_batch(&baseline, &candidate, batch, mode, p, &mut report);
                }
                report.reduced_basis_states_checked += states.len();
            }
        }
    }

    for mode in modes {
        let index = mode_index(mode);
        for p in [false, true] {
            let baseline =
                build_harness(false, mode, p, PRODUCTION_HIGH_WIDTH, PRODUCTION_WORK_WIDTH);
            let candidate =
                build_harness(true, mode, p, PRODUCTION_HIGH_WIDTH, PRODUCTION_WORK_WIDTH);
            report.production_baseline_ops[index] += baseline.builder.ops.len();
            report.production_candidate_ops[index] += candidate.builder.ops.len();
            report.production_baseline_toffoli[index] += emitted_toffoli(&baseline.builder);
            report.production_candidate_toffoli[index] += emitted_toffoli(&candidate.builder);
            let states = (0..64u64)
                .map(|shot| InputState {
                    phase1: (shot & 1) != 0,
                    phase2: (shot & 2) != 0,
                    high: (splitmix64(shot ^ 0x6c73_7061_7269_7479) as usize)
                        & ((1usize << PRODUCTION_HIGH_WIDTH) - 1),
                    work: 0,
                    work_seed: Some(shot ^ 0x7168_6967_682d_7769),
                })
                .collect::<Vec<_>>();
            prove_batch(&baseline, &candidate, &states, mode, p, &mut report);
            report.production_sample_states_checked += states.len();
        }
    }

    report
}
