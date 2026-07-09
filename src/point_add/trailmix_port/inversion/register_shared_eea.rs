//! Classical oracle for the register-sharing EEA of Luo et al. (2026).
//!
//! This module models the published optimized state transition. It is an
//! analysis oracle, not a reversible circuit implementation. In particular,
//! passing this oracle does not establish a challenge-valid qubit count.

use alloy_primitives::U256;
use ruint::aliases::U512;

pub const SECP256K1_BITS: usize = 256;
pub const REGISTER_SHARED_WORK_BITS: usize = SECP256K1_BITS + 3;
pub const REGISTER_SHARED_REFERENCE_STEPS: usize = 1_479;
pub const REGISTER_SHARED_PAPER_INVERSION_QUBITS: usize = 3 * SECP256K1_BITS + 4 * 8 + 20;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
struct Signed512 {
    negative: bool,
    magnitude: U512,
}

impl Signed512 {
    const ZERO: Self = Self {
        negative: false,
        magnitude: U512::ZERO,
    };

    fn unsigned(magnitude: U512) -> Self {
        Self {
            negative: false,
            magnitude,
        }
    }

    fn normalized(mut self) -> Self {
        if self.magnitude.is_zero() {
            self.negative = false;
        }
        self
    }

    fn negated(self) -> Self {
        Self {
            negative: !self.negative,
            magnitude: self.magnitude,
        }
        .normalized()
    }

    fn add(self, other: Self) -> Self {
        if self.negative == other.negative {
            let (magnitude, overflow) = self.magnitude.overflowing_add(other.magnitude);
            assert!(!overflow, "register-sharing signed addition overflow");
            Self {
                negative: self.negative,
                magnitude,
            }
            .normalized()
        } else if self.magnitude >= other.magnitude {
            Self {
                negative: self.negative,
                magnitude: self.magnitude - other.magnitude,
            }
            .normalized()
        } else {
            Self {
                negative: other.negative,
                magnitude: other.magnitude - self.magnitude,
            }
            .normalized()
        }
    }

    fn sub(self, other: Self) -> Self {
        self.add(other.negated())
    }

    fn shifted(self, shift: usize) -> Self {
        assert!(
            bit_len(self.magnitude) + shift <= 512,
            "register-sharing signed shift overflow"
        );
        Self {
            negative: self.negative,
            magnitude: self.magnitude << shift,
        }
        .normalized()
    }

    fn bit_len(self) -> usize {
        bit_len(self.magnitude)
    }

    fn to_mod(self, modulus: U512) -> U512 {
        let residue = self.magnitude % modulus;
        if self.negative && !residue.is_zero() {
            modulus - residue
        } else {
            residue
        }
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
struct Work1 {
    t: Signed512,
    q: U512,
    r: Signed512,
    l_t: usize,
    l_q: usize,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
struct Work2 {
    t_prime: Signed512,
    r_prime: Signed512,
    l_r_prime: usize,
    l_s: usize,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
struct Control {
    phase1: bool,
    phase2: bool,
    iter: bool,
    sign: bool,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
struct State {
    n: usize,
    work1: Work1,
    work2: Work2,
    control: Control,
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct RegisterSharedStepObservation {
    pub step: usize,
    pub r_window: Option<(usize, usize)>,
    pub swap_index: Option<usize>,
    pub t_window_end: Option<usize>,
    pub coefficient_active: bool,
    pub coefficient_sub_enabled: bool,
    pub coefficient_target_above_t: bool,
    pub coefficient_less_than: bool,
    pub coefficient_add_only: bool,
    pub coefficient_t_prime_length_before: usize,
    pub coefficient_shifted_t_length: usize,
    pub coefficient_t_prime_length_after: usize,
    pub length_update: bool,
    pub work1_used: usize,
    pub work2_used: usize,
    pub l_t: usize,
    pub l_q: usize,
    pub l_r_prime_before: usize,
    pub l_r_prime: usize,
    pub transient_l_r_prime: usize,
    pub accepted_remainder_update: bool,
    pub accepted_remainder_strictly_decreased: bool,
    pub accepted_length_nonincreasing: bool,
    pub terminal_padding_length_update: bool,
    pub l_s: usize,
    pub max_intermediate_width: usize,
    pub terminated: bool,
}

/// Canonical packed boundary state for reduced-width gate-level replay.
///
/// This is intentionally limited to small classical oracles whose complete
/// Work1 and Work2 registers fit in `u64`. It lets the reversible port compare
/// whole-step output against the independently implemented integer model.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct RegisterSharedPackedSnapshot {
    pub step: usize,
    pub n: usize,
    pub work1: u64,
    pub work2: u64,
    pub l_t: usize,
    pub l_q: usize,
    pub l_s: usize,
    pub l_r_prime: usize,
    pub phase1: bool,
    pub phase2: bool,
    pub iteration_parity: bool,
    pub sign: bool,
    pub t: u64,
    pub q: u64,
    pub r: u64,
    pub t_prime: u64,
    pub r_prime: u64,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct RegisterSharedFactorAudit {
    pub steps: usize,
    pub euclidean_swaps: usize,
    pub first_termination_step: Option<usize>,
    pub terminal_padding_steps: usize,
    pub maximum_work1_used: usize,
    pub maximum_work2_used: usize,
    pub maximum_shift: usize,
    pub maximum_intermediate_width: usize,
    pub initial_reflected_l_r_prime: usize,
    pub maximum_boundary_l_r_prime: usize,
    pub maximum_transient_l_r_prime: usize,
    pub accepted_remainder_updates: usize,
    pub terminal_padding_length_updates: usize,
    pub accepted_remainder_strict_decrease_failures: usize,
    pub accepted_length_increase_failures: usize,
    pub boundary_l_r_prime_256_observations: usize,
    pub transient_l_r_prime_256_observations: usize,
    pub final_l_t: usize,
    pub final_l_q: usize,
    pub final_l_r_prime: usize,
    pub final_l_s: usize,
    pub inverse_identity_holds: bool,
    pub terminal_gcd_state_holds: bool,
    pub terminal_layout_holds: bool,
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct RegisterSharedNoTransientLengthProofReport {
    pub widths_checked: usize,
    pub moduli_checked: usize,
    pub inputs_checked: usize,
    pub schedule_steps_checked: usize,
    pub accepted_remainder_updates: usize,
    pub terminal_padding_length_updates: usize,
    pub maximum_initial_length: usize,
    pub maximum_boundary_length: usize,
    pub maximum_transient_length: usize,
    pub strict_decrease_failures: usize,
    pub length_increase_failures: usize,
    pub out_of_range_failures: usize,
}

fn bit_len(value: U512) -> usize {
    if value.is_zero() {
        0
    } else {
        512 - value.leading_zeros() as usize
    }
}

fn widen(value: U256) -> U512 {
    let limbs = value.as_limbs();
    U512::from_limbs([limbs[0], limbs[1], limbs[2], limbs[3], 0, 0, 0, 0])
}

fn secp256k1_modulus() -> U512 {
    (U512::from(1u64) << 256) - (U512::from(1u64) << 32) - U512::from(977u64)
}

fn set_bit(value: &mut U512, index: usize, bit: bool) {
    assert!(index < 512);
    let mask = U512::from(1u64) << index;
    let old = !(*value & mask).is_zero();
    match (old, bit) {
        (false, true) => *value |= mask,
        (true, false) => *value -= mask,
        _ => {}
    }
}

fn small_u64(value: U512) -> u64 {
    let limbs = value.as_limbs();
    assert!(limbs[1..].iter().all(|&limb| limb == 0));
    limbs[0]
}

fn rotate_low_small(value: u64, width: usize, amount: usize) -> u64 {
    assert!(width > 0 && width < 64);
    let amount = amount % width;
    let mask = (1u64 << width) - 1;
    if amount == 0 {
        value & mask
    } else {
        ((value >> amount) | (value << (width - amount))) & mask
    }
}

impl State {
    fn new(modulus: U512, input: U512, n: usize) -> Self {
        assert!(!input.is_zero() && input < modulus);
        let half = modulus >> 1;
        let reflected = input > half;
        let adjusted = if reflected { modulus - input } else { input };
        assert!(
            bit_len(adjusted) < n,
            "reflected denominator length must fit in n-1 bits"
        );
        let state = Self {
            n,
            work1: Work1 {
                t: Signed512::unsigned(U512::from(1u64)),
                q: U512::ZERO,
                r: Signed512::unsigned(modulus),
                l_t: 1,
                l_q: 0,
            },
            work2: Work2 {
                t_prime: Signed512::ZERO,
                r_prime: Signed512::unsigned(adjusted),
                l_r_prime: bit_len(adjusted),
                l_s: 0,
            },
            control: Control {
                phase1: false,
                phase2: false,
                iter: reflected,
                sign: false,
            },
        };
        state.assert_boundary_invariants();
        state
    }

    fn assert_boundary_invariants(&self) -> (usize, usize) {
        assert!(!self.work1.t.negative);
        assert!(!self.work1.r.negative);
        assert!(!self.work2.t_prime.negative);
        assert!(!self.work2.r_prime.negative);
        assert_eq!(self.work1.l_t, self.work1.t.bit_len());
        assert_eq!(self.work2.l_r_prime, self.work2.r_prime.bit_len());
        let work1_used = self.work1.l_t + 1 + self.work1.l_q + self.work1.r.bit_len();
        let work2_used = self.work2.t_prime.bit_len() + self.work2.l_r_prime;
        assert!(
            work1_used <= self.n + 3,
            "Work1 packing overflow: {work1_used} > {}",
            self.n + 3
        );
        assert!(
            work2_used <= self.n + 3,
            "Work2 packing overflow: {work2_used} > {}",
            self.n + 3
        );
        assert!(self.work2.l_s < (1usize << 10));
        (work1_used, work2_used)
    }

    fn packed_snapshot(&self, step: usize) -> RegisterSharedPackedSnapshot {
        let width = self.n + 3;
        assert!(width < 64);
        self.assert_boundary_invariants();
        assert!(!self.work1.t.negative);
        assert!(!self.work1.r.negative);
        assert!(!self.work2.t_prime.negative);
        assert!(!self.work2.r_prime.negative);

        let t = small_u64(self.work1.t.magnitude);
        let q = small_u64(self.work1.q);
        let r = small_u64(self.work1.r.magnitude);
        let t_prime = small_u64(self.work2.t_prime.magnitude);
        let r_prime = small_u64(self.work2.r_prime.magnitude);

        let mut work1 = t;
        let q_width = bit_len(self.work1.q);
        if self.work1.l_q > 0 {
            assert!(q_width >= self.work1.l_q);
            let active_q = q >> (q_width - self.work1.l_q);
            for physical_index in 0..self.work1.l_q {
                let source_index = self.work1.l_q - 1 - physical_index;
                if ((active_q >> source_index) & 1) != 0 {
                    work1 |= 1u64 << (self.work1.l_t + 1 + physical_index);
                }
            }
        }
        for bit in 0..bit_len(self.work1.r.magnitude) {
            if ((r >> bit) & 1) != 0 {
                work1 |= 1u64 << (width - 1 - bit);
            }
        }

        let mut work2_raw = t_prime;
        for bit in 0..self.work2.l_r_prime {
            if ((r_prime >> bit) & 1) != 0 {
                work2_raw |= 1u64 << (width - 1 - bit);
            }
        }
        let work2 = rotate_low_small(work2_raw, width, self.work2.l_s);

        RegisterSharedPackedSnapshot {
            step,
            n: self.n,
            work1,
            work2,
            l_t: self.work1.l_t,
            l_q: self.work1.l_q,
            l_s: self.work2.l_s,
            l_r_prime: self.work2.l_r_prime,
            phase1: self.control.phase1,
            phase2: self.control.phase2,
            iteration_parity: self.control.iter,
            sign: self.control.sign,
            t,
            q,
            r,
            t_prime,
            r_prime,
        }
    }

    fn step(&mut self, step: usize) -> RegisterSharedStepObservation {
        let mut observation = RegisterSharedStepObservation {
            step,
            l_r_prime_before: self.work2.l_r_prime,
            transient_l_r_prime: self.work2.l_r_prime,
            ..Default::default()
        };

        if !self.control.phase1 {
            self.work2.l_s += 1;
        }
        if !self.control.phase1 && self.control.phase2 {
            assert!(self.work2.l_s >= 2);
            self.work2.l_s -= 2;
        }

        if !self.control.phase1 && self.work2.l_r_prime > 0 {
            let start = self.work1.l_t + self.work1.l_q + 2;
            let end = self.n + 3 - self.work2.l_s;
            assert!(start <= end);
            observation.r_window = Some((start, end));
            let shifted = self.work2.r_prime.shifted(self.work2.l_s);
            self.work1.r = self.work1.r.sub(shifted);
            observation.max_intermediate_width = observation
                .max_intermediate_width
                .max(self.work1.r.bit_len());
            self.control.sign ^= self.work1.r.negative;
        }
        if !self.control.phase1 && self.control.phase2 && self.work2.l_r_prime > 0 {
            self.control.sign ^= true;
        }
        if !self.control.phase1
            && self.work2.l_r_prime > 0
            && (!self.control.phase2 || !self.control.sign)
        {
            let shifted = self.work2.r_prime.shifted(self.work2.l_s);
            self.work1.r = self.work1.r.add(shifted);
        }

        self.control.phase2 ^= self.control.phase1;
        if self.control.phase2 {
            let index = self.work1.l_t + self.work1.l_q;
            assert!(index < self.n + 3);
            observation.swap_index = Some(index + 1);
            let quotient_bit = !((self.work1.q >> self.work2.l_s) & U512::from(1u64)).is_zero();
            let sign = self.control.sign;
            set_bit(&mut self.work1.q, self.work2.l_s, sign);
            self.control.sign = quotient_bit;
            if self.control.phase1 {
                assert!(self.work1.l_q > 0);
                self.work1.l_q -= 1;
            } else {
                self.work1.l_q += 1;
            }
        }
        self.control.phase2 ^= self.control.phase1;

        if self.control.phase1 {
            let shifted_t = self.work1.t.shifted(self.work2.l_s);
            observation.coefficient_active = true;
            observation.coefficient_sub_enabled = self.control.phase2 || !self.control.sign;
            observation.coefficient_add_only = !observation.coefficient_sub_enabled;
            observation.coefficient_t_prime_length_before = self.work2.t_prime.bit_len();
            observation.coefficient_shifted_t_length = self.work1.l_t + self.work2.l_s;
            observation.coefficient_target_above_t =
                self.work2.t_prime.bit_len() > self.work1.l_t + self.work2.l_s;
            observation.coefficient_less_than = observation.coefficient_sub_enabled
                && self.work2.t_prime.magnitude < shifted_t.magnitude;
        }
        if self.control.phase1 && (self.control.phase2 || !self.control.sign) {
            self.work2.t_prime = self.work2.t_prime.sub(self.work1.t.shifted(self.work2.l_s));
            observation.max_intermediate_width = observation
                .max_intermediate_width
                .max(self.work2.t_prime.bit_len());
        }
        if self.control.phase1 {
            observation.t_window_end = Some(self.work1.l_t + 1);
            self.control.sign ^= true;
            self.control.sign ^= self.work2.t_prime.negative;
            self.work2.t_prime = self.work2.t_prime.add(self.work1.t.shifted(self.work2.l_s));
            observation.coefficient_t_prime_length_after = self.work2.t_prime.bit_len();
        }

        if self.control.phase1 {
            self.work2.l_s += 1;
        }
        if self.control.phase1 && self.control.phase2 {
            assert!(self.work2.l_s >= 2);
            self.work2.l_s -= 2;
        }

        if self.work1.l_q == 0 && self.work2.l_r_prime > 0 {
            self.control.phase2 ^= self.control.sign ^ self.control.phase1;
            self.control.sign ^= self.control.phase2;
        }
        if self.work2.l_s == 0 {
            self.control.phase1 ^= true;
            self.control.phase2 ^= true;
        }

        if self.work1.l_q == 0 && self.work2.l_s == 0 {
            observation.length_update = true;
            assert!(self.work1.q.is_zero());
            assert!(!self.work1.r.negative);
            assert!(!self.work2.r_prime.negative);
            let old_r_prime = self.work2.r_prime;
            let old_l_r_prime = self.work2.l_r_prime;
            let next_l_r_prime = self.work1.r.bit_len();
            observation.transient_l_r_prime = next_l_r_prime;
            observation.accepted_remainder_update = old_l_r_prime > 0;
            observation.terminal_padding_length_update = old_l_r_prime == 0;
            observation.accepted_remainder_strictly_decreased =
                old_l_r_prime == 0 || self.work1.r.magnitude < old_r_prime.magnitude;
            observation.accepted_length_nonincreasing =
                old_l_r_prime == 0 || next_l_r_prime <= old_l_r_prime;
            std::mem::swap(&mut self.work1.t, &mut self.work2.t_prime);
            std::mem::swap(&mut self.work1.r, &mut self.work2.r_prime);
            self.work1.l_t = self.work1.t.bit_len();
            self.work2.l_r_prime = self.work2.r_prime.bit_len();
            self.control.iter ^= true;
        }

        let (work1_used, work2_used) = self.assert_boundary_invariants();
        observation.work1_used = work1_used;
        observation.work2_used = work2_used;
        observation.l_t = self.work1.l_t;
        observation.l_q = self.work1.l_q;
        observation.l_r_prime = self.work2.l_r_prime;
        observation.l_s = self.work2.l_s;
        observation.max_intermediate_width = observation
            .max_intermediate_width
            .max(self.work1.r.bit_len())
            .max(self.work2.t_prime.bit_len());
        observation.terminated = self.work2.l_r_prime == 0;
        observation
    }
}

/// Generate an exact reduced-width packed trace from the classical oracle.
#[must_use]
pub fn register_shared_small_packed_trace(
    input: u64,
    modulus: u64,
    n: usize,
    steps: usize,
) -> Vec<RegisterSharedPackedSnapshot> {
    assert!(n > 0 && n + 3 < 64);
    assert!(input > 0 && input < modulus);
    let mut state = State::new(U512::from(modulus), U512::from(input), n);
    let mut trace = Vec::with_capacity(steps + 1);
    trace.push(state.packed_snapshot(0));
    for step in 1..=steps {
        state.step(step);
        trace.push(state.packed_snapshot(step));
    }
    trace
}

fn run_with_modulus<F>(
    input: U512,
    modulus: U512,
    n: usize,
    steps: usize,
    mut observe: F,
) -> RegisterSharedFactorAudit
where
    F: FnMut(RegisterSharedStepObservation),
{
    let mut state = State::new(modulus, input, n);
    let mut swaps = 0usize;
    let mut first_termination_step = None;
    let mut maximum_work1_used = 0usize;
    let mut maximum_work2_used = 0usize;
    let mut maximum_shift = 0usize;
    let mut maximum_intermediate_width = 0usize;
    let initial_reflected_l_r_prime = state.work2.l_r_prime;
    let mut maximum_boundary_l_r_prime = initial_reflected_l_r_prime;
    let mut maximum_transient_l_r_prime = initial_reflected_l_r_prime;
    let mut accepted_remainder_updates = 0usize;
    let mut terminal_padding_length_updates = 0usize;
    let mut accepted_remainder_strict_decrease_failures = 0usize;
    let mut accepted_length_increase_failures = 0usize;
    let mut boundary_l_r_prime_256_observations =
        usize::from(initial_reflected_l_r_prime == SECP256K1_BITS);
    let mut transient_l_r_prime_256_observations =
        usize::from(initial_reflected_l_r_prime == SECP256K1_BITS);
    for step in 1..=steps {
        let observation = state.step(step);
        swaps += usize::from(observation.length_update);
        if observation.terminated && first_termination_step.is_none() {
            first_termination_step = Some(step);
        }
        maximum_work1_used = maximum_work1_used.max(observation.work1_used);
        maximum_work2_used = maximum_work2_used.max(observation.work2_used);
        maximum_shift = maximum_shift.max(observation.l_s);
        maximum_intermediate_width =
            maximum_intermediate_width.max(observation.max_intermediate_width);
        maximum_boundary_l_r_prime = maximum_boundary_l_r_prime
            .max(observation.l_r_prime_before)
            .max(observation.l_r_prime);
        maximum_transient_l_r_prime =
            maximum_transient_l_r_prime.max(observation.transient_l_r_prime);
        accepted_remainder_updates += usize::from(observation.accepted_remainder_update);
        terminal_padding_length_updates += usize::from(observation.terminal_padding_length_update);
        accepted_remainder_strict_decrease_failures += usize::from(
            observation.accepted_remainder_update
                && !observation.accepted_remainder_strictly_decreased,
        );
        accepted_length_increase_failures += usize::from(
            observation.accepted_remainder_update && !observation.accepted_length_nonincreasing,
        );
        boundary_l_r_prime_256_observations += usize::from(
            observation.l_r_prime_before == SECP256K1_BITS
                || observation.l_r_prime == SECP256K1_BITS,
        );
        transient_l_r_prime_256_observations +=
            usize::from(observation.transient_l_r_prime == SECP256K1_BITS);
        observe(observation);
    }

    let signed_inverse = if state.control.iter {
        state.work2.t_prime
    } else {
        state.work2.t_prime.negated()
    };
    let inverse = signed_inverse.to_mod(modulus);
    let inverse_identity_holds = (input * inverse) % modulus == U512::from(1u64);
    let terminal_gcd_state_holds = state.work1.r == Signed512::unsigned(U512::from(1u64))
        && state.work2.r_prime == Signed512::ZERO
        && state.work1.q.is_zero()
        && state.work1.l_q == 0
        && state.work2.l_r_prime == 0;
    let terminal_layout_holds = terminal_gcd_state_holds
        && state.work1.t == Signed512::unsigned(modulus)
        && state.work1.l_t == bit_len(modulus)
        && !state.work2.t_prime.negative
        && !state.control.phase1
        && !state.control.phase2
        && !state.control.sign;
    RegisterSharedFactorAudit {
        steps,
        euclidean_swaps: swaps,
        first_termination_step,
        terminal_padding_steps: first_termination_step
            .map_or(0, |termination| steps.saturating_sub(termination)),
        maximum_work1_used,
        maximum_work2_used,
        maximum_shift,
        maximum_intermediate_width,
        initial_reflected_l_r_prime,
        maximum_boundary_l_r_prime,
        maximum_transient_l_r_prime,
        accepted_remainder_updates,
        terminal_padding_length_updates,
        accepted_remainder_strict_decrease_failures,
        accepted_length_increase_failures,
        boundary_l_r_prime_256_observations,
        transient_l_r_prime_256_observations,
        final_l_t: state.work1.l_t,
        final_l_q: state.work1.l_q,
        final_l_r_prime: state.work2.l_r_prime,
        final_l_s: state.work2.l_s,
        inverse_identity_holds,
        terminal_gcd_state_holds,
        terminal_layout_holds,
    }
}

pub fn audit_secp256k1_factor<F>(factor: U256, observe: F) -> RegisterSharedFactorAudit
where
    F: FnMut(RegisterSharedStepObservation),
{
    run_with_modulus(
        widen(factor),
        secp256k1_modulus(),
        SECP256K1_BITS,
        REGISTER_SHARED_REFERENCE_STEPS,
        observe,
    )
}

/// Exhaustively prove the reflected-denominator and accepted-remainder length
/// bounds for every reduced-width modulus and nonzero input through eight bits.
#[must_use]
pub fn exhaustive_no_transient_remainder_length_check() -> RegisterSharedNoTransientLengthProofReport
{
    let mut report = RegisterSharedNoTransientLengthProofReport::default();
    for n in 2usize..=8 {
        report.widths_checked += 1;
        let lower = (1u64 << (n - 1)) + 1;
        let upper = 1u64 << n;
        for modulus in (lower..upper).filter(|value| value & 1 == 1) {
            report.moduli_checked += 1;
            for input in 1..modulus {
                let steps = 12 * n;
                let audit =
                    run_with_modulus(U512::from(input), U512::from(modulus), n, steps, |_| {});
                report.inputs_checked += 1;
                report.schedule_steps_checked += steps;
                report.accepted_remainder_updates += audit.accepted_remainder_updates;
                report.terminal_padding_length_updates += audit.terminal_padding_length_updates;
                report.maximum_initial_length = report
                    .maximum_initial_length
                    .max(audit.initial_reflected_l_r_prime);
                report.maximum_boundary_length = report
                    .maximum_boundary_length
                    .max(audit.maximum_boundary_l_r_prime);
                report.maximum_transient_length = report
                    .maximum_transient_length
                    .max(audit.maximum_transient_l_r_prime);
                report.strict_decrease_failures +=
                    audit.accepted_remainder_strict_decrease_failures;
                report.length_increase_failures += audit.accepted_length_increase_failures;
                report.out_of_range_failures += usize::from(
                    audit.initial_reflected_l_r_prime >= n
                        || audit.maximum_boundary_l_r_prime >= n
                        || audit.maximum_transient_l_r_prime >= n,
                );
            }
        }
    }
    report
}

pub fn register_shared_eea_selftest() {
    let seven = Signed512::unsigned(U512::from(7u64));
    let eleven = Signed512::unsigned(U512::from(11u64));
    assert_eq!(seven.sub(eleven).add(eleven), seven);
    assert_eq!(
        eleven.sub(seven).sub(seven),
        Signed512::unsigned(U512::from(3u64)).negated()
    );

    let audit = run_with_modulus(U512::from(13u64), U512::from(37u64), 6, 36, |_| {});
    assert!(audit.inverse_identity_holds);
    assert!(audit.terminal_gcd_state_holds);
    assert!(audit.maximum_work1_used <= 9);
    assert!(audit.maximum_work2_used <= 9);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn explicit_selftest_passes() {
        register_shared_eea_selftest();
    }
}
