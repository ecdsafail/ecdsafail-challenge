/// This file contains code for simulating kickmix circuits.
use crate::circuit::{BitId, Op, OperationType, QubitId, NO_BIT};
use ruint::aliases::U256;

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct SimStats {
    pub clifford_gates: u64,
    pub toffoli_gates: u64,
}

pub struct Simulator<'a, R: sha3::digest::XofReader> {
    pub phase: u64,
    pub qubits: Vec<u64>,
    pub bits: Vec<u64>,
    pub num_qubits: usize,
    pub num_bits: usize,
    pub xof: &'a mut R,
    pub stats: SimStats,
}

impl<'a, R: sha3::digest::XofReader> Simulator<'a, R> {
    pub fn new(num_qubits: usize, num_bits: usize, xof: &'a mut R) -> Self {
        let qubits = vec![0; num_qubits];
        let bits = vec![0; num_bits];

        Self {
            phase: 0,
            qubits,
            bits,
            num_qubits,
            num_bits,
            xof,
            stats: SimStats::default(),
        }
    }

    #[inline(always)]
    pub fn qubit(&self, id: QubitId) -> u64 {
        self.qubits[id.0 as usize]
    }

    #[inline(always)]
    pub fn qubit_mut(&mut self, id: QubitId) -> &mut u64 {
        &mut self.qubits[id.0 as usize]
    }

    #[inline(always)]
    pub fn bit(&self, id: BitId) -> u64 {
        self.bits[id.0 as usize]
    }

    #[inline(always)]
    pub fn bit_mut(&mut self, id: BitId) -> &mut u64 {
        &mut self.bits[id.0 as usize]
    }

    pub fn clear_for_shot(&mut self) {
        for e in &mut self.qubits {
            *e = 0;
        }
        for e in &mut self.bits {
            *e = 0;
        }
        self.phase = 0;
    }

    pub fn apply_iter<'b>(&mut self, ops: impl Iterator<Item = &'b Op>) {
        let mut condition_stack = Vec::new();
        let mut current_base_condition = u64::MAX;

        for op in ops {
            let mut cond = current_base_condition;
            if op.c_condition != NO_BIT {
                cond &= self.bit(op.c_condition);
            }

            let executed_shots = cond.count_ones() as u64;

            match op.kind {
                OperationType::CCZ | OperationType::CCX => {
                    self.stats.toffoli_gates += executed_shots;
                }
                OperationType::CX
                | OperationType::CZ
                | OperationType::Swap
                | OperationType::R
                | OperationType::Hmr => {
                    self.stats.clifford_gates += executed_shots;
                }
                // Note: X and Z are not considered Clifford gates in the
                // stats because they can be tracked in the classical control system.
                // They don't need to cause something to happen on the quantum computer.
                _ => {}
            }

            match op.kind {
                OperationType::CCX => {
                    let v = cond & self.qubit(op.q_control1) & self.qubit(op.q_control2);
                    *self.qubit_mut(op.q_target) ^= v;
                }
                OperationType::CX => {
                    let v = cond & self.qubit(op.q_control1);
                    *self.qubit_mut(op.q_target) ^= v;
                }
                OperationType::Swap => {
                    let mut q_c1 = self.qubit(op.q_control1);
                    let mut q_t = self.qubit(op.q_target);
                    q_c1 ^= q_t;
                    q_t ^= cond & q_c1;
                    q_c1 ^= q_t;
                    *self.qubit_mut(op.q_control1) = q_c1;
                    *self.qubit_mut(op.q_target) = q_t;
                }
                OperationType::X => {
                    *self.qubit_mut(op.q_target) ^= cond;
                }
                OperationType::CCZ => {
                    let v = cond
                        & self.qubit(op.q_target)
                        & self.qubit(op.q_control1)
                        & self.qubit(op.q_control2);
                    self.phase ^= v;
                }
                OperationType::CZ => {
                    let v = cond & self.qubit(op.q_target) & self.qubit(op.q_control1);
                    self.phase ^= v;
                }
                OperationType::Z => {
                    let v = cond & self.qubit(op.q_target);
                    self.phase ^= v;
                }
                OperationType::Neg => {
                    self.phase ^= cond;
                }
                OperationType::Hmr => {
                    let mut buf = [0u8; 8];
                    self.xof.read(&mut buf);
                    let rng_val = u64::from_le_bytes(buf);
                    *self.bit_mut(op.c_target) &= !cond;
                    *self.bit_mut(op.c_target) ^= rng_val & cond;
                    self.phase ^= self.qubit(op.q_target) & rng_val & cond;
                    *self.qubit_mut(op.q_target) &= !cond;
                }
                OperationType::R => {
                    let mut buf = [0u8; 8];
                    self.xof.read(&mut buf);
                    let rng_val = u64::from_le_bytes(buf);
                    self.phase ^= self.qubit(op.q_target) & rng_val & cond;
                    *self.qubit_mut(op.q_target) &= !cond;
                }
                OperationType::BitInvert => {
                    *self.bit_mut(op.c_target) ^= cond;
                }
                OperationType::BitStore0 => {
                    *self.bit_mut(op.c_target) &= !cond;
                }
                OperationType::BitStore1 => {
                    *self.bit_mut(op.c_target) |= cond;
                }
                OperationType::AppendToRegister
                | OperationType::Register
                | OperationType::DebugPrint => {}
                OperationType::PushCondition => {
                    condition_stack.push(current_base_condition);
                    current_base_condition &= self.bit(op.c_condition);
                }
                OperationType::PopCondition => {
                    if let Some(val) = condition_stack.pop() {
                        current_base_condition = val;
                    }
                }
            }
        }
    }

    /// Writes an integer into the qubits/bits of a register.
    ///
    /// Args:
    ///     reg: The qubits and bits making up the register, in little endian order.
    ///         CAUTION: Writes are unchecked!
    ///             Only pass in bits and qubits consistent with num_bits and num_qubits!
    ///         Caution: if a qubit or bit appears multiple times, the write to the more
    ///             significant bit position will overwrite prior writes.
    ///     val: The value to write into the bits/qubits.
    ///     shot_idx: The simulator tracks 64 shots in parallel. This is which shot to write to.
    pub fn set_register(&mut self, reg: &[crate::circuit::QubitOrBit], val: U256, shot_idx: usize) {
        for (i, item) in reg.iter().enumerate() {
            let bit_val = val.bit(i);
            match item {
                crate::circuit::QubitOrBit::Qubit(id) => {
                    if bit_val {
                        *self.qubit_mut(*id) |= 1 << shot_idx;
                    } else {
                        *self.qubit_mut(*id) &= !(1 << shot_idx);
                    }
                }
                crate::circuit::QubitOrBit::Bit(id) => {
                    if bit_val {
                        *self.bit_mut(*id) |= 1 << shot_idx;
                    } else {
                        *self.bit_mut(*id) &= !(1 << shot_idx);
                    }
                }
            }
        }
    }

    /// Reads the qubits/bits of a register as an integer.
    ///
    /// Args:
    ///     reg: The qubits and bits making up the register, in little endian order.
    ///         CAUTION: Reads are unchecked!
    ///             Only pass in bits and qubits consistent with num_bits and num_qubits!
    ///     shot_idx: The simulator tracks 64 shots in parallel. This is which shot to read from.
    ///
    /// Returns:
    ///     The requested integer.
    pub fn get_register(&self, reg: &[crate::circuit::QubitOrBit], shot_idx: usize) -> U256 {
        let mut v = U256::ZERO;
        for (i, item) in reg.iter().enumerate() {
            let bit_val = match item {
                crate::circuit::QubitOrBit::Qubit(id) => (self.qubit(*id) >> shot_idx) & 1,
                crate::circuit::QubitOrBit::Bit(id) => (self.bit(*id) >> shot_idx) & 1,
            };
            v.set_bit(i, bit_val != 0);
        }
        v
    }
}
