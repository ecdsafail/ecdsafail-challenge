/// This file contains code for working with kickmix circuit files.

use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;
use std::num::ParseIntError;

fn parse_u64_below_max(s: &str) -> Result<u64, ParseIntError> {
    let result: u64 = s.parse().unwrap();
    if result == u64::MAX {
        return Err("".parse::<u64>().unwrap_err());
    }
    return Ok(result);
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum OperationType {
    /// Global phase flip.
    Neg = 0,
    /// Ensure a register exists.
    Register = 1,
    /// Annotates that a classical bit or qubit is part of a register.
    AppendToRegister = 2,
    /// Inverts a bit.
    BitInvert = 3,
    /// Writes 0 to a bit.
    BitStore0 = 4,
    /// Writes 1 to a bit.
    BitStore1 = 5,
    /// NOT gate.
    X = 6,
    /// Phase-flip states where a qubit is 1.
    Z = 7,
    /// CNOT gate.
    CX = 8,
    /// Phase-flips states where two qubits are both 1.
    CZ = 9,
    /// Exchanges two qubits.
    Swap = 10,
    /// Reset. Equivalent to HMR with an ignored measurement result.
    R = 11,
    /// X-basis measurement combined with demolition into the 0 state.
    Hmr = 12,
    /// Toffoli gate.
    CCX = 13,
    /// Phase-flips states where three qubits are all 1.
    CCZ = 14,
   /// Pushes a bit onto the condition stack.
    /// (Operations other than PUSH_CONDITION/POP_CONDITION do not
    /// occur unless all values on the condition stack are True.)
    PushCondition = 15,
    /// Pops a bit off of the condition stack.
    /// (Operations other than PUSH_CONDITION/POP_CONDITION do not
    /// occur unless all values on the condition stack are True.)
    PopCondition = 16,
    /// No effect on the simulation. Hints that a value should be
    /// printed, for debugging purposes.
    DebugPrint = 17,
}

impl OperationType {
    pub fn from_name(name: &str) -> Option<Self> {
        match name {
            "NEG" => Some(Self::Neg),
            "REGISTER" => Some(Self::Register),
            "APPEND_TO_REGISTER" => Some(Self::AppendToRegister),
            "BIT_INVERT" => Some(Self::BitInvert),
            "BIT_STORE0" => Some(Self::BitStore0),
            "BIT_STORE1" => Some(Self::BitStore1),
            "X" => Some(Self::X),
            "Z" => Some(Self::Z),
            "CX" => Some(Self::CX),
            "CZ" => Some(Self::CZ),
            "SWAP" => Some(Self::Swap),
            "R" => Some(Self::R),
            "HMR" => Some(Self::Hmr),
            "CCX" => Some(Self::CCX),
            "CCZ" => Some(Self::CCZ),
            "PUSH_CONDITION" => Some(Self::PushCondition),
            "POP_CONDITION" => Some(Self::PopCondition),
            "DEBUG_PRINT" => Some(Self::DebugPrint),
            _ => None,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct QubitId(pub u64);
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct BitId(pub u64);
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct RegisterId(pub u64);

pub const NO_QUBIT: QubitId = QubitId(u64::MAX);
pub const NO_BIT: BitId = BitId(u64::MAX);
pub const NO_REG: RegisterId = RegisterId(u64::MAX);

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum QubitOrBit {
    Qubit(QubitId),
    Bit(BitId),
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Op {
    pub kind: OperationType,
    pub q_control2: QubitId,
    pub q_control1: QubitId,
    pub q_target: QubitId,
    pub c_target: BitId,
    pub c_condition: BitId,
    pub r_target: RegisterId,
}

impl Op {
    pub fn empty() -> Self {
        Self {
            kind: OperationType::Neg,
            q_control2: NO_QUBIT,
            q_control1: NO_QUBIT,
            q_target: NO_QUBIT,
            c_target: NO_BIT,
            c_condition: NO_BIT,
            r_target: NO_REG,
        }
    }

    pub fn validate(&self) {
        // Check for qubit aliasing.
        if self.q_target == self.q_control1 && self.q_target != NO_QUBIT {
            panic!("kind={:?} and q_target==q_control1==q{}", self.kind, self.q_target.0);
        }
        if self.q_target == self.q_control2 && self.q_target != NO_QUBIT {
            panic!("kind={:?} and q_target==q_control2==q{}", self.kind, self.q_target.0);
        }
        if self.q_control1 == self.q_control2 && self.q_control1 != NO_QUBIT {
            panic!("kind={:?} and q_control1==q_control2==q{}", self.kind, self.q_control1.0);
        }

        const BANNED: u8 = 0;
        const ALLOWED: u8 = 1;
        const REQUIRED: u8 = 2;

        let mut q_target_flag = BANNED;
        let mut q_control1_flag = BANNED;
        let mut q_control2_flag = BANNED;
        let mut c_target_flag = BANNED;
        let mut r_target_flag = BANNED;
        let mut c_condition_flag = BANNED;

        match self.kind {
            OperationType::DebugPrint => return,
            OperationType::Register => {
                r_target_flag = REQUIRED;
            }
            OperationType::AppendToRegister => {
                if (self.q_target == NO_QUBIT) == (self.c_target == NO_BIT) {
                    panic!("kind={:?} needs exactly one qubit target or bit target", self.kind);
                }
                c_target_flag = ALLOWED;
                q_target_flag = ALLOWED;
                r_target_flag = REQUIRED;
            }
            OperationType::CCX | OperationType::CCZ => {
                c_condition_flag = ALLOWED;
                q_target_flag = REQUIRED;
                q_control1_flag = REQUIRED;
                q_control2_flag = REQUIRED;
            }
            OperationType::CX | OperationType::CZ | OperationType::Swap => {
                c_condition_flag = ALLOWED;
                q_target_flag = REQUIRED;
                q_control1_flag = REQUIRED;
            }
            OperationType::X | OperationType::Z | OperationType::R => {
                c_condition_flag = ALLOWED;
                q_target_flag = REQUIRED;
            }
            OperationType::Neg => {
                c_condition_flag = ALLOWED;
            }
            OperationType::Hmr => {
                c_condition_flag = ALLOWED;
                q_target_flag = REQUIRED;
                c_target_flag = REQUIRED;
            }
            OperationType::BitInvert | OperationType::BitStore0 | OperationType::BitStore1 => {
                c_condition_flag = ALLOWED;
                c_target_flag = REQUIRED;
            }
            OperationType::PushCondition => {
                c_condition_flag = REQUIRED;
            }
            OperationType::PopCondition => {}
        }

        if c_condition_flag == REQUIRED && self.c_condition == NO_BIT {
            panic!("kind={:?} but c_condition == NO_BIT", self.kind);
        } else if c_condition_flag == BANNED && self.c_condition != NO_BIT {
            panic!("kind={:?} but c_condition != NO_BIT", self.kind);
        }

        if q_target_flag == REQUIRED && self.q_target == NO_QUBIT {
            panic!("kind={:?} but q_target == NO_QUBIT", self.kind);
        } else if q_target_flag == BANNED && self.q_target != NO_QUBIT {
            panic!("kind={:?} but q_target != NO_QUBIT", self.kind);
        }

        if q_control1_flag == REQUIRED && self.q_control1 == NO_QUBIT {
            panic!("kind={:?} but q_control1 == NO_QUBIT", self.kind);
        } else if q_control1_flag == BANNED && self.q_control1 != NO_QUBIT {
            panic!("kind={:?} but q_control1 != NO_QUBIT", self.kind);
        }

        if q_control2_flag == REQUIRED && self.q_control2 == NO_QUBIT {
            panic!("kind={:?} but q_control2 == NO_QUBIT", self.kind);
        } else if q_control2_flag == BANNED && self.q_control2 != NO_QUBIT {
            panic!("kind={:?} but q_control2 != NO_QUBIT", self.kind);
        }

        if c_target_flag == REQUIRED && self.c_target == NO_BIT {
            panic!("kind={:?} but c_target == NO_BIT", self.kind);
        } else if c_target_flag == BANNED && self.c_target != NO_BIT {
            panic!("kind={:?} but c_target != NO_BIT", self.kind);
        }

        if r_target_flag == REQUIRED && self.r_target == NO_REG {
            panic!("kind={:?} but r_target == NO_REG", self.kind);
        } else if r_target_flag == BANNED && self.r_target != NO_REG {
            panic!("kind={:?} but r_target != NO_REG", self.kind);
        }
    }

    pub fn from_text(line: &str) -> Option<Self> {
        let words: Vec<&str> = line.split_whitespace().collect();
        if words.is_empty() || words[0].starts_with('#') {
            return None;
        }

        let mut out = Self::empty();

        if let Some(kind) = OperationType::from_name(words[0]) {
            out.kind = kind;
        } else {
            panic!("Unrecognized operation type '{}'", words[0]);
        }

        let mut cur_word = 1;

        if cur_word < words.len() && words[cur_word].starts_with('q') {
            out.q_target.0 = parse_u64_below_max(&words[cur_word][1..]).unwrap();
            cur_word += 1;

            if cur_word < words.len() && words[cur_word].starts_with('q') {
                out.q_control1 = out.q_target;
                out.q_target.0 = parse_u64_below_max(&words[cur_word][1..]).unwrap();
                cur_word += 1;
            }

            if cur_word < words.len() && words[cur_word].starts_with('q') {
                out.q_control2 = out.q_control1;
                out.q_control1 = out.q_target;
                out.q_target.0 = parse_u64_below_max(&words[cur_word][1..]).unwrap();
                cur_word += 1;
            }
        }

        if cur_word < words.len() && words[cur_word].starts_with('b') {
            out.c_target.0 = parse_u64_below_max(&words[cur_word][1..]).unwrap();
            cur_word += 1;
        }
        if cur_word < words.len() && words[cur_word].starts_with('r') {
            out.r_target.0 = parse_u64_below_max(&words[cur_word][1..]).unwrap();
            cur_word += 1;
        }
        if cur_word + 1 < words.len()
            && words[cur_word] == "if"
            && words[cur_word + 1].starts_with('b')
        {
            out.c_condition.0 = parse_u64_below_max(&words[cur_word + 1][1..]).unwrap();
            cur_word += 2;
        }

        if cur_word < words.len() && words[cur_word].starts_with('#') {
            // Ignore trailing comments
        } else if cur_word != words.len() {
            panic!("Failed to parse line '{}'", line);
        }

        out.validate();

        Some(out)
    }
}

pub struct Circuit {
    pub num_qubits: u64,
    pub num_bits: u64,
    pub num_registers: u64,
    pub operations: Vec<Op>,
    pub registers: Vec<Vec<QubitOrBit>>,
}

impl Circuit {
    pub fn from_text(text: &str) -> Self {
        let mut operations = Vec::new();
        for line in text.lines() {
            if let Some(op) = Op::from_text(line) {
                operations.push(op);
            }
        }
        let (num_qubits, num_bits, num_registers, registers) = analyze_ops(operations.iter());
        Self {
            num_qubits,
            num_bits,
            num_registers,
            operations,
            registers,
        }
    }

    pub fn from_kmx<P: AsRef<Path>>(path: P) -> std::io::Result<Self> {
        let file = File::open(path)?;
        let reader = BufReader::new(file);

        let mut operations = Vec::new();
        for line in reader.lines() {
            let line = line?;
            if let Some(op) = Op::from_text(&line) {
                operations.push(op);
            }
        }

        let (num_qubits, num_bits, num_registers, registers) = analyze_ops(operations.iter());
        Ok(Self {
            num_qubits,
            num_bits,
            num_registers,
            operations,
            registers,
        })
    }
}




#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct DepthStats {
    /// Longest chain of non-Clifford (CCX/CCZ) gates respecting data hazards —
    /// the fault-tolerance-relevant "T-depth" / magic-state depth that sets the
    /// reaction-limited runtime.
    pub toffoli_depth: u64,
    /// Longest chain of quantum gates (any qubit-touching op) respecting hazards.
    pub gate_depth: u64,
}

/// Critical-path depth of the op stream under read/write (RAW/WAR/WAW) hazards on
/// qubits and classical bits.
///
/// A wire's availability is the completion time of the last op that *wrote* it; a
/// read waits for the last write (RAW); a write waits for the last write (WAW)
/// and the last read (WAR). `c_condition` and every bit on the condition stack
/// are treated as READS — a controlled gate depends on the bit but does not
/// advance it — so gates that merely share a condition are not serialized behind
/// it. This keeps the deferred phase corrections of the measurement-based
/// uncompute in the layer their real dependencies allow.
///
/// `num_qubits`/`num_bits` must bound every id in the stream (use `analyze_ops`).
pub fn analyze_depth<'b>(
    ops: impl Iterator<Item = &'b Op>,
    num_qubits: usize,
    num_bits: usize,
) -> DepthStats {
    // last-write / last-read completion time per wire, one set per metric.
    let mut qwt = vec![0u64; num_qubits];
    let mut qrt = vec![0u64; num_qubits];
    let mut bwt = vec![0u64; num_bits];
    let mut brt = vec![0u64; num_bits];
    let mut qwg = vec![0u64; num_qubits];
    let mut qrg = vec![0u64; num_qubits];
    let mut bwg = vec![0u64; num_bits];
    let mut brg = vec![0u64; num_bits];

    let mut stack: Vec<u64> = Vec::new(); // condition-bit ids currently pushed
    let mut max_t = 0u64;
    let mut max_g = 0u64;

    for op in ops {
        match op.kind {
            OperationType::PushCondition => {
                if op.c_condition != NO_BIT {
                    stack.push(op.c_condition.0);
                }
                continue;
            }
            OperationType::PopCondition => {
                stack.pop();
                continue;
            }
            OperationType::Register
            | OperationType::AppendToRegister
            | OperationType::DebugPrint
            | OperationType::Neg => continue,
            _ => {}
        }

        // Read/write qubit sets by kind (validate() guarantees no intra-op aliasing).
        let mut rq = [0u64; 3];
        let mut nrq = 0usize;
        let mut wq = [0u64; 2];
        let mut nwq = 0usize;
        match op.kind {
            OperationType::CCX => {
                rq[0] = op.q_control1.0;
                rq[1] = op.q_control2.0;
                nrq = 2;
                wq[0] = op.q_target.0;
                nwq = 1;
            }
            OperationType::CX => {
                rq[0] = op.q_control1.0;
                nrq = 1;
                wq[0] = op.q_target.0;
                nwq = 1;
            }
            OperationType::Swap => {
                wq[0] = op.q_control1.0;
                wq[1] = op.q_target.0;
                nwq = 2;
            }
            OperationType::X | OperationType::Hmr | OperationType::R => {
                wq[0] = op.q_target.0;
                nwq = 1;
            }
            OperationType::Z => {
                rq[0] = op.q_target.0;
                nrq = 1;
            }
            OperationType::CZ => {
                rq[0] = op.q_target.0;
                rq[1] = op.q_control1.0;
                nrq = 2;
            }
            OperationType::CCZ => {
                rq[0] = op.q_target.0;
                rq[1] = op.q_control1.0;
                rq[2] = op.q_control2.0;
                nrq = 3;
            }
            _ => {} // BitInvert / BitStore*: classical-bit writes only
        }
        let wb: Option<u64> = if op.c_target != NO_BIT {
            Some(op.c_target.0)
        } else {
            None
        };
        let cost_t = matches!(op.kind, OperationType::CCX | OperationType::CCZ) as u64;
        let cost_g = matches!(
            op.kind,
            OperationType::CCX
                | OperationType::CX
                | OperationType::Swap
                | OperationType::X
                | OperationType::Z
                | OperationType::CZ
                | OperationType::CCZ
                | OperationType::Hmr
                | OperationType::R
        ) as u64;

        macro_rules! metric {
            ($qw:ident, $qr:ident, $bw:ident, $br:ident, $cost:expr, $max:ident) => {{
                let mut start = 0u64;
                for k in 0..nrq {
                    let i = rq[k] as usize;
                    if $qw[i] > start { start = $qw[i]; }
                }
                for k in 0..nwq {
                    let i = wq[k] as usize;
                    if $qw[i] > start { start = $qw[i]; }
                    if $qr[i] > start { start = $qr[i]; }
                }
                if op.c_condition != NO_BIT {
                    let i = op.c_condition.0 as usize;
                    if $bw[i] > start { start = $bw[i]; }
                }
                for &c in &stack {
                    let i = c as usize;
                    if $bw[i] > start { start = $bw[i]; }
                }
                if let Some(b) = wb {
                    let i = b as usize;
                    if $bw[i] > start { start = $bw[i]; }
                    if $br[i] > start { start = $br[i]; }
                }
                let comp = start + $cost;
                for k in 0..nrq {
                    let i = rq[k] as usize;
                    if comp > $qr[i] { $qr[i] = comp; }
                }
                for k in 0..nwq {
                    let i = wq[k] as usize;
                    $qw[i] = comp;
                }
                if op.c_condition != NO_BIT {
                    let i = op.c_condition.0 as usize;
                    if comp > $br[i] { $br[i] = comp; }
                }
                for &c in &stack {
                    let i = c as usize;
                    if comp > $br[i] { $br[i] = comp; }
                }
                if let Some(b) = wb {
                    let i = b as usize;
                    $bw[i] = comp;
                }
                if comp > $max { $max = comp; }
            }};
        }
        metric!(qwt, qrt, bwt, brt, cost_t, max_t);
        metric!(qwg, qrg, bwg, brg, cost_g, max_g);
    }

    DepthStats {
        toffoli_depth: max_t,
        gate_depth: max_g,
    }
}

pub fn analyze_ops<'b>(ops: impl Iterator<Item = &'b Op>) -> (u64, u64, u64, Vec<Vec<QubitOrBit>>) {
    let mut registers: Vec<Vec<QubitOrBit>> = Vec::new();
    let mut num_qubits = 0u64;
    let mut num_bits = 0u64;
    let mut num_registers = 0u64;

    for native_op in ops {
        if native_op.q_control2 != NO_QUBIT {
            num_qubits = num_qubits.max(native_op.q_control2.0 + 1);
        }
        if native_op.q_control1 != NO_QUBIT {
            num_qubits = num_qubits.max(native_op.q_control1.0 + 1);
        }
        if native_op.q_target != NO_QUBIT {
            num_qubits = num_qubits.max(native_op.q_target.0 + 1);
        }
        if native_op.c_target != NO_BIT {
            num_bits = num_bits.max(native_op.c_target.0 + 1);
        }
        if native_op.c_condition != NO_BIT {
            num_bits = num_bits.max(native_op.c_condition.0 + 1);
        }
        if native_op.r_target != NO_REG {
            num_registers = num_registers.max(native_op.r_target.0 + 1);
            while registers.len() <= native_op.r_target.0 as usize {
                registers.push(Vec::new());
            }
        }
        if native_op.kind == OperationType::AppendToRegister {
            if native_op.q_target != NO_QUBIT {
                registers[native_op.r_target.0 as usize].push(QubitOrBit::Qubit(native_op.q_target));
            }
            if native_op.c_target != NO_BIT {
                registers[native_op.r_target.0 as usize].push(QubitOrBit::Bit(native_op.c_target));
            }
        }
    }

    (num_qubits, num_bits, num_registers, registers)
}

#[cfg(test)]
mod depth_tests {
    use super::*;

    fn ccx(c1: u64, c2: u64, t: u64) -> Op {
        let mut o = Op::empty();
        o.kind = OperationType::CCX;
        o.q_control1 = QubitId(c1);
        o.q_control2 = QubitId(c2);
        o.q_target = QubitId(t);
        o
    }
    fn cx(c1: u64, t: u64) -> Op {
        let mut o = Op::empty();
        o.kind = OperationType::CX;
        o.q_control1 = QubitId(c1);
        o.q_target = QubitId(t);
        o
    }
    fn push(b: u64) -> Op {
        let mut o = Op::empty();
        o.kind = OperationType::PushCondition;
        o.c_condition = BitId(b);
        o
    }
    fn pop() -> Op {
        let mut o = Op::empty();
        o.kind = OperationType::PopCondition;
        o
    }

    #[test]
    fn dependent_toffolis_stack_in_depth() {
        // Two CCX writing the same target must serialize (WAW): depth 2.
        let ops = [ccx(0, 1, 2), ccx(0, 1, 2)];
        let d = analyze_depth(ops.iter(), 3, 0);
        assert_eq!(d.toffoli_depth, 2);
        assert_eq!(d.gate_depth, 2);
    }

    #[test]
    fn disjoint_toffolis_share_a_layer() {
        // Disjoint qubits => same layer => depth 1.
        let ops = [ccx(0, 1, 2), ccx(3, 4, 5)];
        let d = analyze_depth(ops.iter(), 6, 0);
        assert_eq!(d.toffoli_depth, 1);
        assert_eq!(d.gate_depth, 1);
    }

    #[test]
    fn critical_path_through_cnot() {
        // CCX(0,1,2); CCX(3,4,5); CX(2,3); CCX(0,1,2)
        // toffoli critical path: CCX#1 -> (CX carries 2->3 dep) -> CCX#4 via WAW on 2 = 2.
        // gate critical path: CCX#1 (reads/writes 2) -> CX reads 2 (layer2) -> CCX#4 writes 2
        //   must follow the CX read of 2 (WAR) => layer 3.
        let ops = [ccx(0, 1, 2), ccx(3, 4, 5), cx(2, 3), ccx(0, 1, 2)];
        let d = analyze_depth(ops.iter(), 6, 0);
        assert_eq!(d.toffoli_depth, 2);
        assert_eq!(d.gate_depth, 3);
    }

    #[test]
    fn shared_condition_does_not_serialize() {
        // Both CCX are controlled by bit 0 (a READ). Disjoint qubits => depth 1,
        // i.e. the condition bit does not force them into separate layers.
        let ops = [push(0), ccx(0, 1, 2), ccx(3, 4, 5), pop()];
        let d = analyze_depth(ops.iter(), 6, 1);
        assert_eq!(d.toffoli_depth, 1);
    }
}


