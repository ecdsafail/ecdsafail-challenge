use crate::circuit::{Op, OperationType, NO_BIT, NO_QUBIT};
use std::collections::HashMap;

const NO_INDEX: usize = usize::MAX;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct FanoutWitness {
    pub(crate) first_index: usize,
    pub(crate) blocker_index: usize,
    pub(crate) second_index: usize,
    pub(crate) control_a: u64,
    pub(crate) control_b: u64,
    pub(crate) old_target: u64,
    pub(crate) new_target: u64,
    pub(crate) condition: u64,
}

#[derive(Clone, Copy, Debug, Hash, PartialEq, Eq)]
struct GateKey {
    control_a: u64,
    control_b: u64,
    target: u64,
}

#[derive(Clone, Copy, Debug)]
struct Candidate {
    index: usize,
    snapshot: [u64; 8],
}

struct Epochs {
    x_targets: Vec<u64>,
    x_controls: Vec<u64>,
    z_touches: Vec<u64>,
    hard_touches: Vec<u64>,
    swap_touches: Vec<u64>,
    swap_pairs: HashMap<(u64, u64), u64>,
    last_x_control: Vec<usize>,
}

impl Epochs {
    fn new(wire_count: usize) -> Self {
        Self {
            x_targets: vec![0; wire_count],
            x_controls: vec![0; wire_count],
            z_touches: vec![0; wire_count],
            hard_touches: vec![0; wire_count],
            swap_touches: vec![0; wire_count],
            swap_pairs: HashMap::new(),
            last_x_control: vec![NO_INDEX; wire_count],
        }
    }

    fn swap_pair(&self, a: u64, b: u64) -> u64 {
        *self.swap_pairs.get(&sorted_pair(a, b)).unwrap_or(&0)
    }
}

fn sorted_pair(a: u64, b: u64) -> (u64, u64) {
    if a <= b {
        (a, b)
    } else {
        (b, a)
    }
}

fn ccx_key(op: &Op) -> Option<GateKey> {
    if op.kind != OperationType::CCX || op.c_condition != NO_BIT {
        return None;
    }
    let (control_a, control_b) = sorted_pair(op.q_control1.0, op.q_control2.0);
    Some(GateKey {
        control_a,
        control_b,
        target: op.q_target.0,
    })
}

fn x_controls(op: &Op) -> Option<([u64; 2], usize)> {
    match op.kind {
        OperationType::X => Some(([NO_QUBIT.0; 2], 0)),
        OperationType::CX => Some(([op.q_control1.0, NO_QUBIT.0], 1)),
        OperationType::CCX => Some(([op.q_control1.0, op.q_control2.0], 2)),
        _ => None,
    }
}

fn quantum_support(op: &Op) -> ([u64; 3], usize) {
    match op.kind {
        OperationType::X | OperationType::Z | OperationType::R | OperationType::Hmr => {
            ([op.q_target.0, NO_QUBIT.0, NO_QUBIT.0], 1)
        }
        OperationType::CX | OperationType::CZ | OperationType::Swap => {
            ([op.q_control1.0, op.q_target.0, NO_QUBIT.0], 2)
        }
        OperationType::CCX | OperationType::CCZ => {
            ([op.q_control2.0, op.q_control1.0, op.q_target.0], 3)
        }
        _ => ([NO_QUBIT.0; 3], 0),
    }
}

fn max_wire(ops: &[Op]) -> usize {
    ops.iter()
        .flat_map(|op| {
            let (support, count) = quantum_support(op);
            support.into_iter().take(count)
        })
        .max()
        .unwrap_or(0) as usize
}

fn snapshot(key: GateKey, epochs: &Epochs) -> [u64; 8] {
    let swap_touches = epochs.swap_touches[key.control_a as usize]
        + epochs.swap_touches[key.control_b as usize]
        + epochs.swap_touches[key.target as usize];
    let swap_blockers = swap_touches - 2 * epochs.swap_pair(key.control_a, key.control_b);
    [
        epochs.x_targets[key.control_a as usize],
        epochs.x_targets[key.control_b as usize],
        epochs.x_controls[key.target as usize],
        epochs.z_touches[key.target as usize],
        epochs.hard_touches[key.control_a as usize],
        epochs.hard_touches[key.control_b as usize],
        epochs.hard_touches[key.target as usize],
        swap_blockers,
    ]
}

fn advance_epochs(op: &Op, index: usize, epochs: &mut Epochs) -> bool {
    if matches!(
        op.kind,
        OperationType::PushCondition | OperationType::PopCondition
    ) {
        return true;
    }
    if let Some((controls, count)) = x_controls(op) {
        epochs.x_targets[op.q_target.0 as usize] += 1;
        for &control in &controls[..count] {
            epochs.x_controls[control as usize] += 1;
            epochs.last_x_control[control as usize] = index;
        }
        return false;
    }
    match op.kind {
        OperationType::Z | OperationType::CZ | OperationType::CCZ => {
            let (support, count) = quantum_support(op);
            for &wire in &support[..count] {
                epochs.z_touches[wire as usize] += 1;
            }
        }
        OperationType::Swap => {
            let (a, b) = sorted_pair(op.q_control1.0, op.q_target.0);
            epochs.swap_touches[a as usize] += 1;
            epochs.swap_touches[b as usize] += 1;
            *epochs.swap_pairs.entry((a, b)).or_insert(0) += 1;
        }
        OperationType::R | OperationType::Hmr => {
            epochs.hard_touches[op.q_target.0 as usize] += 1;
        }
        _ => {}
    }
    false
}

fn validate_protected_tail(ops: &[Op], protected: usize) -> Result<Vec<Op>, String> {
    if protected > ops.len() || protected % 2 != 0 {
        return Err("invalid protected-tail length".to_owned());
    }
    let tail = &ops[ops.len() - protected..];
    for (pair_index, pair) in tail.chunks_exact(2).enumerate() {
        if pair[0] != pair[1]
            || pair[0].kind != OperationType::X
            || pair[0].c_condition != NO_BIT
        {
            return Err(format!(
                "protected nonce pair {pair_index} is not unconditional X/X"
            ));
        }
    }
    Ok(tail.to_vec())
}

pub(crate) fn rewrite_first_target_fanout(
    ops: Vec<Op>,
    protected_tail_ops: usize,
) -> Result<(Vec<Op>, FanoutWitness), String> {
    let protected_tail = validate_protected_tail(&ops, protected_tail_ops)?;
    let prefix_len = ops.len() - protected_tail_ops;
    let mut epochs = Epochs::new(max_wire(&ops) + 1);
    let mut candidates = HashMap::<GateKey, Candidate>::new();

    for index in 0..prefix_len {
        let op = ops[index];
        if let Some(key) = ccx_key(&op) {
            let current_snapshot = snapshot(key, &epochs);
            if let Some(prior) = candidates.get(&key).copied() {
                let mut deltas = [0u64; 8];
                let monotonic = deltas
                    .iter_mut()
                    .zip(current_snapshot.into_iter().zip(prior.snapshot))
                    .all(|(delta, (current, old))| {
                        if let Some(value) = current.checked_sub(old) {
                            *delta = value;
                            true
                        } else {
                            false
                        }
                    });
                let blocker_index = epochs.last_x_control[key.target as usize];
                let blocker = (blocker_index != NO_INDEX).then(|| ops[blocker_index]);
                if monotonic
                    && deltas == [0, 0, 1, 0, 0, 0, 0, 0]
                    && prior.index < blocker_index
                    && blocker_index < index
                    && blocker.is_some_and(|blocker| {
                        blocker.kind == OperationType::CX
                            && blocker.q_control1.0 == key.target
                            && blocker.q_target.0 != key.control_a
                            && blocker.q_target.0 != key.control_b
                            && blocker.q_target.0 != key.target
                    })
                {
                    let blocker = blocker.unwrap();
                    let mut replacement = Op::empty();
                    replacement.kind = OperationType::CCX;
                    replacement.q_control2.0 = key.control_a;
                    replacement.q_control1.0 = key.control_b;
                    replacement.q_target = blocker.q_target;
                    replacement.c_condition = blocker.c_condition;
                    let witness = FanoutWitness {
                        first_index: prior.index,
                        blocker_index,
                        second_index: index,
                        control_a: key.control_a,
                        control_b: key.control_b,
                        old_target: key.target,
                        new_target: blocker.q_target.0,
                        condition: blocker.c_condition.0,
                    };
                    let mut rewritten = Vec::with_capacity(ops.len() - 1);
                    for (op_index, stream_op) in ops.into_iter().enumerate() {
                        if op_index == prior.index || op_index == index {
                            continue;
                        }
                        rewritten.push(stream_op);
                        if op_index == blocker_index {
                            rewritten.push(replacement);
                        }
                    }
                    if rewritten.len() + 1 != prefix_len + protected_tail_ops {
                        return Err("single-fanout rewrite changed the wrong op count".to_owned());
                    }
                    if rewritten[rewritten.len() - protected_tail_ops..] != protected_tail {
                        return Err("single-fanout rewrite changed the nonce suffix".to_owned());
                    }
                    return Ok((rewritten, witness));
                }
            }
            candidates.insert(
                key,
                Candidate {
                    index,
                    snapshot: current_snapshot,
                },
            );
        }
        if advance_epochs(&op, index, &mut epochs) {
            candidates.clear();
        }
    }
    Err("no target-fanout conjugation found".to_owned())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::circuit::{BitId, QubitId};

    fn x(target: u64) -> Op {
        let mut op = Op::empty();
        op.kind = OperationType::X;
        op.q_target = QubitId(target);
        op
    }

    fn cx(control: u64, target: u64) -> Op {
        let mut op = Op::empty();
        op.kind = OperationType::CX;
        op.q_control1 = QubitId(control);
        op.q_target = QubitId(target);
        op
    }

    fn ccx(a: u64, b: u64, target: u64) -> Op {
        let mut op = Op::empty();
        op.kind = OperationType::CCX;
        op.q_control2 = QubitId(a);
        op.q_control1 = QubitId(b);
        op.q_target = QubitId(target);
        op
    }

    fn nonce_tail() -> Vec<Op> {
        (0..48).flat_map(|_| [x(0), x(0)]).collect()
    }

    fn eval(ops: &[Op], mut state: u8, condition: bool) -> u8 {
        for op in ops {
            if op.c_condition != NO_BIT && !condition {
                continue;
            }
            match op.kind {
                OperationType::CX => {
                    if ((state >> op.q_control1.0) & 1) != 0 {
                        state ^= 1 << op.q_target.0;
                    }
                }
                OperationType::CCX => {
                    if ((state >> op.q_control1.0) & 1) != 0
                        && ((state >> op.q_control2.0) & 1) != 0
                    {
                        state ^= 1 << op.q_target.0;
                    }
                }
                OperationType::X => state ^= 1 << op.q_target.0,
                _ => {}
            }
        }
        state
    }

    #[test]
    fn first_fanout_rewrite_is_exact_and_tail_stable() {
        let mut blocker = cx(2, 3);
        blocker.c_condition = BitId(7);
        let before_prefix = vec![ccx(0, 1, 2), blocker, ccx(1, 0, 2)];
        let mut before = before_prefix.clone();
        let tail = nonce_tail();
        before.extend(tail.clone());
        let (after, witness) = rewrite_first_target_fanout(before, 96).unwrap();
        assert_eq!(witness.first_index, 0);
        assert_eq!(witness.blocker_index, 1);
        assert_eq!(witness.second_index, 2);
        assert_eq!(witness.condition, 7);
        assert_eq!(&after[after.len() - 96..], tail.as_slice());
        for condition in [false, true] {
            for state in 0..16 {
                assert_eq!(
                    eval(&before_prefix, state, condition),
                    eval(&after[..2], state, condition)
                );
            }
        }
    }

    #[test]
    fn condition_stack_is_a_hard_barrier() {
        let mut push = Op::empty();
        push.kind = OperationType::PushCondition;
        push.c_condition = BitId(9);
        let mut pop = Op::empty();
        pop.kind = OperationType::PopCondition;
        let mut ops = vec![ccx(0, 1, 2), push, cx(2, 3), pop, ccx(0, 1, 2)];
        ops.extend(nonce_tail());
        assert!(rewrite_first_target_fanout(ops, 96).is_err());
    }
}
