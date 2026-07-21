#!/usr/bin/env python3
"""Exact basis-state checks for the Q824 compact metadata/dirty12 route."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
UPSTREAM = Path("/private/tmp/paper2607-upstream")
MODEL = Path("/private/tmp/paper2607-model-agent")
sys.path[:0] = [str(HERE), str(UPSTREAM), str(MODEL), "/private/tmp"]

import eea_circuit_s835_exactwidth_dirty12 as eea
import test_eea_strict_main as test
from algorithm3_model import execute, preprocess, transition

P = 2**256 - 2**32 - 977
WITNESS_1500 = int(
    "5DB3D742C265539D92BA16B83C5C1DC492EC1A6629ED23CC63905323D8E62784", 16
)
WITNESS_1524 = int(
    "5DB3D742C265539D92BA16B83C5C1DC492EC1A6629ED23CC63905323D96EFAEF", 16
)


def apply_instruction(inst, global_qids, lanes, all_cases, *, inverse: bool) -> None:
    name = inst.name.lower()
    if name == "x":
        lanes[global_qids[0]] ^= all_cases
        return
    if name in {"cx", "cnot"}:
        lanes[global_qids[1]] ^= lanes[global_qids[0]]
        return
    if name in {"ccx", "tof", "toffoli"}:
        lanes[global_qids[2]] ^= lanes[global_qids[0]] & lanes[global_qids[1]]
        return
    if name in test.IGNORED:
        return
    definition = inst.definition
    if definition is None:
        raise ValueError(f"unsupported leaf {inst.name}")
    items = list(test._iter_items(definition))
    if inverse:
        items.reverse()
    for subinst, subqubits in items:
        sub_qids = [global_qids[test._qindex(definition, q)] for q in subqubits]
        apply_instruction(subinst, sub_qids, lanes, all_cases, inverse=inverse)


def apply_circuit(qc, lanes, all_cases, *, inverse: bool) -> None:
    items = list(test._iter_items(qc))
    if inverse:
        items.reverse()
    for inst, qargs in items:
        qids = [test._qindex(qc, q) for q in qargs]
        apply_instruction(inst, qids, lanes, all_cases, inverse=inverse)


def rotl(bits: str, amount: int) -> str:
    amount %= len(bits)
    return bits[amount:] + bits[:amount]


def work1_bits(state) -> str:
    t = (bin(state.t)[2:][::-1] if state.t else "") + "0"
    q = bin(state.q)[2:2 + state.l_q] if state.q else ""
    r = bin(state.r)[2:].zfill(state.p.bit_length() + 2 - state.l_t - state.l_q)
    out = t + q + r
    assert len(out) == state.p.bit_length() + 3, (state, out)
    return out


def work2_bits(state) -> str:
    n = state.p.bit_length()
    t = bin(state.t_prime)[2:].zfill(n + 3 - state.l_r_prime)[::-1]
    r = bin(state.r_prime)[2:] if state.r_prime else ""
    out = t + r
    assert len(out) == n + 3, (state, out)
    return rotl(out, state.l_shift)


def enc_lt(value: int) -> int:
    if not 1 <= value <= 256:
        raise AssertionError(f"l_t domain: {value}")
    return value - 1


def enc_lq(value: int) -> int:
    if not 0 <= value <= 256:
        raise AssertionError(f"l_q domain: {value}")
    return (value - 1) % (1 << eea.LQ_WIDTH)


def enc_lrp(value: int) -> int:
    if not 0 <= value <= 255:
        raise AssertionError(f"l_rp domain: {value}")
    return eea.LRP_ZERO if value == 0 else value - 1


def enc_ls(value: int) -> int:
    return (value - 1) % eea.LS_MODULUS


def reg_value(lanes, qc, name: str) -> int:
    reg = test._get_qreg(qc, name)
    return sum((lanes[test._qindex(qc, q)] & 1) << bit for bit, q in enumerate(reg))


def initialize(qc, state, dirty_seed: int) -> list[int]:
    initial: dict[int, int] = {}
    reg = lambda name: test._get_qreg(qc, name)
    for name, value in zip(("Phase1", "Phase2", "Iter", "Sign"), state.controls()):
        if value:
            initial[test._qindex(qc, reg(name)[0])] = 1
    test.set_bits_lr(initial, qc, reg("Work1"), work1_bits(state))
    test.set_bits_lr(initial, qc, reg("Work2"), work2_bits(state))
    for name, value in (
        ("l_t", enc_lt(state.l_t)),
        ("l_q", enc_lq(state.l_q)),
        ("l_s", enc_ls(state.l_shift)),
        ("l_rp", enc_lrp(state.l_r_prime)),
    ):
        test.set_reg_int_le(initial, qc, reg(name), value)
    for bit, qubit in enumerate(reg("DirtyPassenger")):
        if (dirty_seed >> bit) & 1:
            initial[test._qindex(qc, qubit)] = 1
    lanes = [0] * qc.num_qubits
    for qid, value in initial.items():
        lanes[qid] = value
    return lanes


def assert_state(label: str, qc, lanes, expected_state, dirty_before: int) -> None:
    reg = lambda name: test._get_qreg(qc, name)
    controls = tuple(lanes[test._qindex(qc, reg(name)[0])] & 1
                     for name in ("Phase1", "Phase2", "Iter", "Sign"))
    got_lengths = (
        reg_value(lanes, qc, "l_t"),
        reg_value(lanes, qc, "l_q"),
        reg_value(lanes, qc, "l_s"),
        reg_value(lanes, qc, "l_rp"),
    )
    want_lengths = (
        enc_lt(expected_state.l_t), enc_lq(expected_state.l_q),
        enc_ls(expected_state.l_shift), enc_lrp(expected_state.l_r_prime),
    )
    if controls != expected_state.controls() or got_lengths != want_lengths:
        raise AssertionError(
            f"{label}: controls={controls}/{expected_state.controls()} "
            f"lengths={got_lengths}/{want_lengths}"
        )
    if test.get_reg_bits_lr(lanes, qc, reg("Work1")) != work1_bits(expected_state):
        raise AssertionError(f"{label}: Work1 mismatch")
    if test.get_reg_bits_lr(lanes, qc, reg("Work2")) != work2_bits(expected_state):
        raise AssertionError(f"{label}: Work2 mismatch")
    if not test.clean_reg(lanes, qc, reg("Aux")):
        raise AssertionError(f"{label}: Aux not clean")
    dirty_after = reg_value(lanes, qc, "DirtyPassenger")
    if dirty_after != dirty_before:
        raise AssertionError(f"{label}: dirty changed {dirty_before:#x}->{dirty_after:#x}")


def check_step(label: str, x: int, step: int, dirty: int) -> None:
    before = execute(preprocess(P, x).state, step - 1)
    after = transition(before)
    qc = eea.build_step_circuit(
        256, step, T_max=1616, aux_size=eea.CLEAN_AUX_SIZE,
        measurement_uncompute=False,
    )
    if qc.num_qubits != 578:
        raise AssertionError(f"{label}: referenced width {qc.num_qubits}")
    lanes = initialize(qc, before, dirty)
    initial = lanes.copy()
    apply_circuit(qc, lanes, 1, inverse=False)
    assert_state(label, qc, lanes, after, dirty)
    apply_circuit(qc, lanes, 1, inverse=True)
    if lanes != initial:
        changed = [i for i, (a, b) in enumerate(zip(initial, lanes)) if a != b]
        raise AssertionError(f"{label}: reverse mismatch {changed[:12]}")
    print(f"PASS step={step} label={label} dirty={dirty:#x} reverse=exact")


def check_borrowed_c3x() -> None:
    from qiskit import QuantumCircuit, QuantumRegister
    q = QuantumRegister(6, "q")
    qc = QuantumCircuit(q)
    eea._borrowed_c3x(qc, q[0], q[1], q[2], q[3], q[4])
    for raw in range(1 << 5):
        lanes = [(raw >> bit) & 1 for bit in range(6)]
        initial = lanes.copy()
        apply_circuit(qc, lanes, 1, inverse=False)
        if lanes[3] != (initial[3] ^ (initial[0] & initial[1] & initial[2])):
            raise AssertionError("borrowed C3X target")
        if lanes[4] != initial[4]:
            raise AssertionError("borrowed C3X restoration")
        apply_circuit(qc, lanes, 1, inverse=True)
        if lanes != initial:
            raise AssertionError("borrowed C3X inverse")
    print("PASS borrowed_c3x states=32 dirty=restored phase=none")


def check_dirty_mcx_ladder() -> None:
    from qiskit import QuantumCircuit, QuantumRegister

    control_count = 9
    dirty_count = control_count - 2
    controls = QuantumRegister(control_count, "controls")
    target = QuantumRegister(1, "target")
    dirty = QuantumRegister(dirty_count, "dirty")
    qc = QuantumCircuit(controls, target, dirty)
    eea._mcx_dirty_ladder(qc, controls, target[0], dirty)
    counts = eea._e.count_circuit_ops_recursive(qc)
    if counts != {"ccx": 4 * control_count - 8}:
        raise AssertionError(f"dirty MCX primitive count: {counts}")

    width = qc.num_qubits
    states = 1 << width

    def lane_pattern(bit: int) -> int:
        span = 1 << bit
        period = span << 1
        repeats = ((1 << states) - 1) // ((1 << period) - 1)
        return repeats * ((1 << span) - 1) << span

    lanes = [lane_pattern(bit) for bit in range(width)]
    initial = lanes.copy()
    apply_circuit(qc, lanes, (1 << states) - 1, inverse=False)
    predicate = initial[0]
    for lane in initial[1:control_count]:
        predicate &= lane
    expected_target = initial[control_count] ^ predicate
    if lanes[control_count] != expected_target:
        raise AssertionError("dirty MCX target truth table")
    if lanes[:control_count] != initial[:control_count]:
        raise AssertionError("dirty MCX controls changed")
    if lanes[control_count + 1:] != initial[control_count + 1:]:
        raise AssertionError("dirty MCX lenders changed")
    apply_circuit(qc, lanes, (1 << states) - 1, inverse=True)
    if lanes != initial:
        raise AssertionError("dirty MCX inverse")
    print(
        f"PASS dirty_mcx controls={control_count} lenders={dirty_count} "
        f"states={states} ccx={counts['ccx']} reverse=exact phase=none"
    )


def check_mod259() -> None:
    from qiskit import QuantumCircuit, QuantumRegister
    ctrl = QuantumRegister(1, "ctrl")
    reg = QuantumRegister(eea.LS_WIDTH, "reg")
    scratch = QuantumRegister(eea.LS_WIDTH - 1, "scratch")
    inc = QuantumCircuit(ctrl, reg, scratch)
    eea.inc_mod259_1ctrl(inc, ctrl[0], reg, scratch)
    dec = QuantumCircuit(ctrl, reg, scratch)
    eea.dec_mod259_1ctrl(dec, ctrl[0], reg, scratch)
    for c in (0, 1):
        for value in range(1 << eea.LS_WIDTH):
            lanes = [c] + [(value >> bit) & 1 for bit in range(eea.LS_WIDTH)] + [0] * len(scratch)
            initial = lanes.copy()
            apply_circuit(inc, lanes, 1, inverse=False)
            got = sum(lanes[1 + bit] << bit for bit in range(eea.LS_WIDTH))
            expected = value if not c else ((value + 1) % eea.LS_MODULUS if value < eea.LS_MODULUS else None)
            if expected is not None and got != expected:
                raise AssertionError(f"mod259 inc {value}->{got}, expected {expected}")
            if any(lanes[1 + eea.LS_WIDTH:]):
                raise AssertionError("mod259 inc scratch")
            apply_circuit(dec, lanes, 1, inverse=False)
            if lanes != initial:
                raise AssertionError(f"mod259 inverse value={value} ctrl={c}")
    print("PASS mod259 valid=259 invalid=253 global_inverse=512x2")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--quick", action="store_true")
    args = parser.parse_args()
    check_borrowed_c3x()
    check_dirty_mcx_ladder()
    check_mod259()
    cases = [
        ("half-prime-lrp", P // 2, 8, 0x155),
        ("w1500-lrp", WITNESS_1500, 240, 0x2AA),
        ("w1500-r", WITNESS_1500, 1389, 0x3A5),
        ("w1500-swap", WITNESS_1500, 1470, 0x17C),
        ("w1524-lt", WITNESS_1524, 1472, 0x2D3),
        ("w1524-terminal", WITNESS_1524, 1524, 0x0F3),
        ("x1-pad258", 1, 1282, 0x3FF),
        ("x1-pad259", 1, 1283, 0x001),
        ("x1-pad260", 1, 1284, 0x2A6),
        ("x1-pad518", 1, 1542, 0x199),
        ("x1-pad519", 1, 1543, 0x266),
        ("x1-pad592", 1, 1616, 0x155),
    ]
    if args.quick:
        cases = [cases[0], cases[2], cases[5], cases[7], cases[-1]]
    for case in cases:
        check_step(*case)


if __name__ == "__main__":
    main()
