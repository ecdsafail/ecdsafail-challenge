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
MIDTAIL_COUNTEREXAMPLE = int(
    "3388404F41DACAE921DD05E202DD17CCF8EEEB35849727E593938548D173053A", 16
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


def check_clean_c3x_mbu() -> None:
    from qiskit import QuantumCircuit, QuantumRegister

    q = QuantumRegister(5, "q")
    qc = QuantumCircuit(q)
    eea._clean_c3x_mbu(qc, q[0], q[1], q[2], q[3], q[4])
    counts = qc.count_ops()
    markers = sum(value for name, value in counts.items()
                  if name.lower() == "clean_c3x_mbu")
    if markers != 1:
        raise AssertionError(f"clean C3X marker count: {counts}")

    for raw in range(1 << 4):
        lanes = [(raw >> bit) & 1 for bit in range(4)] + [0]
        initial = lanes.copy()
        apply_circuit(qc, lanes, 1, inverse=False)
        if lanes[3] != (initial[3] ^ (initial[0] & initial[1] & initial[2])):
            raise AssertionError(f"clean C3X target raw={raw:#x}")
        if lanes[4] != 0:
            raise AssertionError(f"clean C3X temporary raw={raw:#x}")
        apply_circuit(qc, lanes, 1, inverse=True)
        if lanes != initial:
            raise AssertionError(f"clean C3X inverse raw={raw:#x}")
    print("PASS clean_c3x_mbu states=16 temp=zero reverse=exact phase=discharged")


def check_kg_equality_clean_hmr() -> None:
    from qiskit import QuantumCircuit, QuantumRegister

    q = QuantumRegister(5, "q")
    qc = QuantumCircuit(q)
    eea._kg_toggle_equality(
        qc, base=[q[0], q[1]], c0=q[2], flag=q[3], clean_temp=q[4],
    )
    counts = qc.count_ops()
    markers = sum(value for name, value in counts.items()
                  if name.lower() == "clean_c3x_mbu")
    if markers != 1:
        raise AssertionError(f"KG equality marker count: {counts}")

    for raw in range(1 << 4):
        lanes = [(raw >> bit) & 1 for bit in range(4)] + [0]
        initial = lanes.copy()
        apply_circuit(qc, lanes, 1, inverse=False)
        expected = initial[3] ^ (initial[0] & initial[1] & initial[2])
        if lanes[3] != expected or lanes[4] != 0:
            raise AssertionError(f"KG equality raw={raw:#x}: {lanes}")
        apply_circuit(qc, lanes, 1, inverse=True)
        if lanes != initial:
            raise AssertionError(f"KG equality inverse raw={raw:#x}")
    print("PASS kg_equality_clean_hmr states=16 temp=zero reverse=exact phase=discharged")


def check_r_fused_mode_cell() -> None:
    from qiskit import QuantumCircuit, QuantumRegister

    # mode, ctrl, addend, target, carry, arbitrary reference dirty, clean temp
    q = QuantumRegister(7, "q")
    fused = QuantumCircuit(q)
    eea._apply_r_fused_second_cell_clean_hmr(
        fused, mode=q[0], ctrl=q[1], addend=q[2], target=q[3],
        carry=q[4], clean_temp=q[6],
    )
    finish_sub = QuantumCircuit(q)
    eea._apply_cell_borrowed(
        finish_sub, "sub", "second", q[1], q[2], q[3], q[4], q[5],
    )
    undo_first = QuantumCircuit(q)
    eea._apply_cell_borrowed(
        undo_first, "add", "second", q[1], q[2], q[3], q[4], q[5],
    )
    counts = fused.count_ops()
    markers = sum(value for name, value in counts.items()
                  if name.lower() == "clean_c3x_mbu")
    if markers != 1 or counts.get("ccx", 0) != 4:
        raise AssertionError(f"fused R mode-cell primitive count: {counts}")

    for raw in range(1 << 6):
        initial = [(raw >> bit) & 1 for bit in range(6)] + [0]
        got = initial.copy()
        want = initial.copy()
        apply_circuit(fused, got, 1, inverse=False)
        apply_circuit(undo_first if initial[0] else finish_sub, want, 1, inverse=False)
        if got != want:
            raise AssertionError(f"fused R mode cell raw={raw:#x}: {got} != {want}")
        if got[5] != initial[5]:
            raise AssertionError(f"fused R mode cell dirty changed raw={raw:#x}")
        if got[6] != 0:
            raise AssertionError(f"fused R mode cell clean temp changed raw={raw:#x}")
        apply_circuit(fused, got, 1, inverse=True)
        if got != initial:
            raise AssertionError(f"fused R mode cell reverse raw={raw:#x}")
    print(
        "PASS r_fused_mode_cell states=64 executed_t=6 dirty=restored "
        "temp=zero reverse=exact phase=discharged"
    )


def check_r_fused_one_cell_equivalence() -> None:
    """Exhaust the old four-pass and fused two-pass maps on valid controls."""
    from qiskit import QuantumCircuit, QuantumRegister

    q = QuantumRegister(9, "q")
    ctrl, phase2, phase1, sign, addend, target, carry, dirty, clean = q

    old = QuantumCircuit(q)
    eea._apply_cell_borrowed(old, "sub", "first", ctrl, addend, target, carry, dirty)
    eea._apply_cell_borrowed(old, "sub", "second", ctrl, addend, target, carry, dirty)
    old.ccx(ctrl, phase2, sign)
    old.x(phase1)
    eea._borrowed_c3x(old, phase1, phase2, sign, ctrl, dirty)
    old.x(phase1)
    eea._apply_cell_borrowed(old, "add", "first", ctrl, addend, target, carry, dirty)
    eea._apply_cell_borrowed(old, "add", "second", ctrl, addend, target, carry, dirty)
    old.x(phase1)
    eea._borrowed_c3x(old, phase1, phase2, sign, ctrl, dirty)
    old.x(phase1)

    fused = QuantumCircuit(q)
    eea._apply_cell_clean_hmr(
        fused, "sub", "first", ctrl, addend, target, carry, clean,
    )
    fused.ccx(ctrl, phase2, sign)
    fused.x(phase1)
    fused.ccx(phase2, sign, phase1)
    eea._apply_r_fused_second_cell_clean_hmr(
        fused, mode=phase1, ctrl=ctrl, addend=addend,
        target=target, carry=carry, clean_temp=clean,
    )
    fused.ccx(phase2, sign, phase1)
    fused.x(phase1)

    tested = 0
    for raw in range(1 << 8):
        initial = [(raw >> bit) & 1 for bit in range(8)] + [0]
        c, p2, p1, s = initial[:4]
        valid_control = (
            (c == 1 and p1 == 0)
            or (c == 0 and p1 == 1)
            or (c == 0 and p1 == 0 and p2 == 0 and s == 0)
        )
        if not valid_control:
            continue
        got = initial.copy()
        want = initial.copy()
        apply_circuit(fused, got, 1, inverse=False)
        apply_circuit(old, want, 1, inverse=False)
        if got != want:
            raise AssertionError(
                f"fused R one-cell equivalence raw={raw:#x}: {got} != {want}"
            )
        if got[7] != initial[7]:
            raise AssertionError(f"fused R one-cell dirty changed raw={raw:#x}")
        if got[8] != 0:
            raise AssertionError(f"fused R one-cell clean temp changed raw={raw:#x}")
        apply_circuit(fused, got, 1, inverse=True)
        if got != initial:
            raise AssertionError(f"fused R one-cell reverse raw={raw:#x}")
        tested += 1
    print(
        f"PASS r_fused_one_cell_equivalence states={tested} "
        "old_four_pass=fused_two_pass dirty=restored temp=zero "
        "reverse=exact phase=discharged"
    )


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


def _basis_case_lanes(qc, cases: list[dict[str, int]]) -> tuple[list[int], int]:
    all_cases = (1 << len(cases)) - 1
    lanes = [0] * qc.num_qubits
    for case_index, values in enumerate(cases):
        case_mask = 1 << case_index
        for name, value in values.items():
            reg = test._get_qreg(qc, name)
            for bit, qubit in enumerate(reg):
                if (value >> bit) & 1:
                    lanes[test._qindex(qc, qubit)] |= case_mask
    return lanes, all_cases


def check_midtail_range_scan() -> None:
    from qiskit import QuantumCircuit, QuantumRegister

    cases = [
        {"Ctrl": ctrl, "Boundary": boundary}
        for ctrl in (0, 1)
        for boundary in range(259)
    ]
    for order in ("inc", "dec"):
        ctrl = QuantumRegister(1, "Ctrl")
        boundary = QuantumRegister(eea.LS_WIDTH, "Boundary")
        range_acc = QuantumRegister(1, "RangeAcc")
        path = QuantumRegister(8, "Path")
        output = QuantumRegister(259, "Output")
        qc = QuantumCircuit(ctrl, boundary, range_acc, path, output)

        def leaf(label, boundary_control, _clean_temp) -> None:
            qc.cx(boundary_control, output[label])

        eea._range_scan_259_nine(
            qc, boundary=boundary, ctrl=ctrl[0], range_acc=range_acc[0],
            path=path, leaf_fn=leaf, order=order,
        )
        lanes, all_cases = _basis_case_lanes(qc, cases)
        initial = lanes.copy()
        apply_circuit(qc, lanes, all_cases, inverse=False)
        for label, qubit in enumerate(output):
            expected = 0
            for case_index, case in enumerate(cases):
                if case["Ctrl"] and label <= case["Boundary"]:
                    expected |= 1 << case_index
            if lanes[test._qindex(qc, qubit)] != expected:
                raise AssertionError(f"midtail range {order} label={label}")
        if lanes[test._qindex(qc, range_acc[0])] != 0:
            raise AssertionError(f"midtail range {order} accumulator")
        if any(lanes[test._qindex(qc, qubit)] for qubit in path):
            raise AssertionError(f"midtail range {order} path")
        apply_circuit(qc, lanes, all_cases, inverse=True)
        if lanes != initial:
            raise AssertionError(f"midtail range {order} inverse")
    print("PASS midtail_range_scan valid_boundaries=259 controls=2 orders=2 phase=clean")


def check_midtail_upper_zero_map() -> None:
    from qiskit import QuantumCircuit, QuantumRegister

    ctrl = QuantumRegister(1, "Ctrl")
    boundary = QuantumRegister(eea.LS_WIDTH, "Boundary")
    bits = QuantumRegister(259, "Bits")
    dirty = QuantumRegister(259, "Dirty")
    scratch = QuantumRegister(9, "Scratch")
    qc = QuantumCircuit(ctrl, boundary, bits, dirty, scratch)
    eea._upper_zero_map_midpoint_nine(
        qc, ctrl=ctrl[0], boundary_B=boundary, bits=bits,
        dirty_map=dirty, scratch=scratch,
    )

    mask = (1 << 259) - 1
    boundaries = (0, 1, 127, 128, 255, 256, 257, 258)
    cases = []
    for case_index, boundary_value in enumerate(boundaries):
        for control in (0, 1):
            seed = (case_index + 1) * 0x9E3779B97F4A7C15 + control
            data = (seed * 0xD1342543DE82EF95) & mask
            data ^= ((1 << 259) - 1) // 3
            dirty_value = (seed * 0x94D049BB133111EB) & mask
            cases.append({
                "Ctrl": control,
                "Boundary": boundary_value,
                "Bits": data,
                "Dirty": dirty_value,
            })
    lanes, all_cases = _basis_case_lanes(qc, cases)
    initial = lanes.copy()
    apply_circuit(qc, lanes, all_cases, inverse=False)
    for case_index, case in enumerate(cases):
        suffix = 1
        expected = case["Dirty"]
        for label in range(258, -1, -1):
            data_bit = (case["Bits"] >> label) & 1
            in_range = case["Ctrl"] and label <= case["Boundary"]
            suffix &= 1 ^ (in_range & data_bit)
            if suffix:
                expected ^= 1 << label
        got = sum(
            ((lanes[test._qindex(qc, qubit)] >> case_index) & 1) << bit
            for bit, qubit in enumerate(dirty)
        )
        if got != expected:
            raise AssertionError(
                f"midtail upper map case={case_index} B={case['Boundary']}"
            )
    if any(lanes[test._qindex(qc, qubit)] for qubit in scratch):
        raise AssertionError("midtail upper map scratch")
    apply_circuit(qc, lanes, all_cases, inverse=False)
    if lanes != initial:
        raise AssertionError("midtail upper map involution")
    print(
        "PASS midtail_upper_zero_map cases=16 dirty=arbitrary "
        "scratch=clean involution=exact phase=clean"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--quick", action="store_true")
    args = parser.parse_args()
    check_borrowed_c3x()
    check_clean_c3x_mbu()
    check_kg_equality_clean_hmr()
    check_r_fused_mode_cell()
    check_r_fused_one_cell_equivalence()
    check_dirty_mcx_ladder()
    check_mod259()
    check_midtail_range_scan()
    check_midtail_upper_zero_map()
    cases = [
        ("midtail-counterexample", MIDTAIL_COUNTEREXAMPLE, 7, 0x155),
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
