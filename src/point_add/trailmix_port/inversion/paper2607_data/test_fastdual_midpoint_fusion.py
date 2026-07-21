#!/usr/bin/env python3
"""Exact permutation, cleanup, phase, and resource tests for midpoint fusion."""

from __future__ import annotations

from collections import Counter

from qiskit import QuantumCircuit, QuantumRegister

import eea_circuit_s835_fastdual_aux22 as eea
import test_eea_strict_main as strict
from test_fastdual_aux22_candidate import apply_circuit


def _qreg(qc: QuantumCircuit, name: str):
    return strict._get_qreg(qc, name)


def _set_int(lanes: list[int], qc: QuantumCircuit, reg, value: int, case: int) -> None:
    marker = 1 << case
    for bit, qubit in enumerate(reg):
        if (value >> bit) & 1:
            lanes[strict._qindex(qc, qubit)] |= marker


def _build_pair(n: int, *, measurement_uncompute: bool) -> QuantumCircuit:
    cfg = eea.get_n_config(n)
    len_width = int(cfg["len_width"])
    shift_width = eea.fixed_schedule_shift_width(
        n, int(cfg["shift_width"]), 4 * n + 8
    )
    work_size = n + 3
    eea.set_measurement_uncompute(measurement_uncompute)
    sub = eea.t_prefix_addsub_midpoint_gate(
        n=n, k=1, K=work_size, len_width=len_width,
        shift_width=shift_width, mode="sub", name=f"MID_SUB_N{n}",
    )
    add = eea.t_prefix_addsub_midpoint_gate(
        n=n, k=1, K=work_size, len_width=len_width,
        shift_width=shift_width, mode="add", name=f"MID_ADD_N{n}",
    )
    fixed = 4 + 2 * work_size + 2 * len_width + shift_width
    scratch_size = max(sub.num_qubits - fixed, add.num_qubits - fixed)

    ArithCtrl = QuantumRegister(1, "ArithCtrl")
    TailCtrl = QuantumRegister(1, "TailCtrl")
    Sign = QuantumRegister(1, "Sign")
    Tail = QuantumRegister(1, "Tail")
    Neg = QuantumRegister(1, "Neg")
    Work1 = QuantumRegister(work_size, "Work1")
    Work2 = QuantumRegister(work_size, "Work2")
    l_t = QuantumRegister(len_width, "l_t")
    l_s = QuantumRegister(shift_width, "l_s")
    l_rp = QuantumRegister(len_width, "l_rp")
    Scratch = QuantumRegister(scratch_size, "Scratch")
    qc = QuantumCircuit(
        ArithCtrl, TailCtrl, Sign, Tail, Neg, Work1, Work2,
        l_t, l_s, l_rp, Scratch, name=f"MIDPOINT_PAIR_N{n}",
    )
    args = (
        [ArithCtrl[0], TailCtrl[0], Tail[0], Neg[0]]
        + list(Work1) + list(Work2) + list(l_t) + list(l_s) + list(l_rp)
        + list(Scratch)
    )
    eea._e._append_with_optional_clbits(qc, sub, args[:sub.num_qubits])
    qc.cx(TailCtrl[0], Sign[0])
    qc.ccx(Tail[0], Neg[0], Sign[0])
    eea._e._append_with_optional_clbits(qc, add, args[:add.num_qubits])
    return qc


def _state_int(lanes: list[int], case: int) -> int:
    return sum(((lane >> case) & 1) << qid for qid, lane in enumerate(lanes))


def verify_small_permutations() -> None:
    total_cases = 0
    for n in range(3, 9):
        qc = _build_pair(n, measurement_uncompute=False)
        strict.assert_toffoli_network(qc, f"midpoint pair n={n}")
        len_width = len(_qreg(qc, "l_t"))
        lanes = [0] * qc.num_qubits
        expected_sign_toggle = 0
        combinations: set[tuple[int, int]] = set()
        case = 0

        # Exhaust the Boolean borrow/tail predicate space at every legal A.
        # Propagating all-one addends exercise every lower carry-chain length;
        # every tail bit pattern and both dirty patterns are used.
        for A in range(1, n + 1):
            for borrow in (0, 1):
                addend = (1 << A) - 1
                target = 0 if borrow else addend
                for tail_mask in range(1 << (n - A + 1)):
                    tail_zero = int(tail_mask == 0)
                    for sign in (0, 1):
                        for dirty_pattern in (0, 1):
                            marker = 1 << case
                            lanes[strict._qindex(qc, _qreg(qc, "ArithCtrl")[0])] |= marker
                            lanes[strict._qindex(qc, _qreg(qc, "TailCtrl")[0])] |= marker
                            if sign:
                                lanes[strict._qindex(qc, _qreg(qc, "Sign")[0])] |= marker
                            work1 = addend
                            if dirty_pattern:
                                work1 |= ((1 << (n + 3 - A)) - 1) << A
                            work2 = target | (tail_mask << A)
                            if dirty_pattern:
                                work2 |= 0b11 << (n + 1)
                            _set_int(lanes, qc, _qreg(qc, "Work1"), work1, case)
                            _set_int(lanes, qc, _qreg(qc, "Work2"), work2, case)
                            _set_int(
                                lanes, qc, _qreg(qc, "l_t"),
                                strict.enc_len(A - 1, len_width), case,
                            )
                            # Actual l_s=l_r'=1 gives B=n exactly.
                            _set_int(lanes, qc, _qreg(qc, "l_s"), 0, case)
                            _set_int(lanes, qc, _qreg(qc, "l_rp"), 0, case)
                            if 1 ^ (tail_zero & borrow):
                                expected_sign_toggle |= marker
                            combinations.add((borrow, tail_zero))
                            case += 1

        all_cases = (1 << case) - 1
        initial = lanes.copy()
        initial_states = {_state_int(initial, i) for i in range(case)}
        if len(initial_states) != case:
            raise AssertionError(f"n={n}: duplicate test-domain inputs")
        apply_circuit(qc, lanes, all_cases, inverse=False)
        expected = initial.copy()
        sign_qid = strict._qindex(qc, _qreg(qc, "Sign")[0])
        expected[sign_qid] ^= expected_sign_toggle
        if lanes != expected:
            bad = [qid for qid, (got, want) in enumerate(zip(lanes, expected)) if got != want]
            raise AssertionError(f"n={n}: forward mismatch lanes={bad[:12]}")
        output_states = {_state_int(lanes, i) for i in range(case)}
        if len(output_states) != case:
            raise AssertionError(f"n={n}: permutation collision on promised domain")
        apply_circuit(qc, lanes, all_cases, inverse=True)
        if lanes != initial:
            bad = [qid for qid, (got, want) in enumerate(zip(lanes, initial)) if got != want]
            raise AssertionError(f"n={n}: inverse/ancilla cleanup lanes={bad[:12]}")
        if combinations != {(0, 0), (0, 1), (1, 0), (1, 1)}:
            raise AssertionError(f"n={n}: incomplete predicate combinations {combinations}")
        total_cases += case
        print(
            f"PASS midpoint permutation n={n} cases={case} "
            "borrow_tail=4 inverse_clean=true ancilla_clean=true"
        )
    print(f"PASS midpoint exhaustive total_cases={total_cases}")


def verify_phase_cleanup() -> None:
    # The counted KMX lowering recognizes CLEAN_C3X_MBU and emits
    #   tmp ^= a*b; target ^= c*tmp; clear_and(tmp,a,b).
    # For every H-basis measurement outcome m, demolition contributes phase
    # m*a*b and clear_and's classically conditioned CZ contributes the same
    # phase, so the net debt is identically zero.  The reset returns tmp to 0.
    macro = eea.clean_c3x_mbu_gate()
    macro_names = [item.operation.name for item in macro.definition.data]
    if macro_names != ["ccx", "ccx", "ccx"]:
        raise AssertionError(f"unexpected CLEAN_C3X_MBU definition {macro_names}")
    phase_cases = 0
    for a in (0, 1):
        for b in (0, 1):
            for c in (0, 1):
                for target in (0, 1):
                    for measured in (0, 1):
                        tmp = a & b
                        target_after = target ^ (c & tmp)
                        phase_debt = (measured & tmp) ^ (measured & a & b)
                        if target_after != (target ^ (a & b & c)) or phase_debt:
                            raise AssertionError(
                                f"clean-C3 contract failure a={a} b={b} c={c} "
                                f"target={target} measured={measured}"
                            )
                        phase_cases += 1

    for n in (3, 4, 5, 6, 7, 8, 256):
        qc = _build_pair(n, measurement_uncompute=False)
        primitive_ops = eea._e.PRIMITIVE_OPS
        added = "CLEAN_C3X_MBU" not in primitive_ops
        primitive_ops.add("CLEAN_C3X_MBU")
        try:
            counts = Counter(eea._e.count_circuit_ops_recursive(qc))
        finally:
            if added:
                primitive_ops.remove("CLEAN_C3X_MBU")
        unexpected = set(counts) - {"x", "cx", "ccx", "CLEAN_C3X_MBU"}
        if unexpected:
            raise AssertionError(f"n={n}: phase-bearing leaves {sorted(unexpected)}")
        print(
            f"PASS phase cleanup n={n} clean_c3={counts['CLEAN_C3X_MBU']} "
            "unitary_leaves=x/cx/ccx phase_debt=0"
        )
    print(
        f"PASS CLEAN_C3X_MBU branch-phase contract cases={phase_cases} "
        "reset_clean=true phase_debt=0"
    )


def _executed_toffoli(obj) -> int:
    primitive_ops = eea._e.PRIMITIVE_OPS
    added = "CLEAN_C3X_MBU" not in primitive_ops
    primitive_ops.add("CLEAN_C3X_MBU")
    try:
        circuit = obj if hasattr(obj, "data") else obj.definition
        counts = Counter(eea._e.count_circuit_ops_recursive(circuit))
    finally:
        if added:
            primitive_ops.remove("CLEAN_C3X_MBU")
    return counts["ccx"] + 2 * counts["CLEAN_C3X_MBU"]


def verify_per_step_delta() -> None:
    eea.set_measurement_uncompute(False)
    n, len_width, shift_width = 256, 9, 10
    tail = _executed_toffoli(
        eea.t_tail_zero_toggle_gate(n=n, len_width=len_width, shift_width=shift_width)
    )
    borrow = _executed_toffoli(
        eea.t_lower_borrow_toggle_gate(n=n, len_width=len_width)
    )
    if (tail, borrow) != (6788, 1319):
        raise AssertionError(f"reference costs changed: tail={tail}, borrow={borrow}")
    for M in (1, 2, 5, 17, 65, 129, 259):
        k, K = 1, M
        old_sub = eea.lc_prefix_addsub_unary_gate(
            k=k, K=K, len_width=len_width, mode="sub", sign_update=False,
            target="work2", name=f"OLD_SUB_M{M}",
        )
        old_add = eea.lc_prefix_addsub_unary_gate(
            k=k, K=K, len_width=len_width, mode="add", sign_update=False,
            target="work2", name=f"OLD_ADD_M{M}",
        )
        new_sub = eea.t_prefix_addsub_midpoint_gate(
            n=n, k=k, K=K, len_width=len_width, shift_width=shift_width,
            mode="sub", name=f"NEW_SUB_M{M}",
        )
        new_add = eea.t_prefix_addsub_midpoint_gate(
            n=n, k=k, K=K, len_width=len_width, shift_width=shift_width,
            mode="add", name=f"NEW_ADD_M{M}",
        )
        old = 2 * tail + 2 * borrow + _executed_toffoli(old_sub) + _executed_toffoli(old_add)
        new = _executed_toffoli(new_sub) + _executed_toffoli(new_add) + 1
        expected = 2685 - 2 * M
        if old - new != expected:
            raise AssertionError(
                f"M={M}: old={old}, new={new}, delta={old-new}, expected={expected}"
            )
        print(f"PASS per-step delta M={M} old={old} new={new} saved={expected}")


def main() -> None:
    verify_small_permutations()
    verify_phase_cleanup()
    verify_per_step_delta()


if __name__ == "__main__":
    main()
