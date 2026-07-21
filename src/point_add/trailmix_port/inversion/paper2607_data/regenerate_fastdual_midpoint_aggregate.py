#!/usr/bin/env python3
"""Regenerate and prove the exact Aux=22 midpoint-fusion aggregate count."""

from __future__ import annotations

from collections import Counter
import hashlib
import json
from pathlib import Path

import eea_circuit_s835_fastdual_aux22 as eea


OLD_EXECUTED_T_PER_TRAVERSAL = 59_599_489
OLD_FOUR_TRAVERSAL_T = 238_397_956
OLD_FULL_PA_T = 249_076_062
EXPECTED_STEPS = 1616
EXPECTED_SUM_M = 282_924


def executed_toffoli(obj) -> int:
    primitive_ops = eea._e.PRIMITIVE_OPS
    added_macro = "CLEAN_C3X_MBU" not in primitive_ops
    primitive_ops.add("CLEAN_C3X_MBU")
    try:
        circuit = obj if hasattr(obj, "data") else obj.definition
        counts = Counter(eea._e.count_circuit_ops_recursive(circuit))
    finally:
        if added_macro:
            primitive_ops.remove("CLEAN_C3X_MBU")
    return counts["ccx"] + 2 * counts["CLEAN_C3X_MBU"]


def main() -> None:
    eea.set_measurement_uncompute(False)
    primitive_ops = eea._e.PRIMITIVE_OPS
    added_macro = "CLEAN_C3X_MBU" not in primitive_ops
    primitive_ops.add("CLEAN_C3X_MBU")
    counts: Counter = Counter()
    per_step_t: list[int] = []
    try:
        for step in range(1, EXPECTED_STEPS + 1):
            circuit = eea.build_step_circuit(
                256, step, T_max=EXPECTED_STEPS, aux_size=22,
                measurement_uncompute=False,
            )
            if circuit.num_qubits != 581:
                raise AssertionError(
                    f"step {step}: width={circuit.num_qubits}, expected 581"
                )
            step_counts = Counter(eea._e.count_circuit_ops_recursive(circuit))
            counts.update(step_counts)
            per_step_t.append(
                step_counts["ccx"] + 2 * step_counts["CLEAN_C3X_MBU"]
            )
    finally:
        if added_macro:
            primitive_ops.remove("CLEAN_C3X_MBU")

    rows = eea._certified_window_table["rows"]
    windows = []
    for row in rows:
        window = row["safe"]["t_addsub"]
        windows.append((1, 1) if window is None else (int(window[0]), int(window[1])))
    widths = [K - k + 1 for k, K in windows]
    sum_m = sum(widths)
    if sum_m != EXPECTED_SUM_M:
        raise AssertionError(f"sum M={sum_m}, expected {EXPECTED_SUM_M}")
    expected_savings = 2685 * EXPECTED_STEPS - 2 * sum_m
    tail_t = executed_toffoli(
        eea.t_tail_zero_toggle_gate(n=256, len_width=9, shift_width=10)
    )
    borrow_t = executed_toffoli(
        eea.t_lower_borrow_toggle_gate(n=256, len_width=9)
    )
    exact_delta_by_window = {}
    for k, K in sorted(set(windows)):
        M = K - k + 1
        old_sub = eea.lc_prefix_addsub_unary_gate(
            k=k, K=K, len_width=9, mode="sub", sign_update=False,
            target="work2", name=f"COUNT_OLD_SUB_{k}_{K}",
        )
        old_add = eea.lc_prefix_addsub_unary_gate(
            k=k, K=K, len_width=9, mode="add", sign_update=False,
            target="work2", name=f"COUNT_OLD_ADD_{k}_{K}",
        )
        new_sub = eea.t_prefix_addsub_midpoint_gate(
            n=256, k=k, K=K, len_width=9, shift_width=10,
            mode="sub", name=f"COUNT_NEW_SUB_{k}_{K}",
        )
        new_add = eea.t_prefix_addsub_midpoint_gate(
            n=256, k=k, K=K, len_width=9, shift_width=10,
            mode="add", name=f"COUNT_NEW_ADD_{k}_{K}",
        )
        old_component = (
            2 * tail_t + 2 * borrow_t
            + executed_toffoli(old_sub) + executed_toffoli(old_add)
        )
        new_component = executed_toffoli(new_sub) + executed_toffoli(new_add) + 1
        delta = old_component - new_component
        if delta != 2685 - 2 * M:
            raise AssertionError(
                f"window {k}..{K}: delta={delta}, expected {2685 - 2*M}"
            )
        exact_delta_by_window[f"{k}:{K}"] = delta
    per_step_deltas = [exact_delta_by_window[f"{k}:{K}"] for k, K in windows]
    if sum(per_step_deltas) != expected_savings:
        raise AssertionError(
            f"per-step delta sum={sum(per_step_deltas)}, expected {expected_savings}"
        )
    expected_new_t = OLD_EXECUTED_T_PER_TRAVERSAL - expected_savings
    executed_t = counts["ccx"] + 2 * counts["CLEAN_C3X_MBU"]
    if executed_t != expected_new_t:
        raise AssertionError(
            f"aggregate T={executed_t}, expected {expected_new_t}; "
            f"delta={OLD_EXECUTED_T_PER_TRAVERSAL-executed_t}"
        )

    records = sum(counts.values())
    emitted_ops = records + 3 * counts["CLEAN_C3X_MBU"]
    source = Path(eea.__file__).resolve()
    certificate = source.with_name("active_windows_1616.json")
    result = {
        "schema": "paper2607-eea-midpoint-fusion-count-v1",
        "source_module": source.name,
        "source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
        "window_certificate_sha256": hashlib.sha256(certificate.read_bytes()).hexdigest(),
        "field_width": 256,
        "schedule_steps": EXPECTED_STEPS,
        "aux_size": 22,
        "local_width": 581,
        "integrated_pa_qubits": 837,
        "sum_t_window_widths": sum_m,
        "distinct_t_windows": len(exact_delta_by_window),
        "per_step_delta": "2685 - 2*M",
        "per_step_delta_min": min(per_step_deltas),
        "per_step_delta_max": max(per_step_deltas),
        "old_executed_toffoli_per_traversal": OLD_EXECUTED_T_PER_TRAVERSAL,
        "new_executed_toffoli_per_traversal": executed_t,
        "executed_toffoli_saved_per_traversal": expected_savings,
        "old_four_traversal_executed_toffoli": OLD_FOUR_TRAVERSAL_T,
        "new_four_traversal_executed_toffoli": 4 * executed_t,
        "old_full_pa_toffoli": OLD_FULL_PA_T,
        "new_full_pa_toffoli": OLD_FULL_PA_T - 4 * expected_savings,
        "records_per_traversal": records,
        "emitted_ops_per_traversal": emitted_ops,
        "primitive_counts": {
            "ccx": counts["ccx"],
            "clean_c3x_mbu": counts["CLEAN_C3X_MBU"],
            "cx": counts["cx"],
            "x": counts["x"],
        },
        "per_step_executed_toffoli_min": min(per_step_t),
        "per_step_executed_toffoli_max": max(per_step_t),
        "assertions": {
            "exact_width": True,
            "exact_per_step_delta": True,
            "aggregate_matches_delta_sum": True,
        },
    }
    output = source.with_name("fastdual_aux22_midpoint_aggregate.json")
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, sort_keys=True))
    print(f"wrote {output}")


if __name__ == "__main__":
    main()
