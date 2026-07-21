#!/usr/bin/env python3
"""Independent regression checks for the secp256k1 active-window table."""

from __future__ import annotations

import argparse
import hashlib
import json
import random
from pathlib import Path

import derive_active_windows as derive


def contains(window: list[int] | None, required: list[int]) -> bool:
    return window is not None and window[0] <= required[0] and window[1] >= required[1]


def verify_exact_input(rows: list[dict[str, object]], x: int) -> int:
    checked = 0
    for step, required_by_block in derive.exact_trace_requirements(x):
        if step > derive.SAFE_STEPS:
            raise AssertionError(f"x={hex(x)} exceeds certified schedule at step {step}")
        safe = rows[step - 1]["safe"]
        for block, required in required_by_block.items():
            if not contains(safe[block], required):
                raise AssertionError(
                    f"x={hex(x)} step={step} block={block} required={required} safe={safe[block]}"
                )
            checked += 1
    return checked


def main() -> None:
    here = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser()
    parser.add_argument("--table", type=Path, default=here / "active_windows_1616.json")
    parser.add_argument("--random-cases", type=int, default=10_000)
    args = parser.parse_args()

    encoded = args.table.read_bytes()
    table = json.loads(encoded)
    if table["schema"] != "luo-secp256k1-active-windows-v2":
        raise AssertionError("wrong table schema")
    if len(table["rows"]) != derive.SAFE_STEPS:
        raise AssertionError("wrong row count")
    if table["certificate"]["sha256"] != hashlib.sha256(
        (here / "certificate.json").read_bytes()
    ).hexdigest():
        raise AssertionError("certificate hash mismatch")

    rebuilt = derive.build_table(here / "certificate.json")
    if rebuilt != table:
        raise AssertionError("table is not the deterministic derivation output")

    for row in table["rows"]:
        if int(row["step"]) < 1 or int(row["step"]) > derive.SAFE_STEPS:
            raise AssertionError("bad step index")
        for block, window in row["safe"].items():
            if window is None:
                continue
            maximum = derive.N + 2 if block == "quotient_swap" else derive.WORK_SIZE
            if not (1 <= window[0] <= window[1] <= maximum):
                raise AssertionError(f"invalid safe window step={row['step']} block={block}: {window}")

    known = [
        1,
        2,
        3,
        derive.P // 2,
        int("5DB3D742C265539D92BA16B83C5C1DC492EC1A6629ED23CC63905323D8E62784", 16),
        int("5DB3D742C265539D92BA16B83C5C1DC492EC1A6629ED23CC63905323D96EFAEF", 16),
    ]
    rng = random.Random(0x260713816)
    inputs = known + [rng.randrange(1, derive.P // 2 + 1) for _ in range(args.random_cases)]
    checks = sum(verify_exact_input(table["rows"], x) for x in inputs)

    counterexamples = table["concrete_paper_counterexamples"]
    if len(counterexamples) != 6 or any(row["paper_contains_required"] for row in counterexamples):
        raise AssertionError("paper counterexample set changed")
    late_r = next(row for row in counterexamples if row["step"] == 1389 and row["block"] == "r_addsub")
    if late_r["one_lane_r_repair_contains_required"]:
        raise AssertionError("the one-lane R repair unexpectedly covers the exact late witness")

    shift = table["shift_register_requirement"]
    if shift["maximum_terminal_padding_steps"] != 592 or shift["required_counter_bits"] != 10:
        raise AssertionError("shift-register schedule requirement changed")

    print(f"table_sha256={hashlib.sha256(encoded).hexdigest()}")
    print(f"rows={len(table['rows'])}")
    print(f"exact_inputs={len(inputs)} exact_window_checks={checks}")
    print("paper_counterexamples=6")
    print("shift_counter_bits=10 max_padding=592")


if __name__ == "__main__":
    main()
