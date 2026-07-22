#!/usr/bin/env python3
"""Validate the independent bit-sliced serialized-stream probe output."""

from __future__ import annotations

import argparse
from pathlib import Path
import re


P = 2**256 - 2**32 - 977
SCHEDULE_STEPS = 1616
TRACE = re.compile(
    r"x=([0-9a-f]+) iter=([01]) padding=(\d+) work2=([0-9a-f]+)$"
)


def exact_steps(x: int) -> int:
    x = min(x, P - x)
    cost = 0
    previous = P
    while x:
        quotient, remainder = divmod(previous, x)
        cost += quotient.bit_length()
        previous, x = x, remainder
    return 4 * cost


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path)
    args = parser.parse_args()

    traces = []
    reverse_pass = False
    for raw_line in args.output.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        match = TRACE.fullmatch(line)
        if match:
            traces.append(match.groups())
        elif line == "PASS cases=9 forward_reverse=exact":
            reverse_pass = True
        elif line:
            raise AssertionError(f"unexpected probe output: {line}")

    if len(traces) != 9 or not reverse_pass:
        raise AssertionError(
            f"incomplete probe: traces={len(traces)} reverse_pass={reverse_pass}"
        )
    for x_hex, iteration_text, padding_text, work2_hex in traces:
        x = int(x_hex, 16)
        iteration = int(iteration_text)
        padding = int(padding_text)
        work2 = int(work2_hex, 16)
        if work2 >> 256:
            raise AssertionError(f"x={x_hex}: nonzero terminal padding lanes")
        inverse = work2 if iteration else (P - work2) % P
        expected_inverse = pow(x, -1, P)
        if inverse != expected_inverse:
            raise AssertionError(
                f"x={x_hex}: inverse={hex(inverse)} expected={hex(expected_inverse)}"
            )
        expected_padding = SCHEDULE_STEPS - exact_steps(x)
        if padding != expected_padding:
            raise AssertionError(
                f"x={x_hex}: padding={padding} expected={expected_padding}"
            )

    print(f"PASS traces={len(traces)} inverses=exact padding=exact reverse=exact")


if __name__ == "__main__":
    main()
