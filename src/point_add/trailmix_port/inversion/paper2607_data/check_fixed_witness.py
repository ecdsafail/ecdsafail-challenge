#!/usr/bin/env python3
"""Dependency-free checker for an Algorithm-3 fixed-schedule witness."""

from __future__ import annotations

import argparse
import hashlib
import json
from math import gcd
from pathlib import Path


P = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F


def euclid_word(numerator: int, denominator: int) -> list[int]:
    result = []
    while denominator:
        quotient, remainder = divmod(numerator, denominator)
        result.append(quotient)
        numerator, denominator = denominator, remainder
    return result


def continuant(word: list[int]) -> int:
    previous, current = 0, 1
    for quotient in word:
        assert quotient >= 1
        previous, current = current, quotient * current + previous
    return current


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "certificate",
        nargs="?",
        type=Path,
        default=Path(__file__).with_name("counterexample_cost_384.json"),
    )
    args = parser.parse_args()
    raw = args.certificate.read_bytes()
    record = json.loads(raw)
    assert record["schema"] == "secp256k1-algorithm3-schedule-counterexample-v1"
    assert int(record["p_hex"], 16) == P
    x = int(record["x_hex"], 16)
    assert int(record["x_decimal"]) == x
    assert 1 <= x <= P // 2 and gcd(P, x) == 1
    word = euclid_word(P, x)
    assert word == [int(q) for q in record["quotients"]]
    assert word[0] >= 2 and word[-1] >= 2
    assert continuant(word) == P
    cost = sum(q.bit_length() for q in word)
    assert cost == int(record["weighted_cost"]) == 384
    assert len(word) == int(record["quotient_count"])
    assert 4 * cost == int(record["algorithm3_steps"]) == 1536
    assert int(record["disproves_weighted_cost_at_most"]) == 383
    assert int(record["disproves_fixed_steps_at_most"]) == 1532
    print(
        f"PASS cost={cost} steps={4 * cost} quotients={len(word)} "
        f"sha256={hashlib.sha256(raw).hexdigest()}"
    )


if __name__ == "__main__":
    main()
