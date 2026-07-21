#!/usr/bin/env python3
"""Independent checker for the Luo Algorithm-3 secp256k1 schedule proof."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


SECP256K1_P = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F
EXACT_BOUNDARY_COST = 405


def k_value(word: list[int]) -> int:
    left, right = 0, 1
    for digit in word:
        assert digit >= 1
        left, right = right, digit * right + left
    return right


def quotient_word(numerator: int, denominator: int) -> list[int]:
    result: list[int] = []
    while denominator != 0:
        digit = numerator // denominator
        numerator, denominator = denominator, numerator - digit * denominator
        result.append(digit)
    return result


def minimum_for_length(length: int) -> int:
    if length == 1:
        return 2
    return k_value([2] + [1] * (length - 2) + [2])


def recompute_dp(
    top_cost: int,
) -> tuple[list[list[tuple[int, int]]], list[int | None]]:
    """Recompute all nondominated continuant states without importing prover."""

    states: list[list[tuple[int, int]]] = [[] for _ in range(top_cost + 1)]
    endpoint_minimum: list[int | None] = [None] * (top_cost + 1)
    quotient_width_cap = SECP256K1_P.bit_length()

    for total in range(2, top_cost + 1):
        candidates: list[tuple[int, int]] = []
        endpoint_candidates: list[int] = []

        if total <= quotient_width_cap:
            singleton = 2 ** (total - 1)
            candidates.append((1, singleton))
            endpoint_candidates.append(singleton)

        for final_width in range(1, min(quotient_width_cap, total - 2) + 1):
            smallest_digit = 2 ** (final_width - 1)
            for older, newer in states[total - final_width]:
                completed = smallest_digit * newer + older
                if total > EXACT_BOUNDARY_COST and completed > SECP256K1_P:
                    continue
                candidates.append((newer, completed))
                if final_width >= 2:
                    endpoint_candidates.append(completed)

        if total <= EXACT_BOUNDARY_COST:
            assert endpoint_candidates
        endpoint_minimum[total] = (
            min(endpoint_candidates) if endpoint_candidates else None
        )

        candidates.sort()
        nondominated: list[tuple[int, int]] = []
        smallest_seen_second: int | None = None
        previous_pair: tuple[int, int] | None = None
        for pair in candidates:
            if pair == previous_pair:
                continue
            previous_pair = pair
            if smallest_seen_second is None or pair[1] < smallest_seen_second:
                nondominated.append(pair)
                smallest_seen_second = pair[1]
        states[total] = nondominated
        if total == EXACT_BOUNDARY_COST:
            for earlier in range(2, total + 1):
                states[earlier] = [
                    pair for pair in states[earlier] if pair[1] <= SECP256K1_P
                ]

    return states, endpoint_minimum


def state_hash(states: list[list[tuple[int, int]]]) -> str:
    result = hashlib.sha256()
    for cost, row in enumerate(states):
        if not row:
            continue
        result.update((str(cost) + "\n").encode())
        for first, second in row:
            result.update((str(first) + "," + str(second) + "\n").encode())
    return result.hexdigest()


def check_minimizer(cost: int, record: dict[str, object], expected: int) -> None:
    widths = [int(value) for value in record["bit_lengths"]]
    digits = [int(value) for value in record["quotients"]]
    assert sum(widths) == cost
    assert widths[0] >= 2 and widths[-1] >= 2
    assert digits == [2 ** (width - 1) for width in widths]
    assert sum(digit.bit_length() for digit in digits) == cost
    assert k_value(digits) == expected == int(record["numerator"])


def check_witness(record: dict[str, object]) -> None:
    x = int(record["x_hex"], 16)
    x_used = min(x, SECP256K1_P - x)
    assert 1 <= x_used <= SECP256K1_P // 2
    assert int(record["x_used_hex"], 16) == x_used
    digits = quotient_word(SECP256K1_P, x_used)
    assert digits == [int(value) for value in record["quotients"]]
    assert digits[0] >= 2 and digits[-1] >= 2
    assert k_value(digits) == SECP256K1_P
    weighted = sum(digit.bit_length() for digit in digits)
    assert weighted == int(record["weighted_cost"])
    assert 4 * weighted == int(record["algorithm3_steps"])
    assert len(digits) == int(record["quotient_count"])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "certificate",
        nargs="?",
        type=Path,
        default=Path(__file__).with_name("certificate.json"),
    )
    args = parser.parse_args()
    certificate = json.loads(args.certificate.read_text(encoding="utf-8"))

    assert certificate["schema"] == "luo-algorithm3-fixed-schedule-bound-v1"
    assert int(certificate["p_decimal"]) == SECP256K1_P
    assert int(certificate["p_hex"], 16) == SECP256K1_P

    coarse = certificate["coarse_finite_cap"]
    length_cap = int(coarse["maximum_quotient_count"])
    assert minimum_for_length(length_cap) <= SECP256K1_P
    assert minimum_for_length(length_cap + 1) > SECP256K1_P
    assert int(coarse["minimum_numerator_at_cap"]) == minimum_for_length(length_cap)
    assert int(coarse["minimum_numerator_at_next_length"]) == minimum_for_length(
        length_cap + 1
    )

    log_product_cap = SECP256K1_P.bit_length() - 1
    assert int(coarse["sum_floor_log2_quotients_cap"]) == log_product_cap
    total_cap = length_cap + log_product_cap
    assert int(coarse["weighted_cost_cap"]) == total_cap

    states, minima = recompute_dp(total_cap)
    encoded_minima = [
        str(minima[cost]) for cost in range(2, EXACT_BOUNDARY_COST + 1)
    ]
    pareto = certificate["pareto_dp"]
    assert encoded_minima == pareto[
        "canonical_minima_decimal_through_first_excluded"
    ]
    assert [
        minima[cost] is not None and minima[cost] <= SECP256K1_P
        for cost in range(2, total_cap + 1)
    ] == pareto["canonical_feasible_through_cap"]
    assert [len(states[cost]) for cost in range(2, total_cap + 1)] == pareto[
        "frontier_sizes"
    ]
    assert state_hash(states) == pareto["frontier_sha256"]

    result = certificate["result"]
    bound = int(result["weighted_cost_upper_bound"])
    assert minima[bound] is not None and minima[bound] <= SECP256K1_P
    assert minima[bound + 1] is not None and minima[bound + 1] > SECP256K1_P
    assert all(minima[cost] is None for cost in range(bound + 2, total_cap + 1))
    assert int(result["first_excluded_weighted_cost"]) == bound + 1
    assert int(result["first_excluded_minimum_numerator"]) == minima[bound + 1]
    assert int(result["safe_fixed_schedule_steps"]) == 4 * bound
    assert int(result["all_higher_costs_checked_through"]) == total_cap

    for text_cost, record in pareto["minimizers"].items():
        cost = int(text_cost)
        assert minima[cost] is not None
        check_minimizer(cost, record, minima[cost])
    for witness in certificate["secp_witnesses"]:
        check_witness(witness)

    largest_witness = max(
        int(witness["algorithm3_steps"]) for witness in certificate["secp_witnesses"]
    )
    print(
        f"PASS p={hex(SECP256K1_P)} weighted_cost<={bound} "
        f"fixed_steps<={4 * bound} largest_witness={largest_witness}"
    )
    print(f"checked_costs=2..{total_cap} frontier_sha256={state_hash(states)}")


if __name__ == "__main__":
    main()
