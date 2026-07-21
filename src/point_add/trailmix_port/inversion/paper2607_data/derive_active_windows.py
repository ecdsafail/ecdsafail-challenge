#!/usr/bin/env python3
"""Derive conservative per-step windows for repaired Luo Algorithm 3.

The derivation is intentionally an over-approximation.  It uses the certified
weighted quotient-cost bound, continuant lower bounds, and the exact four-phase
layout of Algorithm 3.  No sampled input is used to construct a window.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path


P = 2**256 - 2**32 - 977
N = 256
WORK_SIZE = N + 3
MAX_WEIGHTED_COST = 404
SAFE_STEPS = 4 * MAX_WEIGHTED_COST
MAX_QUOTIENT_WEIGHT = P.bit_length()
C_EEA = 1.0 / math.log2((math.sqrt(5.0) + 1.0) / 2.0)


def fibonacci_table(limit: int) -> list[int]:
    values = [0, 1]
    while len(values) <= limit:
        values.append(values[-1] + values[-2])
    return values


FIBONACCI = fibonacci_table(MAX_WEIGHTED_COST + 3)


def prefix_continuant_lower(cost: int) -> int | None:
    """Lower-bound a prefix continuant with weighted cost ``cost``.

    A nonempty Euclidean prefix has quotient weights w_0 >= 2 and w_i >= 1,
    sum(w_i)=cost.  For a prefix of m quotients,

      K >= product(q_i) >= 2**(cost-m)
      K >= continuant(2,1,...,1) = Fibonacci[m+2].

    Taking the minimum over every possible m preserves a universal lower
    bound.  Cost zero denotes the initial coefficient 1; cost one is
    impossible because the first quotient is at least two.
    """
    if cost == 0:
        return 1
    if cost < 2:
        return None
    return min(
        max(1 << (cost - count), FIBONACCI[count + 2])
        for count in range(1, cost)
    )


PREFIX_LOWER = [prefix_continuant_lower(cost) for cost in range(MAX_WEIGHTED_COST + 1)]


def prefix_length_bounds(cost: int, current_weight: int) -> tuple[int, int] | None:
    """Bound bit_length(t) before a current quotient of given weight.

    The lower edge comes from ``prefix_continuant_lower``.  The upper edge uses
    both K_prefix < 2**cost and q*K_prefix < p for the current quotient.
    Returning None proves that this relaxed prefix/current pair is impossible.
    """
    lower_value = PREFIX_LOWER[cost]
    if lower_value is None:
        return None
    q_min = 1 << (current_weight - 1)
    if q_min * lower_value >= P:
        return None
    if cost == 0:
        return (1, 1)
    upper_value = min((1 << cost) - 1, (P - 1) // q_min)
    if upper_value < lower_value:
        return None
    return lower_value.bit_length(), upper_value.bit_length()


def feasible_weight_interval(prefix_cost: int, local_step: int) -> tuple[int, int] | None:
    """Return the contiguous relaxed interval of possible current weights."""
    lower_value = PREFIX_LOWER[prefix_cost]
    if lower_value is None:
        return None
    minimum_weight = max(1, (local_step + 3) // 4)
    if prefix_cost == 0:
        minimum_weight = max(2, minimum_weight)
    quotient_room = (P - 1) // lower_value
    maximum_by_product = quotient_room.bit_length()
    maximum_weight = min(
        MAX_QUOTIENT_WEIGHT,
        MAX_WEIGHTED_COST - prefix_cost,
        maximum_by_product,
    )
    if minimum_weight > maximum_weight:
        return None
    return minimum_weight, maximum_weight


def phase_features(weight: int, local_step: int) -> dict[str, int | str]:
    """Return exact logical endpoints used by one Algorithm-3 microstep."""
    if not 1 <= local_step <= 4 * weight:
        raise ValueError("local step outside current quotient")
    if local_step <= weight:
        # Phase A.  Pre-shift has already incremented l_s.
        return {"phase": "A", "l_q": 0, "l_s_r": local_step}
    if local_step <= 2 * weight:
        # Phase B.  R arithmetic and quotient swap precede l_q += 1.
        j = local_step - weight
        return {"phase": "B", "l_q": j - 1, "l_s_r": weight - j, "swap_q": j - 1}
    if local_step <= 3 * weight:
        # Phase C.  Quotient swap precedes l_q -= 1; T arithmetic follows it.
        j = local_step - 2 * weight
        return {"phase": "C", "swap_q": weight - j + 1, "l_s_t": j - 1}
    # Phase D.  T arithmetic precedes the post-shift decrement.
    j = local_step - 3 * weight
    return {"phase": "D", "l_s_t": weight - j + 1}


def envelope(values: list[tuple[int, int]]) -> list[int] | None:
    if not values:
        return None
    return [min(lo for lo, _ in values), max(hi for _, hi in values)]


def arithmetic_windows(step: int) -> tuple[dict[str, list[int] | None], dict[str, int]]:
    r_ranges: list[tuple[int, int]] = []
    swap_ranges: list[tuple[int, int]] = []
    t_ranges: list[tuple[int, int]] = []
    candidates = 0
    phase_counts = {phase: 0 for phase in "ABCD"}

    def intersect(lo: int, hi: int, phase_lo: int, phase_hi: int) -> tuple[int, int] | None:
        lo = max(lo, phase_lo)
        hi = min(hi, phase_hi)
        return (lo, hi) if lo <= hi else None

    for prefix_cost in range(0, min(MAX_WEIGHTED_COST - 1, (step - 1) // 4) + 1):
        if prefix_cost == 1:
            continue
        local_step = step - 4 * prefix_cost
        interval = feasible_weight_interval(prefix_cost, local_step)
        if interval is None:
            continue
        weight_lo, weight_hi = interval
        ell_t_min = PREFIX_LOWER[prefix_cost].bit_length()

        # A: 1 <= o <= w.
        phase_range = intersect(weight_lo, weight_hi, local_step, MAX_QUOTIENT_WEIGHT)
        if phase_range is not None:
            lo, hi = phase_range
            count = hi - lo + 1
            candidates += count
            phase_counts["A"] += count
            r_ranges.append((ell_t_min + 1, WORK_SIZE - local_step))

        # B: w < o <= 2w.
        phase_range = intersect(
            weight_lo, weight_hi, (local_step + 1) // 2, local_step - 1
        )
        if phase_range is not None:
            lo, hi = phase_range
            count = hi - lo + 1
            candidates += count
            phase_counts["B"] += count
            # j=o-w.  Both required edges decrease as w increases.
            r_ranges.append((ell_t_min + local_step - hi, WORK_SIZE + local_step - 2 * lo))
            ell_t_max = prefix_length_bounds(prefix_cost, lo)[1]
            swap_ranges.append(
                (
                    ell_t_min + local_step - hi,
                    min(N + 2, ell_t_max + local_step - lo),
                )
            )

        # C: 2w < o <= 3w.
        phase_range = intersect(
            weight_lo,
            weight_hi,
            (local_step + 2) // 3,
            (local_step - 1) // 2,
        )
        if phase_range is not None:
            lo, hi = phase_range
            count = hi - lo + 1
            candidates += count
            phase_counts["C"] += count
            ell_t_max_hi = prefix_length_bounds(prefix_cost, hi)[1]
            ell_t_max_lo = prefix_length_bounds(prefix_cost, lo)[1]
            # J=ell_t+3w-o+2.  The lower edge grows with w.  The upper
            # edge also grows: increasing w adds three while ell_t_max can
            # fall by at most one.
            swap_ranges.append(
                (
                    ell_t_min + 3 * lo - local_step + 2,
                    min(N + 2, ell_t_max_hi + 3 * hi - local_step + 2),
                )
            )
            t_ranges.append((1, min(N + 1, ell_t_max_lo + 1)))

        # D: 3w < o <= 4w.
        phase_range = intersect(
            weight_lo,
            weight_hi,
            (local_step + 3) // 4,
            (local_step - 1) // 3,
        )
        if phase_range is not None:
            lo, hi = phase_range
            count = hi - lo + 1
            candidates += count
            phase_counts["D"] += count
            ell_t_max = prefix_length_bounds(prefix_cost, lo)[1]
            t_ranges.append((1, min(N + 1, ell_t_max + 1)))

    return {
        "r_addsub": envelope(r_ranges),
        "quotient_swap": envelope(swap_ranges),
        "t_addsub": envelope(t_ranges),
    }, {"relaxed_candidates": candidates, **{f"phase_{k}": v for k, v in phase_counts.items()}}


def relaxed_prefix_states_at_cost(cost: int) -> bool:
    return 2 <= cost <= MAX_WEIGHTED_COST and PREFIX_LOWER[cost] is not None and PREFIX_LOWER[cost] <= P


def length_windows(step: int) -> dict[str, list[int] | None]:
    if step % 4:
        return {"len_update_lt": None, "len_update_lrp": None}
    cost = step // 4
    if not relaxed_prefix_states_at_cost(cost):
        return {"len_update_lt": None, "len_update_lrp": None}

    # Before the last quotient, at least max(0,cost-256) weighted cost has
    # already been consumed.  Scan down to the smallest possible prior
    # coefficient length and up through the largest possible new coefficient.
    prior_cost_floor = max(0, cost - MAX_QUOTIENT_WEIGHT)
    possible_prior_lengths = [
        PREFIX_LOWER[c].bit_length()
        for c in range(prior_cost_floor, cost)
        if PREFIX_LOWER[c] is not None
    ]
    k_lt = min(possible_prior_lengths)
    K_lt = min(N, cost)

    # highest_position_xor_write scans a dynamic boundary label as well as
    # the nonzero coefficient lanes.  If (t,t_next) is the coefficient pair
    # after weighted prefix cost c, induction on the quotient bit lengths
    # gives t+t_next <= 2**c.  The Euclidean invariant
    #
    #     p = r*t_next + r_next*t < r*(t_next+t)
    #
    # therefore gives r > p/2**c.  The prepared decoder label is
    # B = n+3-bit_length(r), so this exact integer lower bound on r supplies
    # a universal upper bound on B.  Omitting B is unsound even when every
    # nonzero coefficient lane itself lies inside the scan window.
    minimum_current_remainder = (P >> cost) + 1
    maximum_boundary_b = N + 3 - minimum_current_remainder.bit_length()
    K_lt = max(K_lt, maximum_boundary_b)

    # After this boundary at most 404-cost quotient-weight units remain.  A
    # suffix continuant of cost d is <2**d; the terminal old remainder is 1.
    remaining_cost = MAX_WEIGHTED_COST - cost
    maximum_remainder_length = max(1, min(N, remaining_cost))
    data_k_lrp = N + 4 - maximum_remainder_length

    # right_length_xor_write also requires the prepared decoder endpoint
    # A=bit_length(t_next)+2 to occur in the scan.  A suffix of weighted cost
    # d has continuant r < 2**d (and r=1 for an empty suffix).  Combining
    # p < 2*r*t_next with x<=p/2 gives the conservative lower bound below.
    # Taking the minimum with the data-lane bound covers both obligations.
    if remaining_cost == 0:
        maximum_current_remainder = 1
    else:
        maximum_current_remainder = min(P // 2, (1 << min(N, remaining_cost)) - 1)
    minimum_next_coefficient = P // (2 * maximum_current_remainder) + 1
    minimum_boundary_a = max(4, minimum_next_coefficient.bit_length() + 2)
    k_lrp = min(data_k_lrp, minimum_boundary_a)
    return {
        "len_update_lt": [k_lt, K_lt],
        "len_update_lrp": [k_lrp, WORK_SIZE],
    }


def ceil_safe(value: float, eps: float = 1e-12) -> int:
    return math.ceil(value - eps)


def floor_safe(value: float, eps: float = 1e-12) -> int:
    return math.floor(value + eps)


def paper_windows(step: int) -> dict[str, list[int]]:
    k1 = max(ceil_safe((step - (N + 2)) / (4.0 * C_EEA - 1.0)), 1) + 2
    k2 = max(ceil_safe((step - 3.0 * (N + 2)) / (4.0 * C_EEA - 3.0)), 1) + 1
    K2 = min(floor_safe(step / 2.0) + 2, N + 2)
    K3 = min(ceil_safe(step / 4.0) + 1, N + 1)
    k4 = max(ceil_safe((step - 4.0 * (N + 2)) / (4.0 * C_EEA - 4.0)), 1)
    K4 = min(floor_safe(step / 4.0 + 3.0), N + 3)
    k5 = ceil_safe(step / (4.0 * C_EEA))
    K5 = min(floor_safe(step / 4.0 + 4.0), N + 3)
    return {
        "r_addsub": [k1, N + 3],
        "quotient_swap": [k2, K2],
        "t_addsub": [1, K3],
        "len_update_lt": [k4, K4],
        "len_update_lrp": [k5, K5],
    }


def contains(outer: list[int], inner: list[int] | None) -> bool:
    return inner is None or (outer[0] <= inner[0] and outer[1] >= inner[1])


def exact_trace_requirements(x: int):
    """Yield exact required ranges for one secp input without using sampling."""
    r_previous, r = P, x
    t_previous, t = 0, 1
    prefix_cost = 0
    while r:
        quotient, r_next = divmod(r_previous, r)
        weight = quotient.bit_length()
        ell_t = t.bit_length()
        t_next = t_previous + quotient * t
        for local_step in range(1, 4 * weight + 1):
            step = 4 * prefix_cost + local_step
            features = phase_features(weight, local_step)
            phase = str(features["phase"])
            required: dict[str, list[int]] = {}
            if phase in "AB":
                ell_q = int(features["l_q"])
                ell_s = int(features["l_s_r"])
                required["r_addsub"] = [ell_t + ell_q + 1, WORK_SIZE - ell_s]
            if phase in "BC":
                selector = ell_t + int(features["swap_q"]) + 1
                required["quotient_swap"] = [selector, selector]
            if phase in "CD":
                required["t_addsub"] = [1, ell_t + 1]
            yield step, required

        boundary = 4 * (prefix_cost + weight)
        # The length decoders need their dynamic boundary labels to occur in
        # the scanned interval; covering only nonzero Work positions is not
        # sufficient for range_scan_leq/range_scan_geq to toggle correctly.
        boundary_b = N + 3 - r.bit_length()
        coefficient_positions = [t.bit_length(), t_next.bit_length(), boundary_b]
        remainder_positions = [N + 4 - r.bit_length()]
        if r_next:
            remainder_positions.append(N + 4 - r_next.bit_length())
        boundary_a = t_next.bit_length() + 2
        remainder_positions.append(boundary_a)
        yield boundary, {
            "len_update_lt": [min(coefficient_positions), max(coefficient_positions)],
            "len_update_lrp": [min(remainder_positions), max(remainder_positions)],
        }
        r_previous, r = r, r_next
        t_previous, t = t, t_next
        prefix_cost += weight


def concrete_paper_counterexamples() -> list[dict[str, object]]:
    witnesses = {
        "x_one": 1,
        "half_prime": P // 2,
        "schedule_1500": int(
            "5DB3D742C265539D92BA16B83C5C1DC492EC1A6629ED23CC63905323D8E62784", 16
        ),
        "schedule_1524": int(
            "5DB3D742C265539D92BA16B83C5C1DC492EC1A6629ED23CC63905323D96EFAEF", 16
        ),
    }
    wanted = {
        ("x_one", 1, "r_addsub"),
        ("half_prime", 8, "len_update_lrp"),
        ("schedule_1500", 240, "len_update_lrp"),
        ("schedule_1500", 1389, "r_addsub"),
        ("schedule_1500", 1470, "quotient_swap"),
        ("schedule_1524", 1472, "len_update_lt"),
    }
    found: list[dict[str, object]] = []
    for name, x in witnesses.items():
        for step, required_by_block in exact_trace_requirements(x):
            if step > 1476:
                continue
            paper = paper_windows(step)
            for block, required in required_by_block.items():
                key = (name, step, block)
                if key not in wanted:
                    continue
                raw = paper[block]
                repaired = [max(1, raw[0] - 1), raw[1]] if block == "r_addsub" else raw
                found.append({
                    "witness": name,
                    "x_hex": hex(x),
                    "step": step,
                    "block": block,
                    "required": required,
                    "paper": raw,
                    "paper_contains_required": contains(raw, required),
                    "one_lane_r_repair": repaired if block == "r_addsub" else None,
                    "one_lane_r_repair_contains_required": contains(repaired, required) if block == "r_addsub" else None,
                })
    if len(found) != len(wanted):
        missing = sorted(wanted - {(r["witness"], r["step"], r["block"]) for r in found})
        raise AssertionError(f"missing concrete paper counterexamples: {missing}")
    return sorted(found, key=lambda row: (int(row["step"]), str(row["block"])))


def build_table(certificate_path: Path) -> dict[str, object]:
    certificate_bytes = certificate_path.read_bytes()
    certificate = json.loads(certificate_bytes)
    result = certificate["result"]
    if int(result["weighted_cost_upper_bound"]) != MAX_WEIGHTED_COST:
        raise AssertionError("certificate weighted-cost bound changed")
    if int(result["safe_fixed_schedule_steps"]) != SAFE_STEPS:
        raise AssertionError("certificate fixed schedule changed")
    if int(certificate["p_decimal"]) != P:
        raise AssertionError("certificate field prime changed")

    rows = []
    paper_first_empty: dict[str, int] = {}
    for step in range(1, SAFE_STEPS + 1):
        arithmetic, counts = arithmetic_windows(step)
        safe = {**arithmetic, **length_windows(step)}
        paper = paper_windows(step)
        for block, window in paper.items():
            if window[0] > window[1] and block not in paper_first_empty:
                paper_first_empty[block] = step
        rows.append({
            "step": step,
            "safe": safe,
            "paper": paper,
            "paper_contains_proved_envelope": {
                block: contains(paper[block], safe[block]) for block in safe
            },
            "proof_state_counts": counts,
        })

    return {
        "schema": "luo-secp256k1-active-windows-v2",
        "field": "secp256k1",
        "p_hex": hex(P),
        "n": N,
        "work_size": WORK_SIZE,
        "weighted_cost_bound": MAX_WEIGHTED_COST,
        "fixed_schedule_steps": SAFE_STEPS,
        "shift_register_requirement": {
            "minimum_exact_steps": 1024,
            "minimum_exact_steps_reason": "continuant p is below 2^sum_weights, so sum_weights>=256; x=1 attains 256",
            "maximum_terminal_padding_steps": SAFE_STEPS - 1024,
            "maximum_terminal_padding_witness_x": "0x1",
            "required_counter_bits": (SAFE_STEPS - 1024).bit_length(),
            "warning": "the existing 9-bit l_s wraps after 511 and cannot canonicalize 592 physical rotations modulo 259",
        },
        "certificate": {
            "path": certificate_path.name,
            "sha256": hashlib.sha256(certificate_bytes).hexdigest(),
        },
        "semantics": {
            "ranges": "inclusive 1-based physical Work labels; null means block is unreachable at that step",
            "r_addsub": "[min(L-1),max(R)], L=ell_t+ell_q+2, R=n+3-ell_s after pre-shift",
            "quotient_swap": "selector J=ell_t+ell_q+1; gate additionally exposes Work[J+1]",
            "t_addsub": "[1,max(ell_t+1)]",
            "len_update_lt": "covers both coefficient fields and decoder label B=n+3-bit_length(r) at an iteration boundary",
            "len_update_lrp": "covers both remainder fields and decoder label A=bit_length(t_next)+2 at an iteration boundary",
        },
        "proof_relaxations": [
            "every canonical quotient word with continuant p has weighted cost at most 404",
            "prefix K >= max(2^(cost-count), Fibonacci[count+2])",
            "prefix K < 2^cost and 2^(current_weight-1)*K < p",
            "suffix continuant of remaining weighted cost d is below 2^d",
            "prefix coefficient sum t+t_next is at most 2^cost, hence B is explicitly bounded",
            "p < 2*r*t_next and the suffix bound give an explicit lower bound on decoder label A",
            "all relaxed prefix lengths, current quotient weights, and four phase positions are enumerated",
        ],
        "paper_first_empty_step": paper_first_empty,
        "concrete_paper_counterexamples": concrete_paper_counterexamples(),
        "rows": rows,
    }


def main() -> None:
    here = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser()
    parser.add_argument("--certificate", type=Path, default=here / "certificate.json")
    parser.add_argument("--out", type=Path, default=here / "active_windows_1616.json")
    parser.add_argument("--tail-csv", type=Path, default=here / "active_windows_1477_1616.csv")
    args = parser.parse_args()
    table = build_table(args.certificate)
    encoded = json.dumps(table, indent=2, sort_keys=True) + "\n"
    args.out.write_text(encoded, encoding="utf-8")
    with args.tail_csv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        blocks = ["r_addsub", "quotient_swap", "t_addsub", "len_update_lt", "len_update_lrp"]
        writer.writerow(["step", *[f"{block}_lo" for block in blocks], *[f"{block}_hi" for block in blocks]])
        for row in table["rows"][1476:]:
            safe = row["safe"]
            writer.writerow(
                [row["step"]]
                + [(safe[block][0] if safe[block] is not None else "") for block in blocks]
                + [(safe[block][1] if safe[block] is not None else "") for block in blocks]
            )
    print(f"wrote={args.out}")
    print(f"sha256={hashlib.sha256(encoded.encode()).hexdigest()}")
    print(f"tail_csv={args.tail_csv}")
    print(f"rows={len(table['rows'])}")
    print(f"paper_first_empty_step={table['paper_first_empty_step']}")
    for row in table["concrete_paper_counterexamples"]:
        print(
            f"counterexample step={row['step']} block={row['block']} "
            f"required={row['required']} paper={row['paper']} x={row['x_hex']}"
        )


if __name__ == "__main__":
    main()
