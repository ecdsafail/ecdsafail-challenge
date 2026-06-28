#!/usr/bin/env python3
"""
z_calibrate.py - Calibrate Muqatta-Z defense parameters from ecdsa.fail scores.

Usage:
    python z_calibrate.py                     # Uses ../../../score.json
    python z_calibrate.py --score 1571592960   # Manual score input
    python z_calibrate.py --track              # Append score to history

Author: Peter Anari Otuke, University of Nairobi
"""
import argparse
import json
import os
from datetime import datetime, timezone

import numpy as np
from scipy import linalg

# Canonical Z-matrix (seed=42). Symmetric 5x5 coupling matrix over the
# zones: ECDLP, Nonce, Verify, SideCh, PQC.
np.random.seed(42)
Z_BASE = np.array([
    [8.50, 1.20, 0.80, 1.50, 2.30],
    [1.20, 7.80, 0.60, 2.10, 1.80],
    [0.80, 0.60, 6.20, 0.40, 0.90],
    [1.50, 2.10, 0.40, 9.10, 1.70],
    [2.30, 1.80, 0.90, 1.70, 10.40],
], dtype=np.float64)
ZONES = ["ECDLP", "Nonce", "Verify", "SideCh", "PQC"]

# Score recorded at the time the framework was calibrated.
C_BASELINE = 10_758_874_395
GAMMA = 0.5


def eigendecompose(z):
    vals, vecs = linalg.eigh(z)
    idx = np.argsort(vals)[::-1]
    return vals[idx], vecs[:, idx]


def calibrate_z_matrix(current_score):
    """Scale the ECDLP<->PQC coupling by the leaderboard improvement ratio."""
    improvement_ratio = C_BASELINE / current_score
    z15_multiplier = improvement_ratio ** GAMMA

    z_cal = Z_BASE.copy()
    z_cal[0, 4] *= z15_multiplier
    z_cal[4, 0] *= z15_multiplier
    z_cal[1, 4] *= z15_multiplier ** 0.5
    z_cal[4, 1] *= z15_multiplier ** 0.5
    return z_cal, z15_multiplier, improvement_ratio


def full_analysis(current_score):
    z_cal, z15_mult, imp_ratio = calibrate_z_matrix(current_score)
    vals, vecs = eigendecompose(z_cal)
    l1, l2, l5 = vals[0], vals[1], vals[-1]

    kappa = l1 / l5
    alpha = l2 / l1
    beta = 0.5
    cr = l2 / (beta * np.sqrt(l1 * l5))

    v1 = vecs[:, 0]
    wstar = np.abs(v1) / np.sum(np.abs(v1))

    eps = 1e-6
    k = int(np.ceil(np.log(eps) / np.log(alpha)))

    w_e, w_nw, w_pq = 71, 43, 3293  # witness-size inputs (bytes)
    t_target = 0.80
    d_cont = t_target * w_pq / (w_e + (1 - t_target) * w_nw)
    d_star = int(np.ceil(d_cont))
    dd_dw = t_target / (w_e + (1 - t_target) * w_nw)

    gershgorin_margins = [
        z_cal[i, i] - (np.sum(np.abs(z_cal[i, :])) - z_cal[i, i])
        for i in range(z_cal.shape[0])
    ]

    return {
        "score": current_score,
        "improvement_ratio": imp_ratio,
        "z15_multiplier": z15_mult,
        "eigenvalues": vals,
        "kappa": kappa,
        "alpha": alpha,
        "cr": cr,
        "K": k,
        "wstar": dict(zip(ZONES, wstar)),
        "D_star": d_star,
        "D_cont": d_cont,
        "dD_dW": dd_dw,
        "gershgorin_min": min(gershgorin_margins),
    }


def load_score_json(path):
    try:
        with open(path) as f:
            return json.load(f)["score"]
    except (FileNotFoundError, KeyError, json.JSONDecodeError):
        return None


def print_report(result):
    print("=" * 72)
    print("MUQATTA-Z DEFENSE CALIBRATION - ecdsa.fail Integration")
    print("=" * 72)
    print(f"\n  INPUT (ecdsa.fail leaderboard):")
    print(f"    Current score:      {result['score']:,.0f}")
    print(f"    Baseline score:     {C_BASELINE:,.0f}")
    print(f"    Improvement ratio:  {result['improvement_ratio']:.2f}x")
    print(f"    Z15 multiplier:     {result['z15_multiplier']:.4f}x (gamma={GAMMA})")

    vs = result["eigenvalues"]
    print(f"\n  SPECTRAL (calibrated Z-matrix):")
    print(f"    lambda = [{', '.join(f'{v:.6f}' for v in vs)}]")
    print(f"    kappa = {result['kappa']:.6f}")
    print(f"    alpha = {result['alpha']:.6f}")
    print(f"    CR = {result['cr']:.6f} ({'PASS' if result['cr'] >= 1.0 else 'FAIL'})")
    print(f"    K = {result['K']} iterations")
    print(f"    Gershgorin min = {result['gershgorin_min']:.6f}")

    print(f"\n  OPTIMAL ALLOCATION w*:")
    for zone, w in result["wstar"].items():
        bar = "#" * int(w * 40)
        print(f"    {zone:8s}: {w*100:6.2f}%  {bar}")

    print(f"\n  SEGWIT-PQ:")
    print(f"    D* = {result['D_star']} (continuous: {result['D_cont']:.4f})")
    print(f"    dD*/dW_pq = {result['dD_dW']:.6f} per byte")
    print("=" * 72)


def track_score(score, history_file="defense_history.jsonl"):
    entry = {
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "score": score,
        "improvement_ratio": C_BASELINE / score,
    }
    with open(history_file, "a") as f:
        f.write(json.dumps(entry) + "\n")
    print(f"  Tracked: score={score:,.0f} at {entry['timestamp']}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Calibrate Muqatta-Z defense from ecdsa.fail scores"
    )
    parser.add_argument("--score", type=int, default=None,
                         help="Circuit score (Toffoli x qubits)")
    parser.add_argument("--track", action="store_true",
                         help="Append score to tracking history")
    args = parser.parse_args()

    score = args.score
    if score is None:
        here = os.path.dirname(os.path.abspath(__file__))
        score = load_score_json(os.path.join(here, "..", "..", "..", "score.json"))
    if score is None:
        score = C_BASELINE
        print(f"No score.json found; using baseline score {score:,.0f}")

    result = full_analysis(score)
    print_report(result)

    if args.track:
        track_score(score)
