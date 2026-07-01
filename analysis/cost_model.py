#!/usr/bin/env python3
"""Map the abstract challenge score to a physical fault-tolerant cost estimate.

The challenge scores  score = round(avg_toffoli_per_shot) * qubits
(src/bin/eval_circuit.rs:434-435), where `toffoli` counts CCX+CCZ executions
(src/sim.rs:86) and `qubits` = max allocated qubit id + 1 (src/circuit.rs:356).
That product is a proxy; it says nothing physical on its own. This script turns
the two real metrics into surface-code resource estimates under explicitly
stated assumptions, so the number means something in qubit-seconds.

EVERYTHING numeric is either (a) read from score.json, or (b) an assumption
printed in the ASSUMPTIONS block below and applied deterministically. No number
is invented; change an assumption and re-run.

Physical model references (assumptions, not repo facts):
  - Fowler, Mariantoni, Martinis, Cleland 2012 (surface codes) -- patch = 2 d^2.
  - Gidney & Ekera 2021, "How to factor 2048-bit RSA in 8 hours" (arXiv:1905.09749)
    -- reaction-limited runtime, t_react ~ 10 us, d ~ 27 at p=1e-3.
  - Gidney 2018 (arXiv:1805.03662) -- Toffoli via measurement = 4 T (repo uses
    measurement-based uncompute, so 4 T/Toffoli is the apt convention; 7 T is the
    Clifford+T textbook upper bound).
  - Roetteler, Naehrig, Svore, Lauter 2017 (arXiv:1706.06752) -- full n-bit ECDLP
    is O(n) point additions; the multiplier below is the (assumed) ladder factor.
"""
import json
import math
import os

HERE = os.path.dirname(os.path.abspath(__file__))
SCORE = os.path.normpath(os.path.join(HERE, "..", "score.json"))

# ----------------------------- REAL INPUTS ---------------------------------
with open(SCORE) as f:
    sj = json.load(f)
TOFFOLI = sj["metrics"]["toffoli"]     # avg CCX+CCZ per shot, rounded
QUBITS = sj["metrics"]["qubits"]       # logical width = max qubit id + 1
SCORE_VAL = sj["score"]

# Optional measured depth from `depth_report` (src/bin/depth_report.rs -> depth.json).
DEPTH = os.path.normpath(os.path.join(HERE, "..", "depth.json"))
TOFFOLI_DEPTH = None
GATE_DEPTH = None
if os.path.exists(DEPTH):
    with open(DEPTH) as f:
        dj = json.load(f)
    TOFFOLI_DEPTH = dj.get("toffoli_depth")
    GATE_DEPTH = dj.get("gate_depth")

# ----------------------------- ASSUMPTIONS ---------------------------------
A = {
    "p_phys": 1e-3,          # physical gate/measurement error rate
    "p_th": 1e-2,            # surface-code threshold (~1%)
    "t_cycle_us": 1.0,       # surface-code cycle time (superconducting)
    "t_react_us": 10.0,      # feed-forward reaction time (Gidney-Ekera)
    "T_per_toffoli": 4,      # measurement-based Toffoli (repo technique); 7 = textbook
    "phys_per_patch": lambda d: 2 * d * d,   # physical qubits per logical patch
    "factory_routing_overhead": 2.0,         # x logical patches for factories+routing
    "distances": [21, 25, 27],
    "ecdlp_point_additions": 1600,  # ASSUMED full-attack ladder factor ~ O(n), n=256
    "target_fail_prob": 0.01,
}


def logical_err_per_cycle(d, p, p_th):
    # standard phenomenological fit: p_L ~ 0.1 (p/p_th)^((d+1)/2)
    return 0.1 * (p / p_th) ** ((d + 1) / 2)


def section(t):
    print("\n" + t + "\n" + "-" * len(t))


print("=" * 68)
print(" ecdsafail-challenge  ->  physical fault-tolerant cost estimate")
print("=" * 68)

section("REAL INPUTS (score.json)")
print(f"  Toffoli (CCX+CCZ, avg/shot, rounded) : {TOFFOLI:,}")
print(f"  Logical qubits (max id + 1)          : {QUBITS:,}")
print(f"  Challenge score (Toffoli x qubits)   : {SCORE_VAL:,}")
print("  NOTE: this circuit = ONE elliptic-curve point addition.")
if TOFFOLI_DEPTH is not None:
    print(f"  Toffoli-depth (measured, depth.json) : {TOFFOLI_DEPTH:,}")
    print(f"  Gate depth   (measured, depth.json)  : {GATE_DEPTH:,}")
    print(f"  Toffoli parallelism (gates/depth)    : {TOFFOLI/TOFFOLI_DEPTH:.2f}x")
else:
    print("  (no depth.json found; run `cargo run --release --bin depth_report`")
    print("   after build_circuit to get measured depth. Runtimes are sequential UBs.)")

section("ASSUMPTIONS (edit + re-run)")
for k, v in A.items():
    if callable(v):
        v = "2*d^2"
    print(f"  {k:26s} = {v}")

section("LOGICAL NON-CLIFFORD VOLUME")
for tpt in (4, 7):
    print(f"  T-count @ {tpt} T/Toffoli : {TOFFOLI * tpt:,}"
          + ("   <- measurement-based (repo)" if tpt == 4 else "   <- Clifford+T textbook"))

section("PER-POINT-ADDITION SURFACE-CODE RESOURCES")
print(f"  {'d':>3} | {'p_L/cycle':>10} | {'phys/patch':>10} | {'phys qubits (incl. factories+routing)':>38}")
for d in A["distances"]:
    pl = logical_err_per_cycle(d, A["p_phys"], A["p_th"])
    per_patch = A["phys_per_patch"](d)
    phys = int(QUBITS * per_patch * A["factory_routing_overhead"])
    print(f"  {d:>3} | {pl:>10.2e} | {per_patch:>10,} | {phys:>38,}")

section("RUNTIME (reaction-limited: one t_react per non-Clifford LAYER)")
# One Toffoli/CCZ consumes one magic state; with enough factories the wall-clock
# is set by the non-Clifford critical path (toffoli-depth), not the total count.
t_seq_s = TOFFOLI * A["t_react_us"] * 1e-6
print(f"  sequential upper bound (all Toffolis) : {t_seq_s:,.1f} s   [count x t_react]")
d_hi = A["distances"][-1]
phys_hi = int(QUBITS * A["phys_per_patch"](d_hi) * A["factory_routing_overhead"])
if TOFFOLI_DEPTH is not None:
    t_meas_s = TOFFOLI_DEPTH * A["t_react_us"] * 1e-6
    print(f"  MEASURED (toffoli-depth x t_react)    : {t_meas_s:,.2f} s   [depth.json]")
    # spacetime volume = physical qubits held x wall-clock, in qubit-seconds.
    vol = phys_hi * t_meas_s
    print(f"  spacetime volume @ d={d_hi}             : {vol:,.3e} physical-qubit-seconds")
    print(f"  (= {phys_hi:,} physical qubits x {t_meas_s:.2f} s)")
else:
    print("  MEASURED runtime unavailable (no depth.json).")
print("  NOTE: depth is the reaction-limited critical path; true wall-clock also")
print("        depends on having enough magic-state factories to feed each layer.")

section("EXTRAPOLATION TO A FULL secp256k1 ECDLP BREAK (order-of-magnitude)")
mult = A["ecdlp_point_additions"]
tof_full = TOFFOLI * mult
print(f"  ASSUMED point additions in the full Shor-ECDLP ladder : ~{mult:,}  (O(n), n=256)")
print(f"  => full-attack Toffoli count : ~{tof_full:,.0f}  (~{tof_full:.1e})")
print(f"  => full-attack T-count @4    : ~{tof_full*4:.1e}")
print(f"  physical qubits @ d={d_hi}        : ~{phys_hi:,}  (single-addition width; the full")
print(f"     algorithm needs a wider register file -> not derivable from this repo)")
if TOFFOLI_DEPTH is not None:
    full_rt_h = TOFFOLI_DEPTH * mult * A["t_react_us"] * 1e-6 / 3600
    print(f"  reaction-limited runtime     : ~{full_rt_h:,.1f} h  (depth-based, if additions are serial)")
else:
    print(f"  reaction-limited runtime UB  : ~{tof_full*A['t_react_us']*1e-6/3600:,.1f} h  (count-based)")
print("\n  Caveat: the multiplier and full-algorithm width are ASSUMPTIONS pending a")
print("  full-circuit build; only the per-point-addition figures are measured.")
print("  See Roetteler et al. 2017 for the exact ladder structure.")
print("=" * 68)
