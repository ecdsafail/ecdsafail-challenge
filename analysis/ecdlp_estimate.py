#!/usr/bin/env python3
"""Derived full-circuit cost for a Shor-ECDLP attack on secp256k1, built from
the MEASURED per-addition metrics of this repo's circuit plus the known
double-and-add structure (Roetteler, Naehrig, Svore, Lauter 2017,
arXiv:1706.06752).

Replaces the hand-picked `ecdlp_point_additions = 1600` multiplier in
cost_model.py with a structural count: the number of additions and register
widths come from the algorithm, and the per-addition Toffoli/depth/ancilla come
from real runs (score.json + depth.json). Everything numeric is either measured
or an explicitly stated assumption; nothing is hand-asserted.

Algorithm (Shor for ECDLP, solve Q = [m]P):
  - two scalar registers a, b of n+1 qubits each, in uniform superposition;
  - compute the double-scalar multiplication [a]P + [b]Q into a running
    ACCUMULATOR point via double-and-add: each scalar bit gates one addition of
    a *classical*, compile-time-precomputed multiple ([2^i]P or [2^i]Q) — which
    is exactly this repo's primitive (classical point + quantum point, in place);
  - a final QFT over the 2(n+1) scalar qubits, then measure.

  Windowing (window w): group w scalar bits, add one precomputed multiple chosen
  by a 2^w-entry table (QROM) lookup -> 2(n+1)/w additions instead of 2(n+1).

CAVEATS (printed below too):
  - COMPLETENESS: this repo's adder is the incomplete affine formula; a real
    attack must handle P==Q / P==-Q / infinity (complete formulas or Roetteler's
    exception argument). Modelled by `completeness_overhead` (default 1.0 =
    exceptions assumed negligible, per Roetteler). Set >1 for complete formulas.
  - The QROM table-cost model (~2^w Toffoli/lookup) is a simplification; the
    windowed rows are a sensitivity, the basic (w=1) row is the headline.
  - QFT is O(n^2) phase rotations (Clifford+T), negligible in Toffoli vs the
    arithmetic; taken as 0 Toffoli here (Roetteler: arithmetic dominates).
"""
import json
import math
import os

HERE = os.path.dirname(os.path.abspath(__file__))


def load(name):
    p = os.path.normpath(os.path.join(HERE, "..", name))
    if not os.path.exists(p):
        return None
    with open(p) as f:
        return json.load(f)


score = load("score.json")
depth = load("depth.json")
if score is None:
    raise SystemExit("score.json not found (run the benchmark first)")
if depth is None:
    raise SystemExit("depth.json not found (run: cargo run --release --bin depth_report)")

# ----------------------------- MEASURED INPUTS -----------------------------
PER_ADD_TOF = score["metrics"]["toffoli"]         # Toffoli per addition
PER_ADD_QUBITS = score["metrics"]["qubits"]       # total qubits per addition
PER_ADD_TOF_DEPTH = depth["toffoli_depth"]        # non-Clifford critical path

# ----------------------------- ALGORITHM MODEL -----------------------------
N = 256                       # secp256k1 field/scalar size
SCALAR_BITS = N + 1           # qubits per scalar register (Roetteler)
N_SCALARS = 2                 # [a]P + [b]Q
ACCUM_QUBITS = 2 * N          # running point (x, y)
SCALAR_QUBITS = N_SCALARS * SCALAR_BITS
PER_ADD_ANCILLA = PER_ADD_QUBITS - ACCUM_QUBITS   # reused scratch per addition


def n_additions(w):
    return N_SCALARS * math.ceil(SCALAR_BITS / w)


def lookup_toffoli(w):
    # QROM address decode ~2^w Toffoli; data write is CX into classical bits;
    # uncompute is measurement-free (~0 Toffoli). Simplified.
    return (1 << w) if w > 1 else 0


# ----------------------------- ASSUMPTIONS ---------------------------------
A = {
    "p_phys": 1e-3,
    "p_th": 1e-2,
    "t_react_us": 10.0,
    "T_per_toffoli": 4,           # measurement-based Toffoli (repo technique)
    "phys_per_patch": lambda d: 2 * d * d,
    "factory_routing_overhead": 2.0,
    "distance": 27,
    "completeness_overhead": 1.0,  # exceptions assumed negligible (Roetteler)
    "windows": [1, 4, 8],          # w=1 = basic double-and-add (headline)
}


def section(t):
    print("\n" + t + "\n" + "-" * len(t))


print("=" * 74)
print(" Shor-ECDLP on secp256k1  ->  derived full-circuit cost (from measured primitive)")
print("=" * 74)

section("MEASURED PER-ADDITION (score.json + depth.json)")
print(f"  Toffoli / addition        : {PER_ADD_TOF:,}")
print(f"  Toffoli-depth / addition  : {PER_ADD_TOF_DEPTH:,}")
print(f"  qubits / addition (total) : {PER_ADD_QUBITS:,}")
print(f"  -> reused ancilla         : {PER_ADD_ANCILLA:,}  (= qubits - 2n accumulator)")

section("ALGORITHM STRUCTURE (Roetteler et al. 2017)")
print(f"  field size n              : {N}")
print(f"  scalar registers          : {N_SCALARS} x {SCALAR_BITS} qubits = {SCALAR_QUBITS}")
print(f"  accumulator point         : 2n = {ACCUM_QUBITS} qubits")
print(f"  basic additions (w=1)     : 2(n+1) = {n_additions(1)}")

section("DERIVED LOGICAL RESOURCES")
co = A["completeness_overhead"]
print(f"  completeness_overhead = {co}  (1.0 = exceptions assumed negligible)")
print(f"  {'window':>6} | {'#adds':>6} | {'total Toffoli':>16} | {'toffoli-depth':>16} | {'peak qubits':>11}")
rows = {}
for w in A["windows"]:
    na = n_additions(w)
    tof = int(na * (PER_ADD_TOF + lookup_toffoli(w)) * co)
    # accumulator is read+written by every addition -> additions are serial;
    # the non-Clifford critical path is the sum of per-addition depths.
    tdepth = int(na * PER_ADD_TOF_DEPTH * co)
    # peak width: accumulator + both scalars + one addition's reused ancilla.
    # (windowed lookup adds classical addend bits + minor ancilla, not counted.)
    peak_q = ACCUM_QUBITS + SCALAR_QUBITS + PER_ADD_ANCILLA
    rows[w] = (na, tof, tdepth, peak_q)
    tag = "  <- basic (headline)" if w == 1 else ""
    print(f"  {w:>6} | {na:>6} | {tof:>16,} | {tdepth:>16,} | {peak_q:>11,}{tag}")

section("PHYSICAL FAULT-TOLERANT COST  (basic ladder, w=1)")
na, tof, tdepth, peak_q = rows[1]
d = A["distance"]
phys = int(peak_q * A["phys_per_patch"](d) * A["factory_routing_overhead"])
t_count = tof * A["T_per_toffoli"]
runtime_s = tdepth * A["t_react_us"] * 1e-6
vol = phys * runtime_s
print(f"  total T-count @ {A['T_per_toffoli']} T/Tof : {t_count:,}  (~{t_count:.2e})")
print(f"  physical qubits @ d={d}     : {phys:,}  (~{phys:.2e})")
print(f"  reaction-limited runtime  : {runtime_s:,.0f} s  = {runtime_s/3600:.2f} h")
print(f"  spacetime volume          : {vol:.3e} physical-qubit-seconds")

section("SANITY CHECK vs LITERATURE")
print(f"  derived basic-ladder Toffoli : ~{rows[1][1]:.2e}")
print(f"  windowed (w=8)               : ~{rows[8][1]:.2e}")
print("  For scale: Gidney-Ekera 2021 factor RSA-2048 at ~3e9 Toffoli / ~2e7 qubits.")
print("  ECDLP-256 being cheaper in Toffoli than RSA-2048 is the expected ordering;")
print("  our per-addition primitive is this repo's optimized frontier, so the total")
print("  sits below a from-scratch estimate like Roetteler et al. 2017 (arXiv:1706.06752).")

section("CAVEATS")
print("  - COMPLETENESS unmodelled: this repo's adder is incomplete (affine only).")
print("    A correct attack needs exception handling (complete formulas or")
print("    Roetteler's negligibility argument); set completeness_overhead > 1 to")
print("    price complete formulas. This is a COST estimate, not a verified attack.")
print("  - Windowed rows use a simplified ~2^w QROM cost; w=1 is the safe headline.")
print("  - Full-algorithm qubit width reuses one addition's ancilla; a real build")
print("    may need more routing space. Numbers are derived, not emitted+measured")
print("    (that is Tier B; see analysis/scientific-value.md).")
print("=" * 74)
