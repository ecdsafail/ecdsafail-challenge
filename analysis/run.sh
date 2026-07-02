#!/usr/bin/env bash
# Run the full scientific-rigor suite: formal proofs + physical cost model.
set -euo pipefail
cd "$(dirname "$0")"

echo "### 1/4  Solinas modular-reduction proof (z3) ###"
python3 verify/solinas_reduction.py
echo
echo "### 2/4  Peephole / adder / comparator proofs (z3) ###"
python3 verify/peephole_identities.py
echo
echo "### 3/4  Physical fault-tolerant cost model ###"
python3 cost_model.py
echo
echo "### 4/4  Derived full-ECDLP cost (measured primitive x ladder structure) ###"
python3 ecdlp_estimate.py
