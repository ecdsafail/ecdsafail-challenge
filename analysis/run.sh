#!/usr/bin/env bash
# Run the full scientific-rigor suite: formal proofs + physical cost model.
set -euo pipefail
cd "$(dirname "$0")"

echo "### 1/5  Solinas modular-reduction proof (z3) ###"
python3 verify/solinas_reduction.py
echo
echo "### 2/5  Peephole / adder / comparator proofs (z3) ###"
python3 verify/peephole_identities.py
echo
echo "### 3/5  Reference kickmix adder validation (source-paper artifacts) ###"
python3 verify/validate_reference_adders.py
echo
echo "### 4/5  Physical fault-tolerant cost model ###"
python3 cost_model.py
echo
echo "### 5/5  Derived full-ECDLP cost (measured primitive x paper's ladder) ###"
python3 ecdlp_estimate.py
