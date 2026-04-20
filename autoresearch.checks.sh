#!/bin/bash
set -euo pipefail

tmp=$(mktemp)
cleanup() {
  rm -f "$tmp"
}
trap cleanup EXIT

cargo run --release -- --note "autoresearch-checks" >"$tmp" 2>&1 || {
  tail -80 "$tmp"
  exit 1
}

# Qubit budget guard: program.md hard cap is 3700.
qubits=$(awk -F: '/qubits/{gsub(/ /, "", $2); print $2; exit}' "$tmp")
if [[ -n "${qubits:-}" ]] && (( qubits > 3700 )); then
  echo "CHECKS FAIL: peak qubits $qubits exceeds hard cap 3700"
  exit 1
fi
