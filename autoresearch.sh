#!/bin/bash
set -euo pipefail

note=${AUTORESEARCH_NOTE:-}
if [[ -z "$note" && -f autoresearch.note ]]; then
  note=$(tr '\n' ' ' < autoresearch.note)
fi
if [[ -z "$note" ]]; then
  note=autoresearch
fi

cargo build --release >/dev/null 2>&1

if [[ -f autoresearch.extra_test ]]; then
  extra_test=$(tr '\n' ' ' < autoresearch.extra_test | xargs)
  if [[ -n "$extra_test" ]]; then
    cargo test --release "$extra_test" -- --nocapture
  fi
fi

tmp=$(mktemp)
cleanup() {
  rm -f "$tmp"
}
trap cleanup EXIT

cargo run --release -- --note "$note" 2>&1 | tee "$tmp"

toffoli=$(awk -F: '/avg executed Toffoli/{gsub(/ /, "", $2); print $2; exit}' "$tmp")
clifford=$(awk -F: '/avg executed Clifford/{gsub(/ /, "", $2); print $2; exit}' "$tmp")
qubits=$(awk -F: '/qubits/{gsub(/ /, "", $2); print $2; exit}' "$tmp")
ops=$(awk -F: '/emitted ops/{gsub(/ /, "", $2); print $2; exit}' "$tmp")

if [[ -n ${toffoli:-} ]]; then
  echo "METRIC avg_toffoli=$toffoli"
fi
if [[ -n ${clifford:-} ]]; then
  echo "METRIC avg_clifford=$clifford"
fi
if [[ -n ${qubits:-} ]]; then
  echo "METRIC qubits=$qubits"
fi
if [[ -n ${ops:-} ]]; then
  echo "METRIC emitted_ops=$ops"
fi
if grep -q '=== experiment OK ===' "$tmp"; then
  echo "METRIC correctness_ok=1"
fi
