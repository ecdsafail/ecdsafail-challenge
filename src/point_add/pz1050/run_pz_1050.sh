#!/usr/bin/env bash
# Decompress the bundled shrunken-PZ kmx and run the official benchmark.
# Result: qubits 1050, all 9024 shots OK (clean).
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$HERE/../../.." && pwd)"
KMX="${KMX:-/tmp/ec_shrunken_pz.kmx}"

if [ ! -f "$KMX" ]; then
  echo "decompressing bundled kmx -> $KMX ..."
  zstd -d --long=27 -f "$HERE/ec_shrunken_pz.kmx.zst" -o "$KMX"
fi

cd "$REPO_ROOT"
echo "running: POINT_ADD_FROM_KMX=$KMX ecdsafail run"
POINT_ADD_FROM_KMX="$KMX" ecdsafail run
