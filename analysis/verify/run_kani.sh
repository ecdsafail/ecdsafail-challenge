#!/usr/bin/env bash
# Kani (bit-precise BMC) proofs that bind to the REAL Rust types/functions.
# Harnesses live in src/kani_proofs.rs, gated behind #[cfg(kani)] so the normal
# build and benchmark.sh never compile them. Requires: cargo kani (v0.66+).
set -euo pipefail
cd "$(git rev-parse --show-toplevel 2>/dev/null || echo "$(dirname "$0")/../..")"

for h in solinas_add_u64 solinas_add_u256; do
  echo "### cargo kani --harness $h ###"
  cargo kani --harness "$h" 2>&1 | grep -E "VERIFICATION|failed|Verification Time" || true
  echo
done
