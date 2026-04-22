# Phase-fix strategy notes

## Confirmed facts
- The strict harness failures for many `k` are dominated by **phase-garbage batches**.
- A targeted probe on the first failing batch for `k = 4` showed the generic and
  experimental circuits already differ in phase **before `lam` is freed`.
- The same divergence persists after `lam` free, so `lam` free itself is not the cause.
- The local specialized step, local forward/backward composition, and isolated
  inverse identity tests do not show the same phase divergence on their own.

## One attempted fix that failed fast
I tried a surgical hypothesis test:
- replace the explicit `kaliski_backward` inside `with_kal_inv_raw` with the
  exact gate-inverse of `kaliski_forward` under an environment switch.

Result:
- this cannot be used directly because `kaliski_forward` contains HMR-based
  measurement uncomputation internally,
- and the local `emit_inverse` helper rejects non-unitary ops like `Hmr`.

So the simple “just use emit_inverse on the full forward block” fallback is not
available without a larger structural rewrite.

## Current likely direction
The phase bug probably lives in a mismatch between:
- the specialized prefix replacement,
- the explicit measurement-based backward strategy,
- and the way those interact with the later top-level scaffold.

The next realistic fix avenue is therefore not a one-line emit-inverse swap, but
carefully matching the specialized prefix's phase history to the generic
explicit backward / uncompute regime.
