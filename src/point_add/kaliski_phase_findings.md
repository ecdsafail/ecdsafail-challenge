# Kaliski bulk-prefix phase findings

## Strongest current localization result
A small targeted probe for the first strict phase-failing batch at `k = 4`
showed:

- the first failing batch under the full `main.rs` harness occurs at batch index **10**,
- and for that exact 64-shot batch, the generic and experimental circuits
  already differ in phase **before `lam` is freed**.

Measured phase masks at the cut immediately before `b.free_vec(&lam)`:

- generic: `0x0000040000000000`
- experimental: `0x0000000000000000`

The same masks persist after `lam` is freed, so the divergence is **not caused
by freeing `lam`** itself.

## Interpretation
This means the phase bug enters **earlier** than the final `lam` cleanup.
More specifically:
- it is not in the isolated specialized step,
- not in the isolated inverse identity,
- not in the pair1-only probe,
- not in the pair1+pair2 probe on ad hoc sampled points,
- but it **does** appear on the first actual circuit-seeded failing batch once
  the full top-level scaffold is used up to the late cut before `lam` free.

So the current best hypothesis is:

> the phase defect is introduced somewhere in the full point-add scaffold before
> `lam` is freed, and it is tied to the actual circuit-seeded batch family that
> triggers the strict harness failure.

## Immediate next target
The next high-value step is to bisect the top-level scaffold between:
- after pair1,
- after the `tx <- Rx - Qx` algebra,
- after `mul3_between_pair`,
- during pair2,

using the **same failing batch** rather than fresh ad hoc samples.
