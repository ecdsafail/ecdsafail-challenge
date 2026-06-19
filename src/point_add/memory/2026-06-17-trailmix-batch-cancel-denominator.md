# TrailMix shrunken-PZ cancel batching check

Context: the sub-1000 TrailMix route pays two shrunken-PZ divide-style passes:
`ec3.inv_fwd` computes `lambda = dy/dx`, and `ec3.alt.cancel` recomputes
`lambda = new_dy/new_dx` to erase the live slope ancilla after the output
coordinates have been formed.

## Algebra result

The second denominator is algebraically available before the first inversion, but
only after clearing the hidden `dx^-2` denominator:

```text
dx     = Qx - Px
dy     = Qy - Py
lambda = dy / dx
Rx     = lambda^2 - Px - Qx
new_dx = Qx - Rx

new_dx = 3*Qx - dx - dy^2/dx^2
D      = (3*Qx - dx)*dx^2 - dy^2
D      = new_dx * dx^2
```

So a Montgomery-style batch can invert `dx` and `new_dx` by inverting
`dx*D` once:

```text
(dx*D)^-1       = batch_inv
dx^-1           = D * batch_inv
new_dx^-1       = dx^2 * D^-1 = dx^3 * batch_inv
new_dy/new_dx   = new_dy * dx^3 * batch_inv = lambda
```

This was checked with a standalone Python field-arithmetic script on 64
deterministic secp256k1 generic additions. A Rust `#[cfg(test)]` guard was added
in `trailmix_port/ec/point_add.rs`, but the repository-wide test target is
currently blocked by pre-existing test-only compile errors unrelated to this
check.

## Low-qubit feasibility result

The algebra is not enough for the current q988 architecture. The active-qubit
trace shows the shrunken-PZ inversion peaks at 988q in both `ec3.inv_fwd` and
`ec3.alt.cancel`. Non-EEA arithmetic phases peak around 789q, so temporary field
registers outside the inversion are affordable, but any extra full 257-bit
cofactor/witness crossing the inversion peak breaks the sub-1000 target.

The straightforward batch needs enough live material across the single inversion
to recover both:

- `lambda = dy * D * (dx*D)^-1`
- `new_dx^-1 = dx^3 * (dx*D)^-1`

That means carrying `D`/`dy*D` and `dx^2`/`dx^3`-class material in addition to the
single denominator/inverse state, or recomputing it later. Carrying it costs at
least one full field register at the 988q inversion peak; recomputing it costs
multiple standalone modular multiplies and still needs a clean way to restore
the overwritten `dx`/`dy` information.

## Decision

For the current two-register sub-1000 TrailMix layout, "batch away the second
cancel inversion" is algebraically valid but not yet a viable circuit cut. It
would need a new representation that transforms the two existing field registers
into an invertible `(denominator, passenger)` pair while retaining enough
cofactor information to derive both slope consumers. Without that representation,
the second inversion is the clean low-qubit way to recover the inverse-add
relation from `(new_dx, new_dy)`.

The next practical Toffoli work should therefore attack the cost of the cancel
pass itself: phase-specific shrunken-PZ schedules, cheaper bitlen/shift/compare
windows for the `new_dx` distribution, or a genuinely new two-register fused
divide representation that keeps the multiply consumers fused instead of
splitting into standalone modular multiplications.
