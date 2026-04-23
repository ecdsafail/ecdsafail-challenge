# HMR phase-fix plan

## Situation
We now know:
- matching the HMR count was enough to repair some strict phase failures,
- but the HMR operand sequence still differs massively,
- and that residual mismatch likely explains the remaining coherent phase bug.

## Practical fix direction
The direct path to a fully phase-correct change is now:

1. **Do not keep the hand-written specialized step body** as the main design.
2. Rebuild the bulk-prefix optimization by preserving the generic measurement
   skeleton exactly.
3. Specialize only the arithmetic subpieces that do not affect the HMR / CZ /
   feed-forward ordering.

In other words, the real inverse does not just need the same classical state
map; it needs the same measurement-phase protocol. The remaining bug is now
narrow enough that trying to keep a fully hand-written alternative body is no
longer the efficient route.
