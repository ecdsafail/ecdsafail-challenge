# BY/compact-DIV SOTA architecture pivot

This is the current non-micro route.  Threshold tuning and local Solinas/Kaliski
polish cannot close a ~1.4M Toffoli gap.  A SOTA-shaped exact affine point-add
has to delete the two current Kaliski-sized inversion objects and replace them
with a compact in-place DIV / product-clean DIV primitive.

## Target architecture

Use a scaled Bernstein–Yang tagged-DIV microprogram instead of Kaliski:

1. Generate a 560-step branch-pattern history for denominator `x` without a
   full 560-bit denominator pair.
2. Decode each 16-bit pattern window into the A-controls.
3. Run the scaled modular microprogram on `(r,s)`:
   - tagged DIV: `(0, y+x) -> (sign(f)*q, 0)` and recover `y/x = sign(f)*q - 1`;
   - inverse/product-clean: `(sign frame q, 0) -> (0, q*x)`.
4. Delete pair1 Kaliski and pair1 schoolbook numerator multiplications; use the
   product-clean inverse shape for pair2 cleanup.

Measured/validated pieces already in the repo:

```text
scaled controlled microstep:      2,046 CCX
560-step replay:                  1,145,760 CCX
pattern+delta decoder:            ~62k CCX for all 35 windows
pattern-history replay peak:      ~1,861q with raw 560-bit pattern history
inverse product-clean replay:     same 1,145,760 CCX in -r frame
```

Budget shape (ignoring selector generation) is therefore in the Google range:

```text
non-inversion affine scaffold after deleting pair1 muls: ~0.64M
one fast BY DIV replay + decode margin:                  ~1.30M
projected fast tagged-DIV point-add:                     ~1.9M-2.3M
low-scratch/vented replay variant:                       ~2.65M
```

So the blocker is no longer the modular replay body.  It is the selector / branch-history generator.

## New selector finding: streaming limb generator

A pure `h = g/f mod 2^16` state is not forward-complete; next-window `h` needs
higher 2-adic denominator data.  A more exact streaming representation is:

```text
(f_j, g_j) = A_j · (p >> 16j, x >> 16j) + c_j
```

where `A_j` is the product of prior 16-step BY matrices and `c_j` is the
low-limb carry.  The test
`streaming_limb_selector_is_exact_but_state_heavy` proves this representation
can generate all 35 branch windows exactly on sampled secp256k1 denominators
without carrying the full denominator pair.

Result:

```text
8 limbs per A/c entry:  fails on all 64 sampled denominators
12 limbs per A/c entry: exact on 64 sampled denominators
state = 4*A entries + 2*c entries = 6 * 12 * 16 = 1152 bits
```

A first structural compression folds the constant `p` column into the carry,
leaving only the two coefficients of the quantum `x` tail plus two carry rows:

```text
16 limbs per folded entry: fails on all 64 sampled denominators
17 limbs per folded entry: exact on 64 sampled denominators
state = 2*b_x entries + 2*c entries = 4 * 17 * 16 = 1088 bits
```

This is constructive but not yet SOTA-shaped.  It replaces a full 560-bit
2-adic denominator pair with a selector state, and the constant-column fold
shows the state can move in the right direction, but 1088 selector bits is still
too large for the ~600 extra-qubit target.

## Next architecture work

Do not spend more time on local Kaliski thresholds.  The next useful work is to
compress or avoid the `A_j` state:

- Recompute `A_j mod 2^k` reversibly from compressed pattern history only when a
  window needs it, then uncompute it before modular replay.
- Factor `A_j` using the already-observed small Hermite/unimodular window
  structure so the live selector state is `O(2 rows + carry)` instead of four
  12-limb matrix entries.
- Try entropy-coded matrix/pattern history as the persistent object, with a
  small rolling carry sidecar; target is `<480` history bits + `<150` selector
  scratch.
- If those fail, the scaled BY route is Toffoli-viable but not low-qubit viable;
  pivot to a compact Montgomery/windowed DIV core rather than returning to
  current Kaliski polishing.
