# secp256k1 Algorithm-3 active-window certificate

## Result

`derive_active_windows.py` produces inclusive, one-based physical Work-register
windows for every repaired Luo Algorithm-3 microstep from 1 through 1616.

Artifacts:

```text
active_windows_1616.json       complete 1616-row machine table
active_windows_1477_1616.csv   compact late-schedule table
derive_active_windows.py       proof-table generator
check_active_windows.py        deterministic reconstruction and trace checks
```

The JSON uses `null` when a block is unreachable for every state in the proof
over-approximation.  Such a block can be omitted rather than instantiated with
an artificial singleton window.

Selected late rows are:

| T | R add/sub | quotient swap | T add/sub | len(t) update | len(r') update |
|---:|:---:|:---:|:---:|:---:|:---:|
| 1477 | 149..258 | 112..258 | 1..256 | - | - |
| 1480 | 149..259 | 113..257 | 1..257 | 48..258 | 224..259 |
| 1481 | 149..258 | 113..258 | 1..256 | - | - |
| 1482 | 149..259 | 114..257 | 1..256 | - | - |
| 1500 | 152..259 | 121..257 | 1..257 | 50..258 | 229..259 |
| 1524 | 155..259 | 131..257 | 1..257 | 52..258 | 235..259 |
| 1600 | 166..259 | 162..257 | 1..257 | 60..258 | 254..259 |
| 1612 | 168..259 | 167..257 | 1..257 | 61..258 | 257..259 |
| 1616 | - | - | 1..257 | 62..258 | 257..259 |

The aggregate scanned widths are 73.6% of full for R add/sub, 69.6% for the
quotient selector, 68.1% for T add/sub, 64.3% for the coefficient length scan,
and 67.8% for the remainder length scan.  These are lane-width ratios, not
claimed Toffoli ratios; exact gate savings require regenerating and counting the
corresponding blocks.

## Proof

Every `1 <= x <= p/2` has a canonical Euclidean quotient word

```text
q_0 >= 2, q_i >= 1, q_last >= 2
```

whose continuant numerator is the secp256k1 prime.  `certificate.json` proves

```text
C = sum(bit_length(q_i)) <= 404,
```

so 1616 microsteps suffice.

For a prefix of weighted cost `c` and quotient count `m`, its continuant `K`
satisfies both

```text
K >= product(q_i) >= 2^(c-m)
K >= continuant(2,1,...,1) = Fibonacci(m+2).
```

The generator minimizes the maximum of those two exact lower bounds over every
possible `m`.  For a current quotient of bit length `w`, it also uses

```text
2^(w-1) * K < p,
K < 2^c.
```

It enumerates every relaxed `(prefix cost, current quotient weight)` pair that
can contain each fixed microstep.  Relaxing away the exact equality
`continuant(word)=p` only adds states, so extrema over this set are conservative
for every real secp256k1 input.

Within a quotient of weight `w`, the exact four phase positions give:

```text
R interval:       L = ell_t + ell_q + 2, R = n + 3 - ell_s
R emitted range:  L - 1 through R        (includes carry/sign extension)
quotient selector J = ell_t + ell_q + 1  (the gate also exposes Work[J+1])
T interval:       1 through ell_t + 1
```

At a quotient boundary, prefix and suffix continuants bound the two coefficient
highest positions and the two remainder lowest positions.  The certificate
also includes the dynamic labels consumed by the range decoders themselves:

```text
B = n + 3 - bit_length(r)
A = bit_length(t_next) + 2
```

This is mandatory: a decoder endpoint can sit one or more lanes outside all
nonzero data while still being needed to toggle the range accumulator.  The
prefix bound `t+t_next <= 2^c` and Euclidean invariant bound `B`; the suffix
continuant bound together with `p < 2*r*t_next` bounds `A`.  A suffix of
remaining weighted cost `d` is below `2^d`; this is what shrinks the late
remainder scans to lanes near 259.

## Unsafe paper windows

The original Section 4.5 formulas are not safe even before they become empty.
Exact secp256k1 counterexamples reproduced by the generator include:

| T | block | required | paper | witness |
|---:|:---|:---:|:---:|:---|
| 1 | R add/sub | 2..258 | 3..259 | `x=1` |
| 8 | len(r') | 5..259 | 2..6 | `x=floor(p/2)` |
| 240 | len(r') | 41..43 | 42..64 | 1500-step witness |
| 1389 | R add/sub | 238..258 | 240..259 | 1500-step witness |
| 1470 | quotient swap | 252 | 254..258 | 1500-step witness |
| 1472 | len(t) | 245..247 | 250..259 | 1524-step witness |

The existing one-lane R widening changes the step-1389 paper window only to
239..259, so it still misses required lane 238.  The raw paper R window first
becomes empty at step 1482; quotient swap at 1484, len(t) at 1489, and len(r')
at 1493.

## Separate shift-width obligation

A 1616-step fixed schedule cannot retain the existing 9-bit `l_s`.  Every full
quotient word has weighted cost at least 256, and `x=1` terminates at exactly
1024 steps.  It therefore takes 592 terminal padding rotations.  Since

```text
592 mod 512 != 592 mod 259,
```

the wrapped 9-bit pointer no longer identifies the physical rotation of the
259-lane Work2 register.  The integration must use a 10-bit shift counter or an
equivalent exact terminal rotation counter before a 1616-step circuit can be
promoted.

## Verification

```sh
python3 derive_active_windows.py
python3 check_active_windows.py --random-cases 10000
```

The checker reconstructs the complete table, verifies the pinned schedule
certificate hash, checks all ranges, and exercises the six concrete paper
failures.  Its exact trace obligations include the hidden `A` and `B` decoder
labels.  The random trace pass is regression evidence; universality comes from
the relaxed continuant derivation above.
