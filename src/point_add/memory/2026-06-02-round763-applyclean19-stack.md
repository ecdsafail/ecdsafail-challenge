# 2026-06-02 round763 + apply-clean-19 stack

Route: promoted chunked apply-F base (`3bd3aa5` / `550895c`) plus two older
Toffoli-only wins that were not present in that head.

Changes stacked:

- Re-enabled the round763 reachable-support compressor lever with
  `DIALOG_GCD_ROUND763_COMPRESS_LEVER=1`.
- Tightened the apply clean boundary comparator from
  `DIALOG_GCD_APPLY_CLEAN_COMPARE_BITS=20` to `19`.
- Kept the promoted head's `DIALOG_GCD_COMPARE_BITS=61`,
  `DIALOG_GCD_APPLY_CHUNKED_F_BLOCKS=2`, and
  `DIALOG_GCD_APPLY_CHUNKED_F_CUT=70`.

Old islands did not transplant:

- `4/15` (promoted chunked route) failed phase garbage.
- `6/3` (old round763 lever route) failed on the stacked route.
- `0/23` (old apply-clean-19 route) failed on the stacked route.
- `1/0` and `3/18` also failed.

Validated island:

- `DIALOG_REROLL=6`
- `DIALOG_POST_SUB_REROLL=28`
- Trace peak: `1567` qubits at
  `round84_fused_square_xtail_dx_sub_lam_square_lowq`.
- Eval: all `9024` shots OK, `0` classical mismatches,
  `0` phase-garbage batches, `0` ancilla-garbage batches.
- Metrics: `1,685,515.000` average executed Toffoli, score
  `2,640,201,005`.

The stack is peak-neutral versus the promoted chunked route and saves `3,990`
average executed Toffoli.
