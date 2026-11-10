# Full-cube order-three decomposition: materiality and its two paths

This deposit resolves one measurement question about the 576-cell composition
record (`composition-audit/records/figure-ratio.json` and its cell archives):
how many order-three components of the measured error are material, and why the
tensor and direct-transform paths of the earlier decomposition reported
different counts.

## Inputs

The two cell archives are the record's own, unmodified:

- `composition-audit/records/cells/cells-msvc-fused.csv.xz`
- `composition-audit/records/cells/cells-msvc-separate.csv.xz`

They are boys v5.0.0, revision `482fd3388c276091834a6c12a642d9de2481c9f2`,
494,208 rows per arm, 0 refusals, as `figure-ratio.json` states.

`fullcube_summary.json` is the earlier decomposition's summary, carried here so
that the verification script can compare its tensor and direct counts against
that record.

## Commands

From this directory:

```sh
python -I verify_transform_resolution.py   # exit 0, about 56 s
python -I exact_field.py                   # exit 0, about 1 s
```

`verify_transform_resolution.py` imports `fullcube_decomposition.py` (co-located,
unmodified) and reads the cell archives; it writes `verification.json` and the
summary lines of `verify.log`. `exact_field.py` proves the borderline share in
exact arithmetic over Q(sqrt 2, sqrt 3, sqrt 6); `exact_field.log` is its
output.

## What the scripts establish

Both paths use the same definition of materiality: a component is material at a
point when its share of the point's non-mean total (the sum of squares of all
non-mean coefficients) strictly exceeds 10 per cent. The earlier note's
0/1/4 counts were a misquote; the decomposition summary's direct counts are
1/1/5/6.

- Double lane, every non-far point with a non-zero centred total (792 of 825):
  **0 of 44** order-three components are material on either route. The largest
  single-component share is 0.0408 fused and 0.0413 separate, both at
  (n=4, x=4).
- Single lane: **4 of 44** separate (largest 0.3403 at n=2, x=1.125) and
  **3 of 44** fused under a strict "exceeds" reading (largest 0.3019 at
  n=18, x=1.125). A fourth fused component at (n=2, x=4) is exactly one tenth:
  `exact_field.py` proves `10 c^2 - D0` is exactly zero there, so an "at least"
  reading gives 4.
- The extra direct-path components occur at far or bit-identical points, where
  the centred total is exactly zero and the quotient is contrast-matrix
  rounding, not data.

## Scope

The analysis is read-only with respect to the record. It says
nothing about order-three structure outside the measured cube; it only fixes
what the deposited cube supports.
