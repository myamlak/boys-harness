# Taxonomy demonstration

A synthetic testbed that reproduces the two demonstrations in the paper's
section on what a row can and cannot carry, and its Table S11. It is
constructed, not found: it models two evaluation paths and the row types that
judge them, then plants a small drift in one path to show which rows can see
it.

## The model

- The reference is the Boys function `F0(x)` on `[0, 11.899848152108484]`,
  evaluated from `math.erf`.
- The scalar and fast paths share a baseline error of `9e-14 * F0(x)` on
  `[9, 10]`. Both paths agree everywhere before the plant, and the baseline
  error is inside the lane accuracy bound.
- The lane accuracy row bounds the worst error of either path against the
  reference by `1e-13 * F0(x)`.
- The band row bounds the fast path on `[2, 3]` by `1e-15 * F0(x)`.
- The consistency row bounds `|fast - scalar|` by the sum of the two paths'
  half-lane bounds, `2e-15 * F0(x)`.
- The exactness row checks `F0(0) == 1` on both paths.
- The ordinary test checks only the lane accuracy bound.

The plant adds `5e-15 * F0(x)` to the fast path on `[2, 3]`. It is inside the
lane accuracy bound, above the band bound and above the consistency
allowance.

## Expected result

Before the plant every row passes. After it:

| row | before | after |
| --- | --- | --- |
| lane accuracy (worst error) | 9.009803e-14, pass | 9.009803e-14, pass |
| band bound | pass | fail at 5.104755e-15 |
| consistency | pass | fail at 5.104755e-15 |
| exactness | pass | pass |
| ordinary test | pass | pass |

The lane accuracy row's worst error is unchanged by the plant: the baseline
error on `[9, 10]` is larger than the drift, so the row cannot see it. The
consistency row and the band row can, because they compare a path with a path
and a band with its own bound.

## Running

```sh
python -I taxonomy_demo.py --json taxonomy_demo.json
```

The script prints the before and after rows and exits 0 when every verdict
matches the table above. `run.log` and `taxonomy_demo.json` are the captured
output of that command; `provenance.json` records the script's digest, the
Python version, the command and the exit status.

## Scope

The testbed is a construction, not a measurement of the released library. It
shows that the property type decides which defects a row can see; it says
nothing about how often a real codebase contains such a defect. The audit's
original scripts were not preserved, and this testbed reproduces the
construction the paper describes rather than the audit's original run.
