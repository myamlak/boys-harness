# The 24 changed cells re-measured at v5.3.1

The composition record was measured at kernel `v5.0.0`. The shipped fp64
region-B rational seed changed at `v5.3.1`, and its coefficients read by the
fp64, rational-minimax, shipped-partition route over region B. Those 24 cells
(2 schemes x 2 packs x 3 division forms x 2 region-B exponentials) were the
only measured cells whose arithmetic changed. Their v5.0.0 errors were stale.

This record re-measures the same 24 cells at `v5.3.1` with the same driver, the
same deposited arguments and the same deposited 60-digit reference, and
compares them cell by cell with the v5.0.0 archive.

## Result

- All 7,128 values at the 24 cells move (297 values per cell).
- The per-cell maxima fall by about a factor of three: the order-0 maximum
  falls from `1.6098653907154834e-14` to `5.6053401067827709e-15`, and the
  order-above-zero maxima fall by factors 3.0 to 3.3.
- The new maximum, `5.6053401067827709e-15`, is below the new delivered figure
  `6.04405276991530667e-15`, as the old maximum was below the old figure.
- The order-0 maximum moves from `x=12` to `x=24`; individual values can rise,
  the largest absolute increase being at `x=30.3, n=0`
  (`3.956883051108498e-17` to `5.594813019461584e-15`). The largest ratio
  anywhere is at `x=29, n=7` (905x, absolutely negligible).
- The other 487,080 rows of the 494,208-row record are bit-identical to the
  v5.0.0 archive.
- A clean v5.0.0 worktree with the same harness and toolchain reproduces the
  archived `msvc-fused` cells byte for byte.

The v5.3.1 to v5.3.3 changes in `include/` and `src/` are the version header
and comments only, so the v5.3.1 numbers carry to v5.3.3.

## Files

| file | what it is |
| --- | --- |
| `comparison-24-summary.csv` | one row per cell: v5.0.0 and v5.3.1 maxima, locations, ratios, movement counts |
| `comparison-24-detail.csv` | 7,128 rows: cell, argument, order, both values, ratio, moved flag |
| `comparison-24.json` | provenance, per-cell summary and all-row movement |
| `commands.json` | the build, run and comparison commands, with neutral placeholders |
| `summary-v531.json` | the full v5.3.1 run summary, 576/576 cells |
| `manifest-v531.sha256` | the run manifest |

The full 494,208-row v5.3.1 cell file is 157 MB and is not deposited. The
detail file carries every changed row; the v5.0.0 archive in
`../records/cells/cells-msvc-fused.csv.xz` carries the unchanged rows; the
commands and the deposited reference regenerate the full file.

## Scope

This is a re-measurement of the changed cells, not a new composition study.
The measured instance, its design and its axes are unchanged. The record is
pinned at v5.3.1; the v5.3.3 reading is the same because the later changes are
version and comments only.
