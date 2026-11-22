# The 576-call composition factorial re-measured at v5.3.3

The composition record in `../records/` was measured at kernel `v5.0.0`. It
carries 494,208 rows per configuration: 576 joint calls (four lanes crossed
with the 144-cell within-lane design) at 26 arguments and 33 orders. The
kernel's shipped fp64 region-B rational seed changed after that tag, so the
v5.0.0 errors at the affected cells are stale. This record re-measures the
whole factorial at the program-of-record kernel `v5.3.3`, with the same
driver, the same deposited arguments and the same deposited 60-digit
reference, and compares every row of both multiply-add routes with the v5.0.0
archive.

## Kernel and harness

- Kernel: boys `v5.3.3`. Annotated tag `v5.3.3`, tag object
  `4119217c8ab23d9f370b301a9bd6338e2c804f7e`, commit
  `01bba7b0ed7e0dd387f2af087e992445222e5546`; the worktree was clean.
- Harness: `boys-harness` commit
  `3be00038384f209dd0879feece817c32465a2258`, version `1.0.43`, and the
  `composition-audit` code at that commit.
- The v5.3.1 to v5.3.3 changes in `include/` and `src/` are the version
  header and comments; the arithmetic is v5.3.1's. The fused cells here hash
  `65c223da...`, the digest the v5.3.1 record names for its measured cells
  (`../v531-24cell/`).

## Build

MSVC 19.51, x64, Release, Visual Studio 18 2026. Two builds, one per
multiply-add route. The driver records the route actually in force in each
`environment.txt`; the separate build reports `separate` on both selected
lanes, so the two configurations measure different arithmetic.

```powershell
cmake -S composition-audit -B build-harness-fused `
    -G "Visual Studio 18 2026" -A x64 `
    -DBOYS_SOURCE_DIR=<boys-v5.3.3> -DBOYS_MULADD_SEPARATE=OFF
cmake --build build-harness-fused --target measure-interaction --config Release

cmake -S composition-audit -B build-harness-separate `
    -G "Visual Studio 18 2026" -A x64 `
    -DBOYS_SOURCE_DIR=<boys-v5.3.3> -DBOYS_MULADD_SEPARATE=ON
cmake --build build-harness-separate --target measure-interaction --config Release
```

## Run

The deposited argument grid and the deposited 60-digit reference are reused.
`records/arguments.csv` has 104 rows, 26 per lane (16 exact and 10 generic).
The LF checkout hashes `04eb574b...`; the CRLF working-tree copy in an older
clone hashes `3ce404b5...` and is content-identical. `records/reference.csv`
is the single 60-digit evaluation at every lane's narrowed input: 3,432
values across 104 lane-arguments, `sha256 484ff224...`. Every one of the
988,416 measured rows (two configurations) carries exactly the `ref_hi` and
`ref_lo` pair that reference states, with 0 mismatches.

```powershell
$env:BOYS_CANARY_HINT_ONLY = "1"
build-harness-fused\Release\measure-interaction.exe `
    --arguments composition-audit\records\arguments.csv `
    --reference composition-audit\records\reference.csv `
    --out results-fused\cells.csv --env results-fused\environment.txt `
    --canary-samples 1
build-harness-separate\Release\measure-interaction.exe `
    --arguments composition-audit\records\arguments.csv `
    --reference composition-audit\records\reference.csv `
    --out results-separate\cells.csv --env results-separate\environment.txt `
    --canary-samples 1
```

The driver reports 494,208 rows, 0 refused cells and 0 readings over the
published figure for each configuration.

## Canary

The driver's canary is a timing witness only; the measured error rows are
deterministic and do not depend on it. The host's canary spread across
samples was 9.36 per cent, above the driver's 5 per cent alarm, so a
multi-sample gate would have labelled the run timing-provisional and stopped
it. Each configuration was therefore run with one sample and
`BOYS_CANARY_HINT_ONLY=1`; the environment files record `canary_samples=1`
and a 0.00 per cent spread, which one sample cannot avoid.

The build is validated independently of timing. Both configurations were
re-run from the same binaries and are byte-identical to the deposited cells
(`65c223da...` fused, `7d1748cc...` separate), and 487,080 of the 494,208 rows
are bit-identical to the v5.0.0 archive. The 7,128 rows that move are exactly
the changed seed's cells.

## Result

- 494,208 rows per configuration. 7,128 values move on each route, and
  487,080 rows are identical to the v5.0.0 archive.
- The 24 affected combinations are the fp64 lane, the rational-minimax
  route, the shipped partition, both schemes, both packs, all three division
  forms and both exponentials: 2 x 2 x 3 x 2 = 24. Only the nine arguments
  whose value selects the changed fit (12, 12.7, 16, 17.3, 24, 24.7, 29,
  30.3 and 30.5) move, at all 33 orders: 24 x 9 x 33 = 7,128.
- The largest value movement is `1.1157741397482823e-14` at `n = 0`,
  `x = 12`, in the rational-minimax / split-clenshaw / shipped / orders /
  exact-division / accurate cell.
- Signed-interaction material counts at 5 / 10 / 20 per cent are, per
  resolving lane, fused fp64 `336 / 325 / 323 -> 321 / 309 / 307` and
  separate fp64 `338 / 327 / 325 -> 319 / 315 / 309`; fp32 is unchanged at
  `276 / 276 / 276`. Over both resolving lanes, fused is
  `912 / 901 / 899 -> 897 / 885 / 883` and separate is
  `914 / 903 / 901 -> 895 / 891 / 885`. The verdict remains material at all
  three margins.
- The contrast structure is unchanged: 1,808 pairwise contrasts, 1,648 of
  them informative, 664 per lane, 480 crossing the precision axis, and 44
  order-three components in the six-axis within-lane decomposition.
- The maximum error over the row's whole published statement falls on the
  double lane from `0.629268` (fused) and `0.627098` (separate) to
  `0.47303043999577227` on both. The other lanes are unchanged: fp32
  `0.819125`, fp16 `0.993403`, bf16 `0.998591`. The overall maximum remains
  `0.9985908700338076` at bf16, so no measured row exceeds its published
  statement on either route.
- The maximum absolute error on the double lane falls from
  `1.6098653907154834e-14` (fused) and `1.6043142755923576e-14` (separate) to
  `1.0337893604199906e-14` on both.

## Files

| file | what it is |
| --- | --- |
| `summary.json` | the re-measurement and comparison: design counts, changed cells, material counts, per-lane maxima and provenance |
| `results-fused/cells.csv.xz`, `results-separate/cells.csv.xz` | the 494,208 measured rows per configuration, xz-compressed with `xz -6 -T4`; the decompressed digests are `65c223da...` and `7d1748cc...` |
| `results-fused/environment.txt`, `results-separate/environment.txt` | the driver's own record of the kernel version, the route in force, the canary and the compiler |
| `commands.txt` | the configure, build, run and comparison commands with exit codes |
| `compare_v533.py` | recomputes every comparison number in `summary.json` from the deposited cells; exits 0 only when all reproduce |
| `b1_calls.tsv` | one row per joint call of the v5.0.0 cube: its largest error over the accessor figure, the row statement and six statements |
| `b1-note.md` | what the per-call ratios do and do not define |
| `b7_ratio_distribution.json` | the distribution of the interaction-to-response ratio over the 1,648 informative contrasts of the v5.0.0 cube |
| `b7-note.md` | the materiality saturation and what it means for the counts |
| `b1_calls.py`, `b7_ratio_distribution.py` | regenerate the two reductions from the deposited v5.0.0 cube |
| `provenance.json` | kernel tag, harness commit, source and output digests, commands |
| `MANIFEST.sha256` | SHA-256 of every file in this directory |

## Scope

This is a re-measurement, not a new composition study. The measured
instance, its axes and its reference are unchanged. The cells are host
measurements; this run carries no device lanes. The canary numbers are not a
timing record. The comparison is against the deposited MSVC v5.0.0 cube, and
the material counts are the signed pointwise statistic defined by
`../signed_interaction.py`.
