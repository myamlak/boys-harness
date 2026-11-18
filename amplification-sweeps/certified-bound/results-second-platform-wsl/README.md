# Second-platform replay of the certified-bound record

Data deposit, not a release. Parent: `boys-harness` v1.0.42. The annotated
tag `second-platform-2026-11-18` names this commit. The run re-executes the
certified-bound crossover on a WSL2/GCC/glibc platform against the same
`boys` v5.3.2 kernel and the same committed inputs.

## Platform

- WSL2, Ubuntu 26.04 LTS, kernel 6.18.33.2-microsoft-standard-WSL2, x86_64.
- Intel Core i7-9850H, 12 logical processors, 15 GiB RAM.
- GCC/G++ 15.2.0 (Ubuntu 15.2.0-16ubuntu1), binutils 2.46.
- glibc 2.43 (Ubuntu GLIBC 2.43-2ubuntu2.4); libm
  `/lib/x86_64-linux-gnu/libm.so.6` from libc6 2.43-2ubuntu2.4.
- Python 3.14.4, mpmath 1.4.1.
- Kernel: detached checkout at tag `v5.3.2`, commit
  `6739447d82765ae5f7f68351a99f418cd018dd9a`, clean.
- The sweep is header-only at this API: it compiles against
  `include/boys` and links no library.

## What was run

The full 55,176-entry crossover, not only the accuracy gate. The commands and
exit codes are in `logs/commands.txt`; the environment capture is
`logs/00-environment.txt`. The build line is

```
g++ -std=c++23 -O2 -ffp-contract=off -fno-fast-math -I <v5.3.2>/include
```

for both the exponential probe and the sweep. The comparison baseline is the
deposited repaired record under `results-underflow-fix/`, which is the record
the paper's numbers use.

## Result

- 55,176 entry-cells; 132 band-metric comparisons; 65 commitments; 49
  truncation-agreeing and 16 truncation-opposite; 0 bound violations; worst
  `|r|/C = 0.71875` at `(A1, 0.0, 24, ship)`; 264/264 population band
  intervals inside.
- 70 of 55,176 delivered doubles differ, all on the recursed route; the
  27,588 direct entries are bit-identical. They sit at five grid arguments
  where glibc `exp(-x)` differs from UCRT by one ulp: A1 `1e-08`,
  `1.5590366116438534`, `1.9917591617705526`; A2 `7.0276408694349346`,
  `11.360263861161261`.
- `band_cert.csv`, `model_cells.csv` and `harness_cells.csv` are byte-identical
  to the deposited record. `certify.log` and `harness.log` are text-identical
  after removing the run's own header and exit lines; see
  `logs/13-log-summary-compare.log` and `logs/11-record-compare.log`.
- Loser certified charge/truncation medians: 83.289559387479912 over the 59
  recursed losers and 13.309575283201504 over the six direct losers. Required
  factor over the 49 agreeing comparisons: median 85.225820155890801,
  minimum 8.9998413453301787, maximum 3937.6896051746826. Charge violations:
  certified 0, paper 5.
- Exponential probe: worst `0.988081126378 x 2^-53` at
  `x = 11.084208097312864`; 2,582 arguments above 0.5 unit and none above one;
  Python `math.exp` mismatches 0. The deposited UCRT control gives
  `0.999662996081 x 2^-53`.

## Files

- `band_cert.csv`, `cells_cert.csv`, `harness_cells.csv`, `model_cells.csv`,
  `summary.json` — the certificate records.
- `exp/exp_bits_second_platform.txt` — the 9,541-argument exponential probe.
- `sweep/cells.csv`, `sweep/admission.csv`, `sweep/per_order.csv`,
  `sweep/summary.txt` — the raw sweep records.
- `logs/` — environment, commands, comparisons, statistics and digests.
- `tools/` — the comparison and statistics scripts used here.
- `SHA256SUMS` — digests of every file above.

The 2.8 GB `paths.csv` sweep intermediate is not deposited; it is reproducible
from the command in `logs/commands.txt`. No paper, harness or kernel file was
edited by the run.
