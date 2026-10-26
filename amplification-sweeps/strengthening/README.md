# Paper A strengthening: scripts, inputs and run record

This directory is the deposit copy of the analysis behind `STRENGTHEN.md`, the strengthening
report for Paper A. It holds the ten analysis scripts, the thirteen logs and the twenty out CSVs
that the report cites, in the form in which they were produced. `STRENGTHEN.md` itself lives in
the paper tree, not in this repository; this README is the index that connects its numbers to the
record here.

In the report, `as/` is this directory and `as/raw/` is the generated `inputs/` directory
described below. The report's `deposit/` is the parent of this directory.

## Requirements

- Python 3.14.7 (the version used for the record).
- `mpmath` 1.4.1. On the recording machine it lives in the user site; run Python **without**
  `-I`, or the import fails and several scripts exit 77.
- A Windows host if you want byte-identical logs. The CSVs are CRLF everywhere because the csv
  module writes them; Python's redirected stdout is CRLF on Windows and LF elsewhere. The
  committed logs were produced on Windows.

Commands below assume the current directory is `amplification-sweeps/strengthening`, and that
`python` is the 3.14 interpreter with `mpmath` importable.

## Inputs

### Compressed streams, decompressed into `inputs/`

`prepare_inputs.py` streams these three deposit files into `inputs/` and checks each decompressed
SHA-256. The `.xz` files are opened read-only and are never modified. `inputs/` is generated, is
not committed, and is ignored by `.gitignore`.

| Deposit stream | Packed bytes | Becomes | Decompressed bytes | SHA-256 of the decompressed stream |
| --- | ---: | --- | ---: | --- |
| `product-amplification/results/cells.csv.xz` | 3,418,908 | `inputs/cells.csv` | 21,670,457 | `dd26d0b2a4352113d9f15433ec991adfacc211607ae9586c0bb51cd12b3f3369` |
| `product-amplification/data/reference.csv.xz` | 8,590,988 | `inputs/reference.csv` | 20,586,094 | `1e554a3540c6a15535b80f7f2597360bfc3a75919fa71d1496b284f348deb5b9` |
| `amplification-form/results/paths.csv.xz` | 31,971,200 | `inputs/paths.csv` | 291,827,555 | `81001a06d80c88e5b3ddac1f3dd011cae6f8ce2e806c8bf1f7f26e85dd8c8219` |

The first two are the first sweep's measured cells and its 60-digit reference; the third is the
3,004,442-row chain record used by the transfer checks. Of the deposit's two campaign trees, only
these three files are compressed: `amplification-form/results/cells.csv` and
`amplification-form/data/reference.csv` are plain files read in place.

`probe_transfer.py` also writes `inputs/probe_cell.csv`, a 33-line extract of its chosen cell's
chain (header, the seed row, and orders 3 to 32), on first use. It is generated too.

### Deposit records read in place

No copy of these is made. They are read from `amplification-form/` and
`product-amplification/` under the deposit root (`HERE.parent`):

- `amplification-form/analyze.py` — the model arithmetic, imported rather than retyped.
- `amplification-form/results/cells.csv`, `data/reference.csv` — second sweep, read in place.
- `amplification-form/results/per_order.csv`, `band_summary.csv`, `prediction.csv`,
  `prediction_scan.csv` — the second sweep's deposited results, used for the zero-mismatch
  reproduction checks.
- `amplification-form/results/summary.txt` — positive controls cited in Sections 4.2, 4.4 and
  4.5 of the report.
- `amplification-form/results/transfer_cells.csv` — the deposit's transfer sample used by
  `check_transfers.py`.
- `amplification-form/results/provenance.json` — thresholds and build defaults.
- `product-amplification/results/per_order.csv`, `band_summary.csv`, `provenance.json` — the
  first sweep's deposited results and provenance.
- `surrogate-control/surrogate_control.py` — imported by `intervals.py`.
- `common/` — on `sys.path` for the model imports.

## Commands in order

Run from this directory. The block reproduces the record's logs exactly, except
`verify_provenance.log` (see the caveats). The exit-code lines are part of the record: where the
logged command has one, append it exactly as shown; where the log ends with the script's own last
line, do not append anything.

```sh
python prepare_inputs.py

python score_sweeps.py --sweep first  --out out/first     --margin 1 > out_first.log 2>&1
echo EXIT=$? >> out_first.log

python score_sweeps.py --sweep second --out out/second    --margin 1 > out_second.log 2>&1
echo EXIT=$? >> out_second.log

python score_sweeps.py --sweep first  --out out/first_m0  --margin 0 > out_first_m0.log 2>&1

python score_sweeps.py --sweep second --out out/second_m0 --margin 0 > out_second_m0.log 2>&1

python compare_sweeps.py > compare_sweeps.log 2>&1
echo EXIT=$? >> compare_sweeps.log

python same_fit_region.py > same_fit_region.log 2>&1
echo EXIT=$? >> same_fit_region.log

python intervals.py > intervals.log 2>&1
echo EXIT=$? >> intervals.log

python check_transfers.py > check_transfers.log 2>&1

python probe_transfer.py > probe_transfer.log 2>&1
echo EXIT=$? >> probe_transfer.log

python probe_start_independence.py > probe_start_independence.log 2>&1
echo "EXIT=$? (observed)" >> probe_start_independence.log

python diagnose_start_independence.py > diagnose_start_independence.log 2>&1

python population_winners.py > population_winners.log 2>&1

python verify_provenance.py > verify_provenance.log 2>&1
echo EXIT=$? >> verify_provenance.log
```

Order matters in one direction only: the four `score_sweeps.py` runs must come before
`compare_sweeps.py`, `intervals.py` and `population_winners.py`, which read their out CSVs.
`same_fit_region.py`, `check_transfers.py`, the two probe scripts and the diagnosis read only the
inputs and the deposit. `prepare_inputs.py` must come first.

Re-running overwrites the committed record files in place. They are committed, so recover them
from version control or copy them aside before a rerun.

### Exit-code and byte conventions

- Python's redirected stdout writes CRLF on Windows. A shell `echo` appended to a log writes a
  bare LF, so a log with an appended exit line has one LF-only line at the end. That is the
  committed byte pattern and the reason the exit line is shown separately above.
- Five logs have no exit line: `out_first_m0.log` (7 lines), `out_second_m0.log` (9),
  `check_transfers.log` (63), `diagnose_start_independence.log` (35) and
  `population_winners.log` (62). Their exit code 0 was observed at the runner; it is not in the
  log. `probe_start_independence.log` ends with `EXIT=0 (observed)` for the same reason.
- `out_second_m0.log:7` reports 33 mismatches against the deposited `prediction.csv`. That is
  expected: the margin-0 run is being checked against a margin-1 deposited record, and the
  margin-0 model commits where the record abstains. It is not a failure.
- `verify_provenance.py` is the one script whose output is not reproduced verbatim from this
  deposit. See the caveats at the end.

## What each command writes

| Command | Writes | Log (lines in the record) |
| --- | --- | --- |
| `prepare_inputs.py` | `inputs/cells.csv`, `inputs/reference.csv`, `inputs/paths.csv` | none in the record; prints three size/hash lines and a closing line |
| `score_sweeps.py` first, margin 1 | `out/first/{prediction,scan,scan_all,cell_scan,per_order_ratio}.csv` | `out_first.log` (8) |
| `score_sweeps.py` second, margin 1 | `out/second/{prediction,scan,scan_all,cell_scan,per_order_ratio}.csv` | `out_second.log` (10) |
| `score_sweeps.py` first, margin 0 | `out/first_m0/{prediction,scan,scan_all,cell_scan,per_order_ratio}.csv` | `out_first_m0.log` (7, no exit line) |
| `score_sweeps.py` second, margin 0 | `out/second_m0/{prediction,scan,scan_all,cell_scan,per_order_ratio}.csv` | `out_second_m0.log` (9, no exit line) |
| `compare_sweeps.py` | none | `compare_sweeps.log` (56) |
| `same_fit_region.py` | none | `same_fit_region.log` (217) |
| `intervals.py` | none | `intervals.log` (41) |
| `check_transfers.py` | none | `check_transfers.log` (63, no exit line) |
| `probe_transfer.py` | `inputs/probe_cell.csv` | `probe_transfer.log` (28) |
| `probe_start_independence.py` | none | `probe_start_independence.log` (30) |
| `diagnose_start_independence.py` | none | `diagnose_start_independence.log` (35, no exit line) |
| `population_winners.py` | none | `population_winners.log` (62, no exit line) |
| `verify_provenance.py` | none | `verify_provenance.log` (18 in the record; the deposit-only rerun is shorter) |

Only `prediction.csv` depends on `--margin`. The other four CSVs are the same file for margin 0
and margin 1: `scan.csv`, `scan_all.csv`, `cell_scan.csv` and `per_order_ratio.csv` are identical
between `out/first` and `out/first_m0`, and between `out/second` and `out/second_m0` (the digest
table below shows the repeated hashes). The margin-0 runs exist to produce the margin-0
`prediction.csv` rows used for the tie counts.

## Where each number in STRENGTHEN.md comes from

Line numbers are 1-based and count the header of a CSV as line 1. Each row covers one number or
one table; the cited log range contains every decimal in that table. Where the same values also
appear in an out CSV, the CSV range is given. "Derived" rows are produced by the count recipe in
the next section, not printed by a script.

### Section 0, provenance and commands

| Report | Number or claim | Log | Out CSV |
| --- | --- | --- | --- |
| 0.2 | three inputs byte-identical to the `.xz` streams, with their hashes | `verify_provenance.log:1-3` (original run; labels `raw/`). The deposit-only rerun prints the same three hashes with `inputs/` labels. | - |
| 0.2 | the three model function bodies identical; paper-tree and deposit `analyze.py` hashes | `verify_provenance.log:4-8` | - |
| 0.2 | same boys commit and tag, both trees clean | `verify_provenance.log:9-11` | - |
| 0.2 | compiler, flags, 15 boys digests, machine | `verify_provenance.log:12-15` | - |
| 0.2 | Python 3.14.7 with mpmath 1.4.1 | no log; runner versions, see Requirements | - |
| 0.3 | theta0 = 1.0855252345349333 | `out_first.log:3`, `out_second.log:3`, `out_first_m0.log:3`, `out_second_m0.log:3` (`model_from`) | - |
| 0.3 | theta17 = 10.783587858916762 and B1's upper edge 5.949924076054242 | no analysis log; `amplification-form/results/provenance.json`, `constants.kTierThresholds[17]` and `constants.kRegionA1Edge` | - |
| 0.3 | build defaults Horner / kUniform | no analysis log; `amplification-form/results/provenance.json` `build_defaults` | - |
| 0.4 | commands and exit codes | the Evidence column of the report's Section 0.4 command table: `out_first.log:8`, `out_second.log:10`, `compare_sweeps.log:56`, `same_fit_region.log:217`, `intervals.log:41`, `probe_transfer.log:28`, `probe_start_independence.log:30`, `verify_provenance.log:18`; the five no-exit logs are listed above | - |
| 0.4 | 264 values and 0 mismatches, 0 band rows, 132 prediction rows, 28 scan rows | `out_second.log:4-8` | `out/second/per_order_ratio.csv`, `out/second/scan.csv`, `out/second/prediction.csv`, `out/second/scan_all.csv` |

### Section 1, out-of-sample scoring

| Report | Number or claim | Log | Out CSV |
| --- | --- | --- | --- |
| 1.1 | 9,352 arguments, 308,616 cells, 4,812 A1 and 4,540 A2 | `out_first.log:2` (totals and cells); 4,812 and 4,540 follow from `compare_sweeps.log:4` (shared + first-only per band) | - |
| 1.2 | the seven-margin ladder, first and second sweeps | - | `out/first/scan_all.csv:2-8`; `out/second/scan_all.csv:2-8` |
| 1.2 | the second margin-1 row reproduces the deposited record | `out_second.log:7-8` | `out/second/scan_all.csv:5`; `out/second/scan.csv:14-17` |
| 1.3 | margin-1 commits/right/wrong/abstains by band and metric | - | `out/first/scan.csv:14-17`; `out/second/scan.csv:14-17` |
| 1.3 | margin-0 commits/right/wrong/ties by band and metric | - | `out/first/scan.csv:2-5`; `out/second/scan.csv:2-5` |
| 1.3 | the one first-sweep tie is A2 absolute order 18 | - | `out/first/scan.csv:4` (`ties_measured` 1); `out/first_m0/prediction.csv:86` |
| 1.3 | the three second-sweep ties are A2 absolute orders 17, 18, 19 | - | `out/second/scan.csv:4` (`ties_measured` 3); `out/second_m0/prediction.csv:85-87` |
| 1.4 | predicted recursed and direct, arithmetic 0, choice product, whole-band winner ship | `population_winners.log:4-6` | `out/first/prediction.csv:2` (predicted values, state, choice); `out/first/per_order_ratio.csv:2` (measured whole-band worsts and winner) |
| 1.4 | measured whole-band recursed 1.6610920772963861e-16 at x=0.6623, direct 2.2079726518177838e-16 at x=0.0603 | `population_winners.log:5` | `out/first/per_order_ratio.csv:2` (values, no x) |
| 1.4 | population recursed 1.3755401165396436e-16 and direct 1.3116554139209733e-16, with x | `population_winners.log:6` | - (population rescore is not written to an out CSV) |
| 1.4 | 69 commits, 69 right on population-restricted winners | `population_winners.log:47-48` | - |
| 1.4 | margin 0: 87 right, 11 wrong, 1 tie against 85/14 whole-band | `population_winners.log:49-50,62` | - |
| 1.4 | 14 of 132 winners change, A1 absolute orders 0, 10, 17-19, 21, 23-27, 29-31, with values and x | `population_winners.log:3-45` | - |
| 1.5 | first-sweep winners 42 recursed, 89 direct, 1 tie | derived | `out/first/per_order_ratio.csv:2-133` (`measured_winner`) |
| 1.5 | second-sweep winners 44, 85, 3 | derived | `out/second/per_order_ratio.csv:2-133` (`measured_winner`) |
| 1.5 | always-direct right at 60 of 69 first-sweep commitments; 59 of 67 second-sweep | derived | `out/first/prediction.csv` + `out/first/per_order_ratio.csv`; `out/second/prediction.csv` + `out/second/per_order_ratio.csv` |
| 1.5 | the model is right at 68 of 69 and 67 of 67 | - | `out/first/scan_all.csv:5`; `out/second/scan_all.csv:5`; `population_winners.log:47` |
| 1.6 | 265,287 first-sweep population cells, 18,876 second-sweep | `out_first.log:6` (sum of the two band cell counts); `out_second.log:6` | `out/first/cell_scan.csv:2-15` (`cells`); `out/second/cell_scan.csv:2-15` (`cells`) |
| 1.6 | the cell table at margins 1, 2 and 8 | - | `out/first/cell_scan.csv:8-11,14-15`; `out/second/cell_scan.csv:8-11,14-15` |
| 1.6 | 920 wrong of 4,577 at margin 8 | - | `out/first/cell_scan.csv:15` |
| 1.7 | 8,704 first-only and 188 second-only arguments; 648 shared (460 A1, 188 A2) | `compare_sweeps.log:3-4` | - |
| 1.7 | 359 of 648 digit-identical references; worst 9.98147e-60 | `compare_sweeps.log:5-6` | - |
| 1.7 | 21,384 shared cells, all bit-identical | `compare_sweeps.log:7` | - |
| 1.7 | 411 shared population arguments; 8,039 first and 572 second | `compare_sweeps.log:8` | - |
| 1.7 | eight measured winners differ, listed with their ratios | `compare_sweeps.log:9-17` | - |
| 1.7 | 648 of 836 and 648 of 9,352 | `compare_sweeps.log:3` | - |

### Section 2, the same-fit region

| Report | Number or claim | Log | Out CSV |
| --- | --- | --- | --- |
| 2.1 | the 16-order first-sweep table: recursed/direct absolute and relative, whole-band and population ratio | `same_fit_region.log:6-37` (absolute and relative alternate, order 17 first) | whole-band ratios also `out/first/per_order_ratio.csv:36-67` |
| 2.1 | 15 consecutive orders 17-31 with ratios 0.343904 to 0.545983, factors 1.83 to 2.91; order 32 is 0.950854; relative 1.239357 to 2.862334 | `same_fit_region.log:6-37` | `out/first/per_order_ratio.csv:36-67` |
| 2.2 | absolute population ratio above 1 at 12 orders, 1.050699 to 1.551652; below 1 at 20, 22, 28, 32 | `same_fit_region.log:6-36` (population ratio column) | - |
| 2.2 | relative maxima at x about 5.88 to 5.95 | `same_fit_region.log:39-103` | - |
| 2.2 | whole-band and population absolute ratios differ at all 16 orders | `same_fit_region.log:6-36` (two ratio columns) | - |
| 2.2 | order 24 and order 26 positions | `population_winners.log:25-27,31-33` | - |
| 2.3 | the 16-order second-sweep table | `same_fit_region.log:114-145` | whole-band ratios also `out/second/per_order_ratio.csv:36-67` |
| 2.3 | 10 above 1, 3 tied (22, 27, 29), 3 below (18, 20, 21); relative below at 17 and 32 | `same_fit_region.log:114-145` | - |
| 2.4 | first sweep: 55,984 cells, 8,612 bit-identical, 8,052 rec correctly rounded, 7,790 direct, quantiles, fraction 0.4733 | `same_fit_region.log:105-108` | - |
| 2.4 | second sweep: 4,640 cells, 1,094, 906, 896, quantiles, fraction 0.4562 | `same_fit_region.log:213-216` | - |
| 2.4 | 48,194 and 3,744 cells with non-zero direct error | `same_fit_region.log:107,215` | - |
| 2.5 | theta17, B1 edge, and the build-default note | no analysis log; `amplification-form/results/provenance.json` (constants, build defaults); kernel source lines are in the report | - |

### Section 3, intervals

| Report | Number or claim | Log | Out CSV |
| --- | --- | --- | --- |
| 3.1 | all intervals | `intervals.log:13-40` | the counts come from `out/second/prediction.csv` and `out/first/prediction.csv` |
| 3.2 | 67/67: Wilson [0.9458, 1.0000], Clopper-Pearson [0.9464, 1.0000] | `intervals.log:17` | `out/second/scan_all.csv:5` |
| 3.2 | 17/100: [0.1089, 0.2555], [0.1023, 0.2582]; 14/100: [0.0853, 0.2214], [0.0787, 0.2237] | `intervals.log:18-19` | `out/second/scan_all.csv:2`; `out/second/scan.csv:2-5` |
| 3.2 | 46/67: [0.5680, 0.7849], [0.5616, 0.7944]; bootstrap [0.5231, 0.8361]; 22 of 33 orders | `intervals.log:20,32-33` | - |
| 3.2 | 68/69: [0.9224, 0.9974], [0.9219, 0.9996] | `intervals.log:21` | `out/first/scan_all.csv:5` |
| 3.2 | six-decision gap: 6.000 +/- 2.250, bootstrap [2, 10] | `intervals.log:38-40` | - |
| 3.3 | 100 commitments and 17 not right at margin 0; 3 ties in A2 absolute | `intervals.log:18-19` | `out/second/scan_all.csv:2`; `out/second/scan.csv:2-5`; `out/second_m0/prediction.csv:85-87` |
| 3.3 | first-sweep counterpart: 99 commitments, 14 not right, 13 wrong and 1 tie at A2 absolute order 18 | `population_winners.log:49,62` | `out/first/scan_all.csv:2`; `out/first/scan.csv:2-5`; `out/first_m0/prediction.csv:86` |
| 3.4 | jackknife 0.1700 +/- 0.0418 against binomial 0.0376; bootstrap [0.0943, 0.2604] | `intervals.log:27-28` | - |
| 3.4 | 67/67 degenerate block interval [1.0000, 1.0000]; Clopper-Pearson lower 0.9464 | `intervals.log:17,34` | - |
| 3.4 | six extra commitments at orders 3, 4, 11, 14, 21, 22 | `intervals.log:37` | - |
| 3.4 | margin-0 wrong counts by order; 13 of 33 orders | `intervals.log:24-25` | - |

### Section 4, the record's unexplained checks

| Report | Number or claim | Log | Out CSV |
| --- | --- | --- | --- |
| 4.1 | 3,004,442 paths rows | `check_transfers.log:1` | - |
| 4.1 | the `transfer_ulp` header | `probe_start_independence.log:1` | - |
| 4.1 | the two probe replays are bit-exact | `probe_transfer.log:4-7`; `probe_start_independence.log:1-3` | - |
| 4.2 | 9,505 sampled cells; worst deviation 5.570598620124316; 232 cells above 2 | `check_transfers.log:4-6` | - |
| 4.2 | distribution median 0.208, p90 0.947, p99 2.75 | `check_transfers.log:28` | - |
| 4.2 | worst cell A2 x=11.896188857019089, start 2, degree 2, n 25; measured 2.8125; 6.57x the closed form 0.42804319097896554 | `check_transfers.log:10,22-23`; `probe_transfer.log:1-2,25` | - |
| 4.2 | replay: 0 base and 0 perturbation deltas differ; both transfers 2.8125 | `probe_transfer.log:4-7` | - |
| 4.2 | the k-ULP series from 1 to 16,777,216 | `probe_transfer.log:9-22` | - |
| 4.2 | canonical transfer 4.0723809634881496e-15 (A1) and 2.9303199895784525e-14 (A2) | no analysis log; `amplification-form/results/summary.txt:22-27` | - |
| 4.3 | 132,496 sampled cells; worst 1.398619401680519e183; 34,190 unmeasurable; all 100 worst unmeasurable | `check_transfers.log:33-35,51` | - |
| 4.3 | the 1.4e183 cell: A1 x=3.1622776601683792e-06, start 28, degree 24, n 0; closed form 3.08995e-183, seed -3.21e-18, ratio 2.8e-186 | `check_transfers.log:38` | - |
| 4.3 | worst measurable 264.029 at A2 x=10.910756753077845, start 2, degree 20, n 26 | `check_transfers.log:49`; `diagnose_start_independence.log:6,20-28` | - |
| 4.3 | seed error, numerator, closed form, measured transfer, seed term, charge, charge ratios, delivered relative error | `probe_start_independence.log:5-17`; delivered relative error `diagnose_start_independence.log:27` | - |
| 4.3 | the four-regime diagnosis table: 61,666/0, 9,903/1,913, 25,001/426, 1,736/371 | `diagnose_start_independence.log:15-18` | - |
| 4.3 | 13,899 of 26,737 measurable cells below theta0 violate the bound | `diagnose_start_independence.log:7` | - |
| 4.4 | 26,648 cells; p=-8.892238833209769, q=3.028657826313866; max residual 687.0537228913399; rms 55.64197937662048 | `check_transfers.log:54-56` | - |
| 4.4 | log-transfer range -11.5563 to 690.082; n log x range -690.776 to 79.2488 | `check_transfers.log:57-58` | - |
| 4.4 | per-argument refit: 833 arguments, worst 4.12561 | `check_transfers.log:59-60` | - |
| 4.4 | three-parameter refit: c3=-0.9999317383958976, max 4.1257, rms 0.988478 | `check_transfers.log:61-63` | - |
| 4.4 | the deposit's derived form: 2.930e-14 maximum, 2.679e-15 rms | no analysis log; `amplification-form/results/summary.txt:34-35` | - |
| 4.5 | 61,666 seed-dominated population cells; 17,912-cell admission check (8,679 down, 33 exact, 9,200 up) | `diagnose_start_independence.log:15`; admission check `amplification-form/results/summary.txt:17-27` | - |

### Section 5, not settled and incidental

| Report | Number or claim | Log | Out CSV |
| --- | --- | --- | --- |
| 5.1 | replay and provenance caveats | `probe_transfer.log:4-7`; `probe_start_independence.log:1-3`; `verify_provenance.log:9-15` | - |
| 5.1 | 289 of 648 references not digit-identical; worst 9.98147e-60 | `compare_sweeps.log:5-6` | - |
| 5.1 | the 90-ULP diagnostic | `probe_transfer.log:27` | - |
| 5.2 | B1 relative: direct 26, recursed 7 | derived | `out/second/per_order_ratio.csv:2-67` (A1 rows) |
| 5.2 | first sweep A1 relative direct 28, recursed 5; A1 absolute recursed 20, direct 13 | derived | `out/first/per_order_ratio.csv:2-67` (A1 rows) |
| 5.2 | 14 wrong and 3 tied, and 13 wrong and 1 tied | as Section 3.3 | as Section 3.3 |
| 5.3 | file inventory | this directory; see Layout above | - |

## Deriving the counts that no script prints

Four cited counts are properties of the out CSVs rather than lines in a log. This snippet prints
all of them, and the A1 route counts of Section 5.2:

```python
import csv

def rows(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle))

for sweep in ("first", "second"):
    ratio = rows(f"out/{sweep}/per_order_ratio.csv")
    print(sweep, "winner counts",
          {w: sum(r["measured_winner"] == w for r in ratio)
           for w in ("ship", "product", "tie")})
    pred = rows(f"out/{sweep}/prediction.csv")
    committed = [r for r in pred if r["state"] != "abstain"]
    winner = {(r["band"], r["metric"], r["n"]): r["measured_winner"] for r in ratio}
    print(sweep, "always-direct right at commitments",
          sum(winner[(r["band"], r["metric"], r["order"])] == "product" for r in committed),
          "of", len(committed))
    for metric in ("absolute", "relative"):
        sub = [r for r in ratio if r["band"] == "A1" and r["metric"] == metric]
        print(sweep, "A1", metric,
              {w: sum(r["measured_winner"] == w for r in sub) for w in ("ship", "product")})
```

It prints:

```text
first winner counts {'ship': 42, 'product': 89, 'tie': 1}
first always-direct right at commitments 60 of 69
first A1 absolute {'ship': 20, 'product': 13}
first A1 relative {'ship': 5, 'product': 28}
second winner counts {'ship': 44, 'product': 85, 'tie': 3}
second always-direct right at commitments 59 of 67
second A1 absolute {'ship': 20, 'product': 13}
second A1 relative {'ship': 7, 'product': 26}
```

## Record digests

SHA-256 of every committed record file. The repeated hashes of `scan.csv`, `scan_all.csv`,
`cell_scan.csv` and `per_order_ratio.csv` between the margin-1 and margin-0 directories are
expected: only `prediction.csv` depends on the margin.

| File | SHA-256 |
| --- | --- |
| `check_transfers.log` | `94459e4555a73e89e51155f22edfa1a22ab6561c41827814e0168f619c83103c` |
| `compare_sweeps.log` | `dbd691d85f74df14a0c020778d7fc6959fd7c43d23851dd7a789e63fa7ba93ad` |
| `diagnose_start_independence.log` | `89a2ca4fbc97691128f3f8d8e80f0ed737750a1b9f210d5703260fba03b41b74` |
| `intervals.log` | `16bcd139d0abbe699645988fa2482704bcc153aab5728d8ee7faa6d3ba426832` |
| `out_first.log` | `32ee1e754dcd2a3fddcf0b04e30b8ce987d2af34f478a81c05757c661967551b` |
| `out_first_m0.log` | `21f9fbeb20ccd7951cd402601229c015d9d88392a2fd645e91b0d59a07ac9a63` |
| `out_second.log` | `c083d1224faf76bf3a6cd3fa9e583c13baeed80adfd75c68f532e064f76f208f` |
| `out_second_m0.log` | `5df321f9e2b48984b7fbfebf2482af040ec4b8a4476d828aad4fb58a8bc9d451` |
| `population_winners.log` | `21aaa965238a684ec88f9c40349208d2791e28eb436aa329bf2c985ff483062a` |
| `probe_start_independence.log` | `5bff654ccc7cc16b05fbf96026cccde58299f8e0cf0e6e74d1c43baa2087711a` |
| `probe_transfer.log` | `3899dbd4c73eabf66bb8dadd74b162caa73613e5e851a47fb40db2c2dc2067bb` |
| `same_fit_region.log` | `b2b3e120bdb34b324a98a69f59ee2c57819f25989d889748a6d4c11e926e3131` |
| `verify_provenance.log` | `b0e80921ae01f37dfc3e0553ffed921e323dd9ec7a86030c80645d0ebbcd0357` |
| `out/first/cell_scan.csv` | `50c117e2d367b458d44650fbecc46ea704c094d07c1f051f992bde77a4f3dfce` |
| `out/first/per_order_ratio.csv` | `c8797c7a259b8f8c8bd7b98f347d023a77ae451a75c46bed0c2ceb2079941fa2` |
| `out/first/prediction.csv` | `63ffba89dab012a281c9839dcfffd4906d14a9e0a65e19e91e2b8ee44c4a6c9b` |
| `out/first/scan.csv` | `5436118418232c63a430e339e4b6618db2f63caf6ce260acea366c861f32238a` |
| `out/first/scan_all.csv` | `4179384d0dfcdf0ae5030d74928533e941838a67c4fe27808bdffa68f726930e` |
| `out/first_m0/cell_scan.csv` | `50c117e2d367b458d44650fbecc46ea704c094d07c1f051f992bde77a4f3dfce` |
| `out/first_m0/per_order_ratio.csv` | `c8797c7a259b8f8c8bd7b98f347d023a77ae451a75c46bed0c2ceb2079941fa2` |
| `out/first_m0/prediction.csv` | `8a020fddc562ba0f7865780fc72f65a063a9fd81c0c49a579dd808667f926be2` |
| `out/first_m0/scan.csv` | `5436118418232c63a430e339e4b6618db2f63caf6ce260acea366c861f32238a` |
| `out/first_m0/scan_all.csv` | `4179384d0dfcdf0ae5030d74928533e941838a67c4fe27808bdffa68f726930e` |
| `out/second/cell_scan.csv` | `e2d6fb80d5d12666cb80b60a92e7a4df2777c7582d4202c99d21e62e4c915b2f` |
| `out/second/per_order_ratio.csv` | `4a1fdd094845be446a72e25eebeb49821dde048a8102e6aa303f749c2cbe8953` |
| `out/second/prediction.csv` | `097fbe2943ba30cfb48cbfd827ccfbd424c78b7a6cb603578a1ece69640aafba` |
| `out/second/scan.csv` | `13b4fdf775e709beef8c0d7a278a10c16bb8d03a5fb8df2cfdf7ae1f0b56656d` |
| `out/second/scan_all.csv` | `dcbdd6d1df1651ed65d76c6f8ac71f0905ec6dd605bd17e4bfcd6172397fdfa3` |
| `out/second_m0/cell_scan.csv` | `e2d6fb80d5d12666cb80b60a92e7a4df2777c7582d4202c99d21e62e4c915b2f` |
| `out/second_m0/per_order_ratio.csv` | `4a1fdd094845be446a72e25eebeb49821dde048a8102e6aa303f749c2cbe8953` |
| `out/second_m0/prediction.csv` | `f51a0d50e9a2f1062ec22477b3848520bac4f233720468f3a6af140810cd229d` |
| `out/second_m0/scan.csv` | `13b4fdf775e709beef8c0d7a278a10c16bb8d03a5fb8df2cfdf7ae1f0b56656d` |
| `out/second_m0/scan_all.csv` | `dcbdd6d1df1651ed65d76c6f8ac71f0905ec6dd605bd17e4bfcd6172397fdfa3` |

## Caveats

- **`verify_provenance.py` cannot reproduce its committed log from this deposit alone.** The
  original log also compared the paper-tree copy of `analyze.py` with the deposited copy and
  printed both hashes (`verify_provenance.log:4-8`). A clone has no paper tree. The deposited
  script is a deposit-only variant: it checks the three inputs, the deposited `analyze.py` hash,
  and the two provenance records, and prints `PASS` lines with `inputs/` labels. Its script output
  is thirteen lines (fourteen with the appended exit line), not eighteen. No hash or provenance
  value changes; the committed log remains the record of the original run.
- **`intervals.py` carries a short import shim.** `surrogate-control/surrogate_control.py`
  names the paper tree it was first run from in a module-level `TREE` and imports `analyze`. The
  shim registers this deposit's `amplification-form/analyze.py` under the name `analyze`, points
  `S.TREE` at `amplification-form/`, and aborts if the control imported any other module. The
  deposit's `analyze.py`, `cells.csv`, `reference.csv` and `provenance.json` are byte-identical
  to the paper tree's, so the shim changes no computed value. It is the only script whose control
  flow was changed, and only in the import block.
- **`as/raw/` in `STRENGTHEN.md` is `inputs/` here.** The report was written before this deposit
  existed and uses the paper-tree path names; the numbers are the same files.
- **`prepare_inputs.py` is a deposit convenience, not part of the original run.** There is no
  `prepare_inputs.log` in the record. Its output is the three decompressed size/hash lines, and
  its SHA-256 checks are the authority that `inputs/` is the deposit's own data.
- **The compression layout is not uniform.** Only `product-amplification/results/cells.csv.xz`,
  `product-amplification/data/reference.csv.xz` and `amplification-form/results/paths.csv.xz`
  exist as streams; the `amplification-form` cells and reference are plain deposited files. The
  first sweep needs the two `product-amplification` streams; the second sweep is read in place.
- **Five logs have no exit line, and one report citation is wrong about it.** `STRENGTHEN.md`
  Section 0.4 lists `population_winners.log:62` as the exit evidence, but line 62 is the final
  tie line; the log has no `EXIT=` line. The exit code 0 for that run was observed at the runner
  and is not recoverable from the log. The other four no-exit logs are marked "0 observed" in the
  report; `population_winners` is the one row that is not. All five are listed above.
- **`out_second_m0.log:7` is expected to report 33 mismatches.** It compares margin-0 predictions
  against the deposited margin-1 `prediction.csv`. A fresh run should print the same 33.
- **Byte identity of logs is a Windows property.** On a non-Windows host the script output is the
  same text with LF endings, and the appended exit lines have the same ending as the rest; the
  CSVs remain CRLF. Compare logs textually there, not by SHA-256.
- **`STRENGTHEN.md` is not in this repository.** It is in the Paper A tree at
  `papers/boys-crossover/STRENGTHEN.md`. This deposit reproduces its evidence, not the report.
