# Clean-edge 60-digit control

This directory records the full forty-row clean-edge coefficient-count ladder
re-solved at 60-digit evaluation with the BRASIL independent rational solver.
The ladder is the one deposited one level up in `../clean-edge-2026-11-22/`:
the clean left edge `a = 62`, five widths (`2^-7`, `2^0`, `2^5`, `ship` and
`paper`) and eight relative targets from `1e-8` to `1e-15`, forty cells. The
intervals, targets and target semantics are the record's; this control changes
how a fit is evaluated, not the ladder.

## Settings

- 2,001-node certifying grid (2,000 intervals); an accepted fit is also scored
  on a 4,001-node check grid.
- Rational denominator degree cap `k <= 12`, denominator guard `1e-3`.
- The fit is produced by BRASIL (`baryrat` 2.1.2) in binary64, with the
  default settings and the tightened settings only where the default route
  finds no meeting fit.
- A fit meets at 60 digits when its monomial coefficients, held at 60
  significant digits and evaluated at 60 decimal digits, are within the
  relative target of `F_0` at every node. `F_0` is evaluated at 70 decimal
  digits; the coefficient conversion also uses 70 decimal digits.
- 8 processes with a 3,600 s per-cell timeout; the completing run took
  1,471 s in total.

The evaluation is at 60 digits; the fit itself remains binary64. Holding the
binary64 coefficients at 60 significant digits measures evaluation arithmetic,
not a 60-digit refit.

## Coverage

All 40 cells completed: 40 `ok`, no errors, no timeouts and no cell without a
60-digit count. Twenty-seven cells reproduce the record's count and thirteen
record a lower count; none required an upward count. BRASIL made 690 solves,
398 of which reported `converged=False`, with no exceptions.

## Counts

Each cell reads `record / 60-digit`.

| Width | `1e-8` | `1e-9` | `1e-10` | `1e-11` | `1e-12` | `1e-13` | `1e-14` | `1e-15` |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `2^-7` | 3 / 3 | 3 / 3 | 3 / 3 | 3 / 3 | 3 / 3 | 4 / 4 | 4 / 4 | 4 / 4 |
| `2^0` | 4 / 4 | 5 / 5 | 5 / 5 | 6 / 6 | 6 / 6 | 6 / 6 | 7 / 7 | 7 / 7 |
| `2^5` | 9 / 9 | 10 / 9 | 11 / 10 | 12 / 11 | 13 / 12 | 14 / 13 | 15 / 14 | 22 / 15 |
| `ship` | 8 / 8 | 9 / 8 | 9 / 9 | 10 / 10 | 11 / 11 | 12 / 12 | 13 / 13 | 16 / 13 |
| `paper` | 8 / 7 | 8 / 8 | 9 / 9 | 10 / 10 | 11 / 11 | 12 / 11 | 13 / 12 | 14 / 13 |

The counts are unchanged at every target for `2^-7` and `2^0`. At `2^5` the
count is one lower from `1e-9` onward and seven lower at `1e-15` (22 to 15).
At `ship` it is one lower at `1e-9` and three lower at `1e-15` (16 to 13). At
`paper` four counts are one lower: `1e-8`, `1e-13`, `1e-14` and `1e-15`.

## Verdict

The primary fit uses all eight targets; the secondary fit uses the contract
range `1e-8` to `1e-13`. The rule is unchanged from the record: the rational
fit fails at a width when the free-line residual exceeds 2 or the absolute
second difference is at least 2, and Chebyshev holds when its residual is at
most one stored coefficient.

The record's primary verdict is **mixed**: the rational fit fails at `2^5` and
`ship`. At 60 digits the rational fit holds at all five widths, so the primary
verdict is **transient-specific**. The secondary verdict is
transient-specific in both the record and the control.

Primary free-line fits:

| Width | Record residual | Record d2 | Record | 60-digit residual | 60-digit d2 | 60-digit |
| --- | --- | --- | --- | --- | --- | --- |
| `2^-7` | 0.464 | 1.0 | holds | 0.464 | 1.0 | holds |
| `2^0` | 0.452 | 1.0 | holds | 0.452 | 1.0 | holds |
| `2^5` | 3.500 | 6.0 | fails | 0.583 | 1.0 | holds |
| `ship` | 1.417 | 2.0 | fails | 0.417 | 1.0 | holds |
| `paper` | 0.583 | 1.0 | holds | 0.464 | 1.0 | holds |

The deep-target departures that produced the record's two primary failures,
22 at `2^5`/`1e-15` and 16 at `ship`/`1e-15`, are not reproduced by this
control. Under it the largest primary rational residual is 0.583 at `2^5`,
within one stored coefficient. In the secondary fit the 60-digit rational
residuals are at most 0.476 with second differences at most one.

## Caveats

- **Joint independent-solver and evaluation-precision control, not a pure
  rescore.** Every 60-digit count comes from re-solving the cell with BRASIL.
  The record's fit was made by the kernel's Remez or Lawson solver, so a
  difference from the record can come from the solver as well as from
  evaluation precision; it is not attributable to precision alone.
- **The recorded fits were not re-scored.** The recorded kernel fits store no
  coefficients, so they cannot be evaluated at 60 digits. Each cell record
  carries the record's count and pair, and the pair is re-solved with BRASIL
  for comparison; the record's own fit is not otherwise used.
- **Counts are existence bounds.** BRASIL reports `converged=False` on 398 of
  its 690 solves. A count is the lowest stored count at which a swept
  admissible pair met the target, so it is an existence bound under this
  solver and these routes, not a proof of minimality. A lower count is
  excluded only for the pairs and settings swept.
- **The Chebyshev arm remains binary64.** Only the rational arm is re-solved.
  The Chebyshev stored counts, and the statement that Chebyshev holds at every
  width, are carried from the binary64 record and are not controlled at 60
  digits.

## Files

| File | Contents |
| --- | --- |
| `combined_60dig.json` | the 40-row combined record: per-cell counts, the per-width free-line table, the summary and the input digests |
| `cells/` | the 40 per-cell records, with the BRASIL levels, the accepted fit and its score on the 4,001-node check grid |
| `analysis-full/analyze_edge62.py` | the pre-registered clean-edge analysis, standard library only |
| `analysis-full/edge62_sweep.json` | the 40-row analysis input: the 60-digit rational counts substituted into the clean-edge record, with the binary64 Chebyshev counts unchanged |
| `analysis-full/edge62_analysis.json` | the analysis output |
| `run_cell_60dig.py` | the per-cell runner |
| `run_all_60dig.py` | the launcher that runs the 40 cells and writes the combined record |

## Identity

| Item | SHA-256 |
| --- | --- |
| `combined_60dig.json` | `dd5297ebbc4d4813cc15e898a6bdbf0110a6593d5c5f3fad42c0b5de964eeb0e` |
| `run_cell_60dig.py` | `eb165d1a1bfe1fbbab051ec5ec267cdedd2cad47dc41e786701fde5694d6a9a8` |
| `run_all_60dig.py` | `d2762e23efabc3ad6f92cb6dd341624d2418066c1ea21875afa0b6678a06fdea` |
| `analysis-full/analyze_edge62.py` | `8cf9d4159608168947fdf99af59bf941e7f32e10f6fcdc67febe88ccd65f9f87` |
| `analysis-full/edge62_sweep.json` | `f61b307d7562330a1c20d185d8f5b8b8ad0ba6674725c37323ec207fb7897f0a` |
| `analysis-full/edge62_analysis.json` | `8cc4f6e98dc1f67f400c5d09d34f010bcab469d55ecd46a502dbd43870f82279` |
| `../clean-edge-2026-11-22/edge62_sweep.json` | `963facd9e18c96b547b6543c8713d5dd8294104d2776fb8d71dd30c0d40bff8b` |
| `../control-60dig/control_60dig_brasil.py` | `8ffedf42c97111fd7512c48d85cedf282007b0aee664fc171569aab7c3e0f559` |

The last two rows are the binary64 clean-edge record and the 60-digit control
core, both deposited in this repository. Their digests are also carried in
`combined_60dig.json`.

## Reproducing and verifying

The analysis regenerates its output from the deposited analysis input:

```sh
cd analysis-full
python -I analyze_edge62.py
```

A fresh run rewrites `edge62_analysis.json` byte-identically. The script
imports only the standard library.

The full run was launched with:

```sh
python run_all_60dig.py --jobs 8 --timeout 3600
```

It needs Python 3.10 or newer, `mpmath`, `numpy`, `scipy` and BRASIL
(`baryrat` 2.1.2). It reads the binary64 clean-edge record and the 60-digit
control core, both deposited in this repository as named in the identity
table. The deposited scripts carry the measurement machine's absolute paths
for these two inputs; point `REC` and `CTRL` in `run_all_60dig.py`, and `CTRL`
in `run_cell_60dig.py`, at the deposit before running. The launcher writes
`cells/`, `logs/` and `combined_60dig.json` in this directory.

A fresh combined record differs from the deposited one in its generation time
and wall times; the deposited record's SHA-256 is the one in the identity
table.
