# Corrected rational counts

These are the records the paper's rational counts now come from. They were
produced by an independent solver, not by the kernel's own Remez and Lawson
implementations.

The two original files are copied verbatim from the verification run of
2026-10-31. The 60-digit audit and the derived reported counts were added for
the Q-429 remedy of 2026-11-15:

| file | sha256 |
| --- | --- |
| `corrected-counts.csv` | `e11899412881e23c65756c7949a2e2f9564c8e30ccf33c397b6a5549ba2bc335` |
| `corrected-counts.json` | `9bb8ac0633b10e8a7922684ca8a419135103f2c7825b4d06b9c71603399b7195` |
| `corrected-counts-60dig.json` | `18a41a42798f0cd1fead4a143309ded3112d1bba94727b093bf27df02e3256e7` |
| `reported-counts.json` | `bb54dfbb86f636406721e6074c21207fe88ce63f4ae07adcf49d1054706583d1` |

## What produced them

The solver is `baryrat` 2.1.2, the BRASIL rational minimax algorithm. It is
independent of the kernel's Remez exchange and Lawson iteratively reweighted
least-squares solvers. The target function is `F_0` on `[x0, x0 + W]` with
`x0 = 11.899848152108484`; the rational form is `p(t) / (1 + q1 t + ... + qk
t^k)` in `t = 2(x - x0)/W - 1`, and the stored count is `(m + 1) + k` with
`k <= 12`.

## The acceptance rule, as applied

The paper's own counting rule was used unchanged. A fit is accepted when its
coefficients, converted exactly to monomials and rounded to binary64, evaluated
by Horner in binary64 by both the fused and the separate routes, hold the
worst error over the 2,001-node grid in the target's kind: relative
`1e-13`, absolute `5e-14`, or relative `1e-14`. The error is the worse of the
two routes. The denominator is guarded at `min |den| >= 1e-3 max |den|`. The
2,001-node grid certifies; a 4,001-node grid checks. The reference is `F_0` at
70 decimal digits.

A count is therefore the smallest stored count at which an exhibited
admissible fit meets, under this solver. It is an upper bound, not a proof of
minimality. One ladder cell is higher than the paper's: `2^5` at `10^-15` is
25 `[22/2]` here against the paper's 24 `[21/2]`, whose coefficients are not in
the record and so could not be re-evaluated. That cell is not independently
reproduced.

## The 60-digit audit and the reported count

The audit re-solved every one of the 45 Group A cells with the v1.0.39
60-digit control and scored two things separately for each recorded pair: its
nominal truncation `(m + 1) + k`, and the full stored witness whose coefficient
arrays the solver zero-pads. The paper reports one count per cell, defined
once and applied by `tools/derive_reported_counts.py`:

* the nominal count of the recorded pair, when that truncation meets the cell
  target on both the 2,001-node certifying and 4,001-node checking grids under
  both the 60-digit and binary64 evaluations (40 cells);
* otherwise the nominal count of a genuine re-solve at the control's count
  `count60` whose truncation meets the target on both grids under binary64 or
  60-digit evaluation (1 cell, `2^5` at the absolute target, pair `[7/6]` under
  binary64 or `[6/7]` under 60-digit evaluation);
* otherwise the full stored witness of the recorded pair, because no genuine
  fit exists at that count (4 cells: `2^5` at `10^-13` relative, and `2^4`,
  `2^5` and the paper width at `10^-14` relative).

The derived record `reported-counts.json` carries the rule, every cell's
nominal, full, reported and Chebyshev counts, and the mode that was applied.
The ladder and the kernel's own route are not re-solved by the audit; their
counts remain nominal degree counts of unaudited records, as the paper states.

## Deposit

The two records above and the raw cell files are deposited under the data tag
`count-audit-2026-11-16`, whose tree is v1.0.39 plus this audit and no other
change. The tag is a data deposit, not a release; the harness VERSION in the
tagged tree remains 1.0.39.

| file | sha256 |
| --- | --- |
| `corrected-counts-60dig.json` | `18a41a42798f0cd1fead4a143309ded3112d1bba94727b093bf27df02e3256e7` |
| `reported-counts.json` | `bb54dfbb86f636406721e6074c21207fe88ce63f4ae07adcf49d1054706583d1` |
| `raw-60dig/sweep_cells.py` | `412a234ca39b1d73498f28700830e9fdf269f5aa6d42ab1ae2f8e3a02eb728ae` |
| `raw-60dig/assemble_record.py` | `1cf2a21a2f3668c61fea16244f0c174673bea1034a872b8cfd137a475656d284` |
| `../tools/derive_reported_counts.py` | `4f93c43976309e271084449c54887b1ed6d61f87c7ff7e9627317a1e020d0e52` |

`raw-60dig/` carries the 45 per-cell records and the four supplementary level
sweeps (`cell_g10`, `cell_g26`, `cell_g40`, `cell_g41`) that the two drivers
wrote, with the two drivers beside them. The drivers expect the campaign tree
beside them (`coefficient-count-sweep/control-60dig/` and the corrected
records); run `sweep_cells.py --index N` for one cell and
`assemble_record.py` to rebuild the audit. `derive_reported_counts.py` reads
the audit and writes the derived record; it verifies the audit's digest and
stops if an aggregate changes. The solver is baryrat 2.1.2, the control is
`coefficient-count-sweep/control-60dig/control_60dig_brasil.py` at v1.0.39,
and the absolute-target wrapper is the only addition to the tagged control.


