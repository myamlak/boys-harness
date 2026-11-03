# Certified per-call rounding bounds for the two region-A Boys entries

This directory is the record for the certified-bound section of the crossover
paper. It bounds both region-A entries of `boys` at tag `v5.2.4` operation by
operation, validates the bound on every measured cell, and records what it can
certify.

The direct entry is `BoysRegionAProduct<ProductMode::kFp64>`; the recursed
entry is `BoysAllOrders`. The bound models the map, the Chebyshev basis
recurrence, the products and pairwise sums, the exact division at zero, the
order-32 downward seed below `theta0 = 1.0855252345349333`, the upward
recurrence above it, and the shared `0.5*exp(-x)`. Each operation carries a
running-error radius; truncation is the exact-real route value at 70 digits
minus the 60-digit reference.

## What is here

- `code/`: the model, exact route, certification, tier variants, replay and
  oracle checks, and the provenance writer.
- `explore/`: the paper's 9.6 and pre-exponential 93 diagnostics, the cell-set
  comparison, and the structural checks.
- `inputs/`, `results/`: the 55,176 cell records, band and model summaries,
  logs and `commands.md`.
- `exp/`: the libm exponential probe source and the measured and reproduced
  bit files.
- `provenance.json`: the kernel, sweep, manuscript, compiler flags and file
  digests of the run.

## How to reproduce

Run from this directory with Python 3.14 and mpmath. The exact commands,
arguments and exit codes are in `results/commands.md`; the principal runs are

    python -I code/certify.py
    python -I code/replay.py <paper>/sweeps/amplification_form --both-downward
    python -I code/oracle_check.py
    python -I explore/dps_check.py
    python -I explore/paper93.py
    python -I code/tiers.py
    python -I code/make_provenance.py

The replay path names the paper's private tree; the same measurement data is
deposited at `../amplification-form`. `inputs/x_bits.txt` carries the argument
bits used by the exponential probe.

## What is not here

`exp/exp_probe.exe` and `exp/exp_probe.obj` are MSVC build products and are not
deposited; rebuild `exp_probe.cpp` with the flags in `results/commands.md`. Its
re-execution reproduces `exp/exp_bits_reproduced.txt` byte for byte, and the
comparison with `exp/exp_bits.txt` is in `exp/exp_check.log`.

## Headline

The bound is sound at all 55,176 entry-cells (0 violations, worst
measured-to-certified ratio 0.719) and certifies 0 of 132 comparisons and 0 of
65 commitments. The obstruction is the losing entry's rounding charge, a median
83 times its truncation; 16 of the 65 comparisons at the paper's per-order
cells are truncation-opposite and cannot be certified by any charge reduction.
