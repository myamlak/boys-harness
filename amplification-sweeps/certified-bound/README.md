# Certified per-call rounding bounds for the two region-A Boys entries

This directory is the record for the certified-bound section of the crossover
paper. It bounds both region-A entries of `boys` operation by operation,
validates the bound on every measured cell, and records what it can certify.
The current kernel reference is tag `v5.3.0`. The record was measured at
`v5.2.5`; the stored tables the model parses are byte-identical at `v5.2.6`,
so the certificate applies at that tag as well, and `v526-identity/` carries
the identity check and the replay.

The direct entry is `BoysRegionAProduct<ProductMode::kFp64>`; the recursed
entry is `BoysAllOrders`. The bound models the map, the Chebyshev basis
recurrence, the products and pairwise sums, the exact division at zero, the
order-32 downward seed below `theta0 = 1.0855252345349333`, the upward
recurrence above it, and the shared `0.5*exp(-x)`. Each operation carries a
running-error radius; truncation is the exact-real route value at 70 digits
minus the 60-digit reference.

## The repaired record (2026-11-06)

The first deposit had two defects, both repaired on branch `fix/cert-underflow`
and rerun on the full grid into `results-underflow-fix/`. The old 2026-11-03
record in `results/` is preserved byte for byte.

**1. The underflow radius.** The grid includes `x = 2^-1074` (the least
positive subnormal). For `0 < x < theta0` the downward recurrence seeds at
order 32; at this argument the exact product `x * F_32` is
`7.6010099360191773e-326`, which rounds to zero. The old multiplication rule
charged `u_+|z|` with `z = 0` and the propagated term `x*rho_b` also rounded
to zero, so the model charged a zero radius for a multiplication whose exact
result is nonzero. No final interval was violated, but the operation-level
enclosure was not proved at that node.

The rule now carried by `mul`, `fma`, `mapping`, `add` and `sub` is:

    bracket = u_+|z| + (propagated terms)
    if |z| < 2^-1022 or bracket < 2^-1022:
        bracket += 2^-1071
    radius = bracket * (1 + 2^-40)

The absolute term is derived as follows. IEEE-754 round-to-nearest bounds an
operation's error by `max(u|ab|, 2^-1075)`, where `2^-1075` is half the least
positive subnormal and the invariant absolute floor below the normal range.
`2^-1075` is not a binary64 value: at ties-to-even it rounds to zero. The
rule therefore carries the representable `2^-1071 = 16*2^-1075`: one half-ulp
for the operation's own rounding and fifteen for the absolute roundings of
the radius arithmetic when its terms are subnormal (the longest rule performs
at most eleven roundings: four products, four bracket sums, the allowance
addition, the operation and the `SAFETY` multiply; below the normal range the
allowance addition is exact, so at most ten, and sixteen half-ulps leave a
margin of at least five). The term
is added only when the computed result or the bracket is below the normal
range `2^-1022`; above it, the existing `u_+|z|` term and `SAFETY` cover the
same floor, because `SAFETY` inflates a normal-range bracket by at least
`2^-1022 * 2^-40 = 2^-1062`, thirteen binary orders above `2^-1075`.

For `mul` and `fma` the allowance covers the product's or the fused sum's own
rounding. For the mapping's division it covers the quotient's. For `add` and
`sub` the operation itself is exact in the subnormal range (the exact sum of
two binary64 values is a multiple of `2^-1074`, hence representable when it
lies below `2^-1022`); there the allowance covers the radius arithmetic only.
`scale2` is exact for every finite value without overflow.
`refined_div`'s derivation assumes normal operands; the model asserts that
condition and refuses a call outside it, and every call on the grid satisfies
it (the downward denominator is `l + 1/2` in `[1/2, 63/2]`, the upward one is
`x >= theta0`).

The allowance propagates through the chain. At `x = 2^-1074` the affected
product's radius is `4e-323` against an exact product of `7.601e-326`; the
order-31 step that consumes it has a total radius of `7.05e-18`, so no
recorded `C` changes. The
regression test `code/test_underflow.py` records the exact and computed
intermediate values and the radii at every affected step, and
`results-underflow-fix/underflow_test.csv` is its record.

**2. Portability.** Every path is now resolved by `code/record_paths.py` from
the command line, the documented environment variable, or a default relative
to the record:

| input | option | environment | default |
| --- | --- | --- | --- |
| kernel checkout | `--kernel` | `BOYS_ROOT` | `<repo>/../boys` |
| sweep data | `--sweep` | `CERT_SWEEP` | `../amplification-form` |
| model record | `--model-record` | `CERT_MODEL_RECORD` | `inputs/model-record` |
| exponential bits | `--exp-bits` | `CERT_EXP_BITS` | `exp/exp_bits.txt` |
| deposit | `--deposit` | `CERT_DEPOSIT` | `../apriori-bound` |
| run directory | `--out` | - | required |

No script inserts an external working copy into `sys.path`; scripts insert
their own `code/` or `explore/` directory so that `python -I` can import the
record's modules, and the in-repo deposit directory only while importing its
sibling-importing modules. Every output is written beneath the `--out`
directory.

The model record is `inputs/model-record/{prediction,prediction_record}.csv`,
the paper model's 132 comparisons and its 65 commitments (sha256 `305b2a28...`
and `d146ab20...`). The deposited sweep's own `prediction.csv` is a different
model revision with 67 commitments, so the certificate reads the paper's
record deposited here.

The kernel is identified by a canonical digest of exactly the parsed tables
(`kPieces`, `kPieceStart`, `kTierThresholds`, `kExtendedBcoeffs`,
`kExtendedBX0`, `kX0`) rather than by the whole-file digest of
`boys_coefficients.hpp`. The parsed content is identical in the v5.1.1 sweep
tree and the v5.2.5 and v5.2.6 releases; the v5.2.5 header differs from the
v5.1.1 one only in the rational region-B seed, which the named policy does
not read, and the v5.2.6 header is byte-identical to the v5.2.5 header.
`provenance.json` in the run directory records both whole-file digests and
the parsed digest.

## What is here

- `code/`: the model, exact route, certification, tier variants, replay,
  oracle and underflow checks, the path resolver and the provenance writer.
- `explore/`: the paper's 9.6 and pre-exponential 93 diagnostics, the cell-set
  comparison, and the structural checks.
- `results/`: the 55,176 cell records, band and model summaries, logs and
  `commands.md` of the 2026-11-03 run.
- `inputs/`: the probe argument bits and `model-record/`, the paper model's
  prediction and comparison records.
- `results-underflow-fix/`: the repaired rerun of 2026-11-06: `summary.json`,
  the per-cell and band CSVs, the logs, the underflow record and `commands.md`.
- `exp/`: the libm exponential probe source and the measured and reproduced
  bit files.
- `v526-identity/`: the coefficient-header identity between v5.2.5 and
  v5.2.6, with the replay that reproduces the recorded values from the
  v5.2.6 headers.
- `provenance.json`: the kernel, sweep, compiler flags and file digests of
  the 2026-11-03 run; the repaired run's provenance is
  `results-underflow-fix/provenance.json`.

## How to reproduce

Run from this directory with Python 3.14 and mpmath. The external input is a
boys checkout; the commands select the current release, tag `v5.3.0`. The
record's identity conclusion is the `v5.2.5`/`v5.2.6` coefficient-header
comparison in `v526-identity/`. The default kernel is `<repository>/../boys`,
i.e. a clone beside the harness checkout (else pass `--kernel` or set
`BOYS_ROOT`). From the repository root:

    git clone --branch v5.3.0 https://github.com/myamlak/boys.git ../boys

The full record, with every output beneath `results-underflow-fix/`:

    python -I code/certify.py --out results-underflow-fix
    python -I explore/harnesscells.py --out results-underflow-fix
    python -I code/replay.py --both-downward > results-underflow-fix/replay.log 2>&1
    python -I code/oracle_check.py --out results-underflow-fix
    python -I explore/dps_check.py > results-underflow-fix/dps_check.log 2>&1
    python -I explore/paper93.py > results-underflow-fix/paper93.log 2>&1
    python -I code/tiers.py --out results-underflow-fix
    python -I code/test_underflow.py --out results-underflow-fix
    python -I exp/exp_check.py > results-underflow-fix/exp_check.log 2>&1
    python -I code/make_provenance.py --out results-underflow-fix

Use a fresh directory for a new run; do not overwrite `results/` or
`results-underflow-fix/`. The exact commands, arguments and exit codes of the
2026-11-06 run are in `results-underflow-fix/commands.md`. That run did not
rerun `tiers.py`: its only model-dependent input is the per-order radius `C`,
whose full table is byte-identical, so the 2026-11-03 seven-variant verdict
stands. The command is listed because it remains part of the record's
reproduction.

## Clean-checkout smoke test

From a fresh harness checkout with the kernel clone above, from this
directory (the default kernel path applies; add `--kernel PATH` if the clone
is elsewhere):

    python -I code/test_underflow.py
    python -I code/replay.py --limit 3
    python -m py_compile code/*.py explore/*.py exp/*.py

The first prints `PASS: the underflow allowance encloses the lost product at
x = 2**-1074`; the second prints
`total {'ship': 198, 'product': 198} mismatches {}`. Both use only the
checkout and the kernel clone: the sweep is `../amplification-form`, the model
record is `inputs/model-record`, and the exponential bits are
`exp/exp_bits.txt`.

## What is not here

`exp/exp_probe.exe` and `exp/exp_probe.obj` are MSVC build products and are not
deposited; rebuild `exp_probe.cpp` with the flags in
`results-underflow-fix/commands.md`. Its re-execution reproduces
`exp/exp_bits_reproduced.txt` byte for byte, and the comparison with
`exp/exp_bits.txt` is in `exp/exp_check.log`.

## Headline

The repaired bound is sound at all 55,176 entry-cells (0 violations, worst
measured-to-certified ratio 0.71875 at (A1, 0, 24, ship)) and certifies 0 of
132 comparisons and 0 of 65 commitments, exactly as the 2026-11-03 record did:
the underflow allowance changes no recorded radius. The obstruction is the
losing entry's rounding charge, a median 83 times its truncation; 16 of the
65 comparisons at the paper's per-order cells are truncation-opposite and
cannot be certified by any charge reduction. The operation-level enclosure
that the old multiplication rule failed at `x = 2^-1074` is now recorded by
the regression test.
