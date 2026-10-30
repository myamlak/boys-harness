# Held-out test

The decision rule and its abstention margin were introduced and scored on the
first two sweeps, and a rule scored on the data that suggested it proves
little.  This sweep freezes the rule and the margin exactly as
`../../tables/prediction-record.tex` records them, and applies them, with no
adjustment, to a third argument set that is disjoint from both earlier
sweeps.

## The frozen rule

    B = band maximum of the recursed bound
    D = band maximum of the direct error
    T = band maximum of the accumulated-arithmetic term

    abstain   if |B - D| <= T
    direct    if D < B - T
    recursed  if D > B + T

The margin is one unit of `T`.  It is hard-coded in `analyze.py` as
`FROZEN_MARGIN = 1.0`; the script offers no option that could change it and no
margin ladder, so nothing in this sweep can tune it on the new data.

The model that forms B, D and T is a line-for-line copy of
`../amplification_form/analyze.py`: the same amplified seed term, the same
accumulated-arithmetic term, the same population (arguments at or above
`kTierThresholds[0]`), the same unserved-order rule (both entries read the same
stored fit there, so the model abstains by construction), and the same band
maxima.

The run proves the copy before it scores anything.  `analyze.py` first applies
it to the second sweep's own cells and reference and requires it to reproduce
`../amplification_form/results/prediction_record.csv` field for field, all 132
rows at 17 digits, and then `../../tables/prediction-record.tex` through the
table's own two format steps (six significant figures for a measured value,
the shortest spelling of its argument) and in the same words.  That
check also reproduces the original score: 65 commitments, all 65 right, 67
abstentions.  A mismatch stops the run before any held-out number is written.

## The held-out argument set

`grid.py`, rule H.  Five rules per band, with the band edges read from the
pinned tree (`kRegionA1Edge`, `kRegionAEnd`):

    H1  shifted uniform, 300 points per band:
            f = (3k + 1) / 900,                        k = 0 .. 299
    H2  Chebyshev-Gauss, 136 points per band:
            f = (1 + cos(pi (j + 1/2) / 136)) / 2,     j = 0 .. 135
    H3  base-2 edge offsets, 27 per edge, both edges:
            f = 2^-m and f = 1 - 2^-m,                 m = 12 .. 38
    H4  upper-half window, 96 points per band, quarter-step offset:
            f = 1/2 + (3j + 2) / 576,                  j = 0 .. 95
    H5  lower-quarter window, 64 points per band:
            f = (3j + 2) / 768,                        j = 0 .. 63

    x = left + (right - left) f

650 points per band, 1300 arguments.  None of the rules uses the committed gate
grid, so the held-out set is independent of the arguments both earlier sweeps
share.  The density, the offsets (one-third and two-thirds of a step) and the edge
coverage (base-2, not base-10) all differ from `product_amplification/grid.py`
and `amplification_form/grid.py`.  `grid.py` asserts exact-double disjointness
against both earlier grids and fails the run on a single collision; this run
checked 9541 distinct earlier arguments and found none.  The rule is stated in
full above and in the script's docstring, so it is regenerable.

## Commands

Run from a Visual Studio developer shell, without `python -I` (mpmath lives in
the user site).  The pipeline does not run git.

    powershell -ExecutionPolicy Bypass -File run.ps1

The steps it runs, in order:

    python grid.py --boys-tree C:\Users\myaml\dev\_boys-v511 --out data\grid.csv
    python provenance.py --boys-tree C:\Users\myaml\dev\_boys-v511 \
        --harness sweep_heldout --harness-dir . --flags "<compile line>" \
        --out results\provenance.json
    python ..\common\reference.py --boys-tree C:\Users\myaml\dev\_boys-v511 \
        --grid data\grid.csv --out data\reference.csv \
        --validation data\reference_validation.txt --dps 70 --digits 60
    cl /nologo /std:c++latest /O2 /EHsc /fp:precise /arch:AVX2 \
        /I C:\Users\myaml\dev\_boys-v511\include \
        /Fo build\sweep_heldout.obj /Fe build\sweep_heldout.exe sweep_heldout.cpp
    build\sweep_heldout.exe --grid data\grid.csv --cells results\cells.csv \
        --admission results\admission.csv
    python analyze.py --tree . --provenance results\provenance.json

## Results

Kernel provenance: v5.1.1 at 834bdbb18a30ac8538b519436429f78d41f187ce, proved
by content with `git_used = false`.  All 15 kernel files the second sweep's
record names match digest for digest; 62 kernel files and the reference chain
are digested, and none differs.

Reference: 1300 arguments, 33 orders, 70 working digits, 60 printed.  Every
evaluated route is inside its stated bound: R2 series 4.21e-62, R2b imported
series 8.75e-62, R3 quadrature 2.14e-71, R5 identity 9.47e-67, R6 F0 against
erf 3.62e-71.  R4 and R7 are not evaluated because the held-out grid excludes
the committed gate arguments by design.

Admission: the shipped ladder is bit-identical to the replica on all 1300
arguments and all 24,332 checked cells (upward 1133 arguments and 18,821
cells, downward 167 arguments and 5,511 cells).  The recursed values scored
here are the kernel's own recurrence.

Model population: arguments at or above `kTierThresholds[0]` =
1.0855252345349333: A1 483 of 650, A2 650 of 650.

Held-out score, orders 0-32, at the frozen margin of one unit of T:

| band | metric | commits | right | wrong | abstains |
| --- | --- | --- | --- | --- | --- |
| A1 | absolute | 13 | 13 | 0 | 20 |
| A1 | relative | 8 | 8 | 0 | 25 |
| A2 | absolute | 20 | 20 | 0 | 13 |
| A2 | relative | 17 | 17 | 0 | 16 |
| all | | 58 | 58 | 0 | 74 |

The rule commits 58 times and is right in all 58.  It abstains 74 times, a
rate of 56.1 per cent on this grid.  Measured on the held-out cells, the
direct entry wins 88 of the 132 comparisons, the recursed entry 43, and one is
a tie (A2 absolute, order 18), which the rule abstains on.

The commitments are 53 to direct and 5 to recursed.  The recursed commitments
are A1 relative order 0 and A2 relative orders 5, 6, 7 and 8.  Every other
commitment is direct, which is the same dominant direction the earlier sweeps report.

One asymmetry is inherited from the original analysis and worth naming: the
model's B and D are maxima over the model population, while the measured
winner is decided by band maxima over the whole held-out grid, including the
A1 arguments below `kTierThresholds[0]`.  At A1 absolute order 0 the measured
direct worst is attained at x = 7.09e-07, outside the population; the
commitment is still right.  The second sweep has the same structure.

## What this establishes and what it does not

Establishes, on the pinned v5.1.1 kernel, its own 60-digit reference, and an
argument set disjoint from both earlier sweeps:

* the implementation this sweep scores is the frozen rule: the model
  reproduces the 132-row table field for field and at its printed
  precision, and reproduces the original 65 / 65 / 67 score, before it sees
  a held-out cell;
* applied unchanged at the frozen margin, the rule is right on all 58
  commitments it makes on the held-out set, and declines to commit on 74 of
  132 comparisons;
* the held-out measurement is a bit-identical recurrence, so the winners it is
  scored against are the kernel's own values.

Does not establish:

* that the 56.1 per cent abstention rate is a property of the rule rather than
  of the grid.  The rate depends on density and coverage: the second sweep
  abstained 67 of 132 at the same margin, this grid 74 of 132.  A denser or
  edge-heavier grid would give another rate.
* anything about the first sweep's edge sets as such.  The held-out set is
  disjoint from them, not representative of them; it samples its own five
  rules.
* anything about A1 arguments below the model population, or about the
  unserved-order cells beyond the model's own unserved rule.
* anything about a margin other than the frozen one.  No margin scan is run
  here, by design: the held-out set is scored once, at the value the rule
  prints.
* statistical confidence.  The 132 comparisons are dependent by construction
  (band maxima over a common grid), so "58 of 58" is a record, not a
  confidence interval.

## Files

    grid.py                     rule H, with the disjointness assertion
    sweep_heldout.cpp           both entries at every order; admission
    analyze.py                  the frozen model, the self-check, the score
    provenance.py               revision by content, constants, digests; no git
    run.ps1                     the pipeline above
    data/grid.csv               1300 held-out arguments
    data/reference.csv          the 60-digit reference for them
    data/reference_validation.txt
    results/cells.csv           band,x,n,recursed,direct (42,900 cells)
    results/admission.csv       replica admission of the shipped ladder
    results/prediction_heldout.csv  the 132 scored rows
    results/per_order.csv       measured band maxima and winners per order
    results/band_summary.csv    measured band worsts and the score per band
    results/summary.txt         the self-check, the score, every commitment
    results/provenance.json     kernel v5.1.1 at 834bdbb, every recorded file
                                digest matched; git_used = false
