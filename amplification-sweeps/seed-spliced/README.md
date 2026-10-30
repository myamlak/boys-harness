# Seed-spliced sweep

The per-order comparison in Paper A sets the recursed entry's delivered error
against the direct entry's, but the two entries start from different stored
fits.  A referee's objection is that the gap therefore mixes two things: the
recurrence's amplification of its seed error, and the different accuracy of
the two fit families.  This sweep removes the second: it seeds the upward
recurrence with the direct entry's own order-0 value and repeats the band-level
comparison.  The seed error and the direct error then come from the same fit
family, and what remains is amplification and accumulated arithmetic.

## Design

The grid is the second sweep's 836 arguments, copied unchanged in
`data/grid.csv` (A1 554, A2 282).  The reference is regenerated here by
`../common/reference.py` at 70 working digits, 60 printed, and the run asserts
it is byte-identical to `../amplification_form/data/reference.csv`; the two
sweeps are read against one reference.

For each argument the driver evaluates

* the shipped recursed entry, `boys::BoysAllOrders<PaperPolicy>`;
* the direct matrix-product entry, `boys::BoysRegionAProduct<kFp64>`, once per
  band;
* the spliced ladder: `F_0` is the direct entry's own order-0 value and
  `F_l = boys::detail::UpwardStep<PaperPolicy::kDivision>(l, F_{l-1}, x,
  1/x, exp(-x)/2)` for `l = 1 .. 32`.  The step is the kernel's own body.

Every step is also run through the replica body of `../common/replica.hpp`
and compared bit for bit; the replica admission of the shipped ladder is the
same check the other sweeps record.

The calibration is the decision model's, on its population (arguments at or
above `kTierThresholds[0]` = 1.0855252345349333, where the entry's order-0
value really is a seed error): for every order and every combination of band
and metric the analysis forms

    eps0(x) = |direct_0(x) - F_0 ref|          the same-family seed error
    S(x)    = eps0(x) A_n(x)                   the seed term, A_n = (2n-1)!!/(2x)^n
    D(x)    = |direct_n(x) - F_n ref|          the direct delivered error
    E(x)    = |spliced_n(x) - F_n ref|         the spliced delivered error

and takes the band maximum of each, as the decision rule compares band maxima.  The
seed term orders the comparison

    predicted spliced  if S < D   (the seed term is below the direct error)
    predicted direct   if S > D

and it is scored against the measured winner, spliced if E < D and direct if
D < E.  Relative errors divide each quantity by the reference at the order.

Order 0 is structural: A_0 = 1, so S_0 = eps0 = D_0 = E_0 at every argument,
and the comparison is a tie between the same measured quantity.  The four
order-0 band-metric rows are recorded in `results/per_order.csv` and excluded,
leaving 128 comparisons.

The full 836-argument grid is not measurable and is not scored.  x = 0 is in
the grid; A_n(0) is not finite for n >= 1; outside the population 103 A1 cells
already have a non-finite spliced value.  Every band maximum there would be an
overflow artefact.  The population is the decision rule's own domain.

## Commands

Run from a Visual Studio developer shell, without `python -I` (mpmath lives in
the user site).  The pipeline does not run git.

    powershell -ExecutionPolicy Bypass -File run.ps1

The steps it runs, in order:

    python provenance.py --boys-tree C:\Users\myaml\dev\_boys-v511 \
        --harness sweep_seed --harness-dir . --flags "<compile line>" \
        --out results/provenance.json
    python ../common/reference.py --boys-tree C:\Users\myaml\dev\_boys-v511 \
        --grid data/grid.csv --out data/reference.csv \
        --validation data/reference_validation.txt --dps 70 --digits 60
    cl /nologo /std:c++latest /O2 /EHsc /fp:precise /arch:AVX2 \
        /I C:\Users\myaml\dev\_boys-v511\include \
        /Fo build\sweep_seed.obj /Fe build\sweep_seed.exe sweep_seed.cpp
    build\sweep_seed.exe --grid data/grid.csv --cells results/cells.csv \
        --admission results/admission.csv --steps results/replica_steps.csv
    python analyze.py --tree . --provenance results/provenance.json

`run.ps1` also fails the run if the regenerated reference differs by one byte
from the second sweep's.

## Results

Admission, on this run: upward 572 arguments and 9,200 cells bit-identical,
downward 263 arguments and 8,679 cells bit-identical, exact 1 argument and 33
cells bit-identical.  The spliced ladder's steps agree with the replica body
on all 26,720 steps over all 836 arguments, so the recurrence this sweep ran is
the kernel's own.

Seed-term comparison, orders 1-32, on the population (A1 290 arguments, A2 282):

| band | metric | ordered | ordered correctly |
| --- | --- | --- | --- |
| A1 | absolute | 32 | 32 |
| A1 | relative | 32 | 32 |
| A2 | absolute | 32 | 26 |
| A2 | relative | 32 | 29 |
| all | | 128 | 119 |

The seed term orders every one of the 128 comparisons (no ties above order 0),
and it orders 119 correctly.  The nine it gets wrong are all in A2:
absolute n = 1, 3, 4, 5, 6, 7 and relative n = 1, 2, 14.  In each the seed
term is below the direct error, so the rule predicts the spliced ladder is
better, while the measured spliced band maximum is above the direct one.  This
is the accumulated-arithmetic term the seed-only rule drops: at those orders
the A2 seed term is `2e-19 .. 5e-18`, an order or two below the delivered
spliced error.

At band level A1 is right at every order, but for a reason worth stating
plainly: from order 8 upward the A1 maxima of both S and E are attained at
x = 1.0855, the population's lower edge, where A_32 is about 1e33 and the
spliced ladder delivers 1.5e17.  There the seed term predicts "direct" because
recursing from that seed is catastrophic, and the ladder is indeed
catastrophic.  The A1 agreement is therefore dominated by the divergence at
the small-x edge, not by a delicate ordering.  A2, whose population starts at
the band edge 5.9499, is where the comparison tests something fine.

Per cell the seed term is a lower bound rather than an upper bound on the
delivered spliced error: over 18,304 population cells with both finite
(n >= 1), E/S has quartiles 0.86, 1.83, 4.07 and maximum 919.  The order-level
sign still agrees at cell level in 16,637 of 18,304 cells (90.9 per cent).

## What this establishes and what it does not

Establishes, on the pinned v5.1.1 kernel and its own 60-digit reference,
with the seed error taken from the direct entry's own fit family:

* the seed-term crossing predicts the direction of the band-level comparison
  at every order 1-32 in both metrics, 128 of 128;
* it predicts the direction correctly in 119 of 128, and the nine failures
  are exactly the cases where a band maximum is set by accumulated arithmetic
  rather than by the seed term;
* the A1 score is carried by the small-x edge's divergence; the informative
  part of the result is A2's 55 of 64.

Does not establish:

* that the two stored fit families are equally accurate.  The experiment
  deliberately makes the seed error the direct fit's own error, so it isolates
  amplification; it says nothing about which table is better.  That is what
  the original comparison measured and what this one controls away.
* that the stored tables do not dominate the original comparison.  The
  original sweep's finding stands as a fit-versus-fit measurement; this sweep
  is the complementary measurement with the fit difference removed.
* anything about order 0: A_0 = 1 makes the comparison an identity there.
* anything about the full 836-argument grid, where the recurrence and the seed
  term both diverge as x -> 0.
* the held-out prediction: the rule and margin frozen from the existing record
  are tested in the sibling `../held-out/` sweep.

## Files

    sweep_seed.cpp              the driver (ladder, admission, replica agreement)
    analyze.py                  band maxima, the seed-term comparison, the score
    provenance.py               revision by content, constants, digests; no git
    run.ps1                     the pipeline above
    data/grid.csv               the second sweep's 836 arguments, copied
    data/reference.csv          regenerated, byte-identical to the second sweep's
    data/reference_validation.txt
    results/cells.csv           band,x,n,shipped,direct,spliced
    results/admission.csv       replica admission of the shipped ladder
    results/replica_steps.csv   per-argument kernel-versus-replica step agreement
    results/seed_term_cells.csv per-cell seed term, direct and spliced errors
    results/per_order.csv       the 132 band-metric rows and the score
    results/summary.txt         the run's own account
    results/provenance.json     kernel v5.1.1 at 834bdbb: every kernel file
                                digest matched against the second sweep's record
