# SURROGATE-CONTROL

Referee control for Paper A: can a differently tuned monotone surrogate reach
the true model's coverage with zero errors?

## Answer

**Yes.** The upward-only monotone surrogate that the manuscript uses, with its
abstention margin lowered from 1 to 0.946 of the recursed bound's
accumulated-arithmetic term, commits at exactly 67 of the 132 band-metric order
comparisons and names the measured winner at all 67. It stays error-free up to
68 commitments and first errs at 69.

The 21-versus-67 contrast in the manuscript is therefore a fixed-threshold
result, not a limit of the monotone shape. The qualification is quantitative:
under matched retuning, each model at its own best error-free margin, the true
shape covers 74 comparisons and the monotone shape 68. The shape adds six
error-free decisions, not the 46 that the fixed-threshold comparison suggests.
The second surrogate, with `A_n` held at `1/(2x)`, cannot do this at any
margin: its error-free coverage is 4 commitments.

## Inputs, commands and row counts

The control reads the existing records of
`sweeps/amplification_form`. No record under the sweep was modified.

| Input | Data rows | SHA-256 |
| --- | ---: | --- |
| `results/cells.csv` | 27,588 (836 arguments x 33 orders) | `2888ea9a76e9582bd5433e782983f42040d30684b30e786840fba6c477a296e6` |
| `data/reference.csv` | 836 (60-digit reference) | `514411570e902bc2320932a8ae6200d46e11d950471994c80c2775871f9a528d` |
| `results/prediction.csv` | 132 (2 bands x 2 metrics x 33 orders) | `597963241019d087d62b293ccee70585e7cf99f0cd4043675b6415f06dd3e371` |
| `results/prediction_scan.csv` | 28 (7 margins x 2 x 2) | `fdf681168a27e7f1dd1f06944944db0ee2511eec602fb48db986fe7ac5d1c85d` |
| `results/provenance.json` | constants, 33 dispatch thresholds | `64df907efb091793d2ad7ac0a771d52fb41629230d11b5fe76d91e8949faabea` |

The model population is the arguments at or above the first dispatch threshold,
`x >= 1.0855252345349333`: A1 290 arguments, A2 282, 18,876 population cells
(A1 9,570, A2 9,306). Of the 132 comparisons, 32 are fit-versus-fit and abstain
at every margin, 4 are order 0 and commit at every margin, and 96 carry a
positive arithmetic term.

Commands, run from `C:\Users\myaml\dev\_surrogate-control` (the scripts live
outside the paper tree because this lane may write only this report into the
paper directory):

```
cd C:/Users/myaml/dev/_surrogate-control && python surrogate_control.py > control6.log 2>&1; echo "EXIT=$?"
cd C:/Users/myaml/dev/_surrogate-control && python explore.py > explore.log 2>&1; echo "EXIT=$?"
```

Both exit 0. The environment is Python 3.14.7 (MSC v.1944 64-bit) and mpmath
1.4.1. All reference comparisons run at `mp.dps = 70`, as in `analyze.py`.

## What the existing code implements, and the definitions used

`sweeps/amplification_form/analyze.py` implements the true model only. It has no
surrogate of any kind; its `amplification(n, x)` is `(2n-1)!!/(2x)^n`. The
manuscript's surrogate numbers were therefore produced by code not present in
the sweep tree, and the freeze was implemented for this control. The definition
used is the manuscript's own sentence, not the looser paraphrase.

**True model, for reference.** For a population cell at order `n`:

- when `n = 0` or `x >= kTierThresholds[n]`, the recursed bound is
  `seed(0, x) * A_n(x) + arith(n, x)`, where `seed(0, x)` is the measured
  order-0 error of the recursed value and `arith` is the accumulated-arithmetic
  term below;
- otherwise both entries read the same stored fit and the recursed bound is the
  measured direct error;
- the direct bound is the measured direct error;
- per band and metric the two predicted bounds are maxima over the population,
  and the measured winner is the band-metric winner over the whole band;
- the model abstains iff `|recursed - direct| <= margin * arith`.

**The monotone surrogate used (recorded).**

```
A_s(n, x) = A_n(x)                        for n <= n_min(x)
A_s(n, x) = A_{n_min(x)}(x)               for n >  n_min(x)
n_min(x)  = min(32, max(0, ceil(x - 1/2)))
```

This is the paper's "keeps `A_n` on the falling branch up to `n_min(x)` and
holds `A_{n_min}` at every higher order, so the curve never rises". The
threshold `n_min` is the crossing rule the analysis code already uses for the
measured argmin, `1 + max{n : 2n + 1 < 2x}`, clipped to the order range.

The literal reading "frozen at its minimum for every order", meaning
`A_s(n, x) = A_{n_min(x)}(x)` for all `n` including the falling branch, was
tested and is **not** the manuscript's surrogate: at margin 1 it gives 19
commits, 17 right and 2 wrong (and 19/19/0 if `A_0 = 1` is kept). Neither
variant reproduces 21/21. Only the falling-branch version above reproduces the
manuscript's 21 commits, 21 right and 46 lost commitments, so that is the
definition used throughout.

**The second surrogate.** `A_s(0, x) = 1` and `A_s(n, x) = 1/(2x)` for
`n >= 1`, the first nonzero-order value `A_1`. The variant that also replaces
`A_0` gives 74 commits and 22 wrong, which is not the manuscript's number;
keeping `A_0 = 1` gives 74 commits and 20 wrong, exactly the manuscript's
figure, so that is the definition used.

**The margin convention.** Only `A_n` in the recursed bound's seed term is
replaced. The accumulated-arithmetic term
`arith(n, x) = sum_{l=1..n} 3 * 2^-53 * F_l(x) * A_n(x)/A_l(x)` is computed once
with the true amplification and shared by every variant, which is the
manuscript's "every other input unchanged". The abstention rule is
`|recursed - direct| <= margin * arith`, with the same `arith`; margin 1 is the
true model's own convention. The validation below confirms this convention:
three independent manuscript numbers reproduce exactly.

## Validation of the reconstruction

| Quantity | Manuscript | This run |
| --- | --- | --- |
| True model at margin 1 | 67 commits, right at all 67 | 67 commits, 67 right, 0 wrong |
| Monotone surrogate at margin 1 | 21 commits, right at all 21 | 21 commits, 21 right, 0 wrong |
| True model's 67 commitments turned into abstentions | 46 | 46 |
| Of those 46, naming a wrong winner | 0 | 0 |
| `A_n = 1` control | 98 commits, wrong at 27 | 98 commits, 71 right, 27 wrong |
| Frozen-first surrogate `A_n = 1/(2x)` | 74 commits, wrong at 20 | 74 commits, 54 right, 20 wrong |

The reconstruction also agrees row-for-row with the stored tables: all 132 rows
of `results/prediction.csv` match on state, choice and rightness (0
mismatches), and the true model's frontier reproduces the stored
`results/prediction_scan.csv` totals at all seven pre-registered margins
(0, 0.25, 0.5, 1, 2, 4, 8), including 100/83/17 at margin 0 and 67/67/0 at
margin 1.

The surrogate's 21 commitments at margin 1 are a subset of the true model's 67,
so the 46 lost commitments are exactly `67 - 21`; the true model names the
measured winner at all 46, which is the manuscript's "none of the 46 names a
wrong winner".

## The margin scan at 67 commitments

The surrogate's critical margin is `|recursed - direct| / arith` for each
comparison; it commits iff the margin is strictly below that value. Lowering
the margin from 1 to the interval below raises the coverage from 21 to 67. The
endpoints are the exact 70-digit critical ratios of the neighbouring
comparisons, printed at 17 significant digits; the exact edges were scored
directly, and any margin strictly inside an interval gives that row.

| Margin (x arithmetic term) | Commits | Right | Wrong |
| --- | ---: | ---: | ---: |
| `0.94674209873355647 .. 0.9493424669680654` (upper neighbour) | 66 | 66 | 0 |
| `0.94674209873355647` (exact upper edge) | 66 | 66 | 0 |
| `0.946` (interior, used for the set comparison) | 67 | 67 | 0 |
| `0.94599162174460234` (exact lower edge) | 67 | 67 | 0 |
| `0.93919403467231277 .. 0.94599162174460234` (lower neighbour) | 68 | 68 | 0 |
| `0.93611279643819232 .. 0.93919403467231277` (next step) | 69 | 68 | 1 |

At margin 1 the surrogate is at 21 commits, 21 right, 0 wrong; the margin 1
value lies strictly inside the interval `[0.99999688167734362,
1.0006081127198742)`, so it is not an edge case.

**The 67 are not the true model's 67 one-for-one.** At margin 0.946 the two
sets share 64 comparisons, the surrogate commits on three that the true model
abstains on (A1 absolute n = 11, A2 absolute n = 3, A2 absolute n = 6), and the
true model commits on three that the surrogate abstains on (A1 relative n = 4,
A2 absolute n = 21, A2 relative n = 22). On the 64 common comparisons the named
winner is identical, and all 67 of the retuned surrogate's choices are right.
Per group at margin 0.946: A1 absolute 14/14, A1 relative 9/9, A2 absolute
25/25, A2 relative 19/19.

The manuscript's own robustness check also holds for the retuned set: with both
measured stores restricted to the model population `x >= 1.0855252345349333`,
no order's measured winner changes, the true model stays at 67/67, and the
retuned surrogate stays at 67/67.

## The frontier (commits versus wrong)

Every commit count from 4 to 100 is attainable; 100 is the ceiling because the
32 fit-versus-fit comparisons cannot commit. The table is complete because the
frontier is contiguous and monotone: within a wrong-count plateau each added
commitment is right.

Monotone surrogate, `A_n` frozen at its minimum:

| Wrong | Commits |
| ---: | ---: |
| 0 | 4 .. 68 |
| 1 | 69 .. 76 |
| 2 | 77 |
| 3 | 78 .. 79 |
| 4 | 80 .. 82 |
| 5 | 83 |
| 6 | 84 .. 85 |
| 7 | 86 |
| 8 | 87 .. 88 |
| 9 | 89 .. 90 |
| 10 | 91 .. 92 |
| 11 | 93 |
| 12 | 94 |
| 13 | 95 .. 96 |
| 14 | 97 .. 98 |
| 15 | 99 |
| 16 | 100 |

The error-free ceiling is 68 commitments; the 67 that match the true model's
coverage are inside it, and the first error appears at 69 (68 right, 1 wrong).
The frontier therefore shows a single clean boundary: every retuning that keeps
coverage at or below 68 is error-free, and every retuning above it is not.

Second surrogate, `A_n` held at `A_1 = 1/(2x)` with `A_0 = 1`:

| Wrong | Commits |
| ---: | ---: |
| 0 | 4 |
| 1 .. 20 | 5 .. 24 (one error added per commitment) |
| 20 | 24 .. 90 |
| 21 | 91 |
| 22 | 92 .. 96 |
| 23 | 97 |
| 24 | 98 |
| 25 | 99 .. 100 |

Its error-free ceiling is 4 commitments, the four order-0 comparisons; the
first comparison added below that margin is already wrong, and its 20 errors
are committed at high margins, so no retuning reaches 67. At 67 coverage
(margin `[1.0013738031529327, 1.0017621243446666)`) it is right at 47 and wrong
at 20; its neighbours are 66 coverage at 46/20 and 68 coverage at 48/20. At
margin 1 it is the manuscript's 74 commits, 54 right, 20 wrong.

## What the frontier shows

The referee's question separates into three statements.

1. **Coverage is recoverable.** A monotone surrogate can reach the true model's
   67-commitment coverage with zero errors, and can go one further, to 68. The
   manuscript's "21 of 132" is the count at margin 1 only; the monotone shape
   is not information-limited at 21.

2. **The recovery costs margin.** The surrogate reaches 67 only after its
   threshold drops to 0.946 of its own accumulated-arithmetic term, while the
   true model keeps margin 1. At the surrogate's retuned margin it is
   committing on predicted differences smaller than the arithmetic its own
   bound carries. The true model's zero-error guarantees are therefore stronger
   at equal coverage, even though the coverage itself is attainable.

3. **The shape's residual advantage is six, not forty-six.** Under matched
   retuning the true model's own frontier is error-free from 4 to 74
   commitments (for example 74 commits, 74 right, 0 wrong at margin 0.94, on the
   same margin interval where the surrogate reaches 68/68/0). The true shape
   therefore buys six more error-free decisions than the monotone shape, not
   46. The manuscript's sentence that "what the non-monotone shape adds is the
   46 decisions" is correct as a fixed-threshold measurement but overstates the
   shape's role as an interpretation: at the fixed threshold the 46 are
   threshold artefact, and at matched error-free tuning the difference is six.

The second surrogate supports the same reading from the other side: holding
`A_n` at `1/(2x)` is monotone but not the upward-only freeze, and it is wrong
at 20 of its 67 commitments with an error-free ceiling of 4. The retuning
rescue is specific to the upward-only freeze; it is not a property of monotone
surrogates in general.

One further point of record: the true model's own error-free ceiling of 74 is
new. The manuscript reports the true model at margin 1 (67/67) and the
pre-registered ladder, whose next point is margin 0.5 (93 commits, 81 right,
12 wrong); the 74/74/0 point at margin 0.94 lies between them and follows from
the same stored tables. No manuscript file was edited by this control.

## Records produced by this control

The control script is
`C:\Users\myaml\dev\_surrogate-control\surrogate_control.py`, with
`explore.py` for the definition check. Logs `control6.log` and `explore.log`
and the exported frontier tables (`frontier_frozen_min.csv`,
`frontier_frozen_first.csv`, `frontier_true.csv`) are in the same directory.
Nothing under `sweeps/amplification_form` was written; the only file added to
the paper directory is this report.
