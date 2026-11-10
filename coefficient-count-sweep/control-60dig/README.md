# Evaluation-precision control on the corrected BRASIL record

This directory is the 60-digit evaluation-precision control on the corrected
rational record, `../results-corrected/corrected-counts.json`. The record's
counts are binary64-delivered-error counts: each fit is converted exactly from
its barycentric solution to monomial coefficients at 70 decimal digits, the
coefficients are rounded to binary64, and the fit is evaluated by Horner in
binary64 by both the fused and the separate routes. At the deep targets the
fit's own error and the binary64 rounding are of comparable size, so a stored
count may be a statement about binary64 storage and evaluation rather than
about the fit. This control keeps the record's solver (BRASIL, baryrat 2.1.2),
settings, denominator cap, grid, target semantics and 70-digit F0 reference,
and changes one thing: a fit is accepted at 60-digit evaluation with its
monomial coefficients held at 60 significant digits.

The record stores each corrected cell's pair, count and delivered error but
not its fitted coefficients, so the paper width's \(10^{-15}\) cell had to be
re-solved before it could be re-evaluated at 60 digits. This directory does
that re-solve and the same control for the ship width's \(10^{-15}\) cell, the
\(W=1\) width's \(10^{-15}\) cell, the \(2^5\) width's \(10^{-15}\) and
\(10^{-14}\) cells, and the three \(10^{-8}\) endpoints.

## Identity

| Item | Value |
| --- | --- |
| Paper commit measured | `935eb662761d7c86e66247ec5881dba280c3e35c` (`main.tex` has uncommitted working-tree edits; nothing here reads `main.tex`) |
| Deposited at | this repository, `coefficient-count-sweep/control-60dig/` |
| Corrected record | `../results-corrected/corrected-counts.json`, SHA-256 `9bb8ac0633b10e8a7922684ca8a419135103f2c7825b4d06b9c71603399b7195` |
| Corrected record README | `../results-corrected/README.md`, SHA-256 `711ef926c6944231682ac3feb1e8cf368986c1f0ca9c44671e35f8dc13bb6bcf` |
| Driver | `control_60dig_brasil.py`; its SHA-256 is in every record's `meta.script_sha256` |

## Environment

The record's own verification environment: Python 3.14.7 (MSC v.1944, 64-bit),
mpmath 1.4.1, numpy 2.5.2, scipy 1.18.1, baryrat 2.1.2, gmpy2 2.3.1,
`math.fma` present. The versions are recorded in every cell record's
`meta.versions`. The scripts are run without `python -I`; the project's known
trap is that isolated mode hides the user-site mpmath. A direct import check
before the runs confirmed the versions above.

## Protocol

1. **Reproduction.** The record's own pair is re-solved at the record's stored
   count with the record's default settings (falling back to the tightened
   settings only if the default solve does not reproduce). BRASIL is
   deterministic: the re-solved fit's binary64 delivered error, computed on the
   2,001-node grid with the record's own evaluator, must equal the record's
   `rational_delivered_2001` (or `delivered_2001`) digit for digit. A cell that
   does not reproduce is a blocker and its walk does not run.
2. **60-digit scoring.** For each fit the barycentric solution is converted
   exactly to monomials \(P/Q\) at 70 decimal digits. Three errors are scored on
   the same 2,001-node grid against the same 70-digit \(F_0\) reference:
   - `e60`: coefficients rounded to 60 significant digits, evaluated at 60 dps;
   - `e64r`: coefficients rounded to binary64, evaluated at 60 dps;
   - `e64`: binary64 coefficients and binary64 Horner, both routes, the record's
     counting route.
   A fit meets at 60 digits when `e60` is at most the target and the 60-digit
   denominator guard is at least \(10^{-3}\) of its maximum. The record's
   acceptance rule otherwise applies unchanged.
3. **The walk.** Start at the record's count and confirm its pair. If it meets,
   solve every admissible pair at each count below, on the default settings
   first and on the tightened settings only when the default route found no
   meeting fit and the negative answer has to be certified. The reported count
   is the first at which no pair meets on either route; that full sweep is the
   record's below-count check, and every intervening count is swept. If the
   record's pair misses, sweep every pair at the record's count and, if none
   meets, search upward until one does. A below-count margin is the smallest
   `e60` among the failing level's admissible fits, i.e. those that pass the
   denominator guard; a fit that meets the error target but fails the guard is
   not admissible and does not bound the margin.
4. **Check grid.** The accepted fit is also evaluated on the 4,001-node grid at
   60 digits, which the record uses as a check rather than as the counting
   grid.

A count is an existence bound under this solver, not a proof of minimality,
exactly as in the records it controls.

## Commands

```sh
cd C:/Users/myaml/dev/papers-c-v5/papers/coefficient-counts/results-control-60dig
python control_60dig_brasil.py --cells A            # exit 0
python control_60dig_brasil.py --cells B            # exit 0
python control_60dig_brasil.py --cells endpoints    # exit 0
python control_60dig_brasil.py --cells C1           # exit 0
python control_60dig_brasil.py --cells C2           # exit 0
python control_60dig_brasil.py --cells all          # rebuilds summary.json from the cell records
python make_summary_table.py
python freeline_effect.py                           # model and free-line effects
python verify_cells.py                              # re-derives every accepted fit from its barycentric data
python verify_freeline.py                           # the paper's own ladder_analysis.py on every measured series
```

Run one cell at a time: the machine is shared and the deep cells are heavy.
`--force` re-runs a cell; without it a completed cell record is reused, which is
how `--cells all` rebuilds `summary.json`.

## Results

### Reproduction of the record

Every target cell was re-solved at the record's stored count and settings and
its binary64 delivered error was compared with the record's:

| Cell | Width | Target | Route | Pair | Record delivered | Re-solved `e64` | Match |
| --- | --- | --- | --- | --- | --- | --- | --- |
| A | paper | \(10^{-15}\) | default | [14/3] | 7.850420647915349e-16 | 7.850420647915349e-16 | yes |
| B | ship | \(10^{-15}\) | tightened | [14/3] | 7.874848760489894e-16 | 7.874848760489894e-16 | yes |
| C1 | \(2^5\) | \(10^{-15}\) | default | [22/2] | 6.645812668009123e-16 | 6.645812668009123e-16 | yes |
| C2 | \(2^5\) | \(10^{-14}\) | default | [14/3] | 6.325061236452513e-15 | 6.325061236452513e-15 | yes |
| D | \(2^0\) | \(10^{-15}\) | default | [6/1] | 4.493761148520244e-16 | 4.493761148520244e-16 | yes |
| E-paper | paper | \(10^{-8}\) | default | [3/3] | 2.9259769819269127e-09 | 2.9259769819269127e-09 | yes |
| E-ship | ship | \(10^{-8}\) | default | [3/3] | 6.4417613663785954e-09 | 6.4417613663785954e-09 | yes |
| E-p5 | \(2^5\) | \(10^{-8}\) | default | [3/4] | 3.2773556954891358e-09 | 3.2773556954891358e-09 | yes |

The ship cell's default-route [14/3] solve gives `e64`
\(2.205384276962203\times10^{-15}\), which does not match the record; the
tightened-route solve matches at \(7.874848760489894\times10^{-16}\) with
guard \(0.06900624689468403\). That record's witness is therefore the
tightened-route fit, and the driver's default-then-tightened fallback is what
reproduces it. Every other cell reproduces on the default route.

### The deep cells

| Cell | Width | Target | Record | Pair | Delivered (binary64) | 60-digit count | Accepted pair | `e60` | `e64r` | `e64` | `e64`/target | Below-count best admissible | Below/target |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| A | paper | \(10^{-15}\) | 18 | [14/3] | 7.850420647915349e-16 | 15 | [7/7] default | 8.0174708439909212e-16 | 2.0183482750010927e-15 | 5.314552509131714e-15 | 5.31 | 4.941162354418863e-15 | 4.94 |
| B | ship | \(10^{-15}\) | 18 | [14/3] | 7.874848760489894e-16 | 16 | [8/7] default | 5.7072831999131074e-16 | 2.4385442918623802e-15 | 5.4364695444468656e-15 | 5.44 | 1.8907963705147082e-15 | 1.89 |
| C1 | \(2^5\) | \(10^{-15}\) | 25 | [22/2] | 6.645812668009123e-16 | 18 | [6/11] default | 7.5298738556386665e-16 | 4.225607213136467e-14 | 7.364509168873253e-14 | 73.6 | 1.940395304493386e-15 | 1.94 |
| C2 | \(2^5\) | \(10^{-14}\) | 18 | [14/3] | 6.325061236452513e-15 | 17 | [9/7] default | 1.9403953044933857e-15 | 5.8115824741051313e-15 | 5.0811976954184616e-14 | 5.08 | 3.384867370061682e-14 | 3.38 |
| D | \(2^0\) | \(10^{-15}\) | 8 | [6/1] | 4.493761148520244e-16 | 8 | [6/1] default | 2.8499005022108707e-16 | 3.6017079915113592e-16 | 4.493761148520244e-16 | 0.449 | 1.0970407294044959e-15 | 1.10 |

The accepted fits at A and B also pass the 4,001-node check at 60 digits
(`e60_4001` equals `e60` in both records), as do C1 and C2; D meets at
`e60_4001` \(=2.8592450096427592\times10^{-16}\). The deep counts
fall by 3, 2, 7, 1 and 0 coefficients respectively; the below-count check at
the next count down solves every admissible pair on the default route and
then, where the default route found no meeting fit, every pair on the
tightened route.

The mechanism is storage first and evaluation second, as in the
kernel-record control. Rounding the accepted 60-digit coefficients to binary64
and changing nothing else (`e64r`) multiplies the error by 2.5 at A, 4.3 at B,
56 at C1 and 1.26 at D; doing the evaluation in binary64 as well (`e64`) adds a
further factor 2.6, 2.2, 1.7 and 1.25. The accepted \(2^5\) \(10^{-15}\) fit meets at 60
digits with \(e_{60}=7.53\times10^{-16}\) but scores
\(7.36\times10^{-14}\), 73.6 times the target, on the paper's counting route:
the binary64 counter stores seven coefficients more than the fit needs.

### The endpoints

All three \(10^{-8}\) endpoints reproduce exactly and keep their counts at
60 digits. The below-count checks at stored 6, 6 and 7 miss by 10 to 20 times
the target on both routes.

| Cell | Width | Record | 60-digit count | Accepted pair | `e60` | `e64r` | `e64` | Below-count best admissible |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| E-paper | paper | 7 | 7 | [3/3] default | 2.9259771855337138e-9 | 2.9259771643497051e-9 | 2.9259769819269127e-09 | 1.26e-7 |
| E-ship | ship | 7 | 7 | [3/3] default | 6.4417613644231688e-9 | 6.4417613049245998e-9 | 6.4417613663785954e-09 | 2.04e-7 |
| E-p5 | \(2^5\) | 8 | 8 | [3/4] default | 3.2773557391748863e-9 | 3.2773557086059786e-9 | 3.2773556954891358e-09 | 1.03e-7 |

### The model-based residual

The paper's model-based test compares the seven-decade increment from
\(10^{-8}\) to the deep target with the model's predicted increment
(Table S7: 5.57 at the paper width, 5.71 at the ship width, 6.41 at \(2^5\)).
The 60-digit endpoint counts are the E-* cells above; the 60-digit deep
counts are the control's.

| Width | Span | Model | Record increment | Record residual | 60-digit increment | 60-digit residual | Removed |
| --- | --- | --- | --- | --- | --- | --- | --- |
| \(2^0\) | \(10^{-8}\)--\(10^{-15}\) | 3.05 | 4 | +0.95 | 4 | +0.95 | 0.00 |
| paper | \(10^{-8}\)--\(10^{-15}\) | 5.57 | 11 | +5.43 | 8 | +2.43 | 3.00 |
| ship | \(10^{-8}\)--\(10^{-15}\) | 5.71 | 11 | +5.29 | 9 | +3.29 | 2.00 |
| \(2^5\) | \(10^{-8}\)--\(10^{-15}\) | 6.41 | 17 | +10.59 | 10 | +3.59 | 7.00 |

The W=1 endpoint is the record's binary64 count 4, not a control cell; the
deep count is measured at 8 unchanged, so the W=1 residual is unchanged and
already inside the two-coefficient band. The two-coefficient rounding band is two stored coefficients; the 60-digit
paper-width residual \(+2.43\) is 0.43 outside it, the ship-width \(+3.29\)
is 1.29 outside it, and the \(2^5\) \(+3.59\) is 1.59 outside it. All three
remain outside, but the paper-width margin is now a small fraction of the
binary64 one.

### The free-line test with the measured 60-digit counts

The paper's pre-set free-line test (supplement Table S6) is reproduced from the
corrected record and re-run with the measured 60-digit counts substituted. The
corrected record's series are the eight above-floor rational counts at each
width; \(2^5\)'s \(10^{-14}\) and \(10^{-15}\) entries are both measured by this
control and both substituted.

| Width | Record series | 60-digit series | Record residual / 2nd diff. | 60-digit residual / 2nd diff. | Record verdict | 60-digit verdict |
| --- | --- | --- | --- | --- | --- | --- |
| \(2^0\) | 4, 5, 5, 5, 6, 6, 6, 8 | 4, 5, 5, 5, 6, 6, 6, 8 | 0.8333 / 2 | 0.8333 / 2 | fail | fail |
| paper | 7, 8, 10, 11, 12, 12, 14, 18 | 7, 8, 10, 11, 12, 12, 14, 15 | 1.7500 / 2 | 0.7857 / 2 | fail | fail |
| ship | 7, 8, 10, 11, 12, 13, 15, 18 | 7, 8, 10, 11, 12, 13, 15, 16 | 1.1667 / 1 | 0.4286 / 1 | fail | **pass** |
| \(2^5\) | 8, 11, 12, 13, 13, 15, 18, 25 | 8, 11, 12, 13, 13, 15, 17, 18 | 3.8333 / 4 | 1.0238 / 2 | fail | fail |

The ship row is an exact flip, verified with the paper's own
`ladder_analysis.py` (`verify_freeline.py`, exit 0): on the paper-format
record the record's ship entry is `max_abs_residual` 1.1666666666666679 and
`max_abs_second_difference` 1, `affine_within_tolerance` false; with 16
substituted it is 0.42857142857142705 and 1, true. The W=1 row is measured at
8 unchanged; `ladder_analysis.py` returns 0.8333333333333321 and 2, false, on
both the record and the substituted series. The global verdict stays
C-NOLAW because the paper, \(W=1\) and \(2^5\) rational rows still fail, and
the flips fall from four to three widths. The paper-width row with 15
substituted is 0.7857142857142847 and 2, still false.

### What the control changes in the paper

- The paper-width seven-decade \(+5.43\) is a binary64 number. At 60-digit
  evaluation the same record gives \(+2.43\), still outside the two-coefficient
  band but by 0.43 rather than 3.43. The paper's statement that this cell
  cannot be re-evaluated at 60 digits without re-solving it is now superseded:
  it has been re-solved, and the six-decade \(+2.22\) is not the only controlled
  paper-width statement.
- The ship-width \(+5.29\) becomes \(+3.29\).
- The \(W=1\) \(10^{-15}\) cell is now controlled: its count stays 8, pair
  \([6/1]\), with \(e_{60}=2.8499005022108707\times10^{-16}\),
  \(e_{64r}=3.6017079915113592\times10^{-16}\),
  \(e_{64}=4.493761148520244\times10^{-16}\), 0.449 times the target, and a
  below-count margin of 1.10 times the target. Its free-line row stays fail
  (0.8333/2), verified with `ladder_analysis.py`.
- The ship width's rational free-line row flips from fail (1.1667/1) to pass
  (0.4286/1) at 60-digit evaluation. The manuscript's "four of the five
  rational widths fail" is therefore a binary64 statement: at 60 digits,
  with all four informative widths controlled (\(W=1\), paper, ship, \(2^5\)),
  three of five fail (\(W=1\), paper, \(2^5\)). The
  global C-NOLAW verdict is unchanged. The ship below-count margin is 1.89
  times the target; a stronger solver lowering the ship count from 16 to 15
  would restore the fail (0.6190/2).
- The paper width's free-line row stays fail (0.7857/2); the second difference
  carries it, as in the kernel-record control.
- The \(2^5\) \(10^{-15}\) record count 25 falls to 18 at 60 digits, and the
  model residual \(+10.59\) becomes \(+3.59\): 7.00 of the 10.59 excess is
  binary64 evaluation, the largest absolute fall in the set. The accepted
  fit's binary64 score is 73.6 times the target, so the counter stores seven
  coefficients more than the 60-digit fit needs. The \(2^5\) free-line row
  fails in both bases (1.0238 / 2), so no verdict flips there.
- The \(2^5\) \(10^{-14}\) cell falls from 18 to 17, accepted [9/7] on the
  default route with \(e_{60}=1.9403953044933857\times10^{-15}\),
  \(e_{64r}=5.8115824741051313\times10^{-15}\) and
  \(e_{64}=5.0811976954184616\times10^{-14}\) (5.08 times the target); its
  4,001-node check meets. Its guard is 0.001008, only 0.8 per cent above the
  floor, and the accepted \(2^5\) \(10^{-15}\) fit's guard is 0.001153, 15 per
  cent above it: at \(2^5\) the deep fits live close to the guard, so the
  counts there depend on the guard as well as on the error. The alternative
  meeting pairs at the same counts have guards of the same order.


## Reproducibility notes

- The record deposits each cell's pair, count and delivered error, but its
  witness coefficients are 17-significant-digit renderings of 70-digit values
  (for example `witness_P` in `corrected-counts.json`). Re-evaluating those
  strings does not reproduce the record's delivered error: for \(2^5\) at
  \(10^{-14}\) the stored strings give \(6.554452133531282\times10^{-15}\)
  against the record's \(6.325061236452513\times10^{-15}\). Re-solving the
  recorded pair reproduces it digit for digit. Every coefficient used here
  therefore comes from the solver, and the records store the barycentric
  nodes, weights and values for every fit, so any fit can be reconstructed and
  re-scored independently.
- Each cell record stores the accepted fit's \(P/Q\) at 60 significant digits
  (`P60`/`Q60`, the coefficients used for `e60`) and at 70 digits
  (`P70`/`Q70`, the exact conversion), the barycentric data, the three errors,
  the counts, the settings and the full below-count sweep.

## Files

| File | Contents |
| --- | --- |
| `control_60dig_brasil.py` | the driver |
| `make_summary_table.py` | the summary table and model-residual recomputation |
| `freeline_effect.py` | the free-line test and model residuals with the measured counts substituted; writes `freeline_effect.json` |
| `verify_cells.py` | rebuilds every accepted fit from its stored barycentric data and re-derives `e60`, `e64r` and `e64` |
| `verify_freeline.py` | builds the paper-format record and runs the paper's own `../analysis/ladder_analysis.py` on the record and every substituted series |
| `fix_margins.py` | recomputes a below-count margin from admissible (guard-passing) fits only |
| `cell_<key>.json` | one cell record (reproduction, walk, accepted fit, below-count check, settings, versions) |
| `summary.json` | the compact summary of every cell run |
| `freeline_effect.json` | the free-line and model-residual effects |
| `_run_A.log`, `_run_B.log`, `_run_C1.log`, `_run_C2.log`, `_run_endpoints.log`, `_run_all.log` | the console log of each run, with the exit code on its last line |
| `_verify_cells.log`, `_verify_freeline.log`, `_freeline_run.log` | the verification and effect runs |
