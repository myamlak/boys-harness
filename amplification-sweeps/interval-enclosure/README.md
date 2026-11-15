# Q-417: the band enclosure recomputed from the deposits

This is a measurement record for Paper A's referee objection Q-417. It recomputes the
paper's band-level decision model from the deposited sweep records, replaces the
questioned lower bound `W - T` with the referee's enclosure `[L, U]`, and scores both
rules against the measured winners. Everything here was computed read-only against the
harness deposit; no paper source, no harness file and no git state was changed. All
output lives under this record directory.

## 1. What the harness actually does

The paper model is the analysis-tree revision at

`amplification-sweeps/product-amplification/analysis/sweeps/amplification_form/analyze.py`

sha256 `75745b378424098cb597e653a972040489a71bdc7c8195991fde89bc28dd994e`.
The same functions are copied in `amplification-sweeps/held-out/analyze.py`. A second,
older revision exists at `amplification-sweeps/amplification-form/analyze.py`
(sha256 `4e7454492a53cef13b7ce2b98ff10bde0529ac406c105a570193499a69077f20`); it lacks
the exponential term; its deposited `results/prediction.csv` records 67 commitments.
The certified-bound record states the
deposited model record is the 65-commitment revision, and the formulas in the task
match it, so that revision was followed (see section 7).

Per band A1/A2, metric absolute/relative and order n = 0..32, with `F_l(x)` the 60-digit
reference and `recursed(x,0)` / `direct(x,n)` the cell record:

- `eps0(x) = |recursed(x,0) - F_0(x)|` (the seed error).
- `A_n(x) = (2n-1)!! / (2x)^n`, `A_0 = 1`, with the odd double factorial formed from
  integer factorials in the code.
- `S(n,x) = eps0(x) A_n(x)`.
- `tau_n(x) = 3*2^-53 * sum_{l=1..n} F_l(x) A_n(x)/A_l(x)
   + 2^-52 * exp(-x) * sum_{l=0..n-1} A_n(x)/((2l+1) A_l(x))`, and `tau_0 = 0`.
- Population: `x >= theta_0 = kTierThresholds[0] = 1.0855252345349333`, from the
  sweep's `results/provenance.json` (33 thresholds; `thresholds[32] = 10.783587858916762`).
- Served order: `n == 0 or x >= thresholds[n]`. Inside the population an unserved
  order is not recursed: the deposited code sets the predicted recursed error to the
  measured direct error and the arithmetic charge to zero (the two entries read the
  same stored fit).
- Band maxima over the population cells: `W_n = max(S + tau)` (unserved cells use
  `dir_abs`), `D_n = max dir_abs`, `T_n = max tau` (unserved cells contribute 0).
- Rule: `gap = |W - D|`; abstain if `gap <= margin * T` with `margin = 1.0` (parser
  default); otherwise `ship` (recursed) if `W < D`, `product` (direct) if `D < W`.
  In the task's terms: abstain if `|W-D| <= T`, direct if `D < W - T`,
  recursed if `D > W + T`.
- Measured winner for the record's `right` flag: per (band, metric, n), compare the
  whole-band maximum recursed error with the whole-band maximum direct error over
  every cell of the sweep, including arguments outside the population.

I reproduced these quantities exactly enough to certify the model: on the second
sweep's model record (132 rows) there are 0 state, winner, choice and right
mismatches, with maximum relative differences 4.24e-17 (`W`), 4.74e-17 (`D`) and
3.59e-17 (`T`); on the held-out record, 0 mismatches and maximum relative differences
4.32e-17, 3.49e-17 and 4.37e-17. Each sweep's `per_order.csv` is reproduced (264
values, 0 mismatched, maximum relative difference 4.7e-17 or less). The `tau` and
`amplification` functions agree with the paper file itself to within 1.1e-70 relative.

Two identities follow and hold in every one of the 132 rows of all four experiments:
`W` computed by the deposit is exactly the upper endpoint maximum `U`, and `D` is
exactly the population-restricted direct maximum `Dpop`. What the deposit never
computes is the lower endpoint `L`.

### The referee's quantities

For each row, with the deposit's unserved convention applied at unserved cells (their
interval is the degenerate `[dir_abs, dir_abs]`):

- `L_n = max over population cells of max(0, S(n,x) - tau_n(x))`
- `U_n = max over population cells of (S(n,x) + tau_n(x))` (equals `W_n`)
- `Dpop_n = max over population cells of dir_abs` (equals `D_n`)

Enclosure rule: direct only if `Dpop < L`; recursed only if `U < Dpop`; abstain
otherwise. The deposit's `W - T` is a difference of maxima and is not `L`; it falls
below `L` in 12 of 132 rows on each first-sweep cell set, 13 on the second and 18 on
the held-out sweep. Worst relative shortfall below `L`: 0.925, 0.927 and 0.979
respectively (held-out worst: A2 absolute order 17, `W-T = 2.485e-21` against
`L = 1.192e-19`).

## 2. Inputs

All input paths are relative to the root of the repository that holds this record.

| role | file | sha256 |
| --- | --- | --- |
| paper model | `amplification-sweeps/product-amplification/analysis/sweeps/amplification_form/analyze.py` | `75745b378424098cb597e653a972040489a71bdc7c8195991fde89bc28dd994e` |
| older model revision (not used) | `amplification-sweeps/amplification-form/analyze.py` | `4e7454492a53cef13b7ce2b98ff10bde0529ac406c105a570193499a69077f20` |
| first sweep cells, 9,353 args | `boys-crossover/exact-tag-replay/run/cells-9353.csv` | `2a2c45d1d7f59b74f4d3aeca647b46e4e20e86c2fb1917767a21a725f86d7da6` |
| first sweep cells, designed 9,352 | `boys-crossover/exact-tag-replay/run/cells-9352.csv` | `dd26d0b2a4352113d9f15433ec991adfacc211607ae9586c0bb51cd12b3f3369` |
| first reference | `amplification-sweeps/product-amplification/first-sweep-9353/sweeps/product_amplification/data/reference.csv` | `c88a0af5bec51674fcbee78586911855299b1101f78c6fb8487e16fb39e5589b` |
| first thresholds | `.../first-sweep-9353/sweeps/product_amplification/results/provenance.json` | `15e4c0f3b04a8b3f335785aab9b27c3d98acd42cdbb5bf0eb083521c794abddb` |
| first per-order | `.../first-sweep-9353/sweeps/product_amplification/results/per_order.csv` | `09d2fe8ce891546f2305c129955377138fd7a1781daa7fd6ed4ffdd58c311d57` |
| second cells | `amplification-sweeps/amplification-form/results/cells.csv` | `2888ea9a76e9582bd5433e782983f42040d30684b30e786840fba6c477a296e6` |
| second reference | `amplification-sweeps/amplification-form/data/reference.csv` | `514411570e902bc2320932a8ae6200d46e11d950471994c80c2775871f9a528d` |
| second thresholds | `amplification-sweeps/amplification-form/results/provenance.json` | `64df907efb091793d2ad7ac0a771d52fb41629230d11b5fe76d91e8949faabea` |
| second per-order | `amplification-sweeps/amplification-form/results/per_order.csv` | `e1ea8f51e37367ddf1f5814e3538164fba691d6b6c6a5fb73f700463f8d1b656` |
| second model record | `amplification-sweeps/certified-bound/inputs/model-record/prediction_record.csv` | `d146ab20bfe0bbcded9f241129359dfd0414c022f437d164c59c05ed03c8b218` |
| second model record (wide) | `amplification-sweeps/certified-bound/inputs/model-record/prediction.csv` | `305b2a28392c48a6eae355a9e8cca04289c214899b165464ced917d4243702eb` |
| held-out cells | `amplification-sweeps/held-out/results/cells.csv` | `252bd8df548fe118fcc04a805d20eae480ac032f521bcb694de1e5c4c5320d94` |
| held-out reference | `amplification-sweeps/held-out/data/reference.csv` | `cd8e16761db6e47563550d3f375d0747c1474d10d160ce77388656d097e7c09d` |
| held-out thresholds | `amplification-sweeps/held-out/results/provenance.json` | `91bd7ed9fa51522a012445b3c40583be17050a5245dd8e0918a5c9bd11d6e07d` |
| held-out per-order | `amplification-sweeps/held-out/results/per_order.csv` | `29c1db32cfb2d08cfd281778de97618fe4ae239aeabccd31327b18ce3777274a` |
| held-out model record | `amplification-sweeps/held-out/results/prediction_heldout.csv` | `7cb9786ea255605f3ab934f061bd8017e9564ed68fd6372c765584d513035c8d` |
| charge record | `amplification-sweeps/certified-bound/results-underflow-fix/harness_cells.csv` | `59fd3052b2d0939fb51edbdd5b9c04280f99d76fb5a72fde2b28ca132a009e72` |
| charge log | `amplification-sweeps/certified-bound/results-underflow-fix/harness.log` | `3e9f834e69bf87e3dbca3e5187d7bfbbdc1af30eb088db16cf06515fe29d95a8` |

The designed `cells-9352.csv` digest is the digest the strengthening README records for
the first sweep's deposited cell set, and `cells-9353.csv` / the second cells digest
are recorded as reproduced at tag v5.3.1 in the exact-tag-replay README. Cell counts
and population arguments: first 9,353 - 308,649 cells, 3,499 (A1) and 4,541 (A2)
population arguments; first 9,352 - 308,616 cells, 3,499 and 4,540; second - 27,588
cells, 290 and 282; held-out - 42,900 cells, 483 and 650.

## 3. Commands, environment, exit codes

Run from `amplification-sweeps/interval-enclosure` (the record's own directory)
with `PYTHONDONTWRITEBYTECODE=1 python -B ...`, Python 3.14.7 (MSC v.1944 64-bit),
mpmath 1.4.1; references are read as mpmath `mpf` at 70-digit working precision. The
`PYTHONDONTWRITEBYTECODE`/`-B` setting keeps the read-only harness free of
`__pycache__` (mpmath is in user site, so `-I` is not usable). Timings are wall clock.
`sh reproduce.sh` runs the eight commands below in order. The scripts derive the
repository root from their own location; set `BOYS_HARNESS` when the record is copied
elsewhere. The deposited inputs were produced at boys `v5.3.2`; no checkout of the
boys library is read by these scripts.

| log | command | exit | time |
| --- | --- | --- | --- |
| `logs/first9353.log` | `python -B work/q417_interval.py --case first9353` | 0 | 29 s |
| `logs/first9352.log` | `python -B work/q417_interval.py --case first9352` | 0 | 30 s |
| `logs/second.log` | `python -B work/q417_interval.py --case second` | 0 | 2 s |
| `logs/heldout.log` | `python -B work/q417_interval.py --case heldout` | 0 | 5 s |
| `logs/containment-second.log` | `python -B work/q417_containment.py --case second` | 0 | 2 s |
| `logs/containment-heldout.log` | `python -B work/q417_containment.py --case heldout` | 0 | 2 s |
| `logs/containment-first9353.log` | `python -B work/q417_containment.py --case first9353` | 0 | 16 s |
| `logs/checks.log` | `python -B work/q417_checks.py` | 0 | 1 s |

`logs/commands.log` is the machine-written copy of the same lines.
`logs/python-env.log` records the interpreter and mpmath versions. The scripts are
`work/q417_interval.py` (per-cell quantities and both rules),
`work/q417_containment.py` (the containment diagnostic) and `work/q417_checks.py`
(assembly of `results.csv`, `results.json` and
`outcome_summary.json`).

## 4. Results

Each experiment has 132 rows (2 bands x 2 metrics x 33 orders). "Direct" is the
deposit's `product`, "recursed" its `ship`. The whole-band winner is the deposit's
own winner; the restricted winner is the same comparison restricted to the model
population `x >= theta_0`. There are no ties among commitments.

| experiment | rule | commits | direct | recursed | abstains | right, whole band | wrong, whole band | right, restricted | wrong, restricted |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| first 9,353 | deposited | 68 | 61 | 7 | 64 | 67 | 1 | 68 | 0 |
| first 9,353 | enclosure | 14 | 3 | 11 | 118 | 13 | 1 | 14 | 0 |
| first 9,352 | deposited | 68 | 61 | 7 | 64 | 67 | 1 | 68 | 0 |
| first 9,352 | enclosure | 14 | 3 | 11 | 118 | 13 | 1 | 14 | 0 |
| second | deposited | 65 | 59 | 6 | 67 | 65 | 0 | 65 | 0 |
| second | enclosure | 14 | 2 | 12 | 118 | 14 | 0 | 14 | 0 |
| held-out | deposited | 58 | 53 | 5 | 74 | 58 | 0 | 58 | 0 |
| held-out | enclosure | 13 | 3 | 10 | 119 | 13 | 0 | 13 | 0 |

The single whole-band mistake of the deposited rule (both first-sweep cell sets) is A1
absolute order 0: direct is chosen (`Dpop = 1.3117e-16 < L = U = W = 1.3755e-16`,
`T = 0`), but the whole-band winner is recursed because arguments below `theta_0`
enter its maximum; restricted to the population the direct entry wins and the
commitment is right. The enclosure makes the same commitment on the same row and the
same single whole-band mistake. The model's own `right` flags are scored against the
whole-band winner, so the deposited rule's paper score is the whole-band column.

### Enclosure versus deposit, row by row

| experiment | deposited commits | justified by the enclosure | not justified | enclosure commits where the deposit abstains |
| --- | --- | --- | --- | --- |
| first 9,353 | 68 | 10 | 58 | 4 |
| first 9,352 | 68 | 10 | 58 | 4 |
| second | 65 | 8 | 57 | 6 |
| held-out | 58 | 8 | 50 | 5 |

- Every not-justified deposited commitment is a **direct** commitment at order `n >= 1`
  (58 first, 57 second, 50 held-out; 165 in total). The exact rows are in
  `results.csv` with `dep_justified_by_enc = 0`.
- Every **recursed** deposited commitment is justified: 7, 6 and 5 respectively, 18 in
  total.
- Every enclosure commitment that the deposit abstains on is **recursed**: 4 first
  (A1 relative 1, A2 relative 3, A2 relative 10, A2 relative 17), 6 second (A1 relative
  1, A2 relative 3, A2 relative 4, A2 relative 9, A2 relative 10, A2 relative 17), 5
  held-out (A1 relative 1, A2 relative 4, A2 relative 9, A2 relative 10, A2 relative
  17). All 15 are right against both winner definitions.
- The justified deposited commitments are the four order-0 rows plus 6, 4 and 4 model
  rows: first - A1 absolute 0 direct, A1 relative 0 recursed, A2 absolute 0 direct,
  A2 relative 0 direct, and A2 relative 4..9 recursed; second - A1 absolute 0
  recursed, A1 relative 0 recursed, A2 absolute 0 direct, A2 relative 0 direct, and
  A2 relative 5..8 recursed; held-out - A1 absolute 0 direct, A1 relative 0 recursed,
  A2 absolute 0 direct, A2 relative 0 direct, and A2 relative 5..8 recursed.
- No row has the two rules committing in opposite directions; the difference is always
  "both commit the same way", "enclosure abstains" or "deposit abstains".

### Population-score pairs (the referee's second score)

Right/committed against the whole-band winner and against the same-population winner:

| experiment | deposited, whole band | deposited, restricted | enclosure, whole band | enclosure, restricted |
| --- | --- | --- | --- | --- |
| first 9,353 (control) | 67/68 | 68/68 | 13/14 | 14/14 |
| first 9,352 (control) | 67/68 | 68/68 | 13/14 | 14/14 |
| second | 65/65 | 65/65 | 14/14 | 14/14 |
| held-out | 58/58 | 58/58 | 13/13 | 13/13 |

The first-experiment control pair (68 commitments, 68 correct against restricted
winners and 67 of 68 against whole-band winners) is reproduced on both deposited
first-sweep cell files. Restricting the winner to the population changes 14 winners
on the first sweep (all A1 absolute, orders 0, 10, 17-19, 21, 23-27, 29-31; the list
matches `strengthening/population_winners.log` exactly), 15 on the second (all A1
absolute, with orders 22, 27 and 29 becoming ties) and 8 on the held-out sweep (all
A1 absolute, with order 9 becoming a tie). Every changed winner on the second and
held-out sweeps, and all but A1 absolute order 0 on the first, is in a row where both
rules abstain; no score pair is altered by the change, and there are no ties among
committed rows.

## 5. The five referee-cited numbers

1. **Second sweep, 65 commitments = 4 order-0 measured comparisons + 61 model
   commitments (57 direct, 4 recursed).** Confirmed. The 4 order-0 commitments are A1
   absolute recursed, A1 relative recursed, A2 absolute direct, A2 relative direct;
   the 128 model rows supply 61 commitments, 57 direct and 4 recursed.
2. **Held-out, 58 commitments = 53 direct + 5 recursed, with always-direct already
   correct on 53 of the 58.** Confirmed. The 5 recursed commitments are A1 relative 0
   and A2 relative 5, 6, 7, 8; in all 58 committed rows the measured winner is the
   deposited choice (53 direct-winner rows, 5 recursed-winner rows), so always-direct
   is right on 53 and wrong on exactly the 5 recursed commitments.
3. **32 of the 74 held-out abstentions are structural (same stored fit, orders
   17-32).** Confirmed. The 32 rows are A1 orders 17-32 in both metrics; they have no
   served population cell, both predicted errors equal the same stored fit's error,
   `W = D`, and the deposit abstains. The other 42 abstentions have served cells.
4. **The heuristic charge is exceeded at five of the 130 harness entry-cells.**
   Confirmed from `results-underflow-fix/harness_cells.csv` and the log line
   `charge violations: cert 0 paper 5`. The five: A1 absolute order 6 product,
   A1 absolute order 15 product, A2 absolute order 0 ship, A2 relative order 0 ship,
   A2 relative order 1 ship. (130 = 65 commitments x 2 entry-cells.)
5. **First-experiment control pair 68/68 restricted and 67/68 whole-band.**
   Reproduced, see the population-score table.

## 6. Does the interval contain the measured error?

The model bounds magnitudes, so the check is `|measured recursed error|` against
`[max(0, S - tau), S + tau]` at every served population cell:

| experiment | served population cells | below lower | above upper |
| --- | --- | --- | --- |
| first 9,353 | 133,684 | 0 | 0 |
| second | 9,200 | 0 | 0 |
| held-out | 18,821 | 0 | 0 |

The signed error is negative at 52-58% of the cells and therefore lies below the
nonnegative lower endpoint by sign alone; that is not a failure of a magnitude bound,
and the paper compares magnitudes throughout.

At band level, the measured recursed band maximum lies outside `[L, U]` only in rows
that contain unserved population cells: 9 below / 24 above on the second sweep, 14/21
on held-out, 8/28 on first 9,353; every one of those rows has unserved cells. Two
facts make this airtight: every fully served row's per-cell magnitudes are inside
their intervals, so its maximum is in `[L, U]`; the rows that fail therefore all
contain unserved cells. The cause is the deposit's
unserved convention, which asserts that both entries read the same stored fit
(`rec_abs := dir_abs`). The printed values are bit-identical at only 2,670 of 9,676
unserved population cells on the second sweep, 4,381 of 18,568 on held-out and 31,807
of 131,636 on first 9,353, so the convention is not exact in the records. This is a
model-convention caveat on every row containing unserved cells, not a failure of the
served-cell algebra; the enclosure rule itself was applied with the deposit's
convention, as the task requires.

The enclosure rule's own commitments are still consistent with the measured data at
the endpoint each decision uses, including in rows with unserved cells: for every
recursed commitment the measured recursed band maximum is at or below `U`, and for
every direct commitment it is at or above `L` - 11 of 11 recursed and 3 of 3 direct
on first 9,353, 12 of 12 and 2 of 2 on the second, 10 of 10 and 3 of 3 on held-out.
So the decisions reported in section 4 are consistent with the measured data at the
endpoint each decision uses, even where the full interval does not bracket the
measured band maximum.

## 7. What the deposit does not supply

- **No paper-model first-sweep prediction record is deposited.** The first-sweep
  deposit predates the paper rule; the strengthening deposit
  (`strengthening/out/first`, the older model revision) has 69 commitments, not the
  paper model's 68. The first-sweep control pair in this record is therefore
  recomputed from the deposited cells-9353 / cells-9352 and the deposited reference,
  not copied from a record. Its 14 changed winners agree with the list in
  `strengthening/population_winners.log`; that log scores the older revision (69
  commitments), so only its winner-change list was used as a cross-check.
- **Population-restricted winners are not deposited as a record.** They are
  recomputed from the deposited cells and reference; the whole-band winners and the
  deposited `right` flags are checked against the certified-bound and held-out model
  records (section 1).
- **The paper model exists in two revisions in the harness.** The sweep's own
  `amplification-sweeps/amplification-form/analyze.py` lacks the exponential term and
  gives 67 second-sweep commitments, while the certified-bound model record and the
  analysis-tree `analyze.py` give 65. This record follows the analysis-tree revision
  and the model record, which match the task's `tau` formula. The older revision was
  not used anywhere.
- Everything else the task asks for is present in the deposits and was reproduced
  within `4.8e-17` relative or better.

## 8. Files in this record

- `results.csv` - 528 rows (4 experiments x 132), one row per comparison, with
  population/served/unserved counts, `W, D, T, W-T, L, U, Dpop`, both states, both
  winners, both right flags and the justification flag. Column names are the ones used
  above.
- `results.json` - the same data plus per-case validation (per-order
  reproduction, formula self-check, deposited-record check), the five checks with
  their confirming evidence, the charge detail, the containment diagnostics and the
  input digests.
- `outcome_summary.json` - the compact counts, justification lists, structural
  rows and the five-check booleans.
- `work/` - the three scripts and their per-case outputs (`<case>.json`,
  `<case>.csv`, `<case>.containment.json`).
- `reproduce.sh` - runs the eight commands in section 3 in order.
- `logs/` - every run above, plus `logs/python-env.log`, `logs/probe_model.log` (the
  first reproduction of the deposited model) and `logs/commands.log`.

Nothing outside this record directory was written, and no git command
was run.
