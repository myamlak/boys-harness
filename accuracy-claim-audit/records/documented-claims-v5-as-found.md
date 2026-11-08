# The v5.0.0 as-found claim verdicts

**Historical record.** These are the as-found verdicts at boys `v5.0.0`
(commit `482fd3388c276091834a6c12a642d9de2481c9f2`). The document records that
release's claims book, not the current kernel state; the current kernel
reference is tag `v5.3.0`. The latest deposited enumeration is
`documented-claims-v525.md`, measured at `v5.2.5`.

`documented-claims-v5-as-found.tsv` gives the as-found verdict for each of
the 37 claims the v5.0.0 claims book registers. The current enumeration's
`carried` column is the repaired reading; it is not the as-found data.

## Method

The wording is the v5.0.0 enumeration as deposited at tag `v1.0.20`
(`accuracy-claim-audit/records/documented-claims-v5.tsv`, SHA-256
`5d39f7de62aec69a56af8cf8612822e8270588fb0854d687caa7a7bc3e6e8348`). That
table's texts for P09, P10, P11 and H07 are the audit's own wording; the
release's registration sentences for those four rows are a later repair.

The rows are the claims book of the v5.0.0 gate, tag `v5.0.0` (commit
`482fd3388c276091834a6c12a642d9de2481c9f2`),
`tests/boys_accuracy_gate.cpp` (SHA-256
`ab17902b36686f79bb51d3d2bdf9a33018a3a499458f0025cf0c6fda242bbfb`). That
book has 32 rows, numbered G1--G32 in printed order. The six rows G33--G38
in the release's registration were added after v5.0.0 and are named in the
`basis` column only where a claim's registration cites one.

A claim is **carried** when a v5.0.0 claims-book row's statement restates the
documented sentence, preserving its domain, bound and quantifier.
**Measured-only** means no claims-book statement restates the sentence but
another gate book measures the same object. **No-row** means no v5.0.0 gate
row carries the sentence as stated.

The `registered_by` column is the release's registration from the v1.0.20
table. The `v500_row` column names the G1--G32 rows that were present at
v5.0.0, so a registration to G33--G38 can be told apart from an as-found row.

## Totals

- 16 carried: D01, D07, F01, X01--X06, H01, H08, P01, P03, P06, P08a and
  P08h.
- 2 measured-only: D13 and P04.
- 19 no-row: D02, F02, H02b, H04, H06, H07, P02, P07, P08b--P08g, P08i,
  P09, P10, P11 and A03.

The totals are 16 + 2 + 19 = 37. The audit's aggregate is stated in
`documented-claims-v524.md` and `ERRATUM.md`; the paper's coverage section
names the same sixteen and the same two measured-only claims, and its
Table S14 caption records the 16/2/19 split. The repaired reading has 35
carried and 2 measured-only, and its measured-only pair is P10 and P11, not
D13 and P04.
