# The v5.3.1 documented-claims enumeration

`documented-claims-v531.tsv` is the paper's enumeration of the quantitative
accuracy claims the `boys` documents state at tag `v5.3.1`. It holds 155
claims in 155 data rows. SHA-256:
`09701aece411d1d9139ef6490522b2bc8a5a8cc76fdeec4385612b3997ebb114`.

It is deposited in this harness repository.

## What changed from the archived enumeration

The set is the previous enumeration, re-read against the documents at
`v5.3.1`. The release's accuracy statements are unchanged; the release moves
the version strings, the changelog and the documentation's run-time tier
material. Five rows were re-based or corrected:

- `P09`, `P10` and `P11` carry the current source and header anchors:
  `src/boys_simd.cpp 374-377`, `src/boys_simd.cpp 627-630`, and
  `include/boys/boys_decl.hpp 1992-1996`.
- `P09` is bound to `G22`, the individual vector-body row, and its evidence is
  that row's measurement, the worst 0.999 of budget at `(n=0, x=17.2242)`. The
  previous evidence quoted the pairwise-consistency row's 19,466 cells, which
  measure a different object.
- `H06` states the release's sentence, that the representation term is very
  nearly the whole bound at each format's worst cell. The older `x = 721`
  figure the sentence used to name is a magnitude the document no longer
  states; the row now reads the sweep's worst cell, and the evidence gives that
  cell's share.
- `H07` carries the lane-contract anchor `lane-contract half 658-661`, the
  range the sentence occupies at the release.

No claim was added or removed. Relative to the archived audit record,
seventeen texts are the release's sentences: the sixteen the repair recorded
and H06. The release split remains 35 carried by a claims-book row, 81 matched
(80 other-book plus 1 measured only) and 39 carried by no row.

## Counting rule and registration

A claim enters the enumeration if it states a bound with its domain, a
delivered figure for a named lane, route, partition, scheme or combination, an
accuracy-limiting relation, or an exactness or bit-identity guarantee. A claim
restated in more than one document counts once. The count is a hand-selected
lower bound, not a coverage rate.

The gate's registration table binds 36 documented claims in 37 entries, of
which 1 is measured without restating its sentence. The gate and its committed
run record are the governing record for that count. The gate does not compute
the carried, matched and none split; the audit reads it.

The none class has 39 claims in seven classes: a second arithmetic build (13),
scope and domain restrictions (8), sweep readings the prose quotes (7),
quantities derived on the way to a bound (5), exactness relations (3), the
accessor's older figures (2), and one published bound certified through a
tighter base (1).

## Reproduction

The paper generates its printed table from this file with
`make_enumeration_table.py`, one row per claim in the file's own order.
`accuracy-claim-audit/tools/check_registration.py` reads this file by default.
The file's digest fixes its content. The earlier enumerations are archived in
`archive/`, each with the revision it records.
