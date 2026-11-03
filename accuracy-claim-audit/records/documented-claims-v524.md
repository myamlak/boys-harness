# The v5.2.4 documented-claims enumeration

`documented-claims-v524.tsv` is the paper's enumeration of the quantitative
accuracy claims the `boys` documents state at tag `v5.2.4`. It holds 155
claims in 155 data rows. SHA-256:
`f0e2fbd3b167a52d2a00ac551e8ae4bb5ad08bbfa97f3be54a7b4c3d6e45f850`.

## Counting rule

A claim enters the enumeration if it states a bound with its domain, a
delivered figure for a named lane, route, partition, scheme or combination, an
accuracy-limiting relation, or an exactness or bit-identity guarantee. A claim
restated in more than one document counts once. The count is a hand-selected
lower bound: the contract page states admitted claims the table omits, so the
true denominator is larger than 155.

The set is the v5.0.0 audit's 154-row table with sixteen texts replaced by the
sentences the release states and one relation added:

- Replaced texts: D04, D09, D10, T03, T04, T08, S03, S06, N01--N06, N08, N09.
  The replacement changes a figure or the document that states it, not
  whether a claim of an admitted kind exists.
- Added relation: D18, the relation between the route table's two bases. The
  bar column states an absolute `3e-14` for region A and a `2.57e-14` band
  maximum for region B; the narrow table's six Chebyshev rows are judged
  against the `1e-13` contract, and its rational piece rows against the two
  absolute figures.

Three borderline readings are decided by the rule: D18 is one claim, the
narrow table's basis is that relation restated rather than a second claim, and
the F0-arm ratio `0.0962` printed beside N05 is measurement detail, not a
separate statement. The enumeration is therefore 155, not 156 or 157.

## Release split

The relation column sums to 155:

- 35 claims carried by claims-book rows.
- 2 claims, P10 and P11, measured only by their claims-book rows.
- 80 claims matched by another gate book.
- 38 claims carried by no gate row.

The 82 matched are the 80 other-book claims plus the 2 measured-only claims.
The status column gives 37 claims-book rows, 80 other-book rows and 38 none.
Three of the 155 claims are not statements of the five prose documents: P09
and P10 are source comments and P11 is a header sentence, so the strict
five-document set is 152.

The 38 no-row claims fall into seven classes: a second arithmetic build (13),
scope and domain restrictions (7), sweep readings the prose quotes (7),
quantities derived on the way to a bound (5), exactness relations (3), the
accessor's older figures (2), and one published bound certified through a
tighter base (1).

## Mapping to the gate's book

The gate's registration table is a different object: 38 entries over the 37
claims, naming 29 distinct rows. The extra entry is D01, registered with two
rows. Five of the book's 38 rows lie outside the registration and carry no
documented statement; four more are supporting rows, and two of those four
state no enumerated claim. The active set is 33 rows (29 named plus 4
supporting). The 35 carried claims are carried by 27 active rows;
`native.half.ceiling` carries ten (P07 and P08a--P08i), and the other 25 are
carried by 26 rows. The 29 registered rows are 27 carrying plus the 2 measured
only. The gate's printed RESULT, 38 of 38 claims met, counts row verdicts, not
documented coverage.

## Provenance of the as-found reading

At version 5.0.0 the pre-repair book registered 32 rows, all met, against 37
of the 154 claims. Read back against the sentences at that version, the
registrations carried 16 of the 37, 2 more were measured only by another book
and 19 by no row. Over the 154 claims the as-found reading is 16 carried, 82
matched and 56 with no row. That split is the audit's record, printed in Table
S14's caption of the paper.

It is not recomputable from the deposited v5.0.0 enumeration: that table's
`carried` column is the repaired reading, 35 carried, 82 matched and 37 none
over 154. The release's table adds D18, so its none count is 38 over 155.

## Reproduction

The paper generates its printed table from this file with
`make_enumeration_table.py`, which asserts the 155 row count, the status
counts 37/80/38, the carried and measured-only counts 35/2, the split 35/82/38
and the seven class counts above. The file's digest fixes its content.
