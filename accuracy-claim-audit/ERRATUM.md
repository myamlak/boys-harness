# Erratum: the documented-claims table

The current table is `records/documented-claims-v531.tsv`; its front page is
`records/documented-claims-v531.md`. The earlier enumerations are archived in
`records/archive/` as provenance. This note records the corrections that the
current table carries.

## Document and anchor corrections

An earlier deposit coded `P09`, `P10` and `P11` with the lane-contract page.
They are not statements of that page: `P09` and `P10` are source comments and
`P11` is a sentence in the shipped public header. The current table records:

- `P09`: `S`, `src/boys_simd.cpp 374-377`.
- `P10`: `S`, `src/boys_simd.cpp 627-630`.
- `P11`: `H`, `include/boys/boys_decl.hpp 1992-1996`.

The anchors are the lines the sentence occupies at the submission tag, and the
gate's registration table names the same ranges for the two claims it binds
(`P09` and `P11`).

## Claim corrections

- `P09` is bound to `G22`, the individual vector-body row. Its evidence is that
  row's measurement, the worst 0.999 of budget at `(n=0, x=17.2242)`. An
  earlier table quoted the pairwise-consistency row's 19,466 cells, which
  measure a different object.
- `P10`, the no-vector scalar fallback, is bound to no row. It is read in the
  none class, under scope and domain restrictions.
- `H06` states the release's sentence, that the representation term is very
  nearly the whole bound at each format's worst cell. The older `x = 721`
  figure the sentence used to name is a magnitude the document no longer
  states; the row reads the sweep's worst cell instead, and the evidence gives
  that cell's share.
- `H07` carries the lane-contract sentence it restates, at its current anchor,
  `lane-contract half 658-661`, the range the gate registers.

## Registration and split

The gate's registration table binds 36 documented claims in 37 entries. Of
those, 35 are carried by the row the table binds to them and 1, `P11`, is
measured without restating its sentence. `P10` is bound to no row. The gate
and its committed run record are the governing record for the count; the gate
does not compute the carried, matched and none split.

Over the 155 enumerated claims the split is 35 carried, 81 matched (80
other-book measurements plus the 1 measured-only claim) and 39 with no row.
The as-found split of 16 carried, 2 matched only and 19 with no row over the
37 claims registered then is a record of the earlier rows, not a column of the
current table.

The table also carries a `carried` column. It answers one question for a
claims-book row: does the row's statement restate the documented sentence,
preserving its domain, bound and quantifier? It is blank for rows that are not
claims-book rows.

The table's SHA-256 is
`09701aece411d1d9139ef6490522b2bc8a5a8cc76fdeec4385612b3997ebb114`.
