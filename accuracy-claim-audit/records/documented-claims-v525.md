# The v5.2.5 documented-claims enumeration

`documented-claims-v525.tsv` is the paper's enumeration of the quantitative
accuracy claims the `boys` documents state at tag `v5.2.5`. It holds 155
claims in 155 data rows. SHA-256:
`164cd9d7456a9a74ac2fbe7cd446f1c3b1af6ab60c04adeb21cc3411589ede95`.

## What changed from the v5.2.4 enumeration

The set is `documented-claims-v524.tsv` with four rows changed and nothing
else. T04 and N08, both the rational route's region-B delivered figure in
`lane-contract.md`, state `6.04e-15` where they stated `1.69e-14`; the v5.2.5
release re-fits that seed, and the contract tables print the new delivered
figure against the unchanged `2.57e-14` bar. P09 and P10 were corrected in the
2026-11-08 erratum below. Every other id, document, location and claim text is
unchanged.

The counting rule is that of `documented-claims-v524.md`; the release split is
155 claims, 35 carried by claims-book rows, 81 matched (80 other-book plus 1
measured only) and 39 carried by no row.

## Erratum, 2026-11-08

Two rows were corrected after the v5.2.5 release, against the gate at boys
v5.3.0. P09 was registered by gate row G21; the gate binds it to G22,
`code.half_simd_cells`, the individual vector-body row, and G21 is the retained
pairwise-consistency row. The evidence column still carries the gate's
19,466-cell individual measurement, which the pairwise row prints. P10 was
recorded as matched by G22 as a measured-only row; no x86-64 gate row measures
the no-vector fallback, so it is carried by no row. The split moves from
35 / 82 / 38 to 35 / 81 / 39 over the 155 claims: one row moves from matched
to none, and P09 stays carried. No claim text and no other row changed.

## Reproduction

The paper generates its printed table from this file with
`make_enumeration_table.py`, one row per claim in the file's own order. The
file's digest fixes its content.
