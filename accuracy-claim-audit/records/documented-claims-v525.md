# The v5.2.5 documented-claims enumeration

`documented-claims-v525.tsv` is the paper's enumeration of the quantitative
accuracy claims the `boys` documents state at tag `v5.2.5`. It holds 155
claims in 155 data rows. SHA-256:
`85e0f7d9950f896bfdec7b73444691e7faea475ddd0e66d64db8c8cbd4dc495f`.

## What changed from the v5.2.4 enumeration

The set is `documented-claims-v524.tsv` with two rows changed and nothing else:
T04 and N08, both the rational route's region-B delivered figure in
`lane-contract.md`, now state `6.04e-15` where they stated `1.69e-14`. The
v5.2.5 release re-fits that seed, and the contract tables print the new
delivered figure against the unchanged `2.57e-14` bar. Every id, document,
location, status, registration and carried flag is unchanged.

The counting rule and the release split are therefore those of
`documented-claims-v524.md`, and are not restated here: 155 claims, 35 carried
by claims-book rows, 82 matched (80 other-book plus 2 measured only) and 38
carried by no row.

## Reproduction

The paper generates its printed table from this file with
`make_enumeration_table.py`, one row per claim in the file's own order. The
file's digest fixes its content.
