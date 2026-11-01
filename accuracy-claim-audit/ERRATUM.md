# Erratum: the documented-claims table

The deposit at tag `v1.0.13` (commit `4893447`) codes `P09`, `P10` and `P11`
with the lane-contract page. They are not statements of that page: `P09` and
`P10` are source comments and `P11` is a sentence in the shipped public
header. This deposit corrects the document column:

- `P09`, `P10`: `S`, `src/boys_simd.cpp:370-373`.
- `P11`: `H`, `include/boys/boys_decl.hpp:1913-1915`, with the corrected
  claim sentence.

The table's SHA-256 changed from
`62f81125f65c0520977f0fbd10203851f53f50cceb9c6dc82b64a4491ba7bff4`
(`accuracy-claim-audit/records/documented-claims-v5.tsv` at `v1.0.13`) to
`f014998c16f057c26e908f934d5874f978705b6eacca5c81ea708e1ee49226e4`.
The table also carries revised prose, with no document change, in `H07` and
`K06`.

The as-found aggregate split 16/2/19 is stated in the caption of Table S14 of
the paper; it is not a column of the TSV, so it is not reproducible from the
table alone.

## Follow-up correction at `v1.0.20`

The deposit at tag `v1.0.19` (commit `c380a87`) records `P11` at
`include/boys/boys_decl.hpp:1913-1915`. The claim sentence spans lines
1913-1917 at both `v5.0.0` and `v5.2.3`, so the anchor is corrected to
`1913-1917`; the rest of the table is unchanged from the `v1.0.19` table.

The table's SHA-256 changed from
`f014998c16f057c26e908f934d5874f978705b6eacca5c81ea708e1ee49226e4`
(`accuracy-claim-audit/records/documented-claims-v5.tsv` at `v1.0.19`) to
`5d39f7de62aec69a56af8cf8612822e8270588fb0854d687caa7a7bc3e6e8348`.

The helper `accuracy-claim-audit/headline_numbers.py` is replaced by the
paper's corrected version, whose labels now match the paper: Boost.Math has 2
row-shaped claims under the interpreted reading, and crlibm has 18 strict, 39
interpreted and 20 aggregate row-shaped claims.
