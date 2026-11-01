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
