# v5.2.6 identity

The certificate bounds the two region-A entries through the coefficient
tables. `boys_coefficients.hpp` is byte-identical at `v5.2.5` and `v5.2.6`:
the same working-tree digest, which is the git blob's SHA-256, and the
same blob object id (both recorded in `digests.txt`). The parsed
content digest the model checks is the same at both tags, and the coarse
Chebyshev table is unchanged.

The replay against a `v5.2.6` checkout reproduces both entries' recorded
values bit for bit: 27,588 ship and 27,588 product values, zero mismatches
with the separate downward step, and the fused variant's 208 mismatches
exactly as recorded at `v5.2.5`.

The run's records were produced from the `v5.2.5` checkout. The identity of
the tables and this replay are what carry the certificate to `v5.2.6`.
`digests.txt` records the commands and the digests; `replay-v526.log` is the
captured replay.
