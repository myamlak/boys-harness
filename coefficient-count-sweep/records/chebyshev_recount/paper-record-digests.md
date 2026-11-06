# The paper-record digests

`paper-record-digests.tsv` maps the five `paper_records` entries in
`records/chebyshev_recount.json` to the deposited files.

## Method

The `paper_records` digests are of the measurement-tree originals at Boys
v5.0.0 (commit `482fd3388c276091834a6c12a642d9de2481c9f2`). Those files were
committed to this repository with LF line endings, so their bytes hash
differently. For each file the deposited blob at the commit named in `note`
was read, each LF was rewritten as CRLF, and the result was hashed: all five
reproduce the original digest exactly. That is the content equivalence the
`note` column records.

The `deposited_sha256` column is the digest of that re-serialised blob as it
is stored in the repository. The `deposited_sha256_current` column is the
digest of the same path in the deposit at v1.0.28 (the bytes a fresh clone
gets).

## Current state

The `ladder` copy is unchanged since the recount deposit: it is still the
re-serialised original, and its 13 cells that the recount changed still carry
`old_stored`.

The four core-side files were re-measured at Boys v5.2.5 in commit
`6681f9d` ("Re-measure the records at the v5.2.5 release"). They are
therefore not content-equivalent to the originals: their Chebyshev stored
counts match the recount record's `new_stored` in every referenced cell (40
rows in `sweep_all.json`, 20 in `sweep_rel_1e-13.json`, 20 in
`sweep_abs_5e-14.json` and 15 in `sweep_rel_1e-14.json`). The recount's
`old_stored` differs from them in 35 cell references across those four
files.

The original digests remain the record of the copies the recount re-solved;
the current digests are the copies the deposit carries.
