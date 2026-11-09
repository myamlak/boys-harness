# The Chebyshev recount

`recount_cheb.py` regenerates `records/chebyshev_recount.json` from the
deposited records. It re-solves the Chebyshev side of the coefficient-count
sweep with a named generator and checks every cell against the record. It
never runs the rational solvers.

## Inputs

- `records/core/sweep_all.json`, `records/core/sweep_abs_5e-14.json`,
  `records/core/sweep_rel_1e-13.json`, `records/core/sweep_rel_1e-14.json`
  and `records/ladder/sweep_all.json` carry the widths, targets, degree
  ladders and grid settings.
- `count_sweep_recorded.py` at the campaign root is the recorded sweep
  harness; the supplement's record inventory names it by digest.
- `generators/gen_boys_coefficients_v5.1.1.py` is the corrected generator.
  `generators/gen_boys_coefficients_v5.0.0.py` is the pre-repair control.

## Node convention

The source records store `settings.nodes = 2001`, the number of points on the
certifying grid. The sweep API takes an interval count:

    make_nodes(g, a, b, n)  ->  n + 1 points

so the driver passes `n = 2000`. The record's per-cell `nodes` field stores
that API parameter, not the point count. The driver asserts both numbers.

## Commands

Run from the campaign directory:

```sh
python records/chebyshev_recount/recount_cheb.py \
    --out records/chebyshev_recount/recount-run-v5.1.1.json \
    --verify records/chebyshev_recount.json \
    --report records/chebyshev_recount/verify-record-new.json \
    --metadata records/chebyshev_recount/run-metadata-new.json

python records/chebyshev_recount/recount_cheb.py \
    --generator records/chebyshev_recount/generators/gen_boys_coefficients_v5.0.0.py \
    --out records/chebyshev_recount/recount-control-v5.0.0.json \
    --verify records/chebyshev_recount.json --expect old \
    --report records/chebyshev_recount/verify-record-old.json \
    --metadata records/chebyshev_recount/run-metadata-control.json
```

The first run re-solves the 95 cells with the corrected generator and checks
the record's `new_stored`, degree and four delivered-error fields. The second
re-solves them with the pre-repair generator and checks the record's
`old_stored` counts. Both write a run record, a verification report and a
metadata file naming the exact invocation, the script and input digests, the
output digest and the environment.

`records/chebyshev_recount.json` is stored with LF line endings, and the
metadata records the committed blob's SHA-256. The Windows-written original
has the same content with CRLF and hashes to
`6ee133f9f94637de0384a77a1035c27e4e2dc4ec7a26e64394b9339d5548337ff5`.

## Files

| File | What it is |
| --- | --- |
| `recount_cheb.py` | the driver |
| `recount-run-v5.1.1.json` | the corrected run |
| `recount-control-v5.0.0.json` | the pre-repair control run |
| `verify-record-new.json` | the corrected run's verification report |
| `verify-record-old.json` | the control run's verification report |
| `run-metadata-new.json` | the corrected run's metadata |
| `run-metadata-control.json` | the control run's metadata |
| `generators/` | the two generator blobs and their digests |
| `paper-record-digests.md` | the source-record digest mapping |
| `paper-record-digests.tsv` | the same mapping as a table |
