# Stored-coefficient counts

`count_sweep.py` measures the smallest stored coefficient count each fit
family needs on `[x0, x0+W)` at a target. It imports the library's generator
rather than reimplementing the fits.

## Commands

```sh
python count_sweep.py --plan
python count_sweep.py --static --out static_counts.json
python count_sweep.py --sweep --out-dir results \
    --targets rel:1e-13,abs:5e-14 --widths all
```

`--plan` prints the width grid, the targets and the solver settings without
running a fit. `--static` reads the shipped header and writes what it stores.
`--sweep` measures the smallest count per family; it is resumable and stops
before starting a new row once `--max-seconds` is spent.

The stability pass re-solves every measured row at the boundary and at the
doubled settings:

```sh
python stability_save_verify.py --plan
python stability_save_verify.py --run
```

`cold_ladder.py` re-solves the one row whose verdict moved.

## Record

`records/chebyshev_recount.json` is the Chebyshev side of the sweep re-solved
with the corrected split-Clenshaw evaluator, the transcription the generator
uses to measure a fit's delivered error. Each of its 95 cells carries the
stored count in the earlier record, the corrected count and the delivered
error; the rational counts, which were not re-solved, are carried beside them.
The file states the library revision it was measured at and the content
digests of the generator, the coefficient header and the five records it was
re-solved from, together with the digest of the generator the earlier counts
came from. Run through that earlier generator, the same script reproduces all
135 recorded rows.

## Requirements

Python 3.10 or newer with `mpmath`. The library is found in a `boys` directory
beside this repository; set `BOYS_SOURCE` to use another checkout. Do not run
the scripts with `python -I`.
