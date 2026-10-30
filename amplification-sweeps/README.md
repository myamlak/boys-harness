# Amplification sweeps

Two sweeps measure how the error of the Boys-function recurrence amplifies
with order. `product-amplification/` compares the recursed entry against the
direct entry; `amplification-form/` varies the seed-fit degree and the order at
which the recurrence starts.

## Build

The drivers are single translation units. From each sweep directory:

```sh
g++ -std=c++23 -O2 -I<boys>/include sweep1.cpp -o sweep1
g++ -std=c++23 -O2 -I<boys>/include sweep2.cpp -o sweep2
```

On Windows with MSVC, use `cl /std:c++latest /O2 /EHsc /I<boys>\include` with
the same source file.

## First sweep

From `product-amplification/`:

```sh
python grid.py --out data/grid.csv
python ../common/reference.py --grid data/grid.csv --out data/reference.csv
./sweep1
python analyze.py --tree .
```

`grid.py` writes 9,352 arguments. `sweep1` runs the full grid in about four
seconds (GCC 15, `-O2`, one core) and writes `results/cells.csv`,
`results/admission.csv` and `results/replica_gain.csv`. `reference.py`
computes the 60-digit reference and is the slow step. `analyze.py` writes the
tables in `results/` and `data/documented.csv`; the analysis needs the full
grid.

## Second sweep

From `amplification-form/`:

```sh
python grid.py --out data/grid.csv
python fits.py --out data/fits.csv --validation data/fits_validation.txt
python ../common/reference.py --grid data/grid.csv --out data/reference.csv
./sweep2
python analyze.py --tree .
```

`fits.py` computes the seed fits and is the slow step. `sweep2` writes
`results/paths.csv`, `results/cells.csv` and `results/admission.csv`.

## Deposited records

The measured records behind the published tables are in each sweep's `data/` and `results/`
directories: the grids, the 60-digit references, the per-cell delivered values of both entries, the
per-order and band summaries, the prediction records and the provenances. They were produced at
Boys v5.0.0 (commit `482fd33`, recorded in `results/provenance.json`).

Four files are stored xz-compressed because of their size. Decompress them in place before running
the analysis:

```sh
xz -dk amplification-form/results/paths.csv.xz \
      product-amplification/data/reference.csv.xz \
      product-amplification/results/cells.csv.xz \
      product-amplification/results/replica_gain.csv.xz
```

Every other record is plain CSV. Each analysis reads only the records above; the remaining files in
`results/` are its outputs, deposited as the tables.

## The v5.2.0 execution check

`product-amplification/v5.2.0-check/` is the source and execution check of the
two sweeps at Boys tag v5.2.0, with the v5.2.1 source-identity note
(`v5.2.1-note.md`) beside its `README.md`. It carries the tag, commit and tree
identity, the file-by-file source comparisons against v5.1.1 and v5.2.1, the
compiler and flag line, and the comparisons of the re-run sweep, table and
integral-engine records against the v5.1.1 run's own output. The v5.1.1 sweep
records are deposited in `product-amplification/v5.1.1-9353/`.

Read its `README.md` first; it is the check's own record. The
`*-compare-*.txt`, `*-diff-*.txt` and `source-diff-*.txt` files are the
evidence, each naming the files and lines that differ. `run/run_execution_check.ps1`
holds the thirteen steps verbatim, and the other scripts under `run/`
regenerate the comparison records from the trees they name; `make_check.py`
regenerates the source manifests. The check's console logs, scratch sweep
trees, build products and Python byte-code are not carried; the two provenance
records its run wrote are in `run/sweeps/*/results/`.

## The analysis and the two later experiments

`product-amplification/analysis/` holds the scripts that generate the tables from the deposited
records, and the tables, with the regeneration recipe in its `README.md`.

`seed-spliced/` holds the seed-spliced sweep: the upward recurrence seeded with the direct entry's
own order-0 value, on the second sweep's 836 arguments, with its driver, grid, 60-digit reference,
records, analysis and provenance.

`held-out/` holds the held-out test: the frozen decision rule scored on a disjoint
1,300-argument grid, with its driver, grid, 60-digit reference, records, analysis and provenance.

## Requirements

Python 3.10 or newer with `mpmath`, and a C++23 compiler. The scripts look for
the library in a `boys` directory beside this repository; set `BOYS_SOURCE` to
use another checkout. Do not run the Python scripts with `python -I`: isolated
mode hides a user-site `mpmath`.
