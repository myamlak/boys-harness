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

## Requirements

Python 3.10 or newer with `mpmath`, and a C++23 compiler. The scripts look for
the library in a `boys` directory beside this repository; set `BOYS_SOURCE` to
use another checkout. Do not run the Python scripts with `python -I`: isolated
mode hides a user-site `mpmath`.
