# First amplification sweep: 9,353 arguments

This is the first of the two amplification sweeps. Every argument of the
designed grid is evaluated at all 33 orders by both routes: the recursed
entry and the direct matrix product. A separate recurrence replica is
admitted against the shipped entry bit for bit. The records here are the
measurement behind the per-order, band-summary, gain, growth and
recovered-argument tables.

## Counts

| Quantity | Value |
| --- | ---: |
| Arguments | 9,353 |
| Band A1 | 4,812 |
| Band A2 | 4,541 |
| Designed by the printed grid rules | 9,352 |
| Recovered | 1 |
| Cells (argument x order) | 308,649 |

`data/grid.csv` has a header and 9,353 rows. `results/cells.csv` has a
header and 308,649 rows, one for every argument and order (9,353 x 33).

The recovered argument is

```text
x = 5.949924076054243 = nextafter(a, +inf)
a = kRegionA1Edge    = 5.949924076054242
```

the first double inside band A2 above the shared band boundary. It is the
last row of `data/grid.csv`, band A2, tagged `recovered`. The surviving
pre-release gain table prints it as the A2 `at_x` at fourteen orders
(n = 1, 2, 3, 5, 6 and 8-16); no printed grid rule generates it. Its 60-digit
reference row, 33 cells, admission row and 16 upward gain rows are appended
to the corresponding records. The grid and reference rows print the argument
as `5.949924076054243`; the cells, admission and gain records print it at
seventeen significant digits as `5.9499240760542431`. Both parse to the same
double. The second pre-release argument is named by no surviving record and
is not synthesised.

## Layout

The run directory mirrors the measurement tree, so the deposited scripts
resolve the relative paths they were written with:

```text
v5.1.1-9353/
  make_tables.py  make_endtoend.py  make_frontier.py
  sweeps/
    common/                 boysprovenance.py, reference.py, replica.hpp
    product_amplification/  grid.py, run.ps1, analyze.py, sweep1.cpp
      data/                 grid, 60-digit reference, validation, documented weights
      results/              cells, admission, replica gains, per-order records, provenance
    amplification_form/     the second sweep's deposited records
      data/  results/
  v5.1.1-check/             the v5.0.0 comparison records and its checker
```

`make_tables.py` reads `sweeps/product_amplification/data` and
`sweeps/product_amplification/results`. `make_endtoend.py` and
`make_frontier.py` read `sweeps/amplification_form/results`, which carries
the second sweep's 836 arguments and 27,588 cells. `run.ps1` and the Python
scripts resolve `../common` and `../../make_tables.py` from
`sweeps/product_amplification`, as they do in the measurement tree.

## Regenerating the tables

Run from this directory with Python 3. In PowerShell:

```powershell
$env:BOYS_SWEEP_TREE = (Resolve-Path .\sweeps\product_amplification).Path
python .\make_tables.py
python .\make_endtoend.py
python .\make_frontier.py
```

In a POSIX shell:

```sh
BOYS_SWEEP_TREE=sweeps/product_amplification python make_tables.py
python make_endtoend.py
python make_frontier.py
```

`make_tables.py` writes `tables/band-level.tex`, `band-ratios.tex`,
`gain.tex`, `growth.tex`, `per-order-A1-abs.tex`, `per-order-A1-rel.tex`,
`per-order-A2-abs.tex`, `per-order-A2-rel.tex` and `recovered-b2.tex`.
`make_endtoend.py` writes `tables/endtoend.tex`; `make_frontier.py` writes
`tables/frontier.tex` and `figures/per-order-ratio.dat`. Every number is read
from the records and formatted; none is retyped. `make_tables.py` must run
first because it creates `tables/`. On the deposited records the three
commands reproduce the measurement tree's `tables/*.tex` and
`figures/per-order-ratio.dat` byte for byte.

Do not run the Python scripts with `python -I`: isolated mode hides the
user-site `mpmath`. `make_tables.py` uses only the standard library;
`make_endtoend.py`, `make_frontier.py` and `analyze.py` need `mpmath`.

The raw analysis can be re-derived from the cells and the 60-digit reference:

```powershell
cd .\sweeps\product_amplification
python .\analyze.py --tree .
```

The provenance file supplies the band constants, so no library checkout is
needed. The step rewrites `results/per_order.csv`, `band_summary.csv`,
`gain_up.csv`, `growth.csv`, `decomposition.csv`, `summary.txt` and
`data/documented.csv`. The deposited files are reproduced byte for byte,
except that `summary.txt` names the reference path of the run being reduced.

## Rerunning the sweep

`run.ps1` rebuilds the harness and regenerates every deposited record. From
`sweeps/product_amplification`:

```powershell
powershell -ExecutionPolicy Bypass -File .\run.ps1 -BoysTree <clean checkout>
```

It needs MSVC `cl` and a clean library checkout at the exact revision the
provenance records; it refuses a dirty tree or any other revision. It runs
`common/boysprovenance.py`, `grid.py`, `common/reference.py`, the C++ driver,
`analyze.py` and `make_tables.py`. The deposited records make this step
unnecessary for regenerating the tables.

## The v5.0.0 comparison

```powershell
python .\v5.1.1-check\check_v500_compare.py
```

The checker compares this record with the deposited v5.0.0 records at `../`
(first sweep) and `../../amplification-form/` (second sweep). It folds line
endings and reports every shared file as `SAME` or `PREFIX`. A `PREFIX` file
carries the recovered argument's rows at its end: one grid row, one reference
row, one admission row, 33 cells and 16 upward gain rows. Every shared second
sweep record is `SAME`; its `prediction_record.csv` is new here. The checker
also verifies the counts above, the recovered argument
(`nextafter(a, +inf)`), its position and tag, and its 33 cells, admission row
and gain rows. It exits 0 when every comparison passes.

The copied files in `v5.1.1-check/` are the check's record.
`sweep-results-compare.txt` is the historical comparison of the shared result
CSVs (it predates the recovered row); `bridge-compare-v511.txt` and
`bridge-run-summary-v511.txt` compare the integral-engine runs;
`compiled-tree-compare.txt` and `include-boys-compare-before.txt` /
`include-boys-compare-after.txt` compare the compiled headers; the two
`*-provenance.json` files record the check's own runs.
`check_v500_compare.py` is this bundle's driver for the sweep-record
comparison, not a file from the measurement tree.

## Files

Every file is listed with its byte count and SHA-256. `README.md` is this
file and is the only file not listed. The copied files are byte-for-byte
copies of the measurement tree's files, except `paths.csv.xz` (see below).

### Root scripts

| File | Bytes | SHA-256 |
| --- | ---: | --- |
| `make_tables.py` | 9,206 | `2b1a9de2a550a302d0f5b40cad6fdb33e8b5eb22345a285959179ac963b7cd2a` |
| `make_endtoend.py` | 8,242 | `9fb0e77f621758c53ea4880f023efa0b8f8c837b0d0d7591c56aafbc73418cbc` |
| `make_frontier.py` | 7,518 | `e512c6b8307a13f1d5ab4e08c4e842aa35331f2110f559c3b98f394373fa9866` |

### `sweeps/common`

| File | Bytes | SHA-256 |
| --- | ---: | --- |
| `sweeps/common/boysprovenance.py` | 11,598 | `c9ebc27db6d60a0c8f1d48d738dcfcaf85b16974cd477ac86af6fe35508af86a` |
| `sweeps/common/reference.py` | 9,553 | `472714cac64787d892a90b808b3b8ae660c3842826c7540e19a6e171df0d6bbe` |
| `sweeps/common/replica.hpp` | 5,068 | `e2ae9aaafc826059981b0435d629e09399a1bcbeecc3a64b9d77981d13456f69` |

### `sweeps/product_amplification`

| File | Bytes | SHA-256 |
| --- | ---: | --- |
| `sweeps/product_amplification/grid.py` | 6,195 | `b44134bad3521056bcda37474e71424109b99b2d2f8ff4362acb43000ef23524` |
| `sweeps/product_amplification/run.ps1` | 3,878 | `b7e09e02cbe445a96f98dede5d347b70959075cd9bc9a9dea3e4006e09e20a33` |
| `sweeps/product_amplification/analyze.py` | 18,225 | `7d5caf0ac3bd8c3a0f418451a6e990915b4238ccf6b96657a30b5e20e35ba686` |
| `sweeps/product_amplification/sweep1.cpp` | 8,569 | `309bcb71e29916512be479e2e51953f03f5d7d29cf4da5d74d6f3c53d1bede38` |
| `sweeps/product_amplification/data/grid.csv` | 228,793 | `780f75744e71f34da9e42b22ef359cc56424283e67a4faf7a16ea58847fae488` |
| `sweeps/product_amplification/data/reference.csv` | 20,588,300 | `c88a0af5bec51674fcbee78586911855299b1101f78c6fb8487e16fb39e5589b` |
| `sweeps/product_amplification/data/reference_validation.txt` | 921 | `81bc77ada3883f2f5fe12fda9bbdf63c0441ce78f7c1e09758c4be198260f03e` |
| `sweeps/product_amplification/data/documented.csv` | 757 | `a922b80573e1120a145d36e4241ba373c78cef150e837eedb01e02c73c633485` |
| `sweeps/product_amplification/results/cells.csv` | 21,672,799 | `2a2c45d1d7f59b74f4d3aeca647b46e4e20e86c2fb1917767a21a725f86d7da6` |
| `sweeps/product_amplification/results/admission.csv` | 335,156 | `3003115faf814830ca95948b9cc85e8880207caecbf208529e00a31f65e5dc91` |
| `sweeps/product_amplification/results/replica_gain.csv` | 7,820,122 | `e25ac559a3d5f23164340802c152bc9933b309a70b7ad09f4f5d7cc931e2855f` |
| `sweeps/product_amplification/results/per_order.csv` | 12,878 | `09d2fe8ce891546f2305c129955377138fd7a1781daa7fd6ed4ffdd58c311d57` |
| `sweeps/product_amplification/results/band_summary.csv` | 777 | `3ed145be4d20e0b1788639f32ed4ed1b39dfae9c8086a315d87ce933b9311baf` |
| `sweeps/product_amplification/results/gain_up.csv` | 2,458 | `b36c6b745371156502175a1ba678667e9cb705df79d938c16787d9fbe06e4f48` |
| `sweeps/product_amplification/results/growth.csv` | 2,897 | `5f295e898d8240f012bcb9efbed2b95f031f35d8050f220561d4b86d14f1e391` |
| `sweeps/product_amplification/results/decomposition.csv` | 927 | `1799cab3132df2155a3e6b56798caf0022a4dbd560a5bd044c641dd9c41fa262` |
| `sweeps/product_amplification/results/summary.txt` | 3,455 | `1397f8b449204bf8dffcdcb3127b24aec646e57648b4fdd59aec4e0d0292b26d` |
| `sweeps/product_amplification/results/provenance.json` | 24,406 | `15e4c0f3b04a8b3f335785aab9b27c3d98acd42cdbb5bf0eb083521c794abddb` |

### `sweeps/amplification_form`

| File | Bytes | SHA-256 |
| --- | ---: | --- |
| `sweeps/amplification_form/data/fits.csv` | 93,903 | `ed9d8ea3cebef5bd734962fcb307a36a8729b752feb6eb0fc1a966a397d835e5` |
| `sweeps/amplification_form/data/fits_validation.txt` | 14,504 | `a719511a6fdcaa47a8d9f6a9862015cfcd87bba7f700502f81b47c64a8ca27cd` |
| `sweeps/amplification_form/data/grid.csv` | 20,389 | `5fee3f7930f52e1b982c421e84d2b6e36d97af450b8512508cd9ea5dc7b9f507` |
| `sweeps/amplification_form/data/reference.csv` | 1,825,960 | `514411570e902bc2320932a8ae6200d46e11d950471994c80c2775871f9a528d` |
| `sweeps/amplification_form/data/reference_validation.txt` | 911 | `cc2a708a2ee67c908acf5917317bbef32ad440ad795dd17fa42e41ee3accf22b` |
| `sweeps/amplification_form/results/admission.csv` | 30,166 | `d57b3bf925ae3740b28048139126ad0806fcab624e996e40586bc91503fb069f` |
| `sweeps/amplification_form/results/band_summary.csv` | 581 | `ef8a89da0b33a92ca1b652f4683660f2aebfb5c9449d299652a556f17f443f10` |
| `sweeps/amplification_form/results/cells.csv` | 1,914,945 | `2888ea9a76e9582bd5433e782983f42040d30684b30e786840fba6c477a296e6` |
| `sweeps/amplification_form/results/paths.csv.xz` | 31,971,200 | `9cdbee313de19ce1b0d10a75a9bd1a9bed7acb4f256a413b2eedac35c4905f41` |
| `sweeps/amplification_form/results/per_order.csv` | 12,816 | `e1ea8f51e37367ddf1f5814e3538164fba691d6b6c6a5fb73f700463f8d1b656` |
| `sweeps/amplification_form/results/prediction.csv` | 13,356 | `597963241019d087d62b293ccee70585e7cf99f0cd4043675b6415f06dd3e371` |
| `sweeps/amplification_form/results/prediction_record.csv` | 18,390 | `50c1e0a3b33701f18c30410599078f8f09df0b4617e66a6a31297cbde652f842` |
| `sweeps/amplification_form/results/prediction_scan.csv` | 929 | `fdf681168a27e7f1dd1f06944944db0ee2511eec602fb48db986fe7ac5d1c85d` |
| `sweeps/amplification_form/results/provenance.json` | 24,538 | `a1ed6ea36dacb67163392eba46556ed2504b2eaa65e3f2122114bf9f83ff38d7` |
| `sweeps/amplification_form/results/summary.txt` | 15,594 | `8d49ea1f795860d718952583d61d078bc348a3845b390159953c880e54ffe646` |
| `sweeps/amplification_form/results/transfer_cells.csv` | 2,191,701 | `b69ff7bd14a68386e022937daf93c7b8cf35a0fff0a1f69abefb2bf12ecf907d` |

### `v5.1.1-check`

| File | Bytes | SHA-256 |
| --- | ---: | --- |
| `v5.1.1-check/amplification-form-provenance.json` | 24,538 | `a1ed6ea36dacb67163392eba46556ed2504b2eaa65e3f2122114bf9f83ff38d7` |
| `v5.1.1-check/bridge-compare-v511.txt` | 13,910 | `78bd95d379e719d05c42d071e4e10d4b367c0cfcc746de6431239311fe845bdc` |
| `v5.1.1-check/bridge-run-summary-v511.txt` | 1,571 | `217faafe2a689a165f9ca67318b79f8bb665feb866d549fb709bf5a1fc4f4f58` |
| `v5.1.1-check/check_v500_compare.py` | 10,190 | `646a271c77861849904621b57b4a7a34f70a8af39136842660e14ecf0f1eb859` |
| `v5.1.1-check/compiled-tree-compare.txt` | 2,047 | `81b436de57562076962dea966b28e2030aa34b964645e4b06bacf03c4b47de59` |
| `v5.1.1-check/include-boys-compare-after.txt` | 10,507 | `b1d47ad145ea649bfdcc3fb410a72ea7599f5acc9624bcac43bd98155305f2c4` |
| `v5.1.1-check/include-boys-compare-before.txt` | 14,949 | `1b5d8ed76c997d306c3ee0b30bc60e3b2f81368181e6c1c37e06837d1039dda1` |
| `v5.1.1-check/product-amplification-provenance.json` | 24,406 | `03c5e2f8a844d27f688402567de16ba4242db9d9a9513bc17199324ecb68725d` |
| `v5.1.1-check/sweep-results-compare.txt` | 832 | `25b5f397a9152898d9e038265bb8d0ce05614174e4c68a1501d8b6af5df200ce` |

`paths.csv.xz` is the second sweep's path record, stored compressed as in the
repository's existing second-sweep deposit. It decompresses to 291,827,555
bytes with SHA-256
`81001a06d80c88e5b3ddac1f3dd011cae6f8ce2e806c8bf1f7f26e85dd8c8219`, the
digest of the measurement tree's `paths.csv`. No table reads it.
