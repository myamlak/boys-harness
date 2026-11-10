# Held-out exact-tag replay: Paper A's held-out score at boys v5.3.1

The held-out sweep's deposited record in `../held-out/` is a v5.1.1
measurement: its `results/provenance.json` names that tag, commit and tree.
This directory replays it against the exact v5.3.1 tag tree, so the paper's
held-out score is a result at the submission release and not only at an
earlier revision.

## What it establishes

1. The 1300-argument grid, regenerated from the tag tree's band edges, is
   byte-identical to `../held-out/data/grid.csv` (sha256 `22f95bfb75f3…`) and
   the disjointness check passes against the 9541 distinct earlier arguments.
2. The kernel tree is the tagged one by content: all 62 headers and 14 source
   files hash to the v5.3.1 column of the certified-bound `v531-check` record,
   the column the exact-tag replay proved against the tag tree; the three
   reference-chain files hash to the second sweep's record, so the reference
   chain did not move between the two kernels.
3. The sweep's outputs at the tag are byte-identical to the v5.1.1 record.
   `cells.csv` and `admission.csv` compare equal, so the 21 files that differ
   between v5.1.1 and v5.3.1 change no held-out delivered value.
4. The analysis outputs are byte-identical to the deposit, and the score is
   unchanged: 58 commitments, 58 right, 0 wrong, 74 abstentions, 53 direct and
   5 recursed. `summary.txt` differs only in the line that names the
   reference's own path.

## Identity

| Item | Value |
| --- | --- |
| boys tag | `v5.3.1` (annotated tag object `d1e63f71c6f2501aecc6340e31216a42c2a390b9`) |
| commit | `2f83e40a6504368b7b2f13682e29b4f24ab2a40c` |
| tree | `344df3ecc1c13122f59a8a98a78afbff3c34337b` |
| version header sha256 | `02e6312dfa1a91b07abe98843c484f2331abb8f3f48b313fb6a15722c9067da9` |
| evaluated files | 76 (62 headers, 14 source); v5.3.1 column 76 of 76 |
| reference chain | identical to the v5.1.1 record (3 files) |
| compiler | MSVC 19.51.36260 for x64 (14.51.36231) |
| flags | `/nologo /std:c++latest /O2 /EHsc /fp:precise /arch:AVX2` |
| python, mpmath | 3.14.7, 1.4.1 |
| reference.csv sha256 | `cd8e16761db6e47563550d3f375d0747c1474d10d160ce77388656d097e7c09d` |

## Relation to `../held-out`

The analysis, grid, provenance and driver scripts are byte-identical copies.
Two things differ in what the run does. The provenance reads
`results/kernel_record_v531.json`, which carries the v5.3.1 digest column
instead of the second sweep's v5.1.1 record; `kernel_record_v531.py` writes it
and fails unless every digest checks out in the tree it is pointed at. The
reference is copied from the deposit rather than regenerated, because the
reference chain is byte-identical and its digest is checked; the grid is
regenerated and compared. The runners differ too: `run.ps1` and `build.bat`
rebuild against the tag and compare the outputs.

The record's `kernel_record_v531.json` names the digests file the run read and
that file's SHA-256. The deposited `certified-bound/v531-check/digests.txt`
carries the v1.0.35 correction note added after that run, so the whole file
hashes differently; its v5.3.1 column is identical row for row to the column
the record carries, and that column is what the record proves. A fresh run
with the deposited file writes its own `revision_source` and digest.

The `harness/` digests in `results/provenance.json` name the scripts as they
ran. The deposited copies differ only where a repository path changed (the
paths above, the driver's replica include, and the hyphenated
`amplification-form`); the replay from this deposit reproduces the five
value-bearing outputs byte for byte.

`results/comparison.txt` is the run's own output:

| output | state |
| --- | --- |
| `cells.csv` | identical, `252bd8df548f…` |
| `admission.csv` | identical, `d3110a880d33…` |
| `prediction_heldout.csv` | identical, `7cb9786ea255…` |
| `band_summary.csv` | identical, `be25d4116928…` |
| `per_order.csv` | identical, `29c1db32cfb2…` |
| `summary.txt` | two lines differ: the reference path |

## Commands

Set `BOYS_TREE` to the v5.3.1 tag checkout, then run the pipeline from this
directory in a shell with `python` on PATH. `build.bat` needs `cl` on PATH, so
start from a Visual Studio developer shell. The pipeline does not run git.

    $env:BOYS_TREE = "C:\path\to\the\v5.3.1\tag\checkout"
    powershell -ExecutionPolicy Bypass -File run.ps1

`run.ps1` performs these steps in order: it writes the kernel digest record
from the tag tree and the v531-check digests (`kernel_record_v531.py`), which
fails unless all 76 evaluated files match; regenerates the grid from the tag
tree's band edges (`grid.py`); copies the deposited reference and checks its
digest; records the kernel identity and compile line (`provenance.py`); builds
`sweep_heldout.cpp` with `build.bat`; runs it; scores the frozen rule
(`analyze.py`); and compares every value-bearing output with `../held-out`.

The analyse step reads the deposited second sweep
(`../amplification-form/results/prediction_record.csv`) and the deposited
`../../tables/prediction-record.tex` at the repository root (sha256
`f1234cff92929a30821f075f70990354d046e2166871d4790c189e7c734e2251`), through
the same default path `../held-out` uses. The run stops before writing a
held-out number if either self-check differs.

This record ships in this deposit beside the pre-tag record; the earlier
deposits are unchanged.
