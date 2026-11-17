# Exact-tag replay: the route-comparison study at boys v5.3.2

This is the replay of the route-comparison study's measurements at the submission tag
`v5.3.2` of `boys`. It is the v5.3.2 counterpart of the deposited
`boys-crossover/exact-tag-replay` record, made at `v5.3.1`.

## What it establishes

1. Both sweep programs, rebuilt against the tagged tree, reproduce every
   deposited sweep record byte for byte: the 9,353-argument first-sweep
   cells, admission and replica-gain records; the designed
   9,352-argument records; and the second sweep's cells, admission and
   path records.
2. The certificate record's model replay re-derives all 27,588 recursed and
   27,588 direct values bit for bit with the separate downward step. The
   fused variant mismatches the same 208 recursed values as the recorded
   run, the first at `(A1, 0.029916890803630927, 0)`, one ulp.
3. The full bound validation reports 0 violations at all 55,176
   entry-cells, the recorded worst measured-to-certified ratio 0.71875 at
   `(A1, 0.0, 24, ship)`, 264/264 population band intervals inside, and
   132/132 band comparisons, 65/65 commitments and 67/67 abstentions
   overlapping. All four certificate output files are byte-identical to
   the deposited run.

The v5.3.2 tag differs from the v5.3.1 tag under `include/` and `src/`
only in `include/boys/version.hpp`; the replayed measurements are
therefore expected to be, and are, unchanged.

## Identity

| Item | Value |
| --- | --- |
| boys tag | `v5.3.2` (annotated tag object `058f10f273f37ca3b5769e9f1173c218f38fce6c`) |
| commit | `6739447d82765ae5f7f68351a99f418cd018dd9a` |
| tree | `71be8196d53c31c8ec5d9be2c9fb48453b297c4e` |
| version header sha256 | `e2a1a783004191cd76bfc952f56d3a2bece49935722c54fb8aee64553c89b27e` |
| evaluated files | 76 (62 `include/boys` headers, 14 `src` files) |
| identical to the v5.3.1 tag tree | 61 headers, 14 src, 75 of 76 |
| difference from the v5.3.1 tag tree | `include/boys/version.hpp` only |
| identical to the v5.1.1 tag tree | 45 headers, 10 src, 55 of 76 |

The kernel tree is a `git archive` of tag `v5.3.2` at
`C:\Users\myaml\dev\_a-exact-tag-replay-v532`; `BOYS` below is that
directory. All 718 blobs of the extracted tree were verified against
`git ls-tree -r v5.3.2` with 0 mismatched (`run/tag-tree-verification.log`,
exit 0; the listing is `run/ls-tree-v532.txt`). The evaluated-file
comparison against the v5.3.1 and v5.1.1 tag trees is
`run/identity-comparison.log`; `git diff v5.3.1 v5.3.2 -- include src`
reports one changed file, `include/boys/version.hpp` (2 insertions, 2
deletions).

Compiler: Microsoft C/C++ Optimizing Compiler 19.51.36260 for x64 (MSVC
14.51.36231), flags `/nologo /std:c++latest /O2 /EHsc /fp:precise
/arch:AVX2`. Python 3.14.7, mpmath 1.4.1. The harness is at VERSION
`1.0.38` (tag `v1.0.38`, commit `03ad030e`; working-tree HEAD `8305917`).

## Commands

Build (`build-tag.log`, exit 0 for both), run from `BOYS`:

    cl /nologo /std:c++latest /O2 /EHsc /fp:precise /arch:AVX2 /I %BOYS%\include ^
       /Fobuild\sweep1.obj /Febuild\sweep1.exe ^
       <papers-a-v5>\papers\boys-crossover\sweeps\product_amplification\sweep1.cpp
    cl /nologo /std:c++latest /O2 /EHsc /fp:precise /arch:AVX2 /I %BOYS%\include ^
       /Fobuild\sweep2.obj /Febuild\sweep2.exe ^
       <papers-a-v5>\papers\boys-crossover\sweeps\amplification_form\sweep2.cpp

Runs (`run-sweeps.log`, exit 0 for all three), with the record directory
as the working directory:

    %BOYS%\build\sweep1.exe --grid <paper>\sweeps\product_amplification\data\grid.csv ^
      --cells run\cells-9353.csv --admission run\admission-9353.csv ^
      --gain run\replica_gain-9353.csv --batch full
    %BOYS%\build\sweep1.exe --grid <harness>\amplification-sweeps\product-amplification\data\grid.csv ^
      --cells run\cells-9352.csv --admission run\admission-9352.csv ^
      --gain run\replica_gain-9352.csv --batch full
    %BOYS%\build\sweep2.exe --grid <harness>\amplification-sweeps\amplification-form\data\grid.csv ^
      --fits <harness>\amplification-sweeps\amplification-form\data\fits.csv ^
      --paths run\paths.csv --cells run\cells-sweep2.csv ^
      --admission run\admission-sweep2.csv

Certificate (`run\replay-tag.log`, `run\certify-tag.log`, exit 0 for both),
from `amplification-sweeps\certified-bound`:

    python -I code/replay.py --kernel <BOYS> --both-downward
    python -I code/certify.py --kernel <BOYS> --sweep ../amplification-form ^
      --model-record inputs/model-record --exp-bits exp/exp_bits.txt ^
      --out <record>\run\certificates

The paths written as `<paper>`, `<harness>` and `<record>` were absolute
paths during the run; `digests.txt` and `provenance.json` record the
inputs' digests and the environment.

## Results

Every produced record has the digest of the deposited v5.3.1 record it
replaces:

| Record | sha256 | Match |
| --- | --- | --- |
| first sweep cells, 9,353 arguments | `2a2c45d1d7f59b74f4d3aeca647b46e4e20e86c2fb1917767a21a725f86d7da6` | yes |
| first sweep admission, 9,353 | `3003115faf814830ca95948b9cc85e8880207caecbf208529e00a31f65e5dc91` | yes |
| first sweep replica gain, 9,353 | `e25ac559a3d5f23164340802c152bc9933b309a70b7ad09f4f5d7cc931e2855f` | yes |
| first sweep cells, designed 9,352 | `dd26d0b2a4352113d9f15433ec991adfacc211607ae9586c0bb51cd12b3f3369` | yes |
| first sweep admission, designed | `19f2c6d9b993f4eb880ec124795f90ef4fc1a72f95dacfe6ddb990d5b111d143` | yes |
| first sweep replica gain, designed | `d0401a0f133ae9ca5ccc03b8f827caebf780cecee730f8bd3eff30026c6bdf29` | yes |
| second sweep cells | `2888ea9a76e9582bd5433e782983f42040d30684b30e786840fba6c477a296e6` | yes |
| second sweep admission | `d57b3bf925ae3740b28048139126ad0806fcab624e996e40586bc91503fb069f` | yes |
| second sweep paths | `81001a06d80c88e5b3ddac1f3dd011cae6f8ce2e806c8bf1f7f26e85dd8c8219` | yes |
| certificate cells_cert.csv | `8041c727ac90c6d0ff9bc5552738e91ca9fd7abe09973ea0fad43756ba62fa6a` | yes |
| certificate band_cert.csv | `1ecc53f44a18b1acf867812cf00bf74e0f4842abcbc2d925991d09a734d2e368` | yes |
| certificate model_cells.csv | `d76f2737c80e9922e7e22432134ded5a1dd0eca33f925c2853fca3de25426efb` | yes |
| certificate summary.json | `75c165fed51c91cdbb4a9546a98b5d526304a40114ba37424b1c0469fae0509a` | yes |

The comparison command is `python -I compare_outputs.py`; it exits 0 and
prints `13 of 13 measurement records identical` (`run/output-comparison.log`).
The deposited files that
are still present were hashed directly and checked against the deposited
`digests.txt`. The deposited record removed its raw `run/paths.csv` after
taking the digest, and this record follows the same practice: the reproduced
file is compared by the digest recorded in `digests.txt`
(`81001a06d80c88e5b3ddac1f3dd011cae6f8ce2e806c8bf1f7f26e85dd8c8219`), and the
raw 291,827,555-byte file is not deposited.

`run/replay-tag.log` and `run/certify-tag.log` contain the same result
lines as the deposited certificate logs; their bytes are an exact prefix
of the deposited files, which append a timing tail. The build, sweep and
identity logs name this run's paths and timings, so they are not
byte-identical to the v5.3.1 record's logs.

## Notes

- The designed 9,352-argument run reads the harness copy of the grid. That
  file has LF line endings; the deposited v5.3.1 provenance records
  `da7a6c12...` for the same 9,352 rows taken from the paper tree with CRLF
  endings. The two hold the same rows and parse to the same arguments; the
  produced records are byte-identical.
- The residual lemma probe (`residual_probe.py`) was not rerun. It is not
  part of the replay's Commands section and reads its grids from the
  manuscript workspace, not from this harness checkout; it produces no
  measurement replayed here.
- No route-comparison study source file, no file of the v5.3.1 record, and no delivery
  folder file was modified by this record.

## Files

| File | Role |
| --- | --- |
| `README.md` | this record |
| `provenance.json` | machine-readable identity, inputs, runs and results |
| `digests.txt` | sha256 of the outputs, logs and checks |
| `build-tag.log` | compiler banner and both builds |
| `run-sweeps.log` | the three sweep runs |
| `run/replay-tag.log` | the certificate model replay |
| `run/certify-tag.log` | the full bound validation |
| `run/*.csv` | the reproduced sweep records |
| `run/certificates/*` | the reproduced certificate records |
| `run/identity-comparison.log` | the evaluated-file identity comparison |
| `run/tag-tree-verification.log` | the extracted tree against `git ls-tree -r v5.3.2` |
| `run/output-comparison.log` | the measurement-record digest comparison output |
| `run/ls-tree-v532.txt` | the `git ls-tree -r v5.3.2` listing compared against |
| `compare_outputs.py` | the measurement-record digest comparison |
| `compare_identity.py` | the evaluated-file identity comparison |
| `verify_tag_tree.py` | the extracted-tree verification |
| `build/sweep1.exe`, `build/sweep2.exe` | not deposited (build products); rebuilt by the two `cl` lines above, SHA-256 `fe61a0a6750f416653d8a4cd0407a22dd09477935a966c8cdf2cc9a6f71e1c2e` and `8159be35fc56207d7f00415bdedaad397f161ef07859bca88e57f9b220e1dbfb` |
