# Exact-tag replay at the submission release v5.3.2

This is the replay of Paper A's measurements at the boys submission tag
`v5.3.2`. It sits beside two earlier records: the pre-tag record
`amplification-sweeps/certified-bound/v531-check/`, whose measurements were
made on the measured tree, and `boys-crossover/exact-tag-replay/`, the same
replay at the tag `v5.3.1`.

## What it establishes

- Both sweep programs, rebuilt against the `v5.3.2` tree, reproduce every
  deposited sweep record byte for byte: the 9,353-argument first-sweep cells,
  admission and replica-gain records; the designed 9,352-argument records; the
  second sweep's cells, admission and path records. The path record's digest
  is recorded and the 291,827,555-byte raw file was removed after it was
  verified, as in the earlier replay.
- The certificate record's model replay re-derives all 27,588 recursed and
  27,588 direct values bit for bit with the separate downward step. The fused
  downward variant mismatches the same 208 recursed values as the recorded
  run, the first at `(A1, 0.029916890803630927, 0)`, one ulp.
- The full bound validation reports 0 violations at all 55,176 entry-cells,
  with the recorded worst measured-to-certified ratio 0.71875 at
  `(A1, 0.0, 24, ship)`, 264/264 population band intervals inside, and
  132/132 band comparisons, 65/65 commitments and 67/67 abstentions
  overlapping. All four certificate output files are byte-identical to the
  deposited run.

The evaluated sources are identical between `v5.3.1` and `v5.3.2`: 76 of 76
files. The full trees differ in two files, neither of them compiled into the
sweeps: `CMakeLists.txt` registers one added test inside the
`BOYS_BUILD_TESTS` branch, and `tests/boys_structural_invariants_test.cpp` is
added. No file under `include/` or `src/` differs between the two tags.

## Identity

| Item | Value |
| --- | --- |
| boys tag | `v5.3.2` (annotated tag object `246d5392c0255c469aa42ff518bd833cbaa1e721`) |
| commit | `bb9c5a143c5253313c15f37eac37562b2577033c` |
| tree | `fc0fb26faeb7804c299bf11305b5e75aabb8f056` |
| version header sha256 | `02e6312dfa1a91b07abe98843c484f2331abb8f3f48b313fb6a15722c9067da9` |
| evaluated files | 76 (62 `include/boys` headers, 14 `src` files) |
| identical to the v5.3.1 tag | 62 headers, 14 src, 76 of 76 |
| identical to the measured v5.1.1 tag tree | 45 headers, 10 src, 55 of 76 |
| full-tree difference from v5.3.1 | `CMakeLists.txt`, `tests/boys_structural_invariants_test.cpp` |

The version header is byte-identical between the two tags; it takes the
version number from the top-level `CMakeLists.txt` at compile time.

## Commands

The boys tree is a `git archive` of tag `v5.3.2`; `BOYS` below is that
directory. The two executables are build products and are not deposited.
Rebuild them with the two `cl` lines below; their digests are in
`digests.txt`.

    cl /nologo /std:c++latest /O2 /EHsc /fp:precise /arch:AVX2 /I %BOYS%\include ^
       /Fo build\sweep1.obj /Fe build\sweep1.exe ^
       papers-a-v5\papers\boys-crossover\sweeps\product_amplification\sweep1.cpp
    cl /nologo /std:c++latest /O2 /EHsc /fp:precise /arch:AVX2 /I %BOYS%\include ^
       /Fo build\sweep2.obj /Fe build\sweep2.exe ^
       papers-a-v5\papers\boys-crossover\sweeps\amplification_form\sweep2.cpp

Runs (`run-sweeps.log`, exits 0):

    build\sweep1.exe --grid <paper>\sweeps\product_amplification\data\grid.csv ^
      --cells run\cells-9353.csv --admission run\admission-9353.csv ^
      --gain run\replica_gain-9353.csv --batch full
    build\sweep1.exe --grid <harness>\amplification-sweeps\product-amplification\data\grid.csv ^
      --cells run\cells-9352.csv --admission run\admission-9352.csv ^
      --gain run\replica_gain-9352.csv --batch full
    build\sweep2.exe --grid <harness>\amplification-sweeps\amplification-form\data\grid.csv ^
      --fits <harness>\amplification-sweeps\amplification-form\data\fits.csv ^
      --paths run\paths.csv --cells run\cells-sweep2.csv ^
      --admission run\admission-sweep2.csv

Certificate (`run\replay-tag.log`, `run\certify-tag.log`, from
`amplification-sweeps\certified-bound`):

    python -I code/replay.py --kernel <boys> --both-downward
    python -I code/certify.py --kernel <boys> --sweep ../amplification-form ^
      --model-record inputs/model-record --exp-bits exp/exp_bits.txt ^
      --out <record>\run\certificates

Identity (`run\identity-comparison.log`):

    python ..\exact-tag-replay\compare_identity.py <v5.3.2 tree> ^
      ..\..\amplification-sweeps\certified-bound\v531-check\digests.txt <v5.1.1 tree>
    python compare_tags.py <v5.3.1 tree> <v5.3.2 tree> ^
      ..\..\amplification-sweeps\certified-bound\v531-check\digests.txt

## Results

Every produced record has the digest of the deposited record it replaces:

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

The produced files themselves are not deposited: each digest equals the
deposited file of the same name under `boys-crossover/exact-tag-replay/run/`
in this repository. A reader who reruns the commands above obtains files with
these digests.

## Files

| File | Role |
| --- | --- |
| `README.md` | this record |
| `provenance.json` | machine-readable identity, inputs, runs and results |
| `digests.txt` | sha256 of the outputs, logs and scripts |
| `build-tag.log` | compiler banner and both builds |
| `run-sweeps.log` | the three sweep runs |
| `run/replay-tag.log` | the certificate model replay |
| `run/certify-tag.log` | the full bound validation |
| `run/identity-comparison.log` | the evaluated-file comparisons |
| `compare_tags.py` | the direct two-tag comparison over the evaluated set |

No paper or deposit file was modified by this record.
