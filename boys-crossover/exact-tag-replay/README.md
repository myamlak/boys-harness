# Exact-tag replay: Paper A at boys v5.3.1

This is the replay of Paper A's measurements at the final boys tag `v5.3.1`.
It was staged for the harness `v1.0.35` record that will carry the exact-tag
replay beside the pre-tag `v531-check` record.

## What it establishes

1. Both sweep programs, rebuilt against the tagged tree, reproduce every
   deposited sweep record byte for byte: the 9,353-argument first-sweep cells,
   admission and replica-gain records; the designed 9,352-argument records; and
   the second sweep's cells, admission and path records.
2. The certificate record's model replay re-derives all 27,588 recursed and
   27,588 direct values bit for bit with the separate downward step; the fused
   variant mismatches 208 recursed values, the first at
   `(A1, 0.029916890803630927, 0)`, one ulp.
3. The full bound validation reports 0 violations at all 55,176 entry-cells,
   the recorded worst measured-to-certified ratio 0.71875 at
   `(A1, 0.0, 24, ship)`, 264/264 population band intervals inside, and
   132/132 band comparisons, 65/65 commitments and 67/67 abstentions
   overlapping. All four certificate output files are byte-identical to the
   deposited run.

The pre-tag `v531-check` record states that its build tree differs from the
tagged tree only in `include/boys/version.hpp`. That statement is confirmed
here: the final tag tree matches the record's v5.3.1 digest column for all 76
evaluated files, and a direct comparison with the v5.3.0 tag shows the same
single-file difference. The pre-tag replay therefore remains valid as
supporting evidence; this record is the replay at the cut tag itself.

## The residual lemma check

`residual_probe.py` runs the same model over the certificate grid and the
first experiment and instruments the refined division. On the 836-argument
certificate grid it finds 8,628 upward and 8,416 downward calls, of which 45
upward and 744 downward return an exactly-zero computed residual; no nonzero
residual is subnormal; the smallest nominal numerator is $1.9\times10^{-6}$
upward and $0.1747$ downward; the smallest nonzero residual is
$3.1\times10^{-25}$ upward and $4.3\times10^{-19}$ downward; and the largest
recurrence radius is $3.7739\times10^{-14}$. The first experiment gives the
same minima with 45 and 3,596 zeros. These are the measurements behind
Supplement~S18.4's computed-zero lemma; the lemma is derived, and the probe
checks its premises on the grids. `run/residual-probe.log` is the output. The
script reads the two grids from the manuscript workspace's
`papers/boys-crossover`; the same grids are deposited under
`amplification-sweeps/`, and the recorded output is the check's result.

## Identity

| Item | Value |
| --- | --- |
| boys tag | `v5.3.1` (annotated tag object `d1e63f71c6f2501aecc6340e31216a42c2a390b9`) |
| commit | `2f83e40a6504368b7b2f13682e29b4f24ab2a40c` |
| tree | `344df3ecc1c13122f59a8a98a78afbff3c34337b` |
| version header sha256 | `02e6312dfa1a91b07abe98843c484f2331abb8f3f48b313fb6a15722c9067da9` |
| evaluated files | 76 (62 `include/boys` headers, 14 `src` files) |
| identical to measured v5.1.1 tag | 45 headers, 10 src, 55 of 76 |
| identical to the pre-tag record's v5.3.1 column | 62 headers, 14 src, 76 of 76 |
| difference from the v5.3.0 tag | `include/boys/version.hpp` only |

The 45/10/55 count is a direct byte comparison with the v5.1.1 tag tree
(`run/identity-comparison.log`). The pre-tag record's `digests.txt` first
column also gives 45/10/55 after the correction carried at harness `v1.0.35`.
In the `v1.0.34` deposit that column gave 47/10/57: two rows,
`boys_build_defaults.hpp` and `boys_coefficients.hpp`, carried the v5.3.1
digest in the measured column instead of the measured-tree digest. The two
entries are corrected in `amplification-sweeps/certified-bound/v531-check/`,
whose README and `digests.txt` state the correction; the corrected column
matches the v5.1.1 tag tree for 45 of the 62 headers, 10 of the 14 src files
and 55 of the 76 files. The v5.3.1 column is unchanged and matches the tagged
tree for 76 of 76 evaluated files.

The harness record files used are those of harness `v1.0.34` (commit
`b4f14af6abd70818f0f98402b215cf322d1392ad`); they are unchanged at the
staging HEAD `39980d5` and are deposited at `v1.0.35`.

Compiler: Microsoft C/C++ Optimizing Compiler 19.51.36260 for x64 (MSVC
14.51.36231), flags `/nologo /std:c++latest /O2 /EHsc /fp:precise
/arch:AVX2`. Python 3.14.7, mpmath 1.4.1.

## Commands

The boys tree is a `git archive` of tag `v5.3.1` at
`C:\Users\myaml\dev\_a-exact-tag-replay`; `BOYS` below is that directory.

The two executables are build products and are not deposited. Rebuild them
with the two `cl` lines below; the digests above are the replay run's. The
byte-identity evidence is the output records, each of which has its deposited
record's digest in `digests.txt`.

Build (`build-tag.log`, exits 0):

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

`run/paths.csv` is 291,827,555 bytes. Its digest is recorded above and it is
byte-identical to the deposited compressed record; the raw file was removed
after verification, and `digests.txt` states this.

## Files

| File | Role |
| --- | --- |
| `README.md` | this record |
| `provenance.json` | machine-readable identity, inputs, runs and results |
| `digests.txt` | sha256 of the binaries, outputs and logs |
| `build-tag.log` | compiler banner and both builds |
| `run-sweeps.log` | the three sweep runs |
| `run/replay-tag.log` | the certificate model replay |
| `run/certify-tag.log` | the full bound validation |
| `run/*.csv` | the reproduced sweep records |
| `run/certificates/*` | the reproduced certificate records |
| `build/sweep1.exe`, `build/sweep2.exe` | not deposited (build products); rebuilt by the two `cl` lines below, SHA-256 `42f45ad1faba2cfa69b8fbbd7ee61af4f4376814f30ad8c6c9c52372498e4d79` and `3b6a92cbb845f31c4af1304257a4552d3285a6a658d9626e4e0f7698f9b8c201` |
| `compare_identity.py` | the evaluated-file identity comparison |
| `run/identity-comparison.log` | its output, including the direct v5.1.1 comparison |
| `residual_probe.py`, `run/residual-probe.log` | the S18.4 residual lemma check |

No paper or deposit file was modified by this record.
