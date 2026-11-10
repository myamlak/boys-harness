# Reproduction manifest

Harness release `v1.0.38` (`VERSION`). Kernel release boys `v5.3.2`. This file
names, for Papers A, B and C, the commands that reproduce or verify each
deposited record, the environment each record was made in, which records are
reproduced and which are only re-serialised, and what the deposit does not
carry.

The papers cite harness `v1.0.38` at submission; the tag is cut from the merge
commit.
The recorded values carried from `v1.0.37` are unchanged. This release adds the
full-cube order-three decomposition at `composition-audit/fullcube/`, whose
numbers Paper B reports. The previous release moved
the kernel pin in `.github/workflows/harness.yml` from `v5.3.1` to `v5.3.2`,
a documentation, version-string and test-only change, so that the papers'
submission tags name one kernel. The `v1.0.36` release added the held-out
exact-tag replay at `amplification-sweeps/held-out-v531/` and its TeX fixture
at `tables/prediction-record.tex`.

## Paths

Every command below runs from the root of this repository (the `boys-harness`
clone) unless its block names another working directory. The two clones and
the tag checkouts:

```sh
git clone --branch v5.3.2 https://github.com/myamlak/boys.git boys
git clone --branch v1.0.38 https://github.com/myamlak/boys-harness.git boys-harness  # (cut after this deposit merges)
cd boys-harness
git -C ../boys worktree add ../boys-v5.1.1 v5.1.1
git -C ../boys worktree add ../boys-v5.2.5 v5.2.5
git -C ../boys worktree add ../boys-v5.2.6 v5.2.6
```

The blocks use these names. In `sh` set them with `export`; in `cmd.exe` set
the same names with `set NAME=value`.

| Name | Value | What it is |
| --- | --- | --- |
| `BOYS` | `$PWD/../boys` | the `v5.3.2` tag checkout |
| `BOYS_V511` | `$PWD/../boys-v5.1.1` | the `v5.1.1` tag worktree |
| `BOYS_BUILD` | `$PWD/../boys-build` | a build directory outside the `boys` tree |
| `HARNESS` | `$PWD` | this repository |
| `BOOST_SRC` | path to the Boost.Math 1.92.0 source tree | external, pinned by archive SHA-256 |
| `CRLIBM_SRC` | path to the crlibm 1.0beta5 source tree | external, pinned by archive SHA-256 |
| `RUN` | a fresh directory for a run's outputs, e.g. `$PWD/run/certificates` | never reuse a deposited `results/` directory |

The paper-tree inputs that `extract_and_verify.py` reads are three files, and
all three are deposited; the stand-in is built from the deposit in the Paper A
block below, so no manuscript checkout is needed. Commands that write into the
deposit are a record's own regeneration commands; run them in a copy if the
deposit must stay unchanged.

## Identity

| Item | Value |
| --- | --- |
| Harness | tag `v1.0.38`, `VERSION` 1.0.38. The tag is annotated and names the release commit. A file inside that commit cannot carry the tag object that points to it: resolve the tag object, commit and tree with `git rev-parse v1.0.38`, `git rev-parse 'v1.0.38^{}'` and `git rev-parse 'v1.0.38^{tree}'`. The preceding release, `v1.0.37`, is tag object `a07d808ea12ca4394104bf4054eabd2407ff1986`, commit `74eecb59de7e8d408872083c06a33284759ef6de`, tree `f832c317e1fdeb9076a439e5c80cf0609f1fe20b`. |
| Kernel | boys `v5.3.2`, tag object `058f10f273f37ca3b5769e9f1173c218f38fce6c`, commit `6739447d82765ae5f7f68351a99f418cd018dd9a`, tree `71be8196d53c31c8ec5d9be2c9fb48453b297c4e` |
| Earlier kernel revisions named by records | v5.0.0 (first and second sweeps, bridge provenance); v5.1.1 (9,353-argument record, seed-spliced, held-out, Chebyshev recount); v5.2.5 and v5.2.6 (certificate; Paper C core and sensitivity sweeps); v5.3.1 (exact-tag replay, accuracy gate, coefficient identity, and the release the earlier paper tags named) |

Each record's own `provenance.json` names the revision it was measured at. That
name is the record's provenance, not a current reference.

Recorded SHA-256 values are of the bytes the measurement machine read. Several
deposited CSVs were committed with LF line endings while the measurement tree
carried CRLF, so a fresh clone's copy of one of those files can hash
differently while parsing identically. The checks that compare values parse
the files; a digest comparison uses the record's own `digests.txt` list, which
is the measurement bytes. The R4 record's README names the three input files
this affects.

## Environments

| Paper | Compilers and toolchains | Python and numerical stack | CUDA |
| --- | --- | --- | --- |
| A | MSVC 19.51.36260 for x64 (14.51.36231), `/nologo /std:c++latest /O2 /EHsc /fp:precise /arch:AVX2`; GCC 15, C++23, for the CI smoke build | Python 3.14.7, mpmath 1.4.1 | no record is a GPU run |
| B | MSVC 14.5 (Visual Studio 18) for the Boost gate, C++17; gcc 15.2 on WSL Ubuntu for crlibm; C++23 and CMake 3.20 or newer for the composition audit | Python 3.10 or newer; mpmath 1.4.1; gmpy2 2.3.1 linked to MPFR 4.2.2 | the CUDA source was inspected, not run on a named device; no GPU measurement |
| C | none for the counts | Python 3.14.7 (MSC v.1944, 64-bit); mpmath 1.4.1; numpy 2.5.2; scipy 1.18.1; baryrat 2.1.2; gmpy2 2.3.1; `math.fma` present | none |

The records were made with the packages above. Install the Python stack into
the interpreter the commands use:

```sh
py -3 -m pip install mpmath==1.4.1 numpy==2.5.2 scipy==1.18.1 baryrat==2.1.2
```

`gmpy2` 2.3.1 is installed against MPFR 4.2.2. The Python records run on
Windows 11 (10.0.26200) on AMD64.

A command that uses `python -I` ignores the user site: its packages must be in
the interpreter's own site-packages, so install without `--user`, or run the
script without `-I`. The certificate commands and the Paper B recount use
`python -I`; every other command here runs without it.

## Paper A

Records: `amplification-sweeps/` (the two sweeps, certificate, strengthening,
seed-spliced, held-out, held-out-v531, surrogate-control, apriori-bound,
workload-pair), `chemistry-bridge/`, and `boys-crossover/exact-tag-replay/` and
`boys-crossover/r4-and-9353/`.

Build and run a sweep from its directory (MSVC form): `product-amplification/`
holds `sweep1.cpp` and `amplification-form/` holds `sweep2.cpp`.

```bat
cl /nologo /std:c++latest /O2 /EHsc /fp:precise /arch:AVX2 /I"%BOYS%\include" sweep1.cpp /Fesweep1.exe
python grid.py --out data/grid.csv
python ../common/reference.py --grid data/grid.csv --out data/reference.csv
sweep1.exe
python analyze.py --tree .
```

Four amplification records are stored xz-compressed. From
`amplification-sweeps/`, decompress them before running the analysis:

```sh
xz -dk amplification-form/results/paths.csv.xz \
      product-amplification/data/reference.csv.xz \
      product-amplification/results/cells.csv.xz \
      product-amplification/results/replica_gain.csv.xz
```

The exact-tag replay rebuilds both drivers against the `v5.3.1` headers and
compares every produced record with the deposited one. From
`boys-crossover/exact-tag-replay`, create the two output directories first
(the repository tracks no files beneath them):

```bat
mkdir build run
cl /nologo /std:c++latest /O2 /EHsc /fp:precise /arch:AVX2 /I"%BOYS%\include" /Febuild\sweep1.exe "%HARNESS%\amplification-sweeps\product-amplification\sweep1.cpp"
cl /nologo /std:c++latest /O2 /EHsc /fp:precise /arch:AVX2 /I"%BOYS%\include" /Febuild\sweep2.exe "%HARNESS%\amplification-sweeps\amplification-form\sweep2.cpp"
build\sweep1.exe --grid "%HARNESS%\amplification-sweeps\product-amplification\first-sweep-9353\sweeps\product_amplification\data\grid.csv" --cells run\cells-9353.csv --admission run\admission-9353.csv --gain run\replica_gain-9353.csv --batch full
build\sweep1.exe --grid "%HARNESS%\amplification-sweeps\product-amplification\data\grid.csv" --cells run\cells-9352.csv --admission run\admission-9352.csv --gain run\replica_gain-9352.csv --batch full
build\sweep2.exe --grid "%HARNESS%\amplification-sweeps\amplification-form\data\grid.csv" --fits "%HARNESS%\amplification-sweeps\amplification-form\data\fits.csv" --paths run\paths.csv --cells run\cells-sweep2.csv --admission run\admission-sweep2.csv
```

The certificate record's commands, from `amplification-sweeps/certified-bound`.
`RUN` must be a fresh directory:

```sh
mkdir -p "$RUN"
python -I code/replay.py --kernel "$BOYS" --both-downward
python -I code/certify.py --kernel "$BOYS" --sweep ../amplification-form \
    --model-record inputs/model-record --exp-bits exp/exp_bits.txt --out "$RUN"
```

Compare the produced digests with `boys-crossover/exact-tag-replay/digests.txt`.
The identity check, from `boys-crossover/exact-tag-replay`:

```sh
python compare_identity.py "$BOYS" "$HARNESS/amplification-sweeps/certified-bound/v531-check/digests.txt" "$BOYS_V511"
```

Expected: measured column 45, v5.3.1 column 76, v5.1.1 tag tree 45. The two
input lists, from `boys-crossover/r4-and-9353`: the script's second argument is
the paper tree, and the three files it reads are deposited, so build the
stand-in first.

```sh
mkdir -p paper-inputs/sweeps/amplification_form/data paper-inputs/sweeps/product_amplification/data
cp "$HARNESS/amplification-sweeps/amplification-form/data/grid.csv" \
   "$HARNESS/amplification-sweeps/amplification-form/data/reference.csv" \
   paper-inputs/sweeps/amplification_form/data/
cp "$HARNESS/amplification-sweeps/product-amplification/first-sweep-9353/sweeps/product_amplification/data/grid.csv" \
   paper-inputs/sweeps/product_amplification/data/grid.csv
python extract_and_verify.py "$BOYS" paper-inputs "$HARNESS/amplification-sweeps"
```

The three paper-tree inputs are the same bytes deposited as
`amplification-sweeps/amplification-form/data/{grid,reference}.csv` and
`amplification-sweeps/product-amplification/first-sweep-9353/sweeps/product_amplification/data/grid.csv`.
The script writes its outputs beside itself; run it in a copy if the deposit
must stay unchanged.

The chemistry bridge rebuilds its probe inside the host framework tree, which
is not deposited. The record names the host tree `bladawiec` at commit
`ac8540b7c23d3818c9e67f8d59bf6cfbbef7fa0d` or
later, and pins the water fixture by SHA-256
`ef2262e4d008f8e96bc2b596dc35be8926f8a1d3a449cd09399b729713e2086c`.
`chemistry-bridge/scripts/commands.md` gives the exact commands, including the
cc-pVTZ arms (`qcx-bridge-probe.exe <out>/ccpvtz-default data/basis/cc-pvtz`,
with `<out>` a fresh output directory) that regenerate the undeposited dense
tensors and argument streams; `benzene-pvdz/commands.md` gives the benzene arm
(`qcx-bridge-probe.exe <out-prefix> data/basis/cc-pvdz benzene`).
`scripts/analyze.py`, `benzene-pvdz/compare_arms.py` and
`scripts/above_b2_max.py` regenerate its summaries from the deposited streams.

Reproduced and re-serialised:

- Reproduced byte for byte at `v5.3.1`: the three sweep records (9,353, the
  designed 9,352, and the second sweep), their admission and
  replica-gain records, and the four certificate outputs (exact-tag
  replay). The second sweep's path record is compared by digest only: it is
  291,827,555 bytes and is stored as a digest, not as a file.
- Reconstructed and verified against the deposits: the 73-argument R4 list
  (2,409 cells, worst difference `4.975492884e-46` at the deposited cell) and
  the 9,353rd argument.
- Delivered as measured, not re-run: the original sweep records (v5.0.0 and
  v5.1.1 provenance), the strengthening, seed-spliced, held-out,
  surrogate-control, apriori-bound and workload-pair records, and the
  chemistry-bridge water and benzene records. Their analysis tables and
  summaries are regenerated from the deposited records.
- `run/paths.csv` (291,827,555 bytes) is stored as a digest only, so
  byte-for-byte reproduction cannot be a comparison against a deposited file:
  compare its SHA-256 with `boys-crossover/exact-tag-replay/digests.txt`.

Missing artifacts:

- The pre-release 688: the pre-release run's count. Its selection criterion was
  not preserved, so it records that run, not a re-applicable rule, and cannot
  be re-applied or checked.
- The host-engine geometry: the water geometry is stated (O at the origin; the
  two H at `(±1.430428808474167, 1.107157044080814, 0)` bohr) and the host
  fixture is identified by SHA-256
  `ef2262e4d008f8e96bc2b596dc35be8926f8a1d3a449cd09399b729713e2086c`,
  but the host framework tree that built the probe is not deposited.
- The cc-pVTZ and benzene dense tensors and argument streams: not deposited;
  the probe regenerates them (about twenty minutes per tensor).

Finite-grid versus continuous: all of Paper A's checks are finite-grid. The
certificate validates 55,176 entry-cells and 132 band comparisons; it is
conditional and makes no continuum claim. The residual lemma is derived, and
the probe checks its premises on the two grids. The replay is byte identity on
the finite records.

## Paper B

Records: `accuracy-claim-audit/` (the coverage enumeration, the captured cell
summary, the Boost.Math and crlibm audits and the taxonomy testbed),
`composition-audit/`, and `chemistry-bridge/` where the program paper cites the
bridge. The accuracy gate's records and reference grids are in the `boys`
repository at `v5.3.1`.

Commands:

```sh
cd accuracy-claim-audit
python headline_numbers.py
python tools/check_registration.py --boys "$BOYS"
cd recount && python -I recount.py --boost "$BOOST_SRC" --crlibm "$CRLIBM_SRC" \
    --registry registry/counts.json --out counts.json
```

The Boost.Math and crlibm audits rebuild from the pinned external sources; the
commands, build stanzas and expected outputs are in `accuracy-claim-audit/README.md`
and `claim-audit-crlibm/work/`. The composition audit, from `composition-audit`:

```bat
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --target measure-interaction
powershell -File run_measurement.ps1 -QuietMachine
```

On Linux or WSL use `bash run_measurement.sh --quiet-machine --reference results/reference.csv --arguments-reference results/arguments.csv`.
`signed_interaction.py`, `synthetic_surface.py` and `figure-ratio.json`
recompute from `records/cells/` and the deposited reference.

The accuracy gate replay is separate from `check_registration.py`. From the
`boys` checkout, configure, build and run the gate, then hold the committed
run to the tree by content:

```sh
cmake -S "$BOYS" -B "$BOYS_BUILD" -DCMAKE_BUILD_TYPE=Release \
    -DBOYS_BUILD_TESTS=ON -DBOYS_BUILD_C_API=ON -DBOYS_BUILD_PROBE=ON -DBOYS_FP16_SEAM=ON
cmake --build "$BOYS_BUILD" --target boys-accuracy-gate --config Release
ctest --test-dir "$BOYS_BUILD" -C Release -R '^boys-accuracy-gate$' --output-on-failure
python "$BOYS/tools/check_recorded_run.py" --run "$BOYS/tests/data/boys_accuracy_gate_run.txt"
```

The four options are the top-level defaults; naming them makes the build the
one the record describes. `--config Release` and `-C Release` name the
configuration for a multi-config generator; a single-config generator ignores
them after `-DCMAKE_BUILD_TYPE=Release`. The gate runs against the committed
reference grids `tests/data/boys_accuracy_gate_reference.csv` (25 significant
digits) and `tests/data/boys_reference.csv` (45 digits), both in the `boys`
tree.
`check_recorded_run.py` holds `tests/data/boys_accuracy_gate_run.txt` to the
tree's version and input digests; it does not run the gate. The committed run
is a CUDA build's, and a build without the device lanes certifies the host
lanes alone: the gate prints its device-lane count and the shortfall.

Reproduced and re-serialised:

- Reproduced: the claim census and the counts from the deposited tools; the
  Boost.Math and crlibm audits from the pinned sources; the composition
  audit's interaction records from the compressed cell files; the gate from the
  `boys` tree.
- Delivered records, checked rather than regenerated by one command:
  `accuracy-claim-audit/records/documented-claims-v531.tsv` (checked by
  `check_registration.py`), `records/summarize_cells.txt` (read by
  `headline_numbers.py`), and the earlier as-found record and erratum, which
  are provenance.

Missing artifacts:

- The un-deposited kernel grid: the gate's 25-significant-digit and 45-digit
  reference grids, and the generators that reproduce them, are committed in the
  `boys` repository at `v5.3.1`, not in this harness deposit. The 45-digit grid
  is also the R4 source used by Paper A and Paper C.
- The external audit sources: Boost.Math 1.92.0 (archive SHA-256
  `c4a3b310ddd2472416e091067166b0713be97c63f38c212c484ada022fd296ce`)
  and crlibm 1.0beta5 (archive SHA-256
  `79d626f9515fccdbc724e8b74cbf848331317b1891ef22aaaf6a195e20bbbfe5`)
  are pinned but not deposited.
- The program demonstration record is in the submission package, not the
  harness.

Finite-grid versus continuous: the gate's verdicts are conformance at committed
grid points and make no statement between them. The audit's claims are measured
on the pinned grids. No continuous statement is made.

## Paper C

Records: `coefficient-count-sweep/` (the core, ladder, stability,
re-certification, second-edge, second-instance, F1, third-edge, shipped-seed,
Chebyshev recount, 60-digit control and strong-optional records).

Commands, from `coefficient-count-sweep`:

```sh
python count_sweep.py --static --out "$RUN/static_counts.json"
python records/chebyshev_recount/recount_cheb.py --out records/chebyshev_recount/recount-run-v5.1.1.json \
    --verify records/chebyshev_recount.json --report records/chebyshev_recount/verify-record-new.json \
    --metadata records/chebyshev_recount/run-metadata-new.json
cd control-60dig && python verify_cells.py && python verify_freeline.py
cd ../strong-optional && python root-isolation/root_isolation.py --selftest \
    && python best-vs-best/best_vs_best.py --selftest && python summarize.py
```

New measurements use the portable `count_sweep.py`; the recorded rows use
`count_sweep_recorded.py`. The recorded rows are 135: 40 in
`records/core/sweep_all.json`, 20 in `records/core/sweep_rel_1e-13.json`, 20 in
`records/core/sweep_abs_5e-14.json`, 15 in `records/core/sweep_rel_1e-14.json`
and 40 in `records/ladder/sweep_all.json`. The four core-side files were
re-measured at `v5.2.5`; the ladder file was not re-run.

To run a recorded driver, place `count_sweep_recorded.py` beside it as
`count_sweep.py`; the recorded revision allows `k <= 12` and the degree modes
the portable revision does not. The drivers differ in how they take the kernel
root: `count_sweep_recorded.py` and `recertify_ladder.py` take `--boys-root`;
`rational_unrounded_control.py` takes `--boys`; `cold_ladder.py` and
`stability_save_verify.py` use `count_sweep.DEFAULT_BOYS` and have no override;
`rational_control_symmetric.py` uses its own `RAT_BOYS` and `CHEB_BOYS`
constants. Several read the measurement tree's `results*/` layout; the deposit
supplies the records and the `--source` options, and `coefficient-count-sweep/README.md`
("Harness revisions") maps each driver.

The stability pass runs against the deposited records:

```sh
python stability_save_verify.py --source records/core/sweep_rel_1e-13.json records/core/sweep_abs_5e-14.json
```

The 60-digit control README gives the per-cell invocations and the reproduction
protocol. The Chebyshev recount's identities are in
`records/chebyshev_recount/README.md`: the corrected generator
`gen_boys_coefficients_v5.1.1.py` (SHA-256
`11185fba03abe74d01b30f093818f48972eb8d88413d446b505198b869115180`), the
pre-repair control `gen_boys_coefficients_v5.0.0.py` (SHA-256
`ccc8b011f4170ce0b530bc4781b09670f42839fa18caaa38ea7d3a9ea9ba9bb4`), both
deposited under `records/chebyshev_recount/generators/`, and the coefficient
header the recount read (SHA-256
`8ce5e17157a4be5f621d695f377094d212cb839732e278d47f23aa225b4d9204`). The four
core-side records were re-measured at `v5.2.5`, whose header is
`d978d88bcfbb60727cab3eed32039c9d66ce917a14c7f5a199b6215db8f7ea00` and whose
generator (not deposited) is
`1d8b88ac30c757138a4c87219d68b56756c175c36230818b2360dc0bf7bdd373`; each
record's `meta` block names both.

Reproduced and re-serialised:

- Reproduced: the Chebyshev recount, both the corrected run and the
  pre-repair control; the 60-digit control's five cells and the free-line
  effect; the strong-optional selftests and report; the static counts. The four
  core-side records were re-measured at the revisions their own `boys` blocks
  name.
- Re-serialised: the ladder record is the re-serialised original and is not
  re-run; the other historical records are content-equivalent after
  re-serialisation. The re-serialisation is enumerated file by file in
  `records/chebyshev_recount/paper-record-digests.tsv`, with the comparison
  method in the `.md` beside it: the measurement-tree originals use CRLF, the
  deposit stores LF, and rewriting each LF as CRLF reproduces the original
  digest exactly. The four core-side files were later re-measured at `v5.2.5`
  and are not content-equivalent to the originals; the ladder file is still the
  re-serialised original. The `baryrat` 2.1.2 records are copied, not
  regenerated (the solver is not deposited). The 135 original core and ladder
  rows are not re-run by this deposit.

Missing artifacts:

- The historical 24-coefficient ladder cell: at width `2^5` and target
  `1e-15`, the corrected record stores 25 `[22/2]` against the historical 24
  `[21/2]`. The historical record carries no fitted coefficients, so its 24
  cannot be re-evaluated at 60 digits; the paper records it as an unresolved
  historical comparison.
- The `baryrat` 2.1.2 solver: the independent-solver records are copied
  records, not regenerated.
- The measurement tree's `results*/` layout: some recorded drivers hardcode it;
  the deposit provides the records and the `--source` options.

Finite-grid versus continuous: the counts are solver-relative existence bounds
on the 2,001-node certifying grid, with the 4,001-node check grid for the
stability pass; a dense maximum is not a continuum certificate. The
root-isolation record is the continuous-side check: exact Sturm isolation and
interval-wide denominator bounds for all 98 accepted rational fits in the
corrected x0, second-edge a2 and 60-digit control records. Continuous-error
certification beyond that is not claimed.

## Consumers

The paper-side gates that consume these records are manuscript tools, not part
of this deposit. Their paths and SHA-256 as of this release:

| Paper | Tool | SHA-256 |
| --- | --- | --- |
| A | `tools/check_paper_numbers.py` | `25867469f6eb41e8cf824ec9f1d0ba6b466dafd892b162f764bb75bbb0079c39` |
| B | `tools/check_source_facts.py --boys "$BOYS"` | `2cd153750ac5ad9340ee340cb9eae6ba4398453187ddd6be98f9ca51849b7a7b` |
| C | `tools/check_paper_numbers.py` | `1e7bfe852e16b303921eb28b8dbb1b29d1efa1a672a3df886927308e07268e73` |
| C | `tools/build_ledger.py` | `1b26b3a69abd62d955bb764b32aebaea226ed9595f7432df644f55388ac79644` |

They read each paper's own tree and its `tools/ledger.json`; an outsider
cannot run them without that workspace. The deposited records they check are
named above. The digest names the revision this manifest was written against,
not a file in this deposit.
