# Reproduction manifest

Harness release `v1.0.35` (`VERSION`). Kernel release boys `v5.3.1`. This file
names, for Papers A, B and C, the commands that reproduce or verify each
deposited record, the environment each record was made in, which records are
reproduced and which are only re-serialised, and what the deposit does not
carry. Commands run from the repository root unless a working directory is
stated. `<boys>` is a checkout of the kernel at the tag named; `<harness>` is
this repository.

The papers cite harness `v1.0.34` today and re-pin to `v1.0.35` at submission.
The records carried from `v1.0.34` are unchanged except two corrected digest
rows in `amplification-sweeps/certified-bound/v531-check/`, recorded there and
in `boys-crossover/exact-tag-replay/`.

## Identity

| Item | Value |
| --- | --- |
| Harness | tag `v1.0.35`, `VERSION` 1.0.35; the tag is the archival identifier and names the release commit |
| Kernel | boys `v5.3.1`, tag object `d1e63f71c6f2501aecc6340e31216a42c2a390b9`, commit `2f83e40a6504368b7b2f13682e29b4f24ab2a40c`, tree `344df3ecc1c13122f59a8a98a78afbff3c34337b` |
| Earlier kernel revisions named by records | v5.0.0 (first and second sweeps, bridge provenance); v5.1.1 (9,353-argument record, seed-spliced, held-out, Chebyshev recount); v5.2.5 and v5.2.6 (certificate; Paper C core and sensitivity sweeps); v5.3.1 (exact-tag replay, accuracy gate, coefficient identity) |

```sh
git clone --branch v5.3.1 https://github.com/myamlak/boys.git boys
git clone --branch v1.0.35 https://github.com/myamlak/boys-harness.git boys-harness
```

Each record's own `provenance.json` names the revision it was measured at. That
name is the record's provenance, not a current reference.

## Environments

| Paper | Compilers and toolchains | Python and numerical stack | CUDA |
| --- | --- | --- | --- |
| A | MSVC 19.51.36260 for x64 (14.51.36231), `/nologo /std:c++latest /O2 /EHsc /fp:precise /arch:AVX2`; GCC 15, C++23, for the CI smoke build | Python 3.14.7, mpmath 1.4.1 | no record is a GPU run |
| B | MSVC 14.5 (Visual Studio 18) for the Boost gate, C++17; gcc 15.2 on WSL Ubuntu for crlibm; C++23 and CMake 3.20 or newer for the composition audit | Python 3.10 or newer; mpmath 1.4.1; gmpy2 2.3.1 linked to MPFR 4.2.2 | the CUDA source was inspected, not run on a named device; no GPU measurement |
| C | none for the counts | Python 3.14.7 (MSC v.1944, 64-bit); mpmath 1.4.1; numpy 2.5.2; scipy 1.18.1; baryrat 2.1.2; gmpy2 2.3.1; `math.fma` present | none |

The Python records run on Windows 11 (10.0.26200) on AMD64. Unless a record's
own README says otherwise (the certificate's commands use `python -I`), run the
scripts without `python -I`: isolated mode hides the user-site `mpmath`.

## Paper A

Records: `amplification-sweeps/` (the two sweeps, certificate, strengthening,
seed-spliced, held-out, surrogate-control, apriori-bound, workload-pair),
`chemistry-bridge/`, and `boys-crossover/exact-tag-replay/` and
`boys-crossover/r4-and-9353/`.

Build and run a sweep from its directory (MSVC form):

```bat
cl /nologo /std:c++latest /O2 /EHsc /fp:precise /arch:AVX2 /I<boys>\include sweep1.cpp /Fe sweep1.exe
python grid.py --out data/grid.csv
python ../common/reference.py --grid data/grid.csv --out data/reference.csv
sweep1.exe
python analyze.py --tree .
```

The exact-tag replay rebuilds both drivers against the `v5.3.1` headers and
compares every produced record with the deposited one. From
`boys-crossover/exact-tag-replay`:

```bat
cl /nologo /std:c++latest /O2 /EHsc /fp:precise /arch:AVX2 /I<boys>\include /Fe build\sweep1.exe <harness>\amplification-sweeps\product-amplification\sweep1.cpp
cl /nologo /std:c++latest /O2 /EHsc /fp:precise /arch:AVX2 /I<boys>\include /Fe build\sweep2.exe <harness>\amplification-sweeps\amplification-form\sweep2.cpp
build\sweep1.exe --grid <harness>\amplification-sweeps\product-amplification\first-sweep-9353\sweeps\product_amplification\data\grid.csv --cells run\cells-9353.csv --admission run\admission-9353.csv --gain run\replica_gain-9353.csv --batch full
build\sweep1.exe --grid <harness>\amplification-sweeps\product-amplification\data\grid.csv --cells run\cells-9352.csv --admission run\admission-9352.csv --gain run\replica_gain-9352.csv --batch full
build\sweep2.exe --grid <harness>\amplification-sweeps\amplification-form\data\grid.csv --fits <harness>\amplification-sweeps\amplification-form\data\fits.csv --paths run\paths.csv --cells run\cells-sweep2.csv --admission run\admission-sweep2.csv
```

The certificate record's commands, from `amplification-sweeps/certified-bound`:

```sh
python -I code/replay.py --kernel <boys> --both-downward
python -I code/certify.py --kernel <boys> --sweep ../amplification-form \
    --model-record inputs/model-record --exp-bits exp/exp_bits.txt --out <run>
```

Compare the produced digests with `boys-crossover/exact-tag-replay/digests.txt`.
The identity check, from `boys-crossover/exact-tag-replay`:

```sh
python compare_identity.py <v5.3.1-tree> <harness>/amplification-sweeps/certified-bound/v531-check/digests.txt <v5.1.1-tree>
```

Expected: measured column 45, v5.3.1 column 76, v5.1.1 tag tree 45. The two
input lists, from `boys-crossover/r4-and-9353`:

```sh
python extract_and_verify.py <boys-v5.3.1-tree> <paper-tree>/papers/boys-crossover <harness>/amplification-sweeps
```

The three paper-tree inputs are the same bytes deposited as
`amplification-sweeps/amplification-form/data/{grid,reference}.csv` and
`amplification-sweeps/product-amplification/first-sweep-9353/sweeps/product_amplification/data/grid.csv`.
The chemistry bridge rebuilds its probe inside the host framework tree (not
deposited); `chemistry-bridge/scripts/commands.md` gives the exact commands, and
`scripts/analyze.py`, `benzene-pvdz/compare_arms.py` and
`scripts/above_b2_max.py` regenerate its summaries from the deposited streams.

Reproduced and re-serialised:

- Reproduced byte for byte at `v5.3.1`: the three sweep records (9,353, the
  designed 9,352, and the second sweep), their admission and
  replica-gain/path records, and the four certificate outputs (exact-tag
  replay).
- Reconstructed and verified against the deposits: the 73-argument R4 list
  (2,409 cells, worst difference `4.975492884e-46` at the deposited cell) and
  the 9,353rd argument.
- Delivered as measured, not re-run: the original sweep records (v5.0.0 and
  v5.1.1 provenance), the strengthening, seed-spliced, held-out,
  surrogate-control, apriori-bound and workload-pair records, and the
  chemistry-bridge water and benzene records. Their analysis tables and
  summaries are regenerated from the deposited records.
- `run/paths.csv` (291,827,555 bytes) is stored as a digest only.

Missing artifacts:

- The pre-release 688: the pre-release run's count. Its selection criterion was
  not preserved, so it records that run, not a re-applicable rule, and cannot
  be re-applied or checked.
- The host-engine geometry: the water geometry is stated (O at the origin; the
  two H at `(±1.430428808474167, 1.107157044080814, 0)` bohr) and the host
  fixture is identified by SHA-256 `ef2262e4d008f8e96bc2b596dc35be8926f8a1d3a449cd09399b729713e2086c`,
  but the host framework tree that built the probe and the fixture's commit are
  not deposited.
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
python tools/check_registration.py --boys <boys>
cd recount && python -I recount.py --boost <boost-1.92.0> --crlibm <crlibm-1.0beta5> \
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
recompute from `records/cells/` and the deposited reference. The accuracy gate
runs as part of the `boys` test invocation at `v5.3.1`.

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

Commands:

```sh
cd coefficient-count-sweep
python count_sweep.py --static --out static_counts.json
python records/chebyshev_recount/recount_cheb.py --out records/chebyshev_recount/recount-run-v5.1.1.json \
    --verify records/chebyshev_recount.json --report records/chebyshev_recount/verify-record-new.json \
    --metadata records/chebyshev_recount/run-metadata-new.json
cd control-60dig && python verify_cells.py && python verify_freeline.py
cd ../strong-optional && python root-isolation/root_isolation.py --selftest \
    && python best-vs-best/best_vs_best.py --selftest && python summarize.py
```

New measurements use the portable `count_sweep.py`; the recorded rows use
`count_sweep_recorded.py`. To run a recorded driver, place
`count_sweep_recorded.py` beside it as `count_sweep.py`; the recorded revision
allows `k <= 12` and the degree modes the portable revision does not.
`stability_save_verify.py --source records/core/sweep_rel_1e-13.json records/core/sweep_abs_5e-14.json`
runs the stability pass against the deposited records. The 60-digit control
README gives the per-cell invocations and the reproduction protocol.

Reproduced and re-serialised:

- Reproduced: the Chebyshev recount, both the corrected run and the
  pre-repair control; the 60-digit control's five cells and the free-line
  effect; the strong-optional selftests and report; the static counts. The four
  core-side records were re-measured at the revisions their own `boys` blocks
  name.
- Re-serialised: the ladder record is the re-serialised original and is not
  re-run; the other historical records are content-equivalent after
  re-serialisation. The `baryrat` 2.1.2 records are copied, not regenerated
  (the solver is not deposited). The 135 original core and ladder rows are not
  re-run by this deposit.

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

The paper-side gates that consume these records are `python tools/check_paper_numbers.py`
in Papers A and C, `python tools/build_ledger.py` in Paper C, and
`python tools/check_source_facts.py --boys <boys>` in Paper B. They live in the
manuscript workspaces and are not part of this deposit.
