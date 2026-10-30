# The v5.2.0 source and execution check

The check described here was run at boys v5.2.0. The current release is
v5.2.1; the source-identity confirmation against it is in
`v5.2.1-note.md`. The identities and results below are the v5.2.0 check as
run.

The annotated tag object
`cb1ce4c79eb22a83bc2fa116c6647df9b08a0e2e` names `v5.2.0` and carries commit
`55ca0e316ad2f9f776424bfee9f365dfc322058b` (short `55ca0e3`), whose tree is
`312e72ab11f314d41834286143ce5549bafad72c`; the tag is dated
2026-10-28T18:09:47Z and its message is "boys 5.2.0: the claim-book
registration table, the measured rows and the exhaustive half sweeps". Its
parent is the historical v5.1.1 commit
`834bdbb18a30ac8538b519436429f78d41f187ce`. The tag, commit and tree were read
from the repository's loose ref and object files; git was not invoked.

Two questions are answered here, and they are separate. First, the source
identity: the evaluated code at v5.2.0 against the v5.1.1 code the paper's
numbers were measured with. Second, the execution check: both sweeps and the
integral-engine bridge were actually compiled and run at v5.2.0, on the
deposited grids, fits and 60-digit references, and every record they produced
was compared with the deposited one. A difference is reported as a difference
throughout; none is smoothed.

## Source identity

The execution check compiled against `C:\Users\myaml\dev\_boys-v520-pub`,
which is clean at the tag: HEAD is the tag's commit, all 710 tracked files
match the tag tree's blobs, 0 differing, 0 missing, 0 extra — "CLEAN AT THE
TAG" (`run/pub-tree-identity.txt`). The evaluation comparison additionally
names the `myamlak/accuracy-gate-umbrella` worktree at
`C:\Users\myaml\dev\_boys-umbrella`, head commit
`14c25bf8df04eaf585790aa2674bafcb15d65b8c`, which has the same tree; its 76
evaluated files match the tagged tree's blob SHA-1 one for one
(`evaluation-sources-compare.txt`).

Compared with the v5.1.1 tree:

- `include/boys`: 62 files, 61 identical and one differing. The difference is
  `version.hpp`, and it is only the version constant: minor 1 to 2, patch 1 to
  0, the string `5.1.1` to `5.2.0`. No other line differs
  (`version-hpp.diff`).
- `src`: 14 files, all identical.
- The v5.1.1 provenance record's 15 evaluated-file digests agree on 14; the
  fifteenth is `version.hpp`.

The changes outside the evaluated path are the gate and probe work the release
is for: `tests/boys_accuracy_gate.cpp`, `tests/boys_gate_reference.hpp`, the
two gate data files, the new `tests/data/half-sweep/` tables, `CMakeLists.txt`,
the repository README and docs, the changelog, the citation and
`tools/check_probe_table.py`. None is compiled by the two sweeps or the
integral-engine bridge. The version constant is read only by `VersionString()`
and by the probe dump's protocol line, not by any evaluation arithmetic.

Every measured figure in the paper is therefore the historical v5.1.1 run's
figure, and v5.2.0 compiles the same evaluation code. The execution check
confirms that rather than assuming it.

## Execution check

The check ran in a scratch tree under `run/`. Its two sweep trees are copies of
the paper tree's `sweeps/product_amplification` and
`sweeps/amplification_form`; every input and harness file the run reads — the
grids, the fits, the 60-digit references, the validation files, `SWEEP.md`,
`analyze.py`, `grid.py`, `fits.py`, `run.ps1` and the two `.cpp` sources — is
byte-identical to the paper tree's copy: 20 compared, 20 identical, 0
differing, 0 missing, 0 extra (`scratch-tree-compare-v520.txt`). The grids,
fits and references are the deposited ones; the check re-runs the evaluation,
not the reference generation.

The compiler was the machine's MSVC 14.51.36231:
`C:\Program Files\Microsoft Visual Studio\18\Community\VC\Tools\MSVC\14.51.36231\bin\HostX64\x64\cl.exe`,
banner `Microsoft (R) C/C++ Optimizing Compiler Version 19.51.36260 for x64`
(`run/logs/compiler-path.log`, `run/logs/compiler-version.log`). No other
compiler or flag set was used.

### The commands

`run/run_execution_check.ps1` holds every command verbatim. With
`S1 = run\sweeps\product_amplification`, `S2 = run\sweeps\amplification_form`
and `B = C:\Users\myaml\dev\_boys-v520-pub`, the compile line is the runbook's
own (`SWEEP.md`), with the include path at the v5.2.0 tree:

```
cl /nologo /std:c++latest /O2 /EHsc /fp:precise /arch:AVX2 /I B\include
   /Fobuild\sweep1.obj /Febuild\sweep1.exe sweep1.cpp
```

and the same line with `sweep2` for the second tree. The run and analysis
commands:

```
S1\build\sweep1.exe --grid S1\data\grid.csv --cells S1\results\cells.csv
    --admission S1\results\admission.csv --gain S1\results\replica_gain.csv
    --batch full
python analyze.py --tree S1 --boys-tree B --provenance S1\results\provenance.json
    --validation S1\data\reference_validation.txt

S2\build\sweep2.exe --grid S2\data\grid.csv --fits S2\data\fits.csv
    --paths S2\results\paths.csv --cells S2\results\cells.csv
    --admission S2\results\admission.csv
python analyze.py --tree S2 --boys-tree B --provenance S2\results\provenance.json
    --validation S2\data\reference_validation.txt
```

and the five table scripts, run from `run/` with `BOYS_SWEEP_TREE=S1`: `python
make_tables.py`, `merge_per_order.py`, `make_endtoend.py`, `make_frontier.py`,
`make_prediction_record.py`. All five are byte-identical to the paper tree's
own copies.

The runbooks' `run.ps1` was not used: its guard calls `git describe` and `git
status` and refuses to start. Its steps were run directly instead, and the
provenance recorder's git-reading `revision()` was replaced by an
environment-supplied identity. The change is 32 diff lines in
`boysprovenance.py` (`recorder-diff-v520.txt`); `reference.py` is identical.
The identity values come from `run/verify_pub_tree.py` and are marked in each
record as `git_invoked: false`, `verified_by:
sweeps/v5.2.0-check/run/verify_pub_tree.py`.

### The steps

`run/logs/execution-check.log`, verbatim (decimal commas are the locale's):

```
v5.2.0 execution check, runbook commands
sweep1-provenance    EXIT=0 SECONDS=0,9
sweep1-build         EXIT=0 SECONDS=23,3
sweep1-run           EXIT=0 SECONDS=1,6
sweep1-analyze       EXIT=0 SECONDS=42,9
sweep2-provenance    EXIT=0 SECONDS=0,9
sweep2-build         EXIT=0 SECONDS=28
sweep2-run           EXIT=0 SECONDS=22,3
sweep2-analyze       EXIT=0 SECONDS=75,6
tables               EXIT=0 SECONDS=3,3
merge-per-order      EXIT=0 SECONDS=0,3
endtoend             EXIT=0 SECONDS=0,6
frontier             EXIT=0 SECONDS=0,9
prediction-record    EXIT=0 SECONDS=0,3
```

## What the check reproduced

### Sweep records

22 sweep-tree files were compared (`sweep-results-compare-v520.txt`): 18
byte-identical after folding line endings, 4 differing.

The 18 identical files: product `results/cells.csv`, `admission.csv`,
`replica_gain.csv`, `per_order.csv`, `band_summary.csv`, `gain_up.csv`,
`growth.csv`, `decomposition.csv` and `data/documented.csv`; form
`results/paths.csv`, `cells.csv`, `admission.csv`, `per_order.csv`,
`band_summary.csv`, `prediction.csv`, `prediction_scan.csv`,
`prediction_record.csv` and `transfer_cells.csv`.

The 4 differing files, with the differing lines:

- product `results/summary.txt`, line 2 only: the printed reference path
  (`sweeps\product_amplification\data\reference.csv` against the scratch
  tree's copy of the same file).
- product `results/provenance.json`: the boys identity fields (`commit`,
  `describe`, `tag_at_head`, `tree`), `git_invoked: false` in place of the
  git-derived status line, and two lines longer (653 against 655).
- form `results/summary.txt`, line 2 only: the printed reference path.
- form `results/provenance.json`: the same identity fields, two lines longer
  (657 against 659).

No computed value differs in any of the four. The differing lines are paths
and the run's own identity.

### Tables and figure

All 14 generated files — the 13 `tables/*.tex` and `figures/per-order-ratio.dat`
— are identical after folding CRLF to LF (`tables-compare-v520.txt`). The
deposited files are CRLF, so the record prints the folded size; the raw sizes
agree too (for example `per-order.tex` is 13199 bytes raw and 13062 folded on
both sides, `prediction-record.tex` 16853 and 16719, `per-order-ratio.dat` 1219
and 1185).

The five table scripts were also run against the deposited records in the
paper tree. Their logs are `run/logs/deposited-*.log`; the logs name exactly
the 14 outputs, and the raw SHA-256 of all 14 files is identical before and
after (`tables-before-deposited.sha256`, `tables-after-deposited.sha256`).

### The new provenance records

The two new records are `run/sweeps/*/results/provenance.json`. Against the
deposited v5.1.1-check copies, product differs in 17 fields and form in 15
(`provenance-diff-v520.txt`). Every difference is an identity or stamp field:

- the six boys identity fields (`commit`, `describe`, `git_invoked`,
  `tag_at_head`, `tree`, `verified_by`);
- the `include/boys/version.hpp` digest (the version constant);
- the harness digests of the scripts that ran — product `SWEEP.md` 11150/`a1fbbf2b`,
  `analyze.py` 18225/`7d5caf0a`, `grid.py` 6195/`b44134ba`, `run.ps1`
  3878/`b7e09e02`; form `SWEEP.md` 10205/`b872096c`, `analyze.py`
  43991/`7da9291c`, `run.ps1` 3601/`57100192`;
- the `flags` string (the v5.2.0 include path and the `/Fo`/`/Fe` form);
- `generated_utc`.

The harness digests in the new records are exactly those of the paper tree's
own scripts, which the check read. That is not true of the deposited
v5.1.1-check copies; see the deposit section below.

### The integral-engine bridge

`compiled-tree-compare-v520.txt` establishes what the bridge binaries were
built from. The bridge's `external/boys` copy has 62/62 `include/boys` files
equal to the v5.2.0 worktree and 61/62 equal to v5.1.1 (the exception is
`version.hpp`); all 14 `src` files are equal to both. Five translation units
are compiled on the CPU lane: `boys.cpp`, `boys_c.cpp`, `boys_orders_simd.cpp`,
`boys_probe.cpp`, `boys_simd.cpp` (the CUDA lane is off). The two probe
binaries:

- `build-bridge-default\bridge\Release\qcx-bridge-probe.exe`: 3668992 bytes,
  sha256 `4d99953feeded6a1251a880bef402adcfa4b40160a8f0f35b1ff3d05212049ec`,
  mtime 2026-10-28T18:50:40Z;
- `build-bridge-paper-a\bridge\Release\qcx-bridge-probe.exe`: 3688448 bytes,
  sha256 `6c75c28de012c42624384eaae6e320a4fe0d779571d00906ef32ce6ed5a63ae0`,
  mtime 2026-10-28T19:03:54Z.

Both are newer than the newest `external/boys` file (`tools_tc/shipped_grid.txt`,
2026-10-28T18:15:18Z). All six probe records state version 5.2.0.

The nine bridge steps all finished EXIT=0 (`bridge-run-summary-v520.txt`
quotes each FINISH line and the resolved command lines): six probe runs at
0.3, 0.2, 0.9, 1.2, 54.8 and 51.1 s and three analyses at 2.4, 11.5 and 441.2
s.

`bridge-compare-v520.txt` compares every record against the v5.1.1 run:

- `.args.txt`, `.matrices.txt`, `.eri_max.txt`, `.default.density.bin` and
  `.eri.txt` are byte-identical in all six arms. The cc-pVTZ dense tensors are
  identical: 270819892 bytes, sha256 `fb36b7c4…` and 270813606 bytes, sha256
  `50d7c1c9…`.
- `.summary.txt` is identical in both STO-3G arms. In the other four the only
  changed key is `basis_directory` (the recorded basis path).
- The twelve `.trace.txt` comparisons differ in one field only, `wall`.
- Each `.boys_record.txt` has one differing line of 102: the version line,
  5.1.1 to 5.2.0.

The stream inventory recomputed from the arguments, identical old and new:

| basis | calls | distinct | nmax | largest x |
| --- | --- | --- | --- | --- |
| STO-3G | 10935 | 2764 | 4 | 28.033987699960576 |
| 6-31G* | 51570 | 14163 | 8 | 153.30510853866076 |
| cc-pVTZ | 262484 | 81627 | 12 | 277.20922853176552 |

For STO-3G the largest deviation is |dF| = 3.219647e-15 at order 0, x =
0.058133111316822465 (3.282% of the contract); for 6-31G* 3.219647e-15 at
order 0, x = 1.92508455331846e-08 (3.224%); for cc-pVTZ 3.219647e-15 at order
0, x = 8.069323580776033e-08, and the largest contract fraction is 5.025% at
order 8, x = 2.0505504716548182 (dF = 2.976785e-15). In B2 the largest
deviation is 1.110223e-16 in every basis. For cc-pVTZ order 12 there is one
cell and its deviation is 0.

The old side is the v5.1.1 run's own output. `bridge-old-records-v520.txt`
records its provenance: 52 files, 52 byte-identical, 0 differing, 0 missing;
the two cc-pVTZ tensors not copied are compared by SHA-256 and are identical;
the old records state 5.1.1. Against the harness working tree
(`papers-harness-repo\chemistry-bridge\outputs`) 48 records match by their
`-v511` names and 4 through the plain file, 0 differing, 0 missing. The
v1.0.11 tag carries a different set: its `chemistry-bridge/outputs` is the
v5.0.0 run, whose records state 5.0.0. Against it, 30 files are the same and
22 differ — 6 boys_record version lines, 12 traces at `wall`, and 4 summaries
on census lines (`basis_source=sto-3g-fixture` against `n=7`; `molecule_atoms=3`
against `shells=10`) — with 0 missing. That difference is between library
versions and is reported as such; it is not the old side of this comparison.

## The v1.0.11 deposit

The deposit is harness tag v1.0.11, tag object
`1aa52650755d7445134a105d289029300096e7dd`, commit
`e4513b62b381d49ee5ec1165f25195ea6bf8ff7f`, tree
`f5596acfddb16b572a75f8dadfe915a1720a0646`, 489 files, VERSION 1.0.11
(`harness-deposit-compare-v111.txt`). The earlier tags are v1.0.9 (436 files,
VERSION 1.0.1) and v1.0.10 (488 files, VERSION 1.0.1); v1.0.10 to v1.0.11 adds
`.gitattributes` and changes 41 files, listed in the record.

The bundle carried by the v1.0.11 tag against this paper tree: 46 entries
mapped, 45 identical, 0 differing, 1 missing — the paper tree does not carry
`sweeps/v5.1.1-check/check_v500_compare.py`. Four bundle entries have no
paper-tree path and are not compared: `README.md`, `make_endtoend.py`,
`make_frontier.py`, `make_tables.py` (their paper tree equivalents live in the
check's `run/`).

The tag carries no `v5.2.0-check` directory and no `SWEEP.md` anywhere. Its
only `-v511` names are the two v5.1.1-check comparison texts; the v5.1.1 run's
own outputs are not deposited. This check is not deposited there either.

### The harness digests in the deposited provenance

`deposit-harness-v520.txt` puts the deposited provenance's harness digests
beside the files they name. They are stored as `digests["harness/<name>"]`;
the top-level `harness` field names the sweep script, not a digest group. The
deposit carries two provenance records per tree, and they do not agree with
each other or, for the form tree, with the paper tree:

- Product: the run record (`results/provenance.json`, generated
  2026-10-27T21:38:40Z) records exactly the digests of the copies the deposit
  carries — `analyze.py` 18225/`7d5caf0a`, `grid.py` 6195/`b44134ba`,
  `run.ps1` 3878/`b7e09e02`, `sweep1.cpp` 8569/`309bcb71` — and the paper
  tree's `SWEEP.md` digest, 11150/`a1fbbf2b`. The v5.1.1-check copy
  (`v5.1.1-check/product-amplification-provenance.json`, generated
  2026-10-26T22:23:19Z) is an earlier vintage: it differs in the four script
  digests (SWEEP.md 10133/`4af424e9`, analyze.py 18145/`a5d94476`, grid.py
  4019/`e7f94908`, run.ps1 3827/`b5e41639`) and `generated_utc`; the
  `sweep1.cpp` digest agrees.
- Form: the two records are byte-identical (24538/`a1ed6ea3`), and their
  `SWEEP.md` (10283/`be29f13e`), `analyze.py` (41469/`1c9b5c49`) and `run.ps1`
  (3550/`3ac72c43`) digests name versions other than the paper tree's; their
  `grid.py` (3559/`30f555b1`) and `sweep2.cpp` (11309/`57e91b7e`) digests
  match it. The form deposit carries none of these files.

All four records name a `SWEEP.md`, though the deposit carries none in either
tree. Every data record is reproduced by the v5.2.0 run regardless, so the
inconsistency is confined to the provenance's harness digest section. It is
reported here, not smoothed.

## What this means for the paper

No measured number changes between v5.1.1 and v5.2.0. The v5.2.0 execution
reproduces the deposited sweep records except for the record paths and the
run's own identity; it reproduces all 14 tables and figure files byte for
byte; it reproduces every bridge record and both cc-pVTZ tensors, with the
traces differing only in `wall` and the summaries only in the recorded basis
path. The v5.1.1 run remains the record of the numbers; this directory is the
source and execution check that shows the release does not change them. The
v1.0.11 deposit is the v5.1.1 deposit and does not carry this check.

## Files

Top-level records:

- `README.md` — this file.
- `provenance.json` — machine-readable summary: tag, commit, tree, the source
  trees, the 76 evaluated files, the 13 execution steps, the input, sweep and
  table comparisons, and the deposit bundle counts.
- `make_check.py` — regenerates the source manifests below from the two trees.
- `evaluation-sources-compare.txt` — the 62-header and 14-source tables and
  the 15 v5.1.1 provenance digests.
- `evaluation-sources-v520.sha256` — SHA-256 of the 76 evaluated files.
- `tag-tree-blobs-v520.txt` — the tag tree's blob SHA-1 against the
  working-tree blob for each of the 76 files.
- `version-hpp.diff` — the one differing header.
- `source-diff-v520.txt` — every non-artifact file that differs, by category.
- `carried-figures.sha256` — the v5.1.1 result-file digests.
- `harness-deposit-compare-v111.txt` — the three harness tag identities, the
  v1.0.10 to v1.0.11 change list, and the v1.0.11 bundle comparison.
- `deposit-harness-v520.txt` — the deposited harness digests beside the files
  they name, both records per tree.
- `provenance-diff-v520.txt` — field-by-field provenance difference: 17 fields
  for product, 15 for form.
- `scratch-tree-compare-v520.txt` — the 20 input and harness files against the
  paper tree's copies.
- `recorder-diff-v520.txt` — the check's `boysprovenance.py` against the paper
  tree's copy, and `reference.py`.
- `sweep-results-compare-v520.txt` — the 22 sweep-tree files, with the
  differing lines quoted.
- `tables-compare-v520.txt` — the 14 regenerated tables and figure files.
- `tables-before-deposited.sha256`, `tables-after-deposited.sha256` — raw
  SHA-256 of the 14 targets before and after the five table scripts ran
  against the deposited records.
- `bridge-compare-v520.txt` — the bridge records against the v5.1.1 run, with
  the stream inventory.
- `bridge-old-records-v520.txt` — the old side's provenance, the harness
  working tree, and the v1.0.11 tag's v5.0.0 set.
- `bridge-run-summary-v520.txt` — the nine bridge steps and the resolved
  command lines.
- `compiled-tree-compare-v520.txt` — the bridge copy, the compiled translation
  units, the two binaries and the six records' stated version.

The `run/` directory:

- `run_execution_check.ps1` — the script that ran the 13 steps; every command
  is there verbatim.
- `verify_pub_tree.py` — reads the tag, commit and tree from the object store,
  checks the worktree file by file, and writes `pub-tree-identity.txt`.
- `pub-tree-identity.txt` — HEAD is the tag's commit; 710 tracked, 0
  differing, 0 missing, 0 extra.
- `harness_refs.py` — reads refs, commits and blobs from the repository's
  object store without git.
- `harness_compare.py` — the bundle comparison behind
  `harness-deposit-compare-v111.txt`.
- `compare_scratch_tree.py`, `compare_provenance.py`, `compare_sweeps.py`,
  `compare_tables.py` — the input, provenance, sweep and table comparisons.
- `compare_compiled_tree.py`, `summarize_bridge_run.py`,
  `compare_bridge_old_records.py`, `compare_deposit_harness.py` — the bridge
  and deposit records.
- `show_recorder_diff.py` — the recorder diff.
- `merge_per_order.py`, `make_tables.py`, `make_endtoend.py`,
  `make_frontier.py`, `make_prediction_record.py` — the five table scripts,
  byte-identical to the paper tree's own copies.
- `logs/` — `execution-check.log`, `compiler-path.log`,
  `compiler-version.log`, the thirteen step logs, the five
  `deposited-*.log` table-script logs, and `before-rerun/` (the two provenance
  records from the check's first provenance step, 18:31:28Z and 18:32:36Z;
  they differ from the final records only in `flags` and `generated_utc` — the
  first recorded the compiler arguments in the invocation order, the final
  records the runbook's `/Fo`/`/Fe` form). `sweep1-run.log` and
  `sweep2-run.log` are empty because the two programs print nothing.
- `sweeps/common/` — `boysprovenance.py` (the adapted recorder),
  `reference.py` (identical to the paper tree's copy) and `replica.hpp`.
- `sweeps/product_amplification/` — `SWEEP.md`, `sweep1.cpp`, `analyze.py`,
  `grid.py`, `run.ps1`; `data/` (`documented.csv`, `grid.csv`,
  `reference.csv`, `reference_validation.txt`); `results/` (`cells.csv`,
  `admission.csv`, `replica_gain.csv`, `per_order.csv`, `band_summary.csv`,
  `gain_up.csv`, `growth.csv`, `decomposition.csv`, `summary.txt`,
  `provenance.json`); `build/` (`sweep1.obj`, `sweep1.exe`).
- `sweeps/amplification_form/` — `SWEEP.md`, `sweep2.cpp`, `analyze.py`,
  `grid.py`, `fits.py`, `run.ps1`; `data/` (`fits.csv`,
  `fits_validation.txt`, `grid.csv`, `reference.csv`,
  `reference_validation.txt`); `results/` (`paths.csv`, `cells.csv`,
  `admission.csv`, `per_order.csv`, `band_summary.csv`, `prediction.csv`,
  `prediction_scan.csv`, `prediction_record.csv`, `transfer_cells.csv`,
  `summary.txt`, `provenance.json`); `build/` (`sweep2.obj`, `sweep2.exe`).
- `tables/` — the 13 regenerated `.tex` files; `figures/per-order-ratio.dat` —
  the regenerated figure data.
- `__pycache__/`, `run/__pycache__/`, `run/sweeps/common/__pycache__/` —
  Python bytecode caches, not records.
