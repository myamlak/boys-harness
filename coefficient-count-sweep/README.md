# Stored-coefficient counts

`count_sweep.py` measures the smallest stored coefficient count each fit
family needs on `[x0, x0+W)` at a target. It imports the library's generator
rather than reimplementing the fits. The deposit carries one harness in two
revisions: the portable `count_sweep.py`, the documented default for new runs,
which takes the library from `BOYS_SOURCE` or a `boys` directory beside this
repository; and `count_sweep_recorded.py`, the exact revision the records
name, deposited so they can be regenerated and verified against it. See
"Harness revisions" below.

The current releases are boys tag `v5.3.1` and harness tag `v1.0.34`; the
kernel and harness revisions the records name are provenance, not current
references.

## Commands

```sh
python count_sweep.py --plan
python count_sweep.py --static --out static_counts.json
python count_sweep.py --sweep --out-dir results \
    --targets rel:1e-13,abs:5e-14 --widths all
```

`--plan` prints the width grid, the targets and the solver settings without
running a fit. `--static` reads the shipped header and writes what it stores.
`--sweep` measures the smallest count per family; it is resumable and stops
before starting a new row once `--max-seconds` is spent.

The stability pass checks each row's accepted pair and every candidate pair at
the stored count below it, at the row's own settings, then evaluates every
saved fit on the 4,001-node grid. It runs a doubled-settings solve only for a
pair whose residual is within a factor of two of the target; the 2026-10-22
pass had 15 such pairs. `--plan` prints the source rows, the pair-job count,
the recorded margin flags and the projection without solving; without
`--plan` the script runs the pass. There is no `--run` option, although the
script's docstring names one.

```sh
python stability_save_verify.py --plan --source \
    records/core/sweep_rel_1e-13.json records/core/sweep_abs_5e-14.json
python stability_save_verify.py --source \
    records/core/sweep_rel_1e-13.json records/core/sweep_abs_5e-14.json
```

Without `--source`, the script reads `results/sweep_rel_1e-13.json` and
`results/sweep_abs_5e-14.json`, the measurement tree's paths, which the
deposit does not have. `records/stability/pass_2026-10-22.json` covers 30
rows: fifteen widths at `rel:1e-13` and the same fifteen at `abs:5e-14`. The
five coarse widths in the core records (`2.8567491415463806`,
`3.405194117991412`, `3.7994698932176263`, `4.091271749574425`,
`4.947466945561672`) are not in that pass.

`cold_ladder.py` re-solves the one row whose verdict moved, `rel:1e-13|2^5`.
Its record is `records/stability/followup_cold_ladder_2026-10-20.json`.

## Harness revisions

The records identify the harness as `count_sweep.py` with SHA-256
`a5def761c8f18d585396210b28864dc6c4a108d2cfda4f3788d87a7f6c854f10`
(`meta.reused_harness` in `records/ladder/recertification.json`;
`meta.harness_sha256` in `records/ladder/rational_control_symmetric.json` and
`records/ladder/rational_unrounded_control.json`). That file is deposited
byte-exactly as `count_sweep_recorded.py`, so the recorded runs can be checked
against their own revision; the portable revision stays the documented
default.

The deposited `count_sweep.py` is a different revision, SHA-256
`fcb45a434fd973446c5bdf50c4fa8348d9a54146fa16f72e129d8ab7a36f1887`. It has
no `cheb_degree_mode`, no `--cheb-degrees` and no `--rat-k-max` option, and
its Chebyshev ladder starts at degree 4 (`CHEB_MIN_DEG`). The recorded
revision's `all` mode also tries degrees 0, 1 and 2; the recount records four
cells with `new_m` 2, all in mode `all`. The portable revision also cannot
reproduce the ladder's rational cells with `k` above the generator's
`RAT_K_MAX` of 8 (the ladder records `rat_k_max: 12`). The records that name
a harness name the recorded revision; the portable revision is the documented
default for new runs.

The deposited `stability_save_verify.py` is the revision named by
`records/ladder/rational_control_symmetric.json`
(`meta.stability_script_sha256`), SHA-256
`51c4ec838819967f13560d1a68804c89cbf54c648f0effd8ea0be5a065e6483d`.
Five drivers import `count_sweep` by name: `cold_ladder.py`,
`stability_save_verify.py`, `recertify_ladder.py`,
`rational_control_symmetric.py` and `rational_unrounded_control.py`. To run
them as recorded, place `count_sweep_recorded.py` beside them as
`count_sweep.py`. `rational_control_symmetric.py` also calls
`cs.cheb_degrees`, which only the recorded revision defines.

The recorded revision's default library root is the measurement machine's.
`count_sweep_recorded.py` and `recertify_ladder.py` take `--boys-root`;
`cold_ladder.py` and `stability_save_verify.py` use `count_sweep.DEFAULT_BOYS`
and have no override; `rational_control_symmetric.py` uses its own `RAT_BOYS`
and `CHEB_BOYS` constants; `rational_unrounded_control.py` takes `--boys`.

## Record

`records/chebyshev_recount.json` is the Chebyshev side of the sweep re-solved
with the corrected split-Clenshaw evaluator. Its `summary` gives 95 unique
cells, 34 changed and 61 unchanged, with deltas of -1 and 0.

Each of the 95 cells carries `old_stored` (the count in the earlier record),
`new_stored` (the corrected count), `delivered_rel` and `delivered_abs` (the
delivered error). No cell carries a rational count. The summary's
`group_b_sums` gives, at `rel:1e-13`, `chebyshev_old` 57, `chebyshev_new` 55
and `rational` 46; at `abs:5e-14` it gives the same Chebyshev values and
`rational` 44. `shipped_totals` gives Chebyshev 55 and rational 39. The static
narrow-B count `kNarrowRatBStored` in `records/core/static_counts.json` is
also 39. The recount record does not reconcile the group-B rational sums
(46 and 44) with 39.

The file states the library revision it was measured at in its own `boys`
block (`boys.version` and `describe`) and the digests of the generator it was
measured with
(`11185fba03abe74d01b30f093818f48972eb8d88413d446b505198b869115180`), the
coefficient header
(`8ce5e17157a4be5f621d695f377094d212cb839732e278d47f23aa225b4d9204`), the
five records it was re-solved from (`paper_records`), and the generator the
earlier counts came from
(`ccc8b011f4170ce0b530bc4781b09670f42839fa18caaa38ea7d3a9ea9ba9bb4`).
Both generators are deposited under
`records/chebyshev_recount/generators/`: the revision the recount was measured
with and the revision the earlier counts came from, with their tags, commits
and digests in the README there. The five `paper_records` digests are of the
measurement-tree originals.
`records/chebyshev_recount/paper-record-digests.tsv` maps each to the
re-serialised copy deposited at the recount and to the file in this deposit:
the `ladder` copy is still a re-serialisation of the original, while the four
core-side records were re-measured at the revision their own `boys` blocks
name and carry the recount's corrected counts.

The five records hold 135 rows: 40 in `records/core/sweep_all.json`, 20 in
`records/core/sweep_rel_1e-13.json`, 20 in
`records/core/sweep_abs_5e-14.json`, 15 in
`records/core/sweep_rel_1e-14.json`, and 40 in
`records/ladder/sweep_all.json`. The four core-side records were re-measured
at the revision their own `boys` blocks name; the `ladder` record was not. The
ladder record regenerates only with `count_sweep_recorded.py`; the portable
revision does not reproduce its `k > 8` rational cells.

## Ladder records

Two records disagree on the cell `rel:1e-13|2^5` (width 32).
`records/ladder/recertification.json` is the final record. It was generated on
2026-10-24 (`meta.generated_utc`), and its `scope.claim_cells` and
`ladder_analysis_scenario.cells_applied` use cells 1, 2, 4 and 6. Cell 6 is
`rel:1e-13|2^5`; its `record_stored` and `certified_stored` are both 16, and
its `why` records that the stability pass moved the cell to 17 and that the
`k <= 12` walk decides 16 or 17.
`records/stability/followup_cold_ladder_2026-10-20.json` is the earlier
cold-ladder follow-up to the 2026-10-20 stability pass; its `new_stored` for
the same row is 17. The deposit gives no reconciliation of the follow-up's 17
with the certified 16 beyond that note.

## Deposited record set

The paper's data statement and supplement name these record blocks; the F1
analysis beside the paper names the F0 control. Each block is copied verbatim
from the paper's record set, with the source column giving its path there,
and no file was edited. `baryrat` 2.1.2,
the BRASIL rational minimax solver, produced `results-corrected/` and
`results-second-edge-a2/`; it is independent of the kernel's Remez and Lawson
implementations. The remaining
records were produced by the harness revision deposited as
`count_sweep_recorded.py`: their drivers import `count_sweep.py`, the recorded
revision described above, and the kernel's generator, and every rational count
in them is a kernel Remez or Lawson solve.

| Block | Source in the paper's record set | Produced by |
| --- | --- | --- |
| `results-corrected/` | `results-corrected/` | `baryrat` 2.1.2 |
| `results-second-instance/` | `results-second-instance/` | `analysis/tail_instance_sweep.py` and `analysis/analyze_tail_instance.py` |
| `results-second-instance-f1/` | `results-second-instance-f1/` | `analysis/second_instance_f1_sweep.py` and `analysis/analyze_second_instance_f1.py`, with `analysis/f1_reference.py` |
| `results-second-instance-f1-control/` | `results-second-instance-f1-control/` | `analysis/second_instance_f1_sweep.py` with `--function F0` |
| `results-second-edge-a2/` | `results-second-edge-a2/` | `second_edge_a2_solver.py`, `analyze_second_edge_a2.py` and `verify_witnesses.py`, with `baryrat` 2.1.2 |
| `recertification/cell_5.json` | `_recert/cell_5.json` | `recertify_ladder.py --cell 5`, over the recorded revision |
| `results-second-instance-f1/knee-f1/` | `_knee-f1/` | `knee_walk_f1.py`, `make_corrected_record.py` and the F1 analysis scripts |
| `analysis/` | the paper's record set, top level | the drivers, analysers and reference named above, with `ladder_analysis.py` and `validate_tail_instance.py` |

`results-second-edge-a2/` is the independent solver's record at the second
Boys left edge `a2 = 2 x0 = 23.799696304216968`, for the three targets
(`rel:1e-13`, `abs:5e-14`, `rel:1e-14`) over the fifteen Group A widths.
`baryrat` 2.1.2, the BRASIL rational minimax algorithm, produced every cell
under the paper's acceptance rule unchanged; the `x0` side of every
comparison is `results-corrected/corrected-counts.json`. The run is 45 cells
on six workers: 1359.7 s of production wall, 1405.6 s including the two
probes. All 45 witnesses pass the independent re-evaluation in
`verify_witnesses.py`. A count is an existence upper bound under this solver,
not a minimality proof.

`results-second-instance-f1/knee-f1/` is the paper's `_knee-f1/` block, copied
verbatim. `knee_walk_f1.py` produced the two walk records with
`count_sweep.refine_pairs` at the F1 sweep's recorded settings (2,001 nodes,
Remez 60, Lawson 300, `k <= 8`); `make_corrected_record.py` built
`corrected/sweep_all.json` from the delivered record and the two walks;
`bracket_f1.py` produced the paper-width bracket record. `model_table_f1.py`,
`class_diff_f1.py`, `drift_like_for_like.py`, `verify_record_f1.py` and
`make_report_tables.py` produced the remaining top-level records and logs.
`corrected/analysis/`, `corrected/model_table.json`, `corrected/class_diff.json`
and `corrected/tables.txt` are the same analyses run on the corrected record.
The scripts hardcode the paper tree's paths; `verify_record_f1.py` also
hardcodes this repository's `records/` path on the measurement machine.

`analysis/validate_tail_instance.py` is the pure tail's validator named in
section S1 of the supplement. It was at the top level of the paper's record
set and is deposited beside the tail driver and analyser.

`recertification/cell_5.json` is the later re-certification of the
paper-width `rel:1e-15` cell: the record stores 18, the 17-to-19 walk finds no
pair at 17 or 18, and the certified count is 19 (`cell.verdict` is `moves`).

One ladder cell is not independently reproduced. `2^5` at `10^-15` is 25
`[22/2]` in `results-corrected/corrected-counts.json` against the paper's 24
`[21/2]`, whose coefficients are not in the record and so cannot be
re-evaluated.

`analysis/analyze_tail_instance.py`, `analysis/analyze_second_instance_f1.py`,
`analysis/tail_instance_sweep.py`, `analysis/second_instance_f1_sweep.py` and
`analysis/f1_reference.py` import `count_sweep` by name; the two F1 scripts
also import `f1_reference`. Run them from a directory that holds both.

### File digests

SHA-256 of every deposited file, with paths relative to this directory. The
two `results-second-instance-f1/sweep_all*.json` files hold the same bytes and
hash the same, as do the knee block's `model_table.json` and
`model_table_recheck.json`.

**`results-corrected/`**

| File | SHA-256 |
| --- | --- |
| `results-corrected/README.md` | `711ef926c6944231682ac3feb1e8cf368986c1f0ca9c44671e35f8dc13bb6bcf` |
| `results-corrected/corrected-counts.csv` | `e11899412881e23c65756c7949a2e2f9564c8e30ccf33c397b6a5549ba2bc335` |
| `results-corrected/corrected-counts.json` | `9bb8ac0633b10e8a7922684ca8a419135103f2c7825b4d06b9c71603399b7195` |

**`results-second-instance/`**

| File | SHA-256 |
| --- | --- |
| `results-second-instance/_edge-a2/chebyshev_only.json` | `8e642782bafb7d3cb4fc7995e6b3c2cb6c9da4ab621bc1a454ad134345c0f450` |
| `results-second-instance/_edge-a2/sweep_abs_5e-14.csv` | `f317e124b033ba15b93a05f3242cc730aed6bd374ded37b35eca68d1224845d6` |
| `results-second-instance/_edge-a2/sweep_abs_5e-14.json` | `f61870429e4ae3a93cc81abc8608038a498ed1e18ca1761d140ab521195e3322` |
| `results-second-instance/_edge-a2/sweep_all.json` | `e1cec20cf988cad3dd0fa81d4c4ea86340b01ae381841849794e84cb02e0ce03` |
| `results-second-instance/_edge-a2/sweep_rel_1e-13.csv` | `cde7ce46ad1fa5dee8c39161d229a2889a6065daa1b2fa6bbd8cb7eb774b61f6` |
| `results-second-instance/_edge-a2/sweep_rel_1e-13.json` | `2f1a39e2a491ae2c8cc6a0f51524abe35c71d02d5429c475eb9cc6c96d175f66` |
| `results-second-instance/_edge-a2/sweep_rel_1e-14.csv` | `5b06bfc2d45f593e80390f4de9d0091f1a19cf09f1d690737661025394276628` |
| `results-second-instance/_edge-a2/sweep_rel_1e-14.json` | `44d8200b94fb18702078a602e285149aae25d2c1786bc84adfb9313b8483ef25` |
| `results-second-instance/_edge-x0/chebyshev_only.json` | `894490bbc933b3a65231cc4a36cf12592cd9a4cf5d770bd6d8cbb4d5e1f58d46` |
| `results-second-instance/_edge-x0/sweep_abs_5e-14.csv` | `cbcc6fad4c9e2d642cb8f925f3ce2d3bac4052f851b6172ff486a9d5a069e58e` |
| `results-second-instance/_edge-x0/sweep_abs_5e-14.json` | `20e850d28db75a3818a73b4c97621b9626618c811c4934c987186fa44de45f51` |
| `results-second-instance/_edge-x0/sweep_all.json` | `6c33aee9988f27f2fdd0afef2332113817000114e6e0b2f0e23b3395a067c973` |
| `results-second-instance/_edge-x0/sweep_rel_1e-13.csv` | `76437901b45a4514a9f13da0a821a4a53181088e633f8c579c43b208655c7b24` |
| `results-second-instance/_edge-x0/sweep_rel_1e-13.json` | `b35b5c4b612782603d4170625af30b9a85934397018c329660e01ac3fb2323ec` |
| `results-second-instance/_edge-x0/sweep_rel_1e-14.csv` | `da316a42799d37e8993c5e7a6e9b972bda97a06870c055e2ff3d76e4ce258c1c` |
| `results-second-instance/_edge-x0/sweep_rel_1e-14.json` | `44aac0d32617f952d4e6fb61da24388802b7f5df3958565a845170be5ea35ca5` |
| `results-second-instance/analysis.json` | `44657bbc6aab32c081e9d72c4ef95e1ba660d0701298d68d5a8ed065c027554b` |
| `results-second-instance/chebyshev_only.json` | `e8a7b66d140d875efc88193f6dcd47b1847bbe60099271add159baa0291234b4` |
| `results-second-instance/comparison.csv` | `6627912d29db96aa3e9605af292755df1fea7b5a9814fb030a63720e02980c16` |
| `results-second-instance/control-f0/chebyshev_only.json` | `94b9add10091f369d56bee378f2746feff8477fb5c0dafe0f5851003f84c27d9` |
| `results-second-instance/counts_by_width_target.csv` | `ab55580b600cb5d34b156bd5998577b09c90375a4307b23371ec8a178f6ab3cd` |
| `results-second-instance/edge_a2.csv` | `97e24a6e0d9ec9e0e777e3846d629565f7fdfc10f3a1a97a880f6d8924b14ad0` |
| `results-second-instance/edge_x0.csv` | `4fc79037c2464c07c2f869681e9fd7444a0eb7b609733c287424cb8671d8dfbf` |
| `results-second-instance/logs/_ti-analysis-final.log` | `d8a57b71401dd555cb357b27196c5d1231ad769918e01eff8d28222d098972a0` |
| `results-second-instance/logs/_ti-analysis-nearfull.log` | `c47c204e6090aa75437ca73c870379eec92c9ee676f847ecce34c66ca34cda5f` |
| `results-second-instance/logs/_ti-analysis-partial.log` | `7728f5892e8751d8250794411876d02f282cd68d7043041c6f8dcfc5f7162bbc` |
| `results-second-instance/logs/_ti-cheb-a2.log` | `8439f63f9997099445679df89bcce77c200694e65262aff2b0eef19ef86f7d96` |
| `results-second-instance/logs/_ti-cheb-x0.log` | `aa535dd28a16292ebfc5bae037c2f1478e58ba8dd4908dddc38f11efdcfaac0b` |
| `results-second-instance/logs/_ti-f0-cheb-control.log` | `f6863eee6c55595783d8e38fc47a4ce60e85f97c519c7778530e8b6ad6516fc6` |
| `results-second-instance/logs/_ti-f0-cheb-probe.log` | `4e1cadd3cb2c32965c07cd290619009da4f308b96d1c6230d76f3fae7e4cd9f0` |
| `results-second-instance/logs/_ti-f0-control-compare.log` | `49ef0f3d80409e3dfd4abb898a3e46776bc703d9fc111b899723611469ce9611` |
| `results-second-instance/logs/_ti-fill-x0-abs-ship.log` | `83cb4f697c6fd6264de9d1fa1131b68dce5a4dd4311130258d76688c76181ff3` |
| `results-second-instance/logs/_ti-rat-probe-x0-2m7.log` | `833a64d448e70933b29cf2e65dfc79362427d56a977284352b3704b3e5e286d1` |
| `results-second-instance/logs/_ti-run-a2.log` | `9e1c0621be6844d50cf6f10c31c3b89b69269f2f6f38bb0bba02d963a3b43332` |
| `results-second-instance/logs/_ti-run-x0.log` | `9548c62d677d35a534d5fd8d1df9f3ca1e6883a3eb7d1abb02e780ec3079c43b` |
| `results-second-instance/logs/_ti-validation-recount.log` | `7706e86c718c4fb08b1cd27ed170bba78c8efae342cc4c59948dbafd949525e3` |
| `results-second-instance/sweep_abs_5e-14.csv` | `485e1589dfb22fb08d2261f9f4b16efc66346052078bb734cffd94665875c997` |
| `results-second-instance/sweep_abs_5e-14.json` | `cbf626561a40d441b2a9266cf770f29654b855802fd44cb43185c86f33ad8f82` |
| `results-second-instance/sweep_all.json` | `d0674eb7031d17b919c78e84d083382e4106d52a9af7293660c74aa8d02e3f2b` |
| `results-second-instance/sweep_rel_1e-13.csv` | `0829158cfe1d2e8b5f6571c003db67d9b7c13be33529fcdfdc20e2c1f6eac085` |
| `results-second-instance/sweep_rel_1e-13.json` | `1c45ad3301681ed95ecb0a4132f1be2705c50a58fbe09a450e03becb2b1f13a1` |
| `results-second-instance/sweep_rel_1e-14.csv` | `2194bacc36becab761e1d94ce143e5f0369338b3825612a655ac395eb112679a` |
| `results-second-instance/sweep_rel_1e-14.json` | `43b589d689c2b61a42a30bb55216546bb73e07946139fce50391c51746117a38` |
| `results-second-instance/validation_recount.json` | `a159d481094b563dd019251440ca3fed74a2061075c04ddfc1a4a7241cb0027f` |

**`results-second-instance-f1/`**

| File | SHA-256 |
| --- | --- |
| `results-second-instance-f1/analysis.json` | `aac896b40f71a2c866312da3ca9356b62d1e2563396005f9e437168e7e644c96` |
| `results-second-instance-f1/comparison.csv` | `c8922b3fa6b3a8436dfce17454b8861d5ab8ae8bf724bb8547a36fedb6a0230a` |
| `results-second-instance-f1/counts_by_width_target.csv` | `9c91ce492afd6e7b8e5051b54427948cbde77a9e67a06707a44131a8eb6262d2` |
| `results-second-instance-f1/f1_validation.json` | `03f749180b8538af35fad297ccee14984c3cb55d86d81d9d13391b31ec17f651` |
| `results-second-instance-f1/fitpath_check.json` | `e3dfdfb54f7e74ec2d4e9d703cea36a14089a70b4dfff94fc7a442b7f1aebd87` |
| `results-second-instance-f1/model.csv` | `679185ced0ef1def7ed50917dfb36b086d199360791c3277538a43fd429a0ba8` |
| `results-second-instance-f1/rate_check.json` | `1e81b9565555fdc458387334ecc4443b67c5234e7a3d5a74055022bd41223085` |
| `results-second-instance-f1/sweep_a2.csv` | `2eac065d2f1dd41d0f9e2aa4e3958ccceba4d5e398352a259d2cee9d27180b06` |
| `results-second-instance-f1/sweep_a2_abs_5e-14.csv` | `4805dd429bfc444097db09b8349acfa6dd1bd4295bb4d71f33585a9531737097` |
| `results-second-instance-f1/sweep_a2_abs_5e-14.json` | `85abf038fe0911c65edccdfa72edeaa95750fa3e8872fa1b67bb585ae3f7e5b0` |
| `results-second-instance-f1/sweep_a2_rel_1e-13.csv` | `654f856ccc9ba9cb9cd5157eae7bc4c48dd0269585246b28b6bda60995c459d7` |
| `results-second-instance-f1/sweep_a2_rel_1e-13.json` | `8b5dd769ebdbeb852b343b16df36bf8adc153e5080ba6a4e31204732903ee213` |
| `results-second-instance-f1/sweep_a2_rel_1e-14.csv` | `c7d252fe47bfa75d0e1a8e33f1587d40251602c634262e0b5ea92a5f25d6ffff` |
| `results-second-instance-f1/sweep_a2_rel_1e-14.json` | `2c0a78868d5a48116edd5da1478e168bbe4ac7baabe71938f3c1c3a87ecaf3db` |
| `results-second-instance-f1/sweep_all.backup.json` | `a08cef0933d08d656e8ff46a23d4b9bd5c5f0b54fa42b688dee56407da462fb3` |
| `results-second-instance-f1/sweep_all.json` | `a08cef0933d08d656e8ff46a23d4b9bd5c5f0b54fa42b688dee56407da462fb3` |
| `results-second-instance-f1/sweep_x0.csv` | `55f68125207a7b8bf437570833fd7910aecd076a7986f076f9894389f3b26bf2` |
| `results-second-instance-f1/sweep_x0_abs_5e-14.csv` | `7eda04ef12d5e7ff83660971c9ad4bad83ca88f7beae1d33704c759ea28bbf3b` |
| `results-second-instance-f1/sweep_x0_abs_5e-14.json` | `c684d8f4cbd2334fc4d987be4451198d41a426985154469143452b0aaf340fd3` |
| `results-second-instance-f1/sweep_x0_rel_1e-13.csv` | `5ec0c8b9e8e3e57af23230c9467d8447e4ffc7cbad7dd5b7309f258d0c8e9e20` |
| `results-second-instance-f1/sweep_x0_rel_1e-13.json` | `15032cce7d39be73700132dd209103ba0ca91390a19d4963be601039d314f9b6` |
| `results-second-instance-f1/sweep_x0_rel_1e-14.csv` | `7b840149a5194cf1a8bd1d87e59549ae430597d6a9709695fb6e32c2f02011a8` |
| `results-second-instance-f1/sweep_x0_rel_1e-14.json` | `7487f342c6bf8822f835bb6ff11fc58ae79674d30a13fb035509b94389e55269` |

**`results-second-instance-f1/knee-f1/`**

| File | SHA-256 |
| --- | --- |
| `results-second-instance-f1/knee-f1/bracket_a2-paper-rel13.json` | `56d9a53136728d9bffa917e38b1863948dec2594752dad68673d8455c2356f5a` |
| `results-second-instance-f1/knee-f1/bracket_a2-paper-rel13.log` | `d64c3277b84a725d00f36ab8f4918372cfe60a8c1920306a4c1532da9cdf92fa` |
| `results-second-instance-f1/knee-f1/bracket_f1.py` | `0167525f275a9ddbb9b21a3f1517fd32e8d1f28f07e71cd4b23debfb2fd83377` |
| `results-second-instance-f1/knee-f1/class_diff.json` | `4ee6ce89541145b60a2e169e454e9fccade08dd5451cba17a7faf135629466b5` |
| `results-second-instance-f1/knee-f1/class_diff.log` | `e90fb8b4c82d979ff772f6455c0bcf6e21365ea18dde9d1bb0f3f57e6713c3bf` |
| `results-second-instance-f1/knee-f1/class_diff_f1.py` | `93ade6aa5573b1e66d761fd6261245492700e9c6c0d2f00de0543c1d4e5429ff` |
| `results-second-instance-f1/knee-f1/compare_f0_control.log` | `49ef0f3d80409e3dfd4abb898a3e46776bc703d9fc111b899723611469ce9611` |
| `results-second-instance-f1/knee-f1/corrected/analysis/analysis.json` | `c85c485c2d8ab77229408bc1a78d6ecc2683aa442a514a2a4a95a0f7395b01ad` |
| `results-second-instance-f1/knee-f1/corrected/analysis/comparison.csv` | `a9221e4c51b68ae717c88308251d5cabf4e0c9e3004766a4970677a57a2f733d` |
| `results-second-instance-f1/knee-f1/corrected/analysis/counts_by_width_target.csv` | `054bc69ae67f90cac2c6a1f244581865fd7b1ea54e206131a5693e4acd6360e4` |
| `results-second-instance-f1/knee-f1/corrected/analysis/model.csv` | `1ac7cfc67676a6fdb461bbc94649b38084f9d8897783bec6ac48fb737246e2b6` |
| `results-second-instance-f1/knee-f1/corrected/class_diff.json` | `ff2c72b5e3466658f89d8ce6890b7c975ab87b97cb9ce41ecf6b9490f52b3f86` |
| `results-second-instance-f1/knee-f1/corrected/model_table.json` | `0701cb75ba6523d35da11e461eece440bd95d1847b20bf21ab728fb85068b914` |
| `results-second-instance-f1/knee-f1/corrected/sweep_all.json` | `42b8f689ddb4f20fa6a6a5e60c171641f9a6fb7bd6edcbcda92d484e88d0a173` |
| `results-second-instance-f1/knee-f1/corrected/tables.txt` | `d633ab41b6d8a333ed669e0413fa7de7e64bf6ee754151c2338b15167c91e5a6` |
| `results-second-instance-f1/knee-f1/corrected_analysis.log` | `af8341a23a5d7b10b636f0020e2ee1ddcf5c380530530970d9dbd5d333a42ebd` |
| `results-second-instance-f1/knee-f1/corrected_class_diff.log` | `b8458bee5d51d06a3cb1c8d418130db396c89a6e8b02bbe7206c3ba292eefcaf` |
| `results-second-instance-f1/knee-f1/corrected_model_table.log` | `9401f70452e247b76ead400739692bae0b90a7b876e8735abef989fa09612171` |
| `results-second-instance-f1/knee-f1/corrected_record.log` | `f15f074f9d8bb4e6008d9bd8405d53cea85df141044004631076e32770929dac` |
| `results-second-instance-f1/knee-f1/drift_corrected.json` | `5df5c8c4a950717b979ebe722aca3aac1d83c8d5d6b0aa1b965edbe3132aa967` |
| `results-second-instance-f1/knee-f1/drift_delivered.json` | `ae3b1eb285836ab16b6a9dfcadd4872536482344fe2f3e3e6fd44c3717dfd8bd` |
| `results-second-instance-f1/knee-f1/drift_like_for_like.py` | `0f2c7831c6a2e1f8273a2dd7db3f9ed13326ee982fc0ed4d374f4ee1d666cd68` |
| `results-second-instance-f1/knee-f1/knee_walk.log` | `3ba84ca3850a87b9945ce3e481d7aadd47d3f9ab70155c3cbe3e97ce707f1150` |
| `results-second-instance-f1/knee-f1/knee_walk_a2-2^5-rel14.json` | `9561b6ca6053d3c4e287548e0043fdbe8c96d1db7878a092b63901ae3f65b1f2` |
| `results-second-instance-f1/knee-f1/knee_walk_a2-paper-rel13.json` | `3c50e25c243029ee81b783ce3235d764eb008261702cf09887cec59368b8ca24` |
| `results-second-instance-f1/knee-f1/knee_walk_attempt1_badroot.log` | `91422f3bf268829a82d1ac95f3297c0a457444a44d19fe40035e1f8a7e1af402` |
| `results-second-instance-f1/knee-f1/knee_walk_f1.py` | `e5750e23c77cf619a4665e08b4d62ff681c2c5bf72c659bd9ebdbf4b872f3510` |
| `results-second-instance-f1/knee-f1/make_corrected_record.py` | `137ccec9d1e553322a5e4de858e79872225f4be497621ca3a62435a04ab13085` |
| `results-second-instance-f1/knee-f1/make_report_tables.py` | `ac7332c441c7d287b770b6f81026a85be49d39d9414d562f8341c59e8965c757` |
| `results-second-instance-f1/knee-f1/model_table.json` | `257073d19b99a9785e6fe63648a4e271ec48dedab39e48984da8327d76860174` |
| `results-second-instance-f1/knee-f1/model_table.log` | `1446e578d8ed348bf880598be165e9f74b790510077f3fedf1cc067351703132` |
| `results-second-instance-f1/knee-f1/model_table_f1.py` | `97081c42d63cfe46256e0e074f76c6a771e8ce767f3d87213a82a508a79fc0bc` |
| `results-second-instance-f1/knee-f1/model_table_recheck.json` | `257073d19b99a9785e6fe63648a4e271ec48dedab39e48984da8327d76860174` |
| `results-second-instance-f1/knee-f1/model_table_recheck.log` | `5705d41d7331190a9363cc8f6b8dcbd3ec55f05fe12eddaf140a838115ebe737` |
| `results-second-instance-f1/knee-f1/tables_delivered.txt` | `15e5fa65ec5bca2e5562ae4d201874b9d5fcc27e47cd839ae9b61c63007283a8` |
| `results-second-instance-f1/knee-f1/verify_record.json` | `8f3351d40936a733aa7b7db0df65527826b02da7537d65ea6530183e0b8a6a5f` |
| `results-second-instance-f1/knee-f1/verify_record.log` | `28946404f33e79047cc2e11b724e940369ccb52a96e5b2a23fd5ed67d819a01f` |
| `results-second-instance-f1/knee-f1/verify_record_f1.py` | `855b4630a30d2418c698a183f798c9173e5edc06a9303f9ba7df8d38a35a874f` |

**`results-second-instance-f1-control/`**

| File | SHA-256 |
| --- | --- |
| `results-second-instance-f1-control/sweep_all.json` | `e316596c308fe42edd83f3a6ab1df6f0da1751d4c85c52e872907ff62b52d80e` |
| `results-second-instance-f1-control/sweep_x0.csv` | `3bc6ca4f7099cee9ee844caf6287a4924d878b2bf3b2419fcc935f9c808863c1` |
| `results-second-instance-f1-control/sweep_x0_abs_5e-14.csv` | `d538fd3fa8d1ddc2750aaaee30dffc935db17f6c66743a5002d1cbfa26be2736` |
| `results-second-instance-f1-control/sweep_x0_abs_5e-14.json` | `5a60207a6fdaf8fdb7f1ab775dc858de15666ba0e476e31d18a67f2ab05c591f` |
| `results-second-instance-f1-control/sweep_x0_rel_1e-13.csv` | `a26986f42c315fc2a04e7a1aa1c82bff5f1eba098f6948d8089ac7eb04adad3c` |
| `results-second-instance-f1-control/sweep_x0_rel_1e-13.json` | `32a730e4d11e44009707806fd47e6b78d7276c3d4ecbf9028b06215cc37fc675` |
| `results-second-instance-f1-control/sweep_x0_rel_1e-14.csv` | `bb426d6109b5bcfce583381d5f7f29cd7b397395a8015e9ed779c7660c9be1e2` |
| `results-second-instance-f1-control/sweep_x0_rel_1e-14.json` | `cb7225b0a211cbe1b34164d5c946b9cff88b6c81cbac316adc04e2194bd0f796` |

**`results-second-edge-a2/`**

| File | SHA-256 |
| --- | --- |
| `results-second-edge-a2/a2_run_summary.json` | `fbdc935ace63fe6fda9867b0e5e57f36ddb35dd285869bcb72c0cfeb6f34149b` |
| `results-second-edge-a2/analysis.json` | `f49f8aa67c845d427ad43a1b0ded76002638996a0fdf43ec18d2ba8c12cd6b4a` |
| `results-second-edge-a2/analyze_second_edge_a2.py` | `2de78c6158544e05b1628398dd6ca26f1dd1d526f622172ead219655c9d1df0f` |
| `results-second-edge-a2/cells/A2_abs514_m1.final.json` | `d9b51b9db28d38286991bb2d511708dc3633204ec3431741e9479dc687b85eda` |
| `results-second-edge-a2/cells/A2_abs514_m1.json` | `9a0596c2c6f4487a8c7989f2a50a0463b7997177ad364a8ef56508531f17e40a` |
| `results-second-edge-a2/cells/A2_abs514_m2.final.json` | `f3567d749277517d7d78c326eba60eefa6c2db4da5cdce33b47436f1db779f93` |
| `results-second-edge-a2/cells/A2_abs514_m2.json` | `4ef1c681ca6160862a3baf35f81a6778f66f35cd7c7ce8d2d5105851beca86df` |
| `results-second-edge-a2/cells/A2_abs514_m3.final.json` | `77d2cd9206b9ac14d31c951e66f5c7a20a8645fc6703c8de0befee4e62545373` |
| `results-second-edge-a2/cells/A2_abs514_m3.json` | `82bd3d572350af46f66307ef481f5b8c83563a6c492d50cd43bfcddff80d4d1f` |
| `results-second-edge-a2/cells/A2_abs514_m4.final.json` | `cbea23cb23c26f9e4b7d3d77c852346b866a2bcb32b342b9216cbd43adf309a5` |
| `results-second-edge-a2/cells/A2_abs514_m4.json` | `7a9b5c8d20261497c65f2d2c82065ac1295114d4c23e656110784d899a8e9d54` |
| `results-second-edge-a2/cells/A2_abs514_m5.final.json` | `555dba18c58735bd333d178d6acd62bf2754d8398c2abd259e00b92243e6c7f6` |
| `results-second-edge-a2/cells/A2_abs514_m5.json` | `a2c77e47a37056303051efeb6aafe84d0a6e77ef8bba4e3669e71c8fd65e653c` |
| `results-second-edge-a2/cells/A2_abs514_m6.final.json` | `537d28c573318d91f935e3514d0a134a689ebe1a7e5d1be02c82fc429233eba3` |
| `results-second-edge-a2/cells/A2_abs514_m6.json` | `e6c55173945e0c21f22e211014eb785edf9f587e28468de24420cc1fa6140698` |
| `results-second-edge-a2/cells/A2_abs514_m7.final.json` | `5329c818d761296650ff1bc1d0e1c8f4f8f65df195212683fc5e06f4423bbba0` |
| `results-second-edge-a2/cells/A2_abs514_m7.json` | `0e084f25dfe18c17f46bd7cbf58f0144a618f4ff80b1c0323b4bef6a706364f2` |
| `results-second-edge-a2/cells/A2_abs514_p0.final.json` | `6de9d2a32128556e259943bfaadd1ec2c1b323b4e4427e4ac4930774e7111ea6` |
| `results-second-edge-a2/cells/A2_abs514_p0.json` | `f9cac70104380118aec919acb1f44a1488ef0761b63dd4f28ff4741cfc9f18dc` |
| `results-second-edge-a2/cells/A2_abs514_p1.final.json` | `947a26e8a7d0db61dc7639ada90b99eff9379ffa89a82cbecd6d510e30cc5d81` |
| `results-second-edge-a2/cells/A2_abs514_p1.json` | `028b75b49c5d30e2bb9b364f523313a9ccac8358e8c91dadb7328109df88fcde` |
| `results-second-edge-a2/cells/A2_abs514_p2.final.json` | `efde211736d5858b05c8fadfea7a84728b33a3d2f53edd7f156e7851fb27c051` |
| `results-second-edge-a2/cells/A2_abs514_p2.json` | `0fbc1790f3aa29697078b3a8d52ef5b7e945bd79c3ac60108d2a00360c28a3ca` |
| `results-second-edge-a2/cells/A2_abs514_p3.final.json` | `f17677164e3a72bad9769d47fea1ea1f8623beb5a23932399e6e631808632d45` |
| `results-second-edge-a2/cells/A2_abs514_p3.json` | `05ca8cddc0bbc55fed59c46bed9f0e65e953134f75ce63f1e3bdc529b58794b8` |
| `results-second-edge-a2/cells/A2_abs514_p4.final.json` | `c444d29db1c1c83b1b80cd4f15176b482118c03744a28ade9007cd5b2a31496f` |
| `results-second-edge-a2/cells/A2_abs514_p4.json` | `48f040e86666ddbbe527a68e4dc6449223cd46841a8d4afd71783be051ffa6ef` |
| `results-second-edge-a2/cells/A2_abs514_p5.final.json` | `fb3b48a5853415abef28c5523e6d2f47bf7834623a99cc35a4ae2ec6f0ace54e` |
| `results-second-edge-a2/cells/A2_abs514_p5.json` | `4927b030a28fe5a44a39f9cd6518358d56a852f70f93c33e6e6086d8dd74102a` |
| `results-second-edge-a2/cells/A2_abs514_paper.final.json` | `5512a09971291123957f2da3fff3c86620748356bebd148f053d16487053a9d1` |
| `results-second-edge-a2/cells/A2_abs514_paper.json` | `04f760b0cf8f0952e03f19c52f37ae1c647b48f381b797ab0c3e80054927b28f` |
| `results-second-edge-a2/cells/A2_abs514_ship.final.json` | `a69ad5ac0f0323bf5a842faa5831ae890f664ca620cd65f8b71885ac6bd56d89` |
| `results-second-edge-a2/cells/A2_abs514_ship.json` | `fe39116fc0ea805e2e02b7d11f928b3b2c116a357a5f2391bc73d2262f8032ec` |
| `results-second-edge-a2/cells/A2_rel13_m1.final.json` | `56f246bbbb5dbd8a0b88e608f2661f3e7bb2a3d8f54dcb846272102cb6f77e78` |
| `results-second-edge-a2/cells/A2_rel13_m1.json` | `8892caf98808e3e8736fc45e08b6aaa9fbff9520b6001ed311b2321173800ad7` |
| `results-second-edge-a2/cells/A2_rel13_m2.final.json` | `117282fd7ebbe7fa347f19632e66940d043069ed3df76d948e2f1306a0d27597` |
| `results-second-edge-a2/cells/A2_rel13_m2.json` | `aa06ff4ed9c1b2336f9e055ca43c5486cc7f26572eadb0ee530e0e1469a7ee8b` |
| `results-second-edge-a2/cells/A2_rel13_m3.final.json` | `fb0df4612f8c3d5e6af1e4809061874b5f60424bc6cb2f52eb463f22dfc382eb` |
| `results-second-edge-a2/cells/A2_rel13_m3.json` | `9fa8f075a49fe71f6e74f2229ea93d4584c9922e4772d7d1e3c5f80e1293c5bd` |
| `results-second-edge-a2/cells/A2_rel13_m4.final.json` | `2c915b35ba47425910e7d689c4038ce595cf42bd143faefb03f4aed7a82c306b` |
| `results-second-edge-a2/cells/A2_rel13_m4.json` | `3693f8a65cff7ddeb0c5f7f79315cc7373fc3834b8b7e639ae833a90f6456a8b` |
| `results-second-edge-a2/cells/A2_rel13_m5.final.json` | `1763628edc6ea7de025636bc754859dd4e32689d280c934feab54bdcfeaca869` |
| `results-second-edge-a2/cells/A2_rel13_m5.json` | `59c08984816fdcb9934db304a6922cd6cb43c063ece9a720a80575d5ea1eba7c` |
| `results-second-edge-a2/cells/A2_rel13_m6.final.json` | `298caa0ebbfcea0b1edf3c688d3546cde1f195ac6a1ab4fecde1aa3194dffd4e` |
| `results-second-edge-a2/cells/A2_rel13_m6.json` | `784bf3d564ac31ef7c2235f89b28e7cc80ad660c9024e2414ebf48e2b138c3f6` |
| `results-second-edge-a2/cells/A2_rel13_m7.final.json` | `90dddd2b7bf940175e4f6932e02e878a4c5243b468f5a65cfc3c39ed22080dc1` |
| `results-second-edge-a2/cells/A2_rel13_m7.json` | `859091300336102ef462c6b8b5caa7c0c12e44777f2da9375a07d3fb8aaf3bff` |
| `results-second-edge-a2/cells/A2_rel13_p0.final.json` | `ebbee9c081586a0245cb14d68ae5a7dcb78dfb3b89557653d0c449fdcf250f93` |
| `results-second-edge-a2/cells/A2_rel13_p0.json` | `9cd2023c0ddbe0239cc69ffe1a2c19c1bf6958a970eaf1d7df33d1d29c231d94` |
| `results-second-edge-a2/cells/A2_rel13_p1.final.json` | `15e7c07670a1c593abad6044ef18c0ca05f866a369fdcbb266a04e4bc4253648` |
| `results-second-edge-a2/cells/A2_rel13_p1.json` | `8ac1c50056e34c1f7ee652e31f656ceacd7a1bae36d5b6e18a5191f74296c2a6` |
| `results-second-edge-a2/cells/A2_rel13_p2.final.json` | `b82e1d3a7e113b0f23d7ddd23203b7902e53a2e5614f4b32bbd124043b888eca` |
| `results-second-edge-a2/cells/A2_rel13_p2.json` | `c34892088beb2b77d21ab8eeceb7b2823da908adfb8bc506091820a09a23def4` |
| `results-second-edge-a2/cells/A2_rel13_p3.final.json` | `830c4dbf1dea48704b4006a09d901238471261f9eb752a8baa2be389e1907e35` |
| `results-second-edge-a2/cells/A2_rel13_p3.json` | `65f834ca875f4891f2894c048ba74a200440a500b3e62a7d7365d2e77422018a` |
| `results-second-edge-a2/cells/A2_rel13_p4.final.json` | `1d11e8b7bc0a42d9370ce336e3b199808bafcbb28069e666e9439bafd75f2acd` |
| `results-second-edge-a2/cells/A2_rel13_p4.json` | `3d48b28dd1e7436f132c177bff3414be2fbad9d1d9bb80a6ca05b0a5ee24fd24` |
| `results-second-edge-a2/cells/A2_rel13_p5.final.json` | `6a11f64a7ecfc7162be1b97ceda92d68f7ceb91aa9088b69a47282e8569747b9` |
| `results-second-edge-a2/cells/A2_rel13_p5.json` | `ad2ecc14648d60803261b060bf0cb384d5084574956c1c9bea9c8bfd702364c0` |
| `results-second-edge-a2/cells/A2_rel13_paper.final.json` | `d93159d20dad87e004cff31851923193a12ef307039faaa83b7330ea8462353e` |
| `results-second-edge-a2/cells/A2_rel13_paper.json` | `053274afc944f9bd5d92e3c482dbc86eda182bff2c869b1309fbebb066d0c822` |
| `results-second-edge-a2/cells/A2_rel13_ship.final.json` | `5d9d0842c36e40323f29f4e03fe78a4383b8482f5611e898c15d7492e73142f4` |
| `results-second-edge-a2/cells/A2_rel13_ship.json` | `7c872b04bbc1095b20d7470b7e9c0ecb16cfe29607184f9f47c8f7440f74a03e` |
| `results-second-edge-a2/cells/A2_rel14_m1.final.json` | `2ce7f2da8a227e1f68f9be757f26606481456a6e189ce97457a663999625718b` |
| `results-second-edge-a2/cells/A2_rel14_m1.json` | `45c924463a7d78132032a8d29a6347c65d358f13b0de1590db4190145e8756d8` |
| `results-second-edge-a2/cells/A2_rel14_m2.final.json` | `a13b7602a1f5bb25185f8a1ddcfdf989a8306f6521252709c4e8b86cabc0f6a3` |
| `results-second-edge-a2/cells/A2_rel14_m2.json` | `a23f04fc93783607475520808d095bcdbdd00c0ad9bda2a406a6bc9c9d8494f0` |
| `results-second-edge-a2/cells/A2_rel14_m3.final.json` | `93bd42b78ddf1a5a92b46f8ec590d79a6ed112b17db82bc44e7202ae0996255c` |
| `results-second-edge-a2/cells/A2_rel14_m3.json` | `bc0c0bc6beaa625b7835d2c2a3871db6b9d717ecf81fab47995bfaeaadf9a0dd` |
| `results-second-edge-a2/cells/A2_rel14_m4.final.json` | `461dfc3c25309a48bd98ae7e178e32ee9d7f5d63e8d0266851c763ae344ef72a` |
| `results-second-edge-a2/cells/A2_rel14_m4.json` | `0b763fcfc3178ed4244d7a17aa69991a212d4d55b62136411a1e1d824d7f4e4e` |
| `results-second-edge-a2/cells/A2_rel14_m5.final.json` | `02efe331223e968a5654d1753f0f2a9d0ff86793ae2cf268d6d7084872ae8234` |
| `results-second-edge-a2/cells/A2_rel14_m5.json` | `afe2fbea0967e389287078caebaad48f27fbd035c504c78a753c9b443a37344d` |
| `results-second-edge-a2/cells/A2_rel14_m6.final.json` | `56d6eaedd28a239e92bf49b5806964208b8252913d55d4bf02bf688eb9f89580` |
| `results-second-edge-a2/cells/A2_rel14_m6.json` | `2603d92ed3d01e59b874500574af8745150538b4f550b912d2233c06f949b67a` |
| `results-second-edge-a2/cells/A2_rel14_m7.final.json` | `6a1c3e4238984475b4e6a92c0a421e8f8e68661dd958be2c966df32fe3ae2164` |
| `results-second-edge-a2/cells/A2_rel14_m7.json` | `ef60b68f287b7077c4fa56a44c97510e39b6b46461946922fc695fdca980fb26` |
| `results-second-edge-a2/cells/A2_rel14_p0.final.json` | `dee5735ded1a7139e360d1f543d80293804b6a315c4d27faa41295895ae20840` |
| `results-second-edge-a2/cells/A2_rel14_p0.json` | `a2b50f999b7928e128a0b6f11481cd9666e445023c9cf20c24d96da68ca81854` |
| `results-second-edge-a2/cells/A2_rel14_p1.final.json` | `4883d586a26a4887db302fa759bfa91d73ff7bef3721a81fbe8da4ea2f72cffd` |
| `results-second-edge-a2/cells/A2_rel14_p1.json` | `cb5a1f532b89d96b7bf036f7ae73bf8485d452a73c58bd21b30e2329297cdf9b` |
| `results-second-edge-a2/cells/A2_rel14_p2.final.json` | `b8c837984e4a8c54ebe074564d8fd354a00badbead59d77fce241e0150c0909a` |
| `results-second-edge-a2/cells/A2_rel14_p2.json` | `40406f0bd17abe777416010e8e6ab1738f14fd58b5896db92feef7e01ce9e099` |
| `results-second-edge-a2/cells/A2_rel14_p3.final.json` | `9f34e34396bda558143f26e651d73923b6753939ca9aa227f0dfe20d61c5257d` |
| `results-second-edge-a2/cells/A2_rel14_p3.json` | `78432affba6927bf1809d52a3f41e1ddc28a9bf5c5b7d973d5f41a2f2a29f1b5` |
| `results-second-edge-a2/cells/A2_rel14_p4.final.json` | `a154e7eeed2d1cbbdf4c88b1cfc51715762ec2edd2e927aba1e4717bea3d1a8a` |
| `results-second-edge-a2/cells/A2_rel14_p4.json` | `11a259e1f9e37678e67a01e1bf258d3da7bdf1372f22a14a25d47924c289da77` |
| `results-second-edge-a2/cells/A2_rel14_p5.final.json` | `21f61910489c9883105e1c388943d11decd02240c24a7cc68d86a9789a72a807` |
| `results-second-edge-a2/cells/A2_rel14_p5.json` | `cd8853593ef1baf15854bf21bebf5183625ea93dd88f80b47c267397db715f5c` |
| `results-second-edge-a2/cells/A2_rel14_paper.final.json` | `9d74e8e98e7d4f67538986d890191cce7679e302d0d713cae045cacbf6e4c157` |
| `results-second-edge-a2/cells/A2_rel14_paper.json` | `44e6438da01a2426932195d5fc6c7a0df96127a63fe0066d03e774e5776cba88` |
| `results-second-edge-a2/cells/A2_rel14_ship.final.json` | `59fbbecd7d432184d2ace9e6953138d9d68fd18d3b4e96d510164f5cf849749b` |
| `results-second-edge-a2/cells/A2_rel14_ship.json` | `ae93ce051261a6a7bdb984a49512b651a80872fd0b04d04852e376c684b5fd9c` |
| `results-second-edge-a2/checks/A_abs514_p1.final.json` | `15f5c7722acdbb45e9deca85ae86f7f626526bbc9d18937c434fae5fee58c3a7` |
| `results-second-edge-a2/checks/A_abs514_p1.json` | `1ca181d9de48060e7b3ffb324086888d77865bf1a05d8f08fef6427ee00cf71e` |
| `results-second-edge-a2/checks/x0_check_summary.json` | `778d1c68a87506fbd0ec25165c2fe44854c41ea3c776213ba0c1f701d334e980` |
| `results-second-edge-a2/checks/x0_fidelity_comparison.json` | `b739793b8b9f6e1553c114b9093cf928d758e958dadd528db6a33a30ede6c06a` |
| `results-second-edge-a2/comparison.csv` | `3e6a6a41bb2960aeb7e19fc4713be48f6a5f9af957ad4e5085f549b154bb92c6` |
| `results-second-edge-a2/counts_by_width_target.csv` | `222e7121c60947f52dffacbfdb8d725b4713af4ec597ab322c245c740ae05c90` |
| `results-second-edge-a2/logs/analysis.log` | `1eaad2b9431334ef2f070e59c60d6171adaa763b0b8f14c42b2c04a474860932` |
| `results-second-edge-a2/logs/chunk1_summary.json` | `ee199cdc20d0afc6873e74938d9f5767fdba11f61fcb2233056b2ec5050fc35c` |
| `results-second-edge-a2/logs/chunk2_summary.json` | `98600ae4cd29093638f4b7e9b356f020c96a96af9587b310f022c0e408394cc1` |
| `results-second-edge-a2/logs/chunk3a_summary.json` | `66ebf8c6054794e135bdcfb41b60891f6d788937055776500ee1f88e8a876768` |
| `results-second-edge-a2/logs/chunk3b_summary.json` | `b40172771fadb8932ed64b31640056513b7b770a71ff8344e54d41658ce350bd` |
| `results-second-edge-a2/logs/chunk3c_summary.json` | `fbdc935ace63fe6fda9867b0e5e57f36ddb35dd285869bcb72c0cfeb6f34149b` |
| `results-second-edge-a2/logs/probe_a2_abs514_p0.log` | `38c89477dd6faa444ec266dbc5cbb56ef48c8d9d91981f1da7f91702fd2715a3` |
| `results-second-edge-a2/logs/probe_x0_abs514_p1.log` | `bab2310083e230274ba820791357b8500e43e8994a8488fe027657aa573c1d3d` |
| `results-second-edge-a2/logs/run_chunk1.log` | `9113e8fc7f7345423e4c4d52fe5e8a21a6b5053cbbbc7ba6cb03a14b157291d0` |
| `results-second-edge-a2/logs/run_chunk2.log` | `f6e87e525b19b56b839682113ede403b90fb550ff16d64469ec3c168eb90bc44` |
| `results-second-edge-a2/logs/run_chunk3a_p5.log` | `0a0be685264e2bbeefec2cee8abd5bbc6525da39315def618343826c3a4b7abb` |
| `results-second-edge-a2/logs/run_chunk3b_p4.log` | `891644febe2b687216926be8036544ae64287e24764f0b25138f8382b40a2dad` |
| `results-second-edge-a2/logs/run_chunk3c_ship_paper.log` | `d716e0468d54a7c66a16d70f58396e35809f27aa6213e470c1e23c57f27c2f63` |
| `results-second-edge-a2/logs/verify_witnesses.log` | `8a2e265c2fb7d516ab9eaa2b3babd60b222385245e2d51ce939c38f8c048e4c7` |
| `results-second-edge-a2/README.md` | `708d7df3961cd930bb1d1911de178b5609078ebaf1b3ff2443253a5f5bd92a87` |
| `results-second-edge-a2/second_edge_a2_solver.py` | `f181b939d3a7091ecb0bb222806fb54d07a6174a8c998e4af23c58d41b13b4b7` |
| `results-second-edge-a2/second-edge-a2-counts.csv` | `f0ddede7107c6cd815f8b5b0ba0a0f557645d4041f89e78512a6060f46b6dcbc` |
| `results-second-edge-a2/second-edge-a2-counts.json` | `0b98c5d33f3cf068bb8eea3182c42e08cc1e831cedd222af81369a2031e7e792` |
| `results-second-edge-a2/verify_witnesses.py` | `a15f2c7af0eda7cdb7ce01cd8373b82c47a4639ca5bd1d2a32f9b6deaea4b381` |

**`recertification/` and `analysis/`**

| File | SHA-256 |
| --- | --- |
| `recertification/cell_5.json` | `3c8144ff86d6ad2e32fed8a43323e9903030cf6e2d587824d48ffdbf344795d7` |
| `analysis/analyze_second_instance_f1.py` | `0b20257d74ee3c1ea248fc1eaa55c703a846228b0303c8cfde1943a484525c07` |
| `analysis/analyze_tail_instance.py` | `d977073baf6dd9576c8c1f45e7084c1cd5a8540328fd9659317fd4d518c69592` |
| `analysis/f1_reference.py` | `711d02e580d7b4fc5d258f21304ea02d37c2fff6bb1b6ba6c208d18da914e56e` |
| `analysis/ladder_analysis.py` | `b10917ccd8bf96b6810279d185fdf964b7e201408d46220bd019b17b3427a3a0` |
| `analysis/second_instance_f1_sweep.py` | `7f682420e74bd6cd85fb1cfc57b7b235edeebc258c493ef24562ca184d471ff4` |
| `analysis/tail_instance_sweep.py` | `a966005882646cf81f25b33d3bf200bd25d819c0fdc6347ec9f1db11955f8e76` |
| `analysis/validate_tail_instance.py` | `06a224e7ef82a6672e15f0131a388998c4a0d0f77be1836c45a1bca3c7f429bd` |

## Not deposited

`recount_cheb.py` is absent. No record names it by digest, and
`records/chebyshev_recount.json` does not name the script that produced it;
it is also outside the measurement tree beside the records. `ladder_analysis.py`
is deposited under `analysis/`, but no record or record metadata names it by
digest either: `recertify_ladder.py` calls it by path and records only its own
digest, so the digest above cannot be checked against a record.

## The 60-digit evaluation control

`control-60dig/` is the evaluation-precision control of the corrected BRASIL
ladder. Five deep cells are re-solved at 60-digit evaluation — the paper width
at `1e-15`, the ship width at `1e-15`, `2^5` at `1e-15` and at `1e-14`, and
`W=1` at `1e-15` — with the record's solver, settings, denominator cap,
2,001-node certifying grid, target semantics and reference, changing one
thing: a fit is accepted at 60-digit evaluation with its coefficients held at
60 significant digits. Each cell first reproduces the record digit for digit,
then walks the counts below the record's; a count is an existence bound under
this solver, in the same sense as the records it controls.

`control_60dig_brasil.py` is the driver; `cell_*.json` are the per-cell
records, `summary.json` the assembled record, `freeline_effect.json` the
free-line test's effect, and the `_run_*.log` and `_verify_*.log` files the
run and verification transcripts. `README.md` there gives the commands, the
reproduction checks and the measured counts — 15, 16, 18, 17 and 8 against
the record's 18, 18, 25, 18 and 8 — and `verify_cells.py` and
`verify_freeline.py` reproduce them.

## Requirements

Python 3.10 or newer with `mpmath`. The portable `count_sweep.py` takes the
library from `BOYS_SOURCE` or a `boys` directory beside this repository. The
recorded scripts hardcode the measurement machine's library paths, and several
read the measurement tree's `results*/` layout; "Harness revisions" names their
options, and the stability command above passes the deposited records with
`--source`. Do not run the scripts with `python -I`.
