# Stored-coefficient counts

`count_sweep.py` measures the smallest stored coefficient count each fit
family needs on `[x0, x0+W)` at a target. It imports the library's generator
rather than reimplementing the fits. The deposited `count_sweep.py` is
portable: it takes the library from `BOYS_SOURCE` or a `boys` directory beside
this repository. It is not the harness the records identify; see
"Recorded harness" below.

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

## Recorded harness

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
harness's `all` mode also tries degrees 0, 1 and 2; the recount records four
cells with `new_m` 2, all in mode `all`. The portable harness also cannot
reproduce the ladder's rational cells with `k` above the generator's
`RAT_K_MAX` of 8 (the ladder records `rat_k_max: 12`). The records that name
a harness name the recorded revision, not this one.

The deposited `stability_save_verify.py` is the revision named by
`records/ladder/rational_control_symmetric.json`
(`meta.stability_script_sha256`), SHA-256
`51c4ec838819967f13560d1a68804c89cbf54c648f0effd8ea0be5a065e6483d`.
Five drivers import `count_sweep` by name: `cold_ladder.py`,
`stability_save_verify.py`, `recertify_ladder.py`,
`rational_control_symmetric.py` and `rational_unrounded_control.py`. To run
them as recorded, place `count_sweep_recorded.py` beside them as
`count_sweep.py`. `rational_control_symmetric.py` also calls
`cs.cheb_degrees`, which only the recorded harness defines.

The recorded harness's default library root is the measurement machine's.
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

The file states the library revision (`boys.version` 5.1.1, `describe`
v5.1.0-1-g3b187d5) and the digests of the generator it was measured with
(`11185fba03abe74d01b30f093818f48972eb8d88413d446b505198b869115180`), the
coefficient header
(`8ce5e17157a4be5f621d695f377094d212cb839732e278d47f23aa225b4d9204`), the
five records it was re-solved from (`paper_records`), and the generator the
earlier counts came from
(`ccc8b011f4170ce0b530bc4781b09670f42839fa18caaa38ea7d3a9ea9ba9bb4`).
Neither generator is in this deposit; the earlier one is the Boys
repository's `tools/gen_boys_coefficients.py`. The five `paper_records`
digests are of the measurement-tree originals. The deposited copies parse to
the same content but are re-serialised, so their bytes hash differently.

The five records hold 135 rows: 40 in `records/core/sweep_all.json`, 20 in
`records/core/sweep_rel_1e-13.json`, 20 in
`records/core/sweep_abs_5e-14.json`, 15 in
`records/core/sweep_rel_1e-14.json`, and 40 in
`records/ladder/sweep_all.json`. The deposit does not re-run them and does not
claim that the deposited harness reproduces them.

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

## Not deposited

`ladder_analysis.py` and `recount_cheb.py` are absent. No record or record
metadata names either file by digest: `recertify_ladder.py` calls
`ladder_analysis.py` by path but records only its own digest, and
`records/chebyshev_recount.json` does not name the script that produced it.
Their digests cannot be checked against the records, so they are not
deposited. `recount_cheb.py` is also outside the measurement tree beside the
records.

## Requirements

Python 3.10 or newer with `mpmath`. The portable `count_sweep.py` takes the
library from `BOYS_SOURCE` or a `boys` directory beside this repository. The
recorded scripts hardcode the measurement machine's library paths, and several
read the measurement tree's `results*/` layout; "Recorded harness" names their
options, and the stability command above passes the deposited records with
`--source`. Do not run the scripts with `python -I`.
