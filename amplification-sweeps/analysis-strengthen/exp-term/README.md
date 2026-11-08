# Exponential-term reruns

The records of the exponential-term reruns of the strengthening, surrogate-control and a-priori
experiments under `amplification-sweeps/`. They extend the pre-exponential runs deposited beside
them and are the source of the per-order scores and the a-priori variant results that use the
exponential term.

## Contents

- `apriori/` — the a-priori bound, coefficient, dominance and exact-term scripts and their CSV
  outputs.
- `surrogate/` — the surrogate-control rerun.
- `out/` — the run logs and the generated records the analyses read: `analyze-exp.log`,
  `apriori-*-exp.log`, `intervals-exp.log`, `make-*-exp.log`, `score-*-exp.log`, and the
  `v5.2.0-exp0/` run directory.
- `../out/first/` and `../out/second/` — the new model's prediction, ladder and per-cell records
  (`prediction.csv`, `scan_all.csv`, `scan.csv`, `cell_scan.csv`, `per_order_ratio.csv`). They are
  the inputs `summarize.py` reads.
- `canon_context.py`, `check_module.py`, `dump_stale.py` — the context and the checks the rerun
  scripts share.

`out/v5.2.0-exp0/` is the previous model's record. `summarize.py` reads it only as the comparison
baseline for the row-by-row changes; it is not a current record. The current model's records are
`../out/first/` and `../out/second/`. Running `python summarize.py` from this directory reproduces
`out/summarize.log`.

## Provenance

The reruns were produced in the analysis tree and are deposited here verbatim, except for
byte-compiled caches and paper-editing helpers. The scripts read the cells deposited under
`amplification-sweeps/amplification-form/` and `amplification-sweeps/product-amplification/`; the
logs name the commands and the versions they ran under.
