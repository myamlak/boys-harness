# Analysis and tables

The scripts that generate the tables from the deposited records, and the
tables themselves. `make_tables.py`, `make_endtoend.py`, `make_frontier.py`,
`make_prediction_record.py` and `merge_per_order.py` produce the printed
tables; `make_heldout_score.py`, `make_seed_fit_ratio.py`,
`make_seedfit_score.py` and `sweeps/amplification_form/seed_fit_ratio.py`
produce the five experiment tables.

## What is here

    make_tables.py             band-level, band-ratios, gain, growth, recovered-b2
    merge_per_order.py         per-order.tex from the four per-order intermediates
    make_endtoend.py           endtoend.tex
    make_frontier.py           frontier.tex and figures/per-order-ratio.dat
    make_prediction_record.py  prediction-record.tex
    make_heldout_score.py      heldout-score.tex and heldout-commitments.tex
    make_seed_fit_ratio.py     seed-fit-ratio.tex
    make_seedfit_score.py      seedfit-score.tex and seedfit-misses.tex
    sweeps/amplification_form/analyze.py
                               the first sweep's prediction, per-order and band records
    sweeps/amplification_form/seed_fit_ratio.py
                               seed_fit_ratio.csv, the input to seed-fit-ratio.tex
    sweeps/product_amplification/analyze.py
                               the second sweep's per-order, band, growth and gain records

`tables/` holds the tables; fourteen are printed. The other four,
`per-order-A1-{abs,rel}.tex` and `per-order-A2-{abs,rel}.tex`, are written by
`make_tables.py` and consumed by `merge_per_order.py`, and are carried here
because the generation chain writes them into the same directory.

## The v5.2.0 check

The table step recorded in `../v5.2.0-check/tables-after-deposited.sha256` was
run with the scripts of that check's own date. The scripts here are the ones
that generate the tables carried here. Three tables (`endtoend.tex`, `frontier.tex`,
`prediction-record.tex`) changed after the check, when `analyze.py` gained the
exponential term in the accumulated-arithmetic model. The check's own record
is not altered.

## Regenerating from the deposit

The scripts resolve their inputs relative to their own directory, at the paths
they occupy in the measurement tree, so assemble a scratch tree and run them there.
From this directory:

    scratch/
      *.py                                  the eight table scripts, from here
      sweeps/amplification_form/analyze.py  from here
      sweeps/amplification_form/seed_fit_ratio.py   from here
      sweeps/product_amplification/analyze.py       from here
      sweeps/common/                        from ../v5.1.1-9353/sweeps/common/
      sweeps/amplification_form/{data,results}
                                            from ../v5.1.1-9353/sweeps/amplification_form/
      sweeps/product_amplification/{data,results}
                                            from ../v5.1.1-9353/sweeps/product_amplification/
      sweeps/seed-spliced/                  from ../../seed-spliced/
      sweeps/held-out/                      from ../../held-out/
      tables/  figures/                     empty

Decompress `sweeps/amplification_form/results/paths.csv.xz` in place before
running anything; the other records in the v5.1.1 tree are plain CSV. (The four
xz archives named in the top-level README belong to the top-level v5.0.0
records, not to this tree.) Then, from `scratch/`:

    python sweeps/amplification_form/analyze.py --tree sweeps/amplification_form
    python sweeps/amplification_form/seed_fit_ratio.py
    python sweeps/product_amplification/analyze.py --tree sweeps/product_amplification
    BOYS_SWEEP_TREE=sweeps/product_amplification python make_tables.py
    python merge_per_order.py
    python make_endtoend.py
    python make_frontier.py
    python make_prediction_record.py
    python sweeps/seed-spliced/analyze.py --tree sweeps/seed-spliced
    python sweeps/held-out/analyze.py --tree sweeps/held-out
    python make_seedfit_score.py
    python make_seed_fit_ratio.py
    python make_heldout_score.py

On PowerShell, set `$env:BOYS_SWEEP_TREE` before `make_tables.py`. The held-out
analysis and `make_heldout_score.py` read the first sweep's
`results/prediction_record.csv` and `tables/prediction-record.tex`, so run the
first sweep's analysis and `make_prediction_record.py` before the held-out
step. Run without `python -I`: isolated mode hides the user-site `mpmath`.

All eighteen `tables/*.tex` and `figures/per-order-ratio.dat` then reproduce
byte-for-byte, as do every regenerated analysis record except the
`summary.txt` files, whose second line embeds the absolute reference path. The
tables and the figure data use CRLF; the record CSVs use LF. On Windows the
regenerated files are byte-identical with no folding; folding CRLF to LF on
another platform preserves the content and is the only adjustment needed.
