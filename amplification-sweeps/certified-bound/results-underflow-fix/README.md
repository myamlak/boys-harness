# Repaired certificate run, 2026-11-06

This directory is the full rerun of the certified-bound record after the
underflow-radius repair on branch `fix/cert-underflow`. It supersedes
`../results/` for the repaired model; the 2026-11-03 record in `../results/`
is preserved unchanged.

Everything recorded is identical to the old record, file for file:

| file | old vs new |
| --- | --- |
| `summary.json` | byte-identical |
| `cells_cert.csv` | byte-identical (55,176 rows) |
| `band_cert.csv` | byte-identical |
| `model_cells.csv` | byte-identical |

The certificate result is 0 violations at 55,176 entry-cells, worst
measured-to-certified ratio 0.71875 at (A1, x=0, n=24, ship), 264/264
population band intervals inside, 0 of 132 comparisons and 0 of 65
commitments certified. The per-order medians are certified charge over the
loser's truncation 83.27 and certified combined charge over the gap 85.23.

The underflow repair changes no recorded radius: every total radius on the
grid is at least 1.7e-17, so the absolute allowance 2^-1071 carried at the
underflowed product propagates into much larger terms and rounds away.
`underflow_test.log` and `underflow_test.csv` record the operation-level
repair at x = 2^-1074: the exact product 7.6010099360191773e-326, the
computed zero, the radius 4e-323 and every affected downward step.

`commands.md` lists the exact commands and exit codes, including the one
diagnostic not rerun here (`tiers.py`; its only model-dependent input is
byte-identical). `provenance.json` digests the kernel, the sweep, the paper
model record, the deposit, the exponential probe files, the code and every
file in this directory.

`harness_cells.csv` differs from the old file only in the trailing digits of
the long `mpmath` decimal expansions of `need_paper` and `need_cert`
(an mpmath version difference); the medians quoted by the record are the same
to every printed digit. `dps_check.log` likewise differs only in the last
printed digits of a 1e-65 quantity.
