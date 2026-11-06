# Commands and exit codes

## The repaired rerun (2026-11-06)

The repaired record was rerun with the portable command set below, from this
directory, with a boys checkout at v5.2.5 (`--kernel`, elsewhere
`$BOYS_ROOT`, default `<repo>/../boys`).  Every output is written beneath the
run directory.  The exact commands, exit codes and outputs of that run are in
`../results-underflow-fix/commands.md`; the 2026-11-03 table further down is
preserved as the historical record of the first run.

    python -I code/certify.py --out results-underflow-fix
    python -I explore/harnesscells.py --out results-underflow-fix
    python -I code/replay.py --both-downward > results-underflow-fix/replay.log 2>&1
    python -I code/oracle_check.py --out results-underflow-fix
    python -I explore/dps_check.py > results-underflow-fix/dps_check.log 2>&1
    python -I explore/paper93.py > results-underflow-fix/paper93.log 2>&1
    python -I code/tiers.py --out results-underflow-fix
    python -I code/test_underflow.py --out results-underflow-fix
    python -I exp/exp_check.py > results-underflow-fix/exp_check.log 2>&1
    python -I code/make_provenance.py --out results-underflow-fix

## The archived run (2026-11-03)

Run on 2026-11-03.  All Python runs use `-I` (isolated mode) and the working
directory shown.  The record's own scripts are under `_a-cert-run`; nothing in
the deposited harness or the paper was modified.  The paths below are
historical: the repaired scripts resolve them with `record_paths.py` and do
not hard-code them.

| # | Command (working directory) | Exit | Key output |
| --- | --- | --- | --- |
| 1 | `python -I code/certify.py > results/certify.log 2>&1` (_a-cert-run) | 0 | 55,176 entry-cells; validation violations 0, worst \|r\|/C 0.71875 at (A1, 0.0, 24, ship), population bands 264/264; 132/132 band comparisons overlap, 65/65 commitments overlap, 67/67 abstentions overlap |
| 2 | `python -I explore/harnesscells.py > results/harness.log 2>&1` (_a-cert-run) | 0 | 65 commitments; exact-truncation ordering 49 agreeing / 16 opposite; certified charge/loser truncation median 83.27 (paper charge 50.22); over the 49: certified combined charge/gap median 85.23 (paper convention 57.39); loser delivered/truncation median 9.573 over all 65 (the paper's 9.6) and 6.274 over the 49; charge violations: certified 0, paper 5; exact-CSV cross-check 65 cells max 1.167e-24 |
| 3 | `python -I code/replay.py /c/Users/myaml/dev/papers-a-v5/papers/boys-crossover/sweeps/amplification_form --both-downward > results/replay.log 2>&1` (_a-cert-run) | 0 | 27,588 ship and 27,588 product values replay bit for bit with separate downward steps (0 mismatches); the fused downward step mismatches 208 ship values |
| 4 | `python -I code/oracle_check.py > results/oracle_check.log 2>&1` (_a-cert-run) | 0 | 60 cells against the deposit's identity oracle; max \|diff\| 7.9872e-24 at A1, x=4.0666245819265257e-4, n=24, product |
| 5 | `python -I explore/dps_check.py > results/dps_check.log 2>&1` (_a-cert-run) | 0 | at that cell the route truncation is stable at dps 70/100/140, agrees with an independent hyp1f1 F_n to 4.4e-65, and the oracle converges to the route value at dps 140 (the oracle's dps-70 quadrature is the inaccurate side) |
| 6 | `python -I explore/paper93.py > results/paper93.log 2>&1` (_a-cert-run) | first 1, rerun 0 | first attempt: `KeyError: 'predicted_ship_at_x'` (read prediction.csv; the at_x columns are in prediction_record.csv); rerun: model band-max cells differ from per_order.csv at 102/132 comparisons; 50 agreeing / 15 opposite; loser \|r\|/t median 3.268 at the model cells (the paper's 9.6 belongs to the per_order cells, row 2); paper charge/gap median 57.30, certified 84.96 |
| 7 | `cut -d' ' -f1 exp_bits.txt > xbits_in.txt` (exp) | 0 | recovers the probe input bits from the deposited probe output |
| 8 | `./exp_probe.exe xbits_in.txt exp_bits_reproduced.txt` (exp) | 0 | rerun of the deposited probe binary |
| 9 | `cmp exp_bits.txt exp_bits_reproduced.txt` (exp) | 0 | byte-identical |
| 10 | `python -I exp_check.py > exp_check.log 2>&1` (exp) | 0 | 9,541 arguments; worst \|delta\|/ref 0.999662996081 * 2^-53 at x=11.084208097312864; 2,585 above 0.5 unit, 0 above 1 unit; Python `math.exp` mismatches 0; half-scaling mismatches 0; conditional eta_E <= 2^-52 supported |
| 11 | `python -I code/tiers.py > results/tiers.log 2>&1` (_a-cert-run) | first 1, rerun 0 | first attempt: `AttributeError: module 'bounds' has no attribute 'certify_magnitude'`; fixed to `bounds.certify` (absolute) and `bounds.certify_relative(..., signed_method="coef")` (relative), matching the deposit's run_all.py verdict; rerun: under all seven variants, both with the deposit's own charge and with the certified C replacing it, 0 of 132 comparisons certified, 132/132 overlap, 0 agreeing, 0 opposite; committed 0/65 certified, 65/65 overlap |
| 12 | `python -I code/make_provenance.py` (_a-cert-run) | 0 | writes `provenance.json` covering every input, script and result file, with digests |

The exponential probe build used
`/nologo /std:c++latest /O2 /EHsc /fp:precise /arch:AVX2` with MSVC
14.51.36231 in the measuring session; `provenance.json` records the compiler,
flags, binary digest and output digest.  The probe run is re-executed above
and reproduces the deposited output byte for byte, so the evidence does not
rest on an unlogged build.

Diagnostic one-liners used to count the commitment composition and the
delivered/truncation ratios are not part of the pipeline; their numbers are
recomputed by `harnesscells.py`, `paper93.py` and `replay.py` and appear in
those logs.
