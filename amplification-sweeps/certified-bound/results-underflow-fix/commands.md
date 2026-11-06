# Commands and exit codes: repaired certificate run, 2026-11-06

Working directory: `amplification-sweeps/certified-bound`.  All Python runs
use `-I` (isolated mode).  Runtime: Python 3.14.7, mpmath 1.4.1, Windows 11
(recorded in `provenance.json`).  Kernel: a boys checkout at
`/c/Users/myaml/dev/_boys-tags`, tag `v5.2.5`, commit
`8b7a289dfce7558c8498e89c17c05f20fa7bb7c8`; the checkout is selected with
`--kernel` (the default is `$BOYS_ROOT`, else `<repository>/../boys`).  The
sweep is the deposited `../amplification-form`, the model record is
`inputs/model-record`, the deposit is the in-repo `../apriori-bound`, and the
exponential bits are `exp/exp_bits.txt`.  Every output is written beneath
`--out results-underflow-fix`.

| # | Command | Exit | Key output |
| --- | --- | --- | --- |
| 1 | `python -I code/certify.py --kernel /c/Users/myaml/dev/_boys-tags --sweep ../amplification-form --model-record inputs/model-record --exp-bits exp/exp_bits.txt --out results-underflow-fix > results-underflow-fix/certify.log 2>&1` | 0 | 55,176 entry-cells; validation violations 0, worst \|r\|/C 0.71875 at (A1, 0.0, 24, ship), population bands 264/264; 132/132 band comparisons overlap, 65/65 commitments overlap, 67/67 abstentions overlap |
| 2 | `python -I explore/harnesscells.py --kernel /c/Users/myaml/dev/_boys-tags --sweep ../amplification-form --model-record inputs/model-record --deposit ../apriori-bound --exp-bits exp/exp_bits.txt --out results-underflow-fix > results-underflow-fix/harness.log 2>&1` | 0 | 65 commitments; exact-truncation ordering 49 agreeing / 16 opposite; certified charge/loser truncation median 83.27064776166068 (paper charge 50.22395288436263); over the 49: certified combined charge/gap median 85.225820155890797 (paper convention 57.393556785409458); loser delivered/truncation median 9.57302184326235 over all 65 and 6.273856942285998 over the 49; charge violations: certified 0, paper 5; exact-CSV cross-check 65 cells max 1.1671903103190006e-24 |
| 3 | `python -I code/replay.py --kernel /c/Users/myaml/dev/_boys-tags --sweep ../amplification-form --exp-bits exp/exp_bits.txt --both-downward > results-underflow-fix/replay.log 2>&1` | 0 | 27,588 ship and 27,588 product values replay bit for bit with separate downward steps (0 mismatches); the fused downward step mismatches 208 ship values, first at (A1, 0.029916890803630927, 0), 1 ulp |
| 4 | `python -I code/oracle_check.py --kernel /c/Users/myaml/dev/_boys-tags --sweep ../amplification-form --deposit ../apriori-bound --exp-bits exp/exp_bits.txt --out results-underflow-fix > results-underflow-fix/oracle_check.log 2>&1` | 0 | 60 cells against the deposit's identity oracle; max \|diff\| 7.9872022e-24 at A1, x=4.0666245819265257e-4, n=24, product |
| 5 | `python -I explore/dps_check.py --kernel /c/Users/myaml/dev/_boys-tags --sweep ../amplification-form --deposit ../apriori-bound --exp-bits exp/exp_bits.txt > results-underflow-fix/dps_check.log 2>&1` | 0 | at that cell the route truncation is stable at dps 70/100/140 (5.57277020162e-19), agrees with an independent hyp1f1 F_n to -4.4e-65, and the oracle converges to the route value at dps 140 (the oracle's dps-70 quadrature is the inaccurate side) |
| 6 | `python -I explore/paper93.py --kernel /c/Users/myaml/dev/_boys-tags --sweep ../amplification-form --model-record inputs/model-record --deposit ../apriori-bound --exp-bits exp/exp_bits.txt > results-underflow-fix/paper93.log 2>&1` | 0 | model band-max cells differ from per_order.csv at 102/132 comparisons; 50 agreeing / 15 opposite; loser \|r\|/t median 3.2679993240833323 at the model cells; paper charge/gap median 57.30090179062915, certified 84.96085481278541 |
| 7 | `python -I code/test_underflow.py --kernel /c/Users/myaml/dev/_boys-tags --exp-bits exp/exp_bits.txt --out results-underflow-fix > results-underflow-fix/underflow_test.log 2>&1` | 0 | at x=2**-1074 the exact product x*F_32 is 7.6010099360191773e-326, the computed product is 0.0, the repaired radius is 4e-323; all 32 downward product nodes covered; the whole chain's worst \|model - exact route\| is 1.85037e-17 inside the certified radii; the normal control x=1e-12 gets no allowance |
| 8 | `python -I exp/exp_check.py --exp-bits exp/exp_bits.txt > results-underflow-fix/exp_check.log 2>&1` | 0 | 9,541 arguments; worst \|delta\|/ref 0.999662996081 * 2^-53 at x=11.084208097312864; 2,585 above 0.5 unit, 0 above 1 unit; Python `math.exp` mismatches 0; half-scaling mismatches 0; conditional eta_E <= 2^-52 supported |
| 9 | `python -m py_compile code/*.py explore/*.py exp/*.py` | 0 | all certificate Python compiles |
| 10 | `python -I code/make_provenance.py --kernel /c/Users/myaml/dev/_boys-tags --sweep ../amplification-form --model-record inputs/model-record --deposit ../apriori-bound --exp-bits exp/exp_bits.txt --manuscript /c/Users/myaml/dev/papers-a-v5/papers/boys-crossover/supplement/main.tex --out results-underflow-fix` | 0 | writes `provenance.json` covering the kernel, sweep, model record, deposit, exponential files, code and every file in this directory |

The repaired run's recorded outputs are byte-identical to the 2026-11-03
record: `summary.json`, `cells_cert.csv`, `band_cert.csv` and
`model_cells.csv` compare equal with `cmp`.  `harness_cells.csv` differs only
in the trailing digits of the long mpmath decimal expansions of `need_paper`
and `need_cert`; `dps_check.log` differs only in the last printed digits of a
1e-65 quantity.  The medians quoted by the record are unchanged.

`code/tiers.py` (the seven a priori variants) was not rerun in this session.
Its only model-dependent input is the per-order cell radius `C`, and the full
`cells_cert.csv` that carries every `C` is byte-identical to the old record;
the nominal values are replay-identical (row 3), and the sweep and model
records are the same files.  Its 2026-11-03 result therefore stands: under all
seven variants, with the deposit's own charge and with the certified `C`,
0 of 132 comparisons and 0 of 65 commitments certify.  The script is
portable (row 9 compiles it); it was stopped after 50 minutes when a
host-side job left it at about one tenth of a core.
