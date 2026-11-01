# Audit replay

This directory is the end-to-end replay of the reconstructed audit package in
`../` (Boost.Math 1.92.0 and crlibm 1.0beta5). It was run on 26–27 October
2026 from the pinned sources, and it is the second derivation of the counts
that the supplement previously read from the session record `counts.json`.

## Sources and provenance

- Boost.Math 1.92.0 archive, SHA-256
  `c4a3b310ddd2472416e091067166b0713be97c63f38c212c484ada022fd296ce`
  (235,384,471 bytes). The tree's `boost/version.hpp` is
  `cf992d9d2c4f2294b69d8819e334927644a48b978bfea3c525e3f9c6dbc7cb9f`.
- crlibm 1.0beta5, commit `eb3063791aa75bc9705b49283bf14250465220a7`,
  archive SHA-256 `79d626f9515fccdbc724e8b74cbf848331317b1891ef22aaaf6a195e20bbbfe5`
  (2,432,556 bytes). The WSL working tree carries CRLF line endings; after
  removing the carriage returns its `csh_fast.h` hashes to
  `0917e43ed1c539176c3a8e09e55980320ee3a705ee0f9221970ae873417f3de0`,
  identical to the archive, so the only difference is the line ending.
- crlibm driver `~/crdrv`, SHA-256
  `7ae0fa4587943300709ec0efcbab295fe61bda38ac600955b3aba2f95240dee6`,
  built 2026-09-24 against the tree above in the **double-extended**
  configuration. The session registry's description says
  "non-double-extended"; `config.log` has `USE_HARDWARE_DE_TRUE=''` (the
  conditional is true) and `libcrlibm.a` contains `log-de.o` with no
  `log-td.o`, so the build is double-extended and the registry's description
  is wrong. The counts are measurements of the double-extended library.
- Python 3.14.7, gmpy2 2.3.1 linked to MPFR 4.2.2 and MPC 1.4.0, mpmath
  1.4.1. The full provenance block is `logs/provenance.log`.

## Commands and result

Every command, its exit status and its wall time are in `logs/commands.tsv`;
each also has its own log in `logs/`. The measurement commands are:

```
python m1_testdata.py
python m2_claims.py        # M2_PERQUERY=<dir> writes the per-query output
python m3_defects.py
python validate_ref.py     # plus the 5,999-trial check in tools/
python tools/count_claims.py
python headline_numbers.py
```

The replay reproduces the session's counts:

| count | `counts.json` | replay |
| --- | --- | --- |
| m1 crlibm vectors | 63,617 | 63,617 |
| m1 disagreements | 9 | 9, all NaN bit patterns; 0 on finite values |
| m2 queries | 361,286 | 361,286 |
| m2 rounding-mode results | 1,445,144 | 1,445,144 |
| m2 correct-rounding failures | 12 (4 sinh, 8 atanpi) | 12 (4 sinh, 8 atanpi) |
| m2 undecided reference points | 44 | 44 |
| Boost documented statements | 320 (298 cells, 22 prose; 9 sections) | 320 (298, 22; 9) |
| Boost row-shaped, interpreted | 2 | 2 |
| crlibm documented claims | 42 | 42 |
| crlibm strictly row-shaped | 18 | 18 |
| crlibm row-shaped if interpreted accepted | 39 of 42 | 39 of 42 |
| false records | 3 (Boost 2, crlibm 1) | 3 (Boost 2, crlibm 1) |
| Boost test-suite bounds | 499 | 499 |
| crlibm vector assertions | 63,617 over 21 files | 63,617 over 21 files |
| crlibm Gappa obligations | 791 over 46 files | 791 over 46 files |
| crlibm source-comment bounds | 41 | 41 |
| aggregate counts compared | 20 | 20 agree, 0 disagree |
| crlibm own-vector points decided | about 48,000 | 42,692 decided of 42,824 checked (132 undecided), 0 disagreements |
| reference validation | 0 mismatches in 5,999 trials, including 72 overflow-boundary cases | 0 |

The Boost registry's `row_shaped` marker is the interpreted reading: both
entries are the one Airy sentence, and the paper counts no Boost claim
strictly row-shaped. The aggregate of 20 row-shaped claims is therefore 2
Boost interpreted plus 18 crlibm strict.

The m2 failures are the four `sinh` modes at `0x408633ce8fb9f87e` and its
negative mirror, and the eight `atanpi` modes at and above
\(2^{54} = 1.8014\times10^{16}\), where the bisection in
`logs/atanpi_onset.log` places the first failing double (`0x4350000000000000`,
the last passing one being `0x434fffffffffffff`); `m3_defects.py` reproduces
both, and `cosh` at the same input is the control and is correct. The session
record's "about \(3\times10^{16}\)" is a rounded threshold; the replay's
bisection is the exact onset. The defect persists to the largest finite
double, not only to \(10^{300}\): the large-point check fails at \(10^{17}\),
\(10^{18}\), \(10^{20}\), \(10^{100}\), \(1.0000000000000001\times10^{300}\)
and `DBL_MAX`.

## Files

- `out/m2_perquery.tsv.xz`: 1,445,144 rows, one per query and rounding mode.
  Columns: `function, mode, input_hex, crlibm_output_hex, reference_status,
  reference_output_hex, result`.
- `out/m2_perinput.tsv.xz`: 361,286 rows, one per query, with the four modes
  and the four reference values side by side, the category, and the failure
  count and failing modes.
- `out/m1_vectors.tsv.xz`: 63,617 rows, one per crlibm test vector. Columns:
  `func, mode, input, expected, got, result, expected_nan, got_nan`.
- `out/m1_results.tsv`: the nine disagreements.
- `out/m2_results.json`: the 18-function summary and the twelve failure rows.
- `logs/`: every command log, `commands.tsv`, `provenance.log`, `recount.log`,
  the bisection log `atanpi_onset.log`, and the five repair diffs.
- `recount/counts.json`: the recount re-run on the pinned sources, with all
  20 aggregate counts compared.
- `scripts/`: the repaired audit scripts, `claims_registry.tsv` and
  `driver.c`.

## Repairs made to the reconstructed package

The reconstruction carried two kinds of defect, both recorded as diffs in
`logs/`:

1. Three scripts (`extract_claims.py`, `parse_testdata.py`,
   `validate_ref.py`) pointed at the dead working directories
   `_claim-audit-crlibm/src/crlibm/{docs,tests}`. They now point at the
   pinned crlibm tree.
2. The measurement scripts wrote only summaries: `m1_testdata.py` now writes
   one row per vector, and `m2_claims.py` writes one row per query and one
   per query-and-mode. The measurement and the reference are unchanged; the
   added code is output only. The nine m1 disagreements are confirmed as NaN
   on both sides by an explicit NaN test.

No count was adjusted to match `counts.json`. The replay was written first and
compared afterwards, and the table above is the comparison.

## What is not replayed

The Boost gate's raw `--cells` stream was not saved in the session and is not
regenerated here; its summariser output is `../records/summarize_cells.txt`,
and the two Boost claims it flags are re-derived argument by argument in the
witness rerun deposited at `accuracy-claim-audit/witness-rerun/` in the
harness repository. The Boost-versus-mpmath 2.7 eps cross-check was not
re-run; the witness rerun states that it does not re-run the Boost arm, and
the audit README gives the `independent_check.py` command that produces it.
The Gappa obligations were
counted, not re-run, as the paper says. The crlibm `static` routines'
worst-case rows remain unsettled for the reason the paper gives.
