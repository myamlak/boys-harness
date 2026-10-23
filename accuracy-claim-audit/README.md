# Accuracy-claim audit

Scripts that audit documented accuracy claims in two external libraries:
Boost.Math 1.92.0 and crlibm 1.0beta5. Each documented statement is classified
as a row, a domain and a bound a program can evaluate, and the row is measured
against an independent MPFR reference.

## Contents

- `boost_accuracy_rows/` — the Boost.Math audit. `src/gate.cpp` applies the row
  table in `claims/rows.csv` and reports one verdict per row. `src/probe.cpp`
  checks a few library values directly. `claims/claims_audit.tsv` classifies
  the documented statements in the nine audited documentation sections. The
  Python tools extract the rendered documentation, the machine-checked bounds
  from the test suite, and the published test data; `independent_check.py`
  compares the gate's reference values with mpmath.
- `claim-audit-crlibm/work/` — the crlibm audit. `ref.py` is the MPFR Ziv-loop
  reference and `validate_ref.py` checks it against an exact Fraction oracle
  and against crlibm's own test vectors. `run_crlibm.py` sends queries to
  `driver.c`, which is built against crlibm. `m1_testdata.py`,
  `m2_claims.py` and `m3_defects.py` are the three measurements;
  `verify_findings.py` checks the two measured defects against mpmath.
  `claims_registry.tsv` is the claim classification.

## Headline numbers

`headline_numbers.py` prints the five counts that follow from the two claim
registries, from any checkout of this directory:

```sh
python headline_numbers.py
```

The three false claims need the captured cell summary
`records/summarize_cells.txt`, which is supplied with the corresponding paper
and is not in this repository; the script says so when it is absent.

## Dependencies

- Python 3.10 or newer with gmpy2 2.3.1 linked to MPFR 4.2.2 and GMP 6.3.0,
  and mpmath 1.4.1. Check the versions with:

  ```sh
  python -c "import gmpy2, mpmath; print(gmpy2.version(), gmpy2.mpfr_version(), mpmath.__version__)"
  ```

- A C++17 compiler. `boost_accuracy_rows/build.bat` uses MSVC 14.5 (Visual
  Studio 18) and expects the Boost tree at
  `boost_accuracy_rows/boost_1_92_0`.
- Boost.Math 1.92.0, including its source, test data and rendered Math
  documentation under `libs/math/doc/html/math_toolkit`.
- crlibm 1.0beta5, built on x86-64 Linux in the non-double-extended
  configuration. The audit uses WSL Ubuntu, autotools, gcc 15.2 and
  `--enable-sse2 --disable-mpfr`.
- A Linux or WSL environment for the crlibm driver.

Set the local paths at the top of these scripts before running: `VSDIR` in
`build.bat`, `DOC` in `extract_claims.py`, `TD` in `parse_testdata.py`, and
`datadir` in `validate_ref.py`. `run_crlibm.py` invokes `wsl.exe` and expects the driver at
`~/crdrv`; on a native Linux host, use `["bash", "-lc", "~/crdrv"]` instead.

## Running the Boost.Math audit

Place Boost 1.92.0 at `boost_accuracy_rows/boost_1_92_0`, then from
`boost_accuracy_rows`:

```bat
build.bat gate
python tools\emit_data.py boost_1_92_0\libs\math\test data
build\gate.exe --rows=claims\rows.csv
```

The gate exits non-zero when a judged row is unmet; that is a measured result.
Add `--cells` for the per-cell record, and
`--refout=claims\reference_dump.csv` to write the reference values, which
`python tools\independent_check.py claims\reference_dump.csv` then checks
against mpmath.

To rebuild the claim census and the machine-checked bounds:

```bat
if not exist extracted_accuracy mkdir extracted_accuracy
python tools\accuracy_text.py boost_1_92_0\libs\math\doc\html\math_toolkit extracted_accuracy > extracted_accuracy\INDEX.tsv
python tools\expected_results.py boost_1_92_0\libs\math\test claims\expected_bounds.csv
python tools\survey.py domains boost_1_92_0\libs\math\test
```

`census.py` takes the extracted section files and prints the count of cells
that report this library. `summarize_cells.py` reports per-row judged and
degenerate cell counts from a gate run with `--cells`.

## Running the crlibm audit

Build crlibm from its source root in WSL or Linux, with `driver.c` from
`claim-audit-crlibm/work/` copied into that root:

```sh
autoreconf -i
./configure --enable-sse2 --disable-mpfr
make -C scs_lib
make libcrlibm.a
gcc -O2 -std=c99 -msse2 -mfpmath=sse -I. -o ~/crdrv driver.c libcrlibm.a -lm
```

Then, from `claim-audit-crlibm/work` with the paths above set:

```sh
python extract_claims.py
python parse_testdata.py
python validate_ref.py
python m1_testdata.py
MQ=20000 MP=10 python m2_claims.py
python m3_defects.py
python verify_findings.py
```

`m1_testdata.py` writes `m1_results.tsv` and `m2_claims.py` writes
`m2_results.json`; the remaining scripts print their results. Running
`validate_ref.py` with function names restricts its second part to those
functions.
