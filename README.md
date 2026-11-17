# Boys-function measurement harnesses

This repository is a deposit of measurement campaigns. The scripts behind them
are here with their measured records, so that each record is regenerated from
the deposit, verified against it, or named as a limitation. There is no current
harness beside a historical one: a record measured at an earlier kernel
revision is provenance inside this deposit, and each record's own provenance
names the revision it was measured at.

The kernel is the Boys library at <https://github.com/myamlak/boys>. An earlier
revision named in a provenance file is that record's provenance, not a current
reference.

`REPRODUCTION-MANIFEST.md` is the map: it names, for each campaign study, the
records, the commands that reproduce or verify them, the measured environments
and the limitations. Read it first.

## What regenerates and what verifies

Each campaign README names the commands and, per record, whether they
regenerate it from the deposit or verify it against the deposit. The division
below is a summary; the per-record statement is the campaign README's and the
limitation column names what neither reproduces nor verifies.

| Campaign | Regenerated from the deposit | Verified against the deposit | Limitation |
| --- | --- | --- | --- |
| `amplification-sweeps/` | the two sweeps from the deposited grids, fits, 60-digit references and drivers; the strengthening summary (`analysis-strengthen/exp-term/summarize.py`); the cutoff summary (`amplification-form/cutoff_summary.py`); the interval/enclosure record (`interval-enclosure/`) from its three scripts | the certificate record and its replays; the coefficient-header identity check | the first sweep's 9,353-argument record is an earlier revision's run; four records are xz-compressed and must be decompressed first; the interval/enclosure record reads only deposited records and no boys checkout |
| `coefficient-count-sweep/` | the 60-digit control (`verify_cells.py`, `verify_freeline.py`) and the stability pass (`stability_save_verify.py --source`); the records that name the harness revision, with `count_sweep_recorded.py` placed beside the drivers | the recorded counts and file digests; the portable `count_sweep.py` checks the static counts | the `baryrat` 2.1.2 blocks are copied records, not regenerated (the solver is not deposited); the portable `count_sweep.py` cannot reproduce the ladder rational cells with `k > 8`; some recorded drivers read the measurement tree's `results*/` layout; one ladder cell has no independent re-evaluation |
| `composition-audit/` | the interaction records from the compressed deposited cell files (`signed_interaction.py`, `figure-ratio.json`); the measurement from the driver and the reference generator; the second fractional-design instance (`second-instance/`) from its driver | the self-checks | a C++23 compiler and CMake 3.20 or newer are needed; the second-instance record additionally needs claimgate 0.1.0 (tag `v0.1.0`, an external C++17 package) |
| `accuracy-claim-audit/` | the claim census and the counts from the deposited tools; the crlibm and Boost audits from the pinned sources; the libint2 accuracy comparison (`comparison-libint2/`) from its probe | `headline_numbers.py` and `recount.py` against the deposited counts and registry | Boost.Math 1.92.0 sources and documentation, and crlibm 1.0beta5 (`VERSION`; archive SHA-256 in `witness-rerun/`), are external; the crlibm driver needs Linux or WSL; the libint2 comparison needs libint2 2.13.1 and a Linux or WSL g++ |
| `chemistry-bridge/` | the probe regenerates the dense tensors from the deposited argument streams and the library; the summaries, matrices, traces and the above-B2 statistic (`scripts/above_b2_max.py`) from the deposited records | `compare_arms.py`, `analyze.py`, the probe's own captures | the dense ERI dumps and the benzene argument streams are not deposited; the probe regenerates them (about twenty minutes per benzene arm) |
| `boys-crossover/` | the two route-comparison study input lists (`r4-and-9353/`) from the deposited grids and references | the exact-tag replays at v5.3.1 (`exact-tag-replay/`) and v5.3.2 (`exact-tag-replay-v5.3.2/`): both sweeps and the certificate, byte for byte against the deposited records | the two rebuilt sweep executables are build products and are not deposited (compiler, build line and SHA-256 are); the 291,827,555-byte path record is deposited as a digest only |

## Layout

- `amplification-sweeps/` — the per-order accuracy sweeps, two C++ drivers, their analysis and
  the interval/enclosure measurement.
- `coefficient-count-sweep/` — the stored-coefficient counts, with the stability pass, its cold
  ladder and the 60-digit evaluation control.
- `composition-audit/` — the full-factorial audit of composed accuracy figures and the second
  fractional-design instance.
- `accuracy-claim-audit/` — the documentation-claim audit against two external libraries and
  the libint2 accuracy comparison.
- `boys-crossover/` — the route-comparison study's exact-tag replays at v5.3.1 and v5.3.2 and the two
  reconstructed input lists.

Each directory has a README with its build and run commands.

## Requirements

- A C++23 compiler and CMake 3.20 or newer for the C++ harnesses.
- The Boys source, checked out beside this repository or named by `BOYS_SOURCE` for the Python
  scripts, by `BOYS_TREE` for the two held-out runners, and by `-DBOYS_SOURCE_DIR` for CMake.
- Python 3.10 or newer with `mpmath`. The manifest's environment table names the
  measured stack for each campaign study; the claim audit additionally needs MPFR
  4.2.2 and the documentation of the two audited libraries.
- The libint2 comparison needs a libint2 2.13.1 install tree (`LIBINT_ROOT`), g++ 15.2 on Linux
  or WSL, and the boys `v5.3.2` checkout named by `BOYS_SOURCE`.
- The second-instance record needs claimgate 0.1.0, tag `v0.1.0`, in a clone beside this
  repository or named by `CLAIMGATE_REPO`, with CMake 3.20 or newer and a C++17 compiler.
