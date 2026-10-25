# Boys-function measurement harnesses

The scripts behind four measurement campaigns on the Boys function kernel, kept so that the
measurements can be re-run and checked. The measured records for each campaign are in the `data/`
and `results/` directories beside its scripts.

The kernel is the Boys library at <https://github.com/myamlak/boys>.

## Layout

- `amplification-sweeps/` — the per-order accuracy sweeps, two C++ drivers and their analysis.
- `coefficient-count-sweep/` — the stored-coefficient counts, with the stability pass and its
  cold ladder.
- `composition-audit/` — the full-factorial audit of composed accuracy figures.
- `accuracy-claim-audit/` — the documentation-claim audit against two external libraries.

Each directory has a README with its build and run commands.

## Requirements

- A C++23 compiler and CMake 3.20 or newer for the C++ harnesses.
- The Boys source, checked out beside this repository or named by `BOYS_SOURCE` for the Python
  scripts and `-DBOYS_SOURCE_DIR` for CMake.
- Python 3.10 or newer with `mpmath`. The claim audit additionally needs MPFR 4.2.2 and the
  documentation of the two audited libraries.
