# Boys-function measurement harnesses

The scripts behind four measurement campaigns on the Boys function kernel, kept so that the
measurements can be re-run and checked. The measured records accompany the corresponding
publications.

The kernel is the Boys library at <https://github.com/myamlak/boys>.

## Layout

- `amplification-sweeps/` — the per-order accuracy sweeps. `sweep1.cpp` and `sweep2.cpp` build
  against the Boys source and emit the cells; `analyze.py` turns the cells into the per-order,
  band and prediction tables. `common/` holds the shared provenance, reference and kernel-replica
  helpers.
- `coefficient-count-sweep/` — the stored-coefficient counts. `count_sweep.py` imports the Boys
  generator from the library's `tools/gen_boys_coefficients.py` and measures the smallest stored
  count each fit family needs; `stability_save_verify.py` and `cold_ladder.py` are the stability
  pass and its follow-up.
- `composition-audit/` — the full-factorial audit of how a library composes its per-axis accuracy
  figures. `measure_interaction.cpp` measures the cube against the Boys source; the Python
  scripts generate the reference, summarise the cells, compare configurations and run the
  self-checks. `CMakeLists.txt` shows the build; pass `-DBOYS_SOURCE_DIR=<boys source>`.
- `accuracy-claim-audit/` — the documentation-claim audit against two external libraries. Its
  README lists the dependencies of each audit and the commands that run it.

## Requirements

- C++23 compiler and CMake 3.20 or newer for the C++ harnesses.
- The Boys source, checked out beside the harness or passed through `BOYS_SOURCE_DIR`.
- Python 3.10 or newer, with `mpmath` for the coefficient-count sweep and the reference
  generator. The claim audit additionally needs MPFR 4.2.2 and the documentation of the two
  audited libraries.
