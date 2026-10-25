# Composition audit

The interaction measurement: how a library composes per-axis accuracy figures
under a joint call. The driver measures every cell of a six-axis, four-lane
design against a 60-digit reference.

## Build

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --target measure-interaction
```

The build uses the `boys` checkout beside this repository by default. Pass
`-DBOYS_SOURCE_DIR=<path>` for another checkout.

## Run

The run is timing-sensitive: use a quiet machine, alone. There are two targets.

Windows (MSVC):

```powershell
powershell -File run_measurement.ps1 -QuietMachine
```

Linux or WSL (GCC), reusing the primary run's reference:

```sh
bash run_measurement.sh --quiet-machine \
    --reference results/reference.csv \
    --arguments-reference results/arguments.csv
```

The driver narrows each nominal argument to its lane's format, emits the
per-lane arguments, and measures the cells against the reference. The
measurement itself takes about a second per configuration; the build and the
60-digit reference are the slow steps. The scripts run the self-checks, the
build, the reference and the measurement, and record the library revision, the
build type and the canary readings beside the results. Use `-MulAddRoute
separate` or `--route separate` for the second multiply-add route.

The Windows fused run writes to `results/`, the Windows separate run to
`results-msvc-separate/`, and the GCC runs to `results-gcc-fused/` and
`results-gcc-separate/`.

## Records

`records/signed-interaction.json` and `records/signed-interaction.tsv` are the
measured factorial audit's summary: every cell of the six-axis design in each
of the four compiler and multiply-add configurations, with the signed and
absolute-value statistics and the interaction sums.

`records/synthetic-surface.json` is the controlled study that separates the two
statistics: 5,200 known-additive surfaces from 200 seeds, with the material
counts each statistic reports on them. `synthetic_surface.py` produces it;
`signed_interaction.py` reduces the measured cells.

## Requirements

A C++23 compiler, CMake 3.20 or newer, the Boys source, and Python 3.10 or
newer with `mpmath` for the reference generator.
