# Call-time probe for the Boys/libint2 comparison

This directory is the timing record cited as Table 4 of the claims-as-rows
paper. It is the companion to the accuracy comparison in the parent directory;
that record measures accuracy only and times nothing.

## Files

- `bench_boys_arms.cpp` — the probe.
- `run-gate.txt` — the recorded output.
- `environment.txt` — host, compiler, flags, grid and protocol.

## What is measured

Four arms, each over the 1,797 unique arguments of the committed gate grid
(`tests/data/boys_accuracy_gate_reference.csv` at tag `v5.3.3`), orders 0 to 32:

- `boys_ladder`: `boys::BoysAllOrders`
- `boys_batch`: `boys::BoysAllN` with a workspace
- `cheb7`: `libint2::FmEval_Chebyshev7<double>` (libint2 2.13.1)
- `reference2`: `libint2::FmEval_Reference2<double>`

Each arm is warmed up twice and timed five times, and the median is printed.
The checksum is printed so that a rerun can confirm it evaluated the same
values.

## Build and run

The record was produced in WSL2 with g++ 15.2.0 at `-O3 -mavx2 -mfma`, linked
against the installed libint2 2.13.1. The include and library paths depend on
the build environment; the source requires the boys and libint2 headers. The
grid loader reads the second CSV column and keeps the unique arguments.

```
g++ -O3 -mavx2 -mfma bench_boys_arms.cpp \
    -I<boys>/include -I<libint2>/include -L<libint2>/lib -lint2 \
    -o bench_boys_arms
./bench_boys_arms <boys>/tests/data/boys_accuracy_gate_reference.csv
```

## Scope

One host and one build. The numbers are microbenchmark medians, not a
library-quality claim, and no accuracy claim is made here.
