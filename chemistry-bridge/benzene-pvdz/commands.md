# Benzene/cc-pVDZ: commands, records and what is not deposited

The chemistry bridge's second molecule. The same dense ERI tensor of benzene is
built twice through the two Boys routes — the library's default class row and the
paper's named recursed policy — and the two tensors are compared value by value
and in the total energy. `probe.cpp` is the probe that produced the records
here, `bridge-CMakeLists.txt` its build stanza, and `compare_arms.py` the
reduction that compares the two arms' dense dumps.

## Scratch tree and edits

The runs were made in `C:\Users\myaml\dev\bladawiec-bridge`, a copy of the
framework tree at version 0.2.0 carrying boys 5.2.0, with the bridge's own
`QCX_BRIDGE_BOYS_POLICY` (`default`/`paper-a`) and `QCX_BRIDGE_BOYS_CAPTURE`
switches (see `../scripts/commands.md` for those edits). This record adds four:

* `bridge/probe.cpp`: a third argument selecting the molecule (`water` or
  `benzene`); the benzene fixture's basis and element filter (`{1, 6}`); and the
  dense-tensor admission cap raised to 4 GiB for benzene, whose estimate is
  2.52 GiB against the 2 GiB default.
* `bridge/CMakeLists.txt`: the benzene fixture source in the executable.
* `integrals/src/internal/md_boys.hpp`: the capture sink deduplicates as it
  records — one entry per distinct `(nmax, x)` — instead of a per-call log. A
  molecule of this size forms 23.8 million arguments, and the per-call vector
  alone is 6.7 GB. The captured records, their counts and the ladder
  comparisons are the same either way; the change moves memory, not figures.
* `scf/src/rhf.cpp`, `scf/src/uhf.cpp`: the J/K supermatrices are captured into
  the Fock lambda by reference, not by value (a second copy of the two n^4
  supermatrices, 2.7 GB at this size).

A water/cc-pVDZ run from the edited probe reproduces the earlier deposit's
`args.txt`, `eri.txt`, `eri_max.txt` and `matrices.txt` byte for byte, which is
what makes the edits safe to state as measurement-neutral.

## Build and run

Two builds of one tree, one per arm (`QCX_BRIDGE_BOYS_POLICY=default` and
`=paper-a`); each is configured with `QCX_BUILD_BRIDGE_PROBE=ON`,
`QCX_BRIDGE_BOYS_CAPTURE=ON` and Release.

```
qcx-bridge-probe.exe <out-prefix> data/basis/cc-pvdz benzene
```

| Arm | Peak working set | Wall |
| --- | --- | --- |
| default | 5.80 GB | 1238 s |
| paper-a | 5.80 GB | 711 s |

Both arms form the same stream: 23 771 671 calls over 5 330 629 distinct
arguments, with no repeated argument returning a different ladder.

## The comparison

```
python compare_arms.py benzene-pvdz-default.eri.txt benzene-pvdz-paper-a.eri.txt \
    benzene-pvdz-default.matrices.txt seed_density
```

`outputs/benzene-pvdz-compare.txt` is that command's own output:

* 168 896 016 values compared (114^4); 71 492 684 of them (42.3%) differ;
* the largest single difference is 3.99680288865056355e-15;
* the two-electron energy at the shared seed density is identical in double
  precision (below half an ulp of 268.8, 2.8e-14).

The converged totals are each arm's own SCF: −230.72206322481188 (default gate)
against −230.72206322481151, a difference of 3.70e-13 hartree, and −230.72206322481205
against −230.72206322481193 at the tight gate, 1.14e-13.

## What is not deposited

The two dense ERI dumps (`*.eri.txt`, 4.4 GB each) and the two argument streams
(`*.args.txt`, 968 MB each) are not deposited. `probe.cpp` regenerates both from
the library and the basis directory, and `compare_arms.py` reduces the dumps;
the summaries, matrices, density dumps and SCF traces deposited beside this file
are what the paper quotes.
