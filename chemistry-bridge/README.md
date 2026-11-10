# Chemistry bridge

The integral-engine bridge of Paper B: the same host tensor built through two
Boys evaluation routes and compared at the floating-point contract and in the
total energy. `scripts/probe.cpp` builds the tensors through each route,
`scripts/analyze.py` reduces the records, and `scripts/bridge-CMakeLists.txt`
is the build stanza; `scripts/commands.md` holds the exact commands and the
harness settings of each run.

## Records

`outputs/` holds, for each basis and each route, the captured Boys argument
streams, the per-order reference comparisons, the integral matrices, the
density dumps and the energy summaries. The two cc-pVTZ argument streams are
stored xz-compressed; decompress them in place with `xz -dk` before running the
analysis.

The two dense cc-pVTZ integral tensors (about 260 MB each) are not deposited.
`scripts/probe.cpp` regenerates them from the deposited argument streams and
the library. The probe commands are in `scripts/commands.md`; the cc-pVTZ arms
are run from the host scratch tree as
`qcx-bridge-probe.exe <out>/ccpvtz-default data/basis/cc-pvtz` and the same
with `paper-a` for the second arm. The benzene arms are in
`benzene-pvdz/commands.md` as
`qcx-bridge-probe.exe <out-prefix> data/basis/cc-pvdz benzene`.

## Above-B2 route difference

The deposited argument streams in `outputs/` carry the class-row and paper-a
values at every distinct `(nmax, x)`. Above B2 (`x >= 11.899848152108484`) the
largest difference between the two routes is `5.5511151231257827e-17` at order
0, the same in all three bases: at `x = 11.974444364155501` in STO-3G,
`12.017923166818386` in 6-31G* and `11.90786820012136` in cc-pVTZ. For
comparison, the largest difference over B2 is `1.1102230246251565e-16` at
order 0, and over all arguments `3.219646771412954e-15` at order 0.
`scripts/above_b2_max.py` prints these values from the deposited streams.

## The water geometry

The molecule is the host repository's fixture, `tests/fixtures/h2o_sto3g.cpp`
in `bladawiec` (SHA-256
`ef2262e4d008f8e96bc2b596dc35be8926f8a1d3a449cd09399b729713e2086c`): O at
the origin and the two H at
`(±1.430428808474167, 1.107157044080814, 0)` bohr, in the xy plane. With the
CODATA 2018 conversion 1 bohr = 0.529177210903 A, the coordinates give
r(OH) = 1.8088458464917874 bohr = 0.9572 A and H-O-H = 104.52 degrees: the
gap to the documented geometry is 0.0 A and 0.0 degrees.

The bridge record names the host tree `bladawiec` at branch
`myamlak/phase2-pr3`, at commit `ac8540b7c23d3818c9e67f8d59bf6cfbbef7fa0d` or
later, for the scratch copy that ran the probe; the host tree itself is not
deposited. The fixture `tests/fixtures/h2o_sto3g.cpp` has SHA-256
`ef2262e4d008f8e96bc2b596dc35be8926f8a1d3a449cd09399b729713e2086c`, is
byte-identical at `ac8540b7`, and its content last changed in commit
`55b9a8251c155d19550d795c2298422b8279fe0f` (`v0.1.0`).

## Benzene/cc-pVDZ

`benzene-pvdz/` is the second molecule: benzene, 12 atoms and 114 basis
functions, with 23 771 671 Boys calls over 5 330 629 distinct arguments. It
holds the probe that produced the records (`probe.cpp`), the comparison script
(`compare_arms.py`), the exact commands (`commands.md`) and, in `outputs/`, both
arms' summaries, integral matrices, density dumps and SCF traces plus the
comparison's own output. Its two dense ERI dumps and two argument streams are
not deposited; the probe regenerates them in about twenty minutes per arm.
