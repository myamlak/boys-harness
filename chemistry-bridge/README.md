# Chemistry bridge

The integral-engine bridge of the accuracy paper: the same host tensor built
through two Boys evaluation routes and compared at the floating-point contract
and in the total energy. `scripts/probe.cpp` builds the tensors through each
route, `scripts/analyze.py` reduces the records, and `scripts/bridge-CMakeLists.txt`
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
the library, and they are available from the author.

## Benzene/cc-pVDZ

`benzene-pvdz/` is the second molecule: benzene, 12 atoms and 114 basis
functions, with 23 771 671 Boys calls over 5 330 629 distinct arguments. It
holds the probe that produced the records (`probe.cpp`), the comparison script
(`compare_arms.py`), the exact commands (`commands.md`) and, in `outputs/`, both
arms' summaries, integral matrices, density dumps and SCF traces plus the
comparison's own output. Its two dense ERI dumps and two argument streams are
not deposited; the probe regenerates them in about twenty minutes per arm.
