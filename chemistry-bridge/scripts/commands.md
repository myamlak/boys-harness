# Chemistry bridge: commands and exit codes

Lane: Paper-A chemistry bridge.
Scratch tree: `C:\Users\myaml\dev\bladawiec-bridge`, a copy of
`C:\Users\myaml\dev\bladawiec` at branch `myamlak/phase2-pr3`, HEAD `ac8540b7` or
later, with no `.git` and no build directories. The boys content is v5.0.0
(`boys::kVersionString = "5.0.0"` read from the copied header; the content hash
`482fd33` is the pin stated in the lane brief, and the copy carries no git
metadata, which is why `BoysBuildRecord` prints `submodule revision unknown`).

No source file in `C:\Users\myaml\dev\bladawiec` was edited and no git command was
run there. The only writes into the real tree are the files under this
`.claude/lane-status/chemistry-bridge/` directory, which the lane brief names as
the output location.

## 1. Scratch copy and edits

The copy was made before this record's first command (pre-compaction). The first
attempt through the Bash tool exited 16 because MSYS rewrote `/E` as a path; the
same robocopy through PowerShell produced the tree. The copy is verified here,
not assumed:

| Purpose | Command | Exit |
| --- | --- | --- |
| Verify the scratch tree | `test -f bridge/probe.cpp && test -f integrals/src/internal/md_boys.hpp && test -f external/boys/include/boys/boys.hpp && test -f tests/fixtures/h2o_sto3g.cpp` | 0 |
| Read the boys version | `grep kVersionString external/boys/include/boys/version.hpp` | 0, prints `5.0.0` |

Edits in the scratch copy only (the diffs against the real tree are in
`scratch-edits/`):

* `integrals/src/internal/md_boys.hpp`: the
  `QCX_BRIDGE_BOYS_POLICY_PAPER_A` policy switch inside `MdBoysBatch`'s fp64 arm
  (default unchanged), the `PaperABoysPolicy` alias, the coarse Chebyshev
  coefficient-header include ahead of `qcx/integrals/boys.hpp`, and the
  `QCX_BRIDGE_BOYS_CAPTURE` argument/ladder capture.
* `integrals/CMakeLists.txt`: the `QCX_BRIDGE_BOYS_POLICY` cache variable
  (`default`/`paper-a`) and `QCX_BRIDGE_BOYS_CAPTURE`, both set as compile
  definitions on the `qcx-boys-route` INTERFACE target, the same seam the
  existing `QCX_BOYS_SHARE_INSTRUMENT` uses.
* `CMakeLists.txt`: the `QCX_BUILD_BRIDGE_PROBE` option and the `bridge`
  subdirectory.
* `bridge/CMakeLists.txt`, `bridge/probe.cpp`: the probe itself (new files).

## 2. Configure, build and run

Both build directories use the same generator and flags; the only difference is
the policy cache variable.

| Purpose | Command | Exit |
| --- | --- | --- |
| Configure default arm | `cmake -S C:/Users/myaml/dev/bladawiec-bridge -B .../build-bridge-default -G "Visual Studio 18 2026" -A x64 -DCMAKE_TOOLCHAIN_FILE=C:/vcpkg/scripts/buildsystems/vcpkg.cmake -DQCX_BUILD_TESTS=OFF -DQCX_BUILD_BENCHMARKS=OFF -DQCX_ENABLE_IO=OFF -DQCX_INTEGRALS_LMAX=2 -DQCX_BUILD_BRIDGE_PROBE=ON -DQCX_BRIDGE_BOYS_CAPTURE=ON -DQCX_BRIDGE_BOYS_POLICY=default` | 0 |
| Configure paper-a arm | the same with `-B .../build-bridge-paper-a` and `-DQCX_BRIDGE_BOYS_POLICY=paper-a` | 0 |
| Build default arm | `cmake --build .../build-bridge-default --config Release --target qcx-bridge-probe` | 0 (background task `bwlx9kctb`) |
| Build paper-a arm | `cmake --build .../build-bridge-paper-a --config Release --target qcx-bridge-probe` | 0 (background task `b4pcw34lb`) |
| Rebuild default arm after the probe fix | `cmake --build .../build-bridge-default --config Release --target qcx-bridge-probe` | 0 (only `probe.cpp` recompiled) |
| Run default arm | `qcx-bridge-probe.exe C:/Users/myaml/dev/bladawiec-bridge/bridge-runs/default` | 0 |
| Run paper-a arm | `qcx-bridge-probe.exe C:/Users/myaml/dev/bladawiec-bridge/bridge-runs/paper-a` | 0 |
| Offline analysis | `python -I analyze.py <default-prefix> <paper-a-prefix> analysis.txt` | 0 |

The configure logs show the arm's own status line: the paper-a log prints
`qcx-integrals: the Paper-A bridge Boys policy is in force`; the default log does
not. Both print the capture line. The probe runs need the vcpkg `bin` directory
on `PATH` (`mkl_sequential.2.dll`); the paper-a build also copies that DLL beside
the executable.

## 3. One correction made during the run

The first default-arm probe run wrote a Coulomb-only fixed-density two-electron
energy, and its own trace-identity check exposed it: `Tr[DH] + E_2e` missed
`HfResult.electronicEnergy` by `9.1048138174668338` Ha, exactly the exchange
term. `probe.cpp` was corrected to evaluate the full closed-shell expression

`E_2e = 1/2 sum D_uv D_ls (uv|ls) - 1/4 sum D_uv D_ls (ul|sv)`,

matching the SCF's own Coulomb and exchange supermatrices (`scf_common.hpp`,
`E[(mu*n+nu),(l*n+s)] = eri(mu,nu,l,s)` and `eri(mu,l,s,nu)`). After the fix the
identity residual is 0.0 at the default-gate density and at most 2.8e-14 Ha
(different summation orders) at the tight-gate density. The SCF totals were never
affected; only the fixed-density reporting was. Both arms were rebuilt with the
corrected probe before the runs recorded above.

## 4. Outputs copied here

* `logs/`: both configure logs, the three build logs, both probe run logs.
* `outputs/`: for each arm, the summary, ERI tensor, ERI maximum, matrices,
  captured argument stream, Boys build record, density dump and both SCF traces.
* `scratch-edits/`: the diffs and the new probe sources, plus the H2O/STO-3G
  fixture the probe ran.
* `analysis.txt`: the analysis output (numbers quoted in the lane report).

## 5. Run facts

* Input: H2O / STO-3G (the repository's own test fixture), n = 7 orbitals,
  RHF, `E_nuc = 9.194964854210685` Ha.
* Boys argument stream: 10,935 engine calls, 2,764 distinct `(nmax, x)` pairs,
  `nmax` in {0,1,2,3,4}, x in [0, 28.033987699960576].
* The two arms' captured streams are identical `(nmax, x bits, multiplicity)`;
  each arm's engine ladder is reproduced bit for bit by the policy it was built
  with; no repeated call on the same argument returned a different ladder.

## 6. Realistic basis: H2O/6-31G*

The probe was extended in the scratch copy only: its second argument is an
optional NWChem basis directory, parsed with
`ParseNwchemDirectoryFiltered(directory, {1, 8})` (the water molecule's
elements) against the same fixture geometry. The summary now also carries
`basis_source`, `basis_directory`, `shells`, `shell_pairs`,
`canonical_quartets`, `screened_quartets` (the dense builder's default
`kNormal` Schwarz screen) and per-order `calls_nmax_N` / `unique_nmax_N`
histograms. The probe was rebuilt; the library objects were not touched, so
both arms are the same frozen LMAX=2 binaries that produced the STO-3G run.
LMAX=2 is sufficient for 6-31G* because the engine canonicalises every task
to `lBra <= lKet` (`CanonicalTaskFormOf`) and the 2e triangle at LMAX=2
covers every canonical class up to `2 * 2 = 4`; the d shell reaches exactly
class (4,4).

| Purpose | Command | Exit | Notes |
| --- | --- | --- | --- |
| Rebuild default probe | `cmake --build build-bridge-default --config Release --target qcx-bridge-probe` | 0 | only `probe.cpp` recompiled |
| Rebuild paper-a probe | the same in `build-bridge-paper-a` | 0 | only `probe.cpp` recompiled |
| Run default arm | `qcx-bridge-probe.exe <out>/631gstar-default data/basis/6-31g-star` | 0 | 2 s wall |
| Run paper-a arm | `qcx-bridge-probe.exe <out>/631gstar-paper-a data/basis/6-31g-star` | 0 | 1 s wall |
| Analyse | `python -I analyze.py <default> <paper-a> analysis-631gstar.txt` | 0 | 11 s wall |

Basis and stream facts (`outputs/631gstar-default.summary.txt`,
`analysis-631gstar.txt`):

* `data/basis/6-31g-star`, NWChem format, `SPHERICAL`; n = 18, 10 shells,
  55 shell pairs, 1,540 canonical quartets, all 1,540 survive the default
  `kNormal` Schwarz screen (the dense stream is therefore the screened one,
  and here the screen drops nothing).
* 51,570 engine calls over 14,163 distinct `(nmax, x)` arguments; `nmax` in
  0..8. Calls by order: 0:23,436; 1:14,040; 2:9,627; 3:3,012; 4:1,199;
  5:196; 6:54; 7:4; 8:2.
* x per order (max over distinct arguments requesting it): order 0..4 reach
  x up to 121.7/121.7/121.7/74.5/73.5; order 5 up to 26.77; order 6 up to
  4.450; orders 7 and 8 each have exactly one distinct argument and it is
  x = 0. Those two are the same-centre (d,d)|(d,d) and (p,d)|(d,d) quartets:
  the d shell reaches order 8, but only where both routes are exact.
* Both arms' captured streams are identical and each arm reproduces its own
  policy's ladder bit for bit; no repeated-call mismatch.

The numbers are in `analysis-631gstar.txt` and quoted in the lane report.

## 7. Larger realistic basis: H2O/cc-pVTZ

The f shell on oxygen makes order 12 reachable and extends order above 6 to
x > 0 for the first time (6-31G* reached order 6 at x > 0, and orders 7-8 only
at x = 0). LMAX=3 is the minimum that admits f shells: the engine canonicalises
every task to `lBra <= lKet`, so the 2e triangle needs every canonical class up
to `2 * 3 = 6`, which is the 46-class LMAX=3 matrix (24 classes at LMAX=2).
Both arms were reconfigured from LMAX=2 to LMAX=3; that forces a rebuild of
every library, not just the probe.

| Purpose | Command | Exit | Notes |
| --- | --- | --- | --- |
| Reconfigure default arm | `cmake -S . -B build-bridge-default -DQCX_INTEGRALS_LMAX=3` | 0 | 1.5 s; capture line present |
| Build default arm | `cmake --build build-bridge-default --config Release --target qcx-bridge-probe` | 0 | 981 s (background task `bhtoxkfzq`; only `probe.cpp` was left) |
| Reconfigure paper-a arm | `cmake -S . -B build-bridge-paper-a -DQCX_INTEGRALS_LMAX=3` | 0 | 5 s; paper-a policy line present |
| Build paper-a arm | `cmake --build build-bridge-paper-a --config Release --target qcx-bridge-probe` | 0 | 1170 s (background task `bgbb7cdjy`) |
| Run default arm | `qcx-bridge-probe.exe <out>/ccpvtz-default data/basis/cc-pvtz` | 0 | 69 s |
| Run paper-a arm | `qcx-bridge-probe.exe <out>/ccpvtz-paper-a data/basis/cc-pvtz` | 0 | 70 s |
| Analyse | `python -I analyze.py <default> <paper-a> analysis-ccpvtz.txt` | 0 | 389 s; lockstep streaming path |
| Rebuild probes with the census lines | `cmake --build <arm> --config Release --target qcx-bridge-probe` | 0 | 36 s / 41 s; only `probe.cpp` |
| Re-run both arms | the two run commands above | 0 | 67 s / 61 s; SCF totals bit-identical |
| Re-analyse all three bases | the same analysis command for each | 0 | 389 s / 12 s / 3 s |

The builds are a one-off cost of the LMAX change; the measurement itself is
about one minute per arm, far inside the bounded hour.

Basis and stream facts (`outputs/ccpvtz-default.summary.txt`,
`analysis-ccpvtz.txt`):

* `data/basis/cc-pvtz`, NWChem format, `SPHERICAL`; n = 58.
* The engine sees **10 contracted shells**, 55 shell pairs and 1,540 canonical
  quartets, and all 1,540 survive the default `kNormal` Schwarz screen. The
  shell count equals 6-31G*'s by coincidence, not by stale data: an NWChem
  general-contraction block is one engine shell with several contraction rows.
  O is 4 shell objects (S with 4 rows, P with 3, D with 2, F with 1) carrying
  30 functions; each H is 3 shell objects carrying 14. The summary now records
  `molecule_atoms=3`, `pair_function_count=58` and the per-element counts, so
  the census can be checked against n directly.
* 262,484 engine calls over 81,627 distinct `(nmax, x)` arguments; `nmax` in
  0..12; x in [0, 277.20922853176552]. Calls by order: 0:63,125; 1:49,500;
  2:61,300; 3:40,060; 4:27,624; 5:12,372; 6:5,825; 7:1,848; 8:651; 9:136;
  10:37; 11:4; 12:2.
* Order 12 is two calls on one distinct argument, and it is x = 0: the
  diagonal same-centre (f,f)|(f,f) quartet. Orders 9-11 do occur at x > 0
  (max x 7.960, 3.975, 0.7866); orders 7-8 reach 53.61 and 10.19. This is the
  first basis where orders 7-11 are exercised away from x = 0.
* Both arms' captured streams are identical and each arm reproduces its own
  policy's ladder bit for bit; no repeated-call mismatch.

### Two corrections found by this run

1. `analyze.py`'s streaming parser never converted the right-hand indices to
   `int`, so the lockstep assertion compared a string tuple with an int tuple
   and could never pass; it would have aborted the cc-pVTZ analysis on its
   first data line. Found with a synthetic 3^4 tensor and fixed; the streamed
   Coulomb/exchange contraction now matches a brute-force contraction of the
   same tensors to 4e-31 (summation-order noise), with the entry statistics
   equal.
2. The 80-digit reference used the upward recurrence above x = 1e-6 at every
   order; the validation exposed 0.011 error at x = 1.0017e-6 and order 12.
   The series branch is now used below x = 1e-2, where the recurrence still
   carries about 40 digits at order 12; the validation reads at most 2.7e-67
   at the sampled arguments and at most 2.6e-82 at the series checks. The
   per-order tables for all three bases are unchanged: the reported maxima sat
   outside the unstable window. The three analysis files carry the corrected
   validation.

### Numbers

* ERI: 11,316,496 entries; 5,078,162 bitwise differing; max |ERI| 4.741181;
  max |ΔERI| 1.554312e-15 at (3,3,3,6); RMS 3.711410e-17; ratio 3.278e-16.
* Fixed density: ΔE_2e at the shared seed = +1.707278e-15 Ha (1.071e-12 x
  chemical accuracy; Coulomb +1.963e-15, exchange -2.559e-16); at the
  converged and tight densities +2.081275e-15 Ha (1.306e-12). The probe's two
  printed seed values are identical because 1.7e-15 is below half an ulp of
  39.93 (ulp 7.1e-15); the streamed per-entry contraction resolves it.
  ΔE_2e/(1e-13 |E_2e(seed)|) = 4.276e-04.
* SCF, default gate: 13 iterations both; -76.057168514878114 vs
  -76.057168514878143 Ha; Δ = -2.842171e-14 Ha (3.737e-03 of the 1e-13 |E|
  scale), achieved dE 5.032e-09 both. Tight gate: 52 iterations both;
  -76.057168514878029 vs -76.057168514878072 Ha; Δ = -4.263256e-14 Ha
  (5.605e-03), achieved dE 8.242e-13 / 9.663e-13.
* S and core are bit-identical; the seed density is bit-identical;
  default_density differs by at most 1.566e-12 and tight_density by 8.555e-13.
* Per-order: class-row worst error 3.512e-02 x bound (order 4, x = 4.195e-03);
  paper-a worst 4.791e-02 (order 8, x = 2.0506); largest route difference
  3.219647e-15 (order 0, x = 8.069e-08) = 3.225e-02 x 1e-13 F_0; order 12
  route difference exactly 0.
