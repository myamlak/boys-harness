# Gate mutation runs, 2026-11-22

The accuracy gate is a falsifier: it must fail when the library's arithmetic
changes outside the published bound, and when the sentence the bound is
published in no longer matches the tree. This record deposits runs that make
three kinds of change at kernel `v5.3.3` and shows which book fails in each.

Kernel: boys `v5.3.3`, annotated tag object
`4119217c8ab23d9f370b301a9bd6338e2c804f7e`, commit
`01bba7b0ed7e0dd387f2af087e992445222e5546`, worktree restored to clean after
each run. Harness deposit: `boys-harness` commit
`3be00038384f209dd0879feece817c32465a2258`.

## Build and commands

MSVC 19.51, x64, Release, Visual Studio 18 2026, `BUILD_CUDA=OFF`. The build
is host-only: the checker below says what that costs.

```powershell
cmake -S <boys> -B build-gate -G "Visual Studio 18 2026" -A x64 -DBUILD_CUDA=OFF
cmake --build build-gate --target boys-accuracy-gate --config Release
build-gate\Release\boys-accuracy-gate.exe --strict
ctest --test-dir build-gate -C Release -R "^boys-accuracy-gate$" --output-on-failure
```

The clean tree passes: `clean-gate.txt`, exit 0, `PASS: every claim this gate
judges is met at this revision`.

## 1. Anchored sentence: `1e-13 -> 1e-14`

A one-character edit to `README.md` line 292, the bound's first paragraph:

```diff
-returned value is at or below **`1e-13 · F_0(x)`** — the ladder's first rung scaled by `1e-13`, not
+returned value is at or below **`1e-14 · F_0(x)`** — the ladder's first rung scaled by `1e-13`, not
```

No rebuild: the gate reads the source tree's `README.md` at run time.
`boys-accuracy-gate --strict` exits 2. `anchored-sentence-gate.txt` shows two
`claim D01` registration failures, one for each anchor that quotes the
sentence (`README.double.single` and `README.double.batch`), followed by
`gate: the registration check failed, so no verdict table is reported`. The
gate refuses to report a verdict table when the anchored sentence and the
tree disagree.

## 2. Default path: `kFlatMonoCoeffs`

A one-character edit to `include/boys/boys_coefficients_d_uniform_mono.hpp`
line 14:

```diff
-    9.89680267367013866e-01,  -1.02235164453333850e-02, ...
+    8.89680267367013866e-01,  -1.02235164453333850e-02, ...
```

This is `detail::kFlatMonoCoeffs`, the monomial table the shipped default
evaluation path reads. Rebuild; `boys-accuracy-gate --strict` exits 1.
`default-path-gate.txt` shows the lane book failing
`README.double.single`, `README.double.batch` and
`header.double.region_A_and_band`: the double-single region-A row moves from
`3.11e-15` to `0.102` against its `1e-13` bound (ratio `1.02e12`) at
`n = 0, x = 0.0582832`. Combination coverage falls from 96 to 95 of 192 and
the gate ends `FAIL (exit status 1)`. Under `ctest` the selected test fails:
`default-path-ctest.txt`, `ctest` exits 8.

## 3. Option paths

These one-character edits leave the default-path rows passing and trip the
book that reads the named option. Each was rebuilt and run the same way; each
gate run exits 1. The default `README.double.single` and
`README.double.batch` rows stay `[verified at this revision]` in all three.

- **`kBcoeffs`**, `include/boys/boys_coefficients.hpp` line 198:
  `{1.99358578668927033e-01,` -> `{9.99358578668927033e-01,` (first character
  1 -> 9). `option-kBcoeffs-gate.txt`: the scheme book loses
  `split-clenshaw / region B`; the stored seed row moves from
  `5.55e-17` to `0.8` (ratio `7.21e15`), and `shipped / stored seed, region
  B`, the orders entry and the single/fixed-order/all-orders entries print
  `EXCEEDED`. Scheme result 28 of 29; combination coverage 94 of 192; option
  space 13 of 16 with the scheme, the stored fit and the shipped granularity
  named as the worklist.
- **`kFlatCoeffs`**, `include/boys/boys_coefficients_d_uniform_cheb.hpp`
  line 14: `9.89728020325774338e-01` -> `8.89728020325774338e-01` (first
  character 9 -> 8). `option-kFlatCoeffs-gate.txt`: the uniform partition's
  shipped row counts 1,470 cells outside its `5.5e-14` bar (worst `n = 0`,
  `x = 4.21697e-08`, delivering `0.1`), and the run ends in the combination
  coverage check at 95 of 192 certified. The option space itself is complete
  (16 of 16), so the failure is a certified row that stopped being certified,
  not a missing combination.
- **`kNarrowBcoeffs`**, `include/boys/boys_coefficients.hpp` line 591:
  `2.43275733172444691e-01,` -> `9.43275733172444691e-01,` (first character
  2 -> 9). `option-kNarrowBcoeffs-gate.txt`: the narrow stored seed row moves
  from `5.55e-17` to `0.7` (ratio `2.72e13`); the scheme result is 23 of 29;
  option space 14 of 16 with split-clenshaw and the narrow granularity named;
  combination coverage 94 of 192.

The three option defects are reported by the scheme, granularity and
combination books, and none of them moves the anchored D01 row. That is the
intended separation: a change under an option is the option's book to catch,
while the default path's own rows stay the default path's evidence.

## Checker and scope

The build is host-only (`BUILD_CUDA=OFF`), so 0 of the 4 device lanes the
build carries were measured on a card. `tools/check_recorded_run.py` in the
kernel tree, run with `--run clean-gate.txt`, labels that run `STALE` for
exactly that reason (`checker-host-only-run.txt`). The committed run's checker
is clean (`checker-committed-run.txt`) because that run was taken from a build
with all four device lanes. These mutation runs exercise the host books; a
device-lane mutation is not carried here.

The gate's own output is deterministic. The sentence run and the default-path
gate run were re-issued on this tree and match the deposited transcripts byte
for byte, and each option-path gate output was reproduced from the
first-coefficient change in `commands.txt` and matches its transcript byte for
byte. Only the ctest timing line differs between ctest invocations.

## Files

| file | what it is |
| --- | --- |
| `clean-gate.txt` | the clean tree's `--strict` run, exit 0 |
| `anchored-sentence-gate.txt` | the edited-sentence run, exit 2 |
| `default-path-gate.txt` | the mutated default-path gate run, exit 1 |
| `default-path-ctest.txt` | the same mutation under ctest, exit 8, test #10 failed |
| `option-kBcoeffs-gate.txt` | the `kBcoeffs` mutation, exit 1 |
| `option-kFlatCoeffs-gate.txt` | the `kFlatCoeffs` mutation, exit 1 |
| `option-kNarrowBcoeffs-gate.txt` | the `kNarrowBcoeffs` mutation, exit 1 |
| `checker-committed-run.txt` | the checker on the committed four-device-lane run |
| `checker-host-only-run.txt` | the checker on `clean-gate.txt`, naming the host-only build STALE |
| `commands.txt` | the exact diffs, build and run commands, and exit codes |
| `MANIFEST.sha256` | SHA-256 of every file in this directory |
