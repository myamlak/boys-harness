# Gate mutation runs, 2026-11-28

A controlled comparison for Paper B: does the accuracy gate catch the historical
half-precision zero-return defect, and does the ordinary test suite catch it?

## Kernel and build

- Kernel: `boys` v5.3.4, annotated tag object
  `e0b6240af077381b0e3f574004f4a47dba1de277`, commit
  `bd8fadbcb5b6f33bf7df76834da0cdd2af464266`.
- Host: Windows 11 Pro, x64, MSVC 19.51.36260.0 (Visual Studio 18 2026),
  CMake 4.4.2.
- Configure: `cmake -S boys-v5.3.4 -B boys-v5.3.4/build-gate -G "Visual Studio 18 2026" -A x64 -DBUILD_CUDA=OFF`
- Build: `cmake --build boys-v5.3.4/build-gate --target boys-accuracy-gate --config Release`
- The library sources were copied from a `git archive` of the tag. The original
  repository was read only; nothing was pushed to the library.

## Clean baseline

- Gate: `boys-accuracy-gate.exe --strict`, exit 0. Final line: PASS, 38 of 39
  judged rows met.
- Ordinary CTest: `ctest -R "^(boys-tests|boys-no-vector-fallback)$"`, exit 0.

## Mutation: the historical reciprocal-form defect

The vector half routine in `src/boys_simd.cpp` originally multiplied by a
rounded reciprocal in the upward recurrence. v2.0.0 replaced that with division.
The mutation restores the multiply without changing the comment's line count, so
no registration anchor moves.

- Diff: `mut1.diff`.
- Gate: exit 1. Claim book 36 of 39 met, 2 exceeded:
  `code.half_simd_cells` and `code.half_domain_stop`. Output: `mut1-gate.txt`.
- Gate through CTest: exit 8, 11.80 s. `ctest-gate-mut1.txt`.
- Ordinary suite on the same tree: exit 0, both tests passed in 77 s.
  `ctest-ordinary-mut1.txt`. The ordinary suite does not catch the defect.
- Direct vector-versus-scalar tests:
  `boys-tests.exe --gtest_filter=BoysTest.SimdMatchesScalarF16:BoysTest.SimdMatchesScalarBf16`,
  exit 0, 2 tests passed. `gtest-simd-mut1.txt`. They fix order 8 and draw
  region-B arguments from `[11.95, 28.95)`, while the defect is at order 32 and
  the stored half argument is 11.9453125, one step below their sample.

A separate run changed the comment's line count as well as the expression. That
shifted the P10 anchor and made the gate exit 2 on the registration check before
reporting verdicts. It is recorded as `mut1-anchor.diff` and
`mut1-gate-comment-anchor.txt` in the experiment directory; the record above is
the line-count-neutral mutation, where the numerical rows decide.

## Files

`REPORT.txt` gives the full transcript and commands. The remaining files are the
diff, gate and suite outputs, build logs, source hashes and the probed point.
`MANIFEST.sha256` covers every file in this directory except itself.
