# Half-lane exhaustive sweep at v5.3.3

Re-derivation of the two published half-lane margins at tag `v5.3.3`, commit
`01bba7b0ed7e0dd387f2af087e992445222e5546` (annotated tag object
`4119217c8ab23d9f370b301a9bd6338e2c804f7e`), detached worktree, working tree clean, default
build options. Run on 21 November 2026 as a referee-requested check of the published half-lane margins.
The harness and oracle are byte-copies of the v5.1.1 deposit, with the hashes recorded in
`run-record.json`.

The comparison revision is the deposited record's: tag `v5.1.1`, commit
`834bdbb18a30ac8538b519436429f78d41f187ce`, in the control worktree `boys-v5.1.1`.

## What is measured

- **Lane entries.** `boys::BoysSingleF16(n, boys::F16)` and `boys::BoysSingleBf16(n, boys::Bf16)`,
  the scalar store-half entries, called once per `(n, x)` by the unchanged harness.
- **Bound.** `HalfBound(got, mantissaBits, minNormalExp) = 1e-7 + 0.5 * UlpOf(got, ...)`, the
  gate's slot, computed by the harness exactly as in the v5.1.1 record. The constants are
  unchanged at v5.3.3: `kBoundHalfBase = 1e-7`, `(10, -14)` for binary16 in
  `tests/boys_gate_reference.hpp`, `(7, -126)` for bfloat16 in `tests/boys_accuracy_gate.cpp`.
- **Arguments and orders.** Every positive finite bit pattern: `0x0001..0x7BFF` (binary16) and
  `0x0001..0x7F7F` (bfloat16), orders 0..32.

One default-policy cell did change between the tags: the fp16 single row moves from
`DivisionForm::kPlainReciprocal` (v5.1.1) to `DivisionForm::kRefinedReciprocal` (v5.3.3); its
granularity was already `kUniform`. The v5.3.3 changelog therefore states that the default fp16
entry has changed since the v5.1.1 record and deposits no rerun. The probe in
`build-v533.txt` static-asserts the compiled row and measures the cell's effect: the two division
forms return identical fp16 bits at all 1,047,519 cells. The bf16 row is unchanged, and the
refined form would differ from it at 2 of its 1,077,087 cells, neither of them reached by the
default.

## Counts

| format | arguments | cells |
| --- | --- | --- |
| binary16 | 31,743 | 1,047,519 |
| bfloat16 | 32,639 | 1,077,087 |
| total | 64,382 | 2,124,606 |

## Reference

The oracle is the deposited one, unchanged: 45 decimal digits, the boys generator's
`boys_ref(n, x, RAT_TAIL_FLOOR, RAT_TERMS)` with the gamma identity above x = 250. The reference
block of the two results files is identical, including the route split (271,887 binary16 and
507,045 bfloat16 cells via the identity) and the series-versus-identity overlap figures.

## Result

Worst delivered error as a fraction of the lane's budget:

| format | worst ratio | margin | worst cell |
| --- | --- | --- | --- |
| binary16 | 0.999299207482777 | 0.07007925172229923 per cent | n=0, x=0x1B82 (0.003665924072265625), delivered 0.9990234375, true 0.9987793680370024, bound 0.000244240625 |
| bfloat16 | 0.9999150781501572 | 0.00849218498427673 per cent | n=0, x=0x3FA4 (1.28125), delivered 0.69921875, true 0.6972656908714802, bound 0.001953225 |

Both worst cells are correctly rounded returns landing just above a rounding midpoint, as in the
v5.1.1 record. Each maximum is unique to 1e-9. There are no non-finite returns: 0 in either
format. Zero returns: 369,692 binary16, 456,092 bfloat16.

## Bound check

The largest ratio in either format is below 1 (0.999299207482777 and 0.9999150781501572), so no
cell exceeds its bound. There is no cell to list under a ratio-greater-than-one heading. The
published half-lane bound passes.

## Top ten ratios per format

binary16:

| # | n | x bits | x | delivered | bound | ratio |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | 0 | 0x1B82 | 0.003665924072265625 | 0.9990234375 | 0.000244240625 | 0.999299207482777 |
| 2 | 0 | 0x3FBE | 1.935546875 | 0.60546875 | 0.000244240625 | 0.999225144839333 |
| 3 | 0 | 0x1D42 | 0.00513458251953125 | 0.99853515625 | 0.000244240625 | 0.999221904928736 |
| 4 | 1 | 0x3780 | 0.46875 | 0.25341796875 | 0.0001221703125 | 0.999164866289538 |
| 5 | 0 | 0x1201 | 0.0007328987121582031 | 0.99951171875 | 0.000244240625 | 0.999159676901074 |
| 6 | 0 | 0x421C | 3.0546875 | 0.50048828125 | 0.000244240625 | 0.999097804253875 |
| 7 | 1 | 0x142B | 0.0010175704956054688 | 0.333251953125 | 0.0001221703125 | 0.999096647164537 |
| 8 | 0 | 0x3928 | 0.64453125 | 0.8212890625 | 0.000244240625 | 0.999082315987290 |
| 9 | 0 | 0x48A6 | 9.296875 | 0.29052734375 | 0.0001221703125 | 0.999007645679590 |
| 10 | 0 | 0x1881 | 0.0021991729736328125 | 0.9990234375 | 0.000244240625 | 0.998966593577407 |

bfloat16:

| # | n | x bits | x | delivered | bound | ratio |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | 0 | 0x3FA4 | 1.28125 | 0.69921875 | 0.001953225 | 0.999915078150157 |
| 2 | 5 | 0x3D31 | 0.043212890625 | 0.087890625 | 0.000244240625 | 0.999592093804728 |
| 3 | 7 | 0x3991 | 0.0002765655517578125 | 0.06689453125 | 0.000244240625 | 0.999551702399791 |
| 4 | 0 | 0x4073 | 3.796875 | 0.451171875 | 0.0009766625 | 0.999441474945192 |
| 5 | 6 | 0x3994 | 0.00028228759765625 | 0.07666015625 | 0.000244240625 | 0.999439972022973 |
| 6 | 0 | 0x4084 | 4.125 | 0.43359375 | 0.0009766625 | 0.999341220564905 |
| 7 | 3 | 0x3BE3 | 0.006927490234375 | 0.1416015625 | 0.00048838125 | 0.999293755278465 |
| 8 | 4 | 0x399C | 0.00029754638671875 | 0.111328125 | 0.000244240625 | 0.999261099799825 |
| 9 | 6 | 0x3993 | 0.0002803802490234375 | 0.0771484375 | 0.000244240625 | 0.999220672516887 |
| 10 | 4 | 0x399D | 0.0002994537353515625 | 0.11083984375 | 0.000244240625 | 0.999210278414804 |

## Comparison with the v5.1.1 record

The old raw sweep was available at
`<previous run root>\half_sweep.bin` (46,741,352 bytes), so the comparison is
cell by cell, not only figure by figure.

| figure | v5.1.1 | v5.3.3 | difference |
| --- | --- | --- | --- |
| binary16 worst ratio | 0.999299207482777 | 0.999299207482777 | 0.0 |
| binary16 margin | 0.07007925172229923 per cent | 0.07007925172229923 per cent | 0.0 |
| bfloat16 worst ratio | 0.9999150781501572 | 0.9999150781501572 | 0.0 |
| bfloat16 margin | 0.00849218498427673 per cent | 0.00849218498427673 per cent | 0.0 |
| cells with changed returned bits | - | 0 binary16, 0 bfloat16 | 0 |
| cells with changed bounds | - | 0 binary16, 0 bfloat16 | 0 |

The worst cells are the same cells: n=0, x=0x1B82 and n=0, x=0x3FA4. The two results files'
`formats`, `orders` and `reference` blocks are identical; only `wall_seconds` differs. The sweep
files themselves are byte-identical (same SHA256,
`6F0F4CC23CDEDE6F7A116E804D5FEE9DCE6DC3D887B280A10EEC9D1FB53327A7`).

A control rebuild of the v5.1.1 headers with today's compiler reproduces the old raw binary with
0 changed bits, so the identity is not an artefact of a changed pipeline.

Revision note: the oracle hardcodes `library_tag: v5.1.1` and the v5.1.1 commit in its output.
`results-v533.json` is its verbatim output; the sweep it scores is v5.3.3 / `01bba7b`, as
`configure.log` and `run-record.json` record.

## H04

The sentence still holds, with the same figures to the digits it states and to every digit
measured:

| figure | H04 | measured at v5.3.3 |
| --- | --- | --- |
| binary16 worst ratio | 0.999299207 | 0.999299207482777 |
| binary16 margin | 0.0700793 per cent | 0.07007925172229923 per cent |
| bfloat16 worst ratio | 0.999915078 | 0.9999150781501572 |
| bfloat16 margin | 0.00849218 per cent | 0.00849218498427673 per cent |

It does not change and does not fail. The v5.3.3 defaults pass the published half-lane bound.

## Runs

Toolchain: cl 19.51.36260.0 (Visual Studio 18, toolset 14.51.36231), Windows SDK 10.0.26100.0,
CMake 4.4.2, generator "Visual Studio 18 2026"; Python 3.14.7, mpmath 1.4.1, numpy 2.5.2. The
deposit's build.txt records cl 19.51.36256.0; the toolset directory is the same. Origin is
`https://github.com/myamlak/boys.git` (the brief said `git@github.com`; the URL differs, the tag
and commit do not). Exact commands are in `build-v533.txt`.

| step | exit | wall time |
| --- | --- | --- |
| configure | 0 | 18.86 s |
| library (`--config Release --target boys`) | 0 | 249.61 s |
| harness compile and link | 0 | 14.77 s |
| sweep | 0 | 0.547 s |
| reference (`--jobs 8`) | 0 | 155.4 s (reference 154.2 s) |
| v5.1.1 control compile / sweep | 0 / 0 | 15.53 s / under 1 s |
| policy probe compile / run | 0 / 0 | 17.99 s / 0.564 s |

The oracle reference pass was faster than the deposited 308.1 s; the machine, Python, mpmath,
numpy and `--jobs 8` are otherwise as recorded.

## Files

- `half_sweep.cpp` - the harness, byte-copy of the deposit
  (`B969BDDC3AD5D7AAAD66559C4DE56CE29414837AB7EF65D49DAFB4581C064B0F`).
- `half_sweep_oracle.py` - the 45-digit oracle, byte-copy of the deposit
  (`24E7BE9A2C0F74132F06DB8854717F9AEB08FE68BF80A90F9E50EB074159087E`).
- `half_sweep.bin` - the v5.3.3 sweep, 46,741,352 bytes
  (`6F0F4CC23CDEDE6F7A116E804D5FEE9DCE6DC3D887B280A10EEC9D1FB53327A7`), byte-identical to the
  v5.1.1 raw sweep.
- `results-v533.json` - the oracle's output for the v5.3.3 sweep
  (`9EC27E0D1B9492FA323CE6B416538725EE2FDFC155DC563DF62CA25DB8DC0BFC`); its revision block is
  the oracle's hardcoded v5.1.1 string.
- `results-v511.json` - the deposited v5.1.1 results, byte-copy
  (`E665159525BD546E5BDDE2242EB1A196FD5C7742917180741994FADB38B35372`).
- `compare-v511-v533.json` - cell-by-cell change counts and the changed-cell lists (empty).
- `run-record.json` - provenance, commands, exit codes, durations and digests.
- `build-v533.txt` - exact build and run commands, as used.
- Logs as run: `configure.log`, `build-boys.log`, `compile-harness.log`, `oracle-v533.log`,
  `probe-policy.log`, `compile-v511-control.log`.
- `probe_policy.cpp` - the default-row and division-form probe.
- `half_sweep-v511-rebuild.bin` - the v5.1.1 control sweep, byte-identical to the deposited raw
  binary.
