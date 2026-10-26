# Recount of the two-library audit counts

This deposit re-derives, from the pinned source trees, every count the paper's
audit pilot attributes to Boost.Math 1.92.0 and crlibm 1.0beta5, and checks the
measurement totals against the audit's deposited session registry. It is the
independent recount cited by the paper. The script reads the trees; it never
builds or runs anything from them and writes only its own output.

Bottom line: every number in the paper holds. Two of the source counts (499
Boost bounds, 63,617 crlibm vectors) are also reproducible with tools deposited
in `../boost_accuracy_rows/tools` and `../claim-audit-crlibm/work`. Two (791
Gappa obligations, 41 source-comment bounds) are reproducible only under a
counting rule this deposit had to infer, because the registry carries the
numbers without a counter. The 12-failure split, the 44 undecided points and
the Table S9 totals are readable from the registry but cannot be recounted: the
per-query output was never saved, as `../README.md` also states.

## Contents

- `recount.py` — the counter. Standard library only; reads the pinned trees and
  the registry, writes `counts.json`.
- `counts.json` — this recount's output, produced from this directory by the
  command below. It contains the provenance checks, per-file breakdowns, every
  counted record, the registry values, and a claim-by-claim comparison.
- `expected_bounds.csv` — output of the audit's deposited
  `../boost_accuracy_rows/tools/expected_results.py`, run as a cross-check of
  the 499 Boost records. It is not an input to `recount.py`.
- `registry/counts.json` — a copy of the audit's session registry
  (`../claim-audit-crlibm/work/counts.json`), so the comparison is
  self-contained. It is an input, not this recount's output.

## Pinned sources

- Boost.Math 1.92.0 release archive, SHA-256
  `c4a3b310ddd2472416e091067166b0713be97c63f38c212c484ada022fd296ce`.
- crlibm commit `eb3063791aa75bc9705b49283bf14250465220a7`; its `VERSION` file
  reads `1.0beta5`. The release tarball's top directory names the commit.
- The same two pins are recorded in `../witness-rerun/README.md`.

The Boost witness tree in the audit record is headers-only: it has no
`libs/math/test`, so the 499 test-suite bounds cannot be counted from it. Pass
the full 1.92.0 source tree as `--boost`, and the witness tree as
`--boost-headers`; the script then checks that every file under `boost/math` is
byte-identical between the two trees (521 of 521) and that both report
`BOOST_VERSION 109200`. The full tree used for the recorded run was
`C:/Users/myaml/dev/_boost-smoke/boost-1_92-src/boost_1_92_0`; the witness tree
was `C:/Users/myaml/dev/_witness/third_party/boost_1_92_0`.

## Reproduce

Run from `accuracy-claim-audit/recount`. The two deposit-relative paths are
`registry/counts.json` (input) and `counts.json` (output). The five paths to the
pinned sources are external and must be supplied by the reader; the values below
are those used for the recorded run.

```sh
cd accuracy-claim-audit/recount
python -I recount.py \
  --boost          C:/Users/myaml/dev/_boost-smoke/boost-1_92-src/boost_1_92_0 \
  --boost-headers  C:/Users/myaml/dev/_witness/third_party/boost_1_92_0 \
  --boost-archive  C:/Users/myaml/dev/_witness/downloads/boost_1_92_0.tar.gz \
  --crlibm         C:/Users/myaml/dev/_witness/third_party/crlibm-1.0beta5 \
  --crlibm-archive C:/Users/myaml/dev/_witness/downloads/crlibm-eb306379.tar.gz \
  --registry       registry/counts.json \
  --out            counts.json
```

`--boost`, `--crlibm` and `--registry` are required. `--boost-headers`,
`--boost-archive` and `--crlibm-archive` add provenance blocks and are optional;
leaving them out changes only the provenance section of the output, never a
count. On Windows use Git Bash; the same block works on Linux with the tree
paths substituted.

The recorded run prints:

```
wrote counts.json
Boost 109200: 522 test files, 499 bounds (plain variant 499)
crlibm 1.0beta5: 63617 vectors over 21 files; Gappa 791 arrow lines over 46 files; source comments 41
registry: 12 failures (4 sinh + 8 atanpi), 44 undecided, 361286 queries; per-query results deposited: False
comparison: 8 engine-recounted (8 agree), 12 aggregate-only (12 agree), disagreements: none
```

The deposited `counts.json` is the output of that command run from this
directory. A fresh filesystem copy of `accuracy-claim-audit` without `.git` or
`__pycache__` was then made at `C:/Users/myaml/dev/_recount-clone`, and the same
command was run from its `recount/` directory: it exited 0 and rewrote
`counts.json` with SHA-256
`d6dc88d9227fbf289e2f1ae8b1892f0d2c7bb0127571332cde18e4cfdc7e19d3`, byte-identical
to the deposited file and deeply equal to it with no field ignored. The
cross-check below also exited 0 and rewrote `expected_bounds.csv` with SHA-256
`9028cd07017be2682b24a704c07f66ff471bd4d0899868d83ee570f3ca1c52c2`, also
byte-identical.

The 499 Boost records can also be produced by the audit's own tool:

```sh
python -I ../boost_accuracy_rows/tools/expected_results.py \
  C:/Users/myaml/dev/_boost-smoke/boost-1_92-src/boost_1_92_0/libs/math/test \
  expected_bounds.csv
```

It prints `files scanned : 522` and `bound records : 499`. The crlibm vector
count can be checked with `../claim-audit-crlibm/work/parse_testdata.py`; calling
its `counts(path)` on the pinned `tests/` directory totals 63,617 over 21 files.

## Counting rules

**B — Boost test-suite bounds.** In `libs/math/test/test_*.cpp`, one bound is one
textual `add_expected_result(...)` call with at least eight top-level arguments.
Comments are removed; preprocessor branches inside an argument list are resolved
with `BOOST_MATH_ENABLE_SYCL` undefined and other macros treated as undefined,
the configuration the audit reports. A variant that keeps both branches and only
drops directive lines gives the same 499.

**V — crlibm vector assertions.** In `tests/*.testdata`, one vector is one
non-comment row whose first token is `N`, `P`, `M` or `Z` and which carries four
hexadecimal words (six for `pow`). A line whose first token is not a mode is a
title or header and is skipped; a malformed mode line would be reported. One
file is one function.

**G — crlibm Gappa obligations.** One obligation is one line containing `->` in
`gappa/**/*.gappa`. This is the rule that reproduces 791. The alternatives are
reported in `counts.json`: 791 arrow occurrences, 783 non-comment arrow lines,
707 non-comment arrow lines that end in `;`, 46 files, and 6 scripted `$GAPPA`
invocations in two shell scripts. A Gappa obligation in the strict sense is a
rewrite rule or hint; 791 is the arrow-line count, not the number of goals Gappa
would report.

**C — crlibm numeric bounds in source comments.** In top-level `*.c` and `*.h`
(so `tests/`, `gappa/`, `maple/`, `docs/` and `scs_lib/` are out), one bound is
one physical comment line containing both a power-of-two literal (`2^n`,
`2^(n)`, `2**(-n)`) and one of the words error, accurate, accuracy, precision or
bound. This is the rule that reproduces 41; it is a line count.

**Registry.** The 12 failures and their split, the 44 undecided points, the
361,286 queries and the 1,445,144 rounding-mode results are read from
`registry/counts.json`, not recounted. `m2_claims.py` writes a per-query
`m2_results.json`; it is not deposited, and neither is `m1_results.tsv`. The
audit's own README states that the per-query output was not saved.

## Results

Source-recounted:

| Paper number | Paper | Recount | Rule | Agrees |
| --- | ---: | ---: | --- | --- |
| Boost.Math suite bounds | 499 | 499 | B | yes |
| Boost.Math test files | 522 | 522 | B | yes |
| crlibm vector assertions | 63,617 | 63,617 | V | yes |
| crlibm functions | 21 | 21 | V | yes |
| crlibm Gappa obligations | 791 | 791 | G | yes |
| crlibm Gappa files | 46 | 46 | G | yes |
| crlibm source-comment bounds | 41 | 41 | C | yes, rule-sensitive |

Aggregate registry only, no per-query rows:

| Paper number | Paper | Registry | Agrees |
| --- | ---: | ---: | --- |
| correct-rounding failures | 12 | 12 | yes |
| — sinh | 4 | 4 | yes |
| — atanpi | 8 | 8 | yes |
| undecided reference points | 44 | 44 | yes |
| queries | 361,286 | 361,286 | yes |
| rounding-mode results | 1,445,144 | 1,445,144 | yes |
| Table S9 oracle | 0 in 5,999, 72 overflow-boundary | same | yes |
| Table S9 own vectors | about 48,000, 0 disagreements | (~48,000), 0 | yes |
| Table S9 vector count | 63,617, 9 disagreements | 63,617, 9 | yes |

The registry also holds an unrelated 44
(`documented_claims.internal_numeric_bound_sentences`); the paper's 44 undecided
points is the `measurements.m2_crlibm_vs_independent_reference` field.

## Caveats

`791` and `41` have no deposited counter; both appear only in
`registry/counts.json`. The rules above reproduce them exactly, but `41` is not
uniquely determined by the record:

- broad line rule (the rule above): **41**
- accuracy word only, dropping the word "bound": **38**
- value bound only, also dropping the 7 "We must prove extra accuracy in the
  interval 0..sin(2^-18)" lines, whose number is a domain endpoint: **31**

The three-line difference between 41 and 38 is `pow.c:581,584,641`, where an
integer exponent is "bounded by 2^11": a magnitude bound, not an accuracy bound.
A reader applying "numeric accuracy bounds" strictly gets 31–38, not 41. The
number holds under the broad rule; the record does not say which rule was used.

For `791`, the alternative semantic counts (783 non-comment arrows, 707 arrows
that terminate the rule on the same line, 46 files, 6 scripted invocations) are
all in `counts.json`. The paper describes the count as "counted, not re-run",
which is consistent with the arrow-line reading.

The 12-failure split, the 44 undecided points and every Table S9 total match
`registry/counts.json`, but they are not independently recounted here. The
registry places the 8 atanpi failures above about 3e16 up to 1e300, outside the
manual's proven interval [9.4e-8, 0.098], and the 4 sinh failures at
`0x408633ce8fb9f87e` and its negative mirror. The per-query rows that would let
a reader count them are not in the deposited record.

The `registry` block of `counts.json` records the local registry copy
(`registry/counts.json`) and what sits beside it. The audit's scripts are one
directory up, in `../claim-audit-crlibm/work/`; because `m2_claims.py` is not
copied beside the registry, the `m2_script_writes_m2_results` check reports
null. The script is present one directory up and does write `m2_results.json`;
that file, and `m1_results.tsv`, are not in the deposit, which is why
`per_query_results_deposited` is false.

## Scope

The script reads only. It does not write into either source tree, does not
build or execute anything from them, and does not run git. No paper or
repository file is edited by the recount.
