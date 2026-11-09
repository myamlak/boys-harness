# R4 and the 9,353rd first-sweep argument

Two Paper A input lists were not carried as argument lists at harness
`v1.0.34`: the 73-argument R4 committed-grid subset, and the 9,353rd argument
of the first amplification sweep. This record reconstructs both from the
deposited files at the submission tags, verifies them against the deposited
counts and summary values, and deposits the lists with their digests.

## R4: the committed 45-digit grid subset

The R4 check compares the paper's 60-digit reference with the kernel's
committed 45-digit grid, `tests/data/boys_reference.csv` at boys `v5.3.1`
(sha256 `e25cf343be108c21842f5590aa98cf73d32faa6bf772862cd8afcf39772f17e8`).

The committed file carries 88 arguments at 33 orders (2,904 cells). The R4
subset is the 73 arguments below the region-A upper edge,
`kRegionAEnd = 11.899848152108484`; the remaining 15 arguments lie above it
and are outside the evaluated bands.

| Quantity | Value |
| --- | ---: |
| committed arguments | 88 |
| R4 arguments (x < `kRegionAEnd`) | 73 |
| orders per argument | 33 |
| R4 cells | 2,409 |
| R4 arguments present in the first sweep's 9,353-argument grid | 73 of 73 |
| R4 arguments present in the second sweep's 836-argument grid | 73 of 73 |

Verification: over the 2,409 cells, the largest absolute difference between
the deposited 60-digit reference and the committed 45-digit value is
`4.975492884e-46` at `(n=2, x=3.162277660168379e-06)`, which is the deposited
summary's `4.975493e-46` at that cell. This reproduces the R4 row from the
list itself and confirms that the 73 arguments are the grid of the check.

Files: `r4-arguments.csv` (73 rows: band, x, C99 hexadecimal spelling, orders,
cells) and `r4-cells.csv` (2,409 rows: band, n, x, the committed 45-digit
value). Band labels are taken from the deposited second-sweep grid.

## The 9,353rd first-sweep argument

The designed first-sweep grid has 9,352 arguments. The paper's first
experiment has 9,353: the last row of the 9,353-argument record is a recovered
B2 argument that no printed grid rule generates.

| Quantity | Value |
| --- | ---: |
| designed grid rows (harness `product-amplification/data/grid.csv`) | 9,352 |
| first-sweep record rows (paper tree and `first-sweep-9353` deposit) | 9,353 |
| relation | the 9,352-row designed record is the 9,353 grid's first 9,352 rows |
| recovered argument | `x = 5.949924076054243 = 0x1.7ccb8e5a1580ep+2` |
| rule | `nextafter(kRegionA1Edge, +inf)`, `kRegionA1Edge = 5.949924076054242` |
| band and tag | A2, `recovered` |
| reference and cells | 33 reference values and 33 cells present in the deposits |

Verification: the last row of the deposited 9,353-argument grid is the
recovered argument; it is exactly `nextafter(kRegionA1Edge, +inf)`; its C99
hexadecimal spelling is `0x1.7ccb8e5a1580ep+2`; no other row carries the
`recovered` tag; the designed 9,352-row record is its prefix.

Files: `argument-9353.csv` (one row with the position, band, both spellings,
tag and rule) and `verification.txt` (the run's output).

## Invocation

From this directory:

    python extract_and_verify.py <boys-v5.3.1-tree> <paper-tree> <harness-root>

with `<boys-v5.3.1-tree>` a `git archive` of boys tag `v5.3.1` (commit
`2f83e40a6504368b7b2f13682e29b4f24ab2a40c`), `<paper-tree>` the
`papers-a-v5/papers/boys-crossover` directory, and `<harness-root>` the
harness's `amplification-sweeps` directory. Exit 0; `run.log` is the captured
output.

## Digests

| Item | sha256 |
| --- | --- |
| committed 45-digit grid at v5.3.1 | `e25cf343be108c21842f5590aa98cf73d32faa6bf772862cd8afcf39772f17e8` |
| second sweep grid (836) | `5fee3f7930f52e1b982c421e84d2b6e36d97af450b8512508cd9ea5dc7b9f507` |
| second sweep 60-digit reference | `514411570e902bc2320932a8ae6200d46e11d950471994c80c2775871f9a528d` |
| first sweep grid (9,353) | `780f75744e71f34da9e42b22ef359cc56424283e67a4faf7a16ea58847fae488` |
| first sweep designed grid (9,352) | `da7a6c12090b0367054f561fdce749fb8147879dbae3a2544902d5b51667414c` |
| `r4-arguments.csv` | `d8c1e1f431596bd2b8f5a311664ccbba1b0133be8ea1add203f539a5fde2e8fd` |
| `r4-cells.csv` | `bb25a9941c2db265b2d3b4dd0bb71c25506fa0d66560a1e29239b6cb5fb3b596` |
| `argument-9353.csv` | `f4a8e10b1f3534a0142c6ea8950b004ee57c9e73a6747ebf7690480d0b0d8224` |
| `extract_and_verify.py` | `4add0ad6fc3c6130fde303b0664dfbeb91c47b4ddddb984385e9cdba3992de23` |

`digests.txt` is the generated digest list; `provenance.json` is the
machine-readable record.

## Result

Both lists reconstructed and verified: R4 is 73 arguments and 2,409 cells,
reproducing the deposited worst disagreement at the deposited cell; the
9,353rd argument is the recovered `nextafter(kRegionA1Edge, +inf)` row of the
first-sweep record.
