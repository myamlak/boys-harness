# The shipped rational seeds at the library's own bar

Six re-fits that test whether the rational seed the library stores in region B holds more
coefficients than it needs. One record covers the coarsest region-B seed on `[x0, 31]`; five cover
the narrow rational pieces on the sub-intervals of that band.

## Method

The solver is `../results-second-edge-a2/second_edge_a2_solver.py`, unchanged: baryrat 2.1.2
(BRASIL) under the same acceptance rule the deposited solver applies, binary64 Horner on the fused
and separate routes, 2001/4001-node grids, denominator guard 1e-3, retry protocol and witness
selection as deposited. `shipped_seed_check.py` fixes the target to an absolute
2.56905715673095443e-14, the region-B bar the library documents, on each interval.

## Records

| Record | Left edge | Width | First meeting stored count |
| --- | --- | --- | --- |
| `coarsestB.json` | 11.899848152108484 | 19.100151847891517 | 13 |
| `narrowB_p0.json` | 11.899848152108484 | 2.8567491415463806 | 9 |
| `narrowB_p1.json` | 14.756597293654865 | 3.405194117991412 | 7 |
| `narrowB_p2.json` | 18.161791411646277 | 4.091271749574425 | 8 |
| `narrowB_p3.json` | 22.253063161220702 | 4.947466945561672 | 8 |
| `narrowB_p4.json` | 27.200530106782374 | 3.7994698932176263 | 7 |

The narrow pieces total 39 stored coefficients, 9 + 7 + 8 + 8 + 7, and the coarse seed 13. Each
record's `meets` list is the set of admissible pairs at its first meeting stored count; every row
below that count misses the bar, so no pair one stored coefficient lower meets it on any of the six
intervals. The solver therefore returns the shipped count everywhere: the seeds are not padded.

The records were made on 2026-11-01 between 20:37Z and 20:48Z. Each record carries its own
`generated_utc`, `levels_complete` and per-pair binary64 figures.

## Reproducing

```
python shipped_seed_check.py
```

The script needs `baryrat` 2.1.2 and `mpmath` importable and writes the six records beside itself.
