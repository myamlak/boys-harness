# Independent rational counts for F1 and tail-type targets

Records of the independent minimax runs that give the rational counts of the second instance for
`F1(x) = (F0(x) - exp(-x)) / (2x)`, and the control runs that place the second edge where the
exponential term is far below every target.

## Solver

The numeric core is `../results-second-edge-a2/second_edge_a2_solver.py`, unchanged: baryrat 2.1.2
(BRASIL) minimax fits under the paper's acceptance rule. A fit is admissible when its binary64
Horner evaluation, on both the fused and separate multiply-add routes and on the 2001/4001-node
grids, stays inside the cell target, with the denominator guard at 1e-3 and the retry protocol and
witness selection as deposited. `rational_instance_solver.py` substitutes two things only: the target
function (F0, F1 or the pure tail `x^-1/2`) and the left edge, with the fifteen Group A widths scaled
by `a/x0` so `W/a`, and with it the cross-ratio, is matched at every edge. Every cell is solved
inside the same search space (k <= 12) and every count is the smallest stored count at which an
exhibited admissible fit meets, an existence upper bound exactly as the deposited solver's counts
are.

## Runs

Each run covers 45 cells: the fifteen Group A widths at three targets (relative 1e-13, absolute
5e-14, relative 1e-14), with max stored 30.

| Record | Function | Left edge | Cross-ratio matched | Workers | Wall (s) | UTC |
| --- | --- | --- | --- | --- | --- | --- |
| `run_f1_x0.json` | F1 | 11.899848152108484 | yes | 8 | 3577.9 | 2026-11-01T18:45:02Z |
| `run_f1_a2.json` | F1 | 23.79969630421697 | yes | 6 | 2082.3 | 2026-11-01T19:21:00Z |
| `run_f1_a62.json` | F1 | 62.0 | yes | 6 | 1364.1 | 2026-11-01T21:52:06Z |
| `run_f1_a124.json` | F1 | 124.0 | yes | 6 | 840.0 | 2026-11-01T22:06:07Z |

`x0` and `a2` are the paper's two second edges. At 62 and 124 the factor 1 - g1 of F1 is 1.06e-26
and 1.8e-53, below every target at both edges, so those two runs are the transient-free control.

## Comparison with the kernel's delivered record

`compare_counts.py --kind f1 --independent run_f1_x0.json run_f1_a2.json` reads the kernel's
corrected F1 record and prints the per-cell pair. Over the 90 cells of the two paper edges, 18
differ: 16 lower in the independent run and two higher, each by one coefficient, with the sum -18
(x0: 12 of 45 differ, sum -13; a2: 6 of 45 differ, sum -5). Where the independent count is lower, an
admissible fit is exhibited at that count, so the kernel's count is not tight in those cells.

## Control

At the two control edges the relative targets agree in all fifteen widths (15 of 15 at 1e-13 and
15 of 15 at 1e-14), and the absolute target agrees in 13 of 15, the two differing cells one
coefficient lower at the larger edge. The second-edge dependence of the F1 counts is therefore not
an artefact of the exponential transient.

## Neighbouring records

`../results-shipped-seed/` holds the check of the library's shipped rational seeds at its own bar,
and `../results-third-edge/` the F0 third-edge control at a left edge of 31.

## Reproducing

```
python rational_instance_solver.py --function f1 --left 11.899848152108484 --tag x0 \
    --max-stored 30 --workers 8 --out-dir out_f1_x0
```

`--only <cell_id>` runs one cell. Records are written under `<out-dir>/cells/`; each run summary is
the JSON beside this file, verbatim as the run wrote it. `compare_counts.py` reads the kernel
records already in this repository.

The runs need `baryrat` 2.1.2 and `mpmath` importable.
