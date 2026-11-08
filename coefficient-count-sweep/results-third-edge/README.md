# The third edge: F0 at a left edge of 31

The control that places the left edge where the F0 transient, `erfc(sqrt(x))`, is below every
target measured here. At 31 it is 3.4e-15. At the second edge `a2` it is 5.2e-12 and at `x0`
1.1e-6, both above at least one target, so a comparison between `a2` and 31 separates the
transient's effect from the edge's.

## Records

- `a31_run_summary.json`: the rational counts at a left edge of 31, over the fifteen Group A widths
  scaled by 31/x0 = 2.605075258418927 so that `W/a` matches the x0 edge, at the three targets
  (relative 1e-13, absolute 5e-14, relative 1e-14). 45 cells, made 2026-11-01T18:28:45Z with six
  workers in 2876.1 s, under the same acceptance rule, search space (k <= 12) and solver as
  `../results-second-edge-a2/`.
- `cheb31.json`: the Chebyshev counts at the same edge and widths, from the recorded sweep harness.
  45 rows, 129.4 s. The three `2^5` rows do not converge at this argument and carry no count.

## Result

The rational counts at 31 and at `a2` agree in 43 of the 45 cells, the two exceptions one
coefficient apart. Between `a2` and 31 the transient stands at most 523 times the target at `a2`
and below every target at 31, so where it is at most that multiple of the target the stored count is
the same function of `W/a` to within one coefficient. The control is consistent with the reading
that the `x0` edge dependence is the transient; it does not locate the crossover, which lies between
5e-12 and 1e-6.

The Chebyshev counts agree with the `a2` record in all 42 cells that carry a count.

The driver that produced the Chebyshev rows imports the recorded revision `count_sweep.py`; the
records are deposited here so that the counts do not depend on a paper's local tree.
