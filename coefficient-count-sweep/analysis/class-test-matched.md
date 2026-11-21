# The matched-solver class test

The alpha-independence test compares the fitted width-law slope of `F_1` with
that of `F_0` at the same edge and target. The comparison is clean only when
both functions' counts come from the same solver. The earlier mixed variant
used the independent solver for `F_1` and the kernel record for `F_0`; it
compares solvers and is not an alpha test.

This record builds the matched comparison:

- `F_1` at `x0` and `a2`: `results-f1-independent/run_f1_x0.json` and
  `run_f1_a2.json`.
- `F_0` at `x0`: `results-corrected/corrected-counts.json`.
- `F_0` at `a2`: `results-second-edge-a2/analysis.json`.

The statistic is the paired slope difference over the thirteen scaling widths.
The script writes `results-f1-independent/class_test_matched.json`.

## Result

Five of the eight relative blocks sit above two paired standard errors. The
rational differences run from `+0.00977` to `+0.02033` per octave; the
Chebyshev blocks are kernel-to-kernel and unchanged. The absolute arm is not
resolved.

The mixed variant is carried as a sensitivity in the same JSON
(`mixed_sensitivity`). It gives the earlier differences and paired `z`; because
it compares the independent solver with the kernel, it is not used as the
alpha test.

## Command

From `coefficient-count-sweep`:

```sh
python analysis/class_test_matched.py
```

The script uses only the Python standard library. It reads the deposited
records beside it and writes the JSON in `results-f1-independent/`.
