# Clean-edge coefficient-count ladder

This directory records the eight-target coefficient-count ladder at the clean
left edge `a = 62`. It asks whether the rational family's departure from the
affine target law, measured at the x0 edge in the coefficient-count paper,
persists where the analytic transient is below every target.

The record holds 40 rows: five widths at eight relative targets from `1e-8` to
`1e-15`. The pre-registration and its amendment are in `PRE-REGISTRATION.md`.
The analysis is in `edge62_analysis.json`; its verdict is stated below.

## Design and scale

The x0 edge is `11.899848152108484`, and `a = 62` is `5.210150516837854` times
it. Every width is multiplied by `a/x0 = 5.210150516837854`, so `W/a` matches
the x0 ladder's design and each fit is measured over the same relative interval
at the new edge. The transient factor `1 - g_1` is below every target, so the
edge is clean at the smallest target, `1e-15`.

| Width | Base width at x0 | Width at `a = 62` | Right endpoint |
| --- | --- | --- | --- |
| `2^-7` | 0.0078125 | 0.040704300912795735 | 62.040704300912796 |
| `2^0` | 1 | 5.210150516837854 | 67.210150516837854 |
| `2^5` | 32 | 166.72481653881133 | 228.72481653881133 |
| `ship` | 19.100151847891516 | 99.514666021973475 | 161.51466602197347 |
| `paper` | 17.089489586712256 | 89.038813002703985 | 151.03881300270398 |

The eight targets are the relative targets `1e-8`, `1e-9`, `1e-10`, `1e-11`,
`1e-12`, `1e-13`, `1e-14` and `1e-15`. The ladder is the full `5 x 8` grid.
The driver writes a row after each completion and reuses finished rows on a
re-run; the completing run ended with `done: 40 rows`. `edge62_run2.log` is its
console transcript, one line per completed row plus the final count.

## Solver settings

The sweep uses the recorded harness revision deposited one level up as
`../count_sweep_recorded.py`, unchanged, and the `v5.0.0` kernel's coefficient
generator. Every row records:

- binary64 (`fp64`) evaluation on a 2,001-node certifying grid (2,000
  intervals);
- coarse stride 5;
- rational scan outer 100, rational refine outer 300 and 60 Remez iterations;
- Chebyshev degree mode `all`;
- rational denominator degree cap `k <= 12`;
- 12 workers.

The rational solver in each row is Remez or Lawson, as the row records. A
stored count is an existence bound under this solver, not a proof of
minimality.

## Pre-registered rule

For each width the analysis fits the stored count against the target exponent
`-log10(target)` on a free line and reports the maximum absolute residual and
the maximum absolute second difference. Two fits are reported:

- **Primary:** all eight targets, `1e-8` to `1e-15`. This fit chooses the
  verdict.
- **Secondary:** the paper's contract range, `1e-8` to `1e-13`, with the same
  thresholds. This fit is a diagnostic and does not choose the verdict.

The rule is:

- **Chebyshev control.** The affine target law holds at a width when the
  free-line residual is at most one stored coefficient. Chebyshev must hold at
  every width for the clean-edge test to be interpretable.
- **Rational failure at a width.** The residual exceeds the two-coefficient
  rounding band (greater than 2), or the absolute second difference is at
  least 2.
- **Clean-edge verdict.** *Failure confirmed at the clean edge* if rational
  fails at at least three of the five widths while Chebyshev holds at every
  width; *transient-specific* if rational holds at at least four of the five
  widths; *mixed* otherwise.

The pre-registration was written at 37 of 40 rows. Its amendment, made at 38 of
40 rows before the last two rows were read, made the all-eight-target fit
primary and the contract-range fit secondary; the thresholds and the three-way
verdict rule are unchanged.

## Verdict

**Primary, all eight targets: mixed.** Rational fails at `2^5` and `ship`, and
holds at `2^-7`, `2^0` and `paper`; Chebyshev holds at every width.

**Secondary, contract range: uninterpretable for clean-edge inference.** The
Chebyshev contract-range residual exceeds one stored coefficient at every
width, so the pre-registered interpretability condition is not met.

Primary residuals and second differences:

| Width | Chebyshev residual | Chebyshev d2 | Chebyshev | Rational residual | Rational d2 | Rational |
| --- | --- | --- | --- | --- | --- | --- |
| `2^-7` | 0.810 | 2.0 | holds | 0.464 | 1.0 | holds |
| `2^0` | 0.810 | 2.0 | holds | 0.452 | 1.0 | holds |
| `2^5` | 0.464 | 1.0 | holds | 3.500 | 6.0 | fails |
| `ship` | 0.929 | 2.0 | holds | 1.417 | 2.0 | fails |
| `paper` | 0.750 | 2.0 | holds | 0.583 | 1.0 | holds |

In the secondary fit the Chebyshev residuals are 2.562, 3.857, 9.705, 8.743
and 7.857, in the width order above; all exceed one stored coefficient. The
secondary rational fit fails at `2^0`, `2^5`, `ship` and `paper` and holds at
`2^-7`, but the failed Chebyshev control makes that pattern uninterpretable for
clean-edge inference.

The clean-edge result is mixed, not a universal confirmation. The rational
failure is not a pure x0 transient artefact: it reappears at the clean edge at
`2^5` and `ship`, so it is not confined to the x0 transient regime. It is,
however, width-dependent at the clean edge, where three of the five widths
hold. The headline claim is therefore limited to what the pattern supports, and
the title keeps the scoped wording "where the affine target law fails" without
a clean-edge generalisation.

## Files

| File | Contents |
| --- | --- |
| `edge62_sweep.py` | the sweep driver |
| `edge62_sweep.json` | the 40-row measurement record |
| `edge62_run2.log` | the console transcript of the completing run |
| `analyze_edge62.py` | the pre-registered analysis, standard library only |
| `edge62_analysis.json` | the analysis output |
| `PRE-REGISTRATION.md` | the pre-registration and its amendment |
| `MANIFEST.sha256` | SHA-256 of every file in this directory except itself |

## Reproducing and verifying

The analysis regenerates the deposited output from the deposited record:

```sh
python -I analyze_edge62.py
sha256sum -c MANIFEST.sha256
```

`analyze_edge62.py` reads `edge62_sweep.json` and rewrites
`edge62_analysis.json`. It imports only the standard library and runs in
isolated mode. A fresh run is byte-identical to the deposited
`edge62_analysis.json`.

Regenerating the 40 rows needs a boys `v5.0.0` checkout for
`tools/gen_boys_coefficients.py`, the recorded harness
`../count_sweep_recorded.py`, and `mpmath`. The driver takes the kernel
directory and the harness directory as arguments and writes the record after
each row; run it without `-I`, because isolated mode hides the user-site
`mpmath`.
