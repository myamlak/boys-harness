# d-second-instance: two fractional designs certified against claimgate 0.1.0

A real, deterministic second instance of the design argument, in two
sizes. Independent series approximations are summed by one joint call; the
error of that call decomposes into independent per-axis errors, so the widest
dependency cone is `w_p = 1`, a strength-2 fractional design exists, and it
certifies the joint accuracy bound at half the cells of the full factorial.

- **three-axis instance**: three series terms, full factorial 8 cells, fraction
  `OA(4,3,2,2)` 4 cells, 4 cells saved (50%).
- **four-axis instance**: the same three plus an independent fourth series
  term, full factorial 16 cells, fraction `OA(8,4,2,2)` 8 cells, 8 cells saved
  (50%).

Both certificates are **CERTIFIED**; the verdicts, coefficients, and
`claimgate` book results below are from the recorded run.

The record was produced against claimgate 0.1.0, tag `v0.1.0`
(tag object `92298bb7ad9623d560fa31d2a9d61e2c04e68a0d`, commit
`97c4105c21d60f2b0e151fb48edce309630dda6e`). The claimgate repository was
not modified and nothing was pushed.

## The shared construction

For `n` axes (`n = 3` or `4`) the joint call is

```
Y(a) = T1(a1) + T2(a2) + ... + Tn(an).
```

Each axis selects one of two truncation levels of an independent alternating
series. Each term is evaluated in `long double` and rounded to the common
fixed-point quantum `Q = 2^-40`, so the delivered term is `n_i Q` with integer
`n_i`. The joint call adds the integer numerators and scales by `Q`:
`N = n1+...+nn`, `Y = N Q`. That addition and scaling are exact in `double`
for these magnitudes.

| axis | term | series | level 0 | level 1 | cone |
| --- | --- | --- | --- | --- | --- |
| A1 | exp(-1) | `sum_{k>=0} (-1)^k / k!` | 5 terms | 8 terms | `{A1}` |
| A2 | ln 2 | `sum_{k>=1} (-1)^(k+1) / k` | 8 terms | 12 terms | `{A2}` |
| A3 | pi/4 | `sum_{k>=0} (-1)^k / (2k+1)` | 6 terms | 10 terms | `{A3}` |
| A4 (four-axis only) | cos(1) | `sum_{k>=0} (-1)^k / (2k)!` | 4 terms | 6 terms | `{A4}` |

The exact reference is the platform's `long double` value for each term and
the sum of those values for the joint call. The per-axis bound is the
alternating-series remainder bound at the coarser level, plus `Q/2` for the
fixed-point rounding, plus a 1e-15 evaluation/reference margin:

| quantity | three-axis | four-axis |
| --- | --- | --- |
| `b_A1` | 0.0083333333337890798 | 0.0083333333337890798 |
| `b_A2` | 0.11111111111156685 | 0.11111111111156685 |
| `b_A3` | 0.076923076923532674 | 0.076923076923532674 |
| `b_A4` | - | 2.4801587757334653e-5 |
| `B = sum b_i` | 0.19636752136888860 | 0.19639232295664594 |
| allowance `eta` | 1e-12 | 1e-12 |
| certified joint bound `B+eta` | 0.19636752136988861 | 0.19639232295764594 |

The per-axis errors are `m_i = delivered_i - exact_i`, and the joint error is
`e = joint delivered - joint exact = m1+...+mn` up to the floating-point floor
of forming the sums (measured largest residual 1.1102230246251565e-16 in both
instances).

**Why `w_p = 1`.** Each series evaluation and its fixed-point rounding is a
local error term that depends only on its own axis: cone `{A_i}`. The final
integer addition and power-of-two scaling have a local error that is
identically zero, so they add no axis to the cone. That is the same reading
this audit applies to an exact final store in the fp64 lane, whose operand may
depend on every live axis but which adds no axis because it has no local
error. The upstream terms already carry the operand dependence. The error
decomposition is therefore exactly the additive `m1+...+mn`, with no
first-order term whose cone has more than one axis. Hence `w_p = 1`,
`U_p = {A1,...,An}`, and the required strength is `2 w_p = 2`.

The two-level contrast basis uses the non-mean row `(1,-1)`, so `G_w = 1`.
With `w_p = 1`, `N_w = 0` and the theorem's non-additive bound is `eta`.

## Designs

Full factorial: all `2^n` cells. Fractional design: the even-parity half
fraction with defining relation `I = A1 ... An`, an `OA(2^(n-1),n,2,2)` of
strength 2. Every column is balanced and every pair of columns contains all
four level combinations equally often; the driver's strength-2 check passes
in both instances.

Three-axis fraction, 4 cells:

| cell | full cell | A1 | A2 | A3 |
| --- | --- | --- | --- | --- |
| 0 | 0 | 0 | 0 | 0 |
| 1 | 3 | 0 | 1 | 1 |
| 2 | 5 | 1 | 0 | 1 |
| 3 | 6 | 1 | 1 | 0 |

Four-axis fraction, 8 cells:

| cell | full cell | A1 | A2 | A3 | A4 |
| --- | --- | --- | --- | --- | --- |
| 0 | 0 | 0 | 0 | 0 | 0 |
| 1 | 3 | 0 | 0 | 1 | 1 |
| 2 | 5 | 0 | 1 | 0 | 1 |
| 3 | 6 | 0 | 1 | 1 | 0 |
| 4 | 9 | 1 | 0 | 0 | 1 |
| 5 | 10 | 1 | 0 | 1 | 0 |
| 6 | 12 | 1 | 1 | 0 | 0 |
| 7 | 15 | 1 | 1 | 1 | 1 |

The full matrices are in `results/three-axis/design-matrices.csv` and
`results/four-axis/design-matrices.csv`; each row carries `full_cell` so the
fraction rows can be mapped back to the full-factorial cell number.

The alias structure is complement pairs: `omega[alpha][alpha] = 1` and
`omega[alpha][alpha ^ (2^n - 1)] = 1`, every other entry zero.

- three-axis pairs: `{} <-> {A1,A2,A3}`, `{A1} <-> {A2,A3}`,
  `{A2} <-> {A1,A3}`, `{A1,A2} <-> {A3}`;
- four-axis pairs: `{} <-> {A1,A2,A3,A4}`, `{A1} <-> {A2,A3,A4}`,
  `{A2} <-> {A1,A3,A4}`, `{A1,A2} <-> {A3,A4}`, `{A3} <-> {A1,A2,A4}`,
  `{A1,A3} <-> {A2,A4}`, `{A2,A3} <-> {A1,A4}`,
  `{A1,A2,A3} <-> {A4}`.

The full `8x8` and `16x16` matrices are in the two `alias-matrix.csv` files.

## Results: three-axis instance

Full-factorial coefficients of the measured error `e` (8 cells):

| mask | support | coefficient |
| --- | --- | --- |
| 0 | `{}` | -0.078893245461092665 |
| 1 | `{A1}` | 0.0035714285713766003 |
| 2 | `{A2}` | -0.0093434343434637412 |
| 4 | `{A3}` | -0.0082241803606848407 |
| 3, 5, 6, 7 | every order >= 2 | 0 |

The four interaction coefficients are exactly zero in the run, the cone
structure's prediction: no operation's cone contains more than one axis.

Fractional estimates against the same coefficient definitions (4 cells):

| support | fractional | full factorial | difference |
| --- | --- | --- | --- |
| `{}` | -0.078893245461092665 | -0.078893245461092665 | 0 |
| `{A1}` | 0.0035714285713766003 | 0.0035714285713766003 | 0 |
| `{A2}` | -0.0093434343434637412 | -0.0093434343434637412 | 0 |
| `{A3}` | -0.0082241803606848407 | -0.0082241803606848407 | 0 |

The aliased estimates are the complementary coefficients: for example the
fractional `{A1,A2}` estimate is `-0.0082241803606848407`, the full `{A3}`
coefficient. `results/three-axis/fractional-contrasts.csv` gives the identity
`estimate = full(alpha) + full(complement)` and its residual on every row.

Certificate: **CERTIFIED**.

- per-axis bounds hold on all full and fractional cells, worst ratio
  `|m_i|/b_i = 0.85446705938018952`;
- the joint bound `B+eta` holds, worst ratio `0.50941361401711416`;
- every unsupported order is zero by the cone structure: the largest
  full-factorial interaction coefficient is 0 over 4 unsupported components,
  within `eta = 1e-12`;
- the fractional mean and main effects match the full factorial within `eta`
  (both differences 0);
- the alias matrix is the complement-pair structure; the strength-2 check
  passes; both levels of every axis are live; the measured additive residual
  1.1102230246251565e-16 is within `eta`;
- `delta`, the largest supported non-mean contrast, is 0.0093434343434637412;
  the non-additive bound is `eta` itself because `N_w = 0`.

## Results: four-axis instance

Full-factorial coefficients of the measured error `e` (16 cells):

| mask | support | coefficient |
| --- | --- | --- |
| 0 | `{}` | -0.078905510544393209 |
| 1 | `{A1}` | 0.0035714285713766003 |
| 2 | `{A2}` | -0.0093434343434637412 |
| 4 | `{A3}` | -0.0082241803606848407 |
| 8 | `{A4}` | -1.2263006738066906e-5 |
| all other masks (11 of them) | every order >= 2 | 0 |

All eleven order >= 2 coefficients are exactly zero in the run, again the
cone structure's prediction. The A4 main coefficient is small because the
difference between its two truncations is the difference between a 4-term and
a 6-term cosine partial sum.

Fractional estimates against the full-factorial coefficients (8 cells):

| support | fractional | full factorial | difference |
| --- | --- | --- | --- |
| `{}` | -0.078905510544393209 | -0.078905510544393209 | 0 |
| `{A1}` | 0.0035714285713766003 | 0.0035714285713766003 | 0 |
| `{A2}` | -0.0093434343434637412 | -0.0093434343434637412 | 0 |
| `{A3}` | -0.0082241803606848407 | -0.0082241803606848407 | 0 |
| `{A4}` | -1.2263006738066906e-5 | -1.2263006738066906e-5 | 0 |

The mean and all four main effects are exactly equal in the two designs. The
six order-2 fractional estimates are all `0`, each the sum of two order-2
coefficients that are both zero by the cone structure (for example the
`{A1,A2}` estimate is `full({A1,A2}) + full({A3,A4})`). The four order-3
estimates equal the complementary main effects, and the order-4 estimate
equals the mean. `results/four-axis/fractional-contrasts.csv` gives the
identity `estimate = full(alpha) + full(complement)` and its residual on every
one of the 16 rows.

Certificate: **CERTIFIED**.

- per-axis bounds hold on all full and fractional cells, worst ratio
  `|m_i|/b_i = 0.98897257218371304` (the tight one is A4 at its coarser
  level; the A1 ratio is 0.85446705938018952);
- the joint bound `B+eta` holds, worst ratio `0.50947417556762009`;
- every unsupported order is zero by the cone structure: the largest
  full-factorial interaction coefficient is 0 over 11 unsupported components,
  within `eta = 1e-12`;
- the fractional mean and main effects match the full factorial within `eta`
  (both differences 0);
- the alias matrix is the complement-pair structure; the strength-2 check
  passes; both levels of every axis are live; the measured additive residual
  1.1102230246251565e-16 is within `eta`;
- `delta`, the largest supported non-mean contrast, is 0.0093434343434637412;
  the non-additive bound is `eta` itself because `N_w = 0`.

## Instance-parameterised summary

`results/instances.csv` carries one row per instance. The run values are:

| quantity | three-axis | four-axis |
| --- | --- | --- |
| axes | 3 | 4 |
| full cells | 8 | 16 |
| fractional cells | 4 | 8 |
| cells saved | 4 (50%) | 8 (50%) |
| `w_p` | 1 | 1 |
| strength | 2 | 2 |
| `B` | 0.19636752136888860 | 0.19639232295664594 |
| `B+eta` | 0.19636752136988861 | 0.19639232295764594 |
| `delta` | 0.0093434343434637412 | 0.0093434343434637412 |
| largest full interaction | 0 | 0 |
| largest main difference | 0 | 0 |
| mean difference | 0 | 0 |
| worst per-axis ratio | 0.85446705938018952 | 0.98897257218371304 |
| worst joint ratio | 0.50941361401711416 | 0.50947417556762009 |
| additive residual | 1.1102230246251565e-16 | 1.1102230246251565e-16 |
| unsupported components | 4 | 11 |
| certificate | CERTIFIED | CERTIFIED |

## claimgate 0.1.0: what it judged, and what it cannot express

The per-axis accuracy claims, the joint additive-bound claim, and the
fractional-contrast claims are expressible as claimgate rows and were measured
cell by cell:

| book | three-axis | four-axis |
| --- | --- | --- |
| full factorial accuracy | 4/4 met, 32 cells, exit 0 | 5/5 met, 80 cells, exit 0 |
| fractional accuracy | 4/4 met, 16 cells, exit 0 | 5/5 met, 40 cells, exit 0 |
| fractional contrasts (mean + main effects, bound `eta`) | 4/4 met, exit 0 | 5/5 met, exit 0 |
| unsupported-order records (zero target, unjudged) | 4 records, exit 0 | 11 records, exit 0 |

The library has no axes, levels, dependency cones, `w_p`/`U_p`, orthogonal
arrays, strength or resolution, contrast or alias matrices, allowance term, or
design-saving model. All of that is computed by the instance driver; claimgate
judges only the per-cell accuracy comparisons. This is the API gap, stated
rather than hidden.

One consequence affects both records: a claim that an interaction coefficient
equals zero cannot be a *judged* claimgate claim. A cell with `ref = 0` always
has `bound >= |ref|`, so claimgate counts it as non-discriminating; a clean
judged row would be `Vacuous`, not met, even though a nonzero coefficient over
the bound would still be caught (`failures > 0`). Those rows are therefore
registered as records, and the same inequality (`|c| <= eta`) is checked
numerically in the driver's certificate. The contrast book's per-lane table
prints `-` for its worst cell because the worst cell is recorded only when the
error ratio is strictly positive, and every contrast difference here is
exactly zero; the verdicts and the `real` (domain) counts are unaffected.

**The four-axis case adds an alias consequence that the three-axis case does
not have.** The four-axis fraction aliases every order-2 coefficient with
another order-2 coefficient (`{A1,A2}` with `{A3,A4}`, and so on). Both are
zero by the cone structure, so the fractional order-2 estimates are zero, but
that zero is the sum of two coefficients, not a measurement of either one. In
the three-axis fraction an order-2 estimate aliases with a main effect, so the
estimate is visibly nonzero (the `{A1,A2}` estimate is the `{A3}` main
coefficient) and cannot be mistaken for the isolated order-2 coefficient.
claimgate has no alias term in its model, so it cannot express either
relation; in the four-axis case the zero is the more deceptive one, because a
zero-target row on the fractional order-2 estimate would be vacuous and would
not distinguish `full({A1,A2})` from `full({A3,A4})`. The driver therefore
checks every order-2 coefficient in the full factorial and uses the fraction
only for the mean and main effects. This is a design-specific consequence of
the existing gap, not a separate API feature that appeared at four axes; the
four-axis books themselves (16 cells, 5 claims, 11 records) are handled
without any further limitation.

## Reproduce

`reproduce.sh` snapshots the `v0.1.0` tag with `git archive`, builds and
tests claimgate, installs the package, builds the driver against it, and runs
`instance all`, which runs both instances. Every command and its exit code is
in `commands.txt`, written with paths relative to this directory; the
per-step output is in `logs/`, and the run output (claimgate reports for both
instances plus both certificates) is `results/run-all.log`. The final run had
20 commands and 0 failures; claimgate's own `ctest` passed 3/3, including the
failing-example test.

Requirements: Git; CMake 3.20 or newer; a C++17 compiler (the recorded run
used MSVC 14.5); and a claimgate clone carrying the public tag `v0.1.0`.
claimgate 0.1.0 is a C++17 CMake package. Set `CLAIMGATE_REPO` when the clone
is not `../claimgate` beside this repository.

```
bash reproduce.sh
```

The driver is C++17 against the installed `claimgate::claimgate` package; the
CMake project is `CMakeLists.txt` and the source is `src/instance.cpp`. On a
multi-config generator the build uses `--config Release`; `CMAKE_BUILD_TYPE`
is reported as unused by CMake, which is expected. `instance all --out DIR`
writes `DIR/three-axis/`, `DIR/four-axis/`, and `DIR/instances.csv`; the
narrower modes `three-axis`, `four-axis`, `full`, `fractional`, and `certify`
are also available.

## Files

| file | contents |
| --- | --- |
| `src/instance.cpp` | both instances, designs, coefficient/alias algebra, claimgate books, certificates |
| `CMakeLists.txt` | driver build against the installed claimgate 0.1.0 package |
| `reproduce.sh` | end-to-end rebuild and rerun; writes `commands.txt` and `logs/` |
| `commands.txt` | every command, its working directory and exit code |
| `results/instances.csv` | one parameterised row per instance: sizes, saving, bounds, ratios, verdict |
| `results/three-axis/axis-cone-table.csv` | axes, levels, cones, remainder bounds, documented bounds, references |
| `results/three-axis/design-matrices.csv` | the 8 full cells and 4 fraction cells, with `full_cell` and coded levels |
| `results/three-axis/full-cells.csv` | the 8 full-factorial cells with numerators, delivered terms, errors and residuals |
| `results/three-axis/fractional-cells.csv` | the 4 fraction cells, same columns |
| `results/three-axis/coefficients.csv` | full and fractional coefficients and their difference |
| `results/three-axis/fractional-contrasts.csv` | contrast estimates, alias partners, alias identity, residual |
| `results/three-axis/alias-matrix.csv` | the 8x8 alias matrix |
| `results/three-axis/bounds.csv` | bounds, allowance, certificate values, cell counts, saving |
| `results/three-axis/certificate.txt` | the certificate checks and verdict |
| `results/four-axis/...` | the same nine files for the 16-cell / 8-cell instance |
| `results/run-all.log` | the full run: both instances' claimgate reports and certificates |
| `logs/` | per-command output used by `commands.txt` |
