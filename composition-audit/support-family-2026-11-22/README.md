# Support-family-aware design certificates

This record deposits the design algebra behind a family-aware reading of the
composition audit's accuracy figure. The audit measures a six-axis cube of
calling choices and decomposes the signed error in an orthonormal basis. The
design question is how many cells a follow-up run needs in order to estimate
the components the path can carry, without assuming that every other component
is absent.

The record contains two exact designs (21 runs on the four-axis pathwise cube,
75 runs on the six-axis pooled cube), their file-backed certificates, the joint
contamination constant that charges the structural zeros once, the reductions
of the deposited 576-cell record that supply the measurement allowances, the
exhaustive search that shows why a smaller exactly orthogonal design does not
exist, and a five-run star design for the four-axis second instance. The design
algebra is exact over Q(sqrt 2, sqrt 3); the numerical allowances and frontier
ratios are record-specific.

`DERIVATION.md` is the corrected derivation: the theorem with the measurement
term and the (T-full) hypothesis, the joint contamination constant, the
no-replication convention, and the provenance argument for the allowance. The
scripts named in section 10 of that document are all in this directory.

## The support family

Write the cube as a product of axes with 2 or 3 levels. The basis is the
product of per-axis orthonormal rows: the mean, a contrast on a two-level axis,
and a linear and a quadratic row on a three-level axis. A component `alpha` is
a choice of one row on every axis; its support is the set of axes on which the
row is not the mean. `G_alpha` is the supremum of the component over the cube,
the product of the per-axis row suprema.

A path cone `C` is a set of axes on which a local error term can depend. The
permitted family is the union of the component sets carried by the path cones:

```
P = { alpha : support(alpha) is contained in some path cone C }.
```

`P` is not the power set of a single cone. On a three-level axis both the
linear and the quadratic row are live, so a cone on `k` axes carries the
product of the per-axis row counts, not `2^k`. The counts below use that
reading.

### Pathwise four-axis cube

Axes: `partition` (3), `pack` (2), `form` (3), `exp` (2); 36 components. Cones:

- upward step: `{partition, form, exp}` (18 components);
- packed sum: `{partition, pack}` (6 components);
- their intersection, the three `partition` rows, is counted once.

So `P` has 18 + 6 - 3 = **21 components**: the mean, six non-mean main
effects, and **14 targets** (permitted components of order 2 or 3). The
complement has **15 structural zeros**. The targets are every permitted
interaction channel, so the (T-full) hypothesis holds with `T` the fourteen
targets.

### Pooled six-axis cube

Axes: `route` (2), `scheme` (2), `partition` (3), `pack` (2), `form` (3),
`exp` (2); 144 components. Cones:

- upward step: `{route, scheme, partition, form, exp}` (72 components);
- packed sum: `{partition, pack}` (6 components);
- their intersection, the three `partition` rows, is counted once.

So `P` has 72 + 6 - 3 = **75 components**: the mean, eight non-mean main
effects, and **66 targets** (permitted components of order 2 to 5). The
complement has **69 structural zeros**. Again `T` is every permitted
interaction channel.

The components outside `P` are structural zeros: their support is contained in
neither cone, and the path decomposition gives them no first-order local error
term. They are not assumed to have zero coefficient. Their measured coefficient vector is allowed for by the
structural-zero allowance

```
eta_Q = sum over beta in Q of G_beta |c_beta|,
```

where `c` is the full-cube coefficient vector. The allowance is measured on
the full cube that the smaller design would replace; section "What the
allowance means" states what that implies.

## The certificate theorem

Let `D` be a design, a set of `N` cells of the cube. Write `phi_P(c)` for the
vector of permitted basis components at cell `c`, and define

```
Omega_PP = (1/N) sum_{c in D} phi_P(c) phi_P(c)^T,
Omega_PQ = (1/N) sum_{c in D} phi_P(c) phi_Q(c)^T,
B        = Omega_PP^{-1} Omega_PQ.
```

The permitted block is estimable when `Omega_PP` is invertible. The theorem
assumes:

1. `Omega_PP` is invertible;
2. the cone table is complete, so the first-order part of `c_beta` vanishes for
   `beta` in `Q`;
3. `eta_Q` bounds the weighted coefficients outside `P`,
   `eta_Q >= sum_{beta in Q} G_beta |c_beta(rho)|`;
4. the adjusted measurement noise `nu_hat_P = Omega_PP^{-1} nu_P` is bounded by
   `delta_alpha >= |nu_hat_alpha|` for `alpha` in `T`, where `nu` is the
   measurement residual;
5. (T-full): `T` is every permitted interaction channel,
   `T = { alpha in P : |support(alpha)| >= 2 }`. If `T` were a strict subset,
   the omitted permitted interactions would have to be added to `T` or bounded
   separately.

With the design estimate

```
chat_P = Omega_PP^{-1} (1/N) sum_{c in D} r(c) phi_P(c),
```

the exact alias identity is

```
chat_P = c_P + B c_Q + nu_hat_P.
```

For a target `alpha`, set

```
kappa_alpha = max over beta in Q of |B_{alpha beta}| / G_beta.
```

Then the componentwise bound is

```
|c_alpha| <= |chat_alpha| + delta_alpha + kappa_alpha eta_Q.
```

The earlier per-target sum of suprema charged each target separately. The
joint statement replaces that sum by one constant. Define

```
kappa_joint = max over beta in Q of
              sup over a in the cube of
              |Phi_beta(a) - sum_{alpha in T} B_{alpha beta} Phi_alpha(a)| / G_beta.
```

Then for every cell `a` of the cube the full non-additive part

```
N(a) = sum_{alpha in T} c_alpha Phi_alpha(a)
       + sum_{beta in Q} c_beta Phi_beta(a)
```

obeys

```
|N(a)|
    <= sum_{alpha in T} G_alpha |chat_alpha|
       + kappa_joint eta_Q
       + sum_{alpha in T} G_alpha delta_alpha.
```

The bound covers the permitted interactions and the structural zeros together,
which is the quantity the composition rule needs.

The charge `kappa_joint eta_Q` is paid once, not once per target. The proof
uses the identity

```
sum_T c_alpha Phi_alpha
    = sum_T chat_alpha Phi_alpha + sum_Q c_beta r_beta
      - sum_T nu_hat_alpha Phi_alpha,
```

with `r_beta(a) = Phi_beta(a) - sum_T B_{alpha beta} Phi_alpha(a)`, and the
definition of `kappa_joint`. The first display now carries the measurement
term the proof produces; a display without `delta_alpha` is the zero-noise
case and is false under assumption 4 with nonzero noise.

**Orthogonal-array special case.** If `D` is an orthogonal array of strength
`t` with `2 w_p <= t <= |U_p|`, then `Omega_PP = I` and `Omega_PQ = 0` on the
surviving supports. With `B = 0`, `r_beta = Phi_beta` and
`kappa_joint = 1`, and the theorem becomes

```
|c_alpha| <= |chat_alpha| + delta + G_alpha eta   (alpha in T),
```

which recovers the orthogonal-array theorem. The full factorial gives the exact
special case. The recovery needs a qualification: `eta_Q` bounds every
component outside the family, including order-two and order-three structural
zeros, whereas the orthogonal-array allowance bounds only orders above `w_p`.
On an orthogonal array those low-order zeros do not alias into the targets
(their unions are at most `2 w_p <= t`), so the orthogonal-array theorem never
needs them in the allowance. The single-allowance bound `kappa_alpha eta_Q` is
therefore in general looser than `delta + G_w eta` on an orthogonal array.
Exact recovery holds either with an alias-support-restricted allowance or at
the proof's intermediate step, where the sum over surviving aliases is bounded
by `G_alpha eta`. The theorem is a strict generalisation of the design
condition, not of the packaged allowance.

**No-replication convention.** "Within the 36-cell cube the only strength-2
array is the full cube" holds for a simple design, in which each cell is used
at most once. With replication, other 36-run strength-2 arrays exist; they use
36 runs and do not lower the minimum. The minimality statement is unaffected.

## The 21-run pathwise design D21

The permitted family has 21 components, so at least 21 runs are needed. Cell
order is `partition, pack, form, exp`, with index
`partition*12 + pack*6 + form*2 + exp`. D21 is:

- the full 3 x 3 `partition x form` grid at `(pack, exp) = (0, 0)`;
- the full 3 x 3 `partition x form` grid at `(pack, exp) = (1, 1)`;
- the three cells at `(pack, exp) = (0, 1)` with `form = 0` and
  `partition = 0, 1, 2`.

The cells are

```
0, 1, 2, 4, 7, 9, 11, 12, 13, 14, 16, 19, 21, 23, 24, 25, 26, 28, 31, 33, 35.
```

Certificate (`certify_d21.py`, `certify_d21.json`, `certify_d21.log`):

- exact rank of the 21 x 21 permitted block over Q(sqrt 2, sqrt 3): 21 of 21;
- `kappa_max = sqrt(2/3) = 0.816496580927726`;
- all 14 target values `G_alpha kappa_alpha` equal 1, so the per-target sum is
  14;
- Gram condition 12, so the model-matrix condition is
  `sqrt(12) = 2 sqrt 3 = 3.464101615`; threshold amplification
  `max_alpha sum_gamma |(Omega_PP^{-1})_{alpha gamma}| G_gamma = 14.024284`.

The equality `G_alpha kappa_alpha = 1` is the design's balance property: the
alias amplification exactly offsets the component's supremum, so every target
carries the same per-target alias term. The joint constant is smaller:

- `kappa_joint(D21) = 4` exactly; the per-beta ratios take the values 1, 2, 3
  and 4, with the maximum 4 attained on the two channels with support
  `{partition, pack, form}` (with and without `exp`; the full breakdown is in
  `kappa_joint.json`);
- measurement factors: max over targets 8.485281, target sum 74.111111;
- taking the residual scale as 1e-29 (fp64) or 1e-22 (fp32), the propagated
  measurement term is at most 2.5e-27 or 2.5e-20, negligible against the bound.

`kappa_joint.py`, `kappa_joint.json` and `kappa_joint.log` contain the
per-beta breakdown. A different valid 21-cell design has the same target
`kappa_max` but per-target sum 11.25 and a worse amplification of 20.830; it is
not the design cited.

## The 75-run pooled design D75

Cell order is `route, scheme, partition, pack, form, exp`, with index
`route*72 + scheme*36 + partition*12 + pack*6 + form*2 + exp`. D75 is:

- all 72 cells with `pack = arguments` (the full five-axis factorial on
  `route x scheme x partition x form x exp`);
- the three cells with `pack = orders` at `partition = 0, 1, 2` and
  `route = scheme = form = exp = 0`.

The three added cells are 6, 18 and 30; the full cell list is in
`certify_d75.json`.

Rank proof. The 72 pack-mean columns restricted to the 72 `pack = arguments`
rows are the five-axis orthonormal basis and are invertible. Subtracting the
matching pack-mean column from each of the three pack-contrast columns leaves a
matrix supported on the three `pack = orders` rows with the 3 x 3 matrix
`[1, linear(p), quadratic(p)]` at `p = 0, 1, 2`; its determinant is
`24 sqrt 3 = 41.569219` in the unnormalised rows (`3 sqrt 3` with the
normalised rows). The rank is also checked numerically and modulo the primes
1009, 2017 and 3001.

Minimality excludes 73 and 74 runs. The 72 pack-mean columns are functions of
the five-axis projection, so 72 distinct projections are necessary. At
`N = 72` each projection appears once and the pack-contrast columns lie in the
pack-mean span. At `N = 73` or `N = 74` the pack-contrast residuals after
regression are supported on one or two duplicated pairs and have rank at most
1 or 2, less than 3. Only at `N = 75`, with three duplicated pairs at distinct
partition levels, is the residual 3 x 3 block invertible.

Certificate (`certify_d75.py`, `certify_d75.json`, `certify_d75.log`):

- rank of the 75 x 75 permitted block: 75 of 75 numerically, 75 of 75 modulo
  each of 1009, 2017 and 3001, and exactly by the block proof;
- `kappa_max = 1.0` (to floating-point precision);
- all 66 target values `G_alpha kappa_alpha` equal 1 within 1.3e-15, so the
  per-target sum is 66;
- Gram condition 48, so the model-matrix condition is
  `sqrt(48) = 4 sqrt 3 = 6.928203230`; threshold amplification
  `max_alpha sum_gamma |(Omega_PP^{-1})_{alpha gamma}| G_gamma = 54.506148`;
- `kappa_joint(D75) = 3` exactly;
- measurement factors: max over targets 33.116168, target sum 253.722222.

`kappa_joint.py`, `kappa_joint.json` and `kappa_joint.log` contain the
per-beta breakdown.

## The exact obstruction

`exact_proofs.py` and `enumerate_designs.py` establish why the designs are not
smaller.

- A strength-2 orthogonal array on levels `(3, 2, 3, 2)` has every pair
  marginal balanced, so `N` is divisible by 9, 6 and 4, hence by
  `lcm(9, 6, 4) = 36`. Within the cube the only simple strength-2 array is the
  full cube; strength 3 and above implies strength 2. With the no-replication
  convention, no proper orthogonal array of any strength exists below 36 runs.
- Condition A, the exact-orthogonality condition of the theorem (target-target
  inner products diagonal, target-main-effect inner products zero, and the
  normalization), has a 36-row constraint matrix of exact rank 36 over
  Q(sqrt 2, sqrt 3). With non-negative weights and the normalization, the only
  solution is the uniform full cube. So there is no 0/1 exact-orthogonality
  design below 36 runs.
- The exhaustive 0/1 search over all designs satisfying
  `Omega_{T,R0} = 0` finds exactly three solutions: the two 18-cell designs and
  the 36-cell cube. The two 18-cell designs have target rank 12, not 14. Their
  null relations are the pack/exp alias relations:

```
V0:  [partition,exp](1,0,0,1) + [partition,pack](1,1,0,0) = 0,
     [partition,exp](2,0,0,1) + [partition,pack](2,1,0,0) = 0;

V6: -[partition,exp](1,0,0,1) + [partition,pack](1,1,0,0) = 0,
    -[partition,exp](2,0,0,1) + [partition,pack](2,1,0,0) = 0.
```

The 21-run and 75-run designs trade exact orthogonality for an invertible
permitted block and a measured alias term, charged once through
`kappa_joint`. That is the whole content of the family-aware certificate.

## Running the design on the deposited record

`pooled_d75_record_check.py` applies D75 to all 858 points of each deposited
arm in the fp64 and fp32 resolving lanes. For every point it computes the
full-cube coefficient vector, the design estimate, the alias identity residual
`chat_P - c_P - B c_Q`, the target estimate, and the alias-only bound
`kappa_alpha eta_Q`.

| arm | lane | max abs(chat - c) | identity residual | max target estimate | points above alias-only bound |
| --- | --- | --- | --- | --- | --- |
| fused | fp64 | 5.905906e-17 | 6.13e-29 | 1.329474e-14 | 857 / 858 |
| fused | fp32 | 5.701968e-08 | 3.90e-22 | 6.584066e-08 | 384 / 858 |
| separate | fp64 | 8.927247e-17 | 1.09e-28 | 1.331010e-14 | 856 / 858 |
| separate | fp32 | 5.375867e-08 | 3.90e-22 | 7.133651e-08 | 386 / 858 |

The identity residual is at the double-precision floor: the alias identity is
exact algebra, so this is the check that the reduction is wired correctly. The
target estimate exceeds the alias-only bound wherever the permitted signal or
the measurement term is larger than the alias allowance; the bound in that
column is the alias term alone, not the full certified bound. The point counts
are close to 858 in fp64 and about 45 per cent in fp32 because the fp32 alias
allowance is larger than the other terms at the remaining points. On fp64 the
target estimates are about 34 times the allowance, so the 75-run design rejects
additivity at the same 857/858 and 856/858 points as the full cube.

The structural-zero allowances that enter the frontier are:

| cube | lane | max eta_Q |
| --- | --- | --- |
| pathwise | fp64 | 2.770933e-16 (both arms) |
| pathwise | fp32 | 1.710590e-07 (both arms) |
| pooled | fp64 | 3.954108e-16 fused, 3.954120e-16 separate |
| pooled | fp32 | 2.375820e-07 fused, 2.177113e-07 separate |

They come from `pathwise_summary.json` and `pooled_allowance.json`.

## Corrected frontier

The family-aware bound is

```
family bound = figure + M_hat + kappa_joint eta_Q,
```

where `M_hat = max_a |sum_{alpha in T} chat_alpha Phi_alpha(a)|` is the
measured reconstruction of the permitted interaction signal and the figure is
the record's published fp64 or fp32 accuracy figure. The full-cube reference is

```
full bound = figure + N_full,
N_full = max_a |sum_{alpha in T} c_alpha Phi_alpha(a)
                + sum_{beta in Q} c_beta Phi_beta(a)|,
```

the figure plus the full cube's measured non-additive part. The figures used
are 2.558314305293831e-14 (fp64) and 1.5e-07 (fp32). Taking the residual scale
as 1e-29 (fp64) or 1e-22 (fp32), the measurement term is below 2.5e-27 or
2.5e-20 and is omitted from the table.

| design / lane | full non-add. | full bound | eta_Q max | kappa_joint | charge | M_hat | family bound | ratio |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| fused fp64 path | 2.9895e-16 | 2.5882e-14 | 2.7709e-16 | 4 | 1.1084e-15 | 2.9895e-16 | 2.6990e-14 | 1.0428 |
| separate fp64 path | 3.4467e-16 | 2.5928e-14 | 2.7709e-16 | 4 | 1.1084e-15 | 3.3367e-16 | 2.6990e-14 | 1.0410 |
| fused fp64 pooled | 2.0975e-14 | 4.6558e-14 | 3.9541e-16 | 3 | 1.1862e-15 | 2.0975e-14 | 4.6558e-14 | 1.000000 |
| separate fp64 pooled | 2.0903e-14 | 4.6486e-14 | 3.9541e-16 | 3 | 1.1862e-15 | 2.0903e-14 | 4.6486e-14 | 1.000000 |
| fused fp32 path | 2.0907e-07 | 3.5907e-07 | 1.7106e-07 | 4 | 6.8424e-07 | 2.4709e-07 | 1.0813e-06 | 3.0114 |
| separate fp32 path | 1.7106e-07 | 3.2106e-07 | 1.7106e-07 | 4 | 6.8424e-07 | 2.2808e-07 | 1.0623e-06 | 3.3088 |
| fused fp32 pooled | 2.2213e-07 | 3.7213e-07 | 2.3758e-07 | 3 | 7.1275e-07 | 2.6609e-07 | 1.1288e-06 | 3.0335 |
| separate fp32 pooled | 2.1456e-07 | 3.6456e-07 | 2.1771e-07 | 3 | 6.5313e-07 | 2.6609e-07 | 1.0291e-06 | 2.8228 |

The `eta_Q`, `charge` and `M_hat` columns are maxima over the record points;
the family bound and ratio use the pointwise maximum of
`M_hat + kappa_joint eta_Q`. On the pooled fp64 rows the point that maximises
`M_hat` has a small pointwise allowance, so the bound coincides with the full
bound. Charging the maximum allowance at every point instead gives the
global-allowance ratio.

The corrected ratios are 1.04 on the fp64 pathwise cube, 1.000000 on the fp64
pooled cube pointwise (1.000255 on either arm with one global allowance, 0.025
per cent), and about 3 on the fp32 cube. The fp64 family-aware bound therefore
essentially matches the full factorial on the pooled cube and costs about four
per cent on the pathwise cube. The fp32 ratios are about 3: the allowance is at
the interaction-signal scale, so the reduced design does not save in fp32.

The earlier per-target sum charged 15 `eta_Q` on the pathwise cube and 67
`eta_Q` on the pooled cube. The pooled row reconciles once the full cube's own
measured non-additive part is subtracted: 67 `eta_Q` less 53.05 `eta_Q` leaves
the apparent 14-`eta_Q` increment. The joint constant removes that double
count and charges 4 `eta_Q` (pathwise) and 3 `eta_Q` (pooled) once.

The per-term breakdown is in `corrected_frontier.json`,
`corrected_frontier.log`, `reconcile_table.json` and `reconcile_table.log`.

## What the allowance means

`eta_Q` is measured on the full cube that the smaller design would replace.
The retrospective check validates the algebra and the estimator relation; it
does not by itself establish a prospective saving. A prospective use needs
either a mathematical bound on the Q components or **amortised
re-certification**: measure the full cube once, then use the smaller design for
new points or arguments under a drift condition that no error mass moves into
Q. The drift condition is the same kernel revision, arithmetic tier and region,
with the allowance refreshed whenever drift is detected.

The fp64 pathwise `eta_Q = 2.77e-16` is the same order as the permitted
interaction signal there (`1.33e-16`); the pooled `eta_Q = 3.95e-16` is about
34 times below the permitted pooled signal (`1.33e-14`). It is a measured
extraction or noise envelope, not a proved physical Taylor remainder. The fp32
`eta_Q` is at the interaction-signal scale, so the reduced design is a valid
but loose bound there and the full cube still decides the verdict.

For a prospective operated protocol with fresh measurements, the measurement
term enters through the amplification: at most 8.485281 over the targets for
D21 and 33.116168 for D75 on the design's cell average, or 14.024284 and
54.506148 in the certificate's supremum form. The 34-fold fp64 margin over the
allowance can therefore shrink. Any protocol that uses the reduced design on
new measurements must state and monitor that amplified floor; the smaller
design does not remove the amplification.

## The four-axis second-instance star

The second instance has four independent two-level axes with singleton cones
`{A1}`, ..., `{A4}`. The permitted family is the mean plus the four main
effects; the eleven order-two to order-four interactions are structural zeros.
A saturated star design estimates the five permitted coefficients in five
runs: the baseline cell (all axes at level 0) plus one cell with a single axis
at level 1, cells `0, 8, 4, 2, 1`. The even-parity half fraction needs eight
runs and the full factorial sixteen.

`second_instance_star.py` reads the recorded four-axis full-factorial
coefficients and errors, builds the star, half-fraction and full design
matrices, and exhibits the alias operator. The star has `kappa_max = 3.0` and
condition number 6.342329. With the record's allowance `eta = 1e-12`, the
certified alias contamination is `kappa_max eta = 3e-12` against the joint
bound `B = 0.19639232295664594` (`B + eta = 0.19639232295764594`). The measured
star estimate differs from the full-factorial coefficient by at most
1.39e-16, well inside the allowance. Its structural remainder is exactly zero
in that arithmetic.

## Inputs, kernel revisions and environment

The cell archives are the composition audit's own, unmodified:

| archive | sha256 |
| --- | --- |
| `composition-audit/records/cells/cells-msvc-fused.csv.xz` | `4af5694810099d07ac9eee268f20073f08d63242731d08a06ebc3d1328c0638d` |
| `composition-audit/records/cells/cells-msvc-separate.csv.xz` | `f3493c4ca6aceaf642bcb8e5ef4b4d330c1214f64e467e284c7526ed1094b5bb` |

The archives hold 494,208 rows per arm (144 cells per lane, 858 points per
lane). Their identities are the ones named in
`composition-audit/fullcube/fullcube_summary.json`.

Kernel provenance:

- measured record: Boys library tag `v5.0.0`, revision
  `482fd3388c276091834a6c12a642d9de2481c9f2`;
- source read: Boys library tag `v5.3.1`, revision
  `2f83e40a6504368b7b2f13682e29b4f24ab2a40c`.

Harness provenance: the record's inputs were deposited at harness tag
`v1.0.42` (commit `32db80209da8bdefa5d6900d1254bb4db362fae8`). The two cell
archives and `fullcube_decomposition.py` are byte-identical at `v1.0.42`,
`v1.0.43` and the current main branch.

Measured environment: Windows 11 Pro 10.0.26200, Python 3.14.7 (MSC v.1944,
64-bit), numpy 2.5.2, scipy 1.18.1, mpmath 1.4.1. No C++ compiler is needed
for this record; the cell archives are already measured. The `.xz` archives
are read with Python's standard `lzma` module.

## Commands

Run from this directory:

```
python certify_d21.py
python certify_d75.py
python kappa_joint.py
python corrected_frontier.py
python reconcile_table.py
python pooled_d75_record_check.py
```

The first two write their JSON beside the script and print the certificate.
`kappa_joint.py` writes the per-beta joint constant, `corrected_frontier.py`
writes the ratio table, and `reconcile_table.py` writes the per-term old/new
breakdown. `pooled_d75_record_check.py` reads both cell archives and writes
`pooled_d75_record_check.json`; it takes under a minute on the recorded
machine. The supporting reductions are:

```
python pooled_allowance.py
python pooled_d75.py
python pooled_d75_exact.py
python pathwise_summary.py
python exact_verify.py
python exact_proofs.py
python enumerate_designs.py
python explore.py
python general_search.py
python second_instance_star.py
```

The record-check, allowance, frontier and second-instance scripts locate the
harness repository two levels above this directory. If the record is extracted
elsewhere, set `BOYS_HARNESS_REPO` to the repository root and those scripts
then find the deposited inputs.

## The premise that matters

The packed-sum cone is `{partition, pack}` once the budget is fixed. The
packed single-precision sum's local error depends on how the sum is split
(`partition`) and on the packing decision (`pack`); after the budget is fixed
it does not depend on `route`, `scheme`, `form` or `exp`. This is the one
modelling premise that sets the 75-run count. The pathwise results do not
depend on it.

Widening the cone changes the design. If the packed-sum cone admitted, say,
`form` or `exp`, then `P` would contain more than 75 components, the permitted
block would need more than 75 independent rows, and the 75-run design would no
longer be estimable. The 75 is `|P|`, not a search artefact. The same reading
gives the pathwise cones `{partition, form, exp}` and `{partition, pack}` and
the 21-run count. Any claim that uses these run counts inherits this premise.

## Files

| file | role |
| --- | --- |
| `DERIVATION.md` | corrected derivation: theorem, joint constant, specialisation, provenance and premise |
| `instance.py` | four-axis pathwise instance: axes, cones, basis, `G_alpha`, `P`, `T`, `Q` |
| `field.py` | exact arithmetic in Q(sqrt 2, sqrt 3) |
| `exact_verify.py` | exact design evaluator (rank, `Omega_PP`, `B`, `kappa`), random screening; writes `exact_designs.json` |
| `certify_d21.py` | file-backed D21 certificate; writes `certify_d21.json` |
| `kappa_joint.py` | exact D21 and numeric D75 joint contamination constant; writes `kappa_joint.json` |
| `exact_proofs.py` | condition-A rank 36 and the 18-cell null relations |
| `enumerate_designs.py` | exhaustive 0/1 MILP search for `Omega_{T,R0} = 0` |
| `pooled_allowance.py` | six-axis pooled cones and the record's `eta_Q`, target and order-6 allowances; writes `pooled_allowance.json` |
| `pooled_d75.py` | D75 construction and candidate evaluation; writes `pooled_d75.json` |
| `pooled_d75_exact.py` | modular rank proof and alias spectrum; writes `pooled_d75_exact.json` |
| `certify_d75.py` | file-backed D75 certificate; writes `certify_d75.json` |
| `corrected_frontier.py` | corrected frontier with `kappa_joint`; writes `corrected_frontier.json` |
| `reconcile_table.py` | per-term old/new frontier breakdown; writes `reconcile_table.json` |
| `pooled_d75_record_check.py` | applies D75 to every record point in both arms and resolving lanes; writes `pooled_d75_record_check.json` |
| `pathwise_summary.py` | per-path reductions of the four-axis cube; writes `pathwise_summary.json` |
| `second_instance_star.py` | four-axis second-instance star design; writes `second_instance_star.json` |
| `explore.py` | polytope exploration for designs orthogonal to the mean and main effects |
| `general_search.py` | single-removal, pair and simulated-annealing design search |

Each `.log` file is the output of the corresponding script at the recorded
environment. `MANIFEST.sha256` covers every file in this directory except
itself.

## Limitations

- The design algebra is exact; the allowances `eta_Q` and the frontier ratios
  are measured on one deposited record, one kernel revision and one platform.
- `eta_Q` is measured on the cube being replaced. The record supports
  amortised re-certification under the drift condition above, not a
  prospective saving without that condition.
- The fp64 `eta_Q` is a measured extraction or noise envelope, not a proved
  Taylor remainder.
- The fp16 and bf16 lanes are present in the archives but are not used.
- The separate arm differs slightly from the fused arm in the fp64 target
  maximum and in the point counts; both are reported.
- The second-instance star is evaluated on the recorded full-factorial
  coefficients and errors, not on a new measurement run.
- The "points above alias-only bound" column is expected to be large because
  the bound in that check is the alias term alone. It is not a failure of the
  certified bound, which also carries `delta_alpha`.
