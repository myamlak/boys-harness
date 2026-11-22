# Support-family-aware design certificates

This record deposits the design algebra behind a family-aware reading of the
composition audit's accuracy figure. The audit measures a six-axis cube of
calling choices and decomposes the signed error in an orthonormal basis. The
design question is how many cells a follow-up run needs in order to estimate
the components the path can carry, without assuming that every other component
is absent.

The record contains two exact designs (21 runs on the four-axis pathwise cube,
75 runs on the six-axis pooled cube), their file-backed certificates, the
reductions of the deposited 576-cell record that supply the measurement
allowances, the exhaustive search that shows why a smaller exactly orthogonal
design does not exist, and a five-run star design for the four-axis second
instance. The design algebra is exact over Q(sqrt 2, sqrt 3); the numerical
allowances and frontier bounds are record-specific.

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
complement has **15 structural zeros**.

### Pooled six-axis cube

Axes: `route` (2), `scheme` (2), `partition` (3), `pack` (2), `form` (3),
`exp` (2); 144 components. Cones:

- upward step: `{route, scheme, partition, form, exp}` (72 components);
- packed sum: `{partition, pack}` (6 components);
- their intersection, the three `partition` rows, is counted once.

So `P` has 72 + 6 - 3 = **75 components**: the mean, eight non-mean main
effects, and **66 targets** (permitted components of order 2 to 5). The
complement has **69 structural zeros**.

The components outside `P` are structural zeros: the path decomposition gives
them no local error term. They are not assumed to have zero coefficient. Their
measured coefficient vector is allowed for by the structural-zero allowance

```
eta_Q = sum over beta in Q of G_beta |c_beta|,
```

where `c` is the full-cube coefficient vector.

## The certificate theorem

Let `D` be a design, a set of `N` cells of the cube. Write `phi_P(c)` for the
vector of permitted basis components at cell `c`, and define

```
Omega_PP = (1/N) sum_{c in D} phi_P(c) phi_P(c)^T,
Omega_PQ = (1/N) sum_{c in D} phi_P(c) phi_Q(c)^T,
B        = Omega_PP^{-1} Omega_PQ.
```

The permitted block is estimable when `Omega_PP` is invertible. If it is, the
design estimate of the permitted coefficients is

```
chat_P = Omega_PP^{-1} (1/N) sum_{c in D} r(c) phi_P(c),
```

and the exact alias identity is

```
chat_P = c_P + B c_Q + (measurement term).
```

For a target `alpha`, set

```
kappa_alpha = max over beta in Q of |B_{alpha beta}| / G_beta.
```

Then the certified bound is

```
|chat_alpha - c_alpha| <= delta_alpha + kappa_alpha eta_Q,
```

where `delta_alpha` is the measurement-threshold term propagated through
`Omega_PP^{-1}`. The first term is the measurement floor; the second is the
worst-case contamination from the structural zeros.

This recovers the orthogonal-array theorem as the zero-alias special case. If
`B = 0` (for example, an orthogonal array of sufficient strength), then
`kappa_alpha = 0` for every target and the bound reduces to the measurement
term alone. Here exact orthogonality is replaced by estimability plus a finite,
measured `kappa_alpha` and a measured `eta_Q`. The certificate is only as good
as the allowance; the design algebra itself is exact.

## The 21-run pathwise design D21

Cell order is `partition, pack, form, exp`, with index
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
- all 14 target values `G_alpha kappa_alpha` equal 1, and their sum is 14;
- `cond(F_P) = 3.464101615`; threshold amplification
  `max_alpha sum_gamma |(Omega_PP^{-1})_{alpha gamma}| G_gamma = 14.024284`.

The equality `G_alpha kappa_alpha = 1` is the design's balance property: the
alias amplification exactly offsets the component's supremum, so every target
carries the same certified alias term.

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
`[1, linear(p), quadratic(p)]` at `p = 0, 1, 2`, whose determinant is
`-3 sqrt 3`. The rank is also checked numerically and modulo the primes 1009,
2017 and 3001.

Certificate (`certify_d75.py`, `certify_d75.json`, `certify_d75.log`):

- rank of the 75 x 75 permitted block: 75 of 75 numerically, 75 of 75 modulo
  each of 1009, 2017 and 3001, and exactly by the block proof;
- `kappa_max = 1.0` (to floating-point precision);
- all 66 target values `G_alpha kappa_alpha` equal 1, and their sum is 66;
- `cond = 6.928203230`; threshold amplification
  `max_alpha sum_gamma |(Omega_PP^{-1})_{alpha gamma}| G_gamma = 54.506148`.

## The exact obstruction

`exact_proofs.py` and `enumerate_designs.py` establish why the designs are not
smaller.

- Condition A, the exact-orthogonality condition of the design theorem
  (target-target inner products diagonal, target-main-effect inner products
  zero, and the normalization), has a 36-row constraint matrix of exact rank
  36 over Q(sqrt 2, sqrt 3). With non-negative weights and the normalization,
  the only solution is the uniform full cube. There is no orthogonal array of
  any strength below 36 runs on levels `(3, 2, 3, 2)`, and no 0/1 design below
  36 runs that is exactly orthogonal in the theorem's sense.
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
permitted block and a measured alias term. That is the whole content of the
family-aware certificate.

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
allowance is larger than the other terms at the remaining points.

The structural-zero allowances that enter the frontier bound are:

| cube | lane | max eta_Q |
| --- | --- | --- |
| pathwise | fp64 | 2.770933e-16 (both arms) |
| pathwise | fp32 | 1.710590e-07 (both arms) |
| pooled | fp64 | 3.954108e-16 fused, 3.954120e-16 separate |
| pooled | fp32 | 2.375820e-07 fused, 2.177113e-07 separate |

They come from `pathwise_summary.json` and `pooled_allowance.json`.

## The fp64/fp32 frontier

The certified joint bound combines the design's propagated measurement
threshold with its structural-zero allowance. The full-factorial reference is
the same bound with the alias term absent. At the deposited record:

| cube | lane | design bound | full-factorial reference |
| --- | --- | --- | --- |
| pathwise (D21, 21 runs) | fp64 | 2.97e-14 | 2.59e-14 |
| pathwise (D21, 21 runs) | fp32 | 2.72e-06 | 2.81e-07 |
| pooled (D75, 75 runs) | fp64 | 7.27e-14 | 4.66e-14 |
| pooled (D75, 75 runs) | fp32 | 1.61e-05 | 3.72e-07 |

The increment over the full-factorial reference is the target count times the
lane's structural-zero allowance: 14 `eta_Q` on the four-axis cube and 66
`eta_Q` on the six-axis cube. For the pooled fp64 row the increment is
66 x 3.954108e-16 = 2.61e-14, so the bound is
4.66e-14 + 2.61e-14 = 7.27e-14. The total allowance charge, including the one
`eta_Q` already in the full-factorial reference, is 67 `eta_Q` = 2.65e-14.

In fp64 the pathwise design costs 15 per cent more bound than the full
factorial at 21 of 36 cells; the pooled design costs 56 per cent more at 75 of
144 cells. In fp32 the alias allowance dominates: the reduced-design bound is
about ten times the reference on the pathwise cube and about forty times on the
pooled cube. The reduced designs are fp64 designs; the record does not claim an
fp32 saving.

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
1.39e-16, well inside the allowance.

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
python pooled_d75_record_check.py
```

The first two write their JSON beside the script and print the certificate.
The third reads both cell archives and writes `pooled_d75_record_check.json`;
it takes under a minute on the recorded machine. The supporting reductions are:

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

The record-check, allowance and second-instance scripts locate the harness
repository two levels above this directory. If the record is extracted
elsewhere, set `BOYS_HARNESS_REPO` to the repository root and those scripts
then find the deposited inputs.

## The premise that matters

The packed-sum cone is `{partition, pack}` once the budget is fixed. The
packed single-precision sum's local error depends on how the sum is split
(`partition`) and on the packing decision (`pack`); after the budget is fixed
it does not depend on `route`, `scheme`, `form` or `exp`. This is the one
modelling premise that sets the 75-run count.

Widening the cone changes the design. If the packed-sum cone admitted, say,
`form` or `exp`, then `P` would contain more than 75 components, the permitted
block would need more than 75 independent rows, and the 75-run design would no
longer be estimable. The 75 is `|P|`, not a search artefact. The same reading
gives the pathwise cones `{partition, form, exp}` and `{partition, pack}` and
the 21-run count. Any claim that uses these run counts inherits this premise.

## Files

| file | role |
| --- | --- |
| `instance.py` | four-axis pathwise instance: axes, cones, basis, `G_alpha`, `P`, `T`, `Q` |
| `field.py` | exact arithmetic in Q(sqrt 2, sqrt 3) |
| `exact_verify.py` | exact design evaluator (rank, `Omega_PP`, `B`, `kappa`), random screening; writes `exact_designs.json` |
| `certify_d21.py` | file-backed D21 certificate; writes `certify_d21.json` |
| `exact_proofs.py` | condition-A rank 36 and the 18-cell null relations |
| `enumerate_designs.py` | exhaustive 0/1 MILP search for `Omega_{T,R0} = 0` |
| `pooled_allowance.py` | six-axis pooled cones and the record's `eta_Q`, target and order-6 allowances; writes `pooled_allowance.json` |
| `pooled_d75.py` | D75 construction and candidate evaluation; writes `pooled_d75.json` |
| `pooled_d75_exact.py` | modular rank proof and alias spectrum; writes `pooled_d75_exact.json` |
| `certify_d75.py` | file-backed D75 certificate; writes `certify_d75.json` |
| `pooled_d75_record_check.py` | applies D75 to every record point in both arms and resolving lanes; writes `pooled_d75_record_check.json` |
| `pathwise_summary.py` | per-path reductions of the four-axis cube; writes `pathwise_summary.json` |
| `second_instance_star.py` | four-axis second-instance star design; writes `second_instance_star.json` |
| `explore.py` | polytope exploration for designs orthogonal to the mean and main effects |
| `general_search.py` | single-removal, pair and simulated-annealing design search |

Each `.log` file is the output of the corresponding script at the recorded
environment. `MANIFEST.sha256` covers every file in this directory except
itself.

## Limitations

- The design algebra is exact; the allowances `eta_Q` and the frontier bounds
  are measured on one deposited record, one kernel revision and one platform.
- The fp16 and bf16 lanes are present in the archives but are not used.
- The separate arm differs slightly from the fused arm in the fp64 target
  maximum and in the point counts; both are reported.
- The second-instance star is evaluated on the recorded full-factorial
  coefficients and errors, not on a new measurement run.
- The "points above alias-only bound" column is expected to be large because
  the bound in that check is the alias term alone. It is not a failure of the
  certified bound, which also carries `delta_alpha`.
