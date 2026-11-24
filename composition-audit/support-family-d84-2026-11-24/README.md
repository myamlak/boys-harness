# The corrected 84-run pooled design

This record deposits the certificate for the 84-run pooled design D84 and the
retrospective record check on the composition audit's 576-cell measurement. It
corrects the packed-sum cone. The packed single-precision sum depends on
`route` and `scheme` as well as on `partition` and `pack`, so the cone is
`{route, scheme, partition, pack}`. The permitted set is therefore 84 channels,
not the 75 of the earlier record, and the pooled design is D84.

The D84 certificate and record check in this directory are the record at tag
`support-family-d84-2026-11-24`. The earlier D75 records remain at tag
`support-family-2026-11-22-fix3`, unchanged.

## The permitted set

The pooled cube has axes and levels `route` (2), `scheme` (2), `partition`
(3), `pack` (2), `form` (3) and `exp` (2). The cell index is

```
route*72 + scheme*36 + partition*12 + pack*6 + form*2 + exp.
```

The two path cones are:

- upward step `{route, scheme, partition, form, exp}`, with
  `2*2*3*3*2 = 72` channels;
- packed sum `{route, scheme, partition, pack}`, with
  `2*2*3*2 = 24` channels;
- their intersection `{route, scheme, partition}`, with `2*2*3 = 12` channels,
  which is counted once.

So the permitted set is `72 + 24 - 12 = 84` channels: the mean, eight non-mean
main effects and 75 targets (every permitted component of order 2 to 5). The
complement has 60 structural zeros. The support family has
`2^5 + 2^4 - 2^3 = 40` distinct supports.

## The design

D84 is all 72 cells with `pack = arguments` (pack level 0), plus the 12 cells
with `pack = orders` (pack level 1) at `route` 0/1, `scheme` 0/1,
`partition` 0/1/2 and `form = exp = 0`. The twelve added cells are

```
6, 18, 30, 42, 54, 66, 78, 90, 102, 114, 126, 138.
```

The full cell list is in `certify_d84.json`.

## The certificate

- The 84 x 84 permitted block has rank 84 of 84 numerically and rank 84 of 84
  modulo each of 1009, 2017 and 3001.
- The Gram matrix has condition number 12, so the model-matrix condition number
  is `sqrt(12) = 2 sqrt 3 = 3.464102`.
- The per-target alias constant is 1: every one of the 75 targets has
  `G_alpha kappa_alpha = 1` within 1.3e-15, so the per-target sum is 75. The
  recorded value is 1.000000000000001.
- The theorem's joint contamination constant is `kappa_joint = 4`. The 60
  structural-zero ratios take the values 1, 2, 2.5, 3 and 4, and the maximum is
  4. **The per-target alias constant (1) and the joint contamination constant
  (4) are different quantities; the per-target constant is not the joint
  constant.**
- The threshold amplification, maximised over the 84 permitted components, is
  `max_alpha sum_gamma |(Omega_PP^-1)_{alpha gamma}| G_gamma = 14.024284`.
- The measurement factors are 9.0 at the maximum target and 357.222222 as the
  `G_alpha`-weighted sum over the targets.
- The 12 x 12 Kronecker block determinant is
  `2^24 * 3^6 = 12,230,590,464`.

### Exact rank proof

The 72 pack-mean columns restricted to the 72 `pack = arguments` rows are the
five-axis orthonormal basis on `route x scheme x partition x form x exp`, so
they are invertible. Subtracting the matching pack-mean column from each of the
twelve pack-contrast columns leaves a matrix supported on the twelve
`pack = orders` rows. Its 12 x 12 block is the Kronecker product of the route
(2 x 2), scheme (2 x 2) and partition (3 x 3) orthonormal bases evaluated at the
selected levels, with the factor -2 from the pack subtraction on each of the
twelve rows. Each factor is invertible, so the block is invertible and the
permitted block has rank 84. The block determinant is

```
(-2)^12 (det R)^6 (det S)^6 (det P)^4 = 2^24 * 3^6 = 12,230,590,464,
```

because `det R = det S = -2` and `det P = 3 sqrt 3`.

## The frontier

The frontier is the maximum over the record's 858 points of the pointwise
family/full ratio, with `kappa_joint = 4`. At each point the family bound is
the pointwise figure plus the measured target reconstruction `M_hat` plus
`4 eta_Q`, and the full bound is the same figure plus the full cube's measured
non-additive part `N_full`; the maxima of the two bounds need not occur at the
same point.

| lane | fused | separate |
| --- | --- | --- |
| fp64 | 1.009916 | 1.009916 |
| fp32 | 1.635887 | 1.705403 |

The values are the `frontier` block of `certify_d84.json`.
`d84_frontier_recompute.json` holds the per-arm term breakdown and the
record-check columns; its `ratio_pointwise_new` fields charge the allowance
once, not with the joint constant of 4.

## The record check

`pooled_d84_record_check.py` applies D84 to all 858 points of each deposited
arm in the fp64 and fp32 resolving lanes. At every point it computes the
full-cube coefficient vector, the design estimate `chat`, the alias identity
residual `chat - c_P - B c_Q`, the largest target estimate, and whether that
estimate exceeds the alias-only bound `kappa_alpha eta_Q`.

| arm | lane | max abs(chat - c) | alias identity residual | max target estimate | points above alias-only bound |
| --- | --- | --- | --- | --- | --- |
| fused | fp64 | 4.3541e-17 | 1.447e-29 | 1.329474e-14 | 857 / 858 |
| fused | fp32 | 5.375867e-8 | 8.603e-23 | 6.584066e-8 | 820 / 858 |
| separate | fp64 | 4.3541e-17 | 1.447e-29 | 1.331010e-14 | 857 / 858 |
| separate | fp32 | 5.375867e-8 | 8.285e-23 | 7.133651e-8 | 821 / 858 |

The alias identity is exact algebra, so its residual is at the
double-precision floor; this is the check that the reduction is wired
correctly. The alias-only bound is the alias term alone, not the full certified
bound, so the large point counts are expected: the target estimate exceeds the
alias term wherever the permitted signal or the measurement term is larger.
The same values are in `d84_frontier_recompute.json`.

## Commands

Run from this directory:

```
python certify_d84.py
python pooled_d84_record_check.py
```

`certify_d84.py` prints the core certificate and writes `certify_d84.json`.
The script replaces that file, so run it in a copy if the deposit must stay
unchanged; the deposited JSON also carries the independently verified
constants and the frontier block, which the script does not write.

`pooled_d84_record_check.py` reads the two deposited cell archives and writes
`pooled_d84_record_check.json`. It takes about 35 seconds on the recorded
machine.

`verify_cert.py` is an independent check built from the claim-gate worked
check's primitives. It prints the modular ranks, the permitted and complement
counts, the per-target constant, the joint constant and its per-beta values,
the Gram and model conditions, the amplification, the measurement factors and
the block determinant. It writes no file. The checker module
`design_rule_check.py` is external to this deposit; set
`DESIGN_RULE_CHECK_DIR` to the directory that holds it:

```
DESIGN_RULE_CHECK_DIR=/path/to/worked-check python verify_cert.py
```

In PowerShell:

```
$env:DESIGN_RULE_CHECK_DIR = 'C:\path\to\worked-check'; python verify_cert.py
```

## Requirements and provenance

Python 3.10 or newer with `numpy` is enough; the record check uses `lzma` from
the standard library. No C++ compiler is needed. The record check and the
frontier read the composition audit's own cell archives, with SHA-256

```
4af5694810099d07ac9eee268f20073f08d63242731d08a06ebc3d1328c0638d  cells-msvc-fused.csv.xz
f3493c4ca6aceaf642bcb8e5ef4b4d330c1214f64e467e284c7526ed1094b5bb  cells-msvc-separate.csv.xz
```

The archives are the measured record at Boys tag `v5.0.0`, revision
`482fd3388c276091834a6c12a642d9de2481c9f2`. The D84 design algebra is exact
over `Q(sqrt 2, sqrt 3)`. The recorded environment is Windows 11 Pro
10.0.26200, Python 3.14.7 and numpy 2.5.2.

## Files

| file | role |
| --- | --- |
| `certify_d84.py` | D84 construction and core certificate; writes `certify_d84.json` |
| `certify_d84.json` | cells, rank, modular ranks, alias spectrum, amplification, joint constant, measurement factors, block determinant and frontier |
| `verify_cert.py` | independent verification with the external claim-gate worked check; prints its checks |
| `pooled_d84_record_check.py` | applies D84 to all 858 points of both arms and resolving lanes; writes `pooled_d84_record_check.json` |
| `d84_frontier_recompute.json` | per-arm frontier terms and record-check values at the corrected cone |
| `pooled_allowance.py` | six-axis basis, supports, `G_alpha` and the full-cube decomposition; the D84 scripts override its packed cone |
| `pooled_d75.py` | cell index and five-axis cell constructors used by the D84 construction |
| `instance.py` | four-axis pathwise basis used by the shared allowance module |
| `README.md` | this file |

## Limitations

- The design algebra is exact; the frontier and record-check values are
  measured on one deposited record, one kernel revision and one platform.
- `eta_Q` is measured on the full cube that the smaller design would replace.
  The record check is retrospective: it validates the algebra and the estimator
  relation, not a prospective saving on new measurements.
- `certify_d84.py` regenerates the core certificate fields but not the
  independently verified constants or the frontier block; those are recorded
  in the deposited `certify_d84.json`.
- `verify_cert.py` needs the external `design_rule_check.py` module and is not
  self-contained.
