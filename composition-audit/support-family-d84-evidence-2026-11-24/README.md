# D84 prospective evidence

This record re-runs the three prospective checks of the claim-gate
article's Section 4.6 on the corrected 84-run pooled design D84,
computes the D84 frontier ratio and its ratio distribution from the
measured record, and settles the joint contamination constant in exact
arithmetic. It is the D84 counterpart of the prospective package that
accompanies the article; the D75 numbers in that package are not reused.

The D84 certificate itself is the record at tag
`support-family-d84-2026-11-24`, in
`composition-audit/support-family-d84-2026-11-24/`. This directory
deposits the evidence that the certificate's constants produce on the
measured record.

## The design

The pooled cube has axes route(2), scheme(2), partition(3), pack(2),
form(3) and exp(2). The corrected packed-sum cone is
`{route, scheme, partition, pack}`. The permitted set is the 72
pack-free channels of the upward cone plus the 24 packed-cone channels
less their 12-channel intersection: **84 permitted, 60 complement and 75
targets** (every permitted channel of order 2 to 5). The joint
contamination constant is **kappa_joint = 4**.

D84 is the 72-cell five-axis factorial at `pack = arguments` plus the 12
cells at `pack = orders` with route, scheme in {0, 1}, partition in
{0, 1, 2} and form = exp = 0. The cell list is in the certificate.

`d84_fixture.py` builds the design from the harness basis and verifies it
against `certify_d84.json`: cell list, counts, rank 84, `G_alpha`
spectrum and sum (75), per-target alias constant 1, amplification
14.024284, measurement factors 9.0 and 357.222222, and the joint constant
values 1, 2, 2.5, 3 and 4. All checks pass.

## The measured record

The checks read the composition audit's two deposited cell archives:

```
4af5694810099d07ac9eee268f20073f08d63242731d08a06ebc3d1328c0638d  cells-msvc-fused.csv.xz
f3493c4ca6aceaf642bcb8e5ef4b4d330c1214f64e467e284c7526ed1094b5bb  cells-msvc-separate.csv.xz
```

The archives are the record at Boys tag `v5.0.0`, revision
`482fd3388c276091834a6c12a642d9de2481c9f2`. Each arm has 858 points
(33 values of `n` times 26 values of `x`) and 144 cells per point. The
scaled figure is `max_cells figure(cell,t) / F0(x)`.

## A. Split-sample transport

Each of 100 random permutations of the 858 points is cut into two
429-point halves. The scalar allowance is the maximum over the estimating
half of the row's pointwise allowance, and the certified inequality

```
max_a |N(a)| <= S(t) + kappa_joint * eta_est + delta_meas,
S(t) = sum_T G_alpha |chat_alpha(t)|
```

is checked at every held-out point. The design coefficients are the
full-record ones, so the check transports the allowance and does not
refit the certificate. The cell-tight variant replaces `S(t)` by
`|sum_T chat_alpha Phi_alpha(a)|` at the tested cell; it is a shortcut,
not a certificate.

| arm | lane | design | coverage | worst ratio | median | p95 | cell-tight coverage | cell-tight worst |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| fused | fp64 | path (D21) | 1 | 0.518949 | 6.33e-5 | 0.081603 | 0.999860 | 1.309 |
| fused | fp64 | pooled (D84) | 1 | 0.988824 | 0.055639 | 0.623794 | 0.999977 | 1.003 |
| fused | fp32 | path (D21) | 1 | 0.691234 | 0.001021 | 0.059471 | 0.999720 | 3.371 |
| fused | fp32 | pooled (D84) | 1 | 0.625241 | 0.002706 | 0.100365 | 0.999814 | 2.870 |
| separate | fp64 | path (D21) | 1 | 0.510094 | 6.01e-5 | 0.071519 | 1 | 0.897 |
| separate | fp64 | pooled (D84) | 1 | 0.980739 | 0.053313 | 0.605286 | 0.999814 | 1.018 |
| separate | fp32 | path (D21) | 1 | 0.691234 | 0.000947 | 0.059810 | 0.999790 | 3.371 |
| separate | fp32 | pooled (D84) | 1 | 0.622855 | 0.002670 | 0.101820 | 0.999907 | 2.365 |

Every row has full coverage over the 42 900 held-out trials. The path
rows use the unchanged D21 design; only the pooled rows change with D84.
The cell-subsample variant (72 of 144 cells) also has full coverage on
every row, with worst ratios 0.093205, 0.224269, 0.093389 and 0.224269 on
the D84 pooled rows.

## B. Planted channel

The injected channel is component 8, the order-two channel supported on
pack and form, with channel choices (1, 1) and `G = sqrt(3/2) =
1.224744871391589`. It lies in neither cone. The null allowance is
`delta0 = kappa_joint * max_t eta_Q(t)` separately for each design and
lane. The injection is swept over `|delta|/delta0` in
`[1e-3, 1e3]` on a 400-point geometric grid, plus zero, at both signs.

Four allowance variants are measured: self-consistent exact (the
planted channel's fitted coefficient re-enters the allowance), self-
consistent upper, global frozen (the null allowance held at its record
maximum) and pointwise frozen (the row's pointwise null allowance).

**Exact alias result.** The exact-arithmetic computation below shows that
all 75 target aliases of this channel are **exactly zero** under D84: the
pooled design is exactly blind to it in its targets, and the leakage is
visible only in the bound. The floating-point value in the tables
(4.44e-16) is rounding noise. The D21 path design, by contrast, maps the
channel into a target at coefficient 1.0000000000000002.

| arm | lane | delta0 | frozen threshold (global) | ratio at delta0, global frozen | ratio at delta0, self-consistent | pointwise frozen threshold | ratio at 1000 delta0, self-consistent |
| --- | --- | --- | --- | --- | --- | --- | --- |
| fused | fp64 | 1.108373e-15 | 0.816497 | 1.224745 | 0.781147 | 2.63e-12 | 0.252636 |
| fused | fp32 | 6.842362e-7 | 0.816497 | 1.224745 | 0.267430 | 4.26e-14 | 0.250018 |
| separate | fp64 | 1.108373e-15 | 0.816497 | 1.224745 | 0.781148 | 2.63e-12 | 0.252627 |
| separate | fp32 | 6.842362e-7 | 0.816497 | 1.224745 | 0.270213 | 4.26e-14 | 0.250021 |

The frozen threshold is the bisection estimate of the crossing between
the last holding and first failing grid points; the global frozen
allowance fails at `delta/delta0 = 0.816497` on all four pooled rows and
at both signs. The self-consistent allowance holds to `10^3 delta0` on
every row. The ratio at `delta0` under the exact re-measured allowance is
at most 0.782 on fp64 and 0.271 on fp32. The path rows are unchanged from
the D21 record: the global frozen allowance does not fail up to
`10^3 delta0`, and the self-consistent ratio at `10^3 delta0` is 0.200.

## C. Ratio distribution

The pointwise family/full ratio

```
(tilde_figure(t) + M_hat(t) + kappa_joint eta_Q(t))
---------------------------------------------------
          (tilde_figure(t) + N_full(t))
```

is evaluated at each of the 858 points with `kappa_joint = 4`.

| arm | lane | design | median | p95 | maximum | within 1% | within 5% |
| --- | --- | --- | --- | --- | --- | --- | --- |
| fused | fp64 | path (D21) | 1.000000 | 1.000079 | 1.011051 | 99.88% | 100% |
| fused | fp64 | pooled (D84) | 1.000000 | 1.000080 | 1.009916 | 100% | 100% |
| fused | fp32 | path (D21) | 1.000000 | 1.014543 | 1.619582 | 94.29% | 98.37% |
| fused | fp32 | pooled (D84) | 1.000000 | 1.013176 | 1.635887 | 94.29% | 98.37% |
| separate | fp64 | path (D21) | 1.000000 | 1.000085 | 1.011051 | 99.88% | 100% |
| separate | fp64 | pooled (D84) | 1.000000 | 1.000086 | 1.009916 | 100% | 100% |
| separate | fp32 | path (D21) | 1.000000 | 1.014543 | 1.657322 | 94.29% | 98.37% |
| separate | fp32 | pooled (D84) | 1.000000 | 1.013130 | 1.705403 | 94.29% | 98.37% |

The medians are 1 plus at most 7.3e-14. All eight maxima are at
`x = 11.5`.

## The frontier and its percentiles

The same pointwise ratio, maximised over the 858 points, is the
certificate's frontier. `frontier_d84.py` recomputes it from the record
and confirms the certificate block to six decimal places:

| arm | lane | frontier | certificate block | median | p95 | within 1% | within 5% |
| --- | --- | --- | --- | --- | --- | --- | --- |
| fused | fp64 | 1.009916 | 1.009916 | 1.000000 | 1.000080 | 100% | 100% |
| fused | fp32 | 1.635887 | 1.635887 | 1.000000 | 1.013176 | 94.29% | 98.37% |
| separate | fp64 | 1.009916 | 1.009916 | 1.000000 | 1.000086 | 100% | 100% |
| separate | fp32 | 1.705403 | 1.705403 | 1.000000 | 1.013130 | 94.29% | 98.37% |

The per-term maxima are:

| arm | lane | N_full | eta_Q max | charge 4 eta_Q | M_hat max | family max | full max |
| --- | --- | --- | --- | --- | --- | --- | --- |
| fused | fp64 | 2.097505e-14 | 2.770933e-16 | 1.108373e-15 | 2.097505e-14 | 1.209751e-13 | 1.209751e-13 |
| fused | fp32 | 2.221286e-7 | 1.710590e-7 | 6.842362e-7 | 2.660919e-7 | 1.906960e-6 | 1.595769e-6 |
| separate | fp64 | 2.090273e-14 | 2.770933e-16 | 1.108373e-15 | 2.090273e-14 | 1.209027e-13 | 1.209027e-13 |
| separate | fp32 | 2.145560e-7 | 1.710590e-7 | 6.842362e-7 | 2.660919e-7 | 1.906960e-6 | 1.595769e-6 |

The measurement term moves the ratio by at most 2.5e-14 absolute and is
invisible at the printed precision. The ratio of the two maxima (family
maximum over full maximum, a different quantity from the pointwise
maximum) is 1.0000000000000004 on fp64 and 1.1950096979 on fp32.

`d84_frontier_recompute.json` in the certificate directory reports a
one-allowance ratio (the allowance charged once, not with the joint
constant of 4). The corresponding values here are 1.001678,
1.195658, 1.001912 and 1.246466; they reproduce that file and are **not**
the certified frontier. The certified frontier is the kappa_joint = 4
value in the table above.

## Exact joint contamination constant

`kappa_joint_d84.py` repeats the D21 exact computation for D84 over
`Q(sqrt 2, sqrt 3)`. It builds the exact 144-cell basis, forms the
84 x 84 permitted block, inverts it exactly (no floating-point inverse),
forms the exact alias matrix, and computes the residual
`Phi_beta - sum_T B_{alpha,beta} Phi_alpha` for every beta in the
60-channel complement at every one of the 144 cells. The supremum is
decided exactly: the sign of `k^2 G_beta^2 - r_beta^2` is obtained from
rational interval bounds on `sqrt2` and `sqrt3`, refined until the
interval excludes zero, so no floating-point comparison enters the
verdict.

The result:

- **kappa_joint = 4 exactly.**
- All 60 complement channels have an exact rational supremum; the exact
  per-channel values are 1, 2, 5/2, 3 and 4, matching the floating-point
  spectrum but without rounding.
- The exact alias matrix agrees with the floating-point one to
  1.33e-15.
- The planted channel's 75 target aliases are exactly zero.

## Relationship to the D75 prospective package

The prospective package that accompanies the article computes its pooled
rows with the earlier 75-run design (D75, `kappa_joint = 3`, delta0
1.186e-15 fp64 and 7.127e-7 / 6.531e-7 fp32). The D84 values differ:

| quantity | D75 | D84 |
| --- | --- | --- |
| pooled transport worst ratio, fused fp64 / fp32 | 0.975262 / 0.795951 | 0.988824 / 0.625241 |
| pooled transport worst ratio, separate fp64 / fp32 | 0.971149 / 0.655768 | 0.980739 / 0.622855 |
| pooled delta0, fp64 (fused / separate) | 1.186232e-15 / 1.186236e-15 | 1.108373e-15 / 1.108373e-15 |
| pooled delta0, fp32 (fused / separate) | 7.127460e-7 / 6.531338e-7 | 6.842362e-7 / 6.842362e-7 |
| global frozen threshold, pooled | 0.796162 (fp64); 0.712398 fused / 0.711854 separate (fp32) | 0.816497 / 0.816497 |
| self-consistent ratio at 1000 delta0, pooled | 0.336112 / 0.333352 (fused) | 0.252636 / 0.250018 |
| largest target alias of the planted channel | ~7.4e-16 (rounding floor) | exactly 0 |
| pooled distribution maximum, fused fp64 / fp32 | 1.010661 / 1.660344 | 1.009916 / 1.635887 |
| pooled distribution maximum, separate fp64 / fp32 | 1.011082 / 1.641661 | 1.009916 / 1.705403 |

The D84 distribution maxima are the values in the article's current
distribution table and the certificate's frontier block; the D75
prospective values are different numbers. The D84 transport and planted
rows above are the ones that belong to D84.

## Self-check against the D75 package

`prospective_d84.py --design d75` runs the same code with the earlier
75-run cone. Its transport and distribution CSV files reproduce the
article package's `transport.csv` and `distribution.csv` byte for byte,
and its planted row counts, delta0 values and frozen thresholds agree
with the package's `planted.csv`; the only differences are last-bit
rounding in quantities at the 1e-30 level and tie-breaking in the
reported worst point. This establishes that the D84 changes are the
design substitution, not a change of protocol.

## Reproduction

Set `BOYS_HARNESS_REPO` to the root of this repository (the scripts also
derive it from their location), then run from this directory:

```
python -I d84_fixture.py
python -I prospective_d84.py
python -I frontier_d84.py
python -I kappa_joint_d84.py
```

`d84_fixture.py` writes `d84_fixture.json`; `prospective_d84.py` writes
`prospective_d84_results.json`, `prospective_d84_transport.csv`,
`prospective_d84_planted.csv`,
`prospective_d84_planted_selfconsistent.csv` and
`prospective_d84_distribution.csv`; `frontier_d84.py` writes
`frontier_d84.json`; `kappa_joint_d84.py` writes
`kappa_joint_d84.json`. The prospective script takes about 2 minutes on
the recorded machine, the fixture and frontier scripts take seconds, and
the exact computation takes about 2.5 minutes. The exact computation is
deterministic and rewrites its JSON identically apart from the recorded
runtime.

`MANIFEST.sha256` lists the digests of the files in this directory.

## Requirements and limitations

- Python 3.10 or newer with `numpy`; `lzma` and `csv` from the standard
  library. No C++ compiler is needed.
- The design algebra is exact. The check values, the frontier and the
  distribution are measured on one deposited record, one kernel revision
  and one platform.
- The allowance is measured on the full cube that the smaller design
  would replace. The checks are retrospective: they validate the algebra
  and the estimator relation, not a prospective saving on new
  measurements.

## Files

| file | role |
| --- | --- |
| `d84_fixture.py` | builds D84 from the harness basis and verifies it against `certify_d84.json` |
| `record_figures.py` | reads the record's figure column and forms the scaled figure |
| `prospective_d84.py` | the three prospective checks; `--design d75` reproduces the earlier package |
| `frontier_d84.py` | the frontier ratio, per-term maxima and ratio percentiles |
| `kappa_joint_d84.py` | the exact joint constant over `Q(sqrt2, sqrt3)` |
| `*_results.json`, `*.csv` | the measured outputs listed under Reproduction |
| `frontier_d84.json` | frontier and percentile record |
| `kappa_joint_d84.json` | exact per-channel constants and the planted alias |
| `*.log` | the printed output of each run |
| `MANIFEST.sha256` | file digests |
