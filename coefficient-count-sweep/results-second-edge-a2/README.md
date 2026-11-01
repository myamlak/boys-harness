# Independent rational counts at the second edge

These are the paper's rational counts at the second Boys left edge
`a2 = 2 x0 = 23.799696304216968`, measured with the same independent solver
that produced `results-corrected`: `baryrat` 2.1.2, the BRASIL rational
minimax algorithm, not the kernel's Remez or Lawson implementations. Every
Group A width is doubled, so `W/a` and the cross-ratio are matched exactly.
`kappa(a2, 2W) = kappa(x0, W)` to a maximum relative difference of `8.3e-59`
at 60 decimal digits, and the interval midpoint `b = a + W/2` doubles exactly
(`b_a2 = 2 b_x0`).

The `x0` side of every comparison is `../results-corrected/corrected-counts.json`
(sha256 `9bb8ac0633b10e8a7922684ca8a419135103f2c7825b4d06b9c71603399b7195`).
Both edges therefore come from this one solver. All 45 cells completed.

## Method

Target function `F_0(x) = sqrt(pi)/(2 sqrt(x)) erf(sqrt(x))` on `[a2, a2+W]`,
where `W` is twice the `x0` Group A width. The rational form is
`p(t)/(1 + q1 t + ... + qk t^k)` in `t = 2(x-a2)/W - 1`; the stored count is
`(m+1)+k` with `k <= 12`.

The paper's acceptance rule was applied unchanged. The barycentric fit is
converted exactly to monomials at 70 decimal digits and rounded to binary64;
it is evaluated by Horner in binary64 by both the fused (`math.fma`) and the
separate route, and the error is the worse of the two. The target's kind is
relative (`max |fhat-f|/f`: `1e-13` or `1e-14`) or absolute
(`max |fhat-f|`: `5e-14`). The denominator is guarded at
`min |den| >= 1e-3 max |den|`. The 2001-node uniform-`t` grid certifies, and
the 4001-node grid checks. The reference is `F_0` at 70 decimal digits.

The default BRASIL settings (`tol=1e-4`, `maxiter=1000`, `npi=-30`,
`init_steps=100`, `max_step_size=0.1`) are scanned upward in stored count from
3 to the `x0` count plus 3. If the first meeting level is `s`, level `s-1` is
re-solved with the tightened settings (`tol=1e-8`, `maxiter=5000`, `npi=-40`,
`init_steps=300`, `max_step_size=0.05`); if it meets, the count drops and the
next level below is retried. The tightened retry ran in every cell and never
met, so no count was lowered. The witness is the admissible fit that also
passes the 4001 grid and its guard, preferring a converged fit and then the
smallest 4001 error. A count is an existence upper bound under this solver,
not a proof of minimality.

The absolute target is carried through `L = log(F_0(a)/eps)` at each edge, as
the paper's second-edge record does. The values are `29.933606208922594` for
`rel:1e-13` at both edges, `29.267707262998931` at `x0` against
`28.921134741707567` at `a2` for `abs:5e-14` (shift `-0.34657252129136455`),
and `32.23619130191664` for `rel:1e-14` at both edges.

## Counts

Stored counts, `x0` from the corrected record and `a2` from `cells/`:

| width | W_x0 | W_a2 | rel:1e-13 | abs:5e-14 | rel:1e-14 |
| --- | --- | --- | --- | --- | --- |
| m7 | 0.0078125 | 0.015625 | 4 -> 4 | 4 -> 4 | 4 -> 4 |
| m6 | 0.015625 | 0.03125 | 4 -> 4 | 4 -> 4 | 4 -> 4 |
| m5 | 0.03125 | 0.0625 | 4 -> 4 | 4 -> 4 | 4 -> 4 |
| m4 | 0.0625 | 0.125 | 4 -> 4 | 4 -> 4 | 5 -> 5 |
| m3 | 0.125 | 0.25 | 5 -> 5 | 5 -> 5 | 5 -> 5 |
| m2 | 0.25 | 0.5 | 5 -> 5 | 5 -> 5 | 5 -> 6 |
| m1 | 0.5 | 1.0 | 6 -> 6 | 6 -> 6 | 6 -> 6 |
| p0 | 1.0 | 2.0 | 6 -> 6 | 6 -> 6 | 6 -> 7 |
| p1 | 2.0 | 4.0 | 7 -> 7 | 7 -> 7 | 9 -> 8 |
| p2 | 4.0 | 8.0 | 9 -> 8 | 9 -> 8 | 10 -> 9 |
| p3 | 8.0 | 16.0 | 11 -> 10 | 10 -> 9 | 11 -> 10 |
| p4 | 16.0 | 32.0 | 12 -> 11 | 12 -> 11 | 14 -> 12 |
| p5 | 32.0 | 64.0 | 15 -> 13 | 14 -> 12 | 18 -> 14 |
| ship | 19.100151847891517 | 38.200303695783035 | 13 -> 12 | 12 -> 11 | 15 -> 12 |
| paper | 17.089489586712254 | 34.17897917342451 | 12 -> 11 | 12 -> 11 | 14 -> 12 |

Count shifts over the fifteen widths, `N(a2) - N(x0)`:

| target | differing cells | sum | largest absolute shift | parity flips | pair changes |
| --- | --- | --- | --- | --- | --- |
| rel:1e-13 | 6 of 15 | -7 | 2 | 5 | 12 |
| abs:5e-14 | 6 of 15 | -7 | 2 | 5 | 12 |
| rel:1e-14 | 9 of 15 | -12 | 4 | 6 | 13 |

For `rel:1e-13` and `abs:5e-14` no shift is positive; the wide widths are
lower at `a2` by up to 2. For `rel:1e-14` the shifts are `+1` at `m2` and
`p0` and negative from `p1` upward, reaching `-4` at `p5`.

## Slopes

The paper's estimator: least squares of `log2 N` on `log2 W` over the
thirteen dyadic widths `2^-7..2^5`, floor rows included, standard error
`sqrt(SSR/(n-2)/SXX)`. Both edges are fitted with the same code and the same
`n = 13`. The paired test carries the covariance of the two edges' residuals
at the same widths:

`se_paired = sqrt(se_x0^2 + se_a2^2 - 2 cov/SXX)`.

| target | x0 slope | a2 slope | difference | unpaired s.e. | unpaired z | paired s.e. | paired z |
| --- | --- | --- | --- | --- | --- | --- | --- |
| rel:1e-13 | 0.16455 +/- 0.01384 | 0.14847 +/- 0.01117 | -0.01608 | 0.01779 | -0.90 | 0.00372 | -4.32 |
| abs:5e-14 | 0.15824 +/- 0.01257 | 0.14132 +/- 0.00991 | -0.01692 | 0.01601 | -1.06 | 0.00391 | -4.33 |
| rel:1e-14 | 0.18051 +/- 0.01572 | 0.15483 +/- 0.00791 | -0.02568 | 0.01760 | -1.46 | 0.01063 | -2.41 |

The unpaired test sees no slope difference at any target. The paired test
sees the one-sided shift pattern at 4.3, 4.3 and 2.4 sigma. The effect is
0.016 to 0.026 in slope, 10 to 14 per cent of the `x0` slope. The endpoint
secants move the same way: `0.1589 -> 0.1417`, `0.1506 -> 0.1321`,
`0.1808 -> 0.1506`.

With the `x0` constant carried across the matched cross-ratio, the model's
real-valued prediction at `a2` is `N(x0)` for the relative targets and
`N(x0) - 0.346573/mu` for the absolute target. The observed residual
`N(a2) - N_model` ranges `-2.00` to `0.00` for `rel:1e-13`, `-1.86` to
`+0.075` for `abs:5e-14`, and `-4.00` to `+1.00` for `rel:1e-14`. The
count-implied `C_r` values are in `analysis.json`.

## Cross-check against the kernel's partial a2 record

The paper's second-edge record has five completed rational cells at `a2`, all
`rel:1e-13` at the widths below `2^-2`. This independent solver reproduces
all five exactly in stored count and in `[m/k]` pair: `4 [1/2]`, `4 [1/2]`,
`4 [1/2]`, `4 [1/2]`, `5 [2/2]` for `2^-7` through `2^-3`. The other forty
cells are new. No count or pair in the existing record is contradicted.

## Verification

`verify_witnesses.py` re-evaluates every stored witness from its coefficients
with a separate implementation, rebuilding each grid from the record's own
`a` and `W`, at 70 decimal digits, on both grids and by both Horner routes.
All 45 cells pass, and every recomputed error and guard ratio agrees with the
record to within four binary64 ulps:

```
python verify_witnesses.py
45 cells, 0 failures, 35.0s          exit 0
```

The adapted driver was also run once at `x0` on `A_abs514_p1` and compared
with the verification harness record `_q330-verify/grid/A_abs514_p1.final.json`
(the full path is in `checks/x0_fidelity_comparison.json`): 29
common fields, 0 differing, all witness fields identical, including the
17-digit `P` and `Q`. The three fields the adapted driver adds (`edge`, `a`,
`x0_count_used_for_bound`) have no counterpart in the old record. The
comparison is in `checks/x0_fidelity_comparison.json`.

## Commands

Working directory is this directory. The probe and the fidelity check:

```
python second_edge_a2_solver.py --edge x0 --only A_abs514_p1 --out-dir checks 2>&1 | tee logs/probe_x0_abs514_p1.log
    driver exit 0, 1 cell, wall 31.4 s (checks/x0_check_summary.json)

python second_edge_a2_solver.py --edge a2 --only A2_abs514_p0 2>&1 | tee logs/probe_a2_abs514_p0.log
    driver exit 0, 1 cell, wall 14.5 s, 10 solves, 8.79 s in BRASIL
```

The production run was split into five chunks, one heavy job at a time, six
workers per chunk, each chunk logged to `logs/`. Each chunk wrote
`a2_run_summary.json`, which was copied to `logs/chunk<n>_summary.json`. The
cell sets are the ones in the logs; the driver selects a set, so the order of
the ids does not affect the outcome:

```
python second_edge_a2_solver.py --edge a2 --workers 6 --only A2_rel13_m7,A2_abs514_m7,A2_rel14_m7,A2_rel13_m6,A2_abs514_m6,A2_rel14_m6,A2_rel13_m5,A2_abs514_m5,A2_rel14_m5,A2_rel13_m4,A2_abs514_m4,A2_rel14_m4,A2_rel13_m3,A2_abs514_m3,A2_rel14_m3,A2_rel13_m2,A2_abs514_m2,A2_rel14_m2,A2_rel13_m1,A2_abs514_m1,A2_rel14_m1,A2_rel13_p0,A2_abs514_p0,A2_rel14_p0,A2_rel13_p1,A2_abs514_p1,A2_rel14_p1 2>&1 | tee logs/run_chunk1.log
    driver exit 0, 27 cells, wall 169.8 s (logs/chunk1_summary.json)

python second_edge_a2_solver.py --edge a2 --workers 6 --only A2_rel13_p2,A2_abs514_p2,A2_rel14_p2,A2_rel13_p3,A2_abs514_p3,A2_rel14_p3 2>&1 | tee logs/run_chunk2.log
    driver exit 0, 6 cells, wall 182.9 s (logs/chunk2_summary.json)

python second_edge_a2_solver.py --edge a2 --workers 6 --only A2_rel13_p5,A2_abs514_p5,A2_rel14_p5 2>&1 | tee logs/run_chunk3a_p5.log
    driver exit 0, 3 cells, wall 330.3 s (logs/chunk3a_summary.json)

python second_edge_a2_solver.py --edge a2 --workers 6 --only A2_rel13_p4,A2_abs514_p4,A2_rel14_p4 2>&1 | tee logs/run_chunk3b_p4.log
    driver exit 0, 3 cells, wall 270.8 s (logs/chunk3b_summary.json)

python second_edge_a2_solver.py --edge a2 --workers 6 --only A2_rel13_ship,A2_abs514_ship,A2_rel14_ship,A2_rel13_paper,A2_abs514_paper,A2_rel14_paper 2>&1 | tee logs/run_chunk3c_ship_paper.log
    driver exit 0, 6 cells, wall 405.9 s (logs/chunk3c_summary.json)
```

The driver returns 0 exactly when no cell records an error, and every log
ends `ALLDONE cells=<n> errors=0`. The production wall is 1359.7 s (22 min
40 s) for 45 cells; with both probes the solver time is 1405.6 s (23 min
26 s). Cost grows with width: a `p0` cell takes 14.5 s alone, while the
widest chunk (`ship` and `paper`, 6 cells on 6 workers) takes 405.9 s, so
the slowest cell in that chunk is about 400 s. The analysis is rerun with:

```
python analyze_second_edge_a2.py     exit 0
```

It reads `cells/A2_*.final.json` and `../results-corrected/corrected-counts.json`
and writes `analysis.json`, `second-edge-a2-counts.json` and the three CSVs.
A rerun stamps a new `generated_utc` in `analysis.json` and the assembled
record, so the hashes above bind the deposited run, not a fresh rerun.

## Files

| file | sha256 |
| --- | --- |
| `second_edge_a2_solver.py` | `f181b939d3a7091ecb0bb222806fb54d07a6174a8c998e4af23c58d41b13b4b7` |
| `verify_witnesses.py` | `a15f2c7af0eda7cdb7ce01cd8373b82c47a4639ca5bd1d2a32f9b6deaea4b381` |
| `analyze_second_edge_a2.py` | `2de78c6158544e05b1628398dd6ca26f1dd1d526f622172ead219655c9d1df0f` |
| `second-edge-a2-counts.json` | `0b98c5d33f3cf068bb8eea3182c42e08cc1e831cedd222af81369a2031e7e792` |
| `second-edge-a2-counts.csv` | `f0ddede7107c6cd815f8b5b0ba0a0f557645d4041f89e78512a6060f46b6dcbc` |
| `analysis.json` | `f49f8aa67c845d427ad43a1b0ded76002638996a0fdf43ec18d2ba8c12cd6b4a` |
| `comparison.csv` | `3e6a6a41bb2960aeb7e19fc4713be48f6a5f9af957ad4e5085f549b154bb92c6` |
| `counts_by_width_target.csv` | `222e7121c60947f52dffacbfdb8d725b4713af4ec597ab322c245c740ae05c90` |
| `checks/x0_check_summary.json` | `778d1c68a87506fbd0ec25165c2fe44854c41ea3c776213ba0c1f701d334e980` |
| `checks/x0_fidelity_comparison.json` | `b739793b8b9f6e1553c114b9093cf928d758e958dadd528db6a33a30ede6c06a` |

`cells/` holds the 45 raw scan records (`A2_*.json`, every tried `[m/k]`
with its deviation) and the 45 finalized records (`A2_*.final.json`, the
accepted pair, delivered errors, guard, flags and the 17-digit witness
`P`/`Q`). The assembled record `second-edge-a2-counts.json` repeats the
witness coefficients, so it is self-contained. `comparison.csv` carries every
compared field, including `kappa`, `L`, pairs, flags, `C_r` and the model
residual.

## Limits

A count is the smallest stored count at which an exhibited admissible fit
meets, under this solver; it is an upper bound and not a minimality proof.
Nineteen of the 45 cells are flagged `margin_limited`: 18 witnesses did not
converge, and `abs:5e-14` at `m7` has the best fit one coefficient lower
within 1.25 times the target (`below_ratio` 1.004). All 45 still meet the
4001 check, pass its denominator guard, and were re-evaluated independently.
All 45 witnesses were selected by the default route; none needed the
tightened retry. The slope comparison is between this solver's `x0` counts in
`results-corrected` and its `a2` counts here; no other record was changed.
