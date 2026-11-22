# Support-family-aware design certification: corrected derivation

This document derives a generalisation of the design rule studied here. It supersedes the
earlier version: the joint bound now charges the structural-zero allowance once through a
contamination kernel, the theorem statement carries the measurement term, and the frontier
table is corrected. Every number below is reproducible from the deposited composition cube;
the scripts and certificates are named in section 10.

## 1. Setting and the permitted support family

Fix a path p of a configured routine. Let U_p be its free axes, with levels L_i, and X_p the
full factorial cube. The response is expanded in an orthonormal basis {Phi_alpha} under the
uniform measure on X_p, with Phi_alpha(a) = prod_i phi_{i,alpha_i}(a_i) and
G_alpha = sup_a |Phi_alpha(a)|. Write e = f + rho, where f = sum_o lambda_o eps_o is the
first-order part.

The first-order cone lemma says that the first-order part of c_alpha vanishes unless S(alpha) lies
inside some operation's dependency cone. Let C_p be the path cones and define the
**permitted support family**

    P_p = { S subset U_p : exists C in C_p with S subset C } = union_{C in C_p} 2^C,

with permitted index set I_p = { alpha : S(alpha) in P_p } and Q_p its complement.

For the Boys resolving fp64/fp32 path, the four free axes are partition (3 levels), pack
(2), division form (3) and region-B exponential (2), with cones

    C1 = {partition, form, exp}   (the region-B upward step),
    C2 = {partition, pack}        (the packed single-precision sum, budget fixed by the lane).

Hence P_p = 2^{G,F,E} union 2^{G,K}, with 10 supports and 21 components (1 mean, 6 mains, 14
interactions). The fourteen targets are the interaction channels of order 2 or 3; the
fifteen structural zeros are the channels whose support is contained in neither cone:
{K,F}, {K,E}, {G,K,F}, {G,K,E}, {K,F,E} and {G,K,F,E} with their channels.

Pooled over route and scheme, the upward cone widens to five axes while the packed-sum cone
stays {partition, pack}. With axes R, S, G, K, F, E at levels (2,2,3,2,3,2),

    P_pooled = 2^{R,S,G,F,E} union 2^{G,K},

with 75 components (the 72-component five-axis basis with pack at its mean row plus the
three pack-contrast components), 66 targets and 69 structural zeros.

## 2. The certificate theorem

Assume:

(i) the permitted block is estimable: Omega_{II} is invertible;
(ii) the cone table is complete, so c_alpha(f) = 0 for alpha in Q (the cone lemma);
(iii) the structural-zero allowance eta_Q bounds the weighted coefficients outside the
    permitted family: eta_Q >= sum_{beta in Q} G_beta |c_beta(rho)|;
(iv) the adjusted measurement noise nu_hat_I = Omega_{II}^{-1} nu_I satisfies
    delta_alpha >= |nu_hat_alpha| for alpha in T. For a per-cell error |nu(a)| <= eps this
    is the deterministic bound
    delta_alpha = eps sum_gamma |(Omega_{II}^{-1})_{alpha gamma}| N^{-1} sum_{a in D} |Phi_gamma(a)|;
(T-full) T is every permitted interaction channel, T = { alpha in I : |S(alpha)| >= 2 }.
    If T is a strict subset, the omitted permitted interactions must be added to T or bounded
    separately.

Put B = Omega_{II}^{-1} Omega_{IQ}, kappa_alpha = max_{beta in Q} |B_{alpha beta}| / G_beta,
and the **joint contamination constant**

    kappa_joint = max_{beta in Q} sup_{a in X_p} |Phi_beta(a) - sum_{alpha in T} B_{alpha beta} Phi_alpha(a)| / G_beta.

**Theorem.** Under (i)-(iv) and (T-full), for alpha in T,

    |c_alpha| <= |c_hat_alpha| + delta_alpha + kappa_alpha eta_Q,

and for every a in X_p the full non-additive part
N(a) = sum_{alpha in T} c_alpha Phi_alpha(a) + sum_{beta in Q} c_beta Phi_beta(a) obeys,
with r_beta(a) = Phi_beta(a) - sum_{alpha in T} B_{alpha beta} Phi_alpha(a),

    |N(a)|
        <= sum_{alpha in T} G_alpha |c_hat_alpha| + kappa_joint eta_Q
           + sum_{alpha in T} G_alpha delta_alpha.

**Proof.** Substituting the expansion of e into the raw contrast gives
g_I = Omega_{II} c_I + Omega_{IQ} c_Q + nu_I, hence
c_hat_I = c_I + B c_Q + Omega_{II}^{-1} nu_I. The componentwise bound follows from the
triangle inequality, the definition of kappa_alpha and (iii), (iv). For the joint bound,
add the structural-zero expansion to the permitted one and use the identity

    N(a) = sum_T c_hat_alpha Phi_alpha + sum_Q c_beta r_beta
           - sum_T nu_hat_alpha Phi_alpha,

then |r_beta(a)| <= kappa_joint G_beta by definition and |Phi_alpha(a)| <= G_alpha. QED.

The first display carries the measurement term the proof produces; a display without
delta_alpha is the nu = 0 case and is false under (iv) with nonzero noise. The second
display bounds the full non-additive part, permitted interactions and structural zeros
together, which is the quantity the composition rule needs.

## 3. Specialisation to the orthogonal-array theorem

If D is an orthogonal array of strength t with 2w_p <= t <= |U_p|, then Omega_{II} = I and
Omega_{IQ} = 0 for every beta with support union at most t. A surviving alias has
|S(beta)| > t - |S(alpha)| >= w_p, so beta lies in Q and is covered by eta. With B = 0,
r_beta = Phi_beta and kappa_joint = 1; kappa_alpha <= G_alpha. The theorem becomes

    |c_alpha| <= |c_hat_alpha| + delta + G_alpha eta  (alpha in T),

which recovers the orthogonal-array theorem, and the full factorial gives the exact special
case. The family-aware theorem is a generalisation, but the recovery needs a qualification:
eta_Q bounds every component outside the family, including order-two and order-three
structural zeros, whereas the orthogonal-array allowance bounds only orders above w_p. On an
orthogonal array those low-order zeros do not alias into the targets (their unions are at
most 2w_p <= t), so the orthogonal-array theorem never needs them in the allowance. The
single-allowance bound kappa_alpha eta_Q is therefore in general looser than
delta + G_wp eta on an orthogonal array, and exact recovery holds either with an
alias-support-restricted allowance or at the proof's intermediate step
(sum over surviving aliases <= G_alpha eta). The theorem is a strict generalisation of the
design condition, not of the packaged allowance.

## 4. What is assumed and what is not

Assumed: orthonormality of the basis; the first-order Taylor model with remainder rho;
completeness of the cone table; the allowance eta_Q; the adjusted measurement threshold;
estimability of the permitted block; and (T-full) for the joint conclusion. The theorem does
not assume an orthogonal array, and it allows aliasing with structural zeros, charging the
alias once through kappa_joint.

**No-replication convention.** "Within the 36-cell cube the only strength-2 array is the
full cube" is true only for a simple design, each cell used at most once. With replication,
other 36-run strength-2 arrays exist; they use 36 runs and do not lower the minimum. The
minimality statement is unaffected.

## 5. The Boys path: no orthogonal or exact-orthogonality design below 36

On levels (3,2,3,2), any strength-2 orthogonal array has every pair marginal balanced, so N
is divisible by 9, 6 and 4, hence by lcm(9,6,4) = 36; the only simple such design inside the
cube is the full cube. Strength 3 and above implies strength 2, so no proper orthogonal
array of any strength exists.

The exact-orthogonality route is also empty below 36. Requiring target estimates to be
unaliased with the mean and mains gives 189 alias equations (98 target-to-mean/main plus 91
target-to-target) plus one normalisation row, a 190 by 36 system of exact rank 36 over
Q(sqrt2, sqrt3); its unique normalised solution is the uniform full cube. Requiring only
target-to-mean/main orthogonality leaves a three-dimensional weight polytope whose eight
vertices generate the feasible cone; an exhaustive 0/1 enumeration returns exactly three
designs: the full cube and two 18-cell designs. Both 18-cell designs have target rank 12,
with the exact null relations

    {partition, exp}(linear channel) + {partition, pack} = 0,
    {partition, exp}(quadratic channel) + {partition, pack} = 0,

up to sign. They make pack and exp complementary on their cells, so {G,E} and {G,K} become
the same function. The blocking feature is the target alias between {partition, exp} and
{partition, pack}, the two cones meeting only in partition.

## 6. D21: a certified 21-run pathwise design

The permitted family has 21 components, so at least 21 runs are needed. The structured design

    D21 = the full 3x3 partition x form grid at (pack, exp) = (0,0);
          the full 3x3 grid at (1,1);
          the three cells at (0,1) with partition 0,1,2 and form 0

has cells

    0, 1, 2, 4, 7, 9, 11, 12, 13, 14, 16, 19, 21, 23, 24, 25, 26, 28, 31, 33, 35.

Its 21x21 permitted block has exact rank 21 over Q(sqrt2, sqrt3). The Gram condition is
12, so the model-matrix condition is sqrt(12) = 2 sqrt3 = 3.464102. The target-restricted
kappa_max is sqrt(2/3) = 0.816496580927726 (over all 21 permitted rows it is 1.0), and every
one of the fourteen G_alpha kappa_alpha equals 1, so the per-target sum is 14. The joint
contamination constant is **kappa_joint = 4 exactly**, with per-beta ratios 1,1,2,3,4
({K,E} and {K,F,E}: 1; {K,F}: 2; {G,K,E}: 3; {G,K,F}: 4). The measurement factors are
max over targets 8.485281 and target sum 74.111111; with eps about 1e-29 (fp64) or 1e-22
(fp32), the measurement terms are 7e-28/2.5e-27 and 7e-21/2.5e-20, negligible against the
bounds. The certificate is reproduced by certify_d21.py; the constant by kappa_joint.py.

A different valid 21-cell design has the same target kappa_max but per-target sum 11.25 and
a worse amplification of 20.830; it is not the design cited.

## 7. D75: a certified 75-run pooled design

The pooled permitted family has 75 components: the 72 pack-mean components require a design
projecting onto all 72 five-axis cells, and the three pack-contrast components need three
further cells. The design

    D75 = all 72 cells with pack = arguments (the full five-axis factorial),
          plus the three cells with pack = orders at partition 0,1,2 and
          route = scheme = form = exp = 0

attains N = 75. Rank 75 is verified numerically, modulo 1009, 2017 and 3001, and by a Schur
block: subtracting the matching pack-mean column from each pack-contrast column leaves a
3x3 block with determinant 24 sqrt3 = 41.569219 on the three pack=orders rows. Minimality
excludes 73 and 74 runs: the 72 pack-mean columns are functions of the five-axis projection,
so 72 distinct projections are necessary; with N = 72 each projection appears once, the
pack-contrast columns lie in the pack-mean span; at N = 73 or 74 the pack-contrast residuals
after regression are supported on one or two duplicated pairs and have rank at most 1 or 2,
less than 3; only at N = 75 with three duplicated pairs at distinct partition levels is the
residual 3x3 invertible. The target-restricted kappa_max is 1.0, all sixty-six
G_alpha kappa_alpha equal 1 within 1.3e-15, and the per-target sum is 66. The Gram condition
is 48, so the model-matrix condition is sqrt(48) = 4 sqrt3 = 6.928203. The joint
contamination constant is **kappa_joint = 3 exactly**. The measurement factors are max over
targets 33.116168 and target sum 253.722222. The certificate is reproduced by
certify_d75.py; the constant by kappa_joint.py.

## 8. Corrected frontier and provenance

The structural-zero allowance is measured by decomposing each fixed-(route,scheme) path or
the pooled cube in the record's basis. The corrected joint bound is

    sum_T G_alpha |c_hat_alpha| + kappa_joint eta_Q + sum_T G_alpha delta_alpha,

with the measured reconstruction M_hat = max_a |sum_T c_hat_alpha Phi_alpha(a)|. The ratio
is the family-aware bound over the full-cube bound, figure + the measured non-additive part.

| design/lane | full non-add. | full bound | eta_Q max | kappa_joint | charge | M_hat | family bound | ratio |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| fused fp64 path | 2.9895e-16 | 2.5882e-14 | 2.7709e-16 | 4 | 1.1084e-15 | 2.9895e-16 | 2.6990e-14 | 1.0428 |
| separate fp64 path | 3.4467e-16 | 2.5928e-14 | 2.7709e-16 | 4 | 1.1084e-15 | 3.3367e-16 | 2.6990e-14 | 1.0410 |
| fused fp64 pooled | 2.0975e-14 | 4.6558e-14 | 3.9541e-16 | 3 | 1.1862e-15 | 2.0975e-14 | 4.6558e-14 | 1.000000 |
| separate fp64 pooled | 2.0903e-14 | 4.6486e-14 | 3.9541e-16 | 3 | 1.1862e-15 | 2.0903e-14 | 4.6486e-14 | 1.000000 |
| fused fp32 path | 2.0907e-7 | 3.5907e-7 | 1.7106e-7 | 4 | 6.8424e-7 | 2.4709e-7 | 1.0813e-6 | 3.0114 |
| separate fp32 path | 1.7106e-7 | 3.2106e-7 | 1.7106e-7 | 4 | 6.8424e-7 | 2.2808e-7 | 1.0623e-6 | 3.3088 |
| fused fp32 pooled | 2.2213e-7 | 3.7213e-7 | 2.3758e-7 | 3 | 7.1275e-7 | 2.6609e-7 | 1.1288e-6 | 3.0335 |
| separate fp32 pooled | 2.1456e-7 | 3.6456e-7 | 2.1771e-7 | 3 | 6.5313e-7 | 2.6609e-7 | 1.0291e-6 | 2.8228 |

The old row did use 67 eta_Q; the apparent 14 eta_Q increment was 67 minus the full cube's
own measured non-additive part of 53.05 eta_Q. With the joint constant the fp64 pathwise
inflation is 4.1-4.3 per cent, the pooled fp64 inflation is 0.000 per cent pointwise
(0.025 per cent with one global allowance), and the fp32 ratios are about 3. The retrospective
run of D75 reproduces the alias identity to 6.1e-29 (fp64) and 3.9e-22 (fp32), with
max|chat - c| 5.9e-17/8.9e-17 (fp64) and 5.7e-8/5.4e-8 (fp32); on fp64 the target estimates
are 34 times the allowance, so the 75-run design rejects additivity at 857/858 and 856/858
points, the same verdict as the 144-run cube.

**Provenance.** eta_Q is measured on the full cube that the smaller design would replace.
The retrospective check validates the algebra and the estimator relation; it does not
establish a prospective saving. A prospective use needs either a mathematical bound on the Q
components or amortised re-certification: measure the full cube once, use the smaller design
for new points or arguments under a drift condition that no error mass moves into Q (same
kernel revision, arithmetic tier and region, with the allowance refreshed when drift is
detected). The fp64 pathwise eta_Q = 2.77e-16 is the same order as the permitted interaction
signal there (1.33e-16); the pooled eta_Q = 3.95e-16 is 34 times below the permitted pooled
signal (1.33e-14). It is a measured noise or extraction envelope, not a proved physical
Taylor remainder. The fp32 eta_Q is at the interaction-signal scale, which is why the
reduced design is a valid but loose bound there and the full cube still decides the
verdict.

For a prospective operated protocol with fresh measurements, the measurement term enters
through the amplification (max over targets 8.485281 for D21 and 33.116168 for D75), so the
34-fold fp64 margin over the allowance can shrink. Any protocol that uses the reduced
design on new measurements must state and monitor that amplified floor.

The second-instance four-axis star design replaces the 8-run half fraction with 5 runs at
ratio 1.00, with contamination at most 3e-12 against B = 0.1964; its structural remainder is
exactly zero in that arithmetic.

## 9. Premise

The pathwise results do not depend on the packed-sum cone. The pooled 75-run count assumes
the packed-sum cone is {partition, pack} after the budget is fixed, as the source audit's
cone table states; if that cone widens to route or scheme, the pooled counts and design
change. Pathwise D21 is unaffected.

## 10. Artifacts and commands

The scripts and their outputs are in this directory: `certify_d21.py`, `certify_d75.py`,
`kappa_joint.py`, `corrected_frontier.py`, `reconcile_table.py`,
`pooled_d75_record_check.py`, `exact_proofs.py` and `second_instance_star.py`, with JSON
and log outputs. Record inputs: the composition-audit full-cube cell archives named in
`fullcube_summary.json`; measured kernel v5.0.0 revision `482fd338`; source read at
`v5.3.1`; harness tag `v1.0.42`. The record is deposited in pull request #59 and is tagged
`support-family-2026-11-22`; this derivation is included in the record.

Run from this directory:

```
python certify_d21.py
python certify_d75.py
python kappa_joint.py
python corrected_frontier.py
python reconcile_table.py
python pooled_d75_record_check.py
```
