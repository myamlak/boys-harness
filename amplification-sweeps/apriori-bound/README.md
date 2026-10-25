# An a priori error bound for the two Boys entries at the paper's cells — measured verdict

**Verdict: no.** With the stored coefficients and exp(−x) alone, no ordering can be certified at
the paper's cells: 0 of the 67 committed comparisons, 0 of all 132, 0 of the 66 band maxima, in
both experiments, under seven bound variants. The reason is structural, not a constant factor: at
the losing entry the delivered error is a median 9.6 times its own truncation error (more than 2×
at 62 of 67), so the ordering is carried by evaluation rounding. A sign-blind bound built from the
stored tables can charge that rounding only by magnitude, and a charge large enough to be correct
contains zero. Even an exact per-piece truncation oracle (accurate to 3.3e-23) certifies 0/67 with
the paper's arithmetic charge: at 16 of 67 the truncation ranks the entries the other way, and at
the other 51 the charge would have to be smaller by a median factor of 92 (paper convention) or
561 (strict scheme charge). The paper's "would need a per-piece truncation bound" is true but
incomplete: the missing object is a certified per-call rounding bound, not a sharper truncation
bound.

## Method (verified)

Identity: 2x F_n'(x) + (2n+1) F_n(x) = e^{−x}; for p with residual R = 2x p' + (2n+1)p − e^{−x}
and E = p − F_n, E(x) = ρE(a) + (1−ρ)/(2n+1)·[R_min,R_max], ρ=(a/x)^{n+1/2}; on [0,x], |E| ≤
max|R|/(2n+1). Four anchors: band 0 at a=0; band 1 at the seam b1 with Δ = p_{n,1}(b1)−p_{n,0}(b1)
exact plus the band-0 error; the extended seed at θ0 with δ = s(θ0)−p_{0,0}(θ0); and propagation
E_n = A_n E_0 upward, E_n = (A_n/A_32)E_32 downward (both checked algebraically and against the
kernel's bodies). Validation: F_n recomputed from the identity against the 60-digit reference
agrees to 3.30e-23 (second experiment) and 3.38e-23 (first).

Bounds from stored coefficients: Chebyshev l1 bounds max|p| ≤ Σ|c_k|, max|p'| ≤ h^{-1}Σk²|c_k|;
the consultant's symmetric envelope |R| ≤ 2b max|p'| + (2n+1)max|p| + 1; a sharper symmetric route
expands R itself in Chebyshev coefficients (the exponential via Bessel coefficients
e^{−m}[I_0(h)+2Σ(−1)^k I_k(h)T_k], tail beyond N = ⌈|h|⌉+80 bounded); a signed residual interval
from R's coefficients plus a sampled-with-slope-correction interval; charges: the paper's
3uΣ|c| per evaluation and 3u s_l A_n/A_l per step (s_l = |p_l|+envelope or 1/(2l+1)), and a strict
scheme charge ((3⌊d/2⌋+3)uΣ|c| scalar; (2+⌈log2(d+1)⌉)uΣ|c| product). Entries built by their
routing (downward fit, upward seed+recursion, per-order fit), relative metric via a certified F_n
denominator, and the comparison certifies only when one entry's magnitude upper edge is below the
other's lower edge.

## Measurements

- **Certifications: 0/67 committed, 0/132 all, 0/66 band maxima, both experiments, seven variants**
  (sampled100/400/1600 signed residual, bound step size, strict charge, whole-piece coefficient
  interval, crude envelope); five relative variants abstain 66/66 each. Sampling density is not the
  obstruction.
- **Looseness versus delivered error:** paper charge + sampled signed residual: median 10.8×
  (second), 6.3× (first); strict charge: median 25.2× / 27.4×; sharp symmetric: 13.6× / 9.2×. The
  paper's charge is not a certified bound: delivered error exceeds it at 7/132 and 19/132 entries
  (worst 1.30× at B1 order 15, product); the strict charge contains all 264 entries.
- **Truncation-only:** the coefficient envelope is a median 43.5× the exact truncation; the sampled
  signed interval a median 15.3×; minimum ≈ 1.
- **Why the lower bound is always zero:** at 67/67 the loser's charged lower bound is 0; its
  truncation interval straddles zero at 66/67; both entry intervals are one-signed at 0/132. The
  loser's delivered error ÷ its exact truncation: min 0.985, median 9.57, q90 48.4, max 133.1
  (>2× at 62/67). Structural argument: if |e_l| > 2|t_l| then any sound sign-blind charge C_l ≥
  |r_l| exceeds |t_l| and the interval contains zero; this holds at 62/67 and is independent of the
  charge formula, the sampling density and the truncation bound's sharpness.
- **Exact-truncation oracle:** truncation rank agrees with delivered rank at 51/67 and disagrees at
  16 (one of them unresolved at 5.97e-23, below twice the oracle's accuracy). Exact truncation with
  zero charge orders 51/67; with the paper's or strict charge, 0/67. The combined charge reduction
  needed at the 51 is a median 92× (paper) / 561× (strict). A sign-blind bound knowing the exact
  rounding magnitudes could order only 3/67.
- **One cell's detail:** B1 absolute order 0 is the only truncation interval one-signed; its charge
  is 42.8× the interval and the combined charge would need a 97.3× reduction; the product's
  delivered error is 20.9× its truncation there.
- **What would change it:** a certified per-call rounding bound with a known sign (higher-precision
  evaluation or error-free transformations — no longer a table-read bound); a rounding budget below
  the truncation (roughly 4× at q10, 9.6× at median, 48× at q90, 133× at the worst cell), which
  would change the measured winners and require re-measurement; or a different comparison (a lower
  bound on the gap is exactly what fails).

## Provenance

Kernel boys 5.1.0 headers (the CHANGELOG states the tables/entries/contract are those of 5.0.0:
commit 482fd33): boys_coefficients.hpp sha256 8ce5e171… (same as the paper's), boys_fit_core.hpp
40b8273a…, boys_fit_rational.hpp b1e79eed…, boys_impl.hpp 4a384f8d…, version.hpp 0e2bc8b0…. Sweeps
read from sweeps/amplification_form and sweeps/product_amplification. The harness is deposited at
github.com/myamlak/boys-harness, amplification-sweeps/apriori-bound/ (scripts, logs and the two
exact-truncation CSVs); all arithmetic is mpmath at 60 digits, stored doubles parsed exactly.

## Suggested paper text (adopted)

Prediction section: "A fully a priori variant was attempted and does not exist at double precision.
The stored coefficients do certify each fit's truncation error through the anchored
Green's-function identity, and the arithmetic can be charged by the convention above; but at the
cells where the model commits the losing entry's delivered error is dominated by its evaluation
rounding (a median 9.6 times its truncation error), and a bound that knows only the magnitude of
that rounding must include zero. Even an exact per-piece truncation error, computed from the
identity to 3·10⁻²³, leaves 16 of the 67 comparisons ordered the other way by truncation, and at
the remaining 51 the arithmetic charge would have to be smaller by a median factor of 92. A
genuine a priori ordering therefore needs a certified per-call rounding bound, not a sharper
stored truncation bound, and none is claimed here."

Limitations: "The score is an across-order consistency test on measured inputs, not a prediction
from stored table bounds. A fully a priori variant is not merely unimplemented: at these cells the
delivered ordering is carried by the evaluation rounding, which a bound built from the stored
coefficients can charge only by magnitude; a charge large enough to be correct then contains zero.
An exact per-piece truncation oracle plus the model's arithmetic charge certifies none of the 67
commitments, so storing a per-piece truncation bound would not change the answer. What would is a
certified per-call rounding bound, which the double-precision entries do not have."
