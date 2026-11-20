# Paper C optimality certificate deposit

Self-contained certificate artifacts for one cell of Paper C, "Stored
coefficient counts for the Boys function F_0: rational against Chebyshev
fits".  Nothing here is inserted into the paper.  This is a deposit for review.

## Correction 1 (binary64 parse and DCT record)

This is the first correction deposit.  Its parent is the immutable data tag
`count-optimality-2026-11-18`.  Two records change:

* `verify_certificate.py` now parses the stored binary64 coefficient strings
  as binary64 values (`Fraction(float(s))`), matching how the record stores
  and evaluates them, instead of as exact decimal literals.  The regenerated
  `verify_certificate.log` reports `9.769346e-15` for the binary64
  coefficients (exactly `9.76934599056004e-15`, the value recorded in
  `certificate_rational_paper_cell.json` under
  `stored_candidate_binary64.exact_grid_check_paper_W.max_error_upper`).
  The parent deposit's log reported `9.594201e-15`, the result of the
  decimal-literal parse; the paper quoted that number.
* `dct/` holds the independent DCT interpolation record:
  `dct/dct_interpolation_check.py`, its JSON output
  `dct/dct_interpolation_paper_cell.json`, and its log
  `dct/dct_interpolation_paper_cell.log`.  They reproduce the degree-17
  miss (`5.599782e-14`) and the degree-18 meeting (`9.879860e-15`) from
  F_0 at 60-digit precision on the same cell and grid as the certificate.
  The parent deposit recorded this comparison only as a README sentence.

The rational certificate JSON, the polynomial certificate JSON and their
generation logs are unchanged.  The parent tag and its commit are not
modified.

## Correction 2 (continuum upper bound and 40-digit precision)

This second correction's parent is the immutable data tag
`count-optimality-2026-11-18-fix`.  Three things change:

* The degree-17 upper bound is now certified over the whole interval, not
  only on the 2,001-node grid.  `continuum_bound.py` computes an exact
  rational interval Taylor enclosure; `continuum_bound_paper_cell.json`
  records the bounds `4.265856e-14` (stored coefficients) and `4.266879e-14`
  (binary64), both below `5e-14`; `verify_continuum_bound.py` recomputes and
  checks them.
* The polynomial certificate stores 40 significant digits
  (`mp.nstr(c, 40)`), not 70.  The 70 is the audit's decimal precision.  This
  README and Paper C now say 40 significant digits for the stored coefficient
  strings.
* `continuum_bound.py` now asserts that the certificate's cell equals the
  hard-coded cell, and `verify_polynomial.py` labels the stored strings
  "40-digit stored".  The deposited `certificate_polynomial_run.log` and
  `verify_polynomial.log` retain the historical "70-digit" label; that label
  describes the generation precision, not the stored format, and the logged
  values are unchanged.

The rational and polynomial certificate JSONs and their generation logs are
unchanged.  The parent tags and their commits are not modified.

## Case

* Function: `F_0(x) = sqrt(pi)/(2 sqrt(x)) erf(sqrt(x))`.
* Interval: the paper cell `[x0, x1]`, with
  `x0 = 11.899848152108484` and `x1 = 28.989337738820740`
  (prior-work width `W_prior = x1 - x0 = 17.089489586712256`).
* Target: `5e-14` absolute.  The counts are checked on the paper's 2,001-node
  certifying grid (2,000 in the record's `nodes` field, 2,001 grid points);
  the degree-17 polynomial's upper bound is also checked over the whole
  interval by the continuum bound.
* Stored rational candidate: the corrected record's nominal `[5/6]`
  truncation, 12 stored coefficients (6 numerator + 7 denominator).  Its
  full witness arrays are 7/7; the nominal truncation is the reported
  representation because it meets the target.
* The audit's width field is `17.089489586712254`, one ulp below
  `x1 - x0`.  The exact grid check is reported on both intervals; the
  certificate nodes lie inside the paper interval.

## Exact statements proved

The proofs are computer-assisted: every inequality below was checked in exact
rational arithmetic.  `F_0` is enclosed by its alternating Taylor series with
an exact rational remainder bound, and every reference denominator is checked
pole-free with an exact Sturm sequence.  The only non-formal ingredients are
the classical de la Vallee Poussin (dVP) argument and the use of Python's
`fractions` arithmetic.

### 1. The rational stored count is exactly 12

No rational function of type `[m/k]` with `m+1+k <= 11` meets `5e-14`
uniformly on `[x0, x1]`.  The proof certifies every type at stored count 11
(`m+k = 10`, `k = 1..9`, the paper's admissible set) with `m+k+2 = 12`
exactly verified alternating nodes, and the search-excluded endpoint types
`[10/0]` and `[0/10]` separately.  Every lower-count pair is embedded by
padding the numerator into one of these types.

| type | stored | certified lower bound | × target |
| --- | --- | --- | --- |
| `[1/9]` | 11 | 1.380585e-11 | 276.12 |
| `[2/8]` | 11 | 1.276121e-11 | 255.22 |
| `[3/7]` | 11 | 1.121254e-12 | 22.43 |
| `[4/6]` | 11 | 4.812850e-13 | 9.63 |
| `[5/5]` | 11 | 4.073395e-13 | 8.15 (weakest) |
| `[6/4]` | 11 | 6.779747e-13 | 13.56 |
| `[7/3]` | 11 | 2.482697e-12 | 49.65 |
| `[8/2]` | 11 | 1.393663e-11 | 278.73 |
| `[9/1]` | 11 | 1.377386e-11 | 275.48 |
| `[10/0]` (extra) | 11 | 3.847030e-09 | 76940.60 |
| `[0/10]` (extra) | 11 | 1.245162e-10 | 2490.32 |

The stored 12-coefficient fit meets the target: exact maximum error
`1.023324e-14` with the 70-digit coefficients, `9.769346e-15` with the
binary64-stored coefficients (the strings evaluated as binary64 values,
exactly `9.76934599056004e-15`), against the record's binary64 Horner value
`9.936496e-15`.  Hence 12 is both achievable and necessary.

At this cell the certified 12 also equals the off-diagonal theorem's even
count `2 ceil((L + log 4)/(2 mu)) = 12`, so theorem-backed existence and the
certificate's necessity meet.

### 2. The polynomial-family minimum is 18 stored coefficients

A degree-17 Remez polynomial meets `5e-14` over the whole cell.  Exact
evaluation on the 2,001-node grid gives `4.042273e-14` (stored
40-significant-digit Chebyshev coefficients) and `4.043317e-14` (binary64).
The exact rational interval Taylor bound in `continuum_bound.py` gives
`4.265856e-14` and `4.266879e-14` over `[x0, x1]`, leaving at least
`7.33e-15`.  A dVP certificate with a degree-16 reference and 18 exactly
verified alternating nodes gives the lower bound `2.212570e-13` (4.43×
target), so no degree-16 polynomial meets.  The all-polynomial minimum is
therefore exactly 18.

The paper's Chebyshev-interpolation count of 19 is the interpolation solver's
minimum, not the polynomial-family minimum: the independent DCT check in
`dct/` gives a degree-17 miss at `5.599782e-14` and a degree-18 meeting at
`9.879860e-15`, on the same cell and 2,001-node grid.  The paper already
labels this count solver-relative.

### 3. The stored candidate's own type (step (a))

The stored `[5/6]` error curve has the required 13 = `m+k+2` alternating
extrema.  The exact dVP bound at those nodes is `9.050746e-15`; the
candidate's exact sup there is `9.770200e-15` (ratio 1.079491).  The fit is
near-minimax but not exactly equioscillating, so the equality form of the
alternation theorem does not certify it as the exact type-`(5,6)` minimax.
This step bounds the fit's own type and does not by itself bound any stored
count.

## Numerical evidence only

* The reference fits are numerical objects.  The stored-11 references are the
  60-digit audit's own default-pass fits (`P70`/`Q70` truncated to the
  nominal type); the endpoint references are pilot BRASIL 2.1.2 solves.
  The certificate is exact for the reference it names; it does not depend on
  how the reference was found.
* The node locations are found numerically.  The certificate verifies the
  nodes exactly; the search only proposes them.
* The stored fit's near-equioscillation (the 7.9% gap in item 3) is numerical
  evidence about the fit, not a proof of type minimality.
* The interval bound is for the coefficient vector, evaluated exactly.  It
  does not by itself bound a particular evaluation route's rounding; the
  binary64 Clenshaw route is a separate delivery question.  The relevant
  margin for the coefficient vector is `7.33e-15`.
* The high-degree scaling probe is numerical: a default BRASIL `[15/12]` solve
  at the ship width returned `converged=False`, deviation 6.0, and a
  denominator pole on the interval; the certificate code rejected it
  (`pole_free=False`, zero nodes).  See the scale limits below.

## Scope limits

* One cell only: the paper cell at `5e-14` absolute.  No other width, target
  or function is certified.
* No rate law, no Stahl modulus, no diagonal constant `C_r = 4` is proved.
  The certificate concerns the actual minimum error at this cell.
* No charged-cell minimality.  The four charged cells (reported counts 21 and
  29) need references at count 20 or 28 (total degree up to 27), which the
  audit does not provide.  The default high-degree BRASIL probe failed as
  described above, so the failing step there is feasible pole-free
  reference construction, not the dVP argument or the exact arithmetic.
  A charged count may also be non-minimal if a lower-degree fit can be padded
  into the count-below type.
* The certificate is computer-assisted, not a formal proof in a proof
  assistant.  The dVP theorem, exact rational arithmetic, the alternating
  series remainder bound and Sturm's theorem are trusted mathematical
  ingredients.

## Method

For a rational reference `r = P/Q` of type `(m,k)`, if `F_0 - r` alternates in
sign at `m+k+2` points with magnitude at least `delta`, then every rational
`R` of type `(m,k)` has uniform error at least `delta`: `R - r` would change
sign `m+k+1` times but its numerator has degree at most `m+k`.  The affine map
between `x` and the reference's `t` variable preserves the type.  The
reference denominators are pole-free on the interval by an exact Sturm count.
`F_0` values are enclosed by the exact alternating series
`F_0(x) = sum (-1)^k x^k/(k! (2k+1))`; the first omitted term bounds the tail
once the terms decrease (checked in code).  The series formula is checked
formally through order 120 via `(2k+1) c_k = e_k`, the coefficient form of
`2x F_0'(x) + F_0(x) = exp(-x)`.

For the degree-17 polynomial the upper bound is an exact rational interval
Taylor expansion.  On 300 subintervals, each of width `2/300` in `t`, the
error `e(t) = p(t) - F_0(x(t))` is expanded to degree 5 about the left
endpoint.  The coefficients are intervals:
`e^(k)(t_i) = p^(k)(t_i) - (W/2)^k (-1)^k F_k(x_i)`, with `F_0` from its
alternating series and `F_1..F_6` from the Boys recurrence.  The degree-5
polynomial is sampled at 101 points per subinterval; the bound adds the
between-sample mean-value slack and the degree-6 Taylor remainder.  Every
inequality is exact rational arithmetic.  `verify_continuum_bound.py`
recomputes the whole bound and checks that the published upper decimal is at
least the recomputed value and below the target.

## Commands and exit codes

Run from this directory unless stated otherwise.  `python` is Python 3.14.7.

| command | exit | result |
| --- | --- | --- |
| `python -I inputs/extract_cells.py` | 0 | extracts the two one-cell inputs, source digests verified |
| `python -I make_certificate.py --audit inputs/paper_cell_abs_audit.json --corrected inputs/paper_cell_abs_corrected.json --out ../regen.json --log ../regen.log` | 0 | 124.8 s; regenerates the rational certificate; conclusion, all nine type bounds, the grid result and the two extras match the deposited certificate exactly (only the recorded input paths/digests and the timestamp differ) |
| `python -I verify_certificate.py > verify_certificate.log 2>&1; echo VERIFY_EXIT=$? >> verify_certificate.log` | 0 | PASS in 189.4 s; re-checks every exact inequality, the Sturm poles, the grid meeting and the conclusion; the corrected binary64 parse reports `9.769346e-15`, meets=True |
| `python -I certify_polynomial.py` | 0 | Remez degrees 17 and 16, degree-16 dVP `2.212570e-13`, degree-17 grid checks; see `certificate_polynomial_run.log` |
| `python -I verify_polynomial.py` | 0 | PASS in 50.4 s; re-checks the 18 degree-16 nodes and the two degree-17 grid checks |
| `python -I continuum_bound.py --certificate certificate_polynomial_paper_cell.json --n-sub 300 --out continuum_bound_paper_cell.json` | 0 | exact rational interval Taylor bound; `4.265856e-14` stored and `4.266879e-14` binary64; see `continuum_bound_paper_cell.log` |
| `python -I verify_continuum_bound.py --certificate certificate_polynomial_paper_cell.json --bound continuum_bound_paper_cell.json` | 0 | PASS; recomputes the exact bound and checks the published upper decimal is at least it and below `5e-14` |
| `python -I dct/dct_interpolation_check.py > dct/dct_interpolation_paper_cell.log 2>&1; echo DCT_EXIT=$? >> dct/dct_interpolation_paper_cell.log` | 0 | PASS; degree-17 DCT miss `5.599782e-14`, degree-18 meeting `9.879860e-15`, 60-digit values with an 80-digit stability check; writes `dct/dct_interpolation_paper_cell.json` and the redirected log |
| `python -I alternation_stored_candidate.py --corrected inputs/paper_cell_abs_corrected.json --out alternation_stored_candidate.json` | 0 | 3.5 s; 13 alternating nodes; type-`(5,6)` dVP bound `9.050746e-15` |
| `python -I bench_high_degree.py` | 1 | negative scaling result: BRASIL `[15/12]` `converged=False`, deviation 6.0; `certify_type` `pole_free=False`, zero nodes |

The deposited `certificate_rational_paper_cell.log` records the authoritative
generation from the full Paper C records:
`python -I make_certificate.py`, exit 0 (166.7 s in that run).  Development
runs before the deposited ones fixed three tooling faults: a Fraction/string
conversion bug, Python's 4300-digit integer-string limit, and a verifier
cross-check tolerance.  They do not affect the deposited results; the final
commands above are the ones that produced and re-checked this deposit.  This
correction deposit fixes a fourth, verifier-only fault: the binary64
coefficient strings were parsed as exact decimal literals instead of binary64
values; see "Correction to the parent deposit" above.

## Input provenance

* Paper C 60-digit audit:
  `C:\Users\myaml\dev\papers-c-v5\papers\coefficient-counts\results-corrected\corrected-counts-60dig.json`,
  SHA-256 `18a41a42798f0cd1fead4a143309ded3112d1bba94727b093bf27df02e3256e7`
  (the digest the supplement names).
* Paper C corrected BRASIL record:
  `C:\Users\myaml\dev\papers-c-v5\papers\coefficient-counts\results-corrected\corrected-counts.json`,
  SHA-256 `9bb8ac0633b10e8a7922684ca8a419135103f2c7825b4d06b9c71603399b7195`.
* The extracted one-cell copies are `inputs/paper_cell_abs_audit.json`
  (SHA-256 `a1f46c56...a536`) and `inputs/paper_cell_abs_corrected.json`
  (SHA-256 `b9b0ba70...b784`); `inputs/extract_cells.py` regenerates them
  and refuses on a source-digest mismatch.

## Software

Python 3.14.7, mpmath 1.4.1, numpy 2.5.2, scipy 1.18.1, baryrat 2.1.2.  The
exact certificate uses only the Python standard library `fractions`; mpmath
is used for node location and for bug-check cross-checks, never for a
certified inequality.

## Artifact digests

SHA-256, computed over the files as deposited.  `SHA256SUMS` covers every
file including this README; the README cannot list its own digest.

| file | SHA-256 |
| --- | --- |
| `alternation_stored_candidate.json` | `7018bcf3934989fef0f9e222198fa6024890c58f082557df9e6aaab8237a6026` |
| `alternation_stored_candidate.log` | `96c5db4b8fb310fab326ba189759680a6d301d63b33d78dbd61a6af2c205ac94` |
| `alternation_stored_candidate.py` | `b95dd192145af03c4ebb4a000e3d0f6a3a2ced5a61b36cf083049f20326e7494` |
| `bench_high_degree.log` | `a93633f4c376cad8b0a89eb42adb2f27be28320be9e949ff9aaa57d2032490b4` |
| `bench_high_degree.py` | `159eb9c3deb85145f4772039eb9cd6999bf8e9042eb3303e31ac7404753511d9` |
| `certificate_polynomial_paper_cell.json` | `b398c185f73c10142a7082ba344ea8dbae8145e74fa9ad1615c3af49b3c93792` |
| `certificate_polynomial_run.log` | `de139edd4df4472d62881deac0b13bd1d7732e0448d56fafc7334c67797a376a` |
| `certificate_rational_paper_cell.json` | `a637b30929623e4651910c0e397e6dc81d85b9c611732127a5a53bdf8629870e` |
| `certificate_rational_paper_cell.log` | `fe85442c720154aa4878bf11268ec63caae7b95e4f4a11f08cb299e68ce097dc` |
| `certify_polynomial.py` | `7ebe9cbfcce30199ff094f35b2b2d4fb02768d0cfc1bb6c0216ded17a9bb65c4` |
| `continuum_bound.py` | `eca008a6987dce19baa3cf552dc971ceccbe247fc57539cdc1423be336335519` |
| `continuum_bound_paper_cell.json` | `014ed7246f19454b535e0adf85cc1db42349cd5037c5592f7adda179fa15b4b7` |
| `continuum_bound_paper_cell.log` | `a017bd35330dc265da768d6261e8937a497f10906072e95973b0a8af0626885a` |
| `continuum_bound_verify.log` | `2d5271aa26a93b57eb4f497dd963a6ab5c7a5a5d6fc7909621d28f04347e298a` |
| `dct/dct_interpolation_check.py` | `7def3973cb1299b8283d75356382e84ebc1b5941050ea9ee9e8ecf36bab399fe` |
| `dct/dct_interpolation_paper_cell.json` | `db494d9231afc02a051af81a88d13e5948c9dbf146a504e4c35befcdf4531d98` |
| `dct/dct_interpolation_paper_cell.log` | `681d6bdda159344e98f40e9388de6c2d97f8a8427c22eaf30ba710ca88994413` |
| `f0exact.py` | `f85114e00bc1506e062e9b8b7bf70341774ed7089d6adbfc924e8ddb2c0777cb` |
| `inputs/extract_cells.py` | `8adc2fc789ed23fa309f7a2aa56cbd74a208aee232a1f1c81b6a183b0a16f3f3` |
| `inputs/paper_cell_abs_audit.json` | `a1f46c56cc82480c36db6d948270be745e36b6d1a00958f1745b70936cb7a536` |
| `inputs/paper_cell_abs_corrected.json` | `b9b0ba70b321f909eba1d47313705cb9130eb4e3360cddd525b9a6e4060bb784` |
| `make_certificate.py` | `1b57288237ac2a2cabbf63edf6cf487d70ca2611172248fc7b203ac6e5b54ad5` |
| `verify_certificate.log` | `f131dfc8ac5e667d289cce25081a6026f71d9bbf6973806c3a18015b2909b352` |
| `verify_certificate.py` | `3a2bef1bbbb71b9fec315bffb1301f21fd4b442b5179457a858a6d0dc1b0f495` |
| `verify_continuum_bound.py` | `f0f75abe3f374a925ce36fc6962d4b9477a3bf2b6961e001b265a71c35c6b07b` |
| `verify_polynomial.log` | `9d37c16243f9ace9e65a3c3369eefd775025a04eec8c29894f2d2bc974b2c2d0` |
| `verify_polynomial.py` | `209e83c4504dfe99035ca46be9d171564209fa5a6b34e85da8a9da64ce5b4512` |

## Proposed harness deposit

Repository `https://github.com/myamlak/boys-harness` (Public, BSD 3-Clause),
data tag `count-optimality-2026-11-20` (parent `count-optimality-2026-11-18-fix`),
directory `coefficient-count-sweep/results-optimality/`, holding this
directory.  It is a data deposit, not a release, in the pattern of the
existing `count-audit-2026-11-16` audit deposit.  The parent tags and their
commits are not modified.
