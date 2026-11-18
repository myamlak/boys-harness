# Paper C optimality certificate deposit

Self-contained certificate artifacts for one cell of Paper C, "Stored
coefficient counts for the Boys function F_0: rational against Chebyshev
fits".  Nothing here is inserted into the paper.  This is a deposit for review.

## Case

* Function: `F_0(x) = sqrt(pi)/(2 sqrt(x)) erf(sqrt(x))`.
* Interval: the paper cell `[x0, x1]`, with
  `x0 = 11.899848152108484` and `x1 = 28.989337738820740`
  (prior-work width `W_prior = x1 - x0 = 17.089489586712256`).
* Target: `5e-14` absolute on the paper's 2,001-node certifying grid
  (2,000 in the record's `nodes` field, 2,001 grid points).
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
`1.023324e-14` with the 70-digit coefficients, `9.594201e-15` with the
binary64-stored coefficients, against the record's binary64 Horner value
`9.936496e-15`.  Hence 12 is both achievable and necessary.

At this cell the certified 12 also equals the off-diagonal theorem's even
count `2 ceil((L + log 4)/(2 mu)) = 12`, so theorem-backed existence and the
certificate's necessity meet.

### 2. The polynomial-family minimum is 18 stored coefficients

A degree-17 Remez polynomial meets `5e-14` with exact maximum error
`4.042273e-14` (70-digit Chebyshev coefficients) and `4.043317e-14`
(binary64 Chebyshev coefficients), so 18 stored coefficients suffice.  A dVP
certificate with a degree-16 reference and 18 exactly verified alternating
nodes gives the lower bound `2.212570e-13` (4.43× target), so no degree-16
polynomial meets.  The all-polynomial minimum is therefore exactly 18.

The paper's Chebyshev-interpolation count of 19 is the interpolation solver's
minimum, not the polynomial-family minimum: the degree-17 DCT interpolant
misses at `5.60e-14` in the pilot's independent DCT check, while the degree-18
interpolant meets.  The paper already labels this count solver-relative.

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
* The degree-17 binary64 Chebyshev coefficients meet under exact evaluation;
  the actual binary64 Clenshaw rounding is covered by a conservative budget,
  not by interval arithmetic.  The margin is about `0.96e-14`.
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

## Commands and exit codes

Run from this directory unless stated otherwise.  `python` is Python 3.14.7.

| command | exit | result |
| --- | --- | --- |
| `python -I inputs/extract_cells.py` | 0 | extracts the two one-cell inputs, source digests verified |
| `python -I make_certificate.py --audit inputs/paper_cell_abs_audit.json --corrected inputs/paper_cell_abs_corrected.json --out ../regen.json --log ../regen.log` | 0 | 124.8 s; regenerates the rational certificate; conclusion, all nine type bounds, the grid result and the two extras match the deposited certificate exactly (only the recorded input paths/digests and the timestamp differ) |
| `python -I verify_certificate.py` | 0 | PASS in 202.0 s; re-checks every exact inequality, the Sturm poles, the grid meeting and the conclusion |
| `python -I certify_polynomial.py` | 0 | Remez degrees 17 and 16, degree-16 dVP `2.212570e-13`, degree-17 grid checks; see `certificate_polynomial_run.log` |
| `python -I verify_polynomial.py` | 0 | PASS in 50.4 s; re-checks the 18 degree-16 nodes and the two degree-17 grid checks |
| `python -I alternation_stored_candidate.py --corrected inputs/paper_cell_abs_corrected.json --out alternation_stored_candidate.json` | 0 | 3.5 s; 13 alternating nodes; type-`(5,6)` dVP bound `9.050746e-15` |
| `python -I bench_high_degree.py` | 1 | negative scaling result: BRASIL `[15/12]` `converged=False`, deviation 6.0; `certify_type` `pole_free=False`, zero nodes |

The deposited `certificate_rational_paper_cell.log` records the authoritative
generation from the full Paper C records:
`python -I make_certificate.py`, exit 0 (166.7 s in that run).  Development
runs before the deposited ones fixed three tooling faults: a Fraction/string
conversion bug, Python's 4300-digit integer-string limit, and a verifier
cross-check tolerance.  They do not affect the deposited results; the final
commands above are the ones that produced and re-checked this deposit.

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
| `f0exact.py` | `f85114e00bc1506e062e9b8b7bf70341774ed7089d6adbfc924e8ddb2c0777cb` |
| `inputs/extract_cells.py` | `8adc2fc789ed23fa309f7a2aa56cbd74a208aee232a1f1c81b6a183b0a16f3f3` |
| `inputs/paper_cell_abs_audit.json` | `a1f46c56cc82480c36db6d948270be745e36b6d1a00958f1745b70936cb7a536` |
| `inputs/paper_cell_abs_corrected.json` | `b9b0ba70b321f909eba1d47313705cb9130eb4e3360cddd525b9a6e4060bb784` |
| `make_certificate.py` | `1b57288237ac2a2cabbf63edf6cf487d70ca2611172248fc7b203ac6e5b54ad5` |
| `verify_certificate.log` | `c4e2830fa7f8dbfcf08de8ed0b3a66cb5a28768c59826bb9a323115bbfde62a5` |
| `verify_certificate.py` | `bc5c2d570dbc093cd46b3b3842129499c3d7036a2a0a14c63bd8e2d16d89950e` |
| `verify_polynomial.log` | `9d37c16243f9ace9e65a3c3369eefd775025a04eec8c29894f2d2bc974b2c2d0` |
| `verify_polynomial.py` | `9ae214c7aff7b9723cdcfa8a6b833cd73b6aed69a1790cbfca92e85f75cc3687` |

## Proposed harness deposit

Repository `https://github.com/myamlak/boys-harness` (Public, BSD 3-Clause),
data tag `count-optimality-2026-11-18`, directory
`coefficient-count-sweep/results-optimality/`, holding this directory.  It is
a deposit, not a release, in the pattern of the existing
`count-audit-2026-11-16` audit deposit.  The harness working tree was not
touched and the tag is proposed, not created.
