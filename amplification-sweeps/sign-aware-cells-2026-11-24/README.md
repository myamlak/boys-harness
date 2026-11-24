# Sign-aware certificate at two amplification reversal cells

This directory deposits the certificate that resolves the sign of the rounding
correction at two of Paper A's amplification cells. At both cells the record's
sign-blind rounding charge is too wide to decide the ordering of the delivered
errors. The truncation centre says one entry has the larger error magnitude;
the delivered errors, once their rounding is taken into account, order the
other way. Exact sign-aware arithmetic decides the sign, and therefore the
ordering, at both cells.

The two cells are A1 absolute order 2 and A2 absolute order 23. This is a
two-cell certificate, not a population theorem. The scope section states what
it does and does not cover.

The *chosen* entry is the direct/product lane at `x_prod`; the *other* entry is
the recursed/ship lane at `x_ship`, following the `state` column of the
certified-bound record's `results/harness_cells.csv`. For an entry with
delivered double `d`, exact route value `R` and true Boys value `F_n`, the
delivered error is `e = d - F_n`, split as `e = t + r` with truncation
`t = R - F_n` and rounding `r = d - R`.

## What is certified

### A1 absolute order 2

Chosen: direct/product at `x = 8.254041852680183e-06`. Other: recursed/ship at
`x = 1.1679101370108986`.

- Delivered errors: `e_c = -5.023109722036405e-17` and
  `e_o = -7.804759823427790e-17`, both negative.
- Certified error ordering: `e_c - e_o` lies in
  `[2.781650101391385e-17, 2.781650101391385e-17]`. It is positive, so
  `e_c > e_o`; with both errors negative, `|e_c| < |e_o|`. The record's
  same-sign test is therefore certified at this cell.
- Certified rounding correction: `r_c - r_o` lies in
  `[3.139723003095558e-17, 3.139723003095558e-17]`, positive.
- Truncation centre: `t_c - t_o` lies in
  `[-3.580729017041734e-18, -3.580729017041734e-18]`. It is negative, so the
  truncation says the chosen magnitude is the larger. The certified rounding
  correction is larger than this separation and reverses the ordering: the
  total error says the chosen magnitude is the smaller.
- The reversal is recorded as `overturns_truncation_centre: true`.

### A2 absolute order 23

Chosen: direct/product at `x = 6.3134226241640246`. Other: recursed/ship at
`x = 10.865022762913437`.

- Delivered errors: `e_c = 8.454164305640834e-20` (positive) and
  `e_o = -5.537929113818018e-19` (negative), mixed sign.
- Certified error ordering: `e_c + e_o` lies in
  `[-4.692512683253935e-19, -4.692512683253935e-19]`. It is negative, so for
  opposite signs `|e_c| < |e_o|`. The signed separation `e_c - e_o` is also
  certified positive.
- Certified rounding correction: `r_c + r_o` lies in
  `[-4.693109582884852e-19, -4.693109582884852e-19]`, negative.
- Truncation centre: `t_c + t_o` lies in
  `[5.968996309173406e-23, 5.968996309173406e-23]`. It is positive, so the
  truncation says the chosen magnitude is the larger. The certified correction
  is larger than this sum and reverses the ordering: the total error says the
  chosen magnitude is the smaller.
- The reversal is recorded as `overturns_truncation_centre: true`.

At both cells the certified correction is the sign-aware quantity that the
record's sign-blind charge cannot supply.

## The certificate constructions

Four constructions are deposited. Each is computed from the stored binary64
coefficients, the arguments and the delivered doubles; none takes the
deposited error columns as a certificate input.

1. **Exact rational intervals** (`exact_arith.py`, `certify_cells.py`). Every
   stored coefficient, argument and mapping constant is a binary64 value,
   hence a dyadic rational. A route value is therefore an affine function of
   the shared exponential, `R = B + D * eh` with `eh = exp(-x)/2`, and `B`, `D`
   exact `Fraction`s; `D = 0` on the direct lane. The values of `exp(-x)/2`
   and `F_n` are enclosed by the alternating exponential and Boys series with
   the alternating-remainder theorem, in exact rational arithmetic. The
   delivered double is replayed with IEEE binary64 operations in the kernel's
   order and checked bit for bit against `model.py` and the delivered cells.
   The resulting intervals for `t`, `r` and `e` are exact-rational. The
   certified separation intervals have widths `1.651e-65` at A1 and
   `3.492e-64` at A2; the record's pair charge is `3.9e49` and `2.1e46` times
   those widths.

2. **EFT composition** (`eft_compose.py`). The operation graph is replayed
   carrying, at every node, both the exact real value at the delivered inputs
   and the delivered double. Each local rounding residual is checked against an
   error-free transformation: `TwoSum` and `TwoDiff` for addition and
   subtraction, `TwoProd` for multiplication, `TwoProd` plus `TwoSum` for the
   fused multiply-add, and the S18.4 identity for the refined division. The
   node values compose exactly, so the final residual is the negative of the
   fixed-input rounding residual. The refined-division identity holds at every
   call: 2 calls at A1 and 23 at A2. The local residuals cancel strongly:
   `sum |local residual| / |net rounding|` is 22.8 and 20.8 at A1 and `1.82e5`
   and `4.98e3` at A2 for the chosen and other entries. A one-sided per-node
   half-ulp band would be `5.34e-15` at A1 and `6.84e-15` at A2, 192 and
   `1.46e4` times the net margin, so a naive EFT composition cannot decide
   these cells; the exact composition can.

3. **Double-double** (`eft_compose.py`). The same route is evaluated in
   double-double and compared with the exact rational route. The differences
   are `2.32e-34` and `5.76e-33` at A1 and `1.98e-37` and `9.25e-36` at A2 for
   the chosen and other entries, many orders below the `2.78e-17` and
   `4.69e-19` margins. The conservative double-double bounds are `5.29e-29`,
   `4.58e-29`, `4.81e-29` and `9.55e-29`; the observed errors sit far inside
   them.

4. **Fixed-point intervals** (`fixed_point_cert.py`). The route is evaluated
   with integer endpoints at scale `2^-P` and outward rounding, for
   `P = 64, 80, 96, 112, 128, 160`. A1 certifies at `P = 64`; A2 first
   certifies at `P = 80`, because its `P = 64` correction interval straddles
   zero. At `P = 160` the correction intervals agree with the exact-rational
   values: `3.139723003095558e-17` at A1 and `-4.693109582884852e-19` at A2.
   The interval widths at `P = 160` are `2.81e-47` and `1.25e-45`.

**Agreement.** All four constructions give the same correction sign at both
cells. Every fixed-point interval at every tested precision contains the
exact-rational correction. The double-double evaluation matches the exact
route far inside the margin. The exact-rational certificate is the narrow one;
the others are independent evaluations of the same route and act as
cross-checks.

## Containment checks

`certify_cells.py` computes 24 checks, 12 per cell: four replay and model bit
checks, two comparisons of the certified route midpoint with the record's
70-digit route (tolerance `1e-40`), and six containment checks of the record's
measured error, truncation and rounding columns for the chosen and other
entries. The containment checks allow the record's `1e-45` tier margin. All 24
pass (`checks_pass: 24`, `checks_total: 24` in `cell_certificate.json`).

The measured columns are evaluated at 100 decimal digits. The largest
difference between a certificate midpoint and the record's measured value is
`9.05e-62` at A1 and `1.52e-62` at A2, consistent with the 60-significant-digit
reference's own rounding. `check_reference.py` also reproduces the record's
own sign-blind signed-residual test: at both cells it intersects (neither
above nor below, and no same-sign or magnitude decision), matching the 0/16
outcome of the certified-bound deposit.

## The sign-only shortfall

The record's rounding charge is sign-blind. Its pair charge `C_c + C_o` is
`6.452827e-16` at A1 and `7.416298e-18` at A2, 23.2 and 15.8 times the net
margins `2.781650e-17` and `4.692513e-19`. Knowing only the sign of each
rounding error and not its magnitude, the deposit's sign-known,
magnitude-blind figure halves that charge: `3.226413e-16` at A1 and
`3.708149e-18` at A2, still 11.6 and 7.9 times the margin. The sign-only bound
therefore falls short of certifying the ordering by factors **11.6** and
**7.9**.

The exact sign-aware correction clears the truncation threshold at both
cells: `r_c - r_o = 3.139723e-17` against `3.580729e-18` at A1, 8.77 times the
threshold and leaving a net margin of `2.781650e-17`; and
`r_c + r_o = -4.693110e-19` against `5.968996e-23` at A2, 7862 times the
threshold and leaving `4.692513e-19`. Its sign is certified, and with it the
error ordering.

## Inputs

The certificate reads the two amplification records as deposited in this
repository:

| input | role |
| --- | --- |
| `amplification-sweeps/certified-bound/code/model.py` | route model and operation sequence; the replay target |
| `amplification-sweeps/certified-bound/code/exactroute.py` | 70-digit route values and the 60-digit reference loader |
| `amplification-sweeps/certified-bound/exp/exp_bits.txt` | recorded `exp(-x)` bits, including the two shared `exp(-x)/2` values |
| `amplification-sweeps/certified-bound/results/harness_cells.csv` | the two cells' `x_ship`, `x_prod` and delivery state |
| `amplification-sweeps/certified-bound/results/cells_cert.csv` | record `C`, `eps` and measured error columns |
| `amplification-sweeps/amplification-form/results/cells.csv` | delivered recursed and direct doubles at the cells |
| `amplification-sweeps/amplification-form/data/reference.csv` | 60-significant-digit reference values |

The certificate scripts load the kernel's stored coefficient tables and
thresholds from a `boys` checkout at tag **`v5.2.5`**, commit
`8b7a289dfce7558c8498e89c17c05f20fa7bb7c8`, the revision named in
`certified-bound/provenance.json`. The amplification-form record's own
provenance names `v5.0.0` for its sweep, with `boys_coefficients.hpp` SHA-256
`8ce5e17157a4be5f621d695f377094d212cb839732e278d47f23aa225b4d9204`; the
certified-bound provenance records that this header differs from the `v5.2.5`
header only in the rational region-B seed, which the modelled policy does not
read, and that the coarse Chebyshev table the policy evaluates is byte-identical
in the two trees. The delivered cells are used as deposited.

## Environment

The recorded run used Windows 11 Pro 10.0.26200, CPython 3.14.7 (64-bit,
MSC v.1944) and mpmath 1.4.1. The scripts use `math.fma`, so Python 3.13 or
newer is required. No C++ compiler is needed.

The recorded commands use `python -I`. On the recording machine mpmath is
installed in the interpreter's own site-packages, so isolated mode finds it.
Isolated mode hides a user-site install; run without `-I` if `mpmath` is only
available there.

## Commands

From this directory, the recorded run is:

```
python -I cell_data.py          exit 0
python -I exact_arith.py        exit 0
python -I certify_cells.py      exit 0
python -I eft_compose.py        exit 0
python -I fixed_point_cert.py   exit 0
python -I check_reference.py    exit 0
```

`run_all.log` records the six commands and exit codes. The per-script logs are
`final_cell_data.log`, `final_exact_arith.log`, `final_certify_cells.log`,
`final_eft_compose.log`, `final_fixed_point_cert.log` and
`final_check_reference.log`. `certify_cells.py` writes `cell_certificate.csv`
and `cell_certificate.json`; `eft_compose.py`, `fixed_point_cert.py` and
`check_reference.py` write `eft_compose.json`, `fixed_point_cert.json` and
`reference_check.json`. `cell_data.py` prints the cell data, and
`exact_arith.py` runs its self-tests.

The deposited scripts are the recorded scripts with one change: the three
absolute input paths are resolved from the script's location and from the
`BOYS_ROOT` and `CERT_SWEEP` environment variables, so a fresh clone can run
them. The arithmetic is unchanged. Rerunning all six commands against the
records in this repository and a `v5.2.5` kernel reproduced every deposited
JSON, CSV and log file byte for byte.

To reproduce, check out the kernel beside the harness repository (the scripts
default to `<repo>/../boys`) and run the six commands above:

```
git clone --branch v5.2.5 https://github.com/myamlak/boys.git ../boys
```

Set `BOYS_ROOT` if the checkout is elsewhere, and `CERT_SWEEP` if the
amplification-form record is elsewhere.

## Scope and assumptions

- **Two cells, not a population theorem.** The certificate covers A1 absolute
  order 2 and A2 absolute order 23 only. It says nothing about the other cells
  or orders, and it does not turn the record's sign-blind bound into a
  sign-aware bound for the grid.
- **The record's H1-H5 assumptions remain.** The certificate is conditional on
  the same hypotheses as the certified-bound record:
  - **H1 (mode).** Binary64 with round-to-nearest, no flush-to-zero, no
    denormals-are-zero, and no overflow at the tested cells.
  - **H2 (no compiler contraction).** The evaluated arithmetic is the model's
    sequence: a separate multiply and subtract in the recurrences and a
    separate coefficient product in the product form.
  - **H3 (exponential).** The C library's `exp` satisfies
    `|exp_hat(x)/exp(x) - 1| <= 2^-52`; the record probes 9,541 arguments and
    the two reference arguments are among them.
  - **H4 (reference).** The reference differs from the true `F_n` by at most
    `5e-60` relative and absolute.
  - **H5 (route).** The 70-digit route differs from the true route by less
    than `1e-59`.
- The exact-rational and fixed-point enclosures do not use H4 or H5; those
  enter the comparisons with the record's measured columns and the
  paper-framework combination. The delivered-value replay uses H1 (binary64),
  H2 (the modelled operation sequence) and H3 (the recorded libm exponential
  bits).
- The certificate does not re-measure the kernel or the delivered cells. It
  uses the deposited delivered doubles and the stored coefficient tables. A
  different kernel table, compiler build or libm would need its own replay.
- The sign-only shortfall is a statement about this record's charge at these
  two cells, not a claim that the record's bound is unsound. The record's
  bound is sound but sign-blind; the certificate resolves the sign with exact
  arithmetic.
- The deposit is tagged `sign-aware-cells-2026-11-24`; that tag carries this
  README, the scripts, the outputs, the logs and the input records at the
  versions used.

## Files

| file | bytes | role |
| --- | ---: | --- |
| `cell_data.py` | 4268 | loads the two cells; prints arguments, delivered bits, model agreement and record columns |
| `exact_arith.py` | 16079 | exact rational arithmetic, EFT self-tests, route affine forms, delivered replay, transcendental enclosures |
| `certify_cells.py` | 17063 | builds the exact-rational certificate and the 24 checks; writes the certificate CSV and JSON |
| `eft_compose.py` | 21130 | EFT composition, S18.4 division identity and double-double cross-check; writes `eft_compose.json` |
| `fixed_point_cert.py` | 8048 | fixed-point outward-rounded intervals at `P = 64` to 160; writes `fixed_point_cert.json` |
| `check_reference.py` | 10028 | reference consistency, the record's sign-blind test and the paper-framework combination; writes `reference_check.json` |
| `cell_certificate.csv` | 6150 | per-cell certified intervals, signs, margins and check results |
| `cell_certificate.json` | 9473 | the same certificate with the 24 checks and the sign-only factors |
| `eft_compose.json` | 17755 | per-node EFT counts, cancellation, division identities, double-double agreement and pair margins |
| `fixed_point_cert.json` | 4077 | correction intervals at each fixed-point precision |
| `reference_check.json` | 8905 | measured columns at 100 digits, containment, certified signs and the paper-framework combination |
| `run_all.log` | 249 | the six commands and exit codes of the recorded run |
| `final_cell_data.log` | 2679 | per-script log of `cell_data.py` |
| `final_exact_arith.log` | 735 | per-script log of `exact_arith.py` |
| `final_certify_cells.log` | 4848 | per-script log of `certify_cells.py` |
| `final_eft_compose.log` | 3093 | per-script log of `eft_compose.py` |
| `final_fixed_point_cert.log` | 1668 | per-script log of `fixed_point_cert.py` |
| `final_check_reference.log` | 1614 | per-script log of `check_reference.py` |

SHA-256 of every deposited file except this README:

```
228550823769cbc5520797c7914e3335789dfc983684d8610ca7fa19f83d9672  cell_data.py
9e0b5d586368b22c3e911126c3aec93fb11f1bfb38f51ab4bb9862c04311fe46  exact_arith.py
3e05db18997552048274904a63a64e9e3556386b71b09a41b71825dd26287c8b  certify_cells.py
102ca4ef9105634036adf89602469d2af1904671bcd765f9afab70b8480a268e  eft_compose.py
804f40051349a62b0aec6869cf9d9861e56287eab26eb5dfdf6ba95a8bb5bb02  fixed_point_cert.py
fadc88253fb9ccb88c3ee77f3af8f86d67cc7e3705c22e3d736d426f9c1bdd20  check_reference.py
4c45f1d3600d3900eb79ba1dabd084629188274af04f254e1590f16fe29d859d  cell_certificate.csv
827605d11d6b4531337f2f4ed86f43f162ce6919c8fd6b71755bbf239350e03f  cell_certificate.json
46aa02ca4e50a3c5daa0cd34ec13ec8f4d41217f0dccb77454ab7569d9cf6ee9  eft_compose.json
e4f8878a3284c38045f3c8f2f4d11a93121f62deccc3e8b0bab4cc8ae867fbb0  fixed_point_cert.json
4f4575a30ee72142af7c372ee0cb1d5b71ba05415c7239382d72700a116e31ae  reference_check.json
813aeb7fbaf86fb188157ccdb5b2572edb4ac33aa163af26be35398ab470f579  run_all.log
80d1665e8252e981f36c2b3d4e94a1f2ebee59defd1af01b7bd241db2fad2974  final_cell_data.log
3794a5aa13fe34f039bb190f37c1ecc742eeb5d60ce14daca6c6b8d606530009  final_exact_arith.log
d2f79e5431f3458a497f498ede452372b079a5bda727a5c7a4cf20483ba4a86e  final_certify_cells.log
87e3d7178ac5d6f987d7117d14b877caaee4bcdf0557f6fa33baf5e413e0a921  final_eft_compose.log
299a324055277a0b80e3a1fd2d23326ac014a479a946069913043ebb2a63cbca  final_fixed_point_cert.log
8737bbd3cc9732fea298f5b79baccaf7725846e2069f435cee3fa98a501fb409  final_check_reference.log
```

The whole directory is about 152 KiB.
