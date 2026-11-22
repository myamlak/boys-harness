This is a copy of the pre-registration and its amendment for the clean-edge
ladder. The pre-registration is dated 2026-11-22 and was written at 37 of 40
rows; its amendment was made the same day at 38 of 40 rows, before the last two
rows (the relative `1e-15` `ship` and `paper` rows) were read. The text below is
unchanged.

---

# Paper C clean-edge decision rule — pre-registered before the last three rows

2026-11-22. Written while the ladder has 37 of 40 rows; the missing rows are the relative
1e-15 target at widths 2^5, ship and paper. This rule is fixed before those rows are read.

## Question

Does the rational family's failure of the affine target law persist at the clean edge
a = 62, where the transient factor 1 - g_1 is below every target, or is it specific to the
x0 transient regime?

## Rule

Use the same analysis as the paper's x0 ladder: for each width, fit the stored count against
the target exponent on a free line over the contract range, and compute the residual and the
second difference.

- **Chebyshev control**: the affine target law holds at a width if its free-line residual is
  at most one stored coefficient. The Chebyshev family must hold at all five widths for the
  clean-edge test to be interpretable.
- **Rational failure at a width**: the free-line residual exceeds the two-coefficient
  rounding band, or the second difference is at least 2 in absolute value.
- **Clean-edge verdict**:
  - *Failure confirmed at the clean edge* if the rational family fails at at least three of
    the five widths while Chebyshev holds at every width. The title and abstract may then
    state that the affine target law fails for rational fits at the clean edge, not only at
    x0.
  - *Transient-specific failure* if the rational family holds at at least four of the five
    widths. The title and headline must be scoped to the x0 transient regime, and the
    clean-edge result becomes a boundary of the failure.
  - *Mixed* otherwise: report the exact widths and residuals; the headline claims only what
    the pattern supports, and the title keeps the scoped "where the affine target law fails"
    without a clean-edge generalisation.

## Reporting

The 40 rows, residuals, second differences and the verdict are recorded in the paper's
supplement and the review ledger. No row is dropped. If the 1e-15 wide-width solves fail to
converge within the stated caps, that is reported as a solver limitation, not counted as a
pass or a failure.

## Why pre-register

The clean-edge result decides whether the paper's title and headline generalise beyond x0.
Fixing the rule before the last rows protects the conclusion from post-hoc choice of width,
target or tolerance.

## Amendment, before the last two rows (2026-11-22, 38/40)

The ladder's full target range is 10^-8 to 10^-15, and the deep targets are where the
rational departure is expected to show. The primary fit therefore uses all eight targets
(the ladder's full range), with the free-line residual and second difference computed on
those eight counts. The paper's narrower contract range 10^-8 to 10^-13 is reported as a
secondary fit with the same thresholds, not used to choose the verdict. Thresholds and the
three-way verdict rule are unchanged. This amendment is recorded while the 1e-15 ship and
paper rows are still missing.
