# Corrected rational counts

These are the records the paper's rational counts now come from. They were
produced by an independent solver, not by the kernel's own Remez and Lawson
implementations.

The two files are copied verbatim from the verification run of 2026-10-31:

| file | sha256 |
| --- | --- |
| `corrected-counts.csv` | `e11899412881e23c65756c7949a2e2f9564c8e30ccf33c397b6a5549ba2bc335` |
| `corrected-counts.json` | `9bb8ac0633b10e8a7922684ca8a419135103f2c7825b4d06b9c71603399b7195` |

## What produced them

The solver is `baryrat` 2.1.2, the BRASIL rational minimax algorithm. It is
independent of the kernel's Remez exchange and Lawson iteratively reweighted
least-squares solvers. The target function is `F_0` on `[x0, x0 + W]` with
`x0 = 11.899848152108484`; the rational form is `p(t) / (1 + q1 t + ... + qk
t^k)` in `t = 2(x - x0)/W - 1`, and the stored count is `(m + 1) + k` with
`k <= 12`.

## The acceptance rule, as applied

The paper's own counting rule was used unchanged. A fit is accepted when its
coefficients, converted exactly to monomials and rounded to binary64, evaluated
by Horner in binary64 by both the fused and the separate routes, hold the
worst error over the 2,001-node grid in the target's kind: relative
`1e-13`, absolute `5e-14`, or relative `1e-14`. The error is the worse of the
two routes. The denominator is guarded at `min |den| >= 1e-3 max |den|`. The
2,001-node grid certifies; a 4,001-node grid checks. The reference is `F_0` at
70 decimal digits.

A count is therefore the smallest stored count at which an exhibited
admissible fit meets, under this solver. It is an upper bound, not a proof of
minimality. One ladder cell is higher than the paper's: `2^5` at `10^-15` is
25 `[22/2]` here against the paper's 24 `[21/2]`, whose coefficients are not in
the record and so could not be re-evaluated. That cell is not independently
reproduced.
