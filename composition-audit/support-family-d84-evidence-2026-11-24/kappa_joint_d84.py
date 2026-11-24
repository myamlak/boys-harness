"""Exact joint contamination constant kappa_joint for D84.

For beta in Q the residual is
    r_beta(a) = Phi_beta(a) - sum_{alpha in T} B_{alpha,beta} Phi_alpha(a),
and kappa_joint = max_{beta in Q} sup_{a in cube} |r_beta(a)| / G_beta.

D84's basis, design and alias matrix lie in Q(sqrt 2, sqrt 3), so the
whole computation is done in that field with the four-coordinate
arithmetic of `field.py` (the same module the D21 exact computation
uses).  The alias matrix is the exact inverse of the 84 x 84 permitted
block; no floating-point inverse enters.

The supremum is decided exactly.  For each beta and cell the script
forms the exact squared ratio r^2 / G^2 and compares it with a rational
square k^2 by the sign of k^2 G^2 - r^2.  The sign of a field element
a + b sqrt2 + c sqrt3 + d sqrt6 is decided by rational interval bounds
on sqrt2 and sqrt3, refined until the interval excludes zero (or the
element is exactly zero).  When a cell attains equality for k, that
beta's exact supremum is k.

The script also reports the exact alias of the planted channel
(component 8, support {pack, form}) into the targets, to settle whether
that alias is zero or merely at the rounding floor.

It writes `kappa_joint_d84.json` beside itself.

Usage:

    python -I kappa_joint_d84.py
"""

from __future__ import annotations

import json
import math
import os
import sys
import time
from fractions import Fraction
from math import isqrt

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import d84_fixture as fx             # noqa: E402
from field import F4, Mat, basis_rows, ONE, ZERO   # noqa: E402

pa = fx.pa
d75 = fx.d75

LEVELS = pa.LEVELS
N = pa.N
KAPPA_CANDIDATES = (Fraction(1), Fraction(2), Fraction(5, 2), Fraction(3),
                    Fraction(4))


# --------------------------------------------------------------------------
# Exact field helpers
# --------------------------------------------------------------------------


def exact_g(choices):
    """Exact G_alpha = product of the per-axis row sup norms."""
    value = ONE
    for count, choice in zip(LEVELS, choices):
        if count == 2:
            sup = F4(1)                      # rows [1,1], [1,-1]
        elif choice == 0:
            sup = F4(1)                      # mean row
        elif choice == 1:
            sup = F4(0, 0, 0, Fraction(1, 2))  # sqrt6/2
        else:
            sup = F4(0, 1)                   # sqrt2
        value = value * sup
    return value


def build_exact_phi():
    rows = [basis_rows(c) for c in LEVELS]
    phi = [[ZERO] * N for _ in range(N)]
    for a in range(N):
        # Decode the cell's level on each axis from its cell index.
        lv = tuple((a // s) % c for c, s in zip(LEVELS, pa.STR))
        for idx in range(N):
            v = ONE
            for axis, choice in enumerate(pa.CHOICES[idx]):
                row = rows[axis][choice]
                v = v * row[lv[axis]]
            phi[a][idx] = v
    return Mat(phi)


_sqrt_cache = {}


def sqrt_bounds(digits):
    """Rational bounds on sqrt2 and sqrt3 at 10^-digits."""
    if digits not in _sqrt_cache:
        scale = 10 ** digits
        s2 = isqrt(2 * scale * scale)
        s3 = isqrt(3 * scale * scale)
        _sqrt_cache[digits] = (
            Fraction(s2, scale), Fraction(s2 + 1, scale),
            Fraction(s3, scale), Fraction(s3 + 1, scale),
        )
    return _sqrt_cache[digits]


def interval(value, digits):
    """Interval [lo, hi] containing the real value of a field element."""
    lo2, hi2, lo3, hi3 = sqrt_bounds(digits)
    a, b, c, d = value.a, value.b, value.c, value.d

    def scaled(coef, lo, hi):
        return (coef * lo, coef * hi) if coef >= 0 else (coef * hi,
                                                         coef * lo)

    b_lo, b_hi = scaled(b, lo2, hi2)
    c_lo, c_hi = scaled(c, lo3, hi3)
    products = (lo2 * lo3, lo2 * hi3, hi2 * lo3, hi2 * hi3)
    if d >= 0:
        d_lo, d_hi = d * min(products), d * max(products)
    else:
        d_lo, d_hi = d * max(products), d * min(products)
    return a + b_lo + c_lo + d_lo, a + b_hi + c_hi + d_hi


def sign_exact(value, start_digits=40, max_digits=20000):
    """Exact sign of a + b sqrt2 + c sqrt3 + d sqrt6, or 0."""
    if value.is_zero():
        return 0
    digits = start_digits
    while digits <= max_digits:
        lo, hi = interval(value, digits)
        if lo > 0:
            return 1
        if hi < 0:
            return -1
        digits *= 2
    raise ArithmeticError(
        f"sign undecided for {value!r} at {max_digits} digits")


# --------------------------------------------------------------------------
# D84 exact design
# --------------------------------------------------------------------------


def main():
    started = time.perf_counter()
    cells = [int(c) for c in fx.D84]
    perm = [any(s <= {0, 1, 2, 4, 5} or s <= {0, 1, 2, 3}
                for s in [pa.SUPPORTS[i]]) for i in range(N)]
    P_IDX = [i for i in range(N) if perm[i]]
    Q_IDX = [i for i in range(N) if not perm[i]]
    T_IDX = [i for i in P_IDX if 2 <= pa.ORDERS[i] <= pa.W_P]
    print(f"cells {len(cells)}, permitted {len(P_IDX)}, complement "
          f"{len(Q_IDX)}, targets {len(T_IDX)}")
    if (len(cells), len(P_IDX), len(Q_IDX), len(T_IDX)) != (84, 84, 60, 75):
        raise SystemExit("unexpected D84 counts")

    print("building the exact 144-cell basis", flush=True)
    E = build_exact_phi()
    print(f"  exact vs float Phi max deviation: "
          f"{np.abs(E.to_float() - pa.PHI).max():.3e} "
          f"({time.perf_counter() - started:.1f} s)", flush=True)

    G = [exact_g(pa.CHOICES[i]) for i in range(N)]
    G_float = np.array([float(g) for g in G])
    print(f"  exact G vs float G max deviation: "
          f"{np.abs(G_float - pa.G_ALPHA).max():.3e}")

    F = E.submatrix(cells, range(N))
    FP = F.submatrix(range(len(cells)), P_IDX)
    FQ = F.submatrix(range(len(cells)), Q_IDX)
    invN = F4(1) / F4(len(cells))
    print("forming the exact Gram matrix", flush=True)
    OmPP = FP.transpose().matmul(FP).scale(invN)
    OmPQ = FP.transpose().matmul(FQ).scale(invN)
    print("inverting the 84 x 84 permitted block exactly", flush=True)
    invPP = OmPP.inverse()
    B = invPP.matmul(OmPQ)
    print(f"  exact inverse and alias matrix done "
          f"({time.perf_counter() - started:.1f} s)", flush=True)

    # Agreement with the floating-point fixture.
    d84 = fx.build("d84")
    B_float = B.to_float()
    print(f"  exact B vs float B max deviation: "
          f"{np.abs(B_float - d84['B']).max():.3e}")

    posT = [P_IDX.index(t) for t in T_IDX]
    print(f"computing {len(Q_IDX)} x {N} exact residuals", flush=True)
    rows = []
    worst_overall = None
    for j, beta in enumerate(Q_IDX):
        g2 = G[beta] * G[beta]
        best_ratio = -1.0
        best = None
        res2 = []
        for a in range(N):
            r = E.rows[a][beta]
            for tr, tcomp in zip(posT, T_IDX):
                r = r - B.rows[tr][j] * E.rows[a][tcomp]
            r2 = r * r
            res2.append(r2)
            ratio = math.sqrt(float(r2) / float(g2))
            if ratio > best_ratio:
                best_ratio = ratio
                best = a
        r2 = res2[best]
        exact_k = None
        for k in KAPPA_CANDIDATES:
            k2 = F4(k * k)
            if (r2 - k2 * g2).is_zero():
                exact_k = k
                break
        # Prove no cell exceeds the candidate.
        if exact_k is not None:
            bound2 = F4(exact_k * exact_k) * g2
            bad = sum(1 for s in res2 if sign_exact(bound2 - s) < 0)
            if bad:
                exact_k = None
        rows.append({
            "beta": int(beta),
            "support": [pa.AXES[i] for i in sorted(pa.SUPPORTS[beta])],
            "choices": [int(c) for c in pa.CHOICES[beta]],
            "G_exact": [str(G[beta].a), str(G[beta].b), str(G[beta].c),
                        str(G[beta].d)],
            "G_float": float(G[beta]),
            "sup_ratio_float": best_ratio,
            "argmax_cell": int(best),
            "argmax_levels": [int((best // s) % c)
                              for c, s in zip(LEVELS, pa.STR)],
            "exact_kappa": (int(exact_k) if exact_k is not None
                            and exact_k.denominator == 1
                            else (str(exact_k) if exact_k else None)),
            "ratio_squared_exact": [str(r2.a), str(r2.b), str(r2.c),
                                    str(r2.d)],
        })
        if exact_k is not None and (worst_overall is None
                                    or exact_k > worst_overall):
            worst_overall = exact_k
    rows.sort(key=lambda r: -r["sup_ratio_float"])
    print("\nexact kappa_joint per beta (largest 10 shown):")
    for row in rows[:10]:
        print(f"  beta {row['beta']:3d} support {row['support']} "
              f"G={row['G_float']:.6f} sup={row['sup_ratio_float']:.12f} "
              f"exact={row['exact_kappa']} at cell {row['argmax_cell']}")
    non_null = [r for r in rows if r["exact_kappa"] is not None]
    print(f"\n{len(non_null)}/{len(rows)} betas have an exact rational "
          f"supremum")
    values = sorted({float(Fraction(r["exact_kappa"]))
                     for r in non_null})
    print(f"exact kappa_joint values: {values}")

    # Exact planted-channel alias into the targets and the residual
    # components.  Component 8 is support {pack, form}.
    plant = 8
    j = Q_IDX.index(plant)
    alias_T = [B.rows[tr][j] for tr in posT]
    alias_T_float = [float(v) for v in alias_T]
    exact_zero_T = [v.is_zero() for v in alias_T]
    print(f"\nplanted channel {plant} support "
          f"{[pa.AXES[i] for i in sorted(pa.SUPPORTS[plant])]}:")
    print(f"  exact zero target aliases: {sum(exact_zero_T)}/"
          f"{len(exact_zero_T)}")
    print(f"  max |target alias| (float): "
          f"{max(abs(v) for v in alias_T_float):.3e}")
    nonzero = [(T_IDX[k], alias_T_float[k]) for k in range(len(posT))
               if not exact_zero_T[k]]
    print(f"  nonzero targets: {nonzero[:5]}")

    out = {
        "field": "Q(sqrt2, sqrt3)",
        "design": "D84",
        "counts": {"cells": len(cells), "permitted": len(P_IDX),
                   "complement": len(Q_IDX), "targets": len(T_IDX)},
        "exact_inverse_exists": True,
        "exact_alias_matrix": True,
        "exact_B_vs_float_max_deviation": float(
            np.abs(B_float - d84["B"]).max()),
        "rows": rows,
        "exact_kappa_joint_values": values,
        "exact_kappa_joint": (
            str(worst_overall) if worst_overall is not None else None),
        "all_betas_exact": len(non_null) == len(rows),
        "planted_channel": {
            "component": plant,
            "support": [pa.AXES[i] for i in sorted(pa.SUPPORTS[plant])],
            "target_aliases_exactly_zero": int(sum(exact_zero_T)),
            "target_count": len(posT),
            "max_abs_target_alias_float": float(
                max(abs(v) for v in alias_T_float)),
            "nonzero_targets": [
                {"component": int(T_IDX[k]),
                 "support": [pa.AXES[i]
                             for i in sorted(pa.SUPPORTS[T_IDX[k]])],
                 "alias": float(alias_T_float[k])}
                for k in range(len(posT)) if not exact_zero_T[k]],
        },
        "runtime_seconds": None,
    }
    out["runtime_seconds"] = time.perf_counter() - started
    with open(os.path.join(HERE, "kappa_joint_d84.json"), "w") as handle:
        json.dump(out, handle, indent=1, default=float)
    print(f"\nruntime {out['runtime_seconds']:.1f} s")
    exact4 = worst_overall is not None and worst_overall == 4
    print("kappa_joint = 4 exactly" if exact4
          else "kappa_joint not shown to be exactly 4")
    return 0 if exact4 else 1


if __name__ == "__main__":
    sys.exit(main())
