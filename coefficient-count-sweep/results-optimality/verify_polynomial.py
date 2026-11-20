#!/usr/bin/env python3
"""Independent re-check of the polynomial-family pilot certificate.

Verifies:
  * degree-16 dVP: 18 alternating nodes with exact lower bounds > 5e-14 for
    the stored monomial coefficients;
  * degree-17 witness: exact grid error <= 5e-14 for the 40-digit Chebyshev
    coefficients and for their binary64 roundings, on the 2,001-node paper
    grid.

Usage: python -I verify_polynomial.py [certificate.json]
Exit 0 iff every check passes.
"""

from __future__ import annotations

import json
import os
import sys
import time
from fractions import Fraction

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import f0exact as fx  # noqa: E402

A = Fraction("11.899848152108484")
B = Fraction("28.989337738820740")
W = B - A
TARGET = fx.TARGET
N_TERMS = 120


def cheb_eval_exact(cs, t):
    if not cs:
        return Fraction(0)
    T0, T1 = Fraction(1), t
    total = cs[0]
    if len(cs) == 1:
        return total
    total += cs[1] * T1
    for j in range(1, len(cs) - 1):
        T2 = 2 * t * T1 - T0
        total += cs[j + 1] * T2
        T0, T1 = T1, T2
    return total


def exact_grid_poly(cs):
    worst_hi, worst_lo = Fraction(0), Fraction(0)
    for i in range(2001):
        t = Fraction(2 * i, 2000) - 1
        x = A + W * (t + 1) / 2
        lo, hi = fx.f0_bounds(x, n_terms=N_TERMS)
        p = cheb_eval_exact(cs, t)
        e_hi = max(abs(lo - p), abs(hi - p))
        e_lo = Fraction(0) if lo <= p <= hi else min(abs(lo - p), abs(hi - p))
        worst_hi = max(worst_hi, e_hi)
        worst_lo = max(worst_lo, e_lo)
    return worst_lo, worst_hi


def main(argv):
    path = argv[1] if len(argv) > 1 else \
        "certificate_polynomial_paper_cell.json"
    print(f"verifying {path}")
    cert = json.load(open(path, encoding="utf-8"))
    ok = True
    t0 = time.time()

    d16 = cert["degree_16"]
    P = [Fraction(s) for s in d16["P_monomial"]]
    Q = [Fraction(s) for s in d16["Q"]]
    nodes = d16["nodes"]
    needed = d16["needed_alternations"]
    if len(nodes) < needed:
        print(f"FAIL: {len(nodes)} nodes < {needed}")
        ok = False
    signs, delta = [], None
    for node in nodes:
        x = Fraction(int(node["x_num"]), int(node["x_den"]))
        tt = Fraction(int(node["t_num"]), int(node["t_den"]))
        lo, hi = fx.f0_bounds(x, n_terms=N_TERMS)
        r, _ = fx.rat_eval(P, Q, tt)
        dlo, dhi = lo - r, hi - r
        if dlo > 0 and dhi > 0:
            sign, mag = 1, dlo
        elif dlo < 0 and dhi < 0:
            sign, mag = -1, -dhi
        else:
            print(f"FAIL: indeterminate sign at x={x}")
            ok = False
            break
        if sign != node["sign"]:
            print(f"FAIL: stored sign mismatch at x={x}")
            ok = False
            break
        if mag <= TARGET:
            print(f"FAIL: |F_0-r|={float(mag):.3e} <= target at x={x}")
            ok = False
            break
        signs.append(sign)
        delta = mag if delta is None else min(delta, mag)
    alt = all(signs[i] != signs[i + 1] for i in range(len(signs) - 1))
    if not alt:
        print("FAIL: degree-16 signs do not alternate")
        ok = False
    if delta is not None:
        print(f"degree 16 dVP: nodes={len(nodes)} alternating={alt} "
              f"delta={float(delta):.6e} "
              f"({float(delta)/float(TARGET):.2f} x target)")

    d17 = cert["degree_17"]
    cs = [Fraction(s) for s in d17["chebyshev_coefficients"]]
    cs64 = [Fraction(float(s)) for s in d17["chebyshev_coefficients"]]
    for label, cc in (("40-digit stored", cs), ("binary64", cs64)):
        lo, hi = exact_grid_poly(cc)
        meets = hi <= TARGET
        print(f"degree 17 grid ({label}): max|err| in "
              f"[{float(lo):.6e}, {float(hi):.6e}] meets={meets}")
        ok &= meets

    conc = cert["conclusion"]
    if conc.get("polynomial_minimum_stored_count") == 18 and not (
            ok and conc.get("degree_17_meets")
            and conc.get("degree_16_certified_to_fail")):
        print("FAIL: conclusion inconsistent")
        ok = False
    print(f"{'PASS' if ok else 'FAIL'} in {time.time()-t0:.1f}s")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
