"""Bounded-precision rigorous interval evaluation of the two cells' routes.

This is the "higher-precision evaluation of the recurrence with a rigorous
error bound" tool, with a definite precision.  Intervals are fixed-point
integers at scale 2^-P; every operation is performed with outward rounding in
exact Python integers, so the result is rigorous by construction, with a
width that is explicit in P.  The shared exponential enters as the exact
rational enclosure of exact_arith.enclose_exp_neg_half, converted outward.

For each entry the script returns an interval for the route R and hence for
the rounding residual r = d - R, and tests the sign of the rounding
correction that the ordering needs:
    A1 absolute 2   r_c - r_o  (both errors negative; need > |t_c - t_o|)
    A2 absolute 23  r_c + r_o  (mixed sign; need < -|t_c + t_o|)
It reports the certified sign if the interval excludes zero, the interval
width, and the precision P at which the test first succeeds.

Run:  python -I fixed_point_cert.py
Writes: fixed_point_cert.json
"""
import csv
import json
import os
import sys
from fractions import Fraction
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(HERE))
import exact_arith as ea                     # noqa: E402

RECORD = REPO / "amplification-sweeps" / "certified-bound"
SWEEP = Path(os.environ.get("CERT_SWEEP", REPO / "amplification-sweeps" / "amplification-form"))
KERNEL = Path(os.environ.get("BOYS_ROOT", REPO.parent / "boys"))
sys.path.insert(0, str(RECORD / "code"))
import model                                 # noqa: E402


# ---------------------------------------------------------------------------
# fixed-point intervals: integer endpoints at scale 2^-P
# ---------------------------------------------------------------------------

def f_floor(fr, P):
    num = fr.numerator << P
    return num // fr.denominator


def f_ceil(fr, P):
    num = fr.numerator << P
    return -((-num) // fr.denominator)


def point(fr, P):
    return f_floor(fr, P), f_ceil(fr, P)


def from_float(x, P):
    return point(Fraction(x), P)


def from_fraction_interval(lo, hi, P):
    return f_floor(lo, P), f_ceil(hi, P)


def iadd(a, b):
    return a[0] + b[0], a[1] + b[1]


def isub(a, b):
    return a[0] - b[1], a[1] - b[0]


def imul(a, b, P):
    vals = (a[0] * b[0], a[0] * b[1], a[1] * b[0], a[1] * b[1])
    return min(vals) >> P, -((-max(vals)) >> P)


def idiv(a, b, P):
    assert b[0] > 0, "positive divisor expected on these routes"
    S = 1 << P
    lo = None
    hi = None
    for n in a:
        for d in b:
            v = n * S
            q0 = v // d
            q1 = -((-v) // d)
            lo = q0 if lo is None else min(lo, q0)
            hi = q1 if hi is None else max(hi, q1)
    return lo, hi


def iscale2(a):
    return 2 * a[0], 2 * a[1]


def width(a, P):
    return Fraction(a[1] - a[0], 1 << P)


# ---------------------------------------------------------------------------
# the routes (same algebraic expressions as exact_arith)
# ---------------------------------------------------------------------------

def map_iv(x, a, b, P):
    t1 = isub(x, a)
    t2 = iscale2(t1)
    den = isub(b, a)
    q = idiv(t2, den, P)
    return isub(q, from_float(1.0, P))


def cheb_iv(c, t, P):
    T1 = from_float(1.0, P)
    total = imul(from_float(c[0], P), T1, P)
    if len(c) == 1:
        return total
    T0 = t
    total = iadd(total, imul(from_float(c[1], P), T0, P))
    for k in range(2, len(c)):
        Tk = isub(imul(iscale2(t), T0, P), T1)
        total = iadd(total, imul(from_float(c[k], P), Tk, P))
        T1, T0 = T0, Tk
    return total


def route_product(x, band, n, P):
    a, b, c0 = model.piece(0, band)
    c = model.piece(n, band)[2]
    t = map_iv(from_float(x, P), from_float(a, P), from_float(b, P), P)
    return cheb_iv(c, t, P)


def route_ship(x, n, eh_iv, P):
    a, b, c = model.ext_piece()
    t = map_iv(from_float(x, P), from_float(a, P), from_float(b, P), P)
    f = cheb_iv(c, t, P)
    xv = from_float(x, P)
    for l in range(1, n + 1):
        p = imul(from_float(float(l) - 0.5, P), f, P)
        a_ = isub(p, eh_iv)
        f = idiv(a_, xv, P)
    return f


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def main():
    model.configure(kernel=KERNEL, exp_bits=RECORD / "exp" / "exp_bits.txt")
    harness = {}
    with open(RECORD / "results" / "harness_cells.csv", newline="") as h:
        for r in csv.DictReader(h):
            harness[(r["band"], r["metric"], int(r["n"]))] = r
    cells = {}
    with open(SWEEP / "results" / "cells.csv", newline="") as h:
        for r in csv.DictReader(h):
            cells[(r["band"], float(r["x"]), int(r["n"]))] = (
                float(r["recursed"]), float(r["direct"]))
    exp = model.load_exp()

    out = {"cells": []}
    for band, metric, n in (("A1", "absolute", 2), ("A2", "absolute", 23)):
        h = harness[(band, metric, n)]
        bi = 0 if band == "A1" else 1
        xs, xp = float(h["x_ship"]), float(h["x_prod"])
        d_o = cells[(band, xs, n)][0]
        d_c = cells[(band, xp, n)][1]
        rec = dict(cell=f"{band} {metric} {n}", precisions=[])
        print("=" * 82)
        print(f"{band} {metric} order {n}")
        for P in (64, 80, 96, 112, 128, 160):
            eh_lo, eh_hi, _, _ = ea.enclose_exp_neg_half(
                Fraction(xs), Fraction(1, 10 ** 62))
            eh_iv = from_fraction_interval(eh_lo, eh_hi, P)
            Rc = route_product(xp, bi, n, P)
            Ro = route_ship(xs, n, eh_iv, P)
            rc = isub(from_float(d_c, P), Rc)     # d - R interval
            ro = isub(from_float(d_o, P), Ro)
            corr = isub(rc, ro) if band == "A1" else iadd(rc, ro)
            corr_w = width(corr, P)
            # truncation threshold from the exact-rational certificate
            Fb_c = ea.enclose_boys(n, Fraction(xp))
            Fb_o = ea.enclose_boys(n, Fraction(xs))
            Bc, Dc = ea.product_affine(model, xp, bi, n)
            Bo, Do = ea.ship_affine(model, xs, n)
            Ehi = (Fraction(eh_lo), Fraction(eh_hi))
            Tc = ea.isub((Bc, Bc), (Fb_c[0], Fb_c[1]))
            Ro_exact = ea.iadd((Bo, Bo), ea.iscale(Ehi, Do))
            To = ea.isub(Ro_exact, (Fb_o[0], Fb_o[1]))
            t_test = ea.isub(Tc, To) if band == "A1" else ea.iadd(Tc, To)
            threshold = max(abs(t_test[0]), abs(t_test[1]))
            # certified sign of the correction
            sgn = (corr[0] > 0) - (corr[1] < 0)
            if band == "A1":
                ok = corr[0] > threshold     # r_c - r_o > |t_c - t_o|
                margin = Fraction(corr[0], 1 << P) - threshold
            else:
                ok = corr[1] < -threshold    # r_c + r_o < -|t_c + t_o|
                margin = -threshold - Fraction(corr[1], 1 << P)
            rec["precisions"].append(dict(
                P=P, corr_lo=float(Fraction(corr[0], 1 << P)),
                corr_hi=float(Fraction(corr[1], 1 << P)),
                corr_width=float(corr_w), sign=sgn,
                threshold=float(threshold),
                certifies=bool(ok),
                margin=float(margin),
                width_over_margin=(float(corr_w / margin) if margin > 0
                                   else None),
            ))
            print(f"  P={P:3d}: correction in "
                  f"[{float(Fraction(corr[0], 1 << P)):.6e}, "
                  f"{float(Fraction(corr[1], 1 << P)):.6e}] "
                  f"width {float(corr_w):.3e} "
                  f"sign {'+' if sgn > 0 else '-' if sgn < 0 else '0'} "
                  f"threshold {float(threshold):.3e} "
                  f"certifies={ok}")
        out["cells"].append(rec)

    with open(HERE / "fixed_point_cert.json", "w") as fh:
        json.dump(out, fh, indent=1)
    print("=" * 82)
    print("wrote fixed_point_cert.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
