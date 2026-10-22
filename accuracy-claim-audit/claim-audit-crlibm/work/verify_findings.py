"""Independent verification of the two measured defects.

Nothing here uses work/ref.py for the decision: the expected values come from
mpmath at 60 digits (a different implementation and different algorithms from
MPFR), and from exact IEEE-754 reasoning about the overflow region.
"""
import sys
import os
import struct
import math

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import run_crlibm
import ref

import mpmath
mpmath.mp.dps = 400

DBL_MAX = float.fromhex("0x1.fffffffffffffp+1023")
TWO1024 = mpmath.mpf(2) ** 1024
RNDN_MID = mpmath.mpf(2) ** 1024 - mpmath.mpf(2) ** 970


def hex2d(h):
    return struct.unpack(">d", bytes.fromhex(h))[0]


def d2hex(x):
    return struct.pack(">d", x).hex()


def ieee_answer(exact, mode):
    """Correctly rounded double for an exact mpmath value, IEEE rules."""
    m = mpmath.mpf
    if abs(exact) <= m(DBL_MAX):
        # in range: use MPFR-free rounding via mpmath -> float then adjust
        f = float(exact)                      # nearest (Python rounds correctly)
        if mode == "rn":
            return f
        exact_f = m(f)
        if mode == "ru":
            return f if exact_f >= exact else math.nextafter(f, math.inf)
        if mode == "rd":
            return f if exact_f <= exact else math.nextafter(f, -math.inf)
        if mode == "rz":
            if exact_f == exact:
                return f
            return math.nextafter(f, 0.0)
    neg = exact < 0
    a = abs(exact)
    # IEEE 754-2019 7.4: overflow gives an infinity when the rounding
    # direction points away from zero, else the largest finite magnitude
    if mode == "ru":
        return -DBL_MAX if neg else math.inf
    if mode == "rd":
        return -math.inf if neg else DBL_MAX
    if mode == "rz":
        return -DBL_MAX if neg else DBL_MAX
    return (math.inf * (-1 if neg else 1)) if a >= RNDN_MID else (
        -DBL_MAX if neg else DBL_MAX)


def check_sinh():
    print("== sinh overflow, directed rounding ==")
    xs = [0x408633ce8fb9f87d, 0x408633ce8fb9f87e, 0x408633ce8fb9f87f,
          0x4087000000000000, 0x4090000000000000, 0x40a0000000000000,
          0x7fefffffffffffff]
    queries, label = [], []
    for h in xs:
        for m in ("rn", "ru", "rd", "rz"):
            queries.append(("sinh", m, d2hex(hex2d("%016x" % h))))
            label.append((h, m))
    got = run_crlibm.evaluate(queries)
    nbad = 0
    for (h, m), g in zip(label, got):
        x = hex2d("%016x" % h)
        if abs(mpmath.mpf(x)) > 1000:
            exp_val = mpmath.mpf(2) ** (x / math.log(2))  # not used; see below
        exact = mpmath.sinh(mpmath.mpf(x))
        want = ieee_answer(exact, m)
        wanthex = d2hex(want)
        flag = "" if wanthex == g else "  <-- DISAGREES"
        if wanthex != g:
            nbad += 1
        print("  x=%016x sinh_%s: crlibm=%s  ieee(expected)=%s%s"
              % (h, m, g, wanthex, flag))
    return nbad


def check_atanpi():
    print("== atanpi near 1/2, directed rounding ==")
    xs = [0x7fefffffffffffff, 0x7feffffffffffffe, 0x4350000000000000,
          0x43c0000000000000, 0x4400000000000000, 0x4360000000000000]
    queries, label = [], []
    for h in xs:
        for m in ("rn", "ru", "rd", "rz"):
            queries.append(("atanpi", m, "%016x" % h))
            label.append((h, m))
    got = run_crlibm.evaluate(queries)
    nbad = 0
    for (h, m), g in zip(label, got):
        x = mpmath.mpf(hex2d("%016x" % h))
        exact = mpmath.atan(x) / mpmath.pi
        want = ieee_answer(exact, m)
        wanthex = d2hex(want)
        if wanthex != g:
            nbad += 1
            print("  x=%016x atanpi_%s: crlibm=%s  ieee(expected)=%s  <-- DISAGREES"
                  % (h, m, g, wanthex))
        else:
            print("  x=%016x atanpi_%s: crlibm=%s  ok" % (h, m, g))
    return nbad


def check_atanpi_boundary():
    """Smallest x for which atanpi_rd is wrong = where 1/pi*atan(x) > 0.5-2^-54."""
    lo, hi = 1.0, float.fromhex("0x1p+60")
    t = mpmath.mpf(1) / 2 - mpmath.mpf(2) ** -54
    # solve atan(x)/pi = t  ->  x = tan(pi*t)
    xb = mpmath.tan(mpmath.pi * t)
    print("== atanpi boundary ==")
    print("  exact boundary x* = tan(pi*(1/2 - 2^-54)) =", mpmath.nstr(xb, 20))
    print("  as double:", float(xb).hex())
    print("  documented proven interval is [tan(2^-25 pi), tan(2^-5 pi)] =",
          mpmath.nstr(mpmath.tan(mpmath.mpf(2) ** -25 * mpmath.pi), 12), "..",
          mpmath.nstr(mpmath.tan(mpmath.mpf(2) ** -5 * mpmath.pi), 12))
    return 0


if __name__ == "__main__":
    a = check_sinh()
    print("sinh disagreements:", a)
    b = check_atanpi()
    print("atanpi disagreements:", b)
    check_atanpi_boundary()
