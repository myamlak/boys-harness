#!/usr/bin/env python3
"""Exact-rational arithmetic for the Paper C optimality pilot.

The certificate this module supports uses only:
  * F_0(x) = sum_{k>=0} (-1)^k x^k / (k! (2k+1)), an alternating series whose
    remainder after N terms is bounded by the first omitted term once the
    terms decrease in magnitude (checked in code);
  * a rational reference r(t) = P(t)/Q(t) evaluated by Horner in exact
    Fraction arithmetic.

Nothing here trusts binary floating point for a certified number.  mpmath is
used only to locate candidate extrema; every inequality in a certificate is
rechecked with Fractions.

The formula is checked formally in `series_ode_check`: the power series
coefficients satisfy (2n+1)c_n = e_n, where e_n are the coefficients of
exp(-x); equivalently 2x F_0'(x) + F_0(x) = exp(-x), with F_0(0) = 1, which
follows from F_0'(x) = -F_1(x) and 2x F_1(x) + F_0(x) = exp(-x).
"""

from __future__ import annotations

import sys
from fractions import Fraction
from math import factorial

# Exact intervals can carry integers longer than Python's default 4300-digit
# int/str conversion limit; all inputs here are local and trusted.
sys.set_int_max_str_digits(0)

# Target of the Paper C cell: 5 x 10^-14 absolute.
TARGET = Fraction(5, 10**14)


def decimal_fraction(s: str) -> Fraction:
    """Exact Fraction from a decimal string, including exponent notation."""
    return Fraction(s)


def f0_bounds(x: Fraction, n_terms: int = 120):
    """Rigorous rational interval [lo, hi] containing F_0(x) for x > 0.

    F_0(x) = sum_{k=0}^{n-1} (-1)^k x^k/(k!(2k+1)) + R_n,
    |R_n| <= a_n = x^n/(n!(2n+1)), provided a_{k+1} <= a_k for every k >= n.
    The ratio a_{k+1}/a_k = x(2k+1)/((k+1)(2k+3)) is checked < 1 at k = n and
    at k = n+1; the ratio is decreasing for k >= 2, so this suffices.
    """
    if x <= 0:
        raise ValueError("x must be positive")
    term = Fraction(1)      # a_0
    total = Fraction(1)     # k = 0 term
    for k in range(1, n_terms):
        term = term * x * (2 * k - 1) / (k * (2 * k + 1))
        total += -term if k % 2 else term
    a_n = term * x * (2 * n_terms - 1) / (n_terms * (2 * n_terms + 1))
    a_n1 = a_n * x * (2 * n_terms + 1) / ((n_terms + 1) * (2 * n_terms + 3))
    if not (a_n1 < a_n < term):
        raise ArithmeticError("alternating-series decrease check failed")
    return total - a_n, total + a_n


def f0_closed_form_mp(x):
    """F_0 through the closed form, for cross-checks only (mpmath mpf)."""
    import mpmath as mp
    return mp.sqrt(mp.pi) / (2 * mp.sqrt(x)) * mp.erf(mp.sqrt(x))


def rat_eval(P, Q, t: Fraction):
    """Exact Horner evaluation of P(t)/Q(t).  Returns (value, denominator)."""
    num = Fraction(0)
    for c in reversed(P):
        num = num * t + c
    den = Fraction(0)
    for c in reversed(Q):
        den = den * t + c
    if den == 0:
        raise ZeroDivisionError("reference denominator vanishes")
    return num / den, den


def dvp_check_mapped(P, Q, nodes_xt, target: Fraction = TARGET,
                     n_terms: int = 120):
    """Verify the de la Vallee Poussin inequalities at mapped nodes.

    `nodes_xt` is an ordered list of (x, t) pairs: x is the point in the
    function's argument, t is the variable of the reference P(t)/Q(t)
    (an affine image of x).  The dVP argument applies because the affine map
    preserves the numerator and denominator degrees.

    Returns (ok, delta, rows).  Each row carries x, t, the exact rational
    sign, and an exact rational lower bound on |F_0(x) - r(x)|.
    """
    rows = []
    for x, t in nodes_xt:
        lo, hi = f0_bounds(x, n_terms=n_terms)
        r, den = rat_eval(P, Q, t)
        dlo, dhi = lo - r, hi - r
        if dlo > 0 and dhi > 0:
            sign, mag = 1, dlo
        elif dlo < 0 and dhi < 0:
            sign, mag = -1, -dhi
        else:
            rows.append(dict(x=decimal_string(x), t=decimal_string(t),
                             sign=0, lower="0", indeterminate=True))
            return False, Fraction(0), rows
        rows.append(dict(x=decimal_string(x), t=decimal_string(t),
                         x_num=str(x.numerator), x_den=str(x.denominator),
                         t_num=str(t.numerator), t_den=str(t.denominator),
                         sign=sign, lower=decimal_string(mag),
                         lower_num=str(mag.numerator),
                         lower_den=str(mag.denominator)))
    signs = [row["sign"] for row in rows]
    alternate = all(signs[i] != signs[i + 1] for i in range(len(signs) - 1))
    delta = min(decimal_fraction(row["lower"]) for row in rows)
    return bool(alternate and delta > target), delta, rows


def poly_eval(coeffs, x: Fraction) -> Fraction:
    """Horner evaluation of a polynomial given ascending coefficients."""
    acc = Fraction(0)
    for c in reversed(coeffs):
        acc = acc * x + c
    return acc


def poly_trim(coeffs):
    """Drop trailing zero coefficients; return at least [0]."""
    out = list(coeffs)
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out


def sturm_root_count(coeffs, lo: Fraction, hi: Fraction):
    """Number of distinct real roots of a rational polynomial in (lo, hi).

    Sturm's theorem, exact rational arithmetic.  The caller must check the
    endpoints separately (poly_eval(coeffs, lo) and at hi are not zero).
    """
    p = poly_trim(coeffs)
    if len(p) <= 1:
        return 0
    seq = [p, [i * p[i] for i in range(1, len(p))]]
    while True:
        a = seq[-2]
        b = seq[-1]
        if len(b) <= 1 and b[0] == 0:
            break
        # remainder of a divided by b, then negate
        r = list(a)
        db = len(b) - 1
        while len(r) - 1 >= db and not (len(r) == 1 and r[0] == 0):
            d = len(r) - 1
            factor = r[-1] / b[-1]
            shift = d - db
            for i, c in enumerate(b):
                r[shift + i] -= factor * c
            r = poly_trim(r)
            if len(r) == 1 and r[0] == 0:
                break
        if len(r) == 1 and r[0] == 0:
            break
        seq.append([-c for c in r])
        if len(seq) > 200:
            raise ArithmeticError("Sturm sequence did not terminate")
    def variations(x):
        signs = []
        for p in seq:
            v = poly_eval(p, x)
            if v != 0:
                signs.append(1 if v > 0 else -1)
        return sum(1 for i in range(len(signs) - 1) if signs[i] != signs[i + 1])
    return variations(lo) - variations(hi)


def decimal_string(q: Fraction, digits: int = 40) -> str:
    """Truncated decimal rendering of q at `digits` fractional digits.

    The rendering is exact to the stated digits (it truncates toward zero);
    for a positive rational it is a lower bound.
    """
    if q == 0:
        return "0"
    sign = "-" if q < 0 else ""
    q = abs(q)
    scaled = q.numerator * 10**digits // q.denominator
    whole, frac = divmod(scaled, 10**digits)
    return f"{sign}{whole}.{frac:0{digits}d}"


def series_ode_check(order: int = 80):
    """Formal power-series check of the F_0 formula.

    With F_0 = sum c_k x^k and exp(-x) = sum e_k x^k, verify
    (2k+1)c_k = e_k coefficientwise through x^order, and c_0 = 1.
    Exact rationals; no floating point.
    """
    c = [Fraction((-1) ** k, factorial(k) * (2 * k + 1))
         for k in range(order + 1)]
    e = [Fraction((-1) ** k, factorial(k)) for k in range(order + 1)]
    if c[0] != 1:
        return False, 0
    for n in range(order + 1):
        if (2 * n + 1) * c[n] != e[n]:
            return False, n
    return True, order
