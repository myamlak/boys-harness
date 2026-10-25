"""Certified a priori error bounds for the region-A Boys fits, from the
stored double coefficients alone.

Everything here is derived from three things and nothing else:

  * the stored Chebyshev coefficients (exact binary doubles),
  * the differential equation 2x F_n'(x) + (2n+1) F_n(x) = exp(-x),
  * exp(-x), evaluated at working precision.

No high-precision Boys value enters any bound.
"""

import csv
import math
from mpmath import mp, mpf, chebyt, diff, quad, exp, besseli, inf

import coeffs

mp.dps = 60

TABLE = coeffs.load()

BAND_NAMES = ("A1", "A2")


# --------------------------------------------------------------------------
# piece selection and evaluation, exactly as the kernel does it
# --------------------------------------------------------------------------

def piece_index(order, x):
    return coeffs.find_piece(TABLE, order, mpf(x) if not isinstance(x, mpf) else x)


def piece_of(order, band_index):
    """The piece an order reads inside a band (0 = A1, 1 = A2)."""
    return TABLE["piece_start"][order] + band_index


def coef_list(index):
    p = TABLE["pieces"][index]
    return [mpf(c) for c in TABLE["coeffs"][p["off"] : p["off"] + p["deg"] + 1]]


def mapping(piece):
    """x = m + h t on t in [-1, 1]."""
    a, b = mpf(piece["a"]), mpf(piece["b"])
    return (a + b) / 2, (b - a) / 2


def t_of(piece, x):
    m, h = mapping(piece)
    return (mpf(x) - m) / h


def poly_value(c, t):
    """Clenshaw summation of sum c_k T_k(t) at working precision."""
    b1 = mpf(0)
    b2 = mpf(0)
    two_t = 2 * t
    for ck in reversed(c[1:]):
        b1, b2 = two_t * b1 - b2 + ck, b1
    return t * b1 - b2 + c[0]


def poly_deriv_coeffs(c):
    """Chebyshev coefficients of dp/dt."""
    d = len(c) - 1
    g = [mpf(0)] * d
    for k in range(1, d + 1):
        ck = c[k]
        if k % 2 == 0:
            for j in range(0, (k - 2) // 2 + 1):
                g[k - 1 - 2 * j] += 2 * k * ck
        else:
            for j in range(0, (k - 3) // 2 + 1):
                g[k - 1 - 2 * j] += 2 * k * ck
            g[0] += k * ck
    return g


def poly_deriv_value(c, t):
    return poly_value(poly_deriv_coeffs(c), t)


# --------------------------------------------------------------------------
# the residual and its Chebyshev coefficients
# --------------------------------------------------------------------------

def residual_cheb(piece_index_, order, n_terms=100):
    """Chebyshev coefficients of R(t) = 2x p'(x) + (2n+1)p(x) - exp(-x).

    Exact algebra on the stored doubles plus the Bessel series of the
    exponential.  R is formally infinite; the returned list holds n_terms and
    the last entry is a rigorous bound on the discarded tail's l1 mass.
    """
    piece = TABLE["pieces"][piece_index_]
    c = coef_list(piece_index_)
    m, h = mapping(piece)
    d = len(c) - 1
    g = poly_deriv_coeffs(c)          # dp/dt
    # P(t) = 2 x(t) p_x'(x) + (2n+1) p(t) = 2 (m + h t) g(t)/h + (2n+1) p(t)
    # = (2m/h) g(t) + 2 t g(t) + (2n+1) p(t)
    P = [mpf(0)] * n_terms

    scale = 2 * m / h
    for j, gj in enumerate(g):
        P[j] += scale * gj
    # t * g(t) in the Chebyshev basis
    if g:
        P[1] += 2 * g[0]
        for j, gj in enumerate(g):
            if j >= 1:
                P[j + 1] += gj
                P[j - 1] += gj
    for j, cj in enumerate(c):
        P[j] += (2 * order + 1) * cj

    # exp(-x(t)) = exp(-m) exp(-h t); exp(-h t) = I0(h) + 2 sum (-1)^k I_k(h) T_k(t)
    em = exp(-m)
    R = list(P)
    R[0] -= em * besseli(0, h)
    for j in range(1, n_terms):
        R[j] -= em * 2 * (-1) ** j * besseli(j, h)

    tail = em * 2 * sum(besseli(j, h) for j in range(n_terms, n_terms + 60))
    return R, tail


def residual_l1(piece, order):
    R, tail = residual_cheb(piece, order)
    return sum(abs(v) for v in R) + tail


def residual_signed_interval(piece, order):
    """Rigorous [min, max] of R over the piece from its Chebyshev coefficients."""
    R, tail = residual_cheb(piece, order)
    s = sum(abs(v) for v in R[1:]) + tail
    return R[0] - s, R[0] + s


def residual_deriv_l1(piece, order):
    R, tail = residual_cheb(piece, order)
    g = poly_deriv_coeffs(R)
    # dR/dt = g(t), so dR/dx = g(t)/h; the l1 mass of g bounds max|dR/dt|.
    return sum(abs(v) for v in g) + 2 * tail


# --------------------------------------------------------------------------
# residue of the consultant's bound: the crude coefficient route
# --------------------------------------------------------------------------

def crude_piece_bounds(piece, order):
    """max|p| <= sum|c|, max|p'| <= sum k^2|c|, max|p''| <= sum k^2(k^2-1)/3 |c|."""
    c = coef_list(piece)
    m, h = mapping(piece)
    s0 = sum(abs(v) for v in c)
    s1 = sum(mpf(k) ** 2 * abs(c[k]) for k in range(len(c)))
    s2 = sum(mpf(k) ** 2 * (mpf(k) ** 2 - 1) / 3 * abs(c[k]) for k in range(len(c)))
    return s0, s1 / h, s2 / (h * h)


def crude_residual_bound(piece, order):
    """The consultant's bound: |R(s)| <= 2s max|p'| + (2n+1) max|p| + exp(-a)."""
    s0, s1, _ = crude_piece_bounds(piece, order)
    p = TABLE["pieces"][piece]
    return 2 * mpf(p["b"]) * s1 + (2 * order + 1) * s0 + exp(-mpf(p["a"]))
