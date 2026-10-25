"""Verification of the consultant's mathematics and the parser, at 60 digits.

The high-precision Boys values here are a *check*, never an input to a bound.
"""

import csv
from mpmath import mp, mpf, exp, quad, inf

import coeffs
import core

mp.dps = 60

ROOT = r"C:/Users/myaml/dev/papers-a-v5/papers/boys-crossover/sweeps/amplification_form"


def boys(n, x):
    """F_n(x) by the defining integral, at working precision (verification only)."""
    return quad(lambda t: t ** (2 * n) * exp(-x * t * t), [0, 1])


def read_cells():
    rows = []
    with open(f"{ROOT}/results/cells.csv", newline="") as handle:
        for row in csv.DictReader(handle):
            rows.append((row["band"], float(row["x"]), int(row["n"]),
                         float(row["recursed"]), float(row["direct"])))
    return rows


def check_polynomial_against_delivered(cells):
    worst_abs = mpf(0)
    worst_rel = mpf(0)
    where = None
    for band, x, n, rec, direct in cells:
        bi = 0 if band == "A1" else 1
        piece = core.piece_of(n, bi)
        c = core.coef_list(piece)
        t = core.t_of(core.TABLE["pieces"][piece], x)
        value = core.poly_value(c, t)
        d = abs(value - mpf(direct))
        rel = d / abs(mpf(direct)) if direct != 0 else d
        if rel > worst_rel:
            worst_rel, where = rel, (band, x, n, direct, float(value))
        worst_abs = max(worst_abs, d)
    print("polynomial vs direct cell: worst abs %.3e, worst rel %.3e at %r"
          % (float(worst_abs), float(worst_rel), where))


def check_derivative(c):
    g = core.poly_deriv_coeffs(c)
    worst = mpf(0)
    for j in range(1, 40):
        t = mpf(j) / 41 - 1
        exact = core.poly_value(g, t)
        h = mpf(10) ** -12
        fd = (core.poly_value(c, t + h) - core.poly_value(c, t - h)) / (2 * h)
        worst = max(worst, abs(exact - fd))
    return worst


def residual_at(order, s):
    idx = core.piece_index(order, s)
    cp = core.coef_list(idx)
    pp = core.TABLE["pieces"][idx]
    ts = core.t_of(pp, s)
    m, h = core.mapping(pp)
    return 2 * s * core.poly_deriv_value(cp, ts) / h + (2 * order + 1) * core.poly_value(cp, ts) - exp(-s)


def check_residual_identity(order, x):
    """E(x) = p(x) - F_n(x) against (1/2) x^{-(n+1/2)} int_0^x s^{n-1/2} R(s) ds."""
    piece = core.piece_index(order, x)
    p = core.TABLE["pieces"][piece]
    t = core.t_of(p, x)
    pc = core.poly_value(core.coef_list(piece), t)
    E = pc - boys(order, x)
    integral = (mpf(1) / 2) * x ** (-(order + mpf(1) / 2)) * quad(
        lambda s: s ** (order - mpf(1) / 2) * residual_at(order, s), [0, x]
    )
    # and the residual identity itself, at the argument
    Rc, _ = core.residual_cheb(piece, order)
    series = core.poly_value(Rc, t)
    direct = residual_at(order, x)
    print("order %2d x %-6s |E - integral| = %.3e   |R_series - R| = %.3e"
          % (order, mp.nstr(x, 6), float(abs(E - integral)), float(abs(series - direct))))


def check_weighted_average(order, x, samples=600):
    """min R/(2n+1) <= E(x) <= max R/(2n+1), on the sampled residual."""
    vals = [residual_at(order, x * mpf(j) / samples) for j in range(samples + 1)]
    lo, hi = min(vals) / (2 * order + 1), max(vals) / (2 * order + 1)
    piece = core.piece_index(order, x)
    pc = core.poly_value(core.coef_list(piece), core.t_of(core.TABLE["pieces"][piece], x))
    E = pc - boys(order, x)
    print("order %2d x %-6s E=%9.3e  sampled [%9.3e, %9.3e]  inside=%s"
          % (order, mp.nstr(x, 6), float(E), float(lo), float(hi), bool(lo <= E <= hi)))


if __name__ == "__main__":
    print("== derivative algebra ==")
    for order, bi in ((0, 0), (5, 1), (16, 0), (32, 1)):
        c = core.coef_list(core.piece_of(order, bi))
        print("order", order, "band", bi, "worst |d - fd| =", float(check_derivative(c)))

    print("== the piece polynomial is the direct entry's value ==")
    check_polynomial_against_delivered(read_cells())

    print("== consultant identity ==")
    for order, x in ((0, mpf("3.5")), (4, mpf("2.0")), (16, mpf("4.9")), (2, mpf("8.5"))):
        check_residual_identity(order, x)

    print("== positive-kernel weighted average ==")
    for order, x in ((0, mpf("3.5")), (4, mpf("2.0")), (16, mpf("4.9")), (2, mpf("8.5"))):
        check_weighted_average(order, x)

    print("== the Boys identity itself: 2x F_n' + (2n+1)F_n = exp(-x) ==")
    for order in (0, 5, 16):
        x = mpf("2.4")
        h = mpf(10) ** -20
        deriv = (boys(order, x + h) - boys(order, x - h)) / (2 * h)
        print("order", order, "residual", mp.nstr(2 * x * deriv + (2 * order + 1) * boys(order, x) - exp(-x), 5))
