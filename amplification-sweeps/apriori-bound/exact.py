"""The stored fit's exact error, from the anchored ODE identity by quadrature.

This evaluates the error; it does not bound it.  Used to check the signed
intervals and, later, to show that a reference-free ordering is possible as a
computation even though it is not possible as a bound.
"""

from mpmath import mp, mpf, quad

import bounds

mp.dps = 70


def fit_exact_error(order, band_index, x):
    x = mpf(x)
    poly = bounds.fit_poly(order, band_index)
    d = 2 * order + 1

    def integrand(p, s):
        return s ** (order - mpf(1) / 2) * p.residual(s)

    if band_index == 0:
        first = poly
        assert x <= first.b
        return (mpf(1) / 2) * x ** (-(order + mpf(1) / 2)) * quad(
            lambda s: integrand(first, s), [0, x])

    first = bounds.fit_poly(order, 0)
    second = poly
    e0 = (mpf(1) / 2) * first.b ** (-(order + mpf(1) / 2)) * quad(
        lambda s: integrand(first, s), [0, first.b])
    delta = second.value(first.b) - first.value(first.b)
    rho = (first.b / x) ** (order + mpf(1) / 2)
    e1 = (mpf(1) / 2) * x ** (-(order + mpf(1) / 2)) * quad(
        lambda s: integrand(second, s), [first.b, x])
    return rho * (delta + e0) + e1


def exact_Fn_from_fit(order, band_index, x):
    poly = bounds.fit_poly(order, band_index)
    return poly.value(x) - fit_exact_error(order, band_index, x)


_e_a1 = {}


def seed_exact_error(x):
    """Exact E(x) of the extended order-0 fit on [theta0, kX0)."""
    x = mpf(x)
    theta0 = bounds.EXT_LO
    poly = bounds.ext_poly()
    if "v" not in _e_a1:
        _e_a1["v"] = fit_exact_error(0, 0, theta0)
    delta = bounds.ext_delta()
    e_int = (mpf(1) / 2) * x ** (-mpf(1) / 2) * quad(
        lambda s: s ** (-mpf(1) / 2) * poly.residual(s), [theta0, x])
    rho = (theta0 / x) ** (mpf(1) / 2)
    return rho * (delta + _e_a1["v"]) + e_int


def entry_truncation_exact(band, n, method, x):
    """Exact truncation error of the value the entry delivers at x."""
    x = mpf(x)
    bi = 0 if band == "A1" else 1
    if method == "product":
        return fit_exact_error(n, bi, x)
    if x < bounds.THETA0:
        transfer = bounds.amplification(n, x) / bounds.amplification(bounds.NMAX, x)
        return transfer * fit_exact_error(bounds.NMAX, 0, x)
    if n == 0 or x >= mpf(bounds.TABLE["thresholds"][n]):
        return bounds.amplification(n, x) * seed_exact_error(x)
    return fit_exact_error(n, bi, x)
