"""Certified a priori error bounds for the region-A Boys fits, from the stored
double coefficients and exp(-x) alone.

The kernel's region A is piecewise:

    order n, band 0: piece [0, b1]          (degree 20)
    order n, band 1: piece [b1, kX0)        (degree 18)
    the recursed entry's seed: one extended F0 fit on [theta0, kX0)  (degree 24)

The error E = p - F_n of a single polynomial over an interval that starts at
its own left end is bounded by the consultant's identity

    E(x) = (1/2) x^(-(n+1/2)) int_0^x s^(n-1/2) R(s) ds,
    R = 2x p'(x) + (2n+1) p(x) - exp(-x),

so |E| <= max|R|/(2n+1).  Applied naively across a seam, or to a fit read
outside its own interval (the extended seed at arguments below theta0, the
A2 piece below b1), the identity is still exact but the residual there is
O(1) and the bound is useless.  The correct form anchors at the left end of
the polynomial's own interval:

    E(x) = rho^(n+1/2) E(a) + (1/2) x^(-(n+1/2)) int_a^x s^(n-1/2) R(s) ds,
    rho = a/x,

and the anchor E(a) is either zero (first piece of the order, a = 0) or the
jump between the two stored pieces plus the first piece's own certified
error.  Every bound below is of this anchored form, which is what makes it
both valid and tight.
"""

import csv
from mpmath import mp, mpf, exp, besseli

import coeffs
import core

mp.dps = 60

TABLE = coeffs.load()
U = mpf(2) ** -53
THETA0 = mpf(TABLE["thresholds"][0])
X0 = mpf(TABLE["ext_hi"])
EXT_LO = mpf(TABLE["ext_lo"])
EXT_C = [mpf(v) for v in TABLE["ext_coeffs"]]
NMAX = 32


# ---------------------------------------------------------------------------
# polynomials and their residuals
# ---------------------------------------------------------------------------

def mapping_of(a, b):
    a, b = mpf(a), mpf(b)
    return (a + b) / 2, (b - a) / 2


def piece_of(order, band_index):
    return TABLE["piece_start"][order] + band_index


def piece_bounds(index):
    p = TABLE["pieces"][index]
    return mpf(p["a"]), mpf(p["b"])


def coefs(index):
    p = TABLE["pieces"][index]
    return [mpf(v) for v in TABLE["coeffs"][p["off"]: p["off"] + p["deg"] + 1]]


def t_of(a, b, x):
    m, h = mapping_of(a, b)
    return (mpf(x) - m) / h


def p_value(c, t):
    return core.poly_value(c, t)


def p_deriv(c, t):
    return core.poly_deriv_value(c, t)


def residual_cheb(c, a, b, order, n_terms=0):
    """Chebyshev coefficients of R(t) and a bound on the discarded tail's l1
    mass and of the mass of dR/dt.

    R(t) = 2 x(t) p'(x) + (2n+1) p(t) - exp(-x(t)),  x = m + h t.
    The polynomial part is exact algebra on the stored doubles; the
    exponential is the Bessel series  exp(-m)(I0(h) + 2 sum (-1)^k I_k(h) T_k).
    """
    m, h = mapping_of(a, b)
    d = len(c) - 1
    g = core.poly_deriv_coeffs(c)
    need = int(abs(h)) + 80
    n_terms = max(n_terms, need)
    P = [mpf(0)] * n_terms

    scale = 2 * m / h
    for j, gj in enumerate(g):
        P[j] += scale * gj
    if g:
        P[1] += 2 * g[0]
        for j, gj in enumerate(g):
            if j >= 1:
                P[j + 1] += gj
                P[j - 1] += gj
    for j, cj in enumerate(c):
        P[j] += (2 * order + 1) * cj

    em = exp(-m)
    R = list(P)
    R[0] -= em * besseli(0, h)
    for j in range(1, n_terms):
        R[j] -= em * 2 * (-1) ** j * besseli(j, h)

    tail = mpf(0)
    tail_d = mpf(0)
    for j in range(n_terms, n_terms + 80):
        ij = besseli(j, h)
        tail += 2 * em * ij
        tail_d += 2 * em * ij * mpf(j) ** 2
    return R, tail, tail_d


class Poly:
    """One stored Chebyshev polynomial with its certified residual data."""

    def __init__(self, name, c, a, b, order):
        self.name = name
        self.c = c
        self.a, self.b = mpf(a), mpf(b)
        self.order = order
        self.m, self.h = mapping_of(a, b)
        self.R, self.tail, self.tail_d = residual_cheb(c, a, b, order)
        self.l1 = sum(abs(v) for v in self.R) + self.tail
        s = sum(abs(v) for v in self.R[1:]) + self.tail
        self.coef_interval = (self.R[0] - s, self.R[0] + s)
        g = core.poly_deriv_coeffs(self.R)
        self.deriv_l1 = sum(abs(v) for v in g) + self.tail_d
        self.coef_l1 = sum(abs(v) for v in c)
        s0 = sum(abs(v) for v in c)
        s1 = sum(mpf(k) ** 2 * abs(c[k]) for k in range(len(c)))
        self.max_p = s0
        self.max_dp = s1 / abs(self.h)
        # the consultant's coefficient route: |R| <= 2 x max|p'| + (2n+1) max|p| + 1
        self.crude = 2 * self.b * self.max_dp + (2 * order + 1) * s0 + 1

    def t(self, x):
        return (mpf(x) - self.m) / self.h

    def value(self, x):
        return p_value(self.c, self.t(x))

    def deriv_value(self, x):
        return p_deriv(self.c, self.t(x)) / self.h

    def residual(self, x):
        s = mpf(x)
        return 2 * s * self.deriv_value(s) + (2 * self.order + 1) * self.value(s) - exp(-s)

    def sampled_interval(self, x_lo, x_hi, samples):
        """[min, max] of R on [x_lo, x_hi] within this polynomial, rigorous."""
        t_lo, t_hi = self.t(x_lo), self.t(x_hi)
        if t_hi < t_lo:
            t_lo, t_hi = t_hi, t_lo
        if t_hi == t_lo:
            r = p_value(self.R, t_lo)
            return r - self.tail, r + self.tail
        n = max(2, int(samples))
        step = (t_hi - t_lo) / n
        vals = [p_value(self.R, t_lo + step * j) for j in range(n + 1)]
        correction = self.deriv_l1 * step / 2 + self.tail
        return min(vals) - correction, max(vals) + correction

    def interval(self, x_lo, x_hi, method, samples):
        """Rigorous [min, max] of R on [x_lo, x_hi] by the chosen route."""
        if method == "coef":
            return self.coef_interval
        if method == "crude":
            return -self.crude, self.crude
        return self.sampled_interval(x_lo, x_hi, samples)


_polys = {}


def fit_poly(order, band_index):
    key = ("fit", order, band_index)
    if key not in _polys:
        index = piece_of(order, band_index)
        a, b = piece_bounds(index)
        _polys[key] = Poly("F%d[%d]" % (order, band_index), coefs(index), a, b, order)
    return _polys[key]


def ext_poly():
    if "ext" not in _polys:
        _polys["ext"] = Poly("F0ext", EXT_C, EXT_LO, X0, 0)
    return _polys["ext"]


# ---------------------------------------------------------------------------
# the anchored error interval of one stored fit
# ---------------------------------------------------------------------------

def fit_error(order, band_index, x, method="sampled", samples=400):
    """[lo, hi] on E(x) = p_n(x) - F_n(x) for the stored fit of this band.

    Band 0 reads the first piece from its own left end; band 1 anchors at the
    seam, where the two stored pieces are two polynomials for the same F_n.
    """
    x = mpf(x)
    poly = fit_poly(order, band_index)
    d = 2 * order + 1
    if band_index == 0:
        lo, hi = poly.interval(0, x, method, samples)
        return lo / d, hi / d
    first = fit_poly(order, 0)
    lo0, hi0 = first.interval(0, first.b, method, samples)
    delta = poly.value(first.b) - first.value(first.b)
    rho = (first.b / x) ** (order + mpf(1) / 2)
    lo1, hi1 = poly.interval(first.b, x, method, samples)
    lo = rho * (delta + lo0 / d) + (1 - rho) * lo1 / d
    hi = rho * (delta + hi0 / d) + (1 - rho) * hi1 / d
    return lo, hi


def fit_envelope(order, band_index, x, method="sampled", samples=400):
    """Symmetric certified envelope max|E(x)| for the stored fit."""
    lo, hi = fit_error(order, band_index, x, method, samples)
    return max(abs(lo), abs(hi))


def fit_sym(order, band_index, x, kind="sharp"):
    """The consultant's symmetric envelope, max|R| over the pieces read.

    kind='sharp' uses the l1 mass of R's own Chebyshev coefficients;
    kind='crude' the coefficient route max|p|, max|p'|.  Both are anchored at
    the seam exactly as the signed interval is.
    """
    x = mpf(x)
    d = 2 * order + 1
    poly = fit_poly(order, band_index)
    sym = (poly.l1 if kind == "sharp" else poly.crude) / d
    if band_index == 0:
        return sym
    first = fit_poly(order, 0)
    delta = abs(poly.value(first.b) - first.value(first.b))
    rho = (first.b / x) ** (order + mpf(1) / 2)
    sym0 = (first.l1 if kind == "sharp" else first.crude) / d
    return rho * (delta + sym0) + (1 - rho) * sym


_delta_cache = {}


def ext_delta():
    """s(theta0) - p_A1(theta0) on the stored doubles (exact algebra)."""
    if "v" not in _delta_cache:
        _delta_cache["v"] = ext_poly().value(EXT_LO) - fit_poly(0, 0).value(EXT_LO)
    return _delta_cache["v"]


def seed_sym(x, kind="sharp"):
    ex = ext_poly()
    rho = (EXT_LO / mpf(x)) ** mpf("0.5")
    sym_edge = fit_sym(0, 0, EXT_LO, kind)
    sym_ext = ex.l1 if kind == "sharp" else ex.crude
    return rho * (abs(ext_delta()) + sym_edge) + (1 - rho) * sym_ext


# ---------------------------------------------------------------------------
# the recursed entry's seed: the extended F0 fit on [theta0, kX0)
# ---------------------------------------------------------------------------

def ext_anchor(method="sampled", samples=400):
    """Certified [lo, hi] on E_ext(theta0) = s(theta0) - F0(theta0).

    E_ext(theta0) = [s(theta0) - p_A1(theta0)] + [p_A1(theta0) - F0(theta0)],
    the first bracket exact on the stored doubles, the second the order-0
    A1 fit's own certified error at theta0 (its piece covers [0, b1]).
    """
    s = ext_poly()
    p = fit_poly(0, 0)
    delta = s.value(EXT_LO) - p.value(EXT_LO)
    lo, hi = fit_error(0, 0, EXT_LO, method, samples)
    return delta + lo, delta + hi, delta


def seed_error(x, method="sampled", samples=400):
    """[lo, hi] on E_ext(x) for x in [theta0, kX0)."""
    x = mpf(x)
    ex = ext_poly()
    lo_a, hi_a, delta = ext_anchor(method, samples)
    rho = (EXT_LO / x) ** mpf("0.5")
    lo_r, hi_r = ex.interval(EXT_LO, x, method, samples)
    lo = rho * lo_a + (1 - rho) * lo_r
    hi = rho * hi_a + (1 - rho) * hi_r
    return lo, hi, delta


def seed_envelope(x, method="sampled", samples=400):
    lo, hi, _ = seed_error(x, method, samples)
    return max(abs(lo), abs(hi))


# ---------------------------------------------------------------------------
# amplification and the paper's accumulated-arithmetic convention
# ---------------------------------------------------------------------------

_odf = {}


def odd_double_factorial(n):
    import math
    if n not in _odf:
        _odf[n] = math.factorial(2 * n) // (2 ** n * math.factorial(n))
    return _odf[n]


def amplification(n, x):
    if n == 0:
        return mpf(1)
    return mpf(odd_double_factorial(n)) / (2 * mpf(x)) ** n


def amplification_ratio(n, l, x):
    if n == l:
        return mpf(1)
    return amplification(n, x) / amplification(l, x)


def _step_size(order, x, size_kind):
    """A size for the value the step produces."""
    if size_kind == "bound":
        return mpf(1) / (2 * order + 1)
    if size_kind == "apriori":
        # F_l = p_l - E_l, so F_l <= |p_l| + the fit's certified envelope
        bi = 0 if x < fit_poly(order, 0).b else 1
        poly = fit_poly(order, bi)
        return abs(poly.value(x)) + fit_sym(order, bi, x, "sharp")
    if size_kind == "fit":
        bi = 0 if x < fit_poly(order, 0).b else 1
        return abs(fit_poly(order, bi).value(x))
    raise ValueError(size_kind)


def step_charge(n, x, size_kind="apriori"):
    """Three roundings at 2^-53 per upward step, carried to order n.

    The paper's convention charges 3u F_l(x) A_n/A_l with the reference F_l.
    A priori, F_l is replaced by the certified size from the stored fit
    (size_kind='apriori') or by the rigorous 1/(2l+1) (size_kind='bound').
    """
    x = mpf(x)
    total = mpf(0)
    for l in range(1, n + 1):
        total += 3 * U * _step_size(l, x, size_kind) * amplification_ratio(n, l, x)
    return total


def downward_charge(n, x, size_kind="apriori"):
    x = mpf(x)
    total = mpf(0)
    for l in range(n, NMAX):
        total += 3 * U * _step_size(l, x, size_kind) * amplification_ratio(n, l, x)
    return total


def charge_uniform(n, band, method, kind="paper"):
    """A single upper bound on the arithmetic charge for the entry over its
    whole population, so a band-maximum envelope stays an upper bound."""
    if method == "product" or kind != "paper":
        poly = fit_poly(n, 0 if band == "A1" else 1)
        return clenshaw_charge(poly) if kind != "paper" else 3 * U * poly.coef_l1
    ch_same = 3 * U * fit_poly(n, 0 if band == "A1" else 1).coef_l1
    # recursed: the step sizes shrink with x, so the charge is largest at theta_n
    x = mpf(TABLE["thresholds"][n])
    ch_rec = 3 * U * ext_poly().coef_l1
    for l in range(1, n + 1):
        ch_rec += 3 * U * _step_size(l, x, "apriori") * amplification_ratio(n, l, x)
    # downward: for x < theta0 and l > n, (2x)^(l-n) grows with x
    x = THETA0
    ch_dn = 3 * U * fit_poly(NMAX, 0).coef_l1
    for l in range(n, NMAX):
        ch_dn += 3 * U * _step_size(l, x, "apriori") * amplification_ratio(n, l, x)
    return max(ch_same, ch_rec, ch_dn)


# ---------------------------------------------------------------------------
# one entry at one argument
# ---------------------------------------------------------------------------

def clenshaw_charge(poly):
    """A worst-case split-Clenshaw summation charge: ~3 roundings per
    even/odd recurrence step against the coefficient l1 mass."""
    d = len(poly.c) - 1
    return (3 * (d // 2) + 3) * U * poly.coef_l1


def product_charge(poly):
    """A pairwise-summation charge: log2(d+1) combination levels plus the
    product and the final accumulation."""
    import math
    d = len(poly.c) - 1
    return (2 + math.ceil(math.log2(d + 1))) * U * poly.coef_l1


def entry_bounds(band, n, method, x, signed_method="sampled", samples=400,
                 size_kind="apriori", charge_kind="paper"):
    """Certified [lo, hi] on the entry's absolute error, with its mode.

    lo/hi come from the signed residual interval (signed_method) plus the
    paper's accumulated-arithmetic convention.  sym_sharp and sym_crude are
    the truncation-only symmetric envelopes of the two coefficient routes.
    """
    x = mpf(x)
    bi = 0 if band == "A1" else 1

    if method == "product":
        poly = fit_poly(n, bi)
        lo, hi = fit_error(n, bi, x, signed_method, samples)
        charge = 3 * U * poly.coef_l1 if charge_kind == "paper" else product_charge(poly)
        return dict(lo=lo - charge, hi=hi + charge, mode="fit", arith=charge,
                    trunc=max(abs(lo), abs(hi)),
                    sym_sharp=fit_sym(n, bi, x, "sharp"),
                    sym_crude=fit_sym(n, bi, x, "crude"),
                    A=mpf(1), charge=charge)

    # the recursed entry
    if x < THETA0:
        # the kernel seeds the downward recursion at order 32 and divides down
        poly32 = fit_poly(NMAX, 0)
        lo32, hi32 = fit_error(NMAX, 0, x, signed_method, samples)
        transfer = amplification(n, x) / amplification(NMAX, x)
        if charge_kind == "paper":
            charge = 3 * U * poly32.coef_l1
        else:
            charge = clenshaw_charge(poly32)
        charge += downward_charge(n, x, size_kind)
        return dict(lo=transfer * lo32 - charge, hi=transfer * hi32 + charge,
                    mode="downward", arith=charge,
                    trunc=transfer * max(abs(lo32), abs(hi32)),
                    sym_sharp=transfer * fit_sym(NMAX, 0, x, "sharp"),
                    sym_crude=transfer * fit_sym(NMAX, 0, x, "crude"),
                    A=transfer, charge=charge)

    if n == 0 or x >= mpf(TABLE["thresholds"][n]):
        lo_s, hi_s, _ = seed_error(x, signed_method, samples)
        A = amplification(n, x)
        if charge_kind == "paper":
            charge = 3 * U * ext_poly().coef_l1
        else:
            charge = clenshaw_charge(ext_poly())
        charge += step_charge(n, x, size_kind)
        return dict(lo=A * lo_s - charge, hi=A * hi_s + charge, mode="recursed",
                    arith=charge, trunc=A * max(abs(lo_s), abs(hi_s)),
                    sym_sharp=A * seed_sym(x, "sharp"),
                    sym_crude=A * seed_sym(x, "crude"),
                    A=A, charge=charge)

    # inside the population but below this order's threshold: the per-order fit
    poly = fit_poly(n, bi)
    lo, hi = fit_error(n, bi, x, signed_method, samples)
    charge = 3 * U * poly.coef_l1 if charge_kind == "paper" else clenshaw_charge(poly)
    return dict(lo=lo - charge, hi=hi + charge, mode="same-fit", arith=charge,
                trunc=max(abs(lo), abs(hi)),
                sym_sharp=fit_sym(n, bi, x, "sharp"),
                sym_crude=fit_sym(n, bi, x, "crude"),
                A=mpf(1), charge=charge)


def entry_envelope(band, n, method, x, kind="sharp", charge_kind="paper"):
    """(truncation envelope, charge) without sampling: cheap, for band maxima.

    The charge is the pointwise paper-convention value for this cell.
    """
    x = mpf(x)
    bi = 0 if band == "A1" else 1
    if method == "product":
        poly = fit_poly(n, bi)
        charge = 3 * U * poly.coef_l1 if charge_kind == "paper" else clenshaw_charge(poly)
        return fit_sym(n, bi, x, kind), charge
    if x < THETA0:
        transfer = amplification(n, x) / amplification(NMAX, x)
        charge = 3 * U * fit_poly(NMAX, 0).coef_l1
        charge += downward_charge(n, x, "apriori")
        return transfer * fit_sym(NMAX, 0, x, kind), charge
    if n == 0 or x >= mpf(TABLE["thresholds"][n]):
        A = amplification(n, x)
        charge = 3 * U * ext_poly().coef_l1 + step_charge(n, x, "apriori")
        return A * seed_sym(x, kind), charge
    poly = fit_poly(n, bi)
    charge = 3 * U * poly.coef_l1 if charge_kind == "paper" else clenshaw_charge(poly)
    return fit_sym(n, bi, x, kind), charge


def magnitude_interval(lo, hi):
    upper = max(abs(lo), abs(hi))
    lower = mpf(0) if lo <= 0 <= hi else min(abs(lo), abs(hi))
    return lower, upper


# ---------------------------------------------------------------------------
# relative error: divide by a certified interval on F_n itself
# ---------------------------------------------------------------------------

def f_n_interval(band, n, x, bound_kind="sharp", signed_method="sampled", samples=400):
    """Certified positive [lo, hi] on F_n(x), from the stored fit of the band."""
    bi = 0 if band == "A1" else 1
    poly = fit_poly(n, bi)
    value = poly.value(x)
    lo, hi = fit_error(n, bi, x, signed_method, samples)
    # F_n = p - E, so F in [p - hi, p - lo]; F_n > 0
    return max(mpf(0), value - hi), value - lo


def relative_interval(band, n, x, entry, samples=400, signed_method="sampled"):
    """Certified [lo, hi] on |value - F_n| / F_n at this cell."""
    err_lo, err_hi = magnitude_interval(entry["lo"], entry["hi"])
    f_lo, f_hi = f_n_interval(band, n, x, signed_method=signed_method, samples=samples)
    if f_lo <= 0:
        return mpf(0), mp.inf
    return err_lo / f_hi, err_hi / f_lo


# ---------------------------------------------------------------------------
# the 132 comparisons
# ---------------------------------------------------------------------------

def certify(lo_s, hi_s, lo_p, hi_p):
    ms = magnitude_interval(lo_s, hi_s)
    mp_ = magnitude_interval(lo_p, hi_p)
    if ms[1] < mp_[0]:
        return "ship"
    if mp_[1] < ms[0]:
        return "product"
    return ""


def certify_relative(band, n, x_s, x_p, ship, product, samples=400,
                     signed_method="sampled"):
    ls, hs = relative_interval(band, n, x_s, ship, samples, signed_method)
    lp, hp = relative_interval(band, n, x_p, product, samples, signed_method)
    if hs < lp:
        return "ship"
    if hp < ls:
        return "product"
    return ""
