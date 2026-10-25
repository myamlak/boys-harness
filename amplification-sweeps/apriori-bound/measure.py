"""The a priori bounds at the paper's cells.

Two certified truncation bounds per (order, argument):

  crude   the consultant's coefficient route: max|p| <= sum|c|,
          max|p'| <= sum k^2|c|, |R| <= 2x max|p'| + (2n+1) max|p| + exp(-a).
  sharp   the l1 mass of R's own Chebyshev coefficients, R expanded exactly
          from the stored doubles and the Bessel series of exp(-x).

and two signed intervals:

  coef    [r0 - S, r0 + S] from R's Chebyshev coefficients.
  sampled a dense sample of R per piece with a rigorous slope correction
          max|R'| * h / 2, max|R'| from the l1 mass of R's derivative
          coefficients.

No Boys reference value enters any of them.
"""

import csv
import math
from mpmath import mp, mpf, exp

import core

mp.dps = 60

ROOT = r"C:/Users/myaml/dev/papers-a-v5/papers/boys-crossover/sweeps/amplification_form"


def band_index(band):
    return 0 if band == "A1" else 1


def odd_double_factorial(n):
    """(2n-1)!! as an exact integer."""
    return math.factorial(2 * n) // (2 ** n * math.factorial(n))


def amplification(n, x):
    """A_n(x) = (2n-1)!! / (2x)^n, A_0 = 1."""
    if n == 0:
        return mpf(1)
    return mpf(odd_double_factorial(n)) / (2 * x) ** n


# ---------------------------------------------------------------------------
# per-piece certified quantities (cached by piece index and order)
# ---------------------------------------------------------------------------

_sharp = {}
_coef_interval = {}
_deriv_l1 = {}


def sharp_l1(piece, order):
    key = (piece, order)
    if key not in _sharp:
        _sharp[key] = core.residual_l1(piece, order)
    return _sharp[key]


def coef_signed_interval(piece, order):
    key = (piece, order)
    if key not in _coef_interval:
        _coef_interval[key] = core.residual_signed_interval(piece, order)
    return _coef_interval[key]


def deriv_l1(piece, order):
    key = (piece, order)
    if key not in _deriv_l1:
        _deriv_l1[key] = core.residual_deriv_l1(piece, order)
    return _deriv_l1[key]


def pieces_up_to(order, x):
    """Every piece of this order that [0, x] touches, with its upper limit."""
    out = []
    for index in range(core.TABLE["piece_start"][order], core.TABLE["piece_start"][order + 1]):
        piece = core.TABLE["pieces"][index]
        lo, hi = mpf(piece["a"]), mpf(piece["b"])
        if lo >= x:
            break
        out.append((index, lo, min(hi, x)))
    return out


def residual_at_piece(piece, order, s):
    p = core.TABLE["pieces"][piece]
    t = core.t_of(p, s)
    m, h = core.mapping(p)
    c = core.coef_list(piece)
    # 2s p'(s) + (2n+1) p(s) - exp(-s)
    return 2 * s * core.poly_deriv_value(c, t) / h + (2 * order + 1) * core.poly_value(c, t) - exp(-s)


def sampled_interval(piece, order, lo, hi, samples):
    """[min, max] of R on [lo, hi] within this piece, slope-corrected."""
    if hi <= lo:
        r = residual_at_piece(piece, order, lo)
        return r, r
    n = max(2, int(samples))
    step = (hi - lo) / n
    vals = [residual_at_piece(piece, order, lo + step * j) for j in range(n + 1)]
    correction = deriv_l1(piece, order) * step / 2
    return min(vals) - correction, max(vals) + correction


def truncation_bound(order, x, kind="sharp"):
    """max|R| over [0, x] for this order, and the resulting |E| envelope."""
    if kind == "crude":
        total = mpf(0)
        for index, _, hi in pieces_up_to(order, x):
            total = max(total, core.crude_residual_bound(index, order))
        return total / (2 * order + 1)
    total = mpf(0)
    for index, _, _ in pieces_up_to(order, x):
        total = max(total, sharp_l1(index, order))
    return total / (2 * order + 1)


def signed_bounds(order, x, method="sampled", samples=400):
    """[lo, hi] of E(x) over [0, x] for this order's fit."""
    lo = mpf(0)
    hi = mpf(0)
    first = True
    for index, plo, phi in pieces_up_to(order, x):
        if method == "coef":
            a, b = coef_signed_interval(index, order)
        else:
            a, b = sampled_interval(index, order, plo, phi, samples)
        if first:
            lo, hi, first = a, b, False
        else:
            lo, hi = min(lo, a), max(hi, b)
    d = 2 * order + 1
    return lo / d, hi / d


# ---------------------------------------------------------------------------
# the paper's arithmetic convention (not certified; stated, not proved)
# ---------------------------------------------------------------------------

def read_reference():
    ref = {}
    with open(f"{ROOT}/data/reference.csv", newline="") as handle:
        reader = csv.reader(handle)
        next(reader)
        for row in reader:
            ref[(row[0], float(row[1]))] = [mpf(v) for v in row[2:]]
    return ref


REFERENCE = read_reference()


def recurrence_arithmetic(n, x, values):
    """The paper's charge: 3 * 2^-53 per upward step, carried by A_n / A_l."""
    total = mpf(0)
    for l in range(1, n + 1):
        ratio = amplification(n, x) / amplification(l, x)
        total += 3 * mpf(2) ** -53 * values[l] * ratio
    return total


def read_per_order(path=None):
    path = path or f"{ROOT}/results/per_order.csv"
    rows = []
    with open(path, newline="") as handle:
        for row in csv.DictReader(handle):
            rows.append(dict(band=row["band"], n=int(row["n"]), method=row["method"],
                             worst_abs=mpf(row["worst_abs"]), at_x_abs=float(row["at_x_abs"]),
                             worst_rel=mpf(row["worst_rel"]), at_x_rel=float(row["at_x_rel"])))
    return rows


# ---------------------------------------------------------------------------
# one entry's certified error at (n, x)
# ---------------------------------------------------------------------------

def entry_bounds(band, n, method, x, bound_kind="sharp", signed_method="sampled",
                 include_arithmetic=True, samples=400):
    """Certified [lo, hi] on the truncation error (plus convention) of one entry.

    half is the certified symmetric truncation envelope, before any convention.
    """
    x = mpf(x)
    bi = band_index(band)
    if method == "product":
        lo, hi = signed_bounds(n, x, signed_method, samples)
        half = truncation_bound(n, x, bound_kind)
        charge = mpf(0)
        if include_arithmetic:
            # the paper's direct charge: the stored coefficient scale
            scale = sum(abs(v) for v in core.coef_list(core.piece_of(n, bi)))
            charge = 3 * mpf(2) ** -53 * scale
        return dict(lo=lo - charge, hi=hi + charge, half=half, mode="fit",
                    arithmetic=charge, A=mpf(0))

    if x < core.TABLE["thresholds"][0]:
        return None  # outside the model's population: the entry seeds downward

    if n == 0 or x >= core.TABLE["thresholds"][n]:
        seed_lo, seed_hi = signed_bounds(0, x, signed_method, samples)
        seed_half = truncation_bound(0, x, bound_kind)
        A = amplification(n, x)
        arith = mpf(0)
        if include_arithmetic and n >= 1:
            values = REFERENCE[(band, float(x))]
            arith = recurrence_arithmetic(n, x, values)
        return dict(lo=A * seed_lo - arith, hi=A * seed_hi + arith,
                    half=seed_half * A, mode="recursed", arithmetic=arith, A=A)
    lo, hi = signed_bounds(n, x, signed_method, samples)
    half = truncation_bound(n, x, bound_kind)
    return dict(lo=lo, hi=hi, half=half, mode="same-fit", arithmetic=mpf(0), A=mpf(0))


def magnitude_interval(lo, hi):
    """[min|E|, max|E|] from a signed interval."""
    upper = max(abs(lo), abs(hi))
    lower = mpf(0) if lo <= 0 <= hi else min(abs(lo), abs(hi))
    return lower, upper
