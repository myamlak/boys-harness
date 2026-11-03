"""Exact-real value of each entry's evaluated route, at 70-digit precision.

The route's real expression is the stored polynomial (or the recurrence
started from it) at the exact mapping of the exact argument, with the exact
exponential e^-x/2 as the recurrence's shared constant.  Its difference from
the Boys function is the entry's truncation error; the delivered double
differs from it by the route's arithmetic, which model.py bounds.

This evaluation is a computation, not a bound; the record states its error
bound (70-digit working precision on 60-digit references, so <= 1e-45 in the
record's tier T) and it is cross-checked against the deposit's identity-based
oracle, which agrees to 4e-56 on the checked cells.
"""
import csv
import sys

from mpmath import mp, mpf

sys.path.insert(0, r"C:/Users/myaml/dev/_a-cert-run/code")
import model

mp.dps = 70

THETA0 = model.THETA0
B1 = model.piece(0, 0)[1]
_EMPTY = object()

_pieces = {}
_ext = None


def _piece(order, band_index):
    global _ext
    key = (order, band_index)
    p = _pieces.get(key)
    if p is None:
        a, b, c = model.piece(order, band_index)
        p = (mpf(a), mpf(b), [mpf(v) for v in c])
        _pieces[key] = p
    return p


def _ext_piece():
    global _ext
    if _ext is None:
        a, b, c = model.ext_piece()
        _ext = (mpf(a), mpf(b), [mpf(v) for v in c])
    return _ext


def _cheb(coeffs, t):
    """Sum c_k T_k(t) by the T recurrence, in the exact real arithmetic."""
    b1 = mpf(1)
    b0 = t
    total = coeffs[0] * b1 + coeffs[1] * b0
    for k in range(2, len(coeffs)):
        bk = 2 * t * b0 - b1
        total += coeffs[k] * bk
        b1, b0 = b0, bk
    return total


def _map(x, a, b):
    return 2 * (x - a) / (b - a) - 1


def fit_value(order, band_index, x):
    a, b, c = _piece(order, band_index)
    return _cheb(c, _map(mpf(x), a, b))


def ext_value(x):
    a, b, c = _ext_piece()
    return _cheb(c, _map(mpf(x), a, b))


def served_count(x):
    served = 0
    while served < 32 and x >= model.TABLE["thresholds"][served + 1]:
        served += 1
    return served


def product_values(x):
    """Exact-real values of the direct entry's 33 orders at x."""
    xm = mpf(x)
    bi = 0 if xm < B1 else 1
    a, b, c0 = _piece(0, bi)
    t = _map(xm, a, b)
    return [_cheb(_piece(n, bi)[2], t) for n in range(33)]


def ship_values(x, exp_half):
    """Exact-real values of the recursed entry's 33 orders at x.

    exp_half is the route's shared constant; passing e^-x/2 gives the route
    whose error from F_n is the truncation.  The recurrence is evaluated in
    the exact real arithmetic, so the only input errors it carries are the
    polynomial's own (none: the coefficients are exact doubles).
    """
    xm = mpf(x)
    if xm == 0:
        return [mpf(1) / (2 * n + 1) for n in range(33)]
    eh = mpf(exp_half)
    if xm < THETA0:
        f = [None] * 33
        f[32] = fit_value(32, 0, xm)
        for l in range(31, -1, -1):
            f[l] = (xm * f[l + 1] + eh) / (l + mpf("0.5"))
        return f
    served = served_count(xm)
    bi = 0 if xm < B1 else 1
    out = [None] * 33
    f = ext_value(xm)
    out[0] = f
    for l in range(1, served + 1):
        f = ((l - mpf("0.5")) * f - eh) / xm
        out[l] = f
    for n in range(served + 1, 33):
        out[n] = fit_value(n, bi, xm)
    return out


def load_reference(root):
    ref = {}
    with open(root + "/data/reference.csv", newline="") as handle:
        for row in csv.DictReader(handle):
            ref[(row["band"], float(row["x"]))] = [
                mpf(row[c]) for c in list(row)[2:]]
    return ref


if __name__ == "__main__":
    import exact
    root = r"C:/Users/myaml/dev/papers-a-v5/papers/boys-crossover/sweeps/amplification_form"
    ref = load_reference(root)
    cases = [("A1", 1.2141318994825876, 5, "product"),
             ("A2", 5.953583371143639, 5, "product"),
             ("A2", 11.0, 5, "ship"), ("A2", 11.0, 20, "ship")]
    for band, x, n, method in cases:
        if (band, x) not in ref:
            print("skip", band, x)
            continue
        if method == "product":
            val = product_values(x)[n]
        else:
            val = ship_values(x, mp.exp(-mpf(x)) / 2)[n]
        tr = val - ref[(band, x)][n]
        orc = exact.entry_truncation_exact(band, n, method, x)
        print(band, x, n, method, "route", mp.nstr(tr, 12),
              "oracle", mp.nstr(orc, 12), "diff", mp.nstr(tr - orc, 6))
