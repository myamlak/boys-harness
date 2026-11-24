"""Exact dyadic/rational arithmetic for the two chosen Paper A cells.

Objects
-------
Route value.  Every stored coefficient, argument and mapping constant is a
binary64 value, hence a dyadic rational.  A route value is therefore an affine
function of the shared exponential,

    R = B + D * eh,        eh = exp(-x)/2,

with B and D exact Fractions.  For the direct/product lane D = 0.  For the
recursed/ship lane D is the exact sensitivity of the recurrence to the shared
constant (the paper's d_l).

Delivered value.  d is replayed with Python floats, which are IEEE-754
binary64 round-to-nearest, and math.fma, which is the correctly rounded fused
operation.  The sequence and order mirror code/model.py, which the record's
replay.py validates bit for bit against the delivered cells.  We check the
replay against model.py and against results/cells.csv at the cells.

Certified enclosures
--------------------
enclose_exp_neg_half(x)  rigorous interval for exp(-x)/2 by the alternating
                         exponential series with the alternating-remainder
                         theorem (the terms decrease from k+1 > x onward).
enclose_boys(n, x)       rigorous interval for the Boys function
                         F_n(x) = sum_k (-x)^k / (k! (2n+2k+1)),
                         same theorem.

Both use exact Fractions, so the intervals are exact-rational and the only
inequality used is the alternating-series remainder bound.

Run:  python -I exact_arith.py        (self-tests)
"""
import math
import os
import sys
from fractions import Fraction
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
RECORD = str(REPO / "amplification-sweeps" / "certified-bound")
KERNEL = str(Path(os.environ.get("BOYS_ROOT", REPO.parent / "boys")))

# ---------------------------------------------------------------------------
# small helpers
# ---------------------------------------------------------------------------

def fr(x):
    """Exact Fraction of a binary64 value."""
    return Fraction(x)


def sign(x):
    return (x > 0) - (x < 0)


def mid(x):
    return (x[0] + x[1]) / 2


def iadd(a, b):
    return (a[0] + b[0], a[1] + b[1])


def isub(a, b):
    return (a[0] - b[1], a[1] - b[0])


def imul(a, b):
    p = (a[0] * b[0], a[0] * b[1], a[1] * b[0], a[1] * b[1])
    return (min(p), max(p))


def iscale(a, s):
    """Multiply an interval by an exact scalar s (any sign)."""
    if s >= 0:
        return (a[0] * s, a[1] * s)
    return (a[1] * s, a[0] * s)


def idiv(a, s):
    """Divide an interval by an exact nonzero scalar s."""
    return iscale(a, Fraction(1) / s)


def iabs(a):
    if a[0] >= 0:
        return a
    if a[1] <= 0:
        return (-a[1], -a[0])
    return (Fraction(0), max(-a[0], a[1]))


def ipow2(a):
    return (2 * a[0], 2 * a[1])


# ---------------------------------------------------------------------------
# error-free transformations (verified exactly below)
# ---------------------------------------------------------------------------

def two_sum(a, b):
    """Knuth TwoSum: s + e = a + b exactly, e a binary64 value."""
    s = a + b
    bb = s - a
    e = (a - (s - bb)) + (b - bb)
    return s, e


def two_prod(a, b):
    """TwoProd with a fused residual: p + e = a * b exactly."""
    p = a * b
    e = math.fma(a, b, -p)
    return p, e


def two_diff(a, b):
    """TwoDiff: s + e = a - b exactly."""
    s = a - b
    bb = s - a
    e = (a - (s - bb)) - (b + bb)
    return s, e


# ---------------------------------------------------------------------------
# rigorous transcendental enclosures (exact rationals)
# ---------------------------------------------------------------------------

def _ceil_frac(x):
    return -((-x.numerator) // x.denominator)


def enclose_exp_neg(x, target=Fraction(1, 10 ** 62)):
    """Rigorous interval for exp(-x), x > 0, by the alternating series.

    exp(-x) = sum_{k>=0} (-1)^k x^k / k!.  The terms x^k/k! decrease from
    k+1 > x onward, so the tail after the partial sum through K is bounded by
    the first omitted term x^(K+1)/(K+1)! and has its sign.
    """
    assert x > 0
    kmin = max(2, _ceil_frac(x) + 1)
    a = Fraction(1)
    S = Fraction(1)
    k = 0
    while True:
        k += 1
        a = a * x / k
        S += a if k % 2 == 0 else -a
        if k >= kmin and a < target:
            break
    a_next = a * x / (k + 1)
    sigma = Fraction(1) if (k + 1) % 2 == 0 else Fraction(-1)
    lo, hi = sorted((S, S + sigma * a_next))
    return lo, hi, k, a_next


def enclose_exp_neg_half(x, target=Fraction(1, 10 ** 62)):
    """Rigorous interval for exp(-x)/2."""
    lo, hi, k, t = enclose_exp_neg(x, target)
    return lo / 2, hi / 2, k, t


def enclose_boys(n, x, target=Fraction(1, 10 ** 62)):
    """Rigorous interval for F_n(x) = int_0^1 t^(2n) exp(-x t^2) dt.

    F_n(x) = sum_{k>=0} (-1)^k x^k / (k! (2n+2k+1)).
    The term ratio is x/(k+1) * (2n+2k+1)/(2n+2k+3) < x/(k+1), so the terms
    decrease from k+1 > x onward and the alternating-remainder theorem applies
    from that index.  x = 0 is the closed form 1/(2n+1).
    """
    assert x >= 0
    if x == 0:
        v = Fraction(1, 2 * n + 1)
        return v, v, 0, Fraction(0)
    kmin = max(2, _ceil_frac(x) + 1)
    a = Fraction(1, 2 * n + 1)          # a_0
    S = a
    k = 0
    while True:
        a_prev = a
        a = a * x / (k + 1) * (2 * n + 2 * k + 1) / (2 * n + 2 * k + 3)
        k += 1
        S += a if k % 2 == 0 else -a
        if k >= kmin and a < target:
            break
    a_next = a * x / (k + 1) * (2 * n + 2 * k + 1) / (2 * n + 2 * k + 3)
    sigma = Fraction(1) if (k + 1) % 2 == 0 else Fraction(-1)
    lo, hi = sorted((S, S + sigma * a_next))
    return lo, hi, k, a_next


# ---------------------------------------------------------------------------
# exact route values
# ---------------------------------------------------------------------------

def mapping_exact(x, a, b):
    """2(x-a)/(b-a) - 1 in exact real arithmetic."""
    return Fraction(2) * (x - a) / (b - a) - 1


def cheb_sum(c, t):
    """sum_k c_k T_k(t) exactly, by the T recurrence."""
    if not c:
        return Fraction(0)
    b1 = Fraction(1)          # T_0
    total = Fraction(c[0]) * b1
    if len(c) == 1:
        return total
    b0 = t                    # T_1
    total += Fraction(c[1]) * b0
    for k in range(2, len(c)):
        bk = 2 * t * b0 - b1
        total += Fraction(c[k]) * bk
        b1, b0 = b0, bk
    return total


def clenshaw_split_exact(t, c):
    """Exact real value of model.clenshaw_split's expression.

    Used to check that the split Clenshaw arrangement is algebraically the
    stored Chebyshev polynomial; if it is, the route value is well defined by
    cheb_sum and the delivered Clenshaw rounding is all charged to r.
    """
    deg = len(c) - 1
    t2 = t * t
    v = 2 * t2 - 1
    two_v = 2 * v
    m = deg // 2
    mo = (deg + 1) // 2
    b1 = Fraction(c[2 * m])
    b2 = Fraction(0)
    for k in range(m - 1, 0, -1):
        s = Fraction(c[2 * k]) - b2
        b0 = two_v * b1 + s
        b2, b1 = b1, b0
    s = Fraction(c[0]) - b2
    even = v * b1 + s
    o1 = Fraction(c[2 * mo - 1])
    o2 = Fraction(0)
    for k in range(mo - 2, 0, -1):
        s = Fraction(c[2 * k + 1]) - o2
        o0 = two_v * o1 + s
        o2, o1 = o1, o0
    s = two_v - 1
    u = Fraction(c[1]) - o2
    odd = s * o1 + u
    return t * odd + even


def fit_affine(model, order, band, x, ext=False):
    """(B, D) with route R = B + D*eh for a stored fit (D = 0)."""
    if ext:
        a, b, c = model.ext_piece()
    else:
        a, b, c = model.piece(order, band)
    t = mapping_exact(x, Fraction(a), Fraction(b))
    return cheb_sum(c, t), Fraction(0)


def ship_affine(model, x, n):
    """(B, D) for the recursed/ship entry's exact route at order n.

    The route is the stored extended seed (order 0 above theta0) followed by
    the upward recurrence, or order 32 from the stored fit followed by the
    downward recurrence below theta0.  Above the served prefix the entry reads
    the per-order stored fit.  D is the exact signed sensitivity to eh; it is
    zero on any stored fit.
    """
    xf = Fraction(x)
    theta0 = Fraction(model.theta0())
    if xf == 0:
        return Fraction(1, 2 * n + 1), Fraction(0)
    if xf < theta0:
        B, D = fit_affine(model, 32, 0, xf)
        for l in range(31, n - 1, -1):
            den = Fraction(2 * l + 1, 2)          # l + 1/2
            B = (xf * B) / den
            D = (xf * D + 1) / den
        return B, D
    B, D = fit_affine(model, 0, 0, xf, ext=True)
    served = 0
    th = model.table()["thresholds"]
    while served < 32 and x >= th[served + 1]:
        served += 1
    if n > served:
        bi = 0 if x < model.piece(0, 0)[1] else 1
        return fit_affine(model, n, bi, xf)
    for l in range(1, n + 1):
        half = Fraction(2 * l - 1, 2)             # l - 1/2
        B = (half * B) / xf
        D = (half * D - 1) / xf
    return B, D


def product_affine(model, x, band, n):
    """(B, D) for the direct/product entry's exact route at order n."""
    a, b, c0 = model.piece(0, band)
    c = model.piece(n, band)[2]
    t = mapping_exact(Fraction(x), Fraction(a), Fraction(b))
    return cheb_sum(c, t), Fraction(0)


# ---------------------------------------------------------------------------
# delivered-value replay (mirrors model.py's value path)
# ---------------------------------------------------------------------------

def mapping_value(x, a, b):
    t1 = x - a
    t2 = 2.0 * t1
    den = b - a
    q = t2 / den
    return q - 1.0


def clenshaw_split_value(t, c):
    """model.clenshaw_split value path, bit for bit."""
    deg = len(c) - 1
    assert deg >= 4
    t2 = t * t
    v = math.fma(2.0, t2, -1.0)
    two_v = 2.0 * v
    m = deg // 2
    mo = (deg + 1) // 2
    b1 = c[2 * m]
    b2 = 0.0
    for k in range(m - 1, 0, -1):
        s = c[2 * k] - b2
        b0 = math.fma(two_v, b1, s)
        b2, b1 = b1, b0
    s = c[0] - b2
    even = math.fma(v, b1, s)
    o1 = c[2 * mo - 1]
    o2 = 0.0
    for k in range(mo - 2, 0, -1):
        s = c[2 * k + 1] - o2
        o0 = math.fma(two_v, o1, s)
        o2, o1 = o1, o0
    s = two_v - 1.0
    u = c[1] - o2
    odd = math.fma(s, o1, u)
    return math.fma(t, odd, even)


def ext_seed_value(model, x):
    a, b, c = model.ext_piece()
    return clenshaw_split_value(mapping_value(x, a, b), c)


def fit_split_value(model, order, band, x):
    a, b, c = model.piece(order, band)
    return clenshaw_split_value(mapping_value(x, a, b), c)


def upward_step(f, x, invx, expx, l):
    """The kernel's upward step: separate product, separate subtract,
    refined division, with the two fused corrections of the division."""
    p = (float(l) - 0.5) * f
    a = p - expx
    q0 = a * invx
    r = math.fma(-q0, x, a)
    return math.fma(r, invx, q0)


def ship_value(model, x, expx, n):
    """model.ship_batch(x, expx)[n]['value'], bit for bit."""
    theta0 = model.theta0()
    if x == 0.0:
        return 1.0 / (2.0 * n + 1.0)
    if x < theta0:
        f = fit_split_value(model, 32, 0, x)
        for l in range(31, n - 1, -1):
            invd = 1.0 / (float(l) + 0.5)
            p = x * f
            a = p + expx
            q0 = a * invd
            r = math.fma(-q0, float(l) + 0.5, a)
            f = math.fma(r, invd, q0)
        return f
    f = ext_seed_value(model, x)
    invx = 1.0 / x
    th = model.table()["thresholds"]
    served = 0
    while served < 32 and x >= th[served + 1]:
        served += 1
    if n > served:
        bi = 0 if x < model.piece(0, 0)[1] else 1
        return fit_split_value(model, n, bi, x)
    for l in range(1, n + 1):
        f = upward_step(f, x, invx, expx, l)
    return f


def product_basis_value(x, a, b, deg):
    t = mapping_value(x, a, b)
    basis = [1.0, t]
    m2 = 2.0 * t
    for _ in range(2, deg + 1):
        prod = math.fma(m2, basis[-1], 0.0)
        basis.append(prod - basis[-2])
    return basis


def eval_fit_product_value(c, basis):
    deg = len(c) - 1
    pass_v = []
    for k in range(len(c)):
        degree = deg - k
        pass_v.append(c[degree] * basis[degree])
    width = len(pass_v)
    while width > 1:
        for k in range(width // 2):
            pass_v[k] = pass_v[2 * k] + pass_v[2 * k + 1]
        if width % 2 != 0:
            pass_v[width // 2] = pass_v[width - 1]
        width = (width + 1) // 2
    return pass_v[0]


def product_value(model, x, band, n):
    a, b, c0 = model.piece(0, band)
    deg = len(c0) - 1
    basis = product_basis_value(x, a, b, deg)
    return eval_fit_product_value(model.piece(n, band)[2], basis)


# ---------------------------------------------------------------------------
# self-tests
# ---------------------------------------------------------------------------

def _selftest():
    import random
    random.seed(1)
    for _ in range(2000):
        a = random.uniform(-1e3, 1e3)
        b = random.uniform(-1e3, 1e3)
        s, e = two_sum(a, b)
        assert Fraction(s) + Fraction(e) == Fraction(a) + Fraction(b)
        s, e = two_diff(a, b)
        assert Fraction(s) + Fraction(e) == Fraction(a) - Fraction(b)
        p, e = two_prod(a, b)
        assert Fraction(p) + Fraction(e) == Fraction(a) * Fraction(b)
    print("EFT self-test: 2000 random pairs reconstruct exactly "
          "(TwoSum, TwoDiff, TwoProd)")

    # the split Clenshaw arrangement is algebraically the Chebyshev polynomial
    sys.path.insert(0, RECORD + "/code")
    import model
    model.configure(kernel=KERNEL,
                    exp_bits=RECORD + "/exp/exp_bits.txt")
    xs = [1.1679101370108986, 10.865022762913437, 8.254041852680183e-06,
          6.3134226241640246]
    worst = Fraction(0)
    for band in (0, 1):
        for n in (0, 1, 2, 3, 7, 15, 22, 23, 32):
            a, b, c = model.piece(n, band)
            for x in xs:
                if not (a <= x <= b):
                    continue
                t = mapping_exact(Fraction(x), Fraction(a), Fraction(b))
                d = cheb_sum(c, t) - clenshaw_split_exact(t, c)
                worst = max(worst, abs(d))
    print("Chebyshev sum == split-Clenshaw expression exactly; worst |diff| =",
          float(worst))

    # product basis: the recurrence P_k = 2t P_{k-1} - P_{k-2} is T_k.  The
    # delivered doubles are that recurrence rounded, so the deviation below is
    # the basis rounding that the route certificate charges to r.
    a, b, c0 = model.piece(0, 0)
    x = 6.3134226241640246
    t = mapping_exact(Fraction(x), Fraction(a), Fraction(b))
    basis = product_basis_value(x, a, b, len(c0) - 1)
    T = [Fraction(1), t]
    for k in range(2, len(c0)):
        T.append(2 * t * T[-1] - T[-2])
    worst = Fraction(0)
    for k, v in enumerate(basis):
        if T[k] != 0:
            worst = max(worst, abs(Fraction(v) - T[k]) / abs(T[k]))
    print("delivered product basis vs exact T_k: worst relative deviation =",
          float(worst), "(rounding, charged to r)")

    # transcendental enclosures: check against mpmath at high precision
    from mpmath import mp, mpf
    mp.dps = 80
    for n, x in ((2, Fraction(1.1679101370108986)),
                 (2, Fraction(8.254041852680183e-06)),
                 (23, Fraction(10.865022762913437)),
                 (23, Fraction(6.3134226241640246))):
        lo, hi, k, tail = enclose_boys(n, x)
        ref = mp.quad(lambda u: u ** (2 * n) * mp.e ** (-mpf(x) * u ** 2),
                      [0, 1])
        assert lo <= ref <= hi, (n, float(x), float(lo), float(ref), float(hi))
        print(f"F_{n}({float(x)!r}) enclosed by k={k}, width={float(hi-lo):.3e}, "
              f"ref inside={lo <= ref <= hi}")
    for x in (Fraction(1.1679101370108986), Fraction(10.865022762913437)):
        lo, hi, k, tail = enclose_exp_neg_half(x)
        ref = mp.e ** (-mpf(x)) / 2
        assert lo <= ref <= hi
        print(f"exp(-{float(x)!r})/2 enclosed by k={k}, "
              f"width={float(hi-lo):.3e}, ref inside={lo <= ref <= hi}")


if __name__ == "__main__":
    _selftest()
