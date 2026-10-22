"""Reference implementation for the crlibm claim audit.

Reference = MPFR (through gmpy2), used with a Ziv loop:

    for increasing precision p:
        lo = f(x) evaluated at precision p with rounding toward -inf
        hi = f(x) evaluated at precision p with rounding toward +inf
        dl = round_double_down(lo), du = round_double_up(hi)
        if dl == du: the exact value is that double exactly, and all four
                     IEEE rounding modes return it -> decided

MPFR guarantees correct rounding of its own functions, so lo <= exact <= hi
holds at every precision. Directed rounding here means the value of the
rounding-mode function of the *real* argument x, which is what the crlibm
documentation claims for the corresponding *_rd/_ru/_rz entry point.

Nothing here uses crlibm.
"""

import gmpy2
from gmpy2 import mpfr, mpz

RNDN, RNDD, RNDU, RNDZ = (gmpy2.RoundToNearest, gmpy2.RoundDown,
                          gmpy2.RoundUp, gmpy2.RoundToZero)

CONVENTION = "__convention__"        # sentinel: value depends on convention

INF = float("inf")
DBL_MAX_VAL = float.fromhex("0x1.fffffffffffffp+1023")

# ---------------------------------------------------------------- conversion


def _ctx(prec, rnd):
    """MPFR context with the exponent range wide open (exact evaluation)."""
    return gmpy2.context(precision=prec, round=rnd,
                         emax=1073741823, emin=-1073741823)


def _ctx53(rnd):
    """MPFR context that *is* the IEEE-754 double format.

    MPFR normalizes to 0.5 <= m < 1, so a double with exponent 1023
    (1.f * 2**1023) has MPFR exponent 1024: emax = 1024, and the smallest
    subnormal 2**-1074 is 2**(emin-1), so emin = -1073.
    """
    return gmpy2.context(precision=53, round=rnd, emax=1024, emin=-1073,
                         subnormalize=True, trap_overflow=False)


DBL_MAX_V = mpfr(2) ** 1024 - mpfr(2) ** 971      # largest finite double
OVF_MID_V = mpfr(2) ** 1024 - mpfr(2) ** 970      # RNDN overflow threshold


def to_double(v, rnd):
    """Correctly round an mpfr value to an IEEE double under `rnd`.

    Outside the overflow region an MPFR context whose precision and exponent
    range are the double format's gives the IEEE rounding directly.  Inside it
    the overflow rule is applied explicitly: MPFR's own threshold is not the
    IEEE one, and this is exactly the region where cosh/sinh/exp vectors live.
    """
    if gmpy2.is_nan(v):
        return float("nan")
    # Compatibility note: gmpy2 2.3.1 exposes is_infinite, not is_inf.
    if gmpy2.is_infinite(v):
        return -INF if v < 0 else INF
    a = abs(v)
    if a > DBL_MAX_V:
        neg = v < 0
        if rnd == RNDU:
            return -DBL_MAX_VAL if neg else INF
        if rnd == RNDD:
            return -INF if neg else DBL_MAX_VAL
        if rnd == RNDZ:
            return -DBL_MAX_VAL if neg else DBL_MAX_VAL
        # RNDN: the tie at the threshold goes to 2**1024 (even significand)
        return INF if (neg == False and a >= OVF_MID_V) else (
            -INF if (neg and a >= OVF_MID_V) else
            (-DBL_MAX_VAL if neg else DBL_MAX_VAL))
    with _ctx53(rnd):
        return float(mpfr(v, 53))


# ---------------------------------------------------------------- evaluators
# Each evaluator returns (lo, hi) at precision p with lo <= exact <= hi.
# `x` is an exact mpfr (a double, or the documented worst-case argument).

def _pi_pair(p):
    with _ctx(p, RNDD):
        lo = gmpy2.const_pi(p)
    with _ctx(p, RNDU):
        hi = gmpy2.const_pi(p)
    if lo == hi:                          # pi exactly representable (never here)
        return lo, hi
    return lo, hi


def _bracket(fn, x, p):
    with _ctx(p, RNDD):
        lo = fn(x)
    with _ctx(p, RNDU):
        hi = fn(x)
    return lo, hi


def _make_simple(fn):
    def ev(x, p):
        return _bracket(fn, x, p)
    return ev


def _bracket_pi_x(x, pp):
    """Valid bracket [tlo, thi] for pi*x, x a nonzero exact mpfr."""
    plo, phi = _pi_pair(pp)
    if x >= 0:
        with _ctx(pp, RNDD):
            tlo = plo * x
        with _ctx(pp, RNDU):
            thi = phi * x
    else:
        with _ctx(pp, RNDD):
            tlo = phi * x
        with _ctx(pp, RNDU):
            thi = plo * x
    return tlo, thi


def _critical_straddle(kind, tlo, thi, p):
    """Which critical points of sin/cos/tan the bracketed argument spans."""
    with _ctx(p + 64, RNDD):
        plo = gmpy2.const_pi(p + 64)
        if kind == "cos":                      # critical points: k*pi
            lo_k = gmpy2.floor(tlo / plo)
            hi_k = gmpy2.floor(thi / plo)
        else:                                  # critical points: pi/2 + k*pi
            half = plo / 2
            lo_k = gmpy2.floor((tlo - half) / plo)
            hi_k = gmpy2.floor((thi - half) / plo)
        k_lo, k_hi = int(lo_k), int(hi_k)
    return k_lo, k_hi


def _apply_trig(kind, tlo, thi, p):
    """Bracket sin/cos/tan over [tlo, thi], honouring critical points."""
    fn = {"sin": gmpy2.sin, "cos": gmpy2.cos, "tan": gmpy2.tan}[kind]
    with _ctx(p, RNDD):
        v1 = fn(tlo)
    with _ctx(p, RNDD):
        v2 = fn(thi)
    with _ctx(p, RNDU):
        w1 = fn(tlo)
    with _ctx(p, RNDU):
        w2 = fn(thi)
    lo, hi = min(v1, v2), max(w1, w2)
    # a critical point inside the interval makes the endpoints an invalid
    # bracket; widen it by the extreme value the function attains there
    k_lo, k_hi = _critical_straddle(kind, tlo, thi, p)
    if k_lo != k_hi:
        with _ctx(p, RNDU):
            hi = max(hi, gmpy2.mpfr(1))
        with _ctx(p, RNDD):
            lo = min(lo, gmpy2.mpfr(-1))
    return lo, hi


def _eval_trigpi(kind):
    """sin/cos/tan of pi*x, with pi bracketed.

    pi*x lands exactly on a zero or a critical point of the function when x is
    an integer or a half-integer; those arguments are decided analytically in
    reference() before this evaluator is reached.
    """
    def ev(x, p):
        pp = p + 64
        tlo, thi = _bracket_pi_x(x, pp)
        return _apply_trig(kind, tlo, thi, pp)
    return ev


def _special_pi(func, xd):
    """Analytic value for the degenerate arguments the manuals enumerate.

    Returns a 16-hex-digit string, or None when the argument is not special.
    """
    import math
    if not math.isfinite(xd):
        # documented: sinpi/cospi/tanpi return NaN on infinities and NaN
        return _double_to_hex(float("nan"))
    if xd != int(xd) and abs(xd * 2) != int(abs(xd * 2)):
        return None
    half = (abs(xd) * 2) % 2 == 1
    neg = xd < 0
    if func == "sinpi":
        if half:
            v = -1.0 if (int(abs(xd) - 0.5) % 2) else 1.0
            return _double_to_hex(-v if neg else v)
        return _double_to_hex(-0.0 if neg else 0.0)
    if func == "cospi":
        if half:
            # cos(pi(n+1/2)) = 0 exactly; the sign of the returned zero is a
            # convention, so these arguments are excluded from the comparison
            return CONVENTION
        n = int(abs(xd))
        v = 1.0 if n % 2 == 0 else -1.0
        return _double_to_hex(v)
    if func == "tanpi":
        if half:
            # tan at a pole is undefined; any infinity is a convention
            return CONVENTION
        return _double_to_hex(-0.0 if neg else 0.0)
    return None


def _eval_atanpi(x, p):
    # atanpi(x) = arctan(x)/pi
    with _ctx(p, RNDD):
        a_lo = gmpy2.atan(x)
    with _ctx(p, RNDU):
        a_hi = gmpy2.atan(x)
    plo, phi = _pi_pair(p)
    # division brackets: min/max over the four combinations, by sign of a
    vals = []
    with _ctx(p, RNDD):
        vals.append(a_lo / plo)
        vals.append(a_lo / phi)
        vals.append(a_hi / plo)
        vals.append(a_hi / phi)
    lo = min(vals)
    with _ctx(p, RNDU):
        vals = [a_lo / plo, a_lo / phi, a_hi / plo, a_hi / phi]
    hi = max(vals)
    return lo, hi


EVALUATORS = {
    "exp": _make_simple(gmpy2.exp),
    "log": _make_simple(gmpy2.log),
    "log2": _make_simple(gmpy2.log2),
    "log10": _make_simple(gmpy2.log10),
    "log1p": _make_simple(gmpy2.log1p),
    "expm1": _make_simple(gmpy2.expm1),
    "sin": _make_simple(gmpy2.sin),
    "cos": _make_simple(gmpy2.cos),
    "tan": _make_simple(gmpy2.tan),
    "asin": _make_simple(gmpy2.asin),
    "acos": _make_simple(gmpy2.acos),
    "atan": _make_simple(gmpy2.atan),
    "sinh": _make_simple(gmpy2.sinh),
    "cosh": _make_simple(gmpy2.cosh),
    "pow": lambda xy, p: _bracket(lambda t: t[0] ** t[1], xy, p),
    "sinpi": _eval_trigpi("sin"),
    "cospi": _eval_trigpi("cos"),
    "tanpi": _eval_trigpi("tan"),
    "atanpi": _eval_atanpi,
}


# ---------------------------------------------------------------- Ziv loop

def reference(func, xhex, max_prec=8192):
    """Return (rn, ru, rd, rz, status) as 16-hex-digit strings.

    status is 'decided' or 'undecided'.
    """
    if isinstance(xhex, tuple):
        x = tuple(mpfr(_hex_to_double(h)) for h in xhex)
    else:
        xd = _hex_to_double(xhex)
        if xd != xd:
            return None, None, None, None, "nan-arg"
        x = mpfr(xd)
    if func in ("sinpi", "cospi", "tanpi") and not isinstance(xhex, tuple):
        h = _special_pi(func, _hex_to_double(xhex))
        if h is CONVENTION:
            return None, None, None, None, "convention"
        if h is not None:
            return h, h, h, h, "special"
    ev = EVALUATORS[func]
    p = 64
    while p <= max_prec:
        lo, hi = ev(x, p)
        # For each mode the exact value's result is pinned once both ends of
        # the bracket round to the same double in that mode, because every
        # rounding map is monotone and exact lies in [lo, hi].
        res = {}
        for name, rnd in (("rn", RNDN), ("ru", RNDU), ("rd", RNDD),
                          ("rz", RNDZ)):
            with _ctx53(rnd):
                ra = float(mpfr(lo, 53))
                rb = float(mpfr(hi, 53))
            if _same_double(ra, rb):
                res[name] = ra
        if len(res) == 4:
            return (_double_to_hex(res["rn"]), _double_to_hex(res["ru"]),
                    _double_to_hex(res["rd"]), _double_to_hex(res["rz"]),
                    "decided")
        p *= 2
    return None, None, None, None, "undecided"


def _same_double(a, b):
    import struct
    return struct.pack(">d", a) == struct.pack(">d", b)


def _hex_to_double(h):
    import struct
    return struct.unpack(">d", bytes.fromhex(h))[0]


def _double_to_hex(d):
    import struct
    return struct.pack(">d", d).hex()


def signed_zero(x):
    import struct
    return (struct.unpack(">Q", struct.pack(">d", x))[0]) >> 63
