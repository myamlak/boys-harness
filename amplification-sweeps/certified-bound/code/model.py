"""The evaluated arithmetic of the v5.2.5 / v5.2.6 region-A entries, nominal
values and running-error radii.

Every nominal operation is the IEEE-754 double operation the kernel performs,
in the kernel's order; math.fma is the correctly rounded fused operation the
backend calls.  Each value carries an absolute radius e such that the exact
real value of the same expression at the exact input lies in
[nominal - e, nominal + e].  The radius rules are stated in the record.

The model is validated bit for bit against the sweep's cells.csv by replay.py.

The kernel checkout and the exponential bit file are selected with
configure(), from the command line or the environment (record_paths.py).  The
stored tables are read when first needed, so importing this module does not
require a checkout.
"""
import hashlib
import math
import os
import re
import struct
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import record_paths

U = 2.0 ** -53
U_ADD = U / (1.0 - U)
MIN_NORMAL = 2.0 ** -1022

# IEEE-754 round-to-nearest bounds the error of an underflowed result by half
# the least positive subnormal, 2**-1075.  That value is not a binary64
# number: it rounds to zero at ties-to-even, so a radius cannot carry it
# literally.  The rules below carry 2**-1071 = 16 * 2**-1075 instead: one
# half-ulp for the operation's own rounding and fifteen for the absolute
# errors of the radius arithmetic when its terms are themselves subnormal.
# The longest rule performs at most eleven roundings (four products, four
# bracket sums, the allowance addition, the operation and the SAFETY
# multiply); below the normal range the allowance addition is exact, so at
# most ten, and sixteen half-ulps leave a margin of at least five.
#
# The term is added when the computed result or the radius bracket lies below
# the normal range, 2**-1022.  Above it, the relative u_+|z| term and SAFETY
# cover the same invariant half-ulp: SAFETY's inflation of a normal-range
# bracket is at least 2**-1022 * 2**-40 = 2**-1062, thirteen binary orders
# above 2**-1075.  A computed result in the subnormal range with a large
# bracket (a large propagated radius) is likewise covered by SAFETY.
UNDERFLOW = 2.0 ** -1071

# Conservative expansion of every accumulated radius: the radius arithmetic
# itself runs in double, and this covers its rounding and the small quadratic
# terms dropped by inspection.
SAFETY = 1.0 + 2.0 ** -40

# Contraction of the bare recurrence expressions, as measured by replay.py.
UPWARD_FUSED = False      # (l - 0.5) * f - expx
DOWNWARD_FUSED = None     # x * f + expx: decided by the replay
PRODUCT_PASS_FUSED = False  # pass[k] = c * basis folded into the pairwise sum

# boys_coefficients.hpp differs between v5.1.1/v5.2.4 and v5.2.5 in the
# rational region-B seed, which the named policy does not read.  The model
# therefore identifies the file by a digest of exactly the fields it parses,
# which is identical in those trees and in v5.2.6, not the whole-file digest.
KNOWN_FIT_DIGESTS = {
    "8ce5e17157a4be5f621d695f377094d212cb839732e278d47f23aa225b4d9204":
        "v5.1.1 / v5.2.4",
    "d978d88bcfbb60727cab3eed32039c9d66ce917a14c7f5a199b6215db8f7ea00":
        "v5.2.5 / v5.2.6",
}
FIT_CONTENT_DIGEST = (
    "74fe52545be4540cfd7abea1b19e6bd5d4e681b656a1be6603585bb97dafa99c")
FIT_D_DIGEST = "7af73817161c34144f149f67a577a275c42e10c6daeb2bdb2245d1367556cf8e"

_KERNEL = None
_EXP_BITS = None
_TABLE = None


def configure(kernel=None, exp_bits=None):
    """Select the kernel checkout and the exponential bit file.

    ``kernel`` is the checkout root; the headers are read from its
    include/boys.  Paths are resolved by record_paths.kernel/exp_bits, so
    None means the documented environment variable or default.
    """
    global _KERNEL, _EXP_BITS, _TABLE
    _KERNEL = record_paths.kernel(kernel)
    _EXP_BITS = record_paths.exp_bits(exp_bits)
    _TABLE = None


def _kernel_root():
    global _KERNEL
    if _KERNEL is None:
        _KERNEL = record_paths.kernel()
    return _KERNEL


def _exp_bits_path():
    global _EXP_BITS
    if _EXP_BITS is None:
        _EXP_BITS = record_paths.exp_bits()
    return _EXP_BITS


def _boys_dir():
    return _kernel_root() / "include" / "boys"


def _read(name):
    with open(_boys_dir() / name, "r", encoding="utf-8") as handle:
        return handle.read()


def _digest(name):
    with open(_boys_dir() / name, "rb") as handle:
        return hashlib.sha256(handle.read()).hexdigest()


def _array_body(text, name):
    start = text.index(name)
    open_token = text.index("({", start)
    depth = 0
    for i in range(open_token + 1, len(text)):
        ch = text[i]
        if ch in "({":
            depth += 1
        elif ch in ")}":
            depth -= 1
            if depth == 0:
                return text[open_token + 2: i]
    raise ValueError(f"unterminated array {name}")


def _doubles(body):
    return [float(tok) for tok in re.findall(r"[-+]?\d+\.\d+e[+-]\d+", body)]


def _ints(body):
    return [int(tok) for tok in re.findall(r"-?\d+", body)]


def _content_digest(pieces, piece_start, thresholds, ext_coeffs, ext_lo, x0):
    """Digest of exactly the stored tables the model reads, in a canonical
    form (float.hex, so the digest does not depend on Python's float repr)."""
    h = hashlib.sha256()

    def line(text):
        h.update(text.encode())
        h.update(b"\n")

    line("pieces")
    for p in pieces:
        line(f"{p['a'].hex()} {p['b'].hex()} {p['deg']} {p['off']}")
    line("piece_start")
    line(" ".join(str(v) for v in piece_start))
    line("thresholds")
    line(" ".join(v.hex() for v in thresholds))
    line("ext_coeffs")
    line(" ".join(v.hex() for v in ext_coeffs))
    line("ext_lo")
    line(ext_lo.hex())
    line("x0")
    line(x0.hex())
    return h.hexdigest()


def load_tables():
    master = _read("boys_coefficients.hpp")
    cheb = _read("boys_coefficients_d_coarse_cheb.hpp")
    assert _digest("boys_coefficients_d_coarse_cheb.hpp") == FIT_D_DIGEST, (
        "the coarse Chebyshev table is not the sweep's; expected "
        f"{FIT_D_DIGEST}")

    pieces_text = _array_body(master, "inline constexpr auto kPieces")
    rows = re.findall(
        r"\{\s*([-+0-9.eE]+)\s*,\s*([-+0-9.eE]+)\s*,\s*(\d+)\s*,\s*(\d+)\s*\}",
        pieces_text)
    pieces = [dict(a=float(a), b=float(b), deg=int(d), off=int(o))
              for a, b, d, o in rows]
    piece_start = _ints(_array_body(master, "inline constexpr auto kPieceStart"))
    thresholds = _doubles(_array_body(master, "inline constexpr auto kTierThresholds"))
    coeffs = _doubles(_array_body(cheb, "kCoeffs = std::to_array<double>"))
    ext_coeffs = _doubles(_array_body(master, "inline constexpr auto kExtendedBcoeffs"))
    ext_lo = float(re.search(
        r"inline constexpr double kExtendedBX0 = ([-+0-9.eE]+);", master).group(1))
    x0 = float(re.search(
        r"inline constexpr double kX0 = ([-+0-9.eE]+);", master).group(1))
    assert len(pieces) == 66 and len(piece_start) == 34
    assert len(coeffs) == 1320 and len(thresholds) == 33
    assert len(ext_coeffs) == 25

    digest = _content_digest(pieces, piece_start, thresholds, ext_coeffs,
                             ext_lo, x0)
    assert digest == FIT_CONTENT_DIGEST, (
        "the parsed stored tables are not the sweep's: content digest "
        f"{digest}, expected {FIT_CONTENT_DIGEST}")
    whole = _digest("boys_coefficients.hpp")
    assert whole in KNOWN_FIT_DIGESTS, (
        f"boys_coefficients.hpp digest {whole} is not a known release")
    return dict(pieces=pieces, piece_start=piece_start, thresholds=thresholds,
                coeffs=coeffs, ext_coeffs=ext_coeffs, ext_lo=ext_lo, ext_hi=x0,
                ext_deg=24, fit_digest=whole,
                release=KNOWN_FIT_DIGESTS[whole])


def table():
    """The stored tables, read on first use."""
    global _TABLE
    if _TABLE is None:
        _TABLE = load_tables()
    return _TABLE


def theta0():
    return table()["thresholds"][0]


def piece_index(order, band_index):
    return table()["piece_start"][order] + band_index


def piece(order, band_index):
    p = table()["pieces"][piece_index(order, band_index)]
    return p["a"], p["b"], table()["coeffs"][p["off"]: p["off"] + p["deg"] + 1]


def ext_piece():
    t = table()
    return t["ext_lo"], t["ext_hi"], t["ext_coeffs"]


# ---------------------------------------------------------------------------
# primitive operations with radii
# ---------------------------------------------------------------------------

def _allowance(z, bracket):
    """The absolute underflow allowance of the rule's comment.

    Added when the computed result or the radius bracket is below the normal
    range; zero otherwise, which leaves every normal-range radius unchanged.
    """
    if abs(z) < MIN_NORMAL or bracket < MIN_NORMAL:
        return UNDERFLOW
    return 0.0


def add(a, ea, b, eb):
    """Addition.  The exact sum of two binary64 values is a multiple of
    2**-1074; when it lies below 2**-1022 it is therefore exactly
    representable, so the operation itself never loses a subnormal result.
    The allowance covers the radius arithmetic's own subnormal roundings."""
    z = a + b
    bracket = U_ADD * abs(z) + ea + eb
    return z, (bracket + _allowance(z, bracket)) * SAFETY


def sub(a, ea, b, eb):
    """Subtraction, with add's exactness property."""
    z = a - b
    bracket = U_ADD * abs(z) + ea + eb
    return z, (bracket + _allowance(z, bracket)) * SAFETY


def mul(a, ea, b, eb):
    """Multiplication.  |RN(ab) - ab| <= max(u|ab|, 2**-1075); the first term
    is u_+|z| for a normal z and the allowance carries the second below the
    normal range, where u_+|z| no longer dominates it."""
    z = a * b
    bracket = U_ADD * abs(z) + abs(b) * ea + abs(a) * eb + ea * eb
    return z, (bracket + _allowance(z, bracket)) * SAFETY


def fma(a, ea, b, eb, c, ec):
    """Fused multiply-add, with mul's underflow bound applied to ab + c."""
    z = math.fma(a, b, c)
    bracket = U_ADD * abs(z) + abs(b) * ea + abs(a) * eb + ec + ea * eb
    return z, (bracket + _allowance(z, bracket)) * SAFETY


def scale2(a, ea):
    """2*a, exact for every finite a without overflow: it changes only the
    exponent, so the radius doubles exactly."""
    return 2.0 * a, 2.0 * ea


def refined_div(a, ea, x, invx):
    """DivideByReciprocal: q0 = a*invx; r = fma(-q0,x,a); q = fma(r,invx,q0).

    The refinement cancels the reciprocal's and the product's first-order
    errors, so the result is the correctly rounded quotient up to a second
    order term (derived in the record): for exact numerator A,
    |q - A/x| <= u(1+8u)|A/x| + ea/|x|.

    The derivation assumes every rounding it charges is relative, which
    requires the denominator and the quotient (and the two intermediates) to
    be normal or exactly zero.  The model's calls satisfy this: the downward
    denominator is l + 1/2 in [1/2, 63/2] and the upward one is x >= theta0,
    and the delivered values never approach the subnormal range.  The check
    below refuses a call outside that domain rather than undercharge it.
    """
    q0 = a * invx
    r = math.fma(-q0, x, a)
    q = math.fma(r, invx, q0)
    if (abs(x) < MIN_NORMAL
            or (q != 0.0 and abs(q) < MIN_NORMAL)
            or (q0 != 0.0 and abs(q0) < MIN_NORMAL)
            or (r != 0.0 and abs(r) < MIN_NORMAL)
            or (q == 0.0 and a != 0.0)):
        raise ArithmeticError(
            "refined division outside its derived range: the quotient and "
            "its intermediates must be normal or exactly zero; got "
            f"a={a!r} x={x!r} q0={q0!r} r={r!r} q={q!r}")
    ud = U * (1.0 + 8.0 * U) / (1.0 - U * (1.0 + 8.0 * U))
    return q, (ud * (abs(q) + ea / abs(x)) + ea / abs(x)) * SAFETY


# ---------------------------------------------------------------------------
# mapping
# ---------------------------------------------------------------------------

def mapping(x, a, b):
    """t = 2.0 * (x - a) / (b - a) - 1.0 with its radius."""
    t1, e1 = sub(x, 0.0, a, 0.0)
    t2 = 2.0 * t1
    e2 = 2.0 * e1
    den, ed = sub(b, 0.0, a, 0.0)
    q = t2 / den
    # |q - T2/DEN| <= u|q| + e2/|den| + |t2| ed/|den|^2, with the division's
    # own underflow bound carried as for mul.
    bracket = U_ADD * abs(q) + e2 / abs(den) + abs(t2) * ed / (den * den)
    eq = (bracket + _allowance(q, bracket)) * SAFETY
    t, et = sub(q, eq, 1.0, 0.0)
    return t, et


# ---------------------------------------------------------------------------
# split Clenshaw, the kernel's ClenshawSplit
# ---------------------------------------------------------------------------

def clenshaw_split(t, et, c):
    """ClenshawSplit<ScalarFp64> at deg = len(c)-1 >= 4.  Returns value, radius."""
    deg = len(c) - 1
    assert deg >= 4
    t2, e_t2 = mul(t, et, t, et)
    v, e_v = fma(2.0, 0.0, t2, e_t2, -1.0, 0.0)
    two_v, e_two_v = scale2(v, e_v)

    m = deg // 2
    mo = (deg + 1) // 2
    b1, e_b1 = c[2 * m], 0.0
    b2, e_b2 = 0.0, 0.0
    for k in range(m - 1, 0, -1):
        s, e_s = sub(c[2 * k], 0.0, b2, e_b2)
        b0, e_b0 = fma(two_v, e_two_v, b1, e_b1, s, e_s)
        b2, e_b2 = b1, e_b1
        b1, e_b1 = b0, e_b0

    s, e_s = sub(c[0], 0.0, b2, e_b2)
    even, e_even = fma(v, e_v, b1, e_b1, s, e_s)

    o1, e_o1 = c[2 * mo - 1], 0.0
    o2, e_o2 = 0.0, 0.0
    for k in range(mo - 2, 0, -1):
        s, e_s = sub(c[2 * k + 1], 0.0, o2, e_o2)
        o0, e_o0 = fma(two_v, e_two_v, o1, e_o1, s, e_s)
        o2, e_o2 = o1, e_o1
        o1, e_o1 = o0, e_o0

    s, e_s = sub(two_v, e_two_v, 1.0, 0.0)
    u, e_u = sub(c[1], 0.0, o2, e_o2)
    odd, e_odd = fma(s, e_s, o1, e_o1, u, e_u)
    ret, e_ret = fma(t, et, odd, e_odd, even, e_even)
    return ret, e_ret


# ---------------------------------------------------------------------------
# per-order fit entry points
# ---------------------------------------------------------------------------

def eval_fit_split(order, band_index, x):
    a, b, c = piece(order, band_index)
    t, et = mapping(x, a, b)
    return clenshaw_split(t, et, c)


def eval_ext_seed(x):
    a, b, c = ext_piece()
    t, et = mapping(x, a, b)
    return clenshaw_split(t, et, c)


def eval_fit_product(order, band_index, x, basis=None):
    """The direct lane's product form for one order.  basis may be shared."""
    a, b, c = piece(order, band_index)
    if basis is None:
        basis = product_basis(x, a, b, len(c) - 1)
    deg = len(c) - 1
    pass_v = []
    pass_e = []
    for k in range(len(c)):
        degree = deg - k
        z, ez = mul(c[degree], 0.0, basis[degree][0], basis[degree][1])
        pass_v.append(z)
        pass_e.append(ez)
    width = len(pass_v)
    while width > 1:
        for k in range(width // 2):
            z, ez = add(pass_v[2 * k], pass_e[2 * k],
                        pass_v[2 * k + 1], pass_e[2 * k + 1])
            pass_v[k], pass_e[k] = z, ez
        if width % 2 != 0:
            pass_v[width // 2] = pass_v[width - 1]
            pass_e[width // 2] = pass_e[width - 1]
        width = (width + 1) // 2
    return pass_v[0], pass_e[0]


def product_basis(x, a, b, deg):
    t, et = mapping(x, a, b)
    basis = [(1.0, 0.0), (t, et)]
    m2, e_m2 = scale2(t, et)
    for _ in range(2, deg + 1):
        upper_v, upper_e = basis[-1]
        lower_v, lower_e = basis[-2]
        prod, e_prod = fma(m2, e_m2, upper_v, upper_e, 0.0, 0.0)
        next_v, e_next = sub(prod, e_prod, lower_v, lower_e)
        basis.append((next_v, e_next))
    return basis


# ---------------------------------------------------------------------------
# the recursed lane's batch
# ---------------------------------------------------------------------------

def _upward(f, ef, x, invx, expx, l):
    if UPWARD_FUSED:
        a, ea = fma(float(l) - 0.5, 0.0, f, ef, -expx, 0.0)
    else:
        p, ep = mul(float(l) - 0.5, 0.0, f, ef)
        a, ea = sub(p, ep, expx, 0.0)
    return refined_div(a, ea, x, invx)


def _downward(f, ef, x, expx, invd, l):
    if DOWNWARD_FUSED:
        a, ea = fma(x, 0.0, f, ef, expx, 0.0)
    else:
        p, ep = mul(x, 0.0, f, ef)
        a, ea = add(p, ep, expx, 0.0)
    return refined_div(a, ea, float(l) + 0.5, invd)


def exp_sensitivity_upward(n, x):
    """|dF_n/d(expx)| for the upward recurrence, exact linear coefficient."""
    d = 0.0
    for l in range(1, n + 1):
        d = (float(l) - 0.5) * d - 1.0
        # divide by x; written as d/x, one rounding
        d = d / x
    return abs(d)


def exp_sensitivity_downward(n, x):
    """|dF_n/d(expx)| for the downward recurrence from order 32."""
    d = 0.0
    for l in range(31, n - 1, -1):
        d = (x * d + 1.0) / (float(l) + 0.5)
    return abs(d)


def ship_batch(x, expx):
    """All 33 recursed orders at one argument, with radii and mode.

    Returns list of dicts: value, radius, mode, c_seed, c_rec, c_exp,
    ratio_seed_amp (amplified seed radius), min_abs.
    """
    out = []
    if x == 0.0:
        for n in range(33):
            # The kernel returns 1.0 / (2.0 * l + 1.0): one correctly rounded
            # division, so the delivered value differs from the exact 1/(2n+1)
            # by at most u/(1-u) times itself.
            v = 1.0 / (2.0 * n + 1.0)
            out.append(dict(value=v, radius=U_ADD * abs(v) * SAFETY,
                            mode="exact", c_seed=0.0,
                            c_rec=U_ADD * abs(v) * SAFETY, c_exp=0.0,
                            seed_amp=0.0, min_abs=v))
        return out

    if x < theta0():
        seed_v, seed_e = eval_fit_split(32, 0, x)
        # the seed's amplified contribution at order n = (A_n/A_32) * seed_e.
        # A_n/A_32 = prod_{l=n+1..32} (l-0.5)/x.
        amp = [1.0] * 33
        for n in range(31, -1, -1):
            amp[n] = amp[n + 1] * (x / (float(n) + 0.5))
        f, ef = seed_v, seed_e
        vals = [0.0] * 33
        rads = [0.0] * 33
        mins = [0.0] * 33
        vals[32], rads[32], mins[32] = f, ef, abs(f)
        for l in range(31, -1, -1):
            invd = 1.0 / (float(l) + 0.5)
            f, ef = _downward(f, ef, x, expx, invd, l)
            vals[l], rads[l], mins[l] = f, ef, abs(f)
        # The radius in the batch carries the whole chain.  The exp part is
        # analytic and linear: subtract it from the propagated radius to
        # report it separately; add it back for the certified total.
        for n in range(33):
            # SAFETY covers the sensitivity's own double-precision evaluation
            # (n roundings at u, well inside 2^-40).
            ce = exp_sensitivity_downward(n, x) * expx * (2.0 ** -52) * SAFETY
            out.append(dict(value=vals[n], radius=rads[n] + ce,
                            mode="downward", c_seed=seed_e * amp[n],
                            c_rec=max(0.0, rads[n] - seed_e * amp[n]),
                            c_exp=ce, seed_amp=seed_e * amp[n],
                            min_abs=mins[n]))
        return out

    source_f, source_e = eval_ext_seed(x)
    invx = 1.0 / x
    # served orders, exactly the kernel's prefix
    nmax = 32
    thresholds = table()["thresholds"]
    served = 0
    while served < nmax and x >= thresholds[served + 1]:
        served += 1
    # ladder with radii (exp not in the radii; added analytically)
    vals = [0.0] * (served + 1)
    rads = [0.0] * (served + 1)
    mins = [0.0] * (served + 1)
    f, ef = source_f, source_e
    vals[0], rads[0], mins[0] = f, ef, abs(f)
    for l in range(1, served + 1):
        f, ef = _upward(f, ef, x, invx, expx, l)
        vals[l], rads[l], mins[l] = f, ef, abs(f)
    for n in range(33):
        if n <= served:
            # SAFETY as in the downward branch: the sensitivity is evaluated
            # in double, so its own rounding is inside 2^-40.
            ce = exp_sensitivity_upward(n, x) * expx * (2.0 ** -52) * SAFETY
            out.append(dict(value=vals[n], radius=rads[n] + ce,
                            mode="recursed" if n else "seed",
                            c_seed=source_e * _amp_up(n, x),
                            c_rec=max(0.0, rads[n] - source_e * _amp_up(n, x)),
                            c_exp=ce, seed_amp=source_e * _amp_up(n, x),
                            min_abs=mins[n]))
        else:
            v, e = eval_fit_split(n, 0 if x < table()["pieces"][0]["b"] else 1, x)
            out.append(dict(value=v, radius=e, mode="same-fit",
                            c_seed=0.0, c_rec=e, c_exp=0.0, seed_amp=0.0,
                            min_abs=abs(v)))
    return out


def _amp_up(n, x):
    """A_n = (2n-1)!!/(2x)^n."""
    a = 1.0
    for l in range(1, n + 1):
        a *= (2.0 * l - 1.0) / (2.0 * x)
    return a


def product_batch(x, band_index):
    a, b, c0 = piece(0, band_index)
    deg = len(c0) - 1
    basis = product_basis(x, a, b, deg)
    out = []
    for n in range(33):
        v, e = eval_fit_product(n, band_index, x, basis)
        out.append(dict(value=v, radius=e, mode="product", min_abs=abs(v)))
    return out


# ---------------------------------------------------------------------------
# probe exponentials
# ---------------------------------------------------------------------------

def load_exp():
    exp = {}
    with open(_exp_bits_path(), "r") as handle:
        for line in handle:
            xb_s, eb_s, hb_s = line.split()
            xb = int(xb_s, 16)
            x = struct.unpack(">d", struct.pack(">Q", xb))[0]
            exp[x] = (struct.unpack(">d", struct.pack(">Q", int(eb_s, 16)))[0],
                      struct.unpack(">d", struct.pack(">Q", int(hb_s, 16)))[0])
    return exp


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    record_paths.add_kernel(parser)
    parser.add_argument("--exp-bits", metavar="FILE", help="probe output file")
    args = parser.parse_args()
    configure(kernel=args.kernel, exp_bits=args.exp_bits)
    print("tables", len(table()["pieces"]), table()["ext_lo"], table()["ext_hi"],
          "release", table()["release"])
    e = load_exp()
    print("exp table", len(e))
    x = 0.88881931955209315
    print("mapping", mapping(x, 0.0, table()["pieces"][0]["b"]))
    print("ship x<theta", ship_batch(x, e[x][1])[0])
