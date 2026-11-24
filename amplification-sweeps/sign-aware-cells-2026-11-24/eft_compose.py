"""EFT composition and double-double cross-check at the two chosen cells.

Purpose
-------
certify_cells.py shows the sign of the rounding is certifiable from exact
rational arithmetic.  This script answers the implementation question: what
arithmetic primitive would a library need at runtime?

1. EFT composition.  It replays each entry's operation graph carrying, at
   every node, both the exact real value V of the node at the delivered
   double inputs x and expx (a Fraction) and the delivered double d.  At every
   node the local rounding residual loc = exact_op(d_children) - d_node is
   computed exactly and checked against an error-free transformation:
   TwoSum/TwoDiff for add/sub, TwoProd for mul, TwoProd+TwoSum for fma, and
   the S18.4 identity for the refined division.  The node values compose
   exactly, so the final Delta = V - d is the negative of the fixed-input
   rounding residual.  This is the "EFT composed along the recurrence and the
   seed fit".

2. Sign-aware refined division.  For every refined-division call it verifies
   the extension of the S18.4 lemma,

       q - a/x = (a/x) * [delta3 - theta1*theta2*(1+delta3)],

   with delta0..delta3 the exact per-operation relative residuals.  The
   bracket is an exactly computable two-term quantity; its sign decides the
   division error's sign.  The magnitude bound of S18.4 cannot.

3. Double-double route.  It evaluates the same exact route in double-double
   and compares with the exact rational route, showing the DD error against
   the margin each cell must resolve.

The deposited columns are used only for the bit checks and the final
comparison, never as certificate inputs.

Run:  python -I eft_compose.py
Writes: eft_compose.json
"""
import csv
import json
import math
import os
import sys
from fractions import Fraction
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(HERE))
import exact_arith as ea                     # noqa: E402

RECORD = REPO / "amplification-sweeps" / "certified-bound"
SWEEP = Path(os.environ.get("CERT_SWEEP", REPO / "amplification-sweeps" / "amplification-form"))
KERNEL = Path(os.environ.get("BOYS_ROOT", REPO.parent / "boys"))
sys.path.insert(0, str(RECORD / "code"))
import model                                 # noqa: E402
import exactroute as er                      # noqa: E402

U = 2.0 ** -53


class Val:
    """A node: exact real value V at the delivered inputs, delivered double d."""
    __slots__ = ("V", "d")

    def __init__(self, V, d):
        self.V = V
        self.d = d


class Machine:
    """The operation graph with exact route values and EFT residual checks."""

    def __init__(self):
        self.count = {}
        self.sum_abs_loc = Fraction(0)
        self.sum_half_ulp = 0.0
        self.max_abs_loc = Fraction(0)
        self.div_calls = []
        self.max_out = 0.0

    def _node(self, kind, V, d, loc, d_args):
        self.count[kind] = self.count.get(kind, 0) + 1
        self.sum_abs_loc += abs(loc)
        self.max_abs_loc = max(self.max_abs_loc, abs(loc))
        self.sum_half_ulp += 0.5 * math.ulp(d)
        self.max_out = max(self.max_out, abs(d))
        for a in d_args:
            if math.isfinite(a):
                self.max_out = max(self.max_out, abs(a))
        return Val(V, d)

    # -- primitives --------------------------------------------------------
    def add(self, x, y):
        V = x.V + y.V
        z = x.d + y.d
        loc = (Fraction(x.d) + Fraction(y.d)) - Fraction(z)
        s, e = ea.two_sum(x.d, y.d)
        assert s == z and Fraction(e) == loc
        return self._node("add", V, z, loc, (x.d, y.d))

    def sub(self, x, y):
        V = x.V - y.V
        z = x.d - y.d
        loc = (Fraction(x.d) - Fraction(y.d)) - Fraction(z)
        s, e = ea.two_diff(x.d, y.d)
        assert s == z and Fraction(e) == loc
        return self._node("sub", V, z, loc, (x.d, y.d))

    def mul(self, x, y):
        V = x.V * y.V
        z = x.d * y.d
        loc = Fraction(x.d) * Fraction(y.d) - Fraction(z)
        p, e = ea.two_prod(x.d, y.d)
        assert p == z and Fraction(e) == loc
        return self._node("mul", V, z, loc, (x.d, y.d))

    def fma(self, u, v, w):
        V = u.V * v.V + w.V
        z = math.fma(u.d, v.d, w.d)
        loc = Fraction(u.d) * Fraction(v.d) + Fraction(w.d) - Fraction(z)
        p, ep = ea.two_prod(u.d, v.d)           # p + ep = u.d*v.d
        s, es = ea.two_sum(p, w.d)              # s + es = p + w.d
        assert Fraction(s) + Fraction(es) + Fraction(ep) - Fraction(z) == loc
        return self._node("fma", V, z, loc, (u.d, v.d, w.d))

    def scale2(self, x):
        V = 2 * x.V
        z = 2.0 * x.d
        assert Fraction(z) == 2 * Fraction(x.d)   # exact exponent change
        return self._node("scale2", V, z, Fraction(0), (x.d,))

    def div(self, x, y):
        V = x.V / y.V
        z = x.d / y.d
        loc = Fraction(x.d) / Fraction(y.d) - Fraction(z)
        return self._node("div", V, z, loc, (x.d, y.d))

    def refined_div(self, a, x, invx):
        """RefinedDivide(a, x, invx), with the S18.4 residual identity checked."""
        V = a.V / x.V
        q0 = a.d * invx.d
        r = math.fma(-q0, x.d, a.d)
        q = math.fma(r, invx.d, q0)
        loc = Fraction(a.d) / Fraction(x.d) - Fraction(q)

        delta0 = Fraction(invx.d) * Fraction(x.d) - 1
        ag = Fraction(a.d) * Fraction(invx.d)
        delta1 = Fraction(q0) / ag - 1
        theta1 = (1 + delta0) * (1 + delta1) - 1
        resid = Fraction(a.d) - Fraction(q0) * Fraction(x.d)
        if resid == 0:
            delta2 = delta3 = delta_div = None
            theta2 = None
            exact_identity = (Fraction(q0) == Fraction(a.d) / Fraction(x.d)
                              and loc == 0)
        else:
            delta2 = Fraction(r) / resid - 1
            theta2 = (1 + delta0) * (1 + delta2) - 1
            sum_rg = Fraction(r) * Fraction(invx.d) + Fraction(q0)
            delta3 = Fraction(q) / sum_rg - 1
            delta_div = delta3 - theta1 * theta2 * (1 + delta3)
            # loc is exact-minus-delivered, so loc = -(a/x)*delta_div.
            exact_identity = (Fraction(a.d) / Fraction(x.d)) * (-delta_div) == loc
        self.div_calls.append(dict(
            a=a.d, x=x.d, invx=invx.d, q0=q0, r=r, q=q,
            delta0=float(delta0),
            delta1=None if delta1 is None else float(delta1),
            delta2=None if delta2 is None else float(delta2),
            delta3=None if delta3 is None else float(delta3),
            abs_delta_div=None if delta_div is None else float(abs(delta_div)),
            loc=float(loc),
            exact_identity=bool(exact_identity),
            # sign of the delivered error q - a/x is sign(delta_div)
            q_error_sign=(0 if delta_div is None
                          else (delta_div > 0) - (delta_div < 0)),
        ))
        assert exact_identity, (a.d, x.d, invx.d, q0, r, q, loc,
                                delta0, delta1, delta2, delta3, delta_div,
                                Fraction(a.d) / Fraction(x.d))
        return self._node("refined_div", V, q, loc, (a.d, x.d, invx.d))


# ---------------------------------------------------------------------------
# operation graphs (mirror model.py; identical to exact_arith's value paths)
# ---------------------------------------------------------------------------

def mapping(m, x, a, b):
    t1 = m.sub(x, Val(Fraction(a), a))
    t2 = m.scale2(t1)
    den = m.sub(Val(Fraction(b), b), Val(Fraction(a), a))
    q = m.div(t2, den)
    return m.sub(q, Val(Fraction(1), 1.0))


def clenshaw_split(m, t, c):
    deg = len(c) - 1
    t2 = m.mul(t, t)
    v = m.fma(Val(Fraction(2), 2.0), t2, Val(Fraction(-1), -1.0))
    two_v = m.scale2(v)
    mm = deg // 2
    mo = (deg + 1) // 2
    co = [Val(Fraction(c[k]), c[k]) for k in range(len(c))]
    b1 = co[2 * mm]
    b2 = Val(Fraction(0), 0.0)
    for k in range(mm - 1, 0, -1):
        s = m.sub(co[2 * k], b2)
        b0 = m.fma(two_v, b1, s)
        b2, b1 = b1, b0
    s = m.sub(co[0], b2)
    even = m.fma(v, b1, s)
    o1 = co[2 * mo - 1]
    o2 = Val(Fraction(0), 0.0)
    for k in range(mo - 2, 0, -1):
        s = m.sub(co[2 * k + 1], o2)
        o0 = m.fma(two_v, o1, s)
        o2, o1 = o1, o0
    s = m.sub(two_v, Val(Fraction(1), 1.0))
    u = m.sub(co[1], o2)
    odd = m.fma(s, o1, u)
    return m.fma(t, odd, even)


def product_basis(m, x, a, b, deg):
    t = mapping(m, x, a, b)
    basis = [Val(Fraction(1), 1.0), t]
    m2 = m.scale2(t)
    for _ in range(2, deg + 1):
        prod = m.fma(m2, basis[-1], Val(Fraction(0), 0.0))
        basis.append(m.sub(prod, basis[-2]))
    return basis


def eval_fit_product(m, c, basis):
    deg = len(c) - 1
    pass_v = []
    for k in range(len(c)):
        degree = deg - k
        pass_v.append(m.mul(Val(Fraction(c[degree]), c[degree]),
                            basis[degree]))
    width = len(pass_v)
    while width > 1:
        for k in range(width // 2):
            pass_v[k] = m.add(pass_v[2 * k], pass_v[2 * k + 1])
        if width % 2 != 0:
            pass_v[width // 2] = pass_v[width - 1]
        width = (width + 1) // 2
    return pass_v[0]


def ship_machine(x, expx, n):
    m = Machine()
    xv = Val(Fraction(x), x)
    ev = Val(Fraction(expx), expx)
    a, b, c = model.ext_piece()
    f = clenshaw_split(m, mapping(m, xv, a, b), c)
    invx = 1.0 / x
    for l in range(1, n + 1):
        p = m.mul(Val(Fraction(2 * l - 1, 2), float(l) - 0.5), f)
        aa = m.sub(p, ev)
        f = m.refined_div(aa, xv, Val(Fraction(invx), invx))
    return m, f


def product_machine(x, band, n):
    m = Machine()
    a, b, c0 = model.piece(0, band)
    xv = Val(Fraction(x), x)
    basis = product_basis(m, xv, a, b, len(c0) - 1)
    return m, eval_fit_product(m, model.piece(n, band)[2], basis)


# ---------------------------------------------------------------------------
# double-double
# ---------------------------------------------------------------------------

def quick_two_sum(a, b):
    s = a + b
    e = b - (s - a)
    return s, e


class DD:
    __slots__ = ("hi", "lo")

    def __init__(self, hi, lo=0.0):
        self.hi = hi
        self.lo = lo


def dd(x):
    return DD(float(x), 0.0)


def dd_add(a, b):
    s1, s2 = ea.two_sum(a.hi, b.hi)
    t1, t2 = ea.two_sum(a.lo, b.lo)
    s2 += t1
    s1, s2 = quick_two_sum(s1, s2)
    s2 += t2
    s1, s2 = quick_two_sum(s1, s2)
    return DD(s1, s2)


def dd_sub(a, b):
    return dd_add(a, DD(-b.hi, -b.lo))


def dd_mul(a, b):
    p1, p2 = ea.two_prod(a.hi, b.hi)
    p2 += a.hi * b.lo + a.lo * b.hi
    p1, p2 = quick_two_sum(p1, p2)
    return DD(p1, p2)


def dd_div(a, b):
    q1 = a.hi / b.hi
    r = dd_sub(a, dd_mul(dd(q1), b))
    q2 = r.hi / b.hi
    q = dd_add(dd(q1), dd(q2))
    r = dd_sub(a, dd_mul(q, b))
    q3 = r.hi / b.hi
    return dd_add(q, dd(q3))


def dd_half(fr):
    """DD approximation of a Fraction; returns (dd, exact remainder bound)."""
    hi = float(fr)                   # CPython int division is correctly rounded
    rem = fr - Fraction(hi)
    lo = float(rem)
    rem2 = rem - Fraction(lo)
    return DD(hi, lo), abs(rem2)


def mapping_dd(x, a, b):
    xv = dd(x)
    av = dd(a)
    bv = dd(b)
    return dd_sub(dd_div(dd_mul(dd(2.0), dd_sub(xv, av)), dd_sub(bv, av)),
                  dd(1.0))


def cheb_sum_dd(c, t):
    T1 = dd(1.0)
    total = dd_mul(dd(c[0]), T1)
    if len(c) == 1:
        return total
    T0 = t
    total = dd_add(total, dd_mul(dd(c[1]), T0))
    for k in range(2, len(c)):
        Tk = dd_sub(dd_mul(dd(2.0), dd_mul(t, T0)), T1)
        total = dd_add(total, dd_mul(dd(c[k]), Tk))
        T1, T0 = T0, Tk
    return total


def route_dd_ship(x, n, eh_dd, maxabs):
    a, b, c = model.ext_piece()
    f = cheb_sum_dd(c, mapping_dd(x, a, b))
    maxabs[0] = max(maxabs[0], abs(f.hi))
    xv = dd(x)
    for l in range(1, n + 1):
        p = dd_mul(dd(l - 0.5), f)
        num = dd_sub(p, eh_dd)
        f = dd_div(num, xv)
        maxabs[0] = max(maxabs[0], abs(f.hi))
    return f


def route_dd_product(x, band, n, maxabs):
    a, b, c0 = model.piece(0, band)
    c = model.piece(n, band)[2]
    deg = len(c0) - 1
    t = mapping_dd(x, a, b)
    T = [dd(1.0), t]
    m2 = dd_mul(dd(2.0), t)
    for k in range(2, deg + 1):
        prod = dd_mul(m2, T[-1])
        T.append(dd_sub(prod, T[-2]))
        maxabs[0] = max(maxabs[0], abs(T[-1].hi))
    total = dd(0.0)
    for k in range(len(c)):
        total = dd_add(total, dd_mul(dd(c[k]), T[k]))
        maxabs[0] = max(maxabs[0], abs(total.hi))
    return total


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def main():
    model.configure(kernel=KERNEL, exp_bits=RECORD / "exp" / "exp_bits.txt")
    harness = {}
    with open(RECORD / "results" / "harness_cells.csv", newline="") as h:
        for r in csv.DictReader(h):
            harness[(r["band"], r["metric"], int(r["n"]))] = r
    cells = {}
    with open(SWEEP / "results" / "cells.csv", newline="") as h:
        for r in csv.DictReader(h):
            cells[(r["band"], float(r["x"]), int(r["n"]))] = (
                float(r["recursed"]), float(r["direct"]))
    csv_cert = {}
    with open(RECORD / "results" / "cells_cert.csv", newline="") as h:
        for r in csv.DictReader(h):
            csv_cert[(r["band"], float(r["x"]), int(r["n"]),
                      r["method"])] = r
    exp = model.load_exp()

    out = {"cells": []}
    for band, metric, n in (("A1", "absolute", 2), ("A2", "absolute", 23)):
        h = harness[(band, metric, n)]
        bi = 0 if band == "A1" else 1
        xs, xp = float(h["x_ship"]), float(h["x_prod"])
        d_o = cells[(band, xs, n)][0]
        d_c = cells[(band, xp, n)][1]
        expx = exp[xs][1]

        print("=" * 82)
        print(f"{band} {metric} order {n}")

        mc, fc = product_machine(xp, bi, n)
        mo, fo = ship_machine(xs, expx, n)
        assert fc.d == d_c, "product replay mismatch"
        assert fo.d == d_o, "ship replay mismatch"

        Bc, Dc = ea.product_affine(model, xp, bi, n)
        Bo, Do = ea.ship_affine(model, xs, n)
        eh_lo, eh_hi, _, _ = ea.enclose_exp_neg_half(Fraction(xs))
        eh_mid = (eh_lo + eh_hi) / 2

        rec = dict(cell=f"{band} {metric} {n}", entries=[])
        for role, mach, node, d, xval, B, D, eh_on in (
                ("chosen-product", mc, fc, d_c, xp, Bc, Dc, False),
                ("other-ship", mo, fo, d_o, xs, Bo, Do, True)):
            exp_used = expx if eh_on else 0.0
            S_fixed = B + D * Fraction(exp_used)
            assert node.V == S_fixed, (role, "route mismatch")
            r_fixed = Fraction(d) - node.V          # d - S(x, expx)
            R_mid = B + D * eh_mid
            r_mid = Fraction(d) - R_mid             # d - R at the true eh
            exp_corr = D * (Fraction(exp_used) - eh_mid) if eh_on else Fraction(0)
            assert r_fixed + exp_corr == r_mid

            # DD route at the true exponential
            eh_dd, eh_rem = dd_half(eh_mid)
            maxabs = [0.0]
            if eh_on:
                R_dd = route_dd_ship(xval, n, eh_dd, maxabs)
                n_dd_ops = 8 + 8 * (len(model.ext_piece()[2]) // 2) + 6 * n
            else:
                R_dd = route_dd_product(xval, bi, n, maxabs)
                n_dd_ops = 8 + 4 * len(model.piece(0, bi)[2]) + \
                    2 * len(model.piece(n, bi)[2])
            r_dd = dd_sub(DD(float(d), 0.0), R_dd)
            r_dd_v = r_dd.hi + r_dd.lo
            dd_err = abs(Fraction(r_dd.hi) + Fraction(r_dd.lo) - r_mid)
            dd_bound = (Fraction(n_dd_ops) * 8 * Fraction(U) ** 2
                        * Fraction(max(maxabs[0], 1.0)) * 4
                        + eh_rem * abs(D))
            dmin = min((c["abs_delta_div"] for c in mach.div_calls
                        if c["abs_delta_div"] is not None), default=None)
            rec["entries"].append(dict(
                role=role, x=xval, d=d,
                ops=mach.count, n_ops=sum(mach.count.values()),
                sum_abs_loc=float(mach.sum_abs_loc),
                max_abs_loc=float(mach.max_abs_loc),
                sum_half_ulp=mach.sum_half_ulp,
                r_fixed=float(r_fixed), r_mid=float(r_mid),
                exp_corr=float(exp_corr),
                cancellation=float(mach.sum_abs_loc / abs(r_fixed))
                if r_fixed else None,
                r_dd=float(r_dd_v), dd_err=float(dd_err),
                dd_bound=float(dd_bound),
                dd_bound_over_err=(float(dd_bound / dd_err) if dd_err else None),
                max_intermediate=maxabs[0],
                div_calls=mach.div_calls,
                div_identity_all=all(c["exact_identity"]
                                     for c in mach.div_calls),
                min_abs_delta_div=dmin,
            ))
            print(f"  {role} at x={xval!r}")
            print(f"    nodes {sum(mach.count.values())} {mach.count}")
            print(f"    sum |local residual| = {float(mach.sum_abs_loc):.6e}; "
                  f"net r_fixed = {float(r_fixed):.6e}; cancellation = "
                  f"{float(mach.sum_abs_loc/abs(r_fixed)):.3g}")
            print(f"    sum half-ulp (one-sided, sign-known) = "
                  f"{mach.sum_half_ulp:.6e}")
            print(f"    refined divisions {len(mach.div_calls)}; S18.4 sign "
                  f"identity all: {all(c['exact_identity'] for c in mach.div_calls)}; "
                  f"min |delta_div| = "
                  f"{(dmin if dmin is not None else float('nan')):.3e}")
            print(f"    DD route: |r_dd - r_exact| = {float(dd_err):.3e}; "
                  f"conservative DD bound {float(dd_bound):.3e}; margin "
                  f"{float(dd_bound/dd_err) if dd_err else float('inf'):.1f}x")
            print(f"    exp correction D*(expx-eh) = {float(exp_corr):.6e}")

        # pair margins from the certified truncation intervals; the ordering
        # test differs by cell:
        #   A1 both-negative: |e_c| < |e_o| <=> e_c > e_o, so the threshold is
        #     the signed separation t_c - t_o;
        #   A2 mixed sign:    |e_c| < |e_o| <=> e_c + e_o < 0, so the threshold
        #     is the truncation magnitude sum t_c + t_o.
        fc_box = ea.enclose_boys(n, Fraction(xp))
        fo_box = ea.enclose_boys(n, Fraction(xs))
        Fc = (fc_box[0], fc_box[1])
        Fo = (fo_box[0], fo_box[1])
        Eh = (eh_lo, eh_hi)
        Rc = (Bc, Bc)
        Ro = ea.iadd((Bo, Bo), ea.iscale(Eh, Do))
        Tc = ea.isub(Rc, Fc)
        To = ea.isub(Ro, Fo)
        Ro_mid = Bo + Do * eh_mid
        r_c_pair = Fraction(d_c) - Bc
        r_o_pair = Fraction(d_o) - Ro_mid
        if band == "A1":
            t_test = ea.isub(Tc, To)
            corr_test = r_c_pair - r_o_pair       # r_c - r_o
            test = "signed separation (both errors negative)"
        else:
            t_test = ea.iadd(Tc, To)
            corr_test = r_c_pair + r_o_pair       # r_c + r_o
            test = "magnitude sum (mixed sign)"
        threshold = float(max(abs(t_test[0]), abs(t_test[1])))
        corr_v = float(abs(corr_test))
        margin = corr_v - threshold
        C_c = float(csv_cert[(band, xp, n, "product")]["C"])
        C_o = float(csv_cert[(band, xs, n, "ship")]["C"])
        rec["truncation_gap"] = threshold
        rec["corr_exact"] = float(corr_test)
        rec["corr_over_gap"] = corr_v / threshold if threshold else None
        rec["ordering_test"] = test
        rec["margin"] = margin
        rec["record_C_sum"] = C_c + C_o
        rec["record_C_over_margin"] = (C_c + C_o) / margin if margin else None
        rec["half_record_C_over_margin"] = ((C_c + C_o) / 2) / margin \
            if margin else None
        sum_half = (rec["entries"][0]["sum_half_ulp"]
                    + rec["entries"][1]["sum_half_ulp"])
        rec["eft_half_ulp_pair_band"] = sum_half
        rec["eft_half_ulp_over_margin"] = sum_half / margin if margin else None
        print(f"  ordering test: {test}")
        print(f"  truncation threshold {threshold:.6e}; exact rounding "
              f"correction {corr_v:.6e}; correction/threshold = "
              f"{(corr_v/threshold if threshold else float('nan')):.6g}")
        print(f"  net margin {margin:.6e}")
        print(f"  record sign-blind C_c+C_o = {C_c+C_o:.6e} "
              f"({(C_c+C_o)/margin if margin else float('nan'):.3g}x margin); "
              f"half of it (sign-known, magnitude-blind) = "
              f"{(C_c+C_o)/2:.6e} "
              f"({((C_c+C_o)/2)/margin if margin else float('nan'):.3g}x)")
        print(f"  EFT one-sided per-node band (sum half-ulps) = "
              f"{sum_half:.6e} ({sum_half/margin if margin else float('nan'):.3g}x)")
        out["cells"].append(rec)

    with open(HERE / "eft_compose.json", "w") as fh:
        json.dump(out, fh, indent=1)
    print("=" * 82)
    print("wrote eft_compose.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
