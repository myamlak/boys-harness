"""Regression test for the downward chain at x = 2**-1074.

At the least positive subnormal argument the order-32 seed is F_32(x) ~ 1/65
and the first downward product x*F_32 has exact value about 7.601e-326,
below half the least positive subnormal (2**-1075), so the computed product
is zero.  Before the repair the model charged a zero radius for that
multiplication; the exact result is nonzero and the operation-level
enclosure was therefore false.  The repaired rule carries the absolute
allowance 2**-1071 whenever the result or the radius bracket is below the
normal range.

This test records the exact and computed intermediate values and the radii
at every affected step, checks |computed - exact| <= radius at each product
node and at the end of the chain, and checks that a normal product (x=1e-300)
does not receive the allowance.

Exit status 0 on pass, 1 on failure.  Run with --kernel / BOYS_ROOT.
"""
import argparse
import csv
import math
import os
import sys

from mpmath import mp, mpf

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import exactroute as er
import model
import record_paths

mp.dps = 70
X_UNDERFLOW = math.ldexp(1.0, -1074)   # 2**-1074, the least positive subnormal
X_CONTROL = 1e-12                      # a normal bracket at the same chain


def check(condition, message):
    if not condition:
        print("FAIL:", message)
        raise SystemExit(1)


def exact_chain(x, seed_exact, exp_half_exact, n=33):
    """Exact-real downward values at x, from the exact seed."""
    f = seed_exact
    out = {32: f}
    for l in range(31, -1, -1):
        f = (mpf(x) * f + exp_half_exact) / (l + mpf("0.5"))
        out[l] = f
    return out


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    record_paths.add_kernel(parser)
    record_paths.add_inputs(parser, exp_bits_arg=True)
    parser.add_argument("--out", metavar="DIR",
                        help="optional run directory for underflow_test.csv")
    args = parser.parse_args(argv)
    model.configure(kernel=args.kernel, exp_bits=args.exp_bits)
    out = record_paths.out_dir(args.out) if args.out else None

    model.table()
    exp = model.load_exp()
    check(X_UNDERFLOW in exp,
          f"the exponential probe has no entry at x={X_UNDERFLOW!r}")
    expx = exp[X_UNDERFLOW][1]
    half_exact = mp.exp(-mpf(X_UNDERFLOW)) / 2

    print("x =", repr(X_UNDERFLOW), "= 2**-1074")
    print("least positive subnormal =", repr(math.ldexp(1.0, -1074)))
    print("2**-1075 as a double =", repr(2.0 ** -1075),
          "(not representable; the rule carries 2**-1071 =",
          repr(2.0 ** -1071), ")")
    print("probe exp(-x)/2 =", repr(expx),
          " exact e^-x/2 =", mp.nstr(half_exact, 25))

    # ---- the seed ---------------------------------------------------------
    seed_v, seed_e = model.eval_fit_split(32, 0, X_UNDERFLOW)
    seed_exact = er.fit_value(32, 0, X_UNDERFLOW)
    seed_err = abs(seed_exact - mpf(seed_v))
    print("\norder-32 seed")
    print("  computed", repr(seed_v))
    print("  exact   ", mp.nstr(seed_exact, 25))
    print("  |error| ", mp.nstr(seed_err, 6), " radius", repr(seed_e))
    check(seed_err <= mpf(seed_e), "the seed radius does not cover the seed")

    # ---- the affected product and the downward chain ----------------------
    rows = []
    f, ef = seed_v, seed_e
    f_exact = seed_exact
    exact_values = exact_chain(X_UNDERFLOW, seed_exact, half_exact)
    first_affected = None
    for l in range(31, -1, -1):
        p, ep = model.mul(X_UNDERFLOW, 0.0, f, ef)
        exact_prod = mpf(X_UNDERFLOW) * f_exact
        err = abs(mpf(p) - exact_prod)
        rows.append(dict(l=l, computed_product=p,
                         exact_product=mp.nstr(exact_prod, 17),
                         exact_product_float=float(exact_prod),
                         error=mp.nstr(err, 17), radius=repr(ep),
                         covered=err <= mpf(ep)))
        if first_affected is None and p == 0.0 and exact_prod > 0:
            first_affected = (l, exact_prod, ep, err)
        check(err <= mpf(ep),
              f"order {l}: |computed - exact product| {mp.nstr(err, 6)} "
              f"exceeds radius {ep!r}")
        check(ep > 0.0, f"order {l}: the product radius is zero")
        a, ea = model.add(p, ep, expx, 0.0)
        invd = 1.0 / (float(l) + 0.5)
        f, ef = model.refined_div(a, ea, float(l) + 0.5, invd)
        f_exact = exact_values[l]

    check(first_affected is not None,
          "no underflowed product was found in the downward chain")
    l0, exact0, radius0, err0 = first_affected
    print("\nfirst affected product (order", l0, ")")
    print("  computed x*F =", repr(0.0))
    print("  exact  x*F  =", mp.nstr(exact0, 17), "(about 7.6e-326)")
    print("  |error|     =", mp.nstr(err0, 17))
    print("  radius      =", repr(radius0))
    check(mpf(radius0) >= exact0,
          "the affected product radius does not cover the exact product")

    # ---- the end of the chain --------------------------------------------
    ship = model.ship_batch(X_UNDERFLOW, expx)
    route = er.ship_values(X_UNDERFLOW, half_exact)
    worst = mpf(0)
    for n in range(33):
        err = abs(route[n] - mpf(ship[n]["value"]))
        worst = max(worst, err)
        check(err <= mpf(ship[n]["radius"]),
              f"order {n}: the delivered chain error {mp.nstr(err, 6)} "
              f"exceeds the total radius {ship[n]['radius']!r}")
    print("\nwhole chain: 33 orders, worst |model - exact route| =",
          mp.nstr(worst, 6), " all inside the certified radii")

    # ---- the control: a normal product keeps the first-order radius -------
    seed_v2, seed_e2 = model.eval_fit_split(32, 0, X_CONTROL)
    seed_exact2 = er.fit_value(32, 0, X_CONTROL)
    p2, ep2 = model.mul(X_CONTROL, 0.0, seed_v2, seed_e2)
    err2 = abs(mpf(p2) - mpf(X_CONTROL) * seed_exact2)
    bracket = model.U_ADD * abs(p2) + abs(seed_v2) * 0.0 + abs(X_CONTROL) * seed_e2
    print("\ncontrol x =", X_CONTROL)
    print("  computed product", repr(p2), " nonzero:", p2 != 0.0)
    print("  |error|", mp.nstr(err2, 6), " radius", repr(ep2))
    print("  allowance applied:", model._allowance(p2, bracket) != 0.0)
    check(p2 != 0.0, "the control product underflowed unexpectedly")
    check(err2 <= mpf(ep2), "the control product radius does not cover it")
    check(model._allowance(p2, bracket) == 0.0,
          "the allowance changed a normal-range product radius")

    if out is not None:
        path = out / "underflow_test.csv"
        with open(path, "w", newline="") as h:
            w = csv.writer(h)
            w.writerow(["order_l", "computed_product", "exact_product",
                        "error", "radius", "covered"])
            for r in rows:
                w.writerow([r["l"], repr(r["computed_product"]),
                            r["exact_product"], r["error"], r["radius"],
                            1 if r["covered"] else 0])
        print("\nwrote", path)

    print("\nPASS: the underflow allowance encloses the lost product at",
          "x = 2**-1074")
    return 0


if __name__ == "__main__":
    sys.exit(main())
