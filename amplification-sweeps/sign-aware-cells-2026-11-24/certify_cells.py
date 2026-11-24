"""Certified sign-aware certificate at the two chosen reversal cells.

Cell 1: A1 absolute order 2   (both delivered errors negative)
Cell 2: A2 absolute order 23  (chosen positive, other negative)

At each cell the chosen entry is the direct/product lane at x_prod and the
other entry is the recursed/ship lane at x_ship (harness_cells.csv state).

The certificate is built from the arithmetic, not from the deposited error
columns:

  * the delivered double is replayed with IEEE binary64 operations in the
    kernel's order (mirror of code/model.py, checked against model.py and
    against results/cells.csv bit for bit);
  * the exact route value R = B + D*eh is computed in exact rational
    arithmetic from the stored binary64 coefficients and arguments, with
    eh = exp(-x)/2 carried symbolically;
  * exp(-x)/2 and F_n(x) are enclosed by exact-rational alternating series
    with the alternating-remainder bound;
  * the delivered error is e = d - F_n, and its decomposition
    e = t + r, t = R - F_n, r = d - R is reported with exact intervals.

The deposited reference and the measured error columns are used only in the
checks: each certificate interval is compared with the measured value, and
the record's own signed-residual test is reproduced.

Run:  python -I certify_cells.py [--target 1e-60]
Writes: cell_certificate.csv, cell_certificate.json (in this directory).
"""
import argparse
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
import exact_arith as ea

RECORD = REPO / "amplification-sweeps" / "certified-bound"
SWEEP = Path(os.environ.get("CERT_SWEEP", REPO / "amplification-sweeps" / "amplification-form"))
KERNEL = Path(os.environ.get("BOYS_ROOT", REPO.parent / "boys"))

sys.path.insert(0, str(RECORD / "code"))
import model                      # noqa: E402
import exactroute as er           # noqa: E402
from mpmath import mp, mpf        # noqa: E402

mp.dps = 90

CELLS = [("A1", "absolute", 2), ("A2", "absolute", 23)]


def load_csv(path):
    with open(path, newline="") as h:
        return list(csv.DictReader(h))


def m(fr):
    return mpf(fr.numerator) / fr.denominator


def s(fr, digits=25):
    return mp.nstr(m(fr), digits)


def sgn_interval(I):
    if I[0] > 0:
        return 1
    if I[1] < 0:
        return -1
    return 0


def mag_interval(I):
    if I[0] > 0:
        return I
    if I[1] < 0:
        return (-I[1], -I[0])
    return (Fraction(0), max(-I[0], I[1]))


def width(I):
    return I[1] - I[0]


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", default="1e-62",
                    help="tail target for the alternating-series enclosures")
    args = ap.parse_args(argv)
    target = Fraction(args.target)

    model.configure(kernel=KERNEL, exp_bits=RECORD / "exp" / "exp_bits.txt")
    harness = {(r["band"], r["metric"], int(r["n"])): r
               for r in load_csv(RECORD / "results" / "harness_cells.csv")}
    cert = {(r["band"], float(r["x"]), int(r["n"]), r["method"]): r
            for r in load_csv(RECORD / "results" / "cells_cert.csv")}
    cells = {}
    for r in load_csv(SWEEP / "results" / "cells.csv"):
        cells[(r["band"], float(r["x"]), int(r["n"]))] = (
            float(r["recursed"]), float(r["direct"]))
    ref = er.load_reference(SWEEP)
    exp = model.load_exp()

    rows = []
    checks = {}
    for band, metric, n in CELLS:
        h = harness[(band, metric, n)]
        bi = 0 if band == "A1" else 1
        xs, xp = float(h["x_ship"]), float(h["x_prod"])
        d_c = cells[(band, xp, n)][1]          # chosen: product at x_prod
        d_o = cells[(band, xs, n)][0]          # other: ship at x_ship
        expx = exp[xs][1]

        # ---- replay and bit checks -------------------------------------
        rep_c = ea.product_value(model, xp, bi, n)
        rep_o = ea.ship_value(model, xs, expx, n)
        mod_c = model.product_batch(xp, bi)[n]["value"]
        mod_o = model.ship_batch(xs, expx)[n]["value"]
        checks[f"{band}-{n}-prod-replay"] = rep_c == d_c
        checks[f"{band}-{n}-ship-replay"] = rep_o == d_o
        checks[f"{band}-{n}-prod-model"] = mod_c == d_c
        checks[f"{band}-{n}-ship-model"] = mod_o == d_o

        # ---- exact route, affine in eh = exp(-x)/2 ---------------------
        Bc, Dc = ea.product_affine(model, xp, bi, n)
        Bo, Do = ea.ship_affine(model, xs, n)
        fc = ea.enclose_boys(n, Fraction(xp), target)
        fo = ea.enclose_boys(n, Fraction(xs), target)
        eh = ea.enclose_exp_neg_half(Fraction(xs), target)
        Fc = (fc[0], fc[1])
        Fo = (fo[0], fo[1])
        Eh = (eh[0], eh[1])
        Rc = ea.iadd((Bc, Bc), ea.iscale(Eh, Dc))
        Ro = ea.iadd((Bo, Bo), ea.iscale(Eh, Do))
        Tc = ea.isub(Rc, Fc)
        To = ea.isub(Ro, Fo)
        Rc_c = ea.isub((Fraction(d_c), Fraction(d_c)), Rc)
        Ro_r = ea.isub((Fraction(d_o), Fraction(d_o)), Ro)
        Ec = ea.isub((Fraction(d_c), Fraction(d_c)), Fc)
        Eo = ea.isub((Fraction(d_o), Fraction(d_o)), Fo)

        # ---- measured check values (70-digit route vs 60-digit ref) ----
        Fc_ref = ref[(band, xp)][n]
        Fo_ref = ref[(band, xs)][n]
        route_c = er.product_values(xp)[n]
        route_o = er.ship_values(xs, mp.exp(-mpf(xs)) / 2)[n]
        e_c_meas = mpf(d_c) - Fc_ref
        e_o_meas = mpf(d_o) - Fo_ref
        t_c_meas = route_c - Fc_ref
        t_o_meas = route_o - Fo_ref
        r_c_meas = mpf(d_c) - route_c
        r_o_meas = mpf(d_o) - route_o
        # our route at the midpoint of the eh enclosure, and the exact-route
        # values at 70 digits, are two independent evaluations of the same R
        eh_mid = (Eh[0] + Eh[1]) / 2
        Rc_mid = Bc + Dc * eh_mid
        Ro_mid = Bo + Do * eh_mid
        checks[f"{band}-{n}-route-c-70d"] = abs(m(Rc_mid) - route_c) < mpf("1e-40")
        checks[f"{band}-{n}-route-o-70d"] = abs(m(Ro_mid) - route_o) < mpf("1e-40")

        # containment of the measured columns in the certified intervals,
        # allowing only the reference's 60-digit rounding and the record's
        # 1e-45 tier margin
        def inside(I, v, tol=mpf("1e-45")):
            return m(I[0]) - tol <= v <= m(I[1]) + tol

        checks[f"{band}-{n}-e-c-contains"] = inside(Ec, e_c_meas)
        checks[f"{band}-{n}-e-o-contains"] = inside(Eo, e_o_meas)
        checks[f"{band}-{n}-t-c-contains"] = inside(Tc, t_c_meas)
        checks[f"{band}-{n}-t-o-contains"] = inside(To, t_o_meas)
        checks[f"{band}-{n}-r-c-contains"] = inside(Rc_c, r_c_meas)
        checks[f"{band}-{n}-r-o-contains"] = inside(Ro_r, r_o_meas)

        # ---- record columns -------------------------------------------
        C_c = float(cert[(band, xp, n, "product")]["C"])
        C_o = float(cert[(band, xs, n, "ship")]["C"])
        eps_c = float(cert[(band, xp, n, "product")]["eps"])
        eps_o = float(cert[(band, xs, n, "ship")]["eps"])
        h_sum = mpf(C_c) + mpf(C_o) + mpf(eps_c) + mpf(eps_o)

        # ---- derived quantities ---------------------------------------
        s_t = ea.isub(Tc, To)                  # truncation-centre separation
        corr = ea.isub(Rc_c, Ro_r)             # sign-aware rounding correction
        e_sep = ea.isub(Ec, Eo)                # signed separation of errors
        e_sum = ea.iadd(Ec, Eo)                # sum of errors
        mag_c = mag_interval(Ec)
        mag_o = mag_interval(Eo)

        signed_above = Ec[0] > Eo[1]           # e_c > e_o certified
        signed_below = Ec[1] < Eo[0]           # e_c < e_o certified
        magnitude = mag_c[1] < mag_o[0]        # |e_c| < |e_o| certified
        # same-sign signed certificate (the paper's T test)
        same_sign = ((Ec[1] < 0 and Eo[1] < 0 and Ec[0] > Eo[1])
                     or (Ec[0] > 0 and Eo[0] > 0 and Ec[1] < Eo[0]))
        sgn_e_c = sgn_interval(Ec)
        sgn_e_o = sgn_interval(Eo)
        sgn_sep = sgn_interval(e_sep)
        sgn_sum = sgn_interval(e_sum)
        sgn_corr = sgn_interval(corr)
        sgn_tsep = sgn_interval(s_t)
        t_sum = ea.iadd(Tc, To)
        sgn_tsum = sgn_interval(t_sum)

        # margins for the two orderings
        # signed: need e_c - e_o > 0; worst-case margin
        signed_margin = Ec[0] - Eo[1]
        # magnitude at A1 (both negative): |e_c| < |e_o| <=> e_c > e_o
        mag_margin = mag_o[0] - mag_c[1]
        # rounding correction over the truncation-centre separation
        corr_minus_gap = corr[0] - abs(s_t[0] if s_t[0] > 0 else s_t[1])
        # for the mixed cell the magnitude centre is t_c + t_o
        corr_over_tsum = -(t_sum[1]) - corr[1] if t_sum[1] > 0 else None

        # required width from the record's sign-blind charge
        gap = abs(float(s_t[0] if s_t[0] > 0 else s_t[1]))
        signblind_factor = float(h_sum / mpf(gap)) if gap else None
        # the exact correction's size and the required reduction factor
        corr_mid = (corr[0] + corr[1]) / 2
        cert_width = width(e_sep)
        record_width = h_sum

        rec = dict(
            cell=f"{band} {metric} {n}",
            band=band, metric=metric, n=n,
            x_ship=xs, x_prod=xp, chosen="product", other="ship",
            d_chosen=d_c, d_other=d_o, expx=expx,
            B_chosen=s(Bc), D_chosen=s(Dc), B_other=s(Bo), D_other=s(Do),
            eh_lo=s(Eh[0]), eh_hi=s(Eh[1]), eh_width=float(width(Eh)),
            F_chosen_lo=s(Fc[0]), F_chosen_hi=s(Fc[1]),
            F_other_lo=s(Fo[0]), F_other_hi=s(Fo[1]),
            t_chosen_lo=s(Tc[0]), t_chosen_hi=s(Tc[1]),
            t_other_lo=s(To[0]), t_other_hi=s(To[1]),
            r_chosen_lo=s(Rc_c[0]), r_chosen_hi=s(Rc_c[1]),
            r_other_lo=s(Ro_r[0]), r_other_hi=s(Ro_r[1]),
            e_chosen_lo=s(Ec[0]), e_chosen_hi=s(Ec[1]),
            e_other_lo=s(Eo[0]), e_other_hi=s(Eo[1]),
            t_sep_lo=s(s_t[0]), t_sep_hi=s(s_t[1]), t_sep_sign=sgn_tsep,
            t_sum_lo=s(t_sum[0]), t_sum_hi=s(t_sum[1]), t_sum_sign=sgn_tsum,
            corr_lo=s(corr[0]), corr_hi=s(corr[1]),
            corr_width=float(width(corr)), corr_sign=sgn_corr,
            e_sep_lo=s(e_sep[0]), e_sep_hi=s(e_sep[1]),
            e_sep_width=float(width(e_sep)), e_sep_sign=sgn_sep,
            e_sum_lo=s(e_sum[0]), e_sum_hi=s(e_sum[1]), e_sum_sign=sgn_sum,
            sign_e_chosen=sgn_e_c, sign_e_other=sgn_e_o,
            signed_above=bool(signed_above), signed_below=bool(signed_below),
            same_sign=bool(same_sign), magnitude=bool(magnitude),
            signed_margin=float(signed_margin),
            magnitude_margin=float(mag_margin),
            corr_minus_gap=float(corr_minus_gap),
            C_chosen=C_c, C_other=C_o,
            eps_chosen=eps_c, eps_other=eps_o, h_sum=float(h_sum),
            signblind_factor=signblind_factor,
            required_cert_width_vs_record=float(record_width / mpf(cert_width)),
            measured_e_chosen=str(e_c_meas), measured_e_other=str(e_o_meas),
            measured_t_chosen=str(t_c_meas), measured_t_other=str(t_o_meas),
            measured_r_chosen=str(r_c_meas), measured_r_other=str(r_o_meas),
            measured_e_sep=str(e_c_meas - e_o_meas),
            measured_corr=str(r_c_meas - r_o_meas),
            measured_corr_plus=str(r_c_meas + r_o_meas),
            reference_F_chosen=str(Fc_ref), reference_F_other=str(Fo_ref),
        )
        # truncation-centre direction and whether the certificate overturns it
        if s_t[0] > 0:
            trunc_dir = "chosen above other"
        elif s_t[1] < 0:
            trunc_dir = "chosen below other"
        else:
            trunc_dir = "straddles"
        rec["trunc_dir"] = trunc_dir
        # A1: signed order.  A2: magnitude order (sign of |t_chosen|-|t_other|
        # is the sign of t_chosen + t_other for opposite signs).
        if metric == "absolute" and band == "A1":
            rec["trunc_centre_ordering"] = (
                "chosen magnitudes larger (t_c < t_o)" if s_t[1] < 0
                else ("chosen magnitudes smaller" if s_t[0] > 0 else "straddles"))
            rec["overturns_truncation_centre"] = bool(
                s_t[1] < 0 and Ec[0] > Eo[1])
        else:
            rec["trunc_centre_ordering"] = (
                "chosen magnitude larger (t_c + t_o > 0)" if t_sum[0] > 0
                else ("chosen magnitude smaller" if t_sum[1] < 0
                      else "straddles"))
            rec["overturns_truncation_centre"] = bool(
                t_sum[0] > 0 and magnitude)
        rows.append(rec)

        # ---- print ------------------------------------------------------
        print("=" * 82)
        print(f"{band} {metric} order {n}: chosen=product at x={xp!r}, "
              f"other=ship at x={xs!r}")
        print(f"  delivered chosen {ea.Fraction(d_c)} = {float.hex(d_c)}")
        print(f"  delivered other  {ea.Fraction(d_o)} = {float.hex(d_o)}")
        print(f"  expx = {float.hex(expx)}; exact eh enclosed in "
              f"[{s(Eh[0], 18)}, {s(Eh[1], 18)}], width {float(width(Eh)):.2e}")
        print(f"  route chosen R = {s(Bc, 18)} (D=0)")
        print(f"  route other  R = {s(Bo, 18)} + ({s(Do, 18)}) * eh")
        print(f"  F_{n} chosen in [{s(Fc[0], 18)}, {s(Fc[1], 18)}] "
              f"(width {float(width(Fc)):.2e})")
        print(f"  F_{n} other  in [{s(Fo[0], 18)}, {s(Fo[1], 18)}] "
              f"(width {float(width(Fo)):.2e})")
        print("  truncation t:")
        print(f"    chosen [{s(Tc[0], 18)}, {s(Tc[1], 18)}]")
        print(f"    other  [{s(To[0], 18)}, {s(To[1], 18)}]")
        print("  rounding r = d - R:")
        print(f"    chosen [{s(Rc_c[0], 18)}, {s(Rc_c[1], 18)}] "
              f"(width {float(width(Rc_c)):.2e})")
        print(f"    other  [{s(Ro_r[0], 18)}, {s(Ro_r[1], 18)}] "
              f"(width {float(width(Ro_r)):.2e})")
        print("  delivered error e = d - F:")
        print(f"    chosen [{s(Ec[0], 18)}, {s(Ec[1], 18)}] sign "
              f"{'+' if sgn_e_c > 0 else '-' if sgn_e_c < 0 else '0'}")
        print(f"    other  [{s(Eo[0], 18)}, {s(Eo[1], 18)}] sign "
              f"{'+' if sgn_e_o > 0 else '-' if sgn_e_o < 0 else '0'}")
        print("  certified roundings bounds:")
        print(f"    r_c - r_o in [{s(corr[0], 18)}, {s(corr[1], 18)}] "
              f"sign {'+' if sgn_corr > 0 else '-' if sgn_corr < 0 else '0'}, "
              f"width {float(width(corr)):.2e}")
        print(f"    r_c + r_o in [{s(ea.iadd(Rc_c, Ro_r)[0], 18)}, "
              f"{s(ea.iadd(Rc_c, Ro_r)[1], 18)}]")
        print(f"  truncation-centre separation t_c - t_o in "
              f"[{s(s_t[0], 18)}, {s(s_t[1], 18)}]")
        print(f"  truncation-centre magnitude sum t_c + t_o in "
              f"[{s(t_sum[0], 18)}, {s(t_sum[1], 18)}]")
        print(f"  certified e_c - e_o in [{s(e_sep[0], 18)}, "
              f"{s(e_sep[1], 18)}] sign "
              f"{'+' if sgn_sep > 0 else '-' if sgn_sep < 0 else '0'}")
        print(f"  certified e_c + e_o in [{s(e_sum[0], 18)}, "
              f"{s(e_sum[1], 18)}] sign "
              f"{'+' if sgn_sum > 0 else '-' if sgn_sum < 0 else '0'}")
        print(f"  CERTIFIED signed above: {signed_above}; "
              f"signed below: {signed_below}; same-sign: {same_sign}; "
              f"magnitude |e_c|<|e_o|: {magnitude}")
        print(f"  signed margin {float(signed_margin):.6e}; magnitude margin "
              f"{float(mag_margin):.6e}; correction over gap "
              f"{float(corr_minus_gap):.6e}")
        print(f"  record sign-blind h_sum {float(h_sum):.6e}; "
              f"truncation gap {gap:.6e}; factor {signblind_factor:.6g}")
        print(f"  certified interval width {float(cert_width):.3e}; record "
              f"width / certified width {float(record_width / mpf(cert_width)):.3e}")
        print(f"  truncation centre: {rec['trunc_centre_ordering']}; "
              f"overturned: {rec['overturns_truncation_centre']}")
        print(f"  measured check: e_c={mp.nstr(e_c_meas, 12)} "
              f"e_o={mp.nstr(e_o_meas, 12)} "
              f"(in certified intervals: {checks[f'{band}-{n}-e-c-contains']}, "
              f"{checks[f'{band}-{n}-e-o-contains']})")

    print("=" * 82)
    bad = [k for k, v in checks.items() if not v]
    print(f"checks: {len(checks) - len(bad)}/{len(checks)} pass")
    for k in bad:
        print("  FAIL", k)

    keys = list(rows[0])
    with open(HERE / "cell_certificate.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(keys)
        for r in rows:
            w.writerow([r[k] for k in keys])
    with open(HERE / "cell_certificate.json", "w") as fh:
        json.dump(dict(
            target=str(target), checks=checks, cells=rows,
            checks_pass=len(checks) - len(bad), checks_total=len(checks),
        ), fh, indent=1)
    print("wrote cell_certificate.csv, cell_certificate.json")
    return 0 if not bad else 1


if __name__ == "__main__":
    sys.exit(main())
