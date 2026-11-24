"""Reference consistency and comparison with the record's signed-residual test.

For the two chosen cells this script
  * recomputes the certificate intervals exactly,
  * evaluates the record's measured quantities at 100 digits
    (delivered - reference, exact route - reference, delivered - exact route),
  * reports the difference between the certificate midpoint and each measured
    quantity (the reference's own rounding),
  * reproduces the record's sign-blind signed-residual test at the two cells
    (it intersects, as the 0/16 deposit result says),
  * reports the certified sign, magnitude and truncation-centre outcome.

The deposited reference enters only in this check, never in the certificate.

Run:  python -I check_reference.py
Writes: reference_check.json
"""
import csv
import json
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
from mpmath import mp, mpf                   # noqa: E402

mp.dps = 100


def m(fr):
    return mpf(fr.numerator) / fr.denominator


def mid(I):
    return (m(I[0]) + m(I[1])) / 2


def sgn(I):
    if I[0] > 0:
        return 1
    if I[1] < 0:
        return -1
    return 0


def mag(I):
    if I[0] > 0:
        return I
    if I[1] < 0:
        return (-I[1], -I[0])
    return (Fraction(0), max(-I[0], I[1]))


def main():
    model.configure(kernel=KERNEL, exp_bits=RECORD / "exp" / "exp_bits.txt")
    harness = {}
    with open(RECORD / "results" / "harness_cells.csv", newline="") as h:
        for r in csv.DictReader(h):
            harness[(r["band"], r["metric"], int(r["n"]))] = r
    cert = {}
    with open(RECORD / "results" / "cells_cert.csv", newline="") as h:
        for r in csv.DictReader(h):
            cert[(r["band"], float(r["x"]), int(r["n"]), r["method"])] = r
    cells = {}
    with open(SWEEP / "results" / "cells.csv", newline="") as h:
        for r in csv.DictReader(h):
            cells[(r["band"], float(r["x"]), int(r["n"]))] = (
                float(r["recursed"]), float(r["direct"]))
    ref = er.load_reference(SWEEP)
    exp = model.load_exp()

    out = {"cells": [], "note": (
        "diffs are certificate midpoint minus the record's measured value; "
        "the reference enters only here")}
    for band, metric, n in (("A1", "absolute", 2), ("A2", "absolute", 23)):
        h = harness[(band, metric, n)]
        bi = 0 if band == "A1" else 1
        xs, xp = float(h["x_ship"]), float(h["x_prod"])
        d_c = cells[(band, xp, n)][1]
        d_o = cells[(band, xs, n)][0]
        expx = exp[xs][1]

        Bc, Dc = ea.product_affine(model, xp, bi, n)
        Bo, Do = ea.ship_affine(model, xs, n)
        eh = ea.enclose_exp_neg_half(Fraction(xs))
        Fc = ea.enclose_boys(n, Fraction(xp))
        Fo = ea.enclose_boys(n, Fraction(xs))
        Eh = (eh[0], eh[1])
        Rc = (Bc, Bc)
        Ro = ea.iadd((Bo, Bo), ea.iscale(Eh, Do))
        Tc = ea.isub(Rc, (Fc[0], Fc[1]))
        To = ea.isub(Ro, (Fo[0], Fo[1]))
        Rc_r = ea.isub((Fraction(d_c), Fraction(d_c)), Rc)
        Ro_r = ea.isub((Fraction(d_o), Fraction(d_o)), Ro)
        Ec = ea.isub((Fraction(d_c), Fraction(d_c)), (Fc[0], Fc[1]))
        Eo = ea.isub((Fraction(d_o), Fraction(d_o)), (Fo[0], Fo[1]))

        # record measured quantities at 100 digits
        Fc_ref = ref[(band, xp)][n]
        Fo_ref = ref[(band, xs)][n]
        route_c = er.product_values(xp)[n]
        route_o = er.ship_values(xs, mp.exp(-mpf(xs)) / 2)[n]
        meas = dict(
            e_c=mpf(d_c) - Fc_ref, e_o=mpf(d_o) - Fo_ref,
            t_c=route_c - Fc_ref, t_o=route_o - Fo_ref,
            r_c=mpf(d_c) - route_c, r_o=mpf(d_o) - route_o,
        )
        boxes = dict(e_c=Ec, e_o=Eo, t_c=Tc, t_o=To, r_c=Rc_r, r_o=Ro_r)
        diffs = {k: mid(boxes[k]) - meas[k] for k in boxes}
        # containment: the rigorous interval plus the record's 1e-45 tier margin
        contains = {k: (m(boxes[k][0]) - mpf("1e-45") <= meas[k]
                        <= m(boxes[k][1]) + mpf("1e-45")) for k in boxes}

        # record sign-blind signed-residual test at this cell
        C_c = mpf(cert[(band, xp, n, "product")]["C"])
        C_o = mpf(cert[(band, xs, n, "ship")]["C"])
        eps_c = mpf(cert[(band, xp, n, "product")]["eps"])
        eps_o = mpf(cert[(band, xs, n, "ship")]["eps"])
        h_sum = C_c + C_o + eps_c + eps_o
        lo_c, hi_c = meas["t_c"] - h_sum, meas["t_c"] + h_sum
        lo_o, hi_o = meas["t_o"] - h_sum, meas["t_o"] + h_sum
        record_above = bool(lo_c > hi_o)
        record_below = bool(hi_c < lo_o)
        record_same_sign = bool(
            (hi_c < 0 and hi_o < 0 and lo_c > hi_o)
            or (lo_c > 0 and lo_o > 0 and hi_c < lo_o))
        def mag_mp(lo, hi):
            if lo > 0:
                return lo, hi
            if hi < 0:
                return -hi, -lo
            return mpf(0), max(-lo, hi)

        record_mag = bool(mag_mp(lo_c, hi_c)[1] < mag_mp(lo_o, hi_o)[0])
        s_t = ea.isub(Tc, To)
        corr = ea.isub(Rc_r, Ro_r)
        e_sep = ea.isub(Ec, Eo)
        e_sum = ea.iadd(Ec, Eo)
        mc, mo = mag(Ec), mag(Eo)
        certified = dict(
            sign_e_c=sgn(Ec), sign_e_o=sgn(Eo),
            sign_e_sep=sgn(e_sep), sign_e_sum=sgn(e_sum),
            signed_above=bool(Ec[0] > Eo[1]),
            same_sign=bool((Ec[1] < 0 and Eo[1] < 0 and Ec[0] > Eo[1])
                           or (Ec[0] > 0 and Eo[0] > 0 and Ec[1] < Eo[0])),
            magnitude=bool(mc[1] < mo[0]),
            sign_corr_sep=sgn(corr),
        )
        # the paper-framework combination: the record's truncation centre
        # (70-digit route minus reference, tier margin eps) plus the certified
        # rounding correction
        eps_pair = eps_c + eps_o
        if band == "A1":
            base = meas["t_c"] - meas["t_o"]
            r_iv = corr
            paper_lo = base + m(r_iv[0]) - eps_pair
            paper_hi = base + m(r_iv[1]) + eps_pair
            paper_ok = bool(paper_lo > 0)
        else:
            base = meas["t_c"] + meas["t_o"]
            r_iv = ea.iadd(Rc_r, Ro_r)
            paper_lo = base + m(r_iv[0]) - eps_pair
            paper_hi = base + m(r_iv[1]) + eps_pair
            paper_ok = bool(paper_hi < 0)
        rec_paper = dict(
            base=str(base),
            rounding_lo=str(m(r_iv[0])), rounding_hi=str(m(r_iv[1])),
            eps_pair=float(eps_pair),
            combined_lo=str(paper_lo), combined_hi=str(paper_hi),
            certifies=paper_ok,
        )
        rec = dict(
            cell=f"{band} {metric} {n}", band=band, n=n,
            measured={k: str(v) for k, v in meas.items()},
            certificate_midpoint_minus_measured={k: str(v) for k, v in diffs.items()},
            max_abs_diff=float(max(abs(v) for v in diffs.values())),
            contains_with_tier_margin=contains,
            record_interval=dict(h_sum=float(h_sum),
                                 truncation_sep=float(abs(s_t[0] if s_t[0] > 0
                                                          else s_t[1])),
                                 h_sum_over_sep=float(
                                     h_sum / abs(s_t[0] if s_t[0] > 0
                                                 else s_t[1]))),
            record_test=dict(above=record_above, below=record_below,
                             same_sign=record_same_sign,
                             magnitude=record_mag),
            paper_framework=rec_paper,
            certified=certified,
            certified_correction=dict(
                r_c_minus_r_o_lo=str(m(corr[0])),
                r_c_minus_r_o_hi=str(m(corr[1])),
                r_c_plus_r_o_lo=str(m(ea.iadd(Rc_r, Ro_r)[0])),
                r_c_plus_r_o_hi=str(m(ea.iadd(Rc_r, Ro_r)[1])),
                e_c_minus_e_o_lo=str(m(e_sep[0])),
                e_c_minus_e_o_hi=str(m(e_sep[1])),
                e_c_plus_e_o_lo=str(m(e_sum[0])),
                e_c_plus_e_o_hi=str(m(e_sum[1])),
            ),
        )
        out["cells"].append(rec)

        print("=" * 78)
        print(f"{band} {metric} order {n}")
        print("  reference check (certificate midpoint - record value):")
        for k in ("e_c", "e_o", "t_c", "t_o", "r_c", "r_o"):
            print(f"    {k}: {mp.nstr(diffs[k], 5)}")
        print(f"  max |diff| = {float(rec['max_abs_diff']):.3e} "
              f"(60-digit reference rounding)")
        print(f"  containment with the record's 1e-45 tier margin: "
              f"{all(contains.values())}")
        print(f"  record sign-blind test: above={record_above} "
              f"below={record_below} same_sign={record_same_sign} "
              f"magnitude={record_mag} (intersects, as the 0/16 deposit says)")
        print(f"  certified: signs e_c={certified['sign_e_c']:+d} "
              f"e_o={certified['sign_e_o']:+d}; "
              f"e_c-e_o sign {certified['sign_e_sep']:+d}; "
              f"e_c+e_o sign {certified['sign_e_sum']:+d}; "
              f"signed_above={certified['signed_above']} "
              f"same_sign={certified['same_sign']} "
              f"magnitude={certified['magnitude']}")
        print(f"  paper framework (record t + certified rounding): "
              f"combined in [{mp.nstr(paper_lo, 8)}, {mp.nstr(paper_hi, 8)}], "
              f"certifies={paper_ok}")

    with open(HERE / "reference_check.json", "w") as fh:
        json.dump(out, fh, indent=1)
    print("=" * 78)
    print("wrote reference_check.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
