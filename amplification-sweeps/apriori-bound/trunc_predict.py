"""What an a priori *evaluation* of the truncation alone could predict.

The ODE identity gives each route's truncation error exactly, without a
reference.  This script asks how far that a priori information goes: if one
knew the exact truncation errors and nothing about the rounding, would the
measured ordering follow?  It also measures the truncation gap against the
delivered gap at the committed cells.
"""

import csv
from collections import Counter

from mpmath import mp, mpf

mp.dps = 50

ROOT = (r"C:/Users/myaml/dev/papers-a-v5/papers/boys-crossover/sweeps/"
        "amplification_form")


def main():
    pred = {}
    with open(f"{ROOT}/results/prediction.csv", newline="") as handle:
        for r in csv.DictReader(handle):
            pred[(r["band"], r["metric"], int(r["order"]))] = r
    exact = {}
    with open("exact_amplification_form.csv", newline="") as handle:
        for r in csv.DictReader(handle):
            exact[(r["band"], int(r["n"]), r["metric"], r["method"])] = r
    order = {}
    with open(f"{ROOT}/results/per_order.csv", newline="") as handle:
        for r in csv.DictReader(handle):
            order[(r["band"], int(r["n"]), r["method"])] = r

    committed = [k for k, r in pred.items() if r["state"] in ("ship", "product")]
    agree = Counter()
    gaps = []
    rows = []
    for band, metric, n in sorted(committed):
        es = exact[(band, n, metric, "ship")]
        ep = exact[(band, n, metric, "product")]
        ts, tp = mpf(es["route_error"]), mpf(ep["route_error"])
        ds = mpf(order[(band, n, "ship")]["worst_abs" if metric == "absolute"
                                           else "worst_rel"])
        dp = mpf(order[(band, n, "product")]["worst_abs" if metric == "absolute"
                                            else "worst_rel"])
        delivered_winner = "ship" if ds < dp else "product" if dp < ds else "tie"
        trunc_winner = ("ship" if abs(ts) < abs(tp)
                        else "product" if abs(tp) < abs(ts) else "tie")
        agree[trunc_winner == delivered_winner] += 1
        gap_del = abs(ds - dp)
        gap_tr = abs(abs(ts) - abs(tp))
        # How much of the delivered separation is explained by truncation?
        if gap_del > 0:
            gaps.append(float(gap_tr / gap_del))
        rows.append(dict(band=band, metric=metric, n=n,
                         measured=pred[(band, metric, n)]["measured_winner"],
                         delivered_winner=delivered_winner,
                         trunc_winner=trunc_winner,
                         ts=mp.nstr(ts, 4), tp=mp.nstr(tp, 4),
                         ds=mp.nstr(ds, 4), dp=mp.nstr(dp, 4),
                         gap_tr=mp.nstr(gap_tr, 3), gap_del=mp.nstr(gap_del, 3)))
    print(f"committed comparisons: {len(committed)}")
    print("truncation-only predictor agrees with delivered winner:",
          dict(agree))
    ratios = sorted(gaps)
    if ratios:
        print("truncation gap / delivered gap: min %.3g median %.3g max %.3g"
              % (ratios[0], ratios[len(ratios) // 2], ratios[-1]))
    for r in rows:
        print("%(band)s %(metric)-8s n=%(n)2d measured=%(measured)-7s "
              "delivered=%(delivered_winner)-7s trunc=%(trunc_winner)-7s "
              "ts=%(ts)-10s tp=%(tp)-10s ds=%(ds)-9s dp=%(dp)-9s "
              "gap_tr=%(gap_tr)-9s gap_del=%(gap_del)s" % r)


if __name__ == "__main__":
    main()
