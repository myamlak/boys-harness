"""Representative committed comparisons with their certified intervals.

Prints the cell-level evidence behind APRIORI-BOUND.md: a compact table of
eight of the paper's 67 committed comparisons (the extremes, the one case
where the loser's truncation interval is one-signed, and cells inside both
dispatch regimes), plus the aggregate loser/winner and gap statistics.
"""

import csv
from collections import Counter

from mpmath import mp, mpf

import bounds

mp.dps = 60

ROOT = (r"C:/Users/myaml/dev/papers-a-v5/papers/boys-crossover/sweeps/"
        "amplification_form")
THETA0 = bounds.THETA0


def main():
    pred = {}
    with open(f"{ROOT}/results/prediction.csv", newline="") as handle:
        for r in csv.DictReader(handle):
            pred[(r["band"], r["metric"], int(r["order"]))] = r
    order = {}
    with open(f"{ROOT}/results/per_order.csv", newline="") as handle:
        for r in csv.DictReader(handle):
            order[(r["band"], int(r["n"]), r["method"])] = r
    committed = sorted(k for k, r in pred.items()
                       if r["state"] in ("ship", "product"))
    rows = []
    below = Counter()
    for band, metric, n in committed:
        xa = "at_x_abs" if metric == "absolute" else "at_x_rel"
        da = "worst_abs" if metric == "absolute" else "worst_rel"
        xs = float(order[(band, n, "ship")][xa])
        xp = float(order[(band, n, "product")][xa])
        ds = mpf(order[(band, n, "ship")][da])
        dp = mpf(order[(band, n, "product")][da])
        bs = bounds.entry_bounds(band, n, "ship", xs)
        bp = bounds.entry_bounds(band, n, "product", xp)
        ms = bounds.magnitude_interval(bs["lo"], bs["hi"])
        mpp = bounds.magnitude_interval(bp["lo"], bp["hi"])
        loser = "ship" if ds > dp else "product"
        win, lose = (ds, dp) if loser == "product" else (dp, ds)
        ratio = lose / win if win else mp.inf
        if xs < THETA0:
            below["ship_cell_below_theta0"] += 1
        if xp < THETA0:
            below["product_cell_below_theta0"] += 1
        rows.append(dict(band=band, metric=metric, n=n, xs=xs, xp=xp,
                         ds=ds, dp=dp, ratio=ratio,
                         los_lo=ms[0] if loser == "ship" else mpp[0],
                         los_hi=ms[1] if loser == "ship" else mpp[1],
                         win_hi=ms[1] if loser == "product" else mpp[1],
                         loser=loser, ms=ms, mp=mpp))
    rows.sort(key=lambda r: r["ratio"])
    print("committed", len(rows), "below", dict(below))
    print("loser/winner delivered: min %.3g (%s %s n=%d) max %.3g (%s %s n=%d)"
          % (rows[0]["ratio"], rows[0]["band"], rows[0]["metric"], rows[0]["n"],
             rows[-1]["ratio"], rows[-1]["band"], rows[-1]["metric"],
             rows[-1]["n"]))
    gaps = sorted(abs(r["ds"] - r["dp"]) for r in rows)
    print("delivered gap: min %.3g median %.3g max %.3g"
          % (gaps[0], gaps[len(gaps) // 2], gaps[-1]))
    one_signed = [r for r in rows if r["los_lo"] > 0]
    print("loser truncation one-signed:", len(one_signed),
          [(r["band"], r["metric"], r["n"]) for r in one_signed])

    picks = [rows[0], rows[-1],
             next(r for r in rows if (r["band"], r["metric"], r["n"])
                  == ("A1", "absolute", 0)),
             next(r for r in rows if (r["band"], r["metric"], r["n"])
                  == ("A1", "absolute", 15)),
             next(r for r in rows if (r["band"], r["metric"], r["n"])
                  == ("A2", "absolute", 32)),
             next(r for r in rows if (r["band"], r["metric"], r["n"])
                  == ("A2", "relative", 32)),
             next(r for r in rows if (r["band"], r["metric"], r["n"])
                  == ("A2", "absolute", 23)),
             next(r for r in rows if (r["band"], r["metric"], r["n"])
                  == ("A1", "relative", 0))]
    print("band metric n x_ship x_prod delivered_ship delivered_prod "
          "loser loser_lo loser_hi winner_hi")
    for r in picks:
        print("%s %s %2d %.6g %.6g %s %s %s %s %s %s"
              % (r["band"], r["metric"], r["n"], r["xs"], r["xp"],
                 mp.nstr(r["ds"], 6), mp.nstr(r["dp"], 6), r["loser"],
                 mp.nstr(r["los_lo"], 6), mp.nstr(r["los_hi"], 6),
                 mp.nstr(r["win_hi"], 6)))


if __name__ == "__main__":
    main()
