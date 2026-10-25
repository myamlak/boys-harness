"""Rounding dominance, with the metrics kept on one scale.

At the absolute cells: delivered absolute error over the evaluated route
truncation error.  At the paper's 67 committed comparisons: the losing
entry's delivered error over its own truncation error (the relative rows
divided by F_n at that entry's argument before the ratio).  The raw
exact_<sweep>.log ratio line pools absolute and relative rows and is not
dimensionally valid for the relative ones; this script is.
"""

import csv
from collections import Counter

from mpmath import mp, mpf

mp.dps = 60

SWEEPS = {
    "amplification_form": r"C:/Users/myaml/dev/papers-a-v5/papers/boys-crossover/sweeps/amplification_form",
    "product_amplification": r"C:/Users/myaml/dev/papers-a-v5/papers/boys-crossover/sweeps/product_amplification",
}


def quantiles(vals):
    vals = sorted(vals)
    n = len(vals)
    return dict(n=n, min=vals[0], q10=vals[n // 10], median=vals[n // 2],
                q90=vals[9 * n // 10], max=vals[-1])


def main():
    for name, root in SWEEPS.items():
        order = {}
        with open(f"{root}/results/per_order.csv", newline="") as handle:
            for r in csv.DictReader(handle):
                order[(r["band"], int(r["n"]), r["method"])] = r
        exact = {}
        at_x = {}
        with open(f"exact_{name}.csv", newline="") as handle:
            for r in csv.DictReader(handle):
                key = (r["band"], int(r["n"]), r["metric"], r["method"])
                exact[key] = mpf(r["route_error"])
                at_x[key] = float(r["x"])
        ref = {}
        with open(f"{root}/data/reference.csv", newline="") as handle:
            for r in csv.DictReader(handle):
                ref[(r["band"], float(r["x"]))] = [mpf(r[c])
                                                   for c in list(r)[2:]]
        vals = []
        dominance = Counter()
        for band in ("A1", "A2"):
            for n in range(33):
                for method in ("ship", "product"):
                    delivered = mpf(order[(band, n, method)]["worst_abs"])
                    trunc = exact[(band, n, "absolute", method)]
                    if trunc == 0:
                        continue
                    q = float(abs(delivered / trunc))
                    vals.append(q)
                    dominance["rounding>10x" if q > 10
                              else "truncation>10x" if q < 0.1
                              else "mixed"] += 1
        print(name, "absolute cells", quantiles(vals), dict(dominance))

        pred = {}
        try:
            with open(f"{root}/results/prediction.csv", newline="") as handle:
                for r in csv.DictReader(handle):
                    pred[(r["band"], r["metric"], int(r["order"]))] = r
        except FileNotFoundError:
            continue
        committed = sorted(k for k, r in pred.items()
                           if r["state"] in ("ship", "product"))
        loser_ratio = []
        over2 = 0
        for band, metric, n in committed:
            ts = exact[(band, n, metric, "ship")]
            tp = exact[(band, n, metric, "product")]
            if metric == "relative":
                ts /= ref[(band, at_x[(band, n, metric, "ship")])][n]
                tp /= ref[(band, at_x[(band, n, metric, "product")])][n]
            da = "worst_abs" if metric == "absolute" else "worst_rel"
            ds = mpf(order[(band, n, "ship")][da])
            dp = mpf(order[(band, n, "product")][da])
            if ds > dp:
                ratio = float(ds / abs(ts)) if ts else float("inf")
            else:
                ratio = float(dp / abs(tp)) if tp else float("inf")
            loser_ratio.append(ratio)
            if ratio > 2:
                over2 += 1
        print(name, "committed losers delivered/trunc",
              quantiles(loser_ratio), ">2x", over2, "of", len(loser_ratio))


if __name__ == "__main__":
    main()
