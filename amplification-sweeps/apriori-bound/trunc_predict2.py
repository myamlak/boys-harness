"""Truncation-only a priori prediction, metric-aware.

For the relative metric the exact truncation error is divided by F_n at the
same argument (from the sweep's 60-digit reference, used here only to put the
two metrics on one scale; the predictor itself uses only the identity).
"""

import csv
from collections import Counter

from mpmath import mp, mpf

mp.dps = 50

ROOT = (r"C:/Users/myaml/dev/papers-a-v5/papers/boys-crossover/sweeps/"
        "amplification_form")


def main():
    ref = {}
    with open(f"{ROOT}/data/reference.csv", newline="") as handle:
        for r in csv.DictReader(handle):
            ref[(r["band"], float(r["x"]))] = [mpf(r[col])
                                               for col in list(r)[2:]]
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
    out = {}
    for metric in ("absolute", "relative"):
        agree = Counter()
        gaps = []
        for band, met, n in sorted(committed):
            if met != metric:
                continue
            es = exact[(band, n, met, "ship")]
            ep = exact[(band, n, met, "product")]
            ts, tp = mpf(es["route_error"]), mpf(ep["route_error"])
            if metric == "relative":
                fs = ref[(band, float(es["x"]))][n]
                fp = ref[(band, float(ep["x"]))][n]
                ts, tp = ts / fs, tp / fp
            ds = mpf(order[(band, n, "ship")]["worst_abs" if metric == "absolute"
                                               else "worst_rel"])
            dp = mpf(order[(band, n, "product")]["worst_abs" if metric == "absolute"
                                                else "worst_rel"])
            delivered_winner = ("ship" if ds < dp else "product" if dp < ds else "tie")
            trunc_winner = ("ship" if abs(ts) < abs(tp)
                            else "product" if abs(tp) < abs(ts) else "tie")
            agree[trunc_winner == delivered_winner] += 1
            gap_del = abs(ds - dp)
            gap_tr = abs(abs(ts) - abs(tp))
            if gap_del > 0:
                gaps.append(float(gap_tr / gap_del))
        gaps.sort()
        out[metric] = dict(
            agree=dict(agree), n=sum(agree.values()),
            gap_min=gaps[0], gap_median=gaps[len(gaps) // 2], gap_max=gaps[-1],
            agreement_rate=agree[True] / max(1, sum(agree.values())))
        print(metric, out[metric])


if __name__ == "__main__":
    main()
