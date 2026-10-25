"""The best a stored per-piece truncation bound could do, with the charge.

The ODE identity gives each route's truncation error from the stored
coefficients alone, to quadrature accuracy -- an oracle for any per-piece
truncation bound the kernel could store.  This script asks whether that
oracle plus the accumulated-arithmetic charge can certify the paper's 67
committed comparisons: loser lower bound max(0, |trunc_loser| - C_loser)
against winner upper bound |trunc_winner| + C_winner.
"""

import csv
from collections import Counter

from mpmath import mp, mpf, inf

import bounds

mp.dps = 60

ROOT = (r"C:/Users/myaml/dev/papers-a-v5/papers/boys-crossover/sweeps/"
        "amplification_form")


def main():
    pred = {}
    with open(f"{ROOT}/results/prediction.csv", newline="") as handle:
        for r in csv.DictReader(handle):
            pred[(r["band"], r["metric"], int(r["order"]))] = r
    order = {}
    with open(f"{ROOT}/results/per_order.csv", newline="") as handle:
        for r in csv.DictReader(handle):
            order[(r["band"], int(r["n"]), r["method"])] = r
    exact = {}
    at_x = {}
    with open("exact_amplification_form.csv", newline="") as handle:
        for r in csv.DictReader(handle):
            exact[(r["band"], int(r["n"]), r["metric"], r["method"])] = mpf(
                r["route_error"])
            at_x[(r["band"], int(r["n"]), r["metric"], r["method"])] = float(
                r["x"])
    ref = {}
    with open(f"{ROOT}/data/reference.csv", newline="") as handle:
        for r in csv.DictReader(handle):
            ref[(r["band"], float(r["x"]))] = [mpf(r[c]) for c in list(r)[2:]]

    committed = sorted(k for k, r in pred.items()
                       if r["state"] in ("ship", "product"))
    out = {}
    for kind in ("paper", "worst"):
        counts = Counter()
        ratios = []
        trunc_over_charge = 0
        for band, metric, n in committed:
            da = "worst_abs" if metric == "absolute" else "worst_rel"
            xs = float(order[(band, n, "ship")]["at_x_abs" if metric == "absolute"
                                              else "at_x_rel"])
            xp = float(order[(band, n, "product")]["at_x_abs" if metric == "absolute"
                                                else "at_x_rel"])
            ts = exact[(band, n, metric, "ship")]
            tp = exact[(band, n, metric, "product")]
            cs = bounds.entry_bounds(band, n, "ship", xs,
                                     charge_kind=kind)["charge"]
            cp = bounds.entry_bounds(band, n, "product", xp,
                                     charge_kind=kind)["charge"]
            if metric == "relative":
                fs = ref[(band, at_x[(band, n, metric, "ship")])][n]
                fp = ref[(band, at_x[(band, n, metric, "product")])][n]
                ts, tp, cs, cp = ts / fs, tp / fp, cs / fs, cp / fp
            ms = bounds.magnitude_interval(ts - cs, ts + cs)
            mpp = bounds.magnitude_interval(tp - cp, tp + cp)
            cert = ""
            if ms[1] < mpp[0]:
                cert = "ship"
            elif mpp[1] < ms[0]:
                cert = "product"
            counts[cert or "abstain"] += 1
            # how short the measured loser's exact truncation is of its charge
            ds = mpf(order[(band, n, "ship")][da])
            dp = mpf(order[(band, n, "product")][da])
            loser = "ship" if ds > dp else "product"
            if loser == "ship":
                ratios.append(float(abs(ts) / cs) if cs else inf)
                if cs and abs(ts) > cs:
                    trunc_over_charge += 1
            else:
                ratios.append(float(abs(tp) / cp) if cp else inf)
                if cp and abs(tp) > cp:
                    trunc_over_charge += 1
        ratios.sort()
        out[kind] = dict(counts=dict(counts),
                         ratio=dict(min=ratios[0],
                                    median=ratios[len(ratios) // 2],
                                    max=ratios[-1]),
                         trunc_over_charge=trunc_over_charge)
        print(kind, out[kind])


if __name__ == "__main__":
    main()
