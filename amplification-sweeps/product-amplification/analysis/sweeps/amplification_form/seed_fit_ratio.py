"""Distribution of the seed-against-fit ratio over the upward cells.

For every cell (x, n >= 1) at or above its dispatch threshold, the order-0
seed value at the same argument is compared with the direct entry's value at
that order.  The ratio reported here is the relative error of the seed value
divided by the relative error of the direct entry:

    r = (|rec(0) - F_0(x)| / F_0(x)) / (|dir(n) - F_n(x)| / F_n(x))

Inputs are the v5.1.1 deposit (cells.csv, reference.csv) and its constants.
Output: results/seed_fit_ratio.csv, one row per band and order, plus a pooled
row per band.  Quartiles use linear interpolation between order statistics
(the type-7 definition).
"""

import csv
import json
import os
import sys
from collections import defaultdict

import mpmath as mp

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
CELLS = os.path.join(ROOT, "product_amplification", "results", "cells.csv")
REFERENCE = os.path.join(ROOT, "product_amplification", "data", "reference.csv")
PROVENANCE = os.path.join(ROOT, "product_amplification", "results", "provenance.json")
OUT = os.path.join(HERE, "results", "seed_fit_ratio.csv")

mp.mp.dps = 60


def quantile(values, p):
    if not values:
        return float("nan")
    values = sorted(values)
    if len(values) == 1:
        return values[0]
    pos = p * (len(values) - 1)
    low = int(pos)
    high = min(low + 1, len(values) - 1)
    frac = pos - low
    return values[low] + (values[high] - values[low]) * frac


def main():
    constants = json.load(open(PROVENANCE))["constants"]
    theta = constants["kTierThresholds"]

    reference = {}
    with open(REFERENCE, newline="") as handle:
        for row in csv.DictReader(handle):
            key = (row["band"], float(row["x"]))
            reference[key] = [mp.mpf(row["f{}".format(n)]) for n in range(33)]

    seed_value = {}
    cells = []
    with open(CELLS, newline="") as handle:
        for row in csv.DictReader(handle):
            band = row["band"]
            x = float(row["x"])
            n = int(row["n"])
            rec = mp.mpf(row["recursed"])
            direct = mp.mpf(row["direct"])
            if n == 0:
                seed_value[(band, x)] = rec
            cells.append((band, x, n, rec, direct))

    ratios = defaultdict(list)
    infinite = defaultdict(int)
    for band, x, n, rec, direct in cells:
        if n == 0 or x < theta[n]:
            continue
        values = reference[(band, x)]
        f0 = values[0]
        fn = values[n]
        seed_rel = abs(seed_value[(band, x)] - f0) / f0
        fit_rel = abs(direct - fn) / fn
        if fit_rel == 0:
            infinite[band] += 1
            continue
        r = float(seed_rel / fit_rel)
        ratios[(band, n)].append(r)
        ratios[(band, "all")].append(r)

    rows = []
    for band in ("A1", "A2"):
        for n in list(range(1, 33)) + ["all"]:
            values = ratios.get((band, n), [])
            rows.append({
                "band": band,
                "order": n,
                "cells": len(values),
                "q25": "{:.10g}".format(quantile(values, 0.25)),
                "median": "{:.10g}".format(quantile(values, 0.50)),
                "q75": "{:.10g}".format(quantile(values, 0.75)),
                "infinite": infinite[band] if n == "all" else "",
            })
    with open(OUT, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    for row in rows:
        if row["order"] == "all":
            print(row)
    total = sum(row["cells"] for row in rows if row["order"] == "all")
    print("cells", total)


if __name__ == "__main__":
    sys.exit(main())
