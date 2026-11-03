"""The four commitments a sound sign-blind charge could possibly order."""
import csv
import os
import sys

sys.path.insert(0, r"C:/Users/myaml/dev/_a-cert-run/code")
import model

from mpmath import mp, mpf

mp.dps = 60
ROOT = r"C:/Users/myaml/dev/papers-a-v5/papers/boys-crossover/sweeps/amplification_form"
HERE = r"C:/Users/myaml/dev/_harness-v1023/amplification-sweeps/apriori-bound"

pred = {}
for r in csv.DictReader(open(ROOT + "/results/prediction.csv")):
    pred[(r["band"], r["metric"], int(r["order"]))] = r
order = {}
for r in csv.DictReader(open(ROOT + "/results/per_order.csv")):
    order[(r["band"], int(r["n"]), r["method"])] = r
exact = {}
for r in csv.DictReader(open(HERE + "/exact_amplification_form.csv")):
    exact[(r["band"], int(r["n"]), r["metric"], r["method"])] = (
        mpf(r["route_error"]), float(r["x"]))
ref = {}
for r in csv.DictReader(open(ROOT + "/data/reference.csv")):
    ref[(r["band"], float(r["x"]))] = [mpf(r[c]) for c in list(r)[2:]]
cells = {}
for r in csv.DictReader(open(ROOT + "/results/cells.csv")):
    cells[(r["band"], float(r["x"]), int(r["n"]))] = (float(r["recursed"]), float(r["direct"]))
exp = model.load_exp()

for band, metric, n in sorted(k for k, v in pred.items()
                              if v["state"] in ("ship", "product")):
    t_s, xs = exact[(band, n, metric, "ship")]
    t_p, xp = exact[(band, n, metric, "product")]
    f_s = ref[(band, xs)][n]
    f_p = ref[(band, xp)][n]
    if not (abs(t_p) > 0 and abs(mpf(cells[(band, xp, n)][1]) - f_p - t_p)
            < abs(t_p)):
        continue
    ship = model.ship_batch(xs, exp[xs][1])[n]
    prod = model.product_batch(xp, 0 if band == "A1" else 1)[n]
    e_s = mpf(cells[(band, xs, n)][0]) - f_s
    e_p = mpf(cells[(band, xp, n)][1]) - f_p
    print(band, metric, n, "state", pred[(band, metric, n)]["state"])
    print("  ship x=%r t=%s e=%s C=%s mode=%s" %
          (xs, mp.nstr(t_s, 6), mp.nstr(e_s, 6), mp.nstr(ship["radius"], 6), ship["mode"]))
    print("  prod x=%r t=%s e=%s C=%s" %
          (xp, mp.nstr(t_p, 6), mp.nstr(e_p, 6), mp.nstr(prod["radius"], 6)))
    ts = abs(t_s)
    tp = abs(t_p)
    rel_s = (mpf(max(0.0, float(ts - ship["radius"]))), (ts + ship["radius"]) / f_s)
    rel_p = (mpf(max(0.0, float(tp - prod["radius"]))), (tp + prod["radius"]) / f_p)
    print("  relative band: ship [%s, %s] product [%s, %s]" %
          (mp.nstr(rel_s[0] / f_s, 4), mp.nstr(rel_s[1] / f_s, 4),
           mp.nstr(rel_p[0] / f_p, 4), mp.nstr(rel_p[1] / f_p, 4)))
