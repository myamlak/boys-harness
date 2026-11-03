"""Necessary condition for any sound sign-blind charge to order a commitment.

Reads the deposited sweep records and the exact-truncation CSV. Writes stdout.
"""
import csv
from collections import Counter
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
    cells[(r["band"], float(r["x"]), int(r["n"]))] = (
        float(r["recursed"]), float(r["direct"]))

comm = sorted(k for k, r in pred.items() if r["state"] in ("ship", "product"))
print("committed", len(comm), Counter(pred[k]["state"] for k in comm))
rows = []
for band, metric, n in comm:
    out = {}
    for method, valcol in (("ship", 0), ("product", 1)):
        t, x = exact[(band, n, metric, method)]
        f = ref[(band, x)][n]
        v = cells[(band, x, n)][valcol]
        e = mpf(v) - f                       # signed absolute delivered error
        r = e - t                            # signed arithmetic rounding
        out[method] = dict(x=x, t=t, e=e, r=r, f=f)
    # the model's claimed loser: the entry it says is larger
    state = pred[(band, metric, n)]["state"]
    loser = "product" if state == "ship" else "ship"   # model says the other is smaller
    winner = "ship" if loser == "product" else "product"
    L, W = out[loser], out[winner]
    rows.append(dict(band=band, metric=metric, n=n, state=state,
                     loser=loser, tw=W["t"], tl=L["t"], rw=W["r"], rl=L["r"],
                     ew=W["e"], el=L["e"], xw=W["x"], xl=L["x"]))

def q(vals, p):
    v = sorted(vals); return v[min(len(v)-1, max(0, int(round(p*(len(v)-1)))))]

ratio = [float(abs(r["rl"] / r["tl"])) for r in rows if r["tl"] != 0]
print("loser |rounding|/|truncation|: n=%d min %.4g med %.4g q90 %.4g max %.4g"
      % (len(ratio), min(ratio), q(ratio,.5), q(ratio,.9), max(ratio)))
print("count |r_l| < |t_l| (necessary for a one-signed loser interval):",
      sum(1 for x in ratio if x < 1.0))
for r in rows:
    if r["tl"] != 0 and abs(r["rl"]) < abs(r["tl"]):
        print("  candidate", r["band"], r["metric"], "n=%d" % r["n"], r["state"],
              "loser", r["loser"], "t_l=%.3e r_l=%.3e ratio=%.3f" %
              (float(r["tl"]), float(r["rl"]), float(abs(r["rl"]/r["tl"]))))
# winner's interval would also need to separate
print("loser t sign vs e sign mismatches:",
      sum(1 for r in rows if (r["tl"]>0) != (r["el"]>0)))
# per state
print("loser entry among committed:", Counter(r["loser"] for r in rows))
