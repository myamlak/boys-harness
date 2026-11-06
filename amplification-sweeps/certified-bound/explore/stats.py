"""Exploratory: size the existing charge against delivered and truncation.

Reads the deposited harness modules read-only; writes only to stdout.
"""
import argparse
import csv
import sys
from collections import Counter
from pathlib import Path

from mpmath import mp, mpf

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "code"))
import record_paths

mp.dps = 60

_parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
record_paths.add_kernel(_parser)
record_paths.add_inputs(_parser, model_record_arg=True, deposit_arg=True)
_args = _parser.parse_args()
ROOT = record_paths.sweep(_args.sweep)
MODEL = record_paths.model_record(_args.model_record)
HERE = record_paths.deposit(_args.deposit)
with record_paths.deposit_on_path(HERE):
    import bounds


def q(vals, p):
    v = sorted(vals)
    return v[min(len(v) - 1, max(0, int(round(p * (len(v) - 1)))))]


order = {}
for r in csv.DictReader(open(ROOT / "results" / "per_order.csv")):
    order[(r["band"], int(r["n"]), r["method"])] = dict(
        worst_abs=mpf(r["worst_abs"]), at_x_abs=float(r["at_x_abs"]),
        worst_rel=mpf(r["worst_rel"]), at_x_rel=float(r["at_x_rel"]))

pred = {}
for r in csv.DictReader(open(MODEL / "prediction.csv")):
    pred[(r["band"], r["metric"], int(r["order"]))] = r

exact = {}
for r in csv.DictReader(open(HERE / "exact_amplification_form.csv")):
    exact[(r["band"], int(r["n"]), r["metric"], r["method"])] = (
        mpf(r["route_error"]), float(r["x"]))

ref = {}
for r in csv.DictReader(open(ROOT / "data" / "reference.csv")):
    ref[(r["band"], float(r["x"]))] = [mpf(r[c]) for c in list(r)[2:]]


def lookup_ref(band, x):
    key = (band, x)
    if key in ref:
        return ref[key]
    best = min((abs(kx - x), kx) for b, kx in ref if b == band)
    return ref[(band, best[1])]


modes = Counter()
rows = []
for band in ("A1", "A2"):
    for metric in ("absolute", "relative"):
        for n in range(33):
            k = (band, metric, n)
            st = pred[k]["state"]
            xa = "at_x_abs" if metric == "absolute" else "at_x_rel"
            da = "worst_abs" if metric == "absolute" else "worst_rel"
            xs = order[(band, n, "ship")][xa]
            xp = order[(band, n, "product")][xa]
            ds = order[(band, n, "ship")][da]
            dp = order[(band, n, "product")][da]
            ts = exact[(band, n, metric, "ship")][0]
            tp = exact[(band, n, metric, "product")][0]
            bs = bounds.entry_bounds(band, n, "ship", xs, signed_method="sampled",
                                     samples=400, charge_kind="paper")
            bp = bounds.entry_bounds(band, n, "product", xp, signed_method="sampled",
                                     samples=400, charge_kind="paper")
            cs, cp = bs["charge"], bp["charge"]
            if metric == "relative":
                fs = lookup_ref(band, xs)[n]
                fp = lookup_ref(band, xp)[n]
                ts, tp = ts / fs, tp / fp
                cs, cp = cs / fs, cp / fp
            modes[(st, bs["mode"], bp["mode"])] += 1
            rows.append(dict(band=band, metric=metric, n=n, state=st,
                             xs=xs, xp=xp, ds=ds, dp=dp, ts=ts, tp=tp,
                             cs=cs, cp=cp, ms=bs["mode"], mp=bp["mode"]))

comm = [r for r in rows if r["state"] in ("ship", "product")]
print("comparisons", len(rows), "committed", len(comm), Counter(r["state"] for r in comm))
print("modes all", modes)
print("modes committed", Counter((r["ms"], r["mp"]) for r in comm))

for meth, dk, tk, ck in (("ship", "ds", "ts", "cs"), ("product", "dp", "tp", "cp")):
    ratio_dt = [float(abs(r[dk] / r[tk])) if r[tk] != 0 else 0.0 for r in comm]
    ratio_cd = [float(abs(r[ck] / r[dk])) if r[dk] != 0 else 0.0 for r in comm]
    print(meth, "delivered/|trunc| min med max",
          min(ratio_dt), q(ratio_dt, .5), max(ratio_dt),
          "| charge/delivered min med max",
          min(ratio_cd), q(ratio_cd, .5), max(ratio_cd))

loser_ratio = []
one_signed = 0
need = []
unfixable = 0
for r in comm:
    loser = "ship" if r["ds"] > r["dp"] else "product"
    d = r["ds"] if loser == "ship" else r["dp"]
    t = r["ts"] if loser == "ship" else r["tp"]
    c = r["cs"] if loser == "ship" else r["cp"]
    loser_ratio.append(float(abs(d / t)) if t != 0 else 0.0)
    if abs(t) > c:
        one_signed += 1
    tw = r["ts"] if loser == "product" else r["tp"]
    cw = r["cs"] if loser == "product" else r["cp"]
    if abs(t) <= abs(tw):
        unfixable += 1
        continue
    slack = abs(t) - abs(tw)
    need.append(float((c + cw) / slack))
print("loser delivered/|trunc| min med q90 max",
      min(loser_ratio), q(loser_ratio, .5), q(loser_ratio, .9), max(loser_ratio))
print("loser interval one-signed under its own paper charge:", one_signed, "of", len(comm))
print("truncation cannot order (|t_l|<=|t_w|):", unfixable)
print("charge reduction needed min med max", min(need), q(need, .5), max(need))

# same-fit cells: identical coefficients, different x? check truncation equality
same = [r for r in comm if r["ms"] == "same-fit" and r["mp"] == "fit"]
print("committed same-fit", len(same),
      "equal trunc", sum(1 for r in same if r["ts"] == r["tp"]))
