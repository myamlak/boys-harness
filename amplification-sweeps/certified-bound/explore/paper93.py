"""The paper's section 5.6 statistics at the model's own decision cells.

For each of the 132 comparisons, takes the model's band-max arguments
(prediction.csv's predicted_ship_at_x / predicted_product_at_x), computes
the exact-real truncation with the fast route, the paper-convention
arithmetic charge and the certified radius, and reports

  - whether those arguments differ from per_order.csv's at_x cells,
  - the exact-truncation ordering (16/49 or otherwise),
  - the combined charge / gap needed factor under both charges,
  - the loser's rounding-to-truncation ratio at these cells; the paper's
    9.6 is a per_order-cell quantity and is checked by harnesscells.py.
"""
import argparse
import csv
import math
import os
import statistics as st
import sys
from pathlib import Path

from mpmath import mp, mpf

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "code"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import exactroute as er
import model
import record_paths
import harnesscells as hc

mp.dps = 70


def med(vals):
    vals = [v for v in vals if v is not None and math.isfinite(v)]
    return st.median(vals) if vals else None


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    record_paths.add_kernel(parser)
    record_paths.add_inputs(parser, model_record_arg=True, deposit_arg=True)
    args = parser.parse_args(argv)
    model.configure(kernel=args.kernel, exp_bits=args.exp_bits)
    root = record_paths.sweep(args.sweep)
    model_record = record_paths.model_record(args.model_record)
    dep = record_paths.deposit(args.deposit)
    with record_paths.deposit_on_path(dep):
        import bounds
    hc.bounds = bounds

    ref = er.load_reference(root)
    exp = model.load_exp()
    cells = {}
    with open(root / "results" / "cells.csv", newline="") as h:
        for r in csv.DictReader(h):
            cells[(r["band"], float(r["x"]), int(r["n"]))] = (
                float(r["recursed"]), float(r["direct"]))
    pred = {}
    with open(model_record / "prediction_record.csv", newline="") as h:
        for r in csv.DictReader(h):
            pred[(r["band"], r["metric"], int(r["order"]))] = r
    order = {}
    with open(root / "results" / "per_order.csv", newline="") as h:
        for r in csv.DictReader(h):
            order[(r["band"], int(r["n"]), r["method"])] = r

    differ = 0
    rows = []
    for (band, metric, n), pr in sorted(pred.items()):
        xs = float(pr["predicted_ship_at_x"])
        xp = float(pr["predicted_product_at_x"])
        key = "at_x_abs" if metric == "absolute" else "at_x_rel"
        if (float(order[(band, n, "ship")][key]) != xs
                or float(order[(band, n, "product")][key]) != xp):
            differ += 1
        sv = er.ship_values(xs, mp.exp(-mpf(xs)) / 2)
        pv = er.product_values(xp)
        fs = ref[(band, xs)][n]
        fp = ref[(band, xp)][n]
        ts = sv[n] - fs
        tp = pv[n] - fp
        cs = mpf(model.ship_batch(xs, exp[xs][1])[n]["radius"])
        cp = mpf(model.product_batch(xp, 0 if band == "A1" else 1)[n]["radius"])
        ps = hc.paper_charge(band, n, "ship", xs, ref)
        pp = hc.paper_charge(band, n, "product", xp, ref)
        rs = abs(mpf(cells[(band, xs, n)][0]) - sv[n])
        rp = abs(mpf(cells[(band, xp, n)][1]) - pv[n])
        ts, tp, cs, cp, ps, pp, fs, fp, rs, rp = (
            mpf(v) for v in (ts, tp, cs, cp, ps, pp, fs, fp, rs, rp))
        if metric == "relative":
            ts, tp, cs, cp, ps, pp = ts / fs, tp / fp, cs / fs, cp / fp, ps / fs, pp / fp
            rs, rp = rs / fs, rp / fp
        state = pr["state"]
        trunc = "ship" if abs(ts) < abs(tp) else "product"
        if state not in ("ship", "product"):
            rows.append(dict(band=band, metric=metric, n=n, state=state,
                             trunc=trunc, t_l=None))
            continue
        loser = "product" if state == "ship" else "ship"
        tl, tw = (abs(ts), abs(tp)) if loser == "ship" else (abs(tp), abs(ts))
        cl, cw = (cs, cp) if loser == "ship" else (cp, cs)
        pl, pw = (ps, pp) if loser == "ship" else (pp, ps)
        rl, rw = (rs, rp) if loser == "ship" else (rp, rs)
        gap = tl - tw
        rows.append(dict(band=band, metric=metric, n=n, state=state,
                         trunc=trunc, tl=tl, tw=tw, cl=cl, cw=cw, pl=pl, pw=pw,
                         gap=gap,
                         need_paper=(pl + pw) / gap if gap > 0 else mp.inf,
                         need_cert=(cl + cw) / gap if gap > 0 else mp.inf,
                         ratio_paper=pl / tl, ratio_cert=cl / tl,
                         round_ratio=rl / tl))

    committed = [r for r in rows if r["state"] in ("ship", "product")]
    agree = [r for r in committed if r["trunc"] == r["state"]]
    opp = [r for r in committed if r["trunc"] != r["state"]]
    print("cell sets differ at", differ, "of 132 comparisons")
    print("committed", len(committed), "trunc agree", len(agree),
          "opposite", len(opp))
    print("median loser |r|/t over all 65 at the model cells:",
          med([float(r["round_ratio"]) for r in committed]))
    print("median paper charge/loser t:", med([float(r["ratio_paper"]) for r in committed]))
    print("median cert charge/loser t:", med([float(r["ratio_cert"]) for r in committed]))
    print("over the", len(agree), "agreeing cells:")
    print("  median paper combined/gap:", med([float(r["need_paper"]) for r in agree]))
    print("  median cert combined/gap:", med([float(r["need_cert"]) for r in agree]))
    print("  min paper", min(float(r["need_paper"]) for r in agree if math.isfinite(float(r["need_paper"]))))


if __name__ == "__main__":
    sys.exit(main())
