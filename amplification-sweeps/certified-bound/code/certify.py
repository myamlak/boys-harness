"""Certified total-error intervals for both Boys entries, and the scoring.

Tier T: the truncation is the exact-real route value minus the 60-digit
reference, with a computed error bound of 1e-45 absolute (70-digit mpmath on
60-digit references).  The deposit's independent identity-based oracle
cross-check agrees to 8e-50 on the non-outlier checked cells; one cell,
recorded in oracle_check.log, resolves at 7.99e-24 and sits inside the tier's
bound.  The arithmetic is the running-error radius of model.py.  The total
interval at a cell is [t - eps - C, t + eps + C].

Outputs, all beneath the run directory selected with --out:
  cells_cert.csv      one row per (band, x, order, method), 55 176 rows
  model_cells.csv     the 132 comparisons at the model's band-max arguments
  band_cert.csv       the 132 comparisons on the population-band intervals
  summary.json        counts, validation, obstruction statistics

No paper or deposit file is written.  The sweep, the model record, the
kernel and the exponential bits are selected with --sweep, --model-record,
--kernel and --exp-bits (record_paths.py).
"""
import argparse
import csv
import json
import math
import os
import sys
from collections import Counter, defaultdict
from pathlib import Path

from mpmath import mp, mpf

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import exactroute as er
import model
import record_paths

mp.dps = 70
EPS_BASE = mpf("1e-45")
F_REL = mpf("1e-55")


def load(root, model_record):
    root = Path(root)
    model_record = Path(model_record)
    ref = er.load_reference(root)
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
    exp = model.load_exp()
    return ref, cells, pred, order, exp


def sub(a, b):
    return a - b


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    record_paths.add_kernel(parser)
    record_paths.add_inputs(parser, model_record_arg=True)
    record_paths.add_out(parser)
    args = parser.parse_args(argv)
    model.configure(kernel=args.kernel, exp_bits=args.exp_bits)
    root = record_paths.sweep(args.sweep)
    model_record = record_paths.model_record(args.model_record)
    out = record_paths.out_dir(args.out)
    theta0 = model.theta0()

    ref, cells, pred, order, exp = load(root, model_record)

    # ---- the per-cell pass ------------------------------------------------
    rec = {}          # (band, x, n, method) -> dict
    band = {}         # (band, metric, n) -> dict of L/U for both entries
    for b in ("A1", "A2"):
        for met in ("absolute", "relative"):
            for n in range(33):
                band[(b, met, n)] = dict(
                    sLo=None, sHi=None, pLo=None, pHi=None,
                    sAt=None, pAt=None, sLoT=None, sHiT=None)

    viol = []
    worst_ratio = mpf(0)
    worst_ratio_at = None
    modes = Counter()
    min_abs_value = None
    max_abs_value = mpf(0)
    max_C = mpf(0)
    n_cells = 0

    for b in ("A1", "A2"):
        bi = 0 if b == "A1" else 1
        xs = sorted(x for bb, x in ref if bb == b)
        for x in xs:
            expx = exp[x][1]
            sv = er.ship_values(x, mp.exp(-mpf(x)) / 2)
            pv = er.product_values(x)
            sC = model.ship_batch(x, expx)
            pC = model.product_batch(x, bi)
            for n in range(33):
                F = ref[(b, x)][n]
                ts = sv[n] - F
                tp = pv[n] - F
                for meth, vals, C in (("ship", sv, sC[n]), ("product", pv, pC[n])):
                    delivered = cells[(b, x, n)][0 if meth == "ship" else 1]
                    route = vals[n]
                    e = mpf(delivered) - F
                    r = mpf(delivered) - route
                    t = ts if meth == "ship" else tp
                    tf = float(t)
                    eps = EPS_BASE + mpf(2) ** -51 * abs(tf)
                    c = float(C["radius"])
                    assert math.isfinite(c) and c >= 0.0, (b, x, n, meth, c)
                    ratio = abs(r) / mpf(c) if c > 0 else mpf(0)
                    if ratio > worst_ratio:
                        worst_ratio, worst_ratio_at = ratio, (b, x, n, meth)
                    inside = abs(e - t) <= eps + mpf(c)
                    if not inside:
                        viol.append((b, x, n, meth, float(abs(r)),
                                     float(mpf(c)), float(t), float(e)))
                    rec[(b, x, n, meth)] = dict(
                        t=t, tf=tf, eps=eps, C=c, F=F, e=e, r=r,
                        mode=C["mode"], c_seed=C.get("c_seed", 0.0),
                        c_rec=C.get("c_rec", 0.0), c_exp=C.get("c_exp", 0.0),
                        inside=inside)
                    modes[C["mode"]] += 1
                    av = abs(mpf(delivered))
                    if av > max_abs_value:
                        max_abs_value = av
                    if av > 0 and (min_abs_value is None or av < min_abs_value):
                        min_abs_value = av
                    max_C = max(max_C, mpf(c))
                    n_cells += 1
                # population band intervals
                if mpf(x) < theta0:
                    continue
                for met in ("absolute", "relative"):
                    Flo = F * (1 - F_REL)
                    Fhi = F * (1 + F_REL)
                    for meth, t, c in (("ship", ts, float(sC[n]["radius"])),
                                       ("product", tp, float(pC[n]["radius"]))):
                        tf = float(t)
                        eps = float(EPS_BASE + mpf(2) ** -51 * abs(tf))
                        lo = max(0.0, abs(tf) - eps - c)
                        hi = abs(tf) + eps + c
                        if met == "relative":
                            lo = lo / float(Fhi)
                            hi = hi / float(Flo)
                        key = (b, met, n)
                        row = band[key]
                        pl, ph, pat = ("sLo", "sHi", "sAt") if meth == "ship" \
                            else ("pLo", "pHi", "pAt")
                        if row[pl] is None or lo > row[pl]:
                            row[pl] = lo
                            row[pat] = x
                        if row[ph] is None or hi > row[ph]:
                            row[ph] = hi

    # ---- write the per-cell table -----------------------------------------
    with open(out / "cells_cert.csv", "w", newline="") as h:
        w = csv.writer(h)
        w.writerow(["band", "x", "n", "method", "t", "eps", "C", "e", "r",
                    "F", "mode", "inside", "c_seed", "c_rec", "c_exp"])
        for (b, x, n, meth), d in sorted(rec.items()):
            w.writerow([b, repr(x), n, meth, mp.nstr(d["t"], 17),
                        mp.nstr(d["eps"], 6), repr(d["C"]),
                        mp.nstr(d["e"], 17), mp.nstr(d["r"], 17),
                        mp.nstr(d["F"], 17), d["mode"],
                        1 if d["inside"] else 0,
                        repr(d["c_seed"]), repr(d["c_rec"]), repr(d["c_exp"])])

    # ---- helpers ----------------------------------------------------------
    def mag(d):
        lo = d["tf"] - float(d["eps"]) - d["C"]
        hi = d["tf"] + float(d["eps"]) + d["C"]
        if lo <= 0 <= hi:
            return 0.0, max(abs(lo), abs(hi)), lo, hi
        return min(abs(lo), abs(hi)), max(abs(lo), abs(hi)), lo, hi

    def relmag(d):
        lo, hi, rlo, rhi = mag(d)
        Flo = float(d["F"]) * (1 - float(F_REL))
        Fhi = float(d["F"]) * (1 + float(F_REL))
        cand = [rlo / Flo, rlo / Fhi, rhi / Flo, rhi / Fhi]
        lo2, hi2 = min(cand), max(cand)
        if lo2 <= 0 <= hi2:
            return 0.0, max(abs(lo2), abs(hi2)), lo2, hi2
        return min(abs(lo2), abs(hi2)), max(abs(lo2), abs(hi2)), lo2, hi2

    def compare(s, p, met):
        if met == "absolute":
            ls, hs, _, _ = mag(s)
            lp, hp, _, _ = mag(p)
        else:
            ls, hs, _, _ = relmag(s)
            lp, hp, _, _ = relmag(p)
        if hs < lp:
            return "ship", ls, hs, lp, hp
        if hp < ls:
            return "product", ls, hs, lp, hp
        return "", ls, hs, lp, hp

    def order_of(s, p, met):
        """Truncation-only ordering, in the metric's own units."""
        if met == "relative":
            a = abs(s["t"]) / s["F"]
            b = abs(p["t"]) / p["F"]
        else:
            a = abs(s["t"])
            b = abs(p["t"])
        return "ship" if a < b else "product"

    # ---- the 132 comparisons at the model's band-max arguments ------------
    model_rows = []
    for (b, met, n), pr in sorted(pred.items()):
        xs = float(pr["predicted_ship_at_x"])
        xp = float(pr["predicted_product_at_x"])
        s = rec[(b, xs, n, "ship")]
        p = rec[(b, xp, n, "product")]
        cert, ls, hs, lp, hp = compare(s, p, met)
        state = pr["state"]
        to = order_of(s, p, met)
        model_rows.append(dict(band=b, metric=met, n=n, state=state,
                               measured=pr["measured_winner"], cert=cert,
                               xs=xs, xp=xp, trunc_order=to,
                               s=s, p=p, ls=ls, hs=hs, lp=lp, hp=hp,
                               predicted_ship=float(pr["predicted_ship"]),
                               predicted_product=float(pr["predicted_product"])))

    # ---- the 132 comparisons on the population band intervals -------------
    band_rows = []
    for (b, met, n), pr in sorted(pred.items()):
        row = band[(b, met, n)]
        ws_lo, ws_hi = row["sLo"], row["sHi"]
        wp_lo, wp_hi = row["pLo"], row["pHi"]
        if ws_hi < wp_lo:
            cert = "ship"
        elif wp_hi < ws_lo:
            cert = "product"
        else:
            cert = ""
        # truncation-only band ordering (C = 0, eps ignored), in the metric's
        # own units: max_x |t(x)| or max_x |t(x)|/F(x)
        def trunc_at(x, meth):
            d = rec[(b, x, n, meth)]
            return abs(d["t"]) / d["F"] if met == "relative" else abs(d["t"])
        ts = max(trunc_at(x, "ship") for x in
                 {xx for bb, xx in ref if bb == b and xx >= theta0})
        tp = max(trunc_at(x, "product") for x in
                 {xx for bb, xx in ref if bb == b and xx >= theta0})
        band_rows.append(dict(band=b, metric=met, n=n, state=pr["state"],
                              measured=pr["measured_winner"], cert=cert,
                              sLo=ws_lo, sHi=ws_hi, pLo=wp_lo, pHi=wp_hi,
                              sAt=row["sAt"], pAt=row["pAt"],
                              trunc_order=("ship" if ts < tp else "product"),
                              sTrunc=float(ts), pTrunc=float(tp)))

    with open(out / "model_cells.csv", "w", newline="") as h:
        w = csv.writer(h)
        w.writerow(["band", "metric", "n", "state", "measured", "cert",
                    "x_ship", "x_prod", "t_ship", "C_ship", "e_ship", "r_ship",
                    "F_ship", "t_prod", "C_prod", "e_prod", "r_prod", "F_prod",
                    "ls", "hs", "lp", "hp", "trunc_order"])
        for r in model_rows:
            s, p = r["s"], r["p"]
            w.writerow([r["band"], r["metric"], r["n"], r["state"], r["measured"],
                        r["cert"], repr(r["xs"]), repr(r["xp"]),
                        mp.nstr(s["t"], 17), repr(s["C"]), mp.nstr(s["e"], 17),
                        mp.nstr(s["r"], 17), mp.nstr(s["F"], 17),
                        mp.nstr(p["t"], 17), repr(p["C"]), mp.nstr(p["e"], 17),
                        mp.nstr(p["r"], 17), mp.nstr(p["F"], 17),
                        r["ls"], r["hs"], r["lp"], r["hp"], r["trunc_order"]])

    with open(out / "band_cert.csv", "w", newline="") as h:
        w = csv.writer(h)
        w.writerow(["band", "metric", "n", "state", "measured", "cert",
                    "sLo", "sHi", "pLo", "pHi", "sAt", "pAt",
                    "trunc_order", "sTrunc", "pTrunc"])
        for r in band_rows:
            w.writerow([r["band"], r["metric"], r["n"], r["state"], r["measured"],
                        r["cert"], r["sLo"], r["sHi"], r["pLo"], r["pHi"],
                        repr(r["sAt"]) if r["sAt"] is not None else "",
                        repr(r["pAt"]) if r["pAt"] is not None else "",
                        r["trunc_order"], r["sTrunc"], r["pTrunc"]])

    # ---- counts and obstruction statistics --------------------------------
    def score(rows):
        out = Counter()
        for r in rows:
            st = r["state"]
            c = r["cert"]
            out["all_" + (c or "overlap")] += 1
            if st in ("ship", "product"):
                out["committed"] += 1
                out["committed_" + (c or "overlap")] += 1
                if c == st:
                    out["committed_agree"] += 1
                elif c:
                    out["committed_opposite"] += 1
                if r["trunc_order"] == st:
                    out["trunc_agree"] += 1
                else:
                    out["trunc_opposite"] += 1
            else:
                out["abstain"] += 1
                out["abstain_" + (c or "overlap")] += 1
        return dict(out)

    mscore = score(model_rows)
    bscore = score(band_rows)

    # obstruction at the committed model cells.  The loser is the entry the
    # model did *not* commit (the complement of the state); the comparison is
    # in the metric's own units, so a relative row divides every quantity by
    # its own entry's denominator F.
    def scaled(d, met):
        f = 1.0 / float(d["F"]) if met == "relative" else 1.0
        return (abs(float(d["t"])) * f, float(d["C"]) * f,
                float(d["eps"]) * f, abs(float(d["r"])) * f, f)

    obstruction = []
    for r in model_rows:
        if r["state"] not in ("ship", "product"):
            continue
        s, p = r["s"], r["p"]
        winner_name = r["state"]
        loser_name = "product" if winner_name == "ship" else "ship"
        L = s if loser_name == "ship" else p
        W = p if loser_name == "ship" else s
        tl, cl, el, rl, fl = scaled(L, r["metric"])
        tw, cw, ew, rw, fw = scaled(W, r["metric"])
        # the gap a combined charge must fit under to separate the intervals
        gap = tl - tw
        need = (cl + cw + el + ew) / gap if gap > 0 else float("inf")
        obstruction.append(dict(
            band=r["band"], metric=r["metric"], n=r["n"], state=r["state"],
            loser=loser_name, winner=winner_name, mode=L["mode"],
            t_l=tl, C_l=cl, r_l=rl,
            t_w=tw, C_w=cw,
            c_seed=float(L["c_seed"]) * fl, c_rec=float(L["c_rec"]) * fl,
            c_exp=float(L["c_exp"]) * fl,
            loser_zero=bool(cl + el >= tl),
            r_gt_t=bool(rl > tl), gap=gap, need_factor=need,
            charge_over_t=(cl / tl if tl else float("inf")),
            trunc_agree=bool(gap > 0),
            overlap=float(max(r["hs"], r["hp"]) - min(r["ls"], r["lp"]))))

    def med(vals):
        vals = sorted(v for v in vals if v is not None and math.isfinite(v))
        if not vals:
            return None
        return vals[len(vals) // 2]

    # measured population band max versus the certified band interval
    pop_inside = 0
    pop_total = 0
    for b in ("A1", "A2"):
        for met in ("absolute", "relative"):
            for n in range(33):
                row = band[(b, met, n)]
                for meth in ("ship", "product"):
                    vals = []
                    for bb, x in ref:
                        if bb != b or mpf(x) < theta0:
                            continue
                        d = rec[(b, x, n, meth)]
                        if met == "absolute":
                            vals.append(abs(float(d["e"])))
                        else:
                            vals.append(abs(float(d["e"] / d["F"])))
                    m = max(vals) if vals else None
                    if m is None:
                        continue
                    lo = row["sLo"] if meth == "ship" else row["pLo"]
                    hi = row["sHi"] if meth == "ship" else row["pHi"]
                    pop_total += 1
                    if lo <= m <= hi:
                        pop_inside += 1

    summary = dict(
        cells=n_cells,
        modes=dict(modes),
        eps_base=float(EPS_BASE),
        max_abs_delivered=float(max_abs_value),
        min_abs_delivered=float(min_abs_value) if min_abs_value else None,
        max_C=float(max_C),
        validation=dict(violations=len(viol),
                        worst_ratio=float(worst_ratio),
                        worst_ratio_at=worst_ratio_at,
                        population_band_inside=[pop_inside, pop_total]),
        model_cell_counts=mscore,
        band_counts=bscore,
        committed=len(obstruction),
        committed_trunc_agree=sum(1 for o in obstruction if o["trunc_agree"]),
        committed_trunc_opposite=sum(1 for o in obstruction if not o["trunc_agree"]),
        loser_zero_count=sum(1 for o in obstruction if o["loser_zero"]),
        loser_zero_agree=sum(1 for o in obstruction
                             if o["loser_zero"] and o["trunc_agree"]),
        r_gt_t_count=sum(1 for o in obstruction if o["r_gt_t"]),
        r_gt_t_agree=sum(1 for o in obstruction
                         if o["r_gt_t"] and o["trunc_agree"]),
        need_factor_median=med([o["need_factor"] for o in obstruction
                                if math.isfinite(o["need_factor"])]),
        charge_over_t_median=med([o["charge_over_t"] for o in obstruction
                                  if o["trunc_agree"]]),
        obstruction=obstruction,
        band_rows=band_rows,
        model_rows=[{k: (v if k not in ("s", "p") else None)
                     for k, v in r.items()} for r in model_rows],
    )
    with open(out / "summary.json", "w") as h:
        json.dump(summary, h, indent=1)

    print("cells", n_cells, "modes", dict(modes))
    print("validation: violations", len(viol),
          "worst |r|/C", mp.nstr(worst_ratio, 6), "at", worst_ratio_at,
          "population band inside", pop_inside, "/", pop_total)
    print("model-cell counts", mscore)
    print("band counts", bscore)
    print("committed", len(obstruction),
          "trunc agree/opposite", summary["committed_trunc_agree"],
          "/", summary["committed_trunc_opposite"],
          "loser zero", summary["loser_zero_count"],
          "|r|>|t|", summary["r_gt_t_count"],
          "need median", summary["need_factor_median"],
          "C/|t| median", summary["charge_over_t_median"])
    if viol:
        print("VIOLATIONS:")
        for v in viol[:10]:
            print(" ", v)
    # the four smallest need factors and the largest
    finite = [o for o in obstruction if math.isfinite(o["need_factor"])]
    finite.sort(key=lambda o: o["need_factor"])
    print("smallest need factors:")
    for o in finite[:5]:
        print(" ", o["band"], o["metric"], o["n"], o["state"], o["loser"],
              "t_l=%.3e C_l=%.3e r_l=%.3e t_w=%.3e gap=%.3e need=%.1f"
              % (o["t_l"], o["C_l"], o["r_l"], o["t_w"], o["gap"],
                 o["need_factor"]))
    return 1 if viol else 0


if __name__ == "__main__":
    sys.exit(main())
