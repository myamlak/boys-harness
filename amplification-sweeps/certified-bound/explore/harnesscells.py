"""The 65 commitments at the sweep's per-order worst cells.

Reproduces the paper's own cell choice (per_order at_x) and reports the
exact-truncation ordering, the certified charge's verdict, and the paper
charge's needed reduction, for comparison with the paper's 16 / 49 / 93
and the deposit's 0/67.

The loser is the entry the model did *not* commit.  The paper's charge is
its own convention: 3u per upward step with the reference F_l carried by
A_n/A_l, plus eta_E = 2^-52 times the exponential sensitivity; the direct
entry is charged 3u times the coefficient l1 mass.
"""
import csv
import math
import sys

from mpmath import mp, mpf

sys.path.insert(0, r"C:/Users/myaml/dev/_a-cert-run/code")
sys.path.insert(0, r"C:/Users/myaml/dev/_harness-v1023/amplification-sweeps/apriori-bound")
import exactroute as er
import model
import bounds

mp.dps = 70
ROOT = r"C:/Users/myaml/dev/papers-a-v5/papers/boys-crossover/sweeps/amplification_form"
EXACT_CSV = r"C:/Users/myaml/dev/_harness-v1023/amplification-sweeps/apriori-bound/exact_amplification_form.csv"
EPS = mpf("1e-45")
U = mpf(2) ** -53


def load():
    ref = er.load_reference(ROOT)
    cells = {}
    with open(ROOT + "/results/cells.csv", newline="") as h:
        for r in csv.DictReader(h):
            cells[(r["band"], float(r["x"]), int(r["n"]))] = (
                float(r["recursed"]), float(r["direct"]))
    pred = {}
    with open(ROOT + "/results/prediction.csv", newline="") as h:
        for r in csv.DictReader(h):
            pred[(r["band"], r["metric"], int(r["order"]))] = r
    order = {}
    with open(ROOT + "/results/per_order.csv", newline="") as h:
        for r in csv.DictReader(h):
            order[(r["band"], int(r["n"]), r["method"])] = r
    exact = {}
    try:
        with open(EXACT_CSV, newline="") as h:
            for r in csv.DictReader(h):
                exact[(r["band"], int(r["n"]), r["metric"], r["method"])] = (
                    mpf(r["route_error"]), float(r["x"]))
    except FileNotFoundError:
        pass
    exp = model.load_exp()
    return ref, cells, pred, order, exact, exp


def paper_charge(band, n, method, x, ref):
    """The paper's accumulated-arithmetic charge, absolute."""
    bi = 0 if band == "A1" else 1
    if method == "product":
        return 3 * U * bounds.fit_poly(n, bi).coef_l1
    xm = mpf(x)
    values = ref[(band, x)]
    total = mpf(0)
    if xm < bounds.THETA0:
        for l in range(n, bounds.NMAX + 1):
            total += 3 * U * values[l] * bounds.amplification_ratio(n, l, xm)
        return total
    for l in range(1, n + 1):
        total += 3 * U * values[l] * bounds.amplification_ratio(n, l, xm)
    if n > 0:
        carried = mpf(0)
        for l in range(n):
            carried += bounds.amplification_ratio(n, l, xm) / (2 * l + 1)
        total += mpf(2) ** -52 * mp.exp(-xm) * carried
    return total


def median(vals):
    vals = sorted(v for v in vals if v is not None and math.isfinite(v))
    if not vals:
        return None
    return vals[len(vals) // 2]


def main():
    ref, cells, pred, order, exact, exp = load()
    rows = []
    for (band, metric, n), pr in sorted(pred.items()):
        if pr["state"] not in ("ship", "product"):
            continue
        key = "at_x_abs" if metric == "absolute" else "at_x_rel"
        xs = float(order[(band, n, "ship")][key])
        xp = float(order[(band, n, "product")][key])
        bi = 0 if band == "A1" else 1
        sv = er.ship_values(xs, mp.exp(-mpf(xs)) / 2)
        pv = er.product_values(xp)
        fs = ref[(band, xs)][n]
        fp = ref[(band, xp)][n]
        ts = sv[n] - fs
        tp = pv[n] - fp
        sc = model.ship_batch(xs, exp[xs][1])[n]
        pc = model.product_batch(xp, bi)[n]
        cs = mpf(sc["radius"])
        cp = mpf(pc["radius"])
        ps = paper_charge(band, n, "ship", xs, ref)
        pp = paper_charge(band, n, "product", xp, ref)
        # validation of the delivered arithmetic at these cells
        # delivered minus the exact-real route value: this *is* the rounding
        es = mpf(cells[(band, xs, n)][0]) - sv[n]
        ep = mpf(cells[(band, xp, n)][1]) - pv[n]
        rs, rp = abs(es), abs(ep)
        es_signed, ep_signed = es, ep
        # exact-CSV cross-check
        xcheck = None
        if (band, n, metric, "ship") in exact:
            t0, x0 = exact[(band, n, metric, "ship")]
            if abs(x0 - xs) < 1e-300:
                xcheck = abs(t0 - ts)
        # metric units
        if metric == "relative":
            ts, tp, cs, cp, ps, pp = (ts / fs, tp / fp, cs / fs, cp / fp,
                                      ps / fs, pp / fp)
            rs, rp = rs / fs, rp / fp
            es_signed, ep_signed = es_signed / fs, ep_signed / fp
        # certified intervals
        ms = bounds.magnitude_interval(ts - EPS - cs, ts + EPS + cs)
        mpp = bounds.magnitude_interval(tp - EPS - cp, tp + EPS + cp)
        cert = ""
        if ms[1] < mpp[0]:
            cert = "ship"
        elif mpp[1] < ms[0]:
            cert = "product"
        # the confirmed ordering is the model's state; the loser is the other
        state = pr["state"]
        trunc_order = "ship" if abs(ts) < abs(tp) else "product"
        loser = "product" if state == "ship" else "ship"
        tl, cl, pl = ((abs(ts), cs, ps) if loser == "ship"
                      else (abs(tp), cp, pp))
        tl_signed, rl_signed = ((ts, es_signed) if loser == "ship"
                                else (tp, ep_signed))
        tw, cw, pw = ((abs(ts), cs, ps) if loser != "ship"
                      else (abs(tp), cp, pp))
        gap = tl - tw
        need_cert = (cl + cw) / gap if gap > 0 else float("inf")
        need_paper = (pl + pw) / gap if gap > 0 else float("inf")
        rows.append(dict(band=band, metric=metric, n=n, state=state, cert=cert,
                         trunc_order=trunc_order, xs=xs, xp=xp,
                         ts=ts, tp=tp, cs=cs, cp=cp, ps=ps, pp=pp,
                         loser=loser, gap=gap,
                         loser_ratio_paper=float(pl / tl) if tl else float("inf"),
                         loser_ratio_cert=float(cl / tl) if tl else float("inf"),
                         loser_t_signed=tl_signed,
                         loser_r_signed=rl_signed,
                         need_paper=need_paper, need_cert=need_cert,
                         r_s=rs, r_p=rp, xcheck=xcheck))

    from collections import Counter
    cnt = Counter()
    for r in rows:
        cnt["committed"] += 1
        cnt["cert_" + (r["cert"] or "overlap")] += 1
        if r["cert"] == r["state"]:
            cnt["cert_agree"] += 1
        elif r["cert"]:
            cnt["cert_opposite"] += 1
        if r["trunc_order"] == r["state"]:
            cnt["trunc_agree"] += 1
        else:
            cnt["trunc_opposite"] += 1
    print("counts", dict(cnt))
    agree = [r for r in rows if r["trunc_order"] == r["state"]]
    opp = [r for r in rows if r["trunc_order"] != r["state"]]
    print("trunc agree", len(agree), "opposite", len(opp))
    print("paper charge/loser |t|: median",
          median([r["loser_ratio_paper"] for r in rows]),
          " cert charge/loser |t|: median",
          median([r["loser_ratio_cert"] for r in rows]))
    print("over the", len(agree), "agreeing cells:")
    print("  paper combined charge / gap: median",
          median([r["need_paper"] for r in agree]),
          "  certified combined charge / gap: median",
          median([r["need_cert"] for r in agree]))
    def ratio(r, blind):
        t = r["loser_t_signed"]
        if t == 0:
            return float("inf")
        return float((abs(t) + abs(r["loser_r_signed"])) / abs(t) if blind
                     else abs(t + r["loser_r_signed"]) / abs(t))
    for label, sel in (("all", rows), ("agreeing", agree)):
        print("loser delivered/truncation over", label, "cells:",
              "signed median", median([ratio(r, False) for r in sel]),
              "sign-blind median", median([ratio(r, True) for r in sel]))
    print("charge violations: cert",
          sum(1 for r in rows if r["r_s"] > r["cs"] or r["r_p"] > r["cp"]),
          "paper",
          sum(1 for r in rows if r["r_s"] > r["ps"] or r["r_p"] > r["pp"]))
    xchecks = [float(r["xcheck"]) for r in rows if r["xcheck"] is not None]
    print("exact CSV cross-check:", len(xchecks),
          "max", max(xchecks) if xchecks else None)

    with open(r"C:/Users/myaml/dev/_a-cert-run/results/harness_cells.csv", "w",
              newline="") as h:
        w = csv.writer(h)
        w.writerow(["band", "metric", "n", "state", "loser", "cert",
                    "trunc_order", "x_ship", "x_prod", "t_ship", "t_prod",
                    "C_ship", "C_prod", "paper_ship", "paper_prod",
                    "gap", "ratio_paper", "ratio_cert", "need_paper",
                    "need_cert", "r_ship", "r_prod"])
        for r in rows:
            w.writerow([r["band"], r["metric"], r["n"], r["state"], r["loser"],
                        r["cert"], r["trunc_order"], repr(r["xs"]),
                        repr(r["xp"]), mp.nstr(r["ts"], 17),
                        mp.nstr(r["tp"], 17), mp.nstr(r["cs"], 17),
                        mp.nstr(r["cp"], 17), mp.nstr(r["ps"], 17),
                        mp.nstr(r["pp"], 17), mp.nstr(r["gap"], 8),
                        r["loser_ratio_paper"], r["loser_ratio_cert"],
                        r["need_paper"], r["need_cert"],
                        mp.nstr(r["r_s"], 17), mp.nstr(r["r_p"], 17)])


if __name__ == "__main__":
    main()
