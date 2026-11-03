"""Per-cell certified intervals at the 132 comparisons' model cells.

Uses the exact-real route value (truncation) with the record's tier-T
precision bound, plus the model's running-error radius.  Reports the
certified direction, the loser's interval, and the structural split.
"""
import csv
import math
import sys

sys.path.insert(0, r"C:/Users/myaml/dev/_a-cert-run/code")
import exactroute as er
import model

from mpmath import mp, mpf

mp.dps = 70
ROOT = r"C:/Users/myaml/dev/papers-a-v5/papers/boys-crossover/sweeps/amplification_form"
EPS_BASE = mpf("1e-45")


def load():
    ref = er.load_reference(ROOT)
    cells = {}
    with open(ROOT + "/results/cells.csv", newline="") as h:
        for r in csv.DictReader(h):
            cells[(r["band"], float(r["x"]), int(r["n"]))] = (
                float(r["recursed"]), float(r["direct"]))
    pred = {}
    with open(ROOT + "/results/prediction_record.csv", newline="") as h:
        for r in csv.DictReader(h):
            pred[(r["band"], r["metric"], int(r["order"]))] = r
    order = {}
    with open(ROOT + "/results/per_order.csv", newline="") as h:
        for r in csv.DictReader(h):
            order[(r["band"], int(r["n"]), r["method"])] = r
    exp = model.load_exp()
    return ref, cells, pred, order, exp


def cell_interval(band, n, method, x, ref, cells, exp, cache):
    key = (band, n, method, x)
    if key in cache:
        return cache[key]
    bi = 0 if band == "A1" else 1
    if method == "ship":
        vals = er.ship_values(x, mp.exp(-mpf(x)) / 2)
        C = model.ship_batch(x, exp[x][1])[n]
    else:
        vals = er.product_values(x)
        C = model.product_batch(x, bi)[n]
    F = ref[(band, x)][n]
    t = vals[n] - F
    tf = float(t)
    eps = EPS_BASE + mpf(2) ** -51 * abs(tf)
    es = mpf(cells[(band, x, n)][0]) - F if method == "ship" else \
         mpf(cells[(band, x, n)][1]) - F
    route = vals[n]
    delivered = mpf(cells[(band, x, n)][0] if method == "ship"
                   else cells[(band, x, n)][1])
    r = delivered - route
    out = dict(t=t, tf=tf, eps=eps, C=mpf(C["radius"]), F=F, e=es, r=r,
               mode=C["mode"], c_seed=mpf(C.get("c_seed", 0.0)),
               c_rec=mpf(C.get("c_rec", 0.0)), c_exp=mpf(C.get("c_exp", 0.0)))
    out["lo"] = t - eps - out["C"]
    out["hi"] = t + eps + out["C"]
    cache[key] = out
    return out


def mag(lo, hi):
    if lo <= 0 <= hi:
        return mpf(0), max(abs(lo), abs(hi))
    return min(abs(lo), abs(hi)), max(abs(lo), abs(hi))


def rel_interval(entry, F):
    flo = F * (1 - mpf("1e-55"))
    fhi = F * (1 + mpf("1e-55"))
    vals = [entry["lo"] / flo, entry["lo"] / fhi,
            entry["hi"] / flo, entry["hi"] / fhi]
    return min(vals), max(vals)


def compare(s, p, metric):
    if metric == "absolute":
        ls, hs = mag(s["lo"], s["hi"])
        lp, hp = mag(p["lo"], p["hi"])
    else:
        ls, hs = rel_interval(s, s["F"])
        lp, hp = rel_interval(p, p["F"])
    if hs < lp:
        return "ship", (ls, hs), (lp, hp)
    if hp < ls:
        return "product", (ls, hs), (lp, hp)
    return "", (ls, hs), (lp, hp)


def main():
    ref, cells, pred, order, exp = load()
    cache = {}
    rows = []
    for (band, metric, n), r in sorted(pred.items()):
        xs = float(r["predicted_ship_at_x"])
        xp = float(r["predicted_product_at_x"])
        s = cell_interval(band, n, "ship", xs, ref, cells, exp, cache)
        p = cell_interval(band, n, "product", xp, ref, cells, exp, cache)
        cert, si, pi = compare(s, p, metric)
        state = r["state"]
        loser = "ship" if state == "ship" else ("product" if state == "product" else "")
        rows.append(dict(band=band, metric=metric, n=n, state=state,
                         xs=xs, xp=xp, cert=cert, s=s, p=p, si=si, pi=pi))

    # certified counts
    from collections import Counter
    cnt = Counter()
    committed = [r for r in rows if r["state"] in ("ship", "product")]
    for r in rows:
        cnt["all_" + (r["cert"] or "overlap")] += 1
        if r["state"] in ("ship", "product"):
            cnt["committed_" + (r["cert"] or "overlap")] += 1
            if r["cert"] == r["state"]:
                cnt["committed_agree"] += 1
            elif r["cert"]:
                cnt["committed_opposite"] += 1
        else:
            cnt["abstain_" + (r["cert"] or "overlap")] += 1
    print("counts", dict(cnt))

    # structural split at the committed pairs
    split = Counter()
    details = []
    for r in committed:
        s, p = r["s"], r["p"]
        ls = "ship" if r["state"] == "ship" else "product"
        L = s if ls == "ship" else p
        tl = abs(L["t"])
        Cl = L["C"]
        rl = abs(L["r"])
        split["loser_zero" if Cl + L["eps"] >= tl else "loser_onesigned"] += 1
        split["r_gt_t" if rl >= tl else "r_le_t"] += 1
        details.append((r["band"], r["metric"], r["n"], r["state"], ls,
                        float(tl), float(Cl), float(rl),
                        float(L["c_seed"]), float(L["c_rec"]), float(L["c_exp"])))
    print("split", dict(split))
    details.sort(key=lambda d: d[7] / d[5] if d[5] else 1e300)
    print("loser t, C, |r|, c_seed, c_rec, c_exp (smallest |r|/|t| first)")
    for d in details[:8]:
        print(" ", d[:4], "loser", d[4], "t=%.3e C=%.3e r=%.3e seed=%.3e rec=%.3e exp=%.3e"
              % (d[5], d[6], d[7], d[8], d[9], d[10]))
    print("  ...")
    for d in details[-4:]:
        print(" ", d[:4], "loser", d[4], "t=%.3e C=%.3e r=%.3e"
              % (d[5], d[6], d[7]))

    # validation of the charge at these cells
    viol = [r for r in rows for e in (r["s"], r["p"]) if abs(e["r"]) > e["C"]]
    print("charge violations at model cells:", len(viol))
    worst = max((abs(e["r"]) / e["C"] for r in rows for e in (r["s"], r["p"])),
                default=mpf(0))
    print("worst |r|/C:", mp.nstr(worst, 6))

    with open(r"C:/Users/myaml/dev/_a-cert-run/results/model_cells.csv", "w",
              newline="") as h:
        w = csv.writer(h)
        w.writerow(["band", "metric", "n", "state", "cert", "x_ship", "x_prod",
                    "t_ship", "C_ship", "e_ship", "r_ship", "F_ship",
                    "t_prod", "C_prod", "e_prod", "r_prod", "F_prod"])
        for r in rows:
            s, p = r["s"], r["p"]
            w.writerow([r["band"], r["metric"], r["n"], r["state"], r["cert"],
                        repr(r["xs"]), repr(r["xp"]),
                        mp.nstr(s["t"], 17), mp.nstr(s["C"], 17),
                        mp.nstr(s["e"], 17), mp.nstr(s["r"], 17),
                        mp.nstr(s["F"], 17),
                        mp.nstr(p["t"], 17), mp.nstr(p["C"], 17),
                        mp.nstr(p["e"], 17), mp.nstr(p["r"], 17),
                        mp.nstr(p["F"], 17)])


if __name__ == "__main__":
    main()
