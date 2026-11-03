"""Rigorous-truncation tier and the seven a priori variants, with the
certified per-call rounding radius replacing the heuristic charge.

For every one of the 132 band-metric comparisons, at the sweep's own
per-order cells, this script obtains the deposit's certified truncation
interval for each entry under each of the seven variants (sampled signed
residual at 100/400/1600 points, coefficient and crude envelopes, a priori
or rigorous step sizes, paper or worst charge), then forms

    total interval = [trunc_lo - C, trunc_hi + C]

with C the running-error radius of model.py at that cell, and scores the
pair.  It also scores the deposit's own charge for the same variant, so the
reproduction of its 0-certified count can be checked beside ours.

Outputs results/tiers_cells.csv, one row per (comparison, variant).
"""
import csv
import itertools
import json
import sys

from mpmath import mp, mpf

sys.path.insert(0, r"C:/Users/myaml/dev/_a-cert-run/code")
sys.path.insert(0, r"C:/Users/myaml/dev/_harness-v1023/amplification-sweeps/apriori-bound")
import model
import bounds

mp.dps = 60
ROOT = r"C:/Users/myaml/dev/papers-a-v5/papers/boys-crossover/sweeps/amplification_form"
OUT = r"C:/Users/myaml/dev/_a-cert-run/results"
F_REL = mpf("1e-55")

VARIANTS = {
    "primary_sampled400_paper": dict(signed_method="sampled", samples=400,
                                     size_kind="apriori", charge_kind="paper"),
    "sampled100_paper": dict(signed_method="sampled", samples=100,
                             size_kind="apriori", charge_kind="paper"),
    "sampled1600_paper": dict(signed_method="sampled", samples=1600,
                              size_kind="apriori", charge_kind="paper"),
    "sampled400_paper_strictsize": dict(signed_method="sampled", samples=400,
                                        size_kind="bound", charge_kind="paper"),
    "sampled400_worstcharge": dict(signed_method="sampled", samples=400,
                                   size_kind="apriori", charge_kind="worst"),
    "coef_paper": dict(signed_method="coef", samples=0,
                       size_kind="apriori", charge_kind="paper"),
    "crude_paper": dict(signed_method="crude", samples=0,
                        size_kind="apriori", charge_kind="paper"),
}


def load():
    pred = {}
    with open(ROOT + "/results/prediction.csv", newline="") as h:
        for r in csv.DictReader(h):
            pred[(r["band"], r["metric"], int(r["order"]))] = r
    order = {}
    with open(ROOT + "/results/per_order.csv", newline="") as h:
        for r in csv.DictReader(h):
            order[(r["band"], int(r["n"]), r["method"])] = r
    ref = {}
    with open(ROOT + "/data/reference.csv", newline="") as h:
        for r in csv.DictReader(h):
            ref[(r["band"], float(r["x"]))] = [mpf(r[c]) for c in list(r)[2:]]
    exp = model.load_exp()
    return pred, order, ref, exp


def magnitude(lo, hi):
    if lo <= 0 <= hi:
        return mpf(0), max(abs(lo), abs(hi))
    return min(abs(lo), abs(hi)), max(abs(lo), abs(hi))


def main():
    pred, order, ref, exp = load()
    cache = {}
    c_cache = {}

    def eb(band, n, method, x, variant):
        key = (band, n, method, x, variant)
        v = cache.get(key)
        if v is None:
            v = bounds.entry_bounds(band, n, method, x, **VARIANTS[variant])
            cache[key] = v
        return v

    def radius(band, n, method, x):
        key = (band, n, method, x)
        v = c_cache.get(key)
        if v is None:
            bi = 0 if band == "A1" else 1
            if method == "ship":
                v = model.ship_batch(x, exp[x][1])[n]["radius"]
            else:
                v = model.product_batch(x, bi)[n]["radius"]
            c_cache[key] = v
        return mpf(v)

    rows = []
    for (band, metric, n), pr in sorted(pred.items()):
        key = "at_x_abs" if metric == "absolute" else "at_x_rel"
        xs = float(order[(band, n, "ship")][key])
        xp = float(order[(band, n, "product")][key])
        fs = ref[(band, xs)][n]
        fp = ref[(band, xp)][n]
        cs = radius(band, n, "ship", xs)
        cp = radius(band, n, "product", xp)
        for variant in VARIANTS:
            es = eb(band, n, "ship", xs, variant)
            ep = eb(band, n, "product", xp, variant)
            # the deposit's own verdict with its charge
            if metric == "absolute":
                dep = bounds.certify(es["lo"], es["hi"],
                                     ep["lo"], ep["hi"])
            else:
                # run_all.py's relative verdict uses the coef F_n interval
                dep = bounds.certify_relative(band, n, xs, xp, es, ep,
                                              signed_method="coef")
            # its truncation-only interval, plus our certified radius
            sl = es["lo"] + es["arith"]
            sh = es["hi"] - es["arith"]
            pl = ep["lo"] + ep["arith"]
            ph = ep["hi"] - ep["arith"]
            sl, sh, pl, ph = sl - cs, sh + cs, pl - cp, ph + cp
            ms = magnitude(sl, sh)
            mp_ = magnitude(pl, ph)
            if metric == "relative":
                ms = (ms[0] / (fs * (1 + F_REL)), ms[1] / (fs * (1 - F_REL)))
                mp_ = (mp_[0] / (fp * (1 + F_REL)), mp_[1] / (fp * (1 - F_REL)))
            if ms[1] < mp_[0]:
                cert = "ship"
            elif mp_[1] < ms[0]:
                cert = "product"
            else:
                cert = ""
            rows.append(dict(band=band, metric=metric, n=n,
                             state=pr["state"], variant=variant,
                             dep_cert=dep, cert=cert,
                             sl=ms[0], sh=ms[1], pl=mp_[0], ph=mp_[1],
                             cs=float(cs), cp=float(cp),
                             dch_s=float(es["arith"]), dch_p=float(ep["arith"])))

    with open(OUT + "/tiers_cells.csv", "w", newline="") as h:
        w = csv.writer(h)
        w.writerow(["band", "metric", "n", "state", "variant", "dep_cert",
                    "cert", "ship_lo", "ship_hi", "prod_lo", "prod_hi",
                    "C_ship", "C_prod", "dep_charge_ship", "dep_charge_prod"])
        for r in rows:
            w.writerow([r["band"], r["metric"], r["n"], r["state"],
                        r["variant"], r["dep_cert"], r["cert"],
                        mp.nstr(r["sl"], 8), mp.nstr(r["sh"], 8),
                        mp.nstr(r["pl"], 8), mp.nstr(r["ph"], 8),
                        repr(r["cs"]), repr(r["cp"]),
                        mp.nstr(r["dch_s"], 8), mp.nstr(r["dch_p"], 8)])

    summary = {}
    for variant in VARIANTS:
        rs = [r for r in rows if r["variant"] == variant]
        committed = [r for r in rs if r["state"] in ("ship", "product")]
        def counts(sel, field):
            c = {"certified": sum(1 for r in sel if r[field]),
                 "overlap": sum(1 for r in sel if not r[field])}
            agree = sum(1 for r in sel if r[field] == r["state"])
            opposite = sum(1 for r in sel if r[field] and r[field] != r["state"])
            c.update(agree=agree, opposite=opposite)
            return c
        summary[variant] = dict(
            all_deposit=counts(rs, "dep_cert"),
            all_certified_charge=counts(rs, "cert"),
            committed_deposit=counts(committed, "dep_cert"),
            committed_certified_charge=counts(committed, "cert"))
        print(variant,
              "all: deposit", summary[variant]["all_deposit"],
              "our C", summary[variant]["all_certified_charge"])
        print("   committed: deposit",
              summary[variant]["committed_deposit"],
              "our C", summary[variant]["committed_certified_charge"])
    with open(OUT + "/tiers_summary.json", "w") as h:
        json.dump(summary, h, indent=1)


if __name__ == "__main__":
    main()
