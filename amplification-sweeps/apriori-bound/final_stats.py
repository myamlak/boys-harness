"""The remaining numbers for APRIORI-BOUND.md, in one reproducible log.

1. Strict-charge absolute envelope / delivered at the 132 absolute entries of
   each sweep (the paper-charge version is in join2_*.log).
2. At the 67 committed comparisons: the winner's certified upper over its
   delivered error, paper and strict, and the cell of the paper-charge
   containment minimum.
3. The exact-truncation oracle: how much the combined arithmetic charge
   would have to shrink for the truncation interval to order the delivered
   winner, and how many comparisons no reduction can fix.
"""

import csv
from collections import Counter

from mpmath import mp, mpf, inf

import bounds

mp.dps = 60

SWEEPS = {
    "amplification_form": r"C:/Users/myaml/dev/papers-a-v5/papers/boys-crossover/sweeps/amplification_form",
    "product_amplification": r"C:/Users/myaml/dev/papers-a-v5/papers/boys-crossover/sweeps/product_amplification",
}


def stats(values):
    vals = sorted(float(v) for v in values if v != inf)
    n = len(vals)
    if n == 0:
        return dict(n=0)

    def q(p):
        return vals[min(n - 1, max(0, int(round(p * (n - 1)))))]
    return dict(n=n, min=vals[0], q10=q(0.1), median=q(0.5), q90=q(0.9),
                max=vals[-1])


def read_order(root):
    out = {}
    with open(f"{root}/results/per_order.csv", newline="") as handle:
        for r in csv.DictReader(handle):
            out[(r["band"], int(r["n"]), r["method"])] = dict(
                worst_abs=mpf(r["worst_abs"]), at_x_abs=float(r["at_x_abs"]),
                worst_rel=mpf(r["worst_rel"]), at_x_rel=float(r["at_x_rel"]))
    return out


def read_ref(root):
    out = {}
    with open(f"{root}/data/reference.csv", newline="") as handle:
        for r in csv.DictReader(handle):
            out[(r["band"], float(r["x"]))] = [mpf(r[c]) for c in list(r)[2:]]
    return out


def strict_ratio(name, root):
    order = read_order(root)
    rows = []
    for band in ("A1", "A2"):
        for n in range(33):
            for meth in ("ship", "product"):
                o = order[(band, n, meth)]
                x, delivered = o["at_x_abs"], o["worst_abs"]
                if delivered == 0:
                    continue
                paper = bounds.entry_bounds(band, n, meth, x,
                                            signed_method="sampled", samples=400,
                                            charge_kind="paper")
                strict = bounds.entry_bounds(band, n, meth, x,
                                             signed_method="sampled", samples=400,
                                             charge_kind="worst")
                rows.append(dict(
                    band=band, n=n, meth=meth, x=x, delivered=delivered,
                    paper=max(abs(paper["lo"]), abs(paper["hi"])),
                    strict=max(abs(strict["lo"]), abs(strict["hi"]))))
    paper = [r["paper"] / r["delivered"] for r in rows]
    strict = [r["strict"] / r["delivered"] for r in rows]
    pm = sum(1 for r in rows if r["delivered"] > r["paper"])
    sm = sum(1 for r in rows if r["delivered"] > r["strict"])
    print(name, "absolute entries", len(rows),
          "paper ratio", stats(paper), "paper misses", pm,
          "strict ratio", stats(strict), "strict misses", sm)


def committed(name, root):
    order = read_order(root)
    ref = read_ref(root)
    pred = {}
    with open(f"{root}/results/prediction.csv", newline="") as handle:
        for r in csv.DictReader(handle):
            pred[(r["band"], r["metric"], int(r["order"]))] = r
    exact = {}
    with open(f"exact_{name}.csv", newline="") as handle:
        for r in csv.DictReader(handle):
            exact[(r["band"], int(r["n"]), r["metric"], r["method"])] = dict(
                t=mpf(r["route_error"]), x=float(r["x"]))
    committed = sorted(k for k, r in pred.items()
                       if r["state"] in ("ship", "product"))
    winner_ratio = {"paper": [], "worst": []}
    min_cell = None
    reduce_ratio = {"paper": [], "worst": []}
    unfixable = 0
    one_signed_trunc = 0
    for band, metric, n in committed:
        xa = "at_x_abs" if metric == "absolute" else "at_x_rel"
        da = "worst_abs" if metric == "absolute" else "worst_rel"
        xs = order[(band, n, "ship")][xa]
        xp = order[(band, n, "product")][xa]
        ds = order[(band, n, "ship")][da]
        dp = order[(band, n, "product")][da]
        loser = "ship" if ds > dp else "product"
        winner = "product" if loser == "ship" else "ship"
        xw = xs if winner == "ship" else xp
        dw = ds if winner == "ship" else dp
        for kind in ("paper", "worst"):
            b = bounds.entry_bounds(band, n, winner, xw,
                                    signed_method="sampled", samples=400,
                                    charge_kind=kind)
            u = max(abs(b["lo"]), abs(b["hi"]))
            winner_ratio[kind].append(u / dw)
            if kind == "paper" and (min_cell is None
                                    or u / dw < min_cell[0]):
                min_cell = (float(u / dw), band, metric, n, winner, xw,
                            float(u), float(dw))
        xl = xs if loser == "ship" else xp
        tl = exact[(band, n, metric, loser)]
        tw = exact[(band, n, metric, winner)]
        t_l, t_w = tl["t"], tw["t"]
        cl = bounds.entry_bounds(band, n, loser, xl,
                                 signed_method="sampled", samples=400)
        cw = bounds.entry_bounds(band, n, winner, xw,
                                 signed_method="sampled", samples=400)
        if metric == "relative":
            f_l = ref[(band, tl["x"])][n]
            f_w = ref[(band, tw["x"])][n]
            t_l, t_w = t_l / f_l, t_w / f_w
            cl = dict(cl, charge=cl["charge"] / f_l)
            cw = dict(cw, charge=cw["charge"] / f_w)
        if bounds.magnitude_interval(t_l - cl["charge"],
                                     t_l + cl["charge"])[0] > 0:
            one_signed_trunc += 1
        if abs(t_l) <= abs(t_w):
            unfixable += 1
            continue
        for kind, key in (("paper", "paper"), ("worst", "worst")):
            slack = abs(t_l) - abs(t_w)
            if slack <= 0:
                continue
            c = cl["charge"] + cw["charge"] if kind == "paper" else None
            if kind == "worst":
                clw = bounds.entry_bounds(band, n, loser, xl,
                                          signed_method="sampled", samples=400,
                                          charge_kind="worst")
                cww = bounds.entry_bounds(band, n, winner, xw,
                                          signed_method="sampled", samples=400,
                                          charge_kind="worst")
                clc, cwc = clw["charge"], cww["charge"]
                if metric == "relative":
                    clc, cwc = clc / f_l, cwc / f_w
                c = clc + cwc
            reduce_ratio[key].append(float(c / slack) if slack else inf)
    print(name, "committed", len(committed),
          "winner upper/delivered paper", stats(winner_ratio["paper"]),
          "worst", stats(winner_ratio["worst"]))
    print("   paper-charge winner minimum at", min_cell)
    print("   loser truncation interval one-signed (exact, paper charge):",
          one_signed_trunc, "of", len(committed))
    print("   truncation picks the other entry (no charge reduction helps):",
          unfixable, "of", len(committed))
    print("   combined charge reduction needed, paper:",
          stats(reduce_ratio["paper"]))
    print("   combined charge reduction needed, worst:",
          stats(reduce_ratio["worst"]))


def main():
    for name, root in SWEEPS.items():
        strict_ratio(name, root)
    committed("amplification_form", SWEEPS["amplification_form"])


if __name__ == "__main__":
    main()
