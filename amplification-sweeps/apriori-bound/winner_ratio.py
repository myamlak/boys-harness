"""The winner's certified upper over its delivered error, metric-aware.

Absolute commitments compare the absolute entry bound with worst_abs;
relative commitments compare the relative interval (entry bound divided by a
certified interval on F_n) with worst_rel.  The earlier mixed-scale line in
detail.log is superseded by this.
"""

import csv

from mpmath import mp, mpf, inf

import bounds

mp.dps = 60

ROOT = (r"C:/Users/myaml/dev/papers-a-v5/papers/boys-crossover/sweeps/"
        "amplification_form")


def stats(values):
    vals = sorted(float(v) for v in values if v != inf)
    n = len(vals)

    def q(p):
        return vals[min(n - 1, max(0, int(round(p * (n - 1)))))]
    return dict(n=n, min=vals[0], q10=q(0.1), median=q(0.5), q90=q(0.9),
                max=vals[-1])


def main():
    pred = {}
    with open(f"{ROOT}/results/prediction.csv", newline="") as handle:
        for r in csv.DictReader(handle):
            pred[(r["band"], r["metric"], int(r["order"]))] = r
    order = {}
    with open(f"{ROOT}/results/per_order.csv", newline="") as handle:
        for r in csv.DictReader(handle):
            order[(r["band"], int(r["n"]), r["method"])] = dict(
                worst_abs=mpf(r["worst_abs"]), at_x_abs=float(r["at_x_abs"]),
                worst_rel=mpf(r["worst_rel"]), at_x_rel=float(r["at_x_rel"]))
    committed = sorted(k for k, r in pred.items()
                       if r["state"] in ("ship", "product"))
    out = {}
    for kind in ("paper", "worst"):
        for metric in ("absolute", "relative", "all"):
            ratios = []
            misses = 0
            n_cells = 0
            for band, met, n in committed:
                if metric != "all" and met != metric:
                    continue
                xa = "at_x_abs" if met == "absolute" else "at_x_rel"
                da = "worst_abs" if met == "absolute" else "worst_rel"
                xs = order[(band, n, "ship")][xa]
                xp = order[(band, n, "product")][xa]
                ds = order[(band, n, "ship")][da]
                dp = order[(band, n, "product")][da]
                winner = "ship" if ds < dp else "product"
                xw = xs if winner == "ship" else xp
                dw = ds if winner == "ship" else dp
                kw = dict(signed_method="coef", samples=0) if met == "relative" \
                    else dict(signed_method="sampled", samples=400)
                b = bounds.entry_bounds(band, n, winner, xw, charge_kind=kind,
                                        **kw)
                if met == "relative":
                    lo, hi = bounds.relative_interval(band, n, xw, b,
                                                      signed_method="coef")
                else:
                    lo, hi = bounds.magnitude_interval(b["lo"], b["hi"])
                upper = hi if met == "relative" else max(abs(b["lo"]),
                                                         abs(b["hi"]))
                if dw > upper:
                    misses += 1
                ratios.append(upper / dw)
                n_cells += 1
            out[(kind, metric)] = (stats(ratios), misses, n_cells)
            print(kind, metric, stats(ratios), "misses", misses, "of", n_cells)


if __name__ == "__main__":
    main()
