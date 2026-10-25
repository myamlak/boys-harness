"""Why the loser's certified interval always contains zero.

For each of the paper's 67 committed comparisons, split the loser's interval
into its truncation part (signed, no charge) and the arithmetic charge:

* L_trunc > 0 : the truncation interval alone excludes zero; the charge is
  what drags the certified interval across it, and the ratio charge/L_trunc
  says how much the charge would have to shrink.
* L_trunc = 0 : the truncation interval already straddles zero, so no
  arithmetic tightening can ever certify the ordering.
"""

import csv
import json

from mpmath import mp, mpf

import bounds

mp.dps = 60

ROOT = (r"C:/Users/myaml/dev/papers-a-v5/papers/boys-crossover/sweeps/"
        "amplification_form")


def trunc_interval(band, n, method, x):
    """Signed [lo, hi] on the truncation error, no arithmetic charge."""
    x = mpf(x)
    bi = 0 if band == "A1" else 1
    if method == "product":
        return bounds.fit_error(n, bi, x, "sampled", 400)
    if x < bounds.THETA0:
        transfer = bounds.amplification(n, x) / bounds.amplification(bounds.NMAX, x)
        lo, hi = bounds.fit_error(bounds.NMAX, 0, x, "sampled", 400)
        return transfer * lo, transfer * hi
    if n == 0 or x >= mpf(bounds.TABLE["thresholds"][n]):
        lo, hi, _ = bounds.seed_error(x, "sampled", 400)
        A = bounds.amplification(n, x)
        return A * lo, A * hi
    return bounds.fit_error(n, bi, x, "sampled", 400)


def main():
    pred = {}
    with open(f"{ROOT}/results/prediction.csv", newline="") as handle:
        for r in csv.DictReader(handle):
            pred[(r["band"], r["metric"], int(r["order"]))] = r
    order = {}
    with open(f"{ROOT}/results/per_order.csv", newline="") as handle:
        for r in csv.DictReader(handle):
            order[(r["band"], int(r["n"]), r["method"])] = r
    committed = sorted(k for k, r in pred.items()
                       if r["state"] in ("ship", "product"))
    out = dict(straddle=0, charge_drags=0, rows=[])
    ratios = []
    for band, metric, n in committed:
        xa = "at_x_abs" if metric == "absolute" else "at_x_rel"
        da = "worst_abs" if metric == "absolute" else "worst_rel"
        xs = float(order[(band, n, "ship")][xa])
        xp = float(order[(band, n, "product")][xa])
        ds = mpf(order[(band, n, "ship")][da])
        dp = mpf(order[(band, n, "product")][da])
        loser = "ship" if ds > dp else "product"
        x = xs if loser == "ship" else xp
        lo, hi = trunc_interval(band, n, loser, x)
        L = bounds.magnitude_interval(lo, hi)[0]
        charge = bounds.entry_bounds(band, n, loser, x, signed_method="coef",
                                     samples=0)["charge"]
        row = dict(band=band, metric=metric, n=n, loser=loser, x=x,
                   trunc_lo=mp.nstr(lo, 6), trunc_hi=mp.nstr(hi, 6),
                   L_trunc=mp.nstr(L, 6), charge=mp.nstr(charge, 6))
        if L == 0:
            out["straddle"] += 1
            row["case"] = "straddle"
        else:
            out["charge_drags"] += 1
            ratio = float(charge / L)
            ratios.append(ratio)
            row["case"] = "charge"
            row["charge_over_L"] = ratio
        out["rows"].append(row)
    ratios.sort()
    if ratios:
        out["ratio_min"] = ratios[0]
        out["ratio_median"] = ratios[len(ratios) // 2]
        out["ratio_max"] = ratios[-1]
    print("committed", len(committed), "straddle", out["straddle"],
          "charge-drags", out["charge_drags"])
    if ratios:
        print("charge / L_trunc: min %.3g median %.3g max %.3g"
              % (ratios[0], ratios[len(ratios) // 2], ratios[-1]))
    with open("sign_diag.json", "w") as handle:
        json.dump(out, handle, indent=1)


if __name__ == "__main__":
    main()
