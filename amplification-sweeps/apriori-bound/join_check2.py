"""Containment of the delivered absolute errors, recomputed from the tables.

For each sweep and each absolute-metric row of exact_<sweep>.csv:

* the certified truncation envelope must contain the exact route truncation
  error computed by the identity (no reference);
* the delivered absolute error (per_order.csv, full precision) must lie in
  the certified magnitude interval under the paper's charge and under the
  strict Clenshaw/product charge.

Also prints, per entry, the certified magnitude upper bound over the
delivered error (the conservativeness of the envelope at the reported cell).
"""

import csv
import sys

from mpmath import mp, mpf

import bounds

mp.dps = 60

ROOTS = {
    "amplification_form": r"C:/Users/myaml/dev/papers-a-v5/papers/boys-crossover/sweeps/amplification_form",
    "product_amplification": r"C:/Users/myaml/dev/papers-a-v5/papers/boys-crossover/sweeps/product_amplification",
}


def magnitude(lo, hi):
    return bounds.magnitude_interval(lo, hi)


def per_order(root):
    out = {}
    with open(f"{root}/results/per_order.csv", newline="") as handle:
        for r in csv.DictReader(handle):
            out[(r["band"], int(r["n"]), r["method"])] = (
                mpf(r["worst_abs"]), float(r["at_x_abs"]))
    return out


def main(which):
    for name in which:
        root = ROOTS[name]
        order = per_order(root)
        rows = []
        with open(f"exact_{name}.csv", newline="") as handle:
            for r in csv.DictReader(handle):
                if r["metric"] == "absolute":
                    rows.append(r)
        paper_miss = strict_miss = trunc_miss = 0
        worst_paper = (mpf(0), None)
        worst_strict = (mpf(0), None)
        ratio = []
        for r in rows:
            band, n, method = r["band"], int(r["n"]), r["method"]
            x = float(r["x"])
            delivered, _ = order[(band, n, method)]
            exact = mpf(r["route_error"])
            pack = bounds.entry_bounds(band, n, method, x,
                                       signed_method="sampled", samples=400,
                                       size_kind="apriori", charge_kind="paper")
            trunc = pack["trunc"]
            if abs(exact) > trunc * (1 + mpf(10) ** -20):
                trunc_miss += 1
                print("   truncation miss", band, n, method, "exact",
                      mp.nstr(exact, 6), "trunc", mp.nstr(trunc, 6))
            lo, hi = pack["lo"], pack["hi"]
            lo_m, hi_m = magnitude(lo, hi)
            if not (lo_m <= delivered <= hi_m):
                paper_miss += 1
                over = max(lo_m - delivered, delivered - hi_m)
                if over > worst_paper[0]:
                    worst_paper = (over, (band, n, method, x))
            strict = bounds.entry_bounds(band, n, method, x,
                                         signed_method="coef", samples=0,
                                         charge_kind="worst")["charge"]
            lo_s = lo + pack["charge"] - strict
            hi_s = hi - pack["charge"] + strict
            lo_s, hi_s = magnitude(lo_s, hi_s)
            if not (lo_s <= delivered <= hi_s):
                strict_miss += 1
                over = max(lo_s - delivered, delivered - hi_s)
                if over > worst_strict[0]:
                    worst_strict = (over, (band, n, method, x))
            if delivered != 0:
                ratio.append(float(hi_m / delivered))
        ratio.sort()
        print(f"{name}: {len(rows)} absolute cells; truncation misses "
              f"{trunc_miss}; paper-charge containment misses {paper_miss}; "
              f"strict-charge containment misses {strict_miss}")
        print(f"   worst paper overshoot {mp.nstr(worst_paper[0], 4)} at "
              f"{worst_paper[1]}")
        print(f"   worst strict overshoot {mp.nstr(worst_strict[0], 4)} at "
              f"{worst_strict[1]}")
        if ratio:
            print(f"   certified upper / delivered: min {ratio[0]:.3g} "
                  f"q10 {ratio[len(ratio) // 10]:.3g} "
                  f"median {ratio[len(ratio) // 2]:.3g} "
                  f"q90 {ratio[9 * len(ratio) // 10]:.3g} "
                  f"max {ratio[-1]:.3g}")


if __name__ == "__main__":
    main(tuple(sys.argv[1:]) or ("amplification_form",))
