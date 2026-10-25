"""Relative-metric certification with the sampled signed interval.

run_all.py certifies the relative comparisons with the coarse coefficient
interval (signed_method='coef') for both the entry error and the F_n divisor,
so its variant table never exercises sampling on the relative metric.  This
script repeats those 66 comparisons per sweep per variant with the sampled
method so the sampling question is answered on both metrics.
"""

import csv
from collections import Counter

from mpmath import mp, mpf

import bounds

mp.dps = 60

SWEEPS = {
    "amplification_form": r"C:/Users/myaml/dev/papers-a-v5/papers/boys-crossover/sweeps/amplification_form",
    "product_amplification": r"C:/Users/myaml/dev/papers-a-v5/papers/boys-crossover/sweeps/product_amplification",
}

VARIANTS = {
    "sampled100_paper": dict(signed_method="sampled", samples=100,
                             size_kind="apriori", charge_kind="paper"),
    "sampled400_paper": dict(signed_method="sampled", samples=400,
                             size_kind="apriori", charge_kind="paper"),
    "sampled1600_paper": dict(signed_method="sampled", samples=1600,
                              size_kind="apriori", charge_kind="paper"),
    "sampled400_paper_strictsize": dict(signed_method="sampled", samples=400,
                                        size_kind="bound", charge_kind="paper"),
    "sampled400_worstcharge": dict(signed_method="sampled", samples=400,
                                   size_kind="apriori", charge_kind="worst"),
}


def per_order(root):
    out = {}
    with open(f"{root}/results/per_order.csv", newline="") as handle:
        for r in csv.DictReader(handle):
            out[(r["band"], int(r["n"]), r["method"])] = r
    return out


def main():
    for name, root in SWEEPS.items():
        order = per_order(root)
        for vname, kw in VARIANTS.items():
            counts = Counter()
            for band in ("A1", "A2"):
                for n in range(33):
                    xs = float(order[(band, n, "ship")]["at_x_rel"])
                    xp = float(order[(band, n, "product")]["at_x_rel"])
                    bs = bounds.entry_bounds(band, n, "ship", xs, **kw)
                    bp = bounds.entry_bounds(band, n, "product", xp, **kw)
                    ls, hs = bounds.relative_interval(
                        band, n, xs, bs, **{k: v for k, v in kw.items()
                                            if k in ("samples", "signed_method")})
                    lp, hp = bounds.relative_interval(
                        band, n, xp, bp, **{k: v for k, v in kw.items()
                                            if k in ("samples", "signed_method")})
                    if hs < lp:
                        counts["ship"] += 1
                    elif hp < ls:
                        counts["product"] += 1
                    else:
                        counts["abstain"] += 1
            print(name, vname, dict(counts), flush=True)


if __name__ == "__main__":
    main()
