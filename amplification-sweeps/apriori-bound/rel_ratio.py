"""Relative-metric envelope/delivered ratios at the per-order worst cells.

The absolute-metric ratios are in run_all.py's table.  This computes the
matching relative ratios with the certificate route actually tested in
run_all.py (coefficient intervals for both the entry error and F_n), so the
ratio is between two relative-error magnitudes.
"""

import csv

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

    def q(p):
        return vals[min(n - 1, max(0, int(round(p * (n - 1)))))]
    return dict(n=n, min=vals[0], q10=q(0.1), median=q(0.5),
                q90=q(0.9), max=vals[-1])


def main():
    for name, root in SWEEPS.items():
        order = {}
        with open(f"{root}/results/per_order.csv", newline="") as handle:
            for r in csv.DictReader(handle):
                order[(r["band"], int(r["n"]), r["method"])] = r
        ratios = []
        contains = 0
        misses = 0
        for band in ("A1", "A2"):
            for n in range(33):
                for meth in ("ship", "product"):
                    r = order[(band, n, meth)]
                    x = float(r["at_x_rel"])
                    delivered = mpf(r["worst_rel"])
                    eb = bounds.entry_bounds(band, n, meth, x,
                                             signed_method="coef", samples=0)
                    lo, hi = bounds.relative_interval(band, n, x, eb,
                                                      signed_method="coef")
                    if delivered != 0:
                        ratios.append(hi / delivered)
                    if lo <= delivered <= hi:
                        contains += 1
                    else:
                        misses += 1
        print(name, "relative envelope/delivered:", stats(ratios),
              "containment", contains, "misses", misses)


if __name__ == "__main__":
    main()
