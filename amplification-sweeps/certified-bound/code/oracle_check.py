"""Cross-check the fast exact-real route evaluator against the deposit's
independent identity oracle on a deterministic 60-cell sample.

The oracle is exact.entry_truncation_exact (the anchored ODE identity with
quadrature, mp.dps=70); the fast evaluator is exactroute.py.  Agreement is
evidence that the truncation t at every cell is the same quantity the
deposit certifies; the record's EPS_BASE=1e-45 is far above the observed
difference.
"""
import csv
import random
import sys

from mpmath import mp, mpf

sys.path.insert(0, r"C:/Users/myaml/dev/_a-cert-run/code")
sys.path.insert(0, r"C:/Users/myaml/dev/_harness-v1023/amplification-sweeps/apriori-bound")
import exactroute as er
import exact

mp.dps = 70
ROOT = r"C:/Users/myaml/dev/papers-a-v5/papers/boys-crossover/sweeps/amplification_form"
OUT = r"C:/Users/myaml/dev/_a-cert-run/results/oracle_check.csv"


def main():
    ref = er.load_reference(ROOT)
    pool = []
    for (band, x), vals in ref.items():
        for n in (0, 1, 2, 5, 8, 13, 17, 24, 31, 32):
            if n < len(vals):
                pool.append((band, x, n, "ship"))
                pool.append((band, x, n, "product"))
    rng = random.Random(20261103)
    cells = rng.sample(pool, 60)
    rows = []
    worst = mpf(0)
    for band, x, n, method in cells:
        if method == "ship":
            route = er.ship_values(x, mp.exp(-mpf(x)) / 2)[n]
        else:
            route = er.product_values(x)[n]
        ours = route - ref[(band, x)][n]
        oracle = exact.entry_truncation_exact(band, n, method, x)
        diff = ours - oracle
        worst = max(worst, abs(diff))
        rows.append((band, x, n, method, mp.nstr(ours, 17),
                     mp.nstr(oracle, 17), mp.nstr(diff, 6)))
        print(band, repr(x), n, method, "diff", mp.nstr(diff, 6))
    with open(OUT, "w", newline="") as h:
        w = csv.writer(h)
        w.writerow(["band", "x", "n", "method", "route_truncation",
                    "oracle_truncation", "diff"])
        w.writerows(rows)
    print("cells", len(rows), "max |diff|", mp.nstr(worst, 8))


if __name__ == "__main__":
    main()
