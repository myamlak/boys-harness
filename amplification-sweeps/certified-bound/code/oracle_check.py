"""Cross-check the fast exact-real route evaluator against the deposit's
independent identity oracle on a deterministic 60-cell sample.

The oracle is exact.entry_truncation_exact (the anchored ODE identity with
quadrature, mp.dps=70); the fast evaluator is exactroute.py.  Agreement is
evidence that the truncation t at every cell is the same quantity the
deposit certifies; the record's EPS_BASE=1e-45 is far above the observed
difference.

Paths: --sweep selects the sweep data, --deposit the in-repo a priori-bound
deposit, --out the run directory (record_paths.py).
"""
import argparse
import csv
import os
import random
import sys

from mpmath import mp, mpf

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import exactroute as er
import model
import record_paths

mp.dps = 70


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    record_paths.add_kernel(parser)
    record_paths.add_inputs(parser, deposit_arg=True)
    record_paths.add_out(parser)
    args = parser.parse_args(argv)
    model.configure(kernel=args.kernel, exp_bits=args.exp_bits)
    root = record_paths.sweep(args.sweep)
    dep = record_paths.deposit(args.deposit)
    out = record_paths.out_dir(args.out)

    with record_paths.deposit_on_path(dep):
        import exact
    mp.dps = 70

    ref = er.load_reference(root)
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
    with open(out / "oracle_check.csv", "w", newline="") as h:
        w = csv.writer(h)
        w.writerow(["band", "x", "n", "method", "route_truncation",
                    "oracle_truncation", "diff"])
        w.writerows(rows)
    print("cells", len(rows), "max |diff|", mp.nstr(worst, 8))
    return 0


if __name__ == "__main__":
    sys.exit(main())
