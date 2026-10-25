"""A priori evaluation of the errors by the ODE identity, no reference.

Two things are measured at every reported cell:

1. The identity check: F_n computed as [stored fit value] - [its exact error
   by the anchored integral] agrees with the sweep's 60-digit reference to
   far below the delivered error.  This uses only the stored coefficients
   and exp(-x).
2. The route check: the exact truncation error of the value each entry
   delivers at its worst cell, against the delivered (reference-measured)
   error at the same cell.  The difference is the evaluation rounding.
"""

import csv
import time
from collections import Counter

from mpmath import mp, mpf

import bounds
import exact

mp.dps = 70

ROOTS = {
    "amplification_form": r"C:/Users/myaml/dev/papers-a-v5/papers/boys-crossover/sweeps/amplification_form",
    "product_amplification": r"C:/Users/myaml/dev/papers-a-v5/papers/boys-crossover/sweeps/product_amplification",
}


def reference(root):
    ref = {}
    with open(f"{root}/data/reference.csv", newline="") as handle:
        for row in csv.DictReader(handle):
            ref[(row["band"], float(row["x"]))] = [
                mpf(row[col]) for col in list(row)[2:]]
    return ref


def lookup(ref, band, x):
    try:
        return ref[(band, x)]
    except KeyError:
        best = min((abs(kx - x), kx) for b, kx in ref if b == band)
        if best[0] > 1e-14 * max(1.0, abs(x)):
            raise
        return ref[(band, best[1])]


def per_order(root):
    out = {}
    with open(f"{root}/results/per_order.csv", newline="") as handle:
        for r in csv.DictReader(handle):
            out[(r["band"], int(r["n"]), r["method"])] = dict(
                worst_abs=mpf(r["worst_abs"]), at_x_abs=float(r["at_x_abs"]),
                worst_rel=mpf(r["worst_rel"]), at_x_rel=float(r["at_x_rel"]))
    return out


def main(which):
    for name in which:
        root = ROOTS[name]
        ref = reference(root)
        order = per_order(root)
        rows = []
        id_worst = (mpf(0), None)
        id_rel = (mpf(0), None)
        t0 = time.time()
        for band in ("A1", "A2"):
            bi = 0 if band == "A1" else 1
            for n in range(33):
                for metric in ("absolute", "relative"):
                    for method in ("ship", "product"):
                        row = order[(band, n, method)]
                        x = row["at_x_abs"] if metric == "absolute" else row["at_x_rel"]
                        delivered = (row["worst_abs"] if metric == "absolute"
                                     else row["worst_rel"])
                        # identity check on the band's own per-order fit
                        e_fit = exact.fit_exact_error(n, bi, x)
                        p = bounds.fit_poly(n, bi).value(x)
                        f_id = p - e_fit
                        fref = lookup(ref, band, x)[n]
                        d = abs(f_id - fref)
                        if d > id_worst[0]:
                            id_worst = (d, (band, n, metric, method, x))
                        rel = d / fref if fref != 0 else mp.inf
                        if rel > id_rel[0]:
                            id_rel = (rel, (band, n, metric, method, x))
                        # route truncation error and its ratio to delivered
                        e_route = exact.entry_truncation_exact(band, n, method, x)
                        rows.append(dict(
                            band=band, n=n, metric=metric, method=method, x=repr(x),
                            f_ref=mp.nstr(fref, 12),
                            identity_delta=mp.nstr(d, 6),
                            fit_error=mp.nstr(e_fit, 8),
                            route_error=mp.nstr(e_route, 8),
                            delivered=str(delivered),
                            ratio=(mp.nstr(abs(delivered / e_route), 6)
                                   if e_route != 0 else ""),
                            rounding=mp.nstr(delivered - e_route, 6)))
                print(f"   {band} n={n:2d} at {time.time() - t0:6.1f}s", flush=True)
        print(f"== {name}: {len(rows)} cells in {time.time() - t0:.1f}s")
        print(f"   identity: max |F_identity - F_ref| = {mp.nstr(id_worst[0], 6)} "
              f"at {id_worst[1]}")
        print(f"   identity: max relative = {mp.nstr(id_rel[0], 6)} at {id_rel[1]}")
        ratios = [float(r["ratio"]) for r in rows if r["ratio"]]
        ratios.sort()
        if ratios:
            print("   |delivered| / |route truncation|: min %.3g q10 %.3g "
                  "median %.3g q90 %.3g max %.3g"
                  % (ratios[0], ratios[len(ratios) // 10],
                     ratios[len(ratios) // 2], ratios[9 * len(ratios) // 10],
                     ratios[-1]))
        dominant = Counter()
        for r in rows:
            if not r["ratio"]:
                continue
            q = float(r["ratio"])
            if q > 10:
                dominant["rounding (>10x)"] += 1
            elif q < 0.1:
                dominant["truncation (>10x)"] += 1
            else:
                dominant["mixed"] += 1
        print("   dominance:", dict(dominant))
        with open(f"C:/Users/myaml/AppData/Local/Temp/apriori-bound/exact_{name}.csv",
                  "w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
        print(f"   wrote exact_{name}.csv")


if __name__ == "__main__":
    import sys
    which = tuple(sys.argv[1:]) or ("amplification_form",)
    main(which)
