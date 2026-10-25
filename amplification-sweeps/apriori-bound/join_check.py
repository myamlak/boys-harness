"""Join the exact route errors with the certified intervals: containment.

For every reported cell:

1. The exact truncation error of the route (computed from the ODE identity,
   no reference) must lie inside the certified symmetric truncation envelope.
2. The delivered error (reference-measured absolute error, a magnitude) must
   lie inside the certified magnitude interval under the paper's arithmetic
   convention; a miss means the paper's convention is not a certified bound
   on the evaluation rounding at that cell.  Only the absolute-metric rows
   carry an absolute delivered error; the relative rows are used for the
   truncation check alone.
3. The same containment under a strict Clenshaw/product charge, obtained by
   widening the paper-charge interval by the difference of the two charges
   (the truncation interval itself is unchanged).
"""

import csv
import json

from mpmath import mp, mpf

import bounds

mp.dps = 60


def magnitude(lo, hi):
    return bounds.magnitude_interval(lo, hi)


def strict_charge(band, n, method, x):
    entry = bounds.entry_bounds(band, n, method, x, signed_method="coef",
                                samples=0, charge_kind="worst")
    return entry["charge"]


def main():
    with open("detail.json") as handle:
        detail = json.load(handle)
    for name, rows in detail.items():
        path = f"exact_{name}.csv"
        try:
            exact_rows = list(csv.DictReader(open(path, newline="")))
        except FileNotFoundError:
            print(f"{name}: no {path}")
            continue
        bad_trunc = 0
        bad_total = 0
        bad_strict = 0
        worst_over = (mpf(0), None)
        worst_strict_over = (mpf(0), None)
        count = 0
        for r in exact_rows:
            band, order, metric, method = r["band"], int(r["n"]), r["metric"], r["method"]
            key = f"{band}|{metric}|{order}"
            d = rows[key]
            pack = d["ship"] if method == "ship" else d["product"]
            exact = mpf(r["route_error"])
            delivered = mpf(r["delivered"])
            trunc = mpf(pack["trunc"])
            lo, hi = mpf(pack["lo"]), mpf(pack["hi"])
            count += 1
            if abs(exact) > trunc * (1 + mpf(10) ** -25):
                bad_trunc += 1
                print("   trunc miss", key, method, "exact", mp.nstr(exact, 4),
                      "trunc", mp.nstr(trunc, 4), "ratio",
                      mp.nstr(abs(exact) / trunc, 4))
            mag_lo, mag_hi = magnitude(lo, hi)
            if metric != "absolute":
                continue
            if not (mag_lo <= delivered <= mag_hi):
                bad_total += 1
                over = max(mag_lo - delivered, delivered - mag_hi)
                if over > worst_over[0]:
                    worst_over = (over, (key, method))
            charge = mpf(pack["charge"])
            sc = strict_charge(band, order, method, float(d["x_ship"] if method == "ship"
                                                         else d["x_prod"]))
            lo_s = lo + charge - sc
            hi_s = hi - charge + sc
            slo, shi = magnitude(lo_s, hi_s)
            if not (slo <= delivered <= shi):
                bad_strict += 1
                over = max(slo - delivered, delivered - shi)
                if over > worst_strict_over[0]:
                    worst_strict_over = (over, (key, method))
        print(f"{name}: {count} cells; truncation misses {bad_trunc}; "
              f"paper-charge magnitude misses {bad_total}; "
              f"strict-charge magnitude misses {bad_strict}")
        print(f"   worst paper-charge overshoot {mp.nstr(worst_over[0], 4)} "
              f"at {worst_over[1]}; worst strict-charge overshoot "
              f"{mp.nstr(worst_strict_over[0], 4)} at {worst_strict_over[1]}")


if __name__ == "__main__":
    main()
