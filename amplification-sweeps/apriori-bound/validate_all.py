"""Containment and tightness of the anchored bounds at every reported cell."""

import csv
import bisect
from mpmath import mp, mpf

import bounds

mp.dps = 60

SWEEPS = {
    "amplification_form": r"C:/Users/myaml/dev/papers-a-v5/papers/boys-crossover/sweeps/amplification_form",
    "product_amplification": r"C:/Users/myaml/dev/papers-a-v5/papers/boys-crossover/sweeps/product_amplification",
}


def load_reference(root):
    ref = {}
    with open(f"{root}/data/reference.csv", newline="") as handle:
        reader = csv.reader(handle)
        next(reader)
        for row in reader:
            ref[(row[0], float(row[1]))] = [mpf(v) for v in row[2:]]
    return ref


def load_cells(root):
    cells = {}
    with open(f"{root}/results/cells.csv", newline="") as handle:
        for row in csv.DictReader(handle):
            key = (row["band"], int(row["n"]))
            cells.setdefault(key, []).append(
                (float(row["x"]), mpf(row["recursed"]), mpf(row["direct"])))
    for key in cells:
        cells[key].sort()
    return cells


def lookup(cells, band, n, x):
    table = cells[(band, n)]
    xs = [row[0] for row in table]
    i = bisect.bisect_left(xs, float(x))
    best = None
    for j in (i - 1, i, i + 1):
        if 0 <= j < len(table):
            d = abs(table[j][0] - float(x))
            if best is None or d < best[0]:
                best = (d, table[j])
    return best[1]


def per_order(root):
    rows = []
    with open(f"{root}/results/per_order.csv", newline="") as handle:
        for row in csv.DictReader(handle):
            rows.append((row["band"], int(row["n"]), row["method"],
                         float(row["at_x_abs"]), mpf(row["worst_abs"]),
                         float(row["at_x_rel"]), mpf(row["worst_rel"])))
    return rows


def main(samples=150):
    for sweep, root in SWEEPS.items():
        ref = load_reference(root)
        cells = load_cells(root)
        worst_env = (mpf(0), None)
        worst_trunc = (mpf(0), None)
        violations = []
        count = 0
        modes = {}
        for band, n, method, xa, wa, xr, wr in per_order(root):
            for metric, x, delivered in (("abs", xa, wa), ("rel", xr, wr)):
                xv, rec, direct = lookup(cells, band, n, x)
                if abs(xv - x) > 1e-15 * max(1.0, abs(x)):
                    raise SystemExit(f"no cell for {band} {n} {x!r}: nearest {xv!r}")
                value = rec if method == "ship" else direct
                exact = value - ref[(band, float(xv))][n]
                b = bounds.entry_bounds(band, n, method, xv, samples=samples,
                                        size_kind="apriori")
                modes[b["mode"]] = modes.get(b["mode"], 0) + 1
                count += 1
                if not (b["lo"] <= exact <= b["hi"]):
                    violations.append((band, n, method, metric, xv,
                                       mp.nstr(exact, 6), mp.nstr(b["lo"], 6),
                                       mp.nstr(b["hi"], 6), b["mode"]))
                if exact != 0:
                    env = max(abs(b["lo"]), abs(b["hi"])) / abs(exact)
                    tr = b["trunc"] / abs(exact)
                    if env > worst_env[0]:
                        worst_env = (env, (band, n, method, metric, xv))
                    if tr > worst_trunc[0]:
                        worst_trunc = (tr, (band, n, method, metric, xv))
        print(f"== {sweep}: {count} cell-metric rows, modes {modes}")
        print("   containment violations:", len(violations))
        for v in violations[:12]:
            print("   ", v)
        print("   worst total envelope / delivered: %.3g at %r" % (float(worst_env[0]), worst_env[1]))
        print("   worst truncation envelope / delivered: %.3g at %r" % (float(worst_trunc[0]), worst_trunc[1]))


if __name__ == "__main__":
    main()
