"""The B1 orders 17-32 region: whole-band maxima against the model population.

    python same_fit_region.py > same_fit_region.log 2>&1

The referee's claim is a whole-band statement: at B1 orders 17-32 the recursed
entry has the smaller band worst in absolute error at 15 consecutive orders by
about a factor two, and the larger relative error at all 16 orders.  The
model's population is x >= theta0; below theta0 the recursed entry is a
downward recursion seeded at order 32 and the model decides nothing.  This
script prints, per order and per metric, the band maxima over the whole band
and over the population alone, so the claim can be read against the region the
model actually covers.  It also prints the per-cell distribution inside the
population, because a band maximum is a tail, not a typical cell.

Inputs are the deposited records (the first sweep's inputs/cells.csv and
inputs/reference.csv are decompressed from the deposit's xz files by
prepare_inputs.py; the second sweep is read in place).
"""

import csv
import statistics
from pathlib import Path

from mpmath import mp, mpf

mp.dps = 70

HERE = Path(__file__).resolve().parent
DEPOSIT = HERE.parent
SECOND = DEPOSIT / "amplification-form"
THETA0 = 1.0855252345349333
BAND = "A1"
ORDERS = range(17, 33)


def read_cells(path):
    cells = {}
    with open(path, newline="") as handle:
        for row in csv.DictReader(handle):
            if row["band"] != BAND:
                continue
            cells[(float(row["x"]), int(row["n"]))] = (
                float(row["recursed"]),
                float(row["direct"]),
            )
    return cells


def read_reference(path):
    reference = {}
    with open(path, newline="") as handle:
        reader = csv.reader(handle)
        next(reader)
        for row in reader:
            if row[0] == BAND:
                reference[float(row[1])] = [mpf(v) for v in row[2:]]
    return reference


def band_maxima(cells, reference, population_only):
    """(metric, order, entry) -> (worst, x_text)."""
    maxima = {}
    for (x, n), (rec, direct) in cells.items():
        if n not in ORDERS:
            continue
        if population_only and x < THETA0:
            continue
        expected = reference[x][n]
        rec_error = abs(mpf(rec) - expected)
        dir_error = abs(mpf(direct) - expected)
        for metric, rec_value, dir_value in (
            ("absolute", rec_error, dir_error),
            ("relative", rec_error / expected, dir_error / expected),
        ):
            for entry, value in (("rec", rec_value), ("dir", dir_value)):
                key = (metric, n, entry)
                if key not in maxima or value > maxima[key][0]:
                    maxima[key] = (value, x)
    return maxima


def per_cell(cells, reference):
    stats = []
    for (x, n), (rec, direct) in cells.items():
        if n not in ORDERS or x < THETA0:
            continue
        expected = float(reference[x][n])
        rec_ok = 1 if rec == expected else 0
        dir_ok = 1 if direct == expected else 0
        denom = abs(direct - expected)
        ratio = abs(rec - expected) / denom if denom else None
        stats.append((n, rec == direct, rec_ok, dir_ok, ratio))
    return stats


def run(label, cells_path, reference_path):
    cells = read_cells(cells_path)
    reference = read_reference(reference_path)
    print("=" * 78)
    print(label, "sweep, band", BAND, "orders 17..32")
    print("cells:", len(cells), "arguments:", len(reference))
    full = band_maxima(cells, reference, population_only=False)
    pop = band_maxima(cells, reference, population_only=True)
    print()
    print("order  metric    full-band rec/dir  ratio   population rec/dir  ratio")
    for n in ORDERS:
        for metric in ("absolute", "relative"):
            fa = full[(metric, n, "rec")][0] / full[(metric, n, "dir")][0]
            pa = pop[(metric, n, "rec")][0] / pop[(metric, n, "dir")][0]
            print("%5d  %-8s  %.9g / %.9g  %.6f  %.9g / %.9g  %.6f"
                  % (n, metric,
                     full[(metric, n, "rec")][0], full[(metric, n, "dir")][0], fa,
                     pop[(metric, n, "rec")][0], pop[(metric, n, "dir")][0], pa))
    print()
    print("population worsts at (order, metric, entry): x")
    for n in ORDERS:
        for metric in ("absolute", "relative"):
            for entry in ("rec", "dir"):
                print("  %2d %-8s %s at x = %s"
                      % (n, metric, entry, mp.nstr(mpf(pop[(metric, n, entry)][1]), 8)))

    stats = per_cell(cells, reference)
    ratios = sorted(ratio for n, same, rec_ok, dir_ok, ratio in stats
                    if ratio is not None)

    def q(p):
        return ratios[int(p * (len(ratios) - 1))]

    print()
    print("per-cell distribution inside the population, orders 17..32:")
    print("  cells %d, bit-identical %d, rec correctly rounded %d, dir correctly rounded %d"
          % (len(stats), sum(1 for s in stats if s[1]),
             sum(s[2] for s in stats), sum(s[3] for s in stats)))
    print("  ratio quantiles over %d cells with a non-zero direct error:"
          " q10 %.6g q25 %.6g median %.6g q75 %.6g q90 %.6g"
          % (len(ratios), q(0.10), q(0.25), q(0.50), q(0.75), q(0.90)))
    print("  fraction of cells with ratio < 1: %.4f"
          % (sum(1 for r in ratios if r < 1.0) / len(ratios)))


def main():
    run("first", HERE / "inputs" / "cells.csv", HERE / "inputs" / "reference.csv")
    run("second", SECOND / "results" / "cells.csv", SECOND / "data" / "reference.csv")


if __name__ == "__main__":
    main()
