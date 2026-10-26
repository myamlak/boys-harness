"""Compare the two sweeps: the shared arguments, the label differences, and the
per-cell structure of the same-fit region.

    python compare_sweeps.py > compare_sweeps.log 2>&1

Reads only deposited or decompressed records:

    inputs/cells.csv                                first sweep, 308616 cells
    inputs/reference.csv                            first sweep, 9352 arguments
    amplification-form/results/cells.csv            second sweep (deposited, in place)
    amplification-form/data/reference.csv
    out/first/per_order_ratio.csv                   first-sweep band-level ratios
    out/second/per_order_ratio.csv                  second-sweep band-level ratios

The first sweep's cells are decompressed from the deposit's
product-amplification/results/cells.csv.xz (308616 data rows) by
prepare_inputs.py; the second sweep is read in place from the deposit's
amplification-form/.  Nothing here writes a record.
"""

import csv
import math
import statistics
from pathlib import Path

from mpmath import mp, mpf

mp.dps = 90

HERE = Path(__file__).resolve().parent
DEPOSIT = HERE.parent
SECOND = DEPOSIT / "amplification-form"
THETA0 = 1.0855252345349333


def read_cells(path):
    """(band, x, n) -> (recursed, direct)."""
    cells = {}
    with open(path, newline="") as handle:
        for row in csv.DictReader(handle):
            cells[(row["band"], float(row["x"]), int(row["n"]))] = (
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
            reference[(row[0], float(row[1]))] = [float(v) for v in row[2:]]
    return reference


def quantile(values, q):
    values = sorted(values)
    return values[int(q * (len(values) - 1))]


def main():
    first_cells = read_cells(HERE / "inputs" / "cells.csv")
    second_cells = read_cells(SECOND / "results" / "cells.csv")
    first_ref = read_reference(HERE / "inputs" / "reference.csv")
    second_ref = read_reference(SECOND / "data" / "reference.csv")
    print("cells read: first %d, second %d" % (len(first_cells), len(second_cells)))
    print("arguments read: first %d, second %d" % (len(first_ref), len(second_ref)))

    # ---- argument overlap -------------------------------------------------
    first_args = set(first_ref)
    second_args = set(second_ref)
    shared = first_args & second_args
    print("shared arguments: %d of the second's %d (first %d, first-only %d, second-only %d)"
          % (len(shared), len(second_args), len(first_args), len(first_args - second_args),
             len(second_args - first_args)))
    print("  by band: " + "  ".join(
        "%s shared %d, first-only %d, second-only %d" % (
            band,
            sum(1 for a in shared if a[0] == band),
            sum(1 for a in first_args - second_args if a[0] == band),
            sum(1 for a in second_args - first_args if a[0] == band),
        ) for band in ("A1", "A2")))

    # ---- the reference values on the shared arguments ---------------------
    def reference_text(path):
        out = {}
        with open(path, newline="") as handle:
            reader = csv.reader(handle)
            next(reader)
            for row in reader:
                out[(row[0], float(row[1]))] = tuple(row[2:])
        return out

    first_text = reference_text(HERE / "inputs" / "reference.csv")
    second_text = reference_text(SECOND / "data" / "reference.csv")
    same_reference = sum(1 for a in shared if first_text[a] == second_text[a])
    print("shared arguments whose 33 reference values agree digit for digit: %d of %d"
          % (same_reference, len(shared)))
    worst_reference = 0.0
    worst_where = None
    for a in shared:
        for n, (ta, tb) in enumerate(zip(first_text[a], second_text[a])):
            va, vb = mpf(ta), mpf(tb)
            if va != 0:
                relative = abs(va - vb) / abs(va)
                if relative > worst_reference:
                    worst_reference = relative
                    worst_where = (a, n)
    print("worst relative reference difference on a shared argument: %s at %s order %d"
          % (mp.nstr(worst_reference, 6), worst_where[0], worst_where[1]))

    # ---- shared cells bit-identical ---------------------------------------
    identical = 0
    compared = 0
    for (band, x) in sorted(shared):
        for n in range(33):
            a = first_cells.get((band, x, n))
            b = second_cells.get((band, x, n))
            if a is None or b is None:
                continue
            compared += 1
            identical += 1 if a == b else 0
    print("shared cells compared: %d, bit-identical (both entries): %d, differing: %d"
          % (compared, identical, compared - identical))

    # ---- model population overlap -----------------------------------------
    first_pop = {(band, x) for (band, x) in first_args if x >= THETA0}
    second_pop = {(band, x) for (band, x) in second_args if x >= THETA0}
    print("population arguments (x >= theta0): first %d, second %d, shared %d,"
          " first-only %d, second-only %d"
          % (len(first_pop), len(second_pop), len(first_pop & second_pop),
             len(first_pop - second_pop), len(second_pop - first_pop)))

    # ---- measured winners: where the two sweeps disagree ------------------
    def winners(path):
        out = {}
        for row in csv.DictReader(open(path)):
            key = (row["band"], row["metric"], int(row["n"]))
            out[key] = (row["measured_winner"], float(row["rec_over_dir"]))
        return out

    w1 = winners(HERE / "out" / "first" / "per_order_ratio.csv")
    w2 = winners(HERE / "out" / "second" / "per_order_ratio.csv")
    differing = [(key, w1[key][0], w2[key][0], w1[key][1], w2[key][1])
                 for key in sorted(w1) if w1[key][0] != w2[key][0]]
    print("measured winners differing between sweeps: %d of %d comparisons"
          % (len(differing), len(w1)))
    for key, a, b, ra, rb in differing:
        print("  %s %s order %2d: first %-7s (ratio %.6f)   second %-7s (ratio %.6f)"
              % (key[0], key[1], key[2], a, ra, b, rb))

    # ---- same-fit per-cell structure, B1 orders 17..32 --------------------
    def per_cell_table(label, arguments, reference, cells, winners_map):
        print()
        print(label)
        print("order  cells  identical  rec_ok  dir_ok  median   q25      q75"
              "      frac<1   band rec/dir")
        for n in range(17, 33):
            ratios = []
            same = rec_ok = dir_ok = 0
            for (band, x) in sorted(arguments):
                if band != "A1":
                    continue
                ref = reference[(band, x)][n]
                rec, direc = cells[(band, x, n)]
                if rec == direc:
                    same += 1
                rec_ok += 1 if rec == ref else 0
                dir_ok += 1 if direc == ref else 0
                denom = abs(direc - ref)
                if denom:
                    ratios.append(abs(rec - ref) / denom)
            band_ratio = winners_map[("A1", "absolute", n)][1]
            print("%5d  %5d  %6d  %6d  %6d  %.6g  %.6g  %.6g  %.4f   %.6f"
                  % (n, len(ratios), same, rec_ok, dir_ok,
                     quantile(ratios, 0.5), quantile(ratios, 0.25), quantile(ratios, 0.75),
                     sum(1 for r in ratios if r < 1.0) / len(ratios), band_ratio))

    per_cell_table("same-fit region B1 orders 17..32, first-sweep arguments;"
                   " ratio = |recursed-ref| / |direct-ref| over cells whose direct"
                   " entry is not the correctly rounded value",
                   first_args, first_ref, first_cells, w1)
    per_cell_table("the same on the second sweep's arguments",
                   second_args, second_ref, second_cells, w2)


if __name__ == "__main__":
    main()
