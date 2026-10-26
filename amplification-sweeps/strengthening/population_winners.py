"""The first sweep scored against measured winners restricted to the model's population.

    python population_winners.py > population_winners.log 2>&1

The model's population is x >= theta0.  Its measured winner, however, is the
entry with the smaller band-level worst over the whole band, and a band worst
may sit below theta0, where the kernel does not run the upward recurrence.
This script recomputes both band maxima twice, once over the whole band and
once over x >= theta0, and rescoring the first-sweep commitments against the
restricted winners shows how much of the first-sweep error is a band-max
argument the model does not cover (main.tex, the paragraph after the model
definition) rather than a wrong decision on the population.

Reads inputs/cells.csv and inputs/reference.csv (decompressed from the
deposit's product-amplification streams) and the first-sweep prediction files
written by score_sweeps.py.  Cell values enter as doubles, the reference as
60-digit mpf, as in the deposited analysis.
"""

import csv
from pathlib import Path

from mpmath import mp, mpf

mp.dps = 70

HERE = Path(__file__).resolve().parent
THETA0 = 1.0855252345349333


def read_reference(path):
    reference = {}
    with open(path, newline="") as handle:
        reader = csv.reader(handle)
        next(reader)
        for row in reader:
            reference[(row[0], float(row[1]))] = [mpf(v) for v in row[2:]]
    return reference


def winner(rec, dir_):
    if rec == dir_:
        return "tie"
    return "ship" if rec < dir_ else "product"


def main():
    reference = read_reference(HERE / "inputs" / "reference.csv")
    best = {}
    cells = 0
    with open(HERE / "inputs" / "cells.csv", newline="") as handle:
        for row in csv.DictReader(handle):
            cells += 1
            band = row["band"]
            x = float(row["x"])
            n = int(row["n"])
            expected = reference[(band, x)][n]
            rec = abs(mpf(float(row["recursed"])) - expected)
            dir_ = abs(mpf(float(row["direct"])) - expected)
            rec_rel = rec / expected if expected != 0 else mpf("inf")
            dir_rel = dir_ / expected if expected != 0 else mpf("inf")
            scopes = (("whole", True), ("pop", x >= THETA0))
            for metric, r_value, d_value in (("absolute", rec, dir_),
                                             ("relative", rec_rel, dir_rel)):
                for scope, keep in scopes:
                    if not keep:
                        continue
                    for entry, value in (("ship", r_value), ("product", d_value)):
                        key = (band, metric, n, scope, entry)
                        if key not in best or value > best[key][0]:
                            best[key] = (value, x)
    print("cells read:", cells, "reference arguments:", len(reference))
    print()

    def winner_of(band, metric, n, scope):
        return winner(best[(band, metric, n, scope, "ship")][0],
                      best[(band, metric, n, scope, "product")][0])

    changes = []
    for band in ("A1", "A2"):
        for metric in ("absolute", "relative"):
            for n in range(33):
                whole = winner_of(band, metric, n, "whole")
                pop = winner_of(band, metric, n, "pop")
                if whole != pop:
                    changes.append((band, metric, n, whole, pop))
    print("comparisons whose measured winner changes when the band maxima are"
          " restricted to x >= theta0: %d of 132" % len(changes))
    for band, metric, n, whole, pop in changes:
        print("  %s %-8s order %2d: whole %-7s -> population %-7s" % (band, metric, n, whole, pop))
        for scope in ("whole", "pop"):
            ship_value, ship_x = best[(band, metric, n, scope, "ship")]
            prod_value, prod_x = best[(band, metric, n, scope, "product")]
            print("      %-5s ship %.17g at x=%.17g   product %.17g at x=%.17g"
                  % (scope, ship_value, ship_x, prod_value, prod_x))
    print()

    for label, path in (("margin 1", HERE / "out" / "first" / "prediction.csv"),
                        ("margin 0", HERE / "out" / "first_m0" / "prediction.csv")):
        rows = list(csv.DictReader(open(path)))
        committed = [r for r in rows if r["state"] != "abstain"]
        whole_right = sum(1 for r in committed if r["right_where_committed"] == "1")
        population_right = 0
        population_wrong = []
        for r in committed:
            band, metric, n = r["band"], r["metric"], int(r["order"])
            choice = r["choice"]
            if winner_of(band, metric, n, "pop") == "tie":
                continue
            if winner_of(band, metric, n, "pop") == choice:
                population_right += 1
            else:
                population_wrong.append((band, metric, n, choice,
                                         winner_of(band, metric, n, "pop"),
                                         r["measured_winner"]))
        print("%s: commits %d, right against the whole-band winners %d, wrong %d"
              % (label, len(committed), whole_right, len(committed) - whole_right))
        print("  against the population-restricted winners: right %d, wrong %d,"
              " abstain-at-tie %d"
              % (population_right, len(population_wrong),
                 len(committed) - population_right - len(population_wrong)))
        for row in population_wrong:
            print("    wrong: %s %s order %d chose %s, population winner %s"
                  " (whole-band winner %s)" % row)
        for r in committed:
            band, metric, n = r["band"], r["metric"], int(r["order"])
            if winner_of(band, metric, n, "pop") == "tie":
                print("    tie: %s %s order %d, choice %s, whole-band winner %s"
                      % (band, metric, n, r["choice"],
                         winner_of(band, metric, n, "whole")))


if __name__ == "__main__":
    main()
