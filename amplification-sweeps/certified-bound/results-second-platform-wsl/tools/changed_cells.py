"""Margin at the cells whose certified record differs from the deposit.

Usage: changed_cells.py DEPOSIT_RUN RUN
Reports how many rows are touched by any value change, the largest |r|/C
among them, the largest increase in |r|/C, and whether every cell is inside.
"""
import csv
import sys
from pathlib import Path


def load(path):
    return {(r["band"], r["x"], r["n"], r["method"]): r
            for r in csv.DictReader(open(path))}


def main():
    a = load(Path(sys.argv[1]) / "cells_cert.csv")
    b = load(Path(sys.argv[2]) / "cells_cert.csv")
    touched = 0
    worst = (0.0, None)
    worst_delta = (0.0, None)
    for k, rb in b.items():
        ra = a[k]
        if (rb["r"], rb["C"]) == (ra["r"], ra["C"]):
            continue
        touched += 1
        c = float(rb["C"])
        ratio = abs(float(rb["r"])) / c if c else 0.0
        ca = float(ra["C"])
        ratio_a = abs(float(ra["r"])) / ca if ca else 0.0
        if ratio > worst[0]:
            worst = (ratio, k)
        if ratio - ratio_a > worst_delta[0]:
            worst_delta = (ratio - ratio_a, k)
    print("rows touched by any r or C change", touched)
    print("max |r|/C among them %.17g at %s" % worst)
    print("max increase in |r|/C %.17g at %s" % worst_delta)
    print("all inside", all(r["inside"] == "1" for r in b.values()))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
