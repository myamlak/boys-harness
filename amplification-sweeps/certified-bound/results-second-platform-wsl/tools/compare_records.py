"""Compare second-platform record CSVs with the deposited repaired record.

Usage: compare_records.py DEPOSIT_DIR RUN_DIR
Compares cells_cert.csv, band_cert.csv, model_cells.csv and harness_cells.csv
by their natural keys.  For each file prints row counts, key equality, the
number of differing rows and, per column, the number of differing values and
the largest relative difference of the numeric ones.
"""
import csv
import math
import sys
from pathlib import Path

KEYS = {
    "cells_cert.csv": ("band", "x", "n", "method"),
    "band_cert.csv": ("band", "metric", "n"),
    "model_cells.csv": ("band", "metric", "n"),
    "harness_cells.csv": ("band", "metric", "n"),
}


def load(path):
    rows = {}
    order = []
    with open(path, newline="") as h:
        reader = csv.DictReader(h)
        fields = reader.fieldnames
        for r in reader:
            key = tuple(r[k] for k in KEYS[Path(path).name])
            rows[key] = r
            order.append(key)
    return rows, order, fields


def num(s):
    try:
        return float(s)
    except (TypeError, ValueError):
        return None


def main():
    dep, run = Path(sys.argv[1]), Path(sys.argv[2])
    for name, keys in KEYS.items():
        a, oa, fa = load(dep / name)
        b, ob, fb = load(run / name)
        print("==", name)
        print("  rows", len(a), len(b), "keys and order equal", oa == ob,
              "fields equal", fa == fb)
        diff_rows = 0
        col_diff = {f: 0 for f in fa}
        col_maxrel = {f: 0.0 for f in fa}
        col_examples = {f: [] for f in fa}
        for key in oa:
            ra, rb = a[key], b.get(key)
            if rb is None:
                continue
            differed = False
            for f in fa:
                va, vb = ra[f], rb[f]
                if va == vb:
                    continue
                differed = True
                col_diff[f] += 1
                na, nb = num(va), num(vb)
                if na is not None and nb is not None:
                    denom = max(abs(na), abs(nb))
                    if denom:
                        rel = abs(na - nb) / denom
                        col_maxrel[f] = max(col_maxrel[f], rel)
                if len(col_examples[f]) < 3:
                    col_examples[f].append((key, va, vb))
            if differed:
                diff_rows += 1
        print("  differing rows", diff_rows)
        for f in fa:
            if col_diff[f]:
                ex = col_examples[f][0]
                print("   %-12s diff=%d maxrel=%.3e first %s: %s -> %s"
                      % (f, col_diff[f], col_maxrel[f], ex[0], ex[1], ex[2]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
