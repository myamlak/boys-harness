"""Recompute the paper's certified-bound statistics from a run directory.

Usage: paper_stats.py RUN_DIR
Reads harness_cells.csv, cells_cert.csv, band_cert.csv and model_cells.csv
and prints the counts, the midpoint medians, the required-factor range and
the validation totals.
"""
import csv
import math
import sys
from pathlib import Path


def med(vals):
    vals = sorted(v for v in vals if v is not None and math.isfinite(v))
    n = len(vals)
    if n == 0:
        return None
    return vals[n // 2] if n % 2 else (vals[n // 2 - 1] + vals[n // 2]) / 2.0


def main():
    run = Path(sys.argv[1])
    rows = list(csv.DictReader(open(run / "harness_cells.csv")))
    print("harness rows", len(rows))
    print("states", {s: sum(1 for r in rows if r["state"] == s)
                     for s in ("ship", "product")})
    agree = [r for r in rows if r["trunc_order"] == r["state"]]
    opp = [r for r in rows if r["trunc_order"] != r["state"]]
    print("truncation agree", len(agree), "opposite", len(opp))
    ship = [float(r["ratio_cert"]) for r in rows if r["loser"] == "ship"]
    prod = [float(r["ratio_cert"]) for r in rows if r["loser"] == "product"]
    print("loser recursed (ship) n=%d median C/|t| %.17g" % (len(ship), med(ship)))
    print("loser direct (product) n=%d median C/|t| %.17g" % (len(prod), med(prod)))
    need = [float(r["need_cert"]) for r in agree]
    print("over the %d agreeing: required factor median %.17g min %.17g max %.17g"
          % (len(need), med(need), min(need), max(need)))
    cert_verdicts = {v: sum(1 for r in rows if r["cert"] == v)
                     for v in ("", "ship", "product")}
    print("harness cert verdicts", cert_verdicts)

    cells = list(csv.DictReader(open(run / "cells_cert.csv")))
    inside = [r for r in cells if r["inside"] == "1"]
    worst = (0.0, None)
    for r in cells:
        c = float(r["C"])
        if c > 0:
            ratio = abs(float(r["r"])) / c
            if ratio > worst[0]:
                worst = (ratio, (r["band"], r["x"], r["n"], r["method"]))
    print("cells_cert rows", len(cells), "inside", len(inside),
          "violations", len(cells) - len(inside))
    print("max |r|/C from CSV %.17g at %s" % worst)

    band = list(csv.DictReader(open(run / "band_cert.csv")))
    print("band rows", len(band),
          "cert verdicts", {v: sum(1 for r in band if r["cert"] == v)
                            for v in ("", "ship", "product")})
    committed = [r for r in band if r["state"] in ("ship", "product")]
    print("band committed", len(committed),
          "trunc agree", sum(1 for r in committed if r["trunc_order"] == r["state"]),
          "opposite", sum(1 for r in committed if r["trunc_order"] != r["state"]))
    model = list(csv.DictReader(open(run / "model_cells.csv")))
    print("model rows", len(model),
          "cert verdicts", {v: sum(1 for r in model if r["cert"] == v)
                            for v in ("", "ship", "product")})
    committed = [r for r in model if r["state"] in ("ship", "product")]
    print("model committed", len(committed),
          "trunc agree", sum(1 for r in committed if r["trunc_order"] == r["state"]),
          "opposite", sum(1 for r in committed if r["trunc_order"] != r["state"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
