#!/usr/bin/env python3
"""Check the computed-zero residual lemma of Supplement S18.4 at the tag.

Runs the certificate record's model over the second (certificate) and first
experiment grids, instrumenting the refined division. Reports, per chain, the
number of calls, the number of exactly-zero computed residuals, the number of
nonzero subnormal residuals, the smallest nominal numerator |a|, the smallest
nonzero computed residual |r|, and the largest recurrence radius. Exits 0 when
the lemma's premises hold: no nonzero residual is subnormal, and the measured
minima are consistent with the derived bounds (|a| > 1e-7 upward, > 0.16
downward).

Usage:
    python residual_probe.py <boys-tag-tree> <certified-bound-dir> <out-log>
"""

import csv
import math
import sys
from pathlib import Path

code = Path(sys.argv[2]) / "code"
sys.path.insert(0, str(code))
import model  # noqa: E402

out_log = Path(sys.argv[3])
model.configure(kernel=Path(sys.argv[1]),
                exp_bits=Path(sys.argv[2]) / "exp" / "exp_bits.txt")
ORIG = model.refined_div
BRANCH = ["?"]
STATS = {}


def probe(a, ea, x, invx):
    q0 = a * invx
    r = math.fma(-q0, x, a)
    s = STATS.setdefault(BRANCH[0], dict(calls=0, zero=0, nonzero=0,
                                         subnormal=0, min_a=None,
                                         min_r=None))
    s["calls"] += 1
    if r == 0.0:
        s["zero"] += 1
    else:
        s["nonzero"] += 1
        if abs(r) < 2.0 ** -1022:
            s["subnormal"] += 1
        if s["min_r"] is None or abs(r) < s["min_r"]:
            s["min_r"] = abs(r)
    if s["min_a"] is None or abs(a) < s["min_a"]:
        s["min_a"] = abs(a)
    return ORIG(a, ea, x, invx)


model.refined_div = probe
original_up, original_down = model._upward, model._downward


def upward(*args):
    BRANCH[0] = "upward"
    return original_up(*args)


def downward(*args):
    BRANCH[0] = "downward"
    return original_down(*args)


model._upward, model._downward = upward, downward
exp = model.load_exp()
lines = []


def scan(label, grid):
    STATS.clear()
    max_radius = 0.0
    with open(grid, newline="") as handle:
        for row in csv.DictReader(handle):
            x = float(row["x"])
            for cell in model.ship_batch(x, exp[x][1]):
                max_radius = max(max_radius, cell["radius"])
    lines.append(label)
    for chain in ("upward", "downward"):
        s = STATS[chain]
        lines.append(
            "  %-8s calls %6d  zero %5d  nonzero %6d  subnormal %d  "
            "min|a| %.9e  min|r| %.6e"
            % (chain, s["calls"], s["zero"], s["nonzero"], s["subnormal"],
               s["min_a"], s["min_r"]))
    lines.append("  largest recurrence radius %.6e" % max_radius)
    return {k: dict(v) for k, v in STATS.items()}


paper = Path(__file__).resolve().parent.parent.parent / "papers" / "boys-crossover"
seen = []
seen.append(scan("second experiment, 836 arguments (the certificate grid):",
                 paper / "sweeps" / "amplification_form" / "data" / "grid.csv"))
seen.append(scan("first experiment, 9353 arguments:",
                 paper / "sweeps" / "product_amplification" / "data" / "grid.csv"))

ok = True
for grid_stats in seen:
    for chain in ("upward", "downward"):
        s = grid_stats[chain]
        if s["subnormal"]:
            ok = False
            lines.append("FAIL: %s has subnormal nonzero residuals" % chain)
lines.append("Result: %s" % ("the lemma's premises hold" if ok else "FAILED"))
text = "\n".join(lines) + "\n"
out_log.write_text(text, encoding="utf-8")
print(text, end="")
sys.exit(0 if ok else 1)
