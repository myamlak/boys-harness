#!/usr/bin/env python3
"""Print the measured route difference in the chemistry bridge, by band region.

The deposited argument streams carry, for every distinct (nmax, x), the engine
value and both routes' values at each order 0..nmax. This script reports the
largest |class-row - paper-a| difference over all arguments, over band B2, and
over the arguments above B2, for the three bases. The B2 upper edge is
11.899848152108484.

Usage:
    python above_b2_max.py
"""

import lzma
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
OUTPUTS = HERE.parent / "outputs"
B2_LOW = 5.949924076054242
B2_HIGH = 11.899848152108484
FILES = (
    ("STO-3G", OUTPUTS / "default.args.txt"),
    ("6-31G*", OUTPUTS / "631gstar-default.args.txt"),
    ("cc-pVTZ", OUTPUTS / "ccpvtz-default.args.txt.xz"),
)


def rows(path):
    if path.suffix == ".xz":
        handle = lzma.open(path, "rt", encoding="utf-8")
    else:
        handle = open(path, "r", encoding="utf-8")
    with handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            yield line.split()


def value(text):
    return float(text.replace("nan(ind)", "nan"))


def maximum(path, region):
    best = (0.0, None, None, None)
    for fields in rows(path):
        nmax = int(fields[0])
        count = int(fields[1])
        x = value(fields[2])
        width = nmax + 1
        class_row = [value(fields[3 + width + n]) for n in range(width)]
        paper_a = [value(fields[3 + 2 * width + n]) for n in range(width)]
        for n, (a, b) in enumerate(zip(class_row, paper_a)):
            if not (math.isfinite(a) and math.isfinite(b)):
                continue
            if region == "B2" and not (B2_LOW <= x < B2_HIGH):
                continue
            if region == "above-B2" and x < B2_HIGH:
                continue
            difference = abs(a - b)
            if difference > best[0]:
                best = (difference, n, x, count)
    return best


def main():
    print("basis              region      max |class-row - paper-a|  order  x  count")
    overall = {}
    for region in ("all", "B2", "above-B2"):
        overall[region] = (0.0, None, None, None, None)
        for name, path in FILES:
            best = maximum(path, region)
            print("%-18s %-10s %.17g  %d  %.17g  %d"
                  % (name, region, best[0], best[1], best[2], best[3]))
            if best[0] > overall[region][0]:
                overall[region] = (best[0], name, best[1], best[2], best[3])
    print()
    for region, best in overall.items():
        print("all bases %-10s %.17g  %s order %d x=%.17g count %d"
              % (region, best[0], best[1], best[2], best[3], best[4]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
