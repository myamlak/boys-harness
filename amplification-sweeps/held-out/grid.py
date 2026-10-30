#!/usr/bin/env python3
"""The held-out sweep's argument list, rule H.

This is the third Paper A argument set and it is disjoint from both earlier
sweeps by construction and by an exact-double assertion.  Neither the
committed gate grid (shared by both earlier sweeps) nor any set of the first
sweep (U, C, L, R, G, T, E) appears here.  Five rules, applied per band:

    H1  shifted uniform, 300 points per band:
            f = (3k + 1) / 900,          k = 0 .. 299
    H2  Chebyshev-Gauss, 136 points per band (a multiple of four, so the
        centre is not a node; both earlier sweeps contain the centre):
            f = (1 + cos(pi (j + 1/2) / 136)) / 2,   j = 0 .. 135
    H3  base-2 edge offsets, 27 per edge, both edges:
            f = 2^-m and f = 1 - 2^-m,   m = 12 .. 38
    H4  upper-half window, 96 points per band, two-thirds-step offset:
            f = 1/2 + (3j + 2) / 576,    j = 0 .. 95
    H5  lower-quarter window, 64 points per band, two-thirds-step offset:
            f = (3j + 2) / 768,          j = 0 .. 63

    x = left + (right - left) f

with left and right the band's own edges read from the pinned tree
(kRegionA1Edge, kRegionAEnd).  The A1 band is [0, edge), the A2 band
[edge, end).  650 points per band, 1300 arguments; no rule needs the gate
grid, so the held-out set is independent of the committed gate arguments.

Why the rules are disjoint from the earlier sets, before the numeric check:

  * H1, H4, H5 have numerators that are never divisible by 3 (3k+1, 3j+290
    and 3j+2 are 1, 2 and 2 mod 3), so the reduced denominator keeps a factor
    3 or 5 that no dyadic rational has: none can equal a U point (k/2048) or
    an H3 point (a power of two).
  * H2 nodes have rational cosine only at cos = 0, +-1/2, +-1 (Niven); with
    the angle (2j+1) pi / 272 that needs 2j + 1 in {272/3, 136, 544/3}, none
    of them an odd integer, and every other family's fraction is rational.
  * H3 offsets 2^-m with m >= 12 are closer to the edge than U's first step
    1/2048 = 2^-11 and are not equal to it, and a power of two is not a
    power of ten (the first sweep's L and R).  m stops at 38: at the region
    end, the m = 39 offset lands on one of the first sweep's R arguments
    (11.899848152097661), so the stated rule stops one octave short of it.

main() asserts every one of these against the two earlier grids by exact
double equality and fails the run if a single argument coincides.

Usage:
    python grid.py --boys-tree C:\\Users\\myaml\\dev\\_boys-v511 --out data/grid.csv
"""

import argparse
import csv
import math
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "common"))
import boysprovenance  # noqa: E402

DEFAULT_TREE = boysprovenance.DEFAULT_TREE
EARLIER_GRIDS = [
    Path(__file__).resolve().parents[1] / "product_amplification/data/grid.csv",
    Path(__file__).resolve().parents[1] / "amplification_form/data/grid.csv",
]


def add(store, fraction, tag):
    value = store.get(fraction)
    if value is None:
        store[fraction] = [tag]
    elif tag not in value:
        value.append(tag)


def band_fractions():
    points = {}

    for k in range(300):
        add(points, (3 * k + 1) / 900.0, "H1")

    for j in range(136):
        add(points, 0.5 * (1.0 + math.cos(math.pi * (j + 0.5) / 136.0)), "H2")

    for m in range(12, 39):
        add(points, 2.0**-m, "H3")
        add(points, 1.0 - 2.0**-m, "H3")

    for j in range(96):
        add(points, 0.5 + (3 * j + 2) / 576.0, "H4")

    for j in range(64):
        add(points, (3 * j + 2) / 768.0, "H5")

    return points


def earlier_arguments():
    seen = {}
    for path in EARLIER_GRIDS:
        if not path.is_file():
            raise SystemExit("the earlier grid is missing, cannot assert disjointness: %s" % path)
        with open(path, newline="") as handle:
            for row in csv.DictReader(handle):
                seen.setdefault(float(row["x"]), set()).add(path.parent.parent.name)
    return seen


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--boys-tree", default=os.environ.get("BOYS_TREE", DEFAULT_TREE))
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    tree = Path(args.boys_tree).resolve()
    constants = boysprovenance.parse_constants(tree)
    edge = constants["kRegionA1Edge"]
    end = constants["kRegionAEnd"]

    fractions = band_fractions()
    earlier = earlier_arguments()

    rows = []
    collisions = []
    per_band = {}
    for band, left, right in (("A1", 0.0, edge), ("A2", edge, end)):
        points = {}
        for fraction, tags in fractions.items():
            value = left + (right - left) * fraction
            if not math.isfinite(value) or not (left <= value < right):
                continue
            add(points, value, "+".join(tags))
        per_band[band] = len(points)
        for value in sorted(points):
            if value in earlier:
                collisions.append((band, value))
            rows.append((band, value, points[value]))

    if collisions:
        raise SystemExit(
            "the held-out set is not disjoint from the earlier sweeps: %s"
            % ", ".join("%s %r" % collision for collision in collisions[:10])
        )

    for band in ("A1", "A2"):
        print("%s: %d arguments" % (band, per_band[band]))
    print(
        "disjointness: %d arguments, none equal to any argument of the two earlier grids "
        "(%d distinct earlier arguments checked)" % (len(rows), len(earlier))
    )

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["band", "x", "sets"])
        for band, value, tags in rows:
            writer.writerow([band, repr(value), tags])
    print("wrote %s: %d arguments (A1 %d, A2 %d)" % (out, len(rows), per_band["A1"], per_band["A2"]))


if __name__ == "__main__":
    main()
