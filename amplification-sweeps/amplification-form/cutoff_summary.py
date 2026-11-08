#!/usr/bin/env python3
"""Regenerate cutoff-summary.csv from the deposited per-cell amplification values.

The resolution rule is the second sweep's own (`amplification_form/analyze.py`):
an argument is resolved when all 33 transfer magnitudes are measurable and the
nearer neighbouring order differs from the minimum by more than the separation
cutoff; two orders whose values agree within the tie tolerance are a tie, and a
tie is a match when the discrete formula's order is one of the two. The record's
bands are A1/A2; the summary uses the paper's B1/B2 labels.

Usage:
    python cutoff_summary.py
"""

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
CUTOFFS = (1e-12, 1e-9)
NMAX = 32
BAND_LABEL = {"A1": "B1", "A2": "B2"}


def read_cells(path):
    values = defaultdict(dict)
    with open(path, newline="") as handle:
        for row in csv.DictReader(handle):
            values[(row["band"], float(row["x"]))][int(row["n"])] = (
                float(row["transfer"]),
                row["measurable"] == "1",
            )
    return values


def counts_at(values, cutoff, tie_tolerance):
    counts = {band: [0, 0] for band in ("A1", "A2")}
    for (band, x), cells in values.items():
        if len(cells) != NMAX + 1:
            continue
        if any(not measured or not math.isfinite(transfer)
               for transfer, measured in cells.values()):
            continue
        magnitudes = [abs(cells[n][0]) for n in range(NMAX + 1)]
        m = min(range(NMAX + 1), key=lambda n: magnitudes[n])
        if magnitudes[m] == 0:
            continue
        nmin = min(NMAX, max(0, math.ceil(x - 0.5)))
        gaps = []
        if m > 0:
            gaps.append(abs(magnitudes[m - 1] - magnitudes[m]) / magnitudes[m])
        if m < NMAX:
            gaps.append(abs(magnitudes[m + 1] - magnitudes[m]) / magnitudes[m])
        gap = min(gaps)
        tie_low = m > 0 and gaps[0] <= tie_tolerance
        tie_high = m < NMAX and abs(magnitudes[m + 1] - magnitudes[m]) / magnitudes[m] <= tie_tolerance
        tie = tie_low or tie_high
        tied = [m - 1, m] if tie_low else ([m, m + 1] if tie_high else [m])
        if gap <= cutoff and not tie:
            continue
        counts[band][0] += 1
        if m == nmin or (tie and nmin in tied):
            counts[band][1] += 1
    return counts


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cells", default=str(HERE / "results" / "transfer_cells.csv"))
    parser.add_argument("--out", default=str(HERE / "cutoff-summary.csv"))
    parser.add_argument("--tie-tolerance", type=float, default=1e-12)
    args = parser.parse_args()

    values = read_cells(args.cells)
    with open(args.out, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["cutoff", "resolved", "matching", "resolved_B1",
                         "resolved_B2", "matching_B1", "matching_B2"])
        for cutoff in CUTOFFS:
            counts = counts_at(values, cutoff, args.tie_tolerance)
            resolved = counts["A1"][0] + counts["A2"][0]
            matching = counts["A1"][1] + counts["A2"][1]
            writer.writerow(["%.0e" % cutoff, resolved, matching,
                             counts["A1"][0], counts["A2"][0],
                             counts["A1"][1], counts["A2"][1]])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
