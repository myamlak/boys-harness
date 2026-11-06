#!/usr/bin/env python3
"""Emit the threshold-free frontier table and the per-order ratio data.

Reads the second experiment's deposited record, beside this file under
``sweeps/amplification_form/results``:

    prediction.csv        one row per band, metric and order at the scored margin
    prediction_scan.csv   the aggregate score over the seven-margin ladder
    per_order.csv         worst error per band, order and method

and writes

    tables/frontier.tex           the frontier table
    figures/per-order-ratio.dat   the per-order ratio series, for pgfplots

The margin rule is the one the analysis used: with ``gap = |p_ship - p_direct|``
and ``sigma`` the recursed bound's accumulated-arithmetic term at the order,
the model abstains when ``gap <= margin * sigma`` and otherwise names the
smaller predicted error.  The aggregate ladder in prediction_scan.csv is
reproduced from prediction.csv as a check; a mismatch exits non-zero and
writes nothing.

The per-order comparison's worst error is the largest, over the orders at
which the model commits, of the named evaluation's per-order worst.  The
better fixed route's worst is the smaller of the two fixed routes' largest
per-order worst over the same orders.

Run without ``python -I``: isolated mode hides the user-site mpmath.
"""

import csv
import math
import os
import sys
from pathlib import Path

if sys.flags.isolated:
    print("SKIP: isolated mode hides user-site mpmath; run without python -I")
    sys.exit(77)

from mpmath import mpf  # noqa: E402

HERE = Path(__file__).resolve().parent
BASE = HERE / "sweeps" / "amplification_form" / "results"
MARGINS = [0.0, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0]
NMAX = 32
BAND_LABEL = {"A1": "B1", "A2": "B2"}


def rows(name):
    with open(BASE / name, newline="") as handle:
        return list(csv.DictReader(handle))


def sci(value, digits=3):
    """A compact LaTeX power of ten with the requested significant digits."""
    if value == 0.0:
        return "0"
    exponent = int(math.floor(math.log10(value)))
    mantissa = value / 10.0**exponent
    text = "%.*f" % (digits - 1, mantissa)
    if float(text) >= 10.0:
        exponent += 1
        text = "%.*f" % (digits - 1, mantissa / 10.0)
    return r"$%s\times10^{%d}$" % (text, exponent)


def state(predicted_ship, predicted_product, arithmetic, margin):
    gap = abs(predicted_ship - predicted_product)
    if gap <= mpf(margin) * arithmetic:
        return None
    return "ship" if predicted_ship < predicted_product else "product"


def main():
    prediction = rows("prediction.csv")
    per_order = rows("per_order.csv")
    scan = rows("prediction_scan.csv")
    if len(prediction) != 132 or len(per_order) != 132:
        raise SystemExit("expected 132 prediction and per-order rows")

    worst = {
        (row["band"], int(row["n"]), row["method"]):
            (float(row["worst_abs"]), float(row["worst_rel"]))
        for row in per_order
    }

    stored = {}
    for row in scan:
        key = float(row["margin"])
        total = stored.setdefault(key, [0, 0, 0, 0])
        total[0] += int(row["commits"])
        total[1] += int(row["right_where_committed"])
        total[2] += int(row["wrong_where_committed"])
        total[3] += int(row["abstains"])

    table_rows = []
    for margin in MARGINS:
        committed = []
        for row in prediction:
            ps = mpf(row["predicted_ship"])
            pd = mpf(row["predicted_product"])
            sigma = mpf(row["arithmetic"])
            choice = state(ps, pd, sigma, margin)
            if choice is None:
                continue
            committed.append((row["band"], row["metric"], int(row["order"]), choice))
        commits = len(committed)
        right = sum(
            1 for band, metric, order, choice in committed
            if choice == next(
                r["measured_winner"] for r in prediction
                if r["band"] == band and r["metric"] == metric
                and int(r["order"]) == order
            )
        )
        wrong = commits - right
        abstains = 132 - commits

        check = [commits, right, wrong, abstains]
        if check != stored[margin]:
            raise SystemExit(
                "margin %s: reproduced %s, stored %s" % (margin, check, stored[margin])
            )

        columns = []
        for metric, index in (("absolute", 0), ("relative", 1)):
            subset = [c for c in committed if c[1] == metric]
            selected = max(worst[(band, order, choice)][index]
                           for band, _, order, choice in subset)
            fixed_ship = max(worst[(band, order, "ship")][index]
                             for band, _, order, _ in subset)
            fixed_product = max(worst[(band, order, "product")][index]
                                for band, _, order, _ in subset)
            better = min(fixed_ship, fixed_product)
            columns.append((selected, better, selected / better))

        table_rows.append((margin, commits, right, wrong, abstains, columns))

    # The per-order ratio series, on the same record.
    series = {}
    for band in ("A1", "A2"):
        for metric, index, label in (
            ("absolute", 0, "abs"),
            ("relative", 1, "rel"),
        ):
            key = "%s%s" % (band, label)
            series[key] = [
                worst[(band, n, "ship")][index] / worst[(band, n, "product")][index]
                for n in range(NMAX + 1)
            ]

    out_tables = HERE / "tables"
    out_figures = HERE / "figures"
    out_tables.mkdir(exist_ok=True)
    out_figures.mkdir(exist_ok=True)

    lines = [
        "% generated by make_frontier.py from",
        "% sweeps/amplification_form/results/{prediction,prediction_scan,per_order}.csv",
        r"\begin{tabular}{lrrrr rr rr rr}",
        r"\toprule",
        r" & & & & & \multicolumn{3}{c}{absolute error} & \multicolumn{3}{c}{relative error}\\",
        r"\cmidrule(lr){6-8}\cmidrule(lr){9-11}",
        r"margin & commits & right & not right & abstains & selected & fixed & ratio & selected & fixed & ratio\\",
        r"\midrule",
    ]
    for margin, commits, right, wrong, abstains, columns in table_rows:
        margin_text = ("%g" % margin)
        cells = ["$%s$" % margin_text if margin else "$0$"]
        cells += [str(commits), str(right), str(wrong), str(abstains)]
        for selected, better, ratio in columns:
            cells += [sci(selected), sci(better), "%.3f" % ratio]
        lines.append(" & ".join(cells) + r"\\")
    lines += [r"\bottomrule", r"\end{tabular}", ""]
    (out_tables / "frontier.tex").write_text("\n".join(lines))

    data = ["n B1abs B1rel B2abs B2rel"]
    for n in range(NMAX + 1):
        data.append("%d %s" % (
            n,
            " ".join("%.6g" % series[key][n]
                    for key in ("A1abs", "A1rel", "A2abs", "A2rel")),
        ))
    (out_figures / "per-order-ratio.dat").write_text("\n".join(data) + "\n")

    print("wrote tables/frontier.tex and figures/per-order-ratio.dat")
    for margin, commits, right, wrong, abstains, columns in table_rows:
        text = ["margin %-4g commits %3d right %3d not right %2d abstains %3d"
                % (margin, commits, right, wrong, abstains)]
        for metric, (selected, better, ratio) in zip(("abs", "rel"), columns):
            text.append("  %s selected %.6g fixed %.6g ratio %.4f"
                        % (metric, selected, better, ratio))
        print("\n".join(text))


if __name__ == "__main__":
    main()
