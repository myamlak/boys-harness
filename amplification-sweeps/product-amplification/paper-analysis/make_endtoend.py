#!/usr/bin/env python3
"""Emit the end-to-end fallback table (Table 8) of the main text.

Reads the second experiment's deposited record, beside this file under
``sweeps/amplification_form/results``:

    prediction.csv   one row per band, metric and order at the scored margin
    per_order.csv    worst error per band, order and method

Three rules are scored over all 132 band-metric-order comparisons:

    always-direct      the direct (matrix-product) route at every comparison
    always-recursed    the recursed route at every comparison
    model+fallback     the route the model names at a commitment, and the
                       named fallback at an abstention

The model is scored at its own margin (one times the recursed bound's
accumulated-arithmetic term), the rung of Table 7 that commits at 65 of the
132 comparisons and is right at every one.  Both fallbacks are reported,
because the answer depends on which one is named.

A choice's error at a comparison is the named evaluation's worst error in
that metric over that band at that order; a choice's worst is the maximum over
the 66 comparisons of one metric.  Reduction is (fixed - model) / fixed, so a
positive value is an improvement on that fixed choice.

Writes

    tables/endtoend.tex

Run without ``python -I``: isolated mode hides the user-site mpmath.
"""

import csv
import math
import sys
from pathlib import Path

if sys.flags.isolated:
    print("SKIP: isolated mode hides user-site mpmath; run without python -I")
    sys.exit(77)

from mpmath import mpf  # noqa: E402

HERE = Path(__file__).resolve().parent
BASE = HERE / "sweeps" / "amplification_form" / "results"
MARGIN = 1.0
METRICS = (("absolute", 0), ("relative", 1))
ROUTES = {"ship": "recursed", "product": "direct"}


def rows(name):
    with open(BASE / name, newline="") as handle:
        return list(csv.DictReader(handle))


def state(predicted_ship, predicted_product, arithmetic, margin):
    gap = abs(predicted_ship - predicted_product)
    if gap <= mpf(margin) * arithmetic:
        return None
    return "ship" if predicted_ship < predicted_product else "product"


def sci(value):
    """Three significant digits as a LaTeX power of ten."""
    exponent = int(math.floor(math.log10(value)))
    mantissa = value / 10.0**exponent
    text = "%.2f" % mantissa
    if float(text) >= 10.0:
        exponent += 1
        text = "%.2f" % (mantissa / 10.0)
    return r"$%s\times10^{%d}$" % (text, exponent)


def pct(value):
    """Three significant digits, or four near 100 per cent."""
    value = 100.0 * value
    if abs(value) < 5e-3:
        return "0.00"
    if abs(value) >= 99.9:
        return "%.2f" % value
    return "%.1f" % value if abs(value) >= 10 else "%.2f" % value


def main():
    prediction = rows("prediction.csv")
    per_order = rows("per_order.csv")
    if len(prediction) != 132 or len(per_order) != 132:
        raise SystemExit("expected 132 prediction and per-order rows")

    worst = {
        (row["band"], int(row["n"]), row["method"]):
            (float(row["worst_abs"]), float(row["worst_rel"]))
        for row in per_order
    }

    choice = {}
    for row in prediction:
        key = (row["band"], row["metric"], int(row["order"]))
        choice[key] = state(
            mpf(row["predicted_ship"]),
            mpf(row["predicted_product"]),
            mpf(row["arithmetic"]),
            MARGIN,
        )

    committed = sum(1 for value in choice.values() if value is not None)
    per_metric_commits = {
        metric: sum(1 for (band, m, order), value in choice.items()
                    if m == metric and value is not None)
        for metric, _ in METRICS
    }
    if committed != 65:
        raise SystemExit("expected 65 commitments at margin 1, got %d" % committed)

    table = {}
    for metric, index in METRICS:
        keys = [k for k in choice if k[1] == metric]

        # A fixed choice names its route at every comparison; the model names
        # its own route at a commitment and the fallback at an abstention.
        def fixed(route):
            return max(worst[(key[0], key[2], route)][index] for key in keys)

        def fixed_at(route):
            return max(keys, key=lambda k: worst[(k[0], k[2], route)][index])

        def model_with(fallback):
            total = 0.0
            for key in keys:
                band, _, order = key
                pick = choice[key] or fallback
                total = max(total, worst[(band, order, pick)][index])
            return total

        def model_at(fallback):
            return max(((key, choice[key] or fallback) for key in keys),
                       key=lambda kc: worst[(kc[0][0], kc[0][2], kc[1])][index])

        direct = fixed("product")
        recursed = fixed("ship")
        model_direct = model_with("product")
        model_recursed = model_with("ship")
        table[metric] = {
            "commits": per_metric_commits[metric],
            "abstains": len(keys) - per_metric_commits[metric],
            "direct": direct,
            "recursed": recursed,
            "model_direct": model_direct,
            "model_recursed": model_recursed,
            "direct_at": fixed_at("product"),
            "recursed_at": fixed_at("ship"),
            "model_direct_at": model_at("product"),
            "model_recursed_at": model_at("ship"),
        }

    lines = [
        "% generated by make_endtoend.py from",
        "% sweeps/amplification_form/results/{prediction,per_order}.csv",
        r"\begin{tabular}{llrrr}",
        r"\toprule",
        r" & choice & worst error & vs always-direct & vs always-recursed\\",
        r"\midrule",
    ]
    for metric, _ in METRICS:
        entry = table[metric]
        lines.append(
            r"\multicolumn{5}{@{}l}{\emph{%s error, %d comparisons: "
            r"%d commits, %d abstentions}}\\"
            % (metric, entry["commits"] + entry["abstains"],
               entry["commits"], entry["abstains"]))
        for label, key in (
            ("always-direct", "direct"),
            ("always-recursed", "recursed"),
            ("model, direct fallback", "model_direct"),
            ("model, recursed fallback", "model_recursed"),
        ):
            value = entry[key]
            if key in ("direct", "recursed"):
                vs_direct = vs_recursed = "--"
            else:
                vs_direct = pct((entry["direct"] - value) / entry["direct"])
                vs_recursed = pct((entry["recursed"] - value) / entry["recursed"])
                vs_direct = ("$%s$" % vs_direct) if vs_direct != "0.0" else "0.0"
                vs_recursed = ("$%s$" % vs_recursed) if vs_recursed != "0.0" else "0.0"
            lines.append(r" & %s & %s & %s & %s\\"
                         % (label, sci(value), vs_direct, vs_recursed))
        if metric != METRICS[-1][0]:
            lines.append(r"\midrule")
    lines += [r"\bottomrule", r"\end{tabular}", ""]
    (HERE / "tables" / "endtoend.tex").write_text("\n".join(lines))

    for metric, _ in METRICS:
        entry = table[metric]
        print(metric, {k: ("%.6g" % v if isinstance(v, float) else v)
                       for k, v in entry.items()})
        print("  reduction vs direct: %.2f%%"
              % (100.0 * (entry["direct"] - entry["model_direct"]) / entry["direct"]))
        print("  reduction vs recursed: %.2f%%"
              % (100.0 * (entry["recursed"] - entry["model_direct"]) / entry["recursed"]))
        print("  direct fallback worst %.6g  (always-direct %.6g)"
              % (entry["model_direct"], entry["direct"]))
        print("  recursed fallback worst %.6g" % entry["model_recursed"])
        print("  always-direct worst at %s" % (entry["direct_at"],))
        print("  always-recursed worst at %s" % (entry["recursed_at"],))
        print("  model(direct) worst at %s" % (entry["model_direct_at"],))
        print("  model(recursed) worst at %s" % (entry["model_recursed_at"],))
    print("wrote tables/endtoend.tex")


if __name__ == "__main__":
    main()
