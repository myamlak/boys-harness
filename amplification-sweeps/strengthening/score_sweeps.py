"""Score the crossing model on both deposited sweeps.

The model arithmetic is imported from the deposit's analyse script
(amplification-form/analyze.py) rather than retyped, so the first sweep's
score is computed by the same code that produced the second sweep's deposited
prediction.csv and prediction_scan.csv.

Run from this directory after prepare_inputs.py:
    python score_sweeps.py --sweep second --out out/second --margin 1
    python score_sweeps.py --sweep first  --out out/first  --margin 1

Every output is written under out/; no deposited record is touched.  The first
sweep's cells and reference are read from inputs/, decompressed from the
deposit's product-amplification streams; the second sweep is read in place
from amplification-form/.
"""

import argparse
import csv
import importlib.util
import json
import math
import os
import sys
from pathlib import Path

from mpmath import mp, mpf

HERE = Path(__file__).resolve().parent
DEPOSIT = HERE.parent
NMAX = 32
MARGINS = [0.0, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0]


def load_analysis():
    """Import the deposited analyse module without running main()."""
    spec = importlib.util.spec_from_file_location(
        "paper_a_analyze", DEPOSIT / "amplification-form" / "analyze.py"
    )
    module = importlib.util.module_from_spec(spec)
    sys.path.insert(0, str(DEPOSIT / "common"))
    sys.path.insert(0, str(DEPOSIT / "amplification-form"))
    spec.loader.exec_module(module)
    return module


def read_cells(path):
    rows = []
    with open(path, newline="") as handle:
        for row in csv.DictReader(handle):
            rows.append(
                (
                    row["band"],
                    row["x"],
                    float(row["x"]),
                    int(row["n"]),
                    float(row["recursed"]),
                    float(row["direct"]),
                )
            )
    return rows


def update(store, key, value, x_text):
    current = store.get(key)
    if current is None or value > current[0]:
        store[key] = (value, x_text)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sweep", choices=["first", "second"], required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--margin", type=float, default=1.0)
    args = parser.parse_args()

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    if args.sweep == "second":
        cells_path = DEPOSIT / "amplification-form" / "results" / "cells.csv"
        reference_path = DEPOSIT / "amplification-form" / "data" / "reference.csv"
        provenance_path = DEPOSIT / "amplification-form" / "results" / "provenance.json"
        deposited_per_order = DEPOSIT / "amplification-form" / "results" / "per_order.csv"
        deposited_band = DEPOSIT / "amplification-form" / "results" / "band_summary.csv"
        deposited_prediction = DEPOSIT / "amplification-form" / "results" / "prediction.csv"
        deposited_scan = DEPOSIT / "amplification-form" / "results" / "prediction_scan.csv"
    else:
        cells_path = HERE / "inputs" / "cells.csv"
        reference_path = HERE / "inputs" / "reference.csv"
        provenance_path = DEPOSIT / "product-amplification" / "results" / "provenance.json"
        deposited_per_order = DEPOSIT / "product-amplification" / "results" / "per_order.csv"
        deposited_band = DEPOSIT / "product-amplification" / "results" / "band_summary.csv"
        deposited_prediction = None
        deposited_scan = None

    module = load_analysis()
    mp.dps = 70
    reference = module.read_reference(str(reference_path))
    cells = read_cells(cells_path)
    constants = json.load(open(provenance_path, encoding="utf-8"))["constants"]
    thresholds = [float(v) for v in constants["kTierThresholds"]]
    model_from = thresholds[0]

    print("sweep:", args.sweep)
    print("cells:", len(cells), "reference arguments:", len(reference))
    print("model_from:", repr(model_from))

    # ---- measured per-order worsts and band winners (the paper's winner) --
    ship_abs, ship_rel, prod_abs, prod_rel = {}, {}, {}, {}
    band_ship = {"A1": {}, "A2": {}}
    band_prod = {"A1": {}, "A2": {}}
    error0 = {}
    winner = {}
    per_order = {}

    for band, x_text, x, n, recursed, direct in cells:
        values = reference[(band, x)]
        expected = values[n]
        rec_error = abs(mpf(recursed) - expected)
        dir_error = abs(mpf(direct) - expected)
        rec_relative = rec_error / expected if expected != 0 else mpf("inf")
        dir_relative = dir_error / expected if expected != 0 else mpf("inf")

        update(ship_abs, (band, n), rec_error, x_text)
        update(prod_abs, (band, n), dir_error, x_text)
        update(ship_rel, (band, n), rec_relative, x_text)
        update(prod_rel, (band, n), dir_relative, x_text)

        for metric, value_now in (("absolute", rec_error), ("relative", rec_relative)):
            entry = band_ship[band].get(metric)
            if entry is None or value_now > entry[0]:
                band_ship[band][metric] = (value_now, n, x_text)
        for metric, value_now in (("absolute", dir_error), ("relative", dir_relative)):
            entry = band_prod[band].get(metric)
            if entry is None or value_now > entry[0]:
                band_prod[band][metric] = (value_now, n, x_text)

        if n == 0:
            error0[(band, x)] = rec_error

    for band in ("A1", "A2"):
        for n in range(NMAX + 1):
            for metric, ship_store, prod_store in (
                ("absolute", ship_abs, prod_abs),
                ("relative", ship_rel, prod_rel),
            ):
                s = ship_store[(band, n)][0]
                p = prod_store[(band, n)][0]
                if s < p:
                    state = "ship"
                elif p < s:
                    state = "product"
                else:
                    state = "tie"
                winner[(band, metric, n)] = state
                per_order[(band, n, metric)] = (s, p, ship_store[(band, n)][1], prod_store[(band, n)][1])

    # ---- validate the per-order record ------------------------------------
    mismatches = 0
    checked = 0
    if deposited_per_order.is_file():
        with open(deposited_per_order, newline="") as handle:
            for row in csv.DictReader(handle):
                band, n = row["band"], int(row["n"])
                for metric, col, at_col in (
                    ("absolute", "worst_abs", "at_x_abs"),
                    ("relative", "worst_rel", "at_x_rel"),
                ):
                    value = float(row[col])
                    mine, at_x = per_order[(band, n, metric)][0 if row["method"] == "ship" else 1], \
                        per_order[(band, n, metric)][2 if row["method"] == "ship" else 3]
                    checked += 1
                    if float(at_x) != float(row[at_col]) or abs(value - float(mine)) / value > 1e-12:
                        mismatches += 1
    print("per_order reproduction: %d values checked, %d mismatched" % (checked, mismatches))

    band_mismatches = 0
    if deposited_band.is_file():
        with open(deposited_band, newline="") as handle:
            for row in csv.DictReader(handle):
                band, metric = row["band"], row["metric"]
                s = band_ship[band][metric]
                p = band_prod[band][metric]
                ok = (
                    abs(float(row["ship_worst"]) - float(s[0])) / float(s[0]) < 1e-12
                    and int(row["ship_at_n"]) == s[1]
                    and abs(float(row["prod_worst"]) - float(p[0])) / float(p[0]) < 1e-12
                    and int(row["prod_at_n"]) == p[1]
                )
                if not ok:
                    band_mismatches += 1
    print("band_summary reproduction: %d mismatched rows" % band_mismatches)

    # ---- the model's band-level inputs ------------------------------------
    predicted_rec_abs = {"A1": [mpf(0)] * (NMAX + 1), "A2": [mpf(0)] * (NMAX + 1)}
    predicted_rec_rel = {"A1": [mpf(0)] * (NMAX + 1), "A2": [mpf(0)] * (NMAX + 1)}
    predicted_dir_abs = {"A1": [mpf(0)] * (NMAX + 1), "A2": [mpf(0)] * (NMAX + 1)}
    predicted_dir_rel = {"A1": [mpf(0)] * (NMAX + 1), "A2": [mpf(0)] * (NMAX + 1)}
    predicted_arith_abs = {"A1": [mpf(0)] * (NMAX + 1), "A2": [mpf(0)] * (NMAX + 1)}
    predicted_arith_rel = {"A1": [mpf(0)] * (NMAX + 1), "A2": [mpf(0)] * (NMAX + 1)}
    model_arguments = {"A1": 0, "A2": 0}
    population_cells = {"A1": 0, "A2": 0}
    model_seen = set()
    amp_cache = {}

    def amp(n, x):
        key = (n, x)
        if key not in amp_cache:
            amp_cache[key] = module.amplification(n, x)
        return amp_cache[key]

    for band, x_text, x, n, recursed, direct in cells:
        if x < model_from:
            continue
        if (band, x) not in model_seen:
            model_seen.add((band, x))
            model_arguments[band] += 1
        seed = error0.get((band, x))
        if seed is None:
            continue
        population_cells[band] += 1
        expected = reference[(band, x)][n]
        dir_abs = abs(mpf(direct) - expected)
        dir_rel = dir_abs / expected if expected != 0 else mpf("inf")

        if n == 0 or x >= thresholds[n]:
            arith_abs = module.recurrence_arithmetic(n, mpf(x), reference[(band, x)]) if n else mpf(0)
            rec_abs = seed * amp(n, mpf(x)) + arith_abs
            rec_rel = rec_abs / expected if expected != 0 else mpf("inf")
            arith_rel = arith_abs / expected if expected != 0 else mpf("inf")
        else:
            arith_abs = mpf(0)
            arith_rel = mpf(0)
            rec_abs = dir_abs
            rec_rel = dir_rel

        if rec_abs > predicted_rec_abs[band][n]:
            predicted_rec_abs[band][n] = rec_abs
        if rec_rel > predicted_rec_rel[band][n]:
            predicted_rec_rel[band][n] = rec_rel
        if dir_abs > predicted_dir_abs[band][n]:
            predicted_dir_abs[band][n] = dir_abs
        if dir_rel > predicted_dir_rel[band][n]:
            predicted_dir_rel[band][n] = dir_rel
        if arith_abs > predicted_arith_abs[band][n]:
            predicted_arith_abs[band][n] = arith_abs
        if arith_rel > predicted_arith_rel[band][n]:
            predicted_arith_rel[band][n] = arith_rel

    print("model population arguments:", model_arguments, "cells:", population_cells)

    def prediction_rows(band, metric, margin):
        rows = []
        if metric == "absolute":
            ship, product, sigma = predicted_rec_abs[band], predicted_dir_abs[band], predicted_arith_abs[band]
        else:
            ship, product, sigma = predicted_rec_rel[band], predicted_dir_rel[band], predicted_arith_rel[band]
        for n in range(NMAX + 1):
            ps, pd = ship[n], product[n]
            gap = abs(ps - pd)
            if gap <= mpf(margin) * sigma[n]:
                state, choice = "abstain", ""
            else:
                state = "ship" if ps < pd else "product"
                choice = state
            measured = winner[(band, metric, n)]
            right = choice == measured if state != "abstain" else False
            rows.append((n, ps, pd, sigma[n], state, measured, choice, right))
        return rows

    # margin-1 prediction record
    with open(out / "prediction.csv", "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["band", "metric", "order", "predicted_rec", "predicted_dir", "arithmetic", "state", "measured_winner", "choice", "right_where_committed"])
        for band in ("A1", "A2"):
            for metric in ("absolute", "relative"):
                for n, ps, pd, sigma, state, measured, choice, right in prediction_rows(band, metric, args.margin):
                    writer.writerow([band, metric, n, mp.nstr(ps, 17), mp.nstr(pd, 17), mp.nstr(sigma, 17), state, measured, choice, int(right)])

    # scan ladder, band-metric order comparisons
    scan_rows = []
    with open(out / "scan.csv", "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["margin", "band", "metric", "commits", "right", "wrong", "abstains", "ties_measured", "orders"])
        for margin in MARGINS:
            for band in ("A1", "A2"):
                for metric in ("absolute", "relative"):
                    rows = prediction_rows(band, metric, margin)
                    commits = sum(1 for r in rows if r[4] != "abstain")
                    right = sum(1 for r in rows if r[4] != "abstain" and r[7])
                    abstains = len(rows) - commits
                    ties = sum(1 for r in rows if r[4] != "abstain" and r[5] == "tie")
                    writer.writerow([margin, band, metric, commits, right, commits - right, abstains, ties, len(rows)])
                    scan_rows.append((margin, band, metric, commits, right, commits - right, abstains))

    # aggregate row per margin
    with open(out / "scan_all.csv", "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["margin", "commits", "right", "wrong", "abstains", "of_which_measured_ties"])
        for margin in MARGINS:
            rows = []
            for band in ("A1", "A2"):
                for metric in ("absolute", "relative"):
                    rows.extend(prediction_rows(band, metric, margin))
            commits = sum(1 for r in rows if r[4] != "abstain")
            right = sum(1 for r in rows if r[4] != "abstain" and r[7])
            abstains = len(rows) - commits
            ties = sum(1 for r in rows if r[4] != "abstain" and r[5] == "tie")
            writer.writerow([margin, commits, right, commits - right, abstains, ties])

    # validate against the deposited prediction.csv and prediction_scan.csv
    if deposited_prediction is not None:
        mismatch = 0
        count = 0
        with open(deposited_prediction, newline="") as handle:
            for row in csv.DictReader(handle):
                band, metric, n = row["band"], row["metric"], int(row["order"])
                rows = prediction_rows(band, metric, args.margin)
                state, choice, right = rows[n][4], rows[n][6], rows[n][7]
                count += 1
                if state != row["state"] or choice != (row["choice"] or "") or int(row["right_where_committed"]) != int(right):
                    mismatch += 1
        print("deposited prediction.csv: %d rows, %d mismatches on state/choice/right" % (count, mismatch))
    if deposited_scan is not None:
        mismatch = 0
        count = 0
        with open(deposited_scan, newline="") as handle:
            for row in csv.DictReader(handle):
                margin = float(row["margin"])
                band, metric = row["band"], row["metric"]
                rows = prediction_rows(band, metric, margin)
                commits = sum(1 for r in rows if r[4] != "abstain")
                right = sum(1 for r in rows if r[4] != "abstain" and r[7])
                abstains = len(rows) - commits
                count += 1
                if (commits, right, commits - right, abstains) != (
                    int(row["commits"]), int(row["right_where_committed"]),
                    int(row["wrong_where_committed"]), int(row["abstains"]),
                ):
                    mismatch += 1
        print("deposited prediction_scan.csv: %d rows, %d mismatches" % (count, mismatch))

    # ---- cell-level scoring ----------------------------------------------
    # The paper's model is a band-level statement.  This block applies the
    # same per-cell formulae and the same abstention rule at each measured
    # (argument, order) cell, which is a strictly finer test; it is reported
    # separately and is not the paper's scoring unit.
    cell_rows = []
    cell_amp = {}

    def camp(n, x):
        key = (n, x)
        if key not in cell_amp:
            cell_amp[key] = module.amplification(n, x)
        return cell_amp[key]

    for band, x_text, x, n, recursed, direct in cells:
        if x < model_from:
            continue
        seed = error0.get((band, x))
        if seed is None:
            continue
        expected = reference[(band, x)][n]
        rec_error = abs(mpf(recursed) - expected)
        dir_abs = abs(mpf(direct) - expected)
        dir_rel = dir_abs / expected if expected != 0 else mpf("inf")
        measured = "tie" if rec_error == dir_abs else ("ship" if rec_error < dir_abs else "product")
        measured_rel = "tie" if rec_error == dir_rel else ("ship" if rec_error < dir_rel else "product")
        if n == 0 or x >= thresholds[n]:
            arith_abs = module.recurrence_arithmetic(n, mpf(x), reference[(band, x)]) if n else mpf(0)
            rec_abs = seed * camp(n, mpf(x)) + arith_abs
            arith_rel = arith_abs / expected if expected != 0 else mpf("inf")
        else:
            arith_abs = mpf(0)
            arith_rel = mpf(0)
            rec_abs = dir_abs
        rec_rel = rec_abs / expected if expected != 0 else mpf("inf")
        cell_rows.append((band, x_text, n, rec_abs, dir_abs, arith_abs, measured,
                          rec_rel, dir_rel, arith_rel, measured_rel))

    with open(out / "cell_scan.csv", "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["margin", "metric", "commits", "right", "wrong", "ties", "abstains", "cells",
                         "commits_A1", "right_A1", "commits_A2", "right_A2"])
        for margin in MARGINS:
            for metric in ("absolute", "relative"):
                commits = right = ties = 0
                per_band = {"A1": [0, 0], "A2": [0, 0]}
                for band, x_text, n, ra, da, aa, meas_abs, rr, dr, ar, meas_rel in cell_rows:
                    if metric == "absolute":
                        rec_v, dir_v, arith_v, meas = ra, da, aa, meas_abs
                    else:
                        rec_v, dir_v, arith_v, meas = rr, dr, ar, meas_rel
                    gap = abs(rec_v - dir_v)
                    if gap <= mpf(margin) * arith_v:
                        continue
                    commits += 1
                    per_band[band][0] += 1
                    if rec_v < dir_v:
                        choice = "ship"
                    elif dir_v < rec_v:
                        choice = "product"
                    else:
                        choice = "tie"
                    if choice == meas:
                        right += 1
                        per_band[band][1] += 1
                    elif meas == "tie":
                        ties += 1
                writer.writerow([margin, metric, commits, right, commits - right - ties, ties,
                                 len(cell_rows) - commits, len(cell_rows),
                                 per_band["A1"][0], per_band["A1"][1],
                                 per_band["A2"][0], per_band["A2"][1]])

    # ---- per-order ratio table (item 2 helper) ----------------------------
    with open(out / "per_order_ratio.csv", "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["band", "n", "metric", "rec_worst", "dir_worst", "rec_over_dir", "measured_winner"])
        for band in ("A1", "A2"):
            for n in range(NMAX + 1):
                for metric in ("absolute", "relative"):
                    s, p = per_order[(band, n, metric)][0], per_order[(band, n, metric)][1]
                    writer.writerow([band, n, metric, mp.nstr(s, 17), mp.nstr(p, 17),
                                     mp.nstr(s / p, 17), winner[(band, metric, n)]])

    print("wrote", out)


if __name__ == "__main__":
    main()
