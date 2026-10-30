#!/usr/bin/env python3
"""Held-out test: the frozen prediction rule on a disjoint argument set.

The decision rule and its abstention margin are taken from Paper A unchanged.
The script proves that before it scores anything, in two steps:

  1. the model below is a line-for-line copy of the one in
     ../amplification_form/analyze.py; applied to the second sweep's own cells
     and reference it must reproduce ../amplification_form/results/
     prediction_record.csv field for field (132 rows, 17-digit values);
  2. the same rows must match tables/prediction-record.tex, the table the paper
     prints, through the same two formatting steps the table's generator uses
     (six significant figures for a measured value, the shortest spelling of
     the argument), and in the paper's own state and winner words.

Only then is the identical model applied to this sweep's cells on the held-out
grid (rule H, grid.py), and its commitments are scored against the measured
winners on that set.  The margin is the frozen 1.0 unit of the model's own
arithmetic term; nothing in this script can change it, because the point of
the exercise is that it was fixed before the new data existed.

The rule, as frozen:

    B = band maximum of the recursed bound (seed term + arithmetic term on the
        served orders; the direct error where the order is unserved, because
        both entries read the same fit there)
    D = band maximum of the direct error
    T = band maximum of the arithmetic term

    abstain   if |B - D| <= T
    direct    if D < B - T
    recursed  if D > B + T

Reported, under results/:

    prediction_heldout.csv  the 132 rows: B, D, T, state, measured winner,
                            the commitment and whether it was right
    per_order.csv           the measured band maxima and the winner per order
    band_summary.csv        the measured band worsts and the score per band
    summary.txt             the self-check, the score and every commitment

mpmath lives in the user site.  Do not run with `python -I`.
"""

import argparse
import csv
import json
import math
import sys
from pathlib import Path

if sys.flags.isolated:
    print("SKIP: isolated mode hides user-site mpmath; run without python -I")
    sys.exit(77)

from mpmath import mp, mpf  # noqa: E402

NMAX = 32
FROZEN_MARGIN = 1.0


def text(value):
    return mp.nstr(value, 17)


def csv_float(text):
    """float(), tolerating MSVC's spelling of a non-finite value."""
    return float(text.strip().replace("nan(ind)", "nan"))


def read_reference(path):
    reference = {}
    with open(path, newline="") as handle:
        reader = csv.reader(handle)
        header = next(reader)

        if header[:2] != ["band", "x"]:
            raise SystemExit("unexpected reference header: %r" % (header[:2],))

        for row in reader:
            reference[(row[0], float(row[1]))] = [mpf(value) for value in row[2:]]

    return reference


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
                    csv_float(row["recursed"]),
                    csv_float(row["direct"]),
                )
            )
    return rows


# ---- the frozen model, copied from ../amplification_form/analyze.py ---------


def amplification(n, x):
    if n == 0:
        return mpf(1)
    odd = math.factorial(2 * n) // (2**n * math.factorial(n))
    return mpf(odd) / (2 * x) ** n


def amplification_ratio(n, l, x):
    """The gain from step l's rounding to order n: A_n(x) / A_l(x)."""
    odd_n = math.factorial(2 * n) // (2**n * math.factorial(n))
    odd_l = math.factorial(2 * l) // (2**l * math.factorial(l))
    return mpf(odd_n) / mpf(odd_l) / (2 * x) ** (n - l)


def recurrence_arithmetic(n, x, values):
    """The paper's accumulated-arithmetic charge, absolute."""
    total = mpf(0)
    for l in range(1, n + 1):
        total += 3 * mpf(2) ** -53 * values[l] * amplification_ratio(n, l, x)

    if n > 0:
        carried = mpf(0)
        for l in range(n):
            carried += amplification_ratio(n, l, x) / (2 * l + 1)
        total += mpf(2) ** -52 * mp.exp(-x) * carried

    return total


class FrozenModel:
    """The paper's band-level model, on one cells file and one reference."""

    def __init__(self, cells, reference, thresholds):
        if len(thresholds) <= NMAX:
            raise SystemExit("the tree's threshold table stops before order %d" % NMAX)

        self.reference = reference
        self.thresholds = thresholds
        self.model_from = thresholds[0]

        ship_abs = {}
        ship_rel = {}
        prod_abs = {}
        prod_rel = {}
        error0 = {}
        winner = {}

        for band, x_text, x, n, recursed, direct in cells:
            values = reference[(band, x)]
            expected = values[n]
            rec_error = abs(mpf(recursed) - expected)
            dir_error = abs(mpf(direct) - expected)
            rec_relative = rec_error / expected if expected != 0 else mpf("inf")
            dir_relative = dir_error / expected if expected != 0 else mpf("inf")

            self._update(ship_abs, (band, n), rec_error, x_text)
            self._update(prod_abs, (band, n), dir_error, x_text)
            self._update(ship_rel, (band, n), rec_relative, x_text)
            self._update(prod_rel, (band, n), dir_relative, x_text)

            if n == 0:
                error0[(band, x)] = rec_error

        for band in ("A1", "A2"):
            for n in range(NMAX + 1):
                for metric, ship_store, prod_store in (
                    ("absolute", ship_abs, prod_abs),
                    ("relative", ship_rel, prod_rel),
                ):
                    ship_value = ship_store[(band, n)][0]
                    prod_value = prod_store[(band, n)][0]
                    if ship_value < prod_value:
                        state = "ship"
                    elif prod_value < ship_value:
                        state = "product"
                    else:
                        state = "tie"
                    winner[(band, metric, n)] = state

        self.ship_abs = ship_abs
        self.ship_rel = ship_rel
        self.prod_abs = prod_abs
        self.prod_rel = prod_rel
        self.winner = winner
        self.error0 = error0

        def zeros():
            return {"A1": [mpf(0) for _ in range(NMAX + 1)], "A2": [mpf(0) for _ in range(NMAX + 1)]}

        def nones():
            return {"A1": [None] * (NMAX + 1), "A2": [None] * (NMAX + 1)}

        predicted_rec_abs = zeros()
        predicted_rec_rel = zeros()
        predicted_dir_abs = zeros()
        predicted_dir_rel = zeros()
        predicted_arith_abs = zeros()
        predicted_arith_rel = zeros()
        predicted_rec_abs_at = nones()
        predicted_rec_rel_at = nones()
        predicted_dir_abs_at = nones()
        predicted_dir_rel_at = nones()
        model_arguments = {"A1": 0, "A2": 0}
        model_seen = set()

        for band, x_text, x, n, recursed, direct in cells:
            if x < self.model_from:
                continue
            if (band, x) not in model_seen:
                model_seen.add((band, x))
                model_arguments[band] += 1

            seed = error0.get((band, x))
            if seed is None:
                continue
            expected = reference[(band, x)][n]
            dir_abs = abs(mpf(direct) - expected)
            dir_rel = dir_abs / expected if expected != 0 else mpf("inf")

            if n == 0 or x >= thresholds[n]:
                arith_abs = recurrence_arithmetic(n, mpf(x), reference[(band, x)]) if n else mpf(0)
                rec_abs = seed * amplification(n, mpf(x)) + arith_abs
                rec_rel = rec_abs / expected if expected != 0 else mpf("inf")
                arith_rel = arith_abs / expected if expected != 0 else mpf("inf")
            else:
                arith_abs = mpf(0)
                arith_rel = mpf(0)
                rec_abs = dir_abs
                rec_rel = dir_rel

            if rec_abs > predicted_rec_abs[band][n]:
                predicted_rec_abs[band][n] = rec_abs
                predicted_rec_abs_at[band][n] = x_text
            if rec_rel > predicted_rec_rel[band][n]:
                predicted_rec_rel[band][n] = rec_rel
                predicted_rec_rel_at[band][n] = x_text
            if dir_abs > predicted_dir_abs[band][n]:
                predicted_dir_abs[band][n] = dir_abs
                predicted_dir_abs_at[band][n] = x_text
            if dir_rel > predicted_dir_rel[band][n]:
                predicted_dir_rel[band][n] = dir_rel
                predicted_dir_rel_at[band][n] = x_text
            if arith_abs > predicted_arith_abs[band][n]:
                predicted_arith_abs[band][n] = arith_abs
            if arith_rel > predicted_arith_rel[band][n]:
                predicted_arith_rel[band][n] = arith_rel

        self.predicted_rec_abs = predicted_rec_abs
        self.predicted_rec_rel = predicted_rec_rel
        self.predicted_dir_abs = predicted_dir_abs
        self.predicted_dir_rel = predicted_dir_rel
        self.predicted_arith_abs = predicted_arith_abs
        self.predicted_arith_rel = predicted_arith_rel
        self.predicted_rec_abs_at = predicted_rec_abs_at
        self.predicted_rec_rel_at = predicted_rec_rel_at
        self.predicted_dir_abs_at = predicted_dir_abs_at
        self.predicted_dir_rel_at = predicted_dir_rel_at
        self.model_arguments = model_arguments

    @staticmethod
    def _update(store, key, value, x_text):
        current = store.get(key)
        if current is None or value > current[0]:
            store[key] = (value, x_text)

    def stores(self, band, metric):
        if metric == "absolute":
            return (
                self.predicted_rec_abs[band],
                self.predicted_dir_abs[band],
                self.predicted_arith_abs[band],
                self.ship_abs,
                self.prod_abs,
            )
        return (
            self.predicted_rec_rel[band],
            self.predicted_dir_rel[band],
            self.predicted_arith_rel[band],
            self.ship_rel,
            self.prod_rel,
        )

    def prediction_rows(self, band, metric, margin=FROZEN_MARGIN):
        rows = []
        ship, product, sigma, _, _ = self.stores(band, metric)

        for n in range(NMAX + 1):
            ps = ship[n]
            pd = product[n]
            gap = abs(ps - pd)
            if gap <= mpf(margin) * sigma[n]:
                state = "abstain"
                choice = ""
            else:
                state = "ship" if ps < pd else "product"
                choice = state
            measured = self.winner[(band, metric, n)]
            right = choice == measured if state != "abstain" else False
            rows.append((n, ps, pd, sigma[n], state, measured, choice, right))

        return rows

    @staticmethod
    def score(rows):
        commits = sum(1 for row in rows if row[4] != "abstain")
        abstains = len(rows) - commits
        right_committed = sum(1 for row in rows if row[4] != "abstain" and row[7])
        return commits, right_committed, abstains


# ---- the self-check: the paper's own record ---------------------------------


def model_record_rows(model):
    rows = []
    for band in ("A1", "A2"):
        for metric in ("absolute", "relative"):
            if metric == "absolute":
                ship_at = model.predicted_rec_abs_at[band]
                product_at = model.predicted_dir_abs_at[band]
            else:
                ship_at = model.predicted_rec_rel_at[band]
                product_at = model.predicted_dir_rel_at[band]
            for n, ps, pd, sigma, state, measured, choice, right in model.prediction_rows(band, metric):
                rows.append(
                    {
                        "band": band,
                        "metric": metric,
                        "order": str(n),
                        "predicted_ship": text(ps),
                        "predicted_ship_at_x": ship_at[n],
                        "predicted_product": text(pd),
                        "predicted_product_at_x": product_at[n],
                        "arithmetic": text(sigma),
                        "state": state,
                        "measured_winner": measured,
                        "choice": choice,
                        "right_where_committed": "1" if right else "0",
                    }
                )
    return rows


def read_prediction_record(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle))


def read_tex_rows(path):
    rows = []
    for line in Path(path).read_text(encoding="ascii").splitlines():
        line = line.strip()
        if not line or line.startswith("%") or "&" not in line:
            continue
        rows.append([field.strip() for field in line.rstrip("\\").split("&")])
    return rows


def self_check(second_sweep, tex_path, thresholds):
    reference = read_reference(second_sweep / "data/reference.csv")
    cells = read_cells(second_sweep / "results/cells.csv")
    model = FrozenModel(cells, reference, thresholds)
    mine = model_record_rows(model)

    recorded = read_prediction_record(second_sweep / "results/prediction_record.csv")
    if len(recorded) != len(mine):
        raise SystemExit("self-check: %d recorded rows, model produced %d" % (len(recorded), len(mine)))

    fields = [
        "band",
        "metric",
        "order",
        "predicted_ship",
        "predicted_ship_at_x",
        "predicted_product",
        "predicted_product_at_x",
        "arithmetic",
        "state",
        "measured_winner",
        "choice",
        "right_where_committed",
    ]
    for index, (got, want) in enumerate(zip(mine, recorded)):
        for field in fields:
            if got[field] != want[field]:
                raise SystemExit(
                    "self-check: row %d field %s: model %r, record %r"
                    % (index + 1, field, got[field], want[field])
                )

    tex = read_tex_rows(tex_path)
    if len(tex) != len(mine):
        raise SystemExit("self-check: %d tex rows, model produced %d" % (len(tex), len(mine)))

    band_word = {"A1": "B1", "A2": "B2"}
    entry_word = {"ship": "recursed", "product": "direct", "tie": "tie", "abstain": "abstain"}
    for index, (row, fields_tex) in enumerate(zip(mine, tex)):
        if len(fields_tex) != 11:
            raise SystemExit("self-check: tex row %d has %d fields: %r" % (index + 1, len(fields_tex), fields_tex))
        # The two format steps of make_prediction_record.py, the table's own
        # generator: a measured value at six significant figures, an argument
        # as the shortest spelling of its double.
        expected = [
            band_word[row["band"]],
            row["metric"],
            row["order"],
            "%.6g" % float(row["predicted_ship"]),
            repr(float(row["predicted_ship_at_x"])),
            "%.6g" % float(row["predicted_product"]),
            repr(float(row["predicted_product_at_x"])),
            "%.6g" % float(row["arithmetic"]),
            entry_word[row["state"]],
            entry_word[row["measured_winner"]],
            "---" if row["state"] == "abstain" else ("yes" if row["right_where_committed"] == "1" else "no"),
        ]
        for field, got, want in zip(fields, fields_tex, expected):
            if got != want:
                raise SystemExit(
                    "self-check: tex row %d (%s %s n=%s): %r is %r, model says %r"
                    % (index + 1, row["band"], row["metric"], row["order"], field, got, want)
                )

    rows = []
    for band in ("A1", "A2"):
        for metric in ("absolute", "relative"):
            rows += model.prediction_rows(band, metric)
    commits, right, abstains = FrozenModel.score(rows)

    return (
        "self-check passed: %d rows reproduced field for field from "
        "amplification_form/results/prediction_record.csv and from "
        "tables/prediction-record.tex; on the second sweep's own cells the "
        "frozen rule commits %d, is right %d, abstains %d"
        % (len(mine), commits, right, abstains)
    )


# ---- the held-out scoring ---------------------------------------------------


def held_out(tree, thresholds, second_sweep, tex_path):
    check_line = self_check(second_sweep, tex_path, thresholds)

    reference = read_reference(tree / "data/reference.csv")
    cells = read_cells(tree / "results/cells.csv")
    model = FrozenModel(cells, reference, thresholds)

    results = tree / "results"
    results.mkdir(parents=True, exist_ok=True)

    entry_word = {"ship": "recursed", "product": "direct", "tie": "tie", "abstain": "abstain"}

    with open(results / "prediction_heldout.csv", "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "band",
                "metric",
                "order",
                "predicted_recursed",
                "predicted_recursed_at_x",
                "predicted_direct",
                "predicted_direct_at_x",
                "arithmetic",
                "state",
                "measured_winner",
                "choice",
                "right_where_committed",
            ]
        )
        for band in ("A1", "A2"):
            for metric in ("absolute", "relative"):
                if metric == "absolute":
                    rec_at = model.predicted_rec_abs_at[band]
                    dir_at = model.predicted_dir_abs_at[band]
                else:
                    rec_at = model.predicted_rec_rel_at[band]
                    dir_at = model.predicted_dir_rel_at[band]
                for n, ps, pd, sigma, state, measured, choice, right in model.prediction_rows(band, metric):
                    writer.writerow(
                        [
                            band,
                            metric,
                            n,
                            text(ps),
                            rec_at[n],
                            text(pd),
                            dir_at[n],
                            text(sigma),
                            entry_word[state],
                            entry_word[measured],
                            entry_word[choice] if choice else "",
                            "1" if right else "0",
                        ]
                    )

    with open(results / "per_order.csv", "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "band",
                "metric",
                "order",
                "measured_recursed_worst",
                "measured_recursed_at_x",
                "measured_direct_worst",
                "measured_direct_at_x",
                "measured_winner",
                "state",
                "choice",
                "right_where_committed",
            ]
        )
        for band in ("A1", "A2"):
            for metric in ("absolute", "relative"):
                _, _, _, ship_store, prod_store = model.stores(band, metric)
                for n, ps, pd, sigma, state, measured, choice, right in model.prediction_rows(band, metric):
                    writer.writerow(
                        [
                            band,
                            metric,
                            n,
                            text(ship_store[(band, n)][0]),
                            ship_store[(band, n)][1],
                            text(prod_store[(band, n)][0]),
                            prod_store[(band, n)][1],
                            entry_word[measured],
                            entry_word[state],
                            entry_word[choice] if choice else "",
                            "1" if right else "0",
                        ]
                    )

    with open(results / "band_summary.csv", "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "band",
                "metric",
                "measured_recursed_worst",
                "measured_recursed_at_n",
                "measured_recursed_at_x",
                "measured_direct_worst",
                "measured_direct_at_n",
                "measured_direct_at_x",
                "ratio_of_band_worsts",
                "commits",
                "right_where_committed",
                "abstains",
            ]
        )
        for band in ("A1", "A2"):
            for metric in ("absolute", "relative"):
                _, _, _, ship_store, prod_store = model.stores(band, metric)
                ship_worst, ship_n, ship_x = max(
                    ((ship_store[(band, n)][0], n, ship_store[(band, n)][1]) for n in range(NMAX + 1)),
                    key=lambda entry: entry[0],
                )
                prod_worst, prod_n, prod_x = max(
                    ((prod_store[(band, n)][0], n, prod_store[(band, n)][1]) for n in range(NMAX + 1)),
                    key=lambda entry: entry[0],
                )
                commits, right, abstains = FrozenModel.score(model.prediction_rows(band, metric))
                writer.writerow(
                    [
                        band,
                        metric,
                        text(ship_worst),
                        ship_n,
                        ship_x,
                        text(prod_worst),
                        prod_n,
                        prod_x,
                        text(ship_worst / prod_worst if prod_worst > 0 else mpf("inf")),
                        commits,
                        right,
                        abstains,
                    ]
                )

    totals = [0, 0, 0]
    band_lines = []
    commitments = []
    for band in ("A1", "A2"):
        for metric in ("absolute", "relative"):
            rows = model.prediction_rows(band, metric)
            commits, right, abstains = FrozenModel.score(rows)
            totals[0] += commits
            totals[1] += right
            totals[2] += abstains
            band_lines.append(
                "  %s %-8s commits %2d, right %2d, wrong %2d, abstains %2d of 33"
                % (band, metric, commits, right, commits - right, abstains)
            )
            for n, ps, pd, sigma, state, measured, choice, right_flag in rows:
                if state == "abstain":
                    continue
                commitments.append(
                    "    %s %s n=%2d: %s (B %s, D %s, T %s; measured %s)%s"
                    % (
                        band,
                        metric,
                        n,
                        entry_word[choice],
                        text(ps),
                        text(pd),
                        text(sigma),
                        entry_word[measured],
                        "" if right_flag else "   WRONG",
                    )
                )

    validation_lines = []
    validation_path = tree / "data/reference_validation.txt"
    if validation_path.is_file():
        validation_lines = validation_path.read_text(encoding="utf-8").splitlines()

    admission_lines = []
    admission_totals = {}
    admission_mismatches = []
    admission_path = tree / "results/admission.csv"
    if admission_path.is_file():
        with open(admission_path, newline="") as handle:
            for row in csv.DictReader(handle):
                entry = admission_totals.setdefault(row["direction"], [0, 0, 0])
                entry[0] += 1
                entry[1] += int(row["checked"])
                entry[2] += int(row["matched"])
                if row["checked"] != row["matched"]:
                    admission_mismatches.append((row["band"], row["x"], row["direction"]))
    for direction in sorted(admission_totals):
        entries, checked, matched = admission_totals[direction]
        admission_lines.append(
            "  %-9s arguments %4d, cells %6d, bit-identical %6d%s"
            % (direction, entries, checked, matched, "" if checked == matched else "  MISMATCH")
        )

    headcount = {}
    with open(tree / "data/grid.csv", newline="") as handle:
        for row in csv.DictReader(handle):
            headcount[row["band"]] = headcount.get(row["band"], 0) + 1

    commits, right, abstains = totals
    summary = [
        "held-out summary",
        "reference: %s" % (tree / "data/reference.csv").resolve(),
        "grid: %d arguments (A1 %d, A2 %d), rule H in grid.py"
        % (sum(headcount.values()), headcount.get("A1", 0), headcount.get("A2", 0)),
        "model population: arguments at or above kTierThresholds[0] = %r: A1 %d, A2 %d"
        % (model.model_from, model.model_arguments["A1"], model.model_arguments["A2"]),
        "",
        check_line,
        "The margin is FROZEN_MARGIN = 1.0 unit of T, as in the paper's record;",
        "this script has no option that changes it, and the grid was fixed before",
        "this run.",
        "",
        "admission:",
    ]
    summary += admission_lines
    for mismatch in admission_mismatches:
        summary.append("  mismatch: %s" % (mismatch,))
    summary += [
        "",
        "the frozen rule: abstain iff |B - D| <= T; direct iff D < B - T;",
        "recursed iff D > B + T.  B = band max of the recursed bound (seed term",
        "plus arithmetic where the order is served, the direct error where it is",
        "not); D = band max of the direct error; T = band max of the arithmetic",
        "term.  The winner is measured on this grid.",
        "",
        "held-out score, orders 0-32:",
    ]
    summary += band_lines
    summary += [
        "  all                commits %2d, right %2d, wrong %2d, abstains %2d of 132"
        % (commits, right, commits - right, abstains),
        "  abstention rate: %d/132 = %.1f per cent" % (abstains, 100.0 * abstains / 132),
        "",
        "every commitment:",
    ]
    summary += commitments
    summary += ["", "reference validation:"]
    summary += ["  " + line for line in validation_lines]

    with open(results / "summary.txt", "w", encoding="utf-8") as handle:
        handle.write("\n".join(summary) + "\n")

    print("\n".join(summary))
    return 0


def main():
    here = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tree", default=str(here))
    parser.add_argument("--provenance", default=None)
    parser.add_argument("--second-sweep", default=str(here.parent / "amplification_form"))
    parser.add_argument("--tex", default=str(here.parents[1] / "tables/prediction-record.tex"))
    args = parser.parse_args()

    tree = Path(args.tree).resolve()
    second_sweep = Path(args.second_sweep).resolve()
    tex_path = Path(args.tex).resolve()

    provenance_path = Path(args.provenance) if args.provenance else tree / "results/provenance.json"
    if not provenance_path.is_file():
        raise SystemExit("provenance.json is missing; it carries kTierThresholds")
    with open(provenance_path, encoding="utf-8") as handle:
        constants = json.load(handle)["constants"]

    thresholds = [float(value) for value in constants["kTierThresholds"]]

    mp.dps = 70
    return held_out(tree, thresholds, second_sweep, tex_path)


if __name__ == "__main__":
    sys.exit(main())
