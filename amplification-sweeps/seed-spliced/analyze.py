#!/usr/bin/env python3
"""Reduce the seed-spliced sweep: does the seed term order the comparison?

sweep_seed.cpp writes one row per argument and order: the shipped recursed
entry, the direct matrix-product entry, and a spliced ladder whose order-0
value is the direct entry's own and whose higher orders come from the
kernel's own upward step.  The spliced ladder is therefore the same initial
question the paper asks - what does the recurrence deliver from a seed error
of this size? - with the seed error taken from the direct entry's own fit
family instead of from a different stored table.

For every order this script forms, on the model's population (arguments at
or above kTierThresholds[0], the same population the second sweep's model
uses):

    eps0(x)  = |direct_0(x) - F_0_ref(x)|          the seed error, measured
    S_n(x)   = eps0(x) A_n(x)                      the seed term
    E_n(x)   = |spliced_n(x) - F_n_ref(x)|         the spliced delivered error
    D_n(x)   = |direct_n(x)  - F_n_ref(x)|         the direct delivered error

and takes the band maximum of each, exactly as the paper's model compares
band maxima.  The seed-term comparison is

    predicted spliced  if S < D,   predicted direct if S > D,

and it is scored against the measured winner (spliced if E < D, direct if
D < E).  Relative errors divide each by the reference value at the order.

Order 0 is structural: A_0 = 1, so S_0 = eps0 = D_0 = E_0 at every argument
and the comparison is a tie between the same measured quantity.  The four
order-0 band-metric comparisons are reported and excluded, leaving 128.

The full 836-argument grid is also attempted and reported as not measurable:
x = 0 is in the grid, A_n(0) is not finite for n >= 1, and as x -> 0 both the
seed term and the spliced ladder overflow, so every band maximum is an
overflow artefact rather than a measured quantity.  The population is the
paper's own domain and the only place the comparison is defined.

Reported, under results/:

    per_order.csv        one row per band, metric and order: the three band
                         maxima with their arguments, the prediction, the
                         measured winner, and whether the prediction was right
    seed_term_cells.csv  one row per population argument and order: the seed
                         term, the spliced error and the direct error, absolute
                         and relative
    summary.txt          counts, admission, replica agreement, the score

mpmath lives in the user site.  Do not run with `python -I`.
"""

import argparse
import csv
import json
import math
import os
import sys
from pathlib import Path

if sys.flags.isolated:
    print("SKIP: isolated mode hides user-site mpmath; run without python -I")
    sys.exit(77)

from mpmath import mp, mpf  # noqa: E402

NMAX = 32


def text(value):
    if value is None:
        return ""
    if isinstance(value, float):
        return "%.17g" % value
    if not mp.isfinite(value):
        return "inf" if value > 0 else ("-inf" if value < 0 else "nan")
    return mp.nstr(value, 17)


def cell_float(text_value):
    """The double a cells.csv field spells; non-finite if it does not spell one."""
    try:
        return float(text_value)
    except ValueError:
        return float("nan")


def amplification(n, x):
    if n == 0:
        return mpf(1)
    odd = math.factorial(2 * n) // (2**n * math.factorial(n))
    return mpf(odd) / (2 * x) ** n


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
                    cell_float(row["shipped"]),
                    cell_float(row["direct"]),
                    cell_float(row["spliced"]),
                )
            )
    return rows


def read_grid(path):
    counts = {"A1": 0, "A2": 0}
    zero = 0
    with open(path, newline="") as handle:
        for row in csv.DictReader(handle):
            counts[row["band"]] = counts.get(row["band"], 0) + 1
            if float(row["x"]) == 0.0:
                zero += 1
    return counts, zero


def read_admission(path):
    totals = {}
    mismatches = []
    with open(path, newline="") as handle:
        for row in csv.DictReader(handle):
            entry = totals.setdefault(row["direction"], [0, 0, 0])
            entry[0] += 1
            entry[1] += int(row["checked"])
            entry[2] += int(row["matched"])
            if row["checked"] != row["matched"]:
                mismatches.append((row["band"], row["x"], row["direction"]))
    return totals, mismatches


def read_steps(path):
    arguments = 0
    steps = 0
    matched = 0
    bad = []
    with open(path, newline="") as handle:
        for row in csv.DictReader(handle):
            arguments += 1
            steps += int(row["steps"])
            matched += int(row["matched"])
            if row["steps"] != row["matched"]:
                bad.append((row["band"], row["x"], row["steps"], row["matched"]))
    return arguments, steps, matched, bad


def error_of(value, expected):
    """|delivered - reference|; a non-finite delivered value is an infinite error."""
    if not math.isfinite(value):
        return mpf("inf")
    return abs(mpf(value) - expected)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tree", default=str(Path(__file__).resolve().parent))
    parser.add_argument("--provenance", default=None)
    args = parser.parse_args()

    tree = Path(args.tree).resolve()
    provenance_path = Path(args.provenance) if args.provenance else tree / "results/provenance.json"
    if not provenance_path.is_file():
        raise SystemExit("provenance.json is missing; it carries kTierThresholds[0]")
    with open(provenance_path, encoding="utf-8") as handle:
        constants = json.load(handle)["constants"]

    thresholds = [float(value) for value in constants["kTierThresholds"]]
    if len(thresholds) <= NMAX:
        raise SystemExit("the threshold table stops before order %d" % NMAX)
    model_from = thresholds[0]

    mp.dps = 70
    reference = read_reference(tree / "data/reference.csv")
    cells = read_cells(tree / "results/cells.csv")
    admission_totals, admission_mismatches = read_admission(tree / "results/admission.csv")
    step_arguments, step_total, step_matched, step_bad = read_steps(tree / "results/replica_steps.csv")
    grid_counts, zero_arguments = read_grid(tree / "data/grid.csv")
    grid_total = grid_counts["A1"] + grid_counts["A2"]

    # ---- per-cell errors on the model population ---------------------------
    #
    # The band maxima are the objects the paper compares; their arguments can
    # differ between the seed term, the spliced error and the direct error.
    # Every store is indexed by band: an A1 maximum and an A2 maximum are
    # different objects even when they share an order.
    BANDS = ("A1", "A2")

    def empties():
        return ({band: [mpf(0) for _ in range(NMAX + 1)] for band in BANDS},
                {band: [None] * (NMAX + 1) for band in BANDS})

    seed_term_abs, seed_term_abs_at = empties()
    seed_term_rel, seed_term_rel_at = empties()
    spliced_abs, spliced_abs_at = empties()
    spliced_rel, spliced_rel_at = empties()
    direct_abs, direct_abs_at = empties()
    direct_rel, direct_rel_at = empties()
    population = {"A1": set(), "A2": set()}
    outside = {"A1": set(), "A2": set()}
    outside_nonfinite = {"A1": 0, "A2": 0}

    # The direct entry's order-0 value per argument: the spliced ladder's seed.
    direct0 = {}
    for band, x_text, x, n, shipped, direct, spliced in cells:
        if n == 0:
            direct0[(band, x_text)] = direct

    cell_rows = []
    seed_ratio = []  # E / S over population cells with both finite

    for band, x_text, x, n, shipped, direct, spliced in cells:
        if x < model_from:
            outside[band].add((x_text, x))
            if not math.isfinite(spliced):
                outside_nonfinite[band] += 1
            continue
        population[band].add((x_text, x))

        values = reference[(band, x)]
        expected = values[n]
        direct_error = error_of(direct, expected)
        spliced_error = error_of(spliced, expected)
        if n == 0:
            # The seed error itself; A_0 = 1.
            seed_value = direct_error
        else:
            seed_error = error_of(direct0[(band, x_text)], reference[(band, x)][0])
            seed_value = seed_error * amplification(n, mpf(x))

        direct_relative = direct_error / expected if expected != 0 else mpf("inf")
        spliced_relative = spliced_error / expected if expected != 0 else mpf("inf")
        seed_relative = seed_value / expected if expected != 0 else mpf("inf")

        for store, at, value in (
            (seed_term_abs, seed_term_abs_at, seed_value),
            (seed_term_rel, seed_term_rel_at, seed_relative),
            (spliced_abs, spliced_abs_at, spliced_error),
            (spliced_rel, spliced_rel_at, spliced_relative),
            (direct_abs, direct_abs_at, direct_error),
            (direct_rel, direct_rel_at, direct_relative),
        ):
            bucket = store[band]
            where = at[band]
            if value > bucket[n] or (where[n] is None and value == bucket[n]):
                bucket[n] = value
                where[n] = x_text

        if n >= 1 and mp.isfinite(seed_value) and seed_value > 0 and mp.isfinite(spliced_error):
            seed_ratio.append(spliced_error / seed_value)

        cell_rows.append(
            (
                band,
                x_text,
                n,
                text(seed_value),
                text(seed_relative),
                text(direct_error),
                text(direct_relative),
                text(spliced_error),
                text(spliced_relative),
            )
        )

    # ---- the band-metric comparison ---------------------------------------
    row_records = []
    for band in BANDS:
        for metric in ("absolute", "relative"):
            if metric == "absolute":
                seed_store, seed_at = seed_term_abs, seed_term_abs_at
                spliced_store, spliced_at = spliced_abs, spliced_abs_at
                direct_store, direct_at = direct_abs, direct_abs_at
            else:
                seed_store, seed_at = seed_term_rel, seed_term_rel_at
                spliced_store, spliced_at = spliced_rel, spliced_rel_at
                direct_store, direct_at = direct_rel, direct_rel_at

            for n in range(NMAX + 1):
                seed_value = seed_store[band][n]
                spliced_value = spliced_store[band][n]
                direct_value = direct_store[band][n]
                if n == 0:
                    state = "structural"
                    predicted = ""
                    measured = "tie"
                    correct = ""
                else:
                    measured = (
                        "spliced"
                        if spliced_value < direct_value
                        else ("direct" if direct_value < spliced_value else "tie")
                    )
                    if seed_value == direct_value:
                        state = "tie"
                        predicted = ""
                        correct = ""
                    else:
                        state = "ordered"
                        predicted = "direct" if seed_value > direct_value else "spliced"
                        correct = "1" if predicted == measured else "0"
                row_records.append(
                    (
                        band,
                        metric,
                        n,
                        text(seed_value),
                        seed_at[band][n] or "",
                        text(direct_value),
                        direct_at[band][n] or "",
                        text(spliced_value),
                        spliced_at[band][n] or "",
                        state,
                        predicted,
                        measured,
                        correct,
                    )
                )

    def tally(records):
        orderable = [row for row in records if row[9] == "ordered"]
        ties = [row for row in records if row[9] == "tie"]
        structural = [row for row in records if row[9] == "structural"]
        right = sum(1 for row in orderable if row[12] == "1")
        return len(records), len(orderable), len(ties), len(structural), right

    def line_for(label, records):
        total, orderable, ties, structural, right = tally(records)
        return (
            "  %-18s comparisons %3d: structural ties %d, seed ties %d, "
            "ordered %3d, ordered correctly %3d"
            % (label, total, structural, ties, orderable, right)
        )

    # ---- the full-grid attempt --------------------------------------------
    spliced_nonfinite = {"A1": 0, "A2": 0}
    for band, x_text, x, n, shipped, direct, spliced in cells:
        if not math.isfinite(spliced):
            spliced_nonfinite[band] += 1

    outside_counts = {band: len(outside[band]) for band in ("A1", "A2")}
    population_counts = {band: len(population[band]) for band in ("A1", "A2")}

    ratio_sorted = sorted(seed_ratio)
    def quartile(fraction):
        if not ratio_sorted:
            return None
        index = min(len(ratio_sorted) - 1, max(0, int(fraction * (len(ratio_sorted) - 1))))
        return ratio_sorted[index]

    # Cell-level sign agreement, for the robustness line.
    cell_ordered = 0
    cell_right = 0
    for band, x_text, x, n, shipped, direct, spliced in cells:
        if x < model_from or n == 0:
            continue
        values = reference[(band, x)]
        direct_error = error_of(direct, values[n])
        spliced_error = error_of(spliced, values[n])
        seed_error = error_of(direct0[(band, x_text)], values[0])
        seed_value = seed_error * amplification(n, mpf(x))
        if not (mp.isfinite(seed_value) and mp.isfinite(spliced_error) and mp.isfinite(direct_error)):
            continue
        if seed_value == direct_error:
            continue  # the same fit at an unserved order would be a tie; none expected here
        cell_ordered += 1
        predicted = "direct" if seed_value > direct_error else "spliced"
        measured = "spliced" if spliced_error < direct_error else ("direct" if direct_error < spliced_error else "tie")
        cell_right += 1 if predicted == measured else 0

    # ---- write ------------------------------------------------------------
    results = tree / "results"
    results.mkdir(parents=True, exist_ok=True)

    with open(results / "per_order.csv", "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "band",
                "metric",
                "order",
                "seed_term",
                "seed_term_at_x",
                "direct_worst",
                "direct_at_x",
                "spliced_worst",
                "spliced_at_x",
                "state",
                "predicted",
                "measured_winner",
                "right_where_ordered",
            ]
        )
        writer.writerows(row_records)

    with open(results / "seed_term_cells.csv", "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "band",
                "x",
                "n",
                "seed_term",
                "seed_term_relative",
                "direct_error",
                "direct_error_relative",
                "spliced_error",
                "spliced_error_relative",
                "in_population",
            ]
        )
        writer.writerows(cell_rows)

    validation_lines = []
    validation_path = tree / "data/reference_validation.txt"
    if validation_path.is_file():
        validation_lines = validation_path.read_text(encoding="utf-8").splitlines()

    summary = [
        "seed-spliced summary",
        "reference: %s" % (tree / "data/reference.csv").resolve(),
        "grid: %d arguments (A1 %d, A2 %d); x = 0 arguments: %d"
        % (grid_total, grid_counts["A1"], grid_counts["A2"], zero_arguments),
        "",
        "population: arguments at or above kTierThresholds[0] = %r"
        % model_from,
        "  A1 %d of %d, A2 %d of %d; outside the population %d (A1 %d, A2 %d)"
        % (
            population_counts["A1"],
            grid_counts["A1"],
            population_counts["A2"],
            grid_counts["A2"],
            outside_counts["A1"] + outside_counts["A2"],
            outside_counts["A1"],
            outside_counts["A2"],
        ),
        "  outside the population the entry below kTierThresholds[0] seeds at order 32 and",
        "  recurses downward; its order-0 value is not a seed error and A_n does not",
        "  describe that path.  The spliced ladder is defined there, but its band maxima",
        "  are dominated by x -> 0: A_n(x) diverges as x^-n, and the ladder overflows.",
        "  Non-finite spliced cells outside the population: A1 %d, A2 %d."
        % (outside_nonfinite["A1"], outside_nonfinite["A2"]),
        "  The full 836-argument comparison is therefore not measurable: with x = 0 in the",
        "  grid, A_n(0) is not finite for n >= 1, and every band maximum is an overflow",
        "  artefact.  The score below is on the paper's own population.",
        "",
        "admission:",
    ]
    for direction in sorted(admission_totals):
        entries, checked, matched = admission_totals[direction]
        summary.append(
            "  %-9s arguments %4d, cells %6d, bit-identical %6d%s"
            % (direction, entries, checked, matched, "" if checked == matched else "  MISMATCH")
        )
    for mismatch in admission_mismatches:
        summary.append("  mismatch: %s" % (mismatch,))
    summary += [
        "",
        "replica step agreement: %d steps over %d arguments, bit-identical %d, mismatched rows %d"
        % (step_total, step_arguments, step_matched, len(step_bad)),
    ]
    for bad in step_bad:
        summary.append("  mismatch: %s" % (bad,))
    summary += [
        "  the spliced ladder is the kernel's own UpwardStep body; the replica row shows",
        "  the two bodies agree everywhere the step was taken.",
        "",
        "order 0 is structural: A_0 = 1, so the seed term is the direct entry's own order-0",
        "error, the spliced ladder's order 0 is that same value, and the comparison is a tie",
        "between the same measured quantity.  The four order-0 band-metric rows are recorded",
        "and excluded from the score.",
        "",
        "seed-term comparison, orders 1-32, on the population:",
    ]
    summary.append(line_for("all", [row for row in row_records if row[2] >= 1]))
    for band in ("A1", "A2"):
        for metric in ("absolute", "relative"):
            rows = [row for row in row_records if row[0] == band and row[1] == metric and row[2] >= 1]
            summary.append(line_for("%s %s" % (band, metric), rows))
    summary.append("  predicted spliced if seed term < direct worst; predicted direct otherwise;")
    summary.append("  right where the measured winner (spliced if its error < direct's) agrees.")
    summary.append("")

    wrong = [row for row in row_records if row[2] >= 1 and row[9] == "ordered" and row[12] == "0"]
    summary.append("rows the seed term orders incorrectly: %d" % len(wrong))
    for row in wrong:
        summary.append(
            "  %s %-8s n=%2d seed %s at %s; direct %s at %s; spliced %s at %s; predicted %s"
            % (row[0], row[1], row[2], row[3], row[4], row[5], row[6], row[7], row[8], row[10])
        )
    ties = [row for row in row_records if row[9] == "tie"]
    if ties:
        summary.append("seed ties at orders above 0: %d" % len(ties))
        for row in ties:
            summary.append("  %s %s n=%d seed %s direct %s" % (row[0], row[1], row[2], row[3], row[5]))
    summary += [
        "",
        "spliced error over seed term, population cells with both finite (n >= 1): %d cells"
        % len(ratio_sorted),
    ]
    if ratio_sorted:
        summary.append(
            "  quartiles %s, %s, %s; maximum %s"
            % (
                text(quartile(0.25)),
                text(quartile(0.5)),
                text(quartile(0.75)),
                text(ratio_sorted[-1]),
            )
        )
        summary.append("  a ratio below one means the delivered spliced error beat the seed term's magnitude.")
    summary += [
        "",
        "cell-level sign agreement (population, n >= 1, all three finite):",
        "  %d cells ordered, %d right where ordered" % (cell_ordered, cell_right),
        "",
        "reference validation:",
    ]
    summary += ["  " + line for line in validation_lines]

    with open(results / "summary.txt", "w", encoding="utf-8") as handle:
        handle.write("\n".join(summary) + "\n")

    print("\n".join(summary))
    return 0


if __name__ == "__main__":
    sys.exit(main())
