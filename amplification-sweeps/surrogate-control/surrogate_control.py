#!/usr/bin/env python3
"""Referee control: can a margin-retuned monotone surrogate match the true
model's 67 commitments with zero errors?

Reads the existing amplification_form records.  Writes nothing.

The prediction block is a verbatim replication of analyze.py's main(): the
model population is x >= kTierThresholds[0], the recursed bound is
seed * A_n(x) + accumulated arithmetic, the direct bound is the measured
direct error, the abstention threshold is the recursed bound's own
accumulated-arithmetic term.  Only A_n(x) is replaced, for the surrogate
runs; every other input (seed, direct bound, arithmetic term, measured
winner, population) is shared.

mpmath lives in the user site.  Do not run with `python -I`.
"""

import csv
import json
import math
import sys
from pathlib import Path

if sys.flags.isolated:
    print("SKIP: isolated mode hides user-site mpmath; run without python -I")
    sys.exit(77)

from mpmath import mp, mpf

TREE = Path(r"C:\Users\myaml\dev\papers-a-v5\papers\boys-crossover\sweeps\amplification_form")
sys.path.insert(0, str(TREE))
import analyze as A  # noqa: E402

mp.dps = 70
NMAX = A.NMAX


def text(value):
    return mp.nstr(value, 17)


def build_measured():
    """The measured stores and the band-metric measured winners (analyze.py)."""
    reference = A.read_reference(TREE / "data/reference.csv")
    cells = A.read_cells(TREE / "results/cells.csv")

    ship_abs, ship_rel, prod_abs, prod_rel = {}, {}, {}, {}
    error0 = {}

    for band, x_text, x, n, recursed, direct in cells:
        values = reference[(band, x)]
        expected = values[n]
        rec_error = abs(mpf(recursed) - expected)
        dir_error = abs(mpf(direct) - expected)
        rec_relative = rec_error / expected if expected != 0 else mpf("inf")
        dir_relative = dir_error / expected if expected != 0 else mpf("inf")

        A.update(ship_abs, (band, n), rec_error, x_text)
        A.update(prod_abs, (band, n), dir_error, x_text)
        A.update(ship_rel, (band, n), rec_relative, x_text)
        A.update(prod_rel, (band, n), dir_relative, x_text)

        if n == 0:
            error0[(band, x)] = rec_error

    winner = {}
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

    return reference, cells, error0, winner, ship_abs, ship_rel, prod_abs, prod_rel


def build_cells(reference, cells, error0):
    """One entry per population cell, with the branch decision of analyze.py.

    For x >= thresholds[n] (or n == 0) the recursed bound is
    seed * A_n(x) + arithmetic; above the served prefix both entries read the
    same fit and the recursed bound is the direct error.  The arithmetic term
    is computed once with the true recurrence gains and reused by every
    amplification variant.
    """
    constants = json.loads((TREE / "results/provenance.json").read_text(encoding="utf-8"))[
        "constants"
    ]
    thresholds = [float(value) for value in constants["kTierThresholds"]]
    if len(thresholds) <= NMAX:
        raise SystemExit("threshold table stops before order %d" % NMAX)
    model_from = thresholds[0]

    entries = {band: [] for band in ("A1", "A2")}
    seen = set()
    for band, x_text, x, n, recursed, direct in cells:
        if x < model_from:
            continue
        seed = error0.get((band, x))
        if seed is None:
            continue
        if (band, x) not in seen:
            seen.add((band, x))

        expected = reference[(band, x)][n]
        dir_abs = abs(mpf(direct) - expected)
        dir_rel = dir_abs / expected if expected != 0 else mpf("inf")

        if n == 0 or x >= thresholds[n]:
            arith_abs = A.recurrence_arithmetic(n, mpf(x), reference[(band, x)]) if n else mpf(0)
            arith_rel = arith_abs / expected if expected != 0 else mpf("inf")
            entries[band].append(
                ("rec", n, mpf(x), seed, expected, dir_abs, dir_rel, arith_abs, arith_rel)
            )
        else:
            entries[band].append(
                ("fit", n, mpf(x), seed, expected, dir_abs, dir_rel, mpf(0), mpf(0))
            )

    # Count distinct arguments as analyze.py does.
    args = {"A1": set(), "A2": set()}
    for band, x_text, x, n, recursed, direct in cells:
        if x >= model_from:
            args[band].add(x)
    population = {band: len(args[band]) for band in args}
    return entries, population, model_from


def predicted(entries, Afun, which):
    """The band-metric predicted bounds for one amplification rule.

    `which` selects (ship_store, product_store, arith_store): absolute or
    relative.  The arithmetic term is not re-derived from the surrogate: the
    paper replaces A_n only and keeps every other input unchanged.
    """
    rec = {"A1": [mpf(0)] * (NMAX + 1), "A2": [mpf(0)] * (NMAX + 1)}
    prod = {"A1": [mpf(0)] * (NMAX + 1), "A2": [mpf(0)] * (NMAX + 1)}
    arith = {"A1": [mpf(0)] * (NMAX + 1), "A2": [mpf(0)] * (NMAX + 1)}
    for band in ("A1", "A2"):
        for kind, n, x, seed, expected, dir_abs, dir_rel, arith_abs, arith_rel in entries[band]:
            if kind == "fit":
                ps, pd, sg = dir_abs, dir_abs, mpf(0)
                if which == "relative":
                    ps = pd = dir_rel
            else:
                sg = arith_abs if which == "absolute" else arith_rel
                if which == "absolute":
                    ps = seed * Afun(n, x) + arith_abs
                    pd = dir_abs
                else:
                    ps_mass = seed * Afun(n, x) + arith_abs
                    ps = ps_mass / expected if expected != 0 else mpf("inf")
                    pd = dir_rel
            if ps > rec[band][n]:
                rec[band][n] = ps
            if pd > prod[band][n]:
                prod[band][n] = pd
            if sg > arith[band][n]:
                arith[band][n] = sg
    return rec, prod, arith


def rows_of(rec, prod, arith, winner, band, metric, margin):
    rows = []
    for n in range(NMAX + 1):
        ps, pd, sigma = rec[band][n], prod[band][n], arith[band][n]
        gap = abs(ps - pd)
        if sigma == 0:
            abstain = gap <= 0
        else:
            abstain = gap <= mpf(margin) * sigma
        if abstain:
            state, choice, right = "abstain", "", False
        else:
            state = "ship" if ps < pd else "product"
            choice = state
            right = choice == winner[(band, metric, n)]
        rows.append(
            dict(
                n=n,
                ps=ps,
                pd=pd,
                sigma=sigma,
                gap=gap,
                measured=winner[(band, metric, n)],
                state=state,
                choice=choice,
                right=right,
            )
        )
    return rows


def score(rows):
    commits = sum(1 for r in rows if r["state"] != "abstain")
    right = sum(1 for r in rows if r["state"] != "abstain" and r["right"])
    return commits, right, commits - right


def winners_over(reference, cells, model_from=None):
    """Measured band-metric winners, optionally over the model population only.

    The manuscript's own robustness check: the measured winning cell may lie at
    an argument below the fit's region-A start, but the winner should not
    change when both measured stores are restricted to x >= kTierThresholds[0].
    """
    ship_abs, ship_rel, prod_abs, prod_rel = {}, {}, {}, {}
    for band, x_text, x, n, recursed, direct in cells:
        if model_from is not None and x < model_from:
            continue
        expected = reference[(band, x)][n]
        rec_error = abs(mpf(recursed) - expected)
        dir_error = abs(mpf(direct) - expected)
        rec_relative = rec_error / expected if expected != 0 else mpf("inf")
        dir_relative = dir_error / expected if expected != 0 else mpf("inf")
        A.update(ship_abs, (band, n), rec_error, x_text)
        A.update(prod_abs, (band, n), dir_error, x_text)
        A.update(ship_rel, (band, n), rec_relative, x_text)
        A.update(prod_rel, (band, n), dir_relative, x_text)

    winner = {}
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
    return winner


def commit_map(table):
    out = {}
    for key, rows in table.items():
        for row in rows:
            if row["state"] != "abstain":
                out[(key[0], key[1], row["n"])] = row["choice"]
    return out



def all_rows(rec_abs, prod_abs, arith_abs, rec_rel, prod_rel, arith_rel, winner, margin):
    out = {}
    for band in ("A1", "A2"):
        for metric in ("absolute", "relative"):
            if metric == "absolute":
                out[(band, metric)] = rows_of(rec_abs, prod_abs, arith_abs, winner, band, metric, margin)
            else:
                out[(band, metric)] = rows_of(rec_rel, prod_rel, arith_rel, winner, band, metric, margin)
    return out


def score_all(table):
    total = [0, 0, 0]
    per = {}
    for key, rows in table.items():
        s = score(rows)
        per[key] = s
        total = [total[i] + s[i] for i in range(3)]
    return total, per


def frontier(rec_abs, prod_abs, arith_abs, rec_rel, prod_rel, arith_rel, winner):
    """The exact commit frontier as the margin falls from +infinity to 0.

    A comparison with sigma == 0 and gap > 0 commits at every margin; one with
    gap == 0 abstains at every margin.  For sigma > 0 the critical margin is
    r = gap / sigma: the comparison commits iff margin < r.  Grouping the
    critical values in descending order gives every attainable commit count
    and its right/wrong split; no margin between two critical values changes
    the set.
    """
    items = []
    always_commit = 0
    always_abstain = 0
    always_right = 0
    for band in ("A1", "A2"):
        for metric, rec, prod, arith in (
            ("absolute", rec_abs, prod_abs, arith_abs),
            ("relative", rec_rel, prod_rel, arith_rel),
        ):
            for n in range(NMAX + 1):
                ps, pd, sigma = rec[band][n], prod[band][n], arith[band][n]
                gap = abs(ps - pd)
                choice = "ship" if ps < pd else "product"
                right = choice == winner[(band, metric, n)]
                if sigma == 0:
                    if gap > 0:
                        always_commit += 1
                        if right:
                            always_right += 1
                    else:
                        always_abstain += 1
                    continue
                r = gap / sigma
                if gap > 0:
                    items.append((r, band, metric, n, choice, right))
                else:
                    always_abstain += 1

    items.sort(key=lambda item: item[0], reverse=True)
    blocks = []
    for item in items:
        if blocks and item[0] == blocks[-1][0][0]:
            blocks[-1].append(item)
        else:
            blocks.append([item])

    steps = []  # (margin_low, margin_high, commits, right, wrong)
    commits = always_commit
    right_total = always_right
    high = mpf("inf")
    for block in blocks:
        low = block[0][0]
        steps.append((low, high, commits, right_total, commits - right_total))
        for _r, _b, _m, _n, _c, right in block:
            commits += 1
            if right:
                right_total += 1
        high = low
    steps.append((mpf(0), high, commits, right_total, commits - right_total))
    return steps, always_commit, always_abstain, blocks


def describe_frontier(steps, always_commit, always_abstain):
    print("  always-commit comparisons (sigma = 0, gap > 0): %d" % always_commit)
    print("  always-abstain comparisons (gap = 0): %d" % always_abstain)
    print("  attainable counts run from %d to %d" % (steps[0][2], steps[-1][2]))


def find_coverage(steps, target):
    """Locate a target commit count in the frontier.

    Returns (kind, info): kind 'exact' with the bracket (low, high) of margins
    that give exactly `target`; or 'skipped' with the two neighbouring
    margin values and the counts either side.
    """
    for low, high, commits, right, wrong in steps:
        if commits == target:
            return "exact", (low, high, commits, right, wrong)
    for index in range(len(steps) - 1):
        if steps[index][2] > target > steps[index + 1][2]:
            return "skipped", (steps[index], steps[index + 1])
    return "absent", None


def report_frontier(name, steps, always_commit, always_abstain, target=67, full=False):
    print()
    print("=== %s ===" % name)
    describe_frontier(steps, always_commit, always_abstain)
    kind, info = find_coverage(steps, target)
    print("  target %d:" % target, end=" ")
    if kind == "exact":
        low, high, commits, right, wrong = info
        print("attained on margin [%s, %s)" % (text(low), text(high)))
        print("    at that margin: commits %d, right %d, wrong %d" % (commits, right, wrong))
        idx = next(i for i, s in enumerate(steps) if s[2] == target)
        if idx + 1 < len(steps):
            lo2, hi2, c2, r2, w2 = steps[idx + 1]
            print(
                "    next higher-coverage margin [%s, %s): commits %d, right %d, wrong %d"
                % (text(lo2), text(hi2), c2, r2, w2)
            )
        if idx > 0:
            lo1, hi1, c1, r1, w1 = steps[idx - 1]
            print(
                "    next lower-coverage margin [%s, %s): commits %d, right %d, wrong %d"
                % (text(lo1), text(hi1), c1, r1, w1)
            )
    elif kind == "skipped":
        above, below = info
        print("skipped between two blocks")
        print(
            "    margin at %s: commits %d, right %d, wrong %d"
            % (text(above[1]), above[2], above[3], above[4])
        )
        print(
            "    margin just below %s: commits %d, right %d, wrong %d"
            % (text(above[1]), below[2], below[3], below[4])
        )
    else:
        print("not reachable")

    if full:
        print("  full frontier (commits, right, wrong) per margin interval:")
        for low, high, commits, right, wrong in reversed(steps):
            print(
                "    [%s, %s) commits %3d right %3d wrong %3d"
                % (text(low), text(high), commits, right, wrong)
            )
    else:
        print("  frontier sampled (commits, right, wrong) as the margin falls:")
        keep = sorted({s[2] for s in steps})
        for count in keep:
            entry = next(s for s in steps if s[2] == count)
            print(
                "    commits %3d right %3d wrong %3d  (margin in [%s, %s))"
                % (entry[2], entry[3], entry[4], text(entry[0]), text(entry[1]))
            )


def main():
    reference, cells, error0, winner, _sa, _sr, _pa, _pr = build_measured()
    entries, population, model_from = build_cells(reference, cells, error0)
    print("records:")
    print("  cells.csv rows: %d" % len(cells))
    print("  reference.csv rows: %d" % len(reference))
    print("  model population x >= %r: A1 %d, A2 %d" % (model_from, population["A1"], population["A2"]))
    print("  population cells (entries): A1 %d, A2 %d" % (len(entries["A1"]), len(entries["A2"])))

    # ---- the true model, verbatim -----------------------------------------
    true_abs = predicted(entries, A.amplification, "absolute")
    true_rel = predicted(entries, A.amplification, "relative")
    table = all_rows(*true_abs, *true_rel, winner, mpf(1))
    total, per = score_all(table)
    print()
    print("true model at margin 1 (A_n unchanged):")
    for key in sorted(per):
        print("  %s %-8s commits %d, right %d, wrong %d" % (key[0], key[1], *per[key]))
    print("  total commits %d, right %d, wrong %d" % tuple(total))

    # ---- cross-check against the stored prediction.csv --------------------
    stored = {}
    with open(TREE / "results/prediction.csv", newline="") as handle:
        for row in csv.DictReader(handle):
            stored[(row["band"], row["metric"], int(row["order"]))] = (
                row["state"],
                row["choice"],
                row["right_where_committed"],
            )
    mismatches = 0
    for key, rows in table.items():
        for row in rows:
            s = stored[(key[0], key[1], row["n"])]
            mine = (
                row["state"],
                row["choice"],
                "1" if row["right"] else "0",
            )
            if s != mine:
                mismatches += 1
                if mismatches <= 10:
                    print("  MISMATCH %s n=%d stored %r mine %r" % (key, row["n"], s, mine))
    print("  stored prediction.csv: %d rows; %d mismatches" % (len(stored), mismatches))

    # ---- the true model's own margin frontier, for context ----------------
    true_steps, true_ac, true_aa, _true_blocks = frontier(*true_abs, *true_rel, winner)
    true_zero = [s[2] for s in true_steps if s[4] == 0]
    print()
    print("true model's own retuned frontier (context, not the control):")
    print("  zero-error commit range: %d .. %d" % (min(true_zero), max(true_zero)))
    for low, high, commits, right, wrong in true_steps:
        if commits in (67, 68, 80, 90, 100):
            print(
                "    commits %3d right %3d wrong %3d  margin in [%s, %s)"
                % (commits, right, wrong, text(low), text(high))
            )
    with open(Path(__file__).resolve().parent / "frontier_true.csv", "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["margin_low", "margin_high", "commits", "right", "wrong"])
        for low, high, commits, right, wrong in true_steps:
            writer.writerow([text(low), text(high), commits, right, wrong])
    for label, m in (
        ("margin 1", mpf(1)),
        ("margin 0.94", mpf("0.94")),
        ("upper edge of 74 interval", mpf("0.94599162174460234")),
    ):
        total, _per = score_all(all_rows(*true_abs, *true_rel, winner, m))
        print(
            "    direct %-28s margin %-22s commits %3d right %3d wrong %3d"
            % (label, text(m), *total)
        )

    # ---- the monotone surrogate: A_n frozen at its minimum ---------------
    def nmin_of(x):
        return min(NMAX, max(0, math.ceil(x - 0.5)))

    cache = {}

    def frozen_min(n, x):
        m = nmin_of(x)
        # The paper's upward-only surrogate: keep A_n on the falling branch
        # up to n_min(x), then hold A_{n_min} at every higher order.
        if n <= m:
            return A.amplification(n, x)
        key = (m, x)
        if key not in cache:
            cache[key] = A.amplification(m, x)
        return cache[key]

    sur_abs = predicted(entries, frozen_min, "absolute")
    sur_rel = predicted(entries, frozen_min, "relative")
    sur_table = all_rows(*sur_abs, *sur_rel, winner, mpf(1))
    sur_total, sur_per = score_all(sur_table)
    print()
    print("monotone surrogate (A_n frozen at the minimum A_nmin(x)) at margin 1:")
    for key in sorted(sur_per):
        print("  %s %-8s commits %d, right %d, wrong %d" % (key[0], key[1], *sur_per[key]))
    print("  total commits %d, right %d, wrong %d" % tuple(sur_total))

    # the 67 -> 46 accounting
    true_set = {
        (key[0], key[1], row["n"])
        for key, rows in table.items()
        for row in rows
        if row["state"] != "abstain"
    }
    sur_set = {
        (key[0], key[1], row["n"])
        for key, rows in sur_table.items()
        for row in rows
        if row["state"] != "abstain"
    }
    lost = true_set - sur_set
    gained = sur_set - true_set
    print("  true commitments %d, surrogate commitments %d" % (len(true_set), len(sur_set)))
    print("  true commitments lost to abstention: %d" % len(lost))
    print("  surrogate commitments outside the true set: %d" % len(gained))
    wrong_lost = [
        key
        for key in lost
        if next(r for r in table[(key[0], key[1])] if r["n"] == key[2])["right"] is False
    ]
    print("  of the lost, true model names a wrong winner at: %d" % len(wrong_lost))
    for key in sorted(lost):
        row = next(r for r in table[(key[0], key[1])] if r["n"] == key[2])
        print(
            "    lost %s n=%2d true choice %-7s measured %-7s right %s"
            % (key[0], key[2], row["choice"], row["measured"], row["right"])
        )

    steps, always_commit, always_abstain, blocks = frontier(*sur_abs, *sur_rel, winner)
    report_frontier(
        "monotone surrogate, A_n frozen at the minimum",
        steps,
        always_commit,
        always_abstain,
        target=67,
        full=True,
    )

    # ---- direct verification at the critical margins ----------------------
    print()
    print("direct scoring at the critical margins (the frontier's own check):")

    def direct(margin, winner_map=winner):
        t = all_rows(*sur_abs, *sur_rel, winner_map, mpf(margin))
        return score_all(t)

    idx67 = next(i for i, s in enumerate(steps) if s[2] == 67)
    low67, high67 = steps[idx67][0], steps[idx67][1]
    low68, high68 = steps[idx67 + 1][0], steps[idx67 + 1][1]
    low66, high66 = steps[idx67 - 1][0], steps[idx67 - 1][1]
    checks = [
        ("margin 1", mpf(1)),
        ("lower neighbour, commits 68", (low68 + high68) / 2),
        ("exact lower edge of 67", low67),
        ("interior of 67", (low67 + high67) / 2),
        ("exact upper edge of 67", high67),
        ("upper neighbour, commits 66", (low66 + high66) / 2),
    ]
    for label, m in checks:
        total, per = direct(m)
        print(
            "  %-30s margin %-24s commits %3d right %3d wrong %3d"
            % (label, text(m), *total)
        )
    with open(Path(__file__).resolve().parent / "frontier_frozen_min.csv", "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["margin_low", "margin_high", "commits", "right", "wrong"])
        for low, high, commits, right, wrong in steps:
            writer.writerow([text(low), text(high), commits, right, wrong])

    # ---- the retuned set against the true model's set ---------------------
    s67 = commit_map(all_rows(*sur_abs, *sur_rel, winner, mpf("0.946")))
    print()
    print("surrogate at margin 0.946 vs the true model at margin 1:")
    print("  surrogate commitments %d, true commitments %d" % (len(s67), len(true_set)))
    print("  common %d" % len(set(s67) & true_set))
    print("  surrogate only %d" % len(set(s67) - true_set))
    print("  true only %d" % len(true_set - set(s67)))
    flipped = []
    for k in sorted(set(s67) & true_set):
        row = next(r for r in table[(k[0], k[1])] if r["n"] == k[2])
        if s67[k] != row["choice"]:
            flipped.append(k)
    print("  common commitments where the named winner differs: %d" % len(flipped))
    for k in sorted(flipped):
        row = next(r for r in table[(k[0], k[1])] if r["n"] == k[2])
        print(
            "    %s %-8s n=%2d surrogate names %-7s true names %-7s measured %-7s"
            % (k[0], k[1], k[2], s67[k], row["choice"], row["measured"])
        )
    print("  comparisons the retuned surrogate commits on but the true model abstains:")
    for k in sorted(set(s67) - true_set):
        print("    %s %-8s n=%2d surrogate names %s" % (k[0], k[1], k[2], s67[k]))
    print("  comparisons the true model commits on but the retuned surrogate abstains:")
    for k in sorted(true_set - set(s67)):
        row = next(r for r in table[(k[0], k[1])] if r["n"] == k[2])
        print("    %s %-8s n=%2d true names %s" % (k[0], k[1], k[2], row["choice"]))
    retuned_per = score_all(all_rows(*sur_abs, *sur_rel, winner, mpf("0.946")))[1]
    print("  per group at margin 0.946:")
    for key in sorted(retuned_per):
        print("    %s %-8s commits %d, right %d, wrong %d" % (key[0], key[1], *retuned_per[key]))

    # ---- the manuscript's population-restricted winner robustness ---------
    winner_pop = winners_over(reference, cells, model_from)
    t_pop, _ = score_all(all_rows(*true_abs, *true_rel, winner_pop, mpf(1)))
    s_pop, _ = direct(mpf("0.946"), winner_pop)
    differing = [
        key
        for key in true_set
        if winner_pop[(key[0], key[1], key[2])] != winner[(key[0], key[1], key[2])]
    ]
    print()
    print("measured winners restricted to x >= %r:" % model_from)
    print("  true model at margin 1: commits %d right %d wrong %d" % tuple(t_pop))
    print("  surrogate at margin 0.946: commits %d right %d wrong %d" % tuple(s_pop))
    print("  orders where the restricted winner differs: %d" % len(differing))
    for key in sorted(differing):
        print(
            "    %s %-8s n=%2d band %-7s population %-7s"
            % (key[0], key[1], key[2], winner[(key[0], key[1], key[2])], winner_pop[(key[0], key[1], key[2])])
        )

    # ---- the frozen-first surrogate: A_n = A_1 = 1/(2x) -------------------
    def frozen_first(n, x):
        return mpf(1) / (2 * x)

    def frozen_first_n0(n, x):
        return mpf(1) if n == 0 else mpf(1) / (2 * x)

    for label, fun in (
        ("frozen-first, every order including n = 0", frozen_first),
        ("frozen-first, n = 0 kept at A_0 = 1", frozen_first_n0),
    ):
        f_abs = predicted(entries, fun, "absolute")
        f_rel = predicted(entries, fun, "relative")
        f_table = all_rows(*f_abs, *f_rel, winner, mpf(1))
        f_total, f_per = score_all(f_table)
        print()
        print("%s at margin 1:" % label)
        for key in sorted(f_per):
            print("  %s %-8s commits %d, right %d, wrong %d" % (key[0], key[1], *f_per[key]))
        print("  total commits %d, right %d, wrong %d" % tuple(f_total))
        steps2, ac2, aa2, _blocks2 = frontier(*f_abs, *f_rel, winner)
        report_frontier(label, steps2, ac2, aa2, target=67, full=False)
        if fun is frozen_first_n0:
            print("  direct scoring at the critical margins:")
            for lab, mm in (
                ("margin 1", mpf(1)),
                ("inside 67 interval", (steps2[next(i for i, s in enumerate(steps2) if s[2] == 67)][0]
                                        + steps2[next(i for i, s in enumerate(steps2) if s[2] == 67)][1]) / 2),
                ("largest zero-error margin", mpf("100")),
            ):
                t = all_rows(*f_abs, *f_rel, winner, mm)
                tot, _p = score_all(t)
                print(
                    "    %-26s margin %-24s commits %3d right %3d wrong %3d"
                    % (lab, text(mm), *tot)
                )
            with open(
                Path(__file__).resolve().parent / "frontier_frozen_first.csv", "w", newline=""
            ) as handle:
                writer = csv.writer(handle)
                writer.writerow(["margin_low", "margin_high", "commits", "right", "wrong"])
                for low, high, commits, right, wrong in steps2:
                    writer.writerow([text(low), text(high), commits, right, wrong])

    # ---- the A_n = 1 control, for pipeline validation ---------------------
    def unit(n, x):
        return mpf(1)

    u_abs = predicted(entries, unit, "absolute")
    u_rel = predicted(entries, unit, "relative")
    u_table = all_rows(*u_abs, *u_rel, winner, mpf(1))
    u_total, u_per = score_all(u_table)
    print()
    print("A_n = 1 control at margin 1 (manuscript: 98 commits, 27 wrong):")
    for key in sorted(u_per):
        print("  %s %-8s commits %d, right %d, wrong %d" % (key[0], key[1], *u_per[key]))
    print("  total commits %d, right %d, wrong %d" % tuple(u_total))

    return 0


if __name__ == "__main__":
    sys.exit(main())
