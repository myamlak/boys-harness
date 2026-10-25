#!/usr/bin/env python3
"""Analyse the Paper-A workload pair probe (scratch copy only).

Reads the six arms' ``.args.txt`` and ``.summary.txt`` files written by
``qcx-pair-probe``, rebuilds the model's margin-1 states from the paper's own
``prediction.csv``, checks every emitted ladder against the selection it should
be, and compares the emitted values with an 80-digit reference.

The comparison unit is the paper's: a (band, order) cell, and inside it the
worst absolute and relative error over the workload's distinct arguments. The
six arms are read in lockstep, one argument at a time, so the memory does not
grow with the basis.
"""

import csv
import sys
from pathlib import Path

if sys.flags.isolated:
    print("SKIP: isolated mode hides user-site mpmath; run without python -I")
    sys.exit(77)

from mpmath import mp, mpf, sqrt, erf, exp  # noqa: E402

mp.dps = 80

PAPER = Path(r"C:/Users/myaml/dev/papers-a-v5/papers/boys-crossover")
PREDICTION = PAPER / "sweeps" / "amplification_form" / "results" / "prediction.csv"

A1_EDGE = 5.94992407605424223
A_END = 11.899848152108484
THETA0 = 1.08552523453493333
MARGIN = mpf(1)

ARMS = ["paper-a", "direct", "model-abs-band", "model-abs-pop", "model-rel-band",
        "model-rel-pop"]
BAND_ARMS = {
    "paper-a": ("paper-a", False),
    "direct": ("direct", False),
    "model-abs-band": ("absolute", False),
    "model-abs-pop": ("absolute", True),
    "model-rel-band": ("relative", False),
    "model-rel-pop": ("relative", True),
}


def model_states():
    """The paper's margin-1 state per band, metric and order."""
    table = {}
    with open(PREDICTION, newline="") as handle:
        for row in csv.DictReader(handle):
            ship = mpf(row["predicted_ship"])
            product = mpf(row["predicted_product"])
            arithmetic = mpf(row["arithmetic"])
            gap = abs(ship - product)
            if gap <= MARGIN * arithmetic:
                state = "abstain"
            else:
                state = "ship" if ship < product else "product"
            stored = row["state"]
            if stored and stored != state:
                raise SystemExit("state mismatch at %s %s %s: %s vs %s"
                                 % (row["band"], row["metric"], row["order"], state, stored))
            table[(row["band"], row["metric"], int(row["order"]))] = state
    if len(table) != 132:
        raise SystemExit("expected 132 prediction rows")
    return table


def reference_ladder(nmax, x_float):
    """F_0(x)..F_nmax(x) at 80 digits.

    The series below x = 2 (no cancellation there); above it the stable upward
    recurrence seeded by F_0 = sqrt(pi)/(2 sqrt(x)) erf(sqrt(x)).
    """
    x = mpf(x_float)
    values = []

    if x < 2:
        for n in range(nmax + 1):
            term = mpf(1) / (2 * n + 1)
            total = term
            k = 0
            while True:
                k += 1
                term = term * (-x) * (2 * n + 2 * k - 1) / (k * (2 * n + 2 * k + 1))
                total += term
                if abs(term) <= mpf("1e-110") * max(abs(total), mpf("1e-300")):
                    break
                if k > 100000:
                    raise SystemExit("series did not converge")
            values.append(total)
        return values

    f0 = sqrt(mp.pi) / (2 * sqrt(x)) * erf(sqrt(x))
    values.append(f0)
    e = exp(-x)
    for n in range(nmax):
        values.append(((2 * n + 1) * values[-1] - e) / (2 * x))
    return values


def reference_by_quad(nmax, x_float, orders):
    """The same ladder from the integral definition, by quadrature.

    F_n(x) = integral_0^1 t^(2n) exp(-x t^2) dt. Used only to validate the
    series/recurrence ladder; the production reference is the ladder.
    """
    x = mpf(x_float)
    out = {}
    for n in orders:
        out[n] = mp.quad(lambda t, n=n: t ** (2 * n) * exp(-x * t * t), [0, 1])
    return out


def validate_reference():
    """Check the ladder against quadrature at sample arguments and orders."""
    sample_x = [0.0, 1e-6, 1e-2, 0.5, 1.0, 2.0, 6.0, 12.0, 50.0,
                277.20922853176552]
    sample_orders = (0, 3, 6, 9, 12)
    worst = mpf(0)
    where = None
    for x in sample_x:
        ladder = reference_ladder(12, x)
        quads = reference_by_quad(12, x, sample_orders)
        for order, value in quads.items():
            deviation = abs(ladder[order] - value) / abs(value)
            if deviation > worst:
                worst = deviation
                where = (x, order)
    return worst, where


def band_of(x):
    if x < A1_EDGE:
        return "A1"
    if x < A_END:
        return "A2"
    return "above"


def stream(path):
    with open(path) as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            tokens = line.split()
            nmax = int(tokens[0])
            count = int(tokens[1])
            x = float(tokens[2])
            width = nmax + 1
            paper_a = [float(t) for t in tokens[3:3 + width]]
            direct = [float(t) for t in tokens[3 + width:3 + 2 * width]]
            engine = [float(t) for t in tokens[3 + 2 * width:3 + 3 * width]]
            yield nmax, count, x, paper_a, direct, engine


def selected_column(arm, x, order, states):
    band = band_of(x)
    mode, population = BAND_ARMS[arm]

    if mode == "paper-a":
        return "paper-a"
    if mode == "direct":
        return "direct"

    # Outside region A the direct entry has no lane: every arm emits the named
    # policy's value.
    if band == "above":
        return "paper-a"
    if population and x < THETA0:
        return "direct"

    state = states[(band, mode, order)]
    if state == "ship":
        return "paper-a"
    if state == "product":
        return "direct"
    return "direct"  # the named fallback


def main():
    if len(sys.argv) != 3:
        raise SystemExit("usage: pair_analyze.py <probe-prefix> <out-file>")
    prefix = sys.argv[1]
    out_path = Path(sys.argv[2])

    states = model_states()
    validation = validate_reference()

    handles = {arm: stream("%s.%s.args.txt" % (prefix, arm)) for arm in ARMS}
    cells = {arm: {} for arm in ARMS}
    above = {arm: [mpf(0), None] for arm in ARMS}
    arguments = 0
    calls = 0
    bands = {}
    mismatches = {arm: 0 for arm in ARMS}

    while True:
        rows = {}
        exhausted = []
        for arm in ARMS:
            try:
                rows[arm] = next(handles[arm])
            except StopIteration:
                exhausted.append(arm)
        if exhausted:
            if len(exhausted) != len(ARMS):
                raise SystemExit("arms with fewer arguments: " + ", ".join(exhausted))
            break

        base = rows["paper-a"]
        nmax, count, x = base[0], base[1], base[2]
        for arm in ARMS[1:]:
            row = rows[arm]
            if row[0] != nmax or row[2] != x:
                raise SystemExit("argument stream differs in arm " + arm)
            if row[3] != base[3] or row[4] != base[4]:
                raise SystemExit("route columns differ in arm " + arm)

        arguments += 1
        calls += count
        band = band_of(x)
        bands[band] = bands.get(band, 0) + 1
        ref = reference_ladder(nmax, x)

        for arm in ARMS:
            nmax_a, count_a, x_a, paper_a, direct, engine = rows[arm]
            for order in range(nmax_a + 1):
                column = selected_column(arm, x_a, order, states)
                expected = paper_a[order] if column == "paper-a" else direct[order]
                if expected != engine[order]:
                    mismatches[arm] += 1
                value = mpf(engine[order])
                error = abs(value - ref[order])
                relative = error / abs(ref[order])
                if band == "above":
                    if error > above[arm][0]:
                        above[arm] = [error, (order, x_a, relative)]
                    continue
                entry = cells[arm].setdefault((band, order), [mpf(0), mpf(0), None, None])
                if error > entry[0]:
                    entry[0] = error
                    entry[2] = x_a
                if relative > entry[1]:
                    entry[1] = relative
                    entry[3] = x_a

    for arm in ARMS:
        if mismatches[arm]:
            raise SystemExit("%d emitted values in arm %s are not the selected route"
                             % (mismatches[arm], arm))

    lines = []

    def emit(text=""):
        lines.append(text)
        print(text)

    emit("# Workload pair analysis")
    emit()
    emit("reference: 80-digit ladder, series below x = 2 and upward recurrence above;")
    emit("  validated against mp.quad at 10 arguments x 5 orders, worst relative")
    emit("  deviation %s at %s" % (mp.nstr(validation[0], 6), validation[1]))
    emit()
    emit("distinct arguments: %d, engine calls: %d" % (arguments, calls))
    emit("distinct arguments per band: " + ", ".join("%s=%d" % item
                                                     for item in sorted(bands.items())))
    emit()

    for arm in ARMS:
        emit("arm %s: above-region worst abs %s at %s"
             % (arm, mp.nstr(above[arm][0], 6), above[arm][1]))
    emit()

    for arm in ARMS:
        emit("## " + arm)
        emit("  cell  worst_abs (at x)  worst_rel (at x)")
        for key in sorted(cells[arm]):
            entry = cells[arm][key]
            emit("  %s n=%2d  %s (%.6g)  %s (%.6g)"
                 % (key[0], key[1], mp.nstr(entry[0], 8), entry[2],
                    mp.nstr(entry[1], 8), entry[3]))
        emit()

    def worst_over(arm, metric, drop=None):
        index = 0 if metric == "abs" else 1
        values = [(arm, key, entry[index]) for key, entry in cells[arm].items()
                  if key != drop]
        return max(values, key=lambda item: item[2])

    def bar(model_arm, metric):
        direct_arm = "direct"
        worst_direct = worst_over(direct_arm, metric)
        worst_model = worst_over(model_arm, metric)
        margin = (worst_direct[2] - worst_model[2]) / worst_direct[2]
        emit("  %s %s: direct %s at %s; model %s at %s; margin %s"
             % (model_arm, metric, mp.nstr(worst_direct[2], 8), worst_direct[1],
                mp.nstr(worst_model[2], 8), worst_model[1], mp.nstr(margin, 6)))
        loo = []
        for key in sorted(cells[direct_arm]):
            d = worst_over(direct_arm, metric, drop=key)
            m = worst_over(model_arm, metric, drop=key)
            loo.append((key, (d[2] - m[2]) / d[2], d[2], m[2]))
        worst_loo = min(loo, key=lambda item: item[1])
        emit("    leave-one-out smallest margin %s after dropping %s (direct %s, model %s)"
             % (mp.nstr(worst_loo[1], 6), worst_loo[0],
                mp.nstr(worst_loo[2], 8), mp.nstr(worst_loo[3], 8)))
        no_order0 = [key for key in cells[direct_arm] if key[1] != 0]
        if no_order0:
            index = 0 if metric == "abs" else 1
            dm = max(cells[direct_arm][key][index] for key in no_order0)
            mm = max(cells[model_arm][key][index] for key in no_order0)
            emit("    orders >= 1 only: direct %s, model %s, margin %s"
                 % (mp.nstr(dm, 8), mp.nstr(mm, 8), mp.nstr((dm - mm) / dm, 6)))
        deciding = []
        index = 0 if metric == "abs" else 1
        for key in sorted(cells[direct_arm]):
            d = cells[direct_arm][key][index]
            m = cells[model_arm][key][index]
            if m < d:
                deciding.append((key, (d - m) / d))
        emit("    cells where the model beats always-direct: "
             + (", ".join("%s n=%d (%s)" % (key[0], key[1], mp.nstr(marg, 4))
                          for key, marg in deciding) if deciding else "none"))
        emit()

    for metric in ("abs", "rel"):
        emit("## referee's bar: " + ("absolute" if metric == "abs" else "relative"))
        for model_arm in ("model-abs-band", "model-abs-pop", "model-rel-band", "model-rel-pop"):
            bar(model_arm, metric)
        emit()

    # Tensor and energy comparisons, streamed for the same reason.
    pairs = [("paper-a", "direct"), ("paper-a", "model-abs-band"),
             ("direct", "model-abs-band"), ("paper-a", "model-abs-pop"),
             ("direct", "model-abs-pop"), ("paper-a", "model-rel-band"),
             ("direct", "model-rel-band"), ("paper-a", "model-rel-pop"),
             ("direct", "model-rel-pop")]
    eri_handles = {arm: open("%s.%s.eri.txt" % (prefix, arm)) for arm in ARMS}
    counts = {}
    for arm in ARMS:
        counts[arm] = int(eri_handles[arm].readline())
    if len(set(counts.values())) != 1:
        raise SystemExit("the arms disagree on the basis size")
    n_basis = counts[ARMS[0]]
    stats = {pair: {"max": 0.0, "where": None, "sumsq": 0.0, "differing": 0}
             for pair in pairs}
    max_eri = 0.0
    entries = 0

    while True:
        rows = {}
        done = []
        for arm in ARMS:
            line = eri_handles[arm].readline()
            if not line:
                done.append(arm)
                continue
            parts = line.split()
            rows[arm] = (int(parts[0]), int(parts[1]), int(parts[2]), int(parts[3]),
                         float(parts[4]))
        if done:
            if len(done) != len(ARMS):
                raise SystemExit("the arms' ERI dumps differ in length")
            break
        entries += 1
        magnitude = abs(rows["paper-a"][4])
        if magnitude > max_eri:
            max_eri = magnitude
        for pair in pairs:
            left, right = pair
            delta = rows[left][4] - rows[right][4]
            entry = stats[pair]
            if delta != 0.0:
                entry["differing"] += 1
            size = abs(delta)
            if size > entry["max"]:
                entry["max"] = size
                entry["where"] = rows[left][:4]
            entry["sumsq"] += delta * delta

    for handle in eri_handles.values():
        handle.close()
    if entries != n_basis ** 4:
        raise SystemExit("an ERI dump is short")

    emit("## tensors: n=%d, %d entries, max |ERI| = %.6g" % (n_basis, entries, max_eri))
    for pair in pairs:
        left, right = pair
        entry = stats[pair]
        rms = (entry["sumsq"] / entries) ** 0.5
        emit("  %-14s vs %-14s max|d| %.6e at %s rms %.6e differing %d"
             % (left, right, entry["max"], entry["where"], rms, entry["differing"]))
    emit()

    def summary(arm):
        values = {}
        with open("%s.%s.summary.txt" % (prefix, arm)) as handle:
            for line in handle:
                key, _, value = line.strip().partition("=")
                values[key] = value
        return values

    summaries = {arm: summary(arm) for arm in ARMS}
    emit("## fixed-density and SCF energies")
    for arm in ARMS:
        entry = summaries[arm]
        emit("  %-14s seed_e2e %s default_e2e %s tight_e2e %s scf_default %s scf_tight %s"
             % (arm, entry["seed_e2e"], entry["default_e2e"], entry["tight_e2e"],
                entry["scf_default_total_energy"], entry["scf_tight_total_energy"]))

    base_summary = summaries["paper-a"]
    for arm in ARMS[1:]:
        entry = summaries[arm]
        deltas = []
        for key in ("seed_e2e", "default_e2e", "tight_e2e",
                    "scf_default_total_energy", "scf_tight_total_energy"):
            deltas.append("%s %+.6e" % (key, float(entry[key]) - float(base_summary[key])))
        emit("  %-14s vs paper-a: %s" % (arm, ", ".join(deltas)))
    emit()

    emit("## E_2e at the shared seed density, delta vs paper-a")
    for arm in ARMS:
        emit("  %-14s %+.6e" % (arm, float(summaries[arm]["seed_e2e"]) -
                                float(base_summary["seed_e2e"])))

    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
