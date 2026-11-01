#!/usr/bin/env python3
"""Score the crossing model with the exponential's error in the arithmetic term.

The recurrence body evaluates 0.5*exp(-x) once per call and shares it across
the upward steps (boys_impl.hpp, RegionBHalfExp).  A relative error eta_E in
exp(-x) perturbs that shared value by eta_E*exp(-x)/2, and the perturbed step
output is carried to order n by the recurrence's own gain.  The kernel states
the sensitivity in its own comment: with t = exp(-x)/2 and

    T(n, x) = (1/(2x)) sum_{l=0}^{n-1} prod_{j=l+1}^{n-1} (2j+1)/(2x),

dF_n/dt = -2 T(n, x), checked symbolically in the kernel.  The contribution
to the recursed error at order n is therefore eta_E*exp(-x)*T(n, x), the
closed form this script evaluates as

    eta_E * exp(-x) * sum_{l=0}^{n-1} A_n(x)/A_l(x) / (2l+1),

which the script checks against the product form at a sample of orders and
arguments.  eta_E is taken at one unit in the last place of a binary64,
2^-52 relative, the largest relative ulp of a normal double; the C library's
exponential accuracy is environment-bound and is not certified here.

The script reads the paper tree's own measured records (cells, reference,
provenance) and writes nothing outside analysis-strengthen/exp-term/out/.
It first runs the model at eta_E = 0 and compares its states with the
recorded predictions, so the rescaling is validated against the scored
record before the exponential term is switched on.

Run from this directory:
    python score_exp.py
"""

import csv
import importlib.util
import json
import sys
from pathlib import Path
from mpmath import mp, mpf

HERE = Path(__file__).resolve().parent
STRENGTHEN = HERE.parent
PAPER = STRENGTHEN.parent
SWEEPS = PAPER / "sweeps"
OUT = HERE / "out"
NMAX = 32
MARGINS = [0.0, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0]
ETA = mpf(2) ** -52
THREE_U = 3 * mpf(2) ** -53


def load_analysis():
    spec = importlib.util.spec_from_file_location(
        "paper_a_analyze", SWEEPS / "amplification_form" / "analyze.py"
    )
    module = importlib.util.module_from_spec(spec)
    sys.path.insert(0, str(SWEEPS / "common"))
    sys.path.insert(0, str(SWEEPS / "amplification_form"))
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


class Model:
    def __init__(self, module, cells_path, reference_path, provenance_path, eta):
        self.eta = eta
        self.reference = module.read_reference(str(reference_path))
        self.cells = read_cells(cells_path)
        constants = json.load(open(provenance_path, encoding="utf-8"))["constants"]
        self.thresholds = [float(v) for v in constants["kTierThresholds"]]
        self.model_from = self.thresholds[0]
        self.tables = {}

    def table(self, band, x):
        key = (band, x)
        entry = self.tables.get(key)
        if entry is None:
            xm = mpf(x)
            gains = [mpf(1)]
            for l in range(NMAX):
                gains.append(gains[-1] * (2 * l + 1) / (2 * xm))
            entry = (gains, mp.exp(-xm))
            self.tables[key] = entry
        return entry

    def arithmetic(self, n, band, x):
        """The rounding charge and the exponential term, absolute."""
        gains, ex = self.table(band, x)
        refs = self.reference[(band, x)]
        rounding = mpf(0)
        for l in range(1, n + 1):
            rounding += THREE_U * refs[l] * gains[n] / gains[l]
        carried = mpf(0)
        for l in range(n):
            carried += gains[n] / (gains[l] * (2 * l + 1))
        eta = self.eta(x) if callable(self.eta) else self.eta
        return rounding, eta * ex * carried

    def band_inputs(self):
        error0 = {}
        for band, x_text, x, n, recursed, direct in self.cells:
            if n == 0:
                expected = self.reference[(band, x)][0]
                error0[(band, x)] = abs(mpf(recursed) - expected)

        rec_abs = {"A1": [(mpf(0), None)] * (NMAX + 1), "A2": [(mpf(0), None)] * (NMAX + 1)}
        rec_rel = {"A1": [(mpf(0), None)] * (NMAX + 1), "A2": [(mpf(0), None)] * (NMAX + 1)}
        dir_abs = {"A1": [(mpf(0), None)] * (NMAX + 1), "A2": [(mpf(0), None)] * (NMAX + 1)}
        dir_rel = {"A1": [(mpf(0), None)] * (NMAX + 1), "A2": [(mpf(0), None)] * (NMAX + 1)}
        ar_abs = {"A1": [mpf(0)] * (NMAX + 1), "A2": [mpf(0)] * (NMAX + 1)}
        ar_rel = {"A1": [mpf(0)] * (NMAX + 1), "A2": [mpf(0)] * (NMAX + 1)}
        for band, x_text, x, n, recursed, direct in self.cells:
            if x < self.model_from:
                continue
            seed = error0.get((band, x))
            if seed is None:
                continue
            expected = self.reference[(band, x)][n]
            da = abs(mpf(direct) - expected)
            dr = da / expected if expected != 0 else mpf("inf")
            if n == 0 or x >= self.thresholds[n]:
                rounding, expt = self.arithmetic(n, band, x)
                arith = rounding + expt
                gains, _ = self.table(band, x)
                ra = seed * gains[n] + arith
                rr = ra / expected if expected != 0 else mpf("inf")
                a_rel = arith / expected if expected != 0 else mpf("inf")
            else:
                arith = mpf(0)
                ra, rr, a_rel = da, dr, mpf(0)
            self._push(rec_abs, ar_abs, (band, n), ra, arith, x_text)
            self._push(rec_rel, ar_rel, (band, n), rr, a_rel, x_text)
            self._push(dir_abs, None, (band, n), da, None, x_text)
            self._push(dir_rel, None, (band, n), dr, None, x_text)
        return {
            "rec_abs": rec_abs, "rec_rel": rec_rel,
            "dir_abs": dir_abs, "dir_rel": dir_rel,
            "ar_abs": ar_abs, "ar_rel": ar_rel,
        }

    @staticmethod
    def _push(store, extra, key, value, other, x_text):
        band, n = key
        current = store[band][n]
        if current[0] is None or value > current[0]:
            store[band][n] = (value, x_text)
        if extra is not None and other > extra[band][n]:
            extra[band][n] = other

    def rows(self, inputs, band, metric, margin, winner):
        if metric == "absolute":
            ship, product = inputs["rec_abs"][band], inputs["dir_abs"][band]
            sigma = inputs["ar_abs"][band]
        else:
            ship, product = inputs["rec_rel"][band], inputs["dir_rel"][band]
            sigma = inputs["ar_rel"][band]
        rows = []
        for n in range(NMAX + 1):
            ps, pd = ship[n][0], product[n][0]
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


def winners(model):
    store_abs, store_rel = {}, {}
    winner = {}
    for band, x_text, x, n, recursed, direct in model.cells:
        expected = model.reference[(band, x)][n]
        rec_error = abs(mpf(recursed) - expected)
        dir_error = abs(mpf(direct) - expected)
        rec_relative = rec_error / expected if expected != 0 else mpf("inf")
        dir_relative = dir_error / expected if expected != 0 else mpf("inf")
        for metric, rv, dv in (("absolute", rec_error, dir_error),
                               ("relative", rec_relative, dir_relative)):
            key = (band, n, metric)
            r, d = store_abs.get(key, (None, None))
            if r is None or rv > r:
                r = rv
            if d is None or dv > d:
                d = dv
            store_abs[key] = (r, d)
    for band in ("A1", "A2"):
        for n in range(NMAX + 1):
            for metric in ("absolute", "relative"):
                r, d = store_abs[(band, n, metric)]
                winner[(band, metric, n)] = (
                    "ship" if r < d else ("product" if d < r else "tie")
                )
    return winner


def score(model, inputs, winner, margin):
    aggregate = {}
    per = {}
    for band in ("A1", "A2"):
        for metric in ("absolute", "relative"):
            rows = model.rows(inputs, band, metric, margin, winner)
            commits = sum(1 for r in rows if r[4] != "abstain")
            right = sum(1 for r in rows if r[4] != "abstain" and r[7])
            per[(band, metric)] = (commits, right, len(rows) - commits, rows)
            aggregate.setdefault("commits", 0)
            aggregate["commits"] += commits
            aggregate["right"] = aggregate.get("right", 0) + right
            aggregate["abstains"] = aggregate.get("abstains", 0) + (len(rows) - commits)
    return aggregate, per


def recorded_states(path):
    states = {}
    with open(path, newline="") as handle:
        for row in csv.DictReader(handle):
            states[(row["band"], row["metric"], int(row["order"]))] = (
                row["state"], row["choice"]
            )
    return states


def main():
    mp.dps = 70
    module = load_analysis()
    OUT.mkdir(parents=True, exist_ok=True)
    log = []

    check = verify_closed_form(module)
    log.append(check)

    trees = {
        "second": (
            SWEEPS / "amplification_form" / "results" / "cells.csv",
            SWEEPS / "amplification_form" / "data" / "reference.csv",
            SWEEPS / "amplification_form" / "results" / "provenance.json",
            STRENGTHEN / "out" / "second" / "prediction.csv",
        ),
        "first": (
            STRENGTHEN / "raw" / "cells.csv",
            STRENGTHEN / "raw" / "reference.csv",
            SWEEPS / "product_amplification" / "results" / "provenance.json",
            STRENGTHEN / "out" / "first" / "prediction.csv",
        ),
    }

    for name, (cells_path, reference_path, provenance_path, recorded) in trees.items():
        log.append("")
        log.append("=== %s experiment" % name)
        base = Model(module, cells_path, reference_path, provenance_path, mpf(0))
        base_inputs = base.band_inputs()
        base_winner = winners(base)
        base_scores = {m: score(base, base_inputs, base_winner, m) for m in MARGINS}
        log.append("cells: %d, model arguments: %d" % (len(base.cells), len(base.tables)))

        states = recorded_states(recorded)
        mismatched = 0
        for band in ("A1", "A2"):
            for metric in ("absolute", "relative"):
                rows = base_scores[1.0][1][(band, metric)][3]
                for n, ps, pd, sigma, state, measured, choice, right in rows:
                    if states[(band, metric, n)] != (state, choice):
                        mismatched += 1
        log.append("validation at eta=0: %d of %d states differ from %s"
                   % (mismatched, 4 * (NMAX + 1), recorded.name))

        ext = Model(module, cells_path, reference_path, provenance_path, ETA)
        ext_inputs = ext.band_inputs()
        ext_winner = winners(ext)
        ext_scores = {m: score(ext, ext_inputs, ext_winner, m) for m in MARGINS}

        log.append("score at margin 1.0, eta=0   : commits %d, right %d, abstains %d"
                   % (base_scores[1.0][0]["commits"], base_scores[1.0][0]["right"],
                      base_scores[1.0][0]["abstains"]))
        log.append("score at margin 1.0, eta=2^-52: commits %d, right %d, abstains %d"
                   % (ext_scores[1.0][0]["commits"], ext_scores[1.0][0]["right"],
                      ext_scores[1.0][0]["abstains"]))

        log.append("per band-metric at margin 1.0 (eta=0 -> eta=2^-52):")
        for key in sorted(base_scores[1.0][1]):
            b = base_scores[1.0][1][key]
            e = ext_scores[1.0][1][key]
            changed = [(b[3][n][0], b[3][n][4], e[3][n][4])
                       for n in range(NMAX + 1) if b[3][n][4] != e[3][n][4]]
            log.append("  %s %-8s %2d/%2d/%2d -> %2d/%2d/%2d; changed rows %d %s"
                       % (key[0], key[1], b[0], b[1], b[2], e[0], e[1], e[2],
                          len(changed), changed))

        log.append("margin ladder, aggregate commits/right/abstains (eta=0 -> eta=2^-52):")
        for m in MARGINS:
            b = base_scores[m][0]
            e = ext_scores[m][0]
            log.append("  margin %-4s %3d/%3d/%3d -> %3d/%3d/%3d"
                       % (m, b["commits"], b["right"], b["abstains"],
                          e["commits"], e["right"], e["abstains"]))

        log.append("arithmetic widening sigma_new/sigma_old over the 33 orders:")
        for band in ("A1", "A2"):
            for metric in ("absolute", "relative"):
                sigma_old = base_inputs["ar_abs" if metric == "absolute" else "ar_rel"][band]
                sigma_new = ext_inputs["ar_abs" if metric == "absolute" else "ar_rel"][band]
                ratios = []
                for n in range(NMAX + 1):
                    if sigma_old[n] > 0:
                        ratios.append((sigma_new[n] / sigma_old[n], n,
                                       sigma_old[n], sigma_new[n]))
                if not ratios:
                    log.append("  %s %-8s no nonzero terms" % (band, metric))
                    continue
                lo = min(ratios)
                hi = max(ratios)
                log.append("  %s %-8s min %.9f at n=%d, max %.9f at n=%d"
                           % (band, metric, lo[0], lo[1], hi[0], hi[1]))

    OUT.mkdir(exist_ok=True)
    text = "\n".join(log) + "\n"
    (OUT / "score-exp.log").write_text(text, encoding="utf-8", newline="\n")
    sys.stdout.write(text)


def verify_closed_form(module):
    """The A-ratio sum against the kernel's product form, at a sample."""
    worst = mpf(0)
    where = None
    for n in (1, 2, 5, 9, 16, 24, 32):
        for x in (mpf("0.3"), mpf("1.0855252345349333"), mpf("4.9"),
                  mpf("11.2"), mpf("11.202305053129882"), mpf("30")):
            product = mpf(0)
            for l in range(n):
                inner = mpf(1)
                for j in range(l + 1, n):
                    inner *= (2 * j + 1) / (2 * x)
                product += inner / (2 * x)
            ratios = mpf(0)
            for l in range(n):
                ratios += module.amplification_ratio(n, l, x) / (2 * l + 1)
            if product == 0:
                continue
            dev = abs(product - ratios) / abs(ratios)
            if dev > worst:
                worst, where = dev, (n, x)
    return ("closed form vs the kernel product form: worst relative deviation %s at n=%s x=%s"
            % (mp.nstr(worst, 5), where[0], mp.nstr(where[1], 10)))


if __name__ == "__main__":
    main()
