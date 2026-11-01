#!/usr/bin/env python3
"""Analysis of the independent rational counts at the second left edge a2.

The x0 side of every comparison is the independent baryrat 2.1.2 record in
results-corrected/corrected-counts.json.  The a2 side is this directory's
cells/A2_*.final.json, produced by second_edge_a2_solver.py with the same
solver, the same acceptance rule and the same retry protocol; only the left
edge and the doubled widths differ.

Statistics, per target, over the paper's fifteen Group A widths:

  count shifts   N(a2) - N(x0): valid cells, cells differing, sum, largest
                 absolute shift, parity flips, [m/k] pair changes, per cell.
  slopes         the paper's estimator: OLS of log2 N on log2 W over the
                 thirteen dyadic widths 2^-7..2^5, all rows, with the
                 paper's standard error sqrt(SSR/(n-2)/SXX).  The two edges
                 are fitted separately; the difference carries the unpaired
                 standard error and the paired one, the latter from the
                 covariance of the two edge residuals at the same widths
                 (analyze_second_edge.py's paired_slope_test).
  L              the absolute target is carried through
                 L(a) = log(F_0(a)/eps) at each edge, natural log, exactly as
                 results-second-edge/analysis.json does: 29.267707262998931
                 at x0 and 28.921134741707567 at a2.  Relative targets have
                 L = log(1/eps) at both edges.
  C_r and model  the count-implied constant C_r = exp((N-1) mu - L), with
                 mu = (pi/2) K(1-kappa^2)/K(kappa^2) and the exactly matched
                 cross-ratio kappa, and the model prediction at a2 carrying
                 the x0 constant across; residual = N(a2) - N_model.

Outputs in this directory:
  second-edge-a2-counts.json   the 45-cell record (corrected-counts schema)
  second-edge-a2-counts.csv    the same cells, one row each
  counts_by_width_target.csv   N(x0), N(a2) and the shift, both edges
  comparison.csv               every compared field
  analysis.json                the summary: count shifts, slopes, paired z,
                               L, C_r and the model residuals
"""
from __future__ import annotations

import csv
import importlib.metadata
import json
import math
import sys
import time
from pathlib import Path

import mpmath as mp

HERE = Path(__file__).resolve().parent
PAPER = HERE.parent
NAMES15 = ["m7", "m6", "m5", "m4", "m3", "m2", "m1",
           "p0", "p1", "p2", "p3", "p4", "p5", "ship", "paper"]
DYAD = ["m7", "m6", "m5", "m4", "m3", "m2", "m1",
        "p0", "p1", "p2", "p3", "p4", "p5"]
TARGETS = ["rel:1e-13", "abs:5e-14", "rel:1e-14"]
TAG_TO_TARGET = {"rel13": "rel:1e-13", "abs514": "abs:5e-14",
                 "rel14": "rel:1e-14"}
X0_TEXT = "11.899848152108484"
mp.mp.dps = 50
X0 = mp.mpf(X0_TEXT)
A2 = 2 * X0


def f0(a):
    return mp.sqrt(mp.pi) * mp.erf(mp.sqrt(a)) / (2 * mp.sqrt(a))


def target_L(target, a):
    kind, value = target.split(":", 1)
    eps = mp.mpf(value)
    if kind == "rel":
        return mp.log(1 / eps), "rel", eps
    return mp.log(f0(a) / eps), "abs", eps


def kappa(a, width):
    s = mp.sqrt(a)
    b = mp.sqrt(a + width)
    return (b - s) / (b + s)


def mu_of(kappa_value):
    return ((mp.pi / 2) * mp.ellipk(1 - kappa_value ** 2)
            / mp.ellipk(kappa_value ** 2))


def ols(xs, ys):
    n = len(xs)
    mx = sum(xs) / n
    my = sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    beta = sxy / sxx
    alpha = my - beta * mx
    resid = [y - (alpha + beta * x) for x, y in zip(xs, ys)]
    s2 = sum(r * r for r in resid) / (n - 2)
    se = math.sqrt(s2 / sxx)
    return {"alpha": alpha, "beta": beta, "se": se,
            "residuals": resid, "n": n, "sxx": sxx}


def paired_slope_test(fit1, fit2):
    """Slope difference (edge 2 minus edge 1), unpaired and paired standard
    errors, exactly analyze_second_edge.py's paired_slope_test."""
    diff = fit2["beta"] - fit1["beta"]
    se_unpaired = math.sqrt(fit1["se"] ** 2 + fit2["se"] ** 2)
    n = fit1["n"]
    cov = sum(r1 * r2 for r1, r2 in
              zip(fit1["residuals"], fit2["residuals"])) / (n - 2)
    se_paired = math.sqrt(max(0.0, fit1["se"] ** 2 + fit2["se"] ** 2
                              - 2 * cov / fit1["sxx"]))
    return {"diff": diff,
            "se_unpaired": se_unpaired,
            "se_paired": se_paired,
            "z_unpaired": diff / se_unpaired if se_unpaired else None,
            "z_paired": diff / se_paired if se_paired else None,
            "residual_cov": cov}


def load_x0():
    doc = json.loads((PAPER / "results-corrected" / "corrected-counts.json")
                     .read_text(encoding="utf-8"))
    cells = {(c["target"], c["width_name"]): c for c in doc["cells"]}
    widths = {name: float(doc["problem"]["widths"][name]) for name in NAMES15}
    return doc, cells, widths


def load_a2():
    cells = {}
    for path in sorted((HERE / "cells").glob("A2_*.final.json")):
        rec = json.loads(path.read_text(encoding="utf-8"))
        if rec.get("error"):
            raise SystemExit(f"error record {path}: {rec['error']}")
        if rec.get("corrected_count") is None:
            raise SystemExit(f"no meeting fit {path}")
        target = TAG_TO_TARGET[rec["target_tag"]]
        cells[(target, rec["width_name"])] = rec
    missing = [(t, n) for t in TARGETS for n in NAMES15
               if (t, n) not in cells]
    if missing:
        raise SystemExit(f"missing a2 cells: {missing}")
    return cells


def versions():
    out = {"python": sys.version.split()[0]}
    for mod in ("baryrat", "numpy", "scipy", "mpmath"):
        try:
            out[mod] = importlib.metadata.version(mod)
        except Exception:  # noqa: BLE001
            out[mod] = None
    return out


def main() -> int:
    x0_doc, x0, widths = load_x0()
    a2 = load_a2()

    cells_out = {}
    counts_rows = []
    comp_rows = []
    for target in TARGETS:
        L1, kind, eps = target_L(target, X0)
        L2, _, _ = target_L(target, A2)
        for name in NAMES15:
            w1 = widths[name]
            w2 = 2.0 * w1
            k1 = kappa(X0, mp.mpf(repr(w1)))
            k2 = kappa(A2, mp.mpf(repr(w2)))
            ell = mp.log(1 / k1)
            mu = mu_of(k1)
            c1 = x0[(target, name)]
            c2 = a2[(target, name)]
            n1 = c1["corrected_count"]
            n2 = c2["corrected_count"]
            d = n2 - n1
            cr1 = mp.exp((n1 - 1) * mu - L1)
            cr2 = mp.exp((n2 - 1) * mu - L2)
            model_real = 1 + (L2 + mp.log(cr1)) / mu
            entry = {
                "width": repr(w1),
                "width_a2": repr(w2),
                "kappa_x0": mp.nstr(k1, 17),
                "kappa_a2": mp.nstr(k2, 17),
                "kappa_rel_diff": mp.nstr(abs(k2 / k1 - 1), 6),
                "ell": mp.nstr(ell, 17),
                "mu": mp.nstr(mu, 17),
                "L_x0": mp.nstr(L1, 17),
                "L_a2": mp.nstr(L2, 17),
                "F0_x0": mp.nstr(f0(X0), 17),
                "F0_a2": mp.nstr(f0(A2), 17),
                "count_x0": n1,
                "count_a2": n2,
                "diff": d,
                "pair_x0": c1["pair"],
                "pair_a2": c2["pair"],
                "pair_changed": c1["pair"] != c2["pair"],
                "parity_x0": "odd" if n1 % 2 else "even",
                "parity_a2": "odd" if n2 % 2 else "even",
                "delivered_2001_x0": c1["delivered_2001"],
                "delivered_2001_a2": c2["delivered_2001"],
                "guard_2001_x0": c1["guard_2001"],
                "guard_2001_a2": c2["guard_2001"],
                "meets_4001_a2": c2["meets_4001"],
                "witness_converged_a2": c2["witness_converged"],
                "route_a2": c2["route"],
                "below_ratio_a2": c2["below_ratio"],
                "Cr_x0": mp.nstr(cr1, 17),
                "Cr_a2": mp.nstr(cr2, 17),
                "Cr_ratio": mp.nstr(cr2 / cr1, 17),
                "Cr_ratio_expected": mp.nstr(mp.exp(L1 - L2), 17),
                "rat_model_real_a2": mp.nstr(model_real, 17),
                "rat_model_ceil_a2": int(mp.ceil(model_real)),
                "rat_model_residual": mp.nstr(n2 - model_real, 17),
            }
            cells_out[f"{target}/{name}"] = entry
            counts_rows.append({
                "target": target, "width_name": name,
                "width_x0": entry["width"], "width_a2": entry["width_a2"],
                "count_x0": n1, "count_a2": n2, "diff": d,
                "pair_x0": f"[{c1['pair'][0]}/{c1['pair'][1]}]",
                "pair_a2": f"[{c2['pair'][0]}/{c2['pair'][1]}]",
            })
            comp_rows.append({"target": target, "width_name": name, **entry})

    # Per-target summaries.
    analysis = {
        "schema": "q330-second-edge-a2-analysis-1",
        "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "x0": X0_TEXT,
        "a2": mp.nstr(A2, 17),
        "versions": versions(),
        "method": {
            "solver": "baryrat 2.1.2, BRASIL",
            "problem": "F_0(x) on [a, a+W], t = 2(x-a)/W - 1, "
                       "p(t)/(1+q1 t+...+qk t^k), stored count (m+1)+k, "
                       "k <= 12",
            "second_edge": "a2 = 2*x0; every Group A width doubled, so the "
                           "cross-ratio and W/a are matched",
            "acceptance": "binary64 coefficients, exact barycentric-to-"
                          "monomial conversion at 70 digits, Horner by both "
                          "the fused and separate routes, worst error; "
                          "denominator guard min|den| >= 1e-3 max|den|; "
                          "2001-node certify grid, 4001-node check; "
                          "reference F_0 at 70 digits",
            "count": "smallest stored count with an exhibited admissible fit "
                     "(upper bound, not a minimality proof)",
            "x0_source": "results-corrected/corrected-counts.json "
                         "(the same independent solver)",
            "slopes": "OLS of log2 N on log2 W, thirteen dyadic widths, "
                      "standard error sqrt(SSR/(n-2)/SXX)",
            "paired_z": "slope difference with the covariance of the two "
                        "edges' residuals at the same widths; "
                        "se_paired^2 = se1^2 + se2^2 - 2 cov / SXX",
            "L": "L = log(F_0(a)/eps) for the absolute target, natural log, "
                 "as the paper's second-edge record does; "
                 "L = log(1/eps) for the relative targets",
        },
        "targets": {},
        "cells": cells_out,
    }

    for target in TARGETS:
        L1, _, eps = target_L(target, X0)
        L2, _, _ = target_L(target, A2)
        diffs = []
        parity_flips = 0
        pair_changes = 0
        max_abs = 0
        for name in NAMES15:
            e = cells_out[f"{target}/{name}"]
            diffs.append({"width": name, "x0": e["count_x0"],
                          "a2": e["count_a2"], "diff": e["diff"]})
            max_abs = max(max_abs, abs(e["diff"]))
            if e["parity_x0"] != e["parity_a2"]:
                parity_flips += 1
            if e["pair_changed"]:
                pair_changes += 1
        shifts = {
            "valid_cells": len(diffs),
            "differing_cells": sum(1 for d in diffs if d["diff"] != 0),
            "max_abs_diff": max_abs,
            "sum_diff": sum(d["diff"] for d in diffs),
            "parity_flips": parity_flips,
            "pair_changes": pair_changes,
            "cells": diffs,
        }

        xs = [-7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5]
        ys1 = [math.log2(cells_out[f"{target}/{n}"]["count_x0"])
               for n in DYAD]
        ys2 = [math.log2(cells_out[f"{target}/{n}"]["count_a2"])
               for n in DYAD]
        fit1 = ols(xs, ys1)
        fit2 = ols(xs, ys2)
        first, last = 0, len(xs) - 1
        slopes = {
            "common_widths": DYAD,
            "x0": {"alpha": fit1["alpha"], "beta": fit1["beta"],
                   "se": fit1["se"], "n": fit1["n"],
                   "endpoint": (ys1[last] - ys1[first]) / (xs[last] - xs[first])},
            "a2": {"alpha": fit2["alpha"], "beta": fit2["beta"],
                   "se": fit2["se"], "n": fit2["n"],
                   "endpoint": (ys2[last] - ys2[first]) / (xs[last] - xs[first])},
            "difference": paired_slope_test(fit1, fit2),
        }

        cr_x0 = [float(cells_out[f"{target}/{n}"]["Cr_x0"]) for n in DYAD]
        cr_a2 = [float(cells_out[f"{target}/{n}"]["Cr_a2"]) for n in DYAD]
        cr_ratio = [float(cells_out[f"{target}/{n}"]["Cr_ratio"])
                    for n in DYAD]
        resid = [float(cells_out[f"{target}/{n}"]["rat_model_residual"])
                 for n in NAMES15]
        cr_block = {
            "x0_min": min(cr_x0), "x0_max": max(cr_x0),
            "x0_median": sorted(cr_x0)[len(cr_x0) // 2],
            "a2_min": min(cr_a2), "a2_max": max(cr_a2),
            "a2_median": sorted(cr_a2)[len(cr_a2) // 2],
            "ratio_min": min(cr_ratio), "ratio_max": max(cr_ratio),
            "ratio_median": sorted(cr_ratio)[len(cr_ratio) // 2],
            "expected_ratio": float(mp.exp(L1 - L2)),
            "model_residual_min": min(resid),
            "model_residual_max": max(resid),
            "model_cells_above_zero": sum(1 for r in resid if r > 0),
        }
        analysis["targets"][target] = {
            "count_shifts": shifts,
            "slopes": slopes,
            "L": {"x0": mp.nstr(L1, 17), "a2": mp.nstr(L2, 17),
                  "shift": mp.nstr(L2 - L1, 17),
                  "kind": "relative" if target.startswith("rel")
                  else "absolute"},
            "Cr_and_model": cr_block,
        }

    # Assembled second-edge record, corrected-counts schema.
    a2_cells = []
    for target in TARGETS:
        for name in NAMES15:
            rec = a2[(target, name)]
            e = cells_out[f"{target}/{name}"]
            a2_cells.append({
                "target": target, "width": e["width_a2"],
                "width_name": name,
                "corrected_count": rec["corrected_count"],
                "pair": rec["pair"],
                "delivered_2001": rec["delivered_2001"],
                "delivered_4001": rec["delivered_4001"],
                "guard": min(rec["guard_2001"], rec["guard_4001"]),
                "guard_2001": rec["guard_2001"],
                "guard_4001": rec["guard_4001"],
                "best_below": rec["best_below"],
                "below_ratio": rec["below_ratio"],
                "best_below_pair": rec["best_below_pair"],
                "route": rec["route"],
                "witness_converged": rec["witness_converged"],
                "meets_4001": rec["meets_4001"],
                "count_lowered_by_tightened":
                    rec["count_lowered_by_tightened"],
                "margin_limited": rec["margin_limited"],
                "x0_independent_count": e["count_x0"],
                "witness_P": rec["witness"]["P"],
                "witness_Q": rec["witness"]["Q"],
            })
    record = {
        "schema": "q330-second-edge-a2-counts-1",
        "generated_utc": analysis["generated_utc"],
        "versions": analysis["versions"],
        "problem": {
            "family": "F_0(x) = sqrt(pi)/(2 sqrt(x)) erf(sqrt(x))",
            "x0": X0_TEXT,
            "a2": mp.nstr(A2, 17),
            "widths_x0": {k: repr(v) for k, v in widths.items()},
            "widths_a2": {k: repr(2.0 * v) for k, v in widths.items()},
            "pair": "p(t)/(1+q1 t+...+qk t^k), stored count = (m+1)+k",
            "t_grid": "2001 nodes uniform in t = 2(x-a2)/W2 - 1",
            "check_grid": "4001 nodes uniform in t",
            "k_max": 12,
            "stored_range": [3, "scan bound per cell = x0 count + 3"],
        },
        "acceptance": {
            "evaluation": "binary64 coefficients, Horner, worse of fused "
                          "(std::fma) and separate routes",
            "reference": "F_0 at 70 decimal digits (mpmath)",
            "denominator_guard": "min|den| >= max|den| * 1e-3",
            "targets": {"rel:1e-13": "max |fhat-f|/f <= 1e-13",
                        "abs:5e-14": "max |fhat-f| <= 5e-14",
                        "rel:1e-14": "max |fhat-f|/f <= 1e-14"},
            "primary_grid": 2001, "check_grid": 4001,
        },
        "solver": {
            "library": "baryrat", "algorithm": "BRASIL",
            "default_settings": {"tol": 1e-4, "maxiter": 1000, "npi": -30,
                                 "init_steps": 100, "max_step_size": 0.1},
            "tightened_settings": {"tol": 1e-8, "maxiter": 5000, "npi": -40,
                                   "init_steps": 300, "max_step_size": 0.05},
            "retry_protocol": "one level below the default first meeting "
                              "count is re-solved with the tightened "
                              "settings; if it meets, the count drops and the "
                              "next level below is retried",
        },
        "cells": a2_cells,
        "notes": [
            "A count is the smallest stored count at which an exhibited "
            "admissible [m/k] fit meets on the 2001 grid and passes the "
            "denominator guard; the 4001 grid is the check.",
            "delivered_2001/4001 are in the target's kind (relative for "
            "rel:*, absolute for abs:*).",
            "guard is the smaller of the 2001 and 4001 denominator ratios of "
            "the selected witness.",
            "best_below is the smallest delivered error at count-1 among "
            "guard-passing pairs; below_ratio = best_below/target.",
            "x0_independent_count is the same solver's count at x0 from "
            "results-corrected/corrected-counts.json.",
            "No count is a proof of minimality; it is an existence bound "
            "under this solver.",
        ],
    }

    (HERE / "analysis.json").write_text(
        json.dumps(analysis, indent=2, default=str) + "\n", encoding="utf-8")
    (HERE / "second-edge-a2-counts.json").write_text(
        json.dumps(record, indent=1, default=str) + "\n", encoding="utf-8")

    def write_csv(path, fields, rows):
        with path.open("w", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(fh, fieldnames=fields)
            writer.writeheader()
            writer.writerows(rows)

    write_csv(HERE / "second-edge-a2-counts.csv",
              ["target", "width", "width_name", "corrected_count", "pair",
               "delivered_2001", "delivered_4001", "guard", "below_ratio",
               "route", "x0_independent_count"],
              [{"target": c["target"], "width": c["width"],
                "width_name": c["width_name"],
                "corrected_count": c["corrected_count"],
                "pair": f"[{c['pair'][0]}/{c['pair'][1]}]",
                "delivered_2001": repr(c["delivered_2001"]),
                "delivered_4001": repr(c["delivered_4001"]),
                "guard": repr(c["guard"]), "below_ratio": repr(c["below_ratio"]),
                "route": c["route"],
                "x0_independent_count": c["x0_independent_count"]}
               for c in a2_cells])
    write_csv(HERE / "counts_by_width_target.csv",
              ["target", "width_name", "width_x0", "width_a2", "count_x0",
               "count_a2", "diff", "pair_x0", "pair_a2"], counts_rows)
    comp_fields = ["target", "width_name", "width", "width_a2", "kappa_x0",
                   "kappa_a2", "kappa_rel_diff", "ell", "mu", "L_x0", "L_a2",
                   "F0_x0", "F0_a2", "count_x0", "count_a2", "diff",
                   "pair_x0", "pair_a2", "pair_changed", "parity_x0",
                   "parity_a2", "delivered_2001_x0", "delivered_2001_a2",
                   "guard_2001_x0", "guard_2001_a2", "meets_4001_a2",
                   "witness_converged_a2", "route_a2", "below_ratio_a2",
                   "Cr_x0", "Cr_a2", "Cr_ratio", "Cr_ratio_expected",
                   "rat_model_real_a2", "rat_model_ceil_a2",
                   "rat_model_residual"]
    write_csv(HERE / "comparison.csv", comp_fields, comp_rows)

    # Terminal report.
    print(f"x0 = {X0_TEXT}; a2 = {mp.nstr(A2, 17)}")
    for target in TARGETS:
        block = analysis["targets"][target]
        s = block["slopes"]
        d = s["difference"]
        sh = block["count_shifts"]
        print(f"{target:>10}: x0 {s['x0']['beta']:.3f}+-{s['x0']['se']:.3f}"
              f" (n={s['x0']['n']}), a2 {s['a2']['beta']:.3f}"
              f"+-{s['a2']['se']:.3f}, diff {d['diff']:+.3f} "
              f"(unpaired se {d['se_unpaired']:.3f}, z {d['z_unpaired']:+.2f};"
              f" paired se {d['se_paired']:.4f}, z {d['z_paired']:+.2f}); "
              f"counts differ {sh['differing_cells']}/{sh['valid_cells']}, "
              f"sum {sh['sum_diff']:+d}, max |d| {sh['max_abs_diff']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
