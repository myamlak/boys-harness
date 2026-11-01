#!/usr/bin/env python3
"""Analyse the F_1 second-instance counts against Paper C's model.

Inputs
  results-second-instance-f1/sweep_all.json
        the F_1 measurement at both Group A edges (x0 and a2 = 2 x0), written
        by second_instance_f1_sweep.py with the paper's own harness.
  the paper's own F_0 records (x0: the published re-count plus the core
        sweep; a2: results-second-edge/sweep_all.json)
        read only to print the F_0 bands beside the F_1 bands.

The estimators are the paper's, copied from analyze_second_edge.py:
  * kappa(W/a) = (sqrt(a+W)-sqrt(a))/(sqrt(a+W)+sqrt(a)), ell = log(1/kappa);
  * the Stahl/Groetzsch modulus mu(kappa) = (pi/2) K(1-kappa^2)/K(kappa^2);
  * the count laws N = 1 + (L + log C)/rate, L = log(1/eps) for a relative
    target and log(f(a)/eps) for an absolute one;
  * the slope estimator: least squares of log2 N on log2 W over the thirteen
    dyadic widths 2^-7..2^5, floors included, with the paper's paired and
    unpaired slope-difference errors.

For F_1 the tail is alpha = 3/2, F_1 = (sqrt(pi)/4) x^(-3/2) g(x), g -> 1.
The tail-class proposition gives the Chebyshev rate ell for every alpha.  The
rational rate for alpha = 3/2 is the Stahl modulus of the corresponding
sheet problem.  Because x^(-3/2) = x^(-1/2) * (1/x) and 1/x is analytic and
non-vanishing on [a, b] with a > 0, analytic-factor invariance (the sandwich
argument) gives the same exponential rate as for alpha = 1/2:
mu_{3/2}(kappa) = mu(kappa), the Groetzsch ring function, with the
alpha-dependence confined to the constant C_r.  The sheet count of x^{3/2} is
2 (the denominator of 3/2 in lowest terms), not 3.  No direct numerical
solution of an alpha-sheeted equilibrium was performed; the rational-rate
prediction used here is that derived rate, and the measurement is reported
against it, and also on its own in case a reader declines the derivation.

Writes, in results-second-instance-f1/:
  analysis.json             the summary, slopes, bands and model residuals
  counts_by_width_target.csv   every measured count, both edges
  comparison.csv            the matched-W/a cell comparison
  model.csv                 per-cell implied constants and model residuals
  f0_bands.json             the same comparison computed for F_0
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path

import mpmath as mp

PAPER = Path(__file__).resolve().parent
sys.path.insert(0, str(PAPER))

NAMES = [f"2^{j}" for j in range(-7, 6)] + ["ship", "paper"]
DYAD = [f"2^{j}" for j in range(-7, 6)]
TARGETS = ["rel:1e-13", "abs:5e-14", "rel:1e-14"]
X0_TEXT = "11.899848152108484"
X1_TEXT = "28.989337738820740"

DEFAULT_F1 = PAPER / "results-second-instance-f1" / "sweep_all.json"
DEFAULT_F0_RECOUNT = Path(
    r"C:\Users\myaml\dev\papers-harness-repo\coefficient-count-sweep"
    r"\records\chebyshev_recount.json")
DEFAULT_F0_CORE = Path(
    r"C:\Users\myaml\dev\papers-harness-repo\coefficient-count-sweep"
    r"\records\core")
DEFAULT_F0_A2 = PAPER / "results-second-edge" / "sweep_all.json"
DEFAULT_BOYS = Path(r"C:\Users\myaml\dev\_mainonly")

mp.mp.dps = 50
X0 = mp.mpf(X0_TEXT)


# ---------------------------------------------------------------------------
# the paper's estimators


def kappa(a, width):
    s = mp.sqrt(a)
    b = mp.sqrt(a + width)
    return (b - s) / (b + s)


def mu_of(kappa_value):
    return (mp.pi / 2) * mp.ellipk(1 - kappa_value ** 2) / mp.ellipk(kappa_value ** 2)


def width_of(name: str):
    if name.startswith("2^"):
        return mp.mpf(2) ** int(name.split("^")[1])
    if name == "ship":
        return mp.mpf("31") - X0
    if name == "paper":
        return mp.mpf(X1_TEXT) - X0
    raise ValueError(name)


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
            "residuals": resid, "n": n, "sxx": sxx,
            "resid_ss": sum(r * r for r in resid)}


def paired_slope_test(fit1, fit2):
    diff = fit2["beta"] - fit1["beta"]
    se_unpaired = math.sqrt(fit1["se"] ** 2 + fit2["se"] ** 2)
    n = fit1["n"]
    cov = sum(r1 * r2 for r1, r2 in zip(fit1["residuals"], fit2["residuals"])) / (n - 2)
    se_paired = math.sqrt(max(0.0, fit1["se"] ** 2 + fit2["se"] ** 2
                              - 2 * cov / fit1["sxx"]))
    return {"diff": diff, "se_unpaired": se_unpaired, "se_paired": se_paired,
            "z_unpaired": diff / se_unpaired if se_unpaired else None,
            "z_paired": diff / se_paired if se_paired else None,
            "residual_cov": cov}


def fmt(value, digits=6):
    if value is None:
        return "--"
    if isinstance(value, (int, str)):
        return str(value)
    return f"{value:.{digits}g}"


# ---------------------------------------------------------------------------
# loading


def load_f1(path: Path):
    if not path.exists():
        raise SystemExit(f"F_1 sweep record missing: {path}")
    doc = json.loads(path.read_text(encoding="utf-8"))
    rows = {}
    for row in doc["rows"]:
        rows[(row["edge"], row["target"], row["width_name"])] = row
    return doc, rows


def f1_edge_maps(rows, edge):
    out = {"cheb": {}, "rat": {}}
    for (e, target, name), row in rows.items():
        if e != edge or name not in NAMES:
            continue
        for fam in ("cheb", "rat"):
            entry = row[fam]
            out[fam][(target, name)] = {
                "stored": entry.get("stored"),
                "m": entry.get("m"),
                "k": entry.get("k"),
                "solver": entry.get("solver"),
                "delivered_rel": entry.get("delivered_rel"),
                "delivered_abs": entry.get("delivered_abs"),
                "below_meets": entry.get("below_meets"),
                "status": entry.get("status"),
            }
    return out


def load_f0_x0(recount_path: Path, core_dir: Path):
    doc = json.loads(recount_path.read_text(encoding="utf-8"))
    cheb = {}
    for cell in doc["cells"]:
        if cell.get("mode") != "records" or cell.get("nodes") != 2000:
            continue
        if cell["target"] not in TARGETS:
            continue
        names = [n for n in cell.get("width_names", []) if n in NAMES]
        if not names:
            continue
        cheb[(cell["target"], names[0])] = {
            "stored": cell["new_stored"], "m": cell["new_m"], "k": None,
            "solver": None, "delivered_rel": cell["delivered_rel"],
            "delivered_abs": cell["delivered_abs"], "below_meets": None,
        }
    rat = {}
    for fname in ("sweep_all.json", "sweep_rel_1e-14.json"):
        path = core_dir / fname
        if not path.exists():
            continue
        part = json.loads(path.read_text(encoding="utf-8"))
        for row in part["rows"]:
            if row["width_name"] not in NAMES or row["target"] not in TARGETS:
                continue
            r = row["rat"]
            rat[(row["target"], row["width_name"])] = {
                "stored": r.get("stored"), "m": r.get("m"), "k": r.get("k"),
                "solver": r.get("solver"),
                "delivered_rel": r.get("delivered_rel"),
                "delivered_abs": r.get("delivered_abs"),
                "below_meets": r.get("below_meets"),
            }
    return {"cheb": cheb, "rat": rat}


def load_f0_a2(path: Path):
    doc = json.loads(path.read_text(encoding="utf-8"))
    out = {"cheb": {}, "rat": {}}
    for row in doc["rows"]:
        if row["width_name"] not in NAMES or row["target"] not in TARGETS:
            continue
        for fam in ("cheb", "rat"):
            entry = row[fam]
            out[fam][(row["target"], row["width_name"])] = {
                "stored": entry.get("stored"), "m": entry.get("m"),
                "k": entry.get("k"), "solver": entry.get("solver"),
                "delivered_rel": entry.get("delivered_rel"),
                "delivered_abs": entry.get("delivered_abs"),
                "below_meets": entry.get("below_meets"),
            }
    return out


def width_table(rows):
    """The Group A widths as the sweep stored them, with a formula fallback.

    The measurement's own width text is the authority on kappa; width_of is
    the fallback for a cell the sweep has not reached yet.
    """
    out = {}
    for (edge, _target, name), row in rows.items():
        if name in NAMES and edge == "x0" and row.get("width"):
            out.setdefault(name, mp.mpf(row["width"]))
    for (_edge, _target, name), row in rows.items():
        if name in NAMES and name not in out and row.get("width_x0"):
            out.setdefault(name, mp.mpf(row["width_x0"]))
    for name in NAMES:
        out.setdefault(name, width_of(name))
    return out


def f1_at(g, a):
    """F_1(a) at 60 digits from the paper's reference routine."""
    import f1_reference as f1r
    with g.mp.workdps(g.RAT_DPS):
        value = f1r.f1_from_f0(g, g.mpf(a))
    return mp.mpf(g.mp.nstr(value, 45))


def target_L(target, a, fval):
    kind, value = target.split(":", 1)
    eps = mp.mpf(value)
    if kind == "rel":
        return mp.log(1 / eps)
    return mp.log(fval / eps)


# ---------------------------------------------------------------------------
# comparison and model machinery


def require(maps, label):
    missing = []
    for fam in ("cheb", "rat"):
        for target in TARGETS:
            for name in NAMES:
                key = (target, name)
                entry = maps[fam].get(key)
                if entry is None or entry.get("stored") is None:
                    missing.append(f"{fam}/{target}/{name}")
    return missing


def compare(m1, m2, L1, L2, widths, label):
    """The paper's matched-W/a comparison for one function.

    L1/L2 map target -> model L at the two edges; ``widths`` maps width name
    to the Group A width (the a2 edge's width is twice it).  Returns cells
    keyed target/width and per-target/family summaries.  The valid-cell rule
    is the paper's: a rational cell is dropped when the below-count check met
    (below_meets true at either edge).
    """
    cells = {}
    for target in TARGETS:
        for name in NAMES:
            w1 = widths[name]
            w2 = w1 * 2
            k1 = kappa(X0, w1)
            k2 = kappa(X0 * 2, w2)
            ell = mp.log(1 / k1)
            mu = mu_of(k1)
            entry = {
                "width": mp.nstr(w1, 17),
                "width_a2": mp.nstr(w2, 17),
                "kappa_x0": mp.nstr(k1, 17),
                "kappa_a2": mp.nstr(k2, 17),
                "kappa_rel_diff": mp.nstr(abs(k2 / k1 - 1), 6),
                "ell": mp.nstr(ell, 17),
                "mu": mp.nstr(mu, 17),
                "L_x0": mp.nstr(L1[target], 17),
                "L_a2": mp.nstr(L2[target], 17),
            }
            for fam, const in (("cheb", "Cc"), ("rat", "Cr")):
                e1 = m1[fam].get((target, name), {})
                e2 = m2[fam].get((target, name), {})
                n1, n2 = e1.get("stored"), e2.get("stored")
                entry[f"{fam}_x0"] = n1
                entry[f"{fam}_a2"] = n2
                entry[f"{fam}_diff"] = (None if n1 is None or n2 is None
                                        else n2 - n1)
                entry[f"{fam}_m_x0"] = e1.get("m")
                entry[f"{fam}_m_a2"] = e2.get("m")
                entry[f"{fam}_k_x0"] = e1.get("k")
                entry[f"{fam}_k_a2"] = e2.get("k")
                entry[f"{fam}_solver_x0"] = e1.get("solver")
                entry[f"{fam}_solver_a2"] = e2.get("solver")
                entry[f"{fam}_deliv_rel_x0"] = e1.get("delivered_rel")
                entry[f"{fam}_deliv_rel_a2"] = e2.get("delivered_rel")
                entry[f"{fam}_deliv_abs_x0"] = e1.get("delivered_abs")
                entry[f"{fam}_deliv_abs_a2"] = e2.get("delivered_abs")
                entry[f"{fam}_below_x0"] = e1.get("below_meets")
                entry[f"{fam}_below_a2"] = e2.get("below_meets")
                entry[f"{fam}_parity_x0"] = (None if n1 is None
                                             else ("odd" if n1 % 2 else "even"))
                entry[f"{fam}_parity_a2"] = (None if n2 is None
                                             else ("odd" if n2 % 2 else "even"))
                rate = ell if fam == "cheb" else mu
                if n1 is not None:
                    entry[f"{const}_x0"] = mp.nstr(
                        mp.exp((n1 - 1) * rate - L1[target]), 17)
                if n2 is not None:
                    entry[f"{const}_a2"] = mp.nstr(
                        mp.exp((n2 - 1) * rate - L2[target]), 17)
                if n1 is not None and n2 is not None:
                    entry[f"{const}_ratio"] = mp.nstr(
                        mp.mpf(entry[f"{const}_a2"]) / mp.mpf(entry[f"{const}_x0"]), 17)
                    entry[f"{const}_ratio_expected"] = mp.nstr(
                        mp.exp(L1[target] - L2[target]), 17)
                # the model's prediction at a2 from the x0 count's constant
                if n1 is not None:
                    c1 = mp.exp((n1 - 1) * rate - L1[target])
                    real = 1 + (L2[target] + mp.log(c1)) / rate
                    entry[f"{fam}_model_real_a2"] = mp.nstr(real, 17)
                    entry[f"{fam}_model_ceil_a2"] = int(mp.ceil(real))
                    entry[f"{fam}_model_resid_a2"] = n2 - real if n2 is not None else None
            cells[f"{target}/{name}"] = entry

    targets = {}
    for target in TARGETS:
        block = {"cheb": {}, "rat": {}}
        for fam in ("cheb", "rat"):
            diffs = []
            parity_flips = 0
            pair_changes = 0
            max_abs = 0
            for name in NAMES:
                e = cells[f"{target}/{name}"]
                n1, n2 = e[f"{fam}_x0"], e[f"{fam}_a2"]
                if n1 is None or n2 is None:
                    continue
                if fam == "rat" and (e["rat_below_x0"] or e["rat_below_a2"]):
                    continue
                d = n2 - n1
                diffs.append({"width": name, "x0": n1, "a2": n2, "diff": d})
                max_abs = max(max_abs, abs(d))
                if e[f"{fam}_parity_x0"] != e[f"{fam}_parity_a2"]:
                    parity_flips += 1
                if (e[f"{fam}_m_x0"], e.get(f"{fam}_k_x0")) != (
                        e[f"{fam}_m_a2"], e.get(f"{fam}_k_a2")):
                    pair_changes += 1
            block[fam] = {
                "valid_cells": len(diffs),
                "differing_cells": sum(1 for d in diffs if d["diff"] != 0),
                "max_abs_diff": max_abs,
                "sum_diff": sum(d["diff"] for d in diffs),
                "parity_flips": parity_flips,
                "pair_changes": pair_changes,
                "cells": diffs,
            }
        targets[target] = block
    return cells, targets


def slope_block(m1, m2, fam, target, rate_of):
    js = []
    ys1, ys2 = [], []
    for j in range(-7, 6):
        name = f"2^{j}"
        e1 = m1[fam].get((target, name), {})
        e2 = m2[fam].get((target, name), {})
        n1, n2 = e1.get("stored"), e2.get("stored")
        if n1 is None or n2 is None:
            continue
        if fam == "rat" and (e1.get("below_meets") or e2.get("below_meets")):
            continue
        js.append(j)
        ys1.append(math.log2(n1))
        ys2.append(math.log2(n2))
    block = {"common_widths": [f"2^{j}" for j in js], "n": len(js)}
    if len(js) < 10:
        block["status"] = "incomplete"
        return block
    xs = [float(j) for j in js]
    fit1 = ols(xs, ys1)
    fit2 = ols(xs, ys2)
    xs_span = xs[-1] - xs[0]
    block.update({
        "x0": {"alpha": fit1["alpha"], "beta": fit1["beta"],
               "se": fit1["se"], "n": fit1["n"],
               "endpoint": (ys1[-1] - ys1[0]) / xs_span},
        "a2": {"alpha": fit2["alpha"], "beta": fit2["beta"],
               "se": fit2["se"], "n": fit2["n"],
               "endpoint": (ys2[-1] - ys2[0]) / xs_span},
        "difference": paired_slope_test(fit1, fit2),
    })
    return block


def model_slope(fam, target, a, widths, L_edge, const=1.0):
    """The model's own width-law slope at a chosen constant C."""
    xs, ys = [], []
    for j in range(-7, 6):
        name = f"2^{j}"
        w = widths[name]
        k1 = kappa(a, w)
        rate = mp.log(1 / k1) if fam == "cheb" else mu_of(k1)
        n = 1 + (L_edge[target] + mp.log(mp.mpf(const))) / rate
        if n <= 0:
            continue
        xs.append(float(j))
        ys.append(float(mp.log(n) / mp.log(2)))
    fit = ols(xs, ys)
    return {"alpha": fit["alpha"], "beta": fit["beta"], "n": fit["n"]}


def constants_band(m1, fam, target, a, widths, L_edge):
    """The implied constant over the dyadic widths at one edge."""
    vals = []
    for name in DYAD:
        entry = m1[fam].get((target, name), {})
        n = entry.get("stored")
        if n is None:
            continue
        if fam == "rat" and entry.get("below_meets"):
            continue
        w = widths[name]
        k1 = kappa(a, w)
        rate = (mp.log(1 / k1) if fam == "cheb" else mu_of(k1))
        vals.append(float(mp.exp((n - 1) * rate - L_edge[target])))
    if not vals:
        return None
    vals_sorted = sorted(vals)
    return {"min": min(vals), "max": max(vals),
            "median": vals_sorted[len(vals_sorted) // 2], "n": len(vals),
            "values": vals}


def width_law_residuals(m1, fam, target, a, widths, L_edge):
    """Residuals of N against N = 1 + (L + log C_med)/rate, dyadic widths."""
    band = constants_band(m1, fam, target, a, widths, L_edge)
    if band is None:
        return None
    c_med = mp.mpf(band["median"])
    out = []
    for name in DYAD:
        entry = m1[fam].get((target, name), {})
        n = entry.get("stored")
        if n is None:
            continue
        if fam == "rat" and entry.get("below_meets"):
            continue
        w = widths[name]
        k1 = kappa(a, w)
        rate = (mp.log(1 / k1) if fam == "cheb" else mu_of(k1))
        real = 1 + (L_edge[target] + mp.log(c_med)) / rate
        out.append({"width": name, "n": n, "real": float(real),
                    "residual": n - float(real)})
    res = [abs(o["residual"]) for o in out]
    return {"C_median": float(c_med), "cells": out,
            "max_abs_residual": max(res) if res else None,
            "n": len(out)}


def target_increments(m1, fam, target_pairs, a, widths, L_edge, floor):
    """Measured vs model count increment between two targets at one width."""
    out = []
    for name in NAMES:
        for lo, hi in target_pairs:
            elo = m1[fam].get((lo, name), {})
            ehi = m1[fam].get((hi, name), {})
            nlo, nhi = elo.get("stored"), ehi.get("stored")
            if nlo is None or nhi is None or nlo <= floor or nhi <= floor:
                continue
            if fam == "rat" and (elo.get("below_meets")
                                 or ehi.get("below_meets")):
                continue
            w = widths[name]
            k1 = kappa(a, w)
            rate = (mp.log(1 / k1) if fam == "cheb" else mu_of(k1))
            pred = float((L_edge[hi] - L_edge[lo]) / rate)
            out.append({"width": name, "from": lo, "to": hi,
                        "n_from": nlo, "n_to": nhi,
                        "measured": nhi - nlo, "predicted": pred,
                        "residual": (nhi - nlo) - pred})
    res = [abs(o["residual"]) for o in out]
    return {"cells": out, "n": len(out),
            "max_abs_residual": max(res) if res else None,
            "mean_residual": (sum(o["residual"] for o in out) / len(out))
            if out else None,
            "outside_one": sum(1 for r in res if r > 1.0),
            "outside_two": sum(1 for r in res if r > 2.0),
            "signs": ["+" if o["residual"] > 0
                      else ("-" if o["residual"] < 0 else "0")
                      for o in out]}


def ladder_test(m, targets, L_edge, floors=None):
    """The paper's recorded free-line test, applied to the three measured arms.

    At each width and family, above the solver floor, fit the stored count
    against L with a free line and require the largest residual and every
    second difference to be at most one stored coefficient; no model rate
    enters.  The three Group A targets are not equally spaced in L, so the
    second difference is corrected by the line's own spacing term,
    ``slope*(L2 - 2 L1 + L0)``; with equal spacing the correction vanishes
    and the raw second difference is the paper's.  A width with fewer than
    three above-floor arms has no free line and is reported as untestable,
    not as a pass.
    """
    floors = {"cheb": 5, "rat": 3} if floors is None else floors
    out = {}
    for fam in ("cheb", "rat"):
        rows_out = []
        failed = []
        for name in NAMES:
            pts = []
            for target in targets:
                e = m[fam].get((target, name), {})
                n = e.get("stored")
                if n is None:
                    continue
                if fam == "rat" and e.get("below_meets"):
                    continue
                pts.append((float(L_edge[target]), n))
            entry = {"width": name, "points_total": len(pts),
                     "counts": [n for _, n in pts]}
            used = sorted((L, n) for L, n in pts if n > floors[fam])
            entry["points_above_floor"] = len(used)
            if len(used) < 3:
                entry["status"] = "fewer than three above-floor targets"
                rows_out.append(entry)
                continue
            Ls = [L for L, _ in used]
            ns = [float(n) for _, n in used]
            fit = ols(Ls, ns)
            resid = [n - (fit["alpha"] + fit["beta"] * L)
                     for L, n in used]
            raw2 = [ns[i + 1] - 2 * ns[i] + ns[i - 1]
                    for i in range(1, len(ns) - 1)]
            spacing = [fit["beta"] * (Ls[i + 1] - 2 * Ls[i] + Ls[i - 1])
                       for i in range(1, len(Ls) - 1)]
            corr2 = [r - s for r, s in zip(raw2, spacing)]
            entry.update({
                "slope": fit["beta"], "slope_se": fit["se"],
                "intercept": fit["alpha"],
                "L": [mp.nstr(mp.mpf(L), 17) for L in Ls],
                "residuals": resid,
                "second_differences": raw2,
                "second_difference_spacing_terms": spacing,
                "second_differences_corrected": corr2,
                "max_abs_residual": max(abs(r) for r in resid),
                "max_abs_second_difference_corrected":
                    max(abs(s) for s in corr2),
            })
            entry["pass"] = (entry["max_abs_residual"] <= 1.0
                             and entry["max_abs_second_difference_corrected"]
                             <= 1.0)
            if not entry["pass"]:
                failed.append(name)
            rows_out.append(entry)
        tested = [r for r in rows_out if "pass" in r]
        out[fam] = {
            "rows": rows_out,
            "tested_widths": [r["width"] for r in tested],
            "n_tested": len(tested),
            "failed": failed,
            "max_abs_residual": max(
                (r["max_abs_residual"] for r in tested), default=None),
            "max_abs_second_difference": max(
                (r["max_abs_second_difference_corrected"] for r in tested),
                default=None),
            "verdict": ("pass" if tested and not failed
                        else ("fail" if failed else "not testable")),
        }
    return out


def width_affine(m, fam, target, a, widths, floor):
    """Free line N-1 = p + q/rate over the dyadic widths, floors excluded.

    The model says q = L + log C with C constant; the test is whether one
    line fits the thirteen widths to within the counting quantum.  Integer
    rounding is the only noise, so a residual over one coefficient is a real
    departure, not scatter.
    """
    xs, ys, names = [], [], []
    for name in DYAD:
        e = m[fam].get((target, name), {})
        n = e.get("stored")
        if n is None or n <= floor:
            continue
        if fam == "rat" and e.get("below_meets"):
            continue
        w = widths[name]
        k1 = kappa(a, w)
        rate = (mp.log(1 / k1) if fam == "cheb" else mu_of(k1))
        xs.append(float(1 / rate))
        ys.append(n - 1)
        names.append(name)
    if len(xs) < 3:
        return {"status": "too few rows above the floor", "n": len(xs),
                "widths": names}
    fit = ols(xs, ys)
    resid = [y - (fit["alpha"] + fit["beta"] * x) for x, y in zip(xs, ys)]
    return {"n": len(xs), "widths": names,
            "q": fit["beta"], "q_se": fit["se"],
            "intercept": fit["alpha"],
            "residuals": [{"width": nm, "residual": r}
                          for nm, r in zip(names, resid)],
            "max_abs_residual": max(abs(r) for r in resid)}


def matched_model_residuals(cells, target, fam):
    """The a2 count against the model carried over from the matched x0 cell."""
    vals = []
    for name in NAMES:
        e = cells[f"{target}/{name}"]
        r = e.get(f"{fam}_model_resid_a2")
        if r is None:
            continue
        if fam == "rat" and (e["rat_below_x0"] or e["rat_below_a2"]):
            continue
        vals.append({"width": name, "residual": r})
    if not vals:
        return None
    mags = [abs(v["residual"]) for v in vals]
    return {"n": len(vals), "cells": vals,
            "max_abs_residual": max(mags),
            "mean_residual": sum(v["residual"] for v in vals) / len(vals),
            "within_one": sum(1 for m in mags if m <= 1.0),
            "within_two": sum(1 for m in mags if m <= 2.0)}


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sweep", default=str(DEFAULT_F1))
    parser.add_argument("--f0-recount", default=str(DEFAULT_F0_RECOUNT))
    parser.add_argument("--f0-core", default=str(DEFAULT_F0_CORE))
    parser.add_argument("--f0-a2", default=str(DEFAULT_F0_A2))
    parser.add_argument("--boys-root", default=str(DEFAULT_BOYS))
    parser.add_argument("--out-dir",
                        default=str(PAPER / "results-second-instance-f1"))
    args = parser.parse_args(argv)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    import count_sweep as cs
    import f1_reference as f1r

    doc, rows = load_f1(Path(args.sweep))
    meta = doc.get("meta", {})
    boys = Path(args.boys_root)
    g = cs.load_generator(boys)
    A2 = X0 * 2
    f1_x0 = f1_at(g, X0)
    f1_a2 = f1_at(g, A2)
    L1 = {t: target_L(t, X0, f1_x0) for t in TARGETS}
    L2 = {t: target_L(t, A2, f1_a2) for t in TARGETS}

    maps1 = f1_edge_maps(rows, "x0")
    maps2 = f1_edge_maps(rows, "a2")
    missing1 = require(maps1, "x0")
    missing2 = require(maps2, "a2")
    widths = width_table(rows)
    widths_a2 = {name: w * 2 for name, w in widths.items()}

    cells, per_target = compare(maps1, maps2, L1, L2, widths, "F1")

    summary = {
        "function": "F1",
        "tail_alpha": "3/2",
        "identity": "F_1(x) = (F_0(x) - exp(-x))/(2x)",
        "x0": X0_TEXT,
        "a2": mp.nstr(A2, 17),
        "F1_at_x0": mp.nstr(f1_x0, 25),
        "F1_at_a2": mp.nstr(f1_a2, 25),
        "targets": {t: {"L_x0": mp.nstr(L1[t], 17), "L_a2": mp.nstr(L2[t], 17)}
                    for t in TARGETS},
        "target_semantics": {
            "rel": "relative to F_1(x) at each node",
            "abs": ("plain absolute 5e-14 on F_1 (not F_0-scaled); its model "
                    "L is log(F_1(a)/5e-14)"),
        },
        "coverage": {"x0_missing": missing1, "a2_missing": missing2,
                     "rows_on_disk": len(rows)},
        "cells": cells,
        "per_target": per_target,
        "slopes": {},
        "model": {},
        "C_bands": {},
        "f0_bands": {},
        "chebyshev_rate_model": {
            "statement": (
                "the Chebyshev rate is log(1/kappa) for every alpha (the "
                "paper's tail-class proposition): the branch point at zero "
                "and the band fix kappa, and a non-vanishing analytic factor "
                "moves only the constant. F_1/F_0 = g/(2 x erf(sqrt(x))) is "
                "analytic and positive on [a, a+W] for a > 0, so F_1 has "
                "F_0's Chebyshev rate at matched W/a; the free constant is "
                "C_c."),
            "used_rate": "ell = log(1/kappa(W/a))",
        },
        "rational_rate_model": {
            "statement": (
                "F_1 = (sqrt(pi)/4) x^(-3/2) g(x) with g analytic and "
                "non-vanishing on the band; x^(-3/2) = x^(-1/2) * (1/x) and "
                "1/x is analytic and non-vanishing for a > 0, so the "
                "exponential rational rate is the alpha = 1/2 rate: "
                "mu_{3/2}(kappa) = mu(kappa), the Groetzsch ring function. "
                "The alpha-dependence is in the constant C_r, not the rate; "
                "the sheet count of x^(3/2) is 2, not 3."),
            "alpha": "3/2",
            "p_over_q": "3/2",
            "sheets": 2,
            "subsidiary_argument": (
                "F_1/F_0 is analytic and non-vanishing on the band; "
                "multiplying by such a factor leaves the nth-root rate "
                "unchanged, so both families inherit F_0's rates at matched "
                "W/a."),
            "direct_sheeted_equilibrium_solved": False,
            "fallback_if_derivation_declined": (
                "the rational-rate prediction is not established "
                "independently; the rational counts and their slopes stand "
                "as measurements"),
            "used_rate": "mu(kappa) = (pi/2) K(1-kappa^2)/K(kappa^2)",
        },
    }

    # ---- slopes, both families, both edges, all targets
    kappa_diffs = [float(e["kappa_rel_diff"]) for e in cells.values()]
    summary["matched_kappa"] = {
        "max_rel_diff": max(kappa_diffs) if kappa_diffs else None,
        "statement": ("the a2 edge doubles every Group A width, so W/a and "
                      "kappa(W/a) are matched; the sweep's stored width "
                      "text is what fixes kappa"),
    }
    for target in TARGETS:
        summary["slopes"][target] = {}
        for fam in ("cheb", "rat"):
            block = slope_block(maps1, maps2, fam, target, None)
            b1 = constants_band(maps1, fam, target, X0, widths, L1)
            b2 = constants_band(maps2, fam, target, X0 * 2, widths_a2, L2)
            block["model_x0"] = {
                "C1": model_slope(fam, target, X0, widths, L1, 1.0),
                "C_med": model_slope(fam, target, X0, widths, L1,
                                     b1["median"]) if b1 else None,
            }
            block["model_a2"] = {
                "C1": model_slope(fam, target, X0 * 2, widths_a2, L2, 1.0),
                "C_med": model_slope(fam, target, X0 * 2, widths_a2, L2,
                                     b2["median"]) if b2 else None,
            }
            summary["slopes"][target][fam] = block

    # ---- count-implied constants at each edge
    for target in TARGETS:
        summary["C_bands"][target] = {}
        for fam in ("cheb", "rat"):
            b1 = constants_band(maps1, fam, target, X0, widths, L1)
            b2 = constants_band(maps2, fam, target, X0 * 2, widths_a2, L2)
            summary["C_bands"][target][fam] = {"x0": b1, "a2": b2}

    # The paper's C_c block: the matched-cell ratio of count-implied
    # constants over the dyadic widths, against its no-drift value.
    summary["C_ratios"] = {}
    for target in TARGETS:
        band = {}
        for fam, const in (("cheb", "Cc"), ("rat", "Cr")):
            vals_x0, vals_a2, ratios = [], [], []
            for name in DYAD:
                e = cells[f"{target}/{name}"]
                if e.get(f"{const}_x0") and e.get(f"{const}_a2"):
                    vals_x0.append(float(e[f"{const}_x0"]))
                    vals_a2.append(float(e[f"{const}_a2"]))
                    ratios.append(float(e[f"{const}_ratio"]))
            if vals_x0:
                band[fam] = {
                    "x0_min": min(vals_x0), "x0_max": max(vals_x0),
                    "x0_median": sorted(vals_x0)[len(vals_x0) // 2],
                    "a2_min": min(vals_a2), "a2_max": max(vals_a2),
                    "a2_median": sorted(vals_a2)[len(vals_a2) // 2],
                    "ratio_min": min(ratios), "ratio_max": max(ratios),
                    "ratio_median": sorted(ratios)[len(ratios) // 2],
                    "expected_ratio": float(mp.exp(L1[target] - L2[target])),
                    "n": len(ratios),
                }
        summary["C_ratios"][target] = band

    # ---- width-law residuals and the free affine test
    for target in TARGETS:
        summary["model"][target] = {}
        for fam in ("cheb", "rat"):
            floor = 5 if fam == "cheb" else 3
            summary["model"][target][fam] = {
                "width_law_x0": width_law_residuals(
                    maps1, fam, target, X0, widths, L1),
                "width_law_a2": width_law_residuals(
                    maps2, fam, target, X0 * 2, widths_a2, L2),
                "width_affine_x0": width_affine(
                    maps1, fam, target, X0, widths, floor),
                "width_affine_a2": width_affine(
                    maps2, fam, target, X0 * 2, widths_a2, floor),
                "matched_model_a2": matched_model_residuals(cells, target, fam),
            }

    # ---- target-affine increments, floor rows excluded
    pairs = [("abs:5e-14", "rel:1e-13"), ("rel:1e-13", "rel:1e-14")]
    summary["target_affine"] = {}
    for fam, floor in (("cheb", 5), ("rat", 3)):
        summary["target_affine"][fam] = {
            "x0": target_increments(maps1, fam, pairs, X0, widths, L1, floor),
            "a2": target_increments(maps2, fam, pairs, X0 * 2, widths_a2, L2,
                                    floor),
        }

    # ---- the paper's recorded free-line test on the three measured arms
    summary["ladder"] = {
        "x0": ladder_test(maps1, TARGETS, L1),
        "a2": ladder_test(maps2, TARGETS, L2),
        "note": ("the paper's free-line rule (largest residual and every "
                 "second difference at most one stored coefficient, no model "
                 "rate) applied to the three Group A targets; second "
                 "differences carry the unequal-L spacing correction"),
    }

    # ---- the F_0 bands, computed by the same code on the paper's records
    try:
        f0x0 = load_f0_x0(Path(args.f0_recount), Path(args.f0_core))
        f0a2 = load_f0_a2(Path(args.f0_a2))
        f0_missing1 = require(f0x0, "F0 x0")
        f0_missing2 = require(f0a2, "F0 a2")

        def f0_of(a):
            return mp.sqrt(mp.pi) * mp.erf(mp.sqrt(a)) / (2 * mp.sqrt(a))

        G1 = {t: target_L(t, X0, f0_of(X0)) for t in TARGETS}
        G2 = {t: target_L(t, A2, f0_of(A2)) for t in TARGETS}
        f0_cells, f0_per_target = compare(f0x0, f0a2, G1, G2, widths, "F0")
        f0_slopes = {}
        f0_cbands = {}
        f0_model = {}
        for target in TARGETS:
            f0_slopes[target] = {}
            f0_cbands[target] = {}
            f0_model[target] = {}
            for fam in ("cheb", "rat"):
                floor = 5 if fam == "cheb" else 3
                f0_slopes[target][fam] = slope_block(f0x0, f0a2, fam, target,
                                                     None)
                f0_cbands[target][fam] = {
                    "x0": constants_band(f0x0, fam, target, X0, widths, G1),
                    "a2": constants_band(f0a2, fam, target, X0 * 2, widths_a2, G2),
                }
                f0_model[target][fam] = {
                    "width_affine_x0": width_affine(
                        f0x0, fam, target, X0, widths, floor),
                    "width_affine_a2": width_affine(
                        f0a2, fam, target, X0 * 2, widths_a2, floor),
                    "matched_model_a2": matched_model_residuals(
                        f0_cells, target, fam),
                }
        summary["f0_bands"] = {
            "coverage": {"x0_missing": f0_missing1, "a2_missing": f0_missing2},
            "per_target": f0_per_target,
            "slopes": f0_slopes,
            "C_bands": f0_cbands,
            "model": f0_model,
            "ladder": {
                "x0": ladder_test(f0x0, TARGETS, G1),
                "a2": ladder_test(f0a2, TARGETS, G2),
            },
            "note": ("the paper's own x0 records (published re-count for "
                     "Chebyshev, core sweep for rational) against "
                     "results-second-edge; the same estimators as the F_1 "
                     "side, so the bands are like-for-like"),
        }
    except Exception as exc:
        summary["f0_bands"] = {"status": f"not computed: {type(exc).__name__}: {exc}"}

    # ---- write outputs
    (out_dir / "analysis.json").write_text(
        json.dumps(summary, indent=2, default=str) + "\n", encoding="utf-8")

    counts_fields = [
        "target", "width_name", "width", "width_a2", "kappa_x0", "kappa_a2",
        "kappa_rel_diff", "ell", "mu", "L_x0", "L_a2",
        "cheb_x0", "cheb_m_x0", "cheb_a2", "cheb_m_a2", "cheb_diff",
        "cheb_deliv_rel_x0", "cheb_deliv_rel_a2",
        "rat_x0", "rat_m_x0", "rat_k_x0", "rat_solver_x0", "rat_below_x0",
        "rat_a2", "rat_m_a2", "rat_k_a2", "rat_solver_a2", "rat_below_a2",
        "rat_diff", "rat_deliv_rel_x0", "rat_deliv_rel_a2",
        "Cc_x0", "Cc_a2", "Cr_x0", "Cr_a2",
    ]
    model_fields = [
        "target", "width_name", "fam", "edge", "n", "implied_C",
    ]
    comp_fields = counts_fields + [
        "cheb_model_real_a2", "cheb_model_ceil_a2", "cheb_model_resid_a2",
        "rat_model_real_a2", "rat_model_ceil_a2", "rat_model_resid_a2",
    ]
    with (out_dir / "counts_by_width_target.csv").open(
            "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=counts_fields)
        writer.writeheader()
        for target in TARGETS:
            for name in NAMES:
                e = cells[f"{target}/{name}"]
                writer.writerow({"target": target, "width_name": name,
                                 **{k: e.get(k) for k in counts_fields
                                    if k not in ("target", "width_name")}})
    with (out_dir / "comparison.csv").open(
            "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=comp_fields)
        writer.writeheader()
        for target in TARGETS:
            for name in NAMES:
                e = cells[f"{target}/{name}"]
                writer.writerow({"target": target, "width_name": name,
                                 **{k: e.get(k) for k in comp_fields
                                    if k not in ("target", "width_name")}})
    with (out_dir / "model.csv").open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=model_fields)
        writer.writeheader()
        for target in TARGETS:
            for name in NAMES:
                e = cells[f"{target}/{name}"]
                for fam, const in (("cheb", "Cc"), ("rat", "Cr")):
                    for edge in ("x0", "a2"):
                        n = e[f"{fam}_{edge}"]
                        if n is None:
                            continue
                        writer.writerow({
                            "target": target, "width_name": name, "fam": fam,
                            "edge": edge, "n": n,
                            "implied_C": e.get(f"{const}_{edge}"),
                        })

    # ---- terminal report
    print(f"F_1 second instance: x0 = {X0_TEXT}, a2 = {mp.nstr(A2, 17)}")
    print(f"F_1(x0) = {mp.nstr(f1_x0, 12)}; F_1(a2) = {mp.nstr(f1_a2, 12)}")
    print(f"L(x0): " + ", ".join(
        f"{t} {mp.nstr(L1[t], 8)}" for t in TARGETS))
    print(f"L(a2): " + ", ".join(
        f"{t} {mp.nstr(L2[t], 8)}" for t in TARGETS))
    print(f"rows on disk: {len(rows)}; missing x0 {len(missing1)}, "
          f"missing a2 {len(missing2)}")
    if missing1:
        print("  missing x0:", ", ".join(missing1))
    if missing2:
        print("  missing a2:", ", ".join(missing2))
    print("\ncounts (cheb x0->a2; rat x0->a2)")
    for target in TARGETS:
        for name in NAMES:
            e = cells[f"{target}/{name}"]
            print(f"{target:>10} {name:>5}  cheb {fmt(e['cheb_x0'])}->"
                  f"{fmt(e['cheb_a2'])}  rat {fmt(e['rat_x0'])}->"
                  f"{fmt(e['rat_a2'])}")
    print("\nper family and target")
    for target in TARGETS:
        for fam in ("cheb", "rat"):
            b = per_target[target][fam]
            print(f"{target:>10} {fam:4}: differing {b['differing_cells']}"
                  f"/{b['valid_cells']}, max |diff| {b['max_abs_diff']}, "
                  f"sum {b['sum_diff']}, parity flips {b['parity_flips']}, "
                  f"pair changes {b['pair_changes']}")
    print("\nslopes (paper estimator, thirteen widths, floors in)")
    for target in TARGETS:
        for fam in ("cheb", "rat"):
            s = summary["slopes"][target][fam]
            if s.get("status") == "incomplete":
                print(f"{target:>10} {fam:4}: INCOMPLETE")
                continue
            d = s["difference"]
            z1 = "--" if d["z_unpaired"] is None else f"{d['z_unpaired']:+.2f}"
            z2 = "--" if d["z_paired"] is None else f"{d['z_paired']:+.2f}"
            m1b = s["model_x0"].get("C_med")
            m2b = s["model_a2"].get("C_med")
            model_text = ("--" if not m1b or not m2b else
                          f"{m1b['beta']:.3f}/{m2b['beta']:.3f}")
            print(f"{target:>10} {fam:4}: x0 {s['x0']['beta']:.3f}"
                  f"+-{s['x0']['se']:.3f} (n={s['x0']['n']}), "
                  f"a2 {s['a2']['beta']:.3f}+-{s['a2']['se']:.3f}, "
                  f"diff {d['diff']:+.3f} (z {z1}/{z2}); "
                  f"model {model_text}; "
                  f"endpoints {s['x0']['endpoint']:.4f}/"
                  f"{s['a2']['endpoint']:.4f}")
    print("\ntarget-affine residuals in coefficients (floor rows excluded)")
    for fam in ("cheb", "rat"):
        for edge in ("x0", "a2"):
            b = summary["target_affine"][fam][edge]
            print(f"{fam:4} {edge}: n={b['n']}, "
                  f"max |resid| {fmt(b['max_abs_residual'], 3)}, "
                  f"mean {fmt(b['mean_residual'], 3)}, "
                  f"outside 1 {b['outside_one']}, outside 2 {b['outside_two']}")
    print("\nfree-line ladder over the three arms (paper rule, floors out)")
    for edge in ("x0", "a2"):
        for fam in ("cheb", "rat"):
            b = summary["ladder"][edge][fam]
            print(f"{fam:4} {edge}: tested {b['n_tested']}/{len(NAMES)} widths, "
                  f"failed {b['failed']}, max |resid| "
                  f"{fmt(b['max_abs_residual'], 3)}, max |2nd diff| "
                  f"{fmt(b['max_abs_second_difference'], 3)}, "
                  f"verdict {b['verdict']}")
            print(f"          tested widths: {b['tested_widths']}")
    print("\nwidth-law residuals (median-constant model)")
    for target in TARGETS:
        for fam in ("cheb", "rat"):
            for edge in ("x0", "a2"):
                w = summary["model"][target][fam][f"width_law_{edge}"]
                if w is None:
                    continue
                print(f"{target:>10} {fam:4} {edge}: C "
                      f"{w['C_median']:.4g}, max |resid| "
                      f"{w['max_abs_residual']:.3f} over {w['n']}")
    print("\nfree width law: N-1 affine in 1/rate (floors out)")
    for target in TARGETS:
        for fam in ("cheb", "rat"):
            for edge in ("x0", "a2"):
                w = summary["model"][target][fam][f"width_affine_{edge}"]
                if w.get("status"):
                    print(f"{target:>10} {fam:4} {edge}: {w['status']}")
                    continue
                print(f"{target:>10} {fam:4} {edge}: n={w['n']}, "
                      f"q={w['q']:.2f}+-{w['q_se']:.2f}, "
                      f"max |resid| {w['max_abs_residual']:.3f}")
    print("\nmatched-cell model residuals at a2 (coefficients)")
    for target in TARGETS:
        for fam in ("cheb", "rat"):
            w = summary["model"][target][fam]["matched_model_a2"]
            if w is None:
                continue
            print(f"{target:>10} {fam:4}: n={w['n']}, max "
                  f"{w['max_abs_residual']:.3f}, mean {w['mean_residual']:+.3f}, "
                  f"|r|<=1 {w['within_one']}, |r|<=2 {w['within_two']}")
    print("\ncount-implied constants over the dyadic widths, a2/x0 ratio")
    for target in TARGETS:
        for fam in ("cheb", "rat"):
            b = summary["C_ratios"][target].get(fam)
            if not b:
                continue
            print(f"{target:>10} {fam:4}: x0 {b['x0_min']:.3g}-"
                  f"{b['x0_max']:.3g}, a2 {b['a2_min']:.3g}-"
                  f"{b['a2_max']:.3g}, ratio {b['ratio_min']:.3g}-"
                  f"{b['ratio_max']:.3g} (median {b['ratio_median']:.3g}, "
                  f"expected {b['expected_ratio']:.3g})")
    if summary["f0_bands"].get("per_target"):
        print("\nF_0 bands from the paper's records (same code)")
        for target in TARGETS:
            for fam in ("cheb", "rat"):
                b = summary["f0_bands"]["per_target"][target][fam]
                print(f"{target:>10} {fam:4}: differing {b['differing_cells']}"
                      f"/{b['valid_cells']}, max |diff| {b['max_abs_diff']}, "
                      f"sum {b['sum_diff']}")
        for target in TARGETS:
            for fam in ("cheb", "rat"):
                s = summary["f0_bands"]["slopes"][target][fam]
                if s.get("status") == "incomplete":
                    print(f"{target:>10} {fam:4}: INCOMPLETE")
                    continue
                print(f"{target:>10} {fam:4} F0: x0 {s['x0']['beta']:.3f}"
                      f"+-{s['x0']['se']:.3f}, a2 {s['a2']['beta']:.3f}"
                      f"+-{s['a2']['se']:.3f}, diff "
                      f"{s['difference']['diff']:+.3f}")
        print("\nF_0 free-line ladder over the three arms (same rule)")
        for edge in ("x0", "a2"):
            for fam in ("cheb", "rat"):
                b = summary["f0_bands"]["ladder"][edge][fam]
                print(f"{fam:4} {edge}: tested {b['n_tested']}/"
                      f"{len(NAMES)} widths, failed {b['failed']}, "
                      f"max |resid| {fmt(b['max_abs_residual'], 3)}, "
                      f"max |2nd diff| "
                      f"{fmt(b['max_abs_second_difference'], 3)}, "
                      f"verdict {b['verdict']}")
        print("\nF_0 free width law and matched-cell residuals")
        for target in TARGETS:
            for fam in ("cheb", "rat"):
                for edge in ("x0", "a2"):
                    w = summary["f0_bands"]["model"][target][fam][
                        f"width_affine_{edge}"]
                    if w.get("status"):
                        print(f"{target:>10} {fam:4} {edge}: {w['status']}")
                    else:
                        print(f"{target:>10} {fam:4} {edge}: n={w['n']}, "
                              f"max |resid| {w['max_abs_residual']:.3f}")
                w = summary["f0_bands"]["model"][target][fam][
                    "matched_model_a2"]
                if w:
                    print(f"{target:>10} {fam:4} matched a2: max "
                          f"{w['max_abs_residual']:.3f}, mean "
                          f"{w['mean_residual']:+.3f}, |r|<=1 "
                          f"{w['within_one']}/{w['n']}")
    else:
        print("\nF_0 bands: " + str(summary["f0_bands"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
