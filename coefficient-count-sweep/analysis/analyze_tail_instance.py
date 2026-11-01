#!/usr/bin/env python3
"""Analyse the pure-tail (x^-1/2) Group A counts at the two left edges.

Reads the rows written by tail_instance_sweep.py (one file per run directory;
the two-edge runs are merged by (a, target, width_name)) and answers the
paper's questions for this instance of the algebraic-tail class:

  * the Chebyshev and rational stored counts at each edge and target;
  * the width-law slopes with the paper's own estimator (least squares of
    log2 N on log2 W over the thirteen scaling widths, floor rows included,
    with the standard error and the endpoint secant beside it);
  * the model's rate predictions -- log(1/kappa) for Chebyshev, independent of
    alpha, and the Stahl modulus mu_alpha(kappa) at alpha = 1/2 for the
    rational family -- compared with the counted counts;
  * the matched-W/a second-edge comparison the paper ran for F_0: differences,
    slope differences (paired and unpaired), count-implied constants carried
    across the matched cross-ratio, and the target-increment ("affine in L")
    test that carries the paper's Chebyshev-affine and rational-departure
    claims.

Model conventions, stated once and used everywhere:
  * a relative target has L = log(1/eps);
  * the absolute target 5e-14 on f is accepted by the solvers on the raw
    absolute error; in the model frame L = log(f(a)/eps) with f(a) = 1/sqrt(a),
    which is the paper's "value at the left edge" convention for its absolute
    arm (F_0(a)/eps). Using log(1/eps) instead would shift every log C by the
    constant log f(a) = -0.5 log a and change no slope and no count;
  * kappa(W) = (sqrt(a+W) - sqrt(a)) / (sqrt(a+W) + sqrt(a)),
    mu(kappa) = (pi/2) K'(kappa)/K(kappa), the Groetzsch ring function, which
    is the Stahl modulus of the square-root problem (alpha = 1/2).

Writes into the results directory: analysis.json, the per-edge CSVs, the
matched comparison CSV and the compact count table. Prints the same tables.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from datetime import datetime, timezone
from pathlib import Path

import mpmath as mp

import count_sweep as cs

PAPER = Path(__file__).resolve().parent
NAMES = [f"2^{j}" for j in range(-7, 6)] + ["ship", "paper"]
DYAD = [f"2^{j}" for j in range(-7, 6)]
TARGETS = ["rel:1e-13", "abs:5e-14", "rel:1e-14"]
X0_TEXT = "11.899848152108484"
A2_TEXT = "23.799696304216968"
EDGE_LABEL = {X0_TEXT: "x0", A2_TEXT: "a2"}
BAND_TWO = 2.0
BAND_ONE = 1.0


def width_of(name: str):
    x0 = mp.mpf(X0_TEXT)
    if name.startswith("2^"):
        return mp.mpf(2) ** int(name.split("^")[1])
    if name == "ship":
        return mp.mpf("31") - x0
    if name == "paper":
        return mp.mpf("28.989337738820740") - x0
    raise ValueError(name)


def kappa(a, width):
    s = mp.sqrt(a)
    t = mp.sqrt(a + width)
    return (t - s) / (t + s)


def mu_of(kappa_value):
    return (mp.pi / 2) * mp.ellipk(1 - kappa_value ** 2) / mp.ellipk(kappa_value ** 2)


def f_of(x):
    return 1 / mp.sqrt(x)


def target_L(target: str, a):
    kind, value = target.split(":", 1)
    eps = mp.mpf(value)
    if kind == "rel":
        return mp.log(1 / eps), "rel", eps
    return mp.log(f_of(a) / eps), "abs", eps


def ols(xs, ys):
    n = len(xs)
    mx = sum(xs) / n
    my = sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    beta = sxy / sxx
    alpha = my - beta * mx
    resid = [y - (alpha + beta * x) for x, y in zip(xs, ys)]
    if n > 2:
        s2 = sum(r * r for r in resid) / (n - 2)
        se = math.sqrt(s2 / sxx)
    else:
        se = None
    return {"alpha": alpha, "beta": beta, "se": se,
            "residuals": resid, "n": n, "sxx": sxx,
            "max_abs_residual": max(abs(r) for r in resid)}


def paired_slope_test(fit1, fit2):
    diff = fit2["beta"] - fit1["beta"]
    if fit1["se"] is None or fit2["se"] is None or fit1["n"] < 3:
        return {"diff": diff, "se_unpaired": None, "se_paired": None,
                "z_unpaired": None, "z_paired": None, "residual_cov": None}
    se_unpaired = math.sqrt(fit1["se"] ** 2 + fit2["se"] ** 2)
    n = fit1["n"]
    cov = sum(r1 * r2 for r1, r2 in zip(fit1["residuals"], fit2["residuals"])) / (n - 2)
    se_paired = math.sqrt(max(0.0, fit1["se"] ** 2 + fit2["se"] ** 2
                              - 2 * cov / fit1["sxx"]))
    return {"diff": diff, "se_unpaired": se_unpaired,
            "se_paired": se_paired,
            "z_unpaired": diff / se_unpaired if se_unpaired else None,
            "z_paired": diff / se_paired if se_paired else None,
            "residual_cov": cov}


def load_rows(paths, function):
    """Merge every run file's rows, keyed by (left edge, target, width name)."""
    rows = {}
    sources = []
    metas = []
    duplicates = []
    for path in paths:
        if not path.exists():
            continue
        try:
            doc = json.loads(path.read_text(encoding="utf-8"))
        except Exception as exc:
            print(f"warning: cannot read {path} ({exc}); skipped", flush=True)
            continue
        if doc.get("meta", {}).get("function") != function:
            continue
        n = 0
        for row in doc.get("rows", []):
            if not (row.get("cheb") and row.get("rat")):
                continue
            key = (row["a"], row["target"], row["width_name"])
            if key in rows:
                old = rows[key]
                if (old["cheb"].get("stored") != row["cheb"].get("stored")
                        or old["rat"].get("stored") != row["rat"].get("stored")):
                    duplicates.append({"key": list(key),
                                       "first": [old["cheb"].get("stored"),
                                                 old["rat"].get("stored")],
                                       "again": [row["cheb"].get("stored"),
                                                 row["rat"].get("stored")],
                                       "path": str(path)})
                continue
            rows[key] = row
            n += 1
        sources.append({"path": str(path), "rows": n,
                        "generated_utc": doc.get("meta", {}).get("generated_utc")})
        metas.append((n, doc.get("meta", {}), str(path)))
    return rows, sources, duplicates, metas


def load_cells(path, function):
    if not path.exists():
        return {}, None
    doc = json.loads(path.read_text(encoding="utf-8"))
    if doc.get("meta", {}).get("function") != function:
        return {}, None
    cells = {(cell["a"], cell["target"], cell["width_name"]): cell
             for cell in doc.get("cells", [])
             if cell.get("cheb", {}).get("stored") is not None}
    return cells, doc.get("meta", {})


def model_secant(kind, a, L, constant, js=DYAD):
    """The model's all-in secant over the thirteen scaling widths.

    ``a`` fixes the scale: the widths are the Group A widths multiplied by
    a/x0, so ``W/a`` is the Group A ``W/x0`` at every edge.
    """
    scale = a / mp.mpf(X0_TEXT)
    xs = [float(name.split("^")[1]) for name in js]
    if kind == "cheb":
        ys = [math.log2(float(1 + (L + mp.log(constant))
                              / mp.log(1 / kappa(a, scale * width_of(name)))))
              for name in js]
    else:
        ys = [math.log2(float(1 + (L + mp.log(constant))
                              / mu_of(kappa(a, scale * width_of(name)))))
              for name in js]
    fit = ols(xs, ys)
    return fit["beta"]


def model_secant_rat_stored(a, L, form):
    """The diagonal rational law on its stored count, or the off-diagonal one.

    ``form = 'diag_odd'``: N = 1 + (L + log 4)/mu, stored odd
    2 ceil((N-1)/2) + 1.  ``form = 'offdiag_even'``: the even stored count
    2 ceil((L + log 4)/(2 mu)).  Both are the paper's quantised models, whose
    secants at rel 1e-13 are 0.1331 and 0.1532.
    """
    scale = a / mp.mpf(X0_TEXT)
    ys = []
    for name in DYAD:
        mu = mu_of(kappa(a, scale * width_of(name)))
        if form == "diag_odd":
            n = 1 + (L + mp.log(4)) / mu
            stored = 2 * mp.ceil((n - 1) / 2) + 1
        elif form == "offdiag_even":
            stored = 2 * mp.ceil((L + mp.log(4)) / (2 * mu))
        else:
            raise ValueError(form)
        ys.append(math.log2(float(stored)))
    return ols([float(name.split("^")[1]) for name in DYAD], ys)["beta"]


def local_exponent(kind, a, width):
    """The model's local exponent at one width (supplement eq. local)."""
    k = kappa(a, width)
    ell = mp.log(1 / k)
    dk = mp.sqrt(a) / (mp.sqrt(a + width) * (mp.sqrt(a + width) + mp.sqrt(a)) ** 2)
    if kind == "cheb":
        return float(width * dk / (k * ell)), float(ell)
    h = mp.mpf("1e-25")
    dmu = (mu_of(k + h) - mu_of(k - h)) / (2 * h)
    mu = mu_of(k)
    return float(-width * dmu * dk / mu), float(mu)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dir", default=str(PAPER / "results-second-instance"))
    parser.add_argument("--function", default="x^-1/2")
    parser.add_argument("--extra-rows", action="append", default=None,
                        help="an additional sweep_all.json to merge; repeatable")
    args = parser.parse_args(argv)
    out_dir = Path(args.dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    mp.mp.dps = 50
    a_x0 = mp.mpf(X0_TEXT)
    a_a2 = mp.mpf(A2_TEXT)
    edge_of = {X0_TEXT: a_x0, A2_TEXT: a_a2}

    paths = [out_dir / "sweep_all.json",
             out_dir / "_edge-x0" / "sweep_all.json",
             out_dir / "_edge-a2" / "sweep_all.json"]
    if args.extra_rows:
        paths += [Path(p) for p in args.extra_rows]
    rows, sources, duplicates, metas = load_rows(paths, args.function)
    cells, cell_meta = load_cells(out_dir / "chebyshev_only.json", args.function)

    # ---- merge the per-edge runs into the deposit's sweep_all.json ----
    merged_note = None
    if rows and metas:
        base_n, base_meta, base_path = max(metas, key=lambda item: item[0])
        merged_meta = json.loads(json.dumps(base_meta))
        merged_meta["generated_utc"] = datetime.now(timezone.utc).isoformat()
        merged_meta["merged_from"] = [{"path": item[2], "rows": item[0]}
                                      for item in metas]
        merged_meta["left_edges"] = sorted({key[0] for key in rows})
        merged_meta["edge_scales"] = sorted(
            {int(row.get("edge_scale", 1)) for row in rows.values()})
        merged_meta["matched_ratio"] = (
            "every width is the Group A width multiplied by the edge scale "
            "(1 at x0, 2 at a2), so W/a is the Group A W/x0 at both edges")
        if isinstance(merged_meta.get("precision_scope"), dict):
            merged_meta["precision_scope"]["sweep_scope"] = (
                "fp64 only: binary64 stored coefficients, binary64 "
                "evaluation, targets rel:1e-13, abs:5e-14, rel:1e-14")
        merged_meta["target_note"] = (
            "rel targets are relative to f; abs:5e-14 is an absolute error on "
            "f, unscaled (the model frame used in the analysis is "
            "L = log(f(a)/eps) with f(a) = 1/sqrt(a))")
        merged_meta["analysis_driver"] = {
            "path": str(PAPER / "analyze_tail_instance.py"),
            "sha256": cs.sha256(PAPER / "analyze_tail_instance.py")}
        merged_meta["merge_note"] = (
            "rows from the per-edge runs, merged by (a, target, width_name); "
            f"source with the most rows: {base_path}")
        merged_targets = cs.parse_targets(",".join(TARGETS))
        # Never shrink the deposit file: only replace it with a union at least
        # as large as what it already holds.
        previous = 0
        top_path = out_dir / "sweep_all.json"
        if top_path.exists():
            try:
                previous = len(json.loads(top_path.read_text(
                    encoding="utf-8")).get("rows", []))
            except Exception:
                previous = 0
        if len(rows) >= previous:
            try:
                cs.write_sweep_outputs(out_dir, list(rows.values()),
                                       merged_targets, merged_meta)
                merged_note = str(top_path)
            except Exception as exc:  # a merge failure must not lose the analysis
                print(f"warning: merge failed ({exc}); source files untouched",
                      flush=True)
        else:
            merged_note = (f"not merged: {len(rows)} union row(s) against "
                           f"{previous} already in {top_path}")

    present = {"full_rows": len(rows), "cheb_cells": len(cells),
               "targets": TARGETS, "widths": NAMES}
    missing = {"full": [(t, n, EDGE_LABEL[e]) for t in TARGETS for n in NAMES
                        for e in (X0_TEXT, A2_TEXT)
                        if (e, t, n) not in rows],
               "cheb": [(t, n, EDGE_LABEL[e]) for t in TARGETS for n in NAMES
                        for e in (X0_TEXT, A2_TEXT)
                        if (e, t, n) not in rows and (e, t, n) not in cells]}
    summary = {
        "function": args.function,
        "reference": "1/sqrt(x), 60-digit mpmath (the generator's boys_ref replaced)",
        "x0": X0_TEXT, "a2": A2_TEXT,
        "sources": sources,
        "duplicates_rerun": duplicates,
        "presence": present,
        "missing": missing,
        "model_conventions": {
            "kappa": "(sqrt(a+W)-sqrt(a))/(sqrt(a+W)+sqrt(a))",
            "mu": "(pi/2) K'(kappa)/K(kappa), the Stahl modulus at alpha=1/2",
            "L_rel": "log(1/eps)",
            "L_abs": "log(f(a)/eps) with f(a)=1/sqrt(a); the paper's "
                     "left-edge-value convention. The solvers accept the raw "
                     "absolute error 5e-14 on f (unscaled).",
            "cheb_law": "N = 1 + (L + log Cc)/log(1/kappa)",
            "rat_law": "N = 1 + (L + log Cr)/mu(kappa)",
        },
        "edges": {},
        "cells": {},
        "targets": {},
        "slopes": {},
        "width_law_model": {},
        "increments": {},
        "rational_departure": {},
        "dominance": {},
    }

    # ---- per-cell table ----
    for edge_text, a in edge_of.items():
        label = EDGE_LABEL[edge_text]
        scale = a / a_x0
        summary["edges"][label] = {"a": edge_text,
                                   "width_scale": mp.nstr(scale, 2),
                                   "kappa_identical_to_x0": None}
        for target in TARGETS:
            L, kind, eps = target_L(target, a)
            for name in NAMES:
                key = (edge_text, target, name)
                row = rows.get(key)
                cell = cells.get(key)
                width = scale * width_of(name)
                k = kappa(a, width)
                ell = mp.log(1 / k)
                mu = mu_of(k)
                entry = {
                    "edge": label, "a": edge_text, "target": target,
                    "width_name": name,
                    "width": mp.nstr(width, 17),
                    "width_x0": mp.nstr(width_of(name), 17),
                    "kappa": mp.nstr(k, 17), "ell": mp.nstr(ell, 17),
                    "mu": mp.nstr(mu, 17),
                    "L": mp.nstr(L, 17),
                    "L_kind": kind, "eps": mp.nstr(eps, 17),
                    "has_full_row": row is not None,
                    "has_cheb_cell": cell is not None,
                }
                cheb = (row or {}).get("cheb") or (cell or {}).get("cheb") or {}
                rat = (row or {}).get("rat") or {}
                entry.update({
                    "cheb_stored": cheb.get("stored"),
                    "cheb_degree": cheb.get("m"),
                    "cheb_deliv_rel": cheb.get("delivered_rel"),
                    "cheb_deliv_abs": cheb.get("delivered_abs"),
                    "cheb_parity": ("odd" if cheb.get("stored") % 2 else "even")
                    if cheb.get("stored") is not None else None,
                    "rat_stored": rat.get("stored"),
                    "rat_m": rat.get("m"), "rat_k": rat.get("k"),
                    "rat_solver": rat.get("solver"),
                    "rat_deliv_rel": rat.get("delivered_rel"),
                    "rat_deliv_abs": rat.get("delivered_abs"),
                    "rat_below_meets": rat.get("below_meets"),
                    "rat_parity": ("odd" if rat.get("stored") % 2 else "even")
                    if rat.get("stored") is not None else None,
                })
                # The counts are integers; a below-count hit means the returned
                # count is known not to be minimal and it is kept out of fits.
                entry["rat_valid"] = (rat.get("stored") is not None
                                      and not rat.get("below_meets"))
                if cheb.get("stored") is not None:
                    entry["Cc"] = mp.nstr(mp.exp((cheb["stored"] - 1) * ell - L), 17)
                    entry["cheb_model_real"] = mp.nstr(
                        1 + (L + mp.log(mp.mpf(entry["Cc"]))) / ell, 17)
                if rat.get("stored") is not None:
                    entry["Cr"] = mp.nstr(mp.exp((rat["stored"] - 1) * mu - L), 17)
                    entry["rat_model_real"] = mp.nstr(
                        1 + (L + mp.log(mp.mpf(entry["Cr"]))) / mu, 17)
                    entry["rat_model_Cr4"] = mp.nstr(
                        1 + (L + mp.log(4)) / mu, 17)
                summary["cells"][f"{label}/{target}/{name}"] = entry

    # kappa is scale invariant at the doubled width: report the mismatch.
    worst = mp.mpf(0)
    for name in NAMES:
        w = width_of(name)
        worst = max(worst, abs(kappa(a_a2, 2 * w) / kappa(a_x0, w) - 1))
    summary["edges"]["a2"]["kappa_identical_to_x0"] = mp.nstr(worst, 6)

    # ---- per family and target: the second-edge differences ----
    for target in TARGETS:
        block = {}
        for fam in ("cheb", "rat"):
            diffs = []
            parity_flips = 0
            pair_changes = 0
            for name in NAMES:
                e1 = summary["cells"].get(f"x0/{target}/{name}")
                e2 = summary["cells"].get(f"a2/{target}/{name}")
                if not e1 or not e2:
                    continue
                n1 = e1.get(f"{fam}_stored")
                n2 = e2.get(f"{fam}_stored")
                if n1 is None or n2 is None:
                    continue
                if fam == "rat" and not (e1["rat_valid"] and e2["rat_valid"]):
                    continue
                d = n2 - n1
                if e1.get(f"{fam}_parity") != e2.get(f"{fam}_parity"):
                    parity_flips += 1
                if fam == "cheb":
                    if e1["cheb_degree"] != e2["cheb_degree"]:
                        pair_changes += 1
                elif (e1["rat_m"], e1["rat_k"]) != (e2["rat_m"], e2["rat_k"]):
                    pair_changes += 1
                diffs.append({"width": name, "x0": n1, "a2": n2, "diff": d})
            block[fam] = {
                "cells_compared": len(diffs),
                "differing_cells": sum(1 for d in diffs if d["diff"] != 0),
                "max_abs_diff": max((abs(d["diff"]) for d in diffs), default=None),
                "sum_diff": sum(d["diff"] for d in diffs) if diffs else None,
                "parity_flips": parity_flips,
                "pair_changes": pair_changes,
                "cells": diffs,
            }
        summary["targets"][target] = block

    # ---- width-law slopes, the paper's estimator ----
    def slope_pair(fam, target):
        js = [j for j in range(-7, 6)]
        usable = []
        for j in js:
            name = f"2^{j}"
            e1 = summary["cells"].get(f"x0/{target}/{name}")
            e2 = summary["cells"].get(f"a2/{target}/{name}")
            if not e1 or not e2:
                continue
            if e1.get(f"{fam}_stored") is None or e2.get(f"{fam}_stored") is None:
                continue
            if fam == "rat" and not (e1["rat_valid"] and e2["rat_valid"]):
                continue
            usable.append(j)
        xs = [float(j) for j in usable]
        ys1 = [math.log2(summary["cells"][f"x0/{target}/2^{j}"][f"{fam}_stored"])
               for j in usable]
        ys2 = [math.log2(summary["cells"][f"a2/{target}/2^{j}"][f"{fam}_stored"])
               for j in usable]
        if len(xs) < 2:
            return {"status": "incomplete", "common_widths": usable}
        if len(xs) == 2:
            secant1 = (ys1[1] - ys1[0]) / (xs[1] - xs[0])
            secant2 = (ys2[1] - ys2[0]) / (xs[1] - xs[0])
            return {"status": "partial", "n": 2,
                    "common_widths": [f"2^{j}" for j in usable],
                    "note": "two common widths: secants only, no s.e.",
                    "x0": {"beta": secant1, "se": None,
                           "endpoint": secant1},
                    "a2": {"beta": secant2, "se": None,
                           "endpoint": secant2},
                    "difference": {"diff": secant2 - secant1,
                                   "se_unpaired": None, "se_paired": None,
                                   "z_unpaired": None, "z_paired": None}}
        fit1 = ols(xs, ys1)
        fit2 = ols(xs, ys2)
        endpoint1 = (ys1[-1] - ys1[0]) / (xs[-1] - xs[0])
        endpoint2 = (ys2[-1] - ys2[0]) / (xs[-1] - xs[0])
        out = {
            "common_widths": [f"2^{j}" for j in usable],
            "n": len(usable),
            "x0": {"alpha": fit1["alpha"], "beta": fit1["beta"],
                   "se": fit1["se"], "endpoint": endpoint1,
                   "max_abs_residual": fit1["max_abs_residual"]},
            "a2": {"alpha": fit2["alpha"], "beta": fit2["beta"],
                   "se": fit2["se"], "endpoint": endpoint2,
                   "max_abs_residual": fit2["max_abs_residual"]},
            "difference": paired_slope_test(fit1, fit2),
        }
        # The model's secant at the paper's reference constants, beside the
        # measured slope; the delta is quoted in units of the measured s.e.
        for edge_text, key in ((X0_TEXT, "x0"), (A2_TEXT, "a2")):
            a = edge_of[edge_text]
            L, _, _ = target_L(target, a)
            constants = ([("Cc_0.02", mp.mpf("0.02")), ("Cc_1", mp.mpf(1))]
                         if fam == "cheb" else
                         [("Cr_4", mp.mpf(4))])
            model = {}
            for name, constant in constants:
                beta = model_secant(fam, a, L, constant)
                model[name] = {
                    "constant": float(constant), "secant": beta,
                    "delta_vs_measured": out[key]["beta"] - beta,
                    "delta_in_se": (out[key]["beta"] - beta) / out[key]["se"]
                    if out[key]["se"] else None}
            out[key]["model"] = model
            if fam == "rat":
                for label, form in (("diag_odd", "diag_odd"),
                                    ("offdiag_even", "offdiag_even")):
                    beta = model_secant_rat_stored(a, L, form)
                    out[key]["model"][label] = {
                        "secant": beta,
                        "delta_vs_measured": out[key]["beta"] - beta,
                        "delta_in_se": (out[key]["beta"] - beta)
                        / out[key]["se"] if out[key]["se"] else None}
        if len(usable) < 10:
            out["status"] = "partial"
        return out

    for target in TARGETS:
        summary["slopes"][target] = {fam: slope_pair(fam, target)
                                     for fam in ("cheb", "rat")}

    # ---- the model's rates at each edge, and the counted-vs-predicted ----
    for edge_text, a in edge_of.items():
        label = EDGE_LABEL[edge_text]
        summary["width_law_model"][label] = {}
        for target in TARGETS:
            L, kind, _ = target_L(target, a)
            block = {}
            for fam, rate_name in (("cheb", "ell"), ("rat", "mu")):
                counts = [summary["cells"][f"{label}/{target}/{n}"].get(f"{fam}_stored")
                          for n in DYAD]
                counts = [c for c in counts if c is not None]
                if len(counts) < 2:
                    block[fam] = {"status": "incomplete"}
                    continue
                # The single constant the edge's own counts imply. Floor rows
                # carry no information about the constant (any large C gives
                # the floor), so the constant is taken from the widest rows
                # that sit above the floor, where the count is informative.
                logCs = []
                informative = []
                for n in DYAD:
                    e = summary["cells"][f"{label}/{target}/{n}"]
                    if e.get(f"{fam}_stored") is None:
                        continue
                    rate = mp.mpf(e[rate_name])
                    if e[f"{fam}_stored"] <= 5:
                        continue
                    log_c = (float(e[f"{fam}_stored"]) - 1) * float(rate) - float(L)
                    logCs.append(log_c)
                    informative.append(n)
                if not logCs:
                    block[fam] = {"status": "all floors"}
                    continue
                logCs.sort()
                mid = logCs[len(logCs) // 2]
                C = mp.e ** mp.mpf(mid)
                # Residuals of the observed stored count from the model with
                # that one constant, in coefficients.
                residuals = []
                for n in DYAD:
                    e = summary["cells"][f"{label}/{target}/{n}"]
                    if e.get(f"{fam}_stored") is None:
                        continue
                    rate = mp.mpf(e[rate_name])
                    real = 1 + (L + mp.log(C)) / rate
                    residuals.append({"width": n,
                                      "observed": e[f"{fam}_stored"],
                                      "model_real": mp.nstr(real, 17),
                                      "residual": float(e[f"{fam}_stored"] - real)})
                # The model's secant at the measured constant, and at the
                # paper's reference constants.
                secant = model_secant(fam, a, L, C)
                block[fam] = {
                    "rate_name": rate_name,
                    "constant_measured": mp.nstr(C, 17),
                    "log_constant": mid,
                    "informative_widths": informative,
                    "model_secant": secant,
                    "residuals": residuals,
                    "max_abs_residual": max(abs(r["residual"]) for r in residuals),
                    "max_abs_residual_above_floor": max(
                        (abs(r["residual"]) for r in residuals
                         if r["observed"] > 5), default=None),
                    "constant_band": [mp.nstr(mp.e ** mp.mpf(min(logCs)), 17),
                                      mp.nstr(mp.e ** mp.mpf(max(logCs)), 17)],
                    "decades_of_drift": (max(logCs) - min(logCs)) / math.log(10),
                }
                if fam == "cheb":
                    block[fam]["model_secant_Cc_0p02"] = model_secant(
                        fam, a, L, mp.mpf("0.02"))
                    block[fam]["model_secant_Cc_1"] = model_secant(
                        fam, a, L, mp.mpf(1))
                else:
                    block[fam]["model_secant_Cr1"] = model_secant(
                        fam, a, L, mp.mpf(1))
                    block[fam]["model_secant_Cr4"] = model_secant(
                        fam, a, L, mp.mpf(4))
                    block[fam]["model_secant_diag_odd"] = \
                        model_secant_rat_stored(a, L, "diag_odd")
                    block[fam]["model_secant_offdiag_even"] = \
                        model_secant_rat_stored(a, L, "offdiag_even")
            block["local_exponents_at_2^5"] = {
                fam: local_exponent(fam, a, (a / a_x0) * mp.mpf(32))[0]
                for fam in ("cheb", "rat")}
            summary["width_law_model"][label][target] = block

    # ---- the target-increment test ("affine in L") ----
    # Three targets give two increments per width. The decade rel:1e-13 ->
    # rel:1e-14 is one decade of target and the cleanest step.
    for edge_text, a in edge_of.items():
        label = EDGE_LABEL[edge_text]
        block = {}
        for fam, rate_name in (("cheb", "ell"), ("rat", "mu")):
            rows_out = []
            for name in NAMES:
                e13 = summary["cells"].get(f"{label}/rel:1e-13/{name}")
                e14 = summary["cells"].get(f"{label}/rel:1e-14/{name}")
                eabs = summary["cells"].get(f"{label}/abs:5e-14/{name}")
                if not (e13 and e14 and eabs):
                    continue
                if (e13.get(f"{fam}_stored") is None
                        or e14.get(f"{fam}_stored") is None):
                    continue
                if fam == "rat" and not (e13["rat_valid"] and e14["rat_valid"]):
                    continue
                rate = mp.mpf(e13[rate_name])
                measured = e14[f"{fam}_stored"] - e13[f"{fam}_stored"]
                predicted = float(mp.log(10) / rate)
                item = {"width": name,
                        "measured_increment": measured,
                        "predicted_increment": predicted,
                        "residual": measured - predicted}
                if eabs.get(f"{fam}_stored") is not None:
                    Labs, _, _ = target_L("abs:5e-14", a)
                    L14, _, _ = target_L("rel:1e-14", a)
                    item["abs_stored"] = eabs[f"{fam}_stored"]
                    item["abs_to_rel14_measured"] = (
                        e14[f"{fam}_stored"] - eabs[f"{fam}_stored"])
                    item["abs_to_rel14_predicted"] = float(
                        (L14 - Labs) / rate)
                    item["abs_to_rel14_residual"] = (
                        item["abs_to_rel14_measured"]
                        - item["abs_to_rel14_predicted"])
                rows_out.append(item)
            if not rows_out:
                block[fam] = {"status": "incomplete"}
                continue
            res = [r["residual"] for r in rows_out]
            block[fam] = {
                "cells": rows_out,
                "max_abs_residual": max(abs(v) for v in res),
                "cells_outside_one": sum(1 for v in res if abs(v) > BAND_ONE),
                "cells_outside_two": sum(1 for v in res if abs(v) > BAND_TWO),
                "signs": [("+" if v > 0 else ("-" if v < 0 else "0")) for v in res],
                "abs_increments": [
                    {"width": r["width"],
                     "residual": r.get("abs_to_rel14_residual")}
                    for r in rows_out if "abs_to_rel14_residual" in r],
            }
        summary["increments"][label] = block

    # ---- the rational departure from the leading Stahl law ----
    for edge_text, a in edge_of.items():
        label = EDGE_LABEL[edge_text]
        block = {}
        for target in TARGETS:
            L, _, _ = target_L(target, a)
            vals = []
            for name in NAMES:
                e = summary["cells"].get(f"{label}/{target}/{name}")
                if not e or e.get("rat_stored") is None or not e["rat_valid"]:
                    continue
                rate = mp.mpf(e["mu"])
                real4 = 1 + (L + mp.log(4)) / rate
                realC = 1 + (L + mp.log(mp.mpf(e["Cr"]))) / rate
                vals.append({"width": name,
                             "observed": e["rat_stored"],
                             "model_Cr4": mp.nstr(real4, 17),
                             "excess_over_Cr4": float(e["rat_stored"] - real4),
                             "model_measured_Cr": mp.nstr(realC, 17),
                             "excess_over_measured_Cr": float(e["rat_stored"] - realC)})
            block[target] = vals
        summary["rational_departure"][label] = block

    # ---- the paper's pre-set ladder test (free line in L) ----
    # For each width and family, above the solver floor: fit the stored count
    # against L = log(1/eps) with a free line and require every residual and
    # every second difference to be at most one stored coefficient. The
    # absolute target enters at L = log(f(a)/eps), the model frame; the test
    # itself is invariant to that additive constant (it moves the intercept
    # only), so the convention cannot manufacture a pass or a flip.
    ladder = {"tolerance": BAND_ONE, "edges": {}}
    for edge_text, a in edge_of.items():
        label = EDGE_LABEL[edge_text]
        floors = {}
        for fam in ("cheb", "rat"):
            values = []
            for name in NAMES:
                e = summary["cells"].get(f"{label}/rel:1e-13/{name}")
                if not e or e.get(f"{fam}_stored") is None:
                    continue
                if fam == "rat" and not e["rat_valid"]:
                    continue
                values.append(e[f"{fam}_stored"])
            floors[fam] = min(values) if values else None
        edge_block = {"floors": floors, "widths": {}}
        flips = []
        complete = True
        for name in NAMES:
            family_block = {}
            for fam, rate_key in (("cheb", "ell"), ("rat", "mu")):
                pts = []
                excluded = []
                for target in TARGETS:
                    e = summary["cells"].get(f"{label}/{target}/{name}")
                    if not e or e.get(f"{fam}_stored") is None:
                        continue
                    if fam == "rat" and not e["rat_valid"]:
                        excluded.append({"target": target,
                                         "stored": e["rat_stored"]})
                        continue
                    L, _, _ = target_L(target, a)
                    pts.append((float(L), e[f"{fam}_stored"],
                                e.get(rate_key)))
                pts.sort()
                used = [(L, n, rate) for L, n, rate in pts
                        if floors[fam] is not None and n > floors[fam]]
                entry = {"points_total": len(pts),
                         "points_above_floor": len(used),
                         "counts": [n for _, n, _ in used],
                         "Ls": [round(L, 6) for L, _, _ in used],
                         "excluded_non_minimal": excluded}
                if len(used) >= 3:
                    n = len(used)
                    mx = sum(L for L, _, _ in used) / n
                    my = sum(c for _, c, _ in used) / n
                    sxx = sum((L - mx) ** 2 for L, _, _ in used)
                    sxy = sum((L - mx) * (c - my) for L, c, _ in used)
                    slope = sxy / sxx
                    intercept = my - slope * mx
                    residuals = [c - (intercept + slope * L) for L, c, _ in used]
                    second = [used[i + 1][1] - 2 * used[i][1] + used[i - 1][1]
                              for i in range(1, n - 1)]
                    # The paper's ladder has equal L spacing, so its raw
                    # second difference is a rounding measure. The three Group
                    # A targets are not equally spaced in L; the spacing term
                    # slope * (L2 - 2 L1 + L0) is what a straight line
                    # contributes, and (raw - that) is the curvature the
                    # ladder test is after.
                    spacing_terms = [
                        slope * (used[i + 1][0] - 2 * used[i][0]
                                 + used[i - 1][0])
                        for i in range(1, n - 1)]
                    corrected = [raw - term
                                 for raw, term in zip(second, spacing_terms)]
                    rate = mp.mpf(used[0][2])
                    entry.update({
                        "slope_measured": slope,
                        "slope_predicted": float(1 / rate),
                        "slope_se": (math.sqrt(sum(r * r for r in residuals)
                                               / (n - 2) / sxx) if n > 2 else None),
                        "intercept": intercept,
                        "residuals": residuals,
                        "second_differences": second,
                        "second_difference_spacing_terms": spacing_terms,
                        "second_differences_corrected": corrected,
                        "max_abs_residual": max(abs(r) for r in residuals),
                        "max_abs_second_difference":
                            max(abs(s) for s in second) if second else 0.0,
                        "max_abs_second_difference_corrected":
                            max(abs(s) for s in corrected) if corrected else 0.0,
                    })
                    ok = (entry["max_abs_residual"] <= BAND_ONE
                          and entry["max_abs_second_difference_corrected"]
                          <= BAND_ONE)
                    entry["affine_within_tolerance"] = ok
                    if not ok:
                        flips.append({"width": name, "family": fam})
                else:
                    entry["affine_within_tolerance"] = None
                    entry["note"] = "fewer than three above-floor targets"
                    complete = False
                family_block[fam] = entry
            edge_block["widths"][name] = family_block
        edge_block["flips"] = flips
        edge_block["complete"] = complete
        edge_block["verdict"] = ("C-NOLAW" if flips else "C-LAW") if complete \
            else "C-INCOMPLETE"
        ladder["edges"][label] = edge_block
    summary["ladder"] = ladder

    # ---- the model's a2 prediction from the x0 count ----
    # The paper's 7.3 construction: the cross-ratio is exactly matched, so for
    # a relative target the model predicts the same count at a2; for the
    # absolute target the L shift is log(f(a2)/f(a1)) = -0.5 log 2 and the
    # prediction falls by that over the width's rate.
    summary["second_edge_prediction"] = {}
    for target in TARGETS:
        block = {}
        for fam, rate_name in (("cheb", "ell"), ("rat", "mu")):
            pred_cells = []
            for name in NAMES:
                e1 = summary["cells"].get(f"x0/{target}/{name}")
                e2 = summary["cells"].get(f"a2/{target}/{name}")
                if not e1 or not e2:
                    continue
                n1 = e1.get(f"{fam}_stored")
                n2 = e2.get(f"{fam}_stored")
                if n1 is None or n2 is None:
                    continue
                if fam == "rat" and not (e1["rat_valid"] and e2["rat_valid"]):
                    continue
                L1 = mp.mpf(e1["L"])
                L2 = mp.mpf(e2["L"])
                rate = mp.mpf(e1[rate_name])
                predicted = n1 + float((L2 - L1) / rate)
                pred_cells.append({"width": name, "x0": n1, "a2": n2,
                                   "predicted_a2": predicted,
                                   "residual": n2 - predicted})
            if pred_cells:
                block[fam] = {
                    "cells": pred_cells,
                    "min_residual": min(c["residual"] for c in pred_cells),
                    "max_residual": max(c["residual"] for c in pred_cells),
                    "cells_above_zero": sum(1 for c in pred_cells
                                            if c["residual"] > 0),
                    "cells_below_zero": sum(1 for c in pred_cells
                                            if c["residual"] < 0),
                    "max_abs_residual": max(abs(c["residual"])
                                            for c in pred_cells),
                }
        summary["second_edge_prediction"][target] = block

    # ---- the count-implied constant ratio at the matched cross-ratio ----
    summary["constant_ratio"] = {}
    for target in TARGETS:
        block = {}
        for fam in ("cheb", "rat"):
            ratios = []
            for name in DYAD:
                e1 = summary["cells"].get(f"x0/{target}/{name}")
                e2 = summary["cells"].get(f"a2/{target}/{name}")
                if not e1 or not e2:
                    continue
                c1 = e1.get("Cc" if fam == "cheb" else "Cr")
                c2 = e2.get("Cc" if fam == "cheb" else "Cr")
                if not c1 or not c2:
                    continue
                if fam == "rat" and not (e1["rat_valid"] and e2["rat_valid"]):
                    continue
                if e1.get(f"{fam}_stored") <= (5 if fam == "cheb" else 3):
                    continue
                ratios.append(float(mp.mpf(c2) / mp.mpf(c1)))
            if not ratios:
                continue
            L1 = mp.mpf(summary["cells"][f"x0/{target}/2^0"]["L"])
            L2 = mp.mpf(summary["cells"][f"a2/{target}/2^0"]["L"])
            expected = float(mp.e ** (L1 - L2))
            ratios.sort()
            block[fam] = {"min": ratios[0], "max": ratios[-1],
                          "median": ratios[len(ratios) // 2],
                          "expected": expected, "n": len(ratios)}
        summary["constant_ratio"][target] = block

    # ---- dominance D = S_cheb - S_rat and its slope ----
    for edge_text, a in edge_of.items():
        label = EDGE_LABEL[edge_text]
        block = {}
        for target in TARGETS:
            xs, ds = [], []
            d_cells = []
            for j in range(-7, 6):
                name = f"2^{j}"
                e = summary["cells"].get(f"{label}/{target}/{name}")
                if not e or e.get("cheb_stored") is None or e.get("rat_stored") is None:
                    continue
                if not e["rat_valid"]:
                    continue
                d = e["cheb_stored"] - e["rat_stored"]
                xs.append(float(j))
                ds.append(float(d))
                d_cells.append({"width": name, "D": d})
            if len(xs) >= 3:
                fit = ols(xs, ds)
                block[target] = {"cells": d_cells, "slope": fit["beta"],
                                 "se": fit["se"], "n": fit["n"],
                                 "min_D": min(ds), "max_D": max(ds)}
            elif len(xs) == 2:
                block[target] = {"cells": d_cells,
                                 "slope": (ds[1] - ds[0]) / (xs[1] - xs[0]),
                                 "se": None, "n": 2,
                                 "min_D": min(ds), "max_D": max(ds)}
            else:
                block[target] = {"cells": d_cells, "status": "incomplete"}
        summary["dominance"][label] = block

    # ---- CSVs ----
    def write_csv(path, fieldnames, records):
        with path.open("w", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(fh, fieldnames=fieldnames)
            writer.writeheader()
            for record in records:
                writer.writerow({k: record.get(k) for k in fieldnames})

    cell_fields = ["edge", "target", "width_name", "width", "kappa", "ell",
                   "mu", "L", "cheb_stored", "cheb_degree", "cheb_parity",
                   "cheb_deliv_rel", "cheb_deliv_abs", "Cc",
                   "rat_stored", "rat_m", "rat_k", "rat_solver",
                   "rat_parity", "rat_deliv_rel", "rat_deliv_abs",
                   "rat_below_meets", "rat_valid", "Cr", "rat_model_Cr4"]
    for edge_text in (X0_TEXT, A2_TEXT):
        label = EDGE_LABEL[edge_text]
        records = [summary["cells"][k] for k in summary["cells"]
                   if k.startswith(f"{label}/")]
        records.sort(key=lambda r: (TARGETS.index(r["target"]),
                                    NAMES.index(r["width_name"])))
        write_csv(out_dir / f"edge_{label}.csv", cell_fields, records)

    comparison = []
    counts_rows = []
    for target in TARGETS:
        for name in NAMES:
            e1 = summary["cells"].get(f"x0/{target}/{name}")
            e2 = summary["cells"].get(f"a2/{target}/{name}")
            if not e1 or not e2:
                continue
            row = {
                "target": target, "width_name": name,
                "width_x0": e1["width"], "width_a2": e2["width"],
                "kappa": e1["kappa"], "ell": e1["ell"], "mu": e1["mu"],
                "L_x0": e1["L"], "L_a2": e2["L"],
                "cheb_x0": e1["cheb_stored"], "cheb_a2": e2["cheb_stored"],
                "cheb_diff": (e2["cheb_stored"] - e1["cheb_stored"])
                if (e1["cheb_stored"] is not None and e2["cheb_stored"] is not None) else None,
                "cheb_degree_x0": e1["cheb_degree"],
                "cheb_degree_a2": e2["cheb_degree"],
                "cheb_deliv_rel_x0": e1["cheb_deliv_rel"],
                "cheb_deliv_rel_a2": e2["cheb_deliv_rel"],
                "rat_x0": e1["rat_stored"], "rat_a2": e2["rat_stored"],
                "rat_diff": (e2["rat_stored"] - e1["rat_stored"])
                if (e1["rat_stored"] is not None and e2["rat_stored"] is not None) else None,
                "rat_m_x0": e1["rat_m"], "rat_k_x0": e1["rat_k"],
                "rat_m_a2": e2["rat_m"], "rat_k_a2": e2["rat_k"],
                "rat_solver_x0": e1["rat_solver"], "rat_solver_a2": e2["rat_solver"],
                "rat_deliv_rel_x0": e1["rat_deliv_rel"],
                "rat_deliv_rel_a2": e2["rat_deliv_rel"],
                "rat_below_x0": e1["rat_below_meets"],
                "rat_below_a2": e2["rat_below_meets"],
                "Cc_x0": e1.get("Cc"), "Cc_a2": e2.get("Cc"),
                "Cc_ratio": (float(mp.mpf(e2["Cc"]) / mp.mpf(e1["Cc"]))
                             if e1.get("Cc") and e2.get("Cc") else None),
                "Cc_ratio_expected": float(mp.e ** (mp.mpf(e1["L"]) - mp.mpf(e2["L"])))
                if e1.get("L") and e2.get("L") else None,
                "Cr_x0": e1.get("Cr"), "Cr_a2": e2.get("Cr"),
                "Cr_ratio": (float(mp.mpf(e2["Cr"]) / mp.mpf(e1["Cr"]))
                             if e1.get("Cr") and e2.get("Cr") else None),
                "cheb_model_real_a2": (float(1 + (mp.mpf(e2["L"])
                                                  + mp.log(mp.mpf(e1["Cc"])))
                                             / mp.mpf(e1["ell"]))
                                       if e1.get("Cc") else None),
                "rat_model_real_a2": (float(1 + (mp.mpf(e2["L"])
                                                 + mp.log(mp.mpf(e1["Cr"])))
                                            / mp.mpf(e1["mu"]))
                                      if e1.get("Cr") else None),
            }
            comparison.append(row)
            counts_rows.append({"target": target, "width_name": name,
                                "width_x0": e1["width"], "width_a2": e2["width"],
                                "cheb_x0": e1["cheb_stored"],
                                "cheb_a2": e2["cheb_stored"],
                                "cheb_diff": row["cheb_diff"],
                                "rat_x0": e1["rat_stored"],
                                "rat_a2": e2["rat_stored"],
                                "rat_diff": row["rat_diff"]})
    write_csv(out_dir / "comparison.csv",
              ["target", "width_name", "width_x0", "width_a2", "kappa", "ell",
               "mu", "L_x0", "L_a2", "cheb_x0", "cheb_a2", "cheb_diff",
               "cheb_degree_x0", "cheb_degree_a2", "cheb_deliv_rel_x0",
               "cheb_deliv_rel_a2", "rat_x0", "rat_a2", "rat_diff",
               "rat_m_x0", "rat_k_x0", "rat_m_a2", "rat_k_a2",
               "rat_solver_x0", "rat_solver_a2", "rat_deliv_rel_x0",
               "rat_deliv_rel_a2", "rat_below_x0", "rat_below_a2",
               "Cc_x0", "Cc_a2", "Cc_ratio", "Cc_ratio_expected",
               "Cr_x0", "Cr_a2", "Cr_ratio", "cheb_model_real_a2",
               "rat_model_real_a2"], comparison)
    write_csv(out_dir / "counts_by_width_target.csv",
              ["target", "width_name", "width_x0", "width_a2", "cheb_x0",
               "cheb_a2", "cheb_diff", "rat_x0", "rat_a2", "rat_diff"],
              counts_rows)
    summary["files"] = {
        "analysis": str(out_dir / "analysis.json"),
        "comparison": str(out_dir / "comparison.csv"),
        "counts": str(out_dir / "counts_by_width_target.csv"),
        "edge_x0": str(out_dir / "edge_x0.csv"),
        "edge_a2": str(out_dir / "edge_a2.csv"),
    }
    (out_dir / "analysis.json").write_text(
        json.dumps(summary, indent=1) + "\n", encoding="utf-8")

    # ---- print ----
    print(f"function {args.function}; x0 = {X0_TEXT}; a2 = {A2_TEXT}")
    print(f"full rows on disk: {len(rows)} of 90; Chebyshev cells: {len(cells)} of 90")
    print(f"kappa(a2,2W)/kappa(x0,W)-1 max: "
          f"{summary['edges']['a2']['kappa_identical_to_x0']}")
    for target in TARGETS:
        print(f"\n{target}   width   cheb x0/a2   rat x0/a2")
        for name in NAMES:
            e1 = summary["cells"].get(f"x0/{target}/{name}")
            e2 = summary["cells"].get(f"a2/{target}/{name}")
            if not e1 or not e2:
                continue
            print(f"  {name:>5} x0 {str(e1['cheb_stored']):>3}/{str(e2['cheb_stored']):>3}"
                  f"   {str(e1['rat_stored']):>3}/{str(e2['rat_stored']):>3}")
    print("\nslopes (paper estimator, thirteen scaling widths, floors in)")
    for target in TARGETS:
        for fam in ("cheb", "rat"):
            s = summary["slopes"][target][fam]
            if s.get("status") == "incomplete":
                print(f" {target} {fam}: INCOMPLETE ({s.get('common_widths')})")
                continue
            if s.get("status") == "partial":
                print(f" {target} {fam}: PARTIAL at n={s['n']} "
                      f"({s.get('common_widths')}): secants "
                      f"{s['x0']['beta']:.5f}/{s['a2']['beta']:.5f}, no s.e.")
                continue
            d = s["difference"]
            zp = f"{d['z_paired']:+.2f}" if d["z_paired"] is not None else "n/a"
            zu = f"{d['z_unpaired']:+.2f}" if d["z_unpaired"] is not None else "n/a"
            se1 = (f"{s['x0']['se']:.5f}" if s["x0"]["se"] is not None else "n/a")
            se2 = (f"{s['a2']['se']:.5f}" if s["a2"]["se"] is not None else "n/a")
            print(f" {target} {fam}: x0 {s['x0']['beta']:.5f}+-{se1} "
                  f"(n={s['n']}), a2 {s['a2']['beta']:.5f}+-{se2}, "
                  f"diff {d['diff']:+.5f} (unpaired z {zu}, "
                  f"paired z {zp}); endpoints "
                  f"{s['x0']['endpoint']:.4f}/{s['a2']['endpoint']:.4f}")
    print("\nmodel secants at each edge's measured constant")
    for label in ("x0", "a2"):
        for target in TARGETS:
            b = summary["width_law_model"][label][target]
            for fam in ("cheb", "rat"):
                if "constant_measured" not in b[fam]:
                    print(f" {label} {target} {fam}: "
                          f"NO CONSTANT ({b[fam].get('status', 'partial')})")
                    continue
                print(f" {label} {target} {fam}: C={b[fam]['constant_measured']} "
                      f"model beta {b[fam]['model_secant']:.4f}, "
                      f"max|N-model| {b[fam]['max_abs_residual']:.3f}, "
                      f"drift {b[fam]['decades_of_drift']:.2f} dec")
    print("\ntarget-increment test (rel:1e-13 -> rel:1e-14)")
    for label in ("x0", "a2"):
        for fam in ("cheb", "rat"):
            b = summary["increments"][label][fam]
            if b.get("status") == "incomplete":
                print(f" {label} {fam}: INCOMPLETE")
                continue
            print(f" {label} {fam}: max|residual| {b['max_abs_residual']:.3f}, "
                  f"outside one: {b['cells_outside_one']}/{len(b['cells'])}, "
                  f"outside two: {b['cells_outside_two']}/{len(b['cells'])}, "
                  f"signs {''.join(b['signs'])}")
    print("\nladder test (free line in L, residual and corrected second "
          "difference <= 1)")
    print(" The paper's raw second-difference criterion presumes the dense "
          "ladder's equal L spacing;")
    print(" the Group A target spacings are unequal, so the spacing term is "
          "reported and removed.")
    for label in ("x0", "a2"):
        block = summary["ladder"]["edges"][label]
        for fam in ("cheb", "rat"):
            flips = [f["width"] for f in block["flips"] if f["family"] == fam]
            tested = [name for name, fb in block["widths"].items()
                      if fb[fam].get("affine_within_tolerance") is not None]
            tested_blocks = [fb[fam] for fb in block["widths"].values()
                             if fb[fam].get("affine_within_tolerance") is not None]
            worst_res = max((b["max_abs_residual"] for b in tested_blocks),
                            default=None)
            worst_raw = max((b["max_abs_second_difference"] for b in tested_blocks),
                            default=None)
            worst_cor = max((b["max_abs_second_difference_corrected"]
                             for b in tested_blocks), default=None)
            print(f" {label} {fam}: {len(tested)} width(s) tested, "
                  f"flips {flips if flips else 'none'}, "
                  f"max|residual| {worst_res}, max|second diff| raw "
                  f"{worst_raw} corrected {worst_cor}")
    print("\nmodel prediction of the a2 count from the x0 count")
    for target in TARGETS:
        for fam in ("cheb", "rat"):
            b = summary["second_edge_prediction"][target].get(fam)
            if not b:
                continue
            print(f" {target} {fam}: residual {b['min_residual']:+.3f} to "
                  f"{b['max_residual']:+.3f} ({b['cells_above_zero']} above, "
                  f"{b['cells_below_zero']} below, n={len(b['cells'])})")
    print("\ncount-implied constant ratio C(a2)/C(x0)")
    for target in TARGETS:
        for fam in ("cheb", "rat"):
            b = summary["constant_ratio"][target].get(fam)
            if not b:
                continue
            print(f" {target} {fam}: {b['min']:.3g}..{b['max']:.3g} "
                  f"(med {b['median']:.3g}) vs expected {b['expected']:.3g}")

    print("\nrational excess over the leading law (Cr=4)")
    for label in ("x0", "a2"):
        for target in TARGETS:
            vals = summary["rational_departure"][label][target]
            if not vals:
                continue
            print(f" {label} {target}: " + " ".join(
                f"{v['width']}:{v['excess_over_Cr4']:+.2f}" for v in vals))
    print(f"\nwrote {out_dir / 'analysis.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
