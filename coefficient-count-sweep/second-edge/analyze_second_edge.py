#!/usr/bin/env python3
"""Compare Paper C's Group A counts at x0 and at the second left edge a2.

The x0 side is the paper's current Group A table: the corrected split-Clenshaw
Chebyshev counts from the published re-count record, and the rational counts
from the core and sensitivity records (the rational solvers and their
evaluator are unchanged between v5.0.0 and 5.1.1; the generator diff is the
odd-degree Clenshaw fix alone).

The a2 side is results-second-edge/sweep_all.json, written by
second_edge_sweep.py with the paper's own count_sweep.py machinery.

The script writes
  results-second-edge/counts_by_width_target.csv   the two edges' counts
  results-second-edge/comparison.csv               every compared field
  results-second-edge/analysis.json                the summary and the slopes
and prints the same tables.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import mpmath as mp

PAPER = Path(__file__).resolve().parent
NAMES = [f"2^{j}" for j in range(-7, 6)] + ["ship", "paper"]
DYAD = [f"2^{j}" for j in range(-7, 6)]
TARGETS = ["rel:1e-13", "abs:5e-14", "rel:1e-14"]
X0_TEXT = "11.899848152108484"
mp.mp.dps = 50
X0 = mp.mpf(X0_TEXT)

RECOUNT = PAPER.parent / "records" / "chebyshev_recount.json"
RESULT_FILES = (PAPER.parent / "records" / "core" / "sweep_all.json",
                PAPER.parent / "records" / "core" / "sweep_rel_1e-14.json")


def width_of(name: str):
    if name.startswith("2^"):
        return mp.mpf(2) ** int(name.split("^")[1])
    if name == "ship":
        return mp.mpf("31") - X0
    if name == "paper":
        return mp.mpf("28.989337738820740") - X0
    raise ValueError(name)


def label_for(target: str) -> str:
    kind, value = target.split(":", 1)
    return f"{kind}:{value}"


def load_baseline(path: Path, result_files=RESULT_FILES):
    """The paper's x0 Group A counts, keyed by (target, width name)."""
    if not path.exists():
        raise SystemExit(f"baseline record missing: {path}")
    doc = json.loads(path.read_text(encoding="utf-8"))
    cells = {}
    for cell in doc["cells"]:
        if cell.get("mode") != "records" or cell.get("nodes") != 2000:
            continue
        if cell["target"] not in TARGETS:
            continue
        names = [n for n in cell.get("width_names", []) if n in NAMES]
        if not names:
            continue
        name = names[0]
        cells[(cell["target"], name)] = {
            "stored": cell["new_stored"],
            "degree": cell["new_m"],
            "delivered_rel": cell["delivered_rel"],
            "delivered_abs": cell["delivered_abs"],
            "at_x_rel": cell["at_x_rel"],
            "at_x_abs": cell["at_x_abs"],
        }
    rational = {}
    for path in result_files:
        doc = json.loads(path.read_text(encoding="utf-8"))
        for row in doc["rows"]:
            if row["width_name"] not in NAMES or row["target"] not in TARGETS:
                continue
            rational[(row["target"], row["width_name"])] = {
                "stored": row["rat"].get("stored"),
                "m": row["rat"].get("m"),
                "k": row["rat"].get("k"),
                "solver": row["rat"].get("solver"),
                "delivered_rel": row["rat"].get("delivered_rel"),
                "delivered_abs": row["rat"].get("delivered_abs"),
                "at_x_rel": row["rat"].get("at_x_rel"),
                "at_x_abs": row["rat"].get("at_x_abs"),
                "below_meets": row["rat"].get("below_meets"),
            }
    missing_cheb = [(t, n) for t in TARGETS for n in NAMES
                    if (t, n) not in cells]
    missing_rat = [(t, n) for t in TARGETS for n in NAMES
                   if (t, n) not in rational]
    if missing_cheb or missing_rat:
        raise SystemExit(f"baseline incomplete: cheb {missing_cheb} "
                         f"rat {missing_rat}")
    return cells, rational


def load_second(path: Path):
    if not path.exists():
        raise SystemExit(f"second-edge record missing: {path}")
    doc = json.loads(path.read_text(encoding="utf-8"))
    rows = {(row["target"], row["width_name"]): row for row in doc["rows"]}
    return doc, rows


def kappa(a, width):
    s = mp.sqrt(a)
    b = mp.sqrt(a + width)
    return (b - s) / (b + s)


def mu_of(kappa_value):
    return (mp.pi / 2) * mp.ellipk(1 - kappa_value ** 2) / mp.ellipk(kappa_value ** 2)


def f0(a):
    return mp.sqrt(mp.pi) * mp.erf(mp.sqrt(a)) / (2 * mp.sqrt(a))


def target_L(target: str, a):
    kind, value = target.split(":", 1)
    eps = mp.mpf(value)
    if kind == "rel":
        return mp.log(1 / eps), "rel", eps
    return mp.log(f0(a) / eps), "abs", eps


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
    """Slope difference at the two edges, unpaired and paired standard errors."""
    diff = fit2["beta"] - fit1["beta"]
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


def fmt(value, digits=6):
    if value is None:
        return "--"
    if isinstance(value, (int, str)):
        return str(value)
    return f"{value:.{digits}g}"


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--second-edge",
                        default=str(PAPER / "sweep_all.json"))
    parser.add_argument("--chebyshev-only",
                        default=str(PAPER / "chebyshev_only.json"))
    parser.add_argument("--recount", default=str(RECOUNT))
    parser.add_argument("--baseline", action="append", default=None,
                        help="a sweep record holding the first edge's counts; repeatable. "
                             "Default: the deposit's records/core/.")
    parser.add_argument("--out-dir",
                        default=str(PAPER))
    args = parser.parse_args(argv)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    mp.mp.dps = 50
    cheb0, rat0 = load_baseline(
        Path(args.recount),
        [Path(p) for p in args.baseline] if args.baseline else RESULT_FILES)
    second_doc, second = load_second(Path(args.second_edge))
    cheb_only = {}
    cheb_only_path = Path(args.chebyshev_only)
    if cheb_only_path.exists():
        only_doc = json.loads(cheb_only_path.read_text(encoding="utf-8"))
        if only_doc.get("meta", {}).get("left_edge") != mp.nstr(X0 * 2, 17):
            raise SystemExit(f"{cheb_only_path} is not at the second left edge")
        cheb_only = {(c["target"], c["width_name"]): c
                     for c in only_doc.get("cells", [])}

    a2 = X0 * 2
    rows_out = []
    counts_rows = []
    summary = {"x0": X0_TEXT, "a2": mp.nstr(a2, 17),
               "targets": {}, "cells": {}, "slopes": {}, "C_c": {},
               "presence": {}}

    for target in TARGETS:
        L1, kind, eps = target_L(target, X0)
        L2, _, _ = target_L(target, a2)
        for name in NAMES:
            width1 = width_of(name)
            width2 = width1 * 2
            k1 = kappa(X0, width1)
            k2 = kappa(a2, width2)
            ell = mp.log(1 / k1)
            mu = mu_of(k1)
            entry = {
                "width": mp.nstr(width1, 17),
                "width_a2": mp.nstr(width2, 17),
                "kappa_x0": mp.nstr(k1, 17),
                "kappa_a2": mp.nstr(k2, 17),
                "kappa_rel_diff": mp.nstr(abs(k2 / k1 - 1), 6),
                "ell": mp.nstr(ell, 17),
                "mu": mp.nstr(mu, 17),
                "L_x0": mp.nstr(L1, 17),
                "L_a2": mp.nstr(L2, 17),
                "F0_x0": mp.nstr(f0(X0), 17),
                "F0_a2": mp.nstr(f0(a2), 17),
            }
            present = (target, name) in second
            only = cheb_only.get((target, name))
            summary["presence"][f"{target}/{name}"] = {
                "full_row": present, "cheb_only": only is not None}
            if present:
                row2 = second[(target, name)]
                cheb2 = row2["cheb"]
                rat2 = row2["rat"]
                entry["a2_width"] = row2["width"]
                entry["a2_b"] = row2["b"]
                entry["cheb_source"] = "sweep"
            elif only is not None:
                cheb2 = only["cheb"]
                rat2 = {}
                entry["a2_width"] = only["width"]
                entry["a2_b"] = only["b"]
                entry["cheb_source"] = "chebyshev_only"
            else:
                cheb2 = rat2 = {}
                entry["a2_width"] = None
                entry["a2_b"] = None
                entry["cheb_source"] = None
            c1 = cheb0[(target, name)]
            r1 = rat0[(target, name)]
            c2s = cheb2.get("stored")
            r2s = rat2.get("stored")
            entry.update({
                "cheb_x0": c1["stored"], "cheb_a2": c2s,
                "cheb_diff": None if c2s is None else c2s - c1["stored"],
                "cheb_degree_x0": c1["degree"],
                "cheb_degree_a2": cheb2.get("m"),
                "cheb_parity_x0": "odd" if c1["stored"] % 2 else "even",
                "cheb_parity_a2": (None if c2s is None
                                   else ("odd" if c2s % 2 else "even")),
                "cheb_deliv_rel_x0": c1["delivered_rel"],
                "cheb_deliv_rel_a2": cheb2.get("delivered_rel"),
                "cheb_deliv_abs_x0": c1["delivered_abs"],
                "cheb_deliv_abs_a2": cheb2.get("delivered_abs"),
                "rat_x0": r1["stored"], "rat_a2": r2s,
                "rat_diff": None if r2s is None else r2s - r1["stored"],
                "rat_m_x0": r1["m"], "rat_k_x0": r1["k"],
                "rat_m_a2": rat2.get("m"), "rat_k_a2": rat2.get("k"),
                "rat_solver_x0": r1["solver"], "rat_solver_a2": rat2.get("solver"),
                "rat_deliv_rel_x0": r1["delivered_rel"],
                "rat_deliv_rel_a2": rat2.get("delivered_rel"),
                "rat_deliv_abs_x0": r1["delivered_abs"],
                "rat_deliv_abs_a2": rat2.get("delivered_abs"),
                "rat_below_x0": r1["below_meets"],
                "rat_below_a2": rat2.get("below_meets"),
                "rat_parity_x0": ("odd" if r1["stored"] % 2 else "even")
                if r1["stored"] is not None else None,
                "rat_parity_a2": ("odd" if r2s % 2 else "even")
                if r2s is not None else None,
            })
            # Count-implied constants, at each edge.
            if c1["stored"] is not None:
                entry["Cc_x0"] = mp.nstr(mp.exp((c1["stored"] - 1) * ell - L1), 17)
            if c2s is not None:
                entry["Cc_a2"] = mp.nstr(mp.exp((c2s - 1) * ell - L2), 17)
            if entry.get("Cc_x0") and entry.get("Cc_a2"):
                entry["Cc_ratio"] = mp.nstr(mp.mpf(entry["Cc_a2"])
                                            / mp.mpf(entry["Cc_x0"]), 17)
                entry["Cc_ratio_expected"] = mp.nstr(mp.exp(L1 - L2), 17)
            if r1["stored"] is not None:
                entry["Cr_x0"] = mp.nstr(mp.exp((r1["stored"] - 1) * mu - L1), 17)
            if r2s is not None:
                entry["Cr_a2"] = mp.nstr(mp.exp((r2s - 1) * mu - L2), 17)
            if entry.get("Cr_x0") and entry.get("Cr_a2"):
                entry["Cr_ratio"] = mp.nstr(mp.mpf(entry["Cr_a2"])
                                            / mp.mpf(entry["Cr_x0"]), 17)
            # The model's own prediction for the count at a2, with the
            # constant measured at x0 and the cross-ratio carried over.
            if c1["stored"] is not None:
                cc1 = mp.exp((c1["stored"] - 1) * ell - L1)
                real = 1 + (L2 + mp.log(cc1)) / ell
                entry["cheb_model_real_a2"] = mp.nstr(real, 17)
                entry["cheb_model_ceil_a2"] = int(mp.ceil(real))
            if r1["stored"] is not None:
                cr1 = mp.exp((r1["stored"] - 1) * mu - L1)
                real = 1 + (L2 + mp.log(cr1)) / mu
                entry["rat_model_real_a2"] = mp.nstr(real, 17)
                entry["rat_model_ceil_a2"] = int(mp.ceil(real))
            summary["cells"][f"{target}/{name}"] = entry
            rows_out.append((target, name, entry))
            counts_rows.append({
                "target": target,
                "width_name": name,
                "width_x0": mp.nstr(width1, 17),
                "width_a2": mp.nstr(width2, 17),
                "cheb_x0": c1["stored"], "cheb_a2": c2s,
                "cheb_diff": entry["cheb_diff"],
                "rat_x0": r1["stored"], "rat_a2": r2s,
                "rat_diff": entry["rat_diff"],
            })

    # Per family and target summaries over the valid cells.
    for target in TARGETS:
        block = {"cheb": {}, "rat": {}}
        for fam in ("cheb", "rat"):
            diffs = []
            parity_flips = 0
            pair_changes = 0
            max_abs = 0
            for name in NAMES:
                e = summary["cells"][f"{target}/{name}"]
                n1 = e[f"{fam}_x0"]
                n2 = e[f"{fam}_a2"]
                if n1 is None or n2 is None:
                    continue
                if fam == "rat" and (e["rat_below_x0"] or e["rat_below_a2"]):
                    continue
                d = n2 - n1
                diffs.append({"width": name, "x0": n1, "a2": n2, "diff": d})
                max_abs = max(max_abs, abs(d))
                if e[f"{fam}_parity_x0"] != e[f"{fam}_parity_a2"]:
                    parity_flips += 1
                if fam == "cheb":
                    if (e["cheb_degree_x0"] != e["cheb_degree_a2"]):
                        pair_changes += 1
                else:
                    if (e["rat_m_x0"], e["rat_k_x0"]) != (e["rat_m_a2"], e["rat_k_a2"]):
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
        summary["targets"][target] = block

    # Width-law slopes over the thirteen scaling widths, floor rows included,
    # on the paper's own estimator: least squares of log2 N on log2 W.
    def slope_block(get1, get2, valid, js):
        xs = [float(j) for j in js]
        ys1 = [math.log2(get1(f"2^{j}")) for j in js]
        ys2 = [math.log2(get2(f"2^{j}")) for j in js]
        fit1 = ols(xs, ys1)
        fit2 = ols(xs, ys2)
        first, last = 0, len(xs) - 1
        endpoint1 = (ys1[last] - ys1[first]) / (xs[last] - xs[first])
        endpoint2 = (ys2[last] - ys2[first]) / (xs[last] - xs[first])
        return {"x0": {"alpha": fit1["alpha"], "beta": fit1["beta"],
                       "se": fit1["se"], "n": fit1["n"],
                       "endpoint": endpoint1},
                "a2": {"alpha": fit2["alpha"], "beta": fit2["beta"],
                       "se": fit2["se"], "n": fit2["n"],
                       "endpoint": endpoint2},
                "difference": paired_slope_test(fit1, fit2)}

    for target in TARGETS:
        summary["slopes"][target] = {}
        for fam in ("cheb", "rat"):
            def get1(n, fam=fam, target=target):
                return summary["cells"][f"{target}/{n}"][f"{fam}_x0"]

            def get2(n, fam=fam, target=target):
                return summary["cells"][f"{target}/{n}"][f"{fam}_a2"]

            def valid(n, fam=fam, target=target):
                e = summary["cells"][f"{target}/{n}"]
                if e[f"{fam}_x0"] is None or e[f"{fam}_a2"] is None:
                    return False
                if fam == "rat" and (e["rat_below_x0"] or e["rat_below_a2"]):
                    return False
                return True

            js = [j for j in range(-7, 6) if valid(f"2^{j}")]
            block = {"common_widths": [f"2^{j}" for j in js]}
            if len(js) < 10:
                block["status"] = "incomplete"
                summary["slopes"][target][fam] = block
                continue
            block.update(slope_block(get1, get2, valid, js))
            summary["slopes"][target][fam] = block

    # Count-implied C_c bands, and the C_c ratio against its no-drift value.
    for target in TARGETS:
        band = {}
        for fam, const in (("cheb", "Cc"), ("rat", "Cr")):
            vals = {"x0": [], "a2": []}
            ratios = []
            for name in DYAD:
                e = summary["cells"][f"{target}/{name}"]
                if e.get(f"{const}_x0") and e.get(f"{const}_a2"):
                    vals["x0"].append(float(e[f"{const}_x0"]))
                    vals["a2"].append(float(e[f"{const}_a2"]))
                    ratios.append(float(e[f"{const}_ratio"]))
            if vals["x0"]:
                band[fam] = {
                    "x0_min": min(vals["x0"]), "x0_max": max(vals["x0"]),
                    "x0_median": sorted(vals["x0"])[len(vals["x0"]) // 2],
                    "a2_min": min(vals["a2"]), "a2_max": max(vals["a2"]),
                    "a2_median": sorted(vals["a2"])[len(vals["a2"]) // 2],
                    "ratio_min": min(ratios), "ratio_max": max(ratios),
                    "ratio_median": sorted(ratios)[len(ratios) // 2],
                    "expected_ratio": float(mp.exp(
                        target_L(target, X0)[0] - target_L(target, a2)[0])),
                }
        summary["C_c"][target] = band

    # Write outputs.
    (out_dir / "analysis.json").write_text(
        json.dumps(summary, indent=2, default=str) + "\n", encoding="utf-8")
    counts_fields = ["target", "width_name", "width_x0", "width_a2",
                     "cheb_x0", "cheb_a2", "cheb_diff",
                     "rat_x0", "rat_a2", "rat_diff"]
    with (out_dir / "counts_by_width_target.csv").open("w", newline="",
                                                       encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=counts_fields)
        writer.writeheader()
        writer.writerows(counts_rows)
    comp_fields = ["target", "width_name", "width", "width_a2", "kappa_x0",
                   "kappa_a2", "kappa_rel_diff", "ell", "mu", "L_x0", "L_a2",
                   "cheb_x0", "cheb_a2", "cheb_diff", "cheb_degree_x0",
                   "cheb_degree_a2", "cheb_parity_x0", "cheb_parity_a2",
                   "cheb_deliv_rel_x0", "cheb_deliv_rel_a2",
                   "cheb_deliv_abs_x0", "cheb_deliv_abs_a2",
                   "rat_x0", "rat_a2", "rat_diff", "rat_m_x0", "rat_k_x0",
                   "rat_m_a2", "rat_k_a2", "rat_solver_x0", "rat_solver_a2",
                   "rat_deliv_rel_x0", "rat_deliv_rel_a2",
                   "rat_deliv_abs_x0", "rat_deliv_abs_a2",
                   "rat_below_x0", "rat_below_a2",
                   "Cc_x0", "Cc_a2", "Cc_ratio", "Cc_ratio_expected",
                   "Cr_x0", "Cr_a2",
                   "cheb_model_real_a2", "cheb_model_ceil_a2",
                   "rat_model_real_a2", "rat_model_ceil_a2"]
    with (out_dir / "comparison.csv").open("w", newline="",
                                           encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=comp_fields)
        writer.writeheader()
        for target, name, entry in rows_out:
            record = {key: entry.get(key) for key in comp_fields}
            record["target"] = target
            record["width_name"] = name
            writer.writerow(record)

    # Terminal report.
    print(f"x0 = {X0_TEXT}; a2 = {mp.nstr(a2, 17)}")
    full = sum(1 for v in summary["presence"].values() if v["full_row"])
    cheb_present = sum(1 for v in summary["presence"].values()
                       if v["full_row"] or v["cheb_only"])
    print(f"full rows on disk: {full} of {len(summary['presence'])}; "
          f"Chebyshev cells: {cheb_present} of {len(summary['presence'])}")
    missing_full = [k for k, v in summary["presence"].items()
                    if not v["full_row"]]
    missing_cheb = [k for k, v in summary["presence"].items()
                    if not (v["full_row"] or v["cheb_only"])]
    if missing_full:
        print("missing rational rows:", ", ".join(missing_full))
    if missing_cheb:
        print("missing Chebyshev cells:", ", ".join(missing_cheb))
    print("\ncounts by target and width (cheb x0 -> a2; rat x0 -> a2)")
    for target in TARGETS:
        for name in NAMES:
            e = summary["cells"][f"{target}/{name}"]
            print(f"{target:>10} {name:>5}  cheb {fmt(e['cheb_x0'])}->"
                  f"{fmt(e['cheb_a2'])}  rat {fmt(e['rat_x0'])}->"
                  f"{fmt(e['rat_a2'])}")
    print("\nper family and target")
    for target in TARGETS:
        for fam in ("cheb", "rat"):
            b = summary["targets"][target][fam]
            print(f"{target:>10} {fam:4}: differing {b['differing_cells']}"
                  f"/{b['valid_cells']}, max |diff| {b['max_abs_diff']}, "
                  f"sum {b['sum_diff']}, parity flips {b['parity_flips']}, "
                  f"pair changes {b['pair_changes']}")
    print("\nslopes (paper estimator, thirteen scaling widths, floors in)")
    for target in TARGETS:
        for fam in ("cheb", "rat"):
            s = summary["slopes"][target][fam]
            if s.get("status") == "incomplete":
                print(f"{target:>10} {fam:4}: INCOMPLETE")
                continue
            d = s["difference"]
            z1 = ("--" if d["z_unpaired"] is None
                  else f"{d['z_unpaired']:+.2f}")
            z2 = ("--" if d["z_paired"] is None else f"{d['z_paired']:+.2f}")
            print(f"{target:>10} {fam:4}: x0 {s['x0']['beta']:.3f}"
                  f"+-{s['x0']['se']:.3f} (n={s['x0']['n']}), "
                  f"a2 {s['a2']['beta']:.3f}+-{s['a2']['se']:.3f}, "
                  f"diff {d['diff']:+.3f} "
                  f"(unpaired se {d['se_unpaired']:.3f}, "
                  f"z {z1}; paired se {d['se_paired']:.3f}, "
                  f"z {z2}); "
                  f"endpoints {s['x0']['endpoint']:.4f}/{s['a2']['endpoint']:.4f}")
    print("\ncount-implied constants (dyadic widths)")
    for target in TARGETS:
        for fam, key in (("cheb", "Cc"), ("rat", "Cr")):
            b = summary["C_c"][target].get(fam)
            if not b:
                continue
            print(f"{target:>10} {key}: x0 {b['x0_min']:.4g}"
                  f"..{b['x0_max']:.4g} (med {b['x0_median']:.4g}); "
                  f"a2 {b['a2_min']:.4g}..{b['a2_max']:.4g} "
                  f"(med {b['a2_median']:.4g}); ratio "
                  f"{b['ratio_min']:.4g}..{b['ratio_max']:.4g} "
                  f"(med {b['ratio_median']:.4g}) vs expected "
                  f"{b['expected_ratio']:.4g}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
