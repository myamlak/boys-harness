#!/usr/bin/env python3
"""Build the summary table and the affine-residual effects from the cell records.

Reads cell_*.json in this directory (written by control_60dig_brasil.py) and
prints a Markdown table plus the model-based residual recomputation.  The
number of decades and the model's predicted increments are the paper's
(Table S7 of the supplement: paper 0.796/decade over 7 decades -> 5.57,
ship 0.816 -> 5.71, 2^5 0.916 -> 6.41); the counts here are the 60-digit
counts of the control.

Usage: python make_summary_table.py [--json out.json]
"""

import argparse
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent

DEEP = [("A", "paper", "10^{-15}", 1e-15),
        ("B", "ship", "10^{-15}", 1e-15),
        ("C1", "2^5", "10^{-15}", 1e-15),
        ("C2", "2^5", "10^{-14}", 1e-14),
        ("D", "2^0", "10^{-15}", 1e-15)]
ENDPOINTS = [("E-paper", "paper", "10^{-8}"),
             ("E-ship", "ship", "10^{-8}"),
             ("E-p5", "2^5", "10^{-8}")]
MODEL = {"2^0": (0.435, 3.05), "paper": (0.796, 5.57),
         "ship": (0.816, 5.71), "2^5": (0.916, 6.41)}


def load(key):
    p = HERE / f"cell_{key}.json"
    return json.loads(p.read_text()) if p.exists() else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--json", default=None)
    args = ap.parse_args()
    out = {"cells": [], "residuals": []}
    print("| Cell | Width | Target | Record | Record pair | Delivered (binary64) "
          "| 60-digit count | Control pair | e60 | e64r | e64 | e64/target "
          "| Below-count best e60 | Below/target |")
    print("| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- "
          "| --- | --- | --- |")
    for key, wname, tlabel, target in DEEP:
        d = load(key)
        if d is None or d.get("count60") is None:
            print(f"| {key} | {wname} | {tlabel} | missing | | | | | | | | | | |")
            continue
        rec = d["reproduction"]["record"]
        acc = d["accepted"]
        fail = d.get("failing_level") or {}
        fail_rows = fail.get("rows") or []
        fail_best = min((float(r["e60"]) for r in fail_rows
                         if r.get("e60")
                         and float(r.get("guard60", 0.0)) >= 1e-3),
                        default=None)
        ratio = float(acc["e64"]) / target
        fail_ratio = (fail_best / target) if fail_best is not None else None
        print(f"| {key} | {wname} | ${tlabel}$ | {rec['stored']} "
              f"| [{rec['pair'][0]}/{rec['pair'][1]}] "
              f"| {rec['delivered_2001']!r} | {d['count60']} "
              f"| [{acc['m']}/{acc['k']}] {acc['settings']} "
              f"| {acc['e60']} | {acc['e64r']} | {acc['e64']} "
              f"| {ratio:.3g} "
              f"| {('%.6g' % fail_best) if fail_best is not None else '-'} "
              f"| {('%.3g' % fail_ratio) if fail_ratio is not None else '-'} |")
        out["cells"].append({"key": key, "width": wname, "target": tlabel,
                             "record": rec, "count60": d["count60"],
                             "below_best_e60": fail_best,
                             "below_ratio": fail_ratio,
                             "accepted": {k: acc.get(k) for k in
                                          ("m", "k", "settings", "e60",
                                           "e64r", "e64", "meets60")}})
    print()
    counts = {key: (load(key) or {}).get("count60") for key, _, _, _ in DEEP}
    ends = {key: (load(key) or {}).get("count60") for key, _, _ in ENDPOINTS}
    for key, wname, _ in ENDPOINTS:
        d = load(key)
        if d:
            rec = d["reproduction"]["record"]
            print(f"endpoint {key}: record {rec['stored']} -> 60-digit "
                  f"{d['count60']} (e60 {d['accepted']['e60']}, "
                  f"e64 {d['accepted']['e64']})")
    print()
    pairs = [("paper", "A", "E-paper"), ("ship", "B", "E-ship"),
             ("2^5", "C1", "E-p5")]
    for wname, deep, endpoint in pairs:
        if counts.get(deep) is None or ends.get(endpoint) is None:
            continue
        per_decade, predicted = MODEL[wname]
        increment = counts[deep] - ends[endpoint]
        residual = increment - predicted
        residual_exact = increment - per_decade * 7
        row = {"width": wname, "deep_count": counts[deep],
               "endpoint_count": ends[endpoint], "increment": increment,
               "predicted": predicted, "residual": round(residual, 2),
               "per_decade": per_decade,
               "residual_exact": round(residual_exact, 3)}
        out["residuals"].append(row)
        print(f"{wname}: 60-digit increment {increment} "
              f"(deep {counts[deep]} - endpoint {ends[endpoint]}) "
              f"against predicted {predicted}: residual "
              f"{residual:+.2f} (exact {residual_exact:+.3f}); "
              f"record residual {11 - predicted if wname != '2^5' else 17 - predicted:+.2f}")
    if args.json:
        Path(args.json).write_text(json.dumps(out, indent=1))
        print(f"\nwrote {args.json}")


if __name__ == "__main__":
    main()
