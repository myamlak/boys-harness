"""De la Vallee Poussin lower bounds from the independent SK fits.

For a fit r of type [m/k], let N = m+k+2.  If the weighted error (P/Q - f)/f
has N alternating extrema of magnitude at least L, then every rational of type
[m/k] has weighted error at least L on the grid (the de la Vallee Poussin
lower bound used by Remez).  This script reads the stored SK fits, refines the
error extrema, selects the best N-point alternating subset and reports L.

Usage:
  python -I lower_bounds.py --json out/p5_s14.json --width-name "2^5" \
      --stored 14 --out out/lb_p5_s14.json
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import mpmath as mp

sys.path.insert(0, str(Path(__file__).resolve().parent))
import remez_rational as R

REC = Path(r"C:\Users\myaml\dev\papers-harness-repo\coefficient-count-sweep"
           r"\clean-edge-2026-11-22\edge62_sweep.json")


def record_width(width_name, target_tag="rel:1e-15"):
    doc = json.loads(REC.read_text(encoding="utf-8"))
    row = next(r for r in doc["rows"] if r["width_name"] == width_name
               and r["target"] == target_tag)
    return float(doc["meta"]["edge"]), float(row["width_f"])


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--json", required=True, help="sk_pairs output")
    ap.add_argument("--width-name", required=True)
    ap.add_argument("--stored", type=int, required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args(argv)
    mp.mp.dps = R.DPS
    a, W = record_width(args.width_name)
    prob = R.problem_f0(a, W)
    ts, tsm, fs = R.build_grid(a, W, 2000)
    ws = [prob.weight(t, f) for t, f in zip(tsm, fs)]
    doc = json.loads(Path(args.json).read_text(encoding="utf-8"))
    rows = [r for r in doc["rows"] if r["stored"] == args.stored]
    out_rows = []
    for row in sorted(rows, key=lambda r: (r["m"], r["k"])):
        m, k = row["m"], row["k"]
        N = m + k + 2
        P = [mp.mpf(c) for c in row["P70"]]
        Q = [mp.mpf(c) for c in row["Q70"]]
        errs, qv = R.profile(P, Q, tsm, fs, ws)
        merged = R.refine_extrema(P, Q, ts, errs, prob)
        window = R.choose_window(merged, N)
        gm = max(abs(e) for e in errs)
        g = R.guard_ratio(qv)
        L = (min(abs(R.error_at(P, Q, t, prob)[0]) for t in window)
             if window else None)

        # grid-node bound: the same selection restricted to 2,001-grid nodes,
        # so it bounds the grid maximum of any competitor of the same type
        cand = [0]
        nlast = len(errs) - 1
        for i in range(1, nlast):
            a_, b_, c_ = abs(errs[i - 1]), abs(errs[i]), abs(errs[i + 1])
            if b_ >= a_ and b_ >= c_ and (b_ > a_ or b_ > c_):
                cand.append(i)
        cand.append(nlast)
        merged_grid = []
        for idx in cand:
            e = errs[idx]
            if e == 0:
                continue
            s = 1 if e > 0 else -1
            if merged_grid and (1 if merged_grid[-1][1] > 0 else -1) == s:
                if abs(e) > abs(merged_grid[-1][1]):
                    merged_grid[-1] = (tsm[idx], e)
            else:
                merged_grid.append((tsm[idx], e))
        window_grid = R.choose_window(merged_grid, N)
        L_grid = (min(abs(R.error_at(P, Q, t, prob)[0])
                      for t in window_grid) if window_grid else None)
        entry = {
            "m": m, "k": k, "N": N,
            "extrema": len(merged),
            "max_abs_error_70": mp.nstr(gm, 17),
            "guard_70": mp.nstr(g, 17),
            "lower_bound_L": (mp.nstr(L, 17) if L is not None else None),
            "L_over_target": (mp.nstr(L / mp.mpf("1e-15"), 8)
                              if L is not None else None),
            "bound_above_target": (bool(L is not None
                                        and L > mp.mpf("1e-15"))),
            "lower_bound_grid_L": (mp.nstr(L_grid, 17)
                                   if L_grid is not None else None),
            "grid_L_over_target": (mp.nstr(L_grid / mp.mpf("1e-15"), 8)
                                   if L_grid is not None else None),
            "grid_bound_above_target": (bool(L_grid is not None
                                             and L_grid > mp.mpf("1e-15"))),
            "sk_e60": row.get("e60"),
            "sk_guard60": row.get("guard60"),
        }
        out_rows.append(entry)
        print(f"[{m:2d}/{k:2d}] N={N} extrema={len(merged):2d} "
              f"L={entry['lower_bound_L']} L/target="
              f"{entry['L_over_target']} "
              f"Lgrid={entry['lower_bound_grid_L']} "
              f"Lgrid/target={entry['grid_L_over_target']} "
              f"max={entry['max_abs_error_70']} guard={entry['guard_70']}")
    min_cont = min((mp.mpf(r["lower_bound_L"]) for r in out_rows
                    if r["lower_bound_L"] is not None), default=None)
    min_grid = min((mp.mpf(r["lower_bound_grid_L"]) for r in out_rows
                    if r["lower_bound_grid_L"] is not None), default=None)
    out = {"schema": "clean-edge-de-la-vallee-poussin-1",
           "width_name": args.width_name, "stored": args.stored,
           "a": a, "width": W, "target": 1e-15,
           "rows": out_rows,
           "all_bounds_above_target": all(r["bound_above_target"]
                                          for r in out_rows),
           "all_grid_bounds_above_target": all(
               r["grid_bound_above_target"] for r in out_rows),
           "min_L": (mp.nstr(min_cont, 17)
                     if min_cont is not None else None),
           "min_grid_L": (mp.nstr(min_grid, 17)
                          if min_grid is not None else None)}
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(json.dumps(out, indent=1, default=str),
                              encoding="utf-8")
    print(f"all bounds above target: {out['all_bounds_above_target']}; "
          f"min L={out['min_L']}")
    print(f"all grid-node bounds above target: "
          f"{out['all_grid_bounds_above_target']}; "
          f"min Lgrid={out['min_grid_L']}")
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
