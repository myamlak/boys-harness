"""Polish a stored independent SK fit with the Remez/Newton exchange.

Reads one pair's P70/Q70 from an `sk_pairs.py` output, rebuilds the clean-edge
problem, and runs the Remez exchange from that fit (Newton first, then DC if
Newton fails), with both global and Maehly exchange.  Reports the best
admissible fit, the equioscillation certificate and the 60-digit score.

Usage:
  python -I polish_pair.py --json out/ship_s13.json --m 6 --k 6 \
      --width-name ship --out out/polish_ship_66.json
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import mpmath as mp

sys.path.insert(0, str(Path(__file__).resolve().parent))
import remez_rational as R

REC = Path(r"C:\Users\myaml\dev\papers-harness-repo\coefficient-count-sweep"
           r"\clean-edge-2026-11-22\edge62_sweep.json")


def load_problem(width_name, target_tag="rel:1e-15"):
    doc = json.loads(REC.read_text(encoding="utf-8"))
    row = next(r for r in doc["rows"] if r["width_name"] == width_name
               and r["target"] == target_tag)
    a = float(doc["meta"]["edge"])
    W = float(row["width_f"])
    return a, W, float(row["target_value"])


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--json", required=True)
    ap.add_argument("--m", type=int, required=True)
    ap.add_argument("--k", type=int, required=True)
    ap.add_argument("--width-name", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args(argv)
    mp.mp.dps = R.DPS
    a, W, target = load_problem(args.width_name)
    prob = R.problem_f0(a, W)
    ts, tsm, fs = R.build_grid(a, W, 2000)
    ts4, tsm4, fs4 = R.build_grid(a, W, 4000)
    src = json.loads(Path(args.json).read_text(encoding="utf-8"))
    row = next(r for r in src["rows"]
               if r["m"] == args.m and r["k"] == args.k)
    P = [mp.mpf(c) for c in row["P70"]]
    Q = [mp.mpf(c) for c in row["Q70"]]
    ws = [prob.weight(t, f) for t, f in zip(tsm, fs)]
    N = args.m + args.k + 2

    errs, qv = R.profile(P, Q, tsm, fs, ws)
    gm0 = max(abs(e) for e in errs)
    g0 = R.guard_ratio(qv)
    merged = R.refine_extrema(P, Q, ts, errs, prob)
    print(f"start [{args.m}/{args.k}] SK e={mp.nstr(gm0, 10)} "
          f"guard={mp.nstr(g0, 6)} extrema={len(merged)} N={N}")
    for t, e in merged:
        print(f"    t={mp.nstr(t, 12):>18s} e={mp.nstr(e, 8)}")

    ref0 = R.choose_window(merged, N)
    attempts = []
    for mode in ("global", "maehly"):
        for phase, dc_iters in (("newton", 0), ("dc", 120)):
            t0 = time.perf_counter()
            r = R.solve_pair(args.m, args.k, tsm, fs, prob, target,
                             dc_iters=dc_iters, max_outer=100,
                             exchange_mode=mode,
                             ref0=ref0, init=(P, Q, mp.mpf(0)),
                             verbose=args.verbose)
            r["attempt"] = {"exchange": mode, "phase": phase,
                            "dc_iters": dc_iters}
            r["wall"] = round(time.perf_counter() - t0, 1)
            attempts.append(r)
            print(f"attempt {mode:6s} {phase:6s}: conv={r['converged']} "
                  f"grid_max={r['grid_max']} guard={r['guard']} "
                  f"ratio={r['ratio']} newton={r['newton_note']} "
                  f"({r['wall']}s)", flush=True)
            if r["converged"]:
                break
        if any(x["converged"] for x in attempts):
            break

    # choose the best by 60-digit e60 among admissible fits
    scored = []
    for r in attempts:
        if "P" not in r:
            continue
        Pb = [mp.mpf(c) for c in r["P"]]
        Qb = [mp.mpf(c) for c in r["Q"]]
        s = R.score60(Pb, Qb, a, W, target, grid=(ts, tsm, fs))
        s4 = R.score60(Pb, Qb, a, W, target, grid=(ts4, tsm4, fs4))
        scored.append((mp.mpf(s["e60"]), r, Pb, Qb, s, s4))
    scored.sort(key=lambda x: x[0])
    best_e, r, Pb, Qb, s, s4 = scored[0]
    out = {
        "schema": "clean-edge-independent-polish-1",
        "width_name": args.width_name, "width": W, "a": a,
        "m": args.m, "k": args.k, "stored": args.m + args.k + 1,
        "target": target,
        "sk_start": {"e60": row.get("e60"), "guard60": row.get("guard60"),
                     "extrema": len(merged)},
        "attempts": [{k: r[k] for k in
                      ("converged", "status", "grid_max", "guard", "ratio",
                       "signs_ok", "E", "ref_min", "iterations", "seconds",
                       "newton_note", "meets_raw", "attempt", "wall")}
                     for r in attempts],
        "best": {
            "exchange": r["attempt"], "converged": r["converged"],
            "status": r["status"], "grid_max": r["grid_max"],
            "guard": r["guard"], "ratio": r["ratio"],
            "e60": s["e60"], "e60_abs": s["e60_abs"],
            "e60_at_t": s["e60_at_t"], "guard60": s["guard60"],
            "meets60": s["meets60"], "e60_4001": s4["e60"],
            "guard60_4001": s4["guard60"], "meets60_4001": s4["meets60"],
            "P70": [mp.nstr(c, 70) for c in Pb],
            "Q70": [mp.nstr(c, 70) for c in Qb],
        },
    }
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(json.dumps(out, indent=1, default=str),
                              encoding="utf-8")
    b = out["best"]
    print(f"\nBEST [{args.m}/{args.k}]: conv={b['converged']} "
          f"e60={b['e60']} guard60={b['guard60']} meets={b['meets60']} "
          f"e60_4001={b['e60_4001']} meets4001={b['meets60_4001']}")
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
