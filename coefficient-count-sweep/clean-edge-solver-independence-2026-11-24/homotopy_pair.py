"""Function homotopy for one clean-edge pair.

f_0(t) = A/(B+t) is a rational function contained exactly in [m/k] with a
positive linear denominator.  f_lambda = (1-lambda) f_0 + lambda f is
continued from lambda = 0.1 to 1.0, warm-starting the Remez/Newton solver from
the previous minimax.  A good-basin denominator (B chosen from the target
solution's linear denominator coefficient) keeps the path away from the
spurious near-pole minima.

Usage:
  python -I homotopy_pair.py --width-name ship --m 6 --k 6 --B 2.34 \
      --out out/homotopy_ship_66.json
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


def load_width(width_name, target_tag="rel:1e-15"):
    doc = json.loads(REC.read_text(encoding="utf-8"))
    row = next(r for r in doc["rows"] if r["width_name"] == width_name
               and r["target"] == target_tag)
    return (float(doc["meta"]["edge"]), float(row["width_f"]),
            float(row["target_value"]))


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--width-name", required=True)
    ap.add_argument("--m", type=int, required=True)
    ap.add_argument("--k", type=int, required=True)
    ap.add_argument("--B", default="1.8")
    ap.add_argument("--out", required=True)
    args = ap.parse_args(argv)
    mp.mp.dps = R.DPS
    a, W, target = load_width(args.width_name)
    m, k = args.m, args.k
    N = m + k + 2
    prob_f = R.problem_f0(a, W)
    ts, tsm, fs = R.build_grid(a, W, 2000)
    ts4, tsm4, fs4 = R.build_grid(a, W, 4000)
    f_mid = prob_f.f(mp.mpf(0))
    B = mp.mpf(args.B)
    A0 = f_mid * B

    def f0_eval(t):
        return A0 / (B + t)

    def make_problem(lam):
        def eval_f(t):
            return (1 - lam) * f0_eval(t) + lam * prob_f.f(t)
        return R.Problem(eval_f, lambda t, f: 1 / f)

    def interp_init(prob):
        nodes = [mp.cos(mp.pi * (2 * i + 1) / (2 * (N - 1)))
                 for i in range(N - 1)]
        A_ = []
        b_ = []
        for t in nodes:
            f = prob.f(t)
            row = []
            tp = mp.mpf(1)
            for _ in range(m + 1):
                row.append(tp)
                tp *= t
            tp = mp.mpf(1)
            for _ in range(k):
                tp *= t
                row.append(-f * tp)
            A_.append(row)
            b_.append(f)
        u = mp.lu_solve(mp.matrix(A_), mp.matrix(b_))
        P = [u[j] for j in range(m + 1)]
        Q = [mp.mpf(1)] + [u[m + 1 + j] for j in range(k)]
        return P, Q, mp.mpf(0)

    lams = [mp.mpf(s) for s in
            ("0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8",
             "0.9", "0.95", "1.0")]
    P = Q = E = None
    ref = None
    steps = []
    t_start = time.perf_counter()
    for lam in lams:
        prob = make_problem(lam)
        fs_lam = [prob.f(t) for t in tsm]
        if P is None:
            init = interp_init(prob)
            P, Q, E = init
            ref0 = None
        else:
            init = (P, Q, E)
            ref0 = ref
        t0 = time.perf_counter()
        r = R.solve_pair(m, k, tsm, fs_lam, prob, target,
                         dc_iters=0, max_outer=80, exchange_mode="global",
                         ref0=ref0, init=init)
        steps.append({"lambda": mp.nstr(lam, 5), "converged": r["converged"],
                      "grid_max": r.get("grid_max"),
                      "guard": r.get("guard"), "ratio": r.get("ratio"),
                      "iterations": r.get("iterations"),
                      "newton": r.get("newton_note"),
                      "seconds": round(time.perf_counter() - t0, 1)})
        print(f"lambda={mp.nstr(lam, 5)} conv={r['converged']} "
              f"grid_max={r.get('grid_max')} guard={r.get('guard')} "
              f"ratio={r.get('ratio')} newton={r.get('newton_note')} "
              f"({steps[-1]['seconds']}s)", flush=True)
        P = [mp.mpf(c) for c in r["P"]]
        Q = [mp.mpf(c) for c in r["Q"]]
        E = mp.mpf(r["E"])
        ref = [mp.mpf(t) for t in r["ref"]]

    s = R.score60(P, Q, a, W, target, grid=(ts, tsm, fs))
    s4 = R.score60(P, Q, a, W, target, grid=(ts4, tsm4, fs4))
    out = {"schema": "clean-edge-independent-homotopy-1",
           "width_name": args.width_name, "width": W, "a": a,
           "m": m, "k": k, "stored": m + k + 1, "target": target,
           "B": args.B, "A0": mp.nstr(A0, 17),
           "steps": steps,
           "final": {"e60": s["e60"], "e60_abs": s["e60_abs"],
                     "e60_at_t": s["e60_at_t"], "guard60": s["guard60"],
                     "meets60": s["meets60"], "e60_4001": s4["e60"],
                     "guard60_4001": s4["guard60"],
                     "meets60_4001": s4["meets60"],
                     "P70": [mp.nstr(c, 70) for c in P],
                     "Q70": [mp.nstr(c, 70) for c in Q]},
           "wall_seconds": round(time.perf_counter() - t_start, 1)}
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(json.dumps(out, indent=1, default=str),
                              encoding="utf-8")
    print(f"\nFINAL [{m}/{k}] e60={s['e60']} guard60={s['guard60']} "
          f"meets={s['meets60']} e60_4001={s4['e60']}")
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
