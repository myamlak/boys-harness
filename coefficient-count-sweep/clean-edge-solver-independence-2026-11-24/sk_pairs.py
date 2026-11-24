"""Independent pair sweep for one (width, stored count) of the clean-edge
record.

For every admissible pair [m/k] at the stored count, builds the iterated
weighted-least-squares (Sanathanan-Koerner) rational fit at 70 decimal digits
and scores it exactly as the control does: coefficients at 60 significant
digits, evaluated at 60 decimal digits on the 2,001-node grid, with F_0 at 70
digits and the 1e-3 denominator guard.

This is the independent, non-BRASIL arm.  A fit that meets at 60 digits is a
constructive existence proof at that stored count; the minimax Remez polish is
a separate arm (`remez_polish.py`).

Usage:
  python -I sk_pairs.py --width-name "2^5" --stored 14 --jobs 8 \
      --out out/p5_s14.json
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone
from pathlib import Path

import mpmath as mp

sys.path.insert(0, str(Path(__file__).resolve().parent))
import remez_rational as R

REC = Path(r"C:\Users\myaml\dev\papers-harness-repo\coefficient-count-sweep"
           r"\clean-edge-2026-11-22\edge62_sweep.json")
DPS = 70
EVAL_DPS = 60
DEN_FLOOR = 1e-3
K_MAX = 12
SK_ITERS = 15


def pairs_at(stored, kmax=K_MAX):
    return [(stored - 1 - k, k) for k in range(1, min(kmax, stored - 2) + 1)]


def sk_fit(m, k, tsm, fs, ws, iters=SK_ITERS):
    """Iterated linearised weighted least squares for P/Q, Q(0)=1.

    Minimises sum_i (P_i - f_i Q_i)^2 / (f_i Q_prev_i)^2 -- the weighted
    relative-error residual linearised around the previous denominator.
    """
    nb = m + 1 + k
    Q = [mp.mpf(1)] + [mp.mpf(0)] * k
    grad = []
    best = None
    for it in range(iters):
        M = [[mp.mpf(0)] * nb for _ in range(nb)]
        rhs = [mp.mpf(0)] * nb
        for t, f, w in zip(tsm, fs, ws):
            q = R.poly_eval(Q, t)
            if abs(q) < mp.mpf("1e-40"):
                return None, None, {"status": "denominator-vanished",
                                    "iteration": it}
            rho = (w / q) ** 2
            phi = []
            tp = mp.mpf(1)
            for _ in range(m + 1):
                phi.append(tp)
                tp *= t
            tp = mp.mpf(1)
            for _ in range(k):
                tp *= t
                phi.append(-f * tp)
            for aidx in range(nb):
                rhs[aidx] += rho * phi[aidx] * f
                for bidx in range(aidx, nb):
                    M[aidx][bidx] += rho * phi[aidx] * phi[bidx]
        for aidx in range(nb):
            for bidx in range(aidx):
                M[aidx][bidx] = M[bidx][aidx]
        try:
            u = mp.lu_solve(mp.matrix(M), mp.matrix(rhs))
        except (ZeroDivisionError, ValueError) as ex:
            return None, None, {"status": "sk-singular", "error": repr(ex),
                                "iteration": it}
        P = [u[j] for j in range(m + 1)]
        Q = [mp.mpf(1)] + [u[m + 1 + j] for j in range(k)]
        errs, qv = R.profile(P, Q, tsm, fs, ws)
        gm = max(abs(e) for e in errs)
        g = R.guard_ratio(qv)
        grad.append((it, gm, g))
        if best is None or (g >= DEN_FLOOR and gm < best[0]):
            best = (gm, g, it, list(P), list(Q))
        if len(grad) >= 2:
            prev = grad[-2][1]
            if prev != 0 and abs(gm - prev) <= mp.mpf("1e-6") * prev:
                break
    gm, g, it, P, Q = best
    return P, Q, {"status": "ok", "iterations": it + 1,
                  "grid_max_70": mp.nstr(gm, 17),
                  "guard_70": mp.nstr(g, 17),
                  "history": [(i, mp.nstr(v, 6), mp.nstr(gg, 6))
                              for i, v, gg in grad]}


def worker(task):
    width_f, stored, pairs, target, a, out_dir = task
    mp.mp.dps = DPS
    prob = R.problem_f0(a, width_f)
    ts, tsm, fs = R.build_grid(a, width_f, 2000)
    ts4, tsm4, fs4 = R.build_grid(a, width_f, 4000)
    ws = [prob.weight(t, f) for t, f in zip(tsm, fs)]
    rows = []
    for m, k in pairs:
        t0 = time.perf_counter()
        P, Q, info = sk_fit(m, k, tsm, fs, ws)
        row = {"m": m, "k": k, "stored": m + k + 1, "sk": info}
        if P is None:
            row["seconds"] = round(time.perf_counter() - t0, 2)
            rows.append(row)
            continue
        gm = mp.mpf(info["grid_max_70"])
        g = mp.mpf(info["guard_70"])
        s = R.score60(P, Q, a, width_f, target, grid=(ts, tsm, fs))
        s4 = R.score60(P, Q, a, width_f, target, grid=(ts4, tsm4, fs4))
        row.update({
            "grid_max_70": info["grid_max_70"],
            "guard_70": info["guard_70"],
            "e60": s["e60"], "e60_abs": s["e60_abs"],
            "e60_at_t": s["e60_at_t"], "guard60": s["guard60"],
            "meets60": s["meets60"],
            "e60_4001": s4["e60"], "guard60_4001": s4["guard60"],
            "meets60_4001": s4["meets60"],
            "P70": [mp.nstr(c, 70) for c in P],
            "Q70": [mp.nstr(c, 70) for c in Q],
            "seconds": round(time.perf_counter() - t0, 2),
        })
        rows.append(row)
    return {"width_f": width_f, "stored": stored, "target": target,
            "rows": rows}


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--width-name", required=True)
    ap.add_argument("--target", default="rel:1e-15")
    ap.add_argument("--stored", type=int, required=True)
    ap.add_argument("--jobs", type=int, default=8)
    ap.add_argument("--out", required=True)
    ap.add_argument("--pairs", default=None,
                    help="optional comma list m/k,m/k")
    args = ap.parse_args(argv)
    doc = json.loads(REC.read_text(encoding="utf-8"))
    a = float(doc["meta"]["edge"])
    row = next(r for r in doc["rows"]
               if r["width_name"] == args.width_name
               and r["target"] == args.target)
    width_f = float(row["width_f"])
    target = float(row["target_value"])
    pairs = (pairs_at(args.stored) if args.pairs is None
             else [tuple(int(x) for x in p.split("/"))
                   for p in args.pairs.split(",")])
    # chunk the pairs over workers; each worker builds the grid once
    nchunks = min(args.jobs, len(pairs))
    chunks = [pairs[i::nchunks] for i in range(nchunks)]
    tasks = [(width_f, args.stored, ch, target, a, str(Path(args.out).parent))
             for ch in chunks if ch]
    started = time.perf_counter()
    results = []
    if nchunks == 1:
        results = [worker(t) for t in tasks]
    else:
        with ProcessPoolExecutor(max_workers=nchunks) as pool:
            results = list(pool.map(worker, tasks))
    rows = [r for res in results for r in res["rows"]]
    rows.sort(key=lambda r: (r["m"], r["k"]))
    wall = round(time.perf_counter() - started, 1)
    out = {
        "schema": "clean-edge-independent-sk-1",
        "generated_utc": datetime.now(timezone.utc).strftime(
            "%Y-%m-%dT%H:%M:%SZ"),
        "width_name": args.width_name, "width": width_f, "a": a,
        "stored": args.stored, "target": target,
        "pairs": pairs, "jobs": nchunks, "wall_seconds": wall,
        "settings": {
            "solver": "iterated weighted least squares "
                      "(Sanathanan-Koerner), 70 dps",
            "sk_iters_max": SK_ITERS,
            "grid_nodes": 2001, "check_nodes": 4001,
            "k_max": K_MAX, "den_floor": DEN_FLOOR,
            "evaluation": "coefficients at 60 significant digits, "
                          "evaluated at 60 dps, F_0 at 70 dps",
            "target": "relative to F_0 at each node",
        },
        "rows": rows,
        "meeting_pairs": [[r["m"], r["k"]] for r in rows
                          if r.get("meets60")],
        "best_e60": min((mp.mpf(r["e60"]) for r in rows
                         if r.get("e60")), default=None),
    }
    out["best_e60"] = (mp.nstr(out["best_e60"], 17)
                       if out["best_e60"] is not None else None)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(json.dumps(out, indent=1, default=str),
                              encoding="utf-8")
    print(f"{args.width_name} stored={args.stored} wall={wall}s")
    for r in rows:
        print(f"  [{r['m']:2d}/{r['k']:2d}] e60={r.get('e60')} "
              f"guard={r.get('guard60')} meet={r.get('meets60')} "
              f"e60_4001={r.get('e60_4001')} "
              f"sk={r['sk'].get('status')} {r.get('seconds')}s")
    print(f"meeting pairs: {out['meeting_pairs']}; best e60: "
          f"{out['best_e60']}")
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
