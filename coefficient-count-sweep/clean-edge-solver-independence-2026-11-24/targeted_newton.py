"""Targeted Newton-Remez for one stored SK fit, with few exchanges.

Cheap variant of the polish: Newton on a fixed reference, then one exchange,
repeated a few times.  Prints the grid error and the equioscillation ratio
after every Newton, so the convergence is visible.

Usage:
  python -I targeted_newton.py --json out/ship_s13.json --m 6 --k 6 \
      --width-name ship --out out/targeted_ship_66.json
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
    ap.add_argument("--json", required=True)
    ap.add_argument("--m", type=int, required=True)
    ap.add_argument("--k", type=int, required=True)
    ap.add_argument("--width-name", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--rounds", type=int, default=8)
    ap.add_argument("--newton", type=int, default=50)
    args = ap.parse_args(argv)
    mp.mp.dps = R.DPS
    a, W, target = load_width(args.width_name)
    prob = R.problem_f0(a, W)
    ts, tsm, fs = R.build_grid(a, W, 2000)
    ts4, tsm4, fs4 = R.build_grid(a, W, 4000)
    ws = [prob.weight(t, f) for t, f in zip(tsm, fs)]
    src = json.loads(Path(args.json).read_text(encoding="utf-8"))
    row = next(r for r in src["rows"]
               if r["m"] == args.m and r["k"] == args.k)
    P = [mp.mpf(c) for c in row["P70"]]
    Q = [mp.mpf(c) for c in row["Q70"]]
    m, k = args.m, args.k
    N = m + k + 2
    E = mp.mpf(0)

    def stats(P, Q):
        errs, qv = R.profile(P, Q, tsm, fs, ws)
        gm = max(abs(e) for e in errs)
        g = R.guard_ratio(qv)
        merged = R.refine_extrema(P, Q, ts, errs, prob)
        return gm, g, merged

    gm, g, merged = stats(P, Q)
    ref = R.choose_window(merged, N)
    print(f"start [{m}/{k}] e={mp.nstr(gm, 10)} guard={mp.nstr(g, 6)} "
          f"extrema={len(merged)} ref={len(ref) if ref else None}")
    best = (gm, g, list(P), list(Q), "start")
    rounds = []
    for rnd in range(args.rounds):
        if ref is None:
            print("no reference; stop")
            break
        t0 = time.perf_counter()
        P2, Q2, E2, ok, note, res = R.newton_solve(
            m, k, ref, prob, P, Q, E)
        P, Q, E = P2, Q2, E2
        gm, g, merged = stats(P, Q)
        ref_errs = [abs(R.error_at(P, Q, t, prob)[0]) for t in ref]
        ratio = gm / min(ref_errs)
        print(f"round {rnd}: newton_ok={ok} note={note} "
              f"e={mp.nstr(gm, 10)} guard={mp.nstr(g, 6)} "
              f"ratio={mp.nstr(ratio, 6)} extrema={len(merged)} "
              f"({time.perf_counter()-t0:.0f}s)", flush=True)
        rounds.append({"round": rnd, "newton_ok": ok, "note": note,
                       "grid_max": mp.nstr(gm, 17),
                       "guard": mp.nstr(g, 17),
                       "ratio": mp.nstr(ratio, 17),
                       "extrema": len(merged),
                       "seconds": round(time.perf_counter() - t0, 1)})
        if g >= R.DEN_FLOOR and gm < best[0]:
            best = (gm, g, list(P), list(Q), f"round {rnd}")
        if ok and ratio <= 1 + R.CONV_TOL:
            print("  CERTIFICATE reached")
            break
        new_ref = R.choose_window(merged, N)
        if new_ref is None or new_ref == ref:
            # keep newton-ing on the same reference
            if new_ref == ref:
                # perturb: choose the global alternating window from a finer
                # extrema set is already done; stop if no progress
                pass
        ref = new_ref if new_ref is not None else ref

    gm, g, Pb, Qb, tag = best
    s = R.score60(Pb, Qb, a, W, target, grid=(ts, tsm, fs))
    s4 = R.score60(Pb, Qb, a, W, target, grid=(ts4, tsm4, fs4))
    out = {"schema": "clean-edge-independent-targeted-newton-1",
           "width_name": args.width_name, "width": W, "a": a,
           "m": m, "k": k, "stored": m + k + 1, "target": target,
           "rounds": rounds,
           "best": {"tag": tag, "grid_max": mp.nstr(gm, 17),
                    "guard": mp.nstr(g, 17), "e60": s["e60"],
                    "e60_abs": s["e60_abs"], "e60_at_t": s["e60_at_t"],
                    "guard60": s["guard60"], "meets60": s["meets60"],
                    "e60_4001": s4["e60"], "guard60_4001": s4["guard60"],
                    "meets60_4001": s4["meets60"],
                    "P70": [mp.nstr(c, 70) for c in Pb],
                    "Q70": [mp.nstr(c, 70) for c in Qb]}}
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(json.dumps(out, indent=1, default=str),
                              encoding="utf-8")
    print(f"BEST {tag}: e60={s['e60']} guard60={s['guard60']} "
          f"meets={s['meets60']} e60_4001={s4['e60']}")
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
