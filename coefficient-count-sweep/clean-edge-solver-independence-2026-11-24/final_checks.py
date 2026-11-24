"""Final checks on the two independent meeting fits.

Loads the p5 [6/8] SK fit (out/p5_s15.json) and the ship [6/6] polished fit
(out/targeted_ship_66.json), scores both at 60 digits on the control's
2,001-node grid, its 4,001-node check grid and a denser 20,001-node grid, and
records their alternating-extrema structure and de la Vallee Poussin bounds.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import mpmath as mp

sys.path.insert(0, str(Path(__file__).resolve().parent))
import remez_rational as R

REC = Path(r"C:\Users\myaml\dev\papers-harness-repo\coefficient-count-sweep"
           r"\clean-edge-2026-11-22\edge62_sweep.json")
HERE = Path(__file__).resolve().parent


def width_of(name):
    doc = json.loads(REC.read_text(encoding="utf-8"))
    row = next(r for r in doc["rows"] if r["width_name"] == name
               and r["target"] == "rel:1e-15")
    return float(doc["meta"]["edge"]), float(row["width_f"])


def fit_of(path, m, k):
    doc = json.loads(Path(path).read_text(encoding="utf-8"))
    key = "best" if "best" in doc else None
    if key:
        b = doc["best"]
        return [mp.mpf(c) for c in b["P70"]], [mp.mpf(c) for c in b["Q70"]]
    row = next(r for r in doc["rows"] if r["m"] == m and r["k"] == k)
    return [mp.mpf(c) for c in row["P70"]], [mp.mpf(c) for c in row["Q70"]]


def grid_node_bound(P, Q, prob, tsm, ts, fs, ws, N):
    errs, qv = R.profile(P, Q, tsm, fs, ws)
    cand = [0]
    for i in range(1, len(errs) - 1):
        a, b, c = abs(errs[i - 1]), abs(errs[i]), abs(errs[i + 1])
        if b >= a and b >= c and (b > a or b > c):
            cand.append(i)
    cand.append(len(errs) - 1)
    merged = []
    for idx in cand:
        e = errs[idx]
        if e == 0:
            continue
        s = 1 if e > 0 else -1
        if merged and (1 if merged[-1][1] > 0 else -1) == s:
            if abs(e) > abs(merged[-1][1]):
                merged[-1] = (tsm[idx], e)
        else:
            merged.append((tsm[idx], e))
    w = R.choose_window(merged, N)
    L = (min(abs(R.error_at(P, Q, t, prob)[0]) for t in w)
         if w else None)
    return len(merged), L


def main():
    mp.mp.dps = R.DPS
    report = {}
    p5_path = (HERE / "out/targeted_p5_68.json"
               if (HERE / "out/targeted_p5_68.json").exists()
               else HERE / "out/p5_s15.json")
    cases = [
        ("2^5", "p5", p5_path, 6, 8),
        ("ship", "ship", HERE / "out/targeted_ship_66.json", 6, 6),
    ]
    for name, token, path, m, k in cases:
        a, W = width_of(name)
        prob = R.problem_f0(a, W)
        P, Q = fit_of(path, m, k)
        entry = {"pair": [m, k], "stored": m + k + 1, "grids": {}}
        for n, label in ((2000, "2001"), (4000, "4001"), (20000, "20001")):
            ts, tsm, fs = R.build_grid(a, W, n)
            s = R.score60(P, Q, a, W, 1e-15, grid=(ts, tsm, fs))
            entry["grids"][label] = {
                "e60": s["e60"], "e60_abs": s["e60_abs"],
                "e60_at_t": s["e60_at_t"], "guard60": s["guard60"],
                "meets60": s["meets60"]}
            print(f"[{name} {m}/{k}] {label}: e60={s['e60']} "
                  f"guard={s['guard60']} meets={s['meets60']}")
            if label == "2001":
                ts1, tsm1, fs1 = ts, tsm, fs
                ws1 = [prob.weight(t, f) for t, f in zip(tsm, fs)]
        N = m + k + 2
        M, L = grid_node_bound(P, Q, prob, tsm1, ts1, fs1, ws1, N)
        entry["extrema"] = M
        entry["de_la_vallee_bound_grid"] = (mp.nstr(L, 17)
                                            if L is not None else None)
        print(f"  extrema={M} N={N} de la Vallee bound={entry['de_la_vallee_bound_grid']}")
        report[token] = entry
    (HERE / "out/final_checks.json").write_text(
        json.dumps(report, indent=1), encoding="utf-8")
    print(f"wrote {HERE / 'out/final_checks.json'}")


if __name__ == "__main__":
    main()
