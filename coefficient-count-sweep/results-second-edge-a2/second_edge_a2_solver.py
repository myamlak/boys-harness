#!/usr/bin/env python3
"""Independent rational counts at the second left edge, a2 = 2 x0.

This is the verification harness of 2026-10-31 (`scan_cell.py`,
`scan_grid.py`, `finalize_grid.py` in the q330 verification tree) copied and
adapted to the second Boys left edge.  The numeric core, the acceptance rule
and the retry protocol are unchanged; only the left edge, the doubled widths
and the worklist bounds are new.

Problem
  Target function F_0(x) = sqrt(pi)/(2 sqrt(x)) erf(sqrt(x)) on
  [a, a + W] with a = a2 = 2 x0 and W = 2 W_x0 for each of the fifteen Group A
  widths.  The rational form is p(t)/(1 + q1 t + ... + qk t^k) in
  t = 2(x - a)/W - 1; the stored count is (m + 1) + k with k <= 12.

Acceptance (the paper's rule unchanged)
  Coefficients are converted exactly from the barycentric rational to
  monomials at 70 decimal digits, rounded to binary64, and evaluated by Horner
  in binary64 by both the fused (math.fma) and the separate route; the error
  is the worse of the two.  The target's kind is relative
  (max |fhat - f| / f) or absolute (max |fhat - f|).  The denominator is
  guarded at min |den| >= 1e-3 max |den|.  The 2001-node uniform-t grid
  certifies; the 4001-node grid checks.  The reference is F_0 at 70 decimal
  digits (mpmath).

Counting protocol
  The default BRASIL settings are scanned by stored count from 3 up to a
  per-cell bound (the corrected x0 count + 3).  If the first meeting level is
  s, level s - 1 is re-solved with the tightened settings; if it meets, the
  count drops and the next level below is retried.  If the default scan finds
  no meeting fit, the tightened settings are scanned upward instead.  A
  meeting level is one where at least one admissible [m/k] meets on the 2001
  grid and passes the denominator guard.  The witness is chosen preferring a
  fit that also meets on the 4001 grid and passes its guard, then a converged
  fit, then the smallest 4001 (or 2001) error.  The count is an existence
  upper bound under this solver, not a proof of minimality.

Usage
  python second_edge_a2_solver.py --edge a2 --only A_rel13_m7
  python second_edge_a2_solver.py --edge a2 --workers 6 --out-dir .
  python second_edge_a2_solver.py --edge x0 --only A_abs514_p1   # fidelity check

Per cell it writes cells/<cell_id>.json (the raw scan) and
cells/<cell_id>.final.json (the finalized record), and prints one DONE line.
"""
import argparse
import json
import math
import os
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import mpmath as mp
from scipy.special import erf
from baryrat import brasil

X0 = 11.899848152108484
A2 = 2.0 * X0
KMAX = 12
DEN_FLOOR = 1e-3
DEFAULT_SETTINGS = {"name": "default", "tol": 1e-4, "maxiter": 1000,
                    "npi": -30, "init_steps": 100, "max_step_size": 0.1}
TIGHT_SETTINGS = {"name": "tightened", "tol": 1e-8, "maxiter": 5000,
                  "npi": -40, "init_steps": 300, "max_step_size": 0.05}
TARGETS = [("rel13", "rel:1e-13", "rel", 1e-13),
           ("abs514", "abs:5e-14", "abs", 5e-14),
           ("rel14", "rel:1e-14", "rel", 1e-14)]
NAMES15 = ["m7", "m6", "m5", "m4", "m3", "m2", "m1",
           "p0", "p1", "p2", "p3", "p4", "p5", "ship", "paper"]

mp.mp.dps = 70


# --------------------------------------------------------------------------
# Numeric core, copied verbatim from scan_cell.py except for the left edge.
# --------------------------------------------------------------------------
def f0_np(x):
    x = np.asarray(x, dtype=float)
    return np.sqrt(np.pi) / (2.0 * np.sqrt(x)) * erf(np.sqrt(x))


def f0_mp(x):
    return mp.sqrt(mp.pi) / (2 * mp.sqrt(x)) * mp.erf(mp.sqrt(x))


def pairs_at(stored, kmax=KMAX):
    return [(stored - 1 - k, k) for k in range(1, min(kmax, stored - 2) + 1)]


def poly_mul(A, B):
    C = [mp.mpf(0)] * (len(A) + len(B) - 1)
    for i, ai in enumerate(A):
        for j, bj in enumerate(B):
            C[i + j] += ai * bj
    return C


def poly_add(A, B):
    n = max(len(A), len(B))
    C = [mp.mpf(0)] * n
    for i, a in enumerate(A):
        C[i] += a
    for i, b in enumerate(B):
        C[i] += b
    return C


def bary_to_monomial(r):
    """Exact conversion of the barycentric rational to monomial P/Q, Q0=1."""
    z = [mp.mpf(float(v)) for v in r.nodes]
    w = [mp.mpf(float(v)) for v in r.weights]
    f = [mp.mpf(float(v)) for v in r.values]
    n = len(z)
    P = [mp.mpf(0)]
    Q = [mp.mpf(0)]
    for j in range(n):
        L = [mp.mpf(1)]
        for k in range(n):
            if k != j:
                L = poly_mul(L, [-z[k], mp.mpf(1)])
        Q = poly_add(Q, [w[j] * c for c in L])
        P = poly_add(P, [f[j] * w[j] * c for c in L])
    q0 = Q[0]
    if q0 == 0:
        return None
    P = [c / q0 for c in P]
    Q = [c / q0 for c in Q]
    return P, Q


def horner_mp(cs, t):
    acc = mp.mpf(0)
    for c in reversed(cs):
        acc = acc * t + c
    return acc


def horner_fma(cs, t):
    acc = 0.0
    for c in cs[::-1]:
        acc = math.fma(acc, float(t), float(c))
    return acc


def horner_sep(cs, t):
    acc = 0.0
    for c in cs[::-1]:
        acc = acc * t + float(c)
    return acc


def evaluate_stored_binary64(P, Q, ts, fs_mp):
    """Rounded monomial p/q, both Horner routes, error against high-precision f."""
    pd = [float(c) for c in P]
    qd = [float(c) for c in Q[1:]]
    wr = 0.0
    wa = 0.0
    dmin = math.inf
    dmax = 0.0
    for fused in (True, False):
        for t, fm in zip(ts, fs_mp):
            t = float(t)
            num = horner_fma(pd, t) if fused else horner_sep(pd, t)
            hq = horner_fma(qd, t) if fused else horner_sep(qd, t)
            den = math.fma(hq, t, 1.0) if fused else hq * t + 1.0
            v = num / den
            f = float(fm)
            e = abs(v - f)
            if e > wa:
                wa = e
            if e / f > wr:
                wr = e / f
            d = abs(den)
            if d < dmin:
                dmin = d
            if d > dmax:
                dmax = d
    return {"bin64_rel": wr, "bin64_abs": wa, "den_ratio": dmin / dmax}


# --------------------------------------------------------------------------
# Cell geometry and one solve, adapted from scan_grid.py.
# --------------------------------------------------------------------------
def build_grid(a, W):
    n = 2000
    ts = np.array([2.0 * i / n - 1.0 for i in range(n + 1)])
    xs = a + W * (ts + 1.0) / 2.0
    fs_mp = [f0_mp(mp.mpf(float(x))) for x in xs]
    return ts, xs, fs_mp


def solve_one(a, kind, W, m, k, settings, ts, xs, fs_mp):
    f = lambda t: f0_np(a + W * (np.asarray(t) + 1.0) / 2.0)
    row = {"W": W, "m": m, "k": k, "stored": m + k + 1, "kind": kind,
           "settings": settings["name"]}
    t0 = time.time()
    try:
        kw = {key: settings[key] for key in
              ("tol", "maxiter", "npi", "init_steps", "max_step_size")
              if key in settings}
        r, info = brasil(f, (-1.0, 1.0), (m, k), info=True, **kw)
    except Exception as ex:  # noqa: BLE001 - record, never abort the scan
        row["exception"] = repr(ex)
        row["seconds"] = round(time.time() - t0, 3)
        return row
    try:
        vals = r(ts)
        fs = f0_np(xs)
        e = np.abs(vals - fs)
        row["brut_rel"] = float(np.max(e / fs))
        row["brut_abs"] = float(np.max(e))
    except Exception as ex:  # noqa: BLE001
        row["eval_exception"] = repr(ex)
    row["converged"] = bool(info.converged)
    row["deviation"] = float(info.deviation)
    row["info_error"] = float(info.error)
    try:
        P, Q = bary_to_monomial(r)
    except Exception as ex:  # noqa: BLE001
        row["convert_exception"] = repr(ex)
        row["seconds"] = round(time.time() - t0, 3)
        return row
    if Q is not None:
        row.update(evaluate_stored_binary64(P, Q, ts, fs_mp))
        row["P"] = [mp.nstr(c, 17) for c in P]
        row["Q"] = [mp.nstr(c, 17) for c in Q]
    row["seconds"] = round(time.time() - t0, 3)
    return row


def meets(row, kind, target):
    if row.get("exception") or row.get("eval_exception"):
        return False
    if row.get("den_ratio", 0.0) < DEN_FLOOR:
        return False
    if "bin64_rel" not in row:
        return False
    return (row["bin64_rel"] <= target if kind == "rel"
            else row["bin64_abs"] <= target)


def run_cell(spec):
    """Default-settings scan from stored 3 to the cell bound, stopping at the
    first stored level with a meeting admissible fit (scan_grid.run)."""
    mp.mp.dps = 70
    a, kind, target, W = spec["a"], spec["kind"], spec["target"], spec["W"]
    out = spec["out"]
    os.makedirs(os.path.dirname(out) or ".", exist_ok=True)
    ts, xs, fs_mp = build_grid(a, W)
    rows = []
    status = {"cell_id": spec["cell_id"], "edge": spec["edge"],
              "kind": kind, "target": target, "W": W, "a": a,
              "max_stored": spec["max_stored"], "settings": DEFAULT_SETTINGS,
              "rows": rows,
              "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ",
                                             time.gmtime())}
    first_meet = None
    for stored in range(3, spec["max_stored"] + 1):
        if stored - 2 > KMAX:
            pairs = [(stored - 1 - k, k) for k in range(1, KMAX + 1)]
        else:
            pairs = pairs_at(stored, KMAX)
        for m, k in pairs:
            row = solve_one(a, kind, W, m, k, DEFAULT_SETTINGS,
                            ts, xs, fs_mp)
            rows.append(row)
            print(f"[{spec['cell_id']}] S={stored} [{m}/{k}] "
                  f"conv={row.get('converged')} "
                  f"brut_rel={row.get('brut_rel', 0):.2e} "
                  f"d_rel={row.get('bin64_rel', 0):.2e} "
                  f"d_abs={row.get('bin64_abs', 0):.2e} "
                  f"den={row.get('den_ratio', 0):.1e} "
                  f"meet={meets(row, kind, target)} {row['seconds']}s",
                  flush=True)
        level = [r for r in rows if r["stored"] == stored]
        status["levels_complete"] = stored
        status["meets"] = [{"m": r["m"], "k": r["k"],
                            "bin64_rel": r.get("bin64_rel"),
                            "bin64_abs": r.get("bin64_abs"),
                            "den_ratio": r.get("den_ratio")}
                           for r in level if meets(r, kind, target)]
        with open(out, "w") as fh:
            json.dump(status, fh, indent=1)
        if status["meets"]:
            first_meet = stored
            break
    status["first_meet_default"] = first_meet
    with open(out, "w") as fh:
        json.dump(status, fh, indent=1)
    return {"cell_id": spec["cell_id"], "first_meet_default": first_meet,
            "levels_complete": status.get("levels_complete")}


_grid4001 = {}


def grid_4001(a, W):
    key = (a, W)
    if key not in _grid4001:
        n = 4000
        ts = np.array([2.0 * i / n - 1.0 for i in range(n + 1)])
        xs = a + W * (ts + 1.0) / 2.0
        fs = [f0_mp(mp.mpf(float(x))) for x in xs]
        _grid4001[key] = (ts, fs)
    return _grid4001[key]


def witness_from_row(a, row, W):
    P = [mp.mpf(s) for s in row["P"]]
    Q = [mp.mpf(s) for s in row["Q"]]
    ts, fs = grid_4001(a, W)
    r4 = evaluate_stored_binary64(P, Q, ts, fs)
    return {
        "m": row["m"], "k": row["k"], "stored": row["stored"],
        "route": row.get("settings", "default"),
        "converged": row.get("converged"),
        "deviation": row.get("deviation"),
        "rel_2001": row.get("bin64_rel"),
        "abs_2001": row.get("bin64_abs"),
        "den_2001": row.get("den_ratio"),
        "rel_4001": r4["bin64_rel"],
        "abs_4001": r4["bin64_abs"],
        "den_4001": r4["den_ratio"],
        "P": row["P"], "Q": row["Q"],
    }


def err(w, kind, grid):
    return w.get(f"{kind}_{grid}")


def ok(w, kind, target, grid):
    e = err(w, kind, grid)
    d = w.get(f"den_{grid}")
    return (e is not None and d is not None
            and e <= target and d >= DEN_FLOOR)


def choose(witnesses, kind, target):
    both = [w for w in witnesses if ok(w, kind, target, 2001)
            and ok(w, kind, target, 4001)]
    prim = [w for w in witnesses if ok(w, kind, target, 2001)]
    pool = both if both else prim
    if not pool:
        return None, None
    both_ids = {id(w) for w in both}
    pool.sort(key=lambda w: (0 if w["converged"] else 1,
                             err(w, kind, 4001) if id(w) in both_ids
                             else err(w, kind, 2001)))
    return pool[0], ("both" if both else "2001_only")


def solve_level(spec, stored, settings, ts, xs, fs_mp):
    a, kind, target, W = spec["a"], spec["kind"], spec["target"], spec["W"]
    if stored - 2 > KMAX:
        pairs = [(stored - 1 - k, k) for k in range(1, KMAX + 1)]
    else:
        pairs = pairs_at(stored, KMAX)
    rows = []
    for m, k in pairs:
        row = solve_one(a, kind, W, m, k, settings, ts, xs, fs_mp)
        rows.append(row)
        print(f"[{spec['cell_id']}] RETRY S={stored} [{m}/{k}] "
              f"conv={row.get('converged')} "
              f"d_rel={row.get('bin64_rel', 0):.2e} "
              f"d_abs={row.get('bin64_abs', 0):.2e} "
              f"den={row.get('den_ratio', 0):.1e} "
              f"meet={meets(row, kind, target)} {row['seconds']}s",
              flush=True)
    return rows


def finalize(spec):
    """One cell: default scan, tightened descent below the first meet, witness
    selection and the record (finalize_grid.finalize, retry protocol intact)."""
    mp.mp.dps = 70
    a, kind, target, W = spec["a"], spec["kind"], spec["target"], spec["W"]
    default = run_cell(spec)
    with open(spec["out"]) as fh:
        data = json.load(fh)
    rows = data["rows"]
    by_stored = {}
    for r in rows:
        by_stored.setdefault(r["stored"], []).append(r)
    default_meet = None
    for s0 in sorted(by_stored):
        if any(meets(r, kind, target) for r in by_stored[s0]):
            default_meet = s0
            break
    ts, xs, fs_mp = build_grid(a, W)
    retry_rows = []
    tightened_levels = []
    s = default_meet
    if s is None:
        for level in range(3, spec["max_stored"] + 1):
            retry_rows += solve_level(spec, level, TIGHT_SETTINGS,
                                      ts, xs, fs_mp)
            tightened_levels.append(level)
            lv = [r for r in retry_rows if r["stored"] == level]
            if any(meets(r, kind, target) for r in lv):
                s = level
                break
    else:
        while s > 3:
            level = s - 1
            if level in tightened_levels:
                break
            new = solve_level(spec, level, TIGHT_SETTINGS, ts, xs, fs_mp)
            retry_rows += new
            tightened_levels.append(level)
            if any(meets(r, kind, target) for r in new):
                s = level
            else:
                break
    all_by_stored = {}
    for r in rows + retry_rows:
        all_by_stored.setdefault(r["stored"], []).append(r)
    record = {
        "cell_id": spec["cell_id"], "edge": spec["edge"], "group": "A",
        "kind": kind, "target": target, "W": W, "a": a,
        "width_name": spec.get("width_name"),
        "target_tag": spec.get("target_tag"),
        "x0_count_used_for_bound": spec.get("x0_count_used_for_bound"),
        "max_stored": spec["max_stored"],
        "levels_complete": data.get("levels_complete"),
        "first_meet_default": default_meet,
        "tightened_levels": tightened_levels,
        "default_settings": DEFAULT_SETTINGS,
        "tightened_settings": TIGHT_SETTINGS,
    }
    if s is None:
        record.update({"corrected_count": None, "no_meet": True,
                       "pair": None, "witness": None})
        return record
    level_rows = list(all_by_stored[s])
    witnesses = [witness_from_row(a, r, W) for r in level_rows if "P" in r]
    wch, how = choose(witnesses, kind, target)
    if wch is not None and not (ok(wch, kind, target, 2001)
                                and ok(wch, kind, target, 4001)) \
            and s not in tightened_levels:
        extra = solve_level(spec, s, TIGHT_SETTINGS, ts, xs, fs_mp)
        retry_rows += extra
        tightened_levels.append(s)
        all_by_stored[s] = list(all_by_stored[s]) + extra
        level_rows = list(all_by_stored[s])
        witnesses = [witness_from_row(a, r, W)
                     for r in level_rows if "P" in r]
        wch, how = choose(witnesses, kind, target)
    below = s - 1
    below_rows = [r for r in all_by_stored.get(below, []) if "P" in r]
    admissible = [r for r in below_rows
                  if r.get("den_ratio", 0) >= DEN_FLOOR]
    best_below = None
    best_below_pair = None
    if admissible:
        def e(r):
            return r["bin64_rel"] if kind == "rel" else r["bin64_abs"]
        br = min(admissible, key=e)
        best_below = e(br)
        best_below_pair = [br["m"], br["k"]]
    any_below = None
    if below_rows:
        any_below = min((r["bin64_rel"] if kind == "rel" else r["bin64_abs"])
                        for r in below_rows)
    below_ratio = (best_below / target) if best_below is not None else None
    route = "tightened" if (s not in [r["stored"] for r in rows]
                            or any(w["route"] == "tightened"
                                   for w in witnesses if w is wch)) else "default"
    record.update({
        "corrected_count": s,
        "pair": [wch["m"], wch["k"]] if wch else None,
        "witness": wch,
        "witness_selection": how,
        "delivered_2001": err(wch, kind, 2001) if wch else None,
        "delivered_4001": err(wch, kind, 4001) if wch else None,
        "guard_2001": wch["den_2001"] if wch else None,
        "guard_4001": wch["den_4001"] if wch else None,
        "best_below": best_below,
        "best_below_pair": best_below_pair,
        "best_below_any": any_below,
        "below_ratio": below_ratio,
        "route": route,
        "count_lowered_by_tightened": (default_meet is not None
                                       and s < default_meet)
        or (default_meet is None),
        "meets_4001": bool(wch and ok(wch, kind, target, 4001)),
        "guard_4001_ok": bool(wch and wch["den_4001"] >= DEN_FLOOR),
        "witness_converged": bool(wch and wch["converged"]),
    })
    rec = record
    rec["margin_limited"] = bool(
        rec["corrected_count"] is not None and (
            not rec["witness_converged"] or not rec["meets_4001"]
            or not rec["guard_4001_ok"]
            or rec["count_lowered_by_tightened"]
            or (rec["below_ratio"] is not None and rec["below_ratio"] < 1.25)))
    rec["solid"] = (rec["corrected_count"] is not None
                    and not rec["margin_limited"])
    return rec


def task(spec):
    try:
        rec = finalize(spec)
    except Exception as ex:  # noqa: BLE001 - report, never lose the cell
        import traceback
        rec = {"cell_id": spec.get("cell_id"), "error": repr(ex),
               "traceback": traceback.format_exc(), "spec": spec}
    with open(spec["out"].replace(".json", ".final.json"), "w") as fh:
        json.dump(rec, fh, indent=1)
    return rec


# --------------------------------------------------------------------------
# Worklist.
# --------------------------------------------------------------------------
def load_widths(corrected_counts):
    doc = json.loads(Path(corrected_counts).read_text())
    widths = {name: float(doc["problem"]["widths"][name]) for name in NAMES15}
    counts = {}
    for cell in doc["cells"]:
        counts[(cell["target"], cell["width_name"])] = cell
    return widths, counts


def build_worklist(edge, widths, counts, only):
    a = A2 if edge == "a2" else X0
    work = []
    for tag, tlabel, kind, target in TARGETS:
        for name in NAMES15:
            rec = counts[(tlabel, name)]
            x0_count = rec["corrected_count"]
            W = widths[name] * 2.0 if edge == "a2" else widths[name]
            cell_id = f"{'A2' if edge == 'a2' else 'A'}_{tag}_{name}"
            work.append({
                "cell_id": cell_id, "edge": edge, "group": "A",
                "kind": kind, "target": target, "a": a, "W": W,
                "width_name": name, "target_tag": tag,
                "x0_count_used_for_bound": x0_count,
                "max_stored": x0_count + 3,
                "out": f"cells/{cell_id}.json",
            })
    if only:
        keep = set(only.split(","))
        work = [w for w in work if w["cell_id"] in keep]
    return work


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--edge", choices=("a2", "x0"), default="a2")
    ap.add_argument("--only", default=None,
                    help="comma-separated cell ids")
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--out-dir", default=str(Path(__file__).resolve().parent))
    ap.add_argument("--corrected-counts", default=str(
        Path(__file__).resolve().parent.parent
        / "results-corrected" / "corrected-counts.json"))
    ap.add_argument("--extend", type=int, default=6,
                    help="extra levels if a cell finds no meeting fit")
    args = ap.parse_args()
    os.chdir(args.out_dir)
    widths, counts = load_widths(args.corrected_counts)
    work = build_worklist(args.edge, widths, counts, args.only)
    if not work:
        raise SystemExit("no cells selected")
    results = []
    t0 = time.time()
    if args.workers <= 1:
        for spec in work:
            rec = task(spec)
            results.append(rec)
            print(f"DONE {len(results)}/{len(work)} {rec.get('cell_id')} "
                  f"count={rec.get('corrected_count')} pair={rec.get('pair')} "
                  f"route={rec.get('route')} "
                  f"below_ratio={rec.get('below_ratio')}", flush=True)
    else:
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            futs = {pool.submit(task, w): w for w in work}
            for i, fut in enumerate(as_completed(futs), 1):
                spec = futs[fut]
                try:
                    rec = fut.result()
                except Exception as ex:  # noqa: BLE001 - worker died
                    import traceback
                    rec = {"cell_id": spec.get("cell_id"),
                           "error": repr(ex),
                           "traceback": traceback.format_exc(),
                           "spec": spec}
                results.append(rec)
                print(f"DONE {i}/{len(work)} {rec.get('cell_id')} "
                      f"count={rec.get('corrected_count')} "
                      f"pair={rec.get('pair')} route={rec.get('route')} "
                      f"below_ratio={rec.get('below_ratio')}", flush=True)
    results.sort(key=lambda r: r.get("cell_id", ""))
    summary = {"generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ",
                                              time.gmtime()),
               "edge": args.edge, "a": A2 if args.edge == "a2" else X0,
               "wall_seconds": round(time.time() - t0, 1),
               "workers": args.workers, "cells": results}
    name = "a2_run_summary.json" if args.edge == "a2" else "x0_check_summary.json"
    with open(name, "w") as fh:
        json.dump(summary, fh, indent=1)
    bad = [r for r in results if r.get("error")]
    print(f"ALLDONE cells={len(results)} errors={len(bad)} "
          f"wall={summary['wall_seconds']}s", flush=True)
    return 0 if not bad else 2


if __name__ == "__main__":
    raise SystemExit(main())
