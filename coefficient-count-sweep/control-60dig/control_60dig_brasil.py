#!/usr/bin/env python3
"""Evaluation-precision control on the corrected BRASIL record, at 60 digits.

The corrected record (`results-corrected/corrected-counts.json`) counts every
rational fit by its delivered error in binary64: the coefficients are converted
exactly from the barycentric solution to monomials at 70 decimal digits,
rounded to binary64, and evaluated by Horner in binary64 by both the fused and
the separate routes; the worse of the two errors is compared with the target.
At the deep targets the fit's own error and the binary64 rounding are of
comparable size, so the corrected record's count may be a statement about
binary64 storage and evaluation rather than about the fit.

This script is the corrected record's evaluation-precision control.  It keeps
the record's solver (BRASIL, baryrat 2.1.2), its settings (default and
tightened), its denominator cap (k <= 12), its 2,001-node certifying grid, its
target semantics (relative to F_0 at each node) and its 70-decimal-digit F_0
reference.  It changes one thing: a fit is accepted at 60-digit evaluation,
with its monomial coefficients held at 60 significant digits.

Per cell the search reproduces the record first: the record's own pair is
re-solved at the record's stored count with the record's settings and its
binary64 delivered error is compared with the record digit for digit.  A cell
that does not reproduce is reported as a blocker and the walk does not run.

The 60-digit count is the smallest stored count at which some admissible pair
meets the target at 60-digit evaluation.  The search starts at the record's
count and confirms its pair.  If it meets, every admissible pair is solved at
each count below, on the default settings first and on the tightened settings
only when the default route found no meeting fit and the negative answer has to
be certified; the walk stops at the first count where neither route meets, and
that full sweep is the record's below-count check.  If the record's pair misses
at 60 digits, every pair at the record's count is swept and, failing a meeting
fit there, counts above the record's are searched until one meets.

A count is an existence bound under this solver, not a proof of minimality, in
the same sense as the records it controls.

Usage:
  python control_60dig_brasil.py --cells A
  python control_60dig_brasil.py --cells B,C1,C2
  python control_60dig_brasil.py --cells endpoints
  python control_60dig_brasil.py --cells all --force
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
import time
import traceback
from datetime import datetime, timezone
from pathlib import Path

import mpmath as mp
import numpy as np
from baryrat import brasil
from scipy.special import erf

HERE = Path(__file__).resolve().parent
PAPER = HERE.parent
CORRECTED = PAPER / "results-corrected" / "corrected-counts.json"

CONVERT_DPS = 70          # the record's reference and monomial conversion
EVAL_DPS = 60             # the control's evaluation precision
K_MAX = 12                # the record's denominator cap
DEN_FLOOR = 1e-3          # the record's denominator guard
MAX_STORED = 40           # upward search bound when the record count misses
SCHEMA = "corrected-60dig-control-1"

DEFAULT_SETTINGS = {"name": "default", "tol": 1e-4, "maxiter": 1000,
                    "npi": -30, "init_steps": 100, "max_step_size": 0.1}
TIGHT_SETTINGS = {"name": "tightened", "tol": 1e-8, "maxiter": 5000,
                  "npi": -40, "init_steps": 300, "max_step_size": 0.05}

# cell key -> (width_name, target tag, kind, target value)
CELLS = {
    "A": ("paper", "rel:1e-15", "rel", 1e-15),
    "B": ("ship", "rel:1e-15", "rel", 1e-15),
    "C1": ("p5", "rel:1e-15", "rel", 1e-15),
    "C2": ("p5", "rel:1e-14", "rel", 1e-14),
    "D": ("p0", "rel:1e-15", "rel", 1e-15),
    "E-paper": ("paper", "rel:1e-08", "rel", 1e-08),
    "E-ship": ("ship", "rel:1e-08", "rel", 1e-08),
    "E-p5": ("p5", "rel:1e-08", "rel", 1e-08),
}
DEEP_CELLS = ("A", "B", "C1", "C2", "D")
ENDPOINT_CELLS = ("E-paper", "E-ship", "E-p5")

mp.mp.dps = CONVERT_DPS


def now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


# ---------------------------------------------------------------------------
# Numeric core: copied from the record's verification harness
# (results-second-edge-a2/second_edge_a2_solver.py, whose docstring names
# scan_cell.py / scan_grid.py / finalize_grid.py as the source).  Only the
# 60-digit scoring that the control adds is new.


def f0_np(x):
    x = np.asarray(x, dtype=float)
    return np.sqrt(np.pi) / (2.0 * np.sqrt(x)) * erf(np.sqrt(x))


def f0_mp(x):
    return mp.sqrt(mp.pi) / (2 * mp.sqrt(x)) * mp.erf(mp.sqrt(x))


def pairs_at(stored, kmax=K_MAX):
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
    """The record's counting route: rounded monomial p/q, both Horner routes.

    Returned errors are the worse of the fused and separate routes; the
    denominator ratio is over both routes, as in the record.
    """
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
    return {"rel": wr, "abs": wa, "den_ratio": dmin / dmax}


def build_grid(a, W, n=2000):
    ts = np.array([2.0 * i / n - 1.0 for i in range(n + 1)])
    xs = a + W * (ts + 1.0) / 2.0
    fs_mp = [f0_mp(mp.mpf(float(x))) for x in xs]
    return ts, xs, fs_mp


def mp_horner(cs, t):
    acc = mp.mpf(0)
    for c in reversed(cs):
        acc = acc * t + c
    return acc


def _mp_errors(P, Q, ts, fs):
    """Worst relative/absolute error and denominator ratio at current dps."""
    wr = mp.mpf(0)
    wa = mp.mpf(0)
    at = None
    dmin = mp.inf
    dmax = mp.mpf(0)
    for t, f in zip(ts, fs):
        tm = mp.mpf(float(t))
        fm = mp.mpf(f)
        num = mp_horner(P, tm)
        den = mp_horner(Q, tm)
        v = num / den
        e = abs(v - fm)
        if e > wa:
            wa = e
        rel = e / fm
        if rel > wr:
            wr, at = rel, t
        d = abs(den)
        if d < dmin:
            dmin = d
        if d > dmax:
            dmax = d
    return {"rel": wr, "abs": wa, "at": at, "den_ratio": dmin / dmax}


def score_at_60(P, Q, ts, fs, target, den_floor=DEN_FLOOR):
    """60-digit errors with coefficients held at 60 significant digits.

    e60   coefficients rounded to 60 significant digits, evaluated at 60 dps.
    e64r  coefficients rounded to binary64, evaluated at 60 dps.
    e64   binary64 coefficients and binary64 Horner, both routes (the record's
          counting route).
    """
    with mp.workdps(EVAL_DPS):
        P60 = [mp.mpf(mp.nstr(c, EVAL_DPS)) for c in P]
        Q60 = [mp.mpf(mp.nstr(c, EVAL_DPS)) for c in Q]
        P64 = [mp.mpf(float(c)) for c in P]
        Q64 = [mp.mpf(float(c)) for c in Q]
        vals = _mp_errors(P60, Q60, ts, fs)
        vals64r = _mp_errors(P64, Q64, ts, fs)
    b = evaluate_stored_binary64(P, Q, ts, fs)
    meet60 = bool(vals["rel"] <= target
                  and vals["den_ratio"] >= den_floor)
    return {
        "e60": mp.nstr(vals["rel"], 17),
        "e60_abs": mp.nstr(vals["abs"], 17),
        "e60_at_t": mp.nstr(vals["at"], 17),
        "guard60": mp.nstr(vals["den_ratio"], 17),
        "e64r": mp.nstr(vals64r["rel"], 17),
        "e64r_abs": mp.nstr(vals64r["abs"], 17),
        "e64": repr(b["rel"]),
        "e64_abs": repr(b["abs"]),
        "guard64": repr(b["den_ratio"]),
        "meets60": meet60,
    }


def score_4001(P, Q, ts4, fs4, target):
    """The check grid: the same 60-digit evaluation on 4,001 nodes."""
    with mp.workdps(EVAL_DPS):
        P60 = [mp.mpf(mp.nstr(c, EVAL_DPS)) for c in P]
        Q60 = [mp.mpf(mp.nstr(c, EVAL_DPS)) for c in Q]
        v = _mp_errors(P60, Q60, ts4, fs4)
    return {"e60_4001": mp.nstr(v["rel"], 17),
            "e60_4001_abs": mp.nstr(v["abs"], 17),
            "guard60_4001": mp.nstr(v["den_ratio"], 17),
            "meets60_4001": bool(v["rel"] <= target
                                 and v["den_ratio"] >= DEN_FLOOR)}


# ---------------------------------------------------------------------------
# One solve and one level


def solve_one(a, W, kind, target, m, k, settings, ts, xs, fs):
    """Solve one [m/k], score it on the 2,001-node grid, keep the fit."""
    started = time.perf_counter()
    row = {"m": m, "k": k, "stored": m + k + 1,
           "settings": settings["name"]}
    try:
        f = lambda t: f0_np(a + W * (np.asarray(t) + 1.0) / 2.0)
        kw = {key: settings[key] for key in
              ("tol", "maxiter", "npi", "init_steps", "max_step_size")}
        r, info = brasil(f, (-1.0, 1.0), (m, k), info=True, **kw)
    except Exception as ex:  # noqa: BLE001 - record, never abort the sweep
        row["exception"] = repr(ex)
        row["seconds"] = round(time.perf_counter() - started, 2)
        return row
    row["converged"] = bool(info.converged)
    row["deviation"] = float(info.deviation)
    row["info_error"] = float(info.error)
    try:
        converted = bary_to_monomial(r)
    except Exception as ex:  # noqa: BLE001
        row["convert_exception"] = repr(ex)
        row["seconds"] = round(time.perf_counter() - started, 2)
        return row
    if converted is None:
        row["convert_exception"] = "barycentric denominator constant is zero"
        row["seconds"] = round(time.perf_counter() - started, 2)
        return row
    P, Q = converted
    try:
        row.update(score_at_60(P, Q, ts, fs, target))
    except Exception as ex:  # noqa: BLE001
        row["score_exception"] = repr(ex)
        row["seconds"] = round(time.perf_counter() - started, 2)
        return row
    row["_P"] = P
    row["_Q"] = Q
    row["P60"] = [mp.nstr(c, EVAL_DPS) for c in P]
    row["Q60"] = [mp.nstr(c, EVAL_DPS) for c in Q]
    row["P70"] = [mp.nstr(c, CONVERT_DPS) for c in P]
    row["Q70"] = [mp.nstr(c, CONVERT_DPS) for c in Q]
    row["nodes"] = [repr(float(v)) for v in r.nodes]
    row["weights"] = [repr(float(v)) for v in r.weights]
    row["values"] = [repr(float(v)) for v in r.values]
    row["seconds"] = round(time.perf_counter() - started, 2)
    return row


def public_row(row):
    return {k: v for k, v in row.items() if k not in ("_P", "_Q")}


def meets(row):
    return bool(row.get("meets60"))


def best_row(rows):
    """The meeting row with the smallest 60-digit error."""
    passing = [r for r in rows if meets(r)]
    if not passing:
        return None
    return min(passing, key=lambda r: float(r["e60"]))


def best_e60(rows):
    """The smallest e60 among admissible rows (guard passing).

    A fit whose denominator guard fails is not admissible, so it does not
    bound the below-count margin; the closest admissible miss is what the
    record's below-count check leaves.
    """
    adm = [r for r in rows if r.get("e60")
           and float(r.get("guard60", 0.0)) >= DEN_FLOOR]
    vals = [float(r["e60"]) for r in adm]
    return repr(min(vals)) if vals else None


def level_sweep(cell, stored, settings, ts, xs, fs, verbose=True):
    """Solve every admissible pair at one stored count, one settings route."""
    pairs = pairs_at(stored)
    rows = []
    started = time.perf_counter()
    for m, k in pairs:
        row = solve_one(cell["a"], cell["width"], cell["kind"],
                        cell["target_value"], m, k, settings, ts, xs, fs)
        rows.append(row)
        if verbose:
            print(f"    {settings['name']:9s} S={stored} [{m}/{k}] "
                  f"conv={row.get('converged')} "
                  f"e60={row.get('e60')} e64={row.get('e64')} "
                  f"guard={row.get('guard60')} meet={meets(row)} "
                  f"{row['seconds']}s", flush=True)
    out = {"stored": stored, "settings": settings["name"],
           "pairs_tried": [f"[{m}/{k}]" for m, k in pairs],
           "solved": len([r for r in rows if "e60" in r]),
           "meets60": any(meets(r) for r in rows),
           "rows": rows,
           "seconds": round(time.perf_counter() - started, 1)}
    if verbose:
        best = best_row(rows)
        print(f"  level {stored} ({settings['name']}): "
              f"meets60={out['meets60']} best={best['e60'] if best else None}"
              f" ({out['seconds']} s)", flush=True)
    return out


def serialize_levels(levels):
    return [{**{k: v for k, v in lv.items() if k != "rows"},
             "rows": [public_row(r) for r in lv["rows"]]}
            for lv in levels]


def record_cell(doc, width_name, target):
    """The record's entry for one (width, target), from cells or ladder."""
    for c in doc["cells"]:
        if c["width_name"] == width_name and c["target"] == target:
            return {
                "source": "results-corrected/corrected-counts.json cells",
                "stored": c["corrected_count"],
                "pair": c["pair"],
                "delivered_2001": c["delivered_2001"],
                "delivered_4001": c.get("delivered_4001"),
                "guard": c.get("guard_2001"),
                "route": c.get("route"),
                "best_below": c.get("best_below"),
                "below_ratio": c.get("below_ratio"),
            }
    for c in doc["ladder"]:
        if c["width_name"] == width_name and c["target"] == target:
            return {
                "source": "results-corrected/corrected-counts.json ladder",
                "stored": c["rational_corrected"],
                "pair": c["rational_pair"],
                "delivered_2001": c["rational_delivered_2001"],
                "delivered_4001": None,
                "guard": c.get("rational_guard"),
                "route": None,
                "best_below": None,
                "below_ratio": None,
            }
    raise KeyError(f"no record entry for {width_name} {target}")


# ---------------------------------------------------------------------------
# The walk


def control_cell(key, doc, verbose=True):
    width_name, target_tag, kind, target_value = CELLS[key]
    width = float(doc["problem"]["widths"][width_name])
    a = float(doc["problem"]["x0"])
    rec = record_cell(doc, width_name, target_tag)
    cell = {"key": key, "width_name": width_name, "width": width,
            "a": a, "target": target_tag, "kind": kind,
            "target_value": target_value}
    print(f"[{key}] {width_name} {target_tag} (W={width!r}, "
          f"record {rec['stored']} {rec['pair']}, "
          f"delivered {rec['delivered_2001']!r})", flush=True)

    mp.mp.dps = CONVERT_DPS
    ts, xs, fs = build_grid(a, width)

    # ---- reproduction: the record's pair at the record's count
    rep = solve_one(a, width, kind, target_value, rec["pair"][0],
                    rec["pair"][1], DEFAULT_SETTINGS, ts, xs, fs)
    rep_route = "default"
    if rep.get("e64") != repr(rec["delivered_2001"]):
        rep2 = solve_one(a, width, kind, target_value, rec["pair"][0],
                         rec["pair"][1], TIGHT_SETTINGS, ts, xs, fs)
        if rep2.get("e64") == repr(rec["delivered_2001"]):
            rep, rep_route = rep2, "tightened"
    reproduced = rep.get("e64") == repr(rec["delivered_2001"])
    reproduction = {
        "record": rec,
        "settings": rep_route,
        "e64": rep.get("e64"),
        "e64_abs": rep.get("e64_abs"),
        "guard64": rep.get("guard64"),
        "converged": rep.get("converged"),
        "deviation": rep.get("deviation"),
        "exact_match": reproduced,
        "record_delivered_2001": rec["delivered_2001"],
        "record_guard": rec["guard"],
    }
    print(f"  reproduction: e64={rep.get('e64')} record="
          f"{rec['delivered_2001']!r} match={reproduced}", flush=True)
    if not reproduced:
        return {"schema": SCHEMA, "cell": cell, "reproduction": reproduction,
                "count60": None, "status": "reproduction-failed"}

    top = rec["stored"]
    levels = []
    confirm = {"stored": top, "settings": rep_route,
               "pairs_tried": [f"[{rec['pair'][0]}/{rec['pair'][1]}]"],
               "solved": 1, "meets60": meets(rep),
               "rows": [rep], "seconds": None,
               "label": "confirm-record-pair"}
    levels.append(confirm)
    final = None
    direction = None
    accepted_pool = []
    fail = None
    minimality = None
    if meets(rep):
        final, direction, accepted_pool = top, "equal", [rep]

    # ---- the record's count, every pair, default then tightened
    if final is None:
        lv = level_sweep(cell, top, DEFAULT_SETTINGS, ts, xs, fs, verbose)
        lv["label"] = "record-count-full-sweep-default"
        levels.append(lv)
        pool = [r for r in lv["rows"] if meets(r)]
        if pool:
            final, direction, accepted_pool = top, "equal-full-sweep", pool
        else:
            lv2 = level_sweep(cell, top, TIGHT_SETTINGS, ts, xs, fs, verbose)
            lv2["label"] = "record-count-full-sweep-tightened"
            levels.append(lv2)
            pool = [r for r in lv2["rows"] if meets(r)]
            if pool:
                final, direction, accepted_pool = top, "equal-full-sweep", pool

    # ---- upward, when nothing at the record's count meets
    if final is None:
        for stored in range(top + 1, MAX_STORED + 1):
            lv = level_sweep(cell, stored, DEFAULT_SETTINGS, ts, xs, fs,
                             verbose)
            lv["label"] = "upward-full-sweep-default"
            levels.append(lv)
            pool = [r for r in lv["rows"] if meets(r)]
            if pool:
                final, direction, accepted_pool = stored, "up", pool
                break
            lv2 = level_sweep(cell, stored, TIGHT_SETTINGS, ts, xs, fs,
                              verbose)
            lv2["label"] = "upward-full-sweep-tightened"
            levels.append(lv2)
            pool = [r for r in lv2["rows"] if meets(r)]
            if pool:
                final, direction, accepted_pool = stored, "up", pool
                break
        if final is not None:
            below = final - 1
            failed = [lv for lv in levels
                      if lv["stored"] == below and not lv["meets60"]
                      and lv["label"] != "confirm-record-pair"]
            fail = summarize_failure(below, failed)
            minimality = (
                "not searched below the record count: the record's own pair "
                "and every admissible pair at its count miss at 60 digits, so "
                "the 60-digit count lies above the record's; the full sweep at "
                f"stored {below} is the below-count check")

    # ---- downward, with the below-count check at the first failing count
    if final is not None and direction != "up":
        minimality = (
            "every admissible pair at each count below was solved on the "
            "default route; a count whose default route found no meeting fit "
            "was also swept on the tightened route; the first count where "
            "neither route meets is the below-count check")
        stored = final
        while stored > 3:
            below = stored - 1
            lv = level_sweep(cell, below, DEFAULT_SETTINGS, ts, xs, fs,
                             verbose)
            lv["label"] = "below-full-sweep-default"
            levels.append(lv)
            pool = [r for r in lv["rows"] if meets(r)]
            if pool:
                final, direction, accepted_pool = below, "down", pool
                stored = below
                continue
            lv2 = level_sweep(cell, below, TIGHT_SETTINGS, ts, xs, fs,
                              verbose)
            lv2["label"] = "below-full-sweep-tightened"
            levels.append(lv2)
            pool = [r for r in lv2["rows"] if meets(r)]
            if pool:
                final, direction, accepted_pool = below, "down", pool
                stored = below
                continue
            fail = summarize_failure(below, [lv, lv2])
            break

    accepted = best_row(accepted_pool) if accepted_pool else None
    acc = None
    if accepted is not None:
        ts4, xs4, fs4 = build_grid(a, width, n=4000)
        check = score_4001(accepted["_P"], accepted["_Q"], ts4, fs4,
                           target_value)
        acc = {
            "m": accepted["m"], "k": accepted["k"],
            "stored": accepted["stored"], "settings": accepted["settings"],
            "converged": accepted.get("converged"),
            "deviation": accepted.get("deviation"),
            "e60": accepted.get("e60"),
            "e60_abs": accepted.get("e60_abs"),
            "e60_at_t": accepted.get("e60_at_t"),
            "e64r": accepted.get("e64r"),
            "e64r_abs": accepted.get("e64r_abs"),
            "e64": accepted.get("e64"),
            "e64_abs": accepted.get("e64_abs"),
            "guard60": accepted.get("guard60"),
            "guard64": accepted.get("guard64"),
            "meets60": accepted.get("meets60"),
            "P60": accepted.get("P60"), "Q60": accepted.get("Q60"),
            "P70": accepted.get("P70"), "Q70": accepted.get("Q70"),
            "nodes": accepted.get("nodes"),
            "weights": accepted.get("weights"),
            "values": accepted.get("values"),
            **check,
        }

    result = {
        "schema": SCHEMA,
        "cell": cell,
        "reproduction": reproduction,
        "count60": final,
        "direction": direction,
        "minimality": minimality,
        "delta_vs_record": None if final is None else final - rec["stored"],
        "accepted": acc,
        "failing_level": fail,
        "levels": serialize_levels(levels),
    }
    if final is not None:
        print(f"  COUNT60 {final} ({direction}); accepted "
              f"[{acc['m']}/{acc['k']}] {acc['settings']} "
              f"e60={acc['e60']} e64r={acc['e64r']} e64={acc['e64']}",
              flush=True)
    else:
        print(f"  no 60-digit count found by stored {MAX_STORED}",
              flush=True)
    return result


def summarize_failure(stored, failed_levels):
    if not failed_levels:
        return None
    return {
        "stored": stored,
        "routes": [{"pass": lv["settings"], "solved": lv["solved"],
                    "meets60": lv["meets60"],
                    "best_e60": best_e60(lv["rows"])}
                   for lv in failed_levels],
        "rows": [public_row(r) for lv in failed_levels for r in lv["rows"]],
    }


# ---------------------------------------------------------------------------
# Driver


def _baryrat_version():
    try:
        import baryrat
        return getattr(baryrat, "__version__", "unknown")
    except Exception:  # noqa: BLE001
        return "unknown"


def meta_block():
    return {
        "generated_utc": now(),
        "script": str(Path(__file__).resolve()),
        "script_sha256": sha256(Path(__file__).resolve()),
        "corrected_record": str(CORRECTED),
        "corrected_record_sha256": sha256(CORRECTED),
        "corrected_readme_sha256": sha256(CORRECTED.parent / "README.md"),
        "versions": {"python": sys.version.split()[0],
                     "mpmath": mp.__version__,
                     "numpy": np.__version__,
                     "baryrat": _baryrat_version()},
        "settings": {
            "solver": "BRASIL (baryrat), the corrected record's solver",
            "default": DEFAULT_SETTINGS, "tightened": TIGHT_SETTINGS,
            "k_max": K_MAX, "den_floor": DEN_FLOOR,
            "nodes": 2001,
            "check_nodes": 4001,
            "convert_dps": CONVERT_DPS, "evaluation_dps": EVAL_DPS,
            "reference": "F_0 at 70 decimal digits (mpmath)",
            "acceptance60": (
                "coefficients from the exact monomial conversion at 70 "
                "digits held at 60 significant digits, evaluated at 60 dps "
                "on the 2,001-node uniform-t grid, relative to F_0 at each "
                "node; denominator guard min|den| >= 1e-3 max|den|"),
        },
    }


def main(argv=None):
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cells", default="A",
                    help="A,B,C1,C2 | endpoints | deep | all | comma list")
    ap.add_argument("--out", default=str(HERE / "summary.json"))
    ap.add_argument("--force", action="store_true",
                    help="re-run a cell whose record already exists")
    ap.add_argument("--max-seconds", type=float, default=0.0,
                    help="stop before starting a new cell once spent")
    args = ap.parse_args(argv)

    if args.cells == "all":
        wanted = list(CELLS)
    elif args.cells == "deep":
        wanted = list(DEEP_CELLS)
    elif args.cells == "endpoints":
        wanted = list(ENDPOINT_CELLS)
    else:
        wanted = [s.strip() for s in args.cells.split(",") if s.strip()]
    for k in wanted:
        if k not in CELLS:
            raise SystemExit(f"unknown cell {k!r}; known: {', '.join(CELLS)}")

    doc = json.loads(CORRECTED.read_text())
    if doc.get("schema") != "q330-corrected-counts-1":
        raise SystemExit(f"unexpected corrected record schema "
                         f"{doc.get('schema')!r}")
    started_all = time.perf_counter()
    results = []
    for key in wanted:
        out_path = HERE / f"cell_{key}.json"
        if out_path.exists() and not args.force:
            prev = json.loads(out_path.read_text())
            if prev.get("schema") == SCHEMA and (
                    prev.get("count60") is not None
                    or prev.get("status") == "reproduction-failed"):
                print(f"[{key}] reused from {out_path.name} "
                      f"(count60={prev.get('count60')})", flush=True)
                results.append(prev)
                continue
        if args.max_seconds and time.perf_counter() - started_all > args.max_seconds:
            print("budget spent; stopping before the next cell", flush=True)
            break
        try:
            res = control_cell(key, doc)
        except Exception as ex:  # noqa: BLE001 - record, never lose the run
            res = {"schema": SCHEMA, "cell": {"key": key},
                   "status": "error", "error": repr(ex),
                   "traceback": traceback.format_exc()}
        res["meta"] = meta_block()
        out_path.write_text(json.dumps(res, indent=1))
        results.append(res)
        print(f"  wrote {out_path}", flush=True)

    summary = {
        "schema": SCHEMA,
        "generated_utc": now(),
        "script_sha256": sha256(Path(__file__).resolve()),
        "corrected_record_sha256": sha256(CORRECTED),
        "versions": {"python": sys.version.split()[0],
                     "mpmath": mp.__version__, "numpy": np.__version__,
                     "baryrat": _baryrat_version()},
        "cells": [{k: r.get(k) for k in
                   ("cell", "reproduction", "count60", "direction",
                    "delta_vs_record", "accepted", "failing_level",
                    "status")}
                  for r in results],
        "runtime_seconds": round(time.perf_counter() - started_all, 1),
    }
    Path(args.out).write_text(json.dumps(summary, indent=1))
    print(f"wrote {args.out}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
