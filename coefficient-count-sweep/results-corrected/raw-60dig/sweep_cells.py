#!/usr/bin/env python3
"""60-digit control sweep over the 45 Group A cells of corrected-counts.json.

The solver, grids, reference, conversion precision, Horner routes and
denominator guard are the tagged v1.0.39 control's own, imported from
`coefficient-count-sweep/control-60dig/control_60dig_brasil.py` extracted with
`git archive v1.0.39`.  Two functions are wrapped so that a cell whose target
is absolute is scored in its target's kind: the tagged script was written for
the five relative deep cells and hard-codes the relative metric.  Nothing else
changes.

Usage:
  python sweep_cells.py --index 0          # one cell, by corrected-record index
  python sweep_cells.py --index 0 --force  # ignore an existing raw record
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import sys
import time
import traceback
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parent
CTRL_DIR = ROOT / "coefficient-count-sweep" / "control-60dig"
SWEEP = ROOT / "sweep"
LOG = ROOT / "logs"

CTRL_PATH = CTRL_DIR / "control_60dig_brasil.py"
spec = importlib.util.spec_from_file_location("ctl_v1039", CTRL_PATH)
ctl = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ctl)
mp = ctl.mp

SCHEMA = "c60-sweep-cell-1"
CONVERT_DPS = ctl.CONVERT_DPS
EVAL_DPS = ctl.EVAL_DPS
DEN_FLOOR = ctl.DEN_FLOOR
CURRENT_KIND = "rel"


def now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


# ---------------------------------------------------------------------------
# Target-kind-aware scoring.  The tagged control's score_at_60/score_4001
# always compare the relative error; the corrected record's absolute target
# `abs:5e-14` is compared in the absolute error.  These wrappers reproduce the
# tagged functions exactly for relative targets and select the absolute metric
# for absolute targets; every other ingredient is unchanged.


def score_at_60_kind(P, Q, ts, fs, target, den_floor=DEN_FLOOR):
    kind = CURRENT_KIND
    with mp.workdps(EVAL_DPS):
        P60 = [mp.mpf(mp.nstr(c, EVAL_DPS)) for c in P]
        Q60 = [mp.mpf(mp.nstr(c, EVAL_DPS)) for c in Q]
        P64 = [mp.mpf(float(c)) for c in P]
        Q64 = [mp.mpf(float(c)) for c in Q]
        vals = ctl._mp_errors(P60, Q60, ts, fs)
        vals64r = ctl._mp_errors(P64, Q64, ts, fs)
    b = ctl.evaluate_stored_binary64(P, Q, ts, fs)
    e60 = vals[kind]
    e64r = vals64r[kind]
    e64 = b[kind]
    return {
        "metric": kind,
        "e60": mp.nstr(e60, 17),
        "e60_rel": mp.nstr(vals["rel"], 17),
        "e60_abs": mp.nstr(vals["abs"], 17),
        "e60_at_t": mp.nstr(vals["at"], 17),
        "guard60": mp.nstr(vals["den_ratio"], 17),
        "e64r": mp.nstr(e64r, 17),
        "e64r_rel": mp.nstr(vals64r["rel"], 17),
        "e64r_abs": mp.nstr(vals64r["abs"], 17),
        "e64": repr(e64),
        "e64_rel": repr(b["rel"]),
        "e64_abs": repr(b["abs"]),
        "guard64": repr(b["den_ratio"]),
        "meets60": bool(e60 <= target and vals["den_ratio"] >= den_floor),
    }


def score_4001_kind(P, Q, ts4, fs4, target):
    kind = CURRENT_KIND
    with mp.workdps(EVAL_DPS):
        P60 = [mp.mpf(mp.nstr(c, EVAL_DPS)) for c in P]
        Q60 = [mp.mpf(mp.nstr(c, EVAL_DPS)) for c in Q]
        v = ctl._mp_errors(P60, Q60, ts4, fs4)
    e60 = v[kind]
    return {
        "metric": kind,
        "e60_4001": mp.nstr(e60, 17),
        "e60_4001_rel": mp.nstr(v["rel"], 17),
        "e60_4001_abs": mp.nstr(v["abs"], 17),
        "guard60_4001": mp.nstr(v["den_ratio"], 17),
        "meets60_4001": bool(e60 <= target and v["den_ratio"] >= DEN_FLOOR),
    }


ctl.score_at_60 = score_at_60_kind
ctl.score_4001 = score_4001_kind


# ---------------------------------------------------------------------------
# Measurement helpers


_GRID_CACHE = {}


def grid(a, W, n):
    key = (a, W, n)
    if key not in _GRID_CACHE:
        _GRID_CACHE[key] = ctl.build_grid(a, W, n=n)
    return _GRID_CACHE[key]


def pq_from_row(row):
    P = [mp.mpf(s) for s in row["P70"]]
    Q = [mp.mpf(s) for s in row["Q70"]]
    return P, Q


def ptrim_mpf(P):
    P = list(P)
    while len(P) > 1 and P[-1] == 0:
        P.pop()
    return P


def eval_fit(P, Q, ts, fs, kind, target):
    """The two scoring bases on one grid, in the target's kind."""
    with mp.workdps(EVAL_DPS):
        P60 = [mp.mpf(mp.nstr(c, EVAL_DPS)) for c in P]
        Q60 = [mp.mpf(mp.nstr(c, EVAL_DPS)) for c in Q]
        P64 = [mp.mpf(float(c)) for c in P]
        Q64 = [mp.mpf(float(c)) for c in Q]
        v60 = ctl._mp_errors(P60, Q60, ts, fs)
        v64r = ctl._mp_errors(P64, Q64, ts, fs)
    b = ctl.evaluate_stored_binary64(P, Q, ts, fs)
    e60 = v60[kind]
    e64r = v64r[kind]
    e64 = b[kind]
    return {
        "metric": kind,
        "e60": mp.nstr(e60, 17),
        "e60_rel": mp.nstr(v60["rel"], 17),
        "e60_abs": mp.nstr(v60["abs"], 17),
        "e60_at_t": mp.nstr(v60["at"], 17),
        "e64r": mp.nstr(e64r, 17),
        "e64r_rel": mp.nstr(v64r["rel"], 17),
        "e64r_abs": mp.nstr(v64r["abs"], 17),
        "e64": repr(e64),
        "e64_rel": repr(b["rel"]),
        "e64_abs": repr(b["abs"]),
        "guard60": mp.nstr(v60["den_ratio"], 17),
        "guard64": repr(b["den_ratio"]),
        "meets60": bool(e60 <= target and v60["den_ratio"] >= DEN_FLOOR),
        "meets64": bool(e64 <= target and b["den_ratio"] >= DEN_FLOOR),
    }


def eval_both_grids(P, Q, a, W, kind, target):
    ts, xs, fs = grid(a, W, 2000)
    ts4, xs4, fs4 = grid(a, W, 4000)
    return {
        "grid_2001": eval_fit(P, Q, ts, fs, kind, target),
        "grid_4001": eval_fit(P, Q, ts4, fs4, kind, target),
    }


def fit_block(P, Q, m, k, a, W, kind, target):
    full = eval_both_grids(P, Q, a, W, kind, target)
    Pn = P[: m + 1]
    Qn = Q[: k + 1]
    nominal = eval_both_grids(Pn, Qn, a, W, kind, target)
    return {
        "full": full,
        "nominal": nominal,
        "full_count": len(ptrim_mpf(P)) + len(ptrim_mpf(Q)) - 1,
        "nominal_count": (m + 1) + k,
        "disagree_60_2001": full["grid_2001"]["meets60"]
        != nominal["grid_2001"]["meets60"],
        "disagree_64_2001": full["grid_2001"]["meets64"]
        != nominal["grid_2001"]["meets64"],
        "disagree_60_4001": full["grid_4001"]["meets60"]
        != nominal["grid_4001"]["meets60"],
        "disagree_64_4001": full["grid_4001"]["meets64"]
        != nominal["grid_4001"]["meets64"],
    }


def nominal_check(P, Q, m, k, a, W, kind, target):
    Pn, Qn = P[: m + 1], Q[: k + 1]
    return eval_both_grids(Pn, Qn, a, W, kind, target)


def nominal_checks(rows, a, W, kind, target_value, seen=None):
    """Nominal-truncation verdicts for every solved row, deduplicated."""
    checks = []
    seen = set() if seen is None else seen
    for r in rows:
        if "P70" not in r:
            continue
        rm, rk = r["m"], r["k"]
        sig = (rm, rk, r.get("settings"))
        if sig in seen:
            continue
        seen.add(sig)
        Pr, Qr = pq_from_row(r)
        nc, ncerr = try_call(nominal_check, Pr, Qr, rm, rk, a, W, kind,
                             target_value)
        if ncerr is not None:
            checks.append({
                "pair": [rm, rk],
                "nominal_count": (rm + 1) + rk,
                "settings": r.get("settings"),
                "full_meets60": r.get("meets60"),
                "full_e60": r.get("e60"),
                "nominal": None,
                "nominal_error": ncerr,
                "nominal_meets60_2001": None,
                "nominal_meets64_2001": None,
                "nominal_meets60_4001": None,
                "nominal_meets64_4001": None,
            })
            continue
        checks.append({
            "pair": [rm, rk],
            "nominal_count": (rm + 1) + rk,
            "settings": r.get("settings"),
            "full_meets60": r.get("meets60"),
            "full_e60": r.get("e60"),
            "nominal": nc,
            "nominal_meets60_2001": nc["grid_2001"]["meets60"],
            "nominal_meets64_2001": nc["grid_2001"]["meets64"],
            "nominal_meets60_4001": nc["grid_4001"]["meets60"],
            "nominal_meets64_4001": nc["grid_4001"]["meets64"],
        })
    return checks


_BULK = ("P60", "Q60", "nodes", "weights", "values")


def try_call(fn, *args, **kwargs):
    """Run a measurement, returning (value, error) and never aborting a cell."""
    try:
        return fn(*args, **kwargs), None
    except Exception as ex:  # noqa: BLE001 - recorded as a failed measurement
        return None, repr(ex)


def slim_row(row):
    return {k: v for k, v in row.items() if k not in _BULK}


def slim_control(res):
    """Drop the barycentric node/weight/value dumps; keep every score field.

    The node values are the control's own intermediate representation, not a
    result; P70/Q70, all errors and guards are kept, so every number reported
    here can be recomputed from this file with the tagged script.
    """
    if "levels" in res:
        res["levels"] = [
            {k: v for k, v in lv.items() if k != "rows"}
            | {"rows": [slim_row(r) for r in lv["rows"]]}
            for lv in res["levels"]]
    if isinstance(res.get("accepted"), dict):
        res["accepted"] = slim_row(res["accepted"])
    if isinstance(res.get("failing_level"), dict):
        fl = res["failing_level"]
        if "rows" in fl:
            fl["rows"] = [slim_row(r) for r in fl["rows"]]
    return res


# ---------------------------------------------------------------------------
# One cell


def run_cell(index, force=False):
    global CURRENT_KIND
    corrected = json.loads(
        (ROOT / "coefficient-count-sweep" / "results-corrected"
         / "corrected-counts.json").read_text())
    cheb = json.loads(
        (ROOT / "coefficient-count-sweep" / "records"
         / "chebyshev_recount.json").read_text())
    rootiso = json.loads(
        (ROOT / "coefficient-count-sweep" / "strong-optional"
         / "root-isolation" / "records" / "root_isolation.json").read_text())

    out_path = SWEEP / f"cell_g{index:02d}.json"
    if out_path.exists() and not force:
        print(f"[g{index:02d}] reused {out_path.name}")
        return json.loads(out_path.read_text())

    entry = corrected["cells"][index]
    width_name = entry["width_name"]
    target = entry["target"]
    kind = "rel" if target.startswith("rel") else "abs"
    target_value = float(target.split(":")[1])
    W = float(corrected["problem"]["widths"][width_name])
    a = float(corrected["problem"]["x0"])
    m, k = entry["pair"]
    CURRENT_KIND = kind

    # The record's own entry, as the tagged control reads it.
    rec = ctl.record_cell(corrected, width_name, target)

    # Route the tagged control to this cell under a unique key.
    key = f"g{index:02d}"
    ctl.CELLS = {key: (width_name, target, kind, target_value)}
    print(f"[{key}] index={index} {width_name} {target} W={W!r} "
          f"record {rec['stored']} [{m}/{k}] delivered={rec['delivered_2001']!r}",
          flush=True)

    started = time.perf_counter()
    try:
        res = ctl.control_cell(key, corrected)
    except Exception as ex:  # noqa: BLE001 - record, never lose the sweep
        res = {"schema": ctl.SCHEMA, "cell": {"key": key},
               "status": "error", "error": repr(ex),
               "traceback": traceback.format_exc()}
    control_seconds = round(time.perf_counter() - started, 1)
    mp.mp.dps = CONVERT_DPS

    cell_doc = {
        "schema": SCHEMA,
        "generated_utc": now(),
        "index": index,
        "key": key,
        "width_name": width_name,
        "width": W,
        "a": a,
        "target": target,
        "target_kind": kind,
        "target_value": target_value,
        "recorded_pair": [m, k],
        "nominal_count": (m + 1) + k,
        "corrected_count": entry["corrected_count"],
        "delivered_2001": entry["delivered_2001"],
        "delivered_4001": entry.get("delivered_4001"),
        "guard_2001": entry.get("guard_2001"),
        "route": entry.get("route"),
        "control_seconds": control_seconds,
        "control": res,
        "measurements": {},
    }

    meas = cell_doc["measurements"]
    meas["reproduction"] = res.get("reproduction")
    meas["control_status"] = res.get("status", "ok")
    meas["count60"] = res.get("count60")
    meas["direction"] = res.get("direction")
    meas["delta_vs_record"] = res.get("delta_vs_record")

    # ---- stored witness lengths, from the record and root isolation
    wp = [float(s) for s in entry["witness_P"]]
    wq = [float(s) for s in entry["witness_Q"]]
    while len(wp) > 1 and wp[-1] == 0.0:
        wp.pop()
    while len(wq) > 1 and wq[-1] == 0.0:
        wq.pop()
    ri = next((f for f in rootiso["fits"]
               if f["source"].startswith("results-corrected")
               and f["index"] == index), None)
    meas["stored_witness"] = {
        "P_length": len(entry["witness_P"]),
        "Q_length": len(entry["witness_Q"]),
        "P_trimmed_length": len(wp),
        "Q_trimmed_length": len(wq),
        "full_evaluated_count": len(wp) + len(wq) - 1,
        "root_isolation": (None if ri is None else {
            "P_stored_length": ri["P_stored_length"],
            "Q_stored_length": ri["Q_stored_length"],
            "P_nominal_length": ri["P_nominal_length"],
            "Q_nominal_length": ri["Q_nominal_length"],
            "P_extra_count": ri["P_extra_count"],
            "Q_extra_count": ri["Q_extra_count"],
            "P_extra_max": ri["P_extra_max"],
            "Q_extra_max": ri["Q_extra_max"],
            "stored_count": ri["stored"],
            "stored_count_match": ri["stored_count_match"],
        }),
    }

    # ---- recorded fit: the re-solved record pair from the control's
    # confirm row (the control's own reproduction solve).
    confirm = None
    for lv in res.get("levels", []):
        if lv.get("label") == "confirm-record-pair":
            confirm = lv["rows"][0]
            break
    if (confirm is None or "P70" not in confirm) and (
            res.get("status") == "reproduction-failed"):
        # The control returns before its confirm row is serialized; re-solve
        # the record pair on the route the control tried, so the decisive
        # full-versus-nominal comparison still exists.
        rep = res.get("reproduction") or {}
        settings = (ctl.DEFAULT_SETTINGS if rep.get("settings") != "tightened"
                    else ctl.TIGHT_SETTINGS)
        ts0, xs0, fs0 = grid(a, W, 2000)
        confirm = ctl.solve_one(a, W, kind, target_value, m, k, settings,
                                ts0, xs0, fs0)
    if confirm is not None and "P70" in confirm:
        P, Q = pq_from_row(confirm)
        fb, fberr = try_call(fit_block, P, Q, m, k, a, W, kind,
                             target_value)
        if fberr is not None:
            meas["recorded_fit"] = None
            meas["recorded_fit_error"] = fberr
        else:
            meas["recorded_fit"] = {
                "pair": [m, k],
                "settings": confirm.get("settings"),
                "control_row": {
                    key2: confirm.get(key2) for key2 in
                    ("e60", "e60_rel", "e60_abs", "e60_at_t", "e64r",
                     "e64r_rel", "e64r_abs", "e64", "e64_rel", "e64_abs",
                     "guard60", "guard64", "meets60", "converged", "deviation")
                },
                **fb,
                "resolved_P_length": len(ptrim_mpf(P)),
                "resolved_Q_length": len(ptrim_mpf(Q)),
                "length_match_root_isolation": (
                    ri is not None
                    and len(ptrim_mpf(P)) == ri["P_stored_length"]
                    and len(ptrim_mpf(Q)) == ri["Q_stored_length"]),
                "delivered_2001_match": (
                    confirm.get("e64") == repr(entry["delivered_2001"])),
            }
            # the full-witness binary64 value in the target kind must
            # reproduce the record's delivered value digit for digit
            meas["recorded_fit"]["reproduced_e64_target"] = (
                confirm.get("e64") if kind == "rel" else confirm.get("e64_abs"))
            meas["delivered_rel"] = confirm.get("e64_rel")
            meas["delivered_2001_of_record"] = entry["delivered_2001"]
            meas["delivered_kind"] = (
                "rel" if kind == "rel" else "abs (the record's target kind)")
            meas["guard_2001_of_record"] = entry.get("guard_2001")
            meas["guard_2001_match"] = (
                entry.get("guard_2001") is None
                or confirm.get("guard64") == repr(entry["guard_2001"]))
    else:
        meas["recorded_fit"] = None

    # ---- Chebyshev count for this cell
    width_map = {
        "m7": "0.0078125", "m6": "0.015625", "m5": "0.03125",
        "m4": "0.0625", "m3": "0.125", "m2": "0.25", "m1": "0.5",
        "p0": "1.0", "p1": "2.0", "p2": "4.0", "p3": "8.0",
        "p4": "16.0", "p5": "32.0", "ship": "ship", "paper": "paper",
    }
    hits = []
    for c in cheb["cells"]:
        if c["mode"] != "records" or c["target"] != target:
            continue
        if width_name in ("ship", "paper"):
            if width_name in c["width_names"]:
                hits.append(c)
        elif abs(float(c["width"]) - float(width_map[width_name])) < 1e-12:
            hits.append(c)
    if len(hits) != 1:
        raise RuntimeError(f"chebyshev join for {width_name} {target}: "
                           f"{len(hits)} hits")
    ch = hits[0]
    meas["chebyshev"] = {
        "source": "coefficient-count-sweep/records/chebyshev_recount.json",
        "joined_by": "width name and target",
        "width_record": ch["width"],
        "width_names": ch["width_names"],
        "target": ch["target"],
        "count_old": ch["old_stored"],
        "count_new": ch["new_stored"],
        "new_m": ch.get("new_m"),
        "delivered_rel": ch.get("delivered_rel"),
    }

    # ---- the 60-digit control's accepted fit
    acc = res.get("accepted")
    if acc is not None and "P70" in acc:
        Pa, Qa = pq_from_row(acc)
        am, ak = acc["m"], acc["k"]
        fb, fberr = try_call(fit_block, Pa, Qa, am, ak, a, W, kind,
                             target_value)
        if fberr is not None:
            meas["control_accepted"] = None
            meas["control_accepted_error"] = fberr
        else:
            meas["control_accepted"] = {
                "pair": [am, ak],
                "nominal_count": (am + 1) + ak,
                "settings": acc.get("settings"),
                "converged": acc.get("converged"),
                "deviation": acc.get("deviation"),
                "scores": {
                    key2: acc.get(key2) for key2 in
                    ("e60", "e60_rel", "e60_abs", "e60_at_t", "e64r",
                     "e64r_rel", "e64r_abs", "e64", "e64_rel", "e64_abs",
                     "guard60", "guard64", "meets60",
                     "e60_4001", "e60_4001_rel", "e60_4001_abs",
                     "guard60_4001", "meets60_4001")
                },
                **fb,
                "genuine_nominal_at_its_count": {
                    "nominal_meets60_2001":
                        fb["nominal"]["grid_2001"]["meets60"],
                    "nominal_meets64_2001":
                        fb["nominal"]["grid_2001"]["meets64"],
                    "nominal_meets60_4001":
                        fb["nominal"]["grid_4001"]["meets60"],
                    "nominal_meets64_4001":
                        fb["nominal"]["grid_4001"]["meets64"],
                },
                "resolved_P_length": len(ptrim_mpf(Pa)),
                "resolved_Q_length": len(ptrim_mpf(Qa)),
            }
    else:
        meas["control_accepted"] = None

    # ---- all pairs solved at count60: is any of them a genuine nominal fit?
    # Direction "down"/"up" cells were fully swept at count60 by the control.
    # Direction "equal" cells confirmed only the record pair; solve the whole
    # level here on the same settings route, default first and tightened only
    # if the default level has no meeting fit, exactly the control's fallback.
    count60 = res.get("count60")
    rows60 = []
    extra_levels = []
    if count60 is not None:
        for lv in res.get("levels", []):
            if lv.get("stored") == count60 and lv.get("label"):
                rows60.extend(lv["rows"])
        full_sweep_labels = ("below-full-sweep", "upward-full-sweep",
                             "record-count-full-sweep")
        swept = any(lv.get("stored") == count60
                    and any(lv.get("label", "").startswith(p)
                            for p in full_sweep_labels)
                    for lv in res.get("levels", []))
        if not swept and count60 is not None:
            cell = res["cell"]
            ts, xs, fs = grid(a, W, 2000)
            lv = ctl.level_sweep(cell, count60, ctl.DEFAULT_SETTINGS,
                                 ts, xs, fs, verbose=False)
            lv["label"] = "genuine-nominal-level-sweep-default"
            rows60.extend(lv["rows"])
            extra_levels.append(ctl.serialize_levels([lv])[0])
            if not any(r.get("meets60") for r in lv["rows"]):
                lv2 = ctl.level_sweep(cell, count60, ctl.TIGHT_SETTINGS,
                                      ts, xs, fs, verbose=False)
                lv2["label"] = "genuine-nominal-level-sweep-tightened"
                rows60.extend(lv2["rows"])
                extra_levels.append(ctl.serialize_levels([lv2])[0])
        checks = nominal_checks(rows60, a, W, kind, target_value)
        meas["count60_level"] = {
            "count60": count60,
            "swept_by_control": bool(swept),
            "extra_level_run": bool(extra_levels),
            "rows_solved": len(rows60),
            "checks": checks,
            "any_nominal_meets60_2001": any(
                c["nominal_meets60_2001"] for c in checks),
            "any_nominal_meets64_2001": any(
                c["nominal_meets64_2001"] for c in checks),
            "any_nominal_meets60_4001": any(
                c["nominal_meets60_4001"] for c in checks),
            "any_nominal_meets64_4001": any(
                c["nominal_meets64_4001"] for c in checks),
            "genuine_pairs_60": [c["pair"] for c in checks
                                 if c["nominal_meets60_2001"]],
            "genuine_pairs_64": [c["pair"] for c in checks
                                 if c["nominal_meets64_2001"]],
        }
        if extra_levels:
            res.setdefault("levels", []).extend(extra_levels)

    slim_control(res)
    cell_doc["control"] = res

    cell_doc["meta"] = {
        "generated_utc": now(),
        "script_sha256": sha256(Path(__file__).resolve()),
        "control_script_sha256": sha256(CTRL_PATH),
        "corrected_record_sha256": sha256(
            ROOT / "coefficient-count-sweep" / "results-corrected"
            / "corrected-counts.json"),
        "chebyshev_record_sha256": sha256(
            ROOT / "coefficient-count-sweep" / "records"
            / "chebyshev_recount.json"),
        "root_isolation_sha256": sha256(
            ROOT / "coefficient-count-sweep" / "strong-optional"
            / "root-isolation" / "records" / "root_isolation.json"),
        "versions": {"python": sys.version.split()[0],
                     "mpmath": mp.__version__,
                     "numpy": __import__("numpy").__version__,
                     "baryrat": getattr(__import__("baryrat"),
                                        "__version__", "unknown")},
        "convert_dps": CONVERT_DPS, "evaluation_dps": EVAL_DPS,
        "den_floor": DEN_FLOOR, "k_max": ctl.K_MAX,
        "grids": {"certifying_nodes": 2001, "check_nodes": 4001},
    }

    SWEEP.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(cell_doc, indent=1))
    print(f"[{key}] wrote {out_path} count60={count60} "
          f"status={meas['control_status']} seconds={control_seconds}",
          flush=True)
    return cell_doc


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--index", type=int, required=True)
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args(argv)
    run_cell(args.index, force=args.force)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
