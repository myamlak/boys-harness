"""Run the full forty-row clean-edge 60-digit control in parallel.

Starts one independent process per (width, target) cell of the clean-edge
record, five to eight at a time, using `run_cell_60dig.py`.  Every cell gets
its own JSON and log; the launcher writes one combined record with the
required per-cell fields, the per-width free-line table on the 60-digit counts
and the SHA-256 identity of the inputs.

Usage:
  python run_all_60dig.py [--jobs 8] [--timeout 3600]
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
RUNNER = HERE / "run_cell_60dig.py"
CELLS = HERE / "cells"
LOGS = HERE / "logs"
COMBINED = HERE / "combined_60dig.json"
REC = Path(r"C:\Users\myaml\dev\papers-harness-repo\coefficient-count-sweep"
           r"\clean-edge-2026-11-22\edge62_sweep.json")
CTRL = Path(r"C:\Users\myaml\dev\papers-c-v5\papers\coefficient-counts"
            r"\results-control-60dig\control_60dig_brasil.py")

WIDTHS = ["2^-7", "2^0", "2^5", "ship", "paper"]
TARGETS = ["rel:1e-08", "rel:1e-09", "rel:1e-10", "rel:1e-11",
           "rel:1e-12", "rel:1e-13", "rel:1e-14", "rel:1e-15"]
EXPONENTS = [8, 9, 10, 11, 12, 13, 14, 15]
CONTRACT = 6
WIDTH_TOKEN = {"2^-7": "m7", "2^0": "p0", "2^5": "p5",
               "ship": "ship", "paper": "paper"}
TARGET_TOKEN = {t: t.replace(":", "_") for t in TARGETS}
WCOST = {"2^5": 1.0, "ship": 0.85, "paper": 0.8, "2^0": 0.4, "2^-7": 0.15}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


# --- clean-edge rule, copied from analyze_edge62.py ---

def fit_line(xs, ys):
    n = len(xs)
    xbar = math.fsum(xs) / n
    ybar = math.fsum(ys) / n
    sxx = math.fsum((x - xbar) ** 2 for x in xs)
    sxy = math.fsum((x - xbar) * (y - ybar) for x, y in zip(xs, ys))
    slope = sxy / sxx
    intercept = ybar - slope * xbar
    return intercept, slope


def metrics(ys, count):
    xs = EXPONENTS[:count]
    ys = ys[:count]
    intercept, slope = fit_line(xs, ys)
    residual = max(abs(y - (intercept + slope * x)) for x, y in zip(xs, ys))
    second = max((abs(ys[i - 1] - 2 * ys[i] + ys[i + 1])
                  for i in range(1, len(ys) - 1)), default=0.0)
    return {"n": count, "intercept": intercept, "slope": slope,
            "max_abs_residual": residual,
            "max_abs_second_difference": second}


def rat_verdict(stats):
    if stats["max_abs_residual"] > 2.0 or stats["max_abs_second_difference"] >= 2.0:
        return "fails"
    return "holds"


def cell_key(width, target):
    # Windows-safe: no ':' or '^'
    return f"{WIDTH_TOKEN[width]}_{TARGET_TOKEN[target]}"


def cost(task):
    stored = task["record_count"]
    pairs = max(0, min(12, stored - 2))
    depth = 1.0 + 0.3 * (TARGETS.index(task["target"]) - 4) / 7.0
    return WCOST[task["width"]] * max(pairs, 1) * depth


def build_tasks(rec_doc):
    tasks = []
    for width in WIDTHS:
        for target in TARGETS:
            row = next(r for r in rec_doc["rows"]
                       if r["width_name"] == width and r["target"] == target)
            tasks.append({
                "key": cell_key(width, target),
                "width": width, "width_name": width, "target": target,
                "width_f": float(row["width_f"]),
                "target_value": float(row["target_value"]),
                "kind": row["target_kind"],
                "record_count": int(row["rat"]["stored"]),
                "record_pair": [int(row["rat"]["m"]), int(row["rat"]["k"])],
                "record_solver": row["rat"]["solver"],
            })
    tasks.sort(key=cost, reverse=True)
    return tasks


def run_one(task, timeout):
    out = CELLS / f"{task['key']}.json"
    log = LOGS / f"{task['key']}.log"
    cmd = [sys.executable, str(RUNNER),
           "--width", task["width"], "--target", task["target"],
           "--out", str(out), "--record", str(REC)]
    env = None
    started = time.perf_counter()
    with open(log, "w", encoding="utf-8") as fh:
        try:
            proc = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT,
                                  timeout=timeout)
            code = proc.returncode
            status = "ok" if code == 0 else "error"
            note = None
        except subprocess.TimeoutExpired:
            code = "timeout"
            status = "timeout"
            note = f"exceeded {timeout} s and was killed"
    wall = round(time.perf_counter() - started, 1)
    if out.exists():
        doc = json.loads(out.read_text(encoding="utf-8"))
        doc["launcher_status"] = doc.get("status", status)
        doc["launcher_note"] = note
        doc["launcher_wall_seconds"] = wall
        doc["launcher_exit_code"] = code
        doc["command"] = cmd
    else:
        doc = {"schema": "clean-edge-60dig-cell-1", "status": status,
               "launcher_note": note, "launcher_wall_seconds": wall,
               "launcher_exit_code": code, "command": cmd,
               "cell": {"width_name": task["width"],
                        "target": task["target"]},
               "log_tail": log.read_text(encoding="utf-8",
                                         errors="replace")[-2000:]}
    return task, doc


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--jobs", type=int, default=8)
    ap.add_argument("--timeout", type=float, default=3600.0)
    args = ap.parse_args(argv)
    if not 5 <= args.jobs <= 8:
        raise SystemExit("--jobs must be between 5 and 8")
    CELLS.mkdir(exist_ok=True)
    LOGS.mkdir(exist_ok=True)
    rec_doc = json.loads(REC.read_text(encoding="utf-8"))
    tasks = build_tasks(rec_doc)
    print(f"forty-row clean-edge 60-digit control: {len(tasks)} cells, "
          f"{args.jobs} processes, timeout {args.timeout:g} s/cell", flush=True)
    started_all = time.perf_counter()
    results = {}
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futures = {pool.submit(run_one, t, args.timeout): t for t in tasks}
        for fut in as_completed(futures):
            task, doc = fut.result()
            results[task["key"]] = doc
            acc = doc.get("accepted") or {}
            print(f"DONE {task['key']:18s} status={doc.get('status')} "
                  f"count60={doc.get('count60')} record="
                  f"{task['record_count']} wall={doc.get('launcher_wall_seconds')}s "
                  f"e60={acc.get('e60')} conv={acc.get('converged')}",
                  flush=True)
    total_wall = round(time.perf_counter() - started_all, 1)

    # ---- combined per-cell records, canonical order
    cells = []
    for width in WIDTHS:
        for target in TARGETS:
            key = cell_key(width, target)
            doc = results[key]
            cell = doc.get("cell") or {}
            acc = doc.get("accepted") or {}
            task = next(t for t in tasks if t["key"] == key)
            cells.append({
                "key": key,
                "width": task["width_f"],
                "width_name": width,
                "target": target,
                "record_count": task["record_count"],
                "record_pair": task["record_pair"],
                "record_solver": task["record_solver"],
                "count60": doc.get("count60"),
                "direction": doc.get("direction"),
                "delta_vs_record": doc.get("delta_vs_record"),
                "accepted_pair": ([acc.get("m"), acc.get("k")]
                                  if acc else None),
                "accepted_stored": acc.get("stored"),
                "accepted_settings": acc.get("settings"),
                "e60": acc.get("e60"),
                "guard60": acc.get("guard60"),
                "e64": acc.get("e64"),
                "e64r": acc.get("e64r"),
                "converged": acc.get("converged"),
                "precision": doc.get("precision"),
                "brasil": doc.get("brasil"),
                "status": doc.get("status"),
                "wall_seconds": doc.get("wall_seconds"),
                "launcher_wall_seconds": doc.get("launcher_wall_seconds"),
                "error": doc.get("error"),
                "launcher_note": doc.get("launcher_note"),
                "command": doc.get("command"),
            })

    # ---- per-width free-line table on the 60-digit counts
    freeline = {}
    changed = []
    partial = []
    for width in WIDTHS:
        rec_series = []
        sub_series = []
        have_all = True
        for target in TARGETS:
            c = next(c for c in cells
                     if c["width_name"] == width and c["target"] == target)
            rec_series.append(float(c["record_count"]))
            if c["count60"] is not None and c["status"] == "ok":
                sub_series.append(float(c["count60"]))
            else:
                sub_series.append(float(c["record_count"]))
                have_all = False
        rec_primary = metrics(rec_series, len(TARGETS))
        sub_primary = metrics(sub_series, len(TARGETS))
        rec_secondary = metrics(rec_series, CONTRACT)
        sub_secondary = metrics(sub_series, CONTRACT)
        rp, sp = rat_verdict(rec_primary), rat_verdict(sub_primary)
        rs, ss = rat_verdict(rec_secondary), rat_verdict(sub_secondary)
        entry = {
            "record_series": rec_series,
            "substituted_series": sub_series,
            "complete": have_all,
            "record_primary": rec_primary,
            "substituted_primary": sub_primary,
            "record_primary_verdict": rp,
            "substituted_primary_verdict": sp,
            "record_secondary": rec_secondary,
            "substituted_secondary": sub_secondary,
            "record_secondary_verdict": rs,
            "substituted_secondary_verdict": ss,
            "changed_primary": rp != sp,
            "changed_secondary": rs != ss,
        }
        freeline[width] = entry
        if rp != sp:
            changed.append(width)
        if not have_all:
            partial.append(width)

    rat_fail_record = [w for w in WIDTHS
                       if freeline[w]["record_primary_verdict"] == "fails"]
    rat_fail_sub = [w for w in WIDTHS
                    if freeline[w]["substituted_primary_verdict"] == "fails"]
    cheb_ok = True  # Chebyshev arm is the binary64 record, unchanged
    if not cheb_ok:
        verdict = "uninterpretable-chebyshev"
    elif len(rat_fail_sub) >= 3:
        verdict = "clean-edge-failure-confirmed"
    elif len([w for w in WIDTHS
              if freeline[w]["substituted_primary_verdict"] == "holds"]) >= 4:
        verdict = "transient-specific"
    else:
        verdict = "mixed"
    rec_verdict = ("clean-edge-failure-confirmed" if len(rat_fail_record) >= 3
                   else "transient-specific" if (5 - len(rat_fail_record)) >= 4
                   else "mixed")

    summary = {
        "cells_total": len(cells),
        "cells_ok": sum(1 for c in cells if c["status"] == "ok"),
        "cells_error": sum(1 for c in cells if c["status"] == "error"),
        "cells_timeout": sum(1 for c in cells if c["status"] == "timeout"),
        "cells_no_count60": sum(1 for c in cells if c["count60"] is None),
        "brasil_solves": sum((c.get("brasil") or {}).get("solves", 0)
                             for c in cells),
        "brasil_unconverged": sum((c.get("brasil") or {}).get("unconverged", 0)
                                  for c in cells),
        "brasil_exceptions": sum((c.get("brasil") or {}).get("exceptions", 0)
                                 for c in cells),
        "changed_widths": changed,
        "partial_widths": partial,
        "record_primary_verdict": rec_verdict,
        "substituted_primary_verdict": verdict,
        "record_rational_fail_widths": rat_fail_record,
        "substituted_rational_fail_widths": rat_fail_sub,
    }

    combined = {
        "schema": "clean-edge-60dig-full-1",
        "generated_utc": datetime.now(timezone.utc).strftime(
            "%Y-%m-%dT%H:%M:%SZ"),
        "record": {"path": str(REC), "sha256": sha256(REC)},
        "control_core": {"path": str(CTRL), "sha256": sha256(CTRL)},
        "runner": {"path": str(RUNNER), "sha256": sha256(RUNNER)},
        "launcher": {"path": str(Path(__file__).resolve()),
                     "sha256": sha256(Path(__file__).resolve())},
        "jobs": args.jobs,
        "cell_timeout_seconds": args.timeout,
        "total_wall_seconds": total_wall,
        "settings": {
            "edge": rec_doc["meta"]["edge"],
            "scale": rec_doc["meta"]["scale"],
            "targets": TARGETS,
            "widths": WIDTHS,
            "grid_nodes": 2001,
            "check_nodes": 4001,
            "k_max": 12,
            "den_floor": 1e-3,
            "convert_dps": 70,
            "eval_dps": 60,
            "solver": "BRASIL (baryrat 2.1.2), binary64 fit",
            "default": {"tol": 1e-4, "maxiter": 1000, "npi": -30,
                        "init_steps": 100, "max_step_size": 0.1},
            "tightened": {"tol": 1e-8, "maxiter": 5000, "npi": -40,
                          "init_steps": 300, "max_step_size": 0.05},
        },
        "cells": cells,
        "freeline": freeline,
        "summary": summary,
    }
    COMBINED.write_text(json.dumps(combined, indent=1, default=str),
                        encoding="utf-8")
    print("\nper-width free-line table (primary, all eight targets):",
          flush=True)
    print(f"{'width':6s} {'record r/d2':>16s} {'verdict':>8s} "
          f"{'60-digit r/d2':>16s} {'verdict':>8s} {'changed':>7s} "
          f"{'complete':>8s}", flush=True)
    for w in WIDTHS:
        e = freeline[w]
        print(f"{w:6s} "
              f"{e['record_primary']['max_abs_residual']:7.3f}/"
              f"{e['record_primary']['max_abs_second_difference']:4.1f} "
              f"{e['record_primary_verdict']:>8s} "
              f"{e['substituted_primary']['max_abs_residual']:7.3f}/"
              f"{e['substituted_primary']['max_abs_second_difference']:4.1f} "
              f"{e['substituted_primary_verdict']:>8s} "
              f"{str(e['changed_primary']):>7s} {str(e['complete']):>8s}",
              flush=True)
    print(f"\nsummary: {json.dumps(summary)}", flush=True)
    print(f"total wall {total_wall} s", flush=True)
    print(f"combined SHA-256 {sha256(COMBINED)}", flush=True)
    print(f"wrote {COMBINED}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
