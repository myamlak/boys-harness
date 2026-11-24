"""One clean-edge cell under the paper C 60-digit evaluation control.

Re-solves one (width, target) cell of the clean-edge record at a = 62 with the
paper C control's numeric core, solver settings and 60-digit scoring, and
records the 60-digit count.  This is the per-cell half of the full forty-row
control; `run_all_60dig.py` starts one process per cell and combines the
records.

Protocol, unchanged from the pilot and from
results-control-60dig/control_60dig_brasil.py:
  * grid 2001 nodes, x = a + W (t + 1)/2, F_0 reference at 70 decimal digits;
  * BRASIL (baryrat 2.1.2), default settings, tightened only where the default
    route finds no meeting fit; k <= 12; denominator guard 1e-3;
  * the record's pair and count are re-solved, then every admissible pair at
    each lower stored count is swept, default first and tightened only for a
    count whose default sweep found no meeting fit; the first count where
    neither route meets is the below-count check;
  * a fit meets at 60 digits when the monomial coefficients held at 60
    significant digits and evaluated at 60 dps are within the target relative
    to F_0 at each node and pass the denominator guard.

Usage:
  python run_cell_60dig.py --width "2^5" --target "rel:1e-15" --out cell.json
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import math
import sys
import time
import traceback
from datetime import datetime, timezone
from pathlib import Path

CTRL = Path(r"C:\Users\myaml\dev\papers-c-v5\papers\coefficient-counts"
            r"\results-control-60dig\control_60dig_brasil.py")
REC = Path(r"C:\Users\myaml\dev\papers-harness-repo\coefficient-count-sweep"
           r"\clean-edge-2026-11-22\edge62_sweep.json")
SCHEMA = "clean-edge-60dig-cell-1"

TARGETS = ["rel:1e-08", "rel:1e-09", "rel:1e-10", "rel:1e-11",
           "rel:1e-12", "rel:1e-13", "rel:1e-14", "rel:1e-15"]
EXPONENTS = [8, 9, 10, 11, 12, 13, 14, 15]
CONTRACT = 6

PRECISION = {
    "evaluation": "monomial coefficients held at 60 significant digits and "
                  "evaluated at 60 decimal digits",
    "conversion_dps": 70,
    "reference": "F_0 at 70 decimal digits (mpmath)",
    "solver": "BRASIL (baryrat 2.1.2), binary64 fit",
    "grid_nodes": 2001,
    "check_nodes": 4001,
    "k_max": 12,
    "den_floor": 1e-3,
    "target": "relative to F_0 at each node",
}


def load_ctrl():
    spec = importlib.util.spec_from_file_location("c60", CTRL)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["c60"] = mod
    spec.loader.exec_module(mod)
    return mod


# --- clean-edge rule, copied from analyze_edge62.py (standard library only) ---

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


def slim_levels(levels):
    """Pair-level rows without the heavy coefficient/nodal payload."""
    fields = ("m", "k", "stored", "settings", "converged", "deviation",
              "info_error", "e60", "e60_abs", "e64r", "e64", "guard60",
              "meets60", "seconds", "exception", "convert_exception",
              "score_exception")
    out = []
    for lv in levels:
        rows = [{k: r.get(k) for k in fields} for r in lv["rows"]]
        out.append({"stored": lv["stored"], "settings": lv["settings"],
                    "label": lv.get("label"), "meets60": lv["meets60"],
                    "solved": lv["solved"],
                    "pairs_tried": lv["pairs_tried"],
                    "seconds": lv["seconds"], "rows": rows})
    return out


def brasil_stats(levels):
    solves = unconverged = exceptions = 0
    for lv in levels:
        for r in lv["rows"]:
            if "converged" in r:
                solves += 1
                if r.get("converged") is False:
                    unconverged += 1
            if any(r.get(k) for k in ("exception", "convert_exception",
                                      "score_exception")):
                exceptions += 1
    return {"solves": solves, "unconverged": unconverged,
            "exceptions": exceptions}


def run_cell(width_name, target_tag, record=REC):
    c60 = load_ctrl()
    rec_doc = json.loads(Path(record).read_text(encoding="utf-8"))
    a = float(rec_doc["meta"]["edge"])
    rows = rec_doc["rows"]
    row = next(r for r in rows
               if r["width_name"] == width_name and r["target"] == target_tag)
    width = float(row["width_f"])
    target = float(row["target_value"])
    kind = row["target_kind"]
    rec_stored = int(row["rat"]["stored"])
    rec_pair = (int(row["rat"]["m"]), int(row["rat"]["k"]))
    rec_solver = row["rat"]["solver"]
    rec_delivered = row["rat"]["delivered_rel"]
    cell = {"key": f"{width_name}:{target_tag}", "width_name": width_name,
            "width": width, "a": a, "target": target_tag, "kind": kind,
            "target_value": target, "record_count": rec_stored,
            "record_pair": rec_pair, "record_solver": rec_solver,
            "record_delivered_rel": rec_delivered}
    print(f"[{cell['key']}] a={a} W={width!r} record {rec_stored} "
          f"[{rec_pair[0]}/{rec_pair[1]}] solver={rec_solver} "
          f"delivered={rec_delivered}", flush=True)

    c60.mp.mp.dps = c60.CONVERT_DPS
    ts, xs, fs = c60.build_grid(a, width)

    levels = []

    def cell_proxy():
        return {"a": a, "width": width, "kind": kind,
                "target_value": target}

    def sweep(stored, settings, label):
        lv = c60.level_sweep(cell_proxy(), stored, settings, ts, xs, fs,
                             verbose=True)
        lv["label"] = label
        levels.append(lv)
        return lv

    # confirm the record's pair, re-solved with BRASIL (informational only:
    # the record's fit was produced by the kernel Remez/Lawson solver)
    started = time.perf_counter()
    rep = c60.solve_one(a, width, kind, target, rec_pair[0], rec_pair[1],
                        c60.DEFAULT_SETTINGS, ts, xs, fs)
    rep_route = "default"
    if not c60.meets(rep):
        rep2 = c60.solve_one(a, width, kind, target, rec_pair[0], rec_pair[1],
                             c60.TIGHT_SETTINGS, ts, xs, fs)
        if c60.meets(rep2):
            rep, rep_route = rep2, "tightened"
    print(f"  record pair re-solve ({rep_route}): converged={rep.get('converged')} "
          f"dev={rep.get('deviation')} e60={rep.get('e60')} "
          f"e64={rep.get('e64')} meets60={c60.meets(rep)} "
          f"{time.perf_counter() - started:.1f}s", flush=True)
    confirm = {"stored": rec_stored, "settings": rep_route,
               "pairs_tried": [f"[{rec_pair[0]}/{rec_pair[1]}]"],
               "solved": 1, "meets60": c60.meets(rep), "rows": [rep],
               "seconds": None, "label": "confirm-record-pair"}
    levels.append(confirm)

    final = None
    direction = None
    accepted_pool = []
    fail = None
    minimality = None
    if c60.meets(rep):
        final, direction, accepted_pool = rec_stored, "equal", [rep]

    if final is None:
        lv = sweep(rec_stored, c60.DEFAULT_SETTINGS,
                   "record-count-full-sweep-default")
        pool = [r for r in lv["rows"] if c60.meets(r)]
        if pool:
            final, direction, accepted_pool = rec_stored, "equal-full-sweep", pool
        else:
            lv2 = sweep(rec_stored, c60.TIGHT_SETTINGS,
                        "record-count-full-sweep-tightened")
            pool = [r for r in lv2["rows"] if c60.meets(r)]
            if pool:
                final, direction, accepted_pool = rec_stored, "equal-full-sweep", pool

    if final is None:
        for stored in range(rec_stored + 1, c60.MAX_STORED + 1):
            lv = sweep(stored, c60.DEFAULT_SETTINGS, "upward-full-sweep-default")
            pool = [r for r in lv["rows"] if c60.meets(r)]
            if pool:
                final, direction, accepted_pool = stored, "up", pool
                break
            lv2 = sweep(stored, c60.TIGHT_SETTINGS, "upward-full-sweep-tightened")
            pool = [r for r in lv2["rows"] if c60.meets(r)]
            if pool:
                final, direction, accepted_pool = stored, "up", pool
                break
        if final is not None:
            below = final - 1
            failed = [lv for lv in levels
                      if lv["stored"] == below and not lv["meets60"]
                      and lv["label"] != "confirm-record-pair"]
            fail = c60.summarize_failure(below, failed)
            minimality = "upward walk from the record count"

    if final is not None and direction != "up":
        minimality = ("every admissible pair at each count below was solved on "
                      "the default route; a count whose default route found no "
                      "meeting fit was also swept on the tightened route; the "
                      "first count where neither route meets is the "
                      "below-count check")
        stored = final
        while stored > 3:
            below = stored - 1
            lv = sweep(below, c60.DEFAULT_SETTINGS, "below-full-sweep-default")
            pool = [r for r in lv["rows"] if c60.meets(r)]
            if pool:
                final, direction, accepted_pool = below, "down", pool
                stored = below
                continue
            lv2 = sweep(below, c60.TIGHT_SETTINGS, "below-full-sweep-tightened")
            pool = [r for r in lv2["rows"] if c60.meets(r)]
            if pool:
                final, direction, accepted_pool = below, "down", pool
                stored = below
                continue
            fail = c60.summarize_failure(below, [lv, lv2])
            break

    accepted = c60.best_row(accepted_pool) if accepted_pool else None
    acc = None
    if accepted is not None:
        ts4, xs4, fs4 = c60.build_grid(a, width, n=4000)
        check = c60.score_4001(accepted["_P"], accepted["_Q"], ts4, fs4, target)
        acc = {
            "m": accepted["m"], "k": accepted["k"],
            "stored": accepted["stored"], "settings": accepted["settings"],
            "converged": accepted.get("converged"),
            "deviation": accepted.get("deviation"),
            "e60": accepted.get("e60"), "e60_abs": accepted.get("e60_abs"),
            "e60_at_t": accepted.get("e60_at_t"),
            "e64r": accepted.get("e64r"), "e64r_abs": accepted.get("e64r_abs"),
            "e64": accepted.get("e64"), "e64_abs": accepted.get("e64_abs"),
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

    series = [float(next(r for r in rows
                         if r["width_name"] == width_name and r["target"] == t)
                    ["rat"]["stored"]) for t in TARGETS]
    sub = list(series)
    if final is not None:
        sub[TARGETS.index(target_tag)] = float(final)
    freeline = {
        "record_series": series, "substituted_series": sub,
        "record_primary": metrics(series, len(TARGETS)),
        "substituted_primary": metrics(sub, len(TARGETS)),
        "record_secondary": metrics(series, CONTRACT),
        "substituted_secondary": metrics(sub, CONTRACT),
    }
    for tag in ("primary", "secondary"):
        freeline[f"record_{tag}_verdict"] = rat_verdict(freeline[f"record_{tag}"])
        freeline[f"substituted_{tag}_verdict"] = rat_verdict(
            freeline[f"substituted_{tag}"])

    return {
        "schema": SCHEMA,
        "cell": cell,
        "precision": PRECISION,
        "reproduction": {
            "record_solver": rec_solver,
            "brasil_route": rep_route,
            "brasil_converged": rep.get("converged"),
            "brasil_deviation": rep.get("deviation"),
            "brasil_e64": rep.get("e64"),
            "record_delivered_rel": rec_delivered,
            "exact_match": rep.get("e64") == repr(float(rec_delivered)),
        },
        "count60": final,
        "direction": direction,
        "minimality": minimality,
        "delta_vs_record": None if final is None else final - rec_stored,
        "accepted": acc,
        "failing_level": fail,
        "freeline": freeline,
        "brasil": brasil_stats(levels),
        "levels": slim_levels(levels),
    }


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--width", required=True)
    ap.add_argument("--target", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--record", default=str(REC))
    args = ap.parse_args(argv)
    out = Path(args.out)
    started = time.perf_counter()
    try:
        result = run_cell(args.width, args.target, args.record)
        result["status"] = "ok"
    except Exception as ex:  # noqa: BLE001 - record, never lose the cell
        result = {
            "schema": SCHEMA,
            "status": "error",
            "cell": {"width_name": args.width, "target": args.target},
            "precision": PRECISION,
            "error": repr(ex),
            "traceback": traceback.format_exc(),
        }
    result["wall_seconds"] = round(time.perf_counter() - started, 1)
    result["generated_utc"] = datetime.now(timezone.utc).strftime(
        "%Y-%m-%dT%H:%M:%SZ")
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=1), encoding="utf-8")
    acc = result.get("accepted") or {}
    print(f"RESULT {args.width} {args.target} status={result['status']} "
          f"count60={result.get('count60')} record="
          f"{(result.get('cell') or {}).get('record_count')} "
          f"pair=[{acc.get('m')}/{acc.get('k')}] e60={acc.get('e60')} "
          f"guard60={acc.get('guard60')} e64={acc.get('e64')} "
          f"converged={acc.get('converged')} "
          f"wall={result['wall_seconds']}s", flush=True)
    print(f"  wrote {out}", flush=True)
    return 0 if result["status"] == "ok" else 1


if __name__ == "__main__":
    raise SystemExit(main())
