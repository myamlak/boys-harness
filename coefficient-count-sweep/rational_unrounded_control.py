#!/usr/bin/env python3
"""Referee control for the dense target ladder: count with unrounded evaluation.

The ladder's stored counts are produced by count_sweep.rational_count, which
solves every fit at 60 digits and then measures the delivered error with the
coefficients rounded to binary64 and evaluated by Horner in binary64
(count_sweep.rat_errors -> generator.rat_eval_double). This script keeps the
whole counting machinery and changes one thing: the delivered error of the fit
is evaluated at the solver's own working precision, 60 digits, on the same
certifying grid against the same 60-digit reference.

Everything else is the record's: the grid (2001 nodes uniform in t, endpoints
included), the solvers (generator.rat_remez with 60 iterations and
generator.rat_lawson with 300 outer iterations), the denominator cap k <= 12,
the target semantics (relative to F_0(x) at each node), and the definition of
the count (the smallest stored count (m+1)+k at which some pair meets).

The search is the paper's own below-check made into a walk. It starts at the
record's stored count and confirms that count with the recorded pair. If the
pair meets, every admissible pair is solved at each count below until none
meets; that failing sweep is the predicate rational_count applies when it
certifies its own count as minimal, and only the acceptance arithmetic changes.
If the record's pair misses, every pair at the record's count is swept and,
failing a meeting fit there, counts above the record's are searched until one
meets: the unrounded count can in principle exceed the recorded one, and that
outcome is recorded with its direction rather than reported as a smaller count.

Usage:
  python rational_unrounded_control.py --cells decisive
  python rational_unrounded_control.py --cells all
  python rational_unrounded_control.py --cells ship:1e-15,2^5:1e-15
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

PAPER = Path(__file__).resolve().parent
DEFAULT_BOYS = Path(r"C:\Users\myaml\dev\bladawiec\external\boys")
OUT_DEFAULT = PAPER / "results-ladder" / "rational_unrounded_control.json"
SWEEP = PAPER / "results-ladder" / "sweep_all.json"
RECERT = PAPER / "results-ladder" / "recertification.json"

K_MAX = 12
EVAL_DPS = 60
SCHEMA = "rational-unrounded-control-1"


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import {path}")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


# ---------------------------------------------------------------------------
# the one changed piece: the delivered error at the solver's own precision


def rat_errors_unrounded(g, p, q, a, b, ts, ref):
    """Worst relative and absolute error of p/q at 60 digits on the grid.

    The fit is a function of t; a and b are carried only so the signature
    matches the binary64 path. The grid's own t values are used, which is the
    same argument the solver fitted on and the same one the binary64 path
    recovers from x up to its own rounding. This is deliberately not
    ``cs.rat_errors`` with ``ts`` substituted: that function maps its argument
    from x to t in binary64 and must be given the x grid.
    """
    mp = g.mp
    worst_rel = g.mpf(0)
    worst_abs = g.mpf(0)
    at_rel = None
    at_abs = None
    with mp.workdps(EVAL_DPS):
        for t, r in zip(ts, ref):
            v = g.rat_value(p, q, t)
            e = abs(v - r)
            if e > worst_abs:
                worst_abs, at_abs = e, t
            rel = e / r
            if rel > worst_rel:
                worst_rel, at_rel = rel, t
    return worst_rel, worst_abs, at_rel, at_abs


def score(g, cs, a, b, ts, xs, ref, target, m, k, name, p, q):
    """Score one solved fit both ways on the same grid."""
    with g.mp.workdps(EVAL_DPS):
        wr60, wa60, at60, _ = rat_errors_unrounded(g, p, q, a, b, ts, ref)
        err64 = cs.rat_errors(g, p, q, a, b, xs, ref)
        wr64, wa64 = (err64[0], err64[1]) if err64 is not None else (None, None)
        return {
            "m": m, "k": k, "stored": (m + 1) + k, "solver": name,
            "rel60": g.mp.nstr(wr60, 17),
            "abs60": g.mp.nstr(wa60, 17),
            "at_t60": g.mp.nstr(at60, 17),
            "rel64": None if wr64 is None else g.mp.nstr(wr64, 17),
            "abs64": None if wa64 is None else g.mp.nstr(wa64, 17),
            "meets60": bool(cs.meets(target, wr60, wa60)),
            "meets64": bool(err64 is not None and cs.meets(target, *err64[:2])),
            "_wr60": wr60,
        }


def level_pass(g, cs, a, b, ts, xs, ref, stored, pairs, target, solver,
               verbose=False):
    """Solve every (m,k) at one stored count with one solver, and score it.

    The two solvers are separated so that the cheap minimax route can settle a
    level on its own. A level meets as soon as one fit meets, so the expensive
    Lawson route is run only when the Remez route has failed to find one and
    the level's negative answer is what has to be certified.
    """
    full_idx = list(range(len(ts)))
    started = time.perf_counter()
    rows = []
    if solver == "remez":
        for m, k in pairs:
            try:
                r = g.rat_remez(m, k, full_idx, ts, ref, None,
                                iters=cs.RAT_REMEZ_ITERS)
            except Exception:
                r = None
            if r is not None:
                rows.append(score(g, cs, a, b, ts, xs, ref, target, m, k,
                                  "remez", r[1], r[2]))
    else:
        ts_text = cs.ship_text(g, ts)
        ref_text = cs.ship_text(g, ref)
        specs = [(cs.BOYS_ROOT_STR, m, k, full_idx, ts_text, ref_text,
                  cs.RAT_REFINE_OUTER, cs.RAT_REMEZ_ITERS) for m, k in pairs]
        outs = g.run_jobs(cs.rat_candidate_job, specs)
        for (m, k), cands in zip(pairs, outs):
            for name, p, q in cands:
                if name == "lawson":
                    rows.append(score(g, cs, a, b, ts, xs, ref, target, m, k,
                                      name, p, q))
    solve_seconds = round(time.perf_counter() - started, 1)
    rows.sort(key=lambda r: r["_wr60"])
    meets = any(r["meets60"] for r in rows)
    out = {
        "stored": stored,
        "pass": solver,
        "pairs_tried": [f"[{m}/{k}]" for m, k in pairs],
        "solved": len(rows),
        "meets60": meets,
        "best_rel60": rows[0]["rel60"] if rows else None,
        "best_pair60": (f"[{rows[0]['m']}/{rows[0]['k']}] {rows[0]['solver']}"
                        if rows else None),
        "seconds": solve_seconds,
        "rows": [{k: v for k, v in r.items() if k != "_wr60"} for r in rows],
    }
    if verbose:
        print(f"  stored {stored} ({solver}): meets60={meets} "
              f"best={out['best_rel60']} ({solve_seconds} s)", flush=True)
    return out


def level_below(g, cs, a, b, ts, xs, ref, stored, target, verbose=False):
    """The full pair sweep at one count, Remez first and Lawson if needed."""
    pairs = [(stored - 1 - k, k)
             for k in range(1, min(K_MAX, stored - 2) + 1)]
    if not pairs:
        return None
    passes = [level_pass(g, cs, a, b, ts, xs, ref, stored, pairs, target,
                         "remez", verbose)]
    if not passes[0]["meets60"]:
        passes.append(level_pass(g, cs, a, b, ts, xs, ref, stored, pairs,
                                 target, "lawson", verbose))
    return {
        "stored": stored,
        "pairs_tried": passes[0]["pairs_tried"],
        "meets60": any(p["meets60"] for p in passes),
        "passes": passes,
    }


def accepted_row(levels):
    """The meeting fit at the smallest stored count, then the smallest error."""
    best = None
    for level in levels:
        passes = level.get("passes") or [level]
        for p in passes:
            for r in p.get("rows", []):
                if not r.get("meets60"):
                    continue
                key = (r["stored"], float(r["rel60"]))
                if best is None or key < best[0]:
                    best = (key, r)
    return best[1] if best else None


def failing_level(levels):
    """The last full sweep at a count, when no pair met there."""
    for level in reversed(levels):
        if level.get("passes") is not None and not level["meets60"]:
            return level
    return None


def count_cell(g, cs, cell, width, target, verbose=True):
    """Walk from the record's count until the 60-digit count is settled.

    If the record's own pair meets at 60 digits the walk descends: every pair
    at each count below is swept and the walk stops where none meets. That
    stopping sweep is the below-count check rational_count applies when it
    certifies its own count as minimal, so the count found here rests on the
    same predicate. If the record's pair misses, every pair at the record's
    count is swept; if none meets there either, the unrounded count is above
    the record's and counts upward are searched until one meets.
    """
    a = g.X0
    b = a + width
    ts, xs = cs.make_nodes(g, a, b, cs.N_NODES)
    started = time.perf_counter()
    ref = cs.reference(g, xs)
    ref_seconds = round(time.perf_counter() - started, 1)
    rec = cell["record"]
    top = cell["start_stored"]
    rec_pair = tuple(cell["record_pair"])
    levels = []
    final = None
    direction = None
    minimality = None

    # A meeting fit settles a count; only a miss needs the full pair sweep.
    confirm = level_pass(g, cs, a, b, ts, xs, ref, top, [rec_pair], target,
                         "remez", verbose)
    confirm["pass"] = "confirm-remez"
    levels.append(confirm)
    if confirm["meets60"]:
        final, direction = top, "equal"

    if final is None:
        recorded_solver = cell.get("record_solver") or "lawson"
        if recorded_solver != "remez":
            fallback = level_pass(g, cs, a, b, ts, xs, ref, top, [rec_pair],
                                  target, recorded_solver, verbose)
            fallback["pass"] = "confirm-recorded-solver"
            levels.append(fallback)
            if fallback["meets60"]:
                final, direction = top, "equal"

    if final is None:
        full = level_below(g, cs, a, b, ts, xs, ref, top, target, verbose)
        full["label"] = "record-count-full-sweep"
        levels.append(full)
        if full["meets60"]:
            final, direction = top, "full-sweep-at-record"

    if final is None:
        for stored in range(top + 1, cs.RAT_MAX_STORED + 1):
            up = level_below(g, cs, a, b, ts, xs, ref, stored, target, verbose)
            up["label"] = "upward-full-sweep"
            levels.append(up)
            if up["meets60"]:
                final, direction = stored, "up"
                break
        if final is not None:
            minimality = (
                "not searched below the record count: the record's own pair and "
                "every pair at its count miss at 60 digits, so the unrounded "
                "count is above the record's")

    if direction in ("equal", "full-sweep-at-record"):
        minimality = (
            "every pair swept at each count below until one failed; that "
            "failing sweep is the below-count check rational_count applies")
        stored = top
        while stored > cs.RAT_MIN_STORED:
            stored -= 1
            level = level_below(g, cs, a, b, ts, xs, ref, stored, target,
                                verbose)
            level["label"] = "below-full-sweep"
            levels.append(level)
            if level["meets60"]:
                final, direction = stored, "down"
            else:
                break

    best = accepted_row(levels)
    fail = failing_level(levels)
    return {
        "id": cell["id"],
        "width_name": cell["width_name"],
        "width": g.mp.nstr(width, 17),
        "target": cs.target_tag(target),
        "a": g.mp.nstr(a, 17),
        "b": g.mp.nstr(b, 17),
        "record": rec,
        "start_stored": top,
        "count60": final,
        "status": (None if final is not None
                   else f"not reached by stored {cs.RAT_MAX_STORED}"),
        "direction": direction,
        "minimality": minimality,
        "delta_vs_start": None if final is None else final - top,
        "delta_vs_record": None if final is None else final - rec["stored"],
        "improved": None if final is None else final < top,
        "count60_m": None if best is None else best["m"],
        "count60_k": None if best is None else best["k"],
        "count60_solver": None if best is None else best["solver"],
        "count60_pair": (None if best is None
                         else f"[{best['m']}/{best['k']}] {best['solver']}"),
        "count60_rel60": None if best is None else best["rel60"],
        "count60_abs60": None if best is None else best["abs60"],
        "count60_rel64": None if best is None else best["rel64"],
        "count60_meets64": None if best is None else best["meets64"],
        "failing_level": None if fail is None else {
            "stored": fail["stored"],
            "routes": [{"pass": p["pass"], "best_rel60": p["best_rel60"],
                        "solved": p["solved"]} for p in fail["passes"]],
        },
        "levels": levels,
        "reference_seconds": ref_seconds,
        "seconds": round(time.perf_counter() - started, 1),
    }


# ---------------------------------------------------------------------------
# cells


def load_record_cells():
    sweep = json.loads(SWEEP.read_text())
    recert = json.loads(RECERT.read_text())
    recert_cells = {c["id"]: c for c in recert["cells"]}
    rows = {}
    for r in sweep["rows"]:
        rows[(r["width_name"], r["target"])] = r
    cells = []
    # The decisive subset: the endpoints of the three increments and the one
    # invalid cell the 2^5 fit excludes.
    # Cheapest first: the two 1e-8 cells exercise the confirm, a full pair
    # sweep and a Lawson negative in minutes, before the deep cells spend
    # hours on the same code path.
    wanted = [
        ("paper", "rel:1e-08"),
        ("ship", "rel:1e-08"), ("2^5", "rel:1e-08"),
        ("paper", "rel:1e-14"),
        ("2^5", "rel:1e-14"),
        ("ship", "rel:1e-15"), ("2^5", "rel:1e-15"),
    ]
    ident = 0
    for name, tag in wanted:
        r = rows[(name, tag)]
        ident += 1
        cell = {"id": ident, "width_name": name, "target": tag,
                "record": {"stored": r["rat"]["stored"],
                           "m": r["rat"].get("m"), "k": r["rat"].get("k"),
                           "solver": r["rat"].get("solver"),
                           "delivered_rel": r["rat"].get("delivered_rel"),
                           "below_meets": r["rat"].get("below_meets"),
                           "source": "results-ladder/sweep_all.json"},
                "record_pair": [r["rat"]["m"], r["rat"]["k"]],
                "record_solver": r["rat"].get("solver"),
                "start_stored": r["rat"]["stored"]}
        # The 2^5 1e-14 row is non-minimal in the record; the re-certification
        # found the true binary64 minimum, which is the fair starting point.
        if name == "2^5" and tag == "rel:1e-14":
            c3 = recert_cells[3]
            cell["record_pair"] = [c3["certifying_pair"]["m"],
                                   c3["certifying_pair"]["k"]]
            cell["record_solver"] = c3["certifying_pair"]["solver"]
            cell["start_stored"] = c3["certified_stored"]
            cell["record"]["re_certified_stored"] = c3["certified_stored"]
            cell["record"]["re_certified_pair"] = cell["record_pair"]
            cell["record"]["re_certified_rel"] = c3["certifying_pair"]["delivered_rel"]
            cell["record"]["source"] = (
                "results-ladder/sweep_all.json; results-ladder/recertification.json cell 3")
        cells.append(cell)
    return cells


def all_cells():
    sweep = json.loads(SWEEP.read_text())
    cells = []
    for i, r in enumerate(sweep["rows"], 1):
        st = r["rat"]["stored"]
        if st is None:
            continue
        cell = {"id": i, "width_name": r["width_name"],
                "target": r["target"],
                "record": {"stored": st, "m": r["rat"].get("m"),
                           "k": r["rat"].get("k"),
                           "solver": r["rat"].get("solver"),
                           "delivered_rel": r["rat"].get("delivered_rel"),
                           "below_meets": r["rat"].get("below_meets"),
                           "source": "results-ladder/sweep_all.json"},
                "record_pair": [r["rat"]["m"], r["rat"]["k"]],
                "record_solver": r["rat"].get("solver"),
                "start_stored": st}
        cells.append(cell)
    return cells


def parse_cell_list(spec):
    out = []
    for part in spec.split(","):
        name, tag = part.split(":", 1)
        out.append((name.strip(), "rel:" + tag.strip()))
    return out


def is_done(cell):
    """A cell that reached an answer, as opposed to one that raised."""
    return ("error" not in cell
            and (cell.get("count60") is not None
                 or cell.get("status") is not None))


def load_previous(out_path, settings):
    """Cells from an earlier invocation under the same counting settings."""
    if not out_path.exists():
        return None
    try:
        data = json.loads(out_path.read_text())
    except Exception as exc:
        print(f"warning: cannot read {out_path} ({exc}); starting fresh",
              flush=True)
        return None
    if data.get("meta", {}).get("schema") != SCHEMA:
        print(f"warning: {out_path.name} has schema "
              f"{data.get('meta', {}).get('schema')!r}, expected {SCHEMA}; "
              f"starting fresh", flush=True)
        return None
    have = data.get("meta", {}).get("settings", {})
    bad = [k for k, v in settings.items() if have.get(k) != v]
    if bad:
        print(f"warning: {out_path.name} was written with other settings "
              f"({', '.join(bad)}); starting fresh", flush=True)
        return None
    return data


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cells", default="decisive",
                    help="decisive | all | name:target,...")
    ap.add_argument("--boys", default=str(DEFAULT_BOYS))
    ap.add_argument("--out", default=str(OUT_DEFAULT))
    ap.add_argument("--max-seconds", type=float, default=0.0,
                    help="stop before starting a new cell once spent")
    ap.add_argument("--workers", type=int, default=0,
                    help="solver pool size; 0 uses one per CPU")
    args = ap.parse_args(argv)

    if str(PAPER) not in sys.path:
        sys.path.insert(0, str(PAPER))
    cs = load_module("count_sweep", PAPER / "count_sweep.py")
    g = cs.load_generator(Path(args.boys))
    if args.workers:
        g.JOBS = args.workers

    if args.cells == "decisive":
        cells = load_record_cells()
    elif args.cells == "all":
        cells = all_cells()
    else:
        cells = []
        wanted = parse_cell_list(args.cells)
        for i, (name, tag) in enumerate(wanted, 1):
            src = {c["width_name"] + ":" + c["target"]: c
                   for c in (load_record_cells() + all_cells())}
            c = src[name + ":" + tag]
            c = dict(c)
            c["id"] = i
            cells.append(c)

    widths = {name: w for name, w in cs.width_set(g, "all")}
    targets = {cs.target_tag(t): t for t in cs.parse_targets(
        "rel:1e-08,rel:1e-09,rel:1e-10,rel:1e-11,rel:1e-12,rel:1e-13,"
        "rel:1e-14,rel:1e-15")}

    out_path = Path(args.out)
    settings = {
        "nodes": cs.N_NODES + 1,
        "remez_iters": cs.RAT_REMEZ_ITERS,
        "refine_outer": cs.RAT_REFINE_OUTER,
        "rat_k_max": K_MAX,
        "reference_dps": int(g.RAT_DPS),
        "evaluation_dps": EVAL_DPS,
    }
    prev = load_previous(out_path, settings)
    done = {}
    resume_block = None
    if prev is not None:
        for c in prev.get("cells", []):
            if is_done(c):
                done[(c["width_name"], c["target"])] = c
        resume_block = {
            "previous_generated_utc": prev.get("meta", {}).get("generated_utc"),
            "previous_script_sha256": prev.get("meta", {}).get("script_sha256"),
            "inherited": [{"id": c.get("id"),
                           "width_name": c.get("width_name"),
                           "target": c.get("target"),
                           "count60": c.get("count60")} for c in done.values()],
            "note": ("a reused cell was computed by the earlier script "
                     "revision; the changes since are the completion print, "
                     "this resume, and the upward direction of the search, none "
                     "of which touches the arithmetic of a cell that completed"),
        }

    record = {
        "meta": {
            "schema": SCHEMA,
            "generated_utc": now(),
            "script": str(Path(__file__).resolve()),
            "script_sha256": sha256(Path(__file__).resolve()),
            "harness": str(PAPER / "count_sweep.py"),
            "harness_sha256": sha256(PAPER / "count_sweep.py"),
            "boys_root": str(args.boys),
            "generator_sha256": sha256(Path(args.boys) / "tools/gen_boys_coefficients.py"),
            "header_sha256": sha256(Path(args.boys) / "include/boys/boys_coefficients.hpp"),
            "cells_source_sha256": sha256(SWEEP),
            "settings": {
                **settings,
                "grid": "uniform in t = 2*(x-a)/(b-a) - 1, endpoints included",
                "acceptance": ("the fit's value at the solver's own precision "
                               "(60 digits) on the certifying grid, relative to "
                               "F_0(x) at each node"),
                "search": ("confirm the record's count with its own pair on the "
                           "minimax route; if it misses, sweep every pair at "
                           "that count and then search upward; if it meets, "
                           "solve every pair at each count below on the full "
                           "grid until none meets, and the Lawson route is run "
                           "at a count only when the Remez route found no "
                           "meeting fit and the negative answer has to be "
                           "certified"),
                "workers": g.worker_count(),
            },
            "resume": resume_block,
        },
        "cells": [],
    }

    started_all = time.perf_counter()
    for cell in cells:
        key = (cell["width_name"], cell["target"])
        if key in done:
            record["cells"].append(done[key])
            print(f"[{cell['id']}] {cell['width_name']} {cell['target']}: "
                  f"reused from {out_path.name} "
                  f"(count60={done[key].get('count60')})", flush=True)
            continue
        if args.max_seconds and time.perf_counter() - started_all > args.max_seconds:
            print("budget spent; stopping before the next cell", flush=True)
            break
        width = widths[cell["width_name"]]
        target = targets[cell["target"]]
        print(f"[{cell['id']}] {cell['width_name']} {cell['target']} "
              f"(record {cell['start_stored']})", flush=True)
        try:
            res = count_cell(g, cs, cell, width, target)
        except Exception as exc:
            res = {"id": cell["id"], "width_name": cell["width_name"],
                   "target": cell["target"], "error": repr(exc),
                   "traceback": traceback.format_exc(),
                   "seconds": round(time.perf_counter() - started_all, 1)}
            record["cells"].append(res)
            out_path.write_text(json.dumps(record, indent=1))
            print(f"    ERROR: {exc!r} (recorded; carrying on)", flush=True)
            continue
        record["cells"].append(res)
        out_path.write_text(json.dumps(record, indent=1))
        print(f"    count60={res['count60']} ({res['direction']}; record "
              f"{res['record']['stored']}, delta {res['delta_vs_record']}) "
              f"pair {res['count60_pair']} rel60 {res['count60_rel60']} "
              f"[{res['seconds']} s]", flush=True)

    seen = {(c.get("width_name"), c.get("target")) for c in record["cells"]}
    for key, c in done.items():
        if key not in seen:
            record["cells"].append(c)

    record["runtime_seconds"] = round(time.perf_counter() - started_all, 1)
    out_path.write_text(json.dumps(record, indent=1))
    print(f"wrote {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
