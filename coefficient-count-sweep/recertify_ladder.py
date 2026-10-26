#!/usr/bin/env python3
"""Re-certify six rational ladder cells on the doubled stability grid.

The dense ladder records a rational stored count for each (width, target).
Six of those cells carry the negative result, and the solver's minimality is
least certain there: the accepted pairs sit at a high numerator and a low
denominator, two rows were flagged below-count, and a count is only the
smallest these solvers reach. This script re-measures the six cells at the
settings the stability pass uses:

    4,001 certifying nodes, uniform in the mapped variable t;
    60-digit references;
    Remez 120 iterations, Lawson 600 outer iterations;
    k <= 12.

It changes no fit and refines no count of its own. The grid, the reference,
the target parser, the solver batch and the acceptance rule all come from
count_sweep.py, which in turn imports the Boys generator's own solvers:

    count_sweep.make_nodes / reference   the certifying grid and references
    count_sweep.refine_pairs             the batched full-grid solve
    count_sweep.meets                    the ladder's acceptance rule
    count_sweep.parse_targets            rel:V / abs:V

refine_pairs drops a fit that misses its target, and this record needs the
delivered error of every pair it tries. The batch is therefore run with the
acceptance threshold set to infinity, which changes no solver and no fit, and
the ladder's real target is applied to the returned fits by count_sweep.meets.
The solve itself is unchanged: every pair is a full-grid Remez and Lawson
solve, in one worker pool, exactly as refine_pairs runs it.

The count is the smallest stored count at which some pair meets on the full
grid, and the stored count below it is settled by solving every pair there,
not a coarse scan. The walk starts one stored count below the record's own
count: if a pair meets, the count walks down until a level where none meets;
if none meets and the recorded pair fails, the walk goes up until one does.

The record is results-ladder/recertification.json. Per-cell fragments and
logs live in _recert/, so an interrupted run resumes from the completed cells
and sweep_all.json is never touched.

Usage:
    python recertify_ladder.py --all
    python recertify_ladder.py --cell 3
    python recertify_ladder.py --merge
"""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import count_sweep as cs  # noqa: E402  the harness this record reuses

SCHEMA = "ladder-recertification-1"
RECORD_REL = Path("results-ladder/recertification.json")
FRAG_DIR_REL = Path("_recert")

# The doubled settings. `nodes` is the count_sweep --nodes argument, so the
# grid has nodes + 1 points. The stability pass's plan line is
# --scan-outer 200 --refine-outer 600 --remez-iters 120 --nodes 4000; the
# coarse scan plays no part in a walk that solves every pair on the full grid.
SETTINGS = {
    "nodes": 4000,
    "grid_points": 4001,
    "remez_iters": 120,
    "refine_outer": 600,
    "rat_k_max": 12,
    "reference_dps": None,          # filled from the generator
    "reference_tail_floor": None,   # filled from the generator
    "reference_terms": None,        # filled from the generator
    "grid": "uniform in t = 2*(x-a)/(b-a) - 1, endpoints included",
    "acceptance": "worse of the fused and separate binary64 evaluation routes, "
                  "relative to F_0(x) at each node",
}

# The six cells. `record_stored` is (m + 1) + k from results-ladder/sweep_all.json
# and is the paper's counting rule. `brief_count` is the figure the work order
# prints; for cells 1-5 it is the numerator degree m of the recorded pair, not
# the stored count, and the record's stored count is the one re-certified here.
CELLS = [
    {
        "id": 1,
        "width": "ship",
        "target": "rel:1e-14",
        "record_stored": 16,
        "brief_count": 12,
        "record_pair": [12, 3],
        "record_below_meets": False,
        "record_solver": "remez",
        "record_delivered_rel": "3.1378578105805594e-15",
        "why": "interior ladder point at the ship width",
    },
    {
        "id": 2,
        "width": "ship",
        "target": "rel:1e-15",
        "record_stored": 19,
        "brief_count": 15,
        "record_pair": [15, 3],
        "record_below_meets": False,
        "record_solver": "lawson",
        "record_delivered_rel": "6.9925232633206473e-16",
        "why": "last point of the ship-width rational increment (+6.29)",
    },
    {
        "id": 3,
        "width": "2^5",
        "target": "rel:1e-14",
        "record_stored": 21,
        "brief_count": 17,
        "record_pair": [17, 3],
        "record_below_meets": True,
        "record_solver": "remez",
        "record_delivered_rel": "5.3492608334410015e-15",
        "why": "flagged below_meets; excluded from the 2^5 rational fit",
    },
    {
        "id": 4,
        "width": "2^5",
        "target": "rel:1e-15",
        "record_stored": 24,
        "brief_count": 21,
        "record_pair": [21, 2],
        "record_below_meets": False,
        "record_solver": "remez",
        "record_delivered_rel": "6.9370546860509225e-16",
        "why": "last point of the 2^5 rational increment (+9.59)",
    },
    {
        "id": 5,
        "width": "paper",
        "target": "rel:1e-15",
        "record_stored": 18,
        "brief_count": 14,
        "record_pair": [14, 3],
        "record_below_meets": True,
        "record_solver": "lawson",
        "record_delivered_rel": "8.5936591774929992e-16",
        "why": "flagged below_meets; excluded from the paper-width fit",
    },
    {
        "id": 6,
        "width": "2^5",
        "target": "rel:1e-13",
        "record_stored": 16,
        "brief_count": 16,
        "record_pair": [7, 8],
        "record_below_meets": False,
        "record_solver": "lawson",
        "record_delivered_rel": "9.3570469664909625e-14",
        "why": "the stability pass moved this cell to 17 [8/8]; the k <= 12 "
               "walk decides 16 or 17",
    },
]

# refine_pairs' only use of its target is the acceptance filter. An infinite
# relative threshold keeps every finite fit, so the record carries the
# delivered error of every pair; the real target is applied below.
LOOSE_TARGET = {"kind": "rel", "value": float("inf"), "spec": "rel:inf"}


def now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


class Log:
    """Print a progress line and append it to a cell's own log file."""

    def __init__(self, path: Path | None):
        self.path = path
        self.fh = None
        if path is not None:
            path.parent.mkdir(parents=True, exist_ok=True)
            self.fh = path.open("a", encoding="utf-8", buffering=1)

    def __call__(self, message: str):
        line = f"[{now()}] {message}"
        print(line, flush=True)
        if self.fh is not None:
            self.fh.write(line + "\n")

    def close(self):
        if self.fh is not None:
            self.fh.close()
            self.fh = None


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def command_line(argv) -> str:
    return subprocess.list2cmdline([sys.executable, str(Path(__file__).resolve()),
                                    *argv])


def solver_fits(fits, m, k):
    return [f for f in fits if f["m"] == m and f["k"] == k]


def solve_level(g, a, b, ts, xs, ref, target, stored, k_max,
                refine_outer, remez_iters, log: Log) -> dict:
    """Solve every pair m + 1 + k = stored with k <= k_max on the full grid."""
    pairs = [(stored - 1 - k, k)
             for k in range(1, min(k_max, stored - 2) + 1)]
    started = time.perf_counter()
    log(f"cell level stored={stored}: {len(pairs)} pairs "
        f"{['[%d/%d]' % p for p in pairs]} on {len(ts)} nodes, "
        f"Remez {remez_iters}, Lawson {refine_outer}, workers "
        f"{g.worker_count()}")
    fits = cs.refine_pairs(g, pairs, a, b, ts, xs, ref, LOOSE_TARGET,
                           refine_outer=refine_outer,
                           remez_iters=remez_iters)
    seconds = time.perf_counter() - started
    rows = []
    for m, k in pairs:
        mine = solver_fits(fits, m, k)
        solvers = []
        for fit in mine:
            solvers.append({
                "solver": fit["name"],
                "delivered_rel": g.mp.nstr(fit["rel"], 17),
                "delivered_abs": g.mp.nstr(fit["abs"], 17),
                "at_x_rel": g.mp.nstr(fit["at_rel"], 17),
                "at_x_abs": g.mp.nstr(fit["at_abs"], 17),
                "meets": cs.meets(target, fit["rel"], fit["abs"]),
            })
        # The delivered figure of a pair is its best fit under the target's
        # own criterion, the reading the sweep records for a count.
        if solvers:
            key = "delivered_rel" if target["kind"] == "rel" else "delivered_abs"
            best = min(solvers, key=lambda s: g.mpf(s[key]))
            best = dict(best)
        else:
            best = None
        met = any(s["meets"] for s in solvers)
        rows.append({
            "m": m,
            "k": k,
            "stored": (m + 1) + k,
            "solvers": solvers,
            "best": best,
            "meets": met,
            "no_finite_fit": not solvers,
        })
    meeting = [r for r in rows if r["meets"]]
    meeting.sort(key=lambda r: g.mpf(r["best"]["delivered_rel"]))
    best_all = None
    ranked = [r for r in rows if r["best"] is not None]
    if ranked:
        ranked.sort(key=lambda r: g.mpf(r["best"]["delivered_rel"]))
        best_all = ranked[0]
    log(f"cell level stored={stored}: done in {seconds:.1f}s; "
        f"meets={bool(meeting)}; "
        + (f"best meeting pair [{meeting[0]['m']}/{meeting[0]['k']}] "
           f"{meeting[0]['best']['solver']} rel {meeting[0]['best']['delivered_rel']}"
           if meeting else
           f"best pair [{best_all['m']}/{best_all['k']}] {best_all['best']['solver']} "
           f"rel {best_all['best']['delivered_rel']}" if best_all else
           "no finite fit at any pair"))
    return {
        "stored": stored,
        "pairs_tried": [f"[{m}/{k}]" for m, k in pairs],
        "rows": rows,
        "any_meets": bool(meeting),
        "meeting_pairs": [f"[{r['m']}/{r['k']}]" for r in meeting],
        "best_meeting": meeting[0] if meeting else None,
        "best_pair": best_all,
        "seconds": round(seconds, 1),
    }


def cell_result(g, cell, args, log: Log) -> dict:
    started = time.perf_counter()
    a = g.X0
    _, width = cs.width_set(g, cell["width"])[0]
    b = a + width
    target = cs.parse_targets(cell["target"])[0]
    log(f"cell {cell['id']}: width {cell['width']} = {g.mp.nstr(width, 17)}, "
        f"target {cell['target']}, record stored {cell['record_stored']} "
        f"[{cell['record_pair'][0]}/{cell['record_pair'][1]}], "
        f"brief count {cell['brief_count']}")
    grid_started = time.perf_counter()
    ts, xs = cs.make_nodes(g, a, b, args.nodes)
    ref = cs.reference(g, xs)
    grid_seconds = time.perf_counter() - grid_started
    log(f"cell {cell['id']}: grid {len(ts)} nodes and 60-digit reference "
        f"built in {grid_seconds:.1f}s")

    levels = {}

    def level(stored):
        if stored not in levels:
            levels[stored] = solve_level(
                g, a, b, ts, xs, ref, target, stored, args.rat_k_max,
                args.refine_outer, args.remez_iters, log)
        return levels[stored]

    stored = cell["record_stored"]
    walk = []
    below = level(stored - 1)
    walk.append({"stored": stored - 1, "role": "one below the record count",
                 "any_meets": below["any_meets"]})
    if below["any_meets"]:
        # The record's count is not minimal. Walk down while the count below
        # the candidate still has a meeting pair.
        candidate = stored - 1
        while candidate - 1 >= cs.RAT_MIN_STORED:
            down = level(candidate - 1)
            walk.append({"stored": candidate - 1, "role": "walk down",
                         "any_meets": down["any_meets"]})
            if not down["any_meets"]:
                break
            candidate -= 1
        certified = candidate
        old_status = "non_minimal (a pair meets one stored count below)"
    else:
        at = level(stored)
        walk.append({"stored": stored, "role": "the record count",
                     "any_meets": at["any_meets"]})
        if at["any_meets"]:
            certified = stored
            old_status = "minimal"
        else:
            # No pair at the record's count meets. The count moves up; the
            # level below is already settled.
            old_status = "no_pair_meets (no pair at the record count meets)"
            certified = None
            up = stored
            while up + 1 <= cs.RAT_MAX_STORED:
                up += 1
                found = level(up)
                walk.append({"stored": up, "role": "walk up",
                             "any_meets": found["any_meets"]})
                if found["any_meets"]:
                    certified = up
                    break
        # A record count that no pair reaches is reported under its own
        # heading even though a higher count is then certified.
    at_level = levels.get(certified)
    certifying = None
    if at_level is not None and at_level["best_meeting"] is not None:
        bm = at_level["best_meeting"]
        certifying = {
            "m": bm["m"], "k": bm["k"], "stored": bm["stored"],
            "solver": bm["best"]["solver"],
            "delivered_rel": bm["best"]["delivered_rel"],
            "delivered_abs": bm["best"]["delivered_abs"],
            "at_x_rel": bm["best"]["at_x_rel"],
            "meets": True,
        }
    if certified is None:
        verdict = "no pair meets"
    elif certified == stored:
        verdict = "holds"
    else:
        verdict = "moves"
    below_certified = levels.get(certified - 1) if certified is not None else None
    result = {
        "id": cell["id"],
        "why": cell["why"],
        "width_name": cell["width"],
        "width": g.mp.nstr(width, 17),
        "width_f": float(width),
        "a": g.mp.nstr(a, 17),
        "b": g.mp.nstr(b, 17),
        "target": cell["target"],
        "target_kind": target["kind"],
        "target_value": target["value"],
        "record_stored": stored,
        "record_pair": cell["record_pair"],
        "record_solver": cell["record_solver"],
        "record_delivered_rel": cell["record_delivered_rel"],
        "record_below_meets": cell["record_below_meets"],
        "brief_count": cell["brief_count"],
        "brief_count_note": (
            "the work order's count for this cell; for cells 1-5 it is the "
            "numerator degree m of the recorded pair, not (m+1)+k"
            if cell["brief_count"] != stored else
            "the work order's count equals the record's stored count"),
        "walk": walk,
        "levels": [levels[k] for k in sorted(levels)],
        "old_count_status": old_status,
        "old_count_certifies": certified == stored,
        "no_pair_meets_at_old_count": not levels[stored]["any_meets"]
        if stored in levels else None,
        "certified_stored": certified,
        "certifying_pair": certifying,
        "certified_level_meets": at_level["any_meets"] if at_level else None,
        "level_below_certified_meets": (below_certified["any_meets"]
                                        if below_certified else None),
        "verdict": verdict,
        "grid_seconds": round(grid_seconds, 1),
        "seconds": round(time.perf_counter() - started, 1),
    }
    log(f"cell {cell['id']}: verdict {verdict}; record {stored} -> "
        f"certified {certified}"
        + (f" by [{certifying['m']}/{certifying['k']}] "
           f"{certifying['solver']} rel {certifying['delivered_rel']}"
           if certifying else "")
        + f"; {result['seconds']}s")
    return result


def meta(args) -> dict:
    boys = Path(args.boys_root)
    gen = boys / cs.GENERATOR_REL
    header = boys / cs.HEADER_REL
    g = cs.load_generator(boys)
    settings = dict(SETTINGS)
    settings.update({
        "reference_dps": int(g.RAT_DPS),
        "reference_tail_floor": g.mp.nstr(g.RAT_TAIL_FLOOR, 8),
        "reference_terms": int(g.RAT_TERMS),
        "remez_iters": args.remez_iters,
        "refine_outer": args.refine_outer,
        "rat_k_max": args.rat_k_max,
        "nodes": args.nodes,
        "grid_points": args.nodes + 1,
        "workers": g.worker_count(),
    })
    return {
        "schema": SCHEMA,
        "generated_utc": now(),
        "script": str(Path(__file__).resolve()),
        "script_sha256": sha256(Path(__file__).resolve()),
        "reused_harness": {
            "path": str(HERE / "count_sweep.py"),
            "sha256": sha256(HERE / "count_sweep.py"),
            "functions": [
                "load_generator", "make_nodes", "reference", "refine_pairs",
                "meets", "parse_targets", "width_set", "RAT_MIN_STORED",
                "RAT_MAX_STORED",
            ],
        },
        "boys_root": str(boys),
        "git": cs.git_rev(boys),
        "papers_git": cs.git_rev(HERE.parents[1]),
        "header": {"path": str(header), "sha256": sha256(header)},
        "generator": {"path": str(gen), "sha256": sha256(gen)},
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "settings": settings,
        "precision_scope": cs.PRECISION_SCOPE,
        "method": {
            "count_unit": "(m + 1) numerator coefficients + k denominator "
                          "coefficients, q_0 held at 1",
            "certified": "the smallest stored count at which some (m,k) pair "
                         "with k <= 12 meets the target on the full 4,001-node "
                         "grid; the count below it is settled by solving every "
                         "pair there, not a coarse scan",
            "batch": "count_sweep.refine_pairs with its acceptance threshold "
                     "set to infinity, so every solved pair's delivered error "
                     "is returned; the ladder target is applied to the fits by "
                     "count_sweep.meets",
            "walk": "start at the record's stored count minus one; if a pair "
                    "meets there, walk down while the count below still has a "
                    "meeting pair; otherwise check the record count and, if no "
                    "pair meets, walk up to the first count that does",
            "target_semantics": "relative to F_0(x) at each node, the ladder's "
                                "own reading",
        },
        "cells_source": str(HERE / "results-ladder" / "sweep_all.json"),
        "cells_source_sha256": sha256(HERE / "results-ladder" / "sweep_all.json"),
        "provenance_note": (
            "boys git records the worktree HEAD when a fragment was written; "
            "the generator and header sha256 digests are the content identity "
            "and are carried per cell in content_identity"
        ),
    }


def run_cell(args, cell, log: Log) -> dict:
    boys = Path(args.boys_root)
    g = cs.load_generator(boys)
    g.JOBS = args.jobs if args.jobs and args.jobs > 0 else None
    result = {
        "command": command_line(["--cell", str(cell["id"])]),
        "started_utc": now(),
        "meta": meta(args),
        "cell": cell_result(g, cell, args, log),
        "finished_utc": now(),
    }
    frag = Path(args.frag_dir) / f"cell_{cell['id']}.json"
    frag.parent.mkdir(parents=True, exist_ok=True)
    frag.write_text(json.dumps(cs.json_ready(result), indent=2) + "\n",
                    encoding="utf-8")
    log(f"cell {cell['id']}: fragment written to {frag}")
    return result


# The cells the paper's three residuals and the ladder table use. Cell 3
# (2^5, 1e-14) and cell 5 (paper, 1e-15) are the two rows the paper excludes
# as non-minimal. Cell 3's re-certification completed before the owner stopped
# the run and is carried in the record as an excluded-row measurement; cell 5
# was never started and is recorded as not re-certified. The ladder-analysis
# scenario applies only the claim cells, so the paper's exclusions stay in
# place and the scenario answers what the claims read after re-certification.
CLAIM_CELL_IDS = (1, 2, 4, 6)


def modified_sweep_record(cells, include_ids=CLAIM_CELL_IDS):
    """The committed ladder record with the claim cells at their certified counts.

    This is a scenario for the reading, never a replacement: sweep_all.json is
    read and copied into a scratch directory. A re-certified count is minimal by
    construction (the count below it was solved in full and no pair meets), so
    its below_meets flag is false.
    """
    path = HERE / "results-ladder" / "sweep_all.json"
    data = json.loads(path.read_text(encoding="utf-8"))
    by_key = {(c["width_name"], c["target"]): c
              for c in cells if c["id"] in include_ids}
    for row in data["rows"]:
        cell = by_key.get((row["width_name"], row["target"]))
        if cell is None or cell["certified_stored"] is None:
            continue
        pair = cell["certifying_pair"]
        row["rat"].update({
            "stored": cell["certified_stored"],
            "m": pair["m"],
            "k": pair["k"],
            "solver": pair["solver"],
            "delivered_rel": pair["delivered_rel"],
            "delivered_abs": pair["delivered_abs"],
            "below_meets": False,
        })
    return data


def ladder_analysis_scenario(cells, work_dir: Path):
    """Run the paper's own ladder_analysis.py over the re-certified claim cells."""
    data = modified_sweep_record(cells)
    work_dir.mkdir(parents=True, exist_ok=True)
    (work_dir / "sweep_all.json").write_text(
        json.dumps(data, indent=2) + "\n", encoding="utf-8")
    out = work_dir / "ladder_analysis.json"
    proc = subprocess.run(
        [sys.executable, str(HERE / "ladder_analysis.py"),
         "--dir", str(work_dir), "--out", str(out)],
        cwd=str(HERE), capture_output=True, text=True)
    if proc.returncode != 0:
        return {"ran": False, "returncode": proc.returncode,
                "stderr": proc.stderr[-4000:]}
    report = json.loads(out.read_text(encoding="utf-8"))
    compact = {}
    for width, families in report.get("widths", {}).items():
        compact[width] = {}
        for fam, entry in families.items():
            compact[width][fam] = {
                key: entry.get(key) for key in (
                    "points_total", "points_above_floor", "counts",
                    "excluded_non_minimal", "slope_measured",
                    "slope_predicted", "max_abs_residual",
                    "max_abs_second_difference", "affine_within_tolerance",
                )
            }
    return {
        "ran": True,
        "cells_applied": list(CLAIM_CELL_IDS),
        "input": "results-ladder/sweep_all.json copied to a scratch directory, "
                 "the claim cells (1, 2, 4, 6) replaced by their certified "
                 "counts and pair; the paper's excluded rows stay at their "
                 "committed values",
        "verdict": report.get("verdict"),
        "complete": report.get("complete"),
        "flips": report.get("flips"),
        "widths": compact,
        "stdout_tail": proc.stdout.strip().splitlines()[-1:] if proc.stdout else [],
    }


# The paper's model-based residuals: measured increment minus the prediction,
# both over the span named. Only a cell that is the last point of a span can
# move its number; the first points (7 at paper and ship, 8 at 2^5) are not
# among the six cells.
RESIDUAL_SPANS = [
    {"name": "+2.22 paper width, 1e-8 to 1e-14",
     "width": "paper", "first_target": "rel:1e-8",
     "last_target": "rel:1e-14", "first_count": 7, "decades": 6,
     "predicted_per_decade": 0.796, "residual": 2.22,
     "endpoint_cell": None},
    {"name": "+6.29 ship width, 1e-8 to 1e-15",
     "width": "ship", "first_target": "rel:1e-8",
     "last_target": "rel:1e-15", "first_count": 7, "decades": 7,
     "predicted_per_decade": 0.816, "residual": 6.29,
     "endpoint_cell": 2},
    {"name": "+9.59 2^5 width, 1e-8 to 1e-15",
     "width": "2^5", "first_target": "rel:1e-8",
     "last_target": "rel:1e-15", "first_count": 8, "decades": 7,
     "predicted_per_decade": 0.916, "residual": 9.59,
     "endpoint_cell": 4},
]


ARM_STATS = HERE / "results-sensitivity" / "arm_stats.py"
CORE_ARM_RECORD = HERE / "results" / "sweep_rel_1e-13.json"


def arm_stats_comparison(work_dir: Path) -> dict:
    """Run the paper's arm_stats.py at the core-run count and at the withdrawn 17.

    The paper's core-arm exponents, endpoint secants, local exponents, two-point
    exponents and dominance slope were computed with the stability pass's
    corrected count 17. The re-certification restores the core-run count 16;
    both readings are run through the paper's own script so the replacement
    numbers come from the record, not from a retyping.
    """
    base = json.loads(CORE_ARM_RECORD.read_text(encoding="utf-8"))
    work_dir.mkdir(parents=True, exist_ok=True)
    out = {"script": str(ARM_STATS), "input": str(CORE_ARM_RECORD),
           "readings": {}}
    for tag, stored in (("core_run_16", None), ("withdrawn_17", 17)):
        data = json.loads(json.dumps(base))
        if stored is not None:
            for row in data["rows"]:
                if (abs(row["width_f"] - 32.0) < 1e-12
                        and row["target"] == "rel:1e-13"):
                    row["rat"]["stored"] = stored
        path = work_dir / f"arm_{tag}.json"
        path.write_text(json.dumps(data), encoding="utf-8")
        proc = subprocess.run(
            [sys.executable, str(ARM_STATS), str(path), tag],
            cwd=str(HERE), capture_output=True, text=True)
        lines = [line for line in proc.stdout.splitlines()
                 if line.startswith(("cheb:", "rat:", "D slope"))]
        out["readings"][tag] = {
            "returncode": proc.returncode,
            "lines": lines,
        }
    out["note"] = (
        "the 2^5 rational count at rel:1e-13 is 16 in the core record and in "
        "the re-certification; 17 is shown only to identify the numbers the "
        "withdrawn stability correction produced")
    return out


def consequences(cells):
    """What each of the paper's three residuals reads after re-certification."""
    by_key = {(c["width_name"], c["target"]): c for c in cells}
    out = {"spans": [], "excluded_rows": []}
    for span in RESIDUAL_SPANS:
        entry = dict(span)
        end = by_key.get((span["width"], span["last_target"]))
        if span["endpoint_cell"] is None:
            entry.update({
                "changed": False,
                "new_residual": span["residual"],
                "reading": (
                    "unchanged: the last point of this span is not one of the "
                    "six re-certified cells, and neither is its first point"
                ),
            })
        elif end is None or end["certified_stored"] is None:
            entry.update({"changed": None,
                          "reading": "no certified endpoint in this record"})
        else:
            delta = end["certified_stored"] - end["record_stored"]
            new_residual = span["residual"] + delta
            entry.update({
                "old_endpoint_stored": end["record_stored"],
                "new_endpoint_stored": end["certified_stored"],
                "delta_coefficients": delta,
                "changed": delta != 0,
                "new_residual": round(new_residual, 2),
                "reading": ("unchanged" if delta == 0 else
                            f"moves by {delta:+d} to {new_residual:+.2f}"),
            })
        out["spans"].append(entry)
    # The paper width's 1e-15 row is excluded from +2.22; if the walk now
    # certifies a minimal count there, a seven-decade variant is possible. It
    # is reported beside the paper's six-decade +2.22, not in place of it.
    paper15 = by_key.get(("paper", "rel:1e-15"))
    if paper15 is not None and paper15["certified_stored"] is not None:
        n = paper15["certified_stored"]
        seven = round(n - 7 - 7 * 0.796, 2)
        out["excluded_rows"].append({
            "cell": 5, "width": "paper", "target": "rel:1e-15",
            "record_stored": paper15["record_stored"],
            "certified_stored": n,
            "effect_on_written_+2.22": (
                "none: the written +2.22 is the six-decade span 1e-8 to 1e-14; "
                "cell 5 is outside it and was excluded as non-minimal"),
            "seven_decade_variant_if_the_row_is_now_valid": {
                "measured_increment": n - 7,
                "predicted": round(7 * 0.796, 2),
                "residual": seven,
            },
        })
    else:
        out["excluded_rows"].append({
            "cell": 5, "width": "paper", "target": "rel:1e-15",
            "record_stored": 18,
            "certified_stored": None,
            "status": "not re-certified",
            "effect_on_written_+2.22": (
                "none: the written +2.22 is the six-decade span 1e-8 to 1e-14, "
                "whose endpoints are not cells of this record; cell 5 is "
                "outside that span and is excluded by the paper as non-minimal"),
            "seven_decade_variant_if_the_row_is_now_valid": (
                "not evaluated: cell 5 was not re-certified"),
        })
    return out


def merge(args) -> dict:
    frag_dir = Path(args.frag_dir)
    cells = []
    commands = []
    provenance = []
    for cell in CELLS:
        frag = frag_dir / f"cell_{cell['id']}.json"
        if not frag.exists():
            print(f"warning: no fragment for cell {cell['id']} ({frag})")
            continue
        data = json.loads(frag.read_text(encoding="utf-8"))
        cells.append(data["cell"])
        commands.append(data["command"])
        frag_meta = data["meta"]
        provenance.append({
            "cell": cell["id"],
            "started_utc": data.get("started_utc"),
            "finished_utc": data.get("finished_utc"),
            "git": frag_meta.get("git"),
            "papers_git": frag_meta.get("papers_git"),
            "generator_sha256": frag_meta.get("generator", {}).get("sha256"),
            "header_sha256": frag_meta.get("header", {}).get("sha256"),
            "settings": frag_meta.get("settings"),
        })
    gen_digests = sorted({p["generator_sha256"] for p in provenance})
    header_digests = sorted({p["header_sha256"] for p in provenance})
    content_identity = {
        "generator_sha256": gen_digests,
        "header_sha256": header_digests,
        "same_generator_for_every_cell": len(gen_digests) <= 1,
        "same_header_for_every_cell": len(header_digests) <= 1,
        "note": (
            "the generator and header digests, not the boys git HEAD, identify "
            "the code that produced the fits; a worktree checkout can move the "
            "HEAD without touching either file, and each cell records the HEAD "
            "it saw at its own start"
        ),
    }
    present = {c["id"] for c in cells}
    scope = {
        "claim_cells": list(CLAIM_CELL_IDS),
        "re_certified_and_present": sorted(present),
        "re_certified_but_excluded_by_the_paper": (
            [3] if 3 in present else []),
        "not_re_certified": [
            {
                "id": cell["id"],
                "width": cell["width"],
                "target": cell["target"],
                "record_stored": cell["record_stored"],
                "reason": (
                    "not re-certified: the paper excludes this row as "
                    "non-minimal and does not report its true minimum; the "
                    "owner stopped the doubling run after cell 3"
                    if cell["id"] == 5 else
                    "no fragment present at merge time"
                ),
            }
            for cell in CELLS if cell["id"] not in present
        ],
        "note": (
            "the three residuals use the ship endpoint (cell 2), the 2^5 "
            "endpoint (cell 4), and the ladder's interior ship and 2^5 cells "
            "(1 and 6); cell 3 is a completed re-certification of an excluded "
            "row and is recorded, not used; cell 5 was not re-certified"
        ),
    }
    record = {
        "meta": meta(args),
        "commands": commands,
        "batch_command": command_line(["--all"]),
        "merge_command": command_line(["--merge"]),
        "scope": scope,
        "provenance": provenance,
        "content_identity": content_identity,
        "cells": cells,
        "consequences": consequences(cells),
        "arm_stats_recomputation": arm_stats_comparison(
            Path(args.frag_dir) / "_analysis"),
        "ladder_analysis_scenario": ladder_analysis_scenario(
            cells, Path(args.frag_dir) / "_analysis"),
        "runtime_seconds": round(sum(c["seconds"] for c in cells), 1),
    }
    record["meta"]["generated_utc"] = now()
    out = Path(args.record_out)
    if not out.is_absolute():
        out = HERE / out
    out.write_text(json.dumps(cs.json_ready(record), indent=2) + "\n",
                   encoding="utf-8")
    print(f"record written to {out} ({len(cells)} cells, "
          f"{record['runtime_seconds']}s of cell runtime)")
    return record


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--boys-root", default=str(cs.DEFAULT_BOYS))
    parser.add_argument("--cell", type=int, default=0,
                        help="run one cell by number (1-6)")
    parser.add_argument("--all", action="store_true",
                        help="run the six cells one at a time, then merge")
    parser.add_argument("--merge", action="store_true",
                        help="assemble results-ladder/recertification.json "
                             "from the per-cell fragments")
    parser.add_argument("--nodes", type=int, default=SETTINGS["nodes"])
    parser.add_argument("--remez-iters", type=int,
                        default=SETTINGS["remez_iters"])
    parser.add_argument("--refine-outer", type=int,
                        default=SETTINGS["refine_outer"])
    parser.add_argument("--rat-k-max", type=int, default=SETTINGS["rat_k_max"])
    parser.add_argument("--jobs", type=int, default=0)
    parser.add_argument("--frag-dir", default=str(FRAG_DIR_REL))
    parser.add_argument("--record-out", default=str(RECORD_REL),
                        help="record path for --merge")
    args = parser.parse_args(argv)

    modes = [args.cell != 0, args.all, args.merge]
    if sum(bool(m) for m in modes) != 1:
        parser.error("name exactly one of --cell N, --all, --merge")

    if args.merge:
        merge(args)
        return 0

    if args.all:
        combined = Log(Path(args.frag_dir) / "all.log")
        for cell in CELLS:
            frag = Path(args.frag_dir) / f"cell_{cell['id']}.json"
            if frag.exists():
                try:
                    data = json.loads(frag.read_text(encoding="utf-8"))
                    settings = data["meta"]["settings"]
                    same = (settings["grid_points"] == args.nodes + 1
                            and settings["remez_iters"] == args.remez_iters
                            and settings["refine_outer"] == args.refine_outer
                            and settings["rat_k_max"] == args.rat_k_max)
                except Exception:
                    same = False
                if same:
                    combined(f"cell {cell['id']}: fragment present, skipped")
                    continue
                combined(f"cell {cell['id']}: fragment has other settings, "
                         f"re-running")
            combined(f"cell {cell['id']}: starting")
            run_cell(args, cell, combined)
        combined("all cells done; merging")
        merge(args)
        combined.close()
        return 0

    cell = next((c for c in CELLS if c["id"] == args.cell), None)
    if cell is None:
        parser.error(f"cell must be 1-6, not {args.cell}")
    log = Log(Path(args.frag_dir) / f"cell_{cell['id']}.log")
    try:
        run_cell(args, cell, log)
    finally:
        log.close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
