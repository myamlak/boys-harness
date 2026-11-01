#!/usr/bin/env python3
"""Measure Paper C's stored-coefficient counts on F_1 at both Group A edges.

The paper measures the width law on F_0(x) at the left edge x0.  A referee
holds that the law is measured at one instance of the class.  This driver
runs the paper's own protocol on a second instance, F_1(x) = (F_0(x) -
exp(-x))/(2x), whose algebraic tail is alpha = 3/2:

    edge x0 : a = X0                      widths = the fifteen Group A widths
    edge a2 : a = 2 X0                    widths = every Group A width doubled

so W/a is matched at the two edges and kappa(W/a) is identical.

What is the paper's and what is new
-----------------------------------
* grid, targets, acceptance, search, evaluators and solvers are
  count_sweep.py's (the paper's harness) and the kernel generator's
  (tools/gen_boys_coefficients.py, C:\\Users\\myaml\\dev\\_mainonly);
* the Chebyshev fit source is patched, module-wide, to f1_reference's
  f1_cheb_fit, which is the generator's own cheb_coeffs_of over the F_1
  reference (the generator carries that entry point for fits whose target is
  not F itself);
* the reference handed to the acceptance test and to the rational solvers is
  f1_reference.f1_reference, the paper's boys_ref plus the exact recurrence;
* nothing else is reimplemented: cheb_count, rational_count, scan_level,
  refine_pairs, rational_candidates, rat_errors, meets, make_nodes and
  write_csv are the paper's functions.

Exactness of the cross-target sharing
-------------------------------------
A cell is one (edge, width) pair measured at all three targets.  The rational
scan at a given stored count is target-independent: the solver job is a pure
function of (m, k, node indices, t grid, reference grid, outer settings).
The driver therefore memoizes the generator's run_jobs at the job boundary,
so the second and third targets at the same width reuse the first target's
solves of the same (m, k, level).  No count depends on the memo: an
identical job returns its identical answer, and acceptance is still applied
per target by count_sweep.rational_count.  Hits and misses are logged.

Rows are written in results-second-edge's row schema (schema 2, meta,
precision_scope, rows), extended with function, edge, width_x0 and parity.

Resumable: a completed target row is written to disk immediately; a re-run
with the same arguments reuses rows whose edge, left edge compared as text,
width and settings match, and re-runs everything else.  A solver exception
writes no row, so the cell is retried by the next invocation.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import sys
import time
from pathlib import Path

PAPER = Path(__file__).resolve().parent
sys.path.insert(0, str(PAPER))

import count_sweep as cs  # noqa: E402
import f1_reference as f1r  # noqa: E402

DEFAULT_BOYS = Path(r"C:\Users\myaml\dev\_mainonly")
DEFAULT_OUT = PAPER / "results-second-instance-f1"
DEFAULT_TARGETS = "rel:1e-13,abs:5e-14,rel:1e-14"
EDGES = {"x0": 1, "a2": 2}

# The measuring order.  The narrowest width first at both edges (the cheapest
# end-to-end proof of the pipeline), then the three wide cells where the
# rational departure lives, then the mid and narrow widths.  The second-edge
# lane lost its wide rational cells to the clock; this order cannot.
ORDER_WIDTHS = ["2^-7", "2^5", "ship", "paper", "2^0",
                "2^-6", "2^-5", "2^-4", "2^-3", "2^-2", "2^-1",
                "2^1", "2^2", "2^3", "2^4"]

IDENTITY_KEYS = ("nodes", "coarse_stride", "scan_outer", "refine_outer",
                 "remez_iters", "cheb_degree_mode", "rat_k_max")


def parity_of(entry):
    """The degree parity of a Chebyshev row, and m/k parity of a rational one."""
    out = {"cheb_degree": None, "cheb_degree_parity": None,
           "rat_m": None, "rat_k": None,
           "rat_m_parity": None, "rat_k_parity": None}
    m = entry.get("cheb", {}).get("m")
    if m is not None:
        out["cheb_degree"] = m
        out["cheb_degree_parity"] = "odd" if m % 2 else "even"
    rm = entry.get("rat", {}).get("m")
    rk = entry.get("rat", {}).get("k")
    if rm is not None:
        out["rat_m"] = rm
        out["rat_m_parity"] = "odd" if rm % 2 else "even"
    if rk is not None:
        out["rat_k"] = rk
        out["rat_k_parity"] = "odd" if rk % 2 else "even"
    return out


def kernel_version(boys: Path) -> str:
    path = boys / "include" / "boys" / "version.hpp"
    if not path.exists():
        return "unknown"
    text = path.read_text(encoding="utf-8", errors="replace")
    found = re.search(r'kVersionString\s*=\s*"([^"]+)"', text)
    return found.group(1) if found else "unknown"


def build_meta(boys: Path, args, g) -> dict:
    """Provenance and settings.  Digests only: the lane does not run git."""
    return {
        "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "schema": cs.ROW_SCHEMA,
        "function": args.function,
        "function_definition": (
            "F_1(x) = (F_0(x) - exp(-x))/(2x) from the n=0->1 Boys "
            "recurrence; F_0 = boys_ref(0, x, RAT_TAIL_FLOOR, RAT_TERMS) at "
            "60 digits; the exact tail is (sqrt(pi)/4) x^(-3/2) g(x), "
            "alpha = 3/2"
            if args.function == "F1" else
            "F_0(x) = boys_ref(0, x, RAT_TAIL_FLOOR, RAT_TERMS) at 60 digits "
            "(the paper's function, as an instrument control)"),
        "boys_root": str(boys),
        "kernel_version": kernel_version(boys),
        "header": {"path": str(boys / cs.HEADER_REL),
                   "sha256": cs.sha256(boys / cs.HEADER_REL)},
        "generator": {"path": str(boys / cs.GENERATOR_REL),
                      "sha256": cs.sha256(boys / cs.GENERATOR_REL)},
        "harness": {"path": str(PAPER / "count_sweep.py"),
                    "sha256": cs.sha256(PAPER / "count_sweep.py")},
        "f1_reference": {"path": str(PAPER / "f1_reference.py"),
                         "sha256": cs.sha256(PAPER / "f1_reference.py")},
        "driver": {"path": str(Path(__file__).resolve()),
                   "sha256": cs.sha256(Path(__file__).resolve())},
        "python": sys.version.split()[0],
        "platform": __import__("platform").platform(),
        "precision_scope": cs.PRECISION_SCOPE,
        "grid": {"nodes": args.nodes + 1, "coarse_stride": args.stride},
        "solver": {
            "scan_outer": args.scan_outer,
            "refine_outer": args.refine_outer,
            "remez_iters": args.remez_iters,
            "cheb_deg": [cs.CHEB_MIN_DEG, cs.CHEB_MAX_DEG],
            "cheb_degree_mode": args.cheb_degrees,
            "cheb_degrees": list(cs.cheb_degrees(args.cheb_degrees)),
            "rat_stored": [cs.RAT_MIN_STORED, cs.RAT_MAX_STORED],
            "rat_k_max": int(args.rat_k_max) if args.rat_k_max else int(g.RAT_K_MAX),
        },
        "left_edges": {"x0": g.mp.nstr(g.X0, 17),
                       "a2": g.mp.nstr(2 * g.X0, 17)},
        "width_scaling": (
            "edge a2 uses every Group A width doubled, so W/a and kappa(W/a) "
            "are matched to the x0 edge; edge x0 uses the Group A widths"),
        "targets": {
            "rel:1e-13": "relative to F_1(x) at each grid node, 1e-13",
            "abs:5e-14": ("plain absolute 5e-14 on F_1(x), the same numerical "
                          "target as the paper's F_0 absolute arm (which is "
                          "5e-14 absolute on F_0); for F_1 the model's "
                          "natural scale is F_1(a), not F_0(a)"),
            "rel:1e-14": "relative to F_1(x) at each grid node, 1e-14",
        },
        "lane": {
            "jobs": int(args.jobs),
            "solver_memo": ("the generator's run_jobs is memoized at the job "
                            "boundary across the three targets of one width; "
                            "an identical solve job returns its identical "
                            "answer, so no count changes"),
        },
    }


def row_identity_ok(row, edge, a_text, function, current):
    if row.get("edge") != edge or row.get("a") != a_text:
        return False
    if row.get("function") != function:
        return False
    settings = row.get("settings", {})
    checks = {
        "nodes": settings.get("nodes"),
        "coarse_stride": settings.get("coarse_stride"),
        "scan_outer": settings.get("scan_outer"),
        "refine_outer": settings.get("refine_outer"),
        "remez_iters": settings.get("remez_iters"),
        "cheb_degree_mode": settings.get("cheb_degree_mode", "records"),
        "rat_k_max": int(settings.get("rat_k_max", 8)),
    }
    return checks == current


def record_rows(path: Path, function):
    """The rows of one record file, or [] when it is absent or unreadable.

    A record that cannot be parsed is ignored with a warning rather than
    taken as empty: the other copy of the record is still read.
    """
    if not path.exists():
        return []
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except Exception as exc:
        print(f"warning: cannot read {path} ({exc}); ignoring it", flush=True)
        return []
    if data.get("meta", {}).get("schema") != cs.ROW_SCHEMA:
        print(f"warning: {path} has another schema; ignoring it", flush=True)
        return []
    if data.get("meta", {}).get("function") != function:
        print(f"warning: {path} was written for "
              f"{data.get('meta', {}).get('function')!r}; ignoring it",
              flush=True)
        return []
    return data.get("rows", [])


def load_done(path: Path, current, function):
    """Completed rows from the primary record and its mirror.

    write_outputs mirrors every save to sweep_all.backup.json through a
    temp-and-replace, so a kill during the primary write leaves the mirror
    one row behind and never loses the run.
    """
    backup = path.with_name(path.name.replace(".json", ".backup.json"))
    done = {}
    for row in record_rows(backup, function) + record_rows(path, function):
        if not (row.get("cheb") and row.get("rat")):
            continue
        edge = row.get("edge")
        if edge not in EDGES:
            continue
        a_text = row.get("a")
        if not row_identity_ok(row, edge, a_text, function, current):
            continue
        done[(edge, row["target"], row["width_name"])] = row
    return done


def memoized_run_jobs(g, orig):
    """Memoize the generator's solver dispatch at the job boundary.

    Only rat_candidate_job is memoized, and its job spec determines the solve
    completely: (root, m, k, indices, ts text, ref text, outer, iters).  The
    worker rebuilds the mpf values at RAT_DPS from the text, so two equal
    specs compute the same two candidates.  Results come back in submission
    order, as the generator's runner returns them.  ``orig`` is the
    unpatched runner; the caller restores it after the cell, so wrappers
    never nest.
    """
    cache: dict = {}
    stats = {"calls": 0, "jobs": 0, "hits": 0}

    def key_of(job):
        root, m, k, indices, ts_text, ref_text, outer, iters = job
        return (root, m, k, tuple(indices), tuple(ts_text), tuple(ref_text),
                outer, iters)

    def wrapped(fn, jobs):
        if fn is not cs.rat_candidate_job:
            return orig(fn, jobs)
        stats["calls"] += 1
        stats["jobs"] += len(jobs)
        results = [None] * len(jobs)
        pending = []
        for i, job in enumerate(jobs):
            key = key_of(job)
            if key in cache:
                results[i] = cache[key]
                stats["hits"] += 1
            else:
                pending.append((i, key, job))
        if pending:
            outs = orig(fn, [job for _, _, job in pending])
            for (i, key, _), out in zip(pending, outs):
                cache[key] = out
                results[i] = out
        return results

    g.run_jobs = wrapped
    return stats, cache


def make_row(g, a, width_name, width, target, cheb, rat, edge, function,
             args, k_max):
    b = a + width
    row = {
        "precision": "fp64",
        "function": function,
        "edge": edge,
        "target": cs.target_tag(target),
        "target_kind": target["kind"],
        "target_value": target["value"],
        "width_name": width_name,
        "width": g.mp.nstr(width, 17),
        "width_f": float(width),
        "a": g.mp.nstr(a, 17),
        "b": g.mp.nstr(b, 17),
        "cheb": cheb,
        "rat": rat,
        "saving": (cheb["stored"] - rat["stored"])
        if cheb.get("stored") is not None and rat.get("stored") is not None
        else None,
        "settings": {
            "nodes": args.nodes + 1,
            "coarse_stride": args.stride,
            "scan_outer": args.scan_outer,
            "refine_outer": args.refine_outer,
            "remez_iters": args.remez_iters,
            "cheb_degree_mode": args.cheb_degrees,
            "rat_k_max": int(k_max),
            "workers": g.worker_count(),
            "solver_memo": True,
        },
    }
    return row


def write_csv_ext(path: Path, rows):
    """The paper's CSV columns plus function, edge and width_x0.

    count_sweep.write_csv reads a fixed field list from schema 2; a row from
    this lane carries three more fields, so the writer here repeats the
    paper's column names in the paper's order and appends them.
    """
    import csv as _csv
    fields = [
        "precision", "function", "edge", "target", "target_kind",
        "target_value", "width_name", "width", "width_x0", "width_f",
        "family", "stored", "m", "k", "solver", "delivered_rel",
        "delivered_abs", "at_x_rel", "at_x_abs", "refined", "below_meets",
        "seconds", "status",
        "nodes", "coarse_stride", "scan_outer", "refine_outer", "remez_iters",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = _csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            settings = row.get("settings", {})
            for family in ("cheb", "rat"):
                entry = row[family]
                writer.writerow({
                    "precision": row.get("precision"),
                    "function": row.get("function"),
                    "edge": row.get("edge"),
                    "target": row["target"],
                    "target_kind": row["target_kind"],
                    "target_value": row["target_value"],
                    "width_name": row["width_name"],
                    "width": row["width"],
                    "width_x0": row.get("width_x0"),
                    "width_f": row["width_f"],
                    "family": entry["family"],
                    "stored": entry.get("stored"),
                    "m": entry.get("m"),
                    "k": entry.get("k"),
                    "solver": entry.get("solver"),
                    "delivered_rel": entry.get("delivered_rel"),
                    "delivered_abs": entry.get("delivered_abs"),
                    "at_x_rel": entry.get("at_x_rel"),
                    "at_x_abs": entry.get("at_x_abs"),
                    "refined": entry.get("refined"),
                    "below_meets": entry.get("below_meets"),
                    "seconds": entry.get("seconds"),
                    "status": entry.get("status"),
                    "nodes": settings.get("nodes"),
                    "coarse_stride": settings.get("coarse_stride"),
                    "scan_outer": settings.get("scan_outer"),
                    "refine_outer": settings.get("refine_outer"),
                    "remez_iters": settings.get("remez_iters"),
                })


def write_outputs(out_dir: Path, rows, targets, meta):
    """sweep_all.json plus per-edge and per-edge-per-target CSVs.

    The primary record is mirrored to sweep_all.backup.json through a
    temp-and-replace: a kill during a write never leaves both copies short.
    """
    payload = {"meta": meta, "precision_scope": cs.PRECISION_SCOPE,
               "rows": rows}
    cs.write_json(out_dir / "sweep_all.json", payload)
    tmp = out_dir / "sweep_all.backup.json.tmp"
    cs.write_json(tmp, payload)
    os.replace(tmp, out_dir / "sweep_all.backup.json")
    for edge in EDGES:
        erows = [r for r in rows if r["edge"] == edge]
        if not erows:
            continue
        write_csv_ext(out_dir / f"sweep_{edge}.csv", erows)
        for target in targets:
            tag = cs.target_tag(target).replace(":", "_")
            trows = [r for r in erows if r["target"] == cs.target_tag(target)]
            if trows:
                cs.write_json(out_dir / f"sweep_{edge}_{tag}.json",
                              {"meta": meta,
                               "precision_scope": cs.PRECISION_SCOPE,
                               "rows": trows})
                write_csv_ext(out_dir / f"sweep_{edge}_{tag}.csv", trows)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--boys-root", default=str(DEFAULT_BOYS))
    parser.add_argument("--out-dir", default=str(DEFAULT_OUT))
    parser.add_argument("--function", choices=("F1", "F0"), default="F1",
                        help="F1 measures the second instance; F0 is the "
                             "instrument control against the paper's records")
    parser.add_argument("--targets", default=DEFAULT_TARGETS)
    parser.add_argument("--widths", default="all",
                        help="all, or a comma list of Group A names")
    parser.add_argument("--edges", default="x0,a2",
                        help="comma list of x0 and/or a2")
    parser.add_argument("--jobs", type=int, default=4,
                        help="solver workers; the lane's cap is 4")
    parser.add_argument("--nodes", type=int, default=cs.N_NODES)
    parser.add_argument("--stride", type=int, default=cs.COARSE_STRIDE)
    parser.add_argument("--scan-outer", type=int, default=cs.RAT_SCAN_OUTER)
    parser.add_argument("--refine-outer", type=int, default=cs.RAT_REFINE_OUTER)
    parser.add_argument("--remez-iters", type=int, default=cs.RAT_REMEZ_ITERS)
    parser.add_argument("--cheb-degrees", choices=("records", "all"),
                        default="records")
    parser.add_argument("--rat-k-max", type=int, default=0)
    parser.add_argument("--cells", default="",
                        help="comma list of edge/target/width keys, e.g. "
                             "x0/rel:1e-13/2^-7")
    parser.add_argument("--max-seconds", type=float, default=0,
                        help="stop before starting a new target row after "
                             "this long")
    parser.add_argument("--progress", action="store_true")
    args = parser.parse_args(argv)

    boys = Path(args.boys_root)
    g = cs.load_generator(boys)
    g.JOBS = args.jobs
    cs.BOYS_ROOT_STR = str(boys)
    if args.function == "F1":
        cs.cheb_fit = f1r.f1_cheb_fit

    targets = cs.parse_targets(args.targets)
    all_widths = cs.width_set(g, "all")
    if args.widths == "all":
        widths = all_widths
    else:
        wanted = [t.strip() for t in args.widths.split(",") if t.strip()]
        names = {name for name, _ in all_widths}
        unknown = [t for t in wanted if t not in names]
        if unknown:
            raise SystemExit(f"unknown Group A width names: {unknown}")
        widths = [(n, w) for n, w in all_widths if n in wanted]
    edges = [e.strip() for e in args.edges.split(",") if e.strip()]
    for edge in edges:
        if edge not in EDGES:
            raise SystemExit(f"unknown edge {edge!r}; use x0 or a2")

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    all_path = out_dir / "sweep_all.json"
    current = {"nodes": args.nodes + 1, "coarse_stride": args.stride,
               "scan_outer": args.scan_outer,
               "refine_outer": args.refine_outer,
               "remez_iters": args.remez_iters,
               "cheb_degree_mode": args.cheb_degrees,
               "rat_k_max": int(args.rat_k_max) if args.rat_k_max else int(g.RAT_K_MAX)}
    done = load_done(all_path, current, args.function)
    meta = build_meta(boys, args, g)

    order_names = [n for n in ORDER_WIDTHS if n in {w for w, _ in widths}]
    if args.cells:
        cells = []
        for token in args.cells.split(","):
            token = token.strip()
            if not token:
                continue
            parts = token.split("/")
            if len(parts) != 3:
                raise SystemExit(f"cell must be edge/target/width: {token}")
            edge, target_spec, name = parts
            target = next((t for t in targets if cs.target_tag(t) == target_spec),
                          None)
            if target is None:
                raise SystemExit(f"cell target not in --targets: {target_spec}")
            if edge not in edges:
                raise SystemExit(f"cell edge not requested: {edge}")
            width = next((w for n, w in widths if n == name), None)
            if width is None:
                raise SystemExit(f"cell width not in this grid: {name}")
            cells.append((edge, name, width, target))
        # group by edge+width so the memo sees all requested targets together
        cell_order = []
        seen = set()
        for edge, name, width, target in cells:
            if (edge, name) not in seen:
                seen.add((edge, name))
                cell_order.append((edge, name, width,
                                   [t for e, n, w, t in cells
                                    if e == edge and n == name]))
    else:
        cell_order = []
        for name in order_names:
            width = next(w for n, w in widths if n == name)
            for edge in edges:
                cell_order.append((edge, name, width, targets))

    print(f"function {args.function}; widths {len(widths)}; edges "
          f"{','.join(edges)}; targets {len(targets)}", flush=True)
    print(f"grid: {len(cell_order)} cell(s), {len(done)} row(s) already on "
          f"disk", flush=True)
    print(f"workers: {g.worker_count()}", flush=True)
    started = time.perf_counter()
    remaining = []
    ran = 0
    base_run_jobs = g.run_jobs
    for edge, name, width, cell_targets in cell_order:
        scale = EDGES[edge]
        a = g.X0 * scale
        a_text = g.mp.nstr(a, 17)
        w_edge = width * scale
        todo = []
        for target in cell_targets:
            key = (edge, cs.target_tag(target), name)
            if key in done:
                print(f"reuse {args.function} {edge}/{cs.target_tag(target)}/"
                      f"{name} cheb={done[key]['cheb'].get('stored')} "
                      f"rat={done[key]['rat'].get('stored')}", flush=True)
                continue
            todo.append(target)
        if not todo:
            continue
        if args.max_seconds and time.perf_counter() - started >= args.max_seconds:
            remaining.extend((edge, cs.target_tag(t), name) for t in todo)
            continue
        b = a + w_edge
        try:
            ts, xs = cs.make_nodes(g, a, b, args.nodes)
            if args.function == "F1":
                ref = f1r.f1_reference(g, xs)
            else:
                ref = cs.reference(g, xs)
            coarse_idx = list(range(0, args.nodes + 1, args.stride))
        except Exception as exc:
            import traceback
            print(f"ERROR setup {edge}/{name}: {type(exc).__name__}: {exc}",
                  flush=True)
            traceback.print_exc()
            continue
        stats, cache = memoized_run_jobs(g, base_run_jobs)
        for target in todo:
            tag = cs.target_tag(target)
            cell_started = time.perf_counter()
            try:
                cheb = cs.cheb_count(g, a, b, xs, ref, target,
                                     degrees=cs.cheb_degrees(args.cheb_degrees),
                                     verbose=args.progress)
                rat = cs.rational_count(
                    g, a, b, ts, xs, ref, target, coarse_idx,
                    scan_outer=args.scan_outer,
                    refine_outer=args.refine_outer,
                    remez_iters=args.remez_iters,
                    k_max=current["rat_k_max"],
                    verbose=args.progress)
            except Exception as exc:
                import traceback
                print(f"ERROR {edge}/{tag}/{name}: {type(exc).__name__}: {exc}",
                      flush=True)
                traceback.print_exc()
                continue
            row = make_row(g, a, name, w_edge, target, cheb, rat, edge,
                           args.function, args, current["rat_k_max"])
            row["width_x0"] = g.mp.nstr(width, 17)
            row["parity"] = parity_of(row)
            done[(edge, tag, name)] = row
            ran += 1
            write_outputs(out_dir, list(done.values()), targets, meta)
            print(f"[{args.function} {edge} {tag}] {name:<7} "
                  f"W={row['width_f']:<12.6g} "
                  f"cheb={row['cheb'].get('stored')} "
                  f"(deg={row['cheb'].get('m')}, "
                  f"{row['parity']['cheb_degree_parity']}) "
                  f"rat={row['rat'].get('stored')} "
                  f"([{row['rat'].get('m')}/{row['rat'].get('k')}], "
                  f"{row['rat'].get('solver')}, below="
                  f"{row['rat'].get('below_meets')}) "
                  f"rel={row['rat'].get('delivered_rel')} "
                  f"({row['cheb'].get('seconds')}s + {row['rat'].get('seconds')}s)"
                  f"  memo {stats['hits']}/{stats['jobs']} jobs reused"
                  f"  cell {time.perf_counter() - cell_started:.1f}s",
                  flush=True)
        cache.clear()
        g.run_jobs = base_run_jobs
    if ran:
        write_outputs(out_dir, list(done.values()), targets, meta)
    if remaining:
        print(f"\n{len(remaining)} row(s) remain in this grid: "
              f"{', '.join(f'{e}/{t}/{w}' for e, t, w in remaining)}", flush=True)
        print("re-run the same command; written rows are reused", flush=True)
    print(f"invocation wall time: {time.perf_counter() - started:.1f} s; "
          f"rows on disk: {len(done)}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
