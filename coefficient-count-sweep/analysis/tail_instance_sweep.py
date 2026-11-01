#!/usr/bin/env python3
"""Measure Paper C's Group A counts for the pure algebraic tail f(x) = x^-1/2.

The referees' clearing test for the width law is a second instance of the
tail class: the paper measures Boys F_0, whose tail is x^-1/2 but which also
carries an entire erfc transient.  The pure tail is exactly in the class the
model states (c = 1, alpha = 1/2, an analytic non-vanishing g, nearest
non-analyticity the branch at the origin) with no transient.

Nothing in the fitting, the acceptance test, the grid or the search space is
reimplemented here.  The driver imports the paper's own harness
(count_sweep.py) and the kernel's generator, and replaces exactly one thing:
the generator's reference routine, ``boys_ref``, with the exact tail
1/sqrt(x) at 60 digits (mpmath, computed 20 digits above the working
precision).  The generator's ``cheb_coeffs`` builds the Chebyshev interpolant
from that module-level reference, so the substitution reaches it; the
rational solvers (``rat_remez``, ``rat_lawson``) take the grid and the
reference values as arguments and are untouched, and their workers receive
the reference values as text, so the worker processes solve the same problem
without needing the substitution.

Rows are written in the paper's row schema (schema 2, meta, precision_scope,
rows), with the per-row fields the paper's records carry, extended by
``function``, ``edge``, ``edge_scale``, ``width_x0`` (the undoubled Group A
width a row is paired with) and ``parity``.

Resumable: every completed row is written to disk, and a re-run with the same
arguments reuses rows whose function, left edge and settings match.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path
from types import SimpleNamespace

PAPER = Path(__file__).resolve().parent
sys.path.insert(0, str(PAPER))

import count_sweep as cs  # noqa: E402

DEFAULT_BOYS = Path(r"C:\Users\myaml\dev\_mainonly")
DEFAULT_OUT = PAPER / "results-second-instance"
DEFAULT_TARGETS = "rel:1e-13,abs:5e-14,rel:1e-14"

FUNCTION_TAIL = "x^-1/2"
FUNCTION_F0 = "boys_F0"


def sha256(path: Path) -> str:
    return cs.sha256(path)


def install_tail(g):
    """Replace the generator's Boys reference with the exact tail 1/sqrt(x).

    The replacement is cached exactly as the generator's own routine is, so
    no node's value is computed twice, and it accepts the same four
    arguments the harness passes (F_0's order, the argument, and the series
    floor and term cap, which a closed form does not use).
    """
    mp = g.mp
    dps = g.RAT_DPS

    def tail_ref(n, x, tail_floor=None, max_terms=None):
        del n, tail_floor, max_terms
        with mp.workdps(dps + 20):
            return 1 / mp.sqrt(mp.mpf(x))

    g.boys_ref = g.lru_cache(maxsize=1 << 20)(tail_ref)
    return g.boys_ref


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


def build_meta(boys: Path, g, args, a, function):
    ns = SimpleNamespace(nodes=args.nodes, stride=args.stride,
                         scan_outer=args.scan_outer,
                         refine_outer=args.refine_outer,
                         remez_iters=args.remez_iters,
                         cheb_degrees=args.cheb_degrees,
                         rat_k_max=args.rat_k_max)
    meta = cs.sweep_meta(boys, g, ns)
    meta["function"] = function
    meta["left_edge"] = g.mp.nstr(a, 17)
    meta["left_edges"] = [g.mp.nstr(a, 17)]
    meta["edge_scales"] = [args.scale]
    meta["base_left_edge"] = g.mp.nstr(g.X0, 17)
    meta["width_scale"] = args.scale
    meta["matched_ratio"] = ("every width is the Group A width multiplied by "
                             f"{args.scale}, so W/a is the Group A W/x0")
    meta["reference"] = (
        "1/sqrt(x) at 60-digit working precision (mpmath.sqrt, 20 digits "
        "above the working precision); the generator's boys_ref is replaced "
        "by this exact value and its cheb_coeffs builds the interpolant "
        "from it. The rational solvers are the generator's, unchanged."
        if function == FUNCTION_TAIL else
        "the generator's own boys_ref(0, x, 1e-55, 900); the control run")
    meta["harness"] = {
        "count_sweep": {"path": str(PAPER / "count_sweep.py"),
                        "sha256": sha256(PAPER / "count_sweep.py")},
        "driver": {"path": str(PAPER / "tail_instance_sweep.py"),
                   "sha256": sha256(PAPER / "tail_instance_sweep.py")},
    }
    meta["target_note"] = (
        "rel targets are relative to f; abs:5e-14 is an absolute error on f, "
        "unscaled (the paper's absolute arm is F_0(a)-scaled through L; the "
        "analysis states the model convention it uses)")
    return meta


def row_identity_ok(row, function, current):
    """Whether a stored row is the row this run would produce, at its own edge."""
    if not row.get("a"):
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
    if checks != current:
        return False
    width_x0 = row.get("width_x0")
    try:
        scale = float(row.get("edge_scale", 1))
    except (TypeError, ValueError):
        return False
    if width_x0 is None:
        return False
    if float(width_x0) * scale != float(row["width"]):
        return False
    return True


def load_done(path: Path, function, current):
    """Every completed row already written for this function, at either edge.

    The file holds both edges (the paper's schema carries the left edge in
    every row), so identity is per row: the function, the row's own left
    edge and scale, the settings, and ``width_x0 * scale == width``.
    """
    if not path.exists():
        return {}
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except Exception as exc:
        print(f"warning: cannot read {path} ({exc}); starting fresh", flush=True)
        return {}
    if data.get("meta", {}).get("schema") != cs.ROW_SCHEMA:
        print(f"warning: {path} has another schema; starting fresh", flush=True)
        return {}
    if data.get("meta", {}).get("function") != function:
        print(f"warning: {path} was written for function "
              f"{data.get('meta', {}).get('function')!r}; starting fresh",
              flush=True)
        return {}
    done = {}
    for row in data.get("rows", []):
        if not (row.get("cheb") and row.get("rat")):
            continue
        if not row_identity_ok(row, function, current):
            continue
        done[(row["a"], row["target"], row["width_name"])] = row
    return done


def load_done_cells(path: Path, function, current):
    """Every Chebyshev cell already written for this function, at either edge."""
    if not path.exists():
        return {}
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except Exception as exc:
        print(f"warning: cannot read {path} ({exc}); starting fresh", flush=True)
        return {}
    meta = data.get("meta", {})
    if meta.get("schema") != cs.ROW_SCHEMA or meta.get("function") != function:
        print(f"warning: {path} is another run's; starting fresh", flush=True)
        return {}
    done = {}
    for cell in data.get("cells", []):
        if cell.get("function") != function or not cell.get("a"):
            continue
        settings = cell.get("settings", {})
        if settings.get("nodes") != current["nodes"]:
            continue
        if settings.get("cheb_degree_mode", "records") != current["cheb_degree_mode"]:
            continue
        if cell.get("cheb", {}).get("stored") is None:
            continue
        done[(cell["a"], cell["target"], cell["width_name"])] = cell
    return done


def refresh_meta(meta, done, cells_done, a_text, scale):
    """The file holds every edge measured so far; the meta says which."""
    edges = {row.get("a") for row in done.values()} | \
            {cell.get("a") for cell in cells_done.values()} | {a_text}
    scales = {int(row.get("edge_scale")) for row in done.values()} | \
             {int(cell.get("edge_scale")) for cell in cells_done.values()} | {scale}
    meta["left_edges"] = sorted(edge for edge in edges if edge)
    meta["edge_scales"] = sorted(scales)
    return meta


def write_cells(path: Path, cells, meta):
    doc = {"meta": meta, "precision_scope": cs.PRECISION_SCOPE,
           "cells": cells}
    cs.write_json(path, doc)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--boys-root", default=str(DEFAULT_BOYS))
    parser.add_argument("--out-dir", default=str(DEFAULT_OUT))
    parser.add_argument("--function", choices=("tail", "f0"), default="tail",
                        help="tail: f(x)=x^-1/2 (the measurement); f0: the "
                             "generator's own Boys reference (the control)")
    parser.add_argument("--family", choices=("both", "cheb", "rat"),
                        default="both")
    parser.add_argument("--targets", default=DEFAULT_TARGETS)
    parser.add_argument("--widths", default="all",
                        help="all, or a comma list of Group A names")
    parser.add_argument("--scale", type=int, default=1,
                        help="a = scale * X0 and W' = scale * W (1 is the "
                             "paper's first edge, 2 the second)")
    parser.add_argument("--jobs", type=int, default=4,
                        help="worker cap; the lane's cap is 4")
    parser.add_argument("--nodes", type=int, default=cs.N_NODES)
    parser.add_argument("--stride", type=int, default=cs.COARSE_STRIDE)
    parser.add_argument("--scan-outer", type=int, default=cs.RAT_SCAN_OUTER)
    parser.add_argument("--refine-outer", type=int, default=cs.RAT_REFINE_OUTER)
    parser.add_argument("--remez-iters", type=int, default=cs.RAT_REMEZ_ITERS)
    parser.add_argument("--cheb-degrees", choices=("records", "all"),
                        default="records")
    parser.add_argument("--rat-k-max", type=int, default=0)
    parser.add_argument("--cells", default="",
                        help="comma list of target/width_name keys to run "
                             "(default: every missing row)")
    parser.add_argument("--max-seconds", type=float, default=0,
                        help="stop before starting a new row after this long")
    parser.add_argument("--progress", action="store_true")
    args = parser.parse_args(argv)

    function = FUNCTION_TAIL if args.function == "tail" else FUNCTION_F0
    boys = Path(args.boys_root)
    g = cs.load_generator(boys)
    g.JOBS = args.jobs
    cs.BOYS_ROOT_STR = str(boys)
    if function == FUNCTION_TAIL:
        install_tail(g)
    a = g.X0 * args.scale
    a_text = g.mp.nstr(a, 17)
    targets = cs.parse_targets(args.targets)
    all_widths = cs.width_set(g, "all")
    if args.widths == "all":
        widths = all_widths
    else:
        wanted = [token.strip() for token in args.widths.split(",") if token.strip()]
        names = {name for name, _ in all_widths}
        unknown = [token for token in wanted if token not in names]
        if unknown:
            raise SystemExit(f"unknown Group A width names: {unknown}")
        widths = [(name, width) for name, width in all_widths if name in wanted]
    scaled = [(name, width * args.scale) for name, width in widths]

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    all_path = out_dir / "sweep_all.json"
    cheb_path = out_dir / "chebyshev_only.json"
    current = {"nodes": args.nodes + 1, "coarse_stride": args.stride,
               "scan_outer": args.scan_outer,
               "refine_outer": args.refine_outer,
               "remez_iters": args.remez_iters,
               "cheb_degree_mode": args.cheb_degrees,
               "rat_k_max": int(args.rat_k_max) if args.rat_k_max else int(g.RAT_K_MAX)}
    done = load_done(all_path, function, current)
    cells_done = load_done_cells(cheb_path, function, current)
    meta = build_meta(boys, g, args, a, function)
    print(f"on disk: {len(done)} full row(s) (both edges), "
          f"{len(cells_done)} Chebyshev cell(s) (both edges)", flush=True)

    order = []
    if args.cells:
        for token in args.cells.split(","):
            token = token.strip()
            if not token:
                continue
            target_spec, name = token.rsplit("/", 1)
            target = next((t for t in targets if cs.target_tag(t) == target_spec),
                          None)
            if target is None:
                raise SystemExit(f"cell target not in --targets: {target_spec}")
            width = next((w for n, w in scaled if n == name), None)
            if width is None:
                raise SystemExit(f"cell width not in this grid: {name}")
            order.append((target, name, width))
    else:
        for target in targets:
            for name, width in scaled:
                order.append((target, name, width))

    print(f"function {function}; left edge a = {a_text} "
          f"(scale {args.scale} x X0); widths scaled by {args.scale}", flush=True)
    print(f"grid: {len(scaled)} width(s) x {len(targets)} target(s); "
          f"{len(done)} full row(s) and {len(cells_done)} Chebyshev cell(s) "
          f"already on disk", flush=True)
    print(f"workers: {g.worker_count()}", flush=True)
    started = time.perf_counter()
    remaining = []
    ran = 0
    ran_cells = 0
    if args.family in ("both", "cheb"):
        # One node grid and one reference per width, shared by its targets.
        by_width = {}
        for target, name, width in order:
            by_width.setdefault((name, width), []).append(target)
        for (name, width), ts_targets in by_width.items():
            missing = [t for t in ts_targets
                       if (a_text, cs.target_tag(t), name) not in cells_done]
            if not missing:
                print(f"reuse cheb {name}: "
                      f"{len(ts_targets)} target(s) on disk", flush=True)
                continue
            if args.max_seconds and time.perf_counter() - started >= args.max_seconds:
                remaining.extend((cs.target_tag(t), name) for t in missing)
                continue
            b = a + width
            ts, xs = cs.make_nodes(g, a, b, args.nodes)
            reference_started = time.perf_counter()
            ref = cs.reference(g, xs)
            reference_seconds = time.perf_counter() - reference_started
            for target in missing:
                cheb = cs.cheb_count(g, a, b, xs, ref, target,
                                     degrees=cs.cheb_degrees(args.cheb_degrees),
                                     verbose=args.progress)
                cell = {
                    "function": function,
                    "edge_scale": args.scale,
                    "target": cs.target_tag(target),
                    "target_kind": target["kind"],
                    "target_value": target["value"],
                    "width_name": name,
                    "width_x0": g.mp.nstr(width / args.scale, 17),
                    "width": g.mp.nstr(width, 17),
                    "width_f": float(width),
                    "a": a_text,
                    "b": g.mp.nstr(b, 17),
                    "cheb": cheb,
                    "parity": {"cheb_degree": cheb.get("m"),
                               "cheb_degree_parity": (
                                   "odd" if cheb.get("m") % 2 else "even")
                               if cheb.get("m") is not None else None},
                    "status": cheb.get("status"),
                    "reference_seconds": round(reference_seconds, 3),
                    "settings": {"nodes": args.nodes + 1,
                                 "cheb_degree_mode": args.cheb_degrees},
                }
                cells_done[(a_text, cs.target_tag(target), name)] = cell
                ran_cells += 1
                write_cells(cheb_path, list(cells_done.values()),
                            refresh_meta(meta, done, cells_done, a_text,
                                         args.scale))
                print(f"cheb [{cell['target']}] {name:<7} W={float(width):<12.6g} "
                      f"stored={cheb.get('stored')} (deg={cheb.get('m')}) "
                      f"rel={cheb.get('delivered_rel')} "
                      f"abs={cheb.get('delivered_abs')} ({cheb.get('seconds')}s)",
                      flush=True)
    if args.family in ("both", "rat"):
        for target, name, width in order:
            tag = cs.target_tag(target)
            key = (a_text, tag, name)
            if key in done:
                print(f"reuse {tag}/{name} cheb={done[key]['cheb'].get('stored')} "
                      f"rat={done[key]['rat'].get('stored')}", flush=True)
                continue
            if args.max_seconds and time.perf_counter() - started >= args.max_seconds:
                remaining.append((tag, name))
                continue
            try:
                row = cs.sweep_one(g, a, name, width, target, args.nodes,
                                   args.stride,
                                   scan_outer=args.scan_outer,
                                   refine_outer=args.refine_outer,
                                   remez_iters=args.remez_iters,
                                   cheb_degree_mode=args.cheb_degrees,
                                   k_max=current["rat_k_max"],
                                   verbose=args.progress)
            except Exception as exc:  # one bad cell must not end the sweep
                import traceback
                print(f"ERROR {tag}/{name}: {type(exc).__name__}: {exc}",
                      flush=True)
                traceback.print_exc()
                continue
            row["function"] = function
            row["edge"] = "x0" if args.scale == 1 else f"{args.scale}x0"
            row["edge_scale"] = args.scale
            row["width_x0"] = g.mp.nstr(width / args.scale, 17)
            row["parity"] = parity_of(row)
            done[key] = row
            ran += 1
            cs.write_sweep_outputs(out_dir, list(done.values()), targets,
                                   refresh_meta(meta, done, cells_done,
                                                a_text, args.scale))
            print(f"[{tag}] {name:<7} W={row['width_f']:<12.6g} "
                  f"cheb={row['cheb'].get('stored')} "
                  f"(deg={row['cheb'].get('m')}, "
                  f"{row['parity']['cheb_degree_parity']}) "
                  f"rat={row['rat'].get('stored')} "
                  f"([{row['rat'].get('m')}/{row['rat'].get('k')}], "
                  f"{row['rat'].get('solver')}) "
                  f"rel={row['rat'].get('delivered_rel')} "
                  f"({row['cheb'].get('seconds')}s + "
                  f"{row['rat'].get('seconds')}s)", flush=True)
            # The accepted Chebyshev count is the same object cheb_count
            # produced; record it in the cell file too.
            if (a_text, tag, name) not in cells_done:
                cells_done[(a_text, tag, name)] = {
                    "function": function, "edge_scale": args.scale,
                    "target": tag, "target_kind": target["kind"],
                    "target_value": target["value"],
                    "width_name": name,
                    "width_x0": row["width_x0"], "width": row["width"],
                    "width_f": row["width_f"], "a": row["a"], "b": row["b"],
                    "cheb": row["cheb"],
                    "parity": {"cheb_degree": row["cheb"].get("m"),
                               "cheb_degree_parity":
                                   row["parity"]["cheb_degree_parity"]},
                    "reference_seconds": None,
                    "settings": {"nodes": args.nodes + 1,
                                 "cheb_degree_mode": args.cheb_degrees},
                }
                write_cells(cheb_path, list(cells_done.values()),
                            refresh_meta(meta, done, cells_done, a_text,
                                         args.scale))
    if ran:
        cs.write_sweep_outputs(out_dir, list(done.values()), targets,
                               refresh_meta(meta, done, cells_done, a_text,
                                            args.scale))
    if remaining:
        print(f"\n{len(remaining)} row(s) remain in this grid: "
              f"{', '.join(f'{t}/{w}' for t, w in remaining)}", flush=True)
        print("re-run the same command; written rows are reused", flush=True)
    print(f"invocation wall time: {time.perf_counter() - started:.1f} s; "
          f"full rows on disk: {len(done)}; Chebyshev cells on disk: "
          f"{len(cells_done)}; this invocation: {ran} row(s), "
          f"{ran_cells} cell(s)", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
