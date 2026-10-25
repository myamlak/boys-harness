#!/usr/bin/env python3
"""Re-run Paper C's Group A width sweep at a second left edge.

The referee's clearing test for the W/a-alone claim is to leave the ratio
W/a fixed while moving the left edge.  This driver takes the paper's own
harness (count_sweep.py) and calls its sweep_one with a left edge
``a = SCALE * X0`` and every Group A width multiplied by ``SCALE``.  The
default is SCALE = 2: ``a2 = 2 * x0 = 23.799696304216968``.

Nothing in the fit, the reference, the target or the acceptance test is
reimplemented here.  The rows are written in the shape of
results-ladder/sweep_all.json (meta, precision_scope, rows), with the
per-row fields the paper's records already carry, plus the two provenance
fields a second-edge comparison needs: ``width_x0`` (the undoubled width the
row is paired with) and ``parity``.

Resumable: every completed row is written to disk, and a re-run with the
same arguments reuses rows whose left edge and settings match.
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
DEFAULT_OUT = PAPER / "results-second-edge"
DEFAULT_TARGETS = "rel:1e-13,abs:5e-14,rel:1e-14"

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


def build_meta(boys: Path, g, args, a):
    ns = SimpleNamespace(nodes=args.nodes, stride=args.stride,
                         scan_outer=args.scan_outer,
                         refine_outer=args.refine_outer,
                         remez_iters=args.remez_iters,
                         cheb_degrees=args.cheb_degrees,
                         rat_k_max=args.rat_k_max)
    meta = cs.sweep_meta(boys, g, ns)
    meta["left_edge"] = g.mp.nstr(a, 17)
    meta["base_left_edge"] = g.mp.nstr(g.X0, 17)
    meta["width_scale"] = args.scale
    meta["matched_ratio"] = ("every width is the Group A width multiplied by "
                             f"{args.scale}, so W/a is the Group A W/x0")
    meta["widths_note"] = ("width is the scaled interval length; width_x0 is "
                           "the Group A width the row is paired with")
    return meta


def row_identity_ok(row, a_text, scale, current):
    if row.get("a") != a_text:
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
    if width_x0 is None:
        return False
    if float(width_x0) * scale != float(row["width"]):
        return False
    return True


def load_done(path: Path, a_text, scale, current):
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
    if data.get("meta", {}).get("left_edge") != a_text:
        print(f"warning: {path} was written at left edge "
              f"{data.get('meta', {}).get('left_edge')!r}, not {a_text!r}; "
              f"starting fresh", flush=True)
        return {}
    if str(data.get("meta", {}).get("width_scale")) != str(scale):
        print(f"warning: {path} was written at width scale "
              f"{data.get('meta', {}).get('width_scale')!r}; starting fresh",
              flush=True)
        return {}
    done = {}
    for row in data.get("rows", []):
        if not (row.get("cheb") and row.get("rat")):
            continue
        if not row_identity_ok(row, a_text, scale, current):
            continue
        done[(row["target"], row["width_name"])] = row
    return done


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--boys-root", default=str(DEFAULT_BOYS))
    parser.add_argument("--out-dir", default=str(DEFAULT_OUT))
    parser.add_argument("--targets", default=DEFAULT_TARGETS)
    parser.add_argument("--widths", default="all",
                        help="all, or a comma list of Group A names")
    parser.add_argument("--scale", type=int, default=2,
                        help="a = scale * X0 and W' = scale * W (2 is the "
                             "second left edge; 1 is an x0 control)")
    parser.add_argument("--jobs", type=int, default=4,
                        help="rational-scan workers; the lane's cap is 4")
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

    boys = Path(args.boys_root)
    g = cs.load_generator(boys)
    g.JOBS = args.jobs
    cs.BOYS_ROOT_STR = str(boys)
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
    current = {"nodes": args.nodes + 1, "coarse_stride": args.stride,
               "scan_outer": args.scan_outer,
               "refine_outer": args.refine_outer,
               "remez_iters": args.remez_iters,
               "cheb_degree_mode": args.cheb_degrees,
               "rat_k_max": int(args.rat_k_max) if args.rat_k_max else int(g.RAT_K_MAX)}
    done = load_done(all_path, a_text, args.scale, current)
    meta = build_meta(boys, g, args, a)

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

    print(f"left edge a = {a_text} (scale {args.scale} x X0); "
          f"widths scaled by {args.scale}", flush=True)
    print(f"grid: {len(scaled)} width(s) x {len(targets)} target(s); "
          f"{len(done)} row(s) already on disk", flush=True)
    print(f"workers: {g.worker_count()}", flush=True)
    started = time.perf_counter()
    remaining = []
    ran = 0
    for target, name, width in order:
        tag = cs.target_tag(target)
        key = (tag, name)
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
        row["width_x0"] = g.mp.nstr(width / args.scale, 17)
        row["parity"] = parity_of(row)
        done[key] = row
        ran += 1
        cs.write_sweep_outputs(out_dir, list(done.values()), targets, meta)
        print(f"[{tag}] {name:<7} W={row['width_f']:<12.6g} "
              f"cheb={row['cheb'].get('stored')} "
              f"(deg={row['cheb'].get('m')}, "
              f"{row['parity']['cheb_degree_parity']}) "
              f"rat={row['rat'].get('stored')} "
              f"([{row['rat'].get('m')}/{row['rat'].get('k')}], "
              f"{row['rat'].get('solver')}) "
              f"rel={row['rat'].get('delivered_rel')} "
              f"({row['cheb'].get('seconds')}s + {row['rat'].get('seconds')}s)",
              flush=True)
    if ran:
        cs.write_sweep_outputs(out_dir, list(done.values()), targets, meta)
    if remaining:
        print(f"\n{len(remaining)} row(s) remain in this grid: "
              f"{', '.join(f'{t}/{w}' for t, w in remaining)}", flush=True)
        print("re-run the same command; written rows are reused", flush=True)
    print(f"invocation wall time: {time.perf_counter() - started:.1f} s; "
          f"rows on disk: {len(done)}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
