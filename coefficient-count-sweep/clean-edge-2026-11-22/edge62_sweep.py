"""Run the dense eight-target ladder at a clean left edge, a = 62.

The recorded harness and the v5.0.0 kernel generator are used unchanged; the
only difference from the deposited x0 ladder is that every width is scaled by
a/x0, so W/a matches the Group A design while the analytic transient lies
below every target. Rows are written after each completion, and a re-run
reuses finished rows.

    python -I edge62_sweep.py --boys-root <kernel> --harness <dir> \
        --out <sweep.json> [--edge 62] [--jobs 12]
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

TARGETS = ("rel:1e-8", "rel:1e-9", "rel:1e-10", "rel:1e-11",
           "rel:1e-12", "rel:1e-13", "rel:1e-14", "rel:1e-15")
WIDTHS = ("2^-7", "2^0", "2^5", "ship", "paper")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--boys-root", required=True)
    ap.add_argument("--harness", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--edge", type=str, default="62")
    ap.add_argument("--jobs", type=int, default=12)
    ap.add_argument("--only", default="", help="target:width_name, for a pilot")
    args = ap.parse_args(argv)

    boys = Path(args.boys_root)
    harness = Path(args.harness)
    out = Path(args.out)
    sys.path.insert(0, str(harness))
    import count_sweep_recorded as cs

    g = cs.load_generator(boys)
    g.JOBS = args.jobs
    cs.BOYS_ROOT_STR = str(boys)
    a = g.mpf(args.edge)
    scale = a / g.X0
    base_widths = dict(cs.width_set(g, "all"))

    data = {"meta": {}, "rows": []}
    if out.exists():
        data = json.loads(out.read_text(encoding="utf-8"))
    done = {(row["target"], row["width_name"]) for row in data["rows"]}
    data["meta"] = {
        "schema": "edge-ladder/1",
        "edge": g.mp.nstr(a, 17),
        "scale": g.mp.nstr(scale, 17),
        "jobs": args.jobs,
        "targets": list(TARGETS),
        "widths": list(WIDTHS),
        "settings": {
            "nodes": cs.N_NODES,
            "coarse_stride": cs.COARSE_STRIDE,
            "scan_outer": cs.RAT_SCAN_OUTER,
            "refine_outer": cs.RAT_REFINE_OUTER,
            "remez_iters": cs.RAT_REMEZ_ITERS,
            "cheb_degree_mode": "all",
            "rat_k_max": 12,
        },
        "generator": {
            "path": str(boys / "tools" / "gen_boys_coefficients.py"),
            "x0": g.mp.nstr(g.X0, 17),
        },
    }

    targets = {spec: cs.parse_targets(spec)[0] for spec in TARGETS}
    for spec in TARGETS:
        for name in WIDTHS:
            key = (spec, name)
            if args.only and f"{spec}:{name}" != args.only:
                continue
            if key in done:
                print(f"skip {spec} {name}", flush=True)
                continue
            width = base_widths[name] * scale
            started = time.perf_counter()
            row = cs.sweep_one(
                g, a, name, width, targets[spec], cs.N_NODES, cs.COARSE_STRIDE,
                scan_outer=cs.RAT_SCAN_OUTER,
                refine_outer=cs.RAT_REFINE_OUTER,
                remez_iters=cs.RAT_REMEZ_ITERS,
                cheb_degree_mode="all", k_max=12, verbose=False)
            elapsed = time.perf_counter() - started
            row["width_x0"] = g.mp.nstr(base_widths[name], 17)
            row["edge"] = g.mp.nstr(a, 17)
            row["edge_wall_seconds"] = round(elapsed, 3)
            data["rows"].append(row)
            out.write_text(json.dumps(data, indent=1), encoding="utf-8")
            print(f"row {spec} {name} {elapsed:.1f}s stored={row['cheb']['stored']}"
                  f"/{row['rat']['stored']}", flush=True)
    print(f"done: {len(data['rows'])} rows", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
