#!/usr/bin/env python3
"""Bracket the two corrected F1 knee cells across four solve configurations.

The knee walk found each corrected cell with one solve per pair at the
recorded settings.  The paper's stability section shows that a solve at this
margin is configuration-dependent, so each cell is re-solved here in four
configurations at its accepted level and at the level below:

    base         2001 nodes, Remez 60, Lawson outer 300, the recorded nodes
    seed         base settings, the interval's left edge moved by 3.0e-32
    doubled      4001 nodes, Remez 120, Lawson outer 600
    doubled-seed doubled settings, the left edge moved by 3.0e-32

Every admissible pair is solved; each pair's fitted candidates (Remez and
Lawson) are evaluated on the full grid and the best delivered relative and
absolute errors are recorded whether or not the pair meets the target, so the
level's verdict and its margin are both visible.

Usage:
    python bracket_f1.py --cell a2-paper-rel13 --jobs 4
    python bracket_f1.py --cell a2-2^5-rel14 --jobs 4

Writes _knee-f1/bracket_<cell>.json after every pair.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import sys
import time
from pathlib import Path

PAPER = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PAPER))

import count_sweep as cs  # noqa: E402
import f1_reference as f1r  # noqa: E402

RECORD = PAPER / "results-second-instance-f1" / "sweep_all.json"
BOYS = Path(json.loads(RECORD.read_text(encoding="utf-8"))["meta"]["boys_root"])

CONFIGS = {
    "base": {"nodes": 2000, "refine_outer": 300, "remez_iters": 60, "shift": 0.0},
    "seed": {"nodes": 2000, "refine_outer": 300, "remez_iters": 60, "shift": 3.0e-32},
    "doubled": {"nodes": 4000, "refine_outer": 600, "remez_iters": 120, "shift": 0.0},
    "doubled-seed": {"nodes": 4000, "refine_outer": 600, "remez_iters": 120,
                     "shift": 3.0e-32},
}

CELLS = {
    "a2-paper-rel13": {
        "edge_scale": 2,
        "width_name": "paper",
        "target": "rel:1e-13",
        "accepted": 12,
        "below": 11,
    },
    "a2-2^5-rel14": {
        "edge_scale": 2,
        "width_name": "2^5",
        "target": "rel:1e-14",
        "accepted": 20,
        "below": 19,
    },
}


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def scan_level(g, a, b, ts, xs, ref, target, stored, outer, remez_iters):
    """Solve every (m,k) at one stored count; keep every pair's best error."""
    k_max = min(g.RAT_K_MAX, stored - 2)
    pairs = [(stored - 1 - k, k) for k in range(1, k_max + 1)]
    indices = list(range(len(ts)))
    ts_text = cs.ship_text(g, ts)
    ref_text = cs.ship_text(g, ref)
    specs = [(cs.BOYS_ROOT_STR, m, k, indices, ts_text, ref_text, outer,
              remez_iters) for m, k in pairs]
    results = g.run_jobs(cs.rat_candidate_job, specs)
    rows = []
    for (m, k), out in zip(pairs, results):
        best = None
        solvers = []
        for name, p, q in out:
            err = cs.rat_errors(g, p, q, a, b, xs, ref)
            if err is None:
                continue
            wr, wa, at_rel, at_abs = err
            met = cs.meets(target, wr, wa)
            solvers.append({"solver": name, "rel": g.mp.nstr(wr, 17),
                            "abs": g.mp.nstr(wa, 17), "meets": met,
                            "at_x_rel": g.mp.nstr(at_rel, 17)})
            key = wr if target["kind"] == "rel" else wa
            if best is None or key < best[0]:
                best = (key, name, met, wr, wa, at_rel, at_abs)
        entry = {"m": m, "k": k, "pair": f"[{m}/{k}]",
                 "solvers": solvers,
                 "best_solver": None if best is None else {
                     "solver": best[1],
                     "rel": g.mp.nstr(best[3], 17),
                     "abs": g.mp.nstr(best[4], 17),
                     "at_x_rel": g.mp.nstr(best[5], 17),
                     "meets": best[2],
                 }}
        rows.append(entry)
    return rows


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cell", required=True, choices=sorted(CELLS))
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--configs", default=",".join(CONFIGS))
    ap.add_argument("--smoke", action="store_true",
                    help="small grid and iteration counts, to check the plumbing")
    ap.add_argument("--out", default=None)
    args = ap.parse_args(argv)

    cell = CELLS[args.cell]
    g = cs.load_generator(BOYS)
    g.JOBS = args.jobs
    target = cs.parse_targets(cell["target"])[0]
    widths = dict(cs.width_set(g, "all"))
    width = widths[cell["width_name"]]
    scale = cell["edge_scale"]
    a = g.X0 * scale
    w_edge = width * scale
    b = a + w_edge

    out_path = Path(args.out) if args.out else (
        PAPER / "_knee-f1" / f"bracket_{args.cell}.json")
    doc = {
        "cell": args.cell,
        "script_digest": digest(Path(__file__)),
        "count_sweep_digest": digest(PAPER / "count_sweep.py"),
        "f1_reference_digest": digest(PAPER / "f1_reference.py"),
        "record": str(RECORD),
        "record_digest": digest(RECORD),
        "boys_root": str(BOYS),
        "python": sys.version,
        "platform": platform.platform(),
        "function": "F1",
        "edge": "a2",
        "target": target["spec"],
        "accepted": cell["accepted"],
        "below": cell["below"],
        "configs": [],
        "started_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    }

    def write():
        tmp = out_path.with_suffix(".json.tmp")
        tmp.write_text(json.dumps(doc, indent=1), encoding="utf-8")
        tmp.replace(out_path)

    write()
    print(f"cell {args.cell}: target {target['spec']}; a={g.mp.nstr(a, 17)} "
          f"b={g.mp.nstr(b, 17)}; accepted {cell['accepted']}, "
          f"below {cell['below']}", flush=True)

    for name in args.configs.split(","):
        cfg = dict(CONFIGS[name])
        if args.smoke:
            cfg.update(nodes=200, refine_outer=40, remez_iters=16)
        a_cfg = a + g.mpf(cfg["shift"])
        t0 = time.perf_counter()
        ts, xs = cs.make_nodes(g, a_cfg, b, cfg["nodes"])
        ref = f1r.f1_reference(g, xs)
        entry = {
            "config": name,
            "settings": {"nodes": cfg["nodes"] + 1,
                         "refine_outer": cfg["refine_outer"],
                         "remez_iters": cfg["remez_iters"],
                         "left_edge_shift": repr(cfg["shift"])},
            "a": g.mp.nstr(a_cfg, 17),
            "b": g.mp.nstr(b, 17),
            "reference_seconds": round(time.perf_counter() - t0, 3),
            "levels": [],
        }
        doc["configs"].append(entry)
        write()
        for level in (cell["accepted"], cell["below"]):
            t1 = time.perf_counter()
            rows = scan_level(g, a_cfg, b, ts, xs, ref, target, level,
                              cfg["refine_outer"], cfg["remez_iters"])
            secs = time.perf_counter() - t1
            hits = [r for r in rows
                    if r["best_solver"] is not None and r["best_solver"]["meets"]]
            lvl = {"stored": level, "seconds": round(secs, 3),
                   "n_pairs": len(rows), "meets": bool(hits),
                   "pairs": rows}
            entry["levels"].append(lvl)
            write()
            best = min((r for r in rows if r["best_solver"] is not None),
                       key=lambda r: float(r["best_solver"]["rel"]),
                       default=None)
            report = ("none" if best is None else
                      f"{best['pair']} {best['best_solver']['solver']} "
                      f"rel={best['best_solver']['rel']}")
            print(f"{name} level {level}: {'MEETS' if hits else 'misses'}; "
                  f"best {report} ({secs:.1f} s)", flush=True)

    doc["finished_utc"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    write()
    print(f"wrote {out_path}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
