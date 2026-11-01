#!/usr/bin/env python3
"""Walk the two F1 knee cells down with the paper's own search code.

The F1 sweep recorded two rational rows with below_meets=True: a stored count
whose count below also reaches the target.  The paper's rule makes the count
the smallest stored that reaches on the certifying grid, so each row needs the
levels below it solved.  This script solves every (m,k) pair at each level with
count_sweep.refine_pairs -- the same full-grid search the sweep's below-count
check uses -- at the recorded settings (2001 nodes, refine_outer=300,
remez_iters=60, k <= 8), starting at the recorded stored count minus one and
walking down while a level meets.

The reference is f1_reference.f1_reference, the 60-digit certifying reference
the sweep used.  The acceptance rule is count_sweep.meets inside refine_pairs:
the worse of the fused and separate binary64 evaluations, relative target for
rel:1e-13 and rel:1e-14, absolute target for abs:5e-14.

Usage:
    python knee_walk_f1.py --cell a2-paper-rel13 --jobs 4
    python knee_walk_f1.py --cell a2-2^5-rel14 --jobs 4

Writes _knee-f1/knee_walk_<cell>.json after every level.
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

CELLS = {
    "a2-paper-rel13": {
        "edge_scale": 2,
        "width_name": "paper",
        "target": "rel:1e-13",
        "recorded_stored": 13,
    },
    "a2-2^5-rel14": {
        "edge_scale": 2,
        "width_name": "2^5",
        "target": "rel:1e-14",
        "recorded_stored": 21,
    },
}


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def recorded_row(cell) -> dict:
    doc = json.loads(RECORD.read_text(encoding="utf-8"))
    for row in doc["rows"]:
        if (row["edge"] == "a2" and row["width_name"] == cell["width_name"]
                and row["target"] == cell["target"]):
            return row
    raise SystemExit("record row not found")


def level_pairs(stored: int, k_max: int):
    return [(stored - 1 - k, k)
            for k in range(1, min(k_max, stored - 2) + 1)]


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cell", required=True, choices=sorted(CELLS))
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--nodes", type=int, default=2000)
    ap.add_argument("--refine-outer", type=int, default=300)
    ap.add_argument("--remez-iters", type=int, default=60)
    ap.add_argument("--k-max", type=int, default=None)
    ap.add_argument("--floor", type=int, default=3)
    ap.add_argument("--out", default=None)
    args = ap.parse_args(argv)

    cell = CELLS[args.cell]
    g = cs.load_generator(BOYS)
    g.JOBS = args.jobs
    k_max = g.RAT_K_MAX if args.k_max is None else args.k_max
    target = cs.parse_targets(cell["target"])[0]

    widths = dict(cs.width_set(g, "all"))
    width = widths[cell["width_name"]]
    scale = cell["edge_scale"]
    a = g.X0 * scale
    w_edge = width * scale
    b = a + w_edge
    ts, xs = cs.make_nodes(g, a, b, args.nodes)
    ref = f1r.f1_reference(g, xs)

    row = recorded_row(cell)
    recorded = row["rat"]
    out_path = Path(args.out) if args.out else (
        PAPER / "_knee-f1" / f"knee_walk_{args.cell}.json")

    doc = {
        "cell": args.cell,
        "script_digest": digest(Path(__file__)),
        "count_sweep_digest": digest(PAPER / "count_sweep.py"),
        "f1_reference_digest": digest(PAPER / "f1_reference.py"),
        "record": str(RECORD),
        "record_digest": digest(RECORD),
        "boys_root": str(BOYS),
        "generator_digest": digest(BOYS / cs.GENERATOR_REL),
        "generator_digest_matches_record": (
            digest(BOYS / cs.GENERATOR_REL)
            == json.loads(RECORD.read_text(encoding="utf-8"))["meta"]
            ["generator"]["sha256"]),
        "python": sys.version,
        "platform": platform.platform(),
        "kernel_version": getattr(g, "VERSION", None),
        "function": "F1",
        "edge": "a2",
        "a": g.mp.nstr(a, 17),
        "b": g.mp.nstr(b, 17),
        "width": g.mp.nstr(w_edge, 17),
        "target": target["spec"],
        "settings": {
            "nodes": args.nodes + 1,
            "refine_outer": args.refine_outer,
            "remez_iters": args.remez_iters,
            "k_max": int(k_max),
            "workers": g.worker_count(),
            "acceptance": ("worse of fused and separate binary64 evaluation; "
                           "count_sweep.meets inside refine_pairs"),
        },
        "recorded": {
            "stored": recorded["stored"],
            "m": recorded["m"],
            "k": recorded["k"],
            "solver": recorded["solver"],
            "delivered_rel": recorded["delivered_rel"],
            "delivered_abs": recorded["delivered_abs"],
            "below_meets": recorded.get("below_meets"),
        },
        "levels": [],
        "started_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    }

    def write():
        tmp = out_path.with_suffix(".json.tmp")
        tmp.write_text(json.dumps(doc, indent=1), encoding="utf-8")
        tmp.replace(out_path)

    write()
    print(f"cell {args.cell}: a={doc['a']} b={doc['b']} "
          f"target={target['spec']} recorded stored={recorded['stored']} "
          f"below_meets={recorded.get('below_meets')}", flush=True)
    print(f"workers {g.worker_count()}; nodes {args.nodes + 1}; "
          f"k<={k_max}", flush=True)

    level = recorded["stored"] - 1
    corrected = None
    miss_level = None
    start_all = time.perf_counter()
    while level >= args.floor:
        pairs = level_pairs(level, k_max)
        t0 = time.perf_counter()
        print(f"level {level}: solving {len(pairs)} pair(s) "
              f"{['[%d/%d]' % p for p in pairs]}", flush=True)
        hits = cs.refine_pairs(g, pairs, a, b, ts, xs, ref, target,
                               args.refine_outer, args.remez_iters)
        secs = time.perf_counter() - t0
        hits.sort(key=lambda h: h["rel"] if target["kind"] == "rel"
                  else h["abs"])
        entry = {
            "stored": level,
            "pairs": [f"[{m}/{k}]" for m, k in pairs],
            "n_hits": len(hits),
            "meets": bool(hits),
            "seconds": round(secs, 3),
            "hits": [
                {
                    "m": h["m"], "k": h["k"], "solver": h["name"],
                    "rel": g.mp.nstr(h["rel"], 17),
                    "abs": g.mp.nstr(h["abs"], 17),
                    "at_x_rel": g.mp.nstr(h["at_rel"], 17),
                    "at_x_abs": g.mp.nstr(h["at_abs"], 17),
                }
                for h in hits
            ],
        }
        doc["levels"].append(entry)
        if hits:
            best = entry["hits"][0]
            print(f"level {level}: MEETS ({len(hits)} fit(s)); best "
                  f"[{best['m']}/{best['k']}] {best['solver']} "
                  f"rel={best['rel']} ({secs:.1f} s)", flush=True)
            corrected = level
            level -= 1
            write()
            continue
        print(f"level {level}: misses ({secs:.1f} s)", flush=True)
        miss_level = level
        write()
        break

    doc["walk_seconds"] = round(time.perf_counter() - start_all, 3)
    doc["corrected"] = {
        "stored": corrected,
        "first_missing_level": miss_level,
        "floor": args.floor,
        "truncated_at_floor": (miss_level is None),
        "moved": (None if corrected is None else
                  recorded["stored"] - corrected),
    }
    if corrected is not None:
        best_entry = next(e for e in doc["levels"] if e["stored"] == corrected)
        doc["corrected"]["best"] = best_entry["hits"][0]
        doc["corrected"]["best"]["solver"] = best_entry["hits"][0]["solver"]
    write()
    print(json.dumps(doc["corrected"], indent=1), flush=True)
    print(f"wrote {out_path}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
