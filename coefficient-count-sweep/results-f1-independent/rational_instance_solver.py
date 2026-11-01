#!/usr/bin/env python3
"""Independent rational counts for F0, F1 and the pure tail at any left edge.

This is a thin driver over the deposited independent solver
`coefficient-count-sweep/results-second-edge-a2/second_edge_a2_solver.py`
(baryrat 2.1.2, the paper's acceptance rule, Horner evaluation, denominator
guard, retry protocol and witness selection, unchanged). Two things are
substituted:

  * the target function: F0 (default), F1 = (F0 - exp(-x))/(2x), or the pure
    tail x^(-1/2);
  * the left edge and the width scaling: the fifteen Group A widths are
    multiplied by a/x0, so W/a is matched at every left edge and the
    cross-ratio is the same.

The point of the driver is the review's two open questions: the F1 and pure
tail rational counts currently come from the kernel's own solver, while F0's
headline counts come from this independent one; and the second-edge design
must be tested where the erfc transient is below the target at both edges.

Usage (one function and one left edge per run):

  C_RATIONAL_FUNCTION=f1 python rational_instance_solver.py \
      --function f1 --left 11.899848152108484 --tag x0 \
      --max-stored 30 --workers 6 --out-dir out_f1_x0

  --only f1_x0_rel13_m7   runs one cell (probe).

The environment variable is what the worker processes read; the --function
flag must agree with it. Records are written under <out-dir>/cells/.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

HERE = Path(__file__).resolve().parent
SOLVER = (HERE.parent / "results-second-edge-a2" / "second_edge_a2_solver.py")

FUNCTION = os.environ.get("C_RATIONAL_FUNCTION", "f0")


def _load_solver():
    spec = importlib.util.spec_from_file_location("deposited_solver", SOLVER)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


orig = _load_solver()
import numpy as np  # noqa: E402
import mpmath as mp  # noqa: E402

X0 = orig.X0
A2 = orig.A2
NAMES15 = orig.NAMES15
TARGETS = orig.TARGETS

_PRISTINE_F0_NP = orig.f0_np
_PRISTINE_F0_MP = orig.f0_mp


def f1_np(x):
    x = np.asarray(x, dtype=float)
    return (_PRISTINE_F0_NP(x) - np.exp(-x)) / (2.0 * x)


def f1_mp(x):
    x = mp.mpf(x)
    return (_PRISTINE_F0_MP(x) - mp.e ** (-x)) / (2 * x)


def tail_np(x):
    return 1.0 / np.sqrt(np.asarray(x, dtype=float))


def tail_mp(x):
    return 1 / mp.sqrt(mp.mpf(x))


def apply_function(function):
    if function == "f0":
        orig.f0_np = _PRISTINE_F0_NP
        orig.f0_mp = _PRISTINE_F0_MP
    elif function == "f1":
        orig.f0_np = f1_np
        orig.f0_mp = f1_mp
    elif function == "tail":
        orig.f0_np = tail_np
        orig.f0_mp = tail_mp
    else:
        raise SystemExit(f"unknown function {function!r}")


apply_function(FUNCTION)


def load_widths(corrected_counts):
    doc = json.loads(Path(corrected_counts).read_text())
    return {name: float(doc["problem"]["widths"][name]) for name in NAMES15}


def build_worklist(args, widths):
    a = float(args.left)
    scale = a / X0 if args.match_cross_ratio else 1.0
    left_tag = args.tag or f"a{a:.6f}".replace(".", "p")
    work = []
    for tag, tlabel, kind, target in TARGETS:
        for name in NAMES15:
            cell_id = f"{args.function}_{left_tag}_{tag}_{name}"
            work.append({
                "cell_id": cell_id,
                "edge": left_tag,
                "group": "A",
                "kind": kind,
                "target": target,
                "a": a,
                "W": widths[name] * scale,
                "width_name": name,
                "target_tag": tag,
                "x0_count_used_for_bound": None,
                "max_stored": args.max_stored,
                "out": f"cells/{cell_id}.json",
            })
    if args.only:
        keep = set(args.only.split(","))
        work = [w for w in work if w["cell_id"] in keep]
    return work


def task(spec):
    apply_function(FUNCTION)
    return orig.task(spec)


def main():
    global FUNCTION
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--function", choices=("f0", "f1", "tail"), default=FUNCTION)
    ap.add_argument("--left", default=str(X0),
                    help="left edge a; x0 is 11.899848152108484, a2 is 23.799696304216968")
    ap.add_argument("--tag", default=None, help="label for the left edge")
    ap.add_argument("--match-cross-ratio", action="store_true", default=True,
                    help="scale every width by a/x0 (default on)")
    ap.add_argument("--no-match-cross-ratio", dest="match_cross_ratio",
                    action="store_false")
    ap.add_argument("--max-stored", type=int, default=None)
    ap.add_argument("--only", default=None)
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--out-dir", default=str(HERE))
    ap.add_argument("--corrected-counts", default=str(
        SOLVER.resolve().parent.parent / "results-corrected"
        / "corrected-counts.json"))
    args = ap.parse_args()

    FUNCTION = args.function
    os.environ["C_RATIONAL_FUNCTION"] = FUNCTION
    apply_function(FUNCTION)
    if args.max_stored is None:
        args.max_stored = {"f0": 26, "f1": 30, "tail": 30}[FUNCTION]

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    os.chdir(out_dir)
    widths = load_widths(args.corrected_counts)
    work = build_worklist(args, widths)
    if not work:
        raise SystemExit("no cells selected")

    left_tag = args.tag or f"a{float(args.left):.6f}".replace(".", "p")
    results = []
    t0 = time.time()
    if args.workers <= 1:
        for spec in work:
            rec = task(spec)
            results.append(rec)
            print(f"DONE {len(results)}/{len(work)} {rec.get('cell_id')} "
                  f"count={rec.get('corrected_count')} pair={rec.get('pair')} "
                  f"route={rec.get('route')}", flush=True)
    else:
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            futs = {pool.submit(task, w): w for w in work}
            for i, fut in enumerate(as_completed(futs), 1):
                spec = futs[fut]
                try:
                    rec = fut.result()
                except Exception as ex:  # noqa: BLE001
                    import traceback
                    rec = {"cell_id": spec.get("cell_id"), "error": repr(ex),
                           "traceback": traceback.format_exc(), "spec": spec}
                results.append(rec)
                print(f"DONE {i}/{len(work)} {rec.get('cell_id')} "
                      f"count={rec.get('corrected_count')} "
                      f"pair={rec.get('pair')} route={rec.get('route')}",
                      flush=True)
    results.sort(key=lambda r: r.get("cell_id", ""))
    summary = {
        "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "function": FUNCTION,
        "left": float(args.left),
        "left_tag": left_tag,
        "match_cross_ratio": args.match_cross_ratio,
        "max_stored": args.max_stored,
        "wall_seconds": round(time.time() - t0, 1),
        "workers": args.workers,
        "cells": results,
    }
    name = f"run_{FUNCTION}_{left_tag}.json"
    with open(name, "w") as fh:
        json.dump(summary, fh, indent=1)
    bad = [r for r in results if r.get("error")]
    print(f"ALLDONE cells={len(results)} errors={len(bad)} "
          f"wall={summary['wall_seconds']}s", flush=True)
    return 0 if not bad else 2


if __name__ == "__main__":
    raise SystemExit(main())
