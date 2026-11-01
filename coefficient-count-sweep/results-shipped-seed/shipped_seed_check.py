#!/usr/bin/env python3
"""Test whether the shipped rational seeds store more than they need.

The shipped region-B rational seed (coarsest, 13 stored) and the narrow
rational route's five region-B pieces (39 stored in total) are re-fitted with
the independent baryrat solver at the library's own bars, using the deposited
solver's acceptance rule and retry protocol unchanged.

Targets are absolute: 2.56905715673095443e-14, the region-B bar
(kNarrowRatBBound and kRegionBFitBar), applied on each interval.
"""
import importlib.util, json, time
from pathlib import Path

SOLVER = Path(__file__).resolve().parent.parent / "results-second-edge-a2" / "second_edge_a2_solver.py"
spec = importlib.util.spec_from_file_location("dep", SOLVER)
dep = importlib.util.module_from_spec(spec); spec.loader.exec_module(dep)

OUT = Path(__file__).resolve().parent
OUT.mkdir(exist_ok=True)
BAR = 2.56905715673095443e-14
EDGES = [11.899848152108484, 14.756597293654865, 18.161791411646277,
         22.253063161220702, 27.200530106782374, 31.0]
CURRENT = [9, 7, 8, 8, 7]  # narrow rational region-B pieces, stored

def run(name, a, W, current, extra=2):
    s = {"cell_id": name, "edge": "narrowB", "group": "B", "kind": "abs",
         "target": BAR, "a": a, "W": W, "width_name": name,
         "target_tag": "abs257", "x0_count_used_for_bound": current,
         "max_stored": current + extra,
         "out": str(OUT / f"{name}.json")}
    t = time.time()
    rec = dep.finalize(s)
    print(f"{name}: current {current} -> independent {rec.get('corrected_count')} "
          f"pair {rec.get('pair')} route {rec.get('route')} "
          f"below_ratio {rec.get('below_ratio')} ({time.time()-t:.0f}s)", flush=True)
    return rec

print("narrow rational region-B pieces, bar 2.569e-14 absolute")
results = []
for i in range(5):
    a, b = EDGES[i], EDGES[i+1]
    results.append(run(f"narrowB_p{i}", a, b - a, CURRENT[i]))
print("narrow total:", sum(CURRENT), "->",
      sum((r.get("corrected_count") or 999) for r in results))

print("\ncoarsest rational region-B seed on [x0, 31], bar 2.569e-14 absolute")
r = run("coarsestB", EDGES[0], EDGES[-1] - EDGES[0], 13, extra=2)
print("coarsest:", 13, "->", r.get("corrected_count"))
