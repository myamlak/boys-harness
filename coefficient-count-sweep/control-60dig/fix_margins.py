#!/usr/bin/env python3
"""Recompute the below-count margins from admissible fits only.

A fit whose 60-digit denominator guard is below 1e-3 is not admissible under
the record's acceptance rule, so it does not bound the below-count margin even
if its error is below the target. The driver's first revision reported the
smallest e60 at the failing level without that filter. This script rewrites
`failing_level.routes[*].best_e60` and adds `best_admissible_e60` and
`best_admissible_ratio` from the rows stored in each record. It is idempotent
and does not touch any other field.

Usage: python fix_margins.py [cell_*.json ...]
"""

import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
DEN_FLOOR = 1e-3


def main(argv):
    paths = [Path(p) for p in argv] if argv else sorted(
        HERE.glob("cell_*.json"))
    for p in paths:
        d = json.loads(p.read_text())
        fail = d.get("failing_level")
        if not fail:
            print(f"{p.name}: no failing level")
            continue
        target = d["cell"]["target_value"]
        rows = fail.get("rows") or []
        adm = [r for r in rows if r.get("e60")
               and float(r.get("guard60", 0.0)) >= DEN_FLOOR]
        best = min((float(r["e60"]) for r in adm), default=None)
        for route in fail.get("routes", []):
            shared = [r for r in rows
                      if r.get("settings") == route["pass"]
                      and r.get("e60")
                      and float(r.get("guard60", 0.0)) >= DEN_FLOOR]
            rb = min((float(r["e60"]) for r in shared), default=None)
            route["best_admissible_e60"] = repr(rb) if rb is not None else None
            route["best_admissible_ratio"] = (rb / target
                                               if rb is not None else None)
            route["best_e60_unfiltered"] = route.get("best_e60")
            route["best_e60"] = repr(rb) if rb is not None else None
        fail["best_admissible_e60"] = repr(best) if best is not None else None
        fail["best_admissible_ratio"] = (best / target
                                         if best is not None else None)
        p.write_text(json.dumps(d, indent=1))
        print(f"{p.name}: best admissible e60 "
              f"{fail['best_admissible_e60']} "
              f"({fail['best_admissible_ratio']:.3f}x target) at stored "
              f"{fail['stored']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
