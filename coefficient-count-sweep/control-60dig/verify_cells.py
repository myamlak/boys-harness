#!/usr/bin/env python3
"""Re-evaluate the accepted fits from their stored barycentric data.

For every cell record in this directory, rebuild the barycentric rational from
`accepted.nodes`, `accepted.weights` and `accepted.values` (the solver's own
output, stored as binary64 reprs), reconvert exactly to monomials at 70 digits,
and recompute `e60`, `e64r` and `e64` on the 2,001-node grid.  The recomputed
values must equal the stored ones; a mismatch means the record's coefficients
or arithmetic do not describe the fit it reports.

This is a check of the record, not a new measurement: no solver is run.

Usage: python verify_cells.py
"""

import importlib.util
import json
import sys
from pathlib import Path

import mpmath as mp
import numpy as np

HERE = Path(__file__).resolve().parent


def load_driver():
    spec = importlib.util.spec_from_file_location(
        "c60", HERE / "control_60dig_brasil.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules["c60"] = mod
    spec.loader.exec_module(mod)
    return mod


class R:
    pass


def main():
    drv = load_driver()
    rec = json.loads(drv.CORRECTED.read_text())
    widths = {n: float(rec["problem"]["widths"][n])
              for n in ("paper", "ship", "p5", "p0")}
    x0 = float(rec["problem"]["x0"])
    ok = True
    for p in sorted(HERE.glob("cell_*.json")):
        d = json.loads(p.read_text())
        if not d.get("accepted"):
            print(f"{p.name}: no accepted fit")
            ok = False
            continue
        a = d["accepted"]
        c = d["cell"]
        mp.mp.dps = drv.CONVERT_DPS
        ts, xs, fs = drv.build_grid(x0, widths[c["width_name"]])
        r = R()
        r.nodes = np.array([float(s) for s in a["nodes"]])
        r.weights = np.array([float(s) for s in a["weights"]])
        r.values = np.array([float(s) for s in a["values"]])
        P, Q = drv.bary_to_monomial(r)
        saved = [mp.nstr(v, drv.CONVERT_DPS) for v in P]
        mismatch = any(x != y for x, y in zip(saved, a["P70"])) or \
            len(saved) != len(a["P70"])
        got = drv.score_at_60(P, Q, ts, fs, c["target_value"])
        fields = ("e60", "e64r", "e64", "guard60")
        same = all(got[k] == a[k] for k in fields)
        print(f"{p.name}: pair [{a['m']}/{a['k']}] "
              f"e60 {got['e60']} {'==' if got['e60'] == a['e60'] else '!='} "
              f"{a['e60']}; e64r {'==' if got['e64r'] == a['e64r'] else '!='}; "
              f"e64 {'==' if got['e64'] == a['e64'] else '!='}; "
              f"P70 {'==' if not mismatch else '!='}")
        ok = ok and same and not mismatch
    print("ALL MATCH" if ok else "MISMATCH")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
