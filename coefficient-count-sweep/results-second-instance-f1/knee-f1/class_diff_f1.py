#!/usr/bin/env python3
"""Direct class test: the per-width count difference F1 minus F0.

If the measured width-law rate does not depend on alpha, the F1 and F0 counts
at matched widths and edges differ by a constant (the two tail constants
differ); the difference regressed on log2 W must have slope zero.  This test
does not use the model's constants at all.

F0 counts are the paper's own records: chebyshev_recount.json and the core
sweep for x0, results-second-edge/sweep_all.json for a2.  Rational cells with
below_meets true are excluded at either function, as the paper's analysis
excludes them.

Usage:
    python class_diff_f1.py [--sweep <record>] [--out <json>]
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

PAPER = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PAPER))

import analyze_second_instance_f1 as af  # noqa: E402


def main(argv=None) -> int:
    import argparse
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sweep", default=str(af.DEFAULT_F1))
    ap.add_argument("--out", default=str(PAPER / "_knee-f1" / "class_diff.json"))
    args = ap.parse_args(argv)

    doc, rows = af.load_f1(Path(args.sweep))
    f0_x0 = af.load_f0_x0(af.DEFAULT_F0_RECOUNT, af.DEFAULT_F0_CORE)
    f0_a2 = af.load_f0_a2(af.DEFAULT_F0_A2)
    out = {"sweep": args.sweep, "blocks": {}}

    for edge, f0 in (("x0", f0_x0), ("a2", f0_a2)):
        f1 = af.f1_edge_maps(rows, edge)
        out["blocks"][edge] = {}
        for target in af.TARGETS:
            out["blocks"][edge][target] = {}
            for fam in ("cheb", "rat"):
                xs, ys, used, dropped = [], [], [], []
                y1s, y0s = [], []
                for j, name in enumerate(af.DYAD):
                    e1 = f1[fam].get((target, name))
                    e0 = f0[fam].get((target, name))
                    if e1 is None or e1.get("stored") is None:
                        dropped.append([name, "f1 missing"])
                        continue
                    if e0 is None or e0.get("stored") is None:
                        dropped.append([name, "f0 missing"])
                        continue
                    if fam == "rat" and (e1.get("below_meets")
                                         or e0.get("below_meets")):
                        dropped.append([name, "below_meets"])
                        continue
                    xs.append(float(j))
                    ys.append(e1["stored"] - e0["stored"])
                    y1s.append(af.math.log2(e1["stored"]))
                    y0s.append(af.math.log2(e0["stored"]))
                    used.append(name)
                fit = af.ols(xs, ys)
                f1fit = af.ols(xs, y1s)
                f0fit = af.ols(xs, y0s)
                paired = af.paired_slope_test(f0fit, f1fit)
                mean = sum(ys) / len(ys)
                out["blocks"][edge][target][fam] = {
                    "n": len(ys), "widths": used, "dropped": dropped,
                    "diffs": ys, "mean": mean,
                    "slope": fit["beta"], "se": fit["se"],
                    "z": (fit["beta"] / fit["se"]) if fit["se"] else None,
                    "max_dev_from_mean": max(abs(y - mean) for y in ys),
                    "loglog": {
                        "f1_slope": f1fit["beta"], "f1_se": f1fit["se"],
                        "f0_slope": f0fit["beta"], "f0_se": f0fit["se"],
                        "diff": paired["diff"],
                        "se_unpaired": paired["se_unpaired"],
                        "se_paired": paired["se_paired"],
                        "z_unpaired": paired["z_unpaired"],
                        "z_paired": paired["z_paired"],
                        "residual_cov": paired["residual_cov"],
                    },
                }
    Path(args.out).write_text(json.dumps(out, indent=1), encoding="utf-8")
    for edge in ("x0", "a2"):
        for target in af.TARGETS:
            for fam in ("cheb", "rat"):
                b = out["blocks"][edge][target][fam]
                ll = b["loglog"]
                z = f"{b['z']:+.2f}" if b["z"] is not None else "n/a"
                print(f"{edge} {target:10s} {fam:4s}: n={b['n']:2d} "
                      f"count diff mean {b['mean']:+.2f} slope {b['slope']:+.3f}"
                      f" +/-{b['se']:.3f} ({z} s.e.); "
                      f"loglog F1 {ll['f1_slope']:.4f}+/-{ll['f1_se']:.4f} "
                      f"F0 {ll['f0_slope']:.4f}+/-{ll['f0_se']:.4f} "
                      f"diff {ll['diff']:+.4f} unpaired {ll['z_unpaired']:+.2f} "
                      f"paired {ll['z_paired']:+.2f}; diffs {b['diffs']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
