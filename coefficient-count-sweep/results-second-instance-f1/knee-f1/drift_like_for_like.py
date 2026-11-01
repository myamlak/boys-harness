"""The tail instance's constant-drift statistic, applied to F1.

The tail report's statistic (`analyze_tail_instance.py`): one count-implied
constant per edge, target and family from the rows whose stored count is above
5; the constant of a row is C = exp((N-1)*rate - L); the drift is the range of
log C over those rows, in decades. This script computes the same statistic on
an F1 sweep record so the two instances can be compared like for like.

Usage:
    python _knee-f1/drift_like_for_like.py --sweep <sweep_all.json> [--json out]
"""
import argparse
import json
import math
import sys
from pathlib import Path

import mpmath as mp

HERE = Path(__file__).resolve().parent
PAPER = HERE.parent
sys.path.insert(0, str(PAPER))

from analyze_second_instance_f1 import (  # noqa: E402
    DYAD, TARGETS, f1_edge_maps, kappa, load_f1, mu_of, target_L,
)

mp.mp.dps = 60

EDGES = {"x0": mp.mpf("11.899848152108484"),
         "a2": mp.mpf("23.799696304216968")}
# F_1 at the two edges, from the paper's reference routine at 60 digits
# (section 2 of SECOND-INSTANCE-F1.md).
F1_EDGE = {"x0": mp.mpf("0.010794210184588879"),
           "a2": mp.mpf("0.003816434577613709")}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sweep", required=True)
    ap.add_argument("--json")
    args = ap.parse_args()

    doc, rows = load_f1(Path(args.sweep))
    out = {"sweep": str(Path(args.sweep)),
           "source_record_digest": doc.get("source_record_digest")
           or doc.get("record_digest"),
           "rule": "stored > 5; C = exp((stored-1)*rate - L); drift = range of log C",
           "blocks": {}}
    for edge, a in EDGES.items():
        m = f1_edge_maps(rows, edge)
        fval = F1_EDGE[edge]
        L = {t: target_L(t, a, fval) for t in TARGETS}
        out["blocks"][edge] = {}
        for target in TARGETS:
            out["blocks"][edge][target] = {}
            for fam in ("cheb", "rat"):
                logcs, widths = [], []
                for name in DYAD:
                    e = m[fam].get((target, name), {})
                    n = e.get("stored")
                    if n is None or n <= 5:
                        continue
                    if fam == "rat" and e.get("below_meets"):
                        continue
                    k1 = kappa(a, mp.mpf(str(rows[(edge, target, name)]["width"])))
                    rate = mp.log(1 / k1) if fam == "cheb" else mu_of(k1)
                    logcs.append((float(n) - 1) * float(rate) - float(L[target]))
                    widths.append(name)
                if not logcs:
                    out["blocks"][edge][target][fam] = {"status": "no informative rows"}
                    continue
                logcs.sort()
                mid = logcs[len(logcs) // 2]
                C = math.exp(mid)
                resid = []
                for name, lc in zip(widths, [None] * len(widths)):
                    e = m[fam].get((target, name), {})
                    n = e.get("stored")
                    k1 = kappa(a, mp.mpf(str(rows[(edge, target, name)]["width"])))
                    rate = mp.log(1 / k1) if fam == "cheb" else mu_of(k1)
                    real = 1 + (L[target] + mp.log(mp.mpf(repr(C)))) / rate
                    resid.append((name, n, float(real), float(n - real)))
                block = {
                    "n_informative": len(logcs),
                    "informative_widths": widths,
                    "C_median": C,
                    "C_min": math.exp(logcs[0]),
                    "C_max": math.exp(logcs[-1]),
                    "drift_decades": (logcs[-1] - logcs[0]) / math.log(10),
                    "max_abs_residual": max(abs(r[3]) for r in resid),
                    "inside_paper_band": math.exp(logcs[0]) >= 1e-2
                    and math.exp(logcs[-1]) <= 1e4,
                    "residuals": [{"width": r[0], "stored": r[1],
                                   "real": r[2], "residual": r[3]}
                                  for r in resid],
                }
                out["blocks"][edge][target][fam] = block

    text = json.dumps(out, indent=2)
    if args.json:
        Path(args.json).write_text(text + "\n", encoding="utf-8")
    for edge in EDGES:
        for target in TARGETS:
            for fam in ("cheb", "rat"):
                b = out["blocks"][edge][target][fam]
                if "C_median" not in b:
                    continue
                print(f"{edge:2s} {target:12s} {fam:4s} n={b['n_informative']:2d} "
                      f"C={b['C_median']:.5g} [{b['C_min']:.3g},{b['C_max']:.3g}] "
                      f"drift={b['drift_decades']:.2f} dec "
                      f"maxres={b['max_abs_residual']:.3f} "
                      f"band={b['inside_paper_band']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
