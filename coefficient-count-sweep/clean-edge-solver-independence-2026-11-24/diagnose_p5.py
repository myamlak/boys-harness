"""Diagnose the rational Remez landscape for p5 [6/8].

Compare the error profile of the BRASIL accepted fit (verified 6.48e-16) with
the differential-correction fit found from a Chebyshev initial reference, and
test whether the exact Newton equations converge from the BRASIL fit.
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import mpmath as mp

sys.path.insert(0, str(Path(__file__).resolve().parent))
import remez_rational as R

mp.mp.dps = R.DPS
HERE = Path(__file__).resolve().parent
CELL = Path(r"C:\Users\myaml\dev\_c60-clean-edge-full-2026-11-24\cells"
            r"\p5_rel_1e-15.json")

doc = json.loads(CELL.read_text(encoding="utf-8"))
acc = doc["accepted"]
a = float(doc["cell"]["a"])
W = float(doc["cell"]["width"])
prob = R.problem_f0(a, W)
ts, tsm, fs = R.build_grid(a, W, 2000)
ws = [prob.weight(t, f) for t, f in zip(tsm, fs)]

P = [mp.mpf(c) for c in acc["P70"]]
Q = [mp.mpf(c) for c in acc["Q70"]]


def report(name, P, Q):
    errs, qvals = R.profile(P, Q, tsm, fs, ws)
    gm = max(abs(e) for e in errs)
    at = max(range(len(errs)), key=lambda i: abs(errs[i]))
    g = R.guard_ratio(qvals)
    print(f"\n{name}: grid_max={mp.nstr(gm, 17)} at t={ts[at]} "
          f"guard={mp.nstr(g, 10)}")
    merged = R.refine_extrema(P, Q, ts, errs, prob)
    print(f"  alternating extrema after merge: {len(merged)}")
    for t, e in merged:
        print(f"    t={mp.nstr(t, 12):>18s}  e={mp.nstr(e, 8)}")
    return errs, merged


errs_b, merged_b = report("BRASIL [6/8]", P, Q)
N = 6 + 8 + 2
ref = R.choose_window(merged_b, N)
print(f"\nwindow of {N} extrema from BRASIL profile: {len(ref)}")
print("  " + " ".join(mp.nstr(t, 10) for t in ref))

for label, P0, Q0 in (("BRASIL init", P, Q),):
    Pn, Qn, En, ok, note, res = R.newton_solve(6, 8, ref, prob, P0, Q0,
                                               mp.mpf("0"))
    print(f"\nNewton from {label}: ok={ok} note={note} res={res}")
    if ok:
        errs, qvals = R.profile(Pn, Qn, tsm, fs, ws)
        gm = max(abs(e) for e in errs)
        print(f"  grid_max={mp.nstr(gm, 17)} guard="
              f"{mp.nstr(R.guard_ratio(qvals), 10)} E={mp.nstr(En, 10)}")
        report("Newton fit", Pn, Qn)

# exact residual of BRASIL fit on the window
Rres = []
for i, t in enumerate(ref):
    f = prob.f(t)
    w = prob.weight(t, f)
    s = 1 if i % 2 == 0 else -1
    Pi = R.poly_eval(P, t)
    Qi = R.poly_eval(Q, t)
    Rres.append(w * (Pi / Qi - f) - s * mp.mpf(0))
print("\nBRASIL residual on window with E=0:", mp.nstr(max(abs(r) for r in Rres), 6))
