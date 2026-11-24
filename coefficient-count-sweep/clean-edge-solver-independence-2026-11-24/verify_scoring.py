"""Independent re-score of the BRASIL accepted fits from the 60-digit control.

Building-block check for the independent Remez re-solve: take the P60/Q60 that
the control stored for the accepted row, evaluate the relative error against
F_0 at 60 decimal digits on the control's 2,001-node grid, its 4,001-node
check grid, and a denser 20,001-node grid.  This replicates the control's
score from the stored coefficients and tests whether the acceptance is a grid
artefact.

Run:  python -I verify_scoring.py
"""

from __future__ import annotations

import json
from pathlib import Path

import mpmath as mp

CELLS = Path(r"C:\Users\myaml\dev\_c60-clean-edge-full-2026-11-24\cells")
OUT = Path(__file__).resolve().parent / "out"
EVAL_DPS = 60
REF_DPS = 70


def f0_mp(x):
    return mp.sqrt(mp.pi) / (2 * mp.sqrt(x)) * mp.erf(mp.sqrt(x))


def grid(a, W, n):
    # Replicates the control's build_grid: t, x, F_0(x) at REF_DPS.
    ts = []
    fs = []
    for i in range(n + 1):
        t = 2.0 * i / n - 1.0          # binary64, as in the control
        x = a + W * (t + 1.0) / 2.0    # binary64, as in the control
        ts.append(t)
        fs.append(f0_mp(mp.mpf(float(x))))
    return ts, fs


def mp_horner(cs, t):
    acc = mp.mpf(0)
    for c in reversed(cs):
        acc = acc * t + c
    return acc


def score(P60, Q60, ts, fs, dps=EVAL_DPS):
    """Worst relative error and denominator ratio, exactly as the control."""
    with mp.workdps(dps):
        P = [mp.mpf(c) for c in P60]
        Q = [mp.mpf(c) for c in Q60]
        wr = mp.mpf(0)
        wa = mp.mpf(0)
        at = None
        dmin = mp.inf
        dmax = mp.mpf(0)
        for t, f in zip(ts, fs):
            tm = mp.mpf(float(t))
            fm = mp.mpf(f)
            num = mp_horner(P, tm)
            den = mp_horner(Q, tm)
            v = num / den
            e = abs(v - fm)
            if e > wa:
                wa = e
            rel = e / fm
            if rel > wr:
                wr, at = rel, t
            d = abs(den)
            if d < dmin:
                dmin = d
            if d > dmax:
                dmax = d
    return {"e60": wr, "e60_abs": wa, "at": at,
            "guard": dmin / dmax}


def main():
    mp.mp.dps = REF_DPS
    report = {}
    for name in ("p5_rel_1e-15", "ship_rel_1e-15"):
        doc = json.loads((CELLS / f"{name}.json").read_text(encoding="utf-8"))
        cell = doc["cell"]
        acc = doc["accepted"]
        a = float(cell["a"])
        W = float(cell["width"])
        target = float(cell["target_value"])
        P60 = acc["P60"]
        Q60 = acc["Q60"]
        P70 = acc["P70"]
        Q70 = acc["Q70"]
        # P60 must be the 60-significant-digit rounding of the 70-digit value.
        with mp.workdps(EVAL_DPS):
            P60_ref = [mp.nstr(mp.mpf(c), EVAL_DPS) for c in P70]
            Q60_ref = [mp.nstr(mp.mpf(c), EVAL_DPS) for c in Q70]
        rounding_ok = (P60 == P60_ref and Q60 == Q60_ref)

        row = {"cell": name, "pair": [acc["m"], acc["k"]],
               "stored": acc["stored"], "target": target,
               "stored_e60": acc["e60"], "stored_guard60": acc["guard60"],
               "P60_is_rounding_of_P70": rounding_ok, "grids": {}}
        for n, label in ((2000, "2001"), (4000, "4001"), (20000, "20001")):
            ts, fs = grid(a, W, n)
            s = score(P60, Q60, ts, fs)
            row["grids"][label] = {
                "e60": mp.nstr(s["e60"], 17),
                "e60_abs": mp.nstr(s["e60_abs"], 17),
                "at_t": repr(s["at"]),
                "guard60": mp.nstr(s["guard"], 17),
                "meets": bool(s["e60"] <= target
                              and s["guard"] >= 1e-3),
            }
        report[name] = row
        print(f"[{name}] pair [{acc['m']}/{acc['k']}] stored e60="
              f"{acc['e60']} guard={acc['guard60']} "
              f"P60=round(P70):{rounding_ok}")
        for label, s in row["grids"].items():
            print(f"   {label:>5s} nodes: e60={s['e60']} at t={s['at_t']} "
                  f"guard={s['guard60']} meets={s['meets']}")
    OUT.mkdir(exist_ok=True)
    (OUT / "verify_scoring.json").write_text(json.dumps(report, indent=1),
                                            encoding="utf-8")
    print(f"wrote {OUT / 'verify_scoring.json'}")


if __name__ == "__main__":
    main()
