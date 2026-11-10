#!/usr/bin/env python3
"""Exact root isolation and interval-wide denominator bounds for accepted fits.

Every accepted rational fit in the deposited records stores its binary64
numerator and denominator coefficients.  This script recovers those stored
values, converts each to an exact rational (the binary64 value is a dyadic
rational, so the conversion is exact), and answers, for the denominator
Q(t) = 1 + q1 t + ... + qk t^k on the fit interval t in [-1, 1]:

  * how many distinct real roots Q has on [-1, 1], with exact isolating
    intervals and multiplicities (Sturm sequences over exact rationals,
    plus Yun's square-free factorisation);
  * a certified lower bound on min |Q| over the whole interval and an upper
    bound on max |Q|, hence a certified lower bound on min|Q|/max|Q|;
  * whether any root in [-1, 1] is off the 2,001- or 4,001-node counting
    grid (an "off-grid pole"), and the distance to the nearest grid node.

The method is exact except for the final decimal rendering of root values.
The interval bounds are Bernstein convex-hull bounds on the isolating
intervals of Q' (the critical points) plus the interval endpoints; every
real extremum of Q is either an endpoint or a critical point, so the bounds
are rigorous.

The W=32, 1e-15, 60-digit control (cell C1) is the worked case: the record
gives no real roots on [-1, 1] and an endpoint ratio that matches the grid
guard.  This script reproduces that and applies the same analysis to every
accepted fit in the corrected x0 record, the second-edge a2 record and the
60-digit control, without assuming the control's behaviour generalises.

Usage:
  python root_isolation.py                 # all blocks, writes records/root_isolation.json
  python root_isolation.py --block corrected
  python root_isolation.py --selftest      # exact-arithmetic unit checks
  python root_isolation.py --plan          # list the fits without solving
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
import sys
import time
import traceback
from fractions import Fraction as F
from pathlib import Path

import mpmath as mp
import numpy as np

HERE = Path(__file__).resolve().parent
PAPER = HERE.parent.parent              # coefficient-count-sweep
CORRECTED = PAPER / "results-corrected" / "corrected-counts.json"
SECOND = PAPER / "results-second-edge-a2" / "second-edge-a2-counts.json"
CONTROL_DIR = PAPER / "control-60dig"
CONTROL_SCRIPT = CONTROL_DIR / "control_60dig_brasil.py"
OUT_DEFAULT = HERE / "records" / "root_isolation.json"
SELFTEST_OUT = HERE / "records" / "selftest.json"
SCHEMA = "c-strong-optional-root-isolation-1"
ROOT_TOL = F(1, 10 ** 30)               # isolating-interval width, exact
GRID_2001 = 2000
GRID_4001 = 4000
DEN_FLOOR = F(1, 1000)                  # the record's guard: 1e-3


# ---------------------------------------------------------------------------
# Exact polynomial arithmetic over Q (coefficient lists, low to high)


def ptrim(p):
    p = list(p)
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return p


def pdeg(p):
    p = ptrim(p)
    return -1 if (len(p) == 1 and p[0] == 0) else len(p) - 1


def peval(p, x):
    acc = F(0)
    for c in reversed(p):
        acc = acc * x + c
    return acc


def pderiv(p):
    p = ptrim(p)
    if len(p) <= 1:
        return [F(0)]
    return ptrim([F(i) * p[i] for i in range(1, len(p))])


def padd(A, B):
    n = max(len(A), len(B))
    out = [F(0)] * n
    for i, a in enumerate(A):
        out[i] += a
    for i, b in enumerate(B):
        out[i] += b
    return ptrim(out)


def psub(A, B):
    return padd(A, [-c for c in B])


def pdivmod(A, B):
    """Exact polynomial division with remainder; returns (quotient, remainder)."""
    A = ptrim(A)
    B = ptrim(B)
    if len(B) == 1 and B[0] == 0:
        raise ZeroDivisionError("polynomial division by zero")
    if len(A) < len(B):
        return [F(0)], A
    q = [F(0)] * (len(A) - len(B) + 1)
    r = A[:]
    while len(r) >= len(B) and not (len(r) == 1 and r[0] == 0):
        d = len(r) - len(B)
        c = r[-1] / B[-1]
        q[d] = c
        for i in range(len(B)):
            r[i + d] -= c * B[i]
        r = ptrim(r)
    return ptrim(q), r


def pmonic(p):
    p = ptrim(p)
    if len(p) == 1 and p[0] == 0:
        return p
    lc = p[-1]
    return [c / lc for c in p]


def pgcd(A, B):
    A = ptrim(A)
    B = ptrim(B)
    while not (len(B) == 1 and B[0] == 0):
        _, r = pdivmod(A, B)
        A, B = B, r
    return pmonic(A)


def yun_squarefree(Q):
    """Yun's algorithm: factors a_i with the roots of multiplicity i."""
    Q = ptrim(Q)
    if pdeg(Q) <= 0:
        return []
    a = Q
    b = pgcd(a, pderiv(a))
    c, _ = pdivmod(a, b)
    d, _ = pdivmod(pderiv(a), b)
    d = psub(d, pderiv(c))
    out = []
    i = 1
    while pdeg(c) > 0 and i < 1000:
        ai = pgcd(c, d)
        if pdeg(ai) > 0:
            out.append((pmonic(ai), i))
        ci, _ = pdivmod(c, ai)
        di, _ = pdivmod(d, ai)
        di = psub(di, pderiv(ci))
        c, d = ci, di
        i += 1
    return out


def sturm_sequence(p):
    p = ptrim(p)
    seq = [p, pderiv(p)]
    while True:
        if len(seq[-1]) == 1 and seq[-1][0] == 0:
            seq.pop()
            break
        _, r = pdivmod(seq[-2], seq[-1])
        if len(r) == 1 and r[0] == 0:
            break
        seq.append([-c for c in r])
    return seq


def sign_variations(seq, x):
    """Number of sign changes of the Sturm sequence at x, zeros skipped."""
    prev = 0
    n = 0
    for p in seq:
        v = peval(p, x)
        if v == 0:
            continue
        cur = 1 if v > 0 else -1
        if prev and cur != prev:
            n += 1
        prev = cur
    return n


def sign_variations_inf(seq, sign):
    """Sign variations at sign*infinity."""
    prev = 0
    n = 0
    for p in seq:
        d = len(p) - 1
        v = p[-1]
        if sign < 0 and d % 2 == 1:
            v = -v
        if v == 0:
            continue
        cur = 1 if v > 0 else -1
        if prev and cur != prev:
            n += 1
        prev = cur
    return n


def real_root_count(p):
    """Number of distinct real roots of p on the whole line."""
    p = ptrim(p)
    if pdeg(p) <= 0:
        return 0
    g = pgcd(p, pderiv(p))
    sq, _ = pdivmod(p, g) if pdeg(g) > 0 else (p, None)
    seq = sturm_sequence(sq)
    return sign_variations_inf(seq, -1) - sign_variations_inf(seq, 1)


def deflate(p, r):
    q, rem = pdivmod(p, [-r, F(1)])
    if not (len(rem) == 1 and rem[0] == 0):
        raise ArithmeticError("deflation by a non-root")
    return ptrim(q)


def isolate_squarefree(p, lo, hi, tol=ROOT_TOL):
    """Isolate the distinct real roots of a square-free p in [lo, hi].

    Returns a list of dicts: {'exact': r} or {'lo': a, 'hi': b}.  Endpoint
    roots and exact roots met during subdivision are deflated; the isolating
    interval brackets exactly one root and is narrower than tol.
    """
    p = ptrim(p)
    roots = []
    if pdeg(p) <= 0:
        return roots
    # Deflate endpoint roots.
    while peval(p, lo) == 0:
        roots.append({"exact": lo})
        p = deflate(p, lo)
    while peval(p, hi) == 0:
        roots.append({"exact": hi})
        p = deflate(p, hi)
    p = ptrim(p)
    if pdeg(p) <= 0:
        return roots
    while True:
        seq = sturm_sequence(p)
        stack = [(lo, hi)]
        restart = False
        while stack and not restart:
            a, b = stack.pop()
            n = sign_variations(seq, a) - sign_variations(seq, b)
            if n <= 0:
                continue
            if n == 1:
                aa, bb = a, b
                while bb - aa > tol:
                    mid = (aa + bb) / 2
                    if peval(p, mid) == 0:
                        roots.append({"exact": mid})
                        p = deflate(p, mid)
                        restart = True
                        break
                    if sign_variations(seq, aa) - sign_variations(seq, mid) >= 1:
                        bb = mid
                    else:
                        aa = mid
                if not restart:
                    roots.append({"lo": aa, "hi": bb})
            else:
                mid = (a + b) / 2
                if peval(p, mid) == 0:
                    roots.append({"exact": mid})
                    p = deflate(p, mid)
                    restart = True
                    break
                stack.append((a, mid))
                stack.append((mid, b))
        if not restart:
            break
    return roots


def bernstein_coeffs(p, a, b):
    """Exact Bernstein coefficients of p on [a, b] (u = (t-a)/(b-a))."""
    p = ptrim(p)
    d = pdeg(p)
    if d <= 0:
        return [p[0] if p else F(0)]
    delta = b - a
    R = [F(0)] * (d + 1)
    for j, c in enumerate(p):
        if c == 0:
            continue
        for i in range(j + 1):
            R[i] += c * math.comb(j, i) * (a ** (j - i)) * (delta ** i)
    bs = []
    for i in range(d + 1):
        s = F(0)
        for j in range(i + 1):
            s += F(math.comb(i, j), math.comb(d, j)) * R[j]
        bs.append(s)
    return bs


def interval_bounds(p, a, b):
    bs = bernstein_coeffs(p, a, b)
    return min(bs), max(bs)


def critical_points(p, lo, hi, tol=ROOT_TOL):
    """Isolate the real roots of p' in [lo, hi]."""
    dp = pderiv(p)
    if pdeg(dp) <= 0:
        return []
    g = pgcd(dp, pderiv(dp))
    sq, _ = pdivmod(dp, g) if pdeg(g) > 0 else (dp, None)
    return isolate_squarefree(sq, lo, hi, tol)


def interval_abs_bounds(Q, lo, hi):
    """Certified bounds on min Q and max Q over [lo, hi] when Q has no roots.

    Every real extremum is an endpoint or a critical point; on each
    isolating interval of Q' the Bernstein hull bounds Q.
    """
    qlo = peval(Q, lo)
    qhi = peval(Q, hi)
    mn = min(qlo, qhi)
    mx = max(qlo, qhi)
    crit = critical_points(Q, lo, hi)
    crit_info = []
    for r in crit:
        if "exact" in r:
            v = peval(Q, r["exact"])
            mn = min(mn, v)
            mx = max(mx, v)
            crit_info.append({"exact": str(r["exact"]), "value": float(v)})
        else:
            a, b = r["lo"], r["hi"]
            l, u = interval_bounds(Q, a, b)
            mn = min(mn, l)
            mx = max(mx, u)
            crit_info.append({
                "lo": str(a), "hi": str(b),
                "q_lower": float(l), "q_upper": float(u),
            })
    return mn, mx, crit_info


def root_value_str(r, digits=30):
    if "exact" in r:
        return mp.nstr(mp.mpf(r["exact"].numerator) / r["exact"].denominator, digits)
    mid = (r["lo"] + r["hi"]) / 2
    return mp.nstr(mp.mpf(mid.numerator) / mid.denominator, digits)


def root_float(r):
    if "exact" in r:
        q = r["exact"]
    else:
        q = (r["lo"] + r["hi"]) / 2
    return float(mp.mpf(q.numerator) / q.denominator)


def nearest_grid(r, n):
    """Distance from the root midpoint to the nearest node of the n-interval grid."""
    x = root_float(r)
    idx = round((x + 1.0) * n / 2.0)
    idx = max(0, min(n, idx))
    t = 2.0 * idx / n - 1.0
    return abs(x - t), idx


def q_at_grid(Q, i, n):
    return peval(Q, F(2 * i, n) - 1)


# ---------------------------------------------------------------------------
# Fit recovery from the deposited records


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_control_module():
    spec = importlib.util.spec_from_file_location("c60_ri", CONTROL_SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["c60_ri"] = mod
    spec.loader.exec_module(mod)
    return mod


class R:
    pass


def corrected_fits():
    doc = json.loads(CORRECTED.read_text())
    x0 = float(doc["problem"]["x0"])
    widths = {k: float(v) for k, v in doc["problem"]["widths"].items()}
    for idx, c in enumerate(doc["cells"]):
        yield {
            "block": "corrected-x0",
            "index": idx,
            "source": str(CORRECTED.relative_to(PAPER)).replace("\\", "/"),
            "target": c["target"],
            "width_name": c["width_name"],
            "width": widths[c["width_name"]],
            "a": x0,
            "pair": c["pair"],
            "m": c["pair"][0],
            "k": c["pair"][1],
            "stored": c["corrected_count"],
            "P": [F(float(s)) for s in c["witness_P"]],
            "Q": [F(float(s)) for s in c["witness_Q"]],
            "P_str": list(c["witness_P"]),
            "Q_str": list(c["witness_Q"]),
            "record_guard": c.get("guard"),
            "record_guard_2001": c.get("guard_2001"),
            "record_delivered": c.get("delivered_2001"),
            "record_pair": c.get("pair"),
            "route": c.get("route"),
            "witness_converged": c.get("witness_converged"),
        }


def second_edge_fits():
    doc = json.loads(SECOND.read_text())
    a = float(doc["problem"]["a2"])
    for idx, c in enumerate(doc["cells"]):
        yield {
            "block": "second-edge-a2",
            "index": idx,
            "source": str(SECOND.relative_to(PAPER)).replace("\\", "/"),
            "target": c["target"],
            "width_name": c["width_name"],
            "width": float(c["width"]),
            "a": a,
            "pair": c["pair"],
            "m": c["pair"][0],
            "k": c["pair"][1],
            "stored": c["corrected_count"],
            "P": [F(float(s)) for s in c["witness_P"]],
            "Q": [F(float(s)) for s in c["witness_Q"]],
            "P_str": list(c["witness_P"]),
            "Q_str": list(c["witness_Q"]),
            "record_guard": c.get("guard"),
            "record_guard_2001": c.get("guard_2001"),
            "record_delivered": c.get("delivered_2001"),
            "record_pair": c.get("pair"),
            "route": c.get("route"),
            "witness_converged": c.get("witness_converged"),
        }


def control_fits():
    mp.mp.dps = 70
    for p in sorted(CONTROL_DIR.glob("cell_*.json")):
        d = json.loads(p.read_text())
        a = d.get("accepted")
        if not a:
            continue
        P70 = a["P70"]
        Q70 = a["Q70"]
        yield {
            "block": "control-60dig",
            "index": p.stem,
            "source": str(p.relative_to(PAPER)).replace("\\", "/"),
            "target": d["cell"]["target"],
            "width_name": d["cell"]["width_name"],
            "width": float(d["cell"]["width"]),
            "a": float(d["cell"]["a"]),
            "pair": [a["m"], a["k"]],
            "m": a["m"],
            "k": a["k"],
            "stored": a["stored"],
            "P": [F(float(mp.mpf(s))) for s in P70],
            "Q": [F(float(mp.mpf(s))) for s in Q70],
            "P_str": list(a["P60"]),
            "Q_str": list(a["Q60"]),
            "P70": list(P70),
            "Q70": list(Q70),
            "record_guard": float(a["guard64"]) if a.get("guard64") is not None else None,
            "record_guard_2001": float(a["guard64"]) if a.get("guard64") is not None else None,
            "record_delivered": float(a["e64"]) if a.get("e64") is not None else None,
            "record_pair": [a["m"], a["k"]],
            "route": a.get("settings"),
            "witness_converged": a.get("converged"),
            "control_record": d,
            "nodes": a.get("nodes"),
            "weights": a.get("weights"),
            "values": a.get("values"),
        }


FIT_SOURCES = {
    "corrected": corrected_fits,
    "second-edge-a2": second_edge_fits,
    "control-60dig": control_fits,
}


# ---------------------------------------------------------------------------
# One fit


def analyse_denominator(Q, label):
    """Exact root structure and certified interval bounds for one denominator."""
    Q = ptrim(Q)
    out = {"label": label, "degree": pdeg(Q), "failures": []}
    try:
        yun = yun_squarefree(Q)
        out["yun_factors"] = [
            {"degree": pdeg(f), "multiplicity": m} for f, m in yun]
        roots_on = []
        total_real = 0
        for factor, mult in yun:
            total_real += real_root_count(factor)
            for r in isolate_squarefree(factor, F(-1), F(1)):
                rec = {
                    "multiplicity": mult,
                    "value": root_value_str(r),
                    "exact": r.get("exact") is not None,
                }
                if "exact" in r:
                    rec["exact_value"] = str(r["exact"])
                else:
                    rec["lo"] = str(r["lo"])
                    rec["hi"] = str(r["hi"])
                    rec["width"] = float(r["hi"] - r["lo"])
                d2, i2 = nearest_grid(r, GRID_2001)
                d4, i4 = nearest_grid(r, GRID_4001)
                rec["nearest_grid_2001"] = d2
                rec["nearest_grid_2001_index"] = i2
                rec["nearest_grid_4001"] = d4
                rec["nearest_grid_4001_index"] = i4
                roots_on.append(rec)
        roots_on.sort(key=lambda r: float(r["value"]))
        out["roots_on_interval"] = roots_on
        out["n_roots_on_interval"] = len(roots_on)
        out["n_real_total"] = total_real
        out["n_roots_outside"] = total_real - len(roots_on)

        # --- off-grid pole test ---------------------------------------------
        # A root is on-grid only if Q at that exact grid node is zero, i.e.
        # only if the exact root value is a node of one of the two grids.  The
        # grid nodes are rational, so this is an exact test.  Each interval
        # root is tested separately: with several roots, one can sit on a node
        # while another lies between nodes, and the fit is then off-grid.
        on_grid = []
        if roots_on:
            for i in range(GRID_2001 + 1):
                if q_at_grid(Q, i, GRID_2001) == 0:
                    on_grid.append({"grid": 2001, "index": i})
            for i in range(GRID_4001 + 1):
                if q_at_grid(Q, i, GRID_4001) == 0:
                    on_grid.append({"grid": 4001, "index": i})
        off_grid = False
        for rec in roots_on:
            exact = rec.get("exact_value")
            rec_on_grid = False
            if exact is not None:
                value = F(exact)
                for n in (GRID_2001, GRID_4001):
                    node_index = (value + 1) * n / 2
                    if node_index.denominator == 1:
                        rec_on_grid = True
                        break
            if not rec_on_grid:
                off_grid = True
        out["on_grid_poles"] = on_grid
        out["on_grid_pole"] = bool(on_grid)
        out["off_grid_pole"] = off_grid

        # --- interval-wide bounds -------------------------------------------
        if roots_on:
            out["min_Q_lower"] = 0.0
            out["max_Q_upper"] = float(
                max(peval(Q, F(-1)), peval(Q, F(1))))
            out["ratio_lower"] = 0.0
            out["interval_bounds"] = "zero: Q has a root on [-1, 1]"
            out["guard_pass_interval"] = False
        else:
            mn, mx, crit_info = interval_abs_bounds(Q, F(-1), F(1))
            out["min_Q_lower"] = float(mn)
            out["max_Q_upper"] = float(mx)
            out["ratio_lower"] = float(mn / mx) if mx > 0 else None
            out["interval_bounds"] = (
                "bernstein hull at the endpoints and at each isolated "
                "critical point")
            out["critical_points"] = crit_info
            out["guard_pass_interval"] = bool(mn / mx >= DEN_FLOOR)

        # --- exact endpoint and sampled extrema -----------------------------
        out["Q_at_minus1"] = float(peval(Q, F(-1)))
        out["Q_at_plus1"] = float(peval(Q, F(1)))
        # Float scan for the sampled ratio at the counting grid; exact
        # evaluation at the selected extrema only (the interval bound above
        # is the rigorous statement).
        qd = [float(c) for c in Q][::-1]
        ts = np.linspace(-1.0, 1.0, GRID_2001 + 1)
        vals = np.polyval(qd, ts)
        imin, imax = int(np.argmin(vals)), int(np.argmax(vals))
        out["grid_2001_exact_min"] = float(q_at_grid(Q, imin, GRID_2001))
        out["grid_2001_exact_max"] = float(q_at_grid(Q, imax, GRID_2001))
        out["ratio_2001_exact_extrema"] = (
            out["grid_2001_exact_min"] / out["grid_2001_exact_max"]
            if out["grid_2001_exact_max"] != 0 else None)
    except Exception as ex:  # noqa: BLE001 - record, never abort the run
        out["failures"].append(f"{type(ex).__name__}: {ex}")
        out["traceback"] = traceback.format_exc()
    return out


def analyse_fit(fit, control_module=None):
    started = time.perf_counter()
    out = {
        "block": fit["block"],
        "index": fit["index"],
        "source": fit["source"],
        "target": fit["target"],
        "width_name": fit["width_name"],
        "width": fit["width"],
        "a": fit["a"],
        "pair": fit["pair"],
        "nominal_m": fit["m"],
        "nominal_k": fit["k"],
        "stored": fit["stored"],
        "route": fit.get("route"),
        "record_guard": fit.get("record_guard"),
        "record_delivered": fit.get("record_delivered"),
        "witness_converged": fit.get("witness_converged"),
        "failures": [],
        "verification": {},
    }
    try:
        P_full = ptrim(fit["P"])
        Q_full = ptrim(fit["Q"])
        m, k = fit["m"], fit["k"]
        # The nominal [m/k] fit stores m+1 numerator coefficients and k
        # denominator coefficients (the constant denominator term is the
        # implicit 1).  The monomial conversion of the barycentric solution
        # can carry roundoff coefficients beyond those nominal degrees; they
        # are kept as the stored/evaluated polynomial and certified below.
        P_nom = ptrim(P_full[:m + 1])
        Q_nom = ptrim(Q_full[:k + 1])
        out["P_stored_length"] = len(P_full)
        out["Q_stored_length"] = len(Q_full)
        out["P_nominal_length"] = len(P_nom)
        out["Q_nominal_length"] = len(Q_nom)
        out["P_extra_count"] = max(0, len(P_full) - (m + 1))
        out["Q_extra_count"] = max(0, len(Q_full) - (k + 1))
        out["P_extra_max"] = (
            float(max((abs(c) for c in P_full[m + 1:]), default=0)))
        out["Q_extra_max"] = (
            float(max((abs(c) for c in Q_full[k + 1:]), default=0)))
        out["stored_count_check"] = (m + 1) + k
        out["stored_count_match"] = out["stored_count_check"] == fit["stored"]
        out["Q_nominal_constant"] = float(Q_nom[0])
        if Q_nom[0] != 1:
            out["failures"].append(
                f"nominal denominator constant {Q_nom[0]} != 1")

        # --- record integrity check on the binary64 storage -----------------
        if control_module is not None and fit["block"] in (
                "corrected-x0", "second-edge-a2"):
            mp.mp.dps = 70
            ts, xs, fs = control_module.build_grid(
                fit["a"], fit["width"], n=GRID_2001)
            Pm = [mp.mpf(s) for s in fit["P_str"]]
            Qm = [mp.mpf(s) for s in fit["Q_str"]]
            ev = control_module.evaluate_stored_binary64(Pm, Qm, ts, fs)
            kind = "rel" if fit["target"].startswith("rel") else "abs"
            got = ev[kind]
            rec = fit.get("record_delivered")
            out["verification"]["delivered_recomputed"] = got
            out["verification"]["delivered_record"] = rec
            out["verification"]["delivered_match"] = (got == rec)
            out["verification"]["delivered_abs_difference"] = (
                abs(got - rec) if rec is not None else None)
            out["verification"]["delivered_rel_difference"] = (
                abs(got - rec) / abs(rec)
                if rec not in (None, 0) else None)
            out["verification"]["delivered_match_within_5e-15"] = (
                rec is not None and abs(got - rec) <= 5e-15)
            g1 = fit.get("record_guard_2001")
            out["verification"]["guard_recomputed"] = ev["den_ratio"]
            out["verification"]["guard_record"] = g1
            out["verification"]["guard_match"] = (
                g1 is not None and ev["den_ratio"] == g1)
            out["verification"]["guard_abs_difference"] = (
                abs(ev["den_ratio"] - g1) if g1 is not None else None)
            out["verification"]["guard_match_within_1e-12"] = (
                g1 is not None and abs(ev["den_ratio"] - g1) <= 1e-12)
        elif fit["block"] == "control-60dig" and control_module is not None:
            mp.mp.dps = 70
            r = R()
            r.nodes = np.array([float(s) for s in fit["nodes"]])
            r.weights = np.array([float(s) for s in fit["weights"]])
            r.values = np.array([float(s) for s in fit["values"]])
            Pm, Qm = control_module.bary_to_monomial(r)
            saved = [mp.nstr(v, 70) for v in Pm]
            out["verification"]["bary_to_P70_match"] = (
                saved == fit["P70"])
            savedq = [mp.nstr(v, 70) for v in Qm]
            out["verification"]["bary_to_Q70_match"] = (
                savedq == fit["Q70"])
            ts, xs, fs = control_module.build_grid(
                fit["a"], fit["width"], n=GRID_2001)
            ev = control_module.evaluate_stored_binary64(Pm, Qm, ts, fs)
            out["verification"]["e64_recomputed"] = ev["rel"]
            out["verification"]["e64_record"] = fit["record_delivered"]
            out["verification"]["e64_match"] = (
                ev["rel"] == fit["record_delivered"])
            out["verification"]["e64_abs_difference"] = abs(
                ev["rel"] - fit["record_delivered"])
            out["verification"]["e64_match_within_5e-15"] = (
                abs(ev["rel"] - fit["record_delivered"]) <= 5e-15)
            out["verification"]["guard64_recomputed"] = ev["den_ratio"]
            out["verification"]["guard64_match"] = (
                ev["den_ratio"] == fit["record_guard"])
            out["verification"]["guard64_abs_difference"] = abs(
                ev["den_ratio"] - fit["record_guard"])
            out["verification"]["guard64_match_within_1e-12"] = (
                abs(ev["den_ratio"] - fit["record_guard"]) <= 1e-12)

        # --- nominal denominator: the accepted [m/k] pair's denominator -----
        nom = analyse_denominator(Q_nom, "nominal")
        out["denominator_nominal"] = nom
        out["yun_factors"] = nom.get("yun_factors")
        out["roots_on_interval"] = nom.get("roots_on_interval", [])
        out["n_roots_on_interval"] = nom.get("n_roots_on_interval", 0)
        out["n_real_total"] = nom.get("n_real_total", 0)
        out["n_roots_outside"] = nom.get("n_roots_outside", 0)
        out["on_grid_poles"] = nom.get("on_grid_poles", [])
        out["on_grid_pole"] = nom.get("on_grid_pole", False)
        out["off_grid_pole"] = nom.get("off_grid_pole", False)
        out["min_Q_lower"] = nom.get("min_Q_lower")
        out["max_Q_upper"] = nom.get("max_Q_upper")
        out["ratio_lower"] = nom.get("ratio_lower")
        out["guard_pass_interval"] = nom.get("guard_pass_interval")

        # --- stored denominator: the full binary64 monomial conversion ------
        extra = Q_full[k + 1:]
        delta = sum((abs(c) for c in extra), F(0))
        out["Q_extra_perturbation_bound"] = float(delta)
        if not extra:
            stored = dict(nom, label="stored (identical to nominal)")
            out["denominator_stored"] = stored
            out["stored_analysis"] = "identical to nominal: no extra Q terms"
            out["stored_ratio_lower"] = nom.get("ratio_lower")
            out["stored_guard_pass_interval"] = nom.get("guard_pass_interval")
            out["stored_n_roots_on_interval"] = nom.get("n_roots_on_interval", 0)
            out["stored_off_grid_pole"] = nom.get("off_grid_pole", False)
        elif (nom.get("n_roots_on_interval", 0) == 0
              and nom.get("min_Q_lower", 0.0) > float(delta)):
            # Q_stored = Q_nom + R with |R(t)| <= sum|extra| on [-1,1].
            # The margin is larger than the perturbation, so Q_stored has
            # the same sign and no roots; the interval bounds follow.
            mn = nom["min_Q_lower"] - float(delta)
            mx = nom["max_Q_upper"] + float(delta)
            out["denominator_stored"] = {
                "label": "stored (nominal plus conversion roundoff)",
                "degree": pdeg(Q_full),
                "analysis": "perturbation certificate",
                "perturbation_bound": float(delta),
                "n_roots_on_interval": 0,
                "off_grid_pole": False,
                "on_grid_pole": False,
                "min_Q_lower": mn,
                "max_Q_upper": mx,
                "ratio_lower": mn / mx if mx > 0 else None,
                "guard_pass_interval": bool(mn / mx >= 1e-3) if mx > 0 else False,
                "extra_roots_added": 0,
            }
            out["stored_analysis"] = (
                "perturbation certificate: |Q_stored-Q_nominal| <= "
                "sum|extra| < min Q_nominal")
            out["stored_ratio_lower"] = mn / mx if mx > 0 else None
            out["stored_guard_pass_interval"] = bool(
                mn / mx >= 1e-3) if mx > 0 else False
            out["stored_n_roots_on_interval"] = 0
            out["stored_off_grid_pole"] = False
        else:
            stored = analyse_denominator(Q_full, "stored (full conversion)")
            out["denominator_stored"] = stored
            out["stored_analysis"] = "direct exact analysis"
            out["stored_ratio_lower"] = stored.get("ratio_lower")
            out["stored_guard_pass_interval"] = stored.get("guard_pass_interval")
            out["stored_n_roots_on_interval"] = stored.get(
                "n_roots_on_interval", 0)
            out["stored_off_grid_pole"] = stored.get("off_grid_pole", False)

        # --- record guard comparison ----------------------------------------
        out["record_guard_pass"] = (
            fit.get("record_guard") is not None
            and fit["record_guard"] >= 1e-3)
        if (fit.get("record_guard") is not None
                and out.get("stored_ratio_lower") is not None):
            out["stored_ratio_gap_to_record_guard"] = (
                fit["record_guard"] - out["stored_ratio_lower"])

    except Exception as ex:  # noqa: BLE001 - record, never abort the run
        out["failures"].append(f"{type(ex).__name__}: {ex}")
        out["traceback"] = traceback.format_exc()
    out["seconds"] = round(time.perf_counter() - started, 4)
    return out


# ---------------------------------------------------------------------------
# Selftests


def selftest():
    cases = []

    def check(name, Q, n_on, total, mults=None, ratio=None):
        yun = yun_squarefree(ptrim(Q))
        roots_on = []
        total_got = 0
        for factor, m in yun:
            total_got += real_root_count(factor)
            for r in isolate_squarefree(factor, F(-1), F(1)):
                roots_on.append((root_value_str(r, 12), m))
        ok = (len(roots_on) == n_on and total_got == total)
        if mults is not None:
            ok = ok and sorted(m for _, m in roots_on) == sorted(mults)
        if ratio is not None and not roots_on:
            mn, mx, _ = interval_abs_bounds(ptrim(Q), F(-1), F(1))
            got = float(mn / mx)
            ok = ok and abs(got - ratio) < 1e-14
            ratio_got = got
        else:
            ratio_got = None
        cases.append({
            "name": name, "ok": ok, "roots_on": roots_on,
            "n_roots_on": len(roots_on), "n_real_total": total_got,
            "ratio": ratio_got,
        })
        return ok

    # 1 + t^2: no real roots, min 1, max 2, ratio 1/2.
    check("1+t^2", [F(1), F(0), F(1)], 0, 0, ratio=0.5)
    # t^2 - 1/4: two roots inside at +/-1/2.
    check("t^2-1/4", [F(-1, 4), F(0), F(1)], 2, 2)
    # t^2 - 1: two endpoint roots.
    check("t^2-1", [F(-1), F(0), F(1)], 2, 2)
    # (t-1/2)^2: one double root inside.
    check("(t-1/2)^2", [F(1, 4), F(-1), F(1)], 1, 1, mults=[2])
    # (t-1/4)(t+3/4) = t^2 + t/2 - 3/16: two simple roots inside.
    check("(t-1/4)(t+3/4)", [F(-3, 16), F(1, 2), F(1)], 2, 2, mults=[1, 1])
    # t: exact root at the initial midpoint 0.
    check("t", [F(0), F(1)], 1, 1, mults=[1])
    # t^2 - 1/2: two interior irrational roots.
    check("t^2-1/2", [F(-1, 2), F(0), F(1)], 2, 2)
    # (t-1/2)^2 (t+3/4) = t^3 - t^2/4 - t/2 + 3/16: roots -3/4 and 1/2.
    check("(t-1/2)^2(t+3/4)", [F(3, 16), F(-1, 2), F(-1, 4), F(1)], 2, 2,
          mults=[1, 2])
    # (t-3)^2 + 1 = t^2 - 6t + 10: no real roots; min 1 at t=3 (outside),
    # so on [-1,1] the minimum is at the right endpoint 5, max at the left 17.
    check("(t-3)^2+1", [F(10), F(-6), F(1)], 0, 0, ratio=5.0 / 17.0)
    # (t-1/3)^2 + 1: interior minimum 1 at t=1/3; ratio 1/((4/3)^2+1)=9/25.
    check("(t-1/3)^2+1", [F(10, 9), F(-2, 3), F(1)], 0, 0, ratio=9.0 / 25.0)
    passed = sum(1 for c in cases if c["ok"])
    return {"cases": cases, "passed": passed, "total": len(cases),
            "all_pass": passed == len(cases)}


# ---------------------------------------------------------------------------


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--block", choices=sorted(FIT_SOURCES) + ["all"],
                    default="all")
    ap.add_argument("--out", type=Path, default=OUT_DEFAULT)
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--plan", action="store_true")
    args = ap.parse_args(argv)

    if args.selftest:
        st = selftest()
        SELFTEST_OUT.parent.mkdir(parents=True, exist_ok=True)
        SELFTEST_OUT.write_text(json.dumps(st, indent=2) + "\n")
        for c in st["cases"]:
            print(f"  {'ok  ' if c['ok'] else 'FAIL'} {c['name']}: "
                  f"roots_on={c['n_roots_on']} total={c['n_real_total']} "
                  f"ratio={c['ratio']}")
        print(f"selftest: {st['passed']}/{st['total']} passed; "
              f"record {SELFTEST_OUT}")
        return 0 if st["all_pass"] else 1

    blocks = (sorted(FIT_SOURCES) if args.block == "all"
              else [args.block])
    fit_specs = []
    for b in blocks:
        fit_specs.extend(list(FIT_SOURCES[b]()))
    if args.plan:
        for f in fit_specs:
            print(f"{f['block']:16s} {f['index']!s:10s} "
                  f"{f['width_name']:6s} {f['target']:12s} "
                  f"[{f['pair'][0]}/{f['pair'][1]}] stored={f['stored']}")
        print(f"fits: {len(fit_specs)}")
        return 0

    control_module = load_control_module()
    started = time.perf_counter()
    results = []
    for n, f in enumerate(fit_specs, 1):
        res = analyse_fit(f, control_module)
        results.append(res)
        flag = []
        if res.get("roots_on_interval"):
            flag.append(f"roots={res['n_roots_on_interval']}")
        if res.get("off_grid_pole"):
            flag.append("OFF-GRID")
        if res.get("on_grid_pole"):
            flag.append("ON-GRID")
        if res.get("failures"):
            flag.append("FAIL")
        print(f"[{n:3d}/{len(fit_specs)}] {res['block']:16s} "
              f"{res['index']!s:10s} {res['width_name']:6s} "
              f"{res['target']:12s} pair=[{res['pair'][0]}/{res['pair'][1]}] "
              f"ratio={res.get('ratio_lower')} "
              f"stored={res.get('stored_ratio_lower')} "
              f"extraQ={res.get('Q_extra_count')} {' '.join(flag)}",
              flush=True)

    summary = {
        "fits": len(results),
        "failures": sum(1 for r in results if r.get("failures")),
        "roots_on_interval": sum(r.get("n_roots_on_interval", 0) for r in results),
        "fits_with_roots": sum(1 for r in results
                               if r.get("n_roots_on_interval", 0) > 0),
        "off_grid_poles": sum(1 for r in results if r.get("off_grid_pole")),
        "on_grid_poles": sum(1 for r in results if r.get("on_grid_pole")),
        "fits_guard_pass_interval": sum(
            1 for r in results if r.get("guard_pass_interval")),
        "fits_guard_fail_interval": sum(
            1 for r in results if r.get("guard_pass_interval") is False),
        "min_ratio_lower": min(
            (r["ratio_lower"] for r in results if r.get("ratio_lower") is not None),
            default=None),
        "max_ratio_lower": max(
            (r["ratio_lower"] for r in results if r.get("ratio_lower") is not None),
            default=None),
        "fits_with_extra_Q_coefficients": sum(
            1 for r in results if r.get("Q_extra_count", 0) > 0),
        "max_Q_extra_coefficients": max(
            (r.get("Q_extra_count", 0) for r in results), default=0),
        "stored_fits_with_roots": sum(
            1 for r in results if r.get("stored_n_roots_on_interval", 0) > 0),
        "stored_off_grid_poles": sum(
            1 for r in results if r.get("stored_off_grid_pole")),
        "stored_guard_pass_interval": sum(
            1 for r in results if r.get("stored_guard_pass_interval")),
        "stored_guard_fail_interval": sum(
            1 for r in results if r.get("stored_guard_pass_interval") is False),
        "min_stored_ratio_lower": min(
            (r["stored_ratio_lower"] for r in results
             if r.get("stored_ratio_lower") is not None), default=None),
        "verification_delivered_mismatches": sum(
            1 for r in results
            if r.get("verification", {}).get("delivered_match") is False),
        "verification_guard_mismatches": sum(
            1 for r in results
            if r.get("verification", {}).get("guard_match") is False),
        "verification_delivered_within_tolerance": sum(
            1 for r in results
            if (r.get("verification", {}).get(
                    "delivered_match_within_5e-15")
                or r.get("verification", {}).get(
                    "e64_match_within_5e-15"))),
        "verification_guard_within_tolerance": sum(
            1 for r in results
            if (r.get("verification", {}).get("guard_match_within_1e-12")
                or r.get("verification", {}).get(
                    "guard64_match_within_1e-12"))),
        "verification_strict_mismatch_cells": [
            {
                "block": r["block"], "index": r["index"],
                "width_name": r["width_name"], "target": r["target"],
                "witness_converged": r.get("witness_converged"),
                "delivered_record": r.get("verification", {}).get(
                    "delivered_record", r.get("verification", {}).get(
                        "e64_record")),
                "delivered_recomputed": r.get("verification", {}).get(
                    "delivered_recomputed", r.get("verification", {}).get(
                        "e64_recomputed")),
                "delivered_abs_difference": r.get("verification", {}).get(
                    "delivered_abs_difference", r.get("verification", {}).get(
                        "e64_abs_difference")),
                "guard_abs_difference": r.get("verification", {}).get(
                    "guard_abs_difference", r.get("verification", {}).get(
                        "guard64_abs_difference")),
            }
            for r in results
            if (r.get("verification", {}).get("delivered_match") is False
                or r.get("verification", {}).get("e64_match") is False
                or r.get("verification", {}).get("guard_match") is False)
        ],
        "seconds": round(time.perf_counter() - started, 2),
    }
    doc = {
        "schema": SCHEMA,
        "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "versions": {
            "python": sys.version.split()[0],
            "mpmath": mp.__version__,
            "numpy": np.__version__,
        },
        "script": {
            "path": str(Path(__file__).resolve()),
            "sha256": sha256(Path(__file__).resolve()),
        },
        "inputs": {
            "corrected": {"path": str(CORRECTED), "sha256": sha256(CORRECTED)},
            "second_edge_a2": {"path": str(SECOND), "sha256": sha256(SECOND)},
            "control_cells": {
                p.name: sha256(p) for p in sorted(CONTROL_DIR.glob("cell_*.json"))
            },
        },
        "method": {
            "coefficient_storage": "binary64, exact dyadic rationals",
            "root_isolation": "Sturm sequences over exact rationals, "
                              "Yun square-free factors, exact deflation",
            "interval_bounds": "Bernstein hull on critical-point isolating "
                               "intervals and at the interval endpoints",
            "root_tolerance_width": "1e-30",
            "grids": [GRID_2001, GRID_4001],
            "denominator_guard": "min|Q|/max|Q| >= 1e-3 over [-1,1]",
        },
        "summary": summary,
        "fits": results,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(doc, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    print(f"record {args.out}")
    return 0 if summary["failures"] == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
